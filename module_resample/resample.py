#resample 
import os
import numpy as np
from scipy import ndimage
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from spectres import spectres
from mpdaf.obj import Cube

# --- constants and configuration ---
c_kms = 299792.458                       # speed of light [km/s]
oII_rest = 3728.48                        # effective [OII] doublet center [Angstrom, vac].
                                         # (3727.092 + 3729.875)/2; treat the blend as one line.
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)    # one cosmology object, reused everywhere
flux_unit = 1e-20                        # MUSE BUNIT scale: erg/s/cm^2/Angstrom per spaxel





#helper function to crop or pad a cube to a fixed size, keeping the center fixed
def _center_crop_or_pad(cube, ny_out, nx_out):
    """Center the last two axes of `cube` on a fixed (ny_out, nx_out) 

    Larger inputs are cropped, smaller ones are padded with NaN, always keeping
    the cube center fixed. Works independently per axis (one can crop while the
    other pads).
    """
    nwave, ny_in, nx_in = cube.shape
    out = np.full((nwave, ny_out, nx_out), np.nan, dtype=cube.dtype)


    def _overlap(n_in, n_out):
        # index ranges of the overlapping (shared) region, centered on both grids
        n = min(n_in, n_out)
        in_lo = (n_in - n) // 2
        out_lo = (n_out - n) // 2
        return in_lo, in_lo + n, out_lo, out_lo + n

    iy0, iy1, oy0, oy1 = _overlap(ny_in, ny_out)
    ix0, ix1, ox0, ox1 = _overlap(nx_in, nx_out)
    out[:, oy0:oy1, ox0:ox1] = cube[:, iy0:iy1, ix0:ix1]
    return out





def spatial_resample_to_kpc(cube, pixscale_arcsec, z,
                            target_kpc=2.0, output_kpc=200.0, order=1):
    

    """Resample each wavelength plane to `target_kpc` pixels, then place them on a
    fixed `output_kpc` x `output_kpc` footprint centered on the cube.

    cube            : (nwave, ny, nx) flux-per-spaxel array 
    pixscale_arcsec : native pixel size on sky [arcsec] -> 
    target_kpc      : output pixel size [kpc] -> 
    output_kpc      : output field-of-view on a side [kpc] (crop larger, NaN-pad smaller)
    Returns the resampled cube, still flux-per-spaxel, with flux conserved.
    
    """
    # Physical size of a native pixel at this redshift (proper kpc per pixel).
    kpc_per_arcsec = 1.0 / COSMO.arcsec_per_kpc_proper(z).value
    old_pix_kpc = pixscale_arcsec * kpc_per_arcsec

    # scipy's zoom factor = output_pixels / input_pixels along an axis.
    # To turn `old_pix_kpc` pixels into `target_kpc` pixels we need this factor:
    zoom_factor = old_pix_kpc / target_kpc

    # Interpolate the VALUES only (spectral axis untouched -> factor 1).
    resampled = ndimage.zoom(cube, (1.0, zoom_factor, zoom_factor), order=order)

    # Flux conservation: rescale by the pixel-AREA ratio (area_after / area_before).
    flux_scale = (target_kpc / old_pix_kpc) ** 2
    resampled *= flux_scale

    # Force a common physical footprint so every field lands on identical axes:
    # e.g. 200 kpc / 2 kpc-per-pixel = 100 pixels on a side, centered on the cube.
    n_out = int(round(output_kpc / target_kpc))
    resampled = _center_crop_or_pad(resampled, n_out, n_out)

    return resampled




def build_velocity_grid(vmin=-2500.0, vmax=2500.0, dv=25.0):
    """Common rest-frame velocity grid shared by ALL fields [km/s]."""
    return np.arange(vmin, vmax + dv, dv)


def velocity_to_obs_wavelength(v_grid, z, line_rest=oII_rest):
    """Map rest-frame velocity offsets to OBSERVED wavelengths [Angstrom]."""
    lam_line_obs = line_rest * (1.0 + z)
    return lam_line_obs * (1.0 + v_grid / c_kms)


def spectral_resample(cube, wave_native, wave_new, err_cube=None):
    """Flux-conserving rebin of every spaxel onto `wave_new` using spectres.

    spectres wants the spectral axis LAST, so we move it there and back.
    Returns (flux_new, err_new) with the spectral axis first again.

    err_cube is sigma array, not variance.  Remember to take the square root of 
    the variance before passing it in.
    """
    cube_last = np.moveaxis(cube, 0, -1)                 # (ny, nx, nwave)
    err_last = np.moveaxis(err_cube, 0, -1) if err_cube is not None else None

    out = spectres(wave_new, wave_native, cube_last,
                   spec_errs=err_last, fill=np.nan, verbose=False)

    if err_last is None:
        return np.moveaxis(out, -1, 0), None
    flux_new, err_new = out
    return np.moveaxis(flux_new, -1, 0), np.moveaxis(err_new, -1, 0)




def to_luminosity_surface_density(f_lam_cube, wave_obs_new, z,
                                  target_kpc=2.0, err_cube=None):
    """[1e-20 erg/s/cm^2/A per spaxel] -> [erg/s/(km/s)/kpc^2]."""
    DL_cm = COSMO.luminosity_distance(z).to(u.cm).value
    area_kpc2 = target_kpc ** 2

    # per-channel scale factor (depends on wavelength through lam/c)
    scale = flux_unit * 4.0 * np.pi * DL_cm**2 * (wave_obs_new / c_kms) / area_kpc2

    ell_v = f_lam_cube * scale[:, None, None]            # broadcast over (y, x)
    ell_v_err = err_cube * scale[:, None, None] if err_cube is not None else None
    return ell_v, ell_v_err





def process_field(cube, wave_native, err_cube, pixscale_arcsec, z,
                  target_kpc=2.0, dv=25.0):
    """Full per-field pipeline. Returns (ell_v cube, ell_v error cube, v_grid)."""
    # 1. spatial -> 2 kpc/pixel
    cube_s = spatial_resample_to_kpc(cube, pixscale_arcsec, z, target_kpc)
    err_s = spatial_resample_to_kpc(err_cube, pixscale_arcsec, z, target_kpc)

    # 2. spectral -> common 25 km/s grid (built once, reused for every field)
    v_grid = build_velocity_grid(dv=dv)
    wave_new = velocity_to_obs_wavelength(v_grid, z)
    f_lam, f_err = spectral_resample(cube_s, wave_native, wave_new, err_s)

    # 3. -> redshift-independent luminosity surface density per velocity
    f_v, f_v_err = to_luminosity_surface_density(f_lam, wave_new, z,
                                                     target_kpc, f_err)
    return f_v, f_v_err, v_grid


#main function to run the pipeline on a single field

def resample_main(z,cubepath, pixscale):

  


    lam_line_obs = oII_rest * (1 + z)    # observed [OII] center




    file = fits.open(cubepath)
    cube_data = file[1].data 
    cube_var = file[2].data
    cube_err = np.sqrt(cube_var)  # convert variance to sigma

    header = file[1].header  # Get the header where the coordinate system is stored

    # Extract the WCS  keywords for the spectral axis (Axis 3)
    nwave = header['NAXIS3']   # Number of channels along the wavelength axis
    crpix = header['CRPIX3']   # Reference pixel
    crval = header['CRVAL3']   # Coordinate value at the reference pixel
    cdelt = 1.25

    # get the pixel indices (0-indexed for Python)
    pixel_indices = np.arange(nwave)

    #Generate the native wavelength array
    wave_native = crval + (pixel_indices + 1 - crpix) * cdelt

    f_v, f_v_err, v_grid = process_field(
    cube_data, wave_native, cube_err, pixscale, z)

    primary_hdu = fits.PrimaryHDU(data=f_v)
    primary_hdu.name = "DATA"
    
    # Velocity WCS
    zero_idx = np.argmin(np.abs(v_grid))
    
    primary_hdu.header["CTYPE3"] = "VELO"
    primary_hdu.header["CUNIT3"] = "km/s"
    primary_hdu.header["CRPIX3"] = zero_idx + 1
    primary_hdu.header["CRVAL3"] = v_grid[zero_idx]
    primary_hdu.header["CDELT3"] = v_grid[1] - v_grid[0]
    
    primary_hdu.header["VELREF"] = oII_rest
    primary_hdu.header["VELUNIT"] = "km/s"
    
    # Other metadata
    primary_hdu.header["PIX_KPC"] = 2.0
    primary_hdu.header["DV_KMS"] = 25.0
    primary_hdu.header["Z"] = z

    # Spatial WCS -- spatial_resample_to_kpc/_center_crop_or_pad only ever
    # operate on plain numpy arrays, so without this the output cube has NO
    # spatial WCS at all and downstream tools (mpdaf, DS9, stack_cubes.py)
    # fall back to a meaningless default grid. Linear physical-offset WCS in
    # kpc, centered on the QSO (which _center_crop_or_pad keeps at the array
    # center by construction).
    ny_out, nx_out = f_v.shape[1], f_v.shape[2]
    primary_hdu.header["CTYPE1"] = "LINEAR"
    primary_hdu.header["CTYPE2"] = "LINEAR"
    primary_hdu.header["CUNIT1"] = "kpc"
    primary_hdu.header["CUNIT2"] = "kpc"
    primary_hdu.header["CDELT1"] = 2.0          # kpc/pixel, matches PIX_KPC
    primary_hdu.header["CDELT2"] = 2.0
    primary_hdu.header["CRPIX1"] = nx_out // 2 + 1   # 1-indexed FITS convention
    primary_hdu.header["CRPIX2"] = ny_out // 2 + 1
    primary_hdu.header["CRVAL1"] = 0.0          # 0 kpc = QSO position
    primary_hdu.header["CRVAL2"] = 0.0

    var_v = f_v_err ** 2

    error_hdu = fits.ImageHDU(data=var_v)
    error_hdu.name = "VAR"
    
    vgrid_hdu = fits.ImageHDU(data=v_grid)
    vgrid_hdu.name = "V_GRID"
    
    hdul = fits.HDUList([
        primary_hdu,
        error_hdu,
        vgrid_hdu
    ])
    

 
    # add headers

    hdul[0].header['PIX_KPC'] = 2.0      # Target spatial resolution
    hdul[0].header['DV_KMS'] = 25.0      # Velocity grid resolution
    hdul[0].header['Z'] = z              # Redshift of the field


    # write the HDUList to a new FITS file
    out_dir = os.path.dirname(cubepath)
    out_base = os.path.basename(cubepath).split('.fits')[0] + '_processed_cube.fits'
    output_filename = os.path.join(out_dir, out_base)
    hdul.writeto(output_filename, overwrite=True)


    hdul.close()
    return output_filename 




#example usage 

'''

z = 0.875
pixscale = 0.2   
cubepath = "/Users/charishall/CASSI_26/cubefolder/myPSFsubbed+croppedcube.fits"

f_v, f_v_err, v_grid = process_field(cube_data, wave_native, cube_err, pixscale, z)
main(z, cubepath, pixscale) 


print("output cube shape (nvel, ny, nx):", f_v.shape)





'''





# mask resampling -- same output grid as resample_main (2 kpc/pixel spatial,
# common velocity grid), but conservative: nearest-neighbor spatially and
# OR-downsampled spectrally, so a mask value never becomes a fraction and a
# flagged voxel never gets interpolated away.
#
# NOTE ON CONVENTION: on disk (both the input mask from build_mask_main and
# the output of resample_mask_main), 1/True = GOOD voxel, 0/False = BAD
# (object OR sky) -- this matches build_mask_main's explicit on-disk
# convention. The two helpers below (_spatial_resample_mask_to_kpc and
# _spectral_resample_mask) are written the opposite way, internally, for
# good reason: their padding and OR-combination logic is only conservative
# if True means BAD (pad = bad, "flag bad if any contributing channel was
# bad"). So resample_mask_main inverts the mask right after loading it and
# inverts back right before saving -- the helpers themselves should NOT be
# changed to take True=good input, since that would make the padding and
# OR-combination behave incorrectly.


def _spatial_resample_mask_to_kpc(mask, pixscale_arcsec, z,
                                   target_kpc=2.0, output_kpc=200.0):
    """Same geometry as spatial_resample_to_kpc, but nearest-neighbor
    (order=0) so mask values stay exactly 0 or 1.

    Expects/returns True = BAD (see convention note above resample_mask_main).
    """
    kpc_per_arcsec = 1.0 / COSMO.arcsec_per_kpc_proper(z).value
    old_pix_kpc = pixscale_arcsec * kpc_per_arcsec
    zoom_factor = old_pix_kpc / target_kpc

    resampled = ndimage.zoom(mask.astype(np.uint8),
                              (1.0, zoom_factor, zoom_factor), order=0)

    n_out = int(round(output_kpc / target_kpc))
    nwave, ny_in, nx_in = resampled.shape
    out = np.ones((nwave, n_out, n_out), dtype=np.uint8)  # pad = 1 (bad)

    def _overlap(n_in, n_out_):
        n = min(n_in, n_out_)
        in_lo = (n_in - n) // 2
        out_lo = (n_out_ - n) // 2
        return in_lo, in_lo + n, out_lo, out_lo + n

    iy0, iy1, oy0, oy1 = _overlap(ny_in, n_out)
    ix0, ix1, ox0, ox1 = _overlap(nx_in, n_out)
    out[:, oy0:oy1, ox0:ox1] = resampled[:, iy0:iy1, ix0:ix1]
    return out.astype(bool)


def _spectral_resample_mask(mask, wave_native, wave_new):
    """OR-downsample the mask onto wave_new: a new bin is flagged bad if
    ANY native channel falling in that bin was bad (conservative -- never
    silently un-flags a voxel through interpolation, unlike a flux-style
    rebin would).

    Expects/returns True = BAD (see convention note above resample_mask_main).
    """
    edges = np.concatenate((
        [wave_new[0] - (wave_new[1] - wave_new[0]) / 2],
        (wave_new[:-1] + wave_new[1:]) / 2,
        [wave_new[-1] + (wave_new[-1] - wave_new[-2]) / 2],
    ))
    mask_new = np.zeros((len(wave_new),) + mask.shape[1:], dtype=bool)
    for i in range(len(wave_new)):
        sel = (wave_native >= edges[i]) & (wave_native < edges[i + 1])
        if sel.any():
            mask_new[i] = mask[sel].any(axis=0)
        else:
            idx = np.argmin(np.abs(wave_native - wave_new[i]))
            mask_new[i] = mask[idx]
    return mask_new


def resample_mask_main(z, maskpath, pixscale, dv=25.0, target_kpc=2.0,
                        output_kpc=200.0):
    """Resample an (already cropped) mask onto the same 2 kpc/pixel spatial
    grid and common velocity grid that resample_main puts the science cube
    on, so the two line up voxel-for-voxel. Saves '<mask>_processed_mask.fits'.

    On disk (both input and output): 1 = good voxel, 0 = bad (object OR sky),
    matching build_mask_main's convention. The resampling helpers used here
    internally assume True = bad (so that edge-padding and OR-combination
    are conservative), so the mask is inverted on the way in and inverted
    back on the way out -- the helpers themselves are untouched.
    """
    with fits.open(maskpath) as hdul:
        mask = hdul[0].data.astype(bool)   # on-disk: True = good
        header = hdul[0].header

        nwave = header['NAXIS3']
        crpix = header['CRPIX3']
        crval = header['CRVAL3']
        cdelt = 1.25  # matches resample_main's hardcoded native step

        pixel_indices = np.arange(nwave)
        wave_native = crval + (pixel_indices + 1 - crpix) * cdelt

    # flip to the internal True=bad convention the helpers below expect
    bad = ~mask

    bad_s = _spatial_resample_mask_to_kpc(bad, pixscale, z, target_kpc, output_kpc)

    v_grid = build_velocity_grid(dv=dv)
    wave_new = velocity_to_obs_wavelength(v_grid, z)
    bad_v = _spectral_resample_mask(bad_s, wave_native, wave_new)

    # flip back to the on-disk True=good convention before saving
    mask_v = ~bad_v

    primary_hdu = fits.PrimaryHDU(data=mask_v.astype(np.uint8))
    primary_hdu.name = "MASK"
    primary_hdu.header['COMMENT'] = '1 = good voxel, 0 = bad (object OR sky)'
    primary_hdu.header['PIX_KPC'] = target_kpc
    primary_hdu.header['DV_KMS'] = dv
    primary_hdu.header['Z'] = z

    vgrid_hdu = fits.ImageHDU(data=v_grid)
    vgrid_hdu.name = "V_GRID"

    hdul_out = fits.HDUList([primary_hdu, vgrid_hdu])
    out_dir = os.path.dirname(maskpath)
    out_base = os.path.basename(maskpath).split('.fits')[0] + '_processed_mask.fits'
    output_filename = os.path.join(out_dir, out_base)
    hdul_out.writeto(output_filename, overwrite=True)

    return output_filename



# example usage
'''
mask_final = resample_mask_main(z, mask_crop_path, pixscale)
'''
