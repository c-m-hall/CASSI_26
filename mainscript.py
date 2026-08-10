
import os
import shutil
 
import numpy as np
from astropy.io import fits
from astropy.table import Table
from pathlib import Path
 
#import all modules scripts
from module_download.extractdata import get_cube_paths, get_z, convert_coords
from module_crop.crop_cube import crop_cube
from module_resample.resample import resample_main
from module_convert_wl.convert_wavel import convert_wl_main
from module_psf_sub.PSFSubtraction import psf_sub_main
 
from module_mask.build_mask import build_mask_main
from module_mask.apply_mask import apply_mask_main
from module_crop.crop_cube import crop_cube, crop_mask
from module_resample.resample import resample_main, resample_mask_main
 
C_KMS = 299792.458
OII_REST = 3728.48
 
BASE_DIR = '/carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/'
 
 
# output subdirectories for each pipeline stage
DIR_WLCONV = os.path.join(BASE_DIR, 'converted_wl') + '/'
 
DIR_PSFSUB = os.path.join(BASE_DIR, 'psf_subbed') + '/'
 
DIR_CROP   = os.path.join(BASE_DIR, 'cropped') + '/'
 
DIR_RESAMPLED  = os.path.join(BASE_DIR, 'resampled') + '/'
 
DIR_FINAL  = os.path.join(BASE_DIR, 'masked') + '/'
 
for d in (DIR_WLCONV, DIR_PSFSUB, DIR_CROP, DIR_RESAMPLED, DIR_FINAL):
    os.makedirs(d, exist_ok=True)
 
def move_to(path, dest_dir):
    """Relocate step's output into its staged subdirectory; return new path."""
    os.makedirs(dest_dir, exist_ok=True)
    new_path = os.path.join(dest_dir, os.path.basename(path))
    shutil.move(path, new_path)
    return new_path
 
 
def set_exptime_header(path, exptime):
    """Force-write EXPTIME into a FITS file's primary header in place.
 
    Called after every stage that writes a new file, since none of the
    module functions (convert_wl_main, psf_sub_main, crop_cube,
    resample_main, apply_mask_main) are guaranteed to carry the keyword
    forward on their own.
    """
    if exptime is None:
        return
    with fits.open(path, mode='update') as hdul:
        hdul[0].header['EXPTIME'] = exptime
        hdul.flush()
 
 
def get_exptime_for_cube(pairs, cube_filename):
    """EXPTIME [s] for this cube from the local pairs table (Name,
    cube_filename, exptime), written by download_sample.py."""
    match = pairs[np.asarray(pairs['cube_filename']).astype(str) == cube_filename]
    if len(match) == 0:
        print(f"  [WARNING] No pairs-table match for {cube_filename}, EXPTIME will be missing")
        return None
    return float(match['exptime'][0])
 
 
def main(sample_path='muse_x_milliquas_sample.fits',
         pairs_path='muse_x_milliquas_pairs_local.fits',
         cube_dir='./scratch/', spatialcrop_pix=200, vel_window_kms=5000, pixscale=0.2,
         sex_config='/Users/charishall/CASSI_26/SExtractor_Files/wl_eso.sex', sex_binary='/opt/homebrew/bin/sex', resume_from='download'):
 
    sample = Table.read(sample_path)
    pairs = Table.read(pairs_path)
 
    for row in sample:
        name, ra, dec, z = row['Name'], row['RA'], row['DEC'], row['z']
 
        cube_paths = get_cube_paths(name, pairs, cube_dir)
        if not cube_paths:
            print(f"No cubes found for {name}, skipping.")
            continue
 
        for cubepath in cube_paths:
            cube_base = os.path.basename(cubepath)
            try:
                x, y = convert_coords(ra, dec, cubepath)
                print(f"{name} / {cube_base}: RA={ra:.6f}, DEC={dec:.6f} -> x={x:.2f}, y={y:.2f}")
 
                exptime = get_exptime_for_cube(pairs, cube_base)
 
                if resume_from == 'mask':
                    # skip re-running convert_wl / psf_sub; find their
                    # already-produced outputs on disk instead
                    vac_name = cube_base.replace('.fits', '_vac.fits')
                    wlconv_path = os.path.join(DIR_WLCONV, vac_name)
                    psfsub_path = os.path.join(
                        DIR_PSFSUB, vac_name.replace('.fits', '_PSFSUBBED.fits'))
                    if not (os.path.exists(wlconv_path) and os.path.exists(psfsub_path)):
                        raise FileNotFoundError(
                            f"resume_from='mask' but missing prior output(s): "
                            f"{wlconv_path if not os.path.exists(wlconv_path) else psfsub_path}"
                        )
                    print(f"  resuming from mask stage using {wlconv_path}, {psfsub_path}")
                else:
                    # 1. air-to-vacuum wavelength correction
                    cube_dirname = os.path.dirname(cubepath) + '/'
                    wlconv_path = convert_wl_main(cube_base, cube_dirname)
                    wlconv_path = move_to(wlconv_path, DIR_WLCONV)
                    set_exptime_header(wlconv_path, exptime)
 
                    # 2. PSF subtraction, on the wavelength-corrected cube
                    psfsub_path = psf_sub_main(wlconv_path, x, y, z)
                    psfsub_path = move_to(psfsub_path, DIR_PSFSUB)
                    set_exptime_header(psfsub_path, exptime)
 
                # 2b. spatial (SExtractor) + spectral (sky-residual) mask,
                # combined into one 3D mask and saved alongside the
                # vacuum-wavelength cube -- make spectral mask  BEFORE PSF subtraction, and then spatial part AFTER PSF subtraction so the
                # white-light image SExtractor sees, and the residual sky
                # spectrum used for the spectral mask, are both unsubtracted but the spatial mask does not include original LARGE PSF
                mask_path = build_mask_main(wlconv_path, psfcube=psfsub_path, sex_config=sex_config,
                                             sex_binary=sex_binary)
 
                # 3. crop around the QSO, centered on observed [OII]
                lam_obs = OII_REST * (1 + z)
                dwave = (2 * vel_window_kms / C_KMS) * lam_obs
                crop_path = crop_cube(x, y, psfsub_path, spatialcrop_pix, lam_obs, dwave)
                crop_path = move_to(crop_path, DIR_CROP)
                set_exptime_header(crop_path, exptime)
 
                # 3b. crop the mask to the identical spatial/spectral window
                mask_crop_path = crop_mask(x, y, mask_path, spatialcrop_pix, lam_obs, dwave)
                mask_crop_path = move_to(mask_crop_path, DIR_CROP)
 
                # 4. spatial/spectral resample -> final product
                final_path = resample_main(z, crop_path, pixscale)
                final_path = move_to(final_path, DIR_RESAMPLED)
                set_exptime_header(final_path, exptime)
 
                # 4b. resample the mask onto the identical output grid
                mask_final_path = resample_mask_main(z, mask_crop_path, pixscale)
                mask_final_path = move_to(mask_final_path, DIR_RESAMPLED)
 
                #5. apply the mask to the cube
                masked_cube_final_path = apply_mask_main(final_path, mask_final_path, out_dir=DIR_FINAL)
                set_exptime_header(masked_cube_final_path, exptime)
 
                print(f"Done: {name} / {cube_base} -> {masked_cube_final_path}")
 
            except Exception as e:
                print(f"Failed on {name} / {cube_base}: {e}")
 
 
if __name__ == '__main__':
    main()
 
