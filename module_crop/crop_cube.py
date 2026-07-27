# Crop a box centered on (x_center, y_center), in both spatial and spectral axes.

from astropy.io import fits


def crop_cube(x_center, y_center, cubepath, spatialcrop_pix, lam_obs, dwave):
    sp_crop_size = spatialcrop_pix / 2

    with fits.open(cubepath) as hdul:
        data = hdul[1].data
        var_data = hdul[2].data
        header = hdul[1].header

        y_min = int(max(0, y_center - sp_crop_size))
        y_max = int(min(data.shape[1], y_center + sp_crop_size))
        x_min = int(max(0, x_center - sp_crop_size))
        x_max = int(min(data.shape[2], x_center + sp_crop_size))

        sp_cropped_cube = data[:, y_min:y_max, x_min:x_max]
        var_cropped_cube = var_data[:, y_min:y_max, x_min:x_max]

        # WCS spectral coordinates -- read CDELT3 from the header rather than
        # hardcoding it, since not every cube is guaranteed to be resampled
        # to the same spectral step at this stage.
        crval3 = header["CRVAL3"]
        cdelt3 = header.get("CD3_3", header.get("CDELT3"))
        if cdelt3 is None:
            raise KeyError(
                f"{cubepath}: neither CD3_3 nor CDELT3 found in header; "
                "cannot determine spectral pixel scale."
            )
        crpix3 = header["CRPIX3"]  # 1-indexed, per FITS convention

        # central wavelength (lam_obs, Angstrom) -> 0-indexed pixel channel
        # pixel = ((value - CRVAL) / CDELT) + (CRPIX - 1)
        lambda_0_pixel = ((lam_obs - crval3) / cdelt3) + (crpix3 - 1)
        deltalambda_pixel = abs(dwave / cdelt3)

        num_wave = sp_cropped_cube.shape[0]
        w_min = max(0, int(round(lambda_0_pixel - deltalambda_pixel / 2)))
        w_max = min(num_wave, int(round(lambda_0_pixel + deltalambda_pixel / 2)))

        final_cube = sp_cropped_cube[w_min:w_max, :, :]
        var_new = var_cropped_cube[w_min:w_max, :, :]

        hdul[1].data = final_cube
        hdul[2].data = var_new

        # shift the spectral reference pixel to match the new cropped start
        header["CRPIX3"] = crpix3 - w_min
        hdul[2].header["CRPIX3"] = crpix3 - w_min

        # --- fix: also shift the SPATIAL reference pixels (CRPIX1/CRPIX2).
        # Without this, the WCS in the cropped file no longer points at the
        # QSO, and anything downstream that re-derives x,y from RA/Dec on
        # this file (or trusts the header WCS) will be wrong.
        for hdu_data_header in (header, hdul[2].header):
            if "CRPIX1" in hdu_data_header:
                hdu_data_header["CRPIX1"] -= x_min
            if "CRPIX2" in hdu_data_header:
                hdu_data_header["CRPIX2"] -= y_min

        output_filename = cubepath.replace(".fits", "_fully_cropped.fits")
        hdul.writeto(output_filename, overwrite=True)

    return output_filename


# example usage
"""
import astropy.io.fits as fits
from astropy.cosmology import Planck18 as COSMO

x_center = 194.22
y_center = 180.8
cubepath = "/Users/charishall/CASSI_26/J0015/J0015_vac_SUBBED.fits"
lam_obs = oII_rest * (1 + z)    # observed [OII] center, ~7945 A

spatialcrop_pix = 200
pixscale = 0.2
vel_kms = 5000
c_kms = 299792.458
dwave = (2 * vel_kms / c_kms) * lam_obs

cropped_cube = crop_cube(x_center, y_center, cubepath, spatialcrop_pix, lam_obs, dwave)
"""




def crop_mask(x_center, y_center, maskpath, spatialcrop_pix, lam_obs, dwave):
    """Crop a 3D mask (single-extension, 0/1, no variance) to the exact same
    spatial box + wavelength window that crop_cube uses on the science cube.
    
    Call with the same x_center/y_center/spatialcrop_pix/lam_obs/dwave passed to crop_cube on the corresponding *_MASK3D.fits file, so the crop
    lines up with the cropped science cube.
    """
    sp_crop_size = spatialcrop_pix / 2

    with fits.open(maskpath) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        y_min = int(max(0, y_center - sp_crop_size))
        y_max = int(min(data.shape[1], y_center + sp_crop_size))
        x_min = int(max(0, x_center - sp_crop_size))
        x_max = int(min(data.shape[2], x_center + sp_crop_size))

        sp_cropped = data[:, y_min:y_max, x_min:x_max]

        crval3 = header["CRVAL3"]
        cdelt3 = header.get("CD3_3", header.get("CDELT3"))
        if cdelt3 is None:
            raise KeyError(
                f"{maskpath}: neither CD3_3 nor CDELT3 found in header; "
                "cannot determine spectral pixel scale."
            )
        crpix3 = header["CRPIX3"]

        lambda_0_pixel = ((lam_obs - crval3) / cdelt3) + (crpix3 - 1)
        deltalambda_pixel = abs(dwave / cdelt3)

        num_wave = sp_cropped.shape[0]
        w_min = max(0, int(round(lambda_0_pixel - deltalambda_pixel / 2)))
        w_max = min(num_wave, int(round(lambda_0_pixel + deltalambda_pixel / 2)))

        final_mask = sp_cropped[w_min:w_max, :, :]

        hdul[0].data = final_mask
        header["CRPIX3"] = crpix3 - w_min
        if "CRPIX1" in header:
            header["CRPIX1"] -= x_min
        if "CRPIX2" in header:
            header["CRPIX2"] -= y_min

        output_filename = maskpath.replace(".fits", "_fully_cropped.fits")
        hdul.writeto(output_filename, overwrite=True)

    return output_filename


# example usage
"""
mask_cropped = crop_mask(x_center, y_center, maskpath, spatialcrop_pix, lam_obs, dwave)
"""

