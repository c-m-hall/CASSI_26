#with compiled sample as 'sample.fits'

import os
import shutil

from astropy.table import Table

#import all modules scripts
from module_download.extractdata import get_cube_paths, get_z, convert_coords
from module_crop.crop_cube import crop_cube
from module_resample.resample import resample_main
from module_convert_wl.convert_wavel import convert_wl_main
from module_psf_sub.PSFSubtraction import psf_sub_main

C_KMS = 299792.458
OII_REST = 3728.48

# output subdirectories for each pipeline stage
DIR_WLCONV = '/Users/charishall/CASSI_26/converted_wl/'
DIR_PSFSUB = '/Users/charishall/CASSI_26/psf_subbed/'
DIR_CROP   = '/Users/charishall/CASSI_26/cropped/'
DIR_FINAL  = '/Users/charishall/CASSI_26/resampled/'


def move_to(path, dest_dir):
    """Relocate step's output into its staged subdirectory; return new path."""
    os.makedirs(dest_dir, exist_ok=True)
    new_path = os.path.join(dest_dir, os.path.basename(path))
    shutil.move(path, new_path)
    return new_path

def main(sample_path='muse_x_milliquas_sample.fits',
         pairs_path='muse_x_milliquas_pairs_local.fits',
         cube_dir='./scratch/', spatialcrop_pix=200, vel_window_kms=5000, pixscale=0.2,
         sex_config='/Users/charishall/CASSI_26/SExtractor_Files/wl_eso.sex', sex_binary='sex'):

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

                # 1. air-to-vacuum wavelength correction
                cube_dirname = os.path.dirname(cubepath) + '/'
                wlconv_path = convert_wl_main(cube_base, cube_dirname)
                wlconv_path = move_to(wlconv_path, DIR_WLCONV)

                # 1b. spatial (SExtractor) + spectral (sky-residual) mask,
                # combined into one 3D mask and saved alongside the
                # vacuum-wavelength cube -- BEFORE PSF subtraction, so the
                # white-light image SExtractor sees, and the residual sky
                # spectrum used for the spectral mask, are both unsubtracted.
                mask_path = build_mask_main(wlconv_path, sex_config=sex_config,
                                             sex_binary=sex_binary)
                print(f"  mask -> {mask_path}")

                # 2. PSF subtraction, on the wavelength-corrected cube
                psfsub_path = psf_sub_main(wlconv_path, x, y, z)
                psfsub_path = move_to(psfsub_path, DIR_PSFSUB)

                # 3. crop around the QSO, centered on observed [OII]
                lam_obs = OII_REST * (1 + z)
                dwave = (2 * vel_window_kms / C_KMS) * lam_obs
                crop_path = crop_cube(x, y, psfsub_path, spatialcrop_pix, lam_obs, dwave)
                crop_path = move_to(crop_path, DIR_CROP)

                # 3b. crop the mask to the identical spatial/spectral window
                mask_crop_path = crop_mask(x, y, mask_path, spatialcrop_pix, lam_obs, dwave)
                mask_crop_path = move_to(mask_crop_path, DIR_CROP)

                # 4. spatial/spectral resample -> final product
                final_path = resample_main(z, crop_path, pixscale)
                final_path = move_to(final_path, DIR_FINAL)

                # 4b. resample the mask onto the identical output grid
                mask_final_path = resample_mask_main(z, mask_crop_path, pixscale)
                mask_final_path = move_to(mask_final_path, DIR_FINAL)

                print(f"Done: {name} / {cube_base} -> {final_path} (mask: {mask_final_path})")

            except Exception as e:
                print(f"Failed on {name} / {cube_base}: {e}")


if __name__ == '__main__':
    main()


