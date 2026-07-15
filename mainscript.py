
#with compiled sample as 'sample.fits'


from astropy.table import Table

#import all modules scripts

from module_download.downloadscript import get_cube_paths, get_z, convert_coords
from module_crop.crop_cube import crop_cube
from module_resample.resample import resample_main
from module_convert_wl.convert_wavel import convert_wl_main
from module_psf_sub.PSFSubtraction import psf_sub_main



OII_REST = 3728.48

# output subdirectories for each pipeline stage
DIR_WLCONV = '/Users/charishall/CASSI_26/converted_wl/'
DIR_PSFSUB = '/Users/charishall/CASSI_26/psf_subbed/'
DIR_CROP   = '/Users/charishall/CASSI_26/cropped/'
DIR_FINAL  = '/Users/charishall/CASSI_26/resampled/'




def move_to(path, dest_dir):
    """Relocate step's output into its staged subdirectory; return new path."""
    new_path = os.path.join(dest_dir, os.path.basename(path))
    shutil.move(path, new_path)
    return new_path


def main(sample_path='sample.fits', pairs_path='muse_x_milliquas_pairs.fits',
         cube_dir='.', spatialcrop_pix=200, vel_window_kms=5000, pixscale=0.2):

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

                # 1. air-to-vacuum wavelength correction
                cube_dirname = os.path.dirname(cubepath) + '/'
                wlconv_path = convert_wl_main(cube_base, cube_dirname)
                wlconv_path = move_to(wlconv_path, DIR_WLCONV)

                # 2. PSF subtraction, on the wavelength-corrected cube
                psfsub_path = psf_sub_main(wlconv_path, x, y, z)
                psfsub_path = move_to(psfsub_path, DIR_PSFSUB)

                # 3. crop around the QSO, centered on observed [OII]
                lam_obs = OII_REST * (1 + z)
                dwave = (2 * vel_window_kms / C_KMS) * lam_obs
                crop_path = crop_cube(x, y, psfsub_path, spatialcrop_pix, lam_obs, dwave)
                crop_path = move_to(crop_path, DIR_CROP)

                # 4. spatial/spectral resample -> final product
                final_path = resample_main(z, crop_path, pixscale)
                final_path = move_to(final_path, DIR_FINAL)

                print(f"Done: {name} / {cube_base} -> {final_path}")

            except Exception as e:
                print(f"Failed on {name} / {cube_base}: {e}")


if __name__ == '__main__':
    main()
