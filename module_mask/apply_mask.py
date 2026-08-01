import os
import numpy as np
from astropy.io import fits


def apply_mask_main(cube_path, mask_path, out_dir=None):
    """
    Multiply a cube by its mask and set masked-out (0) voxels to NaN,
    preparing it for stacking.

    both cube and mask datanare stored in extension 0 
    Mask should already be on the same spatial/spectral grid as the cube (post-resample)
    """
    with fits.open(cube_path) as cube_hdul, fits.open(mask_path) as mask_hdul:
        cube_hdu =  cube_hdul[0]
        mask_hdu = mask_hdul[0]

        cube_data = cube_hdu.data
        mask_data = mask_hdu.data

        if cube_data.shape != mask_data.shape:
            raise ValueError(
                f"Shape mismatch: cube {cube_data.shape} vs mask {mask_data.shape} "
                f"for {os.path.basename(cube_path)}"
            )

        # multiply, then zero -> NaN
        masked_data = cube_data * mask_data
        masked_data[mask_data == 0] = np.nan

        cube_hdu.data = masked_data.astype(cube_data.dtype, copy=False)

        out_dir = out_dir or os.path.dirname(cube_path)
        os.makedirs(out_dir, exist_ok=True)
        base = os.path.basename(cube_path).replace('.fits', '_masked.fits')
        out_path = os.path.join(out_dir, base)

        cube_hdul.writeto(out_path, overwrite=True)

    return out_path

#example usage
'''
masked_final_path = apply_mask_main('/Users/charishall/CASSI_26/resampled/ADP.2016-07-12T10.19.53.321_vac_PSFSUBBED_fully_cropped_processed_cube.fits', '/Users/charishall/CASSI_26/resampled/ADP.2016-06-20T07.31.00.321_vac_MASK3D_fully_cropped_processed_mask.fits', out_dir= '/Users/charishall/CASSI_26/resampled/')
'''       