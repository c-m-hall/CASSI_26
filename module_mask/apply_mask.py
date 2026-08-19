import os
import numpy as np
from astropy.io import fits


def apply_mask_main(cube_path, mask_path, out_dir=None):
    """
    Multiply a cube by its mask and set masked-out (0) voxels to NaN,
    preparing it for stacking.
    Cube: DATA in extension 0, VAR in extension 1 (both masked identically
    so their NaN footprints match for stack_cubes.py's n/n^2 normalization).
    Mask data is stored in extension 0 (1 = good, 0 = bad).
    Mask should already be on the same spatial/spectral grid as the cube (post-resample)
    """
    with fits.open(cube_path) as cube_hdul, fits.open(mask_path) as mask_hdul:
        cube_hdu = cube_hdul[0]
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

        # apply the identical mask to VAR (extension 1) so NaN footprints
        # in DATA and VAR match -- stack_cubes.py's n/n^2 normalization
        # relies on this
        if len(cube_hdul) > 1 and cube_hdul[1].name == 'VAR':
            var_hdu = cube_hdul[1]
            var_data = var_hdu.data
            if var_data.shape != mask_data.shape:
                raise ValueError(
                    f"Shape mismatch: VAR {var_data.shape} vs mask {mask_data.shape} "
                    f"for {os.path.basename(cube_path)}"
                )
            masked_var = var_data * mask_data
            masked_var[mask_data == 0] = np.nan
            var_hdu.data = masked_var.astype(var_data.dtype, copy=False)

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
