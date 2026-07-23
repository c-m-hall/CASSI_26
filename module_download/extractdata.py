"""
Bridges the local pairs table written by download_sample.py (Name,
cube_filename) to mainscript.py, which needs real file paths and QSO pixel
coordinates per source.
"""

import os

import numpy as np
from astropy.wcs import WCS


def get_cube_paths(name, pairs, cube_dir):
    """All local cube file paths for a given QSO `name`.

    `pairs` is the table written by download_sample.py (columns: Name,
    cube_filename). Returns a list of full paths under `cube_dir` -- empty
    if this source has no downloaded cubes.
    """
    matches = pairs[np.asarray(pairs['Name']).astype(str) == str(name)]
    return [os.path.join(cube_dir, str(fn)) for fn in matches['cube_filename']]


def get_z(name, sample):
    """Redshift for `name` looked up from the sample table (Name, z, ...).

    Raises if the name isn't found, since a missing redshift would silently
    break every downstream wavelength/velocity calculation.
    """
    matches = sample[np.asarray(sample['Name']).astype(str) == str(name)]
    if len(matches) == 0:
        raise KeyError(f"'{name}' not found in sample table")
    return float(matches['z'][0])


def convert_coords(ra, dec, cubepath):
    """RA/Dec [deg] -> (x, y) pixel coordinates in `cubepath`'s spatial WCS.

    Reads the celestial WCS from extension 1 (the DATA extension), matching
    what crop_cube/psf_sub_main expect as x_center/y_center.
    """
    from astropy.io import fits

    with fits.open(cubepath) as hdul:
        header = hdul[1].header

    # naxis=2 drops the spectral axis, leaving just the RA/Dec pair.
    wcs = WCS(header, naxis=2)
    x, y = wcs.all_world2pix(ra, dec, 0)
    return float(x), float(y)
