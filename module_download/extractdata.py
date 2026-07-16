"""Look up cube paths for a QSO and convert its sky coordinates to pixel
coordinates in a given cube."""

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table


def get_cube_paths(name, pairs_table, cube_dir):
    """
    Look up all cube files associated with a given QSO name, using the
    pairs table. Adjust the column names below to match your actual
    pairs table schema.
    """
    rows = pairs_table[pairs_table["Name"] == name]
    return [Path(cube_dir) / fname for fname in rows["cube_filename"]]


def get_z(path="path/to/sample.fits"):
    """Redshift column from the sample table."""
    t = Table.read(path)
    return np.asarray(t["z"])


def convert_coords(ra, dec, cubepath):
    """Convert a QSO's (ra, dec) [deg] to pixel (x, y) in `cubepath`'s WCS.

    NB: this uses the RA/Dec you pass in -- it does NOT use the cube's own
    WCS reference position. Using the WCS's own CRVAL here (a previous bug)
    would return the same pixel for every cube regardless of the QSO's
    actual position, which silently breaks all downstream cropping/PSF
    centering.
    """
    with fits.open(cubepath) as hdul:
        wcs3d = WCS(hdul[1].header)
        wcs = wcs3d.celestial  # drop the AWAV axis, keep RA/Dec only

    coord = SkyCoord(ra, dec, unit="deg")
    x, y = wcs.world_to_pixel(coord)
    return float(x), float(y)


# example usage:
#   z = get_z('path/to/sample.fits')
#   x, y = convert_coords(ra, dec, 'path/to/cube.fits')
