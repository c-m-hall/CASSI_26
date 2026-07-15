#download data from catalog and extract values from fits file


from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from astropy.table import Table
import numpy as np




def get_z(path ='path/to/sample.fits file'):
    with fits.open(path) as hdul:
        header = hdul[1].header

    z = np.asarray(header['z'])
    return(z)

#example usage: z = get_z('path/to/sample.fits file')

        

def convert_coords(path=cubepath)
    
    # Load WCS
    with fits.open('cubepath') as hdul:
        wcs3d = WCS(hdul[1].header)
        wcs = wcs3d.celestial   # drop the AWAV axis, keep RA/Dec only

    # Reference RA/Dec pulled from the WCS itself
    ra_ref, dec_ref = wcs.wcs.crval
    coord = SkyCoord(ra_ref, dec_ref, unit='deg')

    x, y = wcs.world_to_pixel(coord)
    sky = wcs.pixel_to_world(x, y)

    return(x,y)

#example usage: ra, dec = convert_coords('path/to/sample.fits file')