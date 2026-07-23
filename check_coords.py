from astropy.table import Table
from module_download.extractdata import get_cube_paths, convert_coords

sample = Table.read('muse_x_milliquas_sample.fits')
pairs = Table.read('muse_x_milliquas_pairs_local.fits')

for row in sample:
    name, ra, dec = row['Name'], row['RA'], row['DEC']
    for cubepath in get_cube_paths(name, pairs, './scratch/'):
        x, y = convert_coords(ra, dec, cubepath)
        print(f"{name} / {cubepath}: RA={ra:.6f} DEC={dec:.6f} -> x={x:.2f} y={y:.2f}")
