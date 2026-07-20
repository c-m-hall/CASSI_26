from pathlib import Path

import numpy as np
import pyvo
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time


# Match radius: half of MUSE's wide-field-mode 1'x1' .... An object within this distance of a pointing center is guaranteed to land inside
# that cube 

MATCH_RADIUS = 30 * u.arcsec  # half MUSE WFM's 1'x1' FoV


# Local cache files -- avoids re-hitting the archives on every run.
MILLIQUAS_CACHE = Path('milliquas_v8.fits')
MUSE_CACHE = Path('eso_muse_cubes.fits')



def cached_tap_query(path, tap_url, adql, maxrec=2_000_000):
    """Run a TAP query and cache the result as FITS; reuse the cache if present."""
    if path.exists():
        # If we've already run this exact query before, just load the saved copy.
        return Table.read(path)
    
    # Otherwise hit the remote TAP (Table Access Protocol) service and run
    # the ADQL query string against it.
    result = pyvo.dal.TAPService(tap_url).search(adql, maxrec=maxrec).to_table()


    for col in result.colnames:
        # variable-width string columns from the VOTable can't be written to FITS
        if result[col].dtype.kind == 'O':
            result[col] = result[col].astype(str)
    result.write(path)
    return result








# 1. Milliquas catalogue 

# get name, coordinates, redshift, redshift source (rz), classification, and R/B magnitudes.

milliquas = cached_tap_query(
    MILLIQUAS_CACHE,
    'https://tapvizier.cds.unistra.fr/TAPVizieR/tap',
    'SELECT "Name", "RAJ2000", "DEJ2000", "z", "rz", "Type", "Rmag", "Bmag" FROM "VII/294/catalog"',
)
print(f'{len(milliquas)} Milliquas objects')



# 2. ESO MUSE datacubes 
#Filtering to cube-type products from MUSE only..... obs_creator_did is needed later to tell individual per-exposure cubes apart from SKY-offset and NFM-mode cubes.

muse = cached_tap_query(
    MUSE_CACHE,
    'https://archive.eso.org/tap_obs',
    '''SELECT target_name, s_ra, s_dec, t_exptime, t_min, obs_collection,
              dp_id, obs_creator_did
       FROM ivoa.ObsCore
       WHERE dataproduct_type = 'cube' AND instrument_name = 'MUSE' ''',
)
print(f'{len(muse)} MUSE cubes')




# 3. Cross match 
# Build SkyCoord arrays for both catalogs, then for every MUSE pointing, find every Milliquas object within MATCH_RADIUS
# Note a "match" is a pair -- the same object can match multiple cubes (repeat pointings of the same field).

mq_coord = SkyCoord(milliquas['RAJ2000'], milliquas['DEJ2000'], unit=u.deg)
muse_coord = SkyCoord(muse['s_ra'], muse['s_dec'], unit=u.deg)
 
idx_mq, idx_muse, sep, _ = muse_coord.search_around_sky(mq_coord, MATCH_RADIUS)
print(f'{len(idx_mq)} object-cube pairs: '
      f'{len(np.unique(idx_mq))} unique objects, '
      f'{len(np.unique(idx_muse))} unique cubes')





 
# 4. Build the pair table 
# One row per (object, cube) match. Stitch together the matched columns from both catalogs using the index arrays from the cross-match above.

pairs = Table()

for col in ('Name', 'RAJ2000', 'DEJ2000', 'z', 'rz', 'Type', 'Rmag', 'Bmag'):
    pairs[col] = milliquas[col][idx_mq]

for col in ('target_name', 'obs_collection', 't_exptime', 'dp_id'):
    pairs[col] = muse[col][idx_muse]
pairs['date_obs'] = Time(muse['t_min'][idx_muse], format='mjd').isot



# ORIGFILE: obs_creator_did is a URI ending in the original filename, e.g.
# ..._WFM-NOAO-N_OBJ.fits. Splitting on '?' and taking the last piece extracts just that filename, which encodes the observing mode 
# (WFM-NOAO-N, WFM-AO-E, NFM-AO-N, ...) -- used below to filter out NFM
# and SKY-offset cubes. Combined/deep products carry no mode string at all.



pairs['origfile'] = [d.split('?')[-1] for d in muse['obs_creator_did'][idx_muse].astype(str)]
pairs['sep'] = sep.arcsec.round(2) * u.arcsec
 

# Milliquas' z column can come back masked (missing values); fill those
# with NaN 
if hasattr(pairs['z'], 'filled'):
    pairs['z'] = pairs['z'].filled(np.nan)
 
# Photometric-redshift flag to determine if z is phot or spec 
#photometric redshift is a lot less accurate than spectroscopic redshift, and for our purpose of stacking an emission line, this uncertainty matters
#we look for sources only with spectroscopic redshift, and also kept a record for where it comes from in case in the future we want a reference for robustness check
z10 = np.asarray(pairs['z'], float) * 10
rounded_to_0p1 = np.abs(z10 - np.round(z10)) < 1e-6
cites_gaia = np.char.find(np.char.upper(np.asarray(pairs['rz']).astype(str)), 'GAIA') >= 0
pairs['z_phot'] = rounded_to_0p1 | cites_gaia
 
# individual cubes = obs_collection 'MUSE', minus SKY-offset exposures and NFM-mode cubes

of = np.asarray(pairs['origfile']).astype(str)
col = np.asarray(pairs['obs_collection']).astype(str)
is_indiv = (col == 'MUSE') & (np.char.find(of, '_SKY') < 0) & (np.char.find(of, 'NFM') < 0)
pairs_indiv = pairs[is_indiv]
print(f'{len(pairs)} total pairs -> {len(pairs_indiv)} individual-cube pairs '
      f'after dropping combined products and SKY/NFM cubes')




# 5. One row per object


# NOTE: this groups the FULL `pairs` table, not `pairs_indiv`. That means
# n_cubes/total_exptime/dp_ids below include combined products, SKY-offset
# exposures, and NFM cubes 


def join_unique(col):
    return ';'.join(np.unique(np.asarray(col).astype(str)))

rows = []
for g in pairs.group_by('Name').groups:
    rows.append((
        g['Name'][0], g['RAJ2000'][0], g['DEJ2000'][0], g['z'][0], g['z_phot'][0],
        g['Type'][0], join_unique(g['target_name']), len(g),
        float(g['t_exptime'].sum()), float(g['sep'].min()), join_unique(g['dp_id']),
    ))

fields = Table(
    rows=rows,
    names=('Name', 'RA', 'DEC', 'z', 'z_phot', 'Type', 'muse_names',
           'n_cubes', 'total_exptime', 'min_sep', 'dp_ids'),
    units={'RA': u.deg, 'DEC': u.deg, 'total_exptime': u.s, 'min_sep': u.arcsec},
)
fields.sort('total_exptime', reverse=True)
print(f'{len(fields)} matched objects')



#   science cut
# Type-1 QSOs, 0.2 < z < 1.5, more than 15 min of exposure, and spec redshift (not phot). 

is_qso = np.array([t.startswith('Q') for t in fields['Type'].astype(str)])
z = np.asarray(fields['z'], float)

# z*100 not close to an integer -> third decimal digit is non-zero, i.e.
# finer than the 0.1-rounded / Gaia ~0.01 precision that z_phot flags.


spec_z_fine = np.isfinite(z) & (np.abs(z * 100 - np.round(z * 100)) > 1e-6)

sample = fields[is_qso & (z > 0.2) & (z < 1.5)
                 & (fields['total_exptime'] > 900)
                 & ~fields['z_phot'] & spec_z_fine]
print(f'{len(sample)} fields')


# 7. Save

pairs.write('muse_x_milliquas_pairs.fits', overwrite=True)
fields.write('muse_x_milliquas_fields.fits', overwrite=True)
fields.write('muse_x_milliquas_fields.csv', overwrite=True)


sample.write('muse_x_milliquas_sample.fits', overwrite=True)
sample.write('muse_x_milliquas_sample.csv', overwrite=True)