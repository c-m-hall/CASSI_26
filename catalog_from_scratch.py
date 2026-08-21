"""
Build a local pairs table (mainscript.py's `pairs_path`) from whatever cube
files are actually present in `cube_dir` right now -- decoupled from any one
download_sample.py run's --n-fields/--max-gb selection.

Why this exists: mainscript.py only processes a cube if it has a row in
`pairs_path` (Name, cube_filename, exptime). download_sample.py's own
pairs_local table only lists the cubes it queued in that ONE run, so if you
download over several separate sessions, the processing step only ever sees
the most recent batch -- everything downloaded earlier gets silently
dropped.

This script instead scans `cube_dir` for every '<dp_id>.fits' cube on disk,
looks each dp_id up in the FULL (unbudgeted) crossmatch table written once
by build_sample.py -- muse_x_milliquas_pairs.fits -- to recover its source
Name and exptime, and writes a pairs table covering everything actually
downloaded so far, regardless of which download_sample.py run fetched it or
when.

Usage (the "download now, process later" two-step workflow):

    # 1. download some cubes -- repeatable across separate sessions, each
    #    run adds to scratch/ without disturbing what's already there
    python download_sample.py

    # ... pause here. come back later, download more if you want ...

    # 2. build a catalog of EVERY cube currently in scratch/
    python catalog_from_scratch.py

    # 3. run wavelength-conversion / PSF-sub / crop / resample on all of it
    python mainscript.py --pairs-path muse_x_milliquas_pairs_local.fits

Or do steps 2+3 together with: python run_pipeline.py --catalog-from-scratch
"""

import argparse
import glob
import os

from astropy.table import Table


def build_catalog(cube_dir, full_pairs_path, out_path):
    """(Name, cube_filename, exptime) for every '<dp_id>.fits' cube found in
    `cube_dir`, looked up against the full crossmatch table. Returns
    `out_path`."""
    full_pairs = Table.read(full_pairs_path)

    dp_id_to_name, dp_id_to_exptime = {}, {}
    for row in full_pairs:
        dp_id = str(row['dp_id'])
        if dp_id not in dp_id_to_name:
            dp_id_to_name[dp_id] = str(row['Name'])
            dp_id_to_exptime[dp_id] = float(row['t_exptime'])

    cube_paths = sorted(glob.glob(os.path.join(cube_dir, '*.fits')))

    rows, skipped = [], []
    for path in cube_paths:
        cube_filename = os.path.basename(path)
        dp_id = os.path.splitext(cube_filename)[0]
        if dp_id not in dp_id_to_name:
            skipped.append(cube_filename)
            continue
        rows.append((dp_id_to_name[dp_id], cube_filename, dp_id_to_exptime[dp_id]))

    table = Table(rows=rows, names=('Name', 'cube_filename', 'exptime'))
    table.write(out_path, overwrite=True)

    print(f"found {len(cube_paths)} cube(s) in {cube_dir!r}")
    print(f"-> {len(rows)} matched to a source, written to {out_path!r}")
    if skipped:
        print(f"! {len(skipped)} cube(s) in {cube_dir!r} had no dp_id match in "
              f"{full_pairs_path!r}, skipped: {skipped}")

    return out_path


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--cube-dir', default='./scratch/',
                    help="folder to scan for downloaded cubes (default: ./scratch/)")
    p.add_argument('--full-pairs', default='muse_x_milliquas_pairs.fits',
                    help="the unbudgeted crossmatch table from build_sample.py "
                         "(NOT muse_x_milliquas_pairs_local.fits)")
    p.add_argument('--out', default='muse_x_milliquas_pairs_local.fits',
                    help="output pairs table -- pass this as --pairs-path to "
                         "mainscript.py")
    args = p.parse_args()

    build_catalog(args.cube_dir, args.full_pairs, args.out)
