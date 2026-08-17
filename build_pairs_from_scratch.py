"""
build_pairs_from_scratch.py
 
Build muse_x_milliquas_pairs_local.fits from whatever cube files are
ACTUALLY sitting in the scratch/download directory right now -- no need
to remember/replay the exact --field-indices or --n-fields used for the
original download. Matches purely off the full dp_id -> Name mapping in
muse_x_milliquas_pairs.fits (written by build_sample.py), so it's robust
to a download that was stopped early, partially failed, or whose exact
field selection you don't have handy anymore.
 
Then immediately runs the rest of the pipeline (mainscript.main(), i.e.
convert_wl onward) on exactly the cubes found.
 
Usage:
    python build_pairs_from_scratch.py --scratch-dir ./scratch/
"""
 
import argparse
import os
 
import numpy as np
from astropy.table import Table
 
import mainscript
 
 
def build_pairs_from_scratch(pairs_path, scratch_dir, pairs_out):
    """(Name, cube_filename, exptime) for every cube actually present in
    scratch_dir, matched via the full dp_id -> Name/t_exptime table.
 
    Assumes ESO Phase-3 downloads are saved as <dp_id>.fits (the
    astroquery default), so a file's dp_id is just its filename stem.
    """
    pairs = Table.read(pairs_path)  # full table: has 'dp_id', 'Name', 't_exptime' columns
    dpid_to_name = dict(zip(pairs['dp_id'].astype(str), pairs['Name'].astype(str)))
    dpid_to_exptime = dict(zip(pairs['dp_id'].astype(str), pairs['t_exptime']))
 
    rows = []
    unmatched = []
    for fname in sorted(os.listdir(scratch_dir)):
        if not fname.endswith('.fits'):
            continue
        dp_id = fname[:-len('.fits')]
        name = dpid_to_name.get(dp_id)
        if name is None:
            unmatched.append(fname)
            continue
        exptime = float(dpid_to_exptime.get(dp_id, np.nan))
        rows.append((name, fname, exptime))
 
    table = Table(rows=rows, names=('Name', 'cube_filename', 'exptime'))
    table.write(pairs_out, overwrite=True)
 
    print(f"{len(table)} cube(s) in {scratch_dir!r} matched to a QSO name "
          f"-> wrote {pairs_out!r}")
    if unmatched:
        print(f"! {len(unmatched)} .fits file(s) in {scratch_dir!r} had no "
              f"matching dp_id in {pairs_path!r} (not a MUSE cube from this "
              f"crossmatch? typo in filename?): {unmatched}")
 
    return table
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Build the local pairs table from scratch/ contents, then run the pipeline on it."
    )
    parser.add_argument("--pairs-path", default="muse_x_milliquas_pairs.fits",
                         help="Full dp_id/Name pairs table from build_sample.py")
    parser.add_argument("--scratch-dir", default="./scratch/",
                         help="Directory containing downloaded cube .fits files")
    parser.add_argument("--pairs-out", default="muse_x_milliquas_pairs_local.fits")
    parser.add_argument("--sample-path", default="muse_x_milliquas_sample.fits")
    parser.add_argument("--sex-config", required=True)
    parser.add_argument("--sex-binary", default="sex")
    args = parser.parse_args()
 
    build_pairs_from_scratch(args.pairs_path, args.scratch_dir, args.pairs_out)
 
    print("\nRunning mainscript.main() on the cubes found...")
    mainscript.main(
        sample_path=args.sample_path,
        pairs_path=args.pairs_out,
        cube_dir=args.scratch_dir,
        sex_config=args.sex_config,
        sex_binary=args.sex_binary,
        resume_from='download',
    )
 
 
if __name__ == "__main__":
    main()
 
