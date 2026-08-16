"""
build_pairs_from_existing.py
 
Build muse_x_milliquas_pairs_local.fits from cubes that are ALREADY on
disk in `dest` -- for when a download was stopped early (Ctrl+C, Slurm
time limit, etc.) before download_sample.py reached its final
pairs_table.write() call, so no pairs table exists yet even though the
cube files themselves are safely downloaded.
 
Re-uses download_sample.py's own field-selection and matching logic, so
this produces the exact same pairs table download_sample.py would have
written -- just built from what's on disk right now, rather than
waiting for every queued cube to finish downloading.
 
Usage:
    python build_pairs_from_existing.py --field-indices 97 98 99 127 137 ... \
        --sample muse_x_milliquas_sample.fits --dest ./scratch/ \
        --pairs-out muse_x_milliquas_pairs_local.fits
 
    # or, to match a run that used --n-fields instead of explicit indices:
    python build_pairs_from_existing.py --n-fields 150 \
        --sample muse_x_milliquas_sample.fits --dest ./scratch/
"""
 
import argparse
 
from download_sample import (
    Config,
    select_fields,
    dp_ids_for_fields,
    build_local_pairs_table,
)
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Build the local pairs table from cubes already on disk."
    )
    parser.add_argument("--sample", default="muse_x_milliquas_sample.fits")
    parser.add_argument("--dest", default="./scratch/",
                         help="Directory the cubes were downloaded into")
    parser.add_argument("--pairs-out", default="muse_x_milliquas_pairs_local.fits")
    parser.add_argument("--field-indices", type=int, nargs="+", default=None,
                         help="Same --field-indices used for the original download")
    parser.add_argument("--n-fields", type=int, default=3,
                         help="Same --n-fields used for the original download, "
                              "if --field-indices wasn't used")
    args = parser.parse_args()
 
    cfg = Config(
        sample=args.sample,
        dest=args.dest,
        pairs_out=args.pairs_out,
        field_indices=args.field_indices,
        n_fields=args.n_fields,
    )
 
    chosen = select_fields(cfg)
    dp_ids, name_for_dpid, n_empty = dp_ids_for_fields(chosen)
 
    print(f"{len(chosen)} field(s) selected -> {len(dp_ids)} individual cube(s) "
          f"expected; checking {cfg.dest!r} for what's actually on disk...")
 
    pairs_table, missing = build_local_pairs_table(dp_ids, name_for_dpid, cfg.dest)
    pairs_table.write(cfg.pairs_out, overwrite=True)
 
    print(f"\nfound {len(pairs_table)}/{len(dp_ids)} cube(s) on disk")
    if missing:
        print(f"  {len(missing)} not yet downloaded (fine -- these just aren't "
              f"in the pairs table yet; rerun this script again later once "
              f"more finish downloading)")
    print(f"wrote {cfg.pairs_out!r} -- ready to pass as pairs_path to mainscript.main()")
 
 
if __name__ == "__main__":
    main()
