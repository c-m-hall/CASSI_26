"""
Resume the pipeline at the resample step for EVERY cube that already has a
cropped output on disk -- regardless of run_pipeline.py's --n-fields/--max-gb
selection.

Why this exists: run_pipeline.py always re-runs download_sample.py, which
re-derives `muse_x_milliquas_pairs_local.fits` from scratch every time,
keeping only the top --n-fields fields, trimmed further to fit --max-gb.
mainscript.py only processes a cube if it has a row in that table, so any
run with modest --n-fields/--max-gb silently drops every other field, even
ones that were already cropped by an earlier, larger run.

This script sidesteps that by reading the FULL crossmatch table
(`muse_x_milliquas_pairs.fits`, written once by build_sample.py and never
budget-trimmed) to recover each cropped cube's source Name, then calls
mainscript.main(..., resume_from='resample') directly -- no build_sample,
no download_sample, no re-selection.

Requires the ORIGINAL raw cubes to still be present in --cube-dir: even in
resume_from='resample' mode, mainscript.py opens the raw cube's WCS header
(via convert_coords) before it ever touches the crop. Cubes whose raw file
is missing are skipped with a warning rather than failing the whole run.

Usage:
    python resume_from_cropped.py \
        --sex-config /carnegie/scidata/groups/muse_stacking/CASSI_26/wl_eso.sex
"""

import argparse
import glob
import os

from astropy.table import Table

import mainscript

PSFSUB_SUFFIX = "_vac_PSFSUBBED_fully_cropped.fits"
MASK_SUFFIX = "_vac_MASK3D_fully_cropped.fits"


def build_local_pairs_from_crops(crop_dir, full_pairs_path, cube_dir, out_path):
    """Scan crop_dir for cropped science+mask pairs, look each dp_id's
    source Name up in the full (unbudgeted) crossmatch table, and write a
    (Name, cube_filename, exptime) table mainscript.py can use directly."""
    full_pairs = Table.read(full_pairs_path)
    dp_id_to_name = {}
    dp_id_to_exptime = {}
    for row in full_pairs:
        dp_id = str(row["dp_id"])
        if dp_id not in dp_id_to_name:
            dp_id_to_name[dp_id] = str(row["Name"])
            dp_id_to_exptime[dp_id] = float(row["t_exptime"])

    psfsub_crops = glob.glob(os.path.join(crop_dir, f"*{PSFSUB_SUFFIX}"))

    rows = []
    skipped_no_mask, skipped_no_name, skipped_no_raw = [], [], []

    for path in sorted(psfsub_crops):
        base = os.path.basename(path)
        dp_id = base[: -len(PSFSUB_SUFFIX)]

        mask_path = os.path.join(crop_dir, dp_id + MASK_SUFFIX)
        if not os.path.exists(mask_path):
            skipped_no_mask.append(dp_id)
            continue

        if dp_id not in dp_id_to_name:
            skipped_no_name.append(dp_id)
            continue

        raw_path = os.path.join(cube_dir, dp_id + ".fits")
        if not os.path.exists(raw_path):
            skipped_no_raw.append(dp_id)
            continue

        rows.append((dp_id_to_name[dp_id], dp_id + ".fits", dp_id_to_exptime[dp_id]))

    table = Table(rows=rows, names=("Name", "cube_filename", "exptime"))
    table.write(out_path, overwrite=True)

    print(f"found {len(psfsub_crops)} cropped cube(s) in {crop_dir}")
    print(f"-> {len(rows)} ready to resume from resample, written to {out_path}")
    if skipped_no_mask:
        print(f"! {len(skipped_no_mask)} skipped (missing cropped mask): {skipped_no_mask}")
    if skipped_no_name:
        print(f"! {len(skipped_no_name)} skipped (dp_id not in {full_pairs_path}): {skipped_no_name}")
    if skipped_no_raw:
        print(f"! {len(skipped_no_raw)} skipped (raw cube missing from {cube_dir}, "
              f"needed for WCS x,y): {skipped_no_raw}")

    return out_path


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", default="muse_x_milliquas_sample.fits")
    p.add_argument("--full-pairs", default="muse_x_milliquas_pairs.fits",
                    help="the unbudgeted crossmatch table from build_sample.py "
                         "(NOT muse_x_milliquas_pairs_local.fits)")
    p.add_argument("--crop-dir", default=os.path.join(
        mainscript.BASE_DIR, "cropped") + "/")
    p.add_argument("--cube-dir", default="./scratch/")
    p.add_argument("--out-pairs-local", default="muse_x_milliquas_pairs_local_from_crops.fits")
    p.add_argument("--pixscale", type=float, default=0.2)
    p.add_argument("--sex-config", default="./sextractor.sex")
    p.add_argument("--sex-binary", default="sex")
    args = p.parse_args()

    pairs_local_path = build_local_pairs_from_crops(
        args.crop_dir, args.full_pairs, args.cube_dir, args.out_pairs_local)

    mainscript.main(
        sample_path=args.sample,
        pairs_path=pairs_local_path,
        cube_dir=args.cube_dir,
        pixscale=args.pixscale,
        sex_config=args.sex_config,
        sex_binary=args.sex_binary,
        resume_from="resample",
    )
