"""
runs the full CASSI_26 pipeline end to end:

    1. build_sample.py     Milliquas x ESO MUSE archive crossmatch
    2. download_sample.py  download the selected cubes
    3. mainscript.py       calibrate / PSF-subtract / crop / resample each cube

Each stage writes files the next one reads (see README), so this script
just calls each stage's main() in order and passes the right filenames
between them. Re-running is cheap: build_sample's TAP queries are cached
to disk, and download_sample skips cubes already on disk.

Usage:
    python run_pipeline.py               # dry run: builds sample, shows what would be downloaded, and stops there
    python run_pipeline.py --download    # actually fetch cubes and process them
"""

import argparse

import build_sample
import download_sample
import mainscript


def run(download: bool, n_fields: int, max_gb: float, sex_config: str,
        sex_binary: str, field_indices=None):
    # 1. crossmatch -> muse_x_milliquas_{pairs,fields,sample}.fits/.csv
    print("=== 1/3: build_sample ===")
    fields, pairs, sample = build_sample.main()

    # 2. download the cubes for the chosen sample rows -> local pairs table
    print("\n=== 2/3: download_sample ===")
    dl_cfg = download_sample.Config(
        sample="muse_x_milliquas_sample.fits",
        dest="./scratch/",
        pairs_out="muse_x_milliquas_pairs_local.fits",
        field_indices=field_indices,
        n_fields=n_fields,
        max_gb=max_gb,
        download=download,
    )
    download_sample.main(dl_cfg)

    if not download:
        print("\nDry run only (pass --download to actually fetch cubes and "
              "run the processing pipeline). Stopping here.")
        return

    # 3. process every downloaded cube (wavelength conversion -> mask ->
    #    PSF subtraction -> crop -> resample, mask carried alongside)
    print("\n=== 3/3: mainscript ===")
    mainscript.main(
        sample_path="muse_x_milliquas_sample.fits",
        pairs_path="muse_x_milliquas_pairs_local.fits",
        cube_dir=dl_cfg.dest,
        sex_config=sex_config,
        sex_binary=sex_binary,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--download", action="store_true",
                    help="actually fetch cubes and run the processing pipeline "
                         "(default: dry run, stops after showing what would download)")
    p.add_argument("--n-fields", type=int, default=3,
                    help="how many top fields (by exposure) to download, if "
                         "--field-indices isn't given")
    p.add_argument("--field-indices", type=int, nargs="+", default=None,
                    help="explicit sample row indices to download instead of "
                         "picking the top --n-fields")
    p.add_argument("--max-gb", type=float, default=25.0,
                    help="size budget in GB for the download step")
    p.add_argument("--sex-config", type=str, default="./sextractor.sex",
                    help="path to your SExtractor .sex config file")
    p.add_argument("--sex-binary", type=str, default="sex",
                    help="SExtractor executable name/path (default: 'sex')")
    args = p.parse_args()

    run(download=args.download, n_fields=args.n_fields, max_gb=args.max_gb,
        sex_config=args.sex_config, sex_binary=args.sex_binary,
        field_indices=args.field_indices)
