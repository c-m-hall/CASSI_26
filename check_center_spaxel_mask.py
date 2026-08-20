#!/usr/bin/env python3
"""
check_center_spaxel_mask.py



Note: apply_mask.py converts 0-valued (bad) voxels to NaN in the final
*_masked.fits DATA cube rather than keeping literal 0s -- so to check for
a literal mask value of 0, run this against the MASK files themselves
(*_MASK3D.fits or *_processed_mask.fits), not the final science cube.
If you actually want to check the final masked science cube for a NaN
center spaxel instead, pass --check-nan and point --glob at 
*_masked.fits files (data extension 0 will be checked for NaN instead of
mask==0).

Center spaxel: (nx//2, ny//2), 0-indexed -- matches the CRPIX = n//2 + 1
convention used elsewhere in the pipeline to keep the source centered.

Usage:
    # Default: scan *_processed_mask.fits files in a directory
    python check_center_spaxel_mask.py /path/to/dir -o results.csv

    # Only flag if masked at ANY wavelength plane (not all)
    python check_center_spaxel_mask.py /path/to/dir --mode any

    # Only check a single reference wavelength plane (e.g. the middle one)
    python check_center_spaxel_mask.py /path/to/dir --mode plane

    # Check final *_masked.fits science cubes for NaN center spaxel instead
    python check_center_spaxel_mask.py /path/to/dir --check-nan --glob "*_masked.fits"
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits


def get_center_indices(ny, nx):
    """0-indexed center spaxel, matching CRPIX = n//2 + 1 (1-indexed) """
    return ny // 2, nx // 2


def check_cube(filepath, mode="all", check_nan=False, hdu_index=0):
    """
    Open one FITS file and evaluate its center spaxel.

    mode (mask-checking ):
        "all"   -> flag if EVERY wavelength channel at the center spaxel is 0 (bad)
        "any"   -> flag if ANY wavelength channel at the center spaxel is 0 (bad)
        "plane" -> flag if the MIDDLE wavelength channel at the center spaxel is 0 (bad)

    check_nan:
        If True, checks the data HDU for NaN at the center spaxel instead of
        checking a mask array for literal 0. Use this for *_masked.fits
        science cubes, where masking was already applied and bad voxels
        were converted to NaN rather than kept as 0.

    Returns a dict of results if the center spaxel is masked/NaN, else None.
    """
    try:
        with fits.open(filepath, memmap=False) as hdul:
            hdu = hdul[hdu_index]
            data = hdu.data
            if data is None:
                print(f"  [skip] No data in extension {hdu_index} of {filepath}")
                return None

            if data.ndim == 3:
                nwave, ny, nx = data.shape
                cy, cx = get_center_indices(ny, nx)
                spectrum = data[:, cy, cx]

                if check_nan:
                    flagged = np.all(np.isnan(spectrum)) if mode == "all" else (
                        np.any(np.isnan(spectrum)) if mode == "any"
                        else np.isnan(spectrum[nwave // 2])
                    )
                else:
                    is_bad = np.isclose(spectrum, 0)
                    if mode == "all":
                        flagged = bool(np.all(is_bad))
                    elif mode == "any":
                        flagged = bool(np.any(is_bad))
                    elif mode == "plane":
                        flagged = bool(is_bad[nwave // 2])
                    else:
                        raise ValueError(f"Unknown mode: {mode}")

            elif data.ndim == 2:
                ny, nx = data.shape
                cy, cx = get_center_indices(ny, nx)
                value = data[cy, cx]
                flagged = bool(np.isnan(value)) if check_nan else bool(np.isclose(value, 0))
          

            if not flagged:
                return None

            header = hdu.header
            redshift = header.get("Z", "NA")

            return {
                "filename": os.path.basename(filepath),
                "spaxel": f"({cx},{cy})",
                "redshift": redshift,
                "center_masked_value": "NaN" if check_nan else 0,
            }

    except Exception as exc:
        print(f"  [error] Could not process {filepath}: {exc}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Check center spaxel mask status in a directory of CASSI_26 mask/cube files."
    )
    parser.add_argument("directory", help="Directory containing mask (or masked cube) FITS files")
    parser.add_argument(
        "-o", "--output", default="center_spaxel_masked_zero.csv",
        help="Output CSV filename (default: center_spaxel_masked_zero.csv)"
    )
    parser.add_argument(
        "-g", "--glob", default="*_processed_mask.fits",
        help="Glob pattern for files to check (default: *_processed_mask.fits). "
             "Use *_MASK3D.fits for pre-resample masks, or *_masked.fits with "
             "--check-nan for final science cubes."
    )
    parser.add_argument(
        "--mode", choices=["all", "any", "plane"], default="all",
        help="For 3D (nwave, ny, nx) arrays: 'all' = every channel at the center "
             "spaxel is masked/NaN (default), 'any' = at least one channel is, "
             "'plane' = only the middle wavelength channel is checked."
    )
    parser.add_argument(
        "--check-nan", action="store_true",
        help="Check for NaN at the center spaxel (for final *_masked.fits DATA "
             "cubes) instead of checking a mask array for literal 0."
    )
    parser.add_argument(
        "--hdu", type=int, default=0,
        help="HDU index to read data/header from (default: 0, the PrimaryHDU "
             "used for both MASK and DATA in this pipeline)."
    )
    args = parser.parse_args()

    search_path = os.path.join(args.directory, args.glob)
    files = sorted(glob.glob(search_path))

    if not files:
        print(f"No files matching {search_path} found.")
        sys.exit(1)

    print(f"Found {len(files)} file(s) to check.")

    results = []
    for f in files:
        print(f"Checking {f} ...")
        result = check_cube(f, mode=args.mode, check_nan=args.check_nan, hdu_index=args.hdu)
        if result:
            results.append(result)
            print(f"  -> center spaxel flagged (z={result['redshift']})")

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile, fieldnames=["filename", "spaxel", "redshift", "center_masked_value"]
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone. {len(results)} of {len(files)} file(s) had a flagged center spaxel.")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
