#!/usr/bin/env python3
"""
check_center_spaxel_mask.py

For each cube in a directory, check whether the spatial CENTER spaxel was
masked out, and if so, log it to a CSV with the spaxel coords, filename,
and redshift ('Z' header keyword).

Written to match the CASSI_26 pipeline convention (github.com/c-m-hall/CASSI_26):
  - Mask files (e.g. *_MASK3D.fits, *_processed_mask.fits): single-extension,
    values are 1 = good, 0 = bad (object or sky contaminated).
  - "Masked cubes" (produced by module_mask/apply_mask.py -> *_masked.fits):
    DATA is in ext 0, VAR in ext 1. Voxels where the mask was 0 are set to
    NaN (masked_data[mask_data == 0] = np.nan) -- the literal value 0 is
    NOT preserved in the final masked cube, NaN is the on-disk signature of
    a masked voxel there.
  - Redshift lives in header['Z'] of the primary HDU (set in resample.py for
    both cubes and masks).
  - Spatial center = (nx//2, ny//2) zero-indexed, matching CRPIX1/CRPIX2 =
    nx//2+1, ny//2+1 set in resample.py (this is the QSO position).

Two modes:
  --mode cube  (default): point this at final masked science cubes
                (e.g. *_masked.fits). Center spaxel is flagged if its DATA
                value is NaN (which is what apply_mask.py produces wherever
                the mask was 0).
  --mode mask  : point this at raw 0/1 mask files (e.g. *_MASK3D.fits,
                *_processed_mask.fits). Center spaxel is flagged if its
                literal value is 0.

For 3D cubes, in addition to flagging the center spaxel, the CSV now also
reports which wavelength channels at that spaxel are NOT masked (i.e. still
non-zero / good), so you can see the extent of partial masking rather than
just a pass/fail. New columns:
    n_channels_total     - total number of wavelength channels checked
    n_channels_masked    - how many of those are masked (0, or NaN in cube mode)
    n_channels_nonzero   - how many are NOT masked (good / non-zero)
    nonzero_channel_idx  - ';'-separated list of the good/non-zero channel
                            indices (capped at MAX_LISTED_CHANNELS to keep
                            the CSV readable; a summary count is always given
                            in n_channels_nonzero regardless of the cap)

Usage:
    python check_center_spaxel_mask.py /path/to/masked_cubes -o results.csv
    python check_center_spaxel_mask.py /path/to/masks --mode mask -o results.csv
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits


REDSHIFT_KEY = "Z"
MAX_LISTED_CHANNELS = 25  # cap how many channel indices we print into the CSV cell


def summarize_channels(is_masked_per_channel):
    """
    Given a 1D boolean array (True = masked/bad at that channel), return a
    dict with counts and a capped, semicolon-joined string of the NON-zero
    (i.e. good / unmasked) channel indices.
    """
    is_masked_per_channel = np.asarray(is_masked_per_channel)
    n_total = is_masked_per_channel.size
    nonzero_idx = np.flatnonzero(~is_masked_per_channel)  # indices that are good/non-zero
    n_nonzero = nonzero_idx.size
    n_masked = n_total - n_nonzero

    if n_nonzero > MAX_LISTED_CHANNELS:
        shown = nonzero_idx[:MAX_LISTED_CHANNELS]
        idx_str = ";".join(str(i) for i in shown) + f";...(+{n_nonzero - MAX_LISTED_CHANNELS} more)"
    else:
        idx_str = ";".join(str(i) for i in nonzero_idx)

    return {
        "n_channels_total": n_total,
        "n_channels_masked": n_masked,
        "n_channels_nonzero": n_nonzero,
        "nonzero_channel_idx": idx_str,
    }


def get_center_indices(shape):
    """(row, col) i.e. (y, x) index of the spatial center for a 2D/3D array."""
    if len(shape) == 3:
        _, ny, nx = shape
    elif len(shape) == 2:
        ny, nx = shape
    else:
        raise ValueError(f"Unsupported array shape: {shape}")
    return ny // 2, nx // 2


def check_cube_mode(filepath, wave_check="mid"):
    """
    Final masked science cube (DATA ext0, VAR ext1). A masked voxel is NaN
    on disk (apply_mask.py sets mask==0 voxels to NaN). Flag the center
    spaxel if it's NaN.

    wave_check: "mid" checks a single reference wavelength plane (fast, and
    sufficient here since the spatial/object mask -- which is what flags the
    QSO center spaxel -- is applied identically across all wavelengths).
    Use "all" to require every wavelength plane at the center spaxel be NaN.
    """
    with fits.open(filepath, memmap=False) as hdul:
        data = hdul[0].data
        if data is None:
            print(f"  [skip] No DATA in ext 0 of {filepath}")
            return None

        cy, cx = get_center_indices(data.shape)

        channel_summary = {}
        if data.ndim == 3:
            spectrum_is_masked = np.isnan(data[:, cy, cx])  # per-channel bool, True=masked
            channel_summary = summarize_channels(spectrum_is_masked)
            if wave_check == "all":
                is_masked = spectrum_is_masked.all()
            else:
                mid = data.shape[0] // 2
                is_masked = spectrum_is_masked[mid]
        else:
            is_masked = np.isnan(data[cy, cx])

        if not is_masked:
            return None

        redshift = hdul[0].header.get(REDSHIFT_KEY, "NA")
        result = {
            "filename": os.path.basename(filepath),
            "spaxel": f"({cx},{cy})",
            "redshift": redshift,
            "center_masked_value": 0,
        }
        result.update(channel_summary)
        return result


def check_mask_mode(filepath):
    """
    Raw 0/1 mask file. Flag the center spaxel if its literal value is 0.
    Looks for a 'MASK' extension by name first, falls back to ext 0.
    """
    with fits.open(filepath, memmap=False) as hdul:
        mask_hdu = None
        for hdu in hdul:
            if hdu.name == "MASK" and hdu.data is not None:
                mask_hdu = hdu
                break
        if mask_hdu is None:
            mask_hdu = hdul[0]

        mask_data = mask_hdu.data
        if mask_data is None:
            print(f"  [skip] No mask data found in {filepath}")
            return None

        cy, cx = get_center_indices(mask_data.shape)

        channel_summary = {}
        if mask_data.ndim == 3:
            spectrum = mask_data[:, cy, cx]
            spectrum_is_masked = np.isclose(spectrum, 0)  # per-channel bool, True=masked(bad)
            channel_summary = summarize_channels(spectrum_is_masked)
            center_value = spectrum[spectrum.shape[0] // 2]
        else:
            center_value = mask_data[cy, cx]

        if not np.isclose(center_value, 0):
            return None

        redshift = hdul[0].header.get(REDSHIFT_KEY, "NA")
        result = {
            "filename": os.path.basename(filepath),
            "spaxel": f"({cx},{cy})",
            "redshift": redshift,
            "center_masked_value": 0,
        }
        result.update(channel_summary)
        return result


def main():
    parser = argparse.ArgumentParser(
        description="Check the center spaxel mask status of cubes in a directory."
    )
    parser.add_argument("directory", help="Directory containing the FITS files")
    parser.add_argument(
        "-o", "--output", default="center_spaxel_masked_zero.csv",
        help="Output CSV filename (default: center_spaxel_masked_zero.csv)"
    )
    parser.add_argument(
        "-g", "--glob", default=None,
        help="Glob pattern for files (default: '*_masked.fits' for cube mode, "
             "'*.fits' for mask mode)"
    )
    parser.add_argument(
        "--mode", choices=["cube", "mask"], default="cube",
        help="'cube' (default): final masked science cubes, flags NaN center "
             "spaxels. 'mask' : raw 0/1 mask files, flags literal-0 center "
             "spaxels."
    )
    parser.add_argument(
        "--wave-check", choices=["mid", "all"], default="mid",
        help="(cube mode only) 'mid': check one reference wavelength plane "
             "(default, faster). 'all': require every wavelength plane at "
             "the center spaxel be masked."
    )
    args = parser.parse_args()

    file_glob = args.glob or ("*_masked.fits" if args.mode == "cube" else "*.fits")
    search_path = os.path.join(args.directory, file_glob)
    files = sorted(glob.glob(search_path))

    if not files:
        print(f"No files matching {search_path} found.")
        sys.exit(1)

    print(f"Found {len(files)} file(s) to check (mode={args.mode}).")

    results = []
    for f in files:
        print(f"Checking {f} ...")
        try:
            if args.mode == "cube":
                result = check_cube_mode(f, wave_check=args.wave_check)
            else:
                result = check_mask_mode(f)
        except Exception as exc:
            print(f"  [error] Could not process {f}: {exc}")
            continue

        if result:
            results.append(result)
            print(f"  -> center spaxel masked (0) (z={result['redshift']})")

    fieldnames = [
        "filename", "spaxel", "redshift", "center_masked_value",
        "n_channels_total", "n_channels_masked", "n_channels_nonzero",
        "nonzero_channel_idx",
    ]
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone. {len(results)} of {len(files)} file(s) had a masked center spaxel.")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
