"""
stack_cubes.py

Stack multiple processed/masked MUSE cubes using the median.

Assumes the input cubes have already been cropped/resampled onto the
same spatial and spectral grid (as your pipeline's final output does) -
stacking only makes sense if they line up pixel-for-pixel.


In addition to the stacked (median) data cube, the output FITS file now
also contains:

    - a NSTACK extension: number of cubes contributing (finite, i.e.
      not NaN) at each voxel
    - an EXPTIME extension: total exposure time contributing at each
      voxel (sum of each input cube's EXPTIME header value, over only
      the cubes that are finite/unexcluded at that voxel)
    - a VAR extension: propagated variance, computed as
      sum(var_i) / n^2, where n is the number of cubes with a finite
      variance at that voxel (this is the correct error propagation
      for an *average* of n measurements, which is what the median
      approximates for the purposes of variance bookkeeping)

Usage:
    # List cubes explicitly
    python stack_cubes.py cube1.fits cube2.fits cube3.fits -o stacked.fits

    # Or grab everything in a directory
    python stack_cubes.py --dir ..../CASSI_26/masked --pattern "*.fits" -o stacked.fits
"""

import argparse
import glob
import os
import sys

import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube


def load_cubes(paths):
    cubes = []
    for p in paths:
        print(f"[INFO] Loading {p}")
        cubes.append(Cube(p))
    return cubes


def get_exptime(cube, path):
    """Pull EXPTIME out of a cube's header, warning if it's missing."""
    header = cube.primary_header if cube.primary_header else cube.data_header
    exptime = header.get("EXPTIME", None)
    if exptime is None:
        print(f"[WARNING] No EXPTIME keyword found in {path}, treating as 0.")
        exptime = 0.0
    return float(exptime)


def stack_cubes(cubes, paths):
    """Median-combine data arrays together, treating NaN as the mask
    indicator (bad/excluded voxels are already NaN in the input cubes),
    and track per-voxel spaxel counts, exposure time, and propagated
    variance.

    A voxel ends up NaN in the output only if it's NaN in every single
    input cube; otherwise it's combined over whichever inputs are
    finite at that voxel.
    """
    stacked = cubes[0].copy()

    # stack across the cube (N) axis, NaN = excluded 
    data_stack = np.stack([np.asarray(c.data) for c in cubes], axis=0)  # (N, nz, ny, nx)
    bad_stack = np.isnan(data_stack)  # True = excluded voxel

    all_bad = bad_stack.all(axis=0)

  
    median_data = np.nanmedian(data_stack, axis=0)
    median_data[all_bad] = np.nan
    stacked.data = median_data

    #n_stack: how many cubes actually contributed (finite) at each voxel
    
    n_stack = (~bad_stack).sum(axis=0).astype(np.int32)  # (nz, ny, nx)

    # exptime: sum of EXPTIME for cubes finite at that voxel 
    exptimes = np.array([get_exptime(c, p) for c, p in zip(cubes, paths)])

    # broadcast each cube's scalar exptime, zero it out where that cube is NaN there
    exptime_per_cube = exptimes[:, None, None, None] * (~bad_stack)
    exptime_cube = exptime_per_cube.sum(axis=0)  # (nz, ny, nx)

    # variance: sum(var_i) / n^2, only where every cube has var 
    var_cube = None
    if all(c.var is not None for c in cubes):
        var_stack = np.stack([np.asarray(c.var) for c in cubes], axis=0)
        # treat var as excluded wherever the data was NaN, or the var itself is NaN
        var_bad = bad_stack | np.isnan(var_stack)
        var_stack_filled = np.where(var_bad, 0.0, var_stack)
        var_sum = var_stack_filled.sum(axis=0)  # (nz, ny, nx)

        # recompute n using var-specific validity, in case var has extra NaNs data doesn't
        n_for_var = (~var_bad).sum(axis=0)

        with np.errstate(divide="ignore", invalid="ignore"):
            var_cube = np.where(n_for_var > 0, var_sum / n_for_var**2, np.nan)
        stacked.var = var_cube
    else:
        print("[WARNING] Not every cube has a variance extension - skipping VAR propagation.")

    return stacked, n_stack, exptime_cube, var_cube


def write_output(stacked, n_stack, exptime_cube, var_cube, output_path):
    """Write the mpdaf-stacked cube, then append NSTACK/EXPTIME extensions."""
    stacked.write(output_path)

    with fits.open(output_path, mode="update") as hdul:
        nstack_hdu = fits.ImageHDU(data=n_stack.astype(np.int32), name="NSTACK")
        exptime_hdu = fits.ImageHDU(data=exptime_cube.astype(np.float32), name="EXPTIME")
        var_hdu = fits.ImageHDU(data=var_cube.astype(np.float32), name="VAR")
        hdul.append(var_hdu)
        hdul.append(nstack_hdu)
        hdul.append(exptime_hdu)
        hdul.flush()


def main():
    parser = argparse.ArgumentParser(description="Stack (median) multiple MUSE cubes.")
    parser.add_argument("cubes", nargs="*", help="Paths to cube FITS files to stack")
    parser.add_argument("--dir", help="Directory to glob cubes from, instead of listing paths")
    parser.add_argument("--pattern", default="*.fits", help="Glob pattern used with --dir (default: *.fits)")
    parser.add_argument("-o", "--output", required=True, help="Output FITS path for the stacked cube")
    args = parser.parse_args()

    if args.dir:
        paths = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    else:
        paths = args.cubes

    if len(paths) < 2:
        sys.exit("[ERROR] Need at least 2 cubes to stack.")

    print(f"[INFO] Stacking {len(paths)} cubes (median):")
    for p in paths:
        print(f"    {p}")

    cubes = load_cubes(paths)

    stacked, n_stack, exptime_cube, var_cube = stack_cubes(cubes, paths)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    write_output(stacked, n_stack, exptime_cube, var_cube, args.output)
    print(f"[INFO] Stacked cube written to {args.output}")
    print("[INFO]   DATA extension: median-combined flux")
    print("[INFO]   VAR extension:  sum(var_i) / n^2")
    print("[INFO]   NSTACK extension: number of cubes contributing per voxel")
    print("[INFO]   EXPTIME extension: summed exposure time per voxel")


if __name__ == "__main__":
    main()
