"""
stack_cubes.py
 
Stack multiple processed/masked MUSE cubes by adding them together.
 
Assumes the input cubes have already been cropped/resampled onto the
same spatial and spectral grid (as your pipeline's final output does) -
stacking only makes sense if they line up pixel-for-pixel.
 
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
from mpdaf.obj import Cube
 
 
def load_cubes(paths):
    cubes = []
    for p in paths:
        print(f"[INFO] Loading {p}")
        cubes.append(Cube(p))
    return cubes
 
 
 
def stack_cubes(cubes):
    """Sum data arrays together, respecting masks.
 
    A pixel is masked in the output only if it is masked in every
    single input cube; otherwise it's summed over whichever inputs
    are unmasked at that pixel.
    """
    stacked = cubes[0].copy()
 
    data_stack = np.ma.stack([c.data for c in cubes], axis=0)
    summed_data = data_stack.sum(axis=0)
    all_masked = np.ma.getmaskarray(data_stack).all(axis=0)
    stacked.data = np.ma.array(np.ma.getdata(summed_data), mask=all_masked)
 
    # Sum variances too, if every cube has one (errors add in quadrature -> variances add linearly)
    if all(c.var is not None for c in cubes):
        var_stack = np.ma.stack([c.var for c in cubes], axis=0)
        stacked.var = var_stack.sum(axis=0)
 
    return stacked
 
 
def main():
    parser = argparse.ArgumentParser(description="Stack (sum) multiple MUSE cubes.")
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
 
    print(f"[INFO] Stacking {len(paths)} cubes:")
    for p in paths:
        print(f"    {p}")
 
    cubes = load_cubes(paths)
 
    stacked = stack_cubes(cubes)
 
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    stacked.write(args.output)
    print(f"[INFO] Stacked cube written to {args.output}")
 
 
if __name__ == "__main__":
    main()
