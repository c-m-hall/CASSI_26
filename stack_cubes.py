"""
stack_cubes.py

Stack multiple processed/masked MUSE cubes together.

Masked spaxels/voxels in the input cubes are stored as NaN, so they're
simply ignored (via NaN-aware numpy functions) rather than tracked with
a separate mask array.

Assumes the input cubes have already been cropped/resampled onto the
same spatial and spectral grid (as your pipeline's final output does) -
stacking only makes sense if they line up voxel-for-voxel.

Both a mean stack and a median stack are always computed; --method picks
which one becomes the primary DATA extension in the output file, and the
mean stack is always additionally saved as a MEAN extension (with its
variance as MEAN_VAR) regardless of --method.

Usage:
    # List cubes explicitly
    python stack_cubes.py cube1.fits cube2.fits cube3.fits -o stacked.fits

    # Or grab everything in a directory
    python stack_cubes.py --dir /scratch/$USER/CASSI_26/masked --pattern "*.fits" -o stacked.fits --method median
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


def _combine(data_stack, n_map, method):
    """Combine an (N, ...) stack of arrays along axis 0, NaN-aware."""
    all_nan = n_map == 0
    with np.errstate(invalid='ignore'):
        if method == "sum":
            out = np.nansum(data_stack, axis=0)
        elif method == "mean":
            out = np.nanmean(data_stack, axis=0)
        elif method == "median":
            out = np.nanmedian(data_stack, axis=0)
        else:
            raise ValueError(f"Unknown method '{method}'; expected 'sum', 'mean', or 'median'")
    out[all_nan] = np.nan
    return out


def stack_cubes(cubes, method="median"):
    """Combine data arrays together, ignoring NaNs (masked spaxels are NaN).

    Always computes both a mean stack and a median stack; `method` picks
    which one (or 'sum') is returned as the primary result.

    method : {'sum', 'mean', 'median'}
        Which combination becomes the primary output.

    A voxel is NaN in the output only if it's NaN in every single input
    cube; otherwise it's combined over whichever inputs are finite there.

    Returns
    -------
    primary : mpdaf.obj.Cube
        The `method`-combined cube, with .var set to sum(var) / n**2.
    mean_cube : mpdaf.obj.Cube
        The mean-combined cube (same var formula), for saving as an
        extra extension regardless of what `method` was.
    n_map : numpy.ndarray
        Per-voxel count of how many input cubes contributed (same shape
        as the data).
    """
    data_stack = np.stack(
        [np.ma.filled(c.data, np.nan).astype(float) for c in cubes], axis=0
    )
    finite_stack = ~np.isnan(data_stack)
    n_map = finite_stack.sum(axis=0)
    n_safe = np.maximum(n_map, 1)
    all_nan = n_map == 0

    has_var = all(c.var is not None for c in cubes)
    var_out = None
    if has_var:
        var_stack = np.stack(
            [np.ma.filled(c.var, np.nan).astype(float) for c in cubes], axis=0
        )
        with np.errstate(invalid='ignore'):
            var_out = np.nansum(var_stack, axis=0) / n_safe**2
        var_out[all_nan] = np.nan

    def make_cube(method_):
        c = cubes[0].copy()
        c.data = _combine(data_stack, n_map, method_)
        if has_var:
            c.var = var_out.copy()
        return c

    primary = make_cube(method)
    mean_cube = primary if method == "mean" else make_cube("mean")

    return primary, mean_cube, n_map


def main():
    parser = argparse.ArgumentParser(description="Stack multiple MUSE cubes (sum, mean, or median).")
    parser.add_argument("cubes", nargs="*", help="Paths to cube FITS files to stack")
    parser.add_argument("--dir", help="Directory to glob cubes from, instead of listing paths")
    parser.add_argument("--pattern", default="*.fits", help="Glob pattern used with --dir (default: *.fits)")
    parser.add_argument("--method", default="median", choices=["sum", "mean", "median"],
                         help="Which combination becomes the primary DATA extension "
                              "(default: median). The mean stack is always saved too, "
                              "as an extra MEAN extension.")
    parser.add_argument("-o", "--output", required=True, help="Output FITS path for the stacked cube")
    args = parser.parse_args()

    if args.dir:
        paths = sorted(glob.glob(os.path.join(args.dir, args.pattern)))
    else:
        paths = args.cubes

    if len(paths) < 2:
        sys.exit("[ERROR] Need at least 2 cubes to stack.")

    print(f"[INFO] Stacking {len(paths)} cubes using method='{args.method}':")
    for p in paths:
        print(f"    {p}")

    cubes = load_cubes(paths)

    primary, mean_cube, n_map = stack_cubes(cubes, method=args.method)

    primary.primary_header['NCUBES'] = (len(cubes), 'number of input cubes stacked')
    primary.primary_header['STACKMTH'] = (args.method, 'method used for primary DATA extension')

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    primary.write(args.output)

    with fits.open(args.output, mode='update') as hdul:
        n_hdu = fits.ImageHDU(data=n_map.astype(np.int16), name='NCUBES_MAP')
        n_hdu.header['COMMENT'] = 'number of input cubes contributing at each voxel'
        hdul.append(n_hdu)

        mean_hdu = fits.ImageHDU(data=np.asarray(mean_cube.data), name='MEAN')
        mean_hdu.header['COMMENT'] = 'mean-combined stack (saved regardless of --method)'
        hdul.append(mean_hdu)
        if mean_cube.var is not None:
            meanvar_hdu = fits.ImageHDU(data=np.asarray(mean_cube.var), name='MEAN_VAR')
            hdul.append(meanvar_hdu)

        hdul.flush()

    print(f"[INFO] Stacked cube written to {args.output} "
          f"(primary='{args.method}', plus NCUBES_MAP and MEAN extensions)")


if __name__ == "__main__":
    main()
