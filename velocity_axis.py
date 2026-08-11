
"""
velocity_axis.py

Rewrites a cube's spectral WCS (extension [1] header, where MUSE cube data
lives) from rest-frame wavelength into velocity [km/s], centered at 0 at
the rest wavelength of the reference line (OII_REST = 3728.48 A by
default). This is a linear rescale, since the wavelength axis is already
linear in pixel index (CRVAL3/CDELT3/CRPIX3) -- no resampling of the data
itself is needed, only the header keywords change.

    v = c * (lambda - lambda_rest) / lambda_rest      [non-relativistic,
                                                         fine for the narrow
                                                         windows used here]

Usage as a library:
    from velocity_axis import set_velocity_axis
    set_velocity_axis("some_cube.fits", rest_wave=3728.48)

Usage as a script (rewrites one or more files in place):
    python velocity_axis.py cube1.fits cube2.fits --rest-wave 3728.48
"""

import argparse

from astropy.io import fits

C_KMS = 299792.458


def set_velocity_axis(path, rest_wave, ext=1):
    """Rewrite extension `ext`'s spectral WCS keywords (CTYPE3/CRVAL3/
    CDELT3/CUNIT3) so the spectral axis reads out in km/s, centered at 0
    at `rest_wave` [Angstrom]. CRPIX3 is left untouched -- only the value
    and step size the axis represents change, not which pixel is the
    reference pixel.

    Handles both CDELT3-style and CD3_3-style linear WCS.
    """
    with fits.open(path, mode="update") as hdul:
        header = hdul[ext].header

        crval3 = header["CRVAL3"]
        crpix3 = header["CRPIX3"]

        if "CDELT3" in header:
            delta_key = "CDELT3"
        elif "CD3_3" in header:
            delta_key = "CD3_3"
        else:
            raise KeyError(
                f"{path}: no CDELT3 or CD3_3 keyword found on extension {ext}; "
                "can't determine spectral pixel scale."
            )
        cdelt3 = header[delta_key]

        # linear rescale: velocity at the reference pixel, and velocity
        # step per pixel, both derived from the current wavelength values
        v_ref = C_KMS * (crval3 - rest_wave) / rest_wave
        v_delta = C_KMS * cdelt3 / rest_wave

        header["CRVAL3"] = v_ref
        header[delta_key] = v_delta
        header["CTYPE3"] = "VELO"
        header["CUNIT3"] = "km/s"
        # CRPIX3 unchanged -- same reference pixel, new value/step at it

        hdul.flush()

    print(f"[INFO] {path}: spectral axis (ext {ext}) now velocity, "
          f"0 km/s at {rest_wave} A (v_ref={v_ref:.2f}, v_delta={v_delta:.4f} km/s/pix)")


def main():
    parser = argparse.ArgumentParser(
        description="Rewrite a cube's spectral WCS from rest wavelength to velocity (in place)."
    )
    parser.add_argument("cubes", nargs="+", help="FITS cube path(s) to rewrite in place")
    parser.add_argument("--rest-wave", type=float, default=3728.48,
                         help="Rest wavelength [Angstrom] to center velocity=0 on (default: [OII] 3728.48)")
    parser.add_argument("--ext", type=int, default=1,
                         help="FITS extension holding the spectral WCS (default: 1, the DATA extension)")
    args = parser.parse_args()

    for path in args.cubes:
        set_velocity_axis(path, rest_wave=args.rest_wave, ext=args.ext)


if __name__ == "__main__":
    main()