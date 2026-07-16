import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits


def air2vac(wvl, precision=1.0e-12, maxiter=100):
    """
    Convert wavelength in air to wavelength in vacuum.
    Using Eq. 8 from Morton 2000. Solved iteratively.

    Parameters
    ----------
    wvl : float or array
        Wavelength in air
    precision : float, optional
        Precision beyond which iteration stops; default 1.0e-12
    maxiter : int, optional
        Maximum number of iterations; default 100

    Returns
    -------
    wvl : float or array
        Wavelength in vacuum
    """
    wvl = np.asarray(wvl, dtype="float")
    lvac = wvl.copy()
    count = 0
    while True:
        count += 1
        s = 1.0e4 / lvac
        n = 1.0 + 8.34254e-5 + 2.406147e-2 / (130 - s**2) + 1.5998e-4 / (38.9 - s**2)
        lair = lvac / n
        dl = lair - wvl
        lvac -= dl
        s = 1.0e4 / lvac
        n = 1.0 + 8.34254e-5 + 2.406147e-2 / (130 - s**2) + 1.5998e-4 / (38.9 - s**2)
        dlair = lvac / n - wvl
        if np.abs(np.max(dlair)) < precision:
            return lvac
        if count > maxiter:
            p = np.max(dlair)
            raise ValueError(f"air2vac: max iteration reached, current precision is {p}")


def convert_wl_main(cubename, dir):
    datapath = dir + cubename
    outpath = dir + cubename.replace(".fits", "_vac.fits")

    with fits.open(datapath) as hdul:
        d, v = hdul[1], hdul[2]
        data, dhdr = d.data, d.header
        var, vhdr = v.data, v.header

        npix = dhdr["NAXIS3"]
        lstep = dhdr["CD3_3"]
        wave0 = dhdr["CRVAL3"]
        wave = wave0 + lstep * np.arange(npix)
        wavev = air2vac(wave)

        # common output grid, resampled to a fixed 1.25 A step
        newwave = np.arange(np.min(wavev), np.max(wavev), 1.25)

        # --- vectorized resample: build ONE interpolator over the whole
        # cube (axis=0 broadcasts over every spaxel simultaneously) instead
        # of constructing ny*nx separate interp1d objects in a Python loop.
        # This also fixes the original shape-mismatch bug, since we now
        # allocate output arrays sized for `newwave` up front rather than
        # assigning into the old, differently-sized `data` array in place.
        interp_flux = interp1d(
            wavev, data, axis=0, kind="linear", bounds_error=False, fill_value=np.nan
        )
        interp_var = interp1d(
            wavev, var, axis=0, kind="linear", bounds_error=False, fill_value=np.nan
        )
        new_data = interp_flux(newwave)
        new_var = interp_var(newwave)

        hdul[1].data = new_data
        hdul[2].data = new_var

        dhdr["NAXIS3"] = len(newwave)
        vhdr["NAXIS3"] = len(newwave)
        dhdr["CRVAL3"] = np.min(newwave)
        vhdr["CRVAL3"] = np.min(newwave)
        dhdr["CD3_3"] = 1.25
        vhdr["CD3_3"] = 1.25
        dhdr["CRPIX3"] = 1
        vhdr["CRPIX3"] = 1

        hdul.writeto(outpath, overwrite=True)

    return outpath


# example usage
"""
dir = '/Users/charishall/CASSI_26/'
convert_wl_main('J0015m0751_dr2_zap_wpsf1.corrEBV.fits', dir)
"""
