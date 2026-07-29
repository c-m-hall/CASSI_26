"""
Build a 3D (wave, y, x) bad-voxel mask for a MUSE cube, run after
air-to-vacuum wavelength conversion and before PSF subtraction.

Two masks are built and OR-combined:

  spatial  (y, x)      -- SExtractor segmentation on the cube's white-light
                           image; True where a detected source (anything
                           that isn't blank sky) sits
  spectral (wave,)      -- sky-subtraction residual masking, per-field: robust-sigma
                           thresholding of the (baseline-subtracted) median
                           background-spaxel spectrum, then morphological
                           closing/dilation to patch gaps and cover line wings

The combined 3D mask is True at every voxel that is EITHER object-contaminated
(that spatial pixel, every wavelength) OR sky-contaminated (that wavelength,
every spatial pixel), and is saved as its own FITS file alongside the input
cube.

Requires a command-line SExtractor binary (`sex` by default) and an existing
.sex config file -- for SExtractor 
"""

import os
import subprocess
import tempfile

import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube
from scipy.ndimage import binary_closing, binary_dilation, median_filter


# 
# spatial mask: SExtractor segmentation on the white-light image


def make_whitelight_image(cube):
    """Collapse an mpdaf Cube along wavelength -> mpdaf Image (sum), which
    carries the cube's spatial WCS through automatically."""
    return cube.sum(axis=0)


def run_sextractor(image_path, sex_config, sex_binary='sex', workdir=None):
    """Run SExtractor on `image_path`

    Returns the path to the segmentation map FITS file.
    """
    workdir = workdir or os.path.dirname(image_path) or '.'
    catalog_path = os.path.join(workdir, '_sextractor_tmp.cat')
    segmap_path = os.path.join(workdir, '_sextractor_tmp_segmap.fits')

    cmd = [
        sex_binary, image_path,
        '-c', sex_config,
        '-CATALOG_NAME', catalog_path,
        '-CHECKIMAGE_TYPE', 'segmentation',
        '-CHECKIMAGE_NAME', segmap_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"SExtractor failed on {image_path} (cmd: {' '.join(cmd)}):\n"
            f"{result.stderr}"
        )
    if not os.path.exists(segmap_path):
        raise RuntimeError(
            f"SExtractor exited cleanly but no segmentation map was written "
            f"to {segmap_path}. Check that your config's DETECT settings "
            f"actually find sources in {image_path}."
        )
    return segmap_path


def spatial_mask_from_segmap(segmap_path):
    """SExtractor segmentation map -> boolean 2D mask, True on any detected
    source (segmap value > 0, i.e. not background)."""
    segmap = fits.getdata(segmap_path)
    return segmap > 0


def build_spatial_mask(cube, sex_config, sex_binary='sex', workdir=None):
    """White-light image -> SExtractor -> boolean object mask (ny, nx)."""
    workdir = workdir or tempfile.mkdtemp(prefix='sex_')
    os.makedirs(workdir, exist_ok=True)

    whitelight = make_whitelight_image(cube)
    wl_path = os.path.join(workdir, '_whitelight_tmp.fits')
    whitelight.write(wl_path)

    segmap_path = run_sextractor(wl_path, sex_config, sex_binary=sex_binary,
                                  workdir=workdir)
    return spatial_mask_from_segmap(segmap_path)



# spectral mask: sky-residual masking (same algorithm as build_sky_mask.ipynb)


def residual_sky_spectrum(cube, spatial_mask):
    """Median spectrum across background (non-object) spaxels only -- 
    per-field residual sky spectrum for input
    """
    wave = cube.wave.coord()
    flux = cube.data.data  # (nwave, ny, nx)
    bg = ~spatial_mask
    if not bg.any():
        raise ValueError("Spatial mask flags every spaxel as a source -- "
                          "no background spaxels left to build a sky "
                          "spectrum from. Check the SExtractor detection "
                          "threshold.")
    sky_spec = np.nanmedian(flux[:, bg], axis=1)
    return wave, sky_spec


def build_spectral_mask(wave, flux, window=151, k=3.0, gap=9, grow=2):
    """Baseline-subtract, robust-sigma threshold (two-sided), patch gaps,
    and grow line wings - Uses build_sky_mask.ipynb as a basis 
    """
    good = np.isfinite(flux)
    ffill = np.interp(wave, wave[good], flux[good])

    # Step 1: remove the slow baseline
    base = median_filter(ffill, size=window)
    r = ffill - base

    # Step 2: robust noise estimate (MAD)
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))

    # Step 3: two-sided threshold
    mask = np.abs(r) > k * sig

    # Step 4: patch gaps, grow edges, flag original NaNs
    mask = binary_closing(mask, structure=np.ones(gap))
    mask = binary_dilation(mask, iterations=grow)
    mask |= ~good

    return mask



# combine + save


def combine_masks(spatial_mask, spectral_mask):
    """(ny, nx) spatial + (nwave,) spectral -> (nwave, ny, nx) combined mask,
    True where either mask flags that voxel  """
    return spectral_mask[:, None, None] | spatial_mask[None, :, :]


def build_mask_main(cubepath, sex_config, sex_binary='sex',
                     window=151, k=3.0, gap=9, grow=2, workdir=None):
    cube = Cube(cubepath)

    spatial_mask = build_spatial_mask(cube, sex_config, sex_binary=sex_binary,
                                       workdir=workdir)
    wave, sky_spec = residual_sky_spectrum(cube, spatial_mask)
    spectral_mask = build_spectral_mask(wave, sky_spec, window=window, k=k,
                                         gap=gap, grow=grow)

    bad3d = combine_masks(spatial_mask, spectral_mask)   # True = bad (internal convention)

    # flip to on-disk convention: 1 = good (keep), 0 = bad (masked)
    mask3d = ~bad3d

    outpath = cubepath.replace('.fits', '_MASK3D.fits')
    hdu = fits.PrimaryHDU(mask3d.astype(np.uint8), header=cube.data_header)
    hdu.header['COMMENT'] = 'Combined mask: 1 = good voxel (keep), 0 = bad (object OR sky)'
    hdu.writeto(outpath, overwrite=True)

    return outpath
# example usage
'''
build_mask_main('cubepath_VACUUM.fits', sex_config='./my_config.sex')
'''
