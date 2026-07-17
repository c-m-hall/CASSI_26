# CASSI_26 -> stacking MUSE datacubes



This repo contains the code for preparing MUSE datacubes for stacking. 
CASSI_26 -> stacking MUSE datacubes: 

it finds targets, downloads their cubes, then calibrates, PSF-subtracts,
crops, and resamples each one onto a common grid.

# Structure

One folder per pipeline step:

  build_sample.py  find targets (cross-match Milliquas x ESO archive)

  - module_download/ download the matched cubes

  - module_convert_wl/   air-to-vacuum wavelength calibration

  - module_psf_sub/        PSF / continuum subtraction

  - module_crop/           crop around the source

  - module_resample/       resample onto a common grid

  - mainscript.py           runs all steps, per cube

# Usage

bashpython build_sample.py       # 1. build the target catalogue

python download_sample.py    # 2. download the cubes

python mainscript.py         # 3. process every cube




Each script writes files the next one reads, so you don't have to
re-query the archive or re-download cubes every time you rerun the
processing step.

Before running mainscript.py, edit the output paths at the top of
the file (DIR_WLCONV, DIR_PSFSUB, DIR_CROP, DIR_FINAL).

# Output

Final, stack-ready cubes land in DIR_FINAL. Each has 3 extensions:
DATA (flux), VAR (variance), V_GRID (shared velocity grid).

# Requirements

pip install numpy scipy astropy mpdaf spectres pyvo astroquery matplotlib
