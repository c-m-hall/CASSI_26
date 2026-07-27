# stacking MUSE datacubes

This repo contains the code for preparing MUSE datacubes for stacking. 
It finds targets by cross-matching the Milliquas catalogs against the ESO MUSE archive,
downloads the cubes, then calibrates the wavelength, performs PSF-subtraction, crops the cube spatially and spectrally, and resamples
each cube onto a common grid centered on the QSO.



## Pipeline overview

```
build_sample.py  ─────►  download_sample.py  ─────►  mainscript.py
 (find targets)           (fetch the cubes)            (process each cube)
      │                          │                            │
      ▼                          ▼                            ▼
muse_x_milliquas_        muse_x_milliquas_          converted_wl/ → psf_subbed/
sample.fits               pairs_local.fits            → cropped/ → resampled/
```

`run_pipeline.py` at the repo root chains all three stages together, so you
normally don't call the individual scripts directly.

## Structure

```
build_sample.py          1. cross-match Milliquas x ESO MUSE archive,
                             apply the science cut, write the target catalogue

download_sample.py       2. download the matched individual-exposure cubes,
                             write a local Name -> cube_filename lookup table

mainscript.py             3. runs steps 3a-3d below, once per downloaded cube

module_download/          RA/Dec -> pixel conversion, cube path lookup
  extractdata.py             (get_cube_paths, get_z, convert_coords)

module_convert_wl/        3a. air-to-vacuum wavelength calibration
  convert_wavel.py

module_psf_sub/           3b. PSF / continuum subtraction, centered on the QSO
  PSFSubtraction.py

module_crop/               3c. spatial + spectral crop around the QSO,
  crop_cube.py                 centered on observed [OII]

module_resample/           3d. resample onto a common spatial/velocity grid
  resample.py

run_pipeline.py            orchestrator: runs all three stages in order
```


## How a source flows through the pipeline

1. **`build_sample.py`** queries Milliquas (VizieR TAP) and the ESO MUSE
   archive (ObsCore TAP), cross-matches them within 30″ (half the MUSE WFM
   field of view), and collapses the matches to one row per QSO. Two sets of
   per-object stats are computed:
   
   - `n_cubes` / `total_exptime` / `dp_ids` — **all** matched products,
     including combined stacks, SKY-offset exposures, and NFM cubes
   - `n_cubes_indiv` / `exptime_indiv` / `dp_ids_indiv` — **individual
     per-exposure cubes only**, i.e. the ones actually worth downloading and
     running through the rest of the pipeline

   The chosen science cut (Type-1 QSO, `0.2 < z < 1.5`, >15 min individual exposure,
   spectroscopic — not photometric — redshift) is then applied to produce
   `muse_x_milliquas_sample.fits`, the target catalogue everything downstream
   reads from. `muse_x_milliquas_fields.fits` (pre-cut) and
   `muse_x_milliquas_pairs.fits` (raw matched pairs) are also written, for
   reference/debugging.

3. **`download_sample.py`** reads the sample, picks which fields to fetch
   (top-N by individual exposure time, or explicit row indices), estimates
   total download size against a GB budget, and — only if run with
   `download=True` — pulls the queued `dp_ids_indiv` cubes from the ESO
   archive via `astroquery`. It then writes
   `muse_x_milliquas_pairs_local.fits`, a `(Name, cube_filename)` table
   mapping each QSO to its actual local file(s) on disk.

4. **`mainscript.py`** loops over every row in the sample, looks up that
   QSO's local cube path(s) via `get_cube_paths`, and for each cube:
   - converts the QSO's catalog RA/Dec to that cube's pixel `(x, y)` via
     `convert_coords` (uses that specific cube's WCS, since the same QSO
     lands at a different pixel position in every cube/pointing)
   - **3a.** air-to-vacuum wavelength correction (`convert_wl_main`)
   - **3b.** PSF/continuum subtraction, centered on `(x, y)` (`psf_sub_main`)
   - **3c.** spatial crop around `(x, y)` and spectral crop around observed
     [OII] (`crop_cube`)
   - **3d.** resample onto a common spatial/velocity grid (`resample_main`)

   Each stage's output is moved into its own subdirectory
   (`DIR_WLCONV`/`DIR_PSFSUB`/`DIR_CROP`/`DIR_FINAL`) before the next stage
   runs. A failure on any one cube is caught and logged, so one bad cube
   doesn't stop the rest of the sample from processing.

## Usage

```bash
# dry run: builds the sample, shows the download plan, stops there
python run_pipeline.py

# actually fetch cubes and run the full processing pipeline
python run_pipeline.py --download

# control what gets downloaded
python run_pipeline.py --download --n-fields 5           # top 5 fields (default 3)
python run_pipeline.py --download --field-indices 0 5 12 # specific sample rows
python run_pipeline.py --download --max-gb 50             # size budget (default 25 GB)
```

Each stage writes files the next one reads, so `build_sample.py`'s TAP
queries are cached locally (`milliquas_v8.fits`, `eso_muse_cubes.fits`) and
`download_sample.py` skips cubes already on disk — reruns are cheap. Delete
the cache files if you want a fresh archive query.

Before your first `--download` run, edit the four output paths at the top of
`mainscript.py` (`DIR_WLCONV`, `DIR_PSFSUB`, `DIR_CROP`, `DIR_FINAL`) to
somewhere that exists on your machine.

You can also run any stage standalone:

```bash
python build_sample.py       # 1. build the target catalogue
python download_sample.py    # 2. download the cubes (edit Config at the
                              #    top of the file, or import + call main()
                              #    with your own Config)
python mainscript.py         # 3. process every downloaded cube
```

## Output

Final, stack-ready cubes land in `DIR_FINAL`. Each has 3 extensions:
`DATA` (flux), `VAR` (variance), `V_GRID` (shared velocity grid).

## Requirements

```bash
pip install numpy scipy astropy mpdaf spectres pyvo astroquery matplotlib
```

`.fits`/`.csv` outputs and `__pycache__/` are gitignored — they're
regenerated by running the pipeline, not tracked as source.
