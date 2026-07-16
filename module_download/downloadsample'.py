"""
MUSE Phase-3 downloads for the MILLIQUAS x MUSE sample. (my version))
 
Pick fields two ways:
  1. FIELD_INDICES = [0, 5, 12]   -> those exact rows
  2. FIELD_INDICES = None         -> top N_FIELDS rows by individual-cube exposure
 
DRY RUN by default: prints how many cubes / how many GB it *would* fetch
(cheap TAP metadata query, no cube bytes) and stops. MUSE WFM cubes are
~5.5 GB each, so MAX_GB caps every real run. Downloads are resumable --
astroquery skips files already in DEST unless you force a re-download.
 
    python download_sample_example_simple.py     # dry run
    # then set DOWNLOAD = True (and MAX_GB) and re-run to actually fetch
"""
 
from pathlib import Path
import numpy as np
from astropy.table import Table
from astroquery.eso import Eso
 
SAMPLE = Path("muse_x_milliquas_sample.fits")
DEST = "./scratch/"
ESO_TAP = "https://archive.eso.org/tap_obs"
 
FIELD_INDICES = [41, 45]   # explicit row indices, or None to rank instead
N_FIELDS = 3               # used only when FIELD_INDICES is None
 
MAX_GB = 25.0               # size budget for a real run
DOWNLOAD = True             # set True to actually fetch
 
 
def select_fields():
    """Return the chosen rows, with an `idx` column giving each row's
    position in the sample file."""
    sample = Table.read(SAMPLE)
    sample["idx"] = np.arange(len(sample))
 
    if FIELD_INDICES is not None:
        n = len(sample)
        bad = [i for i in FIELD_INDICES if not -n <= i < n]
        if bad:
            raise IndexError(f"sample has {n} rows; out-of-range indices: {bad}")
        return sample[list(FIELD_INDICES)]
 
    sample.sort("exptime_indiv", reverse=True)
    return sample[:N_FIELDS]
 
 
def dp_ids_for_fields(chosen):
    """De-duplicated, order-preserving list of dp_ids from the chosen rows,
    plus a count of fields with no downloadable individual cube."""
    ids = []
    seen = set()
    for row in chosen["dp_ids"]:
        text = row if row not in (None, "--") else ""
        for dp_id in str(text).split(";"):
            if dp_id and dp_id != "--" and dp_id not in seen:
                seen.add(dp_id)
                ids.append(dp_id)
 
    n_empty = int(np.sum(chosen["n_cubes_indiv"] == 0))
    return ids, n_empty
 
 
def estimated_sizes_gb(dp_ids, chunk=100):
    """{dp_id: size_GB} from ESO ObsCore access_estsize (kB). Metadata only."""
    import pyvo
    svc = pyvo.dal.TAPService(ESO_TAP)
 
    sizes = {}
    for i in range(0, len(dp_ids), chunk):
        batch = dp_ids[i:i + chunk]
        id_list = ",".join(f"'{d}'" for d in batch)
        query = f"SELECT dp_id, access_estsize FROM ivoa.ObsCore WHERE dp_id IN ({id_list})"
        rows = svc.search(query).to_table()
        for r in rows:
            sizes[str(r["dp_id"])] = float(r["access_estsize"]) / 1e6  # kB -> GB
    return sizes
 
 
def apply_budget(dp_ids, sizes, max_gb):
    """Keep cubes in order until adding the next one would exceed max_gb."""
    kept, total = [], 0.0
    for dp_id in dp_ids:
        size = sizes.get(dp_id, 0.0)
        if kept and total + size > max_gb:
            break
        kept.append(dp_id)
        total += size
    return kept, total
 
 
def download_main():
    chosen = select_fields()
    dp_ids, n_empty = dp_ids_for_fields(chosen)
 
    how = f"sample rows {FIELD_INDICES}" if FIELD_INDICES is not None else f"top-{N_FIELDS} fields by exposure"
    hours = chosen["exptime_indiv"].sum() / 3600
    print(f"{len(chosen)} field(s) selected ({how}) -> {len(dp_ids)} individual "
          f"cube(s), {hours:.1f} hr of exposure")
    chosen["idx", "Name", "z", "n_cubes_indiv", "exptime_indiv"].pprint(max_width=-1)
 
    if n_empty:
        print(f"! {n_empty}/{len(chosen)} selected field(s) have no downloadable "
              f"per-OB cube (NFM/offset-sky/combine only) -- nothing to download.")
 
    sizes = estimated_sizes_gb(dp_ids)
    total_gb = sum(sizes.values())
    queue, queue_gb = apply_budget(dp_ids, sizes, MAX_GB)
    print(f"\nfull set ~{total_gb:.1f} GB; MAX_GB={MAX_GB} -> queue {len(queue)} of "
          f"{len(dp_ids)} cube(s), ~{queue_gb:.1f} GB, into {DEST!r}")
 
    if not DOWNLOAD:
        print("\nDRY RUN -- set DOWNLOAD=True (and raise MAX_GB if you want more) "
              "to fetch the queued cubes.")
        return
 

    Path(DEST).mkdir(parents=True, exist_ok=True)
    eso = Eso()  # public Phase-3 cubes download anonymously; no login() needed
 
    # continuation=False (default) skips cubes already in DEST -> resumable
    # re-runs. Note astroquery's flag is inverted: continuation=True FORCES a
    # re-download of files already on disk. Always pass destination=DEST, or
    # retrieve_data silently falls back to ~/.astropy/cache/astroquery/Eso.
    paths = eso.retrieve_data(queue, destination=DEST, continuation=False, unzip=True)
    if not isinstance(paths, (list, tuple)):
        paths = [paths]
 
    print(f"\nsaved {len(paths)} file(s) to {DEST!r}:")
    for p in paths:
        print("  ", p)
 
 
if __name__ == "__main__":
    download_main()