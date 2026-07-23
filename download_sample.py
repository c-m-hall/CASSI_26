"""
Download the cubes listed in `muse_x_milliquas_sample.fits` and write a local
pairs table (Name, cube_filename) that `mainscript.py` can read directly as
`pairs_path` -- this is the missing link between the notebook's catalogue
(which only lists ESO archive dp_ids, never touches disk) and the pipeline
(which needs real local file paths per QSO name).

Pick fields two ways:
  1. FIELD_INDICES = [0, 5, 12]   -> those exact rows of the sample
  2. FIELD_INDICES = None         -> top N_FIELDS rows by individual-cube exposure

DRY RUN by default: prints how many cubes / how many GB it *would* fetch
(cheap TAP metadata query, no cube bytes) and stops. MUSE WFM cubes are
~5.5 GB each, so MAX_GB caps every real run. Downloads are resumable --
astroquery skips files already in DEST unless you force a re-download.

    python download_data.py     # dry run
    # then set DOWNLOAD = True (and MAX_GB) and re-run to actually fetch
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.table import Table

ESO_TAP = "https://archive.eso.org/tap_obs"



#dataclass to specify size budget, field indices(or number of fields), saved sample fits name, the destination folder for the downloads, and the pairs fits file
@dataclass
class Config:
    sample: str = "muse_x_milliquas_sample.fits"
    dest: str = "./scratch/"
    pairs_out: str = "muse_x_milliquas_pairs_local.fits"

    field_indices: Optional[list] = None   # explicit row indices, or None to rank
    n_fields: int = 3                      # used only when field_indices is None

    max_gb: float = 25.0     # size budget for a real run
    download: bool = False   # set True to actually fetch



#Loads muse_x_milliquas_sample.fits, tags each row with its position (idx) in that file, then picks which fields to actually fetch:
#If we gave explicit field_indices (e.g. [0, 5, 12]), it grabs exactly those rows
#Otherwise it sorts by total individual-cube exposure time and takes the top n_fields

def select_fields(cfg: Config):
    """Chosen sample rows, with an `idx` column giving each row's position
    in the sample file."""
    sample = Table.read(cfg.sample)
    sample["idx"] = np.arange(len(sample))

    if cfg.field_indices is not None:
        n = len(sample)
        bad = [i for i in cfg.field_indices if not -n <= i < n]
        if bad:
            raise IndexError(f"sample has {n} rows; out-of-range indices: {bad}")
        return sample[list(cfg.field_indices)]

    sample.sort("exptime_indiv", reverse=True)
    return sample[: cfg.n_fields]


#each chosen row has a dp_ids column — a ;-separated string of ESO archive dataset IDs (one field can have several cubes/pointings)-> function to unpack that
def dp_ids_for_fields(chosen):
    """De-duplicated, order-preserving list of dp_ids from the chosen rows,
    plus {dp_id: Name} so downloaded files can be traced back to their QSO,
    plus a count of fields with no downloadable individual cube."""
    ids = []
    name_for_dpid = {}
    seen = set()
    for row in chosen:
        text = row["dp_ids"] if row["dp_ids"] not in (None, "--") else ""
        for dp_id in str(text).split(";"):
            if dp_id and dp_id != "--" and dp_id not in seen:
                seen.add(dp_id)
                ids.append(dp_id)
                name_for_dpid[dp_id] = row["Name"]

    n_empty = int(np.sum(chosen["n_cubes_indiv"] == 0))
    return ids, name_for_dpid, n_empty

#tells us how big each file is 
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

#Walks the dp_id list in order, adding up sizes, and stops queuing once the next cube would push the running total over max_gb. 
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


def find_local_file(dp_id, dest):
    """ESO Phase 3 downloads are saved as <dp_id>.fits by default; fall back
    to a substring match in case astroquery/the archive names it slightly
    differently (e.g. a stripped or truncated dp_id)."""
    exact = Path(dest) / f"{dp_id}.fits"
    if exact.exists():
        return exact
    matches = list(Path(dest).glob(f"*{dp_id}*"))
    return matches[0] if matches else None



#Once files are on disk, this looks up each dp_id's actual local filename (find_local_file) 
#builds a (Name, cube_filename) table, and writes it to muse_x_milliquas_pairs_local.fits

def build_local_pairs_table(dp_ids, name_for_dpid, dest):
    """(Name, cube_filename) rows for every dp_id actually found on disk --
    this is the table `get_cube_paths()` reads via `pairs_path`."""
    rows = []
    missing = []
    for dp_id in dp_ids:
        local = find_local_file(dp_id, dest)
        if local is None:
            missing.append(dp_id)
            continue
        rows.append((name_for_dpid[dp_id], local.name))

    table = Table(rows=rows, names=("Name", "cube_filename"))
    return table, missing


def main(cfg: Config = Config()):
    chosen = select_fields(cfg)
    dp_ids, name_for_dpid, n_empty = dp_ids_for_fields(chosen)

    how = (f"sample rows {cfg.field_indices}" if cfg.field_indices is not None
           else f"top-{cfg.n_fields} fields by exposure")
    hours = chosen["exptime_indiv"].sum() / 3600
    print(f"{len(chosen)} field(s) selected ({how}) -> {len(dp_ids)} individual "
          f"cube(s), {hours:.1f} hr of exposure")
    chosen["idx", "Name", "z", "n_cubes_indiv", "exptime_indiv"].pprint(max_width=-1)

    if n_empty:
        print(f"! {n_empty}/{len(chosen)} selected field(s) have no downloadable "
              f"per-OB cube (NFM/offset-sky/combine only) -- nothing to download.")

    sizes = estimated_sizes_gb(dp_ids)
    total_gb = sum(sizes.values())
    queue, queue_gb = apply_budget(dp_ids, sizes, cfg.max_gb)
    print(f"\nfull set ~{total_gb:.1f} GB; MAX_GB={cfg.max_gb} -> queue {len(queue)} of "
          f"{len(dp_ids)} cube(s), ~{queue_gb:.1f} GB, into {cfg.dest!r}")

    if not cfg.download:
        print("\nDRY RUN -- set download=True on the Config (and raise max_gb if "
              "you want more) to fetch the queued cubes.")
        return

    from astroquery.eso import Eso
    Path(cfg.dest).mkdir(parents=True, exist_ok=True)
    eso = Eso()  # public Phase-3 cubes download anonymously; no login() needed

    # continuation=False (default) skips cubes already in DEST -> resumable
    # re-runs. Note astroquery's flag is inverted: continuation=True FORCES a
    # re-download of files already on disk. Always pass destination=cfg.dest,
    # or retrieve_data silently falls back to ~/.astropy/cache/astroquery/Eso.
    eso.retrieve_data(queue, destination=cfg.dest, continuation=False, unzip=True)

    pairs_table, missing = build_local_pairs_table(queue, name_for_dpid, cfg.dest)
    pairs_table.write(cfg.pairs_out, overwrite=True)

    print(f"\nsaved {len(pairs_table)} cube(s) to {cfg.dest!r}")
    if missing:
        print(f"! {len(missing)} queued dp_id(s) not found on disk after download "
              f"(check for a failed/partial fetch): {missing}")
    print(f"wrote local pairs table -> {cfg.pairs_out!r} "
          f"(pass this as pairs_path to mainscript.main())")


if __name__ == "__main__":
    # edit values here (or Config() above) for this run, e.g.:
    #   main(Config(field_indices=[41, 45], download=True, max_gb=50))
    main(Config())
