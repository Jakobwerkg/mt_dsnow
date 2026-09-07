"""
Parse a directory of SMET 1.1 ASCII files into an xarray Dataset.

The output mirrors the Win21 / Mag25 dataset structure:

    Dimensions:   (station: N, time: T)
    Coordinates:
        station   (station)  str
        time      (time)     datetime64[ns]
    Data variables:
        HS        (station, time)  float64   [m]
        SWE       (station, time)  float64   [mm w.e.]
    Coords (1-D per station):
        altitude  (station)        float64   [m a.s.l.]
        easting   (station)        float64   [° lon]
        northing  (station)        float64   [° lat]

Usage
-----
    from read_smet import read_smet_dir

    ds = read_smet_dir("data_18_all/raw_alpsolut")
    ds.to_netcdf("alpsolut_all.nc")
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr


# ── Low-level parser ──────────────────────────────────────────────────────────

def _parse_smet(path: Path) -> dict:
    """
    Read one SMET file and return a dict with:
        station_id, station_name, latitude, longitude, altitude,
        nodata, df  (DataFrame with DatetimeIndex, cols = field names)
    """
    path = Path(path)
    with open(path, encoding="utf-8", errors="replace") as fh:
        raw = fh.read()

    header_block, data_block = raw.split("[DATA]", 1)

    # Parse key = value pairs from the header
    meta = {}
    fields = []
    for line in header_block.splitlines():
        line = line.strip()
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip()
        if key == "fields":
            fields = val.split()
        else:
            meta[key] = val

    nodata = float(meta.get("nodata", -999))

    # Parse data rows — first column is always timestamp
    rows = [l.split() for l in data_block.splitlines() if l.strip()]
    if not rows:
        return None

    times = pd.to_datetime([r[0] for r in rows], utc=True).tz_localize(None)
    data  = np.array([[float(v) for v in r[1:]] for r in rows], dtype=np.float64)

    df = pd.DataFrame(data, index=times, columns=fields[1:])
    df.index.name = "time"
    df.replace(nodata, np.nan, inplace=True)

    return {
        "station_id":   meta.get("station_id", path.stem),
        "station_name": meta.get("station_name", meta.get("station_id", path.stem)),
        "latitude":     float(meta["latitude"]),
        "longitude":    float(meta["longitude"]),
        "altitude":     float(meta["altitude"]),
        "df":           df,
    }


# ── Directory reader ──────────────────────────────────────────────────────────

def read_smet_dir(directory: str | Path,
                  hs_col: str = "HS_meas",
                  swe_col: str = "SWE") -> xr.Dataset:
    """
    Parse all *.smet files in *directory* and return an xarray Dataset.

    Parameters
    ----------
    directory : str or Path
    hs_col    : column name for snow depth   (default "HS_meas")
    swe_col   : column name for SWE          (default "SWE")

    Returns
    -------
    xr.Dataset  with dims (station, time) and variables HS, SWE
    """
    directory = Path(directory)
    smet_files = sorted(directory.glob("*.smet"))
    if not smet_files:
        raise FileNotFoundError(f"No .smet files found in {directory}")

    parsed = []
    for f in smet_files:
        result = _parse_smet(f)
        if result is not None:
            parsed.append(result)

    if not parsed:
        raise ValueError("All .smet files were empty or unparseable.")

    # ── Union time axis ───────────────────────────────────────────────────────
    all_times = pd.DatetimeIndex(
        np.unique(np.concatenate([p["df"].index.values for p in parsed]))
    )

    station_names = [p["station_name"] for p in parsed]
    n_stations    = len(parsed)
    n_times       = len(all_times)

    hs_arr  = np.full((n_stations, n_times), np.nan)
    swe_arr = np.full((n_stations, n_times), np.nan)

    for i, p in enumerate(parsed):
        df  = p["df"]
        idx = all_times.get_indexer(df.index)
        valid = idx >= 0

        if hs_col in df.columns:
            hs_arr[i, idx[valid]] = df[hs_col].values[valid]
        if swe_col in df.columns:
            swe_arr[i, idx[valid]] = df[swe_col].values[valid]

    ds = xr.Dataset(
        {
            "HS":  (["station", "time"], hs_arr),
            "SWE": (["station", "time"], swe_arr),
        },
        coords={
            "station":  station_names,
            "time":     all_times,
            "altitude": ("station", [p["altitude"]  for p in parsed]),
            "easting":  ("station", [p["longitude"] for p in parsed]),
            "northing": ("station", [p["latitude"]  for p in parsed]),
        },
    )

    return ds
