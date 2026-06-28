#!/usr/bin/env python3
"""
Collect all per-subset calibration summary CSVs into one combined table.

Gathers every ``opt_results_summary*.csv`` under AA_opt_out (e.g.
``Rain_Gauge/res/opt_results_summary_rain_gauge.csv``,
``dyn_rho_max/res/opt_results_summary_dyn_rho_max_raing_gauge.csv``, ...),
tags each row with the ``subset`` it came from (the parent-of-res directory),
and writes a single combined CSV.

The ``*_validated.csv`` files are skipped — use ``collect_validated.py`` /
the notebook for those (they carry extra SWE/HNW metric columns).

Usage
-----
    python collect_summaries.py
"""

from pathlib import Path


import pandas as pd

HERE    = Path('calibration/AA_opt_out')          # .../AA_opt_out

print(HERE)


OUT_CSV = HERE / "workflow_helper/all_summaries.csv"


def main() -> None:

    # All summary CSVs, excluding the validated variants.
    csv_files = sorted(
        f for f in HERE.rglob("opt_results_summary*.csv")
        if "validated" not in f.name
    )

    if not csv_files:
        raise SystemExit(f"No opt_results_summary*.csv found under {HERE}")

    print(f"Found {len(csv_files)} summary CSV files:")
    frames = []
    for f in csv_files:
        try:
            d = pd.read_csv(f)
        except pd.errors.EmptyDataError:
            print(f"  {f.relative_to(HERE)}  (empty — skipped)")
            continue
        if d.empty:
            print(f"  {f.relative_to(HERE)}  (no rows — skipped)")
            continue

        d.insert(0, "subset", f.parent.parent.name)   # e.g. Rain_Gauge, Win21
        d["source_file"] = str(f.relative_to(HERE))
        frames.append(d)
        print(f"  {f.relative_to(HERE)}  ({len(d)} rows)")

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUT_CSV, index=False)



    print(f"\n{len(combined)} rows from subsets: {sorted(combined['subset'].unique())}")
    print(f"Written: {OUT_CSV}")


if __name__ == "__main__":
    main()
