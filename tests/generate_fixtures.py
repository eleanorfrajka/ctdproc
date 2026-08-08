"""Generate reference CSV fixtures for pipeline regression tests.

Run this once (from the repo root) to create tests/fixtures/reference/*.csv:

    python tests/generate_fixtures.py

Commit the resulting CSVs.  Re-run only when a deliberate change to the
pipeline output is intended — the regression tests compare against these files.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from ctdproc.load import load_raw_data, convert_data
from ctdproc.pipeline import (
    apply_align,
    apply_low_filter,
    apply_wild_edit,
    apply_celltm,
    apply_loop_edit,
    apply_bin,
)

RAW_DIR = ROOT / "tests" / "fixtures" / "raw"
REF_DIR = ROOT / "tests" / "fixtures" / "reference"
REF_DIR.mkdir(exist_ok=True)

CAST = "MIXSED2_000"
SAMPLE_DT = 1 / 24
BAD_FLAG = np.float64(-9.99e-29)

CONFIG_FEATURES = [
    {
        "name": "align",
        "variables": ["conductivity", "conductivity2"],
        "offset": [0.073, 0.073],
    },
    {
        "name": "wild_edit",
        "variables": [
            "temperature",
            "temperature2",
            "conductivity",
            "conductivity2",
            "pressure",
        ],
        "std_pass_1": [2.0, 2.0, 2.0, 2.0, 2.0],
        "std_pass_2": [5.0, 5.0, 5.0, 5.0, 5.0],
        "scans_per_block": [100, 100, 100, 100, 100],
        "distance_to_mean": [0.0, 0.0, 0.0, 0.0, 0.0],
    },
    {
        "name": "low_filter",
        "variables": ["pressure", "temperature", "conductivity"],
        "time_constant": [0.03, 0.03, 0.03],
    },
    {
        "name": "celltm",
        "variables": ["conductivity", "conductivity2"],
        "amplitude": [0.03, 0.03],
        "time_constant": [7.0, 7.0],
    },
    {
        "name": "loop_edit",
        "min_velocity": [0.25],
        "window_size": [3],
        "mean_speed_percent": [20],
        "remove_surface_soak": [True],
        "min_soak_depth": [5],
        "max_soak_depth": [20],
        "use_deck_pressure_offset": [False],
    },
    {
        "name": "bin",
        "bin_type": ["scans", "pressure"],
        "bin_size": [24.0, 1.0],
        "type_profile": ["BOTH", "DOWN"],
    },
]


def save(df: pd.DataFrame, name: str) -> None:
    path = REF_DIR / f"{name}.csv.gz"
    df.to_csv(path, index=False, float_format="%.6g", compression="gzip")
    print(
        f"  saved {path.name}  ({path.stat().st_size // 1024}KB, {len(df)} rows × {len(df.columns)} cols)"
    )


steps = {s["name"]: s for s in CONFIG_FEATURES}

print(f"Loading {CAST} ...")
raw = load_raw_data(RAW_DIR / f"{CAST}.hex")
data = convert_data(raw, RAW_DIR / f"{CAST}.xmlcon", SAMPLE_DT)
save(data, "stage_00_convert")

print("align ...")
data, _ = apply_align(data, steps["align"], SAMPLE_DT, BAD_FLAG)
save(data, "stage_01_align")

print("wild_edit ...")
data, _ = apply_wild_edit(data, steps["wild_edit"], BAD_FLAG)
save(data, "stage_02_wild_edit")

print("low_filter ...")
data, _ = apply_low_filter(data, steps["low_filter"], SAMPLE_DT)
save(data, "stage_03_low_filter")

print("celltm ...")
data, _ = apply_celltm(data, steps["celltm"], SAMPLE_DT)
save(data, "stage_04_celltm")

print("loop_edit ...")
data, _ = apply_loop_edit(data, steps["loop_edit"], SAMPLE_DT, BAD_FLAG)
save(data, "stage_05_loop_edit")

print("bin ...")
bin_results, _, output_names = apply_bin(data, steps["bin"], BAD_FLAG, CAST)
# output_names: ['MIXSED2_000_24scans.cnv', 'downMIXSED2_000_1dbar.cnv']
for bin_df, name in zip(bin_results, output_names):
    label = Path(name).stem  # e.g. MIXSED2_000_24scans
    save(bin_df, f"stage_06_bin_{label}")

print("Done.")
