"""Pipeline regression tests against committed reference fixtures.

Per-stage tests: load the reference output from the previous stage, apply
the next stage, compare to the reference CSV for that stage.

Integration test: run the full CTD_process on the fixture cast and compare
the written .cnv files to the reference .cnv files.

Fixtures live in tests/fixtures/reference/.  To regenerate them after an
intentional change, run:

    python tests/generate_fixtures.py

then re-commit the updated reference files.
"""

import sys
import warnings
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).parent.parent
FIXTURES = Path(__file__).parent / "fixtures"
RAW_DIR = FIXTURES / "raw"
REF_DIR = FIXTURES / "reference"

sys.path.insert(0, str(ROOT / "src"))

from ctdproc.load import convert_data, load_raw_data
from ctdproc.pipeline import (
    apply_align,
    apply_bin,
    apply_celltm,
    apply_loop_edit,
    apply_low_filter,
    apply_wild_edit,
)

CAST = "MIXSED2_000"
SAMPLE_DT = 1 / 24
BAD_FLAG = np.float64(-9.99e-29)

with open(FIXTURES / "config.yaml") as _f:
    FEATURES = yaml.safe_load(_f)["features"]

STEPS = {s["name"]: s for s in FEATURES}


def load_ref(name: str) -> pd.DataFrame:
    return pd.read_csv(REF_DIR / f"{name}.csv.gz")


def assert_df_close(actual: pd.DataFrame, ref: pd.DataFrame, stage: str) -> None:
    assert list(actual.columns) == list(ref.columns), (
        f"{stage}: column mismatch — got {list(actual.columns)}, expected {list(ref.columns)}"
    )
    assert len(actual) == len(ref), (
        f"{stage}: row count mismatch — got {len(actual)}, expected {len(ref)}"
    )
    for col in ref.columns:
        a = actual[col].to_numpy(dtype=float)
        r = ref[col].to_numpy(dtype=float)
        # NaN positions must match
        np.testing.assert_array_equal(
            np.isnan(a),
            np.isnan(r),
            err_msg=f"{stage} column '{col}': NaN positions differ",
        )
        mask = ~np.isnan(r)
        np.testing.assert_allclose(
            a[mask],
            r[mask],
            rtol=1e-5,
            atol=1e-9,
            err_msg=f"{stage} column '{col}' differs from reference",
        )


# ---------------------------------------------------------------------------
# Module-level fixtures (shared across all per-stage tests)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def stage_convert():
    raw = load_raw_data(RAW_DIR / f"{CAST}.hex")
    return convert_data(raw, RAW_DIR / f"{CAST}.xmlcon", SAMPLE_DT)


@pytest.fixture(scope="module")
def stage_align(stage_convert):
    data = stage_convert.copy()
    data, _ = apply_align(data, STEPS["align"], SAMPLE_DT, BAD_FLAG)
    return data


@pytest.fixture(scope="module")
def stage_wild_edit(stage_align):
    data = stage_align.copy()
    data, _ = apply_wild_edit(data, STEPS["wild_edit"], BAD_FLAG)
    return data


@pytest.fixture(scope="module")
def stage_low_filter(stage_wild_edit):
    data = stage_wild_edit.copy()
    data, _ = apply_low_filter(data, STEPS["low_filter"], SAMPLE_DT)
    return data


@pytest.fixture(scope="module")
def stage_celltm(stage_low_filter):
    data = stage_low_filter.copy()
    data, _ = apply_celltm(data, STEPS["celltm"], SAMPLE_DT)
    return data


@pytest.fixture(scope="module")
def stage_loop_edit(stage_celltm):
    data = stage_celltm.copy()
    data, _ = apply_loop_edit(data, STEPS["loop_edit"], SAMPLE_DT, BAD_FLAG)
    return data


# ---------------------------------------------------------------------------
# Per-stage regression tests
# ---------------------------------------------------------------------------


class TestStageOutputs:
    def test_convert(self, stage_convert):
        assert_df_close(stage_convert, load_ref("stage_00_convert"), "convert_data")

    def test_align(self, stage_align):
        assert_df_close(stage_align, load_ref("stage_01_align"), "align")

    def test_wild_edit(self, stage_wild_edit):
        assert_df_close(stage_wild_edit, load_ref("stage_02_wild_edit"), "wild_edit")

    def test_low_filter(self, stage_low_filter):
        assert_df_close(stage_low_filter, load_ref("stage_03_low_filter"), "low_filter")

    def test_celltm(self, stage_celltm):
        assert_df_close(stage_celltm, load_ref("stage_04_celltm"), "celltm")

    def test_loop_edit(self, stage_loop_edit):
        assert_df_close(stage_loop_edit, load_ref("stage_05_loop_edit"), "loop_edit")

    def test_bin_scans(self, stage_loop_edit):
        data = stage_loop_edit.copy()
        bin_results, _, names = apply_bin(data, STEPS["bin"], BAD_FLAG, CAST)
        label = Path(names[0]).stem
        assert_df_close(
            bin_results[0], load_ref(f"stage_06_bin_{label}"), f"bin/{label}"
        )

    def test_bin_pressure(self, stage_loop_edit):
        data = stage_loop_edit.copy()
        bin_results, _, names = apply_bin(data, STEPS["bin"], BAD_FLAG, CAST)
        label = Path(names[1]).stem
        assert_df_close(
            bin_results[1], load_ref(f"stage_06_bin_{label}"), f"bin/{label}"
        )


# ---------------------------------------------------------------------------
# Integration test — full CTD_process → compare .cnv files
# ---------------------------------------------------------------------------


def _data_lines(cnv_path: Path) -> str:
    """Return the data-only portion of a .cnv file (skip * header lines)."""
    lines = cnv_path.read_text(encoding="utf-8").splitlines(keepends=True)
    data_lines = [l for l in lines if not l.startswith("*")]
    return "".join(data_lines)


def _load_cnv_data(cnv_path: Path) -> pd.DataFrame:
    return pd.read_csv(StringIO(_data_lines(cnv_path)))


class TestIntegration:
    def test_full_pipeline_cnv_outputs(self, tmp_path):
        from ctdproc.main import CTD_process

        # Build a config with concrete paths
        cfg_template = (FIXTURES / "config.yaml").read_text()
        cfg_text = cfg_template.replace("__RAW_DIR__", str(RAW_DIR.resolve()))
        cfg_text = cfg_text.replace("__OUT_DIR__", str(tmp_path.resolve()))

        tmp_cfg = tmp_path / "config.yaml"
        tmp_cfg.write_text(cfg_text)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            CTD_process(str(tmp_cfg))

        expected_files = [
            f"{CAST}_24hz.cnv",
            f"{CAST}_24scans.cnv",
            f"down{CAST}_1dbar.cnv",
        ]

        for fname in expected_files:
            out_path = tmp_path / fname
            ref_path = REF_DIR / fname
            assert out_path.exists(), f"pipeline did not produce {fname}"

            actual = _load_cnv_data(out_path)
            ref = _load_cnv_data(ref_path)
            assert_df_close(actual, ref, fname)
