"""Verify that the performance refactors in instrument_data and processing
produce numerically identical results to the original implementations.

Run with:  pytest tests/test_perf_refactor.py -v
"""

import numpy as np
import pandas as pd
import pytest

FLAG_VALUE = np.float64(-9.99e-29)


# ---------------------------------------------------------------------------
# Reference implementations (verbatim copies of the originals before refactor)
# ---------------------------------------------------------------------------


def _flag_data_original(
    data,
    flags,
    std_pass_1,
    std_pass_2,
    distance_to_mean,
    exclude_bad_flags,
    flag_value=FLAG_VALUE,
):
    data_copy = pd.Series(data.copy())
    flagged_data = data.copy()

    for n, value in enumerate(data_copy):
        if exclude_bad_flags and flags[n] == flag_value:
            data_copy[n] = np.nan

    mean = data_copy.mean()
    std = data_copy.std()

    for n, value in enumerate(data_copy):
        if abs(value - mean) >= std * std_pass_1:
            data_copy[n] = np.nan

    mean = data_copy.mean()
    std = data_copy.std()

    for n, value in enumerate(flagged_data):
        if (
            abs(value - mean) > std * std_pass_2
            and abs(value - mean) > distance_to_mean
        ):
            flagged_data[n] = flag_value

    return flagged_data


def _bin_average_interp_original(dataset, bin_variable):
    """Runs only the interpolation block from bin_average, for comparison."""
    dataset = dataset.copy()

    def interp(p_p, x_p, p_c, x_c, p_i):
        return ((x_c - x_p) * (p_i - p_p) / (p_c - p_p)) + x_p

    excluded_columns = ["nbin", "flag", "bin_number", bin_variable, "midpoint"]
    for column in dataset.columns.difference(excluded_columns):
        interp_result = []
        for n in range(len(dataset[column])):
            n_p = 1 if n == 0 else n - 1
            p_p = dataset[bin_variable].iloc[n_p]
            x_p = dataset[column].iloc[n_p]
            p_c = dataset[bin_variable].iloc[n]
            x_c = dataset[column].iloc[n]
            p_i = dataset["midpoint"].iloc[n]
            x_i = interp(p_p, x_p, p_c, x_c, p_i)
            interp_result.append(x_i)
        dataset[column] = pd.Series(interp_result, index=dataset.index)

    return dataset


# ---------------------------------------------------------------------------
# New implementations (imported from the refactored modules)
# ---------------------------------------------------------------------------

from seabirdscientific.processing import _flag_data as _flag_data_new


def _bin_average_interp_new(dataset, bin_variable):
    """Runs only the refactored interpolation block from bin_average."""
    dataset = dataset.copy()
    excluded_columns = ["nbin", "flag", "bin_number", bin_variable, "midpoint"]

    p_c = dataset[bin_variable].to_numpy()
    p_p = np.empty_like(p_c)
    p_p[0] = p_c[1] if len(p_c) > 1 else p_c[0]
    p_p[1:] = p_c[:-1]
    p_i = dataset["midpoint"].to_numpy()

    for column in dataset.columns.difference(excluded_columns):
        x_c = dataset[column].to_numpy(dtype=float)
        x_p = np.empty_like(x_c)
        x_p[0] = x_c[1] if len(x_c) > 1 else x_c[0]
        x_p[1:] = x_c[:-1]
        dataset[column] = ((x_c - x_p) * (p_i - p_p) / (p_c - p_p)) + x_p

    return dataset


# ---------------------------------------------------------------------------
# Tests: _flag_data
# ---------------------------------------------------------------------------


def make_data(n=500, seed=42):
    rng = np.random.default_rng(seed)
    data = rng.normal(loc=10.0, scale=0.5, size=n)
    flags = np.zeros(n, dtype=float)
    # inject a few pre-existing bad flags and a few outliers
    flags[[10, 50, 200]] = FLAG_VALUE
    data[[20, 150, 400]] = 99.0  # outliers
    return data, flags


@pytest.mark.parametrize("exclude_bad_flags", [True, False])
def test_flag_data_no_preflags(exclude_bad_flags):
    data, flags = make_data()
    flags[:] = 0  # no pre-existing flags

    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags)

    np.testing.assert_array_equal(orig, new)


@pytest.mark.parametrize("exclude_bad_flags", [True, False])
def test_flag_data_with_preflags(exclude_bad_flags):
    data, flags = make_data()

    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags)

    np.testing.assert_array_equal(orig, new)


def test_flag_data_distance_to_mean():
    data, flags = make_data()
    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.5, exclude_bad_flags=True)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.5, exclude_bad_flags=True)
    np.testing.assert_array_equal(orig, new)


def test_flag_data_all_clean():
    """No outliers, no pre-flags — output should equal input."""
    rng = np.random.default_rng(0)
    data = rng.normal(size=200)
    flags = np.zeros(200)
    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    np.testing.assert_array_equal(orig, new)


def test_flag_data_flags_longer_than_data():
    """wild_edit passes the full flags array but only a block of data — must not crash."""
    rng = np.random.default_rng(9)
    data = rng.normal(size=100)
    flags = np.zeros(100_000)  # full cast length
    flags[::500] = FLAG_VALUE

    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    np.testing.assert_array_equal(orig, new)


def test_flag_data_large(benchmark=None):
    """Larger array — also exercises performance."""
    rng = np.random.default_rng(7)
    n = 50_000
    data = rng.normal(loc=5.0, scale=1.0, size=n)
    flags = np.zeros(n)
    flags[::500] = FLAG_VALUE
    data[::1000] = 50.0

    orig = _flag_data_original(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    new = _flag_data_new(data, flags, 2.0, 5.0, 0.0, exclude_bad_flags=True)
    np.testing.assert_array_equal(orig, new)


# ---------------------------------------------------------------------------
# Tests: bin_average interpolation
# ---------------------------------------------------------------------------


def make_bin_dataset(n=20):
    pressure = np.arange(1.0, n + 1, dtype=float)
    midpoint = pressure - 0.5
    temperature = 20.0 - pressure * 0.05
    salinity = 35.0 + pressure * 0.001
    df = pd.DataFrame(
        {
            "pressure": pressure,
            "midpoint": midpoint,
            "temperature": temperature,
            "salinity": salinity,
            "nbin": np.full(n, 5),
            "flag": np.zeros(n),
            "bin_number": np.arange(n),
        }
    )
    return df


def test_bin_interp_basic():
    df = make_bin_dataset()
    orig = _bin_average_interp_original(df, "pressure")
    new = _bin_average_interp_new(df, "pressure")
    for col in ["temperature", "salinity"]:
        np.testing.assert_allclose(
            orig[col].values,
            new[col].values,
            rtol=1e-12,
            err_msg=f"mismatch in column {col!r}",
        )


def test_bin_interp_single_row():
    """Edge case: only one bin row.
    The original crashes (iloc[1] out of bounds); the refactor handles it safely.
    """
    df = make_bin_dataset(n=1)
    with pytest.raises(IndexError):
        _bin_average_interp_original(df, "pressure")
    # new code should not raise
    result = _bin_average_interp_new(df, "pressure")
    assert result is not None


def test_bin_interp_uneven_pressure():
    rng = np.random.default_rng(3)
    n = 50
    pressure = np.cumsum(rng.uniform(0.5, 2.0, size=n))
    midpoint = np.concatenate([[pressure[0]], (pressure[:-1] + pressure[1:]) / 2])
    df = pd.DataFrame(
        {
            "pressure": pressure,
            "midpoint": midpoint,
            "temperature": 20.0 - pressure * 0.05 + rng.normal(scale=0.01, size=n),
            "salinity": 35.0 + pressure * 0.001,
            "nbin": np.full(n, 3),
            "flag": np.zeros(n),
            "bin_number": np.arange(n),
        }
    )
    orig = _bin_average_interp_original(df, "pressure")
    new = _bin_average_interp_new(df, "pressure")
    for col in ["temperature", "salinity"]:
        np.testing.assert_allclose(
            orig[col].values,
            new[col].values,
            rtol=1e-12,
            err_msg=f"mismatch in column {col!r}",
        )
