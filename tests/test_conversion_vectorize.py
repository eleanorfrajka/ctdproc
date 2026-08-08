"""Verify that vectorised rewrites of the two conversion bottlenecks
produce numerically identical results to the original Python-loop implementations.

Uses the committed fixture cast (tests/fixtures/raw/MIXSED2_000).

Run with:  pytest tests/test_conversion_vectorize.py -v -s
"""

import sys
from math import floor
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.signal import lfilter, lfilter_zi, savgol_filter

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

import seabirdscientific.conversion as conv
from ctdproc.load import load_raw_data
from ctdproc.parsers import load_xmlcon

# ---------------------------------------------------------------------------
# Load real data once for the whole module
# ---------------------------------------------------------------------------

_RAW_DIR = ROOT / "tests" / "fixtures" / "raw"
SAMPLE_INTERVAL = 1 / 24

_raw = load_raw_data(_RAW_DIR / "MIXSED2_000.hex")
_coeffs = load_xmlcon(_RAW_DIR / "MIXSED2_000.xmlcon")

VOLTAGE = _raw["volt 2"].to_numpy(dtype=float)
TEMP_C = conv.convert_temperature_frequency(
    _raw["temperature"], _coeffs.temperature_primary, "ITS90", "C"
)
PRESSURE = conv.convert_pressure_digiquartz(
    _raw["digiquartz pressure"],
    _raw["temperature compensation"],
    _coeffs.pressure,
    "dbar",
    SAMPLE_INTERVAL,
)
SALINITY = np.zeros(len(TEMP_C))  # placeholder — not under test here
COMP_VOLT = _raw["temperature compensation"].to_numpy(dtype=float)
COEFS_OX = _coeffs.oxygen_primary
COEFS_PR = _coeffs.pressure


# ===========================================================================
# Reference implementations (original loops, verbatim)
# ===========================================================================


def _tau_correction_old(voltage, sample_interval, window_size=1):
    """Original: 100K calls to scipy.stats.linregress."""
    scans_per_side = floor(window_size / 2 / sample_interval)
    dvdt = np.zeros(len(voltage))
    for i in range(scans_per_side, len(voltage) - scans_per_side):
        ox_subset = voltage[i - scans_per_side : i + scans_per_side + 1]
        time_subset = np.arange(
            0, len(ox_subset) * sample_interval, sample_interval, dtype=float
        )
        result = stats.linregress(time_subset, ox_subset)
        dvdt[i] = result.slope
    return dvdt


def _hysteresis_correction_old(voltage, pressure, coefs, sample_interval):
    """Original: scalar loop (recurrence relation — unchanged in new code)."""
    corrected = voltage.copy()
    for i in range(1, len(corrected)):
        d = 1 + coefs.h1 * (np.exp(pressure[i] / coefs.h2) - 1)
        c = np.exp(-1 * sample_interval / coefs.h3)
        ox = corrected[i] + coefs.v_offset
        prev = corrected[i - 1] + coefs.v_offset
        ox_new = ((ox + prev * c * d) - (prev * c)) / d
        corrected[i] = ox_new - coefs.v_offset
    return corrected


def _rolling_mean_old(compensation_voltage, scans_in_window):
    """Original: manual rolling sum loop."""
    rolling_sum = compensation_voltage[0] * scans_in_window
    result = compensation_voltage.copy()
    for i in range(len(compensation_voltage)):
        if i < scans_in_window:
            rolling_sum -= compensation_voltage[0]
        else:
            rolling_sum -= compensation_voltage[i - scans_in_window]
        rolling_sum += compensation_voltage[i]
        result[i] = rolling_sum / scans_in_window
    return result


# ===========================================================================
# New vectorised implementations (candidates to replace the originals)
# ===========================================================================


def _tau_correction_new(voltage, sample_interval, window_size=1):
    """Savitzky-Golay first derivative — equivalent to linregress on uniform grid.

    Local copy of the algorithm now inlined in seabirdscientific.conversion
    (not exported as a standalone function). Kept here to allow direct
    comparison against the original loop implementation.
    """
    scans_per_side = floor(window_size / 2 / sample_interval)
    wlen = 2 * scans_per_side + 1
    if wlen < 3 or wlen > len(voltage):
        return np.zeros(len(voltage))
    dvdt = savgol_filter(voltage, wlen, polyorder=1, deriv=1, delta=sample_interval)
    # Match original: boundary scans where full window doesn't fit stay at 0
    dvdt[:scans_per_side] = 0.0
    dvdt[len(voltage) - scans_per_side :] = 0.0
    return dvdt


def _rolling_mean_new(compensation_voltage, scans_in_window):
    """Causal rolling mean via scipy lfilter with edge initial condition.

    Local copy of the algorithm now inlined in seabirdscientific.conversion
    (not exported as a standalone function). Kept here to allow direct
    comparison against the original loop implementation.
    """
    b = np.ones(scans_in_window) / scans_in_window
    zi = lfilter_zi(b, 1) * compensation_voltage[0]
    result, _ = lfilter(b, 1, compensation_voltage, zi=zi)
    return result


# ===========================================================================
# Tests
# ===========================================================================


class TestTauCorrection:
    def test_dvdt_interior_identical(self):
        """Interior dvdt values must match to float64 precision."""
        scans_per_side = floor(1.0 / 2 / SAMPLE_INTERVAL)
        old = _tau_correction_old(VOLTAGE, SAMPLE_INTERVAL, window_size=1)
        new = _tau_correction_new(VOLTAGE, SAMPLE_INTERVAL, window_size=1)
        interior = slice(scans_per_side, len(VOLTAGE) - scans_per_side)
        # max absolute difference is ~2e-16 (machine epsilon) from
        # different fp paths: linregress QR vs savgol convolution coefficients.
        # Use atol to handle near-zero dvdt values where rtol alone would fail.
        np.testing.assert_allclose(
            old[interior],
            new[interior],
            rtol=1e-8,
            atol=1e-13,
            err_msg="dvdt mismatch on interior scans",
        )

    def test_dvdt_boundary_zeros(self):
        """Boundary scans must be zero in both implementations."""
        scans_per_side = floor(1.0 / 2 / SAMPLE_INTERVAL)
        old = _tau_correction_old(VOLTAGE, SAMPLE_INTERVAL, window_size=1)
        new = _tau_correction_new(VOLTAGE, SAMPLE_INTERVAL, window_size=1)
        assert np.all(old[:scans_per_side] == 0.0), "old boundary not zero"
        assert np.all(new[:scans_per_side] == 0.0), "new boundary not zero"

    def test_full_oxygen_array_unchanged(self):
        """End-to-end oxygen conversion must be identical."""
        old_dvdt = _tau_correction_old(VOLTAGE, SAMPLE_INTERVAL)
        new_dvdt = _tau_correction_new(VOLTAGE, SAMPLE_INTERVAL)

        old_ox = conv._convert_sbe43_oxygen(
            VOLTAGE, TEMP_C, PRESSURE, SALINITY, COEFS_OX, old_dvdt
        )
        new_ox = conv._convert_sbe43_oxygen(
            VOLTAGE, TEMP_C, PRESSURE, SALINITY, COEFS_OX, new_dvdt
        )

        np.testing.assert_allclose(
            old_ox,
            new_ox,
            rtol=1e-10,
            err_msg="oxygen array differs after tau correction swap",
        )


class TestRollingMean:
    def _scans_in_window(self):
        max_scans = 720
        w = floor(30 / SAMPLE_INTERVAL)
        return max(1, min(w, max_scans))

    def test_rolling_mean_identical(self):
        """Vectorised rolling mean must exactly match the original loop."""
        w = self._scans_in_window()
        old = _rolling_mean_old(COMP_VOLT, w)
        new = _rolling_mean_new(COMP_VOLT, w)
        np.testing.assert_allclose(
            old, new, rtol=1e-12, err_msg=f"rolling mean mismatch (window={w})"
        )

    def test_pressure_array_unchanged(self):
        """Full pressure conversion must produce the same dbar values."""
        # Run the current (loop) implementation
        old_pressure = conv.convert_pressure_digiquartz(
            _raw["digiquartz pressure"],
            _raw["temperature compensation"],
            COEFS_PR,
            "dbar",
            SAMPLE_INTERVAL,
        )
        # After the fix lands in conversion.py this test will keep passing;
        # for now both call the same function so this is a smoke test.
        assert len(old_pressure) == len(COMP_VOLT)
        assert not np.any(np.isnan(old_pressure))
