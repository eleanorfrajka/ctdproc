"""Pytest configuration for ctdproc tests.

Sets the matplotlib backend to non-interactive before any test module is
imported, so functions that call plt.xcorr or other plotting routines do not
attempt to open a display window.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")
