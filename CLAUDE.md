# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`ctdproc` processes raw CTD (Conductivity-Temperature-Depth) data from Sea-Bird SBE911+ instruments into quality-controlled oceanographic profiles. It is configuration-driven, reading a YAML file to orchestrate a modular processing pipeline.

## Setup & Installation

**Conda** (Python 3.12, recommended):
```bash
conda env create -f environment.yml
conda activate ctdproc
```

**Python venv** (Python 3.11–3.12):
```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

## Running the Processor

```bash
python -m ctdproc proc/config.local.yaml
```

`proc/config.yaml` is the committed template with placeholder paths. Copy it to `proc/config.local.yaml` (gitignored) and set your actual `data_dir` and `output_dir` there. The YAML config controls which casts to process and which pipeline stages to apply.

## Architecture

### Package Layout

- `src/ctdproc/` — main application package
- `src/seabirdscientific/` — embedded, locally-modified copy of Sea-Bird Scientific Community Toolkit v2.7.9

### Data Flow

1. **Load** (`load.py`) — reads SBE911+ `.hex` binary files and `.xmlcon` calibration files
2. **Parse calibration** (`parsers.py`) — `CTDCoefficients` dataclass extracted from the `.xmlcon` XML
3. **Convert** (`load.py:convert_data`) — raw counts → physical units via `seabirdscientific/conversion.py`
4. **Pipeline** (`pipeline.py`) — applies ordered processing stages from YAML config:
   - `align` → temporal alignment of sensors
   - `wild_edit` → outlier rejection
   - `low_filter` → low-pass Butterworth filtering
   - `celltm` → conductivity cell thermal mass correction
   - `loop_edit` → velocity-based quality control (flags reversals)
   - `bin` → pressure-bin averaging
5. **Orchestration** (`main.py`) — reads config, calls stages in order

### seabirdscientific subpackage

This is a vendored, modified fork. Key modules:
- `instrument_data.py` — hex file I/O; supports SBE911+, SBE19+, SBE16+, SBE37, SBE39+
- `conversion.py` — raw count → physical unit conversions (temperature, conductivity, pressure, oxygen, chlorophyll, turbidity, pH)
- `processing.py` — signal processing algorithms (filtering, thermal mass, loop editing, bin averaging, wild editing)
- `cal_coefficients.py` — dataclasses for 20+ sensor calibration types
- `visualization.py` — Plotly-based T-S diagrams and profile plots

### Configuration (`proc/config.yaml`)

```yaml
constants:
  sample_interval: "1/24"              # 24 Hz
  name: "MIXSED2"                      # Hex file base name (e.g. MIXSED2_030.hex)
  ind: [30, 38]                        # Cast index range
  convert: "yes"
  time_ref: "System UTC"
  data_dir: "/path/to/raw"             # Directory containing .hex and .xmlcon files
  output_dir: "proc_output"            # Output directory (created if absent; relative or absolute)

features:                   # Ordered list of pipeline stages to run
  - name: "align"
  - name: "wild_edit"
  - name: "low_filter"
  - name: "celltm"
  - name: "loop_edit"
  - name: "bin"
```

### Key Dependencies

- `gsw` — Gibbs SeaWater equations (TEOS-10) for derived oceanographic quantities
- `scipy` — signal processing (Butterworth filters, optimization)
- `pandas` / `numpy` — data manipulation
- `plotly` — interactive visualization

## Testing

```bash
pip install -r requirements-dev.txt   # adds pytest
pytest tests/
```

Tests use cast 000 fixtures committed under `tests/fixtures/`. To regenerate reference outputs after an intentional pipeline change:

```bash
python tests/generate_fixtures.py
```

Then commit the updated `tests/fixtures/reference/` files alongside the code change.

## Notes

- The `seabirdscientific` subpackage should be modified locally; it is not installed from PyPI.
- Upstream remote is `https://github.com/ocean-uhh/ctdproc.git`; origin is `https://github.com/eleanorfrajka/ctdproc.git`.
- `proc/config.local.yaml` is gitignored — it holds machine-specific paths and is never committed.
