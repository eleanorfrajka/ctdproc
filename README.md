# ctdproc-MIXSED2

This repository is designed for processing shipboard Sea-Bird SBE 911plus CTD data. It is currently being developed to support the **MIXSED project** campaign workflow. The long-term goal is to evolve this codebase into a fully operational, generic CTD processing repository and eventually merge it into the main `ctdam` repository.

## Repository Structure

* **`src/`**: Contains the core processing packages and execution modules.
  * **`src/ctdproc/load.py`**: **[NOTE]** This module is currently hardcoded and optimized exclusively to handle data for the **MIXSED2** campaign.
* **`seabirdscientific/`**: **[NOTE]** Contains the official *Sea-Bird Scientific Community Toolkit v2.7.9* with **local modifications** to adapt it to specific project campaign processing requirements.

---

## Installation & Setup

You can set up the environment and install the package using either `pip` or `conda`/`mamba` depending on your requirements.

### Option 1: Using pip (requirements.txt)
```bash
pip install -r requirements.txt
pip install -e .
```

### Option 2: Using Conda/Mamba (environment.yml)
```bash
conda env create -f environment.yml
conda activate ctdproc_env
pip install -e .
```

---

## How to Run the Processing Pipeline

To run the processing execution module, you need to create a dedicated working directory called `proc/` inside the root folder.

### 1. Data Directory Layout
The execution command expects your raw data files to be located one level above the working directory in a `data/` folder:

```text
ctdproc-MIXSED2/
├── data/               # Raw input data (.hex files)
├── proc/               # Your execution directory
│   └── config.yaml    # Campaign configuration file
├── src/
├── pyproject.toml
└── README.md
```

### 2. Execution Command
Navigate into your `proc/` directory and run the module by passing the configuration file:

```bash
cd proc
python -m ctdproc config.yaml
```

---

## Configuration File (`config.yaml`)

The processing sequence is entirely driven by a `config.yaml` file located in your execution directory. Here is an explanation of its components:

```yaml
constants:
  sample_interval: 1/24
  name: "MIXSED2"
  ind: [30, 38]        # Loop iteration indices. Expects filenames formatted as name_030.hex to name_038.hex
  convert: 'yes'       # Set to 'yes' to convert raw hex data to nominal 24 Hz resolution (.cnv)
  time_ref: 'System UTC'

features:
# Only the un-commented features listed here are active during execution.
# - name: "align"
#   variables: ["conductivity","conductivity2"]
#   offset: [0.073,0.073]
# - name: "wild_edit"
#   variables: ["temperature","temperature2","conductivity","conductivity2","pressure"]
#   std_pass_1: [2.0,2.0,2.0,2.0,2.0]
#   std_pass_2: [5.0,5.0,5.0,5.0,5.0]
#   scans_per_block: [100,100,100,100,100]
#   distance_to_mean: [0.0,0.0,0.0,0.0,0.0]
# - name: "low_filter"
#   variables:  ["pressure","temperature","conductivity"]
#   time_constant: [0.03,0.03,0.03]
# - name: "celltm"
#   variables: ["conductivity","conductivity2"]
#   amplitude: [0.03,0.03]
#   time_constant: [7.0,7.0]
# - name: "loop_edit"
#   min_velocity: [0.25]
#   window_size: [3]
#   mean_speed_percent: [20]
#   remove_surface_soak: [True]
#   min_soak_depth: [5]
#   max_soak_depth: [20]
#   use_deck_pressure_offset: [False]
# - name: "bin"
#   bin_type: ["scans","pressure"]
#   bin_size: [24.0,1.0]
#   type_profile: ["BOTH",'DOWN']
```

### Key Parameter Definitions:
* **`ind`**: Defines the file iteration range. The pipeline loops through data indexes (e.g., from cast `030` to cast `038`).
* **`convert`**: When enabled (`'yes'`), it translates raw engineering files into nominal 24 Hz `.cnv` profiles using internal calibrations.
* **`features`**: A list of data processing steps. The pipeline will dynamically search for and execute only the modules explicitly listed and uncommented in this section.


---

## Diagnostic Tools

### CTD Alignment Finder (`ctdalign_finder`)
Computes temporal lags between sensors (e.g., temperature and conductivity) across profiles to optimize your pipeline coefficients.

#### 1. Configuration (`proc/config_align.yaml`)
Create an independent configuration file using this layout:
```yaml
constants:
  sample_interval: 1/24
  name: "MIXSED2"
  ind: [30, 38]

alignment_settings:
  variable_lead: ["temperature"]
  variable_lag: ["conductivity"]
  max_lag: [15]
  pressure_upper: [30.0]
  pressure_lower: [500.0]
```

#### 2. Execution
Run the submodule from your `proc/` directory:
```bash
cd proc
python -m ctdproc.ctdalign_finder config_align.yaml
```

#### 3. Output
Generates **`proc/alignment_report.csv`** containing:
* **Detailed Log**: Individual lag offsets (`lag_dn`, `lag_up`) per cast.
* **Campaign Averages**: Global mean values to paste into your production `config.yaml`.
