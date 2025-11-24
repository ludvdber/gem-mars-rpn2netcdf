![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![CF Conventions](https://img.shields.io/badge/CF-1.10-orange.svg)
![Status](https://img.shields.io/badge/status-production-brightgreen.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Windows-lightgrey.svg)
![NetCDF](https://img.shields.io/badge/format-NetCDF--4-blue.svg)
![Mars](https://img.shields.io/badge/planet-Mars-red.svg)

**CLI tools for GEM-Mars atmospheric model:** Convert RPN/FSTD files to CF-compliant NetCDF-4 + compute diurnal mean cycles with Mars Year/Ls indexing.


## 📑 Table of Contents

- [convert_dm_pm_to_nc.py](#convert_dm_pm_to_nc--readme)
  - [What it does](#what-it-does-default-behavior)
  - [Folder layout](#folder-layout-expected)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Selecting variables](#selecting-variables)
  - [Output & logging](#output--logging)
  - [Troubleshooting](#troubleshooting)
- [compute_diurnal_mean.py](#computing-diurnal-mean-cycles)
  - [What it does](#what-it-does)
  - [Installation](#installation-1)
  - [Basic Usage](#basic-usage)
  - [Advanced Options](#advanced-options)
  - [Mars Year Lookup](#mars-year-lookup-advanced)
  - [Command Reference](#command-reference)
- [Related Resources](#related-resources)

---

# convert_dm_pm_to_nc — README

A tiny CLI tool to convert GEM-Mars RPN files (dm + pm) to NetCDF4, with sane defaults for science usage and batch processing across many subfolders.

## What it does (default behavior)

- **Paired file processing**: Reads dm + pm RPN files via fstd2nc
- **Time coordinate**: Computes `time = reftime + leadtime`
- **Unit conversions**:
  - TT: °C → K (+273.15)
  - PX, P0: hPa → Pa (×100)
  - GZ: decametre → m (×10)
  - UU, VV: knots → m s⁻¹ (×0.5144)
  - WW: already in Pa s⁻¹ (unchanged)
- **Vertical interpolation**: 
  - Thermodynamic variables (TT, PX, GZ, etc.) interpolated to 103 common altitude levels
  - Momentum variables (UU, VV) interpolated to 102 common altitude levels
  - GZ set to height above surface (GZ_surface = 0)
- **Optimized output**:
  - NetCDF-4 format with lossless compression
  - float32 precision (reduces file size)
  - CF-1.10 compliant metadata
- **Progress tracking**: Optional tqdm progress bar

## Folder layout (expected)

Your project root (e.g., `hl-b274`) should contain:

```
hl-b274/
├── dm/
│   ├── 000960/
│   │   ├── hl-b274_dm_000000p_ls000.0000
│   │   ├── hl-b274_dm_000001p_ls000.0100
│   │   └── ...
│   ├── 001920/
│   └── ...
├── pm/
│   ├── 000960/
│   │   ├── hl-b274_pm_000000p_ls000.0000
│   │   ├── hl-b274_pm_000001p_ls000.0100
│   │   └── ...
│   ├── 001920/
│   └── ...
└── netcdf/           # outputs will be created here
    ├── 000960/
    │   ├── hl-b274_000000p_ls000_0000.nc
    │   ├── hl-b274_000001p_ls000_0100.nc
    │   └── ...
    └── ...
```

For each `dm/SUBDIR/file`, the script looks for the matching `pm/SUBDIR/file` and writes to `netcdf/SUBDIR/`.

**Output naming**: Keeps the numbers from the RPN extension:
```
hl-b274_dm_000000p_ls000.0000  →  netcdf/000960/hl-b274_000000p_ls000_0000.nc
```

## Installation
- Python 3.11+ recommended.
  - xarray
  - netCDF4
  - fstd2nc
  - tqdm
  - numpy

tqdm is optional (nice progress bar).
```
pip install numpy xarray netCDF4 fstd2nc tqdm
```


Linux system libs (if needed)

Debian/Ubuntu: 
```
sudo apt-get install libhdf5-dev libnetcdf-dev
```

## Usage

Define your ROOT (path to the folder that contains dm/, pm/, netcdf/)

- ### Linux

    ```
    ROOT="$HOME/hl-b274"
    cd /path/to/your/repo  #folder where convert_dm_pm_to_nc.py lives
    ```

  - ### Windows PowerShell

      ``` 
      $ROOT = "$env:USERPROFILE\Desktop\hl-b274"
      Set-Location C:\Users\USERPROFILE\PythonFile
     ``` 
  
  1. All subfolders and all files

  - ### Linux

      ```
      python convert_dm_pm_to_nc.py --root "$ROOT" --all
     ```

  - ### Windows PowerShell

      ``` 
      python .\convert_dm_pm_to_nc.py --root "$ROOT" --all
     ``` 
  

2. One specific subfolder (e.g., 000960)

``` 
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960
```

3. A range of subfolders (inclusive, lexical order)
```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir-range 000960 003840
```

4. Within selected subfolders: only one file index (e.g., index 7. index starts from 0)
```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --one 7
```

5. Within selected subfolders: a range of file indices (inclusive)
```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --range 0 3
```

Notes

“file index” refers to the ```*_dm_000007p_*``` number → 7.

If the matching pm file is missing, the script logs SKIP and continues.

With tqdm installed, you get a progress bar and a ```[i/total] WROTE <path> ``` line per output.

## Selecting variables

By default, the script loads exactly the variables requested in the project spec:

### Default variables (28 total)

**3D Thermodynamic** (103 levels):
- `TT` - Temperature (K)
- `PX` - Pressure (Pa)
- `GZ` - Geopotential height (m above surface)
- `WW` - Vertical velocity (Pa/s)
- `H2O`, `CO2`, `O3`, `CO` - Trace gases
- `T9` - Temperature perturbation
- `DVM1`, `DVM2`, `DVM3` - Dust moments
- `RWIC` - Ice particle radius

**3D Momentum** (102 levels):
- `UU` - Eastward wind (m/s)
- `VV` - Northward wind (m/s)

**2D Surface**:
- `P0` - Surface pressure (Pa)
- `MTSF` - Surface temperature (K)
- `MLOC` - Local time
- `MALO`, `MCZ`, `MH`, `MCO2`, `MSN` - Additional surface fields

### Output dimensions

- `time`: 1 (unlimited dimension)
- `lat`: 45 (latitudes, -88° to 88°)
- `lon`: 91 (longitudes, 0° to 360°)
- `altitudeT`: 103 (thermodynamic vertical levels)
- `altitudeM`: 102 (momentum vertical levels)

You can override:
- Keep all variables present:
```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --all-vars
```
- Keep a custom list (comma-separated, no spaces):
```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --vars "TT,PX,GZ,WW,P0,UU,VV"
```
## Output & logging

NetCDF written to ```netcdf/<SUBDIR>/...nc``` with zlib level 4 (lossless) and shuffle.

Console output is minimal:

With tqdm: ```progress bar + [i/total] WROTE <output_path> ```per file.

Without tqdm: ```[i/total] progress lines with ETA and WROTE``` lines.

## Troubleshooting

Missing pm counterpart → file is SKIPped; check ```pm/<SUBDIR>/``` has the same filenames as ```dm/<SUBDIR>/ (with _pm_)```.

Permissions / disk → ensure netcdf/ is writable and has enough space.

Linux build errors for netCDF4/HDF5 → install system dev libraries (see Install).

---

## Computing Diurnal Mean Cycles

After converting RPN files to NetCDF, you can compute mean diurnal cycles (averaged over N Martian sols) using `compute_diurnal_mean.py`.

### What it does

- Groups NetCDF files by hour of day (48 timesteps: 0.0h, 0.5h, ..., 23.5h)
- Averages each hour across multiple Martian days
- Outputs one file with 48 timesteps representing a typical 24-hour cycle
- Uses optimized numba JIT compilation for fast processing
- Preserves all variables and metadata

### Installation

Requires `numba` in addition to previous dependencies:
```bash
pip install numba
```

### Basic Usage

**1. Single directory, N days:**
```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5
```

**2. All directories, 10 days each:**
```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --n-days 10
```

**3. Range of directories:**
```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 003840 \
  --n-days 20
```

### Advanced Options

#### Limit number of means (`--max-means`)

Create only the first N means (useful for testing or limiting output):

```bash
# Create only 5 means of 5 days each (= 25 days total)
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5 \
  --max-means 5
```

**Output:** 5 files (sols 0-4, 5-9, 10-14, 15-19, 20-24)

#### Filter by Ls range (`--ls-range`)

Process only files within a specific solar longitude (Ls) range:

```bash
# Only process Ls 0 to 15°
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --n-days 10 \
  --ls-range 0 15
```

#### Cross-directory mode (`--cross-dirs`)

Allow means to span across multiple directories (treats all selected directories as one continuous dataset):

```bash
# 30-day means that can span across directory boundaries
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 001920 \
  --n-days 30 \
  --cross-dirs
```

**Without `--cross-dirs`:** Each directory processed separately  
**With `--cross-dirs`:** All files treated as continuous timeline

**Output location:** `netcdf_mean/cross_dirs/`

#### Single mean mode (`--single-mean`)

Create ONE single mean from ALL files in the specified range (ignores `--n-days` and `--max-means`):

```bash
# Create one mean Ls 0 to 90°
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --ls-range 0 90 \
  --single-mean
```

**Output:** One file  
**Output location:** `netcdf_mean/single_mean/`

### Combined Examples

**Example 1: First directory complete + 5 days from second**
```bash
# Process first directory completely
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5

# Process only first mean (5 days) from second directory
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 001920 \
  --n-days 5 \
  --max-means 1
```

**Example 2: Cross-directory with filters**
```bash
# 5 means of 5 days, only Ls 0-30°, spanning multiple directories
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 003840 \
  --n-days 5 \
  --max-means 5 \
  --ls-range 0 30 \
  --cross-dirs
```

### Output naming

**Standard mode:**
```
Input:  netcdf/000960/hl-b274_000000p_ls000_0000.nc
Output: netcdf_mean/000960/hl-b274_000000p_ls000_0000_sol000to004_5days_mean.nc
```

**Cross-directory mode:**
```
Output: netcdf_mean/cross_dirs/hl-b274_000000p_ls000_0000_sol000to029_30days_mean_crossdir.nc
```

**Single mean mode:**
```
Output: netcdf_mean/single_mean/hl-b274_000000p_ls000_0000_to_ls090_0000_sol000to500_501days_single_mean.nc
```

The filename includes:
- Starting file index (`000000p`)
- Ls range (`ls000_0000` or `ls000_0000_to_ls090_0000`)
- Sol range (`sol000to004`)
- Number of days averaged (`5days`)
- Mode indicator (`mean`, `crossdir`, or `single_mean`)

### Command Reference

```
usage: compute_diurnal_mean.py --netcdf PATH --output PATH 
                               [--dir DIR | --all | --dir-range START END]
                               [--n-days N] [--max-means N] 
                               [--ls-range MIN MAX] [--mars-year MY] [--ls-start LS]
                               [--lookup FILE] [--cross-dirs] [--single-mean]

Required arguments:
  --netcdf PATH          Root directory containing NetCDF subdirs
  --output PATH          Output directory for mean files
  
Directory selection (required unless using --mars-year):
  --dir DIR              Process single subdirectory
  --all                  Process all subdirectories
  --dir-range START END  Process range of subdirectories

Optional arguments:
  --n-days N            Number of days to average per mean (default: 1)
                        Ignored with --single-mean
  --max-means N         Maximum number of means to compute (default: all)
                        Ignored with --single-mean
  --ls-range MIN MAX    Filter files by Ls range (e.g., --ls-range 0 15)
                        Example: 0 90 for spring, 90 180 for summer
  --cross-dirs          Allow means to span across directories
                        Treats all selected directories as continuous dataset
  --single-mean         Create ONE mean from ALL files in range
                        Useful for seasonal/period averages

Mars Year arguments (advanced):
  --mars-year MY        Mars Year number (e.g., 34, 35)
                        Automatically enables cross-directory search
                        Makes directory selection optional
  --ls-start LS         Starting Ls for Mars Year mode (e.g., 9.61)
  --lookup FILE         Mars Year lookup file (.xlsx or .csv)
                        Required when using --mars-year
```

### Mars Year Lookup (Advanced)

Automatically find starting directories using Mars Year and Ls!

**Simplified usage** - no need to specify `--all` when using `--mars-year`:

```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --lookup Mars_year_Ls_timestep_list.xlsx \
  --mars-year 34 \
  --ls-start 9.61 \
  --n-days 2 \
  --max-means 1
```

**What it does:**
- Reads Excel/CSV lookup table with Mars Year, Ls ranges, and directories
- Automatically finds the starting directory for your Mars Year + Ls
- **Automatically searches all directories** (no need for `--all`)
- Computes means starting from that Ls
- Adds `MY34` to output filename for easy identification

**Requirements:**
```bash
pip install pandas openpyxl
```

**Lookup file format:**
The lookup file should have columns:
- `MY` - Mars Year (e.g., 34, 35)
- `Ls start` - Starting Ls for this range
- `Ls end` - Ending Ls for this range  
- `timestep start` - First timestep in range
- `timestep end` - Last timestep in range
- `directory start` - Directory containing these timesteps

**Mars Year parameters:**
- `--mars-year MY` - Mars Year number (e.g., 34, 35) - **automatically enables cross-directory search**
- `--ls-start LS` - Starting Ls (e.g., 9.61, 30.0)
- `--lookup FILE` - Path to Excel (.xlsx) or CSV lookup file (required with --mars-year)

**Example outputs:**
```
# Without Mars Year
hl-b274_000000p_ls007_1234_sol000to004_5days_mean.nc

# With Mars Year
hl-b274_000000p_ls007_1234_MY34_sol000to004_5days_mean.nc
```

**Notes:**
- `--mars-year` automatically enables cross-directory search (searches all directories for matching timesteps)
- Can combine with `--max-means` to limit output
- Can combine with `--ls-range` for additional filtering
- `--lookup` is required when using `--mars-year`

### Help

For full help including examples:
```bash
python compute_diurnal_mean.py --help
```


---

## Related Resources

### This Project
- **Source Code:** [github.com/ludvdber/gem-mars-rpn2netcdf](https://github.com/ludvdber/gem-mars-rpn2netcdf)
- **Institution:** [Royal Belgian Institute for Space Aeronomy (BIRA-IASB)](https://www.aeronomie.be/)

### GEM-Mars Model & Tools
- **fstd2nc Library:** [github.com/neishm/fstd2nc](https://github.com/neishm/fstd2nc) - Read RPN/FSTD files
- **xarray Documentation:** [docs.xarray.dev](https://docs.xarray.dev/) - NetCDF manipulation in Python

### NetCDF & CF Conventions
- **CF-1.10 Conventions:** [cfconventions.org](http://cfconventions.org/) - Climate and Forecast metadata standard
- **NetCDF Documentation:** [unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf/) - NetCDF format specs

### Visualization Tools
- **Paraview:** [paraview.org](https://www.paraview.org/) - 3D scientific visualization
- **Panoply:** [giss.nasa.gov/tools/panoply](https://www.giss.nasa.gov/tools/panoply/) - NetCDF/HDF viewer from NASA

### Performance Tools
- **Numba:** [numba.pydata.org](https://numba.pydata.org/) - JIT compiler (used in compute_diurnal_mean.py)
- **Dask:** [dask.org](https://dask.org/) - Parallel computing library
