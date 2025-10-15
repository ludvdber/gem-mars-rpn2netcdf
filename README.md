# convert_dm_pm_to_nc — README

A tiny CLI tool to convert GEM-Mars RPN files (dm + pm) to NetCDF4, with sane defaults for science usage and batch processing across many subfolders.

## What it does (default behavior)

- Reads paired dm/pm RPN files via fstd2nc.
- Adds a time coordinate: time = reftime + leadtime (if both exist).
- Converts units (per spec):
- - TT: °C → K (+273.15)
- - PX, P0: hPa → Pa (×100)
- - GZ: decametre → m (×10)
- - UU, VV: knots → m s⁻¹ (×0.5144)
- - WW: already in Pa s⁻¹ (unchanged)
- Forces PX and GZ to 103 levels (all thermodynamic + surface).
(Internally selects even indices 0..202 plus surface 203 from the 204-level interleaved grid.)
- Sets GZ(surface) = 0 (height above surface).
- Writes NetCDF-4/HDF5 with lossless compression (zlib level 4, shuffle on).
- Adds simple units metadata.

## Folder layout (expected)

Your project root (e.g. hl-b274) contains:

- hl-b274/ 
  - dm/
    - 000960/
      - hl-b274_dm_000000p_ls000.0000
      - hl-b274_dm_000001p_ls000.0100
      - ...
    - 001920/
      - ...
  - pm/
    - 000960/
      - hl-b274_pm_000000p_ls000.0000
      - hl-b274_pm_000001p_ls000.0100
      - ...
    - 001920/
      - ...
  - netcdf/           # outputs will mirror dm/pm structure here

For each dm/SUBDIR/file, the script looks for the matching pm/SUBDIR/file (same name with _ pm _), and writes to netcdf/SUBDIR/.

Output naming keeps the numbers from the RPN extension:

hl-b274_dm_000000p_ls000.0000   →   netcdf/000960/hl-b274_000000p_ls000_0000.nc

## Install
- Python 3.10+ recommended.
- xarray
- netCDF4
- fstd2nc
- tqdm

tqdm is optional (nice progress bar).
```
pip install xarray netCDF4 fstd2nc tqdm
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
```
TT, PX, GZ, WW, H2O, CO2, O3, CO, T9,
P0,
UU, VV,
MLOC, MALO, MCZ, MH, MTSF, MCO2, MSN,
DVM1, DVM2, DVM3, RWIC
```
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
