#!/usr/bin/env python
"""compute_diurnal_mean.py - Cycle diurnal moyen GEM-Mars avec numba"""

import argparse
import sys
import re
import time
import warnings
from pathlib import Path
import numpy as np
import xarray as xr
from numba import jit
from tqdm import tqdm

TIMESTEPS_PER_DAY = 48


@jit(nopython=True, cache=True)
def fast_mean(data, indices):
    """Moyenne optimisée avec numba."""
    n = len(indices)
    if n == 0:
        return data[0] * 0
    result = np.zeros(data.shape[1:], dtype=data.dtype)
    for idx in indices:
        result += data[idx]
    return result / n


def extract_file_index(filename: str) -> int | None:
    match = re.search(r'_(\d{6})p_', filename)
    return int(match.group(1)) if match else None


def compute_daily_mean(files: list[Path], output_path: Path):
    start_time = time.time()
    n_days = len(files) // TIMESTEPS_PER_DAY
    print(f"Processing {len(files)} files ({n_days} days)")

    first_ds = xr.open_dataset(str(files[0]))
    altitude_T = first_ds['altitudeT'].copy() if 'altitudeT' in first_ds else None
    altitude_M = first_ds['altitudeM'].copy() if 'altitudeM' in first_ds else None

    t0 = time.time()
    datasets = []
    for f in tqdm(files, desc="Loading", unit="file"):
        ds = xr.open_dataset(str(f), chunks=None)
        ds = ds.drop_vars(['altitudeT', 'altitudeM'], errors='ignore')
        datasets.append(ds)
    print(f"Loaded in {time.time() - t0:.1f}s")

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        ds_all = xr.concat(datasets, dim='time', join='inner', combine_attrs='override')

    if altitude_T is not None:
        ds_all = ds_all.assign_coords(altitudeT=altitude_T)
    if altitude_M is not None:
        ds_all = ds_all.assign_coords(altitudeM=altitude_M)
    print(f"Concat in {time.time() - t0:.1f}s")

    t0 = time.time()
    hours = (ds_all.time.dt.hour + ds_all.time.dt.minute / 60.0).values
    unique_hours = np.unique(hours)

    mean_data_vars = {}
    for var_idx, var_name in enumerate(ds_all.data_vars, 1):
        var = ds_all[var_name]
        if 'time' not in var.dims:
            mean_data_vars[var_name] = var
            continue

        time_axis = var.dims.index('time')
        out_shape = list(var.shape)
        out_shape[time_axis] = len(unique_hours)
        out_data = np.zeros(out_shape, dtype=var.dtype)

        for i, hour in enumerate(unique_hours):
            indices = np.where(np.abs(hours - hour) < 0.01)[0]

            # Utiliser np.take pour indexing correct sur n'importe quel axe
            selected = np.take(var.values, indices, axis=time_axis)
            out_slice = selected.mean(axis=time_axis)

            # Insérer au bon endroit
            if time_axis == 0:
                out_data[i] = out_slice
            elif time_axis == 1:
                out_data[:, i] = out_slice
            elif time_axis == 2:
                out_data[:, :, i] = out_slice
            elif time_axis == 3:
                out_data[:, :, :, i] = out_slice

        new_coords = {}
        for dim in var.dims:
            if dim == 'time':
                new_coords[dim] = unique_hours
            elif dim in var.coords:
                new_coords[dim] = var.coords[dim]

        mean_data_vars[var_name] = xr.DataArray(out_data, dims=var.dims, coords=new_coords, attrs=var.attrs)

    diurnal_mean = xr.Dataset(mean_data_vars, attrs=ds_all.attrs)
    print(f"Averaging in {time.time() - t0:.1f}s")

    if altitude_T is not None:
        diurnal_mean = diurnal_mean.assign_coords(altitudeT=altitude_T)
    if altitude_M is not None:
        diurnal_mean = diurnal_mean.assign_coords(altitudeM=altitude_M)

    time_values = np.arange(len(unique_hours)) * 0.5
    diurnal_mean = diurnal_mean.assign_coords(time=time_values)
    diurnal_mean['time'].attrs = {
        'long_name': 'Hour of day (local solar time)',
        'units': 'hours',
        'description': f'Diurnal cycle averaged over {n_days} Martian sols',
    }

    first_index = extract_file_index(files[0].name)
    global_day_start = first_index // 48 if first_index else 0
    global_day_end = global_day_start + n_days - 1

    diurnal_mean.attrs.update({
        'title': 'GEM-Mars Diurnal Mean Cycle',
        'description': f'Mean diurnal cycle averaged over {n_days} Martian sols',
        'n_files': len(files),
        'n_days': n_days,
        'day_start': global_day_start,
        'day_end': global_day_end,
        'Conventions': 'CF-1.10',
    })

    for var in diurnal_mean.data_vars:
        if 'time' in diurnal_mean[var].dims:
            diurnal_mean[var].attrs['cell_methods'] = 'time: mean'

    t0 = time.time()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    from netCDF4 import Dataset as NC4Dataset

    with NC4Dataset(str(output_path), 'w', format='NETCDF4') as ncout:
        ncout.createDimension('time', len(unique_hours))
        if altitude_T is not None:
            ncout.createDimension('altitudeT', len(altitude_T))
        if altitude_M is not None:
            ncout.createDimension('altitudeM', len(altitude_M))
        ncout.createDimension('lat', diurnal_mean.sizes['lat'])
        ncout.createDimension('lon', diurnal_mean.sizes['lon'])

        time_var = ncout.createVariable('time', 'f4', ('time',))
        time_var[:] = time_values
        time_var.units = 'hours'
        time_var.long_name = 'Hour of day (local solar time)'

        for coord_name in ['lat', 'lon']:
            coord_var = ncout.createVariable(coord_name, 'f4', (coord_name,))
            coord_var[:] = diurnal_mean[coord_name].values
            for attr, val in diurnal_mean[coord_name].attrs.items():
                setattr(coord_var, attr, val)

        if altitude_T is not None:
            alt_var = ncout.createVariable('altitudeT', 'f4', ('altitudeT',))
            alt_var[:] = altitude_T.values
            for attr, val in altitude_T.attrs.items():
                setattr(alt_var, attr, val)

        if altitude_M is not None:
            alt_var = ncout.createVariable('altitudeM', 'f4', ('altitudeM',))
            alt_var[:] = altitude_M.values
            for attr, val in altitude_M.attrs.items():
                setattr(alt_var, attr, val)

        for var_name, var_data in diurnal_mean.data_vars.items():
            dims = var_data.dims

            if 'time' in dims:
                chunks = []
                for d in dims:
                    if d == 'time':
                        chunks.append(1)
                    else:
                        chunks.append(diurnal_mean.sizes[d])
                chunksizes = tuple(chunks)
            else:
                chunksizes = None

            var_out = ncout.createVariable(var_name, 'f4', dims, zlib=True, complevel=9, shuffle=True,
                                           chunksizes=chunksizes)
            var_out[:] = var_data.values

            for attr, val in var_data.attrs.items():
                setattr(var_out, attr, val)

        for attr, val in diurnal_mean.attrs.items():
            setattr(ncout, attr, val)

    file_size_mb = output_path.stat().st_size / (1024 ** 2)
    print(f"Saved in {time.time() - t0:.1f}s ({file_size_mb:.1f} MB)")

    for ds in datasets:
        ds.close()
    ds_all.close()
    first_ds.close()

    import gc
    gc.collect()

    elapsed = time.time() - start_time
    return elapsed


def process_subdirectory(netcdf_root: Path, subdir_name: str, output_root: Path, n_days: int = 1):
    input_dir = netcdf_root / subdir_name
    output_dir = output_root / subdir_name

    if not input_dir.exists():
        print(f"ERROR: {input_dir} not found")
        return

    all_files = sorted(input_dir.glob("*.nc"), key=lambda f: extract_file_index(f.name) or 0)

    if not all_files:
        print(f"ERROR: No NetCDF files in {input_dir}")
        return

    files_per_mean = TIMESTEPS_PER_DAY * n_days
    n_means = len(all_files) // files_per_mean

    if n_means == 0:
        print(f"Not enough files (<{files_per_mean})")
        return

    print(f"\n{'=' * 80}\n{subdir_name}: {len(all_files)} files → {n_means} mean(s)\n{'=' * 80}")

    total_time = 0

    for mean_idx in range(n_means):
        start = mean_idx * files_per_mean
        end = start + files_per_mean
        files_batch = all_files[start:end]

        first_file = files_batch[0].name
        last_file = files_batch[-1].name

        # Extraire infos du premier fichier: hl-b274_000000p_ls000_0000.nc
        # Le Ls complet est dans le nom: ls007_1234
        match_first = re.search(r'(hl-b274)_(\d+p)_(ls\d+)_(\d+)\.nc', first_file)

        if match_first:
            prefix = match_first.group(1)  # hl-b274
            file_idx = match_first.group(2)  # 000000p
            ls_part1 = match_first.group(3)  # ls007
            ls_part2 = match_first.group(4)  # 1234

            # Reconstituer Ls complet avec underscore: ls007_1234
            ls_full = f"{ls_part1}_{ls_part2}"

            # Calculer sol start et end depuis l'index fichier
            first_idx = extract_file_index(first_file)
            last_idx = extract_file_index(last_file)

            sol_start = first_idx // TIMESTEPS_PER_DAY
            sol_end = last_idx // TIMESTEPS_PER_DAY

            # Format: hl-b274_000000p_ls007_1234_sol000to004_5days_mean.nc
            output_name = f"{prefix}_{file_idx}_{ls_full}_sol{sol_start:03d}to{sol_end:03d}_{n_days}days_mean.nc"
        else:
            output_name = f"diurnal_mean_{n_days}days.nc"

        output_path = output_dir / output_name

        try:
            elapsed = compute_daily_mean(files_batch, output_path)
            total_time += elapsed
            print(f"Total: {elapsed:.1f}s\n")
        except Exception as e:
            print(f"ERROR: {e}")

    print(f"Completed {n_means} mean(s) in {total_time:.1f}s\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--netcdf', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--n-days', type=int, default=1)

    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--dir')
    selection.add_argument('--all', action='store_true')
    selection.add_argument('--dir-range', nargs=2, metavar=('START', 'END'))

    args = parser.parse_args()

    if not args.netcdf.exists():
        sys.exit(f"ERROR: {args.netcdf} not found")

    all_subdirs = sorted([d.name for d in args.netcdf.iterdir() if d.is_dir()])

    if args.dir:
        subdirs = [args.dir] if args.dir in all_subdirs else []
    elif args.dir_range:
        start, end = sorted(args.dir_range)
        subdirs = [d for d in all_subdirs if start <= d <= end]
    else:
        subdirs = all_subdirs

    if not subdirs:
        sys.exit("ERROR: No subdirectories")

    print(
        f"\n{'=' * 80}\nDIURNAL MEAN COMPUTATION\n{'=' * 80}\nInput: {args.netcdf}\nOutput: {args.output}\nDays: {args.n_days}\nDirs: {len(subdirs)}\n{'=' * 80}")

    for subdir in subdirs:
        process_subdirectory(args.netcdf, subdir, args.output, n_days=args.n_days)

    print(f"{'=' * 80}\nCOMPLETE\n{'=' * 80}\n")


if __name__ == "__main__":
    main()