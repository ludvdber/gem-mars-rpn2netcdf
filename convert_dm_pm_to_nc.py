#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert_dm_pm_to_nc.py — Conversion FSTD vers NetCDF avec interpolation verticale
Batch processing sur sous-dossiers dm/pm, sortie miroir dans netcdf/
"""

import argparse
import sys
import re
import time
import warnings
from pathlib import Path

import xarray as xr
import fstd2nc
from netCDF4 import Dataset
import numpy as np

# Progress bar optionnelle
try:
    from tqdm import tqdm

    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False

# Réduire le bruit console
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ==================== CONFIGURATION ====================

# Variables par défaut à extraire
VARS_DEFAULT = [
    "TT", "PX", "GZ", "WW", "H2O", "CO2", "O3", "CO", "T9",
    "P0",
    "UU", "VV",
    "MLOC", "MALO", "MCZ", "MH", "MTSF", "MCO2", "MSN",
    "DVM1", "DVM2", "DVM3", "RWIC",
]

# Dimensions du modèle GEM-Mars
NLEVT = 103  # niveaux thermodynamiques
NLEVM = 102  # niveaux momentum
NLAT = 45  # latitudes
NLON = 91  # longitudes

# Configuration des variables à exporter
# Format: nom_var: (nom_nc, dimensions, unités, standard_name, long_name)
VARIABLE_CONFIG = {
    # === Variables 3D Thermodynamiques (103 niveaux) ===
    'TT': ('TT', ('time', 'altitudeT', 'lat', 'lon'), 'K',
           'air_temperature', 'Air temperature'),
    'PX': ('PX', ('time', 'altitudeT', 'lat', 'lon'), 'Pa',
           'air_pressure', 'Pressure'),
    'GZ': ('GZ', ('time', 'altitudeT', 'lat', 'lon'), 'm',
           'geopotential_height', 'Geopotential height above surface'),
    'WW': ('WW', ('time', 'altitudeT', 'lat', 'lon'), 'Pa s-1',
           'lagrangian_tendency_of_air_pressure', 'Vertical wind (omega)'),
    'H2O': ('H2O', ('time', 'altitudeT', 'lat', 'lon'), '1',
            'mass_fraction_of_water_in_air', 'Water vapour mixing ratio'),
    'CO2': ('CO2', ('time', 'altitudeT', 'lat', 'lon'), '1',
            'mass_fraction_of_carbon_dioxide_in_air', 'Carbon dioxide mixing ratio'),
    'O3': ('O3', ('time', 'altitudeT', 'lat', 'lon'), '1',
           'mass_fraction_of_ozone_in_air', 'Ozone mixing ratio'),
    'CO': ('CO', ('time', 'altitudeT', 'lat', 'lon'), '1',
           'mass_fraction_of_carbon_monoxide_in_air', 'Carbon monoxide mixing ratio'),
    'T9': ('T9', ('time', 'altitudeT', 'lat', 'lon'), '1',
           'mass_fraction_of_water_ice_in_air', 'Water ice mixing ratio (clouds)'),
    'DVM1': ('DVM1', ('time', 'altitudeT', 'lat', 'lon'), '1',
             'mass_fraction_of_dust_bin_1_in_air', 'Dust mixing ratio (0.1 µm)'),
    'DVM2': ('DVM2', ('time', 'altitudeT', 'lat', 'lon'), '1',
             'mass_fraction_of_dust_bin_2_in_air', 'Dust mixing ratio (1.5 µm)'),
    'DVM3': ('DVM3', ('time', 'altitudeT', 'lat', 'lon'), '1',
             'mass_fraction_of_dust_bin_3_in_air', 'Dust mixing ratio (10 µm)'),
    'RWIC': ('RWIC', ('time', 'altitudeT', 'lat', 'lon'), 'micron',
             'effective_radius_of_water_ice_particles', 'Water ice particle effective radius'),

    # === Variables 3D Momentum (102 niveaux) - TODO ===
    # 'UU': ('UU', ('time', 'altitudeM', 'lat', 'lon'), 'm s-1',
    #        'eastward_wind', 'Eastward wind (E-W)'),
    # 'VV': ('VV', ('time', 'altitudeM', 'lat', 'lon'), 'm s-1',
    #        'northward_wind', 'Northward wind (N-S)'),

    # === Variables 2D Surface ===
    'P0': ('P0', ('time', 'lat', 'lon'), 'Pa',
           'surface_air_pressure', 'Surface pressure'),
    'MLOC': ('MLOC', ('time', 'lat', 'lon'), 'hour',
             'local_time', 'Local time'),
    'MALO': ('MALO', ('time', 'lat', 'lon'), '1',
             'surface_albedo', 'Surface albedo'),
    'MCZ': ('MCZ', ('time', 'lat', 'lon'), '1',
            'solar_zenith_angle', 'Cosine of solar zenith angle'),
    'MH': ('MH', ('time', 'lat', 'lon'), 'm',
           'atmosphere_boundary_layer_thickness', 'Planetary boundary layer height'),
    'MTSF': ('MTSF', ('time', 'lat', 'lon'), 'K',
             'surface_temperature', 'Surface temperature'),
    'MCO2': ('MCO2', ('time', 'lat', 'lon'), 'kg m-2',
             'atmosphere_mass_content_of_carbon_dioxide', 'CO2 column'),
    'MSN': ('MSN', ('time', 'lat', 'lon'), 'micron',
            'surface_water_ice_amount', 'Water ice on surface (precipitable microns)'),
}


# ==================== FONCTIONS UTILITAIRES ====================

def open_dataset(dm_path: Path, pm_path: Path | None, vars_keep: list[str] | None) -> xr.Dataset:
    """Charge les fichiers FSTD dm/pm et retourne un xarray Dataset fusionné."""
    files = [str(dm_path)]
    if pm_path:
        files.append(str(pm_path))

    kwargs = {"vars": vars_keep} if vars_keep else {}
    buf = fstd2nc.Buffer(files, **kwargs)
    return buf.to_xarray(fused=True)


def add_time_coordinate(ds: xr.Dataset) -> xr.Dataset:
    """Ajoute une coordonnée temporelle CF-compliant."""
    if "reftime" in ds and "leadtime" in ds:
        ds = xr.decode_cf(ds, decode_timedelta=True)
        time_coord = ds["reftime"].astype("datetime64[ns]") + ds["leadtime"]
        ds = ds.assign_coords(time=time_coord)
        ds["time"].attrs.setdefault("standard_name", "time")
        ds["time"].attrs.setdefault("long_name", "time")
        ds.attrs.setdefault("Conventions", "CF-1.10")
    return ds


def convert_units(ds: xr.Dataset) -> xr.Dataset:
    """Convertit les unités vers le système SI et conventions CF."""
    ds = ds.copy()

    # Températures: C → K
    if "TT" in ds:
        ds["TT"] = ds["TT"] + 273.15
    if "MTSF" in ds:
        ds["MTSF"] = ds["MTSF"] + 273.15

    # Pressions: hPa → Pa
    if "PX" in ds:
        ds["PX"] = ds["PX"] * 100.0
    if "P0" in ds:
        ds["P0"] = ds["P0"] * 100.0

    # Géopotentiel: dam → m
    if "GZ" in ds:
        ds["GZ"] = ds["GZ"] * 10.0

    # Vents: knots → m/s
    if "UU" in ds:
        ds["UU"] = ds["UU"] * 0.5144
    if "VV" in ds:
        ds["VV"] = ds["VV"] * 0.5144

    return ds


def reduce_px_gz_to_thermo(ds: xr.Dataset) -> xr.Dataset:
    """
    Réduit PX et GZ de 204 niveaux (momentum + thermo) à 103 niveaux (thermo uniquement).

    La grille level2 (204) alterne momentum et thermo: M T M T M T ... M T S
    On extrait les indices pairs (thermo) + le dernier (surface).
    """
    if "PX" not in ds or "GZ" not in ds:
        return ds

    # Vérifier si les variables sont sur level2 (204 niveaux)
    if "level2" not in ds["PX"].dims or "level2" not in ds["GZ"].dims:
        return ds

    # Indices thermodynamiques: 0, 2, 4, ..., 202, 203
    thermo_idx = list(range(0, 203, 2)) + [203]

    ds["PX"] = ds["PX"].rename({"level2": "level1"}).isel(level1=thermo_idx)
    ds["GZ"] = ds["GZ"].rename({"level2": "level1"}).isel(level1=thermo_idx)

    return ds


def make_gz_height_above_surface(ds: xr.Dataset) -> xr.Dataset:
    """Transforme GZ en hauteur relative à la surface (GZ_surface = 0)."""
    if "GZ" not in ds or "level1" not in ds["GZ"].dims:
        return ds

    gz_surface = ds["GZ"].isel(level1=-1)
    ds["GZ"] = ds["GZ"] - gz_surface
    return ds


def interp_along_column(data_col, z_col, z_new):
    """Interpole une colonne de données sur une nouvelle grille verticale."""
    # Trier par altitude croissante
    order = np.argsort(z_col)
    z_col_sorted = z_col[order]
    data_col_sorted = data_col[order]

    # Interpolation linéaire
    return np.interp(z_new, z_col_sorted, data_col_sorted,
                     left=data_col_sorted[0], right=data_col_sorted[-1])


def interpolate_to_common_grid(ds: xr.Dataset, var_name: str, gz_mean: xr.DataArray) -> xr.DataArray | None:
    """
    Interpole une variable 3D sur la grille verticale commune (GZ moyen lat/lon).
    Méthode validée avec Paraview.
    """
    if var_name not in ds or 'level1' not in ds[var_name].dims:
        return None

    var_data = ds.variables[var_name][:]
    var_gz = ds.variables['GZ'][:]

    # Interpolation vectorisée sur toutes les colonnes
    var_interp = xr.apply_ufunc(
        interp_along_column,
        var_data,
        var_gz,
        gz_mean,
        input_core_dims=[['level1'], ['level1'], ['level1']],
        output_core_dims=[['level1']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[var_data.dtype],
    )

    # Réorganise: (time, level1, lat, lon)
    return var_interp.transpose("time", "level1", "lat", "lon")


def write_netcdf(ds: xr.Dataset, out_path: Path, dm_path: Path):
    """
    Écrit le fichier NetCDF avec toutes les variables configurées.
    Format: NetCDF4 avec compression zlib.
    """
    # Validation
    if 'GZ' not in ds.variables:
        raise ValueError("GZ variable required for altitude grid")

    # Calcul de la grille verticale moyenne
    gz_mean = ds.variables['GZ'][:].mean(dim=['lat', 'lon'])

    if gz_mean.shape[1] != NLEVT:
        raise ValueError(f"Expected {NLEVT} vertical levels, got {gz_mean.shape[1]}")

    # Création du fichier
    with Dataset(out_path, mode='w', format='NETCDF4') as ncfile:

        # === DIMENSIONS ===
        ncfile.createDimension('lat', NLAT)
        ncfile.createDimension('lon', NLON)
        ncfile.createDimension('time', None)  # unlimited
        ncfile.createDimension('altitudeT', NLEVT)

        # === ATTRIBUTS GLOBAUX ===
        ncfile.title = 'GEM-Mars simulation output'
        ncfile.subtitle = f'Converted from: {dm_path.name}'
        ncfile.institution = "BIRA-IASB"
        ncfile.Conventions = "CF-1.10"
        ncfile.source = "GEM-Mars atmospheric model"
        ncfile.history = f"Created {time.strftime('%Y-%m-%d %H:%M:%S UTC')}"

        # === COORDONNÉES ===

        # Latitude
        lat_var = ncfile.createVariable('lat', np.float32, ('lat',))
        lat_var.units = 'degrees_north'
        lat_var.long_name = 'latitude'
        lat_var.standard_name = 'latitude'
        lat_var[:] = ds.coords['lat'].values if 'lat' in ds.coords else \
            -88. + (180. / NLAT) * np.arange(NLAT)

        # Longitude
        lon_var = ncfile.createVariable('lon', np.float32, ('lon',))
        lon_var.units = 'degrees_east'
        lon_var.long_name = 'longitude'
        lon_var.standard_name = 'longitude'
        lon_values = ds.coords['lon'].values if 'lon' in ds.coords else \
            (360. / NLON) * np.arange(NLON)
        lon_var[:] = lon_values
        lon_var.comment = f'Model grid: {lon_values[0]:.1f} to {lon_values[-1]:.1f} degrees'

        # Temps
        time_var = ncfile.createVariable('time', np.float64, ('time',))
        time_var.units = 'hours since 0001-01-01 00:00:00'
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        time_var.calendar = 'proleptic_gregorian'

        time_data = ds.variables['time'][:]
        time_ref = np.datetime64('0001-01-01T00:00:00')
        time_hours = (time_data - time_ref) / np.timedelta64(1, 'h')
        time_var[:] = time_hours

        # Altitude (thermodynamique)
        alt_var = ncfile.createVariable('altitudeT', np.float64, ('altitudeT',))
        alt_var.units = 'km'
        alt_var.long_name = 'Altitude on thermodynamic levels'
        alt_var.positive = 'up'
        alt_var.comment = 'Height above surface, averaged over lat/lon'
        alt_var[:] = gz_mean.values / 1000.0  # m → km

        # === VARIABLES DE DONNÉES ===

        variables_skipped = []

        for var_name, (nc_name, dims, units, std_name, long_name) in VARIABLE_CONFIG.items():
            if var_name not in ds:
                continue

            try:
                is_3d = len(dims) == 4

                if is_3d:
                    # Variable 3D: interpolation verticale
                    data_interp = interpolate_to_common_grid(ds, var_name, gz_mean)
                    if data_interp is None:
                        variables_skipped.append(var_name)
                        continue
                    data = data_interp.values
                else:
                    # Variable 2D: pas d'interpolation
                    data = ds.variables[var_name][:].values

                # Création variable avec compression
                var = ncfile.createVariable(nc_name, np.float64, dims,
                                            zlib=True, complevel=4, shuffle=True)
                var[:] = data

                # Attributs
                var.units = units
                var.standard_name = std_name
                var.long_name = long_name

            except Exception as e:
                variables_skipped.append(f"{var_name} ({e})")

        # Message de sortie
        if variables_skipped:
            print(f" Skipped: {', '.join(variables_skipped)}")


def convert_one(dm_path: Path, pm_path: Path | None, out_path: Path, vars_keep: list[str] | None):
    """Pipeline complet de conversion d'un fichier FSTD vers NetCDF."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 1. Chargement
    ds = open_dataset(dm_path, pm_path, vars_keep)

    # 2. Prétraitement
    ds = add_time_coordinate(ds)
    ds = convert_units(ds)
    ds = reduce_px_gz_to_thermo(ds)
    ds = make_gz_height_above_surface(ds)

    # 3. Écriture NetCDF
    write_netcdf(ds, out_path, dm_path)


# ==================== UTILITAIRES BATCH ====================

_DM_INDEX_RE = re.compile(r"_dm_(\d+)p_")


def extract_dm_index(dm_name: str) -> int | None:
    """Extrait l'index numérique d'un nom de fichier dm."""
    match = _DM_INDEX_RE.search(dm_name)
    return int(match.group(1)) if match else None


def pm_for_dm(dm_file: Path, pm_subdir: Path) -> Path | None:
    """Trouve le fichier pm correspondant à un fichier dm."""
    pm_candidate = pm_subdir / dm_file.name.replace("_dm_", "_pm_")
    return pm_candidate if pm_candidate.exists() else None


def build_out_name(dm_file: Path) -> str:
    """
    Construit le nom du fichier NetCDF de sortie.
    Exemple: 'hl-b274_dm_000000p_ls000.0000' → 'hl-b274_000000p_ls000_0000.nc'
    """
    extnum = dm_file.suffix.lstrip(".")
    stem = dm_file.stem.replace("_dm_", "_")
    return f"{stem}_{extnum}.nc"


def list_common_subdirs(dm_dir: Path, pm_dir: Path) -> list[Path]:
    """Liste les sous-dossiers communs entre dm/ et pm/."""
    dm_subs = {p.name: p for p in dm_dir.iterdir() if p.is_dir()}
    pm_subs = {p.name: p for p in pm_dir.iterdir() if p.is_dir()}
    common_names = sorted(set(dm_subs) & set(pm_subs))
    return [dm_subs[name] for name in common_names]


# ==================== PROGRAMME PRINCIPAL ====================

def main():
    """Point d'entrée principal avec gestion des arguments."""
    parser = argparse.ArgumentParser(
        description="Conversion FSTD → NetCDF avec interpolation verticale",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples:
  # Fichier unique
  %(prog)s --dm file.fst --pm file_pm.fst --out output.nc

  # Batch - un sous-dossier
  %(prog)s --root /path/to/hl-b274 --dir 000960

  # Batch - plage de sous-dossiers
  %(prog)s --root /path/to/hl-b274 --dir-range 000960 003840

  # Batch - tous les sous-dossiers
  %(prog)s --root /path/to/hl-b274 --all
        """
    )

    # Mode de fonctionnement
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--dm", help="Fichier dm (.fst) - mode fichier unique")
    modes.add_argument("--root", help="Racine contenant dm/, pm/, netcdf/ - mode batch")

    # Arguments mode fichier unique
    parser.add_argument("--pm", help="Fichier pm (.fst)")
    parser.add_argument("--out", help="Fichier NetCDF de sortie")

    # Sélection sous-dossiers (mode batch)
    subs = parser.add_mutually_exclusive_group()
    subs.add_argument("--all", action="store_true", help="Tous les sous-dossiers")
    subs.add_argument("--dir", help="Un seul sous-dossier (ex: 000960)")
    subs.add_argument("--dir-range", nargs=2, metavar=("START", "END"),
                      help="Plage de sous-dossiers (ex: 000960 003840)")

    # Sélection fichiers (mode batch)
    filesel = parser.add_mutually_exclusive_group()
    filesel.add_argument("--one", type=int, help="Un seul index (ex: 0)")
    filesel.add_argument("--range", nargs=2, type=int, metavar=("START", "END"),
                         help="Plage d'index (ex: 0 5)")

    # Variables
    parser.add_argument("--vars", help="Variables à garder (séparées par virgules)")
    parser.add_argument("--all-vars", action="store_true", help="Toutes les variables")

    args = parser.parse_args()

    # Détermination des variables à charger
    if args.all_vars:
        vars_keep = None
    elif args.vars:
        vars_keep = [v.strip() for v in args.vars.split(",") if v.strip()]
    else:
        vars_keep = list(VARS_DEFAULT)

    # ========== MODE FICHIER UNIQUE ==========
    if args.dm:
        dm = Path(args.dm)
        if not dm.exists():
            sys.exit(f"[ERROR] DM file not found: {dm}")

        pm = Path(args.pm) if args.pm else None
        if pm and not pm.exists():
            sys.exit(f"[ERROR] PM file not found: {pm}")

        if not args.out:
            sys.exit("[ERROR] --out required in single file mode")

        t0 = time.perf_counter()
        print(f"[1/1] Converting {dm.name}...")

        try:
            convert_one(dm, pm, Path(args.out), vars_keep)
            dt = time.perf_counter() - t0
            print(f"✅ Success: {args.out} ({dt:.1f}s)")
        except Exception as e:
            print(f"❌ Failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)

        return

    # ========== MODE BATCH ==========
    root = Path(args.root)
    dm_dir = root / "dm"
    pm_dir = root / "pm"
    out_root = root / "netcdf"

    # Validation des dossiers
    if not dm_dir.exists():
        sys.exit(f"[ERROR] DM directory not found: {dm_dir}")
    if not pm_dir.exists():
        sys.exit(f"[ERROR] PM directory not found: {pm_dir}")

    out_root.mkdir(parents=True, exist_ok=True)

    # Sous-dossiers communs
    dm_subdirs = list_common_subdirs(dm_dir, pm_dir)
    if not dm_subdirs:
        sys.exit("[ERROR] No common subdirectories between dm/ and pm/")

    # Filtrage sous-dossiers
    if args.dir:
        dm_subdirs = [p for p in dm_subdirs if p.name == args.dir]
    elif args.dir_range:
        start, end = sorted(args.dir_range)
        dm_subdirs = [p for p in dm_subdirs if start <= p.name <= end]

    if not dm_subdirs:
        print("[INFO] No subdirectories match selection")
        return

    # Collection des tâches
    tasks = []
    for dm_sub in dm_subdirs:
        pm_sub = pm_dir / dm_sub.name
        out_sub = out_root / dm_sub.name

        dm_files = sorted(dm_sub.glob("*_dm_*.*"))

        # Filtrage par index
        if args.one is not None:
            dm_files = [f for f in dm_files if extract_dm_index(f.name) == args.one]
        elif args.range is not None:
            start, end = sorted(args.range)
            dm_files = [f for f in dm_files
                        if (idx := extract_dm_index(f.name)) is not None and start <= idx <= end]

        for dm_file in dm_files:
            pm_file = pm_for_dm(dm_file, pm_sub)
            out_path = out_sub / build_out_name(dm_file)
            tasks.append((dm_file, pm_file, out_path))

    if not tasks:
        print("[INFO] No files match selection")
        return

    # Traitement batch
    total = len(tasks)
    t0 = time.perf_counter()
    iterator = tqdm(tasks, desc="Converting", unit="file") if _HAS_TQDM else tasks

    success_count = 0
    fail_count = 0

    for i, (dm_file, pm_file, out_path) in enumerate(iterator, 1):
        # Affichage sans tqdm
        if not _HAS_TQDM:
            elapsed = time.perf_counter() - t0
            rate = elapsed / i
            eta_sec = rate * (total - i)
            eta_min, eta_sec_rem = divmod(int(eta_sec), 60)

            print(f"\n{'─' * 80}")
            print(f"📄 File {i}/{total}: {dm_file.name}")
            print(f"📁 Folder: {dm_file.parent.name}")
            print(f"⏱️  {i / total * 100:.1f}% • {rate:.1f}s/file • ETA {eta_min:02d}:{eta_sec_rem:02d}")
            print(f"{'─' * 80}")

        # Vérification PM
        if not pm_file or not pm_file.exists():
            msg = "❌ [SKIP] PM file not found"
            (tqdm.write(msg) if _HAS_TQDM else print(msg))
            fail_count += 1
            continue

        # Conversion
        try:
            t_start = time.perf_counter()
            convert_one(dm_file, pm_file, out_path, vars_keep)
            t_elapsed = time.perf_counter() - t_start

            success_count += 1
            msg = f"✅ Wrote: {out_path.name} ({t_elapsed:.1f}s)"
            (tqdm.write(msg) if _HAS_TQDM else print(msg))

        except Exception as e:
            fail_count += 1
            msg = f"❌ [FAIL] {e}"
            (tqdm.write(msg) if _HAS_TQDM else print(msg, file=sys.stderr))

    # Résumé final
    dt = time.perf_counter() - t0
    avg_time = dt / total if total > 0 else 0

    print(f"\n{'═' * 80}")
    print(f"{' CONVERSION COMPLETE':^80}")
    print(f"{'═' * 80}")
    print(f" Success:     {success_count}/{total} files")
    if fail_count > 0:
        print(f" ❌ Failed:      {fail_count}/{total} files")
    print(f" Total time:  {dt:.1f}s ({avg_time:.1f}s per file)")
    print(f"{'═' * 80}\n")


if __name__ == "__main__":
    main()
