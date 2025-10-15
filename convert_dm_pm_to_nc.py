#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert_dm_pm_to_nc.py — batch sur sous-dossiers dm/pm, sortie miroir dans netcdf/
"""

import argparse
import sys
import re
import time
import warnings
from pathlib import Path

import xarray as xr
import fstd2nc

# Progress bar optionnelle
try:
    from tqdm import tqdm
    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False

# Réduire le bruit console
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Variables par défaut
VARS_DEFAULT = [
    "TT", "PX", "GZ", "WW", "H2O", "CO2", "O3", "CO", "T9",
    "P0",
    "UU", "VV",
    "MLOC", "MALO", "MCZ", "MH", "MTSF", "MCO2", "MSN",
    "DVM1", "DVM2", "DVM3", "RWIC",
]

# ----------------- utilitaires physiques / I/O -----------------
def open_dataset(dm_path: Path, pm_path: Path | None, vars_keep: list[str] | None) -> xr.Dataset:
    files = [str(dm_path)] + ([str(pm_path)] if pm_path else [])
    kwargs = {}
    if vars_keep:
        kwargs["vars"] = vars_keep
    buf = fstd2nc.Buffer(files, **kwargs)
    return buf.to_xarray(fused=True)

def add_time_coordinate(ds: xr.Dataset) -> xr.Dataset:
    if ("reftime" in ds) and ("leadtime" in ds):
        ds = xr.decode_cf(ds, decode_timedelta=True)
        time_coord = ds["reftime"].astype("datetime64[ns]") + ds["leadtime"]
        ds = ds.assign_coords(time=time_coord)
        ds["time"].attrs.setdefault("standard_name", "time")
        ds["time"].attrs.setdefault("long_name", "time")
        ds.attrs.setdefault("Conventions", "CF-1.10")
    return ds

def convert_units(ds: xr.Dataset) -> xr.Dataset:
    ds = ds.copy()
    if "TT" in ds: ds["TT"] = ds["TT"] + 273.15
    if "PX" in ds: ds["PX"] = ds["PX"] * 100.0
    if "P0" in ds: ds["P0"] = ds["P0"] * 100.0
    if "GZ" in ds: ds["GZ"] = ds["GZ"] * 10.0
    if "UU" in ds: ds["UU"] = ds["UU"] * 0.5144
    if "VV" in ds: ds["VV"] = ds["VV"] * 0.5144
    return ds

def add_unit_attributes(ds: xr.Dataset) -> xr.Dataset:
    units = {
        "TT": "K", "PX": "Pa", "GZ": "m", "WW": "Pa s-1",
        "H2O": "1", "CO2": "1", "O3": "1", "CO": "1", "T9": "1",
        "UU": "m s-1", "VV": "m s-1",
        "P0": "Pa",
        "MLOC": "hour", "MALO": "1", "MCZ": "1", "MH": "m", "MTSF": "K",
        "MCO2": "kg m-2", "MSN": "micron",
        "DVM1": "1", "DVM2": "1", "DVM3": "1", "RWIC": "micron",
    }
    for v, u in units.items():
        if v in ds.variables:
            ds[v].attrs.setdefault("long_name", v)
            ds[v].attrs["units"] = u
    return ds

def find_vertical_dim(ds: xr.Dataset, var_name: str, wanted_size: int) -> str | None:
    if var_name not in ds:
        return None
    for dim in ds[var_name].dims:
        if ds.sizes.get(dim, -1) == wanted_size:
            return dim
    return None

def reduce_px_gz_to_thermo(ds: xr.Dataset) -> xr.Dataset:
    if ("PX" not in ds) or ("GZ" not in ds):
        return ds
    dim_px_204 = find_vertical_dim(ds, "PX", 204)
    dim_gz_204 = find_vertical_dim(ds, "GZ", 204)
    if (dim_px_204 is None) or (dim_gz_204 is None):
        return ds  # déjà à 103 ou structure non standard
    nlev = ds.sizes[dim_px_204]     # 204
    last = nlev - 1                 # 203 (surface)
    thermo_idx = list(range(0, nlev - 1, 2)) + [last]
    ds["PX"] = ds["PX"].rename({dim_px_204: "lev"}).isel(lev=thermo_idx)
    ds["GZ"] = ds["GZ"].rename({dim_gz_204: "lev"}).isel(lev=thermo_idx)
    return ds

def make_gz_height_above_surface(ds: xr.Dataset) -> xr.Dataset:
    if "GZ" not in ds:
        return ds
    dim = find_vertical_dim(ds, "GZ", 103) or find_vertical_dim(ds, "GZ", 204)
    if dim is None:
        return ds
    gz_surface = ds["GZ"].isel({dim: -1})
    ds["GZ"] = ds["GZ"] - gz_surface
    return ds

def default_compression_encoding(ds: xr.Dataset) -> dict:
    enc: dict = {}
    for v in ds.data_vars:
        enc[v] = {"zlib": True, "complevel": 4, "shuffle": True}
    return enc

def convert_one(dm_path: Path, pm_path: Path | None, out_path: Path, vars_keep: list[str] | None) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds = open_dataset(dm_path, pm_path, vars_keep=vars_keep)
    ds = add_time_coordinate(ds)
    ds = convert_units(ds)
    ds = reduce_px_gz_to_thermo(ds)
    ds = make_gz_height_above_surface(ds)
    ds = add_unit_attributes(ds)
    encoding = default_compression_encoding(ds)
    ds.to_netcdf(str(out_path), format="NETCDF4", encoding=encoding)

# ----------------- utilitaires batch sous-dossiers -----------------
_DM_INDEX_RE = re.compile(r"_dm_(\d+)p_")

def extract_dm_index(dm_name: str) -> int | None:
    m = _DM_INDEX_RE.search(dm_name)
    return int(m.group(1)) if m else None

def pm_for_dm(dm_file: Path, pm_subdir: Path) -> Path | None:
    cand = pm_subdir / dm_file.name.replace("_dm_", "_pm_")
    return cand if cand.exists() else None

def build_out_name(dm_file: Path) -> str:
    """
    'hl-b274_dm_000000p_ls000.0000' -> 'hl-b274_000000p_ls000_0000.nc'
    """
    extnum = dm_file.suffix.lstrip(".")
    stem = dm_file.stem.replace("_dm_", "_")
    return f"{stem}_{extnum}.nc"

def list_common_subdirs(dm_dir: Path, pm_dir: Path) -> list[Path]:
    dm_subs = {p.name: p for p in dm_dir.iterdir() if p.is_dir()}
    pm_subs = {p.name: p for p in pm_dir.iterdir() if p.is_dir()}
    names = sorted(set(dm_subs).intersection(pm_subs))
    return [dm_subs[n] for n in names]  # on renverra pour dm; on prendra pm via pm_dir/n

# ----------------- programme principal -----------------
def main():
    parser = argparse.ArgumentParser(
        description="FSTD -> NetCDF: miroir dm/pm par sous-dossier, PX/GZ=thermo+surface, GZ(surface)=0, compression lossless."
    )
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--dm", help="Chemin du fichier dm (.fst) - mode fichier unique")
    modes.add_argument("--root", help="Dossier racine contenant dm/, pm/, netcdf/ - mode batch")

    # Unique file
    parser.add_argument("--pm", help="Chemin du fichier pm (.fst) (mode fichier unique)")
    parser.add_argument("--out", help="Chemin du NetCDF de sortie (mode fichier unique)")

    # Sélection des sous-dossiers en mode batch
    subs = parser.add_mutually_exclusive_group()
    subs.add_argument("--all", action="store_true", help="Traiter tous les sous-dossiers (par défaut)")
    subs.add_argument("--dir", help="Un seul sous-dossier (ex: 000960)")
    subs.add_argument("--dir-range", nargs=2, metavar=("START", "END"),
                      help="Plage de sous-dossiers (lexico) ex: 000960 003840")

    # Sélection des fichiers dans chaque sous-dossier (par index dm)
    filesel = parser.add_mutually_exclusive_group()
    filesel.add_argument("--one", type=int, help="Un seul index (ex: 7 -> *_dm_000007p_*)")
    filesel.add_argument("--range", nargs=2, type=int, metavar=("DEBUT", "FIN"),
                         help="Plage d'index (ex: 0 3)")
    # Variables
    parser.add_argument("--vars", help="Variables a garder (liste separee par des virgules)")
    parser.add_argument("--all-vars", action="store_true", help="Tout garder (ignore --vars).")

    args = parser.parse_args()

    # Variables à charger
    if args.all_vars:
        vars_keep: list[str] | None = None
    elif args.vars:
        vars_keep = [v.strip() for v in args.vars.split(",") if v.strip()]
    else:
        vars_keep = list(VARS_DEFAULT)

    # --------- Mode fichier unique ---------
    if args.dm:
        dm = Path(args.dm)
        if not dm.exists(): print(f"[ERREUR] dm introuvable: {dm}", file=sys.stderr); sys.exit(2)
        pm = Path(args.pm) if args.pm else None
        if pm and not pm.exists(): print(f"[ERREUR] pm introuvable: {pm}", file=sys.stderr); sys.exit(2)
        if not args.out: print("[ERREUR] --out requis en mode fichier unique.", file=sys.stderr); sys.exit(2)

        t0 = time.perf_counter()
        print(f"[1/1] Converting {dm.name} ...")
        try:
            convert_one(dm, pm, Path(args.out), vars_keep)
            dt = time.perf_counter() - t0
            print(f"[WROTE] {args.out}  ({dt:.1f}s)")
        except Exception as e:
            print(f"[FAIL] {dm.name} : {e}", file=sys.stderr); sys.exit(1)
        return

    # --------- Mode batch (sous-dossiers) ---------
    root = Path(args.root)
    dm_dir, pm_dir, out_root = root / "dm", root / "pm", root / "netcdf"
    if not dm_dir.exists(): print(f"[ERREUR] dm introuvable: {dm_dir}", file=sys.stderr); sys.exit(2)
    if not pm_dir.exists(): print(f"[ERREUR] pm introuvable: {pm_dir}", file=sys.stderr); sys.exit(2)
    out_root.mkdir(parents=True, exist_ok=True)

    # Sous-dossiers communs dm/pm
    dm_subdirs = list_common_subdirs(dm_dir, pm_dir)
    if not dm_subdirs:
        print("[ERREUR] Aucun sous-dossier commun entre dm/ et pm/", file=sys.stderr); sys.exit(2)

    # Filtre sur les sous-dossiers
    if args.dir:
        dm_subdirs = [p for p in dm_subdirs if p.name == args.dir]
    elif args.dir_range:
        a, b = args.dir_range
        lo, hi = (a, b) if a <= b else (b, a)
        dm_subdirs = [p for p in dm_subdirs if lo <= p.name <= hi]
    # sinon --all (ou rien) -> tout

    if not dm_subdirs:
        print("[INFO] Aucun sous-dossier ne correspond à la sélection."); return

    # Pré-collecte des tâches pour barre globale
    tasks: list[tuple[Path, Path, Path]] = []
    for dm_sub in dm_subdirs:
        pm_sub = pm_dir / dm_sub.name
        out_sub = out_root / dm_sub.name
        # Tous les fichiers dm de ce sous-dossier
        dm_files = sorted(dm_sub.glob("*_dm_*.*"))
        if args.one is not None:
            dm_files = [f for f in dm_files if extract_dm_index(f.name) == args.one]
        elif args.range is not None:
            a, b = args.range; lo, hi = (a, b) if a <= b else (b, a)
            dm_files = [f for f in dm_files if (idx := extract_dm_index(f.name)) is not None and lo <= idx <= hi]
        # sinon tous
        for dm_file in dm_files:
            pm_file = pm_for_dm(dm_file, pm_sub)
            if pm_file is None:
                # on notera le SKIP pendant l'exécution
                tasks.append((dm_file, None, out_sub / build_out_name(dm_file)))
            else:
                tasks.append((dm_file, pm_file, out_sub / build_out_name(dm_file)))

    if not tasks:
        print("[INFO] Aucun fichier ne correspond à la sélection."); return

    total = len(tasks)
    t0 = time.perf_counter()
    iterator = tqdm(tasks, desc="Converting", unit="file") if _HAS_TQDM else tasks

    for i, (dm_file, pm_file, out_path) in enumerate(iterator, 1):
        # Log début (sans casser tqdm)
        if not _HAS_TQDM:
            elapsed = time.perf_counter() - t0
            eta = ""
            if i > 1 and elapsed > 0:
                rem = elapsed * (total - i) / i
                m, s = divmod(int(rem + 0.5), 60)
                eta = f"  ETA {m:02d}:{s:02d}"
            print(f"[{i}/{total}] {dm_file.parent.name}\\{dm_file.name} -> {out_path.name}{eta}")

        # PM manquant ?
        if pm_file is None or not pm_file.exists():
            msg = f"[SKIP] PM introuvable pour {dm_file.parent.name}\\{dm_file.name}"
            (tqdm.write(msg) if _HAS_TQDM else print(msg))
            continue

        try:
            convert_one(dm_file, pm_file, out_path, vars_keep)
            msg = f"[{i}/{total}] WROTE {out_path}"
            (tqdm.write(msg) if _HAS_TQDM else print("    " + msg))
        except Exception as e:
            msg = f"[FAIL] {dm_file.parent.name}\\{dm_file.name} : {e}"
            (tqdm.write(msg) if _HAS_TQDM else print(msg, file=sys.stderr))

    if _HAS_TQDM:
        dt = time.perf_counter() - t0
        tqdm.write(f"[DONE] {total} fichier(s) traité(s) en {dt:.1f}s")


if __name__ == "__main__":
    main()

