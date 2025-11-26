![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![CF Conventions](https://img.shields.io/badge/CF-1.10-orange.svg)
![Status](https://img.shields.io/badge/status-production-brightgreen.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Windows-lightgrey.svg)
![NetCDF](https://img.shields.io/badge/format-NetCDF--4-blue.svg)
![Mars](https://img.shields.io/badge/planet-Mars-red.svg)

> 🌍 **[English version](README.md)** | **Version française**

**Outils CLI pour le modèle atmosphérique GEM-Mars :** Conversion de fichiers RPN/FSTD vers NetCDF-4 conforme CF + calcul de cycles moyens diurnes avec indexation Mars Year/Ls.

## 📑 Table des matières

- [convert_dm_pm_to_nc.py](#convert_dm_pm_to_nc--readme)
  - [Fonctionnalités](#fonctionnalités-comportement-par-défaut)
  - [Structure des dossiers](#structure-des-dossiers-attendue)
  - [Installation](#installation)
  - [Utilisation](#utilisation)
  - [Sélection des variables](#sélection-des-variables)
  - [Sortie et journalisation](#sortie-et-journalisation)
  - [Dépannage](#dépannage)
- [compute_diurnal_mean.py](#calcul-des-cycles-moyens-diurnes)
  - [Fonctionnalités](#fonctionnalités)
  - [Installation](#installation-1)
  - [Utilisation de base](#utilisation-de-base)
  - [Options avancées](#options-avancées)
  - [Recherche Mars Year](#recherche-mars-year-avancé)
  - [Référence des commandes](#référence-des-commandes)
- [Ressources associées](#ressources-associées)

---

# convert_dm_pm_to_nc — README

Un outil CLI compact pour convertir les fichiers RPN GEM-Mars (dm + pm) vers NetCDF4, avec des paramètres par défaut adaptés à un usage scientifique et au traitement par lots de nombreux sous-dossiers.

## Fonctionnalités (comportement par défaut)

- **Traitement de fichiers appariés** : Lit les fichiers RPN dm + pm via fstd2nc
- **Coordonnée temporelle** : Calcule `time = reftime + leadtime`
- **Conversions d'unités** :
  - TT : °C → K (+273.15)
  - PX, P0 : hPa → Pa (×100)
  - GZ : décamètre → m (×10)
  - UU, VV : nœuds → m s⁻¹ (×0.5144)
  - WW : déjà en Pa s⁻¹ (inchangé)
- **Interpolation verticale** :
  - Variables thermodynamiques (TT, PX, GZ, etc.) interpolées sur 103 niveaux d'altitude communs
  - Variables de quantité de mouvement (UU, VV) interpolées sur 102 niveaux d'altitude communs
  - GZ défini comme hauteur au-dessus de la surface (GZ_surface = 0)
- **Sortie optimisée** :
  - Format NetCDF-4 avec compression sans perte
  - Précision float32 (réduit la taille des fichiers)
  - Métadonnées conformes CF-1.10
- **Suivi de progression** : Barre de progression tqdm optionnelle

## Structure des dossiers (attendue)

Votre répertoire racine du projet (par exemple, `hl-b274`) doit contenir :

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
└── netcdf/           # les fichiers de sortie seront créés ici
    ├── 000960/
    │   ├── hl-b274_000000p_ls000_0000.nc
    │   ├── hl-b274_000001p_ls000_0100.nc
    │   └── ...
    └── ...
```

Pour chaque `dm/SUBDIR/file`, le script cherche le `pm/SUBDIR/file` correspondant et écrit dans `netcdf/SUBDIR/`.

**Nommage de sortie** : Conserve les numéros de l'extension RPN :

```
hl-b274_dm_000000p_ls000.0000  →  netcdf/000960/hl-b274_000000p_ls000_0000.nc
```

## Installation

- Python 3.11+ recommandé.
  - xarray
  - netCDF4
  - fstd2nc
  - tqdm
  - numpy

tqdm est optionnel (barre de progression agréable).

```
pip install numpy xarray netCDF4 fstd2nc tqdm
```

Bibliothèques système Linux (si nécessaire)

Debian/Ubuntu :

```
sudo apt-get install libhdf5-dev libnetcdf-dev
```

## Utilisation

Définissez votre ROOT (chemin vers le dossier contenant dm/, pm/, netcdf/)

- ### Linux

    ```
    ROOT="$HOME/hl-b274"
    cd /path/to/your/repo  # dossier où se trouve convert_dm_pm_to_nc.py
    ```

  - ### Windows PowerShell

      ```
      $ROOT = "$env:USERPROFILE\Desktop\hl-b274"
      Set-Location C:\Users\USERPROFILE\PythonFile
     ```
  
  1. Tous les sous-dossiers et tous les fichiers

  - ### Linux

      ```
      python convert_dm_pm_to_nc.py --root "$ROOT" --all
     ```

  - ### Windows PowerShell

      ```
      python .\convert_dm_pm_to_nc.py --root "$ROOT" --all
     ```
  
1. Un sous-dossier spécifique (par exemple, 000960)

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960
```

1. Une plage de sous-dossiers (inclusive, ordre lexical)

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir-range 000960 003840
```

1. Dans les sous-dossiers sélectionnés : un seul index de fichier (par exemple, index 7. l'index commence à 0)

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --one 7
```

1. Dans les sous-dossiers sélectionnés : une plage d'indices de fichiers (inclusive)

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --range 0 3
```

Notes

"file index" fait référence au nombre ```*_dm_000007p_*``` → 7.

Si le fichier pm correspondant est manquant, le script enregistre SKIP et continue.

Avec tqdm installé, vous obtenez une barre de progression et une ligne ```[i/total] WROTE <path>``` par sortie.

## Sélection des variables

Par défaut, le script charge exactement les variables demandées dans les spécifications du projet :

### Variables par défaut (28 au total)

**3D Thermodynamique** (103 niveaux) :

- `TT` - Température (K)
- `PX` - Pression (Pa)
- `GZ` - Hauteur géopotentielle (m au-dessus de la surface)
- `WW` - Vitesse verticale (Pa/s)
- `H2O`, `CO2`, `O3`, `CO` - Gaz traces
- `T9` - Perturbation de température
- `DVM1`, `DVM2`, `DVM3` - Moments de poussière
- `RWIC` - Rayon des particules de glace

**3D Quantité de mouvement** (102 niveaux) :

- `UU` - Vent vers l'est (m/s)
- `VV` - Vent vers le nord (m/s)

**2D Surface** :

- `P0` - Pression de surface (Pa)
- `MTSF` - Température de surface (K)
- `MLOC` - Heure locale
- `MALO`, `MCZ`, `MH`, `MCO2`, `MSN` - Champs de surface supplémentaires

### Dimensions de sortie

- `time` : 1 (dimension illimitée)
- `lat` : 45 (latitudes, -88° à 88°)
- `lon` : 91 (longitudes, 0° à 360°)
- `altitudeT` : 103 (niveaux verticaux thermodynamiques)
- `altitudeM` : 102 (niveaux verticaux de quantité de mouvement)

Vous pouvez surcharger :

- Conserver toutes les variables présentes :

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --all-vars
```

- Conserver une liste personnalisée (séparée par des virgules, sans espaces) :

```
python convert_dm_pm_to_nc.py --root "$ROOT" --dir 000960 --vars "TT,PX,GZ,WW,P0,UU,VV"
```

## Sortie et journalisation

NetCDF écrit dans ```netcdf/<SUBDIR>/...nc``` avec zlib niveau 6 (sans perte) et shuffle.

Sortie console minimale :

Avec tqdm : ```barre de progression + [i/total] WROTE <output_path>``` par fichier.

Sans tqdm : lignes de progression ```[i/total] avec ETA et WROTE```.

## Dépannage

Homologue pm manquant → le fichier est ignoré (SKIP) ; vérifiez que ```pm/<SUBDIR>/``` contient les mêmes noms de fichiers que ```dm/<SUBDIR>/ (avec _pm_)```.

Permissions / disque → assurez-vous que netcdf/ est accessible en écriture et dispose de suffisamment d'espace.

Erreurs de compilation Linux pour netCDF4/HDF5 → installez les bibliothèques de développement système (voir Installation).

---

## Calcul des cycles moyens diurnes

Après avoir converti les fichiers RPN en NetCDF, vous pouvez calculer les cycles moyens diurnes (moyennés sur N sols martiens) en utilisant `compute_diurnal_mean.py`.

### Fonctionnalités

- Regroupe les fichiers NetCDF par heure du jour (48 pas de temps : 0.0h, 0.5h, ..., 23.5h)
- Moyenne chaque heure sur plusieurs jours martiens
- Génère un fichier avec 48 pas de temps représentant un cycle typique de 24 heures
- Utilise la compilation JIT optimisée numba pour un traitement rapide
- Préserve toutes les variables et métadonnées

### Installation

Nécessite `numba` en plus des dépendances précédentes :

```bash
pip install numba
```

### Utilisation de base

**1. Répertoire unique, N jours :**

```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5
```

**2. Tous les répertoires, 10 jours chacun :**

```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --n-days 10
```

**3. Plage de répertoires :**

```bash
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 003840 \
  --n-days 20
```

### Options avancées

#### Limiter le nombre de moyennes (`--max-means`)

Créer uniquement les N premières moyennes (utile pour les tests ou pour limiter la sortie) :

```bash
# Créer seulement 5 moyennes de 5 jours chacune (= 25 jours au total)
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5 \
  --max-means 5
```

**Sortie :** 5 fichiers (sols 0-4, 5-9, 10-14, 15-19, 20-24)

#### Filtrer par plage Ls (`--ls-range`)

Traiter uniquement les fichiers dans une plage de longitude solaire (Ls) spécifique :

```bash
# Traiter uniquement Ls 0 à 15°
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --n-days 10 \
  --ls-range 0 15
```

#### Mode inter-répertoires (`--cross-dirs`)

Permettre aux moyennes de s'étendre sur plusieurs répertoires (traite tous les répertoires sélectionnés comme un ensemble de données continu) :

```bash
# Moyennes de 30 jours pouvant s'étendre sur les limites des répertoires
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 001920 \
  --n-days 30 \
  --cross-dirs
```

**Sans `--cross-dirs` :** Chaque répertoire traité séparément  
**Avec `--cross-dirs` :** Tous les fichiers traités comme une chronologie continue

**Emplacement de sortie :** `netcdf_mean/cross_dirs/`

#### Mode moyenne unique (`--single-mean`)

Créer UNE seule moyenne à partir de TOUS les fichiers dans la plage spécifiée (ignore `--n-days` et `--max-means`) :

```bash
# Créer une moyenne Ls 0 à 90°
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --all \
  --ls-range 0 90 \
  --single-mean
```

**Sortie :** Un fichier  
**Emplacement de sortie :** `netcdf_mean/single_mean/`

### Exemples combinés

**Exemple 1 : Premier répertoire complet + 5 jours du deuxième**

```bash
# Traiter complètement le premier répertoire
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 000960 \
  --n-days 5

# Traiter uniquement la première moyenne (5 jours) du deuxième répertoire
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir 001920 \
  --n-days 5 \
  --max-means 1
```

**Exemple 2 : Inter-répertoires avec filtres**

```bash
# 5 moyennes de 5 jours, seulement Ls 0-30°, sur plusieurs répertoires
python compute_diurnal_mean.py \
  --netcdf "$ROOT/netcdf" \
  --output "$ROOT/netcdf_mean" \
  --dir-range 000960 003840 \
  --n-days 5 \
  --max-means 5 \
  --ls-range 0 30 \
  --cross-dirs
```

### Nommage de sortie

**Mode standard :**

```
Entrée :  netcdf/000960/hl-b274_000000p_ls000_0000.nc
Sortie : netcdf_mean/000960/hl-b274_000000p_ls000_0000_sol000to004_5days_mean.nc
```

**Mode inter-répertoires :**

```
Sortie : netcdf_mean/cross_dirs/hl-b274_000000p_ls000_0000_sol000to029_30days_mean_crossdir.nc
```

**Mode moyenne unique :**

```
Sortie : netcdf_mean/single_mean/hl-b274_000000p_ls000_0000_to_ls090_0000_sol000to500_501days_single_mean.nc
```

Le nom de fichier inclut :

- Index de fichier de départ (`000000p`)
- Plage Ls (`ls000_0000` ou `ls000_0000_to_ls090_0000`)
- Plage de sols (`sol000to004`)
- Nombre de jours moyennés (`5days`)
- Indicateur de mode (`mean`, `crossdir`, ou `single_mean`)

### Référence des commandes

```
usage: compute_diurnal_mean.py --netcdf PATH --output PATH 
                               [--dir DIR | --all | --dir-range START END]
                               [--n-days N] [--max-means N] 
                               [--ls-range MIN MAX] [--mars-year MY] [--ls-start LS]
                               [--lookup FILE] [--cross-dirs] [--single-mean]

Arguments requis :
  --netcdf PATH          Répertoire racine contenant les sous-répertoires NetCDF
  --output PATH          Répertoire de sortie pour les fichiers moyens
  
Sélection de répertoire (requis sauf si --mars-year est utilisé) :
  --dir DIR              Traiter un seul sous-répertoire
  --all                  Traiter tous les sous-répertoires
  --dir-range START END  Traiter une plage de sous-répertoires

Arguments optionnels :
  --n-days N            Nombre de jours à moyenner par moyenne (défaut : 1)
                        Ignoré avec --single-mean
  --max-means N         Nombre maximal de moyennes à calculer (défaut : toutes)
                        Ignoré avec --single-mean
  --ls-range MIN MAX    Filtrer les fichiers par plage Ls (ex. : --ls-range 0 15)
                        Exemple : 0 90 pour le printemps, 90 180 pour l'été
  --cross-dirs          Permettre aux moyennes de s'étendre sur les répertoires
                        Traite tous les répertoires sélectionnés comme un ensemble continu
  --single-mean         Créer UNE moyenne à partir de TOUS les fichiers de la plage
                        Utile pour les moyennes saisonnières/périodiques

Arguments Mars Year (avancé) :
  --mars-year MY        Numéro d'année martienne (ex. : 34, 35)
                        Active automatiquement la recherche inter-répertoires
                        Rend la sélection de répertoire optionnelle
  --ls-start LS         Ls de départ pour le mode Mars Year (ex. : 9.61)
  --lookup FILE         Fichier de recherche Mars Year (.xlsx ou .csv)
                        Requis lors de l'utilisation de --mars-year
```

### Recherche Mars Year (avancé)

Trouver automatiquement les répertoires de départ en utilisant l'année martienne et Ls !

**Utilisation simplifiée** - pas besoin de spécifier `--all` lors de l'utilisation de `--mars-year` :

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

**Ce qu'il fait :**

- Lit une table de recherche Excel/CSV avec l'année martienne, les plages Ls et les répertoires
- Trouve automatiquement le répertoire de départ pour votre année martienne + Ls
- **Recherche automatiquement dans tous les répertoires** (pas besoin de `--all`)
- Calcule les moyennes à partir de ce Ls
- Ajoute `MY34` au nom de fichier de sortie pour une identification facile

**Prérequis :**

```bash
pip install pandas openpyxl
```

**Format du fichier de recherche :**
Le fichier de recherche doit contenir les colonnes :

- `MY` - Année martienne (ex. : 34, 35)
- `Ls start` - Ls de départ pour cette plage
- `Ls end` - Ls de fin pour cette plage  
- `timestep start` - Premier pas de temps de la plage
- `timestep end` - Dernier pas de temps de la plage
- `directory start` - Répertoire contenant ces pas de temps

**Paramètres Mars Year :**

- `--mars-year MY` - Numéro d'année martienne (ex. : 34, 35) - **active automatiquement la recherche inter-répertoires**
- `--ls-start LS` - Ls de départ (ex. : 9.61, 30.0)
- `--lookup FILE` - Chemin vers le fichier de recherche Excel (.xlsx) ou CSV (requis avec --mars-year)

**Exemples de sorties :**

```
# Sans Mars Year
hl-b274_000000p_ls007_1234_sol000to004_5days_mean.nc

# Avec Mars Year
hl-b274_000000p_ls007_1234_MY34_sol000to004_5days_mean.nc
```

**Notes :**

- `--mars-year` active automatiquement la recherche inter-répertoires (recherche dans tous les répertoires les pas de temps correspondants)
- Peut être combiné avec `--max-means` pour limiter la sortie
- Peut être combiné avec `--ls-range` pour un filtrage supplémentaire
- `--lookup` est requis lors de l'utilisation de `--mars-year`

### Aide

Pour l'aide complète incluant les exemples :

```bash
python compute_diurnal_mean.py --help
```

---

## Ressources associées

### Ce projet

- **Code source :** [github.com/ludvdber/gem-mars-rpn2netcdf](https://github.com/ludvdber/gem-mars-rpn2netcdf)
- **Institution :** [Institut royal d'Aéronomie Spatiale de Belgique (BIRA-IASB)](https://www.aeronomie.be/)

### Modèle GEM-Mars et outils

- **Bibliothèque fstd2nc :** [github.com/neishm/fstd2nc](https://github.com/neishm/fstd2nc) - Lecture de fichiers RPN/FSTD
- **Documentation xarray :** [docs.xarray.dev](https://docs.xarray.dev/) - Manipulation NetCDF en Python

### NetCDF et conventions CF

- **Conventions CF-1.10 :** [cfconventions.org](http://cfconventions.org/) - Standard de métadonnées Climate and Forecast
- **Documentation NetCDF :** [unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf/) - Spécifications du format NetCDF

### Outils de visualisation

- **Paraview :** [paraview.org](https://www.paraview.org/) - Visualisation scientifique 3D
- **Panoply :** [giss.nasa.gov/tools/panoply](https://www.giss.nasa.gov/tools/panoply/) - Visualiseur NetCDF/HDF de la NASA

### Outils de performance

- **Numba :** [numba.pydata.org](https://numba.pydata.org/) - Compilateur JIT (utilisé dans compute_diurnal_mean.py)
- **Dask :** [dask.org](https://dask.org/) - Bibliothèque de calcul parallèle
