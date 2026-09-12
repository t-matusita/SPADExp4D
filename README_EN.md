# SPADExp4D

English | [日本語](README_JP.md)

Current version: **1.1.0**

SPADExp4D is research software for calculating and visualizing photoemission angular distributions (PADs) for angle-resolved photoemission spectroscopy (ARPES) in four-dimensional `(kx, ky, kz, E)` space using output from the first-principles code OpenMX.

It is a derivative of [SPADExp](https://github.com/Hiroaki-Tanaka-0606/SPADExp), developed by Hiroaki Tanaka and collaborators. SPADExp4D adds first-Brillouin-zone eigenstate caching, three-dimensional momentum grids, energy-dependent PADs, and calculation-result management.

## Main features

- Reusable HDF5 storage of eigenvalues and LCAO coefficients in the first Brillouin zone
- PAD calculations and cuts in `(kx, ky, kz, E)` space
- Support for the conventional OpenMX axis order: `a1/a2` in plane and `a3` along the surface normal
- CIF import with per-structure settings, caches, and result catalogs
- Batched C++ PAD engine and a Python backend for comparison
- Gaussian smoothing, Fermi-level shifts, polarization, Euler angles, and photon-energy controls
- Reuse of completed eigenstates and recovery after interrupted calculations
- Manual exchange workflow between a local computer and an external compute system

## Distribution scope

This repository distributes **source code only**. It does not include:

- compiled executables (`bin/`);
- build products (`.build/`);
- SPADExp/OpenMX HDF5 data;
- OpenMX PAO and VPS datasets; or
- calculation caches, logs, intermediate results, and user settings.

A clone is therefore not immediately runnable. Users must provide the software dependencies and runtime data described below.

## Platform and dependencies

SPADExp4D runs on **Linux and macOS**. It requires a working OpenMX installation, a C++17 build and runtime environment, and a working Python 3.12 environment.

This release has been operationally tested on macOS. The included Makefile detects the operating system and uses libomp with the Apple Accelerate framework on macOS, or the compiler's OpenMP support with BLAS/LAPACK on Linux.

### C/C++

- A C++17 compiler (Clang is currently assumed)
- GNU Make
- HDF5 built with its C++ API
- An OpenMP runtime (normally libomp on macOS)
- BLAS/LAPACK (the Apple Accelerate framework on macOS, or a suitable implementation on Linux)

### Python

- Python 3.12 (tested with 3.12.7)
- PyQt6 (tested with 6.11.0 and Qt 6.11.0)
- Matplotlib 3.11.1 (verified with PyQt6 using the `QtAgg` backend)
- NumPy
- h5py
- pytest, for running the test suite

### External software and runtime data

- The OpenMX executable, `openmx`
- `cif2omx`, included in the OpenMX source distribution and compiled separately
- OpenMX `DFT_DATA19`, or a compatible PAO/VPS dataset
- `data/Atomic_potentials.hdf5`
- `data/PAO_and_AO_after_opt.hdf5`

Obtain the HDF5 files from a legitimate source whose redistribution and usage terms you have verified, such as the original SPADExp distribution.

## Repository layout

```text
SPADExp4D/
├── cpp/             C++ engines and OpenMX conversion tools
├── python/          GUI, CLI, and 4D calculation control
├── structures/      CIF files and initial per-structure settings
├── tests/           Automated tests
├── exchange/        Manual external-compute exchange area
├── data/            Runtime HDF5 and OpenMX data, ignored by Git
├── bin/             Locally compiled executables, ignored by Git
├── .build/          Intermediate build products, ignored by Git
├── Makefile
└── LICENSE
```

## Setup (external software and this package are separate)

### A. Install separately

The following items are not included in this repository:

- OpenMX (compile its bundled `cif2omx` utility as described below)
- A C++17 compiler, GNU Make, HDF5 with its C++ API, OpenMP, and BLAS/LAPACK
- Python 3.12, PyQt6, Matplotlib 3.11.1, NumPy, and h5py (`pytest` for tests)

#### OpenMX and bundled cif2omx

Install OpenMX first. Then compile the bundled `cif2omx` in the OpenMX source directory:

```bash
cd /path/to/openmx/source
gcc -std=c99 -O2 -Wall -Wextra -pedantic \
  -DCIF2OPENMX_MAIN cif2omx.c -lm -o cif2omx
```

Place `openmx` and the generated `cif2omx` on `PATH`; SPADExp4D detects them automatically. Alternatively, specify them explicitly:

```bash
export SPADEXP4D_OPENMX=/path/to/openmx
export SPADEXP4D_CIF2OMX=/path/to/cif2omx
```

If neither an environment override nor a command on `PATH` is available, SPADExp4D reports an actionable error when the tool is required.

#### Python packages

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install PyQt6 'matplotlib==3.11.1' numpy h5py pytest
```

#### HDF5, OpenMP, and BLAS/LAPACK

Install these through the operating system or a package manager. If HDF5 is not found by `pkg-config`, specify its prefix when building:

```bash
make HDF5_PREFIX=/path/to/hdf5 LIBOMP_PREFIX=/path/to/libomp
```

### B. Runtime data supplied to this package

```text
data/
├── Atomic_potentials.hdf5
└── PAO_and_AO_after_opt.hdf5
```

If the OpenMX dataset is not available as `DFT_DATA19/` beside the OpenMX installation, it can be placed under `data/openmx/`.

### C. Build and run SPADExp4D

Run from the repository root:

```bash
make
make check
```

Executables are generated under `bin/`, and intermediate files under `.build/`. Neither directory is tracked by Git.

## Usage

### GUI

```bash
python3 python/_main.py
```

Built-in defaults are used if the initial settings file does not exist. Changes made in the GUI are saved in files such as the selected structure's `parameter.json`.

## Adding a crystal structure

Create `structures/<structure_name>/` and place a CIF file and `parameter.json` directly inside it. Structures can also be added with “Import CIF…” in the GUI.

```text
structures/<structure_name>/
├── structure.cif or another-name.cif
├── parameter.json
├── cache/        generated, ignored by Git
└── results/      generated, ignored by Git
```

Automatic processing does not support partial occupancy, disordered sites, or undefined elements. Complete-band eigenstate storage currently requires `scf.SpinPolarization off`.

## GUI operation

### 1. Select a calculation target and add a structure

After startup, select a structure from `計算対象` (calculation target). To add a new structure, choose `Import CIF…` from the menu, select a CIF file, and enter a structure name. Structure-specific input files and calculation results are saved under `structures/<structure_name>/`.

### 2. Prepare the first Brillouin zone

Set the first-Brillouin-zone divisions with `N1, N2, N3`. Larger divisions can improve sampling, but increase OpenMX time and HDF5 storage. Press `Calc` to calculate the required eigenstates with OpenMX. Existing caches with matching conditions are reused, and only missing portions are calculated.

Confirm the dialog with `Yes` to start. OpenMX logs and progress are shown while running; use `Cancel` to stop when necessary.

### 3. Set PAD display conditions

- `hν`: photon energy
- `V0`: inner potential; `V0 <= 0` uses the automatic setting
- `EF shift`: display-only Fermi-level shift
- `wfnc`: photoelectron wavefunction parameter
- `EB Begin, End, Step`: binding-energy range and step
- `偏光`: `lx`, `ly`, `c+`, `c-`, or `unpolarized`
- `Euler a, b, c`: rotation angles for polarization/crystal coordinates

Changing PAD conditions reuses the first-Brillouin-zone eigenstate cache. Select either `C++ 一括計算` (batched C++) or the comparison `Python EB毎` backend from `PAD計算`.

### 4. Set cuts and smoothing

Select `kx`, `ky`, or `kz` with `断面方向` and set the cut coordinates `a` and `b`. `KMin, KMax, Kstep` define the displayed momentum range and step. Enable `Brillouin zone表示` to overlay first-Brillouin-zone boundaries.

`Gaussian dK` and `Gaussian dE` apply momentum- and energy-direction Gaussian smoothing to saved PAD data. Changing smoothing does not require OpenMX or PAD recalculation.

### 5. Reuse and save results

Calculated structures, meshes, and PAD data are saved in the structure-specific `cache/` and `results/` directories. Pressing `Calc` again with the same conditions automatically reuses available caches. Use `計算結果を再読み込み` to refresh the result list.

Settings are saved in files such as each structure's `parameter.json`. For external-computer processing, create an upload job from the `受け渡し` menu and import returned results from the download area.

## Known limitations

- Linux and macOS are supported; this release has been operationally tested on macOS.
- Some Linux distributions or installation methods require explicit HDF5, OpenMP, or BLAS/LAPACK library names and search paths.
- Compiled binaries and runtime HDF5 files are not distributed.
- OpenMX and its data files must be provided separately.
- Large three-dimensional k-point grids can require substantial computation time and HDF5 storage.
- Automatic CIF conversion does not support partial occupancy, disordered sites, or undefined elements.
- Complete-band 4D caching for spin-polarized calculations is not currently supported.

## Citation

If SPADExp4D contributes to published research, please cite the original SPADExp paper:

H. Tanaka, K. Kuroda, and T. Matsushita, “SPADExp: A photoemission angular distribution simulator directly linked to first-principles calculations,” *Journal of Electron Spectroscopy and Related Phenomena* **264**, 147297 (2023). https://doi.org/10.1016/j.elspec.2023.147297

## Copyright, attribution, and license

Modifications made for SPADExp4D in 2026:

```text
Copyright (C) 2026 Matsushita Tomohiro
Nara Institute of Science and Technology (NAIST)
```

SPADExp4D is a modified work based on SPADExp by Hiroaki Tanaka and collaborators. Preserve the original copyright and license notices and all third-party notices embedded in the source. `sphere_lebedev_rule.cpp` contains code by John Burkardt carrying a GNU LGPL notice.

SPADExp4D is distributed under the same **GNU General Public License version 3** as SPADExp, without warranty. See [LICENSE](LICENSE) for the complete terms.

## Use of ChatGPT

During the 2026 development of SPADExp4D, OpenAI's ChatGPT was used as a supporting tool for code review, refactoring, testing, release preparation, and documentation drafting. ChatGPT is not a copyright holder or maintainer. Responsibility for the released material and research decisions remains with the copyright holder.

## Support and contributions

When reporting a problem, use this repository's issue tracker and include the operating system, Python/compiler/OpenMX versions, structure used, reproduction steps, and relevant error logs. Do not attach large HDF5 files, calculation results, personal information, credentials, or secrets.

Proposed changes should exclude generated files and include relevant tests and an explanation. Contributions must be distributable under this project's GNU General Public License version 3.
