# SPADExp4D

English | [日本語](README_JP.md)

Current version: **1.3.0**

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

Automatic processing does not support partial occupancy, disordered sites, or undefined elements. When a CIF cannot be used, `Import CIF` explains each problem in plain English and does not add the structure: missing or invalid cell values, no atom list or no fractional coordinates, type symbols that are not elements (for example placeholders such as `T` or `X1` from a structure template), partly occupied sites, several atoms at one position, and elements for which OpenMX has no basis set. cif2omx messages are under Show Details.

## GUI operation

### 1. Select a calculation target and add a structure

After startup, select a structure from `計算対象` (calculation target). To add a new structure, choose `Import CIF…` from the menu, select a CIF file, and enter a structure name. Structure-specific input files and calculation results are saved under `structures/<structure_name>/`.

### 2. Prepare the first Brillouin zone

Set the first-Brillouin-zone divisions with `N1, N2, N3`. Larger divisions can improve sampling, but increase OpenMX time and HDF5 storage. Press `Calc` to calculate the required eigenstates with OpenMX. Existing caches with matching conditions are reused, and only missing portions are calculated.

The `OpenMX Parameters` tab edits OpenMX keywords per calculation target:

- `scf.XcType`, `scf.SpinPolarization` (off / on / NC), `scf.SpinOrbit.Coupling` (NC only), `scf.energycutoff`, `scf.Kgrid`, `scf.ElectronicTemperature`, `scf.criterion`, `scf.maxIter`.
- Values start from the CIF-generated `.openmx.base.dat`. Changed values are shown in bold and stored under `"openmx"` in `parameter.json`; the next `Calc` then runs a new SCF and k-space HDF5 (each keyword set has its own SCF restart).
- `Reset to recommended` sets values derived from the cell (`scf.Kgrid` with a k-point spacing of at most 0.2 Å⁻¹, `scf.maxIter` ≥ 200, …), never weaker than the template.
- With `scf.SpinPolarization on` the PAD of both spins is summed; with `NC` both spinor components are summed. Use `NC` with `scf.SpinOrbit.Coupling` for spin-orbit coupling in a paramagnetic state.
- `Stored bands E−EF` limits the bands stored in the k-space HDF5 to an energy window (`All` stores every band). For NC, a window covering the EB range of interest reduces HDF5 size and memory several-fold; states outside it are absent from PAD.

Before OpenMX runs, the dialog “Calculate the eigenstates of <target> with OpenMX now?” shows the k mesh, any k points already saved by an interrupted run, and how many will be calculated. Choose `Yes` to start, or `No` to stop and, for example, change N1, N2, N3 first. OpenMX logs and progress are shown while running; use `Cancel` to stop when necessary.

The `Results` tab shows `K-space HDF5`: the k-space HDF5 file used for the PAD, when it was calculated, the OpenMX and HDF5 conversion times (total and wall time), and its size. The tooltip shows the full path. The times are stored in the k-space HDF5 file itself (attribute `CalculationTiming`), so they are also shown when the cache is reused, after a restart, and for an opened PAD result. A calculation that was interrupted and resumed adds up the time spent before the interruption. Caches made before this version have no recorded times and show only the date and size.

### 3. Set PAD conditions (`PAD` tab)

Below the OpenMX tabs, the settings are split into the `PAD` and `Band Diagram` tabs. The `PAD` tab lists, in order:

- `PAD calculation`: `C++ batch` or the comparison `Python per EB (comparison)` backend
- `hν`: photon energy
- `Polarization`: `lx`, `ly`, `c+`, `c-`, or `unpolarized`
- `Incident Angles (zyz Euler)`: Euler angles a, b, c (z-y-z rotation) of the incident light (a and c: −360 to 360°, b: 0 to 180°). Press the `Incident Angles (zyz Euler, °)…` button to choose a beamline preset (SPring-8 BL25SU: 0, 90, −90; UVSOR PMM: 0, 68, −90) and apply it with OK
- `EF shift`: display-only Fermi-level shift
- `V0`: inner potential; `V0 <= 0` uses the automatic setting
- `wfnc`: photoelectron wavefunction parameter
- `Slice direction`: `kx` or `ky`, the direction of the PAD and band slices
- `ky fixed` (`kx fixed` for a `ky` slice): the fixed coordinate that places the slice (see 4)
- `KMin, KMax, Kstep`: displayed momentum range and step
- `EB Begin, End, Step`: binding-energy range and step
- `Gaussian dK`, `Gaussian dE`: momentum- and energy-direction Gaussian smoothing of saved PAD data; changing them does not require OpenMX or PAD recalculation
- `Maximum angle θmax`: PAD outside this emission angle is not shown
- `Boundary`: overlay first-Brillouin-zone boundaries

Changing PAD conditions reuses the first-Brillouin-zone eigenstate cache.

### 4. Slice coordinates and the figure (`Band Diagram` tab)

The `Band Diagram` tab is reserved for band-diagram settings that are not implemented yet, so it is disabled in this version. It holds only the second fixed coordinate, which keeps the value saved in `parameter.json` (or in an opened result file); the first one is in the `PAD` tab under `Slice direction`. The fixed coordinates of the slice chosen with `Slice direction` are as follows. For a `kx` slice, `ky fixed` (`PAD` tab) sets ky; kz is not fixed but follows the photoemission sphere at each EB, kz = √((hν + V0 − EB)/(ħ²/2m) − kx² − ky²), so `kz fixed` is not used (likewise for a `ky` slice, with `kx fixed`). `Slice direction` offers only `kx` and `ky`; settings with a `kz` slice saved by an earlier version open as `kx`, and result files with a `kz` slice are still displayed.

Choose the figure layout with `Display type` above `Calc` (no recalculation needed).

| Display type | Panels |
|---|---|
| `Normal` (default) | 2×2: top left `Band Diagram EB′` (band contours at the selected EB), top right `PAD: EB′`, bottom left `Band Diagram (PAD cut)`, bottom right `PAD Slice kx-EB` |
| `With Photoelectron kz` | 3×2: `Band Slice EB′` and `PAD: EB′` (kx–ky maps at the selected EB), `Band Slice kx-EB` and `PAD Slice kx-EB` (the slice versus EB), `Band Diagram (PAD cut)` and the kx–kz photoelectron trajectories |

`Band Diagram EB′` draws, as lines, where the bands cross the photoemission sphere that the PAD map samples at the selected EB, kz = √((hν + V0 − EB)/(ħ²/2m) − kx² − ky²), i.e. where EB = −E_n(kx, ky, kz). The line colour is the PAD at that position (after the same Gaussian dK/dE as `PAD: EB′`), on a scale from the minimum to the maximum PAD shown at that EB (rainbow, colour bar `PAD`); lines where there is no PAD (outside the photoemission horizon or θmax) are grey. The contours are computed for every EB with each calculation and stored with the result; older results without them get them from the k-space HDF5, if it is available, when displayed. A thin red dotted line marks where the PAD slice runs: ky = `ky fixed` (horizontal) for a `kx` slice, kx = `kx fixed` (vertical) for a `ky` slice.

In both layouts, `Band Diagram (PAD cut)` is not a band structure at a fixed kz: it shows the bands at the k points the PAD slice actually samples, i.e. where EB = −E_n(kx, ky, kz(EB)), so it can be compared directly with `PAD Slice kx-EB`. Each band is coloured (rainbow) by the PAD at its position in the slice, after the Gaussian dK/dE, on a scale from the minimum to the maximum PAD of the slice in the displayed range (colour bar `PAD`). Bands where the slice has no PAD (outside the photoemission horizon or θmax) are grey. A thin red dotted line marks the current EB. All panels have the width and horizontal position of the kx–ky maps at the top, so the k axes line up in each column.

### 5. Reuse and save results

Calculated structures, meshes, and PAD data are saved in the structure-specific `cache/` and `results/` directories. Pressing `Calc` again with the same conditions automatically reuses available caches. Use `計算結果を再読み込み` to refresh the result list.

Settings are saved in files such as each structure's `parameter.json`. For external-computer processing, create an upload job from the `受け渡し` menu and import returned results from the download area.

#### Restoring the last session

On exit, the calculation target, every parameter in the sidebar, the selected tabs, `Display type`, the EB position and the window size are saved in `python/spadexp_4d.ini.json`. Each finished calculation is also written to `<calculation target>/cache/last_pad_result.h5`. At the next start, the same target is selected and its last PAD result is shown without recalculation. Gaussian dK/dE, EF shift, θmax and boundary changes made after the calculation are kept. The sidebar parameters come back exactly as they were shown at exit, even a combination that cannot be used yet (for example EB End − EB Begin not a whole multiple of EB Step); the reason is then shown in red above `Calc`, which refuses it until fixed, and `parameter.json` keeps the last valid values. Switching targets uses the parameters saved for that target; only choosing an entry under `Calculation result` loads the conditions of that calculation. If a PAD result file was open at exit, that file is shown instead.

#### Saving and opening PAD results (HDF5)

`File > Save PAD result as HDF5…` writes the displayed result to one HDF5 file. `File > Open PAD result (HDF5)…` shows it again, selects its calculation target and restores its parameters. If the target is not in the data folder, it can be added from the crystal structure stored in the file. The file contains:

| Group | Contents |
|---|---|
| `/` (attributes) | `format = "SPADExp4D PAD result"`, `format_version`, `software`, `created`, `structure_name` |
| `Settings` | all GUI and OpenMX settings as JSON (attribute `json`), plus each scalar setting as an attribute |
| `PAD` | `binding_energy` (eV), `k` (Å⁻¹, shared by kx and ky), unbroadened `pad` and `band_weight` (EB, ky, kx; gzip), `nearest_eigenstate_distance`; attributes `v0_eV`, `request_json` |
| `BandCut` | bands along the PAD slice: `k` and `energy` (E − EF, eV; NaN outside the slice) for lines, `points` (k, E − EF) where a band meets one k more than once or for a kz slice |
| `Structure` | `cif`, `openmx_template` (`.openmx.base.dat`), `lattice` (Å), `reciprocal` (Å⁻¹) |
| `BandContours` | band contours at every EB (`Band Diagram EB′`): `binding_energy`, `offsets` (EB i owns segments `offsets[i]:offsets[i+1]`), `segments` (kx, ky of both segment ends, 1/Å) |
| `Calculation` | k-space HDF5 path and size, OpenMX/conversion/PAD times, the k-space HDF5 calculation record (`mesh_timing`, JSON), maximum sampling distance |

Gaussian dK/dE are stored as settings and applied when the file is shown, so the stored PAD stays unbroadened. The k-space HDF5 (eigenstates) is not included; a new `Calc` with changed conditions needs it or recomputes it.

## Known limitations

- Linux and macOS are supported; this release has been operationally tested on macOS.
- Some Linux distributions or installation methods require explicit HDF5, OpenMP, or BLAS/LAPACK library names and search paths.
- Compiled binaries and runtime HDF5 files are not distributed.
- OpenMX and its data files must be provided separately.
- Large three-dimensional k-point grids can require substantial computation time and HDF5 storage.
- Automatic CIF conversion does not support partial occupancy, disordered sites, or undefined elements.
- `scf.SpinPolarization on` keeps the initial spin populations written by cif2omx (equal up and down). Magnetic order needs initial spins or a magnetic cell in `.openmx.base.dat`.

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
