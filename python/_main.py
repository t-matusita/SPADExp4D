#!/usr/bin/env python3
"""4-D (kx, ky, kz, E) viewer using a reusable first-BZ eigenstate database."""

# SPADExp4D 1.2.0: integrated Python, C++ and per-Structure caches.

from __future__ import annotations

__version__ = "1.2.0"

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
import threading
import tempfile
from datetime import datetime

import spadexp_4d_core as core

import h5py
import numpy as np

try:
    from PyQt6 import QtCore, QtGui, QtWidgets
except ImportError as exc:
    raise ImportError(
        "PyQt6 is required. Install PyQt6 in the Python environment used to run SPADExp4D."
    ) from exc

import matplotlib

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

# The application and its C++ helpers live together; large runtime data stays
# at the project root so existing first-BZ and PAD caches remain reusable.
APP_DIR = Path(__file__).resolve().parent
PROJECT = APP_DIR.parent
INITIAL_SETTINGS_FILE = APP_DIR / "spadexp_4d.ini.json"
STRUCTURES = PROJECT / "structures"
DEFAULT_STRUCTURE = STRUCTURES / "graphite"
EXCHANGE = PROJECT / "exchange"
UPLOAD = EXCHANGE / "upload"
DOWNLOAD = EXCHANGE / "download"
DATA = PROJECT / "data"
BUILD = PROJECT / ".build"

SPADEXP = PROJECT / "bin/spadexp_standard_axes"
SPADEXP_4D = PROJECT / "bin/spadexp_4d_plus"
ATOMIC_ORBITALS = DATA / "PAO_and_AO_after_opt.hdf5"
ATOMIC_POTENTIALS = DATA / "Atomic_potentials.hdf5"

ACTIVE_STRUCTURE = DEFAULT_STRUCTURE
CACHE = ACTIVE_STRUCTURE / "cache"
INPUT = CACHE / "tmp"
OUTPUT = CACHE
SCF_H5 = CACHE / "graphite_2d.hdf5"
PAD_VOLUME = CACHE / "pad_volume.npz"

LATTICE = core.LATTICE
B = core.G
PAD_EH = 27.2114  # Upstream PAD engine energy conversion (setup.cpp).
EH = 27.211386245988
H2M = 3.80998211615486  # hbar^2/(2m), eV Angstrom^2

K_HALF, K_LIMIT, DK = 3.0, 10.0, 0.075
POLARIZATIONS = ("lx", "ly", "c+", "c-", "unpolarized")
# Beamline presets for the Euler angles (a, b, c) in degrees, offered by the
# "Incident Angles (zyz Euler)" button.  a and c may be negative (-360 to 360).
EULER_PRESETS = (
    ("SPring-8 BL25SU", (0.0, 90.0, -90.0)),
    ("UVSOR PMM", (0.0, 68.0, -90.0)),
)

DEFAULT_SETTINGS = {
    "hv": 21.2, "v0": 20.0, "ef_shift": 0.0, "wfnc": 4.0,
    "eb_min": 0.0, "eb_max": 10.0, "de": 0.1,
    "direction": "kx", "a": 0.0, "b": 0.0,
    "dk": 0.2, "dE": 0.2, "theta_max": 80.0,
    "polarization": "ly", "alpha": 0.0, "beta": 0.0, "gamma": 0.0,
    "k_min": -3.0, "k_max": 3.0, "k_step": 0.075,
    "show_bz": True, "use_first_bz_only": True,
    "N1": 24, "N2": 24, "N3": 12, "pad_engine": "cpp",
}

SETTING_KEYS = (
    "hv", "v0", "ef_shift", "wfnc", "eb_min", "eb_max", "de", "direction", "a", "b",
    "dk", "dE", "theta_max", "polarization", "alpha", "beta", "gamma",
    "k_min", "k_max", "k_step", "show_bz", "use_first_bz_only", "N1", "N2", "N3",
)

# Application-level keys kept only in INITIAL_SETTINGS_FILE.  They are not PAD
# parameters, so validate_settings() drops them and write_settings() restores
# them when it rewrites the initialization file.
DATA_FOLDER_KEY = "data_folder"
# Window state restored at the next start: calculation target, EB position,
# selected tabs and window geometry (see Window.save_session).
SESSION_KEY = "session"
APP_SETTING_KEYS = (DATA_FOLDER_KEY, SESSION_KEY)
DEFAULT_DATA_FOLDER = PROJECT / "structures"

# Self-contained PAD result files (HDF5).  Every completed calculation is also
# written to <calculation target>/cache/LAST_PAD_RESULT and shown again at the
# next start.
PAD_RESULT_FORMAT = "SPADExp4D PAD result"
PAD_RESULT_VERSION = 1
LAST_PAD_RESULT = "last_pad_result.h5"
PAD_RESULT_FILTER = "SPADExp4D PAD result (*.h5 *.hdf5);;All files (*)"
CALCULATION_KEYS = ("mesh_path", "mesh_bytes", "openmx_seconds", "conversion_seconds",
                    "openmx_hdf5_seconds", "mesh_seconds", "pad_seconds", "pad_kernel_seconds",
                    "pad_computed_slices", "result_record")


def _is_initial_settings_file(path):
    try:
        return Path(path).expanduser().resolve() == Path(INITIAL_SETTINGS_FILE).resolve()
    except OSError:
        return False


def _read_initial_settings_raw():
    """Return the initialization file as a dict, or {} when absent or invalid."""
    try:
        values = json.loads(Path(INITIAL_SETTINGS_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return values if isinstance(values, dict) else {}


def configured_data_folder():
    """Return the data folder stored in the initialization file.

    A relative path is taken relative to the project folder.  Without a stored
    value the bundled ``structures`` folder is used.
    """
    value = _read_initial_settings_raw().get(DATA_FOLDER_KEY)
    if not isinstance(value, str) or not value.strip():
        return DEFAULT_DATA_FOLDER
    path = Path(value.strip()).expanduser()
    return path if path.is_absolute() else PROJECT / path


def save_data_folder(path):
    """Store the data folder in the initialization file, keeping PAD defaults."""
    values = _read_initial_settings_raw()
    if any(key not in values for key in SETTING_KEYS):
        # read_settings() requires a complete parameter set once the file exists.
        values = {**validate_settings(DEFAULT_SETTINGS), **values}
    values[DATA_FOLDER_KEY] = str(Path(path).expanduser().resolve())
    core.json_write(INITIAL_SETTINGS_FILE, values)


def read_session():
    """Return the window state saved at the last exit ({} when absent)."""
    state = _read_initial_settings_raw().get(SESSION_KEY)
    return state if isinstance(state, dict) else {}


def save_session(settings, state):
    """Store the GUI parameters and window state for the next start."""
    write_settings(INITIAL_SETTINGS_FILE, settings)
    values = _read_initial_settings_raw()
    values[SESSION_KEY] = state
    core.json_write(INITIAL_SETTINGS_FILE, values)


def set_data_folder(path):
    """Point structure discovery and every OpenMX/PAD output at one folder."""
    global STRUCTURES, DEFAULT_STRUCTURE
    STRUCTURES = Path(path).expanduser().resolve()
    DEFAULT_STRUCTURE = STRUCTURES / "graphite"
    core.set_structures_root(STRUCTURES)
    configure_runtime_paths(DEFAULT_STRUCTURE)
    return STRUCTURES



def validate_settings(settings):
    """Return a normalized settings dictionary or raise a readable error."""
    if not isinstance(settings, dict):
        raise ValueError("The top level of the settings file must be a JSON object.")
    settings = dict(settings)
    settings.setdefault("ef_shift", DEFAULT_SETTINGS["ef_shift"])
    missing = [key for key in SETTING_KEYS if key not in settings]
    if missing:
        raise ValueError("Missing settings keys: " + ", ".join(missing))
    try:
        normalized = {
            key: float(settings[key])
            for key in (
                "hv", "v0", "ef_shift", "wfnc", "eb_min", "eb_max", "de", "a", "b",
                "dk", "dE", "theta_max", "alpha", "beta", "gamma",
                "k_min", "k_max", "k_step",
            )
        }
    except (TypeError, ValueError) as exc:
        raise ValueError("Check the numeric settings.") from exc
    normalized["direction"] = str(settings["direction"])
    normalized["polarization"] = str(settings["polarization"])
    normalized["show_bz"] = settings["show_bz"]
    normalized["use_first_bz_only"] = settings["use_first_bz_only"]
    if normalized["direction"] not in ("kx", "ky", "kz"):
        raise ValueError("direction must be one of kx, ky, or kz.")
    if normalized["polarization"] not in POLARIZATIONS:
        raise ValueError("polarization must be one of lx, ly, c+, c-, or unpolarized.")
    if not isinstance(normalized["show_bz"], bool):
        raise ValueError("show_bz must be true or false.")
    if not isinstance(normalized["use_first_bz_only"], bool):
        raise ValueError("use_first_bz_only must be true or false.")
    if not all(np.isfinite(normalized[key]) for key in normalized
               if key not in ("direction", "polarization", "show_bz",
                              "use_first_bz_only")):
        raise ValueError("Numeric settings must be finite.")
    if normalized["hv"] <= 0:
        raise ValueError("hv must be positive.")
    if normalized["wfnc"] < 0:
        raise ValueError("wfnc must be non-negative.")
    if normalized["eb_max"] <= normalized["eb_min"]:
        raise ValueError("Set eb_min less than eb_max.")
    if normalized["de"] <= 0:
        raise ValueError("de must be positive.")
    if normalized["dk"] < 0 or normalized["dE"] < 0:
        raise ValueError("Gaussian widths dk and dE must be non-negative.")
    if not 0 < normalized["theta_max"] <= 90:
        raise ValueError("theta_max must be greater than 0 and at most 90 degrees.")
    if not -360 <= normalized["alpha"] <= 360 or not -360 <= normalized["gamma"] <= 360:
        raise ValueError("Euler angles alpha and gamma must be between -360 and 360 degrees.")
    if not 0 <= normalized["beta"] <= 180:
        raise ValueError("Euler angle beta must be between 0 and 180 degrees.")
    if not (-K_LIMIT <= normalized["k_min"] < normalized["k_max"] <= K_LIMIT):
        raise ValueError(f"Set the k range to {-K_LIMIT:g} <= Min < Max <= {K_LIMIT:g}.")
    if normalized["k_step"] <= 0:
        raise ValueError("k step must be positive.")
    k_intervals = (normalized["k_max"] - normalized["k_min"]) / normalized["k_step"]
    if k_intervals < 1 or not np.isclose(k_intervals, round(k_intervals), rtol=0, atol=1e-7):
        raise ValueError("k Max - k Min must be an integer multiple of k step.")
    intervals = (normalized["eb_max"] - normalized["eb_min"]) / normalized["de"]
    if intervals < 1 or not np.isclose(intervals, round(intervals), rtol=0, atol=1e-7):
        raise ValueError("eb_max - eb_min must be an integer multiple of de.")
    normalized.update(zip(("N1", "N2", "N3"), core.divisions(settings)))
    normalized["use_first_bz_only"] = True
    normalized["pad_engine"] = settings.get("pad_engine", "cpp")
    if normalized["pad_engine"] not in ("cpp", "python"):
        raise ValueError("pad_engine must be cpp or python.")
    # OpenMX keyword overrides on top of the target's CIF template (optional).
    normalized["openmx"] = core.validate_openmx_overrides(settings.get("openmx"))
    window = core.validate_band_window(settings.get("band_window"))
    normalized["band_window"] = None if window is None else list(window)
    return normalized

def read_settings(path=INITIAL_SETTINGS_FILE):
    try:
        values = json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        if Path(path) == INITIAL_SETTINGS_FILE:
            return validate_settings(DEFAULT_SETTINGS)
        raise FileNotFoundError(f"Initial settings file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in settings file ({path}): {exc}") from exc
    if not isinstance(values, dict):
        raise ValueError("The top level of the settings file must be a JSON object.")
    # Backward-compatible migration from separate kx/ky display limits.
    if "k_min" not in values:
        values["k_min"] = min(values.get("kx_min", -K_HALF), values.get("ky_min", -K_HALF))
    if "k_max" not in values:
        values["k_max"] = max(values.get("kx_max", K_HALF), values.get("ky_max", K_HALF))
    values.setdefault("k_step", DK)
    values.setdefault("ef_shift", DEFAULT_SETTINGS["ef_shift"])
    values.setdefault("use_first_bz_only", True)
    for key in ("N1", "N2", "N3"):
        values.setdefault(key, DEFAULT_SETTINGS[key])
    return validate_settings(values)

def write_settings(path, settings):
    values = validate_settings(settings)
    if _is_initial_settings_file(path):
        existing = _read_initial_settings_raw()
        values.update({key: existing[key] for key in APP_SETTING_KEYS if key in existing})
    core.json_write(path, values)

def configure_runtime_paths(structure_path=None):
    """Point every transient PAD path at one Structure's private cache."""
    global ACTIVE_STRUCTURE, CACHE, INPUT, OUTPUT, SCF_H5, PAD_VOLUME
    ACTIVE_STRUCTURE = Path(structure_path or DEFAULT_STRUCTURE).resolve()
    CACHE = ACTIVE_STRUCTURE / "cache"
    INPUT = CACHE / "tmp"
    OUTPUT = CACHE
    SCF_H5 = CACHE / "graphite_2d.hdf5"
    PAD_VOLUME = CACHE / "pad_volume.npz"
    return CACHE


def mkdirs():
    for path in (INPUT, CACHE):
        path.mkdir(parents=True, exist_ok=True)


def ensure_app_folders():
    """Create the user-facing structure and transfer folders on startup."""
    for path in (STRUCTURES, UPLOAD, DOWNLOAD):
        path.mkdir(parents=True, exist_ok=True)


def safe_name(value):
    value = str(value).strip()
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)
    return cleaned.strip("._-") or "structure"


def structure_cif_path(structure_path):
    """Return the CIF supplied in a structure folder.

    ``structure.cif`` is the canonical export name used by the exchange
    format, but users should also be able to drop a naturally named file such
    as ``TaS2.cif`` into the folder.  The canonical file wins when both are
    present; otherwise the first CIF in lexical order is used.
    """
    structure_path = Path(structure_path)
    canonical = structure_path / "structure.cif"
    if canonical.is_file():
        return canonical
    candidates = sorted(path for path in structure_path.glob("*.cif") if path.is_file())
    return candidates[0] if candidates else None


def structure_entries():
    ensure_app_folders()
    entries = []
    for path in sorted(STRUCTURES.iterdir()):
        if path.is_dir() and structure_cif_path(path) is not None:
            entries.append(path)
    return entries


def latest_parameter_path():
    """Return the most recently saved structure parameter file."""
    candidates = [path / "parameter.json" for path in structure_entries()
                  if (path / "parameter.json").is_file()]
    if candidates:
        return max(candidates, key=lambda path: path.stat().st_mtime_ns)
    return INITIAL_SETTINGS_FILE


def result_entries(structure_path):
    results = structure_path / "results"
    if not results.is_dir():
        return []
    prune_duplicate_local_results(structure_path)
    return [path for path in sorted(results.iterdir())
            if path.is_dir() and (path / "manifest.json").is_file()]


def removable_local_result(path):
    """Return the OpenMX HDF5 for a metadata-only local result record."""
    path = Path(path); manifest = core.json_read(path / "manifest.json")
    link = path / "eigenstates.hdf5"
    allowed = {"config.json", "eigenstates.hdf5", "manifest.json"}
    if (not manifest or manifest.get("source") != "local" or not link.is_symlink()
            or any(child.name not in allowed for child in path.iterdir())):
        return None
    try:
        return link.resolve(strict=True)
    except OSError:
        return None


def prune_duplicate_local_results(structure_path):
    """Keep one catalog entry per OpenMX HDF5 and remove PAD-only duplicates."""
    results = Path(structure_path) / "results"
    if not results.is_dir(): return 0
    active = core.json_read(results / "active.json") or {}
    active_name = active.get("path")
    groups = {}
    for path in results.iterdir():
        if not path.is_dir(): continue
        mesh = removable_local_result(path)
        if mesh is not None:
            groups.setdefault(str(mesh), []).append(path)
    removed = 0
    for records in groups.values():
        if len(records) < 2: continue
        preferred = next((path for path in records if path.name == active_name), None)
        keep = preferred or max(records, key=lambda path: path.stat().st_mtime_ns)
        for path in records:
            if path != keep:
                shutil.rmtree(path); removed += 1
    return removed


def register_local_result(structure_path, settings, mesh_path, result=None):
    """Register a completed local cache in the structure's result catalog."""
    structure_path = Path(structure_path)
    mesh_path = Path(mesh_path).resolve()
    cif = structure_cif_path(structure_path)
    if cif is None or not mesh_path.is_file():
        return None
    settings = validate_settings(settings)
    cif_hash = sha256_file(cif)
    identity = json.dumps({"structure_name": structure_path.name,
                           "structure_cif_sha256": cif_hash,
                           "mesh": core.stamp(mesh_path)},
                          sort_keys=True, ensure_ascii=False)
    code = hashlib.sha256(identity.encode()).hexdigest()[:10]
    ns = core.divisions(settings)
    calculation_id = f"openmx_N{ns[0]}x{ns[1]}x{ns[2]}_{code}"
    target = structure_path / "results" / calculation_id
    target.mkdir(parents=True, exist_ok=True)
    write_settings(target / "config.json", settings)
    config_text = json.dumps(validate_settings(settings), sort_keys=True, ensure_ascii=False)
    link = target / "eigenstates.hdf5"
    if not link.exists() and not link.is_symlink():
        link.symlink_to(os.path.relpath(mesh_path, target))
    manifest = {
        "schema": "spadexp-4d-result", "schema_version": 1,
        "source": "local", "structure_name": structure_path.name,
        "calculation_id": calculation_id,
        "cache_kind": "openmx-eigenstates",
        "pad_cache_policy": "latest-volume-only",
        "config_sha256": hashlib.sha256(config_text.encode()).hexdigest(),
        "structure_cif_sha256": cif_hash,
        "eigenstates": "eigenstates.hdf5",
        "mesh_path": str(mesh_path), "mesh_bytes": mesh_path.stat().st_size,
    }
    if result:
        manifest["timing_seconds"] = {
            "openmx": float(result.get("openmx_seconds", 0.0)),
            "hdf5_conversion": float(result.get("conversion_seconds", 0.0)),
            "total": float(result.get("openmx_hdf5_seconds", 0.0)),
            "pad": float(result.get("pad_seconds", 0.0)),
        }
    core.json_write(target / "manifest.json", manifest)
    core.json_write(structure_path / "results" / "active.json", {
        "structure_name": structure_path.name,
        "calculation_id": calculation_id,
        "path": calculation_id,
        "eigenstates": "eigenstates.hdf5",
        "manifest": "manifest.json",
    })
    prune_duplicate_local_results(structure_path)
    return target


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def format_file_size(size):
    """Format a byte count with an IEC unit and exact bytes."""
    size = max(int(size), 0)
    value = float(size); unit = 'B'
    for candidate in ('KiB', 'MiB', 'GiB', 'TiB'):
        if value < 1024: break
        value /= 1024; unit = candidate
    compact = f"{int(value)} {unit}" if unit == 'B' else f"{value:.2f} {unit}"
    return f"{compact} ({size:,} bytes)"


def format_duration(seconds):
    """Format an elapsed time compactly for the status and log windows."""
    seconds = max(float(seconds), 0.0)
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    minutes, remainder = divmod(seconds, 60.0)
    if minutes < 60.0:
        return f"{int(minutes)}m {remainder:04.1f}s"
    hours, minutes = divmod(int(minutes), 60)
    return f"{hours}h {minutes}m {remainder:04.1f}s"

def automatic_v0(h5, requested):
    if requested > 0:
        return requested
    with h5py.File(h5) as f:
        ef = float(f["Output"].attrs["EF_Eh"])
        bottoms = [np.min(f[f"Output/{x}"]) for x in ("Band", "BandUp", "BandDn")
                   if x in f["Output"]]
    return (ef - min(bottoms)) * EH

def rotation_z(angle):
    """Active Cartesian rotation about z by an angle in degrees."""
    angle = math.radians(angle); c, s = math.cos(angle), math.sin(angle)
    return np.array(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)))

def rotation_y(angle):
    """Active Cartesian rotation about y by an angle in degrees."""
    angle = math.radians(angle); c, s = math.cos(angle), math.sin(angle)
    return np.array(((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c)))

def light_vectors(settings):
    """Return (complex polarization, propagation) after a z-y-z Euler rotation."""
    root2 = math.sqrt(2.0)
    basis = {
        "lx": np.array((1.0 + 0j, 0j, 0j)),
        "ly": np.array((0j, 1.0 + 0j, 0j)),
        "c+": np.array((1.0 / root2, 1j / root2, 0j)),
        "c-": np.array((1.0 / root2, -1j / root2, 0j)),
    }[settings["polarization"]]
    rotation = (rotation_z(settings["alpha"])
                @ rotation_y(settings["beta"])
                @ rotation_z(settings["gamma"]))
    field = rotation @ basis
    propagation = rotation @ np.array((0.0, 0.0, 1.0))
    return field, propagation

def polarization_components(settings):
    """Return the independently cached PAD components for a GUI selection."""
    selected = settings["polarization"]
    return ("c+", "c-") if selected == "unpolarized" else (selected,)

def sum_pad_components(first, second):
    """Add two polarization intensities while retaining polarization-free bands."""
    if first['eb'] != second['eb'] or first['pad'].shape != second['pad'].shape:
        raise ValueError('PAD grids for polarization components do not match.')
    result = dict(first)
    result['pad'] = first['pad'] + second['pad']
    result['error'] = max(first['error'], second['error'])
    return result

def field_request(field):
    return [[round(float(value.real), 15), round(float(value.imag), 15)]
            for value in field]

# The SPADExp engine reads file names with sscanf("%s"), which stops at the
# first space: a data folder such as "/Volumes/My Disk/..." truncated
# Output_file and the run failed.  The engine therefore runs inside its
# temporary directory and only sees these space-free relative names.
ENGINE_POTENTIALS = "Atomic_potentials.hdf5"
ENGINE_ORBITALS = "PAO_and_AO_after_opt.hdf5"
ENGINE_LOG = "pad.log"


def prepare_engine_directory(directory):
    """Link the shared atomic data into the engine's run directory."""
    for name, target in ((ENGINE_POTENTIALS, ATOMIC_POTENTIALS), (ENGINE_ORBITALS, ATOMIC_ORBITALS)):
        link = Path(directory) / name
        if not link.exists():
            link.symlink_to(Path(target).resolve())


def keep_engine_log(directory, destination):
    """Copy the engine log from the run directory to the structure cache."""
    log = Path(directory) / ENGINE_LOG
    if log.is_file():
        shutil.copyfile(log, destination)


def pad_text(source, destination, logfile, hv, v0, wf, eb_min, eb_max, de,
             gaussian_de, field):
    """PAD input with file names relative to the engine's run directory."""
    for value in (source, destination, logfile):
        if any(ch.isspace() for ch in str(value)):
            raise ValueError(f"SPADExp engine file names must not contain spaces: {value}")
    ex, ey, ez = field
    return f"""&Control
Calculation PAD
Log_file {logfile}
Console_log true
Output_file {destination}
/
&Atomic-wfn
Potential_file {ENGINE_POTENTIALS}
/
&PAD
Input_file {source}
E_min {-eb_max}
E_max {-eb_min}
E_pixel {de}
dE {gaussian_de}
Initial_state AO
Final_state Calc
Final_state_step 0.01
Photon_energy {hv}
Inner_potential {v0}
Work_function {wf}
Electric_field {ex.real} {ex.imag} {ey.real} {ey.imag} {ez.real} {ez.imag}
Atomic_orbitals_file {ENGINE_ORBITALS}
/
&Radial-grid 3
0.0025 40
0.005 40
0.01 40000
/
&Phase-shift
Skip_points 10
Calc_points 10
/
"""

def pad_ok(path, hv, v0, wf, dimension, eb_min, eb_max, de, gaussian_de,
           field):
    try:
        with h5py.File(path) as f:
            return (int(f.attrs["Dimension"]) == dimension
                    and np.isclose(f.attrs["PhotonEnergy_eV"], hv)
                    and np.isclose(f.attrs["InnerPotentialInput_eV"], v0)
                    and np.isclose(f.attrs["WorkFunction_eV"], wf)
                    and np.isclose(f.attrs["Offset"][-1], -eb_max)
                    and np.isclose(f.attrs["Delta"][-1], de)
                    and np.isclose(f.attrs["dE"], gaussian_de)
                    and np.allclose(f.attrs["ElectricFieldReal"], field.real)
                    and np.allclose(f.attrs["ElectricFieldImag"], field.imag)
                    and np.isclose(
                        f.attrs["Offset"][-1]
                        + f.attrs["Delta"][-1] * (f["Dispersion"].shape[-1] - 1),
                        -eb_min,
                    ))
    except (OSError, KeyError):
        return False

def frac(kcart, reciprocal=None):
    """Cartesian (1/A) -> fractional; pass the grid's OpenMX cell for engine input."""
    return np.asarray(kcart) @ np.linalg.inv(B if reciprocal is None else reciprocal)

def gaussian_smooth_axis(values, width, step, axis, mode="edge", stop=None):
    """Convolve one array axis with a normalized Gaussian in physical units."""
    values = np.asarray(values, dtype=float)
    if values.shape[axis] < 2 or width <= 0 or step == 0:
        return values
    sigma_pixels = width / abs(step)
    radius = min(max(1, int(math.ceil(4 * sigma_pixels))), values.shape[axis] - 1)
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / sigma_pixels) ** 2)
    kernel /= kernel.sum()

    padding = [(0, 0)] * values.ndim; padding[axis] = (radius, radius)
    padded = np.pad(values, padding, mode=mode)
    core.check_stop(stop)
    # smoothed[i] = sum_j kernel[j] * padded[i + j] along the axis, done as one
    # banded matrix product (BLAS) instead of a loop over the kernel taps.
    count = values.shape[axis]; rows = np.arange(count)
    band = np.zeros((count, count + 2 * radius))
    for j, weight in enumerate(kernel):
        band[rows, rows + j] = weight
    smoothed = np.tensordot(band, np.moveaxis(padded, axis, 0), axes=(1, 0))
    core.check_stop(stop)
    return np.ascontiguousarray(np.moveaxis(smoothed, 0, axis))


def gaussian_smooth_k_map(values, dk, extent):
    """Apply the momentum-width Gaussian to both axes of a kx-ky map."""
    dx = (extent[1] - extent[0]) / max(values.shape[1] - 1, 1)
    dy = (extent[3] - extent[2]) / max(values.shape[0] - 1, 1)
    smoothed = gaussian_smooth_axis(values, dk, dx, 1)
    return gaussian_smooth_axis(smoothed, dk, dy, 0)

def bz_polygon():
    poly = [np.array(p, float) for p in ((-20,-20),(20,-20),(20,20),(-20,20))]
    for n1 in range(-2, 3):
        for n2 in range(-2, 3):
            if n1 == n2 == 0: continue
            g = n1 * B[0,:2] + n2 * B[1,:2]; limit = np.dot(g,g)/2
            new = []
            for p, q in zip(poly, poly[1:] + poly[:1]):
                dp, dq = np.dot(p,g)-limit, np.dot(q,g)-limit
                if dp <= 1e-12: new.append(p)
                if dp * dq < 0: new.append(p + (q-p)*(-dp)/(dq-dp))
            poly = new
    return np.asarray(poly)

def repeated_bz_polygons(xlim, ylim, margin=1):
    """Yield translated first-BZ cells covering a Cartesian plotting window."""
    # B changes when the selected structure changes.  Build the Wigner-Seitz
    # cell from the current reciprocal lattice as well; retaining the cell
    # computed at module import mixes Graphite's polygon with the selected
    # material's translations.
    cell = bz_polygon()
    corners = np.array(((xlim[0], ylim[0], 0.0), (xlim[0], ylim[1], 0.0),
                        (xlim[1], ylim[0], 0.0), (xlim[1], ylim[1], 0.0)))
    f = frac(corners)
    n1 = range(math.floor(f[:, 0].min()) - margin,
               math.ceil(f[:, 0].max()) + margin + 1)
    n2 = range(math.floor(f[:, 1].min()) - margin,
               math.ceil(f[:, 1].max()) + margin + 1)
    for i in n1:
        for j in n2:
            shift = i * B[0, :2] + j * B[1, :2]
            polygon = cell + shift
            if (polygon[:, 0].max() >= xlim[0] and polygon[:, 0].min() <= xlim[1]
                    and polygon[:, 1].max() >= ylim[0] and polygon[:, 1].min() <= ylim[1]):
                yield i, j, polygon

def bz_crossings_at_y(y_value, xlim):
    """Return repeated-zone boundary intersections with a horizontal k line."""
    crossings = []
    for _, _, poly in repeated_bz_polygons(xlim, (y_value - 1e-6, y_value + 1e-6), 2):
        closed = np.vstack((poly, poly[0]))
        for p, q in zip(closed[:-1], closed[1:]):
            if min(p[1], q[1]) - 1e-10 <= y_value <= max(p[1], q[1]) + 1e-10:
                if abs(q[1] - p[1]) < 1e-12:
                    crossings.extend((p[0], q[0]))
                else:
                    t = (y_value - p[1]) / (q[1] - p[1])
                    if -1e-10 <= t <= 1 + 1e-10:
                        crossings.append(p[0] + t * (q[0] - p[0]))
    return sorted({round(x, 9) for x in crossings if xlim[0] <= x <= xlim[1]})

def bz_crossings_at_x(x_value, ylim):
    """Return repeated-zone boundary intersections with a vertical k line."""
    crossings = []
    for _, _, poly in repeated_bz_polygons((x_value - 1e-6, x_value + 1e-6), ylim, 2):
        closed = np.vstack((poly, poly[0]))
        for p, q in zip(closed[:-1], closed[1:]):
            if min(p[0], q[0]) - 1e-10 <= x_value <= max(p[0], q[0]) + 1e-10:
                if abs(q[0] - p[0]) < 1e-12:
                    crossings.extend((p[1], q[1]))
                else:
                    t = (x_value - p[0]) / (q[0] - p[0])
                    if -1e-10 <= t <= 1 + 1e-10:
                        crossings.append(p[1] + t * (q[1] - p[1]))
    return sorted({round(y, 9) for y in crossings if ylim[0] <= y <= ylim[1]})


def trajectory_bz_boundaries(xlim, kz_max):
    """Return repeated-zone boundaries visible in the positive-kz trajectory panel."""
    kx_bounds = bz_crossings_at_y(0.0, xlim)
    period = np.linalg.norm(B[2])
    count = int(math.ceil(kz_max / period)) + 1 if period > 0 else 0
    kz_bounds = [(n + 0.5) * period for n in range(count)
                 if (n + 0.5) * period <= kz_max]
    return kx_bounds, kz_bounds

class Canvas(FigureCanvasQTAgg):
    def __init__(self): self.figure = Figure(figsize=(16,13), constrained_layout=True); super().__init__(self.figure)

    def print_figure(self, *args, **kwargs):
        # EB-dependent panels are "animated" for fast slider updates, and
        # matplotlib leaves animated axes out of saved images.  Save them all.
        animated = [artist for artist in self.figure.findobj() if artist.get_animated()]
        for artist in animated: artist.set_animated(False)
        try:
            return super().print_figure(*args, **kwargs)
        finally:
            for artist in animated: artist.set_animated(True)

def energies(settings):
    count = int(round((settings['eb_max'] - settings['eb_min']) / settings['de'])) + 1
    return [round(settings['eb_min'] + i * settings['de'], 8) for i in range(count)]


def cartesian_axis(settings):
    return np.linspace(settings['k_min'], settings['k_max'],
                       int(round((settings['k_max'] - settings['k_min']) / settings['k_step'])) + 1)


def sphere_coordinates(settings, v0, eb):
    axis = cartesian_axis(settings)
    x, y = np.meshgrid(axis, axis)
    z2 = (settings['hv'] + v0 - eb) / H2M - x*x - y*y
    cart = np.stack((x, y, np.sqrt(np.maximum(z2, 0))), axis=-1).reshape(-1, 3)
    return cart, z2 >= 0


def pad_request(grid, s, v0, eb, field, binary=SPADEXP):
    request = {'version': 2, 'spectrum': 'unbroadened-bin-density', 'grid': core.stamp(grid.path), 'eb': eb, 'hv': s['hv'], 'v0': v0,
               'wf': s['wfnc'], 'de': s['de'], 'field': field_request(field),
               'k': [s['k_min'], s['k_max'], s['k_step']], 'spadexp': core.stamp(binary),
               'AO': core.stamp(ATOMIC_ORBITALS),
               'potential': core.stamp(ATOMIC_POTENTIALS)}
    return request


def pad_volume_request(grid, settings, v0, requested, components, binary):
    """Describe one complete unbroadened PAD volume; dK and dE stay display-only."""
    component_requests = []
    for polarization in components:
        component_settings = dict(settings, polarization=polarization)
        field, _ = light_vectors(component_settings)
        request = pad_request(grid, component_settings, v0, 0.0, field, binary)
        request.pop('eb')
        component_requests.append({'polarization': polarization, 'request': request})
    return {
        'version': 1, 'layout': 'EB-ky-kx',
        'binding_energies': [float(eb) for eb in requested],
        'components': component_requests,
    }


def load_pad_volume(request, settings):
    """Load the single latest PAD volume when every calculation input matches."""
    expected = json.dumps(request, sort_keys=True)
    requested = np.asarray(request['binding_energies'], dtype=float)
    count = len(cartesian_axis(settings)); shape = (len(requested), count, count)
    try:
        with np.load(PAD_VOLUME, allow_pickle=False) as saved:
            if (str(saved['request']) != expected
                    or not np.array_equal(saved['binding_energies'], requested)
                    or saved['pad'].shape != shape
                    or saved['band_weight'].shape != shape
                    or saved['error'].shape != (len(requested),)):
                return None
            pad = saved['pad'].copy(); weights = saved['band_weight'].copy()
            errors = saved['error'].copy()
    except (OSError, ValueError, KeyError, EOFError):
        return None
    extent = (settings['k_min'], settings['k_max']) * 2
    return {float(eb): {'eb': float(eb), 'pad': pad[i], 'band_weight': weights[i],
                        'error': float(errors[i]), 'extent': extent}
            for i, eb in enumerate(requested)}


def save_pad_volume(request, raw_maps):
    """Atomically replace the latest cache with one compressed NumPy volume."""
    requested = np.asarray(request['binding_energies'], dtype=float)
    maps = [raw_maps[float(eb)] for eb in requested]
    pad = np.stack([item['pad'] for item in maps])
    weights = np.stack([item['band_weight'] for item in maps])
    errors = np.asarray([item['error'] for item in maps], dtype=float)
    tmp = Path(str(PAD_VOLUME) + '.tmp')
    with tmp.open('wb') as out:
        np.savez_compressed(out, binding_energies=requested, pad=pad,
                            band_weight=weights, error=errors,
                            request=json.dumps(request, sort_keys=True))
    tmp.replace(PAD_VOLUME)
    return PAD_VOLUME.stat().st_size


def structure_files(structure_path):
    """Return the CIF and OpenMX template texts of a calculation target."""
    info = {}
    if not structure_path:
        return info
    cif = structure_cif_path(structure_path)
    if cif is not None:
        info["cif"] = cif.read_text(errors="replace"); info["cif_name"] = cif.name
    template = Path(structure_path) / ".openmx.base.dat"
    if template.is_file():
        info["openmx_template"] = template.read_text(errors="replace")
    return info


def pad_result_complete(result):
    """True when every EB slice of the result's settings has been calculated."""
    if not result or "v0" not in result or "s" not in result:
        return False
    raw = result.get("raw_maps", {})
    return all(eb in raw for eb in energies(result["s"]))


def save_pad_result(path, result, structure_path=None):
    """Write a complete PAD result and everything needed to show it to one HDF5 file.

    The file holds the unbroadened PAD and band-weight volumes (EB, ky, kx),
    every GUI/OpenMX setting, the fixed-coordinate band cut, the crystal
    structure (CIF, OpenMX template, lattice) and calculation metadata.
    Gaussian dK/dE are stored as settings and reapplied when the file is read.
    """
    if not pad_result_complete(result):
        raise ValueError("The PAD result is incomplete. Finish the calculation before saving it.")
    s = validate_settings(result["s"]); ebs = energies(s); raw = result["raw_maps"]
    axis = cartesian_axis(s); count = len(axis)
    pad = np.stack([np.asarray(raw[eb]["pad"], dtype=float) for eb in ebs])
    weight = np.stack([np.asarray(raw[eb]["band_weight"], dtype=float) for eb in ebs])
    error = np.asarray([float(raw[eb]["error"]) for eb in ebs])
    structure = structure_files(structure_path) or dict(result.get("structure") or {})
    lattice = np.asarray(result.get("lattice", core.LATTICE), dtype=float)
    reciprocal = np.asarray(result.get("reciprocal", core.G), dtype=float)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    volume = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=(1, count, count))
    text = h5py.string_dtype()
    with h5py.File(tmp, "w") as f:
        f.attrs.update(format=PAD_RESULT_FORMAT, format_version=PAD_RESULT_VERSION,
                       software=f"SPADExp4D {__version__}",
                       created=datetime.now().astimezone().isoformat(timespec="seconds"),
                       structure_name=result.get("structure_name") or "")
        settings = f.create_group("Settings")
        settings.attrs["json"] = json.dumps(s, sort_keys=True, ensure_ascii=False)
        for key, value in s.items():
            if isinstance(value, (bool, int, float, str)):
                settings.attrs[key] = value
        group = f.create_group("PAD")
        group.attrs.update(layout="EB, ky, kx", v0_eV=float(result["v0"]),
                           description="Unbroadened PAD and band weight. Gaussian dK/dE "
                                       "(Settings) are applied when the file is shown.")
        group.create_dataset("binding_energy", data=np.asarray(ebs)).attrs["units"] = "eV"
        group.create_dataset("k", data=axis).attrs["units"] = "1/Angstrom (kx and ky)"
        group.create_dataset("pad", data=pad, **volume)
        group.create_dataset("band_weight", data=weight, **volume)
        group.create_dataset("nearest_eigenstate_distance", data=error).attrs["units"] = "1/Angstrom"
        if result.get("pad_request"):
            group.attrs["request_json"] = json.dumps(result["pad_request"], sort_keys=True)
        line = result.get("line")
        if line:
            cut = f.create_group("BandCut")
            cut.attrs.update(direction=s["direction"], a=s["a"], b=s["b"],
                             error=float(line.get("error", 0.0)), method=str(line.get("method", "")),
                             path=str(line.get("path", "")))
            if line.get("k") is not None and line.get("bands") is not None:
                cut.create_dataset("k", data=np.asarray(line["k"], dtype=float)).attrs["units"] = "1/Angstrom"
                cut.create_dataset("energy", data=np.asarray(line["bands"], dtype=float)).attrs["units"] = "eV, E - EF"
            if line.get("points") is not None:
                cut.create_dataset("points", data=np.asarray(line["points"], dtype=float).reshape(-1, 2)
                                   ).attrs["units"] = "k (1/Angstrom), E - EF (eV)"
        crystal = f.create_group("Structure")
        crystal.create_dataset("lattice", data=lattice).attrs["units"] = "Angstrom (rows a1, a2, a3)"
        crystal.create_dataset("reciprocal", data=reciprocal).attrs["units"] = "1/Angstrom (rows G1, G2, G3)"
        for key in ("cif", "openmx_template"):
            if structure.get(key):
                crystal.create_dataset(key, data=structure[key], dtype=text)
        if structure.get("cif_name"):
            crystal.attrs["cif_name"] = structure["cif_name"]
        calculation = f.create_group("Calculation")
        calculation.attrs["sampling_error"] = float(result.get("sampling_error", error.max(initial=0)))
        for key in CALCULATION_KEYS:
            value = result.get(key)
            if isinstance(value, (int, float, np.integer, np.floating)):
                calculation.attrs[key] = value
            elif value is not None:
                calculation.attrs[key] = str(value)
    tmp.replace(path)
    return path


def load_pad_result(path):
    """Read a file written by save_pad_result() as an undisplayed result."""
    path = Path(path)
    try:
        with h5py.File(path, "r") as f:
            if f.attrs.get("format") != PAD_RESULT_FORMAT:
                raise ValueError(f"Not a SPADExp4D PAD result file: {path}")
            version = int(f.attrs.get("format_version", 0))
            if version > PAD_RESULT_VERSION:
                raise ValueError(f"This file was written by a newer SPADExp4D (format {version}).")
            s = validate_settings(json.loads(f["Settings"].attrs["json"]))
            group = f["PAD"]; ebs = energies(s); count = len(cartesian_axis(s))
            stored = group["binding_energy"][()]
            pad = group["pad"][()]; weight = group["band_weight"][()]
            error = group["nearest_eigenstate_distance"][()]
            shape = (len(ebs), count, count)
            if (len(stored) != len(ebs) or not np.allclose(stored, ebs, rtol=0, atol=1e-6)
                    or pad.shape != shape or weight.shape != shape or error.shape != (len(ebs),)):
                raise ValueError(f"The PAD arrays do not match the stored settings: {path}")
            extent = (s["k_min"], s["k_max"]) * 2
            raw = {eb: {"eb": eb, "pad": pad[i], "band_weight": weight[i],
                        "error": float(error[i]), "extent": extent} for i, eb in enumerate(ebs)}
            result = dict(s=s, raw_maps=raw, maps={}, cuts={}, v0=float(group.attrs["v0_eV"]),
                          structure_name=str(f.attrs.get("structure_name", "")) or None,
                          created=str(f.attrs.get("created", "")), source_file=str(path))
            if "request_json" in group.attrs:
                result["pad_request"] = json.loads(group.attrs["request_json"])
            if "BandCut" in f:
                cut = f["BandCut"]
                result["line"] = {"k": cut["k"][()] if "k" in cut else None,
                                  "bands": cut["energy"][()] if "energy" in cut else None,
                                  "points": cut["points"][()] if "points" in cut else np.zeros((0, 2)),
                                  "error": float(cut.attrs.get("error", 0.0)),
                                  "method": str(cut.attrs.get("method", "")),
                                  "path": str(cut.attrs.get("path", ""))}
            crystal = f["Structure"]
            result["lattice"] = crystal["lattice"][()]; result["reciprocal"] = crystal["reciprocal"][()]
            structure = {key: crystal[key].asstr()[()] for key in ("cif", "openmx_template") if key in crystal}
            if "cif_name" in crystal.attrs:
                structure["cif_name"] = str(crystal.attrs["cif_name"])
            result["structure"] = structure
            calculation = f["Calculation"].attrs
            for key in CALCULATION_KEYS:
                if key in calculation:
                    value = calculation[key]
                    result[key] = value.item() if isinstance(value, np.generic) else value
            result["sampling_error"] = float(calculation.get("sampling_error", error.max(initial=0)))
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read the PAD result file {path}: {exc}") from exc
    return result


def prune_old_pad_artifacts():
    """Remove legacy slice caches and diagnostics from older PAD requests."""
    removed = 0; removed_bytes = 0
    patterns = ('slice_*.npz', 'slice_*.npz.tmp', 'batch_*_PAD.log',
                'batch_*_stdout.log', 'slice_*_PAD.log', 'slice_*_stdout.log')
    for pattern in patterns:
        for path in OUTPUT.glob(pattern):
            try: removed_bytes += path.stat().st_size
            except OSError: pass
            path.unlink(missing_ok=True); removed += 1
    (OUTPUT / 'latest_pad.json').unlink(missing_ok=True)
    return removed, removed_bytes


def ensure_cpp_binary(log=None, stop=None):
    """Build only the private executable; shared SPADExp binaries are read-only."""
    BUILD.mkdir(parents=True, exist_ok=True)
    with core.file_lock(BUILD / 'cpp_build.lock'):
        core.run_command(['make', '-C', PROJECT / 'cpp', '-j2'],
                         BUILD / 'cpp_build.log', PROJECT, expected=SPADEXP_4D,
                         on_output=log, stop=stop)
    return SPADEXP_4D


def write_batch_input(path, grid, s, v0, requested):
    """Small request file links the first-BZ database without copying LCAO."""
    count = len(cartesian_axis(s))
    cart, _ = sphere_coordinates(s, v0, requested[0])
    reciprocal = getattr(grid, 'reciprocal', None)
    with h5py.File(path, 'w') as out:
        inp = out.create_group('Input')
        for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
            inp[name] = h5py.ExternalLink(str(grid.path.resolve()), '/Input/' + name)
        for name in ('Output', 'Mesh'):
            out[name] = h5py.ExternalLink(str(grid.path.resolve()), '/' + name)
        k = inp.create_group('Kpath')
        k.attrs.update(Dimension=np.int32(2), Curved=np.uint8(1),
            Origin=frac((0, 0, math.sqrt((s['hv']+v0-requested[0])/H2M)), reciprocal),
            Xvector=frac((1, 0, 0), reciprocal), Yvector=frac((0, 1, 0), reciprocal),
            Xcount=count, Ycount=count, Count=count*count,
            Xrange=[s['k_min'], s['k_max']], Yrange=[s['k_min'], s['k_max']],
            minN=getattr(grid, 'band_range', core.BANDS)[0], maxN=getattr(grid, 'band_range', core.BANDS)[1])
        k.create_dataset('Coordinates', data=frac(cart, reciprocal))
        group = out.create_group('Batch')
        group.attrs.update(KMin=s['k_min'], KMax=s['k_max'], H2M=H2M, BandHartree_eV=EH, RawSpectrum=1, EnergyStep_eV=s['de'])
        group.create_dataset('BindingEnergies', data=requested)


def pad_batch(grid, s, v0, requested, field, receive, log=None, stop=None, binary=None):
    """Calculate every requested EB surface in one C++ process."""
    core.check_stop(stop)
    binary = binary or ensure_cpp_binary(log, stop)
    count = len(cartesian_axis(s))
    jobs = [(float(eb), pad_request(grid, s, v0, float(eb), field, binary))
            for eb in requested]
    grid_code = core.request_code(core.stamp(grid.path))
    with core.file_lock(OUTPUT / f'cpp_grid_{grid_code}.lock'):
        with tempfile.TemporaryDirectory(prefix='batch_', dir=INPUT) as temporary:
            directory = Path(temporary); source = directory / 'request.hdf5'
            output = directory / 'batch.hdf5'; inp = directory / 'pad.in'
            todo = [item[0] for item in jobs]
            write_batch_input(source, grid, s, v0, todo)
            code = core.request_code([item[1] for item in jobs])
            prepare_engine_directory(directory)
            inp.write_text(pad_text(source.name, output.name, ENGINE_LOG,
                s['hv'], v0, s['wfnc'], min(todo), max(todo), s['de'], s['de']/10, field))
            consumed = 0
            def collect():
                nonlocal consumed
                while consumed < len(jobs):
                    path = Path(f'{output}.{consumed:06d}.h5')
                    if not path.exists(): break
                    eb, _ = jobs[consumed]
                    with h5py.File(path, 'r') as f:
                        if not f.attrs.get('Complete') or not f.attrs.get('RawSpectrum') or float(f.attrs['EB']) != eb:
                            raise RuntimeError('C++ PAD output conditions do not match.')
                        pad = f['PAD'][:]; weights = f['BandWeight'][:]
                        error = float(f.attrs['MaxSamplingDistance_Ainv'])
                    if pad.shape != (count, count) or weights.shape != pad.shape:
                        raise RuntimeError('C++ PAD output grids do not match.')
                    _, valid = sphere_coordinates(s, v0, eb)
                    pad = np.where(valid, pad, 0); weights = np.where(valid, weights, 0)
                    consumed += 1
                    receive(dict(eb=eb, pad=pad, band_weight=weights, error=error,
                                 extent=(s['k_min'], s['k_max'])*2))
                    path.unlink()
            def on_output(line):
                if log: log(line)
                collect()  # Scan all closed slices, even if the log reader skipped lines.
            if log: log(f'[C++ PAD] Calculating {len(todo)} unsmoothed slices in one batch (Gaussian dK/dE applied during display).')
            try:
                elapsed = core.run_command([binary], OUTPUT / f'batch_{code}_stdout.log',
                    directory, stdin=inp, expected=output, on_output=on_output, stop=stop)
            finally:
                keep_engine_log(directory, OUTPUT / f'batch_{code}_PAD.log')
                collect()
            with h5py.File(output, 'r') as f:
                if not f.attrs.get('Complete') or consumed != len(todo):
                    raise RuntimeError('C++ PAD did not complete all EB calculations.')
                kernel = float(f.attrs['KernelSeconds'])
            if log: log(f'[C++ PAD] Total {elapsed:.3f}s / kernel {kernel:.3f}s / {len(todo)} slices')
            return dict(pad_seconds=elapsed, pad_kernel_seconds=kernel, pad_computed_slices=len(todo))


def pad_slice(grid, settings, v0, eb, field, log=None, stop=None):
    """Calculate one EB surface for the Python compatibility engine."""
    core.check_stop(stop)
    s = settings
    sigma = s['de']/10  # Legacy engine tail_index=0; remove its peak normalization below.
    # A tiny positive pixel offset gives the legacy round() the raw half-open
    # bin convention at exact half-bin boundaries, without changing k geometry.
    engine_eb = eb + s['de']*1e-9
    request = pad_request(grid, s, v0, eb, field)
    code = core.request_code(request)
    count = len(cartesian_axis(s))
    cart, valid = sphere_coordinates(s, v0, eb)
    # Origin and Coordinates must share the cell the engine uses (grid file).
    reciprocal = getattr(grid, 'reciprocal', None)
    attrs = dict(Dimension=np.int32(2), Curved=np.uint8(1),
                 Origin=frac((0, 0, math.sqrt(max(s['hv']+v0-eb, 0)/H2M)), reciprocal),
                 Xvector=frac((1, 0, 0), reciprocal), Yvector=frac((0, 1, 0), reciprocal),
                 Xcount=np.int32(count), Ycount=np.int32(count),
                 Xrange=np.array((s['k_min'], s['k_max'])), Yrange=np.array((s['k_min'], s['k_max'])))
    with tempfile.TemporaryDirectory(prefix='pad_', dir=INPUT) as temporary:
        directory = Path(temporary); source = directory / 'states.hdf5'; output = directory / 'pad.hdf5'
        error = grid.write_path(source, cart, attrs)
        inp = directory / 'pad.in'
        prepare_engine_directory(directory)
        inp.write_text(pad_text(source.name, output.name, ENGINE_LOG,
                               s['hv'], v0, s['wfnc'], engine_eb, engine_eb, s['de'], sigma, field))
        with h5py.File(source, 'r') as f:
            bands = np.concatenate([f[path][:] for path in core.band_paths(core.spin_of(f))], axis=1)
            bands = (bands - grid.ef) * PAD_EH
        # Raw densities: each eigenstate occupies one half-open EB bin.
        active = np.floor((bands + eb) / s['de'] + 0.5) == 0
        if not active.any():
            pad = np.zeros((count, count))
            if log: log(f'[PAD skipped] EB={eb:g}: intensity set to 0 because the EB bin contains no states')
        else:
            try:
                core.run_command([SPADEXP], OUTPUT / f'slice_{code}_stdout.log', directory,
                                 stdin=inp, expected=output, on_output=log, stop=stop)
            finally:
                keep_engine_log(directory, OUTPUT / f'slice_{code}_PAD.log')
            if not pad_ok(output, s['hv'], v0, s['wfnc'], 2, engine_eb, engine_eb, s['de'], sigma, field):
                raise RuntimeError('PAD output energy or calculation conditions do not match.')
            with h5py.File(output, 'r') as f:
                if f['Dispersion'].shape != (count, count, 1):
                    raise RuntimeError('PAD output grid point counts do not match.')
                pad = f['Dispersion'][:, :, 0].T * (sigma*math.sqrt(2*math.pi)/s['de'])
        weights = active.sum(axis=1).reshape(count, count) / s['de']
        pad = np.where(valid, pad, 0); weights = np.where(valid, weights, 0)
    return {'eb': eb, 'pad': pad, 'band_weight': weights, 'error': error,
            'extent': (s['k_min'], s['k_max'])*2}


def interpolate_square(values, axis, x, y):
    x, y = np.broadcast_arrays(np.atleast_1d(x), np.atleast_1d(y))
    outside = (x < axis[0]) | (x > axis[-1]) | (y < axis[0]) | (y > axis[-1])
    u = np.clip((x-axis[0])/(axis[-1]-axis[0])*(len(axis)-1), 0, len(axis)-1)
    v = np.clip((y-axis[0])/(axis[-1]-axis[0])*(len(axis)-1), 0, len(axis)-1)
    i = np.floor(u).astype(int); j = np.floor(v).astype(int)
    ii = np.minimum(i+1, len(axis)-1); jj = np.minimum(j+1, len(axis)-1)
    tx = u-i; ty = v-j
    result = (1-tx)*(1-ty)*values[j,i] + tx*(1-ty)*values[j,ii] + (1-tx)*ty*values[jj,i] + tx*ty*values[jj,ii]
    return np.where(outside, np.nan, result)


def cut_slice(data, s, v0):
    axis = cartesian_axis(s); eb = data['eb']
    if s['direction'] == 'kx': x, y = axis, np.full_like(axis, s['a'])
    elif s['direction'] == 'ky': x, y = np.full_like(axis, s['a']), axis
    else: x, y = np.array([s['a']]), np.array([s['b']])
    kp2 = x*x + y*y; z2 = (s['hv']+v0-eb)/H2M - kp2
    k = axis if s['direction'] != 'kz' else np.sqrt(np.maximum(z2, 0))
    result = {'eb': eb, 'k': k}
    for name in ('pad', 'band_weight'):
        result[name] = interpolate_square(data[name], axis, x, y)
    vacuum = (s['hv'] - eb - s['wfnc']) / H2M
    visible = (vacuum > 0) & (kp2 <= vacuum * math.sin(math.radians(s['theta_max']))**2) & (z2 >= 0)
    result['pad'] = np.where(visible, result['pad'], np.nan)
    result['band_weight'] = np.where(z2 >= 0, result['band_weight'], np.nan)
    return result


def slice_path(s):
    """In-plane (kx, ky) of the PAD slice; kz follows the photoemission sphere."""
    axis = cartesian_axis(s)
    if s['direction'] == 'kx': return np.asarray(axis), np.full_like(axis, s['a'])
    if s['direction'] == 'ky': return np.full_like(axis, s['a']), np.asarray(axis)
    return np.array([s['a']], dtype=float), np.array([s['b']], dtype=float)


def pad_cut_bands(grid, s, v0, rows_per_chunk=None):
    """Bands along the PAD slice, i.e. on the photoemission sphere.

    At binding energy EB the PAD slice samples k = (kx, ky, kz) with
    kz = sqrt((hν + V0 − EB)/H2M − kx² − ky²), so kz changes with EB.  Band n
    appears in the slice where EB = −E_n(k(EB)) (E relative to EF).  The sign
    changes of EB + E_n between neighbouring EB rows locate these points, which
    are refined linearly.  Eigenvalues are trilinearly interpolated on the mesh.

    Returns 'k' and 'bands' (E − EF, one column per band, NaN where the band
    is outside the slice) for lines, and 'points' (k, E − EF) where a band
    meets one k column more than once.  For a kz slice every point is in
    'points' with k = kz.
    """
    ebs = np.asarray(energies(s), dtype=float); x, y = slice_path(s)
    kp2 = x * x + y * y; columns = len(x)
    # About 4096 k points per interpolation call keeps memory to tens of MB.
    rows_per_chunk = rows_per_chunk or max(1, 4096 // columns)
    previous = previous_eb = None; error = 0.0; band_count = 0
    found_columns, found_bands, found_roots = [], [], []
    for start in range(0, len(ebs), rows_per_chunk):
        block = ebs[start:start + rows_per_chunk]
        z2 = (s['hv'] + v0 - block[:, None]) / H2M - kp2[None, :]
        rows, cols = np.nonzero(z2 >= 0)
        g = None
        if len(rows):
            cart = np.stack((x[cols], y[cols], np.sqrt(z2[rows, cols])), axis=1)
            energy, distance = grid.interpolated_bands(cart)
            energy = np.asarray(energy, dtype=float); band_count = energy.shape[1]
            error = max(error, float(np.max(distance, initial=0.0)))
            g = np.full((len(block), columns, band_count), np.nan)
            g[rows, cols] = block[rows, None] + energy  # EB + E_n: zero on the band
        for i, eb in enumerate(block):
            current = None if g is None else g[i]
            if previous is not None and current is not None:
                crossing = ((previous < 0) & (current >= 0)) | ((previous >= 0) & (current < 0))
                col, band = np.nonzero(crossing)
                p, c = previous[col, band], current[col, band]
                found_columns.append(col); found_bands.append(band)
                found_roots.append(previous_eb + (eb - previous_eb) * p / (p - c))
            previous, previous_eb = current, eb
    col = np.concatenate(found_columns) if found_columns else np.zeros(0, dtype=int)
    band = np.concatenate(found_bands) if found_bands else np.zeros(0, dtype=int)
    root = np.concatenate(found_roots) if found_roots else np.zeros(0)
    counts = np.zeros((columns, band_count), dtype=int)
    np.add.at(counts, (col, band), 1)
    single = counts[col, band] == 1 if s['direction'] != 'kz' else np.zeros(len(col), dtype=bool)
    lines = np.full((columns, band_count), np.nan)
    lines[col[single], band[single]] = -root[single]
    lines = lines[:, np.isfinite(lines).any(axis=0)]
    rest = ~single
    if s['direction'] == 'kz':
        k_points = np.sqrt(np.maximum((s['hv'] + v0 - root[rest]) / H2M - kp2[col[rest]], 0.0))
    else:
        k_points = cartesian_axis(s)[col[rest]]
    return {'k': cartesian_axis(s) if s['direction'] != 'kz' else None,
            'bands': lines if s['direction'] != 'kz' else None,
            'points': np.column_stack((k_points, -root[rest])),
            'error': error, 'path': 'photoemission sphere along the PAD slice',
            'method': 'EB = -E_n(k(EB)) crossings; periodic trilinear interpolation of eigenvalues'}


def prepare_display(result, settings, stop=None):
    """Broaden the cached (EB, ky, kx) densities once, then extract all views.

    Values outside the cached cube are zero. No matrix elements or eigenstates
    are recomputed; raw arrays remain immutable for later width changes.
    """
    started = time.monotonic(); raw = dict(result.get('raw_maps', {}))
    view = dict(result, s=dict(settings), raw_maps=raw, maps={}, cuts={})
    if not raw: return view
    ebs = energies(settings); count = len(cartesian_axis(settings))
    positions = {eb: i for i, eb in enumerate(ebs)}
    maps = {eb: dict(data) for eb, data in raw.items() if eb in positions}
    for name in ('pad', 'band_weight'):
        cube = np.zeros((len(ebs), count, count))
        for eb, data in maps.items(): cube[positions[eb]] = data[name]
        for axis, width, step in ((0, settings['dE'], settings['de']),
                                  (1, settings['dk'], settings['k_step']),
                                  (2, settings['dk'], settings['k_step'])):
            cube = gaussian_smooth_axis(cube, width, step, axis, mode='constant', stop=stop)
        for eb in maps: maps[eb][name] = cube[positions[eb]]
    view['maps'] = maps
    view['cuts'] = {eb: cut_slice(data, settings, result['v0']) for eb, data in maps.items()}
    view['smoothing_seconds'] = time.monotonic() - started
    return view


def can_reuse_display(result, settings):
    if not result or set(result.get('raw_maps', {})) != set(energies(settings)): return False
    display_keys = {'dk', 'dE', 'ef_shift', 'theta_max', 'show_bz'}
    return ({k: v for k, v in result['s'].items() if k not in display_keys}
            == {k: v for k, v in settings.items() if k not in display_keys})


class Worker(QtCore.QObject):
    done = QtCore.pyqtSignal(object); error = QtCore.pyqtSignal(str); stopped = QtCore.pyqtSignal(object)
    partial = QtCore.pyqtSignal(object); log = QtCore.pyqtSignal(str)
    progress = QtCore.pyqtSignal(int); status = QtCore.pyqtSignal(str)

    def __init__(self, settings, structure_path=None, save_last=False):
        super().__init__(); self.s = settings; self.structure_path = Path(structure_path) if structure_path else None
        self.structure_name = self.structure_path.name if self.structure_path else "graphite"
        self.stop_event = threading.Event()
        # Write the finished result to <target>/cache/LAST_PAD_RESULT so the
        # next start can show it (the application window enables this).
        self.save_last = save_last

    def request_stop(self):
        self.stop_event.set()

    def stage(self, percent, text):
        self.progress.emit(percent); self.status.emit(text); self.log.emit(text)

    def snapshot(self, result):
        self.partial.emit(prepare_display(result, self.s, self.stop_event.is_set))

    @QtCore.pyqtSlot()
    def run(self):
        result = dict(s=self.s, structure_name=self.structure_name, raw_maps={}, maps={}, cuts={},
                      openmx_seconds=0.0, conversion_seconds=0.0, openmx_hdf5_seconds=0.0)
        grid = None
        try:
            configure_runtime_paths(self.structure_path)
            core.configure_structure(self.structure_path)
            overrides = core.openmx_overrides(self.s); full = core.request_band_range(overrides)
            spin = core.effective_openmx(overrides)['scf.SpinPolarization']; window = core.band_window(self.s)
            self.log.emit(f'[OpenMX bands] scf.SpinPolarization {spin}: {full[1] - full[0] + 1} bands'
                          + (' per spin' if spin == 'on' else '')
                          + ('; all stored' if window is None
                             else f'; stored E−EF {window[0]:g} to {window[1]:g} eV'))
            mkdirs(); core.check_stop(self.stop_event.is_set)
            v0 = (automatic_v0(SCF_H5, self.s['v0'])
                  if self.structure_path is None else (self.s['v0'] if self.s['v0'] > 0 else 20.0))
            result['v0'] = v0
            if self.s['hv'] + v0 - self.s['eb_max'] <= 0:
                raise ValueError('The kinetic energy inside the material is not positive at EB End.')
            self.stage(0, 'Preparing all first-BZ k points (PAD conditions do not affect OpenMX recalculation)')
            mesh_started = time.monotonic()
            def conversion_progress(base, batch, current, total, all_total):
                fraction = base + batch * (current / max(total, 1))
                self.stage(min(49, int(50 * fraction / max(all_total, 1))),
                           f'Converting OpenMX output to HDF5: {current:,}/{total:,} k points '
                           f'(overall {int(fraction):,}/{all_total:,})')
            openmx_total = 0.0; conversion_total = 0.0
            def record_timing(openmx_seconds, conversion_seconds):
                nonlocal openmx_total, conversion_total
                openmx_total += openmx_seconds; conversion_total += conversion_seconds
            path, elapsed = core.ensure_mesh(self.s, log=self.log.emit, stop=self.stop_event.is_set,
                progress=lambda done, total: self.stage(int(50*done/total), f'First BZ: {done:,} / {total:,} k points'),
                conversion=conversion_progress, timing=record_timing)
            result['openmx_seconds'] = openmx_total if openmx_total else elapsed
            result['conversion_seconds'] = conversion_total
            result['openmx_hdf5_seconds'] = result['openmx_seconds'] + conversion_total
            result['mesh_seconds'] = time.monotonic() - mesh_started
            result['mesh_path'] = str(path)
            result['mesh_bytes'] = path.stat().st_size
            self.log.emit(f"[K-space HDF5] {path}; {format_file_size(result['mesh_bytes'])}")
            core.check_stop(self.stop_event.is_set)
            grid = core.EigenstateGrid(path)
            self.log.emit(f'[HDF5] {path}; LCAO=float32 + gzip; all first-BZ points complete')
            result['line'] = pad_cut_bands(grid, self.s, v0)
            result['sampling_error'] = result['line']['error']
            self.snapshot(result)
            requested = energies(self.s); components = polarization_components(self.s)
            combined = {}
            if len(components) == 2:
                self.log.emit('[unpolarized] Fetching and adding c+ and c- PAD intensities separately')
            last_snapshot = time.monotonic()
            def receive_final(data):
                nonlocal last_snapshot
                eb = data['eb']; result['raw_maps'][eb] = data
                result['sampling_error'] = max(result['sampling_error'], data['error'])
                now = time.monotonic()
                if now - last_snapshot >= 0.75 or len(result['raw_maps']) == len(requested):
                    self.snapshot(result); last_snapshot = time.monotonic()
            pad_totals = dict(pad_seconds=0.0, pad_kernel_seconds=0.0, pad_computed_slices=0)
            pad_started = time.monotonic()
            pad_binary = (ensure_cpp_binary(self.log.emit, self.stop_event.is_set)
                          if self.s.get('pad_engine', 'cpp') == 'cpp' else SPADEXP)
            with core.file_lock(OUTPUT / 'latest_pad.lock'):
                removed, removed_bytes = prune_old_pad_artifacts()
                if removed:
                    self.log.emit(f'[PAD cleanup] Removed {removed:,} slices/old logs, '
                                  f'{format_file_size(removed_bytes)} deleted')
                volume_request = pad_volume_request(grid, self.s, v0, requested,
                                                    components, pad_binary)
                result['pad_request'] = volume_request
                cached = load_pad_volume(volume_request, self.s)
                if cached is not None:
                    self.log.emit(f'[PAD volume reused] {PAD_VOLUME}; '
                                  f'{format_file_size(PAD_VOLUME.stat().st_size)} '
                                  '(Gaussian dK/dE reapplied)')
                    for index, eb in enumerate(requested, 1):
                        self.stage(50+int(49*index/len(requested)),
                                   f'PAD volume: {index}/{len(requested)} slices, EB={eb:g} eV')
                        receive_final(cached[float(eb)])
                else:
                    for component_index, polarization in enumerate(components):
                        component_settings = dict(self.s, polarization=polarization)
                        field, _ = light_vectors(component_settings)
                        component_done = [0]
                        def receive_component(data, index=component_index, name=polarization):
                            component_done[0] += 1
                            overall = index * len(requested) + component_done[0]
                            self.stage(50+int(49*overall/(len(components)*len(requested))),
                                       f'PAD {name}: {component_done[0]}/{len(requested)} slices, EB={data["eb"]:g} eV')
                            if len(components) == 1:
                                receive_final(data)
                            elif index == 0:
                                combined[data['eb']] = dict(data)
                            else:
                                receive_final(sum_pad_components(combined.pop(data['eb']), data))
                        if self.s.get('pad_engine', 'cpp') == 'cpp':
                            stats = pad_batch(grid, component_settings, v0, requested, field,
                                              receive_component, self.log.emit, self.stop_event.is_set,
                                              binary=pad_binary)
                            for key in pad_totals:
                                pad_totals[key] += stats.get(key, 0.0)
                        else:
                            for eb in requested:
                                core.check_stop(self.stop_event.is_set)
                                receive_component(pad_slice(grid, component_settings, v0, eb, field,
                                                            self.log.emit, self.stop_event.is_set))
                    core.check_stop(self.stop_event.is_set)
                    volume_bytes = save_pad_volume(volume_request, result['raw_maps'])
                    self.log.emit(f'[PAD volume saved] {PAD_VOLUME}; '
                                  f'{format_file_size(volume_bytes)} / {len(requested)} EB slices')
            if self.s.get('pad_engine', 'cpp') == 'cpp':
                result.update(pad_totals)
            else:
                result['pad_seconds'] = time.monotonic() - pad_started
            self.stage(99, 'Applying Gaussian dK/dE to the saved PAD')
            result = prepare_display(result, self.s, self.stop_event.is_set)
            if self.save_last and self.structure_path is not None:
                self.stage(99, 'Saving the PAD result for the next start')
                try:
                    saved = save_pad_result(self.structure_path / 'cache' / LAST_PAD_RESULT,
                                            result, self.structure_path)
                    self.log.emit(f'[PAD result saved] {saved}; {format_file_size(saved.stat().st_size)}')
                except Exception as exc:  # a failed save must not fail the calculation
                    self.log.emit(f'[PAD result not saved] {exc}')
            self.stage(100, 'Complete'); self.done.emit(result)
        except core.Cancelled:
            self.stopped.emit(prepare_display(result, self.s))
        except Exception:
            self.error.emit(traceback.format_exc())
        finally:
            if grid is not None: grid.close()


class DisplayRelay(QtCore.QObject):
    """Delivers a result smoothed in a background thread to the GUI thread."""
    ready = QtCore.pyqtSignal(int, object, object)  # job, display result, error text


class Window(QtWidgets.QMainWindow):
    def __init__(self, initial_settings, session=False):
        super().__init__(); self.setWindowTitle("SPADExp4D — first-BZ eigenstate cache"); self.resize(1500,1050)
        self.initial_settings = validate_settings(initial_settings)
        self.startup_settings = self.initial_settings
        # session=True (the application, not tests): restore the last window
        # state and PAD result, save them again on exit, and keep every
        # finished result in <target>/cache/LAST_PAD_RESULT.
        self.persist_session = session
        self._session = read_session() if session else {}
        self._preferred_structure = self._session.get("structure")
        self.result = self.thread = self.worker = None; self._closing = False
        self._display_job = 0; self._pending_display = None
        self._display_relay = DisplayRelay(self); self._display_relay.ready.connect(self.background_display_ready)
        self.eb_sweep_timer = QtCore.QTimer(self)
        self.eb_sweep_timer.setInterval(250)
        self.eb_sweep_timer.timeout.connect(self.advance_eb_sweep)
        self.current_structure = None; self.build(); self.direction_changed()
        if session:
            self.restore_session()
    def spin(self, value, lo=-100, hi=100):
        x=QtWidgets.QDoubleSpinBox(); x.setRange(lo,hi); x.setDecimals(3); x.setValue(value); x.setSingleStep(.1); return x
    def build(self):
        initial = self.initial_settings
        ensure_app_folders()
        self.build_menus()
        w=QtWidgets.QWidget(); self.setCentralWidget(w)
        outer=QtWidgets.QHBoxLayout(w); outer.setContentsMargins(4,4,4,4); outer.setSpacing(4)
        sidebar=QtWidgets.QWidget(); sidebar.setObjectName("parameterSidebar")
        sidebar.setMinimumWidth(500); sidebar.setMaximumWidth(600)
        side=QtWidgets.QVBoxLayout(sidebar); side.setContentsMargins(6,6,6,6); side.setSpacing(1)
        title=QtWidgets.QLabel("Calculation and display parameters"); title.setStyleSheet("font-weight: bold; padding-bottom: 2pt;")
        side.addWidget(title)
        g=QtWidgets.QGridLayout(); g.setContentsMargins(0,0,0,0); g.setHorizontalSpacing(4); g.setVerticalSpacing(1)
        row=0
        def row_label(label):
            # A label may also be a widget such as a button (incident angles).
            lab=(label if isinstance(label,QtWidgets.QWidget) else QtWidgets.QLabel(label))
            if isinstance(lab,QtWidgets.QLabel): lab.setStyleSheet("margin: 0; padding: 0;")
            return lab
        def add_row(label, widget):
            nonlocal row
            g.addWidget(row_label(label),row,0); g.addWidget(widget,row,1); row+=1
        def form_page(rows):
            """A tab page laid out like the sidebar: one label and value per row."""
            page=QtWidgets.QWidget(); grid=QtWidgets.QGridLayout(page)
            grid.setContentsMargins(6,6,6,6); grid.setHorizontalSpacing(4); grid.setVerticalSpacing(1)
            for line,(label,widget) in enumerate(rows):
                grid.addWidget(row_label(label),line,0); grid.addWidget(widget,line,1)
            grid.setColumnStretch(1,1); grid.setRowStretch(len(rows),1)
            return page
        def add_full_row(widget):
            nonlocal row
            g.addWidget(widget,row,0,1,2); row+=1
        def inline(*widgets):
            holder=QtWidgets.QWidget(); layout=QtWidgets.QHBoxLayout(holder)
            layout.setContentsMargins(0,0,0,0); layout.setSpacing(3)
            for widget in widgets: layout.addWidget(widget,1)
            return holder
        # Only the tail of the path is shown, so a long data folder path never
        # widens the sidebar; the full path is in the tooltip.
        self.data_folder_button = QtWidgets.QPushButton()
        self.data_folder_button.setStyleSheet("text-align: left; padding: 2px 6px;")
        self.data_folder_button.clicked.connect(self.show_preferences)
        add_row("Data folder", self.data_folder_button)
        self.update_data_folder_button()
        self.structure_combo = QtWidgets.QComboBox()
        self.structure_combo.setToolTip("Calculation targets under <data folder>/<structure_name> (Settings > Preferences…)")
        self.structure_combo.currentIndexChanged.connect(self.structure_changed)
        add_row("Calculation target", self.structure_combo)
        self.result_combo = QtWidgets.QComboBox()
        self.result_combo.setToolTip("Saved calculation results for the selected structure")
        self.result_combo.currentIndexChanged.connect(self.result_selection_changed)
        add_row("Calculation result", self.result_combo)
        # OpenMX settings are grouped in tabs: the first-BZ k division, and the
        # OpenMX input keywords (filled in as they become editable).
        self.openmx_tabs=QtWidgets.QTabWidget()
        division=QtWidgets.QWidget(); kg=QtWidgets.QGridLayout(division)
        kg.setContentsMargins(6,6,6,6); kg.setHorizontalSpacing(12); kg.setVerticalSpacing(1)
        self.nboxes=[]; self.glabels=[]
        for i in range(3):
            box=QtWidgets.QSpinBox(); box.setRange(1,256); box.setValue(initial[f"N{i+1}"])
            box.setToolTip("First-BZ divisions. Changing these selects or calculates the corresponding eigenstate cache.")
            self.nboxes.append(box)
            box.valueChanged.connect(self.update_mesh_labels)
        kg.addWidget(QtWidgets.QLabel("N1, N2, N3"),0,0); kg.addWidget(inline(*self.nboxes),0,1)
        for i in range(3):
            label=QtWidgets.QLabel(); label.setWordWrap(False); self.glabels.append(label)
            kg.addWidget(label,1+i,0,1,2)
        self.total_points=QtWidgets.QLabel(); self.total_points.setWordWrap(True)
        kg.addWidget(QtWidgets.QLabel("Total k-point count"),4,0); kg.addWidget(self.total_points,4,1)
        kg.setColumnStretch(1,1)
        self.openmx_tabs.addTab(division,"OpenMX K division")
        self.openmx_tabs.addTab(self.build_openmx_parameters(inline),"OpenMX Parameters")
        self.openmx_tabs.setTabToolTip(1,OPENMX_PARAMETERS_HELP)
        self.load_openmx_parameters(initial["openmx"]); self.set_band_window(initial["band_window"])
        add_full_row(self.openmx_tabs)
        self.mesh_result=QtWidgets.QLabel("Not available"); self.mesh_result.setWordWrap(True)
        self.mesh_result.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        add_row("K-space HDF5",self.mesh_result)
        self.update_mesh_labels()
        self.pad_engine=QtWidgets.QComboBox()
        self.pad_engine.addItem("C++ batch", "cpp"); self.pad_engine.addItem("Python per EB (comparison)", "python")
        self.pad_engine.setCurrentIndex(self.pad_engine.findData(initial["pad_engine"]))
        self.hv=self.spin(initial["hv"],.01,10000); self.v0=self.spin(initial["v0"],-1000,1000)
        self.ef_shift=self.spin(initial["ef_shift"],-1000,1000); self.wfnc=self.spin(initial["wfnc"],0,100)
        self.ef_shift.setToolTip(
            "Display-only shift added to the calculated binding energy. "
            "Moves the band plot and energy slice within the fixed display range.")
        self.eb_min=self.spin(initial["eb_min"],-1000,1000); self.eb_max=self.spin(initial["eb_max"],-1000,1000); self.de=self.spin(initial["de"],.001,100)
        self.eb_min.setToolTip("Negative EB includes eigenstates above the original EF in the PAD calculation.")
        self.direction=QtWidgets.QComboBox(); self.direction.addItems(("kx","ky","kz")); self.direction.setCurrentText(initial["direction"]); self.direction.currentTextChanged.connect(self.direction_changed)
        self.la,self.lb=QtWidgets.QLabel(),QtWidgets.QLabel(); self.a,self.b=self.spin(initial["a"],-20,20),self.spin(initial["b"],-20,20)
        self.polarization=QtWidgets.QComboBox(); self.polarization.addItems(POLARIZATIONS); self.polarization.setCurrentText(initial["polarization"])
        self.polarization.setToolTip("unpolarized adds the c+ and c- PAD intensities")
        self.alpha=self.spin(initial["alpha"],-360,360); self.beta=self.spin(initial["beta"],0,180); self.gamma=self.spin(initial["gamma"],-360,360)
        self.euler_button=QtWidgets.QPushButton("Incident Angles (zyz Euler, °)…")
        self.euler_button.setToolTip("Euler angles a, b, c (z-y-z) of the incident light. "
                                     "Press to choose a beamline preset; the values stay editable.")
        self.euler_button.clicked.connect(self.choose_euler_preset)
        self.euler_button.setMinimumWidth(self.euler_button.sizeHint().width())  # never clip the text
        self.dk=self.spin(initial["dk"],0,100); self.dE=self.spin(initial["dE"],0,100)
        # Gaussian widths affect the displayed PAD, but changing a spin box
        # must not launch repeated processing while the user is typing.  The
        # values are applied by Calc (or the explicit View menu action).
        for box in (self.dk,self.dE):
            box.setToolTip("Gaussian standard deviation applied to saved PAD. 0 disables smoothing; values outside the calculation range are treated as 0.")
        self.theta_max=self.spin(initial["theta_max"],.001,90)
        self.show_bz=QtWidgets.QCheckBox("Show Brillouin zone"); self.show_bz.setChecked(initial["show_bz"]); self.show_bz.toggled.connect(self.redraw)
        self.k_min=self.spin(initial["k_min"],-K_LIMIT,K_LIMIT); self.k_max=self.spin(initial["k_max"],-K_LIMIT,K_LIMIT)
        self.k_step=self.spin(initial["k_step"],.001,2*K_LIMIT)
        # PAD conditions (with the PAD engine) and the band-diagram slice are
        # grouped in tabs.
        self.parameter_tabs=QtWidgets.QTabWidget()
        self.parameter_tabs.addTab(form_page((
            ("PAD calculation",self.pad_engine),
            ("hν (eV)",self.hv),
            ("Polarization",self.polarization),
            (self.euler_button,inline(self.alpha,self.beta,self.gamma)),
            ("EF shift (eV)",self.ef_shift),
            ("V0 (eV; ≤0:auto)",self.v0),
            ("wfnc (eV)",self.wfnc),
            ("Slice direction",self.direction),
            ("KMin, KMax, Kstep (Å⁻¹)",inline(self.k_min,self.k_max,self.k_step)),
            ("EB Begin, End, Step (eV)",inline(self.eb_min,self.eb_max,self.de)),
            ("Gaussian dK (Å⁻¹)",self.dk),
            ("Gaussian dE (eV)",self.dE),
            ("Maximum angle θmax (°)",self.theta_max),
            ("Boundary",self.show_bz))),"PAD")
        self.parameter_tabs.addTab(form_page((
            (self.la,self.a),
            (self.lb,self.b))),"Band Diagram")
        add_full_row(self.parameter_tabs)
        parameters=QtWidgets.QWidget(); parameters.setLayout(g)
        scroll=QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(parameters)
        scroll.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff); side.addWidget(scroll,1)
        self.calc=QtWidgets.QPushButton("Calc"); self.calc.clicked.connect(self.calculate)
        self.calc.setToolTip("Calculate or load the cache for the current target and parameters")
        side.addWidget(self.calc)
        self._slider_min=initial["eb_min"]; self._slider_step=initial["de"]
        self.slider=QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal); self.slider.setRange(0,int(round((initial["eb_max"]-initial["eb_min"])/initial["de"]))); self.slider.valueChanged.connect(self.slider_changed)
        shift=initial["ef_shift"]
        self.eb_input=self.spin(initial["eb_min"]+shift,initial["eb_min"]+shift,initial["eb_max"]+shift)
        self.eb_input.setSingleStep(initial["de"])
        # Do not commit every intermediate text/valueChanged event.  An EB
        # selection is accepted only when the user presses Enter.
        self.eb_input.setKeyboardTracking(False)
        self.eb_input.lineEdit().returnPressed.connect(self.eb_input_changed)
        self.ef_shift.valueChanged.connect(self.ef_shift_changed)
        self.play_eb_button=QtWidgets.QPushButton("▶ Play")
        self.play_eb_button.setToolTip("Display EB from Begin to End every 250 ms")
        self.play_eb_button.clicked.connect(self.start_eb_sweep)
        self.stop_eb_button=QtWidgets.QPushButton("■ Stop")
        self.stop_eb_button.setToolTip("Stop continuous EB display at the current position")
        self.stop_eb_button.clicked.connect(self.stop_eb_sweep)
        self.stop_eb_button.setEnabled(False)
        eb_row=QtWidgets.QHBoxLayout(); eb_row.setContentsMargins(0,0,0,0); eb_row.setSpacing(4)
        eb_row.addWidget(QtWidgets.QLabel("Current EB′ (eV)")); eb_row.addWidget(self.eb_input,1)
        eb_row.addWidget(self.play_eb_button); eb_row.addWidget(self.stop_eb_button)
        side.addLayout(eb_row); side.addWidget(self.slider)
        side.addWidget(QtWidgets.QLabel("Current operation"))
        self.progress=QtWidgets.QProgressBar(); self.progress.setRange(0,100); self.progress.setValue(0); self.progress.setFormat("Idle")
        self.message=QtWidgets.QLabel("Press Calc")
        self.message.setWordWrap(True); side.addWidget(self.progress); side.addWidget(self.message)
        self.logview=QtWidgets.QPlainTextEdit(); self.logview.setReadOnly(True); self.logview.setMinimumHeight(150); self.logview.setPlaceholderText("OpenMX / SPADExp output")
        self.logview.document().setMaximumBlockCount(1000); side.addWidget(self.logview)
        outer.addWidget(sidebar)
        graphics=QtWidgets.QWidget(); graph=QtWidgets.QVBoxLayout(graphics); graph.setContentsMargins(0,0,0,0); graph.setSpacing(1)
        self.canvas=Canvas(); graph.addWidget(NavigationToolbar2QT(self.canvas,self)); graph.addWidget(self.canvas,1)
        # Fast EB display: after each full draw keep the static background and
        # redraw only the EB-dependent artists (see update_eb_view).
        self._eb_artists=None; self._eb_background=None
        self.canvas.mpl_connect('draw_event',self.on_canvas_draw)
        self._eb_view_timer=QtCore.QTimer(self); self._eb_view_timer.setSingleShot(True)
        self._eb_view_timer.setInterval(0); self._eb_view_timer.timeout.connect(self.update_eb_view)
        outer.addWidget(graphics,1)
        self.refresh_structures()
    OPENMX_LABELS = {
        "scf.XcType": ("scf.XcType", "Exchange-correlation functional"),
        "scf.SpinPolarization": ("scf.SpinPolarization",
                                 "off: no spin; on: collinear up/down spins; "
                                 "NC: non-collinear spinors (needed for spin-orbit coupling)"),
        "scf.SpinOrbit.Coupling": ("SOC", "scf.SpinOrbit.Coupling: spin-orbit coupling (requires NC)"),
        "scf.energycutoff": ("scf.energycutoff (Ry)", "Real-space grid cutoff energy"),
        "scf.Kgrid": ("scf.Kgrid", "k-point grid of the SCF calculation (not the first-BZ mesh)"),
        "scf.ElectronicTemperature": ("scf.ElectronicTemperature (K)", "Electronic temperature for occupations"),
        "scf.criterion": ("scf.criterion (Hartree)", "SCF convergence criterion, e.g. 1e-7"),
        "scf.maxIter": ("scf.maxIter", "Maximum number of SCF iterations"),
    }

    def build_openmx_parameters(self, inline):
        """OpenMX keywords of the selected target, edited on top of its CIF template."""
        tab=QtWidgets.QWidget()
        form=self.openmx_parameter_form=QtWidgets.QFormLayout(tab)
        form.setContentsMargins(6,6,6,6); form.setVerticalSpacing(1)
        self.openmx_xc=QtWidgets.QComboBox(); self.openmx_xc.addItems(core.XC_TYPES)
        self.openmx_spin=QtWidgets.QComboBox(); self.openmx_spin.addItems(core.SPIN_MODES)
        self.openmx_soc=QtWidgets.QCheckBox("scf.SpinOrbit.Coupling")
        spin=QtWidgets.QWidget(); spin_row=QtWidgets.QHBoxLayout(spin)
        spin_row.setContentsMargins(0,0,0,0); spin_row.setSpacing(8)
        spin_row.addWidget(self.openmx_spin,1); spin_row.addWidget(self.openmx_soc)
        self.openmx_cutoff=QtWidgets.QDoubleSpinBox(); self.openmx_cutoff.setRange(10,10000)
        self.openmx_cutoff.setDecimals(1); self.openmx_cutoff.setSingleStep(10)
        self.openmx_kgrid=[QtWidgets.QSpinBox() for _ in range(3)]
        for box in self.openmx_kgrid: box.setRange(1,1000)
        self.openmx_temperature=QtWidgets.QDoubleSpinBox(); self.openmx_temperature.setRange(.1,100000)
        self.openmx_temperature.setDecimals(1); self.openmx_temperature.setSingleStep(100)
        self.openmx_criterion=QtWidgets.QLineEdit()
        self.openmx_criterion.setValidator(QtGui.QDoubleValidator(0.0,1.0e6,20))
        self.openmx_maxiter=QtWidgets.QSpinBox(); self.openmx_maxiter.setRange(1,100000)
        self.openmx_reset=QtWidgets.QPushButton("Reset to recommended")
        self.openmx_reset.setToolTip(OPENMX_RECOMMENDED_HELP)
        self.openmx_reset.clicked.connect(lambda: self.set_openmx_values(self.openmx_recommended))
        last=QtWidgets.QWidget(); last_row=QtWidgets.QHBoxLayout(last)
        last_row.setContentsMargins(0,0,0,0); last_row.setSpacing(6)
        last_row.addWidget(self.openmx_maxiter,1); last_row.addWidget(self.openmx_reset)
        editors={"scf.XcType":self.openmx_xc,"scf.SpinPolarization":spin,"scf.energycutoff":self.openmx_cutoff,
                 "scf.Kgrid":inline(*self.openmx_kgrid),"scf.ElectronicTemperature":self.openmx_temperature,
                 "scf.criterion":self.openmx_criterion,"scf.maxIter":last}
        self.openmx_labels={}; self.openmx_base={}; self.openmx_recommended={}
        for key in core.OPENMX_KEYWORDS:
            text,tip=self.OPENMX_LABELS[key]
            if key=="scf.SpinOrbit.Coupling":  # shares the spin row; its checkbox is the label
                self.openmx_labels[key]=self.openmx_soc; self.openmx_soc.setToolTip(tip); continue
            label=QtWidgets.QLabel(text); self.openmx_labels[key]=label
            label.setToolTip(tip); editors[key].setToolTip(tip)
            form.addRow(label,editors[key])
        # Stored band window (not an OpenMX keyword): limits the k-space HDF5.
        self.band_window_all=QtWidgets.QCheckBox("All")
        self.band_window_low=QtWidgets.QDoubleSpinBox(); self.band_window_high=QtWidgets.QDoubleSpinBox()
        for box,value in ((self.band_window_low,-25.0),(self.band_window_high,5.0)):
            box.setRange(-500,500); box.setDecimals(1); box.setSingleStep(1); box.setValue(value)
        window=QtWidgets.QWidget(); window_row=QtWidgets.QHBoxLayout(window)
        window_row.setContentsMargins(0,0,0,0); window_row.setSpacing(4)
        window_row.addWidget(self.band_window_all); window_row.addWidget(self.band_window_low,1)
        window_row.addWidget(QtWidgets.QLabel("to")); window_row.addWidget(self.band_window_high,1)
        window_label=QtWidgets.QLabel("Stored bands E−EF (eV)")
        for widget in (window_label,window): widget.setToolTip(BAND_WINDOW_HELP)
        form.addRow(window_label,window)
        self.band_window_all.toggled.connect(self.update_band_window_widgets)
        self.update_band_window_widgets()
        self.openmx_xc.currentTextChanged.connect(self.update_openmx_markers)
        self.openmx_spin.currentTextChanged.connect(self.spin_mode_changed)
        self.openmx_soc.toggled.connect(self.update_openmx_markers)
        for box in (self.openmx_cutoff,self.openmx_temperature,self.openmx_maxiter,*self.openmx_kgrid):
            box.valueChanged.connect(self.update_openmx_markers)
        self.openmx_criterion.textChanged.connect(self.update_openmx_markers)
        return tab

    @staticmethod
    def canonical_openmx(key, value):
        try:
            return core.normalize_openmx_value(key, value)
        except ValueError:
            return " ".join(str(value).split())

    def spin_mode_changed(self, *_):
        if self.openmx_spin.currentText()!="NC":
            self.openmx_soc.setChecked(False)  # OpenMX spin-orbit coupling needs NC
        self.openmx_soc.setEnabled(self.openmx_spin.currentText()=="NC")
        self.update_openmx_markers()

    def update_band_window_widgets(self, *_):
        for box in (self.band_window_low,self.band_window_high):
            box.setEnabled(not self.band_window_all.isChecked())

    def band_window_setting(self):
        if self.band_window_all.isChecked():
            return None
        return [self.band_window_low.value(),self.band_window_high.value()]

    def set_band_window(self, window):
        self.band_window_all.setChecked(window is None)
        if window is not None:
            self.band_window_low.setValue(float(window[0])); self.band_window_high.setValue(float(window[1]))
        self.update_band_window_widgets()

    def openmx_widget_values(self):
        return {"scf.XcType":self.openmx_xc.currentText(),
                "scf.SpinPolarization":self.openmx_spin.currentText(),
                "scf.SpinOrbit.Coupling":"on" if self.openmx_soc.isChecked() else "off",
                "scf.energycutoff":self.openmx_cutoff.value(),
                "scf.Kgrid":" ".join(str(box.value()) for box in self.openmx_kgrid),
                "scf.ElectronicTemperature":self.openmx_temperature.value(),
                "scf.criterion":self.openmx_criterion.text(),
                "scf.maxIter":self.openmx_maxiter.value()}

    def set_openmx_values(self, values):
        widgets=(self.openmx_xc,self.openmx_spin,self.openmx_soc,self.openmx_cutoff,
                 self.openmx_temperature,self.openmx_criterion,self.openmx_maxiter,*self.openmx_kgrid)
        for widget in widgets: widget.blockSignals(True)
        try:
            xc=values["scf.XcType"]
            if self.openmx_xc.findText(xc)<0: self.openmx_xc.addItem(xc)  # e.g. HSE from a template
            self.openmx_xc.setCurrentText(xc)
            self.openmx_spin.setCurrentText(values["scf.SpinPolarization"])
            self.openmx_soc.setChecked(values["scf.SpinOrbit.Coupling"]=="on")
            self.openmx_soc.setEnabled(self.openmx_spin.currentText()=="NC")
            self.openmx_cutoff.setValue(float(values["scf.energycutoff"]))
            for box,value in zip(self.openmx_kgrid,values["scf.Kgrid"].split()): box.setValue(int(value))
            self.openmx_temperature.setValue(float(values["scf.ElectronicTemperature"]))
            self.openmx_criterion.setText(values["scf.criterion"])
            self.openmx_maxiter.setValue(int(values["scf.maxIter"]))
        finally:
            for widget in widgets: widget.blockSignals(False)
        self.update_openmx_markers()

    def load_openmx_parameters(self, overrides):
        """Show the target's template values with its saved overrides applied."""
        template={}
        if self.current_structure is not None:
            try: template=core.template_keywords()
            except OSError: template={}
        self.openmx_base={key:self.canonical_openmx(key,template.get(key) or core.OPENMX_DEFAULTS[key])
                          for key in core.OPENMX_KEYWORDS}
        try:
            self.openmx_recommended=(core.recommended_openmx_values()
                                     if self.current_structure is not None else dict(self.openmx_base))
        except (OSError, ValueError):
            self.openmx_recommended=dict(self.openmx_base)
        self.set_openmx_values({**self.openmx_base,**(overrides or {})})
        self.openmx_tabs.widget(1).setEnabled(self.current_structure is not None)

    def openmx_overrides(self):
        """Keywords whose widget value differs from the template (invalid input is ignored)."""
        changed={}
        for key,value in self.openmx_widget_values().items():
            try: text=core.normalize_openmx_value(key,value)
            except ValueError: continue
            if text!=self.openmx_base.get(key): changed[key]=text
        return changed

    def update_openmx_markers(self, *_):
        changed=self.openmx_overrides()
        for key,label in self.openmx_labels.items():
            font=label.font(); font.setBold(key in changed); label.setFont(font)
            label.setToolTip(f"{self.OPENMX_LABELS[key][1]}\nTemplate: {self.openmx_base.get(key, '')}"
                             f"\nRecommended: {self.openmx_recommended.get(key, '')}")
        try:
            core.normalize_openmx_value("scf.criterion",self.openmx_criterion.text()); valid=True
        except ValueError:
            valid=False
        self.openmx_criterion.setStyleSheet("" if valid else "color: #c62828;")
        current={key:self.canonical_openmx(key,value) for key,value in self.openmx_widget_values().items()}
        self.openmx_reset.setEnabled(current!=self.openmx_recommended)

    def build_menus(self):
        file_menu = self.menuBar().addMenu("File")
        self.add_menu_action(file_menu, "Import CIF…", self.import_cif)
        self.add_menu_action(file_menu, "Open calculation target folder…", self.open_structure_folder)
        self.add_menu_action(file_menu, "Save parameter.json", self.save_selected_settings)
        file_menu.addSeparator()
        self.open_result_action = self.add_menu_action(file_menu, "Open PAD result (HDF5)…", self.import_pad_result)
        self.save_result_action = self.add_menu_action(file_menu, "Save PAD result as HDF5…", self.export_pad_result)
        file_menu.addSeparator(); self.add_menu_action(file_menu, "Exit", self.close)

        calc_menu = self.menuBar().addMenu("Calculation")
        self.calc_action = self.add_menu_action(calc_menu, "Calc", self.calculate)
        self.add_menu_action(calc_menu, "Show calculation log", self.show_log)
        self.cancel_action = self.add_menu_action(calc_menu, "Cancel calculation", self.stop_calculation)
        self.cancel_action.setEnabled(False)

        exchange_menu = self.menuBar().addMenu("Exchange")
        self.add_menu_action(exchange_menu, "Create upload job…", self.create_upload_job)
        self.add_menu_action(exchange_menu, "Import downloaded results…", self.import_download_results)
        self.add_menu_action(exchange_menu, "Open upload folder", lambda: self.open_path(UPLOAD))
        self.add_menu_action(exchange_menu, "Open download folder", lambda: self.open_path(DOWNLOAD))

        settings_menu = self.menuBar().addMenu("Settings")
        self.preferences_action = self.add_menu_action(settings_menu, "Preferences…", self.show_preferences)
        keys = QtGui.QKeySequence.keyBindings(QtGui.QKeySequence.StandardKey.Preferences)
        self.preferences_action.setShortcuts(keys or [QtGui.QKeySequence("Ctrl+,")])  # ⌘, on macOS
        # Keep the action in this menu; macOS would otherwise move any
        # "Preferences…" item into the application menu by text heuristics.
        self.preferences_action.setMenuRole(QtGui.QAction.MenuRole.NoRole)
        self.add_menu_action(settings_menu, "Open data folder", lambda: self.open_path(STRUCTURES))

        view_menu = self.menuBar().addMenu("View")
        self.add_menu_action(view_menu, "Reload calculation results", self.refresh_results)
        self.add_menu_action(view_menu, "Reapply Gaussian dK/dE", self.refresh_gaussian)
        self.bz_action = self.add_menu_action(view_menu, "Show Brillouin zone boundary", self.toggle_bz)
        self.bz_action.setCheckable(True); self.bz_action.setChecked(self.initial_settings["show_bz"])

        help_menu = self.menuBar().addMenu("Help")
        self.add_menu_action(help_menu, "About", self.show_about)

    @staticmethod
    def add_menu_action(menu, title, callback):
        action = QtGui.QAction(title, menu)
        action.triggered.connect(callback)
        menu.addAction(action)
        return action

    def open_path(self, path):
        ensure_app_folders(); path = Path(path); path.mkdir(parents=True, exist_ok=True)
        QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(path.resolve())))

    def refresh_structures(self):
        entries = structure_entries()
        selected = self.current_structure.name if self.current_structure else self._preferred_structure
        self.structure_combo.blockSignals(True); self.structure_combo.clear()
        if not entries: self.structure_combo.addItem(f"(Import a CIF into {STRUCTURES.name})", None)
        else:
            for path in entries: self.structure_combo.addItem(path.name, str(path))
            index = next((i for i, path in enumerate(entries) if path.name == selected), 0)
            self.structure_combo.setCurrentIndex(index)
        self.structure_combo.blockSignals(False)
        self.structure_changed(self.structure_combo.currentIndex())

    def save_structure_settings(self, structure_path=None):
        """Persist the PAD/view parameters owned by one structure."""
        structure_path = Path(structure_path or self.current_structure) if (structure_path or self.current_structure) else None
        if structure_path is None:
            return
        write_settings(structure_path / "parameter.json", self.settings())

    def structure_changed(self, index):
        global LATTICE, B
        self._stop_eb_sweep()
        value = self.structure_combo.itemData(index)
        next_structure = Path(value) if value else None
        previous_structure = self.current_structure
        if previous_structure is not None and previous_structure != next_structure:
            try:
                self.save_structure_settings(previous_structure)
            except (OSError, ValueError) as exc:
                self.append_log(f"Could not save PAD parameters: {previous_structure.name}: {exc}")
        self.cancel_background_display()
        self.current_structure = next_structure
        self.result = None
        preferred_result = None
        if self.current_structure:
            try:
                configure_runtime_paths(self.current_structure)
                core.configure_structure(self.current_structure)
                LATTICE = core.LATTICE; B = core.G
                settings_path = self.current_structure / "parameter.json"
                # A CIF-only folder is a valid user input.  Create the
                # companion settings file lazily so renaming the CIF is not
                # part of the normal workflow.
                if not settings_path.is_file():
                    # PAD settings carry over; OpenMX keyword overrides are
                    # material-specific, so a new target starts from its CIF.
                    write_settings(settings_path, dict(self.initial_settings, openmx={}, band_window=None))
                selected_settings = read_settings(settings_path)
                self.apply_settings(selected_settings)
                # Migrate a completed local first-BZ cache created before the
                # result catalog was implemented.  Only metadata and a
                # relative symlink are added; the HDF5 file is not copied.
                plan = core.mesh_plan(selected_settings)
                if plan["reuse"]:
                    preferred_result = register_local_result(
                        self.current_structure, selected_settings, plan["path"])
                self.message.setText(f"Calculation target: {self.current_structure.name}")
            except (OSError, ValueError) as exc:
                self.message.setText(f"Could not load settings: {exc}")
        else:
            configure_runtime_paths(DEFAULT_STRUCTURE)
            try:
                core.configure_structure(None)
                LATTICE = core.LATTICE; B = core.G
            except (OSError, ValueError, RuntimeError):
                # An empty data folder has no default graphite structure.
                # Keep the previous lattice until a CIF is imported.
                pass
            self.load_openmx_parameters({})  # disabled until a target exists
        self.refresh_results(preferred_result)

    def apply_settings(self, settings):
        settings = validate_settings(settings); self.initial_settings = settings
        for key, widget in (("hv", self.hv), ("v0", self.v0), ("ef_shift", self.ef_shift), ("wfnc", self.wfnc),
                            ("eb_min", self.eb_min), ("eb_max", self.eb_max), ("de", self.de),
                            ("a", self.a), ("b", self.b), ("dk", self.dk), ("dE", self.dE),
                            ("theta_max", self.theta_max), ("alpha", self.alpha),
                            ("beta", self.beta), ("gamma", self.gamma), ("k_min", self.k_min),
                            ("k_max", self.k_max), ("k_step", self.k_step)):
            widget.blockSignals(True); widget.setValue(settings[key]); widget.blockSignals(False)
        for i, box in enumerate(self.nboxes, 1):
            box.blockSignals(True); box.setValue(settings[f"N{i}"]); box.blockSignals(False)
        self.direction.setCurrentText(settings["direction"]); self.polarization.setCurrentText(settings["polarization"])
        self.pad_engine.setCurrentIndex(self.pad_engine.findData(settings["pad_engine"]))
        self.show_bz.setChecked(settings["show_bz"]); self.bz_action.setChecked(settings["show_bz"])
        self.load_openmx_parameters(settings["openmx"]); self.set_band_window(settings["band_window"])
        self.update_mesh_labels(); self.direction_changed(); self.update_eb_input()

    def refresh_results(self, selected=None):
        previous = self.result_combo.currentData() if self.result_combo.count() else None
        self.result_combo.blockSignals(True); self.result_combo.clear()
        # Keep the new-local-calculation state available even when previous
        # results exist for the selected structure.  ``None`` is the marker
        # consumed by calculate() to select the Local Machine path.
        self.result_combo.addItem("(No calculation results)", None)
        entries = []
        if self.current_structure:
            entries = result_entries(self.current_structure)
            for path in entries:
                self.result_combo.addItem(path.name, str(path))
        preferred = str(selected) if selected else previous
        preferred_index = self.result_combo.findData(preferred) if preferred else -1
        if preferred_index < 0 and self.current_structure:
            active = core.json_read(self.current_structure / "results" / "active.json")
            if active and active.get("path"):
                preferred = str(self.current_structure / "results" / active["path"])
                preferred_index = self.result_combo.findData(preferred)
        if preferred_index >= 0:
            self.result_combo.setCurrentIndex(preferred_index)
        self.result_combo.blockSignals(False)
        if self.result_combo.currentData() is not None and self.result is None:
            value = self.result_combo.currentData()
            # setCurrentIndex() happened while signals were blocked above,
            # so load the selected result explicitly.
            self.load_result_settings(value)

    def load_result_settings(self, value):
        """Use a saved result as the editable starting point for the next Calc."""
        try:
            result_path = Path(value)
            try:
                settings = read_settings(result_path / "config.json")
            except FileNotFoundError:
                # Imported legacy results may only contain eigenstates.hdf5;
                # retain the current display settings and replace its mesh
                # divisions from the cache below.
                settings = validate_settings(self.settings())
            # The HDF5 cache is authoritative for the first-BZ mesh.  It also
            # supports imported results that do not have a config file.
            eigenstates = result_path / "eigenstates.hdf5"
            if eigenstates.is_file():
                with h5py.File(eigenstates, "r") as file:
                    mesh = file.get("Mesh")
                    divisions = None if mesh is None else mesh.attrs.get("Divisions")
                if divisions is not None and len(divisions) == 3:
                    settings.update({f"N{i + 1}": int(value) for i, value in enumerate(divisions)})
                    settings = validate_settings(settings)
            self.apply_settings(settings)
            self.message.setText(f"Saved calculation result: {Path(value).name} (edit conditions and press Calc)")
        except (OSError, ValueError, KeyError) as exc:
            self.message.setText(f"Could not load saved settings: {exc}")

    def result_selection_changed(self, index):
        self._stop_eb_sweep()
        value = self.result_combo.itemData(index)
        self.result = None
        if value:
            self.load_result_settings(value)
        elif self.current_structure:
            self.message.setText(f"No calculation result: calculating {self.current_structure.name} locally")

    def import_cif(self):
        cif, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select crystal structure CIF", str(PROJECT), "CIF files (*.cif);;All files (*)")
        if not cif: return
        source = Path(cif)
        if source.suffix.lower() != ".cif":
            QtWidgets.QMessageBox.warning(self, "Import CIF", "Please select a CIF file."); return
        default_name = safe_name(source.stem)
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Import CIF", f"Folder name to create in the data folder ({STRUCTURES}):",
            QtWidgets.QLineEdit.EchoMode.Normal, default_name)
        if not ok or not name.strip(): return
        target = STRUCTURES / safe_name(name)
        if target.exists():
            QtWidgets.QMessageBox.warning(self, "Cannot create folder", f"Already exists: {target}"); return
        created = False
        try:
            # Make the current GUI state the newest structure parameter file,
            # then use that file as the imported structure's defaults.
            self.save_structure_settings()
            defaults_path = latest_parameter_path()
            defaults = read_settings(defaults_path)
            target.mkdir(parents=True); created = True
            destination = target / (source.stem + ".cif")
            shutil.copy2(source, destination)
            (target / "results" / "openmx").mkdir(parents=True)
            write_settings(target / "parameter.json", defaults)
            # Validate the CIF and prepare the reusable OpenMX base input now,
            # so an incomplete import is not added to the GUI.
            core.configure_structure(target)
        except (OSError, ValueError, RuntimeError) as exc:
            if created:
                shutil.rmtree(target, ignore_errors=True)
            QtWidgets.QMessageBox.critical(self, "Structure creation error", str(exc)); return
        self.refresh_structures()
        index = self.structure_combo.findData(str(target))
        if index >= 0:
            self.structure_combo.setCurrentIndex(index)
        self.message.setText(
            f"CIF imported: {target.name}\nInitial parameters: {defaults_path}")

    def create_structure(self):
        """Backward-compatible alias for the former menu action."""
        self.import_cif()

    def open_structure_folder(self):
        self.open_path(self.current_structure or STRUCTURES)

    def store_session(self):
        """Save the GUI parameters and window state for the next start."""
        if not self.persist_session:
            return
        try:
            settings = validate_settings(self.settings())
        except ValueError as exc:
            self.append_log(f"Could not save the window state: {exc}"); return
        pending = self._pending_display or {}
        state = {
            "structure": self.current_structure.name if self.current_structure else None,
            "eb": self.eb() if self.result else pending.get("eb"),
            # A result opened from a file is shown again instead of the last
            # calculated one.
            "result_file": (self.result or {}).get("source_file") or pending.get("source_file"),
            "parameter_tab": self.parameter_tabs.currentIndex(),
            "openmx_tab": self.openmx_tabs.currentIndex(),
            "geometry": bytes(self.saveGeometry().toBase64()).decode("ascii"),
        }
        try:
            save_session(settings, state)
        except (OSError, ValueError) as exc:
            self.append_log(f"Could not save the window state: {exc}")

    def restore_session(self):
        """Restore the state saved by store_session() and show the last PAD result."""
        state = self._session; self._preferred_structure = None
        name = self.current_structure.name if self.current_structure else None
        same_target = state.get("structure") == name
        if same_target and state:
            # Selecting the target applied its parameter.json and saved result
            # settings; the values shown at exit take precedence.
            self.apply_settings(self.startup_settings)
        for tabs, key in ((self.parameter_tabs, "parameter_tab"), (self.openmx_tabs, "openmx_tab")):
            index = state.get(key)
            if isinstance(index, int) and 0 <= index < tabs.count():
                tabs.setCurrentIndex(index)
        if isinstance(state.get("geometry"), str):
            self.restoreGeometry(QtCore.QByteArray.fromBase64(state["geometry"].encode("ascii")))
        if self.current_structure is None:
            return
        candidates = [self.current_structure / "cache" / LAST_PAD_RESULT]
        if same_target and state.get("result_file"):
            candidates.insert(0, Path(state["result_file"]))
        for path in candidates:
            if not path.is_file():
                continue
            try:
                result = load_pad_result(path)
            except ValueError as exc:
                self.append_log(f"Could not show the saved PAD result: {exc}"); continue
            if result.get("structure_name") not in (None, name):
                continue
            result["structure_name"] = name
            gui = validate_settings(self.settings())
            # Display-only edits made after the calculation (Gaussian widths,
            # EF shift, θmax, boundary) are kept; otherwise show the result as
            # it was calculated.
            display = gui if can_reuse_display(result, gui) else result["s"]
            self.append_log(f"[PAD result restored] {path}")
            self.display_in_background(
                result, display, state.get("eb") if same_target else None,
                f"Restored PAD result of {name} (calculated {result.get('created') or 'earlier'})")
            return

    def display_in_background(self, result, settings, eb=None, message=None):
        """Apply Gaussian dK/dE to a loaded result off the GUI thread, then show it.

        Smoothing a large volume takes seconds, so the window stays usable
        meanwhile.  Calc, another file or another target cancels the job.
        """
        self.cancel_background_display()
        self._display_job += 1; job = self._display_job
        self._pending_display = dict(job=job, target=self.current_structure, eb=eb, message=message,
                                     source_file=result.get("source_file"))
        self.progress.setRange(0, 0); self.progress.setFormat("Loading PAD result…")
        self.message.setText("Applying Gaussian dK/dE to the saved PAD result…")
        relay = self._display_relay
        def run():
            try:
                view, error = prepare_display(result, settings, lambda: job != self._display_job), None
            except core.Cancelled:
                return
            except Exception:
                view, error = None, traceback.format_exc()
            try:
                relay.ready.emit(job, view, error)
            except RuntimeError:  # the window was closed meanwhile
                pass
        threading.Thread(target=run, name="pad-result-display", daemon=True).start()

    def cancel_background_display(self):
        if self._pending_display is None:
            return
        self._display_job += 1; self._pending_display = None
        self.progress.setRange(0, 100); self.progress.setValue(0); self.progress.setFormat("Idle")

    def background_display_ready(self, job, view, error):
        pending = self._pending_display
        if pending is None or pending["job"] != job:
            return
        self._pending_display = None
        self.progress.setRange(0, 100)
        if self.thread is not None or self.current_structure != pending["target"]:
            return
        if error:
            self.progress.setValue(0); self.progress.setFormat("Could not show the PAD result")
            self.message.setText("Could not show the PAD result (see the log)"); self.append_log(error); return
        self.progress.setValue(100); self.progress.setFormat("PAD result loaded")
        self.show_result(view, pending["eb"])
        if pending["message"]:
            self.message.setText(pending["message"])

    def show_result(self, result, eb=None):
        """Display a complete result and match the EB slider to its energies."""
        self._stop_eb_sweep()
        s = result["s"]; count = len(energies(s))
        target = self.eb() if eb is None else float(eb)
        self._slider_min = s["eb_min"]; self._slider_step = s["de"]
        self.slider.blockSignals(True); self.slider.setRange(0, count - 1)
        self.slider.setValue(int(np.clip(round((target - self._slider_min) / self._slider_step), 0, count - 1)))
        self.slider.blockSignals(False); self.update_eb_input()
        self.result = result; self.update_mesh_result(result); self.redraw()

    def export_pad_result(self):
        """Save the displayed PAD result to a self-contained HDF5 file."""
        result = self.result
        if not pad_result_complete(result):
            QtWidgets.QMessageBox.information(
                self, "Save PAD result", "There is no complete PAD result to save. Press Calc first."); return
        s = result["s"]; name = result.get("structure_name") or "pad"
        folder = self.current_structure or STRUCTURES
        default = folder / f"{safe_name(name)}_hv{s['hv']:g}eV_{datetime.now():%Y%m%d-%H%M%S}.h5"
        selected, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save PAD result", str(default), PAD_RESULT_FILTER)
        if not selected:
            return
        path = Path(selected)
        if path.suffix.lower() not in (".h5", ".hdf5"):
            path = path.with_name(path.name + ".h5")
        structure = (self.current_structure if self.current_structure is not None
                     and self.current_structure.name == result.get("structure_name") else None)
        try:
            save_pad_result(path, result, structure)
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.critical(self, "Save PAD result", str(exc)); return
        self.message.setText(f"PAD result saved: {path}\n{format_file_size(path.stat().st_size)}")

    def import_pad_result(self):
        if self.thread:
            QtWidgets.QMessageBox.information(
                self, "Open PAD result", "Wait for the calculation to finish, or cancel it, first."); return
        selected, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open PAD result", str(self.current_structure or STRUCTURES), PAD_RESULT_FILTER)
        if selected:
            self.open_pad_result(Path(selected))

    def open_pad_result(self, path):
        """Show a saved PAD result with its calculation target and parameters."""
        try:
            result = load_pad_result(path)
        except ValueError as exc:
            QtWidgets.QMessageBox.critical(self, "Open PAD result", str(exc)); return False
        name = result.get("structure_name")
        target = STRUCTURES / safe_name(name) if name else None
        stored = (result.get("structure") or {}).get("cif")
        note = ""
        if target is None or structure_cif_path(target) is None:
            if target is None or not stored:
                QtWidgets.QMessageBox.warning(
                    self, "Open PAD result", "The file does not contain its crystal structure (CIF)."); return False
            answer = QtWidgets.QMessageBox.question(
                self, "Open PAD result",
                f"The calculation target “{target.name}” is not in the data folder.\n"
                f"Add it from the crystal structure saved in this file?\n\n{target}",
                QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.Cancel,
                QtWidgets.QMessageBox.StandardButton.Yes)
            if answer != QtWidgets.QMessageBox.StandardButton.Yes:
                return False
            try:
                self.add_target_from_result(target, result)
            except (OSError, ValueError, RuntimeError) as exc:
                QtWidgets.QMessageBox.critical(self, "Open PAD result", str(exc)); return False
            note = f"\nAdded calculation target: {target}"
        elif stored and structure_files(target).get("cif") != stored:
            note = "\nNote: the CIF in the data folder differs from the one saved in the file."
        self.refresh_structures()
        index = self.structure_combo.findData(str(target))
        if index < 0:
            QtWidgets.QMessageBox.warning(self, "Open PAD result", f"Calculation target not found: {target}"); return False
        self.structure_combo.setCurrentIndex(index)
        result["structure_name"] = target.name
        self.apply_settings(result["s"])
        self.display_in_background(result, result["s"], None,
                                   f"PAD result opened: {Path(path).name}\n"
                                   f"{target.name}, calculated {result.get('created') or 'earlier'}{note}")
        return True

    @staticmethod
    def add_target_from_result(target, result):
        """Create a calculation target from the structure stored in a result file."""
        structure = result["structure"]
        cif_name = Path(structure.get("cif_name") or "structure.cif").name
        if not cif_name.lower().endswith(".cif"):
            cif_name = "structure.cif"
        created = not target.exists()
        target.mkdir(parents=True, exist_ok=True)
        try:
            (target / cif_name).write_text(structure["cif"])
            if structure.get("openmx_template"):
                # Written after the CIF so it is not regenerated by cif2omx.
                (target / ".openmx.base.dat").write_text(structure["openmx_template"])
            (target / "results" / "openmx").mkdir(parents=True, exist_ok=True)
            write_settings(target / "parameter.json", result["s"])
            core.configure_structure(target)
        except Exception:
            if created:
                shutil.rmtree(target, ignore_errors=True)
            raise

    def euler_angles(self):
        return tuple(box.value() for box in (self.alpha, self.beta, self.gamma))

    def choose_euler_preset(self):
        """Set the incident angles (Euler a, b, c) from the beamline chosen in EulerPresetDialog."""
        dialog = EulerPresetDialog(self.euler_angles(), self)
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        angles = dialog.angles()
        if angles is not None:
            for box, value in zip((self.alpha, self.beta, self.gamma), angles):
                box.setValue(value)

    def show_preferences(self):
        dialog = PreferencesDialog(STRUCTURES, self)
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        selected = dialog.data_folder()
        if selected is None:
            QtWidgets.QMessageBox.warning(self, "Preferences", "Enter a data folder."); return
        self.change_data_folder(selected)

    def change_data_folder(self, folder):
        """Switch structure discovery and all OpenMX/PAD output to another folder."""
        folder = Path(folder).expanduser()
        if self.thread:
            QtWidgets.QMessageBox.information(
                self, "Preferences",
                "The data folder cannot be changed while a calculation is running."); return False
        if not folder.is_dir():
            if folder.exists():
                QtWidgets.QMessageBox.warning(self, "Preferences", f"Not a folder: {folder}"); return False
            answer = QtWidgets.QMessageBox.question(
                self, "Create data folder", f"The folder does not exist:\n{folder}\n\nCreate it?",
                QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
                QtWidgets.QMessageBox.StandardButton.Yes)
            if answer != QtWidgets.QMessageBox.StandardButton.Yes:
                return False
            try:
                folder.mkdir(parents=True)
            except OSError as exc:
                QtWidgets.QMessageBox.critical(self, "Create data folder", str(exc)); return False
        try:
            save_data_folder(folder)
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.critical(self, "Preferences", f"Could not save the data folder: {exc}"); return False
        if folder.resolve() == STRUCTURES:
            return True
        self._stop_eb_sweep()
        # refresh_structures() saves the previous target's parameter.json in
        # the old folder before selecting a target inside the new one.
        set_data_folder(folder)
        self.update_data_folder_button()
        self.refresh_structures()
        self.message.setText(f"Data folder: {STRUCTURES}\n" + self.message.text())
        return True

    def save_selected_settings(self):
        path = (self.current_structure / "parameter.json") if self.current_structure else None
        if path is None:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save settings", str(INITIAL_SETTINGS_FILE), "JSON files (*.json)")
        if not path: return
        try:
            write_settings(path, self.settings()); self.message.setText(f"Settings saved: {path}")
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.critical(self, "Settings save error", str(exc))

    def active_settings_path(self):
        return self.current_structure / "parameter.json" if self.current_structure else INITIAL_SETTINGS_FILE

    def save_calc_settings(self, settings):
        """Persist the latest Calc parameters for both startup and structure use."""
        write_settings(self.active_settings_path(), settings)
        if self.current_structure:
            write_settings(INITIAL_SETTINGS_FILE, settings)

    def create_upload_job(self):
        if not self.current_structure:
            QtWidgets.QMessageBox.information(self, "No calculation target", "Create or select a calculation target first."); return
        try:
            cif = structure_cif_path(self.current_structure)
            if cif is None:
                raise ValueError(f"CIF not found: {self.current_structure}")
            settings = validate_settings(self.settings()); config_text = json.dumps(settings, sort_keys=True, ensure_ascii=False)
            config_hash = hashlib.sha256(config_text.encode()).hexdigest(); calculation_id = f"calc_{datetime.now():%Y%m%d_%H%M%S}_{config_hash[:8]}"
            job = UPLOAD / f"{safe_name(self.current_structure.name)}_{calculation_id}"; job.mkdir(parents=True, exist_ok=False)
            shutil.copy2(cif, job / "structure.cif")
            write_settings(job / "parameter.json", settings); write_settings(job / "config.json", settings)
            shutil.copy2(APP_DIR / "spadexp_4d_core.py", job / "spadexp_4d_core.py")
            manifest = {"schema": "spadexp-4d-job", "schema_version": 1, "structure_name": self.current_structure.name,
                        "calculation_id": calculation_id, "config_sha256": config_hash,
                        "structure_cif_sha256": sha256_file(job / "structure.cif")}
            (job / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.critical(self, "Job creation error", str(exc)); return
        self.open_path(job); self.message.setText(f"Upload job created: {job.name}")

    def import_download_results(self):
        ensure_app_folders(); imported = []
        for job in sorted(DOWNLOAD.iterdir()):
            manifest_path = job / "manifest.json"
            if not job.is_dir() or not manifest_path.is_file(): continue
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8")); name = safe_name(manifest["structure_name"])
                calculation_id = safe_name(manifest["calculation_id"]); source_structure = STRUCTURES / name
                if not source_structure.is_dir(): raise ValueError(f"Structure folder not found: {name}")
                expected_cif = manifest.get("structure_cif_sha256")
                source_cif = structure_cif_path(source_structure)
                if source_cif is None:
                    raise ValueError(f"No CIF in structure folder: {source_structure}")
                if expected_cif and sha256_file(source_cif) != expected_cif:
                    raise ValueError("CIF SHA-256 does not match.")
                destination = source_structure / "results" / calculation_id
                if destination.exists(): continue
                destination.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(job, destination)
                active = {"structure_name": name, "calculation_id": calculation_id,
                          "path": calculation_id, "eigenstates": "eigenstates.hdf5",
                          "pad_raw": "pad_raw.hdf5", "manifest": "manifest.json"}
                (destination.parent / "active.json").write_text(json.dumps(active, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
                imported.append(f"{name}/{calculation_id}")
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                self.append_log(f"Skipped import: {job.name}: {exc}")
        self.refresh_results(); self.message.setText("Import complete: " + (", ".join(imported) if imported else "no new results"))

    def show_log(self):
        self.logview.setFocus(); self.logview.ensureCursorVisible()

    def toggle_bz(self, checked):
        self.show_bz.setChecked(checked); self.redraw()

    def show_about(self):
        QtWidgets.QMessageBox.about(self, "SPADExp4D", f"SPADExp4D\nVersion {__version__}")
    def direction_changed(self):
        a,b={"kx":("ky fixed (Å⁻¹)","kz fixed (Å⁻¹)"),"ky":("kx fixed (Å⁻¹)","kz fixed (Å⁻¹)"),"kz":("kx fixed (Å⁻¹)","ky fixed (Å⁻¹)")}[self.direction.currentText()]
        self.la.setText(a); self.lb.setText(b)
        # For kx and ky slices kz is set by the photoemission sphere at each EB
        # (PAD slice and band diagram), so the second coordinate is unused.
        sphere=self.direction.currentText()!="kz"
        self.b.setEnabled(not sphere); self.lb.setEnabled(not sphere)
        self.b.setToolTip("Not used: kz follows the photoemission sphere (hν, V0, EB) along the slice" if sphere else "")
    def update_eb_input(self):
        shift=self.ef_shift.value()
        maximum=self._slider_min+self.slider.maximum()*self._slider_step
        self.eb_input.blockSignals(True)
        minimum=self._slider_min+shift
        # Keep the Fermi-level display value editable even when the currently
        # cached EB interval starts above it after an EF shift.  Calc can then
        # be used after extending EB Begin to generate the missing states.
        if shift > 0:
            minimum=min(minimum,0.0)
        self.eb_input.setRange(minimum,maximum+shift)
        self.eb_input.setSingleStep(self._slider_step)
        self.eb_input.setValue(self.eb()+shift)
        self.eb_input.blockSignals(False)
    def slider_changed(self):
        self.eb_input.blockSignals(True); self.eb_input.setValue(self.eb()+self.ef_shift.value())
        self.eb_input.blockSignals(False)
        if self.slider.isSliderDown():
            self._eb_view_timer.start()  # while dragging, draw only the latest position
        else:
            self.update_eb_view()
    def set_eb_sweep_running(self, running):
        self.play_eb_button.setEnabled(not running and self.thread is None)
        self.stop_eb_button.setEnabled(running)
    def start_eb_sweep(self):
        if self.thread is not None:
            self.message.setText("EB playback is unavailable during calculation.")
            return
        if self.result is None:
            self.message.setText("No calculation result is available for playback. Press Calc first.")
            return
        self.eb_sweep_timer.stop()
        self.slider.blockSignals(True)
        self.slider.setValue(self.slider.minimum())
        self.slider.blockSignals(False)
        self.slider_changed()
        if self.slider.minimum() < self.slider.maximum():
            self.set_eb_sweep_running(True)
            self.eb_sweep_timer.start()
        else:
            self.set_eb_sweep_running(False)
    def advance_eb_sweep(self):
        next_index=self.slider.value()+1
        if next_index > self.slider.maximum():
            self._stop_eb_sweep()
            return
        self.slider.setValue(next_index)
        if next_index == self.slider.maximum():
            self._stop_eb_sweep()
    def _stop_eb_sweep(self):
        self.eb_sweep_timer.stop()
        self.set_eb_sweep_running(False)
    def stop_eb_sweep(self):
        self._stop_eb_sweep()
    def eb_input_changed(self, value=None):
        self._stop_eb_sweep()
        self.eb_input.interpretText()
        requested=self.eb_input.value()
        raw_eb=requested-self.ef_shift.value()
        raw_min=self._slider_min
        raw_max=self._slider_min+self.slider.maximum()*self._slider_step
        outside=raw_eb < raw_min-1e-9 or raw_eb > raw_max+1e-9
        index = int(np.clip(
            round((raw_eb-self._slider_min)/self._slider_step),
            self.slider.minimum(), self.slider.maximum(),
        ))
        self.slider.blockSignals(True)
        self.slider.setValue(index)
        self.slider.blockSignals(False)
        if outside:
            # No calculated slice exists for this requested shifted energy;
            # keep the entered value visible and show the nearest cached slice
            # until Calc is run with an expanded EB range.
            self.redraw()
            self.message.setText(
                f"EB′={requested:g} eV is outside the saved range. Extend EB Begin and press Calc.")
        else:
            self.slider_changed()
    def ef_shift_changed(self, value=None):
        self.update_eb_input()
        if self.result is not None:
            self.result['s']=dict(self.result['s'],ef_shift=self.ef_shift.value())
            self.redraw()
    def eb(self): return self._slider_min+self.slider.value()*self._slider_step
    def settings(self):
        return {"hv":self.hv.value(),"v0":self.v0.value(),"ef_shift":self.ef_shift.value(),"wfnc":self.wfnc.value(),
                "eb_min":self.eb_min.value(),"eb_max":self.eb_max.value(),"de":self.de.value(),
                "direction":self.direction.currentText(),"a":self.a.value(),"b":self.b.value(),
                "dk":self.dk.value(),"dE":self.dE.value(),
                "theta_max":self.theta_max.value(),
                "polarization":self.polarization.currentText(),
                "alpha":self.alpha.value(),"beta":self.beta.value(),"gamma":self.gamma.value(),
                "k_min":self.k_min.value(),"k_max":self.k_max.value(),
                "k_step":self.k_step.value(),
                "show_bz":self.show_bz.isChecked(),
                "pad_engine":self.pad_engine.currentData(), "use_first_bz_only":True, **{f"N{i+1}": box.value() for i,box in enumerate(self.nboxes)},
                "openmx":self.openmx_overrides(),"band_window":self.band_window_setting()}
    def save_settings(self):
        selected, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export settings parameter file",
            str(INITIAL_SETTINGS_FILE), "JSON files (*.json);;All files (*)"
        )
        if not selected: return
        try:
            write_settings(selected, self.settings())
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.critical(self, "Settings save error", str(exc))
            return
        self.message.setText(f"Settings exported: {selected}")
    def update_data_folder_button(self):
        self.data_folder_button.setText(path_tail(STRUCTURES, DATA_FOLDER_BUTTON_CHARS))
        self.data_folder_button.setToolTip(
            f"{STRUCTURES}\nClick to change the data folder (Settings > Preferences…)")
    def update_mesh_labels(self):
        ns=[box.value() for box in self.nboxes]
        for i,label in enumerate(self.glabels):
            g=core.G[i]; step=np.linalg.norm(g)/ns[i]
            label.setText(f"G{i+1} ({g[0]:.2f}, {g[1]:.2f}, {g[2]:.2f}) Å⁻¹   |G{i+1}|={np.linalg.norm(g):.2f} Å⁻¹   dK={step:.2f} Å⁻¹")
        self.total_points.setText(f"{ns[0]} × {ns[1]} × {ns[2]} = {math.prod(ns):,} points")
    def refresh_gaussian(self):
        if self.thread: return
        try: settings=validate_settings(self.settings())
        except ValueError: return
        if can_reuse_display(self.result, settings): self.calculate()
    def calculate(self):
        global LATTICE, B
        if self.thread: return
        self._stop_eb_sweep()
        selected_result = self.result_combo.currentData()
        try:
            if self.current_structure:
                configure_runtime_paths(self.current_structure)
                core.configure_structure(self.current_structure)
                LATTICE = core.LATTICE; B = core.G
            settings=validate_settings(self.settings())
            if can_reuse_display(self.result, settings):
                self.save_calc_settings(settings)
                self.result=prepare_display(self.result, settings)
                self.redraw()
                self.progress.setValue(100); self.progress.setFormat("Smoothing complete")
                self.message.setText(f"Gaussian update: {self.result['smoothing_seconds']:.3f}s\nReused saved PAD (no OpenMX or PAD calculation)")
                return
            if self.current_structure is None and structure_cif_path(DEFAULT_STRUCTURE) is None:
                # Never let OpenMX write outside a real structure in the data folder.
                raise ValueError(f"No calculation target in the data folder:\n{STRUCTURES}\n"
                                 "Import a CIF first (File > Import CIF…).")
            v0=(automatic_v0(SCF_H5, settings["v0"])
                if self.current_structure is None else (settings["v0"] if settings["v0"] > 0 else 20.0))
            if settings["hv"]+v0-settings["eb_max"] <= 0:
                raise ValueError("The kinetic energy inside the material is not positive at EB End.")
            plan=core.mesh_plan(settings)
            if not plan["reuse"]:
                answer=QtWidgets.QMessageBox.question(self,"First-BZ precomputation",
                    f"N={core.divisions(settings)}, total {plan['count']:,} k points.\n"
                    "Completed batches will be reused and missing points calculated by OpenMX.\nThis mesh will then be reused even if PAD conditions change.",
                    QtWidgets.QMessageBox.StandardButton.Yes |
                    QtWidgets.QMessageBox.StandardButton.No,
                    QtWidgets.QMessageBox.StandardButton.No)
                if answer != QtWidgets.QMessageBox.StandardButton.Yes: return
            self.save_calc_settings(settings)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self,"Settings/calculation preparation error",str(exc)); return
        old=self.eb(); self._slider_min=settings["eb_min"]; self._slider_step=settings["de"]
        self.slider.blockSignals(True); self.slider.setRange(0,len(energies(settings))-1)
        self.slider.setValue(int(np.clip(round((old-self._slider_min)/self._slider_step),0,self.slider.maximum())))
        self.slider.blockSignals(False); self.update_eb_input()
        self.cancel_background_display()
        self.result=None; self.logview.clear(); self.progress.setValue(0); self.progress.setFormat("%p%")
        if selected_result:
            self.message.setText(f"Preparing PAD with updated conditions using saved eigenstates: {Path(selected_result).name}")
        elif self.current_structure:
            self.message.setText(f"No calculation result: calculating {self.current_structure.name} locally")
        self.calc.setEnabled(False); self.calc_action.setEnabled(False); self.cancel_action.setEnabled(True)
        self.play_eb_button.setEnabled(False)
        self.thread=QtCore.QThread(self)
        self.worker=Worker(settings, self.current_structure, save_last=self.persist_session)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.status.connect(self.message.setText); self.worker.progress.connect(self.progress.setValue)
        self.worker.log.connect(self.append_log); self.worker.partial.connect(self.partial_result)
        self.worker.done.connect(self.finished); self.worker.error.connect(self.failed); self.worker.stopped.connect(self.stopped)
        self.worker.done.connect(self.thread.quit); self.worker.error.connect(self.thread.quit); self.worker.stopped.connect(self.thread.quit)
        self.thread.finished.connect(self.worker.deleteLater); self.thread.finished.connect(self.cleanup)
        self.thread.start()
    def stop_calculation(self):
        if self.worker:
            self.worker.request_stop(); self.cancel_action.setEnabled(False); self.message.setText("Stopping…")
    def stopped(self,result):
        if "v0" in result:
            self.result=result; self.update_mesh_result(result); self.redraw()
        self.progress.setFormat("Stopped"); self.message.setText("Stopped. Press Calc to resume using completed data.")
    def closeEvent(self,event):
        self._stop_eb_sweep()
        try:
            self.save_structure_settings()
        except (OSError, ValueError) as exc:
            self.append_log(f"Could not save PAD parameters on exit: {exc}")
        self.store_session()
        if self.thread:
            self._closing=True; self.stop_calculation(); event.ignore()
        else: event.accept()
    def append_log(self,text):
        self.logview.appendPlainText(text); bar=self.logview.verticalScrollBar(); bar.setValue(bar.maximum())
    def update_mesh_result(self, result):
        size = result.get('mesh_bytes') if result else None
        if size is None: return
        openmx = result.get('openmx_seconds', 0.0)
        conversion = result.get('conversion_seconds', 0.0)
        combined = result.get('openmx_hdf5_seconds', openmx + conversion)
        wall = result.get('mesh_seconds', combined)
        timing = (f"OpenMX {format_duration(openmx)} / HDF5 conversion {format_duration(conversion)}\n"
                  f"Total {format_duration(combined)} (wall time {format_duration(wall)})"
                  if openmx or conversion else f"Cache reused (checked in {format_duration(wall)})")
        path = Path(result.get('mesh_path', ''))
        display_name = str(Path(path.parent.name) / path.name) if path.name else 'unknown'
        lines=[display_name]
        if openmx or conversion:
            lines.append(timing)
        lines.append(format_file_size(size))
        self.mesh_result.setText("\n".join(lines))
        self.mesh_result.setToolTip(str(path))
    def partial_result(self,r):
        self.result=r; self.update_mesh_result(r); self.redraw()
    def cleanup(self):
        self.thread.deleteLater(); self.thread=self.worker=None
        self.calc.setEnabled(True); self.calc_action.setEnabled(True); self.cancel_action.setEnabled(False)
        self.set_eb_sweep_running(False)
        if self._closing: self.close()
    def finished(self,r):
        self.result=r; self.update_mesh_result(r)
        if self.current_structure:
            try:
                record = register_local_result(
                    self.current_structure, r["s"], r["mesh_path"], r)
                if record:
                    r["result_record"] = str(record)
                    self.refresh_results(record)
            except (OSError, ValueError) as exc:
                self.append_log(f"Failed to register calculation result: {exc}")
        self.progress.setValue(100); self.progress.setFormat("Complete: 100%")
        openmx = r.get("openmx_seconds", 0.0); conversion = r.get("conversion_seconds", 0.0)
        combined = r.get("openmx_hdf5_seconds", openmx + conversion)
        message = (f"OpenMX calculation: {format_duration(openmx)} / HDF5 conversion: {format_duration(conversion)}"
                   f" / total: {format_duration(combined)}"
                   if openmx or conversion else "OpenMX/HDF5 conversion: not run (first-BZ cache reused)")
        self.message.setText(f"Complete: {len(r['maps'])} EB slices\n{message}\nHDF5: {Path(r['mesh_path']).name}, {format_file_size(r['mesh_bytes'])}\nPAD: {format_duration(r.get('pad_seconds', 0))}\nMaximum distance to nearest eigenstate: {r['sampling_error']:.5f} Å⁻¹")
        self.redraw()
    def failed(self,text):
        self.message.setText("Calculation failed"); self.progress.setFormat("Calculation failed"); self.append_log(text)
        box=QtWidgets.QMessageBox(self); box.setIcon(QtWidgets.QMessageBox.Icon.Critical); box.setText("Calculation failed"); box.setDetailedText(text); box.exec()
    def map_bz(self,ax):
        if not self.show_bz.isChecked(): return
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        for i,j,poly in repeated_bz_polygons(xlim,ylim):
            p=np.vstack((poly,poly[0])); ax.plot(p[:,0],p[:,1],"w--",lw=.85,alpha=.8)
    def line_bz(self,ax,s):
        if not self.show_bz.isChecked(): return
        xlim=ax.get_xlim()
        if s["direction"]=="kz":
            period=np.linalg.norm(B[2]); n0=math.floor(xlim[0]/period-.5); n1=math.ceil(xlim[1]/period-.5)
            bounds=[(n+.5)*period for n in range(n0,n1+1)]
        else:
            fixed=s["a"]
            bounds=(bz_crossings_at_y(fixed,xlim) if s["direction"]=="kx"
                    else bz_crossings_at_x(fixed,xlim))
        for x in bounds: ax.axvline(x,color="tab:blue",ls="--",lw=.8,alpha=.75)
    def add_ek_axis(self,ax,s,v0):
        shifted_total=s["hv"]+v0+s.get("ef_shift",0.0)
        right=ax.secondary_yaxis("right",functions=(lambda eb_prime: shifted_total-eb_prime,
                                                     lambda ek: shifted_total-ek))
        right.set_ylabel("material Ek = hν+V0−(EB′−EF shift) (eV)")
        return right
    @staticmethod
    def eb_map_values(s, eb, data, name, axis):
        """kx-ky map of one EB slice; PAD outside the emission cone is hidden."""
        values=data[name]
        if name=='pad':
            vacuum=(s['hv']-eb-s['wfnc'])/H2M
            visible=(vacuum>0)&(axis[:,None]**2+axis[None,:]**2 <= vacuum*math.sin(math.radians(s['theta_max']))**2)
            values=np.where(visible,values,np.nan)
        return values

    @staticmethod
    def map_title(name, eb_display):
        """Title of the kx-ky map of the band weight ('band_weight') or the PAD."""
        prefix = "Band Slice " if name == "band_weight" else "PAD: "
        return f"{prefix}EB′={eb_display:g} eV"

    def update_eb_view(self):
        """EB change: update only the EB-dependent artists and blit them.

        The kx-ky maps (with colorbars and titles) and the EB marker lines are
        the only parts that depend on the selected EB; the full figure is
        rebuilt only when the result or display settings change.
        """
        view=self._eb_artists; r=self.result
        if view is None or r is not view['result']:
            self.redraw(); return
        s=r['s']; eb=round(self.eb(),8); data=r['maps'].get(eb)
        if data is None:
            self.redraw(); return
        eb_display=eb+s.get('ef_shift',0.0)
        for item in view['maps']:
            values=self.eb_map_values(s,eb,data,item['name'],view['axis'])
            item['image'].set_data(values)
            top=float(np.nanmax(values)) if np.isfinite(values).any() else 0.0
            item['image'].set_clim(0,top if top>0 else 1.0)  # as imshow(vmin=0) autoscales
            item['ax'].set_title(self.map_title(item['name'],eb_display))
        for line in view['lines']:
            line.set_ydata([eb_display,eb_display])
        if self._eb_background is None:
            self.canvas.draw_idle(); return  # first full draw still pending
        fig=self.canvas.figure
        self.canvas.restore_region(self._eb_background)
        for artist in view['animated']: fig.draw_artist(artist)
        self.canvas.blit(fig.bbox)

    def on_canvas_draw(self, event):
        """After a full draw keep the static background, then add the EB artists."""
        view=self._eb_artists
        if view is None or event.canvas is not self.canvas:
            self._eb_background=None; return
        fig=self.canvas.figure
        self._eb_background=self.canvas.copy_from_bbox(fig.bbox)
        for artist in view['animated']: fig.draw_artist(artist)

    def redraw(self):
        self._eb_artists=None; self._eb_background=None
        if self.result is None: return
        r=self.result
        if self.current_structure is not None and r.get('structure_name') not in (None, self.current_structure.name):
            return
        s=r['s']; shift=s.get('ef_shift',0.0)
        eb=round(self.eb(),8); eb_display=eb+shift; data=r['maps'].get(eb)
        fig=self.canvas.figure; fig.clear(); axes=fig.subplots(3,2).ravel()
        a0,a1,a2,a3,a4,a5=axes
        axis=cartesian_axis(s); extent=(axis[0],axis[-1],axis[0],axis[-1])
        maps=[]; lines=[]
        for ax,name,cmap,title in ((a0,'band_weight','viridis','Band weight'),(a1,'pad','inferno','PAD')):
            item=dict(ax=ax,name=name,title=title,image=None,colorbar=None); maps.append(item)
            if data is not None:
                values=self.eb_map_values(s,eb,data,name,axis)
                im=ax.imshow(values,origin='lower',extent=extent,cmap=cmap,aspect='equal',vmin=0)
                item['image']=im; item['colorbar']=fig.colorbar(im,ax=ax,label=title)
            else:
                ax.text(.5,.5,'Waiting for this EB slice',ha='center',transform=ax.transAxes)
            ax.set(xlim=extent[:2],ylim=extent[2:],xlabel='kx (Å⁻¹)',ylabel='ky (Å⁻¹)',
                   title=self.map_title(name,eb_display))
            self.map_bz(ax)
        ebs=energies(s); display_ebs=np.asarray(ebs)+shift; cuts=r['cuts']
        for ax,key,cmap,label,name in ((a2,'band_weight','viridis','Band weight','Band'),(a3,'pad','inferno','PAD','PAD')):
            im=None; title=f"{name} Slice {s['direction']}-EB"
            if s['direction']=='kz':
                ordered=[cuts[e] for e in sorted(cuts)]
                if ordered:
                    k=np.concatenate([c['k'] for c in ordered]); val=np.concatenate([c[key] for c in ordered])
                    energy=np.array([c['eb']+shift for c in ordered]); valid=np.isfinite(val)
                    if valid.any(): im=ax.scatter(k[valid],energy[valid],c=val[valid],cmap=cmap,s=18,vmin=0)
            else:
                values=np.full((len(ebs),len(axis)),np.nan)
                for i,e in enumerate(ebs):
                    if e in cuts: values[i]=cuts[e][key]
                if np.isfinite(values).any(): im=ax.pcolormesh(axis,display_ebs,values,shading='auto',cmap=cmap,vmin=0)
                ax.set_xlim(axis[0],axis[-1])
            if im is not None: fig.colorbar(im,ax=ax,label=label)
            else: ax.text(.5,.5,'No completed / visible points',ha='center',transform=ax.transAxes)
            # EF shift moves the spectra inside a fixed requested display
            # window.  Moving these limits by the same amount made the shift
            # cancel visually even though the energy coordinates had changed.
            ax.set(ylim=(s['eb_max'],s['eb_min']),xlabel=f"{s['direction']} (Å⁻¹)",
                   ylabel='EB′ = EB + EF shift (eV)',title=title)
            lines.append(ax.axhline(eb_display,color='cyan',ls='--')); self.add_ek_axis(ax,s,r['v0']); self.line_bz(ax,s)
        # Bands at the k points the PAD slice samples (kz on the
        # photoemission sphere at each EB), see pad_cut_bands().
        line=r.get('line') or {}
        if line.get('k') is not None and line.get('bands') is not None:
            for band in np.asarray(line['bands']).T: a4.plot(line['k'],-band+shift,'k',lw=.8)
        points=np.asarray(line.get('points',np.zeros((0,2))),dtype=float).reshape(-1,2)
        if len(points): a4.plot(points[:,0],-points[:,1]+shift,'k.',ms=2)
        if s['direction']!='kz': a4.set_xlim(axis[0],axis[-1])
        a4.set(ylim=(s['eb_max'],s['eb_min']),xlabel=f"{s['direction']} (Å⁻¹)",
               ylabel='EB′ = EB + EF shift (eV)',title='Band Diagram (PAD cut)')
        lines.append(a4.axhline(eb_display,color='red',ls='--')); self.add_ek_axis(a4,s,r['v0']); self.line_bz(a4,s)
        kx=np.linspace(axis[0],axis[-1],401)
        norm=matplotlib.colors.Normalize(s['eb_min'],s['eb_max']); cmap=matplotlib.colormaps['turbo']
        for e in np.linspace(s['eb_min'],s['eb_max'],min(41,len(ebs))):
            z2=(s['hv']+r['v0']-e)/H2M-kx*kx; valid=z2>=0
            a5.plot(kx[valid],np.sqrt(z2[valid]),color=cmap(norm(e+shift)))
        kz_max=math.sqrt(max((s['hv']+r['v0']-s['eb_min'])/H2M,0.0))*1.04
        if self.show_bz.isChecked():
            kx_bounds,kz_bounds=trajectory_bz_boundaries((axis[0],axis[-1]),kz_max)
            for boundary in kx_bounds:
                a5.axvline(boundary,color='.35',ls='--',lw=.7,alpha=.8)
            for boundary in kz_bounds:
                a5.axhline(boundary,color='.35',ls='--',lw=.7,alpha=.8)
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),ax=a5,label='EB′ = EB + EF shift (eV)')
        a5.set(xlim=(axis[0],axis[-1]),ylim=(0,max(kz_max,1e-6)),
               xlabel='kx (Å⁻¹)',ylabel='kz (Å⁻¹)',title='kx-kz photoelectron trajectories')
        if all(item['image'] is not None for item in maps):
            animated=[a for item in maps for a in (item['ax'],item['colorbar'].ax)]+lines
            for artist in animated: artist.set_animated(True)
            self._eb_artists=dict(result=r,maps=maps,lines=lines,animated=animated,axis=axis)
        self.canvas.draw_idle()


def runtime_executable_status():
    """Return configured runtime executables that are missing or not executable."""
    required = {
        "OpenMX": core.OPENMX,
        "cif2omx": core.CONVERTER,
        "OpenMX postprocessor": core.POSTPROC,
    }
    return [(label, Path(path)) for label, path in required.items()
            if not Path(path).is_file() or not os.access(path, os.X_OK)]


def show_runtime_path_dialog(missing):
    """Explain how to configure missing OpenMX tools before showing the GUI."""
    details = "\n".join(f"• {label}: {path}" for label, path in missing)
    QtWidgets.QMessageBox.critical(
        None,
        "OpenMX path configuration required",
        "Please set the path to the OpenMX executable and related tools.\n\n"
        f"Missing or unavailable files:\n{details}\n\n"
        "Set SPADEXP4D_OPENMX and SPADEXP4D_CIF2OMX, or add the tools to PATH, "
        "then restart SPADExp4D.",
    )


DATA_FOLDER_SUMMARY = (
    "The data folder is where SPADExp4D keeps your calculation targets and "
    "everything calculated for them."
)
# (indent level, name, description) rows of the "What is stored" overview.
DATA_FOLDER_LAYOUT = (
    (0, "<data folder>/", ""),
    (1, "graphite/", "one subfolder per calculation target"),
    (2, "structure.cif", "imported crystal structure"),
    (2, "parameter.json", "settings of this target"),
    (2, "cache/", "OpenMX runs and eigenstates (often 100s of MB)"),
    (3, "last_pad_result.h5", "last PAD result, shown again at the next start"),
    (2, "results/", "list of finished calculations"),
)
DATA_FOLDER_SWITCH_NOTE = (
    "Switching does not move, copy or delete any data. The previous folder "
    "stays as it is, and you can switch back at any time."
)


DATA_FOLDER_BUTTON_CHARS = 32  # sidebar button shows only the last N characters
BAND_WINDOW_HELP = (
    "Energy range (relative to EF) of the bands stored in the k-space HDF5. "
    "All keeps every band. A window that covers the EB range of interest "
    "reduces storage and memory several-fold, e.g. for NC spinors; bands "
    "outside it are missing from PAD. Changing it recomputes the k-space HDF5."
)
OPENMX_RECOMMENDED_HELP = (
    "Set recommended SCF values for this cell (never weaker than the template):\n"
    "• scf.Kgrid: k-point spacing of 0.2 Å⁻¹ or finer\n"
    "• scf.energycutoff ≥ 220 Ry, scf.criterion ≤ 1e-6 Hartree, scf.maxIter ≥ 200\n"
    "• scf.ElectronicTemperature 300 K\n"
    "• scf.XcType as in the template (it matches the pseudopotentials)\n"
    "Each keyword's tooltip shows its template and recommended values."
)
OPENMX_PARAMETERS_HELP = (
    "OpenMX keywords of the selected calculation target, starting from its "
    ".openmx.base.dat (generated from the CIF). Bold keywords differ from the "
    "template; they are saved in parameter.json, and the next Calc runs a new "
    "SCF and k-space HDF5 for them."
)


def path_tail(path, limit):
    """Return at most `limit` characters of a path, keeping its end."""
    text = str(path)
    return text if len(text) <= limit else "…" + text[-(limit - 1):]


def calculation_targets(folder):
    """Names of the calculation targets (subfolders with a CIF) in a data folder."""
    return sorted(path.name for path in Path(folder).iterdir()
                  if path.is_dir() and structure_cif_path(path) is not None)


def data_folder_status(path, current=None):
    """Return (usable, level, text) describing a candidate data folder.

    level is "ok", "warning" or "error"; usable is False only when switching
    to the path is impossible.
    """
    if path is None:
        return False, "error", "Enter a folder path or click Browse…."
    path = Path(path).expanduser()
    if path.exists() and not path.is_dir():
        return False, "error", "This path is a file, not a folder."
    if not path.is_dir():
        return True, "warning", ("This folder does not exist yet. It will be created "
                                 f"inside {nearest_existing_folder(path)}.")
    try:
        targets = calculation_targets(path)
        free = shutil.disk_usage(path).free
    except OSError as exc:
        return False, "error", f"This folder cannot be read: {exc.strerror or exc}"
    is_current = current is not None and path.resolve() == Path(current).resolve()
    state = "This is the current data folder." if is_current else "Existing folder."
    if targets:
        shown = ", ".join(targets[:4]) + (", …" if len(targets) > 4 else "")
        noun = "calculation target" if len(targets) == 1 else "calculation targets"
        contents = f"{len(targets)} {noun}: {shown}"
    else:
        contents = "No calculation targets yet (add one with File > Import CIF…)"
    return True, "ok", f"{state}\n{contents}  ·  {free / 1024**3:.1f} GiB free"


def nearest_existing_folder(path):
    """Return the closest existing ancestor of a path, or the home folder."""
    for candidate in (Path(path).expanduser(), *Path(path).expanduser().parents):
        if candidate.is_dir():
            return candidate
    return Path.home()


def choose_data_folder(parent, start):
    selected = QtWidgets.QFileDialog.getExistingDirectory(
        parent, "Choose Data Folder", str(nearest_existing_folder(start)))
    return Path(selected) if selected else None


def prompt_for_data_folder(missing, parent=None):
    """Ask for a data folder when the configured one is unavailable at startup."""
    missing = Path(missing)
    default = DEFAULT_DATA_FOLDER
    offer_default = default.is_dir() and default.resolve() != missing.resolve()
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(QtWidgets.QMessageBox.Icon.Warning)
    box.setWindowTitle("Data folder not found")
    box.setText("SPADExp4D cannot find its data folder.")
    choices = "choose another folder" + (" or use the default folder" if offer_default else "")
    box.setInformativeText(
        f"{missing}\n\n"
        "If it is on an external or network drive, connect the drive and click "
        f"Retry. Otherwise {choices}.\n\n{DATA_FOLDER_SUMMARY}"
        + (f"\n\nDefault folder: {default}" if offer_default else ""))
    choose = box.addButton("Choose Folder…", QtWidgets.QMessageBox.ButtonRole.AcceptRole)
    use_default = (box.addButton("Use Default Folder", QtWidgets.QMessageBox.ButtonRole.ActionRole)
                   if offer_default else None)
    retry = box.addButton("Retry", QtWidgets.QMessageBox.ButtonRole.ActionRole)
    box.addButton("Quit", QtWidgets.QMessageBox.ButtonRole.RejectRole)
    box.setDefaultButton(choose)
    while True:
        box.exec()
        clicked = box.clickedButton()
        if clicked is not None and clicked is retry:
            if missing.is_dir():
                return missing
            box.setText("The data folder is still not available.")
            continue
        if clicked is not None and clicked is use_default:
            return default
        if clicked is not choose:
            return None
        selected = choose_data_folder(parent, missing)
        if selected is not None:
            return selected


def format_euler(angles):
    return ", ".join(f"{name} = {value:g}°" for name, value in zip("abc", angles))


class EulerPresetDialog(QtWidgets.QDialog):
    """List of beamlines whose Euler angles a, b, c can be applied with OK."""

    def __init__(self, current=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Incident angle presets")
        self.setMinimumWidth(360)
        self.list = QtWidgets.QListWidget()
        for name, angles in EULER_PRESETS:
            item = QtWidgets.QListWidgetItem(name)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, angles)
            item.setToolTip(format_euler(angles))
            self.list.addItem(item)
        self.list.setMinimumHeight(self.list.sizeHintForRow(0) * (len(EULER_PRESETS) + 1))
        self.detail = QtWidgets.QLabel()
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.ok_button = buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject)
        self.list.currentRowChanged.connect(self.selection_changed)
        self.list.itemDoubleClicked.connect(lambda _item: self.accept())
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel("Beamline"))
        layout.addWidget(self.list)
        layout.addWidget(self.detail)
        layout.addWidget(buttons)
        # Start on the preset that matches the current angles, if any.
        current_row = next((row for row, (_, angles) in enumerate(EULER_PRESETS)
                            if current is not None and np.allclose(angles, current)), 0)
        self.list.setCurrentRow(current_row)
        self.selection_changed()

    def selection_changed(self, *_):
        angles = self.angles()
        self.ok_button.setEnabled(angles is not None)
        self.detail.setText("" if angles is None else format_euler(angles))

    def angles(self):
        item = self.list.currentItem()
        if item is None:
            return None
        return tuple(float(value) for value in item.data(QtCore.Qt.ItemDataRole.UserRole))


class PreferencesDialog(QtWidgets.QDialog):
    """Application preferences stored in the initialization file."""

    STATUS_SYMBOL = {"ok": "✓", "warning": "⚠", "error": "✗"}
    STATUS_COLOR = {"ok": None, "warning": "#b26a00", "error": "#c62828"}

    def __init__(self, data_folder, parent=None):
        super().__init__(parent)
        self.current = Path(data_folder)
        self.setWindowTitle("Preferences")
        self.setMinimumWidth(660)

        heading = QtWidgets.QLabel("Data folder")
        font = heading.font(); font.setBold(True); font.setPointSizeF(font.pointSizeF() * 1.25)
        heading.setFont(font)
        summary = QtWidgets.QLabel(DATA_FOLDER_SUMMARY); summary.setWordWrap(True)

        self.folder_edit = QtWidgets.QLineEdit(str(data_folder))
        self.folder_edit.setPlaceholderText("Folder path")
        browse = QtWidgets.QPushButton("Browse…"); browse.clicked.connect(self.browse)
        row = QtWidgets.QHBoxLayout(); row.addWidget(self.folder_edit, 1); row.addWidget(browse)
        self.status = QtWidgets.QLabel(); self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)

        contents = QtWidgets.QGroupBox("What is stored in a data folder")
        grid = QtWidgets.QGridLayout(contents)
        grid.setHorizontalSpacing(18); grid.setVerticalSpacing(2); grid.setColumnStretch(1, 1)
        fixed = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
        for line, (level, name, description) in enumerate(DATA_FOLDER_LAYOUT):
            item = QtWidgets.QLabel(name); item.setFont(fixed); item.setIndent(16 * level)
            grid.addWidget(item, line, 0)
            if description:
                grid.addWidget(QtWidgets.QLabel(description), line, 1)

        note = QtWidgets.QLabel(DATA_FOLDER_SWITCH_NOTE); note.setWordWrap(True)
        try:
            ini = INITIAL_SETTINGS_FILE.relative_to(PROJECT)
        except ValueError:
            ini = INITIAL_SETTINGS_FILE
        saved = QtWidgets.QLabel(f"This setting is saved in {ini}.")
        saved.setEnabled(False)  # secondary, greyed text

        reveal = "Show in Finder" if sys.platform == "darwin" else "Open Folder"
        self.open_button = QtWidgets.QPushButton(reveal); self.open_button.clicked.connect(self.open_folder)
        self.default_button = QtWidgets.QPushButton("Use Default")
        self.default_button.setToolTip(str(DEFAULT_DATA_FOLDER))
        self.default_button.clicked.connect(lambda: self.folder_edit.setText(str(DEFAULT_DATA_FOLDER)))
        cancel = QtWidgets.QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        self.ok_button = QtWidgets.QPushButton("Use This Folder")
        self.ok_button.setDefault(True); self.ok_button.clicked.connect(self.accept)
        buttons = QtWidgets.QHBoxLayout()
        for widget in (self.open_button, self.default_button):
            buttons.addWidget(widget)
        buttons.addStretch(1); buttons.addWidget(cancel); buttons.addWidget(self.ok_button)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(heading); layout.addWidget(summary); layout.addSpacing(4)
        layout.addLayout(row); layout.addWidget(self.status); layout.addSpacing(6)
        layout.addWidget(contents); layout.addWidget(note); layout.addWidget(saved)
        layout.addSpacing(6); layout.addLayout(buttons)
        self.folder_edit.textChanged.connect(self.update_status)
        self.update_status()

    def update_status(self):
        path = self.data_folder()
        usable, level, text = data_folder_status(path, self.current)
        color = self.STATUS_COLOR[level]
        self.status.setStyleSheet(f"color: {color};" if color else "")
        self.status.setText(f"{self.STATUS_SYMBOL[level]}  {text}")
        self.ok_button.setEnabled(usable)
        self.ok_button.setText("Close" if usable and path is not None and path.is_dir()
                               and path.resolve() == self.current.resolve() else "Use This Folder")
        self.open_button.setEnabled(path is not None and path.is_dir())
        self.default_button.setEnabled(path is None or path != DEFAULT_DATA_FOLDER)

    def open_folder(self):
        path = self.data_folder()
        if path is not None and path.is_dir():
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(path)))

    def browse(self):
        selected = choose_data_folder(self, self.folder_edit.text().strip() or STRUCTURES)
        if selected is not None:
            self.folder_edit.setText(str(selected))

    def data_folder(self):
        text = self.folder_edit.text().strip()
        if not text:
            return None
        path = Path(text).expanduser()
        return path if path.is_absolute() else PROJECT / path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--smoke-test',action='store_true')
    parser.add_argument('--prepare',action='store_true',help='Compute only the complete first-BZ eigenstate database')
    parser.add_argument('--settings',type=Path,default=INITIAL_SETTINGS_FILE)
    parser.add_argument('--structure',default='graphite',help='Structure folder name for --check/--prepare')
    parser.add_argument('--pad-engine', choices=('cpp', 'python'), help='PAD backend (default: cpp)')
    args=parser.parse_args()
    settings=read_settings(args.settings)
    if args.pad_engine: settings["pad_engine"] = args.pad_engine
    data_folder = configured_data_folder()
    if args.check or args.prepare:
        if not data_folder.is_dir():
            parser.error(f'Data folder not found: {data_folder}\n'
                         f'Choose one in the GUI or set "{DATA_FOLDER_KEY}" in {INITIAL_SETTINGS_FILE}.')
        set_data_folder(data_folder)
        structure_path = STRUCTURES / safe_name(args.structure)
        if not structure_cif_path(structure_path):
            parser.error(f'No CIF found in structure: {structure_path}')
        configure_runtime_paths(structure_path)
        core.configure_structure(structure_path)
    if args.check:
        plan=core.mesh_plan(settings)
        print('4D settings:',args.settings)
        for i,n in enumerate(core.divisions(settings)):
            print(f'G{i+1}={core.G[i]} A^-1 (step={np.linalg.norm(core.G[i])/n:.6f}); N{i+1}={n}')
        print('Total k points:',plan['count'],'Cache reusable:',plan['reuse']); print(plan['path'])
        print('HDF5 size:', format_file_size(plan['path'].stat().st_size) if plan['path'].is_file() else 'not created'); return
    if args.prepare:
        started=time.monotonic(); path,elapsed=core.ensure_mesh(settings,log=print); total=time.monotonic()-started; print(path)
        print('K-space total time:',format_duration(total)); print('OpenMX time:',format_duration(elapsed)); print('HDF5 size:',format_file_size(path.stat().st_size)); return
    app=QtWidgets.QApplication(sys.argv)
    missing = runtime_executable_status()
    if missing:
        show_runtime_path_dialog(missing)
        return 1
    if not data_folder.is_dir():
        selected = prompt_for_data_folder(data_folder)
        if selected is None:
            return 1
        save_data_folder(selected); data_folder = selected
    set_data_folder(data_folder)
    window=Window(settings, session=True); window.show()
    if args.smoke_test:
        window.persist_session=False  # check the restore path without saving
        QtCore.QTimer.singleShot(500,window.close)
    sys.exit(app.exec())


if __name__=='__main__': main()
