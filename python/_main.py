#!/usr/bin/env python3
"""4-D (kx, ky, kz, E) viewer using a reusable first-BZ eigenstate database."""

# SPADExp4D 1.1.0: integrated Python, C++ and per-Structure caches.

from __future__ import annotations

__version__ = "1.1.0"

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
    if not 0 <= normalized["alpha"] <= 360 or not 0 <= normalized["gamma"] <= 360:
        raise ValueError("Euler angles alpha and gamma must be between 0 and 360 degrees.")
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

def pad_text(source, destination, logfile, hv, v0, wf, eb_min, eb_max, de,
             gaussian_de, field):
    ex, ey, ez = field
    return f"""&Control
Calculation PAD
Log_file {logfile}
Console_log true
Output_file {destination}
/
&Atomic-wfn
Potential_file {ATOMIC_POTENTIALS}
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
Atomic_orbitals_file {ATOMIC_ORBITALS}
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

def frac(kcart):
    return np.asarray(kcart) @ np.linalg.inv(B)

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
    smoothed = np.zeros_like(values)
    for i, weight in enumerate(kernel):
        core.check_stop(stop)
        section = [slice(None)] * values.ndim
        section[axis] = slice(i, i + values.shape[axis])
        smoothed += weight * padded[tuple(section)]
    return smoothed


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
    with h5py.File(path, 'w') as out:
        inp = out.create_group('Input')
        for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
            inp[name] = h5py.ExternalLink(str(grid.path.resolve()), '/Input/' + name)
        for name in ('Output', 'Mesh'):
            out[name] = h5py.ExternalLink(str(grid.path.resolve()), '/' + name)
        k = inp.create_group('Kpath')
        k.attrs.update(Dimension=np.int32(2), Curved=np.uint8(1),
            Origin=frac((0, 0, math.sqrt((s['hv']+v0-requested[0])/H2M))),
            Xvector=frac((1, 0, 0)), Yvector=frac((0, 1, 0)),
            Xcount=count, Ycount=count, Count=count*count,
            Xrange=[s['k_min'], s['k_max']], Yrange=[s['k_min'], s['k_max']],
            minN=core.BANDS[0], maxN=core.BANDS[1])
        k.create_dataset('Coordinates', data=cart @ np.linalg.inv(B))
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
            inp.write_text(pad_text(source, output, OUTPUT / f'batch_{code}_PAD.log',
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
    attrs = dict(Dimension=np.int32(2), Curved=np.uint8(1),
                 Origin=frac((0, 0, math.sqrt(max(s['hv']+v0-eb, 0)/H2M))),
                 Xvector=frac((1, 0, 0)), Yvector=frac((0, 1, 0)),
                 Xcount=np.int32(count), Ycount=np.int32(count),
                 Xrange=np.array((s['k_min'], s['k_max'])), Yrange=np.array((s['k_min'], s['k_max'])))
    with tempfile.TemporaryDirectory(prefix='pad_', dir=INPUT) as temporary:
        directory = Path(temporary); source = directory / 'states.hdf5'; output = directory / 'pad.hdf5'
        error = grid.write_path(source, cart, attrs)
        inp = directory / 'pad.in'
        inp.write_text(pad_text(source, output, OUTPUT / f'slice_{code}_PAD.log',
                               s['hv'], v0, s['wfnc'], engine_eb, engine_eb, s['de'], sigma, field))
        with h5py.File(source, 'r') as f:
            bands = (f['Output/Band'][:] - grid.ef) * PAD_EH
        # Raw densities: each eigenstate occupies one half-open EB bin.
        active = np.floor((bands + eb) / s['de'] + 0.5) == 0
        if not active.any():
            pad = np.zeros((count, count))
            if log: log(f'[PAD skipped] EB={eb:g}: intensity set to 0 because the EB bin contains no states')
        else:
            core.run_command([SPADEXP], OUTPUT / f'slice_{code}_stdout.log', directory,
                             stdin=inp, expected=output, on_output=log, stop=stop)
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


def fixed_band(grid, s):
    axis = cartesian_axis(s); coord = {'kx': 0, 'ky': 1, 'kz': 2}[s['direction']]
    cart = np.zeros((len(axis), 3)); cart[:, coord] = axis
    cart[:, [i for i in range(3) if i != coord]] = (s['a'], s['b'])
    bands, error = grid.interpolated_bands(cart)
    return {'k': axis, 'bands': bands, 'error': float(error.max(initial=0)),
            'method': 'periodic trilinear scalar interpolation'}


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

    def __init__(self, settings, structure_path=None):
        super().__init__(); self.s = settings; self.structure_path = Path(structure_path) if structure_path else None
        self.structure_name = self.structure_path.name if self.structure_path else "graphite"
        self.stop_event = threading.Event()

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
            self.log.emit(f'[OpenMX all bands] {core.BANDS[1] - core.BANDS[0] + 1} bands '
                          f'(indices {core.BANDS[0]}–{core.BANDS[1]}, including states above EF)')
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
            result['line'] = fixed_band(grid, self.s)
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
            self.stage(100, 'Complete'); self.done.emit(result)
        except core.Cancelled:
            self.stopped.emit(prepare_display(result, self.s))
        except Exception:
            self.error.emit(traceback.format_exc())
        finally:
            if grid is not None: grid.close()


class Window(QtWidgets.QMainWindow):
    def __init__(self, initial_settings):
        super().__init__(); self.setWindowTitle("SPADExp4D — first-BZ eigenstate cache"); self.resize(1500,1050)
        self.initial_settings = validate_settings(initial_settings)
        self.result = self.thread = self.worker = None; self._closing = False
        self.eb_sweep_timer = QtCore.QTimer(self)
        self.eb_sweep_timer.setInterval(250)
        self.eb_sweep_timer.timeout.connect(self.advance_eb_sweep)
        self.current_structure = None; self.build(); self.direction_changed()
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
        def add_row(label, widget):
            nonlocal row
            lab=(label if isinstance(label,QtWidgets.QLabel) else QtWidgets.QLabel(label)); lab.setStyleSheet("margin: 0; padding: 0;")
            g.addWidget(lab,row,0); g.addWidget(widget,row,1); row+=1
        def add_full_row(widget):
            nonlocal row
            g.addWidget(widget,row,0,1,2); row+=1
        def inline(*widgets):
            holder=QtWidgets.QWidget(); layout=QtWidgets.QHBoxLayout(holder)
            layout.setContentsMargins(0,0,0,0); layout.setSpacing(3)
            for widget in widgets: layout.addWidget(widget,1)
            return holder
        self.structure_combo = QtWidgets.QComboBox()
        self.structure_combo.setToolTip("Calculation targets under structures/<structure_name>")
        self.structure_combo.currentIndexChanged.connect(self.structure_changed)
        add_row("Calculation target", self.structure_combo)
        self.result_combo = QtWidgets.QComboBox()
        self.result_combo.setToolTip("Saved calculation results for the selected structure")
        self.result_combo.currentIndexChanged.connect(self.result_selection_changed)
        add_row("Calculation result", self.result_combo)
        self.nboxes=[]; self.glabels=[]
        for i in range(3):
            box=QtWidgets.QSpinBox(); box.setRange(1,256); box.setValue(initial[f"N{i+1}"])
            box.setToolTip("First-BZ divisions. Changing these selects or calculates the corresponding eigenstate cache.")
            self.nboxes.append(box)
            box.valueChanged.connect(self.update_mesh_labels)
        add_row("N1, N2, N3",inline(*self.nboxes))
        for i in range(3):
            label=QtWidgets.QLabel(); label.setWordWrap(False); self.glabels.append(label)
            add_full_row(label)
        self.total_points=QtWidgets.QLabel(); self.total_points.setWordWrap(True)
        add_row("Total k-point count",self.total_points)
        self.mesh_result=QtWidgets.QLabel("Not available"); self.mesh_result.setWordWrap(True)
        self.mesh_result.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        add_row("K-space HDF5",self.mesh_result)
        self.update_mesh_labels()
        self.pad_engine=QtWidgets.QComboBox()
        self.pad_engine.addItem("C++ batch", "cpp"); self.pad_engine.addItem("Python per EB (comparison)", "python")
        self.pad_engine.setCurrentIndex(self.pad_engine.findData(initial["pad_engine"]))
        add_row("PAD calculation",self.pad_engine)
        self.hv=self.spin(initial["hv"],.01,10000); self.v0=self.spin(initial["v0"],-1000,1000)
        self.ef_shift=self.spin(initial["ef_shift"],-1000,1000); self.wfnc=self.spin(initial["wfnc"],0,100)
        for label,box in (("hν (eV)",self.hv),("V0 (eV; ≤0:auto)",self.v0),
                          ("EF shift (eV)",self.ef_shift),("wfnc (eV)",self.wfnc)): add_row(label,box)
        self.ef_shift.setToolTip(
            "Display-only shift added to the calculated binding energy. "
            "Moves the band plot and energy slice within the fixed display range.")
        self.eb_min=self.spin(initial["eb_min"],-1000,1000); self.eb_max=self.spin(initial["eb_max"],-1000,1000); self.de=self.spin(initial["de"],.001,100)
        self.eb_min.setToolTip("Negative EB includes eigenstates above the original EF in the PAD calculation.")
        add_row("EB Begin, End, Step (eV)",inline(self.eb_min,self.eb_max,self.de))
        self.direction=QtWidgets.QComboBox(); self.direction.addItems(("kx","ky","kz")); self.direction.setCurrentText(initial["direction"]); self.direction.currentTextChanged.connect(self.direction_changed)
        self.la,self.lb=QtWidgets.QLabel(),QtWidgets.QLabel(); self.a,self.b=self.spin(initial["a"],-20,20),self.spin(initial["b"],-20,20)
        add_row("Slice direction",self.direction); add_row(self.la,self.a); add_row(self.lb,self.b)
        self.polarization=QtWidgets.QComboBox(); self.polarization.addItems(POLARIZATIONS); self.polarization.setCurrentText(initial["polarization"])
        self.polarization.setToolTip("unpolarized adds the c+ and c- PAD intensities")
        self.alpha=self.spin(initial["alpha"],0,360); self.beta=self.spin(initial["beta"],0,180); self.gamma=self.spin(initial["gamma"],0,360)
        add_row("Polarization",self.polarization); add_row("Euler a, b, c (°)",inline(self.alpha,self.beta,self.gamma))
        self.dk=self.spin(initial["dk"],0,100); self.dE=self.spin(initial["dE"],0,100)
        # Gaussian widths affect the displayed PAD, but changing a spin box
        # must not launch repeated processing while the user is typing.  The
        # values are applied by Calc (or the explicit View menu action).
        for box in (self.dk,self.dE):
            box.setToolTip("Gaussian standard deviation applied to saved PAD. 0 disables smoothing; values outside the calculation range are treated as 0.")
        add_row("Gaussian dK (Å⁻¹)",self.dk); add_row("Gaussian dE (eV)",self.dE)
        self.theta_max=self.spin(initial["theta_max"],.001,90)
        add_row("Maximum angle θmax (°)",self.theta_max)
        self.show_bz=QtWidgets.QCheckBox("Show Brillouin zone"); self.show_bz.setChecked(initial["show_bz"]); self.show_bz.toggled.connect(self.redraw)
        add_row("Boundary",self.show_bz)
        self.k_min=self.spin(initial["k_min"],-K_LIMIT,K_LIMIT); self.k_max=self.spin(initial["k_max"],-K_LIMIT,K_LIMIT)
        self.k_step=self.spin(initial["k_step"],.001,2*K_LIMIT)
        add_row("KMin, KMax, Kstep (Å⁻¹)",inline(self.k_min,self.k_max,self.k_step))
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
        outer.addWidget(graphics,1)
        self.refresh_structures()
    def build_menus(self):
        file_menu = self.menuBar().addMenu("File")
        self.add_menu_action(file_menu, "Import CIF…", self.import_cif)
        self.add_menu_action(file_menu, "Open calculation target folder…", self.open_structure_folder)
        self.add_menu_action(file_menu, "Save parameter.json", self.save_selected_settings)
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
        entries = structure_entries(); selected = self.current_structure.name if self.current_structure else None
        self.structure_combo.blockSignals(True); self.structure_combo.clear()
        if not entries: self.structure_combo.addItem("(Place a CIF file under structures)", None)
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
                    write_settings(settings_path, self.initial_settings)
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
            core.configure_structure(None)
            LATTICE = core.LATTICE; B = core.G
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
            self, "Import CIF", "Folder name to create under structures:",
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
        self.eb_input.blockSignals(False); self.redraw()
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
                "pad_engine":self.pad_engine.currentData(), "use_first_bz_only":True, **{f"N{i+1}": box.value() for i,box in enumerate(self.nboxes)}}
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
        self.result=None; self.logview.clear(); self.progress.setValue(0); self.progress.setFormat("%p%")
        if selected_result:
            self.message.setText(f"Preparing PAD with updated conditions using saved eigenstates: {Path(selected_result).name}")
        elif self.current_structure:
            self.message.setText(f"No calculation result: calculating {self.current_structure.name} locally")
        self.calc.setEnabled(False); self.calc_action.setEnabled(False); self.cancel_action.setEnabled(True)
        self.play_eb_button.setEnabled(False)
        self.thread=QtCore.QThread(self)
        self.worker=Worker(settings, self.current_structure)
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
    def redraw(self):
        if self.result is None: return
        r=self.result
        if self.current_structure is not None and r.get('structure_name') not in (None, self.current_structure.name):
            return
        s=r['s']; shift=s.get('ef_shift',0.0)
        eb=round(self.eb(),8); eb_display=eb+shift; data=r['maps'].get(eb)
        fig=self.canvas.figure; fig.clear(); axes=fig.subplots(3,2).ravel()
        a0,a1,a2,a3,a4,a5=axes
        axis=cartesian_axis(s); extent=(axis[0],axis[-1],axis[0],axis[-1])
        for ax,name,cmap,title in ((a0,'band_weight','viridis','Band weight'),(a1,'pad','inferno','PAD')):
            if data is not None:
                values=data[name]
                if name=='pad':
                    vacuum=(s['hv']-eb-s['wfnc'])/H2M
                    visible=(vacuum>0)&(axis[:,None]**2+axis[None,:]**2 <= vacuum*math.sin(math.radians(s['theta_max']))**2)
                    values=np.where(visible,values,np.nan)
                im=ax.imshow(values,origin='lower',extent=extent,cmap=cmap,aspect='equal',vmin=0)
                fig.colorbar(im,ax=ax,label=title)
            else:
                ax.text(.5,.5,'Waiting for this EB slice',ha='center',transform=ax.transAxes)
            ax.set(xlim=extent[:2],ylim=extent[2:],xlabel='kx (Å⁻¹)',ylabel='ky (Å⁻¹)',
                   title=f'{title}: EB′={eb_display:g} eV — first-BZ grid reference')
            self.map_bz(ax)
        ebs=energies(s); display_ebs=np.asarray(ebs)+shift; cuts=r['cuts']
        for ax,key,cmap,label in ((a2,'band_weight','viridis','Band weight'),(a3,'pad','inferno','PAD')):
            im=None
            if s['direction']=='kz':
                ordered=[cuts[e] for e in sorted(cuts)]
                if ordered:
                    k=np.concatenate([c['k'] for c in ordered]); val=np.concatenate([c[key] for c in ordered])
                    energy=np.array([c['eb']+shift for c in ordered]); valid=np.isfinite(val)
                    if valid.any(): im=ax.scatter(k[valid],energy[valid],c=val[valid],cmap=cmap,s=18,vmin=0)
                title=f'{label}: sampled kz trajectory'
            else:
                values=np.full((len(ebs),len(axis)),np.nan)
                for i,e in enumerate(ebs):
                    if e in cuts: values[i]=cuts[e][key]
                if np.isfinite(values).any(): im=ax.pcolormesh(axis,display_ebs,values,shading='auto',cmap=cmap,vmin=0)
                ax.set_xlim(axis[0],axis[-1]); title=f'{label} cut from EB-dependent maps'
            if im is not None: fig.colorbar(im,ax=ax,label=label)
            else: ax.text(.5,.5,'No completed / visible points',ha='center',transform=ax.transAxes)
            # EF shift moves the spectra inside a fixed requested display
            # window.  Moving these limits by the same amount made the shift
            # cancel visually even though the energy coordinates had changed.
            ax.set(ylim=(s['eb_max'],s['eb_min']),xlabel=f"{s['direction']} (Å⁻¹)",
                   ylabel='EB′ = EB + EF shift (eV)',title=title)
            ax.axhline(eb_display,color='cyan',ls='--'); self.add_ek_axis(ax,s,r['v0']); self.line_bz(ax,s)
        line=r.get('line')
        if line:
            for band in line['bands'].T: a4.plot(line['k'],-band+shift,'k',lw=.8)
        fixed_names={'kx':('ky','kz'),'ky':('kx','kz'),'kz':('kx','ky')}[s['direction']]
        fixed_text=f"{fixed_names[0]}={s['a']:.3f}, {fixed_names[1]}={s['b']:.3f} Å⁻¹"
        a4.set(xlim=(axis[0],axis[-1]),ylim=(s['eb_max'],s['eb_min']),xlabel=f"{s['direction']} (Å⁻¹)",
               ylabel='EB′ = EB + EF shift (eV)',title=f"Fixed-coordinate band cut: {fixed_text} — scalar interpolation")
        a4.axhline(eb_display,color='red',ls='--'); self.add_ek_axis(a4,s,r['v0']); self.line_bz(a4,s)
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
    if args.check or args.prepare:
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
    window=Window(settings); window.show()
    if args.smoke_test: QtCore.QTimer.singleShot(500,window.close)
    sys.exit(app.exec())


if __name__=='__main__': main()
