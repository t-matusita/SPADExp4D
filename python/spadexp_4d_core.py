"""First-BZ eigenstate database, independent of the 3-D viewer and its caches."""
from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import time

import h5py
import numpy as np

APP_DIR = Path(__file__).resolve().parent
PROJECT = APP_DIR.parent
STRUCTURES = PROJECT / 'structures'
DEFAULT_STRUCTURE = STRUCTURES / 'graphite'
BUILD = PROJECT / '.build'
CACHE = DEFAULT_STRUCTURE / 'cache'
DATA = PROJECT / 'data/openmx'
TEMPLATE = DEFAULT_STRUCTURE / '.openmx.base.dat'
RESTART_NAME = 'graphite_scf'
RESTART = CACHE / 'scf_rst'
RESTART_FINGERPRINT = None
STRUCTURE_RESULTS = None
STRUCTURE_CACHE = None


def configured_executable(environment, command):
    """Resolve a tool from an environment override, then PATH."""
    configured = os.environ.get(environment)
    if configured:
        return Path(configured).expanduser()
    discovered = shutil.which(command)
    return Path(discovered) if discovered else Path(command)


OPENMX = configured_executable('SPADEXP4D_OPENMX', 'openmx')
CONVERTER = configured_executable('SPADEXP4D_CIF2OMX', 'cif2omx')
POSTPROC = PROJECT / 'bin/openmx_postproc_4d_plus'
LATTICE = np.array([[2.462, 0, 0], [-1.231, 2.1322, 0], [0, 0, 6.707]])
G = 2 * np.pi * np.linalg.inv(LATTICE).T
BANDS = (0, 15)  # Replaced from the selected Structure's complete basis.
CHUNK_POINTS = 1024  # Legacy HDF5 checkpoint layout only; NOT an OpenMX batch limit.
SCHEMA = 2  # Version 2 stores every band in the selected atomic basis.
SHIFTS = np.array(list(itertools.product(range(-2, 3), repeat=3)))


def normalized_cif_text(path):
    """Remove CIF oxidation-state suffixes from elemental type symbols.

    cif2omx expects a bare element such as ``P`` or ``Ta``.  Common CIF
    exporters instead emit type symbols such as ``P0+``, ``Ta4+`` or ``S2-``.
    Atom-site labels such as ``P0`` are deliberately left unchanged.
    """
    source = Path(path).read_text(errors='replace')
    charged_element = re.compile(
        r'(?<![A-Za-z])([A-Z][a-z]?)(?:(?:\d+(?:\.\d+)?)[+-]|[+-](?:\d+(?:\.\d+)?)?)(?![A-Za-z0-9])')
    return charged_element.sub(r'\1', source)


def template_has_undefined_species(path):
    try:
        text = Path(path).read_text(errors='replace')
    except OSError:
        return True
    match = re.search(r'<Definition\.of\.Atomic\.Species(.*?)Definition\.of\.Atomic\.Species>',
                      text, flags=re.S)
    return match is None or bool(re.search(r'\bundefined\b', match.group(1), flags=re.I))


def complete_band_range(path):
    """Return every zero-based band represented by an OpenMX basis.

    For a non-spin calculation OpenMX has one band per atomic basis function.
    No occupation or energy-window estimate is used, so unoccupied states above
    EF remain available to PAD after an EF shift.
    """
    text = Path(path).read_text(errors='replace')
    spin_match = re.search(r'(?mi)^\s*scf\.SpinPolarization\s+(\S+)', text)
    spin = spin_match.group(1).lower() if spin_match else 'off'
    if spin != 'off':
        raise ValueError('SPADExp4D full-band storage currently supports scf.SpinPolarization off.'
                         f' Current setting: {spin}')
    species_match = re.search(
        r'<Definition\.of\.Atomic\.Species\s*(.*?)\s*Definition\.of\.Atomic\.Species>',
        text, flags=re.S)
    atoms_match = re.search(
        r'<Atoms\.SpeciesAndCoordinates\s*(.*?)\s*Atoms\.SpeciesAndCoordinates>',
        text, flags=re.S)
    if species_match is None or atoms_match is None:
        raise ValueError(f'Could not read basis or atom list from OpenMX input: {path}')
    degeneracy = {'s': 1, 'p': 3, 'd': 5, 'f': 7}
    orbitals = {}
    for line in species_match.group(1).splitlines():
        fields = line.split('#', 1)[0].split()
        if len(fields) < 2:
            continue
        basis = fields[1].rsplit('-', 1)[-1]
        groups = re.findall(r'([spdf])(\d+)', basis)
        if not groups or ''.join(kind + count for kind, count in groups) != basis:
            raise ValueError(f'Unrecognized OpenMX basis: {fields[0]} {fields[1]}')
        orbitals[fields[0]] = sum(degeneracy[kind] * int(count)
                                  for kind, count in groups)
    atom_labels = []
    for line in atoms_match.group(1).splitlines():
        fields = line.split('#', 1)[0].split()
        if len(fields) >= 2 and fields[0].isdigit():
            atom_labels.append(fields[1])
    missing = sorted(set(atom_labels) - set(orbitals))
    if not atom_labels or missing:
        detail = ', '.join(missing) if missing else 'no atoms found'
        raise ValueError(f'Could not determine basis for atom species: {detail}')
    count = sum(orbitals[label] for label in atom_labels)
    if count < 1:
        raise ValueError('OpenMX contains zero bands.')
    return (0, count - 1)


def _cif_cell(path):
    """Read unit-cell lengths and angles from a CIF."""
    values = {}
    for line in Path(path).read_text(errors='replace').splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        key = fields[0].lower()
        if key in {'_cell_length_a', '_cell_length_b', '_cell_length_c',
                   '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'}:
            try:
                values[key] = float(re.sub(r'\([^)]*\)', '', fields[1]))
            except ValueError:
                pass
    required = ('_cell_length_a', '_cell_length_b', '_cell_length_c',
                '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma')
    if any(key not in values for key in required):
        raise ValueError(f'Could not read lattice parameters from CIF: {path}')
    a, b, c = (values[f'_cell_length_{x}'] for x in 'abc')
    alpha, beta, gamma = (np.deg2rad(values[f'_cell_angle_{x}']) for x in ('alpha', 'beta', 'gamma'))
    avec = np.array((a, 0.0, 0.0))
    bvec = np.array((b * np.cos(gamma), b * np.sin(gamma), 0.0))
    cz = c * np.sqrt(max(0.0, 1.0 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2
                         + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
    cvec = np.array((c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma), cz))
    return np.vstack((avec, bvec, cvec))


def configure_structure(structure_path=None):
    """Select the CIF/OpenMX context used by subsequent mesh calculations."""
    global TEMPLATE, RESTART_NAME, RESTART, RESTART_FINGERPRINT
    global DATA, LATTICE, G, BANDS, CACHE, STRUCTURE_RESULTS, STRUCTURE_CACHE
    if structure_path is None:
        structure_path = DEFAULT_STRUCTURE
    folder = Path(structure_path).resolve()
    canonical = folder / 'structure.cif'
    cifs = [canonical] if canonical.is_file() else sorted(p for p in folder.glob('*.cif') if p.is_file())
    if not cifs:
        raise FileNotFoundError(f'No CIF file in structure folder: {folder}')
    cif = cifs[0]
    STRUCTURE_RESULTS = folder / 'results'
    STRUCTURE_CACHE = folder / 'cache'
    CACHE = STRUCTURE_CACHE
    CACHE.mkdir(parents=True, exist_ok=True)
    template = folder / '.openmx.base.dat'
    if ((not template.is_file()) or template.stat().st_mtime_ns < cif.stat().st_mtime_ns
            or template_has_undefined_species(template)):
        if not CONVERTER.is_file():
            raise FileNotFoundError(
                'cif2omx not found. Add it to PATH or specify the executable '
                'with SPADEXP4D_CIF2OMX.')
        normalized = folder / '.normalized.cif'
        normalized.write_text(normalized_cif_text(cif))
        try:
            completed = subprocess.run([str(CONVERTER), str(normalized), '-o', str(template)],
                                       cwd=folder, capture_output=True, text=True)
        finally:
            normalized.unlink(missing_ok=True)
        if completed.returncode != 0 or not template.is_file():
            raise RuntimeError('Could not generate OpenMX input from CIF: ' + completed.stderr[-2000:])
        if template_has_undefined_species(template):
            raise RuntimeError('Could not determine OpenMX PAO/VPS from CIF element symbols. '
                               'Check _atom_site_type_symbol.')
    TEMPLATE = template
    RESTART_NAME = folder.name + '_scf'
    RESTART = STRUCTURE_CACHE / 'scf_rst'
    # A locally generated restart accelerates later SCF work, but must not
    # change the identity of an already completed eigenstate mesh.
    RESTART_FINGERPRINT = []
    data_root = OPENMX.parent.parent / 'DFT_DATA19'
    DATA = data_root if data_root.is_dir() else PROJECT / 'data/openmx'
    LATTICE = _cif_cell(cif); G = 2 * np.pi * np.linalg.inv(LATTICE).T
    BANDS = complete_band_range(TEMPLATE)
    return {'name': folder.name, 'cif': cif, 'template': template}


class Cancelled(Exception):
    pass


def check_stop(stop=None):
    if stop and stop():
        raise Cancelled()


def run_command(args, log, cwd, stdin=None, expected=None, on_output=None, stop=None):
    """Stream bounded GUI batches; retain full logs and stop silent children."""
    check_stop(stop)
    process = source = None
    started = time.monotonic()
    try:
        source = Path(stdin).open('rb') if stdin else None
        with Path(log).open('w') as output:
            process = subprocess.Popen(list(map(str, args)), cwd=cwd, stdin=source,
                                       stdout=output, stderr=subprocess.STDOUT, start_new_session=True)
            with Path(log).open(errors='replace') as reader:
                while True:
                    check_stop(stop)
                    code = process.poll()
                    chunk = reader.read()
                    if on_output and chunk.strip():
                        on_output('\n'.join(chunk.rstrip().splitlines()[-60:]))
                    if code is not None:
                        break
                    time.sleep(0.1)
    finally:
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        if source:
            source.close()
    if code != 0:
        raise RuntimeError(f'{args[0]} failed ({code}):\n' + Path(log).read_text(errors='replace')[-4000:])
    if expected and (not Path(expected).is_file() or not Path(expected).stat().st_size):
        tail = Path(log).read_text(errors='replace')[-4000:]
        raise RuntimeError(f'Missing output: {expected}\n{tail}')
    return time.monotonic() - started


def json_write(path, value):
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    temporary.replace(path)


def json_read(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def stamp(path):
    p = Path(path); st = p.stat()
    return {'path': str(p.resolve()), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def divisions(settings):
    values = []
    for key in ('N1', 'N2', 'N3'):
        value = settings[key]
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f'{key} must be an integer from 1 to 256.') from exc
        if isinstance(value, bool) or not np.isfinite(numeric) or not numeric.is_integer() or not 1 <= numeric <= 256:
            raise ValueError(f'{key} must be an integer from 1 to 256.')
        values.append(int(numeric))
    return tuple(values)


def mesh_request(settings):
    # No photoemission / plotting parameter enters this identity. These fixed
    # physical fingerprints only guard against reusing a different material.
    if not OPENMX.is_file():
        raise FileNotFoundError(
            'Please set the path to the OpenMX executable. '
            'Use SPADEXP4D_OPENMX or add openmx to PATH.')
    return {'schema': SCHEMA, 'divisions': list(divisions(settings)), 'lattice': LATTICE.tolist(),
            'bands': list(BANDS), 'template': digest(TEMPLATE),
            'data': {str(p.relative_to(DATA)): digest(p) for p in sorted(DATA.rglob('*')) if p.is_file()},
            'restart': (RESTART_FINGERPRINT if RESTART_FINGERPRINT is not None else
                        ([stamp(p) for p in sorted(RESTART.iterdir()) if p.is_file()]
                         if RESTART.is_dir() else [])),
            'openmx': stamp(OPENMX), 'scf_restart': 'pp'}


def request_code(request):
    return hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()[:20]


def fold_first_bz(points):
    """Wigner-Seitz representatives for this hexagonal, normal-z lattice."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    centered = points - np.rint(points)
    best = np.full(len(points), np.inf)
    result = np.empty_like(centered)
    for shift in SHIFTS:
        candidate = centered - shift
        distance = np.sum((candidate @ G)**2, axis=1)
        improve = distance < best - 1e-14
        result[improve] = candidate[improve]; best[improve] = distance[improve]
    return result


def grid_points(ns, start=0, stop=None):
    """Half-open Gamma-centered periodic grid; N3 varies fastest."""
    end = int(np.prod(ns)) if stop is None else stop
    indices = np.stack(np.unravel_index(np.arange(start, end), ns), axis=-1)
    return fold_first_bz(indices / np.asarray(ns))


@contextmanager
def file_lock(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('The same cache is being used by another 4D calculation.') from exc
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def input_text(stem, points):
    text = TEMPLATE.read_text()
    text = re.sub(r'(?m)^System\.Name\s+\S+', f'System.Name {stem}', text)
    text = re.sub(r'(?m)^DATA\.PATH\s+.*$', f'DATA.PATH {DATA}', text)
    restart_file = RESTART / f'{RESTART_NAME}.crst_check'
    restart_block = (f'scf.restart pp\nscf.restart.filename {RESTART_NAME}\norbitalOpt.Force.Skip on'
                     if restart_file.is_file() else 'scf.restart on')
    if re.search(r'(?m)^scf\.restart\s+', text):
        text = re.sub(r'(?m)^scf\.restart\s+.*$', restart_block, text, count=1)
    else:
        text += '\n' + restart_block + '\n'
    # The converter accepts an explicit MO list. Its 1-D envelope only records
    # the batch length; no fictitious line is sent through the preprocessor.
    block = (f'SPADExp.dimension 1\nSPADExp.curved false\nSPADExp.origin 0 0 0\n'
             f'<SPADExp.range\n 1 0 0 0 1 {len(points)}\nSPADExp.range>\n'
             f'SPADExp.minN {BANDS[0]}\nSPADExp.maxN {BANDS[1]}')
    if re.search(r'SPADExp\.dimension.*?SPADExp\.maxN\s+\d+', text, flags=re.S):
        text = re.sub(r'SPADExp\.dimension.*?SPADExp\.maxN\s+\d+', block, text, count=1, flags=re.S)
    else:
        text += '\n' + block + '\n'
    return text + (f'\nMO.fileout on\nMO.Nkpoint {len(points)}\n<MO.kpoint\n'
                   + '\n'.join(' '.join(f'{v:.12f}' for v in p) for p in points) + '\nMO.kpoint>\n')


def lcao_paths(file):
    paths = []
    file['Output/LCAO'].visititems(lambda name, obj: paths.append('Output/LCAO/' + name)
                                 if isinstance(obj, h5py.Dataset) else None)
    return sorted(paths)


def expected_lcao_paths(file):
    """Build the expected radial-orbital dataset list from OpenMX input data."""
    def text(value):
        return value.decode() if isinstance(value, (bytes, np.bytes_)) else str(value)

    species_labels = [text(value) for value in file['Input/Atomic.Species/Labels'][:]]
    species_orbits = [text(value) for value in file['Input/Atomic.Species/Orbits'][:]]
    orbit_by_species = dict(zip(species_labels, species_orbits))
    atom_labels = [text(value) for value in file['Input/Atoms.SpeciesAndCoordinates/Labels'][:]]
    expected = []
    for atom_index, label in enumerate(atom_labels, 1):
        orbit_text = orbit_by_species[label]
        orbit_groups = re.findall(r'([spdf])(\d+)', orbit_text)
        if not orbit_groups or ''.join(kind + count for kind, count in orbit_groups) != orbit_text:
            raise ValueError(f'Unrecognized OpenMX basis: {label} {orbit_text}')
        for kind, count in orbit_groups:
            for radial_index in range(int(count)):
                expected.append(f'Output/LCAO/{atom_index}_{label}/{kind}{radial_index}')
    return sorted(expected)


def states_ok(file, count, compressed=False):
    try:
        band_count = BANDS[1] - BANDS[0] + 1
        paths = lcao_paths(file)
        return (file['Output/Band'].shape == (count, band_count)
                and str(file['Output'].attrs['Spin']) == 'off'
                and int(file['Output'].attrs['BandIndexStartZeroBased']) == BANDS[0]
                and int(file['Output'].attrs['BandIndexEndZeroBased']) == BANDS[1]
                and paths == expected_lcao_paths(file)
                and all(file[p].shape[:2] == (count, band_count) and file[p].shape[-1] == 2
                        and (not compressed or (file[p].dtype == np.float32 and file[p].compression == 'gzip'))
                        for p in paths)
                and 'Output/Potential/V0' in file)
    except (KeyError, ValueError, OSError):
        return False


def cache_ok(path, request):
    try:
        with h5py.File(path, 'r') as f:
            count = int(np.prod(request['divisions']))
            return (bool(f.attrs.get('Complete', False))
                    and f.attrs.get('Request') == json.dumps(request, sort_keys=True)
                    and f['Mesh/FirstBZCoordinates'].shape == (count, 3)
                    and np.all(f['Mesh/CompletedChunks'][:]) and states_ok(f, count, True))
    except (OSError, ValueError, KeyError):
        return False


def mesh_plan(settings):
    request = mesh_request(settings)
    folder = CACHE / ('mesh_' + request_code(request))
    path = folder / 'eigenstates.hdf5'
    return {'request': request, 'folder': folder, 'path': path, 'reuse': cache_ok(path, request),
            'count': int(np.prod(request['divisions']))}


def completed_openmx(folder, stem, text):
    record = json_read(folder / 'openmx.json')
    out = folder / f'{stem}.out'
    if not out.is_file() or not out.stat().st_size:
        return False
    try:
        with (folder / 'openmx.log').open('rb') as log:
            log.seek(max(0, log.seek(0, 2) - 4096))
            normal = b'The calculation was normally finished.' in log.read()
        output_matches = bool(record and record.get('output') == stamp(out))
        input_matches = bool(record and record.get('input') == hashlib.sha256(text.encode()).hexdigest())
        # Saving a newly generated SCF restart changes only the restart lines
        # of the next input.  Reuse the already completed MO output when the
        # original recorded input differs solely in those lines.
        dat = folder / f'{stem}.dat'
        if not input_matches and record and dat.is_file():
            stored = dat.read_text(errors='replace')
            stored_recorded = record.get('input') == hashlib.sha256(stored.encode()).hexdigest()
            def without_restart(value):
                return re.sub(r'(?m)^(?:scf\.restart(?:\.filename)?|orbitalOpt\.Force\.Skip)\s+.*\n?', '', value)
            input_matches = stored_recorded and without_restart(stored) == without_restart(text)
        return normal and output_matches and input_matches
    except OSError:
        return False


def preserve_structure_restart(folder, stem):
    """Keep one reusable SCF restart beside the selected structure."""
    if STRUCTURE_CACHE is None or RESTART.parent != STRUCTURE_CACHE:
        return
    source = Path(folder) / f'{stem}_rst'
    if not source.is_dir() or RESTART.is_dir():
        return
    temporary = RESTART.with_name(RESTART.name + '.tmp')
    if temporary.exists():
        raise RuntimeError(f'A temporary resume folder remains: {temporary}')
    shutil.copytree(source, temporary)
    for path in sorted(temporary.iterdir()):
        if path.name.startswith(stem):
            path.rename(path.with_name(RESTART_NAME + path.name[len(stem):]))
    temporary.replace(RESTART)


def ensure_postprocessor(log=None, stop=None):
    """Build only the private heap converter; never replace the shared tools."""
    BUILD.mkdir(parents=True, exist_ok=True)
    with file_lock(BUILD / 'cpp_build.lock'):
        run_command(['make', '-C', PROJECT / 'cpp', '-j2', POSTPROC],
                    BUILD / 'postprocessor_build.log', PROJECT, expected=POSTPROC,
                    on_output=log, stop=stop)


def calculate_chunk(folder, points, log=None, stop=None, conversion=None):
    ensure_postprocessor(log, stop)
    folder.mkdir(parents=True, exist_ok=True)
    stem = folder.name; dat = folder / f'{stem}.dat'; raw = folder / f'{stem}.hdf5'
    text = input_text(stem, points)
    check_stop(stop)
    elapsed = 0.0
    conversion_elapsed = 0.0
    # A fully converted raw HDF5 is sufficient to resume the compact import.
    # Validate it before checking the mutable SCF restart lines or deleting any
    # OpenMX output, so a postprocessing retry can never trigger recomputation.
    try:
        with h5py.File(raw, 'r') as existing:
            raw_valid = (states_ok(existing, len(points))
                         and np.allclose(existing['Input/Kpath/Coordinates'][:], points,
                                         rtol=0, atol=2e-11))
    except (OSError, KeyError, ValueError):
        raw_valid = False
    if raw_valid:
        return raw, elapsed, conversion_elapsed
    if not completed_openmx(folder, stem, text):
        restart_file = RESTART / f'{RESTART_NAME}.crst_check'
        if restart_file.is_file():
            link = folder / f'{RESTART_NAME}_rst'
            if not link.exists():
                link.symlink_to(RESTART, target_is_directory=True)
            if not link.is_symlink() or link.resolve() != RESTART.resolve():
                raise RuntimeError(f'Unexpected restart path: {link}')
        (folder / 'openmx.json').unlink(missing_ok=True)
        (folder / f'{stem}.out').unlink(missing_ok=True)
        raw.unlink(missing_ok=True)
        dat.write_text(text)
        if not OPENMX.is_file():
            raise FileNotFoundError(
                'openmx not found. Add it to PATH or specify the executable '
                'with SPADEXP4D_OPENMX.')
        elapsed = run_command([OPENMX, dat.name], folder / 'openmx.log', folder,
                              expected=folder / f'{stem}.out', on_output=log, stop=stop)
        # Normal completion must be established before recording reusable output.
        tail = (folder / 'openmx.log').read_text(errors='replace')[-4096:]
        if 'The calculation was normally finished.' not in tail:
            raise RuntimeError('Could not verify normal OpenMX completion.')
        json_write(folder / 'openmx.json', {'input': hashlib.sha256(text.encode()).hexdigest(),
                                           'output': stamp(folder / f'{stem}.out')})
    preserve_structure_restart(folder, stem)
    conversion_record = {'output': stamp(folder / f'{stem}.out'), 'postprocessor': stamp(POSTPROC)}
    valid = False
    try:
        with h5py.File(raw, 'r') as f:
            valid = (states_ok(f, len(points)) and json_read(folder / 'conversion.json') == conversion_record
                     and np.allclose(f['Input/Kpath/Coordinates'][:], points, rtol=0, atol=2e-11))
    except (OSError, KeyError, ValueError):
        pass
    if not valid:
        raw.unlink(missing_ok=True)
        def postproc_output(text):
            if log:
                log(text)
            if conversion:
                # The private postprocessor emits this marker while parsing
                # OpenMX's text output.  Keep the parser tolerant of a tail
                # containing several repeated log snapshots.
                for current, total in re.findall(r'SPADExp conversion: kpoint=(\d+)/(\d+)', text):
                    conversion(int(current), int(total))
        conversion_started = time.monotonic()
        run_command([POSTPROC, dat], folder / 'postproc.log', folder, expected=raw,
                    on_output=postproc_output, stop=stop)
        conversion_elapsed = time.monotonic() - conversion_started
        with h5py.File(raw, 'r') as f:
            if not states_ok(f, len(points)) or not np.allclose(f['Input/Kpath/Coordinates'][:], points, rtol=0, atol=2e-11):
                raise RuntimeError(f'Incomplete eigenstates: {raw}')
        json_write(folder / 'conversion.json', conversion_record)
    return raw, elapsed, conversion_elapsed


def initialize_database(file, raw, count):
    band_count = BANDS[1] - BANDS[0] + 1
    for key, value in raw.attrs.items():
        file.attrs[key] = value
    inp = file.create_group('Input')
    for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
        raw.copy('Input/' + name, inp, name=name)
    out = file.create_group('Output'); out.attrs.update(raw['Output'].attrs)
    raw.copy('Output/Potential', out, name='Potential')
    out.create_dataset('Band', shape=(count, band_count), dtype='f8',
                       chunks=(min(256, count), band_count), compression='gzip', shuffle=True)
    for path in lcao_paths(raw):
        dataset = raw[path]
        shape = (count,) + dataset.shape[1:]
        parent, name = path.rsplit('/', 1)
        file.require_group(parent).create_dataset(name, shape=shape, dtype='f4',
            chunks=(min(128, count),) + shape[1:], compression='gzip', compression_opts=4, shuffle=True, fletcher32=True)


def ensure_mesh(settings, log=None, progress=None, stop=None, conversion=None, timing=None):
    plan = mesh_plan(settings); request = plan['request']; folder = plan['folder']; final = plan['path']
    folder.mkdir(parents=True, exist_ok=True)
    with file_lock(folder / 'mesh.lock'):
        if cache_ok(final, request):
            if log: log(f'[OpenMX skipped] Reusing first-BZ data for {plan["count"]:,} points: {final}')
            return final, 0.0
        partial = folder / 'eigenstates.partial.hdf5'
        count = plan['count']; ns = request['divisions']; batches = (count + CHUNK_POINTS - 1) // CHUNK_POINTS
        reuse_partial = False
        try:
            with h5py.File(partial, 'r') as f:
                reuse_partial = (f.attrs['Request'] == json.dumps(request, sort_keys=True)
                                 and f['Mesh/CompletedChunks'].shape == (batches,)
                                 and f['Mesh/FirstBZCoordinates'].shape == (count, 3)
                                 and ('Output' not in f or states_ok(f, count, True)))
        except (OSError, KeyError):
            pass
        if not reuse_partial:
            with h5py.File(partial, 'w') as f:
                f.attrs.update(Request=json.dumps(request, sort_keys=True), Complete=False, Schema=SCHEMA)
                mesh = f.create_group('Mesh'); mesh.attrs['Divisions'] = ns
                mesh.attrs['ReciprocalVectors_Ainv'] = G
                mesh.attrs['Ordering'] = 'N3-fast; index=(i1*N2+i2)*N3+i3; Gamma included; endpoints excluded'
                mesh.create_dataset('FirstBZCoordinates', shape=(count, 3), dtype='f8')
                mesh.create_dataset('CompletedChunks', data=np.zeros(batches, dtype=bool))
        elapsed = 0.0
        with h5py.File(partial, 'r+') as f:
            def store_raw(raw_path, selected, points):
                # Write bounded slabs, but commit completion only after ALL slabs
                # are durable. A stopped import retries the same OpenMX job/output.
                with h5py.File(raw_path, 'r') as raw:
                    if 'Output' not in f:
                        initialize_database(f, raw, count)
                    if not np.isclose(f['Output'].attrs['EF_Eh'], raw['Output'].attrs['EF_Eh'], rtol=0, atol=1e-8):
                        raise RuntimeError('The SCF Fermi energy does not match the saved data.')
                    row = 0
                    for batch in selected:
                        start = batch * CHUNK_POINTS; end = min(count, start + CHUNK_POINTS)
                        check_stop(stop)
                        f['Output/Band'][start:end] = raw['Output/Band'][row:row+end-start]
                        for path in lcao_paths(raw):
                            check_stop(stop)
                            f[path][start:end] = raw[path][row:row+end-start].astype(np.float32)
                        f['Mesh/FirstBZCoordinates'][start:end] = points[row:row+end-start]
                        row += end-start
                f.flush()
                f['Mesh/CompletedChunks'][selected] = True; f.flush()
                if progress:
                    flags = f['Mesh/CompletedChunks'][:]
                    done = sum(min(CHUNK_POINTS, count-i*CHUNK_POINTS) for i in np.flatnonzero(flags))
                    progress(done, count)

            pending = np.flatnonzero(~f['Mesh/CompletedChunks'][:])
            # Upgrade an interrupted legacy run without repeating any finished
            # OpenMX chunk whose conversion/import had not yet been committed.
            for batch in pending:
                check_stop(stop)
                legacy = folder / f'chunk_{batch:05d}'
                if not (legacy / 'openmx.json').is_file():
                    continue
                start = batch * CHUNK_POINTS; end = min(count, start + CHUNK_POINTS)
                points = grid_points(ns, start, end)
                if completed_openmx(legacy, legacy.name, input_text(legacy.name, points)):
                    if log: log(f'[Reusing old OpenMX output] k points {start+1:,}–{end:,}')
                    completed_before = sum(min(CHUNK_POINTS, count - i * CHUNK_POINTS)
                                           for i in np.flatnonzero(f['Mesh/CompletedChunks'][:]))
                    def conversion_progress(current, total):
                        if conversion:
                            conversion(completed_before, end - start, current, total, count)
                    raw_path, seconds, conversion_seconds = calculate_chunk(
                        legacy, points, log, stop, conversion_progress)
                    elapsed += seconds
                    import_started = time.monotonic()
                    store_raw(raw_path, np.array([batch]), points)
                    conversion_seconds += time.monotonic() - import_started
                    if timing:
                        timing(seconds, conversion_seconds)
            pending = np.flatnonzero(~f['Mesh/CompletedChunks'][:])
            if len(pending):
                pieces = []
                for batch in pending:
                    check_stop(stop)
                    pieces.append(grid_points(ns, batch*CHUNK_POINTS, min(count, (batch+1)*CHUNK_POINTS)))
                points = np.concatenate(pieces)
                del pieces
                # Stable job name also survives conversion failures and stopped
                # HDF5 imports. The material/mesh identity remains unchanged.
                job_code = request_code({'checkpoints': pending.tolist(), 'size': CHUNK_POINTS})
                job = folder / f'batch_{job_code}'
                if log: log(f'[First-BZ batch] Unsaved {len(points):,} / {count:,} k points (at most one OpenMX run)')
                completed_before = sum(min(CHUNK_POINTS, count - i * CHUNK_POINTS)
                                       for i in np.flatnonzero(f['Mesh/CompletedChunks'][:]))
                def conversion_progress(current, total):
                    if conversion:
                        conversion(completed_before, len(points), current, total, count)
                raw_path, seconds, conversion_seconds = calculate_chunk(
                    job, points, log, stop, conversion_progress)
                elapsed += seconds
                import_started = time.monotonic()
                store_raw(raw_path, pending, points)
                conversion_seconds += time.monotonic() - import_started
                if timing:
                    timing(seconds, conversion_seconds)
            f.attrs['Complete'] = True; f.flush()
        if not cache_ok(partial, request):
            raise RuntimeError('First-BZ cache validation failed.')
        partial.replace(final)
        # The final HDF5 and the structure-level SCF restart are sufficient for
        # reuse. Remove the large OpenMX work directories after committing them.
        for chunk in itertools.chain(folder.glob('chunk_*'), folder.glob('batch_*')):
            if chunk.is_dir():
                shutil.rmtree(chunk)
        return final, elapsed


class EigenstateGrid:
    """Periodic nearest-neighbour reference without mixing eigenvector gauges."""
    def __init__(self, path):
        self.path = Path(path)
        self.file = h5py.File(path, 'r')
        self.points = self.file['Mesh/FirstBZCoordinates'][:]
        self.ns = np.asarray(self.file["Mesh"].attrs["Divisions"], dtype=int)
        self.paths = lcao_paths(self.file)
        self.ef = float(self.file['Output'].attrs['EF_Eh'])

    def close(self):
        self.file.close()

    def lookup(self, cartesian):
        target = np.asarray(cartesian).reshape(-1, 3)
        # G1=(gx,gy,0), G2=(0,gy,0), G3=(0,0,gz). Solve the
        # Cartesian closest-point problem on the infinite grid, then wrap the
        # integer indices. This remains correct for very unequal N1/N2 values;
        # rounding all three fractional coordinates independently does not.
        step = G / self.ns[:, None]
        center = np.rint(target[:, 0] / step[0, 0]).astype(np.int64)
        iz = np.rint(target[:, 2] / step[2, 2]).astype(np.int64)
        def candidate(ix):
            iy = np.rint((target[:, 1] - ix * step[0, 1]) / step[1, 1]).astype(np.int64)
            index = np.stack((ix, iy, iz), axis=-1)
            distance2 = np.sum((index @ step - target)**2, axis=1)
            return index, distance2
        best_index, best = candidate(center)
        radius = int(np.ceil(np.sqrt(best.max(initial=0)) / abs(step[0, 0])))
        for shift in range(-radius, radius + 1):
            index, distance2 = candidate(center + shift)
            improve = distance2 < best - 1e-14
            best_index[improve] = index[improve]; best[improve] = distance2[improve]
        indices = np.ravel_multi_index((best_index % self.ns).T, tuple(self.ns))
        return indices, np.sqrt(best)

    def take(self, path, indices):
        unique, inverse = np.unique(indices, return_inverse=True)
        return self.file[path][unique][inverse]

    def bands(self, cartesian):
        indices, error = self.lookup(cartesian)
        return (self.take('Output/Band', indices) - self.ef) * 27.211386245988, error

    def interpolated_bands(self, cartesian):
        """Periodically trilinear-interpolate scalar eigenvalues on the mesh.

        This is only for band visualization. LCAO coefficients retain nearest
        point sampling because their independently chosen gauges cannot be mixed.
        """
        target = np.asarray(cartesian, dtype=float).reshape(-1, 3)
        scaled = (target @ np.linalg.inv(G)) * self.ns
        lower = np.floor(scaled).astype(np.int64)
        fraction = scaled - lower
        values = None
        for bits in itertools.product((0, 1), repeat=3):
            offset = np.asarray(bits, dtype=np.int64)
            index3 = (lower + offset) % self.ns
            indices = np.ravel_multi_index(index3.T, tuple(self.ns))
            weight = np.prod(np.where(offset, fraction, 1.0-fraction), axis=1)
            corner = self.take('Output/Band', indices)
            if values is None:
                values = np.zeros_like(corner, dtype=float)
            values += weight[:, None] * corner
        _, error = self.lookup(target)
        return (values - self.ef) * 27.211386245988, error

    def write_path(self, destination, cartesian, attrs):
        """SPADExp adapter: requested final-state k, paired cached initial state.

        OpenMX coefficients are reciprocal-periodic in this export convention;
        do not apply an extra exp(i G.tau). Keep requested coordinates for PAD's
        atomic phase and store the actual sampled coordinate/error separately.
        """
        destination = Path(destination); indices, errors = self.lookup(cartesian)
        temporary = Path(str(destination) + '.tmp')
        try:
            with h5py.File(temporary, 'w') as out:
                out.attrs.update(self.file.attrs)
                inp = out.create_group('Input')
                for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
                    self.file.copy('Input/' + name, inp, name=name)
                k = inp.create_group('Kpath'); k.attrs.update(attrs)
                k.attrs.update(Count=len(indices), minN=BANDS[0], maxN=BANDS[1], FirstBZFolded=np.uint8(1),
                               UniqueCount=len(np.unique(indices)), StateSampling='periodic nearest neighbour')
                k.create_dataset('Coordinates', data=np.asarray(cartesian) @ np.linalg.inv(G))
                k.create_dataset('SampledFirstBZCoordinates', data=self.points[indices])
                k.create_dataset('MeshIndex', data=indices)
                k.create_dataset('SamplingDistance_Ainv', data=errors)
                group = out.create_group('Output'); group.attrs.update(self.file['Output'].attrs)
                self.file.copy('Output/Potential', group, name='Potential')
                group.create_dataset('Band', data=self.take('Output/Band', indices))
                for path in self.paths:
                    parent, name = path.rsplit('/', 1)
                    out.require_group(parent).create_dataset(name, data=self.take(path, indices),
                                                             compression='gzip', shuffle=True)
            temporary.replace(destination)
        finally:
            temporary.unlink(missing_ok=True)
        return float(errors.max(initial=0))
