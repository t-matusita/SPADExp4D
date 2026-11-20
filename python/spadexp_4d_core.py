"""First-BZ eigenstate database, independent of the 3-D viewer and its caches."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shlex
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


def set_structures_root(path):
    """Use another data folder for structure inputs and every OpenMX output.

    The selected Structure is cleared; call configure_structure() for a folder
    inside the new root before running OpenMX again.
    """
    global STRUCTURES, DEFAULT_STRUCTURE, CACHE, TEMPLATE, RESTART
    global STRUCTURE_RESULTS, STRUCTURE_CACHE
    STRUCTURES = Path(path).expanduser().resolve()
    DEFAULT_STRUCTURE = STRUCTURES / 'graphite'
    CACHE = DEFAULT_STRUCTURE / 'cache'
    TEMPLATE = DEFAULT_STRUCTURE / '.openmx.base.dat'
    RESTART = CACHE / 'scf_rst'
    STRUCTURE_RESULTS = None
    STRUCTURE_CACHE = None
    return STRUCTURES


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


ELEMENT_SYMBOLS = frozenset("""
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br
    Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho
    Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es
    Fm Md No Lr""".split())


def cif_atom_sites(text):
    """Rows of the CIF ``_atom_site_*`` loop as dicts keyed by tag (lower case)."""
    lines = text.splitlines(); sites = []
    i = 0
    while i < len(lines):
        if lines[i].strip().lower() != 'loop_':
            i += 1; continue
        i += 1; tags = []
        while i < len(lines) and lines[i].strip().startswith('_'):
            tags.append(lines[i].split()[0].lower()); i += 1
        if not any(tag in ('_atom_site_label', '_atom_site_type_symbol') for tag in tags):
            continue
        values = []
        while i < len(lines):
            stripped = lines[i].strip()
            if not stripped or stripped.startswith('#'):
                if values: break
                i += 1; continue
            if stripped.lower() == 'loop_' or stripped.startswith('_') or stripped.lower().startswith('data_'):
                break
            try:
                values.extend(shlex.split(stripped))
            except ValueError:
                values.extend(stripped.split())
            i += 1
        for start in range(0, len(values) - len(tags) + 1, len(tags)):
            sites.append(dict(zip(tags, values[start:start + len(tags)])))
    return sites


CELL_TAGS = ('_cell_length_a', '_cell_length_b', '_cell_length_c',
             '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma')
CIF_FIX_ADVICE = ('Correct the CIF in a text editor (or export it again, for example from VESTA) '
                  'and import it again.')


class CifError(ValueError):
    """A CIF that cannot be used, with every problem explained in plain English.

    ``problems`` holds one sentence per problem (what is wrong and how to fix
    it); ``details`` holds technical output such as the cif2omx messages.
    """
    def __init__(self, path, problems, details='', advice=CIF_FIX_ADVICE):
        self.path = Path(path); self.problems = list(problems); self.details = details
        self.advice = advice
        super().__init__(self.message())

    def message(self):
        return '\n'.join([f'The CIF file "{self.path.name}" cannot be used:']
                         + [f'• {problem}' for problem in self.problems] + [self.advice])


def cif_number(value):
    """A CIF number such as 4.196(2); None when it is not a number."""
    try:
        return float(re.sub(r'\([^)]*\)', '', str(value)))
    except ValueError:
        return None


def cif_cell_values(text):
    """Cell lengths and angles of a CIF, and the problems found with them."""
    values, problems = {}, []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0].lower() in CELL_TAGS:
            values[fields[0].lower()] = fields[1]
    missing = [tag for tag in CELL_TAGS if tag not in values]
    if missing:
        problems.append('The unit cell is incomplete: ' + ', '.join(missing) + ' '
                        + ('is' if len(missing) == 1 else 'are') + ' missing. A CIF needs all six '
                        '_cell_length_a/b/c (Å) and _cell_angle_alpha/beta/gamma (degrees).')
    cell = {}
    for tag in CELL_TAGS:
        if tag not in values:
            continue
        number = cif_number(values[tag])
        if number is None:
            problems.append(f'{tag} is "{values[tag]}", which is not a number.')
        elif tag.startswith('_cell_length') and number <= 0:
            problems.append(f'{tag} is {values[tag]}; a cell length must be larger than 0 Å.')
        elif tag.startswith('_cell_angle') and not 0 < number < 180:
            problems.append(f'{tag} is {values[tag]}; a cell angle must be between 0 and 180 degrees.')
        else:
            cell[tag] = number
    return cell, problems


def cif_problems(path):
    """Everything in a CIF that the automatic OpenMX conversion cannot use.

    Returns plain-English sentences (what is wrong and how to fix it): an
    incomplete or invalid unit cell, no atom-site list or no fractional
    coordinates, type symbols that are not chemical elements (after removing
    oxidation states such as Ta4+), partly occupied sites and several atoms at
    one position.  An empty list means the CIF looks usable.
    """
    try:
        text = normalized_cif_text(path)
    except OSError as exc:
        return [f'The file could not be read ({exc.strerror or exc}).']
    if not text.strip():
        return ['The file is empty.']
    _, problems = cif_cell_values(text)
    sites = cif_atom_sites(text)
    if not sites:
        problems.append('No atoms were found. A CIF needs a loop_ that lists _atom_site_label, '
                        '_atom_site_type_symbol and the fractional coordinates _atom_site_fract_x, '
                        '_atom_site_fract_y and _atom_site_fract_z, one line per atom.')
        return problems
    axes = ('_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z')
    missing = [tag for tag in axes if tag not in sites[0]]
    if missing:
        problems.append('The atom list has no ' + ', '.join(missing) + '. Atom positions must be given as '
                        'fractional coordinates (_atom_site_fract_x, _atom_site_fract_y, _atom_site_fract_z).')
    unknown, partial, positions = [], [], {}
    for site in sites:
        label = site.get('_atom_site_label', '?')
        symbol = site.get('_atom_site_type_symbol') or re.match(r'[A-Za-z]*', label).group(0)
        if symbol not in ELEMENT_SYMBOLS:
            hint = symbol[:1].upper() + symbol[1:].lower()
            unknown.append(f'"{symbol}" (write "{hint}")' if hint in ELEMENT_SYMBOLS
                           else f'"{symbol}"' + ('' if symbol == label else f' at site {label}'))
        occupancy = cif_number(site.get('_atom_site_occupancy', '1'))
        if occupancy is not None and occupancy < 0.999:
            partial.append(f'{label} ({occupancy:g})')
        if not missing:
            xyz = [cif_number(site[tag]) for tag in axes]
            if any(value is None for value in xyz):
                problems.append(f'Atom {label} has a coordinate that is not a number: '
                                + ', '.join(f'{tag} = "{site[tag]}"' for tag, value in zip(axes, xyz)
                                            if value is None) + '.')
            else:
                key = tuple(round(value % 1.0, 4) % 1.0 for value in xyz)
                positions.setdefault(key, []).append(label)
    if unknown:
        problems.append('These _atom_site_type_symbol values are not chemical elements: ' + ', '.join(unknown)
                        + '. Each atom needs a real element symbol such as Fe, O or Si; placeholders '
                        'such as T, M or X from a structure template must be replaced.')
    if partial:
        problems.append('These sites are only partly occupied (_atom_site_occupancy): ' + ', '.join(partial)
                        + '. Mixed or partly occupied sites are not supported; keep one fully occupied '
                        'atom per site.')
    for position, labels in positions.items():
        if len(labels) > 1:
            problems.append(f'Atoms {", ".join(labels)} are at the same position '
                            f'({", ".join(f"{value:g}" for value in position)}). Keep only one atom there.')
    return problems


def undefined_template_species(path):
    """Species of an OpenMX template without a PAO/VPS (cif2omx writes 'undefined')."""
    try:
        text = Path(path).read_text(errors='replace')
    except OSError:
        return []
    match = re.search(r'<Definition\.of\.Atomic\.Species(.*?)Definition\.of\.Atomic\.Species>', text, flags=re.S)
    if match is None:
        return []
    return [line.split()[0] for line in match.group(1).splitlines()
            if line.split() and re.search(r'\bundefined\b', line, flags=re.I)]


def template_has_undefined_species(path):
    try:
        text = Path(path).read_text(errors='replace')
    except OSError:
        return True
    match = re.search(r'<Definition\.of\.Atomic\.Species(.*?)Definition\.of\.Atomic\.Species>',
                      text, flags=re.S)
    return match is None or bool(re.search(r'\bundefined\b', match.group(1), flags=re.I))


def complete_band_range(path, spin=None):
    """Return every zero-based band represented by an OpenMX basis.

    OpenMX has one band per atomic basis function for scf.SpinPolarization off,
    the same count per spin for on, and two per basis function (spinors) for
    NC.  `spin` overrides the template's setting.  No occupation estimate is
    used, so unoccupied states above EF remain available after an EF shift.
    """
    text = Path(path).read_text(errors='replace')
    if spin is None:
        spin_match = re.search(r'(?mi)^\s*scf\.SpinPolarization\s+(\S+)', text)
        spin = spin_match.group(1) if spin_match else 'off'
    spin = str(spin).lower()
    if spin not in ('off', 'on', 'nc'):
        raise ValueError(f'Unsupported scf.SpinPolarization: {spin}')
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
    if spin == 'nc':
        count *= 2
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
        raise CifError(path, cif_cell_values(Path(path).read_text(errors='replace'))[1]
                       or ['The unit cell could not be read.'])
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
        problems = cif_problems(cif)
        if problems:
            raise CifError(cif, problems)
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
        output = (completed.stdout + completed.stderr).strip()[-4000:]
        if completed.returncode != 0 or not template.is_file():
            raise CifError(cif, ['cif2omx could not turn the CIF into an OpenMX input. '
                                 'Its messages are shown in the details.'], output)
        if template_has_undefined_species(template):
            species = undefined_template_species(template)
            template.unlink(missing_ok=True)
            names = ', '.join(f'"{name}"' for name in species) or 'some atoms'
            if species and all(name in ELEMENT_SYMBOLS for name in species):
                # Real elements that the OpenMX library does not cover.
                raise CifError(cif, [f'OpenMX has no basis set (PAO/VPS in DFT_DATA19) for {names}, '
                                     'so this structure cannot be calculated.'], output,
                               advice='Use a structure without these elements, or add PAO/VPS files '
                                      'for them to DFT_DATA19.')
            raise CifError(cif, [f'OpenMX has no basis set (PAO/VPS in DFT_DATA19) for {names}. '
                                 'Check that _atom_site_type_symbol names the element of every atom.'], output)
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


# OpenMX keywords editable in the GUI (OpenMX Parameters tab).  Only values
# that differ from the cif2omx template are stored, as OpenMX text under
# "openmx" in the structure's parameter.json, and they are applied on top of
# the template when the OpenMX input is written.  Without overrides the
# input, the SCF restart and every cache identity are exactly as before.
OPENMX_KEYWORDS = ('scf.XcType', 'scf.SpinPolarization', 'scf.SpinOrbit.Coupling',
                   'scf.energycutoff', 'scf.Kgrid', 'scf.ElectronicTemperature',
                   'scf.criterion', 'scf.maxIter')
XC_TYPES = ('LDA', 'LSDA-CA', 'LSDA-PW', 'GGA-PBE')
SPIN_MODES = ('off', 'on', 'NC')  # NC: non-collinear spinors (needed for spin-orbit coupling)
# OpenMX's own defaults (Input_std.c), used when a template omits a keyword.
OPENMX_DEFAULTS = {'scf.XcType': 'LDA', 'scf.SpinPolarization': 'off',
                   'scf.SpinOrbit.Coupling': 'off', 'scf.energycutoff': '220',
                   'scf.Kgrid': '4 4 4', 'scf.ElectronicTemperature': '300',
                   'scf.criterion': '1e-06', 'scf.maxIter': '40'}


def normalize_openmx_value(key, value):
    """Return the canonical OpenMX text of one editable keyword value."""
    if key not in OPENMX_KEYWORDS:
        raise ValueError(f'Unsupported OpenMX keyword: {key}')
    text = ' '.join(str(value).split())
    try:
        if key in ('scf.XcType', 'scf.SpinPolarization', 'scf.SpinOrbit.Coupling'):
            choices = {'scf.XcType': XC_TYPES, 'scf.SpinPolarization': SPIN_MODES,
                       'scf.SpinOrbit.Coupling': ('off', 'on')}[key]
            match = {name.lower(): name for name in choices}.get(text.lower())
            if match:
                return match
        elif key == 'scf.Kgrid':
            numbers = [float(v) for v in text.split()]
            if len(numbers) == 3 and all(v.is_integer() and 1 <= v <= 1000 for v in numbers):
                return ' '.join(str(int(v)) for v in numbers)
        elif key == 'scf.maxIter':
            number = float(text)
            if number.is_integer() and number >= 1:
                return str(int(number))
        else:
            number = float(text)
            if np.isfinite(number) and number > 0:
                return f'{number:.10g}'
    except ValueError:
        pass
    raise ValueError(f'Invalid value for {key}: {value!r}')


def validate_openmx_overrides(values):
    """Return the canonical {keyword: text} overrides, in keyword order."""
    if values is None:
        return {}
    if not isinstance(values, dict):
        raise ValueError('"openmx" must be a JSON object of OpenMX keywords.')
    unknown = sorted(set(values) - set(OPENMX_KEYWORDS))
    if unknown:
        raise ValueError('Unsupported OpenMX keywords: ' + ', '.join(unknown))
    return {key: normalize_openmx_value(key, values[key]) for key in OPENMX_KEYWORDS if key in values}


def openmx_overrides(settings):
    return validate_openmx_overrides((settings or {}).get('openmx'))


def effective_openmx(overrides=None, path=None):
    """Template keywords with overrides applied; checks spin/SOC consistency."""
    values = {}
    for key, text in template_keywords(path).items():
        try:
            values[key] = normalize_openmx_value(key, text if text is not None else OPENMX_DEFAULTS[key])
        except ValueError:
            values[key] = ' '.join(str(text).split())
    values.update(validate_openmx_overrides(overrides))
    if values['scf.SpinOrbit.Coupling'] == 'on' and values['scf.SpinPolarization'] != 'NC':
        raise ValueError('scf.SpinOrbit.Coupling on requires scf.SpinPolarization NC.')
    return values


def request_band_range(overrides=None):
    """All zero-based bands for the effective spin setting.

    Unless the spin mode is overridden this is BANDS, the range configured for
    the selected structure, so existing cache identities are unchanged.
    """
    spin = effective_openmx(overrides)['scf.SpinPolarization']
    if spin == effective_openmx({})['scf.SpinPolarization']:
        return BANDS
    return complete_band_range(TEMPLATE, spin)


def validate_band_window(window):
    """Return None (store every band) or (low, high) in eV relative to EF."""
    if window is None:
        return None
    try:
        low, high = (float(v) for v in window)
    except (TypeError, ValueError) as exc:
        raise ValueError('band_window must be null or [low, high] in eV relative to EF.') from exc
    if not (np.isfinite(low) and np.isfinite(high) and low < high):
        raise ValueError('band_window needs finite low < high (eV relative to EF).')
    return (low, high)


def band_window(settings):
    return validate_band_window((settings or {}).get('band_window'))


def _keyword_line(key):
    return re.compile(rf'(?mi)^[ \t]*{re.escape(key)}[ \t]+([^#\n]*).*$')


def template_keywords(path=None):
    """Values of the editable keywords in an OpenMX template (None when absent)."""
    text = Path(path or TEMPLATE).read_text(errors='replace')
    values = {}
    for key in OPENMX_KEYWORDS:
        match = _keyword_line(key).search(text)
        values[key] = ' '.join(match.group(1).split()) if match else None
    return values


RECOMMENDED_KGRID_SPACING = 0.2  # 1/A, SCF k-point spacing of the recommended scf.Kgrid


def template_cell(path=None):
    """Lattice vectors (rows, Angstrom) of an OpenMX template."""
    text = Path(path or TEMPLATE).read_text(errors='replace')
    block = re.search(r'<Atoms\.UnitVectors[ \t]*\n(.*?)\n[ \t]*Atoms\.UnitVectors>', text, re.S)
    if block is None:
        raise ValueError('Atoms.UnitVectors block not found in the OpenMX template.')
    cell = np.array([[float(v) for v in line.split()[:3]]
                     for line in block.group(1).strip().splitlines() if line.strip()])
    if cell.shape != (3, 3):
        raise ValueError('Atoms.UnitVectors must contain three vectors.')
    unit = re.search(r'(?mi)^[ \t]*Atoms\.UnitVectors\.Unit[ \t]+(\S+)', text)
    return cell * 0.529177 if unit and unit.group(1).lower() == 'au' else cell


def recommended_openmx_values(path=None):
    """Recommended SCF keywords for a template's cell, never weaker than the template.

    scf.Kgrid uses a spacing of at most RECOMMENDED_KGRID_SPACING; the cutoff is
    at least 220 Ry, the criterion at most 1e-6 Hartree and maxIter at least
    200.  scf.XcType stays as in the template, which matches its pseudopotentials.
    """
    template = template_keywords(path)

    def number(key, fallback):
        try:
            return float(normalize_openmx_value(key, template[key]))
        except ValueError:
            return fallback

    reciprocal = 2 * np.pi * np.linalg.inv(template_cell(path)).T
    grid = [max(1, int(np.ceil(np.linalg.norm(g) / RECOMMENDED_KGRID_SPACING - 1e-9))) for g in reciprocal]
    try:
        current = [int(v) for v in normalize_openmx_value('scf.Kgrid', template['scf.Kgrid']).split()]
        grid = [max(new, old) for new, old in zip(grid, current)]
    except ValueError:
        pass
    try:
        xc = normalize_openmx_value('scf.XcType', template['scf.XcType'])
    except ValueError:
        xc = 'GGA-PBE'

    def kept(key):  # spin settings define the physics model: keep the template's
        try:
            return normalize_openmx_value(key, template[key])
        except ValueError:
            return OPENMX_DEFAULTS[key]

    values = {'scf.XcType': xc,
              'scf.SpinPolarization': kept('scf.SpinPolarization'),
              'scf.SpinOrbit.Coupling': kept('scf.SpinOrbit.Coupling'),
              'scf.energycutoff': max(number('scf.energycutoff', 220.0), 220.0),
              'scf.Kgrid': ' '.join(map(str, grid)),
              'scf.ElectronicTemperature': 300,
              'scf.criterion': min(number('scf.criterion', 1e-6), 1e-6),
              'scf.maxIter': max(number('scf.maxIter', 200), 200)}
    return {key: normalize_openmx_value(key, values[key]) for key in OPENMX_KEYWORDS}


def complete_noncollinear_atoms(text):
    """Give every atom line the explicit NC columns OpenMX reads.

    For scf.SpinPolarization NC OpenMX sscanf()s spin/orbital angles, a
    constraint flag and orbital polarization after the up/down populations;
    absent fields would keep undefined values.  cif2omx writes only
    "id species x y z up down", so append "0 0 0 0 0 off" (no constraint).
    """
    match = re.search(r'(<Atoms\.SpeciesAndCoordinates[^\n]*\n)(.*?)(\n[ \t]*Atoms\.SpeciesAndCoordinates>)',
                      text, re.S)
    if match is None:
        return text
    lines = []
    for line in match.group(2).split('\n'):
        fields = line.split('#', 1)[0].split()
        lines.append(f'{line.rstrip()}  0.0 0.0 0.0 0.0 0 off'
                     if len(fields) == 7 and fields[0].isdigit() else line)
    return text[:match.start(2)] + '\n'.join(lines) + text[match.end(2):]


def apply_openmx_overrides(text, overrides):
    """Replace (or append) the overridden keyword lines of an OpenMX input."""
    for key, value in (overrides or {}).items():
        line = f'{key} {value}'
        pattern = _keyword_line(key)
        if pattern.search(text):
            text = pattern.sub(lambda _match: line, text, count=1)
        else:
            text = text.rstrip('\n') + f'\n{line}\n'
    return text


def scf_restart_folder(overrides=None):
    """SCF restart that belongs to the template plus these keyword overrides.

    Without overrides this is the structure's existing restart.  Overrides
    (another cutoff, functional, ...) get their own folder, so a restart made
    with different SCF settings is never reused through scf.restart pp.
    """
    if not overrides:
        return RESTART
    identity = hashlib.sha256((digest(TEMPLATE) + json.dumps(overrides, sort_keys=True)).encode())
    return RESTART.with_name(f'{RESTART.name}_{identity.hexdigest()[:16]}')


def mesh_request(settings):
    # No photoemission / plotting parameter enters this identity. These fixed
    # physical fingerprints only guard against reusing a different material.
    if not OPENMX.is_file():
        raise FileNotFoundError(
            'Please set the path to the OpenMX executable. '
            'Use SPADEXP4D_OPENMX or add openmx to PATH.')
    overrides = openmx_overrides(settings)
    window = band_window(settings)
    restart = scf_restart_folder(overrides)
    request = {'schema': SCHEMA, 'divisions': list(divisions(settings)), 'lattice': LATTICE.tolist(),
               'bands': list(request_band_range(overrides)), 'template': digest(TEMPLATE),
               'data': {str(p.relative_to(DATA)): digest(p) for p in sorted(DATA.rglob('*')) if p.is_file()},
               'restart': (RESTART_FINGERPRINT if RESTART_FINGERPRINT is not None else
                           ([stamp(p) for p in sorted(restart.iterdir()) if p.is_file()]
                            if restart.is_dir() else [])),
               'openmx': stamp(OPENMX), 'scf_restart': 'pp'}
    if overrides:  # absent when unchanged, so existing cache identities stay valid
        request['openmx_keywords'] = overrides
    if window is not None:  # likewise absent when every band is stored
        request['band_window'] = list(window)
    return request


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


MO_CUBE_KEYS = ('num.HOMOs', 'num.LUMOs')


def input_text(stem, points, openmx=None):
    overrides = validate_openmx_overrides(openmx)
    bands = request_band_range(overrides)
    text = apply_openmx_overrides(TEMPLATE.read_text(), overrides)
    if effective_openmx(overrides)['scf.SpinPolarization'] == 'NC':
        text = complete_noncollinear_atoms(text)
    text = re.sub(r'(?m)^System\.Name\s+\S+', f'System.Name {stem}', text)
    text = re.sub(r'(?m)^DATA\.PATH\s+.*$', f'DATA.PATH {DATA}', text)
    # OpenMX switches to the Cluster solver for a Gamma-only SCF grid.  That
    # solver rejects the multiple MO.kpoint entries used by the mesh batches,
    # so keep a non-Gamma SCF grid when the CIF converter emitted 1 1 1.
    if len(points) > 1:
        text = re.sub(r'(?m)^scf\.Kgrid\s+1\s+1\s+1\s*$',
                      'scf.Kgrid 2 2 1', text)
    restart_file = scf_restart_folder(overrides) / f'{RESTART_NAME}.crst_check'
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
             f'SPADExp.minN {bands[0]}\nSPADExp.maxN {bands[1]}')
    if re.search(r'SPADExp\.dimension.*?SPADExp\.maxN\s+\d+', text, flags=re.S):
        text = re.sub(r'SPADExp\.dimension.*?SPADExp\.maxN\s+\d+', block, text, count=1, flags=re.S)
    else:
        text += '\n' + block + '\n'
    # With MO.fileout on, OpenMX also writes num.HOMOs/num.LUMOs cube files
    # (default 1 each, real and imaginary) for every MO k point: ~1.2 MB per
    # point for graphite, tens of GB for a fine mesh.  The converter reads only
    # the .EV LCAO block, which always lists every band, so disable the cubes.
    cube_lines = ''
    for key in MO_CUBE_KEYS:
        pattern = rf'(?mi)^{re.escape(key)}\s+.*$'
        if re.search(pattern, text):
            text = re.sub(pattern, f'{key} 0', text)
        else:
            cube_lines += f'{key} 0\n'
    return text + (f'\n{cube_lines}MO.fileout on\nMO.Nkpoint {len(points)}\n<MO.kpoint\n'
                   + '\n'.join(' '.join(f'{v:.12f}' for v in p) for p in points) + '\nMO.kpoint>\n')


def lcao_paths(file):
    paths = []
    file['Output/LCAO'].visititems(lambda name, obj: paths.append('Output/LCAO/' + name)
                                 if isinstance(obj, h5py.Dataset) else None)
    return sorted(paths)


def spin_of(file):
    """Spin mode of a converted file: 'off', 'on' or 'nc'."""
    value = file['Output'].attrs['Spin']
    return (value.decode() if isinstance(value, (bytes, np.bytes_)) else str(value)).lower()


def band_paths(spin):
    """Eigenvalue datasets: up and down sets for collinear spin, else one set."""
    return ['Output/BandUp', 'Output/BandDn'] if str(spin).lower() == 'on' else ['Output/Band']


def stored_band_range(file):
    output = file['Output'].attrs
    return (int(output['BandIndexStartZeroBased']), int(output['BandIndexEndZeroBased']))


def expected_lcao_paths(file):
    """Build the expected radial-orbital dataset list from OpenMX input data."""
    def text(value):
        return value.decode() if isinstance(value, (bytes, np.bytes_)) else str(value)

    species_labels = [text(value) for value in file['Input/Atomic.Species/Labels'][:]]
    species_orbits = [text(value) for value in file['Input/Atomic.Species/Orbits'][:]]
    orbit_by_species = dict(zip(species_labels, species_orbits))
    atom_labels = [text(value) for value in file['Input/Atoms.SpeciesAndCoordinates/Labels'][:]]
    suffixes = ('Up', 'Dn') if spin_of(file) == 'on' else ('',)
    expected = []
    for atom_index, label in enumerate(atom_labels, 1):
        orbit_text = orbit_by_species[label]
        orbit_groups = re.findall(r'([spdf])(\d+)', orbit_text)
        if not orbit_groups or ''.join(kind + count for kind, count in orbit_groups) != orbit_text:
            raise ValueError(f'Unrecognized OpenMX basis: {label} {orbit_text}')
        for kind, count in orbit_groups:
            for radial_index in range(int(count)):
                for suffix in suffixes:
                    expected.append(f'Output/LCAO/{atom_index}_{label}/{kind}{radial_index}{suffix}')
    return sorted(expected)


def states_ok(file, count, compressed=False, band_range=None, spin=None):
    """Check a converted or cached eigenstate file.

    band_range defaults to BANDS; pass False to accept the stored range (band
    window).  spin is the expected mode ('off', 'on', 'nc'; default: off,
    as before spin support).
    """
    try:
        if band_range is None:
            band_range = BANDS
        expected_spin = 'off' if spin is None else str(spin).lower()
        start, end = stored_band_range(file)
        band_count = end - start + 1
        mode = spin_of(file)
        digit = 4 if mode == 'nc' else 2
        paths = lcao_paths(file)
        return (mode == expected_spin and band_count >= 1
                and (band_range is False or (start, end) == tuple(band_range))
                and all(file[path].shape == (count, band_count) for path in band_paths(mode))
                and paths == expected_lcao_paths(file)
                and all(file[p].shape[:2] == (count, band_count) and file[p].shape[-1] == digit
                        and (not compressed or (file[p].dtype == np.float32 and file[p].compression == 'gzip'))
                        for p in paths)
                and 'Output/Potential/V0' in file)
    except (KeyError, ValueError, OSError):
        return False


def request_spin(request):
    """Spin mode ('off', 'on', 'nc') that a mesh request stores."""
    return effective_openmx(request.get('openmx_keywords'))['scf.SpinPolarization'].lower()


def request_states_ok(file, request, count, compressed=False):
    window = request.get('band_window')
    bands = tuple(request['bands']) if 'bands' in request else None  # None: BANDS
    return states_ok(file, count, compressed, False if window else bands, request_spin(request))


def cache_ok(path, request):
    try:
        with h5py.File(path, 'r') as f:
            count = int(np.prod(request['divisions']))
            return (bool(f.attrs.get('Complete', False))
                    and f.attrs.get('Request') == json.dumps(request, sort_keys=True)
                    and f['Mesh/FirstBZCoordinates'].shape == (count, 3)
                    and np.all(f['Mesh/CompletedChunks'][:]) and request_states_ok(f, request, count, True))
    except (OSError, ValueError, KeyError):
        return False


MESH_TIMING_ATTR = 'CalculationTiming'


def mesh_timing(path):
    """Calculation record stored in a k-space HDF5 file, or None.

    Returns a dict with ``openmx_seconds``, ``conversion_seconds``,
    ``wall_seconds`` and ``completed_at`` (local time) for caches written by
    ensure_mesh(); older caches only give ``created`` from the converter's
    Datetime attribute.  Reading never modifies the file (its stamp keys the
    PAD caches).
    """
    def text(value):
        return value.decode() if isinstance(value, bytes) else str(value)
    try:
        with h5py.File(path, 'r') as f:
            stored = f.attrs.get(MESH_TIMING_ATTR); created = f.attrs.get('Datetime')
    except (OSError, ValueError):
        return None
    record = {}
    if stored is not None:
        try:
            value = json.loads(text(stored))
        except ValueError:
            value = None
        if isinstance(value, dict):
            record.update(value)
    if created is not None:
        record.setdefault('created', text(created))
    return record or None


def mesh_plan(settings):
    request = mesh_request(settings)
    folder = CACHE / ('mesh_' + request_code(request))
    path = folder / 'eigenstates.hdf5'
    return {'request': request, 'folder': folder, 'path': path, 'reuse': cache_ok(path, request),
            'count': int(np.prod(request['divisions']))}


def saved_mesh_points(plan):
    """k points of a planned mesh already saved by an interrupted calculation.

    Counts the committed batches of the partial HDF5 of the same request (the
    ones ensure_mesh() keeps); 0 when there is none.  Never modifies files.
    """
    if plan.get('reuse'):
        return plan['count']
    count = plan['count']; batches = (count + CHUNK_POINTS - 1) // CHUNK_POINTS
    try:
        with h5py.File(plan['folder'] / 'eigenstates.partial.hdf5', 'r') as f:
            if (f.attrs.get('Request') != json.dumps(plan['request'], sort_keys=True)
                    or f['Mesh/CompletedChunks'].shape != (batches,)):
                return 0
            done = np.flatnonzero(f['Mesh/CompletedChunks'][:])
    except (OSError, KeyError, ValueError):
        return 0
    return int(sum(min(CHUNK_POINTS, count - i * CHUNK_POINTS) for i in done))


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
        # of the next input, and num.HOMOs/num.LUMOs only select cube files.
        # Reuse the already completed MO output when the original recorded
        # input differs solely in those lines.
        dat = folder / f'{stem}.dat'
        if not input_matches and record and dat.is_file():
            stored = dat.read_text(errors='replace')
            stored_recorded = record.get('input') == hashlib.sha256(stored.encode()).hexdigest()
            def without_restart(value):
                return re.sub(r'(?mi)^(?:scf\.restart(?:\.filename)?|orbitalOpt\.Force\.Skip'
                              r'|num\.HOMOs|num\.LUMOs)\s+.*\n?', '', value)
            input_matches = stored_recorded and without_restart(stored) == without_restart(text)
        return normal and output_matches and input_matches
    except OSError:
        return False


def preserve_structure_restart(folder, stem, restart=None):
    """Keep one reusable SCF restart (per keyword set) beside the structure."""
    restart = Path(restart) if restart is not None else RESTART
    if STRUCTURE_CACHE is None or restart.parent != STRUCTURE_CACHE:
        return
    source = Path(folder) / f'{stem}_rst'
    if not source.is_dir() or restart.is_dir():
        return
    temporary = restart.with_name(restart.name + '.tmp')
    if temporary.exists():
        raise RuntimeError(f'A temporary resume folder remains: {temporary}')
    shutil.copytree(source, temporary)
    for path in sorted(temporary.iterdir()):
        if path.name.startswith(stem):
            path.rename(path.with_name(RESTART_NAME + path.name[len(stem):]))
    temporary.replace(restart)


EH_EV = 27.211386245988
BAND_WINDOW_MARGIN_EV = 1.0  # extra range when choosing bands from the SCF k grid


def scf_band_extrema(out_path):
    """EF and per-band (min, max) eigenvalues over SCF k points and spins (Hartree)."""
    ef = None; extrema = {}; inside = seen = False
    row = re.compile(r'\s*(\d+)((?:\s+-?\d+\.\d+(?:[eE][-+]?\d+)?)+)\s*$')
    with Path(out_path).open(errors='replace') as stream:
        for line in stream:
            if not inside:
                inside = 'Eigenvalues (Hartree) of SCF KS-eq.' in line
                continue
            if ef is None:
                match = re.search(r'Chemical Potential \(Hartree\)\s*=\s*(\S+)', line)
                if match:
                    ef = float(match.group(1))
                continue
            match = row.match(line)
            if match:
                seen = True
                band = int(match.group(1)) - 1
                values = [float(v) for v in match.group(2).split()]
                low, high = extrema.get(band, (np.inf, -np.inf))
                extrema[band] = (min(low, *values), max(high, *values))
            elif seen and '****' in line:
                break
    if ef is None or not extrema:
        raise RuntimeError(f'SCF eigenvalues not found in {out_path}')
    return ef, extrema


def window_band_indices(out_path, window, full):
    """Zero-based band range whose SCF energies reach the window (plus margin)."""
    ef, extrema = scf_band_extrema(out_path)
    low = ef + (window[0] - BAND_WINDOW_MARGIN_EV) / EH_EV
    high = ef + (window[1] + BAND_WINDOW_MARGIN_EV) / EH_EV
    bands = [band for band, (lo, hi) in extrema.items()
             if hi >= low and lo <= high and full[0] <= band <= full[1]]
    return (min(bands), max(bands)) if bands else tuple(full)


def window_coverage(raw_path, window):
    """(low_ok, high_ok): no band outside the stored range enters the window."""
    with h5py.File(raw_path, 'r') as f:
        attrs = f['Output'].attrs; ef = float(attrs['EF_Eh'])
        below = attrs.get('MaxEnergyBelowWindow_Eh'); above = attrs.get('MinEnergyAboveWindow_Eh')
    return (below is None or (float(below) - ef) * EH_EV < window[0],
            above is None or (float(above) - ef) * EH_EV > window[1])


def ensure_postprocessor(log=None, stop=None):
    """Build only the private heap converter; never replace the shared tools."""
    BUILD.mkdir(parents=True, exist_ok=True)
    with file_lock(BUILD / 'cpp_build.lock'):
        run_command(['make', '-C', PROJECT / 'cpp', '-j2', POSTPROC],
                    BUILD / 'postprocessor_build.log', PROJECT, expected=POSTPROC,
                    on_output=log, stop=stop)


def calculate_chunk(folder, points, log=None, stop=None, conversion=None, openmx=None, window=None):
    ensure_postprocessor(log, stop)
    folder.mkdir(parents=True, exist_ok=True)
    stem = folder.name; dat = folder / f'{stem}.dat'; raw = folder / f'{stem}.hdf5'
    overrides = validate_openmx_overrides(openmx)
    window = validate_band_window(window)
    full = request_band_range(overrides)
    spin = effective_openmx(overrides)['scf.SpinPolarization'].lower()
    restart = scf_restart_folder(overrides)
    text = input_text(stem, points, overrides)

    def converted_ok(f):
        stored = json_read(folder / 'conversion.json') or {}
        return (states_ok(f, len(points), False, False if window else full, spin)
                and (window is None or stored.get('band_window') == list(window))
                and np.allclose(f['Input/Kpath/Coordinates'][:], points, rtol=0, atol=2e-11))
    check_stop(stop)
    elapsed = 0.0
    conversion_elapsed = 0.0
    # A fully converted raw HDF5 is sufficient to resume the compact import.
    # Validate it before checking the mutable SCF restart lines or deleting any
    # OpenMX output, so a postprocessing retry can never trigger recomputation.
    try:
        with h5py.File(raw, 'r') as existing:
            raw_valid = converted_ok(existing)
    except (OSError, KeyError, ValueError):
        raw_valid = False
    if raw_valid:
        return raw, elapsed, conversion_elapsed
    if not completed_openmx(folder, stem, text):
        restart_file = restart / f'{RESTART_NAME}.crst_check'
        if restart_file.is_file():
            link = folder / f'{RESTART_NAME}_rst'
            if not link.exists():
                link.symlink_to(restart, target_is_directory=True)
            if not link.is_symlink() or link.resolve() != restart.resolve():
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
    preserve_structure_restart(folder, stem, restart)
    conversion_record = {'output': stamp(folder / f'{stem}.out'), 'postprocessor': stamp(POSTPROC)}
    if window is not None:
        conversion_record['band_window'] = list(window)
    valid = False
    try:
        with h5py.File(raw, 'r') as f:
            stored = {k: v for k, v in (json_read(folder / 'conversion.json') or {}).items() if k != 'bands'}
            valid = converted_ok(f) and stored == conversion_record
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
        if window is None:
            run_command([POSTPROC, dat], folder / 'postproc.log', folder, expected=raw,
                        on_output=postproc_output, stop=stop)
        else:
            # Convert only the bands that can reach the window.  The converter
            # reports the neighbouring band edges over every k point; widen and
            # convert again (OpenMX is not rerun) until none enters the window.
            low, high = window_band_indices(folder / f'{stem}.out', window, full)
            step = 4
            while True:
                if log:
                    log(f'[Band window] {window[0]:g} to {window[1]:g} eV: storing bands '
                        f'{low}–{high} of {full[0]}–{full[1]}')
                window_dat = folder / f'{stem}.window.dat'
                window_dat.write_text(re.sub(r'(?m)^SPADExp\.maxN\s+\d+', f'SPADExp.maxN {high}',
                                             re.sub(r'(?m)^SPADExp\.minN\s+\d+', f'SPADExp.minN {low}', text)))
                raw.unlink(missing_ok=True)
                run_command([POSTPROC, window_dat], folder / 'postproc.log', folder, expected=raw,
                            on_output=postproc_output, stop=stop)
                low_ok, high_ok = window_coverage(raw, window)
                wider = (low if low_ok else max(full[0], low - step),
                         high if high_ok else min(full[1], high + step))
                if (low_ok and high_ok) or wider == (low, high):
                    break
                low, high = wider; step *= 2
            conversion_record['bands'] = [low, high]
        conversion_elapsed = time.monotonic() - conversion_started
        with h5py.File(raw, 'r') as f:
            if not (states_ok(f, len(points), False, False if window else full, spin)
                    and np.allclose(f['Input/Kpath/Coordinates'][:], points, rtol=0, atol=2e-11)):
                raise RuntimeError(f'Incomplete eigenstates: {raw}')
        json_write(folder / 'conversion.json', conversion_record)
    return raw, elapsed, conversion_elapsed


def initialize_database(file, raw, count):
    start, end = stored_band_range(raw)
    band_count = end - start + 1
    for key, value in raw.attrs.items():
        file.attrs[key] = value
    inp = file.create_group('Input')
    for name in ('Atomic.Species', 'Atoms.SpeciesAndCoordinates', 'UnitCells'):
        raw.copy('Input/' + name, inp, name=name)
    out = file.create_group('Output'); out.attrs.update(raw['Output'].attrs)
    raw.copy('Output/Potential', out, name='Potential')
    for path in band_paths(spin_of(raw)):
        out.create_dataset(path.rsplit('/', 1)[1], shape=(count, band_count), dtype='f8',
                           chunks=(min(256, count), band_count), compression='gzip', shuffle=True)
    for path in lcao_paths(raw):
        dataset = raw[path]
        shape = (count,) + dataset.shape[1:]
        parent, name = path.rsplit('/', 1)
        file.require_group(parent).create_dataset(name, shape=shape, dtype='f4',
            chunks=(min(128, count),) + shape[1:], compression='gzip', compression_opts=4, shuffle=True, fletcher32=True)


def ensure_mesh(settings, log=None, progress=None, stop=None, conversion=None, timing=None):
    plan = mesh_plan(settings); request = plan['request']; folder = plan['folder']; final = plan['path']
    overrides = openmx_overrides(settings)
    window = band_window(settings)
    folder.mkdir(parents=True, exist_ok=True)
    with file_lock(folder / 'mesh.lock'):
        run_started = time.monotonic()
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
                                 and ('Output' not in f or request_states_ok(f, request, count, True)))
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
            def add_timing(openmx_seconds=0.0, conversion_seconds=0.0, wall_seconds=0.0, **extra):
                # Kept in the partial file, so a resumed calculation adds to the
                # time already spent, and renamed with it into the final cache.
                try:
                    record = json.loads(f.attrs.get(MESH_TIMING_ATTR, '{}'))
                except (TypeError, ValueError):
                    record = {}
                for key, value in (('openmx_seconds', openmx_seconds),
                                   ('conversion_seconds', conversion_seconds),
                                   ('wall_seconds', wall_seconds)):
                    record[key] = float(record.get(key, 0.0)) + float(value)
                record.update(extra)
                f.attrs[MESH_TIMING_ATTR] = json.dumps(record, sort_keys=True); f.flush()
                if timing and (openmx_seconds or conversion_seconds):
                    timing(openmx_seconds, conversion_seconds)

            def store_raw(raw_path, selected, points):
                # Write bounded slabs, but commit completion only after ALL slabs
                # are durable. A stopped import retries the same OpenMX job/output.
                with h5py.File(raw_path, 'r') as raw:
                    if 'Output' not in f:
                        initialize_database(f, raw, count)
                        if window is not None:
                            f['Output'].attrs['BandWindow_eV'] = window
                    if (spin_of(f) != spin_of(raw)
                            or stored_band_range(f) != stored_band_range(raw)):
                        raise RuntimeError('Spin mode or stored band range differs between OpenMX batches.')
                    if not np.isclose(f['Output'].attrs['EF_Eh'], raw['Output'].attrs['EF_Eh'], rtol=0, atol=1e-8):
                        raise RuntimeError('The SCF Fermi energy does not match the saved data.')
                    row = 0
                    for batch in selected:
                        start = batch * CHUNK_POINTS; end = min(count, start + CHUNK_POINTS)
                        check_stop(stop)
                        for path in band_paths(spin_of(raw)):
                            f[path][start:end] = raw[path][row:row+end-start]
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
                if completed_openmx(legacy, legacy.name, input_text(legacy.name, points, overrides)):
                    if log: log(f'[Reusing old OpenMX output] k points {start+1:,}–{end:,}')
                    completed_before = sum(min(CHUNK_POINTS, count - i * CHUNK_POINTS)
                                           for i in np.flatnonzero(f['Mesh/CompletedChunks'][:]))
                    def conversion_progress(current, total):
                        if conversion:
                            conversion(completed_before, end - start, current, total, count)
                    raw_path, seconds, conversion_seconds = calculate_chunk(
                        legacy, points, log, stop, conversion_progress, overrides, window)
                    elapsed += seconds
                    import_started = time.monotonic()
                    store_raw(raw_path, np.array([batch]), points)
                    conversion_seconds += time.monotonic() - import_started
                    add_timing(seconds, conversion_seconds)
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
                    job, points, log, stop, conversion_progress, overrides, window)
                elapsed += seconds
                import_started = time.monotonic()
                store_raw(raw_path, pending, points)
                conversion_seconds += time.monotonic() - import_started
                add_timing(seconds, conversion_seconds)
            add_timing(wall_seconds=time.monotonic() - run_started, k_points=int(count),
                       completed_at=datetime.now().astimezone().isoformat(timespec='seconds'))
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


ENGINE_BOHR_ANGSTROM = 0.529177  # au_ang of the SPADExp engine (setup.cpp)


def openmx_reciprocal(file, fallback=None):
    """Reciprocal vectors (1/A) of the cell stored by the OpenMX postprocessor.

    The SPADExp engine rebuilds Cartesian k from this cell.  Fractional k
    coordinates handed to the engine must therefore use the same cell, not the
    CIF-derived LATTICE: cif2omx rounds the cell to 7 decimals, and the
    resulting ~1e-8 change of |k| exceeds the engine's curved-sphere tolerance,
    so those PAD points were silently dropped.
    """
    for name in ('Input/UnitCells/Bands', 'Input/UnitCells/Atoms'):
        if isinstance(file.get(name), h5py.Dataset):
            cell = np.asarray(file[name][:], dtype=float)
            unit = file[name].parent.attrs.get('Unit', 'Ang')
            unit = unit.decode() if isinstance(unit, (bytes, np.bytes_)) else str(unit)
            if unit.lower() != 'ang':
                cell = cell * ENGINE_BOHR_ANGSTROM
            return 2 * np.pi * np.linalg.inv(cell).T
    return G if fallback is None else fallback


class EigenstateGrid:
    """Periodic nearest-neighbour reference without mixing eigenvector gauges."""
    def __init__(self, path):
        self.path = Path(path)
        self.file = h5py.File(path, 'r')
        self.points = self.file['Mesh/FirstBZCoordinates'][:]
        self.ns = np.asarray(self.file["Mesh"].attrs["Divisions"], dtype=int)
        self.paths = lcao_paths(self.file)
        self.ef = float(self.file['Output'].attrs['EF_Eh'])
        # Cell used by OpenMX and by the PAD engine for fractional -> Cartesian.
        self.reciprocal = openmx_reciprocal(self.file)
        self.spin = spin_of(self.file)
        self.band_paths = band_paths(self.spin)
        self.band_range = stored_band_range(self.file)

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

    def band_values(self, indices):
        """Stored eigenvalues (Hartree); up and down sets side by side for spin on."""
        return np.concatenate([self.take(path, indices) for path in self.band_paths], axis=1)

    def bands(self, cartesian):
        indices, error = self.lookup(cartesian)
        return (self.band_values(indices) - self.ef) * 27.211386245988, error

    def band_table(self):
        """All stored eigenvalues in eV relative to EF (k points × bands), read once."""
        if getattr(self, '_band_table', None) is None:
            values = np.concatenate([self.file[path][()] for path in self.band_paths], axis=1)
            self._band_table = (values - self.ef) * EH_EV
        return self._band_table

    def interpolate_band_table(self, cartesian, columns):
        """interpolated_bands() for the selected band columns, from band_table().

        Much faster for many points and few bands; no nearest-point distance.
        """
        table = self.band_table()[:, columns]
        target = np.asarray(cartesian, dtype=float).reshape(-1, 3)
        scaled = (target @ np.linalg.inv(G)) * self.ns
        lower = np.floor(scaled).astype(np.int64)
        fraction = scaled - lower
        values = np.zeros((len(target), table.shape[1]))
        for bits in itertools.product((0, 1), repeat=3):
            offset = np.asarray(bits, dtype=np.int64)
            indices = np.ravel_multi_index(((lower + offset) % self.ns).T, tuple(self.ns))
            weight = np.prod(np.where(offset, fraction, 1.0-fraction), axis=1)
            values += weight[:, None] * table[indices]
        return values

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
            corner = self.band_values(indices)
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
                k.attrs.update(Count=len(indices), minN=self.band_range[0], maxN=self.band_range[1],
                               FirstBZFolded=np.uint8(1),
                               UniqueCount=len(np.unique(indices)), StateSampling='periodic nearest neighbour')
                k.create_dataset('Coordinates',
                                 data=np.asarray(cartesian) @ np.linalg.inv(self.reciprocal))
                k.create_dataset('SampledFirstBZCoordinates', data=self.points[indices])
                k.create_dataset('MeshIndex', data=indices)
                k.create_dataset('SamplingDistance_Ainv', data=errors)
                group = out.create_group('Output'); group.attrs.update(self.file['Output'].attrs)
                self.file.copy('Output/Potential', group, name='Potential')
                for path in self.band_paths:
                    group.create_dataset(path.rsplit('/', 1)[1], data=self.take(path, indices))
                for path in self.paths:
                    parent, name = path.rsplit('/', 1)
                    out.require_group(parent).create_dataset(name, data=self.take(path, indices),
                                                             compression='gzip', shuffle=True)
            temporary.replace(destination)
        finally:
            temporary.unlink(missing_ok=True)
        return float(errors.max(initial=0))
