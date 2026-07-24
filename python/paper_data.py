#!/usr/bin/env python3
"""paper_data.py -- generic, reusable loaders and derived quantities for the
EXHALE method-paper figures.

Design goal: point these at ANY converged EXHALE run directory (a folder holding
``input.inp`` and an ``output/`` subdirectory) and get clean, named arrays back.
Column indices are NEVER hard-coded -- the loaders parse the ``# columns ...``
header line of each output file and build a name->index map, so they tolerate
metals-off, molecular, He 2^3S-off, and future runs that add EXTRA columns
(e.g. the Phase-2 heating breakdown or an explicit HI-Lya cooling split).

Functional API
--------------
    read_input(path)                 -> dict parsed from input.inp
    load_hydro(run_dir, adv=True)    -> {'r','rho','v','p','T','heat','cool'}
    load_ions(run_dir, adv=True)     -> {'r', species_name: array, ...}
    load_cooling_breakdown(run_dir)  -> {'r','T','ne','cool_total','reco',
                                         'coio','coex_HI','coex_HeI','coex_HeII',
                                         'brem', metal_ion: array}
    load_heating_breakdown(run_dir)  -> {'r','T','ne','heat_total','heat_HI',
                                         'heat_HeI','heat_HeII','heat_He23S',
                                         'heat_H2','heat_metals','heat_Hpe',
                                         'heat_Hdx','heat_He_recomb',
                                         'heat_He23S_Penning'} or {}
    load_excited_h(run_dir)          -> {'r','T',...,'n2s','n2p',...} or None
    adiabatic_cooling(r_cm, p, v)    -> |p * div(v)| [erg cm^-3 s^-1]
    metal_ion_names(ion_dict)        -> list of metal-ion column names present

Convenience
-----------
    load_planet(run_dir, name=None)  -> Planet object bundling all of the above.

All densities are cm^-3, rates erg cm^-3 s^-1, T in K, v in cm/s, r in R_p.
By default the *advected* profiles (``*_adv.txt``) are read because those are the
ones the transit post-processor consumes; pass ``adv=False`` for the raw
snapshot.  Rows containing NaN/Inf (e.g. trailing ghost cells beyond the escape
radius) are dropped.
"""
import os
import re
import numpy as np

# Jupiter radius [cm] (matches EXHALE and exhale_transit_lib).
RJ_CM = 7.1492e9

# Species that are NOT trace metals (used to separate metals from H/He and from
# any optional molecular columns when building the metal list generically).
HHE_SPECIES = ('HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR')
MOLECULES = ('H2', 'H2p', 'H3p', 'HeHp')

# Canonical fixed channels at the front of Cooling_breakdown.txt (the trailing
# columns are one per metal ion, in the Ion_species metal order).  The
# collisional-excitation channel is split by absorber: coex_HI is the
# Lyman-alpha-dominated H-line cooling, coex_HeI (which also carries the He
# 2^3S metastable terms) and coex_HeII the helium lines.  These are only the
# fallback names; the loader reads the actual names from the file header.
COOL_FIXED = ('r', 'T', 'ne', 'cool_total', 'reco', 'coio',
              'coex_HI', 'coex_HeI', 'coex_HeII', 'brem')

# Fixed channels of Heating_breakdown.txt (all named; no trailing metal
# columns).  Fallback only -- the loader reads the header when present.
HEAT_FIXED = ('r', 'T', 'ne', 'heat_total', 'heat_HI', 'heat_HeI',
              'heat_HeII', 'heat_He23S', 'heat_H2', 'heat_metals',
              'heat_Hpe', 'heat_Hdx', 'heat_He_recomb', 'heat_He23S_Penning')

# LaTeX-safe labels for the metal cooling channels most often plotted.
METAL_LABEL = {
    'CII': r'C\,{\sc ii}', 'CIII': r'C\,{\sc iii}',
    'OI': r'O\,{\sc i}', 'OII': r'O\,{\sc ii}', 'OIII': r'O\,{\sc iii}',
    'NII': r'N\,{\sc ii}', 'NIII': r'N\,{\sc iii}',
    'MgII': r'Mg\,{\sc ii}', 'SiII': r'Si\,{\sc ii}',
    'CaII': r'Ca\,{\sc ii}', 'NaI': r'Na\,{\sc i}',
    'FeII': r'Fe\,{\sc ii}', 'FeIII': r'Fe\,{\sc iii}',
}

# Staged Phase-1 run directories (tag -> subdir).  Override with EXHALE_FIGS_RUN
# or the --base CLI flag; the loaders themselves take an explicit path.
STAGED_DEFAULT = os.environ.get(
    'EXHALE_FIGS_RUN',
    '/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/'
    '3c17a7ed-9e0d-4907-a6ee-caeb01dbc98d/scratchpad/figs_run')

# name, tag, log10 Mdot [g/s] (final Phase-1 converged values).
PLANETS = [
    ('HD 209458 b', 'hd209',   9.57),
    ('HD 189733 b', 'hd189',   8.72),
    ('WASP-52 b',   'wasp52', 11.80),
    ('WASP-121 b',  'wasp121', 13.20),
]


# --------------------------------------------------------------------------- #
# Low-level helpers
# --------------------------------------------------------------------------- #
def _clean_name(tok):
    """Strip a units bracket and surrounding brackets from a header token.

    'r[Rp]' -> 'r', 'rho[mH/cm3]' -> 'rho', 'r/Rp' -> 'r', '[H2' -> 'H2'.
    """
    tok = tok.strip().strip('[]')
    if tok in ('r', 'r/Rp', 'r[Rp]'):
        return 'r'
    tok = tok.split('[', 1)[0]
    return tok


def _read_matrix(path):
    """Read a whitespace table, skipping '#' comments and NaN/Inf rows."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            low = s.lower()
            if 'nan' in low or 'inf' in low:
                continue
            rows.append([float(x) for x in s.split()])
    return np.array(rows)


def _columns_header(path):
    """Return the token list of the '# columns ...' header line, or None."""
    with open(path) as fh:
        for ln in fh:
            if ln.lstrip().startswith('#') and 'columns' in ln.lower():
                after = ln.split('columns', 1)[1]
                return [t for t in after.split() if t not in ('=',)]
    return None


def read_input(path):
    """Parse input.inp into a dict keyed by the text before the first ':'."""
    out = {}
    with open(path) as fh:
        for ln in fh:
            if ':' in ln:
                k, v = ln.split(':', 1)
                out[k.strip()] = v.strip()
    return out


def _odir(run_dir):
    return os.path.join(run_dir, 'output')


# --------------------------------------------------------------------------- #
# Functional loaders (header-driven, tolerant of missing/extra columns)
# --------------------------------------------------------------------------- #
def load_hydro(run_dir, adv=True):
    """Return {name: array} for Hydro_ioniz[_adv].txt (r,rho,v,p,T,heat,cool)."""
    suf = '_adv' if adv else ''
    path = os.path.join(_odir(run_dir), 'Hydro_ioniz%s.txt' % suf)
    m = _read_matrix(path)
    names = _columns_header(path)
    if names and len(names) == m.shape[1]:
        keys = [_clean_name(t) for t in names]
    else:  # fixed fallback
        keys = ['r', 'rho', 'v', 'p', 'T', 'heat', 'cool'][:m.shape[1]]
    return {k: m[:, i] for i, k in enumerate(keys)}


def load_ions(run_dir, adv=True):
    """Return {species: array} for Ion_species[_adv].txt, names from the header.

    The radius column is keyed 'r'.  Works with or without metals/molecules/
    He 2^3S because the species names are read from the '# columns' line.
    """
    suf = '_adv' if adv else ''
    path = os.path.join(_odir(run_dir), 'Ion_species%s.txt' % suf)
    m = _read_matrix(path)
    names = _columns_header(path)
    if not names or len(names) != m.shape[1]:
        raise ValueError('cannot parse Ion_species header in %s' % path)
    keys = [_clean_name(t) for t in names]
    keys[0] = 'r'
    return {k: m[:, i] for i, k in enumerate(keys)}


def metal_ion_names(ion_dict):
    """List the metal-ion columns present (excluding r, H/He stages, molecules).

    Preserves the file's column order, which is the canonical order shared with
    Cooling_breakdown.txt.
    """
    skip = set(('r',)) | set(HHE_SPECIES) | set(MOLECULES)
    return [k for k in ion_dict.keys() if k not in skip]


def load_cooling_breakdown(run_dir):
    """Return {channel: array} for Cooling_breakdown.txt.

    The 8 fixed channels are parsed from the 'colN <name>' header; the trailing
    columns are one per metal ion, whose names are taken (in order) from the
    Ion_species header so the mapping stays correct if the metal set changes.
    Extra unnamed trailing columns are kept as 'extra_k' rather than dropped, so
    a future output that appends channels still loads.  Returns {} if the file
    is absent.
    """
    path = os.path.join(_odir(run_dir), 'Cooling_breakdown.txt')
    if not os.path.exists(path):
        return {}
    m = _read_matrix(path)
    ncol = m.shape[1]

    # Fixed channels from the 'colN name' header line.
    fixed = list(COOL_FIXED)
    with open(path) as fh:
        for ln in fh:
            if ln.lstrip().startswith('#') and re.search(r'col\d+', ln):
                found = [_clean_name(t)
                         for t in re.findall(r'col\d+\s+(\S+)', ln)]
                if found:
                    fixed = found
                break

    metals = metal_ion_names(load_ions(run_dir))
    keys = list(fixed) + list(metals)
    if len(keys) < ncol:
        keys += ['extra_%d' % k for k in range(ncol - len(keys))]
    keys = keys[:ncol]
    return {k: m[:, i] for i, k in enumerate(keys)}


def load_heating_breakdown(run_dir):
    """Return {channel: array} for Heating_breakdown.txt.

    All channels are named in the 'colN <name>' header (photoionization heating
    per absorber -- H I, He I, He II, He 2^3S, H2, metals -- plus the excited-H
    photoelectric heat_Hpe, the Ly-alpha de-excitation heat_Hdx, the
    He-recombination-driven H heating, and the He 2^3S Penning heating).  The
    channel sum reproduces the heat_total column.  Returns {} if the file is
    absent.  Header-driven, so extra trailing channels still load.
    """
    path = os.path.join(_odir(run_dir), 'Heating_breakdown.txt')
    if not os.path.exists(path):
        return {}
    m = _read_matrix(path)
    ncol = m.shape[1]

    names = list(HEAT_FIXED)
    with open(path) as fh:
        for ln in fh:
            if ln.lstrip().startswith('#') and re.search(r'col\d+', ln):
                found = [_clean_name(t)
                         for t in re.findall(r'col\d+\s+(\S+)', ln)]
                if found:
                    names = found
                break
    if len(names) < ncol:
        names += ['extra_%d' % k for k in range(ncol - len(names))]
    names = names[:ncol]
    return {k: m[:, i] for i, k in enumerate(names)}


def load_excited_h(run_dir):
    """Return {name: array} for Excited_H.txt, or None if absent.

    Column names are parsed from the 'colN <name>' header line.
    """
    path = os.path.join(_odir(run_dir), 'Excited_H.txt')
    if not os.path.exists(path):
        return None
    m = _read_matrix(path)
    names = None
    with open(path) as fh:
        for ln in fh:
            if ln.lstrip().startswith('#') and re.search(r'col\d+', ln):
                names = [_clean_name(t)
                         for t in re.findall(r'col\d+\s+(\S+)', ln)]
                break
    if not names or len(names) != m.shape[1]:
        names = ['c%d' % i for i in range(m.shape[1])]
        names[0] = 'r'
    return {k: m[:, i] for i, k in enumerate(names)}


def adiabatic_cooling(r_cm, p, v):
    """Adiabatic (PdV expansion) cooling per unit volume [erg cm^-3 s^-1].

        Lambda_adiab = p * (1/r^2) d(r^2 v)/dr,   r in cm.

    The compressional-work term of the energy equation; NOT part of the
    radiative Cooling_breakdown.txt.  Centered finite difference on r^2 v; the
    unsigned magnitude is returned for log plotting (div(v) > 0 in an expanding
    wind, so it is a genuine energy sink).
    """
    r_cm = np.asarray(r_cm, float)
    f = r_cm**2 * np.asarray(v, float)
    dfdr = np.gradient(f, r_cm)
    return np.abs(np.asarray(p, float) * dfdr / r_cm**2)


# --------------------------------------------------------------------------- #
# Convenience object
# --------------------------------------------------------------------------- #
class Planet(object):
    """Bundle one converged EXHALE run: hydro, ions, cooling, excited H."""

    def __init__(self, name, run_dir, adv=True):
        self.name = name
        self.run_dir = run_dir
        self.inp = read_input(os.path.join(run_dir, 'input.inp'))
        self.Rp_RJ = float(self.inp.get('Planet radius [R_J]', 'nan'))
        self.Rp_cm = self.Rp_RJ * RJ_CM
        self.r_esc = float(self.inp.get('Escape radius [R_p]', '1e9'))

        h = load_hydro(run_dir, adv=adv)
        for k, v in h.items():
            setattr(self, k, v)          # self.r, self.rho, self.v, ...

        self.ion = load_ions(run_dir, adv=adv)
        self.r_ion = self.ion['r']
        self.metal_ions = metal_ion_names(self.ion)

        self.cool_ch = load_cooling_breakdown(run_dir)
        self.r_cool = self.cool_ch.get('r')

        self.heat_ch = load_heating_breakdown(run_dir)
        self.r_heat = self.heat_ch.get('r')

        self.exc = load_excited_h(run_dir)

        self.adiab = adiabatic_cooling(self.r * self.Rp_cm, self.p, self.v)


def load_planet(run_dir, name=None, adv=True):
    """Load a run directory (input.inp + output/) into a Planet object."""
    if name is None:
        name = read_input(os.path.join(run_dir, 'input.inp')).get(
            'Planet name', os.path.basename(os.path.normpath(run_dir)))
    return Planet(name, run_dir, adv=adv)


def staged_rundir(tag, base=None):
    """Path to a staged Phase-1 run directory for the given tag."""
    return os.path.join(base or STAGED_DEFAULT, tag)


if __name__ == '__main__':
    for nm, tag, lm in PLANETS:
        try:
            d = load_planet(staged_rundir(tag), name=nm)
            phys = d.r <= d.r_esc
            print('%-14s Rp=%.3f RJ  Tmax=%.0f K  vmax=%.1f km/s  '
                  'metals=%d  log10Mdot=%.2f'
                  % (nm, d.Rp_RJ, np.nanmax(d.T[phys]), np.nanmax(d.v) / 1e5,
                     len(d.metal_ions), lm))
        except Exception as exc:  # noqa: BLE001
            print('%-14s FAILED: %s' % (nm, exc))
