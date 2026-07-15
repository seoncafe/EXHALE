"""exhale_io.py -- reusable readers for EXHALE output files.

Every ATES run writes its profiles to an ``output*/`` directory as plain text
(columns described in the user manual, docs/EXHALE_user_manual.tex).  This module
loads those files into named, physical-unit arrays so analysis scripts and the
example notebook do not have to remember column orders.

Typical use
-----------
>>> import exhale_io as aio
>>> run = aio.load_run('WASP-121b/output_caseA2058', 'WASP-121b/input.inp.caseA')
>>> run.r, run.T, run.v_kms            # radius [R_p], T [K], velocity [km/s]
>>> run.ion['MgII'], run.ion['HII']    # ion number densities [cm^-3]
>>> aio.mdot_log10(run)                # log10 steady-state Mdot [g/s]

All densities are cm^-3, velocity cm/s (``v_kms`` is km/s), pressure erg/cm^3,
temperature K, heating/cooling erg/cm^3/s, radius in units of the base R_p.
"""

import numpy as np

# --- physical constants (cgs) ---
RJ = 6.9911e9          # Jupiter radius [cm]
MJ = 1.898e30          # Jupiter mass  [g]
Msun = 1.989e33        # Solar mass    [g]
AU = 1.495978707e13    # Astronomical unit [cm]
mu = 1.673e-24         # hydrogen mass m_H [g] (NOT a mean molecular weight;
                       # converts the m_H-density column n to g/cm^3)
GYR = 3.15576e16       # 1 Gyr [s]

# Hydrogen + helium columns of Ion_species(_adv).txt, in file order (cols 2-7).
HE_IONS = ['HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR']

# The 27 trace-metal ion columns (cols 8-34), in species_table.f90 mion order.
METAL_IONS = [
    'CI', 'CII', 'CIII', 'OI', 'OII', 'OIII', 'NI', 'NII', 'NIII',
    'MgI', 'MgII', 'MgIII', 'SiI', 'SiII', 'SiIII', 'CaI', 'CaII', 'CaIII',
    'NaI', 'NaII', 'KI', 'KII', 'SI', 'SII', 'FeI', 'FeII', 'FeIII',
]
ION_NAMES = HE_IONS + METAL_IONS

# Cooling_breakdown.txt channels: cols 1-4 are r, T, ne, cool_total; then
# 4 H/He channels, then the 27 metal-ion channels (same order as METAL_IONS).
COOL_HHE_CHANNELS = ['rec', 'coll_ion', 'coll_exc', 'brems']


class Run:
    """Container for one ATES run (one output directory)."""

    def __init__(self):
        self.r = None          # radius [R_p]
        self.n = None          # mass density in m_H units [m_H/cm^3] = rho/m_H
                               # (the rho*n0 column; metals included). n*mu -> g/cm^3.
        self.v = None          # radial velocity [cm/s]
        self.p = None          # pressure [erg/cm^3]
        self.T = None          # temperature [K]
        self.heat = None       # heating rate [erg/cm^3/s]
        self.cool = None       # cooling rate [erg/cm^3/s]
        self.ion = {}          # ion name -> number density [cm^-3]
        self.inp = {}          # parsed input.inp (see read_input)

    @property
    def v_kms(self):
        return self.v / 1.0e5

    def x_ion(self, element_stages):
        """Ionization fraction of a given stage, e.g. x_ion(['HI','HII'])['HII']
        returns nHII/(nHI+nHII). Pass the list of stages of one element."""
        tot = sum(self.ion[s] for s in element_stages)
        tot = np.where(tot > 0, tot, np.nan)
        return {s: self.ion[s] / tot for s in element_stages}


def load_hydro(path):
    """Read Hydro_ioniz.txt or Hydro_ioniz_adv.txt -> dict of physical arrays.
    Columns: r[R_p], n[m_H/cm^3] (mass density = rho/m_H, metals included),
    v[cm/s], p[erg/cm^3], T[K], heat, cool."""
    r, n, v, p, T, heat, cool = np.loadtxt(path, unpack=True)
    return dict(r=r, n=n, v=v, p=p, T=T, heat=heat, cool=cool)


def load_ions(path):
    """Read Ion_species.txt(_adv) -> (r, {ion_name: density[cm^-3]}).
    Reads all 34 columns (r + 6 H/He + 27 metals)."""
    d = np.loadtxt(path, unpack=True)
    r = d[0]
    ion = {name: d[i + 1] for i, name in enumerate(ION_NAMES) if i + 1 < d.shape[0]}
    return r, ion


def load_cooling(path):
    """Read Cooling_breakdown.txt -> dict with r, T, ne, cool_total, and a
    'chan' dict of cooling in each channel [erg/cm^3/s] (H/He + metal lines)."""
    d = np.loadtxt(path, unpack=True)
    out = dict(r=d[0], T=d[1], ne=d[2], cool_total=d[3], chan={})
    names = COOL_HHE_CHANNELS + METAL_IONS
    for i, name in enumerate(names):
        col = 4 + i
        if col < d.shape[0]:
            out['chan'][name] = d[col]
    return out


# Excited_H.txt columns (excited_hydrogen.f90 write_excited_H), in order.
EXCITED_H_COLS = [
    'r', 'T', 'nHI', 'ne', 'Jlya', 'n2s', 'n2p', 'Sproton', 'Hpe', 'Hdx',
    'tau_lya', 'S_ground_HI', 'S_coll_ion', 'S_recomb', 'Jint', 'Jstar',
]


def load_excited_H(path):
    """Read Excited_H.txt -> dict keyed by EXCITED_H_COLS (only present cols)."""
    d = np.loadtxt(path, unpack=True)
    return {name: d[i] for i, name in enumerate(EXCITED_H_COLS) if i < d.shape[0]}


def read_input(path):
    """Parse an input.inp into a dict. Captures the numeric planet/star
    parameters and the optional keyword block. Values are floats where possible.
    Keys: Rp_RJ, Mp_MJ, T0, a_AU, r_esc, HeH, Mstar_Msun, LX, LEUV,
    n0_log10, plus any 'Domain mode', 'Outer radius', etc. by label."""
    inp = {}
    with open(path) as f:
        for line in f:
            if ':' not in line:
                continue
            key, val = line.split(':', 1)
            key = key.strip()
            val = val.strip()
            inp[key] = val
    def num(label, default=None):
        for k, v in inp.items():
            if k.startswith(label):
                try:
                    return float(v.split()[0])
                except (ValueError, IndexError):
                    return default
        return default
    return dict(
        raw=inp,
        n0_log10=num('Log10 lower boundary'),
        Rp_RJ=num('Planet radius'),
        Mp_MJ=num('Planet mass'),
        T0=num('Equilibrium temperature'),
        a_AU=num('Orbital distance'),
        r_esc=num('Escape radius'),
        HeH=num('He/H'),
        Mstar_Msun=num('Parent star mass'),
        LEUV=num('Log10 of EUV luminosity'),
        spherical=('Spherical' in inp.get('Domain mode', '')),
    )


def load_run(outdir, inputfile, adv=True):
    """Load a full run: hydro + ions (+ parsed input). Uses the _adv
    (advection-corrected) files by default; set adv=False for the eq files."""
    import os
    suf = '_adv' if adv else ''
    run = Run()
    h = load_hydro(os.path.join(outdir, 'Hydro_ioniz%s.txt' % suf))
    run.r, run.n, run.v = h['r'], h['n'], h['v']
    run.p, run.T, run.heat, run.cool = h['p'], h['T'], h['heat'], h['cool']
    _, run.ion = load_ions(os.path.join(outdir, 'Ion_species%s.txt' % suf))
    if inputfile:
        run.inp = read_input(inputfile)
    return run


# Mdot output-correction factor for every '2D approximate method' the Fortran
# accepts, keyed on the method token (word 4 of the line, i.e. the first word
# of the value, before input_read.f90 rewrites 'Rate/2'->'Rate/2 + Mdot/2' and
# 'Rate/4'->'Rate/4 + Mdot'). Matches EXHALE_main.f90, which subtracts log10(2)
# for 'Rate/2 + Mdot/2' (factor 0.5) and log10(4) for 'Mdot/4' (factor 0.25)
# and leaves the output Mdot unscaled otherwise: 'Rate/4 + Mdot' already dilutes
# the incident flux by 1/4 during the run (set_energy_vectors.f90), and 'alpha'
# scales the attenuation, so neither takes an output correction.
_MDOT_FACTOR = {
    'Mdot/4': 0.25,   # output divided by 4
    'Rate/2': 0.5,    # 'Rate/2 + Mdot/2': output halved
    'Rate/4': 1.0,    # 'Rate/4 + Mdot':  flux already /4, no output correction
    'alpha':  1.0,    # alpha attenuation: no output correction
}


def _mdot_factor(method):
    """Mdot output-correction factor for a '2D approximate method' value.
    `method` is the raw input.inp value (e.g. 'Rate/2 + Mdot/2' or 'Rate/2');
    the factor is chosen by its first token, matching input_read.f90 +
    EXHALE_main.f90. Raises ValueError on an unrecognized method."""
    tok = method.split()[0] if method and method.split() else ''
    try:
        return _MDOT_FACTOR[tok]
    except KeyError:
        raise ValueError(
            "unrecognized '2D approximate method' value %r; expected one of "
            "Mdot/4, Rate/2[ + Mdot/2], Rate/4[ + Mdot], alpha" % method)


def mdot_log10(run, j_from_top=20):
    """log10 of the steady-state mass-loss rate [g/s], 4*pi*rho*v*r^2 evaluated
    near the outer boundary, with the 2D-approximation factor from input.inp's
    '2D approximate method'. Mirrors EXHALE_main.f90."""
    Rp = run.inp.get('Rp_RJ', 1.0) * RJ
    j = len(run.r) - j_from_top
    mdot = 4.0 * np.pi * run.n[j] * mu * run.v[j] * (run.r[j] * Rp) ** 2
    method = run.inp.get('raw', {}).get('2D approximate method', '')
    mdot *= _mdot_factor(method)
    return np.log10(mdot)


def mdot_Mp_per_Gyr(run, j_from_top=20):
    """Steady-state Mdot in units of the planet mass per Gyr."""
    Mp = run.inp.get('Mp_MJ', 1.0) * MJ
    return 10.0 ** mdot_log10(run, j_from_top) * GYR / Mp
