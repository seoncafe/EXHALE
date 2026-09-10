"""exhale_io.py -- reusable readers for EXHALE output files.

Every EXHALE run writes its profiles to an ``output*/`` directory as plain text
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

The readers return the PHYSICAL cells only; the ghost rows the files carry at
both ends are dropped.  Pass ``ghost=True`` to any of them (or to
``load_run``) for the file exactly as written.  See ``physical_cell_rows``.
"""

import re
import warnings

import numpy as np

# --- physical constants (cgs) ---
# RJ_CM IS THE ONE PYTHON DEFINITION of the Jupiter radius: every tool in this
# repository that turns a `Planet radius [R_J]` into a length imports it from
# here instead of writing the number down again, so a run and the analysis of
# that run cannot be about two different planets.  Its value is the IAU 2015
# nominal EQUATORIAL radius R_J^N(eq) = 7.1492e9 cm (Prsa et al. 2016, AJ 152,
# 41, Table 1), the unit transiting-planet radii are quoted in, and it is the
# same number the Fortran defines as RJ in src/modules/init/parameters.f90.
# Keep the two in step.
RJ_CM = 7.1492e9       # Jupiter equatorial radius [cm], IAU 2015 nominal
RJ = RJ_CM             # the name the readers below and the notebooks use
# The IAU 2015 nominal masses (Prsa et al. 2016, AJ 152, 41, Table 1), the
# units `Planet mass [M_J]` and `Parent star mass [M_sun]` are quoted in, and
# the same numbers parameters.f90 defines as MJ and Msun.  Keep in step: a
# tool that scaled a planet mass by 1.898e30 g while the run used 1.8982e30 g
# reported a planet 1.1e-4 lighter than the one EXHALE relaxed.
MJ = 1.8982e30         # Jupiter mass  [g], IAU 2015 nominal
Msun = 1.98842e33      # Solar mass    [g], IAU 2015 nominal
AU = 1.495978707e13    # Astronomical unit [cm]
mu = 1.67353284e-24    # hydrogen ATOM mass m_H [g] (NOT a mean molecular
                       # weight; converts the m_H-density column n to g/cm^3).
                       # Same value as mu in parameters.f90 -- keep in step.
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
# 6 H/He channels, the H3+ infrared channel and the three molecular infrared
# bands, then the 27 metal-ion channels (same order as METAL_IONS).  The
# collisional-excitation channel is written split by absorber (H I, He I,
# He II), so there are six H/He channels, not four; see the 'col5 reco ...'
# header line the writer emits.  'H3p' is zero unless the run tracks the
# molecular network; 'H2_IR', 'H2O_IR' and 'CO_IR' are zero unless
# "Molecular IR bands" is on, and are NET rates that go negative wherever the
# band heats rather than cools.
COOL_GAS_CHANNELS = ['rec', 'coll_ion', 'coex_HI', 'coex_HeI', 'coex_HeII',
                     'brems', 'H3p', 'H2_IR', 'H2O_IR', 'CO_IR']


# --- ghost cells ---------------------------------------------------------
# Every radial profile EXHALE writes -- Hydro_ioniz, Ion_species,
# Cooling_breakdown, Heating_breakdown, Excited_H, Lyman_Werner, OI_levels,
# and their _adv twins -- is written by a `do j = 1-Ng, N+Ng` loop with
# Ng = 2 fixed in parameters.f90, so the file has N + 4 rows and the first
# two and the last two are GHOST cells: the lower pair is the fixed base
# state and the upper pair a zero-gradient / WENO3 extrapolation, both filled
# by Apply_BC.  They are boundary values, not solution cells.  Taking a
# measure over the raw file counts them: the radial spread of rho*v*r^2 over
# r >= 1.2 R_p of one converged HD 189733 b state is 1.05e-2 with the ghost
# rows and 4.65e-3 without, a factor 2.3, and the factor is state-dependent
# (1.25 on WASP-121b, 1.00 on a hot Uranus), so a consumer cannot assume it
# is small.  Section 133.6 / 137 of docs/Update_EXHALE_stage1.md.
NGHOST = 2

_ROWS_HEADER = re.compile(r'physical cells are rows\s+(\d+)\s+to\s+(\d+)')


def physical_cell_rows(path):
    """0-based half-open row range ``(i0, i1)`` of the physical cells.

    Two ways of knowing, tried in this order.

    1. THE HEADER, which is the contract.  Since section 133.6 the writer
       emits ``# rows 504: 2 ghost cells at each end; physical cells are
       rows 3 to 502`` (1-based, inclusive).  It is a comment, so it costs
       no numeric parse.  Only ``Hydro_ioniz(_adv).txt`` carries it today;
       the other profile files fall through to rule 2.

    2. THE LAYOUT, for every file written before that header existed --
       the regression goldens and every saved ``output*/`` directory.
       ``Ng = 2`` is a Fortran ``parameter`` and every writer loops
       ``1-Ng..N+Ng``, so two rows at each end are ghosts unconditionally,
       for every grid type.  What the file can still be asked is whether it
       has that layout at all, and the base radius answers: ``define_grid``
       puts a cell center exactly at ``r = 1`` and it is a GHOST center --
       ``r(1-Ng) = 1`` for the Uniform and Stretched grids, ``r(2-Ng) = 1``
       for Mixed -- so ``r = 1`` lands on row 1 or row 2 and never later.
       A file whose radius column increases monotonically and equals 1 to
       roundoff at row 1 or row 2 therefore has the ghost layout, and two
       rows come off each end.  Every one of the 9214 profile files in this
       working copy satisfies it.

    A file that satisfies neither is returned whole with a warning: that is
    the pre-2026-09 behavior, and it is what a reader handed something else
    should do rather than silently cutting four rows off it.
    """
    with open(path) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            m = _ROWS_HEADER.search(line)
            if m:
                return int(m.group(1)) - 1, int(m.group(2))
    r = np.atleast_1d(np.loadtxt(path, usecols=0))
    n = r.size
    at_base = any(abs(r[i] - 1.0) <= 1.0e-12 for i in range(min(NGHOST, n)))
    if n > 2*NGHOST + 1 and at_base and bool(np.all(np.diff(r) > 0.0)):
        return NGHOST, n - NGHOST
    warnings.warn('%s: no "physical cells are rows" header and the radius '
                  'column does not have the ghost-cell layout; returning all '
                  '%d rows' % (path, n))
    return 0, n


def loadtxt_cells(path, ghost=False, **kw):
    """``np.loadtxt`` on an EXHALE profile file, ghost rows dropped.

    Every ``np.loadtxt`` keyword is passed through.  ``unpack`` is applied
    after the row selection, so a caller still gets the column-first result
    it asked for.  ``ghost=True`` returns the file as written.
    """
    unpack = kw.pop('unpack', False)
    d = np.asarray(np.loadtxt(path, **kw))
    if not ghost:
        i0, i1 = physical_cell_rows(path)
        d = d[i0:i1]
    return d.T if unpack else d


class Run:
    """Container for one EXHALE run (one output directory)."""

    def __init__(self):
        self.r = None          # radius [R_p]
        self.n = None          # mass density in m_H units [m_H/cm^3] = rho/m_H
                               # (the rho*n0 column; metals included). n*mu -> g/cm^3.
        self.v = None          # radial velocity [cm/s]
        self.p = None          # pressure [erg/cm^3]
        self.T = None          # temperature [K]
        self.heat = None       # heating rate [erg/cm^3/s]
        self.cool = None       # cooling rate [erg/cm^3/s]
        # _adv runs only: the validity of each row's temperature and of its
        # composition, in the five values of the two-field schema (0
        # corrected, 1 retained, 2 failed, 3 unsupported, 4 not evaluated),
        # and ADV_UNKNOWN where the file predates the schema.  See
        # load_hydro, load_adv_header and the file's own header block.
        self.adv_T_status = None
        self.adv_comp_status = None
        self.adv_schema = None        # 1 for a legacy single-field file
        self.adv_input_certified = None   # True / False / None (no line)
        self.adv_certification_reason = ''
        self.adv_status = None        # the single field of a legacy file
        # The measure both fields were decided by, row by row, and the
        # fraction a corrected row is accurate to in it: a corrected row is
        # a conditional correction to that fraction of itself, so the two
        # together are what a caller weighs the row by.
        self.adv_mass_row = None
        self.adv_conditional_tol = None
        self.ion = {}          # ion name -> number density [cm^-3]
        self.inp = {}          # parsed input.inp (see read_input)

    @property
    def v_kms(self):
        return self.v / 1.0e5

    # He and H nuclei carried by one particle of each species, i.e. the
    # bsp_nH / bsp_nHe weights of src/modules/init/species_table.f90.  Only
    # the species with a nonzero count appear.  HeITR is absent on purpose:
    # He 2^3S is an excited level of He I and the HeI column is the TOTAL
    # He I density, triplet included (bsp_is_excited_level), so counting the
    # triplet column as well would count those nuclei twice.
    _NUC_H = {'HI': 1, 'HII': 1, 'H2': 2, 'H2p': 2, 'H3p': 3, 'HeHp': 1,
              'OH': 1, 'H2O': 2}
    _NUC_HE = {'HeI': 1, 'HeII': 1, 'HeIII': 1, 'HeHp': 1}
    # Metal nuclei carried by a molecule. With the oxygen chemistry on the
    # OI column means FREE ATOMIC oxygen and the CI column the carbon not
    # locked in CO, so an element total is its ion stages PLUS these.
    _NUC_METAL = {'O': {'OH': 1, 'H2O': 1, 'CO': 1},
                  'C': {'CO': 1}}

    @property
    def heh_profile(self):
        """He/H ELEMENT ratio n_He/n_H per radius, nuclei counted over every
        species that carries them (H2/H2+ two H nuclei, H3+ three, HeH+ one of
        each; the He 2^3S triplet is already inside the HeI column).  Same definition
        as composition.f90 element_ratio_HeH and as the diffusion operator, so
        a run with He_diffusion on can be read against them without a second
        convention.  Species a run does not track are simply absent."""
        nuc_h = sum(w * self.ion[s] for s, w in self._NUC_H.items()
                    if s in self.ion)
        nuc_he = sum(w * self.ion[s] for s, w in self._NUC_HE.items()
                     if s in self.ion)
        return nuc_he / np.where(nuc_h > 0, nuc_h, np.nan)

    def element_density(self, sym, stages):
        """Total nuclei density of element `sym` [cm^-3]: its ion `stages`
        plus every molecule that carries it.  Use this instead of summing the
        stages when the run has the oxygen chemistry on, where the OI and CI
        columns are the FREE ATOMIC densities and about half the oxygen sits
        in OH, H2O and CO."""
        tot = sum(self.ion[s] for s in stages if s in self.ion)
        for m, w in self._NUC_METAL.get(sym, {}).items():
            if m in self.ion:
                tot = tot + w * self.ion[m]
        return tot

    def x_ion(self, element_stages):
        """Ionization fraction of a given stage, e.g. x_ion(['HI','HII'])['HII']
        returns nHII/(nHI+nHII). Pass the list of stages of one element."""
        tot = sum(self.ion[s] for s in element_stages)
        tot = np.where(tot > 0, tot, np.nan)
        return {s: self.ion[s] / tot for s in element_stages}


# Fixed part of the Hydro_ioniz(_adv).txt layout, in file order.  The file
# carries a '# columns' header and may carry more columns than these (the
# _adv file appends the two status fields), so the header is what is parsed
# and this list is only the fallback for a file written before the header
# existed.
HYDRO_COLS = ['r', 'n', 'v', 'p', 'T', 'heat', 'cool']

# The two validity fields of the advection-corrected files, and the value
# that stands for "this file does not say".  ADV_UNKNOWN is NOT one of the
# five schema values: a file written before the schema existed carries a
# single field that named only the first refusal of a row, which cannot be
# translated into a verdict of each field, so its rows are unknown and must not
# be read as corrected.
ADV_STATUS_COLS = ['adv_T_status', 'adv_comp_status']
# The measure both fields were decided by, one value per row: the fractional
# change of the face mass flux across the cell.  A corrected row is a
# CONDITIONAL correction accurate to 'adv_conditional_tol' of itself in the
# mass flux (the file's own header states the fraction), so this column is
# what a caller weighs a corrected row by.
ADV_MASS_ROW_COL = 'adv_mass_row'
ADV_UNKNOWN = -1
ADV_CORRECTED = 0
ADV_RETAINED = 1
ADV_FAILED = 2
ADV_UNSUPPORTED = 3
ADV_NOT_EVALUATED = 4
ADV_STATUS_NAME = {ADV_UNKNOWN: 'unknown', ADV_CORRECTED: 'corrected',
                   ADV_RETAINED: 'retained', ADV_FAILED: 'failed',
                   ADV_UNSUPPORTED: 'unsupported',
                   ADV_NOT_EVALUATED: 'not_evaluated'}


def load_adv_header(path):
    """Read the validity block of an advection-corrected profile file
    (write_adv_validity_header of src/modules/files_IO/write_output.f90) ->
    dict with

      schema              2 for the two-field schema, 1 for a file that
                          carries no '# adv_schema' line (legacy);
      counts              {'T': [5 ints], 'comp': [5 ints]} or {};
      input_certified     True / False / None if the file does not say;
      certification_reason  the text after the T/F, '' when certified;
      stationarity_operator  the operator the input's stationarity was
                          judged by, '' for a legacy file;
      model_restrictions  the closure's own restrictions, one string;
      product             what the file is (a one-way correction);
      conditional_tol     the fraction a corrected row is accurate to in the
                          mass flux, as a float, or None if the file does not
                          say;
      conditional_tol_text  the whole statement of that line.

    A legacy file gives schema 1 and empty text fields: it states none of
    this, and the caller must not fill it in."""
    out = dict(schema=1, counts={}, input_certified=None,
               certification_reason='', stationarity_operator='',
               model_restrictions='', product='', conditional_tol=None,
               conditional_tol_text='')
    key = None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            w = line[1:].split()
            if not w:
                continue
            if w[0] == 'adv_schema':
                out['schema'] = int(w[1])
                key = None
            elif w[0] == 'adv_status_counts':
                # '# adv_status_counts T n n n n n comp n n n n n'
                try:
                    it = w.index('T')
                    ic = w.index('comp')
                    out['counts'] = {'T': [int(x) for x in w[it + 1:ic]],
                                     'comp': [int(x) for x in w[ic + 1:]]}
                except (ValueError, IndexError):
                    pass
                key = None
            elif w[0] == 'adv_input_certified':
                out['input_certified'] = (w[1] == 'T')
                out['certification_reason'] = ' '.join(w[2:])
                key = None
            elif w[0] == 'adv_conditional_tol':
                # '# adv_conditional_tol 1.0E-02 <what it means>'
                key = 'conditional_tol_text'
                out[key] = ' '.join(w[1:])
                try:
                    out['conditional_tol'] = float(w[1])
                except (ValueError, IndexError):
                    out['conditional_tol'] = None
            elif w[0] == 'adv_stationarity_operator':
                key = 'stationarity_operator'
                out[key] = ' '.join(w[1:])
            elif w[0] == 'adv_model_restrictions':
                key = 'model_restrictions'
                out[key] = ' '.join(w[1:])
            elif w[0] == 'adv_product':
                key = 'product'
                out[key] = ' '.join(w[1:])
            elif w[0].startswith('adv_'):
                key = None
            elif key is not None:
                # Continuation line of the multi-line text fields.
                out[key] = out[key] + ' ' + ' '.join(w)
    return out


def load_hydro(path, ghost=False):
    """Read Hydro_ioniz.txt or Hydro_ioniz_adv.txt -> dict of physical arrays.
    Columns: r[R_p], n[m_H/cm^3] (mass density = rho/m_H, metals included),
    v[cm/s], p[erg/cm^3], T[K], heat, cool.

    `Hydro_ioniz_adv.txt` carries three further columns: `adv_T_status` and
    `adv_comp_status`, the validity of the row's temperature and of its
    composition in the five values of the schema (0 corrected, 1 retained,
    2 failed, 3 unsupported, 4 not evaluated), and `adv_mass_row`, the
    measure both were decided by.  A corrected row is a conditional
    correction accurate to `adv_conditional_tol` of itself in the mass flux,
    so the measure is the CONDITION of the row and comes back as a float
    array beside the two verdicts.  Column names come from the '# columns'
    header, so a file with further columns is read as written.

    A file that carries no '# adv_schema' line is LEGACY: its rows have
    unknown validity, and both fields come back filled with ADV_UNKNOWN
    rather than translated from its single `adv_status` column, which named
    only the first refusal of the row.  That column is still returned under
    its own name.  The header's own statements come back as well:
    'adv_schema', 'adv_input_certified', 'adv_certification_reason',
    'adv_stationarity_operator', 'adv_conditional_tol',
    'adv_conditional_tol_text', 'adv_model_restrictions', 'adv_product' and
    'adv_status_counts' (see load_adv_header).

    Physical cells only unless ghost=True."""
    names = None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# columns'):
                # Header names carry their unit, 'r[Rp]', 'heat[erg/cm3/s]'.
                names = [nm.split('[')[0] for nm in line.split()[2:]]
    # The first seven columns keep the keys this loader has always used
    # ('n' for the mass-density column, whose header name is 'rho'); any
    # further column is keyed by its header name.
    keys = HYDRO_COLS + (names[len(HYDRO_COLS):] if names else [])
    d = loadtxt_cells(path, ghost, unpack=True)
    out = {nm: d[i] for i, nm in enumerate(keys) if i < d.shape[0]}
    for nm in ADV_STATUS_COLS + ['adv_status']:
        if nm in out:
            out[nm] = out[nm].astype(int)
    hdr = load_adv_header(path)
    is_adv = 'adv_status' in out or ADV_STATUS_COLS[0] in out
    if is_adv:
        out['adv_schema'] = hdr['schema']
        out['adv_input_certified'] = hdr['input_certified']
        out['adv_certification_reason'] = hdr['certification_reason']
        out['adv_stationarity_operator'] = hdr['stationarity_operator']
        out['adv_conditional_tol'] = hdr['conditional_tol']
        out['adv_conditional_tol_text'] = hdr['conditional_tol_text']
        out['adv_model_restrictions'] = hdr['model_restrictions']
        out['adv_product'] = hdr['product']
        out['adv_status_counts'] = hdr['counts']
        for nm in ADV_STATUS_COLS:
            if nm not in out:
                out[nm] = np.full(out['r'].size, ADV_UNKNOWN, dtype=int)
    return out


def load_ions(path, ghost=False):
    """Read Ion_species.txt(_adv) -> (r, {ion_name: density[cm^-3]}).
    Physical cells only unless ghost=True.
    Column names come from the '# columns' header, so the molecular columns
    (H2, H2p, H3p, HeHp) and the oxygen-chemistry columns (OH, H2O, CO) are
    picked up when the run tracks them; a file without the header falls back
    to the fixed atomic order ION_NAMES."""
    names = None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# columns'):
                names = line.split()[3:]   # drop '#', 'columns', 'r[Rp]'
    if names is None:
        names = ION_NAMES
    d = loadtxt_cells(path, ghost, unpack=True)
    r = d[0]
    ion = {name: d[i + 1] for i, name in enumerate(names) if i + 1 < d.shape[0]}
    # The advection-corrected file carries the two validity fields next to
    # the species columns.  They are not species: drop them here so that a
    # caller summing this dict cannot add a status integer to a density.
    # Read them with load_hydro or load_adv_status.
    for nm in ADV_STATUS_COLS:
        ion.pop(nm, None)
    return r, ion


def load_adv_status(path, ghost=False):
    """Read the two validity fields of an advection-corrected profile file
    (either Hydro_ioniz_adv.txt or Ion_species_adv.txt) -> (T_status,
    comp_status) as integer arrays, together with the header block.

    Returns (T_status, comp_status, header).  For a legacy file, which
    carries no '# adv_schema' line, both arrays are ADV_UNKNOWN throughout:
    unknown validity, never corrected validity.  Physical cells only unless
    ghost=True."""
    names = None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if line.startswith('# columns'):
                # Every column in file order, units stripped, so the two
                # fields are found by name whatever precedes them.
                names = [nm.split('[')[0] for nm in line.split()[2:]]
    hdr = load_adv_header(path)
    d = loadtxt_cells(path, ghost, unpack=True)
    n = d[0].size
    out = []
    for nm in ADV_STATUS_COLS:
        if names is not None and nm in names:
            out.append(d[names.index(nm)].astype(int))
        else:
            out.append(np.full(n, ADV_UNKNOWN, dtype=int))
    return out[0], out[1], hdr


def load_cooling(path, ghost=False):
    """Read Cooling_breakdown.txt -> dict with r, T, ne, cool_total, and a
    'chan' dict of cooling in each channel [erg/cm^3/s] (H/He + metal lines).
    Physical cells only unless ghost=True."""
    d = loadtxt_cells(path, ghost, unpack=True)
    out = dict(r=d[0], T=d[1], ne=d[2], cool_total=d[3], chan={})
    names = COOL_GAS_CHANNELS + METAL_IONS
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


def load_excited_H(path, ghost=False):
    """Read Excited_H.txt -> dict keyed by EXCITED_H_COLS (only present cols).
    Physical cells only unless ghost=True."""
    d = loadtxt_cells(path, ghost, unpack=True)
    return {name: d[i] for i, name in enumerate(EXCITED_H_COLS) if i < d.shape[0]}


# Lyman_Werner.txt columns (write_output.f90), in order.
LYMAN_WERNER_COLS = [
    'r', 'T', 'x_H2', 'nH2', 'NH2', 'f_shield', 'k_LW', 'heat_LW',
]


def load_lyman_werner(path, ghost=False):
    """Read Lyman_Werner.txt -> dict keyed by LYMAN_WERNER_COLS (physical
    cells only unless ghost=True). Written only
    by a molecular run that carries a "Stellar LW flux": the star-ward H2
    column, the Draine & Bertoldi (1996) self-shielding factor, the
    photodissociation rate [1/s] and its heating [erg/cm^3/s]."""
    d = loadtxt_cells(path, ghost, unpack=True)
    return {name: d[i] for i, name in enumerate(LYMAN_WERNER_COLS)
            if i < d.shape[0]}


# OI_levels.txt columns (write_output.f90), in order.
OI_LEVEL_COLS = [
    'r', 'T', 'ne', 'nHI', 'nOI',
    'f_3P2', 'f_3P1', 'f_3P0', 'n_3P2', 'n_3P1', 'n_3P0',
]


def load_OI_levels(path, ghost=False):
    """Read OI_levels.txt / OI_levels_adv.txt -> dict keyed by OI_LEVEL_COLS
    (physical cells only unless ghost=True).

    The fractional populations of the three O I 2p4 3P ground-term levels,
    solved in the same three-level statistical equilibrium as the [O I]
    63.2/145.5/44.1 um cooling (Cool_coeff.f90), plus the level densities
    n_3P2 + n_3P1 + n_3P0 = nOI. Written by metal-bearing runs only; the
    3P2 / 3P1 / 3P0 levels are the lower levels of the O I 1302.168 /
    1304.858 / 1306.029 A resonance triplet."""
    d = loadtxt_cells(path, ghost, unpack=True)
    return {name: d[i] for i, name in enumerate(OI_LEVEL_COLS)
            if i < d.shape[0]}


def load_lower_atmosphere_profile(path):
    """Read a lower-atmosphere profile file (docs/input_schema.md section 2d).

    Returns (header, columns): `header` maps every `# key value` line to its
    value string, `columns` maps every name of the `# columns:` line to its
    column. Columns this file has no consumer for are kept, not dropped, and
    everything is indexed by name -- the producers order their element list
    differently from run to run.
    """
    header, names = {}, None
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                continue
            body = line[1:].strip()
            if not body:
                continue
            parts = body.split(None, 1)
            key = parts[0].rstrip(':')
            val = parts[1].strip() if len(parts) > 1 else ''
            if key == 'columns':
                names = val.split()
            else:
                header[key] = val
    if names is None:
        raise SystemExit('%s: no "# columns:" header' % path)
    # Plain np.loadtxt, not loadtxt_cells: this is an INPUT file handed to a
    # run (written by run_lower.py / a photochemistry code), not one of the
    # profiles write_output.f90 emits, so it has no ghost rows to drop and
    # its first column is a pressure, not a radius.
    data = np.loadtxt(path, comments='#', ndmin=2)
    if data.shape[1] != len(names):
        raise SystemExit('%s: %d columns, %d names'
                         % (path, data.shape[1], len(names)))
    return header, {nm: data[:, i] for i, nm in enumerate(names)}


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


def load_run(outdir, inputfile, adv=True, ghost=False):
    """Load a full run: hydro + ions (+ parsed input). Uses the _adv
    (advection-corrected) files by default; set adv=False for the eq files.
    Physical cells only unless ghost=True."""
    import os
    suf = '_adv' if adv else ''
    run = Run()
    h = load_hydro(os.path.join(outdir, 'Hydro_ioniz%s.txt' % suf), ghost)
    run.r, run.n, run.v = h['r'], h['n'], h['v']
    run.p, run.T, run.heat, run.cool = h['p'], h['T'], h['heat'], h['cool']
    # Present in the _adv file only (see load_hydro); None for the eq file,
    # whose rows are the run's own solution by construction.  A legacy _adv
    # file gives schema 1 and ADV_UNKNOWN in both fields.
    run.adv_T_status = h.get('adv_T_status')
    run.adv_comp_status = h.get('adv_comp_status')
    run.adv_schema = h.get('adv_schema')
    run.adv_input_certified = h.get('adv_input_certified')
    run.adv_certification_reason = h.get('adv_certification_reason', '')
    run.adv_status = h.get('adv_status')
    run.adv_mass_row = h.get('adv_mass_row')
    run.adv_conditional_tol = h.get('adv_conditional_tol')
    _, run.ion = load_ions(os.path.join(outdir, 'Ion_species%s.txt' % suf),
                           ghost)
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
    'Mdot':   1.0,    # full sphere, no 2D approximation: no output correction
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
            "Mdot, Mdot/4, Rate/2[ + Mdot/2], Rate/4[ + Mdot], alpha" % method)


def mdot_log10(run, j_from_top=20):
    """log10 of the steady-state mass-loss rate [g/s], 4*pi*rho*v*r^2 evaluated
    near the outer boundary, with the 2D-approximation factor from input.inp's
    '2D approximate method'. Mirrors EXHALE_main.f90, which evaluates it at the
    PHYSICAL cell j = N - j_from_top (EXHALE_main.f90:1514) -- 0-based index
    len(r) - 1 - j_from_top of a ghost-free array. It used to be indexed
    len(r) - j_from_top on an array that still carried the ghost rows, i.e.
    three cells further out than the run's own Mdot."""
    Rp = run.inp.get('Rp_RJ', 1.0) * RJ
    j = len(run.r) - 1 - j_from_top
    mdot = 4.0 * np.pi * run.n[j] * mu * run.v[j] * (run.r[j] * Rp) ** 2
    method = run.inp.get('raw', {}).get('2D approximate method', '')
    mdot *= _mdot_factor(method)
    return np.log10(mdot)


def mdot_Mp_per_Gyr(run, j_from_top=20):
    """Steady-state Mdot in units of the planet mass per Gyr."""
    Mp = run.inp.get('Mp_MJ', 1.0) * MJ
    return 10.0 ** mdot_log10(run, j_from_top) * GYR / Mp
