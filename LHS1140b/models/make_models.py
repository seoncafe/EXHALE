#!/usr/bin/env python3
"""Write the LHS 1140 b model tree of MODELS.md section 3.

One directory per case, `models/<group>/HeH<value>/`, holding the `input.inp`
exactly as it will be solved and whatever the lower boundary needs beside it
(`base.inp`, or a copy of the photochemical `lower_atmosphere_profile.dat`).
The group name states the physics and is parsed, not repeated: it is

    <chemistry>_<lower boundary>_<spectrum>_<mixing>

and the four fields decide every key that differs between groups
(MODELS.md section 2).  Nothing else in the file changes from one case to
the next, so a key that is wrong is wrong in one place.

A case that already carries `output/Hydro_ioniz.txt` is a solved case and is
left untouched; `--force` overwrites it.  Run from anywhere:

    python3 models/make_models.py [--force] [--only <substring>] [--list]

The runners are `run_case.sh` (a prescribed-composition case) and
`run_closure.sh` (a flux-closure rung); `status.py` reads the tree back.
"""

import argparse
import os
import re
import shutil
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))          # .../LHS1140b/models
LHS = os.path.dirname(HERE)                                # .../LHS1140b
EXHALE = os.path.dirname(LHS)                              # .../EXHALE_v1.00

# --------------------------------------------------------------------------
# The planet, and the parts of input.inp that no group changes
# --------------------------------------------------------------------------

# LHS 1140 b as `LHS1140b/system_parameters.md` fixes it; the scalar lower
# boundary is T = 226 K at R_0 = 0.157692 R_J and p = 1 microbar
# (MODELS.md section 2).
PLANET = dict(
    radius_RJ='0.157692',
    mass_MJ='0.0176220',
    T_eq_K='226.0',
    a_AU='0.0946',
    r_esc_Rp='2.00',
    Mstar_Msun='0.1844',
    Rstar_Rsun='0.2159',
    log_LX='26.372',
    log_LEUV='26.404',
    # log10 n_tot at 1 microbar and 226 K.  Inert wherever "Base BC:
    # pressure" is in force, which is every case here, and carried only
    # because the solutions of record under examples/ carry it.
    log10_n0='13.506',
)

T_BASE_K = '226.0'          # the scalar base temperature
R_BASE_RJ = '0.157692'      # the scalar base radius
P_BASE_BAR = '1.0e-6'       # the scalar base level, 1 microbar
KZZ_BASE = '1.0e9'          # the eddy coefficient the molecular base carries

# The residual tolerance each lower boundary was solved at
# (`examples/scalar_base_kzz1e9`, `examples/scalar_base_cno`,
# `examples/photochem_profile_base`).
RESID_TOL = {'scalar': '5.0e-5', 'scalarCNO': '4.0e-6', 'photochem': '4.0e-6'}

# --------------------------------------------------------------------------
# The molecular base
# --------------------------------------------------------------------------

# Fraction of the hydrogen nuclei bound into H2 in the photochemical column at
# the matching level: f = 2 q_H2 / (2 q_H2 + q_H) from
# `lower_profile/lower_atmosphere_profile.dat` at p = 1e-6 bar
# (q_H2 = 0.19265464926127412, q_H = 6.488704180144029e-06).  The scalar
# molecular cases hold that fraction at every He/H (MODELS.md section 4).
#
# It puts q_H2_base 2.0e-5 below the ceiling 0.5/(0.5 + He/H) the reader
# enforces, and that proximity is NOT what stops the cold march of a
# molecular case: at f = 0.9998, which leaves 2.1e-4 to 3.7e-4 of room, the
# same run stops the same way 197 steps later
# (docs/lhs1140b_stationary_L4c_20260913.md section 3).  The column's own
# value is therefore kept.
H2_NUCLEUS_FRACTION = 0.9999831600354949


def h2_mixing_ratio_base(heh):
    """q_H2 = n_H2/(n_H2 + n_H + n_He) of the inflowing gas at fixed f.

    Per hydrogen nucleus the mixture carries f/2 molecules, (1 - f) atoms and
    He/H helium atoms, so

        q_H2 = (f/2) / ((1 - f) + f/2 + He/H),

    which approaches the ceiling 0.5/(0.5 + He/H) that `input_read.f90`
    enforces (every H nucleus bound) as f -> 1.
    """
    f = H2_NUCLEUS_FRACTION
    return (f/2.0)/((1.0 - f) + f/2.0 + heh)


def h2_mixing_ratio_ceiling(heh):
    return 0.5/(0.5 + heh)


# --------------------------------------------------------------------------
# The group table of MODELS.md section 3
# --------------------------------------------------------------------------

# Group name -> the He/H ladder, written exactly as MODELS.md writes it: the
# string is both the case name suffix and the value of "He/H number ratio",
# so the table and the tree cannot drift apart.
XUV_SCALES = ('0.01', '0.10', '0.15', '0.20', '0.25', '0.30', '0.33')

GROUPS = [
    ('atomic_scalar_gj1132_wellmixed',
     ['0.083', '0.40', '0.42', '0.44', '0.55', '1', '10', '100', '1000']),
    ('atomic_scalar_gj699_wellmixed',
     ['0.042', '0.046', '0.050', '0.083', '1', '1000']),
    ('atomic_scalar_gj1132_kzz0', ['0.55', '2.6', '3.0', '3.5', '3.7', '3.9']),
    ('atomic_scalar_gj1132_kzz1e5', ['2.6', '3.0', '3.35', '3.64', '3.93']),
    ('atomic_scalar_gj1132_kzz1e6', ['0.55', '2.4', '2.8', '3.19', '3.46', '3.74']),
    ('atomic_scalar_gj1132_kzz1e7', ['0.55', '2.70', '2.94', '3.18']),
    ('atomic_scalar_gj1132_kzz1e8', ['0.55', '2.05', '2.23', '2.41']),
    ('atomic_scalar_gj1132_kzz1e9',
     ['0.55', '1.50', '1.60', '1.70', '2.13', '4.0', '9.7']),
    ('atomic_scalar_gj1132_kzz1e10', ['0.55', '1.06', '1.15', '1.29']),
    ('atomic_scalar_gj1132_kzz1e11', ['0.55', '0.795', '0.865', '0.93']),
    ('atomic_scalarCNO_gj1132_kzz1e9', ['2.13']),
    ('atomic_photochem_gj1132_kzzprofile',
     ['2.09', '3', '5', '7', '8', '9', '9.7', '10', '12']),
] + [
    ('atomic_scalar_gj1132x%s_kzz1e9' % x, ['2.13', '9.7'])
    for x in XUV_SCALES
] + [
    ('atomic_photochem_gj1132x%s_kzzprofile' % x, ['9.7'])
    for x in XUV_SCALES
] + [
    ('molecular_scalar_gj1132_wellmixed', ['0.083', '0.55', '2.13']),
    ('molecular_scalar_gj1132_kzz1e9', ['0.083', '0.55', '2.13', '9.7']),
    ('molecular_photochem_gj1132_kzzprofile', ['2.09', '9']),
]

# EVERY CASE THAT STANDS ON A PHOTOCHEMICAL COLUMN STANDS ON A CLOSURE
# RUNG'S OWN CONVERGED COLUMN, and on no other.  The case name is the rung's
# name; the He/H written into input.inp is the value that column carries at
# the matching level (X_He at p = 1e-6 bar, `profile_match_HeH` below), so
# the case and the column cannot disagree on paper or at startup.
def profile_source(heh):
    """The converged iterate's column of the closure rung at this He/H
    (`atomic_photochem_gj1132_kzzprofile/HeH<heh>/k<k_conv>/`), or None when
    that rung has not converged yet.

    CONVERGED, and not merely certified.  Every iterate of a closure is a
    certified wind -- that is what the closure measures the elemental fluxes
    of -- so "the newest certified k" names an iterate of a closure still in
    progress, whose composition is still moving.  The rung's own log states
    which k closed it, and that is the only k whose column a case may stand
    on; while a rung is being re-run, its dependent cases are reported as
    not yet available rather than written from a passing iterate.

    Two kinds of case read it, for the same reason.  An atomic group at a
    scaled XUV normalization holds the fiducial column fixed and re-solves
    the wind on the scaled spectrum: a Photochem climate solve on a spectrum
    scaled as a whole freezes the deep atmosphere (136 K at 17 bar on the
    0.01 grid) and its chemistry fails elemental closure.  A molecular group
    holds the column of its OWN atomic rung: the molecular case is seeded
    from that rung's certified wind, and a seed is admissible only when the
    two states agree on the column and on the reservoir to a part in 1e6
    (load_IC.f90, metadata field "reservoir")."""
    rung = os.path.join(HERE, 'atomic_photochem_gj1132_kzzprofile', 'HeH%s' % heh)
    log = os.path.join(rung, 'closure.log')
    if not os.path.isfile(log):
        return None
    m = None
    for line in open(log, errors='replace'):
        hit = re.search(r'CONVERGED at k=(\d+)', line)
        if hit:
            m = int(hit.group(1))
    if m is None:
        return None
    k = 'k%02d' % m
    prof = os.path.join(rung, k, 'lower_atmosphere_profile.dat')
    cert = os.path.join(rung, k, 'run.log')
    if os.path.isfile(prof) and os.path.isfile(cert) and \
       'CERTIFIED: every active equation' in open(cert, errors='replace').read():
        return prof
    return None


# --------------------------------------------------------------------------
# The numerical flux of a case
# --------------------------------------------------------------------------

# HLLC everywhere, except the seven lowest-XUV atomic cases below.  There the
# HLLC flux does not reach the stationary root on the 500-cell catalog grid:
# the same case solved with HLLC on a grid of twice the cells reaches the
# state the Roe flux reaches on the catalog grid, so what stands between the
# two fluxes at 500 cells is resolution and not a different wind
# (`docs/lhs1140b_stationary_L25_20260916.md`, step 3, sections 3.3 and 3.5;
# user decision of 2026-09-17).  Where both fluxes solve, Roe against HLLC
# moves the mass-loss rate and the He 10830 equivalent width down by 0.3 to
# 0.6 per cent and the first cell by 6 to 8 per cent, colder and denser, with
# the mass-flux spread unchanged (the same memo, section 3.4).
FLUX_DEFAULT = 'HLLC'

ROE_CASES = frozenset([
    'atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7',
    'atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7',
    'atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7',
    'atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7',
    'atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13',
    'atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7',
    'atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7',
])

# The sentence the input carries beside the key, so that a state's own file
# says why it is not on the flux the rest of the catalog is on.
ROE_NOTE = (
    'Numerical flux: ROE, the HLLC flux does not reach the stationary root'
    ' on this grid at this XUV level (L25 step 3); the Roe flux systematic'
    ' against HLLC where both solve is Mdot and the He 10830 equivalent'
    ' width lower by 0.3 to 0.6 per cent, the first cell colder by 6 to 7'
    ' per cent and denser by 6 to 8 per cent, and the mass-flux spread'
    ' unchanged (L25 memo).')


def case_flux(group_name, heh):
    """`Numerical flux:` of one case: ROE for the seven, HLLC for the rest."""
    return 'ROE' if ('%s/HeH%s' % (group_name, heh)) in ROE_CASES \
        else FLUX_DEFAULT


# What each group is, in one sentence, for the two comment lines at the head
# of every input.inp.  Keyed on the chemistry and the lower boundary, which
# is what the physics of a group actually turns on.
PHYSICS_NOTE = {
    ('atomic', 'scalar'):
        'H, He, He(2^3S) and electrons above a prescribed scalar base'
        ' (T = %s K, R_0 = %s R_J, p = %s bar), metal-free.'
        % (T_BASE_K, R_BASE_RJ, P_BASE_BAR),
    ('atomic', 'scalarCNO'):
        'The same atomic wind over the same scalar base, plus the C, N and O'
        ' reservoirs the photochemical column carries at the matching level.',
    ('atomic', 'photochem'):
        'The atomic wind over the Photochem column handed across as a'
        ' profile; the reservoir He/H is solved by the elemental-flux'
        ' closure, so the composition is an output.',
    ('molecular', 'scalar'):
        'H2, H2+, H3+ and HeH+ with the H2 carrier transported, over the'
        ' scalar molecular base (q_H2_base at p = %s bar).' % P_BASE_BAR,
    ('molecular', 'photochem'):
        'The molecular layer solved above the converged column of its own'
        ' atomic closure rung, held fixed; the profile owns q_H2 and K_zz'
        ' at the matching level.',
}


# --------------------------------------------------------------------------
# Parsing a group name into the four fields
# --------------------------------------------------------------------------

class Group(object):
    """One group name resolved into the keys it implies."""

    def __init__(self, name):
        parts = name.split('_')
        if len(parts) != 4:
            raise SystemExit('group %r is not <chemistry>_<lower boundary>'
                             '_<spectrum>_<mixing>' % name)
        self.name = name
        self.chemistry, self.boundary, self.spectrum, self.mixing = parts
        if self.chemistry not in ('atomic', 'molecular'):
            raise SystemExit('%s: unknown chemistry %r'
                             % (name, self.chemistry))
        if self.boundary not in ('scalar', 'scalarCNO', 'photochem'):
            raise SystemExit('%s: unknown lower boundary %r'
                             % (name, self.boundary))

    # ---- the spectrum ---------------------------------------------------
    @property
    def spectrum_file(self):
        """The SED of this group, as `LHS1140b/sed/` names it."""
        if self.spectrum == 'gj1132':
            return 'lhs1140_sed_gj1132_at_b.txt'
        if self.spectrum == 'gj699':
            return 'lhs1140_sed_gj699_at_b.txt'
        if self.spectrum.startswith('gj1132x'):
            scale = self.spectrum[len('gj1132x'):]
            # 0.30 -> xuv0p30; the file names carry two decimals.
            tag = ('%.2f' % float(scale)).replace('.', 'p')
            return 'lhs1140_sed_gj1132_at_b_xuv%s.txt' % tag
        raise SystemExit('%s: unknown spectrum %r' % (self.name,
                                                      self.spectrum))

    # ---- the mixing -----------------------------------------------------
    @property
    def he_diffusion(self):
        return self.mixing != 'wellmixed'

    @property
    def he_kzz(self):
        """The scalar eddy coefficient, or None when the profile supplies
        K_zz(p) (`kzzprofile`) or there is no diffusion at all."""
        if self.mixing in ('wellmixed', 'kzzprofile'):
            return None
        if not self.mixing.startswith('kzz'):
            raise SystemExit('%s: unknown mixing %r' % (self.name,
                                                        self.mixing))
        value = self.mixing[len('kzz'):]
        if value == '0':
            return '0.0'
        if not value.startswith('1e'):
            raise SystemExit('%s: unknown mixing %r' % (self.name,
                                                        self.mixing))
        return '1.0e%d' % int(value[2:])

    # ---- what the case needs beside input.inp ---------------------------
    @property
    def is_closure(self):
        """A rung whose composition the elemental-flux closure solves, rather
        than a case whose composition is prescribed.  Only the ATOMIC profile
        groups are closed; the molecular ones hold the stored column fixed
        (MODELS.md section 3)."""
        return (self.boundary == 'photochem' and self.chemistry == 'atomic'
                and 'x' not in self.spectrum)

    @property
    def uses_profile(self):
        return self.boundary == 'photochem'

    @property
    def resid_tol(self):
        return RESID_TOL[self.boundary]

    @property
    def physics_note(self):
        return PHYSICS_NOTE[(self.chemistry, self.boundary)]


# --------------------------------------------------------------------------
# input.inp
# --------------------------------------------------------------------------

def input_inp(group, heh, for_closure_template=False, heh_value=None):
    """The whole file, in the order the solutions of record carry it.

    Two routes, and the group decides which (MODELS.md section 6, "Recipe
    adopted 2026-09-13").  An ATOMIC case is solved by the partitioned
    stationary route from an archived state of the same physics mapped onto
    the current grid: `Load IC? True`, `Restart intent: stationary`,
    `Solver: Newton`, no marching thresholds.  A MOLECULAR case has no
    archived state to start from -- no run before 2026-09-13 turned the
    network on -- so it marches from the code's own initial condition,
    two-stage PLM then WENO3, with the JFNK hand-off.

    `Well balanced: True` is in both.  At the base Mach number of this
    planet, 3e-7, the plain HLLC contact speed is set by the reconstruction's
    pressure truncation error rather than by the flow, and the energy row
    then measures itself against its own largest term
    (`docs/lhs1140b_stationary_L5c_20260913.md`); the well-balanced operator
    carries the departure from each cell's own hydrostatic equilibrium
    instead, and the same seed, input and binary then converge quadratically.

    `for_closure_template` writes the form `element_flux_closure.py` drives:
    the seed solution is the initial condition and the solve runs directly,
    so the marching keys are left out (the driver drops them anyway), the
    profile name is the driver's to write, and the route's own keys are in
    `closure.json` under `input_keys`.

    `Numerical flux:` is `case_flux(group.name, heh)`: HLLC for the catalog
    and ROE for the seven lowest-XUV atomic cases of `ROE_CASES`, whose head
    comment then carries `ROE_NOTE`.

    `heh_value` is the number written on the "He/H number ratio" line when it
    differs from the case name: a frozen profile owns the reservoir at the
    matching level and the case name is that value rounded.
    """
    if heh_value is None:
        heh_value = heh
    note = group.physics_note
    flux = case_flux(group.name, heh)
    lines = ['# %s' % group.name,
             '# %s' % note]
    if flux == 'ROE':
        lines += ['# %s' % line for line in textwrap.wrap(ROE_NOTE, 72)]
    lines.append('Planet name: LHS1140b_%s_HeH%s' % (group.name, heh))

    # The base level.  "Base BC: pressure" states it for every case; the
    # density key is redundant there and is refused outright beside a handoff
    # that states its own level (input_read.f90, "two different base levels"),
    # which is what base.inp's p_base and a profile's p_match both do.
    if not (group.uses_profile or group.chemistry == 'molecular'):
        lines.append('Log10 lower boundary number density [cm^-3]: %s'
                     % PLANET['log10_n0'])

    lines += [
        'Planet radius [R_J]: %s' % PLANET['radius_RJ'],
        'Planet mass [M_J]: %s' % PLANET['mass_MJ'],
        'Equilibrium temperature [K]: %s' % PLANET['T_eq_K'],
        'Orbital distance [AU]: %s' % PLANET['a_AU'],
        'Escape radius [R_p]: %s' % PLANET['r_esc_Rp'],
        'He/H number ratio: %s' % heh_value,
        '2D approximate method: Mdot',
        'Parent star mass [M_sun]: %s' % PLANET['Mstar_Msun'],
        'Spectrum type: Load from file..',
        # The closure driver runs EXHALE in the iteration directory, one
        # level below the case, so the template reaches the spectrum from
        # there and a prescribed case from the case directory itself.
        'Spectrum file: %s/sed/%s'
        % ('../../../..' if for_closure_template else '../../..',
           group.spectrum_file),
        'Use only EUV? False',
        '[E_low,E_mid,E_high] = [ 13.60 , 123.98 , 1.24e3 ]',
        'Log10 of X-ray luminosity [erg/s]: %s' % PLANET['log_LX'],
        'Log10 of EUV luminosity [erg/s]: %s' % PLANET['log_LEUV'],
        'Grid type: Mixed',
        'Numerical flux: %s' % flux,
    ]

    # The stationary route rebuilds the state under WENO3 whatever this line
    # says (the residual of a WENO3 solution may not be assembled by the
    # other operator), so it carries the single-stage name; a marching case
    # states the two stages it actually walks through.
    #
    # A MOLECULAR CASE TAKES THE SAME ROUTE.  Its state is not marched from
    # the code's own cold initial condition: the runner builds it out of the
    # certified atomic solution of the same name and He/H
    # (`EXHALE_MOLECULAR_SEED`), so the case is solved from a loaded state
    # like every atomic one, and the marching keys would throw that state
    # away -- `run_case.sh` reads `Load IC?` from this file for the wind pass
    # (docs/lhs1140b_stationary_L7c_20260914.md section 7).
    stationary = for_closure_template or group.chemistry in ('atomic',
                                                             'molecular')
    if stationary:
        lines += ['Reconstruction scheme: PLM',
                  'Include He23S? True',
                  'Load IC? True',
                  'Do only PP: False']
    else:
        lines += ['Reconstruction scheme: PLM+WENO3',
                  'Include He23S? True',
                  'Load IC? False',
                  'Do only PP: False']
    lines.append('Force start: False')

    if group.he_diffusion:
        lines.append('He_diffusion: True')

    lines += [
        'Domain mode: Spherical',
        'Outer radius [R_p]: 30.0',
        'Stellar radius [R_sun]: %s' % PLANET['Rstar_Rsun'],
        'Base BC: pressure 1.0',
        'Resid tol: %s' % group.resid_tol,
    ]

    if group.he_kzz is not None:
        lines.append('He_Kzz: %s' % group.he_kzz)

    if not for_closure_template:
        if stationary:
            # The seed is measured as it stands and the partitioned solve is
            # entered at once; "Restart intent: stationary" is refused
            # without "Solver: Newton" (input_schema.md K43).
            #
            # The secondary-ionization coupling is applied from the first
            # evaluation.  Its STAGED default arms only at a marching
            # convergence, which this route never reaches -- it takes no time
            # step -- so a state solved under the default carries
            # `sec_ion=F` while the post-processing pass, which applies the
            # coupling unconditionally, reports the wind with it: on this
            # planet the two differ by a factor 2.6 in Mdot.  The danger the
            # staging was written for is a cold start, and this route starts
            # from a solved state (input_schema.md K14c).
            #
            # A molecular case adds `equilibrate`: the seed carries the
            # composition the conversion put on it cell by cell, and the
            # measurement the stationary restart makes is of a state whose
            # composition is its own fixed point, not of one a single
            # evaluation is still moving.
            intent = ('stationary equilibrate'
                      if group.chemistry == 'molecular' else 'stationary')
            lines += ['Solver: Newton',
                      'Restart intent: %s' % intent,
                      'Secondary_ionization: Immediate']
        else:
            # Two-stage marching (PLM to du < 0.5, WENO3 to du < 1e-3) and
            # then the JFNK finish, the recipe of README_HOWTO.md.
            lines += ['du_th [PLM,WENO3]: 0.5 1.0e-3',
                      'Solver: Newton']

    if not for_closure_template:
        lines.append('Well balanced: True')

    if group.chemistry == 'molecular':
        lines += ['Molecular chemistry: True',
                  'Molecular carrier transport: True']

    if group.uses_profile and not for_closure_template:
        lines.append('Lower atmosphere profile: lower_atmosphere_profile.dat')

    return '\n'.join(lines) + '\n'


# --------------------------------------------------------------------------
# base.inp
# --------------------------------------------------------------------------

CNO_TEMPLATE = os.path.join(LHS, 'examples', 'scalar_base_cno', 'base.inp')


def base_inp_scalar_cno(heh):
    """The scalar family's base with the C, N and O reservoirs of the
    photochemical column added.

    Only the three element ratios come from
    `examples/scalar_base_cno/base.inp`; the base state is the scalar
    family's (T = 226 K, R_0 = 0.157692 R_J), so this case differs from the
    atomic scalar case at the same He/H in the reservoirs alone, which is
    what MODELS.md section 2 defines `scalarCNO` to be.

    `q_H2_base` of the template is NOT carried over.  Without a molecular
    network it changes nothing (`comp_ntot_bc` uses the H2 binding only under
    `molecular_base`), and at He/H = 2.13 the column's 0.19265 stands above
    the attainable ceiling 0.5/(0.5 + He/H) = 0.190114, which `input_read`
    refuses at startup.
    """
    keep = ('C_H_base', 'N_H_base', 'O_H_base')
    values = {}
    with open(CNO_TEMPLATE) as fh:
        for line in fh:
            if line.lstrip().startswith('#') or not line.strip():
                continue
            word = line.split()
            if len(word) >= 2 and word[0] in keep:
                values[word[0]] = word[1]
    missing = [k for k in keep if k not in values]
    if missing:
        raise SystemExit('%s lacks %s' % (CNO_TEMPLATE, ', '.join(missing)))

    out = ["# base.inp -- the scalar family's base (226 K, 0.157692 R_J)",
           '# with the C, N, O reservoirs of the photochemical column at its',
           '# 1 microbar matching level.  The three element ratios are copied',
           '# from LHS1140b/examples/scalar_base_cno/base.inp; everything else',
           "# is the scalar family's, so the only change against the atomic",
           '# scalar case at the same He/H is the reservoirs.',
           'T_base    %s' % T_BASE_K,
           'r_base    %s' % R_BASE_RJ,
           'HeH_base  %s' % heh,
           'C_H_base  %s' % values['C_H_base'],
           'N_H_base  %s' % values['N_H_base'],
           'O_H_base  %s' % values['O_H_base']]
    return '\n'.join(out) + '\n'


def base_inp_molecular(heh, kzz):
    """The molecular base of MODELS.md section 4: the scalar base state with
    the H2 mixing ratio that holds the column's nucleus fraction f."""
    q = h2_mixing_ratio_base(float(heh))
    ceiling = h2_mixing_ratio_ceiling(float(heh))
    if q > ceiling:
        raise SystemExit('q_H2_base %.6f above the ceiling %.6f at He/H = %s'
                         % (q, ceiling, heh))
    out = ['# base.inp -- the molecular lower boundary of MODELS.md'
           ' section 4.',
           '#',
           '# The inflowing gas keeps the hydrogen-nucleus fraction bound',
           '# into H2 that the photochemical column carries at 1 microbar,',
           '# f = %.16f, so' % H2_NUCLEUS_FRACTION,
           '#',
           '#     q_H2_base = (f/2) / ((1 - f) + f/2 + He/H)',
           '#             = %.9f at He/H = %s,' % (q, heh),
           '#',
           '# just below the ceiling 0.5/(0.5 + He/H) = %.9f that'
           % ceiling,
           '# input_read.f90 enforces.  The base state is the scalar',
           '# family\'s, so the only change against the atomic scalar case',
           '# at the same He/H is the chemistry.',
           '#',
           '# p_base is the level of the run; "Base BC: pressure 1.0" in',
           '# input.inp names the same level, and the density key is left',
           '# out beside it.',
           'T_base    %s' % T_BASE_K,
           'r_base    %s' % R_BASE_RJ,
           'HeH_base  %s' % heh,
           'Kzz_base  %s' % kzz,
           'q_H2_base %.17g' % q,
           'p_base    %s' % P_BASE_BAR]
    return '\n'.join(out) + '\n'


# --------------------------------------------------------------------------
# the stored profiles
# --------------------------------------------------------------------------

def profile_match_HeH(path, p_match=1.0e-6):
    """He/H at the matching level of a stored column: the `X_He` column of
    the row whose pressure is p_match.  The profile runs deep to shallow and
    carries the matching level exactly.  The column is found BY NAME on the
    `# columns:` header, as the code that reads the profile does
    (`lap_element_ratio_at_match` takes "X_<El>"), so a schema that gains a
    column cannot silently move the value."""
    best, best_d, j_he = None, None, None
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith('#'):
                word = line.split()
                if len(word) > 1 and word[1].rstrip(':') == 'columns':
                    names = [w.rstrip(':') for w in word[2:]]
                    if 'X_He' not in names:
                        raise SystemExit('%s names no X_He column' % path)
                    j_he = names.index('X_He')
                continue
            if not line.strip():
                continue
            if j_he is None:
                raise SystemExit('%s carries data rows before its'
                                 ' "# columns:" header' % path)
            word = line.split()
            p, x_he = float(word[0]), float(word[j_he])
            d = abs(p/p_match - 1.0)
            if best_d is None or d < best_d:
                best, best_d = x_he, d
    if best is None:
        raise SystemExit('%s carries no data rows' % path)
    return best


# --------------------------------------------------------------------------
# the closure configuration
# --------------------------------------------------------------------------

def closure_json(group, heh, case_dir):
    """`element_flux_closure.py --config`: the fixed part of both command
    lines for one rung.  Only the trial fluxes and the iteration index vary
    from one iteration to the next, and the driver supplies those.

    The wind of each iteration is solved by the partitioned stationary route
    (`input_keys`), not by the coupled JFNK the driver defaults to: on this
    planet the coupled solve stalls on the base layer's energy row and the
    well-balanced partitioned route converges (MODELS.md section 6)."""
    sed = os.path.join(LHS, 'sed', group.spectrum_file)
    adapter_args = [
        '--mp', PLANET['mass_MJ'],
        '--r-ref', PLANET['radius_RJ'],
        '--p-ref', '1.0',
        '--p-match', '1.0e-6',
        '--climate',
        '--climate-p-deep', '20.0',
        '--boa-pressure-factor', '1.0',
        '--stellar-flux', sed,
        '--flux-at-planet',
        '--wavelength-unit', 'A',
        '--toa', '1.0e-2',
        '--atoms', 'H,He,N,O,C',
        '--abundances', 'He=%s' % heh,
        '--kzz-const', '1.0e9',
    ]
    # Photochem is a compiled extension and is installed into this
    # interpreter; it is named here so the rung reproduces on the build it
    # was run on.  The adapter runs with cwd = the iteration directory, which
    # also keeps the repository root's photochem_clima_data/ from shadowing
    # the installed package.
    body = [
        '{',
        '  "comment": "Elemental-flux closure of %s, starting reservoir'
        ' He/H = %s. The starting trial fluxes are the driver\'s'
        ' --phi0-H / --phi0-He, not this file.",' % (group.name, heh),
        '  "python": "/opt/miniconda3/bin/python3",',
        '  "adapter": "%s",' % os.path.join(EXHALE, 'src', 'utils',
                                            'photochem_to_lower_profile.py'),
        '  "adapter_run_dir": ".",',
        '  "adapter_args": [',
    ]
    body.append('    ' + ', '.join('"%s"' % a for a in adapter_args))
    body += [
        '  ],',
        '  "p_top_bar": null,',
        '  "exhale_bin": "%s",' % os.path.join(EXHALE, 'EXHALE.x'),
        '  "input_template": "%s",'
        % os.path.join(case_dir, 'input_template.inp'),
        '  "omp_num_threads": 8,',
        '  "resid_tol": "%s",' % RESID_TOL['photochem'],
        '',
        '  "input_keys_comment": "The partitioned stationary route of'
        ' MODELS.md section 6, the same one the prescribed atomic cases are'
        ' solved by. CFL is inert in the solve, which takes no time step,'
        ' and holds the post-processing pass to the state it reads.",',
        '  "input_keys": ["Well balanced: True", "Restart intent:'
        ' stationary", "Solver: Newton",',
        '                  "Secondary_ionization: Immediate",'
        ' "CFL: 1.0e-12"],',
        '',
        '  "exhale_env": {"EXHALE_PTC_DTAU0": "1.0"}',
        '}',
    ]
    return '\n'.join(body) + '\n'


# --------------------------------------------------------------------------
# writing one case
# --------------------------------------------------------------------------

def write_if_absent(path, text, force, report):
    if os.path.isfile(path) and not force:
        with open(path) as fh:
            if fh.read() == text:
                return
    with open(path, 'w') as fh:
        fh.write(text)
    report.append(os.path.basename(path))


def make_case(group, heh, force):
    case_dir = os.path.join(HERE, group.name, 'HeH%s' % heh)
    solved = os.path.join(case_dir, 'output', 'Hydro_ioniz.txt')
    if os.path.isfile(solved) and not force:
        return 'kept', case_dir
    os.makedirs(case_dir, exist_ok=True)
    written = []

    if group.is_closure:
        # The driver writes each iteration's own input.inp from this template
        # and runs the chemistry in k00/, k01/ ... beside it.
        write_if_absent(os.path.join(case_dir, 'input_template.inp'),
                        input_inp(group, heh, for_closure_template=True),
                        force, written)
        write_if_absent(os.path.join(case_dir, 'closure.json'),
                        closure_json(group, heh, case_dir),
                        force, written)
    else:
        heh_value = heh
        if group.uses_profile:
            # The profile owns the reservoir; write the value it actually
            # carries at the matching level so the two do not disagree on
            # paper.
            src = profile_source(heh)
            if src is None:
                return 'skipped (its closure rung has not converged yet)', case_dir
            heh_value = '%.9f' % profile_match_HeH(src)
        write_if_absent(os.path.join(case_dir, 'input.inp'),
                        input_inp(group, heh, heh_value=heh_value),
                        force, written)

        if group.uses_profile:
            src = profile_source(heh)
            dst = os.path.join(case_dir, 'lower_atmosphere_profile.dat')
            if force or not os.path.isfile(dst):
                shutil.copyfile(src, dst)
                written.append(os.path.basename(dst))
        elif group.boundary == 'scalarCNO':
            write_if_absent(os.path.join(case_dir, 'base.inp'),
                            base_inp_scalar_cno(heh), force, written)
        elif group.chemistry == 'molecular':
            write_if_absent(os.path.join(case_dir, 'base.inp'),
                            base_inp_molecular(heh, KZZ_BASE), force, written)

    return ('wrote ' + ', '.join(written)) if written else 'unchanged', case_dir


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('--force', action='store_true',
                    help='overwrite a case that already carries a solution')
    ap.add_argument('--only', default=None,
                    help='act only on groups whose name contains this')
    ap.add_argument('--list', action='store_true',
                    help='print the table and the case count, write nothing')
    a = ap.parse_args()

    total = 0
    closure = 0
    for name, ladder in GROUPS:
        g = Group(name)
        if a.only and a.only not in name:
            continue
        total += len(ladder)
        if g.is_closure:
            closure += len(ladder)
        if a.list:
            print('%-46s %-9s %-12s %2d  %s'
                  % (name, ('closure' if g.is_closure else 'prescribed'),
                     g.spectrum_file.replace('lhs1140_sed_', '')
                     .replace('.txt', ''), len(ladder), ' '.join(ladder)))
            continue
        for heh in ladder:
            what, case_dir = make_case(g, heh, a.force)
            print('%-60s %s' % (os.path.relpath(case_dir, HERE), what))

    print('')
    print('%d cases: %d prescribed composition, %d flux-closure rungs'
          % (total, total - closure, closure))
    return 0


if __name__ == '__main__':
    sys.exit(main())
