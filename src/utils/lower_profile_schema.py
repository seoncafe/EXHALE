#!/usr/bin/env python3
"""The lower-atmosphere profile schema, shared by the two adapters.

The schema fixes one file format and one fingerprint for the handoff, and
asks for two producers of it:
`photochem_to_lower_profile.py` (production) and `vulcan_to_lower_profile.py`
(the cross-check), behind the same command line.  Everything that is a
property of the FORMAT rather than of either chemistry code lives here, so
that the two adapters cannot drift apart:

  - the command-line options both take (`add_common_arguments`);
  - elemental accounting: El/H nuclei ratios summed over every carrier, and
    the mean molecular weight from the same counts;
  - the hydrostatic radius, integrated from a stated reference level;
  - level insertion, so that the matching pressure is a node of the table and
    the reader returns its value bit for bit rather than interpolating;
  - the fingerprint `solution_id` and the two files that carry it.

Units are CGS throughout, except the two the schema states otherwise: `p` in
bar and `r` in R_J.

Physics note on the element columns.  Photochemistry moves nuclei between
molecules without creating or destroying them, so the El/H nuclei ratio, and
not any molecule's mixing ratio, is the conserved quantity the wind inherits
and the only one `melem_ab` can hold.  Every carrier is counted, hydrogen
bound in H2O/CH4/NH3 included (the lesson `vulcan_to_base.py` already
records).  Condensed carriers are counted only when the caller asks: the gas
EXHALE's base inherits is the gas phase, and the difference between the two
sums IS the cold trap, which section 5 of the design reports separately.
"""
import hashlib
import os
import re
import sys

import numpy as np

# ---- constants ------------------------------------------------------- #
# The Jupiter mass is EXHALE's own value (parameters.f90), so a mass written
# here means the same mass the code reads back; the radius comes from the one
# Python definition of it, imported below.
KB = 1.380649e-16
MAMU = 1.66053906660e-24
GNEWT = 6.67430e-8
MJ = 1.898e30

# The Jupiter radius is defined once, in examples/exhale_io.py (RJ_CM), and
# imported here rather than written down again: a length this file computes
# from `Planet radius [R_J]` has to be the length EXHALE's parameters.f90
# means by the same key.  The directory is APPENDED to sys.path, never
# prepended, so that nothing in it can shadow a standard-library module.
_EXAMPLES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'examples')
if _EXAMPLES_DIR not in sys.path:
    sys.path.append(_EXAMPLES_DIR)
from exhale_io import RJ_CM                                # noqa: E402
RJ = RJ_CM
RSUN = 6.957e10
AU = 1.495979e13
BAR = 1.0e6                      # dyn/cm^2

# IUPAC standard atomic weights, for the elements the gas-giant networks
# carry.  They enter the mean molecular weight and hence the scale height.
ATOMIC_WEIGHT = {
    'H': 1.008, 'He': 4.002602, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'S': 32.06, 'Na': 22.98976928, 'Mg': 24.305, 'Si': 28.085,
    'K': 39.0983, 'Ca': 40.078, 'Fe': 55.845, 'Ar': 39.948, 'Cl': 35.45,
    'P': 30.973761998, 'Ti': 47.867, 'V': 50.9415,
}

# The molecular mixing ratios the profile carries as diagnostics: what the
# EXHALE molecular network knows, plus what the cold-trap statement needs.
DIAGNOSTIC_MOLECULES = ('H2O', 'CO', 'CO2', 'CH4', 'NH3', 'HCN', 'N2', 'OH',
                        'H2S', 'C2H2')

# The elements a profile carries an X_<El> column for, in this order.  EXHALE
# reads them BY NAME (`lap_element_ratio_at_match`), so the order is only for
# the reader's eye; it is fixed here so two adapters produce diffable files.
PROFILE_ELEMENTS = ('He', 'C', 'N', 'O', 'S', 'Na', 'Mg', 'Si', 'K', 'Ca',
                    'Fe')

REQUIRED_COLUMNS = ('p', 'r', 'T', 'n_tot', 'rho', 'Kzz', 'q_H2', 'q_H',
                    'X_He')


# ---- provenance ------------------------------------------------------- #
def sha256_file(path):
    """sha256 of a file, or 'absent' when there is no file to hash."""
    if path is None or not os.path.exists(path):
        return 'absent'
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def solution_fingerprint(entries):
    """The `solution_id` of section 2.1: sha256 over the lower model's
    configuration, its mechanism and thermodynamic data, its stellar flux
    file, the elemental abundance vector and the trial fluxes of this closure
    iteration.  `entries` is a dict; it is hashed in sorted key order so that
    the same configuration gives the same id whatever order it was built in.

    The iteration INDEX is deliberately not part of it: two iterations that
    were handed the same configuration and the same trial fluxes are the same
    solution, and the fingerprint is what says so.
    """
    h = hashlib.sha256()
    for key in sorted(entries):
        h.update(key.encode())
        h.update(b'\x00')
        h.update(repr(entries[key]).encode())
        h.update(b'\x00')
    return h.hexdigest()


# ---- elemental accounting --------------------------------------------- #
_ELEMENT_RE = re.compile(r'([A-Z][a-z]?)(\d*)')


def formula_elements(name):
    """Element counts of a chemistry-code species label, e.g. C2H2 -> {C:2,
    H:2}.  VULCAN decorates labels with an ionization sign (`H3+`), an
    excited-state suffix (`CH2_1`) or the bare electron (`e`); none of them
    change the nuclei, so they are stripped before the formula is parsed.  A
    label that is not a formula returns no nuclei rather than a wrong count.

    This is the VULCAN route.  Photochem states the composition matrix
    itself (`dat.species_composition`) and that is used instead, because
    labels like `O1D` and `N2D` are excited states, not formulas, and would
    parse into nuclei that do not exist.
    """
    core = name.split('_')[0].rstrip('+-')
    if core in ('e', '', 'M', 'hv'):
        return {}
    out = {}
    pos, n = 0, len(core)
    while pos < n:
        m = _ELEMENT_RE.match(core, pos)
        if m is None or m.start() != pos:
            return {}
        el = m.group(1)
        if el not in ATOMIC_WEIGHT:
            return {}
        out[el] = out.get(el, 0) + int(m.group(2) or 1)
        pos = m.end()
    return out


def element_ratios(mixing, counts):
    """El/H nuclei ratio profiles, summed over every carrier in `mixing`.

    `mixing` maps a species name to its volume mixing ratio profile;
    `counts` maps the same names to their element counts.  Species absent
    from `counts` carry no nuclei and are skipped.
    """
    nuc = {}
    for sp, y in mixing.items():
        for el, k in counts.get(sp, {}).items():
            nuc[el] = nuc.get(el, 0.0) + k*np.asarray(y, dtype=float)
    nH = nuc.get('H')
    if nH is None:
        raise SystemExit('no hydrogen-bearing species in the chemistry '
                         'solution: nothing to state an El/H ratio against')
    return {el: y/nH for el, y in nuc.items() if el != 'H'}, nH


def mean_molecular_weight(mixing, counts):
    """Mean molecular weight [amu] from the same element counts, so that the
    mass and the nuclei budget of the file are one statement and not two."""
    mu = 0.0
    for sp, y in mixing.items():
        m = sum(ATOMIC_WEIGHT[el]*k for el, k in counts.get(sp, {}).items())
        if m > 0.0:
            mu = mu + m*np.asarray(y, dtype=float)
    return mu


# ---- hydrostatic radius ------------------------------------------------ #
def hydrostatic_radius(p_bar, T, mu, mp_MJ, r_ref_RJ, p_ref_bar, nsub=64):
    """Integrate dr/dln p = -k T r^2 / (mu m_amu G M_p) from the reference
    level, outward and inward, and return r [R_J] at every level of `p_bar`.

    Gravity varies with radius, as it does in `make_example_profile.py`; the
    integration is RK4 in ln p, on a grid refined `nsub` times per level, so
    the result does not depend on how the chemistry code chose to space its
    levels.  The reference level is where the planet radius is stated: for a
    transit radius that is 1 bar, the convention `vulcan_to_base.py` uses.

    The mean molecular weight enters as mu * m_amu, the mass of a particle of
    mu atomic mass units.  (`vulcan_to_base.py` multiplies the same mu by the
    hydrogen ATOM mass instead, which is 0.8% heavier; that discrepancy is
    reported, not carried over.)
    """
    p_bar = np.asarray(p_bar, dtype=float)
    order = np.argsort(-p_bar)                # deep -> shallow
    x_src = np.log(p_bar[order])
    T_src = np.asarray(T, dtype=float)[order]
    mu_src = np.asarray(mu, dtype=float)[order]

    x_ref = np.log(p_ref_bar)
    if not (x_src[0] >= x_ref >= x_src[-1]):
        raise SystemExit('the reference pressure %.3e bar is outside the '
                         'chemistry solution (%.3e .. %.3e bar): state a '
                         'reference level the solution covers'
                         % (p_ref_bar, p_bar.max(), p_bar.min()))

    def interp(arr, x):
        return float(np.interp(x, x_src[::-1], arr[::-1]))

    mp = mp_MJ*MJ

    def drdx(r, x):
        pq = x
        return -KB*interp(T_src, pq)*r*r/(interp(mu_src, pq)*MAMU*GNEWT*mp)

    def march(x_end):
        """RK4 from the reference level to x_end, recording (x, r)."""
        n = max(4, int(nsub*len(x_src)))
        xs = np.linspace(x_ref, x_end, n + 1)
        rs = np.empty_like(xs)
        rs[0] = r_ref_RJ*RJ
        for i in range(n):
            h = xs[i+1] - xs[i]
            r0 = rs[i]
            k1 = drdx(r0, xs[i])
            k2 = drdx(r0 + 0.5*h*k1, xs[i] + 0.5*h)
            k3 = drdx(r0 + 0.5*h*k2, xs[i] + 0.5*h)
            k4 = drdx(r0 + h*k3, xs[i] + h)
            rs[i+1] = r0 + h*(k1 + 2*k2 + 2*k3 + k4)/6.0
            if rs[i+1] <= 0.0:
                raise SystemExit('the hydrostatic integration reached a '
                                 'non-positive radius: the column is not '
                                 'bound at the reference level given')
        return xs, rs

    x_up, r_up = march(x_src[-1])             # outward, ln p falling
    x_dn, r_dn = march(x_src[0])              # inward, ln p rising
    x_all = np.concatenate([x_dn[::-1], x_up[1:]])
    r_all = np.concatenate([r_dn[::-1], r_up[1:]])
    o = np.argsort(x_all)
    r_of_x = np.interp(np.log(p_bar), x_all[o], r_all[o])
    return r_of_x/RJ


# ---- the table --------------------------------------------------------- #
def eddy_diffusion_coefficient(p_bar, args):
    """The eddy coefficient the handoff states, on the levels `p_bar` [bar].

    Two forms, and neither is a default: with neither stated this returns
    None and the caller keeps whatever K_zz(p) the underlying solution
    carried.

    * `--kzz-const K` -- one value at every level.
    * `--kzz-power ALPHA` with `--kzz-ref K_ref` and `--kzz-ref-bar p_ref` --

          K_zz(p) = K_ref (p/p_ref)^(-ALPHA),

      the saturated gravity-wave form.  An upward-propagating internal wave
      conserves its energy flux, so its velocity amplitude grows as
      rho^(-1/2) ~ p^(-1/2) until it breaks, and the turbulence that holds
      it at saturation is parameterized as an eddy coefficient rising with
      height (Lindzen 1981).  The exoplanet literature writes it as
      K_zz = K_0 P_bar^(-1/2) (Parmentier et al. 2013, HD 209458 b) or
      P_bar^(-0.4) (Charnay et al. 2015, GJ 1214 b), which is this
      expression with p_ref = 1 bar.  The theory that produces the slope
      also bounds it: the wave field is exhausted at some height and the
      coefficient turns over there (docs/eddy_diffusion_kzz.tex,
      "Breaking gravity waves"), so extrapolating far above the level the
      amplitude was set at reads as an upper bound.
    """
    p_bar = np.asarray(p_bar, dtype=float)
    if getattr(args, 'kzz_power', None) is not None:
        if args.kzz_const is not None:
            refuse('--kzz-const and --kzz-power state the eddy coefficient '
                   'two different ways: give one of them')
        if args.kzz_ref is None:
            refuse('--kzz-power needs its amplitude --kzz-ref [cm^2/s], the '
                   'eddy coefficient at --kzz-ref-bar')
        if not (args.kzz_ref > 0.0 and args.kzz_ref_bar > 0.0):
            refuse('--kzz-ref and --kzz-ref-bar must both be positive')
        return args.kzz_ref*(p_bar/args.kzz_ref_bar)**(-args.kzz_power)
    if args.kzz_const is not None:
        return np.full(p_bar.shape, float(args.kzz_const))
    return None


def insert_level(cols, p_target):
    """Return `cols` with an exact node at `p_target` [bar].

    Every intensive quantity is interpolated linearly in log p, the schema's
    own convention (section 2.3) and the reader's.  A node at exactly the
    matching pressure is what makes `value_at_pressure` return the table's
    own number rather than an interpolation of it, so the base state EXHALE
    builds is the state this file states.

    `n_tot` and `rho` are not free columns: the gas is ideal at every node
    the solution was computed on, `n = p/(k_B T)` and `rho = n mu m_u`, and
    interpolating the two densities alongside p and T breaks that identity at
    the inserted node -- the one level a reader takes its base state from.
    They are therefore recomputed here from the interpolated pressure and
    temperature, with the mean molecular weight `rho/n` (which is intensive
    and is the quantity to interpolate) carrying the composition.
    """
    p = cols['p']
    if np.any(p == p_target):
        return cols
    if not (p[0] > p_target > p[-1]):
        return cols
    j = int(np.searchsorted(-p, -p_target))
    x0, x1 = np.log(p[j-1]), np.log(p[j])
    w = (np.log(p_target) - x0)/(x1 - x0)
    out = {}
    for name, y in cols.items():
        v = y[j-1] + w*(y[j] - y[j-1])
        out[name] = np.insert(y, j, v)
    out['p'][j] = p_target
    if 'n_tot' in out and 'T' in out and out['T'][j] > 0.0:
        mu_m = None
        if 'rho' in out and cols['n_tot'][j-1] > 0.0 and cols['n_tot'][j] > 0.0:
            m0 = cols['rho'][j-1]/cols['n_tot'][j-1]
            m1 = cols['rho'][j]/cols['n_tot'][j]
            mu_m = m0 + w*(m1 - m0)
        out['n_tot'][j] = p_target*BAR/(KB*out['T'][j])
        if mu_m is not None:
            out['rho'][j] = out['n_tot'][j]*mu_m
    return out


def clip_table(cols, p_deep_bar, p_top_bar):
    """Keep the levels inside [p_top, p_deep], deep to shallow."""
    p = cols['p']
    keep = (p <= p_deep_bar*(1.0 + 1.0e-12)) & (p >= p_top_bar*(1.0 - 1.0e-12))
    return {name: y[keep] for name, y in cols.items()}


def sort_deep_to_shallow(cols):
    o = np.argsort(-cols['p'])
    return {name: np.asarray(y, dtype=float)[o] for name, y in cols.items()}


def strictly_decreasing(cols):
    """Drop levels that repeat a pressure: the reader refuses a table that is
    not strictly decreasing, and a chemistry grid can carry a duplicate."""
    p = cols['p']
    keep = np.ones(len(p), dtype=bool)
    last = None
    for i in range(len(p)):
        if last is not None and p[i] >= last:
            keep[i] = False
        else:
            last = p[i]
    return {name: y[keep] for name, y in cols.items()}


# ---- the two files ----------------------------------------------------- #
def write_profile(path, cols, header, column_order=None):
    """Write the profile file of section 2.

    `header` is a dict of the machine-readable header fields; `cols` maps a
    column name to its profile, deep to shallow.  Every number is written at
    full double precision, so a round trip through the file reproduces the
    solution rather than a rounding of it.
    """
    names = list(column_order or cols)
    missing = [c for c in REQUIRED_COLUMNS if c not in names]
    if missing:
        raise SystemExit('the profile is missing required columns: %s'
                         % ', '.join(missing))
    n = len(cols[names[0]])
    with open(path, 'w') as f:
        f.write('# EXHALE lower-atmosphere profile\n')
        for key in ('solution_id', 'source_code', 'source_version',
                    'mechanism', 'stellar_flux'):
            f.write('# %s %s\n' % (key, header[key]))
        for key in ('p_match_bar', 'p_top_bar', 'p_deep_bar',
                    'trial_flux_H', 'trial_flux_He'):
            f.write('# %s %.17E\n' % (key, header[key]))
        # Optional, and written only when the producer measured them: the
        # elemental fluxes the lower model itself carries across its top,
        # summed over every carrier by nuclei count (g/s, outward positive).
        # `trial_flux_*` is what the closure ASKED for; these are what the
        # solution DELIVERS, so the driver can check that the boundary
        # condition took before it forms a residual from it.  The EXHALE
        # reader ignores header keys it has no consumer for
        # (`lower_atmosphere_profile.f90`, the `case default` of its header
        # parser), so the pair is additive to the format.
        for key in ('measured_flux_H', 'measured_flux_He'):
            if header.get(key) is not None:
                f.write('# %s %.17E\n' % (key, header[key]))
        f.write('# iteration %d\n' % header['iteration'])
        f.write('# reached_steady_state %s\n'
                % ('T' if header['reached_steady_state'] else 'F'))
        # The reader holds 1000 characters (`lap_notes`), which is what a
        # climate solve's record needs; the file itself is not truncated
        # more tightly than the reader is.
        f.write('# notes %s\n' % header['notes'].replace('\n', ' ')[:1000])
        f.write('# units p[bar] r[R_J] T[K] n_tot[cm^-3] rho[g/cm^3]'
                ' Kzz[cm^2/s] q_*[-] X_*[El/H nuclei] F_*[g/s]\n')
        f.write('# columns: %s\n' % ' '.join(names))
        for i in range(n):
            f.write(' '.join('%25.17E' % cols[nm][i] for nm in names) + '\n')
    return path


def write_base_inp(path, header, extra_comments=()):
    """The minimal `base.inp` that pairs with a profile.

    With a profile in use EXHALE refuses every EOS-boundary,
    elemental-reservoir and boundary-constraint key (section 2.4), so this
    file carries no physics at all: it exists to state, in the one provenance
    comment the code parses, that it belongs to the same lower-atmosphere
    solution.  Writing any scalar key here would stop the run, by design.
    """
    with open(path, 'w') as f:
        f.write('# base.inp -- provenance only.\n')
        f.write('# The lower-atmosphere PROFILE is the single source of the'
                ' base state and the\n')
        f.write('# elemental reservoirs; every scalar physics key is refused'
                ' beside it.\n')
        f.write('# solution_id %s\n' % header['solution_id'])
        f.write('# source_code %s\n' % header['source_code'])
        f.write('# source_version %s\n' % header['source_version'])
        f.write('# p_match_bar %.6E\n' % header['p_match_bar'])
        f.write('# iteration %d\n' % header['iteration'])
        for line in extra_comments:
            f.write('# %s\n' % line)
    return path


# ---- the command line both adapters take ------------------------------- #
def add_common_arguments(ap):
    """The options that are properties of the handoff, not of either code."""
    ap.add_argument('run_dir', help='directory the two files are written to')
    ap.add_argument('--mp', type=float, required=True,
                    help='planet mass [M_J]')
    ap.add_argument('--r-ref', type=float, required=True,
                    help='planet radius [R_J] at the reference pressure')
    ap.add_argument('--p-ref', type=float, default=1.0,
                    help='reference pressure of --r-ref [bar] (default 1)')
    ap.add_argument('--p-match', type=float, default=1.0e-6,
                    help='matching pressure [bar]: where EXHALE places its'
                         ' base (default 1 microbar)')
    ap.add_argument('--p-top', type=float, default=None,
                    help='shallowest level to carry [bar]; default is the'
                         ' top of the chemistry solution')
    ap.add_argument('--p-deep', type=float, default=None,
                    help='deepest level to carry [bar]; default is the'
                         ' bottom of the chemistry solution')
    ap.add_argument('--trial-flux-H', type=float, default=0.0,
                    help='elemental H flux of this closure iteration [g/s],'
                         ' outward positive')
    ap.add_argument('--trial-flux-He', type=float, default=0.0,
                    help='the same for helium [g/s]')
    ap.add_argument('--iteration', type=int, default=0,
                    help='closure iteration index k (default 0)')
    ap.add_argument('--kzz-const', type=float, default=None,
                    help='override the eddy coefficient with one constant'
                         ' [cm^2/s]; default is the profile the chemistry'
                         ' solution used')
    ap.add_argument('--kzz-power', type=float, default=None,
                    help='override the eddy coefficient with the saturated'
                         ' gravity-wave power law K_zz(p) = K_ref'
                         ' (p/p_ref)^-ALPHA; needs --kzz-ref, and'
                         ' --kzz-ref-bar fixes p_ref. ALPHA > 0 makes K_zz'
                         ' rise with height (1/2 for Lindzen 1981'
                         ' saturation, the value Parmentier et al. 2013 fit;'
                         ' Charnay et al. 2015 fit 0.4)')
    ap.add_argument('--kzz-ref', type=float, default=None,
                    help='the eddy coefficient at --kzz-ref-bar [cm^2/s],'
                         ' i.e. the amplitude of --kzz-power')
    ap.add_argument('--kzz-ref-bar', type=float, default=1.0,
                    help='the pressure --kzz-ref is stated at [bar]'
                         ' (default 1 bar, the literature normalization)')
    ap.add_argument('--count-condensates', action='store_true',
                    help='count condensed carriers in the El/H ratios;'
                         ' default is the gas phase, which is what the wind'
                         ' inherits')
    ap.add_argument('--abundance-tol', type=float, default=1.0e-10,
                    help='relative tolerance of the deepest-level elemental'
                         ' check (section 4.2). What that level can test is'
                         ' narrower than conservation: Photochem pins every'
                         ' gas species there with bc_type=press, so the level'
                         ' carries its equilibrium initialization unchanged'
                         ' and no photochemical or transport leak reaches it.'
                         ' The check is therefore a guard on this adapter own'
                         ' carrier summation, and it also reads out how far'
                         ' the equilibrium solver got. Where that solver'
                         ' tests each element against its own abundance, the'
                         ' second contribution is round-off: over 68 LHS'
                         ' 1140 b handoffs spanning He/H = 0.0969 to 15 the'
                         ' largest departure is 3.3e-13 and the median'
                         ' 1.7e-14, the 13-significant-figure floor of the'
                         ' carrier summation itself. 1e-10 is the design'
                         ' value and keeps a factor 305 over that floor,'
                         ' which leaves a real miscount detectable: each'
                         ' element sits in one dominant carrier at this'
                         ' level, so dropping one is an O(1) error, and even'
                         ' the smallest -- dropping N2 -- is 1.0e-5 of X_N in'
                         ' the coldest band and 4.0e-3 where the deep column'
                         ' is hot enough to make N2. A solver that scales'
                         ' every elemental residual by the largest elemental'
                         ' abundance does not reach this level: photochem'
                         ' 0.8.4 leaves up to 2.1e-4 in N/H, and a refusal'
                         ' against it is the intended signal that the'
                         ' equilibrium solve was never converged element by'
                         ' element. Reproducing a result made with such a'
                         ' build needs an explicit --abundance-tol.')
    ap.add_argument('--profile-name', default='lower_atmosphere_profile.dat')
    ap.add_argument('--no-base-inp', action='store_true',
                    help='do not write the paired provenance base.inp')
    return ap


def refuse(why):
    """No file is written when the solution cannot be handed over: section
    4.2 makes every one of these a refusal, because a profile that is not a
    steady solution over the overlap has no elemental flux to hand on."""
    raise SystemExit('REFUSED: %s\n  No profile was written.' % why)


def write_handoff(args, cols, header, column_order, element_input=None,
                  element_measured=None):
    """The part of the two adapters that is identical: the checks of section
    4.2 that are properties of the file, then the two files."""
    cols = strictly_decreasing(sort_deep_to_shallow(cols))
    if len(cols['p']) < 2:
        refuse('fewer than two levels survive the requested pressure window')

    p_deep = float(cols['p'][0])
    p_top = float(cols['p'][-1])
    p_match = float(header['p_match_bar'])
    if not (p_deep >= p_match > p_top):
        refuse('the matching pressure %.3e bar is not inside the table '
               '(%.3e .. %.3e bar)' % (p_match, p_deep, p_top))
    cols = insert_level(cols, p_match)
    header['p_deep_bar'] = float(cols['p'][0])
    header['p_top_bar'] = float(cols['p'][-1])
    if header['p_top_bar'] >= p_match:
        refuse('p_top_bar >= p_match_bar: the file stops at or below the '
               'match, so the two models never overlap')

    # The adapter's own carrier-summation check (section 4.2).  It is stated
    # at the deepest level, because that is the only level where the elemental
    # ratios are the input ones: transport and the upper boundary condition
    # move nuclei with height, which is the whole point of the closure.
    #
    # What it can and cannot see.  Photochem pins every gas species at the
    # model bottom with bc_type='press' (gasgiants._initialize_atmosphere), so
    # that level carries the equilibrium initialization unchanged -- measured
    # identical to 13 significant figures before and after the photochemical
    # solve, and independent of the trial escape flux.  A photochemical or
    # transport leak therefore never reaches this check.  What it does test is
    # (i) that this adapter sums every carrier of an element, which is the
    # intent, and (ii) how far the equilibrium solver got.  Where that solver
    # tests each element against its own abundance, (ii) is round-off and the
    # tolerance is set by (i): 3.3e-13 worst and 1.7e-14 median over 68 LHS
    # 1140 b handoff writes, against which the 1e-10 default holds a factor 305.
    # (The 3.5e-13 / 1.9e-14 pair quoted elsewhere is the 38 saved deep
    # equilibrium states re-solved directly, a different set.)
    # Where it instead scales every elemental residual by the largest
    # elemental abundance, (ii) is the solver's own residual, 4072x more
    # visible in N than in He, C or O because deep nitrogen sits in NH3 alone
    # and one excess N rides on three excess H; a refusal then reports that
    # solver rather than a fault here, and reproducing such a run needs an
    # explicit --abundance-tol.
    if element_input and element_measured:
        worst, worst_el = 0.0, ''
        for el, want in element_input.items():
            if el == 'H' or want <= 0.0:
                continue
            got = element_measured.get(el)
            if got is None:
                continue
            dev = abs(got/want - 1.0)
            if dev > worst:
                worst, worst_el = dev, el
            if dev > args.abundance_tol:
                refuse('the %s/H ratio at the deepest level, %.6e, departs '
                       'from the input %.6e by %.2e (tolerance %.1e): the '
                       'element columns do not sum to the input abundances'
                       % (el, got, want, dev, args.abundance_tol))
        print('  elemental conservation at the deepest level: worst '
              'departure %.2e (%s/H), tolerance %.1e'
              % (worst, worst_el or '-', args.abundance_tol))

    os.makedirs(args.run_dir, exist_ok=True)
    prof = os.path.join(args.run_dir, args.profile_name)
    write_profile(prof, cols, header, column_order)
    print('wrote %s: %d levels, %.3e .. %.3e bar'
          % (prof, len(cols['p']), header['p_deep_bar'],
             header['p_top_bar']))
    print('  solution_id %s' % header['solution_id'])
    im = int(np.argmin(np.abs(cols['p'] - p_match)))
    print('  match p = %.4e bar: r = %.6f R_J, T = %.2f K, q_H2 = %.6e, '
          'q_H = %.6e, X_He = %.6e'
          % (cols['p'][im], cols['r'][im], cols['T'][im], cols['q_H2'][im],
             cols['q_H'][im], cols['X_He'][im]))
    if not args.no_base_inp:
        write_base_inp(os.path.join(args.run_dir, 'base.inp'), header)
        print('wrote %s (provenance only)'
              % os.path.join(args.run_dir, 'base.inp'))
    return cols
