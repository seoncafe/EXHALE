#!/usr/bin/env python3
"""Convert a VULCAN photochemistry output (.vul pickle) into an EXHALE
`base.inp` handoff file (Tier-3 upgrade of run_lower.py; see
docs/lower_atmosphere_coupling.*).

What it does
------------
Reads the VULCAN output (T(p), composition mixing ratios vs pressure) and
1. extracts, at the requested handoff pressure (default 1 ubar):
   T_base, q_H2, q_H, q_He and the mean molecular weight -- with the
   PHOTOCHEMICAL H2/H dissociation state (the physics the analytic
   chemical-equilibrium column of run_lower.py lacks at Teq ~ 1000-2000 K);
2. integrates the hypsometric relation from the 1-bar radius upward using
   VULCAN's own mu(p) and T(p) profiles to get the base radius r_base;
3. writes base.inp: the keys run_lower.py also writes (T_base, r_base,
   HeH_base, Kzz_base) plus the photochemical H2 mixing ratio q_H2_base and
   the handoff level p_base, which EXHALE uses in place of its
   chemical-equilibrium H2 fit (docs/base_composition_handoff_plan.md), and
   the elemental reservoirs `HeH_base` and `<El>_H_base` for the elements
   VULCAN's networks carry (C, N, O, S), summed over every carrier.
   The molecular mixing ratios stay comments: nothing in EXHALE consumes
   them, so they are diagnostic metadata (docs/input_schema.md section 2c).

What it does NOT provide: atomic-metal release fractions (Na/Mg/Ca/Fe).
VULCAN's networks are H/C/N/O(/S) -- metal/alkali chemistry and condensation
are outside its scope (that part of the Lavvas model has no public analogue);
metals.inp stays user-supplied for those elements.

Usage:
  python3 vulcan_to_base.py <output.vul> <run_dir> --mp <M_J> --r1bar <R_J>
                            [--pbase 1e-6] [--kzz 1e9]
"""
import argparse, math, os, pickle, re, sys

# Same constants as parameters.f90 and lower_profile_schema.py.  MAMU is the
# atomic mass unit: `mu_profile` below is a mean molecular weight in amu, built
# from the species masses VULCAN uses (H2 = 2.016, He = 4.003, ...), so the mass
# of one particle is mu * m_amu.  Multiplying it by the hydrogen ATOM mass
# instead, as this file did until 2026-08-27, makes every particle 0.8% heavier
# and the scale height 0.8% shorter (G was 6.67259e-8 until 2026-08-19).
KB, MAMU, G = 1.380649e-16, 1.66053906660e-24, 6.67430e-8
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
BAR = 1.0e6  # dyn/cm^2


def load_vul(path):
    with open(path, 'rb') as f:
        data = pickle.load(f)
    species = list(data['variable']['species'])
    ymix = data['variable']['ymix']          # (nz, nsp) mixing ratios
    p = data['atm']['pco'] / BAR             # dyn/cm^2 -> bar (increasing up?)
    T = data['atm']['Tco']
    return species, ymix, p, T


def profile(species, ymix, name):
    if name in species:
        return ymix[:, species.index(name)]
    return None


ELEMENT_RE = re.compile(r'([A-Z][a-z]?)(\d*)')


def formula_elements(name):
    """Element counts of a VULCAN species label, e.g. C2H2 -> {C:2, H:2}.

    VULCAN decorates labels with an ionization sign (`H3+`), an excited-state
    suffix (`CH2_1`) or the bare electron (`e`); none of them change the
    nuclei, so they are stripped before the formula is parsed.
    """
    core = name.split('_')[0].rstrip('+-')
    if core in ('e', ''):
        return {}
    out = {}
    pos, n = 0, len(core)
    while pos < n:
        m = ELEMENT_RE.match(core, pos)
        if m is None or m.start() != pos:      # label is not a formula
            return {}
        out[m.group(1)] = out.get(m.group(1), 0) + int(m.group(2) or 1)
        pos = m.end()
    return out


def element_ratios(species, ymix):
    """El/H nuclei ratio profiles, summed over every species carrying El.

    Photochemistry moves nuclei between molecules without creating or
    destroying them, so these ratios are the elemental reservoirs the wind
    inherits -- the quantity EXHALE's `<El>_H_base` keys and `melem_ab` hold.
    Counting every carrier matters for He/H as well: H bound in H2O, CH4 and
    NH3 belongs in the hydrogen denominator.
    """
    counts = [formula_elements(sp) for sp in species]
    nuc = {}
    for j, c in enumerate(counts):
        for el, k in c.items():
            nuc[el] = nuc.get(el, 0.0) + k * ymix[:, j]
    nH = nuc.get('H')
    if nH is None:
        raise SystemExit('no hydrogen-bearing species in the VULCAN output')
    return {el: y / nH for el, y in nuc.items() if el != 'H'}


def mu_profile(species, ymix):
    """Mean molecular weight [amu] from all species (VULCAN compo masses)."""
    # minimal mass table for the dominant species; everything else ~ trace
    mass = {'H2': 2.016, 'H': 1.008, 'He': 4.003, 'H2O': 18.02, 'CH4': 16.04,
            'CO': 28.01, 'CO2': 44.01, 'N2': 28.01, 'NH3': 17.03, 'C2H2': 26.04,
            'HCN': 27.03, 'H2S': 34.08, 'S': 32.06, 'SO2': 64.06, 'e': 0.0}
    mu = 0.0
    for sp, m in mass.items():
        y = profile(species, ymix, sp)
        if y is not None:
            mu = mu + m * y
    return mu


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('vulfile')
    ap.add_argument('run_dir')
    ap.add_argument('--mp', type=float, required=True, help='planet mass [M_J]')
    ap.add_argument('--r1bar', type=float, required=True,
                    help='radius of the 1-bar level [R_J] (transit radius)')
    ap.add_argument('--pbase', type=float, default=1e-6, help='handoff p [bar]')
    ap.add_argument('--kzz', type=float, default=1e9, help='K_zz at base [cm2/s]')
    a = ap.parse_args()

    species, ymix, p, T = load_vul(a.vulfile)
    mu = mu_profile(species, ymix)

    # sort by decreasing pressure (deep -> top) for the upward integration
    import numpy as np
    idx = np.argsort(p)[::-1]
    p, T, mu = p[idx], T[idx], mu[idx]
    ymix = ymix[idx]

    def interp(y, pq):
        return float(np.interp(math.log(pq), np.log(p[::-1]), y[::-1]))

    # hypsometric integration from 1 bar to pbase with VULCAN mu(p), T(p)
    Mp = a.mp * MJ
    r = a.r1bar * RJ
    lp1, lp0 = math.log(1.0), math.log(a.pbase)
    n = 2000
    dx = (lp0 - lp1) / n
    x = lp1
    for _ in range(n):
        def f(r_, lpq):
            pq = math.exp(lpq)
            return -KB * interp(T, pq) / (interp(mu, pq) * MAMU * (G * Mp / r_**2))
        k1 = f(r, x); k2 = f(r + 0.5*dx*k1, x + 0.5*dx)
        k3 = f(r + 0.5*dx*k2, x + 0.5*dx); k4 = f(r + dx*k3, x + dx)
        r += dx * (k1 + 2*k2 + 2*k3 + k4) / 6.0
        x += dx

    Tb  = interp(T, a.pbase)
    q2  = interp(profile(species, ymix, 'H2'), a.pbase)
    qh  = interp(profile(species, ymix, 'H'),  a.pbase)
    qhe = interp(profile(species, ymix, 'He'), a.pbase)
    mub = interp(mu, a.pbase)
    # Elemental reservoirs at the base: El/H nuclei ratios summed over every
    # carrier (H bound in H2O/CH4/NH3 counts as hydrogen), so they are the
    # conserved quantities EXHALE's melem_ab holds.
    elrat = element_ratios(species, ymix)
    heh = interp(elrat['He'], a.pbase) if 'He' in elrat else 0.0

    out = os.path.join(a.run_dir, 'base.inp')
    with open(out, 'w') as f:
        f.write('# base.inp -- lower-atmosphere handoff for EXHALE\n')
        f.write('# written by vulcan_to_base.py from %s\n' % os.path.basename(a.vulfile))
        f.write('# PHOTOCHEMICAL base state (VULCAN): q_H2=%.4e q_H=%.4e q_He=%.4f mu=%.3f\n'
                % (q2, qh, qhe, mub))
        for sp in ('H2O', 'CO', 'CH4', 'CO2', 'NH3', 'HCN'):
            y = profile(species, ymix, sp)
            if y is not None:
                f.write('# q_%s = %.3e\n' % (sp, interp(y, a.pbase)))
        if q2 > 0.3:
            f.write('# NOTE: molecular base persists photochemically -- Tier-2\n'
                    '#       (Molecular chemistry: True) is required for this planet.\n')
        f.write('T_base    %.2f\n' % Tb)
        f.write('r_base    %.5f\n' % (r / RJ))
        f.write('HeH_base  %.6f\n' % heh)
        f.write('Kzz_base  %.3e\n' % a.kzz)
        # Photochemical H2 partition and the level it refers to. EXHALE uses
        # q_H2_base in place of its chemical-equilibrium fit in the
        # molecular-base particle count, and evaluates that fit at p_base.
        f.write('q_H2_base %.6e\n' % q2)
        f.write('p_base    %.3e\n' % a.pbase)
        # Elemental reservoirs, for the elements EXHALE solves and VULCAN's
        # networks carry. They override metals.inp for those elements.
        # Mg/Si/Ca/Na/K/Fe are NOT written: metal and alkali chemistry (and
        # condensation) is outside VULCAN's scope, so metals.inp still
        # supplies them.
        f.write('# elemental reservoirs at p_base (all carriers counted)\n')
        for el in ('C', 'N', 'O', 'S'):
            if el in elrat:
                x = interp(elrat[el], a.pbase)
                if x > 0.0:
                    f.write('%-9s %.6e\n' % (el + '_H_base', x))
    print('wrote %s:  r(%.0e bar) = %.4f R_J,  T_base = %.1f K,  q_H2 = %.3e, q_H = %.3e'
          % (out, a.pbase, r / RJ, Tb, q2, qh))


if __name__ == '__main__':
    main()
