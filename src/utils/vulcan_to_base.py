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
3. writes base.inp (same keys as run_lower.py: T_base, r_base, HeH_base,
   Kzz_base), recording the molecular mixing ratios as comments.

What it does NOT provide: atomic-metal release fractions (Na/Mg/Ca/Fe).
VULCAN's networks are H/C/N/O(/S) -- metal/alkali chemistry and condensation
are outside its scope (that part of the Lavvas model has no public analogue);
metals.inp stays user-supplied.

Usage:
  python3 vulcan_to_base.py <output.vul> <run_dir> --mp <M_J> --r1bar <R_J>
                            [--pbase 1e-6] [--kzz 1e9]
"""
import argparse, math, os, pickle, sys

KB, MH, G = 1.380649e-16, 1.6726e-24, 6.67259e-8
RJ, MJ = 6.9911e9, 1.898e30
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
            return -KB * interp(T, pq) / (interp(mu, pq) * MH * (G * Mp / r_**2))
        k1 = f(r, x); k2 = f(r + 0.5*dx*k1, x + 0.5*dx)
        k3 = f(r + 0.5*dx*k2, x + 0.5*dx); k4 = f(r + dx*k3, x + dx)
        r += dx * (k1 + 2*k2 + 2*k3 + k4) / 6.0
        x += dx

    Tb  = interp(T, a.pbase)
    q2  = interp(profile(species, ymix, 'H2'), a.pbase)
    qh  = interp(profile(species, ymix, 'H'),  a.pbase)
    qhe = interp(profile(species, ymix, 'He'), a.pbase)
    mub = interp(mu, a.pbase)
    # elemental He/H ratio at the base (nuclei)
    heh = qhe / max(qh + 2.0*q2, 1e-30)

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
    print('wrote %s:  r(%.0e bar) = %.4f R_J,  T_base = %.1f K,  q_H2 = %.3e, q_H = %.3e'
          % (out, a.pbase, r / RJ, Tb, q2, qh))


if __name__ == '__main__':
    main()
