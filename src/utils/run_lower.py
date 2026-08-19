#!/usr/bin/env python3
"""Tier-3 lower-atmosphere driver: produce a `base.inp` handoff file for EXHALE.

Reads the planet parameters from an EXHALE `input.inp` in the given run
directory, builds a lower/middle-atmosphere column from the 1-bar level to the
1-ubar escape-model base, and writes `base.inp` (consumed by EXHALE's
input_read; see docs/lower_atmosphere_coupling.*).

Temperature models (increasing fidelity):
  --iso              isothermal at Teq (default; = the Fortran Tier-1 column)
  --guillot          Guillot (2010) semi-grey analytic T(p) with an isothermal
                     skin above the tau_IR ~ kappa_IR p / g photosphere
  (--vulcan          upgrade path, NOT implemented here: run the public VULCAN
                     photochemical kinetics code for the composition and a real
                     RC model for T(p), then write the same base.inp keys plus
                     a metals.inp with the atomic-metal release fractions.)

Composition: chemical-equilibrium H2/H/He partition (Visscher fit as in
Koskinen et al. 2022 Eq. 11); the H2->H photochemical dissociation known to
exceed equilibrium at Teq ~ 1000-2000 K is the --vulcan upgrade.

Usage:
  python3 src/utils/run_lower.py <run_dir> --r1bar <R_J> [--iso|--guillot]
                                 [--tint 100] [--kzz 1e9]
Writes <run_dir>/base.inp.
"""
import argparse, math, os, re, sys

# Same constants as parameters.f90 -- this driver integrates the same column as
# lower_column.f90, so the two must not drift. MH is the hydrogen ATOM, because
# it multiplies a dimensionless mean molecular weight (it was the proton mass
# 1.6726e-24 and G was 6.67259e-8 until 2026-08-19).
KB, MH, G = 1.380649e-16, 1.67353284e-24, 6.67430e-8
RJ, MJ = 6.9911e9, 1.898e30


def read_input_inp(run_dir):
    """Pull Mp [g], Teq [K], HeH from an EXHALE input.inp."""
    path = os.path.join(run_dir, 'input.inp')
    txt = open(path).read().splitlines()
    def grab(label, word):
        for ln in txt:
            if label in ln:
                return float(ln.split()[word - 1])
        raise KeyError(label)
    Mp   = grab('Planet mass', 4) * MJ
    Teq  = grab('Equilibrium temperature', 4)
    HeH  = grab('He/H number ratio', 4)
    return Mp, Teq, HeH


def q_h2_equilibrium(p_bar, T):
    """Koskinen+2022 Eq. 11 (Visscher) mixture H2 mixing ratio."""
    u = -23672.0 / T - math.log10(max(p_bar, 1e-30)) + 6.2645
    if u > 30:  return 1.0
    if u < -30: return 0.0
    tenu = 10.0 ** u
    q = (1.9845 + tenu - math.sqrt(tenu * (3.9690 + tenu))) / 2.3670
    return min(max(q, 0.0), 1.0)


def mu_of(p_bar, T, fhe):
    q = q_h2_equilibrium(p_bar, T)
    x2 = min(1.0, 2.0 * q * (1.0 + fhe) / (1.0 + q))
    npart = (1.0 - x2) + 0.5 * x2 + fhe
    return (1.0 + 4.0 * fhe) / npart


def guillot_T(p_bar, Teq, Tint=100.0, kappa_ir=1e-2, gamma=0.4, g_cgs=1e3):
    """Guillot (2010) Eq. 29 semi-grey global-average T(p) profile.
    kappa_ir [cm2/g], gamma = kappa_vis/kappa_ir; standard hot-Jupiter
    defaults (Guillot 2010; Line+2013 typical values)."""
    tau = kappa_ir * (p_bar * 1e6) / g_cgs          # p in dyn/cm^2
    Tirr4 = 4.0 * Teq ** 4                          # Tirr = sqrt(2)*... (mu*=1/sqrt(3) avg)
    f1 = 0.75 * Tint ** 4 * (2.0 / 3.0 + tau)
    x = gamma * tau * math.sqrt(3.0)
    f2 = 0.75 * Teq ** 4 * (2.0 / 3.0 + 1.0 / (gamma * math.sqrt(3.0))
                            + (gamma / math.sqrt(3.0) - 1.0 / (gamma * math.sqrt(3.0)))
                            * math.exp(-x))
    return (f1 + f2) ** 0.25


def integrate_column(Mp, R1bar, Teq, fhe, tmodel, Tint, p_deep=1.0, p_base=1e-6, n=2000):
    """RK4 in x=ln p for dr/dx = -kT/(mu m_H g)."""
    x, dx, r = math.log(p_deep), (math.log(p_base) - math.log(p_deep)) / n, R1bar
    g1 = G * Mp / R1bar ** 2
    def T_of(p):
        return Teq if tmodel == 'iso' else guillot_T(p, Teq, Tint, g_cgs=g1)
    def f(r_, p_):
        T = T_of(p_)
        return -KB * T / (mu_of(p_, T, fhe) * MH * (G * Mp / r_ ** 2))
    for _ in range(n):
        k1 = f(r, math.exp(x))
        k2 = f(r + 0.5 * dx * k1, math.exp(x + 0.5 * dx))
        k3 = f(r + 0.5 * dx * k2, math.exp(x + 0.5 * dx))
        k4 = f(r + dx * k3, math.exp(x + dx))
        r += dx * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        x += dx
    Tb = T_of(p_base)
    q = q_h2_equilibrium(p_base, Tb)
    x2 = min(1.0, 2.0 * q * (1.0 + fhe) / (1.0 + q))
    npart = (1.0 - x2) + 0.5 * x2 + fhe
    return r, Tb, 0.5 * x2 / npart, (1.0 - x2) / npart, fhe / npart


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('run_dir')
    ap.add_argument('--r1bar', type=float, required=True,
                    help='radius of the 1-bar level [R_J] (e.g. transit radius)')
    ap.add_argument('--guillot', action='store_true',
                    help='Guillot (2010) semi-grey T(p) instead of isothermal Teq')
    ap.add_argument('--tint', type=float, default=100.0, help='internal T [K]')
    ap.add_argument('--kzz', type=float, default=1e9, help='K_zz at base [cm2/s]')
    a = ap.parse_args()

    Mp, Teq, HeH = read_input_inp(a.run_dir)
    tmodel = 'guillot' if a.guillot else 'iso'
    rb, Tb, q2, qh, qhe = integrate_column(Mp, a.r1bar * RJ, Teq, HeH, tmodel, a.tint)

    out = os.path.join(a.run_dir, 'base.inp')
    with open(out, 'w') as f:
        f.write('# base.inp -- lower-atmosphere handoff for EXHALE\n')
        f.write('# written by run_lower.py  (T model: %s; R_1bar = %.4f R_J)\n'
                % (tmodel, a.r1bar))
        f.write('# base q_H2 = %.4f  q_H = %.4f  q_He = %.4f  (chem. equilibrium)\n'
                % (q2, qh, qhe))
        if q2 > 0.3:
            f.write('# WARNING: strongly molecular base in equilibrium -- Tier-2\n'
                    '#          molecular physics (or VULCAN photochemistry) needed.\n')
        f.write('T_base    %.2f\n' % Tb)
        f.write('r_base    %.5f\n' % (rb / RJ))
        f.write('HeH_base  %.6f\n' % HeH)
        f.write('Kzz_base  %.3e\n' % a.kzz)
    print('wrote %s:  r(1ubar) = %.4f R_J,  T_base = %.1f K,  q_H2 = %.3f'
          % (out, rb / RJ, Tb, q2))


if __name__ == '__main__':
    main()
