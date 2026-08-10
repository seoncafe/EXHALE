#!/usr/bin/env python3
"""Compare VULCAN and Photochem on HD 189733 b at the pressure EXHALE reads its
base composition from.

Three runs:

  A  VULCAN     NCHO_photo_network      the code the EXHALE Tier-3 pre-step couples to
  B  Photochem  Zahnle set (H/He/N/O/C) Photochem's own chemistry
  C  Photochem  NCHO_photo_network      the same chemistry as A, converted with
                                        photochem.utils.vulcan2yaml

A vs C isolates the code, B vs C isolates the reaction network.

Run A is produced by VULCAN itself (see docs/vulcan_photochem_comparison.md for
the configuration); this script produces B and C and then compares all three.

Usage
-----
  # needs an environment with photochem installed (conda install -c conda-forge photochem)
  python3 docs/compare_vulcan_photochem.py run     [--workdir DIR]
  python3 docs/compare_vulcan_photochem.py compare [--workdir DIR]

Everything that can be held identical between the codes is: the T(p)/Kzz
profile, the stellar spectrum, the surface gravity, the elemental abundances
and the zenith angle all come from the VULCAN configuration.
"""
import argparse
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VULCAN_DIR = os.path.join(HERE, '..', 'VULCAN')
VULCAN_OUT = os.path.join(VULCAN_DIR, 'output', 'HD189.vul')

R_SUN, AU, R_JUP, GRAV = 6.957e10, 1.495979e13, 7.1492e9, 6.674e-8

# --- HD 189733 b, taken from VULCAN's vulcan_cfg so gravity matches exactly ---
RP = 1.138*R_JUP                  # cm
GS = 2140.0                       # cm/s^2
MP = GS*RP**2/GRAV                # g
R_STAR, A_ORB = 0.805*R_SUN, 0.03142*AU
ZENITH_DEG = 48.0

# The elemental abundances VULCAN actually uses. With `use_solar = True` VULCAN
# ignores the O_H/C_H/N_H/He_H lines of vulcan_cfg.py and reads
# fastchem_vulcan/input/solar_element_abundances.dat (Lodders), so these are
# taken from that file, not from the config. Getting this wrong is easy and
# costly: an earlier pass of this comparison used the config values, which are a
# different solar set, and the He/H mismatch alone (0.0838 vs 0.0969) moved
# HeH_base by 16% -- the one quantity the EXHALE handoff is sensitive to.
ABUNDANCES = {'O': 6.06178e-4, 'C': 2.77588e-4, 'N': 8.18465e-5,
              'He': 9.69170e-2, 'H': 1.0}

SPECIES = ['H2', 'H', 'He', 'H2O', 'CO', 'CH4', 'CO2', 'NH3', 'N2', 'HCN', 'OH']
P_TARGETS = [1.0e3, 1.0e2, 1.0e1, 1.0e0]      # dyn/cm^2; 1 dyn/cm^2 = 1 ubar


# --------------------------------------------------------------------------
# shared inputs
# --------------------------------------------------------------------------
def write_flux_file(path):
    """VULCAN gives erg/cm^2/s/nm at the STELLAR SURFACE; Photochem wants
    mW/m^2/nm AT THE PLANET.  The two units are numerically equal, so only the
    (R_star/a)^2 dilution is applied."""
    wl, f_surf = np.loadtxt(
        os.path.join(VULCAN_DIR, 'atm/stellar_flux/sflux-HD189_Moses11.txt'), unpack=True)
    dilution = (R_STAR/A_ORB)**2
    with open(path, 'w') as fh:
        fh.write('Wavelength (nm)               Solar flux (mW/m^2/nm)\n')
        for w, f in zip(wl, f_surf*dilution):
            fh.write(f'{w:.6e}                  {f:.6e}\n')
    return dilution


def read_climate():
    """T(p) and Kzz, the same file VULCAN read.

    Photochem's own grid has to extend ABOVE the supplied climate grid, so the
    profile is cut below Photochem's default top of atmosphere (1e-7 bar =
    0.1 dyn/cm^2). The 1 ubar handoff EXHALE reads stays inside what is kept,
    but the upper boundary then sits closer to it than VULCAN's does -- see the
    caveats in docs/vulcan_photochem_comparison.md."""
    P, T, Kzz = np.loadtxt(os.path.join(VULCAN_DIR, 'atm/atm_HD189_Kzz.txt'),
                           skiprows=2, unpack=True)
    o = np.argsort(P)[::-1]                       # Photochem wants decreasing P
    P, T, Kzz = P[o], T[o], Kzz[o]
    keep = P >= 0.5
    return P[keep], T[keep], Kzz[keep]


# --------------------------------------------------------------------------
# the two Photochem runs
# --------------------------------------------------------------------------
def run_photochem(workdir, mechanism, thermo, data_dir, out_pkl):
    from photochem.extensions import gasgiants

    flux_file = os.path.join(workdir, 'hd189_flux_photochem.txt')
    write_flux_file(flux_file)
    P, T, Kzz = read_climate()

    pc = gasgiants.EvoAtmosphereGasGiant(mechanism, flux_file, MP, RP,
                                         solar_zenith_angle=ZENITH_DEG,
                                         thermo_file=thermo, data_dir=data_dir)
    pc.gdat.verbose = False

    tot = sum(ABUNDANCES.values())
    molfracs = np.array([ABUNDANCES[a]/tot for a in pc.gdat.gas.atoms_names])
    pc.gdat.gas.molfracs_atoms_sun = molfracs
    print('  atoms:', list(pc.gdat.gas.atoms_names))
    print('  species in mechanism:', len(pc.dat.species_names))

    pc.initialize_to_climate_equilibrium_PT(P, T, Kzz, 1.0, 1.0)
    pc.initialize_robust_stepper(pc.wrk.usol)

    for n in range(300):
        give_up = reached = False
        for _ in range(100):
            give_up, reached = pc.robust_step()
            if give_up or reached:
                break
        print(f'  block {n+1}: t = {pc.wrk.tn:.3e} s'
              f'{"  STEADY" if reached else ""}{"  GAVE UP" if give_up else ""}', flush=True)
        if give_up or reached:
            break

    sol = pc.return_atmosphere()
    with open(out_pkl, 'wb') as fh:
        pickle.dump({k: np.asarray(v) for k, v in sol.items()}, fh)
    print(f'  reached_steady_state = {reached}, gave_up = {give_up} -> {out_pkl}')


def do_run(workdir):
    from photochem.utils import zahnle_rx_and_thermo_files, vulcan2yaml

    os.makedirs(workdir, exist_ok=True)

    print('B  Photochem / Zahnle set')
    rx = os.path.join(workdir, 'pc_rxns.yaml')
    th = os.path.join(workdir, 'pc_thermo.yaml')
    zahnle_rx_and_thermo_files(atoms_names=['H', 'He', 'N', 'O', 'C'],
                               rxns_filename=rx, thermo_filename=th,
                               remove_reaction_particles=True)
    run_photochem(workdir, rx, th, None,
                  os.path.join(workdir, 'pc_hd189_solution.pkl'))

    print('C  Photochem / VULCAN NCHO network')
    conv = os.path.join(workdir, 'vul2yaml')
    os.makedirs(conv, exist_ok=True)
    net_txt = os.path.join(conv, 'NCHO_photo_network.txt')
    if not os.path.exists(net_txt):
        import shutil
        shutil.copy(os.path.join(VULCAN_DIR, 'thermo/NCHO_photo_network.txt'), net_txt)
    cwd = os.getcwd()
    try:
        os.chdir(conv)
        vulcan2yaml('NCHO_photo_network.txt', os.path.join(VULCAN_DIR, 'thermo'))
    finally:
        os.chdir(cwd)
    # vulcan2yaml emits a null placeholder reaction "He <=> He" (A = 0) because
    # VULCAN's network lists He as a non-reacting bulk species; Photochem
    # rejects it as a duplicate. Dropping it changes no rate.
    mech = os.path.join(conv, 'NCHO_photo_network.yaml')
    lines = open(mech).read().split('\n')
    out, i = [], 0
    while i < len(lines):
        if lines[i].strip() == '- equation: He <=> He':
            i += 2
            continue
        out.append(lines[i]);  i += 1
    open(mech, 'w').write('\n'.join(out))
    run_photochem(workdir, mech, mech, os.path.join(conv, 'vulcandata'),
                  os.path.join(workdir, 'pc_hd189_vulnet_solution.pkl'))


# --------------------------------------------------------------------------
# comparison
# --------------------------------------------------------------------------
def load_vulcan(path):
    with open(path, 'rb') as f:
        d = pickle.load(f)
    sp = list(d['variable']['species'])
    ymix = np.asarray(d['variable']['ymix'])
    return (np.asarray(d['atm']['pco']), np.asarray(d['atm']['Tco']),
            {s: ymix[:, sp.index(s)] for s in sp})


def load_photochem(path):
    with open(path, 'rb') as f:
        d = pickle.load(f)
    skip = ('pressure', 'temperature', 'alt', 'density')
    return (np.asarray(d['pressure']), np.asarray(d['temperature']),
            {k: np.asarray(v) for k, v in d.items() if k not in skip})


def at_P(P, arr, Pt):
    o = np.argsort(P)
    return float(np.interp(np.log(Pt), np.log(P[o]), np.asarray(arr)[o]))


def do_compare(workdir):
    runs = [('A VULCAN/NCHO', *load_vulcan(VULCAN_OUT)),
            ('B photochem/Zahnle', *load_photochem(
                os.path.join(workdir, 'pc_hd189_solution.pkl'))),
            ('C photochem/NCHO', *load_photochem(
                os.path.join(workdir, 'pc_hd189_vulnet_solution.pkl')))]

    for name, P, T, m in runs:
        print(f'{name:22s}: {len(P):4d} levels, P {P.max():.2e} -> {P.min():.2e} dyn/cm2, '
              f'{len(m):3d} species')

    for Pt in P_TARGETS:
        print(f'\n--- P = {Pt:.1e} dyn/cm^2 ---')
        print(f"  {'':20s} " + ' '.join(f'{s:>10}' for s in ['T [K]'] + SPECIES))
        for name, P, T, m in runs:
            vals = [f'{at_P(P, T, Pt):10.1f}']
            for s in SPECIES:
                vals.append(f'{at_P(P, m[s], Pt):10.2e}' if s in m else f'{"-":>10}')
            print(f'  {name:20s} ' + ' '.join(vals))

    print('\n=== what EXHALE reads (1 ubar) ===')
    for name, P, T, m in runs:
        qH2, qH, qHe = (at_P(P, m[s], 1.0) if s in m else float('nan')
                        for s in ('H2', 'H', 'He'))
        mu = 2.016*qH2 + 1.008*qH + 4.003*qHe
        print(f'  {name:20s} T = {at_P(P, T, 1.0):7.1f} K   q_H2 = {qH2:.4f}   '
              f'q_H = {qH:.4f}   q_He = {qHe:.4f}   mu(H/He) = {mu:.3f}')


# --------------------------------------------------------------------------
# base.inp handoff, written the same way for every source
# --------------------------------------------------------------------------
KB, MH, MJ, RJ = 1.380649e-16, 1.6726e-24, 1.898e30, 6.9911e9
MASS = {'H2': 2.016, 'H': 1.008, 'He': 4.003, 'H2O': 18.02, 'CH4': 16.04,
        'CO': 28.01, 'CO2': 44.01, 'N2': 28.01, 'NH3': 17.03, 'C2H2': 26.04,
        'HCN': 27.03, 'H2S': 34.08, 'S': 32.06, 'SO2': 64.06}


def write_base_inp(path, label, P, T, mix, mp_MJ, r1bar_RJ, pbase_bar=1e-6, kzz=1e9):
    """Hypsometric integration from 1 bar up to pbase with this run's own mu(p)
    and T(p) -- the same procedure src/utils/vulcan_to_base.py applies to a
    VULCAN output, used here for every source so the handoffs differ only
    through the composition."""
    import math
    mu = np.zeros_like(P, dtype=float)
    for sp, m in MASS.items():
        if sp in mix:
            mu = mu + m*np.asarray(mix[sp])
    p_bar = np.asarray(P)/1.0e6
    o = np.argsort(p_bar)[::-1]
    p_bar, Tp, mup = p_bar[o], np.asarray(T)[o], mu[o]

    def interp(y, pq):
        return float(np.interp(math.log(pq), np.log(p_bar[::-1]), y[::-1]))

    Mp, r = mp_MJ*MJ, r1bar_RJ*RJ
    lp1, lp0, n = math.log(1.0), math.log(pbase_bar), 2000
    dx, x = (lp0 - lp1)/n, math.log(1.0)
    for _ in range(n):
        def f(r_, lpq):
            pq = math.exp(lpq)
            return -KB*interp(Tp, pq)/(interp(mup, pq)*MH*(GRAV*Mp/r_**2))
        k1 = f(r, x);                k2 = f(r + 0.5*dx*k1, x + 0.5*dx)
        k3 = f(r + 0.5*dx*k2, x + 0.5*dx);  k4 = f(r + dx*k3, x + dx)
        r += dx*(k1 + 2*k2 + 2*k3 + k4)/6.0
        x += dx

    q2 = interp(np.asarray(mix['H2']), pbase_bar)
    qh = interp(np.asarray(mix['H']), pbase_bar)
    qhe = interp(np.asarray(mix['He']), pbase_bar)
    heh = qhe/max(qh + 2.0*q2, 1e-30)
    Tb, mub = interp(Tp, pbase_bar), interp(mup, pbase_bar)
    with open(path, 'w') as f:
        f.write('# base.inp -- lower-atmosphere handoff for EXHALE\n')
        f.write(f'# source: {label}\n')
        f.write(f'# q_H2={q2:.4e} q_H={qh:.4e} q_He={qhe:.4f} mu={mub:.3f}'
                '   (q_H and mu are comments; the keys below are read)\n')
        f.write(f'T_base    {Tb:.2f}\n')
        f.write(f'r_base    {r/RJ:.5f}\n')
        f.write(f'HeH_base  {heh:.6f}\n')
        f.write(f'Kzz_base  {kzz:.3e}\n')
        # Photochemical H2 partition and the level it refers to. EXHALE uses
        # q_H2_base in place of its chemical-equilibrium fit in the
        # molecular-base particle count, and evaluates that fit at p_base.
        f.write(f'q_H2_base {q2:.6e}\n')
        f.write(f'p_base    {pbase_bar:.3e}\n')
    print(f'  {label:22s} T_base={Tb:7.2f} K  r_base={r/RJ:.5f} R_J  '
          f'HeH_base={heh:.6f}  mu={mub:.3f}  (q_H={qh:.4f})')


def do_base(workdir, mp_MJ, r1bar_RJ):
    runs = [('VULCAN/NCHO', 'base_vulcan.inp', load_vulcan(VULCAN_OUT)),
            ('photochem/Zahnle', 'base_photochem_zahnle.inp',
             load_photochem(os.path.join(workdir, 'pc_hd189_solution.pkl'))),
            ('photochem/NCHO', 'base_photochem_ncho.inp',
             load_photochem(os.path.join(workdir, 'pc_hd189_vulnet_solution.pkl')))]
    print(f'HD 189733 b: M_p = {mp_MJ} M_J, r(1 bar) = {r1bar_RJ} R_J')
    for label, name, (P, T, mix) in runs:
        write_base_inp(os.path.join(workdir, name), label, P, T, mix, mp_MJ, r1bar_RJ)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('action', choices=['run', 'compare', 'base'])
    p.add_argument('--mp', type=float, default=1.237, help='planet mass [M_J]')
    p.add_argument('--r1bar', type=float, default=1.138, help='radius of the 1-bar level [R_J]')
    p.add_argument('--workdir', default=os.path.join(HERE, '..', 'vulcan_work', 'pc_compare'))
    a = p.parse_args()
    workdir = os.path.abspath(a.workdir)
    if a.action == 'run':
        do_run(workdir)
    elif a.action == 'base':
        do_base(workdir, a.mp, a.r1bar)
    else:
        if not os.path.exists(VULCAN_OUT):
            sys.exit(f'missing VULCAN output {VULCAN_OUT}; run VULCAN first '
                     '(see docs/vulcan_photochem_comparison.md)')
        do_compare(workdir)


if __name__ == '__main__':
    main()
