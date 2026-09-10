#!/usr/bin/env python3
"""Phase-P1 photochemistry comparison: matched reaction network AND matched
vertical domain, plus the H2 reaction budget at the EXHALE handoff level.

`docs/compare_vulcan_photochem.py` (2026-08-09) measured a 7.2x spread in the
neutral-hydrogen fraction at 1 microbar between VULCAN and Photochem and split
it 4.0x network / 1.8x code.  That split was contaminated: the Photochem runs
had their climate grid truncated at 0.5 dyn/cm^2, a factor 2 above the
comparison level, while VULCAN's domain top sat two decades above it.  This
script repeats the comparison with Photochem's top of atmosphere placed at
VULCAN's own P_t, so network / domain / code separate.

Arms (Photochem; VULCAN run A is reused from its own output file):

  ncho      the VULCAN NCHO_photo_network, converted with photochem.utils.vulcan2yaml
  zahnle    Photochem's own set restricted to H/He/N/O/C
  zahnle_s  the same set with sulfur, i.e. the gas-giant H/He/N/O/C/S mechanism
            of photochem.extensions.gasgiants

Usage
-----
  python3 docs/p1_matched_comparison.py run     --planet hd189 --arm ncho --toa 1.0e-2
  python3 docs/p1_matched_comparison.py compare --planet hd189
  python3 docs/p1_matched_comparison.py budget  --planet hd189 --arm ncho --toa 1.0e-2

Run this with the interpreter of the Photochem 0.8.4 environment built for the
2026-08-09 comparison; it is not importable from a stock Python.
Shared helpers come from docs/compare_vulcan_photochem.py so the 2026-08-09
procedure is not duplicated.
"""
import argparse
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_vulcan_photochem as cvp          # noqa: E402  (shared helpers)

ROOT = os.path.abspath(os.path.join(HERE, '..'))
# 7.1492e9 cm is the IAU 2015 nominal EQUATORIAL Jupiter radius (R_JUP_E,
# the one parameters.f90 and VULCAN's vulcan_cfg use) and 6.9911e9 cm is the
# VOLUMETRIC mean radius (R_JUP_V, the constant EXHALE carried before it was
# corrected).  Each planet below keeps the radius its run was made with, so
# the comparison still describes those runs; only the names are put right.
R_SUN, AU, R_JUP_E, R_JUP_V, GRAV = 6.957e10, 1.495979e13, 7.1492e9, 6.9911e9, 6.67430e-8

# Planet configurations.  Every entry repeats what the two codes were actually
# given in the runs being compared; the sources are docs/vulcan_photochem_comparison.md
# (HD 189733 b) and docs/base_composition_handoff_plan.md section 11.1 (HD 209458 b).
PLANETS = {
    'hd189': dict(
        vulcan_dir=os.path.join(ROOT, 'VULCAN'),
        vulcan_out=os.path.join(ROOT, 'VULCAN', 'output', 'HD189.vul'),
        atm='atm/atm_HD189_Kzz.txt',
        sflux='atm/stellar_flux/sflux-HD189_Moses11.txt',
        rp=1.138*R_JUP_E, gs=2140.0, r_star=0.805*R_SUN, a_orb=0.03142*AU,
        mp_MJ=1.237, r1bar_RJ=1.138, p_base_bar=1.0e-6),
    'hd209': dict(
        vulcan_dir=os.path.join(ROOT, 'vulcan_work', 'hd209_vulcan'),
        vulcan_out=os.path.join(ROOT, 'vulcan_work', 'hd209_vulcan', 'output', 'HD209.vul'),
        atm='atm/atm_HD209_Kzz.txt',
        sflux='atm/stellar_flux/Gueymard_solar.txt',
        rp=1.36*R_JUP_V, gs=1008.9, r_star=1.155*R_SUN, a_orb=0.0480*AU,
        mp_MJ=0.720, r1bar_RJ=1.36, p_base_bar=1.0e-4),
}

ZENITH_DEG = 48.0
# Lodders 2009 as VULCAN reads it from fastchem_vulcan/input/solar_element_abundances.dat
# (C 8.4434, He 10.9864, N 7.9130, O 8.7826, S 7.12 dex on H = 12).
ABUNDANCES = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5, O=6.06178e-4,
                  S=1.31826e-5)
ARMS = ('ncho', 'zahnle', 'zahnle_s')


def workdir(planet, toa):
    return os.path.join(ROOT, 'vulcan_work', 'pc_compare_p1',
                        f'{planet}_toa{toa:.0e}'.replace('+', '').replace('-0', '-'))


def sol_path(planet, arm, toa):
    return os.path.join(workdir(planet, toa), f'pc_{arm}_solution.pkl')


def budget_path(planet, arm, toa):
    return os.path.join(workdir(planet, toa), f'pc_{arm}_budget.pkl')


def write_flux_file(cfg, path):
    """VULCAN's stellar-flux files are at the stellar surface in
    erg/cm^2/s/nm; build_atm.py dilutes them by (R_star/a)^2.  Photochem wants
    mW/m^2/nm at the planet, numerically the same unit, so only the dilution
    is applied -- the same conversion compare_vulcan_photochem.py uses."""
    wl, f_surf = np.loadtxt(os.path.join(cfg['vulcan_dir'], cfg['sflux']), unpack=True)
    dilution = (cfg['r_star']/cfg['a_orb'])**2
    with open(path, 'w') as fh:
        fh.write('Wavelength (nm)               Solar flux (mW/m^2/nm)\n')
        for w, f in zip(wl, f_surf*dilution):
            fh.write(f'{w:.6e}                  {f:.6e}\n')
    return dilution


def read_climate(cfg, toa):
    """T(p), Kzz on the same file VULCAN read, cut so that Photochem's
    photochemical grid can sit above the climate grid: gasgiants requires
    3*TOA_pressure_avg < P_in[-1]."""
    P, T, Kzz = np.loadtxt(os.path.join(cfg['vulcan_dir'], cfg['atm']),
                           skiprows=2, unpack=True)
    o = np.argsort(P)[::-1]
    P, T, Kzz = P[o], T[o], Kzz[o]
    keep = P >= 3.05*toa
    return P[keep], T[keep], Kzz[keep]


def mechanism_files(arm, wdir, cfg):
    """Return (mechanism, thermo, data_dir) for one arm."""
    from photochem.utils import zahnle_rx_and_thermo_files, vulcan2yaml
    if arm in ('zahnle', 'zahnle_s'):
        atoms = ['H', 'He', 'N', 'O', 'C'] + (['S'] if arm == 'zahnle_s' else [])
        rx = os.path.join(wdir, f'{arm}_rxns.yaml')
        th = os.path.join(wdir, f'{arm}_thermo.yaml')
        if not os.path.exists(rx):
            zahnle_rx_and_thermo_files(atoms_names=atoms, rxns_filename=rx,
                                       thermo_filename=th,
                                       remove_reaction_particles=True)
        return rx, th, None
    conv = os.path.join(wdir, 'vul2yaml')
    os.makedirs(conv, exist_ok=True)
    mech = os.path.join(conv, 'NCHO_photo_network.yaml')
    if not os.path.exists(mech):
        import shutil
        shutil.copy(os.path.join(ROOT, 'VULCAN', 'thermo', 'NCHO_photo_network.txt'),
                    os.path.join(conv, 'NCHO_photo_network.txt'))
        cwd = os.getcwd()
        try:
            os.chdir(conv)
            vulcan2yaml('NCHO_photo_network.txt', os.path.join(ROOT, 'VULCAN', 'thermo'))
        finally:
            os.chdir(cwd)
        # vulcan2yaml emits a null "He <=> He" placeholder (A = 0) that Photochem
        # rejects as a duplicate; dropping it changes no rate.
        lines = open(mech).read().split('\n')
        out, i = [], 0
        while i < len(lines):
            if lines[i].strip() == '- equation: He <=> He':
                i += 2
                continue
            out.append(lines[i]); i += 1
        open(mech, 'w').write('\n'.join(out))
    return mech, mech, os.path.join(conv, 'vulcandata')


def run_arm(planet, arm, toa, blocks=400):
    from photochem.extensions import gasgiants
    cfg = PLANETS[planet]
    wdir = workdir(planet, toa)
    os.makedirs(wdir, exist_ok=True)
    flux_file = os.path.join(wdir, 'flux_photochem.txt')
    write_flux_file(cfg, flux_file)
    P, T, Kzz = read_climate(cfg, toa)
    mech, thermo, data_dir = mechanism_files(arm, wdir, cfg)
    mp = cfg['gs']*cfg['rp']**2/GRAV

    pc = gasgiants.EvoAtmosphereGasGiant(mech, flux_file, mp, cfg['rp'],
                                         solar_zenith_angle=ZENITH_DEG,
                                         thermo_file=thermo, data_dir=data_dir)
    pc.gdat.verbose = False
    pc.gdat.TOA_pressure_avg = toa
    tot = sum(ABUNDANCES[a] for a in pc.gdat.gas.atoms_names)
    pc.gdat.gas.molfracs_atoms_sun = np.array(
        [ABUNDANCES[a]/tot for a in pc.gdat.gas.atoms_names])
    print(f'  planet {planet}  arm {arm}  TOA {toa:.3e} dyn/cm2')
    print(f'  climate grid {P[0]:.3e} -> {P[-1]:.3e} dyn/cm2, {len(P)} levels, '
          f'T_top = {T[-1]:.0f} K')
    print('  atoms:', list(pc.gdat.gas.atoms_names),
          ' species:', len(pc.dat.species_names), flush=True)

    pc.initialize_to_climate_equilibrium_PT(P, T, Kzz, 1.0, 1.0)
    pc.initialize_robust_stepper(pc.wrk.usol)
    give_up = reached = False
    for n in range(blocks):
        for _ in range(100):
            give_up, reached = pc.robust_step()
            if give_up or reached:
                break
        print(f'  block {n+1}: t = {pc.wrk.tn:.3e} s'
              f'{"  STEADY" if reached else ""}{"  GAVE UP" if give_up else ""}',
              flush=True)
        if give_up or reached:
            break

    sol = pc.return_atmosphere()
    with open(sol_path(planet, arm, toa), 'wb') as fh:
        pickle.dump({'reached_steady_state': bool(reached), 'gave_up': bool(give_up),
                     **{k: np.asarray(v) for k, v in sol.items()}}, fh)
    save_budget(pc, planet, arm, toa)
    print(f'  reached_steady_state = {reached}, gave_up = {give_up}')


def save_budget(pc, planet, arm, toa, species=('H2', 'H2O', 'OH', 'H')):
    """Per-reaction production and loss on the model grid, for the species that
    decide whether the H2O/OH cycle drives H2 destruction at the handoff."""
    usol = pc.wrk.usol
    out = {'pressure': np.asarray(pc.wrk.pressure),
           'temperature': np.asarray(pc.var.temperature),
           'alt': np.asarray(pc.var.z)}
    for sp in species:
        if sp not in pc.dat.species_names:
            continue
        pl = pc.production_and_loss(sp, usol)
        out[sp] = {'production': np.asarray(pl.production),
                   'production_rx': list(pl.production_rx),
                   'loss': np.asarray(pl.loss),
                   'loss_rx': list(pl.loss_rx)}
    with open(budget_path(planet, arm, toa), 'wb') as fh:
        pickle.dump(out, fh)
    print(f'  budget -> {budget_path(planet, arm, toa)}')


def show_budget(planet, arm, toa, p_target, ntop=8):
    with open(budget_path(planet, arm, toa), 'rb') as fh:
        b = pickle.load(fh)
    P = b['pressure']
    k = int(np.argmin(np.abs(np.log(P) - np.log(p_target))))
    print(f'\n=== {planet} / {arm} / TOA {toa:.1e} : budget at P = {P[k]:.3e} dyn/cm^2 '
          f'(T = {b["temperature"][k]:.1f} K) ===')
    for sp in ('H2', 'H2O', 'OH'):
        if sp not in b:
            continue
        for kind in ('loss', 'production'):
            rates = b[sp][kind][k]
            names = b[sp][kind + '_rx']
            tot = float(np.sum(rates))
            print(f'  {sp} {kind}: total {tot:.3e} cm^-3 s^-1')
            for i in np.argsort(rates)[::-1][:ntop]:
                if rates[i] <= 0:
                    break
                print(f'    {100*rates[i]/tot:6.2f}%  {rates[i]:.3e}  {names[i]}')
    return b


def load_arm(planet, arm, toa):
    return cvp.load_photochem(sol_path(planet, arm, toa))


def do_compare(planet, toas, p_targets=(1.0e3, 1.0e2, 1.0e1, 1.0e0, 1.0e-1)):
    cfg = PLANETS[planet]
    runs = []
    if os.path.exists(cfg['vulcan_out']):
        runs.append(('A VULCAN/NCHO', *cvp.load_vulcan(cfg['vulcan_out'])))
    for toa in toas:
        for arm in ARMS:
            if os.path.exists(sol_path(planet, arm, toa)):
                runs.append((f'{arm}/TOA{toa:.0e}', *load_arm(planet, arm, toa)))
    for name, P, T, m in runs:
        print(f'{name:24s}: {len(P):4d} levels, P {P.max():.2e} -> {P.min():.2e} '
              f'dyn/cm2, {len(m):3d} species')
    for Pt in p_targets:
        print(f'\n--- P = {Pt:.1e} dyn/cm^2 ---')
        print(f"  {'':24s} " + ' '.join(f'{s:>10}' for s in ['T [K]'] + cvp.SPECIES))
        for name, P, T, m in runs:
            vals = [f'{cvp.at_P(P, T, Pt):10.1f}']
            for s in cvp.SPECIES:
                vals.append(f'{cvp.at_P(P, m[s], Pt):10.2e}' if s in m else f'{"-":>10}')
            print(f'  {name:24s} ' + ' '.join(vals))
    p_hand = cfg['p_base_bar']*1.0e6
    print(f'\n=== what EXHALE reads (p_base = {cfg["p_base_bar"]:.1e} bar) ===')
    for name, P, T, m in runs:
        qH2, qH, qHe = (cvp.at_P(P, m[s], p_hand) if s in m else float('nan')
                        for s in ('H2', 'H', 'He'))
        print(f'  {name:24s} T = {cvp.at_P(P, T, p_hand):7.1f} K   q_H2 = {qH2:.4f}   '
              f'q_H = {qH:.4f}   q_He = {qHe:.4f}')


def do_base(planet, toa, pbase_bar=None):
    cfg = PLANETS[planet]
    pbase_bar = pbase_bar or cfg['p_base_bar']
    wdir = workdir(planet, toa)
    tag = f'p{pbase_bar:.0e}'.replace('-0', '-')
    runs = [('VULCAN/NCHO', f'base_vulcan_{tag}.inp', cvp.load_vulcan(cfg['vulcan_out']))]
    for arm in ARMS:
        if os.path.exists(sol_path(planet, arm, toa)):
            runs.append((f'photochem/{arm}', f'base_{arm}_{tag}.inp',
                         load_arm(planet, arm, toa)))
    print(f'{planet}: M_p = {cfg["mp_MJ"]} M_J, r(1 bar) = {cfg["r1bar_RJ"]} R_J, '
          f'p_base = {pbase_bar:.1e} bar')
    for label, name, (P, T, mix) in runs:
        cvp.write_base_inp(os.path.join(wdir, name), label, P, T, mix,
                           cfg['mp_MJ'], cfg['r1bar_RJ'], pbase_bar=pbase_bar)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('action', choices=['run', 'compare', 'budget', 'base'])
    p.add_argument('--planet', default='hd189', choices=sorted(PLANETS))
    p.add_argument('--arm', default='ncho', choices=ARMS)
    p.add_argument('--toa', type=float, default=1.0e-2,
                   help='Photochem top-of-atmosphere pressure [dyn/cm^2]')
    p.add_argument('--toas', type=float, nargs='*', default=None,
                   help='TOA values to include in `compare`')
    p.add_argument('--pbase', type=float, default=None, help='handoff level [bar]')
    p.add_argument('--ptarget', type=float, default=None,
                   help='budget level [dyn/cm^2]; default is the handoff level')
    a = p.parse_args()
    if a.action == 'run':
        run_arm(a.planet, a.arm, a.toa)
    elif a.action == 'compare':
        do_compare(a.planet, a.toas if a.toas else [a.toa])
    elif a.action == 'budget':
        pt = a.ptarget if a.ptarget else PLANETS[a.planet]['p_base_bar']*1e6
        show_budget(a.planet, a.arm, a.toa, pt)
    else:
        do_base(a.planet, a.toa, a.pbase)


if __name__ == '__main__':
    main()
