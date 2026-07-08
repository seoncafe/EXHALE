#!/usr/bin/env python3
"""vulcan_driver.py -- run the bundled VULCAN photochemistry as an EXHALE
pre-step ("subroutine-style") and write the run's `base.inp`.

Invoked automatically by EXHALE when `input.inp` contains
    Lower atmosphere: vulcan <R_1bar[R_J]>
(see input_read.f90), or manually:
    python3 src/utils/vulcan_driver.py <run_dir> --r1bar <R_J> [--force]
                                       [--sflux FILE] [--atm FILE]

What it does
------------
1. Reads the planet parameters from <run_dir>/input.inp (name, Mp, Teq,
   orbital distance, stellar radius/Teff).
2. Prepares a private VULCAN work tree <run_dir>/vulcan_work/ (a copy of the
   bundled EXHALE/VULCAN, made once and reused), with:
   - a Guillot (2010) semi-grey T(p)+Kzz atmosphere file (unless --atm),
   - a stellar UV flux chosen by the host Teff from the shipped spectra
     (unless --sflux): >6800 K -> 51 Eri (F0, T7250); 5300-6800 -> solar
     (Gueymard); 4300-5300 -> eps Eri (K2V); <=4300 -> GJ 436 (M).
   - a vulcan_cfg.py for each planet (photochemistry ON, no live plotting).
3. Compiles FastChem once if its binary is missing.
4. Runs VULCAN to steady state (hours!) -- SKIPPED if a converged
   output/<name>-photo.vul already exists (use --force to redo).
5. Converts the result to <run_dir>/base.inp via the same machinery as
   vulcan_to_base.py (photochemical q_H2/q_H/q_He at 1 ubar, hypsometric
   base radius with VULCAN's own mu(p), T(p)).

Opt-out: simply omit the `Lower atmosphere:` key (EXHALE then uses its
classic base), or use `Lower atmosphere: analytic <R_1bar>` for the fast
chemical-equilibrium column (run_lower.py) instead of VULCAN.
"""
import argparse, math, os, shutil, subprocess, sys

HERE   = os.path.dirname(os.path.abspath(__file__))
EXROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
BUNDLE = os.path.join(EXROOT, 'VULCAN')

sys.path.insert(0, HERE)
from run_lower import guillot_T, RJ, MJ                     # noqa: E402

G_CGS = 6.674e-8


def read_input(run_dir):
    raw = {}
    for ln in open(os.path.join(run_dir, 'input.inp')):
        if ':' in ln:
            k, _, v = ln.partition(':')
            raw[k.strip()] = v.strip()
    def num(key, word=0):
        return float(raw[key].split()[word])
    out = {
        'name':  raw.get('Planet name', 'planet').split()[0],
        'Mp':    num('Planet mass [M_J]'),
        'Teq':   num('Equilibrium temperature [K]'),
        'a_AU':  num('Orbital distance [AU]'),
        'Rstar': num('Stellar radius [R_sun]') if 'Stellar radius [R_sun]' in raw else 1.0,
        'Teff':  num('Stellar Teff [K]') if 'Stellar Teff [K]' in raw else 5772.0,
    }
    return out


def pick_sflux(teff):
    if teff > 6800:
        return 'atm/stellar_flux/sflux-51-Eri-IUE-T7250.txt'
    if teff > 5300:
        return 'atm/stellar_flux/Gueymard_solar.txt'
    if teff > 4300:
        return 'atm/stellar_flux/sflux-epseri.txt'
    return 'atm/stellar_flux/sflux-GJ436.txt'


def write_guillot_atm(path, Teq, gs):
    import numpy as np
    P = np.logspace(9, -2, 140)                              # dyn/cm2
    T = [guillot_T(p / 1e6, Teq, Tint=100.0, g_cgs=gs) for p in P]
    with open(path, 'w') as f:
        f.write('#(dyne/cm2) (K)\t\t (cm2/s)\nPressure\tTemp  \t Kzz\n')
        for p, t in zip(P, T):
            f.write('%.3E\t%.1f\t %.3E\n' % (p, t, 1.0e9))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('run_dir')
    ap.add_argument('--r1bar', type=float, required=True,
                    help='radius of the 1-bar level [R_J] (transit radius)')
    ap.add_argument('--force', action='store_true', help='rerun VULCAN even if cached')
    ap.add_argument('--sflux', default=None, help='stellar flux file (bundle-relative)')
    ap.add_argument('--atm', default=None, help='T-P/Kzz file (skip Guillot)')
    a = ap.parse_args()
    run_dir = os.path.abspath(a.run_dir)

    P = read_input(run_dir)
    name = P['name']
    gs = G_CGS * P['Mp'] * MJ / (a.r1bar * RJ) ** 2
    print('(vulcan_driver) %s: Mp=%.3f MJ Teq=%.0f K a=%.4f AU Teff=%.0f -> gs=%.0f cm/s2'
          % (name, P['Mp'], P['Teq'], P['a_AU'], P['Teff'], gs))

    # 0. the VULCAN code is NOT shipped in the EXHALE repository; obtain it
    #    (once) with setup_vulcan.sh, which clones it and applies the EXHALE
    #    patches.  Auto-run the setup if the bundle is missing.
    if not os.path.isdir(BUNDLE):
        setup = os.path.join(HERE, 'setup_vulcan.sh')
        print('(vulcan_driver) VULCAN not found at %s -- running %s' % (BUNDLE, setup))
        rc = subprocess.call(['bash', setup])
        if rc != 0 or not os.path.isdir(BUNDLE):
            sys.exit('(vulcan_driver) ERROR: could not obtain VULCAN. Run\n'
                     '    %s\n'
                     'by hand (needs network access to clone exoclime/VULCAN).' % setup)

    # 1. private work tree (copy once)
    work = os.path.join(run_dir, 'vulcan_work')
    if not os.path.isdir(work):
        print('(vulcan_driver) copying VULCAN ->', work)
        shutil.copytree(BUNDLE, work)

    # 2. inputs
    atm_rel = 'atm/atm_%s_Kzz.txt' % name
    if a.atm:
        shutil.copy(a.atm, os.path.join(work, atm_rel))
    else:
        write_guillot_atm(os.path.join(work, atm_rel), P['Teq'], gs)
    sflux_rel = a.sflux or pick_sflux(P['Teff'])
    out_name = '%s-photo.vul' % name

    cfg = open(os.path.join(BUNDLE, 'vulcan_cfg.py')).read()
    subs = {
        "atm_file = 'atm/atm_HD189_Kzz.txt'":  "atm_file = '%s'" % atm_rel,
        "sflux_file = 'atm/stellar_flux/sflux-HD189_Moses11.txt'":
            "sflux_file = '%s'" % sflux_rel,
        "r_star = 0.805": "r_star = %.4f" % P['Rstar'],
        "Rp = 1.138*7.1492E9": "Rp = %.4f*7.1492E9" % a.r1bar,
        "orbit_radius = 0.03142": "orbit_radius = %.5f" % P['a_AU'],
        "gs = 2140.": "gs = %.1f" % gs,
        "out_name =  'HD189-photo.vul'": "out_name =  '%s'" % out_name,
    }
    for k, v in subs.items():
        if k not in cfg:
            sys.exit('(vulcan_driver) ERROR: cfg anchor missing: ' + k)
        cfg = cfg.replace(k, v, 1)
    open(os.path.join(work, 'vulcan_cfg.py'), 'w').write(cfg)

    # 3. FastChem binary
    fastchem = os.path.join(work, 'fastchem_vulcan', 'fastchem')
    if not os.path.exists(fastchem):
        print('(vulcan_driver) compiling FastChem ...')
        subprocess.check_call(['make'], cwd=os.path.join(work, 'fastchem_vulcan'),
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # 4. run VULCAN (cached unless --force)
    vul = os.path.join(work, 'output', out_name)
    if a.force or not os.path.exists(vul):
        print('(vulcan_driver) running VULCAN photochemistry for %s -- this '
              'takes HOURS; progress in %s/vulcan_run.log' % (name, work))
        with open(os.path.join(work, 'vulcan_run.log'), 'w') as log:
            rc = subprocess.call([sys.executable, 'vulcan.py'], cwd=work,
                                 stdout=log, stderr=subprocess.STDOUT)
        if rc != 0 or not os.path.exists(vul):
            sys.exit('(vulcan_driver) ERROR: VULCAN did not produce ' + vul)
    else:
        print('(vulcan_driver) using cached', vul)

    # 5. convert to base.inp
    rc = subprocess.call([sys.executable,
                          os.path.join(HERE, 'vulcan_to_base.py'),
                          vul, run_dir, '--mp', str(P['Mp']),
                          '--r1bar', str(a.r1bar)])
    if rc != 0:
        sys.exit('(vulcan_driver) ERROR: vulcan_to_base failed')
    print('(vulcan_driver) base.inp ready in', run_dir)


if __name__ == '__main__':
    main()
