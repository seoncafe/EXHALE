#!/usr/bin/env python3
"""T3: start the time-dependent RK march from each of the two states.

A JFNK solve returns a root of the discrete steady residual.  A root can be
a steady state the time-dependent equations hold, or one they run away from.
The march is the test: with the steady solver off, each state is handed to
the RK integrator and left to move.  A state the equations hold stays; a
spurious root leaves.

Both cells use the SAME handoff (the cool arm's profile and base.inp), which
T2 showed is not what distinguishes the states.  The reconstruction is PLM,
the scheme the two states were converged under, and no `Solver:` key is
given, so the JFNK finish never runs.  `du_th` is set below anything the run
can reach so the march is stopped by the step cap and not by a convergence
test -- except for the code's own "steady state reached (dtu < dtu_th)"
exit, which is itself a result and is reported when it happens.

Usage:  python3 T3_time_marching.py [n_steps]
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.abspath(os.path.join(HERE, '../../..'))
SCAN = os.path.join(EX, 'LHS1140b/exhale/kzz_profile_scan/kzz1e7')
COOL = os.path.join(SCAN, 'heh3p4/k02')
HOT = os.path.join(SCAN, 'heh3p399/k03')
BIN = os.path.join(EX, 'EXHALE.x')

TEMPLATE = """Planet name: LHS1140b_lower_profile
Log10 lower boundary number density [cm^-3]: 13.506
Planet radius [R_J]: 0.157692
Planet mass [M_J]: 0.0176220
Equilibrium temperature [K]: 226.0
Orbital distance [AU]: 0.0946
Escape radius [R_p]: 2.00
He/H number ratio: 3.4
2D approximate method: Mdot
Parent star mass [M_sun]: 0.1844
Spectrum type: Load from file..
Spectrum file: %(sed)s
Use only EUV? False
[E_low,E_mid,E_high] = [ 13.60 , 123.98 , 1.24e3 ]
Log10 of X-ray luminosity [erg/s]: 26.372
Log10 of EUV luminosity [erg/s]: 26.404
Grid type: Mixed
Numerical flux: HLLC
Reconstruction scheme: PLM
du_th [PLM,WENO3]: 1.0e-30
Include He23S? True
Load IC? True
Do only PP: False
Force start: False
He_diffusion: True
Domain mode: Spherical
Outer radius [R_p]: 30.0
Stellar radius [R_sun]: 0.2159
Base BC: pressure 1.0
Valve eps: 1.0e-4
Resid tol: 4.0e-4
Lower atmosphere profile: lower_atmosphere_profile.dat
"""

SED = os.path.join(EX, 'LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt')


def run_cell(tag, seed, nsteps, profile_dir):
    d = os.path.join(HERE, os.environ.get('T3_ROOT', 'T3_march'), tag)
    out = os.path.join(d, 'output')
    os.makedirs(out, exist_ok=True)
    for f in ('lower_atmosphere_profile.dat', 'base.inp'):
        shutil.copyfile(os.path.join(profile_dir, f), os.path.join(d, f))
    for name in ('Hydro_ioniz.txt', 'Ion_species.txt'):
        shutil.copyfile(os.path.join(seed, 'output', name),
                        os.path.join(out, name.replace('.txt', '_IC.txt')))
    with open(os.path.join(d, 'input.inp'), 'w') as fh:
        fh.write(TEMPLATE % {'sed': SED})

    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = '4'
    env['EXHALE_MAXSTEPS'] = str(nsteps)
    for k in list(env):
        if k.startswith('EXHALE_PTC'):
            del env[k]
    log = os.path.join(d, 'march.log')
    print('=== %s : marching %d steps from %s' % (tag, nsteps, seed))
    with open(log, 'w') as fh:
        rc = subprocess.call([BIN], cwd=d, env=env, stdout=fh,
                             stderr=subprocess.STDOUT)
    print('   rc = %d' % rc)

    # the post-processing pass, so *_adv.txt and Mdot exist for the endpoint
    for name in ('Hydro_ioniz.txt', 'Ion_species.txt'):
        src = os.path.join(out, name)
        if os.path.isfile(src):
            shutil.copyfile(src, os.path.join(
                out, name.replace('.txt', '_IC.txt')))
    inp = os.path.join(d, 'input.inp')
    txt = open(inp).read().replace('Do only PP: False', 'Do only PP: True')
    open(inp, 'w').write(txt)
    pplog = os.path.join(d, 'pp.log')
    with open(pplog, 'w') as fh:
        subprocess.call([BIN], cwd=d, env=env, stdout=fh,
                        stderr=subprocess.STDOUT)
    mdot = re.findall(r'steady-state Mdot\s*=\s*([-\d.EDed+]+)',
                      open(pplog, errors='replace').read())
    mdot = float(mdot[-1].replace('D', 'E')) if mdot else float('nan')
    stop = [ln.strip() for ln in open(log, errors='replace')
            if '-> converged' in ln or '-> stalled' in ln or 'NaN' in ln]
    print('   log10 Mdot = %.3f   stop: %s' % (mdot, stop[-1] if stop else '(step cap)'))
    return d, mdot, (stop[-1] if stop else '(step cap)')


def main():
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    which = sys.argv[2] if len(sys.argv) > 2 else 'both'
    cells = [(t, s) for t, s in (('from_cool', COOL), ('from_hot', HOT))
             if which in ('both', t)]
    res = []
    for tag, seed in cells:
        res.append((tag,) + run_cell(tag, seed, nsteps, COOL))
    print()
    print('%-12s %12s   %s' % ('cell', 'log10 Mdot', 'stop'))
    for tag, d, mdot, stop in res:
        print('%-12s %12.3f   %s' % (tag, mdot, stop))
    print()
    print('reference: cool 7.540, hot 7.870')


if __name__ == '__main__':
    main()
