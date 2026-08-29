#!/usr/bin/env python3
"""Elemental-flux closure between the photochemical lower atmosphere and the
EXHALE escape wind (`docs/phase_e_flux_closure_design.md` section 6).

The physics.  The two models are joined at the microbar match, and the
handoff is one-way in each direction: the chemistry is told an elemental
flux through its upper boundary, and the wind, given the resulting profile
as its base state, carries an elemental flux of its own.  The solution is
the fixed point at which those agree,

    Phi_El = F_El(Phi_El),      El = H, He,

with

    F_He(r_f) = 4 pi r_f^2 ( rho_f v_f X_upwind + J_f )   [g/s]
    F_H (r_f) = 4 pi r_f^2 ( rho_f v_f (1 - X)   - J_f )  [g/s]

measured face by face by `write_element_flux_profile`
(`binary_element_diffusion.f90`) and reduced to one median and one relative
spread per element over a radial window (design section 3.4).  Nothing in a
Picard iteration of two independent solves keeps them from chasing each
other, so the update is under-relaxed,

    Phi^(k+1) = Phi^(k) + omega ( F^(k) - Phi^(k) ),

with `omega = 0.5` halved (floor 0.125) on any iteration whose residual
failed to fall.  The convergence test reads the *undamped* residual
`eps = |F - Phi| / |F|`, so a small `omega` cannot buy a false
convergence.  Both H and He are iterated here; design section 6.1 held He at
zero, and that is superseded.

Which window supplies the number is a physical question, not a detail.  The
overlap window -- between the base sound-wave region and the profile top --
is what the closure is defined on, but on LHS 1140 b it is only ~0.01 R_p
thick and the mass flux across it is not flat (measured spread 60%-700%,
against 2.9e7 g/s in the steady far field versus 4.7e9 g/s at 1.010 R_p).
The steady window (r >= r_esc) is flat to about 7% there.  This driver
therefore reads the overlap window first and falls back on the steady one
when the overlap is unmeasurable or its spread exceeds the tolerance,
recording in the iteration log which window supplied the number and why.
The substitution is never silent.  Section 6.2's honesty rule stands: a
residual smaller than the spread of the window it was measured on is not a
converged closure, and is reported UNRESOLVED rather than converged.

Usage
-----
    python3 element_flux_closure.py <case_dir> --phi0-H <g/s> [--phi0-He <g/s>]
        [--omega 0.5] [--tol 0.05] [--kmax 8] [--seed <dir>] [--config <file>]
        [--resume] [--dry-run]

`--dry-run <case_dir>` exercises only the reading and measurement path
against an existing EXHALE run directory; nothing is executed.
`--print-config-template` writes a commented configuration to stdout.

The configuration holds the fixed part of both command lines -- the planet
constants, the SED, the mechanism, the executable -- so that none of it is
spelled out inside this file.  Only the iteration-dependent options
(`--trial-flux-H`, `--trial-flux-He`, `--iteration`) are supplied here.
Its `python` key is optional: left out, the chemistry runs in the
interpreter this driver is running under, which is expected to be the one
Photochem is installed into (`README_photochem.md`).  Set it to reproduce a
result made on some other build, and the narrative records which of the two
supplied the interpreter.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

# --------------------------------------------------------------------------
# The interpreter that runs the chemistry
# --------------------------------------------------------------------------

# Photochem is a compiled extension with a dependency set of its own, so the
# chemistry step is run as a separate process.  By default it is run by the
# interpreter running this driver, which is the one Photochem is installed
# into (README_photochem.md).  A configuration that names "python"
# explicitly still wins, which is what keeps every stored closure.json --
# each of which carries an absolute interpreter path -- reproducing on the
# interpreter it was run on.


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

CONFIG_TEMPLATE = """{
  "comment": "Fixed part of both command lines for one planet. Everything the closure does not vary between iterations lives here (design section 6.4).",

  "python_comment": "Optional. Omitted, the chemistry runs in the interpreter this driver runs under, which is expected to be the one Photochem is installed into (README_photochem.md). Name an interpreter here only to reproduce a result made on a different Photochem build, and say in the comment which one.",
  "path_comment": "Paths are used as given. <EXHALE> below stands for the absolute path of this repository; substitute it.",
  "adapter": "<EXHALE>/src/utils/photochem_to_lower_profile.py",
  "adapter_run_dir": ".",
  "adapter_args": [
    "--mp", "0.0176220",
    "--r-ref", "0.157692",
    "--p-match", "1.0e-6",
    "--climate",
    "--climate-p-deep", "20.0",
    "--boa-pressure-factor", "1.0",
    "--stellar-flux", "<EXHALE>/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt",
    "--flux-at-planet",
    "--atoms", "H,He,N,O,C",
    "--abundances", "He=2.09",
    "--kzz-const", "1.0e9"
  ],
  "p_top_bar": null,

  "exhale_bin": "<EXHALE>/EXHALE.x",
  "input_template": "<EXHALE>/LHS1140b/exhale/input.inp",
  "omp_num_threads": 4,
  "resid_tol": "1.0e-4",
  "exhale_env": {"EXHALE_PTC": "1", "EXHALE_PTC_JFNK": "1",
                 "EXHALE_PTC_DTAU0": "1.0"}
}
"""

CONFIG_DEFAULTS = {
    'python': sys.executable,
    'adapter': os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'photochem_to_lower_profile.py'),
    'adapter_run_dir': '.',
    'adapter_args': [],
    'p_top_bar': None,
    'exhale_bin': None,
    'input_template': None,
    'omp_num_threads': 4,
    'exhale_env': {'EXHALE_PTC': '1', 'EXHALE_PTC_JFNK': '1',
                   'EXHALE_PTC_DTAU0': '1.0'},
    'profile_name': 'lower_atmosphere_profile.dat',
    'resid_tol': '1.0e-4',
    'python_named_by_config': False,
}


def read_configuration(path):
    """The fixed part of both command lines.  JSON; a YAML file is accepted
    when PyYAML happens to be importable, which is not a new dependency."""
    with open(path) as fh:
        text = fh.read()
    try:
        cfg = json.loads(text)
    except ValueError:
        try:
            import yaml
        except ImportError:
            raise SystemExit(
                'config %s is not JSON and PyYAML is not available' % path)
        cfg = yaml.safe_load(text)
    if not isinstance(cfg, dict):
        raise SystemExit('config %s must be a mapping' % path)
    merged = dict(CONFIG_DEFAULTS)
    merged.update({k: v for k, v in cfg.items()
                   if not k.endswith('comment')})
    # Whether the interpreter was named here or defaulted is worth saying in
    # the narrative: it is the one thing that told the 0.8.4 and the 0.9.0
    # chemistry apart when both environments existed.
    merged['python_named_by_config'] = 'python' in cfg
    return merged


# --------------------------------------------------------------------------
# The resolved configuration of one EXHALE run
# --------------------------------------------------------------------------

# Every key the closure reads, and whether its absence is fatal.  The reader
# degrades gracefully: an old EXHALE_resolved.out carries only the E1 key
# set, and the driver must be able to say which keys were missing rather
# than crash on the first lookup.
OVERLAP_KEYS = ('lower_profile_flux_state',
                'lower_profile_flux_r_lo_Rp',
                'lower_profile_flux_r_lo_source',
                'lower_profile_flux_r_hi_Rp',
                'lower_profile_flux_nface',
                'lower_profile_F_H_median', 'lower_profile_F_H_spread',
                'lower_profile_F_He_median', 'lower_profile_F_He_spread',
                'lower_profile_Mdot_median', 'lower_profile_Mdot_spread')

STEADY_KEYS = ('steady_flux_window_r_lo_Rp', 'steady_flux_window_nface',
               'steady_F_H_median', 'steady_F_H_spread',
               'steady_F_He_median', 'steady_F_He_spread',
               'steady_Mdot_median', 'steady_Mdot_spread')

PROVENANCE_KEYS = ('lower_profile_solution_id', 'lower_profile_iteration',
                   'lower_profile_trial_flux_H', 'lower_profile_trial_flux_He',
                   'HeH_number_ratio')


def resolved_configuration_path(run_dir):
    """`write_resolved_config` (write_setup_report.f90) writes
    `EXHALE_resolved.out` at the run-directory root, not under `output/`,
    unlike the other run products.  The `output/` path is checked first only
    as a defensive fallback in case that ever changes; today it is never
    where the file is found."""
    for rel in ('output/EXHALE_resolved.out', 'EXHALE_resolved.out'):
        p = os.path.join(run_dir, rel)
        if os.path.isfile(p):
            return p
    return None


def read_resolved_configuration(run_dir):
    """`key<whitespace>value`, one per line, `#` comments.  Returns the raw
    string map; the value of a notes key may contain spaces and is kept
    whole."""
    path = resolved_configuration_path(run_dir)
    if path is None:
        raise SystemExit('no EXHALE_resolved.out under %s' % run_dir)
    keys = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split(None, 1)
            if not parts:
                continue
            keys[parts[0]] = parts[1].strip() if len(parts) > 1 else ''
    return keys, path


def real_key(keys, name):
    """A real-valued key, or None when the run did not write it."""
    if name not in keys:
        return None
    try:
        return float(keys[name].replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return None


def missing_keys(keys, wanted):
    return [k for k in wanted if k not in keys]


# --------------------------------------------------------------------------
# Reducing the measured face fluxes to one number per element
# --------------------------------------------------------------------------

class FluxMeasurement(object):
    """One element pair (H, He) reduced over one radial window, with the
    spread that says how well-defined it is."""

    def __init__(self, window, F_H, F_He, spread_H, spread_He,
                 Mdot, Mdot_spread, nface, r_lo, reason):
        self.window = window          # 'overlap' | 'steady'
        self.F_H = F_H
        self.F_He = F_He
        self.spread_H = spread_H
        self.spread_He = spread_He
        self.Mdot = Mdot
        self.Mdot_spread = Mdot_spread
        self.nface = nface
        self.r_lo = r_lo
        self.reason = reason          # why this window, in words

    @property
    def spread(self):
        finite = [s for s in (self.spread_H, self.spread_He)
                  if s is not None and np.isfinite(s)]
        return max(finite) if finite else float('nan')

    def describe(self):
        return ('%s window: F_H = %.6E g/s (spread %.3G), '
                'F_He = %.6E g/s (spread %.3G), Mdot = %.6E g/s '
                '(spread %.3G), %s faces from r = %s R_p'
                % (self.window, self.F_H, _or_nan(self.spread_H),
                   self.F_He, _or_nan(self.spread_He),
                   _or_nan(self.Mdot), _or_nan(self.Mdot_spread),
                   self.nface, self.r_lo))


def _or_nan(x):
    return float('nan') if x is None else x


def overlap_flux_window(keys):
    """The window the closure is defined on: between the base sound-wave
    region and the profile top (design section 3.4).  Returns None when the
    run did not measure it, together with the state it reported."""
    state = keys.get('lower_profile_flux_state', 'unmeasured')
    if state != 'measured':
        return None, state
    F_H = real_key(keys, 'lower_profile_F_H_median')
    F_He = real_key(keys, 'lower_profile_F_He_median')
    if F_H is None or F_He is None:
        return None, 'unmeasured'
    m = FluxMeasurement(
        'overlap', F_H, F_He,
        real_key(keys, 'lower_profile_F_H_spread'),
        real_key(keys, 'lower_profile_F_He_spread'),
        real_key(keys, 'lower_profile_Mdot_median'),
        real_key(keys, 'lower_profile_Mdot_spread'),
        keys.get('lower_profile_flux_nface', '?'),
        keys.get('lower_profile_flux_r_lo_Rp', '?'),
        'measured over the overlap window (r_lo source %s)'
        % keys.get('lower_profile_flux_r_lo_source', '?'))
    return m, state


def steady_flux_window(keys):
    """The far-field window, r >= r_esc, where the wind is steady and the
    mass flux is flat.  It is not where the two models match, so a number
    taken here is an extrapolation of the closure variable to the match and
    the log must say so."""
    F_H = real_key(keys, 'steady_F_H_median')
    F_He = real_key(keys, 'steady_F_He_median')
    if F_H is None or F_He is None:
        return None
    return FluxMeasurement(
        'steady', F_H, F_He,
        real_key(keys, 'steady_F_H_spread'),
        real_key(keys, 'steady_F_He_spread'),
        real_key(keys, 'steady_Mdot_median'),
        real_key(keys, 'steady_Mdot_spread'),
        keys.get('steady_flux_window_nface', '?'),
        keys.get('steady_flux_window_r_lo_Rp', '?'),
        'steady far-field window (r >= r_esc)')


def select_flux_window(keys, tol):
    """Overlap first; the steady window only when the overlap cannot supply
    the number, and never silently.  Returns (measurement, note, missing)."""
    missing = missing_keys(keys, OVERLAP_KEYS + STEADY_KEYS)
    overlap, state = overlap_flux_window(keys)
    if overlap is not None and overlap.spread <= tol:
        return overlap, overlap.reason, missing

    if overlap is None:
        why = ('the overlap window is not measurable on this run '
               '(lower_profile_flux_state = %s)' % state)
    else:
        why = ('the overlap window is measurable but its spread %.3G '
               'exceeds tol = %.3G, so the flux is not flat across it'
               % (overlap.spread, tol))

    steady = steady_flux_window(keys)
    if steady is None:
        return None, why + '; and the steady window was not written either', \
            missing
    steady.reason = why + '; the number comes from the ' + steady.reason
    return steady, steady.reason, missing


# --------------------------------------------------------------------------
# The residual and the damped Picard update
# --------------------------------------------------------------------------

def elemental_flux_residual(measured, trial):
    """The undamped residual of design section 6.1,

        eps = |F_measured - Phi_trial| / |F_measured|,

    normalized on the MEASURED flux and not on the larger of the two: the
    wind's own statement of what it removes is the physical quantity, and a
    trial flux far above it must show as a large residual rather than
    saturate at 1.  Only both fluxes vanishing is exact agreement; a zero
    measured flux against a nonzero trial is a residual with no scale, and is
    reported as such rather than as a number."""
    if measured is None:
        return float('nan')
    if measured == 0.0:
        return 0.0 if trial == 0.0 else float('inf')
    if not np.isfinite(measured):
        return float('nan')
    return abs(measured - trial)/abs(measured)


def damped_picard_flux_update(trial, measured, omega):
    """Phi^(k+1) = Phi^(k) + omega ( F^(k) - Phi^(k) )."""
    return trial + omega*(measured - trial)


def relaxation_after(omega, eps, eps_prev, floor=0.125):
    """omega is halved on any iteration whose residual failed to fall below
    the previous one (design section 6.1); it never rises again and never
    goes below the floor."""
    if eps_prev is not None and np.isfinite(eps) and np.isfinite(eps_prev) \
            and eps >= eps_prev:
        return max(floor, 0.5*omega)
    return omega


# --------------------------------------------------------------------------
# The two solves
# --------------------------------------------------------------------------

def tail_of(path, nlines=40):
    try:
        with open(path, errors='replace') as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        return '(could not read %s: %s)' % (path, exc)
    return '\n'.join(lines[-nlines:])


def solve_photochemical_lower_profile(cfg, iter_dir, k, phi_H, phi_He, log):
    """The chemistry: the adapter, run with cwd = the iteration directory,
    writes `lower_atmosphere_profile.dat` and `base.inp` there.  A nonzero
    exit is a refusal or a non-steady chemistry solution, and neither has an
    elemental flux to hand over (design section 4.2)."""
    cmd = [cfg['python'], cfg['adapter'], cfg['adapter_run_dir']]
    cmd += [str(a) for a in cfg['adapter_args']]
    cmd += ['--trial-flux-H', repr(float(phi_H)),
            '--trial-flux-He', repr(float(phi_He)),
            '--iteration', str(k)]
    if cfg.get('p_top_bar') is not None:
        cmd += ['--p-top-bar', repr(float(cfg['p_top_bar']))]

    log('  chemistry: ' + ' '.join(cmd))
    logpath = os.path.join(iter_dir, 'adapter.log')
    with open(logpath, 'w') as fh:
        rc = subprocess.call(cmd, cwd=iter_dir, stdout=fh,
                             stderr=subprocess.STDOUT)
    if rc != 0:
        log('  chemistry FAILED rc=%d; tail of %s:' % (rc, logpath))
        log(tail_of(logpath))
        raise ClosureStop('the chemistry step exited %d (refusal or a'
                          ' non-steady solution); see %s' % (rc, logpath))

    profile = os.path.join(iter_dir, cfg['profile_name'])
    if not os.path.isfile(profile):
        log(tail_of(logpath))
        raise ClosureStop('the chemistry step exited 0 but wrote no %s'
                          % profile)
    log('  chemistry ok: ' + profile)
    return profile


INPUT_DROP = ('du_th', 'Solver', 'IC mode')

# The solver configuration of `LHS1140b/exhale/finish_case.sh`: load the
# seed, PLM, the valve and residual tolerances, marching keys removed.
#
# `resid_tol` defaults to 1.0e-4 and not the 1.0e-3 that script uses, and
# the reason is measured rather than cautionary.  The closure compares an
# elemental flux against its own radial spread, so the wind has to be flat
# to better than the tolerance before any residual is meaningful.  On the
# LHS 1140 b profile run, the same solution converged at 1.0e-3 leaves the
# steady window's mass flux spread at 7.4 per cent (F_H 7.1, F_He 7.9),
# above the 5 per cent closure tolerance; continuing the same solution to
# 1.0e-4 brings that spread to 0.50 per cent (F_H 0.47, F_He 0.57) and moves
# the mass-loss rate by 6 per cent, log10 Mdot 7.47 -> 7.50.  The 7 per cent
# was under-convergence, not a floor of the discretization.
#
# It is a configuration key and not a constant because the residual a given
# composition can reach is a property of that solution, not of the closure.
# On the helium-rich end of the LHS 1140 b reservoir ladder the JFNK line
# search stalls on the base contact mode at ||R|| = 1.1e-4 and returns its
# best iterate (info = 2) instead of converging, while the steady window it
# is measured on stays flat to well under the closure tolerance.  Raising
# the target for those runs records the looser number; it does not loosen
# the acceptance test, which is the window spread and is checked separately.
INPUT_JFNK = (('Load IC?', 'True'),
              ('Reconstruction scheme:', 'PLM'),
              ('Valve eps:', '1.0e-4'),
              ('Do only PP:', 'False'))


def write_input_keys(path, updates, drop=(), append_missing=True):
    """Rewrite `input.inp` keys in place, dropping the marching keys."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    seen = set()
    out = []
    for line in lines:
        stripped = line.strip()
        if any(stripped.startswith(d + ' ') or stripped.startswith(d + ':')
               or stripped.startswith(d + ' [') for d in drop):
            continue
        replaced = False
        for key, val in updates:
            if stripped.startswith(key):
                out.append('%s %s' % (key, val))
                seen.add(key)
                replaced = True
                break
        if not replaced:
            out.append(line)
    if append_missing:
        for key, val in updates:
            if key not in seen:
                out.append('%s %s' % (key, val))
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')


def solve_escape_wind(cfg, iter_dir, seed_output, log):
    """The wind, driven exactly as `LHS1140b/exhale/finish_case.sh` drives
    it: the seed solution becomes the IC, JFNK under EXHALE_PTC, then a
    post-processing pass on the solved state.  `output/` must exist before
    EXHALE starts."""
    out_dir = os.path.join(iter_dir, 'output')
    os.makedirs(out_dir, exist_ok=True)

    for name in ('Hydro_ioniz.txt', 'Ion_species.txt'):
        src = os.path.join(seed_output, name)
        if not os.path.isfile(src):
            raise ClosureStop('seed %s has no %s' % (seed_output, name))
        shutil.copyfile(src, os.path.join(
            out_dir, name.replace('.txt', '_IC.txt')))

    inp = os.path.join(iter_dir, 'input.inp')
    shutil.copyfile(cfg['input_template'], inp)
    write_input_keys(inp, INPUT_JFNK + (
        ('Resid tol:', str(cfg['resid_tol'])),
        ('Lower atmosphere profile:', cfg['profile_name']),), drop=INPUT_DROP)

    env = dict(os.environ)
    env.update({k: str(v) for k, v in cfg['exhale_env'].items()})
    env['OMP_NUM_THREADS'] = str(cfg['omp_num_threads'])

    runlog = os.path.join(iter_dir, 'run.log')
    log('  wind: JFNK start %s' % time.strftime('%H:%M:%S'))
    with open(runlog, 'w') as fh:
        rc = subprocess.call([cfg['exhale_bin']], cwd=iter_dir, env=env,
                             stdout=fh, stderr=subprocess.STDOUT)
    info = re.findall(r'done info=(\d+)', tail_of(runlog, 10**7))
    info = int(info[-1]) if info else None
    log('  wind: JFNK done rc=%d info=%s' % (rc, info))
    if info != 0:
        log(tail_of(runlog))
        raise ClosureStop('EXHALE did not report `done info=0` (rc=%d, '
                          'info=%s); see %s' % (rc, info, runlog))

    # post-processing pass on the solved state
    for name in ('Hydro_ioniz.txt', 'Ion_species.txt'):
        shutil.copyfile(os.path.join(out_dir, name),
                        os.path.join(out_dir, name.replace('.txt', '_IC.txt')))
    write_input_keys(inp, (('Do only PP:', 'True'),))
    # The post-processing pass runs WITHOUT the PTC variables, exactly as
    # `finish_case.sh` does (there they are set inline on the JFNK command
    # alone).  With EXHALE_PTC=1 still set the binary re-enters the steady
    # solver and stops before writing the advection-corrected profiles, so
    # `*_adv.txt` and the mass-loss rate never appear.
    pp_env = {k: v for k, v in env.items()
              if not k.startswith('EXHALE_PTC')}
    pplog = os.path.join(iter_dir, 'pp.log')
    with open(pplog, 'w') as fh:
        rc = subprocess.call([cfg['exhale_bin']], cwd=iter_dir, env=pp_env,
                             stdout=fh, stderr=subprocess.STDOUT)
    log('  wind: post-processing pass rc=%d' % rc)
    if rc != 0:
        log(tail_of(pplog))
        raise ClosureStop('the post-processing pass exited %d; see %s'
                          % (rc, pplog))
    return out_dir, info, log10_mass_loss_rate(pplog)


def log10_mass_loss_rate(pplog):
    """`Log10 of steady-state Mdot = ... g/s` from the post-processing log."""
    try:
        with open(pplog, errors='replace') as fh:
            text = fh.read()
    except OSError:
        return float('nan')
    hits = re.findall(r'steady-state Mdot\s*=\s*([-\d.EDed+]+)', text)
    if not hits:
        return float('nan')
    try:
        return float(hits[-1].replace('D', 'E').replace('d', 'e'))
    except ValueError:
        return float('nan')


# --------------------------------------------------------------------------
# The iteration record
# --------------------------------------------------------------------------

HISTORY_COLUMNS = ('k trial_F_H trial_F_He meas_F_H meas_F_He eps_H eps_He'
                   ' omega window_used window_spread_H window_spread_He'
                   ' HeH_match log10_Mdot solution_id exhale_info')

HISTORY_HEADER = (
    '# EXHALE elemental-flux closure history'
    ' (docs/phase_e_flux_closure_design.md section 6)\n'
    '# fluxes [g/s], eps and spreads [-], HeH_match = He/H number ratio at'
    ' the match\n'
    '# window_used: overlap = the match window; steady = the far field,'
    ' substituted with the reason recorded in closure.log\n'
    '# columns: ' + HISTORY_COLUMNS + '\n')


def history_path(case_dir):
    return os.path.join(case_dir, 'closure_history.txt')


def append_history_row(case_dir, row):
    path = history_path(case_dir)
    new = not os.path.isfile(path)
    with open(path, 'a') as fh:
        if new:
            fh.write(HISTORY_HEADER)
        fh.write('%4d %18.10E %18.10E %18.10E %18.10E %12.5E %12.5E'
                 ' %8.4f %-8s %12.5E %12.5E %14.7E %8.3f %s %s\n'
                 % (row['k'], row['trial_F_H'], row['trial_F_He'],
                    row['meas_F_H'], row['meas_F_He'],
                    row['eps_H'], row['eps_He'], row['omega'],
                    row['window_used'], _or_nan(row['window_spread_H']),
                    _or_nan(row['window_spread_He']), row['HeH_match'],
                    row['log10_Mdot'], row['solution_id'],
                    row['exhale_info']))


def read_history(case_dir):
    """Completed iterations, for `--resume`.  A row exists only for an
    iteration whose chemistry, wind and measurement all finished."""
    path = history_path(case_dir)
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith('#') or not line.strip():
                continue
            f = line.split()
            if len(f) < 15:
                continue
            rows.append({'k': int(f[0]), 'trial_F_H': float(f[1]),
                         'trial_F_He': float(f[2]), 'meas_F_H': float(f[3]),
                         'meas_F_He': float(f[4]), 'eps_H': float(f[5]),
                         'eps_He': float(f[6]), 'omega': float(f[7]),
                         'window_used': f[8],
                         'window_spread_H': float(f[9]),
                         'window_spread_He': float(f[10]),
                         'HeH_match': float(f[11]),
                         'log10_Mdot': float(f[12]),
                         'solution_id': f[13], 'exhale_info': f[14]})
    return rows


class ClosureStop(Exception):
    """A condition the closure records and stops on.  Nothing is retried
    silently: every one of these means the next iterate would be built on a
    solution that does not exist."""


def open_narrative(case_dir):
    fh = open(os.path.join(case_dir, 'closure.log'), 'a')

    def log(msg):
        stamp = time.strftime('%Y-%m-%d %H:%M:%S')
        for line in str(msg).splitlines() or ['']:
            fh.write('%s  %s\n' % (stamp, line))
        fh.flush()
        print(msg)
        sys.stdout.flush()
    return log, fh


# --------------------------------------------------------------------------
# The loop
# --------------------------------------------------------------------------

def iteration_directory(case_dir, k):
    return os.path.join(case_dir, 'k%02d' % k)


def measure_iteration(iter_dir, tol, log):
    """Read the resolved configuration of a finished wind solution and
    reduce it to one flux per element, saying which window supplied it."""
    keys, path = read_resolved_configuration(iter_dir)
    log('  measurement read from %s' % path)
    lacking = missing_keys(keys, PROVENANCE_KEYS)
    meas, note, lacking_flux = select_flux_window(keys, tol)
    if lacking or lacking_flux:
        log('  keys absent from this EXHALE_resolved.out (an older EXHALE'
            ' writes only the E1 set): ' + ', '.join(lacking + lacking_flux))
    if meas is None:
        raise ClosureStop(
            'neither flux window is usable: %s. The closure is UNRESOLVED on'
            ' this configuration.' % note)
    log('  ' + note)
    log('  ' + meas.describe())
    return keys, meas


def run_closure(case_dir, cfg, phi_H, phi_He, omega, tol, kmax, seed,
                resume, log):
    history = read_history(case_dir) if resume else []
    k0 = 0
    eps_prev = None
    eps_H = eps_He = float('nan')
    if history:
        last = history[-1]
        k0 = last['k'] + 1
        omega = relaxation_after(last['omega'],
                                 max(last['eps_H'], last['eps_He']),
                                 (max(history[-2]['eps_H'],
                                      history[-2]['eps_He'])
                                  if len(history) > 1 else None))
        eps_prev = max(last['eps_H'], last['eps_He'])
        phi_H = damped_picard_flux_update(last['trial_F_H'],
                                          last['meas_F_H'], omega)
        phi_He = damped_picard_flux_update(last['trial_F_He'],
                                           last['meas_F_He'], omega)
        log('resume: %d completed iterations, continuing at k=%d with'
            ' omega=%.4f, Phi_H=%.6E, Phi_He=%.6E'
            % (len(history), k0, omega, phi_H, phi_He))
        if eps_prev <= tol:
            log('resume: the last recorded residual %.3G already meets'
                ' tol=%.3G; nothing to do.' % (eps_prev, tol))
            return 0
        eps_H, eps_He = last['eps_H'], last['eps_He']

    # k_max = 8 means eight iterations, k = 0 .. k_max-1 (design section 6.2).
    for k in range(k0, kmax):
        iter_dir = iteration_directory(case_dir, k)
        if os.path.isdir(iter_dir):
            stale = iter_dir + '.incomplete.%s' % time.strftime('%Y%m%d%H%M%S')
            os.rename(iter_dir, stale)
            log('k=%d: an unrecorded iteration directory was set aside as %s'
                % (k, stale))
        os.makedirs(os.path.join(iter_dir, 'output'))
        log('k=%d: Phi_H = %.6E g/s, Phi_He = %.6E g/s, omega = %.4f'
            % (k, phi_H, phi_He, omega))

        solve_photochemical_lower_profile(cfg, iter_dir, k, phi_H, phi_He, log)

        seed_output = seed if k == k0 and seed else \
            os.path.join(iteration_directory(case_dir, k - 1), 'output')
        if not os.path.isdir(seed_output):
            raise ClosureStop('seed directory %s does not exist'
                              % seed_output)
        log('  seed: %s' % seed_output)
        _, info, log10_mdot = solve_escape_wind(cfg, iter_dir, seed_output,
                                                log)

        keys, meas = measure_iteration(iter_dir, tol, log)
        eps_H = elemental_flux_residual(meas.F_H, phi_H)
        eps_He = elemental_flux_residual(meas.F_He, phi_He)
        eps = max(eps_H, eps_He)
        log('  residuals (undamped): eps_H = %.4G, eps_He = %.4G' %
            (eps_H, eps_He))

        heh = real_key(keys, 'HeH_number_ratio')
        append_history_row(case_dir, {
            'k': k, 'trial_F_H': phi_H, 'trial_F_He': phi_He,
            'meas_F_H': meas.F_H, 'meas_F_He': meas.F_He,
            'eps_H': eps_H, 'eps_He': eps_He, 'omega': omega,
            'window_used': meas.window,
            'window_spread_H': meas.spread_H,
            'window_spread_He': meas.spread_He,
            'HeH_match': heh if heh is not None else float('nan'),
            'log10_Mdot': log10_mdot,
            'solution_id': keys.get('lower_profile_solution_id', 'unknown'),
            'exhale_info': str(info)})

        # Section 6.2's honesty rule.  The flux is only defined to within its
        # own radial spread, so a window whose spread exceeds the tolerance
        # cannot decide a residual at that tolerance whatever the residual
        # comes out to be: the closure is unresolved on this configuration
        # and says so, rather than converging on a number it cannot resolve.
        # This is a statement about the WINDOW, not about the iterate, so it
        # is tested every iteration and independently of how far eps is.
        if np.isfinite(meas.spread) and meas.spread > tol:
            raise ClosureStop(
                'UNRESOLVED on this configuration: the %s window that'
                ' supplied the flux has a radial spread %.3G above tol %.3G'
                ' (F_H spread %.3G, F_He spread %.3G), so a %.3G residual'
                ' cannot be decided on it. The measured fluxes are'
                ' F_H = %.6E, F_He = %.6E g/s.'
                % (meas.window, meas.spread, tol, _or_nan(meas.spread_H),
                   _or_nan(meas.spread_He), tol, meas.F_H, meas.F_He))

        if eps <= tol:
            log('CONVERGED at k=%d: max residual %.4G <= tol %.4G'
                % (k, eps, tol))
            log('  elemental fluxes at the match: F_H = %.6E g/s,'
                ' F_He = %.6E g/s; He/H = %s; log10 Mdot = %.3f'
                % (meas.F_H, meas.F_He, heh, log10_mdot))
            return 0

        omega = relaxation_after(omega, eps, eps_prev)
        eps_prev = eps
        phi_H = damped_picard_flux_update(phi_H, meas.F_H, omega)
        phi_He = damped_picard_flux_update(phi_He, meas.F_He, omega)

    raise ClosureStop(
        'k_max = %d iterations (k = 0 .. %d) exhausted without convergence;'
        ' last residuals eps_H = %.4G, eps_He = %.4G against tol = %.4G.'
        ' The history is %s.'
        % (kmax, kmax - 1, eps_H, eps_He, tol, history_path(case_dir)))


def report_measurement_only(run_dir, tol):
    """`--dry-run`: the reading and measurement path against an existing
    EXHALE run directory.  Nothing is executed."""
    keys, path = read_resolved_configuration(run_dir)
    print('read %s (%d keys)' % (path, len(keys)))
    lacking = missing_keys(keys, PROVENANCE_KEYS + OVERLAP_KEYS + STEADY_KEYS)
    if lacking:
        print('keys absent from this file (an older EXHALE writes only the'
              ' E1 set):')
        for k in lacking:
            print('   %s' % k)
    else:
        print('every key the closure reads is present.')
    for k in PROVENANCE_KEYS:
        if k in keys:
            print('  %-28s %s' % (k, keys[k]))
    meas, note, _ = select_flux_window(keys, tol)
    print('window selection: ' + note)
    if meas is None:
        print('RESULT: neither window is usable; the closure would report'
              ' UNRESOLVED on this run.')
        return 1
    print('RESULT: ' + meas.describe())
    trial_H = real_key(keys, 'lower_profile_trial_flux_H') or 0.0
    trial_He = real_key(keys, 'lower_profile_trial_flux_He') or 0.0
    print('residuals against the trial fluxes of this run:'
          ' eps_H = %.4G, eps_He = %.4G'
          % (elemental_flux_residual(meas.F_H, trial_H),
             elemental_flux_residual(meas.F_He, trial_He)))
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case_dir', nargs='?',
                    help='directory the iteration subdirectories and the'
                         ' records are written to; with --dry-run, an'
                         ' existing EXHALE run directory to read')
    ap.add_argument('--phi0-H', type=float, default=None,
                    help='initial trial elemental H flux [g/s]')
    ap.add_argument('--phi0-He', type=float, default=0.0,
                    help='initial trial elemental He flux [g/s] (default 0)')
    ap.add_argument('--omega', type=float, default=0.5,
                    help='initial under-relaxation (default 0.5; halved on'
                         ' any iteration whose residual fails to fall,'
                         ' floor 0.125)')
    ap.add_argument('--tol', type=float, default=0.05,
                    help='convergence tolerance on the undamped residual'
                         ' (default 0.05, the measured reproducibility'
                         ' floor of design section 6.2)')
    ap.add_argument('--kmax', type=int, default=8,
                    help='iteration limit (default 8)')
    ap.add_argument('--seed', default=None,
                    help='converged output/ directory seeding iteration 0')
    ap.add_argument('--config', default=None,
                    help='JSON with the fixed part of both command lines')
    ap.add_argument('--resume', action='store_true',
                    help='continue from the last completed iteration')
    ap.add_argument('--dry-run', action='store_true',
                    help='read and measure an existing run directory only')
    ap.add_argument('--print-config-template', action='store_true',
                    help='write a configuration template to stdout and exit')
    a = ap.parse_args()

    if a.print_config_template:
        sys.stdout.write(CONFIG_TEMPLATE)
        return 0
    if not a.case_dir:
        ap.error('case_dir is required')
    if a.dry_run:
        return report_measurement_only(a.case_dir, a.tol)

    if a.phi0_H is None:
        ap.error('--phi0-H is required')
    if not a.config:
        ap.error('--config is required: the planet constants and both fixed'
                 ' command lines live there, not in this file'
                 ' (--print-config-template writes one)')
    cfg = read_configuration(a.config)
    for key in ('exhale_bin', 'input_template'):
        if not cfg.get(key):
            raise SystemExit('config lacks %s' % key)
        if not os.path.isfile(cfg[key]):
            raise SystemExit('config %s = %s does not exist'
                             % (key, cfg[key]))
    if a.seed and not os.path.isdir(a.seed):
        raise SystemExit('--seed %s does not exist' % a.seed)
    if not a.resume and not a.seed:
        raise SystemExit('--seed is required for iteration 0')

    os.makedirs(a.case_dir, exist_ok=True)
    log, fh = open_narrative(a.case_dir)
    log('=' * 70)
    log('elemental-flux closure in %s' % os.path.abspath(a.case_dir))
    log('Phi_H(0) = %.6E g/s, Phi_He(0) = %.6E g/s, omega = %.3f,'
        ' tol = %.3G, kmax = %d'
        % (a.phi0_H, a.phi0_He, a.omega, a.tol, a.kmax))
    log('config %s; EXHALE %s' % (os.path.abspath(a.config),
                                  cfg['exhale_bin']))
    log('chemistry interpreter %s (%s)'
        % (cfg['python'],
           'named by the config' if cfg['python_named_by_config']
           else "this driver's own interpreter"))
    try:
        rc = run_closure(a.case_dir, cfg, a.phi0_H, a.phi0_He, a.omega,
                         a.tol, a.kmax, a.seed, a.resume, log)
    except ClosureStop as exc:
        log('STOPPED: %s' % exc)
        rc = 1
    finally:
        fh.close()
    return rc


if __name__ == '__main__':
    sys.exit(main())
