#!/usr/bin/env python3
"""Overlay EXHALE profiles on published escape models and tabulate Mdot.

Benchmarks (plan rev4 step E):
  salz2016   Salz et al. 2016, A&A 586, A75: HD 209458 b, HD 189733 b
             (CDS J/A+A/586/A75 tables, not digitized)
  huang2023  Huang et al. 2023, ApJ 951, 123: WASP-121 b cases A, B, D
             (Figs. 9, 12, 20, 25 digitized; Table 3)
  yan2022    Yan et al. 2022, ApJ 936, 177: WASP-52 b
             (Figs. 9, 12c, 23 digitized; Table 4)

Usage
  MPLBACKEND=Agg python3 benchmark_profile_comparison.py [salz2016|huang2023|yan2022|all]
        [--runs-root DIR] [--exhale-root DIR] [--out DIR]

Run directories: <runs-root>/salz2016_hd209458b, salz2016_hd189733b,
huang2023_wasp121b/case{A,B,D}, yan2022_wasp52b/<case>; the default
runs-root is the directory of this script (EXHALE_v1.00/benchmarks), the
default out is <runs-root>/comparison.

For every benchmark it writes, under <out>/<bench>/:
  <bench>_profiles.pdf            one page per run: published (markers, figure
                                  number in the legend) against EXHALE (solid:
                                  equilibrium file; dotted: *_adv file if present)
  <bench>_mdot_comparison.txt     run status, printed Mdot, 4 pi r^2 rho v at
                                  several radii in the paper's convention, ratio
  <bench>_profile_differences.txt EXHALE/published - 1 of T, v, n at a few radii

The script only reads files.  It works on a running or unconverged run: the
state is taken from output/Hydro_ioniz.txt, or, before the run has written
one, from the restart seed output/Hydro_ioniz_IC.txt (labeled as such), and
the run status is parsed from run.log (still running, final Mdot line or not,
exit status).  A file caught mid-write is read again after a pause.

Mdot conventions (where READ is stated in each benchmark's README.md):
  EXHALE  run.log prints log10(f x 4 pi r^2 rho v) at the physical cell
          N - 20 (EXHALE_main.f90, steady_mass_loss_rate), with f = 1/4 for
          '2D approximate method: Mdot/4', f = 1/2 for 'Rate/2 + Mdot/2',
          f = 1 otherwise ('Rate/4 + Mdot' dilutes the flux by 1/4 inside the
          run and makes no output correction).
  Salz    Mdot_sim = 1/4 x 4 pi R^2 rho v (Sect. 4.4).
  Huang   4 pi r^2 rho v along the substellar line (Sect. 5.3); Table 3 in
          M_p/Gyr with 1 M_p/Gyr = 7.1e13 g/s (Table 3 note d).
  Yan     4 pi r^2 rho v (full sphere; settled by Fig. 9 against Table 4,
          benchmark README).
"""

import argparse
import glob
import os
import re
import sys
import time
import importlib.util

import numpy as np
import matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif',
                     'font.size': 9, 'legend.fontsize': 6.5})

RUNS_ROOT = os.path.dirname(os.path.abspath(__file__))
EXHALE_ROOT = os.path.dirname(RUNS_ROOT)
OUT_ROOT = os.path.join(RUNS_ROOT, 'comparison')

GAMMA_ATOMIC = 5.0 / 3.0


def import_exhale_io(exhale_root):
    path = os.path.join(exhale_root, 'examples', 'exhale_io.py')
    spec = importlib.util.spec_from_file_location('exhale_io', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


aio = None   # set in main()


# ---------------------------------------------------------------------------
# EXHALE run: status from run.log, state from output/
# ---------------------------------------------------------------------------

_STEP_LINE = re.compile(r'^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$')
_MDOT_LINE = re.compile(r'Log10 of steady-state Mdot =\s*(\S+)')
_MACH_LINE = re.compile(r'Mach number of the outermost cell =\s*(\S+)\s+\(r =\s*(\S+)')
_EXIT_LINE = re.compile(r'exits (?:with status )?(\d)')
_METHOD_LINE = re.compile(r'2D approximate method:\s*(.+?)\s*$')
_CRASH = ('forrtl', 'Program received signal', 'ERROR STOP', 'error stop',
          'Segmentation fault')


def process_alive(run_dir):
    """True if an EXHALE binary (argv[0] EXHALE*.x) runs with run_dir as
    its working directory."""
    target = os.path.realpath(run_dir)
    for pid in os.listdir('/proc'):
        if not pid.isdigit():
            continue
        try:
            if os.path.realpath('/proc/%s/cwd' % pid) != target:
                continue
            with open('/proc/%s/cmdline' % pid, 'rb') as fh:
                cmd = fh.read().decode(errors='replace')
        except OSError:
            continue
        exe = os.path.basename(cmd.split('\0')[0])
        if re.match(r'EXHALE.*\.x$', exe):
            return True
    return False


def parse_run_log(run_dir):
    st = dict(log=False, alive=process_alive(run_dir), step=None, du=None,
              results=False, mdot_printed=None, mach_outer=None,
              r_outer=None, no_outflow=False, exit=None, crashed=False,
              method=None, last_stop='', killed=False)
    path = os.path.join(run_dir, 'run.log')
    if not os.path.isfile(path):
        st['status'] = 'RUNNING (no run.log)' if st['alive'] else 'NO run.log'
        st['final_mdot_line'] = False
        return st
    st['log'] = True
    with open(path, errors='replace') as fh:
        for line in fh:
            m = _STEP_LINE.match(line)
            if m:
                try:
                    st['step'] = int(m.group(1))
                    st['du'] = float(m.group(2))
                except ValueError:
                    pass
                continue
            if '----- Results -----' in line:
                st['results'] = True
            m = _MDOT_LINE.search(line)
            if m:
                st['mdot_printed'] = float(m.group(1))
            m = _MACH_LINE.search(line)
            if m:
                st['mach_outer'] = float(m.group(1))
                st['r_outer'] = float(m.group(2))
            if 'no outward mass flux at the evaluation radius' in line:
                st['no_outflow'] = True
            m = _EXIT_LINE.search(line)
            if m:
                st['exit'] = int(m.group(1))
            m = _METHOD_LINE.search(line)
            if m and '--->' in line:
                st['method'] = m.group(1)
            if '[stop] step' in line:
                st['last_stop'] = ' '.join(line.split())
            if any(c in line for c in _CRASH):
                st['crashed'] = True
            if 'SIGTERM' in line or 'SIGKILL' in line or 'process killed' in line:
                st['killed'] = True
    if st['alive']:
        st['status'] = 'RUNNING'
    elif st['killed']:
        st['status'] = 'KILLED (signal; no Results block)'
    elif st['crashed']:
        st['status'] = 'CRASHED'
    elif st['results']:
        st['status'] = 'FINISHED' + ('' if st['exit'] is None
                                     else ' exit %d' % st['exit'])
    else:
        st['status'] = 'STOPPED (no Results block)'
    st['final_mdot_line'] = st['mdot_printed'] is not None
    return st


def _coupling_header(path):
    with open(path) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            if line.startswith('# coupling') or line.startswith('# derived_from'):
                return ' '.join(line[1:].split()[1:])
    return ''


def _load_pair(hydro_path, ion_path):
    last = None
    for attempt in range(6):
        try:
            h = aio.load_hydro(hydro_path)
            r_ion, ion = aio.load_ions(ion_path)
            if h['r'].size != r_ion.size or not np.allclose(h['r'], r_ion):
                raise ValueError('hydro and ion radii differ (mid-write?)')
            if not np.all(np.isfinite(h['r'])):
                raise ValueError('non-finite radius')
            return h, ion
        except Exception as exc:     # a file caught mid-write
            last = exc
            time.sleep(3.0)
    raise RuntimeError('%s: could not be read: %s' % (hydro_path, last))


class ExhaleState:
    """One EXHALE run directory: profiles, Mdot and status."""

    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.status = parse_run_log(run_dir)
        out = os.path.join(run_dir, 'output')
        self.inp = aio.read_input(os.path.join(run_dir, 'input.inp'))
        res = os.path.join(run_dir, aio.RESOLVED_CONFIG_FILE)
        if os.path.isfile(res):
            self.resolved = aio.read_resolved_config(res)
            self.R0_RJ = float(self.resolved['planet_radius_RJ'])
            self.R0_source = 'EXHALE_resolved.out'
        else:
            self.resolved = {}
            self.R0_RJ = self.inp['Rp_RJ']
            self.R0_source = 'input.inp'
        self.method = self.inp['raw'].get('2D approximate method', 'Mdot')
        self.f_out = aio._mdot_factor(self.method)
        self.ok = False
        for suffix, label in (('', 'output/Hydro_ioniz.txt'),
                              ('_IC', 'output/Hydro_ioniz_IC.txt (restart '
                               'seed: no snapshot written yet)')):
            hp = os.path.join(out, 'Hydro_ioniz%s.txt' % suffix)
            ip = os.path.join(out, 'Ion_species%s.txt' % suffix)
            if os.path.isfile(hp) and os.path.isfile(ip):
                self.h, self.ion = _load_pair(hp, ip)
                self.source = label
                self.source_mtime = time.strftime(
                    '%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(hp)))
                self.coupling = _coupling_header(hp)
                self.ok = True
                break
        self.adv = None
        hp = os.path.join(out, 'Hydro_ioniz_adv.txt')
        ip = os.path.join(out, 'Ion_species_adv.txt')
        if self.ok and self.source.startswith('output/Hydro_ioniz.txt') \
                and os.path.isfile(hp) and os.path.isfile(ip):
            try:
                ha, iona = _load_pair(hp, ip)
                self.adv = (ha, iona)
            except RuntimeError:
                self.adv = None
        self.excited = None
        xp = os.path.join(out, 'Excited_H.txt')
        if self.ok and os.path.isfile(xp):
            try:
                self.excited = aio.load_excited_H(xp)
            except Exception:
                self.excited = None

    # ---- derived quantities, all on the eq file unless adv=True ----
    def _hi(self, adv):
        return (self.adv if adv else (self.h, self.ion))

    def r_R0(self, adv=False):
        return self._hi(adv)[0]['r']

    def mdot_full(self, adv=False):
        """4 pi r^2 rho v [g/s] at every physical cell (no 2D factor)."""
        h = self._hi(adv)[0]
        R0 = self.R0_RJ * aio.RJ
        return 4.0 * np.pi * h['n'] * aio.mu * h['v'] * (h['r'] * R0) ** 2

    def mdot_printed_equivalent(self):
        """log10 of what run.log prints, from the loaded state (cell N-20)."""
        md = self.mdot_full()
        j = md.size - 1 - 20
        if md[j] <= 0:
            return np.nan
        return np.log10(self.f_out * md[j])

    def sonic_radius(self):
        h = self.h
        cs = np.sqrt(GAMMA_ATOMIC * h['p'] / (h['n'] * aio.mu))
        idx = np.where(h['v'] >= cs)[0]
        return h['r'][idx[0]] if idx.size else np.nan

    def electron_density(self, adv=False):
        ion = self._hi(adv)[1]
        ne = np.zeros_like(self._hi(adv)[0]['r'])
        for name, n in ion.items():
            q = ion_charge(name)
            if q:
                ne = ne + q * n
        return ne


def ion_charge(name):
    if name in ('H2p', 'H3p', 'HeHp'):
        return 1
    if name in ('HeITR',):
        return 0
    m = re.match(r'^[A-Z][a-z]?(I+)$', name)
    if not m:
        return 0
    return len(m.group(1)) - 1


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def interp_logr(r_src, y_src, r_eval, logy=True):
    """Interpolate y(r) in log r (and log y when logy); nan outside range."""
    r_src = np.asarray(r_src, float)
    y_src = np.asarray(y_src, float)
    good = np.isfinite(r_src) & np.isfinite(y_src) & (r_src > 0)
    if logy:
        good &= y_src > 0
    if good.sum() < 2:
        return np.full(np.size(r_eval), np.nan)
    o = np.argsort(r_src[good])
    x = np.log10(r_src[good][o])
    y = y_src[good][o]
    y = np.log10(y) if logy else y
    xe = np.log10(np.atleast_1d(r_eval).astype(float))
    out = np.interp(xe, x, y, left=np.nan, right=np.nan)
    return 10.0 ** out if logy else out


def fmt(x, spec='%.3e'):
    if x is None:
        return '-'
    try:
        if not np.isfinite(x):
            return 'nan'
    except TypeError:
        return str(x)
    return spec % x


def status_label(ex):
    s = ex.status
    tag = s['status']
    if not s['final_mdot_line']:
        tag += ', UNCONVERGED/NO FINAL MDOT LINE'
    if 'certified=F' in getattr(ex, 'coupling', ''):
        tag += ', state certified=F'
    elif 'certified=T' in getattr(ex, 'coupling', ''):
        tag += ', state certified=T'
    return tag


def tex_escape(s):
    return (s.replace('\\', '/').replace('_', r'\_').replace('%', r'\%')
            .replace('&', r'\&').replace('#', r'\#'))


def exhale_lines(ax, x, y, color, label, adv_xy=None, logy=False):
    if logy:
        y = np.where(y > 0, y, np.nan)
    ax.plot(x, y, '-', color=color, lw=1.2, label=label)
    if adv_xy is not None:
        xa, ya = adv_xy
        if logy:
            ya = np.where(ya > 0, ya, np.nan)
        ax.plot(xa, ya, ':', color=color, lw=1.2)


def log_radius_axis(ax):
    """Log radius axis with plain-number labels at 1, 1.5, 2, 3, 5, 10, 15."""
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullFormatter
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(FixedLocator([1, 1.2, 1.5, 2, 3, 5, 7, 10, 15, 20]))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(NullFormatter())


def page_title(fig, title, ex, unit_note):
    s = ex.status
    sub1 = '%s; step %s, du %s' % (status_label(ex), s['step'],
                                     fmt(s['du'], '%.2e'))
    sub2 = 'profile: %s (%s)' % (ex.source, ex.source_mtime)
    fig.suptitle('\n'.join([tex_escape(title)] +
                           [r'{\footnotesize ' + tex_escape(t) + r'}'
                            for t in (sub1, sub2, unit_note)]), fontsize=9)


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def mdot_block(case, ex, pub_mdot, pub_note, conv_factor, conv_note,
               r_units, x_of_R0, extra_pub=None):
    """Text block for one run; conv_factor turns 4 pi r^2 rho v into the
    paper's convention; x_of_R0 converts r/R0 to the paper's radius unit."""
    s = ex.status
    L = []
    L.append('== %s  (%s)' % (case, ex.run_dir))
    L.append('  run status            : %s' % status_label(ex))
    L.append('  last marching step    : %s   du = %s' %
             (s['step'], fmt(s['du'], '%.3e')))
    if s['last_stop']:
        L.append('  last [stop] line      : %s' % s['last_stop'][:150])
    L.append('  profile read          : %s  (mtime %s)' % (ex.source,
                                                          ex.source_mtime))
    L.append('  state header          : %s' % (ex.coupling or '-'))
    L.append('  2D approximate method : %s  -> printed log Mdot = '
             'log10(%g x 4 pi r^2 rho v) at cell N-20' % (ex.method, ex.f_out))
    L.append('  printed log Mdot      : %s  (run.log, final line %s)' %
             (fmt(s['mdot_printed'], '%.2f'),
              'present' if s['final_mdot_line'] else 'ABSENT'))
    if s['no_outflow']:
        L.append('  run.log               : "no outward mass flux at the '
                 'evaluation radius" (no Mdot printed)')
    L.append('  same, from the profile: %s  (cell N-20 of the profile read)'
             % fmt(ex.mdot_printed_equivalent(), '%.3f'))
    mach_txt = (fmt(s['mach_outer'], '%.3f') + ' at r = %s R0 (run.log)'
                % fmt(s['r_outer'], '%.3f')) if s['mach_outer'] is not None \
        else '-'
    h = ex.h
    gam = GAMMA_ATOMIC
    mach_prof = abs(h['v'][-1]) / np.sqrt(gam * h['p'][-1] / (h['n'][-1] * aio.mu))
    L.append('  outer-cell Mach       : %s; from the profile (gamma 5/3) %.3f at '
             'r = %.3f R0' % (mach_txt, mach_prof, h['r'][-1]))
    L.append('  sonic point (gamma 5/3): %s R0 = %s %s' %
             (fmt(ex.sonic_radius(), '%.3f'),
              fmt(ex.sonic_radius() * x_of_R0, '%.3f'), r_units))
    L.append('  base radius R0        : %.4f R_J (%s); 1 R0 = %.5f %s' %
             (ex.R0_RJ, ex.R0_source, x_of_R0, r_units))
    L.append('  paper convention      : %s' % conv_note)
    L.append('  published Mdot        : %s g/s (log %s)  [%s]' %
             (fmt(pub_mdot), fmt(np.log10(pub_mdot), '%.3f'), pub_note))
    md = ex.mdot_full() * conv_factor
    x = ex.r_R0() * x_of_R0
    L.append('  %-14s %-14s %-12s %-10s %s' % ('r [%s]' % r_units, 'r [R0]',
                                               'EXHALE Mdot', 'log', 'EX/pub'))
    radii = list(extra_pub.get('radii', [])) if extra_pub else []
    for rr in radii:
        val = interp_logr(x, md, rr, logy=False)[0]
        L.append('  %-14.3f %-14.4f %-12s %-10s %s' %
                 (rr, rr / x_of_R0, fmt(val), fmt(np.log10(val), '%.3f')
                  if val > 0 else 'n/a', fmt(val / pub_mdot, '%.3f')))
    j = md.size - 1 - 20
    L.append('  %-14.3f %-14.4f %-12s %-10s %s   <- cell N-20 (the printed one)' %
             (x[j], ex.r_R0()[j], fmt(md[j]),
              fmt(np.log10(md[j]), '%.3f') if md[j] > 0 else 'n/a',
              fmt(md[j] / pub_mdot, '%.3f')))
    if extra_pub and 'pub_profile_mdot' in extra_pub:
        L.append('  published profile itself, same convention:')
        for rr, val in extra_pub['pub_profile_mdot']:
            L.append('    r = %.3f: %s g/s (log %s, /Table = %s)' %
                     (rr, fmt(val), fmt(np.log10(val), '%.3f'),
                      fmt(val / pub_mdot, '%.3f')))
    L.append('')
    return L


def diff_block(case, ex_label, rows, radii, r_units):
    """rows: list of (quantity, x_ex, y_ex, x_pub, y_pub, logy, source)."""
    L = ['== %s  (%s)' % (case, ex_label),
         '  relative difference EXHALE/published - 1 (interpolation in log r;'
         ' log y for v and n, linear for T); nan = outside a curve',
         '  %-26s %-12s' % ('quantity [figure]', 'r [%s] ->' % r_units)
         + ''.join('%10.2f' % r for r in radii)]
    for q, xe, ye, xp, yp, logy, src in rows:
        pe = interp_logr(xe, ye, radii, logy=logy)
        pp = interp_logr(xp, yp, radii, logy=logy)
        d = pe / pp - 1.0
        L.append('  %-39s' % ('%s [%s]' % (q, src))
                 + ''.join('%10s' % fmt(v, '%+.3f') for v in d))
    L.append('')
    return L


# ---------------------------------------------------------------------------
# Salz et al. 2016
# ---------------------------------------------------------------------------

SALZ_COLS = ['r', 'rho', 'v_kms', 'p', 'T', 'mu', 'heat', 'heat_frac',
             'fH', 'fHp', 'fHe', 'fHep', 'fHepp']


def read_salz(path):
    d = np.loadtxt(path, comments='#')
    out = {k: d[:, i] for i, k in enumerate(SALZ_COLS)}
    R = None
    with open(path) as fh:
        for line in fh:
            m = re.search(r'R_pl\s*=\s*([\d.]+)\s*R_jup', line)
            if m:
                R = float(m.group(1))
    out['R_pl_RJ'] = R
    # check the within-element normalization the README states
    out['sumH_dev'] = np.max(np.abs(out['fH'] + out['fHp'] - 1.0))
    out['sumHe_dev'] = np.max(np.abs(out['fHe'] + out['fHep'] + out['fHepp'] - 1.0))
    out['nonmono_r'] = int(np.sum(np.diff(out['r']) <= 0))
    return out


def bench_salz(runs_root, exhale_root, outdir):
    cases = [('HD 209458 b', 'salz2016_hd209458b', 10.27),
             ('HD 189733 b', 'salz2016_hd189733b', 9.61)]
    radii_mdot = [2.0, 3.0, 4.0, 8.0, 12.0]
    radii_diff = [1.05, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0]
    mlines = ['# Salz et al. 2016 (A&A 586, A75) against EXHALE.',
              '# Published log Mdot_sim: Table 3, 10.27 (HD 209458 b), 9.61 '
              '(HD 189733 b), READ (benchmark README).',
              '# Salz convention Mdot_sim = 1/4 x 4 pi R^2 rho v (Sect. 4.4, READ);'
              ' EXHALE values below are 1/4 x 4 pi r^2 rho v of the profile.',
              '# Radius: r/R_pl; EXHALE r/R0 x R0/R_pl (R0 from '
              'EXHALE_resolved.out, R_pl from the CDS header).', '']
    dlines = ['# Salz et al. 2016 profiles (CDS tables) against EXHALE.', '']
    pdf = PdfPages(os.path.join(outdir, 'salz2016_profiles.pdf'))
    for name, sub, logm in cases:
        rd = os.path.join(runs_root, sub)
        pub = read_salz(os.path.join(rd, 'published_salz2016_profiles.txt'))
        ex = ExhaleState(rd)
        if not ex.ok:
            mlines += ['== %s: no profile in %s' % (name, rd), '']
            continue
        k = ex.R0_RJ / pub['R_pl_RJ']
        xe = ex.r_R0() * k
        pub_md = 0.25 * 4 * np.pi * (pub['r'] * pub['R_pl_RJ'] * aio.RJ) ** 2 \
            * pub['rho'] * pub['v_kms'] * 1e5
        pp = [(rr, interp_logr(pub['r'], pub_md, rr, logy=False)[0])
              for rr in (2.0, 3.0, 4.0, 8.0)]
        mlines += mdot_block(name, ex, 10 ** logm, 'Table 3, READ', 0.25,
                             'Mdot_sim = 1/4 x 4 pi R^2 rho v (Salz Sect. 4.4)',
                             'R_pl', k, dict(radii=radii_mdot,
                                             pub_profile_mdot=pp))
        mlines.insert(len(mlines) - 1,
                      '  CDS table checks: max|H+H+ - 1| = %.1e, max|He+He+ +He++ - 1|'
                      ' = %.1e, nonincreasing radius steps = %d' %
                      (pub['sumH_dev'], pub['sumHe_dev'], pub['nonmono_r']))
        h, ion = ex.h, ex.ion
        nH = ion['HI'] + ion['HII']
        nHe = ion['HeI'] + ion['HeII'] + ion['HeIII']
        adv = ex.adv
        xa = adv[0]['r'] * k if adv else None

        def advxy(fun):
            return (xa, fun(*adv)) if adv else None

        fig, axs = plt.subplots(2, 3, figsize=(11, 7))
        ax = axs[0, 0]
        ax.plot(pub['r'], pub['rho'] / aio.mu, '-', lw=3, alpha=0.35, color='k',
                label='Salz 2016 Fig. 5a (CDS)')
        exhale_lines(ax, xe, h['n'], 'C0', 'EXHALE', advxy(lambda a, b: a['n']),
                     logy=True)
        ax.set_yscale('log')
        ax.set_ylabel(r'$\rho/m_{\mathrm H}\,[\mathrm{cm^{-3}}]$')
        ax = axs[0, 1]
        ax.plot(pub['r'], pub['T'], '-', lw=3, alpha=0.35, color='k',
                label='Salz 2016 Fig. 5b (CDS)')
        exhale_lines(ax, xe, h['T'], 'C3', 'EXHALE', advxy(lambda a, b: a['T']))
        ax.set_ylabel(r'$T\,[\mathrm{K}]$')
        ax = axs[0, 2]
        ax.plot(pub['r'], pub['v_kms'], '-', lw=3, alpha=0.35, color='k',
                label='Salz 2016 Fig. 5c (CDS)')
        exhale_lines(ax, xe, h['v'] / 1e5, 'C2', 'EXHALE',
                     advxy(lambda a, b: a['v'] / 1e5), logy=True)
        ax.set_yscale('log')
        ax.set_ylabel(r'$v\,[\mathrm{km\,s^{-1}}]$')
        ax = axs[1, 0]
        for key, en, c, lab in (('fH', 'HI', 'C0', r'H'), ('fHp', 'HII', 'C1', r'H$^+$')):
            ax.plot(pub['r'], pub[key], '-', lw=3, alpha=0.35, color=c,
                    label='Salz Fig. 6 ' + lab)
            exhale_lines(ax, xe, ion[en] / nH, c, 'EXHALE ' + lab,
                         advxy(lambda a, b, en=en: b[en] / (b['HI'] + b['HII'])),
                         logy=True)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 2)
        ax.set_ylabel(r'fraction of H nuclei')
        ax = axs[1, 1]
        for key, en, c, lab in (('fHe', 'HeI', 'C0', r'He'),
                                ('fHep', 'HeII', 'C1', r'He$^+$'),
                                ('fHepp', 'HeIII', 'C2', r'He$^{2+}$')):
            ax.plot(pub['r'], pub[key], '-', lw=3, alpha=0.35, color=c,
                    label='Salz Fig. 6 ' + lab)
            exhale_lines(ax, xe, ion[en] / nHe, c, 'EXHALE ' + lab,
                         advxy(lambda a, b, en=en: b[en] / (b['HeI'] + b['HeII'] + b['HeIII'])),
                         logy=True)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 2)
        ax.set_ylabel(r'fraction of He nuclei')
        ax = axs[1, 2]
        ax.plot(pub['r'], pub_md, '-', lw=3, alpha=0.35, color='k',
                label=r'Salz CDS $\pi r^2\rho v$')
        md = ex.mdot_full() * 0.25
        exhale_lines(ax, xe, md, 'C4', r'EXHALE $\pi r^2\rho v$', logy=True)
        ax.axhline(10 ** logm, color='k', ls='--', lw=0.8,
                   label=r'Salz Table 3: $10^{%.2f}$' % logm)
        ax.set_yscale('log')
        lo = 10 ** logm
        ax.set_ylim(lo / 30, lo * 30)
        ax.set_ylabel(r'$\frac{1}{4}\,4\pi r^2\rho v\,[\mathrm{g\,s^{-1}}]$')
        for a in axs.flat:
            log_radius_axis(a)
            a.set_xlabel(r'$r/R_{\mathrm{pl}}$')
            a.legend(loc='best')
        page_title(fig, 'Salz et al. 2016, %s' % name, ex,
                   'radius r/R_pl: EXHALE r/R0 x %.5f (R0 %.4f R_J resolved, '
                   'R_pl %.3f R_J CDS); solid EXHALE eq file, dotted _adv file; '
                   'markers published' % (k, ex.R0_RJ, pub['R_pl_RJ']))
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        pdf.savefig(fig)
        plt.close(fig)
        rows = [('rho/m_H', xe, h['n'], pub['r'], pub['rho'] / aio.mu, True, 'Fig. 5a'),
                ('T', xe, h['T'], pub['r'], pub['T'], False, 'Fig. 5b'),
                ('v', xe, h['v'] / 1e5, pub['r'], pub['v_kms'], True, 'Fig. 5c'),
                ('x(H+) in H', xe, ion['HII'] / nH, pub['r'], pub['fHp'], True, 'Fig. 6'),
                ('x(He+) in He', xe, ion['HeII'] / nHe, pub['r'], pub['fHep'], True, 'Fig. 6')]
        dlines += diff_block(name, '%s; %s' % (ex.run_dir, status_label(ex)),
                             rows, radii_diff, 'R_pl')
    pdf.close()
    write(outdir, 'salz2016_mdot_comparison.txt', mlines)
    write(outdir, 'salz2016_profile_differences.txt', dlines)


# ---------------------------------------------------------------------------
# Huang et al. 2023
# ---------------------------------------------------------------------------

HUANG_MDOT_GS_PER_MPGYR = 7.1e13   # Table 3 note d, READ


def read_huang_values(path):
    rows = {}
    markers = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('# Figure'):
                m = re.match(r'# Figure (\d+)', line)
                fig = m.group(1)
                d = {}
                for key, pat in (('sonic', r'sonic point ([\d.]+)'),
                                 ('L1', r'L1 ([\d.]+)'),
                                 ('x50', r'50% H ionization ([\d.]+)'),
                                 ('tau1', r'tau = 1(?: at 13.6 eV)? ([\d.]+)')):
                    mm = re.search(pat, line)
                    if mm:
                        d[key] = float(mm.group(1))
                markers[fig] = d
            if line.startswith('#') or not line.strip():
                continue
            w = line.split()
            rows[w[0]] = dict(Mp=float(w[2]), Mdot_MpGyr=float(w[12]))
    return rows, markers


def bench_huang(runs_root, exhale_root, outdir):
    pubdir = os.path.join(exhale_root, 'benchmarks', 'huang2023_wasp121b')
    rows, markers = read_huang_values(os.path.join(pubdir, 'published_values.txt'))
    cases = [('A', 'published_caseA_fig9_T_v.txt', '9', 'R_p', 1.766),
             ('B', 'published_caseB_fig20_T_v.txt', '20', 'R_px', 1.945),
             ('D', 'published_caseD_fig25_T_v.txt', '25', 'R_px', 1.959)]
    fig12 = np.loadtxt(os.path.join(pubdir, 'published_caseA_fig12_densities.txt'))
    radii_diff = [1.2, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0]
    mlines = ['# Huang et al. 2023 (ApJ 951, 123) against EXHALE, WASP-121 b.',
              '# Published Mdot: Table 3 in M_p/Gyr (READ, published_values.txt),'
              ' converted with 1 M_p/Gyr = 7.1e13 g/s (Table 3 note d, READ).',
              '# Huang convention: 4 pi r^2 rho v along the modeled line (Sect. 5.3);'
              ' EXHALE values below are 4 pi r^2 rho v of the profile (case A: 2 x'
              ' the printed value, which carries Mdot/2).',
              '# Radius: case A r/R_p (R_p = 1.766 R_J), B r/R_px (1.945 R_J), D'
              ' r/R_px (1.959 R_J) (Huang Tables 1, 2, READ); EXHALE r/R0 x R0/R_unit.',
              '']
    dlines = ['# Huang et al. 2023 digitized profiles against EXHALE.', '']
    pdf = PdfPages(os.path.join(outdir, 'huang2023_profiles.pdf'))
    for case, fname, fignum, runit, R_unit in cases:
        rd = os.path.join(runs_root, 'huang2023_wasp121b', 'case' + case)
        pub = np.loadtxt(os.path.join(pubdir, fname))
        rp, Tp, vp = pub[:, 0], pub[:, 3], pub[:, 4]
        ex = ExhaleState(rd)
        if not ex.ok:
            mlines += ['== case %s: no profile in %s' % (case, rd), '']
            continue
        k = ex.R0_RJ / R_unit
        xe = ex.r_R0() * k
        row = rows[case]
        pub_md = row['Mdot_MpGyr'] * HUANG_MDOT_GS_PER_MPGYR
        own = row['Mdot_MpGyr'] * row['Mp'] * aio.MJ / aio.GYR
        mk = markers.get(fignum, {})
        mlines += mdot_block(
            'case %s' % case, ex, pub_md,
            'Table 3: %.3g M_p/Gyr x 7.1e13 (note d); with this case\'s M_p = '
            '%.4f M_J: %.3e g/s' % (row['Mdot_MpGyr'], row['Mp'], own),
            1.0, '4 pi r^2 rho v (Huang Sect. 5.3)', runit, k,
            dict(radii=[2.0, 3.0, 4.0, 6.0]))
        mlines.insert(len(mlines) - 1,
                      '  published sonic point (Fig. %s marker): %s %s' %
                      (fignum, mk.get('sonic', 'n/a'), runit))
        h, ion = ex.h, ex.ion
        adv = ex.adv
        xa = adv[0]['r'] * k if adv else None
        fig, axs = plt.subplots(2, 2, figsize=(9, 7.6))
        axs = list(axs.flat)
        if case != 'A':
            axs[3].set_visible(False)
            axs = axs[:3]
        ax = axs[0]
        ax.plot(rp, Tp, 'o', ms=1.6, mfc='none', mew=0.6, color='k',
                label='Huang 2023 Fig. %s' % fignum)
        exhale_lines(ax, xe, h['T'], 'C3', 'EXHALE',
                     (xa, adv[0]['T']) if adv else None)
        ax.set_ylabel(r'$T\,[\mathrm{K}]$')
        ax = axs[1]
        ax.plot(rp, vp, 'o', ms=1.6, mfc='none', mew=0.6, color='k',
                label='Huang 2023 Fig. %s' % fignum)
        exhale_lines(ax, xe, h['v'], 'C2', 'EXHALE',
                     (xa, adv[0]['v']) if adv else None, logy=True)
        ax.set_yscale('log')
        ax.set_ylabel(r'$v\,[\mathrm{cm\,s^{-1}}]$')
        ax = axs[2]
        exhale_lines(ax, xe, ex.mdot_full(), 'C4', r'EXHALE $4\pi r^2\rho v$',
                     logy=True)
        ax.axhline(pub_md, color='k', ls='--', lw=0.8,
                   label=r'Huang Table 3: %.3g $M_p\,\mathrm{Gyr^{-1}}$' % row['Mdot_MpGyr'])
        ax.set_yscale('log')
        ax.set_ylim(pub_md / 30, pub_md * 30)
        ax.set_ylabel(r'$4\pi r^2\rho v\,[\mathrm{g\,s^{-1}}]$')
        for a in axs[:3]:
            for key, ls in (('sonic', '--'), ('L1', ':')):
                if key in mk:
                    a.axvline(mk[key], color='0.5', ls=ls, lw=0.7)
            rs = ex.sonic_radius() * k
            if np.isfinite(rs):
                a.axvline(rs, color='C0', ls='--', lw=0.7)
        axs[0].plot([], [], '--', color='0.5', lw=0.7, label='Huang sonic point')
        if 'L1' in mk:
            axs[0].plot([], [], ':', color='0.5', lw=0.7, label='Huang L1')
        axs[0].plot([], [], '--', color='C0', lw=0.7,
                    label=r'EXHALE sonic point ($\gamma=5/3$)')
        rows_d = [('T', xe, h['T'], rp, Tp, False, 'Fig. %s' % fignum),
                  ('v', xe, h['v'], rp, vp, True, 'Fig. %s' % fignum)]
        if case == 'A':
            ax = axs[3]
            ne = ex.electron_density()
            spec = [(2, 'e', None, 'C0'), (4, 'Mg I', 'MgI', 'C3'),
                    (5, 'Mg II', 'MgII', 'C1'), (6, 'Fe I', 'FeI', 'C2'),
                    (7, 'Fe II', 'FeII', 'C4'), (8, 'Ca II', 'CaII', 'C5')]
            for col, lab, en, c in spec:
                ax.plot(fig12[:, 0], fig12[:, col], 'o', ms=1.4, mfc='none', mew=0.6,
                        color=c, label='Fig. 12 ' + lab)
                y = ne if en is None else ion.get(en)
                if y is not None:
                    exhale_lines(ax, xe, y, c, None, logy=True)
                    rows_d.append(('n(%s)' % lab, xe, y, fig12[:, 0],
                                   fig12[:, col], True, 'Fig. 12'))
            if ex.excited is not None and 'n2p' in ex.excited:
                ax.plot(fig12[:, 0], fig12[:, 3], 'o', ms=1.4, mfc='none', mew=0.6,
                        color='k', label='Fig. 12 H(2p)')
                exhale_lines(ax, ex.excited['r'] * k, ex.excited['n2p'], 'k',
                             None, logy=True)
            ax.set_yscale('log')
            ax.set_ylim(1e-2, 1e12)
            ax.set_ylabel(r'$n\,[\mathrm{cm^{-3}}]$')
        for a in axs:
            log_radius_axis(a)
            a.set_xlabel(r'$r/R_{\mathrm{p}}$' if runit == 'R_p'
                         else r'$r/R_{\mathrm{px}}$')
            a.legend(loc='best')
        page_title(fig, 'Huang et al. 2023, WASP-121 b case %s' % case, ex,
                   'radius r/%s: EXHALE r/R0 x %.5f (R0 %.4f R_J resolved, %s %.3f R_J);'
                   ' solid EXHALE eq file, dotted _adv; markers published'
                   % (runit, k, ex.R0_RJ, runit, R_unit))
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        pdf.savefig(fig)
        plt.close(fig)
        dlines += diff_block('case %s' % case, '%s; %s' % (ex.run_dir,
                             status_label(ex)), rows_d, radii_diff, runit)
    pdf.close()
    write(outdir, 'huang2023_mdot_comparison.txt', mlines)
    write(outdir, 'huang2023_profile_differences.txt', dlines)


# ---------------------------------------------------------------------------
# Yan et al. 2022
# ---------------------------------------------------------------------------

def read_yan_long(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            q, f, r, v = line.split()
            out.setdefault((q, float(f)), []).append((float(r), float(v)))
    return {k: np.array(v) for k, v in out.items()}


def read_yan_table4(path, beta_col=1):
    """Mdot [g/s] at beta_m = 0.22 (column index 1 of the five)."""
    t = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            w = line.split()
            t[(w[0], float(w[1]))] = float(w[2 + beta_col]) * 1e11
    return t


def bench_yan(runs_root, exhale_root, outdir):
    pubdir = os.path.join(exhale_root, 'benchmarks', 'yan2022_wasp52b')
    fig9 = read_yan_long(os.path.join(pubdir, 'published_yan2022_fig9_hhe92_8.txt'))
    fig12 = read_yan_long(os.path.join(pubdir, 'published_yan2022_fig12c_hhe98_2.txt'))
    fig23 = read_yan_long(os.path.join(pubdir, 'published_yan2022_fig23_hhe98_2.txt'))
    tab4 = read_yan_table4(os.path.join(pubdir, 'published_yan2022_table4.txt'))
    Rp_RJ = 1.27   # Yan sec. 2.1, READ (benchmark README)
    runs = [('mol_h92he8_f1', '92/8', 1.0), ('h92he8_f0p5', '92/8', 0.5),
            ('h92he8_f1', '92/8', 1.0), ('h92he8_f2', '92/8', 2.0),
            ('h98he2_f0p5', '98/2', 0.5), ('h98he2_f1', '98/2', 1.0),
            ('h98he2_f2', '98/2', 2.0)]
    radii_diff = [1.05, 1.1, 1.5, 2.0, 3.0, 5.0, 8.0]
    mlines = ['# Yan et al. 2022 (ApJ 936, 177) against EXHALE, WASP-52 b.',
              '# Published Mdot: Table 4 at beta_m = 0.22, 1e11 g/s (READ,'
              ' published_yan2022_table4.txt).',
              '# Yan convention: 4 pi r^2 rho v, full sphere (benchmark README:'
              ' Fig. 9 against Table 4); EXHALE runs use Rate/4 + Mdot, so the'
              ' printed value is 4 pi r^2 rho v itself.',
              '# Radius: r/R_p with R_p = 1.27 R_J (Yan sec. 2.1, READ); EXHALE'
              ' r/R0 x R0/R_p.', '']
    dlines = ['# Yan et al. 2022 digitized profiles against EXHALE.', '']
    pdf = PdfPages(os.path.join(outdir, 'yan2022_profiles.pdf'))
    for sub, hhe, F in runs:
        rd = os.path.join(runs_root, 'yan2022_wasp52b', sub)
        if not os.path.isdir(rd):
            continue
        ex = ExhaleState(rd)
        name = '%s (H/He %s, %g F0)' % (sub, hhe, F)
        if not ex.ok:
            mlines += ['== %s: no profile yet in %s (status %s)' %
                       (name, rd, ex.status.get('status', '-')), '']
            continue
        k = ex.R0_RJ / Rp_RJ
        xe = ex.r_R0() * k
        pub_md = tab4[(hhe, F)]
        mlines += mdot_block(name, ex, pub_md, 'Table 4, beta_m 0.22, READ', 1.0,
                             '4 pi r^2 rho v (full sphere)', 'R_p', k,
                             dict(radii=[2.0, 4.0, 8.0]))
        if hhe == '92/8':
            src = {'T_K': ('Fig. 9a', fig9), 'v_kms': ('Fig. 9b', fig9),
                   'n_H1s_cm3': ('Fig. 9c', fig9), 'n_Hplus_cm3': ('Fig. 9c', fig9)}
        else:
            src = {'T_K': ('Fig. 23a', fig23), 'v_kms': ('Fig. 23b', fig23),
                   'n_He_total_cm3': ('Fig. 23c', fig23),
                   'n_Heplus_cm3': ('Fig. 23d', fig23),
                   'n_He23S_cm3': ('Fig. 23d', fig23)}
            if F == 1.0:
                src['n_H1s_cm3'] = ('Fig. 12c', fig12)
                src['n_Hplus_cm3'] = ('Fig. 12c', fig12)

        def P(q):
            if q not in src or (q, F) not in src[q][1]:
                return None
            a = src[q][1][(q, F)]
            o = np.argsort(a[:, 0])
            return src[q][0], a[o, 0], a[o, 1]

        h, ion = ex.h, ex.ion
        adv = ex.adv
        xa = adv[0]['r'] * k if adv else None
        ex_q = {'T_K': h['T'], 'v_kms': h['v'] / 1e5,
                'n_H1s_cm3': ion['HI'], 'n_Hplus_cm3': ion['HII'],
                'n_He_total_cm3': ion['HeI'] + ion['HeII'] + ion['HeIII'],
                'n_Heplus_cm3': ion['HeII'], 'n_He23S_cm3': ion['HeITR']}
        if adv:
            ha, ia = adv
            adv_q = {'T_K': ha['T'], 'v_kms': ha['v'] / 1e5,
                     'n_H1s_cm3': ia['HI'], 'n_Hplus_cm3': ia['HII'],
                     'n_He_total_cm3': ia['HeI'] + ia['HeII'] + ia['HeIII'],
                     'n_Heplus_cm3': ia['HeII'], 'n_He23S_cm3': ia['HeITR']}
        labels = {'T_K': r'$T$', 'v_kms': r'$v$', 'n_H1s_cm3': r'H(1s)',
                  'n_Hplus_cm3': r'H$^+$', 'n_He_total_cm3': r'He total',
                  'n_Heplus_cm3': r'He$^+$', 'n_He23S_cm3': r'He(2$^3$S)'}
        if 'H2' in ion:
            ex_q['n_H2_cm3'] = ion['H2']
            labels['n_H2_cm3'] = r'H$_2$'
        fig, axs = plt.subplots(2, 3, figsize=(12, 7.6))
        axs = list(axs.flat)
        axs[5].set_visible(False)
        axs = axs[:5]
        panels = [(axs[0], ['T_K'], False, r'$T\,[\mathrm{K}]$'),
                  (axs[1], ['v_kms'], False, r'$v\,[\mathrm{km\,s^{-1}}]$'),
                  (axs[2], ['n_H1s_cm3', 'n_Hplus_cm3'] +
                   (['n_H2_cm3'] if 'n_H2_cm3' in ex_q else []), True,
                   r'$n\,[\mathrm{cm^{-3}}]$'),
                  (axs[3], ['n_He_total_cm3', 'n_Heplus_cm3', 'n_He23S_cm3'],
                   True, r'$n\,[\mathrm{cm^{-3}}]$')]
        rows_d = []
        for ax, qs, logy, ylab in panels:
            for i, q in enumerate(qs):
                c = 'C%d' % i if len(qs) > 1 else ('C3' if q == 'T_K' else 'C2')
                p = P(q)
                if p is not None:
                    ax.plot(p[1], p[2], 'o', ms=1.4, mfc='none', mew=0.6, color=c,
                            label='Yan %s %s' % (p[0], labels[q]))
                    rows_d.append((q, xe, ex_q[q], p[1], p[2],
                                   q != 'T_K', p[0]))
                exhale_lines(ax, xe, ex_q[q], c, 'EXHALE ' + labels[q],
                             (xa, adv_q[q]) if adv and q in adv_q else None,
                             logy=logy)
            if logy:
                ax.set_yscale('log')
            ax.set_ylabel(ylab)
        ax = axs[4]
        exhale_lines(ax, xe, ex.mdot_full(), 'C4', r'EXHALE $4\pi r^2\rho v$',
                     logy=True)
        ax.axhline(pub_md, color='k', ls='--', lw=0.8,
                   label=r'Yan Table 4: $%.2f\times10^{11}$' % (pub_md / 1e11))
        ax.set_yscale('log')
        ax.set_ylim(pub_md / 30, pub_md * 30)
        ax.set_ylabel(r'$4\pi r^2\rho v\,[\mathrm{g\,s^{-1}}]$')
        for a in axs:
            log_radius_axis(a)
            a.set_xlabel(r'$r/R_{\mathrm{p}}$')
            a.legend(loc='best')
        page_title(fig, 'Yan et al. 2022, WASP-52 b, %s' % name, ex,
                   'radius r/R_p: EXHALE r/R0 x %.5f (R0 %.4f R_J resolved, R_p 1.27 R_J);'
                   ' solid EXHALE eq file, dotted _adv; markers published' %
                   (k, ex.R0_RJ))
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        pdf.savefig(fig)
        plt.close(fig)
        dlines += diff_block(name, '%s; %s' % (ex.run_dir, status_label(ex)),
                             rows_d, radii_diff, 'R_p')
    pdf.close()
    write(outdir, 'yan2022_mdot_comparison.txt', mlines)
    write(outdir, 'yan2022_profile_differences.txt', dlines)


# ---------------------------------------------------------------------------

def write(outdir, name, lines):
    stamp = '# written %s by benchmark_profile_comparison.py' % \
        time.strftime('%Y-%m-%d %H:%M %Z')
    with open(os.path.join(outdir, name), 'w') as fh:
        fh.write(stamp + '\n' + '\n'.join(lines) + '\n')


BENCHES = {'salz2016': bench_salz, 'huang2023': bench_huang, 'yan2022': bench_yan}


def main():
    global aio
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('bench', nargs='?', default='all',
                    choices=list(BENCHES) + ['all'])
    ap.add_argument('--runs-root', default=RUNS_ROOT,
                    help='directory holding the benchmark run directories')
    ap.add_argument('--exhale-root', default=EXHALE_ROOT,
                    help='EXHALE_v1.00 (exhale_io.py, published Huang/Yan data)')
    ap.add_argument('--out', default=OUT_ROOT)
    a = ap.parse_args()
    aio = import_exhale_io(a.exhale_root)
    for b in (BENCHES if a.bench == 'all' else [a.bench]):
        od = os.path.join(a.out, b)
        os.makedirs(od, exist_ok=True)
        BENCHES[b](a.runs_root, a.exhale_root, od)
        print('%s -> %s' % (b, od))


if __name__ == '__main__':
    main()
