#!/usr/bin/env python3
"""Read the LHS 1140 b model tree back and print where every case stands.

One row per case: what state it is in, the mass-loss rate, the He I 10830
red-pair equivalent width and the line shape, and the residual the solver
stopped at.  Under each group's rows, the composition at which that ladder
crosses the measured equivalent width, solved on the ladder itself.

    python3 models/status.py [--markdown] [--only <substring>] [--write]

`--write` puts the table into `LHS1140b/MODELS.md` under "## 7. Results" and
the provenance table under "## 8. How each model was reached", replacing
whatever those sections held.  The provenance table is read from the
`REPRODUCE.md` that each run writes in its own case directory.

Every number is read from the case directory, never from a memo:

  log10 Mdot   `Log10 of steady-state Mdot` in pp.log, else run.log
  ||R||        the solver's own " (JFNK|PTC) done info= ... ||R||=" line
  state        none / running / marching-stop / info=N, and whether the
               stationary certification accepted the state
  EW, depth,   the synthesized He I 10830 curve (`tpm_He10830.txt`) and the
  FWHM         three-Gaussian fit beside it (`tpm_He10830_metrics.txt`)
  flux         `Numerical flux:` of the case's own input.inp

A flux-closure rung is read from its last iteration directory `kNN/`, and its
converged reservoir from `closure_history.txt`.
"""

import argparse
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LHS = os.path.dirname(HERE)
EXHALE = os.path.dirname(LHS)
sys.path.insert(0, HERE)

from make_models import GROUPS, Group                       # noqa: E402

# The measurement the ladders bracket: Cherubim et al. (2026), red-pair
# equivalent width over the vacuum window below.
EW_OBS = 1.108

# The window and the air-to-vacuum factor of `LHS1140b/make_memo_figures.py`,
# so that a crossing read here and a crossing read there are the same number.
AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


# --------------------------------------------------------------------------
# reading one case
# --------------------------------------------------------------------------

def red_pair_equivalent_width(path):
    """[%A] over the measurement's vacuum window, from the
    instrument-convolved curve (column 3 of `tpm_He10830.txt`)."""
    if not os.path.isfile(path):
        return None
    s = np.loadtxt(path)
    if s.ndim != 2 or s.shape[1] < 3:
        return None
    lam = s[:, 0]*AIR
    exc = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    if not m.any():
        return None
    return float(np.trapz(exc[m], lam[m]))


def line_metrics(path):
    """`red_depth` [%] and `fwhm_A` [A] of the three-Gaussian fit."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith('#') or not line.strip():
                continue
            word = line.split()
            if len(word) >= 2:
                try:
                    out[word[0]] = float(word[1])
                except ValueError:
                    pass
    return out


DONE_RE = re.compile(r'^ \((?:JFNK|PTC)\) done info=(-?\d+)'
                     r'(?:\s+\|\|R\|\|=\s*([0-9.EDed+-]+))?')
STATIONARY_RE = re.compile(r'stationary solve returned info\s*=\s*(-?\d+)')
MDOT_RE = re.compile(r'steady-state Mdot\s*=\s*([-\d.EDed+]+)')


def solver_verdict(run_log):
    """(info, ||R||) of the steady solve, or (None, None).

    The partitioned stationary route states its outcome once, at the end,
    while each of its outer passes writes a `done info=` line of its own, so
    the last of those is the last hydrodynamic pass and not the solve; where
    the route's own line is present it is the verdict.  `||R||` stays the
    hydrodynamic norm of the last pass, which is what that line carries.
    """
    info = rnorm = None
    stationary = None
    if not os.path.isfile(run_log):
        return info, rnorm
    with open(run_log, errors='replace') as fh:
        for line in fh:
            m = DONE_RE.match(line)
            if m:
                info = int(m.group(1))
                if m.group(2):
                    rnorm = float(m.group(2).replace('D', 'E')
                                  .replace('d', 'e'))
                continue
            m = STATIONARY_RE.search(line)
            if m:
                stationary = int(m.group(1))
    return (stationary if stationary is not None else info), rnorm


def certification_verdict(run_log):
    """'certified', 'uncertified', or None when no verdict was printed."""
    if not os.path.isfile(run_log):
        return None
    verdict = None
    with open(run_log, errors='replace') as fh:
        for line in fh:
            if 'NOT CERTIFIED:' in line:
                verdict = 'uncertified'
            elif line.strip().startswith('CERTIFIED:'):
                verdict = 'certified'
    return verdict


def log10_mass_loss_rate(*logs):
    for path in logs:
        if not os.path.isfile(path):
            continue
        value = None
        with open(path, errors='replace') as fh:
            for line in fh:
                m = MDOT_RE.search(line)
                if m:
                    value = float(m.group(1).replace('D', 'E')
                                  .replace('d', 'e'))
        if value is not None:
            return value
    return None


def last_iterate(case_dir):
    """The highest-numbered completed iteration directory of a rung."""
    ks = sorted(d for d in os.listdir(case_dir)
                if re.fullmatch(r'k\d\d', d)
                and os.path.isdir(os.path.join(case_dir, d)))
    return os.path.join(case_dir, ks[-1]) if ks else None


def closure_reservoir(case_dir):
    """(k, He/H at the match, log10 Mdot) of the last recorded iteration."""
    path = os.path.join(case_dir, 'closure_history.txt')
    if not os.path.isfile(path):
        return None
    row = None
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith('#') or not line.strip():
                continue
            f = line.split()
            if len(f) >= 15:
                row = f
    if row is None:
        return None
    return int(row[0]), float(row[11]), float(row[12])


def numerical_flux(case_dir):
    """`Numerical flux:` the case is solved with, HLLC or ROE.

    Read from the file itself, not from `make_models.py`, so that a case
    whose input was edited by hand is reported as it stands.  A closure rung
    carries `input_template.inp` in place of `input.inp`.
    """
    for name in ('input.inp', 'input_template.inp'):
        path = os.path.join(case_dir, name)
        if not os.path.isfile(path):
            continue
        with open(path, errors='replace') as fh:
            for line in fh:
                if line.lstrip().startswith('#'):
                    continue
                key, _, value = line.partition(':')
                if key.strip().lower() == 'numerical flux':
                    return value.strip().upper()
    return None


class Case(object):
    """One case directory read back."""

    def __init__(self, group, heh):
        self.group = group
        self.heh = heh
        self.dir = os.path.join(HERE, group.name, 'HeH%s' % heh)
        self.exists = os.path.isdir(self.dir)
        self.closure = None
        self.info = self.rnorm = self.mdot = None
        self.ew = self.depth = self.fwhm = None
        self.cert = None
        self.state = 'none'
        self.flux = None
        if self.exists:
            self.flux = numerical_flux(self.dir)
            self._read()

    # ---- where the run products of this case live ------------------------
    @property
    def product_dir(self):
        if self.group.is_closure:
            return last_iterate(self.dir)
        return self.dir

    def _read(self):
        if self.group.is_closure:
            self.closure = closure_reservoir(self.dir)
        d = self.product_dir
        if d is None:
            self.state = 'none'
            return
        run_log = os.path.join(d, 'run.log')
        pp_log = os.path.join(d, 'pp.log')
        self.info, self.rnorm = solver_verdict(run_log)
        self.cert = certification_verdict(run_log)
        self.mdot = log10_mass_loss_rate(pp_log, run_log)
        self.ew = red_pair_equivalent_width(
            os.path.join(d, 'tpm_He10830.txt'))
        m = line_metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
        self.depth = m.get('red_depth')
        self.fwhm = m.get('fwhm_A')
        self.state = self._state(d, run_log)

    def _state(self, d, run_log):
        solved = os.path.join(d, 'output', 'Hydro_ioniz.txt')
        if not os.path.isfile(run_log):
            return 'none'
        if self.info is None:
            # The solver has not reported; either it is still marching or it
            # stopped on the flux criterion without a Newton finish.
            return 'running' if not os.path.isfile(solved) else 'marching-stop'
        state = 'info=%d' % self.info
        if self.cert == 'certified':
            state += ' certified'
        elif self.cert == 'uncertified':
            state += ' uncertified'
        return state


# --------------------------------------------------------------------------
# the crossing of one ladder
# --------------------------------------------------------------------------

def equivalent_width_crossing(heh, ew, target=EW_OBS):
    """He/H at which a ladder reaches `target`, on the log-log chord between
    the two rungs that bracket it.

    Returns (value, lo, hi) with the bracketing rungs, or None when the ladder
    does not bracket the target.  A chord and not a fit: the ladders are three
    or four rungs wide and a fit would put weight on rungs far from the
    crossing.
    """
    pairs = sorted((h, e) for h, e in zip(heh, ew)
                   if h is not None and e is not None and e > 0.0)
    for (h0, e0), (h1, e1) in zip(pairs, pairs[1:]):
        if (e0 - target)*(e1 - target) <= 0.0 and e0 != e1:
            t = ((np.log10(target) - np.log10(e0))
                 / (np.log10(e1) - np.log10(e0)))
            return (10.0**(np.log10(h0) + t*(np.log10(h1) - np.log10(h0))),
                    h0, h1)
    return None


# --------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------

# `flux` is `Numerical flux:` of the case: HLLC for the catalog, ROE for the
# seven lowest-XUV atomic cases, where the HLLC flux does not reach the
# stationary root on this grid (user decision of 2026-09-17;
# `docs/lhs1140b_stationary_L25_20260916.md`, step 3).  Where both fluxes
# solve, Roe against HLLC lowers Mdot and the He 10830 equivalent width by
# 0.3 to 0.6 per cent and leaves the first cell 6 to 7 per cent colder and 6
# to 8 per cent denser, with the mass-flux spread unchanged, so a ROE row is
# not on the same footing as an HLLC one to that accuracy.
COLUMNS = ('group', 'He/H', 'flux', 'state', 'reason', 'log10 Mdot',
           'red EW [%A]', 'red depth [%]', 'FWHM [A]', '||R||')

# A case that does not solve says so in its own directory, in `not_solved.md`:
# which cell and which row refuse it, and what holds it there.  The first line
# of that file is one sentence and is the `reason` column here; the file
# itself is linked from section 8.  A verdict without a reason is a verdict
# nobody can act on, which is why the column exists.
NOT_SOLVED = 'not_solved.md'


def not_solved_reason(case_dir):
    """The one sentence `not_solved.md` opens with, or None."""
    path = os.path.join(case_dir, NOT_SOLVED)
    if not os.path.isfile(path):
        return None
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            line = line.lstrip('#').strip()
            low = line.lower()
            if low.startswith('not solved:'):
                line = line[len('not solved:'):].strip()
            return line or None
    return None


def fmt(x, spec):
    return '--' if x is None else format(x, spec)


def rows_of_group(group, ladder):
    rows = []
    heh_num, ew_num = [], []
    for heh in ladder:
        c = Case(group, heh)
        label = heh
        if c.closure is not None:
            # A rung's composition is an output: show where it converged.
            label = '%s -> %.4f' % (heh, c.closure[1])
        reason = not_solved_reason(c.dir)
        rows.append([group.name, label, (c.flux or '--'), c.state,
                     (reason.replace('|', r'\|') if reason else '--'),
                     fmt(c.mdot, '.3f'),
                     fmt(c.ew, '.4f'), fmt(c.depth, '.3f'),
                     fmt(c.fwhm, '.4f'), fmt(c.rnorm, '.2E')])
        value = c.closure[1] if c.closure is not None else float(heh)
        heh_num.append(value)
        ew_num.append(c.ew)
    return rows, heh_num, ew_num


def build_table(only=None):
    lines = []
    for name, ladder in GROUPS:
        if only and only not in name:
            continue
        g = Group(name)
        rows, heh_num, ew_num = rows_of_group(g, ladder)
        lines.extend(rows)
        cross = equivalent_width_crossing(heh_num, ew_num)
        if cross is None:
            note = 'crossing He/H = no bracket (EW = %.3f %%A not spanned)' \
                % EW_OBS
        else:
            note = ('crossing He/H = %.4f  (log-log chord between the rungs'
                    ' %g and %g)' % cross)
        lines.append(['^' + name, note] + ['']*(len(COLUMNS) - 2))
    return lines


def render(lines, markdown):
    if markdown:
        # A markdown cell cannot carry a bare pipe, and one column is called
        # ||R||.
        head = [c.replace('|', r'\|') for c in COLUMNS]
        out = ['| ' + ' | '.join(head) + ' |',
               '|' + '|'.join(['---']*len(COLUMNS)) + '|']
        for row in lines:
            if row[0].startswith('^'):
                cells = ['**%s**' % row[0][1:], row[1]]
                cells += ['']*(len(COLUMNS) - len(cells))
            else:
                cells = list(row)
            out.append('| ' + ' | '.join(cells) + ' |')
        return '\n'.join(out)

    data = [r for r in lines if not r[0].startswith('^')] + [list(COLUMNS)]
    width = [max(len(str(r[i])) for r in data) for i in range(len(COLUMNS))]
    out = ['  '.join(c.ljust(w) for c, w in zip(COLUMNS, width))]
    out.append('  '.join('-'*w for w in width))
    for row in lines:
        if row[0].startswith('^'):
            out.append('    ' + row[1])
        else:
            out.append('  '.join(str(v).ljust(w)
                                 for v, w in zip(row, width)))
    return '\n'.join(out)


RESULTS_HEADING = '## 7. Results'
PROVENANCE_HEADING = '## 8. How each model was reached'


# --------------------------------------------------------------------------
# the provenance table, read from each case's own REPRODUCE.md
# --------------------------------------------------------------------------

# The seed line pick_seed.py writes, for a state of either origin: an
# archived directory under archive_20260830/, or a certified case of
# models/.  What is kept is the part of the path that names the state.
SEED_RE = re.compile(
    r'(?:archive_20260830|models)/(\S+)/output\s+(\S+)\s+HeH=(\S+)')
VERDICT_RE = re.compile(r'solver verdict: \*\*info = (-?\d+)\*\*')
CERT_RE = re.compile(r'certification of the state written: \*\*([A-Z ]+)\*\*')
REPRO_MDOT_RE = re.compile(r'log10 Mdot \[g/s\] = \*\*([-\d.]+)\*\*')
REPRO_EW_RE = re.compile(r'equivalent width = \*\*([\d.]+)\*\*')
WALL_RE = re.compile(r'wall clock \*\*(\S+)\*\* at `OMP_NUM_THREADS=(\d+)`')
PASS_ROW_RE = re.compile(r'(?m)^\|\s*(\d+)\s*\|\s*-?\d+\s*\|')


def reproduce_summary(case_dir):
    """The one-line provenance of a case, from the REPRODUCE.md its own run
    wrote.  `None` where no run has written one."""
    path = os.path.join(case_dir, 'REPRODUCE.md')
    if not os.path.isfile(path):
        return None
    text = open(path, errors='replace').read()
    seed = SEED_RE.search(text)
    passes = PASS_ROW_RE.findall(text)
    verdict = VERDICT_RE.search(text)
    cert = CERT_RE.search(text)
    mdot = REPRO_MDOT_RE.search(text)
    ew = REPRO_EW_RE.search(text)
    wall = WALL_RE.search(text)
    return dict(
        seed=('%s (%s, He/H %s)' % (seed.group(1), seed.group(2),
                                    seed.group(3))) if seed else 'cold start',
        passes=(max(int(k) for k in passes) if passes else 0),
        verdict=('info=%s' % verdict.group(1)) if verdict else '--',
        cert=(cert.group(1).strip().lower() if cert else '--'),
        mdot=mdot.group(1) if mdot else '--',
        ew=ew.group(1) if ew else '--',
        wall=('%s at %s thread(s)' % (wall.group(1), wall.group(2)))
        if wall else '--',
        path=os.path.relpath(path, LHS))


def provenance_table(only=None):
    rows = ['| case | seed | outer passes | verdict | certification |'
            ' log10 Mdot | EW [%A] | wall clock | record | why not solved |',
            '|---|---|---|---|---|---|---|---|---|---|']
    written = 0
    for name, ladder in GROUPS:
        if only and only not in name:
            continue
        g = Group(name)
        for heh in ladder:
            c = Case(g, heh)
            d = c.product_dir if g.is_closure else c.dir
            summary = reproduce_summary(c.dir)
            ns = os.path.join(c.dir, NOT_SOLVED)
            ns_cell = ('[%s](%s)' % (NOT_SOLVED, os.path.relpath(ns, LHS))
                       if os.path.isfile(ns) else '--')
            if summary is None:
                # A case that never ran writes no REPRODUCE.md; where it has
                # said WHY in `not_solved.md`, that is its row.
                if ns_cell == '--':
                    continue
                written += 1
                rows.append('| `%s/HeH%s` | -- | -- | -- | not solved | -- |'
                            ' -- | -- | -- | %s |' % (name, heh, ns_cell))
                continue
            written += 1
            rows.append('| `%s/HeH%s` | %s | %d | %s | %s | %s | %s | %s |'
                        ' [REPRODUCE.md](%s) | %s |'
                        % (name, heh, summary['seed'], summary['passes'],
                           summary['verdict'], summary['cert'],
                           summary['mdot'], summary['ew'], summary['wall'],
                           summary['path'], ns_cell))
    if not written:
        return ('No case has been run yet: `REPRODUCE.md` is written by'
                ' `run_case.sh` and `run_closure.sh` as each case finishes.')
    return '\n'.join(rows)


def replace_section(text, heading, body):
    """Replace one `## ` section of MODELS.md, keeping what follows it."""
    if heading not in text:
        raise SystemExit('MODELS.md has no "%s" section' % heading)
    head, _, tail = text.partition(heading)
    rest = re.split(r'(?m)^## ', tail, maxsplit=1)
    following = ('## ' + rest[1]) if len(rest) > 1 else ''
    return head + body + following


def write_into_models_md(table, provenance):
    """Replace the bodies of MODELS.md sections 7 and 8."""
    path = os.path.join(LHS, 'MODELS.md')
    with open(path) as fh:
        text = fh.read()
    if PROVENANCE_HEADING not in text:
        # Section 8 is written the first time a run is recorded, after 7.
        text = text.rstrip('\n') + '\n\n' + PROVENANCE_HEADING + '\n'
    text = replace_section(
        text, RESULTS_HEADING,
        '%s\n\nWritten by `models/status.py`; every number is read from the'
        ' case directory.\n\n%s\n\n' % (RESULTS_HEADING, table))
    text = replace_section(
        text, PROVENANCE_HEADING,
        '%s\n\nOne row per case that has been run, read from the'
        ' `REPRODUCE.md` its own run wrote: which state seeded it,'
        ' how many outer passes the stationary solve took, what the solver'
        ' and the certification said, what came out, and how long it took.'
        ' The linked file carries the commands themselves.\n\n%s\n'
        % (PROVENANCE_HEADING, provenance))
    with open(path, 'w') as fh:
        fh.write(text)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('--markdown', action='store_true',
                    help='render as a markdown table')
    ap.add_argument('--only', default=None,
                    help='only groups whose name contains this')
    ap.add_argument('--write', action='store_true',
                    help='write the table into MODELS.md section 7')
    a = ap.parse_args()

    lines = build_table(a.only)
    table = render(lines, a.markdown or a.write)
    print(table)
    if a.write:
        path = write_into_models_md(render(lines, True),
                                    provenance_table(a.only))
        print('\nwritten into %s (sections 7 and 8)'
              % os.path.relpath(path, LHS))
    return 0


if __name__ == '__main__':
    sys.exit(main())
