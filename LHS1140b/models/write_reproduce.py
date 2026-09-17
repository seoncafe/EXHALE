#!/usr/bin/env python3
"""Write `REPRODUCE.md` beside a solved LHS 1140 b case.

The runners call this at the end of a case, so that the record of how a
result was reached is made by the run that reached it and not afterwards
from memory: what the case is, which binary solved it, which archived state
seeded it, the commands in the order they were given with their environment,
what the solver and the certification said pass by pass, and what came out.

    models/write_reproduce.py <case dir> --kind case|closure [options]

Everything the file states is read from the case directory (`input.inp`,
`run.log`, `pp.log`, `EXHALE_setup.out`, the `tpm_*` files) or from the
binary itself; the options carry only what no file records -- the seed line,
the commands, the clock and the thread count.
"""

import argparse
import hashlib
import os
import platform
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LHS = os.path.dirname(HERE)
EXHALE = os.path.dirname(LHS)

PASS_RE = re.compile(
    r'outer pass\s+(\d+):\s+hydro info=(-?\d+),\s+worst gated species row\s+'
    r'([0-9.EDed+-]+)\s+of\s+([0-9.EDed+-]+)\s+at cell\s+(\d+)\s+\(([^)]*)\),'
    r'\s+mass\s+([0-9.EDed+-]+),\s+momentum\s+([0-9.EDed+-]+),'
    r'\s+energy\s+([0-9.EDed+-]+).*?([0-9.]+)\s+s\s*$')
DONE_RE = re.compile(r'^ \((?:JFNK|PTC)\) done info=(-?\d+)'
                     r'(?:\s+\|\|R\|\|=\s*([0-9.EDed+-]+))?')
STATIONARY_RE = re.compile(r'stationary solve returned info\s*=\s*(-?\d+)')
MDOT_RE = re.compile(r'steady-state Mdot\s*=\s*([-\d.EDed+]+)')
MARCH_RE = re.compile(r'final: count=\s*(\d+)\s+du=\s*([0-9.EDed+-]+)')
# The two answers the post-processing pass reports, kept apart: whether the
# stationary claim the solved state was written with reproduces when that
# state is handed back, and what verdict the state the pass itself evaluated
# gets.  A passing work state is not a confirmation of the original claim.
CLAIM_RE = re.compile(r'^\s*original claim:\s*(.+?)\s*$', re.M)
WORK_RE = re.compile(r'^\s*work state verdict:\s*(.+?)\s*$', re.M)

# The measurement's own vacuum window and the air/vacuum factor, as
# `LHS1140b/make_memo_figures.py` and both runners use them.
AIR = 10832.057/10829.09114
EW_LO, EW_HI = 10832.60, 10834.20


def read(path):
    if not os.path.isfile(path):
        return ''
    with open(path, errors='replace') as fh:
        return fh.read()


def md5(path):
    if not os.path.isfile(path):
        return 'absent'
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def shell(cmd, cwd=EXHALE):
    try:
        out = subprocess.check_output(cmd, cwd=cwd,
                                      stderr=subprocess.DEVNULL)
        return out.decode('utf-8', 'replace').strip()
    except (OSError, subprocess.CalledProcessError):
        return ''


def compiler_of(binary):
    """The compiler that built the binary, from its own `.comment` section.

    `EXHALE_setup.out` states the configuration of the run and not the build,
    so this is where the compiler is actually recorded.
    """
    text = shell(['readelf', '-p', '.comment', binary], cwd=os.getcwd())
    names = []
    for line in text.splitlines():
        if ']' not in line:
            continue
        tag = line.split(']', 1)[1].strip()
        if any(word in tag for word in ('GCC', 'Intel', 'clang', 'ifort',
                                        'ifx')):
            names.append(tag)
    return '; '.join(names) if names else 'not recorded in the binary'


def group_fields(case):
    """`<chemistry>_<lower boundary>_<spectrum>_<mixing>` spelled out."""
    group = case.split('/')[0]
    part = group.split('_')
    if len(part) != 4:
        return []
    chem, boundary, spectrum, mixing = part
    chem_text = {'atomic': 'H, He, He(2^3S) and electrons; no molecular network',
                 'molecular': 'H2, H2+, H3+ and HeH+ with the H2 carrier'
                              ' transported'}.get(chem, chem)
    boundary_text = {
        'scalar': 'prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J,'
                  ' p = 1 microbar, metal-free',
        'scalarCNO': 'the same scalar base plus the C, N and O reservoirs of'
                     ' the photochemical column (base.inp)',
        'photochem': 'the Photochem column handed over as'
                     ' `Lower atmosphere profile:`'}.get(boundary, boundary)
    if mixing == 'wellmixed':
        mixing_text = 'no element diffusion; He/H uniform'
    elif mixing == 'kzzprofile':
        mixing_text = 'binary H/He element diffusion with K_zz(p) from the' \
                      ' profile'
    else:
        mixing_text = 'binary H/He element diffusion with He_Kzz = %s cm^2/s' \
                      % mixing[3:]
    return [('chemistry', '%s (`%s`)' % (chem_text, chem)),
            ('lower boundary', '%s (`%s`)' % (boundary_text, boundary)),
            ('spectrum', '`%s`' % spectrum),
            ('mixing', '%s (`%s`)' % (mixing_text, mixing))]


# The sentence a ROE case carries, with the flux-family systematic measured
# where both fluxes solve (`docs/lhs1140b_stationary_L25_20260916.md`, step 3,
# sections 3.3 to 3.5; user decision of 2026-09-17).
ROE_SENTENCE = (
    'The HLLC flux does not reach the stationary root on this grid at this'
    ' XUV level (L25 step 3), and the same case solved with HLLC on a grid'
    ' of twice the cells reaches the state the Roe flux reaches on the'
    ' catalog grid, so the difference is resolution and not a different'
    ' wind. The Roe flux systematic against HLLC where both solve is Mdot'
    ' and the He 10830 equivalent width lower by 0.3 to 0.6 per cent, the'
    ' first cell colder by 6 to 7 per cent and denser by 6 to 8 per cent,'
    ' and the mass-flux spread unchanged (L25 memo). That systematic is'
    ' carried wherever the numbers of this case are quoted.')


def numerical_flux(case_dir):
    """`Numerical flux:` of the case, from the file it is solved with.

    A closure rung has no `input.inp` of its own; its template is what each
    iteration is written from.
    """
    for name in ('input.inp', 'input_template.inp'):
        for line in read(os.path.join(case_dir, name)).splitlines():
            if line.lstrip().startswith('#'):
                continue
            key, _, value = line.partition(':')
            if key.strip().lower() == 'numerical flux':
                return value.strip().upper()
    return None


def outer_passes(run_log_text):
    rows = []
    for line in run_log_text.splitlines():
        m = PASS_RE.search(line)
        if m:
            rows.append(m.groups())
    return rows


def certification(run_log_text):
    """(verdict, [refusing rows]) of the LAST certification block."""
    verdict, refusals = None, []
    for line in run_log_text.splitlines():
        if 'NOT CERTIFIED:' in line:
            verdict, refusals = 'NOT CERTIFIED', []
        elif line.strip().startswith('CERTIFIED:'):
            verdict, refusals = 'CERTIFIED', []
        elif verdict == 'NOT CERTIFIED' and re.match(
                r'^\s{5,}\S.*(above|refuse)', line):
            refusals.append(line.strip())
    return verdict, refusals


def solve_verdict(run_log_text):
    stationary = STATIONARY_RE.findall(run_log_text)
    info = rnorm = None
    for line in run_log_text.splitlines():
        m = DONE_RE.match(line)
        if m:
            info = int(m.group(1))
            if m.group(2):
                rnorm = m.group(2)
    if stationary:
        return int(stationary[-1]), rnorm, 'partitioned stationary route'
    return info, rnorm, 'marching plus the JFNK hand-off'


def equivalent_width(path):
    if not os.path.isfile(path):
        return None
    try:
        import numpy as np
    except ImportError:
        return None
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    exc = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return float(np.trapz(exc[m], lam[m]))


def line_metrics(path):
    out = {}
    for line in read(path).splitlines():
        if line.lstrip().startswith('#') or not line.strip():
            continue
        word = line.split()
        if len(word) >= 2:
            try:
                out[word[0]] = float(word[1])
            except ValueError:
                pass
    return out


def results_block(d, threads, wall, started, finished, run_log='run.log',
                  pp_log='pp.log'):
    text = read(os.path.join(d, run_log))
    out = []
    info, rnorm, route = solve_verdict(text)
    rows = outer_passes(text)
    if rows:
        out.append('The outer passes of the %s, as `%s` records them:\n'
                   % (route, run_log))
        out.append('| pass | hydro info | mass | momentum | energy |'
                   ' worst gated species row | at cell | s |')
        out.append('|---|---|---|---|---|---|---|---|')
        for (k, hinfo, srow, stol, scell, sname,
             mass, mom, energy, secs) in rows:
            out.append('| %s | %s | %s | %s | %s | %s of %s (%s) | %s | %s |'
                       % (k, hinfo, mass, mom, energy, srow, stol, sname,
                          scell, secs))
        out.append('')
    else:
        m = MARCH_RE.search(text)
        if m:
            out.append('The march stopped at step %s with a mass-flux spread'
                       ' `du` = %s.\n' % (m.group(1), m.group(2)))
    out.append('- solver verdict: **info = %s**%s, by the %s'
               % (info, '' if rnorm is None else ' (last `||R||` = %s)' % rnorm,
                  route))
    verdict, refusals = certification(text)
    out.append('- certification of the state written: **%s**'
               % (verdict or 'no verdict printed'))
    for line in refusals:
        out.append('  - %s' % line)
    pp_text = read(os.path.join(d, pp_log))
    mdot = MDOT_RE.findall(pp_text)
    if mdot:
        out.append('- log10 Mdot [g/s] = **%s**, from the post-processing'
                   ' pass' % mdot[-1])
    claim = CLAIM_RE.findall(pp_text)
    work = WORK_RE.findall(pp_text)
    if claim:
        out.append('- the solved state handed back to the post-processing'
                   ' pass: original claim **%s**' % claim[-1])
    if work:
        out.append('- the state that pass evaluated: **%s**' % work[-1])
    ew = equivalent_width(os.path.join(d, 'tpm_He10830.txt'))
    met = line_metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
    if ew is not None:
        out.append('- He I 10830 red-pair equivalent width = **%.4f** %%A'
                   ' over %.2f to %.2f A (air)' % (ew, EW_LO, EW_HI))
    if 'red_depth' in met:
        out.append('- red-pair depth = **%.4f** %%, FWHM = %.4f A'
                   ' (three-Gaussian fit)'
                   % (met['red_depth'], met.get('fwhm_A', float('nan'))))
    out.append('- wall clock **%s** at `OMP_NUM_THREADS=%s` on `%s`, %s to %s'
               % (wall, threads, platform.node(), started, finished))
    return '\n'.join(out)


def header(case, d, binary, threads, kind, logs_dir, note=''):
    setup = read(os.path.join(logs_dir, 'EXHALE_setup.out')).splitlines()
    run = read(os.path.join(logs_dir, 'run.log')).splitlines()
    lapack = [l.strip() for l in run[:40] if 'LAPACK' in l]
    out = ['# How this model was reached',
           '',
           'Written by `models/write_reproduce.py`, called by the runner at'
           ' the end of the run it describes.',
           '',
           '## The case',
           '',
           '`%s`' % case,
           '']
    for name, text in group_fields(case):
        out.append('- %s: %s' % (name, text))
    flux = numerical_flux(d)
    if flux:
        out.append('- numerical flux: `%s`' % flux)
    out.append('')
    if flux == 'ROE':
        out.append(ROE_SENTENCE)
        out.append('')
    if kind == 'closure':
        out.append('The rung has no `input.inp` of its own:'
                   ' `input_template.inp` beside `closure.json` is what each'
                   ' iteration starts from, and the closure driver writes'
                   ' that iteration\'s own `input.inp` into `kNN/` from it,'
                   ' adding the route\'s keys (`input_keys`) and the profile'
                   ' the chemistry step just produced.')
    else:
        out.append('The `input.inp` in this directory is the file the'
                   ' solution was started from; the runner puts it back after'
                   ' the post-processing pass, which needs two of its keys'
                   ' changed.')
    out.append('')
    if note:
        # What this run had to be told, when the run before it did not reach
        # the root on the defaults.  It is the runner's to state, because the
        # logs of the attempt it replaces are no longer in this directory.
        out.append('## Why this run was made the way it was')
        out.append('')
        out.append(note)
        out.append('')
    out.append('## The binary')
    out.append('')
    out.append('| | |')
    out.append('|---|---|')
    out.append('| path | `%s` |' % binary)
    out.append('| md5 | `%s` |' % md5(binary))
    out.append('| repository HEAD | `%s` |' % (shell(['git', 'rev-parse',
                                                      'HEAD']) or 'unknown'))
    dirty = shell(['git', 'status', '--porcelain'])
    out.append('| working tree | %s |'
               % ('dirty (uncommitted changes present)' if dirty else 'clean'))
    out.append('| compiler | %s |' % compiler_of(binary))
    if lapack:
        out.append('| linear algebra | %s |' % lapack[0].lstrip('- '))
    if setup:
        out.append('| run title | `%s` |' % setup[0].strip().strip('#').strip())
    out.append('')
    return '\n'.join(out)


def reproduction_note(threads):
    return '\n'.join([
        '## Reproducing it',
        '',
        'Run the commands above, in this directory, in the order they are'
        ' given; each pass reads what the one before it wrote, so the order'
        ' is the content and not a convenience.',
        '',
        'The OpenMP parallelization of the cell ionization sweep is bitwise'
        ' identical to the serial result (`README.md`, feature list) and the'
        ' BLAS thread pool is pinned to one thread by the run itself, so'
        ' `OMP_NUM_THREADS` is expected to change the wall clock and not the'
        ' answer; the regression harness nonetheless fixes'
        ' `OMP_NUM_THREADS=1` because that is the only setting under which it'
        ' compares outputs bitwise. This run used'
        ' `OMP_NUM_THREADS=%s`. A different binary -- another compiler,'
        ' another BLAS -- reproduces the physics and not the last digits.'
        % threads,
        '',
        'The seed is another solved state, not a fresh one: a certified case'
        ' of `models/` where one of the right physics exists, otherwise an'
        ' archived state out of `archive_20260830/`, which is never'
        ' rewritten. Either way it is interpolated onto the current grid.'
        ' The solve does not depend on which seed of the right physics is'
        ' used -- three seeds four decades apart in the base velocity land on'
        ' the same fixed point to 1e-8'
        ' (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun'
        ' whose `pick_seed.py` chooses differently, because another case has'
        ' been certified since, is still the same solution.',
        '',
        '## Re-measuring this state, and how closely it comes back',
        '',
        'The two files of `output/` can be handed back to the binary as'
        ' `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with'
        ' `Restart intent: stationary evaluate`, which measures the state as'
        ' it stands, writes it back unchanged and exits 0 only if every'
        ' active equation is within its tolerance.',
        '',
        'That is the post-processing pass of this record: it is the route the'
        ' advection-corrected profiles and the mass-loss line above were'
        ' produced on. Three states are named in it and every product says'
        ' which one it describes: the LOADED state, the conserved variables'
        ' the two files carry, which the pass does not write to; the WORK'
        ' state, which holds those conserved variables fixed and derives the'
        ' pressure and the temperature from the composition one equilibrium'
        ' sweep returns, through the caloric equation of state, and which is'
        ' what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the'
        ' `certified=` pair of their header describe; and the'
        ' ADVECTION-DERIVED composition the post-process builds from the work'
        ' state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt`'
        ' carry. The derived pair states the work state\'s certification'
        ' pair on a `# derived_from:` line, as provenance of the state it was'
        ' built from, and makes no certification claim about its own rows.'
        ' One sweep is not a closed chemical and thermal fixed point, so the'
        ' work state reports its own closure defect in `pp.log` and'
        ' "refreshed" never means "closed".',
        '',
        'Two answers are kept apart in `pp.log` and in the list above:'
        ' whether the stationary claim the solved state was written with'
        ' reproduces when that state is handed back, and what verdict the'
        ' work state gets. A passing work state is not a confirmation of a'
        ' claim that did not reproduce.',
        '',
        'What comes back exactly and what does not, MEASURED over the 120'
        ' certified states of `models/`'
        ' (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8).'
        ' The conserved state round-trips to the last bit: the radius,'
        ' velocity and pressure columns are read as written, and since item'
        ' L18 the mass density is read from its own column rather than'
        ' rebuilt from the species, so the state re-measured is the state'
        ' certified. The ROW MEASURES do not round-trip bitwise and cannot:'
        ' the first equilibrium sweep of the re-entry moves the composition'
        ' by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at'
        ' its rounding floor, where one unit in the last place of the density'
        ' moves the cell-wise maximum of the energy row by about 11 per cent.'
        ' Re-evaluated against in-run, the mass row comes back within a'
        ' factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to'
        ' 1.93 (median 1.22); the momentum row sits four to six decades below'
        ' its tolerance and its ratio is the floor itself. A state certified'
        ' at more than about half of its tolerance may therefore be refused'
        ' on re-evaluation, and that refusal is a property of the'
        ' cancellation and not of this file.',
        ''])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('case_dir')
    ap.add_argument('--kind', choices=('case', 'closure'), default='case')
    ap.add_argument('--note', default='',
                    help='one paragraph the runner wants in the record: what'
                         ' this run had to be told, and why')
    ap.add_argument('--case-name', default=None,
                    help='<group>/<case>; taken from the path when omitted')
    ap.add_argument('--binary', default=os.path.join(EXHALE, 'EXHALE.x'))
    ap.add_argument('--threads', default='8')
    ap.add_argument('--seed-line', default='',
                    help="pick_seed.py's own line, or 'cold start'")
    ap.add_argument('--seed-command', default='',
                    help='the map_state_to_grid.py command as it was given')
    ap.add_argument('--seed-rescale', default='',
                    help='"<old> <new> <factor>" when the seed reservoir He/H'
                         ' was carried onto this case\'s')
    ap.add_argument('--grid-target', default='',
                    help='the file whose "# grid" line the mapped seed carries')
    ap.add_argument('--grid-r0', default='',
                    help='R0 [cm] that file states, when the case builds its'
                         ' grid on a lower-atmosphere profile')
    ap.add_argument('--commands', default='',
                    help='the commands of the run, one per line, with their'
                         ' environment, ready to paste')
    ap.add_argument('--started', default='')
    ap.add_argument('--finished', default='')
    ap.add_argument('--wall', default='')
    ap.add_argument('--iterations', default='',
                    help='closure only: the iteration directories, space'
                         ' separated')
    a = ap.parse_args()

    d = os.path.abspath(a.case_dir)
    case = a.case_name or os.path.join(os.path.basename(os.path.dirname(d)),
                                       os.path.basename(d))
    iters = a.iterations.split()
    logs_dir = os.path.join(d, iters[-1]) if (a.kind == 'closure' and iters) \
        else d
    out = [header(case, d, a.binary, a.threads, a.kind, logs_dir, a.note)]

    out.append('## The seed')
    out.append('')
    molecular_seed = bool(a.seed_line) \
        and a.seed_line.lstrip().startswith('molecular seed from')
    if molecular_seed:
        # A molecular case has no archived state of its physics and pick_seed
        # has no list for it: the binary builds the state out of the
        # certified atomic case of the same name and He/H.
        out.append('The binary built this case\'s state out of a certified'
                   ' ATOMIC solution, which is how a molecular case is'
                   ' seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py`'
                   ' carries no molecular state to choose from):')
        out.append('')
        out.append('```')
        out.append(a.seed_line)
        out.append('```')
        out.append('')
        out.append('The fields are the atomic state directory and which'
                   ' extension of the base H2 partition the conversion'
                   ' carried: `local` gives every cell its own'
                   ' chemical-equilibrium q_H2(p,T), `handoff` carries the'
                   ' base x2 to the outer boundary, a number is that'
                   ' fraction everywhere.')
        if a.seed_command:
            out.append('')
            out.append('The conversion, which writes `output/*_IC.txt` and'
                       ' stops:')
            out.append('')
            out.append('```')
            out.append(a.seed_command)
            out.append('```')
    elif a.seed_line and a.seed_line != 'cold start':
        out.append('`models/pick_seed.py` chose, out of the states that carry'
                   ' this case\'s physics:')
        out.append('')
        out.append('```')
        out.append(a.seed_line)
        out.append('```')
        out.append('')
        out.append('The fields are the state directory, the tier it was taken'
                   ' from (`tier0` this same case\'s own most recent'
                   ' certified state -- the one this tree carries, else the'
                   ' one a preserved tree such as'
                   ' `models_20260914_preL21/` holds, which is what a'
                   ' re-solve of the whole catalog continues from --'
                   ' `tier1` a certified case of the same group,'
                   ' `tier2` a certified case of another group with the same'
                   ' physics, `tier3` the same at another XUV normalization,'
                   ' `tier4` and `tier5` the archive) with the group or the'
                   ' code generation it belongs to, the seed reservoir He/H,'
                   ' this case\'s, the distance'
                   ' |log10 He/H_seed - log10 He/H_case| the choice'
                   ' minimizes, and how many states of that tier carried the'
                   ' physics at all.')
        if a.seed_rescale:
            old, new, fac = a.seed_rescale.split()
            out.append('')
            out.append('The seed was solved at He/H = %s and this case is at'
                       ' %s, so `--reservoir He/H %s` carried its helium onto'
                       ' this composition: every helium column of every row'
                       ' was multiplied by %s, which sets the base rows --'
                       ' where the element-diffusion operator holds its'
                       ' Dirichlet He/H -- to %s exactly and leaves the shape'
                       ' of the diffused He/H profile as it was. The pressure'
                       ' is unchanged; the mass density and the temperature'
                       ' were rewritten to follow the new particle count, and'
                       ' the `# reservoir` line of the seed states %s, which'
                       ' is what `load_IC` compares against the input.'
                       % (old, new, new, fac, new, new))
        if a.grid_r0:
            out.append('')
            out.append('The grid the seed is written on is named by'
                       ' `%s`. This case builds its grid on a'
                       ' lower-atmosphere profile, so its radial scale R0 is'
                       ' the radius that profile carries at its matching'
                       ' level times R_J, %s cm, and not the planet radius'
                       ' the shared `current_grid_Hydro_ioniz.txt` states.'
                       ' `load_IC` compares the `grid` metadata field as'
                       ' text, so the seed has to carry this number: the file'
                       ' is that shared one with the single field replaced,'
                       ' written by `src/utils/profile_match_level.py`.'
                       % (a.grid_target or 'target_grid_Hydro_ioniz.txt',
                          a.grid_r0))
        elif a.grid_target:
            out.append('')
            out.append('The grid the seed is written on is named by `%s`.'
                       % a.grid_target)
        if a.seed_command:
            out.append('')
            out.append('It was interpolated onto the cell centers of the'
                       ' current code by')
            out.append('')
            out.append('```')
            out.append(a.seed_command)
            out.append('```')
    else:
        out.append('None: this case starts from the code\'s own initial'
                   ' condition (`Load IC? False`). The archive holds no state'
                   ' of its physics.')
    out.append('')

    out.append('## The commands')
    out.append('')
    out.append('```bash')
    out.append(a.commands.strip() or '(not recorded)')
    out.append('```')
    out.append('')

    out.append('## What came out')
    out.append('')
    if a.kind == 'closure':
        for k in iters:
            out.append('### Iteration `%s`' % k)
            out.append('')
            out.append(results_block(os.path.join(d, k), a.threads,
                                     a.wall if k == iters[-1] else '(see the'
                                     ' whole rung)', a.started, a.finished))
            out.append('')
        history = os.path.join(d, 'closure_history.txt')
        if os.path.isfile(history):
            out.append('### The rung')
            out.append('')
            out.append('`closure_history.txt`, last row: the converged trial'
                       ' fluxes, the reservoir He/H at the matching level and'
                       ' log10 Mdot.')
            out.append('')
            out.append('```')
            rows = [l for l in read(history).splitlines()
                    if l.strip() and not l.startswith('#')]
            out.append(rows[-1] if rows else '(no data row)')
            out.append('```')
            out.append('')
    else:
        out.append(results_block(d, a.threads, a.wall, a.started, a.finished))
        out.append('')

    out.append(reproduction_note(a.threads))

    path = os.path.join(d, 'REPRODUCE.md')
    with open(path, 'w') as fh:
        fh.write('\n'.join(out).rstrip() + '\n')
    print('wrote %s' % path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
