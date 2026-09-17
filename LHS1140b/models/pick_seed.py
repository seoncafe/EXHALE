#!/usr/bin/env python3
"""Choose the solved state that seeds one case of `models/`.

The stationary route of `MODELS.md` section 6 starts from a solved state of
the same physics mapped onto the current grid: the solve is a Newton solve
and needs a state in the basin, not the exact answer (L5c measured three
seeds four decades apart in the base velocity landing on the same fixed
point to 1e-8).  What the seed has to share with the case is therefore the
physics that decides which equations are solved -- the spectrum, the kind of
lower boundary, whether the elements diffuse and at which eddy coefficient --
and among the states that share it the nearest composition is taken.

    models/pick_seed.py <group>/<case> [--archive DIR] [--path-only]
                                       [--list N]

writes one line naming the chosen state's `output/` directory, the tier it
was taken from, its reservoir He/H and the distance to the case's:

    <output dir>  <tier>  HeH=<seed>  target=<case>  dlog10=<d>  candidates=<n>

Tiers, in the order they are searched.  A state this code certified is a
better seed than one the code of 2026-08-30 wrote, whatever its composition:
every archived state is the answer to equations that have since changed
(`MODELS.md` section 5), while a certified case of `models/` is a fixed
point of the binary the new case will be solved with.

  0  THE SAME CASE'S OWN MOST RECENT CERTIFIED STATE: the one in
     `models/<group>/<case>/output` where this tree already carries one (a
     flux-closure rung: its last certified iterate, `k<NN>/output`), and
     otherwise the one the preserved tree `models_20260914_preL21/` holds at
     the same relative path.  The campaign is re-solved whenever a fix
     reaches the equations -- the base boundary condition of plan item L21,
     and the review fixes after it -- and what changes each time is a part of
     the system, not the wind: the state this case last certified solves the
     rest of the equations already.  No other state of any tier is that near,
     so this one is taken wherever it exists and certified.
  1  a CERTIFIED case of the SAME GROUP, nearest in |log10 He/H|.  The
     group fixes the chemistry, the lower boundary, the spectrum and the
     mixing, so the ladder differs in composition alone.
  2  a CERTIFIED case of another group that carries the same physics all
     the same -- same spectrum, same kind of lower boundary, same diffusion
     and the same K_zz.
  3  the same at another XUV normalization of the same star and the same
     spectral shape.
  4  an archived state of the same physics, newest generation of the code
     first (`archive_20260830/README.md`): `refresh_j96` before
     `crossings_j96` before `crossings_pc090` before the `*_gm25` series
     before the rest.
  5  an archived state at another XUV normalization.

Within a tier the key is |log10(He/H_seed) - log10(He/H_case)|, the
composition ladder being logarithmic (0.04 to 1000).  A seed whose
reservoir is not the case's is carried onto it by
`src/utils/map_state_to_grid.py --reservoir He/H`, which the runner adds; an
archived state carries no reservoir field and `load_IC` rescales it itself.

CERTIFIED means what the run itself said: the `- certification of the state
written: **CERTIFIED**` line of the case's `REPRODUCE.md`, or, where that
file is absent, a `run.log` whose last stationary verdict is `info = 0` and
whose last certification block is `CERTIFIED:` and not `NOT CERTIFIED:`.

Exit status 3, with a message and no line, when nothing carries the case's
physics: the molecular groups are the case in point, since no archived run
turns molecular chemistry on and no molecular case is certified yet.  The
runner reads that status and starts cold.
"""

import argparse
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))          # .../LHS1140b/models
LHS = os.path.dirname(HERE)                                # .../LHS1140b
ARCHIVE = os.path.join(LHS, 'archive_20260830')
# Where a case's own earlier state is looked for when `models/` no longer
# carries one. `models/archive_results.sh` writes trees of this shape each
# time the catalog is re-solved; they are searched NEWEST FIRST, so add a new
# one at the head of this list when it is written.
#
#   models_20260915_db87    the re-run on the corrected base boundary of item
#                           L21, binary db87b88d1ce5
#   models_20260914_preL21  the campaign the L21 correction superseded
PRESERVED_TREES = [os.path.join(LHS, 'models_20260915_db87'),
                   os.path.join(LHS, 'models_20260914_preL21')]
# The newest, for the messages that name one tree.
PRESERVED = PRESERVED_TREES[0]

# The generations of `archive_20260830/README.md`, newest first.  A path is
# placed by the first tag it contains.
GENERATIONS = ('refresh_j96', 'crossings_j96', 'crossings_pc090', '_gm25')


def generation_rank(path):
    for i, tag in enumerate(GENERATIONS):
        if tag in path:
            return i
    return len(GENERATIONS)


def generation_name(path):
    for tag in GENERATIONS:
        if tag in path:
            return tag.lstrip('_')
    return 'original'


def read_keys(path):
    """`input.inp` as a dict of the keys this tool compares."""
    keys = {}
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' in line:
                name, value = line.split(':', 1)
                keys[name.strip()] = value.strip()
    return keys


def base_inp_reservoirs(case_dir):
    """True when `base.inp` hands over the C, N and O reservoirs.

    The H2 mixing ratio is NOT read here: a `base.inp` may carry
    `q_H2_base` beside the element ratios and mean nothing by it, since the
    binding is used only under a molecular network, which `input.inp`
    decides.
    """
    path = os.path.join(case_dir, 'base.inp')
    if not os.path.isfile(path):
        return False
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('#') \
               and line.split()[0] == 'C_H_base':
                return True
    return False


class Physics(object):
    """The part of a run's configuration a seed has to share with its case."""

    def __init__(self, case_dir):
        inp = os.path.join(case_dir, 'input.inp')
        if not os.path.isfile(inp):
            inp = os.path.join(case_dir, 'input_template.inp')
        self.dir = case_dir
        keys = read_keys(inp)
        self.spectrum = os.path.basename(keys.get('Spectrum file', ''))
        self.molecular = keys.get('Molecular chemistry', 'False').lower() \
            .startswith('true')
        self.diffusion = keys.get('He_diffusion', 'False').lower() \
            .startswith('true')
        kzz = keys.get('He_Kzz')
        self.kzz = float(kzz) if kzz is not None else None
        self.heh = float(keys.get('He/H number ratio', 'nan'))
        if 'Lower atmosphere profile' in keys \
           or os.path.isfile(os.path.join(case_dir, 'closure.json')):
            # A closure rung's own input.inp is written by the driver at each
            # iteration; the template beside closure.json states everything
            # else, and the boundary is the photochemical column.
            self.boundary = 'profile'
        elif base_inp_reservoirs(case_dir):
            self.boundary = 'scalarCNO'
        else:
            self.boundary = 'scalar'

    @property
    def xuv_scale(self):
        """The XUV normalization the spectrum file name carries, 1 when it
        carries none (`..._xuv0p30.txt` -> 0.30)."""
        cut = self.spectrum.find('_xuv')
        if cut < 0:
            return 1.0
        tag = self.spectrum[cut + 4:].split('.')[0]
        try:
            return float(tag.replace('p', '.'))
        except ValueError:
            return 1.0

    @property
    def spectrum_family(self):
        """The spectrum with its XUV scaling dropped: the same star, the same
        shape, a different normalization."""
        name = self.spectrum
        cut = name.find('_xuv')
        return name[:cut] + '.txt' if cut >= 0 else name

    def matches(self, other, same_normalization=True, same_kzz=True):
        """Same equations, same radiation field, same transport coefficient.

        `same_normalization=False` accepts a seed of the same star and the
        same spectral shape at another XUV level, which is a state in the
        basin and not a claim about the answer; the scaled grids of section 3
        have no archived counterpart at every point.
        """
        if same_normalization:
            if self.spectrum != other.spectrum:
                return False
        elif self.spectrum_family != other.spectrum_family:
            return False
        if self.molecular != other.molecular:
            return False
        if self.boundary != other.boundary:
            return False
        if self.diffusion != other.diffusion:
            return False
        # K_zz only where input.inp states it; a profile carries K_zz(p) and
        # states no scalar.
        # A seed at another eddy coefficient is a state of the same basin
        # (K_zz moves the He/H profile, not the wind), and nearer in
        # composition than a distant certified state of the same K_zz.
        if self.boundary != 'profile' and self.diffusion and same_kzz:
            if self.kzz is None or other.kzz is None:
                return False
            if not math.isclose(self.kzz, other.kzz, rel_tol=1e-12,
                                abs_tol=1e-30):
                return False
        return True


def certified(case_dir):
    """Whether the run in this case directory certified the state it wrote.

    `REPRODUCE.md` is the run's own record and is read first; where it is
    absent the verdict is taken from `run.log` the way the runner takes it:
    the last `stationary solve returned info = N` must be 0 and the last
    certification block must be the certifying one.
    """
    rep = os.path.join(case_dir, 'REPRODUCE.md')
    if os.path.isfile(rep):
        with open(rep, errors='replace') as fh:
            for line in fh:
                if line.startswith('- certification of the state written:'):
                    return 'NOT CERTIFIED' not in line and 'CERTIFIED' in line
    log = os.path.join(case_dir, 'run.log')
    if not os.path.isfile(log):
        return False
    with open(log, errors='replace') as fh:
        text = fh.read()
    info = re.findall(r'stationary solve returned info *= *(-?\d+)', text)
    if not info or info[-1] != '0':
        return False
    verdict = None
    for line in text.splitlines():
        if 'NOT CERTIFIED:' in line:
            verdict = False
        elif line.strip().startswith('CERTIFIED:'):
            verdict = True
    return verdict is True


def holds_state(case_dir):
    """Whether this directory carries the pair `load_IC` reads."""
    out = os.path.join(case_dir, 'output')
    return os.path.isfile(os.path.join(out, 'Hydro_ioniz.txt')) \
        and os.path.isfile(os.path.join(out, 'Ion_species.txt'))


def certified_state_of(case_rel, root):
    """This same case's certified state under `root`, or None.

    A prescribed-composition case is one directory holding `output/`; a
    flux-closure rung is a series of iterates `k00/ ... k<NN>/`, of which the
    last certified one is the rung's answer.  The verdict is read the way
    every other tier reads it, from `REPRODUCE.md` where the run wrote one
    and from `run.log` otherwise.
    """
    d = os.path.join(root, case_rel)
    if not os.path.isdir(d):
        return None
    if holds_state(d):
        return d if certified(d) else None
    ks = sorted(k for k in os.listdir(d) if re.match(r'k\d\d$', k))
    for k in reversed(ks):
        kd = os.path.join(d, k)
        if holds_state(kd) and certified(kd):
            return kd
    return None


def own_latest_certified_state(case_rel):
    """The most recent state this same case certified: this tree's own first,
    then the preserved trees', newest first."""
    for root in [HERE] + PRESERVED_TREES:
        d = certified_state_of(case_rel, root)
        if d is not None:
            return d
    return None


def solved_states(models_dir, exclude):
    """Every case of `models/` that holds a certified state, as
    (group, case directory)."""
    found = []
    for group in sorted(os.listdir(models_dir)):
        gdir = os.path.join(models_dir, group)
        if group.startswith('.') or not os.path.isdir(gdir):
            continue
        for case in sorted(os.listdir(gdir)):
            cdir = os.path.join(gdir, case)
            if not os.path.isdir(cdir) or os.path.abspath(cdir) == exclude:
                continue
            out = os.path.join(cdir, 'output')
            if not os.path.isfile(os.path.join(cdir, 'input.inp')):
                continue
            if not (os.path.isfile(os.path.join(out, 'Hydro_ioniz.txt'))
                    and os.path.isfile(os.path.join(out, 'Ion_species.txt'))):
                continue
            if not certified(cdir):
                continue
            found.append((group, cdir))
    return found


def archived_states(archive):
    """Every archived directory that holds both a configuration and a state."""
    found = []
    for root, dirs, files in os.walk(archive):
        if 'input.inp' not in files:
            continue
        out = os.path.join(root, 'output')
        if not (os.path.isfile(os.path.join(out, 'Hydro_ioniz.txt'))
                and os.path.isfile(os.path.join(out, 'Ion_species.txt'))):
            continue
        found.append(root)
    return found


def rank_solved(case, case_group, models_dir, same_group,
                same_normalization=True, same_kzz=True):
    """The certified cases of `models/` that seed this one, best first."""
    exclude = os.path.abspath(case.dir)
    out = []
    for group, cdir in solved_states(models_dir, exclude):
        if (group == case_group) != same_group:
            continue
        try:
            other = Physics(cdir)
        except (ValueError, IOError):
            continue
        if not case.matches(other, same_normalization, same_kzz):
            continue
        if not (other.heh > 0.0) or not (case.heh > 0.0):
            continue
        d = abs(math.log10(other.heh) - math.log10(case.heh))
        dx = abs(math.log10(other.xuv_scale) - math.log10(case.xuv_scale))
        out.append((0, dx, d, cdir, other.heh))
    out.sort(key=lambda row: (row[1], row[2], row[3]))
    return out


def rank(case, archive, same_normalization=True, same_kzz=True):
    """The archived states that carry the case's physics, best first."""
    out = []
    for path in archived_states(archive):
        try:
            other = Physics(path)
        except (ValueError, IOError):
            continue
        if not case.matches(other, same_normalization, same_kzz):
            continue
        if not (other.heh > 0.0) or not (case.heh > 0.0):
            continue
        d = abs(math.log10(other.heh) - math.log10(case.heh))
        # Where the normalization is allowed to differ, the nearest one is
        # taken first; it is zero for every exact match.
        dx = abs(math.log10(other.xuv_scale) - math.log10(case.xuv_scale))
        out.append((generation_rank(path), dx, d, path, other.heh))
    out.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('case', help='<group>/<case>, or a directory')
    ap.add_argument('--archive', default=ARCHIVE)
    ap.add_argument('--path-only', action='store_true',
                    help='write the seed output directory and nothing else')
    ap.add_argument('--list', type=int, default=0, metavar='N',
                    help='also write the N best candidates, as comments')
    a = ap.parse_args()

    case_dir = a.case if os.path.isdir(a.case) else os.path.join(HERE, a.case)
    if not os.path.isdir(case_dir):
        sys.stderr.write('pick_seed: no such case %s\n' % case_dir)
        return 1
    case = Physics(case_dir)
    case_group = os.path.basename(os.path.dirname(os.path.abspath(case_dir)))
    OWN = own_latest_certified_state(os.path.join(
        case_group, os.path.basename(os.path.abspath(case_dir))))

    # The tiers of the docstring, in order; the first that holds anything is
    # the one the seed comes from.  A state this binary certified stands
    # ahead of every archived state, however near in composition.
    tiers = (
        ('tier0', lambda: ([(0, 0.0, 0.0, OWN, case.heh)] if OWN else []),
         "this same case's own most recent certified state"),
        ('tier1', lambda: rank_solved(case, case_group, HERE, True, True),
         'a certified case of the same group'),
        ('tier2', lambda: rank_solved(case, case_group, HERE, False, True),
         'a certified case of another group with the same physics'),
        ('tier2k', lambda: rank_solved(case, case_group, HERE, False, True,
                                       False),
         'a certified case of the same physics at another K_zz'),
        ('tier3', lambda: rank_solved(case, case_group, HERE, False, False),
         'a certified case at another XUV normalization'),
        ('tier4', lambda: rank(case, a.archive, True), 'an archived state'),
        ('tier4k', lambda: rank(case, a.archive, True, False),
         'an archived state at another K_zz'),
        ('tier5', lambda: rank(case, a.archive, False),
         'an archived state at another XUV normalization'),
    )
    # Every tier is ranked, then the candidates of all tiers are ordered by
    # ONE rule: a composition step of more than SEED_STEP_MAX in log10(He/H)
    # ranks after every nearer candidate whatever its tier (2026-09-14: from
    # the certified He/H 3.5 state rescaled onto 0.55, a factor 6.4, the
    # element relaxation found no admissible advance at outer pass 1, while
    # the archived state of the case's own composition solves); within that,
    # the tier order of the docstring; within a tier, the tier's own ranking.
    # The first line is the seed; --list N adds the next N as plain lines of
    # the same shape, so a runner can walk down the list when a seed fails.
    SEED_STEP_MAX = 0.3
    ranked = []
    for ti, (name, build, what) in enumerate(tiers):
        for row in build():
            grank, dx, d, path, heh = row
            ranked.append(((d > SEED_STEP_MAX), ti, grank, dx, d, path, heh,
                           name, what))
    ranked.sort(key=lambda r: (r[0], r[1], r[2], r[3], r[4], r[5]))
    if not ranked:
        sys.stderr.write(
            'pick_seed: no certified case and no archived state carries this'
            ' physics (spectrum %s, %s boundary, diffusion %s, K_zz %s):'
            ' start cold\n'
            % (case.spectrum, case.boundary, case.diffusion, case.kzz))
        return 3
    candidates = ranked

    def origin(path):
        """Where a candidate comes from: the group of a case of models/, the
        code generation of an archived state."""
        p = os.path.abspath(path)
        if p.startswith(os.path.abspath(HERE) + os.sep):
            rel = os.path.relpath(p, HERE)
            return 'models/' + (rel if far_row_is_tier0(p)
                                else os.path.basename(os.path.dirname(path)))
        for root in PRESERVED_TREES:
            if p.startswith(os.path.abspath(root) + os.sep):
                return '%s/%s' % (os.path.basename(root),
                                  os.path.relpath(p, root))
        return 'archive/' + generation_name(path)

    def far_row_is_tier0(p):
        return OWN is not None and p == os.path.abspath(OWN)

    def line(row):
        far, ti, grank, dx, d, path, heh, name, what = row
        return ('%s  %s:%s  HeH=%g  target=%g  dlog10=%.4f  candidates=%d  (%s)'
                % (os.path.join(path, 'output'), name, origin(path), heh,
                   case.heh, d, len(candidates), what))
    if a.path_only:
        print(os.path.join(candidates[0][5], 'output'))
    else:
        print(line(candidates[0]))
    for row in candidates[1:1 + a.list]:
        print(line(row))
    return 0


if __name__ == '__main__':
    sys.exit(main())
