#!/usr/bin/env python3
"""What a lower-atmosphere profile states AT THE MATCHING LEVEL, and the
radial scale R0 the run builds its grid on from it.

One implementation, two callers: `element_flux_closure.py` (the elemental-flux
closure, seeding one iteration from the previous one) and
`LHS1140b/models/run_case.sh` (a case seeded from another solved state). Both
need the same number, because both hand a state to `load_IC`, which compares
the `grid` metadata field as TEXT (`load_IC.f90::refuse_unequal_field`) and
refuses a state whose R0 is not the run's.

HOW THE RUN FORMS R0. With `Lower atmosphere profile:` in use,
`input_read::apply_lower_atmosphere_profile` sets `R0` to the profile's `r` at
the matching level, through `lap_value_at_match('r')` ->
`lower_atmosphere_profile::value_at_pressure(ic_r, lap_p_match_bar)`: a level
whose own `p` equals `p_match` returns that level's value bit for bit, and
otherwise the value is linear in ln p between the two levels that bracket it
(the levels run deep to shallow). Nothing else in `src/modules/` writes `R0`
until `input_read.f90:1767`, `R0 = R0*RJ`. So R0 in cm is ONE double multiply
on that radius, and the two operations here are the same two, in the same
order and the same precision, with RJ the value `parameters.f90:835` carries.

Without a profile the run takes R0 from `Planet radius [R_J]` instead and this
module has nothing to say; the caller then uses the grid file as it stands.

THE ELEMENTAL RESERVOIRS COME FROM THE SAME LEVEL, and for the same reason:
`load_IC` compares the `# reservoir` line element by element and refuses a
state whose ratio differs from the run's by more than 1e-6 IN ANY ELEMENT,
and with a profile in use the run reads each ratio as the `X_<El>` column at
the matching level (`input_read::lap_element_ratio_at_match`). So a seed
solved under one column and reloaded under another has to be carried onto
every element that moved, not onto He/H alone.

As a command, either

    profile_match_level.py <profile.dat> <grid_file> <out_file>

which writes <out_file>, a copy of <grid_file> whose `# grid` line states the
R0 of <profile.dat>, and prints that R0 as the metadata block writes it
(exit 2 if the profile carries no matching level or no radius column, 3 if
the grid file carries no `# grid` line), or

    profile_match_level.py --reservoir-options <state_file> <profile.dat>

which prints the `--reservoir <El>/H <value>` arguments that carry the state
in <state_file> onto the composition <profile.dat> states, one line per
element that moved, and nothing when none did.
"""
import math
import os
import sys

# Jupiter equatorial radius (cm), IAU 2015 nominal; parameters.f90:835.
RJ_CM = 7.1492e9


def profile_values_at_match(path):
    """Every column of a lower-atmosphere profile at the matching level, as
    {name: value}, by the rule `value_at_pressure` applies (see the module
    docstring). The matching pressure is the `# p_match_bar` header and the
    column names are the `# columns:` header. None when the file carries no
    `# columns:` line, no `# p_match_bar`, no data rows, or no `p` column."""
    names, p_match, rows = None, None, []
    try:
        with open(path) as fh:
            for line in fh:
                if line.lstrip().startswith('#'):
                    word = line.split()
                    if len(word) > 1 and word[1].rstrip(':') == 'columns':
                        names = [w.rstrip(':') for w in word[2:]]
                    elif len(word) > 2 and word[1] == 'p_match_bar':
                        p_match = float(word[2])
                    continue
                if line.strip():
                    rows.append([float(w) for w in line.split()])
    except (OSError, ValueError, IndexError):
        return None
    if not names or p_match is None or not rows or 'p' not in names:
        return None
    jp = names.index('p')
    names = names[:min(len(names), min(len(r) for r in rows))]

    def at(j):
        for row in rows:                      # a level that IS the level
            if row[jp] == p_match:
                return row[j]
        for i in range(len(rows) - 1):        # the bracketing pair
            if rows[i][jp] > p_match > rows[i+1][jp]:
                wgt = ((math.log(p_match) - math.log(rows[i][jp]))
                       / (math.log(rows[i+1][jp]) - math.log(rows[i][jp])))
                return rows[i][j] + wgt*(rows[i+1][j] - rows[i][j])
        return rows[0][j] if p_match >= rows[0][jp] else rows[-1][j]

    return {nm: at(j) for j, nm in enumerate(names)}


def profile_match_ratios(path):
    """The elemental reservoirs the run will ask for, {'He/H': v, ...}, from
    the `X_<El>` columns of `profile_values_at_match`: that is how
    `input_read` reads them (`lap_element_ratio_at_match`, the column
    "X_<El>" at the matching level), so these are the values the restart
    contract compares the seed's `# reservoir` line against."""
    vals = profile_values_at_match(path)
    if vals is None:
        return None
    return {nm[2:] + '/H': v for nm, v in vals.items() if nm.startswith('X_')}


def reservoir_ratios(state_file):
    """The `# reservoir <El>/H <value> ...` header of a state file, as
    {'He/H': v, 'C/H': v, ...}; None when the file carries no metadata block
    (an archived state)."""
    try:
        with open(state_file, errors='replace') as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                word = line.split()
                if len(word) > 3 and word[1] == 'reservoir':
                    out = {}
                    for i in range(2, len(word) - 1, 2):
                        try:
                            out[word[i]] = float(word[i+1])
                        except ValueError:
                            pass
                    return out
    except OSError:
        return None
    return None


# The restart contract's own allowance, `load_IC.f90`'s heh_dev_tol.
RESERVOIR_TOL = 1.0e-6


def reservoir_moves(state_file, ratios, tol=RESERVOIR_TOL):
    """{El/H: (seed value, run value)} for every element whose ratio the run
    asks for and the state does not carry within `tol`."""
    seed = reservoir_ratios(state_file)
    moved = {}
    for key, old in sorted((seed or {}).items()):
        new = (ratios or {}).get(key)
        if new is None or not (old > 0.0) or not (new > 0.0):
            continue
        if abs(new/old - 1.0) > tol:
            moved[key] = (old, new)
    return moved


def grid_R0_from_profile(path):
    """R0 in cm, as a run with this profile builds its grid on; None when the
    profile states no radius at the matching level."""
    vals = profile_values_at_match(path)
    if not vals or 'r' not in vals:
        return None
    return vals['r']*RJ_CM


def meta_num(x):
    """A number as the restart metadata block writes it: `load_IC.f90`'s
    `meta_num`, an `ES23.16` field left-adjusted, 17 significant digits. The
    `grid` field is compared as text, so a target header has to carry these
    exact characters."""
    return '%.16E' % x


def grid_R0_cm(state_file):
    """The R0 in cm that a state or grid file's `# grid` line states; None
    when the file carries no such line."""
    try:
        with open(state_file, errors='replace') as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                word = line.split()
                if len(word) > 2 and word[1] == 'grid' and 'R0[cm]' in word:
                    return float(word[word.index('R0[cm]') + 1])
    except (OSError, ValueError, IndexError):
        return None
    return None


def write_target_grid_header(grid_file, target, R0_cm):
    """The target `map_state_to_grid.py` is pointed at: a copy of `grid_file`
    whose `# grid` line states `R0_cm`.

    The mapper reads the target's radius column and copies its `# grid` line
    into the state it writes, and `load_IC` compares that line against the
    run's as text. R0 is the profile's radius at the matching level, so it is
    the run's own and neither the seed's nor the shared grid file's: a seed
    mapped against either would carry the wrong one and be refused. Everything
    else on the line -- N, r_min, r_max and the grid mode -- is the grid
    file's, because the grid is built in units of R_p and does not move with
    R0. Were one of them to move as well, `load_IC` would refuse the state and
    name both lines.

    Returns the R0 as it was written, or None when the file carries no
    `# grid` line."""
    out, wrote = [], None
    with open(grid_file) as fh:
        for line in fh:
            if line.startswith('# grid ') and 'R0[cm]' in line:
                word = line.split()
                wrote = meta_num(R0_cm)
                word[word.index('R0[cm]') + 1] = wrote
                line = ' '.join(word) + '\n'
            out.append(line)
    if wrote is None:
        return None
    tmp = target + '.tmp'
    with open(tmp, 'w') as fh:
        fh.writelines(out)
    os.replace(tmp, target)
    return wrote


def main():
    if len(sys.argv) == 4 and sys.argv[1] == '--reservoir-options':
        state_file, profile = sys.argv[2:4]
        moved = reservoir_moves(state_file,
                                profile_match_ratios(profile))
        for key, (old, new) in sorted(moved.items()):
            print('--reservoir %s %.17g' % (key, new))
        return
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    profile, grid_file, target = sys.argv[1:4]
    R0 = grid_R0_from_profile(profile)
    if R0 is None:
        sys.stderr.write('profile_match_level: %s states no radius at a '
                         'matching level\n' % profile)
        sys.exit(2)
    wrote = write_target_grid_header(grid_file, target, R0)
    if wrote is None:
        sys.stderr.write('profile_match_level: %s carries no "# grid" line\n'
                         % grid_file)
        sys.exit(3)
    print(wrote)


if __name__ == '__main__':
    main()
