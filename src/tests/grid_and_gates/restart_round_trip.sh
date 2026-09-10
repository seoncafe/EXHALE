#!/bin/bash
# Writer-to-loader round trip of a restart state, with no evolution between
# the write and the read (PLAN_20260909_rev1 item N9, review F13).
#
# WHAT A ROUND TRIP CAN AND CANNOT BE
#   The writer emits dimensional quantities in a text profile; the loader
#   reads them, rebuilds the mass density from the species with the run's mass
#   policy, and divides back into mass fractions.  Serialization and
#   reconstruction are therefore two different statements and are asserted
#   separately.  The rows below are, in order: the stored scalar fields
#   survive the writer's format; the conserved state and the inventories
#   survive the reconstruction to a round-off budget; the mass density the
#   loader hands on IS the mass policy of the species columns it read (the
#   writer's own rho column is not read, and the distance between the two is
#   the writer state's mass-closure defect, measured and printed here); no
#   chemistry, remap, boundary reinitialization or clock advance happens on
#   load; a second round trip adds no drift.
#
# WHAT IS RUN (six short runs; nothing in backup/regression is written to)
#   A  the hot-Uranus molecular gate of backup/regression/roundtrip, 40 steps,
#      cold start.  Its Hydro_ioniz.txt / Ion_species.txt are the reference
#      state: the molecular network has solved and the state carries HeH+, so
#      the composition path with a species shared between two elements is the
#      one under test.
#   B  A's two files handed back as the _IC pair, EXHALE_DUMP_IC=1, which
#      writes the loaded state before the first equilibrium sweep and stops.
#   C  B's output handed back the same way: the second round trip.
#   D  A's state reloaded with EXHALE_RESIDUAL=1 and EXHALE_RELOAD_EQ=0, which
#      evaluates the steady residual of the state as read.
#   E  the same on B's output: the residual of the state after one round trip.
#   F  A's state with the helium of ONE interior cell raised by half, so its
#      He/H departs from the input at that cell and nowhere else.  The loader
#      refuses it (HeH+ cannot be rescaled by one factor) and the refusal has
#      to name the cell that decided the comparison.
#
# THE ROUND-OFF BUDGET
#   ulp_tol below is 1e-14, about forty-five units in the last place of a
#   double.  MEASURED on this fixture the worst distance of any species
#   column or inventory across one round trip is a few units in the last
#   place; the allowance is loose enough that a different summation order in
#   the mass policy cannot fail the row and tight enough that any arithmetic
#   on the state -- a chemistry step, a boundary rebuild, a remap -- does.
#
# Usage: restart_round_trip.sh   (EXHALE_EXE selects the binary)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
WORK="$OUT/restart_round_trip"
CASE="$ROOT/backup/regression/roundtrip"
JCELL=199        # the physical cell whose helium stage F raises

if [ ! -x "$EXE" ]; then
   echo "FAIL restart_round_trip measured=no_binary reference=$EXE tol=0"
   exit 1
fi
if [ ! -f "$CASE/input.inp" ]; then
   echo "FAIL restart_round_trip measured=no_case reference=$CASE tol=0"
   exit 1
fi

rm -rf "$WORK"
for d in A B C D E F; do
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/$d/"
done
for d in B C D E F; do
   sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/$d/input.inp"
done

run_stage() {   # run_stage <dir> <name> <env assignments...>
   local d="$1"; local name="$2"; shift 2
   ( cd "$WORK/$d" && env OMP_NUM_THREADS=1 "$@" "$EXE" > run.log 2>&1 )
   local rc=$?
   # Exit 2 is a stationary claim the certification refused; the outputs are
   # written in full and that is what this gate reads.
   if [ "$name" != "refusal" ] && [ $rc -ne 0 ] && [ $rc -ne 2 ]; then
      echo "FAIL restart_round_trip_stage_$name measured=exit_$rc reference=0 tol=0"
      tail -n 5 "$WORK/$d/run.log"
      return 1
   fi
   RC_LAST=$rc
   return 0
}

run_stage A stageA EXHALE_MAXSTEPS=40 || exit 1
for f in Hydro_ioniz Ion_species; do
   if [ ! -s "$WORK/A/output/$f.txt" ]; then
      echo "FAIL restart_round_trip_stageA measured=no_$f reference=written tol=0"
      exit 1
   fi
done

cp "$WORK/A/output/Hydro_ioniz.txt" "$WORK/B/output/Hydro_ioniz_IC.txt"
cp "$WORK/A/output/Ion_species.txt" "$WORK/B/output/Ion_species_IC.txt"
run_stage B stageB EXHALE_DUMP_IC=1 || exit 1

cp "$WORK/B/output/Hydro_ioniz.txt" "$WORK/C/output/Hydro_ioniz_IC.txt"
cp "$WORK/B/output/Ion_species.txt" "$WORK/C/output/Ion_species_IC.txt"
run_stage C stageC EXHALE_DUMP_IC=1 || exit 1

cp "$WORK/A/output/Hydro_ioniz.txt" "$WORK/D/output/Hydro_ioniz_IC.txt"
cp "$WORK/A/output/Ion_species.txt" "$WORK/D/output/Ion_species_IC.txt"
run_stage D stageD EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 || exit 1

cp "$WORK/B/output/Hydro_ioniz.txt" "$WORK/E/output/Hydro_ioniz_IC.txt"
cp "$WORK/B/output/Ion_species.txt" "$WORK/E/output/Ion_species_IC.txt"
run_stage E stageE EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0 || exit 1

# Stage F's state: one interior cell's helium raised by half.
cp "$WORK/A/output/Hydro_ioniz.txt" "$WORK/F/output/Hydro_ioniz_IC.txt"
python3 - "$WORK/A/output/Ion_species.txt" \
           "$WORK/F/output/Ion_species_IC.txt" "$JCELL" <<'PY'
import sys
src, dst, jcell = sys.argv[1], sys.argv[2], int(sys.argv[3])
labels, out, irow = None, [], 0
for line in open(src):
    s = line.strip()
    if s.startswith('#') or not s:
        # The field of a header line is its FIRST token: 'species_columns' of
        # the restart metadata block also contains 'columns', and reading it
        # as the label line shifts every column index by one.
        if s.lstrip('#').split()[:1] == ['columns']:
            labels = s.replace('#', '', 1).split()[1:]
        out.append(line)
        continue
    irow += 1
    f = [float(x) for x in s.split()]
    if irow == jcell + 2:            # two ghost rows precede physical cell 1
        f[labels.index('HeI')] *= 1.5
    out.append('  ' + '   '.join(repr(x) for x in f) + '\n')
open(dst, 'w').writelines(out)
PY
run_stage F refusal
rcF=$RC_LAST

python3 - "$WORK" "$JCELL" "$rcF" <<'PY'
import re
import sys

import numpy as np

work, jcell, rcF = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
ulp_tol = 1.0e-14
# m(He)/m(H) of the mass policy: m_He_atom / mu of parameters.f90, the ratio
# calc_rho weights every helium-bearing species with.  H2 and H2+ carry two
# hydrogen masses, H3+ three, HeH+ one of each.
M_HE = 6.6464790722e-24 / 1.67353284e-24
n_fail = 0


def read(path):
    labels, rows = None, []
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            if s.lstrip('#').split()[:1] == ['columns']:
                labels = s.replace('#', '', 1).split()[1:]
            continue
        if s:
            rows.append([float(x) for x in s.split()])
    return labels, np.array(rows)


def col(labels, a, name):
    return a[:, labels.index(name)] if name in labels else np.zeros(a.shape[0])


def maxrel(a, b):
    s = np.maximum(np.abs(a), np.abs(b))
    d = np.abs(a - b)
    m = s > 0
    if not m.any():
        return 0.0, 0
    rel = np.zeros_like(d)
    rel[m] = d[m] / s[m]
    j = int(np.argmax(rel))
    return float(rel[j]), j + 1


def verdict(name, measured, reference, tol):
    global n_fail
    ok = measured <= tol
    if not ok:
        n_fail += 1
    print('%s %s measured=%.6e reference=%.6e tol=%.2e'
          % ('PASS' if ok else 'FAIL', name, measured, reference, tol))


def inventories(labels, a):
    c = lambda n: col(labels, a, n)
    nH = (c('HI') + c('HII') + 2 * (c('H2') + c('H2p')) + 3 * c('H3p')
          + c('HeHp') + c('OH') + 2 * c('H2O'))
    nHe = c('HeI') + c('HeII') + c('HeIII') + c('HeHp')
    nO = c('OI') + c('OII') + c('OIII') + c('OH') + c('H2O') + c('CO')
    nC = c('CI') + c('CII') + c('CIII') + c('CO')
    q = (c('HII') + c('HeII') + 2 * c('HeIII') + c('H2p') + c('H3p')
         + c('HeHp'))
    for e in ('C', 'O', 'N', 'Mg', 'Si', 'Ca', 'Na', 'K', 'S', 'Fe'):
        q = q + col(labels, a, e + 'II') + 2 * col(labels, a, e + 'III')
    inv = {'nH': nH, 'nHe': nHe, 'nO': nO, 'nC': nC, 'charge': q}
    inv['He/H'] = np.where(nH > 0, nHe / np.maximum(nH, 1e-300), 0.0)
    inv['x_H2'] = np.where(nH > 0, 2 * c('H2') / np.maximum(nH, 1e-300), 0.0)
    return inv


def mass_policy(labels, a):
    c = lambda n: col(labels, a, n)
    return (c('HI') + c('HII') + M_HE * (c('HeI') + c('HeII') + c('HeIII'))
            + 2 * c('H2') + 2 * c('H2p') + 3 * c('H3p')
            + (1.0 + M_HE) * c('HeHp'))


def metals_present(labels, a):
    tot = 0.0
    for lab in labels:
        if lab in ('r[Rp]', 'HI', 'HII', 'HeI', 'HeII', 'HeIII', 'HeITR',
                   'H2', 'H2p', 'H3p', 'HeHp'):
            continue
        tot += float(np.abs(col(labels, a, lab)).sum())
    return tot > 0.0


def trip(ref_dir, dump_dir, tag):
    """One round trip: the worst distance over every compared quantity."""
    lh_a, HA = read('%s/output/Hydro_ioniz.txt' % ref_dir)
    lh_b, HB = read('%s/output/Hydro_ioniz.txt' % dump_dir)
    li_a, IA = read('%s/output/Ion_species.txt' % ref_dir)
    li_b, IB = read('%s/output/Ion_species.txt' % dump_dir)
    if HA.shape != HB.shape or IA.shape != IB.shape or li_a != li_b:
        print('FAIL restart_round_trip_%s measured=schema_change '
              'reference=same tol=0' % tag)
        return None
    stored, derived, inv_worst = 0.0, 0.0, 0.0
    for name in ('r[Rp]', 'v[cm/s]', 'p[cgs]'):
        d, j = maxrel(col(lh_a, HA, name), col(lh_b, HB, name))
        print('  %-6s %-12s max_rel %.3e at row %d' % (tag, name, d, j))
        stored = max(stored, d)
    d, j = maxrel(col(lh_a, HA, 'T[K]'), col(lh_b, HB, 'T[K]'))
    print('  %-6s %-12s max_rel %.3e at row %d  (rebuilt from p and the '
          'particle count)' % (tag, 'T[K]', d, j))
    derived = max(derived, d)
    for name in li_a:
        if name == 'r[Rp]':
            continue
        d, j = maxrel(col(li_a, IA, name), col(li_b, IB, name))
        derived = max(derived, d)
    ia, ib = inventories(li_a, IA), inventories(li_b, IB)
    for k in ia:
        d, j = maxrel(ia[k], ib[k])
        print('  %-6s inventory %-8s max_rel %.3e at row %d' % (tag, k, d, j))
        inv_worst = max(inv_worst, d)
    return {'stored': stored, 'derived': max(derived, inv_worst),
            'H_ref': HA, 'lh_ref': lh_a, 'I_ref': IA, 'li_ref': li_a,
            'H_dump': HB, 'lh_dump': lh_b, 'I_dump': IB}


# ---- rows 1 and 2: serialization and reconstruction, one round trip -------
print('---- one round trip: A written, B loaded and re-written ----')
t1 = trip(work + '/A', work + '/B', 'trip1')
if t1 is None:
    sys.exit(1)
verdict('restart_stored_scalar_fields_round_trip', t1['stored'], 0.0, 0.0)
verdict('restart_state_and_inventories_within_budget', t1['derived'], 0.0,
        ulp_tol)

# ---- row 3: the loader's mass density IS the mass policy of the species --
lh_a, HA, li_a, IA = t1['lh_ref'], t1['H_ref'], t1['li_ref'], t1['I_ref']
lh_b, HB = t1['lh_dump'], t1['H_dump']
if metals_present(li_a, IA):
    print('  restart_density_is_the_mass_policy: NOT EVALUATED -- the fixture '
          'carries metal or oxygen columns, whose atomic weights this check '
          'does not hold')
else:
    pol = mass_policy(li_a, IA)
    d_loader, j1 = maxrel(col(lh_b, HB, 'rho[mH/cm3]'), pol)
    d_writer, j2 = maxrel(col(lh_a, HA, 'rho[mH/cm3]'), pol)
    print('  the writer state\'s own mass closure: max |rho_col / '
          'mass_policy(species) - 1| = %.3e at row %d' % (d_writer, j2))
    print('  (that distance, and not a serialization loss, is why the two '
          'rho columns differ: the loader does not read the rho column)')
    verdict('restart_density_is_the_mass_policy', d_loader, 0.0, ulp_tol)

# ---- row 4: no evolution on load -----------------------------------------
print('---- the steady residual of the state as read, before and after '
      'the round trip ----')
ROW_RE = re.compile(r'^\s+(mass|momentum|energy)\s+([0-9.E+-]+)\s+(\d+)\s+'
                    r'([0-9.E+-]+)\s+([0-9.E+-]+)')


def residual_rows(path):
    rows = {}
    for line in open(path):
        m = ROW_RE.match(line)
        if m and m.group(1) not in rows:
            rows[m.group(1)] = (float(m.group(2)), int(m.group(3)),
                                float(m.group(5)))
    return rows


rD = residual_rows(work + '/D/run.log')
rE = residual_rows(work + '/E/run.log')
if len(rD) != 3 or len(rE) != 3:
    print('FAIL restart_no_evolution_on_load measured=no_residual_block '
          'reference=3_rows tol=0')
    n_fail += 1
else:
    worst_max, worst_cell, worst_int = 0.0, 0, 0.0
    for k in ('mass', 'momentum', 'energy'):
        mx_d, cell_d, int_d = rD[k]
        mx_e, cell_e, int_e = rE[k]
        rel = abs(mx_e - mx_d) / max(abs(mx_d), abs(mx_e), 1e-300)
        rel_i = abs(int_e - int_d) / max(abs(int_d), abs(int_e), 1e-300)
        print('  %-9s max %.4E cell %3d -> max %.4E cell %3d   '
              'integrated %.4E -> %.4E' % (k, mx_d, cell_d, mx_e, cell_e,
                                           int_d, int_e))
        worst_max = max(worst_max, rel)
        worst_cell += 0 if cell_d == cell_e else 1
        worst_int = max(worst_int, rel_i)
    # The row asserts the maximum of each residual row and the cell it sits
    # in.  The window integrals are DIAGNOSTIC: the assembled residual of a
    # relaxation snapshot is not a Lipschitz function of the state in the
    # odd-even band above the base, so a round-off difference in the
    # composition moves the integral of that window while leaving every
    # maximum where it was (MEASURED: one unit in the last place on every
    # species column of the hot-Uranus gate moves the 1.03-1.10 mass window
    # by 4 percent and flips the sign of the cell residual at r = 1.075).
    print('  DIAGNOSTIC largest relative movement of a window integral: '
          '%.3e' % worst_int)
    verdict('restart_no_evolution_residual_max', worst_max, 0.0, 1.0e-6)
    verdict('restart_no_evolution_residual_cell', float(worst_cell), 0.0, 0.0)

# ---- row 5: a second round trip adds no drift ----------------------------
print('---- the second round trip ----')
t2 = trip(work + '/B', work + '/C', 'trip2')
if t2 is None:
    sys.exit(1)
verdict('restart_second_trip_stored_fields', t2['stored'], 0.0, 0.0)
verdict('restart_second_trip_no_drift',
        max(t2['derived'] - t1['derived'], 0.0), 0.0, ulp_tol)

# ---- row 6: the refusal names the cell that decided the comparison -------
print('---- the refusal message ----')
msg = open(work + '/F/run.log').read()
named = ('farthest cell %d' % jcell) in msg
refused = rcF != 0
for line in msg.splitlines():
    if 'load_IC' in line or 'farthest cell' in line or 'He/H in the' in line \
            or 'relative departure' in line:
        print('  %s' % line.rstrip())
verdict('restart_refusal_is_taken', 0.0 if refused else 1.0, 0.0, 0.0)
verdict('restart_refusal_names_the_deviating_cell',
        0.0 if named else 1.0, 0.0, 0.0)

sys.exit(1 if n_fail else 0)
PY
exit $?
