#!/bin/bash
# The restart contract: one intent selector, one metadata block
# (docs/restart_contract_design_20260909.md, PLAN_20260909_rev1 item N10).
#
# WHAT IS UNDER TEST
#   1. The three statements that contradict themselves are refused at
#      startup: an intent without a loaded state, a trajectory continued by a
#      run that claims no elapsed time, a stationary intent with no
#      stationary solver to enter.
#   2. "Restart intent: stationary evaluate" measures the loaded state as it
#      stands and writes it back: no step is taken, the stored fields survive
#      the round trip, the rebuilt quantities stay within the budget the
#      chemistry solve's own tolerance sets, and the certification made on
#      the state is the one the writer made.
#   3. The metadata block round trips: written, loaded, written again, field
#      for field identical.
#   4. A state file pair with no block is loaded as before and marked
#      provenance unknown, in the log and in every file the run writes.
#   5. A block that disagrees with the run refuses the load and names the
#      field: the elemental reservoir (its ratios and its element set), the
#      physical grid, the constant set, the physics options, and a pair whose
#      two halves do not carry the same block.
#   6. The stationary residual is a function of the state: the residual the
#      evaluation prints, the residual the EXHALE_RESIDUAL diagnostic prints
#      for the same state, and the first residual the stationary solve
#      judges are the same number.
#
# THE BUDGETS
#   The stored fields (r, v, p) are text in the file and are read back
#   directly, so they return to a few units in the last place; ulp_tol is
#   1e-14, as in restart_round_trip.sh.  The mass density is not read at all
#   (the loader rebuilds it from the species with the run's mass policy), so
#   its budget is the writer state's own mass closure, 1e-12.  The
#   temperature, the species columns and the heating and cooling come from
#   ONE equilibrium sweep of the loaded composition, and the cell chemistry
#   solve stops at xtol = sqrt(machine epsilon) ~ 1.5e-8 of its own
#   unknowns, which is why equilibrate_loaded_composition calls a particle
#   count reproduced to 1e-10 a fixed point: sweep_tol below is 1e-8, the
#   floor of that solve, and the measured distances are printed beside it.
#
# Usage: restart_intent_and_metadata.sh   (EXHALE_EXE selects the binary)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
WORK="$OUT/restart_intent_and_metadata"
CASE="$ROOT/backup/regression/wasp_full_newton"
ACASE="$ROOT/backup/regression/atomic_elem_newton"
# The stationary solve of the atomic fixture does not converge (it is the N7
# entry point) and only its FIRST judged rows are read here, so that run is
# stopped as soon as they are printed, by the PID this script started, and in
# any case after SOLVE_SECONDS. The loop below breaks as soon as the line is
# there, so the backstop only has to be larger than the slowest machine: it
# costs nothing on an idle one (MEASURED 23 s at 8 threads, under two minutes
# single-threaded) and keeps the row from failing on a loaded one.
SOLVE_SECONDS=600

if [ ! -x "$EXE" ]; then
   echo "FAIL restart_intent measured=no_binary reference=$EXE tol=0"
   exit 1
fi
for f in "$CASE/IC/Hydro_ioniz_IC.txt" "$ACASE/IC/Hydro_ioniz_IC.txt"; do
   if [ ! -s "$f" ]; then
      echo "FAIL restart_intent measured=no_fixture reference=$f tol=0"
      exit 1
   fi
done

rm -rf "$WORK"
mkdir -p "$WORK"

# stage <dir> <state-source-dir> <state-file-suffix> : a run directory whose
# output/ carries the _IC pair copied from a written state.
make_run() {   # make_run <dir> <hydro file> <ion file> [extra input lines...]
   local d="$1"; local hf="$2"; local nf="$3"; shift 3
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/metals.inp" "$WORK/$d/"
   sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/$d/input.inp"
   cp "$hf" "$WORK/$d/output/Hydro_ioniz_IC.txt"
   cp "$nf" "$WORK/$d/output/Ion_species_IC.txt"
   local line
   for line in "$@"; do printf '%s\n' "$line" >> "$WORK/$d/input.inp"; done
}

run_it() {   # run_it <dir> [env assignments...]
   local d="$1"; shift
   ( cd "$WORK/$d" && env OMP_NUM_THREADS=1 "$@" "$EXE" > run.log 2>&1 )
   RC_LAST=$?
   return 0
}

META_RE='^# (restart_schema|reservoir|species_columns|grid|constants|options|t_phys\[s\]|source) '

# ---- A: the evaluation of a certified state, from the state as written ----
# The block is stripped from A's own copy of the pair, so that A is a
# BLOCK-LESS load by construction: it is the row that checks the
# provenance-unknown mark, and the state the matrix leaves in the case
# directory carries a block as soon as that case has been re-run.
make_run A "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
        'Restart intent: stationary evaluate'
sed -i -E "/$META_RE/d;/^# option_change /d" \
   "$WORK/A/output/Hydro_ioniz_IC.txt" "$WORK/A/output/Ion_species_IC.txt"
run_it A
rcA=$RC_LAST
if [ $rcA -ne 0 ] && [ $rcA -ne 2 ]; then
   echo "FAIL restart_intent_stage_A measured=exit_$rcA reference=0 tol=0"
   tail -n 5 "$WORK/A/run.log"
   exit 1
fi

# ---- B: the same evaluation of A's output: the metadata round trip -------
make_run B "$WORK/A/output/Hydro_ioniz.txt" "$WORK/A/output/Ion_species.txt" \
        'Restart intent: stationary evaluate'
run_it B
rcB=$RC_LAST

# ---- C to H: one field of the block made to disagree ---------------------
# Every case starts from B's input pair, which carries the block, and changes
# ONE field in BOTH files (the two halves of a state must agree).
for c in reservoir grid constants options elements pair; do
   make_run "$c" "$WORK/A/output/Hydro_ioniz.txt" \
                 "$WORK/A/output/Ion_species.txt" \
            'Restart intent: stationary evaluate'
   case $c in
   reservoir) sed -i -E 's|^# reservoir He/H [0-9.E+-]+|# reservoir He/H 9.0000000000000002E-02|' \
                 "$WORK/$c/output/Hydro_ioniz_IC.txt" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   grid)      sed -i -E 's|^# grid N ([0-9]+)|# grid N 999|' \
                 "$WORK/$c/output/Hydro_ioniz_IC.txt" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   constants) sed -i -E 's|kB\[erg/K\] [0-9.E+-]+|kB[erg/K] 1.3806500000000000E-16|' \
                 "$WORK/$c/output/Hydro_ioniz_IC.txt" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   options)   sed -i -E 's|He23S=T|He23S=F|' \
                 "$WORK/$c/output/Hydro_ioniz_IC.txt" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   elements)  sed -i -E 's|^# reservoir He/H ([0-9.E+-]+).*|# reservoir He/H \1|' \
                 "$WORK/$c/output/Hydro_ioniz_IC.txt" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   pair)      sed -i -E "/$META_RE/d" \
                 "$WORK/$c/output/Ion_species_IC.txt" ;;
   esac
   run_it "$c"
   eval "rc_$c=$RC_LAST"
done

# ---- I to K: the three input errors -------------------------------------
mkdir -p "$WORK/e_noload/output"
cp "$CASE/input.inp" "$CASE/metals.inp" "$WORK/e_noload/"
printf 'Restart intent: relaxation\n' >> "$WORK/e_noload/input.inp"
run_it e_noload
rc_noload=$RC_LAST

make_run e_traj "$WORK/A/output/Hydro_ioniz.txt" \
                "$WORK/A/output/Ion_species.txt" \
         'Run mode: init' 'Restart intent: trajectory'
run_it e_traj
rc_traj=$RC_LAST

make_run e_nosolver "$WORK/A/output/Hydro_ioniz.txt" \
                    "$WORK/A/output/Ion_species.txt" \
         'Restart intent: stationary'
sed -i '/^Solver:/d' "$WORK/e_nosolver/input.inp"
run_it e_nosolver
rc_nosolver=$RC_LAST

# ---- L to N: the residual is a function of the state --------------------
# The atomic element-row reload of backup/regression/atomic_elem_newton: the
# same state measured by the evaluation, by the EXHALE_RESIDUAL diagnostic
# with no re-equilibration, and by the stationary solve's own first judged
# rows.
for d in x_eval x_resid x_solve; do
   mkdir -p "$WORK/$d/output"
   cp "$ACASE/input.inp" "$ACASE/metals.inp" "$WORK/$d/"
   sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/$d/input.inp"
   cp "$ACASE/IC/Hydro_ioniz_IC.txt" "$ACASE/IC/Ion_species_IC.txt" \
      "$WORK/$d/output/"
done
printf 'Restart intent: stationary evaluate\n' >> "$WORK/x_eval/input.inp"
printf 'Restart intent: stationary\n'          >> "$WORK/x_solve/input.inp"
run_it x_eval
run_it x_resid EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0
# The solve does not converge on this fixture, so it is bounded here. timeout
# signals the child it started, by PID; no process is selected by name.
(
   cd "$WORK/x_solve" || exit 1
   env OMP_NUM_THREADS=1 "$EXE" > run.log 2>&1 &
   pid=$!
   t0=$(date +%s)
   while kill -0 "$pid" 2>/dev/null; do
      grep -q 'judged rows' run.log 2>/dev/null && break
      [ $(( $(date +%s) - t0 )) -ge $SOLVE_SECONDS ] && break
      sleep 1
   done
   kill "$pid" 2>/dev/null
   wait "$pid" 2>/dev/null
)

python3 - "$WORK" "$rcA" "$rcB" "$rc_reservoir" "$rc_grid" "$rc_constants" \
          "$rc_options" "$rc_elements" "$rc_pair" "$rc_noload" "$rc_traj" \
          "$rc_nosolver" <<'PY'
import re
import sys

import numpy as np

work = sys.argv[1]
(rcA, rcB, rc_res, rc_grid, rc_const, rc_opt, rc_elem, rc_pair,
 rc_noload, rc_traj, rc_nosolver) = [int(x) for x in sys.argv[2:13]]

ulp_tol = 1.0e-14      # stored text fields, a few units in the last place
rho_tol = 1.0e-12      # the writer state's own mass closure
sweep_tol = 1.0e-8     # the cell chemistry solve's xtol = sqrt(eps)
n_fail = 0
META = ('restart_schema', 'reservoir', 'species_columns', 'grid',
        'constants', 'options', 't_phys[s]', 'source')


def verdict(name, measured, reference, tol):
    global n_fail
    ok = measured <= tol
    if not ok:
        n_fail += 1
    print('%s %s measured=%.6e reference=%.6e tol=%.2e'
          % ('PASS' if ok else 'FAIL', name, measured, reference, tol))


def read(path):
    # The physical rows are the ones the header names: the first and last Ng
    # rows are boundary data, written by Apply_BC from the interior, and the
    # rebuild of a loaded state re-derives them by construction.
    labels, rows, jphys = None, [], (0, 0)
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            body = s.lstrip('#').split()
            if body[:1] == ['columns']:
                labels = body[1:]
            m = re.search(r'physical cells are rows (\d+) to (\d+)', s)
            if m:
                jphys = (int(m.group(1)) - 1, int(m.group(2)))
            continue
        if s:
            rows.append([float(x) for x in s.split()])
    return labels, np.array(rows), jphys


def block(path):
    out = {}
    for line in open(path):
        s = line.strip()
        if not s.startswith('#'):
            continue
        body = s.lstrip('#').strip()
        tok = body.split()[:1]
        if tok and tok[0] in META:
            out[tok[0]] = body
    return out


def maxrel(a, b):
    d, s = np.abs(a - b), np.maximum(np.abs(a), np.abs(b))
    m = s > 0
    rel = np.zeros_like(d)
    rel[m] = d[m] / s[m]
    j = int(np.argmax(rel))
    return float(rel[j]), j + 1


def log(d):
    return open('%s/%s/run.log' % (work, d)).read()


# ---- the evaluation took no step ----------------------------------------
print('---- "Restart intent: stationary evaluate" on a certified state ----')
mA = log('A')
m = re.search(r'steps:\s+(\d+) accepted of (\d+) attempted', mA)
if m is None:
    print('FAIL stationary_evaluate_takes_no_step measured=no_counter '
          'reference=0_of_0 tol=0')
    n_fail += 1
else:
    print('  the run counter reports %s accepted of %s attempted steps'
          % (m.group(1), m.group(2)))
    verdict('stationary_evaluate_takes_no_step',
            float(int(m.group(1)) + int(m.group(2))), 0.0, 0.0)
m = re.search(r"sweep's own root: max \|d\(n_tot\+n_e\)\|/\(n_tot\+n_e\) "
              r"=\s*([0-9.E+-]+)", mA)
if m:
    print('  DIAGNOSTIC the loaded composition against the sweep\'s own '
          'root: %s' % m.group(1))

# ---- the state came back --------------------------------------------------
worst = {}
ghost = (0.0, 'none', 0)
for f in ('Hydro_ioniz', 'Ion_species'):
    la, A, jp = read('%s/A/output/%s_IC.txt' % (work, f))
    lb, B, jq = read('%s/A/output/%s.txt' % (work, f))
    if A.shape != B.shape or la != lb or jp != jq:
        print('FAIL stationary_evaluate_schema measured=changed '
              'reference=same tol=0')
        n_fail += 1
        continue
    j0, j1 = jp
    for i, name in enumerate(la):
        d, j = maxrel(A[j0:j1, i], B[j0:j1, i])
        key = ('stored' if name in ('r[Rp]', 'v[cm/s]', 'p[cgs]')
               else 'rho' if name == 'rho[mH/cm3]' else 'rebuilt')
        if d > worst.get(key, (0.0, '', 0))[0]:
            worst[key] = (d, name, j + j0)
        dg, jg = maxrel(np.concatenate((A[:j0, i], A[j1:, i])),
                        np.concatenate((B[:j0, i], B[j1:, i])))
        if dg > ghost[0]:
            ghost = (dg, name, jg)
print('  DIAGNOSTIC the boundary rows, which the rebuild re-derives from '
      'the interior: worst %s %.3e' % (ghost[1], ghost[0]))
for key, tol, row in (('stored', ulp_tol,
                       'stationary_evaluate_stored_fields'),
                      ('rho', rho_tol, 'stationary_evaluate_mass_density'),
                      ('rebuilt', sweep_tol,
                       'stationary_evaluate_rebuilt_within_sweep_budget')):
    d, name, j = worst.get(key, (0.0, 'none', 0))
    print('  worst %-8s %-16s %.3e at row %d' % (key, name, d, j))
    verdict(row, d, 0.0, tol)

# ---- the certification travelled with it --------------------------------
def coupling(path):
    for line in open(path):
        if line.startswith('# coupling:'):
            return line.strip()
    return ''


cw = coupling('%s/A/output/Hydro_ioniz_IC.txt' % work)
cr = coupling('%s/A/output/Hydro_ioniz.txt' % work)
print('  writer:      %s' % cw)
print('  re-evaluated %s' % cr)
same = ('certified=T' in cw) == ('certified=T' in cr)
verdict('stationary_evaluate_certification_unchanged',
        0.0 if same else 1.0, 0.0, 0.0)


def coupling_field(line, key):
    for tok in line.split():
        if tok.startswith(key + '='):
            return tok.split('=', 1)[1]
    return ''


# THE STEP THE COUPLING WAS ARMED AT IS THE STATE'S, NOT THE READING RUN'S.
# A stationary evaluation takes no step and writes the loaded state back
# unchanged, so the state it writes was armed where the file says it was; a
# 0 here would give the state a provenance no run produced.
step_w = coupling_field(cw, 'sec_ion_step')
step_r = coupling_field(cr, 'sec_ion_step')
print('  armed step: writer %s, re-evaluated %s' % (step_w, step_r))
verdict('stationary_evaluate_carries_the_armed_step',
        0.0 if (step_w != '' and step_r == step_w) else 1.0, 0.0, 0.0)
verdict('stationary_evaluate_exit_status', float(rcA), 0.0, 0.0)

# ---- the metadata round trip --------------------------------------------
print('---- the metadata block: written, loaded, written ----')
bad = 0
for f in ('Hydro_ioniz', 'Ion_species'):
    b_in = block('%s/B/output/%s_IC.txt' % (work, f))
    b_out = block('%s/B/output/%s.txt' % (work, f))
    if sorted(b_in) != sorted(META) or sorted(b_out) != sorted(META):
        print('  %s: fields in  %s' % (f, sorted(b_in)))
        print('  %s: fields out %s' % (f, sorted(b_out)))
        bad += 1
        continue
    for k in META:
        if b_in[k] != b_out[k]:
            print('  %s field %s differs:' % (f, k))
            print('    in  %s' % b_in[k])
            print('    out %s' % b_out[k])
            bad += 1
verdict('restart_metadata_round_trip_identical', float(bad), 0.0, 0.0)

# ---- a pair with no block is loaded and marked --------------------------
print('---- a state pair written before the block existed ----')
marked_log = 'provenance unknown' in mA
src = block('%s/A/output/Hydro_ioniz.txt' % work).get('source', '')
print('  the load said "provenance unknown": %s' % marked_log)
print('  the file this run wrote carries: %s' % src)
verdict('restart_legacy_marked_in_the_log',
        0.0 if marked_log else 1.0, 0.0, 0.0)
verdict('restart_legacy_marked_in_the_written_file',
        0.0 if 'restart_input=provenance_unknown' in src else 1.0, 0.0, 0.0)

# ---- each field's mismatch is refused by name ---------------------------
print('---- a block that disagrees with the run ----')
for case, rc, field in (('reservoir', rc_res, 'reservoir'),
                        ('grid', rc_grid, 'grid'),
                        ('constants', rc_const, 'constants'),
                        ('options', rc_opt, 'options'),
                        ('elements', rc_elem, 'reservoir')):
    txt = log(case)
    named = ('metadata field "%s"' % field) in txt
    for line in txt.splitlines():
        if 'load_IC) ERROR' in line:
            print('  %-10s %s' % (case, line.strip()))
    verdict('restart_metadata_refuses_%s' % case,
            0.0 if (rc != 0 and named) else 1.0, 0.0, 0.0)
txt = log('pair')
paired = 'carries a restart metadata block and the other does not' in txt
verdict('restart_metadata_refuses_a_mismatched_pair',
        0.0 if (rc_pair != 0 and paired) else 1.0, 0.0, 0.0)

# ---- the three input errors ---------------------------------------------
print('---- the intent against the other keys ----')
for case, rc, needle, row in (
        ('e_noload', rc_noload, '"Load IC?" is False',
         'restart_intent_refused_without_a_loaded_state'),
        ('e_traj', rc_traj, '"Restart intent: trajectory" with',
         'restart_intent_trajectory_refused_in_init_mode'),
        ('e_nosolver', rc_nosolver, '"Restart intent: stationary" needs',
         'restart_intent_stationary_needs_a_stationary_solver')):
    txt = log(case)
    for line in txt.splitlines():
        if 'input_read) ERROR' in line:
            print('  %-11s %s' % (case, line.strip()))
    verdict(row, 0.0 if (rc != 0 and needle in txt) else 1.0, 0.0, 0.0)

# ---- the residual is a function of the state ----------------------------
print('---- the same state measured by three routes ----')
RN = re.compile(r'\|\|R\|\| = max_k :\s*([0-9.E+-]+)')
JR = re.compile(r'judged rows:.*\|\|R\|\|\s*([0-9.E+-]+)')


def first(rx, text):
    m = rx.search(text)
    return m.group(1) if m else ''


r_eval = first(RN, log('x_eval'))
r_diag = first(RN, log('x_resid'))
r_jfnk = first(JR, log('x_solve'))
print('  the evaluation of the state as loaded   ||R|| = %s' % r_eval)
print('  EXHALE_RESIDUAL=1 EXHALE_RELOAD_EQ=0    ||R|| = %s' % r_diag)
print('  the first rows the stationary solve judged ||R|| = %s' % r_jfnk)
ok = bool(r_eval) and r_eval == r_diag
verdict('stationary_residual_equals_the_diagnostic',
        0.0 if ok else 1.0, 0.0, 0.0)
# The solve prints three digits where the diagnostic prints five, so the
# comparison is on the digits the solve prints.
ok = bool(r_jfnk) and bool(r_eval) and \
    abs(float(r_jfnk) - float(r_eval)) <= 5.0e-4 * abs(float(r_eval))
verdict('stationary_solve_first_rows_are_the_loaded_state',
        0.0 if ok else 1.0, 0.0, 0.0)

sys.exit(1 if n_fail else 0)
PY
exit $?
