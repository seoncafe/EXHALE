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
#   7. The certification pair of the '# coupling:' header, (certified=,
#      cert_reason=), is metadata OF THE IMPORTED STATE: a reload that writes
#      the state back without measuring it carries the pair through
#      unchanged, whatever the reason says, while a run that measures the
#      state writes the verdict it measured and no earlier token.  A reason
#      longer than the field, a pair whose two halves disagree and a key
#      stated twice refuse the load; a field one half omits is taken from the
#      other with a note, and a token this version does not know is kept as
#      provenance text without certifying anything.
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
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" restart_intent_and_metadata
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
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
   printf '%s\n' "$RC_LAST" > "$WORK/$d/rc"
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

# ---- P: the certification pair the two halves of a state state ----------
# The pair is metadata OF THE IMPORTED STATE. Each row below hands the same
# written state back with its '# coupling:' header edited and reloads it with
# EXHALE_DUMP_IC=1, which writes the state as it was loaded, before the first
# equilibrium sweep, and stops: what the run writes is therefore what the
# loader restored and nothing else. The rows that must be REFUSED never reach
# that write, and their exit status and message are what is read.
CERT_H="$CASE/IC/Hydro_ioniz_IC.txt"
CERT_I="$CASE/IC/Ion_species_IC.txt"
# The destination is 32 characters (state_certification_reason, utilities.f90),
# so one token of exactly that length must be carried and one of 33 refused.
LIMIT32='reason_at_the_field_length_32chr'
OVER33='reason_one_character_over_the_32c'

cert_row() {   # cert_row <dir> <sed program for Hydro> [<sed program for Ion>]
   local d="$1"; local sh="$2"; local si="${3:-}"
   make_run "$d" "$CERT_H" "$CERT_I"
   if [ -n "$sh" ]; then sed -i -E "$sh" "$WORK/$d/output/Hydro_ioniz_IC.txt"; fi
   if [ -n "$si" ]; then sed -i -E "$si" "$WORK/$d/output/Ion_species_IC.txt"; fi
   run_it "$d" EXHALE_DUMP_IC=1
}

# The fixture states certified=T with no token, which is the second row as it
# stands; every other row edits that field.
cert_row p_true_token  's/certified=T/certified=T cert_reason=certified_in_wind/'
cert_row p_true_plain  ''
cert_row p_false_claim 's/certified=T/certified=F cert_reason=no_stationary_claim/'
cert_row p_false_entry 's/certified=T/certified=F cert_reason=failing_entries/'
cert_row p_false_plain 's/certified=T/certified=F/'
# The header is a set of key=value fields, so its order carries no meaning.
cert_row p_reordered \
   's|^# coupling: .*|# coupling: cert_reason=failing_entries certified=F mode=init recon=WENO3 sec_ion=T sec_ion_step=2429|'
cert_row p_limit32     "s/certified=T/certified=T cert_reason=$LIMIT32/"
cert_row p_over32      "s/certified=T/certified=T cert_reason=$OVER33/"
# A token this version does not know: provenance text, and the Boolean stays
# the file's own.
cert_row p_unknown     's/certified=T/certified=F cert_reason=foreign_token_x/'
# One key stated twice: a header that says two things about one state.
cert_row p_duplicate   's/certified=T/certified=T certified=T/'
# The two halves of one state, disagreeing and agreeing about the same field.
# Ion_species carries no '# coupling:' line of its own, so these rows write
# one into it after its first header line.
cert_row p_conflict    '' '1a # coupling: certified=F cert_reason=failing_entries'
cert_row p_legacy      's/certified=T/certified=T cert_reason=certified_in_wind/' \
                       '1a # coupling: certified=T'

# A mapped seed: 'map_state_to_grid.py' writes the state's own '# coupling:'
# line as an initialization seed and keeps what the SOURCE state was produced
# under on a '# mapped-from-coupling:' provenance line beside it. That line is
# not a statement about the state in this file, and the pair and the fields of
# the state's own line are the ones that count. The one written here differs
# from the state's line in both the pair and the armed step.
cert_row p_mapped_from \
   '/^# coupling:/a # mapped-from-coupling: sec_ion=T sec_ion_step=99 recon=WENO3 certified=F cert_reason=mapper_source_token mode=init'

# The verdict of a run that MEASURES the state is that run's own. A's output
# is handed back with a certified claim and a token on it, under a stellar EUV
# luminosity twice the one the state was solved at: the state is not
# stationary for that heating, and the evaluation has to answer with what it
# measured instead of carrying the claim through.
make_run p_reeval "$WORK/A/output/Hydro_ioniz.txt" \
                  "$WORK/A/output/Ion_species.txt" \
         'Restart intent: stationary evaluate'
sed -i -E 's/certified=[TF]( cert_reason=[^ ]+)?/certified=T cert_reason=certified_in_wind/' \
   "$WORK/p_reeval/output/Hydro_ioniz_IC.txt"
sed -i -E 's|^Log10 of EUV luminosity \[erg/s\]:.*|Log10 of EUV luminosity [erg/s]: 30.72|' \
   "$WORK/p_reeval/input.inp"
run_it p_reeval

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
rc_block1=$?

python3 - "$WORK" <<'PY'
# The certification pair of the '# coupling:' header.
import os
import sys

work = sys.argv[1]
n_fail = 0


def verdict(name, ok, measured, reference):
    global n_fail
    if not ok:
        n_fail += 1
    print('%s %s measured=%s reference=%s tol=0'
          % ('PASS' if ok else 'FAIL', name, measured, reference))


def rc(d):
    return int(open('%s/%s/rc' % (work, d)).read().strip())


def log(d):
    return open('%s/%s/run.log' % (work, d)).read()


def coupling(path):
    if not os.path.exists(path):
        return ''
    for line in open(path):
        if line.startswith('# coupling:'):
            return line.strip()
    return ''


def pair(path):
    # (certified, cert_reason) as the line states them; '-' for a field the
    # line does not carry, which is not the same statement as an empty one.
    line = coupling(path)
    out = {'certified': '-', 'cert_reason': '-'}
    for tok in line.split():
        k, _, v = tok.partition('=')
        if k in out and _:
            out[k] = v
    return out['certified'], out['cert_reason']


def written(d):
    return pair('%s/%s/output/Hydro_ioniz.txt' % (work, d))


def imported(d):
    return pair('%s/%s/output/Hydro_ioniz_IC.txt' % (work, d))


print('---- the pair a reload carries through, unmeasured ----')
for d, name in (('p_true_token', 'cert_pair_true_with_token'),
                ('p_true_plain', 'cert_pair_true_without_token'),
                ('p_false_claim', 'cert_pair_false_no_stationary_claim'),
                ('p_false_entry', 'cert_pair_false_failing_entries'),
                ('p_false_plain', 'cert_pair_false_without_reason'),
                ('p_reordered', 'cert_pair_independent_of_field_order'),
                ('p_limit32', 'cert_pair_reason_at_the_field_length'),
                ('p_unknown', 'cert_pair_unknown_token_kept')):
    imp, wrt = imported(d), written(d)
    print('  %-14s imported %-40s written %s'
          % (d, '%s / %s' % imp, '%s / %s' % wrt))
    verdict(name, rc(d) == 0 and wrt == imp,
            'certified=%s cert_reason=%s' % wrt,
            'certified=%s cert_reason=%s' % imp)

# The Boolean is the 'certified=' field alone: a token this version does not
# know is provenance text and certifies nothing by itself.
c_unknown, r_unknown = written('p_unknown')
verdict('cert_unknown_token_does_not_certify', c_unknown == 'F',
        'certified=%s cert_reason=%s' % (c_unknown, r_unknown), 'certified=F')

print('---- the headers that are refused ----')
for d, name, needle in (
        ('p_over32', 'cert_reason_longer_than_the_field_refused',
         'longer than the field'),
        ('p_duplicate', 'cert_duplicate_key_refused', 'stated twice'),
        ('p_conflict', 'cert_conflicting_pair_refused',
         'different stationary claims')):
    txt = log(d)
    for line in txt.splitlines():
        if 'load_IC) ERROR' in line:
            print('  %-12s %s' % (d, line.strip()))
    verdict(name, rc(d) != 0 and needle in txt,
            'exit=%d named=%s' % (rc(d), needle in txt), 'exit!=0 named=True')

print('---- a field one half of the state omits ----')
imp, wrt = imported('p_legacy'), written('p_legacy')
noted = 'only one of the two restart files states "cert_reason"' \
        in log('p_legacy')
print('  imported %s / %s, written %s / %s, note printed %s'
      % (imp[0], imp[1], wrt[0], wrt[1], noted))
verdict('cert_omitted_field_taken_from_the_other_half',
        rc('p_legacy') == 0 and wrt == imp and noted,
        'certified=%s cert_reason=%s note=%s' % (wrt[0], wrt[1], noted),
        'certified=%s cert_reason=%s note=True' % imp)

print('---- a mapped seed, whose provenance line is not its header ----')
imp, wrt = imported('p_mapped_from'), written('p_mapped_from')
step_w = ''
for tok in coupling('%s/p_mapped_from/output/Hydro_ioniz.txt' % work).split():
    if tok.startswith('sec_ion_step='):
        step_w = tok.split('=', 1)[1]
print('  the state\'s own line   certified=%s cert_reason=%s' % imp)
print('  the mapped-from line   certified=F cert_reason=mapper_source_token'
      ' sec_ion_step=99')
print('  written                certified=%s cert_reason=%s sec_ion_step=%s'
      % (wrt[0], wrt[1], step_w))
verdict('mapped_from_provenance_is_not_a_coupling_header',
        rc('p_mapped_from') == 0 and wrt == imp and step_w == '2429',
        'exit=%d certified=%s cert_reason=%s sec_ion_step=%s'
        % (rc('p_mapped_from'), wrt[0], wrt[1], step_w),
        'exit=0 certified=%s cert_reason=%s sec_ion_step=2429' % imp)

print('---- the verdict of a run that measures the state ----')
imp, wrt = imported('p_reeval'), written('p_reeval')
print('  the file claimed   certified=%s cert_reason=%s' % imp)
print('  the evaluation put certified=%s cert_reason=%s' % wrt)
# Exit 2 is a claim the certification refused; the state is written in full.
verdict('reevaluation_writes_its_own_pair',
        rc('p_reeval') in (0, 2) and wrt != imp
        and wrt[1] in ('failing_entries', 'no_stationary_claim', '-',
                       'certified_in_wind'),
        'certified=%s cert_reason=%s' % wrt,
        'the pair this run measured, not certified=%s cert_reason=%s' % imp)

sys.exit(1 if n_fail else 0)
PY
rc_block2=$?

if [ $rc_block1 -ne 0 ] || [ $rc_block2 -ne 0 ]; then
   exit 1
fi
exit 0
