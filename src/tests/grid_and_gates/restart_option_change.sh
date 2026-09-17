#!/bin/bash
# A restart that changes named physics options
# (decision 21 of docs/To_be_determined_by_user_20260906.md, option a;
#  PLAN_20260909_rev1 item N10b).
#
# WHAT IS UNDER TEST
#   1. A rung of an option ladder loads: a state written with one option set is
#      restarted under another, and the load is accepted for exactly the
#      tokens "Restart option change:" names.
#   2. The same difference with nothing named still refuses the load, and the
#      refusal states the token that differs.
#   3. The key is checked against the vocabulary of the state file's
#      '# options' line where the input is read: an unknown token stops the
#      run, a token that decides how many unknowns the state has stops the
#      run with that reason, and the key without a loaded state stops the
#      run.
#   4. A named token that does NOT differ is reported and nothing else: the
#      key is a permission and not an instruction.
#   5. The change is written into the state the run produces as one
#      '# option_change' line, and a rung inherits the lines of the rung it
#      was restarted from, so a ladder can be read back from its last state.
#   6. The ROUTE token carrier_newton is not an equation-set token. An
#      alternation state (the transported balances relaxed at a held wind)
#      loads under "Coupled carrier solve: True" and under "On stall" with
#      nothing named, because the rows, the columns and the tolerances are
#      the same and only the algorithm differs; the difference is reported,
#      written as a '# route_change' line, and the state the run writes
#      records the route it came out of.
#   7. EXHALE_setup.out states what was known about the loaded state's
#      configuration and what this run was allowed to change about it. A pair
#      that carries no block at all and a pair whose state DESCENDS from one
#      are both provenance-unknown and are distinguished there: rung 0 loads
#      the first, and every rung after it the second.
#
# THE FIXTURE is backup/regression/atomic_elem_newton/IC, a pinned pair of
# state files written before the configuration block existed, so the first
# run here is also the provenance-unknown load.  Every run takes
# "Restart intent: stationary evaluate", which measures the loaded state and
# writes it back without a step: what is under test is the load, and no
# marching or solving is needed to reach it.  Nothing in backup/ is written.
#
# EXPECTED AT THE ENTRY TEXT OF N10b (advisor binary EXHALE_lwv.x): RED,
# all eleven rows, MEASURED.  There the key is an unrecognized input line (a
# warning, no refusal), so the ladder rung and its unnamed twin are both
# refused by an exact comparison of the whole field, which states the field
# and not the token that differs; no '# option_change' line is written and no
# setup report line exists.
#
# Usage: restart_option_change.sh   (EXHALE_EXE selects the binary)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}"
WORK="$OUT/restart_option_change"
CASE="$ROOT/backup/regression/atomic_elem_newton"

n_fail=0
verdict() {   # verdict <name> <measured> <reference>
   if [ "$2" = "$3" ]; then
      echo "PASS $1 measured=$2 reference=$3 tol=0"
   else
      echo "FAIL $1 measured=$2 reference=$3 tol=0"
      n_fail=$((n_fail+1))
   fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL restart_option_change measured=no_binary reference=$EXE tol=0"
   exit 1
fi
for f in "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt"; do
   if [ ! -s "$f" ]; then
      echo "FAIL restart_option_change measured=no_fixture reference=$f tol=0"
      exit 1
   fi
done

rm -rf "$WORK"
mkdir -p "$WORK"

# make_run <dir> <hydro state> <ion state> [extra input lines...]
make_run() {
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

log() { cat "$WORK/$1/run.log"; }
n_change_lines() { grep -c '^# option_change ' "$1" 2>/dev/null || true; }

# ---- rung 0: the pinned pair, loaded and written back with the block -----
make_run rung0 "$CASE/IC/Hydro_ioniz_IC.txt" "$CASE/IC/Ion_species_IC.txt" \
         'Restart intent: stationary evaluate'
run_it rung0
rc0=$RC_LAST
if [ $rc0 -ne 0 ] && [ $rc0 -ne 2 ]; then
   echo "FAIL restart_option_change_rung0 measured=exit_$rc0 reference=0_or_2 tol=0"
   tail -n 8 "$WORK/rung0/run.log"
   exit 1
fi
H0="$WORK/rung0/output/Hydro_ioniz.txt"
I0="$WORK/rung0/output/Ion_species.txt"
if [ ! -s "$H0" ] || [ ! -s "$I0" ]; then
   echo "FAIL restart_option_change_rung0 measured=no_state reference=written tol=0"
   exit 1
fi
echo "  rung 0 wrote a state carrying the configuration block:"
grep '^# options ' "$H0" | sed 's/^/    /'

# ---- rung 1: one named token differs -> the load is accepted -------------
make_run rung1 "$H0" "$I0" 'Secondary_ionization: False' \
         'Restart option change: sec_ion' \
         'Restart intent: stationary evaluate'
run_it rung1
rc1=$RC_LAST

# ---- the same difference with nothing named -> refused by name -----------
make_run unnamed "$H0" "$I0" 'Secondary_ionization: False' \
         'Restart intent: stationary evaluate'
run_it unnamed
rc_unnamed=$RC_LAST

# ---- a named token that does not differ -> reported, not refused ---------
make_run inert "$H0" "$I0" 'Restart option change: cond' \
         'Restart intent: stationary evaluate'
run_it inert
rc_inert=$RC_LAST

# ---- the three input errors of the key -----------------------------------
make_run unknown "$H0" "$I0" 'Restart option change: LW' \
         'Restart intent: stationary evaluate'
run_it unknown
rc_unknown=$RC_LAST

make_run layout "$H0" "$I0" 'Restart option change: mol' \
         'Restart intent: stationary evaluate'
run_it layout
rc_layout=$RC_LAST

make_run gridtok "$H0" "$I0" 'Restart option change: N' \
         'Restart intent: stationary evaluate'
run_it gridtok
rc_gridtok=$RC_LAST

# The key with no loaded state: a COLD run of the fixture, so it is bounded
# to one marching step with no solver to follow. What is under test is
# refused before anything runs; the bound is there so that a binary which
# does not refuse it does not relax the fixture instead.
mkdir -p "$WORK/noload/output"
cp "$CASE/input.inp" "$CASE/metals.inp" "$WORK/noload/"
sed -i '/^Solver:/d' "$WORK/noload/input.inp"
printf 'Restart option change: sec_ion\n' >> "$WORK/noload/input.inp"
run_it noload EXHALE_MAXSTEPS=1
rc_noload=$RC_LAST

# ---- rung 2: the history is inherited ------------------------------------
# Restarted from rung 1's state, which already carries one option_change
# line, with a second option allowed to differ.
if [ -s "$WORK/rung1/output/Hydro_ioniz.txt" ]; then
   make_run rung2 "$WORK/rung1/output/Hydro_ioniz.txt" \
                  "$WORK/rung1/output/Ion_species.txt" \
            'Secondary_ionization: False' 'Conduction: True' \
            'Restart option change: cond' \
            'Restart intent: stationary evaluate'
   run_it rung2
   rc2=$RC_LAST
else
   rc2=-1
fi

# ---- the verdicts --------------------------------------------------------
echo "---- a rung that changes a named option ----"
grep -E '(load_IC\) (ERROR|the restart CHANGES))' "$WORK/rung1/run.log" \
   | sed 's/^/    /'
ok=no
if { [ $rc1 -eq 0 ] || [ $rc1 -eq 2 ]; } && \
   grep -q 'the restart CHANGES the physics options' "$WORK/rung1/run.log" && \
   ! grep -q 'load_IC) ERROR' "$WORK/rung1/run.log"; then ok=yes; fi
verdict option_change_named_token_loads "$ok" yes

nl=$(n_change_lines "$WORK/rung1/output/Hydro_ioniz.txt")
ni=$(n_change_lines "$WORK/rung1/output/Ion_species.txt")
echo "    the state rung 1 wrote carries:"
grep '^# option_change ' "$WORK/rung1/output/Hydro_ioniz.txt" \
   | sed 's/^/      /' || true
ok=no
if [ "$nl" = "1" ] && [ "$ni" = "1" ] && \
   grep -q '^# option_change sec_ion=T -> sec_ion=F at restart of ' \
        "$WORK/rung1/output/Hydro_ioniz.txt"; then ok=yes; fi
verdict option_change_line_written_in_both_halves "$ok" yes

echo "---- the same difference with nothing named ----"
grep 'load_IC) ERROR' "$WORK/unnamed/run.log" | sed 's/^/    /' || true
ok=no
if [ $rc_unnamed -ne 0 ] && \
   grep -q 'metadata field "options"' "$WORK/unnamed/run.log" && \
   grep -q 'sec_ion: T -> F' "$WORK/unnamed/run.log"; then ok=yes; fi
verdict option_change_unnamed_difference_refused_by_name "$ok" yes

echo "---- a named token that does not differ ----"
grep 'do not differ between' "$WORK/inert/run.log" | sed 's/^/    /' || true
ok=no
if { [ $rc_inert -eq 0 ] || [ $rc_inert -eq 2 ]; } && \
   grep -q 'names 1 token(s) that do not differ' "$WORK/inert/run.log" && \
   [ "$(n_change_lines "$WORK/inert/output/Hydro_ioniz.txt")" = "0" ]; then
   ok=yes
fi
verdict option_change_named_but_equal_is_reported_only "$ok" yes

echo "---- the key against its own vocabulary ----"
for c in unknown layout gridtok noload; do
   grep 'input_read) ERROR' "$WORK/$c/run.log" | sed "s/^/    $c /" || true
done
ok=no
if [ $rc_unknown -ne 0 ] && \
   grep -q 'unknown token "LW"' "$WORK/unknown/run.log"; then ok=yes; fi
verdict option_change_unknown_token_is_an_input_error "$ok" yes
ok=no
if [ $rc_layout -ne 0 ] && \
   grep -q 'decides how many unknowns' "$WORK/layout/run.log"; then ok=yes; fi
verdict option_change_layout_token_is_an_input_error "$ok" yes
ok=no
if [ $rc_gridtok -ne 0 ] && \
   grep -q 'not a physics option but' "$WORK/gridtok/run.log"; then ok=yes; fi
verdict option_change_grid_field_is_an_input_error "$ok" yes
ok=no
if [ $rc_noload -ne 0 ] && \
   grep -q '"Load IC?" is False' "$WORK/noload/run.log"; then ok=yes; fi
verdict option_change_without_a_loaded_state_is_an_input_error "$ok" yes

echo "---- the history a ladder can be read back from ----"
if [ $rc2 -eq -1 ]; then
   echo "    rung 1 wrote no state, so rung 2 could not be run"
   n2=0
else
   grep '^# option_change ' "$WORK/rung2/output/Hydro_ioniz.txt" \
      | sed 's/^/    /' || true
   n2=$(n_change_lines "$WORK/rung2/output/Hydro_ioniz.txt")
fi
ok=no
if [ "$n2" = "2" ] && \
   grep -q '^# option_change sec_ion=T -> sec_ion=F ' \
        "$WORK/rung2/output/Hydro_ioniz.txt" && \
   grep -q '^# option_change cond=F -> cond=T ' \
        "$WORK/rung2/output/Hydro_ioniz.txt"; then ok=yes; fi
verdict option_change_history_inherited_by_the_next_rung "$ok" yes

# ---- the 21st token: the well-balanced pressure/gravity pair ----------
# "Well balanced:" (K46) changes the discretization of the pressure/gravity
# pair and not how many unknowns the state has, so it is a token of the
# options field that a restart may be allowed to change by name.
make_run wbnamed "$H0" "$I0" 'Well balanced: True' \
         'Restart option change: wellbal' \
         'Restart intent: stationary evaluate'
run_it wbnamed
rc_wbnamed=$RC_LAST

make_run wbunnamed "$H0" "$I0" 'Well balanced: True' \
         'Restart intent: stationary evaluate'
run_it wbunnamed
rc_wbunnamed=$RC_LAST

echo "---- the well-balanced token ----"
grep -E 'wellbal' "$WORK/wbnamed/run.log" | sed 's/^/    named /' || true
grep -E 'wellbal' "$WORK/wbunnamed/run.log" | sed 's/^/    unnamed /' || true
ok=no
if { [ $rc_wbnamed -eq 0 ] || [ $rc_wbnamed -eq 2 ]; } && \
   grep -q 'the restart CHANGES the physics options' \
        "$WORK/wbnamed/run.log" && \
   ! grep -q 'load_IC) ERROR' "$WORK/wbnamed/run.log" && \
   grep -q '^# option_change wellbal=F -> wellbal=T at restart of ' \
        "$WORK/wbnamed/output/Hydro_ioniz.txt"; then ok=yes; fi
verdict option_change_well_balanced_token_loads "$ok" yes
ok=no
if [ $rc_wbunnamed -ne 0 ] && \
   grep -q 'wellbal: F -> T' "$WORK/wbunnamed/run.log"; then ok=yes; fi
verdict option_change_well_balanced_unnamed_refused_by_name "$ok" yes
ok=no
if grep -q 'Well balanced: on' "$WORK/wbnamed/EXHALE_setup.out" && \
   grep -q 'Well balanced: off' "$WORK/rung0/EXHALE_setup.out"; then
   ok=yes
fi
verdict well_balanced_key_is_stated_in_the_setup_report "$ok" yes
ok=no
if grep -q '^well_balanced  *T' "$WORK/wbnamed/EXHALE_resolved.out" && \
   grep -q '^well_balanced  *F' "$WORK/rung0/EXHALE_resolved.out"; then
   ok=yes
fi
verdict well_balanced_key_is_in_the_resolved_input "$ok" yes

echo "---- the setup report ----"
grep -E 'Restart (provenance|option change)' "$WORK/rung0/EXHALE_setup.out" \
   | sed 's/^/    rung0 /' || true
grep -E 'Restart (provenance|option change)|DID change a named option' \
     "$WORK/rung1/EXHALE_setup.out" | sed 's/^/    rung1 /' || true
ok=no
if grep -q 'Restart provenance: UNKNOWN\. The state files carry no' \
        "$WORK/rung0/EXHALE_setup.out" && \
   grep -q 'Restart provenance: UNKNOWN by descent' \
        "$WORK/rung1/EXHALE_setup.out"; then
   ok=yes
fi
verdict option_change_setup_report_states_unknown_provenance "$ok" yes
ok=no
if grep -q 'Restart option change: allowed for sec_ion' \
        "$WORK/rung1/EXHALE_setup.out" && \
   grep -q 'DID change a named option' "$WORK/rung1/EXHALE_setup.out"; then
   ok=yes
fi
verdict option_change_setup_report_states_the_change "$ok" yes

# ---- the route token: the same equations, another algorithm ----------
# THE FIXTURE here is backup/regression/carrier_model_a_newton, a molecular
# case with the H2 carrier transported and "Coupled carrier solve: False",
# which is what makes an ALTERNATION state: the carrier balance is relaxed
# at a held wind instead of standing in the Newton unknown vector. Every
# run takes "Restart intent: stationary evaluate", so what is under test is
# the load and no solve is needed to reach it.
CCASE="$ROOT/backup/regression/carrier_model_a_newton"

cmake_run() {   # cmake_run <dir> <hydro state> <ion state> [extra input lines]
   local d="$1"; local hf="$2"; local nf="$3"; shift 3
   mkdir -p "$WORK/$d/output"
   cp "$CCASE/input.inp" "$CCASE/base.inp" "$WORK/$d/"
   cp "$hf" "$WORK/$d/output/Hydro_ioniz_IC.txt"
   cp "$nf" "$WORK/$d/output/Ion_species_IC.txt"
   sed -i -e '/^Coupled carrier solve:/d' -e '/^Restart intent:/d' \
          "$WORK/$d/input.inp"
   local line
   for line in "$@"; do printf '%s\n' "$line" >> "$WORK/$d/input.inp"; done
}

route_ok=yes
for f in "$CCASE/IC/Hydro_ioniz_IC.txt" "$CCASE/IC/Ion_species_IC.txt"; do
   [ -s "$f" ] || route_ok=no
done

if [ "$route_ok" = no ]; then
   echo "FAIL route_change_fixture measured=no_fixture reference=$CCASE/IC tol=0"
   n_fail=$((n_fail+1))
else
   # rung A: the alternation writes a state, which carries carrier_newton=F
   cmake_run routeA "$CCASE/IC/Hydro_ioniz_IC.txt" "$CCASE/IC/Ion_species_IC.txt" \
            'Coupled carrier solve: False' 'Restart intent: stationary evaluate'
   run_it routeA
   rc_rA=$RC_LAST
   HA="$WORK/routeA/output/Hydro_ioniz.txt"
   IA="$WORK/routeA/output/Ion_species.txt"

   # rung B: the same state under the block, with nothing named
   if [ -s "$HA" ]; then
      cmake_run routeB "$HA" "$IA" \
               'Coupled carrier solve: True' 'Restart intent: stationary evaluate'
      run_it routeB
      rc_rB=$RC_LAST
      cmake_run routeC "$HA" "$IA" \
               'Coupled carrier solve: On stall' \
               'Restart intent: stationary evaluate'
      run_it routeC
      rc_rC=$RC_LAST
   else
      rc_rB=-1; rc_rC=-1
   fi

   # rung D: back to the alternation from the block's state, so the history
   # carries both route lines
   HB="$WORK/routeB/output/Hydro_ioniz.txt"
   if [ -s "$HB" ]; then
      cmake_run routeD "$HB" "$WORK/routeB/output/Ion_species.txt" \
               'Coupled carrier solve: False' 'Restart intent: stationary evaluate'
      run_it routeD
      rc_rD=$RC_LAST
   else
      rc_rD=-1
   fi

   echo "---- the route token, with nothing named ----"
   grep -E 'changes the ROUTE|load_IC\) ERROR' "$WORK/routeB/run.log" 2>/dev/null \
      | sed 's/^/    /' || true
   ok=no
   if { [ $rc_rB -eq 0 ] || [ $rc_rB -eq 2 ]; } && \
      grep -q 'changes the ROUTE and not the equations' \
           "$WORK/routeB/run.log" && \
      ! grep -q 'load_IC) ERROR' "$WORK/routeB/run.log"; then ok=yes; fi
   verdict route_change_block_loads_an_alternation_state "$ok" yes

   ok=no
   if grep -q '^# route_change carrier_newton=F -> carrier_newton=T at restart of ' \
           "$WORK/routeB/output/Hydro_ioniz.txt" 2>/dev/null && \
      grep -q '^# route_change carrier_newton=F -> carrier_newton=T at restart of ' \
           "$WORK/routeB/output/Ion_species.txt" 2>/dev/null; then ok=yes; fi
   verdict route_change_line_written_in_both_halves "$ok" yes

   ok=no
   if grep -m1 '^# options ' "$WORK/routeB/output/Hydro_ioniz.txt" 2>/dev/null \
        | grep -q ' carrier_newton=T' && \
      grep -m1 '^# options ' "$HA" 2>/dev/null | grep -q ' carrier_newton=F'; then
      ok=yes
   fi
   verdict route_change_written_state_records_its_own_route "$ok" yes

   echo "---- the route token under On stall ----"
   ok=no
   if { [ $rc_rC -eq 0 ] || [ $rc_rC -eq 2 ]; } && \
      ! grep -q 'load_IC) ERROR' "$WORK/routeC/run.log" && \
      grep -m1 '^# options ' "$WORK/routeC/output/Hydro_ioniz.txt" 2>/dev/null \
        | grep -q ' carrier_newton=F'; then ok=yes; fi
   verdict route_change_on_stall_loads_an_alternation_state "$ok" yes

   echo "---- the route history a ladder can be read back from ----"
   grep '^# route_change ' "$WORK/routeD/output/Hydro_ioniz.txt" 2>/dev/null \
      | sed 's/^/    /' || true
   nr=$(grep -c '^# route_change ' "$WORK/routeD/output/Hydro_ioniz.txt" 2>/dev/null || true)
   ok=no
   if [ "$nr" = "2" ]; then ok=yes; fi
   verdict route_change_history_inherited_by_the_next_rung "$ok" yes
fi

exit $(( n_fail > 0 ? 1 : 0 ))
