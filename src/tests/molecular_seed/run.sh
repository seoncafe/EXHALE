#!/bin/bash
# The molecular seed built from a solved atomic state (docs/input_schema.md
# appendix D.3), tested on the binary.
#
# WHAT IS UNDER TEST. A molecular run carries four unknowns an atomic run does
# not, so the restart contract refuses an atomic state file for it. The seed
# conversion writes those unknowns instead: it reads the atomic pair of
# another directory, forms the molecular state that pair implies, writes it
# under the names a restart reads, and stops. The assertions are the first
# three tests of docs/PLAN_20260913_lhs_stationary.md section 5.4:
#
#   T-L7-1  conversion identity: a zero requested transfer leaves every
#           species column and the kept primitives where they were
#   T-L7-2   equation-of-state closure: the written composition, p, T and
#           energy describe one state under the production caloric EOS, both
#           at a zero transfer and at the handoff's own fraction
#   T-L7-3  loader acceptance: the written pair loads into the target
#           molecular run with no option permitted to differ, and is treated
#           as initialization
#
# and the two refusals that keep the conversion from being a general
# restart-layout exception.
#
# The fixture is the hot Uranus of backup/regression/mol_base_handoff. Its
# atomic counterpart is the same input with "Molecular chemistry: False" and
# the handoff's q_H2_base removed from base.inp, so that neither the particle
# count nor the ghost composition claims molecules the atomic network cannot
# hold (input_read refuses the other combination by name). The two
# configurations then differ in the molecular tokens alone, and their base
# particle count differs by the H nuclei the handoff binds -- which is a
# different base DENSITY at the same 1 microbar level, not a different grid,
# reservoir or constant set. Both runs are capped at a few steps: what is
# under test is the conversion, not the wind.
#
# Every test prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion; lines beginning with two spaces are context. Exit status
# is nonzero if any assertion fails.
#
# Usage: src/tests/molecular_seed/run.sh
#   EXHALE_EXE       binary to test (default $ROOT/EXHALE.x)
#   EXHALE_TEST_OUT  working directory (default $ROOT/build/tests/molecular_seed)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/molecular_seed}"
REG="$ROOT/backup/regression/mol_base_handoff"
FAIL=0

chk_le () { # name measured tol
   if awk -v m="$2" -v t="$3" 'BEGIN{exit !(m+0 <= t+0)}'; then
      echo "PASS $1 measured=$2 reference=<=$3 tol=$3"
   else
      echo "FAIL $1 measured=$2 reference=<=$3 tol=$3"; FAIL=1
   fi
}
chk_eq () { # name measured reference
   if [ "$2" = "$3" ]; then echo "PASS $1 measured=$2 reference=$3 tol=0"
   else echo "FAIL $1 measured=$2 reference=$3 tol=0"; FAIL=1; fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL molecular_seed_binary measured=no_binary reference=$EXE tol=0"; exit 1
fi
if [ ! -f "$REG/input.inp" ]; then
   echo "FAIL molecular_seed_fixture measured=absent reference=$REG tol=0"; exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK"

# -- the atomic source: the same planet, base and grid, without the molecules --
A="$WORK/atomic"; mkdir -p "$A/output"
grep -v '^q_H2_base' "$REG/base.inp" > "$A/base.inp"
sed -e 's/^Molecular chemistry:.*/Molecular chemistry: False/' "$REG/input.inp" > "$A/input.inp"
( cd "$A" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=20 "$EXE" > run.log 2>&1 )
if [ ! -f "$A/output/Hydro_ioniz.txt" ]; then
   echo "FAIL molecular_seed_atomic_source measured=no_state reference=written tol=0"
   tail -n 15 "$A/run.log"; exit 1
fi
for f in Hydro_ioniz Ion_species; do
   \cp -f "$A/output/$f.txt" "$A/output/${f}_IC.txt"
done

mkcase () { # dir
   local d="$WORK/$1"
   rm -rf "$d"; mkdir -p "$d/output"
   cp "$REG/base.inp" "$d/"
   sed -e 's/^Load IC?.*/Load IC? True/' "$REG/input.inp" > "$d/input.inp"
}

# -------------------------------------------------------------------------- #
# T-L7-1 and the zero-transfer half of T-L7-2: nothing requested, nothing
# moved, and the state still closes under the equation of state.
# -------------------------------------------------------------------------- #
mkcase identity
( cd "$WORK/identity" && OMP_NUM_THREADS=1 EXHALE_MOLECULAR_SEED="$A/output" \
     EXHALE_MOLECULAR_SEED_X2=0 "$EXE" > run.log 2>&1 )
rc=$?
chk_eq molecular_seed_identity_exit "$rc" 0
if [ ! -f "$WORK/identity/output/Ion_species_IC.txt" ]; then
   echo "FAIL molecular_seed_identity_written measured=absent reference=present tol=0"
   tail -n 20 "$WORK/identity/run.log"; FAIL=1
else
   # Every column the two files share, over the PHYSICAL cells. The ghost rows
   # are the target run's own boundary, re-applied by the conversion, and are
   # a statement of the boundary rather than of the transfer.
   D=$(python3 "$HERE/compare_columns.py" "$A/output" "$WORK/identity/output")
   chk_le molecular_seed_identity_columns "$(echo "$D" | sed -n 1p)" 1.0e-12
   chk_le molecular_seed_identity_molecular_zero "$(echo "$D" | sed -n 2p)" 0.0
fi
E=$(sed -n 's/.*p round trip): *//p' "$WORK/identity/run.log" | tail -n 1 | tr -d ' ')
chk_le molecular_seed_identity_eos_closure "${E:-1}" 1.0e-12

# -------------------------------------------------------------------------- #
# T-L7-2 at the handoff's own fraction: the transfer conserves the nuclei and
# the mass, and the written state closes under the caloric equation of state.
# -------------------------------------------------------------------------- #
mkcase handoff
( cd "$WORK/handoff" && OMP_NUM_THREADS=1 EXHALE_MOLECULAR_SEED="$A/output" \
     "$EXE" > run.log 2>&1 )
rc=$?
chk_eq molecular_seed_handoff_exit "$rc" 0
C=$(sed -n 's/.*nuclei, mass: *//p' "$WORK/handoff/run.log" | tail -n 1)
chk_le molecular_seed_handoff_H_nuclei  "$(echo "$C" | awk '{print $1}')" 1.0e-12
chk_le molecular_seed_handoff_He_nuclei "$(echo "$C" | awk '{print $2}')" 1.0e-12
chk_le molecular_seed_handoff_mass      "$(echo "$C" | awk '{print $3}')" 1.0e-12
E=$(sed -n 's/.*p round trip): *//p' "$WORK/handoff/run.log" | tail -n 1 | tr -d ' ')
chk_le molecular_seed_handoff_eos_closure "${E:-1}" 1.0e-12
X=$(sed -n 's/.*largest x2 actually transferred *//p' "$WORK/handoff/run.log" | tail -n 1 | tr -d ' ')
chk_le molecular_seed_handoff_x2_le_one "${X:-2}" 1.0
grep -q '^# molecular_partition:' "$WORK/handoff/output/Hydro_ioniz_IC.txt" \
   && a=present || a=absent
grep -q '^# molecular_partition:' "$WORK/handoff/output/Ion_species_IC.txt" \
   && b=present || b=absent
chk_eq molecular_seed_record_in_hydro   "$a" present
chk_eq molecular_seed_record_in_species "$b" present

# -------------------------------------------------------------------------- #
# T-L7-3 loader acceptance: the pair the conversion wrote loads into the
# target molecular run with nothing permitted to differ, and is initialization.
# -------------------------------------------------------------------------- #
L="$WORK/reload"; rm -rf "$L"; mkdir -p "$L/output"
cp "$WORK/handoff/input.inp" "$WORK/handoff/base.inp" "$L/"
\cp -f "$WORK/handoff/output/Hydro_ioniz_IC.txt" \
       "$WORK/handoff/output/Ion_species_IC.txt" "$L/output/"
( cd "$L" && OMP_NUM_THREADS=1 EXHALE_DUMP_IC=1 "$EXE" > run.log 2>&1 )
chk_eq molecular_seed_reload_exit "$?" 0
n=$(sed -n 's/^# species_columns \([0-9]*\) .*/\1/p' \
    "$WORK/handoff/output/Ion_species_IC.txt" | head -n 1)
chk_eq molecular_seed_reload_columns "${n:-0}" 38
if grep -q 'certified=F cert_reason=molecular_seed mode=init' \
        "$WORK/handoff/output/Hydro_ioniz_IC.txt"; then s=init_uncertified; else s=other; fi
chk_eq molecular_seed_coupling_line "$s" init_uncertified

# -------------------------------------------------------------------------- #
# The two refusals. The conversion is not a general restart-layout exception:
# it is refused where it would have nothing to write, and where the partition
# would be the seed's own estimate rather than upstream information.
# -------------------------------------------------------------------------- #
N="$WORK/refuse_atomic"; rm -rf "$N"; mkdir -p "$N/output"
cp "$A/input.inp" "$A/base.inp" "$N/"
sed -i -e 's/^Load IC?.*/Load IC? True/' "$N/input.inp"
( cd "$N" && OMP_NUM_THREADS=1 EXHALE_MOLECULAR_SEED="$A/output" "$EXE" > run.log 2>&1 )
if [ $? -eq 0 ]; then s=accepted; else s=refused; fi
chk_eq molecular_seed_refuses_atomic_target "$s" refused

H="$WORK/refuse_no_handoff"; rm -rf "$H"; mkdir -p "$H/output"
sed -e 's/^Load IC?.*/Load IC? True/' "$REG/input.inp" > "$H/input.inp"
( cd "$H" && OMP_NUM_THREADS=1 EXHALE_MOLECULAR_SEED="$A/output" "$EXE" > run.log 2>&1 )
if [ $? -eq 0 ]; then s=accepted; else s=refused; fi
chk_eq molecular_seed_refuses_without_handoff "$s" refused

# -------------------------------------------------------------------------- #
# T-L7e: what the 'local' partition means since item L7e -- the SMALLER of the
# thermochemical fit and the root of each cell's own H2 carrier row.  The
# rows below are statements about the seed the binary writes, read from its
# own report (EXHALE_CARRIER_DEBUG=1 prints the two side by side) and from the
# state file.
# -------------------------------------------------------------------------- #
mkcase local
( cd "$WORK/local" && OMP_NUM_THREADS=1 EXHALE_CARRIER_DEBUG=1 \
     EXHALE_MOLECULAR_SEED="$A/output" EXHALE_MOLECULAR_SEED_X2=local \
     "$EXE" > run.log 2>&1 )
chk_eq molecular_seed_local_exit "$?" 0
if [ ! -f "$WORK/local/output/Ion_species_IC.txt" ]; then
   echo "FAIL molecular_seed_local_written measured=absent reference=present tol=0"
   tail -n 20 "$WORK/local/run.log"; FAIL=1
else
   # THE ADOPTED FRACTION IS NEVER ABOVE EITHER STATEMENT, cell by cell, on
   # every line of the report.  That is the definition -- the smaller of the
   # thermochemical fit and the root of the cell's own row, in EVERY cell --
   # and it is the one row that would catch the partition being taken from
   # the wrong one of the two.
   over=$(awk '/fit against the row/{f=1;next} f&&/^ *[0-9]+ /{
             fit=$4; root=$5; ad=$6;
             if (ad > fit*(1+1e-12)) n++;
             if (root > 0 && ad > root*(1+1e-12)) n++ }
             END{print n+0}' "$WORK/local/run.log")
   chk_eq molecular_seed_local_is_the_smaller_of_the_two "${over:-1}" 0
   # AND THE LAYER THE HANDOFF STILL DESCRIBES IS REPORTED, with how far it
   # stands from the network's own root there: a finding the seed records
   # and does not act on (plan item L7f).
   jtop=$(sed -n 's/.*layer the handoff still describes: cells 1 to \([0-9]*\),.*/\1/p' \
          "$WORK/local/run.log" | tail -n 1)
   echo "  DIAGNOSTIC layer the handoff describes on this fixture: cells 1..${jtop:-0}"
   # AND IT IS THE SMALLER ONE WHERE THEY DIFFER: at least one line of the
   # report must have the root below the fit, or the row above is vacuous on
   # this fixture and the test says nothing.
   dif=$(awk '/fit against the row/{f=1;next} f&&/^ *[0-9]+ /{
             if ($5 > 0 && $5 < $4*(1-1e-12)) n++ } END{print n+0}' \
             "$WORK/local/run.log")
   if [ "${dif:-0}" -gt 0 ]; then s=yes; else s=no; fi
   chk_eq molecular_seed_local_root_governs_somewhere "$s" yes
   # AND THE NUMBER CALLED A ROOT IS ONE.  The H2 row is not linear in the
   # unknown -- its dominant formation channel is the three-body
   # association k15 n(H I)^2 and atomic hydrogen is closed out of the
   # element budget, so the production falls as n(H2) rises and vanishes
   # where every nucleus is bound.  Production over loss rate at the trial
   # is therefore NOT the root, and evaluated at a fully molecular trial it
   # collapses with n(H I)^2 (item L7f).  The seed reports |P - L|/(P + L)
   # at the density it hands back, and this row is the statement that the
   # row actually balances there, at every radius the report prints.
   worst=$(awk '/fit against the row/{f=1;next} f&&/^ *[0-9]+ /{
              if ($7+0 > m) m=$7+0 } END{printf "%.3e", m+0}'               "$WORK/local/run.log")
   chk_le molecular_seed_local_root_balances_the_row "${worst:-1}" 1.0e-8
   # THE REPORT AND THE STATE FILE SAY THE SAME THING about how far the root
   # governs: the count on the log and the count in the file's own record
   # line are one number.
   nlog=$(sed -n "s/.*root is the smaller in \([0-9]*\) cell(s).*/\1/p" \
          "$WORK/local/run.log" | tail -n 1)
   nfil=$(sed -n 's/^# molecular_partition_local: cells where the row.s own root is the smaller *\([0-9]*\).*/\1/p' \
          "$WORK/local/output/Ion_species_IC.txt" | tail -n 1)
   chk_eq molecular_seed_local_crossover_is_recorded "${nfil:-none}" "${nlog:-none}"
   # AND THE CONVERSION IS STILL A CONVERSION: the nuclei, the helium and the
   # mass come through it unchanged, as in every other mode.
   D=$(sed -n 's/.*max relative change of H nuclei, He nuclei, mass: *//p' \
       "$WORK/local/run.log" | tail -n 1)
   chk_le molecular_seed_local_conserves_H_nuclei "$(echo $D | awk '{print $1}')" 1.0e-12
   chk_le molecular_seed_local_conserves_mass "$(echo $D | awk '{print $3}')" 1.0e-12
fi
E=$(sed -n 's/.*p round trip): *//p' "$WORK/local/run.log" | tail -n 1 | tr -d ' ')
chk_le molecular_seed_local_eos_closure "${E:-1}" 1.0e-12

exit $FAIL
