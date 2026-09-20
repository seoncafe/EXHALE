#!/bin/bash
# THE MOLECULAR SEED OF A CERTIFIED ATOMIC STATE IS BUILT, AND ITS GHOSTS
# CLOSE (item D9fix; the measurements are
# docs/lhs1140b_stationary_D9fix_20260919.md).
#
# WHAT IS UNDER TEST, row by row.
#
#   1  THE BASE HANDOFF CLOSES ON A GHOST WHOSE FIXED-POINT MAP CONTRACTS
#      SLOWLY. The ghost owes x_H2 = x2 (1 - x_ion) with the x_ion of its own
#      ionization balance at fixed T. On the ghost of this state the map
#      x -> x2 (1 - x_ion(x)) has slope 0.74 (MEASURED), so substitution
#      reaches 1.2e-6 in 30 passes against the 1e-10 the closure asks for
#      and the run stops with "base ghost H2 closure failed". The row reads
#      whether that stop is printed.
#
#   2  THE CONVERSION 2 H -> H2 CONSERVES THE HYDROGEN NUCLEI OF THE GHOSTS.
#      load_IC does not read the ghost rows of the atomic pair and fills them
#      with the molecular run's reservoir, which already carries H2; a
#      transfer that partitions H I alone and overwrites f_H2 there destroys
#      those nuclei (MEASURED: relative change 0.99998). The row reads the
#      largest relative change of the H nuclei the seed report prints and
#      holds it to molecular_seed_closure_tol = 1e-12, the conversion's own
#      acceptance.
#
# THE FIXTURE
#   The certified generation of LHS1140b/models/atomic_scalar_gj1132_kzz1e9/
#   HeH2.13 that the catalog refresh of 2026-09-19 converted, read from its
#   published directory (read-only by the publisher, so it does not move),
#   and the input of molecular_scalar_gj1132_kzz1e9/HeH2.13 written below.
#   Nothing in LHS1140b/ is written to.
#
# Usage: molecular_seed_of_an_atomic_state.sh   (EXHALE_EXE selects the
#        binary, EXHALE_TEST_OUT the work directory)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" molecular_seed_of_an_atomic_state
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/boundary_state}"
WORK="$OUT/molecular_seed_of_an_atomic_state"
GEN="$ROOT/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/states/g0003_20260919T003908Z_c1794389"
SED="$ROOT/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt"

n_fail=0
check_eq() {   # check_eq <name> <measured> <reference>
   if [ "$2" = "$3" ]; then
      echo "PASS $1 measured=$2 reference=$3 tol=0"
   else
      echo "FAIL $1 measured=$2 reference=$3 tol=0"
      n_fail=$((n_fail + 1))
   fi
}
check_le() {   # check_le <name> <measured> <tolerance>
   local v
   v=$(awk -v m="${2:-}" -v t="$3" 'BEGIN {
          if (m == "") { print "missing" }
          else if (m + 0 <= t + 0) { print "within" } else { print "above" } }')
   if [ "$v" = "within" ]; then
      echo "PASS $1 measured=${2:-missing} reference=at_most tol=$3"
   else
      echo "FAIL $1 measured=${2:-missing} reference=at_most tol=$3"
      n_fail=$((n_fail + 1))
   fi
}

for f in "$EXE" "$GEN/Hydro_ioniz.txt" "$GEN/Ion_species.txt" "$SED"; do
   if [ ! -e "$f" ]; then
      echo "FAIL molecular_seed_of_an_atomic_state measured=missing reference=$f tol=0"
      exit 1
   fi
done

rm -rf "$WORK"
mkdir -p "$WORK/src" "$WORK/output"
cp "$GEN/Hydro_ioniz.txt" "$WORK/src/Hydro_ioniz_IC.txt"
cp "$GEN/Ion_species.txt" "$WORK/src/Ion_species_IC.txt"
cat > "$WORK/input.inp" <<INP
Planet name: LHS1140b_molecular_scalar_gj1132_kzz1e9_HeH2.13
Planet radius [R_J]: 0.157692
Planet mass [M_J]: 0.0176220
Equilibrium temperature [K]: 226.0
Orbital distance [AU]: 0.0946
Escape radius [R_p]: 2.00
He/H number ratio: 2.13
2D approximate method: Mdot
Parent star mass [M_sun]: 0.1844
Spectrum type: Load from file..
Spectrum file: $SED
Use only EUV? False
[E_low,E_mid,E_high] = [ 13.60 , 123.98 , 1.24e3 ]
Log10 of X-ray luminosity [erg/s]: 26.372
Log10 of EUV luminosity [erg/s]: 26.404
Grid type: Mixed
Base grid [dr,cells]: 1.9999999494757503e-4 50
Numerical flux: HLLC
Reconstruction scheme: PLM
Include He23S? True
Load IC? True
Do only PP: False
Force start: False
He_diffusion: True
Domain mode: Spherical
Outer radius [R_p]: 30.0
Stellar radius [R_sun]: 0.2159
Base BC: pressure 1.0
Resid tol: 5.0e-5
He_Kzz: 1.0e9
Solver: Newton
Restart intent: stationary
Secondary_ionization: Immediate
Well balanced: True
Molecular chemistry: True
Molecular carrier transport: True
Coupled carrier solve: False
INP
cat > "$WORK/base.inp" <<INP
T_base    226.0
r_base    0.157692
HeH_base  2.13
Kzz_base  1.0e9
q_H2_base 0.19011025828621544
p_base    1.0e-6
INP

( cd "$WORK" && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
     EXHALE_MOLECULAR_SEED="$WORK/src" EXHALE_MOLECULAR_SEED_X2=local \
     "$EXE" > seed.log 2>&1 )

if grep -q 'base ghost H2 closure failed' "$WORK/seed.log"; then
   closure=failed
   grep -A3 'STOP: the base handoff partition' "$WORK/seed.log" | sed 's/^/  /'
else
   closure=closed
fi
check_eq ghost_h2_closure_of_the_molecular_seed "$closure" closed

dnH=$(awk '/max relative change of H nuclei/ { print $(NF-2) }' "$WORK/seed.log" | head -n 1)
[ -z "$dnH" ] && dnH=$(awk '/max relative change of the H nuclei/ { print $(NF-2) }' "$WORK/seed.log" | head -n 1)
check_le hydrogen_nuclei_of_the_molecular_seed "$dnH" 1.0e-12
[ -s "$WORK/output/Ion_species_IC.txt" ] && \
   echo "  seed pair written: $WORK/output/{Hydro_ioniz,Ion_species}_IC.txt"

echo ""
if [ $n_fail -gt 0 ]; then
   echo "molecular_seed_of_an_atomic_state: $n_fail assertion(s) failed"
   exit 1
fi
exit 0
