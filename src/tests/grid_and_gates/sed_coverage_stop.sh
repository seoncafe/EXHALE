#!/bin/bash
# A loaded SED that stops above the photon-grid floor an active absorber
# needs makes the run stop.
#
# QUANTITY UNDER TEST
#   The exit status of EXHALE.x, and the text it prints, when
#   "Spectrum type: Load from file.." names a table whose longest
#   wavelength is shorter than the wavelength of the lowest ionization
#   threshold that is active in the run.  read_sed
#   (src/modules/radiation/sed_read.f90) lowers the grid floor e_low to
#   photon_grid_floor_eV, the lowest ionization threshold over the active
#   absorbers:
#     H(n=2)       3.400 eV = 3647 A, whenever the excited-hydrogen
#                  coupling is armed, i.e. whenever both "Stellar Teff
#                  [K]:" and "Stellar radius [R_sun]:" are stated;
#     He 2^3S      4.768 eV = 2600 A, when the metastable is included;
#     a low-IP metal's neutral stage (K I at 4.341 eV = 2856 A is the
#                  lowest one the species table carries).
#   A table that ends above that floor leaves the absorber with no field
#   over part of its own band, so its photoionization rate is the rate of
#   a truncated spectrum and not of the stated one.
#
# REFERENCE
#   docs/development_plan_20260905_rev3.md section 10.5 decision 17: the
#   run STOPS, there is no key to continue, and the message names the
#   file, the missing band in eV and A, the absorbers that need it, and
#   the remedies.  H(n=2) is an absorber like the others (decision 13), so
#   the stop is about the FILE and not about the window the code opens.
#
# WHAT IS RUN (five short runs of benchmarks/wasp52, none of them in the
# benchmark directory itself; each capped at one step)
#   A  the benchmark as it stands -- "Include He23S? True", the two
#      stellar lines, no metals.inp -- with SED = inputdata/
#      scaled_solar_hd189.sed, whose last row is 1899.5 A = 6.527 eV.
#      Floor 3.400 eV.  Must exit nonzero and name the file, He 2^3S and
#      H(n=2), and the remedy wavelength 3647 A.
#   B  the same with "Include He23S? False", no low-IP metal and the two
#      stellar lines removed: no absorber below the H I edge is active, so
#      the floor is the 13.60 eV of the [E_low,E_mid,E_high] line, which
#      the table covers.  Must exit 0 and write output/Hydro_ioniz.txt.
#   C  the same as B plus a metals.inp carrying Na (5.139 eV) and K
#      (4.341 eV).  Floor 4.341 eV.  Must exit nonzero and name both
#      elements; this is the branch no shipped metals.inp reaches, since
#      none of them lists K.
#   D  "Include He23S? False", the two stellar lines kept, and the eps Eri
#      table CUT at 2600 A: the only absorber below the file's lowest
#      photon is then H(n=2).  Must exit nonzero, name H(n=2) and NOT name
#      He 2^3S, and offer the stellar-line remedy.
#   E  the benchmark as it stands, with its own eps Eri table, which
#      reaches 3700 A: every absorber's band is covered.  Must exit 0 and
#      write output/Hydro_ioniz.txt.  This is the WASP-52b configuration.
#
# ASSERTIONS (eleven)
#   sed_coverage_metastable_stops       A exits nonzero
#   sed_coverage_metastable_message     A's message names the file, the
#                                       missing band, He 2^3S and H(n=2)
#   sed_coverage_remedy_wavelength      the wavelength A tells the user to
#                                       supply, against 3647 A, tol 1 A
#   sed_coverage_covered_floor_runs     B exits 0
#   sed_coverage_covered_floor_output   B wrote output/Hydro_ioniz.txt
#   sed_coverage_low_ip_metal_stops     C exits nonzero
#   sed_coverage_low_ip_metal_message   C's message names Na I and K I
#   sed_coverage_excited_h_stops        D exits nonzero
#   sed_coverage_excited_h_message      D names H(n=2) and the stellar-line
#                                       remedy, and does not name He 2^3S
#   sed_coverage_balmer_table_runs      E exits 0
#   sed_coverage_balmer_table_output    E wrote output/Hydro_ioniz.txt
#
# EXPECTED BEFORE THE 2c-N2FLOOR CHANGE: RED on the four H(n=2)
# assertions (D stopped in excited_hydrogen instead of read_sed and named
# no file; E stopped although its file reaches 3647 A) and on
# sed_coverage_covered_floor_runs (B kept the stellar lines and stopped).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${EXHALE_EXE:-$ROOT/EXHALE.x}"
WORK="$ROOT/build/tests/grid_and_gates"
CASE="$ROOT/benchmarks/wasp52"
SED_FILE="$ROOT/inputdata/scaled_solar_hd189.sed"
EPS_FILE="$ROOT/WASP-52b/eps_eri_sed_fxuv1p0.txt"

n_fail=0
verdict() {  # verdict <PASS|FAIL> <name> <measured> <reference> <tol>
   echo "$1 $2 measured=$3 reference=$4 tol=$5"
   [ "$1" = FAIL ] && n_fail=$((n_fail+1))
   return 0
}

if [ ! -x "$EXE" ]; then
   verdict FAIL sed_coverage_binary no_binary "$EXE" 0
   exit 1
fi
for f in "$SED_FILE" "$EPS_FILE"; do
   if [ ! -f "$f" ]; then
      verdict FAIL sed_coverage_sed_file missing "$f" 0
      exit 1
   fi
done

# The truncation the test relies on: the last row of each table, in A.
last_w="$(grep -v '^ *#' "$SED_FILE" | awk 'NF>=2{w=$1} END{print w}')"
last_e="$(grep -v '^ *#' "$EPS_FILE" | awk 'NF>=2{w=$1} END{print w}')"
echo "  short SED file     : $SED_FILE"
echo "  longest wavelength : $last_w A"
echo "  eps Eri SED file   : $EPS_FILE"
echo "  longest wavelength : $last_e A"

for d in sedA sedB sedC sedD sedE; do
   rm -rf "$WORK/$d"
   mkdir -p "$WORK/$d/output"
   cp "$CASE/input.inp" "$CASE/base.inp" "$WORK/$d/"
done
for d in sedA sedB sedC; do
   sed -i "s|^Spectrum file:.*|Spectrum file: $SED_FILE|" "$WORK/$d/input.inp"
done
sed -i 's/^Include He23S?.*/Include He23S? False/' "$WORK/sedB/input.inp" \
        "$WORK/sedC/input.inp" "$WORK/sedD/input.inp"
# B and C carry no absorber below the H I edge at all, so the two stellar
# lines that arm the excited-hydrogen coupling come out as well.
sed -i '/^Stellar Teff/d;/^Stellar radius/d' "$WORK/sedB/input.inp" \
        "$WORK/sedC/input.inp"
# Solar Na (Asplund 2009) and K; both neutral thresholds are below the
# lowest photon the table carries, K I at 4.341 eV being the lowest
# neutral-metal threshold of the species table.
printf 'Na    1.74e-6\nK     1.32e-7\n' > "$WORK/sedC/metals.inp"
# D: the eps Eri table cut at 2600 A, so that its lowest photon (4.77 eV)
# is above the n=2 edge and below every other threshold that is off.
awk 'NF>=2 && $1+0 <= 2600.0' "$EPS_FILE" > "$WORK/sedD/eps_eri_cut2600.txt"
sed -i "s|^Spectrum file:.*|Spectrum file: $WORK/sedD/eps_eri_cut2600.txt|" \
       "$WORK/sedD/input.inp"
sed -i "s|^Spectrum file:.*|Spectrum file: $EPS_FILE|" "$WORK/sedE/input.inp"

run_stage() {  # run_stage <dir>; echoes the exit status
   ( cd "$WORK/$1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" \
        > run.log 2>&1 )
   echo $?
}

rcA=$(run_stage sedA)
rcB=$(run_stage sedB)
rcC=$(run_stage sedC)
rcD=$(run_stage sedD)
rcE=$(run_stage sedE)
echo "  exit status: A(metastable)=$rcA  B(covered floor)=$rcB  C(Na,K)=$rcC"
echo "               D(H(n=2) only)=$rcD  E(eps Eri, covered)=$rcE"

# ---- stage A: the metastable and n=2 floors together ----
if [ "$rcA" -ne 0 ]; then
   verdict PASS sed_coverage_metastable_stops "exit_$rcA" nonzero 0
else
   verdict FAIL sed_coverage_metastable_stops "exit_$rcA" nonzero 0
fi
msgA="$WORK/sedA/run.log"
missA=""
grep -q "scaled_solar_hd189.sed" "$msgA" || missA="$missA file_name"
grep -q "band not covered"       "$msgA" || missA="$missA band"
grep -q "He 2^3S"                "$msgA" || missA="$missA absorber_He23S"
grep -q "H(n=2)"                 "$msgA" || missA="$missA absorber_Hn2"
grep -q "Include He23S? False"   "$msgA" || missA="$missA remedy_key"
grep -q "Stellar Teff"           "$msgA" || missA="$missA remedy_stellar"
grep -qi "supply an SED"         "$msgA" || missA="$missA remedy_sed"
if [ -z "$missA" ]; then
   verdict PASS sed_coverage_metastable_message all_named \
           "file,band,He_2^3S,H(n=2),3_remedies" 0
else
   verdict FAIL sed_coverage_metastable_message "missing:${missA# }" \
           "file,band,He_2^3S,H(n=2),3_remedies" 0
fi
# The wavelength the message asks for is the n=2 edge, hc/e_th_HI_n2 with
# e_th_HI_n2 = e_th_HI/4 = 3.399609 eV, i.e. 3647.02 A.  1 A is well below
# the spacing of any tabulated SED near that wavelength.
w_ask="$(awk '/supply an SED that reaches/{print $(NF-1)}' "$msgA" | tail -n 1)"
w_dev="$(awk -v w="${w_ask:-0}" 'BEGIN{d=w-3647.02; if(d<0)d=-d; print d}')"
if awk -v d="$w_dev" 'BEGIN{exit !(d<=1.0)}'; then
   verdict PASS sed_coverage_remedy_wavelength "${w_ask:-none}" 3647.02 1.0
else
   verdict FAIL sed_coverage_remedy_wavelength "${w_ask:-none}" 3647.02 1.0
fi

# ---- stage B: a floor the table covers ----
# Exit status 2 = the run declared a stationary state that the A2 certification
# refused (docs/a2_certification_contract_20260906.md section 8); the run wrote
# its outputs in full, which is what this gate reads. Only 1 (a Fortran error
# stop) or a signal is a failed run here.
if [ "$rcB" -eq 0 ] || [ "$rcB" -eq 2 ]; then
   verdict PASS sed_coverage_covered_floor_runs "exit_$rcB" exit_0 0
else
   verdict FAIL sed_coverage_covered_floor_runs "exit_$rcB" exit_0 0
   tail -n 5 "$WORK/sedB/run.log" | sed 's/^/     /'
fi
if [ -s "$WORK/sedB/output/Hydro_ioniz.txt" ]; then
   verdict PASS sed_coverage_covered_floor_output written written 0
else
   verdict FAIL sed_coverage_covered_floor_output no_output written 0
fi

# ---- stage C: the low-IP metal floor ----
if [ "$rcC" -ne 0 ]; then
   verdict PASS sed_coverage_low_ip_metal_stops "exit_$rcC" nonzero 0
else
   verdict FAIL sed_coverage_low_ip_metal_stops "exit_$rcC" nonzero 0
fi
msgC="$WORK/sedC/run.log"
missC=""
grep -q "Na I" "$msgC"          || missC="$missC NaI"
grep -q "K I"  "$msgC"          || missC="$missC KI"
grep -qi "4.34100E+00" "$msgC"  || missC="$missC floor_eV"
grep -q "metals.inp" "$msgC"    || missC="$missC remedy_metals"
if [ -z "$missC" ]; then
   verdict PASS sed_coverage_low_ip_metal_message all_named \
           "Na_I,K_I,floor_eV,metals.inp_remedy" 0
else
   verdict FAIL sed_coverage_low_ip_metal_message "missing:${missC# }" \
           "Na_I,K_I,floor_eV,metals.inp_remedy" 0
fi

# ---- stage D: H(n=2) alone below the file's lowest photon ----
if [ "$rcD" -ne 0 ]; then
   verdict PASS sed_coverage_excited_h_stops "exit_$rcD" nonzero 0
else
   verdict FAIL sed_coverage_excited_h_stops "exit_$rcD" nonzero 0
fi
msgD="$WORK/sedD/run.log"
missD=""
grep -q "eps_eri_cut2600.txt" "$msgD" || missD="$missD file_name"
grep -q "H(n=2)"              "$msgD" || missD="$missD absorber_Hn2"
grep -q "Stellar Teff"        "$msgD" || missD="$missD remedy_stellar"
grep -q "He 2^3S"             "$msgD" && missD="$missD named_He23S_wrongly"
if [ -z "$missD" ]; then
   verdict PASS sed_coverage_excited_h_message all_named \
           "file,H(n=2),stellar_remedy,no_He_2^3S" 0
else
   verdict FAIL sed_coverage_excited_h_message "wrong:${missD# }" \
           "file,H(n=2),stellar_remedy,no_He_2^3S" 0
fi

# ---- stage E: a table that reaches the n=2 edge ----
if [ "$rcE" -eq 0 ] || [ "$rcE" -eq 2 ]; then
   verdict PASS sed_coverage_balmer_table_runs "exit_$rcE" exit_0 0
else
   verdict FAIL sed_coverage_balmer_table_runs "exit_$rcE" exit_0 0
   tail -n 5 "$WORK/sedE/run.log" | sed 's/^/     /'
fi
if [ -s "$WORK/sedE/output/Hydro_ioniz.txt" ]; then
   verdict PASS sed_coverage_balmer_table_output written written 0
else
   verdict FAIL sed_coverage_balmer_table_output no_output written 0
fi

exit $(( n_fail > 0 ))
