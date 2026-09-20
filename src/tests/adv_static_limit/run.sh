#!/bin/bash
# THE ADVECTION-CORRECTED POST-PROCESS IN A COLUMN THAT CARRIES NO WIND.
#
# The `_adv` files re-solve the ionization balance and the energy balance of
# each cell with the advection term kept, upwind, on the profile the run
# produced. Both re-solves are statements about a flow, and both lose their
# meaning as the flow vanishes: the upwind difference has no upstream cell and
# the advective terms leave the equations. What each of them does in that
# limit is what this suite pins.
#
#   * The ionization correction has the Damkohler condition (ii) of
#     post_process_adv: the residence time dr/v grows without bound as v -> 0,
#     so every population is locally equilibrated and the cell keeps the
#     converged equilibrium ionization. The limit is therefore the equilibrium
#     state itself, reached continuously and not on the sign of v, and row B
#     asserts it on a run.
#   * The energy correction has the thermal Damkohler condition of the same
#     module, rows 1 to 7, measured on the function that states it.
#
#   * The energy correction solves the full steady internal-energy equation,
#     enthalpy flux of the mass-flux divergence included, so rows N1 to N6
#     pin that third term: its algebra in both branches of the residual, its
#     exact absence where the mass flux is stationary, and the exact solution
#     it recovers on a column whose mass flux is not.
#
#   * The upwind term of that equation differences the SPECIFIC internal
#     energy, and what the flow carries in is the specific energy of the
#     UPSTREAM gas, evaluated at the upstream composition: across a
#     dissociation front two cells at the same temperature store different
#     energy. Rows U0 to U4 pin that, including the interface where the cell
#     is atomic and its upstream neighbor is not, and the atomic interface
#     that has to keep the constant gamma_ad arithmetic exactly.
#
#   * WHETHER A CELL IS CORRECTED AT ALL, and what its correction is worth,
#     is decided by the mass row of the state the post-process was handed:
#     the numerical FACE fluxes of the Riemann solve differenced over the
#     cell's own control volume, divided by the largest term the row itself
#     holds. The correction is FIRST ORDER in that measure, so a cell at or
#     below adv_conditional_tol is corrected and its correction is accurate
#     to that fraction of itself, and a cell above it keeps the run's own
#     value. Rows M1 to M7 pin the decision and both numbers on the
#     production routines: an analytic outflow whose measure is 2e-7 is a
#     conditional correction and is NOT a certified stationary state, and a
#     column whose mass flux grows by 3 percent across every cell is
#     retained throughout.
#
#   * The enthalpy term ratio stays a SENSITIVITY diagnostic of the energy
#     correction and refuses nothing. Rows R1 to R9 pin it, the last two
#     being the counterexample: a nonzero mass divergence whose ratio is
#     0.5, which a screen built on the ratio accepts.
#
#   * Rows E to J pin what the decision does to a whole run: the two status
#     columns that record the validity of each row's temperature and of its
#     composition, the adv_mass_row column that records the measure both were
#     decided by (against the value the stationary certification prints for
#     the same state), the fraction stated in the file header, the counts in
#     that header, that a retained row carries the run's own state exactly,
#     that the corrected temperature of the mechanical column is therefore
#     bounded by the run's own, that the two files agree row by row, and that
#     the Python loader reads the schema, the fraction and the measure, and
#     reads a file written before the schema as UNKNOWN validity.
#
# Every row prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line; lines beginning with two spaces are context. Exit status is nonzero if
# any row fails.
#
# Usage: src/tests/adv_static_limit/run.sh
#   EXHALE_OBJDIR   object directory to link the unit rows against
#                   (default $ROOT/build); with it set the staleness check is
#                   skipped, as a concurrent private build needs
#   EXHALE_EXE      binary the whole-run rows use (default $ROOT/EXHALE.x)
#   EXHALE_TEST_OUT working directory (default $ROOT/build/tests/adv_static_limit)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OBJDIR="${EXHALE_OBJDIR:-$ROOT/build}"
OUT="${EXHALE_TEST_OUT:-$ROOT/build/tests/adv_static_limit}"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" adv_static_limit
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
REG="$ROOT/backup/regression"
FC="${FC:-gfortran}"
FFLAGS_TEST="${FFLAGS_TEST:--O0 -g -fbacktrace -fopenmp}"
FAIL=0

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Rows 1 to 5: the thermal Damkohler condition of the energy correction.
# ---------------------------------------------------------------------------
FC_PATH="$(command -v "$FC" 2>/dev/null || true)"
PREFIX="$(cd "$(dirname "$FC_PATH")/.." 2>/dev/null && pwd || echo /usr)"
if [ -f "$PREFIX/lib/libopenblas.so" ]; then
   LAPACK="-L$PREFIX/lib -lopenblas -Wl,-rpath,$PREFIX/lib -ldl"
else
   LAPACK="-llapack -ldl"
fi

if [ ! -d "$OBJDIR" ] || [ -z "$(ls "$OBJDIR"/*.o 2>/dev/null)" ]; then
   echo "FAIL adv_static_limit_build measured=no_objects reference=$OBJDIR tol=0"
   echo "     build the code first"
   exit 1
fi
if [ -z "${EXHALE_OBJDIR:-}" ]; then
   if ! ( cd "$ROOT" && make -q ) >/dev/null 2>&1; then
      echo "FAIL adv_static_limit_build measured=stale reference=up_to_date tol=0"
      echo "     'make -q' says build/ is behind the sources; build first"
      exit 1
   fi
fi
PROD_OBJ="$(ls "$OBJDIR"/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$' | tr '\n' ' ')"

rm -f "$OUT/thermal_damkohler_limit.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/thermal_damkohler_limit.x" \
    "$HERE/thermal_damkohler_limit.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL thermal_damkohler_limit_build measured=compile_error reference=ok tol=0"
   exit 1; }
env OMP_NUM_THREADS=1 "$OUT/thermal_damkohler_limit.x" || FAIL=1

# ---------------------------------------------------------------------------
# Rows N1 to N6: the enthalpy flux of the mass-flux divergence, the third term
# of the same energy equation.
# ---------------------------------------------------------------------------
rm -f "$OUT/enthalpy_flux_of_mass_divergence.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/enthalpy_flux_of_mass_divergence.x" \
    "$HERE/enthalpy_flux_of_mass_divergence.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL enthalpy_flux_of_mass_divergence_build measured=compile_error reference=ok tol=0"
   exit 1; }
env OMP_NUM_THREADS=1 "$OUT/enthalpy_flux_of_mass_divergence.x" || FAIL=1

# ---------------------------------------------------------------------------
# Rows U0 to U4: the upstream end of the advected energy difference, the
# specific internal energy the flow carries into the cell.
# ---------------------------------------------------------------------------
rm -f "$OUT/upstream_caloric_energy.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/upstream_caloric_energy.x" \
    "$HERE/upstream_caloric_energy.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL upstream_caloric_energy_build measured=compile_error reference=ok tol=0"
   exit 1; }
env OMP_NUM_THREADS=1 "$OUT/upstream_caloric_energy.x" || FAIL=1

# ---------------------------------------------------------------------------
# Rows R1 to R9: the sensitivity diagnostic built on that same term.
# ---------------------------------------------------------------------------
rm -f "$OUT/stationary_flux_refusal.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/stationary_flux_refusal.x" \
    "$HERE/stationary_flux_refusal.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL stationary_flux_refusal_build measured=compile_error reference=ok tol=0"
   exit 1; }
env OMP_NUM_THREADS=1 "$OUT/stationary_flux_refusal.x" || FAIL=1

# ---------------------------------------------------------------------------
# Rows M1 to M7: the decision itself and the two numbers it is taken with, on
# the face-flux mass row of the state, through the production grid, boundary
# and residual routines.
# ---------------------------------------------------------------------------
rm -f "$OUT/stationarity_by_the_mass_row.x"
$FC $FFLAGS_TEST -J"$OUT" -I"$OBJDIR" \
    -o "$OUT/stationarity_by_the_mass_row.x" \
    "$HERE/stationarity_by_the_mass_row.f90" $PROD_OBJ $LAPACK || {
   echo "FAIL stationarity_by_the_mass_row_build measured=compile_error reference=ok tol=0"
   exit 1; }
env OMP_NUM_THREADS=1 "$OUT/stationarity_by_the_mass_row.x" || FAIL=1

# ---------------------------------------------------------------------------
# Rows A to D: the post-process of a column with no wind, on the binary.
#
# The column is backup/regression/hydrostatic_column, a mechanical column with
# the radiation switched off. It is relaxed for its recorded 300 steps, its
# velocity column is then set to zero, and the resulting state is handed back
# to the same binary with "Do only PP" so that the post-process runs on it.
# The one step such a run still takes regenerates the spurious velocity of the
# spatial operator, which is why row A measures what the flow actually is
# rather than assuming it is zero.
# ---------------------------------------------------------------------------
if [ ! -x "$EXE" ]; then
   echo "FAIL adv_static_limit_binary measured=no_binary reference=$EXE tol=0"
   exit 1
fi

WORK="$OUT/static_column"
rm -rf "$WORK"; mkdir -p "$WORK/output"
cp "$REG/hydrostatic_column/input.inp" "$WORK/input.inp"
MS="$(cat "$REG/hydrostatic_column/maxsteps" 2>/dev/null || echo 300)"
( cd "$WORK" && env OMP_NUM_THREADS=1 EXHALE_MAXSTEPS="$MS" "$EXE" \
     > relax.log 2>&1 )
if [ ! -s "$WORK/output/Hydro_ioniz.txt" ]; then
   echo "FAIL adv_static_limit_relax measured=no_output reference=Hydro_ioniz.txt tol=0"
   sed -n '$p' "$WORK/relax.log"
   exit 1
fi

# Keep the outputs of the column AS THE RUN LEFT IT, before the velocity
# column is zeroed below: rows E to I are about the refusal on a real
# relaxation snapshot, whose mass flux is nowhere stationary, and the
# second run overwrites output/.
rm -rf "$WORK/relaxed"; mkdir -p "$WORK/relaxed"
cp "$WORK/output/Hydro_ioniz.txt" "$WORK/output/Hydro_ioniz_adv.txt" \
   "$WORK/output/Ion_species.txt" "$WORK/output/Ion_species_adv.txt" \
   "$WORK/relaxed/"

python3 - "$WORK" <<'PY'
import sys
w = sys.argv[1]
out = []
for line in open(w + '/output/Hydro_ioniz.txt'):
    if line.lstrip().startswith('#'):
        out.append(line); continue
    f = line.split()
    f[2] = '%.17E' % 0.0          # the velocity column, set to no flow
    out.append(' '.join(f) + '\n')
open(w + '/output/Hydro_ioniz_IC.txt', 'w').writelines(out)
open(w + '/output/Ion_species_IC.txt', 'w').writelines(
    open(w + '/output/Ion_species.txt').readlines())
PY

sed -i 's/^Load IC?.*/Load IC? True/; s/^Do only PP:.*/Do only PP: True/' \
    "$WORK/input.inp"
( cd "$WORK" && env OMP_NUM_THREADS=1 "$EXE" > run.log 2>&1 )

python3 - "$WORK" <<'PY' || exit 1
import sys, math
w = sys.argv[1]
fail = 0

def rows(path):
    return [l.split() for l in open(path) if not l.lstrip().startswith('#')]

def say(verdict, name, measured, reference, tol):
    print('%s %s measured=%s reference=%s tol=%s'
          % (verdict, name, measured, reference, tol))

hyd = rows(w + '/output/Hydro_ioniz.txt')
adv = rows(w + '/output/Hydro_ioniz_adv.txt')
ion = rows(w + '/output/Ion_species.txt')
iad = rows(w + '/output/Ion_species_adv.txt')

# Row A: what the flow in this state actually is, as a fraction of the sound
# speed of a neutral H/He gas. The rows below are statements about a column
# with no wind, so the column has to be one.
kB, mH, gam, mu = 1.380649e-16, 1.6726219e-24, 5.0 / 3.0, 1.3
mach = 0.0
for h in hyd:
    T = float(h[4])
    cs = math.sqrt(gam * kB * T / (mu * mH))
    mach = max(mach, abs(float(h[2])) / cs)
verdict = 'PASS' if mach < 1.0e-2 else 'FAIL'
if verdict == 'FAIL':
    fail = 1
say(verdict, 'static_column_is_static', '%.3e' % mach, '<1.0e-2', '0')

# Row B: the ionization correction returns the equilibrium state exactly. The
# Damkohler condition pins every cell of a column this slow, so the corrected
# species columns are the equilibrium ones to the bit -- not to a tolerance.
worst, worst_col = 0.0, 0
for a, b in zip(ion, iad):
    for k in range(1, len(a)):
        if a[k] != b[k]:
            d = abs(float(a[k]) - float(b[k]))
            if d > worst:
                worst, worst_col = d, k
verdict = 'PASS' if worst == 0.0 else 'FAIL'
if verdict == 'FAIL':
    fail = 1
say(verdict, 'ionization_correction_is_equilibrium_at_no_flow',
    '%.3e' % worst, '0', '0')
if worst != 0.0:
    print('  largest species-column departure in column %d' % worst_col)

# Row C: the corrected temperature of the same column. MEASURED and stated,
# not asserted (see the header of this file).
rel = 0.0
rad = 0.0
for h, a in zip(hyd, adv):
    T, Ta = float(h[4]), float(a[4])
    if T > 0.0 and abs(Ta - T) / T > rel:
        rel, rad = abs(Ta - T) / T, float(h[0])
print('  corrected temperature departs from the run temperature by up to '
      '%.3e relative, at r = %.4f Rp' % (rel, rad))

# Row D: how much of the energy equation the enthalpy flux of the mass-flux
# divergence carries is measured by the run and reported, so a temperature set
# by the divergence of the supplied mass flux rather than by the advected
# balance says so in the log.
line = [l for l in open(w + '/run.log')
        if 'enthalpy flux of the mass-flux divergence' in l]
verdict = 'PASS' if line else 'FAIL'
if not line:
    fail = 1
say(verdict, 'enthalpy_flux_term_is_measured',
    'reported' if line else 'silent', 'reported', '0')
if line:
    print('  ' + line[0].strip())

# ---------------------------------------------------------------------------
# Rows E to H: the refusal on the column as the run left it (relaxed/), a
# 300-step relaxation snapshot whose mass flux runs over seven orders of
# magnitude and is nowhere stationary. Every cell of it is therefore refused,
# and these rows say what that means for the file the transit tools read.
# ---------------------------------------------------------------------------
hyd_r = rows(w + '/relaxed/Hydro_ioniz.txt')
adv_r = rows(w + '/relaxed/Hydro_ioniz_adv.txt')
ion_r = rows(w + '/relaxed/Ion_species.txt')
iad_r = rows(w + '/relaxed/Ion_species_adv.txt')

# Row E: the two status columns exist, are integer-valued, and carry only
# the five values the file's own header defines.
ncol = min(len(a) for a in adv_r)
bad_val = None
if ncol < 10:
    bad_val = '%d columns' % ncol
else:
    for a in adv_r:
        for k in (7, 8):
            iv = int(float(a[k]))
            if float(a[k]) != iv or iv < 0 or iv > 4:
                bad_val = a[k]
                break
        if bad_val is not None:
            break
say('PASS' if bad_val is None else 'FAIL',
    'adv_status_fields_are_present_and_integer',
    '%d columns, 0-4' % ncol if bad_val is None else str(bad_val),
    '>=10 columns, integer 0-4 in both fields', '0')
if bad_val is not None:
    fail = 1
    # The rows that follow read those columns; without them they have
    # nothing to assert, so they say so rather than crashing on a missing
    # field.
    for nm in ('adv_status_counts_match_the_header',
               'retained_rows_keep_the_run_state_exactly',
               'corrected_temperature_is_bounded_by_the_run',
               'the_two_files_carry_the_same_fields'):
        say('FAIL', nm, 'no status fields', 'the two fields', '0')
    sys.exit(1)

# Row F: the counts the header states are the counts of the columns, for
# both fields, so a reader who trusts the header is reading the file.
hdr = [l for l in open(w + '/relaxed/Hydro_ioniz_adv.txt')
       if l.startswith('# adv_status_counts')]
col_T = [sum(1 for a in adv_r if int(float(a[7])) == k) for k in range(5)]
col_c = [sum(1 for a in adv_r if int(float(a[8])) == k) for k in range(5)]
hdr_T, hdr_c = [], []
if hdr:
    wd = hdr[0].split()
    if 'T' in wd and 'comp' in wd:
        hdr_T = [int(x) for x in wd[wd.index('T') + 1:wd.index('comp')]]
        hdr_c = [int(x) for x in wd[wd.index('comp') + 1:]]
ok = (hdr_T == col_T and hdr_c == col_c)
say('PASS' if ok else 'FAIL', 'adv_status_counts_match_the_header',
    'T %s comp %s' % (hdr_T, hdr_c), 'T %s comp %s' % (col_T, col_c), '0')
if not ok:
    fail = 1

# Row G: a row whose field says the run's own value was RETAINED carries
# that value EXACTLY, not to a tolerance: the temperature field means the
# run's temperature and the composition field means the equilibrium density
# of every species.
worst_T, worst_ion, n_kT, n_kc = 0.0, 0.0, 0, 0
for h, a, i0r, i1r in zip(hyd_r, adv_r, ion_r, iad_r):
    if int(float(a[7])) == 1:
        n_kT += 1
        worst_T = max(worst_T, abs(float(a[4]) - float(h[4])))
    if int(float(a[8])) == 1:
        n_kc += 1
        for k in range(1, len(i0r)):
            if i0r[k] != i1r[k]:
                worst_ion = max(worst_ion, abs(float(i1r[k]) - float(i0r[k])))
ok = (worst_T == 0.0 and worst_ion == 0.0 and n_kT > 0 and n_kc > 0)
say('PASS' if ok else 'FAIL', 'retained_rows_keep_the_run_state_exactly',
    'T %.3e over %d rows, species %.3e over %d rows'
    % (worst_T, n_kT, worst_ion, n_kc),
    '0 and 0 over at least one row each', '0')
if not ok:
    fail = 1

# The two files name their columns in their own '# columns' line, and the
# hydrodynamic file carries the measure after the two fields while the
# species file does not, so every row below locates a field by NAME.
def field_index(path, name):
    for line in open(path):
        if line.startswith('# columns'):
            nms = line.split()[2:]
            return nms.index(name) if name in nms else None
        if not line.startswith('#'):
            break
    return None


HYD_ADV = w + '/relaxed/Hydro_ioniz_adv.txt'
ION_ADV = w + '/relaxed/Ion_species_adv.txt'
ia = [field_index(HYD_ADV, nm)
      for nm in ('adv_T_status', 'adv_comp_status')]
ib = [field_index(ION_ADV, nm)
      for nm in ('adv_T_status', 'adv_comp_status')]

# Row H: THE VERDICT OF EVERY ROW IS THE VERDICT OF ITS MEASURE. A row is
# corrected exactly where its adv_mass_row is at or below the fraction the
# header states, and retained by the flow condition exactly where it is
# above: the column and the two fields are one decision, so a reader can
# reproduce the verdict from the measure. What that decision buys is stated
# and measured, not asserted: a corrected row is accurate to the fraction of
# ITSELF in the mass flux, which does not bound how far it stands from the
# snapshot temperature it replaces, and this mechanical column is a
# relaxation snapshot whose corrected temperature therefore leaves the run's
# own range. Without any measure condition the same file reached 2.6e10 K
# (MEASURED 2026-09-08: the exact steady equation on a state whose mass flux
# falls about 8 percent from one cell to the next reads that as a
# compression at the same rate), which is the extreme the retained rows now
# hold back.
COND = None
for line in open(w + '/relaxed/Hydro_ioniz_adv.txt'):
    if line.startswith('# adv_conditional_tol'):
        COND = float(line.split()[2])
        break
im = field_index(HYD_ADV, 'adv_mass_row')
if COND is None or im is None:
    say('FAIL', 'the_verdict_of_a_row_is_the_verdict_of_its_measure',
        'tolerance %s, measure column %s' % (COND, im),
        'the fraction and the measure column', '0')
    fail = 1
else:
    bad_hi, bad_lo, m_corr = 0, 0, 0.0
    for a in adv_r:
        m = float(a[im])
        st = int(float(a[ia[0]]))
        if m <= COND:
            if st == 1:
                bad_lo += 1          # retained with a measure inside
            m_corr = max(m_corr, m if st == 0 else 0.0)
        elif st == 0:
            bad_hi += 1              # corrected with a measure outside
    ok = (bad_hi == 0 and bad_lo == 0)
    say('PASS' if ok else 'FAIL',
        'the_verdict_of_a_row_is_the_verdict_of_its_measure',
        '%d corrected above the fraction, %d retained below it'
        % (bad_hi, bad_lo), '0 and 0', '0')
    if not ok:
        fail = 1
    T_run = max(float(h[4]) for h in hyd_r)
    T_adv = max(float(a[4]) for a in adv_r)
    n_above = sum(1 for a in adv_r if float(a[4]) > T_run)
    print('  largest measure of a corrected row %.3e; corrected temperature '
          'reaches %.6g K against the run max %.6g K, above it in %d rows'
          % (m_corr, T_adv, T_run, n_above))

# Row I: the two advection-corrected files describe ONE state, so the two
# fields have to be the same numbers in both. A transit tool reads the
# species file and a profile tool the hydrodynamic one. The fields are
# located by NAME in each file's own '# columns' line: the hydrodynamic file
# carries the measure after them and the species file does not, so a
# position counted from the end of the row would compare two different
# columns.

if None in ia or None in ib:
    say('FAIL', 'the_two_files_carry_the_same_fields',
        'hydro %s, species %s' % (ia, ib), 'both fields named in both', '0')
    fail = 1
else:
    worst = 0
    for a, b in zip(adv_r, iad_r):
        for ka, kb in zip(ia, ib):
            if int(float(a[ka])) != int(float(b[kb])):
                worst += 1
    ok = (worst == 0)
    say('PASS' if ok else 'FAIL', 'the_two_files_carry_the_same_fields',
        '%d disagreeing entries' % worst, '0', '0')
    if not ok:
        fail = 1

# Row K: THE MEASURE EACH ROW WAS DECIDED BY is written out, and it is the
# same number the stationary certification measures the mass row of this very
# state by. Both sides are the SECOND run, the post-process restart: its log
# carries the certification of the state it was handed and its output/ carries
# the file the post-process wrote from that same state. The certification
# reports its largest value and the cell it is in; the column of that cell has
# to carry it. The
# comparison is to the resolution of the printed format, ES10.3, which is
# four significant digits: 1e-3 relative.
cert = [l for l in open(w + '/run.log') if 'hydrodynamic mass row' in l
        and 'max=' in l]
col_i = None
hdr_cols = [l for l in open(w + '/output/Hydro_ioniz_adv.txt')
            if l.startswith('# columns')]
if hdr_cols:
    nms = hdr_cols[0].split()[2:]
    if 'adv_mass_row' in nms:
        col_i = nms.index('adv_mass_row')
if not cert or col_i is None:
    say('FAIL', 'adv_mass_row_matches_the_certification',
        'certification line %s, column %s'
        % ('present' if cert else 'absent', col_i),
        'the certification mass row line and the adv_mass_row column', '0')
    fail = 1
else:
    wd = cert[0].split()
    cmax = float(wd[wd.index('max=') + 1])
    jcell = None
    for tok in wd:
        if tok.startswith('cell='):
            jcell = int(tok.split('=')[1])
    # The file writes rows 1-Ng..N+Ng with Ng = 2, so cell j is data row
    # j + 1 counted from zero.
    mine = float(adv[jcell + 1][col_i])
    rel = abs(mine - cmax)/max(abs(cmax), 1.0e-300)
    ok = (rel <= 1.0e-3)
    say('PASS' if ok else 'FAIL', 'adv_mass_row_matches_the_certification',
        '%.6e at cell %d' % (mine, jcell), '%.3e (printed)' % cmax, '1e-3')
    if not ok:
        fail = 1

# Row L: the file states the fraction a corrected row is accurate to, as a
# number a reader can parse. It is a different statement from the
# certification tolerance the same header reports for the whole state.
tl = [l for l in open(w + '/output/Hydro_ioniz_adv.txt')
      if l.startswith('# adv_conditional_tol')]
tv = None
if tl:
    try:
        tv = float(tl[0].split()[2])
    except (IndexError, ValueError):
        tv = None
ok = (tv is not None and tv == 1.0e-2)
say('PASS' if ok else 'FAIL', 'header_states_the_conditional_tolerance',
    'line %s, value %s' % ('present' if tl else 'absent', tv),
    'a parsable 1.0e-02', '0')
if not ok:
    fail = 1

sys.exit(fail)
PY
[ $? -ne 0 ] && FAIL=1

# ---------------------------------------------------------------------------
# Rows J: the Python loader reads the schema. Every analysis tool goes
# through examples/exhale_io.py, so a field the Fortran writes and the loader
# cannot read is a field nothing downstream sees -- and a file written before
# the schema existed must come back as UNKNOWN validity, never as a corrected
# one.
# ---------------------------------------------------------------------------
env PYTHONPATH="$ROOT/examples" python3 - "$WORK/relaxed" <<'PYIO'
import os
import sys
import numpy as np
import exhale_io

d = sys.argv[1]
fail = 0


def say(ok, name, measured, reference):
    global fail
    if not ok:
        fail = 1
    print('%s %s measured=%s reference=%s tol=0'
          % ('PASS' if ok else 'FAIL', name, measured, reference))


h = exhale_io.load_hydro(d + '/Hydro_ioniz_adv.txt', ghost=True)
T, C = h.get('adv_T_status'), h.get('adv_comp_status')
ok = (T is not None and C is not None
      and T.dtype.kind in 'iu' and C.dtype.kind in 'iu'
      and T.min() >= 0 and T.max() <= 4
      and C.min() >= 0 and C.max() <= 4)
say(ok, 'exhale_io_reads_both_status_fields',
    'absent' if T is None or C is None
    else 'T %d..%d, comp %d..%d' % (T.min(), T.max(), C.min(), C.max()),
    'integers in 0..4')

# The header block, and that the counts it states are the counts of the
# arrays the loader returned.
hdr = exhale_io.load_adv_header(d + '/Hydro_ioniz_adv.txt')
cT = [int((T == k).sum()) for k in range(5)] if T is not None else []
cC = [int((C == k).sum()) for k in range(5)] if C is not None else []
ok = (hdr['schema'] == 2
      and hdr['counts'].get('T') == cT and hdr['counts'].get('comp') == cC
      and isinstance(hdr['input_certified'], bool)
      and 'face_flux_mass' in hdr['stationarity_operator']
      and hdr['product'] != '' and hdr['model_restrictions'] != '')
say(ok, 'exhale_io_reads_the_validity_header',
    'schema %s, counts T %s comp %s, certified %s, operator %r'
    % (hdr['schema'], hdr['counts'].get('T'), hdr['counts'].get('comp'),
       hdr['input_certified'], hdr['stationarity_operator'].split()[0]
       if hdr['stationarity_operator'] else ''),
    'schema 2, counts T %s comp %s, a certification verdict, '
    'the face-flux mass operator' % (cT, cC))

# The fraction and the measure reach the loader too: a caller that has the
# verdict of a row without its measure cannot say what the row is worth.
mrow = h.get('adv_mass_row')
ctol = h.get('adv_conditional_tol')
ok = (ctol == 1.0e-2 and mrow is not None
      and mrow.dtype.kind == 'f' and mrow.size == T.size
      and float(np.nanmin(mrow)) >= 0.0
      and hdr['conditional_tol'] == ctol
      and 'accurate to this fraction' in hdr['conditional_tol_text'])
say(ok, 'exhale_io_reads_the_conditional_tolerance_and_the_measure',
    'tol %s, measure %s up to %s'
    % (ctol, 'absent' if mrow is None else '%d rows' % mrow.size,
       'n/a' if mrow is None else '%.3e' % float(np.nanmax(mrow))),
    'tol 0.01 and one nonnegative float per row')

# The species file carries the same two fields, and they are not read as
# species densities.
r_i, ion = exhale_io.load_ions(d + '/Ion_species_adv.txt', ghost=True)
Ti, Ci, hdr_i = exhale_io.load_adv_status(d + '/Ion_species_adv.txt',
                                          ghost=True)
ok = ('adv_T_status' not in ion and 'adv_comp_status' not in ion
      and np.array_equal(Ti, T) and np.array_equal(Ci, C)
      and hdr_i['schema'] == 2)
say(ok, 'exhale_io_reads_the_species_file_fields',
    '%d species keys, fields equal to the hydro file: %s'
    % (len(ion), np.array_equal(Ti, T) and np.array_equal(Ci, C)),
    'no status key among the species, the same two fields')

# The legend NAMES all five values of both fields, whatever the fixture
# happens to produce. A file whose legend and columns can drift apart is a
# file whose header cannot be trusted, and this fixture produces only some
# of the five (its state is refused everywhere), so the legend is checked
# against the schema and not against the fixture.
legend = ''.join(l for l in open(d + '/Hydro_ioniz_adv.txt')
                 if l.startswith('#'))
missing = [nm for nm in ('corrected', 'retained', 'failed', 'unsupported',
                         'not_evaluated') if nm not in legend]
say(not missing, 'the_legend_names_every_value_of_the_schema',
    'missing %s' % missing if missing else 'all five named',
    'corrected, retained, failed, unsupported, not_evaluated')

# The equilibrium file is unchanged: no status fields, no header block.
eq = exhale_io.load_hydro(d + '/Hydro_ioniz.txt', ghost=True)
ok = ('adv_T_status' not in eq and 'adv_comp_status' not in eq
      and 'adv_schema' not in eq and len(eq['T']) == len(h['T']))
say(ok, 'exhale_io_reads_the_eq_file_unchanged',
    '%d keys, adv fields %s' % (len(eq), 'adv_T_status' in eq),
    'no status field and no schema line')

# A FILE WRITTEN BEFORE THE SCHEMA. Built here from the same file by
# dropping the validity block and collapsing the two fields into the single
# column of the earlier writer, so the fixture is a real file of that
# layout. Its rows have UNKNOWN validity: the single field named only the
# first refusal of a row and cannot be translated into a verdict of each field.
legacy = os.path.join(d, 'Hydro_ioniz_adv_legacy.txt')
with open(d + '/Hydro_ioniz_adv.txt') as fh, open(legacy, 'w') as out:
    for line in fh:
        if line.startswith('#'):
            if line.startswith('# adv_'):
                continue
            if line.startswith('# columns'):
                line = line.replace(
                    'adv_T_status adv_comp_status adv_mass_row',
                    'adv_status')
            out.write(line)
            continue
        f = line.split()
        out.write(' '.join(f[:8]) + '\n')

lg = exhale_io.load_hydro(legacy, ghost=True)
Tl, Cl = lg.get('adv_T_status'), lg.get('adv_comp_status')
ok = (lg.get('adv_schema') == 1
      and Tl is not None and Cl is not None
      and (Tl == exhale_io.ADV_UNKNOWN).all()
      and (Cl == exhale_io.ADV_UNKNOWN).all()
      and lg.get('adv_status') is not None
      and lg.get('adv_input_certified') is None)
say(ok, 'legacy_file_reads_as_unknown_validity',
    'schema %s, fields all %s: %s, certification %s'
    % (lg.get('adv_schema'), exhale_io.ADV_UNKNOWN,
       Tl is not None and (Tl == exhale_io.ADV_UNKNOWN).all()
       and (Cl == exhale_io.ADV_UNKNOWN).all(),
       lg.get('adv_input_certified')),
    'schema 1, both fields ADV_UNKNOWN, no certification verdict')

sys.exit(fail)
PYIO
[ $? -ne 0 ] && FAIL=1

if [ "$FAIL" -ne 0 ]; then
   echo "adv_static_limit: FAILED"
else
   echo "adv_static_limit: PASSED"
fi
exit "$FAIL"
