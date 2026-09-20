#!/bin/bash
# The run-mode contract (docs/a0_run_mode_contract_20260906.md section 7),
# tested on the binary.
#
# WHAT IS UNDER TEST. A run states which of three things it is doing, and the
# state it writes carries that statement:
#   init  initialization / continuation. Local pseudo-time, PTC and a
#         stationary Newton finish are permitted; no elapsed time is claimed.
#   phys  physical integration. One global dt per step, a handoff test on the
#         state the first step is taken from, and a clock that advances only
#         after a complete accepted step.
# The assertions below are the seven rows of the contract's section 7 table.
#
# Every test prints one
#     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
# line per assertion; lines beginning with two spaces are context. Exit status
# is nonzero if any assertion fails.
#
# Usage: src/tests/run_mode/run.sh
#   EXHALE_EXE           binary to test (default $ROOT/EXHALE.x)
#   EXHALE_TEST_OUT      working directory (default $ROOT/build/tests/run_mode)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
# The binary is selected and its identity stated in one place;
# src/tests/exhale_exe.sh carries the policy.
. "$HERE/../exhale_exe.sh"
exhale_select_exe "$ROOT" run_mode
EXE="$EXHALE_RUN_EXE"
exhale_announce_exe
WORK="${EXHALE_TEST_OUT:-$ROOT/build/tests/run_mode}"
REG="$ROOT/backup/regression"
FAIL=0

say () { echo "$@"; }
chk () { # name measured reference
   if [ "$2" = "$3" ]; then say "PASS $1 measured=$2 reference=$3 tol=0"
   else say "FAIL $1 measured=$2 reference=$3 tol=0"; FAIL=1; fi
}

if [ ! -x "$EXE" ]; then
   echo "FAIL run_mode_binary measured=no_binary reference=$EXE tol=0"; exit 1
fi

mkdir -p "$WORK"
mkcase () { # dir  regression-case  [extra lines...]
   local d="$WORK/$1"; shift
   local c="$1"; shift
   rm -rf "$d"; mkdir -p "$d/output"
   local f
   for f in input.inp base.inp metals.inp opacity.inp; do
      [ -f "$REG/$c/$f" ] && cp "$REG/$c/$f" "$d/"
   done
   for f in "$@"; do echo "$f" >> "$d/input.inp"; done
}

# ---------------------------------------------------------------------------
# 1. "Run mode: phys" with "Time stepping: Local" is refused at startup.
#    Local pseudo-time advances neighbouring cells by different intervals, so
#    the material sum across an internal face is not conserved: the update is
#    a relaxation iterate and has no elapsed time.
# ---------------------------------------------------------------------------
mkcase r1 wasp_localdt "Run mode: phys"
( cd "$WORK/r1" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=1 "$EXE" > run.log 2>&1 )
rc=$?
[ $rc -ne 0 ] && v=refused || v=accepted
chk phys_with_local_time_stepping_is_refused "$v" refused
n=$(grep -c 'Time stepping: Local" cannot both be set' "$WORK/r1/run.log")
chk the_refusal_names_both_keys "$n" 1

# ---------------------------------------------------------------------------
# 2. A physical run: the handoff is tested, the header carries the mode and
#    the clock, and the elapsed time equals the sum of the accepted global dt.
#    The sum is MEASURED from the step clock, not asserted from a formula.
# ---------------------------------------------------------------------------
mkcase r2 hydrostatic_column "Run mode: phys"
( cd "$WORK/r2" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=20 EXHALE_STEP_CLOCK=1 \
     "$EXE" > run.log 2>&1 )
[ $? -eq 0 ] && v=ran || v=failed
chk a_physical_run_completes "$v" ran
n=$(grep -c 'handoff accepted' "$WORK/r2/run.log")
chk the_handoff_is_tested_before_the_first_step "$n" 1
n=$(grep -m1 '^# coupling' "$WORK/r2/output/Hydro_ioniz.txt" | grep -c 'mode=phys t_phys=')
chk the_header_carries_the_mode_and_the_clock "$n" 1
v=$(python3 - "$WORK/r2" <<'PY'
import re,sys
log=open(sys.argv[1]+'/run.log').read()
dt=[float(x) for x in re.findall(r'dt=\s*(\S+)\s+t_phys',log)]
hdr=open(sys.argv[1]+'/output/Hydro_ioniz.txt').readline()
for l in open(sys.argv[1]+'/output/Hydro_ioniz.txt'):
    if l.startswith('# coupling'): hdr=l; break
t=float(re.search(r't_phys=(\S+)',hdr).group(1))
s=sum(dt)
print('equal' if (t!=0.0 and abs(s-t)<=1e-14*abs(t)) else 'differs')
PY
)
chk t_phys_is_the_sum_of_the_accepted_dt "$v" equal

# ---------------------------------------------------------------------------
# 3. A rejected trial advances no clock. The rejection is taken through the
#    positivity retry path itself (EXHALE_REJECT_STEP), which discards the
#    attempt and retakes the step at half dt: one more attempt, no more
#    accepted steps, and the clock still the sum of what was accepted.
# ---------------------------------------------------------------------------
mkcase r3 hydrostatic_column "Run mode: phys"
( cd "$WORK/r3" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=20 EXHALE_STEP_CLOCK=1 \
     EXHALE_REJECT_STEP=3 "$EXE" > run.log 2>&1 )
v=$(grep -o 'steps: [0-9]* accepted of [0-9]* attempted' "$WORK/r3/run.log" \
    | awk '{print $2"_"$5}')
chk a_rejected_trial_is_an_attempt_and_not_a_step "$v" 20_21
v=$(python3 - "$WORK/r3" <<'PY'
import re,sys
log=open(sys.argv[1]+'/run.log').read()
dt=[float(x) for x in re.findall(r'dt=\s*(\S+)\s+t_phys',log)]
tp=[float(x) for x in re.findall(r't_phys=\s*(\S+)',log)]
print('equal' if abs(sum(dt)-tp[-1])<=1e-14*abs(tp[-1]) and len(dt)==20 else 'differs')
PY
)
chk the_rejected_trial_advanced_no_time "$v" equal

# ---------------------------------------------------------------------------
# 4. An initialization snapshot restarted as a physical continuation is
#    refused: its steps were relaxation iterates, it carries no elapsed time,
#    and no clock can be invented for it.
# ---------------------------------------------------------------------------
mkcase r4a hydrostatic_column
( cd "$WORK/r4a" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=10 "$EXE" > run.log 2>&1 )
n=$(grep -c 'Run mode: init (default)' "$WORK/r4a/EXHALE_setup.out")
chk a_run_that_says_nothing_is_an_initialization_run "$n" 1
n=$(grep -m1 '^# coupling' "$WORK/r4a/output/Hydro_ioniz.txt" | grep -c 'mode=init')
chk an_initialization_run_writes_mode_init "$n" 1
n=$(grep -m1 '^# coupling' "$WORK/r4a/output/Hydro_ioniz.txt" | grep -c 't_phys=')
chk an_initialization_state_carries_no_clock "$n" 0
# The handoff: that same snapshot loaded by a physical run is a new
# trajectory. It is tested and its clock starts at zero; no time is invented
# and none is taken from a file that has none.
mkcase r4b hydrostatic_column "Run mode: phys"
sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/r4b/input.inp"
cp "$WORK/r4a/output/Hydro_ioniz.txt" "$WORK/r4b/output/Hydro_ioniz_IC.txt"
cp "$WORK/r4a/output/Ion_species.txt" "$WORK/r4b/output/Ion_species_IC.txt"
( cd "$WORK/r4b" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=2 EXHALE_STEP_CLOCK=1 \
     "$EXE" > run.log 2>&1 )
[ $? -eq 0 ] && v=ran || v=failed
chk an_init_snapshot_loaded_by_a_phys_run_is_the_handoff "$v" ran
n=$(grep -c 'physical time origin is t = 0' "$WORK/r4b/run.log")
chk the_handoff_declares_the_time_origin "$n" 1
v=$(python3 - "$WORK/r4b" <<'PY'
import re,sys
log=open(sys.argv[1]+'/run.log').read()
dt=[float(x) for x in re.findall(r'dt=\s*(\S+)\s+t_phys',log)]
tp=[float(x) for x in re.findall(r'step-clock:.*t_phys=\s*(\S+)',log)]
# The clock started at zero iff the time after the first accepted step is
# that step's own dt, and the last is the sum of all of them.
print('from_zero' if len(dt)>=1 and abs(tp[0]-dt[0])<=1e-14*dt[0]
      and abs(tp[-1]-sum(dt))<=1e-14*sum(dt) else 'elsewhere')
PY
)
chk the_handoff_clock_starts_at_zero "$v" from_zero
# What IS refused: a header that claims a trajectory and carries no time on it.
mkcase r4c hydrostatic_column "Run mode: phys"
sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/r4c/input.inp"
# The clock is written TWICE, as a ' t_phys=' token of the '# coupling:' line
# and as the '# t_phys[s]' line of the restart metadata block, and the loader
# reads the second where the first is absent. Both go, and the metadata line
# goes from BOTH halves of the pair: the two files are two halves of one state
# and a block present in one and absent in the other is refused on that ground
# instead, which would make this assertion pass while testing something else.
# (It used to strip the coupling token alone, so the clock reached the binary
# through the metadata block and the state was accepted.)
sed -e 's/ t_phys=[^ ]*//' -e '/^# t_phys\[s\]/d' \
    "$WORK/r2/output/Hydro_ioniz.txt" > "$WORK/r4c/output/Hydro_ioniz_IC.txt"
sed '/^# t_phys\[s\]/d' "$WORK/r2/output/Ion_species.txt" \
    > "$WORK/r4c/output/Ion_species_IC.txt"
( cd "$WORK/r4c" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=2 "$EXE" > run.log 2>&1 )
[ $? -ne 0 ] && v=refused || v=accepted
chk a_phys_header_without_its_clock_is_refused "$v" refused
n=$(grep -c 'its t_phys field is missing' "$WORK/r4c/run.log")
chk a_phys_header_without_its_clock_names_the_clock "$n" 1

# ---------------------------------------------------------------------------
# 5. A physical state restarted as a physical run continues its clock.
# ---------------------------------------------------------------------------
mkcase r5 hydrostatic_column "Run mode: phys"
sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/r5/input.inp"
cp "$WORK/r2/output/Hydro_ioniz.txt" "$WORK/r5/output/Hydro_ioniz_IC.txt"
cp "$WORK/r5/../r2/output/Ion_species.txt" "$WORK/r5/output/Ion_species_IC.txt"
( cd "$WORK/r5" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=5 "$EXE" > run.log 2>&1 )
[ $? -eq 0 ] && v=ran || v=failed
chk a_physical_state_restarts "$v" ran
v=$(python3 - "$WORK/r2" "$WORK/r5" <<'PY'
import re,sys
def tp(d):
    for l in open(d+'/output/Hydro_ioniz.txt'):
        if l.startswith('# coupling'):
            m=re.search(r't_phys=(\S+)',l)
            return float(m.group(1)) if m else None
    return None
a,b=tp(sys.argv[1]),tp(sys.argv[2])
log=open(sys.argv[2]+'/run.log').read()
m=re.search(r'continues from t_phys =\s*(\S+)',log)
print('continues' if (m and a is not None and b is not None and b>a
                      and abs(float(m.group(1))-a)<=1e-6*a) else 'restarts_at_zero')
PY
)
chk the_clock_continues_from_the_header "$v" continues

# ---------------------------------------------------------------------------
# 6. A handoff state the equations do not describe is refused, and the cell is
#    named. The state is built by putting one cell of a written physical state
#    outside the admissible set (negative pressure), which is what the three
#    handoff conditions exist to catch: the cell is inadmissible, the
#    conserved state has no positive internal energy there, and the sweep
#    accepts that cell under the relaxation amnesty (class 4).
# ---------------------------------------------------------------------------
mkcase r6 hydrostatic_column "Run mode: phys"
sed -i 's/^Load IC?.*/Load IC? True/' "$WORK/r6/input.inp"
python3 - "$WORK/r2/output/Hydro_ioniz.txt" "$WORK/r6/output/Hydro_ioniz_IC.txt" <<'PY'
import sys
src,dst=sys.argv[1],sys.argv[2]
out=[];n=0
for line in open(src):
    if line.lstrip().startswith('#'): out.append(line); continue
    n+=1
    if n==42:                      # two lower ghosts, then physical cell 40
        f=line.split()
        f[3]='-'+f[3].lstrip('-'); f[4]='-'+f[4].lstrip('-')
        line=' '+'   '.join(f)+'\n'
    out.append(line)
open(dst,'w').writelines(out)
PY
cp "$WORK/r2/output/Ion_species.txt" "$WORK/r6/output/Ion_species_IC.txt"
( cd "$WORK/r6" && OMP_NUM_THREADS=1 EXHALE_RELOAD_EQ=0 EXHALE_MAXSTEPS=1 \
     "$EXE" > run.log 2>&1 )
[ $? -ne 0 ] && v=refused || v=accepted
chk an_inadmissible_handoff_state_is_refused "$v" refused
n=$(grep -c 'REFUSED: 1 cell(s) are not admissible; first is cell 40' "$WORK/r6/run.log")
chk the_refusal_names_the_inadmissible_cell "$n" 1
n=$(grep -c 'first such cell: 40' "$WORK/r6/run.log")
chk the_refusal_names_the_cell_without_a_chemical_root "$n" 1

# ---------------------------------------------------------------------------
# 7. init is the default in every configuration, and stating it changes
#    nothing: the numeric content of an initialization run does not depend on
#    whether the key is in the file. That is what keeps every existing input,
#    the regression matrix included, the run it was.
# ---------------------------------------------------------------------------
mkcase r7a hydrostatic_column "Run mode: init"
( cd "$WORK/r7a" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=10 "$EXE" > run.log 2>&1 )
n=$(grep -c 'Run mode: init (given)' "$WORK/r7a/EXHALE_setup.out")
chk a_stated_init_run_says_it_was_stated "$n" 1
if diff -q <(grep -v '^ *#' "$WORK/r7a/output/Hydro_ioniz.txt") \
           <(grep -v '^ *#' "$WORK/r4a/output/Hydro_ioniz.txt") > /dev/null
then v=identical; else v=differs; fi
chk the_stated_and_the_defaulted_init_run_agree "$v" identical

# ---------------------------------------------------------------------------
# 8. THE EXIT STATUS FOLLOWS THE STATIONARITY CLAIM, NOT THE RUN MODE.
#    The `init` exemption covers the PHYSICAL-step claims (the clock, the
#    budgets, the histories), which an initialization run does not make.
#    Stationarity is not one of them: a run stopped by the du criterion
#    asserts that the state it writes is stationary, in whichever mode it
#    made that assertion, and the status is what tells a caller whether the
#    certification accepted it.
#    8a is the same case ended by the step cap: a cap declares nothing, so
#    the status is 0 however the certification came out.
#    8b reloads 8a's snapshot with the du stop armed by the reload and a du
#    threshold above the state's own du, so the stop fires at the first
#    evaluation and the run stops on the du criterion. The state is a
#    10-step relaxation snapshot of a hydrostatic column and the
#    certification refuses it. The loose threshold is what makes the claim
#    cheap to reach; the row is about the status of a declared claim and not
#    about the quality of the state, and the state is written in full either
#    way.
# ---------------------------------------------------------------------------
mkcase r8a hydrostatic_column
( cd "$WORK/r8a" && OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=10 "$EXE" > run.log 2>&1 )
rc=$?
chk a_cap_declares_nothing_and_exits_0 "$rc" 0
n=$(grep -c 'NOT CERTIFIED' "$WORK/r8a/run.log")
[ "$n" -ge 1 ] && v=refused || v=certified
chk the_capped_state_is_itself_uncertified "$v" refused

mkcase r8b hydrostatic_column
sed -i -e 's/^Load IC?.*/Load IC? True/' \
       -e 's/^du_th .*/du_th [PLM,WENO3]: 1.0e3 1.0e3/' "$WORK/r8b/input.inp"
cp "$WORK/r8a/output/Hydro_ioniz.txt" "$WORK/r8b/output/Hydro_ioniz_IC.txt"
cp "$WORK/r8a/output/Ion_species.txt" "$WORK/r8b/output/Ion_species_IC.txt"
( cd "$WORK/r8b" && OMP_NUM_THREADS=1 "$EXE" > run.log 2>&1 )
rc=$?
n=$(grep -c 'stopped: mass flux' "$WORK/r8b/run.log")
chk the_du_stop_declares_a_stationary_state "$n" 1
n=$(grep -c 'Run mode: init (default)' "$WORK/r8b/EXHALE_setup.out")
chk the_declaring_run_is_an_initialization_run "$n" 1
n=$(grep -c 'NOT CERTIFIED' "$WORK/r8b/run.log")
[ "$n" -ge 1 ] && v=refused || v=certified
chk the_declared_state_is_refused_by_the_certification "$v" refused
chk a_refused_stationary_claim_exits_2_in_init_mode "$rc" 2
n=$(grep -c 'exits with status 2 because it is not certified' "$WORK/r8b/run.log")
chk the_status_names_its_reason "$n" 1
# The status is a verdict about a complete set of files, never a reason for
# an incomplete one: the profiles and the advection-corrected post-process
# are on disk before the stop.
v=missing
[ -s "$WORK/r8b/output/Hydro_ioniz.txt" ] && [ -s "$WORK/r8b/output/Hydro_ioniz_adv.txt" ] \
   && [ -s "$WORK/r8b/output/Ion_species.txt" ] && v=written
chk the_refused_run_wrote_its_outputs_in_full "$v" written

exit $FAIL
