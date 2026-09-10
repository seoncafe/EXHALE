#!/usr/bin/env python3
"""The convergence functional the run used, recomputed from the state it wrote.

QUANTITY UNDER TEST
    The marching loop stops on the radial spread of the mass flux rho v r^2
    over the window r >= "Escape radius [R_p]" (cells j_min..N):

        spread = ( max_j F_j - min_j F_j ) / | mean_j F_j | ,
        F_j    = rho_j v_j r_j^2 .

    A steady spherical wind carries one mass flux at every radius, so this
    is the statement that the marching has reached one (ATES, Caldiroli et
    al. 2021, stops on the constancy of rho v r^2; CETIMB, Koskinen et al.
    2013a, asks that F_c = rho v r^2 be constant with altitude).  The sign is
    kept: a window in which the flux reverses is not steady, and taking |F|
    before the extrema would read it as flat.  The normalization is the mean
    of the window, not its minimum.

    THIS TEST NO LONGER TRANSCRIBES THE PRODUCTION EXPRESSION.  The
    functional is `mass_flux_spread` in src/EXHALE_main.f90 (a contained
    function of the main program, so it cannot be linked against from here),
    and the run prints the value it used at its last step in its summary:

        mass_flux_spread: window cells <j_min>..<N>  value=<spread>

    That line is what this test reads, together with the window it names, and
    the assertion is that the same functional evaluated on the profile the
    same run wrote gives the same number.  One functional, one number, and
    the number in the run summary is a property of the file the user reads.

WHAT IS RUN
    A copy of backup/regression/mol_base_handoff (its input.inp and base.inp)
    in build/tests/grid_and_gates/flux_spread, capped at 200 steps.  The cap
    matters: the run must end BY MARCHING, because the summary line reports
    the last marching step and a JFNK finish would move the state after it.
    The regression directory is never written to.

ASSERTIONS
    1. mass_flux_spread_recomputed_from_output: the functional evaluated on
       output/Hydro_ioniz.txt over the window the summary names, against the
       value the summary reports.  Tolerance 1e-9 relative: the profile
       columns are written at full double precision and the functional is
       scale free, so the two differ only by the decimal round trip.
    2. mass_flux_spread_of_sign_reversing_window: the definition above,
       evaluated on a synthetic window carrying -F over its inner half and
       +3F over its outer half, is (3F - (-F))/|(-F + 3F)/2| = 4.  This one
       is a property of the definition, transcribed here from the
       specification and not from the code: it pins what the measure must
       say about a window whose flux reverses, which is the state the stop
       must never accept as flat.  Reading the magnitudes first and
       normalizing by the smallest of them -- which is what the stop did --
       gives (3F - F)/F = 2 on the same window: half as far from flat, on a
       window carrying inflow at one radius and outflow at another.

EXPECTED BEFORE section 10.2 item 5: RED -- the run summary carries no such
line, because the stop and the reported flux spread were two different
functionals of the same window and neither was named.
"""

import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
EXE = os.environ.get("EXHALE_EXE", os.path.join(ROOT, "EXHALE.x"))
CASE = os.path.join(ROOT, "backup", "regression", "mol_base_handoff")
WORK = os.path.join(ROOT, "build", "tests", "grid_and_gates", "flux_spread")
MAXSTEPS = "200"

TINY = 1.0e-30
n_fail = 0


def verdict(name, measured, reference, tol):
    global n_fail
    ok = abs(measured - reference) <= tol * max(abs(reference), TINY)
    print("%s %s measured=%.12e reference=%.12e tol=%.2e"
          % ("PASS" if ok else "FAIL", name, measured, reference, tol))
    if not ok:
        n_fail += 1


def spread(mom):
    """The functional under test, on an array of rho v r^2."""
    return ((np.max(mom) - np.min(mom))
            / max(abs(np.mean(mom)), TINY))


def read_hydro(path):
    """r, rho, v and the grid size, from the file's own header."""
    n_cells = None
    n_ghost = None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = re.search(r"rows\s+\d+:\s+(\d+)\s+ghost cells", line)
            if m:
                n_ghost = int(m.group(1))
            m = re.search(r"\bN=(\d+)", line)
            if m:
                n_cells = int(m.group(1))
    if n_cells is None or n_ghost is None:
        raise RuntimeError("%s: no N / ghost count in the header" % path)
    data = np.loadtxt(path, comments="#")
    if data.shape[0] != n_cells + 2 * n_ghost:
        raise RuntimeError("%s: %d rows for N=%d, Ng=%d"
                           % (path, data.shape[0], n_cells, n_ghost))
    return data[:, 0], data[:, 1], data[:, 2], n_cells, n_ghost


# ---- the run -------------------------------------------------------------
if not os.access(EXE, os.X_OK):
    print("FAIL mass_flux_spread_recomputed_from_output measured=no_binary "
          "reference=%s tol=0" % EXE)
    sys.exit(1)

shutil.rmtree(WORK, ignore_errors=True)
os.makedirs(os.path.join(WORK, "output"))
for f in ("input.inp", "base.inp"):
    shutil.copy(os.path.join(CASE, f), WORK)
env = dict(os.environ, OMP_NUM_THREADS="1", EXHALE_MAXSTEPS=MAXSTEPS)
with open(os.path.join(WORK, "run.log"), "w") as log:
    rc = subprocess.call([EXE], cwd=WORK, stdout=log, stderr=subprocess.STDOUT,
                         env=env)
if rc != 0:
    print("FAIL mass_flux_spread_run measured=exit_%d reference=exit_0 tol=0"
          % rc)
    print("     see %s/run.log" % WORK)
    sys.exit(1)
print("PASS mass_flux_spread_run measured=exit_0 reference=exit_0 tol=0")

# ---- assertion 1: the run's own number, from the run's own state ---------
dump = None
with open(os.path.join(WORK, "run.log")) as fh:
    for line in fh:
        m = re.search(r"mass_flux_spread:\s*window cells\s*(\d+)\.\.(\d+)"
                      r"\s*value=\s*(\S+)", line)
        if m:
            dump = (int(m.group(1)), int(m.group(2)), float(m.group(3)))
if dump is None:
    print("FAIL mass_flux_spread_recomputed_from_output measured=no_dump_line "
          "reference=one_line tol=0")
    print("     the run summary carries no 'mass_flux_spread:' line, so the "
          "value the stop used cannot be read from the run")
    sys.exit(1)

j_min, n_last, value = dump
r_all, rho_all, v_all, n_cells, n_ghost = read_hydro(
    os.path.join(WORK, "output", "Hydro_ioniz.txt"))
# Row k (0-based) of the file is cell j = k + 1 - n_ghost.
sel = slice(j_min + n_ghost - 1, n_last + n_ghost)
mom = rho_all[sel] * v_all[sel] * r_all[sel] ** 2
n_sign_change = int(np.sum(np.sign(mom[1:]) * np.sign(mom[:-1]) < 0))
print("  window cells %d..%d (%d cells, r = %.5f to %.5f), sign changes in "
      "it: %d" % (j_min, n_last, len(mom), r_all[sel][0], r_all[sel][-1],
                  n_sign_change))
verdict("mass_flux_spread_recomputed_from_output", spread(mom), value, 1.0e-9)

# ---- assertion 2: what the definition says about a reversing window ------
m = 50
flux = 3.7e11
synthetic = np.concatenate([np.full(m, -flux), np.full(m, 3.0 * flux)])
print("  synthetic window: %d cells at %.3e and %d at %.3e"
      % (m, -flux, m, 3.0 * flux))
print("     the magnitude-first, minimum-normalized form reads %.6f on it"
      % ((np.max(np.abs(synthetic)) - np.min(np.abs(synthetic)))
         / np.min(np.abs(synthetic))))
verdict("mass_flux_spread_of_sign_reversing_window", spread(synthetic), 4.0,
        1.0e-12)

if n_fail:
    print("mass_flux_spread_functional: %d assertion(s) FAILED" % n_fail)
    sys.exit(1)
print("mass_flux_spread_functional: all assertions PASSED")
