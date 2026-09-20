#!/usr/bin/env python3
"""Focused, non-mutating production-code review checks for PLAN_20260919.

Only this audit directory and fresh temporary directories receive products.
No solve, catalog publication, source change, or golden refresh is performed.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
AUDIT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def function_text(source, name):
    match = re.search(r"^" + re.escape(name) + r"\s*\(\)\s*\{.*?^\}",
                      source, re.M | re.S)
    if not match:
        raise ValueError(name)
    return match.group(0)


def inventory():
    models = ROOT / "LHS1140b/models"
    indices = sorted(models.rglob("state_index.json"))
    records = [(p, json.loads(p.read_text())) for p in indices]
    print("INVENTORY", json.dumps({
        "index_count": len(records),
        "latest_certified_count": sum(bool(r.get("latest_certified")) for _, r in records),
        "without_latest_certified": [str(p.parent.relative_to(models))
                                     for p, r in records if not r.get("latest_certified")],
        "binary_md5": digest(ROOT / "EXHALE.x"),
        "make_q_exit": subprocess.run(["make", "-q"], cwd=ROOT, capture_output=True).returncode,
    }, sort_keys=True))
    paths = ["docs/PLAN_20260919.md", "src/EXHALE_main.f90",
             "src/modules/radiation/ionization_equilibrium.f90",
             "src/modules/states/base_boundary.f90",
             "src/modules/time_step/certification.f90",
             "src/modules/time_step/eval_dt.f90",
             "src/modules/time_step/steady_newton.f90",
             "src/modules/radiation/charge_exchange.f90",
             "src/tests/grid_and_gates/run.sh",
             "LHS1140b/models/run_case.sh", "exhale_transit_lib.py",
             "EXHALE_transit.py"]
    print("SOURCE_MD5", json.dumps({p: digest(ROOT / p) for p in paths}, sort_keys=True))
    sys.path.insert(0, str(models))
    import status as catalog_status
    chosen = []
    for name, ladder in catalog_status.GROUPS:
        for heh in ladder:
            case = catalog_status.Case(catalog_status.Group(name), heh)
            if case.product_dir:
                index = Path(case.product_dir) / "state_index.json"
                if index.exists():
                    chosen.append((index, json.loads(index.read_text())))
    print("CATALOG_GROUP_SELECTION", json.dumps({
        "indexed": len(chosen),
        "latest_certified": sum(bool(data.get("latest_certified")) for _, data in chosen),
    }))


def workflow_checks(scratch):
    source = (ROOT / "LHS1140b/models/run_case.sh").read_text()
    funcs = "\n".join(function_text(source, n) for n in
                      ("certification_block", "state_is_nonfinite", "classify_ending"))
    for name, logtext in [
        ("inner_success_without_outer_verdict", " (JFNK) done info=0 ||R||= 1.0E-12\n"),
        ("outer_success_without_state_files", " the stationary solve returned info = 0\n"),
        ("inner_failure_without_outer_verdict", " (JFNK) done info=2 ||R||= 1.0E-1\n"),
    ]:
        log = scratch / (name + ".log")
        log.write_text(logtext)
        proc = subprocess.run(["bash", "-s", "--", str(log), str(scratch / "absent")],
                              input=funcs + '\nclassify_ending "$1" "$2"\n'
                              + 'printf "%s|%s\\n" "$INFO" "$ENDING_CLASS"\n',
                              text=True, capture_output=True, check=True)
        print("CLASSIFICATION", name, proc.stdout.strip())

    tree = ast.parse((ROOT / "exhale_transit_lib.py").read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                and n.name == "transit_state_files")
    namespace = {"os": os}
    exec(compile(ast.Module(body=[node], type_ignores=[]), "transit_state_files", "exec"), namespace)
    print("TRANSIT_SOLUTION", namespace["transit_state_files"]("fixture", "solution"))
    driver = (ROOT / "src/tests/grid_and_gates/run.sh").read_text()
    invocations = re.findall(r"run_one [^\n]+(?:\\\n[^\n]*)*?bash \"\$HERE/([^\"]+\.sh)\"", driver)
    print("HARNESS_SHELL_TARGETS", len(invocations))
    for script in invocations:
        body = (ROOT / "src/tests/grid_and_gates" / script).read_text()
        if "EXHALE_TEST_OUT" not in body:
            print("HARNESS_CHILD_WITHOUT_OUT", script)
    env = dict(os.environ, EXHALE_TEST_OUT=str(scratch / "requested"))
    inherited = subprocess.run(["bash", "-c", 'env EXHALE_EXE=unused bash -c \'printf "%s" "$EXHALE_TEST_OUT"\''],
                               env=env, text=True, capture_output=True, check=True)
    print("HARNESS_INHERITED_OUT", inherited.stdout == env["EXHALE_TEST_OUT"])


def evaluate(scratch):
    case = ROOT / "LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH9.7"
    generation = case / "states/g0004_20260919T103629Z_9b8394a3"
    stage = scratch / "molecular_evaluate"
    (stage / "output").mkdir(parents=True)
    for name in ("base.inp", "metals.inp", "opacity.inp"):
        if (case / name).exists():
            shutil.copyfile(case / name, stage / name)
    lines = []
    for line in (case / "input.inp").read_text().splitlines():
        if line.startswith("Restart intent:"):
            continue
        key, sep, value = line.partition(":")
        value = value.strip()
        if sep and value and " " not in value and not os.path.isabs(value) and (case / value).exists():
            line = key + ": " + str((case / value).resolve())
        if line.startswith("Load IC?"):
            line = "Load IC? True"
        lines.append(line)
    lines.append("Restart intent: stationary evaluate")
    (stage / "input.inp").write_text("\n".join(lines) + "\n")
    for stem in ("Hydro_ioniz", "Ion_species"):
        shutil.copyfile(generation / (stem + ".txt"), stage / "output" / (stem + "_IC.txt"))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               EXHALE_BOUNDARY_TRACE=str(stage / "boundary_trace.txt"))
    logfile = AUDIT / "plan_20260919_evaluate.log"
    with logfile.open("w") as log:
        proc = subprocess.run(["timeout", "150", str(ROOT / "EXHALE.x")],
                              cwd=stage, env=env, stdout=log, stderr=subprocess.STDOUT)
    print("EVALUATE", json.dumps({"exit": proc.returncode, "scratch": str(stage), "log": str(logfile)}))
    for line in logfile.read_text().splitlines():
        if any(s in line for s in ("hydrodynamic mass row", "ghost H2 partition", "solved ghost counts",
                                   "CERTIFIED:", "NOT CERTIFIED:", "work state verdict")):
            print("EVALUATE_ROW", line.strip())
    trace = stage / "output/boundary_trace.txt"
    if trace.exists():
        shutil.copyfile(trace, AUDIT / "plan_20260919_boundary_trace.txt")


def capacity(scratch):
    compiler = Path(shutil.which("gfortran")).resolve()
    objects = [str(p) for p in sorted((ROOT / "build").glob("*.o"))
               if not re.search(r"(EXHALE_main|_tests|_probe)\.o$", p.name)]
    library = compiler.parent.parent / "lib"
    link = (["-L" + str(library), "-lopenblas", "-Wl,-rpath," + str(library), "-ldl"]
            if (library / "libopenblas.so").exists() else ["-llapack", "-ldl"])
    executable = scratch / "cert_capacity.x"
    command = [str(compiler), "-O0", "-g", "-fopenmp", "-I" + str(ROOT / "build"),
               "-J" + str(scratch), "-o", str(executable),
               str(AUDIT / "plan_20260919_cert_capacity.f90")] + objects + link
    subprocess.run(command, check=True, cwd=scratch, capture_output=True)
    proc = subprocess.run([str(executable)], check=True, text=True,
                          capture_output=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
    print("CAPACITY_PRODUCTION_LINK", proc.stdout)


def analyze_evaluation(stage):
    import numpy as np
    original = ROOT / ("LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH9.7/states/"
                       "g0004_20260919T103629Z_9b8394a3")
    h0 = np.loadtxt(original / "Hydro_ioniz.txt")
    h1 = np.loadtxt(stage / "output/Hydro_ioniz.txt")
    i0 = np.loadtxt(original / "Ion_species.txt")
    i1 = np.loadtxt(stage / "output/Ion_species.txt")
    p = slice(2, -2)
    for name, col in [("rho", 1), ("velocity", 2), ("pressure", 3), ("temperature", 4)]:
        difference = np.max(np.abs(h1[p, col] - h0[p, col]) / np.maximum(np.abs(h0[p, col]), 1e-99))
        print("PHYSICAL_ROUND_TRIP", name, format(difference, ".16e"))
    for row in [0, 1]:
        nh = i1[row, 1:6].sum() + i1[row, 34:38].sum()
        ne_atomic = i1[row, 2] + i1[row, 4] + 2 * i1[row, 5]
        ne_molecular = i1[row, 35:38].sum()
        nall = nh + ne_atomic + ne_molecular
        print("GHOST_THERMODYNAMICS", json.dumps({
            "cell": row - 1, "r": float(h1[row, 0]), "T_K": float(h1[row, 4]),
            "p_cgs": float(h1[row, 3]), "particles_per_mass": float(nall / h1[row, 1]),
            "EOS_p_over_nkT_minus_one": float(h1[row, 3] / (nall * 1.380649e-16 * h1[row, 4]) - 1),
            "molecular_fraction_of_electrons": float(ne_molecular / (ne_atomic + ne_molecular)),
            "T_relative_change_from_saved_ghost": float(h1[row, 4] / h0[row, 4] - 1),
        }))
    trace = stage / "output/boundary_trace.txt"
    if trace.exists():
        shutil.copyfile(trace, AUDIT / "plan_20260919_boundary_trace.txt")


def jacobian(scratch):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               EXHALE_EXE=str(ROOT / "EXHALE.x"), EXHALE_TEST_OUT=str(scratch / "jacobian"),
               EXHALE_JAC_TIMEOUT="90")
    logpath = AUDIT / "plan_20260919_coupled_jacobian.log"
    with logpath.open("w") as log:
        proc = subprocess.run(["timeout", "200", "bash", "src/tests/coupled_block_jacobian/run.sh"],
                              cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    print("JACOBIAN", json.dumps({"exit": proc.returncode, "scratch": str(scratch / "jacobian")}))
    print(logpath.read_text())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluate", action="store_true")
    parser.add_argument("--jacobian", action="store_true")
    parser.add_argument("--capacity", action="store_true")
    parser.add_argument("--analyze-evaluation", type=Path)
    args = parser.parse_args()
    scratch = Path(tempfile.mkdtemp(prefix="exhale_plan_20260919_review_"))
    print("SCRATCH", scratch)
    inventory()
    workflow_checks(scratch)
    if args.evaluate:
        evaluate(scratch)
    if args.jacobian:
        jacobian(scratch)
    if args.capacity:
        capacity(scratch)
    if args.analyze_evaluation:
        analyze_evaluation(args.analyze_evaluation)
