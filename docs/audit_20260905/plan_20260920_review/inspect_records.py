"""Read the exact generations used to review PLAN_20260920; never modify a case."""
from pathlib import Path
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parents[3]
MODELS = ROOT / "LHS1140b/models"
CASES = [
    "atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7",
    "atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13",
    "atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7",
    "molecular_scalar_gj1132_kzz1e9/HeH0.083",
    "molecular_scalar_gj1132_wellmixed/HeH0.083",
    "molecular_scalar_gj1132_kzz1e9/HeH0.55",
    "molecular_scalar_gj1132_wellmixed/HeH0.55",
    "molecular_scalar_gj1132_kzz1e9/HeH2.13",
    "molecular_photochem_gj1132_kzzprofile/HeH9",
]
for case in CASES:
    directory = MODELS / case
    print("\nCASE", case)
    index_path = directory / "state_index.json"
    if not index_path.exists():
        print("NO INDEX")
        continue
    index = json.loads(index_path.read_text())
    print("INDEX", index.get("updated_at"), "latest_complete=",
          index.get("latest_complete"), "latest_certified=",
          index.get("latest_certified"))
    selected = next(g for g in index["generations"]
                    if g["generation_id"] == index["latest_complete"])
    state = directory / selected["path"]
    manifest = json.loads((state / "manifest.json").read_text())
    print("BINARY", manifest.get("source_identity"))
    print("BOUNDARY", manifest.get("model_identity", {}).get("boundary_model"))
    print("ENDING", manifest.get("ending"))
    for name in ("Hydro_ioniz.txt", "Ion_species.txt"):
        product = state / name
        print("MD5", str(product.relative_to(ROOT)),
              hashlib.md5(product.read_bytes()).hexdigest())
    certificate = state / "certification.txt"
    if certificate.exists():
        for line in certificate.read_text().splitlines():
            if any(s in line for s in ("hydrodynamic", "carrier balance H2",
                                      "gated at", "out-of-domain", "H3+")):
                print("CERT", line.strip())
    log = directory / "run.log"
    if log.exists():
        lines = [s.strip() for s in log.read_text(errors="replace").splitlines()
                 if "du =" in s or "norm(R)" in s]
        for line in lines[-2:]:
            print("CASE LOG (not generation identity)", line)

print("\nARITHMETIC: M2 mass/tolerance", 0.2254 / 2.2e-11,
      "decades", math.log10(0.2254 / 2.2e-11))
print("TOY COUNTEREXAMPLE, not an EXHALE simulation:")
x, approximate_jacobian = 1.0, -1.0
step = -x / approximate_jacobian
print("F(x)=x, approximate J=-1: linear residual=",
      approximate_jacobian * step + x, "nonlinear residual=", x + step)
value, scale = 0.1, 3.0
print("FLOATING-POINT CONTRACT: (0.1*3)/3 == 0.1:",
      (value * scale) / scale == value,
      "difference=", (value * scale) / scale - value)
steady = ROOT / "src/modules/time_step/steady_newton.f90"
print("SOURCE SHA256", str(steady.relative_to(ROOT)),
      hashlib.sha256(steady.read_bytes()).hexdigest())
print("H3+ reset references (a definition is not a call):")
for source in sorted((ROOT / "src").rglob("*.f90")):
    for number, line in enumerate(source.read_text(errors="replace").splitlines(), 1):
        if re.search(r"\bh3p_reset_domain_records\b", line, re.I):
            print(str(source.relative_to(ROOT)), number, line.strip())
