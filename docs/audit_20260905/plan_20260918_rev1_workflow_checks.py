"""Read-only workflow checks for the September 18 revised-plan review.

Calls the actual seed-selection functions with mocked inputs and extracts
the mapper's actual nested file-selection function. No catalog is edited.
"""
import ast
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import mock_open, patch

sys.dont_write_bytecode = True
project = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "review_seed_selector", project / "LHS1140b/models/pick_seed.py")
selector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selector)

base_keys = {
    "Spectrum file": "same_spectrum.txt",
    "Molecular chemistry": "False",
    "He_diffusion": "True",
    "He_Kzz": "1e9",
    "He/H number ratio": "2.13",
}
with patch.object(selector.os.path, "isfile", return_value=False), \
        patch.object(selector, "base_inp_reservoirs", return_value=False):
    with patch.object(selector, "read_keys", return_value={
            **base_keys, "Ionization transport": "False"}):
        local = selector.Physics("/synthetic/local")
    with patch.object(selector, "read_keys", return_value={
            **base_keys, "Ionization transport": "True"}):
        transported = selector.Physics("/synthetic/transported")
print("D9 actual Physics.matches ignores ionization-key difference:",
      local.matches(transported))

historical_record = "- certification of the state written: **CERTIFIED**\n"
with patch.object(selector.os.path, "isfile", return_value=True), \
        patch("builtins.open", mock_open(read_data=historical_record)):
    accepted = selector.certified("/synthetic/stale_record")
print("D8 actual certified() accepts record without state-identity check:", accepted)

mapper_tree = ast.parse((project / "src/utils/map_state_to_grid.py").read_text())
file_selector = next(node for node in ast.walk(mapper_tree)
                     if isinstance(node, ast.FunctionDef) and node.name == "src_file")
namespace = {"os": os, "src": "/synthetic/mixed_generation"}
exec(compile(ast.Module(body=[file_selector], type_ignores=[]),
             "actual_mapper_src_file", "exec"), namespace)
with patch.object(os.path, "exists", side_effect=lambda path:
                  path.endswith("/Hydro_ioniz.txt")):
    hydro = namespace["src_file"]("Hydro_ioniz")
    species = namespace["src_file"]("Ion_species")
print("D8 actual mapper selection with one missing product:", hydro, species)

# This is a reduced check of the campaign's shell pipeline, not a campaign run.
pipeline = subprocess.run(
    ["bash", "-c", 'r=$(false 2>&1 | tail -n 1); echo "case: $r"'],
    capture_output=True, text=True, check=False)
print("D8 campaign-style pipeline status after failed task:", pipeline.returncode)

print("D9 binding-cell average slope, pass 1 to 25:", (491 - 217) / 24)
print("D9 binding-cell average slope, pass 2 to 25:", (491 - 416) / 23)
print("D9 500 cells / claimed 3.3 cells each pass:", 500 / 3.3)

# A new seed can change the helium budget while preserving HeH+ and H nuclei.
h_total, hehp, target_ratio = 1.0, 0.1, 0.5
atomic_helium = target_ratio * h_total - hehp
print("D5 feasible new-seed atomic helium with nonzero HeH+:", atomic_helium)
print("D5 resulting total He/H:", (atomic_helium + hehp) / h_total)

assert local.matches(transported)
assert accepted
assert hydro.endswith("Hydro_ioniz.txt") and species.endswith("Ion_species_IC.txt")
assert pipeline.returncode == 0
assert atomic_helium >= 0.0
