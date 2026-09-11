"""Print measurements from the archived runs without altering their products."""

import csv
from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parent / "runs"
paths = [Path(arg) for arg in sys.argv[1:]] or sorted(root.iterdir())
for directory in paths:
    print(f"\nRUN {directory.name}")
    log = (directory / "run.log").read_text()
    state_lines = re.findall(r"^STATE .*|^COMPOSITION .*|^EXPERIMENT total_seconds=.*", log, re.M)
    print("\n".join(state_lines))
    table = directory / "rows.tsv"
    if not table.exists():
        continue
    stages = {}
    with table.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            stages.setdefault(row["stage"], []).append(row)
    print("stage mass momentum energy worst_species_wind H2_whole")
    for stage, rows in stages.items():
        hydro = [float(row["maximum"]) for row in rows if row["row"].startswith("hydrodynamic")]
        species = [float(row["gated"]) for row in rows
                   if row["row"].startswith(("elemental transport", "carrier balance"))]
        h2 = [float(row["maximum"]) for row in rows if row["row"] == "carrier balance H2"]
        print(stage, *[f"{value:.9e}" for value in hydro],
              f"{max(species, default=0):.9e}", f"{max(h2, default=0):.9e}")
    steps = directory / "steps.tsv"
    if steps.exists():
        print(steps.read_text())
    status = directory / "exit_status.log"
    if status.exists():
        print(status.read_text().strip())
