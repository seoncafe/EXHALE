"""Compare the legacy log selection with selection inside the last solve.

Usage: python3 handoff_20260917_log_scope_probe.py PATH_TO_RUN_LOG
This reads an existing log and does not execute or modify EXHALE.
"""
import re
import sys
from pathlib import Path

lines = Path(sys.argv[1]).read_text().splitlines()
done = [i for i, line in enumerate(lines)
        if re.match(r" \((JFNK|PTC)\) done info=", line)]
assert done, "No completed steady solve"
norm = re.compile(r"\|\|R\|\|=\s*(\S+)")
best = re.compile(r" \((JFNK|PTC)\) returning best iterate")
iteration = re.compile(r" \((JFNK|PTC)\) it +[0-9]+ +\|\|R\|\|=")

def select(begin, end):
    candidates = [i for i in range(begin, end) if best.match(lines[i])]
    if not candidates:
        candidates = [i for i in range(begin, end) if iteration.match(lines[i])]
    assert candidates, "No iteration record in this solve"
    i = candidates[-1]
    return i + 1, float(norm.search(lines[i]).group(1))

gate = float(norm.search(lines[done[-1]]).group(1))
legacy_line, legacy = select(0, done[-1])
scoped_line, scoped = select(done[-2] + 1 if len(done) > 1 else 0, done[-1])
print(f"Final done: line={done[-1] + 1}, norm={gate:.9e}")
print(f"Legacy selection: line={legacy_line}, norm={legacy:.9e}, ratio={gate / legacy:.9e}")
print(f"Scoped selection: line={scoped_line}, norm={scoped:.9e}, ratio={gate / scoped:.9e}")
print("This diagnoses record selection, not the accuracy of unprinted internal norms.")
