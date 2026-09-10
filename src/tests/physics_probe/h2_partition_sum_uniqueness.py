#!/usr/bin/env python3
"""One Boltzmann sum over the H2 rovibrational ladder in the production source.

The H2 level ladder of Roueff et al. (2019) fixes three things at once: the
internal energy the caloric equation of state carries, the free energy the
H + H <-> H2 equilibrium constant is built from, and the LTE populations the
H2 quadrupole lines of the infrared coolant are weighted by.  They are the
zeroth, first and second moments of ONE Boltzmann sum, so the source must
form that sum once.  A second copy is not wrong the day it is written -- the
one this test was written for summed the same 302 levels and agreed with the
first to the bit, measured over 100 to 8000 K -- it is wrong the day one copy
gains a level, a weight or a guard the other does not.

WHAT IS COMPARED.  Every .f90 file under src/modules is scanned, with
comments and character strings removed, for a reference to the ladder arrays
h2_lev_T or h2_lev_g.

REFERENCE.  Exactly two files may name them: the module that declares and
fills them, and the module that forms the sum, whose h2_partition_function
every other consumer calls.  Tolerance zero.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MODULES = os.path.join(ROOT, "src", "modules")

LADDER = re.compile(r"(?<![A-Za-z0-9_])h2_lev_(?:T|g)(?![A-Za-z0-9_])",
                    re.IGNORECASE)

# The file that declares and fills the ladder, and the file that forms the
# one Boltzmann sum over it.
ALLOWED = (
    os.path.join("src", "modules", "lower_atmosphere",
                 "molecular_infrared_data.f90"),
    os.path.join("src", "modules", "states", "caloric_eos.f90"),
)


def code_only(text):
    """Blank full-line and trailing comments and the contents of strings.

    Lines are blanked rather than dropped so that the reported line number
    is the line number in the file.
    """
    out = []
    for line in text.splitlines():
        if line[:1] in ("!", "C", "c", "*") and not line[:1].isspace():
            out.append("")
            continue
        line = re.sub(r"'[^']*'", "''", line)
        line = re.sub(r'"[^"]*"', '""', line)
        line = line.split("!", 1)[0]
        out.append(line)
    return out


def main():
    failures = []
    for base, dirs, files in os.walk(MODULES):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for name in sorted(files):
            if not name.endswith(".f90"):
                continue
            path = os.path.join(base, name)
            rel = os.path.relpath(path, ROOT)
            if rel in ALLOWED:
                continue
            with open(path, "r", errors="replace") as fh:
                lines = code_only(fh.read())
            for n, line in enumerate(lines, 1):
                if LADDER.search(line):
                    failures.append("%s:%d: %s" % (rel, n, line.strip()))

    name = "h2_partition_sum_uniqueness"
    if failures:
        print("FAIL %s measured=%d reference=0 tol=0 second sums over the "
              "H2 ladder" % (name, len(failures)))
        for f in failures:
            print("  " + f)
        return 1
    print("PASS %s measured=0 reference=0 tol=0 references to the H2 ladder "
          "outside its declaration and its partition function" % name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
