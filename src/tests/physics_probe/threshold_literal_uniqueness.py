#!/usr/bin/env python3
"""One definition of each ionization threshold in the production source.

The photon grid puts every ionization threshold on a bin edge, the
photoelectron energy is h nu - e_th, the collisional-ionization cooling
removes e_th per event, and the photoionization cross sections turn on
there.  Those are the same physical energy, so the source must state it
once.  When it was stated twice the two copies drifted: e_th_HeI = 24.6
against the He I fit's own E_th = 24.59 left the band [24.59, 24.60] eV,
0.115 per cent of the He I photoionization rate of the default power law,
integrated as zero (measured 2026-09-05).

WHAT IS COMPARED.  Every .f90 file under src/, except src/tests, is
scanned for a numeric literal whose VALUE is one of the ionization
potentials of H I, He I, He II, He 2^3S or H2, in any of the roundings the
code has carried.  Comments and the contents of character strings are
removed first: a printed label is text, not a constant.

REFERENCE.  Exactly one such literal per potential, and all of them in the
declaration of the named constant in src/modules/init/parameters.f90.
Tolerance zero: any other occurrence is a failure, unless it appears in
threshold_literal_exceptions.txt with a reason.

An exception marked `open` is a duplicate this test cannot accept as
correct but whose file was not this item's to edit; it is printed as a
DIAGNOSTIC so that it stays visible.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC = os.path.join(ROOT, "src")
EXCEPTIONS = os.path.join(HERE, "threshold_literal_exceptions.txt")

# The ionization potentials, in every rounding this source has carried.
# Matching is on the parsed VALUE, so 4.8d-10 (a different quantity that
# happens to share digits) is not a match and 13.62 (the O I potential) is
# not one either.
POTENTIALS = {
    "H I": [13.6, 13.60, 13.598, 13.5984, 13.5984346, 13.598434599],
    "He I": [24.6, 24.60, 24.59, 24.5874, 24.587389, 24.587389011],
    "He II": [54.4, 54.40, 54.4178, 54.417765, 54.417765486],
    "He 2^3S": [4.8, 4.80, 4.78, 4.7678, 4.767775],
    "H2": [15.4, 15.40, 15.4259, 15.425927, 15.425933],
    "H(n=2)": [3.4, 3.40, 3.3996, 3.399608650],
}
VALUE_OF = {}
for species, values in POTENTIALS.items():
    for v in values:
        VALUE_OF[v] = species

# The one file allowed to write these values down, and the constants whose
# declaration is that definition.
DEFINITION_FILE = os.path.join("src", "modules", "init", "parameters.f90")
DEFINED = ("e_th_HI", "e_th_HeI", "e_th_HeII", "e_th_HeTR", "e_th_H2")

NUMBER = re.compile(r"(?<![A-Za-z0-9_.])([0-9]+\.[0-9]*)([dDeE][+-]?[0-9]+)?")


def code_only(line):
    """The executable text of one Fortran line: no comment, no string body."""
    out = []
    quote = None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
            continue                      # string contents are text, not code
        if ch in "'\"":
            quote = ch
            continue
        if ch == "!":
            break
        out.append(ch)
    return "".join(out)


def parse(token):
    return float(token.replace("d", "e").replace("D", "E"))


def load_exceptions():
    """file basename | code substring | kind (justified|open) | reason"""
    entries = []
    if not os.path.exists(EXCEPTIONS):
        return entries
    for raw in open(EXCEPTIONS):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("|", 3)]
        if len(parts) != 4:
            print("FAIL threshold_literal_exceptions_syntax "
                  "measured=%r reference=4_fields tol=0" % raw.rstrip())
            return None
        entries.append(parts)
    return entries


def main():
    exceptions = load_exceptions()
    if exceptions is None:
        return 1
    used = [False] * len(exceptions)

    hits = []
    for base, _dirs, files in os.walk(SRC):
        parts = base.split(os.sep)
        # editor auto-save copies (.ipynb_checkpoints) are not source
        if "tests" in parts or ".ipynb_checkpoints" in parts:
            continue
        for name in sorted(files):
            if not name.endswith(".f90"):
                continue
            path = os.path.join(base, name)
            rel = os.path.relpath(path, ROOT)
            with open(path, errors="replace") as fh:
                for lineno, raw in enumerate(fh, 1):
                    code = code_only(raw)
                    for match in NUMBER.finditer(code):
                        token = match.group(1) + (match.group(2) or "")
                        value = parse(token)
                        for known, species in VALUE_OF.items():
                            if abs(value - known) <= 1e-12 * max(1.0, abs(known)):
                                hits.append((rel, lineno, token, species,
                                             code.strip()))
                                break
                        else:
                            continue
                        break

    unlisted = []
    diagnostics = []
    for rel, lineno, token, species, code in hits:
        if rel == DEFINITION_FILE and any(d in code for d in DEFINED):
            continue
        matched = None
        for index, (efile, esub, ekind, ereason) in enumerate(exceptions):
            if os.path.basename(rel) == efile and esub in code:
                matched = (index, ekind, ereason)
                break
        if matched is None:
            unlisted.append((rel, lineno, token, species, code))
            continue
        index, ekind, ereason = matched
        used[index] = True
        if ekind == "open":
            diagnostics.append("DIAGNOSTIC open duplicate of the %s potential: "
                               "%s:%d [%s] -- %s" % (species, rel, lineno,
                                                     token, ereason))

    for line in diagnostics:
        print(line)
    for index, entry in enumerate(exceptions):
        if not used[index]:
            print("DIAGNOSTIC exception no longer matches any line: "
                  "%s | %s" % (entry[0], entry[1]))

    for rel, lineno, token, species, code in unlisted:
        print("  %s:%d [%s] is the %s potential: %s"
              % (rel, lineno, token, species, code[:100]))

    ok = not unlisted
    print("%s threshold_literals_written_once measured=%d reference=0 tol=0"
          % ("PASS" if ok else "FAIL", len(unlisted)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
