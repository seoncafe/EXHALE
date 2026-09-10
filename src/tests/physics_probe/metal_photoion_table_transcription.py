#!/usr/bin/env python3
"""The metal photoionization table in cross_sec.f90 against Verner's files.

The seventeen rows installed in src/modules/functions/cross_sec.f90 are a
transcription of two published tables, distributed by their author as
references/verner_photo/photo.dat (Verner, Ferland, Korista & Yakovlev 1996,
ApJ 465, 487, Table 1) and references/verner_photo/table1.dat (Verner &
Yakovlev 1995, A&AS 109, 125, Table 1), together with the subshell rows the
author's own routine references/verner_photo/phfit2.f carries in its PH1 DATA
statements.  src/utils/metal_photoion_table.py
emits that block and, with --check, reads it back and compares every field
with the file it came from, digit for digit, including E_max and the
handover the fits are joined at.  This driver runs that check.

The distribution lives outside the repository (it is third-party material in
the workspace's references/ tree), so its absence is reported as a skip with
the path that was looked for, not as a failure: what is then untested is the
transcription, and every other assertion about these cross sections, in
metal_photoionization_fits.f90, still runs.
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GEN = os.path.join(ROOT, "src", "utils", "metal_photoion_table.py")
REFDIR = os.path.abspath(
    os.path.join(ROOT, "..", "references", "verner_photo"))

needed = ["photo.dat", "table1.dat", "phfit2.f"]
missing = [f for f in needed if not os.path.isfile(os.path.join(REFDIR, f))]
if missing:
    print("SKIP metal_photoion_table_transcription measured=absent "
          "reference=%s tol=exact" % os.path.join(REFDIR, missing[0]))
    sys.exit(0)

r = subprocess.run([sys.executable, GEN, "--check", "--refdir", REFDIR],
                   capture_output=True, text=True)
sys.stdout.write(r.stdout)
sys.stderr.write(r.stderr)
sys.exit(r.returncode)
