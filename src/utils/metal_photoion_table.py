#!/usr/bin/env python3
"""Metal photoionization fit table for cross_sec.f90, from Verner's own files.

WHAT IT EMITS.  The parameter block of the seventeen photo-ionizable metal
stages EXHALE carries, in the column order of metal_photoion_sigma (which is
species_table's mion_iphot).  Each ion gets

  * the threshold that GATES the fit,
  * E_max, the energy at which the OUTER shell leaves the 1996
    Opacity-Project fit for the 1995 subshell fit of the same shell,
  * the seven parameters of the 1996 fit (E_0, sigma_0, y_a, P, y_w, y_0,
    y_1), valid on E_th <= E < E_max,
  * the five parameters of the 1995 fit of the same shell (E_0, sigma_0,
    y_a, P, y_w) and its orbital quantum number l, used at E >= E_max,
  * every SUBSHELL below the outer one that carries electrons in the
    ground state, with its own 1995 five parameters, its l, and the
    energy at which it turns on.  The photoabsorption of the ion is the
    sum of the outer shell and these, which is what phfit2's caller
    obtains by summing the routine over is = 1 ... 7.

SOURCES, all in references/verner_photo/ (the author's own distribution):

  photo.dat   Verner, Ferland, Korista & Yakovlev 1996, ApJ 465, 487,
              Table 1: Z, N, E_th, E_max, E_0, sigma_0, y_a, P, y_w, y_0,
              y_1, one row per (Z, N), for the OUTER shell of the Opacity
              Project elements.  Byte description in photo.txt.
  table1.dat  Verner & Yakovlev 1995, A&AS 109, 125, Table 1: Z, N, n, l,
              E_th, E_0, sigma_0, y_a, P, y_w, one row per subshell.
              Byte description in table1.txt.
  phfit2.f    the author's reference routine, and the source of the
              SUBSHELL parameters: its PH1 DATA statements are the 1995
              table with, in the author's own words, the "inner-shell
              ionization energies of some low-ionized species ...
              slightly improved to fit smoothly the experimental
              inner-shell ionization energies of neutral atoms".  Parsed
              here from that file, so the subshell rows are the ones the
              reference routine itself evaluates.  MEASURED by comparing
              the two files field by field over the seventeen ions: PH1
              and table1.dat differ in the THRESHOLD column of five
              subshell rows (Mg II 1s, Si II 1s, Ca II 1s and 2s,
              Fe II 1s and 2s) and nowhere else.  Its shell bookkeeping
              is transcribed in SHELL BOOKKEEPING below and is the reason
              E_max is the handover energy.

SHELL BOOKKEEPING, read from phfit2.f (subroutine phfit2, lines 26-56 of the
distributed file).  For an ion (nz, ne) it sets

    nout = ntot(ne)                         the outer shell,
    nout = 7  if nz == ne and nz > 18,      the 4s shell of K, Ca, ... Zn
    nout = 7  if nz == ne+1 and nz in {20, 21, 22, 25, 26}
    nint = ninn(ne)                         the outermost INNER shell,
    einn = 0                if nz in {15, 17, 19} or (nz > 20 and nz != 26)
    einn = 1e30             if ne < 3
    einn = ph1(1,nz,ne,nint)                otherwise,

and for the outer shell (is = nout > nint) evaluates the 1995 six-parameter
form when e >= einn and the 1996 nine-parameter form below it.  einn is the
threshold of the outermost inner shell, and it equals the E_max column of
photo.dat for every row (asserted below).  einn = 0 for potassium is why K I
takes the 1995 fit at every energy, and it is why K I has no photo.dat row.

The same lines fix when each SUBSHELL is > 1 contributes.  A true inner
shell (is <= nint) turns on at its own threshold ph1(1,nz,ne,is).  A shell
between nint and nout is part of the valence complex the 1996 outer-shell
fit already represents, so phfit2 returns zero for it below einn and its
1995 fit above einn; that is the line

    if(is.lt.nout.and.is.gt.nint.and.e.lt.einn)return

Since the thresholds of the inner shells all lie above einn, both cases are
one rule: the shell turns on at max(its threshold, E_max).  A shell with no
electrons in the ground state carries a placeholder row with sigma_0 = 0,
which is how phfit2 returns zero for it; those rows are dropped here.

l(is) = (0, 0, 1, 0, 1, 2, 0) for shells 1s, 2s, 2p, 3s, 3p, 3d, 4s enters
Q = 5.5 + l - P/2 of the 1995 form; the 1996 form has no l term.

THRESHOLDS.  The threshold only gates the fit; it does not enter the fitted
expression.  It must be the energy the photon grid uses as the ion's band
edge (species_table's mion_ethr) and charges the photoelectron h nu - E_th
against, so three ions take the NIST ionization potential in place of the
rounded value the fit tables print.  They are listed in ETH_OVERRIDE.

Usage
-----
    python3 metal_photoion_table.py [--refdir DIR]
        write the Fortran parameter block to stdout.

    python3 metal_photoion_table.py --check [--refdir DIR] [--source FILE]
        read the block installed in src/modules/functions/cross_sec.f90 and
        compare it, digit for digit, with photo.dat and table1.dat.  Prints
        one line per array and exits nonzero on any mismatch.
"""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_REFDIR = os.path.abspath(
    os.path.join(ROOT, "..", "references", "verner_photo"))
DEFAULT_SOURCE = os.path.join(ROOT, "src", "modules", "functions",
                              "cross_sec.f90")

# phfit2.f BLOCK DATA BDATA
L_SHELL = [0, 0, 1, 0, 1, 2, 0]
NINN = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
NTOT = [1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7]
# shell number -> (principal quantum number n, orbital quantum number l)
SHELL_NL = {1: (1, 0), 2: (2, 0), 3: (2, 1), 4: (3, 0), 5: (3, 1),
            6: (3, 2), 7: (4, 0)}

# The column order of metal_photoion_sigma, i.e. species_table's mion_iphot.
IONS = [("CI", 6, 6), ("CII", 6, 5), ("OI", 8, 8), ("OII", 8, 7),
        ("NI", 7, 7), ("NII", 7, 6), ("MgI", 12, 12), ("MgII", 12, 11),
        ("SiI", 14, 14), ("SiII", 14, 13), ("CaI", 20, 20),
        ("CaII", 20, 19), ("NaI", 11, 11), ("KI", 19, 19), ("SI", 16, 16),
        ("FeI", 26, 26), ("FeII", 26, 25)]

# Ions whose gate is the NIST ionization potential rather than the value the
# fit table prints.  Value: (Fortran literal, what it is, the printed value).
ETH_OVERRIDE = {
    "MgI":  ("e_th_MgI",  "NIST ASD, 7.646236 eV",  "7.646E+00"),
    "MgII": ("e_th_MgII", "NIST ASD, 15.035271 eV", "1.504E+01"),
    "FeII": ("1.6199d+01", "NIST ASD, 16.19920 eV", "1.619E+01"),
}


def read_photo_dat(refdir):
    """Verner et al. 1996 Table 1, keyed (Z, N).  Fields kept as strings."""
    rows = {}
    with open(os.path.join(refdir, "photo.dat")) as fh:
        for line in fh:
            if not line.strip():
                continue
            z = int(line[0:2])
            n = int(line[3:5])
            rows[(z, n)] = [line[6 + 10 * i:15 + 10 * i].strip()
                            for i in range(9)]
    return rows


def read_table1_dat(refdir):
    """Verner & Yakovlev 1995 Table 1, keyed (Z, N, n, l)."""
    rows = {}
    with open(os.path.join(refdir, "table1.dat")) as fh:
        for line in fh:
            if not line.strip():
                continue
            z = int(line[0:2])
            n = int(line[3:5])
            nq = int(line[6:7])
            ll = int(line[8:9])
            rows[(z, n, nq, ll)] = [line[10 + 11 * i:20 + 11 * i].strip()
                                    for i in range(6)]
    return rows


def read_ph1(refdir):
    """phfit2.f's PH1 DATA statements, keyed (Z, N, shell).  Strings as printed.

    Fixed-form continuation lines carry a nonblank in column 6; they are
    joined onto the statement they continue before the DATA line is read.
    """
    rows = {}
    stmt = []
    with open(os.path.join(refdir, "phfit2.f")) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line[:1] in ("*", "c", "C") or not line.strip():
                continue
            if len(line) > 5 and line[5] not in (" ", "\t") and stmt:
                stmt[-1] += line[6:]
            else:
                stmt.append(line)
    for line in stmt:
        m = re.match(r"\s*DATA\s*\(PH1\(I,\s*(\d+),\s*(\d+),\s*(\d+)\)"
                     r",I=1,6\)\s*/(.*)/", line.strip())
        if m:
            key = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
            vals = [v.strip() for v in m.group(4).split(",")]
            if len(vals) != 6:
                raise RuntimeError("PH1 row %s has %d fields" % (key, len(vals)))
            rows[key] = vals
    if not rows:
        raise RuntimeError("no PH1 DATA statements found in phfit2.f")
    return rows


def outer_shell(z, ne):
    nout = NTOT[ne - 1]
    if z == ne and z > 18:
        nout = 7
    if z == ne + 1 and z in (20, 21, 22, 25, 26):
        nout = 7
    return nout


def inner_threshold_is_zero(z):
    return z in (15, 17, 19) or (z > 20 and z != 26)


def build(refdir):
    """One record per ion, every number a string exactly as its file prints it."""
    photo = read_photo_dat(refdir)
    t1 = read_table1_dat(refdir)
    ph1 = read_ph1(refdir)
    out = []
    for name, z, ne in IONS:
        nout = outer_shell(z, ne)
        nq, ll = SHELL_NL[nout]
        vy95 = t1[(z, ne, nq, ll)]
        if inner_threshold_is_zero(z):
            # phfit2 sets einn = 0: the 1995 fit at every energy, and the
            # 1996 table has no row for this ion.
            if (z, ne) in photo:
                raise RuntimeError(
                    "%s has einn = 0 but a 1996 row exists" % name)
            e_th = vy95[0]
            e_max = "0.000E+00"
            vfky96 = ["0.000E+00"] * 7
        else:
            row = photo[(z, ne)]
            e_th, e_max, vfky96 = row[0], row[1], row[2:9]
            # The handover energy of phfit2 is the inner-shell threshold; it
            # must be the E_max the 1996 table prints, or the two fits are
            # joined at different energies.
            nint = NINN[ne - 1]
            nq_i, l_i = SHELL_NL[nint]
            einn = float(t1[(z, ne, nq_i, l_i)][0])
            if abs(einn - float(e_max)) > 1.0e-12 * max(einn, 1.0):
                raise RuntimeError(
                    "%s: E_max %s is not the inner-shell threshold %g"
                    % (name, e_max, einn))
        # The subshells below the outer one.  Turn-on = max(threshold,
        # E_max), the one rule phfit2's two return lines amount to; a
        # placeholder row of an empty shell (sigma_0 = 0) is dropped.
        sub = []
        for is_ in range(1, nout):
            row = ph1.get((z, ne, is_))
            if row is None:
                raise RuntimeError("%s: no PH1 row for shell %d" % (name, is_))
            if float(row[2]) == 0.0:
                continue
            # Keep whichever of the two the file prints, so the emitted
            # number stays digit for digit what its source carries.
            e_on = row[0] if float(row[0]) >= float(e_max) else e_max
            sub.append(dict(shell=is_, l=SHELL_NL[is_][1],
                            e_on=e_on, vy95=row[1:6]))
        out.append(dict(name=name, z=z, ne=ne, nout=nout, l=ll,
                        e_th=e_th, e_max=e_max, vfky96=vfky96,
                        vy95=vy95[1:6], sub=sub))
    return out


def fortran(x):
    """A file field as a Fortran double literal, digit for digit."""
    return x.replace("E", "d").replace("e", "d")


def undo_fortran(x):
    return x.replace("d", "E").replace("D", "E").upper()


def emit_array(name, values, comment, width=6):
    lines = ["      ! " + comment,
             "      real*8, parameter :: %s(n_metal_photo) = [ &" % name]
    for i in range(0, len(values), width):
        chunk = ", ".join("%-11s" % v for v in values[i:i + width])
        tail = ", &" if i + width < len(values) else " ]"
        lines.append("           " + chunk.rstrip() + tail)
    return lines


def emit(records):
    out = []
    names = ", ".join(r["name"] for r in records)
    out.append("      ! Column order (= species_table's mion_iphot): " )
    for i in range(0, len(records), 6):
        out.append("      !   " +
                   ", ".join(r["name"] for r in records[i:i + 6]))
    out.append("      integer, parameter :: n_metal_photo = %d"
               % len(records))
    eth = []
    for r in records:
        if r["name"] in ETH_OVERRIDE:
            eth.append(ETH_OVERRIDE[r["name"]][0])
        else:
            eth.append(fortran(r["e_th"]))
    out += emit_array("mph_e_th", eth,
                      "Threshold that gates the fit [eV].")
    out += emit_array("mph_e_max", [fortran(r["e_max"]) for r in records],
                      "E_max [eV]: above it the outer shell takes the 1995 fit.")
    lbl96 = ["E_0 [eV]", "sigma_0 [Mb]", "y_a", "P", "y_w", "y_0", "y_1"]
    nm96 = ["mph96_E0", "mph96_s0", "mph96_ya", "mph96_P", "mph96_yw",
            "mph96_y0", "mph96_y1"]
    for j in range(7):
        out += emit_array(nm96[j],
                          [fortran(r["vfky96"][j]) for r in records],
                          "Verner et al. 1996 Table 1, %s." % lbl96[j])
    lbl95 = ["E_0 [eV]", "sigma_0 [Mb]", "y_a", "P", "y_w"]
    nm95 = ["mph95_E0", "mph95_s0", "mph95_ya", "mph95_P", "mph95_yw"]
    for j in range(5):
        out += emit_array(nm95[j],
                          [fortran(r["vy95"][j]) for r in records],
                          "Verner & Yakovlev 1995 Table 1, %s." % lbl95[j])
    out.append("      ! Orbital quantum number of the outer shell, in")
    out.append("      ! Q = 5.5 + l - P/2 of the 1995 form.")
    out.append("      integer, parameter :: mph_l(n_metal_photo) = [ &")
    vals = [str(r["l"]) for r in records]
    out.append("           " + ", ".join(vals) + " ]")
    out += emit_subshells(records)
    return out


def emit_subshells(records):
    """The subshell block: one column of n_sub_max slots per ion."""
    nmax = max(len(r["sub"]) for r in records)
    out = ["      ! ---- The subshells below the outer one ----",
           "      !",
           "      ! Slots per ion, padded with sigma_0 = 0 rows that the sum",
           "      ! never reaches (the loop runs to mph_n_sub).  Shell order",
           "      ! is phfit2's: 1s, 2s, 2p, 3s, 3p, 3d.  Every row is the",
           "      ! 1995 five-parameter form, evaluated above mph_sub_e_on.",
           "      integer, parameter :: n_metal_subshell = %d" % nmax]
    out += emit_array_i("mph_n_sub", [str(len(r["sub"])) for r in records],
                        "Subshells this ion carries electrons in.")
    out.append("      ! Turn-on [eV] = max(shell threshold, E_max).")
    out += emit_matrix("mph_sub_e_on",
                       [[fortran(x["e_on"]) for x in r["sub"]] for r in records],
                       nmax, "0.000d+00", records)
    lbl = ["E_0 [eV]", "sigma_0 [Mb]", "y_a", "P", "y_w"]
    nms = ["mph_sub_E0", "mph_sub_s0", "mph_sub_ya", "mph_sub_P",
           "mph_sub_yw"]
    for j in range(5):
        out.append("      ! Verner & Yakovlev 1995 subshell fit, %s." % lbl[j])
        out += emit_matrix(nms[j],
                           [[fortran(x["vy95"][j]) for x in r["sub"]]
                            for r in records],
                           nmax, "0.000d+00", records)
    out.append("      ! Orbital quantum number l of each subshell.")
    out += emit_matrix_i("mph_sub_l",
                         [[str(x["l"]) for x in r["sub"]] for r in records],
                         nmax, "0", records)
    return out


def emit_matrix(name, cols, nmax, pad, records):
    return _emit_matrix("real*8", name, cols, nmax, pad, records)


def emit_matrix_i(name, cols, nmax, pad, records):
    return _emit_matrix("integer", name, cols, nmax, pad, records)


SHELL_NAME = {1: "1s", 2: "2s", 3: "2p", 4: "3s", 5: "3p", 6: "3d", 7: "4s"}


def _emit_matrix(kind, name, cols, nmax, pad, records):
    lines = ["      %s, parameter :: %s(n_metal_subshell,n_metal_photo) = &"
             % (kind, name),
             "           reshape( [ &"]
    for i, col in enumerate(cols):
        vals = col + [pad] * (nmax - len(col))
        tail = ", &" if i + 1 < len(cols) else " ], &"
        shells = " ".join(SHELL_NAME[x["shell"]] for x in records[i]["sub"])
        lines.append("           " + ", ".join("%-11s" % v for v in vals).rstrip()
                     + tail + "   ! %-5s %s" % (records[i]["name"], shells))
    lines.append("           [n_metal_subshell,n_metal_photo] )")
    return lines


def emit_array_i(name, values, comment, width=17):
    lines = ["      ! " + comment,
             "      integer, parameter :: %s(n_metal_photo) = [ &" % name]
    for i in range(0, len(values), width):
        chunk = ", ".join(values[i:i + width])
        tail = ", &" if i + width < len(values) else " ]"
        lines.append("           " + chunk + tail)
    return lines


ARRAY_RE = r"real\*8,\s*parameter\s*::\s*%s\(n_metal_photo\)\s*=\s*\[(.*?)\]"


def parse_installed(source):
    text = open(source).read()
    # strip Fortran continuation markers and comments
    text = re.sub(r"!.*", "", text)
    text = re.sub(r"&\s*\n\s*", " ", text)
    found = {}
    for m in re.finditer(
            r"(?:real\*8|integer),\s*parameter\s*::\s*(mph_sub_\w+)"
            r"\(n_metal_subshell,n_metal_photo\)\s*=\s*reshape\(\s*\[([^\]]*)\]",
            text):
        found[m.group(1)] = [t.strip() for t in m.group(2).split(",")
                             if t.strip()]
    m = re.search(
        r"integer,\s*parameter\s*::\s*mph_n_sub\(n_metal_photo\)\s*=\s*\[([^\]]*)\]",
        text)
    if m:
        found["mph_n_sub"] = [t.strip() for t in m.group(1).split(",")
                              if t.strip()]
    for m in re.finditer(
            r"real\*8,\s*parameter\s*::\s*(mph\w+)\(n_metal_photo\)\s*=\s*\[([^\]]*)\]",
            text):
        found[m.group(1)] = [t.strip() for t in m.group(2).split(",")
                             if t.strip()]
    m = re.search(
        r"integer,\s*parameter\s*::\s*mph_l\(n_metal_photo\)\s*=\s*\[([^\]]*)\]",
        text)
    if m:
        found["mph_l"] = [t.strip() for t in m.group(1).split(",")
                          if t.strip()]
    return found


def check(refdir, source):
    records = build(refdir)
    installed = parse_installed(source)
    expect = {}
    expect["mph_e_th"] = [ETH_OVERRIDE[r["name"]][0] if r["name"] in
                          ETH_OVERRIDE else fortran(r["e_th"])
                          for r in records]
    expect["mph_e_max"] = [fortran(r["e_max"]) for r in records]
    for j, nm in enumerate(["mph96_E0", "mph96_s0", "mph96_ya", "mph96_P",
                            "mph96_yw", "mph96_y0", "mph96_y1"]):
        expect[nm] = [fortran(r["vfky96"][j]) for r in records]
    for j, nm in enumerate(["mph95_E0", "mph95_s0", "mph95_ya", "mph95_P",
                            "mph95_yw"]):
        expect[nm] = [fortran(r["vy95"][j]) for r in records]
    expect["mph_l"] = [str(r["l"]) for r in records]
    nmax = max(len(r["sub"]) for r in records)
    expect["mph_n_sub"] = [str(len(r["sub"])) for r in records]
    for nm, key, j, pad in [("mph_sub_e_on", "e_on", None, "0.000d+00"),
                            ("mph_sub_E0", "vy95", 0, "0.000d+00"),
                            ("mph_sub_s0", "vy95", 1, "0.000d+00"),
                            ("mph_sub_ya", "vy95", 2, "0.000d+00"),
                            ("mph_sub_P", "vy95", 3, "0.000d+00"),
                            ("mph_sub_yw", "vy95", 4, "0.000d+00"),
                            ("mph_sub_l", "l", None, "0")]:
        flat = []
        for r in records:
            for x in r["sub"]:
                v = x[key] if j is None else x[key][j]
                flat.append(str(v) if nm == "mph_sub_l" else fortran(v))
            flat += [pad] * (nmax - len(r["sub"]))
        expect[nm] = flat

    failures = 0
    for nm in sorted(expect):
        want = expect[nm]
        got = installed.get(nm)
        if got is None:
            print("FAIL metal_photoion_table_%s missing from %s"
                  % (nm, os.path.relpath(source, ROOT)))
            failures += 1
            continue
        if len(got) != len(want):
            print("FAIL metal_photoion_table_%s measured=%d entries "
                  "reference=%d entries tol=exact"
                  % (nm, len(got), len(want)))
            failures += 1
            continue
        if len(want) == len(records):
            tag = [r["name"] for r in records]
        else:
            tag = ["%s_%d" % (r["name"], j + 1) for r in records
                   for j in range(nmax)]
        bad = [(tag[i], got[i], want[i])
               for i in range(len(want))
               if undo_fortran(got[i]) != undo_fortran(want[i])]
        if bad:
            for ion, g, w in bad:
                print("FAIL metal_photoion_table_%s_%s measured=%s "
                      "reference=%s tol=exact" % (nm, ion, g, w))
            failures += len(bad)
        else:
            print("PASS metal_photoion_table_%s measured=%d_fields_match "
                  "reference=%d_fields tol=exact" % (nm, len(want),
                                                     len(want)))
    return failures


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refdir", default=DEFAULT_REFDIR)
    ap.add_argument("--source", default=DEFAULT_SOURCE)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if args.check:
        sys.exit(1 if check(args.refdir, args.source) else 0)
    print("\n".join(emit(build(args.refdir))))


if __name__ == "__main__":
    main()
