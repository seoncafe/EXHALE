#!/usr/bin/env python3
"""Relaxation of a metal-ion subshell vacancy, for cross_sec.f90.

WHAT IT EMITS.  For every subshell slot of the seventeen photo-ionizable
metal stages (the slots of the mph_sub_* arrays of cross_sec.f90, in the
same order), what a photoabsorption in that subshell does beyond ejecting
the photoelectron h nu - E_th,s:

  mph_sub_e_th        E_th,s, the threshold of the subshell itself [eV]
                      (phfit2.f PH1, first field).  mph_sub_e_on is
                      max(E_th,s, E_max), where the fit turns on; the
                      photoelectron is charged E_th,s.
  mph_sub_p_multi     the probability that the vacancy decays by at least
                      one Auger (autoionizing) transition, i.e. that the
                      absorption ejects two or more electrons.
  mph_sub_e_auger     the mean total kinetic energy of the Auger electrons
                      of one absorption [eV].
  mph_sub_e_fluor     the mean energy that leaves as photons (fluorescence
                      lines and the radiative decay of the excited ion the
                      cascade ends in) [eV].
  mph_sub_e_ion_multi the mean ionization energy spent on the electrons
                      beyond the first [eV]: sum over n of P_n times
                      I_(i+1) + ... + I_(i+n-1).

The four close the energy of the vacancy exactly,

    E_th,s = I_i + e_ion_multi + e_auger + e_fluor ,              (*)

with I_i the ionization potential the code charges the ion (mion_ethr, the
mph_e_th of cross_sec), so the photon h nu is shared among the photoelectron
(h nu - E_th,s), the Auger electrons, the escaping photons and the
ionization energy of every electron removed.

SOURCES.

  Kaastra, J. S., & Mewe, R. 1993, A&AS 97, 443, "X-ray emission from thin
  plasmas. I. Multiple Auger ionisation and fluorescence processes for Be
  to Zn".  Table 2: for each ion (st = 1 neutral) and each vacancy
  s = 1..7 (K, L1, L2, L3, M1, M2, M3), the probability, times 1e4, that
  the vacancy ends with n electrons ejected, the photoelectron included
  (columns 7-16, n = 1..10); "only the vacancies are listed for which Auger
  or fluorescent processes occur" (p. 447).  Table 3: each fluorescent line
  of each vacancy with its energy and its yield omega (photons per
  primary vacancy; "all 11732 lines stronger than 5 x 10^-5
  photons/ionisation", p. 449).  Both are read from the CDS copy of the
  paper's tables, J/A+AS/97/443, in references/kaastra_mewe_1993/
  (table2, table3, ReadMe); the probability columns of the Si I and Si II
  K, L1, L2 and L3 rows were compared against the printed Table 2 on p. 451.  Their Eq. (5), I = E_A + E_p + E_b,
  is the vacancy energy shared among Auger electrons, photons and binding;
  their Eq. (6)-(7) state that their tabulated E_A misses it by epsilon I
  and recommend E'_A = E_A - epsilon I, "forcing Eq. (5) to hold".  (*) is
  the same closure on this code's own thresholds: e_auger is the remainder
  E_th,s - I_i - e_ion_multi - e_fluor, and the script prints it beside
  their E'_A.

  Verner, D. A., & Yakovlev, D. G. 1995, A&AS 109, 125, through phfit2.f
  (references/verner_photo/): E_th,s, the first field of each PH1 row.

  NIST Atomic Spectra Database, ionization energies of C, N, O, Na, Mg, Si,
  S, K, Ca, Fe, all stages, retrieved 2026-09-25
  (references/nist_asd/ionization_energies_C_N_O_Na_Mg_Si_S_K_Ca_Fe.csv):
  the potentials of the stages this code does not photoionize.  For the
  stages it does (I_i, and I_(i+1) of a neutral whose next stage is
  photo-ionizable) the code's own mion_ethr is used, so that (*) holds in
  the energies the balance charges.

SHELL MAP.  Verner's 1s, 2s, 3s are Kaastra & Mewe's K, L1, M1.  Verner's
2p and 3p are single subshells, theirs are split by j: 2p = L2 (2p1/2, two
electrons) + L3 (2p3/2, four), 3p = M2 + M3 alike.  Every 2p and 3p that is
a subshell here is full, and the two j components are combined in the
ratio of their occupancies, 1/3 and 2/3 (equal cross section for each
electron; the spin-orbit splitting, 0.3 eV for Mg 2p and 0.4 eV for Fe 3p
in their Table 2, is below the resolution of the photon grid).  Verner's 3d
(Fe) has no row in either table: its vacancy is a valence hole.

A VACANCY WITH NO TABLE 2 ROW does not autoionize: one electron, the
photoelectron, and the ion is left excited by E_th,s - I_i, which it
radiates (Table 3 lists that line for Na I 2p, Mg II 2p, K I 3p, Ca II 3p
and the valence 2s and 3s holes, omega = 1 within their 5e-5 cut, at
their approximate energies; the energy used here is E_th,s - I_i, which
closes (*) with e_auger = 0).  The script asserts that no such vacancy lies
above the double-ionization threshold I_i + I_(i+1).

Usage
-----
    python3 inner_shell_relaxation_table.py            write the block
    python3 inner_shell_relaxation_table.py --report   print the derivation
    python3 inner_shell_relaxation_table.py --check [--source FILE]
        compare the block installed in cross_sec.f90 with the tables.
"""

import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
REF = os.path.abspath(os.path.join(ROOT, "..", "references"))
KM_DIR = os.path.join(REF, "kaastra_mewe_1993")
NIST_IE = os.path.join(REF, "nist_asd",
                       "ionization_energies_C_N_O_Na_Mg_Si_S_K_Ca_Fe.csv")
DEFAULT_SOURCE = os.path.join(ROOT, "src", "modules", "functions",
                              "cross_sec.f90")
PARAMS = os.path.join(ROOT, "src", "modules", "init", "parameters.f90")

sys.path.insert(0, HERE)
import metal_photoion_table as mpt  # noqa: E402

# Verner shell -> list of (Kaastra & Mewe vacancy number, weight).
SHELL_MAP = {1: [(1, 1.0)], 2: [(2, 1.0)], 3: [(3, 1.0 / 3.0), (4, 2.0 / 3.0)],
             4: [(5, 1.0)], 5: [(6, 1.0 / 3.0), (7, 2.0 / 3.0)], 6: []}

# Stage of each photo-ionizable ion whose NEXT stage the code also
# photoionizes (so I_(i+1) is the code's own mion_ethr): the neutrals of the
# three-stage elements.  Value: the code name of that next stage.
NEXT_STAGE_IN_CODE = {"CI": "CII", "OI": "OII", "NI": "NII", "MgI": "MgII",
                      "SiI": "SiII", "CaI": "CaII", "FeI": "FeII"}

ELEMENT_NAME = {6: "C", 7: "N", 8: "O", 11: "Na", 12: "Mg", 14: "Si",
                16: "S", 19: "K", 20: "Ca", 26: "Fe"}


def read_km_table2():
    rows = {}
    with open(os.path.join(KM_DIR, "table2")) as fh:
        for line in fh:
            if not line.strip():
                continue
            z = int(line[0:2])
            st = int(line[2:5])
            s = int(line[5:7])
            i_vac = float(line[7:14])
            e_a = float(line[14:21])
            eps = int(line[21:25]) * 1.0e-3
            p = []
            for k in range(10):
                f = line[26 + 5 * k:31 + 5 * k]
                p.append(int(f) if f.strip() else 0)
            rows[(z, st, s)] = dict(I=i_vac, EA=e_a, eps=eps, p=p)
    return rows


def read_km_table3():
    rows = {}
    with open(os.path.join(KM_DIR, "table3")) as fh:
        for line in fh:
            if not line.strip():
                continue
            z = int(line[0:3])
            st = int(line[3:6])
            s = int(line[6:9])
            delta = int(line[9:12])
            il = int(line[12:15])
            e = float(line[15:23])
            om = float(line[23:30])
            rows.setdefault((z, st, s), []).append((delta, il, e, om))
    return rows


def read_nist_ie():
    """{(Z, charge): I [eV]} from the NIST ASD export."""
    out = {}
    with open(NIST_IE) as fh:
        rd = csv.reader(fh)
        next(rd)
        for r in rd:
            if len(r) < 10:
                continue
            c = [x.strip().strip("=").strip('"') for x in r]
            try:
                z = int(c[0])
                q = int(c[2].replace("+", ""))
                v = float(c[9])
            except ValueError:
                continue
            out[(z, q)] = v
    return out


def code_thresholds(source):
    """mph_e_th of cross_sec.f90, symbolic entries resolved."""
    text = open(source).read()
    ptext = open(PARAMS).read()
    inst = mpt.parse_installed(source)
    vals = []
    for t in inst["mph_e_th"]:
        if re.match(r"^[0-9.]+d[+-]?\d+$", t, re.I):
            vals.append(float(t.lower().replace("d", "e")))
        else:
            m = re.search(r"%s\s*=\s*([0-9.]+)d([+-]?\d+)" % t, ptext, re.I)
            if not m:
                raise RuntimeError("cannot resolve threshold %s" % t)
            vals.append(float(m.group(1)) * 10.0 ** int(m.group(2)))
    del text
    return vals


def build(source):
    records = mpt.build(mpt.DEFAULT_REFDIR)
    ph1 = mpt.read_ph1(mpt.DEFAULT_REFDIR)
    t2 = read_km_table2()
    t3 = read_km_table3()
    ie = read_nist_ie()
    eth_code = code_thresholds(source)
    name_index = {r["name"]: i for i, r in enumerate(records)}
    out = []
    for k, r in enumerate(records):
        z, ne = r["z"], r["ne"]
        q = z - ne                 # charge of the absorbing ion
        st = q + 1                 # Kaastra & Mewe's stage label
        i_i = eth_code[k]
        # Potential of the n-th electron removed, n = 1 is I_i itself.
        pots = [i_i]
        for n in range(1, 11):
            if n == 1 and r["name"] in NEXT_STAGE_IN_CODE:
                pots.append(eth_code[name_index[NEXT_STAGE_IN_CODE[r["name"]]]])
            else:
                pots.append(ie.get((z, q + n), float("nan")))
        subs = []
        for x in r["sub"]:
            shell = x["shell"]
            e_th = ph1[(z, ne, shell)][0]
            e_th_v = float(e_th)
            p = [0.0] * 10
            e_fl = 0.0
            auger = False
            detail = []
            vac = SHELL_MAP[shell]
            if not vac or all((z, st, sk) not in t2 for sk, _ in vac):
                # No autoionization: the photoelectron alone; the excited
                # ion radiates E_th,s - I_i.
                if any((z, st, sk) in t2 for sk, _ in vac):
                    raise RuntimeError("%s shell %d: partial Table 2 rows"
                                       % (r["name"], shell))
                p[0] = 1.0
                e_fl = e_th_v - i_i
                km_line = []
                for sk, w in vac:
                    for (d, il, e, om) in t3.get((z, st, sk), []):
                        km_line.append((sk, d, il, e, om))
                detail.append(("radiative", km_line))
                if e_th_v > i_i + pots[1]:
                    raise RuntimeError(
                        "%s shell %d: no Auger row, but E_th = %g is above "
                        "the double-ionization threshold %g"
                        % (r["name"], shell, e_th_v, i_i + pots[1]))
            else:
                auger = True
                for sk, w in vac:
                    row = t2.get((z, st, sk))
                    if row is None:
                        raise RuntimeError("%s shell %d: vacancy %d missing"
                                           % (r["name"], shell, sk))
                    tot = float(sum(row["p"]))
                    for n in range(10):
                        p[n] += w * row["p"][n] / tot
                    ep = sum(om * e for (d, il, e, om)
                             in t3.get((z, st, sk), []))
                    e_fl += w * ep
                    detail.append((sk, w, row, ep, tot))
            # ENERGETICS FIRST.  n electrons can leave only if the vacancy
            # energy covers their potentials, E_th,s >= I_i + ... +
            # I_(i+n-1).  Kaastra & Mewe decide each Auger transition on
            # Lotz energy levels (their Sect. 2.2, and their abstract: the
            # potentials "sometimes lead to the uncertainty whether a given
            # Auger transition is energetically allowed or not"); where
            # their distribution puts weight on an n that the measured
            # potentials close, that weight stays at the largest open n
            # and the vacancy energy it would have spent is radiated.
            n_open = 1
            for n in range(2, 11):
                if e_th_v >= sum(pots[0:n]):
                    n_open = n
            closed = sum(p[n_open:])
            if closed > 0.0 and n_open > 1:
                # A cascade cut short after some Auger electrons would
                # need the energy of the closed branch placed by a rule
                # the tables do not give; none of the seventeen ions has
                # one, and a new ion that does is stopped here.
                raise RuntimeError("%s shell %d: partially closed cascade"
                                   % (r["name"], shell))
            if closed > 0.0:
                p[n_open - 1] += closed
                for n in range(n_open, 10):
                    p[n] = 0.0
                if n_open == 1:
                    # No Auger decay is open at all: the radiative case.
                    auger = False
                    e_fl = e_th_v - i_i
                detail.append(("closed", n_open, closed))
            p_multi = 1.0 - p[0]
            e_ion_multi = 0.0
            for n in range(2, 11):
                if p[n - 1] > 0.0:
                    s = sum(pots[1:n])
                    if s != s:
                        raise RuntimeError("%s: missing NIST potential for "
                                           "%d electrons" % (r["name"], n))
                    e_ion_multi += p[n - 1] * s
            e_auger = e_th_v - i_i - e_ion_multi - e_fl
            if auger and e_auger <= 0.0:
                raise RuntimeError("%s shell %d: negative Auger energy %g"
                                   % (r["name"], shell, e_auger))
            if not auger:
                e_auger = 0.0
            subs.append(dict(shell=shell, e_th=e_th, e_th_v=e_th_v,
                             p=p, p_multi=p_multi, e_auger=e_auger,
                             e_fluor=e_fl, e_ion_multi=e_ion_multi,
                             auger=auger, detail=detail))
        out.append(dict(name=r["name"], z=z, ne=ne, i_i=i_i, pots=pots,
                        sub=subs))
    return out


def num(x):
    """A derived number as a Fortran double literal, 9 significant digits."""
    if x == 0.0:
        return "0.0d0"
    s = "%.8e" % x
    m, e = s.split("e")
    return "%sd%d" % (m, int(e))


def field_value(kind, x):
    if kind == "mph_sub_e_th":
        return mpt.fortran(x["e_th"])
    return num(x[{"mph_sub_p_multi": "p_multi",
                  "mph_sub_e_auger": "e_auger",
                  "mph_sub_e_fluor": "e_fluor",
                  "mph_sub_e_ion_multi": "e_ion_multi"}[kind]])


ARRAYS = [
    ("mph_sub_e_th", "Threshold of the subshell itself, E_th,s [eV] "
     "(phfit2.f PH1, field 1)."),
    ("mph_sub_p_multi", "Probability that the vacancy autoionizes, i.e. "
     "that two or more electrons leave (Kaastra & Mewe 1993, Table 2)."),
    ("mph_sub_e_auger", "Mean kinetic energy of the Auger electrons of one "
     "absorption [eV], the remainder of (*)."),
    ("mph_sub_e_fluor", "Mean energy leaving as photons [eV] (Kaastra & "
     "Mewe 1993 Table 3, sum of omega E; E_th,s - I_i for a vacancy that "
     "does not autoionize)."),
    ("mph_sub_e_ion_multi", "Mean ionization energy of the electrons beyond "
     "the first [eV] (NIST ASD; mion_ethr for a stage the code "
     "photoionizes)."),
]


def emit(recs):
    nmax = max(len(r["sub"]) for r in recs)
    out = []
    for name, comment in ARRAYS:
        words, line = comment.split(), "      !"
        for w in words:
            if len(line) + 1 + len(w) > 76:
                out.append(line)
                line = "      !"
            line += " " + w
        out.append(line)
        cols = [[field_value(name, x) for x in r["sub"]] for r in recs]
        lines = ["      real*8, parameter :: %s(n_metal_subshell,"
                 "n_metal_photo) = &" % name, "           reshape( [ &"]
        for i, col in enumerate(cols):
            vals = col + ["0.0d0"] * (nmax - len(col))
            tail = ", &" if i + 1 < len(cols) else " ], &"
            lines.append("           " + ", ".join("%-14s" % v for v in vals)
                         .rstrip() + tail + "   ! " + recs[i]["name"])
        lines.append("           [n_metal_subshell,n_metal_photo] )")
        out += lines
    return out


def report(recs):
    kmn = {1: "K", 2: "L1", 3: "L2", 4: "L3", 5: "M1", 6: "M2", 7: "M3"}
    for r in recs:
        print("%s  I_i = %.4f  I_(i+1..) = %s" % (
            r["name"], r["i_i"],
            " ".join("%.3f" % v for v in r["pots"][1:5])))
        for x in r["sub"]:
            ps = " ".join("%.4f" % v for v in x["p"][:6])
            print("   %s E_th=%8.2f p_multi=%.4f P(n=1..6)=%s e_ion_multi="
                  "%8.3f e_fluor=%7.3f e_auger=%8.3f"
                  % (mpt.SHELL_NAME[x["shell"]], x["e_th_v"], x["p_multi"],
                     ps, x["e_ion_multi"], x["e_fluor"], x["e_auger"]))
            for d in x["detail"]:
                if d[0] == "closed":
                    print("        ENERGETICALLY CLOSED: weight %.4f of the "
                          "Table 2 distribution needs more than the vacancy "
                          "energy; held at n = %d" % (d[2], d[1]))
                elif d[0] == "radiative":
                    for (sk, dd, il, e, om) in d[1]:
                        print("        radiative: KM93 Table 3 %s delta=%d "
                              "line %d E=%.1f omega=%.4f (used: E_th,s - I_i"
                              " = %.2f)" % (kmn[sk], dd, il, e, om,
                                            x["e_fluor"]))
                else:
                    sk, w, row, ep, tot = d
                    e_a_prime = row["EA"] - row["eps"] * row["I"]
                    print("        KM93 %s weight %.3f: I=%.1f E_A=%.1f "
                          "eps=%.3f E'_A=E_A-eps I=%.1f sum(omega E)=%.3f "
                          "row sum=%d" % (kmn[sk], w, row["I"], row["EA"],
                                          row["eps"], e_a_prime, ep, tot))


def parse_installed(source):
    text = open(source).read()
    text = re.sub(r"!.*", "", text)
    text = re.sub(r"&\s*\n\s*", " ", text)
    found = {}
    for m in re.finditer(
            r"real\*8,\s*parameter\s*::\s*(mph_sub_\w+)"
            r"\(n_metal_subshell,n_metal_photo\)\s*=\s*reshape\(\s*\[([^\]]*)\]",
            text):
        found[m.group(1)] = [t.strip() for t in m.group(2).split(",")
                             if t.strip()]
    return found


def to_float(t):
    return float(t.lower().replace("d", "e"))


def check(source):
    recs = build(source)
    inst = parse_installed(source)
    nmax = max(len(r["sub"]) for r in recs)
    fails = 0
    for name, _ in ARRAYS:
        want = []
        for r in recs:
            want += [field_value(name, x) for x in r["sub"]]
            want += ["0.0d0"] * (nmax - len(r["sub"]))
        got = inst.get(name)
        if got is None or len(got) != len(want):
            print("FAIL inner_shell_relaxation_%s missing or wrong size" % name)
            fails += 1
            continue
        bad = [i for i in range(len(want))
               if abs(to_float(got[i]) - to_float(want[i]))
               > 1.0e-8 * max(abs(to_float(want[i])), 1.0e-30)]
        if bad:
            for i in bad:
                print("FAIL inner_shell_relaxation_%s[%d] measured=%s "
                      "reference=%s tol=1e-8" % (name, i, got[i], want[i]))
            fails += len(bad)
        else:
            print("PASS inner_shell_relaxation_%s measured=%d_fields_match "
                  "reference=%d_fields tol=1e-8" % (name, len(want),
                                                     len(want)))
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=DEFAULT_SOURCE)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.check:
        sys.exit(1 if check(a.source) else 0)
    recs = build(a.source)
    if a.report:
        report(recs)
    else:
        print("\n".join(emit(recs)))


if __name__ == "__main__":
    main()
