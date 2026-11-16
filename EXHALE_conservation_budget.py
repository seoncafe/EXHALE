#!/usr/bin/env python3
"""Rebuild the discrete conservation balance of an EXHALE state from the
terms the assembly exported, and report where it fails to close.

It reads output/conservation_budget_<nnnn>.txt, written by the module
conservation_budget when EXHALE_CONSERVATION_BUDGET is set, and forms every
row again from the face fluxes, the face areas, the cell volume, the
potential and the sources.  It shares no arithmetic with the Fortran
assembly: it reads the terms and nothing else, and in particular it never
sums the exported residual and calls that a reconstruction.

Two checks are reported separately and must not be confused.

  1. ASSEMBLY CONSISTENCY.  Does the row rebuilt from the exported terms
     equal the row the assembly returned?  Its tolerance comes from
     arithmetic and from the serialization, not from any physical
     criterion, and it is reported in units of the last bit of the largest
     term the row differences.

  2. STATIONARITY.  Does that row meet the certification's own measure?
     Its tolerance is the certification's.  A state can pass the first
     check and fail this one by many decades, and that is the ordinary
     situation for a checkpoint that was never certified.

WHAT THE THREE ROWS ARE.  Total mass, radial momentum and total energy.  No
species appears in them, so no nucleus-count weighting arises here: the mass
row already carries every nucleus with its own mass.  An elemental or
carrier budget, where H2 alone is not conserved through a dissociation front
and species have to be weighted by their nucleus counts, is a different
export and is not read by this tool.

Usage:

    python3 EXHALE_conservation_budget.py <budget file> [more files ...]
                                          [--top N] [--subdomain J] [--csv]

Standard library only.
"""

import sys

# The radius at which the certification stops reporting and starts gating
# the species and carrier rows (cert_regime_wind_r, certification.f90).  It
# is used here only to place one subdomain boundary, so that the flux
# crossing it is reported on its own.
CERT_REGIME_WIND_R = 1.20

# The certification's tolerances for the three hydrodynamic rows
# (cert_tol_momentum, cert_tol_energy, certification.f90), and the FIXED
# part of the continuity row's, which is a floor: the binary raises it cell
# by cell to the cell's own rounding floor (cert_tol_mass_at), a quantity
# that needs the caloric sound speed and is therefore not rebuilt here.  The
# value below is the smallest tolerance the continuity row can be held to,
# so a cell that fails against it may still be inside the binary's own,
# larger, tolerance.
CERT_TOL = {"mass": 3.0e-12, "momentum": 1.0e-8, "energy": 1.0e-6}

EPS = sys.float_info.epsilon

ROWS = ("mass", "momentum", "energy")


def interdiffusion_divergence(cell):
    """The divergence of the interdiffusion enthalpy flux of one cell, the
    Sidf column of a schema-2 file, ADDED to the energy row by the assembly:
    R_energy = dF_energy - S_energy - (heat - cool) - Sene + Sidf.  A
    schema-1 file was written by a code that had no such term, so it is
    zero there."""
    return cell.get("Sidf", 0.0)


# --------------------------------------------------------------------- #
# reading


class Budget(object):
    """One exported assembly: its header statements and its cells."""

    def __init__(self, path):
        self.path = path
        self.comments = []
        self.meta = {}
        self.columns = None
        self.cells = []
        self._read()

    def _read(self):
        with open(self.path, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    self.comments.append(line)
                    self._header(line)
                    continue
                if not line.strip():
                    continue
                if self.columns is None:
                    raise ValueError(
                        "%s: a data row stands before the '# columns' line"
                        % self.path)
                fields = line.split()
                if len(fields) != len(self.columns):
                    raise ValueError(
                        "%s: %d fields against %d column names"
                        % (self.path, len(fields), len(self.columns)))
                cell = {}
                for name, text in zip(self.columns, fields):
                    if name in ("j", "physical"):
                        cell[name] = int(text)
                    else:
                        cell[name] = float(text)
                self.cells.append(cell)
        for want in ("momentum_branch", "rows_kind", "ieq_sweep_state_kind",
                     "reconstruction"):
            if want not in self.meta:
                raise ValueError("%s: the header does not state %s"
                                 % (self.path, want))

    def _header(self, line):
        body = line[1:].strip()
        if body.startswith("columns "):
            self.columns = body.split()[1:]
            return
        for tag in ("assembly ", "flags "):
            if body.startswith(tag):
                for token in body[len(tag):].split():
                    if "=" in token:
                        key, value = token.split("=", 1)
                        self.meta[key] = value
                return
        if body.startswith("EXHALE conservation_budget schema"):
            self.meta["schema"] = body.split()[-1]

    def physical(self):
        return [c for c in self.cells if c["physical"] == 1]

    def identity_lines(self):
        keep = ("# EXHALE conservation_budget schema",
                "# export index", "# provenance:", "# coupling:",
                "# assembly ", "# flags ", "# rec_method", "# state ",
                "# rows ")
        return [c for c in self.comments if c.startswith(keep)]


# --------------------------------------------------------------------- #
# the rows, rebuilt


def rebuilt_rows(cell, branch, wants_transport):
    """The three residual rows of one cell, formed from its exported terms.

    Returns the three values and, beside each, the largest absolute term
    that entered it, which is what the arithmetic floor of the difference is
    measured against.
    """
    a_lo, a_hi, vol, dr = (cell["area_lo"], cell["area_hi"],
                           cell["volume"], cell["dr"])
    fm_lo, fm_hi = cell["face_mass_lo"], cell["face_mass_hi"]
    fp_lo, fp_hi = cell["face_momentum_lo"], cell["face_momentum_hi"]
    fe_lo, fe_hi = cell["face_energy_lo"], cell["face_energy_hi"]
    p_lo, p_hi = cell["face_p_lo"], cell["face_p_hi"]
    q_dn, q_up = cell["face_q_dn_lo"], cell["face_q_up_hi"]
    phi_lo, phi_hi, phi_c = (cell["phi_face_lo"], cell["phi_face_hi"],
                             cell["phi_cell"])
    smom = cell["Smom"] if wants_transport else 0.0
    sene = cell["Sene"] if wants_transport else 0.0

    out = {}

    # mass
    terms = [a_hi * fm_hi / vol, -a_lo * fm_lo / vol, -cell["S_mass"]]
    out["mass"] = ((a_hi * fm_hi - a_lo * fm_lo) / vol - cell["S_mass"],
                   max(abs(t) for t in terms))

    # momentum, in the branch the header names
    if branch == "plm_ordinary":
        dfm = (a_hi * fp_hi - a_lo * fp_lo) / vol
        terms = [a_hi * fp_hi / vol, a_lo * fp_lo / vol]
    elif branch == "weno3_ordinary":
        dfm = (a_hi * fp_hi - a_lo * fp_lo) / vol + (p_hi - p_lo) / dr
        terms = [a_hi * fp_hi / vol, a_lo * fp_lo / vol,
                 p_hi / dr, p_lo / dr]
    elif branch == "plm_well_balanced":
        dfm = (a_hi * (fp_hi + q_up) - a_lo * (fp_lo + q_dn)) / vol
        terms = [a_hi * (fp_hi + q_up) / vol, a_lo * (fp_lo + q_dn) / vol]
    elif branch == "weno3_well_balanced":
        dfm = ((a_hi * fp_hi - a_lo * fp_lo) / vol + (q_up - q_dn) / dr)
        terms = [a_hi * fp_hi / vol, a_lo * fp_lo / vol,
                 q_up / dr, q_dn / dr]
    else:
        raise ValueError("unrecognized momentum branch %r" % branch)
    terms = terms + [cell["S_momentum"], smom]
    out["momentum"] = (dfm - cell["S_momentum"] - smom,
                       max(abs(t) for t in terms))

    # energy: the gravitational work rides on the mass face flux INSIDE the
    # flux difference, and the assembly divides the sum by the volume, so
    # the sum is formed here before the division as well.
    grav = (a_hi * fm_hi * (phi_hi - phi_c)
            - a_lo * fm_lo * (phi_lo - phi_c))
    dfe = (a_hi * fe_hi - a_lo * fe_lo + grav) / vol
    sidf = interdiffusion_divergence(cell)
    terms = [a_hi * fe_hi / vol, a_lo * fe_lo / vol, grav / vol,
             cell["S_energy"], cell["heat"], cell["cool"], sene, sidf]
    out["energy"] = (dfe - cell["S_energy"] - (cell["heat"] - cell["cool"])
                     - sene + sidf,
                     max(abs(t) for t in terms))
    out["_grav_over_volume"] = grav / vol
    out["_dF"] = {"mass": (a_hi * fm_hi - a_lo * fm_lo) / vol,
                  "momentum": dfm, "energy": dfe}
    return out


def row_scales(cell, wants_transport):
    """The certification's own scale of each row, from the exported terms.

    mass      max(|A F|_hi, |A F|_lo)/V           (mass_flux_row_scale)
    momentum  max(|ram|, |pressure|, |gravity|, |Smom|)
    energy    max(|dF_E|, |S_E|, heat, cool, |Sene|, |Sidf|)

    all of steady_residual.f90.  The momentum and energy scales read the
    exported production attribution, which is what the binary's own measure
    reads, so a measure formed here is comparable with the one the run
    printed.
    """
    smom = cell["Smom"] if wants_transport else 0.0
    sene = cell["Sene"] if wants_transport else 0.0
    tiny = sys.float_info.min
    s_mass = max(abs(cell["face_mass_hi"] * cell["area_hi"]),
                 abs(cell["face_mass_lo"] * cell["area_lo"])) / cell["volume"]
    s_mom = max(abs(cell["momentum_ram"]), abs(cell["momentum_pressure"]),
                abs(cell["momentum_gravity"]), abs(smom))
    s_ene = max(abs(cell["dF_energy"]), abs(cell["S_energy"]),
                abs(cell["heat"]), abs(cell["cool"]), abs(sene),
                abs(interdiffusion_divergence(cell)))
    return {"mass": max(s_mass, tiny), "momentum": max(s_mom, tiny),
            "energy": max(s_ene, tiny)}


def energy_with_potential_defect(cell, wants_transport):
    """The identity of section 13.4 of PLAN_20260920_rev9, at one cell:

        R_E + phi_c R_mass
            = [A_R (F_E,R + phi_R F_m,R)
             - A_L (F_E,L + phi_L F_m,L)]/V - phi_c S_mass - Q

    with Q = heat - cool + Sene - Sidf.  The term phi_c R_mass is KEPT: these
    checkpoints are not stationary and their mass rows are not small, so
    dropping it would impose a stationary mass identity on a state that has
    none.  Returns the defect and the largest term it is a difference of.
    """
    a_lo, a_hi, vol = cell["area_lo"], cell["area_hi"], cell["volume"]
    phi_lo, phi_hi, phi_c = (cell["phi_face_lo"], cell["phi_face_hi"],
                             cell["phi_cell"])
    sene = cell["Sene"] if wants_transport else 0.0
    left_hi = a_hi * (cell["face_energy_hi"]
                      + phi_hi * cell["face_mass_hi"]) / vol
    left_lo = a_lo * (cell["face_energy_lo"]
                      + phi_lo * cell["face_mass_lo"]) / vol
    q = cell["heat"] - cell["cool"] + sene - interdiffusion_divergence(cell)
    right = left_hi - left_lo - phi_c * cell["S_mass"] - q
    left = cell["R_energy"] + phi_c * cell["R_mass"]
    terms = [abs(left_hi), abs(left_lo), abs(q), abs(cell["R_energy"]),
             abs(phi_c * cell["R_mass"])]
    return left - right, max(terms)


# --------------------------------------------------------------------- #
# statistics


class Defects(object):
    """Signed sum, sum of magnitudes and largest magnitude of a set of
    local defects, because a signed sum alone hides opposing errors."""

    def __init__(self, label):
        self.label = label
        self.signed = 0.0
        self.absolute = 0.0
        self.largest = 0.0
        self.largest_at = None
        self.largest_ulp = 0.0
        self.largest_ulp_at = None
        self.n = 0

    def add(self, value, scale, where):
        self.n += 1
        self.signed += value
        self.absolute += abs(value)
        if abs(value) > self.largest:
            self.largest = abs(value)
            self.largest_at = where
        ulp = abs(value) / (EPS * scale) if scale > 0.0 else 0.0
        if ulp > self.largest_ulp:
            self.largest_ulp = ulp
            self.largest_ulp_at = where

    def line(self):
        return ("%-34s n=%4d  signed %+.6E  sum|.| %.6E  max|.| %.6E"
                " at cell %s  max ulp %.3E at cell %s"
                % (self.label, self.n, self.signed, self.absolute,
                   self.largest, self.largest_at, self.largest_ulp,
                   self.largest_ulp_at))


# --------------------------------------------------------------------- #
# the report


def audit(budget, top, extra_subdomains, csv):
    branch = budget.meta["momentum_branch"]
    wants_transport = budget.meta.get("transport_active", "F") == "T"
    cells = budget.physical()
    if not cells:
        raise ValueError("%s: no physical cell" % budget.path)

    print("=" * 78)
    print("CONSERVATION BUDGET AUDIT  %s" % budget.path)
    print("=" * 78)
    for line in budget.identity_lines():
        print(line)
    print("")
    print("physical cells        %d, j = %d to %d"
          % (len(cells), cells[0]["j"], cells[-1]["j"]))
    print("momentum branch       %s" % branch)
    print("transport subtracted  %s"
          % ("yes" if wants_transport else "no (both sources are zero)"))
    print("")

    rebuilt = {}
    scales = {}
    for c in cells:
        rebuilt[c["j"]] = rebuilt_rows(c, branch, wants_transport)
        scales[c["j"]] = row_scales(c, wants_transport)

    # ---- check 1 a: the rows -------------------------------------- #
    print("-" * 78)
    print("CHECK 1  ASSEMBLY CONSISTENCY")
    print("  the row rebuilt from the exported faces, geometry, potential")
    print("  and sources, against the row the assembly returned.")
    print("  'ulp' is the defect divided by eps times the largest term the")
    print("  row differences: 1 means one last bit of that term.")
    print("-" * 78)
    consistency = {}
    for row in ROWS:
        d = Defects("row %s" % row)
        for c in cells:
            value, largest = rebuilt[c["j"]][row]
            d.add(value - c["R_%s" % row], largest, c["j"])
        consistency[row] = d
        print(d.line())

    print("")
    print("  the arithmetic floor of each rebuild: eps times the largest")
    print("  term the row differences, which is the smallest defect a")
    print("  binary64 difference of those terms can resolve.")
    for row in ROWS:
        floors = [EPS * rebuilt[c["j"]][row][1] for c in cells]
        floors_sorted = sorted(floors)
        median = floors_sorted[len(floors_sorted) // 2]
        print("    row %-9s floor  min %.3E  median %.3E  max %.3E"
              % (row, floors_sorted[0], median, floors_sorted[-1]))
    print("")

    d = Defects("dF_energy gravitational work")
    for c in cells:
        d.add(rebuilt[c["j"]]["_grav_over_volume"]
              - c["grav_work_over_volume"],
              abs(c["grav_work_over_volume"]), c["j"])
    print(d.line())

    for row in ROWS:
        d = Defects("flux difference dF_%s" % row)
        for c in cells:
            value = rebuilt[c["j"]]["_dF"][row]
            _, largest = rebuilt[c["j"]][row]
            d.add(value - c["dF_%s" % row], largest, c["j"])
        print(d.line())

    # ---- check 1 b: the energy identity with the potential --------- #
    print("")
    print("  the identity of PLAN_20260920_rev9 section 13.4, kept BESIDE")
    print("  the gas-energy row above and never in place of it:")
    print("    R_E + phi_cell R_mass = [A_R (F_E + phi F_m)_R")
    print("                            - A_L (F_E + phi F_m)_L]/V - Q")
    d_id = Defects("energy + potential identity")
    for c in cells:
        value, largest = energy_with_potential_defect(c, wants_transport)
        d_id.add(value, largest, c["j"])
    print(d_id.line())

    # ---- check 1 c: internal faces --------------------------------- #
    print("")
    print("  internal-face agreement of the SHARED CONSERVATIVE fluxes.")
    print("  It does not cover the equilibrium pressure departures: the two")
    print("  values of one face are measured from two different cells'")
    print("  equilibria and differ by construction (Num_Fluxes.f90).")
    for name, hi, lo in (("area", "area_hi", "area_lo"),
                         ("mass flux", "face_mass_hi", "face_mass_lo"),
                         ("momentum flux", "face_momentum_hi",
                          "face_momentum_lo"),
                         ("energy flux", "face_energy_hi", "face_energy_lo"),
                         ("potential", "phi_face_hi", "phi_face_lo")):
        d = Defects("shared face %s" % name)
        for a, b in zip(cells[:-1], cells[1:]):
            d.add(a[hi] - b[lo], max(abs(a[hi]), abs(b[lo])), a["j"])
        print(d.line())

    # ---- check 1 d: telescoping over subdomains -------------------- #
    print("")
    print("  volume-weighted telescoping over subdomains: the sum of V R")
    print("  over a window against the flux the window's two end faces")
    print("  carry.  For the mass row and for the energy row augmented by")
    print("  the potential energy transport, both of which are exact sums")
    print("  of one conservative flux.")
    gate = None
    for c in cells:
        if c["r_cell"] >= CERT_REGIME_WIND_R:
            gate = c["j"]
            break
    windows = [("the whole physical column", cells[0]["j"], cells[-1]["j"])]
    if gate is not None:
        windows.append(("below the certification gate r = %.2f"
                        % CERT_REGIME_WIND_R, cells[0]["j"], gate - 1))
        windows.append(("above the certification gate r = %.2f"
                        % CERT_REGIME_WIND_R, gate, cells[-1]["j"]))
    worst = {}
    for row in ROWS:
        j = max(cells, key=lambda c: abs(c["R_%s" % row]))["j"]
        worst[row] = j
        windows.append(("to the worst %s cell, j = %d" % (row, j),
                        cells[0]["j"], j))
    for j in extra_subdomains:
        windows.append(("to cell j = %d" % j, cells[0]["j"], j))

    for label, ja, jb in windows:
        window = [c for c in cells if ja <= c["j"] <= jb]
        if not window:
            continue
        sum_mass = sum(c["volume"] * c["R_mass"] for c in window)
        flux_mass = (window[-1]["area_hi"] * window[-1]["face_mass_hi"]
                     - window[0]["area_lo"] * window[0]["face_mass_lo"])
        sum_ene = sum(c["volume"] * (c["R_energy"]
                                     + c["phi_cell"] * c["R_mass"])
                      for c in window)
        flux_ene = ((window[-1]["area_hi"]
                     * (window[-1]["face_energy_hi"]
                        + window[-1]["phi_face_hi"]
                        * window[-1]["face_mass_hi"]))
                    - (window[0]["area_lo"]
                       * (window[0]["face_energy_lo"]
                          + window[0]["phi_face_lo"]
                          * window[0]["face_mass_lo"])))
        src_ene = sum(c["volume"] * (c["heat"] - c["cool"]
                                     + (c["Sene"] if wants_transport else 0.0)
                                     - interdiffusion_divergence(c)
                                     + c["phi_cell"] * c["S_mass"])
                      for c in window)
        print("  %-46s cells %4d..%4d" % (label, ja, jb))
        print("      mass    sum V R  %+.10E   faces %+.10E"
              "   defect %+.3E"
              % (sum_mass, flux_mass, sum_mass - flux_mass))
        print("      energy  sum V(R+phi R_m) %+.10E   faces - sources"
              " %+.10E   defect %+.3E"
              % (sum_ene, flux_ene - src_ene,
                 sum_ene - (flux_ene - src_ene)))

    # ---- check 2: stationarity ------------------------------------- #
    print("")
    print("-" * 78)
    print("CHECK 2  STATIONARITY")
    print("  the reconstructed row against the certification's own measure,")
    print("  |R| divided by the largest term of that row, and against the")
    print("  certification's tolerance.  The continuity row's tolerance is")
    print("  the FIXED FLOOR only: the binary raises it cell by cell to the")
    print("  cell's rounding floor, which needs the caloric sound speed and")
    print("  is not rebuilt here, so a refusal below is a lower bound.")
    print("-" * 78)
    for row in ROWS:
        worst_j, worst_m = None, 0.0
        for c in cells:
            value, _ = rebuilt[c["j"]][row]
            m = abs(value) / scales[c["j"]][row]
            if m > worst_m:
                worst_m, worst_j = m, c["j"]
        exported_m = max(abs(c["R_%s" % row]) / scales[c["j"]][row]
                         for c in cells)
        verdict = "MEETS" if worst_m <= CERT_TOL[row] else "REFUSES"
        print("  %-9s row measure %.6E at cell %-4d  tolerance %.1E"
              "   %s" % (row, worst_m, worst_j, CERT_TOL[row], verdict))
        print("            the same measure of the EXPORTED row: %.6E"
              % exported_m)

    # ---- the worst cells ------------------------------------------- #
    print("")
    print("-" * 78)
    print("THE LARGEST LOCAL ASSEMBLY DEFECTS, %d cells per row" % top)
    print("-" * 78)
    for row in ROWS:
        print("  row %s" % row)
        ranked = sorted(
            cells,
            key=lambda c: abs(rebuilt[c["j"]][row][0] - c["R_%s" % row]),
            reverse=True)[:top]
        for c in ranked:
            value, largest = rebuilt[c["j"]][row]
            defect = value - c["R_%s" % row]
            ulp = abs(defect) / (EPS * largest) if largest > 0.0 else 0.0
            print("    j %4d  r %10.6f  R exported %+.10E  rebuilt %+.10E"
                  "  defect %+.3E  %.2f ulp"
                  % (c["j"], c["r_cell"], c["R_%s" % row], value,
                     defect, ulp))

    if csv:
        print("")
        print("# csv j,r,row,R_exported,R_rebuilt,defect,ulp,measure")
        for c in cells:
            for row in ROWS:
                value, largest = rebuilt[c["j"]][row]
                defect = value - c["R_%s" % row]
                ulp = abs(defect) / (EPS * largest) if largest > 0.0 else 0.0
                print("%d,%.17g,%s,%.17g,%.17g,%.17g,%.6g,%.6g"
                      % (c["j"], c["r_cell"], row, c["R_%s" % row], value,
                         defect, ulp, abs(value) / scales[c["j"]][row]))

    return consistency, d_id


def main(argv):
    paths, top, subdomains, csv = [], 6, [], False
    i = 1
    while i < len(argv):
        a = argv[i]
        if a == "--top":
            i += 1
            top = int(argv[i])
        elif a == "--subdomain":
            i += 1
            subdomains.append(int(argv[i]))
        elif a == "--csv":
            csv = True
        elif a in ("-h", "--help"):
            print(__doc__)
            return 0
        else:
            paths.append(a)
        i += 1
    if not paths:
        print(__doc__)
        return 2
    for path in paths:
        audit(Budget(path), top, subdomains, csv)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
