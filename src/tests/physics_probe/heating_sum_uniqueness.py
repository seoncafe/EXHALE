#!/usr/bin/env python3
"""One assembly of the volumetric heating rate in the production source.

The heating of a cell is a sum of physically distinct deposits: the
photoionization of each absorber, the two H(n=2) channels, the photoelectrons
of the He recombination radiation, the three He(2^3S) collision branches, the
two halves of the Lyman-Werner absorption, the collisional reactions of the
H2/He network, the excess energy of the FUV photolysis events and the
collisional oxygen channels.  It is assembled once, in
heating_of_composition (src/modules/radiation/util_ion_eq.f90), and its three
consumers -- the ionization sweep, the Heating_breakdown.txt dump and the
advection-corrected post-process -- call that one routine.

A second copy of the sum is not wrong the day it is written; it is wrong the
day one copy gains a deposit the other does not.  That is what happened: the
associative He(2^3S) branch and the collisional oxygen channels reached the
energy equation through the sweep's copy while the breakdown dump's copy and
the post-process's copy went on without them, and the two named columns of
one state disagreed by 6.9e-3 in the hot-Uranus gate.

WHAT IS COMPARED.  Every .f90 file under src/modules is scanned, with
comments and character strings removed, for the production entities that
FORM a heating deposit.  Each may be named in the module that defines it and
in the module that holds the one assembly, and nowhere else.

REFERENCE.  Zero files outside those two per name.  Tolerance zero.

The scan does not and cannot catch an assembly written out of literal
arithmetic; what it catches is the way all three copies were actually
written, by calling the routines and reading the energies below.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MODULES = os.path.join(ROOT, "src", "modules")

ASSEMBLY = os.path.join("src", "modules", "radiation", "util_ion_eq.f90")

# name -> the module that defines it.  The one assembly is allowed for all
# of them; nothing else is.
DEPOSITS = {
    # the photoionization deposit of a composition
    "photoheating_of_composition": ASSEMBLY,
    # the collisional H2/He network and the collisional oxygen network
    "molecular_chemical_heating": os.path.join(
        "src", "modules", "lower_atmosphere", "molecular_reaction_heat.f90"),
    "oxygen_chemical_heating": os.path.join(
        "src", "modules", "lower_atmosphere", "molecular_reaction_heat.f90"),
    # the excess energy of one FUV photolysis event
    "heat_per_water_dissociation": os.path.join(
        "src", "modules", "lower_atmosphere", "water_photolysis.f90"),
    "heat_per_hydroxyl_dissociation": os.path.join(
        "src", "modules", "lower_atmosphere", "water_photolysis.f90"),
    # the two halves of one Lyman-Werner absorption
    "e_lw_fragment_erg": os.path.join(
        "src", "modules", "lower_atmosphere", "lyman_werner.f90"),
    "h2_energy_per_bound_fluorescence_erg": os.path.join(
        "src", "modules", "lower_atmosphere",
        "h2_vibrational_relaxation.f90"),
    # the fragment kinetic energy of one CO photodissociation
    "heat_per_co_dissociation": os.path.join(
        "src", "modules", "lower_atmosphere", "co_photodissociation.f90"),
}


# NO EXCEPTIONS.  There was one until 2026-09-06: write_output.f90 re-formed
# the energy of one photolysis event, one Lyman-Werner fragment pair and its
# fluorescence, for the column-integrated band ledger of output/FUV_bands.txt
# and the heating column of output/Lyman_Werner.txt.  Those are aggregates
# over the grid at fixed band rather than rewritings of the assembly's sum
# over bands at fixed cell, so neither was a second heating total; but the
# per-event ENERGIES were written a second time there and could drift from
# the energy equation exactly as the three copies of the heating sum did.
# The ledger now lives in utils_ion_eq::fuv_band_absorption_ledger, beside
# the assembly and reading the same imports, and the Lyman-Werner heating
# column is the assembly's own channel array.  The list is empty and a new
# copy anywhere fails the scan.
EXCEPTIONS = {}


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
    patterns = {n: re.compile(r"(?<![A-Za-z0-9_])" + n + r"(?![A-Za-z0-9_])",
                              re.IGNORECASE)
                for n in DEPOSITS}
    failures = []
    for base, dirs, files in os.walk(MODULES):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for name in sorted(files):
            if not name.endswith(".f90"):
                continue
            path = os.path.join(base, name)
            rel = os.path.relpath(path, ROOT)
            with open(path, "r", errors="replace") as fh:
                lines = code_only(fh.read())
            for deposit, definer in DEPOSITS.items():
                if rel in (ASSEMBLY, definer):
                    continue
                if deposit in EXCEPTIONS.get(rel, ()):
                    continue
                pat = patterns[deposit]
                for n, line in enumerate(lines, 1):
                    if pat.search(line):
                        failures.append("%s:%d: %s -- %s"
                                        % (rel, n, deposit, line.strip()))

    name = "heating_sum_uniqueness"
    if failures:
        print("FAIL %s measured=%d reference=0 tol=0 heating deposits formed "
              "outside the one assembly" % (name, len(failures)))
        for f in failures:
            print("  " + f)
        return 1
    print("PASS %s measured=0 reference=0 tol=0 every heating deposit is "
          "formed in one place" % name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
