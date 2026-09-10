#!/usr/bin/env python3
"""Ownership of emitted and absorbed radiation, and the photoionization sum.

Origin: docs/audit_20260905/review3_ledger_checks.py, whose last two prints
show what the two ownership rules give for the same two cells.  This is
arithmetic over the convention itself, so it holds no production symbol and
is green by construction; it is kept so that the convention the Fortran
tests are written against is pinned in one readable place.

The convention (docs/development_plan_20260905_rev3.md section 4.2 item 3):
"escaping" means leaving the local material system of the emitting cell, not
leaving the atmosphere.  Emission is a loss of the emitting cell, absorption
is a gain of the absorbing cell.  A photon emitted by cell A and fully
absorbed by cell B therefore changes the material energy of the pair by
zero.  The rule that subtracts only the emission that leaves the domain
creates that photon's energy out of nothing.

Verdict lines follow the convention of the directory:
    PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
"""

import sys

TOL = 1.0e-12
failures = 0


def check(name, measured, reference, tol=TOL):
    global failures
    ok = abs(measured - reference) <= tol
    if not ok:
        failures += 1
    print("%s %s measured=%.15e reference=%.15e tol=%.15e absolute"
          % ("PASS" if ok else "FAIL", name, measured, reference, tol))


# One photon of 10 eV emitted by cell A and fully absorbed by cell B, with
# none of it leaving the two-cell system.
photon = 10.0
emitting_cell = -photon
absorbing_cell = +photon
check("two_cell_radiation_exchange_material_sum",
      emitting_cell + absorbing_cell, 0.0)

# The same pair under the rule that counts only domain-escaping emission as
# a loss: the emitting cell loses nothing, the absorbing cell still gains.
# The rule is not the convention; this line records what it would create.
print("note: material energy created by the escaping-emission-only rule "
      "[eV] = %.15e" % (0.0 + absorbing_cell))

# Isolated H photoionization: thermal increment plus chemical increment is
# the absorbed photon energy, for any ionization potential.  13.6 eV is the
# threshold parameters.f90 carries, 13.598434599 eV the CODATA value the
# audit probe used; the identity is independent of both.
for chi in (13.6, 13.598434599):
    photon_energy = 20.0
    thermal = photon_energy - chi
    check("h_photoionization_event_energy_sum_chi_%.9f" % chi,
          thermal + chi, photon_energy)

if failures:
    print("radiation_exchange_ledger: %d assertion(s) failed" % failures)
    sys.exit(1)
print("radiation_exchange_ledger: all assertions passed")
