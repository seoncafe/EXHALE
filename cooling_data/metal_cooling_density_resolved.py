"""Density-resolved C/N/O line cooling above the ground term, from CHIANTI.

WHAT THIS BUILDS.  The closed-form coefficients in Cool_coeff.f90 carry the
line cooling of each metal ion in the CORONAL limit: every collisional
excitation is assumed to radiate.  That is the n_e -> 0 limit of statistical
equilibrium and it overstates the cooling of any line whose critical density
is at or below the electron density of the gas.  This script tabulates the
quantity that replaces it, the statistical equilibrium of the same channels at
the gas's own electron density,

    Lambda_SE,rem(T, n_e)   [erg cm^3 s^-1], on a (log T, log n_e) grid,

the radiative loss of one ion divided by n_e, summed over every transition
whose UPPER level lies OUTSIDE the ground term:

    Lambda_SE,rem(T, n_e) = ( sum_{u not in ground term, l} f_u A_ul E_ul )/n_e

with f_u the fractional level populations of ChiantiPy's `populate`, the full
multilevel statistical equilibrium of the ion under electron collisions and
spontaneous decay.  The ground-term transitions are excluded because
Cool_coeff.f90 solves that term exactly itself, with the neutral-hydrogen
de-excitation channel and the line escape probability that CHIANTI has no way
to carry; leaving them in the table would double count them.

The floor column of the table is the coronal limit of those same channels, so
it is directly comparable with the closed-form fit it replaces; that ratio,
printed to stderr as f_cov, is the measure of the fit's error and is reported
in the memo.  It is NOT unity at low temperature: the fits were built from the
.scups THEORETICAL transition energies, while the level populations use the
observed level energies, and the two differ by up to 14 percent in dE for the
O II 4S* - 2D* excitation, which at 2000 K is a factor 14 in exp(-dE/kT).

METHOD.  Taking the line cooling from the level populations at the local
electron density is the method of MoCHII (~/MoCHII/MoCHII_v1.00,
md/COOLING_LOCAL_NE_PLAN.md, the same author's photoionization code, where it
is the default since 2026-07-26), and the argument for a precomputed
(T, n_e) table rather than a live solve is that document's design point D1.
MoCHII keeps its fits as the carrier and tabulates only the suppression
because its own level models are truncated; here the CHIANTI model is the full
one the fits themselves were made from, so the table carries the cooling
outright and the fits remain as the coronal reference and for the legacy
branch.  The table is emitted as Fortran parameter data, which is how this
repository already carries the Fe II table (cooling_data/fe2_cooling.py).

ATOMIC DATA.  CHIANTI v11.0.2, read through ChiantiPy; set XUVTOP.

GRID.  The axes are the ones Cool_coeff.f90 already declares for the Fe II
table, so no new axis is introduced: log10 T from 3.0 to 5.0 in 0.05 dex (41
points) and log10 n_e from 0 to 14 in 0.5 dex (29 points).  The floor
n_e = 1 cm^-3 is four decades below the lowest critical density of any
transition in the sum (the metastable terms of C I, N I and O II, ~1e4 cm^-3),
the fine-structure lines having been excluded.

RUN:  XUVTOP=<dbase> python3 metal_cooling_density_resolved.py [outfile]
      The default outfile is metal_cooling_density_resolved_tables.txt next to
      this script; its contents are pasted into Cool_coeff.f90.
"""

import sys
import numpy as np
import ChiantiPy.core as ch

H_TIMES_C = 1.9864458571489286e-16      # erg cm

# log10 T and log10 n_e axes of Cool_coeff.f90 (cool_logT, cool_logne)
LOGT = np.linspace(3.0, 5.0, 41)
LOGNE = np.linspace(0.0, 14.0, 29)

# CHIANTI name, the Fortran tag, and the number of levels of the ground term.
# The ground term is the lowest levels of the ground configuration sharing one
# LS term; CHIANTI orders levels by energy, so it is always the first n_gt.
#   C I  2p2 3P_0,1,2     C II 2p 2P_1/2,3/2    N I  2p3 4S_3/2
#   N II 2p2 3P_0,1,2     O I  2p4 3P_2,1,0     O II 2p3 4S_3/2
#   Mg I 3s2 1S0          Mg II 3s 2S_1/2       Ca II 4s 2S_1/2   Na I 3s 2S_1/2
# Mg I, Mg II, Ca II and Na I are NOT here: their coefficients are single
# permitted resonance lines whose critical densities are 3.5e13 to 1.1e16
# cm^-3, decades above anything this code reaches, so their suppression is 1
# to within the printing precision and they keep the coronal fit. Fe I has no
# CHIANTI model atom in v11 and its upper levels are permitted; Fe II already
# carries its own two-dimensional statistical-equilibrium table.
IONS = [("c_1", "CI", 3), ("c_2", "CII", 2), ("n_1", "NI", 1),
        ("n_2", "NII", 3), ("o_1", "OI", 3), ("o_2", "OII", 1)]


def statistical_equilibrium_cooling(chianti_name, T, ne, n_ground_term=0):
    """Radiative loss of one ion divided by n_e [erg cm^3 s^-1].

    T and ne are equal-length arrays of the (T, n_e) pairs to evaluate.  Only
    transitions whose UPPER level lies above n_ground_term are summed, so
    n_ground_term = 0 gives the whole ion and the number of ground-term levels
    gives the part Cool_coeff.f90 takes from the table.  This is the one place
    the sum is formed: the table builder below and any verification of it call
    the same function.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ne = np.atleast_1d(np.asarray(ne, dtype=float))
    ion = ch.ion(chianti_name, temperature=T, eDensity=ne)
    ion.populate()
    f = np.asarray(ion.Population["population"])
    upper = np.asarray(ion.Wgfa["lvl2"])
    avalue = np.asarray(ion.Wgfa["avalue"])
    wvl = np.asarray(ion.Wgfa["wvl"])
    keep = (wvl != 0.0) & (avalue > 0.0) & (upper > n_ground_term)
    e_ul = H_TIMES_C / (np.abs(wvl[keep]) * 1.0e-8)
    return (f[:, upper[keep] - 1] * (avalue[keep] * e_ul)[None, :]).sum(axis=1)/ne


def lambda_se_remainder(chianti_name, n_ground_term):
    """The same sum on the (log T, log n_e) grid of the table, shape (41, 29)."""
    T = np.repeat(10.0**LOGT, len(LOGNE))
    ne = np.tile(10.0**LOGNE, len(LOGT))
    lam = statistical_equilibrium_cooling(chianti_name, T, ne, n_ground_term)
    return lam.reshape(len(LOGT), len(LOGNE))


def fortran_array(name, values, shape2d=False):
    """Emit one Fortran parameter array, eight numbers to a line."""
    flat = np.asarray(values).reshape(-1, order="F") if shape2d else np.asarray(values)
    out = []
    dim = ("(NCOOLT,NCOOLNE)" if shape2d else "(NCOOLT)")
    out.append(f"   real*8, parameter :: {name}{dim} = reshape( [ &"
               if shape2d else f"   real*8, parameter :: {name}{dim} = [ &")
    for i in range(0, flat.size, 6):
        chunk = flat[i:i + 6]
        sep = "," if i + 6 < flat.size else ""
        out.append("        " + ", ".join(f"{v:14.6f}d0" for v in chunk) + sep + " &")
    out.append("        ], [NCOOLT,NCOOLNE] )" if shape2d else "        ]")
    return "\n".join(out)


def main():
    import os
    out = (sys.argv[1] if len(sys.argv) > 1 else
           os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "metal_cooling_density_resolved_tables.txt"))
    with open(out, "w") as fh:
        fh.write("   ! GENERATED by "
                 "cooling_data/metal_cooling_density_resolved.py\n")
        fh.write("   ! from CHIANTI v11.0.2. Do not edit by hand.\n")
        for name, tag, n_gt in IONS:
            lam = lambda_se_remainder(name, n_gt)
            floor = lam[:, 0].copy()
            sys.stderr.write(f"{tag:5s} Lambda_rem(1e4 K, floor) "
                             f"{floor[np.argmin(abs(LOGT-4.0))]:.4e}  "
                             f"suppression at 1e9 cm^-3, 1e4 K "
                             f"{lam[np.argmin(abs(LOGT-4.0)), -10]/floor[np.argmin(abs(LOGT-4.0))]:.3e}\n")
            fh.write(f"\n   ! {tag}\n")
            fh.write(fortran_array(f"cool_logLrem_{tag}",
                                   np.log10(np.maximum(lam, 1.0e-300)),
                                   shape2d=True) + "\n")
    sys.stderr.write(f"written {out}\n")


if __name__ == "__main__":
    main()
