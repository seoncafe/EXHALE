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

with f_u the fractional level populations of the full multilevel statistical
equilibrium of the ion under electron collisions and spontaneous decay.  The
ground-term transitions are excluded because Cool_coeff.f90 solves that term exactly itself, with the neutral-hydrogen
de-excitation channel and the line escape probability that CHIANTI has no way
to carry; leaving them in the table would double count them.

The floor column of the table is the coronal limit of those same channels, so
it is directly comparable with the closed-form fit it replaces.  The two are
NOT equal at low temperature: the fits were built from the
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

ATOMIC DATA.  CHIANTI v11.0.2 (.elvlc, .wgfa, .scups), set XUVTOP.  The
populations are the solution of multilevel_statistical_equilibrium.py, and
the effective collision strengths are descaled by chianti_cooling.upsilon,
the one descaling routine of this directory (the natural cubic spline of
CHIANTI's DESCALE_SCUPS.PRO).  ChiantiPy 0.15.2 `populate` solves the same
equilibrium but descales with a not-a-knot spline, and for C II and N II it
adds proton-impact excitation from the .psplups files at a proton density
set by CHIANTI's coronal ionization equilibrium, which is not the proton
density of the gas this table is used for; neither is taken here.
multilevel_statistical_equilibrium.validate() compares the two solvers with
the descaling made equal.

GRID.  The axes are the ones Cool_coeff.f90 already declares for the Fe II
table, so no new axis is introduced: log10 T from 3.0 to 5.0 in 0.05 dex (41
points) and log10 n_e from 0 to 14 in 0.5 dex (29 points).  The floor
n_e = 1 cm^-3 is four decades below the lowest critical density of any
transition in the sum (the metastable terms of C I, N I and O II, ~1e4 cm^-3),
the fine-structure lines having been excluded.  The Mg I 3s3p 3P0 level has
no radiative decay in CHIANTI at all; its population is independent of n_e
at low density, so the Mg I column is flat to 1 percent below 1e3 cm^-3 and
the floor is its low-density limit, which lies below the coronal sum of the
excitation rates (0.92 of it at 1e4 K) because an excitation to 3P0 is not
radiated.

RUN:  XUVTOP=<dbase> python3 metal_cooling_density_resolved.py [outfile]
      The default outfile is metal_cooling_density_resolved_tables.txt next to
      this script; its contents are pasted into Cool_coeff.f90.
"""

import sys
import numpy as np
import multilevel_statistical_equilibrium as mse

# log10 T and log10 n_e axes of Cool_coeff.f90 (cool_logT, cool_logne)
LOGT = np.linspace(3.0, 5.0, 41)
LOGNE = np.linspace(0.0, 14.0, 29)

# CHIANTI name, the Fortran tag, and the number of levels of the ground term.
# The ground term is the lowest levels of the ground configuration sharing one
# LS term; CHIANTI orders levels by energy, so it is always the first n_gt.
#   C I  2p2 3P_0,1,2     C II 2p 2P_1/2,3/2    N I  2p3 4S_3/2
#   N II 2p2 3P_0,1,2     O I  2p4 3P_2,1,0     O II 2p3 4S_3/2
#   Mg I 3s2 1S0          Mg II 3s 2S_1/2       Ca II 4s 2S_1/2   Na I 3s 2S_1/2
# Ca II and Mg I are here with their single ground level as the "ground
# term", so their tables carry every channel of the ion: the Ca II 4s-3d
# excitation at 1.7 eV (the metastable 3d 2D levels, the [Ca II] 7291/7324
# lines and the 3d-4p infrared triplet they feed) and the Mg I 3s3p 3P
# excitation at 2.7 eV (the 4571 A intercombination line, the metastable
# 3P0 and 3P2 levels and the 3s4s 3S 5167-5184 A triplet above them), besides
# the resonance lines.  At 5e3-1e4 K the critical densities are 4e5-1e7
# cm^-3 for Ca II 3d 2D, 6e4-9e4 cm^-3 for Mg I 3P2 and 4e9-5e9 cm^-3 for
# Mg I 3P1 (3P0 has no radiative decay), inside the range of a wind base.  CHIANTI v11.0.2 has no
# electron collision strengths among the Mg I 3P_J fine-structure levels, and
# the table inherits that omission.
# Mg II and Na I are NOT here: every excited level of Mg II decays by a
# permitted line, and the full Mg II statistical equilibrium equals its
# coronal limit to 0.2 percent at n_e <= 1e12 cm^-3 (2 percent at 1e13), so
# Mg II keeps a closed-form coronal fit (magnesium_ii_line_cooling.py).  Na I
# has no CHIANTI model atom in v11.0.2.  Fe I has no CHIANTI model atom
# either; Fe II carries its own two-dimensional statistical-equilibrium table.
IONS = [("c_1", "CI", 3), ("c_2", "CII", 2), ("n_1", "NI", 1),
        ("n_2", "NII", 3), ("o_1", "OI", 3), ("o_2", "OII", 1),
        ("ca_2", "CaII", 1), ("mg_1", "MgI", 1)]


def statistical_equilibrium_cooling(chianti_name, T, ne, n_ground_term=0):
    """Radiative loss of one ion divided by n_e [erg cm^3 s^-1].

    T and ne are the two axes; returns shape (T.size, ne.size).  Only
    transitions whose UPPER level lies above the first n_ground_term levels
    (energy order) are summed, so n_ground_term = 0 gives the whole ion and
    the number of ground-term levels gives the part Cool_coeff.f90 takes
    from the table.  This is the one place the sum is formed: the table
    builder below and any verification of it call the same function.
    """
    elem, stage = chianti_name.split("_")
    model = mse.load_model_ion(elem, int(stage))
    return mse.line_cooling(model, T, ne, n_ground_term)


def lambda_se_remainder(chianti_name, n_ground_term):
    """The same sum on the (log T, log n_e) grid of the table, shape (41, 29)."""
    return statistical_equilibrium_cooling(chianti_name, 10.0**LOGT,
                                           10.0**LOGNE, n_ground_term)


def fortran_array(name, values, shape2d=False):
    """Emit one Fortran parameter array, six numbers to a line."""
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
