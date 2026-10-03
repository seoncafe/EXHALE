#!/usr/bin/env python3
"""Collisional validity of an EXHALE solution: species-resolved Knudsen
numbers, coupling times, exobase and critical point.

Why this exists.  A hydrodynamic wind solution is a continuum solution.  It
is only a solution of the physical problem where the gas is collisional on
the scale over which the flow varies, and above all through its CRITICAL
POINT, which is where the solution's topology -- and therefore its mass flux
-- is fixed.  If the critical point sits above the exobase, the continuum
equations were integrated through a region in which they do not hold: the
result is then reported as UNVALIDATED.  It is not reported as an
overestimate.  Which way a collisionless critical region moves the answer is
not decided by the continuum solution that assumed it away; it takes a
kinetic or transitional-flow calculation.  The static-Maxwellian (Jeans)
expression of section 3 is a formal comparison on the same profile, not a
bound on either.

--------------------------------------------------------------------------
1. Collision model
--------------------------------------------------------------------------

The momentum-transfer collisions are the SAME four limits the element
diffusion operator uses, in the same Chapman-Enskog first approximation, so
that a mean free path quoted here and a diffusion coefficient quoted by the
wind cannot drift apart:

    neutral-neutral    hard-sphere (Banks & Kockarts 1973)      D ~ T^(1/2)/n
    ion-neutral        induced-dipole (Langevin) polarization
                       and the rigid core, 1/D = 1/D_pol + 1/D_hs
                       (an interpolation between the two limits,
                       not a collision integral of the combined
                       potential)                               (both)
    ion-ion, ion-e     screened Coulomb, Spitzer ln(Lambda)     D ~ T^(5/2)/n

Every formula and every constant below is transcribed from
`src/modules/functions/binary_element_diffusion.f90`, with the Fortran
function or parameter named at each definition (names rather than line
numbers, which move as that file is edited).  Python cannot call the Fortran, so this is a second
implementation of the same expressions and nothing else: no constant is
re-chosen here, and the two are meant to be read side by side.

From the binary diffusion coefficient to a collision frequency.  The
friction force density between two species is

    F_s = (n_s n_t k T)/(n D_st) (v_t - v_s) = n_s mu_st nu_st (v_t - v_s),

which is the definition of the binary diffusion coefficient, so

    nu_st = n_t k T / (n mu_st D_st) = n_t k T / (mu_st Dhat_st),
    Dhat_st = n D_st,

and the total carrier density n cancels: only the density of the PARTNER
enters.  The routines below are therefore called with ntot = 1, which
returns Dhat directly and leaves them identical to the Fortran otherwise.

Species.  H I, H II, He I, He II, He III, e, plus the molecular carriers
H2, H2+, H3+ and HeH+ when a run tracks them. The diagnostic adds HeH+ and
electrons to the Fortran H/He friction lists; HeH+ is not a friction carrier
in the production element diffusion operator. The masses use the same
species table: an H2 is one
partner of mass 2 m_H, a helium atom or ion m_He/m_H = 3.9715259 as in the
species table).  He I 2^3S is NOT a separate species -- it is an excited
level inside the He I column (`species_table.f90` bsp_is_excited_level),
and counting it again would double count those atoms.  METALS (and the
oxygen-chemistry molecules OH, H2O, CO) ARE EXCLUDED as collision partners.
They are kept in the electron budget, where charge neutrality has to close
exactly, and their largest mass and particle fractions on the profile read
are computed and reported, so the size of the omission is measured on each
run rather than assumed.

Omitted channels, all in the direction of a LONGER mean free path, i.e. a
LARGER Knudsen number, i.e. a more conservative verdict.  Omitting a
collision channel can only lower a collision frequency (the frictions of
separate partners add), so every mean free path below is an overestimate
of the one with the channel included:

This monotonic statement holds with the retained coefficients fixed. It
does not establish an upper bound on the physical mean free path: the
polarization-plus-core interpolation and the estimates for resonant pairs
still require independent validation.

  * the METALS and the oxygen-chemistry molecules, as partners (above).
  * RESONANT charge exchange (H+ + H, He+ + He) is not a channel.  The
    Fortran excludes it because a binary ELEMENT diffusion coefficient is
    driven by friction between the elements and both partners of a resonant
    pair carry the same element.  Here the like-element ion-neutral pair is
    a real momentum-transfer channel and it is computed, but with the
    non-resonant (polarization + core) coefficient, whose friction is the
    smaller one: charge transfer at thermal energies makes the momentum
    transfer of a resonant pair exceed the non-resonant estimate.  The mean
    free paths of H I and H II in the partially ionized layer are therefore
    upper limits.  (The size of the resonant momentum-transfer cross section
    was not checked against a published table in this repository.)
  * ELECTRON-NEUTRAL momentum transfer is not in the collision model (the
    Fortran never needs it; an electron carries neither element).  The
    electron Knudsen number is therefore meaningful only where the gas is
    ionized, and is reported with that caveat rather than silently.  The
    electron is not the momentum carrier of the flow in any case: the bulk
    scales below are built from the heavy particles.

--------------------------------------------------------------------------
2. Definitions (all stated, none implicit)
--------------------------------------------------------------------------

Mean free path of species s, over its partners t (the harmonic sum, i.e.
the frictions add, which is the same rule as Blanc's law in the operator):

    nu_s = sum_t nu_st,   lambda_s = vbar_s / nu_s,   1/lambda_s = sum_t 1/lambda_st

with vbar_s = sqrt(8 k T / (pi m_s)) the mean thermal speed.  Note nu_st is
the Chapman-Enskog momentum-transfer frequency, which for rigid spheres of
diameter d is (4/3) n_t pi d^2 gbar_st, gbar_st = sqrt(8 k T/(pi mu_st))
the mean relative speed.  In a gas of like particles lambda_s is therefore
4/3 shorter than Maxwell's 1/(sqrt(2) n pi d^2) and 1.9x shorter than the
elementary 1/(n pi d^2).  The verdict threshold below (0.1) carries that
factor with room to spare, and the convention is stated so a number quoted
from here is not compared against a differently defined one.

Local structure scale, the length over which the state the continuum
closure expands about actually varies.  That state is the local Maxwellian,
and a local Maxwellian is fixed by (n, T, u) -- density, temperature, bulk
velocity -- so those three logarithmic gradients, and the spherical
divergence the geometry adds, are the complete set:

    1/L^2 = (d ln rho/dr)^2 + (d ln T/dr)^2 + (|dv/dr|/c_s)^2 + (1/r)^2,
    L     = max( that, dr_cell ).

Why the sum of squares and not a minimum over the individual lengths.  The
Chapman-Enskog expansion is an expansion in the change of the state over one
mean free path; the state is a vector, its change over a distance d is the
vector (d/H_rho, d/H_T, ...), and the length of that vector is the RMS
above.  Two consequences, both wanted: L is never longer than the shortest
individual scale (each term enters 1/L^2 with a positive sign), and a field
that happens to go logarithmically flat contributes zero rather than
removing itself from a minimum.  A minimum rule has the opposite behavior
-- if the field it is currently taken from flattens, L jumps to the next
one, discontinuously and by whatever factor separates them.  That is not
hypothetical: with the pressure scale alone, across the heating peak of the
LHS 1140 b run `LHS1140b/archive_20260830/exhale/heh0p55` (section 2 of
`md/collisional_validity.md`) d ln p/dr passes through a broad
near-zero (rho falls and T rises with nearly the same log slope, so the
front is close to isobaric while the gas state changes fast), H_p ran up to
3.05e8 cm at 1.0486 R_p, 15x longer than the H_T = 1.80e7 and
H_rho = 1.68e7 cm of the same layer, and Kn showed a spurious 5x dip at
1.03-1.07 R_p.  Pressure is not in the list above for the same
reason it was the wrong list on its own: p = n k T is a derived field, and
adding it would count the density and temperature gradients a second time
wherever they do not cancel.

Why the velocity enters as |dv/dr|/c_s and not as |(1/v) dv/dr|.  The
first-order term the continuum closure drops is the viscous stress, whose
size relative to the pressure is ~ lambda |dv/dr| / vbar: a change of the
bulk velocity distorts the distribution function in proportion to the
THERMAL speed, not to the local bulk speed.  Normalizing by v instead
diverges at every stagnation point -- the cell-centered velocity of the
base cells crosses v = 0 repeatedly, through the collocated two-cell
odd-even mode those cells carry (a mode of the cell-centered field, not a
wave of the solution: the Riemann face mass flux there is the wind's own)
-- which is what an ad-hoc Mach floor was patching before.  c_s is used in place of vbar so that L stays one length,
common to all species; the two differ by an O(1) factor.

The floor at the local cell width dr_cell is a statement about what a
discrete solution can carry: no structure exists below one cell, so a
gradient claiming one is measuring the mesh.  It binds only in the base
cells that carry the odd-even mode (3 of 500 cells in one of the four runs of
`md/collisional_validity.md`, all at r < 1.002 R_p, four orders of
magnitude below the verdict threshold) and never in a critical region.

    Kn_s = lambda_s / L.

Bulk: lambda_bulk is the mass-weighted mean of lambda_s over the HEAVY
species (electrons excluded, see above), and Kn_bulk = lambda_bulk / L.

    EXOBASE  r_exo : the radius where Kn_bulk = 1, taken at the first
                     crossing from below.  If Kn_bulk < 1 everywhere in the
                     domain the exobase is above the outer boundary and is
                     reported as such, not extrapolated.

    CRITICAL POINT r_s : the sonic point, v = c_s with c_s = sqrt(gamma p /
                     rho) and gamma = 5/3.  p and rho are the cgs columns
                     of the output file, so the mean molecular weight is
                     whatever the solution carries.  In atomic gas this is
                     the code's own sound speed (`eval_dt.f90`,
                     cs = sqrt(gamma_ad*p/rho), gamma_ad = 5/3 in
                     parameters.f90).  Where H2 is present and the caloric
                     EOS is active the code uses gamma_eff(T, composition)
                     of `caloric_eos.f90` instead, which is smaller than
                     5/3; this tool does not, so in H2-bearing cells its
                     c_s is the monatomic value, too large by
                     sqrt(5/3 / gamma_eff).  A critical point inside the
                     molecular layer would be misplaced by that factor in
                     c_s; the size of the effect on the verdicts was not
                     measured.

    CRITICAL REGION : from the radius of peak volumetric heating to the
                     sonic point -- the region that heats and accelerates
                     the wind and sets its topology.  "The critical region
                     is collisional" means max(Kn) < 0.1 across it, taken
                     over the bulk and over every species that carries at
                     least `TRACE_FRACTION` of the mass.  When the domain
                     holds no critical point the solution is a subsonic
                     breeze, its mass flux is set by the outer boundary
                     condition, and the region runs from peak heating to
                     the outer boundary.

Coupling times: tau_s = 1/nu_s against the flow time r/|v|; and the
electron-ion ENERGY coupling time, which is longer than the momentum one by
the fraction of energy a collision transfers,

    nu^E_st = 2 mu_st/(m_s + m_t) nu_st     (elastic, Maxwellian average),

against the heating time tau_heat = (3/2) n_tot k T / Q with Q the heating
rate column.  Where tau^E_ei > tau_heat the single-temperature assumption of
the energy equation is itself unvalidated, independently of the Knudsen
number.

--------------------------------------------------------------------------
3. Formal static-Maxwellian (Jeans) comparison
--------------------------------------------------------------------------

Where the collisional criterion is not met, the Jeans (1925) escape flux of
a static Maxwellian is evaluated at the diagnosed exobase, or at the outer
boundary (flagged) when the exobase is outside the grid:

    Phi_J,s = n_s vbar_s / 4 (1 + lambda_J,s) exp(-lambda_J,s),
    vbar_s  = sqrt(8 k T / (pi m_s))          (mean speed)
    lambda_J,s = G M_p m_s / (k T r_exo)      (the Jeans escape parameter)

and Mdot_Jeans = 4 pi r_exo^2 sum_s m_s Phi_J,s.  With the most probable
speed v_mp = sqrt(2 k T / m_s) the same prefactor reads n_s v_mp/(2 sqrt(pi)).
It is the outward flux of an isotropic Maxwellian through a surface,
integrated over v_r > 0 and speeds above the escape speed
v_esc = sqrt(2 G M_p / r):

    Phi = n (m/(2 pi k T))^(3/2) 2 pi int_{v_esc}^inf v^3 exp(-m v^2/(2kT)) dv
          int_0^1 mu dmu
        = n vbar/4 (1 + lambda_J) exp(-lambda_J).

The expression assumes zero bulk drift at the evaluation radius, a
Maxwellian, a point-mass potential and independent escape of each species.
The diagnostic applies it formally to all heavy species, ions included,
and so cannot account for the ambipolar potential or for the wind's own
velocity at that radius.  Its ratio to the continuum mass flux compares two
formulas on the same profile; it is not a bound on, nor an error estimate
for, a kinetic solution.  At small lambda_J (<~ 2-3) the escaping part is
no longer a small tail of the distribution and the usual interpretation of
the Jeans rate does not apply; the number is still printed and labeled
formal.

--------------------------------------------------------------------------
Which state this reads
--------------------------------------------------------------------------

The default is the SOLUTION, `output/Hydro_ioniz.txt` and
`output/Ion_species.txt`, and the reason is that every quantity here that
decides the verdict is built on the momentum equation: the sound speed, the
Mach number, the critical point, the structure scale L that the Knudsen
number is divided by, and the mass flux.  Those are properties of the state
that satisfies that equation, which is the solution and only the solution.

The `_adv` pair is a different state.  It carries the solution's own density
and velocity beside a temperature and a composition obtained from a second
closure -- the steady ionization balance with advection, in place of the
local photoionization equilibrium the wind was solved with -- for which the
momentum equation was never re-solved.  A sound speed formed from its
pressure and the solution's density, and a Mach number formed from that and
the solution's velocity, therefore belong to no single state.  Measured on
the LHS 1140 b wind with a 45 R_p outer boundary, `LHS1140b/models/.L8/r45`:
the solution crosses its critical point at 40.06 R_p; the `_adv` pair,
whose temperature is 1.5, 2.5 and 3.7 times the solution's at 10, 20 and
30 R_p, puts a "critical point" at 32.5 R_p that no state has.  When
ionization-stage transport
is on in the wind solve (`ionization_transport T` in EXHALE_resolved.out),
the solution itself carries transported ion fractions and its composition
is not the local photoionization equilibrium.

`--adv` reads the `_adv` pair anyway, which is the right question to ask of
the COMPOSITION: the mean free path is set by how much of the gas is
neutral, and in a wind whose ionization cannot relax over a flow time the
advected composition is the physical one while the equilibrium composition
the solution carries is not.  The two answers bracket the exobase rather
than agreeing: on that same run the solution puts Kn_bulk = 1 above the
outer boundary (Kn_bulk = 0.11 at 43.4 R_p) and the `_adv` composition puts
it at 29.0 R_p.  The report
names which state it read; neither is a self-consistent state above the
radius where the ionization stops relaxing, and that is a limitation of the
run, not of this diagnostic.

--------------------------------------------------------------------------
Usage
--------------------------------------------------------------------------

    python3 collisional_validity.py <case_dir> [<case_dir> ...]
        [--adv]       read the advection-corrected _adv profiles instead of
                      the solution (see "Which state this reads", below)
        [--kn-threshold 0.1]
        [--json FILE] write the full diagnosis of every case as JSON
        [--profile FILE] write the radial profile of the FIRST case
                      (r, T, v, c_s, L, lambda_s, Kn_s) as a text table

A case directory is a run directory: it holds `input.inp` and `output/`.
The planet radius and mass are those the run was solved with, which
`exhale_io.load_run` takes from the run's `EXHALE_resolved.out` (the
lower-atmosphere profile handoff, for one, moves the base radius away from
`input.inp`'s `Planet radius`); the report prints the values used and the
file they came from.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'examples'))
import exhale_io as eio     # noqa: E402  (path set above)

# --------------------------------------------------------------------------
# Constants.  Every one of these is the value the Fortran uses, cited.
# --------------------------------------------------------------------------

kb_erg = 1.380649e-16       # CODATA 2018, parameters.f90 kb_erg
m_H_g = eio.mu              # hydrogen ATOM mass [g], parameters.f90 mu
m_e_g = 9.1093837015e-28    # electron mass [g], CODATA 2018, parameters.f90 m_e
G_cgs = 6.67430e-8          # CODATA 2018, parameters.f90 Gc
# parameters.f90:879 m_He_atom (helium-4 atom [g]) and :884 amu (CODATA 2018 [g]).
# The Fortran's m_He_over_m_H = m_He_atom/mu is the mass of every helium
# stage in the species table (bsp_mass) and in the diffusion operator
# (hecar_m), and amu_over_m_H = amu/mu turns an atomic weight in u into the
# code's hydrogen-atom mass unit.
m_He_over_m_H = 6.6464790722e-24 / m_H_g
amu_over_m_H = 1.66053906660e-24 / m_H_g

# binary_element_diffusion.f90 e_esu -- CODATA 2018 exact coulomb value.
e_esu = 4.803204713e-10
# binary_element_diffusion.f90 bk_hs_pref -- Banks & Kockarts (1973)
# hard-sphere prefactor.  Put through the Chapman-Enskog first
# approximation, D = 3/(8 n d^2) (k T/(2 pi mu))^(1/2), it is a rigid-sphere
# collision diameter d = 2.99 Angstrom (a diameter of 2.7 Angstrom would
# give a prefactor of 1.86e18).
bk_hs_pref = 1.52e18
# binary_element_diffusion.f90 alpha_HI / alpha_HeI / alpha_H2 -- static
# dipole polarizabilities of the neutral partners, in a_0^3 converted with
# a_0^3 = 1.481847e-25 cm^3.
# alpha(H) = 4.5 a_0^3 is the exact nonrelativistic ground-state value; He and
# H2 are Schwerdtfeger & Nagle (2019, Mol. Phys. 117, 1200).
a0cub_cm3 = 1.481847e-25
alpha_HI = 4.500 * a0cub_cm3
alpha_HeI = 1.383 * a0cub_cm3
alpha_H2 = 5.315 * a0cub_cm3

# Species: name -> (mass [m_H units, counted for each collision partner],
# charge, neutral polarizability [cm^3] or 0 for a charged partner).  Masses
# and charges are the carrier lists of binary_element_diffusion.f90
# (hcar_* for hydrogen, hecar_* for the helium stages, whose mass is
# m_He_over_m_H) and HeH+ = 1 + m_He_over_m_H as in species_table.f90
# bsp_mass; the electron is added because the Coulomb channel needs it and
# the Fortran, which transports elements, does not.  He I 2^3S is
# deliberately absent: it is inside the He I column.
COLLIDERS = {
    'HI':    (1.0, 0.0, alpha_HI),
    'HII':   (1.0, 1.0, 0.0),
    'HeI':   (m_He_over_m_H, 0.0, alpha_HeI),
    'HeII':  (m_He_over_m_H, 1.0, 0.0),
    'HeIII': (m_He_over_m_H, 2.0, 0.0),
    'H2':    (2.0, 0.0, alpha_H2),
    'H2p':   (2.0, 1.0, 0.0),
    'H3p':   (3.0, 1.0, 0.0),
    'HeHp':  (1.0 + m_He_over_m_H, 1.0, 0.0),
    'e':     (m_e_g / m_H_g, -1.0, 0.0),
}
HEAVY = [s for s in COLLIDERS if s != 'e']

# Charge of every ion column that can appear, for the electron budget.  The
# metals are trace for the collisions but not for charge neutrality.
METAL_CHARGE = {}
for _m in eio.METAL_IONS:
    METAL_CHARGE[_m] = {'I': 0, 'II': 1, 'III': 2}[
        _m[len(_m.rstrip('I')):] or 'I']

gamma_ad = 5.0 / 3.0        # parameters.f90 gamma_ad, the monatomic value

# A species below this mass fraction is not asked to validate the solution:
# its Knudsen number is reported but does not enter the verdict.
TRACE_FRACTION = 1.0e-3

# Standard atomic weights [u] of the metal elements, a copy of
# species_table.f90:97 melem_A_u (same element order C O N Mg Si Ca Na K S Fe),
# converted below to the code's hydrogen-atom unit exactly as the Fortran
# forms melem_A = amu_over_m_H*melem_A_u.  Keep the two in step.
METAL_ATOMIC_WEIGHT_U = dict(zip(
    ('C', 'O', 'N', 'Mg', 'Si', 'Ca', 'Na', 'K', 'S', 'Fe'),
    (12.011, 15.999, 14.007, 24.305, 28.085, 40.078, 22.990, 39.098,
     32.06, 55.845)))
# Metal nuclei carried by the oxygen-chemistry molecules (species_table.f90
# bsp_nO, bsp_nC); their hydrogen is counted with the hydrogen.
MOLECULE_METAL_NUCLEI = {'OH': {'O': 1}, 'H2O': {'O': 1},
                         'CO': {'C': 1, 'O': 1}}


def omitted_partner_fractions(ion, rho_over_mH, heavy_particles):
    """Mass and particle fractions of the species this collision model
    leaves out: every metal atom and ion, and the oxygen-chemistry molecules
    OH, H2O and CO.  Mass fraction = metal-nucleus mass / rho (the H in
    those molecules is hydrogen mass); particle fraction = omitted particles
    / (omitted + heavy partners).  Returns two radial arrays."""
    mass = np.zeros_like(rho_over_mH)
    particles = np.zeros_like(rho_over_mH)
    for element, weight_u in METAL_ATOMIC_WEIGHT_U.items():
        stages = [element + 'I', element + 'II', element + 'III']
        for s in stages:
            if s in ion:
                mass = mass + amu_over_m_H * weight_u * ion[s]
                particles = particles + ion[s]
    for mol, nuclei in MOLECULE_METAL_NUCLEI.items():
        if mol in ion:
            for element, count in nuclei.items():
                mass = mass + (amu_over_m_H * METAL_ATOMIC_WEIGHT_U[element]
                               * count * ion[mol])
            particles = particles + ion[mol]
    mass_frac = mass / np.maximum(rho_over_mH, 1.0e-300)
    part_frac = particles / np.maximum(particles + heavy_particles, 1.0e-300)
    return mass_frac, part_frac


# --------------------------------------------------------------------------
# Pair diffusion coefficients -- transcribed from
# src/modules/functions/binary_element_diffusion.f90
# --------------------------------------------------------------------------

def hard_sphere_pair_diffusion(TK, ntot, A_s, A_t):
    """NEUTRAL-NEUTRAL.  Banks & Kockarts (1973):
    D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2) / n  [cm^2/s].
    Fortran: binary_element_diffusion.f90 hard_sphere_pair_diffusion."""
    return bk_hs_pref * np.sqrt(1.0 / A_s + 1.0 / A_t) * np.sqrt(TK) / ntot


def polarization_pair_diffusion(TK, ntot, alpha_n, mu_g):
    """ION-NEUTRAL, induced-dipole (Langevin) channel:
    D = k T / (2.21 pi e n (alpha_n mu)^(1/2))  [cm^2/s], the low-energy
    limit -- never used alone, see ion_neutral_pair_diffusion.
    Fortran: binary_element_diffusion.f90 polarization_pair_diffusion."""
    return kb_erg * TK / (2.21 * np.pi * e_esu * ntot
                          * np.sqrt(np.maximum(alpha_n * mu_g, 1.0e-60)))


def coulomb_logarithm(TK, ne, zz):
    """ln(Lambda) = ln(3 k T lambda_D / (Z_s Z_t e^2)) with the electron
    Debye length (Spitzer 1962, section 5.2), floored at 1.
    Fortran: binary_element_diffusion.f90 coulomb_logarithm."""
    lam_D = np.sqrt(kb_erg * TK / (4.0 * np.pi * np.maximum(ne, 1.0)
                                   * e_esu * e_esu))
    return np.maximum(np.log(3.0 * kb_erg * TK * lam_D
                             / (max(zz, 1.0) * e_esu * e_esu)), 1.0)


def coulomb_pair_diffusion(TK, ntot, ne, Z_s, Z_t, mu_g):
    """ION-ION (and ion-electron).  Screened Coulomb through the
    Chapman-Enskog integral:
    D = 3 (kT)^(5/2) / [4 (2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln(Lambda)].
    Fortran: binary_element_diffusion.f90 coulomb_pair_diffusion."""
    zz = abs(Z_s * Z_t)
    return (3.0 * (kb_erg * TK) ** 2.5
            / (4.0 * np.sqrt(2.0 * np.pi * mu_g) * ntot
               * (zz * e_esu * e_esu) ** 2 * coulomb_logarithm(TK, ne, zz)))


def ion_neutral_pair_diffusion(TK, ntot, alpha_n, A_s, A_t, mu_g):
    """ION-NEUTRAL, non-resonant: polarization and rigid core combined as
    1/D = 1/D_pol + 1/D_hs.  This is an interpolation between the two
    limits (it recovers each where the other is negligible), not a
    collision integral of the combined polarization-plus-core potential;
    in the crossover it can overstate the friction by up to a factor of 2
    against the stronger estimate alone.
    Fortran: binary_element_diffusion.f90 ion_neutral_pair_diffusion."""
    Dpol = polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
    Dhs = hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
    return 1.0 / (1.0 / np.maximum(Dpol, 1.0e-99)
                  + 1.0 / np.maximum(Dhs, 1.0e-99))


def pair_diffusion(TK, ntot, ne, s, t):
    """The pair coefficient in whichever of the three limits the pair
    belongs to -- the Python twin of `stage_pair_diffusion`
    (binary_element_diffusion.f90)."""
    m_s, Z_s, al_s = COLLIDERS[s]
    m_t, Z_t, al_t = COLLIDERS[t]
    mu_g = m_s * m_t / (m_s + m_t) * m_H_g
    if Z_s == 0.0 and Z_t == 0.0:
        return hard_sphere_pair_diffusion(TK, ntot, m_s, m_t)
    if Z_s != 0.0 and Z_t != 0.0:
        return coulomb_pair_diffusion(TK, ntot, ne, Z_s, Z_t, mu_g)
    if Z_s == 0.0:
        return ion_neutral_pair_diffusion(TK, ntot, al_s, m_s, m_t, mu_g)
    return ion_neutral_pair_diffusion(TK, ntot, al_t, m_s, m_t, mu_g)


# --------------------------------------------------------------------------
# Collision frequencies, mean free paths, Knudsen numbers
# --------------------------------------------------------------------------

def momentum_transfer_frequency(TK, ne, dens, s, t):
    """nu_st = n_t k T / (mu_st Dhat_st) [s^-1], the momentum-transfer
    frequency defined by F_s = n_s mu_st nu_st (v_t - v_s), from the same
    binary diffusion coefficient the wind's diffusion operator uses.  The
    velocity relaxation rate of s against fixed t is (mu_st/m_s) nu_st;
    the mean-free-path diagnostic uses nu_st itself as its convention. The
    total carrier density cancels, so `pair_diffusion` is called with
    ntot = 1 and returns Dhat = n D."""
    m_s = COLLIDERS[s][0] * m_H_g
    m_t = COLLIDERS[t][0] * m_H_g
    mu_g = m_s * m_t / (m_s + m_t)
    Dhat = pair_diffusion(TK, 1.0, ne, s, t)
    return dens[t] * kb_erg * TK / (mu_g * Dhat)


def mean_thermal_speed(TK, s):
    """vbar_s = sqrt(8 k T / (pi m_s)) [cm/s]."""
    return np.sqrt(8.0 * kb_erg * TK / (np.pi * COLLIDERS[s][0] * m_H_g))


def collision_rates(TK, ne, dens):
    """Per species: total momentum-transfer frequency nu_s = sum_t nu_st,
    mean free path lambda_s = vbar_s/nu_s, and the pair breakdown.
    Returns (nu, lam, nu_pair)."""
    present = [s for s in COLLIDERS if s in dens]
    nu, lam, nu_pair = {}, {}, {}
    for s in present:
        tot = np.zeros_like(TK)
        for t in present:
            if s == 'e' and COLLIDERS[t][1] == 0.0:
                continue        # electron-neutral is not in the model
            if t == 'e' and COLLIDERS[s][1] == 0.0:
                continue
            nst = momentum_transfer_frequency(TK, ne, dens, s, t)
            nu_pair[(s, t)] = nst
            tot = tot + nst
        nu[s] = tot
        lam[s] = np.where(tot > 0.0, mean_thermal_speed(TK, s)
                          / np.where(tot > 0.0, tot, 1.0), np.inf)
    return nu, lam, nu_pair


def energy_equipartition_frequency(TK, ne, dens, s, t):
    """nu^E_st = 2 mu_st/(m_s+m_t) nu_st: the elastic energy-transfer rate,
    the momentum rate times the fraction of energy one collision moves."""
    m_s = COLLIDERS[s][0] * m_H_g
    m_t = COLLIDERS[t][0] * m_H_g
    mu_g = m_s * m_t / (m_s + m_t)
    return (2.0 * mu_g / (m_s + m_t)
            * momentum_transfer_frequency(TK, ne, dens, s, t))


def structure_scale(r_cm, rho, T, v, cs):
    """The length over which the local Maxwellian (n, T, u) varies, plus the
    spherical divergence, combined as the length of the state-vector
    gradient and floored at the cell width (module header):

        1/L^2 = (dln rho/dr)^2 + (dln T/dr)^2 + (|dv/dr|/c_s)^2 + (1/r)^2
        L     = max(1/sqrt(...), dr_cell)

    Returns (L, scales), scales holding the individual lengths for
    reporting: H_rho, H_T, L_v = c_s/|dv/dr|, r, dr_cell."""
    def inv_logscale(x):
        return np.abs(np.gradient(np.log(np.maximum(x, 1.0e-300)), r_cm))

    g_rho = inv_logscale(rho)
    g_T = inv_logscale(T)
    g_v = np.abs(np.gradient(v, r_cm)) / cs
    g_r = 1.0 / r_cm
    dr_cell = np.gradient(r_cm)
    L = 1.0 / np.maximum(np.sqrt(g_rho ** 2 + g_T ** 2 + g_v ** 2
                                 + g_r ** 2), 1.0e-99)
    L = np.maximum(L, dr_cell)
    scales = {'H_rho': 1.0 / np.maximum(g_rho, 1.0e-99),
              'H_T': 1.0 / np.maximum(g_T, 1.0e-99),
              'L_v': 1.0 / np.maximum(g_v, 1.0e-99),
              'r': r_cm, 'dr_cell': dr_cell}
    return L, scales


# --------------------------------------------------------------------------
# Radii with stated definitions
# --------------------------------------------------------------------------

def _first_upcrossing(r, f):
    """First radius where f crosses zero from below, linearly interpolated.
    None if it never does."""
    idx = np.nonzero((f[:-1] < 0.0) & (f[1:] >= 0.0))[0]
    if idx.size == 0:
        return None
    j = idx[0]
    w = -f[j] / (f[j + 1] - f[j])
    return r[j] + w * (r[j + 1] - r[j])


def sonic_point(r, v, cs):
    """Critical point: v = c_s, c_s = sqrt(gamma p/rho), gamma = 5/3 -- the
    code's own sound speed in atomic gas (eval_dt.f90); in H2-bearing cells
    with the caloric EOS active the code uses gamma_eff < 5/3 (module
    header, CRITICAL POINT)."""
    return _first_upcrossing(r, v - cs)


def exobase(r, Kn_bulk):
    """Exobase: Kn_bulk = 1, first crossing from below."""
    return _first_upcrossing(r, Kn_bulk - 1.0)


def jeans_escape(r_exo_cm, TK, dens, Mp_g):
    """Formal static-Maxwellian (Jeans 1925) escape at radius r_exo_cm
    (module header, section 3):

        Phi_J,s = n_s vbar_s/4 (1 + lambda_J) exp(-lambda_J),
        vbar_s  = sqrt(8 k T/(pi m_s)),
        lambda_J,s = G M_p m_s/(k T r_exo).

    Returns (Mdot [g/s], {species: lambda_J}, {species: Phi [cm^-2 s^-1]}).
    Heavy species only, each treated as escaping on its own; the ambipolar
    field that holds the electrons and acts on the ions is not modeled."""
    Mdot, lam_J, phi = 0.0, {}, {}
    for s in HEAVY:
        if s not in dens:
            continue
        m_s = COLLIDERS[s][0] * m_H_g
        lj = G_cgs * Mp_g * m_s / (kb_erg * TK * r_exo_cm)
        f = (dens[s] * mean_thermal_speed(TK, s) / 4.0
             * (1.0 + lj) * np.exp(-lj))
        lam_J[s], phi[s] = lj, f
        Mdot += 4.0 * np.pi * r_exo_cm ** 2 * m_s * f
    return Mdot, lam_J, phi


# --------------------------------------------------------------------------
# One case
# --------------------------------------------------------------------------

def collisional_diagnosis(case_dir, adv=False, kn_threshold=0.1,
                          outdir='output'):
    """Everything above, on one EXHALE run directory.  Reads
    output/Hydro_ioniz{_adv}.txt and output/Ion_species{_adv}.txt; the
    solution is the default, because the sound speed, the Mach number, the
    critical point and the structure scale are properties of the state that
    satisfies the momentum equation ("Which state this reads", above)."""
    inp_path = os.path.join(case_dir, 'input.inp')
    run = eio.load_run(os.path.join(case_dir, outdir), inp_path, adv=adv)

    # The ghost rows are already gone: exhale_io returns the physical cells
    # (exhale_io.physical_cell_rows).  Cutting them again here would take
    # four more solution cells off each end.
    # Rp_RJ and Mp_MJ are the radius and mass the run was solved with:
    # load_run takes them from EXHALE_resolved.out (see exhale_io.load_run).
    Rp_cm = run.inp['Rp_RJ'] * eio.RJ
    Mp_g = run.inp['Mp_MJ'] * eio.MJ
    r = run.r
    r_cm = r * Rp_cm
    T = run.T
    v = run.v
    p = run.p
    rho = run.n * m_H_g
    heat = run.heat

    ion = dict(run.ion)
    dens = {s: ion[s] for s in COLLIDERS if s in ion and s != 'e'}
    heavy_particles = sum((dens[s] for s in HEAVY if s in dens),
                          np.zeros_like(T))
    omit_mass, omit_part = omitted_partner_fractions(ion, run.n,
                                                     heavy_particles)
    # Electron density from charge neutrality over EVERY tracked ion, metals
    # included: the collisions treat the metals as trace, the charge budget
    # must not.
    ne = np.zeros_like(T)
    for s, (_m, Z, _a) in COLLIDERS.items():
        if s != 'e' and Z > 0.0 and s in ion:
            ne = ne + Z * ion[s]
    for s, Z in METAL_CHARGE.items():
        if Z > 0 and s in ion:
            ne = ne + Z * ion[s]
    dens['e'] = ne

    cs = np.sqrt(gamma_ad * p / rho)
    mach = v / cs
    L, scales = structure_scale(r_cm, rho, T, v, cs)

    nu, lam, _nu_pair = collision_rates(T, ne, dens)
    Kn = {s: lam[s] / L for s in lam}

    # Bulk: mass-weighted over the heavy particles.
    m_heavy = np.zeros_like(T)
    lam_bulk = np.zeros_like(T)
    for s in HEAVY:
        if s not in dens:
            continue
        w = COLLIDERS[s][0] * m_H_g * dens[s]
        m_heavy += w
        lam_bulk += w * lam[s]
    lam_bulk /= np.maximum(m_heavy, 1.0e-300)
    Kn_bulk = lam_bulk / L
    mass_frac = {s: float(np.max(COLLIDERS[s][0] * m_H_g * dens[s] / m_heavy))
                 for s in HEAVY if s in dens}

    r_s = sonic_point(r, v, cs)
    r_exo = exobase(r, Kn_bulk)
    # Where the gas stops being collisional on the flow's own scale.  This
    # is the radius the verdict actually turns on when the domain holds no
    # critical point: it says how far out the continuum solution is a
    # solution at all.
    r_kn_thr = _first_upcrossing(r, Kn_bulk - kn_threshold)
    j_heat = int(np.argmax(heat))
    r_heat = float(r[j_heat])
    j_Tmax = int(np.argmax(T))

    # Critical region: peak heating -> sonic point (or the outer boundary,
    # flagged, when the domain holds no sonic point).
    r_hi = r_s if r_s is not None else float(r[-1])
    crit = (r >= min(r_heat, r_hi)) & (r <= max(r_heat, r_hi))
    if not crit.any():
        crit = np.ones_like(r, dtype=bool)

    kn_crit = {s: float(np.max(Kn[s][crit])) for s in Kn}
    kn_crit_bulk = float(np.max(Kn_bulk[crit]))
    kn_crit_validating = max(
        [kn_crit_bulk] + [kn_crit[s] for s in HEAVY
                          if s in kn_crit and mass_frac.get(s, 0.0)
                          >= TRACE_FRACTION])

    def at(rr, arr):
        return float(np.interp(rr, r, arr)) if rr is not None else None

    # Coupling and energy times.
    tau_flow = r_cm / np.maximum(np.abs(v), 1.0e-30)
    tau_heat = 1.5 * sum(dens[s] for s in dens) * kb_erg * T \
        / np.maximum(heat, 1.0e-300)
    nuE_ei = np.zeros_like(T)
    for s in ('HII', 'HeII', 'HeIII', 'H2p', 'H3p', 'HeHp'):
        if s in dens:
            nuE_ei += (dens[s] / np.maximum(ne, 1.0e-300)) \
                * energy_equipartition_frequency(T, ne, dens, 'e', s)
    tau_E_ei = 1.0 / np.maximum(nuE_ei, 1.0e-300)

    def flag(key):
        v = run.resolved.get(key)
        return {'T': True, 'F': False}.get(v) if v is not None else None

    res = dict(
        case=case_dir, adv=adv, kn_threshold=kn_threshold,
        planet_source=run.inp.get('planet_source'),
        Rp_RJ=run.inp['Rp_RJ'], Mp_MJ=run.inp['Mp_MJ'],
        Rp_RJ_input=run.inp.get('Rp_RJ_input'),
        carrier_transport=flag('carrier_transport'),
        ionization_transport=flag('ionization_transport'),
        omitted_mass_fraction_max=float(np.max(omit_mass)),
        omitted_particle_fraction_max=float(np.max(omit_part)),
        Rp_cm=Rp_cm, Mp_g=Mp_g,
        r=r, T=T, v=v, cs=cs, mach=mach, p=p, rho=rho, ne=ne,
        L=L, scales=scales, lam=lam, Kn=Kn,
        lam_bulk=lam_bulk, Kn_bulk=Kn_bulk, nu=nu,
        tau_flow=tau_flow, tau_heat=tau_heat, tau_E_ei=tau_E_ei,
        mass_frac=mass_frac,
        r_sonic=r_s, r_exobase=r_exo, r_heat_peak=r_heat,
        r_kn_threshold=r_kn_thr,
        r_Tmax=float(r[j_Tmax]), T_max=float(T[j_Tmax]),
        T_sonic=at(r_s, T), mach_max=float(np.max(mach)),
        Kn_at_sonic=at(r_s, Kn_bulk),
        Kn_at_Tmax=float(np.interp(r[j_Tmax], r, Kn_bulk)),
        Kn_crit_max=kn_crit_validating, Kn_crit_bulk_max=kn_crit_bulk,
        Kn_crit_species=kn_crit,
        Kn_top=float(Kn_bulk[-1]), r_top=float(r[-1]),
        tau_E_over_tau_heat_crit=float(np.max(
            (tau_E_ei / np.maximum(tau_heat, 1.0e-300))[crit])),
        log10_Mdot=float(eio.mdot_log10(run)),
    )

    # The kinetic scale, evaluated at the exobase when there is one inside
    # the domain, otherwise at the outer boundary (flagged).
    if r_exo is not None:
        r_k, where = r_exo, 'exobase'
    else:
        r_k, where = float(r[-1]), 'outer boundary (no exobase in domain)'
    dens_k = {s: float(np.interp(r_k, r, dens[s])) for s in dens}
    T_k = float(np.interp(r_k, r, T))
    Mdot_J, lam_J, _phi = jeans_escape(r_k * Rp_cm, T_k, dens_k, Mp_g)
    res.update(jeans_at=where, r_jeans=r_k, T_jeans=T_k,
               Mdot_jeans=Mdot_J,
               lambda_J={s: float(x) for s, x in lam_J.items()},
               lambda_J_bulk=float(min(lam_J.values())) if lam_J else None)
    res['validity_statement'] = validity_statement(res)
    return res


def validity_statement(res):
    """The verdict sentence, with the numbers it rests on."""
    r_s, r_exo = res['r_sonic'], res['r_exobase']
    kn = res['Kn_crit_max']
    thr = res['kn_threshold']
    lines = []

    if r_s is None:
        lines.append(
            "NO CRITICAL POINT IN THE DOMAIN: the solution stays subsonic "
            "out to %.2f R_p (max Mach %.3g), so it is not a transonic wind "
            "and its mass flux is set by the outer boundary condition, where "
            "Kn_bulk = %.3g." % (res['r_top'], res['mach_max'],
                                 res['Kn_top']))
        collisional_through_crit = False
    else:
        collisional_through_crit = (kn < thr) and (
            r_exo is None or r_exo > r_s)

    if r_exo is None:
        exo_txt = ("exobase above the outer boundary (Kn_bulk = %.3g at "
                   "%.2f R_p)" % (res['Kn_top'], res['r_top']))
    else:
        exo_txt = "exobase at %.3f R_p" % r_exo

    if r_s is not None:
        lines.append(
            "critical (sonic) point at %.3f R_p, %s; max Kn over the "
            "heating/acceleration region %.3f R_p - %.3f R_p is %.3g "
            "(threshold %.2g)."
            % (r_s, exo_txt, res['r_heat_peak'], r_s, kn, thr))
    else:
        lines.append("%s; max Kn from peak heating (%.3f R_p) to the "
                     "outer boundary (%.3f R_p) is %.3g (threshold %.2g)."
                     % (exo_txt[0].upper() + exo_txt[1:], res['r_heat_peak'],
                        res['r_top'], kn, thr))

    if collisional_through_crit:
        lines.append(
            "COLLISIONAL CRITERION MET through the critical point, in this "
            "collision model: Kn < %.2g across the heating and acceleration "
            "region, and the critical point lies below the exobase. "
            "log10 Mdot = %.3f."
            % (thr, res['log10_Mdot']))
    else:
        why = []
        if r_s is not None and r_exo is not None and r_exo <= r_s:
            why.append("the exobase (%.3f R_p) lies BELOW the critical point "
                       "(%.3f R_p)" % (r_exo, r_s))
        if kn >= thr:
            why.append("Kn reaches %.3g in the heating/acceleration region"
                       % kn)
        if r_s is None:
            why.append("the domain holds no critical point")
        lines.append(
            "COLLISIONAL CRITERION NOT MET: " + "; ".join(why) + ". "
            "The hydrodynamic result is unvalidated, not overestimated: the "
            "sign and size of the error are not determined by this solution. "
            "log10 Mdot = %.3f." % res['log10_Mdot'])
        lam_b = res['lambda_J_bulk']
        ratio = 10.0 ** res['log10_Mdot'] / max(res['Mdot_jeans'], 1e-300)
        lines.append(
            "Formal static-Maxwellian (Jeans) rate at the %s (T = %.0f K, "
            "lambda_J = %.2f for the lightest heavy species): "
            "Mdot_Jeans = %.3g g/s; the continuum Mdot is %.2g times it. "
            "That ratio compares two formulas on one profile; it is neither "
            "a bound on nor an error estimate for a kinetic solution."
            % (res['jeans_at'], res['T_jeans'], lam_b, res['Mdot_jeans'],
               ratio))
        if lam_b is not None and lam_b < 3.0:
            lines.append(
                "At lambda_J = %.2f the escaping part is not a small tail of "
                "the Maxwellian, so the usual reading of the Jeans rate does "
                "not apply either." % lam_b)
        else:
            lines.append(
                "The expression assumes zero bulk drift at that radius and "
                "omits the ambipolar potential of escaping ions; the "
                "physical escape rate needs a kinetic calculation.")
    return "  ".join(lines)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

_SCALAR_KEYS = ['r_sonic', 'r_exobase', 'r_heat_peak', 'r_Tmax', 'T_max',
                'T_sonic', 'mach_max', 'Kn_at_sonic', 'Kn_at_Tmax',
                'r_kn_threshold', 'Kn_crit_max',
                'Kn_crit_bulk_max', 'Kn_top', 'r_top', 'log10_Mdot',
                'jeans_at', 'r_jeans', 'T_jeans', 'Mdot_jeans',
                'lambda_J_bulk', 'tau_E_over_tau_heat_crit',
                'planet_source', 'Rp_RJ', 'Mp_MJ', 'Rp_RJ_input',
                'carrier_transport', 'ionization_transport',
                'omitted_mass_fraction_max', 'omitted_particle_fraction_max']


def report(res):
    if res['adv']:
        state = ('advection-corrected _adv profiles: the solution\'s density'
                 ' and velocity with a temperature and composition from a'
                 ' closure the momentum equation was not re-solved for')
    elif res['ionization_transport'] is None:
        state = ('the solution; whether its ion fractions are transported is'
                 ' not recorded (EXHALE_resolved.out absent or without an'
                 ' ionization_transport line)')
    elif res['ionization_transport']:
        state = ('the solution; its H/He ion fractions are transported by'
                 ' the wind, not local photoionization equilibrium')
    else:
        state = ('the solution; its composition is the local photoionization'
                 ' equilibrium, which the mean free path inherits')
    print('=' * 74)
    print('case: %s   (%s)' % (res['case'], state))
    print('=' * 74)
    src = res['planet_source'] or 'unknown'
    print('  planet radius          %.9g R_J   mass %.9g M_J   (from %s)'
          % (res['Rp_RJ'], res['Mp_MJ'], src))
    if res['Rp_RJ_input'] is not None and res['Rp_RJ_input'] != res['Rp_RJ']:
        print('                         (input.inp Planet radius %.9g R_J,'
              ' ratio %.6f)' % (res['Rp_RJ_input'],
                                res['Rp_RJ'] / res['Rp_RJ_input']))
    print('  omitted partners (metals, OH/H2O/CO): max mass fraction %.3g,'
          ' max particle fraction %.3g'
          % (res['omitted_mass_fraction_max'],
             res['omitted_particle_fraction_max']))
    print('  peak heating at        %8.3f R_p' % res['r_heat_peak'])
    print('  T max                  %8.3f R_p   %.0f K'
          % (res['r_Tmax'], res['T_max']))
    if res['r_sonic'] is not None:
        print('  critical (sonic) point %8.3f R_p   T = %.0f K'
              % (res['r_sonic'], res['T_sonic']))
    else:
        print('  critical (sonic) point   none in domain (max Mach %.3g)'
              % res['mach_max'])
    if res['r_exobase'] is not None:
        print('  exobase (Kn_bulk = 1)  %8.3f R_p' % res['r_exobase'])
    else:
        print('  exobase (Kn_bulk = 1)    above %.2f R_p (Kn_top = %.3g)'
              % (res['r_top'], res['Kn_top']))
    if res['r_kn_threshold'] is not None:
        print('  Kn_bulk = %.2g at        %8.3f R_p   (collisional below)'
              % (res['kn_threshold'], res['r_kn_threshold']))
    else:
        print('  Kn_bulk stays below %.2g over the whole domain'
              % res['kn_threshold'])
    print('  Kn_bulk at T max           %.3g' % res['Kn_at_Tmax'])
    print('  max Kn in critical region  %.3g   (bulk %.3g)'
          % (res['Kn_crit_max'], res['Kn_crit_bulk_max']))
    print('  Kn by species there: ' + ', '.join(
        '%s %.2g' % (s, res['Kn_crit_species'][s])
        for s in COLLIDERS if s in res['Kn_crit_species']))
    print('  max tau_E(e-ion)/tau_heat in critical region  %.3g'
          % res['tau_E_over_tau_heat_crit'])
    print('  log10 Mdot(hydro) = %.3f   formal Mdot_Jeans = %.3g g/s (%s)'
          % (res['log10_Mdot'], res['Mdot_jeans'], res['jeans_at']))
    print()
    print('  ' + res['validity_statement'].replace('  ', '\n  '))
    print()


def write_profile(res, path):
    cols = ['r', 'T', 'v', 'cs', 'L', 'lam_bulk', 'Kn_bulk']
    data = [res['r'], res['T'], res['v'], res['cs'], res['L'],
            res['lam_bulk'], res['Kn_bulk']]
    for s in COLLIDERS:
        if s in res['Kn']:
            cols += ['lam_' + s, 'Kn_' + s]
            data += [res['lam'][s], res['Kn'][s]]
    hdr = ('collisional_validity.py profile for %s\ncolumns '
           % res['case']) + ' '.join(cols) \
        + '\nr[Rp] T[K] v,cs[cm/s] L,lam[cm] Kn[-]'
    np.savetxt(path, np.column_stack(data), header=hdr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('case', nargs='+', help='EXHALE run directory')
    ap.add_argument('--adv', action='store_true',
                    help='read the advection-corrected _adv profiles instead'
                         ' of the solution')
    ap.add_argument('--outdir', default='output',
                    help='output subdirectory of the case (default output)')
    ap.add_argument('--kn-threshold', type=float, default=0.1)
    ap.add_argument('--json', default=None)
    ap.add_argument('--profile', default=None,
                    help='write the radial profile of the first case here')
    a = ap.parse_args(argv)

    out = []
    for c in a.case:
        res = collisional_diagnosis(c, adv=a.adv,
                                    kn_threshold=a.kn_threshold,
                                    outdir=a.outdir)
        report(res)
        if a.profile and not out:
            write_profile(res, a.profile)
        out.append(res)

    if a.json:
        with open(a.json, 'w') as f:
            json.dump([{k: r[k] for k in _SCALAR_KEYS
                        if r.get(k) is not None} | {'case': r['case'],
                       'validity_statement': r['validity_statement'],
                       'Kn_crit_species': r['Kn_crit_species'],
                       'lambda_J': r['lambda_J']}
                       for r in out], f, indent=2)
    return 0


if __name__ == '__main__':
    sys.exit(main())
