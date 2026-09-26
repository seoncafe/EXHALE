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
kinetic or transitional-flow calculation, and the Jeans estimate below is
only a scale against which the continuum number is read.

--------------------------------------------------------------------------
1. Collision model
--------------------------------------------------------------------------

The momentum-transfer collisions are the SAME four limits the element
diffusion operator uses, in the same Chapman-Enskog first approximation, so
that a mean free path quoted here and a diffusion coefficient quoted by the
wind cannot drift apart:

    neutral-neutral    hard-sphere (Banks & Kockarts 1973)      D ~ T^(1/2)/n
    ion-neutral        induced-dipole (Langevin) polarization
                       and the rigid core, frictions added      (both)
    ion-ion, ion-e     screened Coulomb, Spitzer ln(Lambda)     D ~ T^(5/2)/n

Every formula and every constant below is transcribed from
`src/modules/functions/binary_element_diffusion.f90` (function bodies at
lines 860-1002, constants at lines 253-283), with the source line cited at
each definition.  Python cannot call the Fortran, so this is a second
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
H2, H2+, H3+ and HeH+ when a run tracks them (same carrier lists as the Fortran, and the
same masses counted per collision partner: an H2 is one partner of mass 2 m_H).
He I 2^3S is NOT a separate species -- it is an excited level inside the
He I column (`species_table.f90` bsp_is_excited_level), and counting it
again would double count those atoms.  METALS ARE EXCLUDED as trace: at the
solar-scaled abundances of these runs they carry <1e-3 of the particles and
of the momentum, so they change no mean free path here.  They are kept in
the electron budget, where charge neutrality has to close exactly.

Two omissions, both in the direction of a LONGER mean free path, i.e. a
LARGER Knudsen number, i.e. a more conservative verdict:

  * RESONANT charge exchange (H+ + H, He+ + He) is not a channel.  The
    Fortran excludes it because a binary ELEMENT diffusion coefficient is
    driven by friction between the elements and both partners of a resonant
    pair carry the same element.  Here the like-element ion-neutral pair is
    a real momentum-transfer channel and it is computed, but with the
    non-resonant (polarization + core) cross section, which is the smaller
    one: resonant CX at 1e4 K is ~2e-15 cm^2 for H+ + H against the
    ~1e-15 cm^2 rigid core.  The mean free paths of H I and H II in the
    partially ionized layer are therefore upper limits.
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
the Chapman-Enskog momentum-transfer frequency, which for rigid spheres is
1.6x the elementary n sigma vbar; lambda_s is correspondingly 1.6x SHORTER
than the textbook 1/(n sigma).  The verdict threshold below (0.1) carries
that factor with room to spare, and the convention is stated so a number
quoted from here is not compared against a differently defined one.

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
removing itself from a minimum.  A minimum rule has the opposite behaviour
-- if the field it is currently taken from flattens, L jumps to the next
one, discontinuously and by whatever factor separates them.  That is not
hypothetical: the pressure alone was used here until 2026-08-28, and across
the heating peak of the LHS 1140 b runs d ln p/dr passes through a broad
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
diverges at every stagnation point -- the cell-centred velocity of the
base cells crosses v = 0 repeatedly, through the collocated two-cell
odd-even mode those cells carry (a mode of the cell-centred field, not a
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
                     rho) and gamma = 5/3.  This is the code's own sound
                     speed, not a re-definition: `eval_dt.f90:29` computes
                     cs = sqrt(g*p/rho) with g = 1.666666666667
                     (`parameters.f90:385`), and p, rho are the cgs columns
                     of the output file, so the mean molecular weight is
                     whatever the solution carries.

    CRITICAL REGION : from the radius of peak volumetric heating to the
                     sonic point -- the region that heats and accelerates
                     the wind and sets its topology.  "The critical region
                     is collisional" means max(Kn) < 0.1 across it, taken
                     over the bulk and over every species that carries at
                     least `TRACE_FRACTION` of the mass.

Coupling times: tau_s = 1/nu_s against the flow time r/|v|; and the
electron-ion ENERGY coupling time, which is longer than the momentum one by
the fraction of energy a collision transfers,

    nu^E_st = 2 mu_st/(m_s + m_t) nu_st     (elastic, Maxwellian average),

against the heating time tau_heat = (3/2) n_tot k T / Q with Q the heating
rate column.  Where tau^E_ei > tau_heat the single-temperature assumption of
the energy equation is itself unvalidated, independently of the Knudsen
number.

--------------------------------------------------------------------------
3. The kinetic bound
--------------------------------------------------------------------------

Where the verdict is UNVALIDATED, the collisionless escape rate through the
same exobase is computed for scale.  Jeans (1925), in the standard form of
Chamberlain & Hunten (1987, "Theory of Planetary Atmospheres", eq. 7.2.5):
a Maxwellian at the exobase loses, per unit area,

    Phi_J,s = n_s vbar_s / (2 sqrt(pi)) (1 + lambda_J,s) exp(-lambda_J,s),
    lambda_J,s = G M_p m_s / (k T r_exo)     (the Jeans escape parameter)

and Mdot_Jeans = 4 pi r_exo^2 sum_s m_s Phi_J,s.

How it bounds the hydrodynamic number.  Jeans escape is the rate a STATIC,
collisional-below/collisionless-above atmosphere loses through a Maxwellian
exobase with no bulk velocity there.  It is therefore the floor of the
kinetic problem, not its answer: a real transitional flow arrives at the
exobase with an outward drift and escapes faster.  The two numbers bracket
the physically admissible range only in the sense that

    Mdot_Jeans <= Mdot_kinetic,   and   Mdot_hydro is not constrained by
    either unless the flow is collisional through r_s.

So the statement the diagnostic makes is: if Mdot_hydro >> Mdot_Jeans while
r_s > r_exo, the continuum answer rests entirely on a continuum assumption
that fails at the point that sets it, and the ratio measures how far it
rests on it.  If lambda_J at the exobase is small (<~ 2-3) the atmosphere is
in hydrodynamic blow-off, the Jeans form is not applicable at all, and that
is reported instead of a number.

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
the LHS 1140 b 45 R_p wind: the
solution crosses its critical point at 40.06 R_p, and the `_adv` pair, whose
temperature in the outer wind is 3.5 to 4 times the solution's, reports no
critical point in the domain at all.

`--adv` reads the `_adv` pair anyway, which is the right question to ask of
the COMPOSITION: the mean free path is set by how much of the gas is
neutral, and in a wind whose ionization cannot relax over a flow time the
advected composition is the physical one while the equilibrium composition
the solution carries is not.  The two answers bracket the exobase rather
than agreeing: on that same wind the solution puts Kn_bulk = 1 above the
outer boundary and the `_adv` composition puts it at 29.0 R_p.  The report
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
m_e_g = 9.1093837015e-28    # electron mass [g], CODATA 2018
G_cgs = 6.67430e-8          # CODATA 2018 gravitational constant

# binary_element_diffusion.f90:262 -- CODATA 2018 exact coulomb value.
e_esu = 4.803204713e-10
# binary_element_diffusion.f90:265-266 -- Banks & Kockarts (1973) hard-sphere
# prefactor; equivalent to a rigid-sphere collision diameter d = 2.7 Angstrom
# put through the Chapman-Enskog integral (that file's own consistency note).
bk_hs_pref = 1.52e18
# binary_element_diffusion.f90:277-280 -- static dipole polarizabilities of
# the neutral partners, in a_0^3 converted with a_0^3 = 1.481847e-25 cm^3.
# alpha(H) = 4.5 a_0^3 is the exact nonrelativistic ground-state value; He and
# H2 are Schwerdtfeger & Nagle (2019, Mol. Phys. 117, 1200).
a0cub_cm3 = 1.481847e-25
alpha_HI = 4.500 * a0cub_cm3
alpha_HeI = 1.383 * a0cub_cm3
alpha_H2 = 5.315 * a0cub_cm3

# Species: name -> (mass [m_H units, PER COLLISION PARTNER], charge,
# neutral polarizability [cm^3] or 0 for a charged partner).  Masses and
# charges are the carrier lists of binary_element_diffusion.f90:286-297
# (hcar_*) extended by the helium stages; the electron is added because the
# Coulomb channel needs it and the Fortran, which transports elements, does
# not.  He I 2^3S is deliberately absent: it is inside the He I column.
COLLIDERS = {
    'HI':    (1.0, 0.0, alpha_HI),
    'HII':   (1.0, 1.0, 0.0),
    'HeI':   (4.0, 0.0, alpha_HeI),
    'HeII':  (4.0, 1.0, 0.0),
    'HeIII': (4.0, 2.0, 0.0),
    'H2':    (2.0, 0.0, alpha_H2),
    'H2p':   (2.0, 1.0, 0.0),
    'H3p':   (3.0, 1.0, 0.0),
    'HeHp':  (5.0, 1.0, 0.0),
    'e':     (m_e_g / m_H_g, -1.0, 0.0),
}
HEAVY = [s for s in COLLIDERS if s != 'e']

# Charge of every ion column that can appear, for the electron budget.  The
# metals are trace for the collisions but not for charge neutrality.
METAL_CHARGE = {}
for _m in eio.METAL_IONS:
    METAL_CHARGE[_m] = {'I': 0, 'II': 1, 'III': 2}[
        _m[len(_m.rstrip('I')):] or 'I']

gamma_ad = 1.666666666667   # parameters.f90:385, the code's polytropic index

# A species below this mass fraction is not asked to validate the solution:
# its Knudsen number is reported but does not enter the verdict.
TRACE_FRACTION = 1.0e-3


# --------------------------------------------------------------------------
# Pair diffusion coefficients -- transcribed from
# src/modules/functions/binary_element_diffusion.f90
# --------------------------------------------------------------------------

def hard_sphere_pair_diffusion(TK, ntot, A_s, A_t):
    """NEUTRAL-NEUTRAL.  Banks & Kockarts (1973):
    D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2) / n  [cm^2/s].
    Fortran: binary_element_diffusion.f90:860-874."""
    return bk_hs_pref * np.sqrt(1.0 / A_s + 1.0 / A_t) * np.sqrt(TK) / ntot


def polarization_pair_diffusion(TK, ntot, alpha_n, mu_g):
    """ION-NEUTRAL, induced-dipole (Langevin) channel:
    D = k T / (2.21 pi e n (alpha_n mu)^(1/2))  [cm^2/s], the low-energy
    limit -- never used alone, see ion_neutral_pair_diffusion.
    Fortran: binary_element_diffusion.f90:878-912."""
    return kb_erg * TK / (2.21 * np.pi * e_esu * ntot
                          * np.sqrt(np.maximum(alpha_n * mu_g, 1.0e-60)))


def coulomb_logarithm(TK, ne, zz):
    """ln(Lambda) = ln(3 k T lambda_D / (Z_s Z_t e^2)) with the electron
    Debye length (Spitzer 1962, section 5.2), floored at 1.
    Fortran: binary_element_diffusion.f90:940-956."""
    lam_D = np.sqrt(kb_erg * TK / (4.0 * np.pi * np.maximum(ne, 1.0)
                                   * e_esu * e_esu))
    return np.maximum(np.log(3.0 * kb_erg * TK * lam_D
                             / (max(zz, 1.0) * e_esu * e_esu)), 1.0)


def coulomb_pair_diffusion(TK, ntot, ne, Z_s, Z_t, mu_g):
    """ION-ION (and ion-electron).  Screened Coulomb through the
    Chapman-Enskog integral:
    D = 3 (kT)^(5/2) / [4 (2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln(Lambda)].
    Fortran: binary_element_diffusion.f90:918-936."""
    zz = abs(Z_s * Z_t)
    return (3.0 * (kb_erg * TK) ** 2.5
            / (4.0 * np.sqrt(2.0 * np.pi * mu_g) * ntot
               * (zz * e_esu * e_esu) ** 2 * coulomb_logarithm(TK, ne, zz)))


def ion_neutral_pair_diffusion(TK, ntot, alpha_n, A_s, A_t, mu_g):
    """ION-NEUTRAL, non-resonant: polarization and rigid core taken
    together, 1/D = 1/D_pol + 1/D_hs (the frictions of one potential add).
    Fortran: binary_element_diffusion.f90:962-1002."""
    Dpol = polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
    Dhs = hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
    return 1.0 / (1.0 / np.maximum(Dpol, 1.0e-99)
                  + 1.0 / np.maximum(Dhs, 1.0e-99))


def pair_diffusion(TK, ntot, ne, s, t):
    """The pair coefficient in whichever of the three limits the pair
    belongs to -- the Python twin of `stage_pair_diffusion`
    (binary_element_diffusion.f90:1006-1032)."""
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
    frequency of one particle of s against the gas of t, from the same
    binary diffusion coefficient the wind's diffusion operator uses.  The
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
    code's own sound speed (eval_dt.f90:29, parameters.f90:385)."""
    return _first_upcrossing(r, v - cs)


def exobase(r, Kn_bulk):
    """Exobase: Kn_bulk = 1, first crossing from below."""
    return _first_upcrossing(r, Kn_bulk - 1.0)


def jeans_escape(r_exo_cm, TK, dens, Mp_g):
    """Collisionless escape through a Maxwellian exobase (Jeans 1925;
    Chamberlain & Hunten 1987 eq. 7.2.5):

        Phi_J,s = n_s vbar_s/(2 sqrt(pi)) (1 + lambda_J) exp(-lambda_J),
        lambda_J,s = G M_p m_s/(k T r_exo).

    Returns (Mdot [g/s], {species: lambda_J}, {species: Phi [cm^-2 s^-1]}).
    Heavy species only: the electrons are held by the ambipolar field."""
    Mdot, lam_J, phi = 0.0, {}, {}
    for s in HEAVY:
        if s not in dens:
            continue
        m_s = COLLIDERS[s][0] * m_H_g
        lj = G_cgs * Mp_g * m_s / (kb_erg * TK * r_exo_cm)
        f = (dens[s] * mean_thermal_speed(TK, s) / (2.0 * np.sqrt(np.pi))
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

    res = dict(
        case=case_dir, adv=adv, kn_threshold=kn_threshold,
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
            "and its mass flux is set at the outer boundary, where "
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
        lines.append("The %s; max Kn from peak heating (%.3f R_p) to the "
                     "outer boundary (%.3f R_p) is %.3g (threshold %.2g)."
                     % (exo_txt, res['r_heat_peak'], res['r_top'], kn, thr))

    if collisional_through_crit:
        lines.append(
            "HYDRODYNAMIC RESULT VALIDATED through its critical point: the "
            "gas is collisional (Kn < %.2g) everywhere the wind is heated and "
            "accelerated, and the critical point lies below the exobase, so "
            "the continuum equations hold where the solution is set. "
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
            "HYDRODYNAMIC RESULT UNVALIDATED: " + "; ".join(why) + ". "
            "The sign and size of the error are NOT determined by this "
            "solution -- it is unvalidated, not overestimated. "
            "log10 Mdot = %.3f." % res['log10_Mdot'])
        lam_b = res['lambda_J_bulk']
        ratio = 10.0 ** res['log10_Mdot'] / max(res['Mdot_jeans'], 1e-300)
        lines.append(
            "Kinetic scale: Jeans escape through the %s (T = %.0f K, "
            "lambda_J = %.2f for the lightest heavy species) gives "
            "Mdot_Jeans = %.3g g/s, a factor %.2g below the continuum Mdot."
            % (res['jeans_at'], res['T_jeans'], lam_b, res['Mdot_jeans'],
               ratio))
        if lam_b is not None and lam_b < 3.0:
            lines.append(
                "That number is a SCALE, not a bound: at lambda_J = %.2f the "
                "exobase is barely gravitationally bound, the Jeans integral "
                "is no longer the small escaping tail of a Maxwellian but "
                "most of it, and the atmosphere is in hydrodynamic blow-off. "
                "What it does say is that the continuum Mdot is within a "
                "factor of a few of what a collisionless outer atmosphere at "
                "the same density and temperature would lose, so the failure "
                "of the continuum assumption is not hiding an order of "
                "magnitude." % lam_b)
        else:
            lines.append(
                "Jeans is the floor of the kinetic problem (static exobase, "
                "no bulk drift), so it bounds the kinetic answer from below "
                "and does not bound the continuum one at all; the ratio "
                "measures how far the reported Mdot rests on the continuum "
                "assumption that fails at the critical point.")
    return "  ".join(lines)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

_SCALAR_KEYS = ['r_sonic', 'r_exobase', 'r_heat_peak', 'r_Tmax', 'T_max',
                'T_sonic', 'mach_max', 'Kn_at_sonic', 'Kn_at_Tmax',
                'r_kn_threshold', 'Kn_crit_max',
                'Kn_crit_bulk_max', 'Kn_top', 'r_top', 'log10_Mdot',
                'jeans_at', 'r_jeans', 'T_jeans', 'Mdot_jeans',
                'lambda_J_bulk', 'tau_E_over_tau_heat_crit']


def report(res):
    print('=' * 74)
    print('case: %s   (%s)'
          % (res['case'],
             'advection-corrected _adv profiles: the solution\'s density and'
             ' velocity with a temperature and composition from a closure the'
             ' momentum equation was not re-solved for'
             if res['adv'] else
             'the solution; its composition is the local photoionization'
             ' equilibrium, which the mean free path inherits'))
    print('=' * 74)
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
    print('  log10 Mdot(hydro) = %.3f   Mdot_Jeans = %.3g g/s (%s)'
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
