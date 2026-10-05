#!/usr/bin/env python3
"""Steady temperature of the lower-atmosphere column under the heat the
escape wind conducts into it through the matching level.

THE PROBLEM.  The wind solved by EXHALE above the matching level p_match
conducts a heat flux F [erg cm^-2 s^-1] down through its base face
(`base_conductive_flux_cgs` of EXHALE_resolved.out).  The column below is the
climate solution of `radiative_convective_column.py`, whose stratosphere is
an isotherm at the skin temperature and which has no energy equation of its
own that could take that flux.  This module solves the steady energy balance
of the column below p_match with the flux as its upper boundary condition:

    (1/r^2) d/dz ( r^2 kappa(T) dT/dz ) = C(T, p) - C(T0(p), p) ,

    kappa dT/dz = F            at p = p_match ,
    T = T0(p_bottom)           at the bottom of the domain ,

with T0(p) the climate temperature (the state the column is in without the
flux), kappa the thermal conductivity of the He/H2 mixture and C the net
radiative cooling per unit volume.  Only the DIFFERENCE C(T) - C(T0) enters,
so whatever heating holds the climate state (absorbed starlight, upwelling
infrared) is held at its climate value, and F = 0 returns T = T0 exactly.

COORDINATE.  The equation is solved in s = -ln p at fixed pressure levels.
A layer between two levels holds the column mass dp/g whatever its
temperature; its thickness dz = H ds, H = k T/(m g), follows T.  With the
conductive flux q = kappa dT/dz = (kappa/H) dT/ds (positive when heat flows
DOWN, the direction F enters), the balance is

    d(A q)/ds = A (p/(m g)) sum_X q_X [ e_X(T) - e_X(T0) ] ,

A = r^2 the shell area, q_X the number fraction of the emitter X and e_X its
net cooling per molecule [erg s^-1].  Discretized conservatively on nodes
s_0 (bottom) .. s_M (= -ln p_match), the discrete sum of the right-hand side
equals A_M F - A_{1/2} q_{1/2} exactly at the solution, which is the energy
budget the tests check.  r is the reference profile's radius and is held: a
warming of a few hundred kelvin over a few scale heights moves r by ~1e-3 of
itself, which changes A by ~2e-3.

CONDUCTIVITY.  The mixture rule and the coefficients of
`src/modules/time_step/viscous_conduction.f90` (Eqs. 4' and 4''): kappa =
x_He 299 T^0.69 + x_H 379 T^0.69 + x_H2 kappa_H2(T), the number fractions
over all particles; the trace molecules carry no conductivity of their own.
The same conductivity the wind uses at its base face, so the flux leaving the
wind and the flux entering the column are one law.

EMISSION.  Optically thin band emission of CH4, NH3, H2O, C2H2, CO and CO2
from the k-distributions the climate model itself reads
(photochem_clima_data, HELIOS-K tables from HITRAN/HITEMP): in LTE the
emission per molecule in a wavelength bin is 4 pi sigma_bar(T, p) Bbar(T),
sigma_bar = sum_g w_g k_g the bin-mean cross section and Bbar the Planck
function integrated over the bin.  (The product of two bin means is not the
mean of the product; a bin of the table is 50-100 cm^-1 wide in the
mid-infrared, and the error that makes is measured against the two-level band
formula in the report of `band_strength_report`.)

LTE AND NON-LTE.  Every bin is assigned to one band of its species
(`BANDS`).  Rotational bands are in LTE: the NH3 rotation-inversion lines
have A <= 4.8 s^-1 while their pressure-broadening (collision) rate is
1.4e4 s^-1 at 1 microbar (`ROTATIONAL_LTE_NOTE`).  A vibrational band is a two-level system,

    e = eps (e_LTE(T) - a) ,   eps = C10/(C10 + A10) ,
    C10 = sum_M k_M(T) n_M      (M = He, H2) ,

with A10 the band's Einstein coefficient, k_M the published V-T
de-excitation rate coefficient and a the power absorbed per molecule from the
ambient field (stimulated emission folded into the LTE emission the table
carries).  Derivation: with n1 (A10 + C10) = n0 (C01 + R01) and detailed
balance C01 = C10 (g1/g0) exp(-h nu/k T), the net collisional cooling
h nu (n0 C01 - n1 C10) is eps (h nu n0 C01 A10/C10 - h nu n0 R01), i.e. eps
times (LTE emission minus absorption).  The ambient field in a vibrational
band is taken as the upward intensity of the climate column's isothermal
stratosphere, optically thick in the band core and at T0 there:
J = B(T0)/2, so a = e_LTE(T0)/2 (ABSORBED_FRACTION).  The rates are the
published ones (section PUBLISHED RATES): CH4 dyad + He and + H2 measured,
NH3 nu2 + He measured and + H2 ESTIMATED from an analog (stated there).  A
band without rates for both He and H2 stays in LTE (an upper bound on its
emission) and is flagged; 'nlte_partial' uses the measured rates only and
`treatment='rotational'` gives the lower bound on the cooling
(`TREATMENTS`).

OPTICAL DEPTH, TRANSFER 'thin' (default).  The thin treatment needs the
photons emitted upward to escape.  The test is on the fraction of a bin's
emission excess, weighted over the levels that emit it, that the column
above absorbs on the way up.  Photons emitted isotropically inside a layer
escape upward with the probability of the cooling-to-space transfer below,
the layer integral of E_2(tau) over the layer's optical depths and the
k-distribution g-points (Dickinson 1972, Eq. 5 for the transmission to
space; `escape_factors`), so that the absorbed share of the upward photons
is 2 (1 - beta), beta the cooling-to-space escape of the layer.  (The flux
transmission exp(-1.66 tau), the diffusivity form of 2 E_3(tau), is the
transmission of a beam emerging from the column, not the escape of photons
emitted in it: at tau = 0.04 it absorbs 6.4 per cent where 1 - E_2(tau)
absorbs 14.6.)  The thin treatment is VALID while at most 10 per cent of
that excess is absorbed (THIN_ESCAPE_MIN = 0.9), in any bin that carries
more than 1 per cent of the column's emission excess
(`optical_depth_check`).  Past that limit the solution is kept and FLAGGED
(res['thin_layer']: exceeded, the largest absorbed fraction and its bin;
user decision 2026-10-04), not refused; callers record the flag with the
state.  The emission-weighted vertical optical depth of the column above
each level, tau_bar = sum_g w_g k_g tau_g / sum_g w_g k_g, and the largest
single g-point value are reported with it.  Photons emitted downward are
absorbed in the optically thick stratosphere below; that energy is
deposited where the column's own radiative budget (~1e5 erg cm^-2 s^-1 of
absorbed starlight) dwarfs F, and is not followed.

OPTICAL DEPTH, TRANSFER 'cooling_to_space'.  The absorption of the band
emission by the overlying column is followed in the cooling-to-space
approximation (Rodgers & Walshaw 1966, as stated and generalized to a
two-level band by Dickinson 1972, J. Atmos. Sci. 29, 1531, section 3c Eq. 12
and section 3g Eqs. 30-33; the published version was read).  Dickinson's
Eq. (33) for the net rate of a band per molecule,

    R_H = [1 + (1/2)(A/lambda) T(z,inf)]^-1 [R_ext + A E(z)
                                             - (1/2) A T(z,inf) g e^(-h nu/kT)],

with T(z,inf) the angle-integrated band transmission to space (his Eq. 5),
lambda the collisional de-excitation rate and E(z) = (1/2) int (X(z') - X(z))
dT(z,z') the exchange function (his Eq. 30), is applied with E(z) split in
two parts:
  * toward the column BELOW: photons emitted downward are absorbed in the
    optically thick climate column, whose source function is held at its
    climate value X0.  That part of E is (1/2)(X0 - X(z)): it adds the
    downward half 1/2 to the escape and the field from below, absorbed at
    the climate value, to R_ext (a = ABSORBED_FRACTION e_LTE(T0), as in the
    thin treatment);
  * toward the column ABOVE and inside the warmed layer: neglected (the
    cooling-to-space approximation).  Above p_match the gas is the base of
    the wind, at least as warm as the match (the heat flows down), so
    photons it absorbs are not a net loss of the column; inside the domain
    the neglected exchange is measured on every solution (below).
The escape probability of node i in bin b is then

    beta_b,i = 1/2 + sum_g w_g k_g u_g,i / sum_g w_g k_g ,
    u_g,i = [E_3(tau_t) - E_3(tau_b)] / (2 (tau_b - tau_t)) ,

u the fraction of the emission of the node's control layer (vertical optical
depths tau_t < tau_b of its faces, measured from the top of the profile) that
escapes to space, the layer integral of Dickinson's (1/2) T(z,inf) for one
g-point, T = E_2(tau).  Optical depths and weights are those of the climate
state (absorption held at T0, as in the thin treatment).  LTE band:
excess = sum_bins beta (e_bin(T) - e_bin(T0)).  Two-level band: escape beta =
(sum_bins beta e_bin)/e_band, eps_beta = C10/(C10 + beta A10) (Dickinson's
prefactor; a band with beta A10 << C10 thermalizes), excess =
eps_beta(T)(beta e(T) - a) - eps_beta(T0)(beta e(T0) - a).  With no column
(beta = 1) this is the thin treatment exactly.
VALIDITY.  The neglected exchange inside the domain (photons of the warmed
nodes absorbed by other warmed nodes, Dickinson's E(z) restricted to the
domain) is evaluated on the solution with exact slab exchange fractions
(E_3 of the optical-depth differences between faces, every g-point) over the
nodes warmed by more than 1e-6 of the rise at the match, and its effect on
T_match is computed to first order with the solution's Jacobian.  The solve is
REFUSED (ColumnTransferRefusal) when that shift exceeds exchange_shift_max_K.
For a two-level band the exchanged excess is weighted by eps_beta at both
ends (the population excess of a two-level band is eps_beta times the LTE
one when the absorbed field is held); exact for LTE bands.
`solve_with_exchange` solves the balance with the full exchange of the whole
domain (fixed-point iteration on the exchange heating) to test it.
THE OVERLYING GAS.  Both transfers take the column above p_match from the
lower-atmosphere profile, i.e. the Photochem CH4 and NH3 above the match.
Above the match the gas is the wind's, which carries C, N and O as atoms
and is warmer, so that absorber is not in the coupled model.  On LHS 1140 b
(deep He/H 1.8, 1.5 and 1.2, F = 11.4, 13.8 and 16.5 erg cm^-2 s^-1) it
makes 87-88 per cent of the thin test's absorbed fraction (22.2 -> 2.9, 25.2
-> 3.1 and 27.9 -> 3.3 per cent without it) and 95-96 per cent of the
cooling-to-space rise of T_match (+2.42 -> +0.12, +3.08 -> +0.13 and +3.68
-> +0.14 K) (MEASURED, include_gas_above_match = False, update log stage 3).

WHAT IS HELD.  Composition (number fractions at fixed pressure), the
absorption, r(p), gravity.  The caller iterates composition and temperature
(`photochem_to_lower_profile.py --conducted-heat-flux`).
"""
import os

import numpy as np

KB = 1.380649e-16           # erg/K
HP = 6.62607015e-27         # erg s
CL = 2.99792458e10          # cm/s
C2 = HP*CL/KB               # second radiation constant, cm K
AMU = 1.66053906660e-24
GNEWT = 6.67430e-8
BAR = 1.0e6

# --------------------------------------------------------------------------
# Thermal conductivity: viscous_conduction.f90, Eqs. (4') and (4'')
# --------------------------------------------------------------------------
KAPPA_H_COEF = 379.0        # Salz et al. 2015, Eq. 25
KAPPA_HE_COEF = 299.0       # Sutton et al. 2015, Eq. 5
KAPPA_N_EXPO = 0.69
KAPPA_H2_A, KAPPA_H2_S, KAPPA_H2_C, THETA_V_H2 = 272.523, 0.735652, 0.370446, 5987.0


def einstein_heat_capacity(x):
    """x^2 e^-x/(1 - e^-x)^2, the form of viscous_conduction.f90."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 0.05
    x2 = x[small]**2
    out[small] = 1.0 - x2*(1.0/12.0 - x2*(1.0/240.0 - x2/6048.0))
    big = ~small
    em = np.exp(-np.minimum(x[big], 700.0))
    out[big] = np.where(x[big] > 700.0, 0.0, x[big]**2*em/(1.0 - em)**2)
    return out


def mixture_conductivity(T, x_He, x_H2, x_H):
    """erg cm^-1 s^-1 K^-1 of the neutral He/H2/H mixture (number fractions
    over all particles).  The H2 term is valid 200-2000 K (Eq. 4'' of
    viscous_conduction.f90; +7.6 per cent at 150 K)."""
    T = np.asarray(T, dtype=float)
    k_h2 = KAPPA_H2_A*T**KAPPA_H2_S*(1.0 + KAPPA_H2_C
                                     * einstein_heat_capacity(THETA_V_H2/T))
    return (x_He*KAPPA_HE_COEF*T**KAPPA_N_EXPO + x_H*KAPPA_H_COEF
            * T**KAPPA_N_EXPO + x_H2*k_h2)


# --------------------------------------------------------------------------
# The generic solver: conduction against a local net cooling, in s = -ln p
# --------------------------------------------------------------------------

class ColumnGeometry(object):
    """Nodes s_0 (deepest) .. s_M (the matching level) and what the balance
    needs at each: the shell area A_i, the column mass
    w_i = (p_i/g_i) (s_{i+1/2} - s_{i-1/2}) [g cm^-2] of the node's control
    interval (half an interval at the top node), and the factor
    m_i g_i/k that turns the temperature into the scale height."""

    def __init__(self, s, area, mass_weight, mg_over_k):
        self.s = np.asarray(s, dtype=float)
        self.area = np.asarray(area, dtype=float)
        self.w = np.asarray(mass_weight, dtype=float)
        self.mg_over_k = np.asarray(mg_over_k, dtype=float)
        self.M = self.s.size - 1
        if np.any(np.diff(self.s) <= 0.0):
            raise ValueError('the nodes must rise monotonically in s = -ln p')
        self.area_face = 0.5*(self.area[1:] + self.area[:-1])
        self.mg_face = 0.5*(self.mg_over_k[1:] + self.mg_over_k[:-1])
        self.ds = np.diff(self.s)


def control_interval_widths(s):
    """Width in s of each node's control interval, the top node holding the
    half interval below it and the bottom node (a Dirichlet node) the half
    interval above it."""
    s = np.asarray(s, dtype=float)
    half = 0.5*np.diff(s)
    width = np.zeros_like(s)
    width[:-1] += half
    width[1:] += half
    return width


def assemble_balance(geom, F, T, face_coefficient, net_cooling):
    """Residual R (nodes 0..M, R[0] unused), banded Jacobian on T_1..T_M and
    face fluxes q of the discrete balance of `solve_conducted_response`."""
    M = geom.M
    Af, ds = geom.area_face, geom.ds
    Tf = 0.5*(T[1:] + T[:-1])
    c, dc = face_coefficient(Tf)
    dT = T[1:] - T[:-1]
    q = c*dT/ds                                   # faces 1/2 .. M-1/2
    dq_dlo = (-c + 0.5*dc*dT)/ds                  # d q / d T_i
    dq_dhi = (c + 0.5*dc*dT)/ds                   # d q / d T_{i+1}
    D, dD = net_cooling(T)
    R = np.zeros(M + 1)
    # flux into node i from above minus out below
    top = np.zeros(M + 1)
    top[:M] = Af*q
    top[M] = geom.area[M]*F
    bot = np.zeros(M + 1)
    bot[1:] = Af*q
    R[1:] = (top - bot - geom.area*geom.w*D)[1:]
    # banded Jacobian on the unknowns T_1..T_M
    ab = np.zeros((3, M))
    diag = np.zeros(M + 1)
    upper = np.zeros(M + 1)                       # dR_i/dT_{i+1}
    lower = np.zeros(M + 1)                       # dR_i/dT_{i-1}
    # top face of node i (i < M): A q_{i+1/2}
    diag[:M] += Af*dq_dlo
    upper[:M] += Af*dq_dhi
    # bottom face of node i (i >= 1): -A q_{i-1/2}
    diag[1:] -= Af*dq_dhi
    lower[1:] -= Af*dq_dlo
    diag -= geom.area*geom.w*dD
    ab[1, :] = diag[1:]
    ab[0, 1:] = upper[1:M]
    ab[2, :-1] = lower[2:]
    return R, ab, q


def linear_temperature_response(geom, F, T, face_coefficient, net_cooling,
                                heating):
    """First-order change of the node temperatures when a heating
    `heating`_i [the units of A_i w_i D_i] is added to the balance at the
    solution T: delta T = -J^-1 heating, J the Jacobian of the residual."""
    from scipy.linalg import solve_banded
    _, ab, _ = assemble_balance(geom, F, T, face_coefficient, net_cooling)
    out = np.zeros_like(T)
    out[1:] = solve_banded((1, 1), ab, -np.asarray(heating, float)[1:])
    return out


def solve_conducted_response(geom, F, T_bottom, T_guess, face_coefficient,
                             net_cooling, tol=1e-10, max_iter=200,
                             verbose=False):
    """Newton solve of the discrete balance.

    face_coefficient(T_face) -> (c, dc/dT) per face: the conductive flux
        through face i+1/2 is q = c (T_{i+1} - T_i)/ds_i, positive downward
        (c = kappa m g/(k T) for the physical column).
    net_cooling(T) -> (D, dD/dT) per node: the net cooling EXCESS per unit
        column mass, D_i = sum_X (q_X/m) [e_X(T_i) - e_X(T0_i)]
        [erg g^-1 s^-1] for the physical column.

    Residual of node i = 1..M (node 0 is held at T_bottom):
        R_i = A_{i+1/2} q_{i+1/2} - A_{i-1/2} q_{i-1/2} - A_i w_i D_i ,
    with A_{M+1/2} q_{M+1/2} := A_M F at the top.  Returns T and the face
    fluxes; raises RuntimeError if Newton does not converge.
    """
    from scipy.linalg import solve_banded
    M = geom.M
    T = np.array(T_guess, dtype=float)
    T[0] = T_bottom

    def assemble(T):
        return assemble_balance(geom, F, T, face_coefficient, net_cooling)

    scale = max(abs(geom.area[M]*F), 1e-300)
    for it in range(max_iter):
        R, ab, q = assemble(T)
        res = np.max(np.abs(R[1:]))/scale
        if verbose:
            print('    newton %d  max|R|/(A F) = %.3e' % (it, res))
        if res < tol:
            return T, q, it
        dx = solve_banded((1, 1), ab, -R[1:])
        # damping: no node moves by more than 30 per cent of its temperature
        lim = np.max(np.abs(dx)/(0.3*T[1:]))
        if lim > 1.0:
            dx /= lim
        T[1:] += dx
        if np.any(~np.isfinite(T)) or np.any(T[1:] <= 0.0):
            raise RuntimeError('the column solve left the physical range')
    raise RuntimeError('the column solve did not converge in %d Newton '
                       'iterations (max|R|/(A F) = %.3e)' % (max_iter, res))


def discrete_energy_budget(geom, T, q, F, net_cooling):
    """(sum_i A_i w_i D_i over the free nodes, A_{1/2} q_{1/2}, A_M F): the
    emission excess of the column, the flux leaving through its bottom and
    the flux entering at its top.  At a solution the first two sum to the
    third to the Newton tolerance."""
    D, _ = net_cooling(T)
    excess = float(np.sum((geom.area*geom.w*D)[1:]))
    leak = float(geom.area_face[0]*q[0])
    return excess, leak, float(geom.area[-1]*F)


# --------------------------------------------------------------------------
# Band emission from the climate model's k-distributions
# --------------------------------------------------------------------------

EMITTERS = ('CH4', 'NH3', 'H2O', 'C2H2', 'CO', 'CO2')

# Each bin of a species' table is assigned to one band by the wavenumber of
# its centre [cm^-1].  'rot' bands are pure rotational and in LTE; every
# other band is vibrational and gets the two-level treatment when
# NLTE_BANDS carries published parameters for it, otherwise it stays in LTE
# and is FLAGGED as such in every report.  The windows follow the band
# origins (HITRAN band centres): NH3 nu2 932/968, nu4 1627, nu1/nu3
# 3337/3444; CH4 nu4 1306, nu2 1533, nu3 3019; H2O nu2 1595, nu1/nu3
# 3657/3756; C2H2 nu5 729, nu4+nu5 1328, nu3 3289; CO 2143; CO2 nu2 667,
# nu3 2349.  The rotational windows stop below the P branches of the lowest
# fundamentals.
BANDS = {
    'NH3': (('rotational', 0.0, 600.0, 'rot'),
            ('nu2', 600.0, 1250.0, 'vib'),
            ('nu4+2nu2', 1250.0, 2200.0, 'vib'),
            ('nu1+nu3 and above', 2200.0, 1.0e5, 'vib')),
    'CH4': (('rotational', 0.0, 1000.0, 'rot'),
            ('nu4+nu2 dyad', 1000.0, 1750.0, 'vib'),
            ('nu3 pentad and above', 1750.0, 1.0e5, 'vib')),
    'H2O': (('rotational', 0.0, 1100.0, 'rot'),
            ('nu2', 1100.0, 2300.0, 'vib'),
            ('nu1+nu3 and above', 2300.0, 1.0e5, 'vib')),
    'C2H2': (('rotational', 0.0, 545.0, 'rot'),
             ('nu5', 545.0, 1000.0, 'vib'),
             ('nu4+nu5', 1000.0, 1750.0, 'vib'),
             ('nu3 and above', 1750.0, 1.0e5, 'vib')),
    'CO': (('rotational', 0.0, 1000.0, 'rot'),
           ('fundamental and above', 1000.0, 1.0e5, 'vib')),
    'CO2': (('rotational', 0.0, 545.0, 'rot'),
            ('nu2', 545.0, 1000.0, 'vib'),
            ('nu3 and above', 1000.0, 1.0e5, 'vib')),
}


def kdistribution_directory():
    """The k-distribution files the installed climate model reads."""
    import photochem_clima_data
    return os.path.join(photochem_clima_data.DATA_DIR, 'kdistributions')


def planck_bin_integral(nu_lo, nu_hi, T, npts=48):
    """Integral of B_nu~ over [nu_lo, nu_hi] cm^-1, erg s^-1 cm^-2 sr^-1, by
    Gauss-Legendre in wavenumber; nu_lo, nu_hi arrays of bins, T an array.
    Returns (nbins, nT)."""
    x, wq = np.polynomial.legendre.leggauss(npts)
    nu_lo = np.asarray(nu_lo, float)[:, None]
    nu_hi = np.asarray(nu_hi, float)[:, None]
    nu = 0.5*(nu_hi + nu_lo) + 0.5*(nu_hi - nu_lo)*x[None, :]     # (nb, nq)
    T = np.asarray(T, float)
    arg = C2*nu[:, :, None]/T[None, None, :]
    B = 2.0*HP*CL**2*nu[:, :, None]**3/np.expm1(np.minimum(arg, 700.0))
    return 0.5*(nu_hi - nu_lo)*np.einsum('q,bqt->bt', wq, B)


class KDistributionTable(object):
    """One species' table: bin-mean cross sections sigma_bar = sum_g w_g k_g
    and the g-point cross sections, interpolated bilinearly (log10 k against
    T and log10 p) as the climate model's tables are tabulated."""

    def __init__(self, species, kdir=None):
        import h5py
        path = os.path.join(kdir or kdistribution_directory(), species + '.h5')
        with h5py.File(path, 'r') as h:
            self.log10k = np.array(h['log10k'])          # (bins, T, P, g)
            self.T = np.array(h['T'], dtype=float)
            self.log10P = np.array(h['log10P'], dtype=float)
            self.wg = np.array(h['weights'], dtype=float)
            wl = np.array(h['wavelengths'], dtype=float)  # um, bin edges
        self.species = species
        self.nu_hi = 1.0e4/wl[:-1]
        self.nu_lo = 1.0e4/wl[1:]
        self.nu_c = 0.5*(self.nu_hi + self.nu_lo)
        self.dnu = self.nu_hi - self.nu_lo
        # bins where a species has no lines carry sigma_bar = 0; the floor
        # keeps the log-space interpolation finite and returns ~0 there
        self.log10_sigma_bar = np.log10(np.maximum(
            np.sum(self.wg*10.0**self.log10k, axis=-1), 1e-99))  # (bins,T,P)

    def _weights(self, T, log10p):
        T = np.clip(np.asarray(T, float), self.T[0], self.T[-1])
        lp = np.clip(np.asarray(log10p, float), self.log10P[0],
                     self.log10P[-1])
        iT = np.clip(np.searchsorted(self.T, T) - 1, 0, self.T.size - 2)
        iP = np.clip(np.searchsorted(self.log10P, lp) - 1, 0,
                     self.log10P.size - 2)
        fT = (T - self.T[iT])/(self.T[iT+1] - self.T[iT])
        fP = (lp - self.log10P[iP])/(self.log10P[iP+1] - self.log10P[iP])
        return iT, iP, fT, fP

    def sigma_bar(self, T, log10p):
        """(bins, n) bin-mean cross section [cm^2] at n states."""
        iT, iP, fT, fP = self._weights(T, log10p)
        L = self.log10_sigma_bar
        v = ((1-fT)*(1-fP)*L[:, iT, iP] + fT*(1-fP)*L[:, iT+1, iP]
             + (1-fT)*fP*L[:, iT, iP+1] + fT*fP*L[:, iT+1, iP+1])
        return 10.0**v

    def k_g(self, T, log10p):
        """(bins, n, g) g-point cross sections [cm^2]."""
        iT, iP, fT, fP = self._weights(T, log10p)
        L = self.log10k
        fT = fT[None, :, None]
        fP = fP[None, :, None]
        v = ((1-fT)*(1-fP)*L[:, iT, iP, :] + fT*(1-fP)*L[:, iT+1, iP, :]
             + (1-fT)*fP*L[:, iT, iP+1, :] + fT*fP*L[:, iT+1, iP+1, :])
        return 10.0**v

    def band_of_bins(self):
        """Index into BANDS[species] of each bin."""
        out = np.empty(self.nu_c.size, dtype=int)
        for b, (_, lo, hi, _) in enumerate(BANDS[self.species]):
            out[(self.nu_c >= lo) & (self.nu_c < hi)] = b
        return out


class LTEBandEmission(object):
    """LTE optically thin emission per molecule [erg s^-1] of every band of
    every emitter, 4 pi sum_{bins in band} sigma_bar(T, p) Bbar_bin(T),
    tabulated on a fine temperature grid at the table's own pressure nodes
    and interpolated linearly in log10 p between them (the bin-mean cross
    sections vary by less than 2 per cent from 1e-6 to 1e-1 bar)."""

    T_GRID = np.arange(50.0, 2000.0 + 0.25, 0.5)

    def __init__(self, species=EMITTERS, kdir=None):
        self.tables = {}
        self.e_nodes = {}
        Tg = self.T_GRID
        for sp in species:
            tab = KDistributionTable(sp, kdir)
            self.tables[sp] = tab
            Bbar = planck_bin_integral(tab.nu_lo, tab.nu_hi, Tg)  # (bins, nT)
            band = tab.band_of_bins()
            nb = len(BANDS[sp])
            e = np.zeros((nb, tab.log10P.size, Tg.size))
            for iP, lp in enumerate(tab.log10P):
                sb = tab.sigma_bar(Tg, np.full(Tg.size, lp))      # (bins, nT)
                contrib = 4.0*np.pi*sb*Bbar
                for b in range(nb):
                    e[b, iP] = contrib[band == b].sum(axis=0)
            self.e_nodes[sp] = e

    def band_names(self, sp):
        return [b[0] for b in BANDS[sp]]

    def emission(self, sp, b, T, log10p):
        """e and de/dT of band b of species sp at states (T, log10 p)."""
        tab = self.tables[sp]
        lp = np.clip(np.asarray(log10p, float), tab.log10P[0], tab.log10P[-1])
        iP = np.clip(np.searchsorted(tab.log10P, lp) - 1, 0,
                     tab.log10P.size - 2)
        fP = (lp - tab.log10P[iP])/(tab.log10P[iP+1] - tab.log10P[iP])
        Tg = self.T_GRID
        T = np.asarray(T, float)
        if np.any(T < Tg[0]) or np.any(T > Tg[-1]):
            raise ValueError('temperature outside the k-distribution range '
                             '%.0f-%.0f K' % (Tg[0], Tg[-1]))
        e = self.e_nodes[sp][b]
        j = np.clip(np.searchsorted(Tg, T) - 1, 0, Tg.size - 2)
        fT = (T - Tg[j])/(Tg[j+1] - Tg[j])
        lo = (1 - fP)*e[iP, j] + fP*e[iP+1, j]
        hi = (1 - fP)*e[iP, j+1] + fP*e[iP+1, j+1]
        return lo + fT*(hi - lo), (hi - lo)/(Tg[j+1] - Tg[j])


# --------------------------------------------------------------------------
# Two-level non-LTE parameters of the vibrational bands
# --------------------------------------------------------------------------
# One entry per (species, band) whose V-T de-excitation rate coefficients
# are known (section "PUBLISHED RATES" below):
#   'A10'          band Einstein coefficient [s^-1] as a function of T;
#   'k_vt'         {partner: function T [K] -> k10 [cm^3 s^-1]}, MEASURED
#                  rates, each cited with its temperature range;
#   'k_vt_estimated' {partner: function}, a partner with no measurement,
#                  estimated by a stated argument (used by 'nlte', left out
#                  by 'nlte_partial');
#   'source'       one line naming the papers.
NLTE_BANDS = {}


class LandauTellerRate(object):
    """A published rate coefficient tabulated at temperatures T_i, k_i.

    Inside the tabulated range: ln k interpolated linearly in T^(-1/3)
    between neighbouring points, which reproduces every published value
    exactly.  Outside it: the Landau-Teller law ln k = a + b T^(-1/3)
    continued from the nearest published point, with the slope b of the
    least-squares fit to all the published points of that pair.  Basis: the
    T^(-1/3) dependence is the Landau-Teller form of a V-T rate (the form
    Halthore et al. 1994, Eq. 1, fit to CH4 relaxation over 150-500 K); it is
    an EXTRAPOLATION and every entry states how far it is used beyond its
    data.  `relaxation_to_deexcitation(T)` divides a measured relaxation
    rate (the decay rate of the excess population, k10 + k01) by
    1 + sum_u (g_u/g_0) exp(-E_u/kT) to give k10; None when the paper
    already states the elementary de-excitation rate."""

    def __init__(self, T, k, relaxation_to_deexcitation=None):
        o = np.argsort(T)
        self.T = np.asarray(T, float)[o]
        self.lnk = np.log(np.asarray(k, float)[o])
        self.x = self.T**(-1.0/3.0)
        A = np.vstack([np.ones_like(self.x), self.x]).T
        self.b = float(np.linalg.lstsq(A, self.lnk, rcond=None)[0][1])
        self.correction = relaxation_to_deexcitation

    def __call__(self, T):
        T = np.asarray(T, float)
        x = T**(-1.0/3.0)
        xs, ls = self.x[::-1], self.lnk[::-1]          # x ascending
        lnk = np.interp(x, xs, ls)
        lo = x < xs[0]                                  # hotter than data
        hi = x > xs[-1]                                 # colder than data
        lnk = np.where(lo, ls[0] + self.b*(x - xs[0]), lnk)
        lnk = np.where(hi, ls[-1] + self.b*(x - xs[-1]), lnk)
        k = np.exp(lnk)
        if self.correction is not None:
            k = k/self.correction(T)
        return k


def tabulated_in_T(T_tab, v_tab):
    """A function of T, linear between the tabulated values, held at the end
    values outside them."""
    T_tab = np.asarray(T_tab, float)
    v_tab = np.asarray(v_tab, float)
    return lambda T: np.interp(np.asarray(T, float), T_tab, v_tab)


# ---- PUBLISHED RATES ----------------------------------------------------- #
#
# CH4 BENDING DYAD (nu4 1306 cm^-1 g = 3, nu2 1533 cm^-1 g = 2, held in
# Boltzmann ratio by fast intramode V-V exchange; the bins 1000-1750 cm^-1).
# The process in every rate below is dyad -> ground V-T(R), not intramode
# transfer.
#
# + He: Siddles, Wilson & Simpson (1994, Chem. Phys. 188, 99), Table 5,
#   CH4(nu2, nu4) + 4He, 140-295 K, the relaxation of the coupled pair,
#   laser-induced fluorescence: 1.7e-14 (295), 1.1e-14 (265), 9e-15 (240),
#   5.5e-15 (215), 3.6e-15 (190), 2.6e-15 (165), 1.7e-15 (140) cm^3 s^-1
#   (checked on the page image).  A relaxation rate, so divided by
#   1 + 3 exp(-1879 K/T) + 2 exp(-2206 K/T) (<= 0.6 per cent at 295 K).
#   Above 295 K: Landau-Teller continuation (fitted b = -54.6 K^(1/3)).
#   Cross-checks: Menard-Bourcin et al. (2005, J. Phys. Chem. A 109, 3111),
#   Table 4, elementary k_V-T(He) 1.69e-14 (296 K), 0.54e-14 (193 K, 1.4
#   times the Siddles value there); Yardley, Fertig & Moore (1970, J. Chem.
#   Phys. 52, 1450), Table I, p tau(CH4-4He) = 2.3 microsecond atm at room
#   temperature, i.e. 1.75e-14.
#
# + H2: Menard-Bourcin et al. (2005), Table 4, k_V-T(H2) = 1.44e-13 (193 K)
#   and 3.68e-13 (296 K) cm^3 s^-1, the elementary rate of CH4(nu4) + M ->
#   CH4(0) + M (their process 6, taken equal for nu2, process 7) in a kinetic
#   model with detailed balance, so no correction.  The only low-temperature
#   measurement of this pair.  Outside 193-296 K: Landau-Teller
#   continuation with the two-point slope b = -40.8 K^(1/3) (Halthore et al.
#   1994 used -40 from pure CH4).  Cross-check: Yardley et al. 1970 Table I,
#   p tau(CH4-H2) = 0.12 microsecond atm, i.e. 3.4e-13 at room temperature.
#
# A10: the dyad band coefficient with the rotational levels thermalized,
#   from the HITRAN line list (sum of g' A e^(-E'/kT) over the nu2 and nu4
#   upper levels over sum of g' e^(-E'/kT); main isotopologue, hitran.org
#   2026-10-03): 2.04, 1.96, 1.77, 1.55, 1.43 s^-1 at 150, 185, 296, 500,
#   700 K.  Published nu4 alone: 2.12 s^-1, Yelle (1991, ApJ 383, 380)
#   Table 2.
HC_OVER_K = 1.4387769                     # cm K


def _ch4_dyad_population_factor(T):
    T = np.asarray(T, float)
    return (1.0 + 3.0*np.exp(-HC_OVER_K*1306.0/T)
            + 2.0*np.exp(-HC_OVER_K*1533.0/T))


k_vt_ch4_dyad_he = LandauTellerRate(
    [295.0, 265.0, 240.0, 215.0, 190.0, 165.0, 140.0],
    [1.7e-14, 1.1e-14, 9e-15, 5.5e-15, 3.6e-15, 2.6e-15, 1.7e-15],
    relaxation_to_deexcitation=_ch4_dyad_population_factor)
k_vt_ch4_dyad_h2 = LandauTellerRate([193.0, 296.0], [1.44e-13, 3.68e-13])

NLTE_BANDS[('CH4', 'nu4+nu2 dyad')] = dict(
    A10=tabulated_in_T([150.0, 185.0, 296.0, 500.0, 700.0],
                       [2.04, 1.96, 1.77, 1.55, 1.43]),
    k_vt={'He': k_vt_ch4_dyad_he, 'H2': k_vt_ch4_dyad_h2},
    k_vt_estimated={},
    source='He: Siddles et al. 1994 Table 5; H2: Menard-Bourcin et al. '
           '2005 Table 4; A10: HITRAN lines')

# NH3 nu2 (umbrella mode, 932/968 cm^-1 inversion pair, the bins
# 600-1250 cm^-1).  Process: nu2 -> ground V-T,R (Hovis & Moore measure the
# decay of the nu2 fluorescence, both inversion components together).
#
# + He: Hovis & Moore (1980, J. Chem. Phys. 72, 2397), Table I: 1.1e-13
#   (198 K), 2.8e-13 (293 K, from Hovis & Moore 1978, J. Chem. Phys. 69,
#   4947, Table I), 4.1e-13 (398 K) cm^3 s^-1, +-10-20 per cent.  A
#   relaxation rate: divided by 1 + exp(-1367 K/T) (950 cm^-1, g1/g0 = 1;
#   1.03 at 398 K).  Below 198 K and above 398 K: Landau-Teller continuation
#   with the three-point slope (b = -37.5 K^(1/3)); the column spans 185-560
#   K, so the low end is 13 K and the high end 160 K beyond the data.
#
# + H2: NOT MEASURED (no measurement was found).  ESTIMATE, used by 'nlte'
#   and left out by 'nlte_partial': the He rate times the ratio
#   k(H2)/k(He) measured for a bending mode of an analog molecule, chosen
#   by `nh3_h2_rate_analog` (adapter option
#   --nh3-h2-rate-analog; `NH3_H2_RATE_ANALOGS`).  The argument common to
#   both: the H2 enhancement is V-R transfer into the widely spaced H2
#   rotational levels (Siddles et al. section 4.1), which depends on the
#   quantum to be removed and on H2, and less on the vibrator.  The ratio is
#   not a property of NH3, so the estimate is uncertain by a factor of a few;
#   the two analogs and 'nlte_partial' (H2 left out) bracket it.
#
#   'cd4_bending' (default): the V-T relaxation of the CD4 bending modes
#   (996/1054 cm^-1, the same quantum as NH3 nu2 within 5 per cent), Siddles
#   et al. (1994), Tables 1 and 2: 27.5 (295 K), 36 (240 K), 49 (190 K),
#   interpolated linearly in T, held at 27.5 above 295 K and at 49 below
#   190 K.  The CH4 dyad gives 22-27 (Menard-Bourcin et al. Table 4).  NH3
#   relaxes 16 times faster with He than CD4 does.
#
#   'water_bending_overtone': the H2O bending overtone 2nu2 (3151 cm^-1
#   level, relaxing by 1556-1595 cm^-1 steps), 295 K: k(H2) = 2.9e-12
#   cm^3 s^-1, Zittel & Masturzo (1991, J. Chem. Phys. 95, 8005), Table I,
#   H2(16)O; k(He) = 3.2e-13 cm^3 s^-1 (296 K), Hovis & Moore (1980),
#   Table III, H2(18)O(2nu2).  Ratio 9.06.  The two rates belong to
#   different isotopologues (the He one to H2(18)O, whose bending quantum is
#   0.4 per cent lower), and no He rate of H2O(nu2) itself was found
#   (Zittel & Masturzo 1989, J. Chem. Phys. 90, 977, Table III, gives nu2
#   by H2O and Ar only), so the ratio is taken at the 2nu2 level.  The H2
#   rate exists at 295 K only, so the ratio is HELD CONSTANT IN T over the
#   whole column (185-560 K); the He rate alone carries the temperature
#   dependence.  The quantum removed is 1.6-1.7 times that of NH3 nu2, and
#   water relaxes about 280 times faster with itself than with He (Hovis &
#   Moore 1980 Table III, 293-296 K), so this analog is a hydride like NH3
#   but a less close one in quantum than CD4.
# A10: HITRAN line list as for CH4: 13.46, 13.21, 12.56, 11.89, 11.37 s^-1
#   at 150, 185, 296, 500, 700 K.  Published: 15 s^-1, Weaver & Mumma
#   (1984, ApJ 276, 782) Table 4, from a single band intensity.
k_vt_nh3_nu2_he = LandauTellerRate(
    [198.0, 293.0, 398.0], [1.1e-13, 2.8e-13, 4.1e-13],
    relaxation_to_deexcitation=lambda T: 1.0 + np.exp(
        -HC_OVER_K*950.0/np.asarray(T, float)))
H2_OVER_HE_CD4_BENDING = tabulated_in_T(
    [190.0, 240.0, 295.0], [3.7e-13/7.5e-15, 5.1e-13/1.4e-14,
                            7.7e-13/2.8e-14])
H2_OVER_HE_WATER_BENDING_OVERTONE = tabulated_in_T([295.0],
                                                   [2.9e-12/3.2e-13])
NH3_H2_RATE_ANALOGS = {
    'cd4_bending': H2_OVER_HE_CD4_BENDING,
    'water_bending_overtone': H2_OVER_HE_WATER_BENDING_OVERTONE,
}
NH3_H2_RATE_ANALOG_DEFAULT = 'cd4_bending'


def k_vt_nh3_nu2_h2_estimated(T, analog=NH3_H2_RATE_ANALOG_DEFAULT):
    return k_vt_nh3_nu2_he(T)*NH3_H2_RATE_ANALOGS[analog](T)


def _nh3_nu2_h2_estimate_of(analog):
    return lambda T: k_vt_nh3_nu2_h2_estimated(T, analog)


# k_vt_estimated: {partner: {analog: k10(T)}}, the analog chosen by the
# column (`ColumnResponse(nh3_h2_rate_analog=...)`).
NLTE_BANDS[('NH3', 'nu2')] = dict(
    A10=tabulated_in_T([150.0, 185.0, 296.0, 500.0, 700.0],
                       [13.46, 13.21, 12.56, 11.89, 11.37]),
    k_vt={'He': k_vt_nh3_nu2_he},
    k_vt_estimated={'H2': {a: _nh3_nu2_h2_estimate_of(a)
                           for a in NH3_H2_RATE_ANALOGS}},
    source='He: Hovis & Moore 1980 Table I (1978 Table I at 293 K); H2: '
           'estimated, He x k(H2)/k(He) of an analog: CD4 bending (Siddles '
           'et al. 1994 Tables 1-2) or H2O 2nu2 (Zittel & Masturzo 1991 '
           'Table I over Hovis & Moore 1980 Table III); A10: HITRAN lines')

# Bands left in LTE (flagged): NH3 nu4+2nu2 (2.3 per cent of the flux on
# LHS 1140 b; no rates), CH4 nu3 pentad (0.2 per cent), H2O, C2H2, CO, CO2
# vibrational bands (each below 0.5 per cent).
COLLISION_PARTNERS = ('He', 'H2')

# The ambient field in a vibrational band at the climate state: the upward
# intensity of the isothermal stratosphere below, optically thick in the
# band core at T0, so J = B(T0)/2 and a = e_LTE(T0)/2.
ABSORBED_FRACTION = 0.5

# 'nlte': the two-level treatment of every vibrational band whose major
#   collision partners (He, H2) all have a measured or an estimated rate.
# 'nlte_partial': the same with the MEASURED rates only (an estimated
#   partner left out, so C10 and the band's emission are lower bounds).
# 'lte': every band in LTE, the upper bound on the cooling.
# 'rotational': rotational bands only, the lower bound.
TREATMENTS = ('nlte', 'nlte_partial', 'lte', 'rotational')

# How the band emission leaves the column (module docstring):
# 'thin': all of it is lost (the optically thin treatment, flagged past
#   the 10 per cent limit of `optical_depth_check`);
# 'cooling_to_space': the downward half is lost to the deep column, the
#   upward half escapes to space through the column above with the
#   k-distribution transmission (Dickinson 1972 Eq. 33 with the exchange
#   function toward the warmed layer and the gas above neglected).
TRANSFERS = ('thin', 'cooling_to_space')
# Validity of 'cooling_to_space': the first-order shift of T_match by the
# neglected exchange inside the domain may not exceed this [K].
EXCHANGE_SHIFT_MAX_K = 1.0
# Nodes whose rise is below this fraction of the rise at the match carry
# neither emission excess nor weight in T_match worth the exchange sum.
EXCHANGE_NODE_RISE_MIN = 1.0e-6
# Layers optically thinner than this emit in the exchange sum as sheets.
EXCHANGE_POINT_TAU = 1.0e-6

# Why rotational bands are in LTE.  The NH3 ground-state rotation-inversion
# lines carry A = 0.88 s^-1 intensity-weighted (largest 4.8 s^-1) and a
# pressure-broadening coefficient of 0.078 cm^-1 atm^-1 intensity-weighted
# (HITRAN line data, main isotopologue, below 400 cm^-1; air broadening, the
# He value being smaller by a factor of a few).  The collision rate 2 pi c
# gamma p is then 1.4e4 s^-1 at 1 microbar: even if only a tenth of it were
# rotationally inelastic and He broadened four times less, the de-excitation
# rate would stand 400 times above the largest A, and rises with depth.
ROTATIONAL_LTE_NOTE = ('NH3 rotation: A <= 4.8 s^-1, broadening rate '
                       '1.4e4 s^-1 at 1 microbar (HITRAN)')

THIN_ESCAPE_MIN = 0.9      # validity threshold of the thin treatment: upward
                           # escaping fraction of a bin's emission excess
                           # (flag below)
SHARE_THAT_MATTERS = 0.01   # bins carrying more than 1% of the excess


def band_partner_rates(entry, treatment,
                       analog=NH3_H2_RATE_ANALOG_DEFAULT):
    """The rate functions a treatment uses for a band: the measured ones,
    plus, under 'nlte', the estimated ones of the analog named."""
    rates = dict(entry['k_vt'])
    if treatment == 'nlte':
        for M, by_analog in entry.get('k_vt_estimated', {}).items():
            rates[M] = by_analog[analog]
    return rates


def vibrational_excitation_efficiency(entry, rates, T, n_partner,
                                      escape=None, d_escape=None):
    """eps = C10/(C10 + A10) and d eps/dT at fixed pressure.

    rates: {partner: k10(T)}; n_partner: {partner: number density at T}
    (proportional to 1/T at fixed pressure, which the derivative uses).
    escape, d_escape: the escape probability beta of the emitted photons and
    its T derivative; given, the radiative rate is beta A10 (Dickinson 1972,
    Eq. 33, the factor [1 + (1/2)(A/lambda) T(z,inf)]^-1) and
    eps_beta = C10/(C10 + beta A10)."""
    T = np.asarray(T, float)
    C10 = np.zeros_like(T)
    dC10 = np.zeros_like(T)
    h = 1e-4*T
    for M, kfun in rates.items():
        n = n_partner[M]
        k = kfun(T)
        dk = (kfun(T + h) - kfun(T - h))/(2.0*h)
        C10 += k*n
        dC10 += dk*n - k*n/T
    A = entry['A10'](T)
    dA = (entry['A10'](T + h) - entry['A10'](T - h))/(2.0*h)
    if escape is not None:
        dA = escape*dA + d_escape*A
        A = escape*A
    eps = C10/(C10 + A)
    deps = (A*dC10 - C10*dA)/(C10 + A)**2
    return eps, deps, C10


class ColumnResponse(object):
    """The column below p_match, built from a lower-atmosphere profile.

    p_bar, T0, r_cm, mbar_g: levels of the reference profile (any order);
    mix: {species: number fraction} on the same levels, including 'He',
    'H2', 'H' and the emitters; planet_mass_g; p_match_bar;
    nh3_h2_rate_analog: the analog of the estimated NH3 nu2 + H2 rate
    (`NH3_H2_RATE_ANALOGS`), used by treatment 'nlte' only.
    """

    def __init__(self, p_bar, T0, r_cm, mbar_g, mix, planet_mass_g,
                 p_match_bar, emission=None, dlnp=0.01,
                 nh3_h2_rate_analog=NH3_H2_RATE_ANALOG_DEFAULT):
        if nh3_h2_rate_analog not in NH3_H2_RATE_ANALOGS:
            raise ValueError('nh3_h2_rate_analog must be one of %s'
                             % (tuple(NH3_H2_RATE_ANALOGS),))
        self.nh3_h2_rate_analog = nh3_h2_rate_analog
        o = np.argsort(-np.asarray(p_bar, float))
        self.p_prof = np.asarray(p_bar, float)[o]
        self.T0_prof = np.asarray(T0, float)[o]
        self.r_prof = np.asarray(r_cm, float)[o]
        self.m_prof = np.asarray(mbar_g, float)[o]
        self.mix_prof = {k: np.asarray(v, float)[o] for k, v in mix.items()}
        self.Mp = float(planet_mass_g)
        self.p_match = float(p_match_bar)
        self.dlnp = float(dlnp)
        self.em = emission if emission is not None else LTEBandEmission()
        self.emitters = [sp for sp in self.em.tables if sp in self.mix_prof]
        # the column of each emitter above every pressure, N_X(p) =
        # int_0^p x_X dp'/(m g), from the profile (independent of T at fixed
        # pressure); the shallowest level adds x p/(m g) for the gas above it
        lp = np.log(self.p_prof*BAR)
        g = GNEWT*self.Mp/self.r_prof**2
        self.N_above_prof = {}
        for sp in self.emitters:
            f = self.mix_prof[sp]/(self.m_prof*g)
            N = np.zeros_like(lp)
            N[-1] = f[-1]*self.p_prof[-1]*BAR
            for i in range(lp.size - 2, -1, -1):
                pa, pb = self.p_prof[i+1]*BAR, self.p_prof[i]*BAR
                N[i] = N[i+1] + 0.5*(f[i] + f[i+1])*(pb - pa)
            self.N_above_prof[sp] = N
        # Multiplies every optical depth of the escape to space (the transfer
        # 'cooling_to_space' and the thin test); 1 is the physical column, 0
        # its optically thin limit (a test).
        self.transmission_scale = 1.0
        # DIAGNOSTIC: False removes the profile's gas above p_match from the
        # overlying column (tau_above_match is zero, so the escape of both
        # the thin test and the cooling-to-space transfer starts at p_match,
        # and the reported N_above is measured from p_match).  Above the
        # match the gas belongs to the wind, which carries C, N, O as atoms
        # and is warmer; the switch measures how much of the absorption
        # comes from molecules the wind model does not contain.
        self.include_gas_above_match = True
        self._tau_above = {}

    # ---- the domain ----------------------------------------------------- #
    def _profile_at(self, arr, p_bar, log=False):
        x = np.log(self.p_prof)[::-1]
        y = arr[::-1]
        if log:
            return np.exp(np.interp(np.log(p_bar), x, np.log(np.maximum(y, 1e-300))))
        return np.interp(np.log(p_bar), x, y)

    def build_domain(self, depth_H):
        """Nodes from p_match e^depth_H up to p_match, uniform in ln p."""
        n = int(np.ceil(depth_H/self.dlnp))
        s = -np.log(self.p_match*BAR) - np.linspace(n*self.dlnp, 0.0, n + 1)
        p_bar = np.exp(-s)/BAR
        d = dict(s=s, p_bar=p_bar)
        d['T0'] = self._profile_at(self.T0_prof, p_bar)
        d['r'] = self._profile_at(self.r_prof, p_bar)
        d['m'] = self._profile_at(self.m_prof, p_bar)
        d['mix'] = {k: self._profile_at(v, p_bar, log=(k not in ('He', 'H2')))
                    for k, v in self.mix_prof.items()}
        d['N_above'] = {k: self._profile_at(v, p_bar, log=True)
                        for k, v in self.N_above_prof.items()}
        if not self.include_gas_above_match:
            for k in d['N_above']:
                d['N_above'][k] = np.maximum(
                    d['N_above'][k] - d['N_above'][k][-1], 0.0)
        g = GNEWT*self.Mp/d['r']**2
        d['g'] = g
        width = control_interval_widths(s)
        geom = ColumnGeometry(s, d['r']**2, p_bar*BAR/g*width,
                              d['m']*g/KB)
        d['geom'] = geom
        return d

    # ---- the physics on the nodes --------------------------------------- #
    def band_treatment(self, sp, b, treatment):
        """'lte', 'nlte' or 'off' for band b of species sp."""
        kind = BANDS[sp][b][3]
        if treatment == 'rotational':
            return 'lte' if kind == 'rot' else 'off'
        if treatment == 'lte' or kind == 'rot':
            return 'lte'
        entry = NLTE_BANDS.get((sp, BANDS[sp][b][0]))
        if entry is None:
            return 'lte'
        covered = set(band_partner_rates(entry, treatment,
                                         self.nh3_h2_rate_analog))
        if treatment == 'nlte' and not set(COLLISION_PARTNERS) <= covered:
            return 'lte'
        return 'nlte'

    def band_net_cooling(self, d, sp, b, T, how, treatment='nlte',
                         transfer='thin', bin_emission=None):
        """Net cooling per molecule of one band at the node temperatures T,
        and its T derivative at fixed pressure, minus the same at T0.
        transfer 'cooling_to_space' needs bin_emission = (e_bin, de_bin),
        `bin_emission` of the species at T."""
        lp = np.log10(d['p_bar'])
        e, de = self.em.emission(sp, b, T, lp)
        e0, _ = self.em.emission(sp, b, d['T0'], lp)
        if how == 'off':
            z = np.zeros_like(T)
            return z, z
        escaping = transfer == 'cooling_to_space'
        if escaping:
            # e - sum_bins (1 - beta) e_bin: the emission that escapes
            esc = self.escape_factors(d)[sp]
            sel = esc['band'] == b
            omb = 1.0 - esc['beta'][sel]
            eb, deb = bin_emission
            e_c = e - np.sum(omb*eb[sel], axis=0)
            de_c = de - np.sum(omb*deb[sel], axis=0)
            e0_c = e0 - np.sum(omb*esc['e_bin0'][sel], axis=0)
        if how == 'lte':
            if escaping:
                return e_c - e0_c, de_c
            return e - e0, de
        entry = NLTE_BANDS[(sp, BANDS[sp][b][0])]
        rates = band_partner_rates(entry, treatment,
                                   self.nh3_h2_rate_analog)
        a = ABSORBED_FRACTION*e0
        pdyn = d['p_bar']*BAR

        def partners(TT):
            return {M: d['mix'][M]*pdyn/(KB*TT) for M in rates}
        if escaping:
            beta, dbeta = self._band_escape(e, de, e_c, de_c)
            beta0, _ = self._band_escape(e0, 0.0*e0, e0_c, 0.0*e0)
            eps, deps, _ = vibrational_excitation_efficiency(
                entry, rates, T, partners(T), escape=beta, d_escape=dbeta)
            eps0, _, _ = vibrational_excitation_efficiency(
                entry, rates, d['T0'], partners(d['T0']), escape=beta0,
                d_escape=0.0*beta0)
            return (eps*(e_c - a) - eps0*(e0_c - a),
                    deps*(e_c - a) + eps*de_c)
        eps, deps, _ = vibrational_excitation_efficiency(entry, rates, T,
                                                         partners(T))
        eps0, _, _ = vibrational_excitation_efficiency(entry, rates, d['T0'],
                                                       partners(d['T0']))
        return eps*(e - a) - eps0*(e0 - a), deps*(e - a) + eps*de

    @staticmethod
    def _band_escape(e, de, e_c, de_c):
        """beta = e_c/e of a band and its T derivative (1 where e = 0)."""
        ok = e > 0.0
        es = np.where(ok, e, 1.0)
        beta = np.where(ok, e_c/es, 1.0)
        dbeta = np.where(ok, (de_c*es - e_c*de)/es**2, 0.0)
        return beta, dbeta

    def band_epsilon(self, d, sp, b, T, treatment, bin_emission):
        """eps_beta of a two-level band under 'cooling_to_space' at T (1 for
        an LTE band): the share of an absorbed photon thermalized."""
        how = self.band_treatment(sp, b, treatment)
        if how != 'nlte':
            return np.ones_like(T)
        lp = np.log10(d['p_bar'])
        e, de = self.em.emission(sp, b, T, lp)
        esc = self.escape_factors(d)[sp]
        sel = esc['band'] == b
        omb = 1.0 - esc['beta'][sel]
        eb, deb = bin_emission
        e_c = e - np.sum(omb*eb[sel], axis=0)
        de_c = de - np.sum(omb*deb[sel], axis=0)
        beta, dbeta = self._band_escape(e, de, e_c, de_c)
        entry = NLTE_BANDS[(sp, BANDS[sp][b][0])]
        rates = band_partner_rates(entry, treatment, self.nh3_h2_rate_analog)
        pdyn = d['p_bar']*BAR
        n = {M: d['mix'][M]*pdyn/(KB*T) for M in rates}
        eps, _, _ = vibrational_excitation_efficiency(
            entry, rates, T, n, escape=beta, d_escape=dbeta)
        return eps

    def net_cooling_function(self, d, treatment, transfer='thin'):
        if transfer not in TRANSFERS:
            raise ValueError('transfer must be one of %s' % (TRANSFERS,))

        def D(T):
            tot = np.zeros_like(T)
            dtot = np.zeros_like(T)
            for sp in self.emitters:
                x_over_m = d['mix'][sp]/d['m']
                binem = (self.bin_emission(d, sp, T)
                         if transfer == 'cooling_to_space' else None)
                for b in range(len(BANDS[sp])):
                    how = self.band_treatment(sp, b, treatment)
                    if how == 'off':
                        continue
                    c, dc = self.band_net_cooling(d, sp, b, T, how, treatment,
                                                  transfer, binem)
                    tot += x_over_m*c
                    dtot += x_over_m*dc
            return tot, dtot
        return D

    # ---- transfer through the column ------------------------------------ #
    def bin_emission(self, d, sp, T, h=0.01):
        """LTE emission of every bin of species sp, 4 pi sigma_bar Bbar
        [erg s^-1 per molecule], and its T derivative (central difference
        over +-h K), at the nodes (bins, n)."""
        tab = self.em.tables[sp]
        lp = np.log10(d['p_bar'])

        def e_of(TT):
            return (4.0*np.pi*tab.sigma_bar(TT, lp)
                    * planck_bin_integral(tab.nu_lo, tab.nu_hi, TT))
        T = np.asarray(T, float)
        return e_of(T), (e_of(T + h) - e_of(T - h))/(2.0*h)

    def tau_above_match(self, sp):
        """Vertical optical depth (bins, g) of the profile's gas above
        p_match, int_0^p_match k_g(T0(p), p) x/(m g) dp, trapezoidal over
        the profile levels above the match, plus k x p/(m g) at the
        shallowest level for the gas above it (as N_above)."""
        if not self.include_gas_above_match:
            tab = self.em.tables[sp]
            return np.zeros((tab.nu_c.size, tab.wg.size))
        if sp in self._tau_above:
            return self._tau_above[sp]
        tab = self.em.tables[sp]
        above = self.p_prof < self.p_match
        lnp_m = np.log(self.p_match)
        x = np.log(self.p_prof)[::-1]

        def at_match(arr):
            return np.interp(lnp_m, x, arr[::-1])
        p = np.concatenate([[self.p_match], self.p_prof[above]])
        T0 = np.concatenate([[at_match(self.T0_prof)], self.T0_prof[above]])
        r = np.concatenate([[at_match(self.r_prof)], self.r_prof[above]])
        m = np.concatenate([[at_match(self.m_prof)], self.m_prof[above]])
        xs = np.concatenate([[np.exp(at_match(np.log(np.maximum(
            self.mix_prof[sp], 1e-300))))], self.mix_prof[sp][above]])
        g = GNEWT*self.Mp/r**2
        kg = tab.k_g(T0, np.log10(p))                        # (bins, n, g)
        f = kg*(xs/(m*g))[None, :, None]
        pd = p*BAR
        tau = f[:, -1, :]*pd[-1]
        tau = tau + np.sum(0.5*(f[:, 1:, :] + f[:, :-1, :])
                           * (pd[:-1] - pd[1:])[None, :, None], axis=1)
        self._tau_above[sp] = tau
        return tau

    def escape_factors(self, d):
        """The cooling-to-space escape probabilities of the domain d (held
        at the climate state), cached in d.  Per species: 'beta' (bins, n),
        the emission-weighted escape of every node's control layer;
        'tau_face' (bins, n+1, g), the optical depth from the top of the
        profile at the faces (face k is the bottom face of node k, face n
        the matching level); 'band' the band of each bin; 'e_bin0' the bin
        emission at T0; 'u' (bins, n, g), the escape to space alone."""
        if 'escape' in d:
            return d['escape']
        from scipy.special import expn
        out = {}
        lp = np.log10(d['p_bar'])
        pdyn = d['p_bar']*BAR
        ds = np.diff(d['s'])
        for sp in self.emitters:
            tab = self.em.tables[sp]
            kg0 = tab.k_g(d['T0'], lp)                         # (bins, n, g)
            # d tau / ds = k x p/(m g) (s = -ln p, tau growing downward)
            phi = kg0*(d['mix'][sp]*pdyn/(d['m']*d['g']))[None, :, None]
            dtau = 0.5*(phi[:, 1:, :] + phi[:, :-1, :])*ds[None, :, None]
            tau_top = self.tau_above_match(sp)
            tau_n = np.empty_like(phi)
            tau_n[:, -1, :] = tau_top
            tau_n[:, :-1, :] = tau_top[:, None, :] + np.cumsum(
                dtau[:, ::-1, :], axis=1)[:, ::-1, :]
            tau_n *= self.transmission_scale
            n = tau_n.shape[1]
            face = np.empty((tau_n.shape[0], n + 1, tau_n.shape[2]))
            face[:, 0, :] = tau_n[:, 0, :]
            face[:, n, :] = tau_n[:, -1, :]
            face[:, 1:n, :] = 0.5*(tau_n[:, 1:, :] + tau_n[:, :-1, :])
            tt, tb = face[:, 1:, :], face[:, :-1, :]
            dl = tb - tt
            small = dl < 1e-9
            dls = np.where(small, 1.0, dl)
            u = np.where(small, 0.5*expn(2, 0.5*(tt + tb)),
                         (expn(3, tt) - expn(3, tb))/(2.0*dls))
            wk = tab.wg[None, None, :]*kg0
            sw = np.sum(wk, axis=-1)
            beta = np.where(sw > 0.0,
                            0.5 + np.sum(wk*u, axis=-1)/np.where(sw > 0.0,
                                                                  sw, 1.0),
                            1.0)
            e0, _ = self.bin_emission(d, sp, d['T0'])
            out[sp] = dict(beta=beta, tau_face=face, band=tab.band_of_bins(),
                           e_bin0=e0, u=u)
        d['escape'] = out
        return out

    def exchange_heating(self, d, T, treatment, k0=None, bin_share_min=0.0):
        """The heating that the transfer 'cooling_to_space' leaves out: the
        exchange of the emission excess between the control layers of nodes
        k0..M (Dickinson 1972 Eq. 30 restricted to them), in the units of
        A_i w_i D_i.  Exact plane-parallel slab exchange for every g-point:
        the fraction of layer j's emission absorbed in layer i (i != j) is
        [E3|t_j-b_i| - E3|b_j-b_i| - E3|t_j-t_i| + E3|b_j-t_i|]/(2 dtau_j)
        (t, b the optical depths of a layer's top and bottom faces), the self
        absorption 1 - (1 - 2 E3(dtau))/(2 dtau).  The node gains the excess
        absorbed from all layers and loses the excess its own upward photons
        give to the layers above it in the domain (counted as returned by
        the cooling-to-space escape).  Two-level bands: excess weighted by
        eps_beta at the emitter and at the absorber.  Bins whose excess in
        k0..M is below bin_share_min of the total are skipped.  Returns
        (H over all nodes, zero below k0; record)."""
        from scipy.special import expn
        geom = d['geom']
        M = geom.M
        k0 = 1 if k0 is None else max(int(k0), 0)
        lp = np.log10(d['p_bar'])
        esc_all = self.escape_factors(d)
        Aw = geom.area*geom.w
        H = np.zeros(M + 1)
        blk = slice(k0, M + 1)
        jobs = []
        total = 0.0
        for sp in self.emitters:
            tab = self.em.tables[sp]
            esc = esc_all[sp]
            binem = self.bin_emission(d, sp, T)
            kg = tab.k_g(T, lp)
            kg0 = tab.k_g(d['T0'], lp)
            Bt = planck_bin_integral(tab.nu_lo, tab.nu_hi, T)
            B0 = planck_bin_integral(tab.nu_lo, tab.nu_hi, d['T0'])
            xm = d['mix'][sp]/d['m']
            eps = np.ones((len(BANDS[sp]), M + 1))
            on = np.zeros(len(BANDS[sp]), bool)
            for b in range(len(BANDS[sp])):
                on[b] = self.band_treatment(sp, b, treatment) != 'off'
                if on[b]:
                    eps[b] = self.band_epsilon(d, sp, b, T, treatment, binem)
            for ib in range(tab.nu_c.size):
                b = esc['band'][ib]
                if not on[b]:
                    continue
                # emission excess of each g-point [A w D units], eps-weighted
                dphi = (4.0*np.pi*tab.wg[None, :]*(kg[ib]*Bt[ib][:, None]
                                                    - kg0[ib]*B0[ib][:, None])
                        * (Aw*xm*eps[b])[:, None])           # (n, g)
                dphi[0] = 0.0
                tot_b = float(np.sum(np.abs(dphi[blk])))
                total += tot_b
                jobs.append((sp, ib, eps[b], dphi, tot_b))
        for sp, ib, eps_abs, dphi, tot_b in jobs:
            if tot_b <= bin_share_min*total or tot_b == 0.0:
                continue
            face = esc_all[sp]['tau_face'][ib]                  # (n+1, g)
            for ig in range(face.shape[1]):
                ph = face[k0:, ig]                              # faces k0..M+1
                dl = ph[:-1] - ph[1:]                           # layers k0..M
                if np.all(dl <= 0.0):
                    continue
                ex = dphi[k0:, ig]
                # A layer thinner than EXCHANGE_POINT_TAU emits as a sheet
                # at its mid depth (fluxes (1/2) E_2): the slab form divides
                # E_3 differences by dtau and loses them to round-off there.
                slab = dl >= EXCHANGE_POINT_TAU
                v = np.where(slab, ex/(2.0*np.where(slab, dl, 1.0)), 0.0)
                K = v.size
                w = np.zeros(K + 1)
                w[1:] += v
                w[:-1] -= v
                E = expn(3, np.abs(ph[:, None] - ph[None, :]))
                P = E.T @ w
                tm = 0.5*(ph[:-1] + ph[1:])
                sheet = (~slab) & (ex != 0.0)
                if sheet.any():
                    dd = ph[:, None] - tm[None, sheet]          # (faces, sheets)
                    sgn = np.where(dd < 0.0, 1.0, -1.0)         # +1 above
                    P += np.sum(sgn*0.5*ex[sheet][None, :]
                                * expn(2, np.abs(dd)), axis=1)
                absorbed = P[:-1] - P[1:] + ex
                # upward photons of layer i crossing the matching level
                top = ph[-1]
                cross = np.where(slab, (expn(3, ph[1:] - top)
                                        - expn(3, ph[:-1] - top))
                                 / (2.0*np.where(slab, dl, 1.0)),
                                 0.5*expn(2, tm - top))
                H[k0:] += eps_abs[k0:]*(absorbed - ex*(0.5 - cross))
        return H, dict(k0=k0, bins_total=len(jobs))

    def face_coefficient_function(self, d):
        xf = {k: 0.5*(d['mix'].get(k, 0.0*d['s'])[1:]
                      + d['mix'].get(k, 0.0*d['s'])[:-1])
              for k in ('He', 'H2', 'H')}
        mg = d['geom'].mg_face

        def c(Tf):
            kap = mixture_conductivity(Tf, xf['He'], xf['H2'], xf['H'])
            h = 1e-4*Tf
            dkap = (mixture_conductivity(Tf + h, xf['He'], xf['H2'], xf['H'])
                    - mixture_conductivity(Tf - h, xf['He'], xf['H2'], xf['H'])
                    )/(2*h)
            cc = kap*mg/Tf
            return cc, (dkap*mg/Tf - kap*mg/Tf**2)
        return c

    # ---- the solve ------------------------------------------------------- #
    def solve(self, F, treatment='nlte', depth_H=6.0, depth_step_H=3.0,
              p_floor_bar=None, dT_bottom_max=1e-3, leak_max=1e-6,
              verbose=False, transfer='thin',
              exchange_shift_max_K=EXCHANGE_SHIFT_MAX_K, heating=None):
        """Solve for the flux F [erg cm^-2 s^-1] into the column at p_match.

        The domain starts depth_H scale heights deep and is deepened by
        depth_step_H until the response at the node above the bottom is
        below dT_bottom_max [K] and the bottom flux below leak_max of F, or
        the bottom reaches p_floor_bar (the tropopause, where the climate
        isotherm ends; default the profile's deepest level).  transfer: one
        of TRANSFERS.  heating: a function (d, T) -> extra heating in the
        units of A w D (the full-exchange solve).  Returns a dict.
        """
        if treatment not in TREATMENTS:
            raise ValueError('treatment must be one of %s' % (TREATMENTS,))
        if transfer not in TRANSFERS:
            raise ValueError('transfer must be one of %s' % (TRANSFERS,))
        p_floor = (p_floor_bar if p_floor_bar is not None
                   else self.p_prof[0])
        depth_cap = np.log(p_floor/self.p_match)
        depth = min(depth_H, depth_cap)
        T_prev = None
        while True:
            d = self.build_domain(depth)
            geom = d['geom']
            D = self.net_cooling_function(d, treatment, transfer)
            if heating is not None:
                D = self._with_heating(d, D, heating)
            c = self.face_coefficient_function(d)
            if T_prev is None:
                guess = d['T0'].copy()
            else:
                guess = d['T0'] + np.interp(d['s'], T_prev[0],
                                            T_prev[1], left=0.0)
            if F == 0.0:
                T = d['T0'].copy()
                q = np.zeros(geom.M)
                nit = 0
            else:
                T, q, nit = self._solve_with_continuation(
                    geom, F, d['T0'][0], guess, c, D, verbose)
            dT = T - d['T0']
            leak = abs(geom.area_face[0]*q[0])/max(geom.area[-1]*abs(F), 1e-300)
            ok_depth = (abs(dT[1]) <= dT_bottom_max and leak <= leak_max)
            if (ok_depth or depth >= depth_cap - 1e-12 or F == 0.0
                    or depth_step_H <= 0.0):
                break
            T_prev = (d['s'], dT)
            depth = min(depth + depth_step_H, depth_cap)
        res = dict(domain=d, T=T, dT=dT, q=q, F=F, treatment=treatment,
                   transfer=transfer,
                   newton_iterations=nit, depth_H=depth,
                   depth_limited_by_floor=not ok_depth, leak_fraction=leak)
        excess, leak_abs, top = discrete_energy_budget(geom, T, q, F, D)
        res['budget'] = dict(excess=excess, bottom_flux=leak_abs, top=top,
                             relative_error=((excess + leak_abs - top)
                                             / max(abs(top), 1e-300)))
        res['T_match'] = float(T[-1])
        res['dT_match'] = float(dT[-1])
        if dT[-1] != 0.0:
            x = -np.log(d['p_bar']/self.p_match)          # <= 0 below
            f = dT/dT[-1]
            k = np.where(f >= 1.0/np.e)[0][0]
            if k > 0:
                xe = np.interp(1.0/np.e, [f[k-1], f[k]], [x[k-1], x[k]])
            else:
                xe = x[0]
            res['e_folding_H'] = float(-xe)
        else:
            res['e_folding_H'] = 0.0
        res['band_shares'] = self.band_shares(d, T, treatment, transfer)
        res['optical_depth'] = self.optical_depth_check(d, T, dT, res)
        if transfer == 'cooling_to_space':
            res['escape'] = self.escape_record(d, T, treatment)
            if heating is None:
                res['exchange'] = self.exchange_validity(
                    d, T, F, treatment, c, D, exchange_shift_max_K)
        return res

    @staticmethod
    def _with_heating(d, D, heating):
        """D minus a held heating per unit column mass."""
        geom = d['geom']
        Hn = heating(d)/(geom.area*geom.w)

        def D2(T):
            v, dv = D(T)
            return v - Hn, dv
        return D2

    def exchange_validity(self, d, T, F, treatment, c, D, shift_max):
        """The first-order shift of T_match by the exchange inside the
        domain that 'cooling_to_space' neglects, over the nodes warmed by
        more than EXCHANGE_NODE_RISE_MIN of the rise at the match; refused
        (ColumnTransferRefusal) above shift_max [K]."""
        dT = T - d['T0']
        if F == 0.0 or dT[-1] == 0.0:
            return dict(k0=d['geom'].M, shift_K=0.0, shift_max_K=shift_max)
        k0 = int(np.where(np.abs(dT) > EXCHANGE_NODE_RISE_MIN
                          * abs(dT[-1]))[0][0])
        k0 = max(k0, 1)
        H, rec = self.exchange_heating(d, T, treatment, k0=k0)
        dTr = linear_temperature_response(d['geom'], F, T, c, D, H)
        out = dict(k0=k0, p_k0_bar=float(d['p_bar'][k0]),
                   shift_K=float(dTr[-1]), shift_max_K=shift_max,
                   heating_over_F=float(np.sum(H[1:])
                                        / (d['geom'].area[-1]*F)))
        if abs(out['shift_K']) > shift_max:
            raise ColumnTransferRefusal(
                'the cooling-to-space treatment fails: the exchange of the '
                'emission excess inside the warmed layer (nodes from %.3e bar '
                'up) moves T_match by %.3f K to first order, more than %.3g K'
                % (out['p_k0_bar'], out['shift_K'], shift_max))
        return out

    def escape_record(self, d, T, treatment):
        """Per band: the emission-weighted escape beta at the match and its
        smallest value over the warmed nodes, and the escape to space."""
        esc_all = self.escape_factors(d)
        rows = []
        dT = T - d['T0']
        warm = np.abs(dT) > 1e-3*max(abs(dT[-1]), 1e-300)
        warm[0] = False
        for sp in self.emitters:
            esc = esc_all[sp]
            eb = esc['e_bin0']
            for b, bd in enumerate(BANDS[sp]):
                sel = esc['band'] == b
                if not sel.any():
                    continue
                w = eb[sel]
                sw = np.sum(w, axis=0)
                bb = np.where(sw > 0, np.sum(w*esc['beta'][sel], axis=0)
                              / np.where(sw > 0, sw, 1.0), 1.0)
                rows.append(dict(species=sp, band=bd[0],
                                 beta_match=float(bb[-1]),
                                 beta_min_warm=float(bb[warm].min()
                                                     if warm.any() else 1.0)))
        return rows

    def solve_with_exchange(self, F, treatment='lte', tol_K=1e-3,
                            max_iter=30, bin_share_min=1e-6, verbose=False,
                            **kw):
        """The balance with the full exchange of the emission excess between
        all the nodes of the domain added to 'cooling_to_space': fixed-point
        iteration on the exchange heating (`exchange_heating` with k0 = 1),
        each step a Newton solve with the heating held.  The domain is the
        one 'cooling_to_space' converges to.  Returns (result, history)."""
        base = self.solve(F, treatment=treatment, transfer='cooling_to_space',
                          exchange_shift_max_K=np.inf, **kw)
        depth = base['depth_H']
        T = base['T'].copy()
        hist = [base['T_match']]
        for it in range(max_iter):
            H, _ = self.exchange_heating(base['domain'], T, treatment, k0=1,
                                         bin_share_min=bin_share_min)
            res = self.solve(F, treatment=treatment,
                             transfer='cooling_to_space', depth_H=depth,
                             depth_step_H=0.0, heating=lambda dd, H=H: H,
                             **kw)
            T = res['T']
            hist.append(res['T_match'])
            if verbose:
                print('  exchange iteration %d: T_match %.4f K'
                      % (it, res['T_match']))
            if abs(hist[-1] - hist[-2]) < tol_K:
                res['exchange_heating'] = H
                return res, hist
        raise RuntimeError('the exchange iteration did not settle: %s' % hist)

    def solve_with_included_exchange(self, F, treatment='lte', **kw):
        """`cooling_to_space` with the exchange of the emission excess
        between all the warmed nodes of the domain INCLUDED (the transfer
        `cooling_to_space_exchange`): `solve_with_exchange` and the record the
        writers need.  Where `cooling_to_space` neglects that exchange and
        measures its first-order effect on T_match (refusing above
        `exchange_shift_max_K`), this solves it, so the answer has no
        validity bound of that kind.  res['exchange'] carries
        included = True, shift_K = T_match(with exchange) - T_match
        (cooling to space), the heating of the exchange over F and the
        pressure of the first node.  Validity: Dickinson's exchange function
        restricted to the domain, exact slab E_3 fractions for every g-point;
        for a two-level band the exchanged excess is weighted by eps_beta at
        both ends, exact for LTE bands and approximate for two-level bands
        (module docstring, VALIDITY)."""
        kw.setdefault('max_iter', 100)
        res, hist = self.solve_with_exchange(F, treatment=treatment, **kw)
        H = res['exchange_heating']
        geom = res['domain']['geom']
        res['exchange'] = dict(
            included=True, shift_K=float(hist[-1] - hist[0]),
            shift_max_K=float('inf'),
            heating_over_F=float(np.sum(H[1:])/(geom.area[-1]*F))
            if F != 0.0 else 0.0,
            p_k0_bar=float(res['domain']['p_bar'][1]),
            iterations=len(hist) - 1)
        return res

    def _solve_with_continuation(self, geom, F, Tb, guess, c, D, verbose):
        """Newton at F; on failure, continuation in F from zero."""
        try:
            return solve_conducted_response(geom, F, Tb, guess, c, D,
                                            verbose=verbose)
        except RuntimeError:
            pass
        T = guess.copy()
        for frac in np.linspace(0.1, 1.0, 10):
            T, q, nit = solve_conducted_response(geom, frac*F, Tb, T, c, D,
                                                 verbose=verbose)
        return T, q, nit

    # ---- diagnostics ----------------------------------------------------- #
    def band_shares(self, d, T, treatment, transfer='thin'):
        """Column-integrated emission excess of every band, A w x/m [e(T) -
        e(T0)] summed over the free nodes, and its share of the total."""
        geom = d['geom']
        rows = []
        total = 0.0
        for sp in self.emitters:
            x_over_m = d['mix'][sp]/d['m']
            binem = (self.bin_emission(d, sp, T)
                     if transfer == 'cooling_to_space' else None)
            for b, band in enumerate(BANDS[sp]):
                how = self.band_treatment(sp, b, treatment)
                if how == 'off':
                    val = 0.0
                else:
                    cb, _ = self.band_net_cooling(d, sp, b, T, how, treatment,
                                                  transfer, binem)
                    val = float(np.sum((geom.area*geom.w*x_over_m*cb)[1:]))
                rows.append([sp, band[0], how, val])
                total += val
        for r in rows:
            r.append(r[3]/total if total != 0.0 else 0.0)
        return rows

    def optical_depth_check(self, d, T, dT, res):
        """The optically thin test of the module docstring, bin by bin.

        For every bin, at every node: tau_bar, the emission-weighted mean
        vertical optical depth of the column above the node, and the largest
        single g-point value (reported).  The test is on the fraction of the
        bin's upward emission EXCESS (LTE) that the column above absorbs,

            sum_i ex_i 2 (1 - beta_i) / sum_i ex_i ,

        ex_i the bin's emission excess in the control layer of node i and
        beta_i that layer's cooling-to-space escape (`escape_factors`: one
        half down, the layer integral of E_2/2 over the g-points up), so that
        2 beta_i - 1 is the escaping share of its upward photons; past the
        validity limit
        when it exceeds 1 - THIN_ESCAPE_MIN (10 per cent) in a bin carrying
        more than SHARE_THAT_MATTERS of the column's excess, which sets
        res['thin_layer']['exceeded'] (no refusal).  Returns one record per
        band (largest values over its bins)."""
        geom = d['geom']
        lp = np.log10(d['p_bar'])
        out = []
        failures = []
        largest = (0.0, '', '', 0.0, 0.0)
        top = res['budget']['top'] if res['F'] else 0.0
        # the band's own excess under the treatment solved (band_shares);
        # the LTE excess of its bins only distributes it over the bins
        solved = {(r_[0], r_[1]): (r_[2], r_[3]) for r_ in res['band_shares']}
        esc = self.escape_factors(d)
        for sp in self.emitters:
            tab = self.em.tables[sp]
            kg = tab.k_g(T, lp)                       # (bins, n, g)
            N = d['N_above'][sp]
            tau_g = kg*N[None, :, None]
            wk = tab.wg[None, None, :]*kg
            tau_bar = np.sum(wk*tau_g, axis=-1)/np.sum(wk, axis=-1)
            tau_max = tau_g.max(axis=-1)
            band = tab.band_of_bins()
            sb = tab.sigma_bar(T, lp)
            sb0 = tab.sigma_bar(d['T0'], lp)
            Bt = planck_bin_integral(tab.nu_lo, tab.nu_hi, T)
            B0 = planck_bin_integral(tab.nu_lo, tab.nu_hi, d['T0'])
            ex = (4*np.pi*(sb*Bt - sb0*B0)*(d['mix'][sp]/d['m'])[None, :]
                  * (geom.area*geom.w)[None, :])
            ex[:, 0] = 0.0
            bin_excess = ex.sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                absorbed = (np.sum(ex*2.0*(1.0 - esc[sp]['beta']), axis=1)
                            / bin_excess)
            absorbed = np.where(bin_excess != 0.0, absorbed, 0.0)
            for b, bd in enumerate(BANDS[sp]):
                sel = band == b
                if not sel.any():
                    continue
                how, band_excess = solved[(sp, bd[0])]
                lte_sum = float(np.sum(bin_excess[sel]))
                scale = band_excess/lte_sum if lte_sum != 0.0 else 0.0
                share = band_excess/top if top else 0.0
                out.append(dict(treatment=how,
                    species=sp, band=bd[0], share_of_excess=share,
                    tau_bar_at_match=float(tau_bar[sel][:, -1].max()),
                    tau_max_g_at_match=float(tau_max[sel][:, -1].max()),
                    absorbed_fraction_of_excess=float(absorbed[sel].max())))
                if not top:
                    continue
                for k in np.where(np.abs(scale*bin_excess[sel]/top)
                                  > SHARE_THAT_MATTERS)[0]:
                    if absorbed[sel][k] > largest[0]:
                        largest = (float(absorbed[sel][k]), sp, bd[0],
                                   float(tab.nu_lo[sel][k]),
                                   float(tab.nu_hi[sel][k]))
                    if absorbed[sel][k] > 1.0 - THIN_ESCAPE_MIN:
                        failures.append(
                            '%s %s bin %.0f-%.0f cm^-1: %.1f per cent of its '
                            'emission excess absorbed above'
                            % (sp, bd[0], tab.nu_lo[sel][k], tab.nu_hi[sel][k],
                               100*absorbed[sel][k]))
        # Past the limit the thin solution is kept and flagged (user
        # decision 2026-10-04): the caller must carry the flag into every
        # record of the state.
        res['thin_layer'] = dict(
            exceeded=bool(failures), limit=1.0 - THIN_ESCAPE_MIN,
            max_absorbed=largest[0], species=largest[1], band=largest[2],
            bin_cm=(largest[3], largest[4]), failures=failures,
            message=('the optically thin treatment is past its validity '
                     'limit (more than %.0f per cent of the emission excess '
                     'of a bin carrying more than %.0f per cent of the column '
                     'total is absorbed by the column above): %s'
                     % (100*(1 - THIN_ESCAPE_MIN), 100*SHARE_THAT_MATTERS,
                        '; '.join(failures))) if failures else '')
        return out


class ColumnTransferRefusal(RuntimeError):
    """The exchange of the emission excess inside the warmed layer, which
    the cooling-to-space treatment neglects, moves T_match beyond its
    stated bound."""


# --------------------------------------------------------------------------
# Checks of the k-distribution band strengths
# --------------------------------------------------------------------------
# Integrated band intensities at 296 K [cm^-1/(molecule cm^-2)], sums of the
# HITRAN line intensities of the main isotopologue (natural abundance
# included, as HITRAN lists them), hitran.org line-by-line data retrieved
# 2026-10-03, lines selected by their global vibrational quanta; with the
# published values they agree with:
#   NH3 rotation-inversion (ground state, all below 400 cm^-1) 1.718e-17;
#   NH3 nu2 (both inversion components) 2.247e-17, Weaver & Mumma 1984
#       (ApJ 276, 782) Table 3: 2.2e-17;
#   CH4 nu4 5.141e-18 (dyad nu2+nu4 5.195e-18), Yelle 1991 (ApJ 383, 380)
#       Table 2: 5.42e-18.
HITRAN_BAND_INTENSITY_296K = {
    ('NH3', 'rotational'): 1.718e-17,
    ('NH3', 'nu2'): 2.247e-17,
    ('CH4', 'nu4+nu2 dyad'): 5.195e-18,
}
# Band-integrated Einstein coefficients with the rotational levels
# thermalized, from the same HITRAN lines, at 185 K, and the band origins
# and upper-level degeneracies of the two-level emission h c nu0 A10 g1
# exp(-h c nu0/k T) used to check the k-table emission of a fundamental.
TWO_LEVEL_CHECK = {
    ('NH3', 'nu2'): dict(nu0=950.0, A10=13.21, g1=1),
    ('CH4', 'nu4+nu2 dyad'): dict(nu0=1306.0, A10=2.16, g1=3),
}


def band_strength_report(emission=None, T=296.0, log10p=-6.0):
    """k-table band strength S = sum_{bins in band} sigma_bar dnu against
    the HITRAN sums, and the k-table LTE emission of the fundamentals
    against the two-level formula at 185, 300 and 500 K.  Returns lines."""
    em = emission if emission is not None else LTEBandEmission()
    lines = []
    for (sp, bname), S_hit in HITRAN_BAND_INTENSITY_296K.items():
        tab = em.tables[sp]
        b = [x[0] for x in BANDS[sp]].index(bname)
        sel = tab.band_of_bins() == b
        sb = tab.sigma_bar(np.array([T]), np.array([log10p]))[:, 0]
        S_k = float(np.sum(sb[sel]*tab.dnu[sel]))
        lines.append('%s %s: k-table S(%.0f K) = %.4e, HITRAN %.4e, ratio %.4f'
                     % (sp, bname, T, S_k, S_hit, S_k/S_hit))
    for (sp, bname), par in TWO_LEVEL_CHECK.items():
        b = [x[0] for x in BANDS[sp]].index(bname)
        for TT in (185.0, 300.0, 500.0):
            e_k, _ = em.emission(sp, b, np.array([TT]), np.array([log10p]))
            e_2 = (HP*CL*par['nu0']*par['A10']*par['g1']
                   * np.exp(-C2*par['nu0']/TT))
            lines.append('%s %s at %.0f K: k-table LTE emission %.4e erg/s,'
                         ' two-level h c nu0 A g1 exp(-hc nu0/kT) %.4e,'
                         ' ratio %.3f' % (sp, bname, TT, e_k[0], e_2,
                                          e_k[0]/e_2))
    return lines
