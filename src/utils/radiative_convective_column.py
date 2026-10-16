#!/usr/bin/env python3
"""Radiative-convective equilibrium column of the lower atmosphere.

Milestone E3: the step
that makes T(p) a solution rather than an input, and with it fixes the
tropopause the water cold trap sits at.  Photochem's chemistry API takes
`T(p)` and `K_zz(p)` as inputs, so the climate solution has to be computed
first and handed to it; that is what this module does, and
`photochem_to_lower_profile.py --climate` is its only caller.

The model is Photochem's `clima` (`AdiabatClimate`): a multispecies
pseudoadiabat (Graham et al. 2021) carried upward from a deep boundary and
connected to an isothermal stratosphere, with correlated-k gas opacities,
collision-induced absorption, Rayleigh scattering and the MT_CKD water
continuum, and a non-linear solve for the deep temperature that balances
the absorbed stellar flux against the outgoing longwave flux.
`solve_for_T_trop` is left on, so the stratospheric temperature is the skin
temperature of the solution and not a number chosen here.

With an eddy coefficient (`eddy_kzz`; photochem_to_lower_profile.py
--eddy-heat-transport on; off is the default since 2026-10-10) that isotherm is replaced by the
solution of the stratosphere's steady energy balance with clima's own
radiative transfer, molecular conduction and the eddy heat flux rho c_p K_zz
(dT/dz + g/c_p) of the K_zz that mixes the chemistry
(`stratosphere_energy_balance`), the troposphere and the deep boundary
staying clima's.  Where the eddy cooling of that K_zz exceeds what the
stratosphere's radiative heating can return, the solution leaves the
temperature range of the opacity tables and the solve is refused.

The tropopause of such a solution is the level where the adiabat meets that
isothermal stratosphere, i.e. the temperature minimum, and `clima` states it
directly as `P_trop`/`T_trop`.  It is the cold trap: the saturation vapour
pressure of water at `T_trop` divided by `P_trop` is the largest water
mixing ratio that can pass the level, whatever the deep abundance is.  This
module reports it; the chemistry solver that follows applies it, through the
mechanism's H2O condensate particle.

Deep composition -- the one approximation here, with its validity range.
The climate solve needs a composition at the deep boundary before any
chemistry has been run, so the elemental vector is partitioned into the
molecular carriers a cool hydrogen/helium atmosphere holds in chemical
equilibrium: oxygen in H2O, carbon in CH4, nitrogen in N2, the rest of the
hydrogen in H2, helium atomic.  That partition is the equilibrium one only
below the CO/CH4 and N2/NH3 transitions, roughly `T <~ 1000 K` at the
pressures of a deep boundary of order 1-100 bar (Lodders & Fegley 2002;
Moses et al. 2011).  The solved deep temperature is therefore checked
against `t_partition_max` and the solve is refused above it, rather than
returning a column whose opacity carriers are the wrong molecules.  Nothing
downstream inherits the partition: the photochemistry that consumes T(p)
computes its own composition from the same elemental vector.
"""
import os

import numpy as np

BAR = 1.0e6                      # dyn/cm^2

# The species `clima` has gas opacities for in `photochem_clima_data`
# (k-distributions H2O/CH4/CO2/CO/NH3, CIA H2-H2/H2-He/H2-CH4/CO2-H2, and
# Rayleigh for the rest).  A carrier outside this list would be invisible to
# the radiative transfer, so the partition below stays inside it.
CLIMATE_SPECIES = ('H2', 'He', 'H2O', 'CH4', 'CO2', 'CO', 'N2', 'NH3')
CLIMATE_CONDENSATES = ('H2O',)


def deep_molecular_partition(abundances):
    """Volume mixing ratios of `CLIMATE_SPECIES` from an El/H nuclei vector.

    `abundances` maps an element name to its number ratio to hydrogen (the
    vector the adapters carry, `H` included as 1).  Oxygen goes to H2O,
    carbon to CH4, nitrogen to N2, helium stays atomic and the hydrogen that
    is not bound in those carriers is H2 -- the cool-atmosphere equilibrium
    partition whose validity range the module docstring states.
    """
    xhe = abundances.get('He', 0.0)/abundances.get('H', 1.0)
    xo = abundances.get('O', 0.0)/abundances.get('H', 1.0)
    xc = abundances.get('C', 0.0)/abundances.get('H', 1.0)
    xn = abundances.get('N', 0.0)/abundances.get('H', 1.0)
    # per hydrogen nucleus; H2O takes 2 H, CH4 takes 4 H
    h2 = 0.5*max(0.0, 1.0 - 2.0*xo - 4.0*xc)
    per = {'H2': h2, 'He': xhe, 'H2O': xo, 'CH4': xc, 'N2': 0.5*xn,
           'CO': 0.0, 'CO2': 0.0, 'NH3': 0.0}
    tot = sum(per.values())
    if tot <= 0.0:
        raise SystemExit('the elemental vector carries no climate carrier')
    return {sp: v/tot for sp, v in per.items()}


def background_gas(mixing):
    """The most abundant species: the gas whose pressure `clima` solves for
    so that the stated deep pressure is met.  On a hydrogen atmosphere that
    is H2; on the helium-rich atmosphere the LHS 1140 b line implies it is
    He, and choosing it by abundance rather than by name is what lets one
    code path cover both."""
    return max(mixing, key=lambda sp: mixing[sp])


def write_climate_files(workdir, mp_g, rp_cm, surface_albedo, nlayers):
    """The species and settings files `AdiabatClimate` reads."""
    from photochem.utils import (species_file_for_climate,
                                 settings_file_for_climate)
    sp_file = os.path.join(workdir, 'climate_species.yaml')
    st_file = os.path.join(workdir, 'climate_settings.yaml')
    species_file_for_climate(sp_file, list(CLIMATE_SPECIES),
                             list(CLIMATE_CONDENSATES))
    settings_file_for_climate(st_file, float(mp_g), float(rp_cm),
                              float(surface_albedo),
                              number_of_layers=int(nlayers))
    return sp_file, st_file


def solve_radiative_convective_column(workdir, flux_file, mp_g, rp_cm,
                                      abundances, p_deep_bar,
                                      p_top_dyn=1.0e-2, surface_albedo=0.0,
                                      nlayers=60, t_deep_guess=400.0,
                                      t_trop_guess=120.0,
                                      t_partition_max=1000.0, verbose=True,
                                      zenith_u=None, photon_scale_factor=None,
                                      eddy_kzz=None, radius_cm=None,
                                      p_anchor_dyn=None):
    """Solve the column and return the climate solution.

    `eddy_kzz` (a function p [bar] -> K_zz [cm^2 s^-1], the coefficient the
    chemistry is solved with) replaces clima's isothermal stratosphere by the
    solution of its energy balance with radiation, molecular conduction and
    the eddy heat flux of that K_zz (`stratosphere_energy_balance`, with
    `radius_cm(p_bar, T, m_g)` the hydrostatic radius and `p_anchor_dyn` the
    level with no non-radiative flux through it, the matching level); None
    keeps clima's isotherm (the column before 2026-10-10).

    `zenith_u` (cosines, one per zenith angle of the settings file) and
    `photon_scale_factor` replace clima's own stellar geometry, a diurnal
    factor 0.5 (clima_radtran.f90:51) over a Gauss-Legendre quadrature of
    mu on [0, 1] (clima_eqns.f90:26-40), i.e. F/4 on the horizontal; None
    keeps it.  The adapter's face_geometry sets them so that the column
    takes the stellar energy per unit area the wind takes
    (md/matching_face_radiation_python_20261010.md section 2).

    THE TOP OF THE COLUMN IS A RADIATIVE BOUNDARY WITH NO DOWNWARD INFRARED
    FIELD (clima's boundary condition at P_top).  The infrared the escaping
    wind above the matching face sends down is not exchanged, and this is
    the measured reason (md/Q1_radiation_through_matching_face_20261009.md
    section 2, LHS 1140 b k05, MEASURED there): the wind's H3+ emission
    through the face, 1.757 erg cm^-2 s^-1 (99 per cent of its downward
    cooling radiation), is absorbed 7.2e-5 to 8.0e-5 of it between p_match
    and p_match e^0.5, 0.9-1.0 per cent above 1 mbar, 13-14 per cent at
    1 mbar-0.1 bar, 44 per cent at 0.1-1 bar and 41 per cent below 1 bar,
    by CH4, NH3, H2O and H2 collision-induced absorption: 1.76 against
    5.4e5 erg cm^-2 s^-1 of starlight at normal incidence (3e-6), the
    stellar near infrared the troposphere already absorbs.  Recorded only.

    `p_deep_bar` is the deep boundary the pseudoadiabat starts from and is
    the one free parameter of the solve: it is where the model stops looking
    downward, and both the deep temperature and the tropopause pressure
    depend on it.  It is reported with the solution rather than hidden in it.

    Returns a dict: `P_dyn`, `T` (deep to shallow), `T_deep`, `P_trop_bar`,
    `T_trop`, `mixing` (the deep partition), `bg_gas`, and `f_H2O_trop`, the
    water mixing ratio the climate solution itself carries above the
    tropopause.
    """
    from photochem.clima import AdiabatClimate

    mixing = deep_molecular_partition(abundances)
    bg = background_gas(mixing)
    sp_file, st_file = write_climate_files(workdir, mp_g, rp_cm,
                                           surface_albedo, nlayers)
    c = AdiabatClimate(sp_file, st_file, flux_file)
    if zenith_u is not None:
        zu = np.asarray(zenith_u, dtype=float)
        if zu.size != np.asarray(c.rad.zenith_u).size:
            raise SystemExit('the climate settings carry %d zenith angles,'
                             ' the stated geometry %d'
                             % (np.asarray(c.rad.zenith_u).size, zu.size))
        c.rad.zenith_u = zu
    if photon_scale_factor is not None:
        c.rad.photon_scale_factor = float(photon_scale_factor)
    c.solve_for_T_trop = True
    c.T_trop = float(t_trop_guess)
    c.P_top = float(p_top_dyn)
    c.verbose = bool(verbose)

    names = list(c.species_names)
    missing = [sp for sp in mixing if sp not in names]
    if missing:
        raise SystemExit('the climate species file does not carry %s'
                         % ', '.join(missing))
    p_deep_dyn = p_deep_bar*BAR
    P_i = np.array([mixing.get(sp, 0.0)*p_deep_dyn for sp in names])

    T_deep = c.surface_temperature_bg_gas(P_i, p_deep_dyn, bg,
                                          T_guess=float(t_deep_guess))
    if T_deep > t_partition_max:
        raise SystemExit(
            'the radiative-convective solution puts the deep boundary at '
            '%.0f K, above the %.0f K below which the CH4/H2O/N2 partition '
            'this solve was given is the equilibrium one: state a shallower '
            '--climate-p-deep, or give the climate step a composition from '
            'an equilibrium calculation' % (T_deep, t_partition_max))

    P = np.asarray(c.P, dtype=float)
    T = np.asarray(c.T, dtype=float)
    order = np.argsort(-P)                      # deep -> shallow
    P, T = P[order], T[order]
    f_i = np.asarray(c.f_i, dtype=float)
    iH2O = names.index('H2O')
    f_h2o = f_i[:, iH2O][order]
    # the cold trap: the water the solution carries above the tropopause
    above = P <= c.P_trop
    f_h2o_trop = float(f_h2o[above][0]) if above.any() else float(f_h2o[-1])

    strat = None
    if eddy_kzz is not None:
        strat = stratosphere_energy_balance(
            c, eddy_kzz, radius_cm, mp_g,
            p_anchor_dyn if p_anchor_dyn is not None else c.P_top,
            verbose=verbose)
        # clima's troposphere below the first node, the solved nodes, and
        # above the anchor the anchor's temperature (the continuation the
        # radiative transfer of the solve used)
        ps, Ts = strat['P_dyn'], strat['T']
        deep = P > ps[0]
        top = P < ps[-1]
        P = np.concatenate([P[deep], ps, P[top]])
        T = np.concatenate([T[deep], Ts, np.full(int(top.sum()), Ts[-1])])
        strat['T_min'] = float(np.min(Ts))
        strat['P_T_min_bar'] = float(ps[int(np.argmin(Ts))]/BAR)

    return dict(P_dyn=P, T=T, T_deep=float(T_deep),
                P_deep_bar=float(p_deep_bar),
                P_trop_bar=float(c.P_trop/BAR), T_trop=float(c.T_trop),
                f_H2O_deep=float(f_h2o[0]), f_H2O_trop=f_h2o_trop,
                mixing=mixing, bg_gas=bg, nlayers=int(nlayers),
                surface_albedo=float(surface_albedo), stratosphere=strat)


def summary(sol):
    """One line per statement the profile header has to carry."""
    line = ('climate solved: deep boundary %.4g bar at %.1f K, background '
            '%s, tropopause %.4g bar at %.1f K, cold-trap f_H2O %.3e '
            '(deep %.3e)'
            % (sol['P_deep_bar'], sol['T_deep'], sol['bg_gas'],
               sol['P_trop_bar'], sol['T_trop'], sol['f_H2O_trop'],
               sol['f_H2O_deep']))
    st = sol.get('stratosphere')
    if st is not None:
        line += ('; stratosphere in radiative + molecular + eddy energy '
                 'balance (K_zz %s), T %.1f K at %.3g bar (isotherm %.1f K), '
                 'eddy heat into the troposphere %.4g erg/cm2/s, TOA '
                 'OLR - ASR %.4g erg/cm2/s, residual %.1e'
                 % (st['kzz_note'], st['T'][-1], st['P_dyn'][-1]/BAR,
                    sol['T_trop'], st['bottom_flux'], st['toa_imbalance'],
                    st['residual']))
    return line


# --------------------------------------------------------------------------
# The stratosphere in energy balance with eddy (and molecular) conduction
# --------------------------------------------------------------------------
# Element counts of the climate species, for their mean particle mass
# (lower_profile_schema.mean_molecular_weight).
CLIMATE_ELEMENT_COUNTS = {
    'H2': {'H': 2}, 'He': {'He': 1}, 'H2O': {'H': 2, 'O': 1},
    'CH4': {'C': 1, 'H': 4}, 'CO2': {'C': 1, 'O': 2}, 'CO': {'C': 1, 'O': 1},
    'N2': {'N': 2}, 'NH3': {'N': 1, 'H': 3}}
# Spacing in ln p of the nodes of the stratosphere solve; the column energy
# balance (conducted_heat_column_balance.ColumnResponse, dlnp) uses the same
# lattice anchored at the same level, so its nodes ARE nodes of this solve.
STRATOSPHERE_DLNP = 0.01
STRATOSPHERE_TOL = 1.0e-9           # max|R| over the largest term
STRATOSPHERE_MAX_ITER = 60
STRATOSPHERE_FD_STEP_K = 0.05       # radiative Jacobian, forward difference
# The temperature range of the climate model's k-distribution tables and of
# the column's band emission (conducted_heat_column_balance.LTEBandEmission);
# a stratosphere outside it is refused, not clipped.
STRATOSPHERE_T_MIN = 50.0
# Relative tolerance of clima's altitude integration (dop853) during the
# solve.  At clima's own 1e-9 its fluxes are a step function of the
# temperatures: on LHS 1140 b k05 a change of 1e-6 K in every layer moved an
# edge flux by 4.4 erg cm^-2 s^-1 and one of 1e-8 K by 0.12 (MEASURED,
# md/eddy_heat_transport_column_20261010.md), which breaks the finite-
# difference Jacobian; at 1e-12 the response is linear (1.5e-3 and 1.4e-5).
STRATOSPHERE_CLIMA_RTOL = 1.0e-12


def radiative_edges(c):
    """Pressures [dyn/cm^2] of the edges of clima's radiative layers, surface
    first, up to P_top.  clima builds its grid as 2 nz + 1 points uniform in
    log10 P from P_surf to P_top (clima_adiabat_dry.f90, make_profile_dry:
    `linspace(log10(P_surf), log10(P_top), P)`), the layer centres being the
    even points; with the doubled radiative grid each layer is split into two
    radiative layers of equal thickness, so radiative edge j is grid point j
    (in an isotherm exactly, otherwise to the variation of the scale height
    across half a layer).  The two ghost layers above P_top are not
    returned."""
    nz = len(np.asarray(c.P))
    Pe = 10.0**np.linspace(np.log10(c.P_surf), np.log10(c.P_top), 2*nz + 1)
    Pe[0] = c.P_surf
    Pe[-1] = c.P_top
    return Pe


def net_flux_function(P_edges_dyn, F_edges):
    """F_up(p [bar]): the net upward radiative flux, linear in p between the
    radiative edges (a heating uniform per unit mass inside each radiative
    layer), held at the end values outside them."""
    x = np.asarray(P_edges_dyn, dtype=float)[::-1]/BAR
    y = np.asarray(F_edges, dtype=float)[::-1]
    return lambda p_bar: np.interp(np.asarray(p_bar, dtype=float), x, y)


def _clima_net_flux(c, T_centres, f_i):
    """Net upward radiative flux (infrared + stellar) at the radiative edges
    of `radiative_edges`, for clima's grid at the temperatures T_centres of
    its layer centres (surface held), by clima's own transfer
    (`TOA_fluxes_dry`, which rebuilds the column at fixed P and composition
    and integrates its altitude).  Returns (F_edges, ISR, OLR)."""
    nz = len(T_centres)
    P_in = np.concatenate([[c.P_surf], np.asarray(c.P, dtype=float)])
    T_in = np.concatenate([[c.T_surf], T_centres])
    ISR, OLR = c.TOA_fluxes_dry(P_in, T_in, f_i)
    fir = np.asarray(c.rad.wrk_ir.fup_n) - np.asarray(c.rad.wrk_ir.fdn_n)
    fsol = np.asarray(c.rad.wrk_sol.fup_n) - np.asarray(c.rad.wrk_sol.fdn_n)
    return (fir + fsol)[:2*nz + 1].copy(), float(ISR), float(OLR)


def stratosphere_energy_balance_of_held_troposphere(c, kzz, radius_cm, planet_mass_g,
                                p_anchor_dyn, dlnp=STRATOSPHERE_DLNP,
                                tol=STRATOSPHERE_TOL,
                                max_iter=STRATOSPHERE_MAX_ITER, verbose=False,
                                T_guess=None):
    """The stratosphere of the solved climate column `c` (AdiabatClimate) in
    steady energy balance with radiation, molecular conduction and, with
    `kzz`, eddy heat transport; the troposphere is clima's.

    WHY.  clima's stratosphere is an isotherm at the skin temperature, which
    is not a solution of its own radiative transfer (on LHS 1140 b k05, face
    geometry per_model, the net radiative heating of the column above the
    radiative edge at 0.962 bar is 1.78e4 erg cm^-2 s^-1, MEASURED,
    md/eddy_heat_transport_column_20261010.md), and
    AdiabatClimate has no term through which an eddy (or any non-radiative)
    heat flux enters (clima 0.9: its Python interface takes no heating
    profile; the eddy heat diffusion of clima's time-dependent `Climate`
    model is the mixing-length convective one and is not exposed).  The eddy
    coefficient that mixes the chemistry also carries heat (user decision
    2026-10-10), so the reference climate is solved here with it.

    EQUATION.  In s = -ln p, downward-positive fluxes, the discrete balance of
    conducted_heat_column_balance (`balance_terms`, total form) with the face
    flux q = c dT/ds + b of `conduction_face_functions` (molecular, plus
    eddy q_e = rho c_p K (dT/dz + g/c_p), Johnstone et al. 2018, A&A 617,
    A107, Eq. 44) and the radiative heating of each node's control interval
    H_i = A_i [F(s_{i-1/2}) - F(s_{i+1/2})], F the net upward radiative flux
    of clima's transfer at the current temperatures (`_clima_net_flux`,
    linear in p between its edges, `net_flux_function`):

        A_{i+1/2} q_{i+1/2} - A_{i-1/2} q_{i-1/2} + H_i = 0 ,  i = 1..M ,

    BOUNDARIES.  Bottom: the first node at or below clima's tropopause,
    Dirichlet at clima's temperature (the convective adjustment below stays
    clima's).  Top: the anchor level p_anchor_dyn (the matching level when
    the column is coupled to the wind, else clima's P_top) with NO
    non-radiative flux through it; above it the radiative column keeps the
    temperature of the anchor node (an isothermal continuation; that gas is
    the wind's in the coupled model and enters here only through its
    absorption).  Composition: clima's, held.  Geometry: r(p) from
    `radius_cm(p_bar, T, m_g)` of the current temperature (the column's own
    hydrostatic radius), so A, w and g are those of the column.

    METHOD.  Newton on the node temperatures with the conduction Jacobian
    analytic and the radiative one by forward differences of clima's
    transfer in the temperature of each layer centre inside the solve (one
    transfer per centre; the centres above the anchor move together with
    the anchor node), mapped to the nodes by the linear interpolation in s
    that sets the centre temperatures.  Converged when max|R| is below `tol`
    of the largest term (the radiative heating summed over the nodes, or the
    largest A q).

    Returns a dict: P_dyn, T (nodes, deep to the anchor), s, the radiative
    edges P_edges and net flux F_edges of the solution, the face fluxes
    (q_mol, q_eddy, p_face), the radiative heating H and the residual R of
    the nodes, residual (max|R|/scale), iterations, bottom_flux (the
    non-radiative flux into the troposphere, A-weighted per unit area at the
    bottom face), toa_imbalance (OLR - ASR of the new column; zero in clima's
    solution, the measure of holding the troposphere), lapse (dlnT/dlnp at
    the faces) and adiabat (R/c_p at the faces)."""
    import conducted_heat_column_balance as chb
    import lower_profile_schema as sch

    c.rtol = min(float(c.rtol), STRATOSPHERE_CLIMA_RTOL)
    names = list(c.species_names)
    P_c = np.asarray(c.P, dtype=float)
    T_c0 = np.asarray(c.T, dtype=float)
    f_c = np.asarray(c.f_i, dtype=float)
    f_in = np.maximum(np.vstack([np.asarray(c.f_i_surf, dtype=float)[None, :],
                                 f_c]), 1e-30)
    P_trop = float(c.P_trop)
    p_anchor = float(p_anchor_dyn)
    if not P_c[-1] <= p_anchor < P_trop:
        raise ValueError('the anchor level %.3e dyn/cm2 is not between the '
                         'tropopause and the top of the climate column'
                         % p_anchor)
    n = int(np.ceil(np.log(P_trop/p_anchor)/dlnp - 1e-9))
    s = -np.log(p_anchor) - dlnp*np.arange(n, -1, -1)       # bottom first
    p = np.exp(-s)
    # clima's profile with its surface, ascending ln p for interpolation
    lnP_all = np.log(np.concatenate([[c.P_surf], P_c]))[::-1]
    T_all = np.concatenate([[c.T_surf], T_c0])[::-1]
    T_bottom = float(np.interp(np.log(p[0]), lnP_all, T_all))
    # composition on the nodes (log-linear in p) and its particle mass
    lnPc = np.log(P_c)[::-1]
    mix = {}
    for k, sp in enumerate(names):
        mix[sp] = np.exp(np.interp(np.log(p), lnPc,
                                   np.log(np.maximum(f_c[::-1, k], 1e-300))))
    m_g = sch.mean_molecular_weight(mix, CLIMATE_ELEMENT_COUNTS)*sch.MAMU
    mix_face = {sp: 0.5*(v[1:] + v[:-1]) for sp, v in mix.items()}
    sf = 0.5*(s[1:] + s[:-1])
    K_face = (np.asarray(kzz(np.exp(-sf)/BAR), dtype=float)
              if kzz is not None else None)
    width = chb.control_interval_widths(s)
    Mp = float(planet_mass_g)

    # centre temperatures from the node temperatures: T_c = W T + T_held
    inside = (P_c <= p[0]) & (P_c >= p_anchor)
    above = P_c < p_anchor
    W = np.zeros((P_c.size, s.size))
    for j in np.where(inside)[0]:
        x = -np.log(P_c[j])
        k = int(np.clip(np.searchsorted(s, x) - 1, 0, s.size - 2))
        f = (x - s[k])/(s[k+1] - s[k])
        W[j, k] = 1.0 - f
        W[j, k+1] = f
    W[above, -1] = 1.0
    held = ~(inside | above)

    def centres(T):
        Tc = W @ T
        Tc[held] = T_c0[held]
        return Tc

    Pe = radiative_edges(c)
    # linear-in-p interpolation of the edge fluxes onto the node faces
    faces = np.concatenate([sf, [s[-1]]])
    pf = np.exp(-faces)
    xe = Pe[::-1]
    Pmap = np.zeros((faces.size, Pe.size))
    for i, pp in enumerate(pf):
        k = int(np.clip(np.searchsorted(xe, pp) - 1, 0, xe.size - 2))
        f = (pp - xe[k])/(xe[k+1] - xe[k])
        jlo, jhi = Pe.size - 1 - k, Pe.size - 2 - k        # back to Pe order
        Pmap[i, jlo] += 1.0 - f
        Pmap[i, jhi] += f

    deep_levels = P_c > p[0]

    def geometry(T):
        # deep to shallow: clima's levels below the first node, the nodes
        P_prof = np.concatenate([P_c[deep_levels], p])
        T_prof = np.concatenate([T_c0[deep_levels], T])
        m_prof = np.interp(np.log(P_prof), np.log(p)[::-1], m_g[::-1])
        r = radius_cm(P_prof/BAR, T_prof, m_prof)[int(deep_levels.sum()):]
        g = chb.GNEWT*Mp/r**2
        geom = chb.ColumnGeometry(s, r**2, p/g*width, m_g*g/chb.KB)
        cfun, bfun = chb.conduction_face_functions(
            geom, np.exp(-sf), 0.5*(g[1:] + g[:-1]), mix_face, K_face)
        return geom, cfun, bfun

    def zero_cooling(T):
        return np.zeros_like(T), np.zeros_like(T)

    def residual(T, geom, cfun, bfun):
        Fe, ISR, OLR = _clima_net_flux(c, centres(T), f_in)
        Ff = Pmap @ Fe
        H = np.zeros_like(T)
        H[1:] = geom.area[1:]*(Ff[:-1] - Ff[1:])
        t = chb.balance_terms(geom, 0.0, T, cfun, zero_cooling, None, bfun, H)
        return t, H, Fe, ISR, OLR

    if T_guess is None:
        T = np.interp(np.log(p), lnP_all, T_all)
    else:
        # (P_dyn deep to shallow, T): ascending in s = -ln p
        T = np.interp(s, -np.log(np.asarray(T_guess[0], dtype=float)),
                      np.asarray(T_guess[1], dtype=float))
    T[0] = T_bottom
    active = np.where(inside)[0]
    hist = []
    for it in range(max_iter):
        geom, cfun, bfun = geometry(T)
        t, H, Fe, ISR, OLR = residual(T, geom, cfun, bfun)
        scale = max(float(np.sum(np.abs(H[1:]))),
                    float(np.max(np.abs(geom.area_face*t['q']))), 1e-300)
        res = float(np.max(np.abs(t['R'][1:])))/scale
        hist.append(res)
        if verbose:
            print('    stratosphere newton %d: max|R|/scale %.3e, T(top) '
                  '%.2f K' % (it, res, T[-1]), flush=True)
        if res < tol:
            break
        # Jacobian: conduction (tridiagonal) + radiation (dense)
        _, ab, _ = chb.assemble_balance(geom, 0.0, T, cfun, zero_cooling,
                                        None, bfun, H)
        M = s.size - 1
        J = np.zeros((M, M))
        J[np.arange(M), np.arange(M)] = ab[1]
        J[np.arange(M - 1), np.arange(1, M)] = ab[0, 1:]
        J[np.arange(1, M), np.arange(M - 1)] = ab[2, :-1]
        Tc = centres(T)
        dFe_dTc = np.zeros((Pe.size, active.size))
        for jj, j in enumerate(active):
            Tp = Tc.copy()
            Tp[j] += STRATOSPHERE_FD_STEP_K
            dFe_dTc[:, jj] = (_clima_net_flux(c, Tp, f_in)[0] - Fe) \
                / STRATOSPHERE_FD_STEP_K
        dFe_dT = dFe_dTc @ W[active, :]
        if above.any():
            # the centres above the anchor follow the anchor node together
            Tp = Tc.copy()
            Tp[above] += STRATOSPHERE_FD_STEP_K
            dFe_dT[:, -1] += (_clima_net_flux(c, Tp, f_in)[0] - Fe) \
                / STRATOSPHERE_FD_STEP_K
        dFf_dT = Pmap @ dFe_dT
        dH = np.zeros((s.size, s.size))
        dH[1:] = geom.area[1:, None]*(dFf_dT[:-1] - dFf_dT[1:])
        J += dH[1:, 1:]
        dx = np.linalg.solve(J, -t['R'][1:])
        lim = np.max(np.abs(dx)/(0.3*T[1:]))
        if lim > 1.0:
            dx /= lim
        T[1:] += dx
        if np.any(~np.isfinite(T)) or np.min(T) < STRATOSPHERE_T_MIN:
            raise SystemExit(
                'the stratosphere energy balance (K_zz %s) leaves the '
                'temperature range of the opacity tables (%.0f K at %.3e '
                'bar)' % ('off' if kzz is None else 'on', np.min(T),
                          p[int(np.argmin(T))]/BAR))
    else:
        raise SystemExit('the stratosphere energy balance did not converge '
                         'in %d Newton iterations (max|R|/scale %.3e)'
                         % (max_iter, res))
    geom, cfun, bfun = geometry(T)
    t, H, Fe, ISR, OLR = residual(T, geom, cfun, bfun)
    mol, _ = chb.conduction_face_functions(geom, None, None, mix_face)
    Tf = 0.5*(T[1:] + T[:-1])
    grad = (T[1:] - T[:-1])/geom.ds
    cm, _ = mol(Tf)
    q_mol = cm*grad
    q_eddy = t['q'] - q_mol
    cp = chb.specific_heat(Tf, {k: v for k, v in mix_face.items()
                                if k in chb.JANAF_CP
                                or k in chb.CP_OVER_K_FIXED},
                           0.5*(m_g[1:] + m_g[:-1]))
    adiabat = chb.KB/(0.5*(m_g[1:] + m_g[:-1])*cp)
    lapse = -(np.log(T[1:]) - np.log(T[:-1]))/geom.ds
    # A superadiabatic face is convectively unstable: the convective region
    # would reach above clima's tropopause, which this solve holds.  On LHS
    # 1140 b k05 (per_model) that happens for K_zz <= 8e2 cm^2 s^-1 (lapse
    # 1.4 adiabat at 0.43 bar), where the eddy flux no longer carries the
    # heat the layer above the tropopause absorbs (MEASURED, the note).
    unstable = lapse/adiabat > 1.0 + 1e-9
    if np.any(unstable):
        k = int(np.argmax(lapse/adiabat))
        raise SystemExit(
            'the stratosphere energy balance (K_zz %s) is superadiabatic '
            '(lapse rate %.3f of the adiabat at %.3e bar): the convective '
            'region reaches above the tropopause clima placed, which this '
            'solve holds' % ('off' if kzz is None else 'on',
                             lapse[k]/adiabat[k], np.exp(-sf[k])/BAR))
    return dict(P_dyn=p, T=T, s=s, r=np.sqrt(geom.area), m=m_g,
                P_edges=Pe, F_edges=Fe, H=H, R=t['R'],
                q_mol=q_mol, q_eddy=q_eddy, p_face=np.exp(-sf),
                residual=hist[-1], history=hist, iterations=len(hist) - 1,
                bottom_flux=float(t['q'][0]*geom.area_face[0]
                                  / geom.area[0]),
                ISR=ISR, OLR=OLR, toa_imbalance=OLR - ISR,
                lapse=lapse, adiabat=adiabat,
                kzz_note=('off' if kzz is None else 'on'))


def stratosphere_energy_balance(c, kzz, radius_cm, planet_mass_g,
                                p_anchor_dyn, dlnp=STRATOSPHERE_DLNP,
                                tol=STRATOSPHERE_TOL,
                                max_iter=STRATOSPHERE_MAX_ITER, verbose=False,
                                T_guess=None):
    """stratosphere_energy_balance_of_held_troposphere, with clima's
    solution `c.T` restored on exit.  clima's `TOA_fluxes_dry` leaves the
    last profile it evaluated in `c.T` (MEASURED: a +30 K trial left c.T
    30 K off), so a second call on the same climate object held a drifted
    troposphere (md/eddy_dissipation_heating_review_20261010.md, defect
    found 2026-10-10); the climate object's solution is the caller's and is
    handed back unchanged.  `c.T` is read-only in clima's interface, so the
    climate object is put back the way clima itself writes it: one more
    `TOA_fluxes_dry` at the entry profile and composition."""
    T_entry = np.array(c.T, copy=True)
    Ts_entry = float(c.T_surf)
    P_entry = np.concatenate([[c.P_surf], np.asarray(c.P, dtype=float)])
    # the same surface-plus-layers composition the solve hands clima
    f_entry = np.maximum(np.vstack([np.asarray(c.f_i_surf, dtype=float)[None, :],
                                    np.asarray(c.f_i, dtype=float)]), 1e-30)
    try:
        return stratosphere_energy_balance_of_held_troposphere(
            c, kzz, radius_cm, planet_mass_g, p_anchor_dyn, dlnp=dlnp,
            tol=tol, max_iter=max_iter, verbose=verbose, T_guess=T_guess)
    finally:
        c.TOA_fluxes_dry(P_entry, np.concatenate([[Ts_entry], T_entry]),
                         f_entry)
