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

Deep composition -- the approximation here, with its validity range.
The climate solve needs a composition at the deep boundary before any
chemistry has been run.  It is the Gibbs equilibrium of the elemental vector
(H, He, C, N, O) at the deep pressure and the deep temperature, solved
together with that temperature (`DeepEquilibriumComposition`,
`deep_temperature_in_chemical_equilibrium`): the deep temperature of the
column depends on the opacity of its carriers, and the carriers on the
temperature.  No temperature alone fixes the partition: at fixed (p, T) the
CO/CH4 and N2/NH3 balances depend on the elemental vector, and deep
equilibrium puts nitrogen mainly in NH3, not in N2, in the deepest regions
of hydrogen atmospheres (Moses et al. 2011, ApJ 737, 15, section 3.1); at
20 bar on the solar vector NH3 holds 0.9993 of the nitrogen at 600 K and
0.66 at 900 K, and the CO share of the carbon at 1500 K is 0.45 on the
solar vector and 0.98 at He/H 5 (Photochem's equilibrium, MEASURED
2026-10-29).  The prescribed
partition used before 2026-10-29 (O in H2O, C in CH4, N in N2) and its
1000 K refusal were therefore not an equilibrium statement and are gone
(md/climate_deep_temperature_codex_20261029.md;
md/climate_equilibrium_deep_composition_20261029.md; update log stage 3
section 162).

THE APPROXIMATION: chemical equilibrium AT the deep boundary, i.e. no
quenching from deeper, hotter levels.  It is valid where the conversion time
of the slowest transition, N2 <-> NH3, is shorter than the mixing time at
p_deep; where it is not, the deep composition is the one quenched below the
boundary, which this column does not describe (Moses et al. 2011 for the
comparison of equilibrium and quenched profiles).  Elements other than H,
He, C, N and O (sulfur in the adapters' default vector, S/H = 1.3e-5) have
no carrier with climate opacity and are left out of the column.  The
equilibrium is mapped onto `CLIMATE_SPECIES`, the carriers clima has
opacities for, by the least relative change that keeps the El/H ratios of
the vector (`DeepEquilibriumComposition`), and the solve is refused when the
equilibrium puts more than
`DEEP_OMITTED_TOL` of the mixture outside them (other gases or condensates),
or when the mapped mixture departs from the elemental vector by more than
`DEEP_ELEMENT_RATIO_TOL`.  Nothing downstream inherits the deep partition:
the photochemistry that consumes T(p) computes its own composition from the
same elemental vector.
"""
import os

import numpy as np

BAR = 1.0e6                      # dyn/cm^2

# The species `clima` has gas opacities for in `photochem_clima_data`
# (k-distributions H2O/CH4/CO2/CO/NH3, CIA H2-H2/H2-He/H2-CH4/CO2-H2, and
# Rayleigh for the rest).  A carrier outside this list would be invisible to
# the radiative transfer, so the deep composition is mapped onto it and an
# equilibrium with a significant fraction elsewhere is refused.
CLIMATE_SPECIES = ('H2', 'He', 'H2O', 'CH4', 'CO2', 'CO', 'N2', 'NH3')
CLIMATE_CONDENSATES = ('H2O',)


# Element counts of the climate species: the element balance of the deep
# equilibrium mapped onto them, and their mean particle mass
# (lower_profile_schema.mean_molecular_weight).
CLIMATE_ELEMENT_COUNTS = {
    'H2': {'H': 2}, 'He': {'He': 1}, 'H2O': {'H': 2, 'O': 1},
    'CH4': {'C': 1, 'H': 4}, 'CO2': {'C': 1, 'O': 2}, 'CO': {'C': 1, 'O': 1},
    'N2': {'N': 2}, 'NH3': {'N': 1, 'H': 3}}

# The deep equilibrium (module docstring, "Deep composition").
DEEP_EQUILIBRIUM_ATOMS = ('H', 'He', 'C', 'N', 'O')
DEEP_EQUILIBRIUM_THERMO_FILE = 'deep_equilibrium_thermo.yaml'
# Largest mole fraction of the equilibrium mixture (gas and condensates)
# outside CLIMATE_SPECIES; above it the solve is refused.  What the column
# leaves out at 20 bar on the LHS 1140 b vectors (deep He/H 0.65, 1.2, 5;
# MEASURED 2026-10-29, md/climate_equilibrium_deep_composition_20261029.md)
# is below 3e-11 of the mixture at 900 K, and at 1100-1500 K 86-100 per
# cent atomic H, which carries no infrared opacity, then HCN (at most 12 per
# cent of it, at most 4e-2 of NH3) and CH3.  At 1e-6 of the mixture the
# omitted species move the mean particle mass and heat capacity by 1e-6,
# and their nuclei are returned to the carriers by the mapping.  The total
# reaches 1e-6 near 1450 K at deep He/H 0.65 and 1.2 and near 1500 K at 5;
# a condensate or a molecular carrier at that level is refused.  (The
# prototype's 1e-8 refuses deep He/H 0.65 on LHS 1140 b near its root,
# about 1200 K, for 3e-8 of atomic H.)
DEEP_OMITTED_TOL = 1.0e-6
# Largest relative departure of an El/H ratio of the mapped mixture from the
# elemental vector; above it the solve is refused.
DEEP_ELEMENT_RATIO_TOL = 1.0e-6
# The fixed point T_deep = clima(T_deep; composition(p_deep, T_deep)):
# tolerance of the bracketed root, the largest |clima(T) - T| accepted at
# the re-evaluation, and the bracket search.  The root tolerance is above
# the scatter of clima's own deep-temperature solve with its starting value
# (2e-5 K, MEASURED on LHS 1140 b at deep He/H 5, starts 500-750 K at fixed
# composition, 2026-10-29).
DEEP_ROOT_XTOL_K = 1.0e-4
DEEP_FIXED_POINT_TOL_K = 1.0e-3
DEEP_BRACKET_MAX_STEPS = 12
DEEP_BRACKET_FIRST_STEP_K = 20.0
DEEP_BRACKET_MIN_STEP_K = 1.0
# Starting values of clima's deep-temperature solve, as factors of the
# trial temperature, tried in turn after the trial temperature itself and
# the last root clima returned: clima's solve does not converge from a
# start far below its root (at deep He/H 5 on LHS 1140 b it fails from
# 400 K for a root of 918 K and converges from 500 K, MEASURED 2026-10-29).
# The start selects nothing: the root is clima's to its own tolerance.
CLIMA_START_FACTORS = (1.25, 1.5, 2.0, 0.8, 0.6)


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


class DeepCompositionRefused(SystemExit):
    """The deep equilibrium at a stated (p, T) cannot be represented by the
    climate carriers (or did not converge)."""


class DeepEquilibriumComposition:
    """Gibbs equilibrium of H/He/C/N/O at the deep boundary, mapped onto
    `CLIMATE_SPECIES`.

    Photochem's `ChemEquiAnalysis` on the thermodynamic set of
    `generate_zahnle_earth_thermo` restricted to DEEP_EQUILIBRIUM_ATOMS,
    condensates included (so that a condensing carrier is seen and refused,
    not dropped).  The file is written once, in `workdir`.  `abundances`
    maps an element name to its number ratio to hydrogen (`H` included as
    1); elements outside DEEP_EQUILIBRIUM_ATOMS are left out (module
    docstring).  THE APPROXIMATION -- equilibrium at the deep boundary, no
    quenching from below -- and its validity range are the module
    docstring's.

    THE MAPPING.  The equilibrium fractions x_k of the carriers are moved
    to the nearest mixture (least squares in the relative change x/x_k - 1)
    whose El/H ratios are those of the elemental vector:

        x = x_eq + W A^T (A W A^T)^+ (b - A x_eq),   W = diag(x_eq^2),

    A the element counts of the carriers (CLIMATE_ELEMENT_COUNTS), b the
    vector scaled to the hydrogen of x_eq.  It returns to the carriers the
    nuclei of the species the column cannot hold (at most DEEP_OMITTED_TOL
    of the mixture: atomic H, HCN, CH3, C2H6 near 1000 K) and the element
    imbalance of the Gibbs solve itself (its mass tolerance, 1e-6; MEASURED
    1e-7 to 4e-7 on the LHS 1140 b vectors at 900-1000 K).  Dropping those
    species instead moves a trace ratio by omitted/element total: 2.4e-11
    of HCN is 1.6e-6 of the nitrogen at deep He/H 5, 1000 K (MEASURED
    2026-10-29), above DEEP_ELEMENT_RATIO_TOL for no change of the climate.
    The largest relative change of a carrier is kept as
    `last['carrier_change']`."""

    def __init__(self, workdir, abundances):
        from photochem.equilibrate import (ChemEquiAnalysis,
                                           generate_zahnle_earth_thermo)
        self.thermo_file = os.path.join(workdir, DEEP_EQUILIBRIUM_THERMO_FILE)
        generate_zahnle_earth_thermo(self.thermo_file,
                                     atoms_names=list(DEEP_EQUILIBRIUM_ATOMS))
        self.solver = ChemEquiAnalysis(self.thermo_file)
        self.solver.use_prev_guess = False      # a function of (p, T) only
        h = float(abundances.get('H', 1.0))
        if h <= 0.0:
            raise SystemExit('the elemental vector carries no hydrogen')
        self.ratios = {el: float(abundances.get(el, 0.0))/h
                       for el in DEEP_EQUILIBRIUM_ATOMS}
        atoms = np.array([self.ratios[el] for el in self.solver.atoms_names])
        self.molfracs_atoms = atoms/atoms.sum()
        self.elements = [el for el in DEEP_EQUILIBRIUM_ATOMS
                         if self.ratios[el] > 0.0]
        self.counts = np.array([[CLIMATE_ELEMENT_COUNTS[sp].get(el, 0)
                                 for sp in CLIMATE_SPECIES]
                                for el in self.elements], dtype=float)
        self.last = None

    def composition(self, p_dyn, T):
        """Mole fractions of `CLIMATE_SPECIES` in equilibrium at (p_dyn
        [dyn/cm^2], T [K]), mapped as the class docstring states.  Raises
        DeepCompositionRefused (with the numbers) when the equilibrium does
        not converge, when more than DEEP_OMITTED_TOL of the mixture is
        outside CLIMATE_SPECIES, or when an El/H ratio of the mapped mixture
        departs from the vector by more than DEEP_ELEMENT_RATIO_TOL.
        `self.last` keeps the full equilibrium (`species`, every species
        with its mole fraction), the omitted fraction, the element-ratio
        departures before and after the mapping and the largest relative
        change of a carrier of the last call."""
        e = self.solver
        if not e.solve(float(p_dyn), float(T),
                       molfracs_atoms=self.molfracs_atoms):
            raise DeepCompositionRefused(
                'the deep chemical equilibrium did not converge at %.6g bar,'
                ' %.6g K' % (p_dyn/BAR, T))
        full = dict(zip(e.species_names,
                        np.asarray(e.molfracs_species, dtype=float)))
        omitted = float(sum(v for k, v in full.items()
                            if k not in CLIMATE_SPECIES))
        if omitted > DEEP_OMITTED_TOL:
            worst = max((k for k in full if k not in CLIMATE_SPECIES),
                        key=lambda k: full[k])
            raise DeepCompositionRefused(
                'the deep chemical equilibrium at %.6g bar, %.6g K puts '
                '%.3e of the mixture outside the climate carriers (largest '
                '%s %.3e), above the %.1e the column can leave out'
                % (p_dyn/BAR, T, omitted, worst, full[worst],
                   DEEP_OMITTED_TOL))
        x_eq = np.array([float(full.get(sp, 0.0)) for sp in CLIMATE_SPECIES])
        A = self.counts
        dropped = self.ratio_departure(x_eq)
        b = np.array([self.ratios[el] for el in self.elements]) \
            * (A[self.elements.index('H')] @ x_eq)
        W = x_eq**2
        lam = np.linalg.pinv((A*W) @ A.T) @ (b - A @ x_eq)
        x = x_eq + W*(A.T @ lam)
        if np.any(x < 0.0):
            raise DeepCompositionRefused(
                'the element-conserving mapping of the deep equilibrium at '
                '%.6g bar, %.6g K gives a negative carrier' % (p_dyn/BAR, T))
        change = float(np.max(np.abs(x[x_eq > 0.0]/x_eq[x_eq > 0.0] - 1.0)))
        x = x/x.sum()
        departure, el, r, r0 = self.ratio_departure(x, worst=True)
        if departure > DEEP_ELEMENT_RATIO_TOL:
            raise DeepCompositionRefused(
                'the climate carriers of the deep equilibrium at %.6g '
                'bar, %.6g K hold %s/H = %.9e against %.9e of the '
                'elemental vector (departure %.3e, tolerance %.1e)'
                % (p_dyn/BAR, T, el, r, r0, departure,
                   DEEP_ELEMENT_RATIO_TOL))
        self.last = dict(species=full, omitted=omitted,
                         departure_without_mapping=dropped,
                         element_ratio_departure=departure,
                         carrier_change=change,
                         p_dyn=float(p_dyn), T=float(T))
        return dict(zip(CLIMATE_SPECIES, (float(v) for v in x)))

    def ratio_departure(self, x, worst=False):
        """Largest relative departure of the El/H ratios of the carrier
        fractions `x` (CLIMATE_SPECIES order) from the vector."""
        nuclei = self.counts @ x
        h = nuclei[self.elements.index('H')]
        out = (0.0, 'H', 1.0, 1.0)
        for i, el in enumerate(self.elements):
            if el == 'H':
                continue
            r, r0 = nuclei[i]/h, self.ratios[el]
            d = abs(r/r0 - 1.0)
            if d > out[0]:
                out = (d, el, r, r0)
        return out if worst else out[0]


def deep_temperature_in_chemical_equilibrium(c, equilibrium, p_deep_dyn,
                                             t_deep_guess, t_trop_guess):
    """The deep temperature of the column `c` (AdiabatClimate, its settings
    already given) with the deep composition in chemical equilibrium at that
    temperature: the fixed point

        T_deep = clima(T_deep; equilibrium.composition(p_deep, T_deep)),

    clima(T; x) being the deep temperature `surface_temperature_bg_gas`
    solves (absorbed stellar = outgoing longwave flux) for the deep
    composition x held.  Root of R(T) = clima(T) - T by brentq on a bracket
    built outward from `t_deep_guess`: a first step of |R| (at least
    DEEP_BRACKET_FIRST_STEP_K) in the direction R points, doubled each time
    the sign holds, halved back toward the last accepted point when a trial
    temperature's equilibrium cannot be represented by the carriers
    (DeepCompositionRefused; the refusal stands if the step falls below
    DEEP_BRACKET_MIN_STEP_K), and no lower than the temperature floor of
    the opacity tables (STRATOSPHERE_T_MIN).  Each evaluation starts clima
    from the trial temperature and `t_trop_guess` (then, where clima's
    solve fails from there, from its last root and from CLIMA_START_FACTORS
    times the trial temperature), so R is a function of T alone to clima's
    tolerance.  clima is then re-evaluated at the converged composition, so
    `c` holds the solution of that composition, and |R| there must be below
    DEEP_FIXED_POINT_TOL_K.

    Returns (T_deep from the re-evaluation, mixing, evaluations of R,
    the re-evaluation's R, background gas)."""
    from scipy.optimize import brentq

    names = list(c.species_names)
    missing = [sp for sp in CLIMATE_SPECIES if sp not in names]
    if missing:
        raise SystemExit('the climate species file does not carry %s'
                         % ', '.join(missing))
    count = [0]
    held = {}

    def clima_deep_temperature(T):
        mixing = equilibrium.composition(p_deep_dyn, T)
        bg = background_gas(mixing)
        P_i = np.array([mixing.get(sp, 0.0)*p_deep_dyn for sp in names])
        count[0] += 1
        starts = [float(T)] + ([held['root']] if 'root' in held else []) \
            + [float(T)*f for f in CLIMA_START_FACTORS]
        failures = []
        for start in starts:
            c.T_trop = float(t_trop_guess)
            try:
                T_c = float(c.surface_temperature_bg_gas(
                    P_i, p_deep_dyn, bg, T_guess=start))
            except Exception as exc:              # noqa: BLE001
                failures.append('from %.1f K: %s' % (start, exc))
                continue
            held.update(mixing=mixing, bg=bg, root=T_c)
            return T_c
        raise SystemExit('clima found no deep temperature for the '
                         'equilibrium composition at %.6g K (%s)'
                         % (T, '; '.join(failures)))

    def residual(T):
        return clima_deep_temperature(T) - T

    a = float(t_deep_guess)
    ra = residual(a)
    b, rb = a, ra
    if ra != 0.0:
        direction = 1.0 if ra > 0.0 else -1.0
        step = max(DEEP_BRACKET_FIRST_STEP_K, abs(ra))
        expansions = 0
        while True:
            b = max(a + direction*step, STRATOSPHERE_T_MIN)
            try:
                rb = residual(b)
            except DeepCompositionRefused as exc:
                step *= 0.5
                if step < DEEP_BRACKET_MIN_STEP_K:
                    raise DeepCompositionRefused(
                        'the deep temperature lies beyond %.2f K (clima(T) '
                        '- T = %.4g K there), where the deep equilibrium '
                        'cannot be held by the climate carriers: %s'
                        % (a, ra, exc))
                continue
            if rb == 0.0 or (rb > 0.0) != (ra > 0.0):
                break
            if b == STRATOSPHERE_T_MIN:
                raise SystemExit(
                    'the deep temperature of the column is below the %.0f K '
                    'floor of the opacity tables (clima(T) - T = %.4g K '
                    'there)' % (STRATOSPHERE_T_MIN, rb))
            a, ra = b, rb
            step *= 2.0
            expansions += 1
            if expansions >= DEEP_BRACKET_MAX_STEPS:
                raise SystemExit(
                    'no bracket of the deep temperature in %d steps from '
                    '%.1f K (last trial %.1f K, clima(T) - T = %.4g K)'
                    % (DEEP_BRACKET_MAX_STEPS, float(t_deep_guess), b, rb))
    if ra == 0.0 or rb == 0.0:
        root = a if ra == 0.0 else b
    else:
        root = brentq(residual, min(a, b), max(a, b), xtol=DEEP_ROOT_XTOL_K)
    T_deep = clima_deep_temperature(root)
    fixed = T_deep - root
    if abs(fixed) > DEEP_FIXED_POINT_TOL_K:
        raise SystemExit(
            'the deep temperature is not a fixed point of the climate and '
            'the equilibrium composition: clima %.6f K at the composition '
            'of %.6f K (tolerance %.1e K)'
            % (T_deep, root, DEEP_FIXED_POINT_TOL_K))
    return T_deep, held['mixing'], count[0], fixed, held['bg']


def solve_radiative_convective_column(workdir, flux_file, mp_g, rp_cm,
                                      abundances, p_deep_bar,
                                      p_top_dyn=1.0e-2, surface_albedo=0.0,
                                      nlayers=60, t_deep_guess=400.0,
                                      t_trop_guess=120.0, verbose=True,
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

    The deep composition is the chemical equilibrium of `abundances` at
    (p_deep, T_deep), solved together with T_deep
    (`deep_temperature_in_chemical_equilibrium`; the approximation and its
    validity range are the module docstring's).  The equilibrium's
    thermodynamic file is written in `workdir`.

    Returns a dict: `P_dyn`, `T` (deep to shallow), `T_deep`, `P_trop_bar`,
    `T_trop`, `mixing` (the deep composition on the climate carriers),
    `deep_mixing` (the full deep equilibrium, every species of the Gibbs
    solve with its mole fraction), `deep_evaluations` (the evaluations of
    clima(T) - T), `deep_fixed_point_K` (clima(T) - T at the solution),
    `deep_omitted` (the equilibrium fraction outside the carriers),
    `deep_element_departure` (after the mapping), `deep_carrier_change`
    (the mapping's largest relative change of a carrier), `bg_gas`, and
    `f_H2O_trop`, the water mixing
    ratio the climate solution itself carries above the tropopause.
    """
    from photochem.clima import AdiabatClimate

    equilibrium = DeepEquilibriumComposition(workdir, abundances)
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
    p_deep_dyn = p_deep_bar*BAR
    T_deep, mixing, evaluations, fixed, bg = \
        deep_temperature_in_chemical_equilibrium(c, equilibrium, p_deep_dyn,
                                                 t_deep_guess, t_trop_guess)
    equil = equilibrium.last

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
                mixing=mixing, deep_mixing=equil['species'],
                deep_evaluations=int(evaluations),
                deep_fixed_point_K=float(fixed),
                deep_omitted=equil['omitted'],
                deep_element_departure=equil['element_ratio_departure'],
                deep_carrier_change=equil['carrier_change'],
                bg_gas=bg, nlayers=int(nlayers),
                surface_albedo=float(surface_albedo), stratosphere=strat)


def summary(sol):
    """One line per statement the profile header has to carry."""
    line = ('climate solved: deep boundary %.4g bar at %.1f K, background '
            '%s, tropopause %.4g bar at %.1f K, cold-trap f_H2O %.3e '
            '(deep %.3e)'
            % (sol['P_deep_bar'], sol['T_deep'], sol['bg_gas'],
               sol['P_trop_bar'], sol['T_trop'], sol['f_H2O_trop'],
               sol['f_H2O_deep']))
    if 'deep_evaluations' in sol:
        m = sol['mixing']
        line += ('; deep composition in chemical equilibrium at the deep '
                 'boundary (%d evaluations of clima(T) - T, %.1e K at the '
                 'solution; outside the carriers %.1e): %s'
                 % (sol['deep_evaluations'], sol['deep_fixed_point_K'],
                    sol['deep_omitted'],
                    ', '.join('f_%s %.4e' % (sp, m[sp])
                              for sp in CLIMATE_SPECIES)))
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
# Spacing in ln p of the nodes of the stratosphere solve; the column energy
# balance (conducted_heat_column_balance.ColumnResponse, dlnp) uses the same
# lattice anchored at the same level, so its nodes ARE nodes of this solve.
STRATOSPHERE_DLNP = 0.01
# max|R| over the largest term.  The residual left at the anchor is the
# reference imbalance a column built on the solution answers at zero input,
# where the gas is thinnest: on LHS 1140 b (deep He/H 2.09, equilibrium deep
# composition, K_zz 1e5) a stop at 1.8e-10 left -1.8e-3 K at 1e-6 bar
# (MEASURED 2026-10-29, update log stage 3 section 162), so the solve is
# taken one Newton step further (it fell to 7.2e-12 and the response to
# 1.7e-6 K there; the Newton residual falls 25-200 times per step).
STRATOSPHERE_TOL = 1.0e-11
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
