#!/usr/bin/env python3
"""Radiative-convective equilibrium column of the lower atmosphere.

Milestone E3 of `docs/phase_e_flux_closure_design.md` section 5: the step
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
                                      t_partition_max=1000.0, verbose=True):
    """Solve the column and return the climate solution.

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

    return dict(P_dyn=P, T=T, T_deep=float(T_deep),
                P_deep_bar=float(p_deep_bar),
                P_trop_bar=float(c.P_trop/BAR), T_trop=float(c.T_trop),
                f_H2O_deep=float(f_h2o[0]), f_H2O_trop=f_h2o_trop,
                mixing=mixing, bg_gas=bg, nlayers=int(nlayers),
                surface_albedo=float(surface_albedo))


def summary(sol):
    """One line per statement the profile header has to carry."""
    return ('climate solved: deep boundary %.4g bar at %.1f K, background '
            '%s, tropopause %.4g bar at %.1f K, cold-trap f_H2O %.3e '
            '(deep %.3e)'
            % (sol['P_deep_bar'], sol['T_deep'], sol['bg_gas'],
               sol['P_trop_bar'], sol['T_trop'], sol['f_H2O_trop'],
               sol['f_H2O_deep']))
