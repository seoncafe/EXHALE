#!/usr/bin/env python3
"""Photochem -> EXHALE lower-atmosphere profile (milestones E2 and E3).

Solves a gas-giant photochemistry model to steady state on a T(p), K_zz(p)
column -- prescribed (`--tp-file`, milestone E2) or solved for by the
radiative-convective climate model (`--climate`, milestone E3,
`radiative_convective_column.py`) -- and writes the handoff: one
`lower_atmosphere_profile.dat` and the minimal `base.inp` that carries the
same `solution_id`.  The EXHALE key that consumes it is

    Lower atmosphere profile: lower_atmosphere_profile.dat

The trial elemental fluxes ARE imposed on the chemistry: `--trial-flux-H`
and `--trial-flux-He` [g/s, outward positive] become a flux upper boundary
condition on the carriers at the photochemical model top, so that the
composition the file hands over is the composition of a column that is
losing what the wind takes.  A zero trial flux sets no boundary condition at
all, which is the closed-top solution E2 and E3 wrote.  What the solution
then delivers across its own top is measured with `gas_fluxes()` and written
beside the imposed value as `measured_flux_H` / `measured_flux_He`, so the
closure driver can check that the boundary condition took before it forms a
residual from it.  The `F_H` and `F_He` COLUMNS still carry the stated trial
values at every level; measuring the flux level by level over the overlap is
the wind's side of the closure (`element_flux_closure.py`, milestone E4).

With `--climate` the column is a solution and not an input: the climate
step fixes the tropopause, and with it the water cold trap, which the
chemistry solver then applies through the mechanism's H2O condensate
particle.  Two elemental sums are reported in that case, one over the gas
phase and one including the condensed carriers; their difference IS the cold
trap.

THE MODEL TOP IS THE MATCHING SURFACE (user decision 2026-10-07,
md/AN_overlap_structural_design_20261008.md).  Above p_match the gas is the
wind's: EXHALE alone gives T(p), r(p) and the hydrogen partition there, so
the chemistry column ends with its top FACE at p_match, where the escape
fluxes cross, and nothing of it is held isothermal above the match.
`steady_state_with_top_face_at_match` places the face there on the
solution's own grid (Photochem lays its levels out uniformly in altitude, so
one number fixes every level) and re-converges the chemistry on it,
iterating to the fixed point, and refuses if that will not settle.  The
temperature the chemistry runs on is the climate column with the rise the
column energy balance gives it under the conducted heat, up to the face
(`set_chemistry_temperature`, a persistent T(p) of Photochem).

THE HYDROGEN OF THE ESCAPE FLUX IS SPLIT AS THE WIND SPLITS IT.  With
`--trial-flux-H2` [g/s of hydrogen nuclei carried by H2, outward positive]
the H2 flux through the top face is that number and atomic H carries the
rest of the elemental hydrogen flux, negative where the wind mixes atoms
down into the column (`impose_elemental_escape_flux`).  The atoms that come
in recombine in the column, and their binding energy D0/2 per atom is
deposited in its energy balance with the shape of the column's own
H + H + M rate (`recombination_heating_of_inflowing_atoms`).

The driving pattern is the one `docs/p1_matched_comparison.py` validated on
2026-08-26 and is not re-derived here: dilution-only stellar flux, the Zahnle
set restricted to an atom list, the climate grid cut so the photochemical
grid can sit above it, elemental abundances mapped BY NAME, and the robust
stepper in blocks.

Environment: run this under a Python that has Photochem v0.9.3 built from
the source in `photochem/` -- upstream with the EXHALE corrections of
`README_photochem.md` section 4, which says what to download and how to
build it -- and its data photochem_clima_data v0.3.4.  Every profile this
writes records the Photochem version that made it in its
`# source_version` header line.

Example (HD 209458 b, the P1 configuration):

  python3 src/utils/photochem_to_lower_profile.py <run_dir> \\
      --mp 0.720 --r-ref 1.36 \\
      --tp-file vulcan_work/hd209_vulcan/atm/atm_HD209_Kzz.txt \\
      --stellar-flux vulcan_work/hd209_vulcan/atm/stellar_flux/Gueymard_solar.txt \\
      --r-star 1.155 --a-orb 0.0480 --p-match 1e-6
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lower_profile_schema as sch                       # noqa: E402
import radiative_convective_column as rcc                # noqa: E402

# Lodders 2009 as VULCAN reads it from
# fastchem_vulcan/input/solar_element_abundances.dat, i.e. the same elemental
# vector the P1 comparison gave both codes, so the two producers of E2
# differ in their chemistry and not in what they were told the gas is made of.
SOLAR_ABUNDANCES = dict(H=1.0, He=9.69170e-2, C=2.77588e-4, N=8.18465e-5,
                        O=6.06178e-4, S=1.31826e-5)


def parse_abundances(spec):
    """`--abundances He=0.0969,O=6.1e-4` on top of the Lodders vector."""
    out = dict(SOLAR_ABUNDANCES)
    if not spec:
        return out
    for item in spec.split(','):
        if not item.strip():
            continue
        name, _, val = item.partition('=')
        out[name.strip()] = float(val)
    return out


def read_tp_file(path, p_unit):
    """T(p) and, when the file carries a third column, K_zz(p).

    Accepts the two-or-three column tables both codes already keep (VULCAN's
    `atm/atm_*_Kzz.txt` among them): '#' comments and any non-numeric header
    lines are skipped.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                rows.append([float(w) for w in s.split()])
            except ValueError:
                continue                     # a column-name header line
    if len(rows) < 2:
        sch.refuse('the T(p) file %s carries fewer than two levels' % path)
    ncol = min(len(r) for r in rows)
    a = np.array([r[:ncol] for r in rows], dtype=float)
    P = a[:, 0]*(sch.BAR if p_unit == 'bar' else 1.0)      # -> dyn/cm^2
    T = a[:, 1]
    Kzz = a[:, 2] if ncol >= 3 else None
    o = np.argsort(-P)
    return P[o], T[o], (Kzz[o] if Kzz is not None else None)


def stellar_flux_at_planet(args):
    """(wavelength [nm], flux density [erg cm^-2 s^-1 nm^-1], dilution) of
    the stellar spectrum at the planet.

    A spectrum tabulated per nm in erg/cm^2/s/nm needs no numerical
    conversion (1 erg/cm^2/s/nm = 1 mW/m^2/nm, the unit Photochem and clima
    read); one tabulated per angstrom does, in both columns at once: the
    grid divides by 10 and the flux density multiplies by 10, so that the
    integral over wavelength is unchanged.  Only the dilution (R_star/a)^2
    is applied on top, and only when the spectrum is at the stellar surface.
    """
    wl, f = np.loadtxt(args.stellar_flux, unpack=True)
    if args.wavelength_unit == 'A':
        wl, f = wl/10.0, f*10.0
    if args.flux_at_planet:
        dilution = 1.0
    else:
        if args.r_star is None or args.a_orb is None:
            sch.refuse('--r-star and --a-orb are needed to dilute a stellar-'
                       'surface spectrum (or state --flux-at-planet)')
        dilution = (args.r_star*sch.RSUN/(args.a_orb*sch.AU))**2
    return wl, f*dilution, dilution


def write_flux_table(path, wl, flux):
    """The flux file Photochem and clima read: wavelength [nm], flux density
    [mW/m^2/nm]."""
    with open(path, 'w') as fh:
        fh.write('Wavelength (nm)               Solar flux (mW/m^2/nm)\n')
        for w, y in zip(wl, flux):
            fh.write('%.6e                  %.6e\n' % (w, y))


def write_flux_file(args, path):
    """Stellar flux at the planet, in mW/m^2/nm, on a wavelength grid in nm;
    Photochem and clima both read this one file."""
    wl, f, dilution = stellar_flux_at_planet(args)
    write_flux_table(path, wl, f)
    return dilution


# --------------------------------------------------------------------------
# RADIATION THROUGH THE MATCHING FACE
# (md/matching_face_radiation_python_20261010.md)
#
# One geometric convention for every stellar flux that crosses the face.
# The face fluxes the wind hands over (escape fluxes, conducted heat, the
# recombination heat of the inflowing atoms) are per unit face area of the
# wind's own normalization: the wind takes the stellar beam xi_w F at
# normal incidence on every shell, xi_w = dayside_dilution() of
# global_parameters (parameters.f90, function dayside_dilution; 1 for the
# '2D approximate method: Mdot' of the LHS 1140 b runs).  A column under
# that face describes the same horizontal average only if it receives the
# same stellar energy per unit area, xi_w F, with the angular distribution
# that average implies:
#   xi_w = 1    the substellar ray: all of the beam at mu = 1;
#   xi_w = 1/2  the dayside hemisphere: diurnal factor 1, mu uniform on
#               [0, 1] (the projected area of a hemisphere is uniform in
#               mu), mean mu 1/2;
#   xi_w = 1/4  the planet average: diurnal factor 1/2, mu as above.
# Per model (`per_model`, what each code assumes by itself; the default is
# `wind`, see --face-geometry):
#   Photochem: photon_flux * diurnal_fac at one zenith angle
#     (photochem_common.f90:311; diurnal_fac 0.5, photochem_types.f90:399;
#     the zenith is --zenith, 48 deg), i.e. 0.5 cos(48 deg) F = 0.335 F on
#     the horizontal;
#   clima: diurnal_fac 0.5 (clima_radtran.f90:51) over a Gauss-Legendre
#     quadrature of mu on [0, 1] with weights summing to 1
#     (clima_eqns.f90:26-40, number-of-zenith-angles 4 of the settings file),
#     i.e. 0.5 x 1/2 F = F/4.
# `wind` gives both the wind's xi_w.  Photochem carries one angle: it takes
# the insolation-weighted mu of the distribution, <mu^2>/<mu> = 2/3 for mu
# uniform on [0, 1] (48.19 deg), and the diurnal factor that makes its
# horizontal energy xi_w: xi_w/(2/3).
FACE_GEOMETRIES = ('per_model', 'wind')
PHOTOCHEM_DIURNAL_FACTOR = 0.5          # photochem_types.f90:399
CLIMATE_DIURNAL_FACTOR = 0.5            # clima_radtran.f90:51
HEMISPHERE_INSOLATION_WEIGHTED_MU = 2.0/3.0


def climate_quadrature(n):
    """clima's zenith quadrature (clima_eqns.f90 zenith_angles_and_weights):
    Gauss-Legendre nodes mapped onto mu in [0, 1], weights summing to 1."""
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5*x + 0.5, 0.5*w


def face_geometry(name, xi_wind=None, zenith_deg=48.0, n_climate_angles=4):
    """How each column model takes the stellar beam (the convention of the
    block above).

    Returns a dict: photochem_zenith_deg, photochem_diurnal_fac (None: leave
    Photochem's own), climate_zenith_u (None: clima's own quadrature),
    climate_photon_scale (None: 1), and the beams of each model: `*_mu`
    (direction cosines) and `*_beam` (fraction of the beam, so that the
    energy per unit horizontal area is sum(beam mu) F), with the totals
    `*_energy_factor`.
    """
    mu_q, w_q = climate_quadrature(n_climate_angles)
    if name == 'per_model':
        mu_p = np.cos(np.radians(zenith_deg))
        return dict(name=name, xi_wind=xi_wind,
                    photochem_zenith_deg=float(zenith_deg),
                    photochem_diurnal_fac=None,
                    photochem_mu=np.array([mu_p]),
                    photochem_beam=np.array([PHOTOCHEM_DIURNAL_FACTOR]),
                    photochem_energy_factor=PHOTOCHEM_DIURNAL_FACTOR*mu_p,
                    climate_zenith_u=None, climate_photon_scale=None,
                    climate_mu=mu_q,
                    climate_beam=CLIMATE_DIURNAL_FACTOR*w_q,
                    climate_energy_factor=float(
                        CLIMATE_DIURNAL_FACTOR*np.sum(w_q*mu_q)))
    if name != 'wind':
        raise ValueError('face geometry must be one of %s' % (FACE_GEOMETRIES,))
    if xi_wind is None:
        sch.refuse('--face-geometry wind needs the wind\'s dayside dilution:'
                   ' give --wind-dilution')
    xi = float(xi_wind)
    if abs(xi - 1.0) < 1e-9:
        # the substellar ray: Photochem at mu = 1 with no diurnal factor;
        # clima's quadrature collapsed onto mu = 1 and its fixed diurnal 0.5
        # undone by the photon scale
        return dict(name=name, xi_wind=xi, photochem_zenith_deg=0.0,
                    photochem_diurnal_fac=1.0,
                    photochem_mu=np.array([1.0]), photochem_beam=np.array([1.0]),
                    photochem_energy_factor=1.0,
                    climate_zenith_u=np.ones(n_climate_angles),
                    climate_photon_scale=1.0/CLIMATE_DIURNAL_FACTOR,
                    climate_mu=np.array([1.0]), climate_beam=np.array([1.0]),
                    climate_energy_factor=1.0)
    if abs(xi - 0.5) < 1e-9 or abs(xi - 0.25) < 1e-9:
        mu_p = HEMISPHERE_INSOLATION_WEIGHTED_MU
        d_p = xi/mu_p
        scale = xi/(CLIMATE_DIURNAL_FACTOR*float(np.sum(w_q*mu_q)))
        return dict(name=name, xi_wind=xi,
                    photochem_zenith_deg=float(np.degrees(np.arccos(mu_p))),
                    photochem_diurnal_fac=d_p,
                    photochem_mu=np.array([mu_p]), photochem_beam=np.array([d_p]),
                    photochem_energy_factor=xi,
                    climate_zenith_u=None, climate_photon_scale=scale,
                    climate_mu=mu_q,
                    climate_beam=CLIMATE_DIURNAL_FACTOR*scale*w_q,
                    climate_energy_factor=xi)
    sch.refuse('the wind\'s dayside dilution %.6g is none of 1, 1/2, 1/4, the'
               ' three conventions dayside_dilution() returns' % xi)


# The face geometry every run uses (face_geometry), set by main().
FACE_GEOMETRY = {}
# The eddy heat transport of the reference column and of the column energy
# balance (--eddy-heat-transport): 'kzz', the K_zz(p [bar]) the chemistry is
# solved with (lower_profile_schema.eddy_diffusion_coefficient, the one
# definition), and 'radiative_flux', the net upward radiative flux of the
# reference climate column it was solved with
# (radiative_convective_column.net_flux_function); empty when off.
EDDY_HEAT = {}
# number-of-zenith-angles of the climate settings file
# (photochem.utils.settings_file_for_climate default; the k05
# climate_settings.yaml states 4); radiative_convective_column refuses a
# climate object with another count when the convention changes its angles.
CLIMATE_ZENITH_ANGLES = 4


def write_face_radiation_record(args, geometry, xi_wind):
    """`wind_face_radiation.txt` beside the profile: the face convention by
    which the column models take the stellar beam (its solution_id carries
    the geometry, so the record and the profile name one solution)."""
    path = os.path.join(args.run_dir, 'wind_face_radiation.txt')
    g = geometry
    with open(path, 'w') as f:
        f.write('# radiation through the matching face\n')
        f.write('# face_geometry %s\n' % g['name'])
        f.write('# wind_dilution %s\n' % ('%.6g' % xi_wind
                                          if xi_wind is not None else 'none'))
        f.write('# photochem_zenith_deg %.6f diurnal_fac %s energy_factor'
                ' %.6f\n' % (g['photochem_zenith_deg'],
                             'own (0.5)' if g['photochem_diurnal_fac'] is None
                             else '%.6f' % g['photochem_diurnal_fac'],
                             g['photochem_energy_factor']))
        f.write('# climate_zenith_u %s photon_scale %s energy_factor %.6f\n'
                % ('own' if g['climate_zenith_u'] is None else
                   ' '.join('%.6f' % u for u in g['climate_zenith_u']),
                   'own (1)' if g['climate_photon_scale'] is None else
                   '%.6f' % g['climate_photon_scale'],
                   g['climate_energy_factor']))
    print('wrote %s' % path)


def mechanism_files(args, wdir):
    """(mechanism, thermo, data_dir) for the requested network."""
    if args.mechanism is not None:
        return args.mechanism, (args.thermo or args.mechanism), args.data_dir
    from photochem.utils import zahnle_rx_and_thermo_files
    atoms = [a.strip() for a in args.atoms.split(',') if a.strip()]
    tag = 'zahnle_' + ''.join(atoms)
    rx = os.path.join(wdir, tag + '_rxns.yaml')
    th = os.path.join(wdir, tag + '_thermo.yaml')
    if not os.path.exists(rx):
        zahnle_rx_and_thermo_files(atoms_names=atoms, rxns_filename=rx,
                                   thermo_filename=th,
                                   remove_reaction_particles=True)
    return rx, th, None


def model_top_radius(pc):
    """Radius [cm] of the face the upper boundary condition acts on.

    Photochem's own altitude grid states it: `var.top_atmos` is the top of
    the model domain measured from `dat.planet_radius`, i.e. the upper face
    of the topmost cell, which is exactly the surface the flux crosses.  It
    is used rather than the profile's hydrostatic `r` column because the
    boundary condition is applied on photochem's grid, and the two
    integrations need not agree to the last digit.
    """
    return float(pc.dat.planet_radius + pc.var.top_atmos)


def hydrogen_carrier_split(pc, counts):
    """Fraction of the hydrogen nuclei at the model top carried by atomic H.

    The elemental hydrogen flux has to be imposed on carriers, and at the
    top of a cool H2 atmosphere essentially all of it is H2.  The fraction
    is MEASURED here rather than assumed, so that a column whose top has
    dissociated splits the flux instead of putting it all on H2.
    """
    names = list(pc.dat.species_names)[:int(pc.dat.nq)]
    top = np.asarray(pc.wrk.usol, dtype=float)[:, -1]
    nuc_H = nuc_tot = 0.0
    for j, sp in enumerate(names):
        k = counts.get(sp, {}).get('H', 0)
        if k <= 0:
            continue
        nuc_tot += k*top[j]
        if sp == 'H':
            nuc_H += k*top[j]
    if nuc_tot <= 0.0:
        sch.refuse('the model top carries no hydrogen nuclei: an elemental '
                   'hydrogen flux cannot be imposed on it')
    return nuc_H/nuc_tot


# Atomic-H share of the hydrogen nuclei at the model top below which the
# whole elemental hydrogen flux is imposed on H2 (the dominant-carrier
# approximation), and the largest difference allowed between the share the
# boundary condition carries and the share of the converged top: either way
# at most 1% of the hydrogen flux leaves on the wrong carrier.
CARRIER_SHARE_TOL = 0.01
# Re-impositions of the measured split allowed before the adapter refuses.
CARRIER_ROUNDS = 8
# The most recent imposition (set by impose_elemental_escape_flux): the
# carrier decision ('split', 'H_atomic_share_used'), the stated mass fluxes
# and atomic share it was made with ('phi_H_gs', 'phi_He_gs', 'phi_H2_gs',
# 'f_atomic'), and 'r_top', the top-face radius the flux densities were
# computed with (None when no flux was imposed).
LAST_IMPOSITION = {}
# The atomic share a carrier round imposes (key 'H'): the share of the
# CONVERGED top.  None measures it on the composition the imposition sees;
# a re-solve on a new temperature first re-initializes the composition from
# the climate grid, whose top no longer carries the converged share.
CARRIER_SHARE_STATED = {'H': None}


def impose_elemental_escape_flux(pc, counts, phi_H_gs, phi_He_gs,
                                 phi_H2_gs=None, atomic_share=None):
    """Impose the elemental escape fluxes as a flux upper boundary condition.

    Phi_El [g/s, outward positive] is the mass of El nuclei leaving the top
    of the lower-atmosphere column per unit time.  It becomes a nuclei flux
    density at the model top,

        phi_El = Phi_El / m_El / (4 pi r_top^2)   [nuclei/cm^2/s],

    and is then imposed on the carriers of those nuclei.

    SIGN, determined from photochem's own right-hand side rather than
    guessed (`src/evoatmosphere/photochem_evoatmosphere_rhs.f90`, the upper
    boundary block): with `FluxBC` the top cell gets
    `rhs -= var%upper_flux/dz`, so a POSITIVE flux is a loss from the model
    top.  Outward positive is therefore photochem's convention as well, and
    no sign is flipped here.  `gas_fluxes()` returns the top flux with the
    same sign: in steady state its `top_fluxes` equals `upper_flux` for a
    species held at `FluxBC`.

    HYDROGEN CARRIERS STATED BY THE WIND.  With phi_H2_gs [g/s of the
    hydrogen nuclei H2 carries, outward positive] the H2 flux through the
    top face is phi_H2/(2 m_H)/(4 pi r_top^2) molecules and atomic H carries
    the rest of the elemental flux, (phi_H - phi_H2)/m_H/(4 pi r_top^2),
    which is negative (an inflow of atoms) where the wind's eddies mix the
    atomic hydrogen it makes aloft down through the matching surface.  These
    are the wind's own face fluxes (EXHALE_resolved.out
    base_face_H2_flux_cgs and the elemental flux), so the two models
    exchange the species fluxes of H and H2 and not only their sum.

    HYDROGEN CARRIER WITHOUT THE WIND'S SPLIT.  Without phi_H2_gs (no wind
    solution states one yet) the whole elemental hydrogen flux is put on
    H2, at half the nuclei flux because H2 carries two nuclei.  This is the
    dominant-carrier approximation, and it is applied only when the measured
    atomic-H share of the hydrogen nuclei at the model top is below 1%; the
    share is printed with the imposition.  Above 1% the flux is split
    between H and H2 in proportion to the nuclei each carries there, so the
    approximation never silently covers a dissociated top.  The share is
    that of the composition the imposition sees (the initial one on the
    first solve) unless a carrier round states it (CARRIER_SHARE_STATED);
    main() compares the share the boundary condition carries with that of
    the converged top and, where they differ by more than
    CARRIER_SHARE_TOL, re-imposes the measured split and re-solves.

    A zero flux sets no boundary condition at all, so the closed-top
    solution is reproduced bit for bit.

    THE AREA IS THAT OF THE TOP FACE AS IT STANDS NOW.  The stated fluxes
    are mass loss rates [g/s] through the face, and the flux densities are
    those rates over 4 pi r_top^2 of the current grid.  Any later change of
    the grid moves r_top, and the same densities then carry a different
    mass flux; steady_state_with_top_face_at_match therefore re-imposes the
    condition (impose_escape_flux_at_current_top) whenever r_top differs
    from LAST_IMPOSITION['r_top'], and accepts a state only when the two
    are equal.  `atomic_share` is the atomic-H share of the column-split
    branch to use instead of measuring it, so that a re-imposition keeps the
    split it re-imposes.
    """
    r_top = model_top_radius(pc)
    area = 4.0*np.pi*r_top*r_top
    imposed = {}
    if phi_H2_gs is not None:
        phi_H = phi_H_gs/(sch.ATOMIC_WEIGHT['H']*sch.MAMU)/area
        phi_nuc_H2 = phi_H2_gs/(sch.ATOMIC_WEIGHT['H']*sch.MAMU)/area
        pc.set_upper_bc('H2', bc_type='flux', flux=phi_nuc_H2/2.0)
        pc.set_upper_bc('H', bc_type='flux', flux=phi_H - phi_nuc_H2)
        imposed['H'] = phi_H
        imposed['H2_molecules'] = phi_nuc_H2/2.0
        imposed['H_atoms'] = phi_H - phi_nuc_H2
        LAST_IMPOSITION['split'] = 'wind'
        print('  imposed elemental H escape: %.6e g/s -> %.6e nuclei/cm^2/s '
              'at r_top = %.6e cm, split as the wind states it: H2 %.6e '
              'molecules/cm^2/s, atomic H %.6e atoms/cm^2/s'
              % (phi_H_gs, phi_H, r_top, phi_nuc_H2/2.0,
                 phi_H - phi_nuc_H2))
    elif phi_H_gs != 0.0:
        LAST_IMPOSITION['split'] = 'column'
        phi_H = phi_H_gs/(sch.ATOMIC_WEIGHT['H']*sch.MAMU)/area
        if atomic_share is not None:
            f_atomic = atomic_share
        elif CARRIER_SHARE_STATED.get('H') is None:
            f_atomic = hydrogen_carrier_split(pc, counts)
        else:
            f_atomic = CARRIER_SHARE_STATED['H']
        LAST_IMPOSITION['f_atomic'] = f_atomic
        if f_atomic > CARRIER_SHARE_TOL:
            pc.set_upper_bc('H', bc_type='flux', flux=phi_H*f_atomic)
            pc.set_upper_bc('H2', bc_type='flux',
                            flux=phi_H*(1.0 - f_atomic)/2.0)
            how = ('split H/H2 by the measured nuclei share, atomic %.4f'
                   % f_atomic)
        else:
            pc.set_upper_bc('H2', bc_type='flux', flux=phi_H/2.0)
            how = ('all on H2 (dominant carrier; atomic H carries %.3e of '
                   'the hydrogen nuclei at the top)' % f_atomic)
        imposed['H'] = phi_H
        imposed['H_atomic_share'] = f_atomic
        # The atomic share the boundary condition actually carries (0 on the
        # dominant-carrier branch); main() compares it with the converged top.
        LAST_IMPOSITION['H_atomic_share_used'] = (
            f_atomic if f_atomic > CARRIER_SHARE_TOL else 0.0)
        print('  imposed elemental H escape: %.6e g/s -> %.6e nuclei/cm^2/s '
              'at r_top = %.6e cm, %s' % (phi_H_gs, phi_H, r_top, how))
    if phi_He_gs != 0.0:
        phi_He = phi_He_gs/(sch.ATOMIC_WEIGHT['He']*sch.MAMU)/area
        pc.set_upper_bc('He', bc_type='flux', flux=phi_He)
        imposed['He'] = phi_He
        print('  imposed elemental He escape: %.6e g/s -> %.6e nuclei/cm^2/s'
              % (phi_He_gs, phi_He))
    LAST_IMPOSITION.update(phi_H_gs=phi_H_gs, phi_He_gs=phi_He_gs,
                           phi_H2_gs=phi_H2_gs,
                           r_top=(r_top if imposed else None))
    if phi_H2_gs is not None or phi_H_gs == 0.0:
        LAST_IMPOSITION['f_atomic'] = None
    return imposed, r_top


def impose_escape_flux_at_current_top(pc, counts):
    """Re-impose the escape fluxes of the most recent imposition with the
    top-face radius of the current grid, the same stated mass fluxes and the
    same hydrogen split (impose_elemental_escape_flux)."""
    return impose_elemental_escape_flux(
        pc, counts, LAST_IMPOSITION['phi_H_gs'], LAST_IMPOSITION['phi_He_gs'],
        LAST_IMPOSITION['phi_H2_gs'],
        atomic_share=LAST_IMPOSITION.get('f_atomic'))


def require_flux_imposed_at_final_top(imposed, r_top):
    """Refuse unless the escape flux was last imposed at the final top-face
    radius r_top [cm].

    The imposed numbers are flux densities: the stated mass fluxes over
    4 pi r^2 of the face they were imposed at.  The accepted state carries
    the stated mass fluxes only if that face is the final one, which
    steady_state_with_top_face_at_match guarantees for every state it
    accepts; this checks it on the state that is written.  No condition was
    imposed when `imposed` is empty, and nothing is checked then."""
    if not imposed:
        return
    r_bc = LAST_IMPOSITION.get('r_top')
    if r_bc != r_top:
        sch.refuse('the escape flux was last imposed at a top-face radius of'
                   ' %r cm, but the final top face stands at %r cm, so the'
                   ' state does not carry the stated mass fluxes'
                   % (r_bc, r_top))


def escape_flux_state(pc, counts):
    """One log line: the top-face radius now, the radius the boundary
    condition was imposed with, and the stated and carried H, He and H2
    mass fluxes [g/s] of the current solution."""
    measured, r_now = measured_elemental_escape_flux(pc, counts)
    _, top = pc.gas_fluxes()
    h2 = (float(top.get('H2', 0.0))*2.0*4.0*np.pi*r_now*r_now
          * sch.ATOMIC_WEIGHT['H']*sch.MAMU)
    r_bc = LAST_IMPOSITION.get('r_top')
    return ('r_top %.9e cm, imposed at %s; H stated %.9e carried %.9e, He'
            ' stated %.9e carried %.9e, H2 stated %s carried %.9e g/s'
            % (r_now, 'none' if r_bc is None else '%.9e cm' % r_bc,
               LAST_IMPOSITION.get('phi_H_gs', 0.0), measured['H'],
               LAST_IMPOSITION.get('phi_He_gs', 0.0), measured['He'],
               'none' if LAST_IMPOSITION.get('phi_H2_gs') is None
               else '%.9e' % LAST_IMPOSITION['phi_H2_gs'], h2))


def measured_elemental_escape_flux(pc, counts, elements=('H', 'He')):
    """The elemental fluxes the solution carries across its own top [g/s].

    `gas_fluxes()` gives the top-of-atmosphere flux of every long-lived
    species in molecules/cm^2/s, outward positive, as a dict keyed by
    species name; summing by nuclei count and multiplying by
    `4 pi r_top^2 m_El` turns it into the same quantity the closure iterates
    on.  This is the lower model's own statement of what it is losing, and
    it is a measurement, not the imposed number.
    """
    _, top_fluxes = pc.gas_fluxes()
    r_top = model_top_radius(pc)
    area = 4.0*np.pi*r_top*r_top
    out = {}
    for el in elements:
        s = 0.0
        for sp, phi in top_fluxes.items():
            k = counts.get(sp, {}).get(el, 0)
            if k > 0:
                s += k*float(phi)
        out[el] = s*area*sch.ATOMIC_WEIGHT[el]*sch.MAMU
    return out, r_top


def photochemical_steady_state(args, wdir, flux_file, mech, thermo,
                               data_dir, abundances, column=None):
    """Run the model to steady state, truncating the column from above only
    as far as the thermodynamic data force (section 4.2, trap 2).

    `column` is `(P [dyn/cm^2], T [K], K_zz [cm^2/s] or None)` deep to
    shallow when the caller has already solved for it (`--climate`);
    otherwise the prescribed `--tp-file` is read.
    """
    from photochem.extensions import gasgiants
    from photochem import PhotoException

    if column is None:
        P, T, Kzz_in = read_tp_file(args.tp_file, args.p_unit)
    else:
        P, T, Kzz_in = column
    Kzz_stated = sch.eddy_diffusion_coefficient(P/sch.BAR, args)
    if Kzz_stated is not None:
        # A stated eddy coefficient is the one the chemistry is solved on,
        # not only the one the file reports: K_zz sets the quench levels, so
        # a column solved on one profile and handed over with another would
        # carry a composition nobody's K_zz produced.
        Kzz_in = Kzz_stated
    elif Kzz_in is None:
        sch.refuse('the T(p) file carries no K_zz column: state '
                   '--kzz-const, or --kzz-power with --kzz-ref')

    mp = args.mp*sch.MJ
    rp = args.r_ref*sch.RJ
    truncation = ''
    # The pressure the gas-giant stepper holds its top CELL CENTRE near on
    # the first solve: the matching level.  The face is then placed at the
    # matching level by steady_state_with_top_face_at_match.
    toa = args.p_match*sch.BAR

    # The thermodynamic data cap the temperature (NH2 fails above ~1900 K on
    # some sets).  The column is TRUNCATED, never clipped: a clipped T(p) is
    # a profile nobody solved.  Each retry drops the levels above the cap and
    # lowers the model top with them.
    caps = ([args.t_max] if args.t_max
            else [None, 4000.0, 3000.0, 2500.0, 2000.0, 1800.0])
    last_error = None
    for cap in caps:
        Pc, Tc, Kc = P, T, Kzz_in
        if cap is not None:
            keep = Tc <= cap
            # keep the deep part contiguous: truncate from the top down
            if not keep[0]:
                sch.refuse('the deepest level is already above the '
                           'temperature the thermodynamic data admit')
            n = int(np.argmax(~keep)) if (~keep).any() else len(Tc)
            Pc, Tc, Kc = Pc[:n], Tc[:n], Kc[:n]
            if len(Pc) < 10:
                sch.refuse('the temperature truncation leaves fewer than ten '
                           'levels')
        # The photochemical grid must sit above the climate grid:
        # 3*target_pressure < P_climate_top (trap 1).  Set the top AND cut
        # the grid, never one alone -- the stated TOA is kept and the climate
        # column is cut back to it, which is what makes the requirement hold
        # however shallow the supplied T(p) happened to reach.
        toa_try = toa
        keep = Pc >= 3.05*toa_try
        Pc, Tc, Kc = Pc[keep], Tc[keep], Kc[keep]
        if len(Pc) < 10:
            sch.refuse('the requested model top %.3e bar is not reachable: '
                       'cutting the column at 3.05 x TOA (%.3e dyn/cm^2) '
                       'leaves %d levels, fewer than the ten a solution '
                       'needs. Ask for a deeper top, or supply a T(p) or a '
                       'climate column that reaches above it'
                       % (toa_try/sch.BAR, 3.05*toa_try, len(Pc)))
        try:
            pc = gasgiants.EvoAtmosphereGasGiant(
                mech, flux_file, mp, rp, solar_zenith_angle=args.zenith,
                thermo_file=thermo, data_dir=data_dir)
            pc.gdat.verbose = False
            # The model-top pressure the gas-giant initialization places the
            # top at and the robust stepper holds it within a factor
            # pressure_factor (3, set by EvoAtmosphereGasGiant) of.
            pc.var.toa_pressure_maintenance.target_pressure = toa_try
            # Where the photochemical model bottom sits, in units of the
            # pressure of the deepest quench level.  The gas-giant helper's
            # own 5 assumes a hot deep atmosphere in which the chemistry
            # equilibrates below the model bottom; on a cold planet no
            # quench level lies inside any column the climate model can
            # reach (the deepest level of the grid IS returned as the quench
            # level), and the factor has to be 1 or the model asks for a
            # bottom deeper than the column it was given.
            pc.gdat.BOA_pressure_factor = args.boa_pressure_factor
            atoms = list(pc.gdat.gas.atoms_names)
            missing = [a for a in atoms if a not in abundances]
            if missing:
                sch.refuse('no abundance was stated for %s, which the '
                           'mechanism carries' % ', '.join(missing))
            tot = sum(abundances[a] for a in atoms)
            # BY NAME: the atom order differs from run to run (trap 4).
            pc.gdat.gas.molfracs_atoms_sun = np.array(
                [abundances[a]/tot for a in atoms])
            pc.initialize_to_climate_equilibrium_PT(Pc, Tc, Kc, 1.0, 1.0)
            if FACE_GEOMETRY.get('photochem_diurnal_fac') is not None:
                # the face convention's diurnal factor (face_geometry)
                pc.var.diurnal_fac = FACE_GEOMETRY['photochem_diurnal_fac']
        except PhotoException as exc:
            last_error = str(exc)
            print('  photochem refused the column (%s); truncating' % exc)
            continue
        except Exception as exc:
            if ('does not extend deeply enough for the requested '
                    'photochemical lower boundary') not in str(exc):
                raise
            sch.refuse(
                'the photochemical model asks for a bottom %.1fx deeper in '
                'pressure than its deepest quench level, which is below the '
                'climate column. On an atmosphere cold enough that the '
                'chemistry never equilibrates inside the column, the '
                'deepest quench level is the column bottom itself and '
                '--boa-pressure-factor must be 1 (it is %.1f here); on a hot '
                'one, deepen --climate-p-deep instead'
                % (args.boa_pressure_factor, args.boa_pressure_factor))
        if cap is not None:
            truncation = ('T(p) truncated at %.0f K, model top %.3e dyn/cm2'
                          % (cap, Pc.min()))
        toa = toa_try
        break
    else:
        sch.refuse('photochem refused every truncation of the column; last '
                   'error: %s' % last_error)

    print('  climate grid %.3e -> %.3e dyn/cm2, %d levels, T_top = %.0f K'
          % (Pc[0], Pc[-1], len(Pc), Tc[-1]))
    print('  atoms: %s  species: %d'
          % (list(pc.gdat.gas.atoms_names), len(pc.dat.species_names)))
    sys.stdout.flush()

    # The temperature the chemistry runs on, up to its top face: the
    # reference column itself, whole (the gas-giant initialization holds T
    # constant above the climate grid it was handed, which stops at 3 x its
    # top pressure).
    set_chemistry_temperature(pc, chemistry_temperature_profile(
        args, (P, T, Kzz_in)))

    # The closure's boundary condition, imposed before the column is stepped
    # so that the steady state it reaches is the state of a losing column.
    imposed, r_top_imposed = impose_elemental_escape_flux(
        pc, species_element_counts(pc), args.trial_flux_H, args.trial_flux_He,
        args.trial_flux_H2)

    def report(n, give_up, reached):
        print('  block %d: t = %.3e s%s%s'
              % (n + 1, pc.wrk.tn, '  STEADY' if reached else '',
                 '  GAVE UP' if give_up else ''), flush=True)
    give_up, reached = steady_state_with_stall_retry(pc, args.blocks,
                                                     'the first solve', report)

    if give_up or not reached:
        sch.refuse('the chemistry did not reach steady state, nor on the '
                   'retry from the stalled state at rtol %.0e '
                   '(reached_steady_state = %s, gave_up = %s). A non-steady '
                   'solution has no elemental flux to hand over'
                   % (STALL_RETRY_RTOL, reached, give_up))

    steady_state_with_top_face_at_match(pc, args.p_match*sch.BAR,
                                        args.blocks,
                                        passes=args.top_pressure_passes,
                                        solve='the first solve')
    return pc, reached, truncation, toa, imposed, r_top_imposed


# THE STEP BUDGET.  `--blocks` B is B blocks of 100 robust steps, and that is
# the budget every steady-state pass runs under: Photochem's own cap
# `var.nsteps_before_giveup` (100000 as EvoAtmosphereGasGiant sets it) is
# set to 100 B, so
# the stepper gives up on the stated budget and not on a cap the command
# line does not show.
#
# A STALLED PASS IS RETRIED ONCE.  A pass that ends without a steady state
# (the budget spent, or Photochem's `give_up`) is integrated once more from
# the state the stalled integration holds, with the relative tolerance
# STALL_RETRY_RTOL, and is accepted only under the unchanged criteria: the
# stepper's own steady state (longdy below conv_longdy after min_step_conv
# steps, or the equilibrium time reached, with the top pressure in range),
# and after it the top-face tolerance and the measured-flux and elemental
# conservation checks that every solution passes.  A second stall refuses.
# The retry moves the integrator off a period-2 limit cycle of the BDF
# steps across the switch of the condensation term (its Jacobian is
# discontinuous at the saturation ratio RHc), which keeps longdy flat at a
# fixed step without a steady state (md/k05_chemistry_failure_20261008.md:
# on the LHS 1140 b k = 5 column the retry reached steady state in 126 steps,
# with a handed-over state equal to that of every other completed path).
STALL_RETRY_RTOL = 1.0e-4


def steady_state_retries(pc):
    """The record of the retried passes of this solution (a list of pass
    labels), kept on the Photochem object."""
    if not hasattr(pc, 'steady_state_retries'):
        pc.steady_state_retries = []
    return pc.steady_state_retries


def robust_steps_to_steady_state(pc, blocks, report=None):
    """(give_up, reached_steady_state) of at most `blocks` blocks of 100
    robust steps from the state pc holds, Photochem's step cap set to the
    same budget.  `report(n, give_up, reached)` is called after each block."""
    pc.var.nsteps_before_giveup = 100*blocks
    pc.initialize_robust_stepper(pc.wrk.usol)
    give_up = reached = False
    for n in range(blocks):
        for _ in range(100):
            give_up, reached = pc.robust_step()
            if give_up or reached:
                break
        if report is not None:
            report(n, give_up, reached)
        if give_up or reached:
            break
    return give_up, reached


def steady_state_with_stall_retry(pc, blocks, label, report=None):
    """`robust_steps_to_steady_state`, and on a stall the single retry from
    the stalled state at STALL_RETRY_RTOL (the block comment above).  A
    retry is recorded as `label` in `steady_state_retries(pc)` and logged as
    `steady_state_retry rtol=1e-4 at pass <label>`."""
    give_up, reached = robust_steps_to_steady_state(pc, blocks, report)
    if reached and not give_up:
        return give_up, reached
    print('  steady_state_retry rtol=%.0e at pass %s: the pass stalled at'
          ' t = %.3e s (reached_steady_state = %s, gave_up = %s, %d steps'
          ' against a budget of %d); integrating once more from the stalled'
          ' state' % (STALL_RETRY_RTOL, label, float(pc.wrk.tn), reached,
                      give_up, int(pc.wrk.nsteps_total), 100*blocks),
          flush=True)
    rtol_stated = pc.var.rtol
    pc.var.rtol = STALL_RETRY_RTOL
    try:
        give_up, reached = robust_steps_to_steady_state(pc, blocks, report)
    finally:
        pc.var.rtol = rtol_stated
    steady_state_retries(pc).append(label)
    print('  steady_state_retry rtol=%.0e at pass %s: %s after %d steps,'
          ' t = %.3e s' % (STALL_RETRY_RTOL, label,
                           'STEADY' if reached and not give_up
                           else 'stalled again',
                           int(pc.wrk.nsteps_total), float(pc.wrk.tn)),
          flush=True)
    return give_up, reached


# The top face of the chemistry grid sits at the matching level to this
# tolerance in ln p.  The state handed over at the match is the top face's
# (level_at_top_face), so a face 1e-3 in ln p off the match moves it by 1e-3
# of a scale height, below every gradient the closure resolves.
TOP_FACE_TOL = 1.0e-3
TOP_PRESSURE_PASSES = 8


def top_face_pressure(pc):
    """(p_face, H): the pressure at the top FACE of the chemistry grid
    [dyn/cm^2], `var.top_atmos`, continued from the two top cell centres
    with their own logarithmic slope (the grid is uniform in altitude), and
    the pressure scale height there [cm]."""
    p = np.asarray(pc.wrk.pressure_hydro, dtype=float)
    z = np.asarray(pc.var.z, dtype=float)
    slope = (np.log(p[-1]) - np.log(p[-2]))/(z[-1] - z[-2])
    p_face = float(np.exp(np.log(p[-1])
                          + slope*(float(pc.var.top_atmos) - z[-1])))
    return p_face, -1.0/slope


def steady_state_with_top_face_at_match(pc, p_face_dyn, blocks,
                                        tol=TOP_FACE_TOL,
                                        passes=TOP_PRESSURE_PASSES,
                                        solve='the first solve'):
    """Re-solve the steady state on the grid whose top FACE is the matching
    level p_face_dyn [dyn/cm^2].

    Photochem lays its levels out uniformly in ALTITUDE between a pinned
    bottom and `top_atmos` (`vertical_grid`, `photochem/src/photochem/photochem_eqns.f90`),
    so the altitude of the top face fixes the pressure of every level, and
    the converged composition moves the pressure at a fixed altitude through
    the hydrostatic integration.  The map is therefore iterated: measure the
    face pressure of the converged solution (`top_face_pressure`), move the
    face by H ln(p_face/p_match), re-converge, until the face is at the
    match to `tol` in ln p.  Where two successive moves have opposite signs
    the step is halved (the map's slope is then near -1).

    Photochem's TOA-pressure maintenance is disabled for these passes and
    its previous state restored on every exit (success, refusal or
    exception), so this iteration is the only thing that places the top:
    with `toa_pressure_maintenance.enabled` false, initialize_robust_stepper
    skips its band check and relocation
    (photochem/src/photochem/photochem_evoatmosphere_integrate.f90:1010,
    1023-1036), toa_pressure_state returns WithinTol (:1271), so neither the
    extreme nor the post-convergence resynchronization is triggered by the
    top (:1210-1235), and resync_and_restart, reached at the scheduled
    restarts, updates the grid only when maintenance is enabled
    (:1392-1404).  Inside the passes the grid therefore stays fixed between
    the moves made here; the band maintenance would otherwise relocate the
    top to its target whenever the top-centre gas pressure (`wrk.pressure`,
    not the hydrostatic pressure this iteration measures, :1273) left a
    factor `pressure_factor` about the target at a convergence check.  The
    target itself is set to the top-centre hydrostatic pressure after each
    move and kept on exit, for the solves that follow with maintenance on.

    Each re-convergence is one steady-state pass under the step budget and
    the single stall retry of `steady_state_with_stall_retry`; `solve` names
    the solve the passes follow, for the record of a retry.

    THE ESCAPE FLUX FOLLOWS THE FACE.  The flux densities of the escape
    boundary condition are the stated mass fluxes over the area of the face
    they were imposed at (impose_elemental_escape_flux).  The face moves
    here, and it can also have moved before these passes, in the solve
    that ran with the TOA-pressure maintenance on.  Every pass therefore
    compares the current top radius with the radius of the last
    imposition, re-imposes the condition at the current radius where they
    differ, and re-converges; a state is accepted only with the face at the
    match AND the condition imposed at the radius it stands at, so that the
    mass flux it carries is the stated one to the convergence of the
    chemistry.

    FAILURE IS A REFUSAL: a column whose top will not settle at the match is
    not the column the matching surface needs.
    """
    from photochem import PhotoException

    # Everything that can raise before the passes runs before the
    # maintenance is disabled, so the finally below is the one place that
    # restores it.
    counts = species_element_counts(pc)
    maintenance = pc.var.toa_pressure_maintenance
    maintenance_enabled = bool(maintenance.enabled)
    moves = []
    maintenance.enabled = False
    try:
        for attempt in range(1, passes + 1):
            p_face, H = top_face_pressure(pc)
            off = float(np.log(p_face/p_face_dyn))
            r_bc = LAST_IMPOSITION.get('r_top')
            print('  escape flux at pass %d: %s'
                  % (attempt, escape_flux_state(pc, counts)), flush=True)
            print('  top face pass %d: face %.6e dyn/cm2 (%+.3e in ln p from'
                  ' the match), top centre %.6e dyn/cm2 (gas pressure'
                  ' %.6e), TOA target %.6e dyn/cm2, maintenance off'
                  % (attempt, p_face, off, float(pc.wrk.pressure_hydro[-1]),
                     float(pc.wrk.pressure[-1]),
                     float(maintenance.target_pressure)),
                  flush=True)
            at_face = abs(off) <= tol
            if at_face and (r_bc is None or r_bc == model_top_radius(pc)):
                return p_face, attempt
            if not at_face:
                step = H*off
                if len(moves) >= 1 and moves[-1]*step < 0.0:
                    step *= 0.5
                moves.append(step)
                try:
                    pc.update_vertical_grid(
                        TOA_alt=float(pc.var.top_atmos) + step)
                except PhotoException as exc:
                    sch.refuse('the top face could not be moved to the '
                               'matching level: %s' % exc)
                maintenance.target_pressure = float(
                    pc.wrk.pressure_hydro[-1])
            if r_bc is not None and r_bc != model_top_radius(pc):
                impose_escape_flux_at_current_top(pc, counts)
            give_up, reached = steady_state_with_stall_retry(
                pc, blocks, 'top face %d after %s' % (attempt, solve))
            if give_up or not reached:
                sch.refuse('the chemistry did not re-reach steady state with '
                           'its top face moved toward the matching level, '
                           'nor on the retry from the stalled state at rtol '
                           '%.0e (pass %d of %d, reached_steady_state = %s, '
                           'gave_up = %s)'
                           % (STALL_RETRY_RTOL, attempt, passes, reached,
                              give_up))
        p_face, _ = top_face_pressure(pc)
        sch.refuse('the top face will not settle at the matching level: after'
                   ' %d passes it stands %+.3e in ln p from it, against a'
                   ' tolerance of %.1e' % (passes,
                                           np.log(p_face/p_face_dyn), tol))
    finally:
        maintenance.enabled = maintenance_enabled


def chemistry_temperature_profile(args, column, res=None):
    """(P [dyn/cm^2], T [K], K_zz [cm^2/s]) deep to shallow, ending at the
    matching level, for `set_chemistry_temperature`: the reference column
    (P, T, K or None) with, under a conducted flux, the rise of the column
    solution `res` (`temperature_rise_at`) and its nodes added, so that the
    warmed layer under the match is on the profile the chemistry reads."""
    P, T, K = column
    pm = args.p_match*sch.BAR
    lnP = np.log(np.asarray(P, dtype=float))
    x = lnP[lnP > np.log(pm) + 1.0e-9]
    if res is not None:
        xd = np.log(res['domain']['p_bar']*sch.BAR)
        x = np.concatenate([x, xd[(xd > np.log(pm) + 1.0e-9) &
                                  (xd < lnP[0])]])
    x = np.unique(np.concatenate([x, [np.log(pm)]]))
    x = x[np.concatenate([[True], np.diff(x) > 1.0e-9])][::-1]
    Tq = np.interp(x, lnP[::-1], np.asarray(T, dtype=float)[::-1])
    Pq = np.exp(x)
    if np.any(np.diff(Pq) >= 0.0):
        keep = np.concatenate([[True], np.diff(Pq) < 0.0])
        x, Pq, Tq = x[keep], Pq[keep], Tq[keep]
    if res is not None:
        Tq = Tq + temperature_rise_at(res, Pq)
    Kq = sch.eddy_diffusion_coefficient(Pq/sch.BAR, args)
    if Kq is None:
        Kq = np.interp(x, lnP[::-1], np.asarray(K, dtype=float)[::-1])
    return Pq, Tq, np.asarray(Kq, dtype=float)*np.ones_like(Pq)


def set_chemistry_temperature(pc, profile):
    """Hand Photochem the temperature and eddy profiles it runs on as a
    persistent T(p), K_zz(p) (`set_press_temp_edd_profile`), re-mapped onto
    every grid it takes afterwards.  Above the profile's top (the matching
    level) Photochem holds the top value, but no cell centre lies there:
    the top face is the match.

    The profile is mapped continuously (mode 1: at every right-hand-side
    evaluation, so T and K_zz stay functions of pressure as the composition
    moves the hydrostatic pressure), and the robust stepper keeps its
    TOA-pressure target: the call would otherwise reset it to Photochem's
    default 0.1 dyn/cm^2."""
    P, T, K = profile
    try:
        pc.destroy_stepper()
    except Exception:
        pass
    target = float(pc.var.toa_pressure_maintenance.target_pressure)
    pc.set_press_temp_edd_profile(np.ascontiguousarray(P),
                                  np.ascontiguousarray(T),
                                  np.ascontiguousarray(K),
                                  mode=1, maintain_toa_pressure=True,
                                  target_pressure=target)


def species_element_counts(pc):
    """Element counts BY NAME from Photochem's own composition matrix.

    Never parsed from the label: `O1D` and `N2D` are excited states, not
    formulas, and a formula parser would invent nuclei for them.
    """
    atoms = list(pc.dat.atoms_names)
    names = list(pc.dat.species_names)
    comp = np.asarray(pc.dat.species_composition)        # (natoms, nsp)
    counts = {}
    for j, sp in enumerate(names):
        if j >= comp.shape[1]:
            break
        c = {atoms[i]: int(comp[i, j]) for i in range(len(atoms))
             if comp[i, j] > 0}
        if c:
            counts[sp] = c
    return counts


# --------------------------------------------------------------------------
# The heat the escape wind conducts into the column (--conducted-heat-flux)
# --------------------------------------------------------------------------
# The column's own response to that flux is solved by
# conducted_heat_column_balance.py; this part couples it to the chemistry.
#
# Levels inserted into the T(p) column handed to the chemistry, so that the
# temperature structure the response puts below the match (e-folding depth
# 0.07-0.14 pressure scale heights on LHS 1140 b) reaches the chemistry's
# own grid rather than an interpolation across the climate grid's 0.36 H.
CHEMISTRY_LEVEL_SPACING = 0.05      # in ln p
CHEMISTRY_LEVEL_DEPTH = 8.0         # scale heights below p_match
# Levels of the written profile: every node of the column solution whose
# temperature rise exceeds this, so that the bath temperature EXHALE reads
# half a cell below its base face (linear in r between levels) is the
# solved one and not an interpolation across the chemistry grid.
PROFILE_LEVEL_DT_MIN = 1.0e-3       # K
COLUMN_TEMPERATURE_TOL = 1.0        # K, default of --column-energy-tolerance


def column_with_inserted_levels(P, T, Kzz, p_match_bar):
    """(P, T, Kzz) deep to shallow with levels added every
    CHEMISTRY_LEVEL_SPACING in ln p over CHEMISTRY_LEVEL_DEPTH scale heights
    below the match, T and K_zz linear in ln p (the climate stratosphere is
    an isotherm there, so the temperature profile is unchanged)."""
    pm = p_match_bar*sch.BAR
    x_new = np.log(pm) + np.arange(0.0, CHEMISTRY_LEVEL_DEPTH + 1e-9,
                                   CHEMISTRY_LEVEL_SPACING)
    lnP = np.log(P)
    keep = (x_new < lnP[0]) & (x_new > lnP[-1])
    x_new = x_new[keep]
    x_new = np.array([x for x in x_new
                      if np.min(np.abs(lnP - x)) > 1.0e-6])
    x_all = np.sort(np.concatenate([lnP, x_new]))[::-1]
    T_all = np.interp(x_all, lnP[::-1], T[::-1])
    K_all = (np.interp(x_all, lnP[::-1], Kzz[::-1]) if Kzz is not None
             else None)
    return np.exp(x_all), T_all, K_all


def temperature_rise_at(res, p_dyn):
    """The solved rise Delta T at pressures p_dyn [dyn/cm^2]: the column
    solution inside its domain, its value at the match above it (the
    reference column is an isotherm there), zero below it."""
    s = res['domain']['s']
    x = -np.log(np.asarray(p_dyn, dtype=float))
    out = np.interp(x, s, res['dT'], left=0.0, right=res['dT_match'])
    return out


# The top of the climate solve by default [dyn/cm^2]: 1e-8 bar, two decades
# above the matching level; only the climate column below the match is used.
CLIMATE_P_TOP_DYN = 1.0e-2

# The dissociation energy of H2 from its ground level, D0 = 36118.11 cm^-1
# (Huber & Herzberg 1979; the value molecular_reaction_heat.f90 uses for
# the same reaction), in erg.
D0_H2_ERG = 36118.11*6.62607015e-27*2.99792458e10


def production_rates_on_layers(production, nz, nprocesses):
    """Photochem's `ProductionLoss.production` [cm^-3 s^-1] as an (nz,
    nprocesses) array: one row for each layer, one column for each entry of
    `production_rx`, the layout the v0.9.3 interface declares
    (bindings/photochem/cython/ProductionLoss.pyx: "shape (nz,
    nproduction)").  The shape is checked against that declaration rather
    than inferred, so a square array (nz == nprocesses) is read in the
    declared orientation; any other shape is a refusal."""
    prod = np.asarray(production, dtype=float)
    if prod.ndim != 2 or prod.shape != (int(nz), int(nprocesses)):
        sch.refuse('Photochem returned a production array of shape %s; the '
                   'interface declares (nz, nproduction) = (%d, %d)'
                   % (prod.shape, int(nz), int(nprocesses)))
    return prod


def recombination_heating_of_inflowing_atoms(pc, args):
    """The heat the atomic hydrogen the wind mixes down through the matching
    surface releases in the column, as a `chemical_heating` of
    conducted_heat_column_balance (per node, in the units of A w D), and its
    record; (None, record) when no atoms come in.

    TOTAL: an inflow phi_in [atoms cm^-2 s^-1] of atomic H through the top
    face, at steady state, leaves the column as H2 (the elemental hydrogen
    flux is fixed), so by Hess's law the column receives
    F_chem = phi_in D0/2 [erg cm^-2 s^-1] whatever the path; the wind does
    not charge this energy, because the atoms leave its base cell before
    they recombine.  An outflow of atoms has no counterpart here: their
    dissociation was charged where it happened, in the wind.

    SHAPE: the column's own H + H (+ M) -> H2 (+ M) rate per unit mass,
    normalized to F_chem.  VALIDITY: that rate is dominated by the inflowing
    atoms wherever they raise q_H above its photolytic value, which is where
    the deposition matters; the record gives the ratio of the column's whole
    H + H rate times D0 to F_chem (1 when every recombination is of an
    inflowing atom, above 1 where photolytic atoms recombine as well, which
    then shares the shape but not the total)."""
    rec = dict(F_chem=0.0, phi_in=0.0, rate_ratio=float('nan'),
               reaction='')
    if args.trial_flux_H2 is None:
        return None, rec
    r_top = model_top_radius(pc)
    area = 4.0*np.pi*r_top*r_top
    phi_atomic = ((args.trial_flux_H - args.trial_flux_H2)
                  /(sch.ATOMIC_WEIGHT['H']*sch.MAMU)/area)
    if phi_atomic >= 0.0:
        return None, rec
    phi_in = -phi_atomic
    F_chem = phi_in*0.5*D0_H2_ERG
    pl = pc.production_and_loss('H2', pc.wrk.usol)
    names = [str(n).strip() for n in pl.production_rx]
    # the recombination itself, H + H (+ M) -> H2 (+ M); the abstractions
    # H + HX -> H2 + X also make H2 from an atom but release D0 - D(H-X)
    idx = [i for i, n in enumerate(names)
           if n.replace(' ', '').split('=>')[0] in ('H+H', 'H+H+M')]
    if not idx:
        sch.refuse('the mechanism has no H + H -> H2 reaction to place the '
                   'recombination heat of the inflowing atoms with')
    prod = production_rates_on_layers(pl.production, int(pc.var.nz),
                                      len(names))
    rate = np.sum(prod[:, idx], axis=1)                 # cm^-3 s^-1
    z = np.asarray(pc.var.z, dtype=float)
    dz = np.gradient(z)
    p_dyn = np.asarray(pc.wrk.pressure_hydro, dtype=float)
    rho = np.asarray(pc.wrk.density, dtype=float)*np.asarray(
        pc.wrk.mubar, dtype=float)*sch.MAMU
    per_mass = rate*D0_H2_ERG/rho                       # erg g^-1 s^-1
    rec.update(F_chem=float(F_chem), phi_in=float(phi_in),
               rate_ratio=float(np.sum(rate*dz)*D0_H2_ERG/F_chem),
               reaction='; '.join(names[i] for i in idx))
    lnp_asc = np.log(p_dyn)[::-1]
    pm_asc = per_mass[::-1]

    def heating(d):
        h = np.interp(np.log(d['p_bar']*sch.BAR), lnp_asc, pm_asc)
        H = d['geom'].area*d['geom'].w*h
        H[0] = 0.0
        tot = np.sum(H)
        if not tot > 0.0:
            sch.refuse('the column carries no H + H recombination below '
                       'the match to place the inflowing atoms\' heat with')
        return H*d['geom'].area[-1]*F_chem/tot
    return heating, rec


def conducted_heat_column_state(pc, args, counts, emission, T_reference,
                                p_floor_bar):
    """Solve the column's temperature under the conducted flux with the
    composition the chemistry currently holds.  T_reference(p_dyn) is the
    climate temperature: the reference state, never the warmed one."""
    import conducted_heat_column_balance as chb
    sol = pc.return_atmosphere()
    p_dyn = np.asarray(sol['pressure'], dtype=float)
    p_bar = p_dyn/sch.BAR
    nparticles = int(pc.dat.np)
    particles = set(list(pc.dat.species_names)[:nparticles])
    gas = {sp: np.asarray(y, dtype=float) for sp, y in sol.items()
           if sp in counts and sp not in particles}
    mu = sch.mean_molecular_weight(gas, counts)
    mix = {sp: gas[sp] for sp in ('He', 'H2', 'H') + chb.EMITTERS
           if sp in gas}
    if EDDY_HEAT.get('levels_dyn') is not None:
        # The reference column's own nodes are added to the chemistry levels
        # the column is built from, so that T0 at the column's nodes (the
        # same lattice, ColumnResponse dlnp) is the solved T0 itself and not
        # an interpolation of it across the chemistry grid; the composition
        # at the added levels is the column's own interpolation rule
        # (ColumnResponse._profile_at: He and H2 and the particle mass
        # linear in ln p, the rest log-linear), so the column sees the same
        # composition as without them.
        lv = EDDY_HEAT['levels_dyn']
        lv = lv[lv <= p_dyn.max()]
        x_old = np.log(p_dyn)
        o = np.argsort(x_old)
        x_new = np.unique(np.concatenate([x_old, np.log(lv)]))
        x_new = x_new[np.concatenate([[True], np.diff(x_new) > 1e-9])]
        mu = np.interp(x_new, x_old[o], np.asarray(mu, dtype=float)[o])
        mix = {sp: (np.interp(x_new, x_old[o], v[o]) if sp in ('He', 'H2')
                    else np.exp(np.interp(x_new, x_old[o],
                                          np.log(np.maximum(v[o], 1e-300)))))
               for sp, v in mix.items()}
        p_dyn = np.exp(x_new)
        p_bar = p_dyn/sch.BAR
    # The uppermost chemistry value is at a cell center. Continue its
    # composition to the physical matching face and evaluate the reference
    # temperature there before integrating dr/dln(p).
    p_face_dyn = args.p_match*sch.BAR
    if p_face_dyn < np.min(p_dyn):
        top = int(np.argmin(p_dyn))
        p_dyn = np.append(p_dyn, p_face_dyn)
        mu = np.append(mu, mu[top])
        mix = {sp: np.append(values, values[top]) for sp, values in mix.items()}
        order = np.argsort(p_dyn)[::-1]
        p_dyn, mu = p_dyn[order], mu[order]
        mix = {sp: values[order] for sp, values in mix.items()}
        p_bar = p_dyn/sch.BAR
        p_bar[-1] = args.p_match
    T0 = T_reference(p_dyn)
    r_cm = sch.hydrostatic_radius(p_bar, T0, mu, args.mp, args.r_ref,
                                  args.p_ref)*sch.RJ
    col = chb.ColumnResponse(p_bar, T0, r_cm, mu*sch.MAMU, mix,
                             args.mp*sch.MJ, args.p_match, emission=emission,
                             nh3_h2_rate_analog=args.nh3_h2_rate_analog,
                             kzz=EDDY_HEAT.get('kzz'),
                             reference_radiative_flux=EDDY_HEAT.get(
                                 'radiative_flux'))
    chem, chem_rec = recombination_heating_of_inflowing_atoms(pc, args)
    if args.column_band_transfer == 'thin':
        res = col.solve(args.conducted_heat_flux,
                        treatment=args.column_energy_treatment,
                        p_floor_bar=p_floor_bar, chemical_heating=chem)
    elif args.column_band_transfer == 'cooling_to_space_exchange':
        res = col.solve_with_included_exchange(
            args.conducted_heat_flux, treatment=args.column_energy_treatment,
            p_floor_bar=p_floor_bar, chemical_heating=chem)
    else:
        res = col.solve(args.conducted_heat_flux,
                        treatment=args.column_energy_treatment,
                        p_floor_bar=p_floor_bar,
                        transfer=args.column_band_transfer,
                        exchange_shift_max_K=args.column_exchange_shift_max,
                        chemical_heating=chem)
    res['recombination'] = chem_rec
    return res


def chemistry_on_new_temperature(pc, args, counts, Pc, Tc, Kc, profile):
    """Re-solve the chemistry on the temperature Tc of the climate grid Pc,
    starting from the solution it holds (`reinitialize_to_new_climate_PT`
    with the current mixing ratios), with the whole temperature profile up
    to the top face (`profile`, `chemistry_temperature_profile`) set as the
    one the chemistry runs on, the same escape boundary condition and the
    same placing of the top face as the first solve."""
    nparticles = int(pc.dat.np)
    names = list(pc.dat.species_names)
    nsl = int(pc.dat.nsl)
    gas_names = names[nparticles:len(names) - 2 - nsl]
    clim = pc.return_atmosphere_climate_grid()
    mix = {sp: np.asarray(clim[sp], dtype=float) for sp in gas_names}
    pc.reinitialize_to_new_climate_PT(Pc, Tc, Kc, mix)
    set_chemistry_temperature(pc, profile)
    impose_elemental_escape_flux(pc, counts, args.trial_flux_H,
                                 args.trial_flux_He, args.trial_flux_H2)
    n_solve = getattr(pc, 'new_temperature_solves', 0) + 1
    pc.new_temperature_solves = n_solve
    solve = 'new-temperature solve %d' % n_solve
    give_up, reached = steady_state_with_stall_retry(pc, args.blocks, solve)
    print('  chemistry on the warmed column: t = %.3e s%s%s'
          % (pc.wrk.tn, '  STEADY' if reached else '',
             '  GAVE UP' if give_up else ''), flush=True)
    if give_up or not reached:
        sch.refuse('the chemistry did not reach steady state on the column '
                   'warmed by the conducted heat flux, nor on the retry from '
                   'the stalled state at rtol %.0e (reached_steady_state '
                   '= %s, gave_up = %s)' % (STALL_RETRY_RTOL, reached,
                                            give_up))
    steady_state_with_top_face_at_match(pc, args.p_match*sch.BAR,
                                        args.blocks,
                                        passes=args.top_pressure_passes,
                                        solve=solve)


def reference_column_on_chemistry_grid(pc, args, column):
    """(T_reference(p), Pc, Tc, Kc): the reference column (P, T, K) [dyn/cm2,
    K, cm2/s], interpolated in ln p onto the climate grid Pc of the chemistry;
    K from the eddy-coefficient arguments when they state one."""
    P_ref, T_ref = column[0], column[1]
    lnP_ref = np.log(P_ref)[::-1]
    T_ref_asc = T_ref[::-1]

    def T_reference(p_dyn):
        return np.interp(np.log(p_dyn), lnP_ref, T_ref_asc)

    Pc = np.asarray(pc.gdat.P_clima_grid, dtype=float)
    Tc = T_reference(Pc)
    Kc = sch.eddy_diffusion_coefficient(Pc/sch.BAR, args)
    if Kc is None:
        Kc = np.interp(np.log(Pc), lnP_ref, column[2][::-1])
    return T_reference, Pc, Tc, Kc


def conducted_heat_iteration(pc, args, column, p_floor_bar):
    """Temperature and composition of the column under the conducted flux,
    iterated to agreement: column energy balance with the chemistry's
    composition, chemistry re-solved on the new temperature, until the
    matching-level temperature moves by less than
    --column-energy-tolerance [K].  Returns the last column solution and the
    record of the iterations."""
    import conducted_heat_column_balance as chb
    T_reference, Pc, Tc_ref, Kc = reference_column_on_chemistry_grid(
        pc, args, column)
    counts = species_element_counts(pc)
    emission = chb.LTEBandEmission()
    history = []
    previous = None
    for it in range(args.column_energy_max_iterations):
        res = conducted_heat_column_state(pc, args, counts, emission,
                                          T_reference, p_floor_bar)
        history.append(res['T_match'])
        print('  column energy balance %d: F = %.4g erg/cm2/s, treatment %s,'
              ' T_match = %.3f K (rise %.3f K, e-folding %.3f H, domain %.2f H,'
              ' %d Newton iterations)'
              % (it, args.conducted_heat_flux, args.column_energy_treatment,
                 res['T_match'], res['dT_match'], res['e_folding_H'],
                 res['depth_H'], res['newton_iterations']), flush=True)
        tl = res.get('thin_layer')
        if tl is not None and tl['exceeded'] and \
                args.column_band_transfer == 'thin':
            print('  WARNING: ' + tl['message'] + ' (thin_layer_limit_exceeded'
                  ' T; largest absorbed fraction %.4f, %s %s bin %.0f-%.0f'
                  ' cm^-1)' % (tl['max_absorbed'], tl['species'], tl['band'],
                               tl['bin_cm'][0], tl['bin_cm'][1]), flush=True)
        if 'exchange' in res and res['exchange'].get('included'):
            print('    transfer %s: the exchange inside the warmed layer is'
                  ' solved, T_match %.3f K above the cooling-to-space value'
                  % (args.column_band_transfer, res['exchange']['shift_K']),
                  flush=True)
        elif 'exchange' in res:
            print('    transfer %s: neglected exchange shifts T_match by'
                  ' %.3f K to first order (bound %.3g K)'
                  % (args.column_band_transfer, res['exchange']['shift_K'],
                     res['exchange']['shift_max_K']), flush=True)
        if previous is not None and \
                abs(res['T_match'] - previous) < args.column_energy_tolerance:
            res['iterations'] = it + 1
            res['T_match_history'] = history
            res['T_reference_fn'] = T_reference
            return res
        previous = res['T_match']
        Tc_new = Tc_ref + temperature_rise_at(res, Pc)
        chemistry_on_new_temperature(
            pc, args, counts, Pc, Tc_new, Kc,
            chemistry_temperature_profile(args, column, res))
    sch.refuse('the column temperature and the chemistry did not agree to '
               '%.3g K at the match in %d iterations (T_match history %s)'
               % (args.column_energy_tolerance,
                  args.column_energy_max_iterations,
                  ', '.join('%.2f' % t for t in history)))


def column_heat_note(args, res):
    """The notes sentence of a profile written under a conducted flux."""
    shares = sorted(res['band_shares'], key=lambda r: -abs(r[4]))[:3]
    rec = res.get('recombination') or {}
    chem = (' plus %.4g erg/cm2/s recombination heat of the inflowing'
            ' atomic H' % rec['F_chem']) if rec.get('F_chem', 0.0) > 0.0 \
        else ''
    head = ('Column energy balance below p_match under the conducted heat'
            ' flux %.4g erg/cm2/s' % args.conducted_heat_flux) + chem
    if args.eddy_heat_transport == 'on':
        head += ', with the eddy heat flux rho c_p K_zz (dT/dz + g/c_p)'
    return head + (' (conducted_heat_column_balance.py, %s%s):'
                   ' T_match %.2f K, rise %.2f K, e-folding %.3f H; bands %s;'
                   ' %d chemistry-temperature iterations.'
                   % (args.column_energy_treatment,
                      '' if args.column_band_transfer == 'thin'
                      else ', transfer ' + args.column_band_transfer,
                      res['T_match'], res['dT_match'], res['e_folding_H'],
                      ', '.join('%s %s %.0f%%' % (r[0], r[1], 100*r[4])
                                for r in shares), res['iterations']))


def write_column_heat_record(args, res):
    """`conducted_heat_column.txt` beside the profile: the column solution
    on its own nodes, with the band shares, the optical depths and the
    energy budget in the header."""
    d = res['domain']
    path = os.path.join(args.run_dir, 'conducted_heat_column.txt')
    q = res['q']
    q_nodes = np.concatenate([[q[0]], 0.5*(q[1:] + q[:-1]),
                              [args.conducted_heat_flux]])
    with open(path, 'w') as f:
        f.write('# column energy balance under the conducted heat flux'
                ' (conducted_heat_column_balance.py)\n')
        f.write('# conducted_heat_flux_cgs %.10E\n' % args.conducted_heat_flux)
        f.write('# treatment %s\n' % args.column_energy_treatment)
        if args.column_energy_treatment == 'nlte':
            f.write('# nh3_h2_rate_analog %s\n' % args.nh3_h2_rate_analog)
        if args.column_band_transfer == 'thin':
            tl = res['thin_layer']
            f.write('# thin_layer_limit_exceeded %s\n'
                    % ('T' if tl['exceeded'] else 'F'))
            f.write('# thin_layer_max_absorbed_fraction %.6E (limit %.2f) at'
                    ' %s %s bin %.0f-%.0f cm^-1\n'
                    % (tl['max_absorbed'], tl['limit'], tl['species'],
                       tl['band'].replace(' ', '_'), tl['bin_cm'][0],
                       tl['bin_cm'][1]))
        if args.column_band_transfer != 'thin':
            f.write('# transfer %s\n' % args.column_band_transfer)
            x = res['exchange']
            if x.get('included'):
                f.write('# exchange included: T_match shift %.4f K from the'
                        ' exchange (solved, %d iterations), heating/F %.4e,'
                        ' nodes from %.4e bar\n'
                        % (x['shift_K'], x['iterations'],
                           x['heating_over_F'], x['p_k0_bar']))
            else:
                f.write('# exchange neglected: first-order T_match shift'
                        ' %.4f K (bound %.3g K), heating/F %.4e, nodes from'
                        ' %.4e bar\n'
                        % (x['shift_K'], x['shift_max_K'],
                           x.get('heating_over_F', 0.0),
                           x.get('p_k0_bar', 0.0)))
            for e in res['escape']:
                f.write('# escape %-5s %-22s beta(match) %.5f'
                        ' beta_min(warmed) %.5f\n'
                        % (e['species'], e['band'], e['beta_match'],
                           e['beta_min_warm']))
        rec = res.get('recombination')
        if rec is not None and rec['F_chem'] > 0.0:
            f.write('# recombination heat of the inflowing atomic hydrogen:'
                    ' %.6E erg/cm2/s (%.6E atoms/cm2/s in, D0/2 each),'
                    ' deposited with the shape of %s; column rate x D0 over'
                    ' it %.4f\n' % (rec['F_chem'], rec['phi_in'],
                                    rec['reaction'], rec['rate_ratio']))
        f.write('# eddy_heat_transport %s\n' % args.eddy_heat_transport)
        fl = res.get('face_fluxes_below_match')
        if fl is not None:
            f.write('# flux through the face below the match [erg/cm2/s,'
                    ' down positive]: molecular %.6E eddy %.6E (reference T0:'
                    ' molecular %.6E eddy %.6E) at %.6E bar\n'
                    % (fl['mol'][-1], fl['eddy'][-1], fl['mol0'][-1],
                       fl['eddy0'][-1], fl['p_face'][-1]))
        ri = res.get('reference_imbalance')
        if ri is not None:
            f.write('# reference imbalance max|S0| %.6E r^2 erg/cm2/s'
                    ' (zero where T0 solves its own energy balance with this'
                    ' column\'s conduction); held reference radiative'
                    ' heating %.6E\n'
                    % (ri['max_abs_S0'],
                       res['budget'].get('reference_radiative_heating', 0.0)))
        f.write('# T_match_K %.10E\n' % res['T_match'])
        f.write('# dT_match_K %.10E\n' % res['dT_match'])
        f.write('# e_folding_H %.6f\n' % res['e_folding_H'])
        f.write('# domain_depth_H %.4f (limited by the floor: %s)\n'
                % (res['depth_H'], res['depth_limited_by_floor']))
        f.write('# iterations %d, T_match history %s\n'
                % (res['iterations'],
                   ' '.join('%.4f' % t for t in res['T_match_history'])))
        b = res['budget']
        f.write('# budget [r^2 erg cm^-2 s^-1, r in cm]: emission excess %.10E,'
                ' bottom flux %.10E, top flux %.10E, relative error %.3E\n'
                % (b['excess'], b['bottom_flux'], b['top'],
                   b['relative_error']))
        for r in res['band_shares']:
            f.write('# band %-5s %-22s %-4s share %.5f\n'
                    % (r[0], r[1], r[2], r[4]))
        for o in res['optical_depth']:
            f.write('# tau %-5s %-22s tau_bar(match) %.3e tau_max_g(match)'
                    ' %.3e absorbed fraction of the excess %.3e\n'
                    % (o['species'], o['band'], o['tau_bar_at_match'],
                       o['tau_max_g_at_match'],
                       o['absorbed_fraction_of_excess']))
        f.write('# columns: p[bar] T_reference[K] T[K] dT[K]'
                ' q_down[erg/cm2/s]\n')
        for i in range(d['s'].size):
            f.write('%25.17E %25.17E %25.17E %25.17E %25.17E\n'
                    % (d['p_bar'][i], d['T0'][i], res['T'][i], res['dT'][i],
                       q_nodes[i]))
    print('wrote %s' % path)


def level_at_top_face(cols, p_match_bar):
    """`cols` (deep to shallow) with a level at the matching level p_match,
    the top face of the chemistry grid, half a cell above its top centre.

    The face state is continued from the two top cell centres, linearly in
    ln p: ln q for every positive mixing ratio and element ratio (q_H rises
    steeply toward a face the wind's atoms come in through, and a linear q
    could cross zero), the temperature, the radius and K_zz themselves, the
    mean particle mass rho/n_tot; then n_tot = p/(k T) and rho = n_tot times
    that mass, so the gas stays ideal at the new level.  F_ columns are
    constants.  A table that already reaches the match (or above it) gets
    the level by interpolation (insert_level).  Under a conducted flux the
    temperature and the radius are then replaced by the column solution's
    (profile_with_column_temperature), which is defined up to the match."""
    p = cols['p']
    if np.any(p == p_match_bar) or p[-1] < p_match_bar:
        return sch.insert_level(cols, p_match_bar)
    w = (np.log(p_match_bar) - np.log(p[-1]))/(np.log(p[-1]) - np.log(p[-2]))
    out = {}
    for name, y in cols.items():
        y = np.asarray(y, dtype=float)
        a, b = y[-2], y[-1]
        if name == 'p':
            v = p_match_bar
        elif name.startswith('F_'):
            v = b
        elif (name.startswith('q_') or name.startswith('X_')) and a > 0.0 \
                and b > 0.0:
            v = float(np.exp(np.log(b) + w*(np.log(b) - np.log(a))))
        else:
            v = b + w*(b - a)
        out[name] = np.append(y, v)
    if 'n_tot' in out and 'rho' in out and out['T'][-1] > 0.0:
        m_a = cols['rho'][-2]/cols['n_tot'][-2]
        m_b = cols['rho'][-1]/cols['n_tot'][-1]
        out['n_tot'][-1] = p_match_bar*sch.BAR/(sch.KB*out['T'][-1])
        out['rho'][-1] = out['n_tot'][-1]*(m_b + w*(m_b - m_a))
    return out


def profile_with_column_temperature(cols, res, args):
    """The written profile under the conducted flux: levels added at the
    nodes of the column solution whose rise exceeds PROFILE_LEVEL_DT_MIN, the
    temperature of the column solution at every level (the reference plus
    the rise; composition as the chemistry left it, linear in ln p at the
    added levels), n_tot = p/(k T), rho = n_tot mu m_u at the level's mean
    molecular weight, and r re-integrated hydrostatically from the reference
    level on the new temperature."""
    d = res['domain']
    add = d['p_bar'][np.abs(res['dT']) >= PROFILE_LEVEL_DT_MIN]
    for p in add:
        cols = sch.insert_level(cols, float(p))
    p_bar = cols['p']
    mu = cols['rho']/(cols['n_tot']*sch.MAMU)
    T_ref_cols = res['T_reference_fn'](p_bar*sch.BAR)
    T = T_ref_cols + temperature_rise_at(res, p_bar*sch.BAR)
    cols['T'] = T
    cols['n_tot'] = p_bar*sch.BAR/(sch.KB*T)
    cols['rho'] = cols['n_tot']*mu*sch.MAMU
    cols['r'] = sch.hydrostatic_radius(p_bar, T, mu, args.mp, args.r_ref,
                                       args.p_ref)
    return cols


def cold_trap_note(p_bar, ratios_gas, ratios_all, p_match, sol,
                   climate=None):
    """The cold trap as the two elemental sums the design asks for.

    At the matching level, the oxygen the gas phase carries and the oxygen
    every carrier carries.  The first is what the wind inherits; the second
    is what the deep atmosphere holds; the ratio is the fraction of the
    oxygen the condensate keeps below.  The water mixing ratio at the
    tropopause -- the level the climate solve puts the trap at -- is quoted
    beside it, because that is the number the published estimates state.
    """
    p = np.asarray(p_bar, dtype=float)
    im = int(np.argmin(np.abs(p - p_match)))
    parts = []
    for el in ('O', 'C', 'N'):
        if el in ratios_gas and el in ratios_all:
            g, a = float(ratios_gas[el][im]), float(ratios_all[el][im])
            parts.append('%s/H gas %.4e vs all-carrier %.4e' % (el, g, a))
    note = 'cold trap at the match: ' + '; '.join(parts) if parts else ''
    q = sol.get('H2O')
    if q is not None and climate is not None:
        qa = np.asarray(q, dtype=float)
        it = int(np.argmin(np.abs(p - climate['P_trop_bar'])))
        note += ('; chemistry q_H2O at the tropopause (%.4g bar) %.4e, at '
                 'the match %.4e'
                 % (climate['P_trop_bar'], float(qa[it]), float(qa[im])))
    return note + '.'


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sch.add_common_arguments(ap)
    ap.add_argument('--tp-file', default=None,
                    help='prescribed T(p) [and K_zz(p)] column (milestone'
                         ' E2); mutually exclusive with --climate')
    ap.add_argument('--climate', action='store_true',
                    help='solve the radiative-convective column for T(p)'
                         ' instead of prescribing it (milestone E3)')
    ap.add_argument('--climate-p-deep', type=float, default=10.0,
                    help='deep boundary of the climate solve [bar]'
                         ' (default 10)')
    ap.add_argument('--climate-p-top', type=float, default=None,
                    help='top of the climate solve [dyn/cm^2]; default'
                         ' 1e-2 (1e-8 bar), above the matching level where'
                         ' the chemistry grid ends')
    ap.add_argument('--climate-layers', type=int, default=60,
                    help='number of climate layers (default 60)')
    ap.add_argument('--surface-albedo', type=float, default=0.0,
                    help='albedo of the deep boundary (default 0)')
    ap.add_argument('--climate-t-deep-guess', type=float, default=400.0,
                    help='start of the bracket of the deep-temperature'
                         ' fixed point (deep composition in chemical'
                         ' equilibrium) [K]; the root does not depend on it')
    ap.add_argument('--climate-t-trop-guess', type=float, default=120.0,
                    help='initial guess for the stratospheric temperature'
                         ' [K]; it is solved for, not imposed')
    ap.add_argument('--p-unit', default='dyn', choices=('dyn', 'bar'),
                    help='pressure unit of --tp-file (default dyn/cm^2)')
    ap.add_argument('--stellar-flux', required=True,
                    help='stellar spectrum, erg/cm^2/s/nm')
    ap.add_argument('--flux-at-planet', action='store_true',
                    help='the spectrum is already at the planet')
    ap.add_argument('--wavelength-unit', default='nm', choices=('nm', 'A'),
                    help='wavelength unit and flux-density denominator of'
                         ' --stellar-flux (default nm, i.e. erg/cm^2/s/nm)')
    ap.add_argument('--r-star', type=float, default=None,
                    help='stellar radius [R_sun], for the dilution')
    ap.add_argument('--a-orb', type=float, default=None,
                    help='orbital distance [AU], for the dilution')
    ap.add_argument('--atoms', default='H,He,N,O,C',
                    help='atom list of the Zahnle set (default H,He,N,O,C;'
                         ' add S for the gas-giant sulfur mechanism)')
    ap.add_argument('--mechanism', default=None,
                    help='use this mechanism file instead of a Zahnle set')
    ap.add_argument('--thermo', default=None)
    ap.add_argument('--data-dir', default=None)
    ap.add_argument('--abundances', default=None,
                    help='comma-separated El=value overrides of the Lodders'
                         ' 2009 vector, as number ratios to H')
    ap.add_argument('--zenith', type=float, default=None,
                    help='solar zenith angle of Photochem [deg] (default 48);'
                         ' under --face-geometry wind the convention sets it'
                         ' and stating it is refused')
    # DEFAULT `wind`: the column models take the stellar energy the wind's
    # input states (2D approximate method -> dayside_dilution), one
    # convention for every flux crossing the face. `per_model` keeps each
    # column model's own convention.
    ap.add_argument('--face-geometry', default='wind',
                    choices=FACE_GEOMETRIES,
                    help='how the column models take the stellar beam:'
                         ' per_model (Photochem 0.5 F at --zenith,'
                         ' the climate model F/4 over its zenith quadrature,'
                         ' each code\'s own convention) or wind (default; both take'
                         ' the wind\'s dayside dilution xi_w F per unit face'
                         ' area with the angular distribution it implies)')
    ap.add_argument('--wind-dilution', type=float, default=None,
                    help='the wind\'s dayside dilution xi_w (1, 0.5 or'
                         ' 0.25; EXHALE_setup.out "Dayside dilution") for'
                         ' --face-geometry wind')
    ap.add_argument('--toa', type=float, default=None,
                    help='REFUSED: the model top is the matching level'
                         ' (its top face at --p-match). Kept only to say so')
    ap.add_argument('--p-top-bar', type=float, default=None,
                    help='REFUSED, for the reason --toa is')
    ap.add_argument('--trial-flux-H2', type=float, default=None,
                    help='the part of --trial-flux-H that H2 carries through'
                         ' the top face [g/s of hydrogen nuclei, outward'
                         ' positive]; atomic H carries the rest, negative'
                         ' where the wind mixes atoms down. The wind states'
                         ' it (EXHALE_resolved.out base_face_H2_flux_cgs).'
                         ' Without it the whole elemental flux goes on H2'
                         ' (split by the top\'s own atomic share above 1'
                         ' per cent)')
    ap.add_argument('--boa-pressure-factor', type=float, default=5.0,
                    help='the photochemical model bottom sits this factor'
                         ' deeper in pressure than the deepest quench level'
                         ' (photochem default 5); state 1 on an atmosphere'
                         ' too cold to carry a quench level inside the'
                         ' column')
    ap.add_argument('--t-max', type=float, default=None,
                    help='truncate the column above this temperature')
    ap.add_argument('--blocks', type=int, default=400)
    ap.add_argument('--conducted-heat-flux', type=float, default=0.0,
                    help='heat flux the escape wind conducts into the column'
                         ' through the matching level [erg cm^-2 s^-1],'
                         ' positive into the column (EXHALE_resolved.out'
                         ' base_conductive_flux_cgs). Nonzero solves the'
                         ' column energy balance below p_match'
                         ' (conducted_heat_column_balance.py) and iterates'
                         ' it with the chemistry; 0 (default) writes the'
                         ' profile without it')
    ap.add_argument('--column-energy-treatment', default='nlte',
                    choices=('nlte', 'nlte_partial', 'lte', 'rotational'),
                    help='band emission of the column energy balance: nlte'
                         ' (default; vibrational bands with published V-T'
                         ' rates as two-level systems, the rest in LTE),'
                         ' lte (every band in LTE, an upper bound on the'
                         ' cooling), rotational (rotational bands only, a'
                         ' lower bound), nlte_partial (two-level also for a'
                         ' band whose rates miss a partner, a lower bound on'
                         ' its emission)')
    ap.add_argument('--nh3-h2-rate-analog',
                    default='cd4_bending',
                    choices=('cd4_bending', 'water_bending_overtone'),
                    help='the analog whose measured k(H2)/k(He) turns the'
                         ' measured NH3 nu2 + He de-excitation rate into the'
                         ' unmeasured NH3 nu2 + H2 one (treatment nlte only):'
                         ' cd4_bending (default; CD4 bending modes, Siddles'
                         ' et al. 1994 Tables 1-2, 27.5-49), or'
                         ' water_bending_overtone (H2O 2nu2 at 295 K, Zittel'
                         ' & Masturzo 1991 Table I over Hovis & Moore 1980'
                         ' Table III, 9.06, held constant in T)')
    ap.add_argument('--column-band-transfer', default='thin',
                    choices=('thin', 'cooling_to_space',
                             'cooling_to_space_exchange'),
                    help='how the band emission of the column energy balance'
                         ' leaves the column: thin (default; all of it;'
                         ' flagged, not refused, when more than 10 per cent'
                         ' of a bin\'s excess is absorbed above: a WARNING'
                         ' and thin_layer_limit_exceeded T in the header) or'
                         ' cooling_to_space (the'
                         ' upward half escapes through the column above with'
                         ' the k-distribution transmission, Dickinson 1972'
                         ' Eq. 33; refused when the neglected exchange inside'
                         ' the warmed layer moves T_match by more than'
                         ' --column-exchange-shift-max) or'
                         ' cooling_to_space_exchange (cooling_to_space with'
                         ' that exchange solved instead of neglected: no such'
                         ' bound)')
    ap.add_argument('--column-exchange-shift-max', type=float,
                    default=1.0,
                    help='validity bound of cooling_to_space [K] (default 1)')
    ap.add_argument('--top-pressure-passes', type=int,
                    default=TOP_PRESSURE_PASSES,
                    help='largest number of re-pinned solves of the model top'
                         ' (default %d); the agreement they must reach,'
                         ' TOP_PRESSURE_TOL, is not changed by it. A grid whose'
                         ' top converges by alternating passes needs more of'
                         ' them' % TOP_PRESSURE_PASSES)
    ap.add_argument('--column-energy-tolerance', type=float,
                    default=COLUMN_TEMPERATURE_TOL,
                    help='agreement of the matching-level temperature between'
                         ' successive column solutions [K] (default 1)')
    ap.add_argument('--column-energy-max-iterations', type=int, default=8)
    # DEFAULT off (user decision 2026-10-10, after the measurement of
    # md/eddy_heat_transport_column_20261010.md): with the chemistry's
    # K_zz as a heat diffusivity no reference climate exists above
    # K_zz ~ 1.4e6 cm^2 s^-1 (xi_w = 1; 1.8e5 at the per-model geometry):
    # the eddy cooling drives the stratosphere below 50 K, so the K_zz = 1e9
    # of the LHS 1140 b closures cannot be a turbulent heat diffusivity.
    # K_zz is therefore used as an EFFECTIVE TRACER-MIXING coefficient that
    # carries no heat. THE APPROXIMATION AND ITS SIZE: booked as a change
    # over the isothermal T0, the eddy flux would move T_match by -55 K
    # (K_zz = 1e9) to -2 K (1e5) (section 128); in a self-consistent
    # reference at K_zz = 1e5 it moves T_match by +3 K (section 130). Whether
    # turbulent dissipation heating (Hunten 1974, JGR 79, 2533) makes a
    # consistent treatment possible is under review; `on` keeps the
    # implementation for that review and for heat-diffusivity sensitivities.
    ap.add_argument('--eddy-heat-transport', default='off',
                    choices=('on', 'off'),
                    help='eddy heat flux rho c_p K_zz (dT/dz + g/c_p) with'
                         ' the K_zz of --kzz-const/--kzz-power: on'
                         ' solves the climate stratosphere in radiative +'
                         ' molecular + eddy energy balance (needs --climate)'
                         ' and adds the eddy flux to the column energy'
                         ' balance; off (default) keeps clima\'s isothermal'
                         ' stratosphere and molecular conduction only, K_zz'
                         ' a tracer-mixing coefficient')
    ap.add_argument('--workdir', default=None,
                    help='where the mechanism and flux files are written'
                         ' (default <run_dir>/photochem_work)')
    args = ap.parse_args()

    if args.p_match <= 0.0:
        sch.refuse('--p-match must be positive')
    # The model top is the matching level: above it the wind alone
    # describes the gas.
    if args.toa is not None or args.p_top_bar is not None:
        sch.refuse('--toa/--p-top-bar: the model top is the matching level '
                   '(its top face at --p-match %.3e bar); a top above it is '
                   'the overlap the matching surface removes. Remove '
                   'the option' % args.p_match)
    if args.trial_flux_H2 is not None and args.trial_flux_H == 0.0:
        sch.refuse('--trial-flux-H2 without an elemental --trial-flux-H')
    if args.climate == (args.tp_file is not None):
        sch.refuse('state exactly one of --tp-file (a prescribed column) and'
                   ' --climate (a solved one)')
    wdir = args.workdir or os.path.join(args.run_dir, 'photochem_work')
    os.makedirs(wdir, exist_ok=True)
    flux_file = os.path.join(wdir, 'flux_photochem.txt')
    if args.face_geometry == 'wind' and args.zenith is not None:
        sch.refuse('--zenith with --face-geometry wind: the convention sets'
                   ' Photochem\'s zenith angle')
    xi_wind = args.wind_dilution
    geometry = face_geometry(args.face_geometry, xi_wind,
                             48.0 if args.zenith is None else args.zenith,
                             CLIMATE_ZENITH_ANGLES)
    FACE_GEOMETRY.clear()
    FACE_GEOMETRY.update(geometry)
    args.zenith = geometry['photochem_zenith_deg']
    climate_flux_file = flux_file
    dilution = write_flux_file(args, flux_file)
    mech, thermo, data_dir = mechanism_files(args, wdir)
    abundances = parse_abundances(args.abundances)

    column, climate = None, None
    EDDY_HEAT.clear()
    if args.eddy_heat_transport == 'on':
        if not args.climate and args.conducted_heat_flux != 0.0:
            sch.refuse('--eddy-heat-transport on needs the reference column'
                       ' solved with the same eddy heat flux (--climate);'
                       ' a prescribed --tp-file column is not a solution of'
                       ' it. State --eddy-heat-transport off with --tp-file')
        if args.climate:
            if sch.eddy_diffusion_coefficient(np.array([1.0]), args) is None:
                sch.refuse('--eddy-heat-transport on needs the eddy'
                           ' coefficient: state --kzz-const, or --kzz-power'
                           ' with --kzz-ref')
            EDDY_HEAT['kzz'] = (
                lambda p_bar: sch.eddy_diffusion_coefficient(
                    np.asarray(p_bar, dtype=float), args))
            EDDY_HEAT['radius_cm'] = (
                lambda p_bar, T, m_g: sch.hydrostatic_radius(
                    p_bar, T, m_g/sch.MAMU, args.mp, args.r_ref, args.p_ref)
                * sch.RJ)
    if args.climate:
        # The climate column must reach above the matching level, where the
        # chemistry grid ends; by default it reaches 1e-2 dyn/cm^2 (1e-8
        # bar), and only its part below the match is used.
        p_top_dyn = (args.climate_p_top if args.climate_p_top
                     else CLIMATE_P_TOP_DYN)
        climate = rcc.solve_radiative_convective_column(
            wdir, climate_flux_file, args.mp*sch.MJ, args.r_ref*sch.RJ,
            abundances,
            args.climate_p_deep, p_top_dyn=p_top_dyn,
            surface_albedo=args.surface_albedo,
            nlayers=args.climate_layers,
            t_deep_guess=args.climate_t_deep_guess,
            t_trop_guess=args.climate_t_trop_guess,
            zenith_u=geometry['climate_zenith_u'],
            photon_scale_factor=geometry['climate_photon_scale'],
            eddy_kzz=EDDY_HEAT.get('kzz'),
            radius_cm=EDDY_HEAT.get('radius_cm'),
            p_anchor_dyn=args.p_match*sch.BAR)
        st = climate.get('stratosphere')
        if st is not None:
            EDDY_HEAT['radiative_flux'] = rcc.net_flux_function(
                st['P_edges'], st['F_edges'])
            EDDY_HEAT['levels_dyn'] = np.asarray(st['P_dyn'], dtype=float)
        print('  ' + rcc.summary(climate))
        sys.stdout.flush()
        column = (climate['P_dyn'], climate['T'], None)

    conducted = args.conducted_heat_flux != 0.0
    if conducted:
        # The reference column the conducted heat is added to, on levels
        # fine enough below the match for the chemistry to see the response.
        if column is None:
            P0, T0c, K0 = read_tp_file(args.tp_file, args.p_unit)
        else:
            P0, T0c, K0 = column
        column = column_with_inserted_levels(P0, T0c, K0, args.p_match)

    pc, reached, truncation, toa, imposed, r_top_imposed = \
        photochemical_steady_state(args, wdir, flux_file, mech, thermo,
                                   data_dir, abundances, column=column)

    column_heat = None
    if conducted:
        p_floor = climate['P_trop_bar'] if climate is not None else None
        column_heat = conducted_heat_iteration(pc, args, column, p_floor)
        write_column_heat_record(args, column_heat)

    # The hydrogen carrier split of the escape boundary condition, made
    # self-consistent: the split is decided on the composition the
    # imposition sees, and photochemistry dissociates H2 as it runs, so the
    # converged top can carry more atomic H than the boundary condition put
    # the flux on.  Where the two shares differ by more than
    # CARRIER_SHARE_TOL the measured split is re-imposed and the chemistry
    # re-solved (under the conducted heat, the whole column iteration).
    carrier_rounds = 0
    if 'H_atomic_share' in imposed:
        for _ in range(CARRIER_ROUNDS):
            f_now = hydrogen_carrier_split(pc, species_element_counts(pc))
            used = LAST_IMPOSITION['H_atomic_share_used']
            if abs(f_now - used) <= CARRIER_SHARE_TOL:
                break
            carrier_rounds += 1
            # The converged share falls as the share put on atomic H rises
            # (MEASURED at deep He/H 30 on the GJ 699 spectrum: 0 -> 0.098,
            # 0.098 -> 0.025, 0.025 -> 0.079, 0.079 -> 0.039, slope about
            # -0.75), so substituting the converged share alternates about
            # the fixed point; once a split is imposed the next one is the
            # mean of the imposed and the converged share (slope about
            # +0.13 there).
            CARRIER_SHARE_STATED['H'] = (f_now if used == 0.0
                                         else 0.5*(used + f_now))
            print('  carrier round %d: the converged top carries %.4e of its'
                  ' hydrogen nuclei as atomic H, the boundary condition %.4e;'
                  ' the split %.4e is imposed and the chemistry re-solved'
                  % (carrier_rounds, f_now, used, CARRIER_SHARE_STATED['H']),
                  flush=True)
            if conducted:
                column_heat = conducted_heat_iteration(pc, args, column,
                                                       p_floor)
                write_column_heat_record(args, column_heat)
            else:
                ref = column if column is not None else \
                    read_tp_file(args.tp_file, args.p_unit)
                _, Pc, Tc, Kc = reference_column_on_chemistry_grid(pc, args,
                                                                   ref)
                chemistry_on_new_temperature(
                    pc, args, species_element_counts(pc), Pc, Tc, Kc,
                    chemistry_temperature_profile(args, ref))

    sol = pc.return_atmosphere()
    counts = species_element_counts(pc)
    measured, r_top = measured_elemental_escape_flux(pc, counts)
    flux_note = ('measured elemental escape at the model top (r_top = '
                 '%.6e cm): H %.6e g/s, He %.6e g/s; imposed H %.6e, He '
                 '%.6e' % (r_top, measured['H'], measured['He'],
                           args.trial_flux_H, args.trial_flux_He))
    print('  ' + flux_note)
    if 'H_atomic_share' in imposed:
        # The guard after the carrier rounds above: if the share the
        # boundary condition carries still differs from the converged top's
        # by more than CARRIER_SHARE_TOL, more than 1% of the hydrogen flux
        # was drained from the wrong carrier and the solution is not the one
        # that was asked for.
        f_final = hydrogen_carrier_split(pc, counts)
        used = LAST_IMPOSITION['H_atomic_share_used']
        print('  atomic-H share of the hydrogen nuclei at the model top: '
              '%.3e at the imposition, %.3e converged'
              % (imposed['H_atomic_share'], f_final))
        if abs(f_final - used) > CARRIER_SHARE_TOL:
            sch.refuse('the converged model top carries %.3f of its hydrogen'
                       ' nuclei as atomic H, but the boundary condition put'
                       ' %.3f of the elemental flux on atomic H after %d'
                       ' carrier rounds'
                       % (f_final, used, carrier_rounds))
        flux_note += ('; atomic-H share of the hydrogen nuclei at the top '
                      '%.3e imposed, %.3e converged'
                      % (imposed['H_atomic_share'], f_final))
        if carrier_rounds:
            flux_note += ('; the measured carrier split was re-imposed in %d'
                          ' round(s), atomic share of the boundary condition'
                          ' %.4e' % (carrier_rounds, used))
    require_flux_imposed_at_final_top(imposed, r_top)
    nparticles = int(pc.dat.np)
    particles = set(list(pc.dat.species_names)[:nparticles])

    p_dyn = np.asarray(sol['pressure'], dtype=float)
    T = np.asarray(sol['temperature'], dtype=float)
    p_bar = p_dyn/sch.BAR

    # Two elemental sums, always: over the gas phase, which is what the wind
    # inherits and what the file's X_<El> columns carry by default, and over
    # every carrier including the condensates.  Their difference IS the cold
    # trap -- the nuclei the condensed phase holds back at each level -- so
    # both are computed and both are reported (design section 5, item 4).
    all_carriers = {sp: np.asarray(y, dtype=float) for sp, y in sol.items()
                    if sp in counts}
    gas_only = {sp: y for sp, y in all_carriers.items()
                if sp not in particles}
    ratios_gas, _ = sch.element_ratios(gas_only, counts)
    ratios_all, _ = sch.element_ratios(all_carriers, counts)
    mixing = all_carriers if args.count_condensates else gas_only
    ratios = ratios_all if args.count_condensates else ratios_gas
    mu = sch.mean_molecular_weight(mixing, counts)

    Kzz = sch.eddy_diffusion_coefficient(p_bar, args)
    if Kzz is None:
        if 'Kzz' in sol:
            Kzz = np.asarray(sol['Kzz'], dtype=float)
        else:
            sch.refuse('the solution carries no K_zz: state --kzz-const, '
                       'or --kzz-power with --kzz-ref')

    n_tot = p_dyn/(sch.KB*T)
    rho = n_tot*mu*sch.MAMU
    r_RJ = sch.hydrostatic_radius(p_bar, T, mu, args.mp, args.r_ref,
                                  args.p_ref)

    cols = {'p': p_bar, 'r': r_RJ, 'T': T, 'n_tot': n_tot, 'rho': rho,
            'Kzz': Kzz,
            'q_H2': np.asarray(sol.get('H2', np.zeros_like(p_bar))),
            'q_H': np.asarray(sol.get('H', np.zeros_like(p_bar)))}
    order = ['p', 'r', 'T', 'n_tot', 'rho', 'Kzz', 'q_H2', 'q_H']
    for el in sch.PROFILE_ELEMENTS:
        if el in ratios:
            cols['X_' + el] = ratios[el]
            order.append('X_' + el)
    cols['F_H'] = np.full_like(p_bar, args.trial_flux_H)
    cols['F_He'] = np.full_like(p_bar, args.trial_flux_He)
    order += ['F_H', 'F_He']
    for sp in sch.DIAGNOSTIC_MOLECULES:
        if sp in sol:
            cols['q_' + sp] = np.asarray(sol[sp], dtype=float)
            order.append('q_' + sp)

    p_deep = args.p_deep if args.p_deep else p_bar.max()
    # The table ends at the matching level, the top face of the chemistry
    # grid, whose state is continued from the top cell centres
    # (level_at_top_face).
    cols = level_at_top_face(sch.sort_deep_to_shallow(cols), args.p_match)
    cols = sch.clip_table(cols, p_deep, args.p_match)
    if column_heat is not None:
        cols = profile_with_column_temperature(cols, column_heat, args)

    fp = dict(
        source_code='photochem',
        mechanism=os.path.basename(mech), mechanism_sha=sch.sha256_file(mech),
        thermo=os.path.basename(thermo), thermo_sha=sch.sha256_file(thermo),
        stellar_flux=os.path.basename(args.stellar_flux),
        stellar_flux_sha=sch.sha256_file(args.stellar_flux),
        dilution=dilution,
        tp_file=(os.path.basename(args.tp_file) if args.tp_file
                 else 'climate'),
        tp_sha=sch.sha256_file(args.tp_file),
        boa_pressure_factor=args.boa_pressure_factor,
        mp_MJ=args.mp, r_ref_RJ=args.r_ref, p_ref_bar=args.p_ref,
        zenith_deg=args.zenith, toa_dyn=toa, t_max=args.t_max,
        kzz_const=args.kzz_const, count_condensates=args.count_condensates,
        abundances=sorted(abundances.items()),
        p_match_bar=args.p_match,
        trial_flux_H=args.trial_flux_H, trial_flux_He=args.trial_flux_He,
        model_top='top face at p_match',
    )
    if args.trial_flux_H2 is not None:
        fp.update(trial_flux_H2=args.trial_flux_H2)
    # Only off the defaults, so that a run without them keeps its
    # solution_id: the face convention names a different stellar energy
    # per unit area.
    if args.face_geometry != 'per_model':
        fp.update(face_geometry=args.face_geometry, wind_dilution=xi_wind)
    if args.kzz_power is not None:
        # Only when it is in use: a fingerprint entry that is always present
        # would move the solution_id of every constant-K_zz run.
        fp.update(kzz_power=args.kzz_power, kzz_ref=args.kzz_ref,
                  kzz_ref_bar=args.kzz_ref_bar)
    if climate is not None:
        # The climate solve is part of the solution the id names, so its
        # configuration and its two answers enter the fingerprint.
        fp.update(climate_p_deep_bar=args.climate_p_deep,
                  climate_layers=args.climate_layers,
                  surface_albedo=args.surface_albedo,
                  climate_T_deep=climate['T_deep'],
                  climate_P_trop_bar=climate['P_trop_bar'],
                  climate_T_trop=climate['T_trop'])
        if climate.get('stratosphere') is not None:
            # only when on, so that an `off` run keeps its solution_id
            st = climate['stratosphere']
            fp.update(eddy_heat_transport='on',
                      stratosphere_T_anchor=float(st['T'][-1]),
                      stratosphere_T_min=st['T_min'])
    if column_heat is not None:
        # Only with a conducted flux, so that a run without one keeps its
        # solution_id: the flux, the treatment of the bands and the
        # temperature they gave at the match name a different solution.
        fp.update(conducted_heat_flux=args.conducted_heat_flux,
                  column_energy_treatment=args.column_energy_treatment,
                  column_T_match=column_heat['T_match'])
        if args.column_band_transfer != 'thin':
            # only off the default, so a thin solution keeps its id
            fp.update(column_band_transfer=args.column_band_transfer)
    try:
        import photochem
        version = 'photochem %s' % photochem.__version__
    except Exception:
        version = 'photochem (version unknown)'
    fp['source_version'] = version

    notes = 'F_H and F_He columns are the STATED trial fluxes; the header ' \
            'measured_flux_* are what the solution carries across its own ' \
            'top. ' + flux_note + '. '
    if climate is None:
        notes += 'T(p) prescribed, no climate step. '
    else:
        notes += rcc.summary(climate) + '. '
    notes += cold_trap_note(p_bar, ratios_gas, ratios_all, args.p_match,
                            sol, climate)
    notes += ' ' + truncation
    if column_heat is not None:
        notes += ' ' + column_heat_note(args, column_heat)
    if args.face_geometry != 'per_model':
        # The notes line is cut at the 1000 characters the EXHALE reader
        # holds (lower_profile_schema.write_profile), and the schema writes
        # only the header keys it knows: the record of the face geometry
        # goes to its own file beside the profile.
        write_face_radiation_record(args, geometry, xi_wind)
    p_face, _ = top_face_pressure(pc)
    notes += (' Model top: the top face at %.6e bar (%+.2e in ln p from the'
              ' matching level), the state at the match continued from the'
              ' two top cell centres; above the match the wind alone'
              ' describes the gas.' % (p_face/sch.BAR,
                                       np.log(p_face/(args.p_match*sch.BAR))))
    header = dict(
        solution_id=sch.solution_fingerprint(fp),
        source_code='photochem',
        source_version=version + ', photochem_to_lower_profile.py',
        mechanism='%s sha256 %s' % (os.path.basename(mech),
                                    sch.sha256_file(mech)[:16]),
        stellar_flux='%s sha256 %s dilution %.6E'
                     % (os.path.basename(args.stellar_flux),
                        sch.sha256_file(args.stellar_flux)[:16], dilution),
        p_match_bar=args.p_match, p_top_bar=0.0, p_deep_bar=0.0,
        trial_flux_H=args.trial_flux_H, trial_flux_He=args.trial_flux_He,
        measured_flux_H=measured['H'], measured_flux_He=measured['He'],
        iteration=args.iteration, reached_steady_state=bool(reached),
        notes=notes.strip())
    if steady_state_retries(pc):
        # the passes that reached their steady state only on the single
        # retry from the stalled state (steady_state_with_stall_retry)
        header['steady_state_retry'] = '; '.join(
            'rtol=%.0e at pass %s' % (STALL_RETRY_RTOL, label)
            for label in steady_state_retries(pc))
    if args.trial_flux_H2 is not None:
        _, top_fluxes = pc.gas_fluxes()
        header.update(
            trial_flux_H2=args.trial_flux_H2,
            measured_flux_H2=float(top_fluxes.get('H2', 0.0))*2.0
            * 4.0*np.pi*r_top*r_top*sch.ATOMIC_WEIGHT['H']*sch.MAMU)
        if column_heat is not None and 'recombination' in column_heat:
            header['recombination_heat_flux_cgs'] = \
                column_heat['recombination']['F_chem']
    if column_heat is not None:
        header.update(conducted_heat_flux_cgs=args.conducted_heat_flux,
                      column_energy_treatment=args.column_energy_treatment,
                      column_T_match_K=column_heat['T_match'])
        if args.column_energy_treatment == 'nlte':
            header['nh3_h2_rate_analog'] = args.nh3_h2_rate_analog
        if args.column_band_transfer != 'thin':
            header['column_band_transfer'] = args.column_band_transfer
        else:
            tl = column_heat['thin_layer']
            header['thin_layer_limit_exceeded'] = 'T' if tl['exceeded'] \
                else 'F'
            header['thin_layer_max_absorbed_fraction'] = tl['max_absorbed']
            header['thin_layer_max_absorbed_bin'] = '%s_%s_%.0f-%.0f' % (
                tl['species'], tl['band'].replace(' ', '_'),
                tl['bin_cm'][0], tl['bin_cm'][1])

    deepest = {el: float(v[0]) for el, v in ratios.items()}
    want = {el: abundances[el]/abundances['H'] for el in abundances
            if el != 'H'}
    sch.write_handoff(args, cols, header, order, element_input=want,
               element_measured=deepest)


if __name__ == '__main__':
    main()
