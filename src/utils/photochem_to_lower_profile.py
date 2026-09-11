#!/usr/bin/env python3
"""Photochem -> EXHALE lower-atmosphere profile (milestones E2 and E3).

Solves a gas-giant photochemistry model to steady state on a T(p), K_zz(p)
column -- prescribed (`--tp-file`, milestone E2) or solved for by the
radiative-convective climate model (`--climate`, milestone E3,
`radiative_convective_column.py`) -- and writes the handoff of
`docs/phase_e_flux_closure_design.md` section 2: one
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
trap (`docs/phase_e_flux_closure_design.md` section 5).

THE GRID THE COLUMN IS WRITTEN ON is the solution's own, not the one the
stepper happened to stop on.  Photochem accepts a steady state on any grid
whose top is within a factor 3 of the requested one, so an unpinned run comes
back on one of several grids and the converged column differs with the grid.
`steady_state_at_stated_model_top` re-pins the grid at `--toa` and re-converges
the chemistry on it, iterating to the fixed point, and refuses if that will not
settle (`docs/lower_profile_deep_boundary_sensitivity.md` section 8).

The driving pattern is the one `docs/p1_matched_comparison.py` validated on
2026-08-26 and is not re-derived here: dilution-only stellar flux, the Zahnle
set restricted to an atom list, the climate grid cut so the photochemical
grid can sit above it, elemental abundances mapped BY NAME, and the robust
stepper in blocks.

Environment: run this under a Python that has Photochem 0.9.0 built from
the source in `photochem/` -- upstream with the elemental-closure
corrections of `README_photochem.md` section 4, which says what to download
and how to build it.  Photochem 0.8.4 is a different code with a different
elemental closure; every profile this writes records which one made it in
its `# source_version` header line.

Example (HD 209458 b, the P1 configuration):

  python3 src/utils/photochem_to_lower_profile.py <run_dir> \\
      --mp 0.720 --r-ref 1.36 \\
      --tp-file vulcan_work/hd209_vulcan/atm/atm_HD209_Kzz.txt \\
      --stellar-flux vulcan_work/hd209_vulcan/atm/stellar_flux/Gueymard_solar.txt \\
      --r-star 1.155 --a-orb 0.0480 --p-match 1e-6 --toa 1e-2
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


def write_flux_file(args, path):
    """Stellar flux at the planet, in mW/m^2/nm, on a wavelength grid in nm.

    Photochem and clima both read this one file.  A spectrum tabulated per
    nm in erg/cm^2/s/nm needs no numerical conversion (1 erg/cm^2/s/nm =
    1 mW/m^2/nm); one tabulated per angstrom does, in both columns at once:
    the grid divides by 10 and the flux density multiplies by 10, so that
    the integral over wavelength is unchanged.  Only the dilution
    (R_star/a)^2 is applied on top, and only when the spectrum is at the
    stellar surface.
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
    with open(path, 'w') as fh:
        fh.write('Wavelength (nm)               Solar flux (mW/m^2/nm)\n')
        for w, y in zip(wl, f*dilution):
            fh.write('%.6e                  %.6e\n' % (w, y))
    return dilution


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


def impose_elemental_escape_flux(pc, counts, phi_H_gs, phi_He_gs):
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

    HYDROGEN CARRIER.  The whole elemental hydrogen flux is put on H2, at
    half the nuclei flux because H2 carries two nuclei.  This is the
    dominant-carrier approximation, and it is applied only when the measured
    atomic-H share of the hydrogen nuclei at the model top is below 1%; the
    share is printed with the imposition.  Above 1% the flux is split
    between H and H2 in proportion to the nuclei each carries there, so the
    approximation never silently covers a dissociated top.

    A zero flux sets no boundary condition at all, so the closed-top
    solution is reproduced bit for bit.
    """
    r_top = model_top_radius(pc)
    area = 4.0*np.pi*r_top*r_top
    imposed = {}
    if phi_H_gs != 0.0:
        phi_H = phi_H_gs/(sch.ATOMIC_WEIGHT['H']*sch.MAMU)/area
        f_atomic = hydrogen_carrier_split(pc, counts)
        if f_atomic > 0.01:
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
        print('  imposed elemental H escape: %.6e g/s -> %.6e nuclei/cm^2/s '
              'at r_top = %.6e cm, %s' % (phi_H_gs, phi_H, r_top, how))
    if phi_He_gs != 0.0:
        phi_He = phi_He_gs/(sch.ATOMIC_WEIGHT['He']*sch.MAMU)/area
        pc.set_upper_bc('He', bc_type='flux', flux=phi_He)
        imposed['He'] = phi_He
        print('  imposed elemental He escape: %.6e g/s -> %.6e nuclei/cm^2/s'
              % (phi_He_gs, phi_He))
    return imposed, r_top


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
    toa = args.toa

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
        # 3*TOA_pressure_avg < P_climate_top (trap 1).  Set the top AND cut
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
            pc.gdat.TOA_pressure_avg = toa_try
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
        except PhotoException as exc:
            last_error = str(exc)
            print('  photochem refused the column (%s); truncating' % exc)
            continue
        except Exception as exc:
            if 'BOA in photochemical model' not in str(exc):
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

    # The closure's boundary condition, imposed before the column is stepped
    # so that the steady state it reaches is the state of a losing column.
    imposed, r_top_imposed = impose_elemental_escape_flux(
        pc, species_element_counts(pc), args.trial_flux_H, args.trial_flux_He)

    pc.initialize_robust_stepper(pc.wrk.usol)
    give_up = reached = False
    for n in range(args.blocks):
        for _ in range(100):
            give_up, reached = pc.robust_step()
            if give_up or reached:
                break
        print('  block %d: t = %.3e s%s%s'
              % (n + 1, pc.wrk.tn, '  STEADY' if reached else '',
                 '  GAVE UP' if give_up else ''), flush=True)
        if give_up or reached:
            break

    if give_up or not reached:
        sch.refuse('the chemistry did not reach steady state '
                   '(reached_steady_state = %s, gave_up = %s). A non-steady '
                   'solution has no elemental flux to hand over'
                   % (reached, give_up))

    steady_state_at_stated_model_top(pc, toa, args.blocks)
    return pc, reached, truncation, toa, imposed, r_top_imposed


# The grid is the solution's own: two successive re-pinned solves must agree
# on the model top to this relative tolerance.  It is 2500 times below the
# 2.5e-02 spread between the two grids an unpinned run stops on, so what is
# left of the grid dependence is 2e-06 in the elemental O/H at the matching
# level, against the 5.2e-03 it was.
TOP_PRESSURE_TOL = 1.0e-5
TOP_PRESSURE_PASSES = 6


def steady_state_at_stated_model_top(pc, toa, blocks,
                                     tol=TOP_PRESSURE_TOL,
                                     passes=TOP_PRESSURE_PASSES):
    """Re-solve the steady state on the grid the stated model top defines.

    Photochem lays its levels out uniformly in ALTITUDE between a pinned
    bottom and `top_atmos` (`vertical_grid`, `photochem/src/photochem_eqns.f90`),
    so a single number fixes the pressure of every level.  `robust_step`
    (`photochem/photochem/extensions/gasgiants.py`) re-pins `top_atmos` to the
    requested pressure only every `freq_update_TOA = 1000` internal steps, and
    its exit test accepts any top pressure within a FACTOR 3 of it.  A run is
    therefore declared steady on whatever grid the last re-pinning happened to
    leave: measured on LHS 1140 b, two runs of the same problem come back on
    tops 2.5 % apart, and the converged column differs with the grid -- 2.5 %
    in pressure at fixed level index and 5.2e-03 in the elemental O/H handed
    over at the matching level
    (`docs/lower_profile_deep_boundary_sensitivity.md`).

    The grid is therefore re-pinned here, from the converged solution, and the
    chemistry re-converged on it.  The two steps are ITERATED: re-pinning
    perturbs the solution, and the perturbed solution moves `top_atmos` again
    through the hydrostatic integration.  Their common fixed point --
    `top_atmos` is the altitude at which the steady state computed on the grid
    it defines reaches the stated pressure -- is a property of the problem, so
    the column that comes out no longer records which phase of the re-pinning
    cycle the integration happened to exit on.  The test is that two
    successive passes agree on the model top; the offset from the stated
    pressure is not the test, because the reported top is the top CELL CENTRE
    while `top_atmos` is the domain edge, and the half cell between them is
    about 9 % in pressure at this resolution.

    FAILURE IS A REFUSAL, not a fall-back to the unpinned grid.  A column the
    chemistry could not re-converge on, or one whose grid will not settle, is
    not a solution of the stated problem; writing it anyway would put back
    exactly the grid dependence this removes, and put it back invisibly, in a
    file that carries no record of which grid it stopped on.  A
    `PhotoException` out of `update_vertical_grid` is reported the same way,
    with the library's own message, because the adapter has nothing to add to
    it and nothing to write without it.
    """
    from photochem import PhotoException

    previous = None
    for attempt in range(1, passes + 1):
        try:
            pc.update_vertical_grid(TOA_pressure=toa)
        except PhotoException as exc:
            sch.refuse('the model top could not be pinned at the stated '
                       '%.3e dyn/cm^2: %s' % (toa, exc))
        pc.initialize_robust_stepper(pc.wrk.usol)
        give_up = reached = False
        for _ in range(blocks):
            for _ in range(100):
                give_up, reached = pc.robust_step()
                if give_up or reached:
                    break
            if give_up or reached:
                break
        p_top = float(pc.wrk.pressure_hydro[-1])
        moved = float('nan') if previous is None else p_top/previous - 1.0
        print('  pinned pass %d: t = %.3e s, top cell %.6e dyn/cm2, moved '
              '%s from the previous pass%s%s'
              % (attempt, pc.wrk.tn, p_top,
                 'n/a' if previous is None else '%+.3e' % moved,
                 '  STEADY' if reached else '',
                 '  GAVE UP' if give_up else ''), flush=True)
        if give_up or not reached:
            sch.refuse('the chemistry did not re-reach steady state on the '
                       'grid pinned at the stated model top %.3e dyn/cm^2 '
                       '(pass %d of %d, reached_steady_state = %s, gave_up = '
                       '%s)' % (toa, attempt, passes, reached, give_up))
        if previous is not None and abs(moved) <= tol:
            return p_top, attempt
        previous = p_top
    sch.refuse('the grid will not settle: after %d re-pinned solves the model '
               'top still moves by %+.3e between passes, against a tolerance '
               'of %.1e. The column would be written on a grid that is not '
               'the solution\'s own' % (passes, moved, tol))


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
                    help='top of the climate solve [dyn/cm^2]; default is'
                         ' the photochemical --toa, so the climate column'
                         ' always reaches above the chemistry grid')
    ap.add_argument('--climate-layers', type=int, default=60,
                    help='number of climate layers (default 60)')
    ap.add_argument('--surface-albedo', type=float, default=0.0,
                    help='albedo of the deep boundary (default 0)')
    ap.add_argument('--climate-t-deep-guess', type=float, default=400.0,
                    help='initial guess for the deep temperature [K]')
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
    ap.add_argument('--zenith', type=float, default=48.0,
                    help='solar zenith angle [deg] (default 48)')
    ap.add_argument('--toa', type=float, default=None,
                    help='top-of-atmosphere pressure [dyn/cm^2] (default'
                         ' 1e-2, i.e. 1e-8 bar); --p-top-bar states the same'
                         ' level in bar')
    ap.add_argument('--p-top-bar', type=float, default=None,
                    help='the shallowest pressure the profile carries [bar],'
                         ' stated directly: it IS the model top, and it sets'
                         ' --toa to p_top_bar x 1e6 dyn/cm^2 (the two state'
                         ' the same level in different units, so only one of'
                         ' them may be given). The model settles within a'
                         ' factor of about three of the stated top -- that'
                         ' is the window photochem holds it in -- so the'
                         ' written table stops near, not exactly at, this'
                         ' pressure; --p-top clips the written table and'
                         ' does not move the model top. The column handed'
                         ' to the chemistry is cut at 3.05 x TOA, so a top'
                         ' the supplied T(p) or climate column cannot reach'
                         ' with at least ten levels above the cut is'
                         ' refused')
    ap.add_argument('--boa-pressure-factor', type=float, default=5.0,
                    help='the photochemical model bottom sits this factor'
                         ' deeper in pressure than the deepest quench level'
                         ' (photochem default 5); state 1 on an atmosphere'
                         ' too cold to carry a quench level inside the'
                         ' column')
    ap.add_argument('--t-max', type=float, default=None,
                    help='truncate the column above this temperature')
    ap.add_argument('--blocks', type=int, default=400)
    ap.add_argument('--workdir', default=None,
                    help='where the mechanism and flux files are written'
                         ' (default <run_dir>/photochem_work)')
    args = ap.parse_args()

    if args.p_match <= 0.0:
        sch.refuse('--p-match must be positive')
    # The model top: one level, stated either in bar or in dyn/cm^2.
    if args.p_top_bar is not None:
        if args.toa is not None:
            sch.refuse('--p-top-bar and --toa state the same level in '
                       'different units: give one of them')
        if args.p_top_bar <= 0.0:
            sch.refuse('--p-top-bar must be positive')
        if args.p_top_bar >= args.p_match:
            sch.refuse('--p-top-bar %.3e bar is at or below the matching '
                       'pressure %.3e bar: the two models would never '
                       'overlap' % (args.p_top_bar, args.p_match))
        args.toa = args.p_top_bar*sch.BAR
    elif args.toa is None:
        args.toa = 1.0e-2
    if args.climate == (args.tp_file is not None):
        sch.refuse('state exactly one of --tp-file (a prescribed column) and'
                   ' --climate (a solved one)')
    wdir = args.workdir or os.path.join(args.run_dir, 'photochem_work')
    os.makedirs(wdir, exist_ok=True)
    flux_file = os.path.join(wdir, 'flux_photochem.txt')
    dilution = write_flux_file(args, flux_file)
    mech, thermo, data_dir = mechanism_files(args, wdir)
    abundances = parse_abundances(args.abundances)

    column, climate = None, None
    if args.climate:
        # The climate column must reach above the chemistry grid, whose top
        # the adapter cuts at 3.05 x TOA; solving it to the TOA itself leaves
        # that cut room to work with, whatever TOA the caller asked for.
        p_top_dyn = (args.climate_p_top if args.climate_p_top
                     else args.toa)
        climate = rcc.solve_radiative_convective_column(
            wdir, flux_file, args.mp*sch.MJ, args.r_ref*sch.RJ, abundances,
            args.climate_p_deep, p_top_dyn=p_top_dyn,
            surface_albedo=args.surface_albedo,
            nlayers=args.climate_layers,
            t_deep_guess=args.climate_t_deep_guess,
            t_trop_guess=args.climate_t_trop_guess)
        print('  ' + rcc.summary(climate))
        sys.stdout.flush()
        column = (climate['P_dyn'], climate['T'], None)

    pc, reached, truncation, toa, imposed, r_top_imposed = \
        photochemical_steady_state(args, wdir, flux_file, mech, thermo,
                                   data_dir, abundances, column=column)

    sol = pc.return_atmosphere()
    counts = species_element_counts(pc)
    measured, r_top = measured_elemental_escape_flux(pc, counts)
    flux_note = ('measured elemental escape at the model top (r_top = '
                 '%.6e cm): H %.6e g/s, He %.6e g/s; imposed H %.6e, He '
                 '%.6e' % (r_top, measured['H'], measured['He'],
                           args.trial_flux_H, args.trial_flux_He))
    print('  ' + flux_note)
    if 'H_atomic_share' in imposed:
        # The carrier split had to be decided before the column was stepped,
        # so it was measured on the initial composition.  Photochemistry
        # dissociates as it runs; if the CONVERGED top carries more than 1%
        # of its hydrogen nuclei as atomic H, the flux was drained from the
        # wrong carrier and the solution is not the one that was asked for.
        f_final = hydrogen_carrier_split(pc, counts)
        print('  atomic-H share of the hydrogen nuclei at the model top: '
              '%.3e at the imposition, %.3e converged'
              % (imposed['H_atomic_share'], f_final))
        if f_final > 0.01 >= imposed['H_atomic_share']:
            sch.refuse('the converged model top carries %.3f of its hydrogen'
                       ' nuclei as atomic H, but the elemental flux was'
                       ' imposed entirely on H2 because the initial'
                       ' composition carried %.3e. The dominant-carrier'
                       ' approximation does not hold on this solution'
                       % (f_final, imposed['H_atomic_share']))
        flux_note += ('; atomic-H share of the hydrogen nuclei at the top '
                      '%.3e imposed, %.3e converged'
                      % (imposed['H_atomic_share'], f_final))
    if imposed and abs(r_top/r_top_imposed - 1.0) > 5.0e-3:
        # The imposed number is a flux DENSITY; the grid the solver settles
        # on need not be the grid it started from, so the mass flux the
        # boundary condition represents moves with the top radius.  Say so
        # rather than let the driver read the difference as a residual.
        flux_note += ('; the model top moved from %.6e to %.6e cm while '
                      'solving, so the imposed mass flux moved with its '
                      'area' % (r_top_imposed, r_top))
        print('  NOTE: model top moved %.6e -> %.6e cm'
              % (r_top_imposed, r_top))
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
    p_top = args.p_top if args.p_top else p_bar.min()
    cols = sch.clip_table(sch.sort_deep_to_shallow(cols), p_deep, p_top)

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
    )
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

    deepest = {el: float(v[0]) for el, v in ratios.items()}
    want = {el: abundances[el]/abundances['H'] for el in abundances
            if el != 'H'}
    sch.write_handoff(args, cols, header, order, element_input=want,
               element_measured=deepest)


if __name__ == '__main__':
    main()
