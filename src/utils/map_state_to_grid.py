#!/usr/bin/env python3
"""Map an EXHALE state (Hydro_ioniz.txt + Ion_species.txt) onto the cell
centers of another grid, so that a state written on an earlier grid can be
loaded (`Load IC? True`) by a run whose grid construction has since changed
(load_IC refuses a state whose centers differ from the run's beyond 1e-10).

Usage:
  map_state_to_grid.py <src_dir> <target_grid_file> <out_dir> [--ic]
                       [--extrapolate-beyond <r/R_p>]
                       [--reservoir <El>/H <value> ...]

<src_dir>          holds Hydro_ioniz.txt and Ion_species.txt (or the *_IC.txt pair)
<target_grid_file> any Hydro_ioniz*.txt written by a run on the target grid;
                   only its first column (r) is read, ghost rows included
<out_dir>          receives Hydro_ioniz.txt and Ion_species.txt (with --ic,
                   Hydro_ioniz_IC.txt and Ion_species_IC.txt)

THE OUTPUT IS AN INITIALIZATION SEED, not a continuation of the source
trajectory: the `# coupling` line of the output says `mode=init t_phys=0
certified=F` whatever the source said (the source line is kept as
`# mapped-from-coupling:` provenance), because a nonconservative
interpolation is neither the certified state nor the physical clock of the
run that wrote it (finding B6 of the review of 2026-09-12). The `# grid`
line is taken from the target file, because a state describes the grid it is
written on and `load_IC` refuses a file whose grid field is not the run's
(the line is unchanged when the target is on the same grid). The other
header lines (schema, columns, restart_schema, reservoir, provenance) are
copied: they describe the composition and the origin, which the mapping
preserves unless `--reservoir` is given to change the composition on
purpose, and one `# mapped:` line records the mapping.

WHAT IS CHECKED BEFORE ANYTHING IS WRITTEN (finding B5): both source files
and the target carry finite, positive, strictly increasing radii; the two
source files have the same number of rows and the same radii to 1e-12
relative; the `# columns` lines carry the species the hydrogen-nucleus
count needs; the hydrogen-nucleus density is positive in every source row;
and every PHYSICAL target center (rows 3 to n-2) lies inside the source's
support, ghost rows included (the ghost rows are part of the state: the
lower ones carry the inflow reservoir composition of the molecular
handoff, so they are interpolated like every other row, never rebuilt from
the physical cells). A violation is refused with exit status 2 and a
message; nothing is written. Only a target GHOST row may lie outside the
source's range: it is extrapolated linearly in ln r from the two nearest
source rows, and the mapping line says so. That is the only extrapolation
the tool makes unless `--extrapolate-beyond` is given.

`--extrapolate-beyond <r/R_p>` EXTENDS the state above `r` instead of
interpolating it there, so that a solution on a 30 R_p domain can seed a run
whose domain reaches 45 or 60 R_p. Above `r` the state is not the source's;
it is the continuation of the source state at `r`:

    T(r') = T(r)                                   constant temperature
    rho(r') = rho(r) exp[-(phi(r') - phi(r))/(p/rho)(r)]   isothermal
                                                   hydrostatic stratification
    p(r') = p(r) x the same factor                 so p/rho, and with it T,
                                                   is the anchor's
    v(r') = v(r) [rho(r)/rho(r')] [r/r']^2         constant mass flux
    n_s(r')/n_H(r') = n_s(r)/n_H(r)                frozen composition

with phi = -GM/r' the spherical planetary potential and GM taken from the
source state itself, as the gravitational acceleration its own steady
momentum balance carries at the anchor,
`g = -(1/rho) dp/dr - v dv/dr` from one-sided differences of the source's
last interior rows, scaled outward as 1/r'^2. The density rule is the one
`states/Apply_BC.f90 free_outflow_ghost` writes into the outer ghosts, so
the seed and the boundary condition state the same stratification; the
velocity rule is the steady continuity equation, which the density rule does
not imply, so the two together are a seed and not a solution of the steady
equations. The volumetric heating and cooling columns are carried out at
constant rate per unit mass (they are diagnostics: `load_IC` reads r, v, p
and T from this file and the densities from the species file). `r` must lie
inside the source's rows, and rows at or below it are interpolated and
checked exactly as without the option. The `# mapped:` line records the
anchor, the derived GM, and how many rows were filled this way.

`--reservoir <El>/H <value>` MAKES THE STATE A SEED FOR ANOTHER
COMPOSITION. It may be given once per element -- `--reservoir He/H 2.13
--reservoir C/H 2.7780e-4` -- and every element it names is carried
independently. The columns of that element (He I, He II, He III and the
He 2^3S level for helium, the ionization stages El I, El II, El III for one
of the ten metals of `species_table`) are multiplied by
k = value / (the El/H of the source's `# reservoir` line) in every row,
ghosts included; the hydrogen columns and every element the option does not
name are left as they are, so each cell's n_s/n_H is carried exactly as the
mapping above carries it and the ionization split of the element is
untouched. The `# reservoir` line of both output files then states the new
ratios, and the `# mapped:` line records each rescaling, its factor and the
measured El/H of the base rows and of the column.

One factor over the whole column is what `load_IC.f90` does with the
element on a restart. For HELIUM with `He_diffusion` on, the He/H of a
solved state is not the reservoir value but a profile of it, pinned to the
reservoir at the base rows (the two inner ghosts and cell 1, where the
diffusion operator holds a Dirichlet He/H) and falling outward. Multiplying
every row by the one factor sets those base rows to the new reservoir
exactly and leaves He/H(r)/He/H(base) -- the shape the diffusion and the
advection balance produced -- unchanged. For a METAL the loader does the
same thing itself, by the same rule and for the same reason: the whole
column is multiplied by (handoff El/H)/(El/H at cell 1), one factor for
every ionization stage (`load_IC.f90`, the "elements the handoff states"
loop). That renormalization is reached only by a file with no metadata
block, because the restart contract compares the `# reservoir` line first
and refuses a state whose El/H is not the run's; carrying the state here is
what lets a run at another reservoir load it at all.

THE PRESSURE IS KEPT, AND THE DENSITY AND THE TEMPERATURE FOLLOW IT. The
loader reads r, v, p and T from the hydrodynamic file and rebuilds the mass
density from the species file, and the run takes the temperature from the
pressure and the particle count, T = p/((n_tot + n_e) k_B)
(`EXHALE_main.f90`, the ordering comment of the RK3 step). Adding nuclei at
fixed pressure therefore lowers the temperature by the factor by which
(n_tot + n_e) rises, and both columns are written that way here so that the
file is one state and not a state whose density and temperature belong to
the composition it no longer carries:

    rho  += sum over the rescaled columns of (k - 1) x (the mass one
            particle of that column carries), the weights `calc_rho` uses:
            3.9715259 per free helium nucleus, 4.9715259 per HeH+, and
            amu_over_m_H x A_El per metal nucleus, A_El the standard atomic
            weight in u (`species_table` melem_A_u, melem_A)
    T    = T_src x n_src/(n_src + dn),  n_src = p/(k_B T_src) read off the
           source state, dn the particles and electrons the added nuclei
           bring (He I 1, He II 2, He III 3, HeH+ 2; a metal stage El^(k+)
           1 + k)

n_src is taken from the source's own p and T, so every species the file
carries is inside it, and k_B is the value the state's `# constants` line
states. The metal mass enters rho under the `eos_metals` policy, which is
on by default and is what every state this code writes today carries.

The rescaling of an element is refused when the source states no reservoir
(an archived state written before the metadata block: `load_IC` rescales
such a file itself, by its own rule), when the source's `# reservoir` line
does not name that element, when the species file carries no column of it,
when the value is not positive and finite, or when the base rows do not come
out at the new El/H. The last is the test that catches a species carrying
nuclei of two elements at once: HeH+ holds one helium and one hydrogen
nucleus, so multiplying it changes the hydrogen count as well and no single
factor sets both, which is the same reason `load_IC` refuses a molecular
restart at a composition other than its own; the oxygen-chemistry carriers
OH, H2O and CO do the same for O and for C.

HOW THE VALUES ARE MAPPED. The hydrodynamic columns are interpolated in
ln r: density, pressure and the heating/cooling rates in log10 (values at
or below 1e-290, the writer's zero placeholder, are kept as that
placeholder), the velocity and the temperature linearly. The species are
NOT interpolated as densities: each species is carried as its ratio to the
hydrogen-nucleus density of the cell, n_s / n_H, interpolated linearly in
ln r, and n_H at the new center is the mapped mass density divided by the
mapped mass per hydrogen nucleus (rho / n_H of the source, interpolated). A
ratio that is the same in every source cell (He/H, the element abundances)
is then the same in every mapped cell to round-off, which is what load_IC
checks (He/H to 1e-6), and the mass closure of the state holds to the
interpolation error of a nearly constant quantity. The interpolation is
exact where the two grids coincide and its error scales with the grid
difference, so a state mapped across a 4e-4 relative shift in the centers
(the 2026-09-10 grid change against a 2026-09-05 state) moves by less than
the state's own truncation error. It is not a conservative remap: column-
integrated quantities are not preserved to round-off.
"""
import os, sys
import numpy as np

TINY = 1e-290
NGHOST = 2
# The density column is the mass density in units of the hydrogen atom mass,
# rho/m_H (`# columns rho[mH/cm3]`), with m_H the value the code carries,
# parameters.f90:787.  Needed only by --extrapolate-beyond, which puts the
# pressure and the velocity columns into one momentum balance.
M_H = 1.67353284e-24
# Helium-4 atom mass (g), parameters.f90:795; the ratio is the weight
# `calc_rho` gives He I, He II, He III and the He 2^3S level, and HeH+ takes
# that weight plus one hydrogen atom (utilities.f90 calc_rho).
M_HE = 6.6464790722e-24
M_HE_OVER_M_H = M_HE/M_H
H_NUC = {'HI': 1, 'HII': 1, 'H2': 2, 'H2p': 2, 'H3p': 3, 'HeHp': 1, 'OH': 1, 'H2O': 2}
# The helium nuclei a species carries, and its weight in hydrogen masses.
# He 2^3S is a level of He I and is already inside the He I column, so it is
# rescaled but contributes no nucleus and no mass of its own.
HE_NUC = {'HeI': 1, 'HeII': 1, 'HeIII': 1, 'HeITR': 0, 'HeHp': 1}
HE_MASS = {'HeI': M_HE_OVER_M_H, 'HeII': M_HE_OVER_M_H, 'HeIII': M_HE_OVER_M_H,
           'HeITR': 0.0, 'HeHp': M_HE_OVER_M_H + 1.0}
# What one particle of a helium species contributes to n_tot + n_e: itself
# plus the electrons it has released (He II 1, He III 2, HeH+ 1). That sum
# is the one the temperature is taken against, T = p/((n_tot + n_e) k_B).
HE_PARTICLES = {'HeI': 1.0, 'HeII': 2.0, 'HeIII': 3.0, 'HeITR': 0.0,
                'HeHp': 2.0}
# The ten metal elements of species_table.f90 and their standard atomic
# weights in the unified atomic mass unit u (CIAAW conventional values,
# species_table melem_A_u).  The code's mass unit is the hydrogen ATOM, so a
# weight enters a mass sum as amu_over_m_H x A_El (species_table melem_A);
# the two differ by 0.78 per cent.
AMU = 1.66053906660e-24          # unified atomic mass unit (g), parameters.f90:798
AMU_OVER_M_H = AMU/M_H
MELEM_A_U = {'C': 12.011, 'O': 15.999, 'N': 14.007, 'Mg': 24.305,
             'Si': 28.085, 'Ca': 40.078, 'Na': 22.990, 'K': 39.098,
             'S': 32.06, 'Fe': 55.845}
# The ionization stage a column name states, as the number of electrons the
# species has released: El I is neutral, El II singly ionized, and so on.
ROMAN_STAGE = {'I': 0, 'II': 1, 'III': 2, 'IV': 3, 'V': 4}
# The metal nuclei the oxygen-chemistry carriers hold (species_table bsp rows
# 11-13).  They are NOT rescaled with the ionization stages -- one factor
# cannot set two elements at once -- but they are counted in El/H, so a state
# that carries them is refused by the base-row test rather than mis-scaled.
CARRIER_NUC = {'OH': {'O': 1}, 'H2O': {'O': 1}, 'CO': {'C': 1, 'O': 1}}
# Boltzmann constant, CODATA 2018, the value parameters.f90 carries; the
# state's own '# constants' line is used instead where it states one.
KB = 1.380649e-16
REQUIRED_SPECIES = ('HI', 'HII')


def refuse(msg):
    sys.stderr.write('map_state_to_grid: REFUSED: ' + msg + '\n')
    sys.exit(2)


def read(path):
    header, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    if not rows:
        refuse(f'{path}: no data rows')
    a = np.array(rows)
    return header, a


def columns_of(header, path):
    cols = [h.split()[2:] for h in header if h.startswith('# columns')]
    if not cols:
        refuse(f'{path}: no "# columns" header line')
    return cols[0]


def check_radii(r, what):
    if r.size < 2*NGHOST + 2:
        refuse(f'{what}: only {r.size} rows; two ghost rows at each end and at least two physical cells are needed')
    if not np.all(np.isfinite(r)) or not np.all(r > 0):
        refuse(f'{what}: the radii must be finite and positive')
    if not np.all(np.diff(r) > 0):
        refuse(f'{what}: the radii must increase strictly')


def interp_log(x_new, x_old, y):
    y = np.asarray(y, dtype=float)
    small = y <= TINY
    if small.all():
        return np.full(x_new.size, y[0])
    ly = np.log10(np.where(small, TINY, y))
    out = 10.0**np.interp(x_new, x_old, ly)
    if small.any():
        out[out <= 10.0*TINY] = y[small][0]
    return out


def map_column(lr_new, lr_old_phys, y_phys, log=False):
    """Interpolate y (given on every source row) onto every target center;
    targets below or above the source range (ghost rows only, checked by the
    caller) are extrapolated linearly (in the transformed variable) from the
    two nearest source rows."""
    if log:
        small = np.asarray(y_phys) <= TINY
        if small.all():
            return np.full(lr_new.size, y_phys[0])
        yt = np.log10(np.where(small, TINY, y_phys))
    else:
        yt = np.asarray(y_phys, dtype=float)
    out = np.interp(lr_new, lr_old_phys, yt)
    lo = lr_new < lr_old_phys[0]
    hi = lr_new > lr_old_phys[-1]
    if lo.any():
        s = (yt[1] - yt[0])/(lr_old_phys[1] - lr_old_phys[0])
        out[lo] = yt[0] + s*(lr_new[lo] - lr_old_phys[0])
    if hi.any():
        s = (yt[-1] - yt[-2])/(lr_old_phys[-1] - lr_old_phys[-2])
        out[hi] = yt[-1] + s*(lr_new[hi] - lr_old_phys[-1])
    if log:
        out = 10.0**out
        if small.any():
            out[out <= 10.0*TINY] = np.asarray(y_phys)[small][0]
    return out


def gravity_at(r, rho_cgs, p, v, i):
    """The gravitational acceleration the source state carries at its row i,
    from its own steady momentum balance g = -(1/rho) dp/dr - v dv/dr with
    one-sided second-order differences taken INWARD, so that no ghost row --
    which is the source run's outer boundary condition and not its state --
    enters.  A steady solution satisfies this identically; a source that does
    not is refused below by the positivity test.

    r is the file's own radius column, in R_p, so the value returned is
    R_p times the acceleration in cgs; GM built from it as g r^2 divides by
    p/rho in cgs and gives the dimensionless exponent of the stratification,
    which is the only thing it is used for."""
    h1 = r[i] - r[i-1]
    h2 = r[i] - r[i-2]
    def back(y):
        # second-order backward difference on a non-uniform mesh
        return (y[i]*(h1 + h2)/(h1*h2) - y[i-1]*h2/(h1*(h2 - h1))
                + y[i-2]*h1/(h2*(h2 - h1)))
    return -back(p)/rho_cgs[i] - v[i]*back(v)


def hydrostatic_continuation(r_anchor, r, GM, cT2):
    """The isothermal hydrostatic density factor rho(r)/rho(r_anchor) of the
    spherical planetary potential phi = -GM/r, the same expression
    states/Apply_BC.f90 free_outflow_ghost writes into the outer ghosts. The
    radii and GM are in the units gravity_at returns, the isothermal sound
    speed cT2 = p/rho in cgs."""
    return np.exp(-(-GM/r + GM/r_anchor)/cT2)


def reservoir_map(header):
    """The elemental ratios the source state was solved at, from its
    '# reservoir' line ('# reservoir He/H <v> [<El>/H <v> ...]'), as
    {'He/H': v, 'C/H': v, ...}; None when the file carries no such line."""
    for h in header:
        if h.split()[:2] == ['#', 'reservoir']:
            toks = h.split()[2:]
            out = {}
            for i in range(0, len(toks) - 1, 2):
                try:
                    out[toks[i]] = float(toks[i+1])
                except ValueError:
                    pass
            return out
    return None


def reservoir_key_element(key):
    """The element symbol of a '<El>/H' reservoir key, or None when the key
    is not one of the elements this tool can carry."""
    if not key.endswith('/H'):
        return None
    sym = key[:-2]
    return sym if (sym == 'He' or sym in MELEM_A_U) else None


def reservoir_value(val, key):
    """The ratio of a '--reservoir <El>/H <value>' pair."""
    try:
        v = float(val)
    except ValueError:
        refuse(f'--reservoir {key}: "{val}" is not a nucleus ratio')
    if not (v > 0.0) or not np.isfinite(v):
        refuse(f'--reservoir {key} {val}: the ratio must be finite and positive')
    return v


def reservoir_line_at(header, new):
    """The source's '# reservoir' line with the ratios of `new` replaced and
    every other element left as it stands, written with the 17 significant
    digits that round trip a double."""
    for h in header:
        if h.split()[:2] == ['#', 'reservoir']:
            toks = h.split()[2:]
            for i in range(0, len(toks) - 1, 2):
                if toks[i] in new:
                    toks[i+1] = f'{new[toks[i]]:.16E}'
            return '# reservoir ' + ' '.join(toks)
    return ''


def boltzmann_of(header):
    """k_B as the state's '# constants' line states it, or the CODATA value
    when the file carries no such line."""
    for h in header:
        toks = h.split()
        if toks[:2] == ['#', 'constants'] and 'kB[erg/K]' in toks:
            try:
                return float(toks[toks.index('kB[erg/K]') + 1])
            except (ValueError, IndexError):
                return KB
    return KB


def metal_stage_of(name):
    """(element symbol, electrons released) of a metal ionization-stage
    column name such as 'CII', or None when the name is not one. The symbols
    are tried longest first so that 'CaII' is calcium and not carbon."""
    for sym in sorted(MELEM_A_U, key=len, reverse=True):
        if name.startswith(sym) and name[len(sym):] in ROMAN_STAGE:
            return sym, ROMAN_STAGE[name[len(sym):]]
    return None


def element_nuclei_in(name, el):
    """How many nuclei of element `el` one particle of species `name` holds.
    He 2^3S is a level of He I and is not counted again; HeH+ carries one
    helium and one hydrogen nucleus, and the oxygen-chemistry carriers hold
    the O and C nuclei of CARRIER_NUC, exactly as load_IC.f90 counts them."""
    if el == 'He':
        return HE_NUC.get(name, 0)
    stage = metal_stage_of(name)
    if stage is not None and stage[0] == el:
        return 1
    return CARRIER_NUC.get(name, {}).get(el, 0)


def element_stage_columns(el, scols):
    """The indices of the ionization-stage columns of element `el`: the ones
    the rescaling multiplies. A species holding nuclei of two elements at
    once (HeH+, OH, H2O, CO) is not among them for the metals -- no single
    factor sets two element counts -- and for helium the whole HE_NUC set is
    taken, HeH+ included, so that the base-row test refuses it."""
    if el == 'He':
        return [j for j, name in enumerate(scols) if name in HE_NUC]
    return [j for j, name in enumerate(scols)
            if (metal_stage_of(name) or (None,))[0] == el]


def element_mass_weight(name, el):
    """The mass one particle of the column carries, in hydrogen atoms, at the
    weights calc_rho uses."""
    if el == 'He':
        return HE_MASS[name]
    return AMU_OVER_M_H*MELEM_A_U[el]


def element_particle_weight(name, el):
    """What one particle of the column contributes to n_tot + n_e: itself
    plus the electrons it has released."""
    if el == 'He':
        return HE_PARTICLES[name]
    return 1.0 + metal_stage_of(name)[1]


def element_to_h(d, scols, el):
    """The element-to-hydrogen NUCLEUS ratio of every row of a species
    array, counted over every column that holds a nucleus of either."""
    nEl = np.zeros(d.shape[0])
    nH = np.zeros(d.shape[0])
    for j, name in enumerate(scols):
        k = element_nuclei_in(name, el)
        if k:
            nEl += k*d[:, j]
        if name in H_NUC:
            nH += H_NUC[name]*d[:, j]
    return nEl, nH


def rewrite_coupling(header, tgt_path, src_path, shift, n_new, n_ghost_extrap,
                     ext_note='', grid_line='', reservoir_line=''):
    out = []
    for h in header:
        if reservoir_line and h.split()[:2] == ['#', 'reservoir']:
            # The reservoir the state is a state OF. --reservoir has
            # carried the element onto another one, so the line states that
            # one; load_IC compares this field against the run's input and
            # refuses a state whose reservoir is not the run's.
            out.append(reservoir_line)
        elif grid_line and h.startswith('# grid '):
            # The state describes the grid it is WRITTEN on, not the one it
            # came from; load_IC refuses a file whose 'grid' field is not the
            # run's. The line is taken from the target file verbatim, so a
            # target written on the same grid leaves the header unchanged.
            out.append(grid_line)
        elif h.startswith('# rows'):
            out.append(f'# rows {n_new}: 2 ghost cells at each end; physical cells are rows 3 to {n_new-2}')
        elif h.startswith('# coupling:'):
            toks = h.split()[2:]
            kept = [t for t in toks if not (t.startswith('mode=') or t.startswith('t_phys=') or t.startswith('certified='))]
            out.append('# coupling: mode=init t_phys=0 certified=F ' + ' '.join(kept))
            out.append('# mapped-from-coupling: ' + ' '.join(toks))
        else:
            out.append(h)
    out.append(f'# mapped: onto the cell centers of {os.path.abspath(tgt_path)} from {os.path.abspath(src_path)} '
               f'by map_state_to_grid.py; largest relative center shift of the physical cells {shift:.2e}; '
               f'{n_ghost_extrap} ghost row(s) extrapolated linearly in ln r; an initialization seed, not a continuation'
               + ext_note)
    return out


def main():
    argv = sys.argv[1:]
    args, as_ic, r_ext, res_new = [], False, None, {}
    k = 0
    while k < len(argv):
        a = argv[k]
        if a == '--ic':
            as_ic = True
        elif a == '--reservoir':
            if k + 2 >= len(argv):
                refuse('--reservoir needs an element ratio and a value, '
                       'e.g. --reservoir He/H 2.13')
            key, val = argv[k+1], argv[k+2]
            k += 2
            if reservoir_key_element(key) is None:
                refuse(f'--reservoir {key}: the ratio must be written <El>/H with '
                       'El helium or one of the ten metals of species_table '
                       '(' + ' '.join(MELEM_A_U) + ')')
            if key in res_new:
                refuse(f'--reservoir {key}: given twice')
            res_new[key] = reservoir_value(val, key)
        elif a.startswith('--extrapolate-beyond'):
            val = a.split('=', 1)[1] if '=' in a else (argv[k+1] if k+1 < len(argv) else '')
            if '=' not in a:
                k += 1
            try:
                r_ext = float(val)
            except ValueError:
                refuse('--extrapolate-beyond needs a radius in R_p')
        elif a.startswith('--'):
            sys.exit(__doc__)
        else:
            args.append(a)
        k += 1
    if len(args) != 3:
        sys.exit(__doc__)
    src, tgt, out = args

    def src_file(n):
        p = os.path.join(src, n + '.txt')
        return p if os.path.exists(p) else os.path.join(src, n + '_IC.txt')

    # ---- read and validate ----
    hdr_t, t = read(tgt)
    r_new = t[:, 0]
    check_radii(r_new, f'target {tgt}')
    hdr_h, a = read(src_file('Hydro_ioniz'))
    hdr_s, c = read(src_file('Ion_species'))
    r_h, r_s = a[:, 0], c[:, 0]
    check_radii(r_h, f'source {src_file("Hydro_ioniz")}')
    check_radii(r_s, f'source {src_file("Ion_species")}')
    if r_h.size != r_s.size:
        refuse(f'the two source files have {r_h.size} and {r_s.size} rows')
    if np.max(np.abs(r_h - r_s)/r_h) > 1e-12:
        refuse('the two source files do not carry the same radii (largest relative difference '
               f'{np.max(np.abs(r_h - r_s)/r_h):.2e}); they are not one state')
    cols = columns_of(hdr_h, src_file('Hydro_ioniz'))
    scols = columns_of(hdr_s, src_file('Ion_species'))
    if cols[0] != 'r[Rp]' or scols[0] != 'r[Rp]':
        refuse('the first column of both source files must be r[Rp]')
    for name in REQUIRED_SPECIES:
        if name not in scols:
            refuse(f'the species file carries no {name} column')
    if len(cols) != a.shape[1] or len(scols) != c.shape[1]:
        refuse('a source file has a different number of columns from its "# columns" line')
    nH_old = np.zeros(c.shape[0])
    for j, name in enumerate(scols):
        if name in H_NUC:
            nH_old += H_NUC[name]*c[:, j]
    if not np.all(nH_old > 0):
        refuse('the hydrogen-nucleus density is not positive in every source row')
    # support: every source row, ghosts included
    ph_old = slice(0, r_h.size)
    ph_new = slice(NGHOST, r_new.size - NGHOST)
    lr_old = np.log(r_h)
    lr_new = np.log(r_new)
    if r_ext is None:
        if r_new[ph_new][0] < r_h[0]*(1 - 1e-12) or r_new[ph_new][-1] > r_h[-1]*(1 + 1e-12):
            refuse(f'the target physical cells span {r_new[ph_new][0]:.6f} to {r_new[ph_new][-1]:.6f} but the source '
                   f'rows only {r_h[0]:.6f} to {r_h[-1]:.6f}: no physical extension is implemented')
    else:
        if not (r_h[0]*(1 + 1e-12) < r_ext <= r_h[-1]*(1 + 1e-12)):
            refuse(f'--extrapolate-beyond {r_ext:.6f} is not inside the source rows '
                   f'{r_h[0]:.6f} to {r_h[-1]:.6f}: the anchor of the continuation must be a state the source carries')
        inside = r_new[ph_new][r_new[ph_new] <= r_ext*(1 + 1e-12)]
        if inside.size and inside[0] < r_h[0]*(1 - 1e-12):
            refuse(f'the target physical cells begin at {inside[0]:.6f} but the source rows only at {r_h[0]:.6f}: '
                   f'--extrapolate-beyond extends the state outward, never inward')
    n_ghost_extrap = int(np.sum(lr_new < lr_old[0]) + np.sum(lr_new > lr_old[-1]))
    if r_new.size == r_h.size:
        shift = float(np.max(np.abs(r_new[ph_new] - r_h[ph_new])/r_new[ph_new]))
    else:
        shift = float('nan')

    # ---- hydro ----
    b = np.empty((r_new.size, a.shape[1])); b[:, 0] = r_new
    for j in range(1, a.shape[1]):
        lin = cols[j].startswith('v[') or cols[j].startswith('T[')
        b[:, j] = map_column(lr_new, lr_old, a[ph_old, j], log=not lin)
    rho_new = b[:, 1]
    # ---- species as ratios to the hydrogen nuclei ----
    mu_old = a[ph_old, 1]/nH_old[ph_old]
    nH_new = rho_new/map_column(lr_new, lr_old, mu_old)
    d = np.empty((r_new.size, c.shape[1])); d[:, 0] = r_new
    for j in range(1, c.shape[1]):
        col = c[ph_old, j]
        if (col <= TINY).all():
            d[:, j] = col[0]
        else:
            d[:, j] = map_column(lr_new, lr_old, col/nH_old[ph_old])*nH_new

    # ---- the state above the anchor, when the target reaches past the source --
    ext_note = ''
    if r_ext is not None:
        ext = r_new > r_ext*(1 + 1e-12)
        n_ghost_extrap = int(np.sum(lr_new < lr_old[0]) + np.sum((lr_new > lr_old[-1]) & ~ext))
        def col_index(pred, what):
            hit = [j for j in range(len(cols)) if pred(cols[j])]
            if len(hit) != 1:
                refuse(f'the source names {len(hit)} columns that are the {what}: '
                       f'--extrapolate-beyond cannot be applied to it')
            return hit[0]
        j_rho = col_index(lambda c: c.startswith('rho'), 'mass density')
        j_v = col_index(lambda c: c.startswith('v['), 'velocity')
        j_p = col_index(lambda c: c.startswith('p['), 'pressure')
        j_T = col_index(lambda c: c.startswith('T['), 'temperature')
        lr_a = np.array([np.log(r_ext)])

        def at_anchor(y, log=True):
            return float(map_column(lr_a, lr_old, y, log=log)[0])
        rho_a = at_anchor(a[ph_old, j_rho])
        p_a = at_anchor(a[ph_old, j_p])
        v_a = at_anchor(a[ph_old, j_v], log=False)
        T_a = at_anchor(a[ph_old, j_T], log=False)
        # The gravity of the source's own steady momentum balance, at its last
        # row at or below the anchor; the two rows below it must be interior.
        ia = int(np.searchsorted(r_h, r_ext*(1 + 1e-12)) - 1)
        if ia < 2:
            refuse('--extrapolate-beyond: the anchor is within two rows of the '
                   'bottom of the source, so its gravity cannot be measured')
        gR = gravity_at(r_h, a[:, j_rho]*M_H, a[:, j_p], a[:, j_v], ia)
        if not (gR > 0.0):
            refuse(f'--extrapolate-beyond: the source carries an outward net force at r = {r_h[ia]:.4f} '
                   f'(g = {gR:.3e} in R_p units), so it states no stratification to continue')
        GMr = gR*r_h[ia]**2
        cT2_a = p_a/(rho_a*M_H)
        fac = hydrostatic_continuation(r_ext, r_new[ext], GMr, cT2_a)
        b[ext, j_rho] = rho_a*fac
        b[ext, j_p] = p_a*fac
        b[ext, j_T] = T_a
        b[ext, j_v] = v_a/fac*(r_ext/r_new[ext])**2
        for j in range(1, a.shape[1]):
            if j not in (j_rho, j_v, j_p, j_T):
                b[ext, j] = at_anchor(a[ph_old, j])*fac
        mu_a = at_anchor(a[ph_old, j_rho]/nH_old[ph_old], log=False)
        nH_ext = b[ext, j_rho]/mu_a
        for j in range(1, c.shape[1]):
            col = c[ph_old, j]
            if not (col <= TINY).all():
                d[ext, j] = at_anchor(col/nH_old[ph_old], log=False)*nH_ext
        ext_note = (f'; {int(np.sum(ext))} row(s) above r = {r_ext:.6f} R_p are NOT the source state but its '
                    f'continuation: constant T, isothermal hydrostatic rho and p under GM/R_p = {GMr:.6e} '
                    f'(cm/s)^2 R_p measured from the source momentum balance at r = {r_h[ia]:.6f}, '
                    f'constant mass flux v, frozen composition')

    # ---- the elements carried onto another reservoir -----------------------
    res_note, res_h, res_s = '', '', ''
    if res_new:
        res_file_h = reservoir_map(hdr_h)
        res_file_s = reservoir_map(hdr_s)
        if res_file_h is None or res_file_s is None:
            refuse('--reservoir: the source states no "# reservoir" line, so the '
                   'composition it was solved at is not recorded and no factor can be formed '
                   '(an archived state written before the metadata block: load_IC rescales '
                   'such a file itself)')
        # The mass density the state would be written with: the nuclei the
        # rescaling adds, at the weights calc_rho gives them.
        rho_cols = [j for j in range(len(cols)) if cols[j].startswith('rho')]
        if len(rho_cols) != 1:
            refuse(f'--reservoir: the source names {len(rho_cols)} columns that are the '
                   'mass density, so the added mass cannot be put back into it')
        j_T_r = [j for j in range(len(cols)) if cols[j].startswith('T[')]
        j_p_r = [j for j in range(len(cols)) if cols[j].startswith('p[')]
        if len(j_T_r) != 1 or len(j_p_r) != 1:
            refuse('--reservoir: the source does not name exactly one pressure and one '
                   'temperature column, so the state it becomes cannot be written')
        kb = boltzmann_of(hdr_h)
        # The particle count the pressure and the temperature of the source
        # stand for, read off the state itself so that every species the file
        # carries, metals included, is in it.
        T_src = b[:, j_T_r[0]].copy()
        if not np.all(T_src > 0.0):
            refuse('--reservoir: the source temperature is not positive in every row')
        n_src = b[:, j_p_r[0]]/(kb*T_src)
        d_rho = np.zeros(r_new.size)
        d_part = np.zeros(r_new.size)
        carried = {}
        for key, val in res_new.items():
            el = reservoir_key_element(key)
            if key not in res_file_h or key not in res_file_s:
                refuse(f'--reservoir {key}: the source "# reservoir" line does not name '
                       f'{key}, so the ratio the state was solved at is not recorded and '
                       'no factor can be formed')
            old = res_file_h[key]
            if not (old > 0.0):
                refuse(f'--reservoir {key}: the source reservoir is {old:.6E}, not positive')
            if abs(res_file_s[key] - old) > 1e-12*old:
                refuse(f'--reservoir {key}: the two source files state different reservoirs '
                       f'({old:.16E} and {res_file_s[key]:.16E}); they are not one state')
            el_cols = element_stage_columns(el, scols)
            if not el_cols:
                refuse(f'--reservoir {key}: the species file carries no {el} column')
            fac = val/old
            for j in el_cols:
                d_rho += element_mass_weight(scols[j], el)*d[:, j]*(fac - 1.0)
                d_part += element_particle_weight(scols[j], el)*d[:, j]*(fac - 1.0)
                d[:, j] *= fac
            carried[key] = (el, old, val, fac)
        b[:, rho_cols[0]] = b[:, rho_cols[0]] + d_rho
        # The pressure is held and the temperature follows the particle count,
        # which is the relation the run itself carries; writing it here leaves
        # the file one state rather than a state whose T column belongs to the
        # composition it no longer has.
        if not np.all(n_src + d_part > 0.0):
            refuse('--reservoir: the rescaled particle count is not positive in every row')
        b[:, j_T_r[0]] = T_src*n_src/(n_src + d_part)
        # The base rows come out at the new ratio, or the state carries a
        # species holding nuclei of two elements at once and one factor cannot
        # set both counts. For helium the base rows are the two inner ghosts
        # and cell 1, where the element-diffusion operator holds its Dirichlet
        # He/H; for a metal it is cell 1, the row load_IC.f90 forms its own
        # renormalization factor at.
        notes = []
        for key, (el, old, val, fac) in carried.items():
            nEl, nH = element_to_h(d, scols, el)
            q = nEl/np.where(nH > 0.0, nH, 1.0)
            if el == 'He':
                base, where = q[:NGHOST + 1], 'base rows'
            else:
                base, where = q[NGHOST:NGHOST + 1], 'first physical cell'
            dev = float(np.max(np.abs(base - val)))/val
            if dev > 1e-12:
                refuse(f'--reservoir {key}: the {where} come out at {key} {base.min():.6E} '
                       f'to {base.max():.6E} instead of {val:.6E} (relative departure '
                       f'{dev:.2E}). A state carrying a species with nuclei of two elements '
                       'at once is the case in point: HeH+ holds one helium and one hydrogen '
                       'nucleus, and OH, H2O and CO hold oxygen and carbon, so multiplying '
                       'the ionization stages alone does not move the element count by the '
                       'factor asked for')
            phys = q[NGHOST:r_new.size - NGHOST]
            notes.append(f'the {el} was carried from the reservoir {key} {old:.6E} to '
                         f'{val:.6E} by one factor {fac:.16E} in every row (initialization '
                         f'choice: the column keeps its shape and its base rows take the new '
                         f'reservoir), {key} of the physical cells {phys.min():.6E} to '
                         f'{phys.max():.6E}, {where} {base.max():.6E}')
        res_note = ('; ' + '; '.join(notes)
                    + f'; the pressure is unchanged and the temperature follows the particle '
                      f'count, T {T_src[NGHOST]:.6E} -> {b[NGHOST, j_T_r[0]]:.6E} K at the '
                      f'first physical cell')
        res_h = reservoir_line_at(hdr_h, res_new)
        res_s = reservoir_line_at(hdr_s, res_new)

    tgt_grid_line = ''
    for h in hdr_t:
        if h.startswith('# grid '):
            tgt_grid_line = h

    os.makedirs(out, exist_ok=True)
    for n, hdr, arr, p, res in (('Hydro_ioniz', hdr_h, b, src_file('Hydro_ioniz'), res_h),
                                ('Ion_species', hdr_s, d, src_file('Ion_species'), res_s)):
        hdr = rewrite_coupling(hdr, tgt, p, shift, r_new.size, n_ghost_extrap,
                               ext_note + res_note, tgt_grid_line, res)
        dst = os.path.join(out, n + ('_IC' if as_ic else '') + '.txt')
        with open(dst, 'w') as f:
            f.write('\n'.join(hdr) + '\n')
            for row in arr:
                f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
        print('wrote', dst, arr.shape)


if __name__ == '__main__':
    main()
