#!/usr/bin/env python3
"""Map an EXHALE state (Hydro_ioniz.txt + Ion_species.txt) onto the cell
centers of another grid, so that a state written on an earlier grid can be
loaded (`Load IC? True`) by a run whose grid construction has since changed
(load_IC refuses a state whose centers differ from the run's beyond 1e-10).

Usage:
  map_state_to_grid.py <src_dir> <target_grid_file> <out_dir> [--ic]

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
run that wrote it (finding B6 of the review of 2026-09-12). The other
header lines (schema, columns, restart_schema, reservoir, provenance) are
copied: they describe the composition and the origin, which the mapping
preserves, and one `# mapped:` line records the mapping.

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
the tool makes.

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
H_NUC = {'HI': 1, 'HII': 1, 'H2': 2, 'H2p': 2, 'H3p': 3, 'HeHp': 1, 'OH': 1, 'H2O': 2}
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


def rewrite_coupling(header, tgt_path, src_path, shift, n_new, n_ghost_extrap):
    out = []
    for h in header:
        if h.startswith('# rows'):
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
               f'{n_ghost_extrap} ghost row(s) extrapolated linearly in ln r; an initialization seed, not a continuation')
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    as_ic = '--ic' in sys.argv
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
    if r_new[ph_new][0] < r_h[0]*(1 - 1e-12) or r_new[ph_new][-1] > r_h[-1]*(1 + 1e-12):
        refuse(f'the target physical cells span {r_new[ph_new][0]:.6f} to {r_new[ph_new][-1]:.6f} but the source '
               f'rows only {r_h[0]:.6f} to {r_h[-1]:.6f}: no physical extension is implemented')
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

    os.makedirs(out, exist_ok=True)
    for n, hdr, arr, p in (('Hydro_ioniz', hdr_h, b, src_file('Hydro_ioniz')),
                           ('Ion_species', hdr_s, d, src_file('Ion_species'))):
        hdr = rewrite_coupling(hdr, tgt, p, shift, r_new.size, n_ghost_extrap)
        dst = os.path.join(out, n + ('_IC' if as_ic else '') + '.txt')
        with open(dst, 'w') as f:
            f.write('\n'.join(hdr) + '\n')
            for row in arr:
                f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
        print('wrote', dst, arr.shape)


if __name__ == '__main__':
    main()
