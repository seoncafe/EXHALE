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

The hydrodynamic columns are interpolated in ln r: density, pressure and
the heating/cooling rates in log10 (values at or below 1e-290, the writer's
zero placeholder, are kept as that placeholder), the velocity and the
temperature linearly. The species are NOT interpolated as densities: each
species is carried as its ratio to the hydrogen-nucleus density of the cell,
n_s / n_H, interpolated linearly in ln r, and n_H at the new center is the
mapped mass density divided by the mapped mass per hydrogen nucleus
(rho / n_H of the source, interpolated). A ratio that is the same in every
source cell (He/H, the element abundances) is then the same in every mapped
cell to round-off, which is what load_IC checks (He/H to 1e-6), and the
mass closure of the state holds to the interpolation error of a nearly
constant quantity. The header lines of the source are copied and one
comment line records the mapping; the row count is the target's. The
interpolation is exact where the two grids coincide and its error scales
with the grid difference, so a state mapped across a 4e-4 relative shift in
the centers (the 2026-09-10 grid change against a 2026-09-05 state) moves by
less than the state's own truncation error. It is not a conservative remap:
column-integrated quantities are not preserved to round-off, and a run
started from a mapped state is a restart, not a continuation.
"""
import os, sys
import numpy as np

TINY = 1e-290

def read(path):
    header, rows = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    return header, np.array(rows)

def interp_log(x_new, x_old, y):
    y = np.asarray(y, dtype=float)
    small = y <= TINY
    if small.all():
        return np.full(x_new.size, y[0])
    ly = np.log10(np.where(small, TINY, y))
    out = 10.0**np.interp(x_new, x_old, ly)
    out[out <= 10.0*TINY] = y[small][0] if small.any() else out[out <= 10.0*TINY]
    return out

H_NUC = {'HI': 1, 'HII': 1, 'H2': 2, 'H2p': 2, 'H3p': 3, 'HeHp': 1, 'OH': 1, 'H2O': 2}

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    as_ic = '--ic' in sys.argv
    if len(args) != 3:
        sys.exit(__doc__)
    src, tgt, out = args
    def src_file(n):
        p = os.path.join(src, n + '.txt')
        return p if os.path.exists(p) else os.path.join(src, n + '_IC.txt')
    r_new = read(tgt)[1][:, 0]
    lr_new = np.log(r_new)
    os.makedirs(out, exist_ok=True)
    # --- hydro ---
    hdr_h, a = read(src_file('Hydro_ioniz'))
    lr_old = np.log(a[:, 0])
    cols = [h.split()[2:] for h in hdr_h if h.startswith('# columns')][0]
    b = np.empty((r_new.size, a.shape[1])); b[:, 0] = r_new
    for j in range(1, a.shape[1]):
        if cols[j].startswith('v[') or cols[j].startswith('T['):
            b[:, j] = np.interp(lr_new, lr_old, a[:, j])
        else:
            b[:, j] = interp_log(lr_new, lr_old, a[:, j])
    rho_new = b[:, 1]
    # --- species, as ratios to the hydrogen-nucleus density ---
    hdr_s, c = read(src_file('Ion_species'))
    scols = [h.split()[2:] for h in hdr_s if h.startswith('# columns')][0]
    nH_old = np.zeros(c.shape[0])
    for j, name in enumerate(scols):
        if name in H_NUC:
            nH_old += H_NUC[name]*c[:, j]
    mu_old = a[:, 1]/nH_old                       # mass per H nucleus, in the file's units
    nH_new = rho_new/np.interp(lr_new, lr_old, mu_old)
    d = np.empty((r_new.size, c.shape[1])); d[:, 0] = r_new
    for j in range(1, c.shape[1]):
        col = c[:, j]
        if (col <= TINY).all():
            d[:, j] = col[0]
        else:
            d[:, j] = np.interp(lr_new, lr_old, col/nH_old)*nH_new
    shift = np.max(np.abs(r_new - a[:, 0])/r_new) if r_new.size == a.shape[0] else float('nan')
    for n, hdr, arr, p in (('Hydro_ioniz', hdr_h, b, src_file('Hydro_ioniz')),
                           ('Ion_species', hdr_s, d, src_file('Ion_species'))):
        hdr = [h if not h.startswith('# rows') else
               f'# rows {r_new.size}: 2 ghost cells at each end; physical cells are rows 3 to {r_new.size-2}'
               for h in hdr]
        hdr.append(f'# mapped: onto the cell centers of {os.path.abspath(tgt)} from {os.path.abspath(p)} '
                   f'by map_state_to_grid.py; largest relative center shift {shift:.2e}')
        dst = os.path.join(out, n + ('_IC' if as_ic else '') + '.txt')
        with open(dst, 'w') as f:
            f.write('\n'.join(hdr) + '\n')
            for row in arr:
                f.write(' '.join(f'{x:24.16E}' for x in row) + '\n')
        print('wrote', dst, arr.shape)

if __name__ == '__main__':
    main()
