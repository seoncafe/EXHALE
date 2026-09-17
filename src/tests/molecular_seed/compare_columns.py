#!/usr/bin/env python3
"""Two state pairs, column by column, over the physical cells.

Prints two lines: the largest relative difference of every column the two
pairs share, and the largest absolute value of the four molecular columns of
the second pair.  Used by run.sh for the conversion-identity assertion, where
both numbers must be zero to round-off.

The ghost rows are excluded: the conversion re-applies the target run's own
boundary, so the ghost rows state the boundary and not the transfer.  How many
there are is read from the '# rows' line each file carries.
"""
import sys
import numpy as np


def load(path):
    labels = None
    ghosts = 2
    rows = []
    for line in open(path):
        if line.startswith('#'):
            tok = line[1:].split()
            if tok and tok[0] == 'columns':
                labels = tok[1:]
            elif tok and tok[0] == 'rows':
                # '# rows 504: 2 ghost cells at each end; ...'
                for i, t in enumerate(tok):
                    if t == 'ghost':
                        ghosts = int(tok[i - 1])
                        break
            continue
        rows.append([float(x) for x in line.split()])
    return labels, ghosts, np.array(rows)


def worst(a_dir, b_dir):
    wmax = 0.0
    molmax = 0.0
    for base in ('Hydro_ioniz', 'Ion_species'):
        la, ga, A = load(f'{a_dir}/{base}_IC.txt')
        lb, gb, B = load(f'{b_dir}/{base}_IC.txt')
        sl = slice(ga, A.shape[0] - ga)
        for i, name in enumerate(la):
            # The radiative rates of a seed are written as zero by design:
            # no sweep has run, so they are not a column of the state.
            if name.startswith('heat') or name.startswith('cool'):
                continue
            if name not in lb:
                continue
            a = A[sl, i]
            b = B[slice(gb, B.shape[0] - gb), lb.index(name)]
            den = np.maximum(np.abs(a), np.abs(b))
            d = np.where(den > 0.0, np.abs(a - b) / np.where(den > 0.0, den, 1.0), 0.0)
            wmax = max(wmax, float(d.max()))
        if base == 'Ion_species':
            for name in ('H2', 'H2p', 'H3p', 'HeHp'):
                if name in lb:
                    molmax = max(molmax, float(np.abs(B[:, lb.index(name)]).max()))
    return wmax, molmax


if __name__ == '__main__':
    w, m = worst(sys.argv[1], sys.argv[2])
    print(f'{w:.6e}')
    print(f'{m:.6e}')
