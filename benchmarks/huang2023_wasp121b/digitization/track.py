"""Trace one curve of a given color through a figure image by continuity.

Pixels of the color class are grouped per column into vertical clusters; from
a seed (column, row) the tracer walks left and right, taking in each column
the cluster nearest to the last accepted row within max_jump pixels. Columns
with no acceptable cluster are skipped (dashed or dotted lines). Returns
(column, row) pairs of cluster centres.
"""
import numpy as np
from PIL import Image

def load(path):
    return np.asarray(Image.open(path).convert('RGB')).astype(int)

def mask_color(im, kind):
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    if kind == 'black':
        return (R < 90) & (G < 90) & (B < 200)
    if kind == 'red':
        return (R > 170) & (G < 110) & (B < 110)
    if kind == 'blue':
        return (B > 170) & (R < 90) & (G < 90)
    if kind == 'magenta':
        return (R > 150) & (B > 150) & (G < 110)
    if kind == 'olive':
        return (R > 130) & (G > 130) & (B < 90)
    if kind == 'cyan':
        return (G > 150) & (B > 150) & (R < 100)
    raise ValueError(kind)

def clusters(col):
    idx = np.where(col)[0]
    out = []
    for i in idx:
        if out and i - out[-1][-1] <= 1:
            out[-1].append(i)
        else:
            out.append([i])
    return [((g[0] + g[-1]) / 2.0, len(g)) for g in out]

def trace(m, seed_col, seed_row, x0, x1, max_jump=12, max_len=40, pick=None):
    pts = {}
    for direction in (+1, -1):
        last = seed_row
        c = seed_col
        miss = 0
        while x0 <= c <= x1:
            cl = [cc for cc, n in clusters(m[:, c]) if n <= max_len]
            if cl:
                d = [abs(cc - last) for cc in cl]
                lim = max_jump * (1 + miss)
                k = int(np.argmin(d))
                ok = [i for i in range(len(cl)) if d[i] <= lim]
                if pick == 'lowest' and ok:      # largest row = lowest value
                    k = max(ok, key=lambda i: cl[i])
                elif pick == 'highest' and ok:
                    k = min(ok, key=lambda i: cl[i])
                if d[k] <= lim:
                    pts[c] = cl[k]
                    last = cl[k]
                    miss = 0
                else:
                    miss += 1
            else:
                miss += 1
            c += direction
    cols = sorted(pts)
    return np.array(cols, float), np.array([pts[c] for c in cols])
