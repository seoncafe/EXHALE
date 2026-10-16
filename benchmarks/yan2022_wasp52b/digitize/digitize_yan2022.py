"""Digitize Yan et al. (2022) ApJ 936, 177 profile curves by color from a
300 dpi render of the publisher PDF (pdftoppm -r 300). For each pixel column
inside the frame, pixels within a color distance of the legend color are
grouped; one or two tracks are followed column by column (nearest group to
the previous value). Output: r/Rp and value, every column with a detection."""
import sys, json
import numpy as np
from PIL import Image

COL = {'0.25': (255, 0, 41), '0.5': (255, 186, 0), '1': (92, 255, 0),
       '2': (0, 255, 140), '4': (0, 143, 255), '6': (91, 0, 255),
       '8': (255, 0, 191), 'red': (255, 0, 0)}


def load(png):
    return np.array(Image.open(png).convert('RGB')).astype(float)


def groups(idx):
    g, cur = [], []
    for i in idx:
        if cur and i - cur[-1] > 2:
            g.append(cur); cur = []
        cur.append(i)
    if cur:
        g.append(cur)
    return [float(np.mean(c)) for c in g]


def extract(im, frame, color, ntrack=1, tol=70.0, start='right', exclude=()):
    x0, x1, y0, y1 = frame
    xs = np.arange(int(np.ceil(x0)) + 3, int(np.floor(x1)) - 2)
    if start == 'right':
        xs = xs[::-1]
    c = np.array(color, float)
    tracks = [[] for _ in range(ntrack)]
    last = [None] * ntrack
    for x in xs:
        colpix = im[int(y0) + 3:int(y1) - 2, x]
        d = np.sqrt(((colpix - c) ** 2).sum(axis=1))
        idx = np.where(d < tol)[0]
        for (ex0, ex1, ey0, ey1) in exclude:
            if ex0 <= x <= ex1:
                yy = idx + int(y0) + 3
                idx = idx[(yy < ey0) | (yy > ey1)]
        gs = [g + int(y0) + 3 for g in groups(list(idx))]
        if not gs:
            continue
        if ntrack == 1:
            if last[0] is None:
                y = gs[0] if len(gs) == 1 else None
            else:
                y = min(gs, key=lambda g: abs(g - last[0]))
                if abs(y - last[0]) > 40:
                    y = None
            if y is not None:
                tracks[0].append((x, y)); last[0] = y
        else:
            # two tracks: track 0 is the upper curve (smaller pixel y)
            if last[0] is None:
                if len(gs) >= 2:
                    last = [min(gs), max(gs)]
                    tracks[0].append((x, last[0])); tracks[1].append((x, last[1]))
                continue
            used = set()
            for k in range(2):
                cand = [g for g in gs if g not in used]
                if not cand:
                    continue
                pts = tracks[k][-6:]
                if len(pts) >= 3:
                    sl = (pts[-1][1] - pts[0][1]) / (pts[-1][0] - pts[0][0])
                    pred = pts[-1][1] + sl * (x - pts[-1][0])
                else:
                    pred = last[k]
                y = min(cand, key=lambda g: abs(g - pred))
                gap = abs(x - tracks[k][-1][0])
                if abs(y - pred) <= 4 + 0.5 * gap:
                    tracks[k].append((x, y)); last[k] = y; used.add(y)
    return tracks


def to_data(track, frame, ycal, ylog):
    x0, x1, y0, y1 = frame
    out = []
    for x, y in sorted(track):
        r = 10 ** ((x - x0) / (x1 - x0))
        ypix0, val0, ppu = ycal  # pixel of reference value, value, pixels per unit (or decade)
        v = (ypix0 - y) / ppu
        val = 10 ** (np.log10(val0) + v) if ylog else val0 + v
        out.append((r, val))
    return out
