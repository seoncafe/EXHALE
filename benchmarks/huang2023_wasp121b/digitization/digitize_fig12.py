"""Digitize Huang et al. (2023) Figure 12 (case A number densities), r >= 1.17 R_p.

Run from this directory: python3 digitize_fig12.py . ; writes fig12.npy."""
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from track import load, mask_color, trace, clusters
from overlay import overlay

xt = [251.0, 340.0, 415.0, 510.0, 589.5, 678.5, 787.5, 928.0, 1037.0]
xv = [1.0, 1.2, 1.4, 1.7, 2.0, 2.4, 3.0, 4.0, 5.0]
yt = [104.0, 192.0, 280.0, 368.0, 455.5, 543.5, 631.5]
yv = [12, 10, 8, 6, 4, 2, 0]
ax, bx = np.polyfit(xt, np.log10(xv), 1)
ay, by = np.polyfit(yt, yv, 1)
print('tick fit residuals: logR %.4f dex, logn %.3f dex' % (
    np.abs(np.log10(xv) - (ax * np.array(xt) + bx)).max(),
    np.abs(np.array(yv) - (ay * np.array(yt) + by)).max()))
c_lo = int(np.ceil((np.log10(1.17) - bx) / ax))
im = load('huang2023_fig12.png')
seeds = {  # name: (color, seed column, seed row, pick)
    'e_and_p': ('olive', 700, 269.0, None),
    'H2p': ('black', 700, 658.5, 'highest'),
    'MgII': ('red', 700, 514.0, 'highest'),
    'MgI': ('red', 700, 601.0, 'lowest'),
    'FeII': ('blue', 700, 526.0, 'highest'),
    'FeI': ('blue', 700, 665.0, 'lowest'),
    'CaII': ('magenta', 700, 583.0, 'highest'),
}
res, px = {}, []
for name, (kind, sc, sr, pick) in seeds.items():
    m = mask_color(im, kind)
    m[:106, :] = False; m[717:, :] = False
    m[:, :c_lo] = False; m[:, 1037:] = False
    m[315:452, 515:835] = False          # legend
    m[470:530, 215:300] = False          # "Lower atm." label
    c_x = int(round((np.log10(1.33) - bx) / ax))   # H(2p) and Na cross here
    if name == 'H2p':
        # H(2p) and Na (both black) cross near 1.33 R_p; the Na curve is not
        # separable from H(2p) by this tracer, so H(2p) is kept outside the
        # crossing only and Na is not digitized.
        m[:, :c_x] = False
    if name == 'FeI':
        m[:, int((np.log10(3.85) - bx) / ax):] = False
    cols, rows = trace(m, sc, sr, c_lo, 1036, max_jump=8, pick=pick)
    px.append((cols, rows))
    res[name] = (10 ** (ax * cols + bx), 10 ** (ay * rows + by))
    print(name, len(cols), 'R %.3f-%.3f' % (res[name][0].min(), res[name][0].max()))
overlay('huang2023_fig12.png', px, 'huang2023_fig12_trace_check.png')
np.save(sys.argv[1] + '/fig12.npy', {k: np.vstack(v) for k, v in res.items()}, allow_pickle=True)
