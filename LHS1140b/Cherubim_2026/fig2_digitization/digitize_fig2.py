#!/usr/bin/env python3
"""Digitize Fig. 2 of Cherubim et al. (2026, Science) -- the raw He 10830
transmission spectrum of LHS 1140b from the 2024 transit.

WHAT THIS CAN AND CANNOT RECOVER.  The figure is a raster image in the
published PDF, so there is no vector coordinate list to read, and the markers
are filled discs 24 px across drawn about 5 px apart: neighbouring points
overlap by four fifths of their own width.  The union of overlapping opaque
discs does not determine their individual centres, so the observed points
cannot be recovered one by one, and no attempt is made to.

What is recovered, at each column of the plot, is the interval of marker
centres consistent with the ink -- every position where a full 24 px disc
fits -- and the vertical extent of black through that column, which is that
interval widened by the error bars.  Published points lie inside the first
band; the second is the band plus their uncertainties.  Plot them as bands.

Axis calibration is checked against the figure's own annotation: the three
blue vertical lines mark the He I triplet at vacuum 10832.057, 10833.217 and
10833.306 A, and the check must pass to within 0.05 A or the script stops.

Sign: the paper plots absorption NEGATIVE.  Columns 2-3 keep the paper's own
axis; columns 6-7 flip it to absorption-positive, the convention used by
EXHALE_transit.py and by the p-winds oracle.
"""
import numpy as np
from PIL import Image
from scipy import ndimage as ndi

SRC, OUT = 'fig2/p11-000.png', 'lhs1140b_fig2_band.txt'
RAD, FILL, STEP = 12, 0.995, 2

im = np.asarray(Image.open(SRC).convert('RGB')).astype(int)
r, g, b = im[..., 0], im[..., 1], im[..., 2]
blk = (r < 110) & (g < 110) & (b < 110)
H, W = blk.shape
rows = np.where(blk.sum(axis=1) > 0.6*W)[0]
cols = np.where(blk.sum(axis=0) > 0.6*H)[0]
T, B, L, R = rows.min(), rows.max(), cols.min(), cols.max()
x2lam = lambda x: 10825.0 + (x - L)/(R - L)*15.0
y2val = lambda y: 1.0 + (y - T)/(B - T)*(-3.0)

blue = (b - r > 40) & (b - g > 25)
cand = np.where(blue[T+5:B-5, :].sum(axis=0) > 0.7*(B - T - 10))[0]
groups, cur = [], [cand[0]]
for c in cand[1:]:
    if c - cur[-1] <= 2: cur.append(c)
    else: groups.append(cur); cur = [c]
groups.append(cur)
seen = np.array([x2lam(np.mean(gp)) for gp in groups])
TRUE = np.array([10832.057, 10833.217, 10833.306])
err = np.abs(seen - TRUE).max()
print('triplet check: recovered', np.round(seen, 3), 'A   max |err| = %.3f A' % err)
assert err < 0.05, 'axis calibration failed'

m = np.zeros_like(blk)
m[T+4:B-3, L+4:R-3] = blk[T+4:B-3, L+4:R-3]
yy, xx = np.mgrid[-RAD:RAD+1, -RAD:RAD+1]
disc = ((xx**2 + yy**2) <= RAD**2).astype(float); disc /= disc.sum()
corr = ndi.correlate(m.astype(float), disc, mode='constant')

out = []
for x in range(L+RAD+2, R-RAD-1, STEP):
    fit = np.where(corr[:, x] >= FILL)[0]
    if fit.size == 0:
        continue
    ink = np.where(m[:, x])[0]
    out.append((x2lam(x), y2val(fit.min()), y2val(fit.max()),
                y2val(ink.min()), y2val(ink.max())))

out = np.array(out)
lam, chi, clo, ehi, elo = out.T          # 'hi' = less absorption in paper sign
hdr = (__doc__.strip() + '\n\n'
       'col1 lambda[A] vacuum, planetary rest frame\n'
       'col2 marker band upper edge [%]  (paper sign, absorption negative)\n'
       'col3 marker band lower edge [%]  (paper sign)\n'
       'col4 with error bars, upper [%]  (paper sign)\n'
       'col5 with error bars, lower [%]  (paper sign)\n'
       'col6 marker band, absorption POSITIVE, lower edge = -col2\n'
       'col7 marker band, absorption POSITIVE, upper edge = -col3')
np.savetxt(OUT, np.c_[lam, chi, clo, ehi, elo, -chi, -clo], header=hdr,
           fmt='%12.5f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f')
print('wrote %s: %d columns, %.2f-%.2f A' % (OUT, len(lam), lam.min(), lam.max()))
i = np.argmin(clo)
print('deepest published point: %.3f A, %.3f %% (paper sign)' % (lam[i], clo[i]))
c = (lam < 10831) | (lam > 10836)
print('continuum band away from the line: %.3f to %.3f %% (paper sign)'
      % (np.median(clo[c]), np.median(chi[c])))
