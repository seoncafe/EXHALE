import pickle, os, numpy as np
D = os.environ.get('YAN_RENDER_DIR', os.getcwd())
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
res = pickle.load(open(os.path.join(D, 'digitized.pkl'), 'rb'))
F = {'0.5': 0.5, '1': 1.0, '2': 2.0, 'red': 1.0}


def arr(key):
    d = np.array(res.get(key, []))
    return d[np.argsort(d[:, 0])] if len(d) else d


def crossing(a, b):
    # radius where two tracks meet (min |log ratio| on the common r range)
    rr = np.linspace(max(a[0, 0], b[0, 0]), min(a[-1, 0], b[-1, 0]), 400)
    da = np.interp(rr, a[:, 0], np.log10(a[:, 1])) - np.interp(rr, b[:, 0], np.log10(b[:, 1]))
    m = (rr < 3.0) & (rr > 1.2)
    return rr[m][np.argmin(np.abs(da[m]))]


rows9, rows23, rows12 = [], [], []
for c in ['0.5', '1', '2']:
    for q, panel in (('T_K', 'fig9a'), ('v_kms', 'fig9b')):
        for r, v in arr((panel, c, 0))[::3]:
            rows9.append((q, F[c], r, v))
    up, lo = arr(('fig9c', c, 0)), arr(('fig9c', c, 1))
    rc = crossing(up, lo)
    # inside the crossing H(1s) lies above H+, outside it below; each traced
    # point is assigned by its position against the other trace there, and
    # points where the two traces coincide (one curve traced twice) are dropped
    def split(a, b):
        ob = np.interp(a[:, 0], b[:, 0], np.log10(b[:, 1]))
        dl = np.log10(a[:, 1]) - ob
        return a[dl > 0.02], a[dl < -0.02]
    ua, ub = split(up, lo)
    la, lb = split(lo, up)
    hi_pts = np.vstack([ua, la]); lo_pts = np.vstack([ub, lb])
    hi_pts = hi_pts[np.argsort(hi_pts[:, 0])]; lo_pts = lo_pts[np.argsort(lo_pts[:, 0])]
    h1 = np.vstack([hi_pts[hi_pts[:, 0] < rc], lo_pts[lo_pts[:, 0] >= rc]])
    hp = np.vstack([lo_pts[lo_pts[:, 0] < rc], hi_pts[hi_pts[:, 0] >= rc]])
    h1 = h1[np.argsort(h1[:, 0])]; hp = hp[np.argsort(hp[:, 0])]
    for r, v in h1[::3]:
        rows9.append(('n_H1s_cm3', F[c], r, v))
    for r, v in hp[::3]:
        rows9.append(('n_Hplus_cm3', F[c], r, v))
    print('fig9c', c, 'crossing r =', round(rc, 3))
    for q, panel, k in (('T_K', 'fig23a', 0), ('v_kms', 'fig23b', 0), ('n_He_total_cm3', 'fig23c', 0),
                        ('n_Heplus_cm3', 'fig23d', 0), ('n_He23S_cm3', 'fig23d', 1)):
        for r, v in arr((panel, c, k))[::3]:
            rows23.append((q, F[c], r, v))
up, lo = arr(('fig12c', 'red', 0)), arr(('fig12c', 'red', 1))
rc = crossing(up, lo)
print('fig12c crossing', rc)
for r, v in lo[::3]:
    rows12.append(('n_H1s_cm3', 1.0, r, v))
for r, v in up[up[:, 0] >= rc + 0.05][::3]:
    rows12.append(('n_Hplus_cm3', 1.0, r, v))

common = """# Digitized by color tracking on a 300 dpi render of the publisher PDF
#   (pdftoppm -r 300, references/Yan 2022.pdf); axes calibrated on the frame
#   (x: log r/Rp, 1 and 10 on the frame edges) and on the y tick marks.
#   Scripts: benchmarks/yan2022_wasp52b/digitize/ (digitize_yan2022.py,
#   run_digitize.py, write_published.py); every 3rd pixel column kept.
# Reading precision (estimated): +-2 px of a 4 px wide line, i.e.
#   r +-0.5%%, %s.
#   Where curves of other F_XUV overlap the traced one, columns are missing
#   (the traced color is hidden), so coverage differs between curves.
# columns: quantity  F_XUV/F0  r/Rp  value
"""
def write(name, head, rows, prec):
    with open(os.path.join(OUT, name), 'w') as fh:
        fh.write(head + common % prec)
        for q, f, r, v in rows:
            fh.write('%-16s %4.2f  %.4f  %.5g\n' % (q, f, r, v))

write('published_yan2022_fig9_hhe92_8.txt',
"""# Yan et al. (2022), ApJ 936, 177, Figure 9 (page 15): H/He = 92/8, beta_m = 0.22
#   panel (a) temperature [K], panel (b) velocity [km/s], panel (c) number
#   densities [cm^-3] of H(1s) (solid) and H+ (dashed), against r/Rp.
#   Curves F0*0.5, F0*1, F0*2 (the three Table 4 points run here).
#   Panel (c): H(1s) and H+ of one color cross near 1.5 Rp; the two traced
#   envelopes are exchanged at the crossing radius; within about 0.05 Rp of
#   the crossing the assignment is uncertain.
""", rows9, "T +-35 K, v +-0.25 km/s, n +-0.03 dex")
write('published_yan2022_fig23_hhe98_2.txt',
"""# Yan et al. (2022), ApJ 936, 177, Figure 23 (page 29): H/He = 98/2, beta_m = 0.22
#   panel (a) temperature [K], (b) velocity [km/s], (c) total He number
#   density [cm^-3], (d) He+ (dotted) and He(2^3S) (dash-dotted) [cm^-3].
#   Curves F0*0.5, F0*1, F0*2. In (c) the F0*2 curve is hidden beyond 5.2 Rp.
#   Yan's He(2^3S) is a post-processed rate equilibrium without advection
#   (their section 2.3, Eq. 9), not a hydrodynamic variable.
""", rows23, "T +-35 K, v +-0.25 km/s, n +-0.03 dex (panel c +-0.03 dex, panel d +-0.035 dex)")
write('published_yan2022_fig12c_hhe98_2.txt',
"""# Yan et al. (2022), ApJ 936, 177, Figure 12(c) (page 18): F_XUV = F0, beta_m = 0.22,
#   H/He = 98/2 (red): H(1s) (solid) and H+ (dashed) number density [cm^-3].
#   The red H+ curve lies under the 99/1 and 99.5/0.5 curves inside the
#   H(1s)/H+ crossing (about 1.5 Rp) and is given only outside it.
""", rows12, "n +-0.03 dex")
print(len(rows9), len(rows23), len(rows12))
