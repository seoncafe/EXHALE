"""Write the published_*.txt profile files of this benchmark from the
digitized tracks (fig9.npy, fig12.npy, fig20.npy, fig25.npy).

Each curve is resampled on a uniform grid in log10 r (step 0.004 dex) by
linear interpolation, of T in r and of log10 v / log10 n in log10 r, inside
the radius range the tracer covered; NaN outside it.
Run from this directory: python3 write_published_profiles.py
"""
import numpy as np

HDR = {
    'A': dict(fig='Figure 9 (case A, after three iterations; black solid T, red dashed v_z)',
              unit='R_p', unit_RJ=1.766, R0=2.058,
              note='R_p = 1.766 R_J (Huang Table 1); EXHALE R0 = 2.058 R_J (1 microbar, Table 3 log g = 2.84)'),
    'B': dict(fig='Figure 20 (case B, substellar direction; black solid T, red dashed v_z)',
              unit='R_px', unit_RJ=1.945, R0=2.418,
              note='R_px = 1.945 R_J (Huang Table 2, canonical); EXHALE R0 = 2.418 R_J (Table 3 log g = 2.70)'),
    'D': dict(fig='Figure 25 (case D, substellar streamline; black solid T, red dashed v_z)',
              unit='R_px', unit_RJ=1.959, R0=2.581,
              note='R_px = 1.959 R_J (Huang Table 2, better-match set); EXHALE R0 = 2.581 R_J (Table 3 log g = 2.62)'),
}
PREC = ('# reading precision (estimated): radius 0.002 dex (tick-fit residual <= 0.0003 dex plus\n'
        '#   line width); T +-60 K (line width ~3 px at ~20 K/px); log10 v +-0.02 dex\n'
        '#   (~135 px/dex, 2-3 px line). The red v curve is dashed, so gaps are bridged by\n'
        '#   the interpolation. Case B v stops at 8.05 R_px, where the curve enters the legend.\n')


def resample(x, y, grid, logy):
    o = np.argsort(x)
    x, y = np.log10(x[o]), (np.log10(y[o]) if logy else y[o])
    _, u = np.unique(x, return_index=True)
    x, y = x[u], y[u]
    g = np.log10(grid)
    out = np.interp(g, x, y)
    out[(g < x[0]) | (g > x[-1])] = np.nan
    return 10 ** out if logy else out


def tv_file(case, npy, out):
    h = HDR[case]
    d = np.load(npy, allow_pickle=True).item()
    xT, yT = d['T']; xv, yv = d['v']
    lo, hi = min(xT.min(), xv.min()), max(xT.max(), xv.max())
    grid = 10 ** np.arange(np.log10(lo), np.log10(hi) + 1e-9, 0.004)
    T = resample(xT, yT, grid, False)
    v = resample(xv, yv, grid, True) * 100.0          # m/s -> cm/s
    with open(out, 'w') as f:
        f.write('# Huang, Koskinen, Lavvas & Fossati 2023, ApJ 951, 123 (2023ApJ...951..123H)\n')
        f.write('# source: %s\n' % h['fig'])
        f.write('# axes as printed: x = R/%s (log scale), left y = T [K] (linear), right y = v_z [m s^-1] (log)\n' % h['unit'])
        f.write('# %s\n' % h['note'])
        f.write('# digitization: 300 dpi render (pdftoppm) of the published PDF; axis calibration by\n'
                '#   least-squares fit to the tick marks found in the image; curves traced by color\n'
                '#   (digitization/track.py, digitize_fig9_fig20.py); trace overlays in\n'
                '#   digitization/*_trace_check.png. MEASURED by the digitizer, values READ from the figure.\n')
        f.write(PREC)
        f.write('# columns: r/%s  r/R_J  r/R0_EXHALE  T[K]  v[cm s^-1]\n' % h['unit'])
        for g, t, vv in zip(grid, T, v):
            f.write('%9.5f %9.5f %9.5f %10.1f %12.5e\n' % (g, g * h['unit_RJ'], g * h['unit_RJ'] / h['R0'], t, vv))


def n_file(npy, out):
    d = np.load(npy, allow_pickle=True).item()
    names = ['e_and_p', 'H2p', 'MgI', 'MgII', 'FeI', 'FeII', 'CaII']
    lo = min(d[k][0].min() for k in names); hi = max(d[k][0].max() for k in names)
    grid = 10 ** np.arange(np.log10(lo), np.log10(hi) + 1e-9, 0.004)
    cols = [resample(d[k][0], d[k][1], grid, True) for k in names]
    with open(out, 'w') as f:
        f.write('# Huang, Koskinen, Lavvas & Fossati 2023, ApJ 951, 123 (2023ApJ...951..123H)\n'
                '# source: Figure 12 (case A number densities, upper atmosphere part r >= 1.17 R_p only;\n'
                '#   the green-shaded photochemical lower atmosphere r < 1.165 R_p is not digitized)\n'
                '# axes as printed: x = R/R_p (log), y = n [cm^-3] (log, 1e-2 to 1e12)\n'
                '# R_p = 1.766 R_J; EXHALE R0 = 2.058 R_J\n'
                '# curves: e (olive solid; the proton curve, olive dashed, lies on it within the line\n'
                '#   width for r > 1.3 R_p and below it for r < 1.3 R_p, where this column is n_e),\n'
                '#   H(2p) black solid (r >= 1.33 R_p only: inside it the curve crosses Na, which the\n'
                '#   tracer cannot separate), Mg I red solid, Mg II red dashed, Fe I blue solid (to\n'
                '#   3.77 R_p, where it leaves the axis), Fe II blue dashed, Ca II magenta dashed.\n'
                '#   Not digitized: p separately, Ca I, Na, K, H(2p) below 1.33 R_p. Neutral H is not\n'
                '#   plotted in the paper, so the H ionization fraction is not available from it.\n'
                '# digitization: as the T/v files (digitization/digitize_fig12.py); reading precision\n'
                '#   +-0.05 dex in n (44 px/dex, 3 px lines), 0.002 dex in r.\n'
                '# columns: r/R_p  r/R0_EXHALE  n_e  n_H(2p)  n_MgI  n_MgII  n_FeI  n_FeII  n_CaII  [cm^-3]\n')
        for i, g in enumerate(grid):
            f.write('%9.5f %9.5f ' % (g, g * 1.766 / 2.058) + ' '.join('%11.4e' % c[i] for c in cols) + '\n')


if __name__ == '__main__':
    tv_file('A', 'fig9.npy', '../published_caseA_fig9_T_v.txt')
    tv_file('B', 'fig20.npy', '../published_caseB_fig20_T_v.txt')
    tv_file('D', 'fig25.npy', '../published_caseD_fig25_T_v.txt')
    n_file('fig12.npy', '../published_caseA_fig12_densities.txt')
