"""Digitize Huang et al. (2023) Figures 9 (case A), 20 (case B) and 25 (case D): T(r), v(r).

Run from this directory: python3 digitize_fig9_fig20.py . ; writes fig9.npy, fig20.npy, fig25.npy."""
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from track import load, mask_color, trace

def fit(px, val, log):
    v = np.log10(val) if log else np.asarray(val, float)
    a, b = np.polyfit(px, v, 1)
    res = v - (a * np.asarray(px) + b)
    return (lambda p: 10 ** (a * p + b)) if log else (lambda p: a * p + b), np.abs(res).max()

def run(img, xt, xv, tl, tv, vr, vv, frame, seeds, exclude, outstem):
    im = load(img)
    X, ex = fit(xt, xv, True)
    T, et = fit(tl, tv, False)
    V, ev = fit(vr, vv, True)
    print(img, 'tick fit residuals: logR %.4f dex, T %.1f K, logv %.4f dex' % (ex, et, ev))
    c0, c1, r0, r1 = frame
    res = {}
    res_px = {}
    for name, (kind, sc, sr, jump, pick) in seeds.items():
        m = mask_color(im, kind)
        m[:r0, :] = False; m[r1:, :] = False
        m[:, :c0] = False; m[:, c1:] = False
        for (a, b, c, d) in exclude:
            m[c:d, a:b] = False
        cols, rows = trace(m, sc, sr, c0, c1, max_jump=jump, pick=pick)
        res_px[name] = (cols, rows)
        if name.startswith('v'):
            y = V(rows)
        else:
            y = T(rows)
        res[name] = (X(cols), y)
        print(name, len(cols), 'points, R %.3f-%.3f' % (X(cols).min(), X(cols).max()))
    from overlay import overlay
    overlay(img, list(res_px.values()), outstem + '_trace_check.png')
    return res

if __name__ == '__main__':
    out = sys.argv[1]
    # Figure 9: frame cols 173-998, rows 77-624
    r9 = run('huang2023_fig9.png',
             [251.0, 331.5, 417.5, 488.0, 574.0, 685.0, 771.5, 842.0, 953.0],
             [1.3, 1.6, 2.0, 2.4, 3.0, 4.0, 5.0, 6.0, 8.0],
             [127, 228, 330, 431, 533], [12000, 10000, 8000, 6000, 4000],
             [134.0, 269.0, 404.0, 539.5], [1e4, 1e3, 1e2, 1e1],
             (175, 997, 79, 622),
             {'T': ('black', 574, 213, 14, 'lowest'),
              'v': ('red', 574, 135.5, 10, None)},
             [(480, 700, 500, 610), (270, 300, 330, 600), (338, 372, 330, 560),
              (730, 758, 370, 600), (925, 953, 180, 395)],
             'huang2023_fig9')
    np.save(out + '/fig9.npy', {k: np.vstack(v) for k, v in r9.items()}, allow_pickle=True)

    r20 = run('huang2023_fig20.png',
              [220.0, 295.5, 377.0, 443.0, 524.5, 629.0, 710.5, 777.0, 881.5, 963.0],
              [1.3, 1.6, 2.0, 2.4, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0],
              [143.0, 232.0, 321.0, 410.0, 499.5, 588.5, 677.5],
              [12000, 10000, 8000, 6000, 4000, 2000, 0],
              [228.0, 338.5, 449.0, 559.0, 669.0], [1e5, 1e4, 1e3, 1e2, 1e1],
              (174, 996, 145, 688),
              {'T': ('black', 600, 453, 14, None), 'v': ('red', 600, 263.5, 10, None)},
              [(885, 996, 150, 220), (300, 326, 420, 630), (362, 401, 380, 650),
               (426, 456, 380, 610), (850, 880, 340, 540)],
              'huang2023_fig20')
    np.save(out + '/fig20.npy', {k: np.vstack(v) for k, v in r20.items()}, allow_pickle=True)
    r25 = run('huang2023_fig25.png',
              [198.0, 276.0, 360.0, 429.0, 513.0, 621.0, 705.0, 774.0, 882.5, 966.5],
              [1.3, 1.6, 2.0, 2.4, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0],
              [164.5, 247.0, 329.0, 411.0, 493.0, 575.0, 657.0],
              [12000, 10000, 8000, 6000, 4000, 2000, 0],
              [207.0, 330.5, 454.0, 577.0], [1e5, 1e4, 1e3, 1e2],
              (174, 996, 116, 659),
              {'T': ('black', 700, 272.0, 14, None), 'v': ('red', 705, 229.0, 10, None)},
              [(875, 996, 440, 525), (268, 299, 125, 380), (308, 338, 125, 380),
               (590, 623, 375, 645), (643, 675, 375, 610), (848, 880, 290, 522)],
              'huang2023_fig25')
    np.save(out + '/fig25.npy', {k: np.vstack(v) for k, v in r25.items()}, allow_pickle=True)
