# Fig. 2 digitization — superseded, kept as a cross-check

This folder digitizes **Fig. 2** of Cherubim et al. (2026) from the raster in
the published PDF. It was made before the authors' released spectrum was in
hand and is **not** the observational reference any more: use
`../LHS1140b_He10833_Fig3B_spectrum.csv`, reconstructed from the Zenodo
deposit the paper cites (DOI 10.5281/zenodo.20723095).

It is kept because it independently checks that reconstruction, and because
it records what a figure-only route can and cannot give.

## Agreement with the released data

The released file carries the same raw (pre-GP) spectrum that Fig. 2 plots,
so the two are directly comparable:

| | this digitization | released data |
|---|---|---|
| deepest point | -1.488 % | -1.497 % |
| its wavelength | 10833.426 A | 10833.416 A |

0.009 % in depth and 0.010 A in wavelength.

## Why it is a band and not a list of points

The released spectrum samples every 0.0361 A. The plot markers are 24 px
across, which is 0.28 A on that axis, so consecutive points overlap
eight-fold and the union of the discs does not determine their centres. The
digitizer therefore reports, at each column, the interval of marker centres
consistent with the ink, and that interval widened by the error bars. That
was the right product for a figure sampled this way, and the agreement above
is the evidence it worked.

## Files

| file | what it is |
|---|---|
| `lhs1140b_fig2_band.txt` | the band, in both sign conventions |
| `digitize_fig2.py` | the script, run from this directory |
| `fig2/p11-000.png` | Fig. 2 as extracted from the PDF (`pdfimages -f 11`) |
| `band_check.png` | the band drawn back over the figure |

Axis calibration was checked against the figure's own annotation -- the three
blue lines mark the He I triplet at vacuum 10832.057, 10833.217, 10833.306 A,
recovered to within 0.015 A -- and the script stops if that check fails.
