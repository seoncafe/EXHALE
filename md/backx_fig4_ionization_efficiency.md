# Backx figure 4 read back: the H2 ionization efficiency behind `f_n`

*2026-09-03. For the owner of `scratchpad/e1h2` (section 151.10). The neutral
photodissociation fraction `f_n` of section 151 comes from table I of Chung,
Lee, Masuoka & Samson (1993), J. Chem. Phys. 99, 885, whose photoionization
yield column carries the footnote "Photoionization yields extracted from
Ref. 9" -- Backx, Wight & Van der Wiel (1976), J. Phys. B 9, 315, which
publishes that yield only as **figure 4**. This note reads it back
independently. No code was touched.*

**Verdict, in three parts.**

1. **Chung's table I reproduces Backx figure 4 to better than the figure can
   be read.** Over his twelve energies the mean difference in `f_n` is
   `-0.0011`, the RMS `0.0021`, the largest single difference `0.0056`,
   against a reading error of `0.0034`. Window minimum `37.4 +- 0.1 eV`
   against his `37.5`; depth `f_n = 0.074` against his `0.074`. **No
   correction to `f_n` is called for over 33-41 eV.**
2. **The figure did not have to be digitized at all.** Footnote (double
   dagger) of Backx's own table 1 states that its H2+ column was
   *constructed* from the smooth curve of figure 4, so the two published
   tables invert to give that curve back. The inversion has no pixels in it --
   and where the paper normalized, it returns `eta = 0.999 +- 0.002`
   (20-34 eV) and `1.002 +- 0.017` (44-70 eV).
3. **There is a second window Chung's table does not cover, and over part of
   its range it is the bigger one.** Backx's text says the efficiency is below
   unity in *two* regions. The other is 15.4-18 eV, where the table inversion
   gives `f_n = 0.26` at 16 eV, `0.098` at 16.5 and `0.047` at 17.0 -- against
   the `0.074` maximum of the 37 eV window. A `f_n` built from Chung's table
   alone is identically zero there.

Files: `bx_fig4_digitize.py`, `bx_fig4_analyze.py`, `bx_fig4_report.py`,
`bx_table_inversion.py`, `bx_fig4_fig.py`; `docs/backx_fig4_eta.dat` (the
hand-off table); `fig/backx_fig4_efficiency.png`; `dig/bx_overlay2.png` (the
trace drawn back on the scan, which is how it was checked).

---

## 1. What the figure is, in the paper's own words

Read off 600 and 400 dpi renders of the published pages, not off an OCR.

Section 3.3, step (iii), p. 321: "The total GOS for ionization, i.e.
`(1 + H+/H2+)` times the GOS of `H2+`, was divided by the corresponding
quantity for absorption to give an ionization efficiency (figure 4). This
relative efficiency `eta_i` was found to be constant above 18 eV, and for that
reason was normalized to unity in that region (except for a few points around
37 eV, which are lower than the adjacent regions)."

Step (iv): "Given this `eta_i` for photoionization, `f(0) (H+ + H2+)` was
obtained as the product of `eta_i` (**the smooth curve in figure 4**) and the
`f(0)` for absorption."

Figure 4 caption: "Ratio of the GOS for ionization and absorption of `H2` as
obtained from forward-scattering intensities (`K^2 ~ 0.01 au`). The absolute
scale is obtained by normalizing the average of the ratios at energy losses of
**18-34 eV and 44-70 eV** to unity."

Page 322: "There are **two** regions in the ionization efficiency curve where
`eta_i(E)` (figure 4) is lower than unity (**below 18 eV and around 37 eV**)."

Table 1, footnote (double dagger), p. 320: "`f(0)(H2+) = f(0)_abs eta_i /(1 +
H+/H2+)`, where for `eta_i` the smooth curve through the measured data points
of figure 4 was used and `H+/H2+` is given in table 2."

Four things follow. **The smooth curve, not the scatter, is what the paper
adopts**, so it is the smooth curve Chung read and the smooth curve this note
traces. **The scale is normalized on the 18-34 and 44-70 eV averages by
construction**, so "`eta = 1` outside the window" is a check of the reading,
not independent physical evidence. **Table 1's H2+ column is not independent
data** -- it is that same smooth curve, multiplied through. And **the paper
itself names two windows**, only one of which is in Chung's table.

## 2. The y axis is legible, and it checks

The first thing asked. It is: the axis carries labelled ticks at `0`, `0.4`,
`0.8` and `1.2`, whose ink is found beside the axis line at rows 1485.0,
1097.5, 703.5 and 313.0 of the 600 dpi crop. The three intervals are 387.5,
394.0 and 390.5 px, equal to `+-0.8` per cent. The x axis carries seven ticks
at 333.7 px spacing, labelled 10, 30, 50, 70.

Calibration by least squares **on the tick marks, not on the numerals** (which
sit a few pixels to their left):

```
E   = 0.029971 * col + 4.9528       residual RMS 0.061 eV
eta = -0.00102301 * row + 1.520452  residual RMS 0.0014
```

**Independent check of the y scale**: the paper draws a dashed line at
`eta = 1` across the panel, and the calibration puts it at **`eta = 1.0033`**.
That is the y scale good to 0.3 per cent, measured rather than asserted.

## 3. How the curve was traced, and what fought it

Three obstacles, each handled and each stated:

1. **Two inset boxes** (D2 and HD) cover the panel below `eta = 0.834`, over
   20.0-35.9 eV and 42.0-57.9 eV. The code locates their frames from their own
   long vertical edges rather than assuming coordinates, and masks them.
2. **The data circles sit on the curve.** Every line in the panel is 5-6 px at
   600 dpi and every filled circle about 16 px, so a column's ink groups sort
   by vertical extent: up to 9 px is line, thicker is a circle. That keeps the
   curve where a circle lies on it, which erasing a circle mask does not.
3. **Error bars.** A bar's vertical shaft is one tall group and the thickness
   test removes it; a horizontal cap is not, so the trace is walked outward
   from a clean seed column, taking at each step the thin group nearest a
   linear prediction from the last 25 columns, with the tolerance widening
   across a gap.

The data circles were extracted separately by erosion with a disk of radius 4
and 5 at binarization thresholds 100, 128 and 160 -- six runs, matched by
position. **Threshold and radius move a centroid by less than 0.001 in
`eta`**, so repeat-reading scatter is not the limiting term; the calibration
residual and the finite line width are. Half a drawn line is 3 px, hence

```
sigma_eta = 0.0034      sigma_E = 0.109 eV
```

## 4. The window: position, width, depth

From the smooth curve.

| quantity | figure 4, this work | Chung table I |
|---|---|---|
| minimum of `eta` | `0.926 +- 0.003` | `0.926` |
| at | `37.4 +- 0.1 eV` | `37.5 eV` |
| maximum `f_n = 1 - eta` | `0.074` | `0.074` |
| `eta < 0.99` | 30.7 - 40.7 eV (width 10.1 eV) | 33 - 41 eV tabulated |
| `eta < 0.98` | 30.7 - 40.3 eV (9.6 eV) | 34 - 40 eV |
| `eta < 0.95` | 35.6 - 38.6 eV (3.0 eV) | 35.5 - 38 eV |

Point by point, on Chung's own grid:

| E [eV] | `Y` Chung | `eta` fig. 4 | `f_n` Chung | `f_n` fig. 4 | difference |
|---|---|---|---|---|---|
| 33.0 | 0.992 | 0.994 | 0.008 | 0.006 | -0.002 |
| 34.0 | 0.978 | 0.981 | 0.022 | 0.019 | -0.003 |
| 35.0 | 0.960 | 0.959 | 0.040 | 0.041 | +0.001 |
| 35.5 | 0.951 | 0.950 | 0.049 | 0.050 | +0.001 |
| 36.0 | 0.942 | 0.942 | 0.058 | 0.058 | +0.000 |
| 36.5 | 0.934 | 0.936 | 0.066 | 0.064 | -0.002 |
| 37.0 | 0.929 | 0.929 | 0.071 | 0.071 | +0.000 |
| 37.5 | 0.926 | 0.926 | 0.074 | 0.074 | +0.000 |
| 38.0 | 0.927 | 0.927 | 0.073 | 0.073 | +0.000 |
| 39.0 | 0.940 | 0.946 | 0.060 | 0.054 | -0.006 |
| 40.0 | 0.973 | 0.974 | 0.027 | 0.026 | -0.001 |
| 41.0 | 0.992 | 0.994 | 0.008 | 0.006 | -0.002 |

mean `-0.0011`, RMS `0.0021`, largest `|difference|` `0.0056` at 39.0 eV --
`1.6 sigma`, and the only point outside `1 sigma`.

**Width.** Chung tabulates 33-41 eV; the figure puts the departure from unity
a little earlier, near 30.7 eV. That is not a disagreement: over 30.7-33 eV
the traced deficit is 0.004-0.008, one to two times the reading error, and the
same size as Chung's own end values (0.008 at both 33 and 41 eV). **The window
is 33-41 eV to within what either source resolves.**

## 5. The same curve without any pixels: inverting tables 1 and 2

Because table 1's H2+ column was built from figure 4's smooth curve
(section 1),

```
eta_i(E) = [ f(0)(H2+) / f(0)_abs ] * (1 + H+/H2+)
```

recovers it from two printed tables. `H+/H2+` comes from table 2 in units of
`1e-2`, and is zero below 18.0 eV because the dissociative-ionization
threshold is 18.08 eV and table 2 starts there.

**Where the paper normalized, the inversion returns unity**:

| range | `eta` from tables 1+2 | what the paper set |
|---|---|---|
| 20-34 eV | `0.999 +- 0.002` | 18-34 eV average = 1 |
| 44-70 eV | `1.002 +- 0.017` | 44-70 eV average = 1 |

That is a transcription check and a formula check in one, and the 0.002 spread
over 20-34 eV says both tables were transcribed correctly.

**Inside the window the inversion is coarse and must not be used point by
point**: it gives 0.998, 0.961, 0.928, 0.939, 0.869, 0.929, 0.988 at 34-40 eV,
scattered by `+-0.03` about the smooth curve. Two reasons, both in the paper:
table 1 is printed to three digits at values near 1.3, and footnote (dagger)
says "In part of the spectrum more data points were measured than given in the
table (see figure 2); for those regions the table lists an average over the
interval of the table" -- so the absorption column is an interval average while
`eta_i` is a point value, and the window is exactly where `eta_i` varies fast.
What the inversion does establish independently is that the window **exists**,
**sits at 36-39 eV** and is **7-13 per cent deep**.

## 6. Above 41 eV

The traced curve gives `eta = 1.001 +- 0.011` over 41-46 eV and
`1.001 +- 0.027` over 46-57 eV; the 24 data circles between 41 and 70 eV
average `0.993 +- 0.089` (standard error 0.018); the table inversion gives
`1.002 +- 0.017` over 44-70 eV. **So yes, `eta = 1` above 41 eV** -- with the
standing caveat that the paper normalized the 44-70 eV average to unity, so
this is consistency, not proof. What is not circular is the *shape*: the curve
returns to the normalization level immediately above 41 eV and shows no
further structure to 70 eV.

## 7. The window Chung's table does not cover: 15.4-18 eV

Backx names two sub-unity regions and this is the other one. The table
inversion, which needs no figure reading, gives

| E [eV] | `eta` | `f_n = 1 - eta` |
|---|---|---|
| 15.0 | 0.040 | 0.960 |
| 15.5 | 0.375 | 0.625 |
| 16.0 | 0.742 | 0.258 |
| 16.5 | 0.902 | 0.098 |
| 17.0 | 0.953 | 0.047 |
| 17.5 | 1.002 | -0.002 |
| 18.0 | 1.012 | -0.012 |

**Caveat on the first two rows.** The H2 ionization threshold is 15.43 eV, so
`eta` is exactly zero at 15.0 eV; the 0.040 and 0.375 are the paper's 0.5 eV
energy-loss resolution smearing the onset, not physics. From 16.0 eV upward
the values are away from the onset and are usable.

This is the same physical process as the 37 eV window -- absorption into a
state above the ionization threshold that dissociates neutrally rather than
ionizing, which the figure's insets and shaded areas call "the contribution of
superexcited states" -- and over 16-17 eV it is **larger than the 37 eV window
ever gets**. A `f_n` built only from Chung's table I is identically zero
there. Whether that matters is for the owner of section 151 to judge against
how `f_n` is used and how much of the band lands at 15.4-18 eV; it is flagged,
not decided.

## 8. Where the digitization is weak

* **46-57 eV**: the HD inset frame covers the curve and only six columns
  survive; the value there rests on the data circles and on the table
  inversion.
* **57.0-63.5 eV**: the trace follows an error-bar cap, reading 1.05-1.07
  where the drawn curve and every neighbouring datum are at 1.00. **Dropped**,
  and `docs/backx_fig4_eta.dat` says so.
* **Above 63.5 eV**: four columns reading 1.015.
* **27-28 eV and 43-45 eV**: the trace reads 1.013-1.015, 1.5 per cent above
  unity -- four times the quoted reading error but only twice the drawn line
  width, and the table inversion is at 1.000 there. Tracer drift, not
  structure.
* **Below 16.3 eV**: unusable. The curve is nearly parallel to the `eta` axis
  at the ionization onset and a column-wise trace cannot follow it. Use the
  table inversion of section 7 instead.

## 9. What this settles for section 151

1. `f_n` from Chung table I over 33-41 eV **is** Backx figure 4, to `0.006`
   at worst and `0.002` RMS. The 7.4 per cent maximum at 37.5 eV, the window
   position and the window width are all confirmed. **No correction.**
2. `eta = 1` above 41 eV holds, with the caveat that the figure's scale is
   defined that way; the shape adds no structure to 70 eV.
3. **`eta < 1` also over 15.4-18 eV**, reaching `f_n = 0.26` at 16 eV,
   `0.098` at 16.5 and `0.047` at 17.0 -- from the published tables, not from
   a figure. Chung's table does not carry it.
4. `docs/backx_fig4_eta.dat` carries the digitized curve on a 0.5 eV grid with
   its unreliable ranges marked, and the table inversion below 20 eV in its
   header, if section 151 wants a source for `f_n` outside 33-41 eV.
