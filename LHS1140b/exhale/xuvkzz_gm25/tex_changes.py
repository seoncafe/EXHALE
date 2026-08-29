#!/usr/bin/env python3
"""The list of edits the re-solved numbers imply for
docs/lhs1140b_exhale_vs_pwinds.tex.

Nothing here edits the memo.  Each item quotes the current text straight out
of the file at the line range given, so the quotation cannot drift from what
is actually written there, and states what the re-measured value is.  Line
numbers are those of the file as read at the time this ran; check the quoted
text, not the number, when applying.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX_ROOT = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
TEX = os.path.join(EX_ROOT, 'docs', 'lhs1140b_exhale_vs_pwinds.tex')
sys.path.insert(0, os.path.join(HERE, '..', 'crossings_gm25'))
import measure as M                                            # noqa: E402
os.chdir(HERE)

LINES = open(TEX).readlines()


def quote(a, b):
    """The memo's own text, lines a..b inclusive, 1-based."""
    return ''.join('    | ' + LINES[i-1] for i in range(a, b+1))


def item(n, title, a, b, becomes):
    print('-'*70)
    print('%s  %s   [lines %d-%d]' % (n, title, a, b))
    print('  now reads:')
    print(quote(a, b), end='')
    print('  becomes:')
    for l in becomes.strip('\n').split('\n'):
        print('    ' + l)
    print()


def red(d):
    return M.metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))['red_depth']


def ew(d):
    return M.red_ew(os.path.join(d, 'tpm_He10830.txt'))


def mdot(d):
    return M.wind(d)[0]


def cross(xs, ds, level=0.6):
    xs, ds = np.asarray(xs, float), np.asarray(ds, float)
    o = np.argsort(xs); xs, ds = xs[o], ds[o]
    for i in range(len(xs)-1):
        if (ds[i]-level)*(ds[i+1]-level) <= 0 and ds[i+1] != ds[i]:
            return xs[i] + (level-ds[i])*(xs[i+1]-xs[i])/(ds[i+1]-ds[i])
    s = (ds[-1]-ds[-2])/(xs[-1]-xs[-2])
    return xs[-1] + (level-ds[-1])/s


print('='*70)
print('WHAT THIS IMPLIES FOR docs/lhs1140b_exhale_vs_pwinds.tex')
print('='*70)
print("""
Not applied here -- this directory does not edit the memo.  Every entry
below gives the memo's current text and the re-measured replacement, so one
person can work through the list.  Three things are NOT on it and are worth
saying explicitly:

  * Table tab:kzzformmatch (elemental ratios and carriers at the match) is
    chemistry only.  The Penning coefficient is a wind-side rate and the
    photochemical columns were not re-solved, so every entry stands.
  * The composition statements of Sect. sec:outerinherit -- He/H at the
    match to seven digits, the cold-trap O/H of the two 12.01 arms -- are
    likewise chemistry and stand.
  * Table tab:validity and its exobase range 18.8-25.9 R_p are a different
    block of runs; its caption already carries the retired-coefficient note.
    The ladder arms' own exobases are re-measured here and quoted where the
    outer-region argument uses them.
""")

FRAC = [1.00, 0.33, 0.30, 0.25, 0.20, 0.15, 0.10, 0.01]
TAG = {1.00: '1p00', 0.33: '0p33', 0.30: '0p30', 0.25: '0p25',
       0.20: '0p20', 0.15: '0p15', 0.10: '0p10', 0.01: '0p01'}
S = {f: 'S'+TAG[f] for f in FRAC}
C = {f: 'C'+TAG[f] for f in FRAC if f != 0.30}

print()
print('#'*70)
print('# BLOCK 1: the XUV grid   (Sect. sec:xuvgrid, Table tab:xuvgrid)')
print('#'*70)
print()

rows = []
for f in FRAC:
    sd, cd = S[f], C.get(f)
    rows.append((f, red(sd), ew(sd), mdot(sd),
                 red(cd) if cd else None, ew(cd) if cd else None,
                 mdot(cd) if cd else None))

sx = [f for f in FRAC if f < 1.0]
sc = cross(sx, [red(S[f]) for f in sx])
cx = [f for f in FRAC if f < 1.0 and f in C]
cc = cross(cx, [red(C[f]) for f in cx])

def fm(x):
    """Enough digits to keep the small entries significant, as the memo does."""
    return '%.4f' % x if x < 0.01 else '%.3f' % x


tab = []
for f, sr, se, sm, cr, ce, cm in rows:
    bs = '\\mathbf{%s}' % fm(sr) if abs(f - 0.30) < 1e-9 else fm(sr)
    if cr is None:
        tab.append('$%.2f$ & $%s$ & $%s$ & $%.2f$ & --- & --- & --- \\\\'
                   % (f, bs, fm(se), sm))
    else:
        bc = '\\mathbf{%s}' % fm(cr) if abs(f - 0.33) < 1e-9 else fm(cr)
        tab.append('$%.2f$ & $%s$ & $%s$ & $%.2f$ & $%s$ & $%s$ & $%.2f$ \\\\'
                   % (f, bs, fm(se), sm, bc, fm(ce), cm))

item('1.1', 'Table tab:xuvgrid, every data row', 3522, 3529,
     'The eight rows re-measured (bold still marks the row at which each\n'
     'model is nearest the 0.6 per cent limit from below):\n\n'
     + '\n'.join(tab) + """

Note the bold moves: on the current coefficient the scalar model's last
row above 0.6 per cent is still 0.33 and its first below is 0.30, so the
bold stays where it is; in the closed model no row reaches 0.6 per cent
at all, and the bold on its 0.33 row marks the top of the grid rather
than a crossing.""")

item('1.2', 'the crossing sentence -- THE VERDICT CHANGES', 3534, 3544,
     """The crossing moves UP in both models, not down, and the two no longer
land in the same place.  Re-measured, linear in both variables, on the
scaled grid alone (the fiducial is the model, not a grid point), which is
the reading that reproduces the memo's own 0.31 and 0.36:

  scalar model   0.308 -> %.3f   (still bracketed between 0.30 and 0.33)
  flux-closed    0.360 -> %.3f   (still above the grid, top chord 0.25-0.33)

so the replacement text has to say that the scalar model crosses at about
0.32 and the flux-closed one at about 0.42, a factor %.2f apart rather than
1.2, and that on F_XUV = 33 erg/cm2/s the 2025 epoch needs about
%.0f-%.0f erg/cm2/s or less.  The sentence "The 0.6 per cent limit is met
at about a third of the fiducial XUV in both models" no longer holds for
the flux-closed model.  The factors 4.6 in reservoir and 2.5 in escape
rate are unchanged.""" % (sc, cc, cc/sc, 33*sc, 33*cc))

s33, s10 = red(S[0.33]), red(S[0.10])
c33, c10 = red(C[0.33]), red(C[0.10])
c01, c10b = red(C[0.01]), red(C[0.10])
item('1.3', 'the two caveats: the slopes, and the 0.01 point', 3546, 3554,
     """The contrast between the two slopes survives and grows slightly.  Over
0.33 to 0.10 the depth falls by a factor %.0f in the scalar model against
%.1f in the closed one (was 290 against 3.7).

The non-monotonic lowest closed point: %.3f per cent at 0.01 against
%.3f at 0.10, so it %s.""" % (s33/s10, c33/c10, c01, c10b,
                              'is still non-monotonic and the caveat stands'
                              if c01 > c10b else
                              'is now monotonic with the 0.10 point and the '
                              'caveat has to be dropped or restated'))

item('1.4', 'the retired-coefficient note, and the two model EWs', 3495, 3507,
     """Both grids are now on the current coefficient.  The two models'
fiducial equivalent widths become EW = %.3f (scalar, He/H = 2.13) and
EW = %.3f (flux-closed, He/H = 9.71), from EW = 1.129 and 1.040.  The
prediction the note makes -- "every depth in the grid rises by roughly a
fifth, so each model would cross the 0.6 per cent limit at a slightly
lower XUV fraction" -- is contradicted by the measurement and must be
removed, not merely re-pointed: only the fiducial rises (by %.0f and %.0f
per cent); every scaled point falls, by up to a factor %.1f, and both
crossings move to a HIGHER XUV fraction.  Runs:
LHS1140b/exhale/xuvkzz_gm25/S*, C*.""" % (
         ew(S[1.00]), ew(C[1.00]),
         100*(red(S[1.00])/4.152 - 1), 100*(red(C[1.00])/3.491 - 1),
         max(0.726/red(S[0.33]), 0.551/red(S[0.30]), 0.311/red(S[0.25]),
             0.568/red(C[0.33]), 0.481/red(C[0.25]))))

item('1.5', 'the summary bullet', 4075, 4088,
     """The bullet heading "The 2025 non-detection sits at about $0.3$ of the
fiducial XUV" becomes about 0.32 in the scalar-base model and about 0.42
in the flux-closed one, i.e. roughly %.0f-%.0f erg/cm2/s or less, and the
claim that the two models agree on the crossing weakens from a factor 1.2
to a factor %.2f.  The closing sentence -- "a coefficient that raises every
equivalent width by about a fifth moves the crossing to a slightly lower
XUV fraction in both models" -- is wrong in direction and must go.""" %
     (33*sc, 33*cc, cc/sc))

print()
print('#'*70)
print('# BLOCK 2: the K_zz pressure form   (Sect. sec:kzzform, tab:kzzform)')
print('#'*70)
print()

import math
SLOPE = 0.3819
KA = {'11.11': dict(rec='L11p1', ref='P11_ref', A='P11_A', B='P11_B'),
      '2.09': dict(rec='L2p09', ref='P2_ref', A='P2_A', B='P2_B')}
LABTEX = {'rec': 'recorded, other seed  ', 'ref': 'reference, $K_{zz} = 10^{9}$',
          'A': 'A, chemistry only     ', 'B': 'B, chemistry and wind '}


def heh2rp(d):
    a = np.loadtxt(os.path.join(d, 'output', 'element_flux_profile.txt'),
                   usecols=(1, 2))
    return float(np.interp(2.0, a[:, 0], (a[:, 1]/4.0)/(1.0 - a[:, 1])))


ktab, kderived = [], {}
for res in ('11.11', '2.09'):
    ktab.append('\\multicolumn{7}{l}{\\emph{reservoir} $\\mathrm{He/H} = %s$} \\\\'
                % res)
    e = {}
    for k in ('rec', 'ref', 'A', 'B'):
        d = KA[res][k]
        m = M.metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
        w = ew(d)
        e[k] = w
        ktab.append('%s & $%.4f$ & $%.4f$ & $%.4f$ & $%.3f$ & $%.4f$ & $%.3f$ \\\\'
                    % (LABTEX[k], m['red_depth'], m['fwhm_A'], w,
                       w/M.EW_OBS, heh2rp(d), mdot(d)))
    kderived[res] = e

item('2.1', 'Table tab:kzzform, both blocks of rows', 2642, 2652,
     'The eight rows re-measured:\n\n' + '\n'.join(ktab))

d11, d2 = kderived['11.11'], kderived['2.09']


def shift(e):
    return 11.73*math.exp(-math.log(e)/SLOPE)


cr = [shift(d11['A']/d11['ref']), shift(d11['B']/d11['ref']),
      shift(d2['A']/d2['ref']), shift(d2['B']/d2['ref'])]
pc = [100*(x/11.73 - 1) for x in cr]
sp11 = 100*(d11['rec']/d11['ref'] - 1)
sp2 = 100*(d2['rec']/d2['ref'] - 1)
eff = [100*(d11['A']/d11['ref'] - 1), 100*(d11['B']/d11['ref'] - 1),
       100*(d2['A']/d2['ref'] - 1), 100*(d2['B']/d2['ref'] - 1)]

item('2.2', 'the crossing paragraph, and the path spread it is read against',
     2717, 2731,
     """The four crossings become %.2f and %.2f on the 11.11 arm and %.2f and
%.2f on the 2.09 one, i.e. %+.2f to %+.2f in He/H (was -0.19 to +0.05) and
%+.1f to %+.1f per cent in relative terms (was -1.6 to +0.4).

The path spread they are read against becomes %.2f per cent on the 11.11 arm
and %.2f per cent on the 2.09 one (was 1.5 and 0.31), against the
%.2f-%.2f per cent the eddy form produces (was 0.10-0.61).

The verdict is unchanged and the margin is if anything the same: the largest
effect of the eddy form, %.2f per cent, is still under the largest path
spread, %.2f per cent, so "the signal is buried in the spread of the
procedure that would have to measure it" stands.

NOTE on the memo's own "recorded, other seed" row at 11.11 (red 3.5998,
EW 1.0853).  That row was read from flux_closure/heh11p1/k01, and that
directory has since been re-run in place (mtime 2026-08-28 19:04), so it no
longer holds the numbers the memo printed.  The row above is the re-solved
value; the memo's printed one cannot be reproduced from the tree.""" % (
         cr[0], cr[1], cr[2], cr[3],
         min(cr)-11.73, max(cr)-11.73, min(pc), max(pc),
         abs(sp11), abs(sp2), min(abs(x) for x in eff),
         max(abs(x) for x in eff), max(abs(x) for x in eff),
         max(abs(sp11), abs(sp2))))

item('2.3', 'the retired-coefficient note', 2621, 2624,
     """The test has now been repeated on the current coefficient, run by run,
in LHS1140b/exhale/xuvkzz_gm25/{P11_ref,P11_A,P11_B,P2_ref,P2_A,P2_B} with
the two "recorded, other seed" arms in L11p1 and L2p09.  The reasoning the
sentence gives -- that the coefficient is common to both sides of every
entry -- is confirmed rather than overturned: every ratio the test reports
moves by less than the workflow's own path spread.  The note becomes a
statement that the test was repeated and what it gave.""")

print()
print('#'*70)
print('# BLOCK 3: the outer region a converged wind inherits')
print('#            (Sect. sec:outerinherit, Table tab:outerinherit)')
print('#'*70)
print()

sys.path.insert(0, os.path.join(EX_ROOT, 'src', 'utils'))
import collisional_validity as CV                              # noqa: E402

LAD = [('$2.09$     ', 'L2p09'), ('$3.0$      ', 'L3p0'),
       ('$5.0$      ', 'L5p0'), ('$8.0$      ', 'L8p0'),
       ('$9.7$      ', 'L9p7'), ('$10.3$     ', 'L10p3'),
       ('$10.7$     ', 'L10p7'), ('$11.1$     ', 'L11p1'),
       ('$12.0$     ', 'L12p0ctl')]


def outer3(d):
    rd, fo, t12 = M.outer(d)
    return rd, t12, fo


ltab, shares, rdrops = [], [], []
for lab, d in LAD:
    rd, t12, fo = outer3(d)
    ltab.append('%s & $%.2f$ & $%.0f$ & $%.4f$ \\\\' % (lab, rd, t12, fo))
    shares.append(fo); rdrops.append(rd)
jd, jt, jf = outer3('L12p0jump')
ltab.append('\\midrule')
ltab.append('$12.0$ discarded & $\\mathbf{%.2f}$ & $%.0f$ & $\\mathbf{%.4f}$ \\\\'
            % (jd, jt, jf))

item('3.1', 'Table tab:outerinherit, every cell', 3423, 3434,
     'The ten rows re-measured:\n\n' + '\n'.join(ltab))

ctl_rd, ctl_t12, ctl_fo = outer3('L12p0ctl')
item('3.2', 'the caption\'s retired-coefficient note', 3412, 3417,
     """These arms are now on the current coefficient
(LHS1140b/exhale/xuvkzz_gm25/L*).  The comparison the note draws with the
single-parent ladder of Table tab:closureladder still reads r_drop = 20.46
falling to 19.47 R_p and an outer column share of 0.0091-0.0095 -- those
come from exhale/crossings_pc090, which is already on the current binary --
but the range quoted for THIS ladder changes with the table above, to
%.4f-%.4f.""" % (min(shares), max(shares)))

ew_j, ew_c = ew('L12p0jump'), ew('L12p0ctl')
rd_j, rd_c = red('L12p0jump'), red('L12p0ctl')


def resid(d):
    return M.wind(d)[2]


item('3.3', 'the four bullets that separate the two 12.01 solutions',
     3387, 3399,
     """Two of the four bullets move.

"Both converged": the achieved residuals become ||R|| = %.3fe-4 (jump-seeded)
and %.3fe-4 (control), both at info = 0 and both inside the 4.0e-4 target --
so the bullet's point is unchanged, but the numbers 2.529e-4 and 3.864e-4
are replaced, and note the control now reaches the SMALLER residual, which
reverses the memo's remark elsewhere that "the one that is wrong reports the
smaller residual of the two" (see item 3.5).

"Different line": EW = %.4f against %.4f %%A and a red depth of %.3f against
%.3f per cent, a %.0f per cent difference in equivalent width at fixed
composition (was 17 per cent).

The other two bullets -- same reservoir to seven digits, same cold-trap O/H
-- are chemistry and stand unchanged.""" % (
         resid('L12p0jump')*1e4, resid('L12p0ctl')*1e4,
         ew_j, ew_c, rd_j, rd_c, 100*(ew_j/ew_c - 1)))

item('3.4', 'the paragraph reading the table', 3437, 3448,
     """The trend statements survive.  r_drop still moves inward monotonically
as the reservoir rises, from %.1f to %.1f R_p, and the outer share of the
metastable column still stays in a narrow band across the whole ladder,
now %.4f-%.4f (was 0.0104-0.0119).

The discarded solution still departs from both: %.1f R_p where its own
reservoir gives %.1f, and %.2f times the outer column share of its control
(%.2f times the largest share any other rung reaches), so "twice the outer
column share of any other rung" becomes "about twice".

ONE CLAUSE IN THIS PARAGRAPH NO LONGER HOLDS.  "while its $T(12\\,R_p)$ is
the same as its control's to $0.3$ per cent" is now false: T(12 R_p) is
%.0f K against the control's %.0f K, a difference of %.1f per cent.  The
sentence's purpose -- that the two agree on the base and differ only
outside -- is still carried by the composition bullets and by the profile
ratio in item 3.6, but this particular quantity no longer makes the point
and has to be dropped or replaced.""" % (
         max(rdrops), min(rdrops), min(shares), max(shares),
         jd, ctl_rd, jf/ctl_fo, jf/max(shares), jt, ctl_t12,
         100*abs(jt/ctl_t12 - 1)))

item('3.5', 'the mechanism paragraph', 3449, 3455,
     """One clause changes.  "both report info = 0, and the one that is wrong
reports the smaller residual of the two" -- on the current coefficient the
jump-seeded solution reports the LARGER residual (%.3fe-4 against %.3fe-4).
Both still report info = 0 and both are still inside the target, so the
point that ||R|| cannot tell the two states apart stands; the rhetorical
sharpening about which one reports the smaller residual does not, and the
same clause appears in docs/Update_EXHALE.md Sect. 83.""" % (
         resid('L12p0jump')*1e4, resid('L12p0ctl')*1e4))

exo_c = CV.collisional_diagnosis('L12p0ctl', adv=True)['r_exobase']
exo_j = CV.collisional_diagnosis('L12p0jump', adv=True)['r_exobase']


def prof(d):
    o = os.path.join(d, 'output')
    hn = M.cols(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    sn = M.cols(os.path.join(o, 'Ion_species_adv.txt'))
    h = np.loadtxt(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    a = np.loadtxt(os.path.join(o, 'Ion_species_adv.txt'))
    return h[:, 0], h[:, hn.index('T[K]')], a[:, sn.index('HeITR')]


r, Tc, nc = prof('L12p0ctl')
_, Tj, nj = prof('L12p0jump')
rat = nj/np.maximum(nc, 1e-300)
first10 = r[np.where((r > 2.0) & (abs(rat - 1.0) > 0.10))[0][0]]
first50 = r[np.where((r > 2.0) & (abs(rat - 1.0) > 0.50))[0][0]]

# the same measurement on the stored pair, to separate what the coefficient
# changed from what the memo already had wrong
ro, _To, nco = prof('../flux_closure/heh12_ctl/k02')
_, _Tjo, njo = prof('../flux_closure/heh12/k04')
rato = njo/np.maximum(nco, 1e-300)
old10 = ro[np.where((ro > 2.0) & (abs(rato - 1.0) > 0.10))[0][0]]
old50 = ro[np.where((ro > 2.0) & (abs(rato - 1.0) > 0.50))[0][0]]

exo_all_new, exo_all_old = [], []
for _l, _d in LAD + [('12.0 jump', 'L12p0jump')]:
    v = CV.collisional_diagnosis(_d, adv=True)['r_exobase']
    if v is not None:
        exo_all_new.append(v)
SRC_LAD = ['../flux_closure/' + x for x in
           ('hi/k06', 'heh3/k03', 'heh5/k04', 'heh8/k04', 'heh9p7/k03',
            'heh10p3/k01', 'heh10p7/k01', 'heh11p1/k01', 'heh12_ctl/k02',
            'heh12/k04')]
for _d in SRC_LAD:
    v = CV.collisional_diagnosis(_d, adv=True)['r_exobase']
    if v is not None:
        exo_all_old.append(v)

item('3.6', 'the physical judgement -- THE ARGUMENT HOLDS, THE RADII DO NOT',
     3457, 3467,
     """The judgement survives on the current coefficient.  The exobase does not
move with the coefficient: across the ladder it sits at %.1f-%.1f R_p
re-solved against %.1f-%.1f R_p stored, and for the two 12.01 solutions
themselves it is %.2f R_p (control) and %.2f R_p (jump-seeded).  The bulk of
the region in which they differ is collisionless in both, so the ground for
discarding the jump-seeded rung remains continuity of the trend and not a
physical criterion, and the systematic it puts on the whole ladder's absolute
equivalent widths stands.

BUT THE TWO RADII IN THIS SENTENCE ARE BOTH WRONG, AND WERE ALREADY WRONG
BEFORE THE COEFFICIENT CHANGED.  Measured on the metastable profiles:

  * "The $20$--$30\\,R_p$ region in which the two solutions differ".  The two
    agree to %.1f per cent out to %.1f R_p and then separate abruptly: the
    departure passes 10 per cent at %.2f R_p and 50 per cent at %.2f R_p,
    reaching a factor %.1f.  On the stored pair the same measurement gives
    %.2f and %.2f R_p, so this is a pre-existing inaccuracy in the memo and
    not something the re-solve introduced.  The region is %.0f-30 R_p.

  * "where Sect.~\\ref{sec:validity} places the exobase, $18.8$--$25.9\\,R_p$".
    That range is the three-case range of Table tab:validity, whose cases are
    two non-closure runs and the 2.09 closure arm; it is not the range for
    these arms.  For the pair under discussion the exobase is %.1f-%.1f R_p.

Putting the corrected numbers together sharpens the sentence in one place and
qualifies it in another.  The differing region starts at about %.0f R_p, the
jump-seeded solution's own exobase is %.2f R_p, so for THAT solution -- the
one being discarded -- the whole differing region is outside its exobase.  For
its control the exobase is %.2f R_p, so the innermost %.1f R_p of the
differing region is nominally still collisional there.  The honest form of the
claim is therefore that the region is collisionless in the discarded solution
and marginal in the control over a shell about %.1f R_p deep, rather than that
neither solution is valid anywhere they differ.""" % (
         min(exo_all_new), max(exo_all_new), min(exo_all_old), max(exo_all_old),
         exo_c, exo_j,
         100*max(abs(rat[(r > 2) & (r < 11.9)] - 1)), 11.9,
         first10, first50, rat[r > 2].max(), old10, old50, first10,
         min(exo_c, exo_j), max(exo_c, exo_j),
         first50, exo_j, exo_c, exo_c - first50, exo_c - first50))

item('3.7', 'the two other places the 17 per cent appears', 3860, 3864,
     """"Its equivalent width is $17$ per cent higher for that reason alone"
becomes %.0f per cent.  The same figure appears again at line 4039,
"differ by $17$ per cent in equivalent width", and in
docs/Update_EXHALE.md Sect. 83, whose table of the two solutions carries
the full set of numbers replaced by item 3.3.""" % (100*(ew_j/ew_c - 1)))

print('-'*70)
print('END OF LIST')
