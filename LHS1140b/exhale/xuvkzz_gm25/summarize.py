#!/usr/bin/env python3
"""The measurement record of the three blocks re-solved in this directory.

Every arm here is a stored solution re-solved on the current binary with its
lower atmosphere held fixed (`resolve_wind.sh`), so the composition, the XUV
scaling and K_zz are the source arm's and only the wind and the line are new.
The rows are read out of each arm's own directory by ../crossings_gm25/
measure.py; the stored arms this replaces are read alongside, unchanged, so
the two coefficients can be compared entry by entry.
"""
import math, os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'crossings_gm25'))
import measure as M

HERE = os.path.dirname(os.path.abspath(__file__))
EX_ROOT = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
os.chdir(HERE)


def rows(pairs):
    print(M.HDR)
    return {lab: M.row(lab, d) for lab, d in pairs}


def cross_linear(xs, ds, level):
    """Where the depth crosses `level`, linear in both, on the top chord."""
    xs, ds = np.asarray(xs, float), np.asarray(ds, float)
    o = np.argsort(xs)
    xs, ds = xs[o], ds[o]
    for i in range(len(xs)-1):
        if (ds[i]-level)*(ds[i+1]-level) <= 0 and ds[i+1] != ds[i]:
            return (xs[i] + (level-ds[i])*(xs[i+1]-xs[i])/(ds[i+1]-ds[i]),
                    'bracketed [%g, %g]' % (xs[i], xs[i+1]))
    # not bracketed: continue the top chord
    s = (ds[-1]-ds[-2])/(xs[-1]-xs[-2])
    return xs[-1] + (level-ds[-1])/s, ('above the grid, top chord [%g, %g]'
                                       % (xs[-2], xs[-1]))


def wind_extra(d):
    """He/H at 2 and 10 R_p, and the peak temperature."""
    f = os.path.join(d, 'output', 'element_flux_profile.txt')
    h = np.loadtxt(os.path.join(d, 'output', 'Hydro_ioniz.txt'))
    if not os.path.isfile(f):
        return float('nan'), float('nan'), h[:, 4].max()
    a = np.loadtxt(f, usecols=(1, 2))
    r, X = a[:, 0], a[:, 1]
    heh = (X/4.0)/(1.0 - X)
    return (float(np.interp(2.0, r, heh)), float(np.interp(10.0, r, heh)),
            h[:, 4].max())


SEP = '-'*70
print('='*70)
print('LHS 1140 b: the XUV grid, the K_zz pressure-form test and the outer-')
print('region ladder, re-solved on the current Penning coefficient')
print('='*70)
print('binary   EXHALE.x built 2026-08-29 04:56  (Update_EXHALE 87, 88)')
print('method   each stored arm re-solved by ./resolve_wind.sh: 4 JFNK passes')
print('         seeded from its own output with the lower atmosphere held')
print('         fixed, then one post-process pass and one transit synthesis')
print('         at R = 68,000.  Composition, XUV scaling and K_zz unchanged.')
print('target   red-pair EW = 1.108 +/- 0.030 %A over 10832.60-10834.20 A')
print('stored   the source directories were read only and are quoted here as')
print('         "was" for comparison; nothing in them was rewritten.')
print()

# ---------------------------------------------------------------- XUV grid
print(SEP)
print('1. THE XUV GRID  (docs sec:xuvgrid, Table tab:xuvgrid)')
print('   scalar base He/H = 2.13 (S), flux-closed He/H = 9.71 (C).')
print('   Sources: ../heh2p13_diff_kzz1e9 and ../xuv*_heh2p13;')
print('            ../flux_closure/heh9p7/k03 and ../xuv*_closure9p7.')
print(SEP)
FRAC = [('1p00', 1.00), ('0p33', 0.33), ('0p30', 0.30), ('0p25', 0.25),
        ('0p20', 0.20), ('0p15', 0.15), ('0p10', 0.10), ('0p01', 0.01)]
SRC = {'S1p00': '../heh2p13_diff_kzz1e9',
       'C1p00': '../flux_closure/heh9p7/k03'}
new, old = {}, {}
pairs = []
for t, f in FRAC:
    for fam in ('S', 'C'):
        lab = fam + t
        if os.path.isfile(os.path.join(lab, 'tpm_He10830.txt')):
            pairs.append((lab, lab))
print('-- re-solved --')
new.update(rows(pairs))
print()
print('-- as stored (retired coefficient, or pre-section-88) --')
opairs = []
for lab, _ in pairs:
    t = lab[1:]
    src = SRC.get(lab, '../xuv%s_%s' % (t, 'heh2p13' if lab[0] == 'S'
                                        else 'closure9p7'))
    opairs.append((lab, src))
old.update(rows(opairs))
print()
for fam, name in (('S', 'scalar base, He/H = 2.13'),
                  ('C', 'flux-closed, He/H = 9.71')):
    # the crossing is read on the SCALED grid only, as the memo reads it:
    # the fiducial point is the model, not a point of the XUV scan
    xs = [f for t, f in FRAC if fam+t in new and f < 1.0]
    dn = [new[fam+t]['red'] for t, f in FRAC if fam+t in new and f < 1.0]
    do = [old[fam+t]['red'] for t, f in FRAC if fam+t in new and f < 1.0]
    print('%s:' % name)
    print('  fiducial (1.00) red [%%]: was %.4f, now %.4f, ratio %.3f'
          % (old[fam+'1p00']['red'], new[fam+'1p00']['red'],
             new[fam+'1p00']['red']/old[fam+'1p00']['red']))
    print('  %-8s %9s %9s %8s' % ('F/fid', 'red was', 'red now', 'ratio'))
    for x, a, b in zip(xs, do, dn):
        print('  %-8.2f %9.4f %9.4f %8.3f' % (x, a, b, b/a))
    for lab, d in (('was', do), ('now', dn)):
        c, how = cross_linear(xs, d, 0.6)
        print('  0.6%% crossing %s: %.3f x fiducial (%s)' % (lab, c, how))
    i33 = xs.index(0.33) if 0.33 in xs else None
    i10 = xs.index(0.10) if 0.10 in xs else None
    if i33 is not None and i10 is not None:
        print('  0.33 -> 0.10 depth falls by a factor: was %.1f, now %.1f'
              % (do[i33]/do[i10], dn[i33]/dn[i10]))
    print()

# ------------------------------------------------------------- K_zz form
print(SEP)
print('2. THE K_zz PRESSURE FORM  (docs sec:kzzform, Table tab:kzzform)')
print('   Sources: ../flux_closure/heh11p1/k01 and ../flux_closure/hi/k06')
print('   (recorded, other seed), ../kzz_power/{heh11p1,hi}_{ref,chem,full}')
print('   /k00.  The profiles themselves were not re-solved, so')
print('   Table tab:kzzformmatch (elemental ratios and carriers at the')
print('   match) is chemistry only and is untouched by this.')
print(SEP)
KARMS = [('11.11', [('recorded, other seed', 'L11p1',
                     '../flux_closure/heh11p1/k01'),
                    ('reference, K_zz = 1e9', 'P11_ref',
                     '../kzz_power/heh11p1_ref/k00'),
                    ('A, chemistry only', 'P11_A',
                     '../kzz_power/heh11p1_chem/k00'),
                    ('B, chemistry and wind', 'P11_B',
                     '../kzz_power/heh11p1_full/k00')]),
         ('2.09', [('recorded, other seed', 'L2p09', '../flux_closure/hi/k06'),
                   ('reference, K_zz = 1e9', 'P2_ref',
                    '../kzz_power/hi_ref/k00'),
                   ('A, chemistry only', 'P2_A', '../kzz_power/hi_chem/k00'),
                   ('B, chemistry and wind', 'P2_B',
                    '../kzz_power/hi_full/k00')])]
SLOPE = 0.3819          # dlnEW/dln(He/H) over the 11.11-12.01 interval
for res, arms in KARMS:
    print('=== reservoir He/H = %s ===' % res)
    print('%-24s %8s %8s %9s %8s %9s %9s %8s %8s'
          % ('run', 'red[%]', 'FWHM[A]', 'EW[%A]', 'EW/obs', 'He/H@2Rp',
             'He/H@10Rp', 'Tmax[K]', 'lgMdot'))
    ewn, ewo = {}, {}
    for lab, tag, src in arms:
        for which, d, store in (('now', tag, ewn), ('was', src, ewo)):
            if not os.path.isdir(d):
                continue
            mt = M.metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
            ew = M.red_ew(os.path.join(d, 'tpm_He10830.txt'))
            h2, h10, tmax = wind_extra(d)
            md = M.wind(d)[0]
            store[lab] = ew
            print('%-24s %8.4f %8.4f %9.4f %8.4f %9.4f %9.4f %8.1f %8.3f'
                  % (lab + ' (' + which + ')', mt.get('red_depth', np.nan),
                     mt.get('fwhm_A', np.nan), ew, ew/M.EW_OBS, h2, h10,
                     tmax, md))
    for which, store in (('was', ewo), ('now', ewn)):
        if len(store) < 4:
            continue
        e0 = store['reference, K_zz = 1e9']
        print('  [%s] same profile, other seed: EW %+.2f %% -- path spread'
              % (which, 100*(store['recorded, other seed']/e0 - 1.0)))
        for lab in ('A, chemistry only', 'B, chemistry and wind'):
            d = math.log(store[lab]/e0)
            print('  [%s] %-22s EW %+.2f %%, crossing 11.73 -> %.2f (%+.2f %%)'
                  % (which, lab, 100*(store[lab]/e0 - 1.0),
                     11.73*math.exp(-d/SLOPE),
                     100*(math.exp(-d/SLOPE) - 1.0)))
    print()

# -------------------------------------------------------- outer inherit
print(SEP)
print('3. THE OUTER REGION A CONVERGED WIND INHERITS')
print('   (docs sec:outerinherit, Table tab:outerinherit)')
print('   Sources: the earlier closed ladder, ../flux_closure/{hi/k06,')
print('   heh3/k03, heh5/k04, heh8/k04, heh9p7/k03, heh10p3/k01,')
print('   heh10p7/k01, heh11p1/k01, heh12_ctl/k02} and the discarded')
print('   jump-seeded rung ../flux_closure/heh12/k04.')
print(SEP)
LAD = [('2.09', 'L2p09', '../flux_closure/hi/k06'),
       ('3.0', 'L3p0', '../flux_closure/heh3/k03'),
       ('5.0', 'L5p0', '../flux_closure/heh5/k04'),
       ('8.0', 'L8p0', '../flux_closure/heh8/k04'),
       ('9.7', 'L9p7', '../flux_closure/heh9p7/k03'),
       ('10.3', 'L10p3', '../flux_closure/heh10p3/k01'),
       ('10.7', 'L10p7', '../flux_closure/heh10p7/k01'),
       ('11.1', 'L11p1', '../flux_closure/heh11p1/k01'),
       ('12.0 one step', 'L12p0ctl', '../flux_closure/heh12_ctl/k02'),
       ('12.0 one jump', 'L12p0jump', '../flux_closure/heh12/k04')]
print('-- re-solved --')
lnew = rows([(c, t) for c, t, s in LAD])
print()
print('-- as stored (retired coefficient) --')
lold = rows([(c, s) for c, t, s in LAD])
print()
print('%-16s %9s %9s %9s %9s %9s %9s'
      % ('case', 'r_drop was', 'now', 'T12 was', 'now', 'out10 was', 'now'))
for c, t, s in LAD:
    print('%-16s %9.2f %9.2f %9.0f %9.0f %9.4f %9.4f'
          % (c, lold[c]['r_drop'], lnew[c]['r_drop'], lold[c]['T12'],
             lnew[c]['T12'], lold[c]['outer10'], lnew[c]['outer10']))
print()
print("""CAVEAT on the "was" column.  It is what the stored directory holds now,
not necessarily what the memo printed.  Two rows differ:
  * 11.1 -- flux_closure/heh11p1/k01 was re-run in place on 2026-08-28 19:04,
    after the coefficient changed, so its stored r_drop/T12/outer10
    (19.47/1900/0.0099) already differ from the memo's printed
    19.47/1891/0.0109, and its stored EW 1.1682 differs from the 1.0853 the
    memo's kzz_power table printed for the same directory.
  * every other row reproduces the memo's printed value to the digit.
The memo's printed numbers for the 11.1 arm cannot be reproduced from the
tree and should be treated as superseded rather than as a "before".""")
print()
print('%-16s %9s %9s %9s %9s' % ('case', 'EW was', 'EW now', 'red was',
                                 'red now'))
for c, t, s in LAD:
    print('%-16s %9.4f %9.4f %9.4f %9.4f'
          % (c, lold[c]['ew'], lnew[c]['ew'], lold[c]['red'], lnew[c]['red']))

# --------------------------------------------------- the exobase argument
print()
print('The ground for discarding the jump-seeded 12.01 rung is that the two')
print('solutions differ only outside the exobase, where the fluid solution')
print('is not valid.  Re-measured with src/utils/collisional_validity.py on')
print('the advection-corrected profiles (it reads a run directory and')
print('changes nothing in it):')
print()
sys.path.insert(0, os.path.join(EX_ROOT, 'src', 'utils'))
import collisional_validity as CV
print('%-16s %14s %14s' % ('case', 'r_exo stored', 'r_exo re-solved'))
vo, vn = [], []
for c, t, s in LAD:
    ro = CV.collisional_diagnosis(s, adv=True)['r_exobase']
    rn = CV.collisional_diagnosis(t, adv=True)['r_exobase']
    print('%-16s %14s %14s'
          % (c, 'above domain' if ro is None else '%.2f' % ro,
             'above domain' if rn is None else '%.2f' % rn))
    if ro is not None:
        vo.append(ro)
    if rn is not None:
        vn.append(rn)
print('ladder range: stored %.1f-%.1f, re-solved %.1f-%.1f R_p'
      % (min(vo), max(vo), min(vn), max(vn)))
print()
print('Where the two 12.01 solutions differ, re-solved: the ratio of the')
print('jump-seeded profile to its control, on the _adv profiles.')


def prof(d):
    o = os.path.join(d, 'output')
    hn = M.cols(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    sn = M.cols(os.path.join(o, 'Ion_species_adv.txt'))
    h = np.loadtxt(os.path.join(o, 'Hydro_ioniz_adv.txt'))
    a = np.loadtxt(os.path.join(o, 'Ion_species_adv.txt'))
    return h[:, 0], h[:, hn.index('T[K]')], a[:, sn.index('HeITR')]


r, Tc, nc = prof('L12p0ctl')
_, Tj, nj = prof('L12p0jump')
print('%8s %10s %10s' % ('r [R_p]', 'T ratio', 'n23S ratio'))
for rr in (1.5, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0):
    i = int(np.argmin(abs(r - rr)))
    print('%8.2f %10.3f %10.3f'
          % (r[i], Tj[i]/Tc[i], nj[i]/max(nc[i], 1e-300)))

# ------------------------------------------------ which figures go stale
print()
print(SEP)
print('4. THE FIGURES THAT READ THESE ARMS')
print(SEP)
print("""
Read out of LHS1140b/make_memo_figures.py.  Nothing here was regenerated and
no file under docs/figures/ was touched; this is a list of what a later
regeneration would move, and by how much, so that the decision to regenerate
can be taken separately.

  lhs1140b_closure_ladder.pdf
      Panel (b) IS the XUV grid of block 1 -- it reads heh2p13_diff_kzz1e9,
      xuv*_heh2p13, flux_closure/heh9p7/k03 and xuv*_closure9p7 directly
      (XUV_FAMILIES, XUV_GRID).  Every plotted point moves, the two curves
      cross the 0.6 per cent guide line at a different place, and the
      annotation "0.3x fiducial" drawn at x = 0.30 no longer marks either
      crossing.  Panel (a) reads the discarded 12.01 point from
      flux_closure/heh12 (CLOSURE_DISCARDED), which also moves; the ladder
      curve itself reads crossings_pc090/rw_*, which is already current.
      -> the most stale of the set.

  lhs1140b_closure.pdf
      Panel (a) draws the He 10830 line of flux_closure/{ref,lo,hi}.  Only
      the hi arm (hi/k06) is re-solved here, and its red-pair EW moves
      0.4272 -> 0.5218 %A.  The ref and lo arms are not re-solved in this
      directory, so a regeneration needs them too.  Panels (b) and (c) are
      the handed-over column and do not move.

  lhs1140b_thermostat.pdf, lhs1140b_composition_profiles.pdf
      Both read flux_closure/hi/k06 and flux_closure/heh10p3/k01 (plus
      heh2p13_diff_kzz1e9).  The metastable density is what the Penning
      coefficient acts on, so the n(2^3S) curves move; the temperature
      curves move little.

  lhs1140b_knudsen.pdf
      Reads flux_closure/hi/k06 among its four cases.  Its exobase moves
      18.84 -> 18.79 R_p, so this figure is essentially unchanged.

  lhs1140b_photochem_column.pdf, lhs1140b_photochem_column_2.pdf
      Read lower_atmosphere_profile.dat only.  Chemistry, not wind: not
      stale.
""")

# ----------------------------------------------- the pass record of each arm
print()
print(SEP)
print('5. THE PASS RECORD OF EVERY ARM')
print('   Four JFNK passes per arm, seeded from the arm\'s own output; the')
print('   line below each tag is what the solver returned per pass.  Arms')
print('   that stop at info = 2 are recorded as they stopped and were not')
print('   forced -- the stored solutions they replace stopped there too.')
print('   C0p20 and C0p01 were interrupted mid-pass by an unrelated process')
print('   kill and restarted from scratch; both reproduced their earlier')
print('   pass-1 residuals to four digits, so the restart was deterministic.')
print(SEP)
for f in sorted(os.listdir('.')):
    if f.endswith('.log') and os.path.isdir(f[:-4]):
        print('%-10s %s' % (f[:-4],
                            ' | '.join(l.strip().split('] ', 1)[-1]
                                       for l in open(f) if l.strip())))
