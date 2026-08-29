#!/usr/bin/env python3
"""One measured row per arm of the three He/H crossings.

Every number is read out of the arm's own run directory: the He 10830 line
out of `tpm_He10830.txt` and `tpm_He10830_metrics.txt`, the wind out of
`pp.log`/`run.log`, and the three quality indicators out of `output/`:

  step      cells where n(2^3S) departs from its local trend by >10 percent,
            the test of ../metastable_step_diagnosis/step_detect.py, with the
            T = 4000 K seam of the retired Penning fit tagged separately;
  Zmetal    cells whose C+N+O total is exactly zero (the element-projection
            failure of Update_EXHALE section 84);
  r_drop    the first radius outside 1.5 R_p at which T falls below half of
            T(12 R_p), and the share of the metastable radial column that lies
            outside 10 R_p -- the two outer-region measures of section 83.

usage: ./measure.py <label>=<dir> ...
"""
import os, sys
import numpy as np

AIR = 10832.057/10829.09114          # vacuum -> air scaling used by the tool
EW_LO, EW_HI = 10832.60, 10834.20    # the red pair, vacuum
EW_OBS, EW_ERR = 1.108, 0.030


def red_ew(path):
    if not os.path.isfile(path):
        return float('nan')
    s = np.loadtxt(path)
    lam = s[:, 0]*AIR
    dep = (s[:, 2].max() - s[:, 2])/s[:, 2].max()*100.0
    m = (lam >= EW_LO) & (lam <= EW_HI)
    return np.trapz(dep[m], lam[m])


def metrics(path):
    out = {}
    if os.path.isfile(path):
        for line in open(path):
            if line.startswith('#') or not line.strip():
                continue
            k, v = line.split()
            out[k] = float(v)
    return out


def cols(path):
    return open(path).readlines()[1].split()[2:]


def profiles(d):
    """r, n(2^3S)_adv, T_adv, C+N+O_adv, and the same at equilibrium."""
    o = os.path.join(d, 'output')
    hn = cols(os.path.join(o, 'Hydro_ioniz.txt'))
    sn = cols(os.path.join(o, 'Ion_species.txt'))
    out = {}
    for w, suf in (('eq', ''), ('adv', '_adv')):
        f = os.path.join(o, 'Ion_species%s.txt' % suf)
        g = os.path.join(o, 'Hydro_ioniz%s.txt' % suf)
        if not (os.path.isfile(f) and os.path.isfile(g)):
            continue
        a, h = np.loadtxt(f), np.loadtxt(g)
        out[w] = dict(r=a[:, 0], n23s=a[:, sn.index('HeITR')],
                      T=h[:, hn.index('T[K]')], sp=a, spn=sn)
    return out


def step_flags(r, n, T, rmin=1.05, tol=0.10):
    """Departure from the LOCAL trend (../metastable_step_diagnosis)."""
    q = np.full_like(n, np.nan)
    q[1:] = n[1:]/np.maximum(n[:-1], 1e-300)
    hits = []
    for i in range(2, len(n)-2):
        if r[i] < rmin or n[i] <= 0 or n[i-1] <= 0:
            continue
        if not (q[i-1] > 0 and q[i+1] > 0):
            continue
        trend = np.sqrt(q[i-1]*q[i+1])
        if not np.isfinite(trend) or trend <= 0:
            continue
        if abs(q[i]/trend - 1.0) > tol:
            seam = (T[i-1]-4000.0)*(T[i]-4000.0) <= 0
            hits.append((r[i], q[i], 'T4000' if seam else
                         ('outer' if r[i] > 29 else 'other')))
    return hits


def metal_dropout(d):
    """Cells whose C+N+O total number density is exactly zero."""
    p = profiles(d)
    if 'adv' not in p:
        return float('nan'), float('nan')
    out = []
    for w in ('eq', 'adv'):
        sp, sn = p[w]['sp'], p[w]['spn']
        idx = [i for i, c in enumerate(sn) if c.startswith(('CI', 'NI', 'OI'))]
        tot = sp[:, idx].sum(axis=1) if idx else np.zeros(sp.shape[0])
        # a run with no metals at all carries the columns and leaves them 0;
        # only a run that has metals somewhere can lose them in a cell
        out.append(int((tot <= 0.0).sum()) if tot.max() > 0.0 else -1)
    return out[0], out[1]


def outer(d):
    """r_drop, and the share of the metastable column beyond 10 R_p."""
    p = profiles(d)
    if 'adv' not in p:
        return float('nan'), float('nan'), float('nan')
    r, T, n = p['adv']['r'], p['adv']['T'], p['adv']['n23s']
    T12 = np.interp(12.0, r, T)
    below = np.where((r > 1.5) & (T < 0.5*T12))[0]
    rd = r[below[0]] if below.size else r[-1]   # never: the domain edge
    Rp = 1.0
    col = np.concatenate([[0.0], np.cumsum(0.5*(n[1:]+n[:-1])*np.diff(r))])
    frac = 1.0 - np.interp(10.0, r, col)/col[-1] if col[-1] > 0 else float('nan')
    return rd, frac, T12


def wind(d):
    md = float('nan')
    p = os.path.join(d, 'pp.log')
    if os.path.isfile(p):
        for line in open(p):
            if 'steady-state Mdot' in line:
                md = float(line.split('=')[-1].split()[0])
    info, resid, npass = '', float('nan'), 0
    logs = sorted(f for f in os.listdir(d)
                  if f.startswith(('run', 'cont')) and f.endswith('.log')
                  and 'march' not in f)
    for f in logs:
        for line in open(os.path.join(d, f)):
            if 'done info=' in line:
                info = line.split('done info=')[1].split()[0]
            if '||R||=' in line and 'done' in line:
                resid = float(line.split('||R||=')[1].split()[0])
        npass += 1
    return md, info, resid, npass


def flux_spread(d):
    """Fractional radial spread of rho*v*r^2 above the escape radius."""
    o = os.path.join(d, 'output')
    f = os.path.join(o, 'Hydro_ioniz.txt')
    if not os.path.isfile(f):
        return float('nan')
    hn = cols(f)
    h = np.loadtxt(f)
    r = h[:, 0]
    q = h[:, hn.index('rho[mH/cm3]')]*h[:, hn.index('v[cm/s]')]*r**2
    m = r > 2.0
    q = q[m]
    return float((q.max()-q.min())/abs(np.mean(q))) if q.size else float('nan')


HDR = ('%-14s %8s %9s %7s %7s %8s %8s %7s  %5s %6s %9s %6s %5s %7s %7s'
       % ('arm', 'He/H', 'log10Mdot', 'red[%]', 'blue[%]', 'FWHM[A]',
          'EW[%A]', 'EW/obs', 'info', '||R||', 'fluxspread',
          'step', 'T4k', 'Zdrop', 'r_drop'))


def row(label, d, heh=None):
    mt = metrics(os.path.join(d, 'tpm_He10830_metrics.txt'))
    ew = red_ew(os.path.join(d, 'tpm_He10830.txt'))
    md, info, resid, npass = wind(d)
    p = profiles(d)
    if 'adv' in p:
        st = step_flags(p['adv']['r'], p['adv']['n23s'], p['adv']['T'])
    else:
        st = []
    t4 = sum(1 for x in st if x[2] == 'T4000')
    zeq, zadv = metal_dropout(d)
    rd, fo, t12 = outer(d)
    if heh is None:
        heh = float('nan')
        ip = os.path.join(d, 'input.inp')
        if os.path.isfile(ip):
            for line in open(ip):
                if line.startswith('He/H number ratio:'):
                    heh = float(line.split(':')[1])
    print('%-14s %8.4f %9.4f %7.4f %7.4f %8.5f %8.5f %7.4f  %5s %6.1e %9.2e '
          '%6d %5d %7s %7.3f  outer10=%.4f T12=%.0f'
          % (label, heh, md, mt.get('red_depth', float('nan')),
             mt.get('blue_depth', float('nan')), mt.get('fwhm_A', float('nan')),
             ew, ew/EW_OBS, info, resid, flux_spread(d), len(st), t4,
             ('off' if zadv < 0 else '%d/%d' % (zeq, zadv)), rd, fo, t12))
    return dict(label=label, heh=heh, ew=ew, mdot=md, info=info, resid=resid,
                nstep=len(st), t4=t4, zadv=zadv, r_drop=rd, outer10=fo,
                T12=t12,
                red=mt.get('red_depth', float('nan')),
                blue=mt.get('blue_depth', float('nan')),
                fwhm=mt.get('fwhm_A', float('nan')))


if __name__ == '__main__':
    print(HDR)
    rows = []
    for a in sys.argv[1:]:
        lab, d = a.split('=', 1)
        h = None
        if ':' in d:
            d, h = d.split(':'); h = float(h)
        rows.append(row(lab, d, h))
