import sys, os
sys.path.insert(0,'.')
import numpy as np, p38_lib as L
tab = L.load_cloudy_table()
runs = sys.argv[1:]
store = {}
for r in runs:
    name = os.path.basename(r.rstrip('/'))
    T = {'test1':1300.0,'f0100':100.0}.get(name.split('_')[1], None)
    if T is None: T = float(name.split('_')[1].lstrip('gxf'))
    try:
        d = L.pdr_profile(r)
    except Exception as e:
        print('%-14s not readable yet (%s)' % (name, e)); continue
    N, k = d['NH2'], d['kdiss']
    f = k / k[0]
    m = N > 0
    print()
    print('=== %s   T = %g K   npts=%d   k(0) = %.4e s^-1   N_H2(max) = %.3e' %
          (name, T, len(N), k[0], N[-1]))
    print('%10s %11s %11s %11s %9s' % ('N_H2','Meudon','CLOUDY tab','DB96','M/C'))
    b = L.h2_doppler_parameter(T)
    for Nt in (1e19, 3.16e19, 1e20, 3.16e20, 1e21, 3.16e21, 4.4e21):
        if Nt > N[-1]:
            continue
        fm = 10 ** np.interp(np.log10(Nt), np.log10(N[m]), np.log10(f[m]))
        fc = L.cloudy_fshield([Nt], T, 1e13, tab)[0]
        fd = L.h2_self_shielding_draine_bertoldi(np.array([Nt]), b)[0]
        print('%10.2e %11.3e %11.3e %11.3e %9.2f' % (Nt, fm, fc, fd, fc / fm))
    store['N%d' % int(T)] = N[m]; store['f%d' % int(T)] = f[m]
if store:
    np.savez('meudon_curves.npz', **store)
    print('\nwrote meudon_curves.npz')
