#!/usr/bin/env python3
"""Record of the 2026-08-28 diagnosis of the structure scale L in
   src/utils/collisional_validity.py.

   Until 2026-08-28 that tool used L = min(H_p, L_v, r) -- the pressure
   scale height alone wherever the flow is too slow for L_v -- and across
   the heating peak of the LHS 1140 b runs H_p passes through a broad
   maximum, because rho falls and T rises with nearly the same log slope and
   the front is close to isobaric.  Kn showed a spurious dip there.  The
   tool now uses the state-gradient RMS over (rho, T, v) plus 1/r; this
   script rebuilds the OLD scale locally, so it stays a faithful record of
   what was diagnosed, and prints it against what the tool gives now.

   Diagnostic only -- nothing here modifies collisional_validity.py.
"""
import sys, numpy as np
EXR='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
sys.path.insert(0, EXR+'/src/utils'); sys.path.insert(0, EXR+'/examples')
import collisional_validity as CV, exhale_io as eio

MACH_FLOOR_OLD = 1.0e-2      # the floor the old L_v carried


def old_structure_scale(r_cm, p, v, cs):
    """L = min(H_p, L_v, r), the definition in force before 2026-08-28."""
    H_p = 1.0/np.maximum(np.abs(np.gradient(np.log(p), r_cm)), 1e-99)
    L_v = np.abs(v)/np.maximum(np.abs(np.gradient(v, r_cm)), 1e-99)
    L_v = np.where(np.abs(v) >= MACH_FLOOR_OLD*cs, L_v, np.inf)
    return np.minimum(np.minimum(H_p, L_v), r_cm), H_p, L_v


def scales(d, adv=True):
    run = eio.load_run(d+'/output', d+'/input.inp', adv=adv)
    sl = slice(CV.N_GHOST, -CV.N_GHOST)
    Rp = run.inp['Rp_RJ']*eio.RJ
    r = run.r[sl]; r_cm = r*Rp
    p = run.p[sl]; v = run.v[sl]; T = run.T[sl]; rho = run.n[sl]*CV.m_H_g
    cs = np.sqrt(CV.gamma_ad*p/rho)
    L_old, H_p, L_v = old_structure_scale(r_cm, p, v, cs)
    H_T = 1.0/np.maximum(np.abs(np.gradient(np.log(T), r_cm)), 1e-99)
    H_r = 1.0/np.maximum(np.abs(np.gradient(np.log(rho), r_cm)), 1e-99)
    res = CV.collisional_diagnosis(d, adv=adv)     # the current definition
    return dict(r=r, r_cm=r_cm, T=T, Kn=res['Kn_bulk'], Kns=res['Kn'],
                L=res['L'], L_old=L_old, H_p=H_p, H_T=H_T, H_r=H_r, L_v=L_v,
                heat=run.heat[sl], res=res)


for name, d, adv in [('LHS1140b heh0p55', EXR+'/LHS1140b/exhale/heh0p55', True),
                     ('HD209458b control', EXR+'/backup/phase_d_baseline/new_kzz1e9_d3b', True)]:
    s = scales(d, adv)
    f = s['L_old']/s['L']            # how much longer the old scale was
    r = s['r']; Kn = s['Kn']; Kno = Kn/f       # Kn under the old scale
    j0 = int(np.argmax(s['heat']))
    print('== %s   (heating peak r=%.3f)'%(name, r[j0]))
    print('   max Kn_bulk overall : old %.4g at r=%.3f   ->  now %.4g at r=%.3f'
          %(Kno.max(), r[Kno.argmax()], Kn.max(), r[Kn.argmax()]))
    m = r >= r[j0]
    print('   max Kn_bulk above the heating peak : old %.4g -> now %.4g'%(Kno[m].max(), Kn[m].max()))
    print('   old L longer than the current one by up to %.1fx, at r=%.4f'%(f.max(), r[f.argmax()]))
    m2 = (r > 1.0) & (r < 1.15)
    if m2.any():
        k = int(np.where(m2)[0][np.argmin(Kno[m2])])
        print('   local Kn_bulk minimum in 1.0-1.15 R_p, old scale : %.4g at r=%.4f'
              '  (H_p %.3e, H_T %.3e, H_rho %.3e cm) -> now %.4g'
              %(Kno[k], r[k], s['H_p'][k], s['H_T'][k], s['H_r'][k], Kn[k]))
    print()

# --- max Kn in the critical region (heating peak -> sonic point), as the
#     verdict table of docs/collisional_validity.md defines it ---
print('=== max Kn in the critical region, old scale vs the current one ===')
for name, d in [('LHS1140b heh0p55', EXR+'/LHS1140b/exhale/heh0p55'),
                ('HD209458b control', EXR+'/backup/phase_d_baseline/new_kzz1e9_d3b')]:
    s = scales(d); res = s['res']; r = s['r']
    f = s['L_old']/s['L']
    j0 = int(np.argmax(s['heat'])); rs = res['r_sonic']
    rtop = rs if rs else r[-1]
    m = (r >= r[j0]) & (r <= rtop)
    print(' %s  critical region %.3f - %.3f R_p'%(name, r[j0], rtop))
    print('   %-8s %12s %12s'%('species','max Kn old','max Kn now'))
    for sp in ('bulk',)+tuple(sorted(s['Kns'])):
        k = s['Kn'] if sp == 'bulk' else s['Kns'][sp]
        if not np.isfinite(k[m]).any(): continue
        print('   %-8s %12.4g %12.4g'%(sp, np.nanmax((k/f)[m]), np.nanmax(k[m])))
    print()
