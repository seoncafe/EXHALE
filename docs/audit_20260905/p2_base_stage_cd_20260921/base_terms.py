import numpy as np, sys
fn=sys.argv[1]; lab=sys.argv[2]
cols=None; rows=[]
for line in open(fn):
    if line.startswith('# columns'):
        cols=line.split()[2:]; continue
    if line.startswith('#'): continue
    rows.append([float(x) for x in line.split()])
a=np.array(rows); c={n:i for i,n in enumerate(cols)}
j=a[:,c['j']].astype(int); phys=a[:,c['physical']].astype(int)
fm=a[:,c['face_mass_lo']]; fmh=a[:,c['face_mass_hi']]
alo=a[:,c['area_lo']]; ahi=a[:,c['area_hi']]; vol=a[:,c['volume']]
r=a[:,c['r_cell']]
Phi=fm*alo   # r^2-weighted mass flux through the LOWER face of each row
wind=(r>=1.20)&(phys==1)
mw=np.median(Phi[wind])
print("== %s ==" % lab)
print("wind-window median of r^2*face_mass (lower faces, r>=1.2) = %.6E" % mw)
for jj in range(0,6):
    k=np.where(j==jj)[0][0]
    print("  row j=%d phys=%d r=%.6f : r^2*Phi_lo=%.6E (%.4f x wind)  R_mass=%.6E  R_mom=%.6E  R_ene=%.6E"
          % (jj,phys[k],r[k],Phi[k],Phi[k]/mw,a[k,c['R_mass']],a[k,c['R_momentum']],a[k,c['R_energy']]))
    if jj>=1:
        print("      dF_mass=%.6E dF_energy=%.6E grav_work=%.6E heat=%.6E cool=%.6E"
              % (a[k,c['dF_mass']],a[k,c['dF_energy']],a[k,c['grav_work_over_volume']],a[k,c['heat']],a[k,c['cool']]))
kh=np.where(j==1)[0][0]
print("  r^2*Phi through cell 1 upper face = %.6E (%.4f x wind)" % (ahi[kh]*fmh[kh], ahi[kh]*fmh[kh]/mw))
k2=np.where(j==2)[0][0]
print("  r^2*Phi through cell 2 upper face = %.6E (%.4f x wind)" % (ahi[k2]*fmh[k2], ahi[k2]*fmh[k2]/mw))
