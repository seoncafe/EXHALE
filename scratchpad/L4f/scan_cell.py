import os, shutil, subprocess, sys, numpy as np
EX='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/scratchpad/L4f/'
B=EX+'tree/EXHALE_L4f_diag.x'
src,out=sys.argv[1],sys.argv[2]
cells=[int(x) for x in sys.argv[3].split(',')]
W=EX+out
if os.path.isdir(W): shutil.rmtree(W)
os.makedirs(W+'/output')
for f in ['input.inp','base.inp']:
    if os.path.exists(EX+src+'/'+f): shutil.copy(EX+src+'/'+f, W+'/'+f)
shutil.copy(EX+src+'/output/Ion_species_IC.txt', W+'/output/Ion_species_IC.txt')
base=open(EX+src+'/output/Hydro_ioniz_IC.txt').read().splitlines()
hdr=[l for l in base if l.startswith('#')]; dat=[l for l in base if not l.startswith('#')]
def run():
    e=dict(os.environ); e['OMP_NUM_THREADS']='2'; e['EXHALE_RESIDUAL']='1'; e['EXHALE_L4F']='1'
    subprocess.run([B],cwd=W,stdout=open(W+'/resid.log','w'),stderr=subprocess.STDOUT,env=e)
    return np.loadtxt(W+'/output/row_measures.txt')
open(W+'/output/Hydro_ioniz_IC.txt','w').write('\n'.join(hdr+dat)+'\n')
a0=run()
print('cell   r      m_ene(0)     d meas/d ln p     |out|/s3')
for c in cells:
    j=c+1; eps=1e-6
    v=[float(x) for x in dat[j].split()]; v2=list(v); v2[3]=v[3]*(1+eps); v2[4]=v[4]*(1+eps)
    d2=list(dat); d2[j]=' '+' '.join('%24.16E'%x for x in v2)
    open(W+'/output/Hydro_ioniz_IC.txt','w').write('\n'.join(hdr+d2)+'\n')
    a=run()
    dm=(a[c-1,6]-a0[c-1,6])/eps/a0[c-1,9]
    print('%4d %8.3f %11.3e %15.4f %13.2f'%(c,a0[c-1,0],a0[c-1,12],dm,abs(a0[c-1,13])/a0[c-1,9]))
