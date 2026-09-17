import os, shutil, subprocess, sys, numpy as np
EX='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/scratchpad/L4f/'
B=EX+'tree/EXHALE_L4f_diag.x'
src=sys.argv[1]; out=sys.argv[2]; jrow=501   # cell 500
W=EX+out
if os.path.isdir(W): shutil.rmtree(W)
os.makedirs(W+'/output')
for f in ['input.inp','base.inp']:
    if os.path.exists(EX+src+'/'+f): shutil.copy(EX+src+'/'+f, W+'/'+f)
shutil.copy(EX+src+'/output/Ion_species_IC.txt', W+'/output/Ion_species_IC.txt')
base=open(EX+src+'/output/Hydro_ioniz_IC.txt').read().splitlines()
hdr=[l for l in base if l.startswith('#')]
dat=[l for l in base if not l.startswith('#')]
res=[]
for eps in [0.0,1e-6,2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,-2e-6,-4e-6,-6e-6,-8e-6]:
    v=[float(x) for x in dat[jrow].split()]
    v2=list(v); v2[3]=v[3]*(1.0+eps); v2[4]=v[4]*(1.0+eps)
    d2=list(dat); d2[jrow]=' '+' '.join('%24.16E'%x for x in v2)
    open(W+'/output/Hydro_ioniz_IC.txt','w').write('\n'.join(hdr+d2)+'\n')
    e=dict(os.environ); e['OMP_NUM_THREADS']='2'; e['EXHALE_RESIDUAL']='1'; e['EXHALE_L4F']='1'
    subprocess.run([B],cwd=W,stdout=open(W+'/resid.log','w'),stderr=subprocess.STDOUT,env=e)
    a=np.loadtxt(W+'/output/row_measures.txt')
    res.append((eps,a[-1,6],a[-1,9],a[-1,12],a[-1,16],a[-1,17]))
    print('eps=%+9.2e  R3=%16.9e  s3=%12.5e  meas=%12.5e heat=%12.5e cool=%12.5e'%res[-1], flush=True)
r=np.array(sorted(res))
# straight line through the whole scan
c=np.polyfit(r[:,0],r[:,1],1)
print('slope dR3/deps = %.6e ; residual of the linear fit:'%c[0])
for e,R in zip(r[:,0],r[:,1]):
    print('   eps=%+9.2e  R3-fit = %12.5e  (%.2e of R3)'%(e,R-np.polyval(c,e),abs(R-np.polyval(c,e))/abs(R)))
