import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
S='/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/scratchpad/L4f/'
sets=[('res_mf',  'metal free, certified, 30 $R_p$',  'C0','-'),
      ('res_L8r45','metal free, certified, 45 $R_p$', 'C2','-'),
      ('res_cno', 'C/N/O, campaign state, 30 $R_p$',  'C3','-'),
      ('res_pass8','C/N/O, entry of outer pass 8',    'C1','-')]
fig,ax=plt.subplots(1,2,figsize=(10.4,4.1))
for d,lab,c,ls in sets:
    try: a=np.loadtxt(S+d+'/output/row_measures.txt')
    except OSError: continue
    r=a[:,0]; m3=a[:,12]; s3=a[:,9]; eo=a[:,13]
    ax[0].loglog(r,np.maximum(m3,1e-16),ls,color=c,lw=1.2,label=lab)
    ax[1].semilogx(r,np.abs(eo)/s3,ls,color=c,lw=1.2,label=lab)
ax[0].axhline(1e-6,color='k',ls=':',lw=1.0)
ax[0].text(1.15,1.4e-6,r'certification tolerance $10^{-6}$',fontsize=8)
ax[0].set_xlabel(r'$r\ [R_p]$'); ax[0].set_ylabel(r'$|R_3|/s_3$, energy row')
ax[0].set_ylim(1e-16,1.0); ax[0].legend(fontsize=7,loc='lower left')
ax[1].set_xlabel(r'$r\ [R_p]$')
ax[1].set_ylabel(r'$|A_+F_3^{\rm out}/\Delta V|\,/\,s_3$')
ax[1].legend(fontsize=7,loc='upper left')
for a_ in ax: a_.grid(alpha=0.25)
fig.tight_layout()
fig.savefig('/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/docs/figures/lhs1140b_L4f_outer_energy_row.pdf')
print('written')
