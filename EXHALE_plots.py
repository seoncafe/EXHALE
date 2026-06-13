"""Live / static plotter for EXHALE output.

Usage:
    python3 EXHALE_plots.py            # plot the current ./output once
    python3 EXHALE_plots.py --live     # refresh every 4 s
    python3 EXHALE_plots.py --live n   # refresh every n seconds

Updated for the EXHALE "schema 2" output files: the '# columns'
header line is parsed, so the plot adapts to the actual column layout
(H/He + He 2^3S + up to ten trace metals) instead of assuming the
original fixed 7-column format.  Legacy headerless files still load
(first 7 columns).  If metals are present (any nonzero metal column) a
third figure shows the leading metal ion densities.
"""

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- options
animate = False
sec = 4.0
if len(sys.argv) >= 2:
    animate = True
    if len(sys.argv) >= 3:
        sec = float(sys.argv[2])

# ---------------------------------------------------------------- constants
mu = 1.673e-24          # proton mass [g]
gam = 5.0/3.0
RJ = 6.9911e9

# ---------------------------------------------------------------- loaders


def read_columns_header(path):
    """Return the column-name list from a '# columns ...' header line,
    or None for legacy headerless files."""
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s.startswith('#'):
                return None
            if s.lstrip('#').split()[:1] == ['columns']:
                return s.lstrip('#').split()[1:]
    return None


def load_ion(path):
    """Load an Ion_species(.adv) file as {name: column}, name as in the
    header ('HI', 'HII', ..., 'FeIII'); 'r' is the radius column."""
    cols = read_columns_header(path)
    data = np.loadtxt(path, comments='#', unpack=True)
    if cols is None:                       # legacy: r + 5 species (+HeITR)
        cols = ['r[Rp]', 'HI', 'HII', 'HeI', 'HeII', 'HeIII',
                'HeITR'][:len(data)]
    out = {nm.split('[')[0] if nm.startswith('r') else nm: col
           for nm, col in zip(cols, data)}
    out['r'] = out.pop('r', data[0])
    return out


def load_hydro(path):
    r, n, v, p, T = np.loadtxt(path, comments='#', unpack=True,
                               usecols=(0, 1, 2, 3, 4))
    return r, n*mu, v, p, T


# ---------------------------------------------------------------- input.inp
R0 = None
appx_mth = ''
with open('input.inp') as fh:
    for line in fh:
        if line.startswith('Planet radius'):
            R0 = float(line.split(':')[1])*RJ
        if line.startswith('2D approximate method'):
            appx_mth = line.split(':')[1].strip()
if R0 is None:
    raise RuntimeError('Planet radius not found in input.inp')

# ---------------------------------------------------------------- load data
r, rho, v, p, T = load_hydro('./output/Hydro_ioniz.txt')
ion = load_ion('./output/Ion_species.txt')
N = r.size

METALS = ['CI', 'CII', 'OI', 'OII', 'NI', 'NII', 'MgI', 'MgII',
          'CaII', 'NaI', 'FeI', 'FeII']
have_metals = any(ion.get(m, np.zeros(1)).max() > 0 for m in METALS)


def hydro_derived(r, rho, v):
    mom = 4.*np.pi*v*rho*r**2*R0**2
    lgmom = np.where(mom > 0, np.log10(np.maximum(mom, 1e-300)), -20.0)
    return lgmom


def ion_derived(ion):
    nh = ion['HI'] + ion['HII']
    nhe = ion['HeI'] + ion['HeII'] + ion['HeIII']
    f = dict(fhi=ion['HI']/nh, fhii=ion['HII']/nh,
             fhei=ion['HeI']/nhe, fheii=ion['HeII']/nhe,
             fheiii=ion['HeIII']/nhe)
    f['fheiTR'] = ion.get('HeITR', np.zeros_like(nh))/nhe
    return f


lgmom = hydro_derived(r, rho, v)
frac = ion_derived(ion)

vlim = 1.2*v.max()*1e-5 if v.max() > 0 else 15.0
mom_inf = lgmom.min()
if mom_inf == -20.0:
    mom_inf = 0.8*lgmom[min(200, N - 1)]

# correct the outer mass flux for the 3D approximation method
mom_out = lgmom[-20]
if appx_mth.startswith('Rate/2'):
    mom_out += np.log10(0.5)
if appx_mth.startswith('Mdot/4') or appx_mth.startswith('Rate/4'):
    mom_out += np.log10(0.25)

# ---------------------------------------------------------------- main figure
fig, ax = plt.subplots(2, 4)
fig.set_size_inches(11, 7)
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.93,
                    hspace=0.30, wspace=0.31)

rho_line, = ax[0, 0].semilogy(r, rho)
ax[0, 0].set_title('Density [g cm$^{-3}$]', fontdict={'weight': 'bold'})

v_line, = ax[0, 1].semilogy(r, v*1e-5)
ax[0, 1].set_ylim([1.e-3, vlim])
ax[0, 1].set_title('Velocity [km s$^{-1}$]', fontdict={'weight': 'bold'})

p_line, = ax[0, 2].semilogy(r, p)
ax[0, 2].set_title('Pressure [erg cm$^{-3}$]', fontdict={'weight': 'bold'})

T_line, = ax[0, 3].plot(r, T)
ax[0, 3].set_title('Temperature [K]', fontdict={'weight': 'bold'})

mom_line, = ax[1, 0].plot(r, lgmom)
ax[1, 0].set_title('Log10 Momentum [g s$^{-1}$]',
                   fontdict={'weight': 'bold'})
ax[1, 0].set_ylim([mom_inf, 1.2*lgmom.max()])

ION_PANEL = [('HI', '$n_{HI}$', None), ('HII', '$n_{HII}$', None),
             ('HeI', '$n_{HeI}$', None), ('HeII', '$n_{HeII}$', None),
             ('HeIII', '$n_{HeIII}$', None),
             ('HeITR', '$n_{HeI3}$', '#5baca7')]
ion_lines = {}
for nm, lab, c in ION_PANEL:
    if nm in ion:
        ion_lines[nm], = ax[1, 1].semilogy(r, ion[nm], label=lab, color=c)
ax[1, 1].set_title('Ion densities [cm$^{-3}$]', fontdict={'weight': 'bold'})
ax[1, 1].legend(loc='upper right', fontsize=8)

fhi_line, = ax[1, 2].plot(r, frac['fhi'], label='$f_{HI}$')
fhii_line, = ax[1, 2].plot(r, frac['fhii'], label='$f_{HII}$')
ax[1, 2].set_title('H fractions', fontdict={'weight': 'bold'})
ax[1, 2].set_ylim([0, 1])
ax[1, 2].legend(loc='best')

fhei_line, = ax[1, 3].plot(r, frac['fhei'], label='$f_{HeI}$')
fheii_line, = ax[1, 3].plot(r, frac['fheii'], label='$f_{HeII}$')
fheiii_line, = ax[1, 3].plot(r, frac['fheiii'], label='$f_{HeIII}$')
fheiTR_line, = ax[1, 3].plot(r, frac['fheiTR'], label='$f_{HeI3}$',
                             color='#5baca7')
ax[1, 3].set_title('He fractions', fontdict={'weight': 'bold'})
ax[1, 3].set_ylim([0, 1])
ax[1, 3].legend(loc='best')

for a in ax.flat:
    a.set_xlim([r[0], r[-1]])
    a.set_xlabel('r/R$_P$')

# ------------------------------------------------- post-processed overlay
adv_hydro = os.path.join(os.getcwd(), 'output', 'Hydro_ioniz_adv.txt')
adv_ioniz = os.path.join(os.getcwd(), 'output', 'Ion_species_adv.txt')
if os.path.isfile(adv_hydro) and os.path.isfile(adv_ioniz):
    ra, rhoa, va, pa, Ta = load_hydro(adv_hydro)
    iona = load_ion(adv_ioniz)
    fa = ion_derived(iona)
    ax[0, 2].semilogy(ra, pa, '--')
    ax[0, 3].plot(ra, Ta, '--')
    cyc = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for k, (nm, lab, c) in enumerate(ION_PANEL):
        if nm in iona:
            ax[1, 1].semilogy(ra, iona[nm], '--',
                              color=(c or cyc[k % len(cyc)]))
    ax[1, 2].plot(ra, fa['fhi'], '--', color=cyc[0])
    ax[1, 2].plot(ra, fa['fhii'], '--', color=cyc[1])
    ax[1, 3].plot(ra, fa['fhei'], '--', color=cyc[0])
    ax[1, 3].plot(ra, fa['fheii'], '--', color=cyc[1])
    ax[1, 3].plot(ra, fa['fheiii'], '--', color=cyc[2])
    ax[1, 3].plot(ra, fa['fheiTR'], '--', color='#5baca7')

# ------------------------------------------------------ metal-ion figure
if have_metals:
    figm, axm = plt.subplots(figsize=(7.5, 5))
    for nm in METALS:
        if nm in ion and ion[nm].max() > 0:
            axm.semilogy(r, ion[nm], label=nm)
    axm.set_xlim([r[0], r[-1]])
    ymax = max(ion[m].max() for m in METALS if m in ion)
    axm.set_ylim([ymax*1e-12, ymax*3])
    axm.set_xlabel('r/R$_P$')
    axm.set_ylabel('n [cm$^{-3}$]')
    axm.set_title('Metal ion densities', fontdict={'weight': 'bold'})
    axm.legend(ncol=3, fontsize=8)

# ---------------------------------------------------------------- report
print('2D approximation method: ', appx_mth)
print('Log10 of mass-loss-rate = ', mom_out)

if not animate:
    plt.show()
else:
    plt.show(block=False)

# ---------------------------------------------------------------- live loop
while animate:
    try:
        r, rho, v, p, T = load_hydro('./output/Hydro_ioniz.txt')
        ion = load_ion('./output/Ion_species.txt')
    except Exception:
        time.sleep(sec)        # file being rewritten; retry
        continue
    lgmom = hydro_derived(r, rho, v)
    frac = ion_derived(ion)

    rho_line.set_ydata(rho)
    v_line.set_ydata(v*1e-5)
    p_line.set_ydata(p)
    T_line.set_ydata(T)
    mom_line.set_ydata(lgmom)
    for nm, ln in ion_lines.items():
        ln.set_ydata(ion[nm])
    fhi_line.set_ydata(frac['fhi'])
    fhii_line.set_ydata(frac['fhii'])
    fhei_line.set_ydata(frac['fhei'])
    fheii_line.set_ydata(frac['fheii'])
    fheiii_line.set_ydata(frac['fheiii'])
    fheiTR_line.set_ydata(frac['fheiTR'])

    vlim = 1.2*v.max()*1e-5 if v.max() > 0 else 15.0
    mom_inf = lgmom.min()
    if mom_inf == -20.0:
        mom_inf = 0.8*lgmom[min(200, N - 1)]
    ax[0, 1].set_ylim([1.e-3, vlim])
    ax[0, 3].set_ylim(0.8*T.min(), 1.2*T.max())
    ax[1, 0].set_ylim([mom_inf, 1.2*lgmom.max()])

    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(sec)
