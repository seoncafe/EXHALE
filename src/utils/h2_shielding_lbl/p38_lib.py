"""Helpers for P38: read the EXHALE CLOUDY table, the DB96 fit, and Meudon PDR output."""
import os, re
import numpy as np

EXHALE = '/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00'
TABF = os.path.join(EXHALE, 'src/modules/lower_atmosphere/h2_self_shielding_table.f90')

KB_LW = 1.380649e-16
M_H2 = 2.0 * 1.67353284e-24


def h2_doppler_parameter(T):
    return np.sqrt(2.0 * KB_LW * np.maximum(T, 1.0) / M_H2)


def h2_self_shielding_draine_bertoldi(N_H2, b):
    x = np.maximum(np.asarray(N_H2, dtype=float), 0.0) / 5.0e14
    b5 = np.maximum(b, 1.0) / 1.0e5
    s = np.sqrt(1.0 + x)
    return 0.965 / (1.0 + x / b5) ** 2 + 0.035 / s * np.exp(-8.5e-4 * s)


def _fort_array(name, text):
    """Pull the numbers of a Fortran parameter array declaration."""
    i = text.index(name)
    j = text.index('/)', i)
    body = text[i:j]
    body = body[body.index('(/') + 2:]
    body = re.sub(r'&\s*\n', ' ', body)
    body = body.replace('\n', ' ')
    vals = [float(v.replace('d', 'e').replace('D', 'e'))
            for v in re.findall(r'[-+]?\d*\.\d+[dDeE][-+]?\d+', body)]
    return np.array(vals)


def load_cloudy_table():
    txt = open(TABF).read()
    logT = _fort_array('h2_shield_log_temp(9)', txt)
    logn = _fort_array('h2_shield_log_dens(3)', txt)
    logN = _fort_array('h2_shield_log_col(39)', txt)
    ov = _fort_array('h2_shield_log_col_overlap(n_temp_sh,n_dens_sh)', txt)
    if 'h2_shield_log_fshield(' in txt:
        f = _fort_array('h2_shield_log_fshield(n_col_sh,n_temp_sh,n_dens_sh)', txt)
        renorm = False
    else:
        # 2026-09-02 restructure: the table carries log10 sigma_diss [cm^2];
        # f_shield is sigma_diss(N)/sigma_diss(N_bottom), as
        # h2_self_shielding_level_resolved forms it.
        f = _fort_array('h2_shield_log_sigma_diss(n_col_sh,n_temp_sh,n_dens_sh)', txt)
        renorm = True
    assert ov.size == 27, ov.size
    assert f.size == 39 * 9 * 3, f.size
    ov = ov.reshape((3, 9)).T          # Fortran (n_temp, n_dens)
    f = f.reshape((3, 9, 39)).transpose(2, 1, 0)   # (n_col, n_temp, n_dens)
    if renorm:
        f = f - f[0, :, :][None, :, :]
    return dict(logT=logT, logn=logn, logN=logN, log_overlap=ov, log_f=f)


def cloudy_fshield(N_H2, T, n_H, tab=None):
    """Trilinear in (log10 T, log10 n_H, log10 N_H2), clamped, as the Fortran does."""
    if tab is None:
        tab = load_cloudy_table()

    def loc(ax, v):
        v = min(max(v, ax[0]), ax[-1])
        i = np.searchsorted(ax, v) - 1
        i = min(max(i, 0), len(ax) - 2)
        w = (v - ax[i]) / (ax[i + 1] - ax[i])
        return i, w

    it, wt = loc(tab['logT'], np.log10(T))
    idn, wd = loc(tab['logn'], np.log10(n_H))
    N_H2 = np.atleast_1d(np.asarray(N_H2, dtype=float))
    out = np.empty_like(N_H2)
    for k, N in enumerate(N_H2):
        ic, wc = loc(tab['logN'], np.log10(max(N, 1e-30)))
        v = 0.0
        for a, wa in ((0, 1 - wc), (1, wc)):
            for b_, wb in ((0, 1 - wt), (1, wt)):
                for c, wcd in ((0, 1 - wd), (1, wd)):
                    v += wa * wb * wcd * tab['log_f'][ic + a, it + b_, idn + c]
        out[k] = 10.0 ** v
    return out


def cloudy_overlap_column(T, n_H, tab=None):
    if tab is None:
        tab = load_cloudy_table()

    def loc(ax, v):
        v = min(max(v, ax[0]), ax[-1])
        i = np.searchsorted(ax, v) - 1
        i = min(max(i, 0), len(ax) - 2)
        return i, (v - ax[i]) / (ax[i + 1] - ax[i])
    it, wt = loc(tab['logT'], np.log10(T))
    idn, wd = loc(tab['logn'], np.log10(n_H))
    lo = tab['log_overlap']
    v = ((1 - wt) * (1 - wd) * lo[it, idn] + wt * (1 - wd) * lo[it + 1, idn]
         + (1 - wt) * wd * lo[it, idn + 1] + wt * wd * lo[it + 1, idn + 1])
    return 10.0 ** v


# ---------------------------------------------------------------- Meudon PDR
def read_pdr_ascii(path):
    """Read one of the PDR code's ASCII_files_dat tables.

    The first lines are '@key=value' metadata, one field per column, all
    concatenated on the same line; the rest are numbers.
    Returns (list of column names, ndarray).
    """
    names = None
    rows = []
    for line in open(path):
        if line.startswith('@'):
            if line.startswith('@Name='):
                names = [s[len('Name='):] for s in line.rstrip('\n').split('@')[1:]]
            continue
        line = line.strip()
        if not line:
            continue
        try:
            rows.append([float(v) for v in line.split('val=')])
        except ValueError:
            continue
    return names, np.array(rows)


def pdr_profile(rundir):
    """f_shield(N_H2) and the supporting profile from one PDR run directory."""
    d = os.path.join(rundir, 'ASCII_files_dat')
    npos, pos = read_pdr_ascii(os.path.join(d, 'Positions.dat'))
    ncd, cd = read_pdr_ascii(os.path.join(d, 'ColumnDensities.dat'))
    npr, pr = read_pdr_ascii(os.path.join(d, 'PhotoReactions.dat'))
    ncs, cs = read_pdr_ascii(os.path.join(d, 'CloudStructure.dat'))
    nde, de = read_pdr_ascii(os.path.join(d, 'Densities.dat'))
    iav = npos.index('AV')
    iN = ncd.index('N(H2) profile')
    ik = [i for i, n in enumerate(npr)
          if n == 'Photo reaction probability ( H2 + photon > H + H )'][0]
    out = dict(AV=pos[:, iav], NH2=cd[:, iN], kdiss=pr[:, ik],
               names_cs=ncs, cs=cs, names_de=nde, de=de,
               names_cd=ncd, cd=cd)
    return out
