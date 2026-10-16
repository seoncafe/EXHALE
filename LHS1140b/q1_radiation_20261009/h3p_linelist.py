#!/usr/bin/env python3
"""H3+ line list (ExoMol MiZATeP, Bowesman et al. 2023, MNRAS 519, 6333)
restricted to upper levels below 24000 cm^-1, and its LTE emission.

Files (downloaded from exomol.com 2026-10-09, not modified):
  1H3_p__MiZATeP.states.bz2  -> states_small.txt (id, E[cm^-1], g_tot, J)
  1H3_p__MiZATeP.trans.bz2   -> trans_Eup24000.txt (upper id, lower id, A[s^-1])
                                 rows whose upper level lies below 24000 cm^-1
  1H3_p__MiZATeP.pf          partition function (T, Q), nuclear spin included
"""
import os
import numpy as np

C2 = 1.4387769          # cm K
HC = 6.62607015e-27 * 2.99792458e10   # erg cm


class H3pLines:
    def __init__(self, d, trans='trans_Eup24000.txt'):
        cache = os.path.join(d, trans.replace('.txt', '_lines.npz'))
        pf = np.loadtxt(os.path.join(d, '1H3_p__MiZATeP.pf'))
        if os.path.exists(cache):
            z = np.load(cache)
            for k in ('Eu', 'El', 'gu', 'gl', 'A', 'nu', 'iu', 'il'):
                setattr(self, k, z[k])
        else:
            import pandas as pd
            st = np.loadtxt(os.path.join(d, 'states_small.txt'))
            tr = pd.read_csv(os.path.join(d, trans), sep=r'\s+', header=None,
                             dtype={0: np.int64, 1: np.int64, 2: np.float64}
                             ).to_numpy()
            E = np.zeros(int(st[:, 0].max()) + 1)
            g = np.zeros_like(E)
            E[st[:, 0].astype(int)] = st[:, 1]
            g[st[:, 0].astype(int)] = st[:, 2]
            iu = tr[:, 0].astype(int)
            il = tr[:, 1].astype(int)
            self.Eu, self.El, self.gu, self.gl = E[iu], E[il], g[iu], g[il]
            self.iu, self.il = iu, il
            self.A = tr[:, 2]
            self.nu = self.Eu - self.El
            keep = (self.nu > 0) & (self.gu > 0)
            for k in ('Eu', 'El', 'gu', 'gl', 'A', 'nu', 'iu', 'il'):
                setattr(self, k, getattr(self, k)[keep])
            np.savez(cache, Eu=self.Eu, El=self.El, gu=self.gu, gl=self.gl,
                     A=self.A, nu=self.nu, iu=self.iu, il=self.il)
        self.pfT, self.pfQ = pf[:, 0], pf[:, 1]

    def Q(self, T):
        return np.interp(T, self.pfT, self.pfQ)

    def line_power(self, T):
        """LTE emission of each line per molecule into 4 pi [erg s^-1]."""
        return (self.gu * self.A * HC * self.nu
                * np.exp(-C2 * self.Eu / T) / self.Q(T))

    def tau_center_per_column(self, T, m_amu=3.02293):
        """Line-centre absorption cross section per H3+ molecule (all
        levels, LTE at T, stimulated emission included) [cm^2]."""
        c = 2.99792458e10
        k = 1.380649e-16
        nl = self.gl * np.exp(-C2 * self.El / T) / self.Q(T)
        sig_int = (self.gu / self.gl * self.A / (8 * np.pi * c * self.nu**2)
                   * nl * (1 - np.exp(-C2 * self.nu / T)))      # cm^2 cm^-1
        b = np.sqrt(2 * k * T / (m_amu * 1.66053907e-24))
        dnu_D = self.nu * b / c                                   # cm^-1
        return sig_int / (np.sqrt(np.pi) * dnu_D)


def miller2013_E(T):
    """Miller et al. 2013 Table 5 fit (W molecule^-1 sr^-1), segments as
    published (no join ramp), for the comparison only."""
    cA = [-81.9599, 0.886768, -0.0264611, 0.000462693, -4.70108e-6,
          2.84979e-8, -1.03090e-10, 2.13794e-13, -2.26029e-16, 8.66357e-20]
    cB = [-92.2048, 0.298920, -0.000962580, 1.82712e-6, -2.04420e-9,
          1.24970e-12, -3.22212e-16]
    cC = [-62.7016, 0.0526104, -7.22431e-5, 5.93118e-8, -2.83755e-11,
          7.35415e-15, -8.01994e-19]
    cD = [-55.7672, 0.0162530, -7.68583e-6, 1.98412e-9, -2.68044e-13,
          1.47026e-17]
    c = cA if T < 300 else cB if T < 800 else cC if T < 1800 else cD
    return np.exp(sum(ci * T**i for i, ci in enumerate(c)))


if __name__ == '__main__':
    import sys
    L = H3pLines(sys.argv[1])
    print('lines kept: %d' % L.nu.size)
    for T in (500., 1000., 1500., 2000., 2500., 3000., 3500.):
        P = L.line_power(T)
        E_lines = P.sum() / (4 * np.pi) * 1e-7      # W/molecule/sr
        o = np.argsort(P)[::-1]
        c = np.cumsum(P[o]) / P.sum()
        nb = lambda a, b: P[(L.nu >= a) & (L.nu < b)].sum() / P.sum()
        print('T=%5.0f  E(MiZATeP,Eu<24000)=%.4e  Miller13=%.4e  ratio %.4f  '
              'lines for 90/99%%: %d/%d  share <1500:%.3f 1500-3500:%.3f '
              '3500-6000:%.3f >6000:%.3f'
              % (T, E_lines, miller2013_E(T), E_lines / miller2013_E(T),
                 np.searchsorted(c, 0.9) + 1, np.searchsorted(c, 0.99) + 1,
                 nb(0, 1500), nb(1500, 3500), nb(3500, 6000), nb(6000, 1e9)))
