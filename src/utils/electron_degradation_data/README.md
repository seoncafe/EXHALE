# Electron-impact cross sections for `electron_degradation_spencer_fano.py`

Unpacked on 2026-10-07 from the archives kept in `references/` (the archives
stay there; see their `DOWNLOADED_FROM.txt`). Only the files the degradation
tool reads were kept; everything else in the archives (excited initial
states, vibrationally excited H2) was removed after unpacking.

| directory | source archive | data base and citation | files kept |
|---|---|---|---|
| `CCC_H/` | `references/CCC_database/e-H_XSEC_LS.zip`, from https://atom.curtin.edu.au/CCC-DB/H/I/ELECTRON_IMPACT/ | convergent close-coupling (CCC) data base; Bray, Fursa, Kheifets & Stelbovics (2002), J. Phys. B 35, R117 | `README`, `1S.1S` (elastic), `nL.1S` for n = 2-9 and `:L.1S` (n = 10), all L (1s -> nl excitation), `n.1` (summed over l), `TICS.1S` (total ionization), `TCS.1S`; units a0^2 |
| `CCC_He/` | `references/CCC_database/e-He_sig.zip`, from https://atom.curtin.edu.au/CCC-DB/He/I/ELECTRON_IMPACT/ | CCC data base, e-He with up to 314 target states; theory Fursa & Bray (1995), PRA 52, 1279; the archive's README states the results are unpublished as of 13 April 2026 | `README`, every `X.s1S` (from the ground 1^1S) and `X.t2S` (from the metastable 2^3S), with `TICS` and `TCS`; units cm^2; file name `final.initial`, energy = incident energy |
| `MCCC_H2/MCCC-el-H2-X1Sg/` | `references/MCCC_database/MCCC-el-H2.zip`, from https://www.mccc-db.org/ | molecular CCC data base; Scarlett et al. (2021), ADNDT 137, 101361 (electronic excitation, dissociation, ionization); Scarlett et al. (2023), PRA 107, 062804 (rotational, as the file headers state) | `excitation/vi=0/` (X v = 0 to every bound vf and the dissociative excitation `_DE` of each electronic state; the `X1Sg_vf` "excitation + radiative decay" files are not used), `excitation/vibrational_energies/`, `ionisation/vi=0/`, `dissociation/vi=0/`, `rovibrational_excitation/vi=0/Ni=0..31`; units a0^2; thresholds in the file headers |

Data the tool takes from published tables rather than from these files
(typed into the source with their table numbers): Yoon et al. (2008), JPCRD
37, 913, Tables 3, 5, 13, 14 (H2 momentum transfer, v = 0 -> 1, partial
ionization, dissociative attachment); Janev, Reiter & Samm (2003),
Forschungszentrum Juelich report 4105, eq. (74) (H2 v = 0 -> 2); Kim & Rudd
(1994), PRA 50, 3954, eq. (44) and Table I, and Kim, Johnson & Rudd (2000),
PRA 61, 034702, eq. (6) (binary-encounter-dipole secondary spectra);
He I level energies from NIST ASD (`references/nist_asd/`).
