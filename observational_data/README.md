# Exoplanet Transmission Spectroscopy Observational Datasets

This folder contains public observational data files for Lyman-alpha ($\text{Ly}\alpha$), H-alpha ($\text{H}\alpha$), and Helium I ($10830\text{ \AA}$) lines, formatted to easily compare with model outputs (like `TPM.py`).

## Summary of Datasets

1. **Lyman-alpha ($\text{Ly}\alpha$):**
   - File: [vidalmadjar2003_HD209458b_Lya.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/vidalmadjar2003_HD209458b_Lya.txt)
   - Source: Vidal-Madjar et al. (2003) / Ehrenreich et al. (2008)
   - Format: Two columns (Velocity [km/s], relative flux change $dF/F$)
   - Note: The geocoronal/ISM absorption core ($\pm 40\text{ km/s}$) is omitted.

2. **H-alpha ($\text{H}\alpha$):**
   - File: [cauley2015_HD189733b_Ha.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/cauley2015_HD189733b_Ha.txt)
   - Source: Cauley et al. (2015)
   - Format: Two columns (Wavelength [Angstrom], relative flux change $dF/F$)

3. **Helium I ($10830\text{ \AA}$):**
   - File: [nortmann2018_WASP69b_He.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/nortmann2018_WASP69b_He.txt)
   - Source: Nortmann et al. (2018)
   - Format: Two columns (Wavelength [Angstrom], relative flux change $dF/F$)

4. **Helium I ($10830\text{ \AA}$):**
   - File: [salz2018_HD189733b_He.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/salz2018_HD189733b_He.txt)
   - Source: Salz et al. (2018)
   - Format: Two columns (Wavelength [Angstrom], relative flux change $dF/F$)

---

## How to Compare with Models

These files use the same convention as your existing H-alpha reference file (`HD209458b/jensen2012_HD209458b_Ha.txt`). You can easily load them in Python using:

```python
import numpy as np

# For wavelength-dependent spectra (H-alpha and He I)
lam_obs, dF_obs = np.loadtxt("observational_data/cauley2015_HD189733b_Ha.txt", unpack=True)
T_obs = 1.0 + dF_obs # Transmission probability T_lambda

# For velocity-dependent spectra (Lyman-alpha)
vel_obs, dF_obs = np.loadtxt("observational_data/vidalmadjar2003_HD209458b_Lya.txt", unpack=True)
```
