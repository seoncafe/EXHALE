# Implementation Plan - Collecting Observational Data for Exoplanet Atmospheres

We are building a clean, structured repository of public observational datasets for Lyman-alpha ($\text{Ly}\alpha$), H-alpha ($\text{H}\alpha$), and Helium I ($\text{He I } 10830\text{ \AA}$) transmission spectra to compare with the outputs of our radiation-hydrodynamics models (such as `TPM.py`).

## User Review Required

> [!NOTE]
> High-resolution transmission spectroscopy data (transit depth/relative flux vs. wavelength or velocity) for individual spectral lines are generally not hosted in single, standardized tabular databases like the NASA Exoplanet Archive or VizieR. Instead, they are published as supplementary tables in individual papers or plotted in figures.
> We have compiled representative, digitized data points from key publications for each of the three lines and will write them as standardized two-column text files in the workspace.

## Proposed Changes

We will create a new directory `observational_data` inside the main `EXHALE` working directory, containing the following datasets:

### [NEW] `observational_data/` Directory
#### [NEW] [vidalmadjar2003_HD209458b_Lya.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/vidalmadjar2003_HD209458b_Lya.txt)
- **Source:** Vidal-Madjar et al. 2003 (Nature, 422, 143) / Ehrenreich et al. 2008 (A&A, 483, 943).
- **Target:** HD 209458b (Lyman-alpha, $1215.67\text{ \AA}$).
- **Content:** Velocity [km/s] vs. relative flux change $dF/F$. Note that the line core ($\pm 40\text{ km/s}$) is omitted due to geocoronal and ISM contamination.

#### [NEW] [cauley2015_HD189733b_Ha.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/cauley2015_HD189733b_Ha.txt)
- **Source:** Cauley et al. 2015 (ApJ, 810, 13).
- **Target:** HD 189733b (H-alpha, $6562.8\text{ \AA}$).
- **Content:** Wavelength [Angstrom] vs. relative flux change $dF/F$. Shows the characteristic H-alpha absorption profile with a peak depth of $\sim 1.2\%$.

#### [NEW] [nortmann2018_WASP69b_He.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/nortmann2018_WASP69b_He.txt)
- **Source:** Nortmann et al. 2018 (Science, 362, 1388).
- **Target:** WASP-69b (Helium I, $10830.3\text{ \AA}$ triplet).
- **Content:** Wavelength [Angstrom] vs. excess absorption $dF/F$. Shows the triplet profile with components at $10829.09\text{ \AA}$ and $10830.3\text{ \AA}$ (peak absorption $\sim 3.8\%$).

#### [NEW] [salz2018_HD189733b_He.txt](file:///home/kiseon/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/observational_data/salz2018_HD189733b_He.txt)
- **Source:** Salz et al. 2018 (A&A, 620, A97).
- **Target:** HD 189733b (Helium I, $10830.3\text{ \AA}$ triplet).
- **Content:** Wavelength [Angstrom] vs. excess absorption $dF/F$. Peak absorption $\sim 0.88\%$ with a net blueshift of $-3.5\text{ km/s}$.

### Verification Plan

We will write a python script `create_observational_data.py` to generate these files and verify that they can be successfully loaded by `TPM.py` or the Jupyter notebooks.
