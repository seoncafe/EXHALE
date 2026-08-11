# Expected Effects and Results of the Huang et al. Updates in EXHALE

This note collects the changes in atmospheric structure, ionization structure,
mass-loss rate and observed spectrum that are expected once the physical
processes proposed in `Huang_update_plan.md` (the extended metal set, Ly-alpha
radiative transfer, the $H(n=2)$ excited state, the Roche potential and the
boundary-condition update) are implemented in EXHALE.

---

## 1. Changes in the temperature structure $T(r)$

### Additional metal-line cooling (Mg II and Fe II)

* **Current state:** only the C, N, O cooling lines contribute ([O I] 63 $\mu$m,
  [C II] 158 $\mu$m). In the HD 209458 b simulation, the peak temperature with
  metal cooling on is about **6,525 K**, about **1,750 K** below the
  **8,272 K** reached with metal cooling off.
* **Expected change:** in an ultra-hot Jupiter atmosphere (WASP-121 b and
  similar) the temperature is higher still, so **Mg II and Fe II line cooling**
  becomes the dominant cooling mechanism in the lower thermosphere
  ($1.15 \lesssim r/R_p \lesssim 1.4$). Adding it cools the lower atmosphere
  further and makes the atmosphere more compact.

### Balmer heating from the $H(n=2)$ excited state

* **Expected change:** as Ly-alpha radiation and electron collisions raise the
  population of the $H(n=2)$ states (2s, 2p), photoionization heating by
  stellar Balmer-continuum photons switches on. This acts as a strong heat
  source in the lower atmosphere ($r < 2 R_p$), where the stellar XUV does not
  penetrate, and locally raises the lower-thermosphere temperature.

### Adiabatic cooling from raising the boundary altitude

* **Expected change:** coupling the $1\,\mu$bar lower boundary to a
  photochemical model of the lower atmosphere moves it outward (to
  **$1.46\,R_p$** for WASP-121 b). The outflow then accelerates, strong
  **adiabatic expansion cooling** sets in, and the temperature inside the Roche
  lobe drops overall.

---

## 2. Changes in the ionization and chemical structure

### Metal ionization structure

* **Expected change:** under a strong stellar radiation field, Mg, Ca, Fe and
  Si are mostly doubly ionized ($Mg^{2+}$, $Ca^{2+}$, $Fe^{2+}$) in the upper
  thermosphere, while Na and K are singly ionized ($Na^+$, $K^+$).

### Ionization coupling through charge exchange

* **Expected change:** charge-exchange reactions such as
  $Fe + H^+ \leftrightarrow Fe^+ + H$ tie the hydrogen and metal ionization
  states tightly together. In the lower atmosphere these rates greatly exceed
  the ordinary radiative recombination rates and act as the **dominant sink
  returning $H^+$ to neutral hydrogen**.

### Higher electron density in the lower thermosphere

* **Expected change:** because Ly-alpha radiative transfer populates $H(n=2)$,
  which is then ionized by the Balmer continuum, the hydrogen ionization
  fraction and the electron density $n_e$ rise sharply in the lower atmosphere
  near the molecular-to-atomic transition.

---

## 3. Changes in the mass-loss rate $\dot{M}$

* **Metal cooling alone (lower mass-loss rate).** With only metal cooling
  enabled the temperature drops, the atmospheric scale height shrinks, and the
  mass-loss rate falls. For HD 209458 b it goes from
  $\sim 3.0 \times 10^{10}$ g/s with metals off to
  $\sim 3.1 \times 10^{9}$ g/s with metals on, a factor of about 10.
* **Roche potential plus the boundary-condition change (much higher mass-loss
  rate).** For an ultra-hot Jupiter such as WASP-121 b, applying both the
  **Roche potential** and the shift of the $1\,\mu$bar boundary altitude
  reverses the picture. With a spherically symmetric potential (Case A) the
  mass-loss rate is only **$0.052\,M_p/\mathrm{Gyr}$**
  ($\sim 3.7 \times 10^{12}$ g/s), whereas the final matched model (Case D),
  which includes Roche-lobe overflow (RLOF) and the raised boundary, gives
  **$1.03\,M_p/\mathrm{Gyr}$** ($\sim 7.3 \times 10^{13}$ g/s) — about a factor
  of 20 higher.

---

## 4. Comparison of transit depths

The table below compares the spherically symmetric model (Case A) with the
final Roche-potential-plus-metals model (Case D) for the WASP-121 b simulation
of Huang et al. (2023). Units: $R_p/R_\star$.

*(Correction, recorded in `Update_EXHALE` §14: these values are the effective
transit radius $R_{\rm eff}/R_\star$, not $R_p/R_\star$ and not an absorption
percent. The "$\sim 30\%$" wording in the prose below is a loose gloss of
0.30.)*

| Quantity / line | Case A (spherical, canonical) | Case D (final matched model, with RLOF) | Observed |
| :--- | :---: | :---: | :---: |
| **Mass-loss rate ($\dot{M}$)** | $0.052 \, M_p/\text{Gyr}$ | **$1.03 \, M_p/\text{Gyr}$** | - |
| **Mg II $\lambda 2796$ (NUV 4 Å)** | 0.182 | **0.302** | $0.309 \pm 0.036$ |
| **Ca II K $\lambda 3935$** | 0.199 | **0.278** | $0.281 \pm 0.009$ |
| **H$\alpha$ $\lambda 6563$** | 0.201 | **0.185** | $0.186 \pm 0.003$ |
| **H$\beta$ $\lambda 4861$** | 0.174 | **0.135** | $0.143 \pm 0.005$ |
| **Na D2 $\lambda 5890$** | 0.152 | **0.147** | $0.147 \pm 0.002$ |
| **NUV fit quality ($\chi^2/N$)** | 2.08 | **1.20** | - |

* **Metal lines (Mg II, Ca II).** As the atmosphere expands along the Roche
  lobe in the tidal direction and the outflow speeds up, the metal lines are
  strongly Doppler-broadened. The deep, broad NUV/optical lines of the
  observations — $\sim 28\%$ for Ca II and $\sim 30\%$ for Mg II — are then
  reproduced satisfactorily.
* **Hydrogen Balmer lines (H$\alpha$, H$\beta$).** A model without metal
  cooling predicts too high an atmospheric temperature and therefore
  overestimates the Balmer absorption depth. In the final model (Case D), metal
  cooling together with the adjusted stellar Ly-alpha flux suppresses the
  H$\alpha$ depth appropriately, bringing it in line with the observed
  $18.5\%$.
