# Turbulent (non-thermal) broadening of the He 10830 line: what the
# literature under `references/` does

Survey of all 43 PDFs in `../../references/` (2026-08-24), prompted by the
discovery that the Cherubim et al. (2026) p-winds run sets
`turbulence_broadening=True`. Question: is that a p-winds quirk, or the
field's convention?

Method: `pdftotext` over every PDF, grep for `turbulen|microturbulen|
non-thermal broadening`, then read each hit in context. 15 files matched;
most matches are eddy diffusion (Kzz) or plasma/flow turbulence, not line
broadening. The line-broadening results are below.

## The convention, where it is used

**Yan et al. (2024)** states it most explicitly (Sect. 4.2). The thermal
Doppler width is replaced by an effective width

    dnu_D = (nu0/c) * sqrt(2kT/m + v_turb^2)

and the paper records that **Salz et al. (2018), Lampon et al. (2020) and
Czesla et al. (2022) all assume v_turb = sqrt(5kT/3m)**, the monatomic
sound speed -- i.e. turbulence at Mach 1.

**Salz et al. (2016)**, Eq. (5), does exactly this: "b is the Doppler
parameter, v_th the thermal Doppler velocity, and for the microturbulence
v_micro we use the sound speed", applied when synthesizing the transit
profile. (Separately, their atmosphere solver carries a constant 1 km/s
microturbulence at the lower boundary, which is a different quantity and
"has little impact on the overall atmospheric structure".)

**p-winds** implements the same thing in its own convention. Its Gaussian
is sqrt(kT/m + v_turb^2) -- sigma, not the b-parameter -- with
`turbulence_velocity = sqrt(5/6 kT/m)`. That is sqrt(5kT/3m)/sqrt(2), which
in the b-parameter convention is exactly the Lampon/Salz/Yan sound speed:
the two forms agree, and the apparent factor sqrt(2) is the sigma-vs-b
convention.

Consequence for EXHALE, which carries the b-parameter
v_th = sqrt(2kT/m): the same physics is b -> sqrt(11/6) * v_th, a factor
1.354. That is what `EXHALE_TRANSIT_TURB=1` applies.

## Who uses it, and who argues against

**Uses it as a fitted or assumed term**

- Lampon et al. (2020, 2021): the reference implementation; turbulent
  broadening is a free/assumed parameter of the 1D isothermal Parker fits.
- Salz et al. (2016, 2018), Czesla et al. (2022): v_turb = c_s.
- Cherubim et al. (2026): `turbulence_broadening=True` in the released
  script, so the LHS 1140b retrieval carries it.

**Reports needing it, and says so plainly**

- **Taylor et al. (2026)** is the closest analogue to this project: a
  self-consistent model (energy balance, diffusion, non-LTE H(n=2)) that
  reproduces the HD 189733b line *core* but not its *width*. They state
  that "several prior analyses have matched the HD 189733b width by
  introducing ad hoc (micro)turbulent or Gaussian broadening terms; here,
  we show that this step remains necessary even in our self-consistent
  framework", and quantify it as **~12 km/s of additional nonthermal
  broadening**. They speculate it is plasma turbulence from ions trapped
  on the planet's magnetic field lines, strengthening with stellar
  activity, and call the idea speculative.
- Rumenskikh et al. (2022), as reported by Taylor et al.: a 3D multifluid
  model that instead imposed reduced H line cooling plus prescribed zonal
  jets of 10-20 km/s to reach the observed width.

**Argues against, or does without**

- **Murray-Clay et al. (2009)**: on Garcia Munoz's (2007) proposal to
  broaden the Ly-alpha line by wind turbulence -- "to generate the large
  velocities observed, the energy in such turbulence would need to exceed
  the thermal and bulk kinetic energies in the mean flow by a factor of
  ~100. Such energy requirements seem insurmountable."
- **Koskinen et al. (2013b)**: "non-thermal broadening such as that
  proposed by Ben-Jaffel and Hosseini (2010) does not appear to be
  necessary to explain the current observations", though they allow that
  stellar-wind interaction may produce broadening turbulence.
- **Schulik & Owen (2025)** offer a physical alternative: adiabatic
  cooling produces a helium-triplet *population bump* at high velocity,
  which widens the line without a free parameter. They are explicit that
  the isothermal alternative "can only explain the broadness of the data
  under physically inconsistent assumptions ... and an empirical free
  parameter, turbulent broadening (Lampon et al. 2021)", and note that
  T = 10^4 K already maximizes what thermal broadening alone can give.
- **Yan et al. (2024)** keeps turbulence *out* of its nominal models and
  tests it only as a sensitivity, arguing that ISM turbulence
  anti-correlates with temperature so that Mach ~ 1 is likely an
  overestimate at a few thousand K; they conclude "the turbulence effect
  is likely less significant than the one explored in this paper".
- Yan et al. (2019) computes the Doppler width explicitly "assuming no
  turbulence".

## Where this leaves the LHS 1140b comparison

The literature does not have a settled treatment. Turbulent broadening at
v_turb = c_s is a common *convention* in 1D isothermal Parker-wind fits
(Lampon, Salz, Czesla, and hence Cherubim), while self-consistent and 3D
models either need an even larger ad hoc term (Taylor: ~12 km/s), reach
the width by other physics (Schulik & Owen: adiabatic-cooling population
bump), or argue the term is unjustified on energy grounds (Murray-Clay)
or unnecessary (Koskinen).

For this project that means: (i) enabling it in EXHALE is the right move
for a like-for-like comparison against the paper's p-winds run, since the
paper has it on; (ii) it should stay **opt-in**, because it is a modeling
choice the field disputes, not a settled ingredient; and (iii) the width
shortfall we measure is a known, published problem, not an artifact of
this setup -- Taylor et al. (2026) report the same shortfall from a
comparably self-consistent model.

Non-matches worth recording so the sweep is not repeated: Wogan (2025),
Lavvas (2014, 2019), Robeling (2026), Koskinen (2013b, second half) use
"turbulent" for eddy diffusion Kzz; 2607.18193v1 for flow-flow
interaction in a 3D stellar-wind simulation; Jensen (2012) for
photospheric turbulence driving stellar activity; Taylor (2025) for
gravity-wave breaking. `Vidal_2003_Nature.pdf` matches only on its first
page, which carries the tail of the preceding Nature letter (a Crab
pulsar paper); Vidal-Madjar et al. themselves do not discuss turbulence.
