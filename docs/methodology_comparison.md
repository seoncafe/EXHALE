# Methodology comparison: independent models of XUV-driven atmospheric escape, and what EXHALE took from each

*This is the markdown edition of `methodology_comparison.tex`; the two are
two formats of one document and carry the same content.*

## 1. Purpose, and the verdict first

This memo compares the methodology of four independent modelling efforts for XUV-driven
atmospheric escape from close-in exoplanets: the general-purpose 1-D radiation-hydrodynamics
code AIOLOS (Schulik & Booth 2023, MNRAS 523, 286); the dedicated thermosphere-ionosphere
escape model of Taylor et al. (2025, 2026), built on the Koskinen et al. code; the
multi-fluid PLUTO-based model of Xing et al. (2023); and the photochemical-hydrodynamic
model of Frelikh & Murray-Clay (2026), and then records, item by item, what EXHALE has
taken from each and what it has not.

The memo was first written as a forward-looking recommendation list. It is now a record of
outcomes, so every candidate carries one of three verdicts: **adopted** (and in which
module), **adopted then replaced** (and why the first choice did not survive), or **not
adopted** (and whether the reason still holds). Of the fifteen items assessed in section 9,
six are adopted, two were adopted and then replaced, and seven are not adopted. The two
largest gaps the original assessment identified (no radius-dependent He/H separation, and
an ad hoc lower boundary) are both closed, and closed further than the recommendation
asked for.

Frelikh & Murray-Clay is the newest of the four and the closest peer to EXHALE's molecular
layer, so it is the source of four of the seven not-adopted items and of one finding that
belongs to neither code alone: the thermal-dissociation rate pair that both use is not
thermodynamically consistent below about 3000 K (section 7). Three further points where the
two codes diverge, and where the judgement went against EXHALE, are the photoelectron
partition in a molecular gas (section 6.6), the absence of a composition boundary condition
at the base (section 6.9), and the sharpness of the modeled H2→H transition (section 6.10).

Every statement about EXHALE below was checked against `src/` rather than against EXHALE's
own documentation; the file and routine names given are the ones in the tree. The AIOLOS,
Taylor, Xing and Frelikh entries are what those papers say. Numbers measured here from
EXHALE output, rather than quoted from a paper, are marked as such at the point of use.

## 2. The codes considered

- **AIOLOS**: the general-purpose 1-D multi-species radiation-hydrodynamics code of
  Schulik & Booth (2023, MNRAS 523, 286; "SB23"), distributed at
  `github.com/Schulik/aiolos`.
- **Taylor et al. (2025)**, *ApJ* 989, 68: "A Multispecies Atmospheric Escape Model with
  Excited Hydrogen and Helium: Application to HD209458b," built on the Koskinen et al.
  (2013a,b; 2022) thermosphere-ionosphere code.
- **Taylor et al. (2026)**, *ApJ* 999, 214, "Helium Escape in Context: Comparative
  Signatures of Four Close-in Exoplanets," a direct application and extension of Taylor
  (2025) to HD 209458 b, HD 189733 b, HD 149026 b and GJ 1214 b.
- **Xing et al. (2023)**, *ApJ* 953, 166: "The Mass Fractionation of Helium in the Escaping
  Atmosphere of HD 209458b," a *multi-fluid* model built inside PLUTO.
- **Frelikh & Murray-Clay (2026)**, *ApJ* 996, 96: "Efficiency of Hydrodynamic Atmospheric
  Escape in Hot Jupiters and Super-Earths," a 1-D photochemical-hydrodynamic model of H-He
  escape with H3+ and Lyα cooling, heat conduction, tidal gravity and photoelectron
  secondary ionization.

"Taylor" denotes the shared Koskinen-based framework used in both Taylor papers unless a
specific paper is named. "Frelikh" denotes Frelikh & Murray-Clay (2026).

## 3. One-paragraph characterization of each code

**AIOLOS** is a *general-purpose* radiation-hydrodynamics code, not a dedicated escape code.
It solves the time-dependent Euler equations (mass, momentum, energy) for an arbitrary
number of species (each a **separate fluid** with its own density, momentum and internal
energy) coupled by an inter-species **friction/drag** solver. It uses finite-volume
**HLLC** fluxes with **PLM** reconstruction (monotonized-central limiter) and second-order
**strong-stability-preserving** Runge-Kutta time stepping, in Cartesian, cylindrical or
**spherical** geometry. Radiation is treated with **flux-limited diffusion (FLD)** over
arbitrary in/out spectral bands, solved implicitly through a block-tridiagonal system.
Escape is one of many problems it can run.

Two things about AIOLOS have to be attributed to the distributed code rather than to SB23,
because the paper does not carry them. The paper implements a C2Ray-style two-species
photoionization scheme and calls a fuller chemical solver future work; the code has since
acquired a reaction-network class (`c_reaction`) with selectable photochemistry levels, the
Black (1981) heating and cooling rates with free-free and H3+ cooling (`photochem.cpp`,
`chemistry.cpp`), grid spacing that can be linear, logarithmic or bi-logarithmic
(`PARI_GRID_TYPE`), and a tidal-field option (`USE_TIDES`). Statements below about AIOLOS's
chemistry, cooling, grid or tides refer to the code.

**Taylor (2025/2026)** is a *dedicated* 1-D spherical thermosphere-ionosphere and
hydrodynamic-escape model (Koskinen et al. 2013/2022), integrated to steady state, with a
**single bulk momentum equation** but **multi-species molecular and eddy diffusion** that
redistributes species, so He and H can separate diffusively. Its distinguishing strengths
are in the completeness of the *microphysics*: a full He(2³S) and H(n=2) non-LTE network
solved by the KPP kinetic preprocessor, a **Lyα Monte Carlo** radiative-transfer model
(Huang et al. 2017) iterated with the hydro solution for Hα, a high-resolution He(2³S)
photoionization cross section (B-spline K-matrix), temperature-dependent Penning ionization,
and self-consistent coupling to a photochemical lower and middle atmosphere at the μbar
level.

**Xing (2023)** is a *dedicated* 1-D spherical escape model whose distinguishing feature is
that it is genuinely **multi-fluid**: H, H+, He, He+ and e- each carry their **own
velocity**, coupled by resonant and non-resonant collisional drag. It is the only one of the
three that captures **He/H mass fractionation** (heavy He under-dragged by escaping H) in
the dynamics itself rather than through a diffusion coefficient. It is built on PLUTO (HLL
Riemann solver, third-order TVD RK) with a Riemann split that carries the
electron/electric-field pressure term. The He(2³S) 10830 Å line is post-processed, not
solved in the hydro loop.

**Frelikh (2026)** is a *dedicated* 1-D spherical escape model built around a molecular
base, and of the four the closest peer to EXHALE's lower atmosphere. It is **single-fluid**
in velocity (all species travel with the bulk u, and the diffusion the code contains is
turned off for the published H-He runs), but it carries a mass-conservation equation for
each of the eight H and He species with a chemical source term. The hydrodynamics are
unusual for this field: no Riemann problem at all, but the **constrained interpolation
profile** (CIP) advection scheme of Yabe & Aoki (1991) on a **staggered mesh**,
operator-split against explicitly differenced source terms, with a Landshoff-plus-von
Neumann artificial viscosity for stability. The 22-reaction H-He network (their Table 1,
compiled from Yelle 2004 and García-Muñoz 2007) is **integrated in time** as a stiff ODE
system by the semi-implicit Euler stepper `StepperSie` of Press et al. (2007) with an
**analytic Jacobian**. Radiative cooling is Lyα (Black 1981) and H3+ (Miller et al. 2013);
heat conduction follows García-Muñoz (2007). Photoelectron secondary ionization uses
Dalgarno et al. (1999), including the H2 vibrational, triplet-dissociation and
Lyman-Werner channels. The published result is a three-layer atmosphere: an H3+-cooled H2
layer at the base, a Lyα-cooled neutral-H layer above it, and a PdV-cooled ionized outflow.
There is no He metastable, no He↔H charge exchange and no metal in the network.

**EXHALE** is a 1-D spherical **single-fluid** escape code of ATES heritage. The
hydrodynamics are finite-volume, second-order, with PLM or WENO3 reconstruction (usable as
a two-stage PLM→WENO3 sequence) and a selectable LLF/**HLLC**/Roe numerical flux
(`flux/Num_Fluxes.f90`); reconstructed face states are dropped to the piecewise-constant
limit wherever they lose positivity (`positivity_limited_faces`,
`states/Reconstruction.f90`), and a gated fourth-difference dissipation restores the damping
of the contact mode that HLLC loses as the flow stalls (`flux/low_mach_dissipation.f90`, key
`Low-Mach damping`, off by default). The steady wind is reached either by explicit marching on the
mass-flux (ρvr²) spread or by a direct **Jacobian-free Newton-Krylov** solve:
right-preconditioned GMRES(m) with pseudo-transient continuation
(`time_step/steady_newton.f90`, key `Solver: Newton`). Ionization of H, He, ten metal
elements in 27 ion stages, the He(2³S) metastable and, when molecular chemistry is on,
H2/H2+/H3+/HeH+ is solved as *one* coupled nonlinear algebraic system (analytic-Jacobian
Newton, MINPACK `hybrd1` fallback; the `System_HeH_*` family). The photoelectron heating
fraction is not a free parameter: it follows the Shull & van Steenberg (1985) branching
evaluated on each cell's own ionized fraction and neutral composition. Helium is transported
as the second component of a binary mixture with stage-resolved friction
(`functions/binary_element_diffusion.f90`), the Navier-Stokes viscous force and heat
conduction are available (`time_step/viscous_conduction.f90`), and the base can be handed
over from a photochemical lower atmosphere as scalars or as a full profile. Diagnostics are
a non-LTE He(2³S) network, a non-LTE H(n=2) population (Christie et al. 2013), Lyα transfer
through a Neufeld escape-probability closure or an imported field, and the transmission
post-processor `EXHALE_transit.py`. Core size: 75 source files, 16,568 SLOC excluding blank
and comment lines (`src/utils/codesize.py`, measured 2026-08-30), of which 6,404 are in
modules that do not exist in upstream ATES.

## 4. Side-by-side methodology table

The four external columns are as described in their papers (AIOLOS code entries marked
"code"); the EXHALE entries were verified in `EXHALE_v1.00/src/`. EXHALE options marked
*opt-in* default to off. A sixth column does not fit across one A4 text block in the PDF
edition, so the comparison is split by axis: table 4a carries the numerical method and
table 4b the physics and the diagnostics.

**Table 4a: numerical method.**

| Axis | AIOLOS (SB23) | Taylor 2025/2026 | Xing 2023 | Frelikh 2026 | **EXHALE** |
|---|---|---|---|---|---|
| Purpose | General multi-species RHD | Dedicated escape + line diagnostics | Dedicated escape, He/H fractionation | Dedicated escape; radiative-cooling efficiency | Dedicated escape + line diagnostics |
| Dimensions | 1-D (sph./cyl./cart.) | 1-D spherical | 1-D spherical | 1-D spherical, substellar ray | 1-D spherical |
| Route to steady state | Time-dependent → steady | Time-dependent → steady | Time-dependent → steady | Time-dependent → steady | RK marching on the ρvr² spread, optionally finished by a direct JFNK solve |
| Fluid model | **Multi-fluid** + friction | **Single bulk** + multi-species diffusion | **Multi-fluid** | **Single fluid**; one mass equation per species, all at the bulk u | **Single fluid** + binary element diffusion |
| Advection scheme | PLM + MC limiter | Flux-conservative van Leer advection (Koskinen et al. 2013a); not restated in the Taylor papers | PLM (PLUTO) | **CIP** cubic interpolation carrying X and ∂_r X, on a staggered mesh | PLM, WENO3, or two-stage PLM→WENO3, with a face positivity guard |
| Riemann solver | HLLC | - | HLL (+ e⁻-pressure split) | **None**; Landshoff + von Neumann artificial viscosity | LLF, HLLC or Roe (`Numerical flux`) |
| Time integrator | Second-order SSP RK | Operator split, semi-implicit Crank-Nicholson, forward to steady | Third-order TVD RK | Operator split: CIP advection, then explicitly differenced source terms | RK marching; GMRES(m) JFNK with PTC for the steady solve |
| Grid | Linear, log or bi-log (code) | Radial | 1024 points, stretched, 1-10 Rp | Nonuniform, staggered, 1-8.8 Rp, ghost cells | Uniform, stretched or mixed; cell count set at runtime |

**Table 4b: physics and diagnostics.**

| Axis | AIOLOS (SB23) | Taylor 2025/2026 | Xing 2023 | Frelikh 2026 | **EXHALE** |
|---|---|---|---|---|---|
| Radiation | **FLD**, multi-band, implicit | XUV attenuation + Lyα Monte Carlo | Radial XUV attenuation, 53 bins | Radial XUV attenuation, 37 bins (Richards et al. 1994), 12-248 eV | Radial XUV attenuation on a multi-bin energy grid (power law or tabulated SED) + Lyα transfer |
| Photoelectron heating | From the C2Ray energy split | Tuned efficiency (20-40%) | Tuned frequency-averaged η (0.1-0.5) | **Dalgarno et al. (1999)** W₀(1+Cxᵅ), with the H2 vibrational, triplet and Lyman-Werner returns added back | Shull & van Steenberg (1985) branching renormalized on the cell's own neutral composition; **no efficiency parameter** |
| Chemistry solve | C2Ray two-species (paper); reaction network (code) | KPP preprocessor over a stiff ODE network | Semi-implicit source terms | **Stiff ODE integrated in time**: semi-implicit Euler (`StepperSie`), **analytic** Jacobian, LU | **One coupled nonlinear algebraic system** per cell; analytic-Jacobian Newton, or MINPACK `hybrd` under a continuation |
| He/H separation | Via friction/drag | Molecular + eddy diffusion of the trace ratio n_He/n_H | Via multi-fluid dynamics | Implemented but turned off for the published H-He runs | Binary diffusion of the mass fraction X = ρ_He/ρ, friction resolved by ionization stage; *opt-in* |
| Metals | Through the code's chemistry | Optional (solar) | None | **None** | **10 elements, 27 ion stages** in the coupled system, with CHIANTI closed-form line cooling |
| Molecules | Through the code's chemistry | H2, H2+, H3+, HeH+ (2026, for HD 189733 b and GJ 1214 b only) | None | Same four species; 22 reactions (their Table 1) | Same four species inside the coupled system, with Lyman-Werner photodissociation; *opt-in* |
| He↔H charge exchange | Through the code's chemistry | Yes | Yes | **Absent from the network** | Huang et al. (2023) Table 4, rows B1/B2: a non-radiative and a radiative channel, separately |
| Cooling | Recombination + Lyα (paper); Black (1981) line/free-free + H3+ (code) | Recombination + H I lines + H3+ | **Lyα only** | Lyα + H3+ only; recombination cooling checked and left off | Lyα, recombination, free-free, CHIANTI metal lines with fine-structure line trapping, H3+ infrared |
| H3+ cooling form | Black (1981) (code) | Miller et al. (2013) | - | Miller et al. (2013), LTE, emission into vacuum; zero above 5000 K | Miller et al. (2013) with the non-LTE departure factor; net exchange with a diluted base field, *opt-in* |
| He(2³S) 10830 | Not a focus | In-loop non-LTE | Post-processed (Yan et al. 2022) | **Not in the network** | In-loop non-LTE (default on) + `EXHALE_transit.py` |
| H(n=2)/Hα | Not a focus | Non-LTE + Lyα Monte Carlo | Post-processed | Not computed | Non-LTE (Christie et al. 2013) + Lyα escape probability or an imported field |
| Viscosity / conduction | - | Molecular + eddy diffusion, conduction | Collisional drag | Heat conduction (García-Muñoz 2007); no viscosity | Navier-Stokes viscosity and heat conduction (Crank-Nicolson); *opt-in* |
| Tidal/Roche | `USE_TIDES` (code) | Tested (2025); off (2026) | Stellar tidal term | Tidal term (their Appendix C) | Roche potential, spherical or Roche-lobe domain |
| Lower boundary | Configurable | **μbar, coupled to photochemistry** | 1 Rp, fixed n, T = 1500 K | ≈1 μbar: **species densities and p imposed**, u extrapolated; Newtonian relaxation to T₀ below τ=3 | μbar base: fixed T,n; or an analytic column; or a scalar `base.inp`; or a full photochemical profile with K_zz(p) |
| Transmission | - | Ray-traced Voigt | Post-processed 10830 | H3+ emission spectrum (ExoCross); no transit spectrum | `EXHALE_transit.py`: 10830, Lyα, Hα, Hβ, metal doublets, optional triaxial Roche geometry |

## 5. What they have in common

1. **1-D spherical XUV-driven escape**, with photoionization heating balanced by
   Lyα/recombination cooling and adiabatic expansion, is the shared physical backbone.
2. **A steady transonic wind is the target.** The four external codes evolve the
   time-dependent equations forward; EXHALE marches in the same way but can also solve the
   steady equations directly. The physical end state is the same.
3. **Much of the same core H/He atomic physics.** Voronov collisional ionization, radiative
   and dielectronic recombination, and H↔He charge exchange appear in Xing, Taylor and
   EXHALE with near-identical rate coefficients, differing mainly in the source table
   (Storey & Hummer, Badnell, Verner). EXHALE's charge-exchange network is now Huang et al.
   (2023) Table 4 (`radiation/charge_exchange.f90`). Frelikh is the exception on the last of
   these: its Table 1 has no He⁺ + H and no He + H⁺ reaction of any kind, and neither does
   its text (section 6.4).
4. **He(2³S) as the key observable.** Taylor, Xing and EXHALE all target the 10830 Å line
   through a non-LTE 2³S population (recombination, collisions, radiative decay, Penning
   ionization, photoionization). AIOLOS and Frelikh do not; the He metastable is not a
   species in either.
5. **Finite-volume, slope-limited, shock-capturing hydrodynamics** underlies AIOLOS, Xing
   (PLUTO) and EXHALE. Two of the five are outliers, and in the same direction: the Koskinen
   solver beneath Taylor is flux-conservative van Leer advection operator-split against a
   semi-implicit Crank-Nicholson step (Koskinen et al. 2013a), and Frelikh is CIP advection
   operator-split against explicitly differenced source terms. Both are staggered, both use
   an added dissipation in place of an upwind flux, and **neither poses a Riemann problem**.
   Neither Taylor paper restates its numerical scheme; Frelikh gives its scheme in full in
   its Appendix A.

One item that used to be on this list no longer belongs. It read: none of the four derives
the photoelectron heating efficiency from first principles, and all treat it as the primary
free knob. That is no longer true of EXHALE. Its heating rate is the integral of the excess
photon energy over the energy grid weighted by the Shull & van Steenberg (1985) heat
fraction, and that fraction is a function of the local ionized fraction with the ionization
branching renormalized onto the cell's own neutral H and He (`svs85_secondary_branching`,
`radiation/util_ion_eq.f90`); there is no heating-efficiency key in `input.inp` at all.
Taylor and Xing still tune theirs, and AIOLOS takes its split from C2Ray.

## 6. Where they genuinely differ

### 6.1 Fluid model

- **Xing is multi-fluid**, each species carrying its own velocity, which is what reproduces
  **He/H mass fractionation**: because He is four times heavier it is under-dragged by the H
  wind, so the *escaping* He/H (~0.039) is about half the bulk abundance (0.086), a weak
  10830 line without a subsolar He abundance.
- **Taylor keeps a single bulk momentum equation and adds molecular and eddy diffusion**,
  which produces the same qualitative diffusive separation at far lower cost. Two numbers
  are quoted in Taylor (2025) and they belong to different models: the best-fit model's
  elemental He/H falls from 8% at the base of the thermosphere to about 2.5% at high
  altitude, while the model Taylor sets against Xing's and Schulik & Owen's (2025)
  multi-fluid results is Model B, whose H/He goes from 92/8 at the base to about 96/4.
- **AIOLOS is multi-fluid but escape-agnostic**: its friction solver would in principle give
  fractionation, but the code is not tuned for the transonic-wind problem.
- **EXHALE keeps one bulk velocity and diffuses the composition on top of it**, which is
  Taylor's choice of architecture but not Taylor's variable. Where Taylor transports the
  trace ratio n_He/n_H, the binary operator (`binary_element_diffusion.f90`) transports the
  helium *mass fraction* X = ρ_He/ρ of a two-component mixture, which is bounded at both
  ends of the composition axis and therefore does not break down as helium approaches unity.
  Hydrogen's carriers are H I, H II, H2, H2+ and H3+, helium's are He I, He II and He III,
  and the binary coefficient is a Blanc's-law mixture over their stage fractions, with
  hard-sphere (Banks & Kockarts 1973), polarization ion-neutral (Schunk & Nagy) and Coulomb
  ion-ion (Paquette et al. 1986) pair coefficients chosen by charge. That last point is the
  physical content: an ion is held to the protons by Coulomb friction instead of settling at
  a neutral rate.

### 6.2 Radiation transport

AIOLOS is the only one with a full **FLD** thermal-radiation solver, which is necessary for
optically-thick interiors and dust and overkill for the escaping thermosphere. Taylor is the
only one doing **Lyα Monte Carlo** transfer iterated with the hydro solution. Xing and
EXHALE both use radial XUV attenuation: EXHALE on an energy grid built per run from the
ionization thresholds, driven by either a power-law XUV spectrum or a tabulated SED. On top
of that EXHALE closes Lyα trapping either with the escape probability of the static
plane-parallel damping-wing slab (Neufeld 1990 eq. 3.27, Harrington 1973 eq. 40)
evaluated in line every step, or with an externally computed J_Lyα(r) profile read from file
(`radiation/lya_rt.f90`; keys `Jlya escape-prob` and `Jlya RT file`). The second of those is
how a real Monte Carlo field enters: the LaRT scattering-rate field is produced offline
(`examples/exhale_to_lart.py`) and read back in.

### 6.3 Ionization and chemistry solver

Taylor uses the KPP kinetic preprocessor over a large stiff ODE network; AIOLOS uses the
C2Ray two-species scheme of SB23, or the reaction network its code has since grown; Xing
folds semi-implicit source terms into the PLUTO RK stages. EXHALE solves a single coupled
*nonlinear algebraic* system per cell, which is what a steady-state code needs, and
everything lives in it together: H, He, the He(2³S) metastable, the ten metal elements, and
the four molecular species. Adding physics means adding rows to that system rather than a
parallel decoupled solver. The molecular rows are normalized by a turnover-rate scale before
they reach the solver, because the He and H blocks otherwise differ in magnitude by six
orders at He/H = 10³ and `hybrd1` does not scale rows itself.

**Frelikh integrates the network in time; EXHALE solves it for equilibrium.** This is the
sharpest methodological split in this memo, and it is worth stating precisely because the
two codes model the same layer. Frelikh splits the fluid equations into an advective part
and source terms, and solves the chemistry in a substep of its own, of Equation (6) they
write that the chemical source term is one "which we solve for in a separate substep of our
routine." The network is then a stiff initial-value problem, and their Section 2.9 says so:

> "As the rates of the reactions vary by several orders of magnitude, the reaction network
> cannot be solved via a straightforward explicit method: the set of equations is said to be
> 'stiff.' We implement the semi-implicit method of solving stiff systems of equations from
> W. H. Press et al. (2007), the details of which are given in Appendix B."

Their Appendix B gives the whole of it: semi-implicit Euler as in Numerical Recipes section
17.5.2, the specific stepper named (`StepperSie`), the linear system solved by LU
decomposition with forward and back substitution, Δt crossed by a sequence of substeps dt_i
with Bulirsch-Stoer polynomial extrapolation to test each one, and the tolerances quoted
(rtol = 1e-7, atol = 1e-4·rtol, h1 = 1e-6, hmin = 0). The Jacobian is analytic: "The
Jacobian ∂f/∂y is calculated analytically for our problem, as are the time derivatives of
the concentrations on the left-hand side of Equation (B6)." Cells are independent, so the
network runs in parallel across the grid.

Set against that, EXHALE has the more careful *formulation* and the weaker *linear algebra*.
Its rows are steady-state production-loss balances (`System_HeH_mol.f90`); there is no path
in the tree that integrates the network in time.
`nonlinear_system_solver/constrained_chemical_equilibrium.f90` takes the unknowns to be
u_k = ln n_k, so no iterate can go negative, following the CEA method of Gordon & McBride
(1994); it writes element conservation as an explicit residual row per element rather than
as a closure, so that HeH+ satisfies the H and the He budget simultaneously; and it takes
n_e from charge neutrality rather than carrying it as an unknown, so neutrality holds
identically at every iterate. Each rung is then handed to MINPACK `hybrd` with a
*forward-difference* Jacobian: the opposite of Frelikh's choice.

The reason EXHALE needs the extra apparatus is stated in that file's own header: the
molecular cell is bistable, with "a molecular basin in the dense, shielded base and an
atomic basin above the H2 → H front," so a single starting guess either lands in the right
basin or does not. The fix is a natural-parameter continuation in the radiation field,
tracked from six decades below the cell's own field (where the composition is known in
closed form) up to the true field, bisecting any rung that fails.

**The judgement worth recording is that Frelikh does not meet this problem at all, and not
because its solver is better.** A time-integrated network has no uniqueness question in it:
the basin is selected by the previous step's state, not by a guess. Their only defense
against a composition leaving the physical domain is a hard floor, stated in their section
2.5 as "We implement a number density floor of zero for all species," and they record that
the floor leaves a visible mark: their Figure 4 caption attributes a kink in the H2 density
profile to "applying the density floor in a region where the profile drops off steeply."
EXHALE's continuation is, in effect, a reconstruction inside an equilibrium formulation of
the path that time integration supplies for free. Neither approach is wrong; the point is
that the bistability is a cost of choosing equilibrium, not a fact about the physics.

### 6.4 The molecular reaction network

Frelikh's Table 1 and EXHALE's `lower_atmosphere/mol_rates.f90` are the only two molecular
H-He networks among the five that can be compared reaction by reaction, and they largely
agree. Frelikh lists 22 reactions compiled from Yelle (2004) and García-Muñoz (2007); EXHALE
carries 23 transcribed from Koskinen et al. (2022) Table 1, which was checked here against
the published table rather than against the transcription. All 22 of Frelikh's arrows are
one-directional; there is no reversible pair notation in the table.

**Eleven are identical in value and in source**: H⁺ and He⁺ radiative recombination (Storey
& Hummer 1995), H2+ dissociative recombination (Auerbach et al. 1977), H2+ + H2 (Theard &
Huntress 1974), H2+ + H (Karpas et al. 1979), H⁺ + H2(v≥4) (Yelle 2004), He⁺ + H2 → HeH⁺
(Schauer et al. 1989), HeH⁺ + H2 (Bohme et al. 1980), HeH⁺ + H (Karpas et al. 1979), HeH⁺
dissociative recombination (Yousif & Mitchell 1989), and the three-body H3+ association. The
last of these is worth naming because it is easy to assume it is missing: Frelikh's **k21**,
H⁺ + H2 + M → H3+ + M at 3.2e-29 cm⁶ s⁻¹ (Miller et al. 1968), **is present in EXHALE** as
R13, `rk_R13_Hp_H2_M`, at the same value and from the same source.

Yelle (2004) is now held, and each of these ten was compared against his Table 1 rather than
against Frelikh's transcription of it: R15, R6, R8, R9, R10a, R11, R12, R17, R13 and R14 of
that table are the values `mol_rates.f90` carries as R5, R8, R9, R10, R20, R18, R19, R16, R1
and R2. Every one matches, including the exponential of R10 (Yelle prints
1.0e-9 e^(-2.19E4/T), the module 1.0e-9 exp(-21900/T)).

**Six differ**, and the table below gives the sizes. Five are choices of source table. The
sixth is not.

| Reaction | Frelikh | EXHALE | Ratio |
|---|---|---|---|
| H3+ + e → H2 + H | 2.9e-8 (300/T)^0.65 (Sundström 1994) | 2.16e-8 (300/T)^0.65 (0.30 × Larsson 2008) | 1.34 |
| H3+ + e → 3H | 8.6e-8 (300/T)^0.65 (Datz 1995) | 5.04e-8 (300/T)^0.65 (0.70 × Larsson 2008) | 1.71 |
| *H3+ recombination, total* | 1.15e-7 (300/T)^0.65 | 7.20e-8 (300/T)^0.65 | **1.60 at every T** |
| H3+ + H → H2+ + H2 | 2.0e-9, **no barrier** (Yelle 2004, whose own entry reads "Estimated") | 2.1e-9 e^(-20000/T) (Harada 2010) | 4.6e8 (1000 K); 4.6e6 (1300 K); 2.1e4 (2000 K); 748 (3000 K) |
| H + H + M → H2 + M | 8.0e-33 (300/T)^0.6 = 2.45e-31 T^-0.6 (Ham et al. 1970, via Yelle 2004 R5) | 2.8e-31 T^-0.6 (Cohen & Westberg 1983) | 1.14 at every T |
| H2 + M → H + H + M | 1.5e-9 e^(-4.8e4/T) | k_rec/K_eq, detailed balance of R15 (*since 2026-09-01; was 1.5e-9 e^(-48350/T)*) | EXHALE/Frelikh = 0.079 (1000 K); 0.19 (1300 K); 0.51 (2000 K); 0.87 (3000 K). Was 0.705 (1000 K); 0.839 (2000 K) |
| He⁺ + H2 → H⁺ + H + He | 8.8e-14 (Schauer 1989) | 1.0e-9 e^(-5700/T) (Moses & Bass 2000) | 38 (1000 K); 142 (1300 K); 657 (2000 K) |

**H3+ + H is a physics error in Frelikh, not a source choice.** The reaction
H3+ + H → H2+ + H2 is endothermic. From the standard enthalpies of formation:
1488.3 kJ mol⁻¹ for H2+ (which is IE(H2) = 15.426 eV), 1107 for H3+ and 218.0 for H, the
reaction enthalpy is ΔH = 1488.3 − 1107 − 218.0 = 163.3 kJ mol⁻¹ = 1.69 eV, i.e. a barrier
of 19,600 K. The Harada et al. (2010) form that EXHALE carries, e^(-20000/T), matches that
endothermicity; a temperature-independent 2.0e-9 does not, and there is no v-excited label
on the entry as there is on the neighbouring H⁺ + H2(v≥4) row. **The judgement is that
EXHALE is right here and Frelikh is not**, at the four to eight orders of magnitude the
table gives for the molecular layer. Yelle (2004) has since been read, and it removes the
caveat the previous version of this paragraph carried: the transcription is faithful and the
value is not a measurement. Yelle's Table 1 row R7 is `H3+ + H -> H2+ + H2`, `2.0 × 10−9`,
and the reference column of that row reads **"Estimated"**, exactly as his R9 does. Frelikh
copied a barrier-free estimate, and the estimate itself is what the endothermicity above
contradicts.

**The H3+ recombination difference is a later measurement, not a preference.** Neither
2.16e-8 nor 5.04e-8 is printed anywhere as it stands: `mol_rates.f90` builds them from the
two numbers Larsson, McCall & Orel (2008), Chem. Phys. Lett. 462, 145, p. 149 give, the
thermal rate constant α(300 K) = (7.2 ± 1.1)e-8 cm³ s⁻¹ of the Kokoouline and Greene
calculation, "in good agreement with the new storage ring results", and the three-body
branching ratio 0.70 ± 0.07, "in very good agreement with the CRYRING storage ring results".
The temperature index 0.65 is not Larsson's; it is inherited from the Sundström and Datz
fits that Yelle prints, so above about 1000 K the shape of both codes' rates rests on the
same unmeasured extrapolation and only the 300 K normalization separates them. That
normalization is settled in the source: Larsson's p. 149 states that "the early results
obtained at CRYRING [23,24] and ASTRID [27], which gave results just above or at
10−7 cm³ s⁻¹, were slightly too high because of rotational excitations", and his Ref. 24 is
Sundström et al. **Frelikh's pair is the superseded one**, by the rotationally cold storage
ring measurements and the calculation that reproduces them. The assertions that pin this
construction are `src/tests/physics_probe/h3p_recombination_branching.f90`.

**On H2 + M, Baulch et al. (1992) is held too**, and it sharpens both columns. Its
`H2 + H2 -> 2H + H2` entry, p. 550, is `k0 = 1.5e-9 exp(-48350/T)` over **2500-8000 K**,
uncertainty ±0.5 in log10, with a separate argon entry `3.7e-10 exp(-48350/T)`. Yelle's
`e^(-4.8E4/T)` is that exponent rounded, and the value Frelikh carries through him is a
fit whose stated range begins at 2500 K, above the whole molecular layer, and which is
specific to H2 as the collider.

**Nothing in Frelikh's network is absent from EXHALE.** Its k1, k3, k4, k14 and k22 are
photoreactions, and EXHALE has H, He and H2 photoionization (`util_ion_eq.f90`) and
Lyman-Werner photodissociation (`lyman_werner.f90`, opt-in) of its own. It also splits H2
photoabsorption into four mutually exclusive final states in `h2_photo_channels.f90`, so the
dissociative photoionization channel H2 + hν → H + H⁺ + e⁻ is carried separately, with its
own 18.076 eV threshold, rather than folded into H2+ production.

**Six EXHALE reactions have no counterpart in Frelikh**: collisional ionization of H and of
He (Voronov 1997), electron-impact dissociation of H2 (Stibbe & Tennyson 1999), H2 + He⁺
(Barlow 1984), and the two He↔H charge-exchange channels. The last of these is the one that
matters. Frelikh's seven helium reactions (k14-k20) cover photoionization, two He⁺ + H2
channels, two HeH⁺ channels, He⁺ recombination and HeH⁺ dissociative recombination, and
**contain no charge exchange with hydrogen in either direction**. In EXHALE the live
definitions are in `radiation/charge_exchange.f90`, Group B, and the two are not each
other's reverse: B1, He + H⁺ → He⁺ + H, is the non-radiative collisional channel (Kimura et
al. 1993), while B2, He⁺ + H → He + H⁺ + photon, is radiative charge transfer (Stancil et
al. 1998). A photon-emitting process has no collisional reverse, which is why the module
forbids testing the pair against detailed balance.

### 6.5 H3+ infrared cooling

Both codes take the cooling function from Miller, Stallard, Tennyson & Melin (2013), and
EXHALE's 800-1800 K coefficients are the same seven numbers as Frelikh's Equation (13).
Three things differ.

First, **Frelikh radiates into vacuum and EXHALE can exchange with a field**. Frelikh's
volumetric rate is Λ_H3+ = 4π E(T) n_H3+ × 1e7 erg cm⁻³ s⁻¹, with no incident term, no
stimulated emission and no optical-depth factor; downward emission is treated as lost: "We
treat the radiation that is emitted downward (i.e., toward the lower simulation boundary) as
having left the domain." Optical thinness is verified after the fact in their section 3.4
against ExoCross cross sections. EXHALE's `lower_atmosphere/h3p_cooling.f90` carries the
same emission term but its `h3p_net_cooling_rate` returns
Λ_net = n_H3+ · 4π · s(T, n_H2) · [E(T) − W_dil E(T_rad)] × 1e7, which is exactly zero when
T = T_rad and W_dil = 1. **But W_dil is reached only through the key `Base IR field`, which
defaults to off**, so EXHALE's default behavior is Frelikh's.

Second, EXHALE carries the **non-LTE departure factor** s(T, n_H2) of Miller et al.'s Table
6 and all four published temperature segments (30-300, 300-800, 800-1800, 1800-5000 K);
Frelikh is pure LTE and implements only the upper two, so its H3+ cooling is undefined below
800 K and is set to zero above 5000 K.

Third, and recorded here because it will confuse anyone reading the paper against the code:
**Frelikh's Equations (13) and (14) are printed without their logarithm**, as
E(800 < T < 1800 K) = −62.7016 + 0.0526104 T − … carrying units of W sr⁻¹ molecule⁻¹. A
negative polynomial cannot be a power, and the missing operator is a *natural* logarithm,
not a common one. Two independent checks were run here. The two published segments disagree
at their 1800 K junction by −2.35% read as ln E and by −5.34% read as log₁₀ E, and EXHALE's
transcription of the same Table 5 records the mismatch as "2.4% below": the natural-log
value. And a two-level estimate for the ν2 band (A ≈ 100 s⁻¹, E(ν2) = 3628 K, hν at 4 μm)
gives 1.05e-20 W sr⁻¹ molecule⁻¹ at 1000 K, against 3.34e-20 for the natural-log reading and
1.4e-45 for the common-log reading. Their own figures confirm that the code behind the paper
is right and only the printed equations are wrong: their Figure 5 shows
Λ_H3+ ≈ 2e-7 erg cm⁻³ s⁻¹ where Figure 4 shows n_H3+ ≈ a few × 1e3 cm⁻³, which is the
natural-log ratio.

### 6.6 Photoelectron secondary ionization

Section 5 records that EXHALE no longer treats the photoelectron heating fraction as a free
parameter, which separates it from Taylor and Xing. Frelikh does not treat it as free
either, and it uses a different and (for a molecular gas) a better prescription.

EXHALE uses **Shull & van Steenberg (1985)**, with the branching between the H I and He I
secondary channels renormalized from the n(He)/n(H) = 0.1 composition their Monte Carlo
sampled onto the cell's own neutrals (`svs85_secondary_branching`,
`radiation/util_ion_eq.f90`). Frelikh uses **Dalgarno, Yan & Liu (1999)**:
W = W₀(1 + Cxᵅ) from their Tables 4 and 5, a secondary ionization added to the network as an
extra reaction, and (this is the part EXHALE has no equivalent of) the H2 channels put
back explicitly. Their Equations (25)-(35) return the 5.4 eV dissociation heat following
triplet excitation, the thermalized ν = 1 (0.516 eV) and ν = 2 (1.032 eV) vibrational levels
with three separate ionization-fraction weightings, and the Lyman and Werner electronic
states with their ~90% fluorescent / ~10% dissociative branching.

Frelikh also states the range of validity and checks it: "the approximations presented in
A. Dalgarno et al. (1999) are only valid for a gas with an ionization fraction <0.1," which
they justify by locating the relevant τ_ν = 1 surfaces in the neutral layer, with a stated
error of at most 15%.

**The judgement is that for a molecular layer Dalgarno et al. (1999) is the correct
prescription and Shull & van Steenberg (1985) is not.** SvS85 is a Monte Carlo fit for a
pure atomic H + He plasma; it carries no information about where photoelectron energy goes
in a gas that is mostly H2. Checked here, `svs85_secondary_branching` contains no H2 term at
all: the code carries an H2 photoabsorption column separately, but the photoelectron energy
partition has no H2 vibrational thermalization, no dissociation heat return and no H2
secondary ionization. The two prescriptions are complementary rather than competing: SvS85
remains valid to high ionization fraction, where Dalgarno's fits are not, and Dalgarno
covers the molecular gas, where SvS85 has nothing to say. Neither covers the whole domain.

The sign of the effect is also worth recording, since it is not the intuitive one. Frelikh's
section 6 reports that including secondary ionization raised the ionization fraction in the
molecular layer, increased electron-impact destruction of H3+, reduced H3+ cooling and so
*increased* the escape rate, "from 6.4 × 10⁹ to 9 × 10⁹ g s⁻¹ sr⁻¹ for our fiducial hot
Jupiter run, which is counter to the intuition that secondary ionization results in a
decreased efficiency of escape."

### 6.7 Metals and cooling

EXHALE has the richest *closed-form* metal treatment: ten elements (C, O, N, Mg, Si, Ca, Na,
K, S, Fe) in 27 ion stages inside the coupled system, CHIANTI-fitted line cooling,
ground-term fine-structure levels solved in statistical equilibrium at the local (n_e, n_HI)
rather than in the coronal limit, and line trapping in the eight fine-structure lines of
C I, C II, N II and O I. Taylor can include solar-abundance metals but treats the structure
as H/He-dominated. Xing has **no metals** and **only Lyα cooling**. Frelikh has no metals
either, and only Lyα and H3+ cooling; it checked recombination cooling in postprocessing
against the Osterbrock & Ferland rate and left it off as negligible. Its section 5.4 names
the omission as the main limit on extending the model to smaller planets, since H3+ is
destroyed by H2O and CO and metals bring line coolants of their own. SB23 itself carries
only recombination and Lyα cooling; the AIOLOS code adds the Black (1981) rates with
free-free and H3+ cooling, but neither is specialized to the metal-line escape problem.

### 6.8 Line diagnostics

Taylor computes 10830 and Hα in loop with the most complete non-LTE microphysics. Xing
post-processes 10830 only. EXHALE builds the 2³S and H(n=2) populations in the run (the
metastable is on by default) and produces the transmission spectrum (10830, Lyα, Hα, Hβ,
the Mg II, Ca II and Na I doublets) in `EXHALE_transit.py`, with impact-parameter Voigt
integration, instrument and rotation convolution and an optional triaxial Roche geometry.
AIOLOS does not target these lines.

### 6.9 Lower boundary

Taylor couples to a photochemical lower and middle atmosphere at ~1 μbar, which is the most
physically motivated boundary of the five. Xing uses a fixed-density, fixed-temperature
base. AIOLOS leaves the boundary configurable. EXHALE now has the same range as Taylor: a
fixed base, an analytic chemical-equilibrium column, a scalar handoff, or a full profile;
item (D) below.

**Frelikh imposes the same count as EXHALE, and one thing more.** Its base sits at
P₀ = 0.96 dyn cm⁻², i.e. at "the region of the outflow above the 1 μbar level, near which
molecules such as H2 are present," with T₀ = 1000 K and Rp = 1e10 cm for the fiducial hot
Jupiter (their Table 2). Their section 2.8 states what is imposed and why:

> "At the inner boundary, the flow is subsonic, so we specify the boundary conditions for
> the number densities of each species and the pressure. We extrapolate the velocity from
> the grid: u0 = u1 − (r1 − r0)(u2 − u1)/(r2 − r1), a wave can travel from the domain to the
> inner boundary, as the flow is subsonic, so the value of the velocity is part of the
> solution and cannot be imposed as a boundary condition."

EXHALE's `BC_component_constrho` pins ρ and p and copies the ghost velocity from the first
interior cell. **The count of imposed conditions and the choice of which component floats
are identical**, and Frelikh gives the characteristic argument for that choice explicitly.
**The difference is that Frelikh also imposes the composition.** That is the one item on
EXHALE's own outstanding list for the base that Frelikh has already done; it appears as item
(M) in section 9.

Two further points. Frelikh nowhere reports a supersonic base, nor a guard against one.
Searching the paper for that failure mode returns nothing: the subsonic condition is a
*premise* of the section 2.8 argument and the possibility of its being violated is not
discussed. The structural weakness that follows (once the inflow is supersonic all three
characteristics enter the domain and an extrapolated velocity stops being a boundary
condition) applies to both codes; only EXHALE has found it.

Second, Frelikh's answer to the base thermal problem is not a physical model but a
deliberate override, and the paper says so. Their section 2.2 opens by naming the symptom:
"We implement a bolometric heating and cooling term in this region, as the EUV heating is
insufficient to prevent the gas from cooling to unphysically low temperatures", and then
applies, below the τ = 3 surface of the highest-energy bin,
Γ_bol − Λ_IR = κ σ_SB (T₀⁴ − T⁴) with the opacity set by hand:

> "where we have set κopt and κIR to a single inflated value κ = 1 cm² g⁻¹. In addition, in
> our implementation, we set T₀⁴ = 4T_eq⁴, such that Equation (12) has the form of a
> Newtonian cooling term. … our choices effectively force the bolometrically heated layer to
> be isothermal at the temperature given by the input parameter T₀."

For the fiducial planet that region is r < 1.005 Rp. Heat conduction plays no part in it:
Frelikh's section 3.5 finds conduction "subdominant at all fluxes in our grid of models,"
worth ~2% in Ṁ at 0.3 au, and the only stability remark about it in the paper runs the other
way, their section 2.8 warns that "conduction makes the choice of outer boundary conditions
more important, as the conduction term does not explicitly include a flux limiter."

The bearing on EXHALE is that EXHALE already owns a less blunt instrument for the same
problem. `Base IR field` gives the H3+ band a radiative equilibrium fixed point at
T_rad = 1140 K, W_dil = 1/2, i.e. 975.4 K, instead of letting it run the layer down; it is
off by default (section 6.5). The reservation is that the temperature trough measured in the
molecular gate cases reaches 865 K, *below* that fixed point, so H3+ alone does not account
for it and the other channels in that region have not been separated.

### 6.10 The molecular layer, measured side by side

Frelikh's fiducial hot Jupiter and EXHALE's molecular regression cases are close enough in
irradiation to compare directly. Frelikh quotes F_EUV = 2.6e3 erg cm⁻² s⁻¹ at 0.05 au;
`EXHALE_setup.out` for the gate cases reports log₁₀ F_XUV = 3.193, i.e.
1.56e3 erg cm⁻² s⁻¹, a factor 1.7. The planets differ more: 0.7 M_J against a 14.52 M_⊕ hot
Uranus (Rp = 0.49 R_J, T_eq = 1140 K, a = 0.048 au, He/H = 0.0793).

All EXHALE numbers below were measured here from `backup/regression/*/output/` by parsing
the `# columns` headers; no case was re-run. Frelikh's numbers are from its Table 2 and its
Figure 4, 5 and 9 captions, except the Λ_H3+ plateau, which was read off the Figure 5 axis
and is therefore approximate.

| | **Frelikh HJ** | **Frelikh s-Earth** | **EXHALE `mol_metals`** | **EXHALE `mol_base_handoff`** |
|---|---|---|---|---|
| Planet mass | 0.7 M_J | 7.2 M_⊕ | 14.52 M_⊕ | 14.52 M_⊕ |
| Base pressure | 0.96 μbar | 0.47 μbar | 9.0 μbar | 9.0 μbar |
| Base T | 1000 K (forced) | 1000 K (forced) | 1213 K | 1213 K |
| Base H2 fraction of H nuclei | ≈1 | ≈1 | 0.99989 | 0.9996 |
| H2→H crossover | 1.01 Rp | no sharp transition | 1.158 Rp (1166 K) | 1.158 Rp (1155 K) |
| Molecular layer cut-off | 1.04 Rp, sharp | **none; H2 survives throughout** | 1.19 Rp, sharp | 1.19 Rp, sharp |
| Molecular layer T | 1000→4000 K | - | 1213→1319 K | 1213→1305 K |
| T minimum / maximum | - / ~9500 K | - / 4200 K | **865 K** at 1.190 / 2719 K at 1.910 | 874 K / 2647 K |
| n(H3+) maximum | ~2e5 cm⁻³ | - | 46.3 cm⁻³ at 1.124 Rp | 6.7e4 cm⁻³ at 1.031 Rp |
| N(H3+) | 1.0e13 cm⁻² | 1.5e13 cm⁻² | 1.73e10 cm⁻² | 1.79e13 cm⁻² |
| Λ_H3+ maximum | ~2e-7 (read off Fig. 5) | - | 6.145e-10 at 1.113 Rp | 8.898e-7 at 1.031 Rp |
| … as a fraction of local cooling | most of it | dominant | 6.8% | 100% |
| … as a multiple of local heating | ≲1 | - | 0.003 | **6.95** |
| Sonic point | 3.1 Rp | 4.4 Rp | 2.878 Rp | 2.874 Rp |
| Ionized fraction at sonic point | mostly ionized | mostly neutral | 0.41 | 0.41 |
| log₁₀ Ṁ [g s⁻¹] | 11.05 if integrated over 4π | - | 10.582 | 10.575 |
| Peak Mach number below 1.35 Rp | not reported | - | 0.017 | 0.014 |

Three readings of that table matter.

**The three-layer structure reproduces, but on the wrong side of a qualitative divide.**
Frelikh's section 5.4 sets its own super-Earth against Koskinen et al. (2022) and finds them
alike, and both unlike a hot Jupiter:

> "Our super-Earth runs are qualitatively similar to the results of their hot Uranus model,
> in which instead of the sharp transition between the molecular hydrogen layer and the
> atomic layer, as seen in our hot Jupiter simulations … molecules are present throughout
> the outflow. Particularly interesting is the presence of H3+ throughout the outflow,
> instead of being confined to a narrow region at the base."

EXHALE's hot Uranus gate does the opposite. Measured here, the H2 fraction of H nuclei falls
from 0.48 at 1.16 Rp to 2.7e-5 at 1.30 Rp, and n(H3+) from 35 to 1.2e-3 cm⁻³ over the same
interval, four decades each. The temperature is also low: 865-2719 K, against the
4000-5000 K Frelikh quotes as the Koskinen hot Uranus prediction and the 1000-3800 K of its
own super-Earth. **EXHALE's hot Uranus is behaving like Frelikh's hot Jupiter, not like the
Koskinen hot Uranus the gate case is meant to reproduce.**

**The two EXHALE gate cases disagree with each other by three orders of magnitude** in
N(H3+), 1.73e10 against 1.79e13 cm⁻², and only the second is in Frelikh's range. The
difference is the `base.inp` handoff, whose q_H2_base sets the base H3+. In that case
Λ_H3+ is 6.95 times the local heating at 1.031 Rp and supplies 100% of the cooling there,
which is a net cooling excess visible in the numbers.

**The base flow is not supersonic, and that is not the same as saying it is healthy.** The
peak Mach number below 1.35 Rp is 1.4-1.7e-2 in every gate case, and the He/H of these runs
(0.0793) is far outside the composition window where a supersonic base has been seen. What
the same output does show is described in section 8.

## 7. A defect neither code escapes: the H2 thermal-dissociation pair

**Repaired on the EXHALE side, 2026-09-01; see the end of this section.** The diagnosis
below is left as it was written, because it is what identified the defect, and because it
still holds for Frelikh. The table's EXHALE column describes the code as it was before
that date.

Frelikh's k12 and k13 are the forward and reverse of the same process, H2 + M ⇌ H + H + M,
taken from two different measurements: 1.5e-9 e^(-4.8e4/T) from Baulch et al. (1992) for the
dissociation and 8.0e-33 (300/T)^0.6 from Ham et al. (1970) for the three-body
recombination. EXHALE's R12 and R15 were the same pair from the same two sources, differing
only in that Koskinen et al. (2022) print the barrier as 48,350 K where Frelikh rounds it to
4.8e4 K.

Nothing in Frelikh discusses whether the pair is thermodynamically consistent. The paper was
searched here for *detailed balance*, *thermodynamic*, *equilibrium constant*, *reverse
reaction*, *reversible*, *extrapolat* and *range of validity*, across the body, the
appendices, the footnotes and the table notes. The only hits are the linear velocity
extrapolation of section 2.8 and the Bulirsch-Stoer extrapolation of Appendix B, both
unrelated. **There is no statement about the consistency of the pair, and none about the
temperature range over which either fit is valid.** No such record was found on the EXHALE
side either; the caveats in `mol_rates.f90` concern third-body efficiency and the
helium-rich limit, not forward-reverse consistency.

The ratio k_diss/k_rec is the equilibrium constant n_H²/n_H2, with the third body
cancelling, so it can be checked against thermodynamics directly.[^detbal]

[^detbal]: K = (2πμk_B T/h²)^(3/2) · g_H²/z_int(H2) · e^(−D₀/k_B T) with μ = m_H/2, g_H = 2
and D₀ = 36,118.11 cm⁻¹. The internal partition function was summed explicitly over
v = 0-15 and J = 0-79 with the Huber & Herzberg constants (ω_e = 4401.21, ω_e x_e = 121.33,
ω_e y_e = 0.812, B_e = 60.853, α_e = 3.062 cm⁻¹), discarding levels above D₀, with
ortho/para nuclear-spin weights 3 and 1 divided out by (2I+1)² = 4 so that the nuclear spin
cancels consistently against the H atoms. A rigid-rotor harmonic-oscillator evaluation of
the same expression agrees to 2.8% at 1000 K and 1.7% at 2000 K.

Ratios above unity mean the pair drives the gas more atomic than equilibrium allows.
Computed here.

| T [K] | K_thermo [cm⁻³] | Frelikh k12/k13 | ratio | EXHALE R12/R15 (pre-2026-09-01) | ratio |
|---:|---:|---:|---:|---:|---:|
| 500 | 7.07e-22 | 5.17e-19 | **732** | 2.57e-19 | **363** |
| 800 | 7.71e-05 | 2.96e-03 | **38.4** | 1.91e-03 | **24.8** |
| 1000 | 3.80e+01 | 5.50e+02 | **14.5** | 3.88e+02 | **10.2** |
| 1500 | 1.54e+09 | 6.24e+09 | 4.05 | 4.94e+09 | 3.21 |
| 2000 | 9.96e+12 | 2.21e+13 | 2.22 | 1.86e+13 | 1.86 |
| 3000 | 6.38e+16 | 8.40e+16 | 1.32 | 7.48e+16 | 1.17 |
| 4000 | 4.97e+18 | 5.45e+18 | 1.10 | 4.97e+18 | 1.00 |
| 5000 | 6.62e+19 | 6.87e+19 | 1.04 | 6.40e+19 | 0.97 |

**The judgement is that the pair is thermodynamically consistent only above about
3000-4000 K, and that this is a defect of both codes equally.** Below that it drives the
composition toward atomic hydrogen faster than equilibrium permits, by a factor 14.5 in
Frelikh and 10.2 in EXHALE at 1000 K, and by 38 and 25 at 800 K. EXHALE is the less wrong of
the two only because the higher barrier it transcribes partly compensates; the error is the
same error. Both codes place their molecular layer in exactly the 1000-1500 K band where the
inconsistency is largest: Frelikh forces its base to T₀ = 1000 K, and the EXHALE gate cases
measure 1213-1319 K through the layer.

Two qualifications. Neither code sets the H2/H balance from this pair alone:
photodissociation, the Lyman-Werner channel and the ion chemistry all act, so the
composition error is smaller than the tabulated factors and is confined to wherever thermal
dissociation and three-body recombination dominate. And this is a statement about the pair,
not about either published measurement in its own range of validity; the Baulch fit belongs
to shock-tube temperatures far above the molecular base.

### 7.1 The repair, 2026-09-01

EXHALE no longer transcribes the dissociation fit. `rk_R12_H2_thdis` in
`src/modules/lower_atmosphere/mol_rates.f90` now returns `k3b_H_H_to_H2(T) /
keq_H_H_to_H2(T)`, i.e. the three-body recombination coefficient divided by the equilibrium
constant the module builds from the H2 spectroscopic constants. R15 at the same time moved
from the Ham et al. (1970) room-temperature measurement to the Cohen & Westberg (1983)
recommended k₁(H2) = 2.8e-31 T^-0.6, whose stated 50-5000 K range covers the molecular
layer where Ham's 77-300 K does not.

The consequence for the table above is that the EXHALE ratio column becomes **1 by
construction at every temperature**: R12/R15 is now 1/K_eq exactly, and against the
independently computed K_thermo of the footnote it is 0.998-1.000 over 500-5000 K, the
residual being the two evaluations' different treatment of the H2 internal partition
function. The 363, 24.8, 10.2 and 3.21 of the 500-1500 K rows are gone.

Two checks that the construction is right rather than merely self-consistent. Against the
Baulch fit it replaces, the new k_diss is 1.05 at 8000 K, 1.18 at 5000 K, 0.97 at 3000 K and
0.82 at 2500 K (agreement inside Baulch's own stated 2500-8000 K range) and departs only
below it (0.61 at 2000 K, 0.11 at 1000 K), which is where the fit is an extrapolation.
Against the primary shock-tube measurement rather than the evaluation, Breshears & Bird's
own M = H2 coefficient 5.48e-9 exp(-52989/T), the ratio is 1.11 at 3500 K (the bottom of
their measured range), 0.99 at 4000 K, 0.82 at 5000 K and 0.51 at 8000 K.

**The judgement of section 7 therefore stands for Frelikh and no longer for EXHALE.** The
pair is now thermodynamically exact in EXHALE at every temperature, and the code's molecular
layer is no longer driven toward atomic hydrogen faster than equilibrium permits. What the
detailed-balance route assumes instead is the vibrational state the coefficients refer to,
which is recorded at the code site and in `docs/molecular_hydrogen_treatment.tex`.

## 8. What the molecular gate cases show

The comparison in section 6.10 required reading EXHALE's molecular regression output
closely, and five things turned up in it that are not about Frelikh at all. They are
recorded here with the method, so that they can be re-measured. Every one comes from
`backup/regression/<case>/` with no case re-run; the profile quantities are from
`output/Hydro_ioniz.txt` columns r/Rp, ρ [m_H cm⁻³] and v, and the warnings from `run.log`.
All five gate cases are relaxation snapshots pinned at `maxsteps` = 12,000.

1. **The mass flux is not flat at the base.** ρvr² in `mol_metals` runs 1.6563e16 at
   r = 1.00019 and 3.1420e16 at r = 1.00039: a factor **1.90 across one cell**, at the
   second interior cell. Over the molecular layer (1.0 < r/Rp < 1.19) the spread
   (max−min)/mean is **1.865**. In `mol_lyman_werner` the first two rows are exactly zero,
   the third is **negative**, −1.6507e16, and the fourth returns to +1.3357e15: the flux
   changes sign inside the layer. A converged steady state requires ρvr² to be constant in
   radius, and these are not.
2. **Inflowing cells sit inside the molecular layer in one case.** Counting rows with v < 0:
   `mol_metals`, `mol_base_handoff` and `mol_diffusion` each have 28, all at
   r = 1.305-1.426 Rp, i.e. *outside* the molecular front at 1.16-1.19.
   `mol_lyman_werner` has **144, at r = 1.00019-1.12394 Rp, inside the layer**, with the
   first cell's velocity exactly zero: the one-way valve latched shut.
3. **The accepted chemical states that are not roots sit at the temperature minimum.**
   `mol_metals/run.log` carries 24 lines of `(ioniz_eq) NON-ROOT accepted (relaxation
   amnesty)` and `mol_ir_bands` 11; an example reads `cell 260 step 1 r 1.197620 T[K]
   1.197E+03 info 4 viol 1.210E-08 res 1.065E-03`. They cluster at r ≈ 1.19-1.20 Rp, which
   is where the temperature trough bottoms out at 865 K. **The direction of causation was
   not determined.** Whether a runaway in the cooling breaks the chemical solve, or a state
   that is not a root drives the cooling, or a third cause produces both, is not decided by
   this measurement, and no inference is drawn here. What the coincidence does put in
   question is whether the thermal collapse and the non-uniqueness of the chemical solve are
   two problems or one.
4. **Lyman-Werner self-shielding is being evaluated outside its fit.**
   `mol_lyman_werner/run.log`: "the H2 Lyman-Werner lines remove 1.01E+00 of the 912-1110 A
   band, i.e. all of it. The Draine & Bertoldi (1996) self-shielding fit is a fit to a RATE
   and does not bound the photon budget; the band transmission is floored at zero and this
   column is outside the fit." The band is 101% removed and the transmission floors at zero.
   This is the same case as item 2; no connection between the two is asserted.
5. **`mol_ir_bands` exercises one of the three channels its name implies.** With
   `Base IR field: True` and `Molecular IR bands: True`, the H2O and CO infrared columns of
   `output/Cooling_breakdown.txt` (col. 13 and 14) are **exactly zero at every radius**,
   because `Oxygen chemistry` is off in that configuration and the carriers do not exist.
   The only molecular infrared channel actually switched on is the H2 band, at 8.98e-9 at
   the H3+ peak cell and 2.93e-7 erg cm⁻³ s⁻¹ at its own maximum. As a regression gate the
   case does not protect the H2O and CO paths.

A note on quoting Ṁ from these cases, since it is easy to get wrong. `src/EXHALE_main.f90`
evaluates it at j = N − 20, i.e. at r = 3.9231 Rp, and halves it under
`2D approximate method: Rate/2 + Mdot/2`. Recomputing 4π ρ m_H v r² at that index and
halving reproduces the log₁₀ Ṁ = 10.58 that `EXHALE_setup.out` prints. Evaluating at the
outer edge instead, or omitting the factor of two, gives 10.94 and is not the quantity the
code reports. In any case the mass-flux spread above means none of these numbers should be
quoted as a converged Ṁ.

For completeness: the temperature floor these cases could have hit is 0.01 T₀ = 11.4 K,
because `T_trial` in `energy_semi_implicit.f90` is dimensionless. No case activated it
(`grep -c floor` over `run.log` returns zero in all five), and the measured minima (865,
874, 946 and 970 K) differ case by case rather than resting on a common constant.

## 9. Adoption assessment

Each candidate is given one of three verdicts. **Adopted** names the module that implements
it. **Adopted then replaced** records that the element entered the code and was later
superseded, and by what. **Not adopted** states whether the original reason still holds.

### 9.1 Adopted

**(A) He/H diffusive separation: adopted, in a stronger form than recommended.**
The recommendation was Taylor's route: a molecular-plus-eddy diffusion term on top of a
single bulk momentum solve. That is what `functions/binary_element_diffusion.f90` does, with
the two departures described in section 6.1: the transported variable is the helium mass
fraction of a binary mixture rather than a trace ratio, and the binary coefficient is
resolved by ionization stage rather than taken at a neutral rate. The eddy term comes from
the scalar key `He_Kzz` or, when a lower-atmosphere profile is in use, from that profile's
K_zz(p) interpolated onto the grid. One implicit backward-Euler step is taken per hydro step
and converged by Newton on a tridiagonal Jacobian, with a Péclet-switched donor-cell
treatment of the X(1-X) drift product so the flux vanishes where either component runs out.
Metals ride with hydrogen at fixed metal/H unless `He_metal_diffusion` is set, in which case
each metal additionally diffuses as a trace species through the solved hydrogen background.
The operator is available both inside the marching loop and, through an outer composition
relaxation, on the direct steady route. Default off (`He_diffusion`).

**(D) Whole-atmosphere lower boundary: adopted, and carried past Taylor.**
EXHALE now has three routes below the base, in increasing order of what they carry. An
analytic isothermal column from the 1-bar level to the μbar base using the Koskinen et al.
(2022, Eqs. 11-13) chemical-equilibrium H2/H/He partition
(`lower_atmosphere/lower_column.f90`). A scalar handoff file `base.inp`, generated
analytically or from a VULCAN run, which carries the EOS boundary (T, r, p, q_H2), the
elemental reservoirs (He/H and a `<El>_H_base` key for each of the ten elements, which is
what `melem_ab` is then built from) and the eddy coefficient. And a full *profile* named by
`Lower atmosphere profile:` (`files_IO/lower_atmosphere_profile.f90`), which carries p, r,
T, n_tot, ρ, K_zz, q_H2, q_H, X_He and a mixing-ratio column per element, interpolated in
log p. A scalar physics key sitting beside a profile is refused rather than merged, and the
two files must agree on a `solution_id`.

The part that goes beyond Taylor is the closure. Taylor's coupling is one-way: the
photochemistry sets the base and the wind reads it. EXHALE closes the loop on the elemental
fluxes. `binary_element_diffusion.f90` measures F_H and F_He face by face and reduces them
to a median and a spread over a radial window, and `src/utils/element_flux_closure.py`
iterates Φ_El = F_El(Φ_El) to a fixed point with under-relaxation, the chemistry taking Φ as
its upper boundary condition and the wind returning F. The elemental composition of the wind
is then an output of the coupled pair rather than an input to either.
`photochem_to_lower_profile.py` and `vulcan_to_lower_profile.py` produce the profile file; a
residual smaller than the spread of the window it was measured on is reported as unresolved
rather than as converged.

**(E) Molecular chemistry and H3+ cooling: adopted.**
H2, H2+, H3+ and HeH+ are unknowns of the same coupled system as H, He and the metals
(`System_HeH_mol.f90`, `System_HeH_mol_metals.f90`), so the molecular and metal blocks share
one electron density. H3+ infrared cooling follows Miller et al. (2013), in LTE with the
non-LTE departure factor and with absorption of the ambient infrared field so the rate
returned is the net one (`lower_atmosphere/h3p_cooling.f90`). H2 photodissociation in the
Lyman-Werner bands uses Draine & Bertoldi (1996) self-shielding
(`lower_atmosphere/lyman_werner.f90`, key `Stellar LW flux`, zero by default). The
three-body density M of the association reactions is the total gas-particle density from
`calc_ntot`, not a mass density. Default off (`Molecular chemistry`); turning it on requires
helium and forces the molecular particle count at the base.

**(F) Lyα Monte Carlo: adopted as an option, not as the default.**
The recommendation was to keep the escape-probability closure and use Taylor's Monte Carlo
only as a validation benchmark, because a WASP-121 b column with τ₀ = 1.2×10⁸ makes an
acceleration-free Monte Carlo intractable and forbids core skipping. That still stands for
the in-loop field: the default is the in-line static-slab closure. But EXHALE also
accepts a tabulated J_Lyα(r) from any external solver (`Jlya RT file`), and that path is
used with the LaRT Monte Carlo code, including the Lyα emitted in situ within the wind. The
Monte Carlo is therefore available where it is affordable, without being in the iteration.

**(I) Navier-Stokes viscosity and heat conduction: adopted (from Koskinen, i.e. Taylor's
own hydro heritage).**
This was not on the original list. `time_step/viscous_conduction.f90` adds the viscous
momentum force, its dissipation and thermal conduction (the Navier-Stokes level that CETIMB
solves and that an inviscid HLLC scheme lacks) integrated Crank-Nicolson and entering the
steady residual through the same operator. Both are opt-in (`Viscosity`, `Conduction`). Note
that the module does *not* implement Koskinen et al. (2022) Eqs. (B5) and (B6) as printed:
for a radial flow with the Newtonian stress τ_ij = μ(∂_i w_j + ∂_j w_i) - (2/3)μ δ_ij ∇·w,
the viscous force and the dissipation are

    F_mu = (4/3)(1/r²) d/dr( r² μ dw/dr ) - (4/3)(dμ/dr)(w/r) - (8/3) μ w/r²
    q_mu = (4/3) μ ( dw/dr - w/r )²

and the code carries these rather than the printed forms; the derivation is in the module
header.

**(J) AIOLOS's opacity-model dispatcher: adopted.**
Also not on the original list, which treated AIOLOS as a cross-check only.
`radiation/opacity_models.f90` is an AIOLOS-style dispatcher over photoionization opacity
models (analytic, constant, constant with Robinson & Catling pressure broadening, and
tabulated from a two-column `.opa` file) selected at runtime by the presence of
`opacity.inp`. The pressure factor that was a deferred hook in the version it came through
is applied cell by cell here.

### 9.2 Adopted, then replaced

**(B) Temperature-dependent Penning ionization: adopted from Taylor (2025), then retired.**
The original recommendation was the cheapest item on the list, and it was taken: the
constant 5×10⁻¹⁰ cm³ s⁻¹ of Roberge & Dalgarno (1982) was replaced by the two-branch power
law of Taylor et al. (2025) Table 2, unconditionally. It has since been withdrawn. The fit
was transcribed correctly but is not usable: it **jumps by a factor 1.5724 at 4000 K**
(1.5849×10⁻⁹ → 2.4921×10⁻⁹), and a rate coefficient is a continuous function of temperature.
It also disagrees with its own paper's Figure 19, where the Maxwell-Boltzmann average of the
Cohen & Lane (1971) and Morgner & Niehaus (1979) cross sections rises to 1.50×10⁻⁹ near
1400 K and falls to 1.29×10⁻⁹ at 4000 K while the tabulated branch decreases monotonically,
and above 4000 K it exceeds every cross-section determination collected by García Muñoz
(2025).

What is in the code now is the continuous closed form of García Muñoz (2025), *A&A* 698,
A199: the Maxwell-Boltzmann average of the Movre & Meyer (1997) total ionization cross
sections for the atomic partner (`ioniz_HeI23S_H`) and of Cohen & Lane (1977) for the
molecular one (`ioniz_HeI23S_H2`), both in `radiation/Cool_coeff.f90`. Both return the
*total* rate that removes the metastable; the 0.9/0.1 split between the Penning channel,
which leaves a lasting ion and a free electron, and the associative channel, which makes
HeH+, is carried by `f_penning_HeI23S` at the points where the products matter. Taylor's own
Figure 19 agrees with the adopted expression to 10-20%, so the two independent calculations
are consistent and it was the published fit that was the outlier.

**(K) AIOLOS's C/N/O line-cooling fits: adopted, then replaced by CHIANTI.**
The analytic C/N/O metal-cooling fits of AIOLOS's `photochem.cpp` (code, not paper) were
ported into `Cool_coeff.f90` and were the default for a time. They are now the legacy branch
(`cno_cool 0`) and the default is a set of closed-form fits to CHIANTI, accurate to 0.2-2.1%
over 10³-10⁵ K. The reason is physical, not stylistic: the AIOLOS fits contain no nitrogen
cooling at all, and their O I/O II rates deviate from CHIANTI by 40-70%. The AIOLOS branch
is kept so old results remain reproducible, and its remaining tuning constants are flagged
in place.

### 9.3 Not adopted

**(C) The B-spline K-matrix He(2³S) cross section, not adopted; the reason still holds, but
the earlier statement of it was too generous.**
It was claimed that EXHALE was "already most of the way there." What `functions/cross_sec.f90`
actually has is a smooth background: two Verner et al. (1996) single-shell forms joined by a
log-linear bridge, fitted at the few-percent level to the Norcross (1971) / TOPbase
background, with the autoionization region represented by a *resonance-averaged* bump.
Taylor's cross section resolves the individual EUV resonances, which EXHALE's does not and
cannot. The case for adopting it is unchanged: low cost, a data swap; moderate value, since
whether the resonances matter depends on where the star's coronal lines fall, but it should
be stated as a missing capability rather than as a finished job.

**(G) AIOLOS's FLD radiation and general reaction network, not adopted; the reason still
holds, with one clause of it withdrawn.**
FLD matters for optically-thick interiors and dust, not for an optically-thin thermosphere,
and a general implicit reaction matrix is the wrong shape for a steady-state code whose
ionization system is one algebraic solve. Both still stand. The third clause of the old
argument (that HLLC would be no improvement on EXHALE's "approximate Riemann solver") was
simply wrong about the code: EXHALE has had HLLC since ATES, alongside LLF and Roe, and HLLC
is what every example configuration uses. There was nothing to adopt there.

**(H) Xing's electron-pressure Riemann split and full multi-fluid, not adopted; the reason
is now stronger.**
Converting EXHALE to genuine multi-fluid means velocities per species, collisional drag and
the electron-pressure split, which is a rewrite of the hydro core and hard to reconcile with
a single coupled algebraic ionization solve. The original argument was that Taylor's
diffusion approximation captures most of the same effect far more cheaply. That argument is
now backed by what was built: the stage-resolved binary operator reproduces the friction
physics (an ion held to the protons by Coulomb collisions, a neutral settling at a
hard-sphere rate) that motivated the multi-fluid treatment in the first place. Xing remains
a validation reference for the fractionation itself.

**(L) The Dalgarno et al. (1999) photoelectron partition in a molecular gas, not adopted;
this one is a real gap.** Unlike (C), (G) and (H), the reason here does *not* still hold.
EXHALE uses Shull & van Steenberg (1985) throughout, and its branching routine contains no
H2 term of any kind (section 6.6): no vibrational thermalization, no dissociation heat
returned, no H2 secondary ionization. That was harmless while the code had no molecular
layer. It is not harmless now that `Molecular chemistry` exists, because the layer it omits
the physics from is exactly the layer where the highest-energy photons deposit. The
prescriptions are complementary (SvS85 stays valid where Dalgarno's fits do not, above an
ionization fraction of 0.1), so the shape of the fix is a composition switch between the two
rather than a replacement. Cost: moderate, a set of closed-form fits from Dalgarno et al.
(1999) Tables 4 and 5 plus their Equations (25)-(35). Value: it changes the energy budget of
the molecular layer, and Frelikh reports the effect running counter to intuition, raising Ṁ
rather than lowering it.

**(M) A composition boundary condition at the base, not adopted; recommended.** Frelikh
imposes the number density of every species at the inner boundary along with the pressure
(section 6.9); EXHALE pins ρ and p and lets the interior chemistry decide what the gas at
the base is made of. The imposed-condition count is the same either way, so this is not a
well-posedness question: it is a question of whether an inflowing molecular reservoir is
allowed to arrive with a composition of its own. EXHALE's own outstanding list already
carries this item; the point worth recording is that a peer code does it and reports no
difficulty from it.

**(N) Time integration of the network in molecular cells, not adopted; recorded as an
option rather than a recommendation.** The bistability that forces the radiation-field
continuation in `constrained_chemical_equilibrium.f90` is a cost of solving for equilibrium,
and a short time integration of the same network does not have it (section 6.3). Frelikh's
stepper is an off-the-shelf semi-implicit Euler with an analytic Jacobian, which is not
expensive. Against that: EXHALE is a steady-state code, its molecular rows share one
electron density with the metal and He(2³S) rows, and splitting the molecular block out into
a time-integrated substep would break the single-system architecture that the rest of this
memo records as the right one. The honest statement is that the continuation and the substep
solve the same difficulty and neither has been measured against the other here.

**(O) Post-hoc verification of the radiative-transfer approximations, not adopted; cheap,
and worth doing.** Frelikh does two checks that EXHALE asserts rather than measures. Their
Appendix D runs a Monte Carlo of 1000 Lyα photons through the converged profile and reports
78% escaping to space, 11% lost downward into the bolometric region and 11% thermalized,
which is what licenses the optically thin Lyα cooling rate; repeating it from deeper in the
molecular layer gives 36/27/37 and they say so. Their section 3.4 checks H3+ optical
thinness against ExoCross cross sections at the computed column. EXHALE has the ingredients
for both: the LaRT coupling for the first, the Miller et al. data for the second, and uses
neither as a check on its own closures.

## 10. What was not verified

The Frelikh comparison rests in places on sources that were not available, and on EXHALE
paths that were not followed to the end. Two of the six are now closed, and both the closure
and what remains open are listed so that nothing above is read as more settled than it is.

1. ~~Yelle (2004), Larsson et al. (2008) and Baulch et al. (1992) are not held.~~
   **Closed.** All three were read in the published versions and section 6.4 now rests on
   them. Yelle's Table 1 row R7 gives the barrier-free 2.0e-9 with the reference
   "Estimated", so Frelikh's transcription is faithful and the endothermicity argument
   convicts the estimate itself. The H3+ recombination factor of 1.60 is a superseded value
   over the one that supersedes it, in Larsson's own words. Baulch's H2 + H2 dissociation
   entry is `1.5e-9 exp(-48350/T)` over 2500-8000 K, so Yelle's `4.8E4` is that exponent
   rounded and the fit is used below its stated range by both codes.
2. ~~Whether EXHALE implements the dissociative photoionization channel H2 + hν → H + H⁺ + e⁻
   was not established.~~ **Closed: it does.** `h2_photo_channels.f90` resolves σ_H2 into
   four mutually exclusive final states, of which (S) is H2 + hν → H + H⁺ + e⁻ at a 18.076 eV
   threshold, alongside (M) H2+ + e⁻, (D) H⁺ + H⁺ + 2e⁻ above 51.4 eV and (N) neutral
   dissociation; the four cross sections sum identically to σ_H2. The channel is consumed as
   a source term by `System_HeH_mol.f90`, as an energy recipient by `util_ion_eq.f90` and in
   the band tables of `set_energy_vectors.f90`.
3. **The temperature clamp reported elsewhere was not reproduced.** The five molecular gate
   cases do not activate the 11.4 K floor, and their minima differ case by case (section 8).
   Whatever run exhibited a clamped cell, it is not one of these.
4. **The direction of causation between the temperature trough and the non-root chemical
   states** at r ≈ 1.19 Rp is undetermined (section 8, item 3).
5. **Frelikh contradicts itself on the number of evolved species**, and the memo does not
   adjudicate. Its section 2.9 says "the seven atomic/molecular species (+electrons)"; its
   Appendix A says "the eight component species"; its Table 1 contains eight neutral and
   ionic species (H, H⁺, H2, H2+, H3+, He, He⁺, HeH⁺).

One further observation is recorded as an observation only. In Frelikh's own molecular
layer, reading n(H) ≈ 1e11 and n_e ≈ 3e8 cm⁻³ off its Figure 4, the barrier-free H3+ + H
channel gives ≈200 s⁻¹ against ≈12 s⁻¹ for dissociative recombination, which would make it
the dominant H3+ sink. That sits awkwardly against their section 2.3, which says of the
recombination channels that "we verify post facto that they do indeed dominate in our
results." The numbers behind this were read off a published figure, so it is not asserted as
a finding.

## 11. Summary

- **AIOLOS** is a general-purpose, multi-fluid RHD code with the most sophisticated
  *numerics* and the least specialization to the He/Hα escape-diagnostic problem. It has
  contributed less than Taylor to EXHALE, but not nothing: the opacity-model dispatcher came
  from it, and its C/N/O cooling fits were EXHALE's default until CHIANTI replaced them.
- **Taylor (2025/2026)** was and remains the closest match to EXHALE's goals and its richest
  source of upgrades. Four of the six adopted items come from Taylor or from the Koskinen
  code beneath it: diffusive He/H separation, molecular chemistry with H3+ cooling, the
  whole-atmosphere lower boundary, and Navier-Stokes viscosity and conduction. The one
  Taylor item that did not survive is the Penning rate, and it failed on its own internal
  inconsistency rather than on anything about EXHALE.
- **Xing (2023)** still provides the physically purest treatment of He/H mass fractionation,
  and is still not something to import. It is a validation reference.
- **Frelikh (2026)** is the closest peer to EXHALE's molecular layer and, unlike the other
  three, it is a code EXHALE can be checked *against* rather than borrowed from: twelve of
  its 22 reactions are identical to EXHALE's in value and source, and where the two differ
  EXHALE is ahead on the network (it has the He↔H charge exchange Frelikh lacks, and the
  endothermic barrier on H3+ + H that Frelikh's constant rate is missing), ahead on H3+
  cooling (non-LTE departure factor, four temperature segments, an optional exchange with
  the base field), and behind on exactly two things: the photoelectron partition in a
  molecular gas, and a composition boundary condition at the base.
- **Two of the findings belong to no single code.** The thermal-dissociation pair both codes
  take from Baulch (1992) and Ham (1970) is off thermodynamic equilibrium by a factor 10-15
  at 1000 K in the sense of over-dissociating, in both codes (section 7); and the
  extrapolated-velocity inner boundary both codes use stops being a boundary condition if
  the inflow ever goes supersonic, a weakness only EXHALE has so far found in itself.
- **The two gaps this memo was written to name are closed.** EXHALE has radius-dependent
  He/H separation, and it has a lower boundary that a photochemical model sets rather than
  the user. What is left of the original list is one data swap (the resolved He(2³S) cross
  section) and two architectural items that were correctly declined. The Frelikh comparison
  has since opened four more, of which one (the missing H2 channels in the photoelectron
  partition) is a genuine gap rather than a declined option.
