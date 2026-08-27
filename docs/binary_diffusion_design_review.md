# Review of the Binary H/He Diffusion Design

> [2026-08-27] This is the external review as written on 2026-08-25, against
> the design draft and the code of that date. It is kept unchanged as the
> record. What it describes has since moved: the design was revised against
> every finding (disposition table in `binary_diffusion_design.md` section 9),
> the trace kernel `src/modules/functions/species_diffusion.f90` was deleted
> and replaced by `src/modules/functions/binary_element_diffusion.f90`, the
> `He_Kzz` default is now `0` (`src/modules/init/parameters.f90`), and the
> `input_read.f90` rejection of molecular chemistry together with
> `He_diffusion` is gone.

## Overall assessment

`binary_diffusion_design.md` has the right high-level direction, but it is
not ready for implementation. Replacing the trace-He variable with a bounded
mass fraction, constructing equal-and-opposite binary mass fluxes, and testing
both dilute limits are sound choices. The source-code survey is also largely
accurate. However, the proposed flux has a sign error, the discrete advection
step is not coupled consistently to the hydrodynamic continuity update, and
the definitions used for metals and molecular carriers are not internally
consistent.

The settling sign, hydrodynamic coupling, mixture definition, and ambipolar
field treatment should be resolved before implementation begins.

## Findings

### 1. The settling-flux sign is reversed

Equations (2) and (3) in the design give

```text
J = -A dX/dr - B,
```

where the stated definition of `B` is positive for outward-increasing radius
and inward gravitational settling. Section 3 instead writes

```text
J_f = -(A_f + E_f) dX/dr + B_f.
```

With the given definition of `B_f`, this sends gravitational settling in the
wrong direction. Either `B_f` must include the minus sign or the face flux
must be written as

```text
J_f = -(A_f + E_f) dX/dr - B_f.
```

This is a blocking error because it reverses the physical effect that the new
operator is intended to calculate.

### 2. The proposed operator split does not match the current continuity path

The design proposes advecting `rho X` with the hydrodynamic face mass flux and
uses preservation of uniform `X` as an acceptance criterion. That property
requires the scalar update to use the same old and new density states, Runge-
Kutta stages, and face mass fluxes as the hydrodynamic continuity update.

The current call occurs after the hydrodynamic update has already produced a
new `rho`. The diffusion routine receives only the new cell-centered `rho`,
`v`, `T`, and `dt`; it does not receive the face mass flux or the old density.
Independently applying another conservative advection step to `rho X` at this
point does not, in general, preserve uniform `X` in a compressing or expanding
flow.

The design must choose and specify one of the following structures:

- Evolve `rho X` as an additional hydrodynamic conserved variable through the
  same Runge-Kutta stages.
- Pass the exact hydrodynamic face mass flux and old/new density states to a
  consistent scalar update.
- Restrict the first implementation to a steady-state formulation and defer
  physically correct transient advection explicitly.

Test T6 is useful, but the present interface cannot guarantee that it passes.

### 3. The mass-fraction definition is inconsistent when metals contribute to
`rho`

The derivation defines

```text
X = rho_He / rho
```

and treats `1-X` as the hydrogen mass fraction. Later, the write-back step uses
`1-X-X_met` for hydrogen. If metals contribute to `rho`, the following
identities used by the binary derivation no longer hold:

```text
rho_H + rho_He = rho
x = (X/m_He) / (X/m_He + (1-X)/m_H)
```

Consequently, adding only the H and He equations does not reproduce the total
continuity equation, and `J_H = -J_He` alone does not close the diffusive mass
flux of all material.

The design should either define `X` relative to the H+He subdensity,

```text
X_He,HHe = rho_He / (rho_H + rho_He),
```

or define a physically explicit two-component approximation in which trace
metals belong to one component and their flux is included consistently.
Conservation tests should cover both settings of `eos_include_metals`.

### 4. The endpoint-flux and positivity argument is incorrect

The design correctly derives a finite gradient coefficient,

```text
A = [rho_H rho_He / (rho x(1-x))] D_12 dx/dX.
```

This coefficient is generally nonzero at `X=0` and `X=1`. Therefore, a
concentration-gradient flux does not vanish merely because one cell is at an
endpoint. This is physically necessary: helium must be able to diffuse into a
helium-free cell from a neighboring cell.

The settling coefficient vanishes in the appropriate dilute limit, but the
gradient coefficient does not. The statement that both endpoint states are
invariant because all fluxes vanish should be removed.

Positivity and the upper bound must instead be established from the actual
face interpolation, nonlinear coefficient treatment, boundary values, and
tridiagonal matrix signs. The design should provide a discrete maximum-
principle or M-matrix argument for the complete operator, not only its central
gradient part.

### 5. The hydrogen-plasma ambipolar approximation is not valid for the target
He-rich limit

The retained approximation

```text
Delta m_eff/m_H = 3 - 0.5 (Zbar_He - Zbar_H)
```

assumes an electric field characteristic of a hydrogen-dominated plasma. The
target calculations extend to He/H of `1e3` to `1e4`, where the electron
pressure gradient and ambipolar field are controlled primarily by helium.
The fully ionized value therefore cannot be required to approach `2.5`
independently of composition, as T9 currently states.

Decision D3 should be resolved before accepting the H-trace limit. Suitable
choices are:

- Calculate the electric field from the local electron-pressure gradient and
  apply the corresponding charge-dependent force.
- Disable the ambipolar correction for this phase and state a neutral-only
  validity range at the coefficient calculation.

The test suite should include separate analytic limits for H-dominated and
He-dominated ionized plasmas.

### 6. Option A is conceptually preferable, but its molecular closure is not
derived

Transporting element totals continuously through the molecular region is the
better architectural choice. A single local hydrogen-carrier mass inserted
into the He-H formula, however, does not establish a binary Maxwell-Stefan
closure for a mixture containing H, H2, H2+, H3+, and HeH+.

In particular:

- The proposed `x` counts hydrogen nuclei, while the collision problem counts
  gas particles and distinct carriers.
- He-H and He-H2 collision coefficients cannot generally be combined by
  replacing the hydrogen mass with one mean mass.
- The settling and electric forces depend on the actual carrier species.
- HeH+ carries both an H nucleus and a He nucleus.

Option A should be retained as the intended architecture, but the first
implementation should either be restricted to the atomic region or use a
derived effective element mobility obtained from carrier-specific friction.
The claimed approximately 20 percent coefficient change is not sufficient to
validate the proposed closure.

### 7. The proposed HeH+ write-back can violate nonnegativity

The design proposes scaling HeH+ with the helium factor and compensating its
hydrogen nucleus in HI. This may require removing more HI than is available,
especially in a strongly molecular cell. A subsequent `ioniz_eq` call does
not justify passing through a negative or element-inconsistent intermediate
state.

The reconstruction should satisfy both element totals and nonnegativity by
construction. Two possible approaches are:

- Perform a constrained projection of the current species vector onto the
  new H and He totals.
- Expose an ionization-equilibrium interface that accepts the new element
  totals and reconstructs a valid species state directly.

### 8. Several acceptance criteria conflict with the proposed model

#### T1: equilibrium and conservation

T1 combines a Dirichlet reservoir base with conservation of total helium to
round-off. A Dirichlet boundary may supply or remove helium. The test should
use zero flux at both ends, or it should include the integrated boundary flux
in the mass budget.

#### T2: trace-limit regression

The memo states that the current cap is active over a substantial part of the
reference profile, then requires the uncapped finite-composition result to
match that profile within one percent. In addition, He/H = 0.0833 is not an
asymptotically small trace abundance.

The old calculation should remain a comparison diagnostic, not the physical
acceptance gate. A better gate is a sequence of decreasing He/H values for
which the new binary solution converges to the old trace equation at the
expected order. Observable differences at the tutorial abundance should be
reported rather than required to stay below an arbitrary threshold.

#### T7: constant or radius-dependent Kzz

T7 refers to a `K_zz` profile, while decision D5 says that only a constant
value will be implemented. The test should be described as a series of
constant-`K_zz` calculations unless profile input is added to the scope.

#### Top-boundary inflow

Section 4 specifies outflow-only advection. Section 7.2 says that negative top
velocity imports the outer-ghost composition. These are different boundary
conditions. The design must choose either no inflow or a prescribed external
composition and define the associated elemental mass budget.

## Source survey assessment

The implementation survey in the memo is largely accurate:

- `species_diffusion.f90` evolves helium number density against lagged
  hydrogen, imposes a Dirichlet base, and applies the global He/H cap.
- `input_read.f90` rejects molecular chemistry combined with He diffusion.
- The direct PTC/JFNK route does not run the diffusion co-convergence loop.
- `load_IC.f90` restores the global He/H ratio when the loaded composition
  differs and rejects the corresponding HeH+ rescaling case.
- The current drift metric omits HeTR and molecular carriers.
- The present `He_Kzz` default is `1e9 cm^2/s`, while the proposed design
  changes it to zero.

These findings provide a useful integration map. Moving the outer diffusion
iteration into one routine shared by both steady-solver entry paths and using
one metadata-driven element-count definition are appropriate changes.

## Recommended decision set

- **D1:** Retain Option A as the architecture, but do not accept the proposed
  mean-carrier approximation without a molecular transport derivation.
- **D2:** Retain outer steady-state co-convergence for this phase, while
  separately specifying how transient scalar advection shares the hydro mass
  flux.
- **D3:** Do not retain the hydrogen-plasma ambipolar approximation as the
  accepted He-rich model. Calculate the field consistently or limit the phase
  explicitly to neutral diffusion.
- **D4:** Leaving metal momentum coupling out of scope is reasonable, but the
  mass-fraction algebra must still account for metal mass consistently.
- **D5:** Changing the `He_Kzz` default to zero is reasonable if the setup
  report records the resolved value and existing examples state their value
  explicitly.
- **D6:** Include the direct-steady route and restart corrections; they are
  required for the intended workflow.
- **D7:** Capture the old baseline before changing the code, but treat it as a
  regression comparison rather than a physical truth standard.

## Conclusion

The design is a strong preliminary review of the existing implementation and
selects the correct broad state variable and flux-closure direction. It is not
yet an implementation-ready plan. The minimum required revision is to correct
the settling sign, define a discretely consistent connection to the hydro
mass flux, make the H/He state definition compatible with metal mass, and
replace or explicitly restrict the H-plasma ambipolar approximation. The
molecular closure and HeH+ reconstruction then need a nonnegative,
carrier-consistent formulation before the molecular exclusion can safely be
removed.

This review was based on a read-only comparison of the design memo with the
current implementation. No numerical tests were run. The cited published
diffusion equations were not independently verified because the ADS search
did not identify the intended paper unambiguously.
