# A2: the ROE star-state interface and the vacuum branch (approved 2026-09-05)

Decision 3 of `docs/development_plan_20260905_rev3.md` (sections 4.1 and
10.5) asks for the concrete signature before coding. This is it; the user
approved it on 2026-09-05 and the implementation follows it as written. The change
is confined to `src/modules/flux/speed_estimate_ROE.f90`, the ROE branch of
`Num_flux` in `src/modules/flux/Num_Fluxes.f90` (the estimator's only
production caller, line 170), and one refusal in `input_read.f90`.

## 1. What is wrong today (verified, `docs/development_plan_20260905_rev3.md` sections 1 and 10)

- Two-rarefaction estimate: `sqrt(...)` where the exponent is `1/z`,
  `z = (gam-1)/(2 gam)`.
- Two-shock estimate: the compression ratios are assigned to `rhoL_star`,
  `rhoR_star` without the upstream densities.
- Branch selection: the PVRS guess is kept whenever `p_max/p_min <= 2`, with
  no admissibility test, so the equal-pressure expansion `(1,-1,1)/(1,1,1)`
  returns `p_* < 0` and NaN "sound speeds".
- The estimator's outputs feed only the entropy fix (`l1R = v_star - aL_star`,
  `l3L = v_star + aR_star`); NaN there leaves the fix inactive, and the Roe
  flux itself is not positivity preserving for a vacuum-producing expansion
  (`(1,-10,1)/(1,10,1)`: one Euler update leaves negative thermal energy).

## 2. Proposed interface

```fortran
module S_estimate_ROE
   integer, parameter :: ROE_STAR_OK           = 0   ! admissible star state
   integer, parameter :: ROE_STAR_VACUUM       = 1   ! separating flow, vacuum between the gas edges
   integer, parameter :: ROE_STAR_INADMISSIBLE = 2   ! no branch gave p_* > 0 and rho_* > 0

   type roe_star_state
      real*8 :: p_star                 ! star pressure (0 in the vacuum case)
      real*8 :: rho_L_star, rho_R_star ! star densities (0 in the vacuum case)
      real*8 :: u_star                 ! contact velocity (undefined in the vacuum case, set 0)
      real*8 :: c_L_star, c_R_star     ! star sound speeds (0 in the vacuum case)
      real*8 :: v_edge_L, v_edge_R     ! gas-edge speeds v_L + 2 c_L/(gam-1), v_R - 2 c_R/(gam-1);
                                       ! filled in every case, meaningful in the vacuum case
      integer :: branch                ! 1 PVRS, 2 two-rarefaction, 3 two-shock, 0 none
   end type

   subroutine speed_estimate_ROE(WL, WR, gam, star, status)
      real*8,               intent(in)  :: WL(3), WR(3), gam
      type(roe_star_state), intent(out) :: star
      integer,              intent(out) :: status
```

Branch logic, in order:

1. Vacuum test first: if `v_R - v_L >= 2 (c_L + c_R)/(gam - 1)`, return
   `ROE_STAR_VACUUM` with the edge speeds and zero star quantities.
2. PVRS guess (Toro 9.20). Keep it only if `p_min <= p_* <= p_max`, `Q <= 2`
   and `p_* > 0`, `rho_*_L,R > 0`.
3. Otherwise two-rarefaction (Toro 9.32, exponent `1/z`) when `p_* < p_min`,
   two-shock (Toro 9.42, densities times the upstream densities) when
   `p_* > p_max`; each checked for `p_* > 0` and `rho_* > 0`.
4. If no branch is admissible, return `ROE_STAR_INADMISSIBLE` (with the
   edge speeds filled).

Sound speeds are formed only from admissible states, so no `sqrt` of a
negative ratio can occur.

## 3. The caller (`Num_flux`, ROE branch)

- `ROE_STAR_OK`: the entropy fix as today, using `star%u_star`,
  `star%c_L_star`, `star%c_R_star`.
- `ROE_STAR_VACUUM` or `ROE_STAR_INADMISSIBLE`: the Roe flux is not used for
  this face. The face flux is the HLLE flux (Einfeldt 1988; positively
  conservative, Einfeldt et al. 1991) with the wave bounds
  `S_L = min(v_L - c_L, v_avg - a_avg)`, `S_R = max(v_R + c_R, v_avg + a_avg)`
  from the Roe averages already computed in the branch:
  `NF = (S_R F_L - S_L F_R + S_L S_R (U_R - U_L)) / (S_R - S_L)` for
  `S_L < 0 < S_R`, `F_L` for `S_L >= 0`, `F_R` for `S_R <= 0`.
  `p_out = 0.5 (p_L + p_R)` as today. A counter of faces that took the
  HLLE branch is printed in the run summary (zero included).
- Constant-gamma validity stated in both headers; `input_read` refuses
  `Numerical flux: ROE` when the caloric EOS is active (message naming this
  document), because the Roe averages and the estimator are derived for one
  ideal gas.

## 4. Gates (W1's `src/tests/physics_probe/riemann_wave_speeds.f90`, updated to the new signature when A2 lands)

- Normalization invariance of the star sound speeds (rho, p scaled by 10 at
  fixed v and p/rho).
- Equal-pressure expansion `(1,-1,1)/(1,1,1)`: `ROE_STAR_OK`, `p_* > 0`,
  `c_*` within 5% of 0.957661.
- Vacuum-producing expansion `(1,-10,1)/(1,10,1)`: `ROE_STAR_VACUUM`, edge
  speeds `-6.127017`, `+6.127017`; the final `Num_flux` and one conservative
  Euler update of the left cell at dt/dx = 0.01 leave positive density and
  thermal energy.
- Stationary contact, asymmetric shock and rarefaction against exact
  Riemann solutions at the estimator's accuracy (5%).
- Every regression case byte-identical against the scratch baseline (no case
  selects ROE; `positivity_limited_fluxes` uses LLF, not ROE).

## 5. What changes for a user

`Numerical flux: ROE` behaves as before on admissible faces, refuses the
caloric EOS, and no longer produces NaN or negative energy on separating
faces. No input key is added. The signature of `speed_estimate_ROE` changes
(one caller in the code; the retained audit probes in `docs/audit_20260905/`
keep calling the old signature and are frozen as records, not rebuilt).
