module energy_semi_implicit
   ! Implicit thermal-energy source update at fixed composition and heating.
   ! In the coupled temperature-composition step, T_start gives the source
   ! evaluation temperature and u_th_old retains energy at the old
   ! composition; rebuilding it at the new particle count would change
   ! conserved energy without a source.
   !
   ! For each cell, E(T)-E_old = dt [heat-cool(T)], with caloric_eos energy:
   ! (3/2)(n_tot+n_e) k T plus bound-H2 rovibrational energy, or the
   ! monatomic T/(gamma_ad-1) branch. In energy per particle this is
   ! R(T)=u(T)-u(T_old)-a[heat-cool(T)], a=dt/(n_tot+n_e).
   !
   ! A bracketed Newton-bisection solve uses exact caloric heat capacity
   ! and a secant cooling derivative. The bracket begins at the equation-
   ! of-state floor and extends above the heating-only update. Acceptance
   ! checks the returned state for finite T, bracket membership, and
   ! |R| <= energy_res_tol N(T), with the cooling evaluated at that T and N
   ! the size of the step's own energy changes plus the arithmetic
   ! resolution of R (thermal_balance_norm states it). The stored thermal
   ! energy is not in N: a tolerance on it admits the whole net change of
   ! a nearly balanced step (scratch_thermal_acceptance_20261113/
   ! DESIGN.txt). Missing bracket, floor contact, or nonconvergence fails;
   ! the caller restores the attempted step and halves dt. No failed
   ! update modifies p, T, cool, or u. The floor is not a physical heat
   ! reservoir.
   use global_parameters
   use species_table, only: n_mion, mion_fsp,                        &
                            isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                            isp_HeIII, isp_HeTR,                      &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                            isp_OH, isp_H2O, isp_CO
   use utils
   use utils_ion_eq
   use caloric_eos, only: caloric_mixture_active,                      &
                          internal_energy_per_particle,                &
                          internal_energy_rounding_per_particle,       &
                          heat_capacity_per_particle,                  &
                          temperature_from_energy_per_particle,        &
                          energy_density_from_pressure

   implicit none

   ! Status of the temperature solve of one cell, and of the update as a
   ! whole (the worst cell). Reported through the optional `status` argument
   ! of solve_energy_semi_implicit and through energy_update_last_status, so
   ! the physical-step context and the B3a controller can read a verdict
   ! rather than infer one from the returned numbers.
   integer, parameter :: ENERGY_UPDATE_OK         = 0
   integer, parameter :: ENERGY_UPDATE_ITER_CAP   = 1
   integer, parameter :: ENERGY_UPDATE_NO_BRACKET = 2
   integer, parameter :: ENERGY_UPDATE_FLOOR      = 3
   integer, parameter :: ENERGY_UPDATE_NONFINITE  = 4

   ! Reasons, one per way a cell can fail. Kept separate from the status
   ! because two reasons share the ITER_CAP status (the budget ran out, or
   ! the bracket collapsed with the residual still above tolerance) and they
   ! say different things about the cooling function.
   integer, parameter :: ENERGY_REASON_NONE              = 0
   integer, parameter :: ENERGY_REASON_BUDGET            = 1
   integer, parameter :: ENERGY_REASON_BRACKET_COLLAPSED = 2
   integer, parameter :: ENERGY_REASON_CEILING           = 3
   integer, parameter :: ENERGY_REASON_FLOOR_REACHED     = 4
   integer, parameter :: ENERGY_REASON_BELOW_FLOOR_ON_ENTRY = 5
   integer, parameter :: ENERGY_REASON_NONFINITE_INPUT   = 6
   integer, parameter :: ENERGY_REASON_NONFINITE_COOL    = 7
   ! The cooling was finite but the caloric energy, the residual or the norm
   ! of the balance evaluated to a non-finite number, or the norm to zero:
   ! the acceptance test cannot be applied, so the cell is refused (never
   ! accepted against an infinite or zero denominator).
   integer, parameter :: ENERGY_REASON_NONFINITE_BALANCE = 8
   integer, parameter :: ENERGY_REASON_ABOVE_CEILING_ON_ENTRY = 9

   ! Acceptance of the returned state.
   !   energy_res_tol   the returned temperature T is accepted when R(T)
   !                    and N(T) are finite, N(T) > 0 and
   !                        |R(T)| <= energy_res_tol N(T) ,
   !                    R evaluated with the cooling at T, and N the norm of
   !                    thermal_balance_norm:
   !                        N = |u(T_old) - u_old| + a (|heat| + |cool(T)|)
   !                          + E_arith(T)/energy_res_tol .
   !                    The first two terms are the caloric offset and the
   !                    gross radiative exchange (thermal_balance_norm states
   !                    what the offset is); the
   !                    third is a bound of the arithmetic resolution of R
   !                    under the assumptions stated at
   !                    thermal_balance_norm.
   !                    WHAT IS GUARANTEED: the defect of the discrete
   !                    thermal balance of the returned state, as evaluated,
   !                    is at most 1e-9 of the energy the step exchanges
   !                    (gross) plus the arithmetic allowance. The gross sum
   !                    is a declared approximation of the size of the step:
   !                    when heating and cooling cancel to a fraction f of
   !                    their sum the admitted defect is 1e-9/f of the net
   !                    change, so a strongly canceled net increment is NOT
   !                    bounded relatively (it can be dropped entirely when
   !                    f < 1e-9). The temperature error that follows,
   !                    energy_res_tol N/|R'|, is a linearized estimate; it
   !                    is a bound only with a positive lower bound of |R'|
   !                    between the returned T and the root, e.g. the minimum
   !                    of c_v on that interval when cooling does not fall
   !                    with T (DESIGN.txt section 2).
   !   energy_max_iter  iteration budget, in cooling evaluations of the
   !                    grid: a TESTED OPERATIONAL CAP, not a derived worst
   !                    case. One-sided search at most quadruples a trial
   !                    but need not, and two slow Newton steps may precede
   !                    each forced bisection, so no valid bound of the
   !                    number of evaluations to adjacent doubles is stated.
   !                    Measured (src/tests/thermal_acceptance): a cooling
   !                    jump reaches BRACKET_COLLAPSED in 47-48 evaluations,
   !                    the slow-progress rows of row (j) are listed there.
   !                    A cell still active at the cap is refused with the
   !                    named reason BUDGET; it is never returned as a state.
   real*8,  parameter :: energy_res_tol     = 1.0d-9
   integer, parameter :: energy_max_iter    = 100

   ! Ends of the physical bracket, in K.
   !   The lower end. Two constraints meet here: the rovibrational table of
   !   the caloric equation of state starts at 1 K, below which its energy
   !   would be an extrapolation, and the cooling and rate fits are not
   !   defined at arbitrarily low temperature either. 0.01 T0 is at or above
   !   the first of these for every base temperature T0 >= 100 K and is the
   !   threshold this code has always used, so it is kept as the floor rather
   !   than replaced by a new number; the larger of the two is taken.
   real*8, parameter :: T_eos_floor_K    = 1.0d0
   real*8, parameter :: T_floor_code_min = 0.01d0
   !   The upper end, ONE CEILING FOR EVERY CELL. It has only to stand far
   !   above anything a photoionized wind reaches, so that "no bracket"
   !   reports a heating rate no physical cooling can balance rather than a
   !   bracket set too tight.
   !
   !   VALIDITY WHERE H2 IS PRESENT. The caloric equation of state does not
   !   end at the 5e4 K top of the H2 rovibrational table: above it
   !   h2_rovibrational_energy_and_heat_capacity continues linearly in T with
   !   the end-point heat capacity, so (u_rv, c_rv) stay a consistent pair
   !   (du_rv/dT = c_rv exactly) and u_rv stays monotone, which is the
   !   property the energy -> temperature inverse relies on; the continuation
   !   carries no dissociation energy, so extending it double counts none.
   !   What it IS above the table is an extrapolation of the bound
   !   rovibrational ladder into a region where H2 does not survive, and its
   !   weight in the energy of a cell is bounded: the mean energy of a bound
   !   ladder cannot exceed the dissociation energy, D0/k = 5.20e4 K, so
   !   against the (3/2) T of the translational term the molecular term is at
   !   most x2 D0/k/(1.5 T) of the internal energy, with x2 = n(H2)/(n_tot +
   !   n_e). At 7e4 K and x2 = 4e-4, the LHS 1140 b cells that used to be
   !   refused here, that bound is 1.9e-4. A ceiling selected by the boolean
   !   n(H2) > 0 refused the whole step for a term of that size; what has to
   !   remove H2 from a cell the wind has heated past the table is the
   !   chemistry, not the temperature bracket.
   real*8, parameter :: T_ceiling_K = 1.0d7

   ! Verdict of the last call, for a caller that reads it after the fact
   ! rather than through the optional `status` argument.
   integer, save :: energy_update_last_status   = ENERGY_UPDATE_OK
   integer, save :: energy_update_last_reason   = ENERGY_REASON_NONE
   integer, save :: energy_update_last_cell     = 0
   integer, save :: energy_update_last_step     = -1
   integer, save :: energy_update_last_nfail    = 0
   integer, save :: energy_update_last_iters    = 0
   real*8,  save :: energy_update_last_residual = 0.0d0

   ! Activations of the temperature floor of the energy update. A cell that
   ! reaches the floor has NOT converged to a physical temperature: the
   ! balance asks for a temperature below the range where the equation of
   ! state and the chemical network are defined, and at 0.01 T0 almost any
   ! composition satisfies its reaction balance
   ! (docs/Update_EXHALE_stage1.pdf section 113). These now count floor
   ! FAILURES, that is attempts: no accepted state sits on the floor, because
   ! the clamp that used to produce one is gone.
   !   hits        total activations over the run
   !   first/last  first and last step on which the floor was reached
   !   cell_hits   activations of each cell, so the number of DISTINCT cells
   !               that ever touched the floor can be reported
   integer, save :: n_energy_floor_hits = 0
   ! The same activations split by ledger family: index
   ! ledger_family_init counts the activations taken while the run was
   ! reaching a state, index ledger_family_phys those taken inside accepted
   ! physical steps. Two floor activations belonging to different run states
   ! are not one budget, and the total above is their sum.
   integer, save :: n_energy_floor_hits_family(2) = 0
   integer, save :: energy_floor_first_step = -1
   integer, save :: energy_floor_last_step  = -1
   integer, allocatable, save :: energy_floor_cell_hits(:)

   ! Working state of the bracketed solve, one entry per cell. It is a type
   ! rather than a set of locals because the cooling function is evaluated
   ! for the whole grid at once (eval_cool solves the fine-structure line
   ! transfer over the columns before its cell-local part), so the iteration
   ! has to advance every cell together and hand the caller one vector of
   ! temperatures to evaluate the cooling at. The caller drives the loop:
   !   call energy_balance_init(s, ...)
   !   do while (s%n_active > 0)          ! and within the budget
   !      <evaluate cooling at s%T_eval>
   !      call energy_balance_update(s, cool_at_T_eval)
   !   enddo
   !   call energy_balance_finish(s)
   ! which is what lets the test drivers substitute an analytic cooling for
   ! eval_cool without touching the solver.
   type :: energy_balance_state
      ! The problem
      real*8, allocatable :: a(:)          ! dt/(n_tot + n_e)
      real*8, allocatable :: heat(:)       ! frozen heating of the step
      real*8, allocatable :: T_old(:)
      real*8, allocatable :: u_old(:)      ! the energy the step starts from,
                                           ! per particle
      ! u(T_old) - u_old at the composition of the solve, the first term
      ! of the norm: the offset is the caloric energy at this solve's starting
      ! temperature and composition relative to the source-step energy
      ! anchor. On the first sweep it measures composition-induced change
      ! at fixed temperature, up to reconstruction/rounding; on later
      ! sweeps it can also contain the thermal change accumulated by the
      ! coupled iteration.
      real*8, allocatable :: du_caloric(:)
      real*8, allocatable :: T_floor(:), T_ceil(:)
      ! The iteration
      real*8, allocatable :: T_eval(:)     ! where the caller must evaluate cool
      real*8, allocatable :: lo(:), hi(:)
      logical, allocatable :: has_lo(:), has_hi(:)
      real*8, allocatable :: T_prev(:), cool_prev(:), dcool_dT(:)
      real*8, allocatable :: width_prev(:)
      integer, allocatable :: n_slow(:)
      logical, allocatable :: active(:)
      ! The verdict
      integer, allocatable :: status(:), reason(:)
      real*8, allocatable :: T_new(:), cool_new(:)
      real*8, allocatable :: res_scaled(:)   ! |R|/N at the returned (or
                                             ! last) T, N of
                                             ! thermal_balance_norm
      real*8, allocatable :: res_energy(:)   ! R itself, per particle
      ! N(T) of thermal_balance_norm at the returned (or last) temperature,
      ! the denominator of res_scaled. Not a maximum over iterates: the
      ! acceptance is a property of the returned state, not of the path.
      real*8, allocatable :: scale(:)
      integer :: n_active = 0
      integer :: iters    = 0
   end type energy_balance_state

contains


   logical function energy_is_finite(x)
   ! True for a normal or subnormal double, false for a NaN or an infinity.
   real*8, intent(in) :: x
   energy_is_finite = (x .eq. x) .and. (abs(x) .le. huge(1.0d0))
   end function energy_is_finite


   real*8 function thermal_balance_norm(s, j, T, u_T, C) result(nrm)
   ! N(T) of the acceptance test |R(T)| <= energy_res_tol N(T) of cell j,
   ! for R(T) = u(T) - u_old - a (heat - C) with u_T = u(T) as the residual
   ! evaluated it and C the cooling at T. The same norm serves the test at
   ! T_old (C = cool0), the test of every later trial, and the stored
   ! s%scale.
   !
   !   N = S + E_arith/energy_res_tol ,
   !   S = |u(T_old) - u_old| + a (|heat| + |C|) .
   ! S is the size of the step. Its first term is the caloric offset:
   ! the offset is the caloric energy at this solve's starting
   ! temperature and composition relative to the source-step energy
   ! anchor. On the first sweep it measures composition-induced change
   ! at fixed temperature, up to reconstruction/rounding; on later
   ! sweeps it can also contain the thermal change accumulated by the
   ! coupled iteration. It measures the size of the step and is not
   ! a source in the balance, and no formation energy is added to it.
   !
   ! E_arith bounds what the arithmetic cannot resolve in R as evaluated
   ! (eps = epsilon(1d0) = 2^-52; one rounding of r costs (eps/2)|r|):
   !   err_u(T)                 the rounding of u(T) itself
   !                            (internal_energy_rounding_per_particle);
   !   eps max(|u(T)|, |u_old|) the stored-energy difference u(T) - u_old:
   !                            its subtraction and one unit in the last
   !                            place of the larger operand, the finest
   !                            energy difference two stored energies of
   !                            this size can express;
   !   eps a (|heat| + |C|)     the difference heat - C, the product with a
   !                            and the final sum;
   !   eps c_v(T) |T|           the change of u over one unit in the last
   !                            place of T (spacing <= eps T), the closest
   !                            a representable temperature can come to the
   !                            root.
   ! VALIDITY ASSUMPTIONS (not derived here; E_arith is a bound of the
   ! resolution of R only under them):
   !   (A1) the logarithmic slope of the cooling between neighboring
   !        doubles of T is bounded, n = |dln C/dln T| <= energy_res_tol/eps
   !        = 4.5e6, so that the change of a C over one unit in the last
   !        place of T, eps a |C| n, is inside energy_res_tol a |C|. A
   !        cooling table, an interpolation boundary, an escape factor or a
   !        signed sum of terms can violate it; then a continuous R can
   !        still fail at both ends of an adjacent bracket.
   !   (A2) heat and C as handed in are evaluated with a relative error
   !        below energy_res_tol (the error of eval_cool's own arithmetic
   !        and of the cool0 the caller supplies is not bounded here).
   !   (A3) log and the other intrinsics inside u(T) are accurate to one
   !        unit in the last place, and the arithmetic is ordinary IEEE
   !        double with round to nearest (fused multiply-add only removes
   !        roundings).
   ! No constant floor is added: E_arith > 0 for every finite T > 0 with
   ! c_v > 0. The caller checks that the returned N is finite and positive.
   ! OVERFLOW. The separately evaluated exchange products avoid the
   ! overflow of the unscaled heat-plus-cooling sum. A nonfinite term or
   ! accumulated norm causes refusal (balance_is_testable, at the caller).
   ! The two products are VOLATILE: the production flags (-O3, the
   ! compiler's default floating-point model) allow reassociation, and
   ! without it the compiler refactored a|heat| + a|C| back into
   ! a(|heat| + |C|), which overflowed in the unit test of row (i)
   ! (src/tests/thermal_acceptance). A volatile value is stored and
   ! reloaded, so the products cannot be recombined across it.
   type(energy_balance_state), intent(in) :: s
   integer, intent(in) :: j
   real*8,  intent(in) :: T, u_T, C
   real*8  :: eps_over_tol
   real*8, volatile :: heat_exchange, cool_exchange

   eps_over_tol  = epsilon(1.0d0)/energy_res_tol
   heat_exchange = s%a(j)*abs(s%heat(j))
   cool_exchange = s%a(j)*abs(C)
   nrm = abs(s%du_caloric(j)) + heat_exchange + cool_exchange           &
       + internal_energy_rounding_per_particle(j, T)/energy_res_tol     &
       + eps_over_tol*max(abs(u_T), abs(s%u_old(j)))                    &
       + eps_over_tol*heat_exchange + eps_over_tol*cool_exchange        &
       + eps_over_tol*(heat_capacity_per_particle(j, T)*abs(T))
   end function thermal_balance_norm


   logical function balance_is_testable(u_T, residual, nrm)
   ! The acceptance test |R| <= energy_res_tol N can be applied: the
   ! caloric energy and the residual are finite and N is finite and
   ! positive. Required before EVERY acceptance (at T_old and at a trial).
   real*8, intent(in) :: u_T, residual, nrm
   balance_is_testable = energy_is_finite(u_T) .and.                     &
                         energy_is_finite(residual) .and.                &
                         energy_is_finite(nrm) .and. (nrm .gt. 0.0d0)
   end function balance_is_testable


   integer function n_energy_floor_cells()
   ! Number of DISTINCT cells that have reached the temperature floor of the
   ! energy update at least once. Written as a module function because the
   ! marching loop of EXHALE_main carries a local integer named "count",
   ! which shadows the Fortran intrinsic of the same name at that call site.
   integer :: j
   n_energy_floor_cells = 0
   if (.not. allocated(energy_floor_cell_hits)) return
   do j = lbound(energy_floor_cell_hits,1), ubound(energy_floor_cell_hits,1)
      if (energy_floor_cell_hits(j) .gt. 0)                                &
         n_energy_floor_cells = n_energy_floor_cells + 1
   enddo
   end function n_energy_floor_cells


   subroutine energy_balance_init(s, T_old, a, heat, cool0, T_floor,       &
                                  T_ceil, u_old_in)
   ! Set up the solve of R(T) = 0 for every cell and choose the first
   ! temperature at which the caller must evaluate the cooling.
   !
   ! cool0 is the cooling at T_old and at the composition of this step, which
   ! ioniz_eq has already evaluated, so the residual at T_old costs nothing
   ! and a cell whose step needs no temperature change never calls the
   ! cooling function at all.
   !
   ! u_old_in IS THE SOURCE-STEP ENERGY ANCHOR, per (n_tot + n_e) of the
   ! composition the solve is being run at.  It is a separate argument
   ! because the temperature the solve starts from and the energy anchor
   ! need not belong to the same state.  In the coupled source step the
   ! caller retains u_th_old = u(3) - rho v**2/2 once, before its sweep
   ! loop, and hands that same anchor to every sweep (no formation-energy
   ! term is subtracted on that path; the optional du_form of
   ! solve_energy_semi_implicit is not passed there), while T_old is the
   ! temperature of the current sweep.  u(T_old) - u_old is then the caloric
   ! offset: on the first sweep the composition-induced change at fixed
   ! temperature, up to reconstruction/rounding; on later sweeps it can also
   ! contain the thermal change accumulated by the coupled iteration.
   ! Absent, the routine sets u_old = internal_energy_per_particle(T_old),
   ! which is what a caller that advances the temperature at frozen
   ! composition wants (offset zero by construction).
   type(energy_balance_state), intent(out) :: s
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T_old, a, heat, cool0
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T_floor, T_ceil
   real*8, dimension(1-Ng:N+Ng), intent(in), optional :: u_old_in
   real*8  :: residual_at_T_old, sc, u_expl, T_guess, u_T_old
   logical :: root_is_colder
   integer :: j

   allocate(s%a(1-Ng:N+Ng), s%heat(1-Ng:N+Ng), s%T_old(1-Ng:N+Ng),         &
            s%u_old(1-Ng:N+Ng), s%du_caloric(1-Ng:N+Ng),                   &
            s%T_floor(1-Ng:N+Ng), s%T_ceil(1-Ng:N+Ng),                     &
            s%T_eval(1-Ng:N+Ng), s%lo(1-Ng:N+Ng), s%hi(1-Ng:N+Ng),         &
            s%has_lo(1-Ng:N+Ng), s%has_hi(1-Ng:N+Ng),                      &
            s%T_prev(1-Ng:N+Ng), s%cool_prev(1-Ng:N+Ng),                   &
            s%dcool_dT(1-Ng:N+Ng), s%width_prev(1-Ng:N+Ng),                &
            s%n_slow(1-Ng:N+Ng), s%active(1-Ng:N+Ng),                      &
            s%status(1-Ng:N+Ng), s%reason(1-Ng:N+Ng),                      &
            s%T_new(1-Ng:N+Ng), s%cool_new(1-Ng:N+Ng),                     &
            s%res_scaled(1-Ng:N+Ng), s%res_energy(1-Ng:N+Ng),              &
            s%scale(1-Ng:N+Ng))

   s%a       = a
   s%heat    = heat
   s%T_old   = T_old
   s%T_floor = T_floor
   s%T_ceil  = T_ceil
   s%has_lo  = .false.
   s%has_hi  = .false.
   s%lo      = 0.0d0
   s%hi      = 0.0d0
   s%dcool_dT = 0.0d0
   s%width_prev = huge(1.0d0)
   s%n_slow  = 0
   s%active  = .false.
   s%status  = ENERGY_UPDATE_OK
   s%reason  = ENERGY_REASON_NONE
   s%T_new   = T_old
   s%T_eval  = T_old
   s%u_old   = 0.0d0
   s%du_caloric = 0.0d0
   ! A cell refused on entry keeps scale 0 and res_scaled = huge.
   s%scale   = 0.0d0
   s%cool_new = cool0
   s%res_scaled = 0.0d0
   s%res_energy = 0.0d0
   s%iters   = 0
   s%n_active = 0

   do j = 1-Ng, N+Ng
      if (present(u_old_in)) s%u_old(j) = u_old_in(j)
      ! Admissibility of what this cell was handed.
      if (.not. energy_is_finite(T_old(j)) .or. .not. (T_old(j) .gt. 0.0d0) &
          .or. .not. energy_is_finite(a(j)) .or. .not. (a(j) .ge. 0.0d0)    &
          .or. .not. energy_is_finite(heat(j))                              &
          .or. .not. energy_is_finite(cool0(j))                          &
          .or. .not. energy_is_finite(s%u_old(j))) then
         s%status(j) = ENERGY_UPDATE_NONFINITE
         s%reason(j) = ENERGY_REASON_NONFINITE_INPUT
         s%res_scaled(j) = huge(1.0d0)
         ! The cooling is still evaluated for the whole grid while the
         ! remaining cells iterate, so this entry must hold a number.
         s%T_eval(j) = 1.0d0
         cycle
      endif

      if (.not. present(u_old_in))                                        &
         s%u_old(j) = internal_energy_per_particle(j, T_old(j))
      s%T_eval(j) = T_old(j)
      s%T_prev(j) = T_old(j)
      s%cool_prev(j) = cool0(j)

      ! The incoming temperature is already outside the range in which the
      ! equation of state is defined: nothing this routine can return is a
      ! physical state.
      if (T_old(j) .lt. T_floor(j)) then
         s%status(j) = ENERGY_UPDATE_FLOOR
         s%reason(j) = ENERGY_REASON_BELOW_FLOOR_ON_ENTRY
         s%res_scaled(j) = huge(1.0d0)
         cycle
      endif
      ! The same on the other side: a temperature above the ceiling is
      ! outside the declared range, and accepting it at T_old (zero
      ! sources) would return a state outside [T_floor, T_ceil].
      if (T_old(j) .gt. T_ceil(j)) then
         s%status(j) = ENERGY_UPDATE_NO_BRACKET
         s%reason(j) = ENERGY_REASON_ABOVE_CEILING_ON_ENTRY
         s%res_scaled(j) = huge(1.0d0)
         cycle
      endif

      ! R at T_old.  The first two terms cancel identically only when the
      ! energy the step starts from is the energy of T_old at this
      ! composition; with u_old_in they do not, and what is left is the
      ! caloric offset, also the first term of the norm (see
      ! thermal_balance_norm for what it contains on first and later
      ! sweeps).
      ! With heat, cool0 >= 0 the difference cool0 - heat cannot overflow;
      ! a non-finite residual (a genuinely huge term, or heat and cooling
      ! of opposite sign near the overflow threshold) is refused below.
      u_T_old = internal_energy_per_particle(j, T_old(j))
      s%du_caloric(j) = u_T_old - s%u_old(j)
      residual_at_T_old = s%du_caloric(j) + a(j)*(cool0(j) - heat(j))
      sc = thermal_balance_norm(s, j, T_old(j), u_T_old, cool0(j))
      s%scale(j) = sc
      if (.not. balance_is_testable(u_T_old, residual_at_T_old, sc)) then
         s%status(j) = ENERGY_UPDATE_NONFINITE
         s%reason(j) = ENERGY_REASON_NONFINITE_BALANCE
         s%res_scaled(j) = huge(1.0d0)
         s%res_energy(j) = residual_at_T_old
         cycle
      endif

      ! T_old as a trial. cool0 is the cooling at T_old and at the
      ! composition of this solve, so this is the returned-state test of
      ! energy_balance_update applied to T = T_old, with the same norm, at no
      ! cost. It is passed by an exact equilibrium (R = 0), by a cell with
      ! no source and no composition change, and by a cell whose whole
      ! balance at T_old is within energy_res_tol of the energy the step
      ! moves plus rounding; it is NOT passed because the stored thermal
      ! energy is large.
      if (abs(residual_at_T_old) .le. energy_res_tol*sc) then
         s%T_new(j)      = T_old(j)
         s%cool_new(j)   = cool0(j)
         s%res_scaled(j) = abs(residual_at_T_old)/sc
         s%res_energy(j) = residual_at_T_old
         cycle
      endif

      ! One bracket end is known from the sign at T_old: R is the energy the
      ! state has minus the energy the balance allows it, so R > 0 means the
      ! root is colder.
      root_is_colder = residual_at_T_old .gt. 0.0d0
      if (root_is_colder) then
         s%hi(j) = T_old(j);  s%has_hi(j) = .true.
      else
         s%lo(j) = T_old(j);  s%has_lo(j) = .true.
      endif

      ! First trial: the explicit-Euler point, the temperature the step would
      ! land on if the cooling did not follow T. With the cooling frozen at
      ! cool0 it is the root, so it lies on the root's side of T_old unless
      ! the inverse u -> T or the rounding of u_expl puts it on T_old or
      ! behind it. Then the Newton point from T_old with the caloric heat
      ! capacity, and if that too rounds back to T_old, the adjacent double
      ! on the root's side: never a trial far beyond a root that sits next
      ! to T_old.
      u_expl  = s%u_old(j) + a(j)*(heat(j) - cool0(j))
      T_guess = temperature_from_energy_per_particle(j, u_expl)
      if (.not. strictly_on_root_side(T_guess)) then
         T_guess = T_old(j) - residual_at_T_old                            &
                   /heat_capacity_per_particle(j, T_old(j))
         if (.not. strictly_on_root_side(T_guess)) then
            if (root_is_colder) then
               T_guess = nearest(T_old(j), -1.0d0)
            else
               T_guess = nearest(T_old(j),  1.0d0)
            endif
         endif
      endif
      T_guess = min(max(T_guess, T_floor(j)), T_ceil(j))
      s%T_eval(j) = T_guess
      s%active(j) = .true.
      s%n_active  = s%n_active + 1
      ! Clamping brought the trial back to T_old: T_old is on the floor or
      ! the ceiling and the balance asks for a temperature beyond it.
      if (.not. strictly_on_root_side(T_guess)) then
         if (root_is_colder) then
            call energy_balance_fail(s, j, ENERGY_UPDATE_FLOOR,            &
                                     ENERGY_REASON_FLOOR_REACHED,          &
                                     T_old(j), cool0(j),                   &
                                     abs(residual_at_T_old)/sc,            &
                                     residual_at_T_old)
         else
            call energy_balance_fail(s, j, ENERGY_UPDATE_NO_BRACKET,       &
                                     ENERGY_REASON_CEILING,                &
                                     T_old(j), cool0(j),                   &
                                     abs(residual_at_T_old)/sc,            &
                                     residual_at_T_old)
         endif
         ! The cooling is still evaluated for the whole grid while the
         ! remaining cells iterate, so this entry must hold a number.
         s%T_eval(j) = T_old(j)
      endif
   enddo

   contains

      logical function strictly_on_root_side(T_try)
      ! T_try is finite and on the root's side of T_old, not equal to it.
      real*8, intent(in) :: T_try
      if (.not. energy_is_finite(T_try)) then
         strictly_on_root_side = .false.
      else if (root_is_colder) then
         strictly_on_root_side = T_try .lt. T_old(j)
      else
         strictly_on_root_side = T_try .gt. T_old(j)
      endif
      end function strictly_on_root_side

   end subroutine energy_balance_init


   subroutine energy_balance_update(s, cool_at_T_eval)
   ! Absorb the cooling the caller evaluated at s%T_eval, accept or reject
   ! each active cell, and choose the next temperature for the ones that
   ! remain.
   type(energy_balance_state), intent(inout) :: s
   real*8, dimension(1-Ng:N+Ng), intent(in) :: cool_at_T_eval
   real*8  :: T_e, C, residual_at_T_e, sc, cv, dR, T_n, w, u_T_e
   integer :: j

   s%iters = s%iters + 1

   do j = 1-Ng, N+Ng
      if (.not. s%active(j)) cycle
      T_e = s%T_eval(j)
      C   = cool_at_T_eval(j)

      if (.not. energy_is_finite(C)) then
         call energy_balance_fail(s, j, ENERGY_UPDATE_NONFINITE,           &
                                  ENERGY_REASON_NONFINITE_COOL, T_e, C,    &
                                  huge(1.0d0), huge(1.0d0))
         cycle
      endif

      u_T_e = internal_energy_per_particle(j, T_e)
      residual_at_T_e = u_T_e - s%u_old(j) - s%a(j)*(s%heat(j) - C)

      ! The norm of the state that would be returned, with the cooling at
      ! the temperature that would be returned (thermal_balance_norm). The
      ! test below needs a finite residual and a finite positive norm; a
      ! failed evaluation is a named refusal, never an acceptance.
      sc = thermal_balance_norm(s, j, T_e, u_T_e, C)
      s%scale(j) = sc
      if (.not. balance_is_testable(u_T_e, residual_at_T_e, sc)) then
         call energy_balance_fail(s, j, ENERGY_UPDATE_NONFINITE,           &
                                  ENERGY_REASON_NONFINITE_BALANCE, T_e, C, &
                                  huge(1.0d0), residual_at_T_e)
         cycle
      endif

      ! Acceptance: the residual of the state that would be returned, with
      ! the cooling evaluated at the temperature that would be returned.
      if (abs(residual_at_T_e) .le. energy_res_tol*sc) then
         s%T_new(j)      = T_e
         s%cool_new(j)   = C
         s%res_scaled(j) = abs(residual_at_T_e)/sc
         s%res_energy(j) = residual_at_T_e
         s%active(j)     = .false.
         s%n_active      = s%n_active - 1
         cycle
      endif

      ! Bracket update.
      if (residual_at_T_e .lt. 0.0d0) then
         s%lo(j) = T_e;  s%has_lo(j) = .true.
      else
         s%hi(j) = T_e;  s%has_hi(j) = .true.
      endif

      ! Secant estimate of dcool/dT for the Newton step. It steers only the
      ! path; the bracket owns the convergence.
      if (T_e .ne. s%T_prev(j))                                            &
         s%dcool_dT(j) = (C - s%cool_prev(j))/(T_e - s%T_prev(j))
      s%T_prev(j)    = T_e
      s%cool_prev(j) = C

      cv = heat_capacity_per_particle(j, T_e)
      dR = cv + s%a(j)*s%dcool_dT(j)
      if (.not. energy_is_finite(dR) .or. dR .le. 0.0d0) dR = cv
      T_n = T_e - residual_at_T_e/dR

      if (s%has_lo(j) .and. s%has_hi(j)) then
         w = s%hi(j) - s%lo(j)
         if (nearest(s%lo(j), 1.0d0) .ge. s%hi(j)) then
            ! No double lies strictly between the ends, both of which were
            ! evaluated and failed the test: there is no interior
            ! representable temperature left to try. This proves no more
            ! than that. Under the assumptions (A1)-(A3) stated at
            ! thermal_balance_norm a continuous R would have passed at one
            ! end, so the failure means a jump of cool(T) across the
            ! bracket OR a violation of those assumptions (a cooling
            ! steeper than they allow, or an evaluation noisier than
            ! energy_res_tol of itself). Reported, never returned as a
            ! state.
            call energy_balance_fail(s, j, ENERGY_UPDATE_ITER_CAP,         &
                                     ENERGY_REASON_BRACKET_COLLAPSED,      &
                                     T_e, C, abs(residual_at_T_e)/sc,      &
                                     residual_at_T_e)
            cycle
         endif
         if (.not. energy_is_finite(T_n) .or. T_n .le. s%lo(j)             &
             .or. T_n .ge. s%hi(j)) then
            T_n = bracket_midpoint()
            s%n_slow(j) = 0
         else
            if (w .gt. 0.5d0*s%width_prev(j)) then
               s%n_slow(j) = s%n_slow(j) + 1
            else
               s%n_slow(j) = 0
            endif
            if (s%n_slow(j) .ge. 2) then
               T_n = bracket_midpoint()
               s%n_slow(j) = 0
            endif
         endif
         s%width_prev(j) = w
      else if (.not. s%has_hi(j)) then
         ! Still searching upward: the balance wants a hotter cell than any
         ! temperature tried so far. With R(T_e) < 0 and dR > 0 a finite Newton
         ! point at or below T_e is a step below half a unit in the last
         ! place; the next double up is then the trial, not a doubling.
         if (.not. energy_is_finite(T_n)) then
            T_n = 2.0d0*T_e + 1.0d-3
         else if (T_n .le. T_e) then
            T_n = nearest(T_e, 1.0d0)
         endif
         T_n = min(T_n, 4.0d0*T_e + 1.0d-3)
         T_n = min(T_n, s%T_ceil(j))
         if (T_n .le. T_e) then
            ! At the top of the range where the equation of state and the
            ! cooling fits are defined and still short of the balance: no
            ! bracket exists inside the physical range.
            call energy_balance_fail(s, j, ENERGY_UPDATE_NO_BRACKET,       &
                                     ENERGY_REASON_CEILING, T_e, C,        &
                                     abs(residual_at_T_e)/sc, residual_at_T_e)
            cycle
         endif
      else
         ! Still searching downward (the same rounding case mirrored).
         if (.not. energy_is_finite(T_n)) then
            T_n = 0.5d0*T_e
         else if (T_n .ge. T_e) then
            T_n = nearest(T_e, -1.0d0)
         endif
         T_n = max(T_n, 0.25d0*T_e)
         T_n = max(T_n, s%T_floor(j))
         if (T_n .ge. T_e) then
            ! The balance asks for a temperature below the floor of the
            ! equation of state. res_energy is the energy per particle the
            ! cell would have to shed beyond what the floor allows; divided
            ! by a = dt/(n_tot+n_e) it is the volumetric rate, in the units
            ! of `heat`, that a specified external reservoir would have to
            ! supply for the floor to be a physical state instead of a
            ! failure.
            call energy_balance_fail(s, j, ENERGY_UPDATE_FLOOR,            &
                                     ENERGY_REASON_FLOOR_REACHED, T_e, C,  &
                                     abs(residual_at_T_e)/sc, residual_at_T_e)
            cycle
         endif
      endif

      s%T_eval(j) = T_n
   enddo

   contains

      real*8 function bracket_midpoint()
      ! The midpoint of the bracket of cell j, kept strictly inside it: at
      ! a width of two units in the last place across a binade boundary the
      ! rounded midpoint can land on an end, and the next double above lo
      ! is then the interior point (the caller has checked that one exists).
      bracket_midpoint = 0.5d0*(s%lo(j) + s%hi(j))
      if (bracket_midpoint .le. s%lo(j) .or. bracket_midpoint .ge. s%hi(j)) &
         bracket_midpoint = nearest(s%lo(j), 1.0d0)
      end function bracket_midpoint

   end subroutine energy_balance_update


   subroutine energy_balance_fail(s, j, status_in, reason_in, T_e, C,      &
                                  res_scaled_in, res_energy_in)
   ! Record the failure of cell j and take it out of the iteration. T_new and
   ! cool_new keep the last iterate so that a diagnostic can print it; they
   ! are NOT a state, and solve_energy_semi_implicit assembles no state from
   ! a failed update.
   type(energy_balance_state), intent(inout) :: s
   integer, intent(in) :: j, status_in, reason_in
   real*8,  intent(in) :: T_e, C, res_scaled_in, res_energy_in
   s%status(j)     = status_in
   s%reason(j)     = reason_in
   s%T_new(j)      = T_e
   s%cool_new(j)   = C
   s%res_scaled(j) = res_scaled_in
   s%res_energy(j) = res_energy_in
   s%active(j)     = .false.
   s%n_active      = s%n_active - 1
   end subroutine energy_balance_fail


   subroutine energy_balance_finish(s)
   ! Cells still active when the budget ran out.
   type(energy_balance_state), intent(inout) :: s
   integer :: j
   do j = 1-Ng, N+Ng
      if (.not. s%active(j)) cycle
      s%status(j) = ENERGY_UPDATE_ITER_CAP
      s%reason(j) = ENERGY_REASON_BUDGET
      s%active(j) = .false.
      s%n_active  = s%n_active - 1
   enddo
   end subroutine energy_balance_finish


   subroutine energy_balance_reason_text(reason_in, text)
   integer, intent(in) :: reason_in
   character(len=*), intent(out) :: text
   select case (reason_in)
   case (ENERGY_REASON_BUDGET)
      text = 'iteration budget exhausted'
   case (ENERGY_REASON_BRACKET_COLLAPSED)
      text = 'bracket collapsed with the residual above tolerance'
   case (ENERGY_REASON_CEILING)
      text = 'no bracket below the equation-of-state ceiling'
   case (ENERGY_REASON_FLOOR_REACHED)
      text = 'balance asks for a temperature below the floor'
   case (ENERGY_REASON_BELOW_FLOOR_ON_ENTRY)
      text = 'incoming temperature already below the floor'
   case (ENERGY_REASON_NONFINITE_INPUT)
      text = 'non-finite temperature, step or source on entry'
   case (ENERGY_REASON_NONFINITE_COOL)
      text = 'non-finite cooling'
   case (ENERGY_REASON_NONFINITE_BALANCE)
      text = 'non-finite residual or norm, or zero norm'
   case (ENERGY_REASON_ABOVE_CEILING_ON_ENTRY)
      text = 'incoming temperature already above the ceiling'
   case default
      text = 'converged'
   end select
   end subroutine energy_balance_reason_text


   subroutine solve_energy_semi_implicit(u, W, dt, heat, cool, f_sp, step, &
                                         status, T_start, u_th_old, du_form)
      ! Input/Output variables
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: W
      ! Cell-by-cell pseudo-time steps (uniform = global dt unless
      ! "Time stepping: Local"; the cell solve below is local anyway).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dt
      ! Heating of the state this step starts from, from ioniz_eq. It is
      ! held fixed across the temperature update (the photoheating rates do
      ! not follow T within one step).
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat
      ! On entry: the cooling ioniz_eq evaluated at T_old and at the
      ! post-sweep composition -- the same composition and the same
      ! temperature this routine would evaluate it at, so it is taken as the
      ! residual at T_old instead of being recomputed.
      ! On exit: the cooling at the temperature this step lands on, evaluated
      ! there.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      ! Marching step index, for the diagnostics of a failure.
      integer, intent(in) :: step
      ! Verdict of the update as a whole: the worst cell's status. Optional
      ! so that the call site is unchanged; the same value is left in
      ! energy_update_last_status for a caller that reads it later.
      integer, intent(out), optional :: status
      ! THE INNER STEP OF THE COUPLED SOURCE STEP
      ! All three are absent when
      ! the caller advances the temperature at frozen composition, and the
      ! routine then behaves exactly as it did before B3c.
      !   T_start   the temperature the composition and the sources were
      !             evaluated at.  Absent, it is recovered from the pressure
      !             as p/(n_tot + n_e), which is the same number only when
      !             the pressure was rebuilt at the returned composition --
      !             that rebuild is the C1 projection and B3c removes it.
      !   u_th_old  the thermal energy density the cell HAD before the source
      !             step, at its OLD composition, in the code's energy-density
      !             units (the units of u(3,:)).  Absent, the routine forms it
      !             from T_start at the CURRENT composition, which is the C1
      !             projection written as an energy.
      !   du_form   the change of the formation/excitation reservoir
      !             sum_s eps_s (n_s^new - n_s^old) over the composition
      !             change, same units.  The energy row solved is then
      !               u_th(T,c) - u_th_old + du_form = dt [heat - cool(T)] ,
      !             so that u_th + u_form is conserved by the step up to the
      !             external exchange heat and cool carry (T1.5).
      real*8, dimension(1-Ng:N+Ng), intent(in), optional :: T_start
      real*8, dimension(1-Ng:N+Ng), intent(in), optional :: u_th_old
      real*8, dimension(1-Ng:N+Ng), intent(in), optional :: du_form

      ! Local arrays for dimensional calculations
      real*8, dimension(1-Ng:N+Ng) :: T_K, T_old
      real*8, dimension(1-Ng:N+Ng) :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      ! Metal ion densities (canonical species_table order)
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      ! Molecular densities: H2, H2+, H3+, HeH+ -- adimensional (nmol_l, the
      ! units calc_ne is called with here) and cgs (nmol_dim, the units
      ! eval_cool takes alongside nhi..nm)
      real*8, dimension(1-Ng:N+Ng,4) :: nmol_l, nmol_dim
      ! Oxygen carriers OH, H2O, CO -- adimensional (nox_l, the units
      ! calc_ntot is called with here) and cgs (nox_dim, the units eval_cool
      ! needs for the H2O and CO infrared bands)
      real*8, dimension(1-Ng:N+Ng,3) :: nox_l, nox_dim
      real*8, dimension(1-Ng:N+Ng) :: rchiiB, rcheiiB, rcheiiiB
      real*8, dimension(1-Ng:N+Ng) :: a_ion_HI, a_ion_HeI, a_ion_HeII
      ! Metal rates for each ion returned by eval_cool but unused here
      real*8, dimension(1-Ng:N+Ng,n_mion) :: rec_m, aion_m
      real*8, dimension(1-Ng:N+Ng) :: cool_dim, cool_iter
      real*8, dimension(1-Ng:N+Ng) :: rho, v, p
      real*8, dimension(1-Ng:N+Ng) :: zero_arr
      real*8, dimension(1-Ng:N+Ng) :: a_coef, T_floor, T_ceil
      ! Energy the step starts from, per (n_tot + n_e) of the composition it
      ! is solved at, when the caller supplies it (see T_start above).
      real*8, dimension(1-Ng:N+Ng) :: u_old_pp
      logical :: coupled

      ! Adimensional number densities
      real*8, dimension(1-Ng:N+Ng) :: ne_ad, n_tot_ad
      real*8  :: nk
      integer :: j, im, j_bad, n_fail, verdict
      type(energy_balance_state) :: s
      character(len=64) :: reason_text

      ! Zero array for the calls that take no metal ions
      zero_arr = 0.0d0

      ! Extract primitive variables
      rho = W(1,:)
      v   = W(2,:)
      p   = W(3,:)

      ! Dimensional species densities (n0 is the density normalization)
      nhi  = rho*f_sp(:,isp_HI)*n0
      nhii = rho*f_sp(:,isp_HII)*n0
      if (thereis_He) then
         nhei   = rho*f_sp(:,isp_HeI)*n0
         nheii  = rho*f_sp(:,isp_HeII)*n0
         nheiii = rho*f_sp(:,isp_HeIII)*n0
         if (thereis_HeITR) then
            nheiTR = rho*f_sp(:,isp_HeTR)*n0
         else
            nheiTR = 0.0d0
         endif
      else
         nhei   = 0.0d0
         nheii  = 0.0d0
         nheiii = 0.0d0
         nheiTR = 0.0d0
      endif
      ! Metal ion densities in canonical species_table order
      ! (mion_fsp = [7..18], so nm(:,im) reproduces each ion's
      ! rho*f_sp(:,col)*n0 expressions bit-for-bit).
      do im = 1,n_mion
         nm(:,im) = rho*f_sp(:,mion_fsp(im))*n0
      enddo
      ! Molecular densities (adimensional, matching the nm/n0 units used below).
      ! Zero for non-molecular runs (f_sp molecular columns are zero), so
      ! calc_ne/calc_ntot add exactly zero; for molecular runs they restore the
      ! neutral-H2 particle count and the molecular-ion electrons to n_tot/ne
      ! (first-order near a molecular base, where H2 dominates the particle
      ! budget).
      nmol_l(:,1) = rho*f_sp(:,isp_H2)
      nmol_l(:,2) = rho*f_sp(:,isp_H2p)
      nmol_l(:,3) = rho*f_sp(:,isp_H3p)
      nmol_l(:,4) = rho*f_sp(:,isp_HeHp)
      nmol_dim    = nmol_l*n0

      ! Compute adimensional total and electron densities for T calculation
      ! (nm/n0 = adimensional metal densities; adds the metal electrons and
      ! nuclei under the eos_metals policy)
      call calc_ne(rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), ne_ad, nm/n0, nmol_l)
      ! The oxygen-chemistry carriers are gas particles as well, and their
      ! oxygen and carbon nuclei have been removed from the metal ion
      ! columns by the ionization solve; leaving them out would make this
      ! energy solve and ioniz_eq disagree about the particle count of the
      ! same state. They are neutral, so calc_ne is unaffected.
      nox_l = 0.0d0
      if (thereis_oxychem) then
         nox_l(:,1) = rho*f_sp(:,isp_OH)
         nox_l(:,2) = rho*f_sp(:,isp_H2O)
         nox_l(:,3) = rho*f_sp(:,isp_CO)
      endif
      nox_dim = nox_l*n0
      ! The He 2^3S column is inside the HeI column (bsp_is_excited_level),
      ! so it is not passed and the triplet branch disappears with it.
      if (thereis_He) then
         if (thereis_oxychem) then
            call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                        rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), n_tot_ad, nm/n0, nmol_l, nox_l)
         else
         call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), rho*f_sp(:,isp_HeI), &
                        rho*f_sp(:,isp_HeII), rho*f_sp(:,isp_HeIII), n_tot_ad, nm/n0, nmol_l)
         endif
      else
         call calc_ntot(rho*f_sp(:,isp_HI), rho*f_sp(:,isp_HII), zero_arr, &
                        zero_arr, zero_arr, n_tot_ad, nm/n0, nmol_l)
      endif

      ! Old temperature (adimensional). The pressure was formed as
      ! (n_tot + n_e) T from the state ioniz_eq returned and the same T it
      ! was evaluated at, and the particle counts here come from that same
      ! state through the same calc_ne/calc_ntot, so this reproduces the
      ! ioniz_eq temperature to round-off.
      T_old = p / (n_tot_ad + ne_ad)
      coupled = present(u_th_old)
      if (present(T_start)) T_old = T_start

      ! The coefficient of the source term, and the ends of the physical
      ! bracket of each cell.
      do j = 1-Ng, N+Ng
         nk = n_tot_ad(j) + ne_ad(j)
         if (nk .gt. 0.0d0) then
            a_coef(j) = dt(j)/nk
         else
            a_coef(j) = -1.0d0        ! rejected by the admissibility test
         endif
         T_floor(j) = max(T_floor_code_min, T_eos_floor_K/T0)
         T_ceil(j)  = max(T_ceiling_K/T0, T_floor(j))
         if (coupled) then
            if (nk .gt. 0.0d0) then
               u_old_pp(j) = u_th_old(j)/nk
               if (present(du_form)) u_old_pp(j) = u_old_pp(j)             &
                                                 - du_form(j)/nk
            else
               u_old_pp(j) = 0.0d0     ! rejected by the admissibility test
            endif
         endif
      enddo

      ! The bracketed solve. Each pass evaluates the cooling once, for the
      ! whole grid, at the temperature each cell's iteration asks for.
      if (coupled) then
         call energy_balance_init(s, T_old, a_coef, heat, cool, T_floor,   &
                                  T_ceil, u_old_in = u_old_pp)
      else
         call energy_balance_init(s, T_old, a_coef, heat, cool, T_floor,   &
                                  T_ceil)
      endif
      do while (s%n_active .gt. 0 .and. s%iters .lt. energy_max_iter)
         T_K = s%T_eval * T0
         call eval_cool(T_K, nhi, nhii, nhei, nheii, nheiii, nm,           &
                        rchiiB, rcheiiB, rcheiiiB, rec_m,                  &
                        a_ion_HI, a_ion_HeI, a_ion_HeII, aion_m,           &
                        cool_dim, nheiTR = nheiTR, nmol = nmol_dim,        &
                        nox = nox_dim)
         cool_iter = cool_dim / q0
         call energy_balance_update(s, cool_iter)
      enddo
      call energy_balance_finish(s)

      ! Verdict.
      ! The verdict of the update is the status of the first failing cell in
      ! index order; n_fail says how many cells failed.
      n_fail  = 0
      verdict = ENERGY_UPDATE_OK
      j_bad   = 0
      do j = 1-Ng, N+Ng
         if (s%status(j) .eq. ENERGY_UPDATE_OK) cycle
         n_fail = n_fail + 1
         if (j_bad .eq. 0) then
            j_bad   = j
            verdict = s%status(j)
         endif
         if (s%status(j) .eq. ENERGY_UPDATE_FLOOR) then
            if (.not. allocated(energy_floor_cell_hits)) then
               allocate(energy_floor_cell_hits(1-Ng:N+Ng))
               energy_floor_cell_hits = 0
            endif
            n_energy_floor_hits = n_energy_floor_hits + 1
            n_energy_floor_hits_family(ledger_family) =                    &
               n_energy_floor_hits_family(ledger_family) + 1
            energy_floor_cell_hits(j) = energy_floor_cell_hits(j) + 1
            if (energy_floor_first_step .lt. 0)                            &
               energy_floor_first_step = step
            energy_floor_last_step = step
         endif
      enddo

      energy_update_last_status = verdict
      energy_update_last_nfail  = n_fail
      energy_update_last_step   = step
      energy_update_last_iters  = s%iters
      if (j_bad .gt. 0) then
         energy_update_last_cell     = j_bad
         energy_update_last_reason   = s%reason(j_bad)
         energy_update_last_residual = s%res_scaled(j_bad)
      else
         energy_update_last_cell     = 0
         energy_update_last_reason   = ENERGY_REASON_NONE
         energy_update_last_residual = maxval(s%res_scaled)
      endif
      if (present(status)) status = verdict

      if (n_fail .gt. 0) then
         ! The compact report of a refused energy update, in the form the
         ! controller's own "REFUSED at" line uses: what failed, where, by
         ! how much, and the volumetric rate the balance is short of. The
         ! caller turns this status into a refusal of the attempted step, so
         ! this is one line pair per refused attempt and not a block.
         ! `missing source` is the rate, in the units of `heat`, that would
         ! have to be supplied for the last iterate to satisfy the balance.
         ! `|R|/N` is the residual of the last iterate over the norm N of
         ! thermal_balance_norm at that iterate (the energy the step moves
         ! plus the arithmetic allowance), NOT over the stored thermal
         ! energy; acceptance needs |R|/N <= energy_res_tol. huge marks a
         ! cell refused on entry or by a non-finite evaluation.
         ! For a floor failure it is what a specified external reservoir
         ! would owe; B1's target system decides whether such a reservoir
         ! exists, and until then it is a diagnostic and not a source term.
         call energy_balance_reason_text(s%reason(j_bad), reason_text)
         write(*,'(a,i0,a,i0,a,i0,a,es9.2,a)')                             &
              ' energy update FAILURE, step ', step, ', ', n_fail,         &
              ' cell(s), first cell ', j_bad, ' (T ',                      &
              T_old(j_bad)*T0, ' K)'
         if (a_coef(j_bad) .gt. 0.0d0) then
            write(*,'(a,a,a,es11.4,a,es11.4,a,es11.4,a,i0)')               &
                 '   reason: ', trim(reason_text), ', |R|/N ',         &
                 s%res_scaled(j_bad), ', last iterate ',                   &
                 s%T_new(j_bad)*T0, ' K, missing source ',                 &
                 s%res_energy(j_bad)/a_coef(j_bad)*q0,                     &
                 ' erg cm^-3 s^-1, iterations ', s%iters
         else
            write(*,'(a,a,a,es11.4,a,es11.4,a,i0)')                        &
                 '   reason: ', trim(reason_text), ', |R|/N ',         &
                 s%res_scaled(j_bad), ', last iterate ',                   &
                 s%T_new(j_bad)*T0, ' K, iterations ', s%iters
         endif
         ! NOTHING IS ASSEMBLED: p, T, cool and u keep the values the call
         ! found, so the caller that refuses the step restores a state this
         ! routine did not touch.
         return
      endif

      ! Every cell converged: assemble the state.
      cool = s%cool_new
      p = (n_tot_ad + ne_ad) * s%T_new
      W(3,:) = p
      if (caloric_mixture_active) then
         do j = 1-Ng, N+Ng
            u(3,j) = 0.5d0*rho(j)*v(j)**2.0                             &
                   + energy_density_from_pressure(j, rho(j), p(j))
         enddo
      else
      u(3,:) = 0.5d0 * rho * v**2.0 + p / (gamma_ad - 1.0d0)
      endif

   end subroutine solve_energy_semi_implicit

end module energy_semi_implicit
