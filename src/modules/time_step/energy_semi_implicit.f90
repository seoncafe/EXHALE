module energy_semi_implicit
   ! Implicit source step for the thermal energy of a cell, at fixed
   ! composition and fixed heating.
   !
   ! WHERE IT SITS.  Since B3c this is the INNER step of the coupled
   ! temperature-composition source step of
   ! docs/b1_target_system_20260906.md T1.6: the marching loop alternates it
   ! with the composition sweep until the pair stops moving, so the
   ! composition is fixed for one pass and not for the step.  Two things
   ! follow, and both are optional arguments of solve_energy_semi_implicit
   ! so that a caller advancing the temperature alone is unaffected: the
   ! temperature the sources were evaluated at is handed in rather than
   ! recovered from the pressure (T_start), and the energy the step starts
   ! from is the thermal energy the cell HAD at its OLD composition
   ! (u_th_old), not that energy rebuilt at the new particle count.  The
   ! rebuild was the composition projection, and it changed the conserved
   ! energy with no source behind it.
   !
   ! THE EQUATION.  With the composition frozen over the step and the
   ! photoheating rate frozen at the state the step starts from, the energy
   ! of one cell obeys
   !
   !    E(T) - E_old  =  dt [ heat - cool(T) ] ,
   !
   ! with E the internal energy of the cell from the caloric equation of
   ! state the run uses (caloric_eos: (3/2)(n_tot+n_e) k T plus the
   ! rovibrational energy of the bound H2 in the cell, or the monatomic
   ! T/(gamma_ad - 1) when the mixture branch is inactive).  Everything below
   ! is written per (n_tot + n_e), which is the form internal_energy_per_particle
   ! returns, so the residual solved is
   !
   !    R(T) = u(T) - u(T_old) - a [ heat - cool(T) ] ,   a = dt/(n_tot+n_e) .
   !
   ! HOW IT IS SOLVED.  A safeguarded Newton-bisection hybrid on a physical
   ! bracket.  The bracket ends are the floor of the equation-of-state
   ! validity and a ceiling above the heating-only update; the search starts
   ! from the explicit-Euler point, brackets the root as it goes, and once
   ! both ends are known it takes a Newton step whenever that step stays
   ! strictly inside the bracket and the bracket keeps halving, and bisects
   ! otherwise.  The Newton derivative is exact in the equation-of-state term
   ! (heat_capacity_per_particle) and a secant estimate in dcool/dT, which
   ! only affects the path: the safeguard owns the convergence.
   !
   ! WHAT IS TESTED BEFORE THE UPDATE IS RETURNED.  The residual at the
   ! RETURNED temperature, against |R|/scale with
   ! scale = |E_old| + dt(|heat| + |cool|); admissibility (inside the bracket,
   ! finite).  The cooling written back is cool(T_returned), evaluated there
   ! and not at a previous iterate.  A cell that does not meet the tolerance
   ! within the iteration budget, that has no bracket, or that reaches the
   ! equation-of-state floor is a FAILURE with a named reason.  There is no
   ! clamp: an unconverged temperature is never returned as a state.
   !
   ! WHAT A FAILURE DOES.  It is REPORTED, and nothing else happens here:
   ! the routine writes one compact diagnostic line pair and returns its
   ! status.  The caller refuses the whole attempted step on it (the B3a
   ! controller restores the checkpoint and retakes the step at half dt), in
   ! BOTH run modes: the leniency the initialization mode grants to a fixed
   ! point that was not reached is about an iterate that exists, while a
   ! failed energy update has produced no temperature at all, so there is no
   ! state to march on.  NOTHING IS ASSEMBLED FROM A FAILED UPDATE: p, T,
   ! cool and u are left as the call found them, and a caller must not read
   ! them after a status other than ENERGY_UPDATE_OK.
   !
   ! THE FLOOR IS NOT A RESERVOIR.  Reaching the lower bracket end means the
   ! balance asks for a temperature below the range in which the equation of
   ! state and the rate fits are defined.  A floor that represents real
   ! external heating would have to be a physically specified reservoir with
   ! its heating budgeted, and that specification belongs to B1's target
   ! system, not here.  Until it exists, the floor is a failure.
   use global_parameters
   use species_table, only: n_mion, mion_fsp,                        &
                            isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                            isp_HeIII, isp_HeTR,                      &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                            isp_OH, isp_H2O, isp_CO
   use utils
   use utils_ion_eq
   use caloric_eos, only: caloric_mixture_active, molecular_cell,      &
                          internal_energy_per_particle,                &
                          heat_capacity_per_particle,                  &
                          temperature_from_energy_per_particle,        &
                          energy_density_from_pressure

   implicit none

   ! Status of the temperature solve of one cell, and of the update as a
   ! whole (the worst cell). Reported through the optional `status` argument
   ! of solve_energy_semi_implicit and through energy_update_last_status, so
   ! the physical-step context of docs/a2_certification_contract_20260906.md
   ! section 3 and the B3a controller can read a verdict rather than infer
   ! one from the returned numbers.
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

   ! Acceptance of the returned state.
   !   energy_res_tol   |R|/scale of the returned temperature. The residual
   !                    is an energy per (n_tot+n_e) of order T ~ 1 in code
   !                    units and its floating-point resolution is ~1e-16 of
   !                    that, so 1e-9 is far above the arithmetic noise and
   !                    far below any physically meaningful energy error of
   !                    the step.
   !   energy_T_tol     relative width of the bracket below which the root is
   !                    located as precisely as the arithmetic allows. A
   !                    bracket this narrow moves the residual by about
   !                    c_v * energy_T_tol * T, which is below
   !                    energy_res_tol * scale, so a collapsed bracket that
   !                    still fails the residual test means cool(T) is not
   !                    continuous there, and that is reported as such.
   !   energy_max_iter  iteration budget. Bisection alone needs about 37
   !                    halvings to reach energy_T_tol from a bracket as wide
   !                    as the temperature itself.
   real*8,  parameter :: energy_res_tol     = 1.0d-9
   real*8,  parameter :: energy_T_tol       = 1.0d-11
   integer, parameter :: energy_max_iter    = 60
   ! Absolute floor of the residual scale, so that a cell with no energy at
   ! all divides by a positive number.
   real*8,  parameter :: energy_scale_floor = 1.0d-30

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
   !   The upper end. For a cell holding H2 it is the top of that same
   !   rovibrational table (caloric_eos); for a cell without molecules the
   !   caloric equation of state is exact at any temperature and the ceiling
   !   only has to be far above anything a photoionized wind reaches, so that
   !   "no bracket" reports a heating rate no physical cooling can balance
   !   rather than a bracket set too tight.
   real*8, parameter :: T_ceiling_mol_K    = 5.0d4
   real*8, parameter :: T_ceiling_atomic_K = 1.0d7

   ! TEST HOOK, default off, named for what it does. It makes the update
   ! REPORT a failure at one step count, for a given number of leading
   ! calls, so that the controller's handling of a failed energy update can
   ! be exercised on a configuration that does not happen to produce one.
   ! Nothing else changes: the solve runs as it would have, and only the
   ! status of the first physical cell is overwritten, so what the caller
   ! then does is what it would do for a real failure. Set from
   ! EXHALE_ENERGY_FAIL_AT_STEP and EXHALE_ENERGY_FAIL_CALLS
   ! (energy_update_read_environment).
   integer, save :: energy_fail_at_step = -1
   integer, save :: energy_fail_calls   = 1
   integer, save :: energy_fail_served  = 0

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
   ! (docs/Update_EXHALE_stage1.md section 113). These now count floor
   ! FAILURES, that is attempts: no accepted state sits on the floor, because
   ! the clamp that used to produce one is gone.
   !   hits        total activations over the run
   !   first/last  first and last step on which the floor was reached
   !   cell_hits   activations of each cell, so the number of DISTINCT cells
   !               that ever touched the floor can be reported
   integer, save :: n_energy_floor_hits = 0
   ! The same activations split by ledger family
   ! (docs/a0_run_mode_contract_20260906.md section 5): index
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
      real*8, allocatable :: u_old(:)      ! internal energy per particle at T_old
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
      real*8, allocatable :: res_scaled(:)   ! |R|/scale at the returned T
      real*8, allocatable :: res_energy(:)   ! R itself, per particle
      real*8, allocatable :: scale(:)
      integer :: n_active = 0
      integer :: iters    = 0
   end type energy_balance_state

contains

   ! ------------------------------------------------------!

   logical function energy_is_finite(x)
   ! True for a normal or subnormal double, false for a NaN or an infinity.
   real*8, intent(in) :: x
   energy_is_finite = (x .eq. x) .and. (abs(x) .le. huge(1.0d0))
   end function energy_is_finite

   ! ------------------------------------------------------!

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

   ! ------------------------------------------------------!

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
   ! u_old_in IS THE ENERGY THE STEP STARTS FROM, per (n_tot + n_e) of the
   ! composition the solve is being run at.  It is a separate argument
   ! because the temperature the step starts from and the energy it starts
   ! from need not belong to the same composition.  In the coupled source
   ! step of docs/b1_target_system_20260906.md T1.6 they do not: the energy
   ! is the thermal energy the cell HAD, at its old composition, less the
   ! formation-energy change of the composition change, while T_old is the
   ! temperature the composition was solved at.  Absent, the routine falls
   ! back to internal_energy_per_particle(T_old), which is the same number
   ! whenever the composition did not move, and is what a caller that
   ! advances the temperature at frozen composition wants.
   type(energy_balance_state), intent(out) :: s
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T_old, a, heat, cool0
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T_floor, T_ceil
   real*8, dimension(1-Ng:N+Ng), intent(in), optional :: u_old_in
   real*8  :: R0, sc, u_expl, T_guess
   integer :: j

   allocate(s%a(1-Ng:N+Ng), s%heat(1-Ng:N+Ng), s%T_old(1-Ng:N+Ng),         &
            s%u_old(1-Ng:N+Ng), s%T_floor(1-Ng:N+Ng), s%T_ceil(1-Ng:N+Ng), &
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
   s%scale   = energy_scale_floor
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

      ! R at T_old.  The first two terms cancel identically only when the
      ! energy the step starts from is the energy of T_old at this
      ! composition; with u_old_in they do not, and what is left is the
      ! energy the composition change put into or took out of the cell.
      R0 = internal_energy_per_particle(j, T_old(j)) - s%u_old(j)          &
         + a(j)*(cool0(j) - heat(j))
      sc = abs(s%u_old(j)) + a(j)*(abs(heat(j)) + abs(cool0(j)))
      sc = max(sc, energy_scale_floor)
      s%scale(j) = sc

      if (abs(R0) .le. energy_res_tol*sc) then
         s%T_new(j)      = T_old(j)
         s%cool_new(j)   = cool0(j)
         s%res_scaled(j) = abs(R0)/sc
         s%res_energy(j) = R0
         cycle
      endif

      ! One bracket end is known from the sign at T_old: R is the energy the
      ! state has minus the energy the balance allows it, so R > 0 means the
      ! root is colder.
      if (R0 .gt. 0.0d0) then
         s%hi(j) = T_old(j);  s%has_hi(j) = .true.
      else
         s%lo(j) = T_old(j);  s%has_lo(j) = .true.
      endif

      ! First trial: the explicit-Euler point, the temperature the step would
      ! land on if the cooling did not follow T.
      u_expl  = s%u_old(j) + a(j)*(heat(j) - cool0(j))
      T_guess = temperature_from_energy_per_particle(j, u_expl)
      if (.not. energy_is_finite(T_guess)) T_guess = T_old(j)
      T_guess = min(max(T_guess, T_floor(j)), T_ceil(j))
      ! It must stay on the far side of the end already known, or the
      ! iteration would learn nothing from it.
      if (s%has_hi(j) .and. T_guess .ge. s%hi(j))                          &
         T_guess = max(T_floor(j), 0.5d0*s%hi(j))
      if (s%has_lo(j) .and. T_guess .le. s%lo(j))                          &
         T_guess = min(T_ceil(j), 2.0d0*s%lo(j))
      s%T_eval(j) = T_guess
      s%active(j) = .true.
      s%n_active  = s%n_active + 1
   enddo

   end subroutine energy_balance_init

   ! ------------------------------------------------------!

   subroutine energy_balance_update(s, cool_at_T_eval)
   ! Absorb the cooling the caller evaluated at s%T_eval, accept or reject
   ! each active cell, and choose the next temperature for the ones that
   ! remain.
   type(energy_balance_state), intent(inout) :: s
   real*8, dimension(1-Ng:N+Ng), intent(in) :: cool_at_T_eval
   real*8  :: T_e, C, R, sc, cv, dR, T_n, w
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

      R = internal_energy_per_particle(j, T_e) - s%u_old(j)                &
          - s%a(j)*(s%heat(j) - C)
      if (.not. energy_is_finite(R)) then
         call energy_balance_fail(s, j, ENERGY_UPDATE_NONFINITE,           &
                                  ENERGY_REASON_NONFINITE_COOL, T_e, C,    &
                                  huge(1.0d0), huge(1.0d0))
         cycle
      endif

      sc = abs(s%u_old(j)) + s%a(j)*(abs(s%heat(j)) + abs(C))
      sc = max(sc, s%scale(j), energy_scale_floor)
      s%scale(j) = sc

      ! Acceptance: the residual of the state that would be returned, with
      ! the cooling evaluated at the temperature that would be returned.
      if (abs(R) .le. energy_res_tol*sc) then
         s%T_new(j)      = T_e
         s%cool_new(j)   = C
         s%res_scaled(j) = abs(R)/sc
         s%res_energy(j) = R
         s%active(j)     = .false.
         s%n_active      = s%n_active - 1
         cycle
      endif

      ! Bracket update.
      if (R .lt. 0.0d0) then
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
      T_n = T_e - R/dR

      if (s%has_lo(j) .and. s%has_hi(j)) then
         w = s%hi(j) - s%lo(j)
         if (w .le. energy_T_tol*max(abs(s%hi(j)), s%T_floor(j))) then
            ! The root is located to the precision of the arithmetic and the
            ! residual is still above tolerance: cool(T) is not continuous
            ! across this bracket. Reported, never returned as a state.
            call energy_balance_fail(s, j, ENERGY_UPDATE_ITER_CAP,         &
                                     ENERGY_REASON_BRACKET_COLLAPSED,      &
                                     T_e, C, abs(R)/sc, R)
            cycle
         endif
         if (.not. energy_is_finite(T_n) .or. T_n .le. s%lo(j)             &
             .or. T_n .ge. s%hi(j)) then
            T_n = 0.5d0*(s%lo(j) + s%hi(j))
            s%n_slow(j) = 0
         else
            if (w .gt. 0.5d0*s%width_prev(j)) then
               s%n_slow(j) = s%n_slow(j) + 1
            else
               s%n_slow(j) = 0
            endif
            if (s%n_slow(j) .ge. 2) then
               T_n = 0.5d0*(s%lo(j) + s%hi(j))
               s%n_slow(j) = 0
            endif
         endif
         s%width_prev(j) = w
      else if (.not. s%has_hi(j)) then
         ! Still searching upward: the balance wants a hotter cell than any
         ! temperature tried so far.
         if (.not. energy_is_finite(T_n) .or. T_n .le. T_e)                &
            T_n = 2.0d0*T_e + 1.0d-3
         T_n = min(T_n, 4.0d0*T_e + 1.0d-3)
         T_n = min(T_n, s%T_ceil(j))
         if (T_n .le. T_e) then
            ! At the top of the range where the equation of state and the
            ! cooling fits are defined and still short of the balance: no
            ! bracket exists inside the physical range.
            call energy_balance_fail(s, j, ENERGY_UPDATE_NO_BRACKET,       &
                                     ENERGY_REASON_CEILING, T_e, C,        &
                                     abs(R)/sc, R)
            cycle
         endif
      else
         ! Still searching downward.
         if (.not. energy_is_finite(T_n) .or. T_n .ge. T_e)                &
            T_n = 0.5d0*T_e
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
                                     abs(R)/sc, R)
            cycle
         endif
      endif

      s%T_eval(j) = T_n
   enddo

   end subroutine energy_balance_update

   ! ------------------------------------------------------!

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

   ! ------------------------------------------------------!

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

   ! ------------------------------------------------------!

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
      text = 'non-finite cooling or residual'
   case default
      text = 'converged'
   end select
   end subroutine energy_balance_reason_text

   ! ------------------------------------------------------!

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
      ! (docs/b1_target_system_20260906.md T1.6).  All three are absent when
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
      logical :: is_mol_cell
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
         is_mol_cell = .false.
         if (caloric_mixture_active .and. allocated(molecular_cell))       &
            is_mol_cell = molecular_cell(j)
         T_floor(j) = max(T_floor_code_min, T_eos_floor_K/T0)
         if (is_mol_cell) then
            T_ceil(j) = T_ceiling_mol_K/T0
         else
            T_ceil(j) = T_ceiling_atomic_K/T0
         endif
         T_ceil(j) = max(T_ceil(j), T_floor(j))
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

      ! The forced failure of the test hook, if it was asked for.
      if (energy_fail_at_step .ge. 0 .and. step .eq. energy_fail_at_step   &
          .and. energy_fail_served .lt. energy_fail_calls) then
         energy_fail_served = energy_fail_served + 1
         s%status(1)     = ENERGY_UPDATE_NONFINITE
         s%reason(1)     = ENERGY_REASON_NONFINITE_INPUT
         s%res_scaled(1) = huge(1.0d0)
      endif

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
                 '   reason: ', trim(reason_text), ', |R|/scale ',         &
                 s%res_scaled(j_bad), ', last iterate ',                   &
                 s%T_new(j_bad)*T0, ' K, missing source ',                 &
                 s%res_energy(j_bad)/a_coef(j_bad)*q0,                     &
                 ' erg cm^-3 s^-1, iterations ', s%iters
         else
            write(*,'(a,a,a,es11.4,a,es11.4,a,i0)')                        &
                 '   reason: ', trim(reason_text), ', |R|/scale ',         &
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

   ! ------------------------------------------------------!

   subroutine energy_update_read_environment()
   ! The test hook of this module, off unless asked for, so that a run
   ! without the two variables is the run an unguarded build gives.
   character(len=64) :: env
   call get_environment_variable('EXHALE_ENERGY_FAIL_AT_STEP', env)
   if (len_trim(env) .gt. 0) read(env,*) energy_fail_at_step
   call get_environment_variable('EXHALE_ENERGY_FAIL_CALLS', env)
   if (len_trim(env) .gt. 0) read(env,*) energy_fail_calls
   end subroutine energy_update_read_environment

end module energy_semi_implicit
