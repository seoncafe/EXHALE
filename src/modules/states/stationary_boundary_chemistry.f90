      module stationary_boundary_chemistry
      ! Eliminate lower ghost chemistry and thermodynamics at a fixed column.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use global_parameters
      use species_table, only: n_mion,isp_HeTR,isp_HeI
      use Conversion, only: U_to_W, admissible_conserved_cell
      use composition, only: get_species_densities, comp_T_from_p
      use BC_Apply, only: Apply_BC
      use reservoir_composition, only: base_reservoir_composition_row
      use base_boundary, only: state_file_ghost_row, state_file_ghost_row_have
      use ionization_equilibrium, only: ioniz_eq, ioniz_eq_ledger,       &
           ioniz_eq_sources_of_held_composition, held_composition_sources, &
           set_ioniz_eq_state_may_be_refused, ieq_state_may_be_refused,   &
           ieq_sweep_ledger_last, ghost_acc_res_cell, ieq_res_tol,       &
           xuv_field_block_cells,xuv_self_field_passes
      use excited_hydrogen, only: excited_H_update
      implicit none
      private
      public :: stationary_ghost_state, boundary_chemistry_result
      public :: boundary_chemistry_of_state
      public :: ghost_fraction_accuracy, ghost_thermal_accuracy
      public :: boundary_thermodynamic_invalid
      public :: boundary_chemistry_record,hold_boundary_chemistry,restore_boundary_chemistry
      integer, parameter :: boundary_chemistry_closed = 0
      integer, parameter :: boundary_chemical_nonroot = 1
      integer, parameter :: boundary_thermodynamic_invalid = 2
      integer, parameter :: boundary_iteration_limit = 3
      integer, parameter :: boundary_max_iterations = 60
      real*8, parameter :: ghost_fraction_accuracy = 1.0d-11
      real*8, parameter :: ghost_thermal_accuracy = 128.0d0*epsilon(1.0d0)
      type :: boundary_chemistry_result
         logical :: converged = .false.
         integer :: status = boundary_iteration_limit
         integer :: iterations = 0
         real*8 :: composition_change = huge(1.0d0)
         real*8 :: temperature_change = huge(1.0d0)
         real*8 :: density_change = huge(1.0d0)
      end type boundary_chemistry_result
      type :: boundary_chemistry_record
         type(boundary_chemistry_result) :: result
         real*8, allocatable :: conserved(:,:),composition(:,:)
      end type boundary_chemistry_record
      type(boundary_chemistry_result), save :: result_last
      real*8, allocatable, save :: u_last(:,:),f_last(:,:)
      contains

      subroutine stationary_ghost_state(u,f_sp,heat,cool,eta,n_part,sweep, &
                                         result,canonical_seed)
      ! The physical conserved column and composition are held exactly.
      ! A stationary source is evaluated at the composition it accompanies.
      ! A failed closure is returned through the existing open-ghost verdict.
      real*8, intent(inout) :: u(3,1-Ng:N+Ng)
      real*8, intent(inout) :: f_sp(1-Ng:N+Ng,n_species)
      real*8, intent(out) :: heat(1-Ng:N+Ng),cool(1-Ng:N+Ng)
      real*8, intent(out) :: eta(1-Ng:N+Ng),n_part(1-Ng:N+Ng)
      type(ioniz_eq_ledger), intent(out) :: sweep
      type(boundary_chemistry_result), intent(out) :: result
      logical, optional, intent(in) :: canonical_seed
      real*8 :: w(3,1-Ng:N+Ng),rho(1-Ng:N+Ng),t(1-Ng:N+Ng)
      real*8 :: nhi(1-Ng:N+Ng),nhii(1-Ng:N+Ng),nhei(1-Ng:N+Ng)
      real*8 :: nheii(1-Ng:N+Ng),nheiii(1-Ng:N+Ng),nhetr(1-Ng:N+Ng)
      real*8 :: ne(1-Ng:N+Ng),ntot(1-Ng:N+Ng),nm(1-Ng:N+Ng,n_mion)
      real*8 :: u_held(3,1:N),f_held(1:N+Ng,n_species)
      real*8 :: f_before(1-Ng:0,n_species),rho_before(1-Ng:0)
      real*8 :: t_before(1-Ng:0),f_reservoir(n_species),exc_change
      real*8 :: ghost_residual(1-Ng:0,n_species),residual_previous(1-Ng:0,n_species)
      real*8 :: residual_delta(1-Ng:0,n_species),f_trial(1-Ng:0,n_species)
      real*8 :: omega,omega_previous,denominator
      logical :: seed_reservoir,held_before,refusable_before,gas_ok,trace,mixture_ok
      logical :: ghost_roots_only,full_sweeps
      integer :: j,it,iteration_cap,ios
      character(len=64) :: diagnostic_seed,trace_env,cap_env,full_env
      result = boundary_chemistry_result()
      sweep = ioniz_eq_ledger()
      heat=0.0d0; cool=0.0d0; eta=0.0d0; n_part=0.0d0
      u_held=u(:,1:N); f_held=f_sp(1:N+Ng,:)
      ! A Newton vector contains physical cells only. Install a positive
      ! provisional density before the composition-dependent EOS is formed;
      ! Apply_BC will replace this provisional state before any chemistry.
      do j=1-Ng,0
         u(:,j)=u(:,1)
      enddo
      do j=N+1,N+Ng
         u(:,j)=u(:,N)
      enddo
      held_before=held_composition_sources
      refusable_before=ieq_state_may_be_refused
      seed_reservoir=.true.
      if (present(canonical_seed)) seed_reservoir=canonical_seed
      ! EXHALE_STATIONARY_BOUNDARY_SEED: Selects state_file or reservoir boundary composition. Unset uses
      ! reservoir; changes the boundary state.
      call get_environment_variable('EXHALE_STATIONARY_BOUNDARY_SEED', &
                                     diagnostic_seed)
      if (trim(diagnostic_seed).eq.'state_file') then
         if (.not.state_file_ghost_row_have) then
            write(*,'(A)') ' ERROR: stationary boundary seed file unavailable'
            stop 1
         endif
         f_sp(1-Ng:0,:)=state_file_ghost_row
         seed_reservoir=.false.
      else if (len_trim(diagnostic_seed).ne.0 .and. &
               trim(diagnostic_seed).ne.'reservoir') then
         write(*,'(A)') ' ERROR: unknown stationary boundary seed'
         stop 1
      endif
      if (seed_reservoir) then
         call base_reservoir_composition_row(f_sp(1,:),f_reservoir)
         do j=1-Ng,0
            f_sp(j,:)=f_reservoir
         enddo
      endif
      ! EXHALE_BOUNDARY_CLOSURE_TRACE: Set to 1 to print the boundary closure trace. Unset omits the
      ! trace; diagnostic only.
      call get_environment_variable('EXHALE_BOUNDARY_CLOSURE_TRACE',trace_env)
      trace=trim(trace_env).eq.'1'
      ! EXHALE_BOUNDARY_FULL_SWEEPS: Set to 1 to require full boundary chemistry sweeps. Unset follows
      ! the XUV controls; can change the evaluated closure.
      call get_environment_variable('EXHALE_BOUNDARY_FULL_SWEEPS',full_env)
      full_sweeps=trim(full_env).eq.'1'.or.xuv_field_block_cells.ne.0.or. &
                   xuv_self_field_passes.ne.1
      call ioniz_eq_sources_of_held_composition(.true.)
      call set_ioniz_eq_state_may_be_refused(.true.)
      omega_previous=1.0d0
      residual_previous=0.0d0
      iteration_cap=boundary_max_iterations
      ! EXHALE_BOUNDARY_CLOSURE_MAXITER: Overrides the boundary closure iteration cap within its allowed
      ! range. Unset uses boundary_max_iterations; can change closure acceptance.
      call get_environment_variable('EXHALE_BOUNDARY_CLOSURE_MAXITER',cap_env)
      if (len_trim(cap_env).gt.0) then
         read(cap_env,*,iostat=ios) iteration_cap
         if (ios.ne.0.or.iteration_cap.lt.1.or. &
             iteration_cap.gt.boundary_max_iterations) then
            write(*,'(A)') ' ERROR: invalid diagnostic boundary iteration cap'
            stop 1
         endif
      endif
      do it=1,iteration_cap
         result%iterations=it
         call thermodynamic_state
         if (.not.gas_ok) then
            result%status=boundary_thermodynamic_invalid
            exit
         endif
         f_before=f_sp(1-Ng:0,:)
         rho_before=rho(1-Ng:0); t_before=t(1-Ng:0)
         if (use_excited_H) call excited_H_update(t,rho,f_sp,w(2,:),exc_change)
         ghost_roots_only=it.gt.1.and..not.full_sweeps
         call ioniz_eq(t,rho,f_sp,heat,cool,eta,sweep, &
              lower_ghosts_only=ghost_roots_only,advance_composition_history=.false.)
         if (sweep%n_ghost_open.gt.0 .or. sweep%n_nonfinite.gt.0 .or. &
             any(.not.ieee_is_finite(ghost_acc_res_cell)).or. &
             any(ghost_acc_res_cell.gt.ieq_res_tol)) then
            call thermodynamic_state
            result%status=boundary_chemical_nonroot
            if (.not.gas_ok) result%status=boundary_thermodynamic_invalid
            exit
         endif
         call thermodynamic_state
         if (.not.gas_ok) then
            result%status=boundary_thermodynamic_invalid
            exit
         endif
         result%composition_change=maxval(abs(f_sp(1-Ng:0,:)-f_before))
         result%temperature_change=maxval(abs(t(1-Ng:0)-t_before) &
                                         /max(abs(t_before),tiny(1.0d0)))
         result%density_change=maxval(abs(rho(1-Ng:0)-rho_before) &
                                     /max(abs(rho_before),tiny(1.0d0)))
         if (trace) write(*,'(A,I3,3(A,ES12.4))') &
              ' (boundary chemistry) iteration ',it, &
              ' df=',result%composition_change, &
              ' dT/T=',result%temperature_change, &
              ' drho/rho=',result%density_change
         if (result%composition_change.le.ghost_fraction_accuracy .and. &
             result%temperature_change.le.ghost_thermal_accuracy .and. &
             result%density_change.le.ghost_thermal_accuracy) then
            if (ghost_roots_only) then
               ! Refresh every physical source/rate record before returning.
               ! Verify that a complete sweep leaves this boundary closed.
               f_before=f_sp(1-Ng:0,:)
               rho_before=rho(1-Ng:0); t_before=t(1-Ng:0)
               call ioniz_eq(t,rho,f_sp,heat,cool,eta,sweep, &
                              advance_composition_history=.false.)
               call thermodynamic_state
               if (.not.gas_ok) then
                  result%status=boundary_thermodynamic_invalid
                  exit
               endif
               if (sweep%n_ghost_open.gt.0.or.sweep%n_nonfinite.gt.0.or. &
                   any(.not.ieee_is_finite(ghost_acc_res_cell)).or. &
                   any(ghost_acc_res_cell.gt.ieq_res_tol)) then
                  result%status=boundary_chemical_nonroot
                  exit
               endif
               result%composition_change=maxval(abs(f_sp(1-Ng:0,:)-f_before))
               result%temperature_change=maxval(abs(t(1-Ng:0)-t_before) &
                                      /max(abs(t_before),tiny(1.0d0)))
               result%density_change=maxval(abs(rho(1-Ng:0)-rho_before) &
                                      /max(abs(rho_before),tiny(1.0d0)))
            endif
            if (result%composition_change.le.ghost_fraction_accuracy.and. &
                result%temperature_change.le.ghost_thermal_accuracy.and. &
                result%density_change.le.ghost_thermal_accuracy) then
               result%converged=.true.
               result%status=boundary_chemistry_closed
               exit
            endif
         endif
         if (it.lt.iteration_cap) then
            ! Scalar Aitken relaxation of the conservative composition map.
            ! Salz has a measured contraction near 0.90: 60 plain passes
            ! leave dT/T about 7e-11. Mixing two nucleus-conserving rows
            ! preserves their elemental constraints. The raw map increment
            ! above, never the relaxed increment, decides convergence.
            ghost_residual=f_sp(1-Ng:0,:)-f_before
            omega=1.0d0
            if (it.gt.1) then
               residual_delta=ghost_residual-residual_previous
               denominator=sum(residual_delta**2)
               if (denominator.gt.tiny(1.0d0)) omega=max(1.0d0,min(16.0d0, &
                    -omega_previous*sum(residual_previous*residual_delta)/denominator))
            endif
            do while (omega.gt.1.0d0)
               f_trial=f_before+omega*ghost_residual
               mixture_ok=all(ieee_is_finite(f_trial)).and.all(f_trial.ge.0.0d0)
               if (thereis_HeITR) mixture_ok=mixture_ok.and. &
                    all(f_trial(:,isp_HeTR).le.f_trial(:,isp_HeI))
               if (mixture_ok) exit
               omega=1.0d0+0.5d0*(omega-1.0d0)
               if (omega-1.0d0.lt.32.0d0*epsilon(1.0d0)) omega=1.0d0
            enddo
            if (omega.gt.1.0d0) f_sp(1-Ng:0,:)=f_trial
            if (trace) write(*,'(A,ES12.4)') ' (boundary chemistry) omega=',omega
            residual_previous=ghost_residual
            omega_previous=omega
         endif
      enddo
      if (.not.result%converged) then
         sweep%n_ghost_open=max(1,sweep%n_ghost_open)
         ieq_sweep_ledger_last%n_ghost_open=sweep%n_ghost_open
         write(*,'(A,I0,A,I0)') ' (boundary chemistry) refused status=', &
              result%status,' iterations=',result%iterations
      endif
      call ioniz_eq_sources_of_held_composition(held_before)
      call set_ioniz_eq_state_may_be_refused(refusable_before)
      if (any(u(:,1:N).ne.u_held) .or. any(f_sp(1:N+Ng,:).ne.f_held)) then
         write(*,'(A)') ' ERROR: stationary ghost closure changed held column'
         stop 1
      endif
      u_last=u
      f_last=f_sp
      result_last=result
      contains
      subroutine thermodynamic_state
      ! Install this composition before either the boundary or the EOS reads it.
      call get_species_densities(u(1,:),f_sp,nhi,nhii,nhei,nheii,nheiii, &
                                 nhetr,nm,ne,ntot)
      call Apply_BC(u)
      gas_ok=.true.
      do j=1-Ng,N+Ng
         if (.not.admissible_conserved_cell(u(:,j))) gas_ok=.false.
      enddo
      if (.not.gas_ok) return
      call U_to_W(u,w)
      rho=w(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,nheiii, &
                                 nhetr,nm,ne,ntot)
      call comp_T_from_p(w(3,:),ntot,ne,t)
      n_part=ntot+ne
      gas_ok=all(ieee_is_finite(t)).and.all(t.gt.0.0d0)
      end subroutine thermodynamic_state
      end subroutine stationary_ghost_state

      subroutine boundary_chemistry_of_state(u,f_sp,result,known)
      ! A certificate may use a closure only for the exact state it measured.
      real*8, intent(in) :: u(3,1-Ng:N+Ng)
      real*8, intent(in) :: f_sp(1-Ng:N+Ng,n_species)
      type(boundary_chemistry_result), intent(out) :: result
      logical, intent(out) :: known
      result=boundary_chemistry_result()
      known=.false.
      if (.not.allocated(u_last).or..not.allocated(f_last)) return
      if (any(u.ne.u_last).or.any(f_sp.ne.f_last)) return
      result=result_last
      known=.true.
      end subroutine boundary_chemistry_of_state
      subroutine hold_boundary_chemistry(record)
      type(boundary_chemistry_record), intent(out) :: record
      record%result=result_last
      if (allocated(u_last)) record%conserved=u_last
      if (allocated(f_last)) record%composition=f_last
      end subroutine hold_boundary_chemistry
      subroutine restore_boundary_chemistry(record)
      type(boundary_chemistry_record), intent(in) :: record
      result_last=record%result
      if (allocated(record%conserved)) then
         u_last=record%conserved
      else if (allocated(u_last)) then
         deallocate(u_last)
      endif
      if (allocated(record%composition)) then
         f_last=record%composition
      else if (allocated(f_last)) then
         deallocate(f_last)
      endif
      end subroutine restore_boundary_chemistry
      end module stationary_boundary_chemistry
