      module adiabatic_entropy_transport
      ! THE INTERNAL ENERGY OF A COLD, HIGHLY SUPERSONIC EXPANDING FLOW, FROM
      ! ITS ADIABAT (the dual-energy method: Ryu et al. 1993, ApJ 414, 1;
      ! Bryan et al. 1995, Comput. Phys. Commun. 89, 149).
      !
      ! THE PROBLEM. The marching scheme carries the total energy density
      ! E = rho e + rho v^2/2 and recovers the internal energy as the
      ! difference rho e = E - rho v^2/2. Where the flow is highly
      ! supersonic the two terms of the difference are far apart:
      ! rho e/(rho v^2/2) = 2/(gamma (gamma-1) M^2), 1.1e-3 at M = 40. The
      ! difference then carries the TRUNCATION error of the kinetic energy
      ! flux of the scheme (not round-off), amplified by the same ratio, and
      ! the recovered temperature is noise. MEASURED on the Yan 2022
      ! WASP-52 b, He/H 98/2, run (frozen EXHALE_rochefix_c5c106f1, last
      ! accepted state before step 7272): between 9.03 and 9.14 R_p, where
      ! rho e/(rho v^2/2) falls from 7.0e-3 to 8.2e-4, the adiabat p/rho^gamma
      ! DROPS by a factor 7.2 (1.00e-15 to 1.38e-16 in cgs), although the
      ! net source there is positive (heating 1e-16 against cooling 2e-19 erg
      ! cm^-3 s^-1) and an adiabatic flow with a positive source can only
      ! raise it. The temperature then reaches the 0.01 T0 floor and the
      ! energy update refuses the state.
      !
      ! THE METHOD. Between the heating and cooling updates (which are
      ! operator split and applied afterwards to the total energy) the
      ! hydrodynamic step is adiabatic, and away from shocks the adiabat
      ! K = p/rho^gamma is carried with the gas:
      !
      !    d(rho K)/dt + div(rho K v) = 0 .
      !
      ! It is advanced here on the face mass fluxes, volumes, time step and
      ! SSP-RK3 weights of the mass row of each stage, with the donor-cell
      ! value of K on each face (first order; the variable is used only
      ! where the flow is smooth and expanding). At the end of the step,
      ! in a cell where
      !   (i)  rho e < eta E, with eta = thermal_fraction_for_adiabat (the
      !        internal energy is a small difference of the total and the
      !        kinetic energy; eta is set above the relative truncation
      !        error of the scheme's kinetic energy, of order dr/r = 1e-2
      !        on the stretched grids used, so the method takes over only
      !        where the difference can no longer be trusted),
      !   (ii) the flow expands, div v >= 0 (no shock: the adiabat is not
      !        conserved through a shock, and the total-energy form is the
      !        one that is right there),
      !   (iii) the gas has the constant adiabatic index gamma_ad (atomic;
      !        the caloric equation of state of molecular gas has no single
      !        adiabat p/rho^gamma),
      ! the internal energy is taken from the adiabat, rho e = K rho^gamma/
      ! (gamma - 1), and the total energy is reset to it plus the kinetic
      ! energy. The adiabat is re-taken from the state at the start of
      ! every step, so it carries no history of its own: the heating and
      ! cooling of each step enter through the total energy as before.
      ! Total energy is then not conserved in those cells, by an amount
      ! below eta of their total energy, which is the size of the
      ! truncation error it replaces.
      !
      ! VALIDITY. Flows with M >~ 13 (eta = 1e-2) that expand. A run whose
      ! every cell is below that Mach number or compressive is unchanged.

      use global_parameters
      use caloric_eos, only: adiabatic_index_from_state
      use species_advective_transport, only: species_flux_divergence

      implicit none

      ! Thermal fraction of the total energy below which the internal energy
      ! is taken from the adiabat (criterion (i) above).
      real*8, parameter :: thermal_fraction_for_adiabat = 1.0d-2

      ! The adiabat K = p/rho^gamma at the start of the step and after the
      ! last stage. Zero marks a cell whose adiabat is not defined (no
      ! single gamma, or p <= 0).
      real*8, allocatable, dimension(:) :: adiabat_step, adiabat_stage

      ! Cells whose internal energy was taken from the adiabat, summed over
      ! the accepted steps of the run.
      integer :: n_cells_internal_energy_from_adiabat = 0

      contains

      !-----------------------------------------------!

      subroutine adiabat_of_conserved_state(u_in, K)
      ! K = p/rho^gamma_ad of each cell of a conserved state, where the cell
      ! has the constant adiabatic index gamma_ad and a positive pressure;
      ! zero elsewhere.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: K
      real*8  :: rho_c, p_c
      integer :: j
      K = 0.0d0
      do j = 1-Ng, N+Ng
         rho_c = u_in(1,j)
         if (.not. (rho_c .gt. 0.0d0)) cycle
         p_c = (gamma_ad - 1.0d0)*(u_in(3,j) - 0.5d0*u_in(2,j)*u_in(2,j)/rho_c)
         if (.not. (p_c .gt. 0.0d0)) cycle
         if (adiabatic_index_from_state(j, rho_c, p_c) .ne. gamma_ad) cycle
         K(j) = p_c/rho_c**gamma_ad
      enddo
      end subroutine adiabat_of_conserved_state

      !-----------------------------------------------!

      subroutine adiabat_advection_stage(stage, rho_n, u_in, rho_new,     &
                                         Frho, dt_loc)
      ! One SSP-RK3 stage of (rho K), with the same combination of the
      ! step-start and stage-input states as the mass row of the stage
      ! (species_advective_update states it). Cells 1 .. N are advanced; the
      ! ghost values come from the ghost states of the stage input, which
      ! the boundary conditions set.
      integer, intent(in) :: stage
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: rho_n, rho_new
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: Frho, dt_loc

      real*8, dimension(1-Ng:N+Ng) :: K_in, K_ghost, Kf, FK, divFK
      real*8  :: rhoK
      integer :: j

      if (.not. allocated(adiabat_step)) then
         allocate(adiabat_step(1-Ng:N+Ng), adiabat_stage(1-Ng:N+Ng))
      endif
      if (size(adiabat_step) .ne. N + 2*Ng) then
         deallocate(adiabat_step, adiabat_stage)
         allocate(adiabat_step(1-Ng:N+Ng), adiabat_stage(1-Ng:N+Ng))
      endif

      call adiabat_of_conserved_state(u_in, K_ghost)
      if (stage .eq. 1) then
         adiabat_step  = K_ghost
         adiabat_stage = K_ghost
      endif
      K_in = adiabat_stage
      K_in(1-Ng:0)    = K_ghost(1-Ng:0)
      K_in(N+1:N+Ng)  = K_ghost(N+1:N+Ng)

      ! Donor-cell face value of K on each face j (between cells j, j+1).
      Kf = 0.0d0
      do j = 1-Ng, N+Ng-1
         if (Frho(j) .ge. 0.0d0) then
            Kf(j) = K_in(j)
         else
            Kf(j) = K_in(j+1)
         endif
      enddo
      FK = Frho*Kf
      call species_flux_divergence(FK, divFK)

      do j = 1, N
         select case (stage)
         case (1)
            rhoK = rho_n(j)*adiabat_step(j) - dt_loc(j)*divFK(j)
         case (2)
            rhoK = (3.0d0*rho_n(j)*adiabat_step(j) + u_in(1,j)*K_in(j)  &
                    - dt_loc(j)*divFK(j))/4.0d0
         case default
            rhoK = (rho_n(j)*adiabat_step(j)                             &
                    + 2.0d0*(u_in(1,j)*K_in(j) - dt_loc(j)*divFK(j)))/3.0d0
         end select
         if (rho_new(j) .gt. 0.0d0 .and. adiabat_step(j) .gt. 0.0d0 .and. &
             K_in(j) .gt. 0.0d0) then
            adiabat_stage(j) = max(rhoK, 0.0d0)/rho_new(j)
         else
            adiabat_stage(j) = 0.0d0
         endif
      enddo
      adiabat_stage(1-Ng:0)   = K_ghost(1-Ng:0)
      adiabat_stage(N+1:N+Ng) = K_ghost(N+1:N+Ng)

      end subroutine adiabat_advection_stage

      !-----------------------------------------------!

      subroutine internal_energy_from_adiabat(u, n_changed)
      ! At the end of an accepted hydrodynamic step: in the cells that meet
      ! (i)-(iii) of the module header, the internal energy from the
      ! advected adiabat. n_changed is the number of cells reset.
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      integer, intent(out) :: n_changed
      real*8  :: rho_c, e_kin, e_th, div_v, vp, vm
      integer :: j

      n_changed = 0
      if (.not. allocated(adiabat_stage)) return
      do j = 1, N
         if (.not. (adiabat_stage(j) .gt. 0.0d0)) cycle
         rho_c = u(1,j)
         e_kin = 0.5d0*u(2,j)*u(2,j)/rho_c
         e_th  = u(3,j) - e_kin
         if (.not. (e_th .lt. thermal_fraction_for_adiabat*u(3,j))) cycle
         ! div v of the cell from the neighbouring cell velocities, on the
         ! cell's own faces: (r_+^2 v_+ - r_-^2 v_-) with face velocities
         ! the means of the two cells beside each face.
         vp = 0.5d0*(u(2,j)/u(1,j) + u(2,j+1)/u(1,j+1))
         vm = 0.5d0*(u(2,j)/u(1,j) + u(2,j-1)/u(1,j-1))
         div_v = r_edg(j)*r_edg(j)*vp - r_edg(j-1)*r_edg(j-1)*vm
         if (.not. (div_v .ge. 0.0d0)) cycle
         u(3,j) = e_kin + adiabat_stage(j)*rho_c**gamma_ad/(gamma_ad - 1.0d0)
         n_changed = n_changed + 1
      enddo
      n_cells_internal_energy_from_adiabat =                               &
         n_cells_internal_energy_from_adiabat + n_changed

      end subroutine internal_energy_from_adiabat

      end module adiabatic_entropy_transport
