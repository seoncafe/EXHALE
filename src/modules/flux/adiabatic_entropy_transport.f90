      module adiabatic_entropy_transport
      ! Dual-energy update for smooth, cold, supersonic expansion (Ryu et
      ! al. 1993, ApJ 414, 1; Bryan et al. 1995, Comput. Phys. Commun. 89,
      ! 149). Recovering rho*e = E-rho*v^2/2 magnifies truncation error when
      ! rho*e/(rho*v^2/2) = 2/[gamma*(gamma-1)*M^2] is small.
      ! Between split source updates, advect K=p/rho^gamma through
      ! d(rho*K)/dt + div(rho*K*v) = 0 with the mass-row face flux.
      ! Use K*rho^gamma/(gamma-1) for internal energy only where
      ! rho*e < eta*E, div(v) >= 0, and the EOS has constant gamma_ad;
      ! eta=thermal_fraction_for_adiabat=1e-2. Rebuild total energy, and
      ! reset K from the state each step so split sources remain included.
      ! Shocks and molecular EOS retain total-energy evolution. Valid for
      ! smooth expanding flow at roughly M > 13. The measured failure and
      ! correction are in ../EXHALE_v1.00/md/Update_EXHALE_stage3.md section 139.

      use global_parameters
      use caloric_eos, only: adiabatic_index_from_state
      use species_advective_transport, only: species_flux_divergence

      implicit none

      ! Thermal fraction of the total energy below which the internal energy
      ! is taken from the adiabat.
      real*8, parameter :: thermal_fraction_for_adiabat = 1.0d-2

      ! The adiabat K = p/rho^gamma at the start of the step and after the
      ! last stage. Zero marks a cell whose adiabat is not defined (no
      ! single gamma, or p <= 0).
      real*8, allocatable, dimension(:) :: adiabat_step, adiabat_stage

      ! Cells reset after admissible hydrodynamic stages, counted on each call.
      ! The outer attempted-step controller can still reject the state.
      integer :: n_cells_internal_energy_from_adiabat = 0

      contains


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


      subroutine internal_energy_from_adiabat(u, n_changed)
      ! After admissible hydrodynamic stages, reset internal energy from the
      ! advected adiabat under the module-header conditions. n_changed counts
      ! reset cells; final adoption belongs to the outer step controller.
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
