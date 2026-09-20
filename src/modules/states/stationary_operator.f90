      module stationary_operator
      ! THE SPATIAL OPERATOR A STATIONARY STATE IS A STATE OF, AND THE MASS
      ! FLUX IT PUTS THROUGH EVERY FACE.
      !
      ! A stationary state of this code is a state the WENO3 residual
      ! vanishes on: every route that measures or solves a steady state
      ! selects that reconstruction before it evaluates anything, because a
      ! PLM+WENO3 run marches its first stage in PLM and a state that solves
      ! one discretization is not required to solve another.  The selection
      ! is made here, once, so that a measurement of a stationary state and
      ! the solve that produced it are measurements of the same operator.
      !
      ! THE FACE MASS FLUX.  In spherical symmetry the mass equation is
      !
      !     d(rho)/dt + (1/r^2) d(r^2 rho v)/dr = 0 ,
      !
      ! so the quantity the discrete equation transports across face f is
      ! r_f^2 (rho v)_f, with (rho v)_f the first component of the Riemann
      ! flux at that face.  On a stationary state it is the same number at
      ! every face and equal to the wind's own mass flux, which is what
      ! stationary_face_mass_flux returns and what the certification report
      ! prints beside the hydrodynamic rows.
      !
      ! WHY THE ENTRY POINT INSTALLS THE STATE ITSELF.  Three quantities the
      ! flux depends on are not arguments of the reconstruction: the
      ! composition and the caloric state (the pressure map of a molecular
      ! cell reads them), the ghost cells (the boundary states them), and the
      ! reconstruction selection.  A caller that installs two of the three
      ! measures an operator no route runs.  In particular init applies the
      ! boundary to the conserved array and does not copy the refreshed
      ! ghosts into the primitive array it returns, so a caller that converts
      ! that primitive array back overwrites them with the file's.
      use global_parameters
      use species_table,       only: n_mion
      use composition,         only: get_species_densities
      use Conversion,          only: U_to_W
      use BC_Apply,            only: Apply_BC
      use Reconstruction_step, only: Reconstruct
      use RK_integration,      only: RK_rhs, face_flux
      use base_boundary,       only: wind_window_mass_flux

      implicit none

      ! The reconstruction selection a caller held before it asked for the
      ! stationary one, so that a diagnostic leaves the run's configuration
      ! as it found it.
      type :: reconstruction_selection
         character(len=:), allocatable :: method
         logical :: weno3        = .false.
         logical :: plm          = .false.
         logical :: lambda_blend = .false.
      end type reconstruction_selection

      ! r_f^2 (rho v)_f over the faces of one evaluation, with the identity
      ! of the operator that produced it.  min and max are taken over the
      ! faces 0..N, the faces of the interior cells and the base face; base
      ! is face 0.  window_mean is the mean of the CELL-CENTERED product
      ! rho v r^2 over the wind window r >= r_flux (wind_window_mass_flux),
      ! the flux the wind carries; a face flux equals it only to the
      ! difference between a face value and a cell mean, which on a resolved
      ! wind is O(dr^2) and is not a mass leak.
      type :: face_mass_flux_budget
         logical           :: available        = .false.
         logical           :: window_available = .false.
         real*8            :: window_mean      = 0.0d0
         real*8            :: flux_min         = 0.0d0
         real*8            :: flux_max         = 0.0d0
         real*8            :: flux_base        = 0.0d0
         integer           :: j_min            = 0
         integer           :: j_max            = 0
         character(len=32) :: operator_name    = ' '
         character(len=32) :: flux_name        = ' '
         logical           :: well_balanced    = .false.
      end type face_mass_flux_budget

      contains

      ! ------------------------------------------------------!

      subroutine select_stationary_reconstruction(previous)
      ! Install the reconstruction every stationary route evaluates and
      ! solves under: third-order WENO, with no PLM -> WENO3 blend left
      ! armed (a blended right-hand side is neither operator, and
      ! reconstruction_operator_label names it PLM+WENO3(lambda=...)).
      ! `previous` returns what the caller held, for
      ! restore_reconstruction_selection.
      type(reconstruction_selection), intent(out), optional :: previous

      if (present(previous)) then
         previous%method       = rec_method
         previous%weno3        = use_weno3
         previous%plm          = use_plm
         previous%lambda_blend = recon_lambda_on
      endif

      rec_method      = 'WENO3'
      use_weno3       = .true.
      use_plm         = .false.
      recon_lambda_on = .false.

      end subroutine select_stationary_reconstruction

      ! ------------------------------------------------------!

      subroutine restore_reconstruction_selection(previous)
      ! Put back the selection select_stationary_reconstruction returned.
      type(reconstruction_selection), intent(in) :: previous

      if (allocated(previous%method)) rec_method = previous%method
      use_weno3       = previous%weno3
      use_plm         = previous%plm
      recon_lambda_on = previous%lambda_blend

      end subroutine restore_reconstruction_selection

      ! ------------------------------------------------------!

      subroutine stationary_face_mass_flux(u_state, f_sp_state, budget,   &
                                           r2_rho_v)
      ! r_f^2 (rho v)_f on every face of the state (u_state, f_sp_state),
      ! under the operator a stationary state is a state of.
      !
      ! The state is installed in the order each consumer needs it: the
      ! composition and the caloric state first, because the pressure map
      ! and the boundary read them; the boundary second, so the ghosts the
      ! reconstruction sees are the boundary's; then the reconstruction and
      ! the Riemann solve, which leave the face fluxes in face_flux.
      !
      ! r2_rho_v, when asked for, carries the faces 0 to N, face j being the
      ! one at r_edg(j); the entries outside that range are returned zero,
      ! since no face of the domain lies there.
      !
      ! The conserved array of the caller is not written to: the boundary is
      ! applied to a copy, so the ghosts of the caller's state are its own.
      ! The composition-derived module state (the caloric mixture and the
      ! cell-1 particle count the base boundary reads) IS this state's after
      ! the call: that is what "install the state" means.
      real*8, dimension(3,1-Ng:N+Ng),          intent(in)  :: u_state
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)  :: f_sp_state
      type(face_mass_flux_budget),             intent(out) :: budget
      real*8, dimension(1-Ng:N+Ng), optional,  intent(out) :: r2_rho_v

      real*8, dimension(3,1-Ng:N+Ng) :: u_work, W_work, WL, WR, dFq, Sq
      real*8, dimension(1-Ng:N+Ng)   :: rho_s, nhi_s, nhii_s, nhei_s
      real*8, dimension(1-Ng:N+Ng)   :: nheii_s, nheiii_s, nheiTR_s
      real*8, dimension(1-Ng:N+Ng)   :: ne_s, ntot_s, Fface
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm_s
      real*8  :: d_window
      logical :: have_F
      integer :: j
      type(reconstruction_selection) :: entry_selection

      nhei_s   = 0.0d0
      nheii_s  = 0.0d0
      nheiii_s = 0.0d0
      nheiTR_s = 0.0d0

      u_work = u_state
      rho_s  = u_work(1,:)
      call get_species_densities(rho_s, f_sp_state, nhi_s, nhii_s,        &
                                 nhei_s, nheii_s, nheiii_s, nheiTR_s,     &
                                 nm_s, ne_s, ntot_s)
      call Apply_BC(u_work)
      call U_to_W(u_work, W_work)
      call wind_window_mass_flux(W_work, budget%window_mean, d_window,    &
                                 have_F)

      call select_stationary_reconstruction(entry_selection)
      call Reconstruct(u_work, WL, WR)
      call RK_rhs(u_work, WL, WR, dFq, Sq)
      budget%operator_name = reconstruction_operator_label()
      budget%flux_name     = flux
      budget%well_balanced = well_balanced
      call restore_reconstruction_selection(entry_selection)

      Fface = 0.0d0
      do j = 0, N
         Fface(j) = face_flux(1,j)*r_edg(j)*r_edg(j)
      enddo
      budget%flux_base = Fface(0)
      budget%j_min     = 0
      budget%j_max     = 0
      budget%flux_min  = Fface(0)
      budget%flux_max  = Fface(0)
      do j = 0, N
         if (Fface(j) .lt. budget%flux_min) then
            budget%flux_min = Fface(j);  budget%j_min = j
         endif
         if (Fface(j) .gt. budget%flux_max) then
            budget%flux_max = Fface(j);  budget%j_max = j
         endif
      enddo
      budget%window_available = have_F
      budget%available        = .true.
      if (present(r2_rho_v)) r2_rho_v = Fface

      end subroutine stationary_face_mass_flux

      ! End of module
      end module stationary_operator
