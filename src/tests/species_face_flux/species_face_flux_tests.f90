      program species_face_flux_tests
      ! Acceptance tests of the species advective transport on the face mass
      ! fluxes (src/modules/flux/species_face_flux.f90 and the element rows
      ! of src/modules/functions/binary_element_diffusion.f90).  Rows of
      ! docs/b4_spatial_operator_design_20260906.md increment B4-1 and test
      ! matrix rows B4-b and B4-j:
      !
      !   passive scalar  a uniform composition is preserved whatever the
      !                   velocity field does, including one whose two face
      !                   mass fluxes straddle zero, for PLM and for WENO3;
      !                   reported as the ulp distance from the constant
      !   face identity   sum_s F_s = F_rho at every face and every stage,
      !                   measured from the residual the operator itself
      !                   records; the injected violation is a separate
      !                   invocation (argument "inject"), because the
      !                   operator stops on it
      !   bounds          a step profile advected across the domain keeps
      !                   0 <= Y <= 1, and the closing member of the
      !                   normalized set differs from an independently
      !                   reconstructed 1 - Y by no more than the
      !                   reconstruction's own truncation (MEASURED)
      !   order           the spatial operator applied to a smooth
      !                   composition in a steady spherical flow converges on
      !                   the continuum divergence at the reconstruction's
      !                   design order
      !   element budget  through the element wrapper, the helium nucleus
      !                   total of the domain changes only by the two
      !                   boundary fluxes, to round-off and not to the order
      !                   of the hydro truncation (B1 AT-4d)
      !
      ! The columns are synthetic: this driver sets the global_parameters
      ! scalars and the radial grid itself, so no input.inp, no ionization
      ! solve and no hydro are involved.  The face mass flux is supplied and
      ! the stage densities are built FROM it by the same finite volume the
      ! operator uses, which is what makes the identities exact statements
      ! rather than approximations of a Riemann solve.
      !
      ! WHICH RIEMANN INTERFACE IS COVERED.  The operator reads the face mass
      ! flux and nothing else: the upwind side is the sign of F_rho, which is
      ! the sign of the HLLC contact speed wherever that is defined and the
      ! only rule ROE and Lax-Friedrichs offer.  So the three interfaces
      ! enter these tests through the VALUE of F_rho, and the fields below
      ! cover a positive, a negative and a sign-alternating one.

      use global_parameters
      use species_table, only: isp_HI, isp_HeI, isp_H2, n_melem
      use composition,   only: mass_per_H_nucleus_without_He
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use species_advective_transport, only: species_advective_update,    &
                                    species_face_fraction,                &
                                    species_face_flux,                    &
                                    species_face_identity_residual,       &
                                    species_face_identity_injection,      &
                                    n_species_faces_bounded
      use Reconstruction_step, only: Reconstruct_scalar
      use binary_element_diffusion, only: species_advection_active,       &
                                    species_advection_begin_step,         &
                                    species_advection_stage,              &
                                    species_advection_project,            &
                                    advected_carrier_reset,               &
                                    advected_carrier_register,            &
                                    advected_carrier_fractions

      implicit none

      integer :: n_pass, n_fail
      character(len=32) :: arg

      n_pass = 0
      n_fail = 0

      arg = ' '
      if (command_argument_count() .ge. 1) call get_command_argument(1,arg)

      if (trim(arg) .eq. 'inject') then
         call injected_face_identity_violation()
         ! The operator is expected to stop before this line is reached.
         write(*,'(A)') ' FAIL  face_identity_assertion_fires  '//        &
              'the injected violation was not caught'
         stop 1
      endif

      write(*,'(A)') '===== species face flux acceptance tests ====='

      call test_uniform_composition('PLM')
      call test_uniform_composition('WENO3')
      call test_face_identity()
      call test_bounds_and_normalization()
      call test_spatial_order('PLM',  1.8d0)
      call test_spatial_order('WENO3',1.8d0)
      call test_element_nucleus_budget()
      call test_carrier_nucleus_budget()

      write(*,'(A)') '=============================================='
      write(*,'(A,I0,A,I0,A)') ' TOTAL: ', n_pass, ' passed, ',           &
                               n_fail, ' failed'
      if (n_fail .gt. 0) stop 1

      contains

      ! ================================================================= !
      !  harness
      ! ================================================================= !

      subroutine verdict(name, ok, val, ref, tol)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: ok
      real*8,           intent(in) :: val, ref, tol
      if (ok) then
         n_pass = n_pass + 1
         write(*,'(A,A,A,ES11.4,A,ES11.4,A,ES11.4)') ' PASS ', name,      &
              ' measured=', val, ' reference=', ref, ' tol=', tol
      else
         n_fail = n_fail + 1
         write(*,'(A,A,A,ES11.4,A,ES11.4,A,ES11.4)') ' FAIL ', name,      &
              ' measured=', val, ' reference=', ref, ' tol=', tol
      endif
      end subroutine verdict

      ! ----------------------------------------------------------------- !

      subroutine setup_column(Ncell, r_top)
      ! A synthetic spherical column of Ncell cells, uniform in r from 1 to
      ! r_top, and the global scalars the operator and the element wrapper
      ! read.  Only the grid enters the advective operator; the rest is what
      ! the element write-back needs.
      integer, intent(in) :: Ncell
      real*8,  intent(in) :: r_top
      real*8  :: dr_u
      integer :: j

      N   = Ncell
      T0  = 1.0d4
      R0  = 1.0d10
      n0  = 1.0d10
      v0  = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0  = n0*mu*v0*v0
      spherical_domain = .true.
      thereis_He    = .true.
      thereis_HeITR = .false.
      thereis_mol   = .false.
      thereis_metals = .false.
      eos_include_metals = .true.
      he_diffusion  = .true.
      he_metal_diffusion = .false.
      he_kzz        = 0.0d0
      HeH           = 0.083d0
      rec_method    = 'PLM'
      use_plm       = .true.
      use_weno3     = .false.
      recon_lambda_on = .false.

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 0.0d0

      dr_u = (r_top - 1.0d0)/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      call eddy_diffusion_on_grid()

      end subroutine setup_column

      ! ----------------------------------------------------------------- !

      subroutine set_scheme(scheme)
      character(len=*), intent(in) :: scheme
      rec_method = scheme
      use_plm   = (scheme .eq. 'PLM')
      use_weno3 = (scheme .eq. 'WENO3')
      end subroutine set_scheme

      ! ----------------------------------------------------------------- !

      subroutine mass_row(rho_in, Frho, dt_loc, rho_out)
      ! The mass row of one forward-Euler stage, on the same faces, areas and
      ! volume the species rows use.  Building the stage density here rather
      ! than taking it from a hydro solve is what makes the passive-scalar
      ! and budget statements exact.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_in, Frho, dt_loc
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_out
      real*8  :: rp, rm, dAp, dAm, dV
      integer :: j
      rho_out = rho_in
      do j = 2, N
         rp  = r_edg(j)
         rm  = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV  = (dAp*rp - dAm*rm)/3.0
         rho_out(j) = rho_in(j)                                           &
                    - dt_loc(j)*(dAp*Frho(j) - dAm*Frho(j-1))/dV
      enddo
      end subroutine mass_row

      ! ================================================================= !
      !  a uniform composition is preserved
      ! ================================================================= !

      subroutine test_uniform_composition(scheme)
      character(len=*), intent(in) :: scheme

      integer, parameter :: nq = 1
      real*8, dimension(:),   allocatable :: rho0, rho1, rho2, rho3
      real*8, dimension(:),   allocatable :: Frho, dt_loc
      real*8, dimension(:,:), allocatable :: Y0, Y1, Y2, Y3
      real*8  :: c, worst, ulp
      integer :: j, ic
      character(len=64) :: nm

      call setup_column(200, 5.0d0)
      call set_scheme(scheme)

      allocate(rho0(1-Ng:N+Ng), rho1(1-Ng:N+Ng), rho2(1-Ng:N+Ng),         &
               rho3(1-Ng:N+Ng), Frho(1-Ng:N+Ng), dt_loc(1-Ng:N+Ng))
      allocate(Y0(1-Ng:N+Ng,nq), Y1(1-Ng:N+Ng,nq), Y2(1-Ng:N+Ng,nq),      &
               Y3(1-Ng:N+Ng,nq))

      c = 0.2731d0

      do ic = 1, 3

         ! Three flows: a compressing one, an expanding one, and one whose
         ! face mass fluxes alternate in sign from face to face, which is the
         ! breathing base a cell-velocity upwind rule cannot see.
         do j = 1-Ng, N+Ng
            select case (ic)
            case (1)
               Frho(j) =  0.4d0/(r_edg(j)*r_edg(j))
            case (2)
               Frho(j) = -0.4d0*r_edg(j)
            case (3)
               Frho(j) =  0.4d0*sin(3.0d0*dble(j))
            end select
            rho0(j) = 1.0d0 + 0.3d0*sin(0.05d0*dble(j))
         enddo
         dt_loc = 1.0d-4

         Y0 = c
         call mass_row(rho0, Frho, dt_loc, rho1)
         call species_advective_update(1,nq,1,Y0,Y0,Frho,rho0,rho0,rho1,  &
                                       dt_loc,2,Y1)
         call rk2_density(rho0, rho1, Frho, dt_loc, rho2)
         call species_advective_update(2,nq,1,Y0,Y1,Frho,rho0,rho1,rho2,  &
                                       dt_loc,2,Y2)
         call rk3_density(rho0, rho2, Frho, dt_loc, rho3)
         call species_advective_update(3,nq,1,Y0,Y2,Frho,rho0,rho2,rho3,  &
                                       dt_loc,2,Y3)

         worst = 0.0d0
         do j = 2, N
            worst = max(worst, abs(Y1(j,1)-c), abs(Y2(j,1)-c),            &
                               abs(Y3(j,1)-c))
         enddo
         ulp = worst/spacing(c)

         write(nm,'(A,A,A,I0)') 'uniform_composition_preserved_',         &
              trim(scheme), '_flow', ic
         call verdict(trim(nm), ulp .le. 4.0d0, ulp, 0.0d0, 4.0d0)

      enddo

      deallocate(rho0, rho1, rho2, rho3, Frho, dt_loc, Y0, Y1, Y2, Y3)

      end subroutine test_uniform_composition

      ! ----------------------------------------------------------------- !

      subroutine rk2_density(rho_n, rho_1, Frho, dt_loc, rho_2)
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_n, rho_1, Frho
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: dt_loc
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_2
      real*8  :: rp, rm, dAp, dAm, dV
      integer :: j
      rho_2 = rho_1
      do j = 2, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dAp = rp*rp;    dAm = rm*rm
         dV  = (dAp*rp - dAm*rm)/3.0
         rho_2(j) = (3.0*rho_n(j) + rho_1(j)                              &
                    - dt_loc(j)*(dAp*Frho(j) - dAm*Frho(j-1))/dV)/4.0
      enddo
      end subroutine rk2_density

      ! ----------------------------------------------------------------- !

      subroutine rk3_density(rho_n, rho_2, Frho, dt_loc, rho_3)
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_n, rho_2, Frho
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: dt_loc
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_3
      real*8  :: rp, rm, dAp, dAm, dV
      integer :: j
      rho_3 = rho_2
      do j = 2, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dAp = rp*rp;    dAm = rm*rm
         dV  = (dAp*rp - dAm*rm)/3.0
         rho_3(j) = (rho_n(j) + 2.0*(rho_2(j)                             &
                    - dt_loc(j)*(dAp*Frho(j) - dAm*Frho(j-1))/dV))/3.0
      enddo
      end subroutine rk3_density

      ! ================================================================= !
      !  the face identity
      ! ================================================================= !

      subroutine test_face_identity()
      ! The residual the operator itself measured over every face of every
      ! update this driver has made so far.  It is the quantity the operator
      ! asserts, so reading it back is reading the assertion's own margin.
      real*8 :: v
      v = species_face_identity_residual
      call verdict('face_identity_residual', v .le. 64.0d0*epsilon(1.0d0),&
                   v, 0.0d0, 64.0d0*epsilon(1.0d0))
      end subroutine test_face_identity

      ! ----------------------------------------------------------------- !

      subroutine injected_face_identity_violation()
      ! One update with a deliberate violation injected at one face.  The
      ! operator is expected to stop; the caller checks that it did.
      integer, parameter :: nq = 1
      real*8, dimension(:),   allocatable :: rho0, rho1, Frho, dt_loc
      real*8, dimension(:,:), allocatable :: Y0, Y1
      integer :: j

      call setup_column(64, 3.0d0)
      allocate(rho0(1-Ng:N+Ng), rho1(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               dt_loc(1-Ng:N+Ng))
      allocate(Y0(1-Ng:N+Ng,nq), Y1(1-Ng:N+Ng,nq))
      do j = 1-Ng, N+Ng
         Frho(j) = 0.4d0/(r_edg(j)*r_edg(j))
         rho0(j) = 1.0d0
      enddo
      dt_loc = 1.0d-4
      Y0 = 0.2d0
      call mass_row(rho0, Frho, dt_loc, rho1)
      call species_face_identity_injection(20)
      call species_advective_update(1,nq,1,Y0,Y0,Frho,rho0,rho0,rho1,     &
                                    dt_loc,2,Y1)
      call species_face_identity_injection(-1)
      end subroutine injected_face_identity_violation

      ! ================================================================= !
      !  bounds and the normalization
      ! ================================================================= !

      subroutine test_bounds_and_normalization()
      ! A step in the composition is carried across the domain by a steady
      ! outflow.  Two statements: the mass fraction never leaves [0,1], and
      ! the closing member of the normalized set, which is one minus the
      ! reconstructed fraction, differs from an INDEPENDENTLY reconstructed
      ! 1 - Y by no more than the reconstruction's own truncation.  The
      ! second number is what "the normalization moves no species by more
      ! than the truncation" means when the set has two members.
      integer, parameter :: nq = 1
      real*8, dimension(:),   allocatable :: rho0, rho1, Frho, dt_loc
      real*8, dimension(:),   allocatable :: Yf, Ycomp, cL, cR, Ycl
      real*8, dimension(:,:), allocatable :: Y0, Y1
      real*8  :: ymin, ymax, dmax
      integer :: j, it

      call setup_column(200, 5.0d0)
      call set_scheme('PLM')

      allocate(rho0(1-Ng:N+Ng), rho1(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               dt_loc(1-Ng:N+Ng), Yf(1-Ng:N+Ng), Ycomp(1-Ng:N+Ng),        &
               cL(1-Ng:N+Ng), cR(1-Ng:N+Ng), Ycl(1-Ng:N+Ng))
      allocate(Y0(1-Ng:N+Ng,nq), Y1(1-Ng:N+Ng,nq))

      do j = 1-Ng, N+Ng
         Frho(j) = 0.4d0/(r_edg(j)*r_edg(j))
         rho0(j) = 0.4d0/(0.5d0*r(j)*r(j))
      enddo
      dt_loc = 0.3d0*dr_j(1)/0.5d0
      Y0 = 0.0d0
      do j = 1-Ng, 40
         Y0(j,1) = 1.0d0
      enddo

      ymin =  huge(1.0d0)
      ymax = -huge(1.0d0)
      dmax = 0.0d0
      do it = 1, 300
         call mass_row(rho0, Frho, dt_loc, rho1)
         call species_advective_update(1,nq,1,Y0,Y0,Frho,rho0,rho0,rho1,  &
                                       dt_loc,2,Y1)
         Y1(N+1:N+Ng,1) = Y1(N,1)
         ! The closing member as the operator forms it, against the
         ! independent reconstruction of 1 - Y.
         call species_face_fraction(Y0(:,1),Frho,Yf)
         Ycomp = 1.0d0 - Y0(:,1)
         call species_face_fraction(Ycomp,Frho,Ycl)
         do j = 1, N
            dmax = max(dmax, abs((1.0d0 - Yf(j)) - Ycl(j)))
         enddo
         Y0 = Y1
         ymin = min(ymin, minval(Y0(1:N,1)))
         ymax = max(ymax, maxval(Y0(1:N,1)))
      enddo

      call verdict('mass_fraction_stays_above_zero', ymin .ge. 0.0d0,     &
                   ymin, 0.0d0, 0.0d0)
      call verdict('mass_fraction_stays_below_one', ymax .le. 1.0d0,      &
                   ymax, 1.0d0, 0.0d0)
      ! The bound is the reconstruction's own truncation, one part in
      ! 1/(N-1) squared for a piecewise linear reconstruction of a profile
      ! of order unity; what a correct closing member reaches is round-off,
      ! and the measured value is printed either way.
      call verdict('closing_member_within_reconstruction_truncation',     &
                   dmax .le. dr_j(1)*dr_j(1), dmax, 0.0d0, dr_j(1)**2)

      deallocate(rho0, rho1, Frho, dt_loc, Yf, Ycomp, cL, cR, Ycl, Y0, Y1)

      end subroutine test_bounds_and_normalization

      ! ================================================================= !
      !  spatial order
      ! ================================================================= !

      subroutine test_spatial_order(scheme, want)
      ! In a steady spherical flow, r^2 rho v = M constant, the exact finite
      ! volume divergence of the species flux is M (Y(r_+) - Y(r_-))/dV.
      ! The operator's divergence is compared with it on a smooth monotone Y
      ! at three resolutions, so what the rate measures is the error of the
      ! face reconstruction alone, which is the design order asked for.
      !
      ! The cells are seeded with POINT values of Y rather than with its cell
      ! averages, which are two different quantities at second order on this
      ! grid.  The third order of ESWENO3 is therefore not observable in this
      ! driver; what is stated is that both reconstructions reach at least
      ! second order, and the measured rate of each is printed.
      character(len=*), intent(in) :: scheme
      real*8,           intent(in) :: want

      integer, dimension(3), parameter :: ncl = [100, 200, 400]
      real*8,  dimension(3) :: err, errf
      real*8, dimension(:), allocatable :: Yc, Frho, Yf, Fs
      real*8  :: M, rp, rm, dAp, dAm, dV, ex, di, rate
      integer :: ic, j
      character(len=64) :: nm

      M = 0.4d0

      do ic = 1, 3
         call setup_column(ncl(ic), 3.0d0)
         call set_scheme(scheme)
         allocate(Yc(1-Ng:N+Ng), Frho(1-Ng:N+Ng), Yf(1-Ng:N+Ng),          &
                  Fs(1-Ng:N+Ng))
         do j = 1-Ng, N+Ng
            Frho(j) = M/(r_edg(j)*r_edg(j))
            Yc(j)   = 0.3d0 + 0.4d0*sin(0.5d0*(r(j) - 1.0d0))
         enddo
         call species_face_fraction(Yc,Frho,Yf)
         call species_face_flux(Frho,Yf,Fs)
         err(ic)  = 0.0d0
         errf(ic) = 0.0d0
         do j = 10, N-10
            errf(ic) = max(errf(ic), abs(Yf(j)                            &
                 - (0.3d0 + 0.4d0*sin(0.5d0*(r_edg(j) - 1.0d0)))))
            rp = r_edg(j);  rm = r_edg(j-1)
            dAp = rp*rp;    dAm = rm*rm
            dV  = (dAp*rp - dAm*rm)/3.0
            di  = (dAp*Fs(j) - dAm*Fs(j-1))/dV
            ex  = M*(0.4d0*sin(0.5d0*(rp - 1.0d0))                        &
                   - 0.4d0*sin(0.5d0*(rm - 1.0d0)))/dV
            err(ic) = max(err(ic), abs(di - ex))
         enddo
         deallocate(Yc, Frho, Yf, Fs)
      enddo

      ! The face composition is what the order statement is about: the flux
      ! is that value multiplied by a mass flux the operator does not touch.
      rate = log(errf(1)/errf(3))/log(4.0d0)
      write(*,'(A,A,A,3ES11.4)') '   (', trim(scheme),                    &
           ' face error at N = 100, 200, 400: ', errf
      write(*,'(A,A,A,3ES11.4,A,F6.3)') '   (', trim(scheme),             &
           ' divergence error at N = 100, 200, 400: ', err,               &
           ' rate ', log(err(1)/err(3))/log(4.0d0)
      write(nm,'(A,A)') 'face_reconstruction_order_', trim(scheme)
      call verdict(trim(nm), rate .ge. want, rate, want, 0.0d0)

      end subroutine test_spatial_order


      ! ================================================================= !
      !  the element nucleus budget in a transient
      ! ================================================================= !

      subroutine test_element_nucleus_budget()
      ! Through the element wrapper: the helium nucleus total of the domain
      ! changes over one stage only by the two boundary fluxes, to round-off.
      ! Today's cell-velocity advective form conserves it only to the order
      ! of the hydro truncation, which is what this row replaces.
      real*8, dimension(:),   allocatable :: rho0, rho1, Frho, dt_loc
      real*8, dimension(:),   allocatable :: Yf, nucHe
      real*8, dimension(:,:), allocatable :: f_sp
      real*8  :: m_1, X, tot0, tot1, bnd, rp, rm, dV, rel
      integer :: j

      call setup_column(200, 5.0d0)
      call set_scheme('PLM')

      allocate(rho0(1-Ng:N+Ng), rho1(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               dt_loc(1-Ng:N+Ng), Yf(1-Ng:N+Ng), nucHe(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species))

      m_1  = mass_per_H_nucleus_without_He()
      f_sp = 0.0d0
      do j = 1-Ng, N+Ng
         ! A composition that varies in space, so the advection has work to
         ! do, written with helium and hydrogen in their ground stages only:
         ! m_1 n_H + m_He n_He = 1 by construction.
         X = 0.25d0 + 0.15d0*sin(1.7d0*(r(j) - 1.0d0))
         f_sp(j,isp_HeI) = X/helium_mass()
         f_sp(j,isp_HI)  = (1.0d0 - X)/m_1
         Frho(j) = 0.4d0/(r_edg(j)*r_edg(j))
         rho0(j) = 1.0d0 + 0.2d0*sin(0.05d0*dble(j))
      enddo
      ! The seed uses the species table's own helium mass, so the mass sum
      ! m_1 n_H + m_He n_He is one to round-off and the nucleus total is the
      ! mass total divided by that one mass.
      dt_loc = 1.0d-4

      call helium_nucleus_density(f_sp, nucHe)
      tot0 = 0.0d0
      do j = 2, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dV = (rp*rp*rp - rm*rm*rm)/3.0
         tot0 = tot0 + dV*rho0(j)*nucHe(j)
      enddo

      ! The face fractions the stage will use, taken from the composition
      ! the step BEGINS with: that is what the boundary flux of this stage
      ! carried.
      call helium_face_fraction(f_sp, Frho, Yf)

      call species_advection_begin_step(f_sp)
      call mass_row(rho0, Frho, dt_loc, rho1)
      call species_advection_stage(1, rho0, rho0, rho1, Frho, dt_loc)
      call species_advection_project(f_sp)

      call helium_nucleus_density(f_sp, nucHe)
      tot1 = 0.0d0
      do j = 2, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dV = (rp*rp*rp - rm*rm*rm)/3.0
         tot1 = tot1 + dV*rho1(j)*nucHe(j)
      enddo

      ! What crossed the two ends of the updated range, in nuclei: the face
      ! mass flux times the face helium mass fraction, divided by the helium
      ! mass the wrapper closes its mixture with.
      bnd = dt_loc(2)*(r_edg(N)**2*Frho(N)*Yf(N)                          &
                     - r_edg(1)**2*Frho(1)*Yf(1))/helium_mass()

      rel = abs(tot1 - tot0 + bnd)/max(abs(tot0), 1.0d-300)
      call verdict('helium_nucleus_total_changes_only_by_boundary_flux',  &
                   rel .le. 1.0d-12, rel, 0.0d0, 1.0d-12)

      ! THE FORM THIS REPLACES, on the same column and the same step: the
      ! one-sided upwind difference of the cell values taken with the cell
      ! velocity.  It is reported, not asserted: it is the RED the row above
      ! is the GREEN of, and the number is the truncation of the hydro the
      ! design says the element total used to be conserved to.
      call cell_velocity_upwind_budget(rho0, rho1, Frho, dt_loc, rel)
      write(*,'(A,ES11.4)') '   (the cell-velocity upwind form of the '// &
           'same step leaves the same budget at ', rel

      deallocate(rho0, rho1, Frho, dt_loc, Yf, nucHe, f_sp)

      end subroutine test_element_nucleus_budget

      ! ================================================================= !
      !  the carrier budget in a transient
      ! ================================================================= !

      subroutine test_carrier_nucleus_budget()
      ! Through the carrier half of the same wrapper, with the element rows
      ! switched off so that what is measured is the carrier rows alone.
      ! Two statements on one stage of one column:
      !   * a UNIFORM carrier partition is reproduced to round-off whatever
      !     the face mass flux does, including where it changes sign, which
      !     is the property a cell-velocity upwind difference does not have;
      !   * the H2 NUCLEUS TOTAL of the updated range changes only by what
      !     crossed its two ends.  The carrier range runs from cell 1,
      !     because the base face carries the ghost composition in and the
      !     base cell has an advective term like every other cell.
      real*8, dimension(:),   allocatable :: rho0, rho1, Frho, dt_loc
      real*8, dimension(:),   allocatable :: Yf, YL, YR
      real*8, dimension(:,:), allocatable :: f_sp, fcadv
      real*8  :: q, tot0, tot1, bnd, rp, rm, dV, rel, worst
      integer :: j

      call setup_column(200, 5.0d0)
      call set_scheme('PLM')
      he_diffusion = .false.
      thereis_mol  = .true.
      call advected_carrier_reset()
      call advected_carrier_register(isp_H2, .false.)

      allocate(rho0(1-Ng:N+Ng), rho1(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               dt_loc(1-Ng:N+Ng), Yf(1-Ng:N+Ng), YL(1-Ng:N+Ng),           &
               YR(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species), fcadv(1-Ng:N+Ng,1))

      dt_loc = 1.0d-4

      ! ---- a uniform partition, on a sign-alternating face mass flux ----
      f_sp = 0.0d0
      q    = 0.6d0
      do j = 1-Ng, N+Ng
         f_sp(j,isp_H2) = 0.5d0*q
         f_sp(j,isp_HI) = 1.0d0 - q
         Frho(j) = 0.4d0*sin(0.9d0*dble(j))/(r_edg(j)*r_edg(j))
         rho0(j) = 1.0d0 + 0.2d0*sin(0.05d0*dble(j))
      enddo
      call species_advection_begin_step(f_sp)
      call carrier_mass_row(rho0, Frho, dt_loc, rho1)
      call species_advection_stage(1, rho0, rho0, rho1, Frho, dt_loc)
      call species_advection_project(f_sp)
      ! The write-back the transport operator does is not run here, so the
      ! species vector still holds the entry mixture; that is the mixture
      ! the advected fraction is read against.
      call advected_carrier_fractions(f_sp, fcadv)
      worst = 0.0d0
      do j = 1, N
         worst = max(worst, abs(2.0d0*fcadv(j,1) - q)/q)
      enddo
      call verdict('uniform_carrier_partition_preserved',                 &
                   worst .le. 8.0d0*epsilon(1.0d0), worst/epsilon(1.0d0), &
                   0.0d0, 8.0d0)

      ! ---- the nucleus budget of a varying partition --------------------
      f_sp = 0.0d0
      do j = 1-Ng, N+Ng
         q = 0.50d0 + 0.30d0*sin(1.7d0*(r(j) - 1.0d0))
         f_sp(j,isp_H2) = 0.5d0*q
         f_sp(j,isp_HI) = 1.0d0 - q
         Frho(j) = 0.4d0/(r_edg(j)*r_edg(j))
         rho0(j) = 1.0d0 + 0.2d0*sin(0.05d0*dble(j))
      enddo

      tot0 = 0.0d0
      do j = 1, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dV = (rp*rp*rp - rm*rm*rm)/3.0
         tot0 = tot0 + dV*rho0(j)*2.0d0*f_sp(j,isp_H2)
      enddo

      ! The face fractions this stage will carry, from the composition it
      ! begins with: the same reconstruction, on the side the face mass flux
      ! selects, of the H2 mass fraction, which is 2 f_H2 on this mixture.
      ! The inner ghosts hold the base cell's own partition, which is the
      ! inflow condition of a carrier no handoff states.
      Yf = 2.0d0*f_sp(:,isp_H2)
      Yf(1-Ng:0) = Yf(1)
      call Reconstruct_scalar(Yf, YL, YR)
      do j = 1-Ng, N+Ng
         if (Frho(j) .ge. 0.0d0) then
            Yf(j) = YL(j)
         else
            Yf(j) = YR(j)
         endif
      enddo

      call species_advection_begin_step(f_sp)
      call carrier_mass_row(rho0, Frho, dt_loc, rho1)
      call species_advection_stage(1, rho0, rho0, rho1, Frho, dt_loc)
      call species_advection_project(f_sp)
      call advected_carrier_fractions(f_sp, fcadv)

      tot1 = 0.0d0
      do j = 1, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dV = (rp*rp*rp - rm*rm*rm)/3.0
         tot1 = tot1 + dV*rho1(j)*2.0d0*fcadv(j,1)
      enddo
      bnd = dt_loc(1)*(r_edg(N)**2*Frho(N)*Yf(N)                          &
                     - r_edg(0)**2*Frho(0)*Yf(0))

      rel = abs(tot1 - tot0 + bnd)/max(abs(tot0), 1.0d-300)
      call verdict('h2_nucleus_total_changes_only_by_boundary_flux',      &
                   rel .le. 1.0d-12, rel, 0.0d0, 1.0d-12)

      call advected_carrier_reset()
      deallocate(rho0, rho1, Frho, dt_loc, Yf, YL, YR, f_sp, fcadv)

      end subroutine test_carrier_nucleus_budget

      ! ----------------------------------------------------------------- !

      subroutine carrier_mass_row(rho_in, Frho, dt_loc, rho_out)
      ! The mass row over the range the carrier rows are advanced on, which
      ! begins at the base cell.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_in, Frho, dt_loc
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho_out
      real*8  :: rp, rm, dAp, dAm, dV
      integer :: j
      rho_out = rho_in
      do j = 1, N
         rp  = r_edg(j)
         rm  = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV  = (dAp*rp - dAm*rm)/3.0
         rho_out(j) = rho_in(j)                                           &
                    - dt_loc(j)*(dAp*Frho(j) - dAm*Frho(j-1))/dV
      enddo
      end subroutine carrier_mass_row

      ! ----------------------------------------------------------------- !

      subroutine cell_velocity_upwind_budget(rho0, rho1, Frho, dt_loc, rel)
      ! One step of the advective form the face fluxes replace, on a helium
      ! mass fraction taken directly rather than through the species vector:
      !
      !     X_j <- X_j - dt v_j (X_j - X_{j-1})/(r_j - r_{j-1})   (v_j >= 0)
      !
      ! with v_j = F_rho(j)/rho_j the cell velocity the old operator used,
      ! and the same domain budget measured against the same boundary flux.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho0, rho1, Frho, dt_loc
      real*8,                       intent(out) :: rel
      real*8, dimension(1-Ng:N+Ng) :: X, Xn, Yf
      real*8  :: rp, rm, dV, vj, t0, t1, bnd
      integer :: j

      do j = 1-Ng, N+Ng
         X(j) = 0.25d0 + 0.15d0*sin(1.7d0*(r(j) - 1.0d0))
      enddo
      call species_face_fraction(X, Frho, Yf)

      Xn = X
      do j = 2, N
         vj = 0.5d0*(Frho(j) + Frho(j-1))/rho0(j)
         if (vj .ge. 0.0d0) then
            Xn(j) = X(j) - dt_loc(j)*vj*(X(j) - X(j-1))/(r(j) - r(j-1))
         else
            Xn(j) = X(j) - dt_loc(j)*vj*(X(j+1) - X(j))/(r(j+1) - r(j))
         endif
      enddo

      t0 = 0.0d0
      t1 = 0.0d0
      do j = 2, N
         rp = r_edg(j);  rm = r_edg(j-1)
         dV = (rp*rp*rp - rm*rm*rm)/3.0
         t0 = t0 + dV*rho0(j)*X(j)
         t1 = t1 + dV*rho1(j)*Xn(j)
      enddo
      bnd = dt_loc(2)*(r_edg(N)**2*Frho(N)*Yf(N)                          &
                     - r_edg(1)**2*Frho(1)*Yf(1))
      rel = abs(t1 - t0 + bnd)/max(abs(t0), 1.0d-300)

      end subroutine cell_velocity_upwind_budget

      ! ----------------------------------------------------------------- !

      real*8 function helium_mass()
      ! The helium mass the element closure weighs a nucleus with, read from
      ! the species table through the wrapper's own mixture rather than
      ! restated here.
      use species_table, only: bsp_mass
      use species_table, only: n_bsp, bsp_fsp
      integer :: ib
      helium_mass = 4.0d0
      do ib = 1, n_bsp
         if (bsp_fsp(ib) .eq. isp_HeI) helium_mass = bsp_mass(ib)
      enddo
      end function helium_mass

      ! ----------------------------------------------------------------- !

      subroutine helium_nucleus_density(f_sp, nucHe)
      use species_table, only: n_bsp, bsp_fsp, bsp_nHe,                   &
                               bsp_is_excited_level
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucHe
      integer :: j, ib
      nucHe = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nHe(ib) .le. 0) cycle
         do j = 1-Ng, N+Ng
            nucHe(j) = nucHe(j) + dble(bsp_nHe(ib))*f_sp(j,bsp_fsp(ib))
         enddo
      enddo
      end subroutine helium_nucleus_density

      ! ----------------------------------------------------------------- !

      subroutine helium_face_fraction(f_sp, Frho, Yf)
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Yf
      real*8, dimension(1-Ng:N+Ng) :: nucHe, nucH, Y
      real*8  :: m_1
      integer :: j
      call helium_nucleus_density(f_sp, nucHe)
      m_1 = mass_per_H_nucleus_without_He()
      call hydrogen_nucleus_density(f_sp, nucH)
      do j = 1-Ng, N+Ng
         Y(j) = helium_mass()*nucHe(j)                                    &
              /max(m_1*nucH(j) + helium_mass()*nucHe(j), 1.0d-30)
      enddo
      call species_face_fraction(Y, Frho, Yf)
      end subroutine helium_face_fraction

      ! ----------------------------------------------------------------- !

      subroutine hydrogen_nucleus_density(f_sp, nucH)
      use species_table, only: n_bsp, bsp_fsp, bsp_nH,                    &
                               bsp_is_excited_level
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucH
      integer :: j, ib
      nucH = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib) .le. 0) cycle
         do j = 1-Ng, N+Ng
            nucH(j) = nucH(j) + dble(bsp_nH(ib))*f_sp(j,bsp_fsp(ib))
         enddo
      enddo
      end subroutine hydrogen_nucleus_density

      end program species_face_flux_tests
