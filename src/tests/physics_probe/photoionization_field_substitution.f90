      program photoionization_field_substitution
      ! QUANTITY UNDER TEST
      !   the ORDER in which the photoionization field and the composition
      !   that attenuates it are brought into agreement, in the equilibrium
      !   sweep (ionization_equilibrium.f90, xuv_block_H / xuv_block_HHe).
      !
      !   A cell's photoionization rate is set by the column of the cells
      !   OUTSIDE it and by its own, and by nothing inside it.  The sweep
      !   runs outside in, so it can be given the columns of the composition
      !   it has ALREADY solved (a substitution) or the columns of one whole
      !   composition fixed for the traversal (a lagged, Jacobi map).  Both
      !   are iterated to the same statement -- composition in equilibrium
      !   with the field its own opacity makes -- and three things are
      !   asserted about them.
      !
      !     (1) THE COLUMNS ARE THE SAME ONES.  advance_starward_columns,
      !         run over the whole grid in one block, reproduces
      !         calc_column_dens BIT FOR BIT: same rectangle rule, same
      !         opa_pf weight, same order of accumulation.  Without this the
      !         two schemes would not even share a discretization.
      !
      !     (2) THE FIXED POINT IS THE SAME ONE.  Iterated to 1e-12, the two
      !         orderings land on the same ionization profile and the same
      !         field.  This is the statement the production change rests
      !         on: it may reorder the iteration, it may not move the answer.
      !
      !     (3) THE SUBSTITUTION REACHES IT IN FEWER PASSES.  The lagged map
      !         moves information inward one cell of the front per pass; the
      !         substitution carries the whole column in one traversal, and
      !         what is left for it to iterate is only each cell's OWN
      !         optical depth, which no ordering of the column can remove.
      !         The counts are measured here and asserted against bounds.
      !
      ! CONFIGURATION.  The slab of xuv_cell_mean_attenuation, made thick
      ! enough to hold an ionization front and resolved finely enough that
      ! one cell is a fraction of an optical depth: H I alone, thereis_He
      ! false, no metals, no secondary ionization, opa_pf = 1, a_tau = 0, a
      ! MONOCHROMATIC 20 eV spectrum so that one cell's optical depth is one
      ! number.  The chemistry is photoionization equilibrium at a fixed
      ! recombination coefficient,
      !     P_HI (1-x) n_H = alpha n_H^2 x^2 ,
      ! solved in closed form for x in each cell.  It is not EXHALE's
      ! network; it is the simplest composition that responds to the field
      ! the way the real one does, which is all the ordering statement needs.

      use global_parameters
      use species_table, only: n_mion, n_mphot
      use energy_vectors_construct, only: set_energy_vectors
      use utils,        only: calc_column_dens
      use utils_ion_eq, only: advance_starward_columns,                    &
                              photoionization_field_at_cell_H
      use assertion_report

      implicit none

      ! The slab is sized so that it holds an ionization front: the neutral
      ! column across it is many optical depths, and the recombination
      ! coefficient is set from the UNATTENUATED rate so that the outermost
      ! cell sits at x -> 1 while the innermost is neutral.  Both are
      ! properties of the test problem, not of any planet.
      real*8, parameter :: dr_geo   = 2.5d7      ! cell width [cm]
      real*8, parameter :: alpha_B  = 2.59d-13   ! case B at 1e4 K [cm^3/s]
      ! P/(alpha n_H) at the top of the slab: the ionization parameter the
      ! outermost cell sits at.  The density follows from it and from the
      ! unattenuated rate, so the front sits inside the slab whatever
      ! spectrum normalization the shared setup gives.
      real*8, parameter :: ion_par_top = 3.0d1
      real*8 :: n_H_tot                          ! H nuclei [cm^-3]
      ! The pass is finished when no cell's ionized fraction moved by more
      ! than this.  1e-12 is a few decades above the round-off limit cycle
      ! of a fraction of order one, so what the two schemes are compared at
      ! is their fixed point and not a partly converged state.
      real*8, parameter :: tol_fix  = 1.0d-12
      ! The tolerance the coupled source step of EXHALE_main actually stops
      ! at (csm_comp_tol), which is the pass count that costs a run.
      real*8, parameter :: tol_csm  = 1.0d-8
      integer, parameter :: max_pass = 2000

      real*8, dimension(:), allocatable :: x_jac, x_sub, x_ref
      real*8, dimension(:), allocatable :: nhi, xion, N1_face
      real*8, dimension(:), allocatable :: N1_ref, N15_ref, N2_ref, NTR_ref
      real*8, dimension(:), allocatable :: nzero
      real*8, dimension(:,:), allocatable :: nm_zero
      real*8  :: Nm_c(n_mphot)
      real*8  :: N1_c, N15_c, N2_c, NTR_c, NH2_c
      real*8  :: dmax, dev, xnew, P_top
      integer :: j, np_jac, np_sub, np_jac_csm, np_sub_csm

      ! --- the run-wide input, as input_read resolves it ------------------
      is_PL_sed    = .false.
      do_read_sed  = .false.
      is_monochr   = .true.
      e_low        = 20.0d0
      LEUV         = 30.42d0
      LX           = -3.0d2
      a_orb        = 0.02544d0*AU
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      thereis_He   = .false.
      thereis_HeITR       = .false.
      thereis_Xray        = .false.
      thereis_lowIP_metal = .false.
      use_sec_ion    = .false.
      sec_ion_active = .false.

      call set_energy_vectors

      N  = 250
      R0 = 1.0d10
      call allocate_grid_arrays
      dr_j   = dr_geo/R0
      opa_pf = 1.0d0

      allocate(x_jac(1-Ng:N+Ng), x_sub(1-Ng:N+Ng), x_ref(1-Ng:N+Ng))
      allocate(nhi(1-Ng:N+Ng), xion(1-Ng:N+Ng), N1_face(1-Ng:N+Ng))
      allocate(N1_ref(1-Ng:N+Ng), N15_ref(1-Ng:N+Ng))
      allocate(N2_ref(1-Ng:N+Ng), NTR_ref(1-Ng:N+Ng))
      allocate(nzero(1-Ng:N+Ng), nm_zero(1-Ng:N+Ng,n_mion))
      xion    = 0.0d0
      nzero   = 0.0d0
      nm_zero = 0.0d0

      ! ---------------------------------------------------------------- !
      ! (1) the two column integrations are one discretization
      ! ---------------------------------------------------------------- !
      ! A profile with structure, so that a wrong accumulation order cannot
      ! hide behind a uniform slab.
      do j = 1-Ng,N+Ng
         nhi(j) = 1.0d9*(0.2d0 + 0.8d0*dble(j-(1-Ng))/dble(N+2*Ng-1))
      enddo
      call calc_column_dens(nhi, nzero, nzero, nzero,                      &
                            N1_ref, N15_ref, N2_ref, NTR_ref)
      N1_c = 0.0d0; N15_c = 0.0d0; N2_c = 0.0d0
      NTR_c = 0.0d0; NH2_c = 0.0d0; Nm_c = 0.0d0
      call advance_starward_columns(N+Ng, 1-Ng, nhi, nzero, nzero, nzero,  &
               nzero, nm_zero, .false.,                                    &
               N1_c, N15_c, N2_c, NTR_c, NH2_c, Nm_c,                      &
               N1_face = N1_face)
      ! N1_face(j) is the column OUTSIDE cell j, i.e. calc_column_dens at
      ! j+1; the outermost cell has nothing above it.
      do j = 1-Ng,N+Ng-1
         call check_absolute('starward_column_matches_calc_column_dens',   &
                             N1_face(j), N1_ref(j+1), 0.0d0)
      enddo
      call check_absolute('starward_column_of_the_outermost_cell_is_zero', &
                          N1_face(N+Ng), 0.0d0, 0.0d0)
      call check_absolute('starward_column_leaves_the_whole_column',      &
                          N1_c, N1_ref(1-Ng), 0.0d0)

      ! ---------------------------------------------------------------- !
      ! (2) and (3): the two orderings, from the same neutral start
      ! ---------------------------------------------------------------- !
      ! The recombination coefficient that puts the top of the slab at the
      ! stated ionization parameter, from the rate an unattenuated cell of
      ! zero opacity sees.
      nhi = 0.0d0
      call photoionization_field_at_cell_H(N+Ng, 0.0d0, 0.0d0, 0.0d0,      &
               .false., P_top, dev, dmax, xnew)
      n_H_tot = P_top/(alpha_B*ion_par_top)
      nhi = n_H_tot
      write(*,'(a,es10.3,a,es10.3,a,f8.2)')                                &
        '  DIAGNOSTIC unattenuated rate [1/s] = ', P_top,                  &
        ' , n_H [cm^-3] = ', n_H_tot,                                      &
        ' , neutral slab optical depth = ',                                &
        s_hi(1)*1.0d-18*n_H_tot*dr_geo*dble(N+2*Ng)

      x_jac = 0.0d0
      call iterate_lagged(x_jac, tol_fix, np_jac)
      x_sub = 0.0d0
      call iterate_substitution(x_sub, tol_fix, np_sub)

      x_ref = x_jac
      dmax  = 0.0d0
      do j = 1-Ng,N+Ng
         dmax = max(dmax, abs(x_sub(j) - x_ref(j)))
      enddo
      write(*,'(a,es10.3)')                                                &
        '  DIAGNOSTIC largest difference of the two fixed points = ', dmax
      call check_absolute('the_two_orderings_have_one_fixed_point',        &
                          dmax, 0.0d0, 1.0d-9)

      ! The field at each fixed point, through the production routine.
      dev = 0.0d0
      do j = 1-Ng,N+Ng
         dev = max(dev, abs(rate_of(x_sub,j) - rate_of(x_ref,j))           &
                        /max(rate_of(x_ref,j),1.0d-99))
      enddo
      write(*,'(a,es10.3)')                                                &
        '  DIAGNOSTIC largest relative difference of the two fields = ', dev
      call check_absolute('the_two_orderings_have_one_field',              &
                          dev, 0.0d0, 1.0d-6)

      ! The pass counts at the tolerance the coupled source step stops at.
      x_jac = 0.0d0
      call iterate_lagged(x_jac, tol_csm, np_jac_csm)
      x_sub = 0.0d0
      call iterate_substitution(x_sub, tol_csm, np_sub_csm)
      write(*,'(a,i0,a,i0,a,i0,a,i0)')                                     &
        '  DIAGNOSTIC passes to 1e-8: lagged ', np_jac_csm,                &
        ', substitution ', np_sub_csm,                                     &
        ' ; to round-off: lagged ', np_jac, ', substitution ', np_sub
      ! MEASURED on this slab (250 cells, neutral optical depth 90, front
      ! near cell 190): the lagged map takes 24 passes to reach 1e-8 and the
      ! substitution 8, and 29 against 11 to reach round-off.  What is left
      ! for the substitution to iterate is each cell's OWN optical depth,
      ! which no ordering of the column can remove; that is why the ratio is
      ! about three here and not the grid size.  The bounds are set wide
      ! enough for the arithmetic to move and narrow enough that a change
      ! putting the field back on the lagged column would fail the first.
      call check_at_most('substitution_passes_to_the_csm_tolerance',       &
                         dble(np_sub_csm), 12.0d0)
      call check_at_least('lagged_passes_to_the_csm_tolerance',            &
                          dble(np_jac_csm), 18.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'photoionization_field_substitution: ',       &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'photoionization_field_substitution: '//              &
                     'all assertions passed'

      contains

      !--------------!

      double precision function equilibrium_fraction(P) result(x)
      ! The ionized fraction at which photoionization and case-B
      ! recombination balance in a cell of n_H nuclei:
      !   P (1-x) n_H = alpha n_H^2 x^2 ,  x in [0,1].
      real*8, intent(in) :: P
      real*8 :: a, b
      a = alpha_B*n_H_tot
      b = P
      x = (-b + sqrt(b*b + 4.0d0*a*b))/(2.0d0*a)
      end function equilibrium_fraction

      !--------------!

      double precision function rate_of(x, jc) result(P)
      ! The photoionization rate of cell jc in the field the profile x
      ! attenuates, through the production routine.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: x
      integer, intent(in) :: jc
      real*8 :: c1, c15, c2, cTR, cH2, cm(n_mphot)
      real*8 :: h1, hh, qq
      integer :: jj
      c1 = 0.0d0; c15 = 0.0d0; c2 = 0.0d0
      cTR = 0.0d0; cH2 = 0.0d0; cm = 0.0d0
      do jj = 1-Ng,N+Ng
         nhi(jj) = n_H_tot*(1.0d0 - x(jj))
      enddo
      call advance_starward_columns(N+Ng, jc+1, nhi, nzero, nzero, nzero,  &
               nzero, nm_zero, .false., c1, c15, c2, cTR, cH2, cm)
      call photoionization_field_at_cell_H(jc, c1, nhi(jc), 0.0d0,         &
               .false., P, h1, hh, qq)
      end function rate_of

      !--------------!

      subroutine iterate_lagged(x, tol, npass)
      ! The whole grid's field from one composition, then every cell's
      ! composition from that field: the map a whole-grid column integration
      ! gives.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: x
      real*8, intent(in) :: tol
      integer, intent(out) :: npass
      real*8  :: c1, c15, c2, cTR, cH2, cm(n_mphot)
      real*8  :: P, h1, hh, qq, move, xn
      integer :: jj
      npass = 0
      do
         npass = npass + 1
         do jj = 1-Ng,N+Ng
            nhi(jj) = n_H_tot*(1.0d0 - x(jj))
         enddo
         c1 = 0.0d0; c15 = 0.0d0; c2 = 0.0d0
         cTR = 0.0d0; cH2 = 0.0d0; cm = 0.0d0
         call advance_starward_columns(N+Ng, 1-Ng, nhi, nzero, nzero,      &
                  nzero, nzero, nm_zero, .false.,                          &
                  c1, c15, c2, cTR, cH2, cm, N1_face = N1_face)
         move = 0.0d0
         do jj = N+Ng,1-Ng,-1
            call photoionization_field_at_cell_H(jj, N1_face(jj), nhi(jj), &
                     0.0d0, .false., P, h1, hh, qq)
            xn   = equilibrium_fraction(P)
            move = max(move, abs(xn - x(jj)))
            x(jj) = xn
         enddo
         if (move .le. tol .or. npass .ge. max_pass) exit
      enddo
      end subroutine iterate_lagged

      !--------------!

      subroutine iterate_substitution(x, tol, npass)
      ! The same map with the column carried inward as the cells are solved:
      ! each cell's field is built from the composition already found for
      ! the cells outside it, and from its own entry composition for its own
      ! depth, exactly as the block sweep of ioniz_eq does at width one.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: x
      real*8, intent(in) :: tol
      integer, intent(out) :: npass
      real*8  :: c1, c15, c2, cTR, cH2, cm(n_mphot)
      real*8  :: P, h1, hh, qq, move, xn
      integer :: jj
      npass = 0
      do
         npass = npass + 1
         do jj = 1-Ng,N+Ng
            nhi(jj) = n_H_tot*(1.0d0 - x(jj))
         enddo
         c1 = 0.0d0; c15 = 0.0d0; c2 = 0.0d0
         cTR = 0.0d0; cH2 = 0.0d0; cm = 0.0d0
         move = 0.0d0
         do jj = N+Ng,1-Ng,-1
            call photoionization_field_at_cell_H(jj, c1, nhi(jj),          &
                     0.0d0, .false., P, h1, hh, qq)
            xn   = equilibrium_fraction(P)
            move = max(move, abs(xn - x(jj)))
            x(jj) = xn
            nhi(jj) = n_H_tot*(1.0d0 - xn)
            call advance_starward_columns(jj, jj, nhi, nzero, nzero,       &
                     nzero, nzero, nm_zero, .false.,                       &
                     c1, c15, c2, cTR, cH2, cm)
         enddo
         if (move .le. tol .or. npass .ge. max_pass) exit
      enddo
      end subroutine iterate_substitution

      end program photoionization_field_substitution
