      program coupled_block_linear_system
      ! THE LINEAR SYSTEM OF THE COUPLED BLOCK, MEASURED ON ONE FROZEN
      ! ITERATE (item L32 of docs/PLAN_20260917.md).
      !
      ! The block's trust-region step ends with its Krylov subspace
      ! exhausted above the linear tolerance it asked for. Three
      ! explanations are open: the directional action is inaccurate, the
      ! action is accurate and the banded preconditioner is too weak, or a
      ! coupling the band omits dominates. Only the first is a property of
      ! the map itself, and it is what this program measures.
      !
      ! WHAT IS HELD FIXED, and why the program is one process. The
      ! operator of one block iteration is the state Y, the composition
      ! seed f_sp, the column scales D, the row scales Drow, the species
      ! box, the active bounds, the pseudo-time shift and the frozen
      ! reconstruction weights. The radiation and chemistry caches of
      ! ionization_equilibrium, util_ion_eq and diffusive_photochemistry
      ! have no accessor a second process could be put back into, so the
      ! samples that have to share one operator are taken here, between one
      ! initialization and one exit, with the seed named at every call
      ! (eval_residual takes it as an argument and writes its answer to
      ! another).
      !
      ! THE FOUR MEASUREMENTS THIS PROGRAM MAKES
      !   (1) repeated residual and directional-action evaluations at the
      !       checkpoint, and again after an unrelated state has been
      !       evaluated in between. The reference is the control of
      !       src/tests/residual_determinism, 1.6e-12 in row-scale units.
      !   (2) a perturbation ladder for the action at four arc lengths, the
      !       maximum and the RMS of the error reported separately for the
      !       mass, momentum, energy, elemental and carrier rows, with the
      !       cell of each maximum. A forward difference at the arc is
      !       compared with a central difference over the same arc: the
      !       two differ by the curvature of the residual, so a smooth
      !       residual gives an error RISING with the arc, while rounding
      !       and a residual that is not smooth on the arc give one
      !       falling as the reciprocal of the arc. A family whose error
      !       has a minimum at the arc the route probes at is a family the
      !       route's arc is well chosen for.
      !   (3) homogeneity and additivity of the same action, the property
      !       GMRES and the Arnoldi relation are statements of.
      !   (7) on a small grid, the Jacobian assembled column by column from
      !       the same action and solved directly, as the reference the
      !       matrix-free step is read against. A diagnostic reference, not
      !       a proposal for production.
      !
      ! Measurements 4, 5 and 6 of the item are made by the solver's own
      ! diagnostics on the same checkpoint (EXHALE_GM_TRUE_RESIDUAL,
      ! EXHALE_GM_ORTHO, EXHALE_GM_HISTORY, EXHALE_GM_IMAGE_CHECK,
      ! EXHALE_KRYLOV_SIZE_SCAN, EXHALE_GM_CYCLES, EXHALE_GM_M), because
      ! those read the operator WITH the scaling, the preconditioner and
      ! the pseudo-time shift the solve builds, which a program outside the
      ! solve would have to rebuild and could then no longer claim to be
      ! the same operator.
      !
      ! It runs in a directory that holds input.inp, base.inp and the state
      ! in output/, like any run. It takes no step and writes no state.
      ! One line per assertion:  PASS|FAIL <name> measured= reference=
      ! The ladder and the additivity defect are printed to be read and are
      ! not gated: no tolerance for them was stated before they were
      ! measured, and one chosen afterwards would be a tolerance chosen to
      ! make a measurement pass.
      use global_parameters
      use Read_input,             only: input_read
      use Initialization,         only: init
      use lya_rt,                 only: lya_rt_allocate_arrays
      use excited_hydrogen,       only: excited_H_allocate_arrays
      use ionization_equilibrium, only: ioniz_eq_allocate_arrays
      use mol_rates,              only: h2_thermochemistry_init
      use molecular_infrared_cooling, only: molecular_infrared_init
      use steady_newton, only: nvar_jac, neq_newton, pack_U,              &
           pack_species_rows, eval_residual, cell_state_scales,           &
           cell_row_scales, set_transported_species_rows,                 &
           read_species_unknown_space_controls,                           &
           freeze_species_unknown_box, jacobian_action_of_direction,      &
           probe_length_of_the_jacobian_action, name_of_unknown,          &
           largest_step_inside_the_species_box,                           &
           species_unknowns_outside_their_bounds
      implicit none

      ! The five row families the item asks for each error separately on.
      integer, parameter :: n_family = 5
      character(len=10), parameter, dimension(n_family) :: family_name =  &
           (/ 'mass      ', 'momentum  ', 'energy    ', 'elemental ',     &
              'carrier   ' /)
      ! The noise floor the repeated evaluations are read against: the
      ! control of src/tests/residual_determinism, F(Y) twice with nothing
      ! in between, 1.62e-12 in row-scale units (READ, docs/TO_BE_DONE.md
      ! entry (AG)). A spread at or below it is the arithmetic of the
      ! evaluation and not a state carried between samples.
      real*8, parameter :: noise_floor = 1.6d-12
      ! The arcs of the ladder, as multiples of the arc the route itself
      ! probes at, sqrt(eps)(1 + ||Y||).
      integer, parameter :: n_arc = 4
      real*8, parameter, dimension(n_arc) :: arc_factor =                 &
           (/ 1.0d-2, 1.0d-1, 1.0d0, 1.0d1 /)

      real*8, allocatable :: W(:,:), u(:,:), f_sp(:,:), f_work(:,:)
      real*8, allocatable :: Y(:), F0(:), F1(:), Fa(:), Fb(:), Fother(:)
      real*8, allocatable :: D(:), Drow(:), v(:), vhat(:)
      real*8, allocatable :: Jv1(:), Jv2(:), Jv3(:)
      real*8, allocatable :: afwd(:), acen(:)
      real*8, allocatable :: Yother(:), Ytry(:)
      real*8, allocatable :: u1dir(:), u2dir(:), a1(:), a2(:), asum(:)
      real*8, allocatable :: heat(:), cool(:)
      real*8 :: arc, eps0, vn, ynorm, room, room_bwd
      real*8 :: fam_max(n_family), fam_rms(n_family)
      integer :: fam_cell(n_family), fam_n(n_family)
      real*8 :: spread_rep, spread_after, spread_jv, spread_jv_after
      real*8 :: worst_add, worst_hom
      integer :: neq, i, k, idir, iarc, n_fail, n_out_box
      logical :: ok
      character(len=24) :: dirname(3)
      ! Measurement 7, armed by EXHALE_L32_ASSEMBLE=1 and meant for a grid
      ! of about fifty cells: the Jacobian assembled column by column from
      ! the same action, and a direct solve of the scaled shifted system.
      real*8, allocatable :: Amat(:,:), Akeep(:,:), Acol(:), Aimg(:)
      real*8, allocatable :: rhs(:), dsol(:)
      integer, allocatable :: ipv(:)
      real*8  :: idtau, true_rel
      integer :: lpinfo
      character(len=32) :: env

      n_fail = 0

      ! ---- the run this program measures in ------------------------ !
      call input_read
      open(unit = outfile, file = 'EXHALE_setup.out')
      call lya_rt_allocate_arrays
      call excited_H_allocate_arrays
      call ioniz_eq_allocate_arrays
      if (mol_ir_bands) call molecular_infrared_init(T0)
      call h2_thermochemistry_init
      allocate(W(3,1-Ng:N+Ng), u(3,1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(f_work(1-Ng:N+Ng,n_species))
      allocate(heat(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      call init(W,u,f_sp)

      ! ---- the block's own system, and the frozen checkpoint ------- !
      ! The registry is what makes the transported balances unknowns of
      ! the Newton vector; without it nvar_jac is 3 and this program would
      ! measure the alternation's system instead of the block's.
      call set_transported_species_rows(.true.)
      call read_species_unknown_space_controls
      neq = neq_newton()
      allocate(Y(neq), F0(neq), F1(neq), Fa(neq), Fb(neq), Fother(neq))
      allocate(D(neq), Drow(neq), v(neq), vhat(neq))
      allocate(Jv1(neq), Jv2(neq), Jv3(neq), afwd(neq), acen(neq))
      allocate(Yother(neq), Ytry(neq))
      allocate(u1dir(neq), u2dir(neq), a1(neq), a2(neq), asum(neq))
      call pack_U(u, Y)
      call pack_species_rows(u, f_sp, Y)
      call freeze_species_unknown_box(Y)
      call eval_residual(Y, f_sp, f_work, F0, heat, cool)
      call cell_state_scales(Y, D)
      call cell_row_scales(Y, D, Drow, u)
      ynorm = sqrt(sum(Y*Y))
      arc   = probe_length_of_the_jacobian_action(Y)

      write(*,'(A,I0,A,I0,A,I0)') '  the block system: N = ', N,          &
           ', unknowns per cell = ', nvar_jac, ', rows = ', neq
      write(*,'(A,ES14.6,A,ES14.6)') '  ||Y|| = ', ynorm,                 &
           ', the arc the action probes at = ', arc
      write(*,'(A,ES14.6)') '  ||F/Drow||_inf of the checkpoint = ',      &
           maxval(abs(F0/Drow))
      Ytry = Y
      call species_unknowns_outside_their_bounds(Ytry, n_out_box, .false.)
      write(*,'(A,I0)') '  species unknowns outside their box = ', n_out_box

      ! ================================================================ !
      ! (1) IS ONE OPERATOR BEING SAMPLED?                               !
      ! ================================================================ !
      ! The residual at the checkpoint, twice with nothing in between and
      ! once more after an unrelated state has been evaluated through the
      ! same caches. The measure is the largest difference in the units
      ! each row is judged in, which is the unit the noise floor is stated
      ! in.
      call eval_residual(Y, f_sp, f_work, F1, heat, cool)
      spread_rep = maxval(abs(F1 - F0)/max(abs(Drow), tiny(1.0d0)))

      ! An unrelated state: the same column with its internal energy raised
      ! by a thousandth. It is admissible by construction (positive density,
      ! positive pressure, composition untouched) and it is not near the
      ! checkpoint, so whatever it leaves in a cache is left there.
      Yother = Y
      do k = 3, neq, nvar_jac
         Yother(k) = Yother(k)*1.001d0
      enddo
      call eval_residual(Yother, f_sp, f_work, Fother, heat, cool)
      call eval_residual(Y, f_sp, f_work, F1, heat, cool)
      spread_after = maxval(abs(F1 - F0)/max(abs(Drow), tiny(1.0d0)))

      ! The directional action along one fixed direction, the same three
      ! times, with the unrelated state evaluated between the second and
      ! the third.
      call deterministic_direction(1, vhat)
      v  = D*vhat
      call jacobian_action_of_direction(Y, F0, f_sp, v, Jv1, ok)
      call jacobian_action_of_direction(Y, F0, f_sp, v, Jv2, ok)
      call eval_residual(Yother, f_sp, f_work, Fother, heat, cool)
      call jacobian_action_of_direction(Y, F0, f_sp, v, Jv3, ok)
      spread_jv       = relative_spread(Jv1, Jv2)
      spread_jv_after = relative_spread(Jv1, Jv3)

      write(*,'(A)') '  (1) repeated evaluations at the checkpoint'
      write(*,'(A,ES12.4)') '      residual, twice with nothing in'//     &
           ' between (row-scale units): ', spread_rep
      write(*,'(A,ES12.4)') '      residual, with an unrelated state'//   &
           ' in between            : ', spread_after
      write(*,'(A,ES12.4)') '      action,  twice with nothing in'//      &
           ' between (relative)      : ', spread_jv
      write(*,'(A,ES12.4)') '      action,  with an unrelated state'//    &
           ' in between            : ', spread_jv_after
      call verdict(spread_after .le. max(noise_floor, 2.0d0*spread_rep),  &
           'the_residual_is_one_state_function', spread_after,            &
           noise_floor, n_fail)
      call verdict(spread_jv_after .le. max(1.0d-12, 2.0d0*spread_jv),    &
           'the_action_is_one_map', spread_jv_after, spread_jv, n_fail)

      ! ================================================================ !
      ! (2) THE PERTURBATION LADDER, BY ROW FAMILY                       !
      ! ================================================================ !
      ! Three directions: the one the Krylov cycle starts from, which is
      ! the scaled right-hand side of the step; a carrier direction; an
      ! energy direction. The last two are the pair I1 measured, so the
      ! numbers stand beside its means.
      dirname(1) = 'the right-hand side'
      dirname(2) = 'the carrier rows'
      dirname(3) = 'the energy rows'
      write(*,'(A)') '  (2) the action against a central difference of'// &
           ' the same map, by row family'
      do idir = 1, 3
         call direction_of_kind(idir, F0, Drow, vhat)
         v  = D*vhat
         vn = sqrt(sum(v*v))
         if (vn .le. 0.0d0) cycle
         ! THE ROOM THE FEASIBLE SET LEAVES ALONG THIS DIRECTION, on both
         ! sides, so that a rung whose step the species box cut shorter
         ! than the arc asked for is read as such and not as a rung of the
         ! ladder.
         room     = largest_step_inside_the_species_box(Y,  v)
         room_bwd = largest_step_inside_the_species_box(Y, -v)
         write(*,'(A,A,A,ES11.3,A,ES11.3)') '      direction: ',          &
              trim(dirname(idir)), ', room to the box forward ',          &
              room*vn/arc, ' arcs, backward ', room_bwd*vn/arc
         do iarc = 1, n_arc
            eps0 = arc_factor(iarc)*arc/vn
            ! POSITIVITY AND THE SHARED ABUNDANCE BOUNDS. The step is cut
            ! to what the species box and the element budget leave on both
            ! sides, so neither sample leaves the feasible set and the two
            ! sides of the central difference are the same length.
            eps0 = min(eps0, 0.9d0*room, 0.9d0*room_bwd)
            if (eps0 .le. 0.0d0) cycle
            Ytry = Y + eps0*v
            call eval_residual(Ytry, f_sp, f_work, Fa, heat, cool)
            Ytry = Y - eps0*v
            call eval_residual(Ytry, f_sp, f_work, Fb, heat, cool)
            afwd = (Fa - F0)/eps0
            acen = (Fa - Fb)/(2.0d0*eps0)
            call family_errors(afwd, acen, fam_max, fam_rms, fam_cell,    &
                               fam_n)
            write(*,'(A,ES9.2,A,ES11.3)') '        arc asked x',          &
                 arc_factor(iarc), ', arc taken x', eps0*vn/arc
            do k = 1, n_family
               if (fam_n(k) .eq. 0) cycle
               write(*,'(A,A,A,ES11.3,A,ES11.3,A,I4,A,I0,A)')             &
                    '          ',                                         &
                    family_name(k), '  max ', fam_max(k), '  rms ',       &
                    fam_rms(k), '  worst cell ', fam_cell(k), '  (',      &
                    fam_n(k), ' rows)'
            enddo
         enddo
      enddo

      ! ================================================================ !
      ! (3) HOMOGENEITY AND ADDITIVITY                                   !
      ! ================================================================ !
      ! A(c v) = c A(v) is what the probe length is built to give, and
      ! A(v1 + v2) = A(v1) + A(v2) is the property the Arnoldi relation
      ! needs and the probe length does not give by construction.
      call deterministic_direction(11, vhat);  u1dir = D*vhat
      call deterministic_direction(29, vhat);  u2dir = D*vhat
      call jacobian_action_of_direction(Y, F0, f_sp, u1dir, a1, ok)
      call jacobian_action_of_direction(Y, F0, f_sp, 3.0d0*u1dir, a2, ok)
      worst_hom = relative_spread(a2, 3.0d0*a1)
      call jacobian_action_of_direction(Y, F0, f_sp, u2dir, a2, ok)
      call jacobian_action_of_direction(Y, F0, f_sp, u1dir + u2dir, asum, &
                                        ok)
      call jacobian_action_of_direction(Y, F0, f_sp, u1dir, a1, ok)
      worst_add = norm_ratio(asum - a1 - a2, asum)
      write(*,'(A)') '  (3) the action as a linear map'
      write(*,'(A,ES12.4)') '      homogeneity, A(3v) against 3A(v),'//   &
           ' relative     : ', worst_hom
      write(*,'(A,ES12.4)') '      additivity, ||A(v1+v2) - A(v1) -'//    &
           ' A(v2)||/||A(v1+v2)||: ', worst_add
      call verdict(worst_hom .le. 1.0d-12, 'the_action_is_homogeneous',   &
           worst_hom, 0.0d0, n_fail)

      ! ================================================================ !
      ! (7) THE ASSEMBLED JACOBIAN AND A DIRECT SOLVE, ON A SMALL GRID    !
      ! ================================================================ !
      ! A column of the scaled operator is A_z e_i = idtau e_i
      ! + D^-1 J (D e_i) / Drow, sampled with the SAME action the Krylov
      ! cycle uses, so the assembled matrix and the matrix-free operator
      ! are one map by construction and the direct solve says what the
      ! Krylov cycle would reach with a complete subspace. It costs one
      ! residual evaluation per unknown and is a diagnostic reference, not
      ! a proposal for production.
      call get_environment_variable('EXHALE_L32_ASSEMBLE', env)
      if (trim(env) .eq. '1') then
         idtau = 0.0d0
         call get_environment_variable('EXHALE_L32_IDTAU', env)
         if (len_trim(env) .gt. 0) read(env,*) idtau
         allocate(Amat(neq,neq), Akeep(neq,neq), Acol(neq), Aimg(neq),  &
                  rhs(neq), dsol(neq), ipv(neq))
         do i = 1, neq
            vhat = 0.0d0;  vhat(i) = 1.0d0
            call jacobian_action_of_direction(Y, F0, f_sp, D*vhat, Acol,  &
                                              ok)
            Amat(:,i) = Acol/Drow
            Amat(i,i) = Amat(i,i) + idtau
         enddo
         Akeep = Amat
         rhs  = -F0/Drow
         dsol = rhs
         call dgesv(neq, 1, Amat, neq, ipv, dsol, neq, lpinfo)
         if (lpinfo .ne. 0) then
            write(*,'(A,I0)') '  (7) the assembled system is singular,'// &
                 ' dgesv info = ', lpinfo
         else
            ! TWO RESIDUALS OF ONE STEP, and they say different things.
            ! The first is the residual of the linear system that was
            ! factored, b - (assembled A) x: it says whether the direct
            ! solve solved what it was given. The second replaces the
            ! matrix product by one fresh product of the matrix-free
            ! action along the same step: it says whether the assembled
            ! matrix is the action at all, which it is only if the action
            ! is a linear map of its direction.
            Aimg     = matmul(Akeep, dsol)
            true_rel = norm_ratio(rhs - Aimg, rhs)
            call jacobian_action_of_direction(Y, F0, f_sp, D*dsol, Acol,  &
                                              ok)
            Acol = Acol/Drow + idtau*dsol
            write(*,'(A,I0,A,I0,A)') '  (7) direct solve of the '//       &
                 'assembled ', neq, ' by ', neq, ' system'
            write(*,'(A,ES12.4)') '      ||b - (assembled A) x||/||b||'// &
                 '                     : ', true_rel
            write(*,'(A,ES12.4)') '      ||A(x) - (assembled A) x||/'//   &
                 '||(assembled A) x||: ', norm_ratio(Acol - Aimg, Aimg)
            write(*,'(A,ES12.4)') '      ||b - A(x)||/||b||   (the'//     &
                 ' step through the action) : ', norm_ratio(rhs - Acol,   &
                 rhs)
            write(*,'(A,ES12.4)') '      ||x|| of the direct step = ',    &
                 sqrt(sum(dsol*dsol))
            ! THE SAME COMPARISON ALONG A STEP OF THE LENGTH THE SOLVE
            ! TAKES. The direct step is long, and a secant over a fixed
            ! arc departs from the tangent by more the further the
            ! direction reaches, so the defect above would be read as a
            ! property of the step's length alone unless the same
            ! direction is measured at unit length beside it.
            vhat = dsol/max(sqrt(sum(dsol*dsol)), tiny(1.0d0))
            Aimg = matmul(Akeep, vhat)
            call jacobian_action_of_direction(Y, F0, f_sp, D*vhat, Acol,  &
                                              ok)
            Acol = Acol/Drow + idtau*vhat
            write(*,'(A,ES12.4)') '      the same at unit length: '//     &
                 '||A(v) - (assembled A) v||/||Av||: ',                   &
                 norm_ratio(Acol - Aimg, Aimg)
            call verdict(true_rel .le. 1.0d-6,                            &
                 'the_direct_solve_solves_the_assembled_system',          &
                 true_rel, 1.0d-6, n_fail)
         endif
      endif

      close(outfile)
      if (n_fail .ne. 0) then
         write(*,'(A)') 'coupled_block_linear_system: FAILED'
         stop 1
      endif
      write(*,'(A)') 'coupled_block_linear_system: PASSED'

      contains

      ! ------------------------------------------------------!

      subroutine verdict(passed, name, measured, reference, nf)
      ! One assertion line in the form every suite of this tree prints.
      logical,          intent(in)    :: passed
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference
      integer,          intent(inout) :: nf
      if (passed) then
         write(*,'(A,A,A,ES12.4,A,ES12.4)') 'PASS ', name, ' measured=',  &
              measured, ' reference=', reference
      else
         write(*,'(A,A,A,ES12.4,A,ES12.4)') 'FAIL ', name, ' measured=',  &
              measured, ' reference=', reference
         nf = nf + 1
      endif
      end subroutine verdict

      ! ------------------------------------------------------!

      subroutine deterministic_direction(seed, vv)
      ! A reproducible direction of unit length in the SCALED unknowns, so
      ! that two runs of one binary sample the same direction. The entries
      ! come from a linear congruential recurrence mapped onto [-1,1); the
      ! sequence's quality is irrelevant, it is a direction.
      integer,              intent(in)  :: seed
      real*8, dimension(:), intent(out) :: vv
      integer*8 :: st
      integer   :: ii
      real*8    :: nn
      st = int(seed, 8)
      do ii = 1, size(vv)
         st = iand(1664525_8*st + 1013904223_8, 4294967295_8)
         vv(ii) = 2.0d0*(real(st, 8)/4294967296.0d0) - 1.0d0
      enddo
      nn = sqrt(sum(vv*vv))
      if (nn .gt. 0.0d0) vv = vv/nn
      end subroutine deterministic_direction

      ! ------------------------------------------------------!

      subroutine direction_of_kind(kind, Fv, Dr, vv)
      ! The unit direction of the scaled space the ladder samples along.
      !   1  the scaled right-hand side of the Newton step, -F/Drow/Dc:
      !      the direction the Krylov cycle starts from.
      !   2  the carrier unknowns alone.
      !   3  the energy unknowns alone.
      ! Kinds 2 and 3 are the pair the I1 report measured, so the numbers
      ! here stand beside its means.
      integer,                       intent(in)  :: kind
      real*8, dimension(:),          intent(in)  :: Fv, Dr
      real*8, dimension(:),          intent(out) :: vv
      integer :: ii, jj, kk2
      real*8  :: nn
      vv = 0.0d0
      select case (kind)
      case (1)
         vv = -Fv/max(abs(Dr), tiny(1.0d0))
      case (2)
         do ii = 1, size(vv)
            jj  = (ii - 1)/nvar_jac + 1
            kk2 = ii - nvar_jac*(jj - 1)
            if (kk2 .gt. 3) vv(ii) = 1.0d0
         enddo
      case (3)
         do ii = 1, size(vv)
            jj  = (ii - 1)/nvar_jac + 1
            kk2 = ii - nvar_jac*(jj - 1)
            if (kk2 .eq. 3) vv(ii) = 1.0d0
         enddo
      end select
      nn = sqrt(sum(vv*vv))
      if (nn .gt. 0.0d0) vv = vv/nn
      end subroutine direction_of_kind

      ! ------------------------------------------------------!

      subroutine family_errors(af, ac, fmax, frms, fcell, fn)
      ! The error of the forward difference against the central one, split
      ! by the five row families and reported as a maximum and an RMS
      ! separately, with the cell the maximum sits in.
      !
      ! A ROW WHOSE TWO SIDES BOTH STAND AT THE ROUNDING FLOOR CARRIES NO
      ! DERIVATIVE TO COMPARE: the mass row of a cell does not move with
      ! the composition at all, and the quotient of two rounding numbers
      ! reads whatever the arithmetic gives. A row enters the statistics
      ! only where its central difference reaches 1e-3 of the largest
      ! central difference of its own family, which is the threshold the
      ! I1 report used for "a row the direction moves", and the count of
      ! rows that did enter is printed beside the numbers.
      real*8, dimension(:), intent(in)  :: af, ac
      real*8, dimension(n_family), intent(out) :: fmax, frms
      integer, dimension(n_family), intent(out) :: fcell, fn
      real*8, dimension(n_family) :: biggest
      real*8  :: e
      integer :: ii, jj, kk2, ic
      fmax = 0.0d0;  frms = 0.0d0;  fcell = 0;  fn = 0
      biggest = 0.0d0
      do ii = 1, size(ac)
         ic = family_of_unknown(ii)
         biggest(ic) = max(biggest(ic), abs(ac(ii)))
      enddo
      do ii = 1, size(ac)
         jj  = (ii - 1)/nvar_jac + 1
         kk2 = ii - nvar_jac*(jj - 1)
         ic  = family_of_unknown(ii)
         if (biggest(ic) .le. 0.0d0) cycle
         if (abs(ac(ii)) .lt. 1.0d-3*biggest(ic)) cycle
         e = abs(af(ii) - ac(ii))/abs(ac(ii))
         fn(ic)   = fn(ic) + 1
         frms(ic) = frms(ic) + e*e
         if (e .gt. fmax(ic)) then
            fmax(ic)  = e
            fcell(ic) = jj
         endif
      enddo
      do ii = 1, n_family
         if (fn(ii) .gt. 0) frms(ii) = sqrt(frms(ii)/dble(fn(ii)))
      enddo
      end subroutine family_errors

      ! ------------------------------------------------------!

      integer function family_of_unknown(i) result(ic)
      ! Which of the five row families a flat unknown index belongs to.
      ! The species slots split into the elemental and the carrier
      ! balances, which are different equations with different scales and
      ! the item asks for them apart.
      integer, intent(in) :: i
      integer :: j, k
      character(len=48) :: txt
      j = (i - 1)/nvar_jac + 1
      k = i - nvar_jac*(j - 1)
      if (k .le. 3) then
         ic = k
         return
      endif
      call name_of_unknown(i, txt)
      if (index(txt, 'carrier') .gt. 0) then
         ic = 5
      else
         ic = 4
      endif
      end function family_of_unknown

      ! ------------------------------------------------------!

      real*8 function relative_spread(a, b) result(s)
      ! The largest relative difference of two vectors, read against the
      ! largest entry of the first: a comparison entry by entry would be
      ! a quotient of two rounding numbers wherever the entry is small.
      real*8, dimension(:), intent(in) :: a, b
      real*8 :: amax
      amax = maxval(abs(a))
      if (amax .le. 0.0d0) then
         s = 0.0d0
      else
         s = maxval(abs(a - b))/amax
      endif
      end function relative_spread

      ! ------------------------------------------------------!

      real*8 function norm_ratio(a, b) result(s)
      ! ||a||/||b||, zero where b is zero.
      real*8, dimension(:), intent(in) :: a, b
      real*8 :: bn
      bn = sqrt(sum(b*b))
      if (bn .le. 0.0d0) then
         s = 0.0d0
      else
         s = sqrt(sum(a*a))/bn
      endif
      end function norm_ratio

      end program coupled_block_linear_system
