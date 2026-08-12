      module steady_newton
      ! Steady-state residual in vector form, F(Y) = 0, for a direct
      ! (Newton / pseudo-transient-continuation) solve. Increment (ii)-2:
      ! pack/unpack + newton_residual. The Jacobian and the PTC driver are
      ! added in later increments.
      !
      ! Unknowns: the hydro conserved variables on the PHYSICAL cells
      ! j = 1..N (neq = 3*N). The ghost cells (1-Ng..0 lower, N+1..N+Ng
      ! upper) are NOT unknowns: Apply_BC fills them from the interior each
      ! residual evaluation (lower = fixed rho_bc / p = 1+dp_bc with the
      ! v-valve max(v_1,0); upper = zero-gradient / WENO3 extrapolation).
      !
      ! Ionization/temperature are eliminated LOCALLY inside the residual:
      ! given Y, ioniz_eq re-solves the equilibrium fractions and returns the
      ! consistent heat/cool, so F(Y) is a pure-hydro residual of size 3*N
      ! (the "method A" choice from ATES_sundials_solver_plan.md).

      use global_parameters
      use species_table,  only: n_mion
      use Conversion,      only: U_to_W
      use composition,     only: get_species_densities, comp_T_from_p
      use BC_Apply,        only: Apply_BC
      use ionization_equilibrium, only: ioniz_eq
      use excited_hydrogen,       only: excited_H_update
      use steady_residual_mod,    only: assemble_residual

      implicit none
      private
      public :: neq_newton, pack_U, unpack_U, newton_residual,         &
                eval_residual, frozen_residual, build_banded_jac,      &
                band_matvec, kl_jac, ku_jac, ncolor_jac, solve_steady_ptc, &
                jv_product, solve_steady_jfnk, set_base_fix

      ! Banded-Jacobian geometry. WENO3 reconstruction reaches +/-2 cells,
      ! and within a cell all 3 hydro variables couple, so in the flat
      ! ordering Y(3*(j-1)+k) a column couples to rows within
      ! |row-col| <= 3*2 + 2 = 8. Half-bandwidths kl = ku = 8; a graph
      ! coloring with stride ncolor = kl+ku+1 = 17 makes same-color columns'
      ! row supports disjoint, so the FROZEN-radiation residual (strictly
      ! local) is captured exactly by 17 finite-difference probes.
      integer, parameter :: kl_jac = 8, ku_jac = 8
      integer, parameter :: ncolor_jac = kl_jac + ku_jac + 1   ! = 17

      ! Frozen-base option: the residual of the first nfix_base physical
      ! cells is replaced by the anchor row F_j = Y_j - Yfix_j, pinning them
      ! to the warm-start state, so that Newton solves for the WIND on top of
      ! a given base state. Set via set_base_fix before a solve.
      ! NB: it was previously assumed that the cell-1 momentum row cannot be
      ! satisfied at all, the lower BC (hard-pinned rho_bc/p_bc ghosts + the
      ! v-valve max(v,0)) leaving it non-zero at any steady interior solution.
      ! Direct measurement on the HD 209458 b hand-off state does not support
      ! that: R(2,1) is a four-order cancellation of gravity against the
      ! pressure gradient whose remainder is 6.5e-5 of the gravity term, and a
      ! 2.9 ppm change of the ghost pressure drives it to zero
      ! (docs/newton_scaling_and_base_wall.md). The option is kept as an
      ! experiment, not as a remedy for a base that cannot be solved.
      integer :: nfix_base = 0
      real*8, allocatable :: Yfix_base(:)

      ! Explicit interfaces for the external LAPACK banded-LU routines used by
      ! the direct/preconditioned Newton solves below (double-precision,
      ! general band form). The array dummies are assumed-size so that the
      ! call sites, which pass the right-hand side as a rank-1 vector
      ! (nrhs = 1), match without any argument change; this only gives the
      ! compiler kind/rank/intent information and does not alter the call.
      interface
         subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            integer,          intent(in)    :: m, n, kl, ku, ldab
            real*8,           intent(inout) :: ab(ldab,*)
            integer,          intent(out)   :: ipiv(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrf

         subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            character(1),     intent(in)    :: trans
            integer,          intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real*8,           intent(in)    :: ab(ldab,*)
            integer,          intent(in)    :: ipiv(*)
            real*8,           intent(inout) :: b(*)
            integer,          intent(out)   :: info
         end subroutine dgbtrs
      end interface

      contains

      ! ------------------------------------------------------!

      subroutine set_base_fix(Y, nfix)
      ! Anchor the first nfix physical cells at their current values.
      real*8, dimension(3*N), intent(in) :: Y
      integer, intent(in) :: nfix
      nfix_base = max(0, min(nfix, N-10))
      if (.not. allocated(Yfix_base)) allocate(Yfix_base(3*N))
      Yfix_base = Y
      end subroutine set_base_fix

      ! ------------------------------------------------------!

      subroutine cell_state_scales(Y, D)
      ! Characteristic scale of each conserved unknown IN ITS OWN CELL, from
      ! the local state alone:
      !
      !   mass        D = rho_j
      !   momentum    D = |rho v|_j + rho_j c_s,j = rho_j (|v|_j + c_s,j)
      !   energy      D = E_j
      !
      ! with c_s = sqrt(g p/rho) and p = (g-1)(E - (rho v)^2/2rho), so D is a
      ! function of Y only. It sets the finite-difference step sizes, the
      ! scaled Newton system D^-1 J D, and the line-search merit ||D^-1 F||_2.
      !
      ! rho and E are positive definite and are their own scales. The momentum
      ! density is not: it passes through zero wherever the flow reverses, and
      ! at the hand-off from marching the whole sub-sonic region below the
      ! stagnation point is still infalling, so |rho v| vanishes there and is
      ! no measure of how large a momentum residual is. The scale that does
      ! not vanish is rho times the fastest characteristic speed of the Euler
      ! system, |v| + c_s: the momentum density the cell carries when moved at
      ! its own signal speed. This is parameter-free -- no floor fraction to
      ! choose -- and it is local, so no single cell can set the scale of the
      ! whole domain.
      real*8, dimension(3*N), intent(in)  :: Y
      real*8, dimension(3*N), intent(out) :: D
      real*8  :: rho, mom, ene, pgas, cs
      integer :: j, i1, i2, i3
      do j = 1, N
         i1 = 3*(j-1) + 1;  i2 = i1 + 1;  i3 = i1 + 2
         rho = Y(i1);  mom = Y(i2);  ene = Y(i3)
         if (rho .gt. 0.0d0) then
            pgas = (g - 1.0d0)*(ene - 0.5d0*mom*mom/rho)
            cs   = sqrt(g*max(pgas, 0.0d0)/rho)
         else
            ! Not a physical state (the line search rejects these); fall back
            ! to the bare magnitudes rather than taking a root of a negative.
            rho = abs(rho);  cs = 0.0d0
         endif
         D(i1) = max(rho,                    tiny(1.0d0))
         D(i2) = max(abs(mom) + rho*cs,      tiny(1.0d0))
         D(i3) = max(abs(ene),               tiny(1.0d0))
      enddo
      end subroutine cell_state_scales

      ! ------------------------------------------------------!

      subroutine apply_base_fix(Y, Fvec)
      ! Replace the residual rows of the frozen base cells by anchor rows.
      real*8, dimension(3*N), intent(in)    :: Y
      real*8, dimension(3*N), intent(inout) :: Fvec
      integer :: j, k
      do j = 1, nfix_base
         do k = 1, 3
            Fvec(3*(j-1)+k) = Y(3*(j-1)+k) - Yfix_base(3*(j-1)+k)
         enddo
      enddo
      end subroutine apply_base_fix

      ! ------------------------------------------------------!

      integer function neq_newton()
      ! Number of Newton unknowns = 3 hydro variables over physical cells.
      neq_newton = 3*N
      end function neq_newton

      ! ------------------------------------------------------!

      subroutine pack_U(u, Y)
      ! Physical-cell conserved variables -> flat unknown vector.
      ! Y(3*(j-1)+k) = u(j,k), j = 1..N, k = 1..3.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3*N),         intent(out) :: Y
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Y(3*(j-1)+k) = u(k,j)
         enddo
      enddo
      end subroutine pack_U

      ! ------------------------------------------------------!

      subroutine unpack_U(Y, u)
      ! Flat unknown vector -> physical-cell conserved variables.
      ! Ghost cells are left untouched (Apply_BC sets them).
      real*8, dimension(3*N),         intent(in)    :: Y
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            u(k,j) = Y(3*(j-1)+k)
         enddo
      enddo
      end subroutine unpack_U

      ! ------------------------------------------------------!

      subroutine eval_residual(Y, f_sp, Fvec, heat, cool, n_part)
      ! Full steady residual F(Y) with local ionization-equilibrium
      ! elimination, AND the heat/cool it used (so the caller can FREEZE the
      ! radiation when building the banded Jacobian).
      !   unpack Y -> u(1:N); Apply_BC fills ghosts;
      !   (rho,v,p) -> densities, T; refresh excited-H; ioniz_eq -> heat,cool;
      !   R = assemble_residual(u, n_part, heat, cool); pack R(1:N) -> Fvec.
      ! f_sp is updated in place to the equilibrium fractions (intent inout).
      ! n_part = n_tot + n_e is returned too, so a caller that later builds a
      ! FROZEN-radiation residual can hand the same particle count back and
      ! keep T = p/n_part (hence the transport coefficients) consistent.
      real*8, dimension(3*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(inout) :: f_sp
      real*8, dimension(3*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng),            intent(out)   :: heat, cool
      real*8, dimension(1-Ng:N+Ng), optional,  intent(out)   :: n_part

      real*8, dimension(3,1-Ng:N+Ng) :: u, W, R
      real*8, dimension(1-Ng:N+Ng)   :: rho, v, p, T
      real*8, dimension(1-Ng:N+Ng)   :: nhi, nhii, nhei, nheii, nheiii, nheiTR
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng)   :: ne, n_tot, eta
      real*8 :: rel_change

      u = 0.0d0
      call unpack_U(Y, u)
      ! Composition of the INTERIOR of Y, evaluated before the ghosts are
      ! filled: Apply_BC reads n_part_cell1 for the continuous-temperature base
      ! ghost (T(1) = p(1)/n_part_cell1) and that global is written by
      ! get_species_densities. Without this call it still holds the value left
      ! by the PREVIOUS residual evaluation, so F would depend on the previous Y
      ! as well as on Y -- and a finite-difference Jacobian column would then
      ! mix two states. Cost is one array pass; ioniz_eq below dominates.
      ! u(1,:) is already the density, so no U_to_W here: the ghosts are still
      ! zero at this point and U_to_W would divide by them.
      rho = u(1,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call Apply_BC(u)           ! fill ghosts from the interior

      call U_to_W(u, W)
      rho = W(1,:);  v = W(2,:);  p = W(3,:)
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)
      call comp_T_from_p(p, n_tot, ne, T)
      if (use_excited_H) call excited_H_update(T,rho,f_sp,v,rel_change)
      call ioniz_eq(T,rho,f_sp,heat,cool,eta)
      ! Refresh the particle count from the equilibrium fractions, exactly as
      ! the marching loop does before its transport stage, so the residual's
      ! T = p/(n_tot+n_e) is the same temperature the marching step diffuses.
      ! Only the transport terms read it; with them off nothing changes.
      call get_species_densities(rho,f_sp,nhi,nhii,nhei,nheii,         &
                                 nheiii,nheiTR,nm,ne,n_tot)

      call assemble_residual(u, n_tot + ne, heat, cool, R)
      call pack_R(R, Fvec)
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      if (present(n_part)) n_part = n_tot + ne

      end subroutine eval_residual

      ! ------------------------------------------------------!

      subroutine newton_residual(Y, f_sp, Fvec)
      ! Thin wrapper: full residual, discarding the heat/cool diagnostics.
      real*8, dimension(3*N),                  intent(in)    :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(inout) :: f_sp
      real*8, dimension(3*N),                  intent(out)   :: Fvec
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      call eval_residual(Y, f_sp, Fvec, heat, cool)
      end subroutine newton_residual

      ! ------------------------------------------------------!

      subroutine frozen_residual(Y, n_part, heat, cool, Fvec)
      ! FROZEN-radiation residual: the hydro flux/gravity residual at Y with
      ! heat/cool held FIXED (no ioniz_eq, no column-density recompute). This
      ! is strictly local (WENO3 stencil) and hence exactly banded, so the
      ! colored finite-difference Jacobian of THIS residual is exact. It is
      ! the approximate Jacobian / preconditioner for the inexact Newton: the
      ! weakly non-local radiation response is left to the outer iteration.
      real*8, dimension(3*N),       intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: n_part, heat, cool
      real*8, dimension(3*N),       intent(out) :: Fvec
      real*8, dimension(3,1-Ng:N+Ng) :: u, R
      u = 0.0d0
      call unpack_U(Y, u)
      call Apply_BC(u)
      call assemble_residual(u, n_part, heat, cool, R)
      call pack_R(R, Fvec)
      if (nfix_base .gt. 0) call apply_base_fix(Y, Fvec)
      end subroutine frozen_residual

      ! ------------------------------------------------------!

      subroutine pack_R(R, Fvec)
      ! Pack residual array (physical cells) into the flat vector.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: R
      real*8, dimension(3*N),         intent(out) :: Fvec
      integer :: j, k
      do j = 1, N
         do k = 1, 3
            Fvec(3*(j-1)+k) = R(k,j)
         enddo
      enddo
      end subroutine pack_R

      ! ------------------------------------------------------!

      subroutine build_banded_jac(Y, n_part, heat, cool, ab)
      ! Colored finite-difference banded Jacobian of the FROZEN residual,
      ! stored in LAPACK general-band form for dgbtrf/dgbtrs:
      !   ab(kl+ku+1 + i - j, j) = J(i,j),  i in [j-ku, j+kl]
      ! ldab = 2*kl+ku+1. Columns sharing a color (stride ncolor = kl+ku+1)
      ! have disjoint row supports, so one probe per color recovers their
      ! band entries with no cross-contamination.
      real*8, dimension(3*N),               intent(in)  :: Y
      real*8, dimension(1-Ng:N+Ng),         intent(in)  :: n_part, heat, cool
      real*8, dimension(2*kl_jac+ku_jac+1, 3*N), intent(out) :: ab
      real*8, dimension(3*N) :: F0, Fp, Yp, dYc
      integer :: neq, color, jcol, irow, ilo, ihi
      real*8  :: sqeps

      neq   = 3*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0

      call frozen_residual(Y, n_part, heat, cool, F0)

      do color = 1, ncolor_jac
         Yp  = Y
         dYc = 0.0d0
         do jcol = color, neq, ncolor_jac
            dYc(jcol) = sqeps*max(abs(Y(jcol)), 1.0d0)
            Yp(jcol)  = Y(jcol) + dYc(jcol)
         enddo
         call frozen_residual(Yp, n_part, heat, cool, Fp)
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac

      ! ------------------------------------------------------!

      subroutine build_banded_jac_full(Y, f_sp_base, ab)
      ! Colored-FD banded Jacobian of the FULL residual (eval_residual, i.e.
      ! including the local ionization-equilibrium + heat/cool response). The
      ! frozen version omits the dominant local d(heat-cool)/dE coupling, so
      ! Newton has no descent direction even near the solution; this version
      ! captures it. The weakly non-local column-density coupling slightly
      ! contaminates same-color columns (treated as preconditioner error).
      ! Each probe restores f_sp from the base copy (ioniz_eq mutates it).
      real*8, dimension(3*N),                  intent(in) :: Y
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1, 3*N), intent(out) :: ab
      real*8, dimension(3*N) :: F0, Fp, Yp, dYc, Dsc
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      integer :: neq, color, jcol, irow, ilo, ihi
      real*8  :: sqeps

      neq   = 3*N
      sqeps = sqrt(epsilon(1.0d0))
      ab    = 0.0d0

      ! FD column steps RELATIVE to the scale for each unknown (a flat
      ! max(|Y|,1) floor gives wind cells, ~1e-7 of the base in code
      ! units, order-unity relative kicks and garbage columns).
      call cell_state_scales(Y, Dsc)

      fwork = f_sp_base
      call eval_residual(Y, fwork, F0, heat, cool)

      do color = 1, ncolor_jac
         Yp  = Y;  dYc = 0.0d0
         do jcol = color, neq, ncolor_jac
            dYc(jcol) = sqeps*Dsc(jcol)
            Yp(jcol)  = Y(jcol) + dYc(jcol)
         enddo
         fwork = f_sp_base
         call eval_residual(Yp, fwork, Fp, heat, cool)
         do jcol = color, neq, ncolor_jac
            ilo = max(1,   jcol - ku_jac)
            ihi = min(neq, jcol + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jcol, jcol) =                &
                    (Fp(irow) - F0(irow))/dYc(jcol)
            enddo
         enddo
      enddo

      end subroutine build_banded_jac_full

      ! ------------------------------------------------------!

      subroutine band_matvec(ab, x, y)
      ! y = J*x for the band-stored Jacobian ab (same layout as above).
      real*8, dimension(2*kl_jac+ku_jac+1, 3*N), intent(in)  :: ab
      real*8, dimension(3*N),                    intent(in)  :: x
      real*8, dimension(3*N),                    intent(out) :: y
      integer :: neq, irow, jcol, jlo, jhi
      neq = 3*N
      do irow = 1, neq
         y(irow) = 0.0d0
         jlo = max(1,   irow - kl_jac)
         jhi = min(neq, irow + ku_jac)
         do jcol = jlo, jhi
            y(irow) = y(irow) + ab(kl_jac+ku_jac+1 + irow - jcol, jcol)*x(jcol)
         enddo
      enddo
      end subroutine band_matvec

      ! ------------------------------------------------------!

      subroutine resid_relnorm(F, u, rc, rnorm)
      ! Relative residual rc(k) for each component and its max rnorm over the
      ! wind region [j_min:N] (same definition the marching monitor uses).
      real*8, dimension(3*N),         intent(in)  :: F
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8,                         intent(out) :: rnorm
      integer :: j, k
      real*8  :: fmx, umx, w
      do k = 1, 3
         fmx = 0.0d0;  umx = 0.0d0
         if (resid_vol) then
            ! volume-weighted (default): sum_j |F|V / sum_j |u|V, V = r^2 dr
            do j = j_min, N
               w   = r(j)*r(j)*dr_j(j)
               fmx = fmx + abs(F(3*(j-1)+k))*w
               umx = umx + abs(u(k,j))*w
            enddo
         else
            ! legacy L-inf: max over the wind region
            do j = j_min, N
               fmx = max(fmx, abs(F(3*(j-1)+k)))
               umx = max(umx, abs(u(k,j)))
            enddo
         endif
         rc(k) = fmx/max(umx, 1.0d-30)
      enddo
      rnorm = maxval(rc)
      end subroutine resid_relnorm

      ! ------------------------------------------------------!

      subroutine solve_steady_ptc(u, f_sp, resid_tol, maxit, dtau0, info)
      ! Pseudo-transient-continuation inexact Newton solve of the steady
      ! residual F(Y)=0. Per iteration:
      !   F = eval_residual(Y)                       (full, nonlocal radiation)
      !   J = build_banded_jac (frozen radiation)    (banded approx / precond)
      !   (I/dtau + J) dY = -F   via dgbtrf/dgbtrs    (LAPACK banded LU)
      !   Y <- Y + lam*dY        (backtracking line search + positivity)
      !   dtau <- dtau * rnorm_old/rnorm_new (SER ramp; cut on failure)
      ! As dtau->inf this is Newton; small dtau behaves like explicit
      ! relaxation, giving the robust startup PTC is designed for.
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j
      real*8  :: rnorm, rnorm_try, rc(3), dtau, lam, f2, f2_try
      logical :: ok

      neq  = 3*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      dtau = dtau0
      info = 1

      call pack_U(u, Y)
      call eval_residual(Y, f_sp, F, heat0, cool0)
      call resid_relnorm(F, u, rc, rnorm)
      f2 = sqrt(sum(F*F))                ! smooth line-search merit
      write(*,'(A,ES11.3,A,3ES10.2,A,ES10.2)') ' (PTC) start ||R||=',rnorm, &
           '  R(m,p,E)=',rc,'  ||F||2=',f2

      do iter = 1, maxit
         if (rnorm .lt. resid_tol) then
            info = 0;  exit
         endif

         ! Full-residual banded Jacobian (includes the local source
         ! derivative d(heat-cool)/dE), then the PTC system (I/dtau + J).
         call build_banded_jac_full(Y, f_sp, ab)
         abf = ab
         do jc = 1, neq
            abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + 1.0d0/dtau
         enddo
         dY = -F
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         if (lpinfo .ne. 0) then
            write(*,'(A,I0,A)') ' (PTC) dgbtrf info=',lpinfo,            &
                 ' (singular); cutting dtau'
            dtau = max(dtau*0.25d0, dtau0*1.0d-3);  cycle
         endif
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, ldab, ipiv,      &
                     dY, neq, lpinfo)

         ! Backtracking line search on the SMOOTH merit ||F||_2 (the
         ! Newton direction reduces this; the component-wise max-relative
         ! rnorm is non-smooth and rejected valid steps). Armijo condition
         ! with positivity (rho>0, p>0).
         lam = 1.0d0;  ok = .false.
         do ls = 1, 20
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               f_sp_j = f_sp
               call eval_residual(Ytry, f_sp_j, Ftry, heat0, cool0)
               f2_try = sqrt(sum(Ftry*Ftry))
               if (f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2) then
                  ok = .true.;  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            ! SER ramp on the merit ratio, scaled by the accepted step.
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            ! No descent found: shrink the pseudo-time step (more
            ! relaxation-like) but keep it >= the explicit-stable dtau0.
            dtau = max(dtau*0.25d0, dtau0)
         endif

         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES10.2,A,ES9.2)') ' (PTC) it', &
              iter,'  ||R||=',rnorm,'  ||F||2=',f2,'  dtau=',dtau,'  lam=',lam
      enddo

      call unpack_U(Y, u);  call Apply_BC(u)
      call resid_relnorm(F, u, rc, rnorm)
      write(*,'(A,I0,A,ES11.3)') ' (PTC) done info=',info,' ||R||=',rnorm

      deallocate(Y,F,Ftry,dY,Ytry,ab,abf,ipiv)
      end subroutine solve_steady_ptc

      ! ------------------------------------------------------!

      subroutine jv_product(Y, F0, f_sp_base, v, Jv)
      ! Matrix-free Jacobian-vector product J*v by a forward directional
      ! finite difference of the FULL residual (captures the non-local
      ! radiation coupling the banded Jacobian omits):
      !   J*v ~= ( F(Y + eps*v) - F0 ) / eps,   F0 = F(Y).
      ! eps is the standard scaled step. f_sp is restored from the base
      ! copy (eval_residual mutates it).
      real*8, dimension(3*N),                  intent(in)  :: Y, F0, v
      real*8, dimension(1-Ng:N+Ng,n_species),  intent(in)  :: f_sp_base
      real*8, dimension(3*N),                  intent(out) :: Jv
      real*8, dimension(1-Ng:N+Ng,n_species) :: fwork
      real*8, dimension(3*N) :: Fp
      real*8, dimension(1-Ng:N+Ng) :: heat, cool
      real*8 :: vn, eps
      vn = sqrt(sum(v*v))
      if (vn .le. 0.0d0) then
         Jv = 0.0d0;  return
      endif
      eps   = sqrt(epsilon(1.0d0))*(1.0d0 + sqrt(sum(Y*Y)))/vn
      fwork = f_sp_base
      call eval_residual(Y + eps*v, fwork, Fp, heat, cool)
      Jv = (Fp - F0)/eps
      end subroutine jv_product

      ! ------------------------------------------------------!

      subroutine pgmres(Y, F0, f_sp_base, D, abf, ipiv, idtau, b, x, m, &
                        rtol, gm_iters)
      ! Right-preconditioned GMRES(m), single cycle, for the JFNK step, in
      ! the DIAGONALLY SCALED space:
      !   solve  A_z x = b,   A_z v = v*idtau + D^-1 * J*(D v),
      ! (J*(Dv) matrix-free; b and x are scaled quantities; the caller maps
      ! dY = D*x back). Preconditioner M = factored banded SCALED system
      ! (I/dtau + D^-1 J_banded D) supplied in abf.
      real*8, dimension(3*N),  intent(in)  :: Y, F0, b, D
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_base
      real*8, dimension(2*kl_jac+ku_jac+1,3*N), intent(in) :: abf
      integer, dimension(3*N), intent(in)  :: ipiv
      real*8,  intent(in)  :: idtau, rtol
      integer, intent(in)  :: m
      real*8, dimension(3*N), intent(out) :: x
      integer, intent(out) :: gm_iters

      integer :: neq, i, j, kk, lpinfo
      real*8, allocatable :: V(:,:), Hs(:,:), gg(:), cs(:), sn(:), yy(:)
      real*8, allocatable :: z(:), w(:), u(:)
      real*8 :: beta, hij, nrm, denom, tmp

      neq = 3*N
      allocate(V(neq,m+1), Hs(m+1,m), gg(m+1), cs(m), sn(m), yy(m))
      allocate(z(neq), w(neq), u(neq))
      Hs = 0.0d0; gg = 0.0d0; cs = 0.0d0; sn = 0.0d0; x = 0.0d0

      ! r0 = b - A*0 = b
      beta = sqrt(sum(b*b))
      gm_iters = 0
      if (beta .le. 0.0d0) then
         deallocate(V,Hs,gg,cs,sn,yy,z,w,u);  return
      endif
      V(:,1) = b/beta
      gg(1)  = beta

      do j = 1, m
         gm_iters = j
         ! z = M^{-1} V(:,j)   (right preconditioning)
         z = V(:,j)
         call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, 2*kl_jac+ku_jac+1, &
                     ipiv, z, neq, lpinfo)
         ! w = A_z z = z*idtau + D^-1 * J*(D z)
         call jv_product(Y, F0, f_sp_base, D*z, w)
         w = w/D + idtau*z
         ! Modified Gram-Schmidt
         do i = 1, j
            Hs(i,j) = sum(w*V(:,i))
            w = w - Hs(i,j)*V(:,i)
         enddo
         Hs(j+1,j) = sqrt(sum(w*w))
         if (Hs(j+1,j) .gt. 0.0d0) V(:,j+1) = w/Hs(j+1,j)
         ! Apply previous Givens rotations to column j
         do i = 1, j-1
            tmp       =  cs(i)*Hs(i,j) + sn(i)*Hs(i+1,j)
            Hs(i+1,j) = -sn(i)*Hs(i,j) + cs(i)*Hs(i+1,j)
            Hs(i,j)   =  tmp
         enddo
         ! New Givens rotation to zero Hs(j+1,j)
         denom  = sqrt(Hs(j,j)**2 + Hs(j+1,j)**2)
         if (denom .le. 0.0d0) denom = 1.0d0
         cs(j)  = Hs(j,j)/denom
         sn(j)  = Hs(j+1,j)/denom
         Hs(j,j)   = cs(j)*Hs(j,j) + sn(j)*Hs(j+1,j)
         Hs(j+1,j) = 0.0d0
         gg(j+1) = -sn(j)*gg(j)
         gg(j)   =  cs(j)*gg(j)
         if (abs(gg(j+1)) .le. rtol*beta) exit
      enddo

      ! Back-substitute H(1:kk,1:kk) yy = gg(1:kk)
      kk = gm_iters
      do i = kk, 1, -1
         tmp = gg(i)
         do j = i+1, kk
            tmp = tmp - Hs(i,j)*yy(j)
         enddo
         yy(i) = tmp/Hs(i,i)
      enddo
      ! u = V(:,1:kk) yy ;  x = M^{-1} u  (undo right preconditioning)
      u = 0.0d0
      do i = 1, kk
         u = u + yy(i)*V(:,i)
      enddo
      call dgbtrs('N', neq, kl_jac, ku_jac, 1, abf, 2*kl_jac+ku_jac+1,  &
                  ipiv, u, neq, lpinfo)
      x = u

      deallocate(V,Hs,gg,cs,sn,yy,z,w,u)
      end subroutine pgmres

      ! ------------------------------------------------------!

      subroutine solve_steady_jfnk(u, f_sp, resid_tol, maxit, dtau0,    &
                                   gm_m, info)
      ! Pseudo-transient-continuation Jacobian-free Newton-Krylov solve.
      ! Per outer iteration:
      !   F  = eval_residual(Y)                          (full residual)
      !   M  = factored banded (I/dtau + J_banded)       (preconditioner)
      !   solve (I/dtau + J) dY = -F  by right-precond. GMRES (J*v matrix-
      !                                free; captures non-local radiation)
      !   Y <- Y + lam*dY   (||F||_2 line search + positivity)
      !   dtau <- SER ramp
      real*8, dimension(3,1-Ng:N+Ng),         intent(inout) :: u
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8,  intent(in)  :: resid_tol, dtau0
      integer, intent(in)  :: maxit, gm_m
      integer, intent(out) :: info

      integer :: neq, ldab, iter, ls, lpinfo, jc, gmit
      integer :: jj, kk, jworst, kworst, irow, ilo, ihi
      real*8  :: amx
      real*8, allocatable :: Y(:), F(:), Ftry(:), dY(:), Ytry(:)
      real*8, allocatable :: dZ(:), D(:)
      real*8, allocatable :: ab(:,:), abf(:,:)
      integer, allocatable :: ipiv(:)
      real*8, dimension(1-Ng:N+Ng)            :: heat0, cool0
      real*8, dimension(3,1-Ng:N+Ng)          :: utry, Wtry
      real*8, dimension(1-Ng:N+Ng,n_species)  :: f_sp_j, f_sp_best
      real*8  :: rnorm, rc(3), dtau, lam, f2, f2_try, idtau
      real*8  :: f2hist(5), f2ref, rnorm_best
      integer :: n_no_descent, n_window_used
      ! Stagnation limit: consecutive outer iterations in which the line
      ! search found NO acceptable step at all. Measured on the HD 209458 b
      ! hand-off state (docs/newton_scaling_and_base_wall.md): runs that go on
      ! to converge never string more than 4 such iterations together, runs
      ! that are truly stuck string 34 or more.
      integer, parameter :: n_no_descent_max = 12
      logical :: ok
      real*8, allocatable :: Ybest(:)

      neq  = 3*N
      ldab = 2*kl_jac + ku_jac + 1
      allocate(Y(neq), F(neq), Ftry(neq), dY(neq), Ytry(neq))
      allocate(dZ(neq), D(neq), Ybest(neq))
      allocate(ab(ldab,neq), abf(ldab,neq), ipiv(neq))
      dtau = dtau0;  info = 1

      call pack_U(u, Y)
      call eval_residual(Y, f_sp, F, heat0, cool0)
      call resid_relnorm(F, u, rc, rnorm)
      call cell_state_scales(Y, D)
      f2 = sqrt(sum((F/D)**2))     ! merit in the SCALED space
      f2hist = f2                  ! non-monotone line-search memory
      rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp;  n_no_descent = 0
      n_window_used = 0
      write(*,'(A,ES11.3,A,ES10.2)') ' (JFNK) start ||R||=',rnorm,        &
           '  ||Fs||2=',f2

      do iter = 1, maxit
         if (rnorm .lt. resid_tol) then
            info = 0;  exit
         endif
         idtau = 1.0d0/dtau

         ! Diagonal scaling for this outer iteration.
         call cell_state_scales(Y, D)

         ! FREEZE the WENO3 weights at the current iterate: one mode-1
         ! residual evaluation stores the smoothness factors and yields the
         ! frozen-weights F; the inner evaluations that build the Newton model
         ! (Jacobian probes, J*v) then reuse them (mode 2), so the inner
         ! problem excludes the strongly nonlinear weight response (the
         ! standard lagged-weights remedy for FV steady solves). The line
         ! search does NOT use them -- see there.
         weno_mode = 1
         f_sp_j = f_sp
         call eval_residual(Y, f_sp_j, F, heat0, cool0)
         weno_mode = 2
         f2 = sqrt(sum((F/D)**2))     ! TRUE merit of the current iterate

         ! Grippo non-monotone reference: the worst true merit of the last 5
         ! outer iterates, the current one included. Recording it here, and not
         ! only when a step is accepted, is what keeps the reference from
         ! lagging behind the state the line search actually starts from.
         f2hist = (/ f2hist(2:5), f2 /)

         ! Banded preconditioner of the SCALED system:
         ! M = I/dtau + D^-1 J_banded D, factored.
         call build_banded_jac_full(Y, f_sp, ab)
         do jc = 1, neq
            ilo = max(1,   jc - ku_jac)
            ihi = min(neq, jc + kl_jac)
            do irow = ilo, ihi
               ab(kl_jac+ku_jac+1 + irow - jc, jc) =                     &
                    ab(kl_jac+ku_jac+1 + irow - jc, jc)*D(jc)/D(irow)
            enddo
         enddo
         abf = ab
         do jc = 1, neq
            abf(kl_jac+ku_jac+1, jc) = abf(kl_jac+ku_jac+1, jc) + idtau
         enddo
         call dgbtrf(neq, neq, kl_jac, ku_jac, abf, ldab, ipiv, lpinfo)
         if (lpinfo .ne. 0) then
            dtau = max(dtau*0.25d0, dtau0);  cycle
         endif

         ! GMRES solve of the scaled PTC system; unscale the step.
         call pgmres(Y, F, f_sp, D, abf, ipiv, idtau, -F/D, dZ, gm_m,   &
                     1.0d-1, gmit)
         dY = D*dZ

         ! Scaled ||F/D||_2 backtracking line search with positivity, tested
         ! against the non-monotone reference above.
         !
         ! The trial states are evaluated with weno_mode = 0, i.e. with the
         ! smoothness weights recomputed at the trial state, so what decides
         ! acceptance is the residual the solve is driving to zero. Evaluating
         ! the trials with the weights frozen at Y instead -- the quantity the
         ! inner Newton model minimizes -- is a different measure: on the
         ! HD 189733 b hand-off state the two differ by a median factor 2.2 at
         ! lam = 1, and 7 of the 14 steps the frozen test accepted RAISED the
         ! true residual, by up to 9x (docs/newton_scaling_and_base_wall.md
         ! §10). Mode 0 leaves the stored weights alone, so the frozen model of
         ! this outer iteration survives the search.
         f2ref = maxval(f2hist)
         lam = 1.0d0;  ok = .false.
         weno_mode = 0
         do ls = 1, 20
            Ytry = Y + lam*dY
            call unpack_U(Ytry, utry)
            call U_to_W(utry, Wtry)
            if (minval(Wtry(1,1:N)) .gt. 0.0d0 .and.                    &
                minval(Wtry(3,1:N)) .gt. 0.0d0) then
               f_sp_j = f_sp
               call eval_residual(Ytry, f_sp_j, Ftry, heat0, cool0)
               f2_try = sqrt(sum((Ftry/D)**2))
               if (f2_try .lt. (1.0d0 - 1.0d-4*lam)*f2ref) then
                  ok = .true.
                  if (f2_try .ge. (1.0d0 - 1.0d-4*lam)*f2)              &
                       n_window_used = n_window_used + 1
                  exit
               endif
            endif
            lam = 0.5d0*lam
         enddo

         if (ok) then
            Y = Ytry;  F = Ftry;  f_sp = f_sp_j
            call unpack_U(Y, u);  call Apply_BC(u)
            call resid_relnorm(F, u, rc, rnorm)
            dtau = min(dtau*max(lam,0.1d0)*(f2/max(f2_try,1.0d-30)),     &
                       1.0d14*dtau0)
            f2 = f2_try
         else
            dtau = max(dtau*0.25d0, dtau0)
         endif

         ! Keep the best iterate seen (by the convergence measure) so that a
         ! failed solve returns it rather than wherever it stopped.
         if (rnorm .lt. rnorm_best) then
            rnorm_best = rnorm;  Ybest = Y;  f_sp_best = f_sp
         endif

         ! Stagnation: the failure mode of the non-monotone search is that it
         ! stops finding ANY acceptable step and the iterate no longer moves.
         ! Count consecutive such iterations. (A watchdog on rnorm instead
         ! cannot work here: rnorm measures only [j_min:N], while the solver
         ! minimizes the whole-domain merit, and the opening pseudo-transient
         ! legitimately raises rnorm for ~30 iterations while it repairs the
         ! sub-sonic region -- see docs/newton_scaling_and_base_wall.md.)
         if (ok) then
            n_no_descent = 0
         else
            n_no_descent = n_no_descent + 1
         endif
         if (n_no_descent .ge. n_no_descent_max) then
            info = 2
            write(*,'(A,I0,A)') ' (JFNK) line search found no descent '// &
                 'step in ', n_no_descent_max, ' consecutive iterations'// &
                 ' -- aborting'
            exit
         endif

         ! Locate the cell carrying the largest scaled residual |F/D|, i.e.
         ! the term the merit is actually dominated by, over the whole domain
         ! 1..N. (Normalizing by max_j|u(k,j)| instead reports the base cell
         ! almost always, because the base holds the global maximum of |rho v|
         ! while carrying no wind.)
         amx = -1.0d0;  jworst = 1;  kworst = 1
         do jj = 1, N
            do kk = 1, 3
               if (abs(F(3*(jj-1)+kk))/D(3*(jj-1)+kk) .gt. amx) then
                  amx = abs(F(3*(jj-1)+kk))/D(3*(jj-1)+kk)
                  jworst = jj;  kworst = kk
               endif
            enddo
         enddo
         write(*,'(A,I4,A,ES11.3,A,ES10.2,A,ES9.2,A,I3,A,I4,A,F7.3,A,I1)') &
              ' (JFNK) it',iter,'  ||R||=',rnorm,'  ||Fs||2=',f2,        &
              '  lam=',lam,'  gm=',gmit,'  worst j=',jworst,            &
              ' r=',r(jworst),' k=',kworst
      enddo

      weno_mode = 0                 ! restore default reconstruction
      if (info .ne. 0 .and. rnorm_best .lt. rnorm) then
         Y = Ybest;  f_sp = f_sp_best;  rnorm = rnorm_best
         write(*,'(A,ES11.3)') ' (JFNK) returning best iterate, '//      &
              '||R||=', rnorm
      endif
      call unpack_U(Y, u);  call Apply_BC(u)
      write(*,'(A,I0,A,ES11.3,A,I0)') ' (JFNK) done info=',info,        &
           ' ||R||=',rnorm,'  window-only accepts=',n_window_used
      deallocate(Y,F,Ftry,dY,Ytry,dZ,D,Ybest,ab,abf,ipiv)
      end subroutine solve_steady_jfnk

      ! End of module
      end module steady_newton
