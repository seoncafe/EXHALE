      module Reconstruction_step
      ! Collection of reconstruction procedure subroutine
      
      use global_parameters
      use grid_construction, only: spherical_cell_volume
      use BC_Apply
      use PLM_reconstruction
      use Conversion

      implicit none

      ! Stored WENO3 smoothness factors for the frozen-weights mode
      ! (weno_mode in global_parameters; 0 = off/default, byte-identical).
      real*8, allocatable :: S0sav(:,:), S1sav(:,:)

      ! Number of FACE STATES whose reconstruction left the admissible
      ! thermodynamic state (rho > 0, p > 0) and was scaled back toward its own
      ! cell average, summed over the whole run. Counts one per face state, not
      ! one per face pair: the limiter scales the offending side only.
      ! Reported at the end of a run; zero means the high-order reconstruction
      ! was admissible everywhere and the run is the one an unguarded build
      ! would have produced, to the bit.
      integer :: n_faces_positivity_limited = 0

      ! ---- THE WELL-BALANCED OPTION ("Well balanced:", default off) ----
      ! The local hydrostatic equilibrium of cell j is the one of constant
      ! density through its own (rho_j, p_j),
      !     p_eq,j(r) = p_j - rho_j (phi(r) - phi(r_j)),
      ! whose two face values are wb_P_up(j) at r_edg(j) and wb_P_dn(j) at
      ! r_edg(j-1) (Kaeppeli and Mishra 2016, A&A 587, A94, their eq. 16, with
      ! the potential taken at the true face radius since ours is analytic).
      ! What the reconstruction returns is the DEPARTURE from it,
      !     wb_dev_L(f) = p of the left state of face f  - wb_P_up(f)
      !     wb_dev_R(f) = p of the right state of face f - wb_P_dn(f+1),
      ! and wb_dp_eq(f) = wb_P_dn(f+1) - wb_P_up(f) is the mismatch of the two
      ! neighboring equilibria at the shared face, formed from the difference
      ! of the two cell pressures and two terms of the size of the hydrostatic
      ! pressure drop across a cell. It vanishes on the discrete equilibrium
      ! the scheme preserves. The pressure jump of the Riemann problem is
      ! wb_dp_eq + wb_dev_R - wb_dev_L, every term of the size of the
      ! departure and none of the size of the state, which is the whole point:
      ! the flux then carries the last bit of its own value and not of an O(1)
      ! pressure. Read by RK_rhs; written here at every evaluation and never
      ! carried between them.
      real*8, allocatable :: wb_P_up(:), wb_P_dn(:)
      real*8, allocatable :: wb_p_cell(:), wb_rho_cell(:)
      real*8, allocatable :: wb_dev_L(:), wb_dev_R(:), wb_dp_eq(:)
      ! The face pressures the equilibrium/departure reconstruction assembled,
      ! kept so that a face state the boundary condition or the positivity
      ! limiter has rewritten is recognized (its departure is then re-formed
      ! from the state it now carries, at the cost of the O(1) subtraction the
      ! rest of that option avoids).
      real*8, allocatable :: wb_pL_asm(:), wb_pR_asm(:)

      contains

      subroutine Reconstruct(u_in,WL_out,WR_out) 
      
      integer :: j,k
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng) :: WL,WR
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out

      ! WENO3 variables
      real*8, dimension(3,1-Ng:N+Ng) :: W,dW
      real*8, dimension(3) :: dWp,dWm
      real*8, dimension(3) :: b0,b1
      real*8, dimension(3) :: tau
      real*8, dimension(3) :: S0,S1
      real*8, dimension(1-Ng:N+Ng) :: C1,C2,D1,D2

      ! Well-balanced option: the cell's own equilibrium is formed first, the
      ! reconstruction below then works on the departure from it.
      real*8, dimension(3) :: Wc

      if (well_balanced) call hydrostatic_equilibrium_of_each_cell(u_in)

      select case (rec_method)
      
         !-----------------------------------------------!
         
         case ('PLM')      ! Piecewise linear reconstruction

            if (well_balanced) then
               call PLM_rec(u_in,WL,WR,wb_P_up,wb_P_dn,wb_dev_L,wb_dev_R)
            else
               call PLM_rec(u_in,WL,WR)
            endif

         case ('WENO3')    ! Third-order WENO reconstruction
            !
            ! THE NONLINEAR WEIGHTS.  The smoothness factors S0, S1 below
            ! are the weight functions of Yamaleev and Carpenter (2009,
            ! J. Comput. Phys. 228, 3025; cited in full at the end of this
            ! block), their Eqs. (18), (21) and (22): the weights are
            ! w_r = alpha_r / sum_l alpha_l  with
            ! alpha_r = d_r (1 + tau/(eps + beta_r)), where their tau is
            ! ALREADY a square, tau = (u_{j+1} - 2 u_j + u_{j-1})^2 (their
            ! Eq. 22), and beta_0 = (u_{j+1} - u_j)^2,
            ! beta_1 = (u_j - u_{j-1})^2 (their Eq. 20).  S0 and S1 are
            ! that 1 + tau/(eps + beta_r) with the square written out and
            ! eps folded into b0, b1.  Their ideal weights d_0 = 2/3,
            ! d_1 = 1/3 (Eq. 20) enter only through the ratio
            ! d_1/d_0 = 1/2, which is what D1 and D2 carry: exactly 1/2 on
            ! a grid uniform in the volume coordinate, 1/2 + O(h) on a
            ! stretched one.  The property that makes the reconstruction
            ! third order is their Eq. (24), w_r = d_r + O(Dx^2) on smooth
            ! data, so in a smooth region the scheme is its linear
            ! third-order limit and the stencil coefficients below are what
            ! set the order.
            !
            ! THE FLOOR eps.  Their prescription is eps = O(Dxi^2) (Eq. 62,
            ! made scale invariant as Eq. 64 with Dxi = 1/J the
            ! computational spacing), with the coefficient carrying the
            ! SOLUTION scale,
            !     eps = max(||u_0^2||, ||(u_0)_xi^2||) Dxi^2      (Eq. 65)
            ! so that eps and beta_r are the same kind of quantity.  Their
            ! Eq. (51) requires eps >= O(Dx^2) to hold the design order at
            ! a smooth extremum, where eps = 0 gives tau/beta_r = O(1) and
            ! the scheme drops to second order (Eq. 57); their Eq. (61)
            ! requires eps <= O(Dx^2) so that the stencil astride a
            ! discontinuity is still nullified.  Eq. (62) is the two
            ! together.
            !
            ! The floor used here, eps = dr_j^2, has that h^2 scaling and
            ! reads "O(Dx^2)" as the LOCAL cell width, which is the reading
            ! a grid whose spacing varies by a factor of 20 across the
            ! domain admits and which the published uniform-grid analysis
            ! does not distinguish.  What it drops is the solution scale of
            ! Eq. (65).  It is not dimensionally inconsistent (r is in R0
            ! and rho, v, p in n0, v0, p0, so both sides are
            ! dimensionless), but it is not invariant under a rescaling of
            ! the solution, which Eq. (65) is written to secure.  Two
            ! measurements bound what that costs.  In a smooth region it
            ! costs nothing: with eps = (kappa dr_j)^2 the face error and
            ! the convergence rate on the production grid are flat for
            ! kappa from 1e-1 to 1e3 (L_inf 2.99e-8 at N = 800, rate 2.991)
            ! and deteriorate only below kappa near 1e-2
            ! (weno3_reconstruction_order, MEASURED).  Where the choice
            ! bites is the ENO stencil biasing: on the WASP-121b reference
            ! output the density ratio beta_r/eps runs from 8.4e3 at the
            ! base to 2.2e-6 at the top of the domain, and a third of the
            ! cells sit below 1e-2, where S0 and S1 are within a percent of
            ! 1 and no biasing is left.  Eq. (65) with the base density as
            ! the scale would put a single eps = 4e-6 on the whole grid.
            !
            ! WHAT THIS IS NOT.  The ESWENO scheme of that paper is its
            ! Eq. (9) with P, Q and S of Eqs. (31) and (35) to (37)
            ! TOGETHER with these weights; the energy stability comes from
            ! the extra dissipation flux of Eq. (48), built from mu_j of
            ! Eq. (37) and a second parameter d, and that term is not
            ! computed here.  The paper separates the two itself: "Note
            ! that these weights can also be used with the conventional
            ! third-order WENO scheme" (Section 3.1).  So what runs here is
            ! a third-order WENO reconstruction with their weights, and it
            ! inherits their accuracy result (Eq. 24) and not the energy
            ! estimate of their Section 4.1.  Two further reasons that
            ! estimate does not carry over: their operator is a
            ! finite-difference reconstruction of the FLUX on a uniform
            ! grid (their Section 3.1; neither paper treats a non-uniform
            ! grid), while this is a finite-volume reconstruction of the
            ! primitive STATE handed to a Riemann solver on a stretched
            ! radial grid, so the volume-coordinate coefficients below are
            ! this code's own; and for a system their estimate holds for
            ! the characteristic form, "while it is not clear how to obtain
            ! a similar energy estimate when the system of equations is
            ! discretized in the component-wise fashion" (their Section 5),
            ! which is what the loop below does, on primitive variables at
            ! that.
            !
            !   N. K. Yamaleev and M. H. Carpenter, "Third-order Energy
            !   Stable WENO scheme", J. Comput. Phys. 228, 3025 (2009),
            !   doi:10.1016/j.jcp.2009.01.011.  Every equation number
            !   above is from that paper, read in the published version.
            !   Its companion, "A systematic methodology for constructing
            !   high-order energy stable WENO schemes", J. Comput. Phys.
            !   228, 4248 (2009), doi:10.1016/j.jcp.2009.03.002, carries
            !   the same weight form as its Eq. (58) and the same
            !   solution-scaled floor as its Eq. (79), for design orders
            !   four and above; the third-order scheme used here is the
            !   3025 paper's.
			
            ! Convert to primitive variables
            call U_to_W(u_in,W)
         
            ! Evaluate jumps at interfaces
            dW(:,1-Ng:N+Ng-1) = W(:,2-Ng:N+Ng) - W(:,1-Ng:N+Ng-1)
            dW(:,N+Ng) = 0.0
            
            ! Evaluate reconstruction geometry-dependent coefficients
            call weno3_geometry_coefficients(C1,C2,D1,D2)
            
            ! Construct smoothness indicators and reconstructed values.
            ! weno_mode: 0 = fresh (default), 1 = fresh + store, 2 = reuse
            ! stored factors (frozen weights for the steady Newton solver).
            if (weno_mode .gt. 0 .and. .not. allocated(S0sav)) then
               allocate(S0sav(0:N+1,3), S1sav(0:N+1,3))
               S0sav = 1.0d0;  S1sav = 1.0d0
            endif
            ! Cell-local: the reconstruction of cell j reads the two
            ! interface jumps that bound it and writes its own two face
            ! states, WL at face j and WR at face j-1.  No reduction, so
            ! the arithmetic of a cell is unchanged and the result is
            ! bitwise identical at any number of threads.
            !$omp parallel do default(shared) schedule(static)          &
            !$omp   private(j,k,dWp,dWm,b0,b1,tau,S0,S1,Wc)
            do j = 0,N+1

               dWp = dW(:,j)
               dWm = dW(:,j-1)
               Wc  = W(:,j)

               ! WELL-BALANCED OPTION: the pressure component is reconstructed
               ! in the coordinate of the DEPARTURE from cell j's local
               ! hydrostatic equilibrium, whose value at the cell centre is
               ! zero and at a neighboring centre is that cell's pressure
               ! measured against the equilibrium continued through the face
               ! between them, with cell j's density up to the face and the
               ! neighbor's beyond it (Kaeppeli and Mishra 2016, A&A 587,
               ! A94, their eqs. 18 and 26).  The data then vanishes exactly
               ! on the discrete equilibrium and the stencil, the smoothness
               ! indicators and the volume shares act on it unchanged.
               if (well_balanced) then
                  dWp(3) = dWp(3) + Wc(1)   *(Gphi_i(j)   - Gphi_c(j))   &
                                  + W(1,j+1)*(Gphi_c(j+1) - Gphi_i(j))
                  dWm(3) = dWm(3) + Wc(1)   *(Gphi_c(j)   - Gphi_i(j-1)) &
                                  + W(1,j-1)*(Gphi_i(j-1) - Gphi_c(j-1))
                  Wc(3)  = 0.0d0
               endif

               if (weno_mode .eq. 2) then
                  S0 = S0sav(j,:)
                  S1 = S1sav(j,:)
               else
                  do k = 1,3
                     b0(k) = dWp(k)*dWp(k) + dr_j(j)*dr_j(j)
                     b1(k) = dWm(k)*dWm(k) + dr_j(j)*dr_j(j)
                  enddo
                  tau = dWp - dWm
                  S0 = 1.0 + tau*tau/b0
                  S1 = 1.0 + tau*tau/b1
                  if (weno_mode .eq. 1) then
                     S0sav(j,:) = S0
                     S1sav(j,:) = S1
                  endif
               endif

               ! Each candidate is the linear interpolant of two cell
               ! averages in the volume coordinate, evaluated at a face of
               ! cell j, which lies dV(j)/2 from the centre of cell j: the
               ! jump to the neighbor is weighted by cell j's OWN volume
               ! share, dV(j)/(dV(j) + dV(j+1)) = C2(j) for the jump to j+1
               ! and dV(j)/(dV(j-1) + dV(j)) = C1(j-1) for the jump to j-1,
               ! at either face. With the neighbor's share in their place
               ! the face value carries an O(h) coefficient error times an
               ! O(h) jump and the reconstruction is second order, so these
               ! two shares and not the nonlinear weights are what sets the
               ! order (MEASURED, weno3_reconstruction_order: rate 3.000 on
               ! a grid uniform in r and 2.991 on the production grid with
               ! the shares as written, 1.997 and 1.988 with the two
               ! neighbor shares in their place).
               WL(:,j) = Wc	&
                  + (S0*C2(j)*dWp + D1(j)*S1*C1(j-1)*dWm) &
                  /(S0 + D1(j)*S1)

               WR(:,j-1) = Wc  &
                  - (D2(j)*S0*C2(j)*dWp + S1*C1(j-1)*dWm) &
                  /(D2(j)*S0 + S1)

               ! The departure the stencil returned, kept as its own small
               ! number, and the face pressure it belongs to.
               if (well_balanced) then
                  wb_dev_L(j)   = WL(3,j)
                  wb_dev_R(j-1) = WR(3,j-1)
                  WL(3,j)       = wb_P_up(j) + wb_dev_L(j)
                  WR(3,j-1)     = wb_P_dn(j) + wb_dev_R(j-1)
               endif
            enddo
            !$omp end parallel do

         case default

            write(*,*) 'ERROR: unknown reconstruction scheme: ', trim(rec_method)
            write(*,*) '  allowed: PLM, WENO3'
            error stop 1

      end select

      ! The face pressures the equilibrium/departure reconstruction
      ! assembled, rebuilt from their two parts rather than read out of
      ! WL/WR: it is the same number to the bit where the reconstruction
      ! wrote one (the same sum of the same two operands), and at the
      ! outermost faces, which no reconstruction loop writes and which
      ! Rec_BC states below, the departure is still zero and the value here
      ! is the equilibrium alone, so those faces are recognized as rewritten
      ! instead of being compared against an undefined double.
      if (well_balanced) then
         do j = 1-Ng,N+Ng
            wb_pL_asm(j) = wb_P_up(j) + wb_dev_L(j)
            wb_pR_asm(j) = wb_P_dn(min(j+1,N+Ng)) + wb_dev_R(j)
         enddo
      endif

      ! Apply BC to reconstructed variables
      call Rec_BC(WL,WR,WL_out,WR_out)

      ! Last step before the Riemann solver sees the face states: restore
      ! rho > 0 and p > 0 wherever the reconstruction lost them.
      call positivity_limited_faces(u_in,WL_out,WR_out)

      ! Well-balanced option: the equilibrium mismatch of every face, and the
      ! departures of the face states the two steps above rewrote.
      if (well_balanced) call well_balanced_face_departures(WL_out,WR_out)

      ! End of subroutine
      end subroutine Reconstruct

      !-----------------------------------------------!

      subroutine hydrostatic_equilibrium_of_each_cell(u_in)
      ! The two face values of the local hydrostatic equilibrium of every
      ! cell: constant density rho_j through the cell's own pressure p_j,
      !     p_eq,j(r) = p_j - rho_j (phi(r) - phi(r_j)),
      ! evaluated at r_edg(j) and at r_edg(j-1) (Kaeppeli and Mishra 2016,
      ! A&A 587, A94, eq. 16).  It is the mechanical balance
      ! dp/dr = -rho dphi/dr to second order in the cell width and states
      ! nothing about the temperature or the entropy, which is why it
      ! preserves an arbitrary stratification.
      !
      ! The cell averages of the pressure and the density are kept as well:
      ! the equilibrium mismatch of a face is formed from THEM and not from
      ! the two assembled face pressures, so that the difference of the two
      ! cell pressures (exact in binary floating point while neighboring
      ! cells lie within a factor of two of each other) is not first buried
      ! in a sum with an O(1) number.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng) :: W
      integer :: j

      if (.not. allocated(wb_P_up)) then
         allocate(wb_P_up(1-Ng:N+Ng),   wb_P_dn(1-Ng:N+Ng),             &
                  wb_p_cell(1-Ng:N+Ng), wb_rho_cell(1-Ng:N+Ng),         &
                  wb_dev_L(1-Ng:N+Ng),  wb_dev_R(1-Ng:N+Ng),            &
                  wb_dp_eq(1-Ng:N+Ng),                                  &
                  wb_pL_asm(1-Ng:N+Ng), wb_pR_asm(1-Ng:N+Ng))
      endif

      call U_to_W(u_in,W)

      wb_p_cell   = W(3,:)
      wb_rho_cell = W(1,:)

      do j = 1-Ng,N+Ng
         wb_P_up(j) = W(3,j) - W(1,j)*(Gphi_i(j) - Gphi_c(j))
      enddo
      do j = 2-Ng,N+Ng
         wb_P_dn(j) = W(3,j) + W(1,j)*(Gphi_c(j) - Gphi_i(j-1))
      enddo
      ! The lowest ghost cell has no lower face (RK_rhs forms no row for it),
      ! and Gphi_i does not reach below r_edg(1-Ng).
      wb_P_dn(1-Ng) = W(3,1-Ng)

      wb_dev_L  = 0.0d0
      wb_dev_R  = 0.0d0
      wb_dp_eq  = 0.0d0

      end subroutine hydrostatic_equilibrium_of_each_cell

      !-----------------------------------------------!

      subroutine well_balanced_face_departures(WL,WR)
      ! Finish the well-balanced face data: the equilibrium mismatch of every
      ! face,
      !   wb_dp_eq(f) = P_dn(f+1) - P_up(f)
      !               = (p_f+1 - p_f) + rho_f+1 (phi_c(f+1) - phi_i(f))
      !                               + rho_f   (phi_i(f)   - phi_c(f)),
      ! formed in that order so that it is a sum of quantities of the size of
      ! the hydrostatic pressure drop across a cell and vanishes on the
      ! discrete equilibrium; and the departure of any face state the
      ! boundary condition or the positivity limiter rewrote, which has to be
      ! re-formed by subtracting the equilibrium from the O(1) state it now
      ! carries.  The outermost face has no cell to its right and takes the
      ! same cell on both sides, as the flux does.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL,WR
      integer :: j,jr

      do j = 1-Ng,N+Ng

         jr = min(j+1,N+Ng)

         if (jr .gt. j) then
            wb_dp_eq(j) = (wb_p_cell(jr) - wb_p_cell(j))                 &
                        + wb_rho_cell(jr)*(Gphi_c(jr) - Gphi_i(j))       &
                        + wb_rho_cell(j) *(Gphi_i(j)  - Gphi_c(j))
         else
            wb_dp_eq(j) = wb_P_dn(jr) - wb_P_up(j)
         endif

         if (WL(3,j) .ne. wb_pL_asm(j)) wb_dev_L(j) = WL(3,j) - wb_P_up(j)
         if (WR(3,j) .ne. wb_pR_asm(j)) wb_dev_R(j) = WR(3,j) - wb_P_dn(jr)

      enddo

      end subroutine well_balanced_face_departures

      !-----------------------------------------------!

      subroutine weno3_geometry_coefficients(C1,C2,D1,D2)
      ! Grid-dependent coefficients of the WENO3 reconstruction: the two
      ! volume-weighted interpolation weights C1, C2 of each face and the two
      ! stencil weights D1, D2 of each cell.  They depend on the grid alone,
      ! so the scalar reconstruction of a mass fraction and the
      ! reconstruction of the primitive vector take them from here and cannot
      ! disagree about the grid.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: C1,C2,D1,D2
      real*8, dimension(1-Ng:N+Ng) :: dV
      integer :: j

      ! The shell volumes, from the one expression for them
      ! (grid_construction), over the range the weights below read: the
      ! stencils reach one cell past the physical column at each end, and
      ! the two outermost ghost cells have no stencil of their own and take
      ! the volume of their neighbour.  Only ratios enter the weights below,
      ! so the factor 3 between the shell volume and the difference of cubes
      ! is immaterial here; the volume is what the divergence rows divide by,
      ! and one grid has one volume.
      do j = 0, N+1			
         dV(j) = spherical_cell_volume(j)
      enddo
      dV(1-Ng) = dV(2-Ng)
      dV(N+Ng) = dV(N+Ng-1)

      do j = 1-Ng,N+Ng-1
         C1(j) = dV(j+1)/(dV(j) + dV(j+1))
         C2(j) = 1.0 - C1(j)
      enddo

      do j = 2-Ng,N+Ng-1  	
         D1(j) = dV(j+1)/(dV(j) + dV(j-1))
         D2(j) = dV(j-1)/(dV(j) + dV(j+1))
      enddo

      ! End of subroutine
      end subroutine weno3_geometry_coefficients

      !-----------------------------------------------!

      subroutine Reconstruct_scalar(q,qL,qR)
      ! Reconstruct ONE cell-centered scalar to the two sides of every face
      ! with the reconstruction the run uses for the primitive variables, so
      ! that a species mass fraction and the density whose face flux carries
      ! it are extrapolated to a face by the same scheme.
      !
      ! qL(j) is the state on the left of face r_edg(j) and qR(j) the state
      ! on its right, the pairing WL/WR use.
      !
      ! No boundary reconstruction of its own: the ghost values the caller
      ! supplies are what the outermost faces extrapolate from, which for a
      ! mass fraction is the imposed base composition inside and the
      ! zero-gradient copy outside.  A face state is left unbounded here and
      ! is limited by the caller against the composition simplex, because
      ! [0,1] is a property of a mass fraction and not of a reconstruction.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: q
      real*8, dimension(1-Ng:N+Ng), intent(out) :: qL,qR

      real*8, dimension(1-Ng:N+Ng) :: dq,C1,C2,D1,D2
      real*8  :: dqp,dqm,b0,b1,tau,S0,S1
      integer :: j

      select case (rec_method)

         case ('PLM')

            call PLM_rec_scalar(q,qL,qR)

         case ('WENO3')

            ! Same WENO3 weights and same face formula as the primitive
            ! reconstruction in Reconstruct above, on one component; keep
            ! the two in step.  The frozen-weight mode of the steady solver
            ! is not read here: the stored factors belong to the primitive
            ! variables, and a mass fraction has smoothness indicators of
            ! its own.
            dq(1-Ng:N+Ng-1) = q(2-Ng:N+Ng) - q(1-Ng:N+Ng-1)
            dq(N+Ng) = 0.0

            call weno3_geometry_coefficients(C1,C2,D1,D2)

            do j = 1-Ng,N+Ng
               qL(j) = q(j)
               qR(j) = q(min(j+1,N+Ng))
            enddo

            do j = 0,N+1

               dqp = dq(j)
               dqm = dq(j-1)

               b0  = dqp*dqp + dr_j(j)*dr_j(j)
               b1  = dqm*dqm + dr_j(j)*dr_j(j)
               tau = dqp - dqm
               S0  = 1.0 + tau*tau/b0
               S1  = 1.0 + tau*tau/b1

               ! Cell j's own volume share on each jump (see Reconstruct).
               qL(j) = q(j)                                              &
                  + (S0*C2(j)*dqp + D1(j)*S1*C1(j-1)*dqm)                &
                  /(S0 + D1(j)*S1)

               qR(j-1) = q(j)                                            &
                  - (D2(j)*S0*C2(j)*dqp + S1*C1(j-1)*dqm)                &
                  /(D2(j)*S0 + S1)
            enddo

         case default

            write(*,*) 'ERROR: unknown reconstruction scheme: ',          &
                       trim(rec_method)
            error stop 1

      end select

      ! End of subroutine
      end subroutine Reconstruct_scalar

      !-----------------------------------------------!

      subroutine positivity_limited_faces(u_in,WL,WR)

      ! Scale each face state CONTINUOUSLY toward its own cell average, by the
      ! largest factor that keeps rho > 0 and p > 0 there.
      !
      !     W_face  <-  W_avg + theta ( W_rec - W_avg ),   theta in [0,1]
      !     theta   =   min over rho and p of  (q_avg - eps)/(q_avg - q_rec)
      !
      ! ATTRIBUTION, and where this departs from what it cites. The idea is
      ! the linear scaling limiter of Zhang and Shu (2010, J. Comput. Phys.,
      ! 229, 8918; doi:10.1016/j.jcp.2010.08.016; the publisher PDF is
      ! references/Zhang_2010JCP_229_8918.pdf), their Section 2.2. This is
      ! NOT their limiter, and the differences matter:
      !
      !   * They scale the CONSERVED vector w = (rho, m, E), in two steps:
      !     density first about the cell average, their Eqs. (2.4) and (2.5)
      !     with theta_1 = min{(rho_avg - eps)/(rho_avg - rho_min), 1}, then
      !     the whole vector, their Eqs. (2.10) and (2.11), with theta_2 the
      !     root t_eps of the QUADRATIC in their Eq. (2.12),
      !     p[(1 - t) w_avg + t w_rec] = eps. This routine scales the
      !     PRIMITIVE vector (rho, v, p) in one step, so its theta is the
      !     root of a linear equation in p and not of their quadratic.
      !     The two are not the same number. Along their path rho, m and E
      !     are linear in t and p(t) = (gamma-1)(E - m^2/(2 rho)) is CONCAVE
      !     in t, so p(t) is never below the straight line between p_avg and
      !     p_rec that this routine interpolates: at equal eps this theta is
      !     the smaller of the two, that is the more restrictive one, never
      !     the less safe one. Measured at a Mach 60 state with eps = 1e-3
      !     p_avg: theta = 0.85629 here against t_eps = 0.94910 from their
      !     Eq. (2.12) (src/tests/physics_probe/positivity_limiter_scaling.f90).
      !
      !   * Their theta is one number per CELL, the minimum over that cell's
      !     quadrature points, which is what leaves the cell average of the
      !     limited polynomial equal to w_avg (the conservativity property of
      !     their Section 2.2). Here theta is one number per FACE STATE, so a
      !     cell whose two ends are scaled by different factors no longer
      !     reconstructs to its own average. Measured departure on a cell
      !     with both ends limited: 8.3e-2 of the cell average, against
      !     round-off if the cell's minimum theta is used for both ends.
      !
      !   * Their floor is absolute, eps = 10^-13 in their computations and
      !     eps = min_j {10^-13, rho_avg, p(w_avg)} in their implementation
      !     flowchart. Here it is relative, one unit in the last place of the
      !     cell average.
      !
      ! What carries over is their Lemma 2.5: the limited value is a convex
      ! combination of an admissible cell average with the reconstruction, so
      ! it stays admissible, and no quadrature enters that argument. Their
      ! Theorem 2.1 does NOT carry over. It concludes that the next CELL
      ! AVERAGE is admissible, and its hypotheses are the N-point Legendre
      ! Gauss-Lobatto set of their Eq. (1.7) with 2N - 3 >= k (a
      ! finite-volume reconstruction to two face values supplies only N = 2,
      ! which covers the k = 1 of PLM and not the k = 2 of WENO3; their
      ! Remark 2.6 names finite-volume WENO as the open case) and the CFL
      ! condition of their Eq. (2.1), which nothing here imposes.
      !
      ! Pressure positivity of the SCALED state is not at issue in this form:
      ! p is itself one of the scaled variables, so the limited face carries
      ! the pressure the scaling put there and its internal energy
      ! p/(gamma-1) has the same sign. The nonlinear map that forces their
      ! quadratic has no counterpart here.
      !
      ! WHY IT IS CONTINUOUS AND WHY THAT IS THE POINT (docs/Update_EXHALE_stage1.md
      ! section 138; the measurement is P49). This routine used to be a hard
      ! switch: as soon as a reconstructed rho or p crossed zero, the WHOLE
      ! face pair was replaced by the two cell averages. That is a STEP
      ! DISCONTINUITY of the residual F(Y), and it is what stopped the
      ! molecular reload's steady solve. Measured on the hot Uranus hand-off:
      ! the WENO3 left density at the face of cell 273 (r = 1.2324 R_p) sat at
      ! 2.50e-15 against cell averages of 1e-4, eleven orders below and
      ! positive by a hair, so the iterate sat exactly on the switching
      ! surface. Raising the neighboring cell's density by one part in 1e9
      ! tipped it through zero, the face density jumped to 4.17e-6, the mass
      ! residual of cell 273 jumped by a factor 97, and the merit jumped from
      ! 202.5 to 945.1, BY THE SAME FACTOR 4.67 at every step size from 1e-2
      ! to 1e-9, which is a jump and not a slope. Every direction the solver
      ! could build put its largest scaled component at that face, so all of
      ! them ascended: the Newton/PTC step at every damping over six decades
      ! and all ten damped Gauss-Newton steps, by factors 4.7 to 27.8. With
      ! the scaling above the face value crosses zero continuously instead,
      ! and the merit has a slope there rather than a step.
      !
      ! THE FLOOR eps. The switch this replaces tested q > 0, so its floor was
      ! zero; eps = 0 here would scale the face value to EXACTLY zero and hand
      ! the HLLC solver sqrt(gamma p/0). The floor used is therefore the
      ! smallest one that is not a chosen number: one unit in the last place of
      ! the cell average, eps = epsilon(1.0d0)*q_avg. It introduces no
      ! dimensional constant and no tuning, and the residual jump that survives
      ! at the crossing is 2.2e-16 of the cell average, round-off, against
      ! the factor 1.7e9 the hard switch made.
      !
      ! THE FLOOR IS DELIVERED BY A CLAMP AFTER THE UPDATE. eps is one unit
      ! in the last place of q_avg, and q_avg + theta (q_rec - q_avg) is then
      ! a cancellation whose result is decided by the rounding of the product
      ! rather than by theta: MEASURED over 200000 reconstructed values that
      ! cross zero, the bare update is EXACTLY ZERO for 1013 of them and never
      ! negative (src/tests/physics_probe/positivity_limiter_scaling.f90).
      ! Zero is the value this floor exists to avoid: a zero face density
      ! divides in v = m/rho and in the sound speed of Num_Fluxes.f90. The
      ! scaled density and pressure are therefore clamped to eps after the
      ! update, q <- max(q_avg + theta (q_rec - q_avg), eps), the same number
      ! in exact arithmetic (2026-09-07; moves only runs in which the limiter
      ! fired, and those by one ulp except where the zero would have divided).
      !
      ! BYTE-IDENTITY. A face whose reconstruction is admissible has theta = 1
      ! and is NOT rewritten (W_avg + 1*(W - W_avg) is not bitwise W), so every
      ! run in which the guard never fired is unchanged to the bit. Runs in
      ! which it did fire move, and by more than the floor: the old switch
      ! replaced both states of a face and this scales only the offending one.
      !
      ! Validity of the reconstruction, and why the limiter is needed. PLM and
      ! WENO3 both extrapolate the PRIMITIVE variables (rho, v, p) to a face
      ! with slopes taken from the neighboring cells, which is a valid
      ! approximation only while the solution varies smoothly across the
      ! stencil. It has no positivity property of its own: a face value can
      ! cross zero while every cell average on the stencil is positive. The
      ! state that does it here is a hypersonic layer, where the pressure the
      ! scheme recovers as p = (gamma-1)(E - rho v^2/2) is the difference of two
      ! nearly equal numbers: at Mach 60 the thermal pressure is 2e-4 of the
      ! total energy density, so a reconstruction across the neighboring jump
      ! tips it negative and the HLLC sound speed sqrt(gamma p/rho) takes the
      ! square root of it (Num_Fluxes.f90; the abort of TO_BE_DONE item (O)).
      !
      ! The cell averages are the only states the update guarantees admissible,
      ! so the repair is built from them. Scaling the whole increment keeps the
      ! face state a convex combination of the reconstruction and an admissible
      ! state, so it is the same scheme at reduced order, not a change of the
      ! equations.
      !
      ! Face j takes its left state from cell j and its right state from cell
      ! j+1, which is the pairing the two loops below use. The outermost face
      ! has no cell j+1, so it falls back on the last cell average, the
      ! zero-gradient state the outer BC reconstructs to in any case.
      !
      ! If a CELL AVERAGE is itself non-positive the state is unphysical before
      ! any reconstruction and nothing here can repair it; that face is left as
      ! it is, and the NaN detector of the marching loop reports it.

      real*8, dimension(3,1-Ng:N+Ng), intent(in)    :: u_in
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: WL,WR
      real*8, dimension(3,1-Ng:N+Ng) :: W_avg
      logical :: any_bad
      integer :: j,jr
      real*8  :: th

      ! Scan first: the conversion to cell-average primitives is only needed
      ! when a face has actually left rho > 0, p > 0, which on an admissible
      ! run is never. A face value inside (0, eps] would be missed by this
      ! test, and is left alone deliberately: its theta would be 1 - 2.2e-16.
      ! The tests are written as the NEGATION of "strictly positive" so that a
      ! NaN face state, which compares false against everything, is caught
      ! too; it can only come from the reconstruction arithmetic when the cell
      ! averages below are finite, and the same repair applies.
      any_bad = .false.
      do j = 1-Ng,N+Ng
         if (.not. (WL(1,j) .gt. 0.0d0 .and. WL(3,j) .gt. 0.0d0 .and.   &
                    WR(1,j) .gt. 0.0d0 .and. WR(3,j) .gt. 0.0d0)) then
            any_bad = .true.
            exit
         endif
      enddo
      if (.not. any_bad) return

      call U_to_W(u_in,W_avg)

      do j = 1-Ng,N+Ng
         jr = min(j+1, N+Ng)
         th = positivity_scaling(WL(1,j), WL(3,j), W_avg(1,j), W_avg(3,j))
         if (th .lt. 1.0d0) then
            WL(:,j) = W_avg(:,j) + th*(WL(:,j) - W_avg(:,j))
            ! The update is a cancellation of the size of the floor itself
            ! (theta puts the variable ON the floor, eps = one ulp of the
            ! average), so its rounding can land on zero; the clamp delivers
            ! the floor the scaling was solved for. Same number in exact
            ! arithmetic.
            WL(1,j) = max(WL(1,j), epsilon(1.0d0)*W_avg(1,j))
            WL(3,j) = max(WL(3,j), epsilon(1.0d0)*W_avg(3,j))
            n_faces_positivity_limited = n_faces_positivity_limited + 1
         endif
         th = positivity_scaling(WR(1,j), WR(3,j), W_avg(1,jr), W_avg(3,jr))
         if (th .lt. 1.0d0) then
            WR(:,j) = W_avg(:,jr) + th*(WR(:,j) - W_avg(:,jr))
            WR(1,j) = max(WR(1,j), epsilon(1.0d0)*W_avg(1,jr))
            WR(3,j) = max(WR(3,j), epsilon(1.0d0)*W_avg(3,jr))
            n_faces_positivity_limited = n_faces_positivity_limited + 1
         endif
      enddo

      ! End of subroutine
      end subroutine positivity_limited_faces


      !-----------------------------------------------!

      double precision function positivity_scaling(rho_f, p_f, rho_a, p_a)  &
                                                                 result(th)
      ! The largest theta in [0,1] for which rho and p of
      ! W_avg + theta (W_face - W_avg) both stay at or above their floors in
      ! exact arithmetic; see the header for what round-off then delivers.
      ! theta = 1 whenever the face state is already admissible, so the caller
      ! can leave such a face untouched and keep it bitwise unchanged.
      real*8, intent(in) :: rho_f, p_f, rho_a, p_a
      th = 1.0d0
      ! A cell average that is not itself admissible cannot be the anchor of a
      ! convex combination; that face is left as it is (see the header).
      if (.not. (rho_a .gt. 0.0d0 .and. p_a .gt. 0.0d0)) return
      th = min(th, positive_variable_scaling(rho_f, rho_a))
      th = min(th, positive_variable_scaling(p_f,   p_a))
      end function positivity_scaling

      !-----------------------------------------------!

      double precision function positive_variable_scaling(q_f, q_a) result(th)
      ! One variable's share of the scaling: the theta that puts
      ! q_a + theta (q_f - q_a) on the floor eps = epsilon*q_a when the
      ! reconstruction went below it, and 1 when it did not. Exactly on it in
      ! exact arithmetic only: eps is one ulp of q_a, so the caller's update
      ! is a cancellation of the size of the floor itself (header). q_a > 0 is
      ! the caller's precondition.
      real*8, intent(in) :: q_f, q_a
      real*8 :: qeps
      th = 1.0d0
      ! A face value that is not a number carries no scaling; the cell average
      ! is the only admissible state left, which is theta = 0.
      if (q_f .ne. q_f) then
         th = 0.0d0
         return
      endif
      qeps = epsilon(1.0d0)*q_a
      if (q_f .ge. qeps) return
      th = (q_a - qeps)/(q_a - q_f)
      th = max(0.0d0, min(1.0d0, th))
      end function positive_variable_scaling


      ! End of module
      end module Reconstruction_step
