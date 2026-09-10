!  ONE SOURCE TEXT, TWO REAL KINDS: the hydrodynamic rows of the finite-volume
!  operator
!
!      R = dF - S ,   dF = (A_+ F_+ - A_- F_-)/dV ,
!
!  from the conserved state through the reconstruction, the Riemann flux and
!  the geometric and gravitational sources, written once in the kind fv and
!  instantiated at kind(1.0d0) and at the quadruple kind by the two modules
!  that include this file.  The double instantiation reproduces
!  Reconstruction_step/Numerical_Fluxes/source_func/RK_integration bit for bit
!  (asserted by the suite row generic_double_pipeline_equals_the_operator);
!  the quadruple instantiation is the control experiment of the residual's
!  non-smoothness floor.
!
!  Every routine here mirrors a routine of those modules statement by
!  statement.  The physics, the references and the reasons for each formula
!  are stated at the original; this file states only what the two texts must
!  keep in common.  Two texts of one formula can drift apart, which is why the
!  bitwise assertion above is a test row and not a remark.
!
!  The composition is NOT of this kind: the mixture ratios and the H2
!  rovibrational level table are properties of the gas that the ionization
!  sweep produced in double, and they enter here converted exactly.  What is
!  evaluated in fv is the arithmetic of the operator on that composition.

      ! ------------------------------------------------------!
      ! The caloric equation of state, at the composition caloric_eos holds
      ! ------------------------------------------------------!

      subroutine h2_rovibrational_energy_and_heat_capacity_fv(TK, u_rv, c_rv)
      ! Mean rovibrational energy u_rv(T) [K] of one H2 molecule and its
      ! derivative c_rv = du_rv/dT, from the cubic Hermite interpolant of the
      ! table caloric_eos builds from the observed bound ladder of H2
      ! X^1 Sigma_g^+ (Roueff et al. 2019, A&A 630, A58, table 2).  The nodes
      ! are the double table; the interpolation is in fv.
      real(fv), intent(in)  :: TK
      real(fv), intent(out) :: u_rv, c_rv
      real*8  :: lnT_lo8, dlnT8, T_lo8, T_hi8, u_i8, du_i8, u_ip8, du_ip8
      integer :: n_nodes, i
      real(fv) :: lnT_lo, dlnT, T_lo, T_hi, u_i, du_i, u_ip, du_ip
      real(fv) :: x, s, s2, s3, h00, h10, h01, h11
      real(fv) :: dh00, dh10, dh01, dh11

      if (caloric_eos_monatomic) then
         u_rv = 0.0_fv;  c_rv = 0.0_fv
         return
      endif

      call h2_rovibrational_table_grid(n_nodes, lnT_lo8, dlnT8, T_lo8, T_hi8)
      lnT_lo = real(lnT_lo8, fv);  dlnT = real(dlnT8, fv)
      T_lo   = real(T_lo8,   fv);  T_hi = real(T_hi8, fv)

      if (TK .le. T_lo) then
         call h2_rovibrational_table_node(1, u_i8, du_i8)
         c_rv = real(du_i8, fv)/T_lo
         u_rv = real(u_i8, fv) + c_rv*(TK - T_lo)
         if (u_rv .lt. 0.0_fv) u_rv = 0.0_fv
         return
      endif
      if (TK .ge. T_hi) then
         call h2_rovibrational_table_node(n_nodes, u_i8, du_i8)
         c_rv = real(du_i8, fv)/T_hi
         u_rv = real(u_i8, fv) + c_rv*(TK - T_hi)
         return
      endif

      x = (log(TK) - lnT_lo)/dlnT
      i = 1 + int(x)
      if (i .lt. 1)            i = 1
      if (i .gt. n_nodes - 1)  i = n_nodes - 1
      s  = x - real(i-1, fv)
      s2 = s*s;  s3 = s2*s

      h00 =  2.0_fv*s3 - 3.0_fv*s2 + 1.0_fv
      h10 =         s3 - 2.0_fv*s2 + s
      h01 = -2.0_fv*s3 + 3.0_fv*s2
      h11 =         s3 -        s2
      dh00 =  6.0_fv*s2 - 6.0_fv*s
      dh10 =  3.0_fv*s2 - 4.0_fv*s + 1.0_fv
      dh01 = -6.0_fv*s2 + 6.0_fv*s
      dh11 =  3.0_fv*s2 - 2.0_fv*s

      call h2_rovibrational_table_node(i,   u_i8,  du_i8)
      call h2_rovibrational_table_node(i+1, u_ip8, du_ip8)
      u_i  = real(u_i8, fv);   du_i  = real(du_i8, fv)
      u_ip = real(u_ip8, fv);  du_ip = real(du_ip8, fv)

      u_rv = h00*u_i  + h10*dlnT*du_i                                    &
           + h01*u_ip + h11*dlnT*du_ip
      c_rv = ( dh00*u_i  + dh10*dlnT*du_i                                &
             + dh01*u_ip + dh11*dlnT*du_ip )                             &
             /(dlnT*TK)

      end subroutine h2_rovibrational_energy_and_heat_capacity_fv

      ! ------------------------------------------------------!

      function internal_energy_of_mixture_fv(x2, T) result(e_out)
      ! Internal energy per (n_tot + n_e) in units of k T0.
      real(fv), intent(in) :: x2, T
      real(fv) :: e_out, u_rv, c_rv

      if (x2 .le. 0.0_fv) then
         e_out = T/(real(gamma_ad, fv) - 1.0_fv);  return
      endif
      call h2_rovibrational_energy_and_heat_capacity_fv(T*real(T0,fv),   &
                                                        u_rv, c_rv)
      e_out = 1.5_fv*T + x2*u_rv/real(T0, fv)

      end function internal_energy_of_mixture_fv

      ! ------------------------------------------------------!

      function heat_capacity_of_mixture_fv(x2, T) result(cv_out)
      ! C_V/((n_tot + n_e) k) = 1/(gamma_eff - 1).
      real(fv), intent(in) :: x2, T
      real(fv) :: cv_out, u_rv, c_rv

      if (x2 .le. 0.0_fv) then
         cv_out = 1.0_fv/(real(gamma_ad, fv) - 1.0_fv);  return
      endif
      call h2_rovibrational_energy_and_heat_capacity_fv(T*real(T0,fv),   &
                                                        u_rv, c_rv)
      cv_out = 1.5_fv + x2*c_rv

      end function heat_capacity_of_mixture_fv

      ! ------------------------------------------------------!

      function temperature_of_mixture_fv(x2, u_in) result(T_out)
      ! Invert internal_energy_of_mixture_fv for T: f(T) = 1.5 T
      ! + x2 u_rv(T T0)/T0 - u_in is strictly increasing and 0 < T <= u_in/1.5
      ! brackets its root, so Newton safeguarded by bisection on that bracket
      ! converges from the upper end.  The stopping tolerance and the
      ! iteration limit are the kind's, so the root is resolved to the
      ! precision the arm is run at (fv_newton_tol, fv_newton_maxit).
      real(fv), intent(in) :: x2, u_in
      real(fv) :: T_out, Tlo, Thi, T, Tnew, f, df, u_rv, c_rv
      integer  :: it

      if (x2 .le. 0.0_fv .or. .not. (u_in .gt. 0.0_fv)) then
         T_out = u_in*(real(gamma_ad, fv) - 1.0_fv)
         return
      endif

      Tlo = 0.0_fv
      Thi = u_in/1.5_fv
      T   = Thi
      do it = 1, fv_newton_maxit
         call h2_rovibrational_energy_and_heat_capacity_fv(              &
                 T*real(T0,fv), u_rv, c_rv)
         f  = 1.5_fv*T + x2*u_rv/real(T0,fv) - u_in
         df = 1.5_fv + x2*c_rv
         ! An exact zero IS the root; the bracket update would file it as an
         ! endpoint and the safeguard would bisect it away.  The equality
         ! test on a real is deliberate and is not a tolerance in disguise.
         if (f .eq. 0.0_fv) exit
         if (f .gt. 0.0_fv) then
            Thi = T
         else
            Tlo = T
         endif
         Tnew = T - f/df
         if (Tnew .le. Tlo .or. Tnew .ge. Thi) Tnew = 0.5_fv*(Tlo + Thi)
         if (abs(Tnew - T) .le. fv_newton_tol*T) then
            T = Tnew
            exit
         endif
         T = Tnew
      enddo
      T_out = T

      end function temperature_of_mixture_fv

      ! ------------------------------------------------------!

      function pressure_from_energy_density_fv(j, rho_c, rho_e) result(p_out)
      ! p from the internal energy density of a state whose composition is
      ! that of cell j.
      integer,  intent(in) :: j
      real(fv), intent(in) :: rho_c, rho_e
      real(fv) :: p_out, nk, T
      real*8   :: nk_per_rho8, x28
      logical  :: is_mol

      call caloric_cell_mixture(j, nk_per_rho8, x28, is_mol)
      if (.not. is_mol .or. .not. (rho_c .gt. 0.0_fv)) then
         p_out = (real(gamma_ad, fv) - 1.0_fv)*rho_e;  return
      endif
      nk    = rho_c*real(nk_per_rho8, fv)
      T     = temperature_of_mixture_fv(real(x28, fv), rho_e/nk)
      p_out = nk*T

      end function pressure_from_energy_density_fv

      ! ------------------------------------------------------!

      function energy_density_from_pressure_fv(j, rho_c, p) result(e_out)
      integer,  intent(in) :: j
      real(fv), intent(in) :: rho_c, p
      real(fv) :: e_out, nk
      real*8   :: nk_per_rho8, x28
      logical  :: is_mol

      call caloric_cell_mixture(j, nk_per_rho8, x28, is_mol)
      if (.not. is_mol .or. .not. (rho_c .gt. 0.0_fv)) then
         e_out = p/(real(gamma_ad, fv) - 1.0_fv);  return
      endif
      nk    = rho_c*real(nk_per_rho8, fv)
      e_out = nk*internal_energy_of_mixture_fv(real(x28, fv), p/nk)

      end function energy_density_from_pressure_fv

      ! ------------------------------------------------------!

      function adiabatic_index_from_state_fv(j, rho_c, p) result(gam)
      ! gamma_eff = 1 + (n_tot + n_e) k/C_V of a state (rho_c, p) at the
      ! composition of cell j; the sound speeds are sqrt(gamma_eff p/rho).
      integer,  intent(in) :: j
      real(fv), intent(in) :: rho_c, p
      real(fv) :: gam
      real*8   :: nk_per_rho8, x28
      logical  :: is_mol

      call caloric_cell_mixture(j, nk_per_rho8, x28, is_mol)
      if (.not. is_mol .or. .not. (rho_c .gt. 0.0_fv)) then
         gam = real(gamma_ad, fv);  return
      endif
      gam = 1.0_fv + 1.0_fv/heat_capacity_of_mixture_fv(real(x28, fv),   &
               p/(rho_c*real(nk_per_rho8, fv)))

      end function adiabatic_index_from_state_fv

      ! ------------------------------------------------------!
      ! Conservative <-> primitive
      ! ------------------------------------------------------!

      subroutine U_to_W_comp_fv(U_in, W_out, jcell)
      real(fv), intent(in)  :: U_in(3)
      integer,  intent(in)  :: jcell
      real(fv), intent(out) :: W_out(3)

      W_out(1) = U_in(1)
      W_out(2) = U_in(2)/U_in(1)
      W_out(3) = pressure_from_energy_density_fv(jcell, U_in(1),         &
                 U_in(3)-0.5_fv*U_in(2)*U_in(2)/U_in(1))

      end subroutine U_to_W_comp_fv

      ! ------------------------------------------------------!

      subroutine W_to_U_comp_fv(W_in, U_out, jcell)
      real(fv), intent(in)  :: W_in(3)
      integer,  intent(in)  :: jcell
      real(fv), intent(out) :: U_out(3)

      U_out(1) = W_in(1)
      U_out(2) = W_in(1)*W_in(2)
      U_out(3) = 0.5_fv*W_in(1)*W_in(2)**2.0_fv                          &
                 + energy_density_from_pressure_fv(jcell, W_in(1), W_in(3))

      end subroutine W_to_U_comp_fv

      ! ------------------------------------------------------!

      subroutine U_to_W_fv(U_in, W_out)
      real(fv), intent(in)  :: U_in(3,1-Ng:N+Ng)
      real(fv), intent(out) :: W_out(3,1-Ng:N+Ng)
      integer :: j

      W_out(1,:) = U_in(1,:)
      W_out(2,:) = U_in(2,:)/U_in(1,:)
      !$omp parallel do default(shared) private(j) schedule(static)
      do j = 1-Ng, N+Ng
         W_out(3,j) = pressure_from_energy_density_fv(j, U_in(1,j),      &
            U_in(3,j)-0.5_fv*U_in(2,j)*U_in(2,j)/U_in(1,j))
      enddo
      !$omp end parallel do

      end subroutine U_to_W_fv

      ! ------------------------------------------------------!
      ! Reconstruction
      ! ------------------------------------------------------!

      function minmod_mc_slope_fv(a, b, c) result(d)
      ! Generalized minmod (monotonized central) limited slope, theta = 2
      ! (Kappeli 2016): zero unless the three difference quotients share a
      ! sign, which is what keeps a face value between the cell averages.
      real(fv), intent(in) :: a, b, c
      real(fv) :: d, v_arg(3)
      real(fv), parameter :: theta = 2.0_fv

      v_arg(1) = theta*a
      v_arg(2) = theta*b
      v_arg(3) = c

      if (maxval(v_arg).lt.(0.0_fv)) then
         d = maxval(v_arg)
      elseif (minval(v_arg).gt.(0.0_fv)) then
         d = minval(v_arg)
      else
         d = 0.0_fv
      endif

      end function minmod_mc_slope_fv

      ! ------------------------------------------------------!

      subroutine piecewise_linear_faces_fv(u_in, WL_rec, WR_rec,         &
                                           P_up, P_dn, dev_L, dev_R)
      ! WELL-BALANCED ARM.  With the four optional arrays present the
      ! PRESSURE is reconstructed in the coordinate of the departure from the
      ! cell's own local hydrostatic equilibrium; the value at the cell's own
      ! centre is then identically zero and the equilibrium face value is
      ! added back at the end (Kaeppeli and Mishra 2016, A&A 587, A94,
      ! section 2.1.3).  Mirrors PLM_rec.
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(out) :: WL_rec, WR_rec
      real(fv), dimension(1-Ng:N+Ng), intent(in),  optional :: P_up, P_dn
      real(fv), dimension(1-Ng:N+Ng), intent(out), optional :: dev_L, dev_R
      integer  :: j, k
      real(fv) :: x(-1:1), W(3,-1:1)
      real(fv), dimension(3) :: sp, sm, sd, sc
      logical  :: wb

      wb = present(P_up) .and. present(P_dn) .and.                       &
           present(dev_L) .and. present(dev_R)

      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,k,x,W,sp,sm,sd,sc)
      do j = 2-Ng,N+Ng-1

         x = real(r(j-1:j+1), fv)

         do k = j-1,j+1
            call U_to_W_comp_fv(u_in(:,k),W(:,k-j),k)
         enddo

         ! The departure at a neighboring cell centre is measured against
         ! the equilibrium continued to it through the face, that is with
         ! THIS cell's density up to the face and the NEIGHBOUR's beyond it,
         ! which is what makes the data vanish exactly on the discrete
         ! equilibrium the scheme preserves (Kaeppeli and Mishra 2016, A&A
         ! 587, A94, their eqs. 18 and 26).
         ! The difference of the two cell pressures is taken FIRST, while it
         ! is still an exact floating-point operation, and the hydrostatic
         ! terms are added to it afterwards; the reverse order would bury the
         ! small terms in the rounding of an O(1) sum.
         if (wb) then
            W(3,1)  = (W(3,1)  - W(3,0))                                 &
                    + W(1,0)*(real(Gphi_i(j),fv)   - real(Gphi_c(j),fv))  &
                    + W(1,1)*(real(Gphi_c(j+1),fv) - real(Gphi_i(j),fv))
            W(3,-1) = (W(3,-1) - W(3,0))                                 &
                    - W(1,0) *(real(Gphi_c(j),fv) - real(Gphi_i(j-1),fv)) &
                    - W(1,-1)*(real(Gphi_i(j-1),fv)                       &
                             - real(Gphi_c(j-1),fv))
            W(3,0)  = 0.0_fv
         endif

         sp = (W(:,1) - W(:,0))/(x(1)-x(0))
         sm = (W(:,0) - W(:,-1))/(x(0)-x(-1))
         sd = (W(:,1) - W(:,-1))/(x(1)-x(-1))

         do k = 1,3
            sc(k) = minmod_mc_slope_fv(sp(k),sm(k),sd(k))
         enddo

         WL_rec(:,j)   = W(:,0) + 0.5_fv*sc*(x(1)-x(0))
         WR_rec(:,j-1) = W(:,0) - 0.5_fv*sc*(x(0)-x(-1))

         ! The departure the limited slope returned, and the face pressure
         ! it belongs to
         if (wb) then
            dev_L(j)      = WL_rec(3,j)
            dev_R(j-1)    = WR_rec(3,j-1)
            WL_rec(3,j)   = P_up(j) + dev_L(j)
            WR_rec(3,j-1) = P_dn(j) + dev_R(j-1)
         endif

      enddo
      !$omp end parallel do

      end subroutine piecewise_linear_faces_fv

      ! ------------------------------------------------------!

      subroutine weno3_geometry_coefficients_fv(C1, C2, D1, D2)
      ! Volume-weighted interpolation weights and stencil weights of the
      ! WENO3 reconstruction; a function of the grid alone.
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: C1, C2, D1, D2
      real(fv), dimension(1-Ng:N+Ng) :: dV
      real(fv) :: rm, rp
      integer  :: j

      do j = 0, N+1
         rp = real(r_edg(j),   fv)
         rm = real(r_edg(j-1), fv)
         dV(j) = (rp*rp*rp - rm*rm*rm)
      enddo
      dV(1-Ng) = dV(2-Ng)
      dV(N+Ng) = dV(N+Ng-1)

      do j = 1-Ng,N+Ng-1
         C1(j) = dV(j+1)/(dV(j) + dV(j+1))
         C2(j) = 1.0_fv - C1(j)
      enddo

      do j = 2-Ng,N+Ng-1
         D1(j) = dV(j+1)/(dV(j) + dV(j-1))
         D2(j) = dV(j-1)/(dV(j) + dV(j+1))
      enddo

      end subroutine weno3_geometry_coefficients_fv

      ! ------------------------------------------------------!

      subroutine weno3_faces_fv(u_in, WL, WR, P_up, P_dn, dev_L, dev_R)
      ! Third-order WENO reconstruction with the nonlinear weights of
      ! Yamaleev and Carpenter (2009, J. Comput. Phys. 228, 3025), eqs. (18),
      ! (20) to (22), on the primitive variables; the floor of the smoothness
      ! indicators is dr_j^2.  weno_mode 0 fresh, 1 fresh and stored, 2 the
      ! stored factors (frozen weights).  The stored factors are the double
      ! arrays of Reconstruction_step, so the two schemes cannot disagree
      ! about which weights are frozen.
      !
      ! WELL-BALANCED ARM.  With the four optional arrays present the
      ! pressure component is reconstructed in the coordinate of the
      ! departure from cell j's local hydrostatic equilibrium, whose value at
      ! the cell centre is zero and at a neighboring centre is that cell's
      ! pressure measured against the equilibrium continued through the face
      ! between them, with cell j's density up to the face and the
      ! neighbor's beyond it (Kaeppeli and Mishra 2016, A&A 587, A94, their
      ! eqs. 18 and 26).  Mirrors the WENO3 branch of Reconstruct.
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(out) :: WL, WR
      real(fv), dimension(1-Ng:N+Ng), intent(in),  optional :: P_up, P_dn
      real(fv), dimension(1-Ng:N+Ng), intent(out), optional :: dev_L, dev_R
      real(fv), dimension(3,1-Ng:N+Ng) :: W, dW
      real(fv), dimension(1-Ng:N+Ng)   :: C1, C2, D1, D2
      real(fv), dimension(3) :: dWp, dWm, b0, b1, tau, S0, S1, Wc
      real(fv) :: drc
      integer  :: j, k
      logical  :: wb

      wb = present(P_up) .and. present(P_dn) .and.                       &
           present(dev_L) .and. present(dev_R)

      call U_to_W_fv(u_in, W)

      dW(:,1-Ng:N+Ng-1) = W(:,2-Ng:N+Ng) - W(:,1-Ng:N+Ng-1)
      dW(:,N+Ng) = 0.0_fv

      call weno3_geometry_coefficients_fv(C1,C2,D1,D2)

      if (weno_mode .gt. 0 .and. .not. allocated(S0sav)) then
         allocate(S0sav(0:N+1,3), S1sav(0:N+1,3))
         S0sav = 1.0d0;  S1sav = 1.0d0
      endif

      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,k,dWp,dWm,b0,b1,tau,S0,S1,drc,Wc)
      do j = 0,N+1

         dWp = dW(:,j)
         dWm = dW(:,j-1)
         Wc  = W(:,j)

         if (wb) then
            dWp(3) = dWp(3)                                              &
                   + Wc(1)   *(real(Gphi_i(j),fv)   - real(Gphi_c(j),fv))&
                   + W(1,j+1)*(real(Gphi_c(j+1),fv) - real(Gphi_i(j),fv))
            dWm(3) = dWm(3)                                              &
                   + Wc(1)   *(real(Gphi_c(j),fv) - real(Gphi_i(j-1),fv))&
                   + W(1,j-1)*(real(Gphi_i(j-1),fv)                      &
                             - real(Gphi_c(j-1),fv))
            Wc(3)  = 0.0_fv
         endif

         if (weno_mode .eq. 2) then
            S0 = real(S0sav(j,:), fv)
            S1 = real(S1sav(j,:), fv)
         else
            drc = real(dr_j(j), fv)
            do k = 1,3
               b0(k) = dWp(k)*dWp(k) + drc*drc
               b1(k) = dWm(k)*dWm(k) + drc*drc
            enddo
            tau = dWp - dWm
            S0 = 1.0_fv + tau*tau/b0
            S1 = 1.0_fv + tau*tau/b1
            if (weno_mode .eq. 1) then
               S0sav(j,:) = real(S0, kind(1.0d0))
               S1sav(j,:) = real(S1, kind(1.0d0))
            endif
         endif

         WL(:,j) = Wc                                                    &
            + (S0*C2(j)*dWp + D1(j)*S1*C1(j-1)*dWm)                      &
            /(S0 + D1(j)*S1)

         WR(:,j-1) = Wc                                                  &
            - (D2(j)*S0*C2(j)*dWp + S1*C1(j-1)*dWm)                      &
            /(D2(j)*S0 + S1)

         ! The departure the stencil returned, kept as its own small
         ! number, and the face pressure it belongs to.
         if (wb) then
            dev_L(j)   = WL(3,j)
            dev_R(j-1) = WR(3,j-1)
            WL(3,j)    = P_up(j) + dev_L(j)
            WR(3,j-1)  = P_dn(j) + dev_R(j-1)
         endif
      enddo
      !$omp end parallel do

      end subroutine weno3_faces_fv

      ! ------------------------------------------------------!

      subroutine reconstructed_face_boundary_states_fv(plm, WL_in, WR_in, &
                                                       WL_out, WR_out)
      ! Boundary conditions of the reconstructed variables.  The left state
      ! of the base face IS the boundary condition; the outer faces are the
      ! zero-gradient continuation, extrapolated linearly for WENO3.
      logical, intent(in) :: plm
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: WL_in, WR_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out
      integer :: k

      WL_out = WL_in
      WR_out = WR_in

      WL_out(:,0) = real(base_face_W, fv)
      do k = 2,Ng
         WL_out(:,1-k) = real(base_face_lower_W, fv)
         WR_out(:,1-k) = real(base_face_lower_W, fv)
      enddo

      do k = 1,Ng
         WR_out(:,N+k) = WR_out(:,N)
         if (.not. plm) WR_out(:,N+k) = 2.0_fv*WR_out(:,N+k-1)           &
                                      - WR_out(:,N+k-2)
      enddo

      WL_out(:,N+2) = WL_out(:,N+1)
      if (.not. plm) WL_out(:,N+2) = 2.0_fv*WL_out(:,N+1) - WL_out(:,N)

      end subroutine reconstructed_face_boundary_states_fv

      ! ------------------------------------------------------!

      function positive_variable_scaling_fv(q_f, q_a) result(th)
      ! The theta that puts q_a + theta (q_f - q_a) on the floor
      ! eps = epsilon q_a when the reconstruction went below it, 1 when it
      ! did not.  The floor is one unit in the last place OF THE KIND THE
      ! RECONSTRUCTION WAS EVALUATED IN: it exists to keep the face value off
      ! zero, and zero is what the cancellation of that kind produces.
      real(fv), intent(in) :: q_f, q_a
      real(fv) :: th, qeps
      th = 1.0_fv
      if (q_f .ne. q_f) then
         th = 0.0_fv
         return
      endif
      qeps = epsilon(1.0_fv)*q_a
      if (q_f .ge. qeps) return
      th = (q_a - qeps)/(q_a - q_f)
      th = max(0.0_fv, min(1.0_fv, th))
      end function positive_variable_scaling_fv

      ! ------------------------------------------------------!

      function positivity_scaling_fv(rho_f, p_f, rho_a, p_a) result(th)
      real(fv), intent(in) :: rho_f, p_f, rho_a, p_a
      real(fv) :: th
      th = 1.0_fv
      if (.not. (rho_a .gt. 0.0_fv .and. p_a .gt. 0.0_fv)) return
      th = min(th, positive_variable_scaling_fv(rho_f, rho_a))
      th = min(th, positive_variable_scaling_fv(p_f,   p_a))
      end function positivity_scaling_fv

      ! ------------------------------------------------------!

      subroutine positivity_limited_faces_fv(u_in, WL, WR, n_scaled)
      ! Scale each face state continuously toward its own cell average by the
      ! largest factor that keeps rho > 0 and p > 0 there (the linear scaling
      ! limiter of Zhang and Shu 2010, J. Comput. Phys. 229, 8918, sec. 2.2,
      ! on the primitive vector; the departures from their form are stated at
      ! Reconstruction_step.positivity_limited_faces).  A face whose
      ! reconstruction is admissible has theta = 1 and is not rewritten.
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)    :: u_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(inout) :: WL, WR
      integer, intent(out) :: n_scaled
      real(fv), dimension(3,1-Ng:N+Ng) :: W_avg
      logical  :: any_bad
      integer  :: j, jr
      real(fv) :: th

      n_scaled = 0
      any_bad = .false.
      do j = 1-Ng,N+Ng
         if (.not. (WL(1,j) .gt. 0.0_fv .and. WL(3,j) .gt. 0.0_fv .and.  &
                    WR(1,j) .gt. 0.0_fv .and. WR(3,j) .gt. 0.0_fv)) then
            any_bad = .true.
            exit
         endif
      enddo
      if (.not. any_bad) return

      call U_to_W_fv(u_in, W_avg)

      do j = 1-Ng,N+Ng
         jr = min(j+1, N+Ng)
         th = positivity_scaling_fv(WL(1,j), WL(3,j), W_avg(1,j), W_avg(3,j))
         if (th .lt. 1.0_fv) then
            WL(:,j) = W_avg(:,j) + th*(WL(:,j) - W_avg(:,j))
            WL(1,j) = max(WL(1,j), epsilon(1.0_fv)*W_avg(1,j))
            WL(3,j) = max(WL(3,j), epsilon(1.0_fv)*W_avg(3,j))
            n_scaled = n_scaled + 1
         endif
         th = positivity_scaling_fv(WR(1,j), WR(3,j), W_avg(1,jr), W_avg(3,jr))
         if (th .lt. 1.0_fv) then
            WR(:,j) = W_avg(:,jr) + th*(WR(:,j) - W_avg(:,jr))
            WR(1,j) = max(WR(1,j), epsilon(1.0_fv)*W_avg(1,jr))
            WR(3,j) = max(WR(3,j), epsilon(1.0_fv)*W_avg(3,jr))
            n_scaled = n_scaled + 1
         endif
      enddo

      end subroutine positivity_limited_faces_fv

      ! ------------------------------------------------------!

      subroutine hydrostatic_equilibrium_of_each_cell_fv(u_in,P_up,P_dn, &
                       p_cell,rho_cell,dev_L,dev_R,dp_eq,pL_asm,pR_asm)
      ! The two face values of the local hydrostatic equilibrium of every
      ! cell: constant density rho_j through the cell's own pressure p_j,
      !     p_eq,j(r) = p_j - rho_j (phi(r) - phi(r_j)),
      ! evaluated at r_edg(j) and at r_edg(j-1) (Kaeppeli and Mishra 2016,
      ! A&A 587, A94, eq. 16).  The cell averages of the pressure and the
      ! density are kept as well: the equilibrium mismatch of a face is
      ! formed from THEM, so that the difference of the two cell pressures
      ! is not first buried in a sum with an O(1) number.  Mirrors
      ! Reconstruction_step.hydrostatic_equilibrium_of_each_cell.
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: P_up, P_dn
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: p_cell, rho_cell
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: dev_L, dev_R, dp_eq
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: pL_asm, pR_asm
      real(fv), dimension(3,1-Ng:N+Ng) :: W
      integer :: j

      call U_to_W_fv(u_in, W)

      p_cell   = W(3,:)
      rho_cell = W(1,:)

      do j = 1-Ng,N+Ng
         P_up(j) = W(3,j)                                                &
                 - W(1,j)*(real(Gphi_i(j),fv) - real(Gphi_c(j),fv))
      enddo
      do j = 2-Ng,N+Ng
         P_dn(j) = W(3,j)                                                &
                 + W(1,j)*(real(Gphi_c(j),fv) - real(Gphi_i(j-1),fv))
      enddo
      ! The lowest ghost cell has no lower face (the cell loop forms no row
      ! for it), and Gphi_i does not reach below r_edg(1-Ng).
      P_dn(1-Ng) = W(3,1-Ng)

      dev_L  = 0.0_fv
      dev_R  = 0.0_fv
      dp_eq  = 0.0_fv
      ! No face state has been assembled yet, so every face is "rewritten"
      ! until the reconstruction claims it.

      end subroutine hydrostatic_equilibrium_of_each_cell_fv

      ! ------------------------------------------------------!

      subroutine well_balanced_face_departures_fv(WL,WR,P_up,P_dn,       &
                       p_cell,rho_cell,pL_asm,pR_asm,dev_L,dev_R,dp_eq)
      ! The equilibrium mismatch of every face,
      !   dp_eq(f) = P_dn(f+1) - P_up(f)
      !            = (p_f+1 - p_f) + rho_f+1 (phi_c(f+1) - phi_i(f))
      !                            + rho_f   (phi_i(f)   - phi_c(f)),
      ! formed in that order so that it is a sum of quantities of the size
      ! of the hydrostatic pressure drop across a cell and vanishes on the
      ! discrete equilibrium; and the departure of any face state the
      ! boundary condition or the positivity limiter rewrote, which has to
      ! be re-formed by subtracting the equilibrium from the O(1) state it
      ! now carries.  The outermost face has no cell to its right and takes
      ! the same cell on both sides, as the flux does.  Mirrors
      ! Reconstruction_step.well_balanced_face_departures.
      real(fv), dimension(3,1-Ng:N+Ng), intent(in) :: WL, WR
      real(fv), dimension(1-Ng:N+Ng), intent(in) :: P_up, P_dn
      real(fv), dimension(1-Ng:N+Ng), intent(in) :: p_cell, rho_cell
      real(fv), dimension(1-Ng:N+Ng), intent(in) :: pL_asm, pR_asm
      real(fv), dimension(1-Ng:N+Ng), intent(inout) :: dev_L, dev_R
      real(fv), dimension(1-Ng:N+Ng), intent(out)   :: dp_eq
      integer :: j, jr

      do j = 1-Ng,N+Ng

         jr = min(j+1,N+Ng)

         if (jr .gt. j) then
            dp_eq(j) = (p_cell(jr) - p_cell(j))                          &
                     + rho_cell(jr)*(real(Gphi_c(jr),fv)                 &
                                   - real(Gphi_i(j),fv))                 &
                     + rho_cell(j) *(real(Gphi_i(j),fv)                  &
                                   - real(Gphi_c(j),fv))
         else
            dp_eq(j) = P_dn(jr) - P_up(j)
         endif

         if (WL(3,j) .ne. pL_asm(j)) dev_L(j) = WL(3,j) - P_up(j)
         if (WR(3,j) .ne. pR_asm(j)) dev_R(j) = WR(3,j) - P_dn(jr)

      enddo

      end subroutine well_balanced_face_departures_fv

      ! ------------------------------------------------------!

      subroutine reconstructed_faces_fv(plm, u_in, WL_out, WR_out,       &
                                        n_scaled, dev_L, dev_R, dp_eq)
      ! Mirrors Reconstruct: under the well-balanced arm the cell's own
      ! equilibrium is formed first, the reconstruction works on the
      ! departure from it, and the face departures and the equilibrium
      ! mismatch are finished after the boundary states and the positivity
      ! limiter.  dev_L, dev_R and dp_eq are set only under the arm.
      logical, intent(in) :: plm
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out
      integer, intent(out) :: n_scaled
      real(fv), dimension(1-Ng:N+Ng), intent(out) :: dev_L, dev_R, dp_eq
      real(fv), dimension(3,1-Ng:N+Ng) :: WL, WR
      real(fv), dimension(1-Ng:N+Ng) :: P_up, P_dn, p_cell, rho_cell
      real(fv), dimension(1-Ng:N+Ng) :: pL_asm, pR_asm
      integer :: j

      dev_L = 0.0_fv
      dev_R = 0.0_fv
      dp_eq = 0.0_fv

      if (well_balanced)                                                 &
         call hydrostatic_equilibrium_of_each_cell_fv(u_in,P_up,P_dn,    &
                 p_cell,rho_cell,dev_L,dev_R,dp_eq,pL_asm,pR_asm)

      if (plm) then
         if (well_balanced) then
            call piecewise_linear_faces_fv(u_in, WL, WR,                 &
                                           P_up, P_dn, dev_L, dev_R)
         else
            call piecewise_linear_faces_fv(u_in, WL, WR)
         endif
      else
         if (well_balanced) then
            call weno3_faces_fv(u_in, WL, WR, P_up, P_dn, dev_L, dev_R)
         else
            call weno3_faces_fv(u_in, WL, WR)
         endif
      endif

      ! The face pressures the reconstruction assembled, before the boundary
      ! states and the positivity limiter may rewrite some of them.
      if (well_balanced) then
         ! Rebuilt from their two parts, as Reconstruct does: bitwise the same where the reconstruction
         ! wrote one, and the equilibrium alone at the outermost faces,
         ! which no reconstruction loop writes.
         do j = 1-Ng,N+Ng
            pL_asm(j) = P_up(j) + dev_L(j)
            pR_asm(j) = P_dn(min(j+1,N+Ng)) + dev_R(j)
         enddo
      endif

      call reconstructed_face_boundary_states_fv(plm, WL, WR, WL_out, WR_out)
      call positivity_limited_faces_fv(u_in, WL_out, WR_out, n_scaled)

      if (well_balanced)                                                 &
         call well_balanced_face_departures_fv(WL_out,WR_out,P_up,P_dn,  &
                 p_cell,rho_cell,pL_asm,pR_asm,dev_L,dev_R,dp_eq)

      end subroutine reconstructed_faces_fv

      ! ------------------------------------------------------!
      ! The interface flux
      ! ------------------------------------------------------!

      subroutine physical_flux_fv(plm, W, PF, jcell)
      logical,  intent(in)  :: plm
      real(fv), intent(in)  :: W(3)
      integer,  intent(in)  :: jcell
      real(fv), intent(out) :: PF(3)
      real(fv) :: rho, v, p, E

      rho = W(1)
      v   = W(2)
      p   = W(3)
      E   = 0.5_fv*rho*v*v + energy_density_from_pressure_fv(jcell,rho,p)

      PF(1) = rho*v
      PF(2) = rho*v*v
      PF(3) = v*(E+p)

      ! Add pressure for PLM discretization.  Not under the well-balanced
      ! arm: there the pressure force of both discretizations is assembled
      ! in the cell loop from the departure of the face pressure from the
      ! cell's own hydrostatic equilibrium, so the momentum flux carried
      ! through the face is the kinetic one alone.
      if (plm .and. .not. well_balanced) PF(2) = PF(2) + p

      end subroutine physical_flux_fv

      ! ------------------------------------------------------!

      subroutine lax_friedrichs_flux_fv(plm, WL, WR, NF, p_out, jL, jR)
      ! Rusanov flux: the average of the two physical fluxes plus a jump term
      ! carrying the spectral radius |v| + a of the flux Jacobian on the
      ! larger side (Rusanov 1961, J. Comput. Math. Phys. USSR 1, 267).
      logical,  intent(in)  :: plm
      real(fv), intent(in)  :: WL(3), WR(3)
      integer,  intent(in)  :: jL, jR
      real(fv), intent(out) :: NF(3), p_out
      real(fv) :: uL(3), uR(3), FL(3), FR(3)
      real(fv) :: rhoL, vL, pL, aL, rhoR, vR, pR, aR, a1

      rhoL = WL(1);  vL = WL(2);  pL = WL(3)
      aL   = sqrt(adiabatic_index_from_state_fv(jL,rhoL,pL)*pL/rhoL)

      rhoR = WR(1);  vR = WR(2);  pR = WR(3)
      aR   = sqrt(adiabatic_index_from_state_fv(jR,rhoR,pR)*pR/rhoR)

      call physical_flux_fv(plm,WL,FL,jL)
      call physical_flux_fv(plm,WR,FR,jR)

      a1 = max(abs(vL)+aL,abs(vR)+aR)

      call W_to_U_comp_fv(WL,uL,jL)
      call W_to_U_comp_fv(WR,uR,jR)

      NF = 0.5_fv*(FL + FR - a1*(uR-uL))

      p_out = 0.5_fv*(pR + pL)

      end subroutine lax_friedrichs_flux_fv

      ! ------------------------------------------------------!

      subroutine hlle_flux_fv(WL, WR, FL, FR, SL, SR, jL, jR, NF)
      ! Harten-Lax-van Leer flux with Einfeldt's wave bounds: positively
      ! conservative whenever SL and SR bound the physical waves (Einfeldt
      ! 1988, SIAM J. Numer. Anal. 25, 294; Einfeldt, Munz, Roe & Sjogreen
      ! 1991, J. Comput. Phys. 92, 273).
      real(fv), intent(in)  :: WL(3), WR(3), FL(3), FR(3), SL, SR
      integer,  intent(in)  :: jL, jR
      real(fv), intent(out) :: NF(3)
      real(fv) :: uL(3), uR(3)

      if (SL .ge. 0.0_fv) then
         NF = FL
      elseif (SR .le. 0.0_fv) then
         NF = FR
      else
         call W_to_U_comp_fv(WL,uL,jL)
         call W_to_U_comp_fv(WR,uR,jR)
         NF = (SR*FL - SL*FR + SL*SR*(uR-uL))/(SR-SL)
      endif

      end subroutine hlle_flux_fv

      ! ------------------------------------------------------!

      subroutine riemann_star_state_fv(WL, WR, gam, u_star, cL_star,     &
                                       cR_star, status)
      ! Star state of the face Riemann problem, Toro (2009) chapter 9: the
      ! primitive-variable estimate (eq. 9.20), the two-rarefaction estimate
      ! (eq. 9.32) and the two-shock estimate (eq. 9.42), selected by the
      ! pressure ratio, with the vacuum condition of eq. 4.76 tested first.
      ! Valid for one ideal gas with a constant adiabatic index on both
      ! sides; gam is the mean of the two face values.
      real(fv), intent(in)  :: WL(3), WR(3), gam
      real(fv), intent(out) :: u_star, cL_star, cR_star
      integer,  intent(out) :: status
      real(fv) :: rhoL,uL,pL,cL, rhoR,uR,pR,cR
      real(fv) :: rho_bar,c_bar, p_min,p_max,Q
      real(fv) :: p_pv,u_pv,rhoL_pv,rhoR_pv, z,pLR
      real(fv) :: AL,BL,AR,BR,gL,gR,p_guess
      real(fv) :: p_star, rho_L_star, rho_R_star
      logical  :: pvrs_admissible
      real(fv), parameter :: Q_user = 2.0_fv

      rhoL = WL(1);  uL = WL(2);  pL = WL(3)
      rhoR = WR(1);  uR = WR(2);  pR = WR(3)

      p_star     = 0.0_fv
      rho_L_star = 0.0_fv
      rho_R_star = 0.0_fv
      u_star     = 0.0_fv
      cL_star    = 0.0_fv
      cR_star    = 0.0_fv
      status     = ROE_STAR_INADMISSIBLE

      if (.not. (rhoL .gt. 0.0_fv .and. rhoR .gt. 0.0_fv .and.           &
                 pL   .gt. 0.0_fv .and. pR   .gt. 0.0_fv)) return

      cL = sqrt(gam*pL/rhoL)
      cR = sqrt(gam*pR/rhoR)

      if (uR - uL .ge. 2.0_fv*(cL+cR)/(gam-1.0_fv)) then
         status = ROE_STAR_VACUUM
         return
      endif

      p_min = min(pL,pR)
      p_max = max(pL,pR)
      Q     = p_max/p_min
      rho_bar = 0.5_fv*(rhoL+rhoR)
      c_bar   = 0.5_fv*(cL+cR)

      p_pv    = 0.5_fv*(pL+pR) - 0.5_fv*(uR-uL)*rho_bar*c_bar
      u_pv    = 0.5_fv*(uL+uR) - 0.5_fv*(pR-pL)/(rho_bar*c_bar)
      rhoL_pv = rhoL + (p_pv-pL)/(cL*cL)
      rhoR_pv = rhoR + (p_pv-pR)/(cR*cR)

      pvrs_admissible = (Q .le. Q_user)      .and.                       &
                        (p_pv .ge. p_min)    .and.                       &
                        (p_pv .le. p_max)    .and.                       &
                        (p_pv .gt. 0.0_fv)   .and.                       &
                        (rhoL_pv .gt. 0.0_fv) .and.                      &
                        (rhoR_pv .gt. 0.0_fv)

      if (pvrs_admissible) then

         p_star     = p_pv
         u_star     = u_pv
         rho_L_star = rhoL_pv
         rho_R_star = rhoR_pv

      elseif (p_pv .lt. p_min) then

         z   = 0.5_fv*(gam-1.0_fv)/gam
         pLR = (pL/pR)**z

         p_star = ((cL + cR - 0.5_fv*(gam-1.0_fv)*(uR-uL))/              &
                   (cL/pL**z + cR/pR**z))**(1.0_fv/z)
         u_star = (pLR*uL/cL + uR/cR                                     &
                   + 2.0_fv*(pLR-1.0_fv)/(gam-1.0_fv))/                  &
                  (pLR/cL + 1.0_fv/cR)
         rho_L_star = rhoL*(p_star/pL)**(1.0_fv/gam)
         rho_R_star = rhoR*(p_star/pR)**(1.0_fv/gam)

      else

         AL = 2.0_fv/((gam+1.0_fv)*rhoL)
         BL = (gam-1.0_fv)/(gam+1.0_fv)*pL
         AR = 2.0_fv/((gam+1.0_fv)*rhoR)
         BR = (gam-1.0_fv)/(gam+1.0_fv)*pR

         p_guess = max(0.0_fv, p_pv)
         gL = sqrt(AL/(p_guess + BL))
         gR = sqrt(AR/(p_guess + BR))

         p_star = (gL*pL + gR*pR - (uR-uL))/(gL+gR)
         u_star = 0.5_fv*(uL+uR)                                         &
                + 0.5_fv*((p_star-pR)*gR - (p_star-pL)*gL)
         rho_L_star = rhoL*(p_star/pL + (gam-1.0_fv)/(gam+1.0_fv))       &
                      /((gam-1.0_fv)/(gam+1.0_fv)*p_star/pL             &
                        + 1.0_fv)
         rho_R_star = rhoR*(p_star/pR + (gam-1.0_fv)/(gam+1.0_fv))       &
                      /((gam-1.0_fv)/(gam+1.0_fv)*p_star/pR             &
                        + 1.0_fv)

      endif

      if (.not. (p_star     .gt. 0.0_fv .and.                            &
                 rho_L_star .gt. 0.0_fv .and.                            &
                 rho_R_star .gt. 0.0_fv)) then
         u_star = 0.0_fv
         status = ROE_STAR_INADMISSIBLE
         return
      endif

      cL_star = sqrt(gam*p_star/rho_L_star)
      cR_star = sqrt(gam*p_star/rho_R_star)
      status  = ROE_STAR_OK

      end subroutine riemann_star_state_fv

      ! ------------------------------------------------------!

      subroutine interface_flux_fv(plm, WL, WR, NF, p_out, jL, jR,       &
                                   n_hlle, n_llf,                        &
                                   dev_L, dev_R, dp_eq, q_up, q_dn)
      ! Numerical flux of the face between cells jL and jR.  Each side's
      ! total energy and sound speed carry ITS OWN heat capacity; the single
      ! gamma the Roe average and the star-state estimate are derived for is
      ! the arithmetic mean of the two, which is that value exactly wherever
      ! the two sides agree.
      !
      ! WELL-BALANCED ARM.  With the optional arguments present the PRESSURE
      ! JUMP of the Riemann problem is not pR - pL, the difference of two
      ! O(1) face pressures, but
      !     dp = dp_eq + dev_R - dev_L
      ! from the mismatch of the two neighboring hydrostatic equilibria at
      ! this face and the two face states' departures from their own
      ! equilibrium, all three of the size of the departure.  On a discrete
      ! hydrostatic equilibrium the three vanish together and the face is
      ! the stationary contact discontinuity the well-balanced property
      ! requires the flux to resolve exactly (Kaeppeli and Mishra 2016, A&A
      ! 587, A94, their eq. 10).  The two optional outputs are the face
      ! pressure measured from each side's own equilibrium,
      ! q_up = p_out - P_up(jL) and q_dn = p_out - P_dn(jR).
      ! Mirrors Num_flux.
      logical,  intent(in)  :: plm
      real(fv), intent(in)  :: WL(3), WR(3)
      integer,  intent(in)  :: jL, jR
      real(fv), intent(out) :: NF(3), p_out
      integer, intent(inout) :: n_hlle, n_llf
      real(fv), intent(in),  optional :: dev_L, dev_R, dp_eq
      real(fv), intent(out), optional :: q_up, q_dn
      logical  :: wb, hllc_left
      real(fv) :: dp_wb
      real(fv) :: gam_L, gam_R, gam_face
      real(fv) :: uL(3), uR(3), FL(3), FR(3)
      real(fv) :: rhoL,vL,pL,aL,EL,HL,SL,phiL
      real(fv) :: rhoR,vR,pR,aR,ER,HR,SR,phiR
      real(fv) :: S_star, usL(3), usR(3)
      real(fv) :: s, rho_avg, v_avg, H_avg, a_avg
      real(fv) :: l1, l2, l3, l1L, l1R, l3L, l3R
      real(fv) :: v_star, aL_star, aR_star
      real(fv) :: drho, dvel, dp, a1, a2, a3
      real(fv), dimension(3) :: K1, K2, K3
      integer :: roe_status

      wb = present(dev_L) .and. present(dev_R) .and. present(dp_eq)
      hllc_left = .true.
      dp_wb = 0.0_fv
      if (wb) dp_wb = dp_eq + dev_R - dev_L

      rhoL = WL(1)
      vL   = WL(2)
      pL   = WL(3)
      EL   = 0.5_fv*rhoL*vL*vL + energy_density_from_pressure_fv(jL,rhoL,pL)
      HL   = (EL+pL)/rhoL
      gam_L = adiabatic_index_from_state_fv(jL,rhoL,pL)
      aL   = sqrt(gam_L*pL/rhoL)

      rhoR = WR(1)
      vR   = WR(2)
      pR   = WR(3)
      ER   = 0.5_fv*rhoR*vR*vR + energy_density_from_pressure_fv(jR,rhoR,pR)
      HR   = (ER+pR)/rhoR
      gam_R = adiabatic_index_from_state_fv(jR,rhoR,pR)
      aR   = sqrt(gam_R*pR/rhoR)
      gam_face = 0.5_fv*(gam_L + gam_R)

      select case(flux)

      case ('LLF')

         call lax_friedrichs_flux_fv(plm,WL,WR,NF,p_out,jL,jR)

         ! p_out is the average of the two face pressures
         if (wb) then
            q_up = 0.5_fv*(dev_L + dev_R + dp_eq)
            q_dn = q_up - dp_eq
         endif

         n_llf = n_llf + 1

      case('HLLC')

         call W_to_U_comp_fv(WL,uL,jL)
         call W_to_U_comp_fv(WR,uR,jR)

         SL = min(0.0_fv,min(vL-aL, vR-aR))
         SR = max(0.0_fv,max(vL+aL, vR+aR))
         phiL = rhoL*(SL-vL)
         phiR = rhoR*(SR-vR)
         if (wb) then
            S_star = (dp_wb + phiL*vL - phiR*vR)/(phiL - phiR)
         else
            S_star = (pR - pL + phiL*vL - phiR*vR)/(phiL - phiR)
         endif

         usL(1) = 1.0_fv
         usL(2) = S_star
         usL(3) = EL/rhoL + (S_star-vL)*(S_star + pL/phiL)
         usL = rhoL*(SL-vL)/(SL-S_star)*usL

         usR(1) = 1.0_fv
         usR(2) = S_star
         usR(3) = ER/rhoR + (S_star-vR)*(S_star + pR/phiR)
         usR = rhoR*(SR-vR)/(SR-S_star)*usR

         if(SL.ge.(0.0_fv)) then
            call physical_flux_fv(plm,WL,NF,jL)
            p_out = pL
            hllc_left = .true.
         elseif(SL.lt.(0.0_fv).and.S_star.ge.(0.0_fv)) then
            call physical_flux_fv(plm,WL,NF,jL)
            NF = NF + SL * (usL-uL)
            p_out = pL
            hllc_left = .true.
         elseif(S_star.lt.(0.0_fv).and.SR.ge.(0.0_fv)) then
            call physical_flux_fv(plm,WR,NF,jR)
            NF = NF + SR * (usR-uR)
            p_out = pR
            hllc_left = .false.
         else
            call physical_flux_fv(plm,WR,NF,jR)
            p_out = pR
            hllc_left = .false.
         endif

         ! p_out is one side's own face pressure here, so the departure it
         ! carries is that side's.  hllc_left says which side the branch
         ! above took.
         if (wb) then
            if (hllc_left) then
               q_up = dev_L
               q_dn = dev_L - dp_eq
            else
               q_dn = dev_R
               q_up = dev_R + dp_eq
            endif
         endif

      case ('ROE')

         s = sqrt(rhoL)/(sqrt(rhoL)+sqrt(rhoR))

         rho_avg = sqrt(rhoL*rhoR)

         v_avg = s*vL + (1.0_fv-s)*vR

         H_avg = s*HL + (1.0_fv-s)*HR

         a_avg = sqrt((gam_face-1.0_fv)*(H_avg-0.5_fv*v_avg*v_avg))

         l1 = v_avg - a_avg
         l2 = v_avg
         l3 = v_avg + a_avg

         call riemann_star_state_fv(WL,WR,gam_face,v_star,aL_star,       &
                                    aR_star,roe_status)

         if (roe_status .ne. ROE_STAR_OK) then

            ! No admissible star state: the Roe linearization is not
            ! positivity preserving there, so the face takes the HLLE flux,
            ! whose bounds pair each side's own signal speed with the Roe
            ! average and so enclose the waves of the linearized problem.
            SL = min(vL - aL, v_avg - a_avg)
            SR = max(vR + aR, v_avg + a_avg)

            call physical_flux_fv(plm,WL,FL,jL)
            call physical_flux_fv(plm,WR,FR,jR)
            call hlle_flux_fv(WL,WR,FL,FR,SL,SR,jL,jR,NF)

            n_hlle = n_hlle + 1

         else

            l1L = vL - aL
            l1R = v_star - aL_star
            l3L = v_star + aR_star
            l3R = vR + aR

            if (l1L.lt.(0.0_fv).and.l1R.gt.(0.0_fv)) then
               l1 = l1L*(l1R-l1)/(l1R-l1L)
            endif

            if (l3L.lt.(0.0_fv).and.l3R.gt.(0.0_fv)) then
               l3 = l3R*(l3-l3L)/(l3R-l3L)
            endif

            K1(1) = 1.0_fv
            K1(2) = v_avg - a_avg
            K1(3) = H_avg - v_avg*a_avg

            K2(1) = 1.0_fv
            K2(2) = v_avg
            K2(3) = 0.5_fv*v_avg*v_avg

            K3(1) = 1.0_fv
            K3(2) = v_avg + a_avg
            K3(3) = H_avg + v_avg*a_avg

            ! Under the well-balanced arm the pressure jump is the one
            ! built from the departures (see the header), not the
            ! difference of the two O(1) face pressures.
            drho = rhoR - rhoL
            dvel = vR - vL
            dp   = pR - pL
            if (wb) dp = dp_wb

            a1 = 0.5_fv/a_avg**2.0_fv*(dp - rho_avg*a_avg*dvel)
            a2 = drho - dp/a_avg**2.0_fv
            a3 = 0.5_fv/a_avg**2.0_fv*(dp + rho_avg*a_avg*dvel)

            call physical_flux_fv(plm,WL,FL,jL)
            call physical_flux_fv(plm,WR,FR,jR)

            NF = 0.5_fv*(FR+FL)                                          &
               - 0.5_fv*(a1*abs(l1)*K1 + a2*abs(l2)*K2 + a3*abs(l3)*K3)

         endif

         p_out = 0.5_fv*(pR + pL)

         ! p_out is the average of the two face pressures, on the Roe branch
         ! and on the HLLE fallback alike
         if (wb) then
            q_up = 0.5_fv*(dev_L + dev_R + dp_eq)
            q_dn = q_up - dp_eq
         endif

      case default

         write(*,*) 'ERROR: unknown numerical flux: ', trim(flux)
         write(*,*) '  allowed: LLF, HLLC, ROE'
         error stop 1

      end select

      end subroutine interface_flux_fv

      ! ------------------------------------------------------!
      ! The rows
      ! ------------------------------------------------------!

      subroutine geometric_and_gravitational_source_fv(plm,j,dr,dAp,dAm, &
                                                       dV,u,WR,WL,S)
      logical,  intent(in)  :: plm
      integer,  intent(in)  :: j
      real(fv), intent(in)  :: dr,dAp,dAm,dV
      real(fv), intent(in)  :: u(3),WR(3),WL(3)
      real(fv), intent(out) :: S(3)
      real(fv) :: rhoL,rhoR,pC

      rhoL = WR(1)
      rhoR = WL(1)

      if (.not. well_balanced)                                           &
         pC = pressure_from_energy_density_fv(j, u(1),                   &
                 u(3)-0.5_fv*u(2)*u(2)/u(1))

      ! WELL-BALANCED ARM.  The gravitational source and the geometric
      ! pressure term are then not evaluated here at all: they cancel the
      ! equilibrium part of the momentum flux difference IDENTICALLY, and
      ! the cell loop assembles the momentum row from what is left, the
      ! departure of the face pressure from the cell's own hydrostatic
      ! equilibrium (Kaeppeli and Mishra 2014, J. Comput. Phys. 259, 199,
      ! eq. 2.26).  Returning zero here is not the absence of gravity: it
      ! is gravity carried in the flux difference in a form that vanishes
      ! on the discrete equilibrium.
      if (well_balanced) then
         S = 0.0_fv
         return
      endif

      S(1) = 0.0_fv
      S(2) = - 0.5_fv*(rhoL + rhoR)                                      &
             *(real(Gphi_i(j),fv) - real(Gphi_i(j-1),fv))/dr
      S(3) = 0.0_fv

      if (plm) S(2) = S(2) + (dAp-dAm)/dV*pC

      end subroutine geometric_and_gravitational_source_fv

      ! ------------------------------------------------------!

      subroutine flux_difference_and_source_of_scheme_fv(plm,u_in,WL,WR, &
                       dF,S,ff,fpr,n_scaled,n_hlle,n_llf)
      ! One scheme's flux difference and source over the padded range: the
      ! face loop assembles every interface flux, the cell loop differences
      ! the two faces of each cell over its volume and adds the source.  The
      ! lowest ghost cell has no lower face, so it carries no equation and
      ! both outputs are zero there.
      logical, intent(in) :: plm
      real(fv), dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real(fv), dimension(3,1-Ng:N+Ng), intent(out) :: WL,WR,dF,S,ff
      real(fv), dimension(1-Ng:N+Ng),   intent(out) :: fpr
      integer, intent(out) :: n_scaled
      integer, intent(inout) :: n_hlle, n_llf
      integer  :: j, nh_j, nl_j
      real(fv) :: dr, rp, rm, dAp, dAm, dV, dF3p, pL, pR
      real(fv), dimension(3) :: Fp, Fm
      ! WELL-BALANCED ARM.  The departures and the equilibrium mismatch the
      ! reconstruction returned, and the face pressure of every interface
      ! measured from the hydrostatic equilibrium of the cell on each side,
      !   face_q_up(f) = face_p(f) - P_up(f),
      !   face_q_dn(f) = face_p(f) - P_dn(f+1),
      ! both of the size of the departure from equilibrium.  They are what
      ! the momentum row is assembled from when the arm is on: the O(1)
      ! pressure of the cell and the gravitational source cancel in the
      ! ALGEBRA instead of in floating point.
      real(fv), dimension(1-Ng:N+Ng) :: wb_dev_L, wb_dev_R, wb_dp_eq
      real(fv), dimension(1-Ng:N+Ng) :: face_q_up, face_q_dn
      real(fv) :: qp, qm

      call reconstructed_faces_fv(plm, u_in, WL, WR, n_scaled,           &
                                  wb_dev_L, wb_dev_R, wb_dp_eq)

      face_q_up = 0.0_fv
      face_q_dn = 0.0_fv

      dF(:,1-Ng) = 0.0_fv
      S(:,1-Ng)  = 0.0_fv

      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,Fp,pR,nh_j,nl_j,qp,qm) reduction(+:n_hlle,n_llf)
      do j = 1-Ng,N+Ng
         nh_j = 0;  nl_j = 0
         if (well_balanced) then
            call interface_flux_fv(plm,WL(:,j),WR(:,j),Fp,pR,j,          &
                                   min(j+1,N+Ng),nh_j,nl_j,              &
                                   wb_dev_L(j),wb_dev_R(j),wb_dp_eq(j),  &
                                   qp,qm)
            face_q_up(j) = qp
            face_q_dn(j) = qm
         else
            call interface_flux_fv(plm,WL(:,j),WR(:,j),Fp,pR,j,          &
                                   min(j+1,N+Ng),nh_j,nl_j)
         endif
         ff(:,j) = Fp
         fpr(j)   = pR
         n_hlle  = n_hlle + nh_j
         n_llf   = n_llf  + nl_j
      enddo
      !$omp end parallel do

      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,dr,rp,rm,dAp,dAm,dV,Fp,Fm,pL,pR,dF3p)
      do j = 2-Ng,N+Ng

         dr = real(dr_j(j), fv)
         rp = real(r_edg(j), fv)
         rm = real(r_edg(j-1), fv)
         dAp = rp*rp
         dAm = rm*rm
         dV = (dAp*rp - dAm*rm)/3.0_fv

         Fm = ff(:,j-1)
         pL = fpr(j-1)
         Fp = ff(:,j)
         pR = fpr(j)

         call geometric_and_gravitational_source_fv(plm,j,dr,dAp,dAm,dV, &
                     u_in(:,j),WR(:,j-1),WL(:,j),S(:,j))

         dF(1,j) = (dAp*Fp(1) - dAm*Fm(1))/dV
         dF(2,j) = (dAp*Fp(2) - dAm*Fm(2))/dV

         if (.not. plm)  dF(2,j) = dF(2,j) + (pR - pL)/dr

         ! THE MOMENTUM ROW UNDER THE WELL-BALANCED ARM.  The pressure is
         ! not in Fp(2)/Fm(2) here (physical_flux_fv leaves it out), and
         ! what replaces it is the face pressure measured from THIS cell's
         ! own hydrostatic equilibrium.  The two identities that make the
         ! substitution exact, with P_up = p_j - rho_j (phi_i(j) - phi_c(j))
         ! and P_dn = p_j + rho_j (phi_c(j) - phi_i(j-1)):
         !
         !   A+ P_up - A- P_dn - (A+ - A-) p_j
         !        = -rho_j [ A+ (phi_i(j) - phi_c(j))
         !                 + A- (phi_c(j) - phi_i(j-1)) ]        (PLM form)
         !   P_up - P_dn = -rho_j (phi_i(j) - phi_i(j-1))        (WENO3 form)
         !
         ! The left-hand sides are the equilibrium part of the pressure
         ! terms the row carries (the flux difference and, under PLM, the
         ! geometric source); the right-hand sides are the discrete
         ! gravitational source of the well-balanced scheme.  They cancel
         ! here in the algebra, cell by cell, with the cell's own O(1)
         ! pressure never appearing, so the source returns S(2,j) = 0,
         ! NEITHER side is evaluated, and what is left is the departure
         ! alone.  This is the momentum source of Kaeppeli and Mishra (2014,
         ! J. Comput. Phys. 259, 199, their eq. 2.26) in spherical geometry.
         if (well_balanced) then
            if (plm) then
               dF(2,j) = (dAp*(Fp(2) + face_q_up(j))                     &
                        - dAm*(Fm(2) + face_q_dn(j-1)))/dV
            else
               dF(2,j) = (dAp*Fp(2) - dAm*Fm(2))/dV                      &
                       + (face_q_up(j) - face_q_dn(j-1))/dr
            endif
         endif

         dF3p  = dAp*Fp(1)*(real(Gphi_i(j),fv) - real(Gphi_c(j),fv))     &
               - dAm*Fm(1)*(real(Gphi_i(j-1),fv) - real(Gphi_c(j),fv))
         dF(3,j) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV

      enddo
      !$omp end parallel do

      end subroutine flux_difference_and_source_of_scheme_fv

      ! ------------------------------------------------------!

      subroutine hydrodynamic_flux_difference_and_source(u8,WL8,WR8,dF8, &
                       S8,ff8,fp8,n_scaled,n_hlle,n_llf,               &
                       dF_at_working_kind)
      ! THE ENTRY POINT.  Takes the conserved state in double, converts it
      ! exactly to fv, evaluates the reconstruction, the interface flux, the
      ! sources and the flux difference of the scheme (or of the PLM/WENO3
      ! homotopy at the current lambda) entirely in fv, and rounds the result
      ! to double ONCE, on the way out.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u8
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL8,WR8,dF8,S8,ff8
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: fp8
      integer, intent(out) :: n_scaled, n_hlle, n_llf
      ! The rows BEFORE the single rounding to double, so that a test can
      ! measure the arithmetic of this instantiation below the width of a
      ! double.  Nothing on a run's path asks for them.
      real(fv), dimension(3,1-Ng:N+Ng), intent(out), optional ::         &
                                                    dF_at_working_kind
      real(fv), dimension(3,1-Ng:N+Ng) :: u,WL,WR,dF,S,ff
      real(fv), dimension(3,1-Ng:N+Ng) :: WLw,WRw,dFw,Sw,ffw
      real(fv), dimension(1-Ng:N+Ng)   :: fpr,fprw
      real(fv) :: lam, om
      integer  :: ns_w

      u = real(u8, fv)
      n_scaled = 0;  n_hlle = 0;  n_llf = 0

      if (.not. recon_lambda_on) then
         call flux_difference_and_source_of_scheme_fv(use_plm,u,WL,WR,   &
                 dF,S,ff,fpr,n_scaled,n_hlle,n_llf)
      else
         lam = real(recon_lambda, fv)
         om  = 1.0_fv - lam
         if (lam .le. 0.0_fv) then
            call flux_difference_and_source_of_scheme_fv(.true.,u,WL,WR, &
                    dF,S,ff,fpr,n_scaled,n_hlle,n_llf)
         else if (lam .ge. 1.0_fv) then
            call flux_difference_and_source_of_scheme_fv(.false.,u,WL,WR,&
                    dF,S,ff,fpr,n_scaled,n_hlle,n_llf)
         else
            call flux_difference_and_source_of_scheme_fv(.true.,u,WL,WR, &
                    dF,S,ff,fpr,n_scaled,n_hlle,n_llf)
            call flux_difference_and_source_of_scheme_fv(.false.,u,WLw,  &
                    WRw,dFw,Sw,ffw,fprw,ns_w,n_hlle,n_llf)
            n_scaled = n_scaled + ns_w
            dF = om*dF + lam*dFw
            S  = om*S  + lam*Sw
            WL = om*WL + lam*WLw
            WR = om*WR + lam*WRw
            ff = om*ff + lam*ffw
            fpr = om*fpr + lam*fprw
         endif
      endif

      if (present(dF_at_working_kind)) dF_at_working_kind = dF

      WL8 = real(WL, kind(1.0d0))
      WR8 = real(WR, kind(1.0d0))
      dF8 = real(dF, kind(1.0d0))
      S8  = real(S,  kind(1.0d0))
      ff8 = real(ff, kind(1.0d0))
      fp8 = real(fpr, kind(1.0d0))

      end subroutine hydrodynamic_flux_difference_and_source
