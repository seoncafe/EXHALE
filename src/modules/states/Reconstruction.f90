      module Reconstruction_step
      ! Collection of reconstruction procedure subroutine
      
      use global_parameters
      use BC_Apply
      use PLM_reconstruction
      use Conversion

      implicit none

      ! Stored ESWENO3 smoothness factors for the frozen-weights mode
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

      contains

      subroutine Reconstruct(u_in,WL_out,WR_out) 
      
      integer :: j,k
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng) :: WL,WR
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out

      ! ESWENO3 variables
      real*8, dimension(3,1-Ng:N+Ng) :: W,dW
      real*8, dimension(3) :: dWp,dWm
      real*8, dimension(3) :: b0,b1
      real*8, dimension(3) :: tau
      real*8, dimension(3) :: S0,S1
      real*8 :: rm,rp
      real*8, dimension(1-Ng:N+Ng) :: dV,C1,C2,D1,D2
	
      select case (rec_method)
      
         !-----------------------------------------------!
         
         case ('PLM')      ! Piecewise linear reconstruction

            call PLM_rec(u_in,WL,WR)

         case ('WENO3')    ! ESWENO3 Reconstruction
			
            ! Convert to primitive variables
            call U_to_W(u_in,W)
         
            ! Evaluate jumps at interfaces
            dW(:,1-Ng:N+Ng-1) = W(:,2-Ng:N+Ng) - W(:,1-Ng:N+Ng-1)
            dW(:,N+Ng) = 0.0
            
            ! Calculate cell volumes
            do j = 0, N+1			
               rp = r_edg(j)
               rm = r_edg(j-1)
               dV(j) = (rp*rp*rp - rm*rm*rm)
            enddo
            dV(1-Ng) = dV(2-Ng)
            dV(N+Ng) = dV(N+Ng-1)
            
            
            ! Evaluate reconstruction geometry-dependent coefficients
            do j = 1-Ng,N+Ng-1
               C1(j) = dV(j+1)/(dV(j) + dV(j+1))
               C2(j) = 1.0 - C1(j)
            enddo
            
            do j = 2-Ng,N+Ng-1  	
               D1(j) = dV(j+1)/(dV(j) + dV(j-1))
               D2(j) = dV(j-1)/(dV(j) + dV(j+1))
            enddo
            
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
            !$omp   private(j,k,dWp,dWm,b0,b1,tau,S0,S1)
            do j = 0,N+1

               dWp = dW(:,j)
               dWm = dW(:,j-1)

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

               WL(:,j) = W(:,j)	&
                  + (S0*C1(j)*dWp + D1(j)*S1*C1(j-1)*dWm) &
                  /(S0 + D1(j)*S1)
               
               WR(:,j-1) = W(:,j)  &
                  - (D2(j)*S0*C2(j)*dWp + S1*C2(j-1)*dWm) &
                  /(D2(j)*S0 + S1)
            enddo
            !$omp end parallel do

         case default

            write(*,*) 'ERROR: unknown reconstruction scheme: ', trim(rec_method)
            write(*,*) '  allowed: PLM, WENO3'
            error stop 1

      end select

      ! Apply BC to reconstructed variables
      call Rec_BC(WL,WR,WL_out,WR_out)

      ! Last step before the Riemann solver sees the face states: restore
      ! rho > 0 and p > 0 wherever the reconstruction lost them.
      call positivity_limited_faces(u_in,WL_out,WR_out)

      ! End of subroutine
      end subroutine Reconstruct

      !-----------------------------------------------!

      subroutine positivity_limited_faces(u_in,WL,WR)

      ! Scale each face state CONTINUOUSLY toward its own cell average, by the
      ! largest factor that keeps rho > 0 and p > 0 there.
      !
      !     W_face  <-  W_avg + theta ( W_rec - W_avg ),   theta in [0,1]
      !     theta   =   min over rho and p of  (q_avg - eps)/(q_avg - q_rec)
      !
      ! This is the positivity-preserving limiter of Zhang & Shu (2010,
      ! J. Comput. Phys., 229, 8918; doi:10.1016/j.jcp.2010.08.016), in the
      ! form that scales the reconstruction increment of one cell toward that
      ! cell's admissible mean. The publisher PDF is not in references/, so
      ! only the bibliography is given here.
      !
      ! WHY IT IS CONTINUOUS AND WHY THAT IS THE POINT (docs/Update_EXHALE.md
      ! section 138; the measurement is P49). This routine used to be a hard
      ! switch: as soon as a reconstructed rho or p crossed zero, the WHOLE
      ! face pair was replaced by the two cell averages. That is a STEP
      ! DISCONTINUITY of the residual F(Y), and it is what stopped the
      ! molecular arm's steady solve. Measured on the hot Uranus hand-off:
      ! the WENO3 left density at the face of cell 273 (r = 1.2324 R_p) sat at
      ! 2.50e-15 against cell averages of 1e-4 -- eleven orders below and
      ! positive by a hair -- so the iterate sat exactly on the switching
      ! surface. Raising the neighbouring cell's density by one part in 1e9
      ! tipped it through zero, the face density jumped to 4.17e-6, the mass
      ! residual of cell 273 jumped by a factor 97, and the merit jumped from
      ! 202.5 to 945.1 -- BY THE SAME FACTOR 4.67 at every step size from 1e-2
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
      ! at the crossing is 2.2e-16 of the cell average -- round-off, against
      ! the factor 1.7e9 the hard switch made.
      !
      ! BYTE-IDENTITY. A face whose reconstruction is admissible has theta = 1
      ! and is NOT rewritten (W_avg + 1*(W - W_avg) is not bitwise W), so every
      ! run in which the guard never fired is unchanged to the bit. Runs in
      ! which it did fire move, and by more than the floor: the old switch
      ! replaced both states of a face and this scales only the offending one.
      !
      ! Validity of the reconstruction, and why the limiter is needed. PLM and
      ! ESWENO3 both extrapolate the PRIMITIVE variables (rho, v, p) to a face
      ! with slopes taken from the neighbouring cells, which is a valid
      ! approximation only while the solution varies smoothly across the
      ! stencil. It has no positivity property of its own: a face value can
      ! cross zero while every cell average on the stencil is positive. The
      ! state that does it here is a hypersonic layer, where the pressure the
      ! scheme recovers as p = (gamma-1)(E - rho v^2/2) is the difference of two
      ! nearly equal numbers -- at Mach 60 the thermal pressure is 2e-4 of the
      ! total energy density -- so a reconstruction across the neighbouring jump
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
      ! has no cell j+1, so it falls back on the last cell average -- the
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
      ! NaN face state -- which compares false against everything -- is caught
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
            n_faces_positivity_limited = n_faces_positivity_limited + 1
         endif
         th = positivity_scaling(WR(1,j), WR(3,j), W_avg(1,jr), W_avg(3,jr))
         if (th .lt. 1.0d0) then
            WR(:,j) = W_avg(:,jr) + th*(WR(:,j) - W_avg(:,jr))
            n_faces_positivity_limited = n_faces_positivity_limited + 1
         endif
      enddo

      ! End of subroutine
      end subroutine positivity_limited_faces


      !-----------------------------------------------!

      double precision function positivity_scaling(rho_f, p_f, rho_a, p_a)  &
                                                                 result(th)
      ! The largest theta in [0,1] for which rho and p of
      ! W_avg + theta (W_face - W_avg) both stay at or above their floors.
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
      ! q_a + theta (q_f - q_a) exactly on the floor eps = epsilon*q_a when the
      ! reconstruction went below it, and 1 when it did not. q_a > 0 is the
      ! caller's precondition.
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
