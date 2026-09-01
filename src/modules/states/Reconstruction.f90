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

      ! Number of faces at which the reconstruction left the admissible
      ! thermodynamic state (rho > 0, p > 0) and was dropped to first order,
      ! summed over the whole run. Reported at the end of a run; zero means the
      ! high-order reconstruction was admissible everywhere and the run is the
      ! one an unguarded build would have produced.
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

      ! Drop a face pair to the piecewise-constant (first-order) limit of the
      ! same reconstruction wherever the higher-order face state has left the
      ! admissible thermodynamic state, rho > 0 and p > 0.
      !
      ! Validity of the reconstruction, and why the test is needed. PLM and
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
      ! so the repair is to use them: WL(face j) = W(cell j), WR(face j) =
      ! W(cell j+1). That is the same scheme at first order, so it is a local
      ! loss of accuracy at an under-resolved face, not a change of the
      ! equations. It is the standard positivity guard of a high-resolution
      ! finite-volume scheme.
      !
      ! Outermost face: there is no cell j+1, so the right state falls back on
      ! the last cell average -- the zero-gradient state the outer BC
      ! reconstructs to in any case.
      !
      ! If a CELL AVERAGE is itself non-positive the state is unphysical before
      ! any reconstruction and nothing here can repair it; that face is left as
      ! it is, and the NaN detector of the marching loop reports it.

      real*8, dimension(3,1-Ng:N+Ng), intent(in)    :: u_in
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: WL,WR
      real*8, dimension(3,1-Ng:N+Ng) :: W_avg
      logical :: any_bad
      integer :: j,jr

      ! Scan first: the conversion to cell-average primitives is only needed
      ! when a face has actually failed, which on an admissible run is never.
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
         if (WL(1,j) .gt. 0.0d0 .and. WL(3,j) .gt. 0.0d0 .and.          &
             WR(1,j) .gt. 0.0d0 .and. WR(3,j) .gt. 0.0d0) cycle
         jr = min(j+1, N+Ng)
         if (.not. (W_avg(1,j)  .gt. 0.0d0 .and. W_avg(3,j)  .gt. 0.0d0 &
              .and. W_avg(1,jr) .gt. 0.0d0 .and. W_avg(3,jr) .gt. 0.0d0)) cycle
         WL(:,j) = W_avg(:,j)
         WR(:,j) = W_avg(:,jr)
         n_faces_positivity_limited = n_faces_positivity_limited + 1
      enddo

      ! End of subroutine
      end subroutine positivity_limited_faces


      ! End of module
      end module Reconstruction_step
