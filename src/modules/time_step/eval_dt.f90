   module eval_time_step

   use global_parameters
   use grid_construction, only: spherical_cell_volume
   use caloric_eos, only: caloric_mixture_active,                     &
                         adiabatic_index_from_state

   implicit none

   contains

   subroutine eval_dt(W,dt,dt_loc)
   ! THE EXPLICIT-STABLE INTERVAL OF THE STATE, cell by cell and as a
   ! minimum over the cells that are evolved.
   !
   ! dt     = global step (the minimum over the evolved cells 1..N).
   ! dt_loc = the same restriction held cell by cell, the pseudo-time the RK
   !          stages and the source updates advance one cell with.  With
   !          "Time stepping: Local" each cell takes its own interval (the
   !          fixed point dF = S, heating = cooling, is dt-independent);
   !          otherwise dt_loc is uniformly the global dt.
   !
   ! THE RESTRICTION, DERIVED FOR THE SCHEME IN USE.
   !
   ! The hydrodynamic row of cell j is assembled in RK_rhs as the spherical
   ! finite-volume divergence
   !
   !     dF_j = ( A_{j+1/2} F_{j+1/2} - A_{j-1/2} F_{j-1/2} ) / V_j ,
   !     A = r^2 ,   V_j = ( r_{j+1/2}^3 - r_{j-1/2}^3 )/3 ,
   !
   ! and integrated by SSP-RK3, every stage of which is a convex combination
   ! of forward-Euler steps of the same length dt (SSP coefficient 1), so
   ! the restriction is the forward-Euler one of that divergence.
   !
   ! F at a face is a Riemann flux of the two reconstructed states bounding
   ! the face (HLLC, Roe or Rusanov, with PLM or WENO3 reconstruction).
   ! Every one of them is built from the waves of that face's Riemann
   ! problem, whose fastest signal is bounded by
   !
   !     S_{j+1/2} = max( |v| + c )  over the two states bounding the face.
   !
   ! Godunov's condition is that the fan opened at a face has not left the
   ! cell it enters when the step ends.  Written with the geometry the
   ! divergence above actually uses, the face carries signal into cell j at
   ! the rate A_f S_f per unit volume V_j, so
   !
   !     dt_j = CFL * V_j / max( A_{j-1/2} S_{j-1/2} , A_{j+1/2} S_{j+1/2} )
   !
   ! with the Godunov number CFL <= 1.  In the Cartesian uniform limit
   ! A_{j-1/2} = A_{j+1/2} = A and V_j = A dr this is CFL*dr/(|v| + c), the
   ! convention the "CFL:" key carries (default 0.6), so the input number
   ! keeps its meaning.  The stricter condition that the two fans of a cell
   ! may not meet each other inside it replaces the max by the sum of the
   ! two face terms, which is the same restriction at CFL/2.
   !
   ! In spherical geometry V_j/A_{j+1/2} < dr_j < V_j/A_{j-1/2}: the outer
   ! face of a cell is the binding one, and the cell-centred width dr_j is
   ! not an upper bound of the interval on either face.  On the base grids
   ! used here dr/r is 2e-4 and the two differ by that much; on the
   ! stretched outer grid dr/r reaches order unity and they do not.
   !
   ! WHICH CELLS ARE RESTRICTED, AND WHERE THE GHOSTS ENTER.  Only cells
   ! 1..N are evolved: the ghost rows are overwritten by Apply_BC after
   ! every stage, so a ghost's own width is not a stability requirement and
   ! the lowest ghost does not even carry a width of its own (define_grid
   ! copies dr_j(1-Ng) from the cell above it).  The ghosts nevertheless
   ! bound the two boundary faces r_edg(0) and r_edg(N), which are faces of
   ! the evolved cells 1 and N, so the base Riemann problem restricts cell 1
   ! with the reservoir's own signal speed although no ghost cell is
   ! evolved, and the outer one restricts cell N.
   !
   ! VALIDITY.  S_f is formed from the CELL AVERAGES bounding the face, not
   ! from the reconstructed face states, which can carry a larger |v| + c
   ! inside a steep cell; that margin is the one the Godunov number below 1
   ! leaves, and it is the margin this expression has always carried.
   real*8, dimension(3,1-Ng:N+Ng), intent(in) :: W
   real*8, dimension(1-Ng:N+Ng) :: rho,v,p,cs
   real*8, dimension(0:N) :: s_face
   real*8 :: dt_j
   integer :: j
   real*8, intent(out) :: dt
   real*8, dimension(1-Ng:N+Ng), intent(out) :: dt_loc

   ! Extract physical variables
   rho = W(1,:)
   v   = W(2,:)
   p   = W(3,:)

   ! Evaluate sound speed
   if (caloric_mixture_active) then
      do j = 1-Ng, N+Ng
         cs(j) = sqrt(adiabatic_index_from_state(j,rho(j),p(j))         &
                      *p(j)/rho(j))
      enddo
   else
   cs = sqrt(gamma_ad*p/rho)
   endif

   ! The fastest signal of each face of an evolved cell, the two boundary
   ! faces included.
   do j = 0, N
      s_face(j) = max(abs(v(j))   + cs(j),                              &
                      abs(v(j+1)) + cs(j+1))
   enddo

   ! The interval of each evolved cell, and the minimum over them.
   dt = huge(1.0d0)
   do j = 1, N
      dt_j = CFL*spherical_cell_volume(j)                               &
             /max(r_edg(j)*r_edg(j)*s_face(j),                          &
                  r_edg(j-1)*r_edg(j-1)*s_face(j-1))
      dt_loc(j) = dt_j
      dt = min(dt, dt_j)
   enddo

   if (use_local_dt) then
      ! A ghost row is not evolved: it is rebuilt by Apply_BC after every
      ! stage. It is given the interval of the evolved cell it borders so
      ! that a stage formed over the padded range is arithmetic on a
      ! physical interval and not on an undefined one.
      dt_loc(1-Ng:0)   = dt_loc(1)
      dt_loc(N+1:N+Ng) = dt_loc(N)
   else
      dt_loc = dt
   endif

   ! End of subroutine
   end subroutine eval_dt

   ! End of module
   end module eval_time_step
