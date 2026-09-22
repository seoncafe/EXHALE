      module grid_construction
      ! Constructs the radial grid, its faces and its cell sizes, and holds
      ! the one spherical-shell geometry every transport operator on this
      ! grid divides by (spherical_face_area_and_cell_volume below).

      use global_parameters
      
      implicit none
      
      contains

      subroutine define_grid
      integer :: j
      ! Base resolution of the Mixed grid: N_low uniform cells of size drc
      ! [R_p] on the base, then N_up stretched cells out to r_max. Both come
      ! from input.inp ("Base grid [dr,cells]:") and default to the historical
      ! hardcoded values 2.0e-4 and 50; see the scale-height requirement
      ! documented with dr_base in parameters.f90.
      integer :: N_low
      integer :: N_up
      real*8 :: drc
      real*8 :: x0,x1
      real*8 :: f,df
      real*8 :: tol   ! Newton convergence measure of the stretch factor;
                      ! reset before each search (an initializer here would
                      ! imply SAVE and skip the loop on every later call)
      real*8 :: q
      real*8 :: dr

      N_low = N_low_cells
      drc   = dr_base
      N_up  = N - N_low

      if (grid_type .eq. 'Mixed') then
         if (N_low .lt. 2 .or. N_low .gt. N-10) then
            write(*,'(A,I0,A,I0,A)') ' (define_grid.f90) ERROR: "Base grid'// &
               ' cells: ', N_low, '" must lie in [2,', N-10, '] (the '//      &
               'stretched region needs the remaining cells).'
            error stop 1
         endif
         if (drc .le. 0.0d0 .or. 1.0d0 + N_low*drc .ge. r_max) then
            write(*,'(A,ES10.3,A,I0,A,F7.3,A)') ' (define_grid.f90) ERROR: '//&
               'base grid spacing ', drc, ' R_p x ', N_low, ' cells does '//  &
               'not fit inside the domain r_max = ', r_max, ' R_p.'
            error stop 1
         endif
      endif

      select case (grid_type)
      
      case ('Uniform')
         !------ Uniform spaced grid ------!
         
         ! Grid spacing
         dr = (r_max-1.0)/(1.0*N)
         
         ! Lower ghost cells
         r(1-Ng) = 1.0
         
         ! Loop for others cell centers
         do j = 2-Ng,N+Ng
               r(j) = r(j-1) + dr
         enddo
      
      !--------------------------------------------------
       
      case ('Stretched')
        
         !------ Regular stretched grid ------!
         r   = (/ (r_max**((j-1+Ng)*1.0/(N*1.0 + 2.0*Ng - 1.0) ), j = 1-Ng,N+Ng) /)

      !--------------------------------------------------

      case ('Mixed') 
      
         !------ Mixed grid ------!
         ! Constructed with N_low uniform spaced points
         ! and N_up points in a stretched grid
         
         ! Lower ghost cells
         r(1-Ng) = 1.0 - drc
         r(2-Ng) = 1.0
         
         ! Grid centers in the uniform region
         do j = 1,N_low
            r(j) = r(j-1) + drc
         enddo
         
         ! Solve for the region of stretched grid
         ! Solves the equation f(x) = 0 using 
         ! Newton-Raphson method with
         ! 
         !     f(x) = (1-x^N)/(1-x) - (r_up-r_low)/dr
         ! 
         ! x = stretch parameter
         ! r_up, r_low = upper and lower boundary 
         ! of the domain
         ! dr = starting grid dimension
         
         ! Initial guess
         x0 = 1.01
         
         tol = 1.0d0
         do while(tol.ge.(1.0d-8))
               
            ! Evaluate function and its derivative
            
            ! Auxiliary constant
            q  = (r_max-r(N_low))/drc
            
            ! Function
            f  = (1.0-x0**(1.0*N_up))/(1.0-x0) - q
            
            ! Analytic function derivative
            df = (f  + q - N_up*x0**(N_up-1.0))/(1.0-x0)
               
            ! Guess of solution
            x1 = x0 - f/df

            ! Evaluate tolerance
            tol = abs(x1-x0)
            
            ! Update point for the next step
            x0 = x1
               
         ! End of while loop
         enddo
         
         ! Construct stretched grid
         do j = N_low + 1,N
            r(j) = r(j-1) + x0**(1.0*j - N_low -1)*drc
         enddo
         
         ! Add ghost points at the top of the domain
         do j = 1,Ng
            r(N+j) = 2.0*r(N+j-1) - r(N+j-2)
         enddo
       
       !--------------------------------------------------

      case default

         write(*,*) 'ERROR: unknown grid type: ', trim(grid_type)
         write(*,*) '  allowed: Uniform, Stretched, Mixed'
         error stop 1

      end select
      
      !--- Cell edges r_{j+1/2} ---!
      
      ! Cell edges (N+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
      r_edg(1-Ng:N+Ng-1) = 0.5*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      
      !--- Cell dimensions r_{j+1/2} - r_{j-1/2} --- !
      ! The width stored for cell j is the distance between the two faces of
      ! THAT cell,  dr_j(j) = r_edg(j) - r_edg(j-1),  with r_edg(j) = r_{j+1/2}.
      ! This is an exact identity of the finite-volume discretization, not an
      ! approximation: it is the same face pair RK_rhs takes for cell j to form
      ! the volume dV = (r_edg(j)^3 - r_edg(j-1)^3)/3, the same pair Source
      ! divides (Gphi_i(j) - Gphi_i(j-1)) by to get the gravity of cell j, and
      ! the same width the column integrals n(j)*dr_j(j) and the cell optical
      ! depths use with the cell-centred densities. Tested to round-off.
      ! r_edg(-Ng) lies outside the array, so the innermost ghost takes the
      ! width of its neighbour.
      dr_j(2-Ng:N+Ng) = r_edg(2-Ng:N+Ng) - r_edg(1-Ng:N+Ng-1)
      dr_j(1-Ng) = dr_j(2-Ng)
      
      
      ! Do a smoothing of the mixed-type grid
      ! The 1-2-1 pass runs downward from the uniform/stretched junction to
      ! the base, so its start index is N_low (it was written as the literal
      ! 50 when N_low itself was hardcoded to 50).
      if (grid_type .eq. 'Mixed') then

         ! The pass below smooths the widths and then rebuilds the centres
         ! from them as r(j) = r(j-1) + (dr_j(j) + dr_j(j-1))/2, i.e. half of
         ! cell j-1 plus half of cell j: it reads dr_j(j) as the width of
         ! cell j, which is what the definition above stores.
         do j = N_low,2-Ng,-1
            dr_j(j) = 0.25*(dr_j(j-1) + 2.0*dr_j(j) + dr_j(j+1))
         enddo
	      
         r(2-Ng) = 1.0
         r(1-Ng) = r(2-Ng) - 0.5*(dr_j(1-Ng) + dr_j(2-Ng))
   
         do j = 3-Ng,N+Ng
            r(j) = r(j-1) + 0.5*(dr_j(j) + dr_j(j-1))
         enddo  
      
      	! Rescale to [1,r_max]
         r = (r-1)/(r(N+Ng) - 1.0)*(r_max - 1.0) + 1.0
      	
      	! Re-eval edges and cell size
      	
      	! Cell edges (N+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
         r_edg(1-Ng:N+Ng-1) = 0.5*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
         r_edg(N+Ng) = 2.0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
		
         !--- Cell dimensions r_{j+1/2} - r_{j-1/2} --- !
         ! Same identity as above: dr_j(j) = r_edg(j) - r_edg(j-1) is the
         ! width of cell j itself, the innermost ghost taking its neighbour's.
         dr_j(2-Ng:N+Ng) = r_edg(2-Ng:N+Ng) - r_edg(1-Ng:N+Ng-1)
         dr_j(1-Ng) = dr_j(2-Ng)
   
      endif

      !-----------------------------!
      
      ! Select relevant domain for constant momentum
      do j = 1-Ng,N+Ng
         if(r(j).ge.r_esc) goto 111
      enddo

111   j_min = j

      ! Guard: if the escape radius lies outside the (L1-truncated) domain,
      ! NO cell satisfies r>=r_esc and the DO loop above leaves j = N+Ng+1,
      ! so j_min > N and the convergence range [j_min:N] is EMPTY. maxval /
      ! minval over an empty array return -/+HUGE, giving du = Inf and
      ! dtu = -HUGE < dtu_th -> a spurious "steady state" exit at step 0.
      ! This bites deep-RLOF cases whose Roche lobe (r_max) is smaller than
      ! the default escape radius (e.g. WASP-121b Case D: r_max ~ 1.35 R_p <
      ! r_esc = 1.5 R_p). Clamp the convergence range to mid-domain and warn
      ! loudly; the user should set a smaller "Escape radius [R_p]:".
      if (j_min .gt. N) then
         do j = 1-Ng,N+Ng
            if (r(j) .ge. (1.0d0 + 0.5d0*(r_max-1.0d0))) goto 112
         enddo
112      j_min = j
         write(*,'(A,F6.3,A,F6.3,A)')                                       &
            '    (define_grid.f90) WARNING: escape radius r_esc = ', r_esc,  &
            ' R_p >= domain r_max = ', r_max, ' R_p.'
         write(*,*) '       The escape radius is outside the L1-truncated ' //&
                    'domain; the constant-momentum'
         write(*,*) '       convergence range would be empty. Clamping it ' //&
                    'to mid-domain (r >= '
         write(*,'(A,F6.3,A)') '        ', 1.0d0+0.5d0*(r_max-1.0d0),        &
                    ' R_p). Set a smaller "Escape radius [R_p]:" in input.inp.'
      endif
      ! The window is a range of PHYSICAL cells, so its first index is cell 1
      ! even when the escape radius sits at or below the base: the ghosts
      ! (1-Ng .. 0) carry the boundary closure, not a solution, and a spread
      ! taken over them is not a property of the wind. Same clamp as j_flux
      ! below.
      j_min = max(j_min, 1)


      ! Inner edge of the FLUX window: first cell with r >= r_flux, over which
      ! the flux gate measures the spread of rho*v*r^2 (section 133). Same
      ! guard as j_min: if r_flux lies outside the domain the window would be
      ! empty, so clamp it to mid-domain and warn.
      do j = 1-Ng,N+Ng
         if (r(j) .ge. r_flux) goto 121
      enddo
121   j_flux = j
      if (j_flux .gt. N) then
         do j = 1-Ng,N+Ng
            if (r(j) .ge. (1.0d0 + 0.5d0*(r_max-1.0d0))) goto 122
         enddo
122      j_flux = j
         write(*,'(A,F6.3,A,F6.3,A)')                                       &
            '    (define_grid.f90) WARNING: flux-window radius r_flux = ',   &
            r_flux, ' R_p >= domain r_max = ', r_max,                        &
            ' R_p; clamping the flux gate window to mid-domain.'
      endif
      j_flux = max(j_flux, 1)

      ! The faces of every physical cell are strictly ordered, which is what
      ! makes the shell between them a body with a positive volume.  A grid
      ! that does not meet it carries no finite-volume divergence at all, so
      ! it is refused here instead of being floored downstream.
      do j = 1, N
         if (r_edg(j) .le. r_edg(j-1)) then
            write(*,'(A,I0,A,ES13.6,A,ES13.6,A)')                        &
               ' (define_grid.f90) ERROR: cell ', j, ' has faces ',       &
               r_edg(j-1), ' and ', r_edg(j),                             &
               ' R_p, so it has no positive volume.'
            error stop 1
         endif
      enddo

      ! End of subroutine
      end subroutine define_grid

      ! ---------------------------------------------------------------- !

      subroutine spherical_face_area_and_cell_volume(face_area,          &
                                                     cell_volume)
      ! THE ONE SPHERICAL GEOMETRY EVERY TRANSPORT OPERATOR ON THIS GRID
      ! DIVIDES BY.  For the shell between the faces r_- = r_edg(j-1) and
      ! r_+ = r_edg(j),
      !
      !    A(f) = r_edg(f)^2
      !    V(j) = (r_+^3 - r_-^3)/3 = (r_+ - r_-)(r_+^2 + r_+ r_- + r_-^2)/3
      !
      ! so that the conservative divergence
      !
      !    D(F)_j = [ A_+ F_+ - A_- F_- ] / V_j
      !
      ! telescopes exactly down the column: an internal face leaves cell j
      ! carrying the same A F it enters cell j+1 with, and the sum over the
      ! column is the difference of the two boundary face fluxes alone.
      ! Every contribution to one conserved quantity therefore has to divide
      ! by THIS V(j).  Replacing it in one term by r_j^2 (r_+ - r_-) leaves
      ! that term weighted by V_j/(r_j^2 dr_j), a factor that differs
      ! between neighbours on a stretched grid, and the internal faces stop
      ! cancelling: the mixed operator is then the divergence of no single
      ! flux.
      !
      ! The factored second form of V is the one evaluated.  The difference
      ! of cubes cancels its leading digits when the cell is thin against
      ! its radius, and the base cells of the production grids here are
      ! dr/r ~ 2e-4, where it loses about three decimal digits that the
      ! factored form keeps.
      !
      ! UNITS.  Both arrays are dimensionless, in the radius unit of r and
      ! r_edg (the planetary radius R_p).  A module working in cm multiplies
      ! the area by R0^2 and the volume by R0^3.
      !
      ! It is filled from the CURRENT r_edg at every call rather than stored,
      ! so that a caller that built a grid of its own by writing r_edg -- the
      ! refinement rows of the acceptance suites do -- gets the geometry of
      ! the grid it built and not of a previous one.
      !
      ! It covers the PHYSICAL column, faces 0..N and cells 1..N, which is
      ! the range the element and carrier transport operators assemble their
      ! rows on.  A row assembled over the padded range as well -- the
      ! hydrodynamic rows and the species flux divergence are -- reads the
      ! volume of one cell at a time from spherical_cell_volume below, which
      ! is the expression this fills the array with.
      real*8, dimension(0:N), intent(out) :: face_area
      real*8, dimension(1:N), intent(out) :: cell_volume
      integer :: j

      do j = 0, N
         face_area(j) = r_edg(j)*r_edg(j)
      enddo
      do j = 1, N
         cell_volume(j) = spherical_cell_volume(j)
         if (cell_volume(j) .le. 0.0d0) then
            write(*,'(A,I0,A,ES13.6,A,ES13.6,A)')                        &
               ' (define_grid.f90) ERROR: cell ', j, ' has faces ',       &
               r_edg(j-1), ' and ', r_edg(j),                             &
               ' R_p, so its volume is not positive.'
            error stop 1
         endif
      enddo

      end subroutine spherical_face_area_and_cell_volume

      ! ---------------------------------------------------------------- !

      real*8 function spherical_cell_volume(j)
      ! THE ONE EXPRESSION FOR THE SHELL VOLUME OF ONE CELL, the value the
      ! array above is filled with and the value every operator that
      ! transports anything through this grid divides its row by:
      !
      !    V(j) = (r_+^3 - r_-^3)/3 = (r_+ - r_-)(r_+^2 + r_+ r_- + r_-^2)/3
      !
      ! between the faces r_- = r_edg(j-1) and r_+ = r_edg(j).  The factored
      ! second form is the one evaluated: the difference of cubes cancels its
      ! leading digits when the cell is thin against its radius, and the base
      ! cells of the production grids here are dr/r ~ 2e-4, where it loses
      ! about three decimal digits that the factored form keeps.  Two
      ! spellings of this volume leave two operators on the same column
      ! disagreeing at that size about what a cell holds, which is why there
      ! is one.
      !
      ! It takes the cell index rather than an array so that a row assembled
      ! over the padded index range 2-Ng..N+Ng can read it: the ghost cells
      ! lie outside the physical column the array form covers.  It reads the
      ! CURRENT r_edg, for the same reason the array form is filled and not
      ! stored.
      !
      ! UNITS.  Dimensionless, in the radius unit of r and r_edg (the
      ! planetary radius R_p).  A module working in cm multiplies by R0^3.
      !
      ! The positivity of the volume is asserted by the array form over the
      ! physical column and by define_grid over the faces it builds; it is
      ! not restated per call, which would put a branch in the innermost
      ! loop of every flux difference.
      integer, intent(in) :: j
      real*8 :: rm, rpl

      rm  = r_edg(j-1)
      rpl = r_edg(j)
      spherical_cell_volume = (rpl - rm)*(rpl*rpl + rpl*rm + rm*rm)/3.0d0

      end function spherical_cell_volume

      ! ---------------------------------------------------------------- !

      ! Index of the cell whose center lies nearest a given radius [R_p].
      ! minloc counts positions from 1 whatever the declared lower bound of
      ! the array is, while the grid arrays are declared 1-Ng:N+Ng, so the
      ! position it returns becomes a subscript of r only through the lbound
      ! offset. Written once here so that every caller asking "which cell is
      ! at radius x" gets the same answer.
      integer function cell_nearest_radius(r_target) result(j_near)
      real*8, intent(in) :: r_target
      j_near = minloc(abs(r - r_target), dim = 1) + lbound(r,1) - 1
      end function cell_nearest_radius
      
      ! End of module
      end module grid_construction
