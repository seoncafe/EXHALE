      module grid_construction
      ! Constructs radial grid, edges and sized of grids

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
      real*8 :: tol = 1.0
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
         
         do while(tol.ge.(1.0e-8))
               
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
      ! Cell size (N+2*Ng points) - dr(j) = dimension of cell j
      dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)
      dr_j(1-Ng) = dr_j(2-Ng)
      dr_j(N+Ng) = dr_j(N+Ng-1)
      
      
      ! Do a smoothing of the mixed-type grid
      ! The 1-2-1 pass runs downward from the uniform/stretched junction to
      ! the base, so its start index is N_low (it was written as the literal
      ! 50 when N_low itself was hardcoded to 50).
      if (grid_type .eq. 'Mixed') then

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
         ! Cell size (N+2*Ng points) - dr(j) = dimension of cell j
         dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)
         dr_j(1-Ng) = dr_j(2-Ng)
         dr_j(N+Ng) = dr_j(N+Ng-1)
   
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

      ! End of subroutine
      end subroutine define_grid
      
      ! End of module
      end module grid_construction
