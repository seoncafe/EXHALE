      module grid_construction
      ! Constructs the radial grid, its faces and its cell sizes, and holds
      ! the one spherical-shell geometry every transport operator on this
      ! grid divides by (spherical_face_area_and_cell_volume below).

      use global_parameters
      
      implicit none
      
      ! The two stretch ratios of the grid, stated by define_grid for the
      ! setup report: the width ratio dr_j(j+1)/dr_j(j) of the stretched
      ! region of the Mixed grid (0 for the other grid types), and the width
      ! ratio of the outer shells (0 without them). The second is the one the
      ! cell count of "Outer shells" sets, and the first is the ratio it
      ! continues.
      real*8 :: mixed_stretch_ratio      = 0.0d0
      real*8 :: outer_shells_width_ratio = 0.0d0

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
      ! Cells of the CONSTRUCTED grid, the one "Grid type", "Grid cells",
      ! "Base grid" and "Outer radius" define. The outer shells, when
      ! "Outer shells" asks for them, are appended beyond its outer face as
      ! cells nc+1 .. N, so every statement of the construction below runs
      ! on the index range 1-Ng .. nc+Ng and produces the same centers,
      ! faces and widths whether or not shells follow.
      integer :: nc

      nc    = N - n_outer_shells
      N_low = N_low_cells
      drc   = dr_base
      N_up  = nc - N_low

      if (grid_type .eq. 'Mixed') then
         if (N_low .lt. 2 .or. N_low .gt. nc-10) then
            write(*,'(A,I0,A,I0,A)') ' (define_grid.f90) ERROR: "Base grid'// &
               ' cells: ', N_low, '" must lie in [2,', nc-10, '] (the '//     &
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

      mixed_stretch_ratio = 0.0d0

      select case (grid_type)
      
      case ('Uniform')
         !------ Uniform spaced grid ------!
         
         ! Grid spacing
         dr = (r_max-1.0)/(1.0*nc)
         
         ! Lower ghost cells
         r(1-Ng) = 1.0
         
         ! Loop for others cell centers
         do j = 2-Ng,nc+Ng
               r(j) = r(j-1) + dr
         enddo
      
      !--------------------------------------------------
       
      case ('Stretched')
        
         !------ Regular stretched grid ------!
         r(1-Ng:nc+Ng) = (/ (r_max**((j-1+Ng)*1.0/(nc*1.0 + 2.0*Ng - 1.0) ), &
                             j = 1-Ng,nc+Ng) /)

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
         
         mixed_stretch_ratio = x0

         ! Construct stretched grid
         do j = N_low + 1,nc
            r(j) = r(j-1) + x0**(1.0*j - N_low -1)*drc
         enddo
         
         ! Add ghost points at the top of the domain
         do j = 1,Ng
            r(nc+j) = 2.0*r(nc+j-1) - r(nc+j-2)
         enddo
       
       !--------------------------------------------------

      case default

         write(*,*) 'ERROR: unknown grid type: ', trim(grid_type)
         write(*,*) '  allowed: Uniform, Stretched, Mixed'
         error stop 1

      end select
      
      !--- Cell edges r_{j+1/2} ---!
      
      ! Cell edges (nc+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
      r_edg(1-Ng:nc+Ng-1) = 0.5*(r(1-Ng:nc+Ng-1) + r(2-Ng:nc+Ng))
      r_edg(nc+Ng) = 2.0*r_edg(nc+Ng-1) - r_edg(nc+Ng-2)
      
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
      dr_j(2-Ng:nc+Ng) = r_edg(2-Ng:nc+Ng) - r_edg(1-Ng:nc+Ng-1)
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
   
         do j = 3-Ng,nc+Ng
            r(j) = r(j-1) + 0.5*(dr_j(j) + dr_j(j-1))
         enddo  
      
      	! Rescale to [1,r_max]: the center of the outermost ghost cell of
         ! the constructed grid is placed at r_max, so its outer PHYSICAL
         ! face r_edg(nc) lies about one and a half cells inside r_max.
         r(1-Ng:nc+Ng) = (r(1-Ng:nc+Ng)-1)/(r(nc+Ng) - 1.0)*(r_max - 1.0) &
                         + 1.0
      	
      	! Re-eval edges and cell size
      	
      	! Cell edges (nc+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
         r_edg(1-Ng:nc+Ng-1) = 0.5*(r(1-Ng:nc+Ng-1) + r(2-Ng:nc+Ng))
         r_edg(nc+Ng) = 2.0*r_edg(nc+Ng-1) - r_edg(nc+Ng-2)
		
         !--- Cell dimensions r_{j+1/2} - r_{j-1/2} --- !
         ! Same identity as above: dr_j(j) = r_edg(j) - r_edg(j-1) is the
         ! width of cell j itself, the innermost ghost taking its neighbour's.
         dr_j(2-Ng:nc+Ng) = r_edg(2-Ng:nc+Ng) - r_edg(1-Ng:nc+Ng-1)
         dr_j(1-Ng) = dr_j(2-Ng)
   
      endif

      ! Shells beyond the outer face of the constructed grid ("Outer
      ! shells"). They replace its two outer ghost cells by physical cells
      ! and leave cells 1-Ng .. nc and faces 1-Ng .. nc as constructed.
      outer_shells_width_ratio = 0.0d0
      if (n_outer_shells .gt. 0) call append_outer_shells(nc)

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
            if (r(j) .ge. (1.0d0 + 0.5d0*(domain_outer_radius()-1.0d0)))   &
               goto 112
         enddo
112      j_min = j
         write(*,'(A,F6.3,A,F6.3,A)')                                       &
            '    (define_grid.f90) WARNING: escape radius r_esc = ', r_esc,  &
            ' R_p >= domain r_max = ', domain_outer_radius(), ' R_p.'
         write(*,*) '       The escape radius is outside the L1-truncated ' //&
                    'domain; the constant-momentum'
         write(*,*) '       convergence range would be empty. Clamping it ' //&
                    'to mid-domain (r >= '
         write(*,'(A,F6.3,A)') '        ',                                   &
                    1.0d0+0.5d0*(domain_outer_radius()-1.0d0),               &
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
            if (r(j) .ge. (1.0d0 + 0.5d0*(domain_outer_radius()-1.0d0)))   &
               goto 122
         enddo
122      j_flux = j
         write(*,'(A,F6.3,A,F6.3,A)')                                       &
            '    (define_grid.f90) WARNING: flux-window radius r_flux = ',   &
            r_flux, ' R_p >= domain r_max = ', domain_outer_radius(),        &
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

      subroutine append_outer_shells(nc)
      ! SHELLS BEYOND THE OUTER FACE OF THE CONSTRUCTED GRID.
      !
      ! The constructed grid of nc cells ends at the face r_edg(nc). This
      ! appends m = n_outer_shells cells beyond it, cells nc+1 .. N = nc+m,
      ! whose last face r_edg(N) is the radius r_F = r_outer_shells_face of
      ! "Outer shells [r_face,cells]:", and leaves everything at and below
      ! r_edg(nc) as the construction made it: the centers r(1-Ng .. nc),
      ! the faces r_edg(1-Ng .. nc), the widths dr_j(1-Ng .. nc) and with
      ! them every volume below r_edg(nc) are the ones a run without shells
      ! builds, bit for bit.
      !
      ! Every face of this grid is the midpoint of the two centers beside it,
      ! r_edg(j) = (r(j) + r(j+1))/2, which is what the reconstruction reads
      ! (a face value is extrapolated from the center by (r(j+1) - r(j))/2).
      ! The shells keep that rule, so they are built from their centers:
      !
      !   r(nc+1)            = the center of the first outer ghost of the
      !                        constructed grid, unchanged, which keeps
      !                        r_edg(nc) = (r(nc) + r(nc+1))/2 as it was;
      !   r(nc+k) - r(nc+k-1) = s x^(k-1),  k = 2 .. m,  s = r(nc+1) - r(nc);
      !   r(N+1)             = 2 r_F - r(N), so that r_edg(N) = r_F;
      !   r(N+j)             = 2 r(N+j-1) - r(N+j-2), j = 2 .. Ng, the
      !                        linear extrapolation the Mixed grid applies
      !                        to its own outer ghosts.
      !
      ! With center spacings in geometric progression the widths are too,
      ! dr_j(nc+k) = s (1 + x) x^(k-1)/2, so x is the width ratio of the
      ! shells, and it is the root of
      !
      !   S(x) = sum_{i=0}^{m-1} x^i + x^m/2 = (r_F - r(nc))/s,
      !
      ! the statement r_edg(N) = r_F with the ghost spacing r(N+1) - r(N)
      ! continuing the progression, s x^m. S is increasing and convex for
      ! x > 0, so Newton's method started above the root descends to it
      ! monotonically; S(0) = 1, so a positive root exists iff r_F lies
      ! beyond r(nc+1). The cell count m decides x: the width ratio is
      ! continuous across the old outer face when m is chosen so that x is
      ! the stretch ratio of the grid below (mixed_stretch_ratio for the
      ! Mixed grid). Both ratios are written to the setup report.
      integer, intent(in) :: nc
      integer :: j, k, m, it
      real*8  :: s, q, x, xn, sx, dsx, sp, rf

      m  = n_outer_shells
      rf = r_outer_shells_face
      s  = r(nc+1) - r(nc)
      if (.not. (rf .gt. r(nc+1))) then
         write(*,'(A,ES13.6,A,ES13.6,A)') ' (define_grid.f90) ERROR: '//     &
            '"Outer shells" face r = ', rf, ' R_p does not lie beyond'//     &
            ' the first outer ghost center of the grid, ', r(nc+1), ' R_p.'
         error stop 1
      endif
      q = (rf - r(nc))/s

      ! A starting point above the root, then Newton down to it. The
      ! descent stops when a step no longer lowers x: that is the rounding
      ! floor of S(x) - q.
      x = 1.0d0
      call shell_spacing_sum(x, m, sx, dsx)
      do while (sx .lt. q)
         x = 2.0d0*x
         call shell_spacing_sum(x, m, sx, dsx)
      enddo
      do it = 1, 200
         call shell_spacing_sum(x, m, sx, dsx)
         xn = x - (sx - q)/dsx
         if (.not. (xn .lt. x)) exit
         x = xn
      enddo
      outer_shells_width_ratio = x

      sp = s
      do k = 2, m
         sp = sp*x
         r(nc+k) = r(nc+k-1) + sp
      enddo
      ! The ghost center that puts the last face on r_F. One rounding in the
      ! difference and one in the sum can leave the midpoint one unit in the
      ! last place off r_F; the neighbouring double is then taken, so the
      ! face lands on r_F exactly.
      r(N+1) = 2.0d0*rf - r(N)
      do it = 1, 4
         if (0.5d0*(r(N) + r(N+1)) .eq. rf) exit
         r(N+1) = nearest(r(N+1), rf - 0.5d0*(r(N) + r(N+1)))
      enddo
      do j = 2, Ng
         r(N+j) = 2.0d0*r(N+j-1) - r(N+j-2)
      enddo

      ! Faces and widths above r_edg(nc), by the rules of define_grid.
      r_edg(nc+1:N+Ng-1) = 0.5d0*(r(nc+1:N+Ng-1) + r(nc+2:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j(nc+1:N+Ng) = r_edg(nc+1:N+Ng) - r_edg(nc:N+Ng-1)

      if (r_edg(N) .ne. rf) then
         write(*,'(A,ES23.16,A,ES23.16,A)') ' (define_grid.f90) ERROR: '// &
            'the outer face of the shells is ', r_edg(N), ' R_p, not the ', &
            rf, ' R_p "Outer shells" states.'
         error stop 1
      endif

      end subroutine append_outer_shells

      ! ---------------------------------------------------------------- !

      subroutine shell_spacing_sum(x, m, sx, dsx)
      ! S(x) = sum_{i=0}^{m-1} x^i + x^m/2 and its derivative, by Horner's
      ! rule, which has no division by x - 1 and so no special case at the
      ! uniform spacing x = 1.
      real*8,  intent(in)  :: x
      integer, intent(in)  :: m
      real*8,  intent(out) :: sx, dsx
      integer :: i
      sx  = 0.5d0
      dsx = 0.0d0
      do i = 1, m
         dsx = dsx*x + sx
         sx  = sx*x + 1.0d0
      enddo
      end subroutine shell_spacing_sum

      ! ---------------------------------------------------------------- !

      real*8 function domain_outer_radius()
      ! THE OUTER RADIUS OF THE DOMAIN THE RUN SOLVES [R_p]: the outer face
      ! of the last shell when "Outer shells" appends shells, else r_max,
      ! the radius the grid is constructed for ("Outer radius" or the
      ! Roche/Hill radius), which is the center of the outermost ghost cell
      ! of the Mixed and Stretched grids.
      if (n_outer_shells .gt. 0) then
         domain_outer_radius = r_outer_shells_face
      else
         domain_outer_radius = r_max
      endif
      end function domain_outer_radius

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
