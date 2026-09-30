      module Initialization
      ! Initialize the simulation:
      ! - Set global options (thread number, time loop counter, initial
      !                       initial momentum variation)
      ! - Construct spatial grid
      ! - Construct energy grids for ionization equilibrium
      ! - Construct initial conditions (load or set)
      
      use global_parameters
      use grid_construction
      use energy_vectors_construct
      use gravity_grid_construction
      use initial_conditions
      use Conversion
      use BC_Apply
      use omp_lib
      use blas_thread_policy, only: blas_threads_set_policy, blas_threads_report
      use IC_load
      use initial_conditions
      use composition, only: get_species_densities
      use species_table, only: n_mion
      use wae_exhale_bridge, only: wae_generate_ic
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use base_boundary, only: base_reservoir_p
      
      implicit none 
      
      contains
      
      subroutine init(W,u,f_sp)
      ! Initialize the simulation setup

      integer :: n_omp_threads, omp_env_st
      integer :: n_omp_obtained
      character(len=32) :: omp_env
      real*8, dimension(1-Ng:N+Ng)   :: rho,v,p,T
      ! composition scratch: only n_part_cell1 (set inside get_species_densities)
      ! is wanted here, but the single policy point returns the whole set
      real*8, dimension(1-Ng:N+Ng)   :: nhi,nhii,nhei,nheii,nheiii,nheiTR
      real*8, dimension(1-Ng:N+Ng)   :: ne,n_tot
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng,n_species),intent(out) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),intent(out) :: W,u
      logical :: ex_outdir
      integer :: rc_outdir
      
      write(*,*) '(init.f90) Initializing the simulation..'

      ! Every output of the run -- from the IC dump a few calls below to the
      ! final profiles -- goes to ./output, and nothing else creates it, so a
      ! fresh run directory used to die at the first write with a bare
      ! Fortran runtime error that named neither the directory nor the cure.
      ! Create it here; if that fails (permissions, a FILE named output),
      ! stop and say what is missing.
      !
      ! The test is "can a file be created inside it", not an INQUIRE on
      ! the directory: the standard does not define INQUIRE for a
      ! directory, gfortran answers .true. for 'output/.' and ifx answers
      ! .false. (measured 2026-09-05: every ifx run stopped here with the
      ! directory present), and a writable directory is what the run needs.
      ex_outdir = output_directory_writable()
      if (.not. ex_outdir) then
         call execute_command_line('mkdir -p output', exitstat=rc_outdir)
         ex_outdir = output_directory_writable()
         if (rc_outdir .ne. 0 .or. .not. ex_outdir) then
            write(*,*) '(init.f90) ERROR: cannot create the ./output '//  &
                       'directory the run writes to; create it and rerun.'
            error stop 1
         endif
         write(*,*) '(init.f90) Created the ./output directory.'
      endif

      !---- Global options ----!
      
      ! Set the number of threads used. The OpenMP regions are the ionization
      ! sweep (ionization_equilibrium), the rate and cooling tables
      ! (util_ion_eq), the carrier residual and Jacobian
      ! (diffusive_photochemistry), the reconstructions (PLM_rec,
      ! Reconstruction), the conserved/primitive conversions (UW_conversions),
      ! the marching right-hand side (RK_rhs) and the stationary hydrodynamic
      ! rows (hydrodynamic_rows). Measured 2026-09-28 on 2 x Xeon Gold 6154
      ! (36 cores, 72 hardware threads), gfortran 16.2 -O3: the marching of
      ! wasp_full (300 steps) speeds up 4.3x at 8 threads, 5.8x at 18 and
      ! 6.3x at 36, and is slower at 72 (4.8x); one stationary outer pass of
      ! a 500-cell molecular state speeds up 6.5x at 8, 12.2x at 18, 17.5x
      ! at 36 and 17.1x at 72 (LHS1140b/models/.P1/omp_scaling_lart2).
      ! A second hardware thread on a core adds nothing (72 against 36), so
      ! the default is the number of physical cores. Policy: honor an
      ! explicit OMP_NUM_THREADS (any value); otherwise default to the
      ! physical cores available (physical_core_count).
      ! Dynamic team adjustment is switched OFF so that the team size is the
      ! one requested and the threadprivate scratch of the sweeps persists
      ! across regions (OpenMP guarantees that only for a fixed team), and
      ! the size actually obtained is reported. The thread pool of the
      ! LAPACK library is a separate matter (blas_thread_policy).
      call get_environment_variable('OMP_NUM_THREADS', omp_env, status=omp_env_st)
      if (omp_env_st .eq. 0 .and. len_trim(omp_env) .gt. 0) then
         n_omp_threads = omp_get_max_threads()          ! user set it explicitly
      else
         n_omp_threads = physical_core_count()
      endif
      call omp_set_dynamic(.false.)
      call omp_set_num_threads(n_omp_threads)
      n_omp_obtained = 0
!$omp parallel
!$omp single
      n_omp_obtained = omp_get_num_threads()
!$omp end single
!$omp end parallel
      if (omp_env_st .eq. 0 .and. len_trim(omp_env) .gt. 0) then
         write(*,'(A,I3,A)') '    - Using', n_omp_threads,                   &
                             ' OMP threads (from OMP_NUM_THREADS)'
      else
         write(*,'(A,I3,A)') '    - Using', n_omp_threads,                   &
                             ' OMP threads (default: the physical cores)'
      endif
      if (n_omp_obtained .ne. n_omp_threads)                                  &
         write(*,'(A,I0,A)') '      (team obtained: ', n_omp_obtained,        &
                             ' threads)'
      call blas_threads_set_policy()
      call blas_threads_report(6)
      
      ! Loop parameters
      marching_step = 0  
      du  = 1.0
	   dtu = 1.0
      
      !------------------------------------------------!
      
      ! Construction of radial grid
      write(*,*) '    - Constructing the spatial grid..'
      call define_grid         

      ! Eddy diffusion coefficient of every cell. Needs the grid, and needs
      ! he_kzz / the lower-atmosphere profile to be final, so it sits here
      ! rather than in allocate_grid_arrays. Without a profile every cell
      ! takes the scalar he_kzz.
      call eddy_diffusion_on_grid
      
      !------------------------------------------------!
      
      ! Construction of energy grid
      write(*,*) '    - Precalculating energy grid, cross section and flux vector..'
      call set_energy_vectors      
      
      !------------------------------------------------!
      
      ! Pre-evaluate gravity at grid center and edges
      write(*,*) '    - Precalculating the gravitational potential..'
      call set_gravity_grid
      
      !------------------------------------------------!
      
      !---- Initial conditions ----!
      
      if (ic_mode .eq. 4 .and. .not. do_load_IC) then

         ! In-process Wind-AE warm-start IC ("IC mode: windae"): the ported
         ! Wind-AE solver builds an IC on the EXHALE grid and writes
         ! output/*_IC.txt, which load_IC then ingests -- the whole
         ! EXHALE-input -> Wind-AE solve -> IC -> run is one invocation.
         write(*,*) '    - Generating Wind-AE warm-start IC in-process..'
         call wae_generate_ic()
         write(*,*) '    - Loading the generated Wind-AE IC..'
         call load_IC(rho,v,p,T,f_sp,W)
         marching_step = 1

      else if (.not. do_load_IC) then

         ! Set IC to isothermal atmosphere
         write(*,*) '    - Setting the default, isothermal IC..'
      	call set_IC(W,T,f_sp)

      else  ! Load existing initial conditions

         write(*,*) '    - Loading IC from file..'
	      ! Load thermodynamic profiles
	      call load_IC(rho,v,p,T,f_sp,W)

	      ! A seed mapped onto this grid from another one may have its
	      ! lower layer projected onto the discrete hydrostatic equilibrium
	      ! of this grid (measurement key, off unless set).
	      call seed_layer_in_discrete_equilibrium(W)

	      ! Change starting loop counting index
	      marching_step = 1
    	endif
      
      !------------------------------------------------!

      ! Evaluate the composition of the initial state BEFORE the ghosts are
      ! filled. The lower boundary reads n_part_cell1 (the cell-1 particle
      ! count) to turn cell 1's pressure into the temperature and the particle
      ! count per unit mass its compatibility relation needs; that global is
      ! written by get_species_densities, so without this call the first
      ! Apply_BC of a run builds a face state from the input_read placeholder
      ! n_part_cell1 = ntot_bc + dp_bc, i.e. from a DIFFERENT state than the one
      ! it bounds. Every later Apply_BC in the marching loop is preceded by a
      ! composition solve, so this is the one entry point where the value can be
      ! stale, and it is the value the standalone residual/Newton diagnostics
      ! (which stop right after init) see. On HD 189733 b the placeholder made
      ! the old ghost pressure disagree with the interior by 4.4%, which the
      ! base face read as a contact discontinuity.
      nhei = 0.0d0;  nheii = 0.0d0;  nheiii = 0.0d0;  nheiTR = 0.0d0
      call get_species_densities(W(1,:),f_sp,nhi,nhii,nhei,nheii,           &
                                 nheiii,nheiTR,nm,ne,n_tot)

      ! Apply BC to initial condition
      call W_to_U(W,u)
      call Apply_BC(u)

      write(*,*) '(init.f90) Done.'
      
      ! End of subroutine
      end subroutine init

      !----------------------------------!

      subroutine seed_layer_in_discrete_equilibrium(W)
      ! THE LOWER LAYER OF A MAPPED SEED PROJECTED ONTO THIS GRID'S DISCRETE
      ! HYDROSTATIC EQUILIBRIUM. Measurement key
      ! EXHALE_SEED_HYDROSTATIC_CELLS=K (2 <= K < N; off when unset); it
      ! requires "Well balanced: True", whose reconstruction and base
      ! boundary (model ..._discrete_equilibrium_v5) share the equilibrium
      ! built here (md/Update_EXHALE_stage3.md sections 36 and 37).
      !
      ! WHY. A state interpolated onto shifted cell centers
      ! (map_state_to_grid.py) is off the discrete equilibrium of the new
      ! grid by the interpolation error, about 1e-4 to 1e-3 of the pressure
      ! in a lower layer a few cells to a scale height. The Roe pressure-jump
      ! term turns that into a face mass flux of order mismatch/c_s, which
      ! at a wind Mach number of 1e-8 is 1e4 times the wind: the mapped
      ! atomic He/H 9.7 seed at 0.06 of the fiducial XUV loads at ||R|| 1.88
      ! (md/atomic_heh97_reduced_xuv_20260929_review.md, step 4).
      !
      ! WHAT IS DONE, holding the temperature and the mass fractions of
      ! every cell (p/rho of every cell kept):
      !  - cell 1 is set so that its constant-density equilibrium reaches
      !    the reservoir pressure at the level face,
      !        p_1 + rho_1 (phi_c(1) - phi_i(0)) = p_reservoir ;
      !  - cells 2 .. K so that the equilibrium mismatch of faces 1 .. K-1
      !    vanishes (Reconstruction.f90, well_balanced_face_departures),
      !        p_j+1 - p_j + rho_j+1 (phi_c(j+1) - phi_i(j))
      !                    + rho_j   (phi_i(j)   - phi_c(j)) = 0 ;
      !  - every cell above K, the upper ghosts included, is scaled by the
      !    one factor that balances face K as well (a uniform factor on rho
      !    and p scales the mismatch of every face above it by the same
      !    factor, so no join is made);
      !  - the momentum of cells 1 .. K is set to the mass flux of cell K+1,
      !    rho v = F / r^2 with F = (rho v r^2)(K+1) after the scaling, and
      !    above K the velocity is divided by the density factor, so the
      !    mass flux of the seed is kept there.
      ! The momentum flux rho v^2 is left out of the balance: where the key
      ! is meant to act it is below 1e-10 of the pressure (v of order cm/s
      ! against c_s of order km/s). The energy balance is not projected; the
      ! stationary solve that follows does that.
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: W
      ! The projection is built on a copy and adopted only when every
      ! factor, density and pressure it gives is finite and positive: the
      ! constant-density continuation of the pressure crosses zero where a
      ! cell is not thin against the scale height, and a negative or zero
      ! factor is no equilibrium of the state (it would also divide the
      ! velocity by zero above K).
      real*8, dimension(3,1-Ng:N+Ng) :: Wt
      real*8, dimension(1-Ng:N+Ng) :: fac
      character(len=32) :: env
      integer :: K, j, st
      real*8  :: chi_j, rhs_eq, den_eq, flux_K, mom_ratio

      call get_environment_variable('EXHALE_SEED_HYDROSTATIC_CELLS', env,  &
                                    status=st)
      if (st .ne. 0 .or. len_trim(env) .eq. 0) return
      read(env,*,iostat=st) K
      if (st .ne. 0 .or. K .lt. 2 .or. K .ge. N) then
         write(*,*) '(init) ERROR: EXHALE_SEED_HYDROSTATIC_CELLS takes a'//  &
                    ' cell count from 2 to N-1, not "'//trim(env)//'".'
         error stop 1
      endif
      if (.not. well_balanced) then
         write(*,*) '(init) ERROR: EXHALE_SEED_HYDROSTATIC_CELLS builds the'// &
                    ' well-balanced equilibrium and needs "Well balanced:'//  &
                    ' True".'
         error stop 1
      endif

      Wt  = W
      fac = 1.0d0
      chi_j  = Wt(3,1)/Wt(1,1)
      den_eq = chi_j + (Gphi_c(1) - Gphi_i(0))
      if (.not. (den_eq .gt. 0.0d0 .and. base_reservoir_p .gt. 0.0d0)) &
         call projection_refused(1, base_reservoir_p, den_eq)
      fac(1) = base_reservoir_p/den_eq/Wt(1,1)
      Wt(1,1) = Wt(1,1)*fac(1)
      Wt(3,1) = chi_j*Wt(1,1)
      do j = 1, K
         chi_j    = Wt(3,j+1)/Wt(1,j+1)
         rhs_eq   = Wt(3,j) - Wt(1,j)*(Gphi_i(j) - Gphi_c(j))
         den_eq   = chi_j + (Gphi_c(j+1) - Gphi_i(j))
         if (.not. (rhs_eq .gt. 0.0d0 .and. den_eq .gt. 0.0d0))            &
            call projection_refused(j+1, rhs_eq, den_eq)
         fac(j+1) = rhs_eq/den_eq/Wt(1,j+1)
         Wt(1,j+1) = Wt(1,j+1)*fac(j+1)
         Wt(3,j+1) = chi_j*Wt(1,j+1)
      enddo
      do j = K+2, N+Ng
         fac(j) = fac(K+1)
         Wt(1,j) = Wt(1,j)*fac(j)
         Wt(3,j) = Wt(3,j)*fac(j)
      enddo
      do j = 1, N+Ng
         if (.not. (fac(j) .gt. 0.0d0 .and. fac(j) .lt. huge(1.0d0) .and.  &
                    Wt(1,j) .gt. 0.0d0 .and. Wt(1,j) .lt. huge(1.0d0) .and. &
                    Wt(3,j) .gt. 0.0d0 .and. Wt(3,j) .lt. huge(1.0d0)))     &
            call projection_refused(j, Wt(1,j), Wt(3,j))
      enddo
      Wt(2,K+1:N+Ng) = Wt(2,K+1:N+Ng)/fac(K+1:N+Ng)
      flux_K = Wt(1,K+1)*Wt(2,K+1)*r(K+1)**2
      do j = 1, K
         Wt(2,j) = flux_K/(Wt(1,j)*r(j)**2)
      enddo
      W = Wt
      ! The balance leaves the momentum flux out; its size in the layer is
      ! reported, not judged (the stationary rows that follow judge it).
      mom_ratio = maxval(W(1,1:K)*W(2,1:K)**2/W(3,1:K))

      write(*,'(A,I0,A,ES11.3,A,I0,A,ES11.3)') ' (init) seed layer of ', K, &
           ' cells projected onto the discrete hydrostatic equilibrium:'//  &
           ' largest density change ', maxval(abs(fac(1:K+1) - 1.0d0)),    &
           ' (cell ', maxloc(abs(fac(1:K+1) - 1.0d0), 1),                  &
           '); cells above scaled by', fac(K+1) - 1.0d0
      write(*,'(A,ES11.3)') '   (temperature and mass fractions kept; the'// &
           ' layer carries the mass flux of cell K+1;'//                    &
           ' EXHALE_SEED_HYDROSTATIC_CELLS); largest rho v^2 / p left out'// &
           ' of the balance in the layer:', mom_ratio

      ! THE FILE'S STATIONARY CLAIM IS ABOUT THE STATE IT CARRIED, NOT THIS
      ! ONE. The projection moved density, pressure and velocity, so no
      ! certificate of the loaded pair may be written with it (an IC dump)
      ! or held against its evaluation (the original-claim check of a
      ! stationary evaluate): the state is a seed from here on.
      ic_certified   = .false.
      ic_cert_reason = 'projected_seed'

      contains

      subroutine projection_refused(jc, a, b)
      integer, intent(in) :: jc
      real*8,  intent(in) :: a, b
      write(*,'(A,I0,A,2ES12.4)') ' (init) ERROR: the hydrostatic'//      &
           ' projection of EXHALE_SEED_HYDROSTATIC_CELLS is not admissible'// &
           ' at cell ', jc, ': the constant-density pressure continuation'// &
           ' (numerator, denominator; or density, pressure) gives ', a, b
      write(*,'(A)') '   (the cell is not thin against the scale height'//  &
           ' there, or the state is not positive); nothing was adopted.'
      error stop 1
      end subroutine projection_refused

      end subroutine seed_layer_in_discrete_equilibrium

      !----------------------------------!

      integer function physical_core_count()
      ! The physical cores this process may run on: the processors of its
      ! affinity mask (Cpus_allowed_list of /proc/self/status), each one
      ! assigned to its core by the first entry of its Linux list
      ! /sys/devices/system/cpu/cpu<i>/topology/thread_siblings_list
      ! ("0,36" on a machine with two hardware threads a core), and the
      ! distinct cores counted. Where either list cannot be read, every
      ! available processor (omp_get_num_procs) counts as a core.
      integer, parameter :: n_cpu_max = 4096
      integer :: allowed(n_cpu_max), siblings(n_cpu_max), core_of(n_cpu_max)
      integer :: n_allowed, n_siblings, n_cores, k, u, ios
      character(len=4096) :: line
      character(len=96) :: path
      physical_core_count = max(1, omp_get_num_procs())
      n_allowed = 0
      open(newunit=u, file='/proc/self/status', action='read', iostat=ios)
      if (ios .ne. 0) return
      do
         read(u,'(A)', iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:18) .eq. 'Cpus_allowed_list:') then
            call cpu_list_parse(line(19:), allowed, n_allowed)
            exit
         endif
      enddo
      close(u)
      if (n_allowed .le. 0) return
      n_cores = 0
      do k = 1, n_allowed
         write(path,'(A,I0,A)') '/sys/devices/system/cpu/cpu', allowed(k),  &
                                '/topology/thread_siblings_list'
         open(newunit=u, file=trim(path), action='read', iostat=ios)
         if (ios .ne. 0) return
         read(u,'(A)', iostat=ios) line
         close(u)
         if (ios .ne. 0) return
         call cpu_list_parse(line, siblings, n_siblings)
         if (n_siblings .le. 0) return
         if (all(core_of(1:n_cores) .ne. minval(siblings(1:n_siblings)))) then
            n_cores = n_cores + 1
            core_of(n_cores) = minval(siblings(1:n_siblings))
         endif
      enddo
      physical_core_count = max(1, n_cores)
      end function physical_core_count

      !----------------------------------!

      subroutine cpu_list_parse(text, ids, n)
      ! The processor numbers of a Linux CPU list such as "0-3,8,10-11"
      ! (the format of Cpus_allowed_list and thread_siblings_list). n = 0
      ! when the text cannot be read as such a list.
      character(len=*), intent(in) :: text
      integer, intent(out) :: ids(:), n
      integer :: i, i0, a, b, dash, ios, m, last
      character(len=64) :: tok
      n = 0
      last = len_trim(text)
      i0 = 1
      do i = 1, last + 1
         if (i .le. last) then
            if (text(i:i) .ne. ',') cycle
         endif
         tok = adjustl(text(i0:i-1))
         i0 = i + 1
         if (len_trim(tok) .eq. 0) cycle
         dash = index(tok, '-')
         if (dash .gt. 0) then
            read(tok(:dash-1), *, iostat=ios) a
            if (ios .eq. 0) read(tok(dash+1:), *, iostat=ios) b
         else
            read(tok, *, iostat=ios) a
            b = a
         endif
         if (ios .ne. 0 .or. b .lt. a) then
            n = 0
            return
         endif
         do m = a, b
            if (n .ge. size(ids)) return
            n = n + 1
            ids(n) = m
         enddo
      enddo
      end subroutine cpu_list_parse

      !----------------------------------!

      logical function output_directory_writable()
      ! Whether a file can be created in ./output: a probe file is opened
      ! with status='replace' and deleted on close. Portable where an
      ! INQUIRE on the directory itself is not (see the caller).
      integer :: u, ios
      open(newunit=u, file='output/.writable_probe', status='replace',   &
           action='write', iostat=ios)
      output_directory_writable = (ios .eq. 0)
      if (ios .eq. 0) close(u, status='delete')
      end function output_directory_writable
      
      ! End of module
      end module Initialization
