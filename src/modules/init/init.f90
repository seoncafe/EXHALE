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
      use IC_load
      use initial_conditions
      use composition, only: get_species_densities
      use species_table, only: n_mion
      use wae_exhale_bridge, only: wae_generate_ic
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      
      implicit none 
      
      contains
      
      subroutine init(W,u,f_sp)
      ! Initialize the simulation setup

      integer :: n_omp_threads, omp_env_st
      character(len=32) :: omp_env
      real*8, dimension(1-Ng:N+Ng)   :: rho,v,p,T
      ! composition scratch: only n_part_cell1 (set inside get_species_densities)
      ! is wanted here, but the single policy point returns the whole set
      real*8, dimension(1-Ng:N+Ng)   :: nhi,nhii,nhei,nheii,nheiii,nheiTR
      real*8, dimension(1-Ng:N+Ng)   :: ne,n_tot
      real*8, dimension(1-Ng:N+Ng,n_mion) :: nm
      real*8, dimension(1-Ng:N+Ng,n_species),intent(out) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng),intent(out) :: W,u
      
      write(*,*) '(init.f90) Initializing the simulation..'

      !---- Global options ----!
      
      ! Set the number of threads used. The current OpenMP coverage is limited
      ! (only the radiation rate loop + the ionization loop are parallel), so a
      ! ~500-cell grid scales out by ~16 threads and beyond that the fork/join
      ! overhead makes it SLOWER (measured: 16 thr ~ 59 steps/s, 60 thr ~ 56).
      ! Policy: honor an explicit OMP_NUM_THREADS (any value); otherwise default
      ! to min(cores, 16) instead of grabbing every core for ~no gain.
      call get_environment_variable('OMP_NUM_THREADS', omp_env, status=omp_env_st)
      if (omp_env_st .eq. 0 .and. len_trim(omp_env) .gt. 0) then
         n_omp_threads = omp_get_max_threads()          ! user set it explicitly
      else
         n_omp_threads = min(omp_get_max_threads(), 16) ! sensible default
      endif
      call omp_set_num_threads(n_omp_threads)
      write(*,'(A12,I3,A21)') '    - Using',n_omp_threads,                   &
                              ' OMP threads (default)'
      if (omp_env_st .eq. 0 .and. len_trim(omp_env) .gt. 0)                  &
         write(*,'(A)') '      (from OMP_NUM_THREADS)'
      
      ! Loop parameters
      count = 0  
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
         count = 1

      else if (.not. do_load_IC) then

         ! Set IC to isothermal atmosphere
         write(*,*) '    - Setting the default, isothermal IC..'
      	call set_IC(W,T,f_sp)

      else  ! Load existing initial conditions

         write(*,*) '    - Loading IC from file..'
	      ! Load thermodynamic profiles
	      call load_IC(rho,v,p,T,f_sp,W)

	      ! Change starting loop counting index
	      count = 1
    	endif
      
      !------------------------------------------------!

      ! Evaluate the composition of the initial state BEFORE the ghosts are
      ! filled. Apply_BC reads n_part_cell1 (the cell-1 particle count) for the
      ! continuous-temperature base ghost, T(1) = p(1)/n_part_cell1; that global
      ! is written by get_species_densities, so without this call the first
      ! Apply_BC of a run uses the input_read placeholder n_part_cell1 =
      ! ntot_bc + dp_bc, i.e. a ghost built from a DIFFERENT state than the one
      ! it bounds. Every later Apply_BC in the marching loop is preceded by a
      ! composition solve, so this is the one entry point where the value can be
      ! stale, and it is the value the standalone residual/Newton diagnostics
      ! (which stop right after init) see. On HD 189733 b the placeholder made
      ! the ghost pressure disagree with the interior by 4.4%, which the base
      ! face reads as a contact discontinuity (docs/hd189_base_checkerboard.md
      ! section 15). Runs that do not set "Base ghost temperature: continuous"
      ! never read n_part_cell1, so they are unaffected.
      nhei = 0.0d0;  nheii = 0.0d0;  nheiii = 0.0d0;  nheiTR = 0.0d0
      call get_species_densities(W(1,:),f_sp,nhi,nhii,nhei,nheii,           &
                                 nheiii,nheiTR,nm,ne,n_tot)

      ! Same rule for the other state-dependent quantity the lower BC reads:
      ! with "Base velocity: massflux" the ghost velocity is F_c/(rho_bc r^2),
      ! and F_c is the wind mass-flux constant that the marching loop updates
      ! each step. It starts at -1, which Apply_BC reads as "not available yet"
      ! and silently falls back to the valve -- so the first step, and any
      ! diagnostic that stops right after init, evaluate a DIFFERENT velocity
      ! boundary condition than the run uses. Seed it from the initial state,
      ! by the same average over the escape region the loop takes.
      if (base_v_massflux .and. j_min .le. N)                                &
         base_flux_const = sum(W(1,j_min:N)*W(2,j_min:N)                     &
                               *r(j_min:N)*r(j_min:N))/dble(N - j_min + 1)

      ! Apply BC to initial condition
      call W_to_U(W,u)
      call Apply_BC(u)

      write(*,*) '(init.f90) Done.'
      
      ! End of subroutine
      end subroutine init
      
      ! End of module
      end module Initialization
