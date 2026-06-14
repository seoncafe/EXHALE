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
      use wae_exhale_bridge, only: wae_generate_ic
      
      implicit none 
      
      contains
      
      subroutine init(W,u,f_sp)
      ! Initialize the simulation setup

      integer :: n_omp_threads, omp_env_st
      character(len=32) :: omp_env
      real*8, dimension(1-Ng:N+Ng)   :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng,n_species),intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng,3),intent(out) :: W,u    
      
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

      ! Apply BC to initial condition
      call W_to_U(W,u)
      call Apply_BC(u,u)    

      write(*,*) '(init.f90) Done.'
      
      ! End of subroutine
      end subroutine init
      
      ! End of module
      end module Initialization
