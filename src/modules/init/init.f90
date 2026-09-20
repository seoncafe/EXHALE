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
      ! rows (hydrodynamic_rows): a ~500-cell grid scales out by ~16 threads
      ! and beyond that the fork/join overhead makes it SLOWER (measured
      ! 2026-07 with fewer regions: 16 thr ~ 59 steps/s, 60 thr ~ 56; not
      ! remeasured with the present coverage). Policy: honor an explicit
      ! OMP_NUM_THREADS (any value); otherwise default to min(cores, 16).
      ! Dynamic team adjustment is switched OFF so that the team size is the
      ! one requested and the threadprivate scratch of the sweeps persists
      ! across regions (OpenMP guarantees that only for a fixed team), and
      ! the size actually obtained is reported. The thread pool of the
      ! LAPACK library is a separate matter (blas_thread_policy).
      call get_environment_variable('OMP_NUM_THREADS', omp_env, status=omp_env_st)
      if (omp_env_st .eq. 0 .and. len_trim(omp_env) .gt. 0) then
         n_omp_threads = omp_get_max_threads()          ! user set it explicitly
      else
         n_omp_threads = min(omp_get_max_threads(), 16) ! sensible default
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
                             ' OMP threads (default: min(cores, 16))'
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
      ! base face read as a contact discontinuity
      ! (docs/hd189_base_checkerboard.md section 15).
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
