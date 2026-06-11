      module global_parameters
      ! Definition of global parameters and vectors
      
      implicit none
      
      integer, parameter :: outfile = 99  ! Unit number of report file
      integer, parameter :: N = 500       ! Number of computational cells
      integer, parameter :: Ng = 2        ! Number of ghost cells
      integer, parameter :: Nl_fix = 200  ! Number of default energy bins
      ! Number of ion-fraction species carried in f_sp:
      !  1-6   HI,HII,HeI,HeII,HeIII,HeITR
      !  7-9   CI,CII,CIII    10-12 OI,OII,OIII    13-15 NI,NII,NIII
      !  16-18 MgI,MgII,MgIII 19-21 SiI,SiII,SiIII  22-24 CaI,CaII,CaIII
      !  25-26 NaI,NaII       27-28 KI,KII          29-30 SI,SII
      !  31-33 FeI,FeII,FeIII
      ! (Si/Ca/Fe carry three stages like Mg; Na/K/S carry two stages.)
      integer, parameter :: n_species = 33
      integer :: Nl, NlTR                 ! Number of points for energy integrations
      integer :: N_eq                     ! Numbers of equations in NL solver
      integer :: lwa                      ! Working array length for NL solver
      integer :: info                     ! Output info variable of NL solver
      integer :: j_min
      integer :: count
      
      character(len = 9), parameter   :: inp_file = 'input.inp'
      character(len = :), allocatable :: p_name
      character(len = :), allocatable :: grid_type
      character(len = :), allocatable :: flux 
      character(len = :), allocatable :: rec_method 
      character(len = :), allocatable :: appx_mth
      character(len = :), allocatable :: sp_type
      character(len = :), allocatable :: sed_file
      
      logical :: is_mom_const  = .false.  ! Is momentum constant within tolerance
      logical :: is_zero_dt    = .false.  ! Is time derivative really zero 
      logical :: force_start   = .false.  ! Force to do first 1000 iterations
      logical :: do_only_pp    = .false.  ! Do only the post processing
      logical :: do_load_IC    = .false.  ! Load existing IC
      logical :: thereis_He    = .false.  ! Is He included in computations
      logical :: do_read_sed   = .false.  ! Is numerical SED read from file
      logical :: is_monochr    = .false.  ! Is monochromatic radiation selected
      logical :: is_PL_sed     = .false.  ! Is the SED a power law
      logical :: thereis_Xray  = .false.  ! Include only EUV band	
      logical :: thereis_HeITR = .false.  ! Include calculations for He triplet
      logical :: thereis_metals = .false. ! Include trace-metal species
      ! EOS electron/particle policy (composition module). .false. = legacy
      ! behavior: ne and n_tot count H/He only (metals enter ionization but
      ! not the bulk gas EOS). .true. = fully-coupled metal EOS (Phase 3
      ! hook; currently unused, keep .false. for byte-identical results).
      logical :: eos_include_metals = .false.
                                          !  (C/N/O/Mg/Si/Ca/Na/K/S)
      logical :: thereis_lowIP_metal = .false. ! An active metal whose neutral
                                          !  ionization potential lies below the
                                          !  13.6 eV HI edge (e.g. Mg I, 7.646
                                          !  eV) is present; triggers extension
                                          !  of the energy grid/SED below 13.6 eV
      logical :: use_2lev_cool  = .false. ! Two-level fine-structure metal
                                          !  cooling ([O I] 63um, [C II] 158um)
      integer :: pp_metal_mode  = 1       ! Metal treatment in the advection
                                          !  post-process (post_process_adv):
                                          !  0 = metal-free (legacy: metals
                                          !      zeroed in _adv heating/cooling);
                                          !  1 = frozen (use converged eq metal
                                          !      densities, held fixed) (default);
                                          !  2 = re-solve metal ionization at the
                                          !      post-process T/n_e (not yet
                                          !      implemented; warns and falls
                                          !      back to mode 1).
                                          !  Default is 1 (frozen): the frozen
                                          !  _adv profile already tracks the eq
                                          !  solution closely, so the re-solve
                                          !  would change little.
                                          !  Set at runtime via 'pp_metals' in
                                          !  metals.inp. No effect when metals
                                          !  are off (nm = 0 either way).
      logical :: use_weno3     = .false.  ! Use WENO3 reconstruction
      logical :: use_plm       = .false.  ! Use PLM reconstruction
      logical :: is_stalled    = .false.  ! Convergence stalled at a du plateau
      logical :: spherical_domain = .false. ! Domain mode (set via input.inp):
                                          !  .false. = full Roche potential
                                          !   truncated at the Hill/L1 radius
                                          !   (default ATES behavior);
                                          !  .true.  = pure planetary potential
                                          !   (-b0/r), tidal + centrifugal terms
                                          !   dropped, domain extended to
                                          !   r_out_user [R_p] (Huang Case A-like)
      logical :: transonic_ic  = .false.  ! IC type (set via input.inp,
                                          !  "Transonic IC: True"; default off):
                                          !  .false. = isothermal hydrostatic
                                          !   atmosphere + v~0 seed (default);
                                          !  .true.  = transonic isothermal-wind
                                          !   IC (Parker-type, solved from the
                                          !   ATES potential via the Bernoulli
                                          !   integral). Needed for deep-RLOF
                                          !   bases (sonic point near L1) that
                                          !   launch no wind from the hydrostatic
                                          !   IC. Auto-enables force_start so the
                                          !   cold wind heats to its hot steady
                                          !   state before the convergence test
                                          !   (the IC already has rho*v*r^2=const,
                                          !   so du~0 would exit prematurely).
                                          !  Falls back to the hydrostatic IC if
                                          !   no interior sonic point exists.

      ! Hot-Parker IC (set via input.inp "Hot Parker IC: <Twind[K]>"; default
      ! off): a transonic isothermal-wind IC built at a WARM effective wind
      ! temperature T_wind_ic (instead of the cold base T0), with T and H/He
      ! ionization ramped from the cold neutral base to the hot ionized wind, so
      ! v, rho, T, and ionization all start close to the final escaping-wind state
      ! (a p-winds-style Parker initial guess). Auto-enables force_start.
      logical :: hot_parker_ic = .false.
      real*8  :: T_wind_ic     = 1.0d4    ! effective wind temperature [K] for the IC
      ! Cold-base transition radius [R_p] for the hot-Parker IC: below it the IC
      ! is forced to the cold, static, cold-pressure base state (matching the
      ! lower BC), above it the warm Parker wind is ramped in. This removes the
      ! inverted base pressure gradient that otherwise drives a base inflow and
      ! a thermal runaway (see docs/initial_condition_benchmark, root-cause).
      real*8  :: hp_base_rtr   = 1.3d0    ! transition radius [R_p]; >1

      ! Nonlinear-solver selection (set via input.inp; defaults = the upgraded
      ! solvers). use_newton_ieq: ionization equilibrium via the analytic-
      ! Jacobian Newton (+hybrd1 fallback) vs. legacy MINPACK hybrd1
      ! ("Newton solver: False"). use_brent_tsolve: the post-process energy
      ! equation via the bracketing Brent root-finder vs. legacy hybrd1 + the
      ! 2x-band reject ("Brent solver: False"). See docs/Update_ATES_solver.
      logical :: use_newton_ieq  = .true.
      logical :: use_brent_tsolve = .true.

      !------- Global constants -------!
      
      ! Physical constants
      real*8,parameter ::  pi      = 3.1415926536     ! pi   
      real*8,parameter ::  kb_erg  = 1.38e-16         ! Boltzmann constant in CGS units
      real*8,parameter ::  kb_eV   = 8.6167e-05       ! Boltzmann constant (eV/K)
      real*8,parameter ::  mu      = 1.673e-24        ! Hydrogen mass (g)
      real*8,parameter ::  g       = 1.666666666667   ! Polytropic index
      real*8,parameter ::  Gc      = 6.67259e-8       ! Gravitational constant (CGS)   
      real*8,parameter ::  erg2eV  = 6.241509075e11   ! 1 erg measured in eV
      real*8,parameter ::  hp_erg  = 6.62620e-27      ! Planck constant in CGS units
      real*8,parameter ::  hp_eV   = 4.1357e-15       ! Planck constant (eV*s)
      real*8,parameter ::  c_light = 2.99792458e10    ! Speed of light in cm/s
      real*8,parameter ::  parsec  = 3.08567758147e18 ! 1 pc in cm
      real*8,parameter ::  AU      = 1.495978707e13   ! Astronomical unit
      real*8,parameter ::  RJ      = 6.9911e9         ! Jupiter radius (cm)
      real*8,parameter ::  MJ      = 1.898e30         ! Jupiter mass (g)
      real*8,parameter ::  Msun    = 1.989e33         ! Sun mass (g)
      real*8,parameter ::  Rsun    = 6.957e10         ! Sun radius (cm)
      real*8,parameter ::  R_earth = 6.3725e8         ! Earth radius (cm)
      real*8,parameter ::  M_earth = 5.9726e27        ! Earth mass (g)
      real*8,parameter ::  ih      = 1.0              ! Atomic number of Hydrogen
      real*8,parameter ::  ihe     = 2.0              ! Atomic number of Helium
	
	   ! -- - Energy constants

      ! Energy intervals
      real*8 ::  e_top
      real*8 ::  e_mid
      real*8 ::  e_low
      
      ! Threshold energies
      real*8,parameter ::  e_th_HI   = 13.6      ! Threshold for HI ionization
      real*8,parameter ::  e_th_HeI  = 24.6      ! Threshold for HeI ionization
      real*8,parameter ::  e_th_HeII = 54.4      ! Threshold for HeII ionization
      real*8,parameter ::  e_th_HeTR = 4.80      ! Threshold for HeI triplet ionization
      real*8,parameter ::  e_th_MgI  = 7.646     ! Threshold for MgI ionization
      real*8,parameter ::  e_th_MgII = 15.035    ! Threshold for MgII ionization

	   ! Numerical constants
      real*8 ::  CFL    = 0.6         ! CFL number; settable in input.inp via "CFL:" (lower = smaller dt, may damp a numerical limit cycle)
      real*8 ::  du_th     = 1.0d-3   ! final (stage-2 / WENO3) escape-momentum threshold; settable in input.inp via "du_th [PLM,WENO3]:". 1e-3 is the original ATES-Code-main value (ATES-metal had loosened it to 2e-2, accepting ~2% mass-flux spread).
      real*8 ::  du_th_plm = -1.0d0   ! stage-1 (PLM) threshold; if > du_th the run is two-stage: PLM until du<du_th_plm, then switch reconstruction to WENO3 and converge at du<du_th. <=0 => single-stage at du_th.
      real*8,parameter ::  dtu_th = 1.0d-8      ! Threshold variation of time deriv.
      real*8           ::  du                   ! Initial momentum variation
      real*8           ::  dtu                  ! Norm of time derivative

      ! Convergence stall detection (ported from ATES_extended): stop runs
      ! whose du settles on a plateau above du_th (steady state reached but
      ! the strict thresholds physically unreachable for the given setup).
      integer,parameter :: count_max = 1000000  ! hard cap on iterations
      integer :: N_stall   = 2000     ! window of stalled steps; settable via
                                      !   input line "Stall [tol,N]: <tol> <N>"
      real*8  :: stall_tol = 1.0d-6   ! rel. du change defining a stall (same line)

      ! Mass-flux LEVEL stability tolerance: in stage 2 a run only counts as
      ! converged/stalled when the mean |rho*v*r^2| over [j_min:N] changed by
      ! less than lev_th (relative) across the last N_stall steps. du is the
      ! SPREAD of the mass flux and is blind to a uniform drift of the whole
      ! profile; without this gate, runs can stop on spatially flat but still
      ! level-drifting states that differ by tens of percent in Mdot
      ! (path-dependent quasi-steady snapshots; see refactor_stage1_log.md).
      ! "Level tol: <val>" in input.inp; <= 0 disables (DEFAULT, legacy
      ! behavior preserved byte-identically). Kept opt-in until a steady-
      ! state Newton/PTC reference calibrates the production tolerance: the
      ! WASP-121b study showed the true converged state needs >35k marching
      ! steps (du dips below 1e-3 long before, at a premature ~0.08-dex-low
      ! snapshot), so a default-on level gate would change every result
      ! without a validated "truth" to anchor the tolerance to.
      real*8  :: lev_th = -1.0d0

      ! Production steady-solver wiring ("Solver: Newton [R_switch]").
      ! When .true., the normal marching loop (including the automatic
      ! two-stage PLM->WENO3) runs as a WARM-UP; once the periodic steady-
      ! residual monitor reports max_k ||R_k|| < newton_R_switch, the JFNK
      ! steady solver finishes the run to ||R|| < resid_th (default 1e-3 if
      ! "Resid tol:" was not given) and the standard final outputs /
      ! post-processing follow. The switch is residual-based, NOT du-based:
      ! du dips transiently while the energy residual is still O(1-30)
      ! (premature-dip lesson), far outside Newton's basin; warm starts with
      ! ||R|| ~ 2e-2 converge robustly, hence the 5e-2 default. Requires a
      ! smooth base valve: if "Valve eps:" was not set, 1e-4 is adopted.
      logical :: use_newton_solver = .false.
      real*8  :: newton_R_switch   = 5.0d-2

      ! WENO3 weight-freezing mode for the steady Newton solver (lagged /
      ! frozen nonlinear weights, the standard FV steady-solve remedy):
      !   0 = compute weights fresh every call (DEFAULT, byte-identical)
      !   1 = compute AND store the smoothness factors S0/S1
      !   2 = reuse the stored factors (linearizes the reconstruction so the
      !       inner Newton problem is mildly nonlinear; refreshed per outer
      !       iteration by a mode-1 evaluation)
      integer :: weno_mode = 0

      ! Base-valve smoothing ("Valve eps: <v_eps>", code velocity units).
      ! The lower-BC ghost velocity is v_g = max(v_1, 0) (one-way valve); the
      ! kink at v_1 = 0 sits exactly where the breathing base lives (steady
      ! states hug v_1 -> 0+), making the steady residual non-differentiable
      ! there and blocking Newton line searches. With valve_eps > 0 the valve
      ! becomes the softplus 0.5*(v + sqrt(v^2 + eps^2)) - 0.5*eps, which is
      ! smooth, -> max(v,0) as eps -> 0, and = v - O(eps^2/v) for v >> eps.
      ! Default <= 0 keeps the exact legacy valve (byte-identical).
      real*8  :: valve_eps = -1.0d0

      ! Residual-based convergence ("Resid tol: <val>"). When > 0, the run
      ! converges when the finite-volume steady residual
      !   R = dF - S        (mass, momentum)
      !   R = dF_E - S_E - (heat - cool)   (energy)
      ! has max_j |R(j,k)|/max|u(:,k)| over [j_min:N] below resid_th, REPLACING
      ! the du<du_th test (du measures only the mass-flux spread and is blind
      ! to an operator-split energy imbalance: at the premature WASP golden,
      ! du/dtu were tiny while the energy residual was ~30). Evaluated every
      ! N_resid steps (cheap: one Reconstruct + RK_rhs, no extra ioniz_eq).
      ! Default <= 0 disables (legacy du-based stop, byte-identical).
      real*8  :: resid_th = -1.0d0
      integer :: N_resid  = 500

      ! Energy source-term integrator: .true. = semi-implicit backward-Euler
      ! cell solve (ATES-metal default); .false. = original explicit forward
      ! Euler ("Energy solver: Explicit" in input.inp). Runtime switch kept for
      ! solver-component isolation tests.
      logical :: use_semi_implicit_energy = .true.

      ! Local (per-cell) pseudo-time stepping ("Time stepping: Local" in
      ! input.inp). Steady-state acceleration only: each cell marches with its
      ! own CFL step dt_j = CFL*dr_j/(|v|+cs), removing the global-dt
      ! bottleneck set by the smallest base cell. Time accuracy is lost but
      ! the steady fixed point (dF = S, heating = cooling) is unchanged.
      ! Default .false. = original global-dt marching, byte-identical.
      logical :: use_local_dt = .false.
      
      !------- Planetary parameters -------!

      real*8  ::  n0          ! Density at lower boundary (cm^-3)
      real*8  ::  R0          ! Planetary radius (cm)
      real*8  ::  Mp          ! Planet mass (g)
      real*8  ::  T0          ! Temperature at lower boundary (K)
      real*8  ::  LX          ! Log10 X-Ray luminosity (erg/s)
      real*8  ::  LEUV        ! Log10 EUV luminosity (erg/s)
      real*8  ::  J_XUV       ! Bolometric star XUV flux (erg/cm^2*s)
      real*8  ::  Mstar       ! Mass of companion star
      real*8  ::  Mrapp       ! Ratio M_star/M_p
      real*8  ::  atilde      ! Orbital radius in unit of R0 (a/R0)
      real*8  ::  r_esc       ! Escape radius for constant momentum
      real*8  ::  a_orb       ! Orbital distance
      real*8  ::  r_max       ! Maximum radius = Roche Lobe dimension
      real*8  ::  r_out_user  ! User-set outer radius [R_p], spherical mode only
      real*8  ::  HeH         ! He/H ratio
      real*8  ::  X_C         ! C/H number ratio (solar ~ 2.69e-4)
      real*8  ::  X_N         ! N/H number ratio (solar ~ 6.76e-5)
      real*8  ::  X_O         ! O/H number ratio (solar ~ 4.90e-4)
      real*8  ::  X_Mg        ! Mg/H number ratio (solar ~ 3.98e-5)
      real*8  ::  X_Si        ! Si/H number ratio (solar ~ 3.24e-5)
      real*8  ::  X_Ca        ! Ca/H number ratio (solar ~ 2.19e-6)
      real*8  ::  X_Na        ! Na/H number ratio (solar ~ 1.74e-6)
      real*8  ::  X_K         ! K/H  number ratio (solar ~ 1.07e-7)
      real*8  ::  X_S         ! S/H  number ratio (solar ~ 1.32e-5)
      real*8  ::  X_Fe        ! Fe/H number ratio (solar ~ 3.16e-5)
      real*8, allocatable :: melem_ab(:) ! per-element abundance, canonical
                                          !  element order (iel_*); set in
                                          !  input_read so metal code can index
                                          !  abundance by element, not by name
      real*8  ::  rho_bc      ! Adimensional number density at origin
      real*8  ::  a_tau       ! Rate correction coefficient
      real*8  ::  PLind       ! Index of spectral power law
      real*8  ::  dp_bc       ! Boundary condition for pressure
	
      !------- Normalizations -------!

      real*8 :: v0      ! Velocity normalization
      real*8 :: t_s     ! Time normalization
      real*8 :: p0      ! Pressure normalization
      real*8 :: q0      ! Scale normalization 
      real*8 :: b0      ! Jeans parameter at planet surface

      !------- Global vectors -------!
      
      real*8, dimension(:), allocatable :: e_v, de_v
      real*8, dimension(:), allocatable :: s_hi,s_hei,s_heii,s_heiTR
      ! Metal photoionization cross sections, one column per photo-ionizable
      ! metal ion in species_table iphot order (1=CI,2=CII,3=OI,4=OII,
      ! 5=NI,6=NII,7=MgI,8=MgII,9=SiI,10=SiII,11=CaI,12=CaII,13=NaI,14=KI,
      ! 15=SI). Replaces the former loose s_ci..s_mgii.
      real*8, dimension(:,:), allocatable :: sigma_tab

      !------- Opacity-model dispatcher (set via opacity.inp) -------!
      ! See src/modules/radiation/opacity_models.f90 for the dispatcher.
      character(len = 1) :: opacity_model = 'A'  ! 'A','C','P','T'
      ! Constant-model multiplicative factors (default 1.0):
      real*8 :: opa_const_HI    = 1.0d0
      real*8 :: opa_const_HeI   = 1.0d0
      real*8 :: opa_const_HeII  = 1.0d0
      real*8 :: opa_const_HeITR = 1.0d0
      ! Robinson & Catling pressure broadening f(p)=1+a*(p/p_pivot)^n:
      real*8 :: opa_pb_factor   = 0.0d0    ! a
      real*8 :: opa_pb_exponent = 1.0d0    ! n
      real*8 :: opa_pb_pivot    = 1.0d5    ! dyne/cm^2 = 0.1 bar
      ! Tabulated-model (.atesopa) per-species file paths:
      character(len = :), allocatable :: opa_file_HI
      character(len = :), allocatable :: opa_file_HeI
      character(len = :), allocatable :: opa_file_HeII
      character(len = :), allocatable :: opa_file_HeITR
      real*8, dimension(:), allocatable :: F_XUV
      ! Per-cell pressure-broadening multiplier for the opacity ('P'
      ! model). 1.0 everywhere unless opacity_model='P'; set before each
      ! photoheating call and used to weight the opacity column density.
      real*8, dimension(1-Ng:N+Ng) :: opa_pf = 1.0d0
      real*8, dimension(1-Ng:N+Ng) :: r,r_edg,dr_j
      real*8, dimension(1-Ng:N+Ng) :: Gphi_c,Gphi_i

      !------- Phase 3a: excited hydrogen H(n=2) coupling -------!
      ! Christie+2013 / Huang+2017 n=2 (2s/2p) model feeding back into
      ! the H ionization balance and the energy equation. All terms are
      ! ZERO unless use_excited_H = .true. (enabled at runtime when a
      ! stellar T_eff is supplied in input.inp), so the default build is
      ! byte-identical to the Phase-2 result. See excited_hydrogen.f90.
      logical :: use_excited_H   = .false. ! master switch (set if T_star_eff>0)
      logical :: incl_deexc_heat = .false. ! add collisional de-excitation
                                           !  heating (overlaps the existing HI
                                           !  coex cooling; off by default)
      real*8  :: T_star_eff = 0.0d0        ! stellar effective temperature [K]
                                           !  (<=0 disables the Balmer continuum)
      real*8  :: R_star     = 0.0d0        ! stellar radius [cm] (Balmer dilution)
      real*8  :: gamma2_bal = 0.0d0        ! n=2 Balmer photoion. rate [s^-1]
      real*8  :: hpe2_bal   = 0.0d0        ! photoelec. heat per n=2 atom [erg/s]
      ! Ly-alpha mean-intensity source for the n=2 pumping term:
      !  jlya_mode = 0 -> parameterized J_Lya = 0.1*F_LyC/Dnu_D with a
      !                    top-down 1/(1+tau_Lya) line-center attenuation;
      !  jlya_mode = 1 -> read a J_Lya(r) profile from jlya_rt_file (an
      !                    external RT result, e.g. a real Monte Carlo) directly;
      !  jlya_mode = 2 -> in-line escape-probability RT (Neufeld/Harrington wing
      !                    escape), computed every timestep from the current state
      !                    in excited_H_update (see lya_rt.f90).
      integer :: jlya_mode  = 0
      character(len=200) :: jlya_rt_file = 'jlya_rt.txt'
      ! ----- Phase 3b: Ly-alpha radiative transfer ----- !
      !  F_Lya_star -> incident stellar Ly-alpha flux at the planet
      !  [erg cm^-2 s^-1] (the stellar source of J_lya in jlya_mode 2; Case D
      !  scales it x0.35). The "cooling trapping" is realized through the Phase-3a
      !  H(n=2) heating channels (photoelectric + collisional de-excitation), as in
      !  Huang et al. 2023, NOT a separate escape-probability factor on the Ly-alpha
      !  cooling. OFF by default => byte-identical to Phase 3a.
      real*8  :: F_Lya_star = 0.0d0
      ! Stellar Ly-alpha line half-width [km/s] (broad plateau ~+-70 km/s,
      ! Huang+2017); sets how deep the (trapped) stellar beam penetrates via the
      ! line wings. Editable in input.inp ("Lya stellar halfwidth [km/s]:").
      real*8  :: dv_star_lya = 70.0d0
      ! Bounded stellar trapping buildup: the penetrating stellar Ly-alpha photons
      ! scatter and build the mean intensity above free-streaming by E=min(boost,
      ! 1/beta) (-> 1 in the thin outer wind, capped at lya_star_boost in the thick
      ! region; the full 1/beta over-counts). Tuned to Huang+2023 Fig. 11.
      ! Editable: "Lya stellar boost [-]:".
      real*8  :: lya_star_boost = 5.0d0
      ! Per-cell feedback arrays injected into ioniz_eq (zero unless enabled):
      real*8, dimension(1-Ng:N+Ng) :: gph_balmer_HI = 0.0d0 ! extra HI photoion [s^-1]
      real*8, dimension(1-Ng:N+Ng) :: heat_balmer   = 0.0d0 ! extra heat [erg cm^-3 s^-1]
      ! Per-cell diagnostics for output/Excited_H.txt:
      real*8, dimension(1-Ng:N+Ng) :: Jlya_arr  = 0.0d0 ! Lya mean intensity [cgs]
      real*8, dimension(1-Ng:N+Ng) :: n2s_arr   = 0.0d0 ! H(2s) density [cm^-3]
      real*8, dimension(1-Ng:N+Ng) :: n2p_arr   = 0.0d0 ! H(2p) density [cm^-3]
      real*8, dimension(1-Ng:N+Ng) :: Sproton_arr = 0.0d0 ! Balmer proton src [cm^-3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: Hpe_arr   = 0.0d0 ! photoelec. heat [erg cm^-3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: Hdx_arr   = 0.0d0 ! de-excit. heat [erg cm^-3 s^-1]
      ! Ground-state H proton-budget rate coefficients, captured by ioniz_eq on
      ! the converged pass so write_excited_H can compare the n=2 photoionization
      ! proton source against the ground-state channels (Huang Figs. 11/27):
      real*8, dimension(1-Ng:N+Ng) :: gph_ground_HI = 0.0d0 ! ground-state HI photoion. [s^-1]
      real*8, dimension(1-Ng:N+Ng) :: cion_HI       = 0.0d0 ! HI collisional ioniz. [cm^3 s^-1]
      real*8, dimension(1-Ng:N+Ng) :: arec_HII      = 0.0d0 ! HII recombination [cm^3 s^-1]

      ! NL solver vectors
      real*8, dimension(:), allocatable :: sys_sol, sys_x
      real*8, dimension(:), allocatable :: wa
      
      contains

      ! End of module      
      end module global_parameters
