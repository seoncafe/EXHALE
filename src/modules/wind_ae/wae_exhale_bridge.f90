      module wae_exhale_bridge
      ! In-process bridge: build a Wind-AE warm-start IC directly from
      ! EXHALE's already-parsed parameters and grid, writing
      ! output/{Hydro_ioniz,Ion_species}_IC.txt for load_IC. Invoked by
      ! init.f90 when "IC mode: windae" (ic_mode=4). This module is the
      ! ONLY place the Wind-AE tree touches EXHALE's global_parameters, so
      ! the standalone wind_ae_ic.x build (which excludes this file) stays
      ! self-contained.
      use global_parameters, only: Mp, R0, Mstar, a_orb, T0, n0, HeH,    &
                                    J_XUV, b0, T_star_eff, R_star,        &
                                    spherical_domain, r, N, Ng,           &
                                    windae_seed_file, windae_seed_out
      use wae_params,      only: par => wae_par
      use wae_spectrum,    only: wae_load_spectrum
      use wae_rate_coeffs, only: wae_rate_coeffs_init
      use wae_rnew,        only: rnew_init
      use wae_rrec,        only: rrec_init
      use wae_fe,          only: fe_init
      use wae_ion_pots,    only: ion_pots_init
      use wae_lc_cii,      only: cii_init
      use wae_lc_ciii,     only: ciii_init
      use wae_lc_oii,      only: oii_init
      use wae_lc_oiii,     only: oiii_init
      use wae_config,      only: wae_nsp => wae_nspecies, wae_M => wae_m, &
                                 wae_addpts, RHO0 => wae_RHO0, WT0 => wae_T0
      use wae_params,      only: wae_cs0_val
      use wae_grid,        only: wae_x
      use wae_continuation,only: load_seed, setup_indices_scales, ramp_to, &
                                 ramp_tidal, cy, rhoscale, dump_seed,      &
                                 pick_nearest_seed
      use wae_intode,      only: wae_integrate_ode
      use wae_ic_writer,   only: wae_write_ic
      implicit none
      real*8, parameter :: WAE_PI = 3.141592653589793d0
      real*8, parameter :: WAE_SIGSB = 5.6705d-5
      contains

      subroutine wae_generate_ic()
      character(len=256) :: seedf, specf
      real*8  :: Ftot_t, Lstar_t, CS0
      integer :: i, j, row, ntot, tp, rc, ng_grid
      real*8, allocatable :: sr(:), srho(:), sv(:), sT(:), sq(:), sz(:)
      real*8, allocatable :: sYs(:,:), sNcol(:,:)
      real*8, allocatable :: dr(:), drho(:), dv(:), dT(:), dHI(:), dHeI(:)
      real*8, allocatable :: rgrid(:)

      ! Resolve the warm-start seed: an explicit CSV path from input.inp
      ! ("Wind-AE seed:"), or the keyword 'grid'/'auto' to auto-pick the
      ! nearest solution in inputdata/windae_grid/ to this planet's (Mp,Rp,flux)
      ! -- a far shorter, more robust ramp than always starting from the single
      ! fiducial HJ seed.
      if (trim(windae_seed_file) .eq. 'grid' .or.                          &
          trim(windae_seed_file) .eq. 'auto') then
         call pick_nearest_seed(Mp, R0, J_XUV, 'inputdata/windae_grid', seedf)
      else
         seedf = trim(windae_seed_file)
      end if
      specf = 'inputdata/windae_spectrum.inp'

      ! data tables + spectrum + seed (sets wae_par to the seed)
      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call load_seed(trim(seedf))
      ! rho convergence scale from the loaded seed's base density, matching the
      ! reference C code's 10^floor(log10(rho_rmin*0.01)); the self-consistent-BC ramp updates
      ! it as rho_rmin changes. (The old "b0>120 -> 10" heuristic mis-scaled the
      ! seed -- whose rho_rmin ~ 2e4 needs ~100 -- and stalled the first solve.)
      rhoscale = 10.0d0**floor(log10(par%rho_rmin*0.01d0))
      call setup_indices_scales()
      ! EXHALE Domain mode -> tidalforce. The seed is converged WITH tidal
      ! gravity (tidalforce=1); for a spherical domain, ramp it smoothly to 0
      ! (a discrete 1->0 jump would disrupt the seed and stall the first solve).
      if (spherical_domain) rc = ramp_tidal(0.0d0)

      ! targets from EXHALE globals
      Ftot_t = J_XUV                                   ! = (10^LX+10^LEUV)/(4 pi a^2)
      if (T_star_eff .gt. 0.0d0 .and. R_star .gt. 0.0d0) then
         Lstar_t = 4.0d0*WAE_PI*R_star**2*WAE_SIGSB*T_star_eff**4
      else
         Lstar_t = 0.0d0
      end if

      ! Static-BC ramp first: fast for seed-adjacent hot Jupiters.
      write(*,'(A)') '      (wae bridge) ramping seed -> this planet (static-BC ramp)...'
      rc = ramp_to(Ftot_t, Mp, R0, Mstar, a_orb, Lstar_t)
      if (rc .ne. 0) then
         ! static-BC ramp stalled -> strongly-bound / far-from-seed planet (e.g. HD189733b).
         ! Reload the seed and retry with the self-consistent-BC ramp: re-converge the base BCs and
         ! turn the molecular layer off once the base falls inside the wind.
         write(*,'(A)') '      (wae bridge) static-BC ramp stalled; retrying with'// &
                        ' self-consistent-BC ramp...'
         call load_seed(trim(seedf))
         rhoscale = 10.0d0**floor(log10(par%rho_rmin*0.01d0))
         call setup_indices_scales()
         if (spherical_domain) rc = ramp_tidal(0.0d0)
         rc = ramp_to(Ftot_t, Mp, R0, Mstar, a_orb, Lstar_t, static_bcs=.false.)
      end if
      if (rc .ne. 0) then
         write(*,*) '      (wae bridge) RAMP FAILED (code', rc, ').'
         write(*,*) '      This planet may be unreachable from the shipped seed; '// &
                    'try "IC mode: auto" instead.'
         stop '(wae_generate_ic) ramp failed'
      end if

      ! optionally bank the converged relaxation solution as a reusable seed
      ! ("Wind-AE seed out:"), e.g. the first spherical (tidalforce=0) solution
      ! or one close to a hard target -- growing inputdata/windae_grid/.
      if (len_trim(windae_seed_out) .gt. 0) call dump_seed(trim(windae_seed_out))

      ! relax solution is in cy; integrate outward
      tp = wae_M + wae_addpts
      allocate(sr(tp), srho(tp), sv(tp), sT(tp), sq(tp), sz(tp))
      allocate(sYs(tp,wae_nsp), sNcol(tp,wae_nsp))
      do row = 1, wae_M
         sq(row)=wae_x(row); sz(row)=cy(2,row)
         sr(row)=wae_x(row)*cy(2,row)+par%Rmin
         srho(row)=cy(3,row); sv(row)=cy(1,row); sT(row)=cy(4,row)
         do j = 1, wae_nsp
            sYs(row,j)=cy(4+j,row); sNcol(row,j)=cy(4+wae_nsp+j,row)
         end do
      end do
      call wae_integrate_ode(sr, srho, sv, sT, sYs, sNcol, sq, sz, ntot)

      ! dimensionalize (CODE -> cgs) for the IC writer
      CS0 = wae_cs0_val
      allocate(dr(ntot), drho(ntot), dv(ntot), dT(ntot), dHI(ntot), dHeI(ntot))
      do row = 1, ntot
         dr(row)   = sr(row)*R0      ! cm  (R0 = wind-ae Rp = EXHALE radius)
         drho(row) = srho(row)*RHO0
         dv(row)   = sv(row)*CS0
         dT(row)   = sT(row)*WT0
         dHI(row)  = sYs(row,1)
         dHeI(row) = sYs(row,2)
      end do

      ! EXHALE grid (incl. ghosts), in R0 units, ascending
      ng_grid = N + 2*Ng
      allocate(rgrid(ng_grid))
      do i = 1, ng_grid
         rgrid(i) = r(i - Ng)        ! r is indexed 1-Ng : N+Ng
      end do

      ! write the IC onto the EXHALE grid (base anchor = EXHALE n0/T0/Rp/HeH)
      call wae_write_ic(dr, drho, dv, dT, dHI, dHeI, ntot, rgrid, ng_grid, &
                        log10(n0), T0, R0, HeH, par%Mstar/par%Mp,          &
                        par%semimajor/par%Rp, par%tidalforce, 'output')
      write(*,'(A)') '      (wae bridge) wrote output/{Hydro_ioniz,'//    &
                     'Ion_species}_IC.txt'
      end subroutine wae_generate_ic

      end module wae_exhale_bridge
