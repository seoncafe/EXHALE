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
                                    spherical_domain, r, N, Ng
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
                                 cy, rhoscale
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

      seedf = 'inputdata/windae_seed.csv'
      specf = 'inputdata/windae_spectrum.inp'

      ! data tables + spectrum + seed (sets wae_par to the seed)
      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call load_seed(trim(seedf))

      ! EXHALE Domain mode -> tidalforce (discrete, not ramped)
      par%tidalforce = 1.0d0
      if (spherical_domain) par%tidalforce = 0.0d0
      ! convergence scale: strongly-bound planets need RHOSCALE=10
      rhoscale = 100.0d0
      if (b0 .gt. 120.0d0) rhoscale = 10.0d0
      call setup_indices_scales()

      ! targets from EXHALE globals
      Ftot_t = J_XUV                                   ! = (10^LX+10^LEUV)/(4 pi a^2)
      if (T_star_eff .gt. 0.0d0 .and. R_star .gt. 0.0d0) then
         Lstar_t = 4.0d0*WAE_PI*R_star**2*WAE_SIGSB*T_star_eff**4
      else
         Lstar_t = 0.0d0
      end if

      write(*,'(A)') '      (wae bridge) ramping seed -> this planet...'
      rc = ramp_to(Ftot_t, Mp, R0, Mstar, a_orb, Lstar_t)
      if (rc .ne. 0) then
         write(*,*) '      (wae bridge) RAMP FAILED (code', rc, ').'
         write(*,*) '      This planet may be far from the shipped seed; '// &
                    'try "IC mode: auto" instead.'
         stop '(wae_generate_ic) ramp failed'
      end if

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
                        log10(n0), T0, R0, HeH, 'output')
      write(*,'(A)') '      (wae bridge) wrote output/{Hydro_ioniz,'//    &
                     'Ion_species}_IC.txt'
      end subroutine wae_generate_ic

      end module wae_exhale_bridge
