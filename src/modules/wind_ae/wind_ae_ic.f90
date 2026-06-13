      program wind_ae_ic
      ! P2 deliverable: read an EXHALE input.inp, solve the Wind-AE wind
      ! (ramp from a shipped seed with static BCs, then relax + outward
      ! integrate), and write EXHALE Load-IC files. Metal-free H/He
      ! warm-start IC generator -- no Python in the loop.
      !
      ! Usage: wind_ae_ic <EXHALE_input.inp> <seed.csv> <spectrum.inp> \
      !                   <IC_dump_grid> <outdir> [rhoscale]
      use wae_config,      only: nsp => wae_nspecies, m => wae_m,         &
                                 ne => wae_ne, addpts => wae_addpts,      &
                                 RHO0 => wae_RHO0, T0 => wae_T0
      use wae_params,      only: par => wae_par, wae_cs0_val
      use wae_spectrum,    only: wae_load_spectrum
      use wae_exhale_input,only: wae_read_exhale_input, exh_HeH, exh_Teq, &
                                 exh_lognbase, exh_spherical
      use wae_rate_coeffs, only: wae_rate_coeffs_init
      use wae_rnew,        only: rnew_init
      use wae_rrec,        only: rrec_init
      use wae_fe,          only: fe_init
      use wae_ion_pots,    only: ion_pots_init
      use wae_lc_cii,      only: cii_init
      use wae_lc_ciii,     only: ciii_init
      use wae_lc_oii,      only: oii_init
      use wae_lc_oiii,     only: oiii_init
      use wae_grid,        only: wae_x
      use wae_intode,      only: wae_integrate_ode
      use wae_continuation,only: load_seed, setup_indices_scales, ramp_to, &
                                 ramp_tidal, cy, rhoscale, dump_seed
      use wae_ic_writer,   only: wae_write_ic
      implicit none
      character(len=4096) :: exin, seed, specf, gridf, outdir, rhoarg, line, dumpf
      real*8  :: Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t, Rp_cgs, CS0
      integer :: i, j, row, ntot, tp, r, u, ios, ng
      real*8, allocatable :: sr(:), srho(:), sv(:), sT(:), sq(:), sz(:)
      real*8, allocatable :: sYs(:,:), sNcol(:,:), rgrid(:)
      real*8, allocatable :: dr(:), drho(:), dv(:), dT(:), dHI(:), dHeI(:)

      if (command_argument_count() .lt. 5) then
         write(*,*) 'usage: wind_ae_ic <EXHALE_input.inp> <seed.csv> '//  &
                    '<spectrum.inp> <IC_dump_grid> <outdir> '//           &
                    '[rhoscale] [seed_out.csv]'
         stop 2
      end if
      call get_command_argument(1, exin)
      call get_command_argument(2, seed)
      call get_command_argument(3, specf)
      call get_command_argument(4, gridf)
      call get_command_argument(5, outdir)
      ! optional arg 7: dump the converged solution as a reusable seed CSV
      dumpf = ''
      if (command_argument_count() .ge. 7) call get_command_argument(7, dumpf)

      ! 1. EXHALE input -> target params (capture before load_seed overwrites)
      call wae_read_exhale_input(trim(exin))
      Ftot_t=par%Ftot; Mp_t=par%Mp; Rp_t=par%Rp
      Mstar_t=par%Mstar; a_t=par%semimajor; Lstar_t=par%Lstar

      ! 2. spectrum + data tables
      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))

      ! 3. seed (sets wae_par = seed params/BCs/composition + cy)
      call load_seed(trim(seed))
      rhoscale = 10.0d0**floor(log10(par%rho_rmin*0.01d0))  ! from seed base rho
      if (command_argument_count() .ge. 6) then
         call get_command_argument(6, rhoarg); read(rhoarg,*) rhoscale
      end if
      call setup_indices_scales()
      ! tidalforce: the seed has it on (=1); for a spherical EXHALE domain ramp
      ! it smoothly to 0 (a discrete jump would disrupt the tidally-converged
      ! seed and stall the first relaxation).
      if (exh_spherical) r = ramp_tidal(0.0d0)

      ! 4. ramp to EXHALE target params. Stage C-1 first (static base BCs):
      !    fast for seed-adjacent planets. If it stalls (far-from-seed /
      !    strongly-bound), reload the seed and retry with stage C-2.
      write(*,'(A)') ' (wind_ae_ic) ramping seed -> EXHALE planet (stage C-1)...'
      r = ramp_to(Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t)
      if (r .ne. 0) then
         write(*,'(A)') ' (wind_ae_ic) C-1 stalled; retrying with stage C-2...'
         call load_seed(trim(seed))
         rhoscale = 10.0d0**floor(log10(par%rho_rmin*0.01d0))
         call setup_indices_scales()
         if (exh_spherical) r = ramp_tidal(0.0d0)
         r = ramp_to(Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t, static_bcs=.false.)
      end if
      if (r .ne. 0) then
         write(*,*) ' RAMP FAILED (code', r, '). Planet may be unreachable '// &
                    'from the shipped seed.'; stop 1
      end if

      ! optionally bank the converged relaxation solution as a reusable seed
      if (len_trim(dumpf) .gt. 0) call dump_seed(trim(dumpf))

      ! 5. relax solution is in cy; integrate outward
      tp = m + addpts
      allocate(sr(tp), srho(tp), sv(tp), sT(tp), sq(tp), sz(tp))
      allocate(sYs(tp,nsp), sNcol(tp,nsp))
      do row = 1, m
         sq(row)=wae_x(row); sz(row)=cy(2,row)
         sr(row)=wae_x(row)*cy(2,row)+par%Rmin
         srho(row)=cy(3,row); sv(row)=cy(1,row); sT(row)=cy(4,row)
         do j = 1, nsp
            sYs(row,j)=cy(4+j,row); sNcol(row,j)=cy(4+nsp+j,row)
         end do
      end do
      call wae_integrate_ode(sr, srho, sv, sT, sYs, sNcol, sq, sz, ntot)

      ! 6. dimensionalize and write EXHALE IC on the EXHALE grid
      Rp_cgs = par%Rp
      CS0 = wae_cs0_val
      allocate(dr(ntot), drho(ntot), dv(ntot), dT(ntot), dHI(ntot), dHeI(ntot))
      do row = 1, ntot
         dr(row)   = sr(row)*Rp_cgs        ! cm
         drho(row) = srho(row)*RHO0        ! g/cm^3
         dv(row)   = sv(row)*CS0           ! cm/s
         dT(row)   = sT(row)*T0            ! K
         dHI(row)  = sYs(row,1)            ! neutral fractions
         dHeI(row) = sYs(row,2)
      end do
      ! read EXHALE grid (r in col 1)
      ng = 0
      open(newunit=u, file=trim(gridf), status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .ne. '#' .and. len_trim(line) .gt. 0) ng = ng + 1
      end do
      allocate(rgrid(ng)); rewind(u); i = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         i = i + 1; read(line,*) rgrid(i)
      end do
      close(u)

      call wae_write_ic(dr, drho, dv, dT, dHI, dHeI, ntot, rgrid, ng,    &
                        exh_lognbase, exh_Teq, Rp_cgs, exh_HeH,          &
                        par%Mstar/par%Mp, par%semimajor/par%Rp,          &
                        par%tidalforce, trim(outdir))
      write(*,'(A,A)') ' (wind_ae_ic) wrote IC files to ', trim(outdir)
      end program wind_ae_ic
