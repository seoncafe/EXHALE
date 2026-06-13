      program ramp_test
      ! C-1 test: load a seed windsoln, ramp the system params to a target
      ! (static BCs), integrate outward, write the full windsoln.
      ! Usage: ramp_test <seed.csv> <spectrum.inp> <out.csv> <rhoscale> \
      !                  <Ftot> <Mp> <Rp> <Mstar> <a> <Lstar>
      use wae_config,      only: nsp => wae_nspecies, m => wae_m,         &
                                 ne => wae_ne, addpts => wae_addpts
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
      use wae_grid,        only: wae_x
      use wae_intode,      only: wae_integrate_ode
      use wae_continuation,only: load_seed, setup_indices_scales, ramp_to, &
                                 cy, rhoscale
      implicit none
      character(len=4096) :: seed, specf, outcsv, a(7)
      real*8  :: tgt(6)
      integer :: i, j, row, ntot, tp, uo, r
      real*8, allocatable :: sr(:), srho(:), sv(:), sT(:), sq(:), sz(:)
      real*8, allocatable :: sYs(:,:), sNcol(:,:)

      call get_command_argument(1, seed)
      call get_command_argument(2, specf)
      call get_command_argument(3, outcsv)
      do i = 1, 7
         call get_command_argument(3+i, a(i))
      end do
      read(a(1),*) rhoscale
      do i = 1, 6
         read(a(1+i),*) tgt(i)
      end do

      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call load_seed(trim(seed))
      call setup_indices_scales()

      write(*,'(A,6ES12.5)') ' seed params Ftot/Mp/Rp/Ms/a/Ls = ',       &
         par%Ftot, par%Mp, par%Rp, par%Mstar, par%semimajor, par%Lstar
      write(*,'(A,6ES12.5)') ' target                        = ', tgt
      r = ramp_to(tgt(1), tgt(2), tgt(3), tgt(4), tgt(5), tgt(6))
      if (r .ne. 0) then
         write(*,*) 'RAMP FAILED, code', r; stop 1
      end if
      write(*,*) 'ramp complete; integrating outward...'

      ! map relaxed cy -> full-grid storage
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

      open(newunit=uo, file=trim(outcsv), status='replace', action='write')
      write(uo,'(A)') '#vars: r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI,q,z'
      do row = 1, ntot
         write(uo,'(9(ES24.17,A),ES24.17)') sr(row),',',srho(row),',',   &
            sv(row),',',sT(row),',',sYs(row,1),',',sYs(row,2),',',       &
            sNcol(row,1),',',sNcol(row,2),',',sq(row),',',sz(row)
      end do
      close(uo)
      write(*,'(A,I0,A)') ' wrote ', ntot, ' rows.'
      end program ramp_test
