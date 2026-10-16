      program eqn_trace_f
      ! Fortran counterpart of wind-ae src/eqn_trace.c: the residual-layer
      ! gate. Loads a windsoln (first M rows) into y(ne,m)/x(m), calls
      ! wae_linearize_dvdr_crit, then dumps wae_eval_eqn at the same k and
      ! eqn numbers with ymod=0.
      use wae_config,     only: nsp => wae_nspecies, m => wae_m, ne => wae_ne
      use wae_spectrum,   only: wae_load_spectrum, wae_load_params
      use wae_params,     only: wae_load_all_params
      use wae_rate_coeffs,only: wae_rate_coeffs_init
      use wae_rnew,       only: rnew_init
      use wae_rrec,       only: rrec_init
      use wae_fe,         only: fe_init
      use wae_ion_pots,   only: ion_pots_init
      use wae_lc_cii,     only: cii_init
      use wae_lc_ciii,    only: ciii_init
      use wae_lc_oii,     only: oii_init
      use wae_lc_oiii,    only: oiii_init
      use wae_soe,        only: wae_linearize_dvdr_crit
      use wae_eqns,       only: wae_eval_eqn, VEQN, RHOEQN, TEQN, IONEQN, &
                                NCOLEQN, SPVEQN, SPCRITEQN
      implicit none

      character(len=4096) :: wsoln, outcsv, specf, indir, line
      integer :: u, uo, ios, row, t, k
      real*8  :: x(m), y(ne, m), ymod(0:ne, 0:1)
      real*8  :: cval(10)
      integer :: ks(10)
      real*8  :: rv, rr, rt, ri1, ri2, rn1, rn2, spv, spc

      if (command_argument_count() .lt. 4) then
         write(*,*) 'usage: eqn_trace_f <windsoln> <out.csv> <spectrum.inp> <inputs_dir/>'
         stop 2
      end if
      call get_command_argument(1, wsoln)
      call get_command_argument(2, outcsv)
      call get_command_argument(3, specf)
      call get_command_argument(4, indir)

      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call wae_load_params(trim(indir)//'guess.inp')
      call wae_load_all_params(trim(indir))

      ! load windsoln first M rows into x,y
      open(newunit=u, file=trim(wsoln), status='old', action='read')
      row = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         if (row .ge. m) exit
         call parse_csv(line, cval, 10)
         row = row + 1
         ! windsoln cols: r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI,q,z
         x(row)     = cval(9)     ! q
         y(1,row)   = cval(3)     ! v
         y(2,row)   = cval(10)    ! z
         y(3,row)   = cval(2)     ! rho
         y(4,row)   = cval(4)     ! T
         y(5,row)   = cval(5)     ! Ys_HI
         y(6,row)   = cval(6)     ! Ys_HeI
         y(7,row)   = cval(7)     ! Ncol_HI
         y(8,row)   = cval(8)     ! Ncol_HeI
      end do
      close(u)
      if (row .lt. m) then
         write(*,*) 'only', row, 'of', m, 'rows'; stop 4
      end if

      call wae_linearize_dvdr_crit(x, y, m, ne)
      ymod = 0.0d0

      ks = [2, 10, 50, 100, 300, 700, 1100, 1400, 1500, 1501]
      open(newunit=uo, file=trim(outcsv), status='replace', action='write')
      do t = 1, 10
         k = ks(t)
         rv  = wae_eval_eqn(k, x, y, ymod, VEQN, 1, m, ne)
         rr  = wae_eval_eqn(k, x, y, ymod, RHOEQN, 1, m, ne)
         rt  = wae_eval_eqn(k, x, y, ymod, TEQN, 1, m, ne)
         ri1 = wae_eval_eqn(k, x, y, ymod, IONEQN, 1, m, ne)
         ri2 = wae_eval_eqn(k, x, y, ymod, IONEQN, 2, m, ne)
         rn1 = wae_eval_eqn(k, x, y, ymod, NCOLEQN, 1, m, ne)
         rn2 = wae_eval_eqn(k, x, y, ymod, NCOLEQN, 2, m, ne)
         write(uo,'(I0,7(A,ES24.17))') k, ',',rv, ',',rr, ',',rt,        &
               ',',ri1, ',',ri2, ',',rn1, ',',rn2
      end do
      spv = wae_eval_eqn(m, x, y, ymod, SPVEQN, 1, m, ne)
      spc = wae_eval_eqn(m, x, y, ymod, SPCRITEQN, 1, m, ne)
      write(uo,'(I0,2(A,ES24.17))') -1, ',',spv, ',',spc
      close(uo)

      contains
      subroutine parse_csv(s, vals, n)
      character(len=*), intent(in) :: s
      integer, intent(in) :: n
      real*8, intent(out) :: vals(n)
      integer :: p, j, b, L
      character(len=64) :: tok
      p = 1; L = len_trim(s)
      do j = 1, n
         b = index(s(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(s(p:L)); p = L + 1
         else
            tok = adjustl(s(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(j)
      end do
      end subroutine parse_csv
      end program eqn_trace_f
