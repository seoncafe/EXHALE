      module wae_exhale_input
      ! Front-end that maps an EXHALE input.inp to the Wind-AE planet
      ! parameters (port plan §3). Fills the planet/star/flux/composition
      ! fields of wae_par; the base BCs (Rmin, rho_rmin, Ncol_sp, ...) are
      ! NOT set here -- those are inherited from the continuation seed and
      ! ramped, exactly as the Python wrapper does.
      use wae_config, only: nsp => wae_nspecies, T0 => wae_T0
      use wae_types,  only: wae_paramlist
      use wae_params, only: wae_par
      implicit none

      ! astronomical constants (cgs)
      real*8, parameter :: MJ   = 1.8982d30
      real*8, parameter :: RJ   = 7.1492d9
      real*8, parameter :: MSUN = 1.98842d33
      real*8, parameter :: RSUN = 6.957d10
      real*8, parameter :: AU   = 1.49598d13
      real*8, parameter :: SIGSB = 5.6705d-5
      real*8, parameter :: PIc  = 3.141592653589793d0

      ! values read from EXHALE input.inp the caller may also want
      real*8  :: exh_HeH = 0.0d0, exh_LEUV = 0.0d0, exh_LX = 0.0d0
      real*8  :: exh_Teff = 0.0d0, exh_Rstar = 0.0d0, exh_Teq = 0.0d0
      real*8  :: exh_lognbase = 14.0d0
      logical :: exh_spherical = .false.
      real*8  :: exh_rout = 0.0d0
      contains

      subroutine wae_read_exhale_input(fname)
      character(len=*), intent(in) :: fname
      character(len=1024) :: line, val
      integer :: u, ios
      real*8  :: MpMJ, RpRJ, MsMsun, aAU
      MpMJ=0; RpRJ=0; MsMsun=0; aAU=0
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0) cycle
         if      (has(line,'lower boundary number density')) then; call aft(line,val); read(val,*) exh_lognbase
         else if (has(line,'Planet radius'))  then; call aft(line,val); read(val,*) RpRJ
         else if (has(line,'Planet mass'))    then; call aft(line,val); read(val,*) MpMJ
         else if (has(line,'Equilibrium temperature')) then; call aft(line,val); read(val,*) exh_Teq
         else if (has(line,'Orbital distance')) then; call aft(line,val); read(val,*) aAU
         else if (has(line,'He/H number ratio')) then; call aft(line,val); read(val,*) exh_HeH
         else if (has(line,'Parent star mass')) then; call aft(line,val); read(val,*) MsMsun
         else if (has(line,'X-ray luminosity')) then; call aft(line,val); read(val,*) exh_LX
         else if (has(line,'EUV luminosity'))   then; call aft(line,val); read(val,*) exh_LEUV
         else if (has(line,'Stellar Teff'))     then; call aft(line,val); read(val,*) exh_Teff
         else if (has(line,'Stellar radius'))   then; call aft(line,val); read(val,*) exh_Rstar
         else if (has(line,'Domain mode'))      then
            call aft(line,val)
            if (index(val,'Spherical') .gt. 0) exh_spherical = .true.
         else if (has(line,'Outer radius'))     then; call aft(line,val); read(val,*) exh_rout
         end if
      end do
      close(u)

      ! --- map to wind-ae PARAMLIST ---
      wae_par%Mp        = MpMJ*MJ
      wae_par%Rp        = RpRJ*RJ
      wae_par%Mstar     = MsMsun*MSUN
      wae_par%semimajor = aAU*AU
      wae_par%Ftot = (10.0d0**exh_LEUV + 10.0d0**exh_LX) /              &
                     (4.0d0*PIc*wae_par%semimajor**2)
      if (exh_Teff .gt. 0.0d0 .and. exh_Rstar .gt. 0.0d0) then
         wae_par%Lstar = 4.0d0*PIc*(exh_Rstar*RSUN)**2*SIGSB*exh_Teff**4
      else
         wae_par%Lstar = 0.0d0   ! bolo terms negligible in the wind
      end if
      ! composition: He/H number ratio -> H/He mass fractions
      wae_par%HX(1) = 1.0d0/(1.0d0 + 4.0d0*exh_HeH)
      wae_par%HX(2) = 4.0d0*exh_HeH/(1.0d0 + 4.0d0*exh_HeH)
      ! CODATA H-atom and He-atom masses; wind-ae originally passed
      ! 1.6733d-24 and 6.6464790722d-24 here.
      wae_par%atomic_mass(1) = 1.67353284d-24
      wae_par%atomic_mass(2) = 6.6464790722d-24
      wae_par%species(1) = 'HI'
      wae_par%species(2) = 'HeI'
      wae_par%Z(1) = 1; wae_par%Z(2) = 2
      wae_par%N_e(1) = 1; wae_par%N_e(2) = 2
      wae_par%T_rmin = exh_Teq/T0
      wae_par%tidalforce = 1.0d0
      if (exh_spherical) wae_par%tidalforce = 0.0d0
      end subroutine wae_read_exhale_input

      ! ---- helpers ----
      logical function has(line, key)
      character(len=*), intent(in) :: line, key
      has = index(line, key) .gt. 0
      end function has

      subroutine aft(line, val)
      ! text after ':' with trailing '#...' comment stripped
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: val
      integer :: a, h
      a = index(line, ':') + 1
      val = adjustl(line(a:))
      h = index(val, '#')
      if (h .gt. 0) val = val(1:h-1)
      end subroutine aft

      end module wae_exhale_input
