      module wae_lc_oiii
      ! Data table OIII(4,4) ported verbatim from wind-ae.
      implicit none
      integer, parameter :: OIII_nr = 4
      integer, parameter :: OIII_nc = 4
      real*8 :: OIII(4, 4) = 0.0d0
      contains
      subroutine OIII_init()
      OIII(1, :) = [ 5.2000000000d+05, -3.1385200000d-18, 2.7768200000d+02, 2.5493000000d+03 ]
      OIII(2, :) = [ 5.0000000000d+03, -3.3866780000d-14, 2.8728600000d+04, 9.6674100000d+05 ]
      OIII(3, :) = [ 1.6600000000d+02, -6.5997900000d-10, 8.6632400000d+04, 1.4758890000d+10 ]
      OIII(4, :) = [ 8.3500000000d+01, -1.7520500000d+03, 1.7256970000d+05, 5.4059370000d+21 ]
      end subroutine OIII_init
      end module wae_lc_oiii
