      module wae_lc_ciii
      ! Data table CIII(2,4) ported verbatim from wind-ae.
      implicit none
      integer, parameter :: CIII_nr = 2
      integer, parameter :: CIII_nc = 4
      real*8 :: CIII(2, 4) = 0.0d0
      contains
      subroutine CIII_init()
      CIII(1, :) = [ 1.9100000000d+03, -3.8422302400d-10, 7.5460800000d+04, 1.3147895300d+09 ]
      CIII(2, :) = [ 9.7700000000d+02, -1.7905083400d-03, 1.4726390000d+05, 7.1716480000d+14 ]
      end subroutine CIII_init
      end module wae_lc_ciii
