      module wae_lc_oii
      ! Data table OII(4,4) ported verbatim from wind-ae.
      implicit none
      integer, parameter :: OII_nr = 4
      integer, parameter :: OII_nc = 4
      real*8 :: OII(4, 4) = 0.0d0
      contains
      subroutine OII_init()
      OII(1, :) = [ 8.3400000000d+02, -5.7862200000d-04, 1.7242160000d+05, 1.3217400000d+15 ]
      OII(2, :) = [ 2.7410000000d+03, -3.8119800000d-13, 5.8225300000d+04, 4.4877700000d+07 ]
      OII(3, :) = [ 3.7270000000d+03, -4.2990100000d-16, 3.8575000000d+04, 5.3646100000d+03 ]
      OII(4, :) = [ 7.3200000000d+03, -3.7692900000d-13, 5.3063600000d+04, 3.1101800000d+07 ]
      end subroutine OII_init
      end module wae_lc_oii
