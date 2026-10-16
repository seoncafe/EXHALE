      module wae_lc_cii
      ! Data table CII(3,4) ported verbatim from wind-ae.
      implicit none
      integer, parameter :: CII_nr = 3
      integer, parameter :: CII_nc = 4
      real*8 :: CII(3, 4) = 0.0d0
      contains
      subroutine CII_init()
      CII(1, :) = [ 1.5700000000d+06, -1.7829205400d-20, 9.1200000000d+01, 1.3877966800d+01 ]
      CII(2, :) = [ 2.3260000000d+03, -1.2147102300d-10, 6.1853900000d+04, 1.2097763300d+09 ]
      CII(3, :) = [ 1.3340000000d+03, -2.4130450800d-03, 1.0771810000d+05, 3.7400877000d+15 ]
      end subroutine CII_init
      end module wae_lc_cii
