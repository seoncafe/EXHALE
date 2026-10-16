      module wae_grid
      ! The independent-variable grid x(1:M) (= q), shared between the
      ! relax driver and difeq exactly as wind-ae's global x[] is. The
      ! driver fills wae_x before calling solvde; difeq reads it.
      use wae_config, only: wae_m
      implicit none
      real*8 :: wae_x(wae_m) = 0.0d0
      end module wae_grid
