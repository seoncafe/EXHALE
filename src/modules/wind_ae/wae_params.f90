      module wae_params
      ! The parameter list of the Wind-AE port (the PARAMLIST of wind-ae's
      ! io.c, with the C field meanings) and the physical constants it is
      ! built on. wae_exhale_input fills it from an EXHALE input.inp.
      use wae_types,  only: wae_paramlist
      implicit none

      ! Physical constants. Wind-AE's defs.h carries older values; EXHALE
      ! uses the current CODATA ones so that the IC it hands over is built
      ! on the same constants as the rest of the code (2026-08-19).
      ! Original defs.h values, kept here for comparison against the C/Python
      ! oracle: G = 6.67259d-8, K = 1.380658d-16, MH = 1.6733d-24.
      real*8, parameter :: wae_G  = 6.67430d-8      ! CODATA 2018
      real*8, parameter :: wae_K  = 1.380649d-16    ! CODATA exact
      real*8, parameter :: wae_MH = 1.67353284d-24  ! H atom, CODATA 2018
      real*8, parameter :: wae_CS0 = 1.0d0   ! set in init (sqrt(K*T0/MH))

      type(wae_paramlist) :: wae_par
      real*8 :: wae_cs0_val = 0.0d0

      end module wae_params
