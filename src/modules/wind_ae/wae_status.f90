      module wae_status
      ! Shared solve-status flag so the relaxation stack can FAIL GRACEFULLY
      ! (return a code) instead of stopping the program -- required by the
      ! continuation driver, which recovers from a failed step by halving
      ! the ramp delta. 0 = ok; nonzero = a failure mode:
      !   1 pinvs singular matrix
      !   2 get_dvdr erroneous weight
      !   3 linearize_dvdr_crit negative sqrt argument
      !   4 solvde too many iterations
      implicit none
      integer :: wae_err = 0
      end module wae_status
