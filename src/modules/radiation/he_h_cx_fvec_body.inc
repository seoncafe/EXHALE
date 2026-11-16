!  THE He <-> H CHARGE-EXCHANGE PAIR AND He2+ + H CAPTURE ON THE ROWS, ONE SOURCE
!  TEXT FOR TWO ARITHMETICS: included by he_h_cx_fvec_real (real*8) and
!  he_h_cx_fvec_dual (dual densities and rows) of charge_exchange.  The
!  orientation rules are stated above he_h_cx_fvec_real.
      if (.not. he_h_charge_exchange) return
      R1 = k1*n_hei*n_hii       ! He0 + H+ -> He+ + H0
      R2 = k2*n_heii*n_hi       ! He+ + H0 -> He0 + H+
      R3 = k3*n_heiii*n_hi      ! He2+ + H0 -> He+ + H+ + photon
      ! H row (HI->HII positive): B2 and B3 ionize H0 (+R2, +R3), B1
      ! recombines H+ (-R1).
      fvec(1) = fvec(1) + R2 - R1 + R3
      ! He row: B1 ionizes HeI->HeII (+R1), B2 recombines HeII->HeI (-R2), in
      ! the ionization-positive orientation; he_row_sign flips it for the
      ! HeI-gain (TR) orientation.
      fvec(2) = fvec(2) + he_row_sign*(R1 - R2)
      if (heii_row_is_stage_source) fvec(2) = fvec(2) + R3
      fvec(3) = fvec(3) - R3
