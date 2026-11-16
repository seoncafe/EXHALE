!  THE NEUTRAL HELIUM OF A TRIAL STAGE PARTITION, ONE SOURCE TEXT FOR TWO
!  ARITHMETICS: included by carrier_helium_neutral_of_partition_real and
!  carrier_helium_neutral_of_partition_dual of diffusive_photochemistry
!  (the closure is stated above the first).
      n_heiSI = carrier_helium_available_to_stages(j) - n_heii - n_heiii
      if (carrier_solved(ic_HeTR)) n_heiSI = n_heiSI - n_heiTR
      n_hei   = n_heiSI + n_heiTR
