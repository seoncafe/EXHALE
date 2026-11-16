!  THE REACTANT DENSITY TABLE OF THE CHARGE-EXCHANGE SET, ONE SOURCE TEXT FOR
!  TWO ARITHMETICS: included by cx_reactant_densities_real and
!  cx_reactant_densities_dual of charge_exchange.
      dens = 0.0d0
      do i = 1, n_melem
         dens(i,0) = nm0(i)
         dens(i,1) = nm1(i)
         dens(i,2) = nm2(i)
      enddo
      dens(cx_H,0)  = n_hi
      dens(cx_H,1)  = n_hii
      dens(cx_He,0) = n_hei
      dens(cx_He,1) = n_heii
      dens(cx_He,2) = n_heiii
