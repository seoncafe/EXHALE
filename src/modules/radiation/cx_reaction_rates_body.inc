!  THE VOLUMETRIC RATE OF EVERY ACTIVE CHARGE-EXCHANGE REACTION, ONE SOURCE TEXT
!  FOR TWO ARITHMETICS: included by cx_reaction_rates_real and
!  cx_reaction_rates_dual of charge_exchange.
      call cx_require_cell(T_cell, caller)
      do i = 1, cx_nact
         r = cx_act(i)
         rrate(i) = cx_kc(i)*dens(cx_don_el(r), cx_don_stg(r))            &
                            *dens(cx_acc_el(r), cx_acc_stg(r))
      enddo
