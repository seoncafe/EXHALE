!  THE CHARGE-EXCHANGE SET AS STOICHIOMETRIC STAGE SOURCES, ONE SOURCE TEXT FOR
!  TWO ARITHMETICS: included by charge_exchange_stage_sources_real (real*8
!  densities, with the gross production and loss) and
!  charge_exchange_stage_sources_dual (dual densities, value and first
!  derivatives) of charge_exchange.  The rule is stated above
!  charge_exchange_stage_sources_real.
      call cx_reactant_densities(nm0, nm1, nm2, n_hi, n_hii, n_hei,       &
                                 n_heii, n_heiii, dens)
      call cx_reaction_rates(dens, T_cell,                                &
                             'charge_exchange_stage_sources', rrate)

      ssrc  = 0.0d0
      do i = 1, cx_nact
         r  = cx_act(i)
         de = cx_don_el(r);  ds = cx_don_stg(r)
         ae = cx_acc_el(r);  as = cx_acc_stg(r)
         ssrc(de,ds)    = ssrc(de,ds)    - rrate(i)
         ssrc(de,ds+1)  = ssrc(de,ds+1)  + rrate(i)
         ssrc(ae,as)    = ssrc(ae,as)    - rrate(i)
         ssrc(ae,as-1)  = ssrc(ae,as-1)  + rrate(i)
      enddo
      if (maxval(abs(value_of(ssrc(:,3)))) .gt. 0.0d0) then
         write(*,*) '(charge_exchange) ERROR: a reaction of the active',  &
                    ' set promotes a donor above ionization stage 2,'
         write(*,*) '  which the three-stage element layout of this',     &
                    ' module and of the metal systems cannot hold.'
         write(*,*) '  Aborting.'
         error stop 1
      endif

      s_H(0:1)  = ssrc(cx_H,0:1)
      s_He(0:2) = ssrc(cx_He,0:2)
