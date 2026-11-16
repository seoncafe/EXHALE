!  THE ASSEMBLY OF THE CARRIER SOURCES FROM THE EVALUATED ROWS, ONE SOURCE
!  TEXT FOR TWO ARITHMETICS: included by carrier_source and
!  carrier_source_dual of diffusive_photochemistry, after the balance rows
!  (fv), the metal charge exchange (cx_sH, cx_sHe) and the closure
!  (carrier_source_closure.inc) have been evaluated.

      ! THE ROWS OF THE CARRIERS, from the balances evaluated above.  In
      ! molecular gas the H2, molecular-ion and oxygen-carrier rows are the
      ! network's own rows (4)-(7), (9) and (10); an atomic gas has none of
      ! them and their sources stay at zero.
      if (thereis_mol) then
         src(ic_H2)  = fv(4)
         if (carrier_solved(ic_H2p)) src(ic_H2p) = fv(5)
         if (carrier_solved(ic_H3p)) src(ic_H3p) = fv(6)
         if (carrier_solved(ic_HeHp)) src(ic_HeHp) = fv(7)
         src(ic_OH)  = fv(9)
         src(ic_H2O) = fv(10)
      endif
      ! ---- CO: the one-sided destruction model -----------------------
      ! Two channels and no formation term, which is what "one-sided"
      ! means:
      !
      !   D1  He+ + CO -> C+ + O + He, UMIST RATE22 entry 4068, measured,
      !       accuracy better than 25 percent, temperature independent at
      !       the Langevin value 1.6e-9 cm^3 s^-1 (oxygen_rates::
      !       rk_D1_Hep_CO).  n_heii is the frozen He+ of this step.
      !   D2  CO + hv -> C + O on the 912-1201 A beam, with the Visser,
      !       van Dishoeck & Black (2009) shielding function of the
      !       star-ward CO and H2 columns.  cph_kco is that rate, the cell
      !       mean, rebuilt from the CURRENT carrier columns at the top of
      !       the step, so the CO layer shields itself inside the
      !       relaxation.
      !
      ! Both are proportional to the trial nc(ic_CO), so the term is a
      ! diagonal contribution to the CO row and the Jacobian sees it.
      !
      ! WHERE THE NUCLEI GO.  The oxygen is already accounted for: n_o0
      ! above is the closure nO_free - ions - OH - H2O - CO, so an oxygen
      ! nucleus that leaves CO appears as free atomic O in the same
      ! evaluation of the same cell and the oxygen element total does not
      ! move.  The carbon returns to the carbon pool and carrier_write_back
      ! shares it over C I / C II / C III in proportion to what they held;
      ! the ionization sweep that follows re-solves that partition from its
      ! own balance.  The step moves nuclei, the sweep moves charge.
      !
      ! WHERE THE He+ IS CONSUMED.  D1 destroys a helium ion as well as a CO
      ! molecule, and the two operators split the one reaction between them:
      ! this row removes the CO at the frozen He+ of the step, and row (2)
      ! of mol_heh_rows removes the He+ at the frozen CO of the sweep, both
      ! from rk_D1_Hep_CO, so the reaction is complete over the pair and
      ! counted once in each.  What the split leaves is a lag of one
      ! operator, not a missing sink; the domain record counts the cells in
      ! which D1 is the leading He+ loss (carrier_co_domain_record), which
      ! is where that lag is worth reading.
      !
      ! WHERE THE MODEL IS OUT OF DOMAIN.  Below the helium ionization
      ! front, in the shielded molecular layer, both channels are slow
      ! against the residence time and the omitted formation is what would
      ! set the CO abundance.  The rates are still evaluated there -- they
      ! are physics, not a switch -- and the cells are counted; over the
      ! run's own physical time a destruction that slow cannot remove what
      ! an equally omitted formation cannot rebuild, so the transported
      ! value stands, which is the correct behaviour of a transport
      ! operator in a quenched layer.
      src(ic_CO)  = -(rk_D1_Hep_CO()*n_heii + cph_kco(j))*nc(ic_CO)
      ! Row (1) of mol_heh_rows is the complete proton balance and is taken
      ! whole: photoionization of H including the secondaries, the
      ! dissociative and double photoionization branches of H2, the H2+ + H
      ! and He+ + H2 channels and the He(2^3S) Penning term, less radiative
      ! recombination and less the H+ + H2 reactions R10 and R13. There is
      ! no proton chemistry written here that the local solve does not
      ! already have; the difference is only where the row is evaluated.
      src(ic_Hp)  = 0.0d0
      if (ionization_transport) src(ic_Hp) = fv(1) + cx_sH(1)
      ! Rows (2) and (3) of the same block are the complete He+ and He++
      ! balances -- photoionization and collisional ionization of He I and
      ! of He+, the He(2^3S) channels, radiative and dielectronic
      ! recombination, the molecular sinks R17/R23 and the He+ + CO
      ! channel -- taken whole, exactly as row (1) is.  There is no helium
      ! chemistry written here that the local solve does not already have;
      ! the difference is only where the row is evaluated.
      src(ic_HeII)  = 0.0d0
      src(ic_HeIII) = 0.0d0
      if (carrier_solved(ic_HeII)) then
         if (thereis_mol) then
            src(ic_HeII)  = fv(2) + cx_sHe(1)
            src(ic_HeIII) = fv(3) + cx_sHe(2)
         else
            ! heh_tr_rows writes the SUMMED He I balance in row (2), He I
            ! gain positive, and the He++ balance in row (3).  A helium
            ! nucleus that leaves He I enters He+, and one that leaves He+
            ! enters He++, so the He+ balance is minus the He I balance
            ! less the He++ balance.  The three stages of the element then
            ! sum to zero source, which is what makes the transport of the
            ! partition conserve helium nuclei.
            src(ic_HeII)  = -fv(2) - fv(3) + cx_sHe(1)
            src(ic_HeIII) =  fv(3)         + cx_sHe(2)
         endif
      endif
      ! THE He 2^3S ROW IS THE LEVEL BALANCE THE SWEEP SOLVES, taken
      ! whole: tr_triplet_row (ion_residual_core), row (4) of the atomic
      ! rows and row (8) of the molecular ones, the latter with the
      ! He(2^3S) + H2 quench.  Production positive, as the stage rows.
      ! A transition between 2^3S and the singlet moves no nucleus out of
      ! He I, so the summed He I balance of row (2) and hence the He II
      ! source above are the same whichever of the two levels the neutral
      ! sits in; what the level changes there is which of its own channels
      ! (triplet photoionization, collisional ionization) runs, at the
      ! trial n_heiTR.  The metal charge exchange has no 2^3S channel
      ! (its helium reactant is the ground singlet).
      src(ic_HeTR) = 0.0d0
      if (carrier_solved(ic_HeTR)) then
         if (thereis_mol) then
            src(ic_HeTR) = fv(8)
         else
            src(ic_HeTR) = fv(4)
         endif
      endif
