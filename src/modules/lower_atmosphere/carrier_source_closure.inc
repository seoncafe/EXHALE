!  THE DENSITY CLOSURE OF THE CARRIER SOURCE, ONE SOURCE TEXT FOR TWO
!  ARITHMETICS: included by carrier_source (real*8 trial densities: the
!  residual) and carrier_source_dual (dual trial densities: the exact
!  derivatives of the Jacobian) of diffusive_photochemistry.  It forms, from
!  the trial carrier densities nc and the frozen background of the step, the
!  species densities the chemistry rows are evaluated at: the proton, the
!  molecular ions, the He 2^3S level and the helium stages (trial or
!  frozen), the free electrons, atomic hydrogen and free atomic oxygen (the
!  closures), and loads the cell state of the last sweep.

      ! THE PROTON IS EITHER FROZEN BACKGROUND OR THE TRIAL UNKNOWN.
      ! With the ionization state carried, n(H+) is what this Newton is
      ! solving for in this cell, so the row has to be evaluated at the
      ! trial value and not at the value the last sweep left behind --
      ! otherwise the recombination sink and the H closure below both refer
      ! to a different proton than the one the row is moving, and the
      ! Jacobian column for it is zero.
      if (ionization_transport) then
         n_hii = nc(ic_Hp)
      else
         n_hii = cbg_nhii(j)
      endif
      n_h2p   = cbg_nh2p(j)
      n_h3p   = cbg_nh3p(j)
      n_hehp  = cbg_nhehp(j)
      if (carrier_solved(ic_H2p)) n_h2p = nc(ic_H2p)
      if (carrier_solved(ic_H3p)) n_h3p = nc(ic_H3p)
      if (carrier_solved(ic_HeHp)) n_hehp = nc(ic_HeHp)
      ! THE He 2^3S POPULATION IS FROZEN BACKGROUND OR THE TRIAL UNKNOWN,
      ! for the reason the proton is: carried ("He 2^3S transport"), its
      ! row is solved here, and the triplet photoionization, collisional
      ! ionization and Penning terms of the He II, H II and singlet rows
      ! below then see the population the row is moving.
      n_heiTR = cbg_nheiTR(j)
      if (carrier_solved(ic_HeTR)) n_heiTR = nc(ic_HeTR)
      ! THE TWO IONIZED HELIUM STAGES ARE EITHER FROZEN BACKGROUND OR THE
      ! TRIAL UNKNOWNS, for the reason the proton is: with the helium
      ! ionization state carried, n(He+) and n(He++) are what this Newton
      ! is solving for in this cell, and rows evaluated at the last sweep's
      ! values would give the Jacobian columns of those rows a zero
      ! derivative.  The neutral helium then closes the element budget, the
      ! same closure the hydrogen has below: the helium AVAILABLE TO THE
      ! STAGES less the two ionized stages is the ground singlet, and the
      ! neutral stage is that singlet plus the frozen metastable.  The
      ! helium held in HeH+ and the population of the He 2^3S level are
      ! RESERVED from that budget rather than subtracted and clipped
      ! afterwards (carrier_helium_available_to_stages): on the admissible
      ! set of the two stages both differences are nonnegative, and the
      ! projection that defines that set is the one the limiter and every
      ! Newton trial go through.
      if (carrier_solved(ic_HeII)) then
         n_heii  = nc(ic_HeII)
         n_heiii = nc(ic_HeIII)
         call carrier_helium_neutral_of_partition(j, n_heii, n_heiii,    &
                                                  n_heiTR, n_hei,        &
                                                  n_heiSI)
         if (carrier_solved(ic_HeHp)) then
            n_hei = n_hei + cbg_nhehp(j) - n_hehp
            n_heiSI = n_heiSI + cbg_nhehp(j) - n_hehp
         endif
      else
         ! The frozen background's own partition.  Its singlet is the
         ! floored difference of the two columns the sweep wrote, the
         ! floor and its reason stated at he_ground_singlet_density
         ! (a restart whose two columns come from different solves can put
         ! the difference on the wrong side of zero); it is written out
         ! here rather than called because this routine runs inside the
         ! threaded Jacobian assembly and that function keeps a count.
         n_heii  = cbg_nheii(j)
         n_heiii = cbg_nheiii(j)
         n_hei   = cbg_nhei(j)
         n_heiSI = positive_part(n_hei - n_heiTR)
      endif
      ! THE ELECTRONS THAT RECOMBINE WITH THE TRIAL PROTONS ARE THE TRIAL
      ! PROTONS' OWN. The rest of the electron budget -- He+, He++, the
      ! molecular ions, the metals -- is background and stays frozen with
      ! it, but holding the H+ electrons frozen while H+ itself moves would
      ! break charge neutrality inside the solve and, where H+ carries the
      ! electrons (the wind of this configuration), would evaluate the
      ! recombination sink a_hii n_e n_hii at an n_e belonging to a
      ! different ionization state. n_e is a sum over charges, so replacing
      ! one term by its trial value is exact.
      n_e     = cbg_ne(j)
      if (ionization_transport)                                              &
         n_e = positive_part(cbg_ne(j) - cbg_nhii(j) + n_hii)
      ! The helium electrons of the same sum, He+ once and He++ twice, moved
      ! to their trial values for the same reason.  n_e is a sum over
      ! charges, so replacing one term by its trial value is exact.
      if (carrier_solved(ic_HeII))                                           &
         n_e = positive_part(n_e - cbg_nheii(j)  + n_heii                              &
                       - 2.0d0*cbg_nheiii(j) + 2.0d0*n_heiii)

      if (carrier_solved(ic_H2p)) n_e = n_e + n_h2p - cbg_nh2p(j)
      if (carrier_solved(ic_H3p)) n_e = n_e + n_h3p - cbg_nh3p(j)
      if (carrier_solved(ic_HeHp)) n_e = n_e + n_hehp - cbg_nhehp(j)

      ! Atomic H closes the hydrogen budget: what the element has left after
      ! the ions, the molecular ions and the trial carriers.  This is the
      ! coupling that makes the H2 source a function of the unknown.
      ! nH_free is the element headroom, so the frozen ions come off here.
      ! With the proton carried it is one of the carriers this closure
      ! subtracts (one H nucleus each), and hydrogen_available_to_carriers
      ! has stopped subtracting it as frozen background.
      if (ionization_transport) then
         n_hi = positive_part(hydrogen_available_to_carriers(j, nH_free)           &
                    - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O)     &
                    - nc(ic_Hp))
      else
      n_hi = positive_part(hydrogen_available_to_carriers(j, nH_free)              &
                 - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O))
      endif
      if (carrier_solved(ic_H2p)) n_hi = n_hi - 2.0d0*n_h2p
      if (carrier_solved(ic_H3p)) n_hi = n_hi - 3.0d0*n_h3p
      if (carrier_solved(ic_HeHp)) n_hi = n_hi - n_hehp
      if (carrier_solved(ic_HeHp) .and. .not. carrier_solved(ic_HeII)) then
         n_hei = n_hei + cbg_nhehp(j) - n_hehp
         n_heiSI = n_heiSI + cbg_nhehp(j) - n_hehp
      endif
      ! Free atomic oxygen closes the oxygen budget the same way: the
      ! element total less its ionized stages and less the carriers.
      n_o0 = positive_part(nO_free - cbg_nOion(j) - nc(ic_OH) - nc(ic_H2O)         &
                         - nc(ic_CO))

      ! The cell state of the last sweep, with only the photolysis rates
      ! refreshed: n_ofam and n_co are carried over unchanged because the
      ! rows called below take the densities they need as arguments.
      ieq_cell       = bg_cell(j)
