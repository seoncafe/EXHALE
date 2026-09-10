      program e1_h2_channel_check
      ! Standalone check of src/modules/functions/h2_photo_channels.f90,
      ! the E1 deliverable of docs/open_defects_20260903_review.md section 6.
      ! It is NOT part of the EXHALE build; compile it by hand, from
      ! src/tests/e1_h2/:
      !
      !   gfortran -O2 -o e1_h2_channel_check.x \
      !       stub_global_parameters.f90 \
      !       ../../modules/functions/cross_sec_h2_only.f90 \
      !       ../../modules/functions/h2_photo_channels.f90 \
      !       e1_h2_channel_check.f90
      !   ./e1_h2_channel_check.x > e1_h2_channel_check.out
      !
      ! RELATION TO THE PHASE A LEDGER.  Item A3 of the same review asks
      ! for "an H2 photoevent ledger, requiring conservation of H nuclei
      ! and charge for every final-state channel".  Phase A implements it
      ! as test E3, on the TWO-way branching, and its own comment hands the
      ! gate to this file: "These are REPORTED, not gated ... the gate
      ! belongs with the explicit-channel cross sections of phase E1, not
      ! here."  The two are therefore a pair and not a duplicate: E3 sizes
      ! the unresolved channels, test 2 below gates them once resolved.
      ! The invariants are E3's; they are not redefined here, only applied
      ! to four channels instead of two.  E3's two reported numbers are
      ! reconciled against this construction in section 151.3b of
      ! docs/Update_EXHALE_stage1.md -- note that E3's double-ionization figure is
      ! PROTON-weighted (0.20 f_di) where test 5 below is event-weighted,
      ! a factor 2 that is the P31b defect itself and not a disagreement.
      !
      ! What it checks, in order:
      !   1. sum of the four channel cross sections = the total absorption
      !      cross section, at machine precision, on a fine energy grid
      !      that straddles every join.
      !   2. the stoichiometric table conserves H nuclei and charge for
      !      every channel (the A3 ledger invariant).
      !   3. continuity of each channel across the joins of the three
      !      inputs: 15.4/18/85 eV (sigma_H2), 33/41 eV (f_n),
      !      51.4/80/110 eV (the double-ionization model).
      !   4. the default model 'off' reproduces the pre-E1 two-way split
      !      f_di exactly, so a run with the new channels off is unchanged.
      !   5. a rate-weighted ledger over a power-law spectrum: H nuclei and
      !      charge produced per unit time must balance the H2 destroyed.

      use h2_photo_channels
      use Cross_sections, only: sigma_H2, frac_H2_dissociative_ionization,   &
                               sigma_H2_k_hi
      implicit none
      integer, parameter :: dp = kind(1.0d0)
      real(dp) :: E, st, sig(n_h2_channels), s, worst, rel
      real(dp) :: dH2, dH2p, dH, dHp, dele, hn, ch
      integer  :: i, ich, nbad
      character(len=16) :: models(3)

      models = (/ 'off             ', 'chung80         ', 'yan_rho         ' /)

      write(*,'(a)') '=============================================================='
      write(*,'(a)') ' E1: H2 photoabsorption as mutually exclusive channels'
      write(*,'(a)') '=============================================================='
      write(*,*)

      ! ---------------------------------------------------------------
      write(*,'(a)') '1. CHANNEL SUM = TOTAL ABSORPTION CROSS SECTION'
      write(*,'(a)') '   max |sum(sigma_ch)/sigma_H2 - 1| over 15.4-2000 eV'
      do i = 1, 3
         h2_double_ionization_model = models(i)
         worst = 0.0d0
         E = 15.4d0
         do while (E .lt. 2000.0d0)
            st = sigma_H2(E)
            call h2_channel_cross_sections(E, st, sig)
            s = sum(sig)
            if (st .gt. 0.0d0) then
               rel = abs(s/st - 1.0d0)
               if (rel .gt. worst) worst = rel
            endif
            E = E + 0.01d0
         enddo
         write(*,'(4x,a16,es12.3)') models(i), worst
      enddo
      write(*,*)

      ! ---------------------------------------------------------------
      write(*,'(a)') '2. STOICHIOMETRY: H NUCLEI AND CHARGE PER PHOTON EVENT'
      write(*,'(a)') '   channel                    H2   H2+    H    H+    e-'// &
                     '   dH_nuc  dcharge'
      nbad = 0
      do ich = 1, n_h2_channels
         call h2_channel_stoichiometry(ich, dH2, dH2p, dH, dHp, dele)
         hn = 2.0d0*dH2 + 2.0d0*dH2p + dH + dHp
         ch = dH2p + dHp - dele
         write(*,'(4x,a24,5f6.1,2f9.1)') h2_channel_name(ich),            &
              dH2, dH2p, dH, dHp, dele, hn, ch
         if (abs(hn) .gt. 0.0d0 .or. abs(ch) .gt. 0.0d0) nbad = nbad + 1
      enddo
      write(*,'(4x,a,i0)') 'channels violating an invariant: ', nbad

      ! ---------------------------------------------------------------
      write(*,*)
      write(*,'(a)') '3. CONTINUITY ACROSS THE JOINS'
      write(*,'(a)') '   sigma of each channel just below / just above'
      do i = 1, 3
         h2_double_ionization_model = models(i)
         write(*,'(4x,a,a)') 'model ', trim(models(i))
         call join(15.40d0);  call join(18.00d0);  call join(18.076d0)
         call join(33.00d0);  call join(41.00d0);  call join(51.40d0)
         call join(80.00d0);  call join(85.00d0);  call join(110.00d0)
         call join(124.00d0)
      enddo

      ! ---------------------------------------------------------------
      write(*,*)
      write(*,'(a)') '4. MODEL ''off'' REPRODUCES THE PRE-E1 TWO-WAY SPLIT'
      write(*,'(a)') '   outside 33-41 eV the old code had'
      write(*,'(a)') '     sigma_S = f_di*sigma_H2,  sigma_M = (1-f_di)*sigma_H2'
      h2_double_ionization_model = 'off'
      worst = 0.0d0
      E = 15.4d0
      do while (E .lt. 2000.0d0)
         if (E .lt. 32.99d0 .or. E .gt. 41.01d0) then
            st = sigma_H2(E)
            call h2_channel_cross_sections(E, st, sig)
            if (st .gt. 0.0d0) then
               rel = abs(sig(ICH_S) - frac_H2_dissociative_ionization(E)*st)/st
               if (rel .gt. worst) worst = rel
               rel = abs(sig(ICH_M)                                       &
                    - (1.0d0-frac_H2_dissociative_ionization(E))*st)/st
               if (rel .gt. worst) worst = rel
            endif
         endif
         E = E + 0.01d0
      enddo
      write(*,'(4x,a,es12.3)') 'max |new - old|/sigma_H2 = ', worst
      write(*,'(4x,a)') '(inside 33-41 eV they differ by construction: that'
      write(*,'(4x,a)') ' band is where the neutral channel is taken out)'

      ! ---------------------------------------------------------------
      write(*,*)
      write(*,'(a)') '5. RATE-WEIGHTED LEDGER over a power-law spectrum'
      write(*,'(a)') '   F(E) ~ E^-1 over 15.4-1240 eV, 400 log bins.'
      write(*,'(a)') '   Per unit n(H2), per unit time:'
      write(*,'(a)') '   model            R(H2 destroyed)   dH_nuclei    dcharge'
      do i = 1, 3
         h2_double_ionization_model = models(i)
         call ledger(models(i))
      enddo

      ! ---------------------------------------------------------------
      write(*,*)
      write(*,'(a)') '6. sigma_H2 ITSELF IS CONTINUOUS (E1b)'
      write(*,'(a)') '   The three-piece Yan fit stepped by x1.719 at 18 eV and'
      write(*,'(a)') '   x0.619 at 85 eV. sigma_H2 is now measured throughout:'
      write(*,'(a)') '   Backx+1976 Table 1 over 15.4-18 eV, Samson & Haddad'
      write(*,'(a)') '   1994 Table 1 over 18-300 eV, and only the E^-7/2'
      write(*,'(a)') '   tail above 300 eV is a rescaled fit.'
      write(*,'(a)') '   Scan for any step above the local slope:'
      ! The scan starts at 15.5 eV: the 15.4 eV IONIZATION THRESHOLD is a
      ! real step (0 -> 0.173) and is not a join.
      worst = 0.0d0
      rel = 0.0d0
      E = 15.5d0
      s = -1.0d0
      do while (E .lt. 2000.0d0)
         st = sigma_H2(E)
         if (s .gt. 0.0d0 .and. st .gt. 0.0d0) then
            if (abs(st/s - 1.0d0) .gt. worst) then
               worst = abs(st/s - 1.0d0)
               rel   = E
            endif
         endif
         s = st
         E = E + 1.0d-3
      enddo
      write(*,'(4x,a,es12.3,a,f9.3,a)')                                   &
           'largest fractional step per 1e-3 eV = ', worst,               &
           ' at E = ', rel, ' eV'
      write(*,'(4x,a)') '(a smooth cross section gives ~1e-4 here; either'
      write(*,'(4x,a)') ' surviving join would give 4e-1 or 6e-1)'
      write(*,'(4x,a,es12.3)') 'the 15.4 eV threshold step itself     = ',  &
           sigma_H2(15.401d0)
      write(*,'(4x,a,f10.6)') 'k_hi, the 300 eV continuity constant = ',   &
           sigma_H2_k_hi()

      write(*,*)
      write(*,'(a)') 'done.'

      contains

      subroutine join(E0)
      real(dp), intent(in) :: E0
      real(dp) :: sa(n_h2_channels), sb(n_h2_channels), ea, eb
      ea = E0 - 1.0d-6;  eb = E0 + 1.0d-6
      call h2_channel_cross_sections(ea, sigma_H2(ea), sa)
      call h2_channel_cross_sections(eb, sigma_H2(eb), sb)
      write(*,'(8x,a,f9.3,a,4(es11.3),a,4(es11.3))') 'E=', E0,           &
           '  below', sa, '   above', sb
      end subroutine join

      subroutine ledger(mdl)
      character(len=16), intent(in) :: mdl
      integer, parameter :: nb = 400
      real(dp) :: Elo, Ehi, lg, dlg, Ec, w, stt, sg(n_h2_channels)
      real(dp) :: rdes, rhn, rch
      real(dp) :: a2, a2p, aH, aHp, ae
      integer  :: k, jc
      Elo = 15.4d0;  Ehi = 1240.0d0
      dlg = (log(Ehi) - log(Elo))/dble(nb)
      rdes = 0.0d0;  rhn = 0.0d0;  rch = 0.0d0
      do k = 1, nb
         lg = log(Elo) + dlg*(dble(k) - 0.5d0)
         Ec = exp(lg)
         w  = (1.0d0/Ec)*Ec*dlg          ! E^-1 spectrum, log bin
         stt = sigma_H2(Ec)
         call h2_channel_cross_sections(Ec, stt, sg)
         do jc = 1, n_h2_channels
            call h2_channel_stoichiometry(jc, a2, a2p, aH, aHp, ae)
            rdes = rdes - w*sg(jc)*a2                 ! H2 removed
            rhn  = rhn  + w*sg(jc)*(2.0d0*a2 + 2.0d0*a2p + aH + aHp)
            rch  = rch  + w*sg(jc)*(a2p + aHp - ae)
         enddo
      enddo
      write(*,'(4x,a16,es16.6,2es13.3)') mdl, rdes, rhn, rch
      end subroutine ledger

      end program e1_h2_channel_check
