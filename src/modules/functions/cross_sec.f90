      module Cross_sections
	   ! Energy-dependent photoionization cross sections (1e-18 cm^2)
	
      use global_parameters
	
      implicit none

      ! Scaling energy of the Yan, Sadeghpour & Dalgarno (1998), ApJ 496,
      ! 1044, Eq. 19 sum-rule tail of the H2 photoabsorption cross section:
      ! they write the fit in x = E/I with the H2 ionization potential quoted
      ! as 15.4 eV.  It is a COEFFICIENT OF THAT FIT, valid only with the
      ! seven polynomial coefficients as published, and so it is not the
      ! measured threshold e_th_H2 that gates the cross section and charges
      ! the photoelectron.  Used above 300 eV only (the measured table covers
      ! 15.4-300 eV), and by sigma_H2_k_hi, which rescales the tail to meet
      ! that table at 300 eV.
      real*8, parameter :: E_I_yan = 15.4d0

      ! ---- Metal photoionization fits, and the range each is valid in ----
      !
      ! Verner, D. A., Ferland, G. J., Korista, K. T., & Yakovlev, D. G.
      ! 1996, ApJ 465, 487, Table 1 and their Eqs. (1)-(4), and
      ! Verner, D. A., & Yakovlev, D. G. 1995, A&AS 109, 125, Table 1 and
      ! their Eq. (1).  The 1995 form is the 1996 one with y_0 = y_1 = 0
      ! and with the orbital quantum number l of the shell kept in
      ! Q = 5.5 + l - P/2; sigma_VFKY96 and sigma_Verner96 are the two.
      !
      ! WHY TWO FITS PER ION.  The 1996 fit is an OUTER-SHELL fit and its
      ! table carries an E_max, the energy at which the outermost INNER
      ! shell opens.  Above it the fit is outside the range it was made in,
      ! and Verner's own routine phfit2 leaves it there for the 1995
      ! subshell fit of the same outer shell (its variable einn, which is
      ! that inner-shell threshold and equals the E_max column row for
      ! row).  The photon grid of this code runs to 1240 eV while E_max is
      ! 34 to 558 eV, so the outer shell is evaluated well past E_max
      ! unless the handover is made.  Extrapolating the 1996 fit there is
      ! not a small error and it is not always conservative: the Fe I row
      ! has P = 7.923, hence Q = 1.5385, so the fit RISES as E^(2-Q) and
      ! gives 24.9 Mb at 1240 eV against a true Fe I total of 0.51 Mb.
      !
      ! THE SUBSHELLS BELOW THE OUTER ONE ARE IN THE SUM, AND WHAT THE
      ! MODEL DOES WITH THEM.  metal_photoion_sigma returns the shell SUM,
      ! which is what phfit2's caller obtains by summing that routine over
      ! is = 1 ... 7.  Above a K or L edge the cross section is therefore
      ! the whole photoabsorption of the ion; the outer shell alone was
      ! 0.099 per cent (C II) to 7.0 per cent (Mg I) of it at 1240 eV.
      !
      ! WHAT AN ABSORPTION IN EACH SHELL DOES.  An absorption in the outer
      ! shell ejects one electron of h nu - E_th and advances the ion one
      ! stage.  An absorption in a subshell s below it ejects the
      ! photoelectron h nu - E_th,s (E_th,s the threshold of that shell,
      ! mph_sub_e_th) and leaves a vacancy of energy E_th,s - I_i above the
      ! next ion's ground state.  The vacancy relaxes:
      !  * by Auger (autoionizing) transitions where it lies above the
      !    double-ionization threshold I_i + I_(i+1): one or more further
      !    electrons leave, the ion ends two or more stages up, and the
      !    Auger electrons carry what the ionization energies of those
      !    electrons and the fluorescence photons leave of E_th,s - I_i;
      !  * radiatively where it lies below that threshold (Na I 2p,
      !    Mg II 2p, K I 3p, Ca II 3p and 3s, and every valence 2s, 3s, 3d
      !    hole): the ion ends one stage up, excited by E_th,s - I_i, which
      !    leaves as line radiation.
      ! The probability of each number of ejected electrons and the photon
      ! energy of the cascade are Kaastra & Mewe (1993, A&AS 97, 443),
      ! Tables 2 and 3, in the block mph_sub_p_multi ... mph_sub_e_ion_multi
      ! below, with the energy of the vacancy closed exactly (their Eq. 5
      ! on this code's thresholds):
      !     E_th,s = I_i + e_ion_multi + e_auger + e_fluor .
      ! What each part becomes: the photoelectron and the Auger electrons
      ! are the electron energy (electron_energy_degradation divides it
      ! between heat and secondary ionization); I_i and, for a neutral whose
      ! next stage this code photoionizes, I_(i+1) are the ionization energy
      ! the stage balance stores; the photons leave the cell.
      !
      ! STAGES ABOVE THE CARRIED ONES.  species_table carries three stages
      ! of C, N, O, Mg, Si, Ca, Fe and two of Na, K, S, and the top stage is
      ! not photoionized.  A neutral of a three-stage element that loses two
      ! or more electrons goes to its third stage (the balance rows carry
      ! that jump: System_HeH_metals, met_g02).  Every other multiple
      ! ionization (a singly ionized absorber, a neutral of Na, K or S, and
      ! the third and further electrons of any absorber) would leave the
      ! carried stages: the ion is counted in the top carried stage, its
      ! charge there is lower than the true one by the electrons beyond it,
      ! and the ionization energy of those electrons is spent (it is not
      ! heat) but not stored by the balance: physically it returns as the
      ! recombination radiation of the stages this code does not carry.
      ! The electrons the top stage does not count are missing from n_e;
      ! for trace metals that is below the precision n_e is solved to.
      !
      ! FLUORESCENCE AND ITS YIELDS.  The fluorescence of a K vacancy is
      ! a small branch at the energies this grid reaches: READ from Krause
      ! 1979, J. Phys. Chem. Ref. Data 8, 307, Table 3 (p. 315), omega_K =
      ! 2.8e-3, 5.2e-3, 8.3e-3 and 0.023 for C, N, O and Na, the only K
      ! shells below the default e_top = 1240 eV; Kaastra & Mewe normalize
      ! their radiative K rates to the yields of Bambynek et al. (1972)
      ! (their Sect. 2.1).  Their Table 3 is the whole photon output of the
      ! cascade, the valence lines of the final ion included (the 2s-2p
      ! lines of 5-15 eV that follow a C, N or O KLL decay carry most of
      ! e_fluor there); a line weaker than 5e-5 photons per vacancy is not
      ! listed in it and is inside e_auger here.
      !
      ! APPROXIMATIONS, with their range.  (i) The 2p and 3p subshells of
      ! Verner's fits are Kaastra & Mewe's L2 + L3 and M2 + M3 combined in
      ! the ratio 1 : 2 of their occupancies.  (ii) Kaastra & Mewe decide
      ! each Auger transition on Lotz energy levels; the Ca II 3s vacancy,
      ! which their Table 2 lets autoionize, lies at E_th,s - I_i = 48.5 eV,
      ! below the Ca III potential of 50.91 eV (NIST), and is taken as
      ! radiative.  (iii) The Auger energy is the mean over the cascades of
      ! a vacancy; how it is shared among several Auger electrons is not
      ! carried, and the degradation treats it as one electron.  (iv) The
      ! fluorescence photons and the de-excitation lines leave the cell;
      ! their absorption elsewhere (a 50.3 eV Mg III line ionizes H I and
      ! He I) is not followed, like every photon this code does not
      ! transport.  (v) The yields are those of the free ion in its ground
      ! state; Kaastra & Mewe's Sect. 3 compares them with the independent
      ! Fe calculation of Jacobs & Rozsnyai (1986) (their Fig. 1).
      ! src/utils/inner_shell_relaxation_table.py emits the block from the
      ! tables and, with --check, compares it back (the physics probe runs
      ! that check).
      !
      ! HOW MUCH THIS MOVES A RUN.  MEASURED with the incident spectrum as
      ! the weight and no attenuation: above 100 eV hydrogen and helium
      ! carry nearly all of the absorption (the metals were 2.2 per cent of
      ! it summed over the column in wasp_full and in mol_metals before the
      ! subshells were added), and the whole 100-1240 eV band is 0.37 per
      ! cent of the 13.6-1240 eV absorption in either run.
      !
      ! The rows below are the published tables digit for digit, from the
      ! author's own distribution (references/verner_photo/photo.dat for
      ! 1996, table1.dat for the 1995 outer shell, and phfit2.f's own PH1
      ! DATA statements for the subshells, whose thresholds that file
      ! states were "slightly improved" over the printed 1995 table in a
      ! few low-ionized species: MEASURED, the two differ in the threshold
      ! column of Mg II 1s, Si II 1s, Ca II 1s and 2s, and Fe II 1s and 2s,
      ! and nowhere else).  src/utils/metal_photoion_table.py
      ! emits this block and, with --check, compares it back against those
      ! two files; the physics probe runs that check.  E_max = 0 for K I
      ! is phfit2's einn = 0 for potassium, so K I takes the 1995 fit at
      ! every energy and has no 1996 row (the 1996 table has no potassium);
      ! its 1996 columns are zero and are never reached.
      !
      ! THRESHOLDS.  A threshold only GATES its fit and does not enter the
      ! fitted expression, so it is the energy the photon grid uses as the
      ! ion's band edge (species_table's mion_ethr) and charges the
      ! photoelectron h nu - E_th against.  Mg I, Mg II and Fe II therefore
      ! take the NIST ionization potential where the fit tables print a
      ! rounded value (15.04 for Mg II, 16.19 for Fe II).
      ! Column order (= species_table's mion_iphot): 
      !   CI, CII, OI, OII, NI, NII
      !   MgI, MgII, SiI, SiII, CaI, CaII
      !   NaI, KI, SI, FeI, FeII
      integer, parameter :: n_metal_photo = 17
      ! Threshold that gates the fit [eV].
      real*8, parameter :: mph_e_th(n_metal_photo) = [ &
           1.126d+01  , 2.438d+01  , 1.362d+01  , 3.512d+01  , 1.453d+01  , 2.960d+01, &
           e_th_MgI   , e_th_MgII  , 8.152d+00  , 1.635d+01  , 6.113d+00  , 1.187d+01, &
           5.139d+00  , 0.4341d+01 , 1.036d+01  , 7.902d+00  , 1.6199d+01 ]
      ! E_max [eV]: above it the outer shell takes the 1995 fit.
      real*8, parameter :: mph_e_max(n_metal_photo) = [ &
           2.910d+02  , 3.076d+02  , 5.380d+02  , 5.581d+02  , 4.048d+02  , 4.236d+02, &
           5.490d+01  , 6.569d+01  , 1.060d+02  , 1.186d+02  , 3.443d+01  , 4.090d+01, &
           3.814d+01  , 0.000d+00  , 1.700d+02  , 6.600d+01  , 7.617d+01 ]
      ! Verner et al. 1996 Table 1, E_0 [eV].
      real*8, parameter :: mph96_E0(n_metal_photo) = [ &
           2.144d+00  , 4.058d-01  , 1.240d+00  , 1.386d+00  , 4.034d+00  , 6.128d-02, &
           1.197d+01  , 8.139d+00  , 2.317d+01  , 2.556d+00  , 1.278d+01  , 1.553d+01, &
           6.139d+00  , 0.000d+00  , 1.808d+01  , 5.461d-02  , 1.761d-01 ]
      ! Verner et al. 1996 Table 1, sigma_0 [Mb].
      real*8, parameter :: mph96_s0(n_metal_photo) = [ &
           5.027d+02  , 8.709d+00  , 1.745d+03  , 5.967d+01  , 8.235d+02  , 1.944d+00, &
           1.372d+08  , 3.278d+00  , 2.506d+01  , 4.140d+00  , 5.370d+05  , 1.064d+07, &
           1.601d+00  , 0.000d+00  , 4.564d+04  , 3.062d-01  , 4.365d+03 ]
      ! Verner et al. 1996 Table 1, y_a.
      real*8, parameter :: mph96_ya(n_metal_photo) = [ &
           6.216d+01  , 1.261d+02  , 3.784d+00  , 3.175d+01  , 8.033d+01  , 8.163d+02, &
           2.228d-01  , 4.341d+07  , 2.057d+01  , 1.337d+01  , 3.162d-01  , 7.790d-01, &
           6.148d+03  , 0.000d+00  , 1.000d+00  , 2.671d+07  , 6.298d+03 ]
      ! Verner et al. 1996 Table 1, P.
      real*8, parameter :: mph96_P(n_metal_photo) = [ &
           5.101d+00  , 8.578d+00  , 1.764d+01  , 8.943d+00  , 3.928d+00  , 8.773d+00, &
           1.574d+01  , 3.610d+00  , 3.546d+00  , 1.191d+01  , 1.242d+01  , 2.130d+01, &
           3.839d+00  , 0.000d+00  , 1.361d+01  , 7.923d+00  , 5.204d+00 ]
      ! Verner et al. 1996 Table 1, y_w.
      real*8, parameter :: mph96_yw(n_metal_photo) = [ &
           9.157d-02  , 2.093d+00  , 7.589d-02  , 1.934d-02  , 9.097d-02  , 1.043d+01, &
           2.805d-01  , 0.000d+00  , 2.837d-01  , 1.570d+00  , 4.477d-01  , 6.453d-01, &
           0.000d+00  , 0.000d+00  , 6.385d-01  , 2.069d+01  , 1.141d+01 ]
      ! Verner et al. 1996 Table 1, y_0.
      real*8, parameter :: mph96_y0(n_metal_photo) = [ &
           1.133d+00  , 4.929d+01  , 8.698d+00  , 2.131d+01  , 8.598d-01  , 4.280d+02, &
           0.000d+00  , 0.000d+00  , 1.672d-05  , 6.634d+00  , 1.012d-03  , 2.161d-03, &
           0.000d+00  , 0.000d+00  , 9.935d-01  , 1.382d+02  , 9.272d+01 ]
      ! Verner et al. 1996 Table 1, y_1.
      real*8, parameter :: mph96_y1(n_metal_photo) = [ &
           1.607d+00  , 3.234d+00  , 1.271d-01  , 1.503d-02  , 2.325d+00  , 2.030d+01, &
           0.000d+00  , 0.000d+00  , 4.207d-01  , 1.272d-01  , 1.851d-02  , 6.706d-02, &
           0.000d+00  , 0.000d+00  , 2.486d-01  , 2.481d-01  , 1.075d+02 ]
      ! Verner & Yakovlev 1995 Table 1, E_0 [eV].
      real*8, parameter :: mph95_E0(n_metal_photo) = [ &
           0.9435d+01 , 0.1094d+02 , 0.1391d+02 , 0.1745d+02 , 0.1164d+02 , 0.1827d+02, &
           0.9393d+01 , 0.8139d+01 , 0.2212d+02 , 0.2123d+02 , 0.7366d+01 , 0.4155d+01, &
           0.5968d+01 , 0.3824d+01 , 0.2975d+02 , 0.1277d+02 , 0.1014d+02 ]
      ! Verner & Yakovlev 1995 Table 1, sigma_0 [Mb].
      real*8, parameter :: mph95_s0(n_metal_photo) = [ &
           0.1152d+04 , 0.1792d+03 , 0.1220d+06 , 0.5186d+03 , 0.1029d+05 , 0.1724d+03, &
           0.3034d+01 , 0.3278d+01 , 0.1845d+03 , 0.6975d+02 , 0.2373d+01 , 0.2235d+01, &
           0.1460d+01 , 0.7363d+00 , 0.5644d+02 , 0.1468d+01 , 0.1084d+01 ]
      ! Verner & Yakovlev 1995 Table 1, y_a.
      real*8, parameter :: mph95_ya(n_metal_photo) = [ &
           0.5687d+01 , 0.3308d+02 , 0.1364d+01 , 0.1728d+02 , 0.2361d+01 , 0.8893d+02, &
           0.2625d+08 , 0.4341d+08 , 0.3849d+01 , 0.4907d+01 , 0.2082d+03 , 0.1595d+05, &
           0.2557d+08 , 0.2410d+08 , 0.1321d+02 , 0.1116d+06 , 0.2562d+05 ]
      ! Verner & Yakovlev 1995 Table 1, P.
      real*8, parameter :: mph95_P(n_metal_photo) = [ &
           0.6336d+01 , 0.4150d+01 , 0.1140d+02 , 0.4995d+01 , 0.8821d+01 , 0.3348d+01, &
           0.3923d+01 , 0.3610d+01 , 0.9721d+01 , 0.9525d+01 , 0.4841d+01 , 0.4313d+01, &
           0.3789d+01 , 0.4427d+01 , 0.7513d+01 , 0.4112d+01 , 0.4167d+01 ]
      ! Verner & Yakovlev 1995 Table 1, y_w.
      real*8, parameter :: mph95_yw(n_metal_photo) = [ &
           0.4474d+00 , 0.5276d+00 , 0.4103d+00 , 0.2182d-01 , 0.4239d+00 , 0.4209d+00, &
           0.0000d+00 , 0.0000d+00 , 0.2921d+00 , 0.3169d+00 , 0.5841d-03 , 0.3539d+00, &
           0.0000d+00 , 0.2049d-03 , 0.2621d+00 , 0.3238d-01 , 0.1598d-01 ]
      ! Orbital quantum number of the outer shell, in
      ! Q = 5.5 + l - P/2 of the 1995 form.
      integer, parameter :: mph_l(n_metal_photo) = [ &
           1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0 ]
      ! ---- The subshells below the outer one ----
      !
      ! Slots per ion, padded with sigma_0 = 0 rows that the sum
      ! never reaches (the loop runs to mph_n_sub).  Shell order
      ! is phfit2's: 1s, 2s, 2p, 3s, 3p, 3d.  Every row is the
      ! 1995 five-parameter form, evaluated above mph_sub_e_on.
      integer, parameter :: n_metal_subshell = 6
      ! Subshells this ion carries electrons in.
      integer, parameter :: mph_n_sub(n_metal_photo) = [ &
           2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 3, 5, 4, 6, 6 ]
      ! Turn-on [eV] = max(shell threshold, E_max).
      real*8, parameter :: mph_sub_e_on(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           2.910d+02  , 2.910d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           3.076d+02  , 3.076d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           5.380d+02  , 5.380d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           5.581d+02  , 5.581d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           4.048d+02  , 4.048d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           4.236d+02  , 4.236d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           1.311d+03  , 9.400d+01  , 5.490d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           1.320d+03  , 9.881d+01  , 6.569d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           1.846d+03  , 1.560d+02  , 1.060d+02  , 1.060d+02  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           1.848d+03  , 1.619d+02  , 1.186d+02  , 1.186d+02  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           4.043d+03  , 4.425d+02  , 3.523d+02  , 4.830d+01  , 3.443d+01  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           4.047d+03  , 4.445d+02  , 3.638d+02  , 6.037d+01  , 4.090d+01  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           1.079d+03  , 7.084d+01  , 3.814d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           3.614d+03  , 3.843d+02  , 3.014d+02  , 4.080d+01  , 2.466d+01  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           2.477d+03  , 2.350d+02  , 1.700d+02  , 1.700d+02  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           7.124d+03  , 8.570d+02  , 7.240d+02  , 1.040d+02  , 6.600d+01  , 6.600d+01, &   ! FeI   1s 2s 2p 3s 3p 3d
           7.140d+03  , 8.608d+02  , 7.341d+02  , 1.102d+02  , 7.617d+01  , 7.617d+01 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Verner & Yakovlev 1995 subshell fit, E_0 [eV].
      real*8, parameter :: mph_sub_E0(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           8.655d+01  , 1.026d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           9.113d+01  , 2.991d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           1.774d+02  , 1.994d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           1.690d+02  , 1.759d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           1.270d+02  , 1.482d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           1.242d+02  , 1.094d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           2.711d+02  , 4.587d+01  , 4.937d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           3.709d+02  , 5.142d+01  , 4.940d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           5.322d+02  , 7.017d+01  , 7.808d+01  , 1.413d+01  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           4.580d+02  , 6.738d+01  , 7.154d+01  , 1.364d+01  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           6.947d+02  , 1.201d+02  , 1.529d+02  , 3.012d+01  , 4.487d+01  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           7.997d+02  , 1.335d+02  , 1.448d+02  , 3.176d+01  , 4.498d+01  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           4.216d+02  , 4.537d+01  , 3.655d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           1.171d+03  , 1.602d+02  , 2.666d+02  , 2.910d+01  , 4.138d+01  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           8.114d+02  , 1.047d+02  , 9.152d+01  , 1.916d+01  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           8.044d+02  , 5.727d+01  , 2.948d+02  , 4.334d+01  , 7.630d+01  , 1.407d+01, &   ! FeI   1s 2s 2p 3s 3p 3d
           8.931d+02  , 1.431d+02  , 2.281d+02  , 4.663d+01  , 7.750d+01  , 1.933d+01 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Verner & Yakovlev 1995 subshell fit, sigma_0 [Mb].
      real*8, parameter :: mph_sub_s0(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           7.421d+01  , 4.564d+03  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           6.649d+01  , 1.184d+03  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           3.237d+01  , 2.415d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           3.584d+01  , 1.962d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           4.748d+01  , 7.722d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           5.002d+01  , 7.483d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           3.561d+01  , 1.671d+01  , 2.023d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           1.767d+01  , 1.290d+01  , 2.049d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           1.184d+01  , 1.166d+01  , 1.532d+02  , 1.166d+01  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           1.628d+01  , 1.263d+01  , 1.832d+02  , 8.454d+00  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           1.586d+01  , 1.010d+01  , 1.282d+02  , 7.227d+00  , 9.017d+01  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           1.168d+01  , 8.939d+00  , 1.446d+02  , 6.924d+00  , 7.314d+01  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           1.119d+01  , 1.142d+01  , 2.486d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           4.540d+00  , 6.389d+00  , 3.107d+01  , 6.377d+00  , 2.614d+01  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           6.649d+00  , 8.520d+00  , 1.883d+02  , 1.003d+01  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           2.055d+01  , 1.076d+01  , 7.191d+01  , 5.921d+00  , 6.298d+01  , 1.850d+04, &   ! FeI   1s 2s 2p 3s 3p 3d
           1.666d+01  , 9.289d+00  , 1.241d+02  , 5.434d+00  , 4.624d+01  , 3.679d+04 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Verner & Yakovlev 1995 subshell fit, y_a.
      real*8, parameter :: mph_sub_ya(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           5.498d+01  , 1.568d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           9.609d+01  , 3.085d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           3.812d+02  , 3.241d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           1.894d+02  , 4.020d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           1.380d+02  , 2.306d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           9.100d+01  , 2.793d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           2.374d+01  , 2.389d+01  , 1.079d+04  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           1.346d+02  , 5.775d+01  , 4.112d+03  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           2.580d+02  , 4.742d+01  , 5.765d+06  , 2.288d+01  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           7.706d+01  , 3.623d+01  , 3.537d+02  , 7.489d+01  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           2.563d+01  , 2.468d+01  , 2.217d+02  , 1.736d+02  , 1.465d+01  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           3.233d+01  , 2.914d+01  , 1.349d+02  , 5.246d+02  , 1.898d+01  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           5.642d+07  , 2.395d+02  , 3.222d+02  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           6.165d+03  , 1.044d+02  , 7.187d+06  , 2.229d+03  , 2.143d+02  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           3.734d+03  , 9.469d+01  , 7.193d+01  , 3.296d+01  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           3.633d+01  , 2.785d+01  , 3.219d+02  , 5.293d+01  , 1.479d+01  , 4.458d+00, &   ! FeI   1s 2s 2p 3s 3p 3d
           2.381d+01  , 2.026d+01  , 8.058d+01  , 9.271d+01  , 2.155d+01  , 3.564d+00 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Verner & Yakovlev 1995 subshell fit, P.
      real*8, parameter :: mph_sub_P(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           1.503d+00  , 1.085d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           1.338d+00  , 1.480d+01  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           1.083d+00  , 8.037d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           1.185d+00  , 7.999d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           1.252d+00  , 9.139d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           1.335d+00  , 9.956d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           1.952d+00  , 4.742d+00  , 2.960d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           1.225d+00  , 3.903d+00  , 2.995d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           1.102d+00  , 3.933d+00  , 2.639d+00  , 5.334d+00  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           1.385d+00  , 4.172d+00  , 3.133d+00  , 4.676d+00  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           1.966d+00  , 4.592d+00  , 3.087d+00  , 4.165d+00  , 7.498d+00  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           1.744d+00  , 4.269d+00  , 3.300d+00  , 3.771d+00  , 7.152d+00  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           7.736d-01  , 3.380d+00  , 3.570d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           8.392d-01  , 3.159d+00  , 2.067d+00  , 3.587d+00  , 5.631d+00  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           8.646d-01  , 3.346d+00  , 3.633d+00  , 5.038d+00  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           2.118d+00  , 6.635d+00  , 2.837d+00  , 5.129d+00  , 7.672d+00  , 1.691d+01, &   ! FeI   1s 2s 2p 3s 3p 3d
           2.262d+00  , 5.376d+00  , 3.572d+00  , 4.640d+00  , 7.138d+00  , 1.566d+01 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Verner & Yakovlev 1995 subshell fit, y_w.
      real*8, parameter :: mph_sub_yw(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CI    1s 2s
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! CII   1s 2s
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OI    1s 2s
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! OII   1s 2s
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NI    1s 2s
           0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NII   1s 2s
           0.000d+00  , 0.000d+00  , 1.463d-02  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgI   1s 2s 2p
           0.000d+00  , 0.000d+00  , 2.223d-05  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! MgII  1s 2s 2p
           0.000d+00  , 0.000d+00  , 2.774d-04  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! SiI   1s 2s 2p 3s
           0.000d+00  , 0.000d+00  , 2.870d-04  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! SiII  1s 2s 2p 3s
           0.000d+00  , 0.000d+00  , 3.343d-03  , 0.000d+00  , 2.754d-01  , 0.000d+00, &   ! CaI   1s 2s 2p 3s 3p
           0.000d+00  , 0.000d+00  , 3.358d-04  , 0.000d+00  , 2.735d-01  , 0.000d+00, &   ! CaII  1s 2s 2p 3s 3p
           0.000d+00  , 0.000d+00  , 1.465d-01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! NaI   1s 2s 2p
           0.000d+00  , 0.000d+00  , 5.274d-01  , 0.000d+00  , 2.437d-01  , 0.000d+00, &   ! KI    1s 2s 2p 3s 3p
           0.000d+00  , 0.000d+00  , 2.485d-01  , 0.000d+00  , 0.000d+00  , 0.000d+00, &   ! SI    1s 2s 2p 3s
           0.000d+00  , 0.000d+00  , 6.314d-02  , 0.000d+00  , 2.646d-01  , 4.039d-01, &   ! FeI   1s 2s 2p 3s 3p 3d
           0.000d+00  , 0.000d+00  , 1.223d-01  , 0.000d+00  , 2.599d-01  , 4.144d-01 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )
      ! Orbital quantum number l of each subshell.
      integer, parameter :: mph_sub_l(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           0          , 0          , 0          , 0          , 0          , 0, &   ! CI    1s 2s
           0          , 0          , 0          , 0          , 0          , 0, &   ! CII   1s 2s
           0          , 0          , 0          , 0          , 0          , 0, &   ! OI    1s 2s
           0          , 0          , 0          , 0          , 0          , 0, &   ! OII   1s 2s
           0          , 0          , 0          , 0          , 0          , 0, &   ! NI    1s 2s
           0          , 0          , 0          , 0          , 0          , 0, &   ! NII   1s 2s
           0          , 0          , 1          , 0          , 0          , 0, &   ! MgI   1s 2s 2p
           0          , 0          , 1          , 0          , 0          , 0, &   ! MgII  1s 2s 2p
           0          , 0          , 1          , 0          , 0          , 0, &   ! SiI   1s 2s 2p 3s
           0          , 0          , 1          , 0          , 0          , 0, &   ! SiII  1s 2s 2p 3s
           0          , 0          , 1          , 0          , 1          , 0, &   ! CaI   1s 2s 2p 3s 3p
           0          , 0          , 1          , 0          , 1          , 0, &   ! CaII  1s 2s 2p 3s 3p
           0          , 0          , 1          , 0          , 0          , 0, &   ! NaI   1s 2s 2p
           0          , 0          , 1          , 0          , 1          , 0, &   ! KI    1s 2s 2p 3s 3p
           0          , 0          , 1          , 0          , 0          , 0, &   ! SI    1s 2s 2p 3s
           0          , 0          , 1          , 0          , 1          , 2, &   ! FeI   1s 2s 2p 3s 3p 3d
           0          , 0          , 1          , 0          , 1          , 2 ], &   ! FeII  1s 2s 2p 3s 3p 3d
           [n_metal_subshell,n_metal_photo] )

      ! ---- Relaxation of a subshell vacancy (Kaastra & Mewe 1993) ----
      !
      ! The same slots as the mph_sub_* fits above.  Written by
      ! src/utils/inner_shell_relaxation_table.py, whose header states the
      ! sources, the shell map and the energetic rule; the physics is at
      ! the head of this module ("WHAT AN ABSORPTION IN EACH SHELL DOES").
      ! Closure, exact to the printed digits:
      !   mph_sub_e_th = I_i + mph_sub_e_ion_multi + mph_sub_e_auger
      !                  + mph_sub_e_fluor ,
      ! I_i the threshold mph_e_th of the ion.
      ! Threshold of the subshell itself, E_th,s [eV] (phfit2.f PH1, field
      ! 1).
      real*8, parameter :: mph_sub_e_th(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           2.910d+02     , 1.939d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CI
           3.076d+02     , 3.047d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CII
           5.380d+02     , 2.848d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OI
           5.581d+02     , 4.599d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OII
           4.048d+02     , 2.541d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NI
           4.236d+02     , 3.796d+01     , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NII
           1.311d+03     , 9.400d+01     , 5.490d+01     , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgI
           1.320d+03     , 9.881d+01     , 6.569d+01     , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgII
           1.846d+03     , 1.560d+02     , 1.060d+02     , 1.517d+01     , 0.0d0         , 0.0d0, &   ! SiI
           1.848d+03     , 1.619d+02     , 1.186d+02     , 2.240d+01     , 0.0d0         , 0.0d0, &   ! SiII
           4.043d+03     , 4.425d+02     , 3.523d+02     , 4.830d+01     , 3.443d+01     , 0.0d0, &   ! CaI
           4.047d+03     , 4.445d+02     , 3.638d+02     , 6.037d+01     , 4.090d+01     , 0.0d0, &   ! CaII
           1.079d+03     , 7.084d+01     , 3.814d+01     , 0.0d0         , 0.0d0         , 0.0d0, &   ! NaI
           3.614d+03     , 3.843d+02     , 3.014d+02     , 4.080d+01     , 2.466d+01     , 0.0d0, &   ! KI
           2.477d+03     , 2.350d+02     , 1.700d+02     , 2.130d+01     , 0.0d0         , 0.0d0, &   ! SI
           7.124d+03     , 8.570d+02     , 7.240d+02     , 1.040d+02     , 6.600d+01     , 1.470d+01, &   ! FeI
           7.140d+03     , 8.608d+02     , 7.341d+02     , 1.102d+02     , 7.617d+01     , 2.193d+01 ], &   ! FeII
           [n_metal_subshell,n_metal_photo] )
      ! Probability that the vacancy autoionizes, i.e. that two or more
      ! electrons leave (Kaastra & Mewe 1993, Table 2).
      real*8, parameter :: mph_sub_p_multi(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           9.97400000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CI
           9.98100000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CII
           9.90600000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OI
           9.89100000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OII
           9.94000000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NI
           9.92900000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NII
           9.99900000d-1 , 1.00000000d0  , 9.96700000d-1 , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgI
           9.65900000d-1 , 9.99900000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgII
           9.98800000d-1 , 1.00000000d0  , 9.97333333d-1 , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiI
           9.99300000d-1 , 1.00000000d0  , 9.95833333d-1 , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiII
           9.99900000d-1 , 1.00000000d0  , 1.00000000d0  , 1.00000000d0  , 9.95000000d-1 , 0.0d0, &   ! CaI
           9.86000000d-1 , 9.99000000d-1 , 1.00000000d0  , 0.0d0         , 0.0d0         , 0.0d0, &   ! CaII
           9.74000000d-1 , 9.99900000d-1 , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NaI
           9.89000000d-1 , 9.99300000d-1 , 1.00000000d0  , 9.99900000d-1 , 0.0d0         , 0.0d0, &   ! KI
           9.95400000d-1 , 9.99900000d-1 , 9.99566667d-1 , 0.0d0         , 0.0d0         , 0.0d0, &   ! SI
           9.95302349d-1 , 1.00000000d0  , 9.84600000d-1 , 1.00000000d0  , 9.99900000d-1 , 0.0d0, &   ! FeI
           9.95299530d-1 , 1.00000000d0  , 9.84533333d-1 , 1.00000000d0  , 9.99900000d-1 , 0.0d0 ], &   ! FeII
           [n_metal_subshell,n_metal_photo] )
      ! Mean kinetic energy of the Auger electrons of one absorption [eV],
      ! the remainder of (*).
      real*8, parameter :: mph_sub_e_auger(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           2.45856558d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CI
           2.19568517d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CII
           4.79735708d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OI
           4.49126597d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OII
           3.47690530d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NI
           3.32711922d2  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NII
           1.14309526d3  , 2.16873000d1  , 3.21154955d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgI
           1.10178074d3  , 3.63051436d0  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgII
           1.63521362d3  , 9.06880672d1  , 7.41871267d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiI
           1.55348775d3  , 6.66478916d1  , 5.97012375d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiII
           3.26632745d3  , 2.73375158d2  , 2.74709803d2  , 4.57949000d0  , 1.63961833d1  , 0.0d0, &   ! CaI
           3.22583840d3  , 2.42757999d2  , 2.55300933d2  , 0.0d0         , 0.0d0         , 0.0d0, &   ! CaII
           9.46973749d2  , 1.84131786d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NaI
           3.00960138d3  , 2.49977863d2  , 2.33930458d2  , 4.83389250d0  , 0.0d0         , 0.0d0, &   ! KI
           2.14344108d3  , 1.54915246d2  , 1.33561563d2  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SI
           4.58748986d3  , 6.34063545d2  , 6.31366342d2  , 5.12063791d1  , 4.18955566d1  , 0.0d0, &   ! FeI
           4.54217288d3  , 5.72919782d2  , 5.88819011d2  , 1.08323370d1  , 2.93166518d1  , 0.0d0 ], &   ! FeII
           [n_metal_subshell,n_metal_photo] )
      ! Mean energy leaving as photons [eV] (Kaastra & Mewe 1993 Table 3,
      ! sum of omega E; E_th,s - I_i for a vacancy that does not
      ! autoionize).
      real*8, parameter :: mph_sub_e_fluor(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           9.56683000d0  , 8.13000000d0  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CI
           1.58546900d1  , 6.09000000d0  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CII
           9.85442000d0  , 1.48600000d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OI
           1.95166600d1  , 1.08700000d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OII
           1.31570700d1  , 1.08800000d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NI
           1.41796400d1  , 8.36000000d0  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NII
           6.05989700d1  , 4.96317000d1  , 1.53120000d-1 , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgI
           8.68860000d1  , 8.90000000d-3 , 5.06550000d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgII
           1.07694610d2  , 8.33512000d0  , 7.35447333d0  , 7.01800000d0  , 0.0d0         , 0.0d0, &   ! SiI
           1.50030090d2  , 1.77054000d0  , 9.19531667d0  , 6.05000000d0  , 0.0d0         , 0.0d0, &   ! SiII
           5.42251620d2  , 2.75089700d1  , 9.09636667d0  , 2.57375100d1  , 1.10166667d-1 , 0.0d0, &   ! CaI
           5.51529850d2  , 2.33296200d1  , 3.04471333d1  , 4.85000000d1  , 2.90300000d1  , 0.0d0, &   ! CaII
           5.64938600d1  , 6.19000000d-3 , 3.30010000d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! NaI
           4.11980470d2  , 1.85945300d1  , 2.15322067d1  , 3.27000000d-3 , 2.03190000d1  , 0.0d0, &   ! KI
           2.06302280d2  , 1.27039900d1  , 2.75067000d0  , 1.09400000d1  , 0.0d0         , 0.0d0, &   ! SI
           2.25243645d3  , 1.72398100d1  , 1.13988867d1  , 6.35000000d-3 , 5.06333333d-3 , 6.79800000d0, &   ! FeI
           2.29861582d3  , 3.89464100d1  , 1.35467333d1  , 7.23000000d-3 , 6.41333333d-3 , 5.73100000d0 ], &   ! FeII
           [n_metal_subshell,n_metal_photo] )
      ! Mean ionization energy of the electrons beyond the first [eV] (NIST
      ! ASD; mion_ethr for a stage the code photoionizes).
      real*8, parameter :: mph_sub_e_ion_multi(n_metal_subshell,n_metal_photo) = &
           reshape( [ &
           2.43166120d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CI
           4.77967932d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! CII
           3.47898720d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OI
           5.43367426d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! OII
           2.94224000d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NI
           4.71084384d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NII
           9.96597677d1  , 1.50350000d1  , 1.49853845d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgI
           1.16298259d2  , 8.01355856d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! MgII
           9.49397697d1  , 4.88248128d1  , 1.63064000d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiI
           1.28132158d2  , 7.71315684d1  , 3.33534458d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SiII
           2.28307929d2  , 1.35502872d2  , 6.23808300d1  , 1.18700000d1  , 1.18106500d1  , 0.0d0, &   ! CaI
           2.57761750d2  , 1.66542381d2  , 6.61819340d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! CaII
           7.03933906d1  , 4.72816314d1  , 0.0d0         , 0.0d0         , 0.0d0         , 0.0d0, &   ! NaI
           1.88077152d2  , 1.11386607d2  , 4.15963349d1  , 3.16218375d1  , 0.0d0         , 0.0d0, &   ! KI
           1.16896635d2  , 5.70207642d1  , 2.33277669d1  , 0.0d0         , 0.0d0         , 0.0d0, &   ! SI
           2.76171694d2  , 1.97794645d2  , 7.33327716d1  , 4.48852709d1  , 1.61973801d1  , 0.0d0, &   ! FeI
           2.83012302d2  , 2.32734808d2  , 1.15535256d2  , 8.31614330d1  , 3.06479349d1  , 0.0d0 ], &   ! FeII
           [n_metal_subshell,n_metal_photo] )

	   contains
	
	   !----- Hydrogenic atoms -----! 
	
      pure double precision function sigma(E,Z,E_th)
      ! Hydrogenic photoionization cross section of an ion of nuclear charge
      ! Z, in 1e-18 cm^2.
      !
      ! E_th is the MEASURED ionization potential of that ion.  It both gates
      ! the cross section and sets the scale of the fit, so it must be the
      ! same energy the photon grid uses as a band edge and the photoelectron
      ! h nu - E_th is charged against; a caller passes the named global
      ! constant of its ion (e_th_HI for H I, e_th_HeII for He II).  Omitted,
      ! it is the hydrogenic scaling Z^2 e_th_HI, which is exact only for
      ! infinite nuclear mass and without relativistic corrections: for He II
      ! it gives 54.3937 eV against the measured 54.4178 eV, 0.044 per cent
      ! low, which would put a second definition of that threshold in the
      ! code and would misplace the (E_0/E)^4 scale by 0.18 per cent.
      real*8, intent(in) :: Z,E
      real*8, intent(in), optional :: E_th
      real*8 :: E_0,eps

      if (present(E_th)) then
         E_0 = E_th
      else
         E_0 = e_th_HI*Z*Z
      endif

      if (E.lt.E_0) then

         ! Below the threshold there is no bound-free absorption.  The gate is
         ! the threshold itself: a cross section that turned on a fraction
         ! below it (the retired 0.99999*E_0, 1.4e-4 eV for H I) is absorption
         ! by a photon that cannot ionize.
         sigma = 0.0d0

      else if (E.gt.E_0) then      ! For E > E_th

         ! Substitution
         eps = sqrt(E/E_0-1.0d0)               
         
         ! Cross section value
         sigma = 6.3d0/(Z*Z)*(E_0/E)**4.0d0        &
                  *exp(4.0d0-4.0d0*atan(eps)/eps)    &
                  /(1.0d0-exp(-2.0d0*pi/eps))
            
      else
         ! Exactly at the threshold: the eps -> 0 limit of the expression
         ! above, which is finite and equal to the threshold cross section.
         sigma = 6.3d0/(Z*Z) 
      
      endif

      ! End of function
      end function
      
      !----------------------------------------
      
      !----- HeI -----! 
      
      ! He I GROUND state (1^1S, singlet) photoionization.  Default: Verner+1996
      ! (VFKY96) single-shell fit.  ates_photoion_rate=.true. reverts to the legacy
      ! ATES two-term fit.  (The 2^3S metastable TRIPLET is handled separately by
      ! sigma_HeI23S; this routine is the singlet ground state only.)
      pure double precision function sigma_HeI(E)
      real*8,intent(in) :: E
      real*8 :: eth

      if (ates_photoion_rate) then
         ! --- Legacy ATES two-term fit ---
         ! Gated at the threshold itself, the same energy the default branch
         ! and the photon grid use.  (ATES turned this fit on 0.1 per cent
         ! below 24.6 eV so that a grid point sitting exactly on the threshold
         ! would not be dropped; the thresholds are bin EDGES now, so a
         ! sub-threshold turn-on only lets a photon that cannot ionize be
         ! absorbed.)
         eth = e_th_HeI
         if (E.ge.eth) then
            sigma_HeI = 0.6935d0/((E*1.0d-2)**1.82d0+(E*1.0d-2)**3.23d0)
         else
            sigma_HeI = 0.0d0
         endif
      else
         ! --- Default: Verner, Ferland, Korista & Yakovlev 1996 He I 1^1S ---
         ! [E_th,E_0,sigma_0,y_a,P,y_w,y_0,y_1] (sigma_0 in 1e-18 cm^2)
         ! E_th is passed as the global He I ionization threshold, not as the
         ! 24.59 eV of Verner et al.'s Table 1.  E_th only GATES the fit
         ! (sigma = 0 below it); it does not appear in the functional form,
         ! whose parameters E_0, sigma_0, y_a, P, y_w, y_0, y_1 are the
         ! published ones unchanged.  The shift is therefore confined to the
         ! 0.003 eV at the onset, and it removes the band that was integrated
         ! as zero because the grid's He I edge sat above this turn-on.
         sigma_HeI = sigma_VFKY96(E, e_th_HeI, 1.361d1, 9.492d2, 1.469d0, &
                                     3.188d0, 2.039d0, 4.434d-1, 2.136d0)
      endif

      ! End of function
      end function sigma_HeI
      
      !----------------------------------------
      
      ! HeI(23S) triplet photoionization cross section.
      ! Two Verner+1996 (VFKY96) single-shell forms joined by a power-law
      ! (log-linear) bridge, used piecewise:
      !   wing A = threshold -> Cooper minimum   (E <= x3 ~ 34.7 eV),
      !   bridge = log-linear rise to the bump   (x3 < E < x4 ~ 45.6 eV),
      !   wing B = bump -> high-E tail           (E >= x4).
      ! Each wing is a ~4% fit to the former Norcross-1971 broken PL / TOPbase
      ! background; the bridge joins them continuously.  Threshold (E<4.78 eV)
      ! and the high-E tail are handled by wing A/B themselves (no hard cutoff).
      ! See docs/photoion_cross_sections.tex (sec. "He 2^3S as two VFKY96
      ! forms") and docs/Update_EXHALE.  The former Norcross broken-PL +
      ! TOPbase-tail implementation is retained, commented out, below.
      !	Fit of data from Norcross (1971)
      ! Note: Oklopcic calculations include HeITR only in the
      !	range [4.8-13.6] eV
	
      double precision function sigma_HeI23S(E)
      real*8, intent(in) :: E
      real*8 :: logE, x3, x4, y3b, y4b, mb
      ! VFKY96 fit parameters [E_th,E_0,sigma_0,y_a,P,y_w,y_0,y_1] (Mb, eV):
      !   wing A -- former segments 1-2 (threshold -> Cooper min, ~4.8-34.7 eV)
      ! EthA is the global He 2^3S ionization threshold, not a fit parameter:
      ! as in sigma_HeI it only gates where wing A turns on, and it is the
      ! same energy that floors the photon grid when the metastable is
      ! carried, so the grid's lowest bin no longer starts below the
      ! turn-on.  E0A ... y1A are the fitted parameters and are unchanged.
      real*8, parameter :: E0A=2.645d0, s0A=2.08d1,  yaA=1.0d12,             &
                           PA=3.42d0,    ywA=2.681d0, y0A=1.956d0, y1A=2.603d0
      real*8, parameter :: EthA=e_th_HeTR
      !   wing B -- former segments 4-5 (bump -> tail, ~45.6-320 eV and beyond)
      real*8, parameter :: EthB=45.59d0, E0B=49.68d0, s0B=1.052d3, yaB=4.393d-2, &
                           PB=2.941d0,   ywB=1.717d0, y0B=5.488d-5,y1B=1.118d0

      ! Transition (bridge) node energies: Cooper minimum x3 and bump x4
      x3 = log10(hp_eV*c_light/(357.340d0*1d-8))   ! ~34.70 eV
      x4 = log10(hp_eV*c_light/(271.940d0*1d-8))    ! ~45.59 eV

      logE = log10(E)

      if (logE .le. x3) then
         ! Wing A: threshold -> Cooper minimum (VFKY96; =0 below E=4.78 eV)
         sigma_HeI23S = sigma_VFKY96(E, EthA,E0A,s0A,yaA,PA,ywA,y0A,y1A)
      else if (logE .ge. x4) then
         ! Wing B: resonance-averaged bump -> high-E tail (VFKY96)
         sigma_HeI23S = sigma_VFKY96(E, EthB,E0B,s0B,yaB,PB,ywB,y0B,y1B)
      else
         ! Bridge: log-linear (power-law) join of the two wings
         y3b = log10(sigma_VFKY96(10.0d0**x3, EthA,E0A,s0A,yaA,PA,ywA,y0A,y1A))
         y4b = log10(sigma_VFKY96(10.0d0**x4, EthB,E0B,s0B,yaB,PB,ywB,y0B,y1B))
         mb  = (y4b - y3b)/(x4 - x3)
         sigma_HeI23S = 10.0d0**(mb*(logE - x3) + y3b)
      endif

      ! ================================================================
      ! FORMER implementation (Norcross 1971 broken power law + re-aimed
      ! last segment + TOPbase tail), retained for reference:
      ! ----------------------------------------------------------------
      !   real*8 :: logsigma, a1,c1,a2,c2,a3,c3, x1,x2,xJ, y3,y4,m,
      !  &          a3n,c3n,ptail,ctail,tailJ
      !   x1 = log10(hp_eV*c_light/(2593.01*1e-8))*0.9999   ! ~4.78 eV
      !   x2 = log10(hp_eV*c_light/(1655.63*1e-8))          ! ~7.49 eV
      !   xJ = log10(70.0d0)                                ! 70 eV junction
      !   a1=-0.8134; a2=-1.772; a3=-3.039; c1=1.240
      !   c2 = c1 + x2*(a1-a2); c3 = 5.470                  ! continuity
      !   y3 = a2*x3 + c2; y4 = a3*x4 + c3; m = (y4-y3)/(x4-x3)
      !   ptail=-2.9964; ctail=5.4623; tailJ=ptail*xJ+ctail
      !   a3n=(tailJ-y4)/(xJ-x4); c3n=y4-a3n*x4
      !   if (logE .lt. x1) then; sigma_HeI23S=0.0; return; endif
      !   if (logE .le. x2)                    logsigma = a1*logE + c1
      !   if (logE .gt. x2 .and. logE .le. x3) logsigma = a2*logE + c2
      !   if (logE .gt. x3 .and. logE .lt. x4) logsigma = m*(logE-x3) + y3
      !   if (logE .ge. x4 .and. logE .le. xJ) logsigma = a3n*logE + c3n
      !   if (logE .gt. xJ)                    logsigma = ptail*logE + ctail
      !   sigma_HeI23S = 10.0**logsigma
      !   [pre-extension variant: hard cutoff at x5=log10(59.2 eV),
      !    last segment slope a3=-3.039]
      ! ================================================================

      ! End of function
      end function sigma_HeI23S
      
      !----------------------------------------

      !----- Metal photoionization (Verner et al. 1996) -----!
      !
      ! Single-shell analytic fit from Verner, Ferland, Korista &
      ! Yakovlev 1996, ApJ 465, 487, Eqs. (1)-(4). Returns sigma in
      ! units of 1e-18 cm^2 (consistent with rest of this module).
      !
      ! Two forms are provided:
      !
      !  sigma_VFKY96 -- the FULL published Eq. (1) with the (y_0,y_1)
      !  offset/asymptote parameters (photo.dat Table 1, phfit2.f ph2
      !  branch):
      !       x     = E/E_0 - y_0
      !       z     = sqrt(x^2 + y_1^2)
      !       sigma = sigma_0 * ((x-1)^2 + y_w^2) * z^(-Q)
      !                       * (1 + sqrt(z/y_a))^(-P)
      !       Q     = 5.5 - 0.5*P                     [no orbital-l term]
      !  photo.dat / ph2 columns: E_0, sigma_0, y_a, P, y_w, y_0, y_1.
      !
      !  sigma_Verner96 -- the SIMPLER Verner & Yakovlev 1995 (A&AS 109,
      !  125) inner-shell form, retained for ions whose published fit has
      !  y_0=y_1=0 (Mg I/II, Na I) and for the non-OP K I (phfit2 ph1):
      !       y     = E/E_0
      !       Q     = 5.5 + l - 0.5*P
      !       sigma = sigma_0 * ((y-1)^2 + y_w^2) * y^(-Q)
      !                       * (1 + sqrt(y/y_a))^(-P)
      !  The two agree exactly when y_0 = y_1 = 0 and l = 0.

      pure double precision function sigma_VFKY96(E,E_th,E_0,s_0,y_a,P,y_w, &
                                             y_0,y_1)
      real*8, intent(in) :: E,E_th,E_0,s_0,y_a,P,y_w,y_0,y_1
      real*8 :: x,z,Q,Fy

      if (E .lt. E_th) then
         sigma_VFKY96 = 0.0d0
         return
      endif

      x  = E/E_0 - y_0
      z  = sqrt(x*x + y_1*y_1)
      Q  = 5.5d0 - 0.5d0*P
      Fy = ((x-1.0d0)**2 + y_w**2) * z**(-Q)          &
                                   * (1.0d0 + sqrt(z/y_a))**(-P)
      sigma_VFKY96 = s_0 * Fy

      end function sigma_VFKY96

      !--------------!

      double precision function sigma_Verner96(E,E_th,E_0,s_0,y_a,P,y_w,l)
      real*8, intent(in) :: E,E_th,E_0,s_0,y_a,P,y_w
      integer, intent(in) :: l
      real*8 :: y,Q,Fy

      if (E .lt. E_th) then
         sigma_Verner96 = 0.0d0
         return
      endif

      y  = E/E_0
      Q  = 5.5d0 + dble(l) - 0.5d0*P
      Fy = ((y-1.0d0)**2 + y_w**2) * y**(-Q)            &
                                 * (1.0d0 + sqrt(y/y_a))**(-P)
      sigma_Verner96 = s_0 * Fy

      end function sigma_Verner96

      !--------------!

      ! H2 TOTAL PHOTOABSORPTION cross section.
      !
      ! WHAT THIS REPLACED, AND WHY (E1b).  It used to be the three-piece
      ! analytic fit of Yan, Sadeghpour & Dalgarno (1998), ApJ 496, 1044,
      ! Eqs. 17-19 throughout.  Those three pieces DO NOT MEET at the two
      ! energies where the paper joins them; evaluated on either side, in
      ! units of 1e-18 cm^2,
      !
      !     18 eV   Eq. 17 gives 5.693, Eq. 18 gives 9.787   (x 1.719)
      !     85 eV   Eq. 18 gives 0.1281, Eq. 19 gives 0.0793 (x 0.619)
      !
      ! Both formulas were re-derived by hand from the published
      ! coefficients, so this is a property of the fit as published and not
      ! a transcription error.  A 72 per cent step in an opacity is a
      ! physical defect: the JFNK residual differentiates it, and the H2
      ! destruction rate jumps across it.
      !
      ! WHAT IS USED NOW.  Over the range it was measured, the MEASUREMENT
      ! is used, taken from the paper that made it:
      !
      !   Samson, J. A. R., & Haddad, G. N. 1994, J. Opt. Soc. Am. B 11,
      !   277, "Total photoabsorption cross sections of H2 from 18 to
      !   113 eV", their TABLE 1, "Recommended Total Photoabsorption Cross
      !   Sections of H2 from 18 eV to 300 eV", 72 rows.
      !
      ! Their stated accuracy is +/-2% to +/-3% over 18-113 eV, the measured
      ! range, and +/-3% to +/-4% from 113 to 300 eV.  Yan Eq. 18 is not
      ! used at all: the table covers its whole range.
      !
      !   15.4 <= E < 18 eV : k_bx * Backx et al. (1976) Table 1
      !   18 <= E <= 300 eV : Samson & Haddad (1994) Table 1
      !   E > 300 eV        : k_hi * Yan Eq. 19     (E^-7/2 sum-rule tail)
      !
      ! all tables interpolated in log sigma against log E, and
      !   k_bx = 9.85/Backx(18)     = 1.016415
      !   k_hi = 0.00154/Eq19(300)  = 0.885805
      !
      ! NO PIECE OF THE YAN FIT SURVIVES BELOW 300 eV.
      !
      ! WHY THE TABLE IS THIS ONE AND NOT CHUNG'S.  Chung, Lee, Masuoka &
      ! Samson (1993) Table II carries a sigma(abs) column over
      ! 18.076-124 eV, and this code used to take it from there.  It is the
      ! SAME data: their text says "Values for a(abs) were taken from the
      ! data of Samson and Haddad" and their reference 31 is
      ! "J. A. R. Samson and G. N. Haddad, J. Opt. Soc. Am. (submitted)",
      ! the paper above.  Checked row by row at the 50 energies the two
      ! tables share, 49 agree exactly and one differs in its last digit
      ! (62 eV, 0.237 against 0.236, 0.42 per cent).  The original is used
      ! because it reaches 300 eV where Chung's stops at 124, so the
      ! high-energy join moves out of the band the wind actually absorbs
      ! in.  Chung's Table II is still the source of the H+/H2+ BRANCHING
      ! (frac_H2_dissociative_ionization) and of the neutral fraction
      ! (h2_photo_channels), both of which are ratios and are unaffected.
      !
      ! 15.4-18 eV: A TABLE, from the experiment the other two normalize to.
      !
      !   Backx, C., Wight, G. R., & Van der Wiel, M. J. 1976, J. Phys. B 9,
      !   315, "Oscillator strengths (10-70 eV) for absorption, ionization
      !   and dissociation in H2, HD and D2, obtained by an electron-ion
      !   coincidence method", their TABLE 1, column f^0(E).
      !
      ! UNITS AND NORMALIZATION, taken from the paper and not assumed.  Their
      ! Fig. 2 plots f^0(E) in eV^-1 with a x10^-1 axis multiplier and the
      ! table is in 10^-2 eV^-1; the absolute scale is fixed by the TRK sum
      ! rule, "normalized on an integral value of two (TRK sum rule) with a
      ! 1% correction for energies beyond 80 eV" (their Fig. 2 caption).  The
      ! conversion used here is the standard sigma[Mb] = 109.75 df/dE[1/eV].
      ! That constant is CONFIRMED against the other measurement rather than
      ! assumed: it puts Backx at 7.628 Mb at 20 eV and 2.700 at 30 eV
      ! against Samson & Haddad's 7.65 and 2.71, i.e. 0.3 and 0.4 per cent.
      ! (A factor 100 instead would put them 9 per cent low.)
      !
      ! THE JOIN.  Backx gives 9.691 Mb at 18 eV where Samson & Haddad give
      ! 9.850, so the branch is scaled by k_bx = 1.016415 to meet the table
      ! exactly.  THAT 1.6 PER CENT IS A MEASURED DIFFERENCE BETWEEN TWO
      ! INDEPENDENT EXPERIMENTS, not a fitted constant, and it is inside
      ! both quoted accuracies.
      !
      ! HOW BIG THE CORRECTION IS.  The branch this replaced was Yan Eq. 17
      ! scaled by a continuity constant.  Eq. 17 falls to 0.10 (0.17 after
      ! scaling) at 15.4 eV, where the measurement is 14.5: the analytic fit
      ! went to ZERO at the ionization threshold, which no photoionization
      ! cross section does.  The error ran from a factor 85 at threshold to
      ! 1.7 at 18 eV.
      !
      ! AN EARLIER VERSION OF THIS BRANCH was digitized from Fig. 1 of Chan,
      ! Cooper & Brion (1992), Chem. Phys. 168, 375, whose Tables 1-4 give
      ! only the discrete Lyman and Werner transitions.  That reading is kept
      ! in section 151.10 of docs/Update_EXHALE_stage1.pdf as a cross-check -- it
      ! agrees with the table used here to 2 per cent at 16 eV and 6 per cent
      ! at 17 eV -- and is not used in the code.  A table beats a figure.
      !
      ! WHAT IS EXCLUDED.  The table is evaluated only at and above the code
      ! threshold, 15.4 eV.  Its 15.0 eV row is carried so the interpolation
      ! into 15.4 has a left end; it is never returned.  Below the
      ! ionization potential the spectrum is the discrete Lyman and Werner
      ! band system, which this code carries separately as Lyman-Werner
      ! photodissociation and must not add to an ionization-continuum
      ! opacity twice.
      !
      ! THE SPREAD BETWEEN EXPERIMENTS IS THE REAL UNCERTAINTY.  Compared at
      ! the energies all three cover (section 151.10 has the table):
      ! Backx/Samson-Haddad runs 0.98 to 1.03 over 18-31 eV and 1.06 to 1.23
      ! over 32-68 eV; Lee, Carlson & Judge (1976), JQSRT 16, 873 runs 0.92
      ! to 0.96 over 18-31 eV and 1.0 to 1.5 above.  So sigma_H2 is good to
      ! a few per cent below about 30 eV and to 10-20 per cent above it,
      ! whatever any single paper's quoted accuracy says.
      !
      ! WHAT THE HIGH JOIN COSTS, and it is a real tension.  Eq. 19 agrees
      ! with the table to 0.33 per cent at 113 eV (0.033270 against
      ! 0.033160), so over the measured range the sum-rule tail and the
      ! measurement are the same curve.  They part company above it: at
      ! 300 eV the table reads 0.00154 and Eq. 19 gives 0.00174, and Yan et
      ! al. say why in their section 4 -- they RAISED the 300 eV point from
      ! 1.54e-21 to 1.75e-21 cm^2 precisely because the S(2) sum rules were
      ! not satisfied with the measured value.  So k_hi = 0.886 carries the
      ! measured normalization into the tail and gives up 12 per cent of
      ! their sum-rule correction.  Continuity is chosen over the sum rule
      ! here because a step in an opacity is what this whole change exists
      ! to remove, and because 300 eV is four decades below the peak of the
      ! cross section: sigma_H2(300) is 1.5e-3 against 9.85 at 18 eV.
      !
      ! Chung and Samson & Haddad tabulate megabarns, which are already
      ! 1e-18 cm^2; the Yan pieces are in barns and are converted
      ! (1 barn = 1e-6 Mb).
      double precision function sigma_H2(E)
      real*8, intent(in) :: E
      real*8 :: x, EkeV, sb
      integer :: k
      ! Turn-on of the H2 photoabsorption: the global H2 ionization threshold,
      ! the same energy the photoelectron h nu - e_th_H2 is charged against.
      ! It is NOT the 15.4 eV that scales the Yan, Sadeghpour & Dalgarno
      ! (1998) sum-rule tail below (E_I_yan), which is a parameter of their
      ! Eq. 19 and stays as published.
      real*8, parameter :: eth = e_th_H2
      integer, parameter :: n_ab = 72
      ! Samson & Haddad (1994) Table 1, photon energy [eV] ...
      real*8, parameter :: e_ab(n_ab) = (/                                &
              18.0000d0,     19.0000d0,     20.0000d0,     21.0000d0,     22.0000d0,       &
              23.0000d0,     24.0000d0,     25.0000d0,     26.0000d0,     27.0000d0,       &
              28.0000d0,     29.0000d0,     30.0000d0,     31.0000d0,     32.0000d0,       &
              33.0000d0,     34.0000d0,     35.0000d0,     36.0000d0,     37.0000d0,       &
              38.0000d0,     39.0000d0,     40.0000d0,     42.0000d0,     44.0000d0,       &
              46.0000d0,     48.0000d0,     50.0000d0,     52.0000d0,     54.0000d0,       &
              56.0000d0,     58.0000d0,     60.0000d0,     62.0000d0,     64.0000d0,       &
              66.0000d0,     68.0000d0,     70.0000d0,     72.0000d0,     74.0000d0,       &
              76.0000d0,     78.0000d0,     80.0000d0,     85.0000d0,     90.0000d0,       &
              95.0000d0,    100.0000d0,    105.0000d0,    110.0000d0,    115.0000d0,       &
             120.0000d0,    125.0000d0,    130.0000d0,    135.0000d0,    140.0000d0,       &
             145.0000d0,    150.0000d0,    160.0000d0,    170.0000d0,    180.0000d0,       &
             190.0000d0,    200.0000d0,    210.0000d0,    220.0000d0,    230.0000d0,       &
             240.0000d0,    250.0000d0,    260.0000d0,    270.0000d0,    280.0000d0,       &
             290.0000d0,    300.0000d0 /)
      ! ... and their recommended sigma [1e-18 cm^2].
      real*8, parameter :: a_ab(n_ab) = (/                                &
              9.850000d0,     8.680000d0,     7.650000d0,     6.770000d0,     6.030000d0,       &
              5.400000d0,     4.830000d0,     4.310000d0,     3.850000d0,     3.500000d0,       &
              3.200000d0,     2.950000d0,     2.710000d0,     2.410000d0,     2.130000d0,       &
              1.930000d0,     1.740000d0,     1.570000d0,     1.420000d0,     1.290000d0,       &
              1.170000d0,     1.080000d0,     0.990000d0,     0.840000d0,     0.722000d0,       &
              0.626000d0,     0.543000d0,     0.477000d0,     0.420000d0,     0.373000d0,       &
              0.331000d0,     0.294000d0,     0.262000d0,     0.236000d0,     0.212000d0,       &
              0.191000d0,     0.173000d0,     0.157000d0,     0.143000d0,     0.130000d0,       &
              0.118000d0,     0.108000d0,     0.098000d0,     0.079200d0,     0.065800d0,       &
              0.055600d0,     0.047900d0,     0.041300d0,     0.035700d0,     0.031600d0,       &
              0.027100d0,     0.023800d0,     0.021100d0,     0.018600d0,     0.016500d0,       &
              0.014800d0,     0.013300d0,     0.010800d0,     0.008950d0,     0.007550d0,       &
              0.006400d0,     0.005450d0,     0.004700d0,     0.004030d0,     0.003500d0,       &
              0.003080d0,     0.002700d0,     0.002400d0,     0.002130d0,     0.001900d0,       &
              0.001710d0,     0.001540d0 /)
      integer, parameter :: n_lo = 7
      ! Backx et al. (1976) Table 1, photon energy [eV] ...
      real*8, parameter :: e_lo(n_lo) = (/                                &
              15.0000d0,     15.5000d0,     16.0000d0,     16.5000d0,     17.0000d0,       &
              17.5000d0,     18.0000d0 /)
      ! ... and k_bx * 109.75 * f^0(E), their absorption column [1e-18 cm^2].
      real*8, parameter :: a_lo(n_lo) = (/                                &
              14.8364d0,     14.0555d0,     13.0515d0,     12.0476d0,     11.2667d0,       &
              10.6866d0,      9.8500d0 /)

      if (E .lt. eth) then
         sigma_H2 = 0.0d0
         return
      endif

      if (E .lt. e_ab(1)) then
         ! The Backx table, interpolated in log sigma against log E.
         k = 1
         do while (k .lt. n_lo-1 .and. E .gt. e_lo(k+1))
            k = k + 1
         enddo
         sigma_H2 = exp(log(a_lo(k))                                      &
              + (log(a_lo(k+1)) - log(a_lo(k)))                           &
                *(log(E) - log(e_lo(k)))/(log(e_lo(k+1)) - log(e_lo(k))))
         return
      endif

      if (E .le. e_ab(n_ab)) then
         ! The measurement, interpolated in log sigma against log E.  Both
         ! are positive and smooth over the whole table, and a power law
         ! between adjacent rows is what the cross section is between them.
         k = 1
         do while (k .lt. n_ab-1 .and. E .gt. e_ab(k+1))
            k = k + 1
         enddo
         sigma_H2 = exp(log(a_ab(k))                                      &
              + (log(a_ab(k+1)) - log(a_ab(k)))                           &
                *(log(E) - log(e_ab(k)))/(log(e_ab(k+1)) - log(e_ab(k))))
         return
      endif

      ! Sum-rule tail, Eq. 19, rescaled to meet the table at 300 eV.  Yan et
      ! al. scale the photon energy by the H2 ionization potential they quote,
      ! 15.4 eV; that number is a coefficient of their fit, so it is written
      ! as E_I_yan and not as the measured threshold e_th_H2.
      x    = E/E_I_yan
      EkeV = E*1.0d-3
      sb = 45.57d0*(1.0d0 - 2.003d0/sqrt(x) - 4.806d0/x                  &
                    + 50.577d0/x**1.5d0 - 171.044d0/(x*x)                &
                    + 231.608d0/x**2.5d0 - 81.885d0/x**3)                &
           / EkeV**3.5d0
      if (sb .lt. 0.0d0) sb = 0.0d0
      sigma_H2 = sb*1.0d-6*sigma_H2_k_hi()

      end function sigma_H2

      !--------------!

      ! The one surviving continuity constant of sigma_H2: the ratio of the
      ! measured table endpoint at 300 eV to Eq. 19 evaluated there.

      double precision function sigma_H2_k_hi()
      real*8 :: x, EkeV, sb
      real*8, parameter :: E0 = 300.0d0, a0 = 0.00154d0
      x    = E0/E_I_yan
      EkeV = E0*1.0d-3
      sb = 45.57d0*(1.0d0 - 2.003d0/sqrt(x) - 4.806d0/x                  &
                    + 50.577d0/x**1.5d0 - 171.044d0/(x*x)                &
                    + 231.608d0/x**2.5d0 - 81.885d0/x**3)                &
           / EkeV**3.5d0
      sigma_H2_k_hi = a0/(sb*1.0d-6)
      end function sigma_H2_k_hi

      !--------------!

      ! Fraction of an H2 PHOTOIONIZATION that is DISSOCIATIVE,
      !
      !     H2 + hv -> H + H+ + e-      (threshold 18.076 eV),
      !
      ! as opposed to the non-dissociative channel
      !
      !     H2 + hv -> H2+ + e-         (threshold 15.4 eV).
      !
      ! sigma_H2 above is the TOTAL photoionization cross section, so this is
      ! a BRANCHING of a cross section the code already carries and not an
      ! extra opacity: the rate at which photons destroy H2 does not change,
      ! only what the destruction leaves behind. The dissociative channel
      ! puts a proton and an H atom into the gas instead of an H2+ ion, so it
      ! bypasses the H3+ chain (H2+ + H2 -> H3+ + H) altogether.
      !
      ! SOURCE. Chung, Lee, Masuoka & Samson (1993), J. Chem. Phys. 99, 885,
      ! "Dissociative photoionization of H2 from 18 to 124 eV", their
      ! TABLE II, which tabulates sigma(H+) and sigma(H2+) at 69 energies
      ! from the threshold to 124 eV. The table below is
      !
      !     f_di(E) = sigma(H+) / [ sigma(H+) + sigma(H2+) ] ,
      !
      ! the share of the IONIZATIONS that release a proton, linearly
      ! interpolated in E. Their two partial columns reproduce their own
      ! ratio column to better than 0.5% at every row, and their sum
      ! reproduces their sigma(abs) column exactly except at the ten rows
      ! flagged with their footnote b (33-41 eV), where the photoionization
      ! yield is below unity and the difference is their Table I neutral
      ! cross section -- both checks were run on the transcription.
      !
      ! WHY NOT THE NUMBERS QUOTED BY YAN, SADEGHPOUR & DALGARNO (1998).
      ! Their section 4 says the dissociative "ratio to the total
      ! photoionization cross section increases from zero at the threshold of
      ! 18.08 eV to 0.284 at 76 eV and then decreases slowly to 0.258 at
      ! 124 eV". Those two numbers are Chung et al.'s H+/H2+ COLUMN, whose
      ! caption reads "the partial cross sections sigma(H+), sigma(H2+), and
      ! the ratio H+/H2+": they are a ratio to the H2+ partial cross section,
      ! not to the total. The ratio to the total is r/(1+r) -- 0.221 at
      ! 76 eV and 0.205 at 124 eV -- so taking the quoted numbers as a
      ! branching of sigma_H2 overstates the channel by a factor 1.26-1.28
      ! above 70 eV, and, because the true curve rises within 0.03 eV of
      ! threshold to a 0.02 plateau instead of ramping linearly, understates
      ! it by up to 20x below 20 eV. It also misses the Q1 Rydberg
      ! autoionization resonance, which takes f_di to 0.109 at 35 eV and back
      ! down to 0.091 at 37.5 eV.
      !
      ! MEASUREMENT UNCERTAINTY, as the source states it: "An estimate of the
      ! total uncertainty in the data is about +/-4%-5%"; the ionization
      ! yields used between 18 and 28 eV are unity "within an accuracy of
      ! about +/-5%". Against other determinations their H+ data lie 0%-5%
      ! below Backx, Wight & van der Wiel over 18-33 eV (26% at 70 eV), and
      ! within 5% of Kossmann et al. over 25-75 eV (11%-20% below them over
      ! 75-110 eV).
      !
      ! ABOVE 124 eV, the top of the measured range, f_di is held at its
      ! 124 eV value, 0.2048. The source quotes no higher-energy data. The
      ! curve is nearly flat over its last 50 eV (0.2219 at 74 eV to 0.2048
      ! at 124 eV), so holding it is a mild extrapolation, but it is an
      ! extrapolation.
      !
      ! TWO THINGS THIS BRANCHING CARRIES THAT IT DOES NOT NAME, both stated
      ! by the source and both confined to where sigma_H2 is already small.
      !  * DOUBLE IONIZATION. Chung et al.'s eqs. (3)-(4) note that "the
      !    branching ratio H+/(H+ + H2+) is understood to include both single
      !    and double ionization where appropriate", and their section on
      !    Fig. 4 says that "by 80 eV about 20% of sigma(H+) comes from
      !    double ionization. This percentage remains fairly constant towards
      !    higher energies." Because their sigma(H+) counts PROTONS and a
      !    double ionization yields two, using f_di as a share of the H2
      !    ionization rate gives the PROTON source rate correctly; what it
      !    mis-assigns is the co-product, since those events leave no H atom
      !    and no H2+. Above 80 eV that is about 4% of the H2 ionizations.
      !  * NEUTRAL DISSOCIATION over 33-41 eV. sigma_H2 above follows Samson
      !    & Haddad's ABSORPTION cross section there (it is 1.76 against
      !    Chung et al.'s sigma(abs) 1.74 and their ionization sum 1.699 at
      !    34 eV), so the H2 destruction rate the code applies in that band
      !    includes the neutral channel H2 + hv -> H + H, up to 7% of it at
      !    37.5 eV, and the 1 - f_di branch hands that share to H2+. This is
      !    a property of the cross-section fit, not of the branching.
      double precision function frac_H2_dissociative_ionization(E)
      real*8, intent(in) :: E
      integer :: k
      integer, parameter :: n_di = 69
      ! Chung et al. (1993) Table II: photon energy [eV] ...
      real*8, parameter :: e_di_tab(n_di) = (/                            &
          18.076d0,   18.100d0,   18.150d0,   18.200d0,   18.300d0,       &
          18.400d0,   18.500d0,   18.600d0,   18.800d0,   19.000d0,       &
          19.500d0,   20.000d0,   20.500d0,   21.000d0,   21.500d0,       &
          22.000d0,   23.000d0,   24.000d0,   25.000d0,   26.000d0,       &
          27.000d0,   28.000d0,   29.000d0,   30.000d0,   31.000d0,       &
          32.000d0,   33.000d0,   34.000d0,   35.000d0,   35.500d0,       &
          36.000d0,   36.500d0,   37.000d0,   37.500d0,   38.000d0,       &
          39.000d0,   40.000d0,   41.000d0,   42.000d0,   43.000d0,       &
          44.000d0,   45.000d0,   46.000d0,   48.000d0,   50.000d0,       &
          52.000d0,   54.000d0,   56.000d0,   58.000d0,   60.000d0,       &
          62.000d0,   64.000d0,   66.000d0,   68.000d0,   70.000d0,       &
          72.000d0,   74.000d0,   76.000d0,   78.000d0,   80.000d0,       &
          85.000d0,   90.000d0,   95.000d0,  100.000d0,  105.000d0,       &
         110.000d0,  115.000d0,  120.000d0,  124.000d0 /)
      ! ... and sigma(H+)/[sigma(H+) + sigma(H2+)] at those energies.
      real*8, parameter :: f_di_tab(n_di) = (/                            &
         0.00000d0,  0.00219d0,  0.00457d0,  0.00626d0,  0.00853d0,       &
         0.01009d0,  0.01133d0,  0.01268d0,  0.01461d0,  0.01613d0,       &
         0.01867d0,  0.01961d0,  0.02018d0,  0.02082d0,  0.02114d0,       &
         0.02156d0,  0.02222d0,  0.02298d0,  0.02365d0,  0.02454d0,       &
         0.02674d0,  0.03216d0,  0.03938d0,  0.04832d0,  0.06107d0,       &
         0.07512d0,  0.08759d0,  0.09947d0,  0.10845d0,  0.10940d0,       &
         0.10448d0,  0.09524d0,  0.09167d0,  0.09091d0,  0.09225d0,       &
         0.10049d0,  0.10892d0,  0.11641d0,  0.12381d0,  0.12949d0,       &
         0.13194d0,  0.13419d0,  0.13633d0,  0.14470d0,  0.15375d0,       &
         0.16370d0,  0.17203d0,  0.18152d0,  0.19020d0,  0.19817d0,       &
         0.20354d0,  0.20867d0,  0.21342d0,  0.21739d0,  0.22045d0,       &
         0.22160d0,  0.22188d0,  0.22119d0,  0.22037d0,  0.22041d0,       &
         0.21843d0,  0.21581d0,  0.21403d0,  0.21294d0,  0.21123d0,       &
         0.20964d0,  0.20786d0,  0.20635d0,  0.20477d0 /)

      if (E .le. e_di_tab(1)) then
         frac_H2_dissociative_ionization = 0.0d0
      else if (E .ge. e_di_tab(n_di)) then
         frac_H2_dissociative_ionization = f_di_tab(n_di)
      else
         k = 1
         do while (k .lt. n_di-1 .and. E .gt. e_di_tab(k+1))
            k = k + 1
         enddo
         frac_H2_dissociative_ionization = f_di_tab(k)                    &
              + (f_di_tab(k+1) - f_di_tab(k))                             &
                *(E - e_di_tab(k))/(e_di_tab(k+1) - e_di_tab(k))
      endif

      end function frac_H2_dissociative_ionization

      !--------------!
      !----- Metal photoionization: the sum over the shells -----!

      ! TOTAL photoabsorption cross section [1e-18 cm^2 = Mb] of the metal
      ! ion whose photo-table column index is k, i.e. species_table's
      ! mion_iphot(i): the outer shell plus every subshell that is open at
      ! this energy.  THE COLUMN ORDER IS DEFINED IN THE TABLE ABOVE, and
      ! set_energy_vectors fills sigma_tab through this function, so the
      ! ordering exists once.  It is needed off the e_v grid because the He
      ! recombination channels (util_ion_eq, he_rec_coupling) emit at their
      ! own photon energies.  Returns zero for k outside 1..n_metal_photo
      ! and below the ion's threshold.
      !
      ! The outer shell takes the fit that is valid at this energy: below
      ! E_max the 1996 Opacity-Project fit, at and above it the 1995 fit of
      ! the same shell.  This is the handover of Verner's own phfit2 (see
      ! the table above), and it is what keeps a fit from being
      ! extrapolated through the XUV: the Fe I row has Q = 5.5 - P/2 =
      ! 1.5385, so outside its range it RISES as E^0.46 and reaches 24.9 Mb
      ! at 1240 eV, where the true Fe I total is 0.51 Mb.  Each subshell
      ! takes its 1995 fit above its own turn-on.  Summed this way the
      ! function reproduces phfit2 shell for shell.
      !
      ! This is the absorption; what each shell's absorption does to the
      ! ion and to the photon energy is metal_shell_relaxation, shell by
      ! shell (metal_shell_photoion_sigma), and the head of this module.
      double precision function metal_photoion_sigma(k,E)
      integer, intent(in) :: k
      real*8,  intent(in) :: E

      integer :: is
      real*8  :: s

      metal_photoion_sigma = 0.0d0
      if (k .lt. 1 .or. k .gt. n_metal_photo) return
      if (E .lt. mph_e_th(k)) return

      if (E .ge. mph_e_max(k)) then
         s = sigma_Verner96(E, mph_e_th(k),                               &
              mph95_E0(k), mph95_s0(k), mph95_ya(k), mph95_P(k),          &
              mph95_yw(k), mph_l(k))
      else
         s = sigma_VFKY96(E, mph_e_th(k),                                 &
              mph96_E0(k), mph96_s0(k), mph96_ya(k), mph96_P(k),          &
              mph96_yw(k), mph96_y0(k), mph96_y1(k))
      endif

      do is = 1,mph_n_sub(k)
         if (E .lt. mph_sub_e_on(is,k)) cycle
         s = s + sigma_Verner96(E, mph_sub_e_on(is,k),                    &
              mph_sub_E0(is,k), mph_sub_s0(is,k), mph_sub_ya(is,k),       &
              mph_sub_P(is,k), mph_sub_yw(is,k), mph_sub_l(is,k))
      enddo

      metal_photoion_sigma = s

      end function metal_photoion_sigma

      !--------------!

      ! The energies at which the photoabsorption of metal photo-table
      ! column k turns on above its own ionization threshold: the K and L
      ! edges of the subshells, and, for a valence subshell that the 1996
      ! outer-shell fit already represents below E_max, the E_max at which
      ! phfit2 lets it appear.  A step of the cross section has to be a bin
      ! EDGE of the photon grid (set_energy_vectors), so this is the list
      ! that routine adds for an active ion.  Ascending order is not
      ! guaranteed and is not needed: band_sub_boundaries sorts.
      subroutine metal_subshell_turn_on(k, e_on, n_on)
      integer, intent(in)  :: k
      real*8,  intent(out) :: e_on(:)
      integer, intent(out) :: n_on

      integer :: is

      n_on = 0
      if (k .lt. 1 .or. k .gt. n_metal_photo) return
      do is = 1,mph_n_sub(k)
         n_on       = n_on + 1
         e_on(n_on) = mph_sub_e_on(is,k)
      enddo

      end subroutine metal_subshell_turn_on

      !--------------!

      ! PARTIAL photoabsorption cross section [Mb] of ONE shell of metal
      ! photo-table column k: is = 0 is the outer shell, is = 1 ...
      ! mph_n_sub(k) the subshells of the mph_sub_* slots.  The same fits,
      ! turn-ons and handover as metal_photoion_sigma, which is their sum
      ! in the order is = 0, 1, 2, ...; summed in that order the shells
      ! reproduce it to the last bit.  Zero outside 1..n_metal_photo, for
      ! a slot the ion has no electrons in, and below the shell's turn-on.
      double precision function metal_shell_photoion_sigma(k, is, E)
      integer, intent(in) :: k, is
      real*8,  intent(in) :: E

      metal_shell_photoion_sigma = 0.0d0
      if (k .lt. 1 .or. k .gt. n_metal_photo) return
      if (E .lt. mph_e_th(k)) return
      if (is .eq. 0) then
         if (E .ge. mph_e_max(k)) then
            metal_shell_photoion_sigma = sigma_Verner96(E, mph_e_th(k),  &
                 mph95_E0(k), mph95_s0(k), mph95_ya(k), mph95_P(k),       &
                 mph95_yw(k), mph_l(k))
         else
            metal_shell_photoion_sigma = sigma_VFKY96(E, mph_e_th(k),    &
                 mph96_E0(k), mph96_s0(k), mph96_ya(k), mph96_P(k),       &
                 mph96_yw(k), mph96_y0(k), mph96_y1(k))
         endif
      else if (is .le. mph_n_sub(k)) then
         if (E .lt. mph_sub_e_on(is,k)) return
         metal_shell_photoion_sigma = sigma_Verner96(E,                  &
              mph_sub_e_on(is,k), mph_sub_E0(is,k), mph_sub_s0(is,k),     &
              mph_sub_ya(is,k), mph_sub_P(is,k), mph_sub_yw(is,k),        &
              mph_sub_l(is,k))
      endif

      end function metal_shell_photoion_sigma

      !--------------!

      ! What ONE absorption in shell is (0 = outer) of metal photo-table
      ! column k hands on, beyond the photon's own energy [eV]:
      !   e_th     the threshold the photoelectron is charged against, so
      !            the photoelectron carries h nu - e_th;
      !   e_aug    the mean kinetic energy of the Auger electrons;
      !   p_multi  the probability that two or more electrons leave;
      !   e_fluor  the mean energy that leaves as photons;
      !   e_ion_multi  the mean ionization energy of the electrons beyond
      !            the first.
      ! The outer shell ejects one electron and nothing else.  The five
      ! close the vacancy exactly, e_th = I + e_ion_multi + e_aug + e_fluor
      ! with I = mph_e_th(k) (the table block above, and the head of this
      ! module for what each part becomes).
      pure subroutine metal_shell_relaxation(k, is, e_th, e_aug, p_multi, &
                                             e_fluor, e_ion_multi)
      integer, intent(in)  :: k, is
      real*8,  intent(out) :: e_th, e_aug, p_multi, e_fluor, e_ion_multi

      e_th        = 0.0d0
      e_aug       = 0.0d0
      p_multi     = 0.0d0
      e_fluor     = 0.0d0
      e_ion_multi = 0.0d0
      if (k .lt. 1 .or. k .gt. n_metal_photo) return
      if (is .eq. 0) then
         e_th = mph_e_th(k)
      else if (is .le. mph_n_sub(k)) then
         ! The closure is imposed here in the code's own arithmetic, so the
         ! photon energy is shared out to the last bit: the Auger energy of
         ! an autoionizing vacancy is the remainder of the tabulated parts
         ! (mph_sub_e_auger is that remainder to the printed digits), and a
         ! radiative vacancy radiates exactly E_th,s - I.
         e_th        = mph_sub_e_th(is,k)
         p_multi     = mph_sub_p_multi(is,k)
         e_ion_multi = mph_sub_e_ion_multi(is,k)
         if (p_multi .gt. 0.0d0) then
            e_fluor = mph_sub_e_fluor(is,k)
            e_aug   = e_th - mph_e_th(k) - e_ion_multi - e_fluor
         else
            e_fluor = e_th - mph_e_th(k)
         endif
      endif

      end subroutine metal_shell_relaxation

      !----------------------------------------

      ! End of module
      end module Cross_sections
