      program metal_photoionization_fits
      ! Every metal photoionization cross section the code carries, against
      ! the published fits evaluated here from the published parameters,
      ! and each fit only where it is valid.
      !
      ! Production routine exercised: metal_photoion_sigma of
      ! src/modules/functions/cross_sec.f90, which is the one place the
      ! photo-table column order is defined and the one entry through which
      ! set_energy_vectors fills sigma_tab.
      !
      ! REFERENCES, and which ion takes which where.
      !
      !   Verner, D. A., Ferland, G. J., Korista, K. T., & Yakovlev, D. G.
      !   1996, ApJ 465, 487, "Atomic Data for Astrophysics. II. New
      !   Analytic FITS for Photoionization Cross Sections of Atoms and
      !   Ions", their Table 1 and Eqs. (1)-(4):
      !
      !       x     = E/E_0 - y_0
      !       y     = sqrt(x^2 + y_1^2)
      !       sigma = sigma_0 F(y),
      !       F(y)  = ((x-1)^2 + y_w^2) y^(-Q) (1 + sqrt(y/y_a))^(-P),
      !       Q     = 5.5 + l - 0.5 P    with l = 0 for their Table 1,
      !
      !   an OUTER-SHELL fit, on E_th <= E < E_max.  Their Table 1 covers
      !   the ions of H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, S,
      !   Ar, Ca and Fe, so it has no row for potassium.
      !
      !   Verner, D. A., & Yakovlev, D. G. 1995, A&AS 109, 125, "Analytic
      !   FITS for partial photoionization cross sections", their Table 1
      !   and Eq. (1), the same expression with y_0 = 0, y_1 = 0, i.e.
      !   y = E/E_0, and with the orbital quantum number l of the shell
      !   kept in Q: the fit of the SAME outer shell at E >= E_max, and
      !   the fit at every energy for K I, whose E_max is zero.
      !
      ! WHAT E_max IS.  It is the threshold of the outermost inner shell,
      ! i.e. the energy at which that shell opens and the outer-shell fit
      ! of 1996 leaves the range it was made in.  Verner's own routine
      ! phfit2 switches the outer shell there (its variable einn) and the
      ! code does the same, so the assertions below evaluate the 1996 form
      ! under E_max and the 1995 form over it, and check the switch itself
      ! at 1e-6 eV on either side.  phfit2 sets einn = 0 for potassium,
      ! which is why K I takes the 1995 fit everywhere.
      !
      ! The parameter rows below are the published ones, digit for digit,
      ! read from the author's own distribution of the two tables
      ! (references/verner_photo/photo.dat for the 1996 table,
      ! references/verner_photo/table1.dat for the 1995 one), which agree
      ! with the printed Table 1 of the 1996 paper row by row.  They are
      ! typed here independently of the table installed in cross_sec.f90,
      ! so this driver is a second transcription and not a copy of the
      ! first; metal_photoion_table_transcription.py separately compares
      ! the installed table with the distribution files field by field.
      !
      ! WHAT IS ASSERTED.  Each row is a transcription, not a fit of ours,
      ! so the assertion is that the production routine returns the
      ! published expression evaluated with the published constants: the
      ! tolerance is 1e-10 relative, round-off of the double-precision
      ! evaluation and nothing else.  Energies per ion: 1.05, 2 and 8 times
      ! the threshold, so that the onset rise ((x-1)^2 + y_w^2), the
      ! y^(-Q) fall and the (1 + sqrt(y/y_a))^(-P) tail are each exercised;
      ! 1e-6 eV under and over E_max; and 120, 350 and 1240 eV, the last
      ! being the top of the photon grid.
      !
      ! THE SHELL SUM.  metal_photoion_sigma returns the TOTAL
      ! photoabsorption of the ion, the outer shell plus every subshell
      ! open at that energy, which is what a caller of phfit2 obtains by
      ! summing that routine over is = 1 ... 7.  The published expression
      ! this driver compares against is built the same way from the rows
      ! below.  The turn-ons themselves are asserted 1e-6 eV on either
      ! side: under a K or L edge the shell must contribute nothing, over
      ! it exactly its 1995 fit.
      !
      ! WHAT THIS DOES NOT ASSERT.  Nothing here is a statement about what
      ! the absorption is charged to.  The model advances the ion one stage
      ! and hands the degradation cascade one electron of h nu - E_th(outer)
      ! whichever shell absorbed; that under-count, and the neglected
      ! fluorescence, are stated at the table in cross_sec.f90.

      use global_parameters, only: e_th_MgI, e_th_MgII
      use species_table,     only: n_mion, mion_iphot, mion_ethr,         &
                                   mion_name
      use Cross_sections,    only: metal_photoion_sigma
      use assertion_report

      implicit none

      integer, parameter :: n_ion = 17

      character(len=5), parameter :: ion_name(n_ion) = &
           [ 'CI   ','CII  ','OI   ','OII  ','NI   ','NII  ',              &
             'MgI  ','MgII ','SiI  ','SiII ','CaI  ','CaII ',              &
             'NaI  ','KI   ','SI   ','FeI  ','FeII ' ]

      ! The threshold that gates the fit.  Mg I, Mg II and Fe II take the
      ! NIST ionization potential, which is the energy the photon grid uses
      ! as their band edge, where the fit tables print 7.646, 15.04 and
      ! 16.19; a threshold gates the fit and does not enter it.
      real*8, parameter :: E_th(n_ion) = [                                 &
           1.126d1, 2.438d1, 1.362d1, 3.512d1, 1.453d1, 2.960d1,           &
           e_th_MgI, e_th_MgII, 8.152d0, 1.635d1, 6.113d0, 1.187d1,        &
           5.139d0, 4.341d0, 1.036d1, 7.902d0, 16.199d0 ]

      ! Verner et al. (1996) Table 1, E_max: the inner-shell threshold at
      ! which the outer-shell fit is handed over.  Zero for K I.
      real*8, parameter :: E_max(n_ion) = [                                &
           2.910d2, 3.076d2, 5.380d2, 5.581d2, 4.048d2, 4.236d2,           &
           5.490d1, 6.569d1, 1.060d2, 1.186d2, 3.443d1, 4.090d1,           &
           3.814d1, 0.0d0,   1.700d2, 6.600d1, 7.617d1 ]

      ! Verner et al. (1996) Table 1, [E_0, sigma_0, y_a, P, y_w, y_0, y_1]
      ! with E_0 in eV and sigma_0 in Mb, in the column order of
      ! metal_photoion_sigma.  K I (column 14) has no 1996 row.
      real*8, parameter :: E_0(n_ion) = [                                  &
           2.144d0, 4.058d-1, 1.240d0, 1.386d0, 4.034d0, 6.128d-2,         &
           1.197d1, 8.139d0, 2.317d1, 2.556d0, 1.278d1, 1.553d1,           &
           6.139d0, 0.0d0, 1.808d1, 5.461d-2, 1.761d-1 ]
      real*8, parameter :: sigma_0(n_ion) = [                              &
           5.027d2, 8.709d0, 1.745d3, 5.967d1, 8.235d2, 1.944d0,           &
           1.372d8, 3.278d0, 2.506d1, 4.140d0, 5.370d5, 1.064d7,           &
           1.601d0, 0.0d0, 4.564d4, 3.062d-1, 4.365d3 ]
      real*8, parameter :: y_a(n_ion) = [                                  &
           6.216d1, 1.261d2, 3.784d0, 3.175d1, 8.033d1, 8.163d2,           &
           2.228d-1, 4.341d7, 2.057d1, 1.337d1, 3.162d-1, 7.790d-1,        &
           6.148d3, 0.0d0, 1.000d0, 2.671d7, 6.298d3 ]
      real*8, parameter :: P_fit(n_ion) = [                                &
           5.101d0, 8.578d0, 1.764d1, 8.943d0, 3.928d0, 8.773d0,           &
           1.574d1, 3.610d0, 3.546d0, 1.191d1, 1.242d1, 2.130d1,           &
           3.839d0, 0.0d0, 1.361d1, 7.923d0, 5.204d0 ]
      real*8, parameter :: y_w(n_ion) = [                                  &
           9.157d-2, 2.093d0, 7.589d-2, 1.934d-2, 9.097d-2, 1.043d1,       &
           2.805d-1, 0.0d0, 2.837d-1, 1.570d0, 4.477d-1, 6.453d-1,         &
           0.0d0, 0.0d0, 6.385d-1, 2.069d1, 1.141d1 ]
      real*8, parameter :: y_0(n_ion) = [                                  &
           1.133d0, 4.929d1, 8.698d0, 2.131d1, 8.598d-1, 4.280d2,          &
           0.0d0, 0.0d0, 1.672d-5, 6.634d0, 1.012d-3, 2.161d-3,            &
           0.0d0, 0.0d0, 9.935d-1, 1.382d2, 9.272d1 ]
      real*8, parameter :: y_1(n_ion) = [                                  &
           1.607d0, 3.234d0, 1.271d-1, 1.503d-2, 2.325d0, 2.030d1,         &
           0.0d0, 0.0d0, 4.207d-1, 1.272d-1, 1.851d-2, 6.706d-2,           &
           0.0d0, 0.0d0, 2.486d-1, 2.481d-1, 1.075d2 ]

      ! Verner & Yakovlev (1995) Table 1, [E_0, sigma_0, y_a, P, y_w] of
      ! the SAME outer shell (2p for C, N, O; 3s for Na, Mg; 3p for Si, S;
      ! 4s for K, Ca, Fe), and the orbital quantum number l of that shell.
      real*8, parameter :: E_0_vy95(n_ion) = [                             &
           9.435d0, 1.094d1, 1.391d1, 1.745d1, 1.164d1, 1.827d1,           &
           9.393d0, 8.139d0, 2.212d1, 2.123d1, 7.366d0, 4.155d0,           &
           5.968d0, 3.824d0, 2.975d1, 1.277d1, 1.014d1 ]
      real*8, parameter :: sigma_0_vy95(n_ion) = [                         &
           1.152d3, 1.792d2, 1.220d5, 5.186d2, 1.029d4, 1.724d2,           &
           3.034d0, 3.278d0, 1.845d2, 6.975d1, 2.373d0, 2.235d0,           &
           1.460d0, 7.363d-1, 5.644d1, 1.468d0, 1.084d0 ]
      real*8, parameter :: y_a_vy95(n_ion) = [                             &
           5.687d0, 3.308d1, 1.364d0, 1.728d1, 2.361d0, 8.893d1,           &
           2.625d7, 4.341d7, 3.849d0, 4.907d0, 2.082d2, 1.595d4,           &
           2.557d7, 2.410d7, 1.321d1, 1.116d5, 2.562d4 ]
      real*8, parameter :: P_vy95(n_ion) = [                               &
           6.336d0, 4.150d0, 1.140d1, 4.995d0, 8.821d0, 3.348d0,           &
           3.923d0, 3.610d0, 9.721d0, 9.525d0, 4.841d0, 4.313d0,           &
           3.789d0, 4.427d0, 7.513d0, 4.112d0, 4.167d0 ]
      real*8, parameter :: y_w_vy95(n_ion) = [                             &
           4.474d-1, 5.276d-1, 4.103d-1, 2.182d-2, 4.239d-1, 4.209d-1,     &
           0.0d0, 0.0d0, 2.921d-1, 3.169d-1, 5.841d-4, 3.539d-1,           &
           0.0d0, 2.049d-4, 2.621d-1, 3.238d-2, 1.598d-2 ]
      integer, parameter :: l_shell(n_ion) = [ 1,1,1,1,1,1,0,0,1,1,        &
                                               0,0,0,0,1,0,0 ]

      ! The SUBSHELLS below the outer one, from phfit2's own PH1 DATA
      ! statements: their 1995 five parameters, their orbital quantum
      ! number, and the energy at which each turns on.  A shell with no
      ! electrons in the ground state is not listed (phfit2 carries a
      ! sigma_0 = 0 placeholder row for it and returns zero).
      !
      ! THE TURN-ON is max(shell threshold, E_max), the one rule phfit2's
      ! two guards amount to: a true inner shell (is <= nint) opens at its
      ! own threshold, while a valence subshell between nint and nout is
      ! part of what the 1996 outer-shell fit already represents and so is
      ! held at zero until E_max.  For C, N and O the 2s shell is such a
      ! subshell, which is why its turn-on is the K edge and not 19-46 eV.
      !
      ! Five thresholds here are NOT the ones table1.dat prints: phfit2's
      ! header says the inner-shell ionization energies of some low-ionized
      ! species were "slightly improved to fit smoothly the experimental
      ! inner-shell ionization energies of neutral atoms", and PH1 carries
      ! the improved values for Mg II 1s, Si II 1s, Ca II 1s and 2s, and
      ! Fe II 1s and 2s.  phfit2 is what the code reproduces, so PH1 is the
      ! reference here.
      integer, parameter :: n_sub_max = 6
      integer, parameter :: n_sub(n_ion) = [ &
           2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 3, 5, 4, 6, 6 ]
      real*8, parameter :: E_on_sub(n_sub_max,n_ion) = reshape( [ &
           2.910d+02 , 2.910d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           3.076d+02 , 3.076d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           5.380d+02 , 5.380d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           5.581d+02 , 5.581d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           4.048d+02 , 4.048d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           4.236d+02 , 4.236d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           1.311d+03 , 9.400d+01 , 5.490d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           1.320d+03 , 9.881d+01 , 6.569d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           1.846d+03 , 1.560d+02 , 1.060d+02 , 1.060d+02 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           1.848d+03 , 1.619d+02 , 1.186d+02 , 1.186d+02 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           4.043d+03 , 4.425d+02 , 3.523d+02 , 4.830d+01 , 3.443d+01 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           4.047d+03 , 4.445d+02 , 3.638d+02 , 6.037d+01 , 4.090d+01 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           1.079d+03 , 7.084d+01 , 3.814d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           3.614d+03 , 3.843d+02 , 3.014d+02 , 4.080d+01 , 2.466d+01 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           2.477d+03 , 2.350d+02 , 1.700d+02 , 1.700d+02 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           7.124d+03 , 8.570d+02 , 7.240d+02 , 1.040d+02 , 6.600d+01 , 6.600d+01, &  ! FeI   1s 2s 2p 3s 3p 3d
           7.140d+03 , 8.608d+02 , 7.341d+02 , 1.102d+02 , 7.617d+01 , 7.617d+01 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      real*8, parameter :: E_0_sub(n_sub_max,n_ion) = reshape( [ &
           8.655d+01 , 1.026d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           9.113d+01 , 2.991d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           1.774d+02 , 1.994d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           1.690d+02 , 1.759d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           1.270d+02 , 1.482d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           1.242d+02 , 1.094d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           2.711d+02 , 4.587d+01 , 4.937d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           3.709d+02 , 5.142d+01 , 4.940d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           5.322d+02 , 7.017d+01 , 7.808d+01 , 1.413d+01 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           4.580d+02 , 6.738d+01 , 7.154d+01 , 1.364d+01 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           6.947d+02 , 1.201d+02 , 1.529d+02 , 3.012d+01 , 4.487d+01 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           7.997d+02 , 1.335d+02 , 1.448d+02 , 3.176d+01 , 4.498d+01 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           4.216d+02 , 4.537d+01 , 3.655d+01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           1.171d+03 , 1.602d+02 , 2.666d+02 , 2.910d+01 , 4.138d+01 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           8.114d+02 , 1.047d+02 , 9.152d+01 , 1.916d+01 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           8.044d+02 , 5.727d+01 , 2.948d+02 , 4.334d+01 , 7.630d+01 , 1.407d+01, &  ! FeI   1s 2s 2p 3s 3p 3d
           8.931d+02 , 1.431d+02 , 2.281d+02 , 4.663d+01 , 7.750d+01 , 1.933d+01 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      real*8, parameter :: sigma_0_sub(n_sub_max,n_ion) = reshape( [ &
           7.421d+01 , 4.564d+03 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           6.649d+01 , 1.184d+03 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           3.237d+01 , 2.415d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           3.584d+01 , 1.962d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           4.748d+01 , 7.722d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           5.002d+01 , 7.483d+02 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           3.561d+01 , 1.671d+01 , 2.023d+02 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           1.767d+01 , 1.290d+01 , 2.049d+02 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           1.184d+01 , 1.166d+01 , 1.532d+02 , 1.166d+01 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           1.628d+01 , 1.263d+01 , 1.832d+02 , 8.454d+00 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           1.586d+01 , 1.010d+01 , 1.282d+02 , 7.227d+00 , 9.017d+01 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           1.168d+01 , 8.939d+00 , 1.446d+02 , 6.924d+00 , 7.314d+01 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           1.119d+01 , 1.142d+01 , 2.486d+02 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           4.540d+00 , 6.389d+00 , 3.107d+01 , 6.377d+00 , 2.614d+01 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           6.649d+00 , 8.520d+00 , 1.883d+02 , 1.003d+01 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           2.055d+01 , 1.076d+01 , 7.191d+01 , 5.921d+00 , 6.298d+01 , 1.850d+04, &  ! FeI   1s 2s 2p 3s 3p 3d
           1.666d+01 , 9.289d+00 , 1.241d+02 , 5.434d+00 , 4.624d+01 , 3.679d+04 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      real*8, parameter :: y_a_sub(n_sub_max,n_ion) = reshape( [ &
           5.498d+01 , 1.568d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           9.609d+01 , 3.085d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           3.812d+02 , 3.241d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           1.894d+02 , 4.020d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           1.380d+02 , 2.306d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           9.100d+01 , 2.793d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           2.374d+01 , 2.389d+01 , 1.079d+04 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           1.346d+02 , 5.775d+01 , 4.112d+03 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           2.580d+02 , 4.742d+01 , 5.765d+06 , 2.288d+01 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           7.706d+01 , 3.623d+01 , 3.537d+02 , 7.489d+01 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           2.563d+01 , 2.468d+01 , 2.217d+02 , 1.736d+02 , 1.465d+01 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           3.233d+01 , 2.914d+01 , 1.349d+02 , 5.246d+02 , 1.898d+01 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           5.642d+07 , 2.395d+02 , 3.222d+02 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           6.165d+03 , 1.044d+02 , 7.187d+06 , 2.229d+03 , 2.143d+02 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           3.734d+03 , 9.469d+01 , 7.193d+01 , 3.296d+01 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           3.633d+01 , 2.785d+01 , 3.219d+02 , 5.293d+01 , 1.479d+01 , 4.458d+00, &  ! FeI   1s 2s 2p 3s 3p 3d
           2.381d+01 , 2.026d+01 , 8.058d+01 , 9.271d+01 , 2.155d+01 , 3.564d+00 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      real*8, parameter :: P_sub(n_sub_max,n_ion) = reshape( [ &
           1.503d+00 , 1.085d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           1.338d+00 , 1.480d+01 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           1.083d+00 , 8.037d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           1.185d+00 , 7.999d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           1.252d+00 , 9.139d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           1.335d+00 , 9.956d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           1.952d+00 , 4.742d+00 , 2.960d+00 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           1.225d+00 , 3.903d+00 , 2.995d+00 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           1.102d+00 , 3.933d+00 , 2.639d+00 , 5.334d+00 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           1.385d+00 , 4.172d+00 , 3.133d+00 , 4.676d+00 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           1.966d+00 , 4.592d+00 , 3.087d+00 , 4.165d+00 , 7.498d+00 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           1.744d+00 , 4.269d+00 , 3.300d+00 , 3.771d+00 , 7.152d+00 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           7.736d-01 , 3.380d+00 , 3.570d+00 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           8.392d-01 , 3.159d+00 , 2.067d+00 , 3.587d+00 , 5.631d+00 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           8.646d-01 , 3.346d+00 , 3.633d+00 , 5.038d+00 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           2.118d+00 , 6.635d+00 , 2.837d+00 , 5.129d+00 , 7.672d+00 , 1.691d+01, &  ! FeI   1s 2s 2p 3s 3p 3d
           2.262d+00 , 5.376d+00 , 3.572d+00 , 4.640d+00 , 7.138d+00 , 1.566d+01 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      real*8, parameter :: y_w_sub(n_sub_max,n_ion) = reshape( [ &
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CI    1s 2s
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! CII   1s 2s
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OI    1s 2s
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! OII   1s 2s
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NI    1s 2s
           0.000d+00 , 0.000d+00 , 0.0d0     , 0.0d0     , 0.0d0     , 0.0d0, &  ! NII   1s 2s
           0.000d+00 , 0.000d+00 , 1.463d-02 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgI   1s 2s 2p
           0.000d+00 , 0.000d+00 , 2.223d-05 , 0.0d0     , 0.0d0     , 0.0d0, &  ! MgII  1s 2s 2p
           0.000d+00 , 0.000d+00 , 2.774d-04 , 0.000d+00 , 0.0d0     , 0.0d0, &  ! SiI   1s 2s 2p 3s
           0.000d+00 , 0.000d+00 , 2.870d-04 , 0.000d+00 , 0.0d0     , 0.0d0, &  ! SiII  1s 2s 2p 3s
           0.000d+00 , 0.000d+00 , 3.343d-03 , 0.000d+00 , 2.754d-01 , 0.0d0, &  ! CaI   1s 2s 2p 3s 3p
           0.000d+00 , 0.000d+00 , 3.358d-04 , 0.000d+00 , 2.735d-01 , 0.0d0, &  ! CaII  1s 2s 2p 3s 3p
           0.000d+00 , 0.000d+00 , 1.465d-01 , 0.0d0     , 0.0d0     , 0.0d0, &  ! NaI   1s 2s 2p
           0.000d+00 , 0.000d+00 , 5.274d-01 , 0.000d+00 , 2.437d-01 , 0.0d0, &  ! KI    1s 2s 2p 3s 3p
           0.000d+00 , 0.000d+00 , 2.485d-01 , 0.000d+00 , 0.0d0     , 0.0d0, &  ! SI    1s 2s 2p 3s
           0.000d+00 , 0.000d+00 , 6.314d-02 , 0.000d+00 , 2.646d-01 , 4.039d-01, &  ! FeI   1s 2s 2p 3s 3p 3d
           0.000d+00 , 0.000d+00 , 1.223d-01 , 0.000d+00 , 2.599d-01 , 4.144d-01 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )
      integer, parameter :: l_sub(n_sub_max,n_ion) = reshape( [ &
           0         , 0         , 0         , 0         , 0         , 0, &  ! CI    1s 2s
           0         , 0         , 0         , 0         , 0         , 0, &  ! CII   1s 2s
           0         , 0         , 0         , 0         , 0         , 0, &  ! OI    1s 2s
           0         , 0         , 0         , 0         , 0         , 0, &  ! OII   1s 2s
           0         , 0         , 0         , 0         , 0         , 0, &  ! NI    1s 2s
           0         , 0         , 0         , 0         , 0         , 0, &  ! NII   1s 2s
           0         , 0         , 1         , 0         , 0         , 0, &  ! MgI   1s 2s 2p
           0         , 0         , 1         , 0         , 0         , 0, &  ! MgII  1s 2s 2p
           0         , 0         , 1         , 0         , 0         , 0, &  ! SiI   1s 2s 2p 3s
           0         , 0         , 1         , 0         , 0         , 0, &  ! SiII  1s 2s 2p 3s
           0         , 0         , 1         , 0         , 1         , 0, &  ! CaI   1s 2s 2p 3s 3p
           0         , 0         , 1         , 0         , 1         , 0, &  ! CaII  1s 2s 2p 3s 3p
           0         , 0         , 1         , 0         , 0         , 0, &  ! NaI   1s 2s 2p
           0         , 0         , 1         , 0         , 1         , 0, &  ! KI    1s 2s 2p 3s 3p
           0         , 0         , 1         , 0         , 0         , 0, &  ! SI    1s 2s 2p 3s
           0         , 0         , 1         , 0         , 1         , 2, &  ! FeI   1s 2s 2p 3s 3p 3d
           0         , 0         , 1         , 0         , 1         , 2 ], &  ! FeII  1s 2s 2p 3s 3p 3d
           [n_sub_max,n_ion] )

      ! Three energies as multiples of the threshold, and three fixed ones
      ! in the XUV, the last at the top of the photon grid.
      real*8, parameter :: e_scale(3) = [ 1.05d0, 2.0d0, 8.0d0 ]
      real*8, parameter :: e_fixed(3) = [ 1.20d2, 3.50d2, 1.240d3 ]

      integer :: i, k, m
      real*8  :: E, measured, reference, s96, s95
      character(len=48) :: label

      do i = 1,n_ion
         do k = 1,3
            E = e_scale(k)*E_th(i)
            write(label,'(a,a,a,f0.2,a)') 'sigma_', trim(ion_name(i)),    &
                 '_at_', e_scale(k), '_Eth'
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 sigma_published(i,E), 1.0d-10)
         enddo
         do k = 1,3
            E = e_fixed(k)
            if (E .lt. E_th(i)) cycle
            write(label,'(a,a,a,i0,a)') 'sigma_', trim(ion_name(i)),      &
                 '_at_', nint(E), 'eV'
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 sigma_published(i,E), 1.0d-10)
         enddo

         ! Below the threshold there is no bound-free absorption.
         write(label,'(a,a)') 'sigma_below_threshold_', trim(ion_name(i))
         call check_absolute(trim(label),                                 &
              metal_photoion_sigma(i, 0.999d0*E_th(i)), 0.0d0, 0.0d0)

         ! THE HANDOVER ITSELF.  Just under E_max the returned value is the
         ! 1996 fit, just over it the 1995 fit of the same shell.  Where
         ! the two rows differ the assertion has content on both sides; Mg
         ! II is the one ion whose 1996 row IS its 1995 row, so its two
         ! sides agree by construction and the switch is invisible in it.
         if (E_max(i) .gt. 0.0d0) then
            E   = E_max(i) - 1.0d-6
            s96 = sigma_vfky96_published(i,E)
            write(label,'(a,a)') 'sigma_below_Emax_is_1996_',             &
                 trim(ion_name(i))
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 s96, 1.0d-10)
            ! Over E_max the outer shell is the 1995 fit AND every
            ! subshell whose turn-on is E_max itself has opened, so the
            ! reference is the sum.  The two are separated below.
            E   = E_max(i) + 1.0d-6
            s95 = sigma_vy95_published(i,E) + sigma_subshells_published(i,E)
            write(label,'(a,a)') 'sigma_above_Emax_is_1995_',             &
                 trim(ion_name(i))
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 s95, 1.0d-10)
         else
            ! E_max = 0: the 1995 fit at every energy, from the threshold up.
            E = 1.05d0*E_th(i)
            write(label,'(a,a)') 'sigma_is_1995_everywhere_',             &
                 trim(ion_name(i))
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 sigma_vy95_published(i,E), 1.0d-10)
         endif
      enddo

      ! EVERY SUBSHELL TURN-ON, 1e-6 eV on either side.  Under it the shell
      ! contributes nothing and the total is what it was; over it the total
      ! carries that shell's 1995 fit as well.  Both sides are compared
      ! with the published sum built here, so a shell that opened at the
      ! wrong energy, or not at all, fails on one side or the other.
      do i = 1,n_ion
         do m = 1,n_sub(i)
            E = E_on_sub(m,i) - 1.0d-6
            write(label,'(a,a,a,i0)') 'sigma_under_subshell_turn_on_',    &
                 trim(ion_name(i)), '_', m
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 sigma_published(i,E), 1.0d-10)
            E = E_on_sub(m,i) + 1.0d-6
            write(label,'(a,a,a,i0)') 'sigma_over_subshell_turn_on_',     &
                 trim(ion_name(i)), '_', m
            call check_relative(trim(label), metal_photoion_sigma(i,E),   &
                 sigma_published(i,E), 1.0d-10)
            ! The step itself: a shell that opens raises the total.
            write(label,'(a,a,a,i0)') 'sigma_rises_at_subshell_turn_on_', &
                 trim(ion_name(i)), '_', m
            call check_positive(trim(label),                              &
                 metal_photoion_sigma(i,E_on_sub(m,i)+1.0d-6)             &
                 - metal_photoion_sigma(i,E_on_sub(m,i)-1.0d-6))
         enddo
      enddo

      ! A column index outside the table is not a cross section.
      call check_absolute('sigma_column_zero',                            &
           metal_photoion_sigma(0, 5.0d1), 0.0d0, 0.0d0)
      call check_absolute('sigma_column_past_end',                        &
           metal_photoion_sigma(n_ion+1, 5.0d1), 0.0d0, 0.0d0)

      ! ONE definition of each metal threshold.  The photon grid puts
      ! mion_ethr on a bin EDGE and charges the photoelectron
      ! h nu - mion_ethr against it, so a cross section that turns on at
      ! its own copy of the threshold mis-charges the band between the two
      ! copies: above the edge, that band is integrated as zero; below it,
      ! a photon that cannot ionize is absorbed.  Each fit is evaluated
      ! 1e-6 eV on either side of the threshold species_table carries,
      ! which is narrower than any gap this can catch (the one it did
      ! catch was Mg II, 0.005 eV) and wider than the double-precision
      ! spacing of these energies.
      do i = 1,n_mion
         if (mion_iphot(i) .le. 0) cycle
         write(label,'(a,a)') 'sigma_zero_below_ethr_', trim(mion_name(i))
         call check_absolute(trim(label),                                 &
              metal_photoion_sigma(mion_iphot(i),                         &
                                   mion_ethr(i)-1.0d-6), 0.0d0, 0.0d0)
         write(label,'(a,a)') 'sigma_positive_above_ethr_',                &
              trim(mion_name(i))
         call check_positive(trim(label),                                 &
              metal_photoion_sigma(mion_iphot(i),                         &
                                   mion_ethr(i)+1.0d-6))
      enddo

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      double precision function sigma_published(i,E)
      ! The total: the outer shell with the fit that is valid at this
      ! energy (the 1996 one below E_max, the 1995 one of the same shell at
      ! and above it), plus every subshell already open.
      integer, intent(in) :: i
      real*8,  intent(in) :: E
      if (E .lt. E_th(i)) then
         sigma_published = 0.0d0
         return
      else if (E .ge. E_max(i)) then
         sigma_published = sigma_vy95_published(i,E)
      else
         sigma_published = sigma_vfky96_published(i,E)
      endif
      sigma_published = sigma_published + sigma_subshells_published(i,E)
      end function sigma_published

      !--------------!

      double precision function sigma_subshells_published(i,E)
      ! Verner & Yakovlev (1995) Eq. (1) for each subshell that is open.
      integer, intent(in) :: i
      real*8,  intent(in) :: E
      integer :: m
      sigma_subshells_published = 0.0d0
      do m = 1,n_sub(i)
         if (E .lt. E_on_sub(m,i)) cycle
         sigma_subshells_published = sigma_subshells_published            &
              + sigma_one_shell(E, E_0_sub(m,i), sigma_0_sub(m,i),        &
                                y_a_sub(m,i), P_sub(m,i), y_w_sub(m,i),   &
                                l_sub(m,i))
      enddo
      end function sigma_subshells_published

      !--------------!

      double precision function sigma_one_shell(E,E_0s,s_0,y_as,Ps,y_ws,ls)
      ! Verner & Yakovlev (1995) Eq. (1), one shell, written here from the
      ! paper: y = E/E_0, Q = 5.5 + l - P/2.
      real*8,  intent(in) :: E, E_0s, s_0, y_as, Ps, y_ws
      integer, intent(in) :: ls
      real*8 :: y, Q
      y = E/E_0s
      Q = 5.5d0 + dble(ls) - 0.5d0*Ps
      sigma_one_shell = s_0*((y-1.0d0)**2 + y_ws**2)                      &
           *y**(-Q)*(1.0d0 + sqrt(y/y_as))**(-Ps)
      end function sigma_one_shell

      !--------------!

      double precision function sigma_vfky96_published(i,E)
      ! Verner et al. (1996) Eqs. (1)-(4), written here from the paper and
      ! not from the production module: the offset x, the asymptote y, and
      ! Q = 5.5 - 0.5 P (their Table 1 carries no l term).
      integer, intent(in) :: i
      real*8,  intent(in) :: E
      real*8 :: x, y, Q
      x = E/E_0(i) - y_0(i)
      y = sqrt(x*x + y_1(i)*y_1(i))
      Q = 5.5d0 - 0.5d0*P_fit(i)
      sigma_vfky96_published = sigma_0(i)*((x-1.0d0)**2 + y_w(i)**2)      &
           *y**(-Q)*(1.0d0 + sqrt(y/y_a(i)))**(-P_fit(i))
      end function sigma_vfky96_published

      !--------------!

      double precision function sigma_vy95_published(i,E)
      ! Verner & Yakovlev (1995) Eq. (1) in y = E/E_0, with the orbital
      ! quantum number of the shell in Q = 5.5 + l - 0.5 P.
      integer, intent(in) :: i
      real*8,  intent(in) :: E
      real*8 :: y, Q
      y = E/E_0_vy95(i)
      Q = 5.5d0 + dble(l_shell(i)) - 0.5d0*P_vy95(i)
      sigma_vy95_published = sigma_0_vy95(i)*((y-1.0d0)**2 + y_w_vy95(i)**2) &
           *y**(-Q)*(1.0d0 + sqrt(y/y_a_vy95(i)))**(-P_vy95(i))
      end function sigma_vy95_published

      end program metal_photoionization_fits
