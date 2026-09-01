      program a2_m2_g3_run_check
      ! Gate G3 of docs/a2_oxygen_option_design.md, measured on a finished
      ! run rather than algebraically: with every photolysis rate at zero,
      ! the solved O / OH / H2O partition must BE the chemical equilibrium
      ! of the cell's own (T, n_H2, n_H0).  It is the test that makes the
      ! thermodynamically reversed pairs trustworthy -- the reverse rates
      ! are never transcribed, so nothing but this says they are right.
      !
      ! The algebraic half of the same gate is section 2 of
      ! a2_m2_kinetics_check, which puts the two rows on the equilibrium
      ! partition and finds the residuals at round-off.  This program adds
      ! the half that the algebra cannot reach: that the coupled solve, the
      ! row scaling, the closure and the write-back deliver that partition.
      !
      ! Not part of any build.  From src/tests/a2_m2/:
      !
      !   gfortran -O2 -o a2_m2_g3_run_check.x \
      !       ../../modules/lower_atmosphere/oxygen_rates.f90 \
      !       a2_m2_g3_run_check.f90
      !   ./a2_m2_g3_run_check.x <run_dir>
      !
      ! where <run_dir> is a run made with "Oxygen chemistry: True" and
      ! every FUV band flux at zero.
      !
      ! TWO CLASSES OF CELL ARE EXCLUDED, and neither exclusion hides an
      ! error.  A cell with no H2 or no atomic H has no reaction to run at
      ! all.  A cell in which one DIRECTION of the dominant pair is slower
      ! than 1e-10 s^-1 -- 97 K with 0.02 cm^-3 of atomic H, in the
      ! snapshot this was first run on -- satisfies both rows to the
      ! solver's tolerance over a range of partitions: a local steady state
      ! is only determined where something runs in both directions, and
      ! where it is not, the partition is a question for the transport, not
      ! for this gate.  Both counts are printed.

      use oxygen_rates
      implicit none
      integer, parameter :: mxc = 20000
      integer :: j, k, u, ios, nskip, nfroz, nused, nbad
      real*8  :: v(15), r, T, nO, nOII, nOIII, nOH, nH2O, nCO
      real*8  :: nH2(mxc), nH0(mxc), w(41)
      real*8  :: f_oh, f_h2o, fam, dev, worst, wr, kf, kr
      character(len=200000) :: line
      character(len=512)    :: run_dir

      call get_command_argument(1, run_dir)
      if (len_trim(run_dir) == 0) run_dir = '.'

      ! --- n_H2 and n_H0 of every cell, from Ion_species.txt -------------
      ! Layout: r, 6 H/He columns, 27 metal stages, H2 H2p H3p HeHp,
      ! OH H2O CO = 41 fields.  Column 2 is H I and column 35 is H2.
      open(newunit=u, file=trim(run_dir)//'/output/Ion_species.txt',       &
           status='old')
      j = 0
      do
        read(u,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#') cycle
        j = j + 1
        if (j > mxc) stop 'a2_m2_g3_run_check: raise mxc'
        read(line,*) w(1:41)
        nH0(j) = w(2)
        nH2(j) = w(35)
      end do
      close(u)

      ! --- the solved partition, from Oxygen_chemistry.txt ---------------
      open(newunit=u, file=trim(run_dir)//'/output/Oxygen_chemistry.txt',  &
           status='old')
      k = 0; nskip = 0; nfroz = 0; nused = 0; nbad = 0
      worst = 0.0d0; wr = 0.0d0
      do
        read(u,'(A)',iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#') cycle
        k = k + 1
        read(line,*) v(1:15)
        r = v(1); T = v(2); nO = v(3); nOII = v(4); nOIII = v(5)
        nOH = v(6); nH2O = v(7); nCO = v(8)
        ! The equilibrium partition is among the NEUTRAL carriers O, OH and
        ! H2O.  The ionized stages are a separate balance that no reaction
        ! of this network touches, so they are not in the denominator.
        fam = nO + nOH + nH2O
        if (fam <= 0.0d0) cycle
        if (nH2(k) <= 0.0d0 .or. nH0(k) <= 0.0d0) then
          nskip = nskip + 1
          cycle
        end if
        kf = rk_O1_OH_H2_water(T)
        kr = rate_from_detailed_balance(kf, (/ ith_OH, ith_H2 /),         &
                                            (/ ith_H2O, ith_H /), T)
        if (kf*nH2(k) < 1.0d-10 .or. kr*nH0(k) < 1.0d-10) then
          nfroz = nfroz + 1
          cycle
        end if
        nused = nused + 1
        call oxygen_chemical_equilibrium_fractions(T, nH2(k), nH0(k),     &
                                                   f_oh, f_h2o)
        dev = max(abs(nH2O/fam - f_h2o), abs(nOH/fam - f_oh))
        if (dev > 1.0d-6) nbad = nbad + 1
        if (dev > worst) then
          worst = dev
          wr    = r
        end if
      end do
      close(u)

      write(*,'(a,a)')      ' run:   ', trim(run_dir)
      write(*,'(a,i6)')     ' cells: ', k
      write(*,'(a,i6,a,i6,a,i6)') ' used: ', nused,                       &
           '   no H2/H0: ', nskip, '   one-way (frozen): ', nfroz
      write(*,'(a,i6)')     ' cells above 1e-6: ', nbad
      write(*,'(a,es12.4,a,f10.5)')                                       &
           ' G3 max |solved - chemical equilibrium| = ', worst,           &
           '  at r/Rp = ', wr
      end program a2_m2_g3_run_check
