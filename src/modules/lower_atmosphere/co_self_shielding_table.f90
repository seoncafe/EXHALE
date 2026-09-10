      module co_self_shielding_table
      ! Shielding function Theta[N(12CO), N(H2)] of CO photodissociation,
      ! from Visser, van Dishoeck and Black (2009), A&A 503, 323.
      !
      ! WHAT THETA IS.  Their eq. (2) writes the rate of one isotopologue as
      !
      !     k_i = chi k0_i Theta_i exp(-gamma A_V) ,
      !
      ! an unattenuated rate in a stated radiation field, a LINE shielding
      ! function of the star-ward columns, and a continuum term.  Theta
      ! "accounts for self-shielding and shielding by H, H2 and the other CO
      ! isotopologues" (their sec. 3.3).  CO predissociates in 37 discrete
      ! bands between 912.70 and 1076.08 A (their Table 1), so its shielding
      ! is a line problem and not an exp(-tau): this module is the CO
      ! counterpart of h2_self_shielding_table, not of water_photolysis.
      !
      ! Reducing the line-by-line problem to a function of two columns is
      ! the paper's own construction, and it is justified in their sec. 5.1:
      ! "The transition from atomic to molecular hydrogen occurs much closer
      ! to the edge of the cloud than the C+-C-CO transition, so the column
      ! density of atomic H is roughly constant at the depths where
      ! shielding of CO is important.  In addition, H shields CO by only a
      ! few per cent.  Therefore, it is a good approximation to compute the
      ! shielding functions on a grid of CO and H2 column densities, while
      ! taking a constant column of H."  Their sec. 5.2 measures what that
      ! costs against the full integration: "The rate from our approximate
      ! method is within 10% of the 'real' rate in 98.3% of all points
      ! (Fig. 6).  In no cases is the difference between the approximate
      ! rates and the full model more than 40%."
      !
      ! ---------------------------------------------------------------
      ! 1. WHICH TABLE IS SHIPPED, AND THE DOMAIN THAT COMES WITH IT
      !
      ! The paper gives four parameter sets.  The one carried here is
      ! TABLE 6: b(CO) = 0.3 km/s, T_ex(CO) = 50 K, T_ex(H2) = 501.5 K,
      ! N(12CO)/N(13CO) = 69, "Additional rotational lines of CO and H2 were
      ! included as described in Sect. 4.4" (its footnote a).  It is the
      ! warmest of the four sets and therefore the closest to the 800-3000 K
      ! molecular layer this code integrates.  Table 5, the paper's
      ! reference set at T_ex(CO) = 5 K, is carried beside it as
      ! co_self_shielding_texc5 and is NOT called by the rate: it exists so
      ! that the spread between the two excitation temperatures can be
      ! measured from this tree, which is what
      ! src/tests/physics_probe/co_shielding_table.f90 does.
      !
      ! NEITHER TABLE IS AT THE EXCITATION TEMPERATURE OF THIS GAS, and the
      ! reason is a data limit rather than a choice of convenience.  READ,
      ! their sec. 4.4: "T ex (CO) is raised from 4 to 512 K in steps of
      ! factors of two.  The v'' = 1 vibrational level of 12CO lies at
      ! 2143 cm^-1 above the v'' = 0 level, so it starts to be thermally
      ! populated at ~500 K.  No data are available on dissociative
      ! transitions out of this level, so we choose not to go to higher
      ! excitation temperatures."  Their sec. 5.2 states the size of the
      ! error that follows: the Table 5 functions "can easily give
      ! photodissociation rates off by a factor of two when applied to a
      ! high-density, high-temperature PDR".
      !
      ! So co_shield_tex_limit_K() below is a STATED OUT-OF-DOMAIN
      ! TEMPERATURE, not a clamp.  Above it the table is still evaluated --
      ! there is nothing better to evaluate -- and the cells above it are
      ! counted and reported, on the same reading that makes the H3+ cooling
      ! domain record informational.  A factor of two in Theta is a factor
      ! of two in the rate, so it has to be printed.
      !
      ! The gas also sits between the tables in Doppler width and not above
      ! them: b(CO) = sqrt(2kT/m) is 0.42 km/s at 300 K and 1.34 km/s at
      ! 3000 K, against 0.3 km/s for Tables 5 and 6 and 3.0 km/s for
      ! Table 7.  No shipped set is closer on both axes at once, which is
      ! why the excitation temperature is the axis the choice was made on.
      !
      ! ---------------------------------------------------------------
      ! 2. THE GRID, AND WHAT HAPPENS OUTSIDE IT
      !
      !     log10 N(12CO) [cm^-2] :  0, 13, 14, 15, 16, 17, 18, 19
      !     log10 N(H2)   [cm^-2] :  0, 19, 20, 21, 22, 23
      !
      ! Only the 12CO block of the table is transcribed.  EXHALE carries one
      ! carbon and one oxygen species with no isotopic structure, so the
      ! five isotopologue blocks of every table -- which are the paper's own
      ! subject -- can never be called from here.
      !
      ! Interpolation is BILINEAR IN (log10 N_CO, log10 N_H2, log10 Theta).
      ! Theta falls by more than three decades along the CO axis and by more
      ! than six along the H2 axis, so an interpolation linear in Theta
      ! between two rows that differ by five decades would return the upper
      ! row almost everywhere between them.  Interpolating the logarithm is
      ! the same convention h2_self_shielding_table uses on its own column
      ! axis and for the same reason.
      !
      ! Outside the grid the EDGE VALUE is returned and no extrapolation is
      ! made.  Their sec. 5.1 justifies that at the deep end:
      ! "photodissociation at these depths is typically already so slow a
      ! process that it is no longer the dominant destruction pathway for
      ! CO".  Below the bottom of an axis (a column under 1 cm^-2, including
      ! a column of exactly zero) the same edge value is returned, and there
      ! it is the unshielded limit Theta = 1 that the table itself carries.
      !
      ! ---------------------------------------------------------------
      ! References: Visser, van Dishoeck & Black (2009), A&A 503, 323
      ! (publisher PDF, references/Visser_2009A&A_503_323.pdf): sec. 3.3
      ! eq. (2), secs. 4.4, 5.1 and 5.2, Table 6 (Online Material p. 1) for
      ! the shipped set and Table 5 (p. 334) for the spread.

      implicit none
      private
      public :: co_self_shielding, co_self_shielding_texc5,               &
                co_shield_max_co_column, co_shield_max_h2_column,         &
                co_shield_tex_limit_K,                                    &
                n_co_shield_co, n_co_shield_h2,                           &
                co_shield_log_co, co_shield_log_h2

      integer, parameter :: dp = kind(1.0d0)

      integer, parameter :: n_co_shield_co = 8
      integer, parameter :: n_co_shield_h2 = 6

      real(dp), parameter :: co_shield_log_co(n_co_shield_co) =           &
           (/ 0.0d0, 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0,      &
              19.0d0 /)
      real(dp), parameter :: co_shield_log_h2(n_co_shield_h2) =           &
           (/ 0.0d0, 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0 /)

      ! Visser et al. Table 6, the 12CO block, transcribed row by row in the
      ! order the table prints (one row per log N(H2), the eight log N(CO)
      ! entries across).  The block's own header carries the unattenuated
      ! rate k0 = 2.590e-10 s^-1 in the Draine (1978) field; this code does
      ! not use k0 (co_photodissociation.f90 states why), so it is recorded
      ! here and nowhere else.
      real(dp), parameter :: k0_draine_table6 = 2.590d-10
      real(dp), parameter ::                                              &
        theta6(n_co_shield_co, n_co_shield_h2) = reshape( (/              &
          1.000d0,   9.405d-1, 7.046d-1, 4.015d-1,                        &
          9.964d-2,  1.567d-2, 3.162d-3, 4.839d-4,                        &
          7.546d-1,  6.979d-1, 4.817d-1, 2.577d-1,                        &
          6.505d-2,  1.135d-2, 2.369d-3, 3.924d-4,                        &
          5.752d-1,  5.228d-1, 3.279d-1, 1.559d-1,                        &
          3.559d-2,  6.443d-3, 1.526d-3, 2.751d-4,                        &
          2.493d-1,  2.196d-1, 1.135d-1, 4.062d-2,                        &
          7.864d-3,  1.516d-3, 4.448d-4, 9.367d-5,                        &
          1.550d-3,  1.370d-3, 6.801d-4, 2.127d-4,                        &
          5.051d-5,  1.198d-5, 6.553d-6, 3.937d-6,                        &
          8.492d-8,  8.492d-8, 8.492d-8, 8.492d-8,                        &
          8.492d-8,  8.492d-8, 8.488d-8, 8.453d-8 /),                     &
          (/ n_co_shield_co, n_co_shield_h2 /) )

      ! Visser et al. Table 5, the 12CO block: b(CO) = 0.3 km/s,
      ! T_ex(CO) = 5 K, T_ex(H2) = 51.5 K, k0 = 2.592e-10 s^-1.  The paper's
      ! reference set.  NOT CALLED BY THE RATE -- it is here so that the
      ! spread between T_ex(CO) = 5 K and 50 K can be measured from this
      ! tree and reported as the model's own uncertainty on that axis.
      real(dp), parameter :: k0_draine_table5 = 2.592d-10
      real(dp), parameter ::                                              &
        theta5(n_co_shield_co, n_co_shield_h2) = reshape( (/              &
          1.000d0,   8.080d-1, 5.250d-1, 2.434d-1,                        &
          5.467d-2,  1.362d-2, 3.378d-3, 5.240d-4,                        &
          8.176d-1,  6.347d-1, 3.891d-1, 1.787d-1,                        &
          4.297d-2,  1.152d-2, 2.922d-3, 4.662d-4,                        &
          7.223d-1,  5.624d-1, 3.434d-1, 1.540d-1,                        &
          3.515d-2,  9.231d-3, 2.388d-3, 3.899d-4,                        &
          3.260d-1,  2.810d-1, 1.953d-1, 8.726d-2,                        &
          1.907d-2,  4.768d-3, 1.150d-3, 1.941d-4,                        &
          1.108d-2,  1.081d-2, 9.033d-3, 4.441d-3,                        &
          1.102d-3,  2.644d-4, 7.329d-5, 1.437d-5,                        &
          3.938d-7,  3.938d-7, 3.936d-7, 3.923d-7,                        &
          3.901d-7,  3.893d-7, 3.890d-7, 3.875d-7 /),                     &
          (/ n_co_shield_co, n_co_shield_h2 /) )

      contains

      ! Shielding function Theta of CO photodissociation for star-ward
      ! columns N_CO and N_H2 [cm^-2], from the shipped table (Table 6).
      ! Dimensionless, 1 in the unshielded limit.
      double precision function co_self_shielding(N_CO, N_H2) result(theta)
      real(dp), intent(in) :: N_CO, N_H2
      theta = co_shield_of_block(theta6, N_CO, N_H2)
      end function co_self_shielding

      ! The same at the paper's reference excitation temperature,
      ! T_ex(CO) = 5 K (Table 5).  Its only caller is the test that reports
      ! the spread between the two sets; the rate never reads it.
      double precision function co_self_shielding_texc5(N_CO, N_H2)       &
                                result(theta)
      real(dp), intent(in) :: N_CO, N_H2
      theta = co_shield_of_block(theta5, N_CO, N_H2)
      end function co_self_shielding_texc5

      ! Bilinear interpolation of one isotopologue block in
      ! (log10 N_CO, log10 N_H2, log10 Theta), with the edge value outside
      ! the grid on either axis.  One expression, so the shipped table and
      ! the comparison table cannot be interpolated by two different rules.
      double precision function co_shield_of_block(blk, N_CO, N_H2)       &
                                result(theta)
      real(dp), intent(in) :: blk(n_co_shield_co, n_co_shield_h2)
      real(dp), intent(in) :: N_CO, N_H2
      real(dp) :: xc, xh, tx, ty, t00, t10, t01, t11, a0, a1
      integer  :: ic, ih

      call co_shield_locate(co_shield_log_co, n_co_shield_co, N_CO, ic, tx)
      call co_shield_locate(co_shield_log_h2, n_co_shield_h2, N_H2, ih, ty)
      xc = tx
      xh = ty

      t00 = log10(blk(ic,   ih  ))
      t10 = log10(blk(ic+1, ih  ))
      t01 = log10(blk(ic,   ih+1))
      t11 = log10(blk(ic+1, ih+1))

      a0 = t00 + xc*(t10 - t00)
      a1 = t01 + xc*(t11 - t01)
      theta = 10.0d0**(a0 + xh*(a1 - a0))
      end function co_shield_of_block

      ! Interval index i and fractional position t of log10(max(N,1)) on a
      ! monotone axis, clamped to the grid so that the edge value is
      ! returned outside it (header sec. 2).  A column of zero -- the
      ! outermost cell has nothing above it -- gives log10(1) = 0, the
      ! bottom knot of both axes, which is the unshielded entry.
      subroutine co_shield_locate(axis, n, col, i, t)
      real(dp), intent(in)  :: axis(*)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: col
      integer,  intent(out) :: i
      real(dp), intent(out) :: t
      real(dp) :: x
      x = log10(max(col, 1.0d0))
      if (x .le. axis(1)) then
         i = 1
         t = 0.0d0
         return
      endif
      if (x .ge. axis(n)) then
         i = n - 1
         t = 1.0d0
         return
      endif
      i = 1
      do while (i .lt. n-1 .and. x .ge. axis(i+1))
         i = i + 1
      enddo
      t = (x - axis(i))/(axis(i+1) - axis(i))
      end subroutine co_shield_locate

      ! Top of the CO column axis [cm^-2].  Above it the table returns its
      ! edge value rather than a calculated one.
      double precision function co_shield_max_co_column() result(N)
      N = 10.0d0**co_shield_log_co(n_co_shield_co)
      end function co_shield_max_co_column

      ! Top of the H2 column axis [cm^-2], same convention.
      double precision function co_shield_max_h2_column() result(N)
      N = 10.0d0**co_shield_log_h2(n_co_shield_h2)
      end function co_shield_max_h2_column

      ! Highest CO excitation temperature the shielding model was computed
      ! for [K].  Visser et al. sec. 4.4: above about 500 K the v'' = 1
      ! level is populated and no data exist for dissociating transitions
      ! out of it, so their grid stops at 512 K.  A cell hotter than this is
      ! recorded as out of domain; nothing is clamped.
      double precision function co_shield_tex_limit_K() result(T)
      T = 512.0d0
      end function co_shield_tex_limit_K

      ! End of module
      end module co_self_shielding_table
