      module lower_atmosphere_profile
      ! The lower atmosphere's solution over an interval of pressure, handed
      ! to EXHALE as a table instead of the single-level scalars of base.inp.
      !
      ! Physics this module carries
      ! ---------------------------
      ! A lower-atmosphere model (photochemistry + climate) solves the same
      ! gas EXHALE solves, over a pressure interval that reaches ABOVE the
      ! level where EXHALE places its base.  Where the two overlap, both are
      ! defined, and the composition, the temperature and the eddy mixing are
      ! properties of that shared gas rather than of either code.  This module
      ! reads that solution and hands EXHALE
      !
      !   - the state at the matching level p_match (T, r, q_H2, the elemental
      !     nuclei ratios He/H and El/H), which is what the base EOS and the
      !     elemental reservoirs are built from;
      !   - the eddy diffusion coefficient K_zz as a PROFILE.  K_zz is the one
      !     thing in the handoff that is not a single number: the homopause,
      !     the only consequential thing K_zz does, is where K_zz = D_12, a
      !     property of two profiles, and the molecular coefficient rises by
      !     nearly five decades between the base and 1.2 R_p
      !     (docs/eddy_diffusion_kzz.tex, LHS1140b/kzz_decision.md).
      !
      ! Schema: docs/input_schema.md section 2d.  Example file:
      ! examples/17_lower_profile/lower_atmosphere_profile.dat.
      !
      ! Interpolation.  Every intensive quantity is interpolated linearly in
      ! log p, which is the variable both models solve on; the two codes'
      ! radius scales are equal only if their hydrostatic integrations agree.
      ! Extrapolation is never performed.  Above the shallowest level the file
      ! carries, K_zz takes that level's value (the eddy coefficient is a
      ! lower-atmosphere property and the file makes no statement above its
      ! top); nothing else is taken from the profile there.
      !
      ! When no profile is given every entry point here is inert and the
      ! scalar base.inp path of read_base_inp is unchanged, bit for bit.

      use global_parameters

      implicit none

      private

      ! ---- what a run states about its profile -------------------------- !
      public :: lap_in_use, lap_file
      ! ---- header fields ------------------------------------------------ !
      public :: lap_solution_id, lap_source_code, lap_source_version
      public :: lap_mechanism, lap_stellar_flux, lap_notes
      public :: lap_p_match_bar, lap_p_top_bar, lap_p_deep_bar
      public :: lap_trial_flux_H, lap_trial_flux_He
      public :: lap_iteration, lap_steady, lap_iterable
      public :: lap_r_top_RJ
      ! ---- entry points -------------------------------------------------- !
      public :: read_lower_atmosphere_profile
      public :: lap_value_at_match
      public :: lap_element_ratio_at_match
      public :: eddy_diffusion_on_grid
      public :: lap_report_provenance
      ! ---- the elemental flux measured over two radial windows ----------- !
      ! Filled by binary_element_diffusion when a profile is in use; written
      ! by write_resolved_config.  lap_flux_measured stays .false. until a
      ! diffusion step has produced numbers, so an unmeasured run says so.
      !
      ! TWO windows are reported, because they answer two different questions
      ! and only measurement says which one carries a usable number:
      !
      !  (a) the OVERLAP window, the interval both models describe: from the
      !      top of the base sound-wave region up to the radius the profile
      !      reaches.  Its lower edge is measured, not assumed -- the first
      !      face at which the 5-face moving spread of r^2 rho v falls below
      !      10 per cent -- and falls back to 1.02 R_p when no face qualifies.
      !  (b) the STEADY-FLUX window, r >= r_esc, the escape region the solver
      !      itself uses to declare the wind steady.  At a steady state the
      !      elemental flux is radius independent, so a flux measured here IS
      !      the flux through the matching level; that identity is what makes
      !      this window a legitimate statement about the handoff and not a
      !      different quantity.
      public :: lap_flux_measured, lap_flux_window_empty
      public :: lap_flux_r_lo_Rp, lap_flux_r_hi_Rp, lap_flux_r_lo_measured
      public :: lap_FH_median, lap_FH_spread
      public :: lap_FHe_median, lap_FHe_spread
      public :: lap_Mdot_median, lap_Mdot_spread
      public :: lap_flux_nface
      public :: lap_steady_r_lo_Rp, lap_steady_nface
      public :: lap_steady_FH_median, lap_steady_FH_spread
      public :: lap_steady_FHe_median, lap_steady_FHe_spread
      public :: lap_steady_Mdot_median, lap_steady_Mdot_spread

      logical            :: lap_in_use = .false.
      character(len=256) :: lap_file   = ''

      character(len=128) :: lap_solution_id    = ''
      character(len=64)  :: lap_source_code    = ''
      character(len=128) :: lap_source_version = ''
      character(len=192) :: lap_mechanism      = ''
      character(len=192) :: lap_stellar_flux   = ''
      ! The climate solve of milestone E3 states its deep boundary, its
      ! tropopause and the two elemental sums that bracket the cold trap in
      ! this one field, so it is long enough to hold that record whole.
      character(len=1000) :: lap_notes         = ''
      real*8  :: lap_p_match_bar  = -1.0d0
      real*8  :: lap_p_top_bar    = -1.0d0
      real*8  :: lap_p_deep_bar   = -1.0d0
      real*8  :: lap_trial_flux_H  = 0.0d0
      real*8  :: lap_trial_flux_He = 0.0d0
      integer :: lap_iteration    = -1
      logical :: lap_steady       = .false.
      ! .true. only if both `iteration` and `trial_flux_H` were stated: a
      ! hand-written profile is a legitimate one-shot input, it just cannot
      ! be driven by the flux-closure iteration.
      logical :: lap_iterable     = .false.
      real*8  :: lap_r_top_RJ     = -1.0d0

      logical :: lap_flux_measured     = .false.
      logical :: lap_flux_window_empty = .false.
      integer :: lap_flux_nface        = 0
      real*8  :: lap_flux_r_lo_Rp      = -1.0d0
      real*8  :: lap_flux_r_hi_Rp      = -1.0d0
      ! .true. when the lower edge came from the 10 per cent spread rule,
      ! .false. when no face qualified and the 1.02 R_p fallback was taken.
      logical :: lap_flux_r_lo_measured = .false.
      real*8  :: lap_FH_median  = 0.0d0, lap_FH_spread  = 0.0d0
      real*8  :: lap_FHe_median = 0.0d0, lap_FHe_spread = 0.0d0
      real*8  :: lap_Mdot_median = 0.0d0, lap_Mdot_spread = 0.0d0

      integer :: lap_steady_nface   = 0
      real*8  :: lap_steady_r_lo_Rp = -1.0d0
      real*8  :: lap_steady_FH_median   = 0.0d0, lap_steady_FH_spread   = 0.0d0
      real*8  :: lap_steady_FHe_median  = 0.0d0, lap_steady_FHe_spread  = 0.0d0
      real*8  :: lap_steady_Mdot_median = 0.0d0, lap_steady_Mdot_spread = 0.0d0

      ! ---- the table ----------------------------------------------------- !
      integer :: nlev = 0, ncol = 0
      character(len=24), allocatable :: colname(:)
      real*8,            allocatable :: coltab(:,:)     ! (nlev, ncol)
      real*8,            allocatable :: logp(:)         ! ln of column `p`
      integer :: ic_p = 0, ic_r = 0, ic_kzz = 0

      contains

      ! ------------------------------------------------------------------ !

      subroutine read_lower_atmosphere_profile
      ! Read the file named by the "Lower atmosphere profile:" key.  Absent
      ! key: return with lap_in_use = .false. and nothing touched.

      ! Wide enough for a table row of ~150 columns at full double precision.
      character(len=4096) :: raw
      character(len=24)  :: key_s
      character(len=:), allocatable :: val_s
      integer :: iu, ios, ilev, icol, ndat
      logical :: ex

      if (len_trim(lap_file) .eq. 0) return

      inquire(file=trim(lap_file), exist=ex)
      if (.not. ex) then
         write(*,*) '(lower_atmosphere_profile) ERROR: "Lower atmosphere'//&
                    ' profile: '//trim(lap_file)//'"'
         write(*,*) '  was set but the file does not exist in the run'//   &
                    ' directory. Remove the key to run on the scalar'
         write(*,*) '  base.inp handoff, or put the file in place.'
         error stop 1
      endif

      write(*,*) '(lower_atmosphere_profile) Reading '//trim(lap_file)//   &
                 ' (lower-atmosphere solution)..'

      ! ---- pass 1: header lines and the number of data rows -------------- !
      ndat = 0
      if (allocated(coltab)) deallocate(coltab)
      if (allocated(logp))   deallocate(logp)
      open(newunit=iu, file=trim(lap_file), status='old', action='read')
      do
         read(iu,'(A)',iostat=ios) raw
         if (ios .ne. 0) exit
         if (len_trim(raw) .eq. 0) cycle
         if (index(adjustl(raw),'#') .eq. 1) then
            call split_header(raw, key_s, val_s)
            ! A header key with no value states nothing; skip it rather than
            ! read past the end of an empty string.
            if (len_trim(val_s) .eq. 0) cycle
            select case (trim(key_s))
               case ('columns');  call set_column_names(val_s)
               case ('solution_id');       lap_solution_id    = val_s
               case ('source_code');       lap_source_code    = val_s
               case ('source_version');    lap_source_version = val_s
               case ('mechanism');         lap_mechanism      = val_s
               case ('stellar_flux');      lap_stellar_flux   = val_s
               case ('notes');             lap_notes          = val_s
               case ('p_match_bar');   read(val_s,*) lap_p_match_bar
               case ('p_top_bar');     read(val_s,*) lap_p_top_bar
               case ('p_deep_bar');    read(val_s,*) lap_p_deep_bar
               case ('trial_flux_H');  read(val_s,*) lap_trial_flux_H
               case ('trial_flux_He'); read(val_s,*) lap_trial_flux_He
               case ('iteration');     read(val_s,*) lap_iteration
               case ('reached_steady_state')
                  lap_steady = (val_s(1:1) .eq. 'T' .or. val_s(1:1) .eq. 't')
               case default
                  ! Header keys this version has no consumer for are carried
                  ! by the file and ignored here, exactly as unknown columns
                  ! are kept and not dropped.
            end select
            cycle
         endif
         ndat = ndat + 1
      enddo

      if (ncol .le. 0) call refuse('the "# columns:" schema line is missing')
      if (ndat .lt. 2) call refuse('fewer than two levels in the table')

      nlev = ndat
      allocate(coltab(nlev, ncol), logp(nlev))

      ! ---- pass 2: the table --------------------------------------------- !
      rewind(iu)
      ilev = 0
      do
         read(iu,'(A)',iostat=ios) raw
         if (ios .ne. 0) exit
         if (len_trim(raw) .eq. 0) cycle
         if (index(adjustl(raw),'#') .eq. 1) cycle
         ilev = ilev + 1
         read(raw,*,iostat=ios) (coltab(ilev,icol), icol = 1, ncol)
         if (ios .ne. 0) then
            write(*,'(A,I0,A,I0,A)') ' (lower_atmosphere_profile) ERROR:'//&
               ' data row ', ilev, ' does not carry the ', ncol,           &
               ' columns the schema line names.'
            write(*,*) '  row: '//trim(raw)
            error stop 1
         endif
      enddo
      close(iu)

      call locate_required_columns
      call validate_table

      lap_in_use = .true.
      ! A stated iteration index is what makes a profile drivable by the
      ! flux-closure loop; the trial flux is legitimately zero on the first
      ! iteration, so a zero value says nothing.  A hand-written profile with no
      ! iteration index runs, it just cannot be iterated.
      lap_iterable = (lap_iteration .ge. 0)

      end subroutine read_lower_atmosphere_profile

      ! ------------------------------------------------------------------ !

      subroutine split_header(raw, key_s, val_s)
      ! "# key value ..." -> key (trailing ':' stripped) and the rest.
      character(len=*), intent(in)  :: raw
      character(len=24), intent(out) :: key_s
      character(len=:), allocatable, intent(out) :: val_s
      character(len=4096) :: body
      integer :: ib, lk
      body  = adjustl(raw)
      body  = adjustl(body(2:))               ! drop the '#'
      ib    = index(trim(body), ' ')
      if (ib .eq. 0) then
         key_s = trim(body)
         val_s = ''
      else
         key_s = body(1:ib-1)
         val_s = trim(adjustl(body(ib+1:)))
      endif
      lk = len_trim(key_s)
      if (lk .gt. 0) then
         if (key_s(lk:lk) .eq. ':') key_s = key_s(1:lk-1)
      endif
      end subroutine split_header

      ! ------------------------------------------------------------------ !

      subroutine set_column_names(val_s)
      ! Store the names of the "# columns:" line, in file order.  Unknown
      ! names are kept, not dropped: a profile may carry columns this version
      ! of EXHALE has no consumer for, and they must survive a round trip.
      character(len=*), intent(in) :: val_s
      character(len=4096) :: rest
      integer :: ib, icol, nn
      rest = adjustl(val_s)
      nn = 0
      do
         if (len_trim(rest) .eq. 0) exit
         nn = nn + 1
         ib = index(trim(rest), ' ')
         if (ib .eq. 0) exit
         rest = adjustl(rest(ib+1:))
      enddo
      if (nn .eq. 0) return
      if (allocated(colname)) deallocate(colname)
      allocate(colname(nn))
      ncol = nn
      rest = adjustl(val_s)
      do icol = 1, nn
         ib = index(trim(rest), ' ')
         if (ib .eq. 0) then
            colname(icol) = trim(rest)
            rest = ''
         else
            colname(icol) = rest(1:ib-1)
            rest = adjustl(rest(ib+1:))
         endif
      enddo
      end subroutine set_column_names

      ! ------------------------------------------------------------------ !

      integer function column_index(nm)
      ! Position of a named column, 0 if the file does not carry it.  Columns
      ! are found BY NAME everywhere: the producers order their element list
      ! differently from run to run.
      character(len=*), intent(in) :: nm
      integer :: icol
      column_index = 0
      if (.not. allocated(colname)) return
      do icol = 1, ncol
         if (trim(colname(icol)) .eq. trim(nm)) then
            column_index = icol
            return
         endif
      enddo
      end function column_index

      ! ------------------------------------------------------------------ !

      subroutine locate_required_columns
      character(len=8), parameter :: need(*) = [ character(len=8) ::       &
         'p', 'r', 'T', 'n_tot', 'rho', 'Kzz', 'q_H2', 'q_H', 'X_He' ]
      integer :: k
      do k = 1, size(need)
         if (column_index(trim(need(k))) .eq. 0)                           &
            call refuse('the required column "'//trim(need(k))//           &
                        '" is missing from the schema line')
      enddo
      ic_p   = column_index('p')
      ic_r   = column_index('r')
      ic_kzz = column_index('Kzz')
      end subroutine locate_required_columns

      ! ------------------------------------------------------------------ !

      subroutine validate_table
      ! Monotonicity, the overlap requirement, and the match level's coverage.
      ! Each of these makes an interpolation below either meaningless or an
      ! extrapolation, so they are refusals, not warnings.
      integer :: ilev
      real*8  :: p_shallow, p_deepest

      do ilev = 1, nlev
         if (coltab(ilev, ic_p) .le. 0.0d0)                                &
            call refuse('a non-positive pressure in column "p"')
         logp(ilev) = log(coltab(ilev, ic_p))
      enddo
      do ilev = 2, nlev
         if (coltab(ilev, ic_p) .ge. coltab(ilev-1, ic_p))                 &
            call refuse('column "p" is not strictly decreasing (the table'//&
                        ' runs deep to shallow)')
      enddo
      p_deepest = coltab(1,    ic_p)
      p_shallow = coltab(nlev, ic_p)

      if (lap_p_match_bar .le. 0.0d0)                                      &
         call refuse('the header does not state "p_match_bar"')
      if (lap_p_top_bar .le. 0.0d0)                                        &
         call refuse('the header does not state "p_top_bar"')
      if (lap_p_top_bar .ge. lap_p_match_bar)                              &
         call refuse('p_top_bar >= p_match_bar: the file stops at or below'&
                     //' the match, so the two models never overlap')
      if (lap_p_match_bar .gt. p_deepest .or.                              &
          lap_p_match_bar .lt. p_shallow)                                  &
         call refuse('p_match_bar lies outside the pressures the table'//  &
                     ' covers')
      if (lap_p_top_bar .lt. p_shallow)                                    &
         call refuse('p_top_bar lies above the shallowest level in the'//  &
                     ' table')
      if (len_trim(lap_solution_id) .eq. 0)                                &
         call refuse('the header does not state "solution_id"')

      lap_r_top_RJ = value_at_pressure(ic_r, lap_p_top_bar)
      end subroutine validate_table

      ! ------------------------------------------------------------------ !

      subroutine refuse(why)
      character(len=*), intent(in) :: why
      write(*,*) '(lower_atmosphere_profile) ERROR in '//trim(lap_file)//':'
      write(*,*) '  '//trim(why)//'.'
      write(*,*) '  Schema: docs/input_schema.md section 2d.'
      error stop 1
      end subroutine refuse

      ! ------------------------------------------------------------------ !

      real*8 function value_at_pressure(icol, p_bar)
      ! Linear in log p, never extrapolated (the caller has already checked
      ! that p_bar is covered).  A target that falls exactly on a level
      ! returns that level's value bit for bit, which is what makes a
      ! single-level profile reproduce the scalar handoff exactly.
      integer, intent(in) :: icol
      real*8,  intent(in) :: p_bar
      integer :: ilev, ilo
      real*8  :: wgt, xt
      ilo = 0
      do ilev = 1, nlev
         if (coltab(ilev, ic_p) .eq. p_bar) then
            value_at_pressure = coltab(ilev, icol)
            return
         endif
      enddo
      do ilev = 1, nlev-1
         if (p_bar .lt. coltab(ilev, ic_p) .and.                           &
             p_bar .gt. coltab(ilev+1, ic_p)) then
            ilo = ilev
            exit
         endif
      enddo
      if (ilo .eq. 0) then
         ! Outside the table: clamp to the nearer end rather than invent a
         ! value.  validate_table has already refused the cases where this
         ! would matter physically.
         if (p_bar .ge. coltab(1, ic_p)) then
            value_at_pressure = coltab(1, icol)
         else
            value_at_pressure = coltab(nlev, icol)
         endif
         return
      endif
      xt  = log(p_bar)
      wgt = (xt - logp(ilo))/(logp(ilo+1) - logp(ilo))
      value_at_pressure = coltab(ilo, icol)                                &
                          + wgt*(coltab(ilo+1, icol) - coltab(ilo, icol))
      end function value_at_pressure

      ! ------------------------------------------------------------------ !

      subroutine lap_value_at_match(nm, val, found)
      ! Interpolated value of a named column at the matching level.
      character(len=*), intent(in)  :: nm
      real*8,           intent(out) :: val
      logical,          intent(out) :: found
      integer :: icol
      val   = 0.0d0
      found = .false.
      if (.not. lap_in_use) return
      icol = column_index(nm)
      if (icol .eq. 0) return
      val   = value_at_pressure(icol, lap_p_match_bar)
      found = .true.
      end subroutine lap_value_at_match

      ! ------------------------------------------------------------------ !

      subroutine lap_element_ratio_at_match(elname, val, found)
      ! El/H nuclei ratio at the matching level, from the column "X_<El>".
      character(len=*), intent(in)  :: elname
      real*8,           intent(out) :: val
      logical,          intent(out) :: found
      call lap_value_at_match('X_'//trim(elname), val, found)
      end subroutine lap_element_ratio_at_match

      ! ------------------------------------------------------------------ !

      subroutine eddy_diffusion_on_grid
      ! Fill kzz_cell, the eddy diffusion coefficient of every cell, once the
      ! radial grid exists.
      !
      ! Without a profile every cell carries the scalar he_kzz, so the three
      ! sites in binary_element_diffusion that used to read he_kzz see the
      ! same number they saw before: for a uniform array the face average
      ! 0.5*(a+a) is exactly a in IEEE double, and D + K_face is then the same
      ! sum of the same two operands as D + he_kzz.
      !
      ! With a profile, K_zz is interpolated onto the grid through the
      ! profile's own radius column: the cell radius selects the bracketing
      ! levels, and the value is linear in log p between them.  Because the
      ! radius -> log p map inside a bracket is itself taken linearly, this is
      ! the same weight either way; the log-p statement is what fixes the
      ! choice at the level of the schema.  Above the shallowest level the
      ! file carries, the cell takes that level's value -- the eddy
      ! coefficient is a lower-atmosphere property and the file makes no
      ! statement above its top.
      integer :: jc, ilev, ilo
      real*8  :: r_RJ, wgt

      if (allocated(kzz_cell)) then
         if (size(kzz_cell) .ne. N + 2*Ng) deallocate(kzz_cell)
      endif
      if (.not. allocated(kzz_cell)) allocate(kzz_cell(1-Ng:N+Ng))

      if (.not. lap_in_use) then
         kzz_cell = he_kzz
         return
      endif

      do jc = 1-Ng, N+Ng
         r_RJ = r(jc)*R0/RJ
         ilo  = 0
         if (r_RJ .le. coltab(1, ic_r)) then
            kzz_cell(jc) = coltab(1, ic_kzz)
            cycle
         endif
         if (r_RJ .ge. coltab(nlev, ic_r)) then
            kzz_cell(jc) = coltab(nlev, ic_kzz)
            cycle
         endif
         do ilev = 1, nlev-1
            if (r_RJ .ge. coltab(ilev, ic_r) .and.                         &
                r_RJ .le. coltab(ilev+1, ic_r)) then
               ilo = ilev
               exit
            endif
         enddo
         if (ilo .eq. 0) then
            kzz_cell(jc) = coltab(nlev, ic_kzz)
            cycle
         endif
         wgt = (r_RJ - coltab(ilo, ic_r))                                  &
               /max(coltab(ilo+1, ic_r) - coltab(ilo, ic_r), 1.0d-300)
         kzz_cell(jc) = coltab(ilo, ic_kzz)                                &
                        + wgt*(coltab(ilo+1, ic_kzz) - coltab(ilo, ic_kzz))
      enddo

      write(*,'(A,ES10.3,A,ES10.3,A)')                                     &
         ' (lower_atmosphere_profile) K_zz(r) on the grid: ',              &
         minval(kzz_cell(1:N)), ' to ', maxval(kzz_cell(1:N)), ' cm^2/s'

      end subroutine eddy_diffusion_on_grid

      ! ------------------------------------------------------------------ !

      subroutine lap_report_provenance(iu)
      ! Provenance block for EXHALE_resolved.out, so the closure driver and
      ! src/utils/element_budget.py read one authority.
      integer, intent(in) :: iu
      write(iu,'(A,L1)')    'lower_profile_present     ', lap_in_use
      if (.not. lap_in_use) return
      write(iu,'(A)')       'lower_profile_file        '//trim(lap_file)
      write(iu,'(A)')       'lower_profile_solution_id '//trim(lap_solution_id)
      write(iu,'(A)')       'lower_profile_source_code '//trim(lap_source_code)
      write(iu,'(A)')       'lower_profile_source_ver  '//trim(lap_source_version)
      write(iu,'(A,ES23.15E3)') 'lower_profile_p_match_bar ', lap_p_match_bar
      write(iu,'(A,ES23.15E3)') 'lower_profile_p_top_bar   ', lap_p_top_bar
      write(iu,'(A,ES23.15E3)') 'lower_profile_trial_flux_H ', lap_trial_flux_H
      write(iu,'(A,ES23.15E3)') 'lower_profile_trial_flux_He', lap_trial_flux_He
      write(iu,'(A,I0)')    'lower_profile_iteration   ', lap_iteration
      write(iu,'(A,L1)')    'lower_profile_iterable    ', lap_iterable
      write(iu,'(A,L1)')    'lower_profile_steady      ', lap_steady
      write(iu,'(A)')       'lower_profile_notes       '//trim(lap_notes)
      end subroutine lap_report_provenance

      ! End of module
      end module lower_atmosphere_profile
