      module metals_input
      ! Optional metals.inp reader: set trace-metal abundances at runtime
      ! (no recompile). Format is one ion per line:
      !
      !     <ION>  <abundance n_X/n_H by number>
      !
      ! Blank lines and lines starting with '#' are ignored. Recognized
      ! ions (neutral label = element): CI -> X_C, NI -> X_N, OI -> X_O,
      ! MgI -> X_Mg, SiI -> X_Si, CaI -> X_Ca, NaI -> X_Na, KI -> X_K,
      ! SI -> X_S, FeI -> X_Fe (the bare element symbols
      ! C/N/O/Mg/Si/Ca/Na/K/S/Fe are also accepted). A missing metals.inp
      ! leaves the X_* defaults from
      ! input_read.f90 untouched, so the H/He-only default is preserved.
      !
      ! Mirrors the ATES_extended metals.inp format so the two trees
      ! share input files.

      use global_parameters
      use charge_exchange, only: cx_full,  &  ! full-Table-4 toggle
                                 cx_o2p_h_scale  ! Group-E O2+ + H0 dial

      implicit none

      private
      public :: read_metals_input

      character(len = *), parameter :: met_inp_file = 'metals.inp'

      contains

      !----------------------------------------------------------!

      subroutine read_metals_input
      integer :: io, n_set
      character(len = 256) :: line
      character(len = :), allocatable :: trimmed, tok, val
      integer :: sp_pos
      real*8  :: ab
      logical :: file_exists

      inquire(file = met_inp_file, exist = file_exists)
      if (.not. file_exists) then
         write(*,*) '(metals_input) No metals.inp found; using built-in', &
                    ' X_C/X_N/X_O/X_Mg defaults.'
         return
      endif

      write(*,*) '(metals_input) Reading metals.inp..'
      open(unit = 36, file = met_inp_file, status = 'old', action = 'read')

      n_set = 0
      do
         read(36, '(A)', iostat = io) line
         if (io /= 0) exit

         trimmed = trim(adjustl(line))
         if (len(trimmed) == 0) cycle
         if (trimmed(1:1) == '#') cycle

         ! Split on first whitespace: <token> <value>
         sp_pos = scan(trimmed, ' '//achar(9))
         if (sp_pos <= 1) then
            write(*,*) '  WARN: skipping malformed line: ', trim(trimmed)
            cycle
         endif

         tok = trim(adjustl(trimmed(1:sp_pos-1)))
         val = trim(adjustl(trimmed(sp_pos+1:)))
         read(val, *, iostat = io) ab
         if (io /= 0) then
            write(*,*) '  WARN: bad abundance on line: ', trim(trimmed)
            cycle
         endif

         ! Labels are matched CASE-SENSITIVELY against the canonical
         ! chemical symbols, because the symbols themselves are case
         ! sensitive: 'Si' vs 'S', 'Na' vs 'N', 'Ca' vs 'C'. Each metal
         ! accepts its neutral ion label (e.g. 'SiI') and its bare symbol
         ! (e.g. 'Si'). The single-letter sulphur is 'S'/'SI'; silicon is
         ! 'Si'/'SiI'.
         ! 'cx_full <0|1>' toggles the full Huang Table 4 (He-H and
         ! metal-metal charge exchange) instead of the metal-H default.
         if (trim(tok) == 'cx_full' .or. trim(tok) == 'CX_FULL') then
            cx_full = (ab > 0.5d0)
            write(*,'(a,l1)') '   charge-exchange full Table 4 mode = ', &
                              cx_full
            cycle
         endif

         ! 'cx_O2p_H <scale>' rescales the ONE charge-transfer reaction that
         ! is absent from Huang Table 4 and whose absence leaves the O III
         ! profile wrong by decades, O2+ + H0 -> O+ + H+, in units of the
         ! published Barragan et al. (2006) rate. 1 (default) is that rate;
         ! 0 reproduces the Table-4-only reaction set exactly; other values
         ! are bounding experiments. It is deliberately NOT folded into
         ! cx_full, which means "all of Table 4".
         if (trim(tok) == 'cx_O2p_H' .or. trim(tok) == 'cx_o2p_h') then
            cx_o2p_h_scale = max(ab, 0.0d0)
            if (cx_o2p_h_scale > 0.0d0) then
               write(*,'(a,es9.2,a)') '   O2+ + H0 -> O+ + H+ charge '     &
                  // 'transfer ON at ', cx_o2p_h_scale,                    &
                  ' x Barragan+2006 [default 1]'
            else
               write(*,'(a)') '   O2+ + H0 charge transfer OFF '           &
                  // '(Huang Table 4 only)'
            endif
            cycle
         endif

         ! 'cno_cool <0|1>' selects the C/N/O line-cooling source:
         ! 1 (default) = CHIANTI v11 closed-form fits including N I/N II,
         ! 0 = legacy AIOLOS analytic fits (no N cooling; O deviates
         !     40-70% from CHIANTI in the wind region).
         if (trim(tok) == 'cno_cool' .or. trim(tok) == 'CNO_COOL') then
            cno_chianti = (ab > 0.5d0)
            if (cno_chianti) then
               write(*,'(a)') '   C/N/O cooling source = CHIANTI v11' &
                              // ' fits (incl. N I/N II) [default]'
            else
               write(*,'(a)') '   C/N/O cooling source = legacy AIOLOS' &
                              // ' fits (no N cooling)'
            endif
            cycle
         endif

         ! 'eos_metals <0|1>' selects the bulk-gas metal policy:
         ! 1 (default) = metals contribute their mass (rho_bc, calc_rho,
         !     hence mu/v0/p0/b0), electrons (calc_ne, dp_bc) and nuclei
         !     (calc_ntot, ghost pressure) to the gas budget;
         ! 0 = legacy trace approximation (H/He-only budget).
         if (trim(tok) == 'eos_metals' .or. trim(tok) == 'EOS_METALS') then
            eos_include_metals = (ab > 0.5d0)
            if (eos_include_metals) then
               write(*,'(a)') '   EOS metal policy = metals in mass/'   &
                              // 'electron/particle budget [default]'
            else
               write(*,'(a)') '   EOS metal policy = legacy trace'      &
                              // ' approximation (H/He-only budget)'
            endif
            cycle
         endif

         ! 'pp_metals <0|1|2>' selects how metals are treated in the
         ! advection post-process (post_process_adv): 0 metal-free,
         ! 1 frozen eq metals (default), 2 re-solve. 'pp_metal_mode' is
         ! accepted as a synonym.
         if (trim(tok) == 'pp_metals' .or. trim(tok) == 'pp_metal_mode') then
            pp_metal_mode = nint(ab)
            if (pp_metal_mode < 0 .or. pp_metal_mode > 2) then
               write(*,'(a,i0,a)') '  WARN: pp_metals = ', pp_metal_mode, &
                  ' out of range [0,2]; clamping.'
               pp_metal_mode = min(max(pp_metal_mode, 0), 2)
            endif
            write(*,'(a,i0,a)') '   post-process metal mode pp_metals = ', &
               pp_metal_mode, &
               ' (0=metal-free, 1=frozen, 2=re-solve)'
            cycle
         endif

         select case (trim(tok))
            case ('CI',  'C');   X_C  = ab
            case ('NI',  'N');   X_N  = ab
            case ('OI',  'O');   X_O  = ab
            case ('MgI', 'Mg');  X_Mg = ab
            case ('SiI', 'Si');  X_Si = ab
            case ('CaI', 'Ca');  X_Ca = ab
            case ('NaI', 'Na');  X_Na = ab
            case ('KI',  'K');   X_K  = ab
            case ('SI',  'S');   X_S  = ab
            case ('FeI', 'Fe');  X_Fe = ab
            case default
               write(*,*) '  WARN: unrecognized ion ', &
                          '(only C/N/O/Mg/Si/Ca/Na/K/S/Fe, case sensitive): ', &
                          trim(tok)
               cycle
         end select

         n_set = n_set + 1
      enddo

      close(36)

      write(*,'(a,i0,a)') ' (metals_input) Done. ', n_set, &
                          ' metal abundance(s) set.'
      write(*,'(a,4es11.3)') '   X_C, X_N, X_O, X_Mg     = ', &
                             X_C, X_N, X_O, X_Mg
      write(*,'(a,5es11.3)') '   X_Si, X_Ca, X_Na, X_K, X_S = ', &
                             X_Si, X_Ca, X_Na, X_K, X_S
      write(*,'(a,es11.3)')  '   X_Fe                     = ', X_Fe

      end subroutine read_metals_input

      ! End of module
      end module metals_input
