      module IC_load
      ! Module to load previous ICs, stored in the following files:
      ! - Hydro_ioniz_IC.txt
      ! - Ion_species_IC.txt
      !
      ! Two file generations are supported:
      !  * schema-2 files (written by write_output with '# ...' header
      !    lines): species columns are identified by the labels on the
      !    "# columns" line, so ALL species present in the file --
      !    including the metal ions -- are restored. Restarts therefore
      !    preserve the metal ionization state.
      !  * legacy headerless files: the original fixed 7-column read
      !    (r + H/He/HeITR); metals fall back to neutral-from-abundance
      !    (the historical behavior).

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,    &
                               isp_HeIII, isp_HeTR,                    &
                               n_mion, n_melem, mion_fsp, mion_name,   &
                               mion_elem, melem_i0, melem_top

      implicit none

      contains

      subroutine load_IC(rho,v,p,T,f_sp,W)

      ! Integer variables
      integer :: j, k, ios, nlab, c, e, i0

      ! Loaded number densities for every f_sp column (zero = not in file)
      real*8, dimension(1-Ng:N+Ng,n_species) :: nsp_l
      logical :: col_present(n_species), elem_ok
      ! Auxiliary temporary variable
      real*8 :: tmp
      real*8 :: vals(80)

      character(len=8192) :: line
      character(len=16)   :: labels(80)
      integer :: col2fsp(80)
      logical :: has_header

      ! Output variables
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: W


	   !-------------------------------------!

      ! Load thermodynamic variables (skip any '#' header lines)
      open(unit = 1, file = 'output/Hydro_ioniz_IC.txt')
         read(1,'(A)') line
         do while (is_comment(line))
            read(1,'(A)') line
         enddo
         read(line,*) tmp, tmp, v(1-Ng), p(1-Ng), T(1-Ng), tmp, tmp
         do j = 2-Ng,N+Ng
            read(1,*) tmp, tmp, v(j), p(j), T(j), tmp, tmp
         enddo
      close(1)

      ! Adimensionalize
      v = v/v0
      p = p/p0
      T = T/T0

      ! Load ionization profiles
      nsp_l       = 0.0d0
      col_present = .false.

      open(unit = 2, file = 'output/Ion_species_IC.txt')
      read(2,'(A)') line
      has_header = is_comment(line)

      if (has_header) then
         ! ---- schema-2 file: map columns by label ----
         nlab = 0
         do while (is_comment(line))
            if (index(line,'columns') .gt. 0) call parse_labels(line, labels, nlab)
            read(2,'(A)',iostat=ios) line
            if (ios .ne. 0) exit
         enddo
         if (nlab .lt. 2) stop '(load_IC) header found but no "# columns" line'

         ! Build the label -> f_sp column map (0 = ignore). labels(1) is r.
         col2fsp = 0
         do k = 2, nlab
            col2fsp(k) = species_column(labels(k))
            if (col2fsp(k) .gt. 0) col_present(col2fsp(k)) = .true.
         enddo

         ! First data record is already in 'line'
         read(line,*) (vals(k), k = 1,nlab)
         call scatter_row(1-Ng, vals, col2fsp, nlab, nsp_l)
         do j = 2-Ng,N+Ng
            read(2,*) (vals(k), k = 1,nlab)
            call scatter_row(j, vals, col2fsp, nlab, nsp_l)
         enddo

      else
         ! ---- legacy headerless file: fixed 7-column layout ----
         rewind(2)
         do j = 1-Ng,N+Ng
            read(2,*) r(j),                  &
                      nsp_l(j,isp_HI),       &
                      nsp_l(j,isp_HII),      &
                      nsp_l(j,isp_HeI),      &
                      nsp_l(j,isp_HeII),     &
                      nsp_l(j,isp_HeIII),    &
                      nsp_l(j,isp_HeTR)
         enddo
         col_present(isp_HI:isp_HeTR) = .true.
      endif
      close(2)

      ! Construct mass density profile (adimensional).
      ! Same H/He-mass formula as the historical loader (HeITR and trace
      ! metals excluded), so legacy reloads are bit-identical.
      rho = (nsp_l(:,isp_HI) + nsp_l(:,isp_HII)                        &
             + 4.0*(nsp_l(:,isp_HeI) + nsp_l(:,isp_HeII)               &
                    + nsp_l(:,isp_HeIII)))/n0

      ! H/He(+HeITR) fractions
      f_sp(:,isp_HI)    = nsp_l(:,isp_HI)/(rho*n0)
      f_sp(:,isp_HII)   = nsp_l(:,isp_HII)/(rho*n0)
      f_sp(:,isp_HeI)   = nsp_l(:,isp_HeI)/(rho*n0)
      f_sp(:,isp_HeII)  = nsp_l(:,isp_HeII)/(rho*n0)
      f_sp(:,isp_HeIII) = nsp_l(:,isp_HeIII)/(rho*n0)
      f_sp(:,isp_HeTR)  = nsp_l(:,isp_HeTR)/(rho*n0)

      ! Metals, element by element: if every ion stage of the element was
      ! present in the IC file, restore the loaded state; otherwise fall
      ! back to the historical neutral-from-abundance initialization
      ! (n_X/n_tot = X_X/(1+4*HeH), higher stages zero).
      do e = 1, n_melem
         i0 = melem_i0(e)
         elem_ok = .true.
         do k = 0, melem_top(e)
            if (.not. col_present(mion_fsp(i0+k))) elem_ok = .false.
         enddo
         if (elem_ok) then
            do k = 0, melem_top(e)
               c = mion_fsp(i0+k)
               f_sp(:,c) = nsp_l(:,c)/(rho*n0)
            enddo
         else
            do k = 0, melem_top(e)
               c = mion_fsp(i0+k)
               if (k .eq. 0) then
                  f_sp(:,c) = melem_ab(e)/(1.0 + 4.0*HeH)
               else
                  f_sp(:,c) = 0.0
               endif
            enddo
         endif
      enddo

      ! Construct matrix of primitive profiles
      W(:,1) = rho
      W(:,2) = v
      W(:,3) = p

      ! End of subroutine
      end subroutine load_IC

      !-------------------------------------!

      logical function is_comment(line)
      ! True if the (left-adjusted) line starts with '#'
      character(len=*), intent(in) :: line
      character(len=len(line)) :: t
      t = adjustl(line)
      is_comment = (len_trim(t) .gt. 0 .and. t(1:1) .eq. '#')
      end function is_comment

      !-------------------------------------!

      subroutine parse_labels(line, labels, nlab)
      ! Split a "# columns r[Rp] HI HII ..." line into its column labels
      ! ('#' and 'columns' tokens are dropped; labels(1) is the r column).
      character(len=*), intent(in)  :: line
      character(len=16), intent(out) :: labels(:)
      integer, intent(out) :: nlab
      integer :: pos, l, start
      character(len=len(line)) :: rest
      character(len=64) :: tok

      rest = adjustl(line)
      nlab = 0
      do
         l = len_trim(rest)
         if (l .eq. 0) exit
         pos = index(rest, ' ')
         if (pos .le. 1) then
            tok = rest(1:l); rest = ''
         else
            tok = rest(1:pos-1); rest = adjustl(rest(pos:))
         endif
         if (trim(tok) .eq. '#')        cycle
         if (trim(tok) .eq. 'columns')  cycle
         nlab = nlab + 1
         labels(nlab) = trim(tok)
      enddo
      end subroutine parse_labels

      !-------------------------------------!

      integer function species_column(label)
      ! f_sp column for a species label (0 = unknown/ignored, e.g. r[Rp])
      character(len=*), intent(in) :: label
      integer :: i
      select case (trim(label))
         case ('HI');    species_column = isp_HI
         case ('HII');   species_column = isp_HII
         case ('HeI');   species_column = isp_HeI
         case ('HeII');  species_column = isp_HeII
         case ('HeIII'); species_column = isp_HeIII
         case ('HeITR'); species_column = isp_HeTR
         case default
            species_column = 0
            do i = 1, n_mion
               if (trim(label) .eq. trim(mion_name(i))) then
                  species_column = mion_fsp(i)
                  return
               endif
            enddo
      end select
      end function species_column

      !-------------------------------------!

      subroutine scatter_row(j, vals, col2fsp, nlab, nsp_l)
      ! Store one data record: vals(1) is r, vals(k>=2) go to their
      ! mapped f_sp columns.
      integer, intent(in) :: j, nlab
      real*8,  intent(in) :: vals(:)
      integer, intent(in) :: col2fsp(:)
      real*8,  intent(inout) :: nsp_l(1-Ng:N+Ng, n_species)
      integer :: k
      r(j) = vals(1)
      do k = 2, nlab
         if (col2fsp(k) .gt. 0) nsp_l(j,col2fsp(k)) = vals(k)
      enddo
      end subroutine scatter_row

      ! End of module
      end module IC_load
