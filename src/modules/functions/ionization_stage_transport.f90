      module ionization_stage_transport
      ! THE TRANSPORT OF THE IONIZATION STAGES OF ONE ELEMENT, WRITTEN ON
      ! THAT ELEMENT'S OWN NUCLEUS FLUX.
      !
      ! WHAT IS CARRIED.  For an element with stages k = 1 ... n_k + 1, the
      ! fractions of that element's NUCLEI in each stage, x_k = n_k/n_el.
      ! The n_k carried rows are the ionized stages; the neutral stage is the
      ! CLOSING one and is not a row, it is one minus the sum of the others.
      ! Hydrogen carries x(H II) and closes on x(H I); helium carries
      ! x(He II) and x(He III) and closes on x(He I).  The He 2^3S level is a
      ! sublevel inside He I and not a stage of this partition
      ! (species_table.f90 bsp_is_excited_level), so it is not carried.
      !
      ! THE FACE FLUX.  At the face f between cells j and j+1,
      !
      !    F_k(f) = x_k(f) N_el(f)
      !           - n_el(f) K(f) [ x_k(j+1) - x_k(j) ] / dr(f)          (1)
      !
      ! with N_el(f) the element NUCLEUS face flux the element operator forms
      ! (advection plus gradient plus eddy plus settling drift, one object,
      ! binary_element_diffusion.f90 element_nucleus_face_flux), x_k(f) the
      ! face fraction, n_el(f) the face element nucleus density, K(f) the
      ! face eddy coefficient and dr(f) the face spacing.  The second term is
      ! the stage's own eddy term and nothing else: the element's eddy flux
      ! is already inside N_el(f), and adding the whole mixing-ratio eddy
      ! flux -n_tot K d(x_k n_el/n_tot)/dr inside every stage would carry
      ! x_k J_el,eddy twice.
      !
      ! THE IDENTITY.  Summed over ALL stages, neutral included,
      !
      !    sum_k F_k(f) = N_el(f)     to rounding at every face,          (2)
      !
      ! which says that moving charge between stages cannot move a nucleus.
      ! It holds under three conditions:
      !
      !   C1  ONE element face flux multiplies every stage.  A second copy,
      !       rebuilt per stage from its own donor rule, does not return one
      !       element flux.  This module never rebuilds N_el; the caller
      !       reads it from the element operator.
      !   C2  the face fractions close the simplex exactly: the carried
      !       stages are reconstructed and the CLOSING stage is one minus
      !       their sum, never reconstructed on its own.  The donor side is
      !       taken from the sign of N_el and not from any stage's own
      !       gradient, so every stage leaves the cell its element leaves.
      !   C3  n_el(f), K(f) and dr(f) are ONE face value shared by every
      !       stage, and the closing stage's eddy term is MINUS the sum of
      !       the others'.  Forming it from the closing stage's own cell
      !       values is exact in exact arithmetic and loses the gradient to
      !       the rounding of 1 - sum wherever the carried stages are trace
      !       (MEASURED at 2.9e-3 of the term).
      !
      ! THE BOUNDARY FACES.  The element operator carries no diffusive flux
      ! through f = 0 or f = N: cell 1 is the Dirichlet reservoir the gas
      ! enters from and the top face is an outflow whose ghost repeats cell
      ! N.  A stage eddy term at either face would move charge
      ! through a face the element itself cannot cross, so the second term of
      ! (1) is zero there and the boundary face carries x_k(f) N_el(f) alone.
      ! Which cell that x_k(f) comes from is the sign of N_el(f) and nothing
      ! else: with N_el(0) > 0 the reservoir's partition enters, with
      ! N_el(0) < 0 (the breathing base the code admits) the domain's own
      ! leaves; at the top, N_el(N) >= 0 is outflow needing no data and
      ! N_el(N) < 0 brings the outermost cell's own composition back in.
      !
      ! THE DIVERGENCE.  On the one spherical geometry of this grid
      ! (grid_construction spherical_face_area_and_cell_volume),
      !
      !    D_k(j) = [ A_+ F_k(j) - A_- F_k(j-1) ] / V_j                  (3)
      !
      ! in cm^-3 s^-1, with A in cm^2 and V in cm^3.  Because (2) holds face
      ! by face and (3) telescopes on one geometry, the sum over stages of
      ! the cell divergences is the element's own nucleus divergence and the
      ! column sum of V_j sum_k D_k(j) is the difference of the two boundary
      ! nucleus fluxes alone.  Both are acceptance rows.
      !
      ! THE ONE-VELOCITY APPROXIMATION AND ITS VALIDITY RANGE.  Every stage
      ! of an element is given the element's velocity: no drift of an ion
      ! against its own neutral is carried.  What decides it is the drift
      ! w = D_in dln(n_e T)/dr an ion acquires against its neutral in the
      ! ambipolar field, measured against the bulk velocity.
      !
      ! MEASURED on the certified atomic fiducial LHS1140b/models/.L26/
      ! fid_resolve (LHS 1140 b, He_Kzz = 1e9, He_diffusion on, 500 cells to
      ! 29 R_p), with the code's own non-resonant ion-neutral coefficients
      ! and the resonant charge-exchange channel they exclude added at
      ! sigma = 2e-15 cm^2 (docs/collisional_validity.md section 1):
      !
      !   |w/v| stays below 0.1 only inside r = 1.22 R_p (H) and 1.24 R_p
      !   (He); over r >= 1.2 R_p it reaches 0.62 (H) and 0.43 (He) without
      !   the resonant channel and 0.52 and 0.34 with it.  The momentum
      !   transfer itself is fast: t_coll/t_flow is below 0.1 inside
      !   14.5 R_p for hydrogen and 8.2 R_p for helium, so the stages are
      !   collisionally locked in the sense of time scales and still drift,
      !   because the drift is set by the ambipolar field against the
      !   friction and not by the friction alone.
      !
      ! So the approximation is exact only inside about 1.2 R_p, which is
      ! also where the Damkohler number returns the local ionization root
      ! whatever the transport does, and outward of it the neglected drift
      ! is a systematic of order 0.02 to 0.6 of the advective stage flux --
      ! not smaller than the stage transport this module exists to carry.
      ! Above about 13 R_p the continuum stage equation is itself
      ! unvalidated by the run's own Knudsen measure (max Kn 0.88 in the
      ! heating region on this state), and the drift is not the leading
      ! approximation there.  An ambipolar stage drift is a separate physical
      ! term, not a refinement of this one.
      !
      ! RELATION TO THE PROTON CARRIER.  The transported proton of
      ! "Ionization transport: True" (diffusive_photochemistry.f90, carrier
      ! ic_Hp) is the same physical object as the hydrogen row here, written
      ! in the other variable: a fraction per unit MASS carried on the BULK
      ! face mass flux, with the whole mixing-ratio eddy flux inside it.
      ! That is correct exactly while the bulk mass flux carries no eddy term
      ! of its own, and it stops being correct the moment the advective half
      ! becomes the element nucleus flux, which does.  This module is the
      ! replacement, not a second transport:
      ! there is to be ONE spelling of the stage flux, and the carrier
      ! operator's ionization rows are to be formed from (1) here.

      use global_parameters
      use grid_construction, only: spherical_face_area_and_cell_volume
      use species_advective_transport, only: species_face_fraction
      use binary_element_diffusion, only: element_nucleus_face_flux,      &
                                          element_nucleus_counts,         &
                                          m_He_amu

      implicit none
      private

      ! The face fractions and the face fluxes of one element's stages, and
      ! the cell divergence of those fluxes on the one spherical geometry.
      public :: ionization_stage_face_fractions
      public :: ionization_stage_face_flux
      public :: ionization_stage_divergence
      ! What identity (2) can stand at in floating point, and the largest
      ! value the ratio of magnitude sums in it can take. The identity is
      ! algebraic, so this is the whole of what a gate on it clears.
      public :: ionization_stage_sum_rounding_bound
      public :: ionization_stage_sum_scale_ratio_limit
      ! The derivative of the face flux and of the cell row with respect to
      ! the carried cell fractions.  The Jacobian is formed from the SAME
      ! coefficients the flux is, so a residual change cannot leave the rows
      ! behind.
      public :: ionization_stage_face_jacobian
      public :: ionization_stage_row_jacobian
      ! The admissible set of the carried fractions is the simplex, not a
      ! box: an ionization fraction is bounded by the sum rule and not by a
      ! free nucleus density.
      public :: stage_simplex_projection
      ! The nucleus face flux of both elements, their face nucleus densities
      ! and the face eddy coefficient and spacing the stage rows are written
      ! against, all from the element operator's one public flux.
      public :: hydrogen_and_helium_nucleus_face_flux

      ! How far the last projection had to move the carried fractions:
      ! the largest amount by which their sum exceeded the nuclei available
      ! to them, and the largest negative fraction it met, both before the
      ! projection.  Negative means the state was already inside the
      ! simplex.  MEASURED and reported, in the shape he_fraction_over_one
      ! and species_face_excursion already have, rather than asserted.
      real*8, protected, public :: stage_simplex_sum_over_limit = -1.0d0
      real*8, protected, public :: stage_fraction_under_zero    = -1.0d0

      contains

      ! ------------------------------------------------------------------ !

      subroutine ionization_stage_face_fractions(xion, Nel, xf, xclose)
      ! THE FACE FRACTIONS OF ONE ELEMENT'S STAGES, closing the simplex at
      ! the RECONSTRUCTED face and not only at the cell centres (condition
      ! C2 of the header).
      !
      ! Every carried stage goes through species_face_fraction, the same
      ! reconstruction, limiter and donor rule the element mass fractions
      ! take; the donor side is selected by the ELEMENT nucleus flux Nel, so
      ! one stage cannot be donated from one side of a face while another is
      ! donated from the other.  The closing stage is then one minus their
      ! sum, which makes sum_k x_k(f) = 1 to one subtraction's rounding
      ! (MEASURED 2.2e-16).
      !
      ! What this rule buys is the SUM and nothing else.  The carried stages
      ! come back inside [0,1] because species_face_fraction scales them
      ! there; the closing stage is a difference and can still leave it
      ! (MEASURED largest negative face value -8.6e-4 on the simplex-face
      ! column).  Bounding the closing face value as well would break (2),
      ! so the excursion is measured at the CELL state by
      ! stage_simplex_projection and reported, never repaired at the face.
      real*8, dimension(:,1-Ng:),   intent(in)  :: xion
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Nel
      real*8, dimension(:,0:),      intent(out) :: xf
      real*8, dimension(0:N),       intent(out) :: xclose

      real*8, dimension(1-Ng:N+Ng) :: work, wf
      integer :: k, nk, j

      nk     = size(xion,1)
      xclose = 1.0d0
      do k = 1, nk
         work = xion(k,:)
         call species_face_fraction(work, Nel, wf)
         do j = 0, N
            xf(k,j) = wf(j)
            xclose(j) = xclose(j) - wf(j)
         enddo
      enddo

      end subroutine ionization_stage_face_fractions

      ! ------------------------------------------------------------------ !

      subroutine ionization_stage_face_flux(xion, Nel, n_elf, Kf, drf,    &
                                            xf, xclose, Fk, Fclose,       &
                                            donor_value, Ek_out)
      ! EQUATION (1) OF THE HEADER, at the faces f = 0 ... N, in nuclei
      ! cm^-2 s^-1 when Nel is in nuclei cm^-2 s^-1, n_elf in cm^-3, Kf in
      ! cm^2 s^-1 and drf in cm.
      !
      ! Fk(k,f) is the carried stage k, Fclose(f) the closing stage.  The
      ! closing stage's eddy term is MINUS the sum of the carried ones'
      ! (condition C3), and the shared face quantities n_elf, Kf, drf are
      ! used as given: a donor-cell density chosen from each stage's own
      ! gradient breaks the telescoping outright (MEASURED 2.0e-2).
      !
      ! The eddy term is zero at f = 0 and f = N (the header's boundary
      ! paragraph).  No width floor is applied to drf: a zero-width cell is
      ! a grid error the grid constructor refuses.
      !
      ! donor_value = .true. replaces the reconstruction by the donor cell's
      ! own value, which is the face the Jacobian below differentiates.  The
      ! identity (2) holds in both, because it is a statement about the
      ! stages sharing one face rule and not about which rule that is.
      !
      ! Ek_out returns the stage eddy terms themselves, the second term of
      ! (1).  They are what the rounding bound of identity (2) is written
      ! against (ionization_stage_sum_rounding_bound), and they are returned
      ! rather than re-formed by the caller so that the eddy term has ONE
      ! spelling, the header's rule for the stage flux.
      real*8, dimension(:,1-Ng:),   intent(in)  :: xion
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Nel
      real*8, dimension(0:N),       intent(in)  :: n_elf, Kf, drf
      real*8, dimension(:,0:),      intent(out) :: xf
      real*8, dimension(0:N),       intent(out) :: xclose
      real*8, dimension(:,0:),      intent(out) :: Fk
      real*8, dimension(0:N),       intent(out) :: Fclose
      logical, optional,            intent(in)  :: donor_value
      real*8, dimension(:,0:), optional, intent(out) :: Ek_out

      real*8  :: nKdr, Ek
      integer :: j, k, nk, jd
      logical :: donor

      nk    = size(xion,1)
      donor = .false.
      if (present(donor_value)) donor = donor_value

      if (donor) then
         xclose = 1.0d0
         do j = 0, N
            jd = j
            if (Nel(j) .lt. 0.0d0) jd = min(j+1, N+Ng)
            do k = 1, nk
               xf(k,j)   = xion(k,jd)
               xclose(j) = xclose(j) - xion(k,jd)
            enddo
         enddo
      else
         call ionization_stage_face_fractions(xion, Nel, xf, xclose)
      endif

      Fk     = 0.0d0
      Fclose = 0.0d0
      if (present(Ek_out)) Ek_out = 0.0d0
      do j = 0, N
         if (j .eq. 0 .or. j .eq. N) then
            nKdr = 0.0d0
         else
            nKdr = n_elf(j)*Kf(j)/drf(j)
         endif
         Fclose(j) = xclose(j)*Nel(j)
         do k = 1, nk
            Ek        = -nKdr*(xion(k,j+1) - xion(k,j))
            Fk(k,j)   = xf(k,j)*Nel(j) + Ek
            Fclose(j) = Fclose(j) - Ek
            if (present(Ek_out)) Ek_out(k,j) = Ek
         enddo
      enddo

      end subroutine ionization_stage_face_flux

      ! ------------------------------------------------------------------ !

      function ionization_stage_sum_scale_ratio_limit(nk) result(g)
      ! THE LARGEST VALUE g = S/S' CAN TAKE at a face of an element carrying
      ! nk stages, with
      !
      !    S  = |N| + |F_close| + sum_k |F_k| + sum_k |E_k|                (4)
      !    S' = |N| + |F_close| + sum_k |F_k|
      !
      ! the magnitude sum that bounds the rounding of identity (2) and the
      ! magnitude sum the measure of it divides by.  E_k is the stage's own
      ! eddy term of (1); it enters S because the same E_k is added to a
      ! carried stage and subtracted from the closing one, so its own
      ! rounding cancels in the sum and only its magnitude leaks.
      !
      ! DERIVED.  With the face fractions in the simplex,
      ! |x_close| + sum_k x_k(f) = 1, so term by term
      !
      !    |F_close| + sum_k |F_k|  >=  |sum_k E_k| + sum_k |E_k| - |N| ,
      !
      ! hence S' >= |sum_k E_k| + sum_k |E_k| and
      !
      !    g = 1 + sum_k |E_k|/S'
      !      <= 1 + sum_k |E_k| / ( sum_k |E_k| + |sum_k E_k| ) .
      !
      ! ONE carried stage has |sum_k E_k| = sum_k |E_k| and reaches 3/2; two
      ! or more can cancel their eddy terms against one another, and reach 2
      ! in that limit.  The carried face fractions are in [0,1] by
      ! construction (ionization_stage_face_fractions); the closing one is a
      ! difference and leaves it by at most -8.6e-4 (MEASURED, that routine),
      ! which raises the right-hand side by that much of itself.
      !
      ! MEASURED over manufactured columns whose eddy coefficient spans
      ! K_0 = 0 to 1e15 cm^2 s^-1: 1.00 to 1.50.
      integer, intent(in) :: nk
      real*8  :: g
      if (nk .le. 1) then
         g = 1.5d0
      else
         g = 2.0d0
      endif
      end function ionization_stage_sum_scale_ratio_limit

      ! ------------------------------------------------------------------ !

      function ionization_stage_sum_rounding_bound(nk, g) result(dbound)
      ! WHAT IDENTITY (2) CAN STAND AT IN FLOATING POINT for an element
      ! carrying nk stages, in the units the identity is measured in:
      ! relative to max(|N_el(f)|, |F_close(f)| + sum_k |F_k(f)|) at the
      ! face.  This is the one place the number is formed; the certification
      ! gates the identity against it.
      !
      ! The identity is ALGEBRAIC and carries no truncation error, so what a
      ! tolerance on it clears is the rounding of the sums alone.  On the
      ! path from the carried cell fractions to the measured difference
      ! there are 5 nk + 2 rounding events:
      !
      !    nk   forming x_close = 1 - sum_k x_k
      !     1   the product x_close N
      !    nk   subtracting each eddy term from the closing flux
      !   2nk   each carried stage's product and its sum with its eddy term
      !    nk   summing the nk + 1 fluxes in the measure
      !     1   the final difference
      !
      ! each bounded by eps times a quantity no larger than S of (4), and the
      ! measure divides by max(|N|, S' - |N|) >= S'/2, so
      !
      !    d  <=  2 (5 nk + 2) eps g ,      g = S/S' .                     (5)
      !
      ! eps = epsilon(1.0d0) is conservative by two, one rounding being
      ! bounded by the unit roundoff eps/2.
      !
      ! g IS MEASURED at the face the identity was measured at and passed in
      ! (ionization_stage_nucleus_sum returns it beside the measure), so the
      ! bound is the rounding of that face and not of the worst face the
      ! construction admits.  Without it the value returned is the bound at
      ! the algebraic ceiling of g, 4.663e-15 for hydrogen (nk = 1) and
      ! 1.066e-14 for helium (nk = 2), which is what a caller with no
      ! measured g gets.
      !
      ! VALIDATED BY EXECUTION.  On manufactured columns spanning K_0 = 0 to
      ! 1e15 cm^2 s^-1 the largest measure reaches 0.0703 (hydrogen) and
      ! 0.0580 (helium) of the bound at the face's own g, and the production
      ! states reach 0.06 to 0.08 of it (their worst faces carry g = 1 to
      ! five decimals, the advective nucleus flux dominating the eddy term
      ! there); a broken construction of the same identity stands at 1.5e-4
      ! to 3.6e-4, ten decades higher.
      integer,          intent(in) :: nk
      real*8, optional, intent(in) :: g
      real*8 :: dbound, gg
      gg = ionization_stage_sum_scale_ratio_limit(nk)
      if (present(g)) gg = g
      dbound = 2.0d0*dble(5*nk + 2)*epsilon(1.0d0)*gg
      end function ionization_stage_sum_rounding_bound

      ! ------------------------------------------------------------------ !

      subroutine ionization_stage_divergence(Fk, Fclose, Dk, Dclose)
      ! EQUATION (3) OF THE HEADER: the cell divergence of the stage face
      ! fluxes on the one spherical geometry of this grid, in cm^-3 s^-1
      ! for a face flux in cm^-2 s^-1.
      !
      ! The face areas and the shell volumes come from
      ! spherical_face_area_and_cell_volume, which returns them in the
      ! radius unit of r and r_edg; they are put into cm here by R0^2 and
      ! R0^3.  Every contribution to one conserved quantity divides by THIS
      ! volume, which is what makes the internal faces cancel in the column
      ! sum.
      real*8, dimension(:,0:), intent(in)  :: Fk
      real*8, dimension(0:N),  intent(in)  :: Fclose
      real*8, dimension(:,1:), intent(out) :: Dk
      real*8, dimension(1:N),  intent(out) :: Dclose

      real*8, dimension(0:N) :: fa
      real*8, dimension(1:N) :: cv
      real*8  :: R0sq, R0cb, sL, sR, Vj
      integer :: j, k, nk

      nk = size(Fk,1)
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0

      do j = 1, N
         sL = fa(j-1)*R0sq
         sR = fa(j)  *R0sq
         Vj = cv(j)*R0cb
         do k = 1, nk
            Dk(k,j) = (sR*Fk(k,j) - sL*Fk(k,j-1))/Vj
         enddo
         Dclose(j) = (sR*Fclose(j) - sL*Fclose(j-1))/Vj
      enddo

      end subroutine ionization_stage_divergence

      ! ------------------------------------------------------------------ !

      subroutine ionization_stage_face_jacobian(Nel, n_elf, Kf, drf,      &
                                                dFdl, dFdr)
      ! THE DERIVATIVE OF THE FACE FLUX (1) WITH RESPECT TO THE CARRIED CELL
      ! FRACTIONS, taken at the donor-cell face value:
      !
      !    dF_k(f)/dx_k(j)   =  N_el(f) [N_el(f) >= 0]  +  n_el K/dr
      !    dF_k(f)/dx_k(j+1) =  N_el(f) [N_el(f) <  0]  -  n_el K/dr
      !
      ! It is the same for every carried stage, because F_k depends on x_k
      ! alone: the closing stage is where the stages couple, and the closing
      ! stage is not a row.
      !
      ! The limiter and the second-order part of the reconstruction are left
      ! out, exactly as element_advective_face_coefficients leaves them out
      ! of the element row: they are the strongly nonlinear part of the
      ! operator and reach beyond the tridiagonal band.  The residual is the
      ! full operator either way, so this changes what a Newton converges AT,
      ! not what it converges TO.  The acceptance row differentiates the
      ! donor-value face and measures the second-order part separately.
      !
      ! The sign pattern is the M-matrix one: the eddy half puts a positive
      ! entry on the left cell and a negative one on the right, and the
      ! advective half puts the flux on its donor side alone.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Nel
      real*8, dimension(0:N),       intent(in)  :: n_elf, Kf, drf
      real*8, dimension(0:N),       intent(out) :: dFdl, dFdr

      real*8  :: nKdr
      integer :: j

      do j = 0, N
         if (j .eq. 0 .or. j .eq. N) then
            nKdr = 0.0d0
         else
            nKdr = n_elf(j)*Kf(j)/drf(j)
         endif
         if (Nel(j) .ge. 0.0d0) then
            dFdl(j) = Nel(j) + nKdr
            dFdr(j) =        - nKdr
         else
            dFdl(j) =          nKdr
            dFdr(j) = Nel(j) - nKdr
         endif
      enddo

      end subroutine ionization_stage_face_jacobian

      ! ------------------------------------------------------------------ !

      subroutine ionization_stage_row_jacobian(dFdl, dFdr, aa, bb, cc)
      ! THE TRIDIAGONAL ROW OF THE STAGE DIVERGENCE (3), from the face
      ! derivatives above and the SAME areas and volume the divergence
      ! divides by:
      !
      !    aa(j) = d D(j)/d x(j-1) = -A_- dF(j-1)/dx(j-1) / V_j
      !    bb(j) = d D(j)/d x(j)   = [A_+ dF(j)/dx(j)
      !                               - A_- dF(j-1)/dx(j)] / V_j
      !    cc(j) = d D(j)/d x(j+1) =  A_+ dF(j)/dx(j+1) / V_j
      !
      ! The time term of a backward-Euler step, n_el(j)/dt, is the caller's
      ! and is added to bb by the caller: this routine is the derivative of
      ! the transport alone, so that one change of the flux coefficients
      ! moves the residual and these rows together.
      real*8, dimension(0:N), intent(in)  :: dFdl, dFdr
      real*8, dimension(1:N), intent(out) :: aa, bb, cc

      real*8, dimension(0:N) :: fa
      real*8, dimension(1:N) :: cv
      real*8  :: R0sq, R0cb, sL, sR, Vj
      integer :: j

      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0

      do j = 1, N
         sL = fa(j-1)*R0sq
         sR = fa(j)  *R0sq
         Vj = cv(j)*R0cb
         aa(j) = -sL*dFdl(j-1)/Vj
         bb(j) = ( sR*dFdl(j) - sL*dFdr(j-1))/Vj
         cc(j) =   sR*dFdr(j)/Vj
      enddo

      end subroutine ionization_stage_row_jacobian

      ! ------------------------------------------------------------------ !

      subroutine stage_simplex_projection(xion, xsum_max)
      ! THE ADMISSIBLE SET OF THE CARRIED FRACTIONS IS THE SIMPLEX.  Each
      ! x_k lies in [0,1] and they sum to at most the share of the element
      ! the carried stages may hold, the remainder being the closing stage,
      ! so the element headroom of a free density does not apply: an
      ! ionization fraction is bounded by the sum rule.
      !
      ! THAT SHARE IS ONE ONLY WHERE EVERY NUCLEUS OF THE ELEMENT IS IN AN
      ! ATOMIC STAGE.  A nucleus held in a molecule is not partitioned by
      ! the stage rows and is not theirs to take, and neither is a nucleus
      ! in a level this step holds frozen: for helium the two are the HeH+
      ! bond and the He 2^3S metastable, whose reservation makes the bound
      ! x(He II) + x(He III) <= 1 - (n_HeH+ + n_He(2^3S))/n_He,nuc
      ! (carrier_helium_available_to_stages, which forms the numerator from
      ! the species table).  Without that reservation the simplex admits a
      ! partition whose neutral remainder is negative, and the element
      ! inventory of the state written from it exceeds the nuclei the cell
      ! has.  xsum_max carries the bound cell by cell; absent, it is one,
      ! which is the hydrogen case of this code (H II is the only carried
      ! stage of an element whose molecular hydrogen the carriers hold
      ! separately).
      !
      ! The projection is the largest scaling toward the interior that
      ! returns the set: negatives are raised to zero first, then, if the
      ! sum still exceeds the bound, every carried fraction is scaled by
      ! xsum_max/sum.  Scaling is used rather than clipping the largest
      ! stage because it keeps the RATIOS of the carried stages, which are
      ! what the charge sum n_e = sum_k Z_k x_k n_el reads.
      !
      ! What it had to move is MEASURED into stage_simplex_sum_over_limit
      ! and stage_fraction_under_zero and reported; a projection that fires
      ! often is a statement about the step, not a repair of it.
      real*8, dimension(:,1-Ng:), intent(inout) :: xion
      real*8, dimension(1-Ng:), intent(in), optional :: xsum_max

      real*8  :: ssum, sc, smax
      integer :: j, k, nk

      nk = size(xion,1)
      stage_simplex_sum_over_limit = -1.0d0
      stage_fraction_under_zero    = -1.0d0

      do j = 1-Ng, N+Ng
         smax = 1.0d0
         if (present(xsum_max)) smax = max(min(xsum_max(j), 1.0d0), 0.0d0)
         ssum = 0.0d0
         do k = 1, nk
            if (xion(k,j) .lt. 0.0d0)                                     &
               stage_fraction_under_zero =                                &
                  max(stage_fraction_under_zero, -xion(k,j))
            xion(k,j) = max(0.0d0, xion(k,j))
            ssum = ssum + xion(k,j)
         enddo
         if (ssum .gt. smax) then
            stage_simplex_sum_over_limit =                                &
               max(stage_simplex_sum_over_limit, ssum - smax)
            sc = smax/ssum
            do k = 1, nk
               xion(k,j) = xion(k,j)*sc
            enddo
         endif
      enddo

      end subroutine stage_simplex_projection

      ! ------------------------------------------------------------------ !

      subroutine hydrogen_and_helium_nucleus_face_flux(rho, Tcode, f_sp,  &
                            Frho_cgs, N_H, N_He, n_Hf, n_Hef, Kf, drf,    &
                            n_H_cell, n_He_cell)
      ! THE TWO ELEMENT NUCLEUS FACE FLUXES AND THE FACE STATE THE STAGE
      ! ROWS ARE WRITTEN AGAINST, all from the element operator's ONE public
      ! flux (element_nucleus_face_flux).  Nothing here rebuilds a face
      ! coefficient of that operator, which is condition C1 of the header.
      !
      ! Frho_cgs(f) is the face mass flux of the state in g cm^-2 s^-1, so
      ! the advective half comes back in the same units (the units paragraph
      ! of element_nucleus_face_flux), and with the diffusive half J(f) in
      ! g cm^-2 s^-1 the two add directly:
      !
      !    N_He(f) = [ F_rho(f) Y_He(f) + J(f) ] / m_He
      !    N_H (f) = [ F_rho(f) - F_rho(f) Y_He(f) - J(f) ] / m_1(f)
      !
      ! the second line being the binary closure of the element operator,
      ! Y_1 = 1 - Y_He and J_1 = -J_He: component 1 carries the hydrogen
      ! nuclei together with the metals and the heavy nuclei bound in the
      ! molecular carriers, and m_1(f) is its mass per hydrogen nucleus, so
      ! the quotient is a hydrogen nucleus flux.  m_1 is one face value
      ! shared by every hydrogen stage, for the reason n_el(f) is.
      !
      ! n_Hf, n_Hef [cm^-3] are the face nucleus densities of the two
      ! elements, Kf [cm^2 s^-1] the face eddy coefficient and drf [cm] the
      ! face spacing, all on the same arithmetic face rule the element
      ! operator's own face coefficients use.
      !
      ! Faces 0 and N carry no diffusive flux, so N_el there is advective
      ! alone; the stage eddy term vanishes at both by
      ! ionization_stage_face_flux.
      !
      ! WHICH HALVES OF THE ELEMENT FLUX THIS RUN ACTUALLY CARRIES.  The
      ! stage rows ride on the element flux of THIS run and not on the flux
      ! a differently configured run would have:
      !
      !   the advective half is present exactly where the caller passes a
      !   nonzero Frho_cgs, which is where the stage row itself carries the
      !   material advection (the fixed-wind relaxation and the stationary
      !   balance); where the Runge-Kutta stages carry it instead the
      !   caller passes zero and the row is the operator-split remainder,
      !   the same convention element_transport_residual and
      !   carrier_residual follow;
      !
      !   the diffusive half is present exactly where the run transports
      !   the elements against each other, he_diffusion with helium in the
      !   mixture.  With that option off no element crosses a face by
      !   diffusion in this run, so a stage riding on one would be carried
      !   by a flux the elements themselves do not have.  The stage's OWN
      !   eddy term, -n_el K dx/dr of equation (1), is not a part of this
      !   and is carried in every configuration: K_zz is a bulk mixing
      !   coefficient of the gas and is blind to which element it mixes.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho_cgs
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: N_H, N_He
      real*8, dimension(0:N),                 intent(out) :: n_Hf, n_Hef
      real*8, dimension(0:N),                 intent(out) :: Kf, drf
      ! The two elements' nucleus densities at the CELL centres [cm^-3],
      ! the quantity the face densities above are the arithmetic mean of.
      ! A stage row's unknown is a fraction of these and its time term is
      ! n_el(j) dx/dt, so they come from the same count and not from a
      ! second one.
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: n_H_cell
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: n_He_cell

      real*8, dimension(0:N) :: Fadv, Jdif, n_el, m_one, n_one
      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe
      real*8  :: m_He_g
      integer :: j

      call element_nucleus_face_flux(rho, Tcode, f_sp, Frho_cgs,          &
                                     Fadv, Jdif, n_el, m_one, n_one)
      if (.not. (he_diffusion .and. thereis_He)) Jdif = 0.0d0

      m_He_g = m_He_amu*mu

      N_H  = 0.0d0
      N_He = 0.0d0
      do j = 0, N
         N_He(j) = (Fadv(j) + Jdif(j))/m_He_g
         N_H(j)  = (Frho_cgs(j) - Fadv(j) - Jdif(j))                      &
                   /max(m_one(j), 1.0d-300)
         n_Hef(j) = n_el(j)
         n_Hf(j)  = n_one(j)
         Kf(j)    = 0.5d0*(kzz_cell(j) + kzz_cell(j+1))
         drf(j)   = (r(j+1) - r(j))*R0
      enddo

      if (present(n_H_cell) .or. present(n_He_cell)) then
         call element_nucleus_counts(f_sp, nucH, nucHe)
         if (present(n_H_cell))  n_H_cell  = nucH *rho*n0
         if (present(n_He_cell)) n_He_cell = nucHe*rho*n0
      endif

      end subroutine hydrogen_and_helium_nucleus_face_flux

      end module ionization_stage_transport
