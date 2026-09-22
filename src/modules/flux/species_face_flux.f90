      module species_advective_transport
      ! Advective transport of species mass fractions on the SAME faces, with
      ! the SAME mass flux, as the density.
      !
      ! THE IDENTITY THIS MODULE EXISTS FOR.  There is one set of faces,
      ! r_edg(1-Ng:N+Ng), and one mass flux per face, F_rho(j), the one the
      ! Riemann solver returns and RK_rhs stores.  Every species mass flux is
      ! that scalar multiplied by a face composition that sums to one,
      !
      !     F_s(j) = F_rho(j) Y_s^face(j),     sum_s Y_s^face(j) = 1,
      !
      ! so that
      !
      !     sum_s F_s(j) = F_rho(j)
      !
      ! holds by construction and not by tolerance.  Y_s is the species mass
      ! fraction with the electron mass carried on its ion, so the free
      ! electrons contribute no mass flux.  The closing member of the set is
      ! not reconstructed: it is one minus the reconstructed others, which is
      ! what makes the sum exact.
      !
      ! WHY THIS FORM AND NOT A CELL-VELOCITY UPWIND DIFFERENCE.  A cell can
      ! then only lose what the mass row says it loses.  The two failures the
      ! cell-centered forms each worked around separately do not arise:
      !   * a cell whose two FACE velocities straddle zero has an advective
      !     term of its own here, because the term is a difference of two face
      !     fluxes and not a one-sided difference of cell values;
      !   * a cell that is outflowing at both faces cannot be evacuated of a
      !     trace element, because its species outflow is its mass outflow
      !     multiplied by a fraction and the density is drained by exactly the
      !     same fluxes.
      !
      ! THE UPWIND SIDE.  The face composition is taken from the side the face
      ! MASS FLUX selects: F_rho(j) >= 0 takes the left state, F_rho(j) < 0
      ! the right one.  With HLLC the same side is the one the contact wave
      ! selects: in the star region the material speed is S_star and the mass
      ! flux across the face is rho_star S_star with rho_star > 0, so the two
      ! signs agree; where the fan is supersonic the face flux is the upwind
      ! state's own rho v and the signs agree again.  ROE and the
      ! Lax-Friedrichs flux return one flux without a contact state, so the
      ! sign of F_rho is the only rule defined for them, and it is the rule
      ! used for all three.  A face with F_rho = 0 carries no species flux and
      ! its side is immaterial.
      !
      ! NO LOW-MACH DISSIPATION TERM (design decision 4).  The dissipation the
      ! momentum and energy rows carry has an identically zero mass component,
      ! so the species rows, which ride on the mass flux, carry none either.
      !
      ! BOUNDS.  A reconstructed face fraction is scaled back toward its own
      ! donor cell average by the largest factor that keeps it in [0,1].  That
      ! is a property of a mass fraction, not of a reconstruction, which is
      ! why it is imposed here and not in Reconstruct_scalar.  How often it
      ! fires and how far the reconstruction had left the range are counted
      ! and exposed rather than assumed.

      use global_parameters
      use Reconstruction_step, only: Reconstruct_scalar
      use grid_construction, only: spherical_cell_volume

      implicit none
      private

      public :: species_face_fraction
      public :: species_face_flux
      public :: species_flux_divergence
      public :: species_advective_update
      ! Number of FACE STATES whose reconstructed mass fraction left [0,1] and
      ! was scaled back toward its own donor cell average, summed over the
      ! run, and the largest excursion seen before the scaling.  Zero means
      ! the reconstruction stayed on the composition simplex everywhere.
      integer, public :: n_species_faces_bounded = 0
      real*8,  public :: species_face_excursion  = 0.0d0
      ! How far the UPDATED cell mass fraction left [0,1], signed so that a
      ! negative value means it stayed inside.  It is measured and not
      ! clipped: a clipped value cannot tell a bound that holds from one
      ! that is being enforced, and the bound belongs to the scheme.
      real*8,  public :: species_fraction_excursion = -1.0d0
      ! Largest relative face-identity residual seen, kept so that a caller or
      ! a test reads the number the operator itself measured.
      real*8,  public :: species_face_identity_residual = 0.0d0
      public :: species_face_identity_injection

      ! Face at which a deliberate violation of the face identity is injected,
      ! so that the assertion below can be shown to fire.  Unset (a value
      ! below 1-Ng) means no injection, which is every production run.
      integer :: inject_face = -huge(1)

      ! The identity is exact in exact arithmetic; what is left is the
      ! round-off of one multiplication and one sum per member of the set.
      real*8, parameter :: identity_tol = 64.0d0*epsilon(1.0d0)

      contains

      !-----------------------------------------------!

      subroutine species_face_identity_injection(jface)
      ! Arm (jface >= 1-Ng) or disarm (any smaller value) the injected
      ! violation of the face identity.  It exists because an assertion that
      ! no run has ever made fire is not known to fire.
      integer, intent(in) :: jface
      inject_face = jface
      end subroutine species_face_identity_injection

      !-----------------------------------------------!

      subroutine species_face_fraction(Y,Frho,Yf)
      ! Face composition of one species: the reconstruction of its cell mass
      ! fraction, evaluated on the side the face mass flux selects and scaled
      ! back onto [0,1] against its own donor cell average.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Y,Frho
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Yf

      real*8, dimension(1-Ng:N+Ng) :: YL,YR
      real*8  :: Yav,th,exc
      integer :: j,jd

      call Reconstruct_scalar(Y,YL,YR)

      do j = 1-Ng,N+Ng
         if (Frho(j) .ge. 0.0d0) then
            Yf(j) = YL(j)
            jd    = j
         else
            Yf(j) = YR(j)
            jd    = min(j+1,N+Ng)
         endif
         ! The anchor of the scaling below: the donor cell's own average,
         ! held inside [0,1].  A cell average is a mass fraction and belongs
         ! there, but the ionization solve can hand back a trace species at a
         ! small negative density, and an anchor outside [0,1] is not an
         ! admissible face value to scale toward: with it the denominators
         ! below can vanish (Yav = Yf, reached by a negative average and a
         ! negative face value; MEASURED as a SIGFPE under -ffpe-trap=zero on
         ! the metals row).  Clamped, Yf < 0 <= Yav and Yav <= 1 < Yf hold in
         ! the two branches, so each denominator is strictly positive.
         Yav = min(1.0d0, max(0.0d0, Y(jd)))
         if (Yf(j) .lt. 0.0d0 .or. Yf(j) .gt. 1.0d0 .or.                  &
             Yf(j) .ne. Yf(j)) then
            exc = max(-Yf(j), Yf(j) - 1.0d0)
            if (Yf(j) .ne. Yf(j)) exc = huge(1.0d0)
            species_face_excursion = max(species_face_excursion, exc)
            n_species_faces_bounded = n_species_faces_bounded + 1
            ! Largest scaling toward the donor cell average that puts the
            ! face value back on [0,1].  The anchor is in [0,1] by the clamp
            ! above, so the scaling exists and its denominator cannot vanish.
            th = 1.0d0
            if (Yf(j) .ne. Yf(j)) then
               th = 0.0d0
            else if (Yf(j) .lt. 0.0d0) then
               th = Yav/(Yav - Yf(j))
            else
               th = (1.0d0 - Yav)/(Yf(j) - Yav)
            endif
            th = max(0.0d0, min(1.0d0, th))
            Yf(j) = Yav + th*(Yf(j) - Yav)
            Yf(j) = max(0.0d0, min(1.0d0, Yf(j)))
         endif
      enddo

      end subroutine species_face_fraction

      !-----------------------------------------------!

      subroutine species_face_flux(Frho,Yf,Fs)
      ! The species mass flux carried by a face: the face mass flux times the
      ! face composition.  One multiplication, so that the sum over a set of
      ! fractions summing to one returns the mass flux itself.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Frho,Yf
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Fs
      Fs = Frho*Yf
      end subroutine species_face_flux

      !-----------------------------------------------!

      subroutine species_flux_divergence(Fs,divF,divmag)
      ! Divergence of one species face flux, per unit volume, on the faces,
      ! the areas and the volume of the MASS ROW of the same cell:
      !
      !     divF(j) = ( A_+ F_s(j) - A_- F_s(j-1) ) / dV_j ,
      !     A_+ = r_{j+1/2}^2,  A_- = r_{j-1/2}^2,
      !     dV_j = ( r_{j+1/2}^3 - r_{j-1/2}^3 ) / 3 .
      !
      ! THERE IS ONE OF THESE DIVERGENCES IN THE CODE, and every equation
      ! that transports a species on the mass flux reads it here: the
      ! Runge-Kutta stages below, the stationary carrier rows and the
      ! stationary elemental rows.  Two spellings of one term give a marched
      ! state and a Newton state that are different objects, which is what
      ! the stationary rows carried until the cell-velocity upwind
      ! difference was replaced by this call.
      !
      ! The volume itself is grid_construction's spherical_cell_volume, the
      ! one expression for it: the hydrodynamic rows that drain the same cell
      ! of mass divide by that value, and a species row that divided by a
      ! second spelling of it would take a cell's composition away from the
      ! mass it rides on at the size of the difference between the two.
      !
      ! divmag is the same expression with both fluxes taken in magnitude,
      ! A_+|F_s(j)| + A_-|F_s(j-1)| over dV_j: the size of the term, for a
      ! row scale that has to say how large the terms it balances are.
      !
      ! Cell 1-Ng has no left face inside the array and is left at zero.
      real*8, dimension(1-Ng:N+Ng), intent(in)            :: Fs
      real*8, dimension(1-Ng:N+Ng), intent(out)           :: divF
      real*8, dimension(1-Ng:N+Ng), intent(out), optional :: divmag

      real*8  :: rp,rm,dAp,dAm,dV
      integer :: j

      divF = 0.0d0
      if (present(divmag)) divmag = 0.0d0
      do j = 2-Ng,N+Ng
         rp  = r_edg(j)
         rm  = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV  = spherical_cell_volume(j)
         divF(j) = (dAp*Fs(j) - dAm*Fs(j-1))/dV
         if (present(divmag))                                             &
            divmag(j) = (dAp*abs(Fs(j)) + dAm*abs(Fs(j-1)))/dV
      enddo

      end subroutine species_flux_divergence

      !-----------------------------------------------!

      subroutine species_advective_update(stage,nq,n_norm,Y_n,Y_in,Frho,  &
                                          rho_n,rho_in,rho_new,dt_loc,    &
                                          jlo,Y_out)
      ! One SSP-RK3 stage of the species rows, on the faces, the volumes and
      ! the time step of the mass row of the same stage.
      !
      !     (rho Y_s)^stage = combination of (rho Y_s)^n and (rho Y_s)^in
      !                       - dt ( A_+ F_s(j) - A_- F_s(j-1) ) / dV
      !     Y_s^stage       = (rho Y_s)^stage / rho^stage
      !
      ! with rho^stage the density the mass row of this stage returned, never
      ! a density recomputed from the composition.  A uniform composition is
      ! then reproduced to round-off whatever the velocity field does,
      ! because the species update is the mass update multiplied by the same
      ! constant.
      !
      ! nq is the number of species carried; the first n_norm of them form
      ! the set whose face fractions are normalized, and the closing member
      ! of that set (index 0, not stored) takes one minus their sum.  Members
      ! beyond n_norm ride on the same faces without entering the
      ! normalization, which is the discretization of a species whose mass is
      ! already accounted for inside the closing member: the trace limit, and
      ! it is valid exactly where that trace closure is.
      !
      ! Rows 1-Ng .. jlo-1 and N+1 .. N+Ng are not updated: they are the
      ! boundary composition, which the caller imposes.
      integer, intent(in) :: stage,nq,n_norm,jlo
      real*8, dimension(1-Ng:N+Ng,nq), intent(in)  :: Y_n,Y_in
      real*8, dimension(1-Ng:N+Ng),    intent(in)  :: Frho
      real*8, dimension(1-Ng:N+Ng),    intent(in)  :: rho_n,rho_in,rho_new
      real*8, dimension(1-Ng:N+Ng),    intent(in)  :: dt_loc
      real*8, dimension(1-Ng:N+Ng,nq), intent(out) :: Y_out

      real*8, dimension(1-Ng:N+Ng,nq) :: Yf,Fs,divFs
      real*8, dimension(1-Ng:N+Ng)    :: Yclose,Fclose
      real*8  :: rhoY,res,scal
      integer :: j,q

      do q = 1,nq
         call species_face_fraction(Y_in(:,q),Frho,Yf(:,q))
      enddo

      ! The closing member of the normalized set, and with it the exactness
      ! of sum_s F_s = F_rho.
      Yclose = 1.0d0
      do q = 1,n_norm
         Yclose = Yclose - Yf(:,q)
      enddo

      do q = 1,nq
         call species_face_flux(Frho,Yf(:,q),Fs(:,q))
      enddo
      call species_face_flux(Frho,Yclose,Fclose)

      ! The injected violation is put into the FLUX and not into the
      ! fraction, so that the fractions still sum to one and only the
      ! identity below is broken.
      if (inject_face .ge. 1-Ng .and. inject_face .le. N+Ng .and.         &
          nq .ge. 1)                                                      &
         Fs(inject_face,1) = Fs(inject_face,1)                            &
                           + 1.0d-6*max(abs(Frho(inject_face)),1.0d-30)

      ! THE FACE IDENTITY, ASSERTED HERE AND NOT IN A TEST (design row B4-j).
      ! It is a statement about every face of every stage of every step, so
      ! it is checked where it is built.  The faces checked are the ones the
      ! update below uses, jlo-1 to N: the outer ghost faces bound no updated
      ! cell and carry no species row to be consistent with.
      do j = jlo-1,N
         res = Fclose(j)
         do q = 1,n_norm
            res = res + Fs(j,q)
         enddo
         res = abs(res - Frho(j))
         scal = abs(Frho(j))
         if (scal .gt. 0.0d0) then
            res = res/scal
         else if (res .gt. 0.0d0) then
            res = huge(1.0d0)
         endif
         species_face_identity_residual =                                 &
              max(species_face_identity_residual,res)
         if (.not. (res .le. identity_tol)) then
            write(*,'(A,I0,A,ES12.5,A,ES12.5)')                           &
               ' ERROR: the species face fluxes do not sum to the mass'// &
               ' flux at face ', j, ': relative residual ', res,          &
               ' against a tolerance of ', identity_tol
            error stop 1
         endif
      enddo

      Y_out = Y_in

      ! The one flux divergence of this operator, shared with the stationary
      ! rows that balance the same term (species_flux_divergence).
      do q = 1,nq
         call species_flux_divergence(Fs(:,q),divFs(:,q))
      enddo

      do j = jlo,N

         do q = 1,nq

            select case (stage)
            case (1)
               rhoY = rho_n(j)*Y_n(j,q) - dt_loc(j)*divFs(j,q)
            case (2)
               rhoY = (3.0*rho_n(j)*Y_n(j,q) + rho_in(j)*Y_in(j,q)        &
                      - dt_loc(j)*divFs(j,q))/4.0
            case (3)
               rhoY = (rho_n(j)*Y_n(j,q)                                  &
                      + 2.0*(rho_in(j)*Y_in(j,q)                          &
                             - dt_loc(j)*divFs(j,q)))/3.0
            case default
               write(*,*) 'ERROR: species_advective_update: stage ',stage
               error stop 1
            end select

            if (rho_new(j) .gt. 0.0d0) then
               Y_out(j,q) = rhoY/rho_new(j)
            else
               Y_out(j,q) = Y_in(j,q)
            endif

            species_fraction_excursion = max(species_fraction_excursion,  &
                 -Y_out(j,q), Y_out(j,q) - 1.0d0)

         enddo

      enddo

      end subroutine species_advective_update

      ! End of module
      end module species_advective_transport
