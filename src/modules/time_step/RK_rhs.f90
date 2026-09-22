      module RK_integration
      ! Evaluate RK right hand side (convection + source)

      use global_parameters
      use grid_construction, only: spherical_cell_volume
      use Numerical_Fluxes
      use Conversion
      use source_func
      use low_mach_dissipation, only: low_mach_damping_active,          &
                                      contact_mode_dissipation_flux
      ! The well-balanced face data of the reconstruction that produced the
      ! face states this right-hand side is built from ("Well balanced:",
      ! default off).
      use Reconstruction_step, only: wb_dev_L, wb_dev_R, wb_dp_eq,      &
                                     wb_P_up, wb_P_dn

      implicit none

      ! Numerical flux and face pressure at every interface j -- the interface
      ! at r_edg(j), between cells j and j+1 -- as last assembled by RK_rhs.
      ! They are kept so that positivity_limited_fluxes can rebuild the update
      ! of a single cell from a mixture of these fluxes and a first-order
      ! replacement, without re-running the whole right-hand side.
      real*8, dimension(:,:), allocatable :: face_flux
      real*8, dimension(:),   allocatable :: face_p

      ! WELL-BALANCED OPTION.  The face pressure of every interface measured
      ! from the hydrostatic equilibrium of the cell on each side,
      !   face_q_up(f) = face_p(f) - P_up(f),
      !   face_q_dn(f) = face_p(f) - P_dn(f+1),
      ! both of the size of the departure from equilibrium.  They are what
      ! the momentum row is assembled from when that option is on: the O(1)
      ! pressure of the cell and the gravitational source cancel in the
      ! ALGEBRA (see the header of the cell loop) instead of in floating
      ! point.  Filled only when well_balanced is set.
      real*8, dimension(:),   allocatable :: face_q_up, face_q_dn

      ! WELL-BALANCED OPTION.  The magnitude of the pressure force of each
      ! cell's own hydrostatic equilibrium, which is the gravitational
      ! weight the cell carries and is the term the momentum row's residual
      ! is read against under it, since the row itself no longer holds
      ! either of them.  Written out at
      ! equilibrium_pressure_force_of_state, which is the only place it is
      ! formed.  Filled only when well_balanced is set.
      real*8, dimension(:),   allocatable :: equilibrium_pressure_force

      ! THE THREE PHYSICAL TERMS OF EVERY CELL'S MOMENTUM ROW, as the
      ! evaluation that produced dF and S assembled them.  The equation is
      !
      !   d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc
      !                  ^ram          ^pressure ^gravity
      !
      ! and these are its three left-hand terms, gathered out of the pieces
      ! the discretization splits them into:
      !
      !   momentum_ram_divergence(j)     the momentum flux divergence
      !                                  WITHOUT the pressure
      !   momentum_pressure_gradient(j)  the whole spherical pressure
      !                                  gradient the row carries
      !   momentum_gravity(j)            rho dphi/dr as `source` forms it,
      !                                  and the equilibrium pressure force
      !                                  under the well-balanced option, where
      !                                  the row itself carries neither
      !
      ! The three add up to the assembled row dF(2,j) - S(2,j) exactly.
      !
      ! WHAT READS THEM, AND WHY THE DISCRETIZATION'S PIECES ARE NOT THESE
      ! TERMS.  The momentum row's reference scale is the largest PHYSICAL
      ! term of the equation (momentum_row_scale, steady_residual), and
      ! dF(2) and S(2) are not those terms: under PLM the pressure sits
      ! partly in the momentum flux, which Phys_flux gives p, and partly in
      ! the source, which carries the geometric term (A+ - A-) p_c/dV with
      ! the opposite sign, so each of the two holds an O(2 p/r) part that
      ! cancels against the other.  A scale built from them therefore reads
      ! 2 p/r on a uniform pressure at rest, where the physical force is
      ! zero, and under WENO3, where the whole gradient is one number
      ! (p_R - p_L)/dr, it reads |dp/dr| ~ |rho g| in a near-hydrostatic
      ! cell only because the two balance there.  Gathering the pieces back
      ! into the terms of the equation is what these arrays are for.
      ! Nothing here changes dF or S.
      !
      ! THE LOWEST GHOST CELL HAS NO LOWER FACE and carries no equation, so
      ! its entry is zero, as its rows are.
      real*8, dimension(:),   allocatable :: momentum_ram_divergence
      real*8, dimension(:),   allocatable :: momentum_pressure_gradient
      real*8, dimension(:),   allocatable :: momentum_gravity

      ! Number of interfaces whose flux was dropped to first order to keep an
      ! RK stage inside rho > 0, rho e > 0, summed over the whole run, over
      ! the calls that succeeded in restoring the set. A call that returns
      ! repaired = .false. leaves the state to the caller's dt bisection and
      ! adds nothing here, so this counts repaired interfaces and not
      ! attempted ones.
      ! Reported at the end of a run; zero means the high-order fluxes were
      ! admissible everywhere and the run is the one an unguarded build gives.
      integer :: n_faces_flux_positivity_limited = 0

      ! Number of calls to positivity_limited_fluxes, that is the number of
      ! RK stages that left rho > 0, rho e > 0 and were offered the repair,
      ! whether or not the repair succeeded. It is the counter that says
      ! whether a run reached the repair path at all: zero means every stage
      ! stayed inside the admissible set and no first-order flux, and hence
      ! no Lax-Friedrichs signal speed, entered the trajectory.
      integer :: n_calls_flux_positivity_repair = 0

      ! Repaired interfaces belonging to steps the marching loop accepted.
      ! A step retaken at half dt is recomputed from the state at its
      ! beginning, so the repairs of its discarded attempts are not part of
      ! the trajectory the run kept; those are in the total above but not
      ! here. Advanced by the marching loop of EXHALE_main at the point where
      ! it leaves the retry loop with the step accepted.
      integer :: n_faces_flux_positivity_limited_accepted = 0

      contains

      subroutine RK_rhs(u_in,WL,WR,dF,S)
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL,WR
      integer :: j
      real*8 :: dr
      real*8 :: rp,rm
      real*8 :: dAp,dAm
      real*8 :: dV
      real*8, dimension(3) ::  Fp,Fm
      real*8 :: dF3p
      real*8 :: pL,pR
      real*8 :: qp,qm
      ! Gated fourth-difference dissipation of the stagnant-layer contact
      ! mode, added to the numerical flux below so that the marching RHS and
      ! the steady residual (which reaches this routine through
      ! assemble_residual) solve the same equation. Off by default: the array
      ! is then never filled and never read.
      logical :: damp_lowmach
      real*8, dimension(3,1-Ng:N+Ng) :: Ddis
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: dF
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: S

      if (.not. allocated(face_flux)) then
         allocate(face_flux(3,1-Ng:N+Ng), face_p(1-Ng:N+Ng))
         face_flux = 0.0d0
         face_p    = 0.0d0
      endif
      ! The face departures are tested on their own: face_flux and face_p
      ! are also allocated by the kind-generic rows (hydrodynamic_rows),
      ! which do not carry these, and a run that reaches a stationary
      ! evaluation before any marching stage would otherwise find
      ! face_flux allocated and write these unallocated.
      if (.not. allocated(face_q_up)) then
         allocate(face_q_up(1-Ng:N+Ng), face_q_dn(1-Ng:N+Ng))
         face_q_up = 0.0d0
         face_q_dn = 0.0d0
      endif

      ! THE LOWEST GHOST CELL HAS NO LOWER FACE, so no flux difference and
      ! no geometric source are defined for it and its residual row is zero:
      ! the cell loop below starts at 2-Ng, and the column is written here
      ! so that every element of both outputs is defined.  A caller that
      ! evaluates dF - S over the whole padded range (the marching stage
      ! update, whose ghosts Apply_BC then rewrites, and assemble_residual)
      ! must not read an undefined double.
      dF(:,1-Ng) = 0.0d0
      S(:,1-Ng)  = 0.0d0

      if (well_balanced) call equilibrium_pressure_force_of_state(u_in)

      damp_lowmach = low_mach_damping_active()
      if (damp_lowmach) call contact_mode_dissipation_flux(u_in,Ddis)

      ! TWO LOOPS, ONE OVER FACES AND ONE OVER CELLS, AND WHY.  The Riemann
      ! problem at a face is a function of that face's own two reconstructed
      ! states, and the flux difference of a cell is a function of its own
      ! two faces: both are local, and the only thing that tied them into a
      ! single sequential loop was the reuse of the previous cell's
      ! right-hand flux as this cell's left-hand one.  Every flux is stored
      ! in face_flux / face_p in any case -- the positivity repair below
      ! rebuilds a cell from them -- so the face loop writes exactly what the
      ! cell loop reads and the reuse buys nothing.  Split that way both
      ! loops run over independent indices, in ONE parallel region: the
      ! barrier that closes the first !$omp do is what makes every face flux
      ! visible to the cell loop.
      !
      ! BITWISE IDENTICAL TO THE SEQUENTIAL FORM, at any number of threads.
      ! Each face carries the same Num_flux call with the same arguments --
      ! face j lies between cells j and j+1, and the outermost face of the
      ! ghost range has no cell j+1 and so takes the composition of the last
      ! cell that exists -- and each cell the same flux differences of the
      ! same two stored fluxes.  Neither loop contains a reduction, so no
      ! sum changes order.
      !$omp parallel default(shared)                                    &
      !$omp   private(j,dr,rp,rm,dAp,dAm,dV,Fp,Fm,pL,pR,qp,qm,dF3p)

      !$omp do schedule(static)
      do j = 1-Ng,N+Ng
         if (well_balanced) then
            call Num_flux(WL(:,j),WR(:,j),Fp,pR,j,min(j+1,N+Ng),        &
                          wb_dev_L(j),wb_dev_R(j),wb_dp_eq(j),qp,qm)
            face_q_up(j) = qp
            face_q_dn(j) = qm
         else
            call Num_flux(WL(:,j),WR(:,j),Fp,pR,j,min(j+1,N+Ng))
         endif
         if (damp_lowmach) Fp = Fp + Ddis(:,j)
         face_flux(:,j) = Fp
         face_p(j)      = pR
      enddo
      !$omp end do

      !$omp do schedule(static)
      do j = 2-Ng,N+Ng

         ! Substitutions
         dr = dr_j(j)
         rp = r_edg(j)
         rm = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV = spherical_cell_volume(j)

         ! The two interface fluxes of this cell, as the face loop assembled
         ! them
         Fm = face_flux(:,j-1)
         pL = face_p(j-1)
         Fp = face_flux(:,j)
         pR = face_p(j)

         ! Evaluate source
         call source(j,dr,dAp,dAm,dV,    &
                     u_in(:,j),WR(:,j-1),WL(:,j),S(:,j))

         ! Evaluate flux differences
         dF(1,j) = (dAp*Fp(1) - dAm*Fm(1))/dV
         dF(2,j) = (dAp*Fp(2) - dAm*Fm(2))/dV

         ! Correct for WENO3 discretization
         if (use_weno3)  dF(2,j) = dF(2,j) + (pR - pL)/dr

         ! THE MOMENTUM ROW UNDER THE WELL-BALANCED OPTION.  The pressure is
         ! not in Fp(2)/Fm(2) here (Phys_flux leaves it out), and what
         ! replaces it is the face pressure measured from THIS cell's own
         ! hydrostatic equilibrium.  The two identities that make the
         ! substitution exact, with P_up = p_j - rho_j (phi_i(j) - phi_c(j))
         ! and P_dn = p_j + rho_j (phi_c(j) - phi_i(j-1)):
         !
         !   A+ P_up - A- P_dn - (A+ - A-) p_j
         !        = -rho_j [ A+ (phi_i(j) - phi_c(j))
         !                 + A- (phi_c(j) - phi_i(j-1)) ]        (PLM form)
         !   P_up - P_dn = -rho_j (phi_i(j) - phi_i(j-1))        (WENO3 form)
         !
         ! The left-hand sides are the equilibrium part of the pressure
         ! terms the row carries (the flux difference and, under PLM, the
         ! geometric source); the right-hand sides are the discrete
         ! gravitational source of the well-balanced scheme.  They cancel
         ! here in the algebra, cell by cell, with the cell's own O(1)
         ! pressure never appearing, so `source` returns S(2,j) = 0, NEITHER
         ! side is evaluated, and what is left is the departure alone.  This
         ! is the momentum source of Kaeppeli and Mishra (2014, J. Comput.
         ! Phys. 259, 199, their eq. 2.26) in spherical geometry.
         !
         ! What cancels here is still the physics the row balances, and its
         ! magnitude is the row's reference scale: it is formed, from the
         ! right-hand side of whichever identity the assembled row belongs
         ! to, by equilibrium_pressure_force_of_state below.
         if (well_balanced) then
            if (use_plm) then
               dF(2,j) = (dAp*(Fp(2) + face_q_up(j))                    &
                        - dAm*(Fm(2) + face_q_dn(j-1)))/dV
            else
               dF(2,j) = (dAp*Fp(2) - dAm*Fm(2))/dV                     &
                       + (face_q_up(j) - face_q_dn(j-1))/dr
            endif
         endif

         dF3p  = dAp*Fp(1)*(Gphi_i(j) - Gphi_c(j))         &
               - dAm*Fm(1)*(Gphi_i(j-1) - Gphi_c(j))
         dF(3,j) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV

      enddo
      !$omp end do

      !$omp end parallel

      ! The terms of the equation the momentum row of every cell holds, out
      ! of the pieces the loop above assembled it from (the module header
      ! says why the pieces are not the terms).  It reads the stored face
      ! data and S, and writes nothing the rows are built from.
      call momentum_row_terms_of_state(WL,WR,S)

      ! End of subroutine
      end subroutine RK_rhs

      !-----------------------------------------------------------!

      subroutine equilibrium_pressure_force_of_state(u_in)
      ! THE MAGNITUDE OF THE PRESSURE FORCE OF EACH CELL'S OWN HYDROSTATIC
      ! EQUILIBRIUM, and the single definition of it.  With P_up and P_dn
      ! the face values of the constant-density equilibrium through the
      ! cell's own (rho_j, p_j) (Kaeppeli and Mishra 2016, A&A 587, A94,
      ! their eq. 16, in spherical geometry), the identities the momentum
      ! row of the well-balanced option rests on are
      !
      !   A+ P_up - A- P_dn - (A+ - A-) p_j
      !        = -rho_j [ A+ (phi_i(j) - phi_c(j))
      !                 + A- (phi_c(j) - phi_i(j-1)) ]        (PLM form)
      !   P_up - P_dn = -rho_j (phi_i(j) - phi_i(j-1))        (WENO3 form)
      !
      ! and this is the right-hand side of whichever one the assembled row
      ! belongs to, divided by dV or by dr as that row is, and taken as a
      ! magnitude.  It is at once the gravitational weight the cell carries
      ! and the equilibrium part of the pressure terms, which is why the
      ! two cancel in the row and neither appears in it; the physics the
      ! row balances is this size all the same, so it is the term the row's
      ! residual is read against (steady_residual's momentum row scale).
      ! With no gravity it is zero and that scale falls back on the row's
      ! own dynamic terms.
      !
      ! It is a function of the cell's density, of the potential and of the
      ! reconstruction in use, and of nothing the flux assembly does, so it
      ! is the same number whatever arithmetic assembled the rows.
      !
      ! THE LOWEST GHOST CELL HAS NO LOWER FACE and carries no equation, so
      ! its entry is zero, as its rows are.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      integer :: j
      real*8  :: rp, rm

      if (.not. allocated(equilibrium_pressure_force))                   &
         allocate(equilibrium_pressure_force(1-Ng:N+Ng))

      equilibrium_pressure_force(1-Ng) = 0.0d0
      if (use_plm) then
         do j = 2-Ng,N+Ng
            rp = r_edg(j)
            rm = r_edg(j-1)
            equilibrium_pressure_force(j) =                              &
               abs(u_in(1,j)*(rp*rp*(Gphi_i(j)   - Gphi_c(j))            &
                            + rm*rm*(Gphi_c(j)   - Gphi_i(j-1))))        &
               /spherical_cell_volume(j)
         enddo
      else
         do j = 2-Ng,N+Ng
            equilibrium_pressure_force(j) =                              &
               abs(u_in(1,j)*(Gphi_i(j) - Gphi_i(j-1)))/dr_j(j)
         enddo
      endif

      ! End of subroutine
      end subroutine equilibrium_pressure_force_of_state

      !-----------------------------------------------------------!

      subroutine momentum_row_terms_of_cell(dr,dAp,dAm,dV,              &
                                            Fp2,Fm2,pR,pL,qp,qm,        &
                                            S2,weight,w_face_p,         &
                                            ram,pgrad,grav)
      ! ONE CELL'S MOMENTUM ROW SPLIT INTO THE TERMS OF THE PHYSICAL
      ! EQUATION, and the single definition of that split.  The three
      ! outputs add up to the assembled row dF(2) - S(2) exactly, whichever
      ! reconstruction produced it and whether or not the well-balanced
      ! option is on:
      !
      !   PLM    Phys_flux gives the momentum flux the pressure, so the
      !          face value is rho v v + p and the flux difference holds
      !          (A+ p_up - A- p_dn)/dV, while `source` holds the geometric
      !          part (A+ - A-) p_c/dV with the opposite sign.  Both belong
      !          to dp/dr: the pressure piece comes out of the flux
      !          difference and joins the geometric term.
      !   WENO3  Phys_flux leaves the pressure out and RK_rhs adds the face
      !          difference (p_R - p_L)/dr, which is the whole gradient the
      !          row carries; `source` has no geometric term.
      !   WB     Phys_flux leaves the pressure out on either reconstruction
      !          and the row carries the face pressure measured from each
      !          side's own hydrostatic equilibrium, which is the pressure
      !          gradient OF THE DEPARTURE.  The equilibrium part of the
      !          pressure terms and the gravitational source cancel there
      !          analytically, and their common magnitude is the weight the
      !          caller passes (equilibrium_pressure_force).
      !
      ! w_face_p is the weight the caller's assembly gave the face-pressure
      ! difference: one under WENO3, zero under PLM, and recon_lambda where
      ! the positivity repair rebuilds a cell on the PLM to WENO3 homotopy.
      !
      ! S2 AND weight ARE READ ONLY TO RECOVER THE GEOMETRIC TERM, as
      ! p_geom = S2 + weight, which is `source`'s own (A+ - A-) p_c/dV with
      ! the p_c it itself used, and not a second call to the caloric EOS
      ! that would be a second definition of the cell pressure.
      real*8, intent(in)  :: dr,dAp,dAm,dV
      real*8, intent(in)  :: Fp2,Fm2,pR,pL,qp,qm
      real*8, intent(in)  :: S2,weight,w_face_p
      real*8, intent(out) :: ram,pgrad,grav
      real*8 :: p_flux

      grav = weight
      ram  = (dAp*Fp2 - dAm*Fm2)/dV

      if (well_balanced) then
         if (use_plm) then
            pgrad = (dAp*qp - dAm*qm)/dV
         else
            pgrad = (qp - qm)/dr
         endif
         return
      endif

      pgrad = w_face_p*(pR - pL)/dr

      if (use_plm) then
         p_flux = (dAp*pR - dAm*pL)/dV
         ram    = ram   - p_flux
         pgrad  = pgrad + p_flux - (S2 + weight)
      endif

      ! End of subroutine
      end subroutine momentum_row_terms_of_cell

      !-----------------------------------------------------------!

      subroutine momentum_row_terms_of_state(WL,WR,S)
      ! The three terms of every cell's momentum row, from the face data the
      ! flux assembly stored (face_flux, face_p and, under the well-balanced
      ! option, face_q_up / face_q_dn) and the source that assembly returned.
      !
      ! THE GRAVITATIONAL TERM IS THE ONE `source` FORMS, the half-sum of
      ! the two reconstructed face densities times the interface potential
      ! difference over dr (Source.f90); the two must stay in step.  Under
      ! that option `source` returns zero and the weight is the equilibrium
      ! pressure force, which is the same physics in that option's own
      ! discretization (equilibrium_pressure_force_of_state).
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL,WR,S
      integer :: j
      real*8  :: dr,rp,rm,dAp,dAm,dV,wgt,w_face_p

      if (.not. allocated(momentum_ram_divergence))                      &
         allocate(momentum_ram_divergence(1-Ng:N+Ng),                    &
                  momentum_pressure_gradient(1-Ng:N+Ng),                 &
                  momentum_gravity(1-Ng:N+Ng))
      ! The face departures are passed below whether that option is on or not,
      ! and the kind-generic rows allocate face_flux / face_p without them,
      ! so a stationary evaluation through those rows reached before any marching
      ! stage would otherwise pass an unallocated array element.
      if (.not. allocated(face_q_up)) then
         allocate(face_q_up(1-Ng:N+Ng), face_q_dn(1-Ng:N+Ng))
         face_q_up = 0.0d0
         face_q_dn = 0.0d0
      endif

      momentum_ram_divergence(1-Ng)    = 0.0d0
      momentum_pressure_gradient(1-Ng) = 0.0d0
      momentum_gravity(1-Ng)           = 0.0d0

      w_face_p = 0.0d0
      if (use_weno3) w_face_p = 1.0d0

      do j = 2-Ng,N+Ng
         dr  = dr_j(j)
         rp  = r_edg(j)
         rm  = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV  = spherical_cell_volume(j)
         if (well_balanced) then
            wgt = equilibrium_pressure_force(j)
         else
            wgt = 0.5d0*(WR(1,j-1) + WL(1,j))                           &
                  *(Gphi_i(j) - Gphi_i(j-1))/dr
         endif
         call momentum_row_terms_of_cell(dr,dAp,dAm,dV,                 &
                 face_flux(2,j), face_flux(2,j-1),                      &
                 face_p(j), face_p(j-1),                                &
                 face_q_up(j), face_q_dn(j-1),                          &
                 S(2,j), wgt, w_face_p,                                 &
                 momentum_ram_divergence(j),                            &
                 momentum_pressure_gradient(j),                         &
                 momentum_gravity(j))
      enddo

      ! End of subroutine
      end subroutine momentum_row_terms_of_state

      !-----------------------------------------------------------!

      ! Repair an RK stage that left the admissible set rho > 0, rho e > 0 by
      ! replacing, at the offending cells only, the high-order interface
      ! fluxes with the first-order Lax-Friedrichs flux of the neighboring
      ! cell averages, and rebuilding those cells' update from the mixture.
      !
      ! The first-order Lax-Friedrichs update of cell averages is positivity
      ! preserving under dt(|v|+c)/dr <= 1 (Perthame & Shu 1996, Numer. Math.
      ! 73, 119; the LF lemma of Zhang & Shu 2010, J. Comput. Phys. 229,
      ! 3091), and each stage of SSP-RK3 is a convex combination of forward
      ! Euler steps, so it inherits the property. Replacing the flux at a
      ! single interface rather than everywhere is the flux correction of Hu,
      ! Adams & Shu (2013, J. Comput. Phys. 242, 169); the same first-order
      ! interface substitution is standard practice in astrophysical codes
      ! (Stone et al. 2020, ApJS 249, 4, sec. 4.6).
      !
      ! Arguments: stage = 1, 2, 3 of SSP-RK3; u_n = state at the beginning of
      ! the step; u_stage = state this stage was built from (u_n for stage 1);
      ! S = the source term RK_rhs returned for this stage; u_new = the stage
      ! output, repaired in place. repaired is false when the cell averages
      ! themselves are inadmissible or when both interfaces of a cell are
      ! already first order and it is still inadmissible -- no flux choice
      ! repairs that, and the caller falls back on halving dt.
      !
      ! Only cells that are rebuilt are written: every other cell keeps the
      ! bit pattern the uncorrected update gave it.
      !
      ! S is not recomputed. The source depends on the cell state and on the
      ! reconstructed interface states, which positivity_limited_faces has
      ! already made admissible; the correction here changes fluxes only.
      subroutine positivity_limited_fluxes(stage,u_n,u_stage,S,dt_loc,   &
                                           u_new,repaired)

      integer, intent(in) :: stage
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_n,u_stage,S
      real*8, dimension(1-Ng:N+Ng), intent(in) :: dt_loc
      real*8, dimension(3,1-Ng:N+Ng), intent(inout) :: u_new
      logical, intent(out) :: repaired

      real*8, dimension(3,1-Ng:N+Ng) :: W_avg
      real*8, dimension(3,1-Ng:N+Ng) :: flux_lo
      real*8, dimension(1-Ng:N+Ng)   :: p_lo
      real*8, dimension(1-Ng:N+Ng)   :: q_up_lo,q_dn_lo
      logical, dimension(1-Ng:N+Ng)  :: is_first_order
      logical, dimension(1-Ng:N+Ng)  :: rebuild_cell
      integer :: j,jf,sweep,n_new_faces,n_repl
      real*8  :: dr,rp,rm,dAp,dAm,dV,dF3p
      real*8, dimension(3) :: Fp,Fm,dFc
      real*8  :: pL,pR,rho_e
      real*8  :: qp,qm
      real*8  :: wgt,w_face_p

      n_calls_flux_positivity_repair = n_calls_flux_positivity_repair + 1

      ! The momentum row of a rebuilt cell is assembled below, so the term
      ! it is read against belongs here as well.  A replaced face does not
      ! move it: the equilibrium pressure force is a property of the cell's
      ! own density and of the potential, and this stage's state is the one
      ! whose fluxes are being repaired.
      if (well_balanced) call equilibrium_pressure_force_of_state(u_stage)

      repaired       = .false.
      is_first_order = .false.
      flux_lo        = 0.0d0
      p_lo           = 0.0d0
      q_up_lo        = 0.0d0
      q_dn_lo        = 0.0d0
      n_repl         = 0

      ! Cell averages of the state this stage was built from. Their two-cell
      ! pairs are the input of the first-order flux, so the update it produces
      ! is the positivity-preserving one of the references above.
      call U_to_W(u_stage,W_avg)

      ! A sweep marks at least one new interface or finds nothing left to
      ! repair, and there are at most N+1 interfaces bounding cells 1..N, so
      ! the repeat terminates.
      do sweep = 1,N+2

         n_new_faces  = 0
         rebuild_cell = .false.

         ! Locate the cells still outside rho > 0, rho e > 0 and mark their
         ! two interfaces. Each test is the negation of "strictly positive",
         ! so a NaN -- which compares false against everything -- is caught.
         do j = 1,N

            if (u_new(1,j) .gt. 0.0d0) then
               rho_e = u_new(3,j) - 0.5d0*u_new(2,j)*u_new(2,j)/u_new(1,j)
               if (rho_e .gt. 0.0d0) cycle
            endif

            ! Both interfaces already first order: the violation does not come
            ! from the flux discretization (a source term, or a time step past
            ! the positivity bound), and halving dt is the remedy.
            if (is_first_order(j-1) .and. is_first_order(j)) return

            do jf = j-1,j

               if (is_first_order(jf)) cycle

               ! The first-order flux is only positivity preserving if its own
               ! input states are admissible; when they are not, no interface
               ! choice repairs this stage.
               if (.not. (W_avg(1,jf)   .gt. 0.0d0 .and.                &
                          W_avg(3,jf)   .gt. 0.0d0 .and.                &
                          W_avg(1,jf+1) .gt. 0.0d0 .and.                &
                          W_avg(3,jf+1) .gt. 0.0d0)) return

               ! No low-Mach damping term is added on a replaced interface:
               ! the purpose here is to restore positivity, not to put further
               ! numerical dissipation into a cell that already lost it.
               call lax_friedrichs_flux(W_avg(:,jf),W_avg(:,jf+1),      &
                                        flux_lo(:,jf),p_lo(jf),         &
                                        jf,min(jf+1,N+Ng))
               ! Under the well-balanced option the row is assembled from the
               ! face pressure measured against each side's equilibrium, so
               ! a replaced face needs its own pair.  It is formed here by
               ! subtracting two O(1) numbers, which is what that option avoids
               ! everywhere else: a repaired face is not well balanced, as a
               ! first-order Lax-Friedrichs face cannot be (it resolves no
               ! stationary contact).
               if (well_balanced) then
                  q_up_lo(jf) = p_lo(jf) - wb_P_up(jf)
                  q_dn_lo(jf) = p_lo(jf) - wb_P_dn(min(jf+1,N+Ng))
               endif
               is_first_order(jf) = .true.
               n_new_faces = n_new_faces + 1

               if (jf   .ge. 1 .and. jf   .le. N) rebuild_cell(jf)   = .true.
               if (jf+1 .ge. 1 .and. jf+1 .le. N) rebuild_cell(jf+1) = .true.

            enddo

         enddo

         ! Nothing marked means no cell is inadmissible any more.
         if (n_new_faces .eq. 0) then
            repaired = .true.
            exit
         endif
         n_repl = n_repl + n_new_faces

         ! Rebuild the cells that touch a newly replaced interface. The flux
         ! differences below mirror those of RK_rhs above, term by term, with
         ! the replaced interfaces substituted; keep the two in step.
         do j = 1,N

            if (.not. rebuild_cell(j)) cycle

            dr = dr_j(j)
            rp = r_edg(j)
            rm = r_edg(j-1)
            dAp = rp*rp
            dAm = rm*rm
            dV = spherical_cell_volume(j)

            if (is_first_order(j-1)) then
               Fm = flux_lo(:,j-1)
               pL = p_lo(j-1)
               qm = q_dn_lo(j-1)
            else
               Fm = face_flux(:,j-1)
               pL = face_p(j-1)
               qm = face_q_dn(j-1)
            endif

            if (is_first_order(j)) then
               Fp = flux_lo(:,j)
               pR = p_lo(j)
               qp = q_up_lo(j)
            else
               Fp = face_flux(:,j)
               pR = face_p(j)
               qp = face_q_up(j)
            endif

            dFc(1) = (dAp*Fp(1) - dAm*Fm(1))/dV
            dFc(2) = (dAp*Fp(2) - dAm*Fm(2))/dV

            ! Correct for WENO3 discretization. Under the PLM -> WENO3
            ! continuation the stored interface fluxes are already on the
            ! homotopy, so the face-pressure term carries the same weight.
            ! (The first-order replacement flux itself is built with the
            ! flags currently set, i.e. the PLM form; it is reached only on
            ! cells the correction actually repairs, and the counter
            ! n_faces_flux_positivity_limited reports how many those were.)
            if (recon_lambda_on) then
               dFc(2) = dFc(2) + recon_lambda*(pR - pL)/dr
            else if (use_weno3) then
               dFc(2) = dFc(2) + (pR - pL)/dr
            endif

            ! The well-balanced momentum row, as the cell loop of RK_rhs
            ! assembles it; keep the two in step.  Under the PLM to WENO3
            ! continuation the face departures stored here are those of the
            ! LAST of the two evaluations and not their blend (the blend is
            ! formed in reconstruction_continuation_rhs on dF, S, face_flux
            ! and face_p), so a face this repair replaces inside a
            ! continuation ramp carries that evaluation's pressure force rather
            ! than the homotopy's.  A replaced face is not well balanced in
            ! any case (see above), and the repair runs only on a stage that
            ! left rho > 0, rho e > 0.
            if (well_balanced) then
               if (use_plm) then
                  dFc(2) = (dAp*(Fp(2) + qp) - dAm*(Fm(2) + qm))/dV
               else
                  dFc(2) = (dAp*Fp(2) - dAm*Fm(2))/dV + (qp - qm)/dr
               endif
            endif

            ! THE TERMS OF THE MOMENTUM EQUATION THIS ROW HOLDS FOLLOW
            ! THE REBUILD, from the same mixture of fluxes.  The weight is
            ! a property of the cell's own density and of the potential and
            ! no replaced face can move it, so it is the one this stage's
            ! own right-hand side left.
            w_face_p = 0.0d0
            if (recon_lambda_on) then
               w_face_p = recon_lambda
            else if (use_weno3) then
               w_face_p = 1.0d0
            endif
            wgt = momentum_gravity(j)
            call momentum_row_terms_of_cell(dr,dAp,dAm,dV,              &
                    Fp(2),Fm(2),pR,pL,qp,qm,S(2,j),wgt,w_face_p,        &
                    momentum_ram_divergence(j),                         &
                    momentum_pressure_gradient(j),                      &
                    momentum_gravity(j))

            dF3p  = dAp*Fp(1)*(Gphi_i(j) - Gphi_c(j))         &
                  - dAm*Fm(1)*(Gphi_i(j-1) - Gphi_c(j))
            dFc(3) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV

            ! Same stage updates as the marching loop in EXHALE_main
            select case (stage)
            case (1)
               u_new(:,j) = u_n(:,j) - dt_loc(j)*(dFc - S(:,j))
            case (2)
               u_new(:,j) = (3.0*u_n(:,j) + u_stage(:,j)                &
                            - dt_loc(j)*(dFc - S(:,j)))/4.0
            case (3)
               u_new(:,j) = (u_n(:,j) + 2.0*(u_stage(:,j)               &
                            - dt_loc(j)*(dFc - S(:,j))))/3.0
            end select

         enddo

      enddo

      ! THE STORED FACE FLUXES FOLLOW THE REPAIR.  The species rows ride on
      ! face_flux(1,:), so a face whose flux was replaced here has to carry
      ! the replacement there as well: the mass row of every rebuilt cell is
      ! the divergence of exactly these fluxes, and a species flux built on
      ! the discarded high-order flux would not sum to it.  A cell that was
      ! not rebuilt cannot touch a replaced face, because replacing a face
      ! marks both of its cells, so writing them back leaves every cell's
      ! mass update the divergence of the array.
      do jf = 1-Ng,N+Ng
         if (is_first_order(jf)) then
            face_flux(:,jf) = flux_lo(:,jf)
            face_p(jf)      = p_lo(jf)
            if (well_balanced) then
               face_q_up(jf) = q_up_lo(jf)
               face_q_dn(jf) = q_dn_lo(jf)
            endif
         endif
      enddo

      if (repaired)                                                     &
         n_faces_flux_positivity_limited =                              &
            n_faces_flux_positivity_limited + n_repl

      ! End of subroutine
      end subroutine positivity_limited_fluxes

      ! End of module
      end module RK_integration
