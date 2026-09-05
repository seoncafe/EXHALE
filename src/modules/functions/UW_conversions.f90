      module Conversion
      ! Conversion conservative-primitive variables subroutines
      
      use global_parameters
      use caloric_eos, only: caloric_mixture_active,                    &
                             pressure_from_energy_density,              &
                             energy_density_from_pressure
      
      implicit none
      
      contains
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive (one component).
      ! jcell names the cell whose COMPOSITION this state has: for a cell
      ! average it is the cell itself, for a reconstructed face state the
      ! cell the state was reconstructed from.  The caloric EOS needs it
      ! because the heat capacity is a property of the mixture, not of the
      ! three numbers in U.
      subroutine U_to_W_comp(U_in,W_out,jcell)
      real*8, intent(in) :: U_in(3)
      integer, intent(in) :: jcell
      real*8, intent(out) :: W_out(3)
      
      W_out(1) = U_in(1)
      W_out(2) = U_in(2)/U_in(1)
      W_out(3) = pressure_from_energy_density(jcell, U_in(1),            &
                 U_in(3)-0.5*U_in(2)*U_in(2)/U_in(1))
 
      end subroutine U_to_W_comp
      
      !-----------------------------------------------------------!
      
      ! Primitive to conservative (one component); jcell as above.
      subroutine W_to_U_comp(W_in,U_out,jcell)
      real*8, intent(in) :: W_in(3)
      integer, intent(in) :: jcell
      real*8, intent(out) :: U_out(3)
      
      U_out(1) = W_in(1)
      U_out(2) = W_in(1)*W_in(2)
      U_out(3) = 0.5*W_in(1)*W_in(2)**2.0                                &
                 + energy_density_from_pressure(jcell, W_in(1), W_in(3))
 
      end subroutine W_to_U_comp
      
      !-----------------------------------------------------------!
      
      ! Primitive to conservative
      subroutine W_to_U(W_in,U_out)
      real*8, intent(in) :: W_in(3,1-Ng:N+Ng)
      real*8, intent(out) :: U_out(3,1-Ng:N+Ng)
      integer :: j
      
      U_out(1,:) = W_in(1,:)
      U_out(2,:) = W_in(1,:)*W_in(2,:)
      ! The array form is kept verbatim for a gas with no molecules in it,
      ! so that an atomic run reproduces the constant-gamma arithmetic to
      ! the bit; the cell loop is the same statement with the cell's own
      ! heat capacity in place of the constant.
      ! Cell-local: the caloric EOS of a cell is a function of that cell's
      ! own state and its own composition, so the loop carries no
      ! dependence between cells and no reduction.  Each cell's expression
      ! is untouched, so the result is bitwise the sequential one at any
      ! number of threads.  No caller is itself inside a parallel region.
      if (caloric_mixture_active) then
         !$omp parallel do default(shared) private(j) schedule(static)
         do j = 1-Ng, N+Ng
            U_out(3,j) = 0.5*W_in(1,j)*W_in(2,j)**2.0                    &
               + energy_density_from_pressure(j, W_in(1,j), W_in(3,j))
         enddo
         !$omp end parallel do
      else
      U_out(3,:) = 0.5*W_in(1,:)*W_in(2,:)**2.0    &
                   + W_in(3,:)/(gamma_ad-1.0)
      endif
 
      end subroutine W_to_U
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive
      subroutine U_to_W(U_in,W_out)
      real*8, intent(in) :: U_in(3,1-Ng:N+Ng)
      real*8, intent(out) :: W_out(3,1-Ng:N+Ng)
      integer :: j
      
      W_out(1,:) = U_in(1,:)
      W_out(2,:) = U_in(2,:)/U_in(1,:)
      ! Cell-local, and the expensive half of the pair: the pressure of a
      ! molecular cell is the inverse of the caloric EOS and costs a Newton
      ! iteration per cell (temperature_of_mixture), which is why this loop
      ! rather than the array statement below is what a molecular run
      ! spends its reconstruction time in.  Same argument as W_to_U above:
      ! one cell, one composition, no reduction, bitwise unchanged.
      if (caloric_mixture_active) then
         !$omp parallel do default(shared) private(j) schedule(static)
         do j = 1-Ng, N+Ng
            W_out(3,j) = pressure_from_energy_density(j, U_in(1,j),      &
               U_in(3,j)-0.5*U_in(2,j)*U_in(2,j)/U_in(1,j))
         enddo
         !$omp end parallel do
      else
      W_out(3,:) = (gamma_ad-1.0)*                   &
                   (U_in(3,:)-0.5*U_in(2,:)*U_in(2,:)/U_in(1,:))
      endif
 
      end subroutine U_to_W

      !-----------------------------------------------------------!

      ! Is the conservative state thermodynamically admissible, i.e. does
      ! every physical cell carry a positive mass density AND a positive
      ! internal energy rho e = E - rho v^2 / 2?
      !
      ! That set is where the Euler equations are defined and where the sound
      ! speed sqrt(gamma p/rho) exists. A conservative update preserves it only
      ! under a CFL bound -- dt (|v| + c)/dr <= 1/2 for the first-order HLL
      ! family (Einfeldt et al. 1991; Batten et al. 1997), and tighter with a
      ! high-order reconstruction -- which the CFL number of a run does not
      ! enforce. The marching loop tests each RK stage with this function and
      ! retakes the step at half dt when it fails.
      !
      ! Only the physical cells 1..N are tested: they are what the update
      ! writes, and the ghosts are set from them by Apply_BC.
      !
      ! Each test is the negation of "strictly positive", so a NaN -- which
      ! compares false against everything -- is inadmissible too.
      logical function positive_density_and_internal_energy(U_in)
      real*8, intent(in) :: U_in(3,1-Ng:N+Ng)
      integer :: j
      real*8  :: rho_e

      positive_density_and_internal_energy = .false.
      do j = 1,N
         if (.not. (U_in(1,j) .gt. 0.0d0)) return
         rho_e = U_in(3,j) - 0.5d0*U_in(2,j)*U_in(2,j)/U_in(1,j)
         if (.not. (rho_e .gt. 0.0d0)) return
      enddo
      positive_density_and_internal_energy = .true.

      end function positive_density_and_internal_energy

      !-----------------------------------------------------------!

      ! Conservative to primitive on the interior cells 1..N only; the
      ! ghost columns of W_out are set to zero. For a state whose ghosts
      ! are not filled yet -- before Apply_BC, or a line-search trial that
      ! unpacks only the interior -- U_to_W over the whole array would
      ! evaluate 0/0 in the ghosts (harmless when the result is discarded,
      ! fatal under -ffpe-trap=invalid).
      subroutine U_to_W_interior(U_in,W_out)
      real*8, intent(in) :: U_in(3,1-Ng:N+Ng)
      real*8, intent(out) :: W_out(3,1-Ng:N+Ng)
      integer :: j
      
      W_out = 0.0d0
      do j = 1,N
         call U_to_W_comp(U_in(:,j),W_out(:,j),j)
      enddo
      
      end subroutine U_to_W_interior
        
      !-----------------------------------------------------------!
      
      ! End of module
      end module Conversion
