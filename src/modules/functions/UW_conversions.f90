      module Conversion
      ! Conversion conservative-primitive variables subroutines
      
      use global_parameters
      
      implicit none
      
      contains
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive (one component)
      subroutine U_to_W_comp(U_in,W_out)
      real*8, intent(in) :: U_in(3)
      real*8, intent(out) :: W_out(3)
      
      W_out(1) = U_in(1)
      W_out(2) = U_in(2)/U_in(1)
      W_out(3) = (g-1.0)*(U_in(3)-0.5*U_in(2)*U_in(2)/U_in(1))
 
      end subroutine U_to_W_comp
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive (one component)
      subroutine W_to_U_comp(W_in,U_out)
      real*8, intent(in) :: W_in(3)
      real*8, intent(out) :: U_out(3)
      
      U_out(1) = W_in(1)
      U_out(2) = W_in(1)*W_in(2)
      U_out(3) = 0.5*W_in(1)*W_in(2)**2.0 + W_in(3)/(g-1.0)
 
      end subroutine W_to_U_comp
      
      !-----------------------------------------------------------!
      
      ! Primitive to conservative
      subroutine W_to_U(W_in,U_out)
      real*8, intent(in) :: W_in(3,1-Ng:N+Ng)
      real*8, intent(out) :: U_out(3,1-Ng:N+Ng)
      
      U_out(1,:) = W_in(1,:)
      U_out(2,:) = W_in(1,:)*W_in(2,:)
      U_out(3,:) = 0.5*W_in(1,:)*W_in(2,:)**2.0    &
                   + W_in(3,:)/(g-1.0)
 
      end subroutine W_to_U
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive
      subroutine U_to_W(U_in,W_out)
      real*8, intent(in) :: U_in(3,1-Ng:N+Ng)
      real*8, intent(out) :: W_out(3,1-Ng:N+Ng)
      
      W_out(1,:) = U_in(1,:)
      W_out(2,:) = U_in(2,:)/U_in(1,:)
      W_out(3,:) = (g-1.0)*                   &
                   (U_in(3,:)-0.5*U_in(2,:)*U_in(2,:)/U_in(1,:))
 
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
         call U_to_W_comp(U_in(:,j),W_out(:,j))
      enddo
      
      end subroutine U_to_W_interior
        
      !-----------------------------------------------------------!
      
      ! End of module
      end module Conversion
