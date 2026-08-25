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
