      subroutine write_frozen_derivative_steps(y_ref,direction,j_ref,f_ref, &
                                                npart,fcomp,h_ref,c_ref)
      ! A diagnostic of the frozen-source residual, not the live chemical
      ! Jacobian. Both sides stay inside the positive gas domain.
      real*8, intent(in) :: y_ref(:),direction(:),j_ref(:),f_ref(:)
      real*8, intent(in) :: npart(1-Ng:N+Ng), h_ref(1-Ng:N+Ng)
      real*8, intent(in) :: c_ref(1-Ng:N+Ng),fcomp(1-Ng:N+Ng,n_species)
      real*8 :: yp(size(y_ref)),ym(size(y_ref)),fp(size(y_ref)),fm(size(y_ref))
      real*8 :: dp(size(y_ref)),dm(size(y_ref)),dc(size(y_ref))
      real*8 :: hp,hm,h_requested,scale,errp,errm,errc,asymmetry
      integer :: exponent,unit_d,loc(1),cell_d,slot_d
      open(newunit=unit_d,file='output/frozen_derivative_steps.txt', &
           status='replace',action='write')
      write(unit_d,'(A)') '# requested_h h_plus h_minus error_forward'// &
           ' error_backward error_central asymmetry worst_cell worst_slot'
      scale=max(maxval(abs(j_ref)),1.0d-300)
      do exponent=4,10
         h_requested=10.0d0**(-exponent)
         hp=h_requested; hm=h_requested
         call gas_admissible_perturbation(y_ref,direction,hp,yp,       &
                                          'frozen forward difference')
         call gas_admissible_perturbation(y_ref,-direction,hm,ym,     &
                                          'frozen backward difference')
         call frozen_residual(yp,npart,fcomp,h_ref,c_ref,fp)
         call frozen_residual(ym,npart,fcomp,h_ref,c_ref,fm)
         dp=(fp-f_ref)/hp; dm=(f_ref-fm)/hm
         dc=(hm*dp+hp*dm)/(hp+hm)
         errp=maxval(abs(dp-j_ref))/scale
         errm=maxval(abs(dm-j_ref))/scale
         errc=maxval(abs(dc-j_ref))/scale
         asymmetry=maxval(abs(dp-dm))/scale
         loc=maxloc(abs(dc-j_ref))
         cell_d=(loc(1)-1)/nvar_jac+1
         slot_d=mod(loc(1)-1,nvar_jac)+1
         write(unit_d,'(7ES24.15E3,2I8)') h_requested,hp,hm,errp,errm, &
                                           errc,asymmetry,cell_d,slot_d
         write(*,'(A,ES10.3,A,ES10.3,A,I0,A,I0)') ' (jac_sweep) h=', &
              h_requested,' central error=',errc,' cell=',cell_d,   &
              ' slot=',slot_d
      enddo
      close(unit_d)
      end subroutine write_frozen_derivative_steps
