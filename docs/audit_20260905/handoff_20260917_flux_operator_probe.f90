program flux_operator_probe
  ! Review diagnostic: compare reconstruction operators on one restart.
  ! Link current production objects; run only in an isolated input directory.
  use global_parameters
  use Read_input, only: input_read
  use Initialization, only: init
  use lya_rt, only: lya_rt_allocate_arrays
  use excited_hydrogen, only: excited_H_allocate_arrays
  use ionization_equilibrium, only: ioniz_eq_allocate_arrays
  use mol_rates, only: h2_thermochemistry_init
  use Conversion, only: W_to_U
  use BC_Apply, only: Apply_BC
  use Reconstruction_step, only: Reconstruct
  use RK_integration, only: RK_rhs, face_flux
  use base_boundary, only: wind_window_mass_flux
  implicit none
  real*8, allocatable :: W(:,:), u(:,:), f(:,:), WL(:,:), WR(:,:), dF(:,:), S(:,:)
  real*8 :: F0, spread, flux_min, flux_max
  logical :: have
  integer :: mode, j, refresh
  call input_read
  open(unit=outfile, file='EXHALE_setup.out')
  call lya_rt_allocate_arrays
  call excited_H_allocate_arrays
  call ioniz_eq_allocate_arrays
  call h2_thermochemistry_init
  allocate(W(3,1-Ng:N+Ng), u(3,1-Ng:N+Ng), f(1-Ng:N+Ng,n_species))
  allocate(WL(3,1-Ng:N+Ng), WR(3,1-Ng:N+Ng), dF(3,1-Ng:N+Ng), S(3,1-Ng:N+Ng))
  call init(W,u,f)
  call wind_window_mass_flux(W,F0,spread,have)
  write(*,'(A,A)') 'INPUT reconstruction=',trim(rec_method)
  write(*,'(A,ES24.16,A,L1)') 'WINDOW F0=',F0,' available=',have
  do refresh=0,1
    do mode=1,2
      if (mode == 1) then
        rec_method='PLM'; use_plm=.true.; use_weno3=.false.
      else
        rec_method='WENO3'; use_plm=.false.; use_weno3=.true.
      endif
      call W_to_U(W,u)
      if (refresh == 1) call Apply_BC(u)
      call Reconstruct(u,WL,WR)
      call RK_rhs(u,WL,WR,dF,S)
      write(*,'(A,A,A,I1)') 'OPERATOR ',trim(rec_method),' refresh_BC=',refresh
      do j=0,6
        write(*,'(A,I3,A,ES24.16)') 'FACE ',j,' flux/F0=',face_flux(1,j)*r_edg(j)**2/F0
      enddo
      flux_min=minval(face_flux(1,0:N)*r_edg(0:N)**2)/F0
      flux_max=maxval(face_flux(1,0:N)*r_edg(0:N)**2)/F0
      write(*,'(A,2ES24.16)') 'ALL_FACES min,max /F0=',flux_min,flux_max
    enddo
  enddo
  close(outfile)
end program flux_operator_probe
