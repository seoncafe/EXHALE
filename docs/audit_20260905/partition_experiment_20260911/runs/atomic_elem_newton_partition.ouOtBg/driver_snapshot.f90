! Diagnostic driver for the current production modules. No production source
! is replaced. Each invocation initializes one copied fixture without marching.
program transport_wind_experiment
  use global_parameters, only: N, Ng, n_species, outfile, T0, do_load_IC, &
       ionization_transport, mol_ir_bands, rec_method, use_weno3, use_plm, &
       sec_ion_active, use_sec_ion, sec_ion_armed_step, resid_th, he_diffusion, &
       thereis_mol, carrier_transport, thereis_oxychem, r
  use species_table, only: n_mion, n_melem, melem_name, isp_H2, isp_H2p, &
       isp_H3p, isp_HeHp, isp_OH, isp_H2O, isp_CO
  use Read_input, only: input_read
  use Initialization, only: init
  use lya_rt, only: lya_rt_allocate_arrays
  use excited_hydrogen, only: excited_H_allocate_arrays
  use ionization_equilibrium, only: ioniz_eq_allocate_arrays, ieq_sweep_ledger_last
  use molecular_infrared_cooling, only: molecular_infrared_init
  use mol_rates, only: h2_thermochemistry_init
  use Conversion, only: U_to_W
  use BC_Apply, only: Apply_BC
  use composition, only: get_species_densities, comp_T_from_p
  use utils, only: calc_rho
  use steady_newton, only: set_transported_species_rows, neq_newton, &
       pack_U, pack_species_rows, eval_residual, solve_steady_jfnk, &
       read_species_unknown_space_controls
  use steady_residual_mod, only: assemble_residual, face_mass_flux_of_state, &
       n_cells_without_chemical_root
  use binary_element_diffusion, only: relax_element_composition, &
       element_transport_residual
  use diffusive_photochemistry, only: relax_photochemical_composition, &
       carrier_steady_residual, n_carrier_max, n_carrier, carrier_name, carrier_solved
  use certification, only: certification_evaluate, certification_report_write, &
       cert_report, cert_context_stationary, cert_evaluated, cert_regime_wind_r, cert_scale_floor
  use omp_lib, only: omp_get_wtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real*8, allocatable :: au(:,:), aw(:,:), af(:,:), afout(:,:), ar(:,:)
  real*8, allocatable :: ay(:), avec(:), ah(:), ac(:), anpart(:), atemp(:), aflux(:)
  real*8, allocatable :: ani(:), anii(:), aneutral(:), anion(:), anion2(:), atr(:)
  real*8, allocatable :: ane(:), antot(:), amet(:,:), amol(:,:), aoxy(:,:), arho(:)
  real*8, allocatable :: er(:), es(:), trr(:,:), trs(:,:), cr(:,:), cs(:,:)
  logical, allocatable :: carried(:)
  real*8 :: start_wall, step_wall, adrift, aomega, atrust, acmax, avol, closure_change
  integer :: pass, passes, maxit, ainfo, nstep, aj, ai, ak, unit_rows, unit_steps
  logical :: admissible, he_ok, tr_ok
  character(len=40) :: mode, arg, label
  type(cert_report) :: report

  mode = 'partition'; passes = 3; maxit = 40; aomega = 0.5d0; atrust = 1.0d-2
  call get_command_argument(1, arg)
  if (len_trim(arg)>0) mode = arg
  call get_command_argument(2, arg)
  if (len_trim(arg)>0) read(arg,*) passes
  call get_command_argument(3, arg)
  if (len_trim(arg)>0) read(arg,*) maxit
  call get_command_argument(4, arg)
  if (len_trim(arg)>0) read(arg,*) aomega
  call get_command_argument(5, arg)
  if (len_trim(arg)>0) read(arg,*) atrust
  if (mode/='partition' .and. mode/='coupled' .and. mode/='species') error stop 'unknown mode'
  start_wall = omp_get_wtime()
  call input_read
  if (.not.do_load_IC) error stop 'This experiment requires a saved initial state.'
  if (ionization_transport) error stop 'Transported proton outside this experiment scope.'
  open(unit=outfile,file='EXHALE_setup.out',status='replace')
  call lya_rt_allocate_arrays
  call excited_H_allocate_arrays
  call ioniz_eq_allocate_arrays
  allocate(au(3,1-Ng:N+Ng),aw(3,1-Ng:N+Ng),ar(3,1-Ng:N+Ng))
  allocate(af(1-Ng:N+Ng,n_species),afout(1-Ng:N+Ng,n_species))
  allocate(ah(1-Ng:N+Ng),ac(1-Ng:N+Ng),anpart(1-Ng:N+Ng),atemp(1-Ng:N+Ng),aflux(1-Ng:N+Ng))
  allocate(ani(1-Ng:N+Ng),anii(1-Ng:N+Ng),aneutral(1-Ng:N+Ng),anion(1-Ng:N+Ng))
  allocate(anion2(1-Ng:N+Ng),atr(1-Ng:N+Ng),ane(1-Ng:N+Ng),antot(1-Ng:N+Ng),arho(1-Ng:N+Ng))
  allocate(amet(1-Ng:N+Ng,n_mion),amol(1-Ng:N+Ng,4),aoxy(1-Ng:N+Ng,3))
  allocate(er(N),es(N),trr(N,n_melem),trs(N,n_melem),carried(n_melem))
  allocate(cr(N,n_carrier_max),cs(N,n_carrier_max))
  aneutral=0d0; anion=0d0; anion2=0d0; atr=0d0
  if (mol_ir_bands) call molecular_infrared_init(T0)
  call h2_thermochemistry_init
  call init(aw,au,af)
  rec_method='WENO3'; use_weno3=.true.; use_plm=.false.
  sec_ion_active=use_sec_ion
  if (sec_ion_active) sec_ion_armed_step=0
  open(newunit=unit_rows,file='rows.tsv',status='replace')
  write(unit_rows,'(A)') 'stage'//achar(9)//'row'//achar(9)//'status'//achar(9)// &
       'maximum'//achar(9)//'wind'//achar(9)//'reported_region'//achar(9)// &
       'gated'//achar(9)//'tolerance'//achar(9)//'pass'//achar(9)//'cell'
  open(newunit=unit_steps,file='steps.tsv',status='replace')
  write(unit_steps,'(A)') 'stage operation info steps drift seconds'
  write(*,'(A,A,A,I0,A,I0,A,ES12.4,A,ES12.4)') 'EXPERIMENT mode=',trim(mode), &
       ' passes=',passes,' maxit=',maxit,' omega=',aomega,' carrier_trust=',atrust
  call equilibrium_state_report('initial')

  if (mode=='coupled') then
    step_wall=omp_get_wtime()
    call set_transported_species_rows(.true.)
    call solve_steady_jfnk(au,af,resid_th,maxit,1d0,40,ainfo)
    write(unit_steps,*) 'coupled hydro_species ',ainfo,maxit,0d0,omp_get_wtime()-step_wall
    call equilibrium_state_report('coupled_return')
  else
    do pass=1,passes
      if (mode=='partition') then
        write(label,'(A,I2.2)') 'hydro_',pass
        call set_transported_species_rows(.false.)
        step_wall=omp_get_wtime()
        call solve_steady_jfnk(au,af,resid_th,maxit,1d0,40,ainfo)
        write(unit_steps,*) trim(label),' hydro ',ainfo,maxit,0d0,omp_get_wtime()-step_wall
        call equilibrium_state_report(trim(label))
        ! A failed subsolve does not certify this alternation. Continue only
        ! as a bounded diagnostic of whether its following composition update helps.
      endif
      call primitive_state
      call face_mass_flux_of_state(au(1,:),aflux)
      write(label,'(A,I2.2)') 'frozen_before_',pass
      call frozen_transport_report(trim(label))
      step_wall=omp_get_wtime()
      if (he_diffusion) then
        call relax_element_composition(aw(1,:),aw(2,:),atemp,af,aflux,aomega,adrift,nstep)
      else if (thereis_mol .and. carrier_transport) then
        call relax_photochemical_composition(aw(1,:),aw(2,:),af,atrust,adrift,nstep)
      else
        error stop 'No active transport block.'
      endif
      write(label,'(A,I2.2)') 'species_',pass
      write(unit_steps,*) trim(label),' species ',0,nstep,adrift,omp_get_wtime()-step_wall
      write(label,'(A,I2.2)') 'frozen_after_',pass
      call frozen_transport_report(trim(label))
      call composition_record(trim(label))
      flush(unit_steps)
      if (.not.all(ieee_is_finite(af))) error stop 'Nonfinite composition returned by transport'
      write(label,'(A,I2.2)') 'species_',pass
      call equilibrium_state_report(trim(label))
      flush(unit_rows); flush(unit_steps)
    enddo
  endif
  write(*,'(A,ES16.8)') 'EXPERIMENT total_seconds=',omp_get_wtime()-start_wall
  close(unit_rows); close(unit_steps); close(outfile)
contains
  subroutine composition_record(stage)
    character(len=*),intent(in) :: stage
    integer :: us,j
    write(*,'(A,A,A,I0,A,I0,A,ES23.15)') 'COMPOSITION ',stage, &
         ' nonfinite=',count(.not.ieee_is_finite(af)), &
         ' negative=',count(af<0d0),' maximum=',maxval(af)
    call get_species_densities(aw(1,:),af,ani,anii,aneutral,anion,anion2,atr,amet,ane,antot)
    amol=0d0; aoxy=0d0
    if (thereis_mol) then
      amol(:,1)=aw(1,:)*af(:,isp_H2); amol(:,2)=aw(1,:)*af(:,isp_H2p)
      amol(:,3)=aw(1,:)*af(:,isp_H3p); amol(:,4)=aw(1,:)*af(:,isp_HeHp)
    endif
    if (thereis_oxychem) then
      aoxy(:,1)=aw(1,:)*af(:,isp_OH); aoxy(:,2)=aw(1,:)*af(:,isp_H2O)
      aoxy(:,3)=aw(1,:)*af(:,isp_CO)
    endif
    call calc_rho(ani,anii,aneutral,anion,anion2,arho,amet,amol,aoxy)
    write(*,'(A,A,A,ES23.15)') 'COMPOSITION ',stage,' mass_closure=', &
         maxval(abs(arho(1:N)-aw(1,1:N))/max(abs(aw(1,1:N)),1d-300))
    open(newunit=us,file=trim(stage)//'.state',status='replace')
    write(us,'(A)') '# r u(1:3) frozen_T f_sp(1:n_species); before chemical refresh'
    do j=1-Ng,N+Ng
      write(us,'(*(ES24.16,1X))') r(j),au(:,j),atemp(j),af(j,:)
    enddo
    close(us)
    flush(6)
  end subroutine

  subroutine primitive_state
    call get_species_densities(au(1,:),af,ani,anii,aneutral,anion,anion2,atr,amet,ane,antot)
    call Apply_BC(au)
    call U_to_W(au,aw)
    call get_species_densities(aw(1,:),af,ani,anii,aneutral,anion,anion2,atr,amet,ane,antot)
    call comp_T_from_p(aw(3,:),antot,ane,atemp)
  end subroutine

  subroutine equilibrium_state_report(stage)
    character(len=*),intent(in) :: stage
    integer :: i,j,us
    real*8 :: mass_error, hydro_difference
    call set_transported_species_rows(.true.)
    call read_species_unknown_space_controls
    if (allocated(ay)) deallocate(ay,avec)
    allocate(ay(neq_newton()),avec(neq_newton()))
    call pack_U(au,ay)
    call pack_species_rows(au,af,ay)
    call eval_residual(ay,af,afout,avec,ah,ac,anpart,admissible=admissible)
    closure_change=maxval(abs(afout(1:N,:)-af(1:N,:)))
    af=afout
    call primitive_state
    ! Reassemble on the actual adopted composition and boundary, so the
    ! certification receives one state even if the chemical map moved it.
    call assemble_residual(au,antot+ane,ah,ac,ar)
    call face_mass_flux_of_state(au(1,:),aflux)
    hydro_difference=0d0
    do j=1,N
      do i=1,3
        hydro_difference=max(hydro_difference,abs(ar(i,j)-avec((neq_newton()/N)*(j-1)+i)))
      enddo
    enddo
    call certification_evaluate(cert_context_stationary,au,ar,af,resid_th, &
         n_cells_without_chemical_root(ieq_sweep_ledger_last%acc_n),.true.,report)
    call certification_report_write(report,stage)
    do i=1,report%n
      if (report%e(i)%status==0) cycle
      write(unit_rows,'(A,A,A,A,I0,5(A,ES23.15),A,L1,A,I0)') trim(stage),achar(9), &
           trim(report%e(i)%name),achar(9),report%e(i)%status,achar(9),report%e(i)%row_max, &
           achar(9),report%e(i)%row_max_wind,achar(9),report%e(i)%row_max_reported, &
           achar(9),report%e(i)%row_max_gate,achar(9),report%e(i)%tol, &
           achar(9),report%e(i)%within_tol,achar(9),report%e(i)%jworst
    enddo
    amol=0d0; aoxy=0d0
    if (thereis_mol) then
      amol(:,1)=aw(1,:)*af(:,isp_H2); amol(:,2)=aw(1,:)*af(:,isp_H2p)
      amol(:,3)=aw(1,:)*af(:,isp_H3p); amol(:,4)=aw(1,:)*af(:,isp_HeHp)
    endif
    if (thereis_oxychem) then
      aoxy(:,1)=aw(1,:)*af(:,isp_OH); aoxy(:,2)=aw(1,:)*af(:,isp_H2O)
      aoxy(:,3)=aw(1,:)*af(:,isp_CO)
    endif
    call calc_rho(ani,anii,aneutral,anion,anion2,arho,amet,amol,aoxy)
    mass_error=maxval(abs(arho(1:N)-aw(1,1:N))/max(abs(aw(1,1:N)),1d-300))
    write(*,'(A,A,A,L1,A,L1,4(A,ES16.8))') 'STATE ',stage,' admissible=',admissible, &
         ' certified=',report%certified,' mass_closure=',mass_error, &
         ' minimum_fraction=',minval(af(1:N,:)),' closure_map_change=',closure_change, &
         ' hydro_reassembly_difference=',hydro_difference
    if (.not.all(ieee_is_finite(au)) .or. .not.all(ieee_is_finite(af))) error stop 'Nonfinite state'
    open(newunit=us,file=trim(stage)//'.state',status='replace')
    write(us,'(A)') '# r u(1:3) T f_sp(1:n_species); code units'
    do j=1-Ng,N+Ng
      write(us,'(*(ES24.16,1X))') r(j),au(:,j),atemp(j),af(j,:)
    enddo
    close(us)
    open(newunit=us,file=trim(stage)//'.faces',status='replace')
    write(us,'(A)') '# face_index mass_flux; code units, physical faces 0:N'
    do j=0,N
      write(us,'(I0,1X,ES24.16)') j,aflux(j)
    enddo
    close(us)
    flush(unit_rows)
  end subroutine

  subroutine frozen_transport_report(stage)
    character(len=*),intent(in) :: stage
    integer :: i
    if (he_diffusion) then
      call element_transport_residual(aw(1,:),atemp,af,aflux,er,es,he_ok,trr,trs,tr_ok,carried)
      if (he_ok) call print_native(stage,'He',er,es)
      if (tr_ok) then
        do i=1,n_melem
          if (carried(i)) call print_native(stage,trim(melem_name(i)),trr(:,i),trs(:,i))
        enddo
      endif
    endif
    if (thereis_mol .and. carrier_transport) then
      call carrier_steady_residual(aw(1,:),aw(2,:),af,acmax,aj,ai,rvol=avol,res_out=cr,terms_out=cs)
      do i=1,n_carrier
        if (carrier_solved(i)) call print_native(stage,trim(carrier_name(i)),cr(:,i),cs(:,i))
      enddo
    endif
  end subroutine

  subroutine print_native(stage,name,res,scale)
    character(len=*),intent(in) :: stage,name
    real*8,intent(in) :: res(N),scale(N)
    real*8 :: norm(N)
    norm=abs(res)/max(scale,cert_scale_floor)
    write(*,'(A,A,1X,A,2(A,ES23.15))') 'FROZEN ',stage,name, &
         ' whole=',maxval(norm),' wind=',maxval(norm,mask=r(1:N)>=cert_regime_wind_r)
  end subroutine
end program transport_wind_experiment
