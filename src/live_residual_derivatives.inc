      subroutine live_residual_derivatives
      ! A fixed four-pass live thermochemical map, with converged boundary
      ! chemistry at every sample. This diagnoses directional differences;
      ! it does not certify the eliminated physical composition or a root.
      integer, parameter :: chemical_passes=4
      integer :: neq_live,unit_s,unit_v,unit_i,exponent,jl,kl,offset,idx
      real*8, allocatable :: y0(:),yp(:),ym(:),f0(:),fp(:),fm(:),fheld(:)
      real*8, allocatable :: direction(:),dscale(:),rscale(:)
      real*8, allocatable :: dfp(:),dfm(:),dfc(:),previous(:)
      real*8 :: hp,hm,hreq,asym(3),change(3),dnorm(3),reference
      real*8 :: fseed(1-Ng:N+Ng,n_species),fbase(1-Ng:N+Ng,n_species)
      real*8 :: ftrial(1-Ng:N+Ng,n_species)
      logical :: usable_base,usable_p,usable_m,have_previous
      neq_live=neq_newton()
      allocate(y0(neq_live),yp(neq_live),ym(neq_live),f0(neq_live), &
           fp(neq_live),fm(neq_live),fheld(neq_live),direction(neq_live), &
           dscale(neq_live),rscale(neq_live),dfp(neq_live),dfm(neq_live), &
           dfc(neq_live),previous(neq_live))
      call pack_U(u,y0)
      call pack_species_rows(u,f_sp,y0)
      call freeze_species_unknown_box(y0)
      fseed=f_sp
      call eval_residual(y0,fseed,fbase,f0,heat,cool, &
           admissible=usable_base,may_be_adopted=.false., &
           n_eq_sweeps_fixed=chemical_passes)
      if (.not.usable_base) then
         write(*,'(A)') ' ERROR: live derivative baseline is inadmissible'
         stop 1
      endif
      ! Reassemble the held returned state through the final-write operation.
      ! Physical U is unchanged; physical f is the map's returned composition.
      f_sp=fbase
      call sources_and_ghosts_of_the_held_composition
      call assemble_residual(u,n_tot+ne,f_sp,heat,cool,Rres)
      call write_stationary_residual_record('live_newton_baseline')
      call evaluation_state_dump('live_newton_baseline',u,f_sp,p,T)
      fheld=0.0d0
      do jl=1,N
         do kl=1,3
            fheld(nvar_jac*(jl-1)+kl)=Rres(kl,jl)
         enddo
      enddo
      open(newunit=unit_i,file=trim(output_dir)//'live_baseline_identity.txt',status='replace')
      write(unit_i,'(A)') '# cell slot live_R held_R difference (code units)'
      do jl=1,N
         do kl=1,3
            idx=nvar_jac*(jl-1)+kl
            write(unit_i,'(2I8,3ES24.15E3)') jl,kl,f0(idx),fheld(idx), &
                                              f0(idx)-fheld(idx)
         enddo
      enddo
      close(unit_i)
      call cell_state_scales(y0,dscale)
      call cell_row_scales(y0,dscale,rscale,u)
      do jl=1,neq_live
         direction(jl)=sin(0.1d0*dble(jl))
      enddo
      direction=direction/maxval(abs(direction))*dscale
      open(newunit=unit_s,file=trim(output_dir)//'live_derivative_steps.txt',status='replace')
      write(unit_s,'(A,I0)') '# fixed physical composition sweeps: ',chemical_passes
      write(unit_s,'(A)') '# h hp hm usable_plus usable_minus '// &
           'asym_mass asym_mom asym_energy change_mass change_mom '// &
           'change_energy norm_mass norm_mom norm_energy'
      open(newunit=unit_v,file=trim(output_dir)//'live_derivative_vectors.txt',status='replace')
      write(unit_v,'(A)') '# h cell slot central forward backward (row scaled)'
      have_previous=.false.
      do exponent=4,10
         hreq=10.0d0**(-exponent); hp=hreq; hm=hreq
         call gas_admissible_perturbation(y0,direction,hp,yp,'live forward')
         call gas_admissible_perturbation(y0,-direction,hm,ym,'live backward')
         call eval_residual(yp,fseed,ftrial,fp,heat0,cool0, &
              admissible=usable_p,may_be_adopted=.false., &
              n_eq_sweeps_fixed=chemical_passes,state_is_discarded=.true.)
         call eval_residual(ym,fseed,ftrial,fm,heat0,cool0, &
              admissible=usable_m,may_be_adopted=.false., &
              n_eq_sweeps_fixed=chemical_passes,state_is_discarded=.true.)
         if (.not.usable_p.or..not.usable_m) then
            write(*,'(A)') ' ERROR: live derivative has an inadmissible sample'
            stop 1
         endif
         dfp=(fp-f0)/hp/max(rscale,1.0d-300)
         dfm=(f0-fm)/hm/max(rscale,1.0d-300)
         dfc=(hm*dfp+hp*dfm)/(hp+hm)
         do kl=1,3
            dnorm(kl)=maxval(abs(dfc(kl:neq_live:nvar_jac)))
            reference=max(dnorm(kl),1.0d-300)
            asym(kl)=maxval(abs(dfp(kl:neq_live:nvar_jac) &
                              -dfm(kl:neq_live:nvar_jac)))/reference
            change(kl)=-1.0d0
            if (have_previous) change(kl)= &
                 maxval(abs(dfc(kl:neq_live:nvar_jac) &
                            -previous(kl:neq_live:nvar_jac)))/reference
         enddo
         write(unit_s,'(3ES24.15E3,2I4,9ES24.15E3)') hreq,hp,hm, &
              merge(1,0,usable_p),merge(1,0,usable_m),asym,change,dnorm
         do jl=1,N
            do kl=1,3
               offset=nvar_jac*(jl-1)+kl
               write(unit_v,'(ES24.15E3,2I8,3ES24.15E3)') hreq,jl,kl, &
                    dfc(offset),dfp(offset),dfm(offset)
            enddo
         enddo
         previous=dfc; have_previous=.true.
      enddo
      close(unit_s); close(unit_v)
      write(*,'(A)') ' (live_derivative) seven live directional samples completed'
      end subroutine live_residual_derivatives
