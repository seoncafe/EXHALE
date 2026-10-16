      module regional_convergence
      ! Regional stabilization: finite-duration regional stability.
      ! Stationary convergence requires a separate stationary root test.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      private
      public :: regional_read, regional_begin, regional_sample, regional_incomplete
      public :: regional_enabled, regional_passed, regional_exhausted
      public :: regional_lower, regional_upper, regional_window
      public :: regional_profile_every
      logical :: regional_enabled=.false.,regional_passed=.false.
      logical :: regional_exhausted=.false.
      real*8 :: regional_lower=1.2d0,regional_upper=1.2d0
      real*8 :: regional_window=0d0,maximum_seconds=0d0,settling_seconds=0d0
      ! Proposed finite-duration limits: lower, transition, upper.
      ! Prioritize upper stability while permitting larger lower variations.
      real*8 :: mean_limit(3)=(/0.10d0,0.10d0,0.01d0/)
      real*8 :: range_limit(3)=(/0.20d0,0.20d0,0.05d0/)
      real*8 :: inventory_limit(3)=(/0.05d0,0.05d0,0.01d0/)
      real*8 :: flow_limit(3)=(/0.05d0,0.05d0,0.01d0/)
      integer :: validation_windows=4,confirmation_factor=2,minimum_samples=40
      integer,allocatable :: samples(:)
      integer :: regional_profile_every=40
      integer :: nwindow,iw=1,candidate=1,log_unit,accepted_samples=0
      integer,allocatable :: regions(:),fields(:)
      real*8,allocatable :: previous(:),reference(:),scales(:)
      real*8,allocatable :: integral(:,:),minimum(:,:),maximum(:,:)
      real*8 :: elapsed=0d0,within=0d0
      logical :: started=.false.
      contains

      subroutine regional_read(path)
      character(*),intent(in) :: path
      character(1024) :: line,key,value
      integer :: unit,ios,c,hash
      real*8 :: inventory_override(3),flow_override(3),common_limit
      logical :: inventory_given(3),flow_given(3)
      inventory_given=.false.;flow_given=.false.
      inventory_override=0d0;flow_override=0d0
      open(newunit=unit,file=path,status='old',action='read')
      do
         read(unit,'(A)',iostat=ios) line
         if(ios.lt.0) exit
         if(ios.ne.0) error stop 'cannot read regional convergence input'
         line=adjustl(line)
         if(index(line,'Regional ').ne.1) cycle
         hash=index(line,'#');if(hash.gt.0) line=line(:hash-1)
         c=index(line,':');if(c.eq.0) error stop 'regional input requires colon'
         key=trim(line(:c-1));value=trim(adjustl(line(c+1:)))
         ios=0
         select case(trim(key))
         case('Regional convergence')
            if(trim(value).ne.'True'.and.trim(value).ne.'False') &
               error stop 'Regional convergence requires True or False'
            regional_enabled=trim(value).eq.'True'
         case('Regional lower radius');read(value,*,iostat=ios) regional_lower
         case('Regional upper radius');read(value,*,iostat=ios) regional_upper
         case('Regional window seconds');read(value,*,iostat=ios) regional_window
         case('Regional maximum seconds');read(value,*,iostat=ios) maximum_seconds
         case('Regional settling seconds');read(value,*,iostat=ios) settling_seconds
         case('Regional validation windows');read(value,*,iostat=ios) validation_windows
         case('Regional confirmation factor');read(value,*,iostat=ios) confirmation_factor
         case('Regional lower mean limit');read(value,*,iostat=ios) mean_limit(1)
         case('Regional upper mean limit');read(value,*,iostat=ios) mean_limit(3)
         case('Regional transition mean limit');read(value,*,iostat=ios) mean_limit(2)
         case('Regional lower range limit');read(value,*,iostat=ios) range_limit(1)
         case('Regional upper range limit');read(value,*,iostat=ios) range_limit(3)
         case('Regional transition range limit');read(value,*,iostat=ios) range_limit(2)
         case('Regional inventory limit')
            read(value,*,iostat=ios) common_limit
            if(ios.eq.0) inventory_limit=common_limit
         case('Regional mass flow limit')
            read(value,*,iostat=ios) common_limit
            if(ios.eq.0) flow_limit=common_limit
         case('Regional lower inventory limit')
            read(value,*,iostat=ios) inventory_override(1);inventory_given(1)=.true.
         case('Regional transition inventory limit')
            read(value,*,iostat=ios) inventory_override(2);inventory_given(2)=.true.
         case('Regional upper inventory limit')
            read(value,*,iostat=ios) inventory_override(3);inventory_given(3)=.true.
         case('Regional lower mass flow limit')
            read(value,*,iostat=ios) flow_override(1);flow_given(1)=.true.
         case('Regional transition mass flow limit')
            read(value,*,iostat=ios) flow_override(2);flow_given(2)=.true.
         case('Regional upper mass flow limit')
            read(value,*,iostat=ios) flow_override(3);flow_given(3)=.true.
         case('Regional minimum samples');read(value,*,iostat=ios) minimum_samples
         case('Regional profile interval');read(value,*,iostat=ios) regional_profile_every
         case default
            error stop 'unknown Regional input key'
         end select
         if(ios.ne.0) error stop 'invalid Regional input value'
      enddo
      close(unit)
      ! Explicit regional values override a common value in either line order.
      where(inventory_given) inventory_limit=inventory_override
      where(flow_given) flow_limit=flow_override
      if(.not.regional_enabled) return
      if(.not.all(ieee_is_finite((/regional_lower,regional_upper,regional_window, &
            maximum_seconds,settling_seconds,mean_limit,range_limit, &
            inventory_limit,flow_limit/)))) error stop 'nonfinite regional input'
      if(regional_lower.gt.regional_upper.or.regional_lower.le.0d0) &
         error stop 'invalid regional radii'
      if(minval(mean_limit).le.0d0.or.minval(range_limit).le.0d0.or. &
         minval(inventory_limit).le.0d0.or.minval(flow_limit).le.0d0) error stop 'nonpositive regional limit'
      if(validation_windows.lt.4.or.validation_windows.gt.32.or. &
         confirmation_factor.lt.2.or.confirmation_factor.gt.4) &
         error stop 'regional windows require 4..32 and confirmation factor 2..4'
      if(regional_window.lt.0d0.or.maximum_seconds.lt.0d0.or.settling_seconds.lt.0d0.or. &
         regional_profile_every.lt.1.or.minimum_samples.lt.1) error stop 'invalid regional duration or cadence'
      end subroutine

      subroutine regional_begin(values,normalization,region,field,acoustic_seconds)
      real*8,intent(in) :: values(:),normalization(:),acoustic_seconds
      integer,intent(in) :: region(:),field(:)
      integer :: m
      if(.not.regional_enabled) return
      if(started) error stop 'regional convergence already started'
      if(regional_window.eq.0d0) regional_window=2d0*acoustic_seconds
      nwindow=validation_windows*confirmation_factor
      if(maximum_seconds.eq.0d0) &
         maximum_seconds=settling_seconds+2d0*nwindow*regional_window
      if(.not.ieee_is_finite(regional_window).or.regional_window.le.0d0.or. &
         .not.ieee_is_finite(maximum_seconds)) error stop 'invalid regional time scale'
      if(maximum_seconds.lt.settling_seconds+nwindow*regional_window) &
         error stop 'regional maximum time shorter than confirmation horizon'
      m=size(values)
      if(size(region).ne.m.or.size(field).ne.m.or.size(normalization).ne.m) &
         error stop 'regional vector size mismatch'
      if(any(region.lt.1).or.any(region.gt.3).or.any(field.lt.1).or.any(field.gt.49)) &
         error stop 'invalid regional quantity metadata'
      if(.not.any(region.eq.1).or..not.any(region.eq.3)) error stop 'empty lower or upper region'
      if(.not.all(ieee_is_finite(normalization)).or.any(normalization.le.0d0)) &
         error stop 'invalid regional reference scale'
      if(.not.all(ieee_is_finite(values))) error stop 'nonfinite regional initial state'
      allocate(previous(m),reference(m),scales(m),regions(m),fields(m))
      allocate(samples(nwindow))
      allocate(integral(m,nwindow),minimum(m,nwindow),maximum(m,nwindow))
      previous=values;reference=values;scales=normalization
      regions=region;fields=field
      call reset_candidate(values)
      ! Excluded startup states cannot bound the first retained window.
      ! Its extrema begin with the interpolated state at the settling edge.
      if(settling_seconds.gt.0d0) then
         minimum(:,1)=huge(1d0);maximum(:,1)=-huge(1d0)
      endif
      open(newunit=log_unit,file='output/regional_convergence.csv',status='replace')
      write(log_unit,'(A)') 'elapsed_seconds,candidate,windows,region,field,mean_change,range_width,'// &
         'mean_limit,range_limit,passed,minimum_samples_observed,minimum_samples_required'
      flush(log_unit)
      call status_record('collecting')
      started=.true.
      end subroutine

      subroutine reset_candidate(values)
      real*8,intent(in) :: values(:)
      samples=0
      integral=0d0;minimum=huge(1d0);maximum=-huge(1d0)
      minimum(:,1)=values;maximum(:,1)=values;iw=1;within=0d0
      end subroutine

      subroutine regional_sample(dt,values)
      ! Integrate linearly between accepted states, splitting at exact physical
      ! window boundaries. Every accepted step contributes to the extrema.
      real*8,intent(in) :: dt,values(:)
      real*8 :: a,b,h,remaining,end_fraction,x0(size(values)),x1(size(values))
      logical :: ok
      if(.not.regional_enabled) return
      if(.not.started) error stop 'regional convergence not initialized'
      if(.not.ieee_is_finite(dt).or.dt.le.0d0.or..not.all(ieee_is_finite(values))) &
         error stop 'nonfinite regional accepted sample'
      if(regional_passed.or.regional_exhausted) error stop 'sample after regional termination'
      a=0d0
      end_fraction=min(1d0,max(0d0,(maximum_seconds-elapsed)/dt))
      if(regional_window.lt.dt*end_fraction/(2d0*nwindow)) &
         error stop 'regional window is unresolved within an accepted interval'
      if(elapsed.lt.settling_seconds) a=min(1d0,(settling_seconds-elapsed)/dt)
      do while(a.lt.end_fraction)
         remaining=regional_window-within
         h=min(dt*(end_fraction-a),remaining);b=min(end_fraction,a+h/dt)
         if(b.le.a) error stop 'regional window is below numerical time resolution'
         x0=previous+a*(values-previous);x1=previous+b*(values-previous)
         ! The accepted RK mass flow is an interval average, not an endpoint.
         where(fields.eq.49)
            x0=values;x1=values
         end where
         integral(:,iw)=integral(:,iw)+0.5d0*h*(x0+x1)
         minimum(:,iw)=min(minimum(:,iw),x0,x1)
         maximum(:,iw)=max(maximum(:,iw),x0,x1)
         if(b.ge.1d0-1d-12) samples(iw)=samples(iw)+1
         within=within+h;a=b
         if(within.ge.regional_window*(1d0-1d-12)) then
            if(iw.eq.1) then
               reference=integral(:,1)/regional_window
               ! Positive fluid quantities and inventories use their FIRST
               ! window mean. Velocity and populations retain frozen scales.
               where(fields.ne.2.and.(fields.le.4.or.(fields.ge.47.and.fields.le.48))) &
                  scales=abs(reference)
               ! A zero mean flow has no relative scale. Retain a rounding-scale
               ! absolute denominator from the initial acoustic mass flow.
               where(fields.eq.49) scales=max(abs(reference),epsilon(1d0)*scales)
               if(any(scales.le.0d0)) error stop 'zero regional normalization'
            endif
            call evaluate_windows(elapsed+dt*a,ok)
            if(iw.eq.nwindow) then
               ! The final saved state is the adopted endpoint, even when a
               ! step crosses a window edge. Include it in the finite bound.
               minimum(:,iw)=min(minimum(:,iw),values)
               maximum(:,iw)=max(maximum(:,iw),values)
               call evaluate_windows(elapsed+dt,ok)
               ! The confirmation edge must be inside the allowed horizon.
               ! Its adopted endpoint can lie just beyond that edge, as for
               ! the ordinary time-budget stop; its range was checked above.
               if(ok.and.elapsed+dt*a.le.maximum_seconds*(1d0+1d-12)) then
                  regional_passed=.true.;exit
               endif
               candidate=candidate+1
               call reset_candidate(x1)
            else
               iw=iw+1;within=0d0;minimum(:,iw)=x1;maximum(:,iw)=x1
            endif
         endif
      enddo
      elapsed=elapsed+dt;previous=values
      accepted_samples=accepted_samples+1
      if(elapsed.ge.maximum_seconds.and..not.regional_passed) regional_exhausted=.true.
      if(regional_passed) then
         call status_record('regional_limits_satisfied')
      else if(regional_exhausted) then
         call status_record('time_budget_exhausted')
      else if(mod(accepted_samples,regional_profile_every).eq.0) then
         call status_record('collecting')
      endif
      end subroutine

      subroutine evaluate_windows(t,ok)
      real*8,intent(in) :: t
      logical,intent(out) :: ok
      real*8 :: dm(size(scales)),span(size(scales)),ml,rl,dmmax,spmax
      logical :: mask(size(scales)),passed
      integer :: reg,f,w
      dm=0d0
      do w=1,iw
         dm=max(dm,abs(integral(:,w)/regional_window-reference)/scales)
         if(w.gt.1) dm=max(dm,abs(integral(:,w)-integral(:,w-1))/regional_window/scales)
      enddo
      ! Means on opposite sides of the original reference can each remain
      ! within its limit while their total span exceeds that limit.
      dm=max(dm,(maxval(integral(:,:iw),dim=2)-minval(integral(:,:iw),dim=2))/regional_window/scales)
      span=(maxval(maximum(:,:iw),dim=2)-minval(minimum(:,:iw),dim=2))/scales
      ok=all(samples(:iw).ge.minimum_samples)
      do reg=1,3
         do f=1,49
            mask=regions.eq.reg.and.fields.eq.f
            if(.not.any(mask)) cycle
            ml=mean_limit(reg);rl=range_limit(reg)
            if(f.eq.47.or.f.eq.48) then
               ml=inventory_limit(reg);rl=inventory_limit(reg)
            else if(f.eq.49) then
               ml=flow_limit(reg);rl=huge(1d0)
            endif
            dmmax=maxval(dm,mask=mask);spmax=maxval(span,mask=mask)
            passed=dmmax.le.ml.and.spmax.le.rl.and.all(samples(:iw).ge.minimum_samples)
            ok=ok.and.passed
            write(log_unit,'(ES25.16E3,4(",",I0),4(",",ES25.16E3),3(",",I0))') &
               t,candidate,iw,reg,f,dmmax,spmax,ml,rl,merge(1,0,passed),minval(samples(:iw)),minimum_samples
         enddo
      enddo
      flush(log_unit)
      call status_record('collecting',t)
      end subroutine

      subroutine regional_incomplete(reason)
      character(*),intent(in) :: reason
      if(regional_enabled.and.started) call status_record(reason)
      end subroutine

      subroutine status_record(reason,measured_time)
      character(*),intent(in) :: reason
      real*8,intent(in),optional :: measured_time
      real*8 :: time
      integer :: unit
      time=elapsed
      if(present(measured_time)) time=measured_time
      ! A closed status file is read by the plotter. A concurrent reader may
      ! retry an incomplete replacement; a partial file never grants a stop.
      open(newunit=unit,file='output/regional_convergence_status.txt',status='replace')
      write(unit,'(A)') 'criterion=finite_duration_regional_stability'
      write(unit,'(A)') 'target_convergence_level=regional_stabilization'
      if(regional_passed) then
         write(unit,'(A)') 'convergence_level=regional_stabilization'
      else
         write(unit,'(A)') 'convergence_level=none'
      endif
      write(unit,'(A)') 'stationary_root_claim=False'
      write(unit,'(A)') 'status='//trim(reason)
      write(unit,'(A,ES24.16)') 'elapsed_seconds=',time
      write(unit,'(A,ES24.16)') 'window_seconds=',regional_window
      write(unit,'(A,ES24.16)') 'maximum_seconds=',maximum_seconds
      write(unit,'(A,I0)') 'confirmation_windows=',nwindow
      write(unit,'(A,I0)') 'candidate=',candidate
      write(unit,'(A,I0)') 'completed_windows=',merge(nwindow,iw-1,regional_passed)
      write(unit,'(A,ES24.16)') 'current_window_seconds=',within
      write(unit,'(A,I0)') 'required_samples_each_window=',minimum_samples
      write(unit,'(A,*(I0,1X))') 'samples_in_candidate_windows=',samples
      write(unit,'(A,ES24.16)') 'lower_radius_R0=',regional_lower
      write(unit,'(A,ES24.16)') 'upper_radius_R0=',regional_upper
      write(unit,'(A,3(ES24.16,1X))') 'mean_limits_lower_transition_upper=',mean_limit
      write(unit,'(A,3(ES24.16,1X))') 'range_limits_lower_transition_upper=',range_limit
      write(unit,'(A,3(ES24.16,1X))') 'inventory_limits_lower_transition_upper=',inventory_limit
      write(unit,'(A,3(ES24.16,1X))') 'mass_flow_limits_lower_transition_upper=',flow_limit
      close(unit)
      end subroutine
      end module regional_convergence
