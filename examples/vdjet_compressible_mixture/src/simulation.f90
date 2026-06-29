!> Low-Mach variable-temperature vdjet case with a single-component NASG EOS.
!> This uses the public src/variable_density lowmach_class/vdscalar_class API.
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs, ss
   type(timetracker), public :: time

   type(ensight) :: ens_out
   type(event)   :: ens_evt

   type(monitor) :: mfile, cflfile, consfile

   public :: simulation_init, simulation_run, simulation_final

   ! Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU, resV, resW, resSC
   real(WP), dimension(:,:,:), allocatable :: Ui, Vi, Wi

   ! Output/property arrays
   real(WP), dimension(:,:,:), allocatable :: mu_mix, cp_mix, cv_mix, gamma_mix
   real(WP), dimension(:,:,:), allocatable :: speed_of_sound, Mach_mix

   ! Single-component NASG and transport parameters
   real(WP) :: P_inf, T_inf, T_cof_inlet, T_jet_inlet
   real(WP) :: gamma_ref, visc_ref, R_ref, Pref_ref, q_ref, b_ref
   real(WP) :: cv_ref, cp_ref, Pr_ref, kappa_ref

   ! Inlet parameters
   real(WP) :: Djet, Ujet, Ucof, xjet
   real(WP) :: x0_init, inlet_velocity_radius, inlet_ramp_time

   ! Integral of pressure residual
   real(WP) :: int_RP=0.0_WP

contains

   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = (i.eq.pg%imin)
   end function xm_locator

   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = (i.eq.pg%imin-1)
   end function xm_locator_sc

   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = (i.eq.pg%imax+1)
   end function xp_locator

   subroutine nasg_properties(pres,temp,rho,mu,cp,cv,gamma,a)
      implicit none
      real(WP), intent(in)  :: pres, temp
      real(WP), intent(out) :: rho, mu, cp, cv, gamma, a
      real(WP) :: peff, tcell, bulkmod

      peff    = max(pres+Pref_ref,1.0e-20_WP)
      tcell   = max(temp,1.0e-6_WP)
      rho     = 1.0_WP/max(b_ref+R_ref*tcell/peff,1.0e-20_WP)
      mu      = visc_ref
      cp      = cp_ref
      cv      = cv_ref
      gamma   = gamma_ref
      bulkmod = gamma*peff/max(1.0_WP-rho*b_ref,1.0e-20_WP)
      a       = sqrt(max(bulkmod/max(rho,1.0e-20_WP),0.0_WP))
   end subroutine nasg_properties

   function nasg_density(pres,temp) result(rho)
      implicit none
      real(WP), intent(in) :: pres, temp
      real(WP) :: rho, mu, cp, cv, gamma, a
      call nasg_properties(pres,temp,rho,mu,cp,cv,gamma,a)
   end function nasg_density

   ! Update density, viscosity, and diagnostic NASG properties from temperature.
   subroutine update_properties()
      implicit none
      integer :: i,j,k
      real(WP) :: temp, rho, mu, cp, cv, gamma, a

      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               temp = max(sc%SC(i,j,k),1.0e-6_WP)
               call nasg_properties(P_inf,temp,rho,mu,cp,cv,gamma,a)

               sc%rho(i,j,k)          = rho
               sc%diff(i,j,k)         = kappa_ref/max(cp_ref,1.0e-20_WP)
               fs%visc(i,j,k)         = mu
               mu_mix(i,j,k)          = mu
               cp_mix(i,j,k)          = cp
               cv_mix(i,j,k)          = cv
               gamma_mix(i,j,k)       = gamma
               speed_of_sound(i,j,k)  = a
            end do
         end do
      end do

      call sc%cfg%sync(sc%rho)
      call sc%cfg%sync(sc%diff)
      call fs%cfg%sync(fs%visc)
      call fs%cfg%sync(mu_mix)
      call fs%cfg%sync(cp_mix)
      call fs%cfg%sync(cv_mix)
      call fs%cfg%sync(gamma_mix)
      call fs%cfg%sync(speed_of_sound)
   end subroutine update_properties

   subroutine update_mach()
      implicit none
      integer :: i,j,k
      real(WP) :: velmag

      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               velmag = sqrt(Ui(i,j,k)**2+Vi(i,j,k)**2+Wi(i,j,k)**2)
               Mach_mix(i,j,k) = velmag/max(speed_of_sound(i,j,k),1.0e-20_WP)
            end do
         end do
      end do
      call cfg%sync(Mach_mix)
   end subroutine update_mach

   function plane_jet_blend(xloc,rloc) result(blend)
      implicit none
      real(WP), intent(in) :: xloc, rloc
      real(WP) :: blend
      real(WP) :: rho_ref, nu_ref, Uexcess, Mref, xeff
      real(WP) :: arg, arg_edge, amplitude, sech2, sech2_edge, denom, radius

      rho_ref = nasg_density(P_inf,T_cof_inlet)
      nu_ref  = visc_ref/max(rho_ref,1.0e-20_WP)
      Uexcess = max(abs(Ujet-Ucof),1.0e-20_WP)
      Mref    = rho_ref*Uexcess*Uexcess*Djet
      xeff    = max(xloc+x0_init,x0_init)
      radius  = max(inlet_velocity_radius,1.0e-20_WP)

      if (abs(rloc).ge.radius) then
         blend = 0.0_WP
         return
      end if

      arg      = abs(rloc)*(Mref/(48.0_WP*rho_ref*nu_ref*nu_ref))**(1.0_WP/3.0_WP)*xeff**(-2.0_WP/3.0_WP)
      arg_edge = radius   *(Mref/(48.0_WP*rho_ref*nu_ref*nu_ref))**(1.0_WP/3.0_WP)*xeff**(-2.0_WP/3.0_WP)
      amplitude  = (xeff/x0_init)**(-1.0_WP/3.0_WP)
      sech2      = 1.0_WP/cosh(arg     )**2
      sech2_edge = 1.0_WP/cosh(arg_edge)**2
      denom      = max(1.0_WP-sech2_edge,epsilon(1.0_WP))
      blend      = amplitude*(sech2-sech2_edge)/denom
      blend      = max(0.0_WP,min(1.0_WP,blend))
   end function plane_jet_blend

   subroutine get_inlet_profile(j,k,Tcell,uvel)
      implicit none
      integer, intent(in) :: j,k
      real(WP), intent(out) :: Tcell, uvel
      real(WP) :: r_dist, blend, ramp

      r_dist    = sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
      blend     = plane_jet_blend(0.0_WP,r_dist)
      ramp      = min(1.0_WP,time%t/max(inlet_ramp_time,1.0e-20_WP))
      uvel      = Ucof + ramp*(Ujet-Ucof)*blend
      Tcell     = T_cof_inlet + (T_jet_inlet-T_cof_inlet)*blend
   end subroutine get_inlet_profile

   subroutine enforce_scalar_inlet()
      use vdscalar_class, only: bcond
      implicit none
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      real(WP) :: Tcell, udum

      call sc%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         call get_inlet_profile(j,k,Tcell,udum)
         sc%SC(i,j,k) = Tcell
      end do
   end subroutine enforce_scalar_inlet

   subroutine enforce_velocity_inlet()
      use lowmach_class, only: bcond
      implicit none
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      real(WP) :: Tcell, uvel

      call fs%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         call get_inlet_profile(j,k,Tcell,uvel)
         fs%U(i,j,k)    = uvel
         fs%rhoU(i,j,k) = uvel*sum(fs%itpr_x(:,i,j,k)*fs%rho(i-1:i,j,k))
      end do
   end subroutine enforce_velocity_inlet

   subroutine simulation_init
      use param, only: param_read
      implicit none

      read_params: block
         call param_read('NASG gamma',     gamma_ref)
         call param_read('NASG Pref',      Pref_ref)
         call param_read('NASG q',         q_ref)
         call param_read('NASG b',         b_ref)
         call param_read('NASG viscosity', visc_ref)
         call param_read('NASG constant',  R_ref)
         call param_read('NASG Prandtl',   Pr_ref, default=0.72_WP)

         cv_ref    = R_ref/max(gamma_ref-1.0_WP,1.0e-20_WP)
         cp_ref    = cv_ref + R_ref
         kappa_ref = visc_ref*cp_ref/max(Pr_ref,1.0e-20_WP)

         call param_read('U jet',                 Ujet)
         call param_read('Jet diameter',          Djet)
         call param_read('Jet location',          xjet)
         call param_read('Initial jet offset',    x0_init, default=Djet)
         call param_read('Inlet velocity radius', inlet_velocity_radius, default=0.5_WP*cfg%yL)
         call param_read('Coflow velocity',       Ucof)
         call param_read('Inlet ramp time',       inlet_ramp_time, default=1.0e-3_WP)
         call param_read('Static pressure',       P_inf)
         call param_read('Static temperature',    T_inf)
         call param_read('Coflow temperature',    T_cof_inlet, default=T_inf)
         call param_read('Jet temperature',       T_jet_inlet, default=T_inf)
      end block read_params

      create_velocity_solver: block
         use lowmach_class,   only: dirichlet, clipped_neumann
         use hypre_str_class, only: pcg_pfmg2

         fs=lowmach(cfg=cfg,name='NASG property low Mach NS')
         call fs%add_bcond(name='inflow' , type=dirichlet,       face='x', dir=-1, &
                           canCorrect=.false., locator=xm_locator)
         call fs%add_bcond(name='outflow', type=clipped_neumann, face='x', dir=+1, &
                           canCorrect=.true.,  locator=xp_locator)
         ps = hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel = 18
         call param_read('Pressure iteration', ps%maxit)
         call param_read('Pressure tolerance', ps%rcvg)
         vs = ddadi(cfg=cfg,name='Velocity',nst=7)
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_velocity_solver

      create_scalar_solver: block
         use vdscalar_class, only: dirichlet, neumann, quick

         sc=vdscalar(cfg=cfg,scheme=quick,name='Temperature')
         call sc%add_bcond(name='inflow' , type=dirichlet, locator=xm_locator_sc)
         call sc%add_bcond(name='outflow', type=neumann,   locator=xp_locator, dir='+x')
         ss = ddadi(cfg=cfg,name='Scalar',nst=13)
         call sc%setup(implicit_solver=ss)
      end block create_scalar_solver

      allocate_work_arrays: block
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resSC(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
         allocate(mu_mix        (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(cp_mix        (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(cv_mix        (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gamma_mix     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(speed_of_sound(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Mach_mix      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot)
         call param_read('Max timestep size', time%dtmax)
         call param_read('Max cfl number',    time%cflmax)
         call param_read('Max time',          time%tmax)
         time%dt = time%dtmax
         call param_read('Subiterations', time%itmax, default=2)
      end block initialize_timetracker

      initialize_scalar: block
         sc%SC = T_cof_inlet
         call enforce_scalar_inlet()
         call sc%apply_bcond(time%t,time%dt)
         call update_properties()
         sc%rhoold = sc%rho
         sc%SCold = sc%SC
         call sc%rho_multiply()
      end block initialize_scalar

      initialize_velocity: block
         fs%rho = sc%rho
         fs%rhoold = fs%rho
         fs%U = Ucof
         fs%V = 0.0_WP
         fs%W = 0.0_WP
         call fs%apply_bcond(time%t,time%dt)
         call enforce_velocity_inlet()
         call fs%rho_multiply()
         fs%Uold = fs%U
         fs%Vold = fs%V
         fs%Wold = fs%W
         fs%rhoUold = fs%rhoU
         fs%rhoVold = fs%rhoV
         fs%rhoWold = fs%rhoW
         call fs%interp_vel(Ui,Vi,Wi)
         resSC = 0.0_WP
         call fs%get_div(drhodt=resSC)
         call fs%get_mfr()
         call update_mach()
      end block initialize_velocity

      create_ensight: block
         ens_out = ensight(cfg=cfg,name='GasJet_lowmach')
         ens_evt = event(time=time,name='Ensight output')
         call param_read('Ensight output period', ens_evt%tper)
         call ens_out%add_vector('velocity',         Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',         fs%P)
         call ens_out%add_scalar('density',          sc%rho)
         call ens_out%add_scalar('temperature',      sc%SC)
         call ens_out%add_scalar('viscosity',        mu_mix)
         call ens_out%add_scalar('cp_mix',           cp_mix)
         call ens_out%add_scalar('cv_mix',           cv_mix)
         call ens_out%add_scalar('gamma_mix',        gamma_mix)
         call ens_out%add_scalar('speed_of_sound',   speed_of_sound)
         call ens_out%add_scalar('Mach',             Mach_mix)
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      create_monitor: block
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         mfile = monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,        'Step')
         call mfile%add_column(time%t,        'Time')
         call mfile%add_column(time%dt,       'dt')
         call mfile%add_column(time%cfl,      'CFL_max')
         call mfile%add_column(fs%Umax,       'Umax')
         call mfile%add_column(fs%Vmax,       'Vmax')
         call mfile%add_column(fs%Wmax,       'Wmax')
         call mfile%add_column(fs%Pmax,       'Pmax')
         call mfile%add_column(sc%SCmax,      'Tmax')
         call mfile%add_column(sc%SCmin,      'Tmin')
         call mfile%add_column(sc%rhomax,     'RHOmax')
         call mfile%add_column(sc%rhomin,     'RHOmin')
         call mfile%add_column(int_RP,        'Int(RP)')
         call mfile%add_column(fs%divmax,     'divmax')
         call mfile%add_column(fs%psolv%it,   'P_iter')
         call mfile%add_column(fs%psolv%rerr, 'P_err')
         call mfile%write()

         cflfile = monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,    'Step')
         call cflfile%add_column(time%t,    'Time')
         call cflfile%add_column(fs%CFLc_x, 'CFL_conv_x')
         call cflfile%add_column(fs%CFLc_y, 'CFL_conv_y')
         call cflfile%add_column(fs%CFLc_z, 'CFL_conv_z')
         call cflfile%add_column(fs%CFLv_x, 'CFL_visc_x')
         call cflfile%add_column(fs%CFLv_y, 'CFL_visc_y')
         call cflfile%add_column(fs%CFLv_z, 'CFL_visc_z')
         call cflfile%write()

         consfile = monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,        'Step')
         call consfile%add_column(time%t,        'Time')
         call consfile%add_column(sc%SCint,      'T_integral')
         call consfile%add_column(sc%rhoint,     'RHO_integral')
         call consfile%add_column(sc%rhoSCint,   'rhoT_integral')
         call consfile%write()
      end block create_monitor
   end subroutine simulation_init

   subroutine simulation_run
      implicit none

      do while (.not.time%done())
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         sc%rhoold = sc%rho
         sc%SCold  = sc%SC

         fs%rhoold  = fs%rho
         fs%Uold    = fs%U
         fs%Vold    = fs%V
         fs%Wold    = fs%W
         fs%rhoUold = fs%rhoU
         fs%rhoVold = fs%rhoV
         fs%rhoWold = fs%rhoW

         do while (time%it.le.time%itmax)
            ! Temperature equation.
            sc%SC = 0.5_WP*(sc%SC+sc%SCold)
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            resSC = time%dt*resSC - (2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            sc%SC = 2.0_WP*sc%SC - sc%SCold + resSC
            call sc%apply_bcond(time%t,time%dt)
            call enforce_scalar_inlet()

            ! NASG properties at fixed thermodynamic pressure.
            call update_properties()
            call sc%rho_multiply()

            ! Momentum predictor.
            fs%rho  = 0.5_WP*(sc%rho+sc%rhoold)
            fs%U    = 0.5_WP*(fs%U+fs%Uold)
            fs%V    = 0.5_WP*(fs%V+fs%Vold)
            fs%W    = 0.5_WP*(fs%W+fs%Wold)
            fs%rhoU = 0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%rhoV = 0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%rhoW = 0.5_WP*(fs%rhoW+fs%rhoWold)

            call fs%get_dmomdt(resU,resV,resW)
            resU = time%dtmid*resU - (2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV = time%dtmid*resV - (2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW = time%dtmid*resW - (2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            fs%U = 2.0_WP*fs%U - fs%Uold + resU
            fs%V = 2.0_WP*fs%V - fs%Vold + resV
            fs%W = 2.0_WP*fs%W - fs%Wold + resW
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)
            call enforce_velocity_inlet()

            ! Pressure projection with variable density divergence source.
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%correct_mfr(drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs = -fs%cfg%vol*fs%div/time%dtmid
            call cfg%integrate(A=fs%psolv%rhs,integral=int_RP)
            fs%psolv%sol = 0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P = fs%P + fs%psolv%sol
            fs%rhoU = fs%rhoU - time%dtmid*resU
            fs%rhoV = fs%rhoV - time%dtmid*resV
            fs%rhoW = fs%rhoW - time%dtmid*resW
            call fs%rho_divide()

            time%it = time%it + 1
         end do

         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)
         call update_mach()

         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
      end do
   end subroutine simulation_run

   subroutine simulation_final
      implicit none
      deallocate(resU,resV,resW,resSC,Ui,Vi,Wi)
      deallocate(mu_mix,cp_mix,cv_mix,gamma_mix,speed_of_sound,Mach_mix)
   end subroutine simulation_final

end module simulation
