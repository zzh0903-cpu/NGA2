!> Low-Mach two-component heptane/N2 vdjet case with a NASG EOS.
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
   type(vdscalar),    public :: scT, scY
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs, ssT, ssY
   type(timetracker), public :: time

   type(ensight) :: ens_out
   type(event)   :: ens_evt

   type(monitor) :: mfile, cflfile, consfile

   public :: simulation_init, simulation_run, simulation_final

   ! Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU, resV, resW, resT, resY, resRho
   real(WP), dimension(:,:,:), allocatable :: Ui, Vi, Wi

   ! Output/property arrays
   real(WP), dimension(:,:,:), allocatable :: Xmole_heptane
   real(WP), dimension(:,:,:), allocatable :: mu_mix, cp_mix, cv_mix, gamma_mix
   real(WP), dimension(:,:,:), allocatable :: speed_of_sound, Mach_mix

   ! NASG and transport parameters
   real(WP) :: MW_heptane, MW_N2
   real(WP) :: P_inf, T_inf, T_cof_inlet, T_jet_inlet
   real(WP) :: Y_cof_inlet, Y_jet_inlet, D_species
   real(WP) :: gamm_n2_ref, visc_n2_ref, R_n2_ref, Pref_n2_ref, q_n2_ref, b_n2_ref
   real(WP) :: gamm_heptane_ref, visc_heptane_ref, R_heptane_ref, Pref_heptane_ref
   real(WP) :: q_heptane_ref, b_heptane_ref
   real(WP) :: cv_n2_ref, cv_heptane_ref, Pr_n2_ref, Pr_heptane_ref

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

   subroutine strict_mixture_terms(Y_heptane,pres,Av,Bv,Rm,cvm,cpm,gammam,Prefm,qm,bm,mum,Prm)
      implicit none
      real(WP), intent(in)  :: Y_heptane, pres
      real(WP), intent(out) :: Av, Bv, Rm, cvm, cpm, gammam, Prefm, qm, bm, mum, Prm
      real(WP) :: yh, yn2, den_h, den_n2

      yh  = max(0.0_WP,min(1.0_WP,Y_heptane))
      yn2 = 1.0_WP-yh

      den_h  = max(pres+Pref_heptane_ref,1.0e-20_WP)
      den_n2 = max(pres+Pref_n2_ref,     1.0e-20_WP)

      Av = yn2*R_n2_ref/den_n2 + yh*R_heptane_ref/den_h
      Bv = yn2*b_n2_ref        + yh*b_heptane_ref

      Rm     = yn2*R_n2_ref       + yh*R_heptane_ref
      cvm    = yn2*cv_n2_ref      + yh*cv_heptane_ref
      cpm    = cvm + Rm
      gammam = cpm/max(cvm,1.0e-20_WP)
      Prefm  = yn2*Pref_n2_ref    + yh*Pref_heptane_ref
      qm     = yn2*q_n2_ref       + yh*q_heptane_ref
      bm     = Bv
      mum    = yn2*visc_n2_ref    + yh*visc_heptane_ref
      Prm    = yn2*Pr_n2_ref      + yh*Pr_heptane_ref
   end subroutine strict_mixture_terms

   subroutine nasg_properties(pres,temp,Y_heptane,rho,mu,cp,cv,gamma,a,kappa)
      implicit none
      real(WP), intent(in)  :: pres, temp, Y_heptane
      real(WP), intent(out) :: rho, mu, cp, cv, gamma, a, kappa
      real(WP) :: Av, Bv, Rm, Prefm, qm, bm, Prm, bulkmod, v

      call strict_mixture_terms(Y_heptane,pres,Av,Bv,Rm,cv,cp,gamma,Prefm,qm,bm,mu,Prm)
      v       = Bv + max(temp,1.0e-6_WP)*Av
      rho     = 1.0_WP/max(v,1.0e-20_WP)
      kappa   = mu*cp/max(Prm,1.0e-20_WP)
      bulkmod = gamma*max(pres+Prefm,1.0e-20_WP)/max(1.0_WP-rho*bm,1.0e-20_WP)
      a       = sqrt(max(bulkmod/max(rho,1.0e-20_WP),0.0_WP))
   end subroutine nasg_properties

   function nasg_density(pres,temp,Y_heptane) result(rho)
      implicit none
      real(WP), intent(in) :: pres, temp, Y_heptane
      real(WP) :: rho, mu, cp, cv, gamma, a, kappa
      call nasg_properties(pres,temp,Y_heptane,rho,mu,cp,cv,gamma,a,kappa)
   end function nasg_density

   function heptane_mass_to_mole_frac(Y_heptane) result(X)
      implicit none
      real(WP), intent(in) :: Y_heptane
      real(WP) :: X, y, n_h, n_n2

      y    = max(0.0_WP,min(1.0_WP,Y_heptane))
      n_h  = y/max(MW_heptane,1.0e-20_WP)
      n_n2 = (1.0_WP-y)/max(MW_N2,1.0e-20_WP)
      X    = n_h/max(n_h+n_n2,1.0e-20_WP)
   end function heptane_mass_to_mole_frac

   ! Update density, viscosity, and diagnostic NASG properties from T and Y_heptane.
   subroutine update_properties()
      implicit none
      integer :: i,j,k
      real(WP) :: temp, y, rho, mu, cp, cv, gamma, a, kappa

      do k=scT%cfg%kmino_,scT%cfg%kmaxo_
         do j=scT%cfg%jmino_,scT%cfg%jmaxo_
            do i=scT%cfg%imino_,scT%cfg%imaxo_
               temp = max(scT%SC(i,j,k),1.0e-6_WP)
               y    = max(0.0_WP,min(1.0_WP,scY%SC(i,j,k)))
               call nasg_properties(P_inf,temp,y,rho,mu,cp,cv,gamma,a,kappa)

               scT%rho(i,j,k)         = rho
               scY%rho(i,j,k)         = rho
               scT%diff(i,j,k)        = kappa/max(cp,1.0e-20_WP)
               scY%diff(i,j,k)        = rho*D_species
               fs%visc(i,j,k)         = mu
               mu_mix(i,j,k)          = mu
               cp_mix(i,j,k)          = cp
               cv_mix(i,j,k)          = cv
               gamma_mix(i,j,k)       = gamma
               speed_of_sound(i,j,k)  = a
               Xmole_heptane(i,j,k)   = heptane_mass_to_mole_frac(y)
            end do
         end do
      end do

      call scT%cfg%sync(scT%rho)
      call scY%cfg%sync(scY%rho)
      call scT%cfg%sync(scT%diff)
      call scY%cfg%sync(scY%diff)
      call fs%cfg%sync(fs%visc)
      call fs%cfg%sync(mu_mix)
      call fs%cfg%sync(cp_mix)
      call fs%cfg%sync(cv_mix)
      call fs%cfg%sync(gamma_mix)
      call fs%cfg%sync(speed_of_sound)
      call fs%cfg%sync(Xmole_heptane)
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
      real(WP) :: rho_ref, mu_ref, cp_ref, cv_ref, gamma_ref, a_ref, kappa_ref
      real(WP) :: nu_ref, Uexcess, Mref, xeff
      real(WP) :: arg, arg_edge, amplitude, sech2, sech2_edge, denom, radius

      call nasg_properties(P_inf,T_cof_inlet,Y_cof_inlet,rho_ref,mu_ref,cp_ref,cv_ref,gamma_ref,a_ref,kappa_ref)
      nu_ref  = mu_ref/max(rho_ref,1.0e-20_WP)
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

   subroutine get_inlet_profile(j,k,Tcell,Y_heptane,uvel)
      implicit none
      integer, intent(in) :: j,k
      real(WP), intent(out) :: Tcell, Y_heptane, uvel
      real(WP) :: r_dist, blend, ramp

      r_dist    = sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
      blend     = plane_jet_blend(0.0_WP,r_dist)
      ramp      = min(1.0_WP,time%t/max(inlet_ramp_time,1.0e-20_WP))
      uvel      = Ucof + ramp*(Ujet-Ucof)*blend
      Tcell     = T_cof_inlet + (T_jet_inlet-T_cof_inlet)*blend
      Y_heptane = Y_cof_inlet + (Y_jet_inlet-Y_cof_inlet)*blend
   end subroutine get_inlet_profile

   subroutine enforce_scalar_inlet()
      use vdscalar_class, only: bcond
      implicit none
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      real(WP) :: Tcell, Y_h, udum

      call scT%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         call get_inlet_profile(j,k,Tcell,Y_h,udum)
         scT%SC(i,j,k) = Tcell
      end do

      call scY%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         call get_inlet_profile(j,k,Tcell,Y_h,udum)
         scY%SC(i,j,k) = Y_h
      end do
   end subroutine enforce_scalar_inlet

   subroutine enforce_velocity_inlet()
      use lowmach_class, only: bcond
      implicit none
      type(bcond), pointer :: mybc
      integer :: n,i,j,k
      real(WP) :: Tdum, Ydum, uvel

      call fs%get_bcond('inflow',mybc)
      do n=1,mybc%itr%no_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         call get_inlet_profile(j,k,Tdum,Ydum,uvel)
         fs%U(i,j,k)    = uvel
         fs%rhoU(i,j,k) = uvel*sum(fs%itpr_x(:,i,j,k)*fs%rho(i-1:i,j,k))
      end do
   end subroutine enforce_velocity_inlet

   subroutine simulation_init
      use param, only: param_read
      implicit none

      read_params: block
         call param_read('N2 gamma',     gamm_n2_ref)
         call param_read('N2 Pref',      Pref_n2_ref)
         call param_read('N2 q',         q_n2_ref)
         call param_read('N2 b',         b_n2_ref)
         call param_read('N2 viscosity', visc_n2_ref)
         call param_read('N2 constant',  R_n2_ref)
         call param_read('N2 MW',        MW_N2)
         call param_read('N2 Prandtl',   Pr_n2_ref, default=0.72_WP)

         call param_read('Heptane gamma',     gamm_heptane_ref)
         call param_read('Heptane Pref',      Pref_heptane_ref)
         call param_read('Heptane q',         q_heptane_ref)
         call param_read('Heptane b',         b_heptane_ref)
         call param_read('Heptane viscosity', visc_heptane_ref)
         call param_read('Heptane constant',  R_heptane_ref)
         call param_read('Heptane MW',        MW_heptane)
         call param_read('Heptane Prandtl',   Pr_heptane_ref, default=0.85_WP)

         cv_n2_ref      = R_n2_ref/max(gamm_n2_ref-1.0_WP,1.0e-20_WP)
         cv_heptane_ref = R_heptane_ref/max(gamm_heptane_ref-1.0_WP,1.0e-20_WP)

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
         call param_read('Coflow heptane mass fraction', Y_cof_inlet, default=0.0_WP)
         call param_read('Jet heptane mass fraction',    Y_jet_inlet, default=1.0_WP)
         call param_read('Species diffusivity',   D_species, default=1.0e-5_WP)
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

         scT=vdscalar(cfg=cfg,scheme=quick,name='Temperature')
         call scT%add_bcond(name='inflow' , type=dirichlet, locator=xm_locator_sc)
         call scT%add_bcond(name='outflow', type=neumann,   locator=xp_locator, dir='+x')
         ssT = ddadi(cfg=cfg,name='Temperature scalar',nst=13)
         call scT%setup(implicit_solver=ssT)

         scY=vdscalar(cfg=cfg,scheme=quick,name='MassFrac_heptane')
         call scY%add_bcond(name='inflow' , type=dirichlet, locator=xm_locator_sc)
         call scY%add_bcond(name='outflow', type=neumann,   locator=xp_locator, dir='+x')
         ssY = ddadi(cfg=cfg,name='Heptane scalar',nst=13)
         call scY%setup(implicit_solver=ssY)
      end block create_scalar_solver

      allocate_work_arrays: block
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resT  (scT%cfg%imino_:scT%cfg%imaxo_,scT%cfg%jmino_:scT%cfg%jmaxo_,scT%cfg%kmino_:scT%cfg%kmaxo_))
         allocate(resY  (scY%cfg%imino_:scY%cfg%imaxo_,scY%cfg%jmino_:scY%cfg%jmaxo_,scY%cfg%kmino_:scY%cfg%kmaxo_))
         allocate(resRho(scT%cfg%imino_:scT%cfg%imaxo_,scT%cfg%jmino_:scT%cfg%jmaxo_,scT%cfg%kmino_:scT%cfg%kmaxo_))
         allocate(Xmole_heptane(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         scT%SC = T_cof_inlet
         scY%SC = Y_cof_inlet
         call enforce_scalar_inlet()
         call scT%apply_bcond(time%t,time%dt)
         call scY%apply_bcond(time%t,time%dt)
         call update_properties()
         scT%rhoold = scT%rho
         scT%SCold  = scT%SC
         scY%rhoold = scY%rho
         scY%SCold  = scY%SC
         call scT%rho_multiply()
         call scY%rho_multiply()
      end block initialize_scalar

      initialize_velocity: block
         fs%rho = scT%rho
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
         resRho = 0.0_WP
         call fs%get_div(drhodt=resRho)
         call fs%get_mfr()
         call update_mach()
      end block initialize_velocity

      create_ensight: block
         ens_out = ensight(cfg=cfg,name='GasJet_lowmach')
         ens_evt = event(time=time,name='Ensight output')
         call param_read('Ensight output period', ens_evt%tper)
         call ens_out%add_vector('velocity',         Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',         fs%P)
         call ens_out%add_scalar('density',          scT%rho)
         call ens_out%add_scalar('temperature',      scT%SC)
         call ens_out%add_scalar('MassFrac_heptane', scY%SC)
         call ens_out%add_scalar('MoleFrac_heptane', Xmole_heptane)
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
         call scT%get_max()
         call scY%get_max()
         call scT%get_int()
         call scY%get_int()
         mfile = monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,        'Step')
         call mfile%add_column(time%t,        'Time')
         call mfile%add_column(time%dt,       'dt')
         call mfile%add_column(time%cfl,      'CFL_max')
         call mfile%add_column(fs%Umax,       'Umax')
         call mfile%add_column(fs%Vmax,       'Vmax')
         call mfile%add_column(fs%Wmax,       'Wmax')
         call mfile%add_column(fs%Pmax,       'Pmax')
         call mfile%add_column(scT%SCmax,     'Tmax')
         call mfile%add_column(scT%SCmin,     'Tmin')
         call mfile%add_column(scY%SCmax,     'Ymax')
         call mfile%add_column(scY%SCmin,     'Ymin')
         call mfile%add_column(scT%rhomax,    'RHOmax')
         call mfile%add_column(scT%rhomin,    'RHOmin')
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
         call consfile%add_column(scT%SCint,     'T_integral')
         call consfile%add_column(scY%SCint,     'Y_integral')
         call consfile%add_column(scT%rhoint,    'RHO_integral')
         call consfile%add_column(scT%rhoSCint,  'rhoT_integral')
         call consfile%add_column(scY%rhoSCint,  'rhoY_integral')
         call consfile%write()
      end block create_monitor
   end subroutine simulation_init

   subroutine simulation_run
      implicit none

      do while (.not.time%done())
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         scT%rhoold = scT%rho
         scT%SCold  = scT%SC
         scY%rhoold = scY%rho
         scY%SCold  = scY%SC

         fs%rhoold  = fs%rho
         fs%Uold    = fs%U
         fs%Vold    = fs%V
         fs%Wold    = fs%W
         fs%rhoUold = fs%rhoU
         fs%rhoVold = fs%rhoV
         fs%rhoWold = fs%rhoW

         do while (time%it.le.time%itmax)
            ! Heptane mass fraction equation.
            scY%SC = 0.5_WP*(scY%SC+scY%SCold)
            call scY%get_drhoSCdt(resY,fs%rhoU,fs%rhoV,fs%rhoW)
            resY = time%dt*resY - (2.0_WP*scY%rho*scY%SC-(scY%rho+scY%rhoold)*scY%SCold)
            call scY%solve_implicit(time%dt,resY,fs%rhoU,fs%rhoV,fs%rhoW)
            scY%SC = 2.0_WP*scY%SC - scY%SCold + resY
            scY%SC = max(0.0_WP,min(1.0_WP,scY%SC))

            ! Temperature equation.
            scT%SC = 0.5_WP*(scT%SC+scT%SCold)
            call scT%get_drhoSCdt(resT,fs%rhoU,fs%rhoV,fs%rhoW)
            resT = time%dt*resT - (2.0_WP*scT%rho*scT%SC-(scT%rho+scT%rhoold)*scT%SCold)
            call scT%solve_implicit(time%dt,resT,fs%rhoU,fs%rhoV,fs%rhoW)
            scT%SC = 2.0_WP*scT%SC - scT%SCold + resT

            call enforce_scalar_inlet()
            call scT%apply_bcond(time%t,time%dt)
            call scY%apply_bcond(time%t,time%dt)

            ! NASG properties at fixed thermodynamic pressure.
            call update_properties()
            call scT%rho_multiply()
            call scY%rho_multiply()

            ! Momentum predictor.
            fs%rho  = 0.5_WP*(scT%rho+scT%rhoold)
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
            call scT%get_drhodt(dt=time%dt,drhodt=resRho)
            call fs%correct_mfr(drhodt=resRho)
            call fs%get_div(drhodt=resRho)
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
         call scT%get_drhodt(dt=time%dt,drhodt=resRho)
         call fs%get_div(drhodt=resRho)
         call update_mach()

         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         call fs%get_max()
         call scT%get_max()
         call scY%get_max()
         call scT%get_int()
         call scY%get_int()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
      end do
   end subroutine simulation_run

   subroutine simulation_final
      implicit none
      deallocate(resU,resV,resW,resT,resY,resRho,Ui,Vi,Wi)
      deallocate(Xmole_heptane,mu_mix,cp_mix,cv_mix,gamma_mix,speed_of_sound,Mach_mix)
   end subroutine simulation_final

end module simulation
