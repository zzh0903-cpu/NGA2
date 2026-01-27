 !> Various definitions and tools for running an NGA2 simulation
! TEST 2026-01-14
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use matm_class,        only: matm
   use hypre_str_class,   only: hypre_str
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile, cvgfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   integer :: relax_model
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   real(WP) :: P_inf        ! Infinity Pressure
   real(WP) :: T_inf        ! background temperature
   real(WP) :: T_jet        ! jet temperature
   !> Inlet parameters
   real(WP) :: Djet,Dcof
   real(WP) :: Ujet,Ucof
   real(WP) :: thick
   real(WP) :: L_sp        
   real(WP) :: A_sp         ! sponge
   !> Liquid Properties 
   real(WP) :: gas_const_R
   real(WP) :: rho_l_ref    ! Liquid reference density
   real(WP) :: LP,GP
   real(WP) :: xjet
contains

   ! Locate the Circular Jet Inlet
   function jet_inlet_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: r_dist, radius
      isIn = .false.
      if (i == pg%imin) then
         r_dist = sqrt(pg%ym(j)**2 + pg%zm(k)**2)
         radius = 0.5_WP * Djet 
         if (r_dist <= radius) isIn = .true.
      end if
   end function jet_inlet_locator

   ! Locate the Solid Wall/Coflow region
   function wall_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin-1) then
         if (.not. jet_inlet_locator(pg,i,j,k)) isIn = .true.
      end if
   end function wall_locator

   ! X+ (Outflow) boundary locator
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i == pg%imax) isIn = .true. 
   end function xp_locator

   ! Y- (Bottom) boundary locator
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmin) isIn = .true.
   end function ym_locator

   ! Y+ (Top) boundary locator
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmax) isIn = .true.
   end function yp_locator

   ! Z- (Back)boundary locator
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmin) isIn = .true.
   end function zm_locator

   !Z+ (Front) boundary locator
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmax) isIn = .true.
   end function zp_locator
   

   ! subroutine apply_sponges(dt)
   !    use mathtools, only: Pi
   !    implicit none
   !    real(WP), intent(in) :: dt
   !    integer :: i,j,k
   !    logical :: in_sponge
   !    real(WP) :: sigma,dist,weight,damp_factor
   !    ! Target State Variables
   !    real(WP) :: P_target, rho_g_target, rho_l_target, U_target
      
   !    ! Set target values (Must match initialization/ambient values)
   !    P_target     = P_inf
   !    rho_g_target = P_inf / (gas_const_R * T_inf) 
   !    rho_l_target = 1000.0_WP  ! Reference liquid density
   !    U_target     = Ucof       ! Coflow velocity (often 0.0)
      
   !    ! Note: A_sp (Strength) and L_sp (Length) are global variables defined in the module

   !    do k = fs%cfg%kmin_, fs%cfg%kmax_
   !       do j = fs%cfg%jmin_, fs%cfg%jmax_
   !          do i = fs%cfg%imin_, fs%cfg%imax_
               
   !             in_sponge = .false.
   !             sigma = 0.0_WP
               
   !             ! --- Calculate Damping Strength (Quadratic Profile) ---
               
   !             ! X+ Outflow Sponge
   !             if (fs%cfg%xm(i) > (fs%cfg%xL - L_sp)) then
   !                in_sponge = .true.
   !                dist = (fs%cfg%xm(i) - (fs%cfg%xL - L_sp)) / L_sp
   !                sigma = max(sigma, A_sp * dist**2)
   !             end if
               
   !             ! Y Direction Sponge (Top/Bottom)
   !             if (abs(fs%cfg%ym(j)) > (fs%cfg%yL/2.0_WP - L_sp)) then
   !                in_sponge = .true.
   !                dist = (abs(fs%cfg%ym(j)) - (fs%cfg%yL/2.0_WP - L_sp)) / L_sp
   !                sigma = max(sigma, A_sp * dist**2)
   !             end if

   !             ! Z Direction Sponge (Left/Right)
   !             if (abs(fs%cfg%zm(k)) > (fs%cfg%zL/2.0_WP - L_sp)) then
   !                in_sponge = .true.
   !                dist = (abs(fs%cfg%zm(k)) - (fs%cfg%zL/2.0_WP - L_sp)) / L_sp
   !                sigma = max(sigma, A_sp * dist**2)
   !             end if

   !             ! --- Apply Damping ---
   !             if (in_sponge) then
   !                weight = dt * sigma
   !                ! Prevent weight from becoming too large (numerical stability)
   !                ! damp_factor = 1.0 / (1.0 + weight)
                  
   !                ! Relax variables towards target state:
   !                ! New = Old - Weight * (Old - Target)
                  
   !                ! Correct Density (Gas and Liquid separately)
   !                fs%Grho(i,j,k) = fs%Grho(i,j,k) - weight * (fs%Grho(i,j,k) - rho_g_target)
   !                fs%Lrho(i,j,k) = fs%Lrho(i,j,k) - weight * (fs%Lrho(i,j,k) - rho_l_target)
                  
   !                ! Correct Pressure
   !                fs%GP(i,j,k)   = fs%GP(i,j,k)   - weight * (fs%GP(i,j,k)   - P_target)
   !                fs%LP(i,j,k)   = fs%LP(i,j,k)   - weight * (fs%LP(i,j,k)   - P_target)
                  
   !                ! Correct Velocity
   !                fs%Ui(i,j,k)   = fs%Ui(i,j,k)   - weight * (fs%Ui(i,j,k)   - U_target)
   !                fs%Vi(i,j,k)   = fs%Vi(i,j,k)   - weight * (fs%Vi(i,j,k)   - 0.0_WP)
   !                fs%Wi(i,j,k)   = fs%Wi(i,j,k)   - weight * (fs%Wi(i,j,k)   - 0.0_WP)
                  
   !                ! [IMPORTANT] Update Total Energy via EOS to maintain thermodynamic consistency
   !                fs%GrhoE(i,j,k) = matmod%EOS_energy(fs%GP(i,j,k), fs%Grho(i,j,k), &
   !                                                    fs%Ui(i,j,k), fs%Vi(i,j,k), fs%Wi(i,j,k), 'gas')
   !                fs%LrhoE(i,j,k) = matmod%EOS_energy(fs%LP(i,j,k), fs%Lrho(i,j,k), &
   !                                                    fs%Ui(i,j,k), fs%Vi(i,j,k), fs%Wi(i,j,k), 'liquid')
   !             end if
               
   !          end do
   !       end do
   !    end do

   ! end subroutine apply_sponges



   ! Function that defines a level set function for a initial wall_jet
   function levelset_jet_stub(xyz, t) result(G)
      implicit none
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      integer :: n
      G = 0.5_WP*Djet - sqrt( (xyz(1)-xjet)**2 + xyz(2)**2 + xyz(3)**2 )
   end function levelset_jet_stub
   
   

   !> Initialization of problem solver for Liquid Jet
   subroutine simulation_init
      use param, only: param_read, param_getsize
      implicit none
      real(WP) :: radius
      integer :: i,j,k

      initialize_timetracker: block
         time = timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size', time%dtmax)
         call param_read('Max cfl number', time%cflmax)
         call param_read('Max time', time%tmax)
         time%dt = time%dtmax 
         time%itmax = 2 
      end block initialize_timetracker

      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: elvira,VFhi,VFlo,remap
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with elvira reconstruction
         call vf%initialize(cfg=cfg,reconstruction_method=elvira,transport_method=remap,name='VOF')
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_jet_stub,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  ! Round up to fully liquid within wall (accuracy of mdot improves with resolution)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Boundary conditions on VF are built into the mast solver
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set initial interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Create a two-phase flow solver
      create_and_initialize_flow_solver: block
         use mast_class,      only: clipped_neumann, dirichlet,mech_egy_mech_hhz,bc_scope,bcond
         use hypre_str_class, only: hypre_str, pcg_pfmg
         integer :: i,j,k,n
         real(WP) :: gamm_l,gamm_g,visc_g,visc_l,sigma_st,Pref_l 
         type(bcond), pointer :: mybc

         real(WP) :: Ma_set, P_tot, T_tot
         real(WP) :: P_stat, T_stat, c_sound, U_cof_calc, Ujet_val    
         matmod = matm(cfg=cfg, name='Liquid-gas models')

         ! Ideal Gas
         call param_read('Gas gamma', gamm_g)
         call param_read('Gas constant', gas_const_R)
         call param_read('Gas viscosity', visc_g)

         ! Liquid Properties (Stiffened Gas) 
         call param_read('Liquid gamma', gamm_l)
         call param_read('Liquid Pref', Pref_l) 
         call param_read('Liquid viscosity', visc_l)
         call param_read('Liquid density', rho_l_ref)

         ! Register equations of state
         call matmod%register_idealgas('gas', gamm_g)
         call matmod%register_stiffenedgas('liquid', gamm_l, Pref_l)
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)
         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('gas',    fs%Grho, fs%Ui, fs%Vi, fs%Wi, fs%GrhoE, fs%GP)
         call matmod%register_thermoflow_variables('liquid', fs%Lrho, fs%Ui, fs%Vi, fs%Wi, fs%LrhoE, fs%LP)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g, viscconst_liquid=visc_l)
         call param_read('Surface tension coefficient', fs%sigma, default=0.072_WP) 
      
         ! Pressure Solver
         ps = hypre_str(cfg=cfg, name='Pressure', method=pcg_pfmg, nst=7)
         ps%maxlevel = 10 
         call param_read('Pressure iteration', ps%maxit)
         call param_read('Pressure tolerance', ps%rcvg)
         
         ! Implicit Velocity Solver
         vs = hypre_str(cfg=cfg, name='Velocity', method=pcg_pfmg, nst=7)
         call param_read('Implicit iteration', vs%maxit)
         call param_read('Implicit tolerance', vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps, implicit_solver=vs)
         
         ! Initialize liquid jet(s)
         call param_read('Mach number', Ma_set)
         call param_read('Stagnation pressure', P_tot)
         call param_read('Stagnation temperature', T_tot)
         call param_read('U jet', Ujet_val) 
         call param_read('Jet diameter', Djet)
         call param_read('Jet location',xjet)

         ! Calculate Isentropic Relations
         ! Static Temperature: T = T0 / (1 + (g-1)/2 * M^2)
         T_stat = T_tot / (1.0_WP + (gamm_g - 1.0_WP) * 0.5_WP * Ma_set**2)
         
         ! Static Pressure: P = P0 * (T/T0)^(g/(g-1))
         P_stat = P_tot * (T_stat / T_tot)**(gamm_g / (gamm_g - 1.0_WP))
         
         ! Sound Speed: c = sqrt(gamma * R * T)
         c_sound = sqrt(gamm_g * gas_const_R * T_stat)
         
         ! Coflow Velocity: U = M * c
         U_cof_calc = Ma_set * c_sound

         ! Store for global usage
         P_inf = P_stat
         T_inf = T_stat
         Ucof  = U_cof_calc 
         Ujet  = Ujet_val


         ! ! Read Sponge Params (Global variables)
         ! call param_read('Sponge length', L_sp, default=cfg%xL*0.1_WP) 
         ! call param_read('Sponge strength', A_sp, default=5.0_WP)

         ! C. Populate the Field
         do k = fs%cfg%kmino_, fs%cfg%kmaxo_
            do j = fs%cfg%jmino_, fs%cfg%jmaxo_
               do i = fs%cfg%imino_, fs%cfg%imaxo_

                  ! ! VOF
                  ! fs%VF(i,j,k)=vf%VF(i,j,k)
                  ! Velocity
                  if (vf%VF(i,j,k) > 0.5_WP) then
                     fs%Ui(i,j,k) = Ujet ! Liquid core velocity
                  else
                     fs%Ui(i,j,k) = Ucof ! Calculated high-speed gas velocity
                  end if
                  fs%Vi(i,j,k) = 0.0_WP
                  fs%Wi(i,j,k) = 0.0_WP

                  ! Density, Energy
                  fs%Grho(i,j,k)  = P_stat / (gas_const_R * T_stat)
                  fs%GrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Grho(i,j,k), &
                                                      fs%Ui(i,j,k), 0.0_WP, 0.0_WP, 'gas')
                  
                  ! Liquid: Reference Density + Stiffened Gas EOS
                  fs%Lrho(i,j,k)  = rho_l_ref 
                  fs%LrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Lrho(i,j,k), &
                                                      fs%Ui(i,j,k), 0.0_WP, 0.0_WP, 'liquid')
      
               end do
            end do
         end do

         ! Boundary Conditions 
         
         ! Jet Inlet (X- center) 
         call fs%add_bcond(name='jet_inlet', type=dirichlet, locator=jet_inlet_locator, celldir='xm')
         ! Coflow Inlet (X- outside) 
         call fs%add_bcond(name='coflow', type=dirichlet, locator=wall_locator, celldir='xm')
         ! 3. Open Boundaries (Outflow/Sides) -> Uses Sponge Layer
         call fs%add_bcond(name='outflow', type=clipped_neumann, locator=xp_locator, celldir='xp')
         call fs%add_bcond(name='bottom',  type=clipped_neumann, locator=ym_locator, celldir='ym')
         call fs%add_bcond(name='top',     type=clipped_neumann, locator=yp_locator, celldir='yp')
         call fs%add_bcond(name='left',    type=clipped_neumann, locator=zm_locator, celldir='zm')
         call fs%add_bcond(name='right',   type=clipped_neumann, locator=zp_locator, celldir='zp')
         ! Calculate face velocities
         call fs%interp_vel_basic(vf, fs%Ui, fs%Vi, fs%Wi, fs%U, fs%V, fs%W)
         ! Apply face BC - air inflow
         call fs%get_bcond('coflow', mybc) 
         do n = 1, mybc%itr%n_
            i = mybc%itr%map(1,n)
            j = mybc%itr%map(2,n)
            k = mybc%itr%map(3,n)
            fs%U(i:i+1, j, k) = fs%Ui(i, j, k) 
         end do
         ! Apply face BC - water inflow
         call fs%get_bcond('jet_inlet', mybc) 
         do n = 1, mybc%itr%n_
            i = mybc%itr%map(1,n)
            j = mybc%itr%map(2,n)
            k = mybc%itr%map(3,n)
            fs%U(i:i+1, j, k) = fs%Ui(i, j, k)
         end do
         ! Apply face BC - outflows
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt, bc_scope)

         ! RHO = alpha * rho_l + (1-alpha) * rho_g
         fs%RHO   = (1.0_WP - vf%VF) * fs%Grho + vf%VF * fs%Lrho
         fs%rhoUi = fs%RHO * fs%Ui
         fs%rhoVi = fs%RHO * fs%Vi
         fs%rhoWi = fs%RHO * fs%Wi

         ! Initial Relaxation
         relax_model = mech_egy_mech_hhz
         call fs%pressure_relax(vf, matmod, mech_egy_mech_hhz)
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
         fs%psolv%sol = 0.0_WP

         call matmod%update_temperature(vf,fs%Tmptr)
         print *, "DEBUG INIT: T_stat=", T_stat, " P_stat=", P_stat, " U_cof=", U_cof_calc

      end block create_and_initialize_flow_solver
      
      create_ensight: block
         ens_out = ensight(cfg=cfg, name='LiquidJet')
         ens_evt = event(time=time, name='Ensight output')
         call param_read('Ensight output period', ens_evt%tper)
         
         call ens_out%add_vector('velocity', fs%Ui, fs%Vi, fs%Wi)
         call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('Pressure', fs%P)
         call ens_out%add_scalar('Density', fs%RHO)
         call ens_out%add_scalar('VOF', vf%VF)     
         call ens_out%add_scalar('Temperature', fs%Tmptr)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('curvature',vf%curv)
         
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
         
      create_monitor: block
         call fs%get_cfl(time%dt, time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Simulation Log
         mfile = monitor(fs%cfg%amRoot, 'simulation')
         call mfile%add_column(time%n,  'Step')
         call mfile%add_column(time%t,  'Time')
         call mfile%add_column(time%dt, 'dt')
         call mfile%add_column(time%cfl,'CFL_max') 
         call mfile%add_column(fs%Umax, 'U_max')
         call mfile%add_column(fs%Pmax, 'P_max')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         
         ! CFL Breakdown (Watch Acoustic CFL!)
         cflfile = monitor(fs%cfg%amRoot, 'cfl')
         call cflfile%add_column(time%n, 'Step')
         call cflfile%add_column(fs%CFLc_x, 'CFL_conv_x')
         call cflfile%add_column(fs%CFLa_x, 'CFL_acous_x') ! Critical for compressible stability
         call cflfile%write()

         ! Solver Convergence
         cvgfile = monitor(fs%cfg%amRoot, 'cvg')
         call cvgfile%add_column(time%n, 'Step')
         call cvgfile%add_column(time%it,'SubIter')
         call cvgfile%add_column(fs%psolv%it,   'P_iter') 
         call cvgfile%add_column(fs%psolv%rerr, 'P_err')
         call cvgfile%write()
      end block create_monitor

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()
         
         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Insert viscous step here, or possibly incorporate into predictor above
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! ! Perform sponge forcing
            ! call apply_sponges(time%t)
            
            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()


            fs%psolv%sol=fs%PA-fs%P


            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)

            ! Record convergence monitor
            call cvgfile%write()
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
          ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         fs%PA = matmod%EOS_all(vf);
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call fs%get_viz()
         call mfile%write()
         call cflfile%write()
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
     
   end subroutine simulation_final
   

end module simulation