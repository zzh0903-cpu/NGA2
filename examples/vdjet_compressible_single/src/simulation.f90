 !> Various definitions and tools for running an NGA2 simulation
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

   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs
   type(timetracker), public :: time

   type(ensight) :: ens_out
   type(event)   :: ens_evt

   type(monitor) :: mfile,cflfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   integer :: relax_model
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wit
   real(WP), dimension(:,:,:), allocatable :: Ymass_heptane, Ymass_heptane_old
   real(WP), dimension(:,:,:), allocatable :: rhoY_heptane, rhoY_heptane_old
   real(WP), dimension(:,:,:), allocatable :: Xmole_heptane
   real(WP), dimension(:,:,:), allocatable :: R_mix, gamma_mix, Pref_mix, q_mix, b_mix, mu_mix
   real(WP), dimension(:,:,:), allocatable :: nu_mix_offset
   real(WP), parameter :: MW_heptane = 100.2_WP    ! n-heptane
   real(WP), parameter :: MW_N2   = 28.0134_WP
   real(WP) :: x0_init
   real(WP) :: P_inf
   real(WP) :: T_inf
   real(WP) :: T_jet_inlet
   real(WP) :: T_cof_inlet
   !> Inlet parameters
   real(WP) :: Djet,Dcof
   real(WP) :: Ujet,Ucof
   real(WP) :: inlet_velocity_radius
   real(WP) :: L_sp
   real(WP) :: A_sp
   real(WP) :: xjet
   real(WP) :: D_species
   real(WP) :: gamm_n2_ref,visc_n2_ref,R_n2_ref,Pref_n2_ref,q_n2_ref,b_n2_ref
   real(WP) :: gamm_heptane_ref,visc_heptane_ref,R_heptane_ref,Pref_heptane_ref,q_heptane_ref,b_heptane_ref
   real(WP) :: Pr_n2_ref,Pr_heptane_ref
   real(WP) :: cv_n2_ref,cv_heptane_ref
   real(WP) :: k_n2_ref,k_heptane_ref
contains

   function jet_inlet_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      real(WP) :: r_dist, radius

      isIn = .false.
      if (i == pg%imin-1) then
         r_dist = sqrt(pg%ym(j)**2 + pg%zm(k)**2)
         radius = 0.5_WP*Djet
         if (r_dist <= radius) isIn = .true.
      end if
   end function jet_inlet_locator

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

   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i == pg%imax) isIn = .true.
   end function xp_locator

   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmin) isIn = .true.
   end function ym_locator

   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmax) isIn = .true.
   end function yp_locator

   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmin) isIn = .true.
   end function zm_locator

   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmax) isIn = .true.
   end function zp_locator

   ! function levelset_jet_stub(xyz, t) result(G)
   !    implicit none
   !    real(WP), dimension(3), intent(in) :: xyz
   !    real(WP), intent(in) :: t
   !    real(WP) :: G
   !    G = 0.5_WP*Djet - sqrt( (xyz(1)-xjet)**2 + xyz(2)**2 + xyz(3)**2 )
   ! end function levelset_jet_stub

   subroutine get_mixture_coeffs(Y_heptane,Rm,gammam,Prefm,qm,bm,mum,cvm,km)
      implicit none
      real(WP), intent(in)  :: Y_heptane
      real(WP), intent(out) :: Rm,gammam,Prefm,qm,bm,mum,cvm,km
      real(WP) :: y,yn2

      y = max(0.0_WP,min(1.0_WP,Y_heptane))
      yn2 = 1.0_WP-y
      Rm = yn2*R_n2_ref + y*R_heptane_ref
      cvm = yn2*cv_n2_ref + y*cv_heptane_ref
      gammam = 1.0_WP + Rm/max(cvm,1.0e-20_WP)
      Prefm = yn2*Pref_n2_ref + y*Pref_heptane_ref
      qm    = yn2*q_n2_ref    + y*q_heptane_ref
      bm    = yn2*b_n2_ref    + y*b_heptane_ref
      mum   = yn2*visc_n2_ref + y*visc_heptane_ref
      km    = yn2*k_n2_ref    + y*k_heptane_ref
   end subroutine get_mixture_coeffs

   function nasg_density(pres,temp,Rm,Prefm,bm) result(rho)
      implicit none
      real(WP), intent(in) :: pres,temp,Rm,Prefm,bm
      real(WP) :: rho,peff,teff,denom

      if (pres+Prefm.gt.1.0e-20_WP) then
         peff = pres+Prefm
      else
         peff = max(P_inf+Prefm,1.0e-20_WP)
      end if
      if (temp.gt.1.0e-6_WP) then
         teff = temp
      else
         teff = T_inf
      end if
      denom = Rm*teff/peff + bm
      if (denom.gt.1.0e-20_WP) then
         rho = 1.0_WP/denom
      else
         rho = 1.0_WP/1.0e-20_WP
      end if
   end function nasg_density

   function nasg_energy(pres,dens,uvel,vvel,wvel,gammam,Prefm,qm,bm) result(energy)
      implicit none
      real(WP), intent(in) :: pres,dens,uvel,vvel,wvel,gammam,Prefm,qm,bm
      real(WP) :: energy,KE

      KE = 0.5_WP*dens*(uvel**2+vvel**2+wvel**2)
      energy = (pres+gammam*Prefm)*(1.0_WP-dens*bm)/max(gammam-1.0_WP,1.0e-20_WP) &
             + dens*qm + KE
   end function nasg_energy

   function nasg_pressure(dens,temp,Rm,Prefm,bm) result(pres)
      implicit none
      real(WP), intent(in) :: dens,temp,Rm,Prefm,bm
      real(WP) :: pres,rhoeff,teff

      rhoeff = max(dens,1.0e-20_WP)
      teff = max(temp,1.0e-6_WP)
      pres = max(Rm*teff/max(1.0_WP/rhoeff-bm,1.0e-20_WP)-Prefm,1.0e-20_WP)
   end function nasg_pressure

   function nasg_temperature(pres,dens,Rm,Prefm,bm) result(temp)
      implicit none
      real(WP), intent(in) :: pres,dens,Rm,Prefm,bm
      real(WP) :: temp,peff,rhoeff

      peff = max(pres+Prefm,1.0e-20_WP)
      rhoeff = max(dens,1.0e-20_WP)
      temp = max((1.0_WP/rhoeff-bm)*peff/max(Rm,1.0e-20_WP),1.0e-6_WP)
   end function nasg_temperature

   function nasg_temperature_from_energy(energy,dens,uvel,vvel,wvel,cvm,Prefm,qm,bm) result(temp)
      implicit none
      real(WP), intent(in) :: energy,dens,uvel,vvel,wvel,cvm,Prefm,qm,bm
      real(WP) :: temp,KE,rhoeff

      rhoeff = max(dens,1.0e-20_WP)
      KE = 0.5_WP*rhoeff*(uvel**2+vvel**2+wvel**2)
      temp = max((energy-KE-rhoeff*qm-(1.0_WP-rhoeff*bm)*Prefm)/max(rhoeff*cvm,1.0e-20_WP),1.0e-6_WP)
   end function nasg_temperature_from_energy

   function nasg_bulkmod(pres,dens,gammam,Prefm,bm) result(bulkmod)
      implicit none
      real(WP), intent(in) :: pres,dens,gammam,Prefm,bm
      real(WP) :: bulkmod

      bulkmod = gammam*max(pres+Prefm,1.0e-20_WP)/max(1.0_WP-dens*bm,1.0e-20_WP)
   end function nasg_bulkmod

   function nasg_EOS_v(dens,gammam,bm) result(eosv)
      implicit none
      real(WP), intent(in) :: dens,gammam,bm
      real(WP) :: eosv

      eosv = (gammam-1.0_WP)/max(1.0_WP-dens*bm,1.0e-20_WP)
   end function nasg_EOS_v

   function heptane_mass_fraction_to_mole_fraction(Y_heptane) result(X_heptane)
      implicit none
      real(WP), intent(in) :: Y_heptane
      real(WP) :: X_heptane,y,n_heptane,n_N2

      y = max(0.0_WP,min(1.0_WP,Y_heptane))
      n_heptane = y/MW_heptane
      n_N2    = (1.0_WP-y)/MW_N2
      X_heptane = n_heptane/max(n_heptane+n_N2,1.0e-20_WP)
   end function heptane_mass_fraction_to_mole_fraction

   function plane_jet_blend(xloc,rloc) result(blend)
      implicit none
      real(WP), intent(in) :: xloc, rloc
      real(WP) :: blend
      real(WP) :: rho_ref, nu_ref, Uexcess, Mref, xeff, arg, arg_edge
      real(WP) :: amplitude, sech2, sech2_edge, denom, radius

      rho_ref = 1.0_WP / ( (R_n2_ref*T_cof_inlet)/(P_inf + Pref_n2_ref) + b_n2_ref )
      nu_ref  = visc_n2_ref / max(rho_ref,1.0e-20_WP)
      Uexcess = max(abs(Ujet-Ucof), 1.0e-20_WP)
      Mref    = rho_ref * Uexcess * Uexcess * Djet

      xeff = max(xloc + x0_init, x0_init)
      radius = max(inlet_velocity_radius, 1.0e-20_WP)

      if (abs(rloc).ge.radius) then
         blend = 0.0_WP
         return
      end if

      arg      = abs(rloc) * (Mref/(48.0_WP*rho_ref*nu_ref*nu_ref))**(1.0_WP/3.0_WP) * xeff**(-2.0_WP/3.0_WP)
      arg_edge = radius    * (Mref/(48.0_WP*rho_ref*nu_ref*nu_ref))**(1.0_WP/3.0_WP) * xeff**(-2.0_WP/3.0_WP)

      amplitude  = (xeff/x0_init)**(-1.0_WP/3.0_WP)
      sech2      = 1.0_WP / cosh(arg     )**2
      sech2_edge = 1.0_WP / cosh(arg_edge)**2
      denom      = max(1.0_WP-sech2_edge, epsilon(1.0_WP))

      blend = amplitude * (sech2-sech2_edge) / denom
      blend = max(0.0_WP, min(1.0_WP, blend))
   end function plane_jet_blend

   subroutine get_inlet_profile(j,k,Y_heptane,uvel,tcell)
      implicit none
      integer,  intent(in)  :: j,k
      real(WP), intent(out) :: Y_heptane,uvel,tcell
      real(WP) :: r_dist

      r_dist = sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
      uvel = Ucof + (Ujet - Ucof)*plane_jet_blend(0.0_WP, r_dist)
      ! Single-phase debug: pure N2 everywhere, uniform temperature
      Y_heptane = 0.0_WP
      tcell = T_cof_inlet
   end subroutine get_inlet_profile

   subroutine sync_mixture_fields(sync_velocity,sync_energy)
      implicit none
      logical, intent(in), optional :: sync_velocity,sync_energy
      logical :: do_velocity,do_energy

      do_velocity = .false.
      do_energy   = .false.
      if (present(sync_velocity)) do_velocity = sync_velocity
      if (present(sync_energy))   do_energy   = sync_energy

      call fs%cfg%sync(vf%VF)
      call fs%cfg%sync(Ymass_heptane)
      call fs%cfg%sync(rhoY_heptane)
      call fs%cfg%sync(Xmole_heptane)
      call fs%cfg%sync(R_mix)
      call fs%cfg%sync(gamma_mix)
      call fs%cfg%sync(Pref_mix)
      call fs%cfg%sync(q_mix)
      call fs%cfg%sync(b_mix)
      call fs%cfg%sync(mu_mix)
      call fs%cfg%sync(nu_mix_offset)
      call fs%cfg%sync(fs%RHO)
      call fs%cfg%sync(fs%Grho)
      call fs%cfg%sync(fs%Lrho)
      call fs%cfg%sync(fs%P)
      call fs%cfg%sync(fs%GP)
      call fs%cfg%sync(fs%LP)
      call fs%cfg%sync(fs%PA)
      call fs%cfg%sync(fs%GrhoSS2)
      call fs%cfg%sync(fs%LrhoSS2)
      call fs%cfg%sync(fs%RHOSS2)
      call fs%cfg%sync(fs%visc)
      call fs%cfg%sync(fs%therm_cond)
      call fs%cfg%sync(fs%rhoCv)
      call fs%cfg%sync(fs%EOS_v)
      call fs%cfg%sync(fs%Tmptr)
      if (do_velocity) then
         call fs%cfg%sync(fs%Ui)
         call fs%cfg%sync(fs%Vi)
         call fs%cfg%sync(fs%Wi)
         call fs%cfg%sync(fs%rhoUi)
         call fs%cfg%sync(fs%rhoVi)
         call fs%cfg%sync(fs%rhoWi)
      end if
      if (do_energy) then
         call fs%cfg%sync(fs%GrhoE)
         call fs%cfg%sync(fs%LrhoE)
      end if
   end subroutine sync_mixture_fields

   subroutine set_cell_mixture_state(i,j,k,Y_heptane,uvel,vvel,wvel,tcell,pcell)
      implicit none
      integer,  intent(in) :: i,j,k
      real(WP), intent(in) :: Y_heptane,uvel,vvel,wvel,tcell,pcell
      real(WP) :: y,rhocell,cvm,kcell,bulkcell,eosvcell

      y = max(0.0_WP,min(1.0_WP,Y_heptane))
      call get_mixture_coeffs(y,R_mix(i,j,k),gamma_mix(i,j,k),Pref_mix(i,j,k), &
                              q_mix(i,j,k),b_mix(i,j,k),mu_mix(i,j,k),cvm,kcell)

      rhocell = nasg_density(pcell,tcell,R_mix(i,j,k),Pref_mix(i,j,k),b_mix(i,j,k))
      bulkcell = nasg_bulkmod(pcell,rhocell,gamma_mix(i,j,k),Pref_mix(i,j,k),b_mix(i,j,k))
      eosvcell = nasg_EOS_v(rhocell,gamma_mix(i,j,k),b_mix(i,j,k))

      vf%VF(i,j,k) = 0.0_WP
      Ymass_heptane(i,j,k) = y
      rhoY_heptane(i,j,k) = rhocell*y
      Xmole_heptane(i,j,k) = heptane_mass_fraction_to_mole_fraction(y)
      nu_mix_offset(i,j,k) = 0.0_WP

      fs%Ui(i,j,k) = uvel
      fs%Vi(i,j,k) = vvel
      fs%Wi(i,j,k) = wvel
      fs%RHO(i,j,k) = rhocell
      fs%Grho(i,j,k) = rhocell
      fs%Lrho(i,j,k) = 0.0_WP
      fs%P(i,j,k) = pcell
      fs%GP(i,j,k) = pcell
      fs%LP(i,j,k) = 0.0_WP
      fs%PA(i,j,k) = pcell
      fs%GrhoE(i,j,k) = nasg_energy(pcell,rhocell,uvel,vvel,wvel,gamma_mix(i,j,k), &
                                    Pref_mix(i,j,k),q_mix(i,j,k),b_mix(i,j,k))
      fs%LrhoE(i,j,k) = 0.0_WP
      fs%rhoUi(i,j,k) = rhocell*uvel
      fs%rhoVi(i,j,k) = rhocell*vvel
      fs%rhoWi(i,j,k) = rhocell*wvel
      fs%GrhoSS2(i,j,k) = bulkcell
      fs%LrhoSS2(i,j,k) = 0.0_WP
      fs%RHOSS2(i,j,k) = bulkcell
      fs%visc(i,j,k) = mu_mix(i,j,k)
      fs%therm_cond(i,j,k) = kcell
      fs%rhoCv(i,j,k) = rhocell*cvm
      fs%EOS_v(i,j,k) = eosvcell
      fs%Tmptr(i,j,k) = tcell
   end subroutine set_cell_mixture_state

   subroutine apply_inlet_state()
      use mast_class, only: bc_scope
      implicit none
      integer :: i,j,k
      real(WP) :: ycell,uvel,tcell

      if (.not.allocated(Ymass_heptane)) return

      if (fs%cfg%imin_==fs%cfg%imin) then
         i = fs%cfg%imin_
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               call get_inlet_profile(j,k,ycell,uvel,tcell)
               call set_cell_mixture_state(i-1,j,k,ycell,uvel,0.0_WP,0.0_WP,tcell,P_inf)
               call set_cell_mixture_state(i  ,j,k,ycell,uvel,0.0_WP,0.0_WP,tcell,P_inf)
            end do
         end do
      end if

      call sync_mixture_fields(sync_velocity=.true.,sync_energy=.true.)
      call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
      bc_scope = 'velocity'
      call fs%apply_bcond(time%dt,bc_scope)
      call fs%interp_pressure_density(vf)
   end subroutine apply_inlet_state

   subroutine update_mixture_from_species(update_pressure)
      implicit none
      logical, intent(in), optional :: update_pressure
      logical :: set_pressure
      integer :: i,j,k
      real(WP) :: y,pcell,tcell,rhocell,ecell,bulkcell
      real(WP) :: cvcell,kcell,eosvcell

      if (.not.allocated(Ymass_heptane)) return
      set_pressure = .false.
      if (present(update_pressure)) set_pressure = update_pressure

      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               y = max(0.0_WP,min(1.0_WP,Ymass_heptane(i,j,k)))
               Ymass_heptane(i,j,k) = y
               call get_mixture_coeffs(y,R_mix(i,j,k),gamma_mix(i,j,k),Pref_mix(i,j,k), &
                                       q_mix(i,j,k),b_mix(i,j,k),mu_mix(i,j,k),cvcell,kcell)
               rhocell = max(fs%RHO(i,j,k),1.0e-20_WP)
               ecell = fs%GrhoE(i,j,k)
               tcell = nasg_temperature_from_energy(ecell,rhocell,fs%Ui(i,j,k),fs%Vi(i,j,k),fs%Wi(i,j,k), &
                                                    cvcell,Pref_mix(i,j,k),q_mix(i,j,k),b_mix(i,j,k))
               pcell = nasg_pressure(rhocell,tcell,R_mix(i,j,k),Pref_mix(i,j,k),b_mix(i,j,k))
               bulkcell = nasg_bulkmod(pcell,rhocell,gamma_mix(i,j,k),Pref_mix(i,j,k),b_mix(i,j,k))
               eosvcell = nasg_EOS_v(rhocell,gamma_mix(i,j,k),b_mix(i,j,k))

               vf%VF(i,j,k) = 0.0_WP
               fs%RHO(i,j,k) = rhocell
               fs%Grho(i,j,k) = rhocell
               fs%Lrho(i,j,k) = 0.0_WP
               fs%rhoUi(i,j,k) = rhocell*fs%Ui(i,j,k)
               fs%rhoVi(i,j,k) = rhocell*fs%Vi(i,j,k)
               fs%rhoWi(i,j,k) = rhocell*fs%Wi(i,j,k)
               fs%GP(i,j,k) = pcell
               fs%LP(i,j,k) = 0.0_WP
               fs%PA(i,j,k) = pcell
               if (set_pressure) fs%P(i,j,k) = pcell
               fs%GrhoSS2(i,j,k) = bulkcell
               fs%LrhoSS2(i,j,k) = 0.0_WP
               fs%RHOSS2(i,j,k) = bulkcell
               fs%visc(i,j,k) = mu_mix(i,j,k)
               fs%therm_cond(i,j,k) = kcell
               fs%rhoCv(i,j,k) = rhocell*cvcell
               fs%EOS_v(i,j,k) = eosvcell
               nu_mix_offset(i,j,k) = 0.0_WP
               rhoY_heptane(i,j,k) = rhocell*y
               Xmole_heptane(i,j,k) = heptane_mass_fraction_to_mole_fraction(y)
               fs%Tmptr(i,j,k) = tcell
            end do
         end do
      end do

      call sync_mixture_fields(sync_velocity=.true.)
   end subroutine update_mixture_from_species

   subroutine update_heptane_mole_fraction()
      implicit none
      integer :: i,j,k

      do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               Ymass_heptane(i,j,k) = max(0.0_WP,min(1.0_WP,Ymass_heptane(i,j,k)))
               rhoY_heptane(i,j,k) = fs%RHO(i,j,k)*Ymass_heptane(i,j,k)
               Xmole_heptane(i,j,k) = heptane_mass_fraction_to_mole_fraction(Ymass_heptane(i,j,k))
            end do
         end do
      end do
      call fs%cfg%sync(Ymass_heptane)
      call fs%cfg%sync(rhoY_heptane)
      call fs%cfg%sync(Xmole_heptane)
   end subroutine update_heptane_mole_fraction


   subroutine apply_species_bcond()
      implicit none
      integer :: i,j,k
      real(WP) :: ycell,dummy_uvel,dummy_tcell

      if (.not.allocated(Ymass_heptane)) return

      ! jet_inlet
      if (fs%cfg%imin_==fs%cfg%imin) then
         i = fs%cfg%imin_
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               call get_inlet_profile(j,k,ycell,dummy_uvel,dummy_tcell)
               Ymass_heptane(i-1,j,k) = ycell
               rhoY_heptane (i-1,j,k) = fs%RHO(i-1,j,k)*ycell
               Xmole_heptane(i-1,j,k) = heptane_mass_fraction_to_mole_fraction(ycell)
            end do
         end do
      end if
      ! xp_max
      if (fs%cfg%imax_==fs%cfg%imax) then
         i = fs%cfg%imax_
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               Ymass_heptane(i+1,j,k) = Ymass_heptane(i,j,k)
               rhoY_heptane (i+1,j,k) = rhoY_heptane (i,j,k)
               Xmole_heptane(i+1,j,k) = Xmole_heptane(i,j,k)
            end do
         end do
      end if
      ! ym_min and yp_max
      if (.not.fs%cfg%yper) then
         if (fs%cfg%jmin_==fs%cfg%jmin) then
            j = fs%cfg%jmin_
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Ymass_heptane(i,j-1,k) = Ymass_heptane(i,j,k)
                  rhoY_heptane (i,j-1,k) = rhoY_heptane (i,j,k)
                  Xmole_heptane(i,j-1,k) = Xmole_heptane(i,j,k)
               end do
            end do
         end if
         if (fs%cfg%jmax_==fs%cfg%jmax) then
            j = fs%cfg%jmax_
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Ymass_heptane(i,j+1,k) = Ymass_heptane(i,j,k)
                  rhoY_heptane (i,j+1,k) = rhoY_heptane (i,j,k)
                  Xmole_heptane(i,j+1,k) = Xmole_heptane(i,j,k)
               end do
            end do
         end if
      end if
      ! zm_min and zp_max
      if (.not.fs%cfg%zper) then
         if (fs%cfg%kmin_==fs%cfg%kmin) then
            k = fs%cfg%kmin_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Ymass_heptane(i,j,k-1) = Ymass_heptane(i,j,k)
                  rhoY_heptane (i,j,k-1) = rhoY_heptane (i,j,k)
                  Xmole_heptane(i,j,k-1) = Xmole_heptane(i,j,k)
               end do
            end do
         end if
         if (fs%cfg%kmax_==fs%cfg%kmax) then
            k = fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  Ymass_heptane(i,j,k+1) = Ymass_heptane(i,j,k)
                  rhoY_heptane (i,j,k+1) = rhoY_heptane (i,j,k)
                  Xmole_heptane(i,j,k+1) = Xmole_heptane(i,j,k)
               end do
            end do
         end if
      end if

      call fs%cfg%sync(Ymass_heptane)
      call fs%cfg%sync(rhoY_heptane)
      call fs%cfg%sync(Xmole_heptane)
   end subroutine apply_species_bcond

   ! transport of methane mass fraction with an explicit finite-volume scheme
   subroutine advance_heptane_species(dt)
      implicit none
      real(WP), intent(in) :: dt
      integer :: i,j,k
      real(WP) :: mass_flux,Y_face,rhoD_face,gradY
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ,rhoY_new

      call apply_inlet_state()
      call apply_species_bcond()

      allocate(FX(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(FY(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(FZ(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      allocate(rhoY_new(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
      FX = 0.0_WP; FY = 0.0_WP; FZ = 0.0_WP
      rhoY_new = rhoY_heptane
      ! Compute diffusive fluxes with upwind advection and central diffusion.
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_+1
               mass_flux = fs%rho_U(i,j,k)*fs%U(i,j,k)
               if (mass_flux>=0.0_WP) then
                  Y_face = Ymass_heptane(i-1,j,k)
               else
                  Y_face = Ymass_heptane(i,j,k)
               end if
               rhoD_face = max(D_species,0.0_WP)*fs%rho_U(i,j,k)
               gradY = (Ymass_heptane(i,j,k)-Ymass_heptane(i-1,j,k))*fs%cfg%dxmi(i)
               FX(i,j,k) = mass_flux*Y_face - rhoD_face*gradY
            end do
         end do
      end do

      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_+1
            do i=fs%cfg%imin_,fs%cfg%imax_
               mass_flux = fs%rho_V(i,j,k)*fs%V(i,j,k)
               if (mass_flux>=0.0_WP) then
                  Y_face = Ymass_heptane(i,j-1,k)
               else
                  Y_face = Ymass_heptane(i,j,k)
               end if
               rhoD_face = max(D_species,0.0_WP)*fs%rho_V(i,j,k)
               gradY = (Ymass_heptane(i,j,k)-Ymass_heptane(i,j-1,k))*fs%cfg%dymi(j)
               FY(i,j,k) = mass_flux*Y_face - rhoD_face*gradY
            end do
         end do
      end do

      do k=fs%cfg%kmin_,fs%cfg%kmax_+1
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               mass_flux = fs%rho_W(i,j,k)*fs%W(i,j,k)
               if (mass_flux>=0.0_WP) then
                  Y_face = Ymass_heptane(i,j,k-1)
               else
                  Y_face = Ymass_heptane(i,j,k)
               end if
               rhoD_face = max(D_species,0.0_WP)*fs%rho_W(i,j,k)
               gradY = (Ymass_heptane(i,j,k)-Ymass_heptane(i,j,k-1))*fs%cfg%dzmi(k)
               FZ(i,j,k) = mass_flux*Y_face - rhoD_face*gradY
            end do
         end do
      end do

      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               rhoY_new(i,j,k) = rhoY_heptane_old(i,j,k) - dt*( &
                  (FX(i+1,j,k)-FX(i,j,k))*fs%cfg%dxi(i) + &
                  (FY(i,j+1,k)-FY(i,j,k))*fs%cfg%dyi(j) + &
                  (FZ(i,j,k+1)-FZ(i,j,k))*fs%cfg%dzi(k) )
               rhoY_new(i,j,k) = max(0.0_WP,min(fs%RHO(i,j,k),rhoY_new(i,j,k)))
               rhoY_heptane(i,j,k) = rhoY_new(i,j,k)
               Ymass_heptane(i,j,k) = rhoY_heptane(i,j,k)/max(fs%RHO(i,j,k),1.0e-20_WP)
            end do
         end do
      end do

      deallocate(FX,FY,FZ,rhoY_new)
      call update_heptane_mole_fraction()
      call apply_inlet_state()
      call apply_species_bcond()
   end subroutine advance_heptane_species


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read, param_getsize
      implicit none
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
         use vfs_class, only: elvira, remap
         integer :: i,j,k
         call param_read('Jet diameter', Djet)
         call param_read('Jet location', xjet)

         call vf%initialize(cfg=cfg,reconstruction_method=elvira,transport_method=remap,name='VOF')
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Start with pure N2 outside the injected n-heptane jet.
                  vf%VF(i,j,k)      = 0.0_WP
                  vf%Lbary(:,i,j,k) = [vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k) = [vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do
         call vf%update_band()
         call vf%build_interface()
         call vf%set_full_bcond()
         call vf%polygonalize_interface()
         call vf%distance_from_polygon()
         call vf%subcell_vol()
         call vf%get_curvature()
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      create_and_initialize_flow_solver: block
         use mast_class,      only: clipped_neumann, dirichlet,mech_egy_mech_hhz,bc_scope,bcond
         use hypre_str_class, only: hypre_str, pcg_pfmg
         integer :: i,j,k,n
         real(WP) :: gamm_n2,visc_n2,R_n2,Pref_n2,q_n2,b_n2
         real(WP) :: gamm_heptane,visc_heptane,R_heptane,Pref_heptane,q_heptane,b_heptane
         type(bcond), pointer :: mybc
         real(WP) :: P_stat, T_stat, T_cof_val, T_jet_val, U_cof_calc, Ujet_val
         real(WP) :: Pr_n2, Pr_heptane, cp_n2, cp_heptane
         real(WP) :: ycell, uvel, tcell, rho_gas, rho_liq
         matmod = matm(cfg=cfg, name='gas-gas models')

         call param_read('N2 gamma',    gamm_n2)
         call param_read('N2 Pref',     Pref_n2)
         call param_read('N2 q',        q_n2)
         call param_read('N2 b',        b_n2)
         call param_read('N2 viscosity',visc_n2)
         call param_read('N2 constant', R_n2)
         call param_read('N2 Prandtl',  Pr_n2, default=0.72_WP)

         call param_read('Heptane gamma',    gamm_heptane)
         call param_read('Heptane Pref',     Pref_heptane)
         call param_read('Heptane q',        q_heptane)
         call param_read('Heptane b',        b_heptane)
         call param_read('Heptane viscosity',visc_heptane)
         call param_read('Heptane constant', R_heptane)
         call param_read('Heptane Prandtl',  Pr_heptane, default=0.72_WP)
         call param_read('U jet',            Ujet_val)
         call param_read('Jet diameter',     Djet)
         call param_read('Initial jet offset', x0_init, default=Djet)
         call param_read('Inlet velocity radius', inlet_velocity_radius, default=0.5_WP*cfg%yL)
         call param_read('Jet location',     xjet)

         call param_read('Static pressure',  P_stat)
         call param_read('Static temperature',T_stat)
         call param_read('Coflow temperature',T_cof_val,default=T_stat)
         call param_read('Jet temperature',  T_jet_val,default=T_stat)
         call param_read('Coflow velocity',  U_cof_calc)
         call param_read('Species diffusivity', D_species, default=1.0e-5_WP)

         P_inf = P_stat
         T_cof_inlet = T_cof_val
         T_jet_inlet = T_jet_val
         T_inf = T_cof_inlet
         Ucof  = U_cof_calc
         Ujet  = Ujet_val
         gamm_n2_ref  = gamm_n2
         Pref_n2_ref  = Pref_n2
         q_n2_ref     = q_n2
         b_n2_ref     = b_n2
         visc_n2_ref  = visc_n2
         R_n2_ref     = R_n2
         gamm_heptane_ref = gamm_heptane
         Pref_heptane_ref = Pref_heptane
         q_heptane_ref    = q_heptane
         b_heptane_ref    = b_heptane
         visc_heptane_ref = visc_heptane
         R_heptane_ref    = R_heptane
         Pr_n2_ref    = Pr_n2
         Pr_heptane_ref   = Pr_heptane
         cv_n2_ref    = R_n2_ref /max(gamm_n2_ref -1.0_WP,1.0e-20_WP)
         cv_heptane_ref   = R_heptane_ref/max(gamm_heptane_ref-1.0_WP,1.0e-20_WP)
         cp_n2        = cv_n2_ref + R_n2_ref
         cp_heptane       = cv_heptane_ref + R_heptane_ref
         k_n2_ref     = visc_n2_ref * cp_n2 / max(Pr_n2_ref,1.0e-20_WP)
         k_heptane_ref    = visc_heptane_ref* cp_heptane/max(Pr_heptane_ref,1.0e-20_WP)

         call matmod%register_NobleAbelstiffenedgas('gas',    gamm_n2,  Pref_n2,  q_n2,  b_n2)
         call matmod%register_NobleAbelstiffenedgas('liquid', gamm_heptane, Pref_heptane, q_heptane, b_heptane)

         fs=mast(cfg=cfg,name='Single-phase All-Mach',vf=vf)
         call matmod%register_thermoflow_variables('gas',    fs%Grho, fs%Ui, fs%Vi, fs%Wi, fs%GrhoE, fs%GP)
         call matmod%register_thermoflow_variables('liquid', fs%Lrho, fs%Ui, fs%Vi, fs%Wi, fs%LrhoE, fs%LP)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_n2, viscconst_liquid=visc_heptane)
         fs%sigma = 0.0_WP

         ps = hypre_str(cfg=cfg, name='Pressure', method=pcg_pfmg, nst=7)
         ps%maxlevel = 10
         call param_read('Pressure iteration', ps%maxit)
         call param_read('Pressure tolerance', ps%rcvg)

         vs = hypre_str(cfg=cfg, name='Velocity', method=pcg_pfmg, nst=7)
         call param_read('Implicit iteration', vs%maxit)
         call param_read('Implicit tolerance', vs%rcvg)
         call fs%setup(pressure_solver=ps, implicit_solver=vs)
         fs%use_external_temperature = .true.
         fs%use_external_transport   = .true.
         fs%use_external_rhoCv       = .true.
         fs%use_external_EOS_v       = .true.

         rho_gas = 1.0_WP / ( (R_n2  * T_cof_inlet)/(P_stat + Pref_n2)  + b_n2  )
         rho_liq = 1.0_WP / ( (R_heptane * T_jet_inlet)/(P_stat + Pref_heptane) + b_heptane )

         ! Initialize the domain as uniform coflow; the jet profile is only
         ! imposed at the inlet below and in apply_inlet_state().
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  uvel = Ucof
                  fs%Ui(i,j,k) = uvel
                  fs%Vi(i,j,k) = 0.0_WP
                  fs%Wi(i,j,k) = 0.0_WP

                  fs%Grho(i,j,k)  = rho_gas
                  fs%GP(i,j,k)    = P_stat
                  fs%GrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Grho(i,j,k), uvel, 0.0_WP, 0.0_WP, 'gas')

                  fs%Lrho(i,j,k)  = rho_liq
                  fs%LP(i,j,k)    = P_stat
                  fs%LrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Lrho(i,j,k), uvel, 0.0_WP, 0.0_WP, 'liquid')
               end do
            end do
         end do
         fs%P  = P_stat
         fs%PA = P_stat
         ! Boundary conditions
         call fs%add_bcond(name='jet_inlet', type=dirichlet,        locator=jet_inlet_locator, celldir='xm')
         call fs%add_bcond(name='coflow',    type=dirichlet,        locator=wall_locator,      celldir='xm')
         call fs%add_bcond(name='outflow',   type=clipped_neumann,  locator=xp_locator,        celldir='xp')

         call fs%interp_vel_basic(vf, fs%Ui, fs%Vi, fs%Wi, fs%U, fs%V, fs%W)
         call fs%get_bcond('coflow', mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            call get_inlet_profile(j,k,ycell,uvel,tcell)
            fs%U(i:i+1,j,k) = uvel
         end do
         call fs%get_bcond('jet_inlet', mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            call get_inlet_profile(j,k,ycell,uvel,tcell)
            fs%U(i:i+1,j,k) = uvel
         end do
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt, bc_scope)

         fs%RHO   = (1.0_WP - vf%VF)*fs%Grho + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO * fs%Ui
         fs%rhoVi = fs%RHO * fs%Vi
         fs%rhoWi = fs%RHO * fs%Wi

         relax_model = mech_egy_mech_hhz
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  uvel = fs%Ui(i,j,k)
                  fs%P(i,j,k)     = P_stat
                  fs%PA(i,j,k)    = P_stat
                  fs%GP(i,j,k)    = P_stat
                  fs%LP(i,j,k)    = P_stat
                  fs%GrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Grho(i,j,k), uvel, 0.0_WP, 0.0_WP, 'gas')
                  fs%LrhoE(i,j,k) = matmod%EOS_energy(P_stat, fs%Lrho(i,j,k), uvel, 0.0_WP, 0.0_WP, 'liquid')
               end do
            end do
         end do
         call fs%init_phase_bulkmod(vf,matmod)
         fs%psolv%sol = 0.0_WP

      end block create_and_initialize_flow_solver

      create_ensight: block
         ! Allocate mass/mole fraction arrays
         allocate(Ymass_heptane(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ymass_heptane_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rhoY_heptane(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rhoY_heptane_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Xmole_heptane(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(R_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gamma_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Pref_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(q_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(b_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(mu_mix(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(nu_mix_offset(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         Ymass_heptane = 0.0_WP; Ymass_heptane_old = 0.0_WP
         rhoY_heptane = 0.0_WP; rhoY_heptane_old = 0.0_WP
         Xmole_heptane = 0.0_WP
         R_mix = R_n2_ref
         gamma_mix = gamm_n2_ref
         Pref_mix = Pref_n2_ref
         q_mix = q_n2_ref
         b_mix = b_n2_ref
         mu_mix = visc_n2_ref
         nu_mix_offset = 0.0_WP

         ens_out = ensight(cfg=cfg, name='GasJet')
         ens_evt = event(time=time, name='Ensight output')
         call param_read('Ensight output period', ens_evt%tper)

         call ens_out%add_vector('velocity',     fs%Ui, fs%Vi, fs%Wi)
         call ens_out%add_scalar('Pressure',     fs%P)
         call ens_out%add_scalar('PA',           fs%PA)
         call ens_out%add_scalar('Density',      fs%RHO)
         call ens_out%add_scalar('VOF',          vf%VF)
         call ens_out%add_scalar('Temperature',  fs%Tmptr)
         call ens_out%add_scalar('Grho',         fs%Grho)
         call ens_out%add_scalar('Lrho',         fs%Lrho)
         call ens_out%add_scalar('MassFrac_heptane', Ymass_heptane)
         call ens_out%add_scalar('MoleFrac_heptane', Xmole_heptane)
         call ens_out%add_scalar('R_mix',        R_mix)
         call ens_out%add_scalar('gamma_mix',    gamma_mix)
         call ens_out%add_scalar('mu_mix',       mu_mix)
         call ens_out%add_scalar('kappa_mix',    fs%therm_cond)
         call ens_out%add_scalar('rhoCv_mix',    fs%rhoCv)

         call update_heptane_mole_fraction()
         Ymass_heptane_old = Ymass_heptane
         rhoY_heptane_old  = rhoY_heptane
         call apply_inlet_state()
         call apply_species_bcond()
         call update_mixture_from_species(update_pressure=.true.)
         call fs%interp_pressure_density(vf)
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      create_monitor: block
         call fs%get_cfl(time%dt, time%cfl)
         call fs%get_max()
         call vf%get_max()

         mfile = monitor(fs%cfg%amRoot, 'simulation')
         call mfile%add_column(time%n,         'Step')
         call mfile%add_column(time%t,         'Time')
         call mfile%add_column(time%dt,        'dt')
         call mfile%add_column(time%cfl,       'CFL_max')
         call mfile%add_column(fs%Umax,        'U_max')
         call mfile%add_column(fs%Pmax,        'P_max')
         call mfile%add_column(vf%VFmax,       'VOF maximum')
         call mfile%add_column(vf%VFmin,       'VOF minimum')
         call mfile%add_column(vf%VFint,       'VOF integral')
         call mfile%add_column(fs%psolv%it,    'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,  'Pressure error')
         call mfile%write()

         cflfile = monitor(fs%cfg%amRoot, 'cfl')
         call cflfile%add_column(time%n,      'Step')
         call cflfile%add_column(fs%CFLc_x,   'CFL_conv_x')
         call cflfile%add_column(fs%CFLa_x,   'CFL_acous_x')
         call cflfile%write()

         cvgfile = monitor(fs%cfg%amRoot, 'cvg')
         call cvgfile%add_column(time%n,          'Step')
         call cvgfile%add_column(time%it,         'SubIter')
         call cvgfile%add_column(fs%psolv%it,     'P_iter')
         call cvgfile%add_column(fs%psolv%rerr,   'P_err')
         call cvgfile%write()
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      do while (.not.time%done())

         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         call apply_inlet_state()
         call apply_species_bcond()
         call update_mixture_from_species(update_pressure=.true.)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold  = fs%RHO
         fs%Grhoold = fs%Grho;  fs%Lrhoold = fs%Lrho
         fs%GrhoEold= fs%GrhoE; fs%LrhoEold= fs%LrhoE
         fs%GPold   = fs%GP;    fs%LPold   = fs%LP
         Ymass_heptane_old = Ymass_heptane
         rhoY_heptane_old = rhoY_heptane

         call vf%copy_interface_to_old()
         call fs%flow_reconstruct(vf)

         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         call fs%flag_sl(time%dt,vf)

         do while (time%it.le.time%itmax)

            call fs%advection_step(time%dt,vf,matmod)
            call update_mixture_from_species(update_pressure=.false.)
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)
            call update_mixture_from_species(update_pressure=.false.)
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            call fs%psolv%setup()
            fs%psolv%sol = fs%PA - fs%P
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            fs%P = fs%P + fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)
            call update_mixture_from_species(update_pressure=.false.)
            call cvgfile%write()
            time%it = time%it + 1

         end do

         call fs%interp_pressure_density(vf)
         call advance_heptane_species(time%dt)
         call update_mixture_from_species(update_pressure=.true.)
         call fs%interp_pressure_density(vf)
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

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
      if (allocated(Ymass_heptane))     deallocate(Ymass_heptane)
      if (allocated(Ymass_heptane_old)) deallocate(Ymass_heptane_old)
      if (allocated(rhoY_heptane))      deallocate(rhoY_heptane)
      if (allocated(rhoY_heptane_old))  deallocate(rhoY_heptane_old)
      if (allocated(Xmole_heptane))     deallocate(Xmole_heptane)
      if (allocated(R_mix))         deallocate(R_mix)
      if (allocated(gamma_mix))     deallocate(gamma_mix)
      if (allocated(Pref_mix))      deallocate(Pref_mix)
      if (allocated(q_mix))         deallocate(q_mix)
      if (allocated(b_mix))         deallocate(b_mix)
      if (allocated(mu_mix))        deallocate(mu_mix)
      if (allocated(nu_mix_offset)) deallocate(nu_mix_offset)
   end subroutine simulation_final


end module simulation

