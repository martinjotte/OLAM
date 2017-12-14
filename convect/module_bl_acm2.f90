module module_bl_acm2
  
  use consts_coms, only: alvlocp,  & ! L_v / C_p
                         alvlor,   & ! L_v / R_dry
                         eps_vap,  & ! for converting vap press and spec humid
                         eps_virt, & ! used in virtual temperature equation
                         grav,     & ! gravitational acceleration
                         grav2,    & ! 2 * grav
                         vonk,     & ! von Karman's constant
                         cpi         ! 1 / specific heat

  use mem_grid,    only: dzt,      & ! Layer thickness at T (half-layer)
                         dzit,     & ! Inverse Layer thickness at t
                         dzm,      & ! Layer thickness at W (half-layer)
                         dzim,     & ! Inverse Layer thickness at W
                         mza,      & ! Number of model levels
                         zm,       & ! Full layer height (W level)
                         zt          ! Half-layer height (T level)

  implicit none

  real, parameter :: ric    = 0.25  ! critical richardson number
  real, parameter :: mvonk  =      - vonk
  real, parameter :: vonk72 = 0.72 * vonk

! Don't re-export symbols from other modules

  private :: alvlocp, alvlor, eps_vap, eps_virt, grav, grav2, vonk
  private :: dzt, dzit, dzm, dzim, mza, zm, zt

contains

!=======================================================================

  subroutine acm2_scalars( iw, moli, pblh, kpblh, zkh )

    use mem_grid,   only: arw, volti, lpw, lsw, nsw_max
    use mem_turb,   only: frac_sfc, sxfer_tk
    use mem_basic,  only: rho, thil
    use mem_ijtabs, only: itab_w
    use misc_coms,  only: dtlm, io6
    use tridiag,    only: tridv, acm_matrix
    use var_tables, only: num_scalar, scalar_tab
    use mem_tend,   only: thilt
    use mem_para

    implicit none

    integer, intent(in)    :: iw, kpblh
    real,    intent(in)    :: pblh
    real,    intent(in)    :: moli
    real,    intent(in)    :: zkh(:)

    logical :: cnvct
    real :: massflx(mza)
    real :: seddy(mza)
    real :: akodz(mza)

    real :: low(mza)
    real :: dia(mza)
    real :: upp(mza)
    real :: dtom(mza)
    real :: frac_sum(nsw_max)
    real :: aa(mza,nsw_max)
    real :: cbot(mza)

    real :: aflux(mza,num_scalar+1)
    real :: rhs  (mza,num_scalar+1)
    real :: soln (mza,num_scalar+1)

    integer :: k, ks, n, ns, ksm, ksp, kk, ksmax
    real :: hovl, mbar, fnl, dens
    real :: dtl, dtli
    
    integer :: kbot, ktop, nsfc, nlev

    kbot = lpw(iw)
    ktop = mza
    nsfc = lsw(iw)
    nlev = ktop - kbot + 1

    dtl  = dtlm(itab_w(iw)%mrlw)
    dtli = 1. / dtl

    hovl = pblh * moli 

    seddy  (:) = zkh(:)
    massflx(:) = 0.0
    
    cnvct = ((hovl < -0.1) .and. (kpblh - kbot > 2))

    do k = 1, nsfc
       frac_sum(k) = sum(frac_sfc(1:k,iw))
    enddo

    ! IF CONVECTIVE, COMPUTE NON-LOCAL CONVECTIVE MASS FLUX
       
    if (cnvct) then

       dens = 0.5 * (rho(kbot,iw) + rho(kbot+1,iw))
       fnl  = vonk72 / ( vonk72 + (mvonk / hovl) ** 0.33333333)
       mbar = zkh(kbot) / (pblh - dzt(kbot)) * dzit(kbot) * fnl / dens

       do k = kbot, kpblh-1
          seddy(k) = zkh(k) * ( 1.0 - fnl)
       enddo

       do k = kbot, min(kbot + nsfc - 1, kpblh-1)
          dens = 0.5 * (rho(k,iw) + rho(k+1,iw))
          massflx(k) = arw(k,iw) * mbar * (pblh - (zm(k) - zm(kbot-1))) * dens
       enddo

       do k = kbot+nsfc, kpblh-1
          dens = 0.5 * (rho(k,iw) + rho(k+1,iw))
          massflx(k) = arw(kbot+nsfc-1,iw) * mbar * (pblh - (zm(k) - zm(kbot-1))) * dens
       enddo

    endif

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - SCALARS

    do k = kbot, ktop - 1
       akodz(k) = arw(k,iw) * seddy(k) * dzim(k)
    enddo

    akodz(kbot-1) = 0.0
    akodz(ktop  ) = 0.0

    do k = kbot, ktop
       ks = k - kbot + 1

       dtom(k) = dtl * volti(k,iw) / rho(k,iw)

       low(ks) = - dtom(k) * akodz(k-1)
       upp(ks) = - dtom(k) * akodz(k  )
       dia(ks) = 1.0 - low(ks) - upp(ks)
    enddo

    ! Scalar variables with long-timestep forcing included

    do n = 1, num_scalar
       do k = kbot, ktop
          ks = k - kbot + 1
          rhs(ks,n) = scalar_tab(n)%var_p(k,iw)  &
                    + dtl * scalar_tab(n)%var_t(k,iw) / rho(k,iw)
       enddo
    enddo

    n = num_scalar + 1
    do k = kbot, ktop
       ks = k - kbot + 1
       rhs(ks,n) = thil(k,iw) + dtl * thilt(k,iw) / rho(k,iw)
    enddo

    ! IF CONVECTIVE, INCLUDE NONLOCAL TERMS
       
    if (cnvct) then

       do k = kbot+1, kpblh
          ks = k - kbot + 1
          dia(ks)   = dia(ks)   + dtom(k)   * massflx(k-1)
          upp(ks-1) = upp(ks-1) - dtom(k-1) * massflx(k-1)
       enddo

       aa(:,:) = 0.0
       ksmax = min(nsfc, kpblh - kbot)

       do kk = 1, ksmax
          do k = kk + kbot, kpblh
             
             ks  = k - kbot + 1
             ksp = min(ks,  nsfc)
             ksm = min(ks-1,nsfc)

             aa(ks,kk) = dtom(k) * frac_sfc(kk,iw) * &
                       ( massflx(k) / frac_sum(ksp) - massflx(k-1) / frac_sum(ksm) )
          enddo
       enddo

       do k = kbot, min(kbot + nsfc - 1, kpblh)
          ks = k - kbot + 1
          dia(ks) = dia(ks) + dtom(k) * massflx(k) * frac_sfc(ks,iw) / frac_sum(ks)
          low(ks+1) = low(ks+1) + aa(ks+1,ks)
       enddo

       call acm_matrix(aa, dia, low, rhs, upp, soln, nlev, num_scalar+1, ksmax)

    else

       call tridv(low, dia, upp, rhs, soln, 1, nlev, mza, num_scalar+1)

    endif

    ! Now, soln contains future(t+1) values

    ! Compute internal vertical turbulent fluxes due to Kh

    do n = 1, num_scalar+1
       do k = kbot, ktop-1
          ks = k - kbot + 1
          aflux(k,n) = akodz(k) * (soln(ks,n) - soln(ks+1,n))
       enddo
    enddo

    ! Include internal nonlocal vertical turbulent fluxes

    if (cnvct) then

       do n = 1, num_scalar+1

          cbot(1) = soln(1,n)
          do k = 2, nsfc
             cbot(k) = sum( soln(1:k,n)*frac_sfc(1:k,iw) ) / frac_sum(k)
          enddo
          cbot(nsfc+1:) = cbot(nsfc)

          do k = kbot, kpblh
             ks = k - kbot + 1
             aflux(k,n) = aflux(k,n) + massflx(k) * (cbot(ks) - soln(ks+1,n))
          enddo
          
       enddo
    endif

    ! Set bottom and top internal fluxes to zero
    
    aflux(kbot-1,:) = 0.
    aflux(ktop,  :) = 0.

    ! Compute tendencies due to internal turbulent fluxes

    do n = 1, num_scalar
       do k = kbot, ktop
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + volti(k,iw) * (aflux(k-1,n) - aflux(k,n))
       enddo
    enddo

    n = num_scalar + 1
    do k = kbot, ktop
       thilt(k,iw) = thilt(k,iw) + volti(k,iw) * (aflux(k-1,n) - aflux(k,n))
    enddo

  end subroutine acm2_scalars

!=======================================================================

  subroutine acm2_momentum( iw, zkm )

    use mem_grid,   only: arw, volti, lpw, lsw, nsw_max
    use mem_turb,   only: vkm_sfc
    use mem_basic,  only: vxe, vye, vze, rho
    use mem_tend,   only: vmxet, vmyet, vmzet
    use mem_ijtabs, only: itab_w
    use misc_coms,  only: dtlm
    use tridiag,    only: tridv

    implicit none

    integer, intent(in)    :: iw
    real,    intent(in)    :: zkm(:)

    real :: aflux(mza,3)
    real :: rhs  (mza,3)
    real :: soln (mza,3)

    real :: akodz(mza)

    real :: low(mza)
    real :: dia(mza)
    real :: upp(mza)
    real :: dtom(mza)
    real :: fact(nsw_max)

    integer :: k, ks
    integer :: kbot, ktop, nsfc

    kbot = lpw(iw)
    ktop = mza
    nsfc = lsw(iw)

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - MOMENTUM

    akodz(kbot-1) = 0.0

    do k = kbot, ktop-1
       akodz(k) = arw(k,iw) * zkm(k) * dzim(k)
    enddo

    akodz(ktop) = 0.0

    do k = kbot, ktop
       dtom(k)  = dtlm(itab_w(iw)%mrlw) * volti(k,iw) / rho(k,iw)
       low(k)   = -dtom(k) * akodz(k-1)
       upp(k)   = -dtom(k) * akodz(k  )
       dia(k)   = 1.0 - low(k) - upp(k)
       rhs(k,1) = vxe(k,iw) 
       rhs(k,2) = vye(k,iw)
       rhs(k,3) = vze(k,iw)
    enddo

    ! Surface drag
    
    do k = kbot, kbot + nsfc - 1
       ks = k - kbot + 1
       fact(ks) = 2.0 * vkm_sfc(ks,iw) * dzim(k-1) * (arw(k,iw) - arw(k-1,iw))
       dia(k)   = dia(k) + dtom(k) * fact(ks)
    enddo

    ! tridv - low, diag, upper, rhs, soln
               
    call tridv( low, dia, upp, rhs, soln, kbot, ktop, mza, 3)

    ! Now, soln contains future(t+1) values

    ! Compute internal vertical turbulent fluxes

    do k = kbot, mza-1
       aflux(k,1) = akodz(k) * (soln(k,1) - soln(k+1,1))
       aflux(k,2) = akodz(k) * (soln(k,2) - soln(k+1,2))
       aflux(k,3) = akodz(k) * (soln(k,3) - soln(k+1,3))
    enddo

    ! Set bottom and top internal fluxes to zero
    
    aflux(kbot-1,1:3) = 0.
    aflux(mza,   1:3) = 0.

    ! Compute tendencies due to internal turbulent fluxes

    do k = kbot, mza
       vmxet(k,iw) = vmxet(k,iw) + volti(k,iw) * (aflux(k-1,1) - aflux(k,1))
       vmyet(k,iw) = vmyet(k,iw) + volti(k,iw) * (aflux(k-1,2) - aflux(k,2))
       vmzet(k,iw) = vmzet(k,iw) + volti(k,iw) * (aflux(k-1,3) - aflux(k,3))
    enddo

    ! Include tendencies due to surface drag

    do k = kbot, kbot + nsfc - 1
       ks = k - kbot + 1
       vmxet(k,iw) = vmxet(k,iw) - volti(k,iw) * fact(ks) * soln(k,1)
       vmyet(k,iw) = vmyet(k,iw) - volti(k,iw) * fact(ks) * soln(k,2)
       vmzet(k,iw) = vmzet(k,iw) - volti(k,iw) * fact(ks) * soln(k,3)
    enddo

  end subroutine acm2_momentum
    
!=======================================================================

  subroutine acm2_eddyx(iw, moli, ustar, wstar, pblh, kpblh, kbot, ktop, zkm, &
       zkh, vx, vy, vz, qw, qv, ql, thil, theta, thetav, tair, delr, rho)

    ! THREE METHODS FOR COMPUTING KZ:
    !   1. Boundary scaling similar to Holtslag and Boville (1993)
    !   2. Local Kz computed as function of local Richardson # and
    !      vertical wind shear, similar to LIU & CARROLL (1996)
    !   3. Cloud-top radiational cooling scaling, Lock et al. (2000)

    use misc_coms,  only: io6
    use mem_cuparm, only: iactcu, kcubot
    implicit none

    ! INPUT VARIABLES

    real,    intent(in) :: moli      ! inverse Monin-Obukhov length
    real,    intent(in) :: ustar     ! u* surface-layer velocity scale
    real,    intent(in) :: wstar     ! w* convective PBL velocity scale
    real,    intent(in) :: pblh      ! PBL height
    integer, intent(in) :: kpblh     ! level corresponding to PBL height
    integer, intent(in) :: kbot      ! lowest model level
    integer, intent(in) :: ktop      ! highest model level
    real,    intent(in) :: vx    (:) ! Earth-cartesian x-wind profile
    real,    intent(in) :: vy    (:) ! Earth-cartesian y-wind profile
    real,    intent(in) :: vz    (:) ! Earth-cartesian z-wind profile
    real,    intent(in) :: qw    (:) ! total water specific humidity profile
    real,    intent(in) :: qv    (:) ! water vapor specific humidity profile
    real,    intent(in) :: ql    (:) ! condensate specific humidity profile
    real,    intent(in) :: thil  (:) ! ice-liquid potential temperature profile
    real,    intent(in) :: theta (:) ! potential temperature profile
    real,    intent(in) :: thetav(:) ! virtual potential temperature profile
    real,    intent(in) :: tair  (:) ! air temperature profile
    real,    intent(in) :: delr      ! cloud top cooling at PBL top (K m/s)
    real(8), intent(in) :: rho   (:) ! air density
    integer, intent(in) :: iw

    ! OUTPUT VARIABLES

    real, intent(out) :: zkm(:) ! eddy diffusivity for momentum
    real, intent(out) :: zkh(:) ! eddy diffusivity for scalars

    ! LOCAL VARIABLES

    integer :: k, kpblm1, kpblp1, kb, ka
    real    :: zagl, zoh, zsol, dz, dzi
    real    :: buoy, ss, kprof, rlam
    real    :: ri, zk, sql, fh, fm, pr
    real    :: phimi, phihi
    real    :: we, dthv, deltar, wcloud
    logical :: inlay

    real :: edyzh(mza) ! K-profile eddy diffusivity within PBL
    real :: edyzm(mza) ! K-profile eddy diffusivity within PBL 
    real :: edyrh(mza) ! Richardson-number eddy diffusivity 
    real :: edyrm(mza) ! Richardson-number eddy diffusivity 
    real :: kzo  (mza) ! background eddy diffusivity

    ! PARAMETERS

    REAL, PARAMETER :: RLAM_stab =  30.0
    REAL, PARAMETER :: RLAM_unst = 150.0
    REAL, PARAMETER :: GAMH      =  16.0 ! Holtslag and Boville (1993)

    REAL, PARAMETER :: EDYZ0  = 0.0   ! New Min Kz
!   REAL, PARAMETER :: EDYZ0  = 0.01  ! New Min Kz

    kzo  (:) = edyz0
    edyzh(:) = 0.0
    edyzm(:) = 0.0
    edyrh(:) = 0.0
    edyrm(:) = 0.0

!! Holtslag's Eddy-Diffusivity Profile Within the convective PBL

    if (moli <= 0.0) then

       ! Vertical loop over W levels
       do k = kbot, kpblh-1
          zagl     = zm(k) - zm(kbot-1)
          zsol     = min(0.1 * pblh, zagl) * moli
          
          phimi    =     (1.0 - gamh * zsol)**0.25
          phihi    = sqrt(1.0 - gamh * zsol)

          kprof    = vonk * ustar * zagl * (1.0 - zagl/pblh)**2
          edyzh(k) = kprof * phihi
          edyzm(k) = kprof * phimi
       enddo

    endif

!! Lock's Cloud-Top Driven Eddy-Diffusivity Profile

    wcloud  = 0.0
    kpblm1  = max(kpblh-1, kbot)

    if ( ql(kpblh) > 1.e-8 .or. ql(kpblm1) > 1.e-8 .or. &
         (iactcu(iw) > 0 .and. kcubot(iw) <= kpblh) ) then
       deltar = max( delr, 0.0 )
    else
       deltar = 0.0
    endif

    if (deltar > 2.e-7) then

       wcloud = (grav / thetav(kpblh) * pblh * deltar)**0.33333333

       do k = kbot, kpblh-1
          zagl    = zm(k) - zm(kbot-1)
          zoh     = zagl / pblh
               
        ! If cloud effects limited to just cloud layer
        ! (see ECMWF physics technote):
        ! zagl    = max(zm(k) - zm(kbot-1) - zbase, 0.0)
        ! zoh     = zagl / (pblh - zbase)

          kprof    = 0.85 * vonk * wcloud * zagl * zoh * sqrt(1.0 - zoh)
          edyzh(k) = edyzh(k) +        kprof
          edyzm(k) = edyzm(k) + 0.75 * kprof
       enddo

    endif

!! Explicit inclusion of entrainment

    if (deltar >= 1.e-9 .or. moli <= 0.0) then

       kpblp1 = min(kpblh+1, ktop)
       kpblm1 = max(kpblh-1, kbot)

       dthv = max(thetav(kpblp1) - thetav(kpblm1), 0.5)
       we   = 0.2 * (wstar**3 + wcloud**3 + 5.0*ustar**3) * thetav(kpblh) &
            / (grav * pblh * dthv)

       edyzh(kpblh)   =        we * dzm(kpblh)
       edyzm(kpblh)   = 0.75 * we * dzm(kpblh)

       edyzh(kpblm1) = max(       we * dzm(kpblm1), edyzh(kpblm1))
       edyzm(kpblm1) = max(0.75 * we * dzm(kpblm1), edyzm(kpblm1))

    endif

!! Louis's Richardson-number dependent eddy diffusivity; used everywhere
!! except for within the unstable PBL

    edyrh(:) = 0.0
    edyrm(:) = 0.0

    if (moli <= 0.0) then
       kb = kpblh
    else
       kb = kbot
    endif

    inlay = .false.
    ka    = kb
      
    ! Vertical loop over W levels
    do k = kb, ktop-1

       if (.not. inlay) ka = k

       dz  = zt(k+1) - zt(ka)
       dzi = 1.0 / dz
       zk  = vonk * dz

       ss = ( (vx(k+1) - vx(ka))**2 &
            + (vy(k+1) - vy(ka))**2 &
            + (vz(k+1) - vz(ka))**2 ) * dzi * dzi
       ss = max(ss, 1.e-6)

       buoy = (thetav(k+1) - thetav(ka)) * dzi

       ri = grav2 / (thetav(k+1) + thetav(ka)) * buoy / ss
       ri = max(ri, -1.0)

       if (ri > ric) then
          inlay = .false.
       else
          inlay = .true.
       endif

       if (ri > 0.0) then

          if (k < kpblh) then
             rlam = max(0.1 * pblh, rlam_stab)
          else
             rlam = rlam_stab
          endif

          pr = min(1.0 + 3.7 * ri, 3.0)
          fh = 1.0 / (1.0 + ri * ( 10.0 + ri * ( 50.0 + 5000.0 * ri * ri ) ) )
          fm = fh / pr

          sql = (zk * rlam / (rlam + zk))**2       

          edyrm(k) = sqrt(ss) * sql * fm
          edyrh(k) = sqrt(ss) * sql * fh

       else

          sql = (zk * rlam_unst / (rlam_unst + zk))**2

          phimi =     (1.0 - 16.0 * ri)**0.25
          phihi = sqrt(1.0 - 16.0 * ri)

          edyrm(k) = sqrt(ss) * sql * phimi * phimi
          edyrh(k) = sqrt(ss) * sql * phimi * phihi

       endif

    enddo

!! Now compute combined eddy diffusivity from surface and cloud K-profiles
!! and Richardson-number eddy diffusivity

    do k = kbot, mza-1
       zkh(k) = max( edyrh(k), edyzh(k), kzo(k) )
       zkh(k) = min( zkh(k), 1000.0 )
       zkh(k) = 0.5 * (rho(k+1) + rho(k)) * zkh(k)

       zkm(k) = max( edyrm(k), edyzm(k), kzo(k) )
       zkm(k) = min( zkm(k), 1000.0 )
       zkm(k) = 0.5 * (rho(k+1) + rho(k)) * zkm(k)
    enddo

    zkh(1:kbot-1) = 0.0
    zkh(mza)      = 0.0

    zkm(1:kbot-1) = 0.0
    zkm(mza)      = 0.0

  end subroutine acm2_eddyx

!=======================================================================

  subroutine acm2_pblhgt( ustar, wstar, wtv0, kbot, ktop, nsfc, &
                          fracsfc, thv, vx, vy, vz, kpblh, pblh )
    
    ! PURPOSE: TO CALCULATE BOUNDARY LAYER HEIGHT WITH THE METHOD
    !          PUBLISHED BY HOLTSLAG (MWR AUGUST 1990).

    implicit none

    ! INPUT VARIABLES

    real,    intent(in) :: ustar  ! u* surface velocity scale [m/s]
    real,    intent(in) :: wstar  ! w* convective PBL velocity scale [m/s]
    real,    intent(in) :: wtv0   ! surface virtual temperature flux [K-m/s]

    integer, intent(in) :: kbot   ! index of lowest layer
    integer, intent(in) :: ktop   ! index of highest layer
    integer, intent(in) :: nsfc   ! number of levels that intersect ground

    real,    intent(in) :: thv(:) ! virtual temperature profile
    real,    intent(in) :: vx (:) ! earth-cartesian x velocity profile
    real,    intent(in) :: vy (:) ! earth-cartesian y velocity profile
    real,    intent(in) :: vz (:) ! earth-cartesian z velocity profile

    real,    intent(in) :: fracsfc(:) ! fraction of the total surface area that 
                                      ! intersects each nsfc layer

    ! OUTPUT VARIABLES

    integer, intent(out) :: kpblh ! Layer corresponding to PBL height
    real,    intent(out) :: pblh  ! Height of PBL

    ! LOCAL VARIABLES

    real    :: thvsfc ! near-ground virtual temperature
    real    :: fint   ! vertical interpolation factor
    real    :: rib    ! bulk Richardson number
    real    :: zmix   ! height of mixed layer
    integer :: kmix   ! vertical index of mixed layer height
    real    :: wscale ! convective velocity scale
    real    :: wssq   ! wind shear term in the Richardson number 

    real :: vxmix, vymix, vzmix, rib_sav
    integer :: k

    if (kbot+1 >= ktop-1) then
       kpblh = kbot
       pblh  = dzt(kbot)
       return
    endif

    ! COMPUTE AN AVERAGE NEAR-SURFACE VIRTUAL POTENTIAL TEMPERATURE

    if (nsfc == 1) then
       thvsfc = thv(kbot)
    else
       thvsfc = sum( thv(kbot:kbot+nsfc-1) * fracsfc(1:nsfc) )
    endif

    ! IF SURFACE LAYER IS UNSTABLE, COMPUTE CONVECTIVE VELOCITY SCALE 
    ! AND THERMAL EXCESS TEMPERATURE:

    if (wtv0 > 0.0) then
       wscale = (ustar**3 + 0.6 * wstar**3) ** 0.33333333
       thvsfc = thvsfc + 8.5 * wtv0 / wscale
    else
       wscale = 0.0
    endif

    ! FIND FIRST LEVEL ABOVE THE MIXED LAYER, IF IT EXISTS

    k = kbot+1
    do k = kbot+1, ktop-1
       if (thv(k) > thvsfc) exit
    enddo
    kmix = min(k,ktop-1)

    ! IF MIXED-LAYER EXISTS; INTERPOLATE HEIGHT AND WIND VELOCITIES
    ! TO THE POINT WHERE THETAV EQUALS ITS SURFACE-LAYER VALUE

    if (kmix > kbot+1) then
       fint  = (thvsfc - thv(kmix-1)) / (thv(kmix) - thv(kmix-1))
       zmix  = fint * (zt(kmix)-zt(kmix-1)) + zt(kmix-1)
       vxmix = fint * (vx(kmix)-vx(kmix-1)) + vx(kmix-1)
       vymix = fint * (vy(kmix)-vy(kmix-1)) + vy(kmix-1)
       vzmix = fint * (vz(kmix)-vz(kmix-1)) + vz(kmix-1)
    else
       zmix  = zt(kbot)
       vxmix = vx(kbot)
       vymix = vy(kbot)
       vzmix = vz(kbot)
    endif

    ! STARTING FROM THE TOP OF THE MIXED LAYER (OR SURFACE IF NO MIXED LAYER
    ! EXISTS), SEARCH UPWARD UNTIL THE BULK RICHARDSON NUMBER EQUALS ITS
    ! CRITICAL VALUE. THIS WILL BE THE PBL HEIGHT

    rib     = 0.0
    rib_sav = 0.0

    do k = kmix, ktop-1
       rib_sav = rib
       
       wssq = wscale * wscale     &
            + (vx(k) - vxmix)**2  &
            + (vy(k) - vymix)**2  &
            + (vz(k) - vzmix)**2
            
       wssq = max(wssq, 0.1)

       rib = grav2 * (thv(k) - thvsfc) / (thv(k) + thvsfc) &
            * (zt(k) - zmix) / wssq

       if (rib >= ric) exit
    enddo

    kpblh = min(k,ktop-1)

    if (kpblh == ktop-1) then

       pblh = max( zm(kpblh) - zm(kbot-1), dzt(kbot) )

    else

       ! INTERPOLATE BETWEEN LEVELS TO DETERMINE THE PBL HEIGHT

       fint = (ric - rib_sav) / (rib - rib_sav)

       if (fint > 0.5) then
          fint  = fint - 0.5
       else
          kpblh = kpblh - 1
          fint  = fint  + 0.5
       endif

       pblh = fint * (zm(kpblh) - zm(kpblh-1)) + zm(kpblh-1) - zm(kbot-1)

    endif

  end subroutine acm2_pblhgt

END MODULE module_bl_acm2
