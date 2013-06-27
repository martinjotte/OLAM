module module_bl_acm2
  
  use consts_coms, only: alvlocp,  & ! L_v / C_p
                         alvlor,   & ! L_v / R_dry
                         eps_vap,  & ! for converting vap press and spec humid
                         eps_virt, & ! used in virtual temperature equation
                         grav,     & ! gravitational acceleration
                         grav2,    & ! 2 * grav
                         vonk        ! von Karman's constant

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

  subroutine acm2( iw, rhot, moli, ustar, pblh, kpblh, thetav, zkh )

    use mem_grid,   only: arw, volt, volti, lpw, lsw, nsw_max, mwa, zfacm
    use mem_turb,   only: vkm_sfc, frac_sfc, sxfer_rk
    use mem_basic,  only: vxe, vye, vze, rho, sh_w
    use mem_tend,   only: vmxet, vmyet, vmzet
    use mem_ijtabs, only: itab_w
    use misc_coms,  only: dtlm, io6
    use tridiag,    only: tridv, acm_matrix
    use var_tables, only: num_scalar, scalar_tab, sxfer_map, num_sxfer, &
                          emis_map, num_emis
    use emis_defn,  only: emlays

    implicit none

    integer, intent(in)    :: iw, kpblh
    real,    intent(in)    :: pblh
    real,    intent(in)    :: ustar
    real,    intent(in)    :: moli
    real,    intent(in)    :: thetav(:)
    real,    intent(in)    :: zkh(:)
    real,    intent(inout) :: rhot(mza,mwa)

    logical :: cnvct
    real :: mbarks(mza)
    real :: seddy(mza)
    real :: akodz(mza)

    real :: low(mza)
    real :: dia(mza)
    real :: upp(mza)
    real :: dtom(mza)
    real :: frac_sum(nsw_max)
    real :: aa(mza,nsw_max)
    real :: cbot(mza)

    real :: aflux(mza,max(3,num_scalar))
    real :: rhs  (mza,max(3,num_scalar))
    real :: soln (mza,max(3,num_scalar))

    real :: fact(nsw_max)

    integer :: k, ks, n, ns, ksm, ksp, kk, ksmax, kmax
    real :: hovl, mbar, fnl
    real :: dtl, dtli
    
    integer :: kbot, ktop, nsfc, nlev

    kbot = lpw(iw)
    ktop = mza-1
    nsfc = lsw(iw)
    nlev = ktop - kbot + 1

    dtl  = dtlm(itab_w(iw)%mrlw)
    dtli = 1. / dtl

    hovl = pblh * moli 

    seddy (:) = zkh(:)
    mbarks(:) = 0.0
    
    cnvct = ((hovl < -0.1) .and. (kpblh - kbot > 2))

    do k = 1, nsfc
       frac_sum(k) = sum(frac_sfc(1:k,iw))
    enddo

    ! IF CONVECTIVE, COMPUTE NON-LOCAL CONVECTIVE MASS FLUX
       
    if (cnvct) then

       fnl  = vonk72 / ( vonk72 + (mvonk / hovl) ** 0.33333333)
       mbar = zkh(kbot) / (pblh - dzt(kbot)) * dzit(kbot) * fnl

       do k = kbot, kpblh-1
          seddy(k) = zkh(k) * ( 1.0 - fnl)
       enddo
       
       do k = kbot, kpblh-1
          mbarks(k) = arw(k,iw) * mbar * (pblh - zm(k) + zm(kbot-1))
       enddo

       do k = kbot + nsfc, kpblh-1
          mbarks(k) = mbarks(k) * zfacm(kbot+nsfc-1)**2 / zfacm(k)**2
       enddo

    endif

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - MOMENTUM

    akodz(kbot-1)     = 0.0
    aflux(kbot-1,1:3) = 0.0

    do k = kbot, ktop-1
       akodz(k)   = arw(k,iw) * zkh(k) * dzim(k)
       aflux(k,1) = akodz(k) * (vxe(k,iw) - vxe(k+1,iw))
       aflux(k,2) = akodz(k) * (vye(k,iw) - vye(k+1,iw))
       aflux(k,3) = akodz(k) * (vze(k,iw) - vze(k+1,iw))
    enddo

    akodz(ktop    ) = 0.0
    aflux(ktop,1:3) = 0.0

    do k = kbot, ktop
       dtom(k)  = dtl * volti(k,iw) / rho(k,iw)
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
       dia(k)   = dia(k)   + dtom(k) * fact(ks)
    enddo

    ! tridv - low, diag, upper, rhs, soln
               
    call tridv( low, dia, upp, rhs, soln, kbot, ktop, mza, 3)

    ! Now, soln contains future(t+1) values

    ! Compute internal vertical turbulent fluxes

    do k = kbot, mza-2
       aflux(k,1) = akodz(k) * (soln(k,1) - soln(k+1,1))
       aflux(k,2) = akodz(k) * (soln(k,2) - soln(k+1,2))
       aflux(k,3) = akodz(k) * (soln(k,3) - soln(k+1,3))
    enddo

    ! Set bottom and top internal fluxes to zero
    
    aflux(kbot-1,1:3) = 0.
    aflux(mza-1, 1:3) = 0.

    ! Compute tendencies due to internal turbulent fluxes

    do k = kbot, mza-1
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

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - SCALARS

    do k = kbot, ktop - 1
       akodz(k)   = arw(k,iw) * seddy(k) * dzim(k)
    enddo

    akodz(kbot-1) = 0.0
    akodz(ktop  ) = 0.0

    do k = kbot, ktop
       ks = k - kbot + 1

       low(ks)   = - dtom(k) * akodz(k-1)
       upp(ks)   = - dtom(k) * akodz(k  )
       dia(ks)   = 1.0 - low(ks) - upp(ks)
    enddo

    do n = 1, num_scalar
       do k = kbot, ktop
          ks = k - kbot + 1
          rhs(ks,n) = scalar_tab(n)%var_p(k,iw)
       enddo
    enddo

    ! Include surface exchange
    
    if (num_sxfer > 0) then
       kmax = min(kbot + nsfc - 1, ktop-1)
       do ns = 1, num_sxfer
          n = sxfer_map(ns)
          do k = kbot, kmax
             ks = k - kbot + 1
             rhs(ks,n) = rhs(ks,n) + scalar_tab(n)%sxfer(ks,iw) * volti(k,iw) / rho(k,iw)
          enddo
       enddo
    endif

    ! Include emissions (emis units are concentration / sec )

    if (num_emis > 0) then
       kmax = min(kbot + nsfc + emlays - 2, ktop-1)
       do ns = 1, num_emis
          n = emis_map(ns)
          do k = kbot, kmax
             ks = k - kbot + 1
             rhs(ks,n) = rhs(ks,n) + scalar_tab(n)%emis(k,iw) * dtl
          enddo
       enddo
    endif

    ! IF CONVECTIVE, INCLUDE NONLOCAL TERMS
       
    if (cnvct) then

       do k = kbot+1, kpblh
          ks = k - kbot + 1
          dia(ks)   = dia(ks)   + dtom(k)   * mbarks(k-1)
          upp(ks-1) = upp(ks-1) - dtom(k-1) * mbarks(k-1)
       enddo

       aa(:,:) = 0.0
       ksmax = min(nsfc, kpblh - kbot)

       do kk = 1, ksmax
          do k = kk + kbot, kpblh
             
             ks  = k - kbot + 1
             ksp = min(ks,  nsfc)
             ksm = min(ks-1,nsfc)

             aa(ks,kk) = dtom(k) * frac_sfc(kk,iw) * &
                       ( mbarks(k) / frac_sum(ksp) - mbarks(k-1) / frac_sum(ksm) )
          enddo
       enddo
       
       do k = kbot, min(kbot + nsfc - 1, kpblh)
          ks = k - kbot + 1
          dia(ks) = dia(ks) + dtom(k) * mbarks(k) * frac_sfc(ks,iw) / frac_sum(ks)
          low(ks+1) = low(ks+1) + aa(ks+1,ks)
       enddo

       call acm_matrix(aa, dia, low, rhs, upp, soln, nlev, num_scalar, ksmax)

    else

       call tridv(low, dia, upp, rhs, soln, 1, nlev, mza, num_scalar)

    endif

    ! Now, soln contains future(t+1) values

    ! Compute internal vertical turbulent fluxes due to Kh

    do n = 1, num_scalar
       do k = kbot, ktop-1
          ks = k - kbot + 1
          aflux(k,n) = akodz(k) * (soln(ks,n) - soln(ks+1,n))
       enddo
    enddo

    ! Include internal nonlocal vertical turbulent fluxes

    if (cnvct) then

       do n = 1, num_scalar

          cbot(1) = soln(1,n)
          do k = 2, nsfc
             cbot(k) = sum( soln(1:k,n)*frac_sfc(1:k,iw) ) / frac_sum(k)
          enddo
          cbot(nsfc+1:) = cbot(nsfc)

          do k = kbot, kpblh
             ks = k - kbot + 1
             aflux(k,n) = aflux(k,n) + mbarks(k) * (cbot(ks) - soln(ks+1,n))
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

    ! Include tendencies due to surface exchange

    if (num_sxfer > 0) then
       kmax = min(kbot + nsfc - 1, ktop-1)
       do ns = 1, num_sxfer
          n = sxfer_map(ns)
          do k = kbot, kmax
             ks = k - kbot + 1
             scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                       + dtli * volti(k,iw) * scalar_tab(n)%sxfer(ks,iw)
          enddo
       enddo
    endif

    ! Include tendencies due to emissions (emis units are concentration / sec )

    if (num_emis > 0) then
       kmax = min(kbot + nsfc + emlays - 2, ktop-1)
       do ns = 1, num_emis
          n = emis_map(ns)
          do k = kbot, kmax
             scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                       + rho(k,iw) * scalar_tab(n)%emis(k,iw)
          enddo
       enddo
    endif

    ! Apply surface vapor xfer [kg_vap] directly to rhot [kg_air / (m^3 s)]

    do k = kbot, kbot + nsfc - 1
       ks = k - kbot + 1
       rhot(k,iw) = rhot(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)
    enddo

  end subroutine acm2
    
!=======================================================================

  subroutine acm2_eddyx(iw, moli, ustar, pblh, kpblh, kbot, ktop, zkh, &
       vx, vy, vz, qw, qv, ql, thil, theta, thetav, tair, fthrd, rho)

    ! THREE METHODS FOR COMPUTING KZ:
    !   1. Boundary scaling similar to Holtslag and Boville (1993)
    !   2. Local Kz computed as function of local Richardson # and
    !      vertical wind shear, similar to LIU & CARROLL (1996)
    !   3. Cloud-top radiational cooling scaling, Lock et al. (2000)

    implicit none
    
    ! INPUT VARIABLES
    
    real,    intent(in) :: moli      ! inverse Monin-Obukhov length
    real,    intent(in) :: ustar     ! u* surface-layer velocity scale
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
    real,    intent(in) :: fthrd (:) ! radiational heating/cooling profile
    real(8), intent(in) :: rho   (:) ! air density
    integer, intent(in) :: iw

    ! OUTPUT VARIABLES

    real, intent(out) :: zkh(:) ! eddy diffusivity for scalars
!   real, intent(out) :: zkm(:) ! eddy diffusivity for momentum

    ! LOCAL VARIABLES

    INTEGER  :: ILX, KL, KLM, K, I
    real :: zagl, deltar, wcloud, zoh, a, b, buoy
    REAL     :: ZOVL, WT, ZSOL, ZFUNC, DZF, SS
    REAL     :: RI, QMEAN, TMEAN, XLV, ALPH, CHI, ZK, SQL, DENSF, KZO
    REAL     :: FH

    real :: phih(mza) ! M-O similarity nondimensional scalar gradient
    real :: edyz(mza) ! Holtslag's K-profile eddy diffusivity
    real :: edyc(mza) ! Lock et al. (2000) cloud-top driven K-profile
    real :: edyr(mza) ! Richardson-number eddy diffusivity

    real :: alpha(mza)
    real :: beta(mza)

    ! PARAMETERS

    REAL, PARAMETER :: RV     = 461.5
    REAL, PARAMETER :: RC     = 0.25
    REAL, PARAMETER :: RLAM   = 80.0
    REAL, PARAMETER :: GAMH   = 16.0 !15.0  !  Holtslag and Boville (1993)
    REAL, PARAMETER :: BETAH  = 5.0   !  Holtslag and Boville (1993)
    
!   REAL, PARAMETER :: FHMIN = 0.01  ! Minimum value of f(Ri)
    REAL, PARAMETER :: FHMIN = 0.0

    REAL, PARAMETER :: EDYZ0  = 0.0   ! New Min Kz
!   REAL, PARAMETER :: EDYZ0  = 0.01  ! New Min Kz
!   REAL, PARAMETER :: EDYZ0  = 0.1

    ! Constants for the "moist" richardson number

    real, parameter :: c1 = 1.0 + eps_virt
    real, parameter :: c2 = eps_vap * alvlor
    real, parameter :: c3 = c2 * alvlocp

    kzo = edyz0

!! Holtslag's Eddy-Diffusivity Profile Within the PBL

    edyz(:) = 0.0

    ! Vertical loop over W levels
    do k = kbot, kpblh-1
            
       zagl = zm(k) - zm(kbot-1)
       zovl = zagl * moli

       if (zovl < 0.0) then

          ! Convective PBL:
          if (zagl < 0.1 * pblh) then
             phih(k) = 1.0 / sqrt(1.0 - gamh * zovl)
          else
             zsol = 0.1 * pblh * moli
             phih(k) = 1.0 / sqrt(1.0 - gamh * zsol)
          endif

       elseif (zovl > 1.0) then

          ! Very stable PBL:
          phih(k) = 1.0 + betah * zovl

       else

          ! Stable PBL:
          phih(k) = betah + zovl
          
       endif

       zfunc = zagl * (1.0 - zagl / pblh)**2
       edyz(k) = vonk * ustar * zfunc / phih(k)

    enddo

!! Lock's Cloud-Top Driven Eddy-Diffusivity Profile

    edyc(:) = 0.0

    if (kpblh > kbot) then

       deltar = 0.0
    
       if (ql(kpblh) > 1.e-6 .or. ql(kpblh-1) > 1.e-6) then
       
          deltar = fthrd(kpblh-1) * dzt(kpblh-1) * tair(kpblh-1) / theta(kpblh-1) &
                 + fthrd(kpblh  ) * dzt(kpblh  ) * tair(kpblh  ) / theta(kpblh  )
       else
          
          deltar = 0.0

       endif

       if (deltar < -1.e-6) then

          wcloud = (- grav / thetav(kpblh) * pblh * deltar)**0.33333333

          do k = kbot, kpblh-1
             zagl    = zm(k) - zm(kbot-1)
             zoh     = zagl / pblh
               
           ! If cloud effects limited to just cloud layer
           ! (see ECMWF physics technote):
           ! zagl    = max(zm(k) - zm(kbot-1) - zbase, 0.0)
           ! zoh     = zagl / (pblh - zbase)

             edyc(k) = 0.85 * vonk * wcloud * zagl * zoh * sqrt(1.0 - zoh)
          enddo
            
       endif
    endif

!! Vertical loop over T levels to compute the alpha and beta coefficients
!! in the buoyancy term of the moist Richardson number, from Cuipers and
!! Duynkerke (JAS, 1993, pp 3894-3908)

    do k = kbot, ktop

       if (ql(k) < 1.e-6) then

          alpha(k) = 1.0 + eps_virt * qv(k)
          beta (k) = eps_virt * theta(k)

       else

          alpha(k) = (1.0 - qw(k) + c1 * qv(k) * (1.0 + c2/tair(k))) / &
                     (1.0 + c3 * qv(k) / (tair(k) * tair(k)))
          beta(k)  = alpha(k) * alvlocp * theta(k) / tair(k) - theta(k)

       endif

    enddo

!! Louis's Richardson-number dependent eddy diffusivity
!! Vertical loop over W levels
      
    do k = kbot, ktop-1

       ! This uses the moist Richardson number of Brinkop and Roeckner 
       ! (Tellus, 1995, pp 197-222) computed from THIL and total water,
       ! but with the alpha and beta coefficients of Cuipers and Duynkerke
    
       ss = ( (vx(k+1) - vx(k))**2 &
            + (vy(k+1) - vy(k))**2 &
            + (vz(k+1) - vz(k))**2 ) * dzim(k) * dzim(k) + 1.e-9

       a  = 0.5 * (alpha(k+1) + alpha(k))
       b  = 0.5 * (beta (k+1) + beta (k))

       buoy = a * (thil(k+1) - thil(k)) + b * (qw(k+1) - qw(k))

       ri = grav2 / (thetav(k+1) + thetav(k)) * buoy * dzim(k) / ss

       zagl = zm(k) - zm(kbot-1)
       zk   = 0.4 * zagl
            
       if (ri >= 0.0) then

        ! if (zagl < pblh .and. wtv0 < 0.0) then
        !    fh  = max( (1. - zagl/pblh)**2, fhmin) * phih(k)**(-2)
        !    sql = zk**2
        ! else
             fh  = ( max( 1.0 - ri/rc, fhmin) )**2
             sql = ( zk * rlam / (rlam + zk) )**2
        ! endif

             edyr(k) = kzo + sqrt(ss) * fh * sql
       !     pran(k) = pr0 / ( 1.0 - (1.0 - pr0) * min(ri,rc) / rc )

       else

          sql     = ( zk * rlam / (rlam + zk) )**2
          edyr(k) = kzo + sqrt(ss * (1.0 - 25.0 * ri)) * sql
       !  pran(k) = pr0 / ( 1.0 - (1.0 - pr0) * min(ri,rc) / rc )

       endif

    enddo

!! Now computed combined eddy diffusivity from surface and cloud K-profiles
!! and Richardson-number eddy diffusivity

    do k = kbot, mza-2
       zkh(k) = max( edyr(k), edyz(k) + edyc(k), kzo )
       zkh(k) = min( zkh(k), 1000.0 )
       zkh(k) = 0.5 * (rho(k+1) + rho(k)) * zkh(k)
    enddo

    zkh(1:kbot-1) = 0.0
    zkh(mza-1)  = 0.0
    zkh(mza)    = 0.0

  end subroutine acm2_eddyx



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

    do k = kbot+1, ktop-1
       if (thv(k) > thvsfc) exit
    enddo
    kmix = k

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

    rib = 0.0

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

    kpblh = k

    ! INTERPOLATE BETWEEN LEVELS TO DETERMINE THE PBL HEIGHT

    fint = (ric - rib_sav) / (rib - rib_sav)

    if (fint > 0.5) then
       fint  = fint - 0.5
    else
       kpblh = kpblh - 1
       fint  = fint  + 0.5
    endif

    pblh = fint * (zm(kpblh) - zm(kpblh-1)) + zm(kpblh-1) - zm(kbot-1)

  end subroutine acm2_pblhgt

END MODULE module_bl_acm2
