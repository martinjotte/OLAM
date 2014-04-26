MODULE module_cu_gf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This convective parameterization is build to attempt     !
!     a smooth transition to cloud resolving scales as proposed!
!     by Arakawa et al (2011, ACP). It currently does not use  !
!     subsidencespreading as in G3. Difference and details     !
!     will be described in a forthcoming paper by              !
!     Grell and Freitas (2013). The parameterization also      !
!     offers options to couple with aerosols. Both, the smooth !
!     transition part as well as the aerosol coupling are      !
!     experimental. While the smooth transition part is turned !
!     on, nd has been tested dow to a resolution of about 3km  !
!     the aerosol coupling is turned off.                      !
!     More clean-up as well as a direct coupling to chemistry  !
!     will follow for V3.5.1                                   !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mem_para,  only: myrank
  use misc_coms, only: io6
  implicit none

  private :: myrank, io6

  integer, parameter, private :: iens    =  1
  integer, parameter, private :: maxiens =  1
  integer, parameter, private :: maxens  =  1
  integer, parameter, private :: maxens2 =  1
  integer, parameter, private :: maxens3 = 16
  integer, parameter, private :: ensdim  = 16

  integer, parameter, private :: its = 1, ite = 1, itf = 1
  integer, parameter, private :: jts = 1, jte = 1, jtf = 1

  real,    parameter, private :: xmbmax = 1.0

CONTAINS

  subroutine gf_driver(iw, dtlong)

     use mem_turb,    only: kpblh, frac_land, fthpbl, fqtpbl, wtv0, &
                            ustar, sflux_t, sflux_r
     use consts_coms, only: cp, alvl, grav, p00i, rocp, rvap, erad, &
                           gravi, alvlocp, vonk, r8
     use mem_radiate, only: rshort, fthrd_sw, fthrd_lw
     use mem_grid,    only: mza, lpw, arw0, zm, zt, xew, yew, zew, &
                           glatw, glonw, dzt, arw, lpv, arv, volt
     use mem_ijtabs,  only: itab_w
     use mem_basic,   only: wmc, vmc, theta, tair, press, rho, sh_v, &
                           vxe, vye, vze
     use mem_cuparm,  only: thsrc, rtsrc, conprr

     implicit none

     integer, intent(in)  :: iw
     real,    intent(in)  :: dtlong

     real,    parameter :: ccnclean = 250.
     real,    parameter :: beta     = 0.03

     ! CCN-dependent autoconversion of cloud water to rain:
     integer, parameter :: autoconv = 2 ! 1=old (no ccn), 2=berry c0 with ccn

     ! Include CCN effect on precipitation efficiency:
     integer, parameter :: aeroevap = 3 ! 1=old (no ccn), 2=ccn only, 3=average

     integer, parameter :: training = 0
     integer, parameter :: use_excess = 0
     integer, parameter :: use_excess_sh = 0

     real    :: xmb    (1)     ! deep conv convective mass flux
     real    :: zo     (1,mza) ! height
     real    :: outqc  (1,mza) ! output deep conv cloud water tendency
     integer :: j              ! dummy horizontal index
     real    :: aaeq   (1)     ! turn off convection if negative
     real    :: t      (1,mza) ! input temperature
     real    :: q      (1,mza) ! input mixing ratio
     real    :: z1     (1)     ! surface elevation
     real    :: sub_mas(1,mza)
     real    :: tn     (1,mza) ! input forced temperature
     real    :: qo     (1,mza) ! input forced mixing ratio
     real    :: po     (1,mza) ! input forced pressure
     real    :: pre    (1)     ! output precipitation rate
     real    :: p      (1,mza) ! input pressure
     real    :: outt   (1,mza) ! output deep conv temp tendency
     real    :: outq   (1,mza) ! output deep conv water vapor tendency
     real    :: dtime          ! time step
     integer :: ktau           ! random # seed, not used
     real    :: psur   (1)     ! input surface pressure
     real    :: us     (1,mza) ! input zonal wind
     real    :: vs     (1,mza) ! input meridional wind
     real    :: tcrit          ! temperature for cloud water phase
     real    :: ztexec (1)     ! PBL temperature excess
     real    :: zqexec (1)     ! PBL humidity excess
     real    :: ccn    (1)     ! ccn concentration ( #/cm^3 ?)
     real    :: r      (1,mza) ! air density
     real    :: dx             ! grid spacing
     real    :: mconv  (1,    1) ! integrated column moisture convergence
     real    :: omeg   (1,mza,1) ! vertical velocity in pressure coordinates
     integer :: k22    (1)     ! deep convectin updraft source level
     integer :: kbcon  (1)     ! deep convection LCL
     integer :: ktop   (1)     ! deep convection cloud top
     real    :: cupclw (1,mza) ! deep convection cloud water
     real    :: xf_ens (1,1,ensdim) ! mass flux ensembles
     real    :: pr_ens (1,1,ensdim) ! precipitation ensembles
     real    :: xland  (1)     ! 1 = land, 2 = water
     real    :: gsw    (1)     ! surface incoming radiation
     real    :: subt   (1,mza) ! temp tendency due to subsidence
     real    :: subq   (1,mza) ! water vapor tendency due to subsidence
     real    :: xl             ! latent heat of vaporization
     real    :: rv             ! gas constant for water vapor
     real    :: cpd            ! heat capacity of air
     real    :: g              ! gravitational acceleration
     integer :: ichoice        ! flag if only want one closure (usually set to zero)
     integer :: ipr,jpr        ! unused
     integer :: ens4           ! number of horizontal ensembles (not used with GF)

     character(50) :: ierrc(1) ! error messages

     real :: APR_GR   (1,1)    ! precip (mm/hr) for different closures/caps
     real :: APR_W    (1,1)
     real :: APR_MC   (1,1)
     real :: APR_ST   (1,1)
     real :: APR_AS   (1,1)
     real :: APR_CAPMA(1,1)
     real :: APR_CAPME(1,1)
     real :: APR_CAPMI(1,1)

     integer :: kts            ! starting k level
     integer :: ktf            ! ending   k level
     integer :: kte            ! max size of k dimensions

     ! shallow convection variables
     
     real    :: tshall (1,mza) ! PBL forced temperature for shallow convection
     real    :: qshall (1,mza) ! PBL forced water vapor for shallow convection
     real    :: xmbs   (1)     ! shallow conv convective mass flux
     real    :: dhdt   (1,mza) ! moist static energy tendency
     integer :: k22s   (1)     ! deep convectin updraft source level
     integer :: kbcons (1)     ! deep convection LCL
     integer :: ktops  (1)     ! deep convection cloud top
     real    :: cupclws(1,mza) ! deep convection cloud water
     real    :: tscl_kf        ! shallow convection timescale
     integer :: ierr(1)        ! error code
     integer :: kpbl (1)       ! layer corresponding to PBL height
     real    :: outts(1,mza)   ! shallow conv temp tendency
     real    :: outqs(1,mza)   ! shallow conv water vapor tendency

     ! local variables

     integer :: k, ka, kc, npoly, n, iwn, jv, iv
     real    :: raxis, omega_ave, dirv, flx, zws
     real    :: exner(mza)
     real    :: vflux(mza), vflux_the(mza), vflux_vap(mza)
     real    :: hflux, hflux_the, hflux_vap
     real    :: fqvadv, fthadv
     real(r8):: qsum, qav, tsum, tav

     ka    = lpw(iw)
     npoly = itab_w(iw)%npoly
     kpbl  = kpblh(iw) - ka + 1

     kts  = 1
     ktf  = mza - ka
     kte  = mza
     ens4 = 1

     ! Go no higher than 50mb for convective calculations to prevent any
     ! problems when esat gets near ambient pressure
     do k = ka, mza-2
        if (press(k,iw) < 50.e2) then
           ktf = k - ka + 1
           exit
        endif
     enddo

     ! A. Betts for shallow convection: suggestion for the KF timescale < DELTAX / 25 m/s
     tscl_kf = sqrt(arw0(iw)) / 25.

     dx = sqrt(arw0(iw))
     dtime = dtlong
     gsw   = rshort(iw)
     psur  = 0.01 * (press(ka,iw) + (zt(ka)-zm(ka-1))*rho(ka,iw)*grav)
     z1    = zm(ka-1)

     if (allocated(frac_land)) then
        xland = 2.0 - frac_land(iw)
     else
        xland = 1.0
     endif

     ccn = ccnclean
     if (xland(1) < 1.1) then
        ccn = 1000.0
     endif

     if (wtv0(iw) > 0.0 .and. (use_excess>=1 .or. use_excess_sh>=1)) then
        zws = vonk * grav * wtv0(iw) * (zt(ka)-zm(ka-1)) / theta(ka,iw)
        zws = (ustar(iw)**3 + zws) ** 0.33333333
        zws = max(zws,0.1)
        ztexec(1) = max( sflux_t(iw) / zws, 0.0)
        zqexec(1) = max( sflux_r(iw) / zws, 0.0)
     else
        ztexec(1) = 0.0
        zqexec(1) = 0.0
     endif

     raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

     ! Vertical advective theta and water vapor fluxes (W levels)

     do k = ka,mza-2
        vflux(k) = arw(k,iw) * wmc(k,iw)

        ! upwinded
        if (wmc(k,iw) >= 0.0) then
           vflux_vap(k) = vflux(k) * sh_v(k,iw)
           vflux_the(k) = vflux(k) * theta(k,iw)
        else
           vflux_vap(k) = vflux(k) * sh_v(k+1,iw)
           vflux_the(k) = vflux(k) * theta(k+1,iw)
        endif

        ! centered
        ! vflux_the(k) = vflux(k) * 0.5 * (theta(k,iw) + theta(k+1,iw))
        ! vflux_vap(k) = vflux(k) * 0.5 * (sh_v (k,iw) + sh_v (k+1,iw))
     enddo

     vflux    (ka-1)  = 0.
     vflux    (mza-1) = 0.
     vflux_the(ka-1)  = 0.
     vflux_the(mza-1) = 0.
     vflux_vap(ka-1)  = 0.
     vflux_vap(mza-1) = 0.

     mconv(1,1) = 0.0

     ! Loop over T levels

     do kc = 1, ktf
        k  = kc + ka - 1

        exner(k) = tair(k,iw) / theta(k,iw)

        ! Current temp, water vapor, density, and pressure

        t(1,kc) = max( tair(k,iw), 200.0 )
        q(1,kc) = max( sh_v(k,iw), 1.e-8 )
        r(1,kc) = rho(k,iw)
        p(1,kc) = 0.01 * press(k,iw)

        ! Height of T (half) levels

        zo(1,kc) = zt(k)

        ! Horizontal advective mass and water vapor fluxes

        hflux     = 0.
        hflux_the = 0.
        hflux_vap = 0.
      
        do jv = 1, npoly
           iv  = itab_w(iw)%iv(jv)

           if (k >= lpv(iv)) then
              iwn  = itab_w(iw)%iw(jv)
              dirv = itab_w(iw)%dirv(jv)

              flx   = dirv * vmc(k,iv) * arv(k,iv)
              hflux = hflux + flx

              ! upwinded
              if (flx >= 0.0) then
                 hflux_vap = hflux_vap + flx * sh_v (k,iwn)
                 hflux_the = hflux_the + flx * theta(k,iwn)
              else
                 hflux_vap = hflux_vap + flx * sh_v (k,iw)
                 hflux_the = hflux_the + flx * theta(k,iw)
              endif

              ! centered
              ! hflux_the = hflux_the + flx * 0.5 * (theta(k,iw) + theta(k,iwn))
              ! hflux_vap = hflux_vap + flx * 0.5 * (sh_v (k,iw) + sh_v (k,iwn))
           endif
        enddo

        fthadv = ((vflux_the(k-1) - vflux_the(k) + hflux_the) &
               - (vflux(k-1) - vflux(k) + hflux) * theta(k,iw)) &
               / (volt(k,iw) * rho(k,iw))

        fqvadv = ((vflux_vap(k-1) - vflux_vap(k) + hflux_vap) &
               - (vflux(k-1) - vflux(k) + hflux) * sh_v(k,iw)) &
               / (volt(k,iw) * rho(k,iw))

        ! "forced" temp, water vapor, and pressure

        tn(1,kc) = tair(k,iw) + (fthrd_lw(k,iw) + fthrd_sw(k,iw) + &
                                 fthpbl(k,iw) + fthadv) * dtlong * exner(k)
        qo(1,kc) = sh_v(k,iw) + (fqtpbl(k,iw) + fqvadv) * dtlong
        tn(1,kc) = max( tn(1,kc), 200.0 )
        qo(1,kc) = max( qo(1,kc), 1.e-8 )
        po(1,kc) = p(1,kc)

        ! U and V winds
        if (raxis > 1.e3) then
           us(1,kc) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) / raxis

           vs(1,kc) = vze(k,iw) * raxis / erad  &
                - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) / (raxis * erad)  
        else
           us(1,kc) = vxe(k,iw)
           vs(1,kc) = vye(k,iw)
        endif
        
        ! average vertical velocity from current and surrounding cells
        omeg(1,kc,1) = -grav * (wmc(k,iw) + sum(wmc(k,itab_w(iw)%iw(1:npoly)))) / real(npoly+1)

        ! moisture convergence
        mconv(1,1) = mconv(1,1) + omeg(1,kc,1) * (sh_v(k+1,iw) - sh_v(k,iw)) * gravi

     enddo

     mconv(1,1) = max(mconv(1,1), 0.0)

     APR_GR    = 0.0
     APR_W     = 0.0
     APR_MC    = 0.0
     APR_ST    = 0.0
     APR_AS    = 0.0
     APR_CAPMA = 0.0
     APR_CAPME = 0.0
     APR_CAPMI = 0.0

     ipr = iw
     jpr = 0

     xl  = alvl
     rv  = rvap
     cpd = cp
     g   = grav

     outt (1,:) = 0.0
     outq (1,:) = 0.0
     outqc(1,:) = 0.0
     subt (1,:) = 0.0
     subq (1,:) = 0.0
     pre  (1)   = 0.0

     tcrit = 258.0
     aaeq = 1.
     j = 1
     ichoice = 0
     xmb = 0.0
     sub_mas = 0.0
     xf_ens = 0.0
     pr_ens = 0.0
     k22 = 0
     kbcon = 0
     ktop = 0

     call CUP_gf(xmb,zo,OUTQC,J,AAEQ,T,Q,Z1,sub_mas,           &
          TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,PSUR,US,VS,      &
          TCRIT,ztexec,zqexec,ccn,ccnclean,r,dx,mconv,         &
          omeg,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,              &
          APR_CAPMA,APR_CAPME,APR_CAPMI,k22,kbcon,ktop,cupclw, &
          xf_ens,pr_ens,xland,gsw,subt,subq,                   &
          xl,rv,cp,g,ichoice,ipr,jpr,ierrc,ens4,               &
          beta,autoconv,aeroevap,ktf,training,                 &
          use_excess,kts,kte                                   )

     thsrc(:,iw) = 0.0
     rtsrc(:,iw) = 0.0
     conprr (iw) = 0.0

     if (pre(1) > 1.e-16) then

        ! Deep convecton is active. Copy tendencies to model arrays.

        do kc = 1, ktf
           k  = kc + ka - 1

           ! Total water tendency
           outq(1,kc) = outq(1,kc) + subq(1,kc) + outqc(1,kc)

           ! Any cloud condensate is evaporated since we do not feed back
           ! to resolved microphysics
           outt(1,kc) = outt(1,kc) + subt(1,kc) - alvlocp * outqc(1,kc)
        enddo

        ! Slightly modify tendencies to ensure heat and moisture conservation

        qav  = sum(     outq(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) ) + pre(1)
        qsum = sum( abs(outq(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )

        tav  = cp * sum(     outt(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) ) - pre(1) * alvl
        tsum = cp * sum( abs(outt(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
        
        qav = qav / max(qsum, 1.e-20_r8)
        tav = tav / max(tsum, 1.e-20_r8)
        
        do kc = 1, ktf
           k  = kc + ka - 1
           rtsrc(k,iw) =  outq(1,kc) - qav * abs(outq(1,kc))
           thsrc(k,iw) = (outt(1,kc) - tav * abs(outt(1,kc))) / exner(k)
        enddo

        conprr(iw) = pre(1)

     else

        ! If there is no deep convection, check for shallow convection

        outts(1,:) = 0.0
        outqs(1,:) = 0.0
        outqc(1,:) = 0.0
        pre  (1)   = 0.0

        cupclws = 0.
        xmbs = 0.
        kbcons = 0
        ktops = 0
        k22s = 0
        ierr = 0
        ierrc = ""

        ! Shallow convection uses current T and Q plus PBL tendencies

        do kc = 1, ktf
           k  = kc + ka - 1
           tshall(1,kc) = t(1,kc) + fthpbl(k,iw) * dtlong * exner(k)
           qshall(1,kc) = q(1,kc) + fqtpbl(k,iw) * dtlong
           dhdt  (1,kc) = cp * fthpbl(k,iw) + alvl * fqtpbl(k,iw)
        enddo

        call CUP_gf_sh(xmbs,zo,OUTQC,J,AAEQ,T,Q,Z1,      &
              tshall,qshall,PO,PRE,P,OUTTS,OUTQS,DTIME,  &
              ktau,PSUR,US,VS,TCRIT,                     &
              ztexec,zqexec,ccn,ccnclean,r,dx,dhdt,      &
              kpbl,kbcons,ktops,cupclws,k22s,            &
              xland,gsw,tscl_kf,                         &
              xl,rv,cp,g,ichoice,ipr,jpr,ierr,ierrc,     &
              autoconv,ktf,use_excess_sh,kts,kte         )

        if (kbcons(1) > 0 .and. ktops(1) >= kbcons(1)) then

           ! Shallow convection is active; copy tendencies

           ! Slightly modify tendencies to ensure heat and moisture conservation

           qav  = sum(     outqs(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
           qsum = sum( abs(outqs(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )

           tav  = cp * sum(     outts(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
           tsum = cp * sum( abs(outts(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )

           qav = qav / max(qsum, 1.e-20_r8)
           tav = tav / max(tsum, 1.e-20_r8)

           do kc = 1, ktf
              k  = kc + ka - 1
              rtsrc(k,iw) =  outqs(1,kc) - qav * abs(outqs(1,kc))
              thsrc(k,iw) = (outts(1,kc) - tav * abs(outts(1,kc))) / exner(k)
           enddo

        endif

     endif

   end subroutine gf_driver



   SUBROUTINE CUP_gf(xmb_out,zo,OUTQC,J,AAEQ,T,Q,Z1,sub_mas,       &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,PSUR,US,VS,      &
              TCRIT,ztexec,zqexec,ccn,ccnclean,rho,dx,mconv,       &
              omeg,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,              &
              APR_CAPMA,APR_CAPME,APR_CAPMI,k22,kbcon,ktop,cupclw, &
              xf_ens,pr_ens,xland,gsw,subt,subq,                   &
              xl,rv,cp,g,ichoice,ipr,jpr,ierrc,ens4,               &
              beta,autoconv,aeroevap,ktf,training,                 &
              use_excess,kts,kte                                   )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,aeroevap,ktf,ktau,training,use_excess,        &
        kts,kte,ipr,jpr,ens4
     integer, intent (in   )              ::                           &
        j,ichoice
  !
  ! 
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        xf_ens,pr_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout )                  ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,     &
               APR_CAPME,APR_CAPMI
    real, dimension( its:ite , jts:jte )                               &
          :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw

  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        OUTT,OUTQ,OUTQC,subt,subq,sub_mas,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre, xmb_out
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop,k22
!    integer,    dimension (its:ite)                                   &
!       ,intent (in  )                   ::                           &
!       kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,T,PO,P,US,VS,tn
     real,    dimension (its:ite,kts:kte,1:ens4)                       &
        ,intent (inout   )                   ::                           &
        omeg
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ztexec,zqexec,ccn,Z1,PSUR,AAEQ,xland
     real, dimension (its:ite,1:ens4)                                         &
        ,intent (in   )                   ::                           &
        mconv

       
       real                                                            &
        ,intent (in   )                   ::                           &
        beta,dx,ccnclean,dtime,tcrit,xl,cp,rv,g


!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens,subt_ens,subq_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,zuhe
     real,    dimension (its:ite,1:ens4) ::                                   &
        axx
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,entr_rate,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: xxx,xx1,xx2,power_entr,zustart,zufinal,dzm1,dzp1


     integer :: k1,k2,kbegzu,kfinalzu,kstart,jmini,levadj
     logical :: keep_going
     real xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (kts:kte) :: smth

      levadj=5
      power_entr=2. ! 1.2
      zustart=.1
      zufinal=1.
      day=86400.
      do i=its,itf
        xmb_out(i)=0.
        closure_n(i)=16.
        xland1(i)=1.
        if(xland(i).gt.1.5)xland1(i)=0.
        cap_max_increment(i)=25.
        ierrc(i)=" "
!       cap_max_increment(i)=1.
      enddo
!
!--- specify entrainmentrate and detrainmentrate
!--- highly tuneable !
!
      entr_rate=7.e-5
      radius=.2/entr_rate
      frh=3.14*(radius*radius)/dx/dx
      if(frh .gt. 0.7)then
         frh=.7
         radius=sqrt(frh*dx*dx/3.14)
         entr_rate=.2/radius
      endif
      do i=its,itf
         sig(i)=(1.-frh)**2
      enddo
!      sig(:)=1.

!
!--- entrainment of mass
!
      mentrd_rate=entr_rate ! 0.
      xlamdd=mentrd_rate
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        cupclw(i,k)=0.
        cd(i,k)=1.*entr_rate
        cdd(i,k)=xlamdd
        hcdo(i,k)=0.
        qrcdo(i,k)=0.
        dellaqc(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax=1.
      edtmin=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=1000.   ! gg 500
!
!--- maximum depth (mb) of capping inversion that convection
!--- can ovecome (smaller cap_max = less convection)
!
!     cap_maxs=75.
!     cap_maxs=50.
      cap_maxs=40.

      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
        edt(i)=0.
        kstabm(i)=ktf-1
        IERR(i)=0
        IERR2(i)=0
        IERR3(i)=0
      enddo
      do i=its,itf
          cap_max(i)=cap_maxs
        iresult=0

      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1250.   !1000
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           ktf,kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
      endif
      enddo

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22,ierr,ktf,kts,kte)
      DO 36 i=its,itf
         IF(ierr(I).eq.0)THEN
           ! frh=q_cup(i,k22(i))/qes_cup(i,k22(i))
           ! IF(omeg(i,k22(i),1).lt.0. .and. frh.ge.0.99 .and. sig(i).lt.0.091)ierr(i)=1200
           IF(K22(I).GE.KBMAX(i))THEN
             ierr(i)=2
             ierrc(i)="could not find k22"
           ENDIF
         ENDIF
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      do i=its,itf
       IF(ierr(I).eq.0)THEN
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)+1
             hkb(i) =he_cup(i,k22(i)) ! sum(he_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
        else if(use_excess <= 1)then
         hkb(i)=he_cup(i,k22(i)) ! +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
         hkbo(i)=heo_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
        endif  ! excess
       endif ! ierr
      enddo


      call cup_kbcon(ierrc,cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess,ktf,kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,ktf,kts,kte)
      !DO i=its,itf
      !   IF(ierr(I).eq.0)THEN
      !   do k=k22(i),kbcon(i)
      !   frh=q_cup(i,k)/qes_cup(i,k)
      !   if(omeg(i,k,1).lt.-1.e-6 .and. frh.ge.0.99 .and. sig(i).lt.0.091)ierr(i)=1200
      !   enddo
      !   endif
      !enddo
!
! the following section insures a smooth normalized mass flux profile. See Grell
! and Freitas (2013) for a description
!
      DO i=its,itf
         IF(ierr(I).eq.0)THEN
            do k=kts,ktf
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate*(1.3-frh)
            enddo
            zuhe(i)=zustart
            kstart=1
            xx1=z_cup(i,kbcon(i))*z_cup(i,kbcon(i))
            xx2=z_cup(i,kstart)*z_cup(i,kstart)
            frh=(zufinal-zustart)/(xx1-xx2)
            dh=zustart-frh*xx2
            do k=kstart,kbcon(i)-1
             dz=z_cup(i,k+1)-z_cup(i,k)
!            cd(i,k)=entr_rate_2d(i,kbcon(i))
             if(p_cup(i,k).gt. p_cup(i,kstabi(i)))cd(i,k)=1.e-6
             xxx=z_cup(i,k+1)*z_cup(i,k+1)
             entr_rate_2d(i,k)=((frh*xxx+dh)/zuhe(i)-1.+cd(i,k)*dz)/dz
             zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
            enddo
            kbegzu=kstabi(i)+4
            kbegzu=min(kbegzu,ktf-2)
            kfinalzu=kbegzu+1
            !do k=kts,ktf
            !   cd(i,k)=entr_rate_2d(i,kbcon(i))
            !enddo
               do k=kbcon(i),kbegzu
                cd(i,k)=entr_rate_2d(i,kbcon(i))
                if(p_cup(i,k).gt. p_cup(i,kstabi(i)))cd(i,k)=1.e-6
                dz=z_cup(i,k+1)-z_cup(i,k)
                zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
               enddo
         do k=kstabi(i),ktf-2
          if((hkb(i)-hes_cup(i,k)).lt.0.)then
              kfinalzu=k-3
              go to 411
          endif
         enddo
411      continue
             kfinalzu=max(kfinalzu,kbegzu+1)
             kfinalzu=min(kfinalzu,ktf-1)
             xx1=z_cup(i,kfinalzu)*z_cup(i,kfinalzu)
             xx2=z_cup(i,kbegzu)*z_cup(i,kbegzu)
             frh=-(0.2-zuhe(i))/(xx1-xx2)
             dh=zuhe(i)+frh*xx2
               do k=kbegzu+1,kfinalzu
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 xxx=z_cup(i,k+1)*z_cup(i,k+1)
                 cd(i,k)=-((-frh*xxx+dh)/zuhe(i)-1.-entr_rate_2d(i,k)*dz)/dz
                 zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
               enddo
               do k=kfinalzu+1,ktf
                   cd(i,k)=entr_rate_2d(i,k)
               enddo
               do k=kts+1,ktf-2
                 dzm1=z_cup(i,k)-z_cup(i,k-1)
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 dzp1=z_cup(i,k+2)-z_cup(i,k+1)
                 smth(k)=.25*(dzm1*cd(i,k-1)+2.*dz*cd(i,k)+dzp1*cd(i,k+1))
               enddo
               do k=kts+1,ktf-2
                 dzm1=z_cup(i,k)-z_cup(i,k-1)
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 dzp1=z_cup(i,k+2)-z_cup(i,k+1)
                 cd(i,k)=smth(k)/dz ! (.25*(dzm1+2.*dz+dzp1))
               enddo

            smth(:)=0.
            do k=2,ktf-2
                 dzm1=z_cup(i,k)-z_cup(i,k-1)
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 dzp1=z_cup(i,k+2)-z_cup(i,k+1)
              smth(k)=.25*(dzm1*entr_rate_2d(i,k-1)+2.*dz*entr_rate_2d(i,k)+dzp1*entr_rate_2d(i,k+1))
            enddo
            do k=2,ktf-2 
                 dz=z_cup(i,k+1)-z_cup(i,k)
              entr_rate_2d(i,k)=smth(k)/dz
            enddo
            zuhe(i)=zustart
            do k=2,kbegzu 
              dz=z_cup(i,k+1)-z_cup(i,k)
              frh=zuhe(i)
              zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
            enddo
         ENDIF
       enddo
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
         hco(i,k)=0.
         DBYo(I,K)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
       endif ! ierr
      enddo
!
!
      do i=its,itf
         if(ierr(i).eq.0)then
         zu(i,1)=zustart
         zuo(i,1)=zustart
!    mass entrainment and detrinament is defined on model levels
         do k=2,ktf-1
          dz=zo_cup(i,k)-zo_cup(i,k-1)
          up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
          up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
          zuo(i,k)=zuo(i,k-1)+up_massentro(i,k-1)-up_massdetro(i,k-1)
          if(zuo(i,k).lt.0.05)then
             zuo(i,k)=.05
             up_massdetro(i,k-1)=zuo(i,k-1)-.05  + up_massentro(i,k-1)
             cd(i,k-1)=up_massdetro(i,k-1)/dz/zuo(i,k-1)
          endif
          zu(i,k)=zuo(i,k)
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo
         do k=kbcon(i)+1,ktf-1
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo
         do k=kbcon(i)+1,ktf
          if(dbyo(i,k).lt.0.)then
              ktop(i)=k-1
              go to 41
          endif
         enddo
41       continue
         if(ktop(i).lt.kbcon(i)+2)then
            ierr(i)=5
            ierrc(i)='ktop too small'
         endif
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
      endif
      enddo
!
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr,ktf,kts,kte)
      DO 100 i=its,itf
         IF(ierr(I).eq.0)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
         jmini = jmin(i)
         keep_going = .TRUE.
         do while ( keep_going )
           keep_going = .FALSE.
           if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
           if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
           ki = jmini
           hcdo(i,ki)=heso_cup(i,ki)
           DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
           dh=0.
           do k=ki-1,1,-1
             hcdo(i,k)=heso_cup(i,jmini)
             DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
             dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
             if(dh.gt.0.)then
               jmini=jmini-1
               if ( jmini .gt. 5 ) then
                 keep_going = .TRUE.
               else
                 ierr(i) = 9
                 ierrc(i) = "could not find jmini9"
                 exit
               endif
             endif
           enddo
         enddo
         jmin(i) = jmini 
         if ( jmini .le. 5 ) then
           ierr(i)=4
           ierrc(i) = "could not find jmini4"
         endif
       ENDIF
100   continue
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
         IF(ierr(I).eq.0)THEN
            IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
               ierr(i)=6
               ierrc(i)="cloud depth very shallow"
            endif
         endif
      enddo

!
!--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      do k=kts,ktf
      do i=its,itf
       zd(i,k)=0.
       zdo(i,k)=0.
       cdd(i,k)=0.
       dd_massentr(i,k)=0.
       dd_massdetr(i,k)=0.
       dd_massentro(i,k)=0.
       dd_massdetro(i,k)=0.
       hcdo(i,k)=heso_cup(i,k)
       dbydo(i,k)=0.
      enddo
      enddo
      do i=its,itf
          bud(i)=0.
          IF(ierr(I).eq.0)then
            mentrd_rate_2d(i,:)=mentrd_rate
            cdd(i,1:jmin(i))=xlamdd
            cdd(i,jmin(i))=0.
! start from dd origin
            zd(i,jmin(i))=0.2
            zdo(i,jmin(i))=0.2
            xx1=z_cup(i,jmin(i)-levadj)*z_cup(i,jmin(i)-levadj)
            xx2=z_cup(i,jmin(i))*z_cup(i,jmin(i))
            frh=(zdo(i,jmin(i))-1.)/(-xx1+xx2)
            dh=zdo(i,jmin(i))-frh*xx2
            zuhe(i)=zdo(i,jmin(i))
            do ki=jmin(i)-1,jmin(i)-levadj,-1
             cdd(i,ki)=0.
             dz=z_cup(i,ki+1)-z_cup(i,ki)
             xxx=z_cup(i,ki)*z_cup(i,ki)
             mentrd_rate_2d(i,ki)=((frh*xxx+dh)/zuhe(i)-1.)/dz
             zuhe(i)=zuhe(i)+mentrd_rate_2d(i,ki)*dz*zuhe(i)
            enddo
! now we know the max zd, for detrainment we will go back to beta at level 1
            kstart=max(kbcon(i),kdet(i))-1
            kstart=min(jmin(i)-levadj,kstart)
            kstart=max(2,kstart)
            if(kstart.lt.jmin(i)-levadj-1)then
              do ki=jmin(i)-levadj-1,kstart,-1
                dz=z_cup(i,ki+1)-z_cup(i,ki)
                mentrd_rate_2d(i,ki)=mentrd_rate
                cdd(i,ki)=xlamdd
                zuhe(i)=zuhe(i)-cdd(i,ki)*dz*zuhe(i)+mentrd_rate_2d(i,ki)*dz*zuhe(i)
              enddo
            endif
            xx1=z_cup(i,kstart)*z_cup(i,kstart)
            xx2=z_cup(i,1)*z_cup(i,1)
            frh=(zuhe(i)-beta)/(xx1-xx2)
            dh=beta-frh*xx2
            mentrd_rate_2d(i,kstart)=0.
            do ki=kstart+1,1,-1
             mentrd_rate_2d(i,ki)=0.
             dz=z_cup(i,ki+1)-z_cup(i,ki)
             xxx=z_cup(i,ki)*z_cup(i,ki)
             cdd(i,ki)=max(0.,(1.-(frh*xxx+dh)/zuhe(i))/dz)
             zuhe(i)=zuhe(i)-cdd(i,ki)*dz*zuhe(i)
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'low cd ',ki,zuhe(i),cdd(i,ki)
            enddo

! now that we have entrainment and detrainment rates, 
! calculate downdraft mass terms
!
            do ki=jmin(i)-1,1,-1
               mentrd_rate=mentrd_rate_2d(i,ki)
               dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
               dd_massentro(i,ki)=mentrd_rate*dzo*zdo(i,ki+1)
               dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
               zdo(i,ki)=zdo(i,ki+1)+dd_massentro(i,ki)-dd_massdetro(i,ki)
            enddo
! downdraft moist static energy + moisture budget
            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))
            do ki=jmin(i)-1,1,-1
             dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
             hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+ &
                        dd_massentro(i,ki)*heo(i,ki))   /            &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
             bud(i)=bud(i)+dbydo(i,ki)*dzo
            enddo
          endif

        if(bud(i).gt.0.)then
          ierr(i)=7
          ierrc(i)='downdraft is not negatively buoyant '
        endif
      enddo
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_new(ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,xl,ktf,kts,kte)
!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture('deep',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           ccnclean,p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,        &
           ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum,psumh,&
           autoconv,aeroevap,1,ktf,j,ipr,jpr,kts,kte)
      do k=kts,ktf
      do i=its,itf
         cupclw(i,k)=qrco(i,k)
      enddo
      enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,ktf,kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,ktf,kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
!mjo       if(aa1(i).eq.0.)then
           if(aa1(i).lt.1.e-12)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       do i=1,ens4
       axx(:,i)=aa1(:)
       enddo

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwo,ccn,pwevo,edtmax,edtmin,edtc,psum,psumh, &
           ccnclean,rho,aeroevap,ktf,j,ipr,jpr,kts,kte)
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=edtc(i,iedt)
         edto(i)=edtc(i,iedt)
         edtx(i)=edtc(i,iedt)
         if(maxens2.eq.3)then
            edt(i)=edtc(i,maxens2)
            edto(i)=edtc(i,maxens2)
            edtx(i)=edtc(i,maxens2)
         endif
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
!
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellah(i,k)=0.
        dsubt(i,k)=0.
        dsubh(i,k)=0.
        dellaq(i,k)=0.
        dsubq(i,k)=0.
      enddo
      enddo
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp
         dellaq(i,1)=(edto(i)*zdo(i,2)*qrcdo(i,2)   &
                     -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
         dsubt(i,1)=0.
         dsubq(i,1)=0.

         do k=kts+1,ktop(i)
! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
! detrainment and entrainment for fowndrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
!
!         SPECIAL LEVELS
!
            if(k.eq.jmin(i))then
               entdoj=edto(i)*zdo(i,k)
            endif
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
            endif
            totmas=subin-subdown+detup-entup-entdo+ &
             detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
!               print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               write(0,123)k,subin,subdown,detup,entup, &
                           detdo,entdo,entupk,detupk
123     formAT(10X,i2,8E12.4)
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k)=(detup*.5*(HCo(i,K+1)+HCo(i,K)) &
                    +detdo*.5*(HCDo(i,K+1)+HCDo(i,K)) &
                    -entup*heo(i,k) &
                    -entdo*heo(i,k) &
                    +subin*heo_cup(i,k+1) &
                    -subdown*heo_cup(i,k) &
                    +detupk*(hco(i,ktop(i))-heo_cup(i,ktop(i)))    &
                    -entupk*heo_cup(i,k22(i)) &
                    -entdoj*heo_cup(i,jmin(i)) &
                     )*g/dp
            dellaq(i,k)=(detup*.5*(qco(i,K+1)+qco(i,K)-qrco(i,k+1)-qrco(i,k)) &
                    +detdo*.5*(qrcdo(i,K+1)+qrcdo(i,K)) &
                    -entup*qo(i,k) &
                    -entdo*qo(i,k) &
                    +subin*qo_cup(i,k+1) &
                    -subdown*qo_cup(i,k) &
                    +detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))    &
                    -entupk*qo_cup(i,k22(i)) &
                    -entdoj*qo_cup(i,jmin(i)) &
                     )*g/dp
!
! updraft subsidence only
!
           if(k.lt.ktop(i))then
             dsubt(i,k)=(zuo(i,k+1)*heo_cup(i,k+1) &
                    -zuo(i,k)*heo_cup(i,k))*g/dp
             dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1) &
                    -zuo(i,k)*qo_cup(i,k))*g/dp
           endif
!
       enddo   ! k

        endif
      enddo
!
!-- take out cloud liquid water for detrainment
!
      do k=kts,ktf-1
      do i=its,itf
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k)-po_cup(i,k+1))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*up_massdetro(i,k)*.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k)-po_cup(i,k+1))
         endif
         dellaqc(i,k)=max(0.,dellaqc(i,k))
       endif
      enddo
      enddo
!
!--- using dellas, calculate changed environmental profiles
!
      mbdt=mbdt_ens(1)
      do i=its,itf
      xaa0_ens(i,:)=0.
      enddo

      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
!           if(i.eq.ipr.and.j.eq.jpr.and.k.eq.kts)write(0,*)'mbdt = ',mbdt,mbdt_ens,dtime
            dsubh(i,k)=dsubt(i,k)
            XHE(I,K)=(dsubt(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
!           XQ(I,K)=(dsubq(i,k)+DELLAQ(I,K)+dellaqc(i,k))*MBDT+QO(I,K)
            XQ(I,K)=(dsubq(i,k)+DELLAQ(I,K))*MBDT+QO(I,K)
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            dSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
!           XT(I,K)= (DELLAT(I,K)+dsubt(i,k)-dellaqc(i,k)*xl/cp)*MBDT+TN(I,K)
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k))*MBDT+TN(I,K)
!           IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
            IF(XQ(I,K).LT.1.E-08)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      xhkb(i)=hkbo(i)+(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
!     IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      IF(XQ(I,ktf).LT.1.E-08)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,ktf,kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
!     do i=its,itf
!       if(ierr(i).eq.0)then
!         xhkb(i)=xhe(i,k22(i))
!       endif
!     enddo
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
!        if(use_excess == 2) then
!            k1=max(1,k22(i)-1)
!            k2=max(1,min(kbcon(i)-1,k22(i)+1))
!            k1=1
!            k2=k22(i)+1
!            xhkb(i) =sum(xhe_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
!        else if(use_excess <= 1) then
!            xhkb(i)=xhe_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))

!        endif
         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
         k=kbcon(i)
         xhc(i,k)=xhkb(i)
         xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
      xzu(i,:)=zuo(i,:)
      do k=kbcon(i)+1,ktop(i)
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,ktf,kts,kte)
      do 200 nens=1,maxens
      do i=its,itf 
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +edto(i)*pwdo(i,k)             &
                                    +pwo(i,k) 
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +.5*edto(i)*pwdo(i,k)          &
                                 +  pwo(i,k)
                 else
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k) ! +edto(i)*pwdo(i,k)
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,j,nall+7).lt.1.e-6)then
            ierr(i)=18
            ierrc(i)="total normalized condensate too small"
!           if(i.eq.ipr.and.j.eq.jpr)write(0,*)ierr(i),ierrc(i)
            do nens3=1,maxens3
               pr_ens(i,j,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,j,nall+nens3).lt.1.e-4)then
            pr_ens(i,j,nall+nens3)=0.
           endif
         enddo
         endif
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'ierrc = ',ierr(i),ierrc(i)
      enddo
 200  continue
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero, assuming this 
!        ensemble is chosen
!
!
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
         k22x(i)=k22(i)
      enddo
       if(maxens.gt.0)then
!     CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
!          ktf,kts,kte)
      call cup_kbcon(ierrc,cap_max_increment,2,k22x,kbconx,heo_cup, &
           heso_cup,hkbo,ierr2,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess,       &
           ktf,kts,kte)
      call cup_kbcon(ierrc,cap_max_increment,3,k22x,kbconx,heo_cup, &
           heso_cup,hkbo,ierr3,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess,       &
           ktf,kts,kte)
      endif
!
!--- calculate cloud base mass flux
!

      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,iedt,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           ichoice,ipr,jpr,ktf,kts,kte,ens4,ktau)
!
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           subt_ens(i,k,iedt)=dsubt(i,k)
           subq_ens(i,k,iedt)=dsubq(i,k)
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edt(i)*pwdo(i,k)
        else 
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!
       call cup_output_ens_3d(xf_ens,ierr,dellat_ens,dellaq_ens, &
            dellaqc_ens,subt_ens,subq_ens,subt,subq,outt,     &
            outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
            j,'deep',ierr2,ierr3,         &
            pr_ens,                    &
            sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                &
            APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,   &
            weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training, &
            ipr,jpr,ktf,kts,kte  )
      k=1
      do i=its,itf
         if(ierr(i).eq.0) PRE(I)=MAX(PRE(I),0.)
         if(ierr(i).eq.0) xmb_out(i)=xmb(i)
      enddo
!
!---------------------------done------------------------------
!
   END SUBROUTINE CUP_gf


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,edtc,psum2,psumh, &
              ccnclean,rho,aeroevap,ktf,j,ipr,jpr,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        j,ipr,jpr,aeroevap,ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws
     real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
     prop_c=8. !10.386
     alpha3 = 1.9
     beta3  = -1.13
     pefc=0.

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.0.9)pef=0.9
            if(pef.lt.0.1)pef=0.1
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
            if(pefb.gt.0.9)pefb=0.9
            if(pefb.lt.0.1)pefb=0.1
            EDT(I)=1.-.5*(pefb+pef)
            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1.)) !*1.e6
!              if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
!              prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1.)) !*1.e6
!              if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
!!             if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
               if(aeroevap.eq.3)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif


!--- edt here is 1-precipeff!
            if (maxens2==1) then
               edtc(i,1) = edt(i)
            else
               einc=.2*edt(i)
               do k=1,maxens2
                  edtc(i,k)=edt(i)+float(k-2)*einc
               enddo
            endif
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,           &
              ktf,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf,kts,kte
  ! cdd= detrainment function 
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he 
     real                                                              &
        ,intent (in   )                   ::                           &
        xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      DH=HCD(I,k)-HES_cup(I,K)
      if(dh.lt.0.)then
        QRCD(I,K)=(qes_cup(i,k)+(1./XL)*(GAMMA_cup(i,k) &
                  /(1.+GAMMA_cup(i,k)))*DH)
        else
          qrcd(i,k)=qes_cup(i,k)
        endif
      pwd(i,jmin(i))=zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
      qcd(i,k)=qrcd(i,k)
      pwev(i)=pwev(i)+pwd(i,jmin(i))
!
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1)                          &
                  -.5*dd_massdetr(i,ki)*qcd(i,ki+1)+ &
                  dd_massentr(i,ki)*q(i,ki))   /            &
                  (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
!        write(0,*)'qcd in dd_moi = ',qcd(i,ki)

!
!--- to be negatively buoyant, hcd should be smaller than hes!
!--- ideally, dh should be negative till dd hits ground, but that is not always
!--- the case
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)then
          dqeva=0.
          qrcd(i,ki)=qcd(i,ki)
         endif
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_new

   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              ktf,kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z,q
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar
!      real, external :: satvap
!      real :: satvap


      HT(1)=XL/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
       QES(I,K)=.622*E/(100.*P(I,K)-(1.-.622)*E)
!        e=satvap(t(i,k))
!        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
!       IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
        IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.le.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              ktf,kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts,ktf
      do i=its,itf
        qes_cup(i,k)=0.
        q_cup(i,k)=0.
        hes_cup(i,k)=0.
        he_cup(i,k)=0.
        z_cup(i,k)=0.
        p_cup(i,k)=0.
        t_cup(i,k)=0.
        gamma_cup(i,k)=0.
      enddo
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
!       hes_cup(i,1)=hes(i,1)
!       he_cup(i,1)=he(i,1)
        hes_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
        he_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*q(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        z_cup(i,1)=z1(i)
        p_cup(i,1)=psur(i)
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,iedt,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,      &
              icoic,ipr,jpr,ktf,kts,kte,ens4,ktau                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ipr,jpr,ktf,kts,kte,ens4,ktau
     integer, intent (in   )              ::                           &
        j,iedt
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte,1:ens4)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,xland
     real,    dimension (its:ite,1:ens4)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,kclim
     real                                 ::                           &
       fens4,a1,massfld,a_ave,xff0,xff00,xxx,xomg

     integer :: nall2,ixxx
     integer,  dimension (8) :: seed
     real, dimension (its:ite) :: ens_adj
     
     integer, parameter :: irandom = 0

       seed=0
       do i=its,itf
        if(ierr(i).eq.0)then
          seed(1)=int(aa0(i))
          seed(2)=int(aa1(i))
          exit
        endif
       enddo

       fens4=float(ens4)
!
!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
          ens_adj(i)=1.
          if(ierr2(i).gt.0.and.ierr3(i).eq.0)ens_adj(i)=0. ! 2./3.
          if(ierr2(i).gt.0.and.ierr3(i).gt.0)ens_adj(i)=0.
!
!---
!
             if(name.eq.'deeps')then

                a_ave=0.
                do ne=1,ens4
                  a_ave=a_ave+axx(i,ne)
                enddo
                a_ave=max(0.,a_ave/fens4)
                a_ave=min(a_ave,aa1(i))
                a_ave=max(0.,a_ave)
                do ne=1,16
                  xff_ens3(ne)=0.
                enddo
                xff0= (AA1(I)-AA0(I))/DTIME
                xff_ens3(1)=max(0.,(AA1(I)-AA0(I))/dtime)
                xff_ens3(2)=max(0.,(a_ave-AA0(I))/dtime)

                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(3)=max(0.,(axx(i,ixxx)-AA0(I))/dtime)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(13)=max(0.,(axx(i,ixxx)-AA0(I))/dtime)
                else
                   xff_ens3(3)=max(0.,(AA1(I)-AA0(I))/dtime)
                   xff_ens3(13)=max(0.,(AA1(I)-AA0(I))/dtime)
                endif
!   
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(14)=0.
                do ne=1,ens4
                  xff_ens3(14)=xff_ens3(14)-omeg(i,k22(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(14).lt.0.)xff_ens3(14)=0.
                xff_ens3(5)=0.
                do ne=1,ens4
                  xff_ens3(5)=xff_ens3(5)-omeg(i,kbcon(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
!  
! minimum below kbcon
!
                   xff_ens3(4)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.lt.xff_ens3(4))xff_ens3(4)=xomg
                   enddo
                   enddo
                   if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
!
! max below kbcon
                   xff_ens3(6)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
                   enddo
                   if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
!
!--- more like Krishnamurti et al.; pick max and average values
!
                xff_ens3(7)=mconv(i,1)
                xff_ens3(8)=mconv(i,1)
                xff_ens3(9)=mconv(i,1)
                if(ens4.gt.1)then
                   do ne=2,ens4
                      if (mconv(i,ne).gt.xff_ens3(7))xff_ens3(7)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      if (mconv(i,ne).lt.xff_ens3(8))xff_ens3(8)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      xff_ens3(9)=xff_ens3(9)+mconv(i,ne)
                   enddo
                   xff_ens3(9)=xff_ens3(9)/fens4
                endif
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(15)=mconv(i,ixxx)
                else
                   xff_ens3(15)=mconv(i,1)
                endif
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=AA0(i)/(60.*40.)
                xff_ens3(11)=AA0(I)/(60.*40.)
                xff_ens3(16)=AA0(I)/(60.*40.)
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AA0(I)/(60.*40.)
                else
                   xff_ens3(12)=AA1(I)/(60.*40.)
                endif

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(1)
                   if(xk(nens).le.0.0.and.xk(nens).gt.-1.e-2) &
                           xk(nens)=-1.e-2
                   if(xk(nens).gt.0.0.and.xk(nens).lt.1.e-2) &
                           xk(nens)=1.e-2
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
! over water, enforce small cap for some of the closures
!
                if(maxens.gt.0 .and. xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =ens_adj(i)*xff_ens3(1)
                      xff_ens3(2) =ens_adj(i)*xff_ens3(2)
                      xff_ens3(3) =ens_adj(i)*xff_ens3(3)
                      xff_ens3(13) =ens_adj(i)*xff_ens3(13)
                      xff_ens3(10) =ens_adj(i)*xff_ens3(10)
                      xff_ens3(11) =ens_adj(i)*xff_ens3(11)
                      xff_ens3(12) =ens_adj(i)*xff_ens3(12)
                      xff_ens3(16) =ens_adj(i)*xff_ens3(16)
                      xff_ens3(7) =ens_adj(i)*xff_ens3(7)
                      xff_ens3(8) =ens_adj(i)*xff_ens3(8)
                      xff_ens3(9) =ens_adj(i)*xff_ens3(9)
                      xff_ens3(15) =ens_adj(i)*xff_ens3(15)
!                     xff_ens3(7) =0.
!                     xff_ens3(8) =0.
!                     xff_ens3(9) =0.
                 endif
                endif
!
! end water treatment
!
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!

                      if(xff0.gt.0. .and. xk(ne).lt.0.)then
                         if(xff_ens3(1).gt.0.)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))
                         if(xff_ens3(2).gt.0.)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne))
                         if(xff_ens3(3).gt.0.)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne))
                         if(xff_ens3(13).gt.0.)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4))
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5))
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6))
                         xf_ens(i,j,nall+14)=max(0.,xff_ens3(14))
                         a1=max(1.e-2,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9)/a1)
                         a1=max(1.e-2,pr_ens(i,j,nall+15))
                         xf_ens(i,j,nall+15)=max(0.,xff_ens3(15)/a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
                            xf_ens(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))
                            xf_ens(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))
                            xf_ens(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne))
                         endif
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! 16 is a randon pick from the other 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
                endif
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!!!!    NOT USED FOR "NORMAL" APPLICATION (maxens=1)

                if(maxens.gt.1)then
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif
                endif

                   endif
 350            continue
                if(maxens.gt.1)then
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                   endif
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
             enddo
          endif
 100   continue

!mjo    Make sure any individual closures are not going crazy
        xf_ens(:,:,:) = min(xf_ens(:,:,:), 2.*xmbmax)

   END SUBROUTINE cup_forcing_ens_3d


   SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,                         &
              xl,cp,ztexec,zqexec,use_excess,       &
              ktf,kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        use_excess,ktf,kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ztexec,zqexec,cap_max,cap_inc
     real,intent (in   )                  ::                           &
        xl,cp
     real,    dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        hkb
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2
     real                                 ::                           &
        pbcdif,plus,hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)+1
      if(iloop.eq.5)KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
         endif
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=hkb(i) ! HE_cup(I,K22(I))
      if(iloop.eq.5)then
       hetest=HKB(I)
!      do k=1,k22(i)
!        hetest=max(hetest,he_cup(i,k))
!      enddo
      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))then
!       write(0,*)'htest',k22(i),kbcon(i),HETEST,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        GO TO 31
      endif

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).eq.0)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=25.
      if(iloop.eq.5.and.cap_max(i).gt.25.)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
!       write(0,*)'htest',k22(i),kbcon(i),plus,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)+1
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=max(1,min(kbcon(i)-1,k22(i)+1))  !kbcon(i)-1
             k2=k22(i)+1
             hkb(i)=sum(he_cup(i,k1:k2))/float(k2-k1+1)+(xl*zqexec(i)+cp*ztexec(i))/float(k2-k1+1)
        else if(use_excess <= 1)then
             hkb(i)=he_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
        endif  ! excess

        if(iloop.eq.5)KBCON(I)=K22(I)
        IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
         endif
        GO TO 27
      ENDIF
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ierrc,ilo,dby,kbcon,ktop,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
          if(ilo.eq.1)ierrc(i)="problem with defining ktop"
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
            ierrc(i)="kbcon == ktop "
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf,kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              ktf,kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf,kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0



   SUBROUTINE cup_output_ens_3d(xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,ierr2,ierr3,pr_ens,             &
              sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training,  &
	      ipr,jpr,ktf,kts,kte)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ipr,jpr,ktf,kts,kte
     integer, intent (in   )              ::                           &
        j,training
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens
!srf ------
!    real,    dimension (its:ite,jts:jte)                              &
     real,    dimension (its:ite,jts:jte)                              &
         ,intent (inout)                   ::                          &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI 
     real, dimension( its:ite , jts:jte )                      &
         ,intent(in) :: weight_gr,weight_w,weight_mc,weight_st,weight_as
!-srf---
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
     real,   dimension (its:ite)                                      &
         ,intent (in  )                   ::                           &
        sig
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite,kts:kte,1:maxens2)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,prerate,clos_wei,xmbhelp
     real                                 ::                           &
        dtts,dtqs
     real,    dimension (its:ite)         ::                           &
       xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm,wm1,wm2,wm3
     real, dimension (its:ite,5) :: xmb_w

!
      character *(*), intent (in)        ::                           &
       name

!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        subt(i,k)=0.
        subq(i,k)=0.
        sub_mas(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=0.
         xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*maxens*maxens2*maxens3+1,iens*maxens*maxens2*maxens3
           if(pr_ens(i,j,n).le.0.)then
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'pr_ens',n,pr_ens(i,j,n),xf_ens(i,j,n)
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
       xmb_w=0.
!
!-- now do feedback
!
      ddtes=100.
      do i=its,itf
        if(ierr(i).eq.0)then
         k=0
         xmb_ave(i)=0.
         do n=(iens-1)*maxens*maxens2*maxens3+1,iens*maxens*maxens2*maxens3
          k=k+1
          xmb_ave(i)=xmb_ave(i)+xf_ens(i,j,n)
         enddo
         xmb_ave(i)=xmb_ave(i)/float(k)
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
         xmb(i)=sig(i)*xmb_ave(i)
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))

           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
!mjo limit mass flux
           xmb(i) = min(xmb(i),xmbmax)
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
      ENDDO
      DO k=kts,ktf
      do i=its,itf
            dtt =0.
            dtts=0.
            dtq =0.
            dtqs=0.
            dtqc=0.
            dtpw=0.
        IF(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,maxens2
              dtt =dtt  + dellat  (i,k,n)
              dtts=dtts + subt_ens(i,k,n)
              dtq =dtq  + dellaq  (i,k,n)
              dtqs=dtqs + subq_ens(i,k,n)
              dtqc=dtqc + dellaqc (i,k,n)
              dtpw=dtpw + pw      (i,k,n)
           enddo
           OUTTEM(I,K)= XMB(I)* dtt /float(maxens2)
           SUBT  (I,K)= XMB(I)* dtts/float(maxens2)
           OUTQ  (I,K)= XMB(I)* dtq /float(maxens2)
           SUBQ  (I,K)= XMB(I)* dtqs/float(maxens2)
           OUTQC (I,K)= XMB(I)* dtqc/float(maxens2)
	   PRE(I)=PRE(I)+XMB(I)*dtpw/float(maxens2)
           sub_mas(i,k)=zu(i,k)*xmb(i)
!          xf_ens(i,j,:)=sig(i)*xf_ens(i,j,:)*dtpw/float(maxens2)
        endif
       enddo
      enddo

!!      do i=its,itf
!!        if(ierr(i).eq.0)then
!!        do k=(iens-1)*maxens*maxens2*maxens3+1,iens*maxens*maxens2*maxens3
!!          xf_ens(i,j,k)=sig(i)*xf_ens(i,j,k)*xfac1(i)
!!        enddo
!!        endif
!!      ENDDO

!srf-fix for preci
      do i=its,itf
        if(ierr(i).ne. 0)then
            apr_w (i,j)=0.0
	    apr_st(i,j)=0.0
	    apr_gr(i,j)=0.0
	    apr_mc(i,j)=0.0
	    apr_as(i,j)=0.0
        endif
      ENDDO
!srf
   END SUBROUTINE cup_output_ens_3d
!-------------------------------------------------------
   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              ccnclean,p_cup,kbcon,ktop,cd,dby,clw_all,&
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,         &
              ZQEXEC,use_excess,ccn,rho, &
              up_massentr,up_massdetr,psum,psumh,                 &
              autoconv,aeroevap,itest,ktf,j,ipr,jpr,                &
              kts,kte                     )

   IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (maritime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  use_excess,itest,autoconv,aeroevap,ktf,           &
                                  j,ipr,jpr, kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
     real,    dimension (its:ite)                              &
        ,intent (in   )                   ::                           &
        zqexec
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite,kts:kte) ::                           &
        qch,qrcb,pwh,clw_allh
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
!
!  local variables in this routine
!

     integer                              ::                           &
        iounit,iprop,iall,i,k,k1,k2
     real                                 ::                           &
        prop_ave,qrcb_h,bdsp,dp,g,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real,    dimension (kts:kte)         ::                           &
        prop_b
!
        prop_b(kts:kte)=0
        iall=0
        c0=.002
        g=9.81
        bdsp=BDISPM
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')c0=0.
        do i=its,itf
          pwav(i)=0.
          pwavh(i)=0.
          psum(i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          pwh(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
        enddo
        enddo
      if(use_excess < 2 ) then
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=2,kbcon(i)-1
        DZ=Z_cup(i,K)-Z_cup(i,K-1)
        qc(i,k)=qe_cup(i,k22(i))+float(use_excess)*zqexec(i)
        qch(i,k)=qe_cup(i,k22(i))+float(use_excess)*zqexec(i)
        if(qc(i,k).gt.qes_cup(i,kbcon(i)-1))then
            pw(i,k)=zu(i,k)*(qc(i,k)-qes_cup(i,kbcon(i)-1))
            qc(i,k)=qes_cup(i,kbcon(i)-1)
            qch(i,k)=qes_cup(i,kbcon(i)-1)
            PWAV(I)=PWAV(I)+PW(I,K)
            Psum(I)=Psum(I)+pw(I,K)*dz
        endif
      enddo
      endif
      enddo
      else if(use_excess == 2) then
        do i=its,itf
         if(ierr(i).eq.0)then
             k1=max(1,k22(i)-1)
             k2=k22(i)+1
          do k=2,kbcon(i)-1
             DZ=Z_cup(i,K)-Z_cup(i,K-1)
             qc (i,k)=sum(qe_cup(i,k1:k2))/float(k2-k1+1) +zqexec(i)
             qch(i,k)=sum(qe_cup(i,k1:k2))/float(k2-k1+1) +zqexec(i)
             if(qc(i,k).gt.qes_cup(i,kbcon(i)-1))then
                 pw(i,k)=zu(i,k)*(qc(i,k)-qes_cup(i,kbcon(i)-1))
                 qc(i,k)=qes_cup(i,kbcon(i)-1)
                 qch(i,k)=qes_cup(i,kbcon(i)-1)
                 PWAV(I)=PWAV(I)+PW(I,K)
                 Psum(I)=Psum(I)+pw(I,K)*dz
             endif
          enddo !k
         endif  !ierr
        enddo !i
      endif  ! use_excess

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Lt.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         rhoc=.5*(rho(i,k)+rho(i,k-1))
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
         DP=p_cup(i,K)-p_cup(i,K-1)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
       qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
       qch(i,k)= (qch(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*qch(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

        if(qc(i,k).le.qrch)qc(i,k)=qrch
        if(qch(i,k).le.qrch)qch(i,k)=qrch
!
!------- Total condensed water before rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
        clw_allh(i,k)=QCH(I,K)-QRCH
        QRCB(I,K)=(QCH(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
    IF(autoconv.eq.2) then


! 
! normalized berry
!
! first calculate for average conditions, used in cup_dd_edt!
! this will also determine proportionality constant prop_b, which, if applied,
! would give the same results as c0 under these conditions
!
         q1=1.e3*rhoc*qrcb(i,k)  ! g/m^3 ! g[h2o]/cm^3
!        berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCNclean/ &
!               ( q1 * BDSP)  ) ) !/(
!mjo avoid problems with small q1
         q1 = max(q1,1.e-3)
         berryc0 = q1*q1*q1 / (q1*300.0 + 2.1960 * CCNclean/BDSP)
!mjo
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',k,rhoc,rho(i,k)
!         qrcb_h=qrcb(i,k)/(1.+c0*dz)
         qrcb_h=((QCH(I,K)-QRCH)*zu(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
         prop_b(k)=c0*qrcb_h*zu(i,k)/(1.e-3*berryc0)
         pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
         berryc=qrcb(i,k)
         qrcb(i,k)=((QCh(I,K)-QRCH)*zu(i,k)-pwh(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
!        QRCb(I,K) = qrcb(i,k) - pwh(i,k)
         if(qrcb(i,k).lt.0.)then
           berryc0=(qrcb(i,k-1)*(.5*up_massdetr(i,k-1))-(QCh(I,K)-QRCH)*zu(i,k))/zu(i,k)*1.e-3*dz*prop_b(k)
           pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k)
           qrcb(i,k)=0.
         endif
!     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',zu(i,k),pwh(i,k),dz,qrch,qrcb(i,k),clw_allh(i,k)
      QCh(I,K)=QRCb(I,K)+qrch
      PWAVH(I)=PWAVH(I)+pwh(I,K)
      Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
!
! then the real berry
!
          q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
!         berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
!               ( q1 * BDSP)  ) ) !/(
!mjo avoid problems with small q1
          q1 = max(q1,1.e-3)
          berryc0 = q1*q1*q1 / (q1*300.0 + 2.1960 * CCN(i)/BDSP)
!mjo
          berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
          berryc=qrc(i,k)
          qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-zu(i,k)*berryc0-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
          if(qrc(i,k).lt.0.)then
            berryc0=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/zu(i,k)
            qrc(i,k)=0.
          endif
          pw(i,k)=berryc0*zu(i,k)
          QC(I,K)=QRC(I,K)+qrch
!
!  if not running with berry at all, do the following
!
       ELSE       !c0=.002
         qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
         if(qrc(i,k).lt.0.)then
           qrc(i,k)=0.
           pw(i,k)=0.
         endif
!
!
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
        QC(I,K)=QRC(I,K)+qrch
      endif !autoconv
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
         Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz
 100     CONTINUE
       prop_ave=0.
       iprop=0
       do k=kts,kte
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0.)iprop=iprop+1
       enddo
       iprop=max(iprop,1)
!      write(11,*)'prop_ave = ',prop_ave/float(iprop)
!      print *,'pwav = ',pwav(1)

   END SUBROUTINE cup_up_moisture

!--------------------------------------------------------------------

      real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
     &        ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
     &    log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
        satvap = 10 ** eilog
      else
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
     &             (log(tsot) / log(10.))
        ewlog2 = ewlog - 1.3816e-07 * &
     &             (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + .0081328 * &
     &             (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
        satvap = 10 ** ewlog4
      end if
      return
      end function

   SUBROUTINE CUP_gf_sh(xmb_out,zo,OUTQC,J,AAEQ,T,Q,Z1,                    &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,PSUR,US,VS,    &
              TCRIT,                                        &
              ztexec,zqexec,ccn,ccnclean,rho,dx,dhdt,                               &
              kpbl,kbcon,ktop,cupclws,k22,         &   !-lxz
              xland,gsw,tscl_kf,              &
              xl,rv,cp,g,ichoice,ipr,jpr,ierr,ierrc,         &
              autoconv,ktf,use_excess,kts,kte                                &
                                                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,ktf,ktau,use_excess,kts,kte,ipr,jpr
     integer, intent (in   )              ::                           &
        j,ichoice
  !
  ! 
  !
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        cupclws,OUTT,OUTQ,OUTQC
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb_out
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop,k22
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,T,PO,P,US,VS,tn,dhdt
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ztexec,zqexec,ccn,Z1,PSUR,AAEQ,xland
       
       real                                                            &
        ,intent (in   )                   ::                           &
        tscl_kf,dx,ccnclean,dtime,tcrit,xl,cp,rv,g


!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq,subt,subq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       xmbmax,XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,zuhe
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,entr_rate,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd,fsum
      
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: xx1,xx2,xxx,power_entr,zustart,zufinal,dzm1,dzp1


     integer :: tun_lim,jprnt,k1,k2,kbegzu,kfinalzu,kstart,jmini,levadj
     logical :: keep_going
     real xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (kts:kte) :: smth
      zustart=.1
      zufinal=1.
      levadj=4
      power_entr=2.
      day=86400.
      do i=its,itf
        xmb_out(i)=0.
        xland1(i)=1.
        if(xland(i).gt.1.5)xland1(i)=0.
        cap_max_increment(i)=25.
        ierrc(i)=" "
        pre(i)=0.
      enddo
!
!--- initial entrainment rate (these may be changed later on in the
!--- program
!
      entr_rate =.2/200.
      tun_lim=7
      
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        up_massentro(i,k)=0.
        up_massdetro(i,k)=0.
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        qrco(i,k)=0.
        cd(i,k)=1.*entr_rate
        dellaqc(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!
!--- minimum depth (m), clouds must have
!
      depth_min=50.
!
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
!
      cap_maxs=25.
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
      enddo
      do i=its,itf
          cap_max(i)=cap_maxs
        iresult=0
      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           ktf,kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           ktf,kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
      kbmax(i)=min(kbmax(i),ktf-4)
      endif
      enddo

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22,ierr, &
           ktf,kts,kte)
       DO 36 i=its,itf
         if(kpbl(i).gt.5)cap_max(i)=po_cup(i,kpbl(i))
         IF(ierr(I).eq.0)THEN
         IF(K22(I).GT.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
         endif
            if(kpbl(i).gt.5)then
               k22(i)=kpbl(i)
               ierr(i)=0
               ierrc(i)="reset to zero becausof kpbl"
             endif
         else
             ierrc(i)="why here? "
         endif
!      if(j.eq.jpr .and. i.eq.ipr)write(0,*)'initial k22 = ',k22(ipr),kpbl(i)
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0)THEN
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)+1
             hkb(i) =sum(he_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             qkbo(i)=sum(qo_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)
!            write(0,*)sum(heo_cup(i,k1:k2))/float(k2-k1+1),heo_cup(i,k1),heo(i,k1:k2)
        else if(use_excess <= 1) then
             hkb(i)=he_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             hkbo(i)=heo_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             qkbo(i)=qo_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i))
        endif  ! excess
         do k=1,k22(i)
            hkb(i)=max(hkb(i),he_cup(i,k))
            hkbo(i)=max(hkbo(i),heo_cup(i,k))
            qkbo(i)=max(qkbo(i),qo_cup(i,k))
         enddo
       endif ! ierr
      enddo
      call cup_kbcon(ierrc,cap_max_increment,5,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess,       &
           ktf,kts,kte)
!
!--- increase detrainment in stable layers
!
      DO 887 i=its,itf
         IF(ierr(I).eq.0)THEN
            if(kbcon(i).ge.ktf-tun_lim)then
                ierr(i)=231
                go to 887
            endif
            do k=kts,ktf
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate*(1.3-frh)
               cd(i,k)=entr_rate_2d(i,k)
            enddo
            zuhe(i)=zustart
            kstart=1
            xx1=z_cup(i,kbcon(i))*z_cup(i,kbcon(i))
            xx2=z_cup(i,kstart)*z_cup(i,kstart)
            frh=(zufinal-zustart)/(xx1-xx2)
            dh=zuhe(i)-frh*xx2
            do k=kstart,kbcon(i)-1
             dz=z_cup(i,k+1)-z_cup(i,k)
             cd(i,k)=0.
             xxx=z_cup(i,k+1)*z_cup(i,k+1)
             entr_rate_2d(i,k)=((frh*xxx+dh)/zuhe(i)-1.+cd(i,k)*dz)/dz
             zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'entr = ',k,entr_rate_2d(i,k),dh,frh,zuhe(i),dz
            enddo
            xx1=z_cup(i,kbcon(i)+tun_lim)*z_cup(i,kbcon(i)+tun_lim)
            xx2=z_cup(i,kbcon(i)-1)*z_cup(i,kbcon(i)-1)
            frh=-(0.1-zuhe(i))/(xx1-xx2)
            dh=zuhe(i)+frh*z_cup(i,kbcon(i))*z_cup(i,kbcon(i))
               do k=kbcon(i),kbcon(i)+tun_lim
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 xxx=z_cup(i,k+1)*z_cup(i,k+1)
                 cd(i,k)=-((-frh*xxx+dh)/zuhe(i)-1.-entr_rate_2d(i,k)*dz)/dz
                 zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)-cd(i,k)*dz*zuhe(i)
!            if(i.eq.ipr.and.j.eq.jpr)write(0,*)'entr = ',k,entr_rate_2d(i,k),cd(i,k),zuhe(i)
               enddo
               do k=kbcon(i)+tun_lim+1,ktf
                entr_rate_2d(i,k)=0.
                cd(i,k)=0.
               enddo

        ENDIF
 887  enddo
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
         hco(i,k)=0.
         DBYo(I,K)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
            qco(i,k)=qkbo(i)
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         qco(i,k)=qkbo(i)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
         trash=QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
              /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)
         qrco(i,k)=max(0.,qco(i,k)-trash)
       endif ! ierr
      enddo
!
!
      do 42 i=its,itf
         if(ierr(i).eq.0)then
         zu(i,1)=zustart
         zuo(i,1)=zustart
!    mass entrainment and detrinament is defined on model levels
         do k=2,ktf-1 !kbcon(i)+4 ! ktf-1
          dz=zo_cup(i,k)-zo_cup(i,k-1)
          up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
          up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
          zuo(i,k)=zuo(i,k-1)+up_massentro(i,k-1)-up_massdetro(i,k-1)
          if(zuo(i,k).lt.0.05)then
             zuo(i,k)=.05
             up_massdetro(i,k-1)=zuo(i,k-1)-.05  + up_massentro(i,k-1)
             cd(i,k-1)=up_massdetro(i,k-1)/dz/zuo(i,k-1)
          endif
          zu(i,k)=zuo(i,k)
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
!          zu(i,k)=max(0.01,zu(i,k-1)+up_massentr(i,k-1)-up_massdetr(i,k-1))
         enddo
         do k=kbcon(i)+1,ktf-1
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo
         do k=kbcon(i)+1,ktf
          if(dbyo(i,k).lt.0.)then
              ktop(i)=k-1
              go to 41
          endif
         enddo
41       continue
         if(ktop(i).lt.kbcon(i)+1)then
            ierr(i)=5
            ierrc(i)='ktop is less than kbcon+1'
             go to 42
         endif
         if(ktop(i).gt.ktf-2)then
             ierr(i)=5
             ierrc(i)="ktop is larger than ktf-2"
             go to 42
         endif
         do k=kbcon(i)+1,ktop(i)
          trash=QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
              /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)
          qco(i,k)=   (qco(i,k-1)*zuo(i,k-1)-.5*up_massdetr(i,k-1)* qco(i,k-1)+ &
                         up_massentr(i,k-1)*qo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          qrco(i,k)=max(0.,qco(i,k)-trash)
          cupclws(i,k)=qrco(i,k)*.5
         enddo
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           zuo(i,k)=0.
           qco(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
!        if(i.eq.ipr.and.j.eq.jpr)then
!           write(0,*)'hcnew = '
!           do k=1,ktf
!             write(0,*)k,hco(i,k),dbyo(i,k)
!           enddo
!        endif
      endif
42    continue
!     enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellah(i,k)=0.
        dsubt(i,k)=0.
        dsubh(i,k)=0.
        dellaq(i,k)=0.
        dsubq(i,k)=0.
      enddo
      enddo
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
             dsubt(i,1)=0. 
             dsubq(i,1)=0. 
         do k=kts+1,ktop(i)
               subin=0.
               subdown=0.
! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
!
!         SPECIAL LEVELS
!
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               entup=0.
               detup=0.
            endif
            totmas=subin-subdown+detup-entup  &
             -entupk+detupk+zuo(i,k+1)-zuo(i,k)
!               print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               print *,jmin(i),k22(i),kbcon(i),ktop(i)
               write(0,123)k,subin,subdown,detup,entup, &
                           entupk,detupk,zuo(i,k+1),zuo(i,k)
123     formAT(1X,i2,10E12.4)
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k)=(detup*.5*(HCo(i,K+1)+HCo(i,K)) &
                    -entup*heo(i,k) &
                    +subin*heo_cup(i,k+1) &
                    -subdown*heo_cup(i,k) &
                    +detupk*(hco(i,ktop(i))-heo_cup(i,ktop(i)))    &
                    -entupk*heo_cup(i,k22(i)) &
                     )*g/dp
            dellaq(i,k)=(detup*.5*(qco(i,K+1)+qco(i,K)-qrco(i,k+1)-qrco(i,k)) &
                    -entup*qo(i,k) &
                    +subin*qo_cup(i,k+1) &
                    -subdown*qo_cup(i,k) &
                    +detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))    &
                    -entupk*qo_cup(i,k22(i)) &
                     )*g/dp
          
!
! updraft subsidence only
!
           if(k.lt.ktop(i))then
             dsubt(i,k)=(zuo(i,k+1)*heo_cup(i,k+1) &
                    -zuo(i,k)*heo_cup(i,k))*g/dp
             dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1) &
                    -zuo(i,k)*qo_cup(i,k))*g/dp
!          if(i.eq.ipr.and.j.eq.jpr)then
!           write(0,*)'dq3',k,zuo(i,k+1)*heo_cup(i,k+1),zuo(i,k)*heo_cup(i,k)
!          endif
           endif
!
       enddo   ! k

        endif
      enddo
!
!-- take out cloud liquid water for detrainment
!
      do k=kts,ktf-1
      do i=its,itf
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k)-po_cup(i,k+1))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*up_massdetro(i,k)*.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k)-po_cup(i,k+1))
         endif
         if(dellaqc(i,k).lt.0.)write(0,*)'neg della',i,j,k,ktop(i),qrco(i,k), &
              qrco(i,k+1),up_massdetro(i,k),zuo(i,ktop(i))
         dellaqc(i,k)=max(0.,dellaqc(i,k))
       endif
      enddo
      enddo
!
!--- using dellas, calculate changed environmental profiles
!
      mbdt=3.e-4

      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            dsubh(i,k)=dsubt(i,k)
            dellaq(i,k)=dellaq(i,k)+dellaqc(i,k)
            dellaqc(i,k)=0.
            XHE(I,K)=(dsubt(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ(I,K)=(dsubq(i,k)+DELLAQ(I,K))*MBDT+QO(I,K)
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            dSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k))*MBDT+TN(I,K)
!           IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
            IF(XQ(I,K).LE.1.E-08)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      xhkb(i)=hkbo(i)+(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
!     IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      IF(XQ(I,ktf).LE.1.E-08)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           ktf,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           ktf,kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
!     do i=its,itf
!       if(ierr(i).eq.0)then
!         xhkb(i)=xhe(i,k22(i))
!       endif
!     enddo
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
!        if(use_excess == 2) then
!            k1=max(1,k22(i)-1)
!            k2=k22(i)+1
!            xhkb(i) =sum(xhe_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
!        else if(use_excess <= 1) then
!            xhkb(i)=xhe_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
!        endif

         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
      xzu(i,:)=zuo(i,:)
      do k=kbcon(i)+1,ktop(i)
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           ktf,kts,kte)
!
! now for shallow forcing
!
       do i=its,itf
        xmb(i)=0.
        xff_shal(1:9)=0.
        if(ierr(i).eq.0)then
          xmbmax(i)=0.1  
          xkshal=(xaa0(i)-aa1(i))/mbdt
          if(xkshal.ge.0.)xkshal=+1.e6
          if(xkshal.gt.-1.e-4 .and. xkshal.lt.0.)xkshal=-1.e-4
          xff_shal(1)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
          xff_shal(1)=min(xmbmax(i),xff_shal(1))
          xff_shal(2)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
          xff_shal(2)=min(xmbmax(i),xff_shal(2))
          xff_shal(3)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
          xff_shal(3)=min(xmbmax(i),xff_shal(3))
          if(aa1(i).le.0.)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
          if(aa1(i)-aa0(i).le.0.)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
! boundary layer QE (from Saulo Freitas)
          blqe=0.
          trash=0.
          if(k22(i).lt.kpbl(i)+1)then
             do k=1,kbcon(i)-1
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash=max((hc(i,kbcon(i))-he_cup(i,kbcon(i))),1.e1)
             xff_shal(7)=max(0.,blqe/trash)
             xff_shal(7)=min(xmbmax(i),xff_shal(7))
          else
             xff_shal(7)=0.
          endif
          if(xkshal.lt.-1.1e-04)then ! .and.  &
!            ((aa1(i)-aa0(i).gt.0.) .or. (xff_shal(7).gt.0)))then
          xff_shal(4)=max(0.,-aa0(i)/(xkshal*tscl_KF))
          xff_shal(4)=min(xmbmax(i),xff_shal(4))
          xff_shal(5)=xff_shal(4)
          xff_shal(6)=xff_shal(4)
          else
           xff_shal(4)=0.
           xff_shal(5)=0.
           xff_shal(6)=0.
          endif
!         write(0,888)'i0=',i,j,kpbl(i),blqe,xff_shal(7)
!888       format(a3,3(1x,i3),2e12.4)
          xff_shal(8)= xff_shal(7)
          xff_shal(9)= xff_shal(7)
          fsum=0.
          do k=1,9
           xmb(i)=xmb(i)+xff_shal(k)
           fsum=fsum+1.
          enddo
          xmb(i)=min(xmbmax(i),xmb(i)/fsum)
!         if(i.eq.ipr.and.j.eq.jpr)write(0,*)',ierr,xffs',ierr(i),xff_shal(1:9),xmb(i),xmbmax(i)
          if(xmb(i).eq.0.)ierr(i)=22
          if(xmb(i).eq.0.)ierrc(i)="22"
          if(xmb(i).lt.0.)then
             ierr(i)=21
             ierrc(i)="21"
!            write(0,*)'neg xmb,i,j,xmb for shallow = ',i,j,k22(i),ierr(i)
          endif
        endif
        if(ierr(i).ne.0)then
           k22(i)=0
           kbcon(i)=0
           ktop(i)=0
           xmb(i)=0
           do k=kts,ktf
              outt(i,k)=0.
              outq(i,k)=0.
              outqc(i,k)=0.
           enddo
        else if(ierr(i).eq.0)then
!
! got the mass flux, sanity check, first for heating rates
!
          trash=0.
!         kmaxx=0
          do k=2,ktop(i)
           trash=max(trash,86400.*(dsubt(i,k)+dellat(i,k))*xmb(i))
          enddo
          if(trash.gt.100.)then
             xmb(i)=xmb(i)*100./trash
          endif
          trash=0.
          do k=2,ktop(i)
           trash=min(trash,86400.*(dsubt(i,k)+dellat(i,k))*xmb(i))
          enddo
          if(trash.lt.-100.)then
              xmb(i)=-xmb(i)*100./trash
          endif
!
! sanity check on moisture tendencies: do not allow anything that may allow neg
! tendencies
!
          do k=2,ktop(i)
           trash=q(i,k)+(dsubq(i,k)+dellaq(i,k))*xmb(i)*dtime
          if(trash.lt.1.e-12)then
! max allowable tendency over tendency that would lead to too small mix ratios
!
            trash=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
            xmb(i)=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
          endif
          enddo
          xmb_out(i)=xmb(i)
! 
! final tendencies
!
          do k=2,ktop(i)
           outt(i,k)=(dsubt(i,k)+dellat(i,k))*xmb(i)
           outq(i,k)=(dsubq(i,k)+dellaq(i,k))*xmb(i)
          enddo
        endif
       enddo
!      
! done shallow
!--------------------------done------------------------------
!
   END SUBROUTINE CUP_gf_sh

END MODULE module_cu_gf
