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
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use consts_coms, only: cp, xl=>alvl, rv=>rvap, g=>grav
  implicit none

  !-- ichoice_deep:  0 ensemble, 1 grell, 4 omega , 7 MCONV, 10 = KF, 13 = PB
  integer, parameter :: ichoice = 0

  !-- ichoice_shallow:  0 ensemble, 1 Wstar, 4 heat-engine or 7 BLQE 
  integer, parameter :: ichoice_s = 0

  integer, parameter :: iens    =  1
  integer, parameter :: maxens  =  1
  integer, parameter :: maxens2 =  1
  integer, parameter :: maxens3 = 16
  integer, parameter :: ensdim  = maxens * maxens2 * maxens3

  integer, parameter :: its = 1, ite = 1, itf = 1
  integer, parameter :: jts = 1, jte = 1, jtf = 1
  integer, parameter :: kts = 1

  integer, parameter ::  ens4 = 1
  real,    parameter :: fens4 = real(ens4)

  integer, parameter :: j = 1
  integer, parameter :: imid = 0
  real,    parameter :: ccnclean = 250.
  integer, parameter :: autoconv = 2 ! 1=Kessler, 2=Berry with ccn
  integer, parameter :: aeroevap = 3 ! 1=orig, 2=aerosol dep, 3=orig+aerosol

  logical, parameter :: masscon = .false.  !- new mass conserv form flag
  logical, parameter :: dicycle = .false.  !- diurnal cycle flag
  logical, parameter :: entrnew = .false.  !- new entr formulation

  integer, parameter :: use_excess = 0    ! Include PBL temp/humidity excess
                                          ! for shallow Cu
  real,    parameter :: tcrit = 253.16

  private
  public :: gf_driver

CONTAINS

  subroutine gf_driver(iw, dtlong)

    use mem_turb,    only: kpblh, frac_land, fthpbl, fqtpbl, wtv0, &
                           ustar, wstar, sfluxt, sfluxr
    use consts_coms, only: p00i, eradi, grav, gravi, alvl, alvlocp, vonk, r8
    use mem_radiate, only: fthrd_sw, fthrd_lw
    use mem_grid,    only: mza, lpw, arw0, zm, zt, xew, yew, zew, &
                           dzt, arw, lpv, arv, volt
    use mem_ijtabs,  only: itab_w
    use mem_micro,   only: cldnum
    use mem_basic,   only: wmc, vmc, theta, tair, press, rho, sh_v, &
                           vxe, vye, vze
    use oname_coms,  only: nl
    use mem_cuparm,  only: thsrc, rtsrc, conprr, kcutop, kcubot, cbmf, &
                           qwcon, vxsrc, vysrc, vzsrc
    implicit none

    integer, intent(in)  :: iw
    real,    intent(in)  :: dtlong

    real, parameter :: onethird = 1./3.

    integer :: ktf                 ! number of vertical levels
    real    :: dx                  ! grid spacing
    real    :: dtime               ! time step
    integer :: kpbl   (1)          ! layer corresponding to PBL height
    real    :: ccn    (1)          ! ccn concentration ( #/cm^3 ?)
    real    :: zo     (1,mza)      ! height
    real    :: r      (1,mza)      ! air density
    real    :: aaeq   (1)          ! turn off convection if negative
    real    :: t      (1,mza)      ! input temperature
    real    :: q      (1,mza)      ! input mixing ratio
    real    :: z1     (1)          ! surface elevation
    real    :: tn     (1,mza)      ! input forced temperature
    real    :: qo     (1,mza)      ! input forced mixing ratio
    real    :: po     (1,mza)      ! input forced pressure
    real    :: p      (1,mza)      ! input pressure
    real    :: psur   (1)          ! input surface pressure
    real    :: us     (1,mza)      ! input zonal wind
    real    :: vs     (1,mza)      ! input meridional wind
    real    :: zws    (1)          ! wstar
    real    :: dhdt   (1,mza)      ! moist static energy tendency
    integer :: ierr   (1)          ! error code
    real    :: xmb    (1)          ! deep conv convective mass flux
    real    :: sub_mas(1,mza)
    real    :: subt   (1,mza)      ! temp tendency due to subsidence
    real    :: subq   (1,mza)      ! water vapor tendency due to subsidence
    real    :: pre    (1)          ! output precipitation rate
    real    :: outt   (1,mza)      ! output temp tendency
    real    :: outq   (1,mza)      ! output water vapor tendency
    real    :: outqc  (1,mza)      ! output cloud water tendency
    real    :: outu   (1,mza)      ! output conv u wind tendency
    real    :: outv   (1,mza)      ! output conv v wind tendency
    real    :: cupclw (1,mza)      ! cloud water
    integer :: kbcon  (1)          ! cloud base (LCL)
    integer :: ktop   (1)          ! cloud top
    integer :: k22    (1)          ! updraft originating level
    real    :: xf_ens (1,1,ensdim) ! mass flux ensembles
    real    :: pr_ens (1,1,ensdim) ! precipitation ensembles
    real    :: xland  (1)          ! 1 = land, 0 = water
    real    :: mconv  (1,1)        ! integrated column moisture convergence
    real    :: omeg   (1,mza,ens4) ! vertical velocity in pressure coordinates
    real    :: ztexec (1)          ! PBL temperature excess
    real    :: zqexec (1)          ! PBL humidity excess

    character(50) :: ierrc(1)  ! error messages

    integer :: k, ka, kc, npoly, iwn, jv, iv
    real    :: raxis, raxisi, dirv, flx, ws, uvtr
    real    :: exner(mza)
    real    :: vflux(mza), vflux_the(mza), vflux_vap(mza)
    real    :: hflux, hflux_the, hflux_vap
    real    :: fqvadv, fthadv
    real(r8):: qsum, qav, tsum, tav

    ka    = lpw(iw)
    npoly = itab_w(iw)%npoly
    kpbl  = kpblh(iw) - ka + 1

    ktf  = mza + 1 - ka  ! number of ATM levels to process

    ! Go no higher than 50mb for convective calculations to prevent any
    ! problems when esat gets near ambient pressure
    do k = ka, mza-1
       if (press(k,iw) < 50.e2) then
          ktf = k - ka + 1
          exit
       endif
    enddo

    dx    = sqrt(arw0(iw))
    dtime = dtlong
    psur  = 0.01 * (press(ka,iw) + (zt(ka)-zm(ka-1))*rho(ka,iw)*grav)
    z1    = zm(ka-1)
    ccn   = cldnum(iw) * 1.e-6
    zws   = ( ustar(iw)**3 + wstar(iw)**3 )**onethird

    !if (allocated(frac_land)) then
    !   xland = frac_land(iw)
    !else
    !   xland = 1.0
    !endif
    xland = 1.0

    if (wtv0(iw) > 0.0 .and. use_excess>=1) then
       ws = vonk * grav * wtv0(iw) * (zt(ka)-zm(ka-1)) / theta(ka,iw)
       ws = ( ustar(iw)**3 + ws )**onethird
       ws = max(ws,0.1)
       ztexec(1) = max( sfluxt(iw) / ws, 0.0)
       zqexec(1) = max( sfluxr(iw) / ws, 0.0)
    else
       ztexec(1) = 0.0
       zqexec(1) = 0.0
    endif

    raxis  = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis
    raxisi = 1.0 / max(raxis, 1.e-12)

    ! Vertical advective theta and water vapor fluxes (W levels)

    do k = ka, mza-1
       vflux(k) = arw(k,iw) * wmc(k,iw)

       if (wmc(k,iw) >= 0.0) then
          vflux_vap(k) = vflux(k) * sh_v(k,iw)
          vflux_the(k) = vflux(k) * theta(k,iw)
       else
          vflux_vap(k) = vflux(k) * sh_v(k+1,iw)
          vflux_the(k) = vflux(k) * theta(k+1,iw)
       endif
    enddo

    vflux    (ka-1)  = 0.
    vflux    (mza) = 0.
    vflux_the(ka-1)  = 0.
    vflux_the(mza) = 0.
    vflux_vap(ka-1)  = 0.
    vflux_vap(mza) = 0.

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

       ! Only pbl forcing changes moist static energy
       dhdt(1,kc) = cp * exner(k) * fthpbl(k,iw) + alvl * fqtpbl(k,iw)

       ! U and V winds
       if (raxis > 1.e3) then
          us(1,kc) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi

          vs(1,kc) = vze(k,iw) * raxis * eradi  &
               - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi
       else
          us(1,kc) = vxe(k,iw)
          vs(1,kc) = vye(k,iw)
       endif

       omeg(1,kc,1) = -grav * wmc(k,iw)

       ! moisture convergence
       mconv(1,1) = mconv(1,1) + omeg(1,kc,1) * (sh_v(k+1,iw) - sh_v(k,iw)) * gravi

    enddo

    mconv(1,1) = max(mconv(1,1), 0.0)

    aaeq = 1.
    xmb = 0.0
    sub_mas = 0.0
    xf_ens = 0.0
    pr_ens = 0.0
    k22 = 0
    kbcon = 0
    ktop = 0
    cupclw = 0.0

    cbmf(iw) = 0.0

    outt (1,:) = 0.0
    outq (1,:) = 0.0
    outqc(1,:) = 0.0
    outu (1,:) = 0.0
    outv (1,:) = 0.0
    subt (1,:) = 0.0
    subq (1,:) = 0.0
    pre  (1)   = 0.0

    call cup_gf(ktf, dx, dtime, kpbl, ccn, r, mconv, omeg, aaeq, t, q,  &
                z1, xland, tn, qo, zo, po, p, psur, us, vs, zws, dhdt,  &
                ierr, ierrc, xmb, sub_mas, subt, subq, pre, outt, outq, &
                outqc, outu, outv, cupclw, kbcon, ktop, k22, xf_ens, pr_ens )

    if (pre(1) > 1.e-16 .and. kbcon(1) > 0 .and. ktop(1) >= kbcon(1)) then
        
       ! Deep convecton is active. Copy tendencies to model arrays.

       kcutop(iw) = ktop (1) + ka - 1
       kcubot(iw) = kbcon(1) + ka - 1

       cbmf(iw) = xmb(1)

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

       tsum = cp * sum( abs(outt(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
       tav  = cp * sum(     outt(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) ) &
            - pre(1) * alvl
 
       qav = qav / max(qsum, 1.e-20_r8)
       tav = tav / max(tsum, 1.e-20_r8)
        
       do kc = 1, ktf
          k  = kc + ka - 1
          rtsrc(k,iw) =  outq(1,kc) - qav * abs(outq(1,kc))
          thsrc(k,iw) = (outt(1,kc) - tav * abs(outt(1,kc))) / exner(k)
          qwcon(k,iw) = cupclw(1,kc)
       enddo

       conprr(iw) = pre(1)

       ! convective momentum transport

       if (nl%conv_uv_mix > 0) then

          if (raxis > 1.e3) then
             do kc = 1, ktf
                k  = kc + ka - 1
                uvtr = -outv(1,kc) * zew(iw) * eradi
                vxsrc(k,iw) = (-outu(1,kc) * yew(iw) + uvtr * xew(iw)) * raxisi
                vysrc(k,iw) = ( outu(1,kc) * xew(iw) + uvtr * yew(iw)) * raxisi
                vzsrc(k,iw) =   outv(1,kc) * raxis * eradi 
             enddo
          else
             do kc = 1, ktf
                k  = kc + ka - 1
                vxsrc(k,iw) = outu(1,kc)
                vysrc(k,iw) = outv(1,kc)
                vzsrc(k,iw) = 0.0
             enddo
          endif

       endif

    else

       ! If there is no deep convection, check for shallow convection

       outt (1,:) = 0.0
       outq (1,:) = 0.0
       outqc(1,:) = 0.0
       pre  (1)   = 0.0

       cupclw = 0.
       xmb = 0.
       k22 = 0
       kbcon = 0
       ktop = 0
       ierr = 0
       ierrc = ""
       aaeq = 1

       call CUP_gf_sh(ktf, dtime, psur, z1, &
                      ztexec, zqexec, dhdt, zws, zo, aaeq, &
                      t, q, tn, qo, po, p, kpbl, &
                      ierr, ierrc, xmb, outt, outq, outqc, &
                      cupclw, kbcon, ktop, k22)

       if (kbcon(1) > 0 .and. ktop(1) >= kbcon(1)) then

          ! Shallow convection is active; copy tendencies

          kcutop(iw) = ktop (1) + ka - 1
          kcubot(iw) = kbcon(1) + ka - 1

          cbmf(iw) = xmb(1)

          ! Slightly modify tendencies to ensure heat and moisture conservation

          qav  = sum(     outq(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
          qsum = sum( abs(outq(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )

          tav  = cp * sum(     outt(1,1:ktf)  * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )
          tsum = cp * sum( abs(outt(1,1:ktf)) * rho(ka:ktf+ka-1,iw) * dzt(ka:ktf+ka-1) )

          qav = qav / max(qsum, 1.e-20_r8)
          tav = tav / max(tsum, 1.e-20_r8)

          do kc = 1, ktf
             k  = kc + ka - 1
             rtsrc(k,iw) =  outq(1,kc) - qav * abs(outq(1,kc))
             thsrc(k,iw) = (outt(1,kc) - tav * abs(outt(1,kc))) / exner(k)
             qwcon(k,iw) = cupclw(1,kc)
          enddo

       endif

       ! Although deep convection mixes momentum, shallow convecton
       ! module does not. Compute momentum mixing here.

       if ( nl%conv_uv_mix > 0 ) then
          call acmcld_uvmix(iw, dtlong)
       endif

    endif

  end subroutine gf_driver

!========================================================================================

   SUBROUTINE CUP_gf(ktf	       &
     		    !input data
     		    ,dx 	       &
     		    ,DTIME	       &
     		    ,kpbl	       &
     		    ,ccn	       &
     		    ,rho	       &
     		    ,mconv	       &
     		    ,omeg	       &
     		    ,AAEQ	       &
     		    ,T  	       &
     		    ,Q  	       &
     		    ,Z1                &
                    ,xland             &
     		    ,TN 	       &
     		    ,QO 	       &
     		    ,zo 	       &
     		    ,PO 	       &
     		    ,P  	       &
     		    ,PSUR	       &
     		    ,US 	       &
     		    ,VS 	       &
     		    ,zws	       &
     		    ,dhdt	       &
     		    !output data
     		    ,ierr,ierrc        &
     		    ,xmb_out	       &
     		    ,sub_mas ,subt,subq&
     		    ,PRE	       &
     		    ,OUTT	       &
     		    ,OUTQ	       &
     		    ,OUTQC	       &
     		    ,outu	       &
     		    ,outv	       &
     		    ,cupclw	       &
     		    ,kbcon,ktop,k22    &
     		    ,xf_ens,pr_ens)
  
     IMPLICIT NONE
     		    
     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  !
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        xf_ens,pr_ens
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (inout  )                   ::                           &
        outu,outv,OUTT,OUTQ,OUTQC,subt,subq,sub_mas,cupclw
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        dhdt,rho,T,PO,P,US,VS,tn
     real,    dimension (its:ite,kts:ktf,1:ens4)                       &
        ,intent (inout   )                   ::                           &
        omeg
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ccn,Z1,PSUR,xland
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        AAEQ
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        zws
     real, dimension (its:ite,1:ens4)                                         &
        ,intent (in   )                   ::                           &
        mconv


       real  ,intent (in   )                   ::                          &
        dx,dtime
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
     real,    dimension (its:ite,kts:ktf,1:maxens2) ::                 &
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
  ! ichoice     = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:ktf) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,hcot,                                                   &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq,          &
        u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv

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
       cap_max_increment,closure_n,psum,psumh,sig,sigd,zuhe
     real,    dimension (its:ite,1:ens4) ::                                   &
        axx
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult,nvar
     real                                 ::                           &
      day,dz,dzo,mbdt,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd,radiusd,frhd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas,beta,entr_rate
      real :: zusafe,xxx,xx1,xx2,zutop,zustart,zufinal,dzm1,dzp1


     integer :: k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini,levadj
     logical :: keep_going
     real tot_time_hr,blqe
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:ktf) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real dts,fp,fpi,up_massent,up_massdet
     real :: xff_mid(its:ite,2)

      integer :: iversion=1
      real :: umean,T_star
      real, dimension (its:ite)         :: aa1_bl,hkbo_bl,tau_bl,tau_ecmwf,wmean
      real, dimension (its:ite,kts:ktf) :: tn_bl, qo_bl, qeso_bl, heo_bl, heso_bl &
                                          ,qeso_cup_bl,qo_cup_bl, heo_cup_bl,heso_cup_bl&
                                          ,gammao_cup_bl,tn_cup_bl,hco_bl,DBYo_bl
      real, dimension(its:ite) :: xf_dicycle
      real :: C_up, E_dn,G_rain,trash2

      real, parameter :: zustart_dp = 0.3
      real, parameter :: zufinal_dp = 1.0
      real, parameter :: zutop_dp   = 0.2
     
!srf- end
!
!proportionality constant to estimate pressure gradient of updraft (Zhang and Wu, 2003, JAS
!
!- ecmwf formulation
     real, parameter :: lambau=2.
     real, parameter :: pgcon=0.

!    if(imid.eq.1)then
!- SAS formulation
!     real, parameter :: lambau=0.
!     real, parameter :: pgcon=-.55
!    endif

      zustart=zustart_dp
      zufinal=zufinal_dp
      zutop = zutop_dp
      if(imid.eq.1)zutop=.1

      levadj=5
      cap_maxs=150.

      do i=its,itf
        edto(i) = 0.0
        closure_n(i)=16.
        xland1(i)=xland(i) ! 1.
        xmb_out(i)=0.
        cap_max(i)=cap_maxs
        cap_max_increment(i)=20.
        if(imid.eq.1)cap_max_increment(i)=10.
        if(xland(i).gt.1.5 .or. xland(i).lt.0.5)then
            xland1(i)=0.
            if(imid.eq.0)cap_max(i)=cap_maxs-50.
            if(imid.eq.1)cap_max(i)=cap_maxs-50.
        endif
        ierrc(i)=" "
!       cap_max_increment(i)=1.
      enddo
!
!--- initial entrainment rate (these may be changed later on in the

      if(ENTRNEW) then
         entr_rate  = 1.e-3
         mentrd_rate= 0.3*entr_rate
      else
         entr_rate  = 1.e-4 !-21/10/2015 7.0e-5
         mentrd_rate= 2.0*entr_rate
      endif

      radius =.2/entr_rate
      radiusd=.2/mentrd_rate
      frh =3.14*radius*radius/dx/dx
      frhd=3.14*radiusd*radiusd/dx/dx
      if(frh .gt. 0.55)then !srf orig: 0.7
         frh=.55            !srf orig: 0.7
         radius=sqrt(frh*dx*dx/3.14)
         entr_rate=.2/radius
      endif
      if(frhd .gt. 0.7)then
         frhd=.7
         radiusd=sqrt(frhd*dx*dx/3.14)
         mentrd_rate=.2/radiusd
      endif
      do i=its,itf
         sig(i)=(1.-frh)**2
         sigd(i)=1.
!        sigd(i)=sig(i)/(1.-frhd)**2
!        if(imid.eq.0)sig(i)=sig(i)
      enddo
!
!--- entrainment of mass
!
      xlamdd=mentrd_rate
!
!--- initial detrainmentrates
!
      do k=kts,ktf
       do i=its,itf
        hcot(i,k)=0.
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
      if(imid.eq.1)edtmax=.3
      edtmin=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=1000.
      if(imid.eq.1)depth_min=500.
!
!--- maximum depth (mb) of capping
!--- inversion (larger cap = no convection)
!
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
        edt(i)=0.
        kstabm(i)=ktf-1
        IERR2(i)=0
        IERR3(i)=0
      enddo

      iresult=0
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
      if(imid.eq.1)zkbmax=3000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1000.
!     if(imid.eq.1)z_detr=800.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      mbdt_ens(1)=10. !dtime*1.E-02
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
      call cup_env(z,qes,he,hes,t,q,p,z1, psur,ierr,-1,ktf)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1,psur,ierr,-1,ktf)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, ierr,z1,ktf)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,ierr,z1,ktf)
      do i=its,itf
        if(ierr(i).eq.0)then
          u_cup(i,kts)=us(i,kts)
          v_cup(i,kts)=vs(i,kts)
          do k=kts+1,ktf
           u_cup(i,k)=.5*(us(i,k-1)+us(i,k))
           v_cup(i,k)=.5*(vs(i,k-1)+vs(i,k))
          enddo
        endif
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!     if(ierr(i).eq.0)then
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
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      !-srf to increase are coverage change "2" below to "1", the start point
      CALL cup_MAXIMI(HEO_CUP,2,KBMAX,K22,ierr,ktf)
       DO 36 i=its,itf
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GE.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
           ktop(i)=0
           k22(i)=0
           kbcon(i)=0
         endif
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0.)THEN
          hkb(i)=he_cup(i,k22(i))
          hkbo(i)=heo_cup(i,k22(i))
       endif ! ierr
      enddo
      iloop=1
!     if(imid.eq.1)iloop=5
      call cup_kbcon(ierrc,cap_max_increment,iloop,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max,ktf,z_cup,entr_rate,heo)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,ktf)
      DO i=its,itf
         IF(ierr(I).eq.0.)THEN
            kdefi=0
            do k=kts,ktf

               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
                   !
               !-------------------------------------------
               if(ENTRNEW) then
               !- v 2
                  if(k >= kbcon(i)) then
                     entr_rate_2d(i,k)=entr_rate*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,kbcon(i)))**3
                  else
                     entr_rate_2d(i,k)=entr_rate*(1.3-frh)
                  endif
                  cd(i,k)=0.75e-4*(1.6-frh)        
           
               ELSE        
                  !- v 1
                  entr_rate_2d(i,k)=entr_rate*(1.3-frh)
                  cd(i,k)=1.*entr_rate
        
               ENDIF
               !
            enddo
          endif
       enddo
!
! get entrainment and detrainmentrates for updraft
!
     call rates_up_pdf(ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
                    kstabi,k22,kbcon,kts,ktf,zuo,kpbl,deep=.true.)
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         uc(i,k)=0.
         vc(i,k)=0.
         hc(i,k)=0.
         dby (i,k)=0.
         hco (i,k)=0.
         dbyo(i,k)=0.
         up_massentro(i,k)=0.
         up_massdetro(i,k)=0.
         up_massentr (i,k)=0.
         up_massdetr (i,k)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
            uc(i,k)=u_cup(i,k22(i))
            vc(i,k)=v_cup(i,k22(i))
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         uc(i,k)=u_cup(i,k22(i))
         vc(i,k)=v_cup(i,k22(i))

         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
       endif ! ierr
      enddo
!
!
      do i=its,itf
         if(ierr(i).eq.0)then
        
         do k=1,ktop(i)+1
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
         enddo

         !- mass entrainment and detrinament is defined on model levels

         do k=2,maxloc(zuo(i,:),1)
         !=> below location of maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)

           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(up_massentro(i,k-1).lt.0.)then
              up_massentro(i,k-1)=0.
              up_massdetro(i,k-1)=zuo(i,k-1)-zuo(i,k)
              if(zuo(i,k-1).gt.0.) &
                cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
           endif
           if(zuo(i,k-1).gt.0.) &
             entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktop(i)
         !=> above location of maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
           if(up_massdetro(i,k-1).lt.0.)then
              up_massdetro(i,k-1)=0.
              up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)
              if(zuo(i,k-1).gt.0.) &
                entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
           endif

           if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo

         do k=2,ktf-1
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo

!---mass con
!        do k=kbcon(i)+1,ktop(i)  ! original
         do k=k22(i)  +1,ktop(i)  ! mass cons option
!---mass con
        
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*uc(i,k-1) &
                         -lambau*up_massdetr(i,k-1)*uc(i,k-1) +         &
                         (up_massentr(i,k-1)+lambau*up_massdetr(i,k-1))*us(i,k-1) &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1)))  /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*vc(i,k-1)                &
                         -lambau*up_massdetr(i,k-1)*vc(i,k-1) +                        &
                         (up_massentr(i,k-1)+lambau*up_massdetr(i,k-1))*vs(i,k-1)      &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /    &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo

         if(ktop(i).lt.kbcon(i)+2)then
            ierr(i)=5
            ierrc(i)='ktop too small'
         endif
         !
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           UC(i,K)=u_cup(i,k)
           VC(i,K)=v_cup(i,k)
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
!     call cup_minimi(HEso_cup,K22,kstabi,JMIN,ierr, &
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr,ktf)
      DO 100 i=its,ite
         IF(ierr(I).eq.0.)THEN
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
         IF(ierr(I).eq.0.)THEN
            if ( jmin(i) - 1 .lt. kdet(i)   ) kdet(i) = jmin(i)-1
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
       ucd(i,k)=u_cup(i,k)
       vcd(i,k)=v_cup(i,k)
       dbydo(i,k)=0.
      enddo
      enddo
      do i=its,itf

          beta=.05
          if(imid.eq.1)beta=.0
          bud(i)=0.
          IF(ierr(I).eq.0)then

        !- this calls routine to get downdrafts normalized mass flux
        !
        zutop  = 0.2 ! zd at jmin
        zustart= beta  ! zd at surface
        zufinal= 1.  ! zd at jmin-lev_adjust
        mentrd_rate_2d(i,:)=mentrd_rate
        cdd(i,1:jmin(i))=xlamdd
        cdd(i,jmin(i))=0.
        dd_massdetro(i,:)=0.
        dd_massentro(i,:)=0.

        call get_zu_zd_pdf(ierr(i),kdet(i),jmin(i),zdo(i,:),kts,ktf,kpbl(i),draft=3)

        xzd(i,jmin(i))= zdo(i,jmin(i))
        zd (i,jmin(i))= zdo(i,jmin(i))
!        write(92,*)'k,zdo(i,k)'
!       do k=jmin(i),1,-1
!        write(92,*)k,zdo(i,k)
!       enddo

!-srf-21jul2015 mass cons
       !do ki=jmin(i)-1,maxloc(zdo(i,:),1),-1
        do ki=jmin(i)  ,maxloc(zdo(i,:),1),-1
!-srf mass cons
          !=> from jmin to maximum value zd -> change entrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
          dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
          if(dd_massentro(i,ki).lt.0.)then
             dd_massentro(i,ki)=0.
             dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
             if(zdo(i,ki+1) > 0.0)&
               cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1) > 0.0)&
            mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
        enddo
        mentrd_rate_2d(i,1)=0.
        do ki=maxloc(zdo(i,:),1)-1,1,-1
          !=> from maximum value zd to surface -> change detrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massentro(i,ki)=mentrd_rate_2d(i,ki)*dzo*zdo(i,ki+1)
          dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
          if(dd_massdetro(i,ki).lt.0.)then
            dd_massdetro(i,ki)=0.
            dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
            if(zdo(i,ki+1) > 0.0)&
              mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1) > 0.0)&
            cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
        enddo

!         write(92,*)'k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)'
        do k=jmin(i),1,-1
          xzd(i,k)= zdo(i,k)
          zd (i,k)= zdo(i,k)
          dd_massentr(i,k)=dd_massentro(i,k)
          dd_massdetr(i,k)=dd_massdetro(i,k)
!         write(92,*)k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)
        enddo

! downdraft moist static energy + moisture budget
            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))
!-srf-21jul2015 mass cons
           !do ki=jmin(i)-1,1,-1
            do ki=jmin(i)  ,1,-1
!-srf mass cons
             dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
             ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*ucd(i,ki+1)+ &
                        dd_massentro(i,ki)*us(i,ki)          &
                        -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /     &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*vcd(i,ki+1)+ &
                        dd_massentro(i,ki)*vs(i,ki)            &
                        -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+ &
                        dd_massentro(i,ki)*heo(i,ki))   /            &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
             bud(i)=bud(i)+dbydo(i,ki)*dzo
            enddo
          endif

        if(bud(i).gt.0)then
          ierr(i)=7
          ierrc(i)='downdraft is not negatively buoyant '
        endif
      enddo
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_new(ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,ktf)
!
!--- calculate moisture properties of updraft
!
      if(imid.eq.0)then
      call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo, &
           p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,        &
           ccn,rho,up_massentr,up_massdetr,psum,psumh,1,ktf,deep=.true.)
      else if(imid.eq.1)then
      call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo, &
           p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,        &
           ccn,rho,up_massentr,up_massdetr,psum,psumh,1,ktf,deep=.false.)
      endif
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=kts,ktop(i)
         cupclw(i,k)=qrco(i,k)        ! my mod
      enddo
      endif
      enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,ktf)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,ktf)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
           aaeq(i)=aa1(i)
         endif
      enddo

!=================================================================
!-srf- begin
      !--- AA1 from boundary layer (bl) processes only
      aa1_bl    (:) = 0.0
      xf_dicycle(:) = 0.0
      !- way to calculate the fraction of cape consumed by shallow convection
      iversion=1 ! ecmwf
      !iversion=0 ! orig
      !
      ! Bechtold et al 2008 time-scale of cape removal
      DO i=its,itf
            if(ierr(i).eq.0)then
                !- mean vertical velocity
                wmean(i) = 7.0 ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
                !- time-scale cape removal from  Betchold et al. 2008
                tau_ecmwf(i)=( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                tau_ecmwf(i)= tau_ecmwf(i) * (1.0061 + 1.23E-2 * (dx/1000.))! dx must be in meters
            endif
      enddo
      !
      IF(dicycle) then

        DO i=its,itf

            if(ierr(i).eq.0)then
                if(xland(i) > 0.9 ) then
                  !- over water
                  umean= 2.0+sqrt(2.0*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
                  tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
                else
                  !- over land
                  tau_bl(i) =( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                endif

            endif
        ENDDO

        if(iversion == 1) then
        !-- version ecmwf

           !- T_star = temp scale in original paper = 1 K
           ! T_star = 1.0
           T_star = 3.0

           !-- calculate pcape from BL forcing only
            call cup_up_aa1bl(aa1_bl,t,tn,q,qo,dtime, &
                              zo_cup,zuo,kbcon,ktop,ierr,ktf)

            DO i=its,itf
                  
                  if(ierr(i).eq.0)then

                      !- only for convection rooting in the PBL
                      !if(zo_cup(i,kbcon(i))-z1(i) > 500.0) then !- instead 500 -> zo_cup(kpbl(i))
                      !        aa1_bl(i) = 0.0
                      !else
                      !        !- multiply aa1_bl/T_* by the " time-scale" - tau_bl
                         aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i)
                      !endif
                      !print*,'aa0,aa1bl=',aa0(i),aa1_bl(i),aa0(i)-aa1_bl(i),tau_bl(i)!,dtime,xland(i)        
                  endif
            ENDDO

        else
        
          !- version for real cloud-work function
        
          !-get the profiles modified only by bl tendencies
          DO i=its,itf
           tn_bl(i,:)=0.;qo_bl(i,:)=0.
           if ( ierr(i) == 0 )then
            !below kbcon -> modify profiles
            tn_bl(i,1:kbcon(i)) = tn(i,1:kbcon(i))
            qo_bl(i,1:kbcon(i)) = qo(i,1:kbcon(i))
                 !above kbcon -> keep environment profiles
            tn_bl(i,kbcon(i)+1:ktf) = t(i,kbcon(i)+1:ktf)
            qo_bl(i,kbcon(i)+1:ktf) = q(i,kbcon(i)+1:ktf)
           endif
          ENDDO
          !--- calculate moist static energy, heights, qes, ... only by bl tendencies
          call cup_env(zo,qeso_bl,heo_bl,heso_bl,tn_bl,qo_bl,po,z1,psur,ierr,-1,ktf)
          !--- environmental values on cloud levels only by bl tendencies
          call cup_env_clev(tn_bl,qeso_bl,qo_bl,heo_bl,heso_bl,zo,po,qeso_cup_bl,qo_cup_bl, &
                            heo_cup_bl,heso_cup_bl,zo_cup,po_cup,gammao_cup_bl,tn_cup_bl, &
                            psur,ierr,z1,ktf)
          DO i=its,itf
            IF(ierr(I).eq.0.)THEN
               hkbo_bl(i)=heo_cup_bl(i,k22(i))
            endif ! ierr
          ENDDO
          DO k=kts,ktf
           do i=its,itf
             hco_bl (i,k)=0.
             DBYo_bl(i,k)=0.
           enddo
          ENDDO
          DO i=its,itf
            IF(ierr(I).eq.0.)THEN
             do k=1,kbcon(i)-1
              hco_bl(i,k)=hkbo_bl(i)
             enddo
             k=kbcon(i)
             hco_bl (i,k)=hkbo_bl(i)
             DBYo_bl(i,k)=Hkbo_bl(i) - HESo_cup_bl(i,k)
            ENDIF
          ENDDO
!        
!        
          DO i=its,itf
            if(ierr(i).eq.0)then
               do k=kbcon(i)+1,ktop(i)
                    hco_bl(i,k)=(hco_bl(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_bl(i,k-1)+ &
                               up_massentro(i,k-1)*heo_bl(i,k-1))   /               &
                               (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                    dbyo_bl(i,k)=hco_bl(i,k)-heso_cup_bl(i,k)
               enddo
               do k=ktop(i)+1,ktf
                  hco_bl (i,k)=heso_cup_bl(i,k)
                  dbyo_bl(i,k)=0.0
               enddo
            endif
          ENDDO

          !--- calculate workfunctions for updrafts
          call cup_up_aa0(aa1_bl,zo,zuo,dbyo_bl,GAMMAo_CUP_bl,tn_cup_bl, &
                        kbcon,ktop,ierr,ktf)

          DO i=its,itf

            if(ierr(i).eq.0)then
                !- get the increment on AA0 due the BL processes
                aa1_bl(i) = aa1_bl(i) - aa0(i)
                !- only for convection rooting in the PBL
                !if(zo_cup(i,kbcon(i))-z1(i) > 500.0) then !- instead 500 -> zo_cup(kpbl(i))
                !   aa1_bl(i) = 0.0
                !else
                !   !- multiply aa1_bl the "normalized time-scale" - tau_bl/ model_timestep
                   aa1_bl(i) = aa1_bl(i)* tau_bl(i)/ dtime
                !endif
                !print*,'aa0,aa1bl=',aa0(i),aa1_bl(i),aa0(i)-aa1_bl(i),tau_bl(i)!,dtime,xland(i)
            endif
           ENDDO
        ENDIF
     ENDIF  ! version of implementation
!-srf -end
!=================================================================


       do i=1,ens4
       axx(:,i)=aa1(:)
       enddo

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwo,ccn,pwevo,edtmax,edtmin,edtc,psum,psumh, &
           rho,ktf)
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=sigd(i)*edtc(i,iedt)
         edto(i)=sigd(i)*edtc(i,iedt)
         edtx(i)=sigd(i)*edtc(i,iedt)
         if(maxens2.eq.3)then
            edt(i)=sigd(i)*edtc(i,maxens2)
            edto(i)=sigd(i)*edtc(i,maxens2)
            edtx(i)=sigd(i)*edtc(i,maxens2)
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
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellu(i,k)=0.
        dellv(i,k)=0.
        dellah(i,k)=0.
        dsubt(i,k)=0.
        dsubh(i,k)=0.
        dellat(i,k)=0.
        dellaq(i,k)=0.
        dellaqc(i,k)=0.
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
!----------------------------------------------  cloud level 3  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if(.not. masscon) then
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!-------------------------------------------------- form orig
      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellu(i,1)=(edto(i)*zdo(i,2)*ucd(i,2)   &
                     -edto(i)*zdo(i,2)*u_cup(i,2))*g/dp
         dellv(i,1)=(edto(i)*zdo(i,2)*vcd(i,2)   &
                     -edto(i)*zdo(i,2)*v_cup(i,2))*g/dp
         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp
         dellaq(i,1)=(edto(i)*zdo(i,2)*qrcdo(i,2)   &
                     -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
         dsubh(i,1)=0.
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
! updraft originates at k22, only updraft term at k22-1 is (zu-entupk)*he
!           if(k.eq.k22(i)-1)then
!              entupk=zuo(i,k+1)
!           endif
! downdraft originating level, similiar to k22-1 for updraft
!-srf- 21jul2015 mass con            
            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !endif
!-srf- mass con            
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
!              print *,jmin(i),k22(i),kbcon(i),ktop(i)
               write(0,123)k,subin,subdown,detup,entup, &
                           detdo,entdo,entupk,detupk
123     formAT(1X,i2,8E12.4)
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
            dellu(i,k)=((detup+lambau*detup)*.5*(uC(i,K+1)+uC(i,K)) &
                    +detdo*.5*(UCD(i,K+1)+UCD(i,K)) &
                    -(entup+lambau*detup)*us(i,k) &
                    -entdo*us(i,k) &
                    +subin*u_cup(i,k+1) &
                    -subdown*u_cup(i,k) &
                    +detupk*(uc(i,ktop(i))-u_cup(i,ktop(i)))    &
                    -entupk*u_cup(i,k22(i)) &
                    -entdoj*u_cup(i,jmin(i)) &
                     )*g/dp
            dellv(i,k)=((detup+lambau*detup)*.5*(vC(i,K+1)+vC(i,K)) &
                    +detdo*.5*(vCD(i,K+1)+vCD(i,K)) &
                    -(entup+lambau*detup)*vs(i,k) &
                    -entdo*vs(i,k) &
                    +subin*v_cup(i,k+1) &
                    -subdown*v_cup(i,k) &
                    +detupk*(vc(i,ktop(i))-v_cup(i,ktop(i)))    &
                    -entupk*v_cup(i,k22(i)) &
                    -entdoj*v_cup(i,jmin(i)) &
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
             dsubh(i,k)=(zuo(i,k+1)*heo_cup(i,k+1) &
                    -zuo(i,k)*heo_cup(i,k))*g/dp
             dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1) &
                    -zuo(i,k)*qo_cup(i,k))*g/dp
             dellu(i,k)=dellu(i,k)+(zuo(i,k+1)*u_cup(i,k+1) &
                    -zuo(i,k)*u_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)))*g/dp
             dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
                    -zuo(i,k)*v_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
!            dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
!                   -zuo(i,k)*v_cup(i,k)   &
!                   +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
!          elseif (k.eq.ktop(i))then
!            dsubt(i,k)=(-zuo(i,k)*heo_cup(i,k))*g/dp
!            dsubq(i,k)=(-zuo(i,k)*qo_cup(i,k))*g/dp
           endif
!
! in igh res case, subsidence terms are for meighbouring points only. This has to be
! done mass consistent with the della term
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
         if(dellaqc(i,k).lt.0)write(0,*)'neg della',i,j,k,ktop(i),qrco(i,k), &
              qrco(i,k+1),up_massdetro(i,k),zuo(i,ktop(i))
         dellaqc(i,k)=max(0.,dellaqc(i,k))
       endif
      enddo
      enddo
      dellah(:,:)=dellah(:,:)+dsubh(:,:)
      dellaq(:,:)=dellaq(:,:)+dsubq(:,:)
      dsubh(:,:)=0.0
      dsubq(:,:)=0.0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
else

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellu(i,1)=(edto(i)*zdo(i,2)*ucd(i,2)   &
                     -edto(i)*zdo(i,2)*u_cup(i,2))*g/dp
         dellv(i,1)=(edto(i)*zdo(i,2)*vcd(i,2)   &
                     -edto(i)*zdo(i,2)*v_cup(i,2))*g/dp

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
            ! downdraft originating level, similiar to k22-1 for updraft
!-srf- 21jul2015 mass con            
!            if(k.eq.jmin(i))then
!               entdoj=edto(i)*zdo(i,k)
!            endif
!-srf- mass con
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
!              print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               !print *,jmin(i),k22(i),kbcon(i),ktop(i)
               write(0,123)k,subin,subdown,detup,entup, &
                           detdo,entdo,entupk,detupk
               !123     formAT(1X,i2,8E12.4)
               ! call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellu(i,k)=((detup+lambau*detup)*.5*(uC(i,K+1)+uC(i,K)) &
                    +detdo*.5*(UCD(i,K+1)+UCD(i,K)) &
                    -(entup+lambau*detup)*us(i,k) &
                    -entdo*us(i,k) &
                    +subin*u_cup(i,k+1) &
                    -subdown*u_cup(i,k) &
                    +detupk*(uc(i,ktop(i))-u_cup(i,ktop(i)))    &
                    -entupk*u_cup(i,k22(i)) &
                    -entdoj*u_cup(i,jmin(i)) &
                     )*g/dp
            dellv(i,k)=((detup+lambau*detup)*.5*(vC(i,K+1)+vC(i,K)) &
                    +detdo*.5*(vCD(i,K+1)+vCD(i,K)) &
                    -(entup+lambau*detup)*vs(i,k) &
                    -entdo*vs(i,k) &
                    +subin*v_cup(i,k+1) &
                    -subdown*v_cup(i,k) &
                    +detupk*(vc(i,ktop(i))-v_cup(i,ktop(i)))    &
                    -entupk*v_cup(i,k22(i)) &
                    -entdoj*v_cup(i,jmin(i)) &
                     )*g/dp


!
! updraft subsidence only
!
           if(k.lt.ktop(i))then

             dellu(i,k)=dellu(i,k)+(zuo(i,k+1)*u_cup(i,k+1) &
                    -zuo(i,k)*u_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)))*g/dp
             dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
                    -zuo(i,k)*v_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
           endif
       enddo   ! k

        endif
      enddo


      do i=its,itf
        trash  = 0.0
        trash2 = 0.0
        if(ierr(i).eq.0)then
         dsubq(i,:)  =0.0
         dsubh(i,:)  =0.0

         dp=100.*(po_cup(i,1)-po_cup(i,2))

         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp

         dellaqc(i,1)=0.0
         dellaq (i,1)=(edto(i)*zdo(i,2)*qcdo(i,2)   &
                      -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
        
         G_rain=  0.5*(pwo (i,1)+pwo (i,2))*g/dp
         E_dn  = -0.5*(pwdo(i,1)+pwdo(i,2))*g/dp*edto(i)  ! pwdo < 0 and E_dn must > 0
         dellaq(i,1) = dellaq(i,1)+ E_dn-G_rain

         !--- conservation check
         !- water mass balance
         trash = trash  + (dellaq(i,1)+dellaqc(i,1)+G_rain-E_dn)*dp/g         
         !- H  budget
         trash2 = trash2+ (dellah(i,1))*dp/g
        
         !write(0,*) "delH1=",1,dellah(i,k),dellaq(i,k),E_dn
         !write(3,*)'=>H k22 kbcon ktop= ',k22(i),kbcon(i),ktop(i)
         !write(3,*)'=>H= ',1,real(trash2,4),real(dellah(i,1),4)

         !write(4,*)'=>W k22 kbcon ktop= ',k22(i),kbcon(i),ktop(i)
         !write(4,*)'=>W= ',1,real(trash,4),real(dellaq(i,1),4)

         do k=kts+1,ktop(i)
            ! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
            ! detrainment and entrainment for downdrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
            ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            ! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
            ! downdraft originating level, similiar to k22-1 for updraft

            !write(0,*)"down",k,edto(i),detdo,entdo,subin,subdown
        
            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !endif
            !if(k.eq.ktop(i))then
            !   detupk=zuo(i,ktop(i))
            !   subin=0.
            !   subdown=0.
            !   detdo=0.
            !   entdo=0.
            !   entup=0.
            !  detup=0.
            !endif
            !totmas=subin-subdown+detup-entup-entdo+ &
            !       detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)
                        

            !- check H conservation
             trash2 = trash2+ (dellah(i,k))*dp/g
        
        
            !===tmp
            !entdoj=0.0
            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !   !write(0,*)"entdoj=",zdo(i,k),zdo(i,k+1),dd_massentro(i,k),k!;stop 43
            !endif
            !===tmp

            !-- take out cloud liquid water for detrainment
            detup=up_massdetro(i,k)
            dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            
            !--10jul latest form of Georg - to be tested (see cup_up_moisture)
            !dellaqc(i,k) = zuo(i,k)*c1*qrco(i,k)*dz/dp*g !detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            !if(k.eq.ktop(i))dellaqc(i,k)= detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            !
            
            !---
            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
            !
            !write(0,*) "eva=",k,pwdo(i,k),E_dn,zdo(i,k  )
            !
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid water (qrco) + converted to rain
        
            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain
        
            !-- water vapor budget
            !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
            !--   - condensation term + evaporation
            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                         - C_up + E_dn
            !- check water conservation liq+condensed (including rainfall)
             trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

         !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
         !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
        
         enddo   ! k
         !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
         !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
         !  write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
         !  !stop 33
         !endif

        endif

      enddo
endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(1)
      do i=its,itf
      xaa0_ens(i,:)=0.
      enddo

      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            !adding dellaqc to dellaq:
            !dellaq (i,k)= dellaq(i,k)+dellaqc(i,k)
            !dellaqc(i,k)=0.0
            !        
            XHE(I,K)=(dsubh(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ (I,K)=(dsubq(i,k)+DELLAQ(I,K)+DELLAQC(i,k))*MBDT+QO(I,K)

            !- don't feed dellat with dellaqc if
            !- the detrainment of liquid water will be used as
            !- a source for cloud microphysics
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            dSUBT(I,K)=(1./cp)*(dsubh(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k)-xl/cp*dellaqc(i,k))*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,ktf)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,ierr,z1,ktf)
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
!        do k=1,k22(i)
!           xhc(i,k)=xhe_cup(i,k)
!        enddo
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
!ss     xzu(i,:)=zuo(i,:)
!-- mass con
!     do k=kbcon(i)+1,ktop(i)  ! orig
      do k=k22(i)  +1,ktop(i)  ! mass cons option
!-- mass con
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
           kbcon,ktop,ierr,ktf)
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
      IF(maxens.gt.0)then
       call cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr,ktf)
       call cup_kbcon(ierrc,cap_max_increment,2,k22x,kbconx,heo_cup, &
            heso_cup,hkbo,ierr2,kbmax,po_cup,cap_max, &
            ktf,z_cup,entr_rate,heo)
       call cup_kbcon(ierrc,cap_max_increment,3,k22x,kbconx,heo_cup, &
            heso_cup,hkbo,ierr3,kbmax,po_cup,cap_max, &
            ktf,z_cup,entr_rate,heo)
      ENDIF
!
!--- calculate cloud base mass flux
!
      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,'deeps',axx,                 &
           iedt,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           ktf, &
           tau_ecmwf,aa1_bl,xf_dicycle)
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           subt_ens(i,k,iedt)=dsubt(i,k)
           subq_ens(i,k,iedt)=dsubq(i,k)
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edto(i)*pwdo(i,k)
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
       IF(imid.eq.1)THEN
         do i=its,itf
          xff_mid(i,1)=0.
          xff_mid(i,2)=0.
          if(ierr(i).eq.0)then
          blqe=0.
          trash=0.
          if(k22(i).lt.kpbl(i)+1)then
             do k=1,kpbl(i)
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash=max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
             xff_mid(i,1)=max(0.,blqe/trash)
             xff_mid(i,1)=min(0.1,xff_mid(i,1))
          endif
          xff_mid(i,2)=.03*zws(i)
          endif
         enddo
       ENDIF
       call cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat_ens,dellaq_ens, &
            dellaqc_ens,subt_ens,subq_ens,subt,subq,outt,     &
            outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
            'deep',ierr2,ierr3,pr_ens,                        &
            sig,closure_n,ktf,xf_dicycle )
      k=1
      do i=its,itf
          if(ierr(i).eq.0) then
             PRE(I)=MAX(PRE(I),0.)
             xmb_out(i)=xmb(i)
             do k=kts,ktop(i)
               outu(i,k)=dellu(i,k)*xmb(i)
               outv(i,k)=dellv(i,k)*xmb(i)
               sub_mas(i,k)=zuo(i,k)*xmb(i)
             enddo
          else
             ktop(i)=0
          endif
      enddo

   END SUBROUTINE CUP_gf

!==============================================================================

   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,edtc,psum2,psumh, &
              rho,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:ktf)                              &
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
        edtmax,edtmin
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

     real :: pefc, aeroadd, rhoc, prop_c

     real, parameter :: alpha3 =  1.9
     real, parameter :: beta3  = -1.13
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
!              prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1.)) !*1.e6
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
               if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif


!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*pwav(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,q,he,t_cup,iloop,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ktf
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he 
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:ktf)                              &
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
      if(dh.lt.0)then
        QRCD(I,K)=(qes_cup(i,k)+(1./XL)*(GAMMA_cup(i,k) &
                  /(1.+GAMMA_cup(i,k)))*DH)
      else
          qrcd(i,k)=qes_cup(i,k)
      endif
      pwd(i,jmin(i))=zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
      qcd(i,k)=qrcd(i,k)
      pwev(i)=pwev(i)+pwd(i,jmin(i)) ! *dz
!
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
!        QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
!                 +entr*DZ*q(i,Ki) &
!                )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
!        dz=qcd(i,ki)
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
         pwev(i)=pwev(i)+pwd(i,ki) ! *dz
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.0.and.iloop.eq.1)then
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
              psur,ierr,itest,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
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
       i,k
      real, dimension (its:ite,kts:ktf) :: tv
      real :: e,tvbar,pp

!      real, external :: satvap
!      real :: satvap
      real, external :: rsif, rslf


!!      HT(1)=XL/CP
!!      HT(2)=2.834E6/CP
!!      BE(1)=.622*HT(1)/.286
!!      AE(1)=BE(1)/273.+ALOG(610.71)
!!      BE(2)=.622*HT(2)/.286
!!      AE(2)=BE(2)/273.+ALOG(610.71)

      do i=its,itf
      if(ierr(i).eq.0)then
      DO k=kts,ktf
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
!        IPH=1
!        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
!       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
!       QES(I,K)=.622*E/(100.*P(I,K)-E)
!!        e=satvap(t(i,k))
!!        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
!!        
        pp = 100. * p(i,k)
        if (t(i,k) .le. tcrit) then
           qes(i,k) = rsif(pp,t(i,k))
        else
           qes(i,k) = rslf(pp,t(i,k))
        endif

        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
      enddo
      endif
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
         do i=its,itf
           if(ierr(i).eq.0)then
           DO K=kts+1,ktf
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           enddo
           endif
         enddo
      else if(itest.eq.2)then
         do i=its,itf
         if(ierr(i).eq.0)then
            do k=kts,ktf
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
            enddo
         endif
         enddo
      else if(itest.eq.-1)then
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
      do i=its,itf
      if(ierr(i).eq.0)then
      DO k=kts,ktf
         if(itest.le.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
      enddo
      endif
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,ktf)

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
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
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=kts+1,ktf
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
      enddo
      endif
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
              xf_ens,name,axx,iedt,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,      &
              ktf,tau_ecmwf,aa1_bl,xf_dicycle  )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
     integer, intent (in   )              ::                           &
        iedt
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
  ! ichoice     = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:ktf,1:ens4)                              &
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
      character *(*), intent (in)         ::                           &
       name
!-srf begin
      real,    intent(IN)   , dimension (its:ite) :: aa1_bl,tau_ecmwf
      real,    intent(INOUT), dimension (its:ite) :: xf_dicycle
      !- local var
      real  :: xff_dicycle
!-srf end
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       a1,massfld,a_ave,xff0,xxx,xomg,aclim1,aclim2,aclim3,aclim4

     integer :: nall2,ixxx
     integer,  dimension (8) :: seed
     real, dimension (its:ite) :: ens_adj

     integer, parameter :: irandom = 0

       ens_adj=1.
       seed=0
       do i=its,itf
        if(ierr(i).eq.0)then
          seed(1)=int(aa0(i))
          seed(2)=int(aa1(i))
          exit
        endif
       enddo

!      nens=0
!      irandom=0
!      fens4=float(ens4)

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
!
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
                xff_ens3(2)=max(0.,(AA1(I)-AA0(I))/dtime)
!               xff_ens3(2)=max(0.,(a_ave-AA0(I))/dtime)

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
                xff_ens3(5)=xff_ens3(6)
                xff_ens3(4)=xff_ens3(6)
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
!
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
                xff_ens3(10)=AA0(i)/(60.*20.)
                xff_ens3(11)=AA0(I)/(60.*20.)
                xff_ens3(16)=AA0(I)/(60.*20.)
                if(irandom.eq.1)then
                call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AA0(I)/(60.*20.)
                else
                   xff_ens3(12)=AA0(I)/(60.*20.)
                endif
!
!-srf begin
!- more like Bechtold et al. (JAS 2014)
                !IF(DICYCLE) then
                ! xff_ens3(13)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                ! xff_ens3(14)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                ! xff_ens3(15)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                !ENDIF
!-srf implementation of 2014
!                xff_ens3(13)=(AA0(i))/tau_ecmwf(i)
!                xff_ens3(14)=(AA0(i))/tau_ecmwf(i)
!                xff_ens3(15)=(AA0(i))/tau_ecmwf(i)
!-srf implementation of 2015 - AA1 already appplied the forcing
                xff_ens3(13)=(AA1(i))/tau_ecmwf(i)
                xff_ens3(14)=(AA1(i))/tau_ecmwf(i)
                xff_ens3(15)=(AA1(i))/tau_ecmwf(i)
                
                xff_dicycle = AA1_BL(i)/tau_ecmwf(i)
!-srf end
!
!gtest
                if(ichoice.eq.0)then
                if(xff0.lt.0.)then
                     xff_ens3(1)=0.
                     xff_ens3(2)=0.
                     xff_ens3(3)=0.
                     !srf if(.not. dicycle) xff_ens3(13)=0.
                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
                     !-srf begin
                     xff_ens3(13)= 0.
                     xff_ens3(14)= 0.
                     xff_ens3(15)= 0.
                     xff_dicycle = 0.
                     !-srf end
                endif
                  if(xff0.lt.0 .and. xland(i).lt.0.1 .and. imid.eq.0)then
                     xff_ens3(:)=0.
                  endif
                endif

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(1)
                   !if(xk(nens).le.0.and.xk(nens).gt.-1.e-2) &
                   !        xk(nens)=-1.e-2
                   !lixo if(xk(nens).le.0.and.xk(nens).gt.-10.*mbdt(1)) &
                   !lixo        xk(nens)=-10.*mbdt(1)
                   if(xk(nens).le.0.and.xk(nens).gt.-.1*mbdt(1)) &
                           xk(nens)=-.1*mbdt(1)
                   
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-2) &
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
! over water, enfor!e small cap for some of the closures
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
                      !srf
                      xff_dicycle = ens_adj(i)*xff_dicycle
                      !srf end
!                     xff_ens3(7) =0.
!                     xff_ens3(8) =0.
!                     xff_ens3(9) =0.
                 endif
                endif
!
! end water treatment
!
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!
                      if(xff0.gt.0.)then
                         if(xff_ens3(1).gt.0)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))
                         if(xff_ens3(2).gt.0)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne))
                         if(xff_ens3(3).gt.0)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne))
                         !srf if(.not. dicycle) then
                         !srf  if(xff_ens3(13).gt.0)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                         !srf endif

                      else
                         xff_ens3(1)=0
                         xff_ens3(2)=0
                         xff_ens3(3)=0
                         !srf if(.not. dicycle)xff_ens3(13)=0
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4))
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5))
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6))
                         !srf if(.not. dicycle) xf_ens(i,j,nall+14)=max(0.,xff_ens3(14))

                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7)/a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8)/a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9)/a1)
                         !srf a1=max(1.e-3,pr_ens(i,j,nall+15))
                         !srf if(.not. dicycle) xf_ens(i,j,nall+15)=max(0.,xff_ens3(15)/a1)

                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
                            xf_ens(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))
                            xf_ens(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))
                            xf_ens(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne))
                         endif
!srf-begin
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                            xf_ens(i,j,nall+14)=max(0.,-xff_ens3(14)/xk(ne))
                            xf_ens(i,j,nall+15)=max(0.,-xff_ens3(15)/xk(ne))
                            !ask georg if the line below should be here -OR- out of this if statement
                            !xf_dicycle(i)      =  max(0.,-xff_dicycle /xk(ne))
                         endif
                          
                         !testing out  
                         xf_dicycle(i)     =    max(0.,-xff_dicycle /xk(ne))  
                                                   
!srf-end
                      if(ichoice.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+ichoice)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ichoice)
                      endif
                     !print*,"xf_dic2=",xf_dicycle(i),xff_dicycle,maxval(xf_ens(i,j,:));call flush(6)
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
!               else
!                  xf_ens(i,j,nall+16)=xf_ens(i,j,nall+1)
                endif
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!

!    NOT USED FOR "NORMAL" APPLICATION (maxens=1)
!
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
               !srf begin
               xf_dicycle(i) = 0.
               !srf end        
             enddo
          endif
 100   continue

   END SUBROUTINE cup_forcing_ens_3d

   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:ktf)                              &
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


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ktf
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:ktf)                              &
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
              kbcon,ktop,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:ktf)                              &
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

     DO i=its,itf
        IF (ierr(i).eq.0) then
           do k = kbcon(i)+1, ktop(i)
              DZ = Z(I,K) - Z(I,K-1)
              da = zu(i,k)*DZ*9.81/(1004.*T_cup(I,K))*DBY(I,K-1) / &
                   (1.+GAMMA_CUP(I,K))
              IF (K.eq.KTOP(I)) da = max(da,0.0)
              AA0(I)=AA0(I)+da
              if (aa0(i).lt.0.) aa0(i)=0.
           enddo
        endif
     enddo

   END SUBROUTINE cup_up_aa0

!====================================================================

   SUBROUTINE cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              name,ierr2,ierr3,pr_ens,                    &
              sig,closure_n,ktf,xf_dicycle )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in  )                   ::                           &
        zu
     real,   dimension (its:ite)                                      &
         ,intent (in  )                   ::                           &
        sig
     real,   dimension (its:ite,2)                                      &
         ,intent (in  )                   ::                           &
        xff_mid
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n
     real,    dimension (its:ite,kts:ktf,1:maxens2)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!-srf begin
      real,    intent(IN), dimension (its:ite) :: xf_dicycle
!-srf end
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
       xmb_ave
     real, dimension (5) :: weight, wm
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
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
           if(pr_ens(i,j,n).le.0.)then
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!

       xmb_w=0.
!
!-- now do feedback
!
      ddtes=100.
      if(imid.eq.0)then
      do i=its,itf
        if(ierr(i).eq.0)then
         k=0
         xmb_ave(i)=0.
         do n=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
          k=k+1
          xmb_ave(i)=xmb_ave(i)+xf_ens(i,j,n)
         enddo
         xmb_ave(i)=xmb_ave(i)/float(k)
         !srf begin
         !print*,"XMB=",xmb_ave(i),xf_dicycle(i),xmb_ave(i) - xf_dicycle(i);call  flush(6)
         if(DICYCLE) then
            xmb_ave(i)=xmb_ave(i) - xf_dicycle(i)
         endif
         !srf end
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
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
      ENDDO
      else
         do i=its,itf
         IF(ierr(i).eq.0)then
           xmb(i)=.5*sig(i)*(xff_mid(i,1)+xff_mid(i,2))
           if (ichoice.gt.0) then
              xmb(i)=sig(i)*xff_mid(i,max(ichoice,1))
           endif
         endif
         enddo
      endif

      DO i=its,itf
       IF(ierr(i).eq.0)then
         DO k=kts,ktop(i)
           dtt =0.
           dtts=0.
           dtq =0.
           dtqs=0.
           dtqc=0.
           dtpw=0.
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
           xf_ens(i,j,:)=sig(i)*xf_ens(i,j,:)*dtpw/float(maxens2)
           sub_mas(i,k)=zu(i,k)*xmb(i)
         ENDDO
       ENDIF
      ENDDO

      do i=its,itf
        if(ierr(i).eq.0)then
        do k=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
          xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
        enddo
        endif
      ENDDO

   END SUBROUTINE cup_output_ens_3d

!-------------------------------------------------------

   SUBROUTINE cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,   &
              p_cup,kbcon,ktop,cd,dby,clw_all,             &
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,     &
              ccn,rho,up_massentr,up_massdetr,psum,psumh,  &
              itest,ktf,deep)

  IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (maritime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)

  real, parameter :: bdspi = 2.1960 / bdispm

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itest,ktf
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
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
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite,kts:ktf) ::                           &
        qch,qrcb,pwh,clw_allh
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
     logical, intent(in) :: deep
!
!  local variables in this routine
!

     integer                              ::                           &
        iprop,iall,i,k,k1,k2
     real                                 ::                           &
        prop_ave,qrcb_h,dp,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real,    dimension (kts:ktf)         ::                           &
        prop_b

        prop_b(kts:ktf)=0.

        if (deep) then
           c0=.002
           iall=0
        else
           ! no precip for small clouds
           c0=0.
           iall=1
        endif

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
          if(ierr(i).eq.0)qc(i,k)=qe_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qe_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
        enddo
        enddo

      do i=its,itf
      if(ierr(i).eq.0.)then

!-- mass con
!     do k=2,kbcon(i)-1  ! orignal
      do k=2,k22(i)      ! mass cons option
!-- mass con

        DZ=Z_cup(i,K)-Z_cup(i,K-1)
        qc(i,k)=qe_cup(i,k22(i))
        qch(i,k)=qe_cup(i,k22(i))

      enddo
      endif
      enddo

        DO i=its,itf
        IF(ierr(i).eq.0) then
        DO k=k22(i)+1,ktop(i)

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

        if(qc(i,k).le.qrch)then
          qc(i,k)=qrch
!         ierr(i)=232
!         go to 100
        endif
        if(qch(i,k).le.qrch)then
          qch(i,k)=qrch
!         ierr(i)=232
!         go to 100
        endif
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
         q1 = max(q1,0.0)
         berryc0 = q1*q1*q1 / (q1*300.0 + CCNclean * BDSPI)
!        qrcb_h=qrcb(i,k)/(1.+c0*dz)
         qrcb_h=((QCH(I,K)-QRCH)*zu(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
         prop_b(k)=c0*qrcb_h*zu(i,k)/(1.e-3*berryc0+1.e-8)
         pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
         berryc=qrcb(i,k)
         qrcb(i,k)=((QCh(I,K)-QRCH)*zu(i,k)-pwh(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
         if(qrcb(i,k).lt.0.)then
           berryc0=(qrcb(i,k-1)*(.5*up_massdetr(i,k-1))-(QCh(I,K)-QRCH)*zu(i,k))/zu(i,k)*1.e-3*dz*prop_b(k)
           pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k)
           qrcb(i,k)=0.
         endif
      QCh(I,K)=QRCb(I,K)+qrch
      PWAVH(I)=PWAVH(I)+pwh(I,K)
      Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
!
! then the real berry
!
          q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
!          berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
!                ( q1 * BDSP)  ) ) !/(
!mjo avoid problems with small q1
          q1 = max(q1,0.0)
          berryc0 = q1*q1*q1 / (q1*300.0 + CCN(i) * BDSPI)
          berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
          berryc=qrc(i,k)
         qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-zu(i,k)*berryc0-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
          if(qrc(i,k).lt.0.)then
            berryc0=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/zu(i,k)
            qrc(i,k)=0.
          endif
          pw(i,k)=max(0.,berryc0*zu(i,k))
          QC(I,K)=QRC(I,K)+qrch
!
!  if not running with berry at all, do the following
!
       ELSE       !c0=.002
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        else

!---gg 10jul
         QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)      ! original form 
!        QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c1+C0)*DZ) ! new to be tested
!--
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
!        qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
!                  (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
!        PW(i,k)=c0*dz*qrc(i,k)*zu(i,k)
         if(qrc(i,k).lt.0.)then
           qrc(i,k)=0.
           pw(i,k)=0.
         endif
        endif
        QC(I,K)=QRC(I,K)+qrch
      endif !autoconv
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
         Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz

      enddo
      endif
      enddo

       prop_ave=0.
       iprop=0
       do k=kts,ktf
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0)iprop=iprop+1
       enddo
       iprop=max(iprop,1)

   END SUBROUTINE cup_up_moisture

!--------------------------------------------------------------------

   real function satvap(temp2)
     implicit none
 
     real, intent(in) :: temp2

     real :: toot, toto, eilog, tsot, tsto
     real :: ewlog, ewlog2, ewlog3, ewlog4

     real, parameter :: l10i = 1.0 / log(10.)
     real, parameter :: cil1 = 3.56654 * l10i
     real, parameter :: cil2 = log(6.1071) * l10i

     real, parameter :: cwl1 = 5.02808 * l10i
     real, parameter :: cwl4 = log(1013.246) * l10i

     if (temp2 .lt. 253.16) then   !!!! ice saturation

        toot = 273.16 / temp2
        toto = temp2 / 273.16

        eilog = -9.09718 * (toot - 1.) - cil1 * log(toot) &
              + 0.876793 * (1. - toto) + cil2
        satvap = 10. ** eilog

     else

        tsot = 373.16 / temp2
        tsto = temp2 / 373.16

        ewlog = -7.90298 * (tsot - 1.) + cwl1 * log(tsot)
        ewlog2 = ewlog - 1.3816e-07 * (10. ** (11.344 * (1. - tsto)) - 1.)
        ewlog3 = ewlog2 + .0081328 * (10. ** (-3.49149 * (tsot - 1.)) - 1.)
        ewlog4 = ewlog3 + cwl4
        satvap = 10. ** ewlog4

     endif

   end function satvap

!--------------------------------------------------------------------

   SUBROUTINE CUP_gf_sh(ktf	                 &
	     !input data                         &
	     ,DTIME                              &
	     ,PSUR				 &
             ,Z1		                 &
             ,ztexec				 &
	     ,zqexec				 &
	     ,dhdt        			 &
             ,zws				 &
	     ,zo				 &
	     ,AAEQ				 &
	     ,T 				 &
	     ,Q 				 &
             ,TN				 &
	     ,QO				 &
	     ,PO				 &
	     ,P 				 &
             ,kpbl  				 &
	     !output data			 &
	     ,ierr,ierrc   			 &
             ,xmb_out				 &
	     ,OUTT				 &
	     ,OUTQ				 &
	     ,OUTQC				 &
             ,cupclw				 &
	     ,kbcon,ktop,k22                     )  

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  !
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (inout  )                   ::                           &
        OUTT,OUTQ,OUTQC,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        xmb_out
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
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        T,PO,P,tn,dhdt
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        zws,ztexec,zqexec,Z1,PSUR,AAEQ

       
       real                                                            &
        ,intent (in   )                   ::                           &
        dtime
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
  ! ichoice     = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

     real,    dimension (its:ite,kts:ktf) ::                           &
        entr_rate_2d,he,hes,qes,z,                                     &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,                &
        xt_cup,dby,hc,qrc,zu,dbyo,qco,hco,qrco,zuo,      &
        xdby,xhc,xqrc,xzu,tempco,    &
        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq,subt,subq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKBO,                                    &
       xmbmax,XMB,PWAV,PWEV,                                &
       BU,cap_max,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,zuhe,cape
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd,fsum
      
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zufinals,zustart,zufinal,dzm1,dzp1


     integer :: icount,tun_lim,k1,k2,kbegzu,kfinalzu,kstart,jmini,levadj,nvar
     logical :: keep_going
     real xff_shal(2),blqe
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:ktf) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real :: C_up,trash2,entr_rate

     real, parameter :: zustart_sh = 0.1
     real, parameter :: zufinal_sh = 1.0

      zustart=zustart_sh
      zufinal=zufinal_sh
      entr_rate = 1.2e-3
      levadj=4
      icount=0
      day=86400.
      do i=its,itf
        xmb_out(i)=0.
        cap_max_increment(i)=25.
        ierrc(i)=" "
      enddo
!
!--- initial entrainment rate       
!--- initial detrainment rates
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
        cupclw(i,k)=0.
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
      cap_maxs=75.
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
      zkbmax=3000.
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1,psur,ierr,-1,ktf)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1,psur,ierr,-1,ktf)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,ktf)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,ktf)
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
      CALL cup_MAXIMI(HEO_CUP,1,KBMAX,K22,ierr,ktf)
       DO 36 i=its,itf
!srf     k22(i)=kpbl(i)
         k22(i)=max(2,kpbl(i)-2)
!         k22(i)=2
         if(kpbl(i).gt.3)cap_max(i)=po_cup(i,kpbl(i))
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GT.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
         endif
!           if(kpbl(i).gt.3)then
!              k22(i)=kpbl(i)
!              ierr(i)=0
!              ierrc(i)=" to zero becausof kpbl"
!            endif
!        else
!            ierrc(i)="why here? "
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0.)THEN
        if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)
             hkb(i) =sum(he_cup (i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i) 
             qkbo(i)=sum(qo_cup (i,k1:k2))/float(k2-k1+1)+ zqexec(i)
        else if(use_excess == 1) then
             hkb (i)=he_cup (i,1)+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             hkbo(i)=heo_cup(i,1)+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))            
             qkbo(i)=qo_cup (i,1)+float(use_excess)*zqexec(i)
        else
             hkb (i)=he_cup (i,1)
             hkbo(i)=heo_cup(i,1)
             qkbo(i)=qo_cup (i,1)
        endif
       endif ! ierr
      enddo

      do i=its,itf
      do k=kts,ktf
          dbyo(i,k)=hkbo(i)-heso_cup(i,k)
      enddo
      enddo
      do i=its,itf
      do k=kts,ktf
         if(dbyo(i,k).lt.0)then
            ktop(i)=k-1
            go to 441
         endif
      enddo
 441       continue
      enddo

      call cup_kbcon(ierrc,cap_max_increment,5,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max,ktf,z_cup,entr_rate,heo)
!
!--- increase detrainment in stable layers
!
      DO i=its,itf
         IF(ierr(I).eq.0.)THEN
         !print*,"kbcon/k22/kpbl=",kbcon(i),k22(i),kpbl(i);call flush(6)
            if(kbcon(i).gt.ktf-4)then
                ierr(i)=231
            endif
            do k=kts,ktf
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate*(1.3-frh)
               cd(i,k)=entr_rate_2d(i,k)
            enddo
         endif
      enddo
       call rates_up_pdf(ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
                    kstabi,k22,kbcon,kts,ktf,zuo,kpbl,deep=.false.)
      do i=its,itf
         if(ierr(i).eq.0)then
         do k=ktop(i)-1,1,-1
             if(zuo(i,k).lt.1.e-6)then
     !          k22(i)=k !<<<<<<<<<<<<<<<<<<<
               exit
             endif
         enddo

         do k=ktop(i)+1,ktf
           zuo(i,k)=0.
         enddo
      !   k22(i)=max(2,k22(i))!<<<<<<<<<<<<<<<<<<<<
         endif
      enddo
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc  (i,k)=0.
         qco (i,k)=0.
         qrco(i,k)=0.
         DBY (I,K)=0.
         hco (i,k)=0.
         DBYo(I,K)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0.)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
            qco(i,k)=qkbo(i)
            if(qco(i,k).gt.qeso_cup(i,k))then           
               qrco  (i,k)= qco(i,k)-qeso_cup(i,k)
               qco   (i,k)= qeso_cup(i,k)
               cupclw(i,k)= qrco(i,k)               
            endif
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         qco(i,k)=qkbo(i)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
         !- in-cloud saturation value
         trash=QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
                             /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)
         
         !- do not allow qco be sub-saturated at cloud base.
         qco(i,k)=max(trash,qco(i,k))
         
         if(qco(i,k)>=trash) then
           qrco(i,k)= qco(i,k)-trash
           qco (i,k)= trash
         else
           qrco(i,k)= 0.0
         endif          
         cupclw(i,k)=qrco(i,k)
       endif ! ierr
      enddo
!
!
      do 42 i=its,itf
         if(ierr(i).eq.0)then
         !srf zu (i,1)=0.
         zuo(i,1)=0.
         xzu(i,1)= zuo(i,1)
         zu (i,1)= zuo(i,1)

         !- mass entrainment and detrinament is defined on model levels
        
         do k=2,maxloc(zuo(i,:),1)
         !=> below maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
        
           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(zuo(i,k-1).gt.0.) &
              entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktf-1
         !=> above maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
           
           !-special treatment for ktop
           if(k-1==ktop(i))up_massentro(i,k-1)=0.0

           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
                     
           if(zuo(i,k-1).gt.0.) &
             cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo
         
         
         do k=2,ktf-1
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo
         do k=kbcon(i)+1,ktop(i)
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
          
         enddo
 
         do k=kbcon(i)+1,ktop(i)
          !- in-cloud saturation value
          trash =QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
                             /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)
          
          !- total water liq+vapour
          trash2  = qco(i,k-1)+qrco(i,k-1)
          qco (i,k)=   (trash2* ( zuo(i,k-1)-0.5*up_massdetr(i,k-1)) + &
                       up_massentr(i,k-1)*qo(i,k-1))   /            &
                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

!          qco (i,k)=   (qco(i,k-1)*zuo(i,k-1)-.5*up_massdetr(i,k-1)* qco(i,k-1)+ &
!                       up_massentr(i,k-1)*qo(i,k-1))   /            &
!                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!          
!          qrco(i,k)=   (qrco(i,k-1)*zuo(i,k-1)-.5*up_massdetr(i,k-1)* qrco(i,k-1)) /            &
!                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!
          
          if(qco(i,k)>=trash) then
              ! cloud liquid water
              qrco(i,k)= qco(i,k)-trash
              ! cloud water vapor 
              qco (i,k)= trash 
              
          else
              qrco(i,k)= 0.0
          endif         
          cupclw(i,k)=qrco(i,k)
           !print*,"tas=",ktop(i),k,trash,DBYo(I,K),qco (i,k),qrco(i,k)
         enddo
        
         !- 
         do k=ktop(i)+1,ktf-1
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           qco (i,k)=QESo_cup(I,K)
           qrco(i,k)=0.0
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
      endif
42    continue
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,ktf)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,ktf)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!
!--- calculate in-cloud air temperature for CAPE       
!
      do i=its,itf
         tempco(i,:)=t_cup(i,:)
         if(ierr(i)== 0)then
            !print*,"tempco",ktop(i)
           do k=1,ktop(i)+1
            tempco(i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xl*qco(i,k))
            !print*,"tempco",k,tempco(i,k),t_cup(i,k)
           enddo
          endif
      enddo
      call cup_up_cape(cape,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,tempco,qco,qrco,qo_cup,ktf)
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
!
      do i=its,itf
        trash=0.        
        if(ierr(i).eq.0)then
         do k=kts,ktop(i)+1  

           ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)

            totmas=detup-entup+zuo(i,k+1)-zuo(i,k)
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               write(0,*)k22(i),kbcon(i),ktop(i)
               write(0,123)k,detup,entup,detupk,zuo(i,k+1),zuo(i,k)      
               123     format(1X,i2,5E12.4)
            !        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo_cup(i,k+1) )-     &
                           zuo(i,k  )*(hco(i,k  )-heo_cup(i,k  ) ))*g/dp 

            !-- take out cloud liquid water for detrainment
            dellaqc(i,k)=   detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp  

            !-- condensation source term = detrained + flux divergence of 
            !-- cloud liquid water (qrco)
            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp  

            !-- water vapor budget (flux divergence of Q_up-Q_env - condensation term)
            dellaq(i,k) =-(zuo(i,k+1)*(qco(i,k+1)-qo_cup(i,k+1) ) -         &
                           zuo(i,k  )*(qco(i,k  )-qo_cup(i,k  ) ) )*g/dp &
                           - C_up

            !- check water conservation liq+condensed 
            trash=trash+ (dellaq(i,k)+dellaqc(i,k))*dp/g
         enddo   ! k
         !
         if(abs(trash)>1.e-6) then
           write(6,*)'=> not water mass cons for shallow= ',i,trash
         endif

       endif
      enddo
!

!--- using dellas, calculate changed environmental profiles
!
       mbdt=3.e-4      
      !do k=kts,ktf
      !  write(0,*)'zuo,k22 = ',k,zuo(1,k),k22(1)
      !enddo
      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            dsubh(i,k)=0.0
            dsubq(i,k)=0.0
            dsubt(i,k)=0.0
            !adding dellaqc to dellaq:
            !dellaq (i,k)= dellaq(i,k)+dellaqc(i,k)
            !dellaqc(i,k)=0.0
            !                       
            XHE(I,K)=(dsubh(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ(I,K) =(dsubq(i,k)+DELLAQ(I,K)+DELLAQC(i,k))*MBDT+QO(I,K)

            !- don't feed dellat with dellaqc if
            !- the detrainment of liquid water will be used as
            !- a source for cloud microphysics (then microphysics will
            !- evaporate this excess and provide the cooling at the
            !- detrainment region)
            !DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*(DELLAQ(I,K)+dellaqc(i,k)))
             DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*(DELLAQ(I,K)))
            
            DSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+DSUBT(I,K)-XL/CP*DELLAQC(I,K))*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      xhkb(i)=hkbo(i)+(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,ktf)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,ierr,z1,ktf)
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
!     xzu(i,:)=zuo(i,:)
      xzu(i,1:ktf)=zuo(i,1:ktf)        !ss 2/19/14
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
           kbcon,ktop,ierr,ktf)
!
! now for shallow forcing
!
       do i=its,itf

          xmb(i)=0.
          xff_shal(1:2)=0.

          if(ierr(i).eq.0)then

             xmbmax(i)=100.*(p(i,kbcon(i))-p(i,kbcon(i)+1))/(g*dtime)

!- closure from Grant (200X) based on wstar

             xff_shal(1)=.03*zws(i)

!- closure from boundary layer QE (Raymond 1995)

             blqe=0.
             trash=0.
             if (k22(i).lt.kpbl(i)+1) then
                do k=1,kbcon(i)-1
                   blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
                enddo
                trash=max((hc(i,kbcon(i))-he_cup(i,kbcon(i))),1.e1)
                xff_shal(2)=max(0.,blqe/trash)
             else
                xff_shal(2)=0.0
             endif

!- average forcings

             xmb(i) = min( 0.5 * (xff_shal(1) + xff_shal(2)), xmbmax(i) )

             if(xmb(i).eq.0.)ierr(i)=22
             if(xmb(i).eq.0.)ierrc(i)='22'
             if(xmb(i).lt.0.)then
                ierr(i)=21
                ierrc(i)='21'
             endif
          endif
        
          if (ierr(i).ne.0) then
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
                   ! max allowable tendency over tendency that would lead to 
                   ! too small mix ratios
                   trash=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
                   xmb(i)=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
                endif
             enddo
             xmb_out(i)=xmb(i)
! 
! final tendencies
!
             do k=2,ktop(i)
                outt (i,k)=(dsubt(i,k)+dellat(i,k))*xmb(i)
                outq (i,k)=(dsubq(i,k)+dellaq(i,k))*xmb(i)
                outqc(i,k)=(dellaqc(i,k)          )*xmb(i)
             enddo
          endif
       enddo

     END SUBROUTINE CUP_gf_sh

!--------------------------------------------------------

   SUBROUTINE cup_up_aa1bl(aa0,t,tn,q,qo,dtime, &
                           z,zu,kbcon,ktop,ierr,ktf)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        z,zu,t,tn,q,qo
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
     real, intent(in) :: dtime
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
        dz,dA
!
     DO i=its,itf
        AA0(I)=0.
     ENDDO

     DO i=its,itf
        if (ierr(i) .eq. 0) then
           do k = kts+1, kbcon(i)
              DZ = Z(I,K)-Z(I,K-1)
              dA = DZ*9.81*( tn(i,k)-t(i,k) + 0.608*(qo(i,k)-q(i,k)) ) / dtime
              AA0(I)=AA0(I)+dA
           enddo
        endif
     enddo

   END SUBROUTINE cup_up_aa1bl
!---------------------------------------------------------------------- 
 SUBROUTINE rates_up_pdf(ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
                         kstabi,k22,kbcon,kts,ktf,zuo,kpbl,deep)
     implicit none
     integer, intent(in) :: kts,ktf
     real, dimension (its:ite,kts:ktf),intent (inout) :: entr_rate_2d,zuo
     real, dimension (its:ite,kts:ktf),intent (in) ::p_cup, heo,heso_cup,z_cup
     real, dimension (its:ite),intent (in) :: hkbo
     integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon,kpbl
     integer, dimension (its:ite),intent (inout) :: ierr,ktop
     real, dimension (its:ite,kts:ktf) :: hcot
     logical, intent(in) :: deep

     real :: dz,dh
     real :: dby(kts:ktf)
     integer :: i,k,kdefi,kstart,kbegzu,kfinalzu

     real, parameter :: dbythresh=1.0  ! the range of this parameter is 0-1, higher => lower
                                       ! overshoting (cheque AA0 calculation)
                                       ! rainfall is too sensible this parameter
                                       ! for now, keep =1.
     DO i=its,itf
      dby(:)=0.0
      if(ierr(i).eq.0)then
        hcot(i,kbcon(i))=hkbo(i)
        dz=z_cup(i,kbcon(i))-z_cup(i,kbcon(i)-1)
        dby(kbcon(i))=(hcot(i,kbcon(i))-heso_cup(i,kbcon(i)))*dz
        
        do k=kbcon(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1) &
                      + entr_rate_2d(i,k-1)*dz*heo(i,k-1))/ &
                      (1.+0.5*entr_rate_2d(i,k-1)*dz)
           dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz

       enddo

       do k=maxloc(dby(:),1),ktf-2
          if(dby(k).lt.dbythresh*maxval(dby))then
              kfinalzu = k - 1
              ktop(i)  = kfinalzu
              go to 412
          endif
       enddo

       kfinalzu=ktf-2
       ktop(i)=kfinalzu
412    continue

       if ( deep ) then
          if(kfinalzu.le.kbcon(i)+2)then
              ierr(i)=41
              ktop(i)= 0
          else
           call get_zu_zd_pdf(ierr(i),kbcon(i),kfinalzu,zuo(i,kts:ktf),kts,ktf,kpbl(i),draft=1)
          endif
       else
          if(kfinalzu.le.kbcon(i)+1)then
              ierr(i)=41
              ktop(i)= 0
          else
              call get_zu_zd_pdf(ierr(i),k22(i),kfinalzu,zuo(i,kts:ktf),kts,ktf,kpbl(i),draft=2)
          endif
         endif
      endif
     enddo

  END SUBROUTINE rates_up_pdf

!-------------------------------------------------------------------------

  subroutine get_zu_zd_pdf(ierr,kb,kt,zu,kts,ktf,kpbli,draft)

    implicit none

    integer, intent(in) ::kb,kt,kts,ktf,kpbli
    real, intent(inout) :: zu(kts:ktf)
    real  :: zuh(kts:ktf)
    integer, intent(inout) :: ierr
    integer, intent(in) :: draft  !  1='UP', 2='SH', 3='DOWN'

    !- local var
    integer :: kk,i,k,kb_adj,kpbli_adj
    real :: beta, alpha,kratio,tunning,krmax
    !- kb cannot be at 1st level

    !-- fill zu with zeros
    zu=0.0
    zuh=0.0

    IF (draft == 1) then ! "UP"
       kb_adj=max(kb,2)
       !beta=4.  !=> must be larger than 1
                 !=> higher makes the profile sharper
                 !=> around the maximum zu 
       !- 2nd approach for beta and alpha parameters
       !- the tunning parameter must be between 0.5 (low  level max zu)
       !-                                   and 1.5 (high level max zu)
       tunning = 0.9
       beta    = 2.0/tunning
       alpha   = tunning*beta
       !
       !- this alpha constrains the location of the maximun ZU to be at 
       !- "kb_adj" vertical level
       !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
       !
       do k=kts+1,min(ktf,kt+1)
          kratio= float(k)/float(kt+1)
          zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
       enddo

    ELSEIF (draft == 2) then ! "SH"

       kb_adj=kb ! level where mass flux starts
       kpbli_adj=kpbli

       if(kpbli_adj < kb_adj) then 
          kpbli_adj = kb_adj + 1
       endif
   
       !- location of the maximum Zu
       krmax=real(kpbli_adj-kb_adj+1)/real(kt+1)  

       !- this alpha imposes the maximum zu at kpbli
       alpha=1.+5.*krmax/(1.-krmax)

       !- Beta PDF 
       do k=kts+kb_adj-1,min(ktf,kt+1)
          kratio = float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)
          zu(k)  = kratio**(alpha-1.0) * (1.0-kratio)**5
       enddo

    ELSEIF (draft == 3) then ! "DOWN"

       tunning = 0.8
       beta    = 3.0/tunning
       alpha   = tunning*beta
       zuh(:)  = 0.
       do k=kts,min(kt+1,ktf)
          kratio= float(k)/float(kt+1)
          zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
       enddo
       if(maxloc(zuh(:),1).ge.kb)then
          do k=maxloc(zuh(:),1),1,-1
             kk=kb+k-maxloc(zuh(:),1)
             if(kk.gt.1)zu(kk)=zuh(k)
          enddo
          do k=maxloc(zuh(:),1)+1,kt
             kk=kb+k-maxloc(zuh(:),1)
             if(kk.le.kt)zu(kk)=zuh(k)
          enddo
       else
          do k=2,kt ! maxloc(zuh(:),1)
             zu(k)=zuh(k-1)
          enddo
       endif
    ENDIF

    !- normalize ZU
    zu(kts:min(ktf,kt+1))= zu(kts:min(ktf,kt+1)) / maxval(zu(kts:min(ktf,kt+1)))

  end subroutine get_zu_zd_pdf

!-------------------------------------------------------------------------

  SUBROUTINE cup_up_cape(aa0,z,zu,dby,GAMMA_CUP,t_cup,  &
                         kbcon,ktop,ierr,tempco,qco,qrco, qo_cup,ktf)

    IMPLICIT NONE
    integer ,intent (in   )                   ::        &
         ktf
  
    ! aa0 = dummy array for CAPE 
    ! gamma_cup = gamma on model cloud levels
    ! t_cup = temperature (Kelvin) on model cloud levels
    ! dby = buoancy term
    ! zu= normalized updraft mass flux
    ! z = heights of model levels 
    ! ierr = error value, maybe modified in this routine
    ! tempco = in-cloud temperature (Kelvin) on model cloud levels
    ! qco    = in-cloud water vapor mixing ratio on model cloud levels
    ! qo_cup = environ water vapor mixing ratio on model cloud levels
    ! qrco   = in-cloud liquid water mixing ratio on model cloud levels

    real,    dimension (its:ite,kts:ktf) ,intent (in   )    ::        &
         z,zu,gamma_cup,t_cup,dby,tempco,qco,qrco, qo_cup
    integer, dimension (its:ite)                                      &
         ,intent (in   )                   ::                           &
         kbcon,ktop

    ! input and output
    integer, dimension (its:ite),intent (inout)   ::                  &
         ierr
    real,    dimension (its:ite)                                      &
         ,intent (out  )                   ::                           &
         aa0

    ! local variables in this routine
    integer                              ::   i,k
    real                                 ::   dz,daa0

    DO i=its,itf
       AA0(i)=0.
    ENDDO

    DO i=its,itf
       IF(ierr(i) == 0) then

          DO k=max(kbcon(i),2),ktop(i)
             DZ=Z(I,K)-Z(I,K-1)
             daa0=9.81*DZ*( (tempco(i,k)- t_cup(i,k) )/ t_cup(i,k) + &
                    0.608*  (qco   (i,k)-qo_cup(i,k) )/qo_cup(i,k) - &
                             qrco  (i,k)                             )
             AA0(I)=AA0(I)+daa0
          ENDDO

       ENDIF
    ENDDO

  END SUBROUTINE cup_up_cape

 !====================================================================

    SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,ktf,z_cup,entr_rate,heo)

   IMPLICIT NONE

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ktf  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:ktf)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max,cap_inc
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

     real, dimension (its:ite,kts:ktf),intent (in) :: z_cup,heo
     real,intent (in) :: entr_rate
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2
     real                                 ::                           &
        pbcdif,plus,hetest,dz
     real, dimension (its:ite,kts:ktf) :: hcot
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      DO 27 i=its,itf
        kbcon(i)=1
        IF(ierr(I).ne.0)GO TO 27
        KBCON(I)=K22(I)+1
        if(iloop.eq.5)KBCON(I)=K22(I)
 
       !== including entrainment for hetest
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)= ( (1.-0.5*entr_rate*dz)*hcot(i,k-1)   &
                         + entr_rate*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate*dz)
        enddo
       !==
       
        GO TO 32
 31     CONTINUE
        KBCON(I)=KBCON(I)+1
      
        IF(KBCON(I).GT.KBMAX(i)+2)THEN
           if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
           endif
           GO TO 27
        ENDIF
 32     CONTINUE
        
       !== 
       !hetest= hkb(i) 
        hetest= hcot(i,KBCON(I))     
       !==
      
       if(iloop.eq.5)then
           hetest=HKB(I)  
!          do k=1,k22(i)
!          hetest=max(hetest,he_cup(i,k))
!          enddo
        endif
     
        IF(HETEST.LT.HES_cup(I,KBCON(I)))then
           GO TO 31
        ENDIF

!       cloud base pressure and max moist static energy pressure
!       i.e., the depth (in mb) of the layer of negative buoyancy
        if(KBCON(I)-K22(I).eq.1)go to 27
      
        if(iloop.eq.5 .and. (KBCON(I)-K22(I)).le.2)go to 27
      
        PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))

!----------        
        if(iloop.eq.4)plus=cap_max(i)
        !
        ! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
        if(iloop.eq.5)plus=150.
        if(iloop.eq.5.and.cap_max(i).gt.25)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
!----------        
      
        IF(PBCDIF.GT.plus)THEN

            K22  (I)=K22(I)+1
            KBCON(I)=K22(I)+1
            
            !==     including entrainment for hetest
            hcot(i,1:k22(i)) = HKB(I)
            do k=k22(i)+1,KBMAX(i)+3
               dz=z_cup(i,k)-z_cup(i,k-1)

               hcot(i,k)= ( (1.-0.5*entr_rate*dz)*hcot(i,k-1)        &
                                  + entr_rate*dz* heo (i,k-1)        )/ &
                            (1.+0.5*entr_rate*dz)
            enddo
            !==
     
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

END MODULE module_cu_gf
