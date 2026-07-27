module module_bl_acm2

contains

!=======================================================================

  subroutine acm2_eddyx(iw)

    ! THREE METHODS FOR COMPUTING KZ:
    !   1. Boundary scaling similar to Holtslag and Boville (1993)
    !   2. Local Kz computed as function of local Richardson # and
    !      vertical wind shear, similar to LIU & CARROLL (1996)
    !   3. Cloud-top radiational cooling scaling, Lock et al. (2000)

    use mem_cuparm,  only: iactcu, qwcon
    use mem_grid,    only: mza, zm, zt, dzm, dzt, dzimsq, lpw, lsw, dzit_bot
    use consts_coms, only: vonk, grav, grav2, alvl, rvap, eps_virt, &
                           alvlocp, alviocp, eps_virt, cpio2
    use mem_radiate, only: cloud_frac, pbl_cld_forc
    use mem_basic,   only: vxe, vye, vze, rr_v, rr_w, tair, theta, rho
    use mem_turb,    only: frac_sfc, ustar_k, wtv0_k, pblh, kpblh,     &
                           ustar, wstar, vkh, vkm, akm_sfc, frac_sfck, &
                           gam_rkm, gam_rkh, frac_norm, arw_sfc
    use mem_micro,   only: rr_c, rr_p
    use therm_lib,   only: rhovsl
    use buoyancy,    only: comp_buoy
    use oname_coms,  only: nl
    use mem_tend,    only: thilt

    implicit none

    ! INPUT VARIABLES

    integer, intent(in) :: iw

    ! LOCAL VARIABLES

    integer :: k, kpbl, kpblm, kpblp, ks, ka, kpblpp, kpblmm
    real    :: zagl, zoh, zol, zsol, zi
    real    :: dthl, dqt, dql, cfac, qsat, gams, chi_s, dbwet, db
    real    :: kprofm, kprofh, rlam, shear, molki, vels2, dissk
    real    :: zfac, wstar_k, wnonloc, pr, gam
    real    :: zk, sql, ri, dens
    real    :: phimi, phihi, phim
    real    :: we, dthv, wcloud, wboyr3
    real    :: qliqm, qliq, qliqp
    real    :: qicem, qice, qicep
    real    :: tam, ta, tap
    real    :: thlm, thl, thlp
    real    :: qtm, qt, qtp
    real    :: qlm, ql, qlp
    real    :: phim_top, phimi_top

    real :: edyzh(mza) ! K-profile eddy diffusivity within PBL
    real :: edyzm(mza) ! K-profile eddy diffusivity within PBL
    real :: edyrh(mza) ! Richardson-number eddy diffusivity for heat/scalars
    real :: edyrm(mza) ! Richardson-number eddy diffusivity for momentum
    real :: edych(mza)
    real :: edycm(mza)
    real :: kzo  (mza) ! background eddy diffusivity
    real :: gamkm(mza) ! nonlocal term
    real :: gamkh(mza) ! nonlocal term
    real :: ss   (mza) ! shear squared
    real :: dissp(mza) ! dissipation

    real :: thetav(mza), buoy(mza)
    real :: bt(mza), bq(mza), btwet(mza), bqwet(mza)

    ! PARAMETERS

    REAL, PARAMETER :: RLAM_stab = 30.0
    REAL, PARAMETER :: RLAM_unst = 30.0  ! 150.0
    REAL, PARAMETER :: GAMH      = 16.0  ! Holtslag and Boville (1993)
    REAL, PARAMETER :: EDYZ0     =  0.0  ! New Min Kz
    real, parameter :: cs        =  0.2
    real, parameter :: onethird  = 1. / 3.
    real, parameter :: alvlorvap = alvl / rvap
    real, parameter :: mvonk     = -vonk
    real, parameter :: vonk72    = 0.72 * vonk
    real, parameter :: c1_nonloc = 0.6    ! Holtslag and Boville
    real, parameter :: a_nonloc  = 7.2    ! Holtslag and Boville
    real, parameter :: pr_cap    = 1.0    ! Prandl # at top of PBL
    real, parameter :: pri_cap   = 1.0 / pr_cap
    real, parameter :: prmin     = 0.5

    kzo  (:) = edyz0
    edyzh(:) = 0.0
    edyzm(:) = 0.0

    edyrh(:) = 0.0
    edyrm(:) = 0.0

    edych(:) = 0.0
    edycm(:) = 0.0

    gamkm(:) = 0.0
    gamkh(:) = 0.0

    kpblpp = min(mza, kpblh(iw)+2)
    kpblp  = min(mza, kpblh(iw)+1)
    kpbl   = kpblh(iw)
    kpblm  = max(lpw(iw), kpblh(iw)-1)
    kpblmm = max(lpw(iw), kpblh(iw)-2)

    ! Virtual potential temperature

    do k = lpw(iw), mza
       ql = 0.
       if (allocated(rr_c)) ql = ql + rr_c(k,iw)
       if (allocated(rr_p)) ql = ql + rr_p(k,iw)
       thetav(k) = theta(k,iw) * (1.0 + eps_virt * rr_v(k,iw) - ql)
    enddo

!! Compute buoyancy terms

    if (nl%moist_buoy > 0) then
       call comp_buoy(iw, buoy, a=bt, b=bq, awet=btwet, bwet=bqwet)
    else
       call comp_buoy(iw, buoy, a=bt, b=bq)
    endif

!! Holtslag's Eddy-Diffusivity Profile Within the PBL

    phimi_top = 0.

    do ks = 1, lsw(iw)
       ka = lpw(iw) + ks - 1

       ! skip if there are no surface overlaps at this level
       if (frac_sfc(ks,iw) < 1.e-6) cycle

       zi    = max( pblh(iw) + zm(lpw(iw)-1) - zm(ka-1), dzt(ka) )
       molki = -grav * vonk * wtv0_k(ks,iw) / (ustar_k(ks,iw)**3 * thetav(ka))

       if (molki <= 0.) then
          phimi_top = phimi_top + frac_sfc(ks,iw)
       else
          zol = zi * molki
          phim_top = 1.0 + 5.0 * zol
          ! This prevents diffusivities from getting too small
          ! in stable PBL, check if necessary:
          if (zol > 1.0) phim_top = 5.0 + zol
          phimi_top = phimi_top + frac_sfc(ks,iw) / phim_top
       endif

       ! skip K-profile if PBL too shallow
       if (ka > kpbl-1) cycle

       if ( molki >= 0. ) then

          ! Vertical loop over W levels
          do k = ka, kpbl-1
             zagl = zm(k) - zm(ka-1)
             zol  = zagl * molki

             phim = 1.0 + 5.0 * zol
             ! This prevents diffusivities from getting too small
             ! in stable PBL, check if necessary:
             if (zol > 1.0) phim = 5.0 + zol

             phimi = 1.0 / phim
             zfac  = 1. - zagl / zi

             kprofm   = frac_sfc(ks,iw) * vonk * ustar_k(ks,iw) * phimi * zagl * zfac**2
             edyzm(k) = edyzm(k) + kprofm
             edyzh(k) = edyzh(k) + kprofm
          enddo

       else

          wstar_k = ( grav / thetav(ka) * zi * wtv0_k(ks,iw) )**onethird
          wnonloc = ( ustar_k(ks,iw)**3 + c1_nonloc * wstar_k**3 )**onethird

          zsol = 0.1 * zi * molki
          phihi = (1.0 - gamh * zsol)**0.5
          phimi = (1.0 - gamh * zsol)**onethird
          pr = max( phimi / phihi, prmin )

          ! Vertical loop over W levels
          do k = ka, kpbl-1
             zagl = zm(k) - zm(ka-1)
             zsol = min(0.1 * zi, zagl) * molki

             phihi = (1.0 - gamh * zsol)**0.50
             phimi = (1.0 - gamh * zsol)**onethird

             zfac = 1. - zagl / zi

             kprofm   = frac_sfc(ks,iw) * vonk * ustar_k(ks,iw) * phimi * zagl * zfac**2
             kprofh   = kprofm / pr
             edyzm(k) = edyzm(k) + kprofm
             edyzh(k) = edyzh(k) + kprofh

             gam      = a_nonloc * wstar_k / (wnonloc**2 * zi)
             gamkm(k) = gamkm(k) + gam * kprofm
             gamkh(k) = gamkh(k) + gam * kprofm
          enddo

       endif
    enddo

    ka = lpw(iw)

    do ks = 1, lsw(iw)-1
       k  = ka + ks - 1
       edyzh(k) = edyzh(k) * frac_norm(ks,iw)
       edyzm(k) = edyzm(k) * frac_norm(ks,iw)

       gamkh(k) = gamkh(k) * frac_norm(ks,iw)
       gamkm(k) = gamkm(k) * frac_norm(ks,iw)
    enddo

!! Lock's Cloud-Top Driven Eddy-Diffusivity Profile

    wcloud = 0.0

    if ( allocated(cloud_frac) .and. allocated(pbl_cld_forc) ) then

       if ( cloud_frac(kpbl,iw) + cloud_frac(kpblm,iw) > 0.05 .and. &
            pbl_cld_forc(iw) > 1.e-7 .and. kpbl > lpw(iw) ) then

          wcloud = (grav / tair(kpbl,iw) * pblh(iw) * pbl_cld_forc(iw)) ** onethird

          do k = lpw(iw), kpbl-1
             zagl     = zm(k) - zm(lpw(iw)-1)
             zoh      = zagl / pblh(iw)
             kprofm   = vonk * wcloud * pblh(iw) * zoh**3 * sqrt(1.0 - zoh)

             edycm(k) = kprofm
             edych(k) = kprofm * pri_cap
          enddo

       endif
    endif

!! Explicit inclusion of entrainment

!   if (moli(iw) < 0.0 .or. wcloud > 1.e-7) then

       if (allocated(rr_c) .and. nl%moist_buoy > 0) then
          qliqm = rr_c(kpblm,iw)
          qliq  = rr_c(kpbl ,iw)
          qliqp = rr_c(kpblp,iw)
       else
          qliqm = 0.
          qliq  = 0.
          qliqp = 0.
       endif

       if (allocated(rr_p) .and. nl%moist_buoy > 1) then
          qicem = rr_p(kpblm,iw)
          qice  = rr_p(kpbl ,iw)
          qicep = rr_p(kpblp,iw)
       else
          qicem = 0.
          qice  = 0.
          qicep = 0.
       endif

       if (iactcu(iw) > 0 .and. nl%moist_buoy > 2) then
          qliqp = qliqp + qwcon(kpblp,iw)
          qliq  = qliq  + qwcon(kpbl ,iw)
          qliqm = qliqm + qwcon(kpblm,iw)
       endif

       tam = max(tair(kpblm,iw), 253.)
       ta  = max(tair(kpbl ,iw), 253.)
       tap = max(tair(kpblp,iw), 253.)

       thlm = theta(kpblm,iw) * tam / (tam + alvlocp * qliqm + alviocp * qicem)
       thl  = theta(kpbl ,iw) * ta  / (ta  + alvlocp * qliq  + alviocp * qice )
       thlp = theta(kpblp,iw) * tap / (tap + alvlocp * qliqp + alviocp * qicep)

       qtm = rr_w(kpblm,iw)
       qt  = rr_w(kpbl ,iw)
       qtp = rr_w(kpblp,iw)

       qlp = qliqp + qicep
       ql  = qliq  + qice
       qlm = qliqm + qicem

       cfac = max(0.0, max(cloud_frac(kpbl,iw), cloud_frac(kpblm,iw)) &
                     - min(cloud_frac(kpbl,iw), cloud_frac(kpblp,iw)) )

       dthl = max(thl, thlp) - min(thl,thlm)
       dthl = max(dthl, 0.2)

       dqt = max(qt, qtp) - min(qt, qtm)
       dqt = min(dqt, -1.e-8)

       dthv = max(bt(kpbl) * dthl + bq(kpbl) * dqt, 0.0)

       wboyr3 = 0.

       if (nl%moist_buoy > 0) then

          dql = max(ql, qlm) - min(ql,qlp)

!         if ( ql > 1.e-8 .and. cfac > 0.05 .and. dthv > 1.e-9) then
          if (dql > 1.e-8 .and. cfac > 0.05 .and. dthv > 1.e-6) then

             if (pblh(iw) > zt(kpbl) - zm(lpw(iw)-1)) then
                k = kpbl
             else
                k = kpblm
             endif

             qsat   = rhovsl(tair(k,iw)) / real(rho(k,iw))
             gams   = alvlorvap * qsat / tair(k,iw)**2
!            chi_s  =  ql * (1. + alvlocp * gams) / (gams * dthl - dqt)
             chi_s  = dql * (1. + alvlocp * gams) / (gams * dthl - dqt)

             dbwet  = grav * (btwet(k) * dthl + bqwet(k) * dqt) / thetav(k)
             db     = grav * dthv / thetav(k)

             wboyr3 = 0.1 * chi_s**2 * max(0., -dbwet) &
                          * sqrt(db * pblh(iw)) * pblh(iw) * cfac
          endif
       endif

       zi = max(pblh(iw), zm(ka+1) - zm(ka-1))

       we = 0.24 * (wstar(iw)**3 + wcloud**3 + 15.0*(ustar(iw)*phimi_top)**3 + wboyr3) &
                 * thetav(kpblm) / (grav * max(zi,200.) * max(dthv,0.5))

       we = min(we, 0.5)

       edych(kpbl) =          we * dzm(kpbl)
       edycm(kpbl) = pr_cap * we * dzm(kpbl)

       edych(kpblm) = max(         we * dzm(kpblm), edych(kpblm))
       edycm(kpblm) = max(pr_cap * we * dzm(kpblm), edycm(kpblm))

!   endif

!! Louis's Richardson-number dependent eddy diffusivity

    ! Vertical loop over W levels
    do k = lpw(iw), mza-1

       zk = vonk * (zm(k) - zm(lpw(iw) - 1))

       ss(k) = ( (vxe(k+1,iw) - vxe(k,iw))**2 &
               + (vye(k+1,iw) - vye(k,iw))**2 &
               + (vze(k+1,iw) - vze(k,iw))**2 ) * dzimsq(k)

       ri = grav2 * buoy(k) / ( (thetav(k+1) + thetav(k)) * max(ss(k),1.e-5) )

       shear = max(sqrt(ss(k)), 1.e-6)

       if (ri < 0.) then

          rlam = rlam_unst
          sql = (zk * rlam / (rlam + zk))**2

          phihi = (1.0 - 16.0 * ri)**0.50
          phimi = (1.0 - 16.0 * ri)**onethird

          edyrm(k) = shear * sql * phimi * phimi
          edyrh(k) = shear * sql * phimi * phihi

       else ! test NOT limiting ri ! if (ri < 0.5) then

          rlam = rlam_stab
          sql = (zk * rlam / (rlam + zk))**2

          ! ECMWF stability functions, based on Louis et al. (1982)
          ! with modifications by Viterbo et al. (1999)
          edyrm(k) = shear * sql / (1. + 10. * ri / sqrt(1. + ri))
          edyrh(k) = shear * sql / (1. + 10. * ri * sqrt(1. + ri))

       endif

    enddo

    do k = lpw(iw), mza-1
       dens = 0.5 * real(rho(k+1,iw) + rho(k,iw))

       vkh(k,iw) = max( edyrh(k), edyzh(k) + edych(k), kzo(k) )
       vkh(k,iw) = min( vkh(k,iw), 1000.0 )
       vkh(k,iw) = dens * vkh(k,iw)

       vkm(k,iw) = max( edyrm(k), edyzm(k) + edycm(k), kzo(k) )
       vkm(k,iw) = min( vkm(k,iw), 1000.0 )
       vkm(k,iw) = dens * vkm(k,iw)

       gam_rkh(k,iw) = dens * gamkh(k)
       gam_rkm(k,iw) = dens * gamkm(k)

       dissp(k) = ss(k) * vkm(k,iw)
    enddo

    vkm    (mza,iw) = 0.0
    vkh    (mza,iw) = 0.0
    gam_rkh(mza,iw) = 0.0
    gam_rkm(mza,iw) = 0.0

    vkm    (1:ka-1,iw) = 0.0
    vkh    (1:ka-1,iw) = 0.0
    gam_rkh(1:ka-1,iw) = 0.0
    gam_rkm(1:ka-1,iw) = 0.0

    ! Compute surface dissipation

    vels2       = vxe(ka,iw)**2 + vye(ka,iw)**2 + vze(ka,iw)**2
    dissp(ka-1) = akm_sfc(1,iw) * vels2 * dzit_bot(ka)**2 / arw_sfc(1,iw)
    dissp(mza)  = 0.0

    ! Special for shaved cells (if multiple surface overlaps)

    do ks = 2, lsw(iw)
       k = ks + ka - 1
       vels2      = vxe(k,iw)**2 + vye(k,iw)**2 + vze(k,iw)**2
       dissk      = akm_sfc(ks,iw) * vels2 * dzit_bot(k)**2 / max(arw_sfc(ks,iw), 1.e-8)
       dissp(k-1) = (1. - frac_sfck(ks,iw)) * dissp(k-1) &
                  +       frac_sfck(ks,iw)  * dissk
    enddo

    ! Apply dissipative heating to THIL
    do k = lpw(iw), mza
       thilt(k,iw) = thilt(k,iw) + (dissp(k) + dissp(k-1)) * cpio2 * theta(k,iw) / tair(k,iw)
    enddo

  end subroutine acm2_eddyx

!=======================================================================

  subroutine acm2_pblhgt( iw )

    ! PURPOSE: TO CALCULATE BOUNDARY LAYER HEIGHT WITH THE METHOD
    !          PUBLISHED BY HOLTSLAG (MWR AUGUST 1990).

    use mem_grid,    only: mza, dzm, dzt, zm, zt, lpw, lsw
    use consts_coms, only: grav, grav2, eps_virt, alvlocp, alviocp
    use mem_turb,    only: ustar, wstar, wtv0, kpblh, pblh
    use mem_basic,   only: theta, tair, vxe, vye, vze, rr_v
    use mem_radiate, only: pbl_cld_forc
    use mem_micro,   only: rr_c, rr_p

    implicit none

    integer, intent(in) :: iw

    ! LOCAL VARIABLES

    real    :: thlvsfc   ! near-ground virtual temperature
    real    :: fint      ! vertical interpolation factor
    real    :: rib       ! bulk Richardson number
    real    :: zmix      ! height of mixed layer
    integer :: kmix      ! vertical index of mixed layer height
    real    :: wscale    ! convective velocity scale
    real    :: wcld      ! cloud-top cooling velocity scale
    real    :: wssq      ! wind shear term in the Richardson number
    real    :: thlv(mza) ! a 'virtual' ice-liquid pot temp

    real    :: vxmix, vymix, vzmix, rib_sav, wss, zi, thl, ta, ric
    integer :: k, kbot, ktop, nsfc, kpbl
    real    :: ql(mza), qi(mza)

    real, parameter :: onethird = 1. / 3.

    kbot = lpw(iw)
    ktop = mza-1
    nsfc = lsw(iw)

    if (kbot+1 >= ktop) then
       kpblh(iw) = kbot
       pblh (iw) = dzt(kbot)
       return
    endif

    ql(:) = 0.0
    if (allocated(rr_c)) then
       ql(kbot:ktop) = rr_c(kbot:ktop,iw)
    endif

    qi(:) = 0.0
    if (allocated(rr_p)) then
       qi(kbot:ktop) = rr_p(kbot:ktop,iw)
    endif

    do k = kbot, ktop
       ta  = max( tair(k,iw), 253. )

       thl = theta(k,iw) * ta &
           / ( ta + alvlocp * ql(k) + alviocp * qi(k) )

       thlv(k) = thl * (1.0 + eps_virt * rr_v(k,iw) - ql(k) - qi(k))
    enddo

    ! IF SURFACE LAYER IS UNSTABLE, COMPUTE CONVECTIVE VELOCITY SCALE
    ! AND THERMAL EXCESS TEMPERATURE:

    if (wtv0(iw) > 0.) then
       wscale = (ustar(iw)**3 + 0.6 * wstar(iw)**3) ** onethird
       thlvsfc = thlv(kbot) + 9.0 * wtv0(iw) / wscale
       ric     = 0.45
    else
       wscale  = ustar(iw)
       thlvsfc = thlv(kbot)
       ric     = 0.45
    endif

    ! FIND FIRST LEVEL ABOVE THE MIXED LAYER, IF IT EXISTS

    do k = kbot+1, ktop
       if (thlv(k) > thlvsfc) exit
    enddo

    ! EXIT IF NO INVERSION FOUND

    if (k > ktop) then
       kpblh(iw) = ktop
       pblh (iw) = zm(ktop)
       return
    endif

    kmix = k

    ! IF MIXED-LAYER EXISTS; INTERPOLATE HEIGHT AND WIND VELOCITIES
    ! TO THE POINT WHERE THETAV EQUALS ITS SURFACE-LAYER VALUE

    if (kmix > kbot+1) then
       fint  = (thlvsfc - thlv(kmix-1)) / max(thlv(kmix) - thlv(kmix-1), 1.e-7)
       fint  = max(min(fint, .999), 0.0)
       zmix  = fint * dzm(kmix-1) + zt(kmix-1)
       vxmix = fint * (vxe(kmix,iw)-vxe(kmix-1,iw)) + vxe(kmix-1,iw)
       vymix = fint * (vye(kmix,iw)-vye(kmix-1,iw)) + vye(kmix-1,iw)
       vzmix = fint * (vze(kmix,iw)-vze(kmix-1,iw)) + vze(kmix-1,iw)
    else
       zmix  = zt(kbot)
       vxmix = vxe(kbot,iw)
       vymix = vye(kbot,iw)
       vzmix = vze(kbot,iw)
    endif

    ! STARTING FROM THE TOP OF THE MIXED LAYER (OR SURFACE IF NO MIXED LAYER
    ! EXISTS), SEARCH UPWARD UNTIL THE BULK RICHARDSON NUMBER EQUALS ITS
    ! CRITICAL VALUE. THIS WILL BE THE PBL HEIGHT

    wss = wscale**2

    if (allocated(pbl_cld_forc)) then
       if (pbl_cld_forc(iw) > 1.e-10) then
          wcld = 0.6 * (grav / thlv(kmix) * zmix * pbl_cld_forc(iw)) ** onethird
          wss  = wss + wcld**2
       endif
    endif

    rib  = 0.0
    do k = kmix, ktop
       rib_sav = rib

       wssq = wss + (vxe(k,iw) - vxmix)**2  &
                  + (vye(k,iw) - vymix)**2  &
                  + (vze(k,iw) - vzmix)**2

       wssq = max(wssq, 0.1)

       rib = grav2 * (thlv(k) - thlvsfc) * (zt(k) - zmix)  &
           / ( (thlv(k) + thlvsfc) * wssq )

       if (rib > ric) exit
    enddo

    kpbl = min(k,ktop)

    if (kpbl == ktop) then

      zi = max( zm(kpbl) - zm(kbot-1), dzt(kbot) )

    else

       ! ZI IS BETWEEN ZT(KPBL-1) and ZT(KPBL).
       ! INTERPOLATE BETWEEN LEVELS TO DETERMINE THE PBL HEIGHT.

       fint = (ric - rib_sav) / max(rib - rib_sav, 1.e-6)
       fint = max(min(fint, 1.0), 0.0)

       zi   = fint * dzm(kpbl-1) + zt(kpbl-1) - zm(kbot-1)

       if (zi + zm(kbot-1) <= zm(kpbl-1)) then
          kpbl = kpbl - 1
       endif

    endif

    kpblh(iw) = kpbl
    pblh (iw) = zi

  end subroutine acm2_pblhgt

END MODULE module_bl_acm2
