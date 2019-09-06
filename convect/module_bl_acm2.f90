module module_bl_acm2

contains

!=======================================================================

  subroutine acm2_driver( iw )

    use mem_grid, only: mza, lpw
    use mem_turb, only: vkm, vkh

    implicit none

    integer, intent(in) :: iw
    real                :: mflx
    real                :: acm2_km(mza)
    real                :: acm2_kh(mza)
    integer             :: k

    call acm2_eddyx   ( iw, mflx, acm2_kh, acm2_km)
    call acm2_scalars ( iw, mflx, acm2_kh)
    call acm2_momentum( iw,       acm2_km)

    do k = lpw(iw), mza-1
       vkm(k,iw) = acm2_km(k)
       vkh(k,iw) = acm2_kh(k)
    enddo

    vkm(lpw(iw)-1,iw) = 0.0
    vkh(lpw(iw)-1,iw) = 0.0

    vkm(mza,iw) = 0.0
    vkh(mza,iw) = 0.0

  end subroutine acm2_driver

!=======================================================================

  subroutine acm2_scalars( iw, mflx, zkh )

    use mem_grid,   only: mza, arw, volti, lpw, lsw, nsw_max, dzt, zm, dzim
    use mem_turb,   only: frac_sfc, kpblh, pblh, vkh
    use mem_basic,  only: rho, thil
    use mem_ijtabs, only: itab_w
    use misc_coms,  only: dtlm
    use tridiag,    only: tridv8, acm_matrix
    use var_tables, only: num_scalar, scalar_tab
    use mem_tend,   only: thilt
    use consts_coms,only: r8

    implicit none

    integer, intent(in)    :: iw
    real,    intent(in)    :: mflx
    real,    intent(in)    :: zkh(mza)

    real :: akodz(mza), dtom(mza), dtorho(mza)

    real(r8) :: low (mza)
    real(r8) :: dia (mza)
    real(r8) :: upp (mza)
    real(r8) :: rhs (mza,num_scalar+1)
    real(r8) :: soln(mza,num_scalar+1)

    real(r8), allocatable :: frac_sumi(:)
    real(r8), allocatable :: massflx  (:)
    real(r8), allocatable :: aa     (:,:)

    real(r8) :: fracs(nsw_max+1)
    real(r8) :: fsum

    logical :: cnvct
    integer :: k, ks, n, kk, ksmax, kpbl
    integer :: kbot, ktop, nsfc, nlev
    real    :: mbar, dens, zpbl
    real    :: dtl, dti

    kbot = lpw(iw)
    ktop = mza
    nsfc = lsw(iw)
    kpbl = kpblh(iw)
    zpbl = pblh(iw)

    nlev = ktop - kbot + 1
    dtl  = dtlm(itab_w(iw)%mrlw)
    dti  = 1.0 / dtl

    cnvct = (mflx > 1.e-8 .and. kpbl - kbot > 2)

    ! IF CONVECTIVE, COMPUTE NON-LOCAL CONVECTIVE MASS FLUX

    if (cnvct) then

       if (nsfc == 1) then
          fracs(1) = 1.0_r8
       else
          fracs(1:nsfc) = frac_sfc(1:nsfc,iw)
       endif

       if (dzt(kbot) < 0.1 * zpbl .and. fracs(1) > 0.7) then
          if (nsfc < 2) fracs(2) = 0.0
          nsfc = max(2,nsfc)
          fracs(1) = .6 * frac_sfc(1,iw)
          fracs(2) = fracs(2) + .4 * frac_sfc(1,iw)
       endif

       if (nsfc > 1) then
          fracs(nsfc) = 1._r8 - sum(fracs(1:nsfc-1))
       endif

       allocate(frac_sumi(kpbl - kbot + 1))
       allocate(massflx  (kpbl - kbot + 1))

       if (nsfc > 1) then
          fsum = fracs(1)
          frac_sumi(1) = 1._r8 / fsum
       endif

       if (nsfc > 2) then
          do ks = 2, min(nsfc-1, kpbl-kbot+1)
             fsum = fsum + fracs(ks)
             frac_sumi(ks) = 1._r8 / fsum
          enddo
       endif

       do ks = nsfc, kpbl-kbot+1
          frac_sumi(ks) = 1._r8
       enddo

       mbar = mflx / (zpbl - dzt(kbot))

       do k = kbot, min(kbot + nsfc - 1, kpbl-1)
          ks = k - kbot + 1
          dens = 0.5 * (rho(k,iw) + rho(k+1,iw))
          massflx(ks) = arw(k,iw) * mbar * (zpbl - (zm(k) - zm(kbot-1))) * dens
       enddo

       do k = kbot+nsfc, kpbl-1
          ks = k - kbot + 1
          dens = 0.5 * (rho(k,iw) + rho(k+1,iw))
          massflx(ks) = arw(kbot+nsfc-1,iw) * mbar * (zpbl - (zm(k) - zm(kbot-1))) * dens
       enddo

       massflx(kpbl-kbot+1) = 0._r8
    endif

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - SCALARS

    do k = kbot, ktop - 1
       akodz(k) = arw(k,iw) * (vkh(k,iw) + zkh(k)) * dzim(k) * 0.5
    enddo

    akodz(kbot-1) = 0.0
    akodz(ktop  ) = 0.0

    do ks = 1, nlev
       k  = ks + kbot - 1

       dtorho(k) = dtl / real(rho(k,iw))
       dtom (ks) = dtorho(k) * volti(k,iw)

       low(ks) = - dtom(ks) * akodz(k-1)
       upp(ks) = - dtom(ks) * akodz(k  )
       dia(ks) = 1._r8 - low(ks) - upp(ks)
    enddo

    ! Scalar variables with long-timestep forcing included

    do n = 1, num_scalar
       do ks = 1, nlev
          k  = ks + kbot - 1
          rhs(ks,n) = scalar_tab(n)%var_p(k,iw)  &
                    + dtorho(k) * scalar_tab(n)%var_t(k,iw)
       enddo
    enddo

    ! Load potential temperature

    do ks = 1, nlev
       k  = ks + kbot - 1
       rhs(ks,n) = thil(k,iw) + dtorho(k) * thilt(k,iw)
    enddo

    ! IF CONVECTIVE, INCLUDE NONLOCAL TERMS

    if (cnvct) then

       ksmax = min(nsfc, kpbl - kbot + 1)

       allocate(aa(mza,ksmax))
       aa = 0.0_r8

       do ks = 2, kpbl - kbot + 1
          dia(ks)   = dia(ks)   + dtom(ks)   * massflx(ks-1)
          upp(ks-1) = upp(ks-1) - dtom(ks-1) * massflx(ks-1)
       enddo

       do kk = 1, ksmax
          do ks = kk+1, kpbl - kbot + 1
             aa(ks,kk) = dtom(ks) * fracs(kk) * &
                       ( massflx(ks) * frac_sumi(ks) - massflx(ks-1) * frac_sumi(ks-1) )
          enddo
       enddo

       do ks = 1, ksmax
          dia(ks)   = dia(ks) + dtom(ks) * massflx(ks) * fracs(ks) * frac_sumi(ks)
          low(ks+1) = low(ks+1) + aa(ks+1,ks)
       enddo

      call acm_matrix(aa, dia, low, rhs, upp, soln, nlev, mza, num_scalar+1, ksmax)

    else

       call tridv8(low, dia, upp, rhs, soln, 1, nlev, mza, num_scalar+1)

    endif

    !dir$ nounroll_and_jam
    do n = 1, num_scalar

       ! Vertical loop over T levels
       do k = kbot, ktop
          ks = k - kbot + 1

          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + dti * real( (soln(ks,n) - rhs(ks,n)) * rho(k,iw) )
       enddo

    enddo

    n = num_scalar + 1

    ! Vertical loop over T levels

    do k = kbot, ktop
       ks = k - kbot + 1
       thilt(k,iw) = thilt(k,iw) + dti * real( (soln(ks,n) - rhs(ks,n)) * rho(k,iw) )
    enddo

  end subroutine acm2_scalars

!=======================================================================

  subroutine acm2_momentum( iw, zkm )

    use mem_grid,   only: mza, arw, volti, lpw, lsw, nsw_max, dzim
    use mem_turb,   only: vkm_sfc, vkm
    use mem_basic,  only: vxe, vye, vze, rho
    use mem_tend,   only: vmxet, vmyet, vmzet
    use mem_ijtabs, only: itab_w
    use misc_coms,  only: dtlm
    use tridiag,    only: tridv

    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: zkm(mza)

    real :: aflux(mza,3)
    real :: rhs  (mza,3)
    real :: soln (mza,3)

    real :: akodz(mza)

    real :: low(mza)
    real :: dia(mza)
    real :: upp(mza)
    real :: dtom(mza)
    real :: fact(nsw_max)
    real :: dtl

    integer :: k, ks
    integer :: kbot, ktop, nsfc

    kbot = lpw(iw)
    ktop = mza
    nsfc = lsw(iw)
    dtl  = dtlm(itab_w(iw)%mrlw)

    ! EDDY DIFFUSIVITY TERMS FOR SEMI-IMPLICIT SOLVER - MOMENTUM

    akodz(kbot-1) = 0.0

    do k = kbot, ktop-1
       akodz(k) = arw(k,iw) * (vkm(k,iw) + zkm(k)) * dzim(k) * 0.5
    enddo

    akodz(ktop) = 0.0

    do k = kbot, ktop
       dtom(k)  = dtl * volti(k,iw) / real(rho(k,iw))
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
       fact(ks) = 2.0 * dzim(k-1) * vkm_sfc(ks,iw) * (arw(k,iw) - arw(k-1,iw))
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

  subroutine acm2_eddyx(iw, mflx, zkh, zkm)

    ! THREE METHODS FOR COMPUTING KZ:
    !   1. Boundary scaling similar to Holtslag and Boville (1993)
    !   2. Local Kz computed as function of local Richardson # and
    !      vertical wind shear, similar to LIU & CARROLL (1996)
    !   3. Cloud-top radiational cooling scaling, Lock et al. (2000)

    use mem_cuparm,  only: iactcu, qwcon
    use mem_grid,    only: mza, zm, zt, dzm, dzt, dzimsq, lpw, lsw
    use consts_coms, only: vonk, grav, grav2, alvl, rvap, alvlocp, eps_virt
    use mem_radiate, only: cloud_frac, pbl_cld_forc
    use mem_basic,   only: vxe, vye, vze, rr_w, rr_v, thil, tair, theta, rho
    use mem_turb,    only: frac_sfck, ustar_k, wtv0_k, pblh, kpblh, &
                           ustar, wstar, moli
    use therm_lib,   only: rhovsl
    use buoyancy,    only: comp_buoy

    implicit none

    ! INPUT VARIABLES

    integer, intent(in) :: iw

    ! OUTPUT VARIABLES

    real, intent(out) :: zkh(mza) ! eddy diffusivity for scalars
    real, intent(out) :: zkm(mza) ! eddy diffusivity for momentum
    real, intent(out) :: mflx     ! nonlocal mass flux

    ! LOCAL VARIABLES

    integer :: k, kpbl, kpblm, kpblp, ks, ka, kpblpp, kpblmm
    real    :: zagl, zoh, zol, zsol, zi, hovl, fnl
    real    :: dthl, dqt, ql, cfac, qsat, gams, chi_s, dbwet, db
    real    :: ss, kprof, rlam, shear, molkm
    real    :: zk, sql, fh, fm, pr, ri
    real    :: phimi, phihi
    real    :: we, dthv, wcloud, wboyr3
    real    :: rho4

    real :: edyzh(mza) ! K-profile eddy diffusivity within PBL
    real :: edyzm(mza) ! K-profile eddy diffusivity within PBL
    real :: edyrh(mza) ! Richardson-number eddy diffusivity for heat/scalars
    real :: edyrm(mza) ! Richardson-number eddy diffusivity for momentum
    real :: edych(mza)
    real :: edycm(mza)
    real :: kzo  (mza) ! background eddy diffusivity

    real :: thetav(mza), buoy(mza)
    real :: bt(mza), bq(mza), btwet(mza), bqwet(mza)

    ! PARAMETERS

    real, parameter :: ric       =   0.25 ! critical richardson number
    REAL, PARAMETER :: RLAM_stab =  30.0
    REAL, PARAMETER :: RLAM_unst = 150.0
    REAL, PARAMETER :: GAMH      =  16.0  ! Holtslag and Boville (1993)
    REAL, PARAMETER :: EDYZ0     =   0.0  ! New Min Kz
    real, parameter :: cs        = 0.2
    real, parameter :: onethird  = 1. / 3.
    real, parameter :: alvlorvap = alvl / rvap
    real, parameter :: mvonk     = -vonk
    real, parameter :: vonk72    = 0.72 * vonk

    kzo   = edyz0
    edyzh = 0.0
    edyzm = 0.0

    edyrh = 0.0
    edyrm = 0.0

    edych = 0.0
    edycm = 0.0

    kpblpp = min(mza, kpblh(iw)+2)
    kpblp  = min(mza, kpblh(iw)+1)
    kpbl   = kpblh(iw)
    kpblm  = max(lpw(iw), kpblh(iw)-1)
    kpblmm = max(lpw(iw), kpblh(iw)-2)

    ! Virtual potential temperature

    do k = lpw(iw), mza
       ql = rr_w(k,iw) - rr_v(k,iw)
       thetav(k) = theta(k,iw) * (1.0 + eps_virt * rr_v(k,iw) - ql)
    enddo

!! Compute buoyancy terms

    call comp_buoy(iw, buoy=buoy, a=bt, b=bq, awet=btwet, bwet=bqwet)

!! Holtslag's Eddy-Diffusivity Profile Within the convective PBL

    do ks = 1, lsw(iw)
       ka = lpw(iw) + ks - 1

       if (frac_sfck(ks,iw) < 1.e-4) cycle

       zi    = max(pblh(iw) + zm(lpw(iw)-1) - zm(ka-1), zt(ka+1) - zm(ka-1))
       molkm = -grav * vonk * wtv0_k(ks,iw) / (ustar_k(ks,iw)**3 * thetav(ka))

       ! Vertical loop over W levels
       do k = ka, max(ka, kpbl-1)
          zagl     = zm(k) - zm(ka-1)
          zol      = zagl * molkm
          zsol     = min(0.1 * zi, zagl) * molkm

          if (zol <= 0.0) then
             phihi = (1.0 - gamh * zsol)**0.50
             phimi = (1.0 - gamh * zsol)**0.25
          elseif (zol <= 1.0) then
             phihi = 1.0 + 5.0 * zol
             phimi = phihi
          else
             phihi = 5.0 + zol
             phimi = phihi
          endif

          kprof    = vonk * ustar_k(ks,iw) * zagl * (1.0 - zagl/zi)**2
          edyzh(k) = edyzh(k) + frac_sfck(ks,iw) * kprof * phihi
          edyzm(k) = edyzm(k) + frac_sfck(ks,iw) * kprof * phimi

       enddo
    enddo

!! Lock's Cloud-Top Driven Eddy-Diffusivity Profile

    wcloud = 0.0

    if ( cloud_frac(kpbl,iw) + cloud_frac(kpblm,iw) > 0.05 &
         .and. pbl_cld_forc(iw) > 1.e-7 .and. kpbl > lpw(iw) ) then

       wcloud = (grav / tair(kpbl,iw) * pblh(iw) * pbl_cld_forc(iw)) ** onethird

       do k = lpw(iw), kpbl-1
          zagl     = zm(k) - zm(lpw(iw)-1)
          zoh      = zagl / pblh(iw)
          kprof    = vonk * wcloud * pblh(iw) * zoh**3 * sqrt(1.0 - zoh)

          edych(k) =        kprof
          edycm(k) = 0.75 * kprof
       enddo

    endif

!! Explicit inclusion of entrainment

    if (moli(iw) < 0.0 .or. wcloud > 1.e-7) then

       ql = max( rr_w(kpbl, iw) - rr_v(kpbl ,iw), &
                 rr_w(kpblm,iw) - rr_v(kpblm,iw))

       if (iactcu(iw) > 0) ql = ql + max(qwcon(kpbl,iw), qwcon(kpblm,iw))

       cfac = max(0.0, max(cloud_frac(kpbl,iw), cloud_frac(kpblm,iw)) &
                     - min(cloud_frac(kpbl,iw), cloud_frac(kpblp,iw)) )

       dthl = max(thil(kpbl,iw), thil(kpblp,iw), thil(kpblpp,iw)) &
            - min(thil(kpbl,iw), thil(kpblm,iw), thil(kpblmm,iw))

       dqt  = min(rr_w(kpbl,iw), rr_w(kpblp,iw), rr_w(kpblpp,iw)) &
            - max(rr_w(kpbl,iw), rr_w(kpblm,iw), rr_w(kpblmm,iw))

       dthl = max(dthl, 0.2)
       dqt  = min(dqt, -1.e-8)
       dthv = max(bt(kpbl) * dthl + bq(kpbl) * dqt, 0.0)

       wboyr3 = 0.

       if (ql > 1.e-8 .and. cfac > 0.05) then

          if (pblh(iw) > zt(kpbl) - zm(lpw(iw)-1)) then
             k = kpbl
          else
             k = kpblm
          endif

          qsat   = rhovsl(tair(k,iw)) / real(rho(k,iw))
          gams   = alvlorvap * qsat / tair(k,iw)**2
          chi_s  = ql * (1. + alvlocp * gams) / (gams * dthl - dqt)

          dbwet  = grav * (btwet(k) * dthl + bqwet(k) * dqt) / thetav(k)
          db     = grav * dthv / thetav(k)

          wboyr3 = 0.1 * chi_s**2 * max(0., -dbwet) &
                       * sqrt(db * pblh(iw)) * pblh(iw) * cfac
       endif

       zi = max(pblh(iw), zm(ka+1) - zm(ka-1))

       we = 0.24 * (wstar(iw)**3 + wcloud**3 + 15.0*ustar(iw)**3 + wboyr3) &
                 * thetav(kpblm) / (grav * zi * max(dthv,0.5))

       edych(kpbl) =        we * dzm(kpbl)
       edycm(kpbl) = 0.75 * we * dzm(kpbl)

       edych(kpblm) = max(       we * dzm(kpblm), edych(kpblm))
       edycm(kpblm) = max(0.75 * we * dzm(kpblm), edycm(kpblm))

    endif

!! Louis's Richardson-number dependent eddy diffusivity

    ! Vertical loop over W levels
    do k = lpw(iw), mza-1

       zk = vonk * (zm(k) - zm(lpw(iw) - 1))

       ss = ( (vxe(k+1,iw) - vxe(k,iw))**2 &
            + (vye(k+1,iw) - vye(k,iw))**2 &
            + (vze(k+1,iw) - vze(k,iw))**2 ) * dzimsq(k)

       ri = grav2 * buoy(k) / ( (thetav(k+1) + thetav(k)) * max(ss,1.e-5) )

       shear = max(sqrt(ss), 1.e-6)

       if (buoy(k) > 0.0) then

          rlam = min(rlam_stab, cs * dzm(k))

          sql = (zk * rlam / (rlam + zk))**2

          pr = min(1.0 + 3.7 * ri, 3.0)
          fh = 1.0 / (1.0 + ri * ( 10.0 + ri * ( 50.0 + 5000.0 * ri * ri ) ) )
          fm = fh * pr

          edyrm(k) = shear * sql * fm
          edyrh(k) = shear * sql * fh

       else

          rlam = min(rlam_unst, cs * dzm(k))

          sql = (zk * rlam / (rlam + zk))**2

          phimi = (1.0 - 16.0 * ri)**0.25
          phihi = (1.0 - 16.0 * ri)**0.50

          edyrm(k) = shear * sql * phimi * phimi
          edyrh(k) = shear * sql * phimi * phihi

       endif

    enddo

    do k = kpblp, mza-1
       if (buoy(k-1) <= 0.0 .and. buoy(k) > 0.0) then
          edyrm(k) = max(edyrm(k), 0.1 * edyrm(k-1))
          edyrh(k) = max(edyrh(k), 0.1 * edyrh(k-1))
       endif
    enddo

!! Now compute combined eddy diffusivity from surface and cloud K-profiles
!! and Richardson-number eddy diffusivity

    hovl = pblh(iw) * moli(iw)
    if ((hovl < -0.1) .and. (kpbl - lpw(iw) > 2)) then
       fnl  = min(vonk72 / ( vonk72 + (mvonk / hovl) ** onethird) - 0.1, 0.5)
       mflx = vonk * ustar(iw) * sqrt(1.0 - gamh * min(0.1*pblh(iw), dzt(k)) * moli(iw)) * fnl
    else
       fnl  = 0.
       mflx = 0.
    endif

    do k = lpw(iw), mza-1
       rho4 = 0.5 * real(rho(k+1,iw) + rho(k,iw))

       zkh(k) = max( edyrh(k), edyzh(k) + edych(k), kzo(k) )
       zkh(k) = min( zkh(k), 1000.0 )
       zkh(k) = rho4 * zkh(k)

       zkm(k) = max( edyrm(k), edyzm(k) + edycm(k), kzo(k) )
       zkm(k) = min( zkm(k), 1000.0 )
       zkm(k) = rho4 * zkm(k)
    enddo

    zkh(1:lpw(iw)-1) = 0.0
    zkh(mza)         = 0.0

    zkm(1:lpw(iw)-1) = 0.0
    zkm(mza)         = 0.0

  end subroutine acm2_eddyx

!=======================================================================

  subroutine acm2_pblhgt( iw )

    ! PURPOSE: TO CALCULATE BOUNDARY LAYER HEIGHT WITH THE METHOD
    !          PUBLISHED BY HOLTSLAG (MWR AUGUST 1990).

    use mem_grid,    only: mza, dzm, dzt, zm, zt, lpw, lsw
    use consts_coms, only: grav, grav2, eps_virt
    use mem_turb,    only: ustar, wstar, wtv0, kpblh, pblh, frac_sfc
    use mem_basic,   only: thil, vxe, vye, vze, rr_w, rr_v
    use mem_cuparm,  only: iactcu, qwcon
    use mem_radiate, only: pbl_cld_forc

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

    real    :: vxmix, vymix, vzmix, rib_sav, ql, wss, zi
    integer :: k, kbot, ktop, nsfc, kpbl

    real, parameter :: ric = 0.30  ! critical richardson number
    real, parameter :: onethird = 1. / 3.

    kbot = lpw(iw)
    ktop = mza-1
    nsfc = lsw(iw)

    if (kbot+1 >= ktop) then
       kpblh(iw) = kbot
       pblh (iw) = dzt(kbot)
       return
    endif

    do k = kbot, ktop
       ql = rr_w(k,iw) - rr_v(k,iw)
       if (iactcu(iw) > 0) ql = ql + qwcon(k,iw)

       thlv(k) = thil (k,iw) * (1.0 + eps_virt * rr_v(k,iw) - ql)
    enddo

    ! COMPUTE AN AVERAGE NEAR-SURFACE VIRTUAL POTENTIAL TEMPERATURE

    if (nsfc == 1) then
       thlvsfc = thlv(kbot)
    else
       thlvsfc = sum( thlv(kbot:kbot+nsfc-1) * frac_sfc(1:nsfc,iw) )
    endif

    ! IF SURFACE LAYER IS UNSTABLE, COMPUTE CONVECTIVE VELOCITY SCALE
    ! AND THERMAL EXCESS TEMPERATURE:

    if (wtv0(iw) >= 0.0) then
       wscale = (ustar(iw)**3 + 0.6 * wstar(iw)**3) ** onethird
       thlvsfc = thlvsfc + 8.5 * wtv0(iw) / wscale
    else
       wscale = 0.0
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
    if (pbl_cld_forc(iw) > 1.e-10) then
       wcld = 0.6 * (grav / thlvsfc * zmix * pbl_cld_forc(iw)) ** onethird
       wss  = wss + wcld**2
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
