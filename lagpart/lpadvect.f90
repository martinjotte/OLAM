subroutine lpadvect()

  use mem_grid,    only: lpw, wnx, wny, wnz, xew, yew, zew, zt, dzim
  use mem_lp,      only: atp, nlp, lpturb, vxeh, vyeh, vzeh, dtlp, nent, anorm
  use consts_coms, only: erad, eradi
  use mem_ijtabs,  only: itab_w
  use map_proj,    only: de_ps
  use misc_coms,   only: io6

  implicit none (external, type)

  integer :: iw, iwn, iw1, iw2, iw3, iwstrt
  integer :: j, j1, j2
  integer :: k
  integer :: ilp, npoly

  real :: raxis, raxisi, sinwlat, coswlat, sinwlon, coswlon
  real :: dxep, dyep, dzep, factor
  real :: dxe, dye, dze
  real :: dot00(7), dot01(7), dot11(7)
  real :: dot02, dot12
  real :: denomi(7), dn(7), xwn(7), ywn(7)
  real :: wh1, wh2, wh3, wz11, wz21,  wz12, wz22, wz13, wz23
  real :: sigu, sigv, sigw, dsw2
  real :: acu, acv, acw, acu2, acv2, acw2
  real :: tlulp, tlvlp, tlwlp
  real :: uppp, vppp, wppp
  real :: wdrift
  real :: qx, qy
  real :: u, v
  real :: entfac
  real :: vxep, vyep, vzep

  real, parameter :: fuzz = 0.001

  integer, save :: ix1 = 1000, ix2 = 2000, ix3 = 3000
  integer, save :: ia = 129, ib = 1, it = 1048576

  integer :: k11, k12, k13, k21, k22, k23
  real    :: dzp

  external :: lpe_kiw

  entfac = real(nent-1) / real(it)

  ! Loop over all particles

  write(io6,'(A,I0)') " LP particles currently advected: ", count(atp(1:nlp)%ncluster == 1)

  do ilp = 1,nlp

     ! Skip this particle if it is not active

     if (atp(ilp)%ncluster == 0) cycle

     ! Current particle grid cell indices

     iw = atp(ilp)%iw
     k  = atp(ilp)%k

     ! THIS SECTION COULD BE PRECOMPUTED FOR ALL IW COLUMNS, STORING:
     !    coswlat(iw), sinwlat(iw), coswlon(iw), sinwlon(iw),
     !    xwn(j,iw), ywn(j,iw), dot00(j,iw), dot01(j,iw), dot11(j,iw),
     !    denomi(j,iw), dn(j,iw).  HOWEVER, THIS IS PROBABLY NOT WORTH
     !    DOING SINCE PLOTTING IS INFREQUENT COMPARED TO TIME STEPPING.

     ! Prepare iw cell geometry for interpolation to particle location

     raxis = sqrt( xew(iw)**2 + yew(iw)**2 )

     sinwlat = zew(iw) * eradi
     coswlat = raxis   * eradi

     ! For points less than 100 m from Earth's polar axis, make arbitrary
     ! assumption that longitude = 0 deg.  This is just to settle on a PS
     ! planar coordinate system in which to do the algebra.

     if (raxis >= 1.e2) then
        raxisi = 1.0 / raxis
        sinwlon = yew(iw) * raxisi
        coswlon = xew(iw) * raxisi
     else
        sinwlon = 0.
        coswlon = 1.
     endif

     ! Loop over neighbor W points

     npoly = itab_w(iw)%npoly
     do j = 1,npoly
        iwn = itab_w(iw)%iw(j)

        ! Transform neighbor W points to PS coordinates tangent at IW point

        dxe = xew(iwn) - xew(iw)
        dye = yew(iwn) - yew(iw)
        dze = zew(iwn) - zew(iw)

        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xwn(j),ywn(j))

     enddo

     ! Loop over each pair of consecutive neighbor W points, and set up
     ! triangle-check coefficients that depend only on W points

     do j1 = 1,npoly
        j2 = j1 + 1
        if (j1 == npoly) j2 = 1

        dot00(j1) = xwn(j1) * xwn(j1) + ywn(j1) * ywn(j1)
        dot01(j1) = xwn(j1) * xwn(j2) + ywn(j1) * ywn(j2)
        dot11(j1) = xwn(j2) * xwn(j2) + ywn(j2) * ywn(j2)

        denomi(j1) = 1. / (dot00(j1) * dot11(j1) - dot01(j1) * dot01(j1))
        dn    (j1) = 1. / (xwn(j2) * ywn(j1) - xwn(j1) * ywn(j2))
     enddo

     ! END OF SECTION THAT COULD BE PRECOMPUTED FOR ALL IW COLUMNS

     ! Transform particle location to PS coordinates tangent at IW point.

     dxe = atp(ilp)%xep0 - xew(iw)
     dye = atp(ilp)%yep0 - yew(iw)
     dze = atp(ilp)%zep0 - zew(iw)

     call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,qx,qy)

     ! Loop over each pair of consecutive neighbor W points

     do j1 = 1,npoly
        j2 = j1 + 1
        if (j1 == npoly) j2 = 1

        ! Set up triangle-check coefficients that depend on lat/lon points

        dot02 = xwn(j1) * qx + ywn(j1) * qy
        dot12 = xwn(j2) * qx + ywn(j2) * qy

        u = (dot11(j1) * dot02 - dot01(j1) * dot12) * denomi(j1)
        v = (dot00(j1) * dot12 - dot01(j1) * dot02) * denomi(j1)

        if (u > -fuzz .and. v > -fuzz .and. u + v < 1.0 + fuzz) then

           ! particle is inside current triangle; store indices of
           ! both IWN neighbors

           iw1 = iw
           iw2 = itab_w(iw)%iw(j1)
           iw3 = itab_w(iw)%iw(j2)

           ! Compute and store 3 horizontal interpolation weights

           wh2 = dn(j1) * (-ywn(j2) * qx + xwn(j2) * qy)
           wh3 = dn(j1) * ( ywn(j1) * qx - xwn(j1) * qy)
           wh1 = 1. - wh2 - wh3

           exit
        endif
     enddo  ! j1 loop

!! Am considering whether it is correct to only use vxeh, vyeh, vzeh.
!! Particles should also (presumably) not be allowed to cross a closed
!! V face of a grid cell since the air doesn't cross it.  Therefore,
!! should particles instead use vmasc and wmasc?  For now, we try the
!! following approach using vxeh, vyeh, vzeh.

     ! Vertical interpolation weights, which may differ between iw1, iw2, and iw3
     ! when particles are adjacent to closed cells

     if (atp(ilp)%zp <= zt(lpw(iw1))) then
        k21 = lpw(iw1)
        wz21 = 1.0
     elseif (atp(ilp)%zp <= zt(k)) then
        k21 = k
        wz21 = (atp(ilp)%zp - zt(k-1)) * dzim(k-1)
     else
        k21 = k+1
        wz21 = (atp(ilp)%zp - zt(k)) * dzim(k)
     endif
     k11 = max(lpw(iw1),k21-1)
     wz11 = 1.0 - wz21

     if (atp(ilp)%zp <= zt(lpw(iw2))) then
        k22 = lpw(iw2)
        wz22 = 1.0
     elseif (atp(ilp)%zp <= zt(k)) then
        k22 = k
        wz22 = (atp(ilp)%zp - zt(k-1)) * dzim(k-1)
     else
        k22 = k+1
        wz22 = (atp(ilp)%zp - zt(k)) * dzim(k)
     endif
     k12 = max(lpw(iw2),k22-1)
     wz12 = 1.0 - wz22

     if (atp(ilp)%zp <= zt(lpw(iw3))) then
        k23 = lpw(iw3)
        wz23 = 1.0
     elseif (atp(ilp)%zp <= zt(k)) then
        k23 = k
        wz23 = (atp(ilp)%zp - zt(k-1)) * dzim(k-1)
     else
        k23 = k+1
        wz23 = (atp(ilp)%zp - zt(k)) * dzim(k)
     endif
     k13 = max(lpw(iw3),k23-1)
     wz13 = 1.0 - wz23

     ! Interpolate velocity components, sigmas, and Lagrangian time scales
     !    from OLAM grid to particle location

     vxep = wh1 * (wz11 * vxeh(k11,iw1) + wz21 * vxeh(k21,iw1)) &
          + wh2 * (wz12 * vxeh(k12,iw2) + wz22 * vxeh(k22,iw2)) &
          + wh3 * (wz13 * vxeh(k13,iw3) + wz23 * vxeh(k23,iw3))

     vyep = wh1 * (wz11 * vyeh(k11,iw1) + wz21 * vyeh(k21,iw1)) &
          + wh2 * (wz12 * vyeh(k12,iw2) + wz22 * vyeh(k22,iw2)) &
          + wh3 * (wz13 * vyeh(k13,iw3) + wz23 * vyeh(k23,iw3))

     vzep = wh1 * (wz11 * vzeh(k11,iw1) + wz21 * vzeh(k21,iw1)) &
          + wh2 * (wz12 * vzeh(k12,iw2) + wz22 * vzeh(k22,iw2)) &
          + wh3 * (wz13 * vzeh(k13,iw3) + wz23 * vzeh(k23,iw3))

     if (lpturb == 1) then

        ! For now, set turb quantities to zero

        dsw2  = 0. ! = metc(mgp)%ds2
        sigu  = 0. ! = metc(mgp)%sigu
        sigv  = 0. ! = metc(mgp)%sigv
        sigw  = 0. ! = metc(mgp)%sigw
        tlulp = 0. ! = metc(mgp)%tlu
        tlvlp = 0. ! = metc(mgp)%tlv
        tlwlp = 0. ! = metc(mgp)%tlw

        ! Compute turbulent and drift correction velocity based on sigs and tls,
        ! and add to ulp, vlp, and wlp

        ! THE SIGMAS COMPUTED IN THE TURBULENCE ROUTINES ARE ACTUALLY VARIANCES.
        ! CONVERT THEM TO STANDARD DEVIATIONS HERE FOR WHAT FOLLOWS.

        sigu = sqrt(sigu)
        sigv = sqrt(sigv)
        sigw = sqrt(sigw)

        ! ALSO, DSW2 HAS BEEN CORRECTED IN HTURB.F.

        acu = exp(-dtlp / tlulp)
        acv = exp(-dtlp / tlvlp)
        acw = exp(-dtlp / tlwlp)

        acu2 = acu * acu
        acv2 = acv * acv
        acw2 = acw * acw

        ! compute random turbulent velocity components as products of
        ! normally-distributed random numbers gx and the standard deviations
        ! computed in equation (6).  (THEY ARE NOT STANDARD DEVIATIONS).

        ix1 = mod(ia*ix1+ib,it)
        ix2 = mod(ia*ix2+ib,it)
        ix3 = mod(ia*ix3+ib,it)

        uppp = sigu * sqrt(1. - acu2) * anorm(1 + int(entfac * real(ix1)))
        vppp = sigv * sqrt(1. - acv2) * anorm(1 + int(entfac * real(ix2)))
        wppp = sigw * sqrt(1. - acw2) * anorm(1 + int(entfac * real(ix3)))

        ! LADM FORMULA FOR GENERAL INHOMOGENEOUS GAUSSIAN TURBULENCE

!e        wdrift = dtlp * 0.5 * dsw2 * (1.0 + ( atp(l)%wpp / sigw )**2.)

        ! Compute Lagrangian turbulent velocity fluctuation from equation (5).

!e        atp(ilp)%upp = atp(ilp)%upp * acu + uppp
!e        atp(ilp)%vpp = atp(ilp)%vpp * acv + vppp
!e        atp(ilp)%wpp = atp(ilp)%wpp * acw + wppp

     else  ! zero out turbulent velocities and drift velocity

!e        atp(ilp)%upp = 0.
!e        atp(ilp)%vpp = 0.
!e        atp(ilp)%wpp = 0.
        wdrift     = 0.

     endif

     ! ADD UPP, VPP, WPP TO VXEPP, VYEPP, VZEPP
     ! ADD WDRIFT TO VXEP, VYEP, VZEP

     ! RAMS VERSION:   z_new = z_old + (wlp + atp(l)%wpp + wdrift) * dtlp

     ! Compute particle displacement from mean and turbulent flow contributions

     dxep = (vxep + atp(ilp)%vxepp) * dtlp
     dyep = (vyep + atp(ilp)%vyepp) * dtlp
     dzep = (vzep + atp(ilp)%vzepp) * dtlp

     ! Project particle displacement onto local vertical to get change in zp

     dzp = dxep * wnx(iw) + dyep * wny(iw) + dzep * wnz(iw)
     atp(ilp)%zp = atp(ilp)%zp + dzp - atp(ilp)%fallvel * dtlp

     ! Add particle displacement to particle earth coordinates and then
     ! project them vertically onto earth sphere

     factor = erad / (erad + atp(ilp)%zp)
     atp(ilp)%xep0 = (atp(ilp)%xep0 + dxep) * factor
     atp(ilp)%yep0 = (atp(ilp)%yep0 + dyep) * factor
     atp(ilp)%zep0 = (atp(ilp)%zep0 + dzep) * factor

     ! Get particle k,iw from its xep0, yep0, and zep0, and zp

     iwstrt = atp(ilp)%iw
     call lpe_kiw(ilp, iwstrt)

  enddo

end subroutine lpadvect

!==============================================================================

subroutine lpe_kiw(ilp, iwstrt)

  use consts_coms, only:
  use mem_grid,    only: mza, lpw, zm, xew, yew, zew
  use mem_lp,      only: atp
  use mem_ijtabs,  only: itab_w

  implicit none (external, type)

  integer, intent(in) :: ilp     ! Lagrangian particle index
  integer, intent(in) :: iwstrt  ! Most recent LP iw index

  integer :: klo, khi, kmid, iw0, iw, iwn, npoly, j
  real    :: dist, dist0

  ! Particle distance to starting IW center

  iw0   = iwstrt
  dist0 = sqrt((atp(ilp)%xep0 - xew(iw0))**2 + (atp(ilp)%yep0 - yew(iw0))**2 + (atp(ilp)%zep0 - zew(iw0))**2)

  ! Iterate as precaution for high Courant number (or for using source iw at initial LP release)

  iw = 1
  do while (iw /= iw0)
     iw = iw0

     npoly = itab_w(iw)%npoly
     do j = 1,npoly
        iwn = itab_w(iw)%iw(j)

        dist = sqrt((atp(ilp)%xep0 - xew(iwn))**2 + (atp(ilp)%yep0 - yew(iwn))**2 + (atp(ilp)%zep0 - zew(iwn))**2)

        if (dist0 > dist) then
           iw0   = iwn
           dist0 = dist
        endif
     enddo
  enddo

  atp(ilp)%iw = iw0  ! New particle cell IW

  ! Now that iw0 is found, check if particle is below ground in iw0 cell
  ! Apply bottom boundary condition here, but for now, if particle
  ! went below ground or above model top, set ncluster = 0.

  if (atp(ilp)%zp < zm(lpw(iw0)-1)) then
     atp(ilp)%ncluster = 0
  elseif (atp(ilp)%zp > zm(mza)) then
     atp(ilp)%ncluster = 0
  else
     atp(ilp)%ncluster = 1
     klo = lpw(iw0)
     khi = mza
     do while (khi - klo > 1)
        kmid = (klo + khi) / 2
        if (atp(ilp)%zp < zm(kmid)) then
           khi = kmid
        else
           klo = kmid
        endif
     enddo
     atp(ilp)%k = khi
  endif

end subroutine lpe_kiw
