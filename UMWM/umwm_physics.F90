module umwm_physics

contains

!===============================================================================

subroutine umwm_diag(i)

  use consts_coms, only: pi2
  use umwm_module, only: evs, cth, sth, kdk, freq, k, kdk, dth, th, ucurr, vcurr, &
                         mwp, mwd, mss, mwl, dwd, swh, dwp, dcp0, dcg0, dcp, dcg, &
                         opeak, ppeak, dwl, om, pm, cp0, cg0, umwm
  implicit none

  integer, intent(in) :: i

  integer             :: o, opk, ppk, pmax(om)
  real                :: m0, m2, m3, xcomp, ycomp
  real                :: evsum(om), exsum(om), eysum(om), evmax(om)

  if (.not. umwm%iactive(i)) then

     mwp  (i) = 0.0
     mwd  (i) = 0.0
     swh  (i) = 0.0
     mss  (i) = 0.0
     mwl  (i) = 0.0
     dwd  (i) = 0.0
     dwl  (i) = 0.0
     dwp  (i) = 0.0
     dcp0 (i) = 0.0
     dcg0 (i) = 0.0
     dcp  (i) = 0.0
     dcg  (i) = 0.0
     opeak(i) = 0.0
     ppeak(i) = 0.0

     return
  endif

  do o = 1, om
     evsum(o) = sum( evs(1:pm,o,i) )
     exsum(o) = sum( evs(1:pm,o,i) * cth(1:pm) )
     eysum(o) = sum( evs(1:pm,o,i) * sth(1:pm) )
     pmax (o) = maxloc( evs(1:pm,o,i), dim=1 )
     evmax(o) = evs(pmax(o),o,i)
  enddo

  evmax(om) = max(evmax(om),1.e-30)

  evsum(1:om) = evsum(1:om) * kdk(1:om,i)
  evmax(1:om) = evmax(1:om) * kdk(1:om,i)

  m0    = sum( evsum(1:om)                   )
  m2    = sum( evsum(1:om) * freq(1:om)**2 )
  m3    = sum( evsum(1:om) * k  (1:om,i)**2  )
  xcomp = sum( exsum(1:om) * kdk(1:om,i)     )
  ycomp = sum( eysum(1:om) * kdk(1:om,i)     )

  opk   = maxloc( evmax(1:om), dim=1 )
  ppk   = pmax(opk)

  mwp(i) = sqrt( m0 / max(m2,1.e-30) )       ! mean wave period
  mwd(i) = atan2(ycomp,xcomp)                ! mean wave direction
  swh(i) = 4. * sqrt(m0 * dth)               ! significant wave height

  mss(i) = m3 * dth                          ! mean-squared slope
  mwl(i) = pi2 * sqrt( m0 / max(m3,1.e-30) ) ! mean wavelength

  dwd(i) = th(ppk)                         ! dominant wave direction
  dwl(i) = pi2 / k(opk,i)                  ! dominant wave length
  dwp(i) = 1. / freq(opk)                  ! dominant wave period

  dcp0(i) = cp0(opk,i)                     ! dominant phase speed, intrinsic
  dcg0(i) = cg0(opk,i)                     ! dominant group speed, intrinsic

  dcp(i) = dcp0(i) + ucurr(i) * cth(ppk) + vcurr(i) * sth(ppk) ! dominant phase speed
  dcg(i) = dcg0(i) + ucurr(i) * cth(ppk) + vcurr(i) * sth(ppk) ! dominant group speed

  opeak(i) = opk  ! store frequency of peak wave
  ppeak(i) = ppk  ! store direction of peak wave

end subroutine umwm_diag

!======================================================================

subroutine umwm_diag_orig(i)

  use consts_coms, only: pi2
  use umwm_module, only: evs, cth, sth, kdk, freq, k, kdk, dth, th, ucurr, vcurr, &
                         mwp, mwd, mss, mwl, dwd, swh, dwp, dcp0, dcg0, dcp, dcg, &
                         opeak, ppeak, dwl, om, pm, cp0, cg0, k3dk

  implicit none

  integer, intent(in) :: i

  integer :: o, p, opk, ppk

  real :: mag, xcomp, ycomp
  real :: m0, m2
  real :: evskdk(pm,om)
  real :: evskdkmax
! real :: evskdkovcp

  m0       = 0.
  m2       = 0.
! momx (i) = 0.
! momy (i) = 0.
! cgmxx(i) = 0.
! cgmxy(i) = 0.
! cgmyy(i) = 0.

  evskdkmax = -1.e6
  do o = om,1,-1
     do p = 1,pm
        evskdk(p,o) = evs(p,o,i) * kdk(o,i)
        m0 = m0 + evskdk(p,o)
        m2 = m2 + freq(o)**2 * evskdk(p,o)

        if (evskdkmax < evskdk(p,o)) then
           evskdkmax = evskdk(p,o)
           opk = o
           ppk = p
        endif
     enddo
  enddo
  mwp(i) = sqrt(m0 / (m2 + tiny(m2))) ! mean wave period

  ! total wave momentum:
! do o = 1,oc(i)
!    do p = 1,pm
!       evskdkovcp = evskdk(p,o) * invcp0(o,i)
!
!       momx(i) = momx(i) + evskdkovcp * cth(p)
!       momy(i) = momy(i) + evskdkovcp * sth(p)
!
!       cgmxx(i) = cgmxx(i) + cg0(o,i) * evskdkovcp * cth(p)**2
!       cgmxy(i) = cgmxy(i) + cg0(o,i) * evskdkovcp * cth(p) * sth(p)
!       cgmyy(i) = cgmyy(i) + cg0(o,i) * evskdkovcp * sth(p)**2
!    enddo
! enddo
!
! momx (i) = momx (i) * rhosw * dthg
! momy (i) = momy (i) * rhosw * dthg
! cgmxx(i) = cgmxx(i) * rhosw * dthg
! cgmxy(i) = cgmxy(i) * rhosw * dthg
! cgmyy(i) = cgmyy(i) * rhosw * dthg

  ! significant wave height:
  swh(i) = 0.
  do o = 1,om
     do p = 1,pm
        swh(i) = swh(i) + evskdk(p,o)
     enddo
  enddo
  swh(i) = 4. * sqrt(swh(i) * dth)

  ! mean wave direction:
  xcomp = 0.
  ycomp = 0.
  do o = 1, om
     do p = 1, pm
        xcomp = xcomp + mag * cth(p)
        ycomp = ycomp + mag * sth(p)
     enddo
  enddo
  mwd(i) = atan2(ycomp,xcomp)

  ! wavenumber spectrum moments:
  m0 = 0.
  m2 = 0.
  do o = 1,om
     do p = 1,pm
        m0 = m0 + evskdk(p,o)
        m2 = m2 + evs(p,o,i) * k3dk(o,i)
     enddo
  enddo

  mss(i) = m2 * dth ! mean-squared slope

  mwl(i) = pi2 * sqrt(m0 / (m2 + tiny(m2))) ! mean wavelength

  dwd(i) = th(ppk)                       ! dominant wave direction
  dwl(i) = pi2 / k(opk,i)                ! dominant wave length
  dwp(i) = 1. / freq(opk)              ! dominant wave period

  dcp0(i) = cp0(opk,i)                   ! dominant phase speed, intrinsic
  dcg0(i) = cg0(opk,i)                   ! dominant group speed, intrinsic

  dcp(i) = dcp0(i) + ucurr(i) * cth(ppk) + vcurr(i) * sth(ppk) ! dominant phase speed
  dcg(i) = dcg0(i) + ucurr(i) * cth(ppk) + vcurr(i) * sth(ppk) ! dominant group speed

  opeak(i) = opk  ! store frequency of peak wave
  ppeak(i) = ppk  ! store direction of peak wave

end subroutine umwm_diag_orig

!======================================================================

end module umwm_physics
