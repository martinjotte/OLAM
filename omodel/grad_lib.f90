module grad_lib

contains


!=======================================================================

! This routine computes the vertical derivative of a scalar or velocity
! vector component with respect to the local upward Z direction at each
! cell center (T level).

subroutine grad_z(iw, scp, gzps)

  use mem_grid,   only: mza, lpw, dzim

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: scp (mza)
  real,    intent(out) :: gzps(mza)

  integer :: kb, k
  real    :: gwz(mza)

  kb = lpw(iw)

  ! Vertical loop over W levels
  do k = kb, mza-1
     gwz(k) = dzim(k) * (scp(k+1) - scp(k))
  enddo

  ! Constant gradient top and bottom:
  ! gwz(kb-1) = gwz(kb)
  ! gwz(mza)  = gwz(mza-1)

  ! Zero-gradient top and bottom:
  gwz(kb-1) = 0.0
  gwz(mza)  = 0.0

  ! Vertical loop over T levels
  do k = kb, mza
     gzps(k) = 0.5 * (gwz(k-1) + gwz(k))
  enddo

end subroutine grad_z

!=======================================================================

! This routine computes the first and second derivatives of a scalar or
! velocity vector component with respect to the local upward Z direction
! at each cell center (T level). The second derivative is bounded so that
! the integral of the polynomial:
!
! scp + z * gzps + z^2 * gzzps
!
! over a cell equals its mean value scp0.

subroutine grad_z_quad(iw, scp, gzps, gzzps)

  use mem_grid,   only: mza, lpw
  use mem_adv,    only: a_v

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: scp  (mza)
  real,    intent(out) :: gzps (mza)
  real,    intent(out) :: gzzps(mza)

  integer :: k, ka
  real    :: ds(mza)

  ka = lpw(iw)

  ! Loop over W levels
  do k = ka, mza-1
     ds(k) = scp(k+1) - scp(k)
  enddo

  ! Zero gradient bottom
  ds(ka-1) = 0.0

  ! Constant gradient bottom
  ! ds(ka-1) = ds(ka)

  ! Zero gradient top
  ! ds(mza) = 0.0

  ! Constant gradient top
  ds(mza) = ds(mza-1)

  ! Loop over T levels
  do k = ka, mza
     gzzps(k) = ds(k) * a_v(k,2,1) + ds(k-1) * a_v(k,2,2)
     gzps (k) = ds(k) * a_v(k,1,1) + ds(k-1) * a_v(k,1,2)
  enddo

end subroutine grad_z_quad

!=========================================================================

! Computes horizontal gradients in a local rotated polar
! stereographic projection tangent to the local W

subroutine grad_t2d(iw, scp, gxps, gyps)

  use mem_ijtabs, only: itab_w
  use mem_grid,   only: mza, mwa, lpv, gxps_coef, gyps_coef

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: scp (mza,mwa)
  real,    intent(out) :: gxps(mza)
  real,    intent(out) :: gyps(mza)

  integer :: n, iwn, ivn, k
  real    :: dscp

  gxps = 0.
  gyps = 0.

  ! Loop over W neighbors of this W cell

  !dir$ loop count max=7
  do n = 1, itab_w(iw)%npoly
     iwn = itab_w(iw)%iw(n)
     ivn = itab_w(iw)%iv(n)

     do k = lpv(ivn), mza
        dscp    = scp(k,iwn) - scp(k,iw)
        gxps(k) = gxps(k) + gxps_coef(n,iw) * dscp
        gyps(k) = gyps(k) + gyps_coef(n,iw) * dscp
     enddo

  enddo

end subroutine grad_t2d

!=========================================================================

! Computes terms of a horizontal quadratic interpolation in a local
! rotated polar stereographic projection tangent to the local W by
! least-squares fitting. The polynomial is bounded so that the integral
! over the current cell equals its mean value scp0.

subroutine grad_t2d_quad(iw, scp, gxps, gyps, gxxps, gxyps, gyyps)

  use mem_ijtabs, only: itab_w
  use mem_grid,   only: mza, mwa, lpv
  use mem_adv,    only: xy_h

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: scp  (mza,mwa)
  real,    intent(out) :: gxps (mza)
  real,    intent(out) :: gyps (mza)
  real,    intent(out) :: gxxps(mza)
  real,    intent(out) :: gxyps(mza)
  real,    intent(out) :: gyyps(mza)

  integer :: iwn, ivn, k, n
  real    :: sc

  gxps  = 0.
  gyps  = 0.
  gxxps = 0.
  gxyps = 0.
  gyyps = 0.

  ! Loop over neighbors of this W cell

  !dir$ loop count max=7
  do n = 1, itab_w(iw)%npoly
     iwn = itab_w(iw)%iw(n)
     ivn = itab_w(iw)%iv(n)

     do k = lpv(ivn), mza
        sc = scp(k,iwn) - scp(k,iw)

        gxps (k) = gxps (k) + sc * xy_h(1,n,iw)
        gyps (k) = gyps (k) + sc * xy_h(2,n,iw)
        gxxps(k) = gxxps(k) + sc * xy_h(3,n,iw)
        gxyps(k) = gxyps(k) + sc * xy_h(4,n,iw)
        gyyps(k) = gyyps(k) + sc * xy_h(5,n,iw)
     enddo

  enddo

end subroutine grad_t2d_quad

!=======================================================================

! This routine computes gradients of earth-cartesian velocities at each
! cell center in earth-cartesian xe-ye-ze ccordinates using a generalized
! form of the divergence theorum. Cell-centered earth-cartesian velocities
! are averaged to each face to perform the integration. Used for computing
! 3d velocity strain rates for the SGS model.

  subroutine comp_vel_grads_ec(iw, DvxeDxe, DvxeDye, DvxeDze, &
                                   DvyeDxe, DvyeDye, DvyeDze, &
                                   DvzeDxe, DvzeDye, DvzeDze  )
    use mem_ijtabs, only: itab_w
    use mem_basic,  only: vxe, vye, vze
    use mem_grid,   only: mza, lpw, dzit, wnxo2, wnyo2, wnzo2, &
                          vnxo2, vnyo2, vnzo2, arw0i, dnu, lpv
    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: DvxeDxe(mza), DvxeDye(mza), DvxeDze(mza)
    real,    intent(out) :: DvyeDxe(mza), DvyeDye(mza), DvyeDze(mza)
    real,    intent(out) :: DvzeDxe(mza), DvzeDye(mza), DvzeDze(mza)
    integer              :: k, j, iv, iwn
    real                 :: dvxe, dvye, dvze, dsoa, nxds, nyds, nzds

    do k = lpw(iw), mza-1

       if (k == lpw(iw)) then
          dvxe = (vxe(k+1,iw) - vxe(k,iw)) * dzit(k) * 2.
          dvye = (vye(k+1,iw) - vye(k,iw)) * dzit(k) * 2.
          dvze = (vze(k+1,iw) - vze(k,iw)) * dzit(k) * 2.
       else
          dvxe = (vxe(k+1,iw) - vxe(k-1,iw)) * dzit(k)
          dvye = (vye(k+1,iw) - vye(k-1,iw)) * dzit(k)
          dvze = (vze(k+1,iw) - vze(k-1,iw)) * dzit(k)
       endif

       DvxeDxe(k) = dvxe * wnxo2(iw)
       DvxeDye(k) = dvxe * wnyo2(iw)
       DvxeDze(k) = dvxe * wnzo2(iw)

       DvyeDxe(k) = dvye * wnxo2(iw)
       DvyeDye(k) = dvye * wnyo2(iw)
       DvyeDze(k) = dvye * wnzo2(iw)

       DvzeDxe(k) = dvze * wnxo2(iw)
       DvzeDye(k) = dvze * wnyo2(iw)
       DvzeDze(k) = dvze * wnzo2(iw)

    enddo

    DvxeDxe(mza) = DvxeDxe(mza-1)
    DvxeDye(mza) = DvxeDye(mza-1)
    DvxeDze(mza) = DvxeDze(mza-1)

    DvyeDxe(mza) = DvyeDxe(mza-1)
    DvyeDye(mza) = DvyeDye(mza-1)
    DvyeDze(mza) = DvyeDze(mza-1)

    DvzeDxe(mza) = DvzeDxe(mza-1)
    DvzeDye(mza) = DvzeDye(mza-1)
    DvzeDze(mza) = DvzeDze(mza-1)

    do j = 1, itab_w(iw)%npoly
       iv  = itab_w(iw)%iv(j)
       iwn = itab_w(iw)%iw(j)

       dsoa = dnu(iv) * itab_w(iw)%dirv(j) * arw0i(iw)

       nxds = vnxo2(iv) * dsoa
       nyds = vnyo2(iv) * dsoa
       nzds = vnzo2(iv) * dsoa

       do k = lpv(iv), mza
          dvxe = vxe(k,iwn) - vxe(k,iw)
          dvye = vye(k,iwn) - vye(k,iw)
          dvze = vze(k,iwn) - vze(k,iw)

          DvxeDxe(k) = DvxeDxe(k) - dvxe * nxds
          DvxeDye(k) = DvxeDye(k) - dvxe * nyds
          DvxeDze(k) = DvxeDze(k) - dvxe * nzds

          DvyeDxe(k) = DvyeDxe(k) - dvye * nxds
          DvyeDye(k) = DvyeDye(k) - dvye * nyds
          DvyeDze(k) = DvyeDze(k) - dvye * nzds

          DvzeDxe(k) = DvzeDxe(k) - dvze * nxds
          DvzeDye(k) = DvzeDye(k) - dvze * nyds
          DvzeDze(k) = DvzeDze(k) - dvze * nzds
       enddo
    enddo

  end subroutine comp_vel_grads_ec

!=======================================================================

! This routine computes scalar gradients in earth-cartesian xe-ye-ze
! ccordinates using a generalized form of the divergence theorum.
! Cell-centered aarth-cartesian velocities are averaged to each face
! to perform the integration.

  subroutine comp_scp_grads_ec(iw, scp, DsDxe, DsDye, DsDze )

    use mem_ijtabs, only: itab_w
    use mem_grid,   only: mza, mwa, lpw, dzit, wnxo2, wnyo2, wnzo2, &
                          vnxo2, vnyo2, vnzo2, arw0i, dnu, lpv
    implicit none

    integer, intent(in ) :: iw
    real,    intent(in ) :: scp(mza,mwa)
    real,    intent(out) :: DsDxe(mza), DsDye(mza), DsDze(mza)
    integer              :: k, j, iv, iwn
    real                 :: dsc, dsoa, nxds, nyds, nzds

    do k = lpw(iw), mza-1
       dsc = (scp(k+1,iw) - scp(k-1,iw)) * dzit(k)

       DsDxe(k) = dsc * wnxo2(iw)
       DsDye(k) = dsc * wnyo2(iw)
       DsDze(k) = dsc * wnzo2(iw)
    enddo

    DsDxe(mza) = DsDxe(mza-1)
    DsDye(mza) = DsDye(mza-1)
    DsDze(mza) = DsDze(mza-1)

    do j = 1, itab_w(iw)%npoly
       iv  = itab_w(iw)%iv(j)
       iwn = itab_w(iw)%iw(j)

       dsoa = dnu(iv) * itab_w(iw)%dirv(j) * arw0i(iw)

       nxds = vnxo2(iv) * dsoa
       nyds = vnyo2(iv) * dsoa
       nzds = vnzo2(iv) * dsoa

       do k = lpv(iv), mza
          dsc = scp(k,iwn) - scp(k,iw)

          DsDxe(k) = DsDxe(k) - dsc * nxds
          DsDye(k) = DsDye(k) - dsc * nyds
          DsDze(k) = DsDze(k) - dsc * nzds
       enddo
    enddo

  end subroutine comp_scp_grads_ec

!=========================================================================

! Computes the horizontal laplacian of scalar variables at cell centers
! (W points) on the hexagonal grid

subroutine laplacian2d(iw, scp, delsq)

  use mem_ijtabs, only: itab_w
  use mem_grid,   only: mza, mwa, lpv
  use mem_adv,    only: xx_yy

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: scp  (mza,mwa)
  real,    intent(out) :: delsq(mza)

  integer :: iwn, ivn, k, n
  real    :: sc

  delsq = 0.0

  ! Loop over neighbors of this W cell
  do n = 1, itab_w(iw)%npoly

     iwn = itab_w(iw)%iw(n)
     ivn = itab_w(iw)%iv(n)

     do k = lpv(ivn), mza
        sc = scp(k,iwn) - scp(k,iw)
        delsq(k) = delsq(k) + sc * xx_yy(n,iw)
     enddo

  enddo

end subroutine laplacian2d

!=========================================================================

! Computes the horizontal laplacian of scalar variables at cell vertices
! (M points) on the hexagonal grid

subroutine laplacian_m(im, scm, delsq, lpm)

  use mem_ijtabs, only: itab_m
  use mem_grid,   only: mza, mma, lpv
  use mem_adv,    only: xx_yy_m

  implicit none

  integer, intent(in)           :: im
  real,    intent(in)           :: scm(mza,mma)
  integer, intent(in), optional :: lpm(mma)
  real,    intent(out)          :: delsq(mza)

  integer :: imn, k, ka, n
  real    :: sc

  delsq(:) = 0.0

  ! Loop over neighbors of this M point
  do n = 1, 3

     imn = itab_m(im)%im(n)

     if (present(lpm)) then
        ka = max( lpm(im), lpm(imn) )
     else
        ka = lpv( itab_m(im)%iv(n) )
     endif

     do k = ka, mza
        sc = scm(k,imn) - scm(k,im)
        delsq(k) = delsq(k) + sc * xx_yy_m(n,im)
     enddo

  enddo

end subroutine laplacian_m

end module grad_lib
