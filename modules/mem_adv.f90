module mem_adv

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable :: a_v(:,:,:)

  real, allocatable :: xx0(:), xy0(:), yy0(:), xy_h(:,:,:), xx_yy(:,:), xx_yy_m(:,:)

  real, allocatable :: dssq(:)

  real, allocatable :: xx0_v(:,:), xy0_v(:,:), yy0_v(:,:)

  real, allocatable :: gxps_scp (:,:), gyps_scp (:,:)
  real, allocatable :: gxxps_scp(:,:), gxyps_scp(:,:), gyyps_scp(:,:)

  real, allocatable :: dxps_w(:,:), dyps_w(:,:), dzps_w(:,:), dzzps_w(:,:)

  real, allocatable ::  dxps_v(:,:),  dyps_v(:,:),  dzps_v(:,:)
  real, allocatable :: dxxps_v(:,:), dxyps_v(:,:), dyyps_v(:,:)

  real, allocatable :: dxpsw_v (:,:), dypsw_v (:,:)
  real, allocatable :: dxxpsw_v(:,:), dxypsw_v(:,:), dyypsw_v(:,:)

contains

!==============================================================================

  subroutine alloc_adv()

    use mem_grid,    only: mma, mwa, mva, mza, dzm, dztsqo12, dnu, dniv, arw0i, &
                           xew, yew, zew, glatw, glonw, dniu, dnv, arm0
    use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog, itab_v, jtab_v, jtv_wadj, &
                           itab_m, jtab_m, jtm_prog
    use quadrature,  only: hex_quad
    use misc_coms,   only: mdomain, rinit, iparallel
    use consts_coms, only: r8
    use oname_coms,  only: nl
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
    use obnd,        only: lbcopy_w, lbcopy_m
    use map_proj,    only: ec_ps

    implicit none

    integer  :: k, j, iw, iwn, ivn, np, n, iv, iw1, iw2, im, im1, im2, npoly
    real     :: xw(7), yw(7), at(5,5)
    real     :: xe(7), ye(7), ze(7)
    real(r8) :: fint
    real     :: z1(2), z2(2), az1(2,2), az2(2,2), a_h(5,5)
    real     :: dxps1, dyps1, dxps2, dyps2

!!    real, parameter :: wt1 = 1. / 6.
!!    real, parameter :: wt2 = 2. / 3.

    allocate(gxps_scp(mza,mwa))! ; gxps_scp = rinit
    allocate(gyps_scp(mza,mwa))! ; gyps_scp = rinit

    allocate(dxpsw_v(7,mwa))    !; dxpsw_v = rinit
    allocate(dypsw_v(7,mwa))    !; dypsw_v = rinit

    if (nl%horiz_adv_order == 3) then
       allocate(gxxps_scp(mza,mwa)) !; gxxps_scp = rinit
       allocate(gxyps_scp(mza,mwa)) !; gxyps_scp = rinit
       allocate(gyyps_scp(mza,mwa)) !; gyyps_scp = rinit

       allocate(dxxpsw_v(7,mwa))    !; dxxpsw_v = rinit
       allocate(dxypsw_v(7,mwa))    !; dxypsw_v = rinit
       allocate(dyypsw_v(7,mwa))    !; dyypsw_v = rinit
    endif

!    allocate( dxps_w(mza,mwa)) !;  dxps_w = rinit
!    allocate( dyps_w(mza,mwa)) !;  dyps_w = rinit
!   allocate( dzps_w(mza,mwa)) !;  dzps_w = rinit
!   allocate(dzzps_w(mza,mwa)) !; dzzps_w = rinit

!    allocate(dxps_v(mza,mva)) !; dxps_v = rinit
!    allocate(dyps_v(mza,mva)) !; dyps_v = rinit
!    allocate(dzps_v(mza,mva)) !; dzps_v = rinit

    allocate(xx_yy(7,mwa)) !; xx_yy = rinit

    allocate(xx_yy_m(3,mma)) !; xx_yy_m = rinit

!!    if (nl%horiz_adv_order == 3) then
!       allocate(dxxps_v(mza,mva)) ; dxxps_v = rinit
!       allocate(dxyps_v(mza,mva)) ; dxyps_v = rinit
!       allocate(dyyps_v(mza,mva)) ; dyyps_v = rinit

       allocate(xy_h(5,7,mwa)) !; xy_h = rinit

       allocate(xx0(mwa)) !; xx0 = rinit
       allocate(xy0(mwa)) !; xy0 = rinit
       allocate(yy0(mwa)) !; yy0 = rinit

!!       allocate(dxps_m1(2,mva)) ; dxps_m1 = rinit
!!       allocate(dyps_m1(2,mva)) ; dyps_m1 = rinit
!!
!!       allocate(dxps_m2(2,mva)) ; dxps_m2 = rinit
!!       allocate(dyps_m2(2,mva)) ; dyps_m2 = rinit

       allocate(xx0_v(2,mva)) !; xx0_v = rinit
       allocate(xy0_v(2,mva)) !; xy0_v = rinit
       allocate(yy0_v(2,mva)) !; yy0_v = rinit

       if (nl%iscal_monot > 0) allocate(dssq(mwa))
!!    endif

    ! Arrays for finding the vertical quadratic polynomial representation

    allocate(a_v(mza,2,2))

    do k = 2, mza
       z1 = [ -dzm(k-1), dzm(k) ]
       z2 = z1 * z1 - dztsqo12(k)

       az1(1:2,1) = z1
       az1(1:2,2) = z2
       az2 = matmul(transpose(az1), az1)
       call ludcmp(az2, 2)

       a_v(k,2,1) =  (z2(2) - dzm(k)   * az2(2,1)) * az2(2,2)
       a_v(k,2,2) = -(z2(1) + dzm(k-1) * az2(2,1)) * az2(2,2)

       a_v(k,1,1) = az2(1,1) * (dzm(k  ) - az2(1,2) * a_v(k,2,1))
       a_v(k,1,2) = az2(1,1) * (dzm(k-1) - az2(1,2) * a_v(k,2,2))
    enddo

    a_v(1,:,:) = a_v(2,:,:)

    !$omp parallel do private(iw,n,ivn,npoly)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
       npoly = itab_w(iw)%npoly

       do n = 1, npoly
          ivn = itab_w(iw)%iv(n)

         ! Coefficients for horizontal Laplacian at W point
          xx_yy(n,iw) = dniv(ivn) * dnu(ivn) * arw0i(iw)

          if (itab_w(iw)%dirv(n) < 0.) then
             dxpsw_v(n,iw) = itab_v(ivn)%dxps(1)
             dypsw_v(n,iw) = itab_v(ivn)%dyps(1)
          else
             dxpsw_v(n,iw) = itab_v(ivn)%dxps(2)
             dypsw_v(n,iw) = itab_v(ivn)%dyps(2)
          endif

       enddo

       if (nl%horiz_adv_order == 3 .and. nl%iscal_monot > 0) then
          dssq(iw) = minval( dxpsw_v(1:npoly,iw)**2 + dypsw_v(1:npoly,iw)**2 )
       endif

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_w(svara1=xx_yy)
       call mpi_recv_w(svara1=xx_yy)
    endif
    call lbcopy_w(s1=xx_yy)

    !$omp parallel do private(im,n,ivn)
    do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

       do n = 1, 3
          ivn = itab_m(im)%iv(n)

          ! Coefficients for horizontal Laplacian at M point
          xx_yy_m(n,im) = dnv(ivn) * dniu(ivn) / arm0(im)
       enddo

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_m(rvara1=xx_yy_m)
       call mpi_recv_m(rvara1=xx_yy_m)
    endif
    call lbcopy_m(a1=xx_yy_m)

    ! Coefficients for quadratic interpolating polynomical for advection

    if (nl%horiz_adv_order == 3) then

       !$omp parallel do private(iw,np,n,iwn,xw,yw,at,a_h,fint,xe,ye,ze)
       do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
          np  = itab_w(iw)%npoly

          if (mdomain <= 1) then

             do n = 1, np
                iwn = itab_w(iw)%iw(n)
                xe(n) = xew(iwn)
                ye(n) = yew(iwn)
                ze(n) = zew(iwn)
             enddo

             call ec_ps( xe, ye, ze, glatw(iw), glonw(iw), xw, yw, np )

          else

             do n = 1, np
                iwn   = itab_w(iw)%iw(n)
                xw(n) = xew(iwn) - xew(iw)
                yw(n) = yew(iwn) - yew(iw)
             enddo

          endif

          call hex_quad(iw, 2, fint, fxx)
          xx0(iw) = real(fint) * arw0i(iw)

          call hex_quad(iw, 2, fint, fxy)
          xy0(iw) = real(fint) * arw0i(iw)

          call hex_quad(iw, 2, fint, fyy)
          yy0(iw) = real(fint) * arw0i(iw)

          xy_h(1,1:np,iw) = xw(1:np)
          xy_h(2,1:np,iw) = yw(1:np)
          xy_h(3,1:np,iw) = xw(1:np) * xw(1:np) - xx0(iw)
          xy_h(4,1:np,iw) = xw(1:np) * yw(1:np) - xy0(iw)
          xy_h(5,1:np,iw) = yw(1:np) * yw(1:np) - yy0(iw)

          at = matmul( xy_h(:,1:np,iw), transpose(xy_h(:,1:np,iw)) )
          call ludcmp(at, 5)

          a_h = transpose(at)

          xy_h(5,1:np,iw) = xy_h(5,1:np,iw) * a_h(5,5)
          xy_h(4,1:np,iw) = xy_h(4,1:np,iw) * a_h(4,4)
          xy_h(3,1:np,iw) = xy_h(3,1:np,iw) * a_h(3,3)
          xy_h(2,1:np,iw) = xy_h(2,1:np,iw) * a_h(2,2)
          xy_h(1,1:np,iw) = xy_h(1,1:np,iw) * a_h(1,1)

          a_h(1:4,5) = a_h(1:4,5) * a_h(5,5)

          a_h(5  ,4) = a_h(5  ,4) * a_h(4,4)
          a_h(1:3,4) = a_h(1:3,4) * a_h(4,4)
          a_h(4,5  ) = a_h(4,5  ) / a_h(4,4)

          a_h(4:5,3) = a_h(4:5,3) * a_h(3,3)
          a_h(1:2,3) = a_h(1:2,3) * a_h(3,3)
          a_h(3,4:5) = a_h(3,4:5) / a_h(3,3)

          a_h(3:5,2) = a_h(3:5,2) * a_h(2,2)
          a_h(1,  2) = a_h(1  ,2) * a_h(2,2)
          a_h(2,3:5) = a_h(2,3:5) / a_h(2,2)

          a_h(2:5,1) = a_h(2:5,1) * a_h(1,1)
          a_h(1,2:5) = a_h(1,2:5) / a_h(1,1)

          xy_h(2,1:np,iw) = xy_h(2,1:np,iw) - a_h(1,2) * xy_h(1,1:np,iw)

          xy_h(3,1:np,iw) = xy_h(3,1:np,iw) - a_h(1,3) * xy_h(1,1:np,iw) &
                                            - a_h(2,3) * xy_h(2,1:np,iw)

          xy_h(4,1:np,iw) = xy_h(4,1:np,iw) - a_h(1,4) * xy_h(1,1:np,iw) &
                                            - a_h(2,4) * xy_h(2,1:np,iw) &
                                            - a_h(3,4) * xy_h(3,1:np,iw)

          xy_h(5,1:np,iw) = xy_h(5,1:np,iw) - a_h(1,5) * xy_h(1,1:np,iw) &
                                            - a_h(2,5) * xy_h(2,1:np,iw) &
                                            - a_h(3,5) * xy_h(3,1:np,iw) &
                                            - a_h(4,5) * xy_h(4,1:np,iw)

          xy_h(4,1:np,iw) = xy_h(4,1:np,iw) - a_h(5,4) * xy_h(5,1:np,iw)

          xy_h(3,1:np,iw) = xy_h(3,1:np,iw) - a_h(5,3) * xy_h(5,1:np,iw) &
                                            - a_h(4,3) * xy_h(4,1:np,iw)

          xy_h(2,1:np,iw) = xy_h(2,1:np,iw) - a_h(5,2) * xy_h(5,1:np,iw) &
                                            - a_h(4,2) * xy_h(4,1:np,iw) &
                                            - a_h(3,2) * xy_h(3,1:np,iw)

          xy_h(1,1:np,iw) = xy_h(1,1:np,iw) - a_h(5,1) * xy_h(5,1:np,iw) &
                                            - a_h(4,1) * xy_h(4,1:np,iw) &
                                            - a_h(3,1) * xy_h(3,1:np,iw) &
                                            - a_h(2,1) * xy_h(2,1:np,iw)
       enddo
       !$omp end parallel do

       if (iparallel == 1) then
          call mpi_send_w(r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
          call mpi_recv_w(r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
       endif
       call lbcopy_w(v1=xx0, v2=xy0, v3=yy0)

       !$omp parallel
       !$omp do private(iv,iw1,iw2,im1,im2,dxps1,dyps1,dxps2,dyps2)
       do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)
          iw1 = itab_v(iv)%iw(1)
          iw2 = itab_v(iv)%iw(2)

          im1 = itab_v(iv)%im(1)
          im2 = itab_v(iv)%im(2)

          dxps1 = itab_v(iv)%dxps(1)
          dyps1 = itab_v(iv)%dyps(1)

          dxps2 = itab_v(iv)%dxps(2)
          dyps2 = itab_v(iv)%dyps(2)

!!          if (mdomain <= 1) then
!!
!!             call ec_ps( xem(im1),  yem(im1),  zem(im1),  glatw(iw1), glonw(iw1), &
!!                        dxps_m1(1,iv), dyps_m1(1,iv) )
!!
!!             call ec_ps( xem(im2),  yem(im2),  zem(im2),  glatw(iw1), glonw(iw1), &
!!                        dxps_m2(1,iv), dyps_m2(1,iv) )
!!
!!             call ec_ps( xem(im1),  yem(im1),  zem(im1),  glatw(iw2), glonw(iw2), &
!!                        dxps_m1(2,iv), dyps_m1(2,iv) )
!!
!!             call ec_ps( xem(im2),  yem(im2),  zem(im2),  glatw(iw2), glonw(iw2), &
!!                        dxps_m2(2,iv), dyps_m2(2,iv) )
!!
!!          else
!!
!!             dxps_m1(1,iv) = xem(im1) - xew(iw1)
!!             dyps_m1(1,iv) = yem(im1) - yew(iw1)
!!
!!             dxps_m2(1,iv) = xem(im2) - xew(iw1)
!!             dyps_m2(1,iv) = yem(im2) - yew(iw1)
!!
!!             dxps_m1(2,iv) = xem(im1) - xew(iw2)
!!             dyps_m1(2,iv) = yem(im1) - yew(iw2)
!!
!!             dxps_m2(2,iv) = xem(im2) - xew(iw2)
!!             dyps_m2(2,iv) = yem(im2) - yew(iw2)
!!
!!          endif
!!
!!          xx0_v(1,iv) = (dxps_m1(1,iv) * dxps_m1(1,iv) + dxps_m2(1,iv) * dxps_m2(1,iv)) * wt1 &
!!                      + dxps1 * dxps1 * wt2 - xx0(iw1)
!!
!!          xx0_v(2,iv) = (dxps_m1(2,iv) * dxps_m1(2,iv) + dxps_m2(2,iv) * dxps_m2(2,iv)) * wt1 &
!!                      + dxps2 * dxps2 * wt2 - xx0(iw2)
!!
!!          xy0_v(1,iv) = (dxps_m1(1,iv) * dyps_m1(1,iv) + dxps_m2(1,iv) * dyps_m2(1,iv)) * wt1 &
!!                      + dxps1 * dyps1 * wt2 - xy0(iw1)
!!
!!          xy0_v(2,iv) = (dxps_m1(2,iv) * dyps_m1(2,iv) + dxps_m2(2,iv) * dyps_m2(2,iv)) * wt1 &
!!                      + dxps2 * dyps2 * wt2 - xy0(iw2)
!!
!!          yy0_v(1,iv) = (dyps_m1(1,iv) * dyps_m1(1,iv) + dyps_m2(1,iv) * dyps_m2(1,iv)) * wt1 &
!!                      + dyps1 * dyps1 * wt2 - yy0(iw1)
!!
!!          yy0_v(2,iv) = (dyps_m1(2,iv) * dyps_m1(2,iv) + dyps_m2(2,iv) * dyps_m2(2,iv)) * wt1 &
!!                      + dyps2 * dyps2 * wt2 - yy0(iw2)
!!
          xx0_v(1,iv) = dxps1 * dxps1 - xx0(iw1)

          xx0_v(2,iv) = dxps2 * dxps2 - xx0(iw2)

          xy0_v(1,iv) = dxps1 * dyps1 - xy0(iw1)

          xy0_v(2,iv) = dxps2 * dyps2 - xy0(iw2)

          yy0_v(1,iv) = dyps1 * dyps1 - yy0(iw1)

          yy0_v(2,iv) = dyps2 * dyps2 - yy0(iw2)


       enddo
       !$omp end do

       !$omp do private(iw,n,ivn)
       do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
          do n = 1, itab_w(iw)%npoly
             ivn = itab_w(iw)%iv(n)

             if (itab_w(iw)%dirv(n) < 0.) then
                dxxpsw_v(n,iw) = xx0_v(1,ivn)
                dxypsw_v(n,iw) = xy0_v(1,ivn)
                dyypsw_v(n,iw) = yy0_v(1,ivn)
             else
                dxxpsw_v(n,iw) = xx0_v(2,ivn)
                dxypsw_v(n,iw) = xy0_v(2,ivn)
                dyypsw_v(n,iw) = yy0_v(2,ivn)
             endif

          enddo
       enddo
       !$omp end do
       !$omp end parallel

    endif

  end subroutine alloc_adv

!==============================================================================

  real(r8) function fxx(x, y)

    use consts_coms, only: r8
    implicit none

    real(r8), intent(in)::  x, y
    fxx = x * x
  end function fxx

!==============================================================================

  real(r8) function fxy(x, y)

    use consts_coms, only: r8
    implicit none

    real(r8), intent(in)::  x, y
    fxy = x * y

  end function fxy

!==============================================================================

  real(r8) function fyy(x, y)

    use consts_coms, only: r8
    implicit none

    real(r8), intent(in)::  x, y
    fyy = y * y

  end function fyy

!==============================================================================

  subroutine ludcmp(a, n)
    implicit none

    integer, intent(   in) :: n
    real,    intent(inout) :: a(n,n)
    real                   :: dum, sum
    real                   :: vv(n)
    integer                :: i, j, k

    ! ***************************************************************
    ! * Given an N x N matrix A, this routine replaces it by the LU *
    ! * decomposition of a rowwise permutation of itself. A and N   *
    ! * are input. This routine is used with LUBKSB to solve        *
    ! * linear equations or to invert a matrix.                     *
    ! ***************************************************************

    do i = 1, n
       vv(i) = 1.0 / maxval(abs(a(i,1:n)))
    enddo

    do j = 1, n
       do i = 1, j-1
          sum = a(i,j)
          do k = 1, i-1
             sum = sum - a(i,k)*a(k,j)
          enddo
          a(i,j) = sum
       enddo

       do i = j, n
          sum = a(i,j)
          do k = 1, j-1
             sum = sum - a(i,k)*a(k,j)
          enddo
          a(i,j) = sum
       enddo

       if (j /= n) then
          dum = 1.0 / a(j,j)
          do i=j+1,n
             a(i,j) = a(i,j)*dum
          enddo
       endif

    enddo

    do i = 1, n
       a(i,i) = 1.0 / a(i,i)
    enddo

  end subroutine ludcmp

end module mem_adv
