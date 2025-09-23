module mem_adv

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable :: a_v (:,:,:), a_w (:,:,:)
  real, allocatable :: a_vu(:,:,:), a_wu(:,:,:)

  real, allocatable :: xy_h (:,:,:)
  real, allocatable :: xy_hu(:,:,:)

  real, allocatable ::  xx_yy(:,:), xx_yy_m(:,:)

  real, allocatable :: xx0_v (:,:), xy0_v (:,:), yy0_v (:,:)
  real, allocatable :: xx0_vu(:,:), xy0_vu(:,:), yy0_vu(:,:)

contains

!==============================================================================

  subroutine alloc_adv()

    use mem_grid,    only: mma, mwa, mva, mza, dzm, dzt, dnu, dniv, arw0i, &
                           dniu, dnv, arm0, dzt_bot, dzt_top
    use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog, itab_v, jtab_v, jtv_wadj, &
                           itab_m, jtab_m, jtm_prog
    use quadrature,  only: hex_quad
    use misc_coms,   only: iparallel
    use consts_coms, only: r8
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
    use obnd,        only: lbcopy_w, lbcopy_m
    use map_proj,    only: ec_ps

    implicit none

    integer  :: k, j, iw, ivn, np, n, iv, iw1, iw2, im, im1, im2, npoly
    real     :: xw(7), yw(7), at(5,5)
    real(r8) :: fint
    real     :: z1(2), z2(2), az1(2,2), az2(2,2), a_h(5,5)
    real     :: dxps1, dyps1, dxps2, dyps2
    real     :: a0
    real     :: xx0(mwa), yy0(mwa), xy0(mwa)

    allocate( xx_yy  (7,mwa) )
    allocate( xx_yy_m(3,mma) )

    allocate( xy_h (5,7,mwa) )
    allocate( xy_hu(5,7,mwa) )

    allocate( xx0_v (2,mva), xy0_v (2,mva), yy0_v (2,mva) )
    allocate( xx0_vu(2,mva), xy0_vu(2,mva), yy0_vu(2,mva) )

    ! Arrays for finding the vertical quadratic polynomial representation. These
    ! are bounded so that the integral between -dzt_bot and dzt_top equals 0.

    allocate(a_v(mza,2,2))
    allocate(a_w(mza,2,2))

    do k = 2, mza
       a0 =  (dzt_top(k)**3 + dzt_bot(k)**3) / (3. * dzt(k))

       z1 = [ -dzm(k-1), dzm(k) ]
       z2 = z1 * z1 - a0

       az1(1:2,1) = z1
       az1(1:2,2) = z2
       az2 = matmul(transpose(az1), az1)
       call ludcmp(az2, 2)

       a_v(k,2,1) =  (z2(2) - dzm(k)   * az2(2,1)) * az2(2,2)
       a_v(k,2,2) = -(z2(1) + dzm(k-1) * az2(2,1)) * az2(2,2)

       a_v(k,1,1) = az2(1,1) * (dzm(k  ) - az2(1,2) * a_v(k,2,1))
       a_v(k,1,2) = az2(1,1) * (dzm(k-1) - az2(1,2) * a_v(k,2,2))

       a_w(k,1,1) = a_v(k,2,1) * (dzt_bot(k)**2 - a0) - a_v(k,1,1) * dzt_bot(k)
       a_w(k,2,1) = a_v(k,2,2) * (dzt_bot(k)**2 - a0) - a_v(k,1,2) * dzt_bot(k)

       a_w(k,1,2) = a_v(k,2,1) * (dzt_top(k)**2 - a0) + a_v(k,1,1) * dzt_top(k)
       a_w(k,2,2) = a_v(k,2,2) * (dzt_top(k)**2 - a0) + a_v(k,1,2) * dzt_top(k)
    enddo

    a_v(1,:,:) = a_v(2,:,:)
    a_w(1,:,:) = a_w(2,:,:)

    ! Arrays for finding the vertical quadratic polynomial representation.
    ! There are unbounded.

    allocate(a_vu(mza,2,2))
    allocate(a_wu(mza,2,2))

    do k = 2, mza
       z1 = [ -dzm(k-1), dzm(k) ]
       z2 = z1 * z1

       az1(1:2,1) = z1
       az1(1:2,2) = z2
       az2 = matmul(transpose(az1), az1)
       call ludcmp(az2, 2)

       a_vu(k,2,1) =  (z2(2) - dzm(k)   * az2(2,1)) * az2(2,2)
       a_vu(k,2,2) = -(z2(1) + dzm(k-1) * az2(2,1)) * az2(2,2)

       a_vu(k,1,1) = az2(1,1) * (dzm(k  ) - az2(1,2) * a_vu(k,2,1))
       a_vu(k,1,2) = az2(1,1) * (dzm(k-1) - az2(1,2) * a_vu(k,2,2))

       a_wu(k,1,1) = a_vu(k,2,1) * dzt_bot(k)**2 - a_vu(k,1,1) * dzt_bot(k)
       a_wu(k,2,1) = a_vu(k,2,2) * dzt_bot(k)**2 - a_vu(k,1,2) * dzt_bot(k)

       a_wu(k,1,2) = a_vu(k,2,1) * dzt_top(k)**2 + a_vu(k,1,1) * dzt_top(k)
       a_wu(k,2,2) = a_vu(k,2,2) * dzt_top(k)**2 + a_vu(k,1,2) * dzt_top(k)
    enddo

    a_vu(1,:,:) = a_vu(2,:,:)
    a_wu(1,:,:) = a_wu(2,:,:)

    ! Coefficients for horizontal Laplacian at T/W points

    !$omp parallel do private(iw,n,ivn,npoly)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
       npoly = itab_w(iw)%npoly

       do n = 1, npoly
          ivn = itab_w(iw)%iv(n)
          xx_yy(n,iw) = dniv(ivn) * dnu(ivn) * arw0i(iw)
       enddo

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_w(svara1=xx_yy)
       call mpi_recv_w(svara1=xx_yy)
    endif
    call lbcopy_w(s1=xx_yy)

    ! Coefficients for horizontal Laplacian at M/P points

    !$omp parallel do private(im,n,ivn)
    do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

       do n = 1, 3
          ivn = itab_m(im)%iv(n)
          xx_yy_m(n,im) = dnv(ivn) * dniu(ivn) / arm0(im)
       enddo

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_m(rvara1=xx_yy_m)
       call mpi_recv_m(rvara1=xx_yy_m)
    endif
    call lbcopy_m(a1=xx_yy_m)

    ! Coefficients for quadratic interpolating polynomical for advection.
    ! This first set bounds the integral over each W cell.

    !$omp parallel do private(iw,np,n,iv,xw,yw,at,a_h,fint)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
       np  = itab_w(iw)%npoly

       do n = 1, np
          iv = itab_w(iw)%iv(n)

          if (itab_w(iw)%dirv(n) < 0.) then
             xw(n) = itab_v(iv)%dxps(1) * 2.
             yw(n) = itab_v(iv)%dyps(1) * 2.
          else
             xw(n) = itab_v(iv)%dxps(2) * 2.
             yw(n) = itab_v(iv)%dyps(2) * 2.
          endif
       enddo

       ! first compute bounded coefficients

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

       ! now compute unounded coefficients

       xy_hu(1,1:np,iw) = xw(1:np)
       xy_hu(2,1:np,iw) = yw(1:np)
       xy_hu(3,1:np,iw) = xw(1:np) * xw(1:np)
       xy_hu(4,1:np,iw) = xw(1:np) * yw(1:np)
       xy_hu(5,1:np,iw) = yw(1:np) * yw(1:np)

       at = matmul( xy_hu(:,1:np,iw), transpose(xy_hu(:,1:np,iw)) )
       call ludcmp(at, 5)

       a_h = transpose(at)

       xy_hu(5,1:np,iw) = xy_hu(5,1:np,iw) * a_h(5,5)
       xy_hu(4,1:np,iw) = xy_hu(4,1:np,iw) * a_h(4,4)
       xy_hu(3,1:np,iw) = xy_hu(3,1:np,iw) * a_h(3,3)
       xy_hu(2,1:np,iw) = xy_hu(2,1:np,iw) * a_h(2,2)
       xy_hu(1,1:np,iw) = xy_hu(1,1:np,iw) * a_h(1,1)

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

       xy_hu(2,1:np,iw) = xy_hu(2,1:np,iw) - a_h(1,2) * xy_hu(1,1:np,iw)

       xy_hu(3,1:np,iw) = xy_hu(3,1:np,iw) - a_h(1,3) * xy_hu(1,1:np,iw) &
                                           - a_h(2,3) * xy_hu(2,1:np,iw)

       xy_hu(4,1:np,iw) = xy_hu(4,1:np,iw) - a_h(1,4) * xy_hu(1,1:np,iw) &
                                           - a_h(2,4) * xy_hu(2,1:np,iw) &
                                           - a_h(3,4) * xy_hu(3,1:np,iw)

       xy_hu(5,1:np,iw) = xy_hu(5,1:np,iw) - a_h(1,5) * xy_hu(1,1:np,iw) &
                                           - a_h(2,5) * xy_hu(2,1:np,iw) &
                                           - a_h(3,5) * xy_hu(3,1:np,iw) &
                                           - a_h(4,5) * xy_hu(4,1:np,iw)

       xy_hu(4,1:np,iw) = xy_hu(4,1:np,iw) - a_h(5,4) * xy_hu(5,1:np,iw)

       xy_hu(3,1:np,iw) = xy_hu(3,1:np,iw) - a_h(5,3) * xy_hu(5,1:np,iw) &
                                           - a_h(4,3) * xy_hu(4,1:np,iw)

       xy_hu(2,1:np,iw) = xy_hu(2,1:np,iw) - a_h(5,2) * xy_hu(5,1:np,iw) &
                                           - a_h(4,2) * xy_hu(4,1:np,iw) &
                                           - a_h(3,2) * xy_hu(3,1:np,iw)

       xy_hu(1,1:np,iw) = xy_hu(1,1:np,iw) - a_h(5,1) * xy_hu(5,1:np,iw) &
                                           - a_h(4,1) * xy_hu(4,1:np,iw) &
                                           - a_h(3,1) * xy_hu(3,1:np,iw) &
                                           - a_h(2,1) * xy_hu(2,1:np,iw)
    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_w(r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
       call mpi_recv_w(r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
    endif
    call lbcopy_w(v1=xx0, v2=xy0, v3=yy0)

    !$omp parallel do private(iv,iw1,iw2,im1,im2,dxps1,dyps1,dxps2,dyps2)
    do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)
       iw1 = itab_v(iv)%iw(1)
       iw2 = itab_v(iv)%iw(2)

       im1 = itab_v(iv)%im(1)
       im2 = itab_v(iv)%im(2)

       dxps1 = itab_v(iv)%dxps(1)
       dyps1 = itab_v(iv)%dyps(1)

       dxps2 = itab_v(iv)%dxps(2)
       dyps2 = itab_v(iv)%dyps(2)

       xx0_v(1,iv) = dxps1 * dxps1 - xx0(iw1)
       xx0_v(2,iv) = dxps2 * dxps2 - xx0(iw2)

       xx0_vu(1,iv) = dxps1 * dxps1
       xx0_vu(2,iv) = dxps2 * dxps2

       xy0_v(1,iv) = dxps1 * dyps1 - xy0(iw1)
       xy0_v(2,iv) = dxps2 * dyps2 - xy0(iw2)

       xy0_vu(1,iv) = dxps1 * dyps1
       xy0_vu(2,iv) = dxps2 * dyps2

       yy0_v(1,iv) = dyps1 * dyps1 - yy0(iw1)
       yy0_v(2,iv) = dyps2 * dyps2 - yy0(iw2)

       yy0_vu(1,iv) = dyps1 * dyps1
       yy0_vu(2,iv) = dyps2 * dyps2
    enddo
    !$omp end parallel do


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

!==============================================================================

end module mem_adv
