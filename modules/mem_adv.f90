module mem_adv

  real, allocatable :: a_v(:,:,:)

  real, allocatable :: xx0(:), xy0(:), yy0(:), xy_h(:,:,:), xx_yy(:,:)

  real, allocatable :: dxps_m1(:,:), dxps_m2(:,:)
  real, allocatable :: dyps_m1(:,:), dyps_m2(:,:)

  real, allocatable :: xx0_v(:,:), xy0_v(:,:), yy0_v(:,:)

  real, allocatable :: gxps_scp (:,:), gyps_scp (:,:), gzps_scp (:,:)
  real, allocatable :: gxxps_scp(:,:), gxyps_scp(:,:), gyyps_scp(:,:), gzzps_scp(:,:)

  real, allocatable :: gxps_vxe (:,:), gyps_vxe (:,:), gzps_vxe (:,:)
  real, allocatable :: gxxps_vxe(:,:), gxyps_vxe(:,:), gyyps_vxe(:,:), gzzps_vxe(:,:)

  real, allocatable :: gxps_vye (:,:), gyps_vye (:,:), gzps_vye (:,:)
  real, allocatable :: gxxps_vye(:,:), gxyps_vye(:,:), gyyps_vye(:,:), gzzps_vye(:,:)

  real, allocatable :: gxps_vze (:,:), gyps_vze (:,:), gzps_vze (:,:)
  real, allocatable :: gxxps_vze(:,:), gxyps_vze(:,:), gyyps_vze(:,:), gzzps_vze(:,:)

  real, allocatable :: dxps_w(:,:), dyps_w(:,:), dzps_w(:,:), dzzps_w(:,:)

  real, allocatable ::  dxps_v(:,:),  dyps_v(:,:),  dzps_v(:,:)
  real, allocatable :: dxxps_v(:,:), dxyps_v(:,:), dyyps_v(:,:)

contains

  subroutine alloc_adv()

    use mem_grid
    use mem_ijtabs
    use quadrature
    use misc_coms,   only: mdomain, rinit, iparallel
    use consts_coms, only: r8
    use oname_coms,  only: nl
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w
    use obnd,        only: lbcopy_w

    implicit none

    integer  :: k, j, iw, iwn, ivn, np, n, iv, iw1, iw2, im1, im2
    real     :: xw(7), yw(7), at(5,5)
    real(r8) :: fint
    real     :: z1(2), z2(2), az1(2,2), az2(2,2), a_h(5,5)

    allocate(gxps_scp(mza,mwa)) ; gxps_scp = rinit
    allocate(gyps_scp(mza,mwa)) ; gyps_scp = rinit
    allocate(gzps_scp(mza,mwa)) ; gzps_scp = rinit

    allocate(gxps_vxe(mza,mwa)) ; gxps_vxe = rinit
    allocate(gyps_vxe(mza,mwa)) ; gyps_vxe = rinit
    allocate(gzps_vxe(mza,mwa)) ; gzps_vxe = rinit

    allocate(gxps_vye(mza,mwa)) ; gxps_vye = rinit
    allocate(gyps_vye(mza,mwa)) ; gyps_vye = rinit
    allocate(gzps_vye(mza,mwa)) ; gzps_vye = rinit

    allocate(gxps_vze(mza,mwa)) ; gxps_vze = rinit
    allocate(gyps_vze(mza,mwa)) ; gyps_vze = rinit
    allocate(gzps_vze(mza,mwa)) ; gzps_vze = rinit

    if (nl%horiz_adv_order == 3) then
       allocate(gxxps_scp(mza,mwa)) ; gxxps_scp = rinit
       allocate(gxyps_scp(mza,mwa)) ; gxyps_scp = rinit
       allocate(gyyps_scp(mza,mwa)) ; gyyps_scp = rinit

       allocate(gxxps_vxe(mza,mwa)) ; gxxps_vxe = rinit
       allocate(gxyps_vxe(mza,mwa)) ; gxyps_vxe = rinit
       allocate(gyyps_vxe(mza,mwa)) ; gyyps_vxe = rinit

       allocate(gxxps_vye(mza,mwa)) ; gxxps_vye = rinit
       allocate(gxyps_vye(mza,mwa)) ; gxyps_vye = rinit
       allocate(gyyps_vye(mza,mwa)) ; gyyps_vye = rinit

       allocate(gxxps_vze(mza,mwa)) ; gxxps_vze = rinit
       allocate(gxyps_vze(mza,mwa)) ; gxyps_vze = rinit
       allocate(gyyps_vze(mza,mwa)) ; gyyps_vze = rinit
    endif

    allocate(gzzps_scp(mza,mwa)) ; gzzps_scp = rinit
    allocate(gzzps_vxe(mza,mwa)) ; gzzps_vxe = rinit
    allocate(gzzps_vye(mza,mwa)) ; gzzps_vye = rinit
    allocate(gzzps_vze(mza,mwa)) ; gzzps_vze = rinit

    allocate( dxps_w(mza,mwa)) ;  dxps_w = rinit
    allocate( dyps_w(mza,mwa)) ;  dyps_w = rinit
    allocate( dzps_w(mza,mwa)) ;  dzps_w = rinit
    allocate(dzzps_w(mza,mwa)) ; dzzps_w = rinit

    allocate(dxps_v(mza,mva)) ; dxps_v = rinit
    allocate(dyps_v(mza,mva)) ; dyps_v = rinit
    allocate(dzps_v(mza,mva)) ; dzps_v = rinit

    allocate(xx_yy(7,mwa)) ; xx_yy = rinit

    if (nl%horiz_adv_order == 3) then
       allocate(dxxps_v(mza,mva)) ; dxxps_v = rinit
       allocate(dxyps_v(mza,mva)) ; dxyps_v = rinit
       allocate(dyyps_v(mza,mva)) ; dyyps_v = rinit

       allocate(xy_h(5,7,mwa)) ; xy_h = rinit

       allocate(xx0(mwa)) ; xx0 = rinit
       allocate(xy0(mwa)) ; xy0 = rinit
       allocate(yy0(mwa)) ; yy0 = rinit

       allocate(dxps_m1(2,mva)) ; dxps_m1 = rinit
       allocate(dyps_m1(2,mva)) ; dyps_m1 = rinit

       allocate(dxps_m2(2,mva)) ; dxps_m2 = rinit
       allocate(dyps_m2(2,mva)) ; dyps_m2 = rinit

       allocate(xx0_v(2,mva)) ; xx0_v = rinit
       allocate(xy0_v(2,mva)) ; xy0_v = rinit
       allocate(yy0_v(2,mva)) ; yy0_v = rinit
    endif

    ! Arrays for finding the vertical quadratic polynomial representation

    allocate(a_v(mza,2,2))

    do k = 2, mza
       z1 = (/ -dzm(k-1), dzm(k) /)
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

    ! Coefficients for horizontal Laplacian

    !$omp parallel do private(iw,n,ivn)
    do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
       do n = 1, itab_w(iw)%npoly
          ivn = itab_w(iw)%iv(n)
          xx_yy(n,iw) = dniv(ivn) * dnu(ivn) * arw0i(iw)
       enddo
    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_w(1, svara1=xx_yy)
       call mpi_recv_w(1, svara1=xx_yy)
    endif
    call lbcopy_w(1, s1=xx_yy)

    ! Coefficients for quadratic interpolating polynomical for advection

    if (nl%horiz_adv_order == 3) then

       !$omp parallel do private(iw,np,n,iwn,xw,yw,at,a_h,fint)
       do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
          np  = itab_w(iw)%npoly

          do n = 1, np
             iwn = itab_w(iw)%iw(n)
             if (mdomain <= 1) then
                call e_ps( xew(iwn), yew(iwn), zew(iwn), glatw(iw), glonw(iw), &
                           xw(n), yw(n) )
             else
                xw(n) = xew(iwn) - xew(iw)
                yw(n) = yew(iwn) - yew(iw)
             endif
          enddo

          call hex_quad(iw, 2, fint, fxx)
          xx0(iw) = fint / arw0(iw)

          call hex_quad(iw, 2, fint, fxy)
          xy0(iw) = fint / arw0(iw)

          call hex_quad(iw, 2, fint, fyy)
          yy0(iw) = fint / arw0(iw)

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
          call mpi_send_w(1, r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
          call mpi_recv_w(1, r1dvara1=xx0, r1dvara2=xy0, r1dvara3=yy0)
       endif
       call lbcopy_w(1, v1=xx0, v2=xy0, v3=yy0)

       !$omp parallel do private(iv,iw1,iw2,im1,im2)
       do j = 1,jtab_v(jtv_wadj)%jend(1); iv = jtab_v(jtv_wadj)%iv(j)
          iw1 = itab_v(iv)%iw(1)
          iw2 = itab_v(iv)%iw(2)

          im1 = itab_v(iv)%im(1)
          im2 = itab_v(iv)%im(2)

          if (mdomain <= 1) then

             call e_ps( xem(im1),  yem(im1),  zem(im1),  glatw(iw1), glonw(iw1), &
                        dxps_m1(1,iv), dyps_m1(1,iv) )

             call e_ps( xem(im2),  yem(im2),  zem(im2),  glatw(iw1), glonw(iw1), &
                        dxps_m2(1,iv), dyps_m2(1,iv) )

             call e_ps( xem(im1),  yem(im1),  zem(im1),  glatw(iw2), glonw(iw2), &
                        dxps_m1(2,iv), dyps_m1(2,iv) )

             call e_ps( xem(im2),  yem(im2),  zem(im2),  glatw(iw2), glonw(iw2), &
                        dxps_m2(2,iv), dyps_m2(2,iv) )

          else

             dxps_m1(1,iv) = xem(im1) - xew(iw1)
             dyps_m1(1,iv) = yem(im1) - yew(iw1)

             dxps_m2(1,iv) = xem(im2) - xew(iw1)
             dyps_m2(1,iv) = yem(im2) - yew(iw1)

             dxps_m1(2,iv) = xem(im1) - xew(iw2)
             dyps_m1(2,iv) = yem(im1) - yew(iw2)

             dxps_m2(2,iv) = xem(im2) - xew(iw2)
             dyps_m2(2,iv) = yem(im2) - yew(iw2)

          endif

          xx0_v(1,iv) = (dxps_m1(1,iv) * dxps_m1(1,iv) + dxps_m2(1,iv) * dxps_m2(1,iv))/12. - xx0(iw1)
          xx0_v(2,iv) = (dxps_m1(2,iv) * dxps_m1(2,iv) + dxps_m2(2,iv) * dxps_m2(2,iv))/12. - xx0(iw2)

          xy0_v(1,iv) = (dxps_m1(1,iv) * dyps_m1(1,iv) + dxps_m2(1,iv) * dyps_m2(1,iv))/12. - xy0(iw1)
          xy0_v(2,iv) = (dxps_m1(2,iv) * dyps_m1(2,iv) + dxps_m2(2,iv) * dyps_m2(2,iv))/12. - xy0(iw2)

          yy0_v(1,iv) = (dyps_m1(1,iv) * dyps_m1(1,iv) + dyps_m2(1,iv) * dyps_m2(1,iv))/12. - yy0(iw1)
          yy0_v(2,iv) = (dyps_m1(2,iv) * dyps_m1(2,iv) + dyps_m2(2,iv) * dyps_m2(2,iv))/12. - yy0(iw2)

       enddo
       !$omp end parallel do

    endif

  end subroutine alloc_adv


  real(8) function fxx(x, y)
    implicit none
    real(8), intent(in)::  x, y
    fxx = x * x
  end function fxx


  real(8) function fxy(x, y)
    implicit none
    real(8), intent(in)::  x, y
    fxy = x * y
  end function fxy


  real(8) function fyy(x, y)
    implicit none
    real(8), intent(in)::  x, y
    fyy = y * y
  end function fyy


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
