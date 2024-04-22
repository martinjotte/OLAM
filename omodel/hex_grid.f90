subroutine voronoi()

  use mem_ijtabs,   only: mloops, itab_m, itab_v, itab_w, alloc_itabs

  use mem_delaunay, only: itab_md, itab_ud, itab_wd, &
                          xemd, yemd, zemd, nmd, nud, nwd

  use mem_grid,     only: nma, nua, nva, nwa, mma, mua, mva, mwa, &
                          xem, yem, zem, xew, yew, zew, arw0, &
                          alloc_xyzem, alloc_xyzew

  use misc_coms,    only: mdomain
  use consts_coms,  only: pi2, erad

  implicit none

  integer :: im1,im2
  integer :: iw1,iw2,iw3,im,iw
  integer :: iwd,iv,iud,imd,npoly,j
! integer :: iud1,iud2,j1
  real    :: expansion

  ! Interchange grid dimensions

  nma = nwd
  nua = nud
  nva = nud
  nwa = nmd

  mma = nma
  mua = nua
  mva = nva
  mwa = nwa

  ! Allocate Voronoi set of itabs

  call alloc_itabs(nma,nva,nwa)

  ! Allocate XEW,YEW,ZEW arrays, and fill their values from XEMD,YEMD,ZEMD, which
  ! still have the OLD nmad dimension which is the NEW nwa dimension

  call move_alloc(xemd, xew)
  call move_alloc(yemd, yew)
  call move_alloc(zemd, zew)

  allocate(arw0(nwa)) ; arw0 = 0.0

  ! Allocate XEM,YEM,ZEM to NEW nma dimension

  call alloc_xyzem(nma)

  ! Since XEM,YEM,ZEM have just been re-allocated, initialize their values to be
  ! barycenters of Delaunay triangles whose vertices are at XEW,YEW,ZEW

  do iwd = 2, nwd
     im = iwd

  ! Indices of 3 M points surrounding WD point

     if (any(itab_wd(iwd)%im(1:3) < 2)) cycle

     iw1 = itab_wd(iwd)%im(1)
     iw2 = itab_wd(iwd)%im(2)
     iw3 = itab_wd(iwd)%im(3)

     xem(im) = (xew(iw1) + xew(iw2) + xew(iw3)) / 3.
     yem(im) = (yew(iw1) + yew(iw2) + yew(iw3)) / 3.
     zem(im) = (zew(iw1) + zew(iw2) + zew(iw3)) / 3.

     ! If mdomain <= 1, push M point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt(xem(im) ** 2  &
                              + yem(im) ** 2  &
                              + zem(im) ** 2  )

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

     endif

  enddo

  ! Loop over V points

  do iv = 2,nva
     iud = iv

     itab_v(iv)%loop(1:mloops)  = itab_ud(iud)%loop(1:mloops)

     itab_v(iv)%ivp      = itab_ud(iud)%iup
     itab_v(iv)%ivglobe  = iv
     itab_v(iv)%mrlv     = itab_ud(iud)%mrlu

   ! itab_v(iv)%im(1:6)  = itab_ud(iud)%iw(1:6)
   ! itab_v(iv)%iw(1:4)  = itab_ud(iud)%im(1:4)

     itab_v(iv)%im(1:2)  = itab_ud(iud)%iw(1:2)
     itab_v(iv)%iw(1:2)  = itab_ud(iud)%im(1:2)

   ! itab_v(iv)%iv(1:4)  = itab_ud(iud)%iu(1:4)
   ! itab_v(iv)%iv(1:12) = itab_ud(iud)%iu(1:12)

  ! For periodic Cartesian hex domain, compute coordinates for outer M points

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     if (itab_wd(im1)%npoly < 3) then ! itab_m(im1)%npoly not filled yet
        xem(im1) = xew(iw1) + xew(iw2) - xem(im2)
        yem(im1) = yew(iw1) + yew(iw2) - yem(im2)
        zem(im1) = 0.
     elseif (itab_wd(im2)%npoly < 3) then ! itab_m(im2)%npoly not filled yet
        xem(im2) = xew(iw1) + xew(iw2) - xem(im1)
        yem(im2) = yew(iw1) + yew(iw2) - yem(im1)
        zem(im2) = 0.
     endif

     ! Extract information from IMD1 neighbor

   ! imd = itab_ud(iud)%im(1)
   ! npoly = itab_md(imd)%npoly

   ! do j = 1,npoly
   !    j1 = j + 1
   !    if (j == npoly) j1 = 1

   !    iud1 = itab_md(imd)%iu(j)
   !    iud2 = itab_md(imd)%iu(j1)

   !    ! IW(3) and IW(4) neighbors of IV

   !    if (iud2 == iv) then
   !       iw1 = itab_ud(iud1)%im(1)
   !       iw2 = itab_ud(iud1)%im(2)

   !       if (iw1 == imd) then
   !          itab_v(iv)%iw(3) = iw2
   !       else
   !          itab_v(iv)%iw(3) = iw1
   !       endif
   !    endif

   !    if (iud1 == iv) then
   !       iw1 = itab_ud(iud2)%im(1)
   !       iw2 = itab_ud(iud2)%im(2)

   !       if (iw1 == imd) then
   !          itab_v(iv)%iw(4) = iw2
   !       else
   !          itab_v(iv)%iw(4) = iw1
   !       endif
   !    endif

   ! enddo

  enddo

  ! Loop over W points

  do iw = 2,nwa
     imd = iw

     itab_w(iw)%loop(1:mloops) = itab_md(imd)%loop(1:mloops)

     if (mdomain == 0) then
        itab_w(iw)%iwp = iw
     else
        itab_w(iw)%iwp = itab_md(imd)%imp ! Could this be used always?
     endif

     itab_w(iw)%npoly   = itab_md(imd)%npoly
     itab_w(iw)%iwglobe = iw

     itab_w(iw)%mrlw      = itab_md(imd)%mrlm
     itab_w(iw)%mrlw_orig = itab_md(imd)%mrlm_orig
     itab_w(iw)%ngr       = itab_md(imd)%ngr

     npoly = itab_w(iw)%npoly

     ! Loop over IM/IV neighbors of IW

     do j = 1,itab_w(iw)%npoly
        im = itab_md(imd)%iw(j)
        iwd = im
        iv = itab_md(imd)%iu(j)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        itab_w(iw)%im(j) = im
        itab_w(iw)%iv(j) = iv

        if (iw1 == iw) then
           itab_w(iw)%iw(j)   = iw2
           itab_w(iw)%dirv(j) = -1.
        else
           itab_w(iw)%iw(j)   = iw1
           itab_w(iw)%dirv(j) = 1.
        endif

     enddo

  enddo

  ! Loop over M points

  do im = 2,nma
     iwd = im

     itab_m(im)%loop(1:mloops) = itab_wd(iwd)%loop(1:mloops)

     if (mdomain == 0) then
        itab_m(im)%imp = im
     else
        itab_m(im)%imp = itab_wd(iwd)%iwp
     endif

     itab_m(im)%npoly     = itab_wd(iwd)%npoly
     itab_m(im)%imglobe   = im

     itab_m(im)%mrlm      = itab_wd(iwd)%mrlw
     itab_m(im)%ngr       = itab_wd(iwd)%ngr

     itab_m(im)%mrlm_orig = itab_wd(iwd)%mrlw_orig
     itab_m(im)%mrow      = itab_wd(iwd)%mrow

     itab_m(im)%iv(1:3)   = itab_wd(iwd)%iu(1:3)
     itab_m(im)%iw(1:3)   = itab_wd(iwd)%im(1:3)
     itab_m(im)%im(1:3)   = itab_wd(iwd)%iw(1:3)
  enddo

  ! Special: for a cartesian domain, do a lateral boundary copy of MRL

  if (mdomain /= 0) then

     do im = 2, nma
        if (itab_m(im)%imp /= im) then
           itab_m(im)%mrlm = itab_m( itab_m(im)%imp )%mrlm
        endif
     enddo

     do iw = 2, nwa
        if (itab_w(iw)%iwp /= iw) then
           itab_w(iw)%mrlw = itab_w( itab_w(iw)%iwp )%mrlw
        endif
     enddo

     do iv = 2, nva
        if (itab_v(iv)%ivp /= iv) then
           itab_v(iv)%mrlv = itab_v( itab_v(iv)%ivp )%mrlv
        endif
     enddo

  endif

  deallocate(itab_md,itab_ud,itab_wd)

end subroutine voronoi

!===============================================================================

subroutine pcvt()

  ! Iterative procedure for defining centroidal voronoi cells

  use mem_ijtabs,  only: itab_m
  use mem_grid,    only: nma, xem, yem, zem, xew, yew, zew
  use consts_coms, only: erad, eradi
  use misc_coms,   only: mdomain

  implicit none

  integer :: im,iw1,iw2,iw3

  real :: raxis,raxisi,expansion
  real :: sinwlat,coswlat,sinwlon,coswlon
  real :: dxe,dye,dze
  real :: xebc,yebc,zebc
  real :: x1,x2,x3,y1,y2,y3
  real :: dx12,dx13,dx23
  real :: s1,s2,s3
  real :: xcc,ycc

  ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
  ! This establishes W cell as voronoi.

  ! Loop over all M points

  !$omp parallel
  !$omp do private(iw1,iw2,iw3,xebc,yebc,zebc,raxis,raxisi, &
  !$omp            sinwlat,coswlat,sinwlon,coswlon,dxe,dye,dze,x1,y1, &
  !$omp            x2,y2,x3,y3,dx12,dx13,dx23,s1,s2,s3,ycc,xcc)
  do im = 2,nma

     ! Indices of 3 W points surrounding M point

     if (any(itab_m(im)%iw(1:3) < 2)) cycle

     iw1 = itab_m(im)%iw(1)
     iw2 = itab_m(im)%iw(2)
     iw3 = itab_m(im)%iw(3)

     ! These were initialized to be the barycenter of each triangle

     xebc = xem(im)
     yebc = yem(im)
     zebc = zem(im)

     if (mdomain <= 1) then

        ! For global domain, transform from sphere to PS plane

        raxis  = sqrt(xebc ** 2 + yebc ** 2)

        sinwlat = zebc  * eradi
        coswlat = raxis * eradi

        ! For points less than 100 m from Earth's polar axis, make arbitrary
        ! assumption that longitude = 0 deg.  This is just to settle on a PS
        ! planar coordinate system in which to do the algebra.

        if (raxis >= 1.e2) then
           raxisi = 1.0 / raxis
           sinwlon = yebc * raxisi
           coswlon = xebc * raxisi
        else
           sinwlon = 0.
           coswlon = 1.
        endif

        ! Transform 3 W points to PS coordinates

        dxe = xew(iw1) - xebc
        dye = yew(iw1) - yebc
        dze = zew(iw1) - zebc
        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,x1,y1)

        dxe = xew(iw2) - xebc
        dye = yew(iw2) - yebc
        dze = zew(iw2) - zebc
        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,x2,y2)

        dxe = xew(iw3) - xebc
        dye = yew(iw3) - yebc
        dze = zew(iw3) - zebc
        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,x3,y3)

     else

        ! For Cartesian domain, use given planar X,Y coordinates

        x1 = xew(iw1) - xebc
        x2 = xew(iw2) - xebc
        x3 = xew(iw3) - xebc

        y1 = yew(iw1) - yebc
        y2 = yew(iw2) - yebc
        y3 = yew(iw3) - yebc

     endif

     ! Compute intermediate quanties

     dx12 = x2 - x1
     dx13 = x3 - x1
     dx23 = x3 - x2

     s1 = x1**2 + y1**2
     s2 = x2**2 + y2**2
     s3 = x3**2 + y3**2

     ! Algebraic solution for circumcenter Y coordinate

     ycc = .5 * (dx13 * s2 - dx12 * s3 - dx23 * s1) &
              / (dx13 * y2 - dx12 * y3 - dx23 * y1)

     ! Algebraic solution for circumcenter X coordinate

     if (abs(dx12) > abs(dx13)) then
        xcc = (s2 - s1 - ycc * 2. * (y2 - y1)) / (2. * dx12)
     else
        xcc = (s3 - s1 - ycc * 2. * (y3 - y1)) / (2. * dx13)
     endif

     ! For global domain, transform circumcenter from PS to earth coordinates.

     if (mdomain <= 1) then

        call ps_de(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xcc,ycc)

        xem(im) = dxe + xebc
        yem(im) = dye + yebc
        zem(im) = dze + zebc

     else

        xem(im) = xcc + xebc
        yem(im) = ycc + yebc

     endif

  enddo
  !$omp end do nowait

  ! Adjust each M point to the Earth's radius for global domain

  if (mdomain <= 1) then

     !$omp do private(expansion)
     do im = 2, nma

        expansion = erad / sqrt( xem(im) ** 2 &
                               + yem(im) ** 2 &
                               + zem(im) ** 2 )

        xem(im) = xem(im) * expansion
        yem(im) = yem(im) * expansion
        zem(im) = zem(im) * expansion

     enddo
     !$omp end do nowait

  endif
  !$omp end parallel

end subroutine pcvt

!===============================================================================

subroutine grid_geometry_hex()

  use mem_ijtabs,  only: itab_m, itab_v, itab_w
  use mem_grid,    only: nma, nva, nwa, xev, yev, zev, xem, yem, zem, &
                         xew, yew, zew, unx, uny, unz, wnx, wny, wnz, &
                         vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, &
                         arm0, glonm, glatm, glatv, glonv
  use misc_coms,   only: mdomain, nxp
  use consts_coms, only: erad, eradi, piu180, pio2
  use oplot_coms,  only: op
  use consts_coms, only: r8, pio180

  implicit none

  integer            :: im,iv,iw
  integer            :: im1,im2
  integer            :: iw1,iw2
  integer            :: j,npoly
  real               :: expansion
  real               :: dvm1,dvm2
  integer            :: j1,j2
  integer            :: npoly1,npoly2,np
  real               :: xm1,xm2,xv,ym1,ym2,yv,frac
  real               :: xw1,xw2,yw1,yw2
  real               :: xq1, yq1, xq2, yq2, psiz, vsprd
  integer            :: iskip, iwp, ivp, imp
  logical            :: dops
  real               :: quarter_kite(2,nva)
  character(10)      :: string
  real               :: raxis,raxisi,dxe,dye,dze,arwi
  real               :: sinwlat,coswlat,sinwlon,coswlon

  integer               :: lwork, info
  real(r8)              :: a(3,7), b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7)
  real(r8)              :: vrot_x(7), vrot_y(7)
  real(r8), allocatable :: work(:)
  real(r8)              :: wsize(1), vdotw, vmagi, fact

  ! Loop over all M points

  !$omp parallel

  !$omp do private(raxis)
  do im = 2,nma

  ! Latitude and longitude at M points

     if (mdomain <= 1) then
        raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
        glatm(im) = atan2(zem(im),raxis)   * piu180
        glonm(im) = atan2(yem(im),xem(im)) * piu180
     else
        glatm(im) = 0.  ! want it this way?
        glonm(im) = 0.  ! want it this way?
     endif

     ! Fill global index (replaced later if this run is parallel)

     itab_m(im)%imglobe = im

  enddo
  !$omp end do

  ! Loop over all V points

  !$omp do private(im1,im2,iw1,iw2,expansion,raxis,dvm1,dvm2,frac)
  do iv = 2,nva

  ! Fill global index (replaced later if this run is parallel)

     itab_v(iv)%ivglobe = iv

     ! M-point indices of two end points of V segment

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

     ! W-point indices on either side of V segment

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! V point is midway between W points of Voronoi cells

     xev(iv) = .5 * (xew(iw1) + xew(iw2))
     yev(iv) = .5 * (yew(iw1) + yew(iw2))
     zev(iv) = .5 * (zew(iw1) + zew(iw2))

     ! If mdomain <= 1, push V point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt(xev(iv) ** 2 &
                              + yev(iv) ** 2 &
                              + zev(iv) ** 2 )

        xev(iv) = xev(iv) * expansion
        yev(iv) = yev(iv) * expansion
        zev(iv) = zev(iv) * expansion

     endif

     ! Latitude and longitude at V point

     if (mdomain <= 1) then
        raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis
        glatv(iv) = atan2(zev(iv),raxis)   * piu180
        glonv(iv) = atan2(yev(iv),xev(iv)) * piu180
     else
        glatv(iv) = 0. ! want it this way?
        glonv(iv) = 0. ! want it this way?
     endif

     ! Normal distance across U face
     ! Unit vector components of U face
     !x Convert normal distance across U face to geodesic arc length

     dnu(iv) = sqrt( (xem(im1) - xem(im2))**2 &
                   + (yem(im1) - yem(im2))**2 &
                   + (zem(im1) - zem(im2))**2 )

     unx(iv) = (xem(im2) - xem(im1)) / dnu(iv)
     uny(iv) = (yem(im2) - yem(im1)) / dnu(iv)
     unz(iv) = (zem(im2) - zem(im1)) / dnu(iv)

  !x   dnu(iv) = erad2 * asin(dnu(iv) / erad2)
     dniu(iv) = 1. / dnu(iv)

     ! Normal distance across V face
     ! Unit vector components of V face
     !x Convert normal distance across V face to geodesic arc length

     dnv(iv) = sqrt( (xew(iw1) - xew(iw2))**2 &
                   + (yew(iw1) - yew(iw2))**2 &
                   + (zew(iw1) - zew(iw2))**2 )

     vnx(iv) = (xew(iw2) - xew(iw1)) / dnv(iv)
     vny(iv) = (yew(iw2) - yew(iw1)) / dnv(iv)
     vnz(iv) = (zew(iw2) - zew(iw1)) / dnv(iv)

  !x   dnv(iv) = erad2 * asin(dnv(iv) / erad2)
     dniv(iv) = 1. / dnv(iv)

     ! Skip this U point if iw1 < 2 or iw2 < 2

     if (iw1 < 2 .or. iw2 < 2) cycle

     itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw,itab_w(iw2)%mrlw)

     ! Compute IM1 and IM2 values of quarter kite area,
     ! and add to ARM0 and ARW0 arrays

     dvm1 = sqrt((xev(iv) - xem(im1))**2 &
          +      (yev(iv) - yem(im1))**2 &
          +      (zev(iv) - zem(im1))**2)

     dvm2 = sqrt((xev(iv) - xem(im2))**2 &
          +      (yev(iv) - yem(im2))**2 &
          +      (zev(iv) - zem(im2))**2)

     ! Fractional distance along V edge where intersection with U edge is located

     frac = dvm1 * dniu(iv)

     if (im1 > 1 .and. im2 > 1 .and. (frac < .0001 .or. frac > .9999)) then
        print*, 'Non-intersecting U-V edges detected in grid geometry'
        print*, 'FRAC  = ',frac
        print*, 'IW1 = ', iw1, ' IW2 = ', iw2
        print*, 'IV    = ',iv
        print*, 'GLATV = ',glatv(iv)
        print*, 'GLONV = ',glonv(iv)

        print*, 'dnu(iv),dniu(iv) ',dnu(iv),dniu(iv)

        stop 'STOP U-V edges'
     endif

     quarter_kite(1,iv) = .25 * dvm1 * dnv(iv)
     quarter_kite(2,iv) = .25 * dvm2 * dnv(iv)

  enddo
  !$omp end do
  !$omp end parallel

  !dir$ novector
  do iv = 2, nva
     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     arm0(im1) = arm0(im1) + 2. * quarter_kite(1,iv)
     arm0(im2) = arm0(im2) + 2. * quarter_kite(2,iv)

     arw0(iw1) = arw0(iw1) + quarter_kite(1,iv) + quarter_kite(2,iv)
     arw0(iw2) = arw0(iw2) + quarter_kite(1,iv) + quarter_kite(2,iv)
  enddo

  ! Lateral boundary copy of arw0

  do iw = 2,nwa
     iwp = itab_w(iw)%iwp
     if (iw /= iwp) arw0(iw) = arw0(iwp)
  enddo

  ! Lateral boundary copy of arm0

  do im = 2,nma
     imp = itab_m(im)%imp
     if (im /= imp) arm0(im) = arm0(imp)
  enddo

  !$omp parallel private(a,b,wsize,lwork,work,info)
  !$omp do private(raxis)
  do iw = 2,nwa

     ! Fill global index (replaced later if this run is parallel)

     itab_w(iw)%iwglobe = iw

     ! Fill outward unit vector components and latitude and longitude of W point

     if (mdomain <= 1) then

        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

        glatw(iw) = atan2(zew(iw),raxis)   * piu180
        glonw(iw) = atan2(yew(iw),xew(iw)) * piu180

        wnx(iw) = xew(iw) * eradi
        wny(iw) = yew(iw) * eradi
        wnz(iw) = zew(iw) * eradi

     else

        glatw(iw) = 0. ! want it this way?
        glonw(iw) = 0. ! want it this way?

        wnx(iw) = 0.0
        wny(iw) = 0.0
        wnz(iw) = 1.0

     endif
  enddo
  !$omp end do

  !$omp single
  iw = minloc(arw0(2:),1) + 1
  write(*,*)
  write(*,'(A,f0.4,A)')       " Minimum atmos grid spacing is ", .001*sqrt(arw0(iw)), " km"
  write(*,'(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)

  iw = maxloc(arw0(2:),1) + 1
  write(*,*)
  write(*,'(A,f0.4,A)')       " Maximum atmos grid spacing is ", .001*sqrt(arw0(iw)), " km"
  write(*,'(A,I0,2(A,f0.4))') " at iw = ", iw, ", Lat = ", glatw(iw), ", Lon = ", glonw(iw)
  write(*,*)
  !$omp end single

  !$omp do private(npoly,j2,j1,iv,iw1,iw2,dops,npoly1,npoly2,np, &
  !$omp            im1,xw1,xw2,yw1,yw2,xv,yv,raxis,raxisi,vmagi, &
  !$omp            sinwlat,coswlat,sinwlon,coswlon,dxe,dye,dze,arwi)
  do iw = 2,nwa

     ! Number of polygon edges/vertices

     npoly = itab_w(iw)%npoly

     if (mdomain <= 1) then

        raxis  = sqrt(xew(iw)**2 + yew(iw)**2)

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

     endif

     arwi = 1.0 / arw0(iw)

     ! Loop over all polygon edges

     do j2 = 1,npoly
        j1 = j2 - 1
        if (j2 == 1) j1 = npoly

        iv  = itab_w(iw)%iv(j2)
        iw2 = itab_w(iw)%iw(j2)
        iw1 = itab_w(iw)%iw(j1)

        ! Fractional area of arw0(iw) that is occupied by M and V sectors

        if (itab_v(iv)%iw(1) == iw1) then
           itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(1,iv) * arwi
           itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(2,iv) * arwi
        else
           itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(2,iv) * arwi
           itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(1,iv) * arwi
        endif

        itab_w(iw)%farv(j2) = (quarter_kite(1,iv) + quarter_kite(2,iv)) * arwi

        !----------------------------------------
        ! NEW SECTION JULY 2011
        !----------------------------------------

        ! Special - skip gradient calculation if we are at the periodic
        ! domain border and iw1 and iw2 do not share a common vertex

        if ( mdomain <= 1 .or. iw == itab_w(iw)%iwp ) then

           dops = .true.

        else

           dops = .false.
           npoly1 = itab_w(iw1)%npoly
           npoly2 = itab_w(iw2)%npoly

           do np = 1, npoly1
              im1 = itab_w(iw1)%im(np)
              if (im1 == 1) cycle
              if (any( itab_w(iw2)%im(1:npoly2) == itab_w(iw1)%im(np) )) then
                 dops = .true.
                 exit
              endif
           enddo

        endif

        if (dops) then

           ! Evaluate x,y coordinates of IW1 and IW2 points on polar stereographic plane
           ! tangent at IW

           if (mdomain <= 1) then

              dxe = xew(iw1) - xew(iw)
              dye = yew(iw1) - yew(iw)
              dze = zew(iw1) - zew(iw)
              call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xw1,yw1)

              dxe = xew(iw2) - xew(iw)
              dye = yew(iw2) - yew(iw)
              dze = zew(iw2) - zew(iw)
              call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xw2,yw2)

              dxe = xev(iv) - xew(iw)
              dye = yev(iv) - yew(iw)
              dze = zev(iv) - zew(iw)
              call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xv,yv)

           else

              xw1 = xew(iw1) - xew(iw)
              yw1 = yew(iw1) - yew(iw)
              xw2 = xew(iw2) - xew(iw)
              yw2 = yew(iw2) - yew(iw)
              xv  = xev(iv)  - xew(iw)
              yv  = yev(iv)  - yew(iw)

           endif

           ! Coefficients for eastward and northward components of gradient

           vmagi = 1.0 / (xw1 * yw2 - xw2 * yw1)

           itab_w(iw)%gxps1(j1) =  yw2 * vmagi
           itab_w(iw)%gxps2(j1) = -yw1 * vmagi

           itab_w(iw)%gyps1(j1) = -xw2 * vmagi
           itab_w(iw)%gyps2(j1) =  xw1 * vmagi

           !----------------------------------------

           vmagi = 1.0 / sqrt(xw2*xw2 + yw2*yw2)

           if (itab_w(iw)%dirv(j2) < 0.) then

            ! alpha = atan2(yw2,xw2)   ! VC(iv) direction counterclockwise from east
            ! itab_v(iv)%cosv(1) = cos(alpha)
            ! itab_v(iv)%sinv(1) = sin(alpha)

              itab_v(iv)%cosv(1) = xw2 * vmagi
              itab_v(iv)%sinv(1) = yw2 * vmagi

              itab_v(iv)%dxps(1) = xv
              itab_v(iv)%dyps(1) = yv

           else

            ! alpha = atan2(-yw2,-xw2) ! VC(iv) direction counterclockwise from east
            ! itab_v(iv)%cosv(2) = cos(alpha)
            ! itab_v(iv)%sinv(2) = sin(alpha)

              itab_v(iv)%cosv(2) = -xw2 * vmagi
              itab_v(iv)%sinv(2) = -yw2 * vmagi

              itab_v(iv)%dxps(2) = xv
              itab_v(iv)%dyps(2) = yv

           endif

        endif

        !----------------------------------------
        ! END NEW SECTION JULY 2011
        !----------------------------------------

     enddo

     ! Earth-grid components of rotated polar stereographic easterly
     ! (or cartesian positive x) horizontal unit vector

     if (mdomain <= 1) then
        itab_w(iw)%unx_w = -sinwlon
        itab_w(iw)%uny_w =  coswlon
     else
        itab_w(iw)%unx_w = 1.0
        itab_w(iw)%uny_w = 0.0
     endif

     ! Earth-grid components of rotated polar stereographic northerly
     ! (or cartesian positive y) horizontal unit vector

     if (mdomain <= 1) then
        itab_w(iw)%vnx_w = -sinwlat * coswlon
        itab_w(iw)%vny_w = -sinwlat * sinwlon
        itab_w(iw)%vnz_w =  coswlat
     else
        itab_w(iw)%vnx_w = 0.0
        itab_w(iw)%vny_w = 1.0
        itab_w(iw)%vnz_w = 0.0
     endif

  enddo
  !$omp end do

  ! Loop over all V points

  !$omp do private(ivp,iw1,iw2)
  do iv = 2,nva

  ! Let's not do this section on the boundary cells

     ivp = itab_v(iv)%ivp
     if (mdomain > 1 .and. iv /= ivp) cycle

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! FARW(1) and FARW(2) interpolation coefficients for ARW and VOLT
     ! (taking V control volume to be full DNU(IV) * DNV(IV) rectangle)

     itab_v(iv)%farw(1) = 2. * (quarter_kite(1,iv) + quarter_kite(2,iv)) / arw0(iw1)
     itab_v(iv)%farw(2) = 2. * (quarter_kite(1,iv) + quarter_kite(2,iv)) / arw0(iw2)

     itab_v(iv)%mrlv = max(itab_w(iw1)%mrlw,itab_w(iw2)%mrlw)

  enddo  ! IV
  !$omp end do

  ! Now copy back border iv cell values that we skipped

  if (mdomain > 1) then

     !$omp do private(ivp)
     do iv = 2, nva
        ivp = itab_v(iv)%ivp
        if (iv /= ivp) then

           itab_v(iv)%cosv(:) = itab_v(ivp)%cosv(:)
           itab_v(iv)%sinv(:) = itab_v(ivp)%sinv(:)

           itab_v(iv)%dxps(:) = itab_v(ivp)%dxps(:)
           itab_v(iv)%dyps(:) = itab_v(ivp)%dyps(:)

           itab_v(iv)%farw(:) = itab_v(ivp)%farw(:)
           itab_v(iv)%mrlv    = itab_v(ivp)%mrlv
        endif
     enddo
     !$omp end do

  endif

  ! Scale eastward and northward gradient components by farm

  !$omp do private(j)
  do iw = 2, nwa

     ! The gradient components are not computed at the lateral boundaries
     if (iw /= itab_w(iw)%iwp) cycle

     do j = 1, itab_w(iw)%npoly

        itab_w(iw)%gxps1(j) = itab_w(iw)%gxps1(j) * itab_w(iw)%farm(j)
        itab_w(iw)%gyps1(j) = itab_w(iw)%gyps1(j) * itab_w(iw)%farm(j)

        itab_w(iw)%gxps2(j) = itab_w(iw)%gxps2(j) * itab_w(iw)%farm(j)
        itab_w(iw)%gyps2(j) = itab_w(iw)%gyps2(j) * itab_w(iw)%farm(j)

     enddo
  enddo
  !$omp end do

  ! Coefficients for converting earth-cartesian velocity to V and W

  if (mdomain < 2 .or. mdomain == 5) then

     a = 0.0
     b = 0.0
     call dgels( 'N', 3, 7, 1, a, 3, b, 7, wsize, -1, info )

     lwork = nint(wsize(1)) + 1
     allocate(work(lwork))

     !$omp do private(npoly, fo, j, iv, vdotw, vmagi, fact, &
     !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y)
     do iw = 2, nwa

        npoly = itab_w(iw)%npoly

        ! Default coefficients from Perot
        fo(1:npoly) = 2.0_r8 * itab_w(iw)%farv(1:npoly)

        if (mdomain < 2) then

           do j = 1, npoly
              iv = itab_w(iw)%iv(j)

              ! Compute the components of the V unit normals perpendicular to W

              vdotw = vnx(iv)*wnx(iw) + vny(iv)*wny(iw) + vnz(iv)*wnz(iw)

              vnx_ps(j) = vnx(iv) - vdotw * wnx(iw)
              vny_ps(j) = vny(iv) - vdotw * wny(iw)
              vnz_ps(j) = vnz(iv) - vdotw * wnz(iw)

              ! Normalize these new vectors to unit length

              vmagi = 1._r8 / sqrt( vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2 )

              vnx_ps(j) = vnx_ps(j) * vmagi
              vny_ps(j) = vny_ps(j) * vmagi
              vnz_ps(j) = vnz_ps(j) * vmagi

              ! Rotate these new unit normals to a coordinate system with Z aligned with W

              if (wnz(iw) >= 0.0) then

                 fact = ( wny(iw)*vnx_ps(j) - wnx(iw)*vny_ps(j) ) / ( 1._r8 + wnz(iw) )

                 vrot_x(j) = vnx_ps(j)*wnz(iw) - vnz_ps(j)*wnx(iw) + wny(iw)*fact
                 vrot_y(j) = vny_ps(j)*wnz(iw) - vnz_ps(j)*wny(iw) - wnx(iw)*fact

              else

                 fact = ( wny(iw)*vnx_ps(j) - wnx(iw)*vny_ps(j) ) / ( 1._r8 - wnz(iw) )

                 vrot_x(j) = -vnx_ps(j)*wnz(iw) + vnz_ps(j)*wnx(iw) + wny(iw)*fact
                 vrot_y(j) = -vny_ps(j)*wnz(iw) + vnz_ps(j)*wny(iw) - wnx(iw)*fact

              endif

           enddo

        else

           do j = 1, npoly
              iv = itab_w(iw)%iv(j)

              vnx_ps(j) = vnx(iv)
              vny_ps(j) = vny(iv)
              vnz_ps(j) = 0._r8

              vrot_x(j) = vnx_ps(j)
              vrot_y(j) = vny_ps(j)
           enddo

        endif

        a(1,1:npoly) = vrot_x(1:npoly) * vrot_x(1:npoly)
        a(2,1:npoly) = vrot_y(1:npoly) * vrot_y(1:npoly)
        a(3,1:npoly) = vrot_x(1:npoly) * vrot_y(1:npoly)

        b(1) = 1._r8 - sum( fo(1:npoly) * a(1,1:npoly) )
        b(2) = 1._r8 - sum( fo(1:npoly) * a(2,1:npoly) )
        b(3) =       - sum( fo(1:npoly) * a(3,1:npoly) )

        call dgels( 'N', 3, npoly, 1, a, 3, b, 7, work, lwork, info )

        ! Vector b is now the correction to the coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly)

        if (info == 0 .and. all(b(1:npoly) > 0.03_r8) .and. all(b(1:npoly) < 0.7_r8)) then

            itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly)
            itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly)
            itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly)

        else

           write(*,*) "Problem optimizing vector coefficients for iw = ", iw
           write(*,*) glatw(iw), glonw(iw)
           write(*,*) info
           write(*,*) real(b (1:npoly))
           write(*,*) real(fo(1:npoly))
           write(*,*) "Using default coefficients."

           itab_w(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
           itab_w(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
           itab_w(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

        endif

     enddo
     !omp end do

     deallocate(work)

  else

     !$omp do private(npoly)
     do iw = 2, nwa
        npoly = itab_w(iw)%npoly
        itab_w(iw)%ecvec_vx(1:npoly) = 2.0 * itab_w(iw)%farv(1:npoly) &
                                     * vnx(itab_w(iw)%iv(1:npoly))
        itab_w(iw)%ecvec_vy(1:npoly) = 2.0 * itab_w(iw)%farv(1:npoly) &
                                     * vny(itab_w(iw)%iv(1:npoly))
        itab_w(iw)%ecvec_vz(1:npoly) = 2.0 * itab_w(iw)%farv(1:npoly) &
                                     * vnz(itab_w(iw)%iv(1:npoly))
     enddo
     !$omp end do

  endif

  !$omp end parallel

  ! Plot grid lines

  if (.false.) then

     call o_reopnwk()
     call plotback()
     call oplot_set(1)
     psiz = .035 / real(nxp) ! not good with nested grids
     vsprd = .10 * sqrt(arw0(nwa))

     do iv = 2, nva

        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        call oplot_transform(1,xem(im1),yem(im1),zem(im1),xm1,ym1)
        call oplot_transform(1,xem(im2),yem(im2),zem(im2),xm2,ym2)
        call oplot_transform(1,xev(iv),yev(iv),zev(iv),xv,yv)
        call oplot_transform(1,xew(iw2),yew(iw2),zew(iw2),xw2,yw2)
        call oplot_transform(1,xew(iw1),yew(iw1),zew(iw1),xw1,yw1)

        call trunc_segment(xm1,xm2,ym1,ym2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        call o_frstpt (xq1,yq1)
        call o_vector (xq2,yq2)

        if ( xm1 < op%xmin .or.  &
             xm1 > op%xmax .or.  &
             ym1 < op%ymin .or.  &
             ym1 > op%ymax ) cycle

        write(string,'(I0)') im1
        call o_plchlq (xm1,ym1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') im2
        call o_plchlq (xm2,ym2,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iv
        call o_plchlq (xv,yv,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw1
        call o_plchlq (xw1,yw1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw2
        call o_plchlq (xw2,yw2,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') itab_m(im1)%imp
        call o_plchlq (xm1,ym1+vsprd,trim(adjustl(string)),.5*psiz,0.,0.)

        write(string,'(I0)') itab_m(im2)%imp
        call o_plchlq (xm2,ym2+vsprd,trim(adjustl(string)),.5*psiz,0.,0.)

        write(string,'(I0)') itab_v(iv)%ivp
        call o_plchlq (xv,yv+vsprd,trim(adjustl(string)),.5*psiz,0.,0.)

        write(string,'(I0)') itab_w(iw1)%iwp
        call o_plchlq (xw1,yw1+vsprd,trim(adjustl(string)),.5*psiz,0.,0.)

        write(string,'(I0)') itab_w(iw2)%iwp
        call o_plchlq (xw2,yw2+vsprd,trim(adjustl(string)),.5*psiz,0.,0.)

     enddo  ! IV

     call o_frame()
     call o_clswk()

  endif

end subroutine grid_geometry_hex

!===============================================================================

subroutine ctrlvols_hex()

  use mem_ijtabs,  only: jtab_m, jtab_v, jtab_w, itab_m, itab_v, itab_w, &
                         jtm_grid, jtv_grid, jtv_lbcp, jtw_grid, jtw_lbcp, &
                         jtm_lbcp
  use mem_sfcg,    only: nwsfc, sfcg, itab_wsfc
  use misc_coms,   only: io6, mdomain
  use consts_coms, only: r8
  use mem_grid,    only: nsw_max, nza, nwa, lpm, lpv, lpw, lsw, topm, zm, &
                         dzt, zfact, zfacm2, dnu, arw0, arv, arw, volt, lve2, &
                         nve2_max, dzt_bot
  implicit none

  integer  :: j,iw,iwp,iv,ivp,im1,im2,k,km,im,imp,iw1,iw2
  integer  :: kw,npoly,jv,iv1,iv2,iv3,iwsfc
  real     :: hmin, hmax, amin, facw
  real(r8) :: vmin, vmax
  logical  :: docheck
  integer  :: nopen, js, jj, jasfc, kbot

  integer, allocatable :: kworig(:)
  integer, allocatable :: jsorig(:)

  real, parameter :: afrac_min = 0.30  ! minimum area fraction (compared to fully
                                       ! open cell) below which cell face is closed

  write(io6,*) 'Defining control volume areas'

  arw(:,:) = 0.

  ! Loop over all SURFACE cells

  do iwsfc = 2,nwsfc

     ! Loop over ATM cells that couple to current SURFACE cell.  Add coupling
     ! area to all ARW areas of ATM grid at and above the KW index that is
     ! (provisionally) attached to the coupling area.

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw  = itab_wsfc(iwsfc)%iwatm(j)
        kw  = itab_wsfc(iwsfc)%kwatm(j)

        arw(kw:nza-1,iw) = arw(kw:nza-1,iw) + itab_wsfc(iwsfc)%arc(j)

        itab_w(iw)%jsfc2 = itab_w(iw)%jsfc2 + 1
     enddo
  enddo

  ! Compute the mappings between the atmosphere and surface meshes

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)
     allocate( itab_w(iw)%iwsfc( itab_w(iw)%jsfc2 ) )
     allocate( itab_w(iw)%jasfc( itab_w(iw)%jsfc2 ) )
     itab_w(iw)%jsfc2 = 0
  enddo

  do iwsfc = 2, nwsfc
     do j = 1, itab_wsfc(iwsfc)%nwatm
        iw = itab_wsfc(iwsfc)%iwatm(j)

        itab_w(iw)%jsfc2 = itab_w(iw)%jsfc2 + 1
        itab_w(iw)%iwsfc(  itab_w(iw)%jsfc2  ) = iwsfc
        itab_w(iw)%jasfc(  itab_w(iw)%jsfc2  ) = j
     enddo
  enddo

  ! Loop over all ATM cells and close those whose ARW is below a specified limit

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)
     do k = nza-1, 1, -1
        if (arw(k,iw) < afrac_min * arw0(iw)) then
           arw(1:k,iw) = 0.0
           lpw    (iw) = k + 1
           exit
        endif
     enddo
  enddo

  ! Lateral boundary copy of ARW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
     arw(:,iw) = arw(:,iwp)
     lpw  (iw) = lpw  (iwp)
  enddo

  ! Loop over all ATM grid V edges and get FIRST ESTIMATE of ARV, subject to
  ! subsequent adjustment

  !$omp parallel do private(iv,im1,im2,iw1,iw2,k,hmin,hmax,km,kbot)
  do j = 1,jtab_v(jtv_grid)%jend; iv = jtab_v(jtv_grid)%iv(j)
     im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

     ! lowest possible lpv
     kbot = max( lpw(iw1), lpw(iw2) )

     arv(1:kbot-1,iv) = 0.
     lpv         (iv) = kbot

     if (dnu(iv) < 1.e-6) then

        arv(kbot:nza,iv) = 0.
        lpv(iv) = nza + 1

     elseif (itab_w(iw1)%jsfc2 == 1 .and. itab_w(iw2)%jsfc2 == 1) then

        do k = kbot, nza
           arv(k,iv) = dnu(iv) * dzt(k)
        enddo

     else

        if (topm(im2) > topm(im1)) then
           hmin = topm(im1)
           hmax = topm(im2)
        else
           hmin = topm(im2)
           hmax = topm(im1)
        endif

        do k = nza, kbot, -1
           km = k - 1

           if (zm(k) <= hmin) then

              ! Close remaining V faces if below terrain height

              arv(kbot:k,iv) = 0.
              lpv       (iv) = k + 1
              exit

           elseif (zm(km) >= hmax) then

              arv(k,iv) = dnu(iv) * dzt(k)

           elseif (zm(k) < hmax .and. zm(km) < hmin) then

              arv(k,iv) = dnu(iv) * .5 * (zm(k) - hmin)**2 / (hmax - hmin)

           elseif (zm(k) <  hmax .and. zm(km) >=  hmin) then

              arv(k,iv) = dnu(iv) * dzt(k)  &
                        * (.5 * (zm(k) + zm(km)) - hmin) / (hmax - hmin)

           elseif (zm(k) >= hmax .and. zm(km) < hmin) then

              arv(k,iv) = dnu(iv) * (zm(k) - .5 * (hmax + hmin))

           elseif (zm(k) >= hmax .and. zm(km) >= hmin) then

              arv(k,iv) = dnu(iv)  &
                        * ( dzt(k) - .5 * (hmax - zm(km) )**2 / (hmax - hmin))

           endif

           ! Close remaining V faces if ARV is small

           if (arv(k,iv) < afrac_min * dnu(iv) * dzt(k)) then

              arv(kbot:k,iv) = 0.
              lpv       (iv) = k + 1
              exit

           endif

        enddo ! k

     endif

  enddo ! j,iv
  !$omp end parallel do

  ! Lateral boundary copy of ARV

  do j = 1,jtab_v(jtv_lbcp)%jend; iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp
     arv(:,iv) = arv(:,ivp)
     lpv  (iv) = lpv  (ivp)
  enddo

  ! Option for stability: if cell has only one or two lateral faces open,
  ! close entire cell:

  docheck = .true.
  do while (docheck)

     docheck = .false.
     do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

        npoly = itab_w(iw)%npoly
        kbot  = lpw(iw)

        do k = kbot, nza-1
           nopen = count( arv( k,itab_w(iw)%iv(1:npoly) ) > 1.e-8 )

           if (nopen < 3) then

              docheck   = .true.
              arw(k,iw) = 0.0
              lpw  (iw) = k + 1

              arv(k,itab_w(iw)%iv(1:npoly)) = 0.0
              lpv  (itab_w(iw)%iv(1:npoly)) = max(k+1,lpv(itab_w(iw)%iv(1:npoly)))

           else

              exit

           endif
        enddo

     enddo

     ! Lateral boundary copy of ARV

     do j = 1,jtab_v(jtv_lbcp)%jend; iv = jtab_v(jtv_lbcp)%iv(j)
        ivp = itab_v(iv)%ivp
        arv(:,iv) = arv(:,ivp)
        lpv  (iv) = lpv  (ivp)
     enddo

     ! Lateral boundary copy of ARW

     do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
        iwp = itab_w(iw)%iwp
        arw(:,iw) = arw(:,iwp)
        lpw  (iw) = lpw  (iwp)
     enddo

  enddo

  !=============================================================================
  ! At this point, ATM cells with insufficient top ARW and/or only 1 or 2 ARV
  ! sides open have been completely closed (ARV and ARW have been set to 0)
  !
  ! No further complete closures of ARW or ARV are allowed beyond this point,
  !
  ! Now help maintain CFL stability in shaved cells by keeping the ratio of
  ! open lateral face area to top area (and thus cell volume) consistent
  !=============================================================================

  j = maxval( itab_w(2:nwa)%jsfc2 )

  allocate( kworig( j ) )
  allocate( jsorig( j ) )

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

     do js = 1, itab_w(iw)%jsfc2
        iwsfc = itab_w(iw)%iwsfc(js)
        jasfc = itab_w(iw)%jasfc(js)

        kworig(js) = itab_wsfc(iwsfc)%kwatm(jasfc)
        jsorig(js) = js
     enddo

     ! Sort surface patches by kw level

     if (itab_w(iw)%jsfc2 > 1) then
        call isort2( itab_w(iw)%jsfc2, kworig, jsorig )
     endif

     do k = nza-1, lpw(iw), -1

        if (arw(k,iw) > 0.9 * arw0(iw)) cycle

        ! Each V face contributes to open up its fraction farv of the
        ! current cell's top open area:

        facw = 0.0
        do jv = 1, itab_w(iw)%npoly
           iv = itab_w(iw)%iv(jv)
           facw = facw + itab_w(iw)%farv(jv) * arv(k,iv) / (dnu(iv) * dzt(k))
        enddo

        amin = min( facw, 0.9 ) * arw0(iw)

        ! If arw is less than amin, increase arw by moving surface cell patches
        ! down to this level.

        if (arw(k,iw) < amin) then

           do jj = 1, itab_w(iw)%jsfc2
              js    = jsorig(jj)
              iwsfc = itab_w(iw)%iwsfc(js)
              jasfc = itab_w(iw)%jasfc(js)

              if ( itab_wsfc(iwsfc)%kwatm(jasfc) == k+1 ) then
                 itab_wsfc(iwsfc)%kwatm(jasfc) = k
                 arw(k,iw) = arw(k,iw) + itab_wsfc(iwsfc)%arc(jasfc)
              endif

              if (arw(k,iw) >= amin) exit
           enddo

        endif

     enddo
  enddo

  ! Lateral boundary copy of ARW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
     arw(:,iw) = arw(:,iwp)
  enddo

  !=============================================================================
  ! No further adjustments of ARW or ARV are allowed beyond this point.
  !=============================================================================

  do j = 1,jtab_m(jtm_grid)%jend; im = jtab_m(jtm_grid)%im(j)
     iv1 = itab_m(im)%iv(1)
     iv2 = itab_m(im)%iv(2)
     iv3 = itab_m(im)%iv(3)

     lpm(im) = max( lpv(iv1), lpv(iv2), lpv(iv3) )
  enddo

  ! Lateral boundary copy of LPM

  do j = 1,jtab_m(jtm_lbcp)%jend; im = jtab_m(jtm_lbcp)%im(j)
     imp = itab_m(im)%imp
     lpm(im) = lpm(imp)
  enddo

  allocate(sfcg%dzt_bot(nwsfc))

  ! Loop over all SURFACE cells

  !$omp parallel
  !$omp do private(iwsfc,j,iw,kw)
  do iwsfc = 2,nwsfc

     ! Loop over ATM cells that couple to current SURFACE cell

     do j = 1,itab_wsfc(iwsfc)%nwatm
        iw  = itab_wsfc(iwsfc)%iwatm(j)
        kw  = itab_wsfc(iwsfc)%kwatm(j)

        ! If provisional KW of coupling area is less than LPW of this ATM cell,
        ! reset KW to LPW.

        kw = max( kw, lpw(iw) )
        itab_wsfc(iwsfc)%kwatm(j) = kw

        ! Compute ratios of coupling area to ATM and SFC grid cell areas

        !!!!!!!!!!!! error check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if (arw(kw,iw) - arw(kw-1,iw) < 1.) then
        !   print*, 'small arw dif! ',iw,kw,itab_wsfc(iwsfc)%arc(j),arw(kw,iw),arw(kw-1,iw)
        !endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        itab_wsfc(iwsfc)%arcoariw (j) = itab_wsfc(iwsfc)%arc(j) / arw0(iw)
        itab_wsfc(iwsfc)%arcoarkw (j) = itab_wsfc(iwsfc)%arc(j) / max( arw(kw,iw) - arw(kw-1,iw), 1.)
        itab_wsfc(iwsfc)%arcoarsfc(j) = itab_wsfc(iwsfc)%arc(j) / sfcg%area(iwsfc)

        ! Expand coupling area for spherical geometry, and add change to SFC cell area

        if (mdomain < 2) then
           itab_wsfc(iwsfc)%arc(j) = itab_wsfc(iwsfc)%arc(j) * zfacm2(kw-1)
           sfcg%area(iwsfc) = sfcg%area(iwsfc) &
                            + itab_wsfc(iwsfc)%arc(j) * (zfacm2(kw-1) - 1.0)
        endif

        ! Compute SFC grid cell dzt_bot as weighted average over ATM cells

        sfcg%dzt_bot(iwsfc) = sfcg%dzt_bot(iwsfc) + dzt_bot(kw) * itab_wsfc(iwsfc)%arcoarsfc(j)

     enddo

  enddo
  !$omp end do

  ! Loop over all ATM grid columns

  !$omp do private(iw,k,vmin,vmax,jv,iv,js,iwsfc,jasfc,kw)
  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

     ! special at model top
     volt(nza,iw) = dzt(nza) * arw0(iw)

     ! Loop over vertical levels

     do k = nza-1, lpw(iw), -1

        if (arw(k-1,iw) > 0.999 * arw(k,iw)) then

           volt(k,iw) = dzt(k) * arw(k,iw)

        else

           call volt_from_flux(k,iw)

           ! set limits on volt for CFL stability
           vmin = 0.5 * (arw(k,iw)+arw(k-1,iw)) * dzt(k)
           vmax = arw(k,iw) * dzt(k)
           volt(k,iw) = max( vmin, min(vmax, volt(k,iw) ) )

        endif

     enddo

     volt(1:lpw(iw)-1,iw) = 0._r8

     ! Expand ARW and VOLT with height for spherical geometry

     if (mdomain < 2) then
        do k = lpw(iw), nza
           arw (k,iw) = arw (k,iw) * zfacm2(k)
           volt(k,iw) = volt(k,iw) * zfact (k)**2
        enddo
     endif

     ! Compute number of underground v[xyz]e2 levels in this IW column

     lve2(iw) = 0.0

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)
        lve2(iw) = max(lve2(iw), lpv(iv) - lpw(iw))
     enddo

     ! Compute number of levels that intersect with surface

     lsw(iw) = 0.0

     do js = 1, itab_w(iw)%jsfc2
        iwsfc = itab_w(iw)%iwsfc(js)
        jasfc = itab_w(iw)%jasfc(js)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)
        lsw(iw) = max(lsw(iw), kw - lpw(iw) + 1)
     enddo

  enddo
  !$omp end do
  !$omp end parallel

  ! Lateral boundary copy of ARW, VOLT, LPW, and LSW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp

     arw (:,iw) = arw (:,iwp)
     volt(:,iw) = volt(:,iwp)
     lpw   (iw) = lpw   (iwp)
     lsw   (iw) = lsw   (iwp)
     lve2  (iw) = lve2  (iwp)
  enddo

  nsw_max  = maxval( lsw (2:))
  nve2_max = maxval( lve2(2:))

  ! Expand ARV with height for spherical geometry

  if (mdomain < 2) then
     do j = 1,jtab_v(jtv_grid)%jend; iv = jtab_v(jtv_grid)%iv(j)
        do k = lpv(iv), nza
           arv(k,iv) = arv(k,iv) * zfact(k)
        enddo
     enddo
  endif

  ! Lateral boundary copy of ARV

  do j = 1,jtab_v(jtv_lbcp)%jend; iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp

     arv(:,iv) = arv(:,ivp)
     lpv  (iv) = lpv  (ivp)
  enddo

end subroutine ctrlvols_hex

!===============================================================================

subroutine ctrlvols_hex_nosfcg()

  use mem_ijtabs,  only: jtab_m, jtab_v, jtab_w, itab_m, itab_v, itab_w, &
                         jtm_grid, jtv_grid, jtv_lbcp, jtw_grid, jtw_lbcp
  use misc_coms,   only: io6, mdomain
  use mem_grid,    only: nsw_max, nza, nma, nva, nwa, lpm, lpv, lpw, lsw, &
                         zm, dzt, zfact, zfacm2, dnu, &
                         arw0, arv, arw, volt, lve2, nve2_max, &
                         glatw, glonw, xew, yew, zew, topw, &
                         glatm, glonm, xem, yem, zem, topm
  implicit none

  integer  :: j,iw,iwp,iv,ivp,im1,im2,k,km,im,iw1,iw2,iw3
  integer :: kw,npoly,jv,iv1,iv2,iv3, ipat
  real    :: hmin,hmax,topmin, topmax
  logical :: docheck

  integer :: npats(7), kw_pats(nza,7)
  real    :: area_pats(nza,7)

  real :: arw0_check

  ! This subroutine is called instead of ctrlvols_hex when OLAM is not using
  ! a surface grid.  It must therefore define control volume parameters directly
  ! from topography info, which it first obtains from a call to topo_init.

  ! Initialize TOPW and TOPM by calling topo_init

  call topo_init(nwa,topw,glatw,glonw,xew,yew,zew)
  call topo_init(nma,topm,glatm,glonm,xem,yem,zem)

  ! Prevent TOPM from being lower than lowest model level zm(1)

  topm(2:nma) = max(topm(2:nma),zm(1))

  write(io6,*) 'Defining control volume areas'

  arw(:,:) = 0.

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)
     npoly = itab_w(iw)%npoly

     ! Prevent TOPW from being outside the range of surrounding TOPM values

     topmax = maxval(topm(itab_w(iw)%im(1:npoly)))
     topmin = minval(topm(itab_w(iw)%im(1:npoly)))

     topw(iw) = max(topmin,min(topmax,topw(iw)))

     ! Evaluate levels where topographic surface intersects this IW column

     call atm_grid_topocut(iw, npoly, npats(1:npoly), kw_pats(:,1:npoly), &
                                                    area_pats(:,1:npoly))

     arw0_check = 0.

     do jv = 1,npoly
        do ipat = 1,npats(jv)
           kw = kw_pats(ipat,jv)

           arw(kw:nza-1,iw) = arw(kw:nza-1,iw) + area_pats(ipat,jv)
           arw0_check       = arw0_check       + area_pats(ipat,jv)

           lsw(iw) = max(lsw(iw), kw - lpw(iw) + 1)
        enddo
     enddo

     ! Check for equality between arw0 and arw0_check

     if (abs(arw0(iw) - arw0_check) > 1.e-6*min(arw0(iw),arw0_check)) then
        print*, 'arw0_check ',iw,arw0(iw),arw0_check,arw0(iw)/arw0_check
     endif

  enddo

  ! Loop over all ATM cells and close those whose ARW is below a specified limit

GOTO 1

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)
     do k = nza-1, 2, -1
        if (arw(k,iw) < 0.2 * arw0(iw)) then
           arw (1:k,iw) = 0.0
           exit
        endif
     enddo
  enddo

1 CONTINUE

  ! Lateral boundary copy of ARW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
     arw (:,iw) = arw (:,iwp)
  enddo

  ! Loop over all ATM grid V edges and get FIRST ESTIMATE of ARV, subject to
  ! subsequent adjustment

  !$omp parallel do private(iv,im1,im2,iw1,iw2,k,hmin,hmax,km)
  do j = 1,jtab_v(jtv_grid)%jend; iv = jtab_v(jtv_grid)%iv(j)
     im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

     arv(1,iv) = 0.

     if (dnu(iv) < 1.e-6) then

        do k = 2,nza
           arv(k,iv) = 0.
        enddo

     else

        if (topm(im2) > topm(im1)) then
           hmin = topm(im1)
           hmax = topm(im2)
        else
           hmin = topm(im2)
           hmax = topm(im1)
        endif

        do k = nza,2,-1
           km = k - 1

           if (k < nza .and. (arw(k,iw1) <= 1.e-9 .or. arw(k,iw2) <= 1.e-9)) then

              ! close V if top of either T neighbor is completely closed
              arv(k,iv) = 0.

           elseif (zm(k) <= hmin) then

              ! close V if below terrain height
              arv(k,iv) = 0.

           elseif (zm(km) >= hmax) then

              arv(k,iv) = dnu(iv) * dzt(k)

           elseif (zm(k) < hmax .and. zm(km) < hmin) then

              arv(k,iv) = dnu(iv) * .5 * (zm(k) - hmin)**2 / (hmax - hmin)

           elseif (zm(k) <  hmax .and. zm(km) >=  hmin) then

              arv(k,iv) = dnu(iv) * dzt(k)  &
                 * (.5 * (zm(k) + zm(km)) - hmin) / (hmax - hmin)

           elseif (zm(k) >= hmax .and. zm(km) < hmin) then

              arv(k,iv) = dnu(iv) * (zm(k) - .5 * (hmax + hmin))

           elseif (zm(k) >= hmax .and. zm(km) >=  hmin) then

              arv(k,iv) = dnu(iv)  &
                 * (dzt(k) - .5 * (hmax - zm(km)) ** 2 / (hmax - hmin))

           else

              write(io6,*) 'arv option not reached ',k,iv,j,  &
                 zm(k),zm(km),hmax,hmin
              stop 'stop arv defn'

           endif

GOTO 2

           if (arv(k,iv) < 0.2 * dnu(iv) * dzt(k)) then
              arv(1:k,iv) = 0.0
              exit
           endif

2 CONTINUE

        enddo ! k

     endif

  enddo ! j,iv
  !$omp end parallel do

  ! option for stability: if hexagon has only one lateral face open, close entire cell

  docheck = .true.
  do while (docheck)

     docheck = .false.
     do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

        npoly = itab_w(iw)%npoly

        do k = nza-1, 2, -1

           if ( count( arv( k,itab_w(iw)%iv(1:npoly) ) < 1.e-8 ) == npoly-1 ) then
              docheck = .true.
              arw (1:k,iw) = 0.0
              arv (1:k,itab_w(iw)%iv(1:npoly)) = 0.0
              exit
           endif
        enddo
     enddo

  enddo

  ! Lateral boundary copy of ARV

  do j = 1,jtab_v(jtv_lbcp)%jend; iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp
     arv(:,iv) = arv(:,ivp)
  enddo

  ! Lateral boundary copy of ARW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
     arw (:,iw) = arw (:,iwp)
  enddo

!===============================================================================
! At this point, ATM cells with insufficient top ARW and/or only one ARV
! side open have been completely closed (ARV and ARW have been set to 0)
!===============================================================================
! No further complete closures of ARW or ARV are allowed beyond this point,
! although ARW will be recomputed based on closures that have already been done
! {also, ARV may be adjusted downward if it is too large for an adjacent ARW??}
!===============================================================================

  ! At each ATM grid M point, compute upper bound for LPM and for LPW and LPV
  ! of neighbor points using neighbor ARVs and ARWs.  This is final value for
  ! LPM, LPV, and LPW.

  lpm(2:nma) = nza
  lpv(2:nva) = nza
  lpw(2:nwa) = nza

  do j = 1,jtab_m(jtm_grid)%jend; im = jtab_m(jtm_grid)%im(j)

     iv1 = itab_m(im)%iv(1)
     iv2 = itab_m(im)%iv(2)
     iv3 = itab_m(im)%iv(3)

     iw1 = itab_m(im)%iw(1)
     iw2 = itab_m(im)%iw(2)
     iw3 = itab_m(im)%iw(3)

     ! Loop over vertical levels from top to bottom

     do k = nza,2,-1
        if (arv(k,iv1) < .005 * dnu(iv1) * dzt(k)) exit
        if (arv(k,iv2) < .005 * dnu(iv2) * dzt(k)) exit
        if (arv(k,iv3) < .005 * dnu(iv3) * dzt(k)) exit

        if (k < nza) then
           if (arw(k,iw1) < .001 * arw0(iw1)) exit
           if (arw(k,iw2) < .001 * arw0(iw2)) exit
           if (arw(k,iw3) < .001 * arw0(iw3)) exit
        endif

        lpm(im) = k

        lpv(iv1) = min(lpv(iv1),k)
        lpv(iv2) = min(lpv(iv2),k)
        lpv(iv3) = min(lpv(iv3),k)

        lpw(iw1) = min(lpw(iw1),k)
        lpw(iw2) = min(lpw(iw2),k)
        lpw(iw3) = min(lpw(iw3),k)
     enddo

  enddo

  lsw(:)  = 1
  nsw_max = 1

  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

     ! Loop over vertical levels from top to bottom

     do k = nza-1,2,-1

        ! Increase LSW if this W level intersects topography in this cell

        if (arw(k,iw) > 0. .and. arw(k,iw) < .999 * arw0(iw)) then
           lsw(iw) = lsw(iw) + 1
        endif

     enddo  ! k

     nsw_max = max(nsw_max, lsw(iw))
  enddo

  volt(:,:) = 0.
  lve2(:)   = 0.

  ! Loop over all ATM grid columns

  !$omp parallel do private(iw,k,jv,iv)
  do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)

     ! special at model top
     volt(nza,iw) = dzt(nza) * arw0(iw)

     ! Loop over vertical levels

     do k = nza-1, lpw(iw), -1

        if (arw(k-1,iw) > 0.999 * arw(k,iw)) then
           volt(k,iw) = dzt(k) * arw(k,iw)
        else
           call volt_from_flux(k,iw)
        endif

     enddo

     ! Compute number of underground v[xyz]e2 levels in this IW column

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)
        lve2(iw) = max(lve2(iw), lpv(iv) - lpw(iw))
     enddo

  enddo
  !$omp end parallel do

  nve2_max = maxval(lve2(:))

  if (mdomain < 2) then

     ! Expand ARW and VOLT with height for spherical geometry

     do j = 1,jtab_w(jtw_grid)%jend; iw = jtab_w(jtw_grid)%iw(j)
        do k = lpw(iw), nza
           arw (k,iw) = arw (k,iw) * zfacm2(k)
           volt(k,iw) = volt(k,iw) * zfact(k)**2
        enddo
     enddo

     ! Expand ARV with height for spherical geometry

     do j = 1,jtab_v(jtv_grid)%jend; iv = jtab_v(jtv_grid)%iv(j)
        do k = lpv(iv), nza
           arv(k,iv) = arv(k,iv) * zfact(k)
        enddo
     enddo

  endif

  ! Lateral boundary copy of ARW, VOLT, LPW, and LSW

  do j = 1,jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp

     arw  (:,iw) = arw  (:,iwp)
     volt (:,iw) = volt (:,iwp)
     lpw    (iw) = lpw    (iwp)
     lsw    (iw) = lsw    (iwp)
  enddo

  ! Lateral boundary copy of ARV

  do j = 1,jtab_v(jtv_lbcp)%jend; iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp

     arv(:,iv) = arv(:,ivp)
     lpv(iv) = lpv(ivp)
  enddo

end subroutine ctrlvols_hex_nosfcg

!============================================================================

subroutine atm_grid_topocut(iw, npoly, npats, kw_pats, area_pats)

  use mem_grid,   only: nza, topm, topw, dnv, xem, yem, zem, xev, yev, zev
  use mem_ijtabs, only: itab_w

  implicit none

  integer, intent(in)    :: iw, npoly
  integer, intent(inout) :: npats(npoly), kw_pats(nza,npoly)
  real,    intent(inout) :: area_pats(nza,npoly)

  integer :: jv, jm1, jm2, iv, im1, im2
  integer :: t1, t2, t3, ipat
  integer :: jp1, jp2

  real :: z1,z2,z3,x1,x2,x3,y1,y2,y3
  real :: xm1, ym1, xm2, ym2
  real :: area
  real :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer :: tmpat(5,nza)  ! Flag denoting the edge or vertex where a new M
                           ! point is defined
  integer :: kmpat(5,nza)
  integer :: kwpat(nza),lpoly(nza),npat
  real :: dvm1, dvm2

  ! Loop over all polygon edges

  do jv = 1,npoly
     jm1 = jv - 1
     if (jm1 == 0) jm1 = npoly
     jm2 = jv

     iv  = itab_w(iw)%iv(jv)
     im1 = itab_w(iw)%im(jm1)
     im2 = itab_w(iw)%im(jm2)

     ! Each triangular sector of a hexagonal (or pentagonal or heptagonal) grid
     ! cell, formed by the IW point and two adjacent IM points, has topographic
     ! height specified on each of its vertices, so the triangle is treated as
     ! a planar topographic surface.  Find lines on this surface where it is
     ! intersected by each zm level of the atmospheric grid, and construct
     ! polygonal areas bounded by the lines of intersection and triangular
     ! sector edges (using the algorithm in cont3f).  Each of these polygons
     ! defines a separate contribution to the volume and surface area of a
     ! cut cell of the atmosphere grid.

     ! For cont3f algorithm, use local coordinate system with IW point at
     ! origin and IV edge being a line of constant positive y.

     dvm1 = sqrt((xev(iv) - xem(im1))**2 &
          +      (yev(iv) - yem(im1))**2 &
          +      (zev(iv) - zem(im1))**2)

     dvm2 = sqrt((xev(iv) - xem(im2))**2 &
          +      (yev(iv) - yem(im2))**2 &
          +      (zev(iv) - zem(im2))**2)

     x1 = 0.
     x2 = dvm1
     x3 = -dvm2
     y1 = 0.
     y2 = 0.5 * dnv(iv)
     y3 = 0.5 * dnv(iv)
     z1 = topw(iw)
     z2 = topm(im1)
     z3 = topm(im2)
     t1 = 1
     t2 = 2
     t3 = 3

     if     (z1 >= z2 .and. z2 >= z3) then
        call cont3sfc2(z1,z2,z3,x1,x2,x3,y1,y2,y3,t1,t2,t3, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z1 >= z3 .and. z3 >= z2) then
        call cont3sfc2(z1,z3,z2,x1,x3,x2,y1,y3,y2,t1,t3,t2, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     elseif (z2 >= z1 .and. z1 >= z3) then
        call cont3sfc2(z2,z1,z3,x2,x1,x3,y2,y1,y3,t2,t1,t3, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     elseif (z2 >= z3 .and. z3 >= z1) then
        call cont3sfc2(z2,z3,z1,x2,x3,x1,y2,y3,y1,t2,t3,t1, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z3 >= z1 .and. z1 >= z2) then
        call cont3sfc2(z3,z1,z2,x3,x1,x2,y3,y1,y2,t3,t1,t2, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z3 >= z2 .and. z2 >= z1) then
        call cont3sfc2(z3,z2,z1,x3,x2,x1,y3,y2,y1,t3,t2,t1, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     endif

     ! Loop over new patches that have been found in this sector

     do ipat = 1,npat

        ! Compute polygon area

        area = 0.

        do jp1 = 1,lpoly(ipat)
           jp2 = jp1 + 1
           if (jp2 > lpoly(ipat)) jp2 = 1

           xm1 = xmpat(jp1,ipat)
           ym1 = ympat(jp1,ipat)

           xm2 = xmpat(jp2,ipat)
           ym2 = ympat(jp2,ipat)

           area = area + 0.5 * (xm1 * ym2 - xm2 * ym1)
        enddo ! jp1

        npats         (jv) = npat
        kw_pats  (ipat,jv) = kwpat(ipat)
        area_pats(ipat,jv) = area

     enddo  ! ipat

  enddo   ! jv

end subroutine atm_grid_topocut

!===============================================================================

subroutine cont3sfc2(z1,z2,z3,x1,x2,x3,y1,y2,y3,t1,t2,t3,zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)

  use mem_grid, only: nza, zm

  implicit none

  real,    intent(in   ) :: z1,z2,z3,x1,x2,x3,y1,y2,y3
  integer, intent(in   ) :: t1,t2,t3
  real,    intent(inout) :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer, intent(inout) :: tmpat(5,nza), kmpat(5,nza)
  integer, intent(out  ) :: kwpat(nza),lpoly(nza),npat

  integer :: to1, to2, km, iflag
  real :: zo1, zo2, xo1, xo2, yo1, yo2, contlev ! "old" values of xmpat1,xmpat2,ympat1,ympat2

  ! Determine lowest zm level that exceeds z3, the lowest triangle vertex

  km = 1
  do while (z3 >= zm(km) .and. km < nza-1)
     km = km + 1
  enddo
  contlev = zm(km)
  npat = 1
  kmpat(:,:) = 0

  ! If all 3 points are in the same contour interval

  if (contlev >= z1 .or. contlev <= z3) then

     zmpat(1,npat) = z1
     xmpat(1,npat) = x1
     ympat(1,npat) = y1
     tmpat(1,npat) = t1

     zmpat(2,npat) = z2
     xmpat(2,npat) = x2
     ympat(2,npat) = y2
     tmpat(2,npat) = t2

     zmpat(3,npat) = z3
     xmpat(3,npat) = x3
     ympat(3,npat) = y3
     tmpat(3,npat) = t3

     kwpat(npat) = km

     lpoly(npat) = 3

     return
  endif

  ! Draw all contour lines of value less than z1

  iflag = 0

  do while (contlev < z1)

     if (iflag > 0) then
        zo1 = zmpat(1,npat-1)
        xo1 = xmpat(1,npat-1)
        yo1 = ympat(1,npat-1)
        to1 = tmpat(1,npat-1)

        zo2 = zmpat(2,npat-1)
        xo2 = xmpat(2,npat-1)
        yo2 = ympat(2,npat-1)
        to2 = tmpat(2,npat-1)
     endif

     if ( abs(z1-z3) > 1.e-25 ) then
        zmpat(1,npat) = contlev
        xmpat(1,npat) = x3 + (x1-x3) * (contlev-z3) / (z1-z3)
        ympat(1,npat) = y3 + (y1-y3) * (contlev-z3) / (z1-z3)
        tmpat(1,npat) = t2*10
        kmpat(1,npat) = km
     else
        zmpat(1,npat) = z3
        xmpat(1,npat) = x3
        ympat(1,npat) = y3
        tmpat(1,npat) = t3
     endif

     if (contlev < z2) then

        if ( abs(z2-z3) > 1.e-25 ) then
           zmpat(2,npat) = contlev
           xmpat(2,npat) = x3 + (x2-x3) * (contlev-z3) / (z2-z3)
           ympat(2,npat) = y3 + (y2-y3) * (contlev-z3) / (z2-z3)
           tmpat(2,npat) = t1*10
           kmpat(2,npat) = km
        else
           zmpat(2,npat) = z3
           xmpat(2,npat) = x3
           ympat(2,npat) = y3
           tmpat(2,npat) = t3
        endif

        if (iflag == 0) then  ! lowest contour interval: z3 is a node
           zmpat(3,npat) = z3
           xmpat(3,npat) = x3
           ympat(3,npat) = y3
           tmpat(3,npat) = t3

           kwpat(npat) = km
           lpoly(npat) = 3
        else
           zmpat(3,npat) = zo2
           xmpat(3,npat) = xo2
           ympat(3,npat) = yo2
           tmpat(3,npat) = to2
           kmpat(3,npat) = km - 1

           zmpat(4,npat) = zo1
           xmpat(4,npat) = xo1
           ympat(4,npat) = yo1
           tmpat(4,npat) = to1
           kmpat(4,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 4
        endif

        if (km == nza-1) then  ! highest model level: z1 and z2 are nodes
                               ! (should not actually happen)
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = z1
           xmpat(4,npat) = x1
           ympat(4,npat) = y1
           tmpat(4,npat) = t1

           kwpat(npat) = km
           lpoly(npat) = 4

           return
        elseif (zm(km+1) >= z1) then  ! highest contour interval:
                                      ! z1 and z2 are nodes
           zmpat(1,npat+1) = zmpat(2,npat)
           xmpat(1,npat+1) = xmpat(2,npat)
           ympat(1,npat+1) = ympat(2,npat)
           tmpat(1,npat+1) = tmpat(2,npat)
           kmpat(1,npat+1) = km

           zmpat(2,npat+1) = zmpat(1,npat)
           xmpat(2,npat+1) = xmpat(1,npat)
           ympat(2,npat+1) = ympat(1,npat)
           tmpat(2,npat+1) = tmpat(1,npat)
           kmpat(2,npat+1) = km

           zmpat(3,npat+1) = z1
           xmpat(3,npat+1) = x1
           ympat(3,npat+1) = y1
           tmpat(3,npat+1) = t1

           zmpat(4,npat+1) = z2
           xmpat(4,npat+1) = x2
           ympat(4,npat+1) = y2
           tmpat(4,npat+1) = t2

           kwpat(npat+1) = km + 1
           lpoly(npat+1) = 4

           npat = npat + 1

           return
        endif

        iflag = 1

     else

        if ( abs(z1-z2) > 1.e-25 ) then
           zmpat(2,npat) = contlev
           xmpat(2,npat) = x2 + (x1-x2) * (contlev-z2) / (z1-z2)
           ympat(2,npat) = y2 + (y1-y2) * (contlev-z2) / (z1-z2)
           tmpat(2,npat) = t3*10
           kmpat(2,npat) = km
        else
           zmpat(2,npat) = z2
           xmpat(2,npat) = x2
           ympat(2,npat) = y2
           tmpat(2,npat) = t2
        endif

        if (iflag == 0) then  ! lowest contour color: z2 and z3 are nodes
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = z3
           xmpat(4,npat) = x3
           ympat(4,npat) = y3
           tmpat(4,npat) = t3

           kwpat(npat) = km
           lpoly(npat) = 4
        elseif (iflag == 1) then ! switching from contlev < z2 to contlev > z2
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = zo2
           xmpat(4,npat) = xo2
           ympat(4,npat) = yo2
           tmpat(4,npat) = to2
           kmpat(4,npat) = km - 1

           zmpat(5,npat) = zo1
           xmpat(5,npat) = xo1
           ympat(5,npat) = yo1
           tmpat(5,npat) = to1
           kmpat(5,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 5
        else                     ! keeping to contlev > z2
           zmpat(3,npat) = zo2
           xmpat(3,npat) = xo2
           ympat(3,npat) = yo2
           tmpat(3,npat) = to2
           kmpat(3,npat) = km - 1

           zmpat(4,npat) = zo1
           xmpat(4,npat) = xo1
           ympat(4,npat) = yo1
           tmpat(4,npat) = to1
           kmpat(4,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 4
        endif

        if (km == nza-1) then  ! highest model level: z1 is a node
           zmpat(3,npat) = z1  ! (should not actually happen)
           xmpat(3,npat) = x1  ! (should not actually happen)
           ympat(3,npat) = y1
           tmpat(3,npat) = t1

           kwpat(npat) = km
           lpoly(npat) = 3

           return
        elseif (zm(km+1) >= z1) then  ! highest contour interval: z1 is a node
           zmpat(1,npat+1) = zmpat(2,npat)
           xmpat(1,npat+1) = xmpat(2,npat)
           ympat(1,npat+1) = ympat(2,npat)
           tmpat(1,npat+1) = tmpat(2,npat)
           kmpat(1,npat+1) = km

           zmpat(2,npat+1) = zmpat(1,npat)
           xmpat(2,npat+1) = xmpat(1,npat)
           ympat(2,npat+1) = ympat(1,npat)
           tmpat(2,npat+1) = tmpat(1,npat)
           kmpat(2,npat+1) = km

           zmpat(3,npat+1) = z1
           xmpat(3,npat+1) = x1
           ympat(3,npat+1) = y1
           tmpat(3,npat+1) = t1

           kwpat(npat+1) = km + 1
           lpoly(npat+1) = 3

           npat = npat + 1

           return
        endif

        iflag = 2

     endif

     npat = npat + 1
     km = km + 1
     contlev = zm(km)

  enddo

end subroutine cont3sfc2

!==========================================================================

subroutine reverse_polygon(zmpat, xmpat, ympat, tmpat, kmpat, lpoly, npat)

  use mem_grid, only: nza

  implicit none

  real,    intent(inout) :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer, intent(inout) :: tmpat(5,nza), kmpat(5,nza)
  integer, intent(in   ) :: npat, lpoly(nza)

  integer :: ipat, j, jc, istore
  real :: store

  do ipat = 1,npat
     do j = 1,lpoly(ipat) / 2
        jc = lpoly(ipat) + 1 - j

        store = zmpat(j,ipat)
        zmpat(j,ipat) = zmpat(jc,ipat)
        zmpat(jc,ipat) = store

        store = xmpat(j,ipat)
        xmpat(j,ipat) = xmpat(jc,ipat)
        xmpat(jc,ipat) = store

        store = ympat(j,ipat)
        ympat(j,ipat) = ympat(jc,ipat)
        ympat(jc,ipat) = store

        istore = tmpat(j,ipat)
        tmpat(j,ipat) = tmpat(jc,ipat)
        tmpat(jc,ipat) = istore

        istore = kmpat(j,ipat)
        kmpat(j,ipat) = kmpat(jc,ipat)
        kmpat(jc,ipat) = istore
     enddo
  enddo

end subroutine reverse_polygon

!==========================================================================

subroutine volt_from_flux(k,iw)

  use mem_grid,    only: arv, arw, arw0, volt, dnu, dzt, vnx, vny, vnz, &
                         xew, yew, zew
  use mem_ijtabs,  only: itab_w
  use consts_coms, only: pi2, eradi
  use misc_coms,   only: mdomain

  implicit none

  integer, intent(in) :: k, iw

  integer, parameter :: naz = 18
  integer :: iaz, jv, iv, iloop

  real :: az, uzonal, umerid, raxis, vxe, vye, vze, vc, wc, dirv, hflux
  real :: hflux_cut_in, hflux_cut_out, hflux_opn_in, hflux_opn_out
  real :: vflux_cut_in, vflux_cut_out, vflux_opn_in
  real :: del_arw, factor, maxfactor

  real, parameter :: onem = 1.0 - 1.e-6

  maxfactor = 0.
  iloop = 1

  ! Loop over increments of azimuth angle, counterclockwise from east, up to half a circle

  do iaz = 1, naz
     az = pi2 * real(iaz - 1) / real(naz)

     ! Horizontal velocity components assuming 1 m/s total horizontal wind

     uzonal = cos(az)
     umerid = sin(az)

     raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

     if (mdomain < 2 .and. raxis > 1.e3) then
        vxe = -yew(iw) / raxis * uzonal - xew(iw) * zew(iw) * eradi / raxis * umerid
        vye =  xew(iw) / raxis * uzonal - yew(iw) * zew(iw) * eradi / raxis * umerid
        vze =                                                 raxis * eradi * umerid
     else
        vxe = uzonal
        vye = umerid
        vze = 0.
     endif

     ! Project cell-center velocity components onto each hex face, and get
     ! horizontal fluxes for this cut cell and for open cell of same dimensions

     hflux_cut_in  = 0.
     hflux_cut_out = 0.
     hflux_opn_in  = 0.
     hflux_opn_out = 0.

     do jv = 1, itab_w(iw)%npoly
        iv   = itab_w(iw)%iv(jv)
        dirv = itab_w(iw)%dirv(jv)

        vc = vnx(iv) * vxe &
           + vny(iv) * vye &
           + vnz(iv) * vze

        hflux = dirv * vc

        if (hflux > 0.) then
           hflux_cut_in = hflux_cut_in + hflux * arv(k,iv)
           hflux_opn_in = hflux_opn_in + hflux * dnu(iv) * dzt(k)
        else
           hflux_cut_out = hflux_cut_out - hflux * arv(k,iv)
           hflux_opn_out = hflux_opn_out - hflux * dnu(iv) * dzt(k)
        endif
     enddo

     ! Compute vertical velocity that gives zero net flux for cut cell.
     ! For arbitrarily steep terrain slope, limit vertical velocity to
     ! +/- 10 m/s, which would correspond to a slope of 10:1.

     del_arw = max(.01 * arw(k,iw), arw(k,iw) - arw(k-1,iw))
     wc = max(-10., min(10., (hflux_cut_in - hflux_cut_out) / del_arw))

     ! Get vertical fluxes for this cut cell and for open cell of same dimensions.

     if (wc > 0.) then
        vflux_cut_in  = wc * arw(k-1,iw)
        vflux_cut_out = wc * arw(k,iw)
        vflux_opn_in  = wc * arw0(iw)
     else
        vflux_cut_in  = -wc * arw(k,iw)
        vflux_cut_out = -wc * arw(k-1,iw)
        vflux_opn_in  = -wc * arw0(iw)
     endif

     ! Compute volt

     factor = (hflux_cut_in + vflux_cut_in) / (hflux_opn_in + vflux_opn_in)

     if (factor > maxfactor) maxfactor = factor

     ! Verify that cut cell net flux is zero

!     if (iaz == 1) then
!        print*, ' '
!        write(*,'(120a)') '          iaz   iw     k   hfcutin    hfcutout    vfcutin    vfcutout  fsum   factor  maxfactor'
!        print*, ' '
!     endif

!     write(*,'(a,3i6,4f11.1,f8.1,2f8.3)') 'volts ', iaz, iw, k, &
!        hflux_cut_in, hflux_cut_out, vflux_cut_in, vflux_cut_out, &
!        hflux_cut_in - hflux_cut_out + vflux_cut_in - vflux_cut_out, &
!        factor, maxfactor

     if (maxfactor > onem) then
        maxfactor = 1.0
        exit
     endif

  enddo

  volt(k,iw) = arw0(iw) * dzt(k) * maxfactor

end subroutine volt_from_flux

!==========================================================================

subroutine isort2(n,i1,i2)

  ! Sort n integers in each of i1 and i2 into ascending order by i1

  implicit none

  integer, intent(in)    :: n
  integer, intent(inout) :: i1(n), i2(n)
  integer                :: i, j, i0

  do i = 1, n-1
     do j = i+1, n

        if (i1(j) < i1(i)) then
           i0 = i1(i)
           i1(i) = i1(j)
           i1(j) = i0

           i0 = i2(i)
           i2(i) = i2(j)
           i2(j) = i0
        endif

     enddo
  enddo

end subroutine isort2

