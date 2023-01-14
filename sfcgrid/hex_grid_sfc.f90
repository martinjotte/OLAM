subroutine voronoi_sfc()

  use mem_sfcg,     only: nmsfc, nvsfc, nwsfc, itab_msfc, itab_vsfc, &
                          itab_wsfc, sfcg, alloc_sfcgrid1

  use mem_delaunay, only: itab_md, itab_ud, itab_wd, nmd, nud, nwd, &
                          xemd, yemd, zemd

  use misc_coms,    only: mdomain
  use consts_coms,  only: erad, eradi, piu180

  implicit none

  integer :: iw1, iw2, iw3, im, iv, iw
  integer :: imd,iud,iwd,j
  real    :: expansion, raxis, raxisi
  real    :: xebc,yebc,zebc
  real    :: x1,x2,x3,y1,y2,y3
  real    :: dx12,dx13,dx23
  real    :: s1,s2,s3
  real    :: xcc,ycc
  real    :: dxe,dye,dze
  real    :: sinwslat,coswslat
  real    :: sinwslon,coswslon

  ! Interchange grid dimensions

  nmsfc = nwd
  nvsfc = nud
  nwsfc = nmd

  ! Allocate Voronoi set of arrays

  call alloc_sfcgrid1(nmsfc, nvsfc, nwsfc, alloc_xyzew=.false.)

  ! Transfer information from Delaunay to Voronoi arrays

  call move_alloc(xemd, sfcg%xew)
  call move_alloc(yemd, sfcg%yew)
  call move_alloc(zemd, sfcg%zew)

  ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
  ! This establishes W cell as voronoi.

  !$omp parallel
  !$omp do private(iwd,iw1,iw2,iw3,xebc,yebc,zebc,expansion,raxis,raxisi, &
  !$omp            sinwslat,coswslat,sinwslon,coswslon,dxe,dye,dze,x1,y1, &
  !$omp            x2,y2,x3,y3,dx12,dx13,dx23,s1,s2,s3,ycc,xcc)
  do im = 2,nmsfc
     iwd = im

     if (any(itab_wd(iwd)%im(1:3) < 2)) cycle

     ! Indices of 3 M points surrounding WD point

     iw1 = itab_wd(iwd)%im(1)
     iw2 = itab_wd(iwd)%im(2)
     iw3 = itab_wd(iwd)%im(3)

     ! First, compute barycenter of 3 W points

     xebc = (sfcg%xew(iw1) + sfcg%xew(iw2) + sfcg%xew(iw3)) / 3.
     yebc = (sfcg%yew(iw1) + sfcg%yew(iw2) + sfcg%yew(iw3)) / 3.
     zebc = (sfcg%zew(iw1) + sfcg%zew(iw2) + sfcg%zew(iw3)) / 3.

     if (mdomain <= 1) then

        ! If mdomain <= 1, push M point coordinates out to earth radius

        expansion = erad / sqrt( xebc ** 2 &
                               + yebc ** 2 &
                               + zebc ** 2 )

        xebc = xebc * expansion
        yebc = yebc * expansion
        zebc = zebc * expansion

        ! Get latitude and longitude of barycentric point

        raxis  = sqrt(xebc ** 2 + yebc ** 2)

        sinwslat = zebc  * eradi
        coswslat = raxis * eradi

        ! For points less than 100 m from Earth's polar axis, make arbitrary
        ! assumption that longitude = 0 deg.  This is just to settle on a PS
        ! planar coordinate system in which to do the algebra.

        if (raxis >= 1.e2) then
           raxisi = 1.0 / raxis
           sinwslon = yebc * raxisi
           coswslon = xebc * raxisi
        else
           sinwslon = 0.
           coswslon = 1.
        endif

        ! Transform 3 W points to PS coordinates

        dxe = sfcg%xew(iw1) - xebc
        dye = sfcg%yew(iw1) - yebc
        dze = sfcg%zew(iw1) - zebc
        call de_ps(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,x1,y1)

        dxe = sfcg%xew(iw2) - xebc
        dye = sfcg%yew(iw2) - yebc
        dze = sfcg%zew(iw2) - zebc
        call de_ps(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,x2,y2)

        dxe = sfcg%xew(iw3) - xebc
        dye = sfcg%yew(iw3) - yebc
        dze = sfcg%zew(iw3) - zebc
        call de_ps(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,x3,y3)

     else

        ! For Cartesian domain, use given planar X,Y coordinates

        x1 = sfcg%xew(iw1) - xebc
        x2 = sfcg%xew(iw2) - xebc
        x3 = sfcg%xew(iw3) - xebc

        y1 = sfcg%yew(iw1) - yebc
        y2 = sfcg%yew(iw2) - yebc
        y3 = sfcg%yew(iw3) - yebc

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

     ! For global domain, transform circumcenter from PS to earth coordinates

     if (mdomain <= 1) then

        call ps_de(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xcc,ycc)

        sfcg%xem(im) = dxe + xebc
        sfcg%yem(im) = dye + yebc
        sfcg%zem(im) = dze + zebc

     else

        sfcg%xem(im) = xcc + xebc
        sfcg%yem(im) = ycc + yebc

     endif

  enddo
  !$omp end do nowait

  if (mdomain <= 1) then
     !$omp do private(expansion)
     do im = 2, nmsfc

        expansion = erad / sqrt( sfcg%xem(im) ** 2 &
                               + sfcg%yem(im) ** 2 &
                               + sfcg%zem(im) ** 2 )

        sfcg%xem(im) = sfcg%xem(im) * expansion
        sfcg%yem(im) = sfcg%yem(im) * expansion
        sfcg%zem(im) = sfcg%zem(im) * expansion

     enddo
     !$omp end do nowait
  endif

  ! Loop over V points

  !$omp do private(iud)
  do iv = 2,nvsfc
     iud = iv

     itab_vsfc(iv)%imn(1:2)  = itab_ud(iud)%iw(1:2)
     itab_vsfc(iv)%iwn(1:2)  = itab_ud(iud)%im(1:2)
  enddo
  !$omp end do

  ! Loop over WSFC points

  !$omp do private(imd,j,im,iwd,iv,iw1,iw2)
  do iw = 2,nwsfc
     imd = iw

     itab_wsfc(iw)%npoly = itab_md(imd)%npoly

     ! Loop over IM/IV neighbors of IW

     do j = 1,itab_wsfc(iw)%npoly
        im = itab_md(imd)%iw(j)
        iwd = im
        iv = itab_md(imd)%iu(j)

        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)

        itab_wsfc(iw)%imn(j) = im
        itab_wsfc(iw)%ivn(j) = iv

        if (iw1 == iw) then
           itab_wsfc(iw)%iwn(j)  = iw2
           itab_wsfc(iw)%dirv(j) = -1.
        else
           itab_wsfc(iw)%iwn(j)  = iw1
           itab_wsfc(iw)%dirv(j) = 1.
        endif
     enddo

  enddo
  !$omp end do

  ! Loop over MSFC points

  !$omp do private(iwd)
  do im = 2,nmsfc
     iwd = im
     itab_msfc(im)%ivn(1:3) = itab_wd(iwd)%iu(1:3)
     itab_msfc(im)%iwn(1:3) = itab_wd(iwd)%im(1:3)
     itab_msfc(im)%imn(1:3) = itab_wd(iwd)%iw(1:3)
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(itab_md,itab_ud,itab_wd)

end subroutine voronoi_sfc

!===============================================================================

subroutine grid_geometry_hex_sfc()

  use mem_sfcg,    only: nmsfc, nvsfc, nwsfc, itab_msfc, itab_vsfc, itab_wsfc, sfcg
  use misc_coms,   only: mdomain, nxp
  use consts_coms, only: erad, piu180, r8
  use oplot_coms,  only: op

  implicit none

  integer       :: im,iv,iw,ivn
  integer       :: im1,im2
  integer       :: iw1,iw2
  integer       :: j,npoly,j1,j2
  real          :: raxis, expansion
  real          :: dvm1,dvm2
  real          :: xm1,xm2,xv,ym1,ym2,yv,frac,alpha
  real          :: xw1,xw2,yw1,yw2
  real          :: xp1, yp1, xp2, yp2
  real          :: xq1, yq1, xq2, yq2, psiz, vsprd
  integer       :: iskipm, iskipw
  real          :: quarter_kite(2,nvsfc)
  character(10) :: string

  integer               :: lwork, info
  real(r8)              :: a(3,7), b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7)
  real(r8)              :: vrot_x(7), vrot_y(7)
  real(r8), allocatable :: work(:)
  real(r8)              :: wsize(1), vdotw, vmagi, fact

  ! Loop over all M points and compute their latitude and longitude

  !$omp parallel
  !$omp do private(raxis)
  do im = 2,nmsfc
     if (mdomain <= 1) then
        raxis = sqrt(sfcg%xem(im) ** 2 + sfcg%yem(im) ** 2)  ! dist from earth axis
        sfcg%glatm(im) = atan2(sfcg%zem(im),raxis)        * piu180
        sfcg%glonm(im) = atan2(sfcg%yem(im),sfcg%xem(im)) * piu180
     else
        sfcg%glatm(im) = 0.  ! want it this way?
        sfcg%glonm(im) = 0.  ! want it this way?
     endif

     ! Fill global index (replaced later if this run is parallel)

     itab_msfc(im)%imglobe = im

  enddo
  !$omp end do

  ! Loop over all V points

  !$omp do private(im1,im2,iw1,iw2,expansion,dvm1,dvm2,frac)
  do iv = 2,nvsfc

     ! Fill global index (replaced later if this run is parallel)

     itab_vsfc(iv)%ivglobe = iv

     ! M-point indices of two end points of V segment

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     ! W-point indices on either side of V segment

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     ! V point is midway between W points of Voronoi cells

     sfcg%xev(iv) = .5 * (sfcg%xew(iw1) + sfcg%xew(iw2))
     sfcg%yev(iv) = .5 * (sfcg%yew(iw1) + sfcg%yew(iw2))
     sfcg%zev(iv) = .5 * (sfcg%zew(iw1) + sfcg%zew(iw2))

     ! If mdomain <= 1, push V point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt(sfcg%xev(iv) ** 2 &
                              + sfcg%yev(iv) ** 2 &
                              + sfcg%zev(iv) ** 2 )

        sfcg%xev(iv) = sfcg%xev(iv) * expansion
        sfcg%yev(iv) = sfcg%yev(iv) * expansion
        sfcg%zev(iv) = sfcg%zev(iv) * expansion

     endif

     ! Normal distance across U face

     sfcg%dnu(iv) = sqrt( (sfcg%xem(im1) - sfcg%xem(im2))**2 &
                        + (sfcg%yem(im1) - sfcg%yem(im2))**2 &
                        + (sfcg%zem(im1) - sfcg%zem(im2))**2 )

     sfcg%dniu(iv) = 1. / sfcg%dnu(iv)

     sfcg%unx(iv) = (sfcg%xem(im2) - sfcg%xem(im1)) / sfcg%dnu(iv)
     sfcg%uny(iv) = (sfcg%yem(im2) - sfcg%yem(im1)) / sfcg%dnu(iv)
     sfcg%unz(iv) = (sfcg%zem(im2) - sfcg%zem(im1)) / sfcg%dnu(iv)

     ! Normal distance across V face

     sfcg%dnv(iv) = sqrt( (sfcg%xew(iw1) - sfcg%xew(iw2))**2 &
                        + (sfcg%yew(iw1) - sfcg%yew(iw2))**2 &
                        + (sfcg%zew(iw1) - sfcg%zew(iw2))**2 )

     sfcg%dniv(iv) = 1. / sfcg%dnv(iv)

     sfcg%vnx(iv) = (sfcg%xew(iw2) - sfcg%xew(iw1)) / sfcg%dnv(iv)
     sfcg%vny(iv) = (sfcg%yew(iw2) - sfcg%yew(iw1)) / sfcg%dnv(iv)
     sfcg%vnz(iv) = (sfcg%zew(iw2) - sfcg%zew(iw1)) / sfcg%dnv(iv)

     ! Compute IM1 and IM2 values of quarter kite area,
     ! and add to ARM0 and ARW0 arrays

     dvm1 = sqrt((sfcg%xev(iv) - sfcg%xem(im1))**2 &
          +      (sfcg%yev(iv) - sfcg%yem(im1))**2 &
          +      (sfcg%zev(iv) - sfcg%zem(im1))**2)

     dvm2 = sqrt((sfcg%xev(iv) - sfcg%xem(im2))**2 &
          +      (sfcg%yev(iv) - sfcg%yem(im2))**2 &
          +      (sfcg%zev(iv) - sfcg%zem(im2))**2)

     ! Fractional distance along V edge where intersection with U edge is located

     frac = dvm1 * sfcg%dniu(iv)

     if (im1 > 1 .and. im2 > 1 .and. (frac < .0001 .or. frac > .9999)) then
        print*, 'Non-intersecting U-V edges detected in sfcgrid geometry'
        print*, 'FRAC  = ',frac
        print*, 'IW1 = ', iw1, ' IW2 = ', iw2
        print*, 'IV    = ',iv
        print*, 'GLATM = ', sfcg%glatm(im1)
        print*, 'GLONM = ', sfcg%glonm(im1)
        print*, 'SFCG%XEV = ',sfcg%xev(iv)
        print*, 'SFCG%YEV = ',sfcg%yev(iv)
        print*, 'SFCG%ZEV = ',sfcg%zev(iv)

        print*, 'sfcg%dnu(iv),sfcg%dniu(iv) ',sfcg%dnu(iv),sfcg%dniu(iv)

        stop 'STOP U-V edges'
     endif

     quarter_kite(1,iv) = .25 * dvm1 * sfcg%dnv(iv)
     quarter_kite(2,iv) = .25 * dvm2 * sfcg%dnv(iv)

  enddo
  !$omp end do
  !$omp end parallel

  !dir$ novector
  do iv = 2, nvsfc
     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     sfcg%arm0(im1) = sfcg%arm0(im1) + 2. * quarter_kite(1,iv)
     sfcg%arm0(im2) = sfcg%arm0(im2) + 2. * quarter_kite(2,iv)

     sfcg%area(iw1) = sfcg%area(iw1) + quarter_kite(1,iv) + quarter_kite(2,iv)
     sfcg%area(iw2) = sfcg%area(iw2) + quarter_kite(1,iv) + quarter_kite(2,iv)
  enddo

  !$omp parallel private(a,b,wsize,lwork,work,info)

  !$omp do private(raxis,npoly,j1,j2,ivn,iw1,iw2,xw1,yw1,xw2,yw2,xv,yv,alpha)
  do iw = 2,nwsfc

     ! Fill global index (replaced later if this run is parallel)

     itab_wsfc(iw)%iwglobe = iw

     ! Fill outward unit vector components and latitude and longitude of W point

     if (mdomain <= 1) then

        raxis = sqrt(sfcg%xew(iw) ** 2 + sfcg%yew(iw) ** 2)

        sfcg%glatw(iw) = atan2(sfcg%zew(iw),raxis)        * piu180
        sfcg%glonw(iw) = atan2(sfcg%yew(iw),sfcg%xew(iw)) * piu180

        sfcg%wnx(iw) = sfcg%xew(iw) / erad
        sfcg%wny(iw) = sfcg%yew(iw) / erad
        sfcg%wnz(iw) = sfcg%zew(iw) / erad

     else

        sfcg%glatw(iw) = 0. ! want it this way?
        sfcg%glonw(iw) = 0. ! want it this way?

        sfcg%wnx(iw) = 0.0
        sfcg%wny(iw) = 0.0
        sfcg%wnz(iw) = 1.0

     endif

     ! Number of polygon edges/vertices

     npoly = itab_wsfc(iw)%npoly

     ! Loop over all polygon edges

     do j2 = 1,npoly
        j1 = j2 - 1
        if (j2 == 1) j1 = npoly

        ivn = itab_wsfc(iw)%ivn(j2)
        iw2 = itab_wsfc(iw)%iwn(j2)
        iw1 = itab_wsfc(iw)%iwn(j1)

        ! Fractional area of sfcg%area(iw) that is occupied by M and V sectors

        if (itab_vsfc(ivn)%iwn(1) == iw1) then
           itab_wsfc(iw)%farm(j1) = itab_wsfc(iw)%farm(j1) + quarter_kite(1,ivn) / sfcg%area(iw)
           itab_wsfc(iw)%farm(j2) = itab_wsfc(iw)%farm(j2) + quarter_kite(2,ivn) / sfcg%area(iw)
        else
           itab_wsfc(iw)%farm(j1) = itab_wsfc(iw)%farm(j1) + quarter_kite(2,ivn) / sfcg%area(iw)
           itab_wsfc(iw)%farm(j2) = itab_wsfc(iw)%farm(j2) + quarter_kite(1,ivn) / sfcg%area(iw)
        endif
        itab_wsfc(iw)%farv(j2) = (quarter_kite(1,ivn) + quarter_kite(2,ivn)) / sfcg%area(iw)

        ! Evaluate x,y coordinates of IW1 and IW2 points on polar stereographic plane
        ! tangent at IW

        if (mdomain <= 1) then
           call e_ps(sfcg%xew(iw1),sfcg%yew(iw1),sfcg%zew(iw1),sfcg%glatw(iw),sfcg%glonw(iw),xw1,yw1)
           call e_ps(sfcg%xew(iw2),sfcg%yew(iw2),sfcg%zew(iw2),sfcg%glatw(iw),sfcg%glonw(iw),xw2,yw2)
           call e_ps(sfcg%xev(ivn), sfcg%yev(ivn), sfcg%zev(ivn), sfcg%glatw(iw),sfcg%glonw(iw),xv,yv)
        else
           xw1 = sfcg%xew(iw1) - sfcg%xew(iw)
           yw1 = sfcg%yew(iw1) - sfcg%yew(iw)
           xw2 = sfcg%xew(iw2) - sfcg%xew(iw)
           yw2 = sfcg%yew(iw2) - sfcg%yew(iw)
           xv  = sfcg%xev(ivn) - sfcg%xew(iw)
           yv  = sfcg%yev(ivn) - sfcg%yew(iw)
        endif

        ! Coefficients for eastward and northward components of gradient (they apply at M points)

        itab_wsfc(iw)%gxps1(j1) =  yw2 / (xw1 * yw2 - xw2 * yw1)
        itab_wsfc(iw)%gxps2(j1) = -yw1 / (xw1 * yw2 - xw2 * yw1)

        itab_wsfc(iw)%gyps1(j1) = -xw2 / (xw1 * yw2 - xw2 * yw1)
        itab_wsfc(iw)%gyps2(j1) =  xw1 / (xw1 * yw2 - xw2 * yw1)

        if (itab_wsfc(iw)%dirv(j2) < 0.) then
           alpha = atan2(yw2,xw2)   ! VC(ivn) direction counterclockwise from east

           itab_vsfc(ivn)%cosv(1) = cos(alpha)
           itab_vsfc(ivn)%sinv(1) = sin(alpha)

           itab_vsfc(ivn)%dxps(1) = xv
           itab_vsfc(ivn)%dyps(1) = yv
        else
           alpha = atan2(-yw2,-xw2) ! VC(ivn) direction counterclockwise from east

           itab_vsfc(ivn)%cosv(2) = cos(alpha)
           itab_vsfc(ivn)%sinv(2) = sin(alpha)

           itab_vsfc(ivn)%dxps(2) = xv
           itab_vsfc(ivn)%dyps(2) = yv
        endif

     enddo

  enddo
  !$omp end do

  ! Scale eastward and northward gradient components by farm

  !$omp do private(j)
  do iw = 2, nwsfc
     do j = 1, itab_wsfc(iw)%npoly
        itab_wsfc(iw)%gxps1(j) = itab_wsfc(iw)%gxps1(j) * itab_wsfc(iw)%farm(j)
        itab_wsfc(iw)%gyps1(j) = itab_wsfc(iw)%gyps1(j) * itab_wsfc(iw)%farm(j)

        itab_wsfc(iw)%gxps2(j) = itab_wsfc(iw)%gxps2(j) * itab_wsfc(iw)%farm(j)
        itab_wsfc(iw)%gyps2(j) = itab_wsfc(iw)%gyps2(j) * itab_wsfc(iw)%farm(j)
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

     !$omp do private(npoly, fo, j, ivn, vdotw, vmagi, fact, &
     !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y)
     do iw = 2, nwsfc

        npoly = itab_wsfc(iw)%npoly

        ! Default coefficients from Perot
        fo(1:npoly) = 2.0_r8 * itab_wsfc(iw)%farv(1:npoly)

        if (mdomain < 2) then

           do j = 1, npoly
              ivn = itab_wsfc(iw)%ivn(j)

              ! Compute the components of the V unit normals perpendicular to W

              vdotw = sfcg%vnx(ivn)*sfcg%wnx(iw) + sfcg%vny(ivn)*sfcg%wny(iw) + sfcg%vnz(ivn)*sfcg%wnz(iw)

              vnx_ps(j) = sfcg%vnx(ivn) - vdotw * sfcg%wnx(iw)
              vny_ps(j) = sfcg%vny(ivn) - vdotw * sfcg%wny(iw)
              vnz_ps(j) = sfcg%vnz(ivn) - vdotw * sfcg%wnz(iw)

              ! Normalize these new vectors to unit length

              vmagi = 1.0_r8 / sqrt( vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2 )

              vnx_ps(j) = vnx_ps(j) * vmagi
              vny_ps(j) = vny_ps(j) * vmagi
              vnz_ps(j) = vnz_ps(j) * vmagi

              ! Rotate these new unit normals to a coordinate system with Z aligned with W

              if (sfcg%wnz(iw) >= 0.0) then

                 fact = ( sfcg%wny(iw)*vnx_ps(j) - sfcg%wnx(iw)*vny_ps(j) ) / ( 1._r8 + sfcg%wnz(iw) )

                 vrot_x(j) = vnx_ps(j)*sfcg%wnz(iw) - vnz_ps(j)*sfcg%wnx(iw) + sfcg%wny(iw)*fact
                 vrot_y(j) = vny_ps(j)*sfcg%wnz(iw) - vnz_ps(j)*sfcg%wny(iw) - sfcg%wnx(iw)*fact

              else

                 fact = ( sfcg%wny(iw)*vnx_ps(j) - sfcg%wnx(iw)*vny_ps(j) ) / ( 1._r8 - sfcg%wnz(iw) )

                 vrot_x(j) = -vnx_ps(j)*sfcg%wnz(iw) + vnz_ps(j)*sfcg%wnx(iw) + sfcg%wny(iw)*fact
                 vrot_y(j) = -vny_ps(j)*sfcg%wnz(iw) + vnz_ps(j)*sfcg%wny(iw) - sfcg%wnx(iw)*fact

              endif

           enddo

        else

           do j = 1, npoly
              ivn = itab_wsfc(iw)%ivn(j)

              vnx_ps(j) = sfcg%vnx(ivn)
              vny_ps(j) = sfcg%vny(ivn)
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

        call dgels( 'N', 3, npoly, 1, a(:,1:npoly), 3, b, 7, work, lwork, info )

        ! Vector b is now the correction to the coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly)

        if (info == 0 .and. all(b(1:npoly) > 0.03_r8) .and. all(b(1:npoly) < 0.7_r8)) then

           itab_wsfc(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly)
           itab_wsfc(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly)
           itab_wsfc(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly)

        else

           write(6,'(a,i9)')   "Problem optimizing vector coefficients for iw = ", iw
           write(6,'(7f10.4)') sfcg%glatw(iw), sfcg%glonw(iw)
           write(6,'(i5)')     info
           write(6,'(7f10.4)') real(b (1:npoly))
           write(6,'(7f10.4)') real(fo(1:npoly))
           write(6,'(a)')      "Using default coefficients."

           itab_wsfc(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
           itab_wsfc(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
           itab_wsfc(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

        endif

     enddo
     !omp end do

     deallocate(work)

  else

     !$omp do private(npoly)
     do iw = 2, nwsfc
        npoly = itab_wsfc(iw)%npoly
        itab_wsfc(iw)%ecvec_vx(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vnx(itab_wsfc(iw)%ivn(1:npoly))
        itab_wsfc(iw)%ecvec_vy(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vny(itab_wsfc(iw)%ivn(1:npoly))
        itab_wsfc(iw)%ecvec_vz(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vnz(itab_wsfc(iw)%ivn(1:npoly))
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
     vsprd = .10 * sqrt(sfcg%area(nwsfc))

     do iv = 2, nvsfc

        im1 = itab_vsfc(iv)%imn(1)
        im2 = itab_vsfc(iv)%imn(2)

        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)

        call oplot_transform(1,sfcg%xem(im1),sfcg%yem(im1),sfcg%zem(im1),xm1,ym1)
        call oplot_transform(1,sfcg%xem(im2),sfcg%yem(im2),sfcg%zem(im2),xm2,ym2)
        call oplot_transform(1,sfcg%xew(iw2),sfcg%yew(iw2),sfcg%zew(iw2),xw2,yw2)
        call oplot_transform(1,sfcg%xew(iw1),sfcg%yew(iw1),sfcg%zew(iw1),xw1,yw1)

        call trunc_segment(xm1,xm2,ym1,ym2,xq1,xq2,yq1,yq2,iskipm)
        call trunc_segment(xw1,xw2,yw1,yw2,xp1,xp2,yp1,yp2,iskipw)

        if (iskipm == 0) then
           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)
        endif

        if (iskipw == 0) then
           call o_frstpt (xp1,yp1)
           call o_vector (xp2,yp2)
        endif

        if ( xm1 < op%xmin .or.  &
             xm1 > op%xmax .or.  &
             ym1 < op%ymin .or.  &
             ym1 > op%ymax ) cycle

        write(string,'(I0)') im1
        call o_plchlq (xm1,ym1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') im2
        call o_plchlq (xm2,ym2,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw1
        call o_plchlq (xw1,yw1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw2
        call o_plchlq (xw2,yw2,trim(adjustl(string)),psiz,0.,0.)

     enddo  ! IV

     call o_frame()
     call o_clswk()

  endif

end subroutine grid_geometry_hex_sfc

!===============================================================================

subroutine sfc_atm_hex_overlay()

  use mem_grid,   only: nwa, nma, topm, topw, glatw, glonw, glonm, glatm, &
                        xew, yew, zew, xem, yem, zem
  use mem_sfcg,   only: nwsfc
  use ll_bins,    only: latlon_bins, itab_w0
  use mem_ijtabs, only: itab_w
  use misc_coms,  only: io6

  implicit none

  integer :: iw, iwsfc, im, j
  real    :: ds2

  integer, allocatable :: icountw(:)
  integer, allocatable :: icountm(:)
  real,    allocatable :: dsw_max(:)

  write(io6,*) "Creating latlon bins for atm/sfc overlay"

  ! Allocate and fill itab_w0 data structure with ATM GRID values that will be
  ! used inside subroutine latlon_bins

  allocate(itab_w0(nwa))
  do iw = 1,nwa
     itab_w0(iw)%npoly = itab_w(iw)%npoly
     itab_w0(iw)%iw(:) = itab_w(iw)%iw(:)
  enddo

  ! Allocate and fill latlon bins at different resolutions to compartmentalize
  ! all IW points of the ATM grid

  call latlon_bins(nwa, glatw, glonw)

  deallocate(itab_w0)

  ! Allocate arrays to record ATM IW and IM points that will have topography
  ! height filled in overlay_2 subroutine

  allocate(icountw(nwa)) ; icountw = 0
  allocate(icountm(nma)) ; icountm = 0

  ! Loop over all sfcgrid IW cells and for each, compute its overlay with
  ! ATM grid cells.  Also define topography of ATM IW and IM points.

  topw(:) = 0.0
  topm(:) = 0.0

  write(io6,*) "Computing atm/sfc overlaps"

  allocate(dsw_max(nwa))

  !$omp parallel
  !$omp do private(ds2)
  do iw = 2, nwa
     ds2 = 0.
     do j = 1, itab_w(iw)%npoly
        im = itab_w(iw)%im(j)
        ds2 = max(ds2, (xew(iw)-xem(im))**2 &
                     + (yew(iw)-yem(im))**2 &
                     + (zew(iw)-zem(im))**2 )
     enddo
     dsw_max(iw) = 1.0001 * sqrt(ds2)
  enddo
  !$omp end do

  !$omp do schedule(guided)
  do iwsfc = 2, nwsfc
     call sfc_atm_hex_overlay_2( icountw, icountm, iwsfc, dsw_max )
  enddo
  !$omp end do

  !$omp do
  do iw = 2, nwa

     if (icountw(iw) == 0) then

        write(*,*) "topw not set at iw = ", iw, glatw(iw), glonw(iw)
        stop 'error in sfc_atm_hex_overlay'

     elseif (icountw(iw) == 2) then

        topw(iw) = 0.5 * topw(iw)

     elseif (icountw(iw) == 3) then

        topw(iw) = topw(iw) / 3.

     elseif (icountw(iw) > 3) then

        topw(iw) = topw(iw) / real(icountw(iw))

        write(*,*) "illegal value for topw at iw = ", iw, glatw(iw), glonw(iw)
        stop 'error in sfc_atm_hex_overlay'

     endif

  enddo
  !$omp end do

  !$omp do
  do im = 2, nma

     if (icountm(im) == 0) then

        write(*,*) "topm not set at im = ", im, glatm(im), glonm(im)
        stop 'error in sfc_atm_hex_overlay'

     elseif (icountm(im) == 2) then

        topm(im) = 0.5 * topm(im)

     elseif (icountm(im) == 3) then

        topm(im) = topm(im) / 3.

     elseif (icountm(im) > 3) then

        topm(im) = topm(im) / real(icountm(im))

        write(*,*) "illegal value for topm at im = ", im
        stop 'error in sfc_atm_hex_overlay'

     endif

  enddo
  !$omp end do
  !$omp end parallel

end subroutine sfc_atm_hex_overlay

!============================================================================

subroutine sfc_atm_hex_overlay_2( icountw, icountm, iwsfc, dswmax )

  use mem_grid,    only: nwa, nma, xem, yem, zem, xew, yew, zew, arw0, &
                         topm, topw
  use mem_ijtabs,  only: itab_w, itab_m
  use mem_sfcg,    only: itab_wsfc, sfcg
  use consts_coms, only: r8, eradi, pio180
  use ll_bins,     only: latlon_bins, delat, bset

  implicit none

  integer, intent(in)    :: iwsfc
  integer, intent(inout) :: icountw(nwa)
  integer, intent(inout) :: icountm(nma)
  real,    intent(in)    :: dswmax (nwa)

  integer, parameter :: npmax = 7  ! Heptagons are max polygon for SINGLE GRID LEVEL
                                   ! in atm polygon cell

  integer, parameter :: nqmax = 7  ! Land cells can also be up to heptagons at this stage

  integer, parameter :: npqmax  = max(npmax,nqmax)
  real,    parameter :: oneplus = 1.0 + 5. * epsilon(1.)

  integer  :: iw, npoly, nsfcpoly, jmsfc, imsfc, jm, im, nwatm, idum
  real     :: xm0(npmax), xs0(nqmax), xw0, dum
  real     :: ym0(nqmax), ys0(nqmax), yw0

  real(r8) :: xm(npmax), xs(npmax), xw
  real(r8) :: ym(nqmax), ys(nqmax), yw
  real(r8) :: alpham(npmax), alphas(nqmax), alphaw

  real :: area, dsmax
  real :: sinwslat, coswslat
  real :: sinwslon, coswslon
  real :: raxis, raxisi
  real :: dxe(npqmax), dye(npqmax), dze(npqmax)
  real :: dxew, dyew, dzew
  real :: hlatw, hlonw

  integer :: iset, i, j, jw
  integer :: nwbin
  integer, allocatable :: iwbin(:)

  ! Evaluate xs,ys coordinates of SFCGRID iwsfc cell M points on gnomonic plane
  ! tangent at iwsfc

  raxis = sqrt( sfcg%xew(iwsfc) ** 2 + sfcg%yew(iwsfc) ** 2 )

  sinwslat = sfcg%zew(iwsfc) * eradi
  coswslat = raxis           * eradi

  ! For points less than 100 m from Earth's polar axis, make arbitrary
  ! assumption that longitude = 0 deg.  This is just to settle on a
  ! gnomonic tangent planar coordinate system in which to do the algebra.

  if (raxis >= 1.e2) then
     raxisi = 1.0 / raxis
     sinwslon = sfcg%yew(iwsfc) * raxisi
     coswslon = sfcg%xew(iwsfc) * raxisi
  else
     sinwslon = 0.
     coswslon = 1.
  endif

  ! Loop over all neighbor M points of this iwsfc

  nsfcpoly = itab_wsfc(iwsfc)%npoly
  dsmax    = 0.0

  do jmsfc = 1, nsfcpoly
     imsfc = itab_wsfc(iwsfc)%imn(jmsfc)

     dxe(jmsfc) = sfcg%xem(imsfc) - sfcg%xew(iwsfc)
     dye(jmsfc) = sfcg%yem(imsfc) - sfcg%yew(iwsfc)
     dze(jmsfc) = sfcg%zem(imsfc) - sfcg%zew(iwsfc)

     ! Very slightly increase sfc cell distances to avoid precision issues.
     ! This helps to ensure that an atmospheric W/M point that falls exactly
     ! on a surface cell boundary is matched with that surface cell.

     dxe(jmsfc) = oneplus * dxe(jmsfc)
     dye(jmsfc) = oneplus * dye(jmsfc)
     dze(jmsfc) = oneplus * dze(jmsfc)

     dsmax = max(dsmax, dxe(jmsfc)**2 + dye(jmsfc)**2 + dze(jmsfc)**2)
  enddo

  call de_gn_mult(nsfcpoly,dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xs0,ys0)

  xs(1:nsfcpoly) = real(xs0(1:nsfcpoly),r8)
  ys(1:nsfcpoly) = real(ys0(1:nsfcpoly),r8)

  dsmax = sqrt(dsmax)

  ! Get ATM bin indices for this iwsfc point

  hlonw = max(-179.9999,min(179.9999,sfcg%glonw(iwsfc)))
  hlatw = max( -89.9999,min( 89.9999,sfcg%glatw(iwsfc)))

  do iset = 1,4
     j = int(delat(iset) * (hlatw +  90.)) + 1
     i = int(bset(iset)%blat(j)%delon * (hlonw + 180.)) + 1

     ! Check if bins are allocated, starting with smallest

     if (allocated(bset(iset)%blat(j)%bins(i)%iw)) then
        nwbin = bset(iset)%blat(j)%bins(i)%nw
        allocate(iwbin(nwbin))
        iwbin(:) = bset(iset)%blat(j)%bins(i)%iw(:)
        go to 10
     endif
  enddo  ! iset

  ! If iset loop has not branched to statement 10, search over all IW points

  print*, 'No allocated ATM bin found for iwsfc = ', iwsfc, hlatw, hlonw
  print*, 'Conducting search over all IW points '

  nwbin = nwa - 1
  allocate(iwbin(nwbin))
  do iw = 2,nwa
     iwbin(iw-1) = iw
  enddo

  10 continue

  ! Loop over all ATM IW points in bin

  nwatm = 0

  do jw = 1, nwbin
     iw = iwbin(jw)

     dum = sqrt( (xew(iw) - sfcg%xew(iwsfc))**2 &
               + (yew(iw) - sfcg%yew(iwsfc))**2 &
               + (zew(iw) - sfcg%zew(iwsfc))**2 )

     ! Skip this overlap if atm and sfc cells are too far apart

     if (dum > dsmax + dswmax(iw)) cycle

     ! Loop over all neighbor M points of this IW

     npoly = itab_w(iw)%npoly

     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        ! Evaluate xm,ym coordinates of current M point on gnomonic tangent
        ! plane tangent at iwsfc

        dxe(jm) = xem(im) - sfcg%xew(iwsfc)
        dye(jm) = yem(im) - sfcg%yew(iwsfc)
        dze(jm) = zem(im) - sfcg%zew(iwsfc)
     enddo

     call de_gn_mult(npoly,dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xm0,ym0)

     xm(1:npoly) = real(xm0(1:npoly),r8)
     ym(1:npoly) = real(ym0(1:npoly),r8)

     ! Evaluate possible overlap of ATM and SURFACE polygons

     call polygon_overlap2(nsfcpoly,npoly,xs,ys,xm,ym,sfcg%area(iwsfc),arw0(iw),area,alphas,alpham)

     ! Set topm of atmospheric IM points that are inside or on the boundary
     ! of this surface cell

     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        ! only check each im point once
        if (itab_m(im)%iw(1) == iw) then

           if (alpham(jm) > 1.0_r8) then

              !$omp atomic
              icountm(im) = icountm(im) + 1

              !$omp atomic
              topm(im) = topm(im) + sfcg%topw(iwsfc)

           endif
        endif
     enddo

     ! Skip further computation if overlap is small.

     if (area < 1.e-5 * sfcg%area(iwsfc)) cycle

     ! Evaluate x,y coordinates of current IW point on gnomonic plane
     ! tangent at iwsfc

     dxew = xew(iw) - sfcg%xew(iwsfc)
     dyew = yew(iw) - sfcg%yew(iwsfc)
     dzew = zew(iw) - sfcg%zew(iwsfc)

     call de_gn(dxew,dyew,dzew,coswslat,sinwslat,coswslon,sinwslon,xw0,yw0)

     xw = real(xw0,r8)
     yw = real(yw0,r8)

     ! If overlap is positive, even though it may be very small, check whether
     ! this SFCG cell overlaps with the IW point of this ATM cell. Where overlap
     ! is found, set topo height of ATM cell point from TOPW of SURFACE cell

     call inout_check(nsfcpoly,xs,ys,xw,yw,alphaw)

     if (alphaw > 1.0_r8) then

        !$omp atomic
        icountw(iw) = icountw(iw) + 1

        !$omp atomic
        topw(iw) = topw(iw) + sfcg%topw(iwsfc)

     endif

     ! If overlap area is less than 1.e-5 of sfcg cell area or 2.e-6 of
     ! atmospheric cell area, this overlap will not be counted.

     if (area < 2.0e-6 * arw0(iw)) cycle

     ! This iwsfc SURFACE cell overlaps with IW ATM cell.

     nwatm = nwatm + 1

     itab_wsfc(iwsfc)%nwatm = nwatm
     itab_wsfc(iwsfc)%iwatm(nwatm) = iw
     itab_wsfc(iwsfc)%arc  (nwatm) = area

  enddo  ! jw/iwbin

  ! Order multiple overlap areas so that largest is first

  if (nwatm > 1) then
     do j = 2,nwatm
        if (itab_wsfc(iwsfc)%arc(j)  > itab_wsfc(iwsfc)%arc(1)) then
           idum                      = itab_wsfc(iwsfc)%iwatm(j)
           itab_wsfc(iwsfc)%iwatm(j) = itab_wsfc(iwsfc)%iwatm(1)
           itab_wsfc(iwsfc)%iwatm(1) = idum

           dum                     = itab_wsfc(iwsfc)%arc(j)
           itab_wsfc(iwsfc)%arc(j) = itab_wsfc(iwsfc)%arc(1)
           itab_wsfc(iwsfc)%arc(1) = dum
        endif
     enddo
  endif

  end subroutine sfc_atm_hex_overlay_2
