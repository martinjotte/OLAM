!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University;
   ! Colorado State University Research Foundation ; ATMET, LLC

   ! This software is free software; you can redistribute it and/or modify it
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version.

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.

   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
   ! (http://www.gnu.org/licenses/gpl.html)
   !----------------------------------------------------------------------------

!===============================================================================
subroutine voronoi_sfc()

  use mem_sfcg,     only: nmsfc, nvsfc, nwsfc, itab_msfc, itab_vsfc, &
                          itab_wsfc, sfcg, alloc_sfcgrid1, im_orig

  use mem_delaunay, only: itab_md, itab_ud, itab_wd, nmd, nud, nwd, &
                          xemd, yemd, zemd, iwdorig

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

  call move_alloc(iwdorig, im_orig)

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
        raxisi = 1.0 / raxis

        sinwslat = zebc  * eradi
        coswslat = raxis * eradi

        sinwslon = yebc * raxisi
        coswslon = xebc * raxisi

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
  use mem_para,    only: myrank

  implicit none

  integer       :: im,iv,iw,ivn
  integer       :: im1,im2
  integer       :: iw1,iw2
  integer       :: j,npoly,j1,j2
  real          :: raxis, expansion
  real          :: dvm1,dvm2
  real          :: xm1,xm2,xv,ym1,ym2,yv,frac,alpha
  real          :: xw1,xw2,yw1,yw2
  real          :: xq1, yq1, xq2, yq2, psiz, vsprd
  integer       :: iskip
  real          :: quarter_kite(2,nvsfc)
  character(10) :: string

  integer               :: lwork
  integer               :: info
  real(r8)              :: b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7), vrot_x(7), vrot_y(7)
  real(r8), allocatable :: work(:)
  real(r8), allocatable :: a(:,:)
  real(r8)              :: wsize(1), vdotw, vmag, fact

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

  !$omp parallel

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

  !$omp end parallel

  ! Coefficients for converting earth-cartesian velocity to V and W

  if (mdomain < 2 .or. mdomain == 5) then

 !!    !$omp do private(npoly, fo, a, b, work, info, j, ivn, vdotw, vmag, fact, &
 !!    !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y, wsize, lwork)
     do iw = 2, nwsfc
      
        npoly = itab_wsfc(iw)%npoly

        ! Default coefficients from Perot
        fo(1:npoly) = 2.0_r8 * itab_wsfc(iw)%farv(1:npoly)

        if (allocated(a)) then
           if (size(a,2) /= npoly) deallocate(a)
        endif

        if (.not. allocated(a)) allocate(a(3,npoly))

        if (mdomain < 2) then

           do j = 1, npoly
              ivn = itab_wsfc(iw)%ivn(j)

              ! Compute the components of the V unit normals perpendicular to W

              vdotw = sfcg%vnx(ivn)*sfcg%wnx(iw) + sfcg%vny(ivn)*sfcg%wny(iw) + sfcg%vnz(ivn)*sfcg%wnz(iw)

              vnx_ps(j) = sfcg%vnx(ivn) - vdotw * sfcg%wnx(iw)
              vny_ps(j) = sfcg%vny(ivn) - vdotw * sfcg%wny(iw)
              vnz_ps(j) = sfcg%vnz(ivn) - vdotw * sfcg%wnz(iw)

              ! Normalize these new vectors to unit length

              vmag = sqrt( vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2 )

              vnx_ps(j) = vnx_ps(j) / vmag
              vny_ps(j) = vny_ps(j) / vmag
              vnz_ps(j) = vnz_ps(j) / vmag

              ! Rotate these new unit normals to a coordinate system with Z aligned with W

              if (sfcg%wnz(iw) >= 0.0) then

                 fact = ( sfcg%wny(iw)*vnx_ps(j) - sfcg%wnx(iw)*vny_ps(j) ) / ( 1.0 + sfcg%wnz(iw) )

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

        b(1) = 1._r8 - sum( fo(1:npoly) * a(1,:) )
        b(2) = 1._r8 - sum( fo(1:npoly) * a(2,:) )
        b(3) =       - sum( fo(1:npoly) * a(3,:) )

        call dgels( 'N', 3, npoly, 1, a, 3, b, 7, wsize, -1, info )
        lwork = nint(wsize(1)) + 1

        if (allocated(work) .and. size(work) < lwork) deallocate(work)
        if (.not. allocated(work)) allocate(work(lwork))

        call dgels( 'N', 3, npoly, 1, a, 3, b, 7, work, size(work), info )

        ! Vector b is now the correction to the coefficients fo
        b(1:npoly) = b(1:npoly) + fo(1:npoly)

        if (info == 0 .and. all(b(1:npoly) > 0.05_r8) .and. all(b(1:npoly) < 0.7_r8)) then

           ! itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly)   ! This is old ATM version
           ! itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly)   ! This is old ATM version
           ! itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly)   ! This is old ATM version

           fact =  sum( b(1:npoly) * vnx_ps(1:npoly) * sfcg%vnx(itab_wsfc(iw)%ivn(1:npoly)) )
           fact = (1._r8 - sfcg%wnx(iw)**2) / max(fact, 1.e-30_r8)

           itab_wsfc(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly) * fact

           fact = sum( b(1:npoly) * vny_ps(1:npoly) * sfcg%vny(itab_wsfc(iw)%ivn(1:npoly)) )
           fact = (1._r8 - sfcg%wny(iw)**2) / max(fact, 1.e-30_r8)

           itab_wsfc(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly) * fact

           fact = sum( b(1:npoly) * vnz_ps(1:npoly) * sfcg%vnz(itab_wsfc(iw)%ivn(1:npoly)) )
           fact = (1._r8 - sfcg%wnz(iw)**2) / max(fact, 1.e-30_r8)

           itab_wsfc(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly) * fact

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
 !!    !omp end do
   
     if (allocated(a))    deallocate(a)
     if (allocated(work)) deallocate(work)

  else

 !!    !$omp do private(npoly)
     do iw = 2, nwsfc
        npoly = itab_wsfc(iw)%npoly
        itab_wsfc(iw)%ecvec_vx(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vnx(itab_wsfc(iw)%ivn(1:npoly))
        itab_wsfc(iw)%ecvec_vy(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vny(itab_wsfc(iw)%ivn(1:npoly))
        itab_wsfc(iw)%ecvec_vz(1:npoly) = 2.0 * itab_wsfc(iw)%farv(1:npoly) &
                                        * sfcg%vnz(itab_wsfc(iw)%ivn(1:npoly))
     enddo
 !!    !$omp end do

  endif

  ! Plot grid lines

  if (.false.) then

     if (myrank /= 0) return

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

  use mem_grid,    only: nwa, nma, topm, topw, glatw, glonw, glonm, glatm
  use mem_sfcg,    only: nwsfc

  implicit none

  integer :: iw, iwsfc, im

  integer, allocatable :: icountw(:)
  integer, allocatable :: icountm(:)

  allocate(icountw(nwa)) ; icountw = 0
  allocate(icountm(nma)) ; icountm = 0

  !$omp parallel do
  do iwsfc = 2, nwsfc

     call sfc_atm_hex_overlay_2( icountw, icountm, iwsfc )

  enddo
  !$omp end parallel do

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
!       stop 'error in sfc_atm_hex_overlay'

     endif

  enddo

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
!       stop 'error in sfc_atm_hex_overlay'

     endif

  enddo

end subroutine sfc_atm_hex_overlay

!============================================================================

subroutine sfc_atm_hex_overlay_2( icountw, icountm, iwsfc )

  use mem_grid,    only: nwa, nma, xem, yem, zem, xew, yew, zew, arw0, &
                         topm, topw
  use mem_ijtabs,  only: itab_w, itab_m
  use mem_sfcg,    only: itab_wsfc, sfcg, im_orig
  use consts_coms, only: r8, eradi, pio180

  implicit none

  integer, intent(in)    :: iwsfc
  integer, intent(inout) :: icountw(nwa)
  integer, intent(inout) :: icountm(nma)

  integer, parameter :: npmax = 7  ! Heptagons are max polygon for SINGLE GRID LEVEL
                                   ! in atm polygon cell

  integer, parameter :: nqmax = 7  ! Land cells can also be up to heptagons at this stage

  real,    parameter :: oneplus = 1.0 + 5. * epsilon(1.)

  integer :: iw, npoly, nsfcpoly, jmsfc, imsfc, jm, im, nwatm, idum
  real    :: xp0, yp0, xq0, yq0, dum

  real(r8) :: xw, yw, alpha
  real(r8) :: xp(npmax),yp(npmax)
  real(r8) :: xq(nqmax),yq(nqmax)
  real(r8) :: alphap(npmax)
  real(r8) :: alphaq(nqmax)

  real :: area
  real :: sinwslat, coswslat
  real :: sinwslon, coswslon
  real :: raxis, raxisi
  real :: dxe, dye, dze

  integer :: nw, imatm(7), iwatm(8), iwnew, jw, np, j

  np          = itab_wsfc(iwsfc)%npoly
  imatm(1:np) = im_orig( itab_wsfc(iwsfc)%imn(1:np) )

  nw         = 3
  iwatm(1:3) = itab_m( imatm(1) )%iw(1:3)

  do jm = 2, itab_wsfc(iwsfc)%npoly
     if ( any(imatm(jm) == imatm(1:jm-1)) ) cycle

     do jw = 1, 3
        iwnew = itab_m( imatm(jm) )%iw(jw)
        if ( any(iwnew == iwatm(1:nw)) ) cycle

        nw = nw + 1
        if (nw > 8) stop 'too many points!'
        iwatm(nw) = iwnew
     enddo
  enddo

  nwatm = 0
  nsfcpoly = itab_wsfc(iwsfc)%npoly

  raxis  = sqrt( sfcg%xew(iwsfc) ** 2 + sfcg%yew(iwsfc) ** 2 )
  raxisi = 1.0 / raxis

  sinwslat = sfcg%zew(iwsfc) * eradi
  coswslat = raxis           * eradi

  sinwslon = sfcg%yew(iwsfc) * raxisi
  coswslon = sfcg%xew(iwsfc) * raxisi

  ! Loop over all neighbor M points of this iwsfc

  do jmsfc = 1, nsfcpoly
     imsfc = itab_wsfc(iwsfc)%imn(jmsfc)

     ! Evaluate x,y coordinates of LAND cell M points on gnomic
     ! plane tangent at iwsfc

     dxe = sfcg%xem(imsfc) - sfcg%xew(iwsfc)
     dye = sfcg%yem(imsfc) - sfcg%yew(iwsfc)
     dze = sfcg%zem(imsfc) - sfcg%zew(iwsfc)

     ! Very slightly increase sfc cell distances to avoid precision issues.
     ! This helps to ensure that an atmospheric W/M point that falls exactly
     ! on a surface cell boundary is matched with that surface cell.
     dxe = oneplus * dxe
     dye = oneplus * dye
     dze = oneplus * dze

     call de_gn(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xq0,yq0)

     xq(jmsfc) = real(xq0,r8)
     yq(jmsfc) = real(yq0,r8)

  enddo

  do jw = 1, nw
     iw = iwatm(jw)

           ! Loop over all neighbor M points of this IW

           npoly = itab_w(iw)%npoly

           do jm = 1,npoly
              im = itab_w(iw)%im(jm)

              ! Evaluate x,y coordinates of current M point on gnomic
              ! plane tangent at iwsfc

              dxe = xem(im) - sfcg%xew(iwsfc)
              dye = yem(im) - sfcg%yew(iwsfc)
              dze = zem(im) - sfcg%zew(iwsfc)

              call de_gn(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xp0,yp0)

              xp(jm) = real(xp0,r8)
              yp(jm) = real(yp0,r8)
           enddo

           ! Evaluate possible overlap of ATM and SURFACE polygons

           call polygon_overlap2(npoly,nsfcpoly,xp,yp,xq,yq,area,alphap,alphaq)

           ! Set topm of atmospheric IM points that are inside or on the boundary
           ! of this surface cell

           do jm = 1,npoly
              im = itab_w(iw)%im(jm)

              ! only check each im point once
              if (itab_m(im)%iw(1) == iw) then

                 if (alphap(jm) > 1.0_r8) then

                    !$omp atomic
                    icountm(im) = icountm(im) + 1

                    !$omp atomic
                    topm(im) = topm(im) + sfcg%topw(iwsfc)

                 endif
              endif

           enddo

           ! Skip further computation if overlap is zero

           if (area < 1.0e-7) cycle

           ! Evaluate x,y coordinates of current IW point on gnomic plane
           ! tangent at iwsfc

           dxe = xew(iw) - sfcg%xew(iwsfc)
           dye = yew(iw) - sfcg%yew(iwsfc)
           dze = zew(iw) - sfcg%zew(iwsfc)

           call de_gn(dxe,dye,dze,coswslat,sinwslat,coswslon,sinwslon,xp0,yp0)

           xw = real(xp0,r8)
           yw = real(yp0,r8)

           ! If overlap is positive, even though it may be very small, check whether
           ! this SFCG cell overlaps with the IW point of this ATM cell. Where overlap
           ! is found, set topo height of ATM cell point from TOPW of SURFACE cell

           call inout_check(nsfcpoly,xq,yq,xw,yw,alpha)
           if (alpha > 1.0_r8) then

              !$omp atomic
              icountw(iw) = icountw(iw) + 1

              !$omp atomic
              topw(iw) = topw(iw) + sfcg%topw(iwsfc)

           endif

           ! If overlap area is less than 1.e-5 of sfcg cell area or 2.e-6 of
           ! atmospheric cell area, this overlap will not be counted.

           if (area < 1.0e-5 * sfcg%area(iwsfc)) cycle
           if (area < 2.0e-6 * arw0(iw)) cycle

           ! This iwsfc SURFACE cell overlaps with IW ATM cell.

           nwatm = nwatm + 1

           itab_wsfc(iwsfc)%nwatm = nwatm
           itab_wsfc(iwsfc)%iwatm(nwatm) = iw
           itab_wsfc(iwsfc)%arc  (nwatm) = area

!        enddo  ! iw
!     enddo  ! ibin
  enddo  ! jbin

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

!============================================================================

subroutine polygon_overlap2(np,nq,xp,yp,xq,yq,area,alphap,alphaq)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use consts_coms, only: r8

  implicit none

  integer,  intent(in) :: np,nq ! Number of vertices in p and q
  real(r8), intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
  real(r8), intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

  real,     intent(out) :: area        ! area of overlap of p and q
  real(r8), intent(out) :: alphap(np)  ! if any perimeter points of p are on/inside q
  real(r8), intent(out) :: alphaq(nq)  ! if any perimeter points of q are on/inside p

  real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer :: nev      ! # of events
  integer :: nsect    ! # of intersections between strip centerline and p,q edges
  real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
  real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges
  real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges
  real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges
  real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values

  integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real(r8) :: p0,q0,dx,dy,dxtrap

  real(r8) :: alpha    ! angle swept by ray from point to polygon perimeter during circuit

  ! Find vertices in p that are not outside q

  nev = 0
  area = 0.

  do ip = 1,np
     call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

     alphap(ip) = alpha

     if (abs(alpha) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  ! Find vertices in q that are not outside p

  do iq = 1,nq
     call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

     alphaq(iq) = alpha

     if (abs(alpha) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  ! Find intersecting edges of polygons p and q

  do ipa = 1,np
     ipb = ipa + 1
     if (ipa == np) ipb = 1

     do iqa = 1,nq
        iqb = iqa + 1
        if (iqa == nq) iqb = 1

        call intersect(xp(ipa),yp(ipa),xp(ipb),yp(ipb)  &
                      ,xq(iqa),yq(iqa),xq(iqb),yq(iqb),p0,q0,iflag)

        if (iflag == -1) cycle
        if (p0 < -0.000000001_r8 .or. p0 > 1.000000001_r8) cycle
        if (q0 < -0.000000001_r8 .or. q0 > 1.000000001_r8) cycle

  ! Line segments pa-pb and qa-qb intersect; find y-coordinate of intersection

        nev = nev + 1
        yev(nev) = (1.0_r8 - p0) * yp(ipa) + p0 * yp(ipb)
     enddo
  enddo

  ! Sort event list to obtain increasing order of yev

  call fltsort(nev,yev)

  ! Loop over event points

  do iev = 1,nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)

     if (dy < 0.1_r8) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5_r8 * dy
     dx = 0.0_r8

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1,np
        ib = ia + 1
        if (ia == np) ib = 1

        if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

        nsect = nsect + 1
        xmid(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) / (yp(ib) - yp(ia))
        xbot(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) / (yp(ib) - yp(ia))
        xtop(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) / (yp(ib) - yp(ia))

     enddo

     do ia = 1,nq
        ib = ia + 1
        if (ia == nq) ib = 1

        if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

        nsect = nsect + 1
        xmid(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) / (yq(ib) - yq(ia))
        xbot(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) / (yq(ib) - yq(ia))
        xtop(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) / (yq(ib) - yq(ia))
     enddo

     ! Sort xmid values into increasing order

     call fltsort3(nsect,xmid,xbot,xtop)

     ! see if the segment is inside both polygons

     do is = 1,nsect - 1
        xcent = 0.5_r8 * (xmid(is) + xmid(is+1))

        if (xcent == xmid(is)) cycle          ! if segment length is 0

        call inout_check(np,xp,yp,xcent,ymid,alpha)
        if (abs(alpha) < 1.0_r8) cycle

        call inout_check(nq,xq,yq,xcent,ymid,alpha)
        if (abs(alpha) < 1.0_r8) cycle

        dxtrap = xmid(is+1) - xmid(is)

        dx = dx + dxtrap
     enddo

     area = area + dx * dy

  enddo

end subroutine polygon_overlap2
