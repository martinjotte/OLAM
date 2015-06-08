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
subroutine voronoi()

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd, itab_m, itab_v, itab_w, &
                       mrls, alloc_itabs, mloops

use mem_grid,    only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, &
                       xem, yem, zem, xew, yew, zew, &
                       alloc_xyzem, alloc_xyzew
use misc_coms,   only: io6, mdomain
use consts_coms, only: pi2, erad, erador5
use oname_coms,  only: nl

implicit none

integer :: im1,im2
integer :: iw1,iw2,iw3,im,iw
integer :: nmad,nuad,nwad
integer :: iwd,iv,iud,iud1,iud2,imd,npoly,j,j1
real    :: expansion

! Interchange grid dimensions

nmad = nma
nuad = nua
nwad = nwa

nma = nwad
nua = nuad
nva = nuad
nwa = nmad

mma = nma
mua = nua
mva = nva
mwa = nwa

! Allocate Voronoi set of itabs

call alloc_itabs(nma,nva,nwa,0)

! Allocate XEW,YEW,ZEW arrays, and fill their values from XEM,YEM,ZEM, which
! still have the OLD nmad dimension which is the NEW nwa dimension

call alloc_xyzew(nwa)

do iw = 1,nwa
   xew(iw) = xem(iw)
   yew(iw) = yem(iw)
   zew(iw) = zem(iw)
enddo

! Deallocate XEM,YEM,ZEM and re-allocate to NEW nma dimension

deallocate (xem,yem,zem)
call alloc_xyzem(nma)

! Since XEM,YEM,ZEM have just been re-allocated, initialize their values to be
! barycenters of Delaunay triangles whose vertices are at XEW,YEW,ZEW

do iwd = 2,nwad
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
   itab_v(iv)%ivglobe  = itab_ud(iud)%iuglobe
   itab_v(iv)%mrlv     = itab_ud(iud)%mrlu

   itab_v(iv)%im(1:6)  = itab_ud(iud)%iw(1:6)
   itab_v(iv)%iw(1:2)  = itab_ud(iud)%im(1:2)

   itab_v(iv)%iv(1:4)  = itab_ud(iud)%iu(1:4)
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

   imd = itab_ud(iud)%im(1)
   npoly = itab_md(imd)%npoly

   do j = 1,npoly
      j1 = j + 1
      if (j == npoly) j1 = 1

      iud1 = itab_md(imd)%iu(j)
      iud2 = itab_md(imd)%iu(j1)

! IW(3) and IW(4) neighbors of IV   

      if (iud2 == iv) then
         iw1 = itab_ud(iud1)%im(1)
         iw2 = itab_ud(iud1)%im(2)

         if (iw1 == imd) then
            itab_v(iv)%iw(3) = iw2
         else
            itab_v(iv)%iw(3) = iw1
         endif
      endif

      if (iud1 == iv) then
         iw1 = itab_ud(iud2)%im(1)
         iw2 = itab_ud(iud2)%im(2)

         if (iw1 == imd) then
            itab_v(iv)%iw(4) = iw2
         else
            itab_v(iv)%iw(4) = iw1
         endif
      endif

   enddo           

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
   itab_m(im)%imglobe   = itab_wd(iwd)%iwglobe

   itab_m(im)%mrlm      = itab_wd(iwd)%mrlw
   itab_m(im)%ngr       = itab_wd(iwd)%ngr

   itab_m(im)%mrlm_orig = itab_wd(iwd)%mrlw_orig
   itab_m(im)%mrow      = itab_wd(iwd)%mrow

   itab_m(im)%iv(1:3)   = itab_wd(iwd)%iu(1:3)
   itab_m(im)%iw(1:3)   = itab_wd(iwd)%im(1:3)
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

return
end subroutine voronoi

!===============================================================================

subroutine pcvt()

! Iterative procedure for defining centroidal voronoi cells

use mem_ijtabs,  only: itab_m, itab_w
use mem_grid,    only: nma, nwa, xem, yem, zem, xew, yew, zew
use consts_coms, only: erad, piu180
use misc_coms,   only: io6, mdomain

implicit none

integer, parameter :: niter = 1
real, parameter :: fracw = 1., fracm = 1.

integer :: jm,jm1,iw,im,iw1,iw2,iw3,iter,npoly

real :: xm(7),ym(7)
real :: raxis,glatw,glonw,area,xc,yc,xec,yec,zec,expansion
real :: xebc,yebc,zebc
real :: glatbc,glonbc
real :: x1,x2,x3,y1,y2,y3
real :: dx12,dx13,dx23
real :: s1,s2,s3
real :: xcc,ycc

! Main iteration loop

do iter = 1,niter

! Loop over all M points

   do im = 2,nma

! Indices of 3 W points surrounding M point

      if (any(itab_m(im)%iw(1:3) < 2)) cycle

      iw1 = itab_m(im)%iw(1)
      iw2 = itab_m(im)%iw(2)
      iw3 = itab_m(im)%iw(3)

! Fill circumcentric M point coordinates from coordinates of 3 W points.
! This establishes W cell as voronoi.

! For global domain, transform from sphere to PS plane

      if (mdomain <= 1) then

! First, compute barycenter of 3 W points

         xebc = (xew(iw1) + xew(iw2) + xew(iw3)) / 3.
         yebc = (yew(iw1) + yew(iw2) + yew(iw3)) / 3.
         zebc = (zew(iw1) + zew(iw2) + zew(iw3)) / 3.

! Push barycentric point coordinates out to earth radius

         expansion = erad / sqrt(xebc ** 2  &
                               + yebc ** 2  &
                               + zebc ** 2  )

         xebc = xebc * expansion
         yebc = yebc * expansion
         zebc = zebc * expansion

! Get latitude and longitude of barycentric point

         raxis = sqrt(xebc ** 2 + yebc ** 2)

         glatbc = atan2(zebc,raxis) * piu180
         glonbc = atan2(yebc,xebc) * piu180

! Transform 3 W points to PS coordinates

         call e_ps(xew(iw1),yew(iw1),zew(iw1),glatbc,glonbc,x1,y1)
         call e_ps(xew(iw2),yew(iw2),zew(iw2),glatbc,glonbc,x2,y2)
         call e_ps(xew(iw3),yew(iw3),zew(iw3),glatbc,glonbc,x3,y3)

! For Cartesian domain, use given planar X,Y coordinates

      else

! First, compute barycenter of 3 W points

         xebc = (xew(iw1) + xew(iw2) + xew(iw3)) / 3.
         yebc = (yew(iw1) + yew(iw2) + yew(iw3)) / 3.
         zebc = (zew(iw1) + zew(iw2) + zew(iw3)) / 3.

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

! For global domain, transform circumcenter from PS to earth coordinates

      if (mdomain <= 1) then
         call ps_e(xem(im),yem(im),zem(im),glatbc,glonbc,xcc,ycc)
      else
         xem(im) = xcc + xebc
         yem(im) = ycc + yebc
      endif

   enddo

   if (iter == niter) exit

! Loop over all W points

   do iw = 2,nwa

! Determine latitude and longitude of W point

      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

      glatw = atan2(zew(iw),raxis)   * piu180
      glonw = atan2(yew(iw),xew(iw)) * piu180

! Loop over all M points adjacent to current W point

      npoly = itab_w(iw)%npoly

      do jm = 1,npoly
         im = itab_w(iw)%im(jm)

! Determine local PS coordinates of current M point

         if (mdomain <= 1) then
            call e_ps(xem(im),yem(im),zem(im),glatw,glonw,xm(jm),ym(jm))
         else
            xm(jm) = xem(im) - xew(iw)
            ym(jm) = yem(im) - yew(iw)
         endif
      enddo

! Compute Voronoi cell area for current W point

      area = 0.

      do jm = 1,npoly
         jm1 = jm + 1
         if (jm1 > npoly) jm1 = 1

         area = area + .5 * (xm(jm) * ym(jm1) - xm(jm1) * ym(jm))
      enddo

! Compute Voronoi cell centroid (in local PS coordinates) for current W point

      xc = 0.
      yc = 0.

      do jm = 1,npoly
         jm1 = jm + 1
         if (jm1 > npoly) jm1 = 1

         xc = xc + (xm(jm) + xm(jm1)) * (xm(jm) * ym(jm1) - xm(jm1) * ym(jm))
         yc = yc + (ym(jm) + ym(jm1)) * (xm(jm) * ym(jm1) - xm(jm1) * ym(jm))
      enddo

      xc = xc / (6. * area)
      yc = yc / (6. * area)

! Transform local PS coordinates of centroid to Earth coordinates

      if (mdomain <= 1) then
         call ps_e(xec,yec,zec,glatw,glonw,xc,yc)
      else
         xec = xc + xew(iw)
         yec = yc + yew(iw)
         zec = 0.
      endif

! Move current W point (fractionally) to centroid location

      xew(iw) = (1. - fracw) * xew(iw) + fracw * xec
      yew(iw) = (1. - fracw) * yew(iw) + fracw * yec
      zew(iw) = (1. - fracw) * zew(iw) + fracw * zec

! Adjust W point coordinates to earth radius in case required

      expansion = erad / sqrt(xew(iw) ** 2 + yew(iw) ** 2 + zew(iw) ** 2)

      xew(iw) = xew(iw) * expansion
      yew(iw) = yew(iw) * expansion
      zew(iw) = zew(iw) * expansion

   enddo

enddo

return
end subroutine pcvt

!===============================================================================

subroutine grid_geometry_hex()

use mem_ijtabs,  only: itab_m, itab_v, itab_w
use mem_grid,    only: nza, nma, nva, nwa, xev, yev, zev, xem, yem, zem, &
                       xew, yew, zew, unx, uny, unz, wnx, wny, wnz,      &
                       vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, &
                       arm0, glonm, glatm, glatv, glonv
use misc_coms,   only: io6, mdomain, grdlat, grdlon, nxp, rinit
use consts_coms, only: erad, erad2, piu180, eradsq,pio2
use oplot_coms,  only: op
use oname_coms,  only: nl
use mem_para,    only: myrank

implicit none

integer            :: im,iv,iw
integer            :: im1,im2
integer            :: iw1,iw2
integer            :: j,npoly
real               :: expansion
real               :: raxis
real               :: dvm1,dvm2
integer            :: j1,j2
integer            :: npoly1,npoly2,np
real               :: xm1,xm2,xv,ym1,ym2,yv,frac,alpha
real               :: xw1,xw2,yw1,yw2
real               :: xq1, yq1, xq2, yq2, psiz, vsprd
integer            :: iskip, iwp, ivp, imp
logical            :: dops
real               :: quarter_kite(2,nva)
character(10)      :: string
integer, parameter :: lwork = 200
real               :: work(lwork)
integer            :: info
real               :: vdotw, vmag, fact
real               :: b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7), vrot_x(7), vrot_y(7)
real, allocatable  :: a(:,:)

! Loop over all M points

!$omp parallel

!$omp do private(raxis,npoly,j,iw)
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

   if (im1 > 1 .and. im2 > 1 .and. (frac < .001 .or. frac > .999)) then
      print*, 'Non-intersecting U-V edges detected in grid geometry'
      print*, 'FRAC  = ',frac
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

!$omp parallel
!$omp do private(raxis)
do iw = 2,nwa

! Fill global index (replaced later if this run is parallel)

   itab_w(iw)%iwglobe = iw

! Fill outward unit vector components and latitude and longitude of W point

   if (mdomain <= 1) then

      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

      glatw(iw) = atan2(zew(iw),raxis)   * piu180
      glonw(iw) = atan2(yew(iw),xew(iw)) * piu180
      
      wnx(iw) = xew(iw) / erad
      wny(iw) = yew(iw) / erad
      wnz(iw) = zew(iw) / erad

   else

      glatw(iw) = 0. ! want it this way?
      glonw(iw) = 0. ! want it this way?

      wnx(iw) = 0.0
      wny(iw) = 0.0
      wnz(iw) = 1.0

   endif
enddo
!$omp end do

!$omp do private(npoly,j2,j1,iv,iw1,iw2,dops,npoly1,npoly2,np,&
!$omp            im1,xw1,xw2,yw1,yw2,xv,yv,alpha)
do iw = 2,nwa

! Number of polygon edges/vertices

   npoly = itab_w(iw)%npoly

! Loop over all polygon edges

   do j2 = 1,npoly
      j1 = j2 - 1
      if (j2 == 1) j1 = npoly

      iv  = itab_w(iw)%iv(j2)
      iw2 = itab_w(iw)%iw(j2)
      iw1 = itab_w(iw)%iw(j1)

! Fractional area of arw0(iw) that is occupied by M and V sectors

      if (itab_v(iv)%iw(1) == iw1) then
         itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(1,iv) / arw0(iw)
         itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(2,iv) / arw0(iw)
      else
         itab_w(iw)%farm(j1) = itab_w(iw)%farm(j1) + quarter_kite(2,iv) / arw0(iw)
         itab_w(iw)%farm(j2) = itab_w(iw)%farm(j2) + quarter_kite(1,iv) / arw0(iw)
      endif
      itab_w(iw)%farv(j2) = (quarter_kite(1,iv) + quarter_kite(2,iv)) / arw0(iw)

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
         call e_ps(xew(iw1),yew(iw1),zew(iw1),glatw(iw),glonw(iw),xw1,yw1)
         call e_ps(xew(iw2),yew(iw2),zew(iw2),glatw(iw),glonw(iw),xw2,yw2)
         call e_ps(xev(iv),yev(iv),zev(iv),glatw(iw),glonw(iw),xv,yv)
      else
         xw1 = xew(iw1) - xew(iw)
         yw1 = yew(iw1) - yew(iw)
         xw2 = xew(iw2) - xew(iw)
         yw2 = yew(iw2) - yew(iw)
         xv  = xev(iv)  - xew(iw)
         yv  = yev(iv)  - yew(iw)
      endif

! Coefficients for eastward and northward components of gradient

      itab_w(iw)%gxps1(j1) =  yw2 / (xw1 * yw2 - xw2 * yw1)
      itab_w(iw)%gxps2(j1) = -yw1 / (xw1 * yw2 - xw2 * yw1)

      itab_w(iw)%gyps1(j1) = -xw2 / (xw1 * yw2 - xw2 * yw1)
      itab_w(iw)%gyps2(j1) =  xw1 / (xw1 * yw2 - xw2 * yw1)

!----------------------------------------

      if (itab_w(iw)%dirv(j2) < 0.) then
         alpha = atan2(yw2,xw2)   ! VC(iv) direction counterclockwise from east

         itab_v(iv)%cosv(1) = cos(alpha)
         itab_v(iv)%sinv(1) = sin(alpha)

         itab_v(iv)%dxps(1) = xv
         itab_v(iv)%dyps(1) = yv
      else
         alpha = atan2(-yw2,-xw2) ! VC(iv) direction counterclockwise from east

         itab_v(iv)%cosv(2) = cos(alpha)
         itab_v(iv)%sinv(2) = sin(alpha)

         itab_v(iv)%dxps(2) = xv
         itab_v(iv)%dyps(2) = yv
      endif

      endif

      itab_w(iw)%unx_w = -sin(glonw(iw))
      itab_w(iw)%uny_w = cos(glonw(iw))

      itab_w(iw)%vnx_w = -sin(glatw(iw)) * cos(glonw(iw))
      itab_w(iw)%vny_w = -sin(glatw(iw)) * sin(glonw(iw))
      itab_w(iw)%vnz_w = cos(glatw(iw))

!----------------------------------------
! END NEW SECTION JULY 2011
!----------------------------------------

   enddo

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

   !$omp do private(npoly, fo, a, b, work, info, j, iv, vdotw, vmag, fact, &
   !$omp            vnx_ps, vny_ps, vnz_ps, vrot_x, vrot_y)
   do iw = 2, nwa
      
      npoly = itab_w(iw)%npoly

      ! Default coefficients from Perot
      fo(1:npoly) = 2.0 * itab_w(iw)%farv(1:npoly)

      if (allocated(a)) then
         if (size(a,2) /= npoly) deallocate(a)
      endif

      if (.not. allocated(a)) allocate(a(3,npoly))

      if (mdomain < 2) then

         do j = 1, npoly
            iv = itab_w(iw)%iv(j)

            ! Compute the components of the V unit normals perpendicular to W 

            vdotw = vnx(iv)*wnx(iw) + vny(iv)*wny(iw) + vnz(iv)*wnz(iw)

            vnx_ps(j) = vnx(iv) - vdotw * wnx(iw)
            vny_ps(j) = vny(iv) - vdotw * wny(iw)
            vnz_ps(j) = vnz(iv) - vdotw * wnz(iw)

            ! Normalize these new vectors to unit length

            vmag = sqrt( vnx_ps(j)**2 + vny_ps(j)**2 + vnz_ps(j)**2 )

            vnx_ps(j) = vnx_ps(j) / vmag
            vny_ps(j) = vny_ps(j) / vmag
            vnz_ps(j) = vnz_ps(j) / vmag

            ! Rotate these new unit normals to a coordinate system with Z aligned with W

            if (wnz(iw) >= 0.0) then
            
               fact = ( wny(iw)*vnx_ps(j) - wnx(iw)*vny_ps(j) ) / ( 1.0 + wnz(iw) )

               vrot_x(j) = vnx_ps(j)*wnz(iw) - vnz_ps(j)*wnx(iw) + wny(iw)*fact
               vrot_y(j) = vny_ps(j)*wnz(iw) - vnz_ps(j)*wny(iw) - wnx(iw)*fact

            else

               fact = ( wny(iw)*vnx_ps(j) - wnx(iw)*vny_ps(j) ) / ( 1.0 - wnz(iw) )

               vrot_x(j) = -vnx_ps(j)*wnz(iw) + vnz_ps(j)*wnx(iw) + wny(iw)*fact
               vrot_y(j) = -vny_ps(j)*wnz(iw) + vnz_ps(j)*wny(iw) - wnx(iw)*fact

            endif

         enddo

      else

         do j = 1, npoly
            iv = itab_w(iw)%iv(j)

            vnx_ps(j) = vnx(iv)
            vny_ps(j) = vny(iv)
            vnz_ps(j) = 0.0

            vrot_x(j) = vnx_ps(j)
            vrot_y(j) = vny_ps(j)
         enddo

      endif

      a(1,1:npoly) = vrot_x(1:npoly) * vrot_x(1:npoly)
      a(2,1:npoly) = vrot_y(1:npoly) * vrot_y(1:npoly)
      a(3,1:npoly) = vrot_x(1:npoly) * vrot_y(1:npoly)

      b(1) = 1.0 - sum( fo(1:npoly) * a(1,:) )
      b(2) = 1.0 - sum( fo(1:npoly) * a(2,:) )
      b(3) =     - sum( fo(1:npoly) * a(3,:) )

      call sgels( 'N', 3, npoly, 1, a, 3, b, 7, work, lwork, info )

      ! Vector b is now the correction to the coefficients fo
      b(1:npoly) = b(1:npoly) + fo(1:npoly)

      if (info == 0 .and. all(b(1:npoly) > 0.05) .and. all(b(1:npoly) < 0.7)) then

         itab_w(iw)%ecvec_vx(1:npoly) = b(1:npoly) * vnx_ps(1:npoly)
         itab_w(iw)%ecvec_vy(1:npoly) = b(1:npoly) * vny_ps(1:npoly)
         itab_w(iw)%ecvec_vz(1:npoly) = b(1:npoly) * vnz_ps(1:npoly)

      else
      
         write(*,*) "Problem optimizing vector coefficients for iw = ", iw
         write(*,*) "Using default coefficients."

         itab_w(iw)%ecvec_vx(1:npoly) = fo(1:npoly) * vnx_ps(1:npoly)
         itab_w(iw)%ecvec_vy(1:npoly) = fo(1:npoly) * vny_ps(1:npoly)
         itab_w(iw)%ecvec_vz(1:npoly) = fo(1:npoly) * vnz_ps(1:npoly)

      endif

   enddo
   !omp end do
   
   if (allocated(a)) deallocate(a)

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

   if (myrank /= 0) return

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
                         jtm_grid, jtv_grid, jtv_wall, jtv_lbcp, &
                         jtw_grid, jtw_lbcp
  use misc_coms,   only: io6, mdomain, itopoflg
  use consts_coms, only: erad, r8
  use mem_grid,    only: nsw_max, nza, nma, nva, nwa, lpm, lpv, lpw, lsw, &
                         topm, topw, zm, dzt, zt, zfacm, zfact, dnu, dniu, dnv, &
                         arm0, arw0, arv, arw, volt, &
                         xem, yem, zem, xew, yew, zew, glatw, glonw, glatm, glonm
  use leaf_coms,   only: isfcl,nwl
  use sea_coms,    only: nws
  use mem_leaf,    only: land, itab_wl
  use mem_sea,     only: sea, itab_ws
  use land_db,     only: land_database_read
  use olam_mpi_atm,only: olam_stop

  implicit none

  integer :: j,iw,iwp,iv,ivp,im1,im2,k,km,im,iw1,iw2,iw3
  integer :: iws,iwl,kw,ks,npoly,jv,iv1,iv2,iv3
  real    :: hmin,hmax,topc,facw
  real(r8):: area
  integer :: status

  real(r8), allocatable :: area_sum(:,:)

! Loop over all land and sea cells, and sum information from them to compute
! ARW and VOLT

  do iwl = 2,nwl
     iw = itab_wl(iwl)%iw
     kw = itab_wl(iwl)%kw

     area = land%area(iwl)
     topc = land%topw(iwl)

     volt(kw,iw) = volt(kw,iw) + area * (zm(kw) - topc)

     do k = kw+1,nza
        arw(k-1,iw) = arw(k-1,iw) + area
         volt(k,iw) =  volt(k,iw) + area * (zm(k) - zm(k-1))
     enddo
  enddo

  do iws = 2,nws
     iw = itab_ws(iws)%iw
     kw = itab_ws(iws)%kw

     area = sea%area(iws)
     topc = sea%topw(iws)

     volt(kw,iw) = volt(kw,iw) + area * (zm(kw) - topc)

     do k = kw+1,nza
        arw(k-1,iw) = arw(k-1,iw) + area
         volt(k,iw) =  volt(k,iw) + area * (zm(k) - zm(k-1))
     enddo
  enddo

! Lateral boundary copy of ARW, VOLT

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

     arw (:,iw) = arw (:,iwp)
     volt(:,iw) = volt(:,iwp)
  enddo

! ARV

  write(io6,*) 'Defining control volume areas'

!----------------------------------------------------------------------
  !$omp parallel do private(iv,im1,im2,iw1,iw2,k,hmin,hmax,km)
  do j = 1,jtab_v(jtv_grid)%jend(1); iv = jtab_v(jtv_grid)%iv(j)
     im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

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

           if (volt(k,iw1) <= 1.e-9 .or. volt(k,iw2) <= 1.e-9) then

              ! close V if either T neighbor is completely closed
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

        enddo ! k

     endif

  enddo ! j,iv
  !$omp end parallel do

! Set ARV to zero at non-topo walls
! [ONLY FOR ISFCL = 0]

!----------------------------------------------------------------------
  do j = 1,jtab_v(jtv_wall)%jend(1); iv = jtab_v(jtv_wall)%iv(j)
!----------------------------------------------------------------------

     do k = 2,nza
        arv(k,iv) = 0.
     enddo
    lpv(iv) = nza

  enddo

! Lateral boundary copy of ARV

!----------------------------------------------------------------------
  do j = 1,jtab_v(jtv_lbcp)%jend(1); iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------

     arv(:,iv) = arv(:,ivp)

  enddo

!----------------------------------------------------------------------
  !$omp parallel do private(iw,npoly,k,facw,jv,iv)
  do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------

! Number of vertices of IW polygon

     npoly = itab_w(iw)%npoly

! Loop over vertical levels

     do k = 2,nza

! Close top area of T cell if volume is zero

        if (volt(k,iw) < 1.e-9) then

!          print*, 'closing arw(k,iw) since volt(k,iw) is zero ',k,iw

           arw (k,iw) = 0.
           volt(k,iw) = 1.e-9
        endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Option for stability: expand volt if too small relative to any grid cell face

        volt(k,iw) = max(volt(k,iw), &
                     0.5_r8 * real(dzt(k) * (arw(k,iw) + arw(k-1,iw)), r8))

! Loop over faces of IW polygon

        facw = 0.0

        do jv = 1,npoly
           iv = itab_w(iw)%iv(jv)

           ! Each V face contributes to open up its fraction farv
           ! of the current polygon:
           facw = facw + itab_w(iw)%farv(jv) * arv(k,iv) / (dnu(iv) * dzt(k))

        enddo

        facw = max( min(facw,1.0), 0.0)
        volt(k,iw) = max( volt(k,iw), real(arw0(iw) * dzt(k) * facw, r8) )

! Reset arw(k,iw) if we increased volt to avoid hollow box

        if (k < nza .and. volt(k,iw) > 1.e-9) then
           arw(k,iw) = max( arw(k,iw), real( volt(k,iw) / dzt(k)) )
        endif

! Make sure arw increases with height if we have adjusted it

        if (k < nza-1) then
           arw(k+1,iw) = max( arw(k+1,iw), arw(k,iw) )
        endif

! End option for stability
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     enddo  ! k

! Set arw = 0 for bottom (k = 1) and wall-on-top (k = nza) levels

     arw(1,iw) = 0.
     arw(nza,iw) = 0. 

  enddo
  !$omp end parallel do

! Lateral boundary copy of ARW, VOLT

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

     arw (:,iw) = arw (:,iwp)
     volt(:,iw) = volt(:,iwp)

  enddo

! LPM, LPW, LPV

  lpm(2:nma) = nza
  lpv(2:nva) = nza
  lpw(2:nwa) = nza

!----------------------------------------------------------------------
  do j = 1,jtab_m(jtm_grid)%jend(1); im = jtab_m(jtm_grid)%im(j)
!----------------------------------------------------------------------

     iv1 = itab_m(im)%iv(1)
     iv2 = itab_m(im)%iv(2)
     iv3 = itab_m(im)%iv(3)

     iw1 = itab_m(im)%iw(1)
     iw2 = itab_m(im)%iw(2)
     iw3 = itab_m(im)%iw(3)

! Loop over vertical levels from top to bottom

     do k = nza,2,-1

        if (topm(im) >= zm(k)) exit

        if (arv(k,iv1) < .005 * dnu(iv1) * dzt(k)) exit
        if (arv(k,iv2) < .005 * dnu(iv2) * dzt(k)) exit
        if (arv(k,iv2) < .005 * dnu(iv2) * dzt(k)) exit

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

! Close all T cells that are below LPW

  nsw_max = 1

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------

! Number of vertices of IW polygon

     npoly = itab_w(iw)%npoly

     do k = 2,nza

        if (k < lpw(iw)) then
           arw (k,iw) = 0.
           volt(k,iw) = 1.e-9
        endif

     enddo  ! k

! Set arw = 0 for bottom (k = 1) and wall-on-top (k = nza) levels

     arw(1,iw) = 0.   
     arw(nza,iw) = 0.   

! Initialize LSW(IW)

     lsw(iw) = 1

! Loop over vertical levels from top to bottom

     do k = nza-1, lpw(iw), -1

! Increase LSW if K-1 W level intersects topography in this cell

        if (arw(k,iw) < .999 * arw0(iw)) then
           lsw(iw) = lsw(iw) + 1
        else
           arw(k,iw) = arw0(iw)
        endif

    enddo  ! k

! Increase nsw_max if necessary

   if (lsw(iw) > nsw_max) nsw_max = lsw(iw)

  enddo

! In case ARW has been reset to 0 anywhere (because it was nearly zero), 
! transfer the sea and land cell values to KW = LPW(IW).
! Also compute fractional sea and land cell areas per vertical level

  allocate(area_sum(nsw_max,nwa), stat=status)
  if (status /= 0) then
     call olam_stop( "Error allocating memory in hex_grid." )
  endif

  area_sum(:,:) = 0.0

  do iws = 2,nws
     iw = itab_ws(iws)%iw
     kw = itab_ws(iws)%kw

     if (kw < lpw(iw))               itab_ws(iws)%kw = lpw(iw)
     if (kw > lpw(iw) + lsw(iw) - 1) itab_ws(iws)%kw = lpw(iw) + lsw(iw) - 1

     ks = itab_ws(iws)%kw - lpw(iw) + 1
     area_sum(ks,iw) = area_sum(ks,iw) + sea%area(iws)
  enddo

  do iwl = 2,nwl
     iw = itab_wl(iwl)%iw
     kw = itab_wl(iwl)%kw

     if (kw < lpw(iw))               itab_wl(iwl)%kw = lpw(iw)
     if (kw > lpw(iw) + lsw(iw) - 1) itab_wl(iwl)%kw = lpw(iw) + lsw(iw) - 1

     ks = itab_wl(iwl)%kw - lpw(iw) + 1
     area_sum(ks,iw) = area_sum(ks,iw) + land%area(iwl)
  enddo

  do iws = 2,nws
     iw = itab_ws(iws)%iw
     kw = itab_ws(iws)%kw
     ks = kw - lpw(iw) + 1
     sea%area(iws) = sea%area(iws) * (arw(kw,iw) - arw(kw-1,iw)) / area_sum(ks,iw)
  enddo

  do iwl = 2,nwl
     iw = itab_wl(iwl)%iw
     kw = itab_wl(iwl)%kw
     ks = kw - lpw(iw) + 1
     land%area(iwl) = land%area(iwl) * (arw(kw,iw) - arw(kw-1,iw)) / area_sum(ks,iw)
  enddo

  deallocate(area_sum)

! Compute sea and land area fractions and expand surface areas with height
! for spherical geometry

  do iws = 2, nws
     kw = itab_ws(iws)%kw
     iw = itab_ws(iws)%iw
     itab_ws(iws)%arf_iw = sea%area(iws) / arw0(iw)
     itab_ws(iws)%arf_kw = sea%area(iws) / (arw(kw,iw) - arw(kw-1,iw))
     if (mdomain < 2) sea%area(iws) = sea%area(iws) * zfacm(kw-1)**2
  enddo

  do iwl = 2, nwl
     kw = itab_wl(iwl)%kw
     iw = itab_wl(iwl)%iw
     itab_wl(iwl)%arf_iw = land%area(iwl) / arw0(iw)
     itab_wl(iwl)%arf_kw = land%area(iwl) / (arw(kw,iw) - arw(kw-1,iw))
     if (mdomain < 2) land%area(iwl) = land%area(iwl) * zfacm(kw-1)**2
  enddo

! Expand ARW and VOLT with height for spherical geometry

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------

! Loop over vertical levels

     do k = 2,nza

        if (volt(k,iw) > 1.e-9) then
           if (mdomain < 2) then
              arw (k,iw) = arw (k,iw) * zfacm(k)**2
              volt(k,iw) = volt(k,iw) * zfact(k)**2
           endif
        else
           arw (k,iw) = 0.
           volt(k,iw) = 1.e-9
        endif

     enddo  ! k

! Set arw = 0 for bottom (k = 1) and wall-on-top (k = nza) levels

     arw(1,iw) = 0.   
     arw(nza,iw) = 0.   

  enddo

! Lateral boundary copy of ARW, VOLT, LPW, and LSW

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
     iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

     arw  (:,iw) = arw  (:,iwp)
     volt (:,iw) = volt (:,iwp)
     lpw    (iw) = lpw    (iwp)
     lsw    (iw) = lsw    (iwp)

  enddo

! Expand ARV with height for spherical geometry

!----------------------------------------------------------------------
  do j = 1,jtab_v(jtv_grid)%jend(1); iv = jtab_v(jtv_grid)%iv(j)
!----------------------------------------------------------------------

     if (mdomain < 2) then
        do k = 2,nza
           arv(k,iv) = arv(k,iv) * zfact(k)
        enddo
   endif

  enddo

! Lateral boundary copy of ARV

!----------------------------------------------------------------------
  do j = 1,jtab_v(jtv_lbcp)%jend(1); iv = jtab_v(jtv_lbcp)%iv(j)
     ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------

     arv(:,iv) = arv(:,ivp)
     lpv(iv) = lpv(ivp)

  enddo

end subroutine ctrlvols_hex
