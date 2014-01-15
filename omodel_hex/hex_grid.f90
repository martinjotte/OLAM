!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine voronoi()

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd, itab_m, itab_v, itab_w,  &
                       mrls, alloc_itabs, mloops

use mem_grid,    only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa,  &
                       xem, yem, zem, xew, yew, zew,  &
                       alloc_xyzem, alloc_xyzew
use misc_coms,   only: io6, mdomain, meshtype
use consts_coms, only: pi2, erad, erador5

implicit none

integer :: im1,im2,im3
integer :: iw1,iw2,iw3,im,iw
integer :: nmad,nuad,nwad
integer :: iwd,iv,iud,iud1,iud2,imd,imd1,imd2,npoly,j,j1

real :: expansion

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

call alloc_itabs(meshtype,nma,nua,nva,nwa)

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

   itab_v(iv)%loop(:)  = itab_ud(iud)%loop(:)

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

!!! IV(13) and IV(15) neighbors of IV
!!
!!      if (npoly >= 6 .and. iud2 == itab_v(iv)%iv(5)) then
!!         itab_v(iv)%iv(13) = iud1
!!      endif
!!
!!      if (npoly == 7 .and. iud1 == itab_v(iv)%iv(9)) then
!!         itab_v(iv)%iv(15) = iud2
!!      endif

   enddo           

!!! Extract information from IMD2 neighbor
!!
!!   imd = itab_ud(iud)%im(2)
!!   npoly = itab_md(imd)%npoly
!!
!!   do j = 1,npoly
!!      j1 = j + 1
!!      if (j == npoly) j1 = 1
!!
!!      iud1 = itab_md(imd)%iu(j)
!!      iud2 = itab_md(imd)%iu(j1)
!!
!!! IV(14) and IV(16) neighbors of IV
!!
!!      if (npoly >= 6 .and. iud1 == itab_v(iv)%iv(8)) then
!!         itab_v(iv)%iv(14) = iud2
!!      endif
!!
!!      if (npoly == 7 .and. iud2 == itab_v(iv)%iv(12)) then
!!         itab_v(iv)%iv(16) = iud1
!!      endif
!!
!!   enddo           

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

   itab_w(iw)%mrlw = itab_md(imd)%mrlm
   itab_w(iw)%mrlw_orig = itab_md(imd)%mrlm_orig

   npoly = itab_w(iw)%npoly

! Loop over IM/IV neighbors of IW

   do j = 1,itab_w(iw)%npoly
      im = itab_md(imd)%iw(j)
      iwd = im
      iv = itab_md(imd)%iu(j)

      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

!! Set mrlw to max of itab_wd neighbors
!
!      if (itab_w(iw)%mrlw < itab_wd(iwd)%mrlw) then
!          itab_w(iw)%mrlw = itab_wd(iwd)%mrlw
!      endif

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
   itab_m(im)%mrlm_orig = itab_wd(iwd)%mrlw_orig
   itab_m(im)%mrow      = itab_wd(iwd)%mrow
   itab_m(im)%mrowh     = itab_wd(iwd)%mrowh

   itab_m(im)%iv(1:3)   = itab_wd(iwd)%iu(1:3)
   itab_m(im)%iw(1:3)   = itab_wd(iwd)%im(1:3)

! Loop over IW neighbors of IM

!   do j = 1,itab_m(im)%npoly
!      iw = itab_m(im)%iw(j)
!
!! Set mrlm to max of itab_w neighbors
!
!      if (itab_m(im)%mrlm < itab_w(iw)%mrlw) then
!          itab_m(im)%mrlm = itab_w(iw)%mrlw
!      endif
!   enddo
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
real :: onmx,onmy,onmz

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

subroutine grid_geometry_hex(quarter_kite)

use mem_ijtabs,  only: itab_m, itab_v, itab_w
use mem_grid,    only: nza, nma, nva, nwa, xev, yev, zev, xem, yem, zem, &
                       xew, yew, zew, unx, uny, unz, wnx, wny, wnz,      &
                       vnx, vny, vnz, glonw, glatw, dnu, dniu, dnv, dniv, arw0, &
                       arm0, glonm, glatm, glatv, glonv
use misc_coms,   only: io6, mdomain, grdlat, grdlon, nxp, rinit
use consts_coms, only: erad, erad2, piu180, eradsq,pio2
use oplot_coms,  only: op

implicit none

real, intent(out) :: quarter_kite(2,nva)

integer :: im,iv,iw
integer :: im1,im2,im3,im4,im5,im6
integer :: iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8,iv9,iv10
integer :: iv11,iv12,iv13,iv14,iv15,iv16
integer :: iw1,iw2,iw3
integer :: itpn,npoly
integer :: j,jj,jmax,jmaxneg,jminpos,jn
integer :: ivn,iwn

real :: expansion
real :: raxis
real :: hper

real :: xiw,yiw,ziw
real :: xij(6),yij(6),zij(6)
real :: x1,x2,x3,y1,y2,y3
real :: scalprod,scalprod_max
real :: vecjx,vecjy,vecjz
real :: vecmx,vecmy,vecmz
real :: vecmjx,vecmjy,vecmjz
real :: xi,yi,xj(6),yj(6)
real :: vecprodz,vecprodz_maxneg,vecprodz_minpos

real :: ef,x12,y12,z12,x34,y34,z34
real :: b11,b21,b31,b12,b22,b32,b13,b23,b33 &
       ,c11,c21,c31,c12,c22,c32,c13,c23,c33 &
       ,d11,d21,d31,d12,d22,d32,d13,d23,d33
real :: s1, s2, s3

real :: dwm1,dwm2,dwv,dvm1,dvm2,swm1,swm2,swv,svm1,svm2,angw1,angw2,angm1,angm2

integer :: jm,jv,j1,j2

integer :: iv_w
integer :: npoly1,npoly2,npolym,jw_m,iw_m,np

real :: xm1,xm2,xv,xw,ym1,ym2,yv,yw,frac,accum_kite,alpha

real :: xw1,xw2,yw1,yw2
real :: gx1,gx2,gy1,gy2

real :: xem1, xem2, yem1, yem2, zem1, zem2

real :: xq1, yq1, xq2, yq2, psiz, vsprd
integer :: iskip, iwp, ivp
logical :: dops

character(10) :: string

integer, parameter :: lwork = 200
real               :: work(lwork)
integer            :: info
real               :: vdotw, vmag, fact
real               :: b(7), fo(7), vnx_ps(7), vny_ps(7), vnz_ps(7), vrot_x(7), vrot_y(7)
real, allocatable  :: a(:,:)

ef = 1.01  ! radial expansion factor (from earth center) for defining 
           ! auxiliary point for computing unit normal to U face

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

! Number of polygon edges/vertices

!   npoly = itab_m(im)%npoly
!
! Loop over all polygon edges
!
!   do j = 1,npoly
!      iw = itab_m(im)%iw(j)
!
!      itab_m(im)%mrlm = max(itab_m(im)%mrlm, itab_w(iw)%mrlw)
!   enddo

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

! fww is inner product of w and wn unit vectors

      itab_w(iw)%fww(j2) = wnx(iw) * wnx(iw2) &
                         + wny(iw) * wny(iw2) &
                         + wnz(iw) * wnz(iw2)

! The following needed in case of roundoff error

      if (itab_w(iw)%fww(j2) > 1.) itab_w(iw)%fww(j2) = 1.

! fwv from fww using pythagorean theorem, then halved for averaging 2 V values

      itab_w(iw)%fwv(j2) = .5 * sqrt(1. - (itab_w(iw)%fww(j2))**2) &
                         * itab_w(iw)%dirv(j2)

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

! FARW(1) and FARW(2) interpolation coefficients for ARW and VOLV
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

!---------------------------------------------------
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

! End plot grid lines
!---------------------------------------------------

!! THE FOLLOWING COEFFICIENTS ARE UNUSED WITH PEROT'S METHOD FOR HEXAGONS
!!
!!! Loop over all V points
!!
!!do iv = 2,nva
!!
!!! Indices of neighboring M-points
!!
!!   im1 = itab_v(iv)%im(1)
!!   im2 = itab_v(iv)%im(2)
!!   im3 = itab_v(iv)%im(3)
!!   im4 = itab_v(iv)%im(4)
!!   im5 = itab_v(iv)%im(5)
!!   im6 = itab_v(iv)%im(6)
!!
!!! V neighbor indices
!!
!!   iv1  = itab_v(iv)%iv(1)
!!   iv2  = itab_v(iv)%iv(2)
!!   iv3  = itab_v(iv)%iv(3)
!!   iv4  = itab_v(iv)%iv(4)
!!   iv5  = itab_v(iv)%iv(5)
!!   iv6  = itab_v(iv)%iv(6)
!!   iv7  = itab_v(iv)%iv(7)
!!   iv8  = itab_v(iv)%iv(8)
!!   iv9  = itab_v(iv)%iv(9)
!!   iv10 = itab_v(iv)%iv(10)
!!   iv11 = itab_v(iv)%iv(11)
!!   iv12 = itab_v(iv)%iv(12)
!!   iv13 = itab_v(iv)%iv(13)
!!   iv14 = itab_v(iv)%iv(14)
!!   iv15 = itab_v(iv)%iv(15)
!!   iv16 = itab_v(iv)%iv(16)
!!
!!! OR:  ivv(1:16) = itab_v(iv)%iv(1:16)   
!!
!!   iw1 = itab_v(iv)%iw(1)
!!   iw2 = itab_v(iv)%iw(2)
!!
!!! Number of V neighbors of IW1 and IW2
!!
!!   npoly1 = itab_w(iw1)%npoly
!!   npoly2 = itab_w(iw2)%npoly
!!
!!! FUV interpolation coefficients
!!
!!!---------------------------------------------------------------------
!!! Progressing clockwise from IV in IW1
!!
!!   if (itab_v(iv1)%iw(1) == iw1) then
!!      accum_kite = quarter_kite(1,iv) + quarter_kite(2,iv1)
!!      itab_v(iv)%fuv(1) = -(.5 - accum_kite / arw0(iw1))
!!   else
!!      accum_kite = quarter_kite(1,iv) + quarter_kite(1,iv1)
!!      itab_v(iv)%fuv(1) = (.5 - accum_kite / arw0(iw1))
!!   endif
!!
!!   if (itab_v(iv5)%iw(1) == iw1) then
!!      accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv1) &
!!                 + quarter_kite(2,iv1) + quarter_kite(2,iv5)
!!      itab_v(iv)%fuv(5) = -(.5 - accum_kite / arw0(iw1))
!!   else
!!      accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv1) &
!!                 + quarter_kite(2,iv1) + quarter_kite(1,iv5)
!!      itab_v(iv)%fuv(5) =  (.5 - accum_kite / arw0(iw1))
!!   endif
!!
!!   if (iv13 > 1) then
!!
!!      if (itab_v(iv13)%iw(1) == iw1) then
!!         accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv1) &
!!                    + quarter_kite(2,iv1) + quarter_kite(1,iv5) &
!!                    + quarter_kite(2,iv5) + quarter_kite(2,iv13)
!!         itab_v(iv)%fuv(13) = -(.5 - accum_kite / arw0(iw1))
!!      else
!!         accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv1) &
!!                    + quarter_kite(2,iv1) + quarter_kite(1,iv5) &
!!                    + quarter_kite(2,iv5) + quarter_kite(1,iv13)
!!         itab_v(iv)%fuv(13) =  (.5 - accum_kite / arw0(iw1))
!!      endif
!!
!!   endif
!!
!!!---------------------------------------------------------------------
!!! Progressing counterclockwise from IV in IW1
!!
!!   if (itab_v(iv3)%iw(1) == iw1) then
!!      accum_kite = quarter_kite(2,iv) + quarter_kite(1,iv3)
!!      itab_v(iv)%fuv(3) = (.5 - accum_kite / arw0(iw1))
!!   else
!!      accum_kite = quarter_kite(2,iv) + quarter_kite(2,iv3)
!!      itab_v(iv)%fuv(3) = -(.5 - accum_kite / arw0(iw1))
!!   endif
!!
!!   if (itab_v(iv9)%iw(1) == iw1) then
!!      accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv3) &
!!                 + quarter_kite(2,iv3) + quarter_kite(1,iv9)
!!      itab_v(iv)%fuv(9) = (.5 - accum_kite / arw0(iw1))
!!   else
!!      accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv3) &
!!                 + quarter_kite(2,iv3) + quarter_kite(2,iv9)
!!      itab_v(iv)%fuv(9) =  -(.5 - accum_kite / arw0(iw1))
!!   endif
!!
!!   if (iv15 > 1) then
!!   
!!      if (itab_v(iv15)%iw(1) == iw1) then
!!         accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv3) &
!!                    + quarter_kite(2,iv3) + quarter_kite(1,iv9) &
!!                    + quarter_kite(2,iv9) + quarter_kite(1,iv15)
!!         itab_v(iv)%fuv(15) = (.5 - accum_kite / arw0(iw1))
!!      else
!!         accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv3) &
!!                    + quarter_kite(2,iv3) + quarter_kite(1,iv9) &
!!                    + quarter_kite(2,iv9) + quarter_kite(2,iv15)
!!         itab_v(iv)%fuv(15) =  -(.5 - accum_kite / arw0(iw1))
!!      endif
!!
!!   endif
!!
!!!---------------------------------------------------------------------
!!! Progressing counterclockwise from IV in IW2
!!
!!   if (itab_v(iv2)%iw(1) == iw2) then
!!      accum_kite = quarter_kite(1,iv) + quarter_kite(1,iv2)
!!      itab_v(iv)%fuv(2) = -(.5 - accum_kite / arw0(iw2))
!!   else
!!      accum_kite = quarter_kite(1,iv) + quarter_kite(2,iv2)
!!      itab_v(iv)%fuv(2) = (.5 - accum_kite / arw0(iw2))
!!   endif
!!
!!   if (itab_v(iv8)%iw(1) == iw2) then
!!      accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv2) &
!!                 + quarter_kite(2,iv2) + quarter_kite(1,iv8)
!!      itab_v(iv)%fuv(8) = -(.5 - accum_kite / arw0(iw2))
!!   else
!!      accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv2) &
!!                 + quarter_kite(2,iv2) + quarter_kite(2,iv8)
!!      itab_v(iv)%fuv(8) =  (.5 - accum_kite / arw0(iw2))
!!   endif
!!
!!   if (iv14 > 1) then
!!
!!      if (itab_v(iv14)%iw(1) == iw2) then
!!         accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv2) &
!!                    + quarter_kite(2,iv2) + quarter_kite(1,iv8) &
!!                    + quarter_kite(2,iv8) + quarter_kite(1,iv14)
!!         itab_v(iv)%fuv(14) = -(.5 - accum_kite / arw0(iw2))
!!      else
!!         accum_kite = quarter_kite(1,iv)  + quarter_kite(1,iv2) &
!!                    + quarter_kite(2,iv2) + quarter_kite(1,iv8) &
!!                    + quarter_kite(2,iv8) + quarter_kite(2,iv14)
!!         itab_v(iv)%fuv(14) =  (.5 - accum_kite / arw0(iw2))
!!      endif
!!
!!   endif
!!
!!!---------------------------------------------------------------------
!!! Progressing clockwise from IV in IW2
!!
!!   if (itab_v(iv4)%iw(1) == iw2) then
!!      accum_kite = quarter_kite(2,iv) + quarter_kite(2,iv4)
!!      itab_v(iv)%fuv(4) = (.5 - accum_kite / arw0(iw2))
!!   else
!!      accum_kite = quarter_kite(2,iv) + quarter_kite(1,iv4)
!!      itab_v(iv)%fuv(4) = -(.5 - accum_kite / arw0(iw2))
!!   endif
!!
!!   if (itab_v(iv12)%iw(1) == iw2) then
!!      accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv4) &
!!                 + quarter_kite(2,iv4) + quarter_kite(2,iv12)
!!      itab_v(iv)%fuv(12) = (.5 - accum_kite / arw0(iw2))
!!   else
!!      accum_kite = quarter_kite(2,iv)  + quarter_kite(1,iv4) &
!!                 + quarter_kite(2,iv4) + quarter_kite(1,iv12)
!!      itab_v(iv)%fuv(12) = -(.5 - accum_kite / arw0(iw2))
!!   endif
!!
!!   if (iv16 > 1) then
!!
!!      if (itab_v(iv16)%iw(1) == iw2) then
!!         accum_kite = quarter_kite(2,iv)   + quarter_kite(1,iv4)  &
!!                    + quarter_kite(2,iv4)  + quarter_kite(1,iv12) &
!!                    + quarter_kite(2,iv12) + quarter_kite(2,iv16)
!!         itab_v(iv)%fuv(16) = (.5 - accum_kite / arw0(iw2))
!!      else
!!         accum_kite = quarter_kite(2,iv)   + quarter_kite(1,iv4)  &
!!                    + quarter_kite(2,iv4)  + quarter_kite(1,iv12) &
!!                    + quarter_kite(2,iv12) + quarter_kite(1,iv16)
!!         itab_v(iv)%fuv(16) = -(.5 - accum_kite / arw0(iw2))
!!      endif
!!
!!   endif
!!
!!!---------------------------------------------------------------------
!!
!!! Skip this V point if iw1 < 2 or iw2 < 2
!!
!!   if (iw1 < 2 .or. iw2 < 2) cycle
!!
!!! Project neighbor V & W unit vectors onto V unit vector
!!
!!   call matrix_3x3(vnx(iv5),vnx(iv6),xem(im3)/erad  &
!!                  ,vny(iv5),vny(iv6),yem(im3)/erad  &
!!                  ,vnz(iv5),vnz(iv6),zem(im3)/erad  &
!!                  ,vnx(iv),vny(iv),vnz(iv)     &
!!                  ,itab_v(iv)%fvv(5)           &
!!                  ,itab_v(iv)%fvv(6)           &
!!                  ,itab_v(iv)%fvw(1)           )
!!
!!   call matrix_3x3(vnx(iv7),vnx(iv8),xem(im4)/erad  &
!!                  ,vny(iv7),vny(iv8),yem(im4)/erad  &
!!                  ,vnz(iv7),vnz(iv8),zem(im4)/erad  &
!!                  ,vnx(iv),vny(iv),vnz(iv)     &
!!                  ,itab_v(iv)%fvv(7)           &
!!                  ,itab_v(iv)%fvv(8)           &
!!                  ,itab_v(iv)%fvw(2)           )
!!
!!   call matrix_3x3(vnx(iv9),vnx(iv10),xem(im5)/erad  &
!!                  ,vny(iv9),vny(iv10),yem(im5)/erad  &
!!                  ,vnz(iv9),vnz(iv10),zem(im5)/erad  &
!!                  ,vnx(iv),vny(iv),vnz(iv)      &
!!                  ,itab_v(iv)%fvv(9)            &
!!                  ,itab_v(iv)%fvv(10)           &
!!                  ,itab_v(iv)%fvw(3)            )
!!
!!   call matrix_3x3(vnx(iv11),vnx(iv12),xem(im6)/erad  &
!!                  ,vny(iv11),vny(iv12),yem(im6)/erad  &
!!                  ,vnz(iv11),vnz(iv12),zem(im6)/erad  &
!!                  ,vnx(iv),vny(iv),vnz(iv)       &
!!                  ,itab_v(iv)%fvv(11)            &
!!                  ,itab_v(iv)%fvv(12)            &
!!                  ,itab_v(iv)%fvw(4)             )
!!
!!! Divide fvw1-4 by 2 for use with two W points
!!
!!   itab_v(iv)%fvw(1:4) = .5 * itab_v(iv)%fvw(1:4)
!!
!!enddo  ! IV
!!
!!! Loop over all M points
!!
!!do im = 2,nma
!!
!!   iv1 = itab_m(im)%iv(1)
!!   iv2 = itab_m(im)%iv(2)
!!   iv3 = itab_m(im)%iv(3)
!!
!!   iw1 = itab_m(im)%iw(1)
!!   iw2 = itab_m(im)%iw(2)
!!   iw3 = itab_m(im)%iw(3)
!!
!!! ITAB_M(IM)%FMW coefficient for interpolating from IW to IM [not currently used]
!!
!!   if (itab_v(iv1)%im(1) == im) then
!!      itab_m(im)%fmw(2) = itab_m(im)%fmw(2) + quarter_kite(1,iv1) / arm0(im)
!!      itab_m(im)%fmw(3) = itab_m(im)%fmw(3) + quarter_kite(1,iv1) / arm0(im)
!!   else
!!      itab_m(im)%fmw(2) = itab_m(im)%fmw(2) + quarter_kite(2,iv1) / arm0(im)
!!      itab_m(im)%fmw(3) = itab_m(im)%fmw(3) + quarter_kite(2,iv1) / arm0(im)
!!   endif
!!
!!   if (itab_v(iv2)%im(1) == im) then
!!      itab_m(im)%fmw(3) = itab_m(im)%fmw(3) + quarter_kite(1,iv2) / arm0(im)
!!      itab_m(im)%fmw(1) = itab_m(im)%fmw(1) + quarter_kite(1,iv2) / arm0(im)
!!   else
!!      itab_m(im)%fmw(3) = itab_m(im)%fmw(3) + quarter_kite(2,iv2) / arm0(im)
!!      itab_m(im)%fmw(1) = itab_m(im)%fmw(1) + quarter_kite(2,iv2) / arm0(im)
!!   endif
!!
!!   if (itab_v(iv3)%im(1) == im) then
!!      itab_m(im)%fmw(1) = itab_m(im)%fmw(1) + quarter_kite(1,iv3) / arm0(im)
!!      itab_m(im)%fmw(2) = itab_m(im)%fmw(2) + quarter_kite(1,iv3) / arm0(im)
!!   else
!!      itab_m(im)%fmw(1) = itab_m(im)%fmw(1) + quarter_kite(2,iv3) / arm0(im)
!!      itab_m(im)%fmw(2) = itab_m(im)%fmw(2) + quarter_kite(2,iv3) / arm0(im)
!!   endif
!!
!!enddo

end subroutine grid_geometry_hex

!===============================================================================

subroutine ctrlvols_hex(quarter_kite)

use mem_ijtabs,  only: jtab_m, jtab_v, jtab_w, itab_m, itab_v, itab_w, &
                       jtm_grid, jtv_grid, jtv_wall, jtv_lbcp, &
                       jtw_grid, jtw_lbcp
use misc_coms,   only: io6, mdomain, itopoflg
use consts_coms, only: erad
use mem_grid,    only: nsw_max, nza, nma, nva, nwa, lpm, lpv, lpw, lsw,  &
                       topm, topw, zm, dzt, zt, zfacm, zfact, dnu, dniu, dnv, &
                       arm0, arw0, arv, arw, volt, volti, volvi, volwi,  &
                       xem, yem, zem, xew, yew, zew, glatw, glonw, glatm, glonm
use leaf_coms,   only: isfcl,nml
use sea_coms,    only: nms
use mem_sflux,   only: init_fluxcells,nseaflux,nlandflux,seaflux,landflux
use mem_leaf,    only: land
use mem_sea,     only: sea

implicit none

real, intent(in) :: quarter_kite(2,nva)

integer :: j,jw,iw,iwp,iter,iv,ivp,im1,im2,k,km,im  &
   ,im11,im21,im12,im22,iu1,iu2,iw1,iw2,iw3,kp  &
   ,iu1a,iu1b,iu2a,iu2b  &
   ,iuo1a,iuo1b,iuo2a,iuo2b  &
   ,im3,iu3,ka,ks
integer :: isf,ilf,kw

integer, parameter :: npass = 2
integer :: ipass

real :: hmin,hmax,sum1,sumk,w1,w2,t1,t2,t3,hm,qk,dt13,dt32,t13,t32,farw1,farw2
real, allocatable :: area_kw_sum(:,:)

!!!!!!!!!!!!! special quadrature parameters

integer, parameter :: np = 30
integer :: ihalf,iwn,ip,jp,ipair,npoly,jm,jm1,jm2,jv,js1,js2,npolym
integer :: iv1,iv2,iv3,iv4,iv5,iv8,iv9,iv12,iv13,iv14,iv15,iv16

real :: top1,top2,top3,top4
real :: topp,topm13,topm32
real :: facj1,facj2,faci1,faci2
real :: del_arw,v_height

real :: fuv(16)
real(8) :: facw
real :: az,bz,cz,xq(3),yq(3),zq(3),fracz

!!!!!!!!!!!!! end special

! Check whether LEAF and SEA models are being used

if (isfcl == 0) then

! LEAF and SEA models are NOT being used.  Initialize topography without
! reading a surface file.

   if (itopoflg == 1) then

! Read TOPW from dataset

      call topo_database_read(nwa,xew,yew,zew,topw,arw0,glatw,glonw)

   else

! Initialize TOPW from default value (may customize here)

      call topo_init(nwa,topw,glatw,glonw,xew,yew,zew)  

   endif

   !$omp parallel

! Prevent TOPW lower than lowest model level zm(1)

   !$omp workshare
   topw(2:nwa) = max(topw(2:nwa),zm(1))
   !$omp end workshare

! Loop over all Hex W points

   !$omp do private(k,fracz)
   do iw = 2,nwa

! Loop over KM levels

      do k = 1,nza-2

! Check whether topw(iw) is at or above current KM level

         if (zm(k) <= topw(iw) .and. zm(k+1) >= topw(iw)) then

! Check topw(iw) within height range of KT level

            fracz = (topw(iw) - zm(k)) / (zm(k+1) - zm(k))

! If topw(iw) is in lowest or highest 1% of KT level, move it to the limit.      

            if (fracz < .01) then
               topw(iw) = zm(k)
            elseif (fracz > .99) then
               topw(iw) = zm(k+1)
            endif

            exit

         endif

      enddo

   enddo
   !$omp end do

! Fill TOPM by bi-linearly interpolating TOPW values

! Loop over all Hex M points

   !$omp do private(jw,iw,xq,yq,zq,az,bz,cz)
   do im = 2,nma

! Loop over all neighbor W points of this IM

      do jw = 1,3
         iw = itab_m(im)%iw(jw)

! Evaluate x,y coordinates of W point on polar stereographic plane
! tangent at IM

         if (mdomain <= 1) then
            call e_ps(xew(iw),yew(iw),zew(iw),glatm(im),glonm(im),xq(jw),yq(jw))
         else
            xq(jw) = xew(iw) - xem(im)
            yq(jw) = yew(iw) - yem(im)
         endif

! Store topography height of W point in zq array

         zq(jw) = topw(iw)

      enddo

! Find fit coefficients for linear elevation surface of triangle defined
! by 3 IW points

      call matrix_3x3(1.,xq(1),yq(1)    &
                     ,1.,xq(2),yq(2)    &
                     ,1.,xq(3),yq(3)    &
                     ,zq(1),zq(2),zq(3) &
                     ,az,bz,cz          )

! Find height of M point (located at x=0,y=0)

      topm(im) = az

   enddo
   !$omp end do

! Fill ARW and VOLT from topography using quadrature method

   call psub()
!----------------------------------------------------------------------
   !$omp do private(iw,npoly,t1,jv,jm1,jm2,im1,im2,iv,iwn,t2,ihalf,t3,    &
   !$omp            qk,k,dt13,dt32,del_arw,jp,facj1,facj2,ip,faci1,faci2, &   
   !$omp            top1,top2,top3,top4,ipair,topp)
   do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------
   call qsub('W',iw)

! Initialize arw and volt to zero

      arw (1:nza,iw) = 0.
      volt(1:nza,iw) = 0.

! Number of vertices of IW polygon

      npoly = itab_w(iw)%npoly

      t1 = topw(iw)

! Loop over polygon edges

      do jv = 1,npoly

! Polygon vertices adjacent to jv

         jm1 = jv - 1
         jm2 = jv
         if (jv == 1) jm1 = npoly

         im1 = itab_w(iw)%im(jm1)
         im2 = itab_w(iw)%im(jm2)

         iv = itab_w(iw)%iv(jv)

         iwn = itab_w(iw)%iw(jv)

! Topography at V point (t2)

         t2 = .5 * (topw(iw) + topw(iwn))

! Loop over both right-triangle halves of jv sector

         do ihalf = 1,2

! Quarter-kite and topography of 3rd point

            if (ihalf == 1) then
               t3 = topm(im1)

               if (itab_v(iv)%im(1) == im1) then
                  qk = quarter_kite(1,iv)
               else
                  qk = quarter_kite(2,iv)
               endif
            else
               t3 = topm(im2)

               if (itab_v(iv)%im(1) == im1) then
                  qk = quarter_kite(2,iv)
               else
                  qk = quarter_kite(1,iv)
               endif
            endif

! Loop over model levels

            do k = 2,nza-1

               if (t1 < zm(k-1) .and. t2 < zm(k-1) .and. t3 < zm(k-1)) then

! Bottom of this model layer is above all 3 topography points

                  arw (k,iw) = arw (k,iw) + qk
                  volt(k,iw) = volt(k,iw) + qk * dzt(k)

               else

! Interpolate 3 topography heights to np*np sub-triangles in this sector

                  dt13 = t3 - t1
                  dt32 = t2 - t3

                  del_arw = qk / real(np*np)

                  do jp = 1,np

                     facj1 = real(jp-1) / real(np)
                     facj2 = real(jp)   / real(np)

                     do ip = 1,jp

                        faci1 = real(ip-1) / real(np)
                        faci2 = real(ip)   / real(np)

                        top1 = t1 + facj1 * dt13 + faci1 * dt32
                        top2 = t1 + facj1 * dt13 + faci2 * dt32
                        top3 = t1 + facj2 * dt13 + faci1 * dt32
                        top4 = t1 + facj2 * dt13 + faci2 * dt32

                        do ipair = 1,2

                           if (ipair == 1) then
                              topp = (top1 + top3 + top4) / 3.
                           else
                              topp = (top1 + top2 + top4) / 3.
                           endif

                           if (topp < zm(k)) then
                              arw (k,iw) = arw (k,iw) + del_arw
                              volt(k,iw) = volt(k,iw) + del_arw  &
                                         * (zm(k) - max(topp,zm(k-1)))
                           endif

                           if (ip == jp) exit 

                        enddo  ! ipair

                     enddo  ! ip

                  enddo  ! jp

               endif  ! t1,t2,t3 vs zm

            enddo  ! k

         enddo  ! ihalf

      enddo  ! jv (sector)

   enddo  ! j
   !$omp end do

   !$omp end parallel

else  ! isfcl = 1

! LEAF and SEA are being used.  

! Adjust topography information that was read from LANDFILE and SEAFILE,
! if necessary, to prevent values less than lowest model level zm(1)

   land%zm(2:nml) = max(land%zm(2:nml), zm(1))
   sea%zm (2:nms) = max( sea%zm(2:nms), zm(1))

!Fill TOPM and TOPW from surface file topography.
! Determine and initialize flux cells for entire model domain.
! Initialize ARW and VOLT from surface file topography.

   write(io6,'(/,a)') 'ctrvols calling init_fluxcells'
   call init_fluxcells()

endif

! Lateral boundary copy of ARW, VOLT

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   arw (:,iw) = arw (:,iwp)
   volt(:,iw) = volt(:,iwp)

enddo
call rsub('W',jtw_lbcp)

! ARV

write(io6,*) 'Defining control volume areas'

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iv,im1,im2,iw1,iw2,k,hmin,hmax,km)
do j = 1,jtab_v(jtv_grid)%jend(1); iv = jtab_v(jtv_grid)%iv(j)
   im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
call qsub('V',iv)

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
call rsub('V',1)

! Set ARV to zero at non-topo walls
! [ONLY FOR ISFCL = 0]

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_wall)%jend(1); iv = jtab_v(jtv_wall)%iv(j)
!----------------------------------------------------------------------
call qsub('V',iv)

   do k = 2,nza
      arv(k,iv) = 0.
   enddo
   lpv(iv) = nza

enddo
call rsub('V',3)

! Lateral boundary copy of ARV

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_lbcp)%jend(1); iv = jtab_v(jtv_lbcp)%iv(j)
   ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------
call qsub('V',iv)

   arv(:,iv) = arv(:,ivp)

enddo
call rsub('V',jtv_lbcp)

! Topographic adjustments to ARW and VOLT...

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,k,facw,jv,iv)
do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Number of vertices of IW polygon

   npoly = itab_w(iw)%npoly

! Loop over vertical levels

   do k = 2,nza

! Close top area of T cell if volume is zero

      if (volt(k,iw) < 1.e-9) then
         arw (k,iw) = 0.
         volt(k,iw) = 1.e-9
      endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!go to 1
! Option for stability: expand volt if too small relative to any grid cell face

      volt(k,iw) = max(volt(k,iw), 0.5_8 * real(dzt(k) * (arw(k,iw) + arw(k-1,iw)), 8))

! Loop over faces of IW polygon

      facw = 0.0_8

      do jv = 1,npoly
         iv = itab_w(iw)%iv(jv)

         ! new way: each V face contributes to open up its fraction farv 
         ! of the current polygon:
         facw = facw + itab_w(iw)%farv(jv) * arv(k,iv) / (dnu(iv) * dzt(k))

         ! old way: each V face contributes to open up the entire polygon
         ! v_height = arv(k,iv) / dnu(iv)
         ! volt(k,iw) = max( volt(k,iw), real(arw0(iw),8) * v_height)
      enddo

      facw = max( min(facw,1.0_8), 0.0_8)
      volt(k,iw) = max( volt(k,iw), facw * real(arw0(iw),8) * real(dzt(k),8) )

! Reset arw(k,iw) if we increased volt to avoid hollow box

      if (volt(k,iw) > 1.e-9) then
         arw(k,iw) = max( arw(k,iw), real( volt(k,iw) / dzt(k)))
      endif

! Make sure arw increases with height if we have adjusted it

      if (k < nza-2) then
         arw(k+1,iw) = max( arw(k+1,iw), arw(k,iw))
      endif

1 continue

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   enddo  ! k

! Set arw = 0 for bottom (k = 1), wall-on-top (k = nza-1), 
! and top (k = nza) levels

   arw(1,iw) = 0.
   arw(nza-1,iw) = 0.
   arw(nza,iw) = 0. 

enddo
!$omp end parallel do
call rsub('W',3)

! Lateral boundary copy of ARW, VOLT

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   arw (:,iw) = arw (:,iwp)
   volt(:,iw) = volt(:,iwp)

enddo
call rsub('W',jtw_lbcp)

! LPM, LPW, LPV

lpm(2:nma) = nza-1
lpv(2:nva) = nza-1
lpw(2:nwa) = nza-1

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_m(jtm_grid)%jend(1); im = jtab_m(jtm_grid)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)

   iv1 = itab_m(im)%iv(1)
   iv2 = itab_m(im)%iv(2)
   iv3 = itab_m(im)%iv(3)

   iw1 = itab_m(im)%iw(1)
   iw2 = itab_m(im)%iw(2)
   iw3 = itab_m(im)%iw(3)

! Loop over vertical levels from top to bottom

   do k = nza-1,2,-1

      if (topm(im) >= zm(k)) exit

      if (arv(k,iv1) < .005 * dnu(iv1) * dzt(k)) exit
      if (arv(k,iv2) < .005 * dnu(iv2) * dzt(k)) exit
      if (arv(k,iv2) < .005 * dnu(iv2) * dzt(k)) exit

      if (k < nza-1) then
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
call rsub('M',1)

! Close all T cells that are below LPW

nsw_max = 1

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Number of vertices of IW polygon

   npoly = itab_w(iw)%npoly

   do k = 2,nza

      if (k < lpw(iw)) then
         arw (k,iw) = 0.
         volt(k,iw) = 1.e-9
      endif

   enddo  ! k

! Set arw = 0 for bottom (k = 1), wall-on-top (k = nza-1), 
! and top (k = mza) levels

   arw(1,iw) = 0.   
   arw(nza-1,iw) = 0.   
   arw(nza,iw) = 0.   

! Initialize LSW(IW)

   lsw(iw) = 1

! Loop over vertical levels from top to bottom

   do k = nza-1,2,-1

! Increase LSW if K-1 W level intersects topography in this cell

      if (arw(k,iw) > 0. .and. arw(k,iw) < .999 * arw0(iw)) then
         lsw(iw) = lsw(iw) + 1
      endif

   enddo  ! k

! Increase nsw_max if necessary

   if (lsw(iw) > nsw_max) nsw_max = lsw(iw)

enddo
call rsub('W',3)

! In case ARW has been reset to 0 anywhere (because it was nearly zero), 
! transfer the seaflux and landflux cell values to KW = LPW(IW).
! Also compute fractional flux cell areas per vertical level

if (isfcl == 1) then

   allocate(area_kw_sum(nsw_max,nwa))
   area_kw_sum(:,:) = 0.0

   do isf = 2,nseaflux
      iw = seaflux(isf)%iw
      kw = seaflux(isf)%kw

      if (kw < lpw(iw)) seaflux(isf)%kw = lpw(iw)
      if (kw > lpw(iw) + lsw(iw) - 1) seaflux(isf)%kw = lpw(iw) + lsw(iw) - 1

      ks = seaflux(isf)%kw - lpw(iw) + 1
      area_kw_sum(ks,iw) = area_kw_sum(ks,iw) + seaflux(isf)%area
   enddo

   do ilf = 2,nlandflux
      iw = landflux(ilf)%iw
      kw = landflux(ilf)%kw

      if (kw < lpw(iw)) landflux(ilf)%kw = lpw(iw)
      if (kw > lpw(iw) + lsw(iw) - 1) landflux(ilf)%kw = lpw(iw) + lsw(iw) - 1

      ks = landflux(ilf)%kw - lpw(iw) + 1
      area_kw_sum(ks,iw) = area_kw_sum(ks,iw) + landflux(ilf)%area
   enddo

   do isf = 2, nseaflux
      kw = seaflux(isf)%kw
      iw = seaflux(isf)%iw
      ks = kw - lpw(iw) + 1
      seaflux(isf)%arf_kw = seaflux(isf)%area / area_kw_sum(ks,iw)
   enddo

   do ilf = 2, nlandflux
      kw = landflux(ilf)%kw
      iw = landflux(ilf)%iw
      ks = kw - lpw(iw) + 1
      landflux(ilf)%arf_kw = landflux(ilf)%area / area_kw_sum(ks,iw)
   enddo

   deallocate(area_kw_sum)

endif

! Expand ARW and VOLT with height for spherical geometry, and compute VOLTI

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_grid)%jend(1); iw = jtab_w(jtw_grid)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

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

      volti(k,iw) = 1. / volt(k,iw)

   enddo  ! k

! VOLWI from VOLT

   do k = 2,nza
      kp = min(k+1,nza)
      volwi(k,iw) = 2. / (volt(k,iw) + volt(kp,iw)) 
   enddo

! modify volwi for lpw and lpw-1 levels

   ka = lpw(iw)
   volwi(ka,iw)   = 1. / (volt(ka,iw) + .5 * volt(ka+1,iw))
   volwi(ka-1,iw) = 1. / (.5 * volt(ka,iw))

! Set arw = 0 for bottom (k = 1), wall-on-top (k = nza-1), 
! and top (k = nza) levels

   arw(1,iw) = 0.   
   arw(nza-1,iw) = 0.   
   arw(nza,iw) = 0.   

enddo
call rsub('W',3)

! Lateral boundary copy of ARW, VOLT, VOLTI, VOLWI, LSW

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_lbcp)%jend(1); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   arw  (:,iw) = arw  (:,iwp)
   volt (:,iw) = volt (:,iwp)
   volti(:,iw) = volti(:,iwp)
   volwi(:,iw) = volwi(:,iwp)
   lsw    (iw) = lsw    (iwp)

enddo
call rsub('W',jtw_lbcp)

! Expand ARV with height for spherical geometry; Compute VOLVI

volvi(1:nza,1:nva) = 1.e-9

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_grid)%jend(1); iv = jtab_v(jtv_grid)%iv(j)
iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
call qsub('V',iv)

   farw1 = itab_v(iv)%farw(1)
   farw2 = itab_v(iv)%farw(2)

   if (mdomain < 2) then
      do k = 2,nza
         arv(k,iv) = arv(k,iv) * zfact(k)
      enddo
   endif

   do k = 2, nza-1
      volvi(k,iv) = 1. / (farw1 * volt(k,iw1) + farw2 * volt(k,iw2))
   enddo

   volvi(1  ,iv) = 1.e9
   volvi(nza,iv) = 1.e9

enddo
call rsub('V',4)

! Lateral boundary copy of ARV, VOLVI

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_lbcp)%jend(1); iv = jtab_v(jtv_lbcp)%iv(j)
   ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------
call qsub('V',iv)

   arv  (:,iv) = arv  (:,ivp)
   volvi(:,iv) = volvi(:,ivp)

enddo
call rsub('V',jtv_lbcp)

return
end subroutine ctrlvols_hex
