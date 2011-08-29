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
subroutine delaunay()

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd, itab_m, itab_u, itab_w,  &
                       mrls, alloc_itabs

use mem_grid,    only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa,  &
                       xem, yem, zem, xew, yew, zew,  &
                       alloc_xyzem, alloc_xyzew
use misc_coms,   only: io6, mdomain, meshtype
use consts_coms, only: pi2, erad, erador5

implicit none

! Copy grid dimensions

mma = nma
mua = nua
nva = nua
mva = nva
mwa = nwa

! Allocate Delaunay set of itabs

call alloc_itabs(meshtype,nma,nua,nva,nwa)

! Copy itabs

itab_m(1:nma) = itab_md(1:nma)
itab_u(1:nua) = itab_ud(1:nua)
itab_w(1:nwa) = itab_wd(1:nwa)

deallocate(itab_md,itab_ud,itab_wd)

! Allocate XEW,YEW,ZEW arrays

call alloc_xyzew(nwa)

return
end subroutine delaunay

!===============================================================================

subroutine grid_geometry_tri()

use mem_ijtabs,  only: itab_m, itab_u, itab_w
use mem_grid,    only: nza, mma, mua, mwa, xeu, yeu, zeu, xem, yem, zem,  &
                       xew, yew, zew, unx, uny, unz, vnx, vny, vnz,  &
                       wnx, wny, wnz, dnu, dniu, dnv, dniv, arm0, arw0,  &
                       glonw, glatw, glonm, glatm, glatu, glonu
use misc_coms,   only: io6, mdomain, grdlat, grdlon
use consts_coms, only: erad, erad2, piu180

implicit none

integer :: im,iu,iw
integer :: im1,im2,im3
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: itpn

real :: expansion
real :: raxis
real :: hper

real :: xiw,yiw,ziw
real :: xij(6),yij(6),zij(6)
real :: xw1,xw2,xw3,xw4,xw5,xw6,yw1,yw2,yw3,yw4,yw5,yw6
real :: vecjx,vecjy,vecjz
real :: vecmx,vecmy,vecmz
real :: vecmjx,vecmjy,vecmjz
real :: xi,yi,xj(6),yj(6)

real :: ef,x12,y12,z12,x34,y34,z34
real :: b11,b21,b31,b12,b22,b32,b13,b23,b33  &
       ,c11,c21,c31,c12,c22,c32,c13,c23,c33  &
       ,d11,d21,d31,d12,d22,d32,d13,d23,d33
real :: dnupgf
real :: s1, s2, s3
real :: xm1,xm2,ym1,ym2,tux,tuy

real :: xec(mua),yec(mua),zec(mua)
real :: xc0,xc5,xc6,xc7,xc8,xc9,xc10,xc11,xc12
real :: yc0,yc5,yc6,yc7,yc8,yc9,yc10,yc11,yc12
real :: xcw3,xcw4,xcw5,xcw6
real :: ycw3,ycw4,ycw5,ycw6

ef = 1.01  ! radial expansion factor (from earth center) for defining 
           ! auxiliary point for computing unit normal to U face

! Loop over all U points

do iu = 2,mua

! Fill global index (replaced later if this run is parallel)

   itab_u(iu)%iuglobe = iu

! M-point indices of two end points of U segment

   im1 = itab_u(iu)%im(1)
   im2 = itab_u(iu)%im(2)

! Fill U point coordinates from M point coordinates

   xeu(iu) = .5 * (xem(im1) + xem(im2))
   yeu(iu) = .5 * (yem(im1) + yem(im2))
   zeu(iu) = .5 * (zem(im1) + zem(im2))

! If mdomain <= 1, push U point coordinates out to earth radius

   if (mdomain <= 1) then

      expansion = erad / sqrt(xeu(iu) ** 2  &
                            + yeu(iu) ** 2  &
                            + zeu(iu) ** 2  )

      xeu(iu) = xeu(iu) * expansion
      yeu(iu) = yeu(iu) * expansion
      zeu(iu) = zeu(iu) * expansion

   endif

! Fill latitude and longitude of U point

   if (mdomain <= 1) then
      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)

      glatu(iu) = atan2(zeu(iu),raxis)   * piu180
      glonu(iu) = atan2(yeu(iu),xeu(iu)) * piu180
   else
      glatu(iu) = grdlat(1,1)  ! want it this way?
      glonu(iu) = grdlon(1,1)  ! want it this way?
   endif

! Length of U side

   dnv(iu) = sqrt( (xem(im1) - xem(im2))**2  &
                 + (yem(im1) - yem(im2))**2  &
                 + (zem(im1) - zem(im2))**2  )
                 
   dniv(iu) = 1. / dnv(iu)
                 
! Convert to geodesic arc length

!   dnv(iu) = erad2 * asin(dnv(iu) / erad2)

! Fill U face unit normal vector components

   if (mdomain <= 1) then   ! Spherical geometry case

      call unit_normal(xem(im2),yem(im2),zem(im2)           &
                      ,xem(im1),yem(im1),zem(im1)           &
                      ,ef*xem(im1),ef*yem(im1),ef*zem(im1)  &
                      ,unx(iu),uny(iu),unz(iu))   

   else                     ! Cartesian geometry case

      call unit_normal(xem(im2),yem(im2),zem(im2)         &
                      ,xem(im1),yem(im1),zem(im1)         &
                      ,xem(im1),yem(im1),zem(im1) + 1.e3  &
                      ,unx(iu),uny(iu),unz(iu))   

   endif                

! Fill U face unit horizontal tangential vector components

   vnx(iu) = (xem(im2) - xem(im1)) * dniv(iu)
   vny(iu) = (yem(im2) - yem(im1)) * dniv(iu)
   vnz(iu) = (zem(im2) - zem(im1)) * dniv(iu)

enddo

! Loop over all W points

do iw = 2,mwa

! Fill global index (replaced later if this run is parallel)

   itab_w(iw)%iwglobe = iw

! Indices of 3 M points surrounding W point

   im1 = itab_w(iw)%im(1)
   im2 = itab_w(iw)%im(2)
   im3 = itab_w(iw)%im(3)

! Indices of 3 U points surrounding W point

   iu1 = itab_w(iw)%iu(1)
   iu2 = itab_w(iw)%iu(2)
   iu3 = itab_w(iw)%iu(3)

! Fill W face unit normal vector components

   call unit_normal(xem(im1),yem(im1),zem(im1)  &
                   ,xem(im2),yem(im2),zem(im2)  &
                   ,xem(im3),yem(im3),zem(im3)  &
                   ,wnx(iw) ,wny(iw) ,wnz(iw)   )

! Fill W point coordinates from coordinates of 3 M points

! OPTION 1: barycentric coordinates for IW point

   xew(iw) = (xem(im1) + xem(im2) + xem(im3)) / 3.
   yew(iw) = (yem(im1) + yem(im2) + yem(im3)) / 3.
   zew(iw) = (zem(im1) + zem(im2) + zem(im3)) / 3.

! OPTION 2: Circumcentric coordinates for IW point

!   xew(iw) = wnx(iw) * erad
!   yew(iw) = wny(iw) * erad
!   zew(iw) = wnz(iw) * erad

! If mdomain <= 1, push W point coordinates out to earth radius

   if (mdomain <= 1) then

      expansion = erad / sqrt(xew(iw) ** 2  &
                            + yew(iw) ** 2  &
                            + zew(iw) ** 2  )

      xew(iw) = xew(iw) * expansion
      yew(iw) = yew(iw) * expansion
      zew(iw) = zew(iw) * expansion

   endif

! Fill latitude and longitude of W point

   if (mdomain <= 1) then
      raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)

      glatw(iw) = atan2(zew(iw),raxis)   * piu180
      glonw(iw) = atan2(yew(iw),xew(iw)) * piu180
   else
      glatw(iw) = grdlat(1,1)  ! want it this way?
      glonw(iw) = grdlon(1,1)  ! want it this way?
   endif

! IW triangle area at sea level

   hper = .5 * (dnv(iu1) + dnv(iu2) + dnv(iu3))  ! half perimeter
         
   arw0(iw) = sqrt(hper * (hper - dnv(iu1))  &
                        * (hper - dnv(iu2))  &
                        * (hper - dnv(iu3))  )

enddo

! Loop over all M points

do im = 2,mma

! Latitude and longitude at M points

   if (mdomain <= 1) then
      raxis = sqrt(xem(im) ** 2 + yem(im) ** 2)  ! dist from earth axis
      glatm(im) = atan2(zem(im),raxis)   * piu180
      glonm(im) = atan2(yem(im),xem(im)) * piu180
   else
      glatm(im) = grdlat(1,1)  ! want it this way?
      glonm(im) = grdlon(1,1)  ! want it this way?
   endif

! Fill global index (replaced later if this run is parallel)

   itab_m(im)%imglobe = im

! Set area of dual cell surrounding M to zero prior to summation over triangles

   arm0(im) = 0.

! Loop over all U neighbors of M

   do itpn = 1,itab_m(im)%npoly

      iu = itab_m(im)%iu(itpn)

      iw1 = itab_u(iu)%iw(1)
      iw2 = itab_u(iu)%iw(2)

! Contribution to dual-cell area around M point of triangle M-IW1-IW2

      s1 = sqrt( (xew(iw1) - xew(iw2))**2  &
               + (yew(iw1) - yew(iw2))**2  &
               + (zew(iw1) - zew(iw2))**2  )

      s2 = sqrt( (xem(im) - xew(iw2))**2  &
               + (yem(im) - yew(iw2))**2  &
               + (zem(im) - zew(iw2))**2  )
            
      s3 = sqrt( (xew(iw1) - xem(im))**2  &
               + (yew(iw1) - yem(im))**2  &
               + (zew(iw1) - zem(im))**2  )
            
      hper = .5 * (s1 + s2 + s3)  ! half perimeter of triangle
 
      arm0(im) = arm0(im) + sqrt(hper * (hper - s1) * (hper - s2) * (hper - s3))

   enddo

enddo

!! Loop over all U points

do iu = 2,mua
   iw1 = itab_u(iu)%iw(1)
   iw2 = itab_u(iu)%iw(2)

   xec(iu) = (arw0(iw1) * xew(iw1) + arw0(iw2) * xew(iw2)) / (arw0(iw1) + arw0(iw2))
   yec(iu) = (arw0(iw1) * yew(iw1) + arw0(iw2) * yew(iw2)) / (arw0(iw1) + arw0(iw2))
   zec(iu) = (arw0(iw1) * zew(iw1) + arw0(iw2) * zew(iw2)) / (arw0(iw1) + arw0(iw2))

! If mdomain <= 1, push ec point coordinates out to earth radius

   if (mdomain <= 1) then

      expansion = erad / sqrt(xec(iu) ** 2  &
                            + yec(iu) ** 2  &
                            + zec(iu) ** 2  )

      xec(iu) = xec(iu) * expansion
      yec(iu) = yec(iu) * expansion
      zec(iu) = zec(iu) * expansion

   endif

enddo

! Loop over all U points

do iu = 2,mua

   im1 = itab_u(iu)%im(1)
   im2 = itab_u(iu)%im(2)

   iu1  = itab_u(iu)%iu(1)
   iu2  = itab_u(iu)%iu(2)
   iu3  = itab_u(iu)%iu(3)
   iu4  = itab_u(iu)%iu(4)
   iu5  = itab_u(iu)%iu(5)
   iu6  = itab_u(iu)%iu(6)
   iu7  = itab_u(iu)%iu(7)
   iu8  = itab_u(iu)%iu(8)
   iu9  = itab_u(iu)%iu(9)
   iu10 = itab_u(iu)%iu(10)
   iu11 = itab_u(iu)%iu(11)
   iu12 = itab_u(iu)%iu(12)

   iw1 = itab_u(iu)%iw(1)
   iw2 = itab_u(iu)%iw(2)
   iw3 = itab_u(iu)%iw(3)
   iw4 = itab_u(iu)%iw(4)
   iw5 = itab_u(iu)%iw(5)
   iw6 = itab_u(iu)%iw(6)

! Project neighbor U & W unit vectors onto U horizontal tangential vector

   if (iw1 > 1)                                       &
   call matrix_3x3(unx(iu1),unx(iu2),wnx(iw1)         &
                  ,uny(iu1),uny(iu2),wny(iw1)         &
                  ,unz(iu1),unz(iu2),wnz(iw1)         &
                  ,vnx(iu),vny(iu),vnz(iu)            &
                  ,b11,b12,b13                        )

   if (iw2 > 1)                                       &
   call matrix_3x3(unx(iu3),unx(iu4),wnx(iw2)         &
                  ,uny(iu3),uny(iu4),wny(iw2)         &
                  ,unz(iu3),unz(iu4),wnz(iw2)         &
                  ,vnx(iu),vny(iu),vnz(iu)            &
                  ,c11,c12,c13                        )

   if (iw1 > 1 .and. iw2 > 1) then
      itab_u(iu)%tuu(1) = .5 * b11
      itab_u(iu)%tuu(2) = .5 * b12
      itab_u(iu)%tuu(3) = .5 * c11
      itab_u(iu)%tuu(4) = .5 * c12
   elseif (iw1 > 1) then
      itab_u(iu)%tuu(1) = b11
      itab_u(iu)%tuu(2) = b12
   elseif (iw2 > 1) then
      itab_u(iu)%tuu(3) = c11
      itab_u(iu)%tuu(4) = c12
   endif

! Defaults:

   itab_u(iu)%crossmm = .5
   itab_u(iu)%crossww = .5

   itab_u(iu)%gcf36 = 0.
   itab_u(iu)%gcf45 = 0.

   if (iw1 > 1 .and. iw2 > 1) then

      call e_ps(xem(im1),yem(im1),zem(im1),glatu(iu),glonu(iu),xm1,ym1)
      call e_ps(xem(im2),yem(im2),zem(im2),glatu(iu),glonu(iu),xm2,ym2)
      call e_ps(xew(iw1),yew(iw1),zew(iw1),glatu(iu),glonu(iu),xw1,yw1)
      call e_ps(xew(iw2),yew(iw2),zew(iw2),glatu(iu),glonu(iu),xw2,yw2)

      tux = (xm2 - xm1) / sqrt((xm2 - xm1)**2 + (ym2 - ym1)**2) 
      tuy = (ym2 - ym1) / sqrt((xm2 - xm1)**2 + (ym2 - ym1)**2) 

! Crossmm and crossww values

!      itab_u(iu)%crossmm = (ym1 * xw2 - yw2 * xm1) &
!                         / (yw2 * xm2 - yw2 * xm1 - xw2 * ym2 + xw2 * ym1)

!      itab_u(iu)%crossww = (xm2 * ym1 - xm1 * ym2) &
!                         / (yw2 * xm2 - yw2 * xm1 - xw2 * ym2 + xw2 * ym1)

      itab_u(iu)%crossmm &
         = ((xw2 - xw1) * ym1 + (xm1 - xw2) * yw1 + (xw1 - xm1) * yw2) &
         / ((xm2 - xm1) * (yw2 - yw1) + (xw2 - xw1) * (ym1 - ym2))
      
      itab_u(iu)%crossww &
         = ((xm1 - xm2) * yw1 + (xm2 - xw1) * ym1 + (xw1 - xm1) * ym2) &
         / ((xm2 - xm1) * (yw2 - yw1) + (xw2 - xw1) * (ym1 - ym2))
      
! The above code is the same as...

!tryit   call intersect(xm1,ym1,xm2,ym2,xw1,yw1,xw2,yw2, &
!tryit      itab_u(iu)%crossmm,itab_u(iu)%crossww,iflag)


! Transverse gradient coefficients 

      if (iw3 > 1 .and. iw4 > 1) then

         call e_ps(xew(iw3),yew(iw3),zew(iw3),glatu(iu),glonu(iu),xw3,yw3)
         call e_ps(xew(iw4),yew(iw4),zew(iw4),glatu(iu),glonu(iu),xw4,yw4)
   
         call matrix_2x2(xw3-xw1,xw4-xw1,yw3-yw1,yw4-yw1,tux,tuy, &
            itab_u(iu)%guw(1),itab_u(iu)%guw(2))

      endif

      if (iw5 > 1 .and. iw6 > 1) then

         call e_ps(xew(iw5),yew(iw5),zew(iw5),glatu(iu),glonu(iu),xw5,yw5)
         call e_ps(xew(iw6),yew(iw6),zew(iw6),glatu(iu),glonu(iu),xw6,yw6)
   
         call matrix_2x2(xw5-xw2,xw6-xw2,yw5-yw2,yw6-yw2,tux,tuy, &
            itab_u(iu)%guw(3),itab_u(iu)%guw(4))

      endif

! 1 Jan 2010 - interpolation coefficients to get U flux from U momentum

! Horizontal center of mass of U control volume
! (Horizontal center of U edge is at origin: (xu,yu) = 0)

      call e_ps(xec(iu),yec(iu),zec(iu),glatu(iu),glonu(iu),xc0,yc0)

      call e_ps(xec(iu5),yec(iu5),zec(iu5),glatu(iu),glonu(iu),xc5,yc5)
      call e_ps(xec(iu6),yec(iu6),zec(iu6),glatu(iu),glonu(iu),xc6,yc6)
      call e_ps(xec(iu7),yec(iu7),zec(iu7),glatu(iu),glonu(iu),xc7,yc7)
      call e_ps(xec(iu8),yec(iu8),zec(iu8),glatu(iu),glonu(iu),xc8,yc8)
      call e_ps(xec(iu9),yec(iu9),zec(iu9),glatu(iu),glonu(iu),xc9,yc9)
      call e_ps(xec(iu10),yec(iu10),zec(iu10),glatu(iu),glonu(iu),xc10,yc10)
      call e_ps(xec(iu11),yec(iu11),zec(iu11),glatu(iu),glonu(iu),xc11,yc11)
      call e_ps(xec(iu12),yec(iu12),zec(iu12),glatu(iu),glonu(iu),xc12,yc12)

!mod
      call e_ps(xeu(iu5),yeu(iu5),zeu(iu5),glatu(iu),glonu(iu),xc5,yc5)
      call e_ps(xeu(iu6),yeu(iu6),zeu(iu6),glatu(iu),glonu(iu),xc6,yc6)
      call e_ps(xeu(iu7),yeu(iu7),zeu(iu7),glatu(iu),glonu(iu),xc7,yc7)
      call e_ps(xeu(iu8),yeu(iu8),zeu(iu8),glatu(iu),glonu(iu),xc8,yc8)
      call e_ps(xeu(iu9),yeu(iu9),zeu(iu9),glatu(iu),glonu(iu),xc9,yc9)
      call e_ps(xeu(iu10),yeu(iu10),zeu(iu10),glatu(iu),glonu(iu),xc10,yc10)
      call e_ps(xeu(iu11),yeu(iu11),zeu(iu11),glatu(iu),glonu(iu),xc11,yc11)
      call e_ps(xeu(iu12),yeu(iu12),zeu(iu12),glatu(iu),glonu(iu),xc12,yc12)

      xcw3 = .5 * (xc5 + xc6)
      xcw4 = .5 * (xc7 + xc8)
      xcw5 = .5 * (xc9 + xc10)
      xcw6 = .5 * (xc11 + xc12)

      ycw3 = .5 * (yc5 + yc6)
      ycw4 = .5 * (yc7 + yc8)
      ycw5 = .5 * (yc9 + yc10)
      ycw6 = .5 * (yc11 + yc12)

      if (iw3 > 1 .and. iw4 > 1 .and. iw5 > 1 .and. iw6 > 1) then

! gcf36 and gcf45 were used in test simulations in version 4.0i, but are 
! not used in version 4.0j.

         call matrix_2x2(xcw6-xcw3,xcw5-xcw4,ycw6-ycw3,ycw5-ycw4,-xc0,-yc0, &
            itab_u(iu)%gcf36,itab_u(iu)%gcf45)

      endif

   endif

! Skip this U point if iw1 < 2 or iw2 < 2

   if (iw1 < 2 .or. iw2 < 2) cycle

   itab_u(iu)%mrlu = max(itab_w(iw1)%mrlw,itab_w(iw2)%mrlw)

! Project neighbor U & W unit vectors onto U unit vector

   if (iw3 > 1)                                &
   call matrix_3x3(unx(iu5),unx(iu6),wnx(iw3)  &
                  ,uny(iu5),uny(iu6),wny(iw3)  &
                  ,unz(iu5),unz(iu6),wnz(iw3)  &
                  ,unx(iu),uny(iu),unz(iu)     &
                  ,itab_u(iu)%fuu(5)           &
                  ,itab_u(iu)%fuu(6)           &
                  ,itab_u(iu)%fuw(3)           )

   if (iw4 > 1)                                &
   call matrix_3x3(unx(iu7),unx(iu8),wnx(iw4)  &
                  ,uny(iu7),uny(iu8),wny(iw4)  &
                  ,unz(iu7),unz(iu8),wnz(iw4)  &
                  ,unx(iu),uny(iu),unz(iu)     &
                  ,itab_u(iu)%fuu(7)           &
                  ,itab_u(iu)%fuu(8)           &
                  ,itab_u(iu)%fuw(4)           )

   if (iw5 > 1)                                 &
   call matrix_3x3(unx(iu9),unx(iu10),wnx(iw5)  &
                  ,uny(iu9),uny(iu10),wny(iw5)  &
                  ,unz(iu9),unz(iu10),wnz(iw5)  &
                  ,unx(iu),uny(iu),unz(iu)      &
                  ,itab_u(iu)%fuu(9)            &
                  ,itab_u(iu)%fuu(10)           &
                  ,itab_u(iu)%fuw(5)            )

   if (iw6 > 1)                                  &
   call matrix_3x3(unx(iu11),unx(iu12),wnx(iw6)  &
                  ,uny(iu11),uny(iu12),wny(iw6)  &
                  ,unz(iu11),unz(iu12),wnz(iw6)  &
                  ,unx(iu),uny(iu),unz(iu)       &
                  ,itab_u(iu)%fuu(11)            &
                  ,itab_u(iu)%fuu(12)            &
                  ,itab_u(iu)%fuw(6)             )

! Divide fuw3, fuw4, fuw5, and fuw6 by 2 for use with two W points

   if (iw3 > 1) itab_u(iu)%fuw(3) = .5 * itab_u(iu)%fuw(3)
   if (iw4 > 1) itab_u(iu)%fuw(4) = .5 * itab_u(iu)%fuw(4)
   if (iw5 > 1) itab_u(iu)%fuw(5) = .5 * itab_u(iu)%fuw(5)
   if (iw6 > 1) itab_u(iu)%fuw(6) = .5 * itab_u(iu)%fuw(6)

! Find earth velocity coefficients for Coriolis force

   call matinv3x3(unx(iu1),uny(iu1),unz(iu1)          &
                 ,unx(iu2),uny(iu2),unz(iu2)          &
                 ,wnx(iw1),wny(iw1),wnz(iw1)          &
                 ,b11,b21,b31,b12,b22,b32,b13,b23,b33 )

   call matinv3x3(unx(iu3),uny(iu3),unz(iu3)          &
                 ,unx(iu4),uny(iu4),unz(iu4)          &
                 ,wnx(iw2),wny(iw2),wnz(iw2)          &
                 ,c11,c21,c31,c12,c22,c32,c13,c23,c33 )

   itab_u(iu)%vxu_u(1) = b11
   itab_u(iu)%vxu_u(2) = b21
   itab_u(iu)%vxw_u(1) = b31 / 2.

   itab_u(iu)%vyu_u(1) = b12
   itab_u(iu)%vyu_u(2) = b22
   itab_u(iu)%vyw_u(1) = b32 / 2.

   itab_u(iu)%vxu_u(3) = c11
   itab_u(iu)%vxu_u(4) = c21
   itab_u(iu)%vxw_u(2) = c31 / 2.

   itab_u(iu)%vyu_u(3) = c12
   itab_u(iu)%vyu_u(4) = c22
   itab_u(iu)%vyw_u(2) = c32 / 2.

! U point PGF coefficients for full 6-point stencil

   x12 = xew(iw2) - xew(iw1)
   y12 = yew(iw2) - yew(iw1)
   z12 = zew(iw2) - zew(iw1)

   dnu(iu) = unx(iu) * x12 + uny(iu) * y12 + unz(iu) * z12
   dniu(iu) = 1. / dnu(iu)

   dnupgf = (arw0(iw1) + arw0(iw2)) / dnv(iu)


!!!!!!!!! ALTERNATE FORM
   dnupgf = 1.
!!!!!!!!!!!!!!!!!!!!!!!!   

   if (iw3 > 1 .and. iw4 > 1 .and. iw5 > 1 .and. iw6 > 1) then

      x34 = xew(iw3) + xew(iw5) - xew(iw4) - xew(iw6)
      y34 = yew(iw3) + yew(iw5) - yew(iw4) - yew(iw6)
      z34 = zew(iw3) + zew(iw5) - zew(iw4) - zew(iw6)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif
                 
      itab_u(iu)%pgc12 = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc45 = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)
      itab_u(iu)%pgc63 = itab_u(iu)%pgc45
      
   endif

! U point PGF coefficients for (1-2-4-5) 4-point stencil

   if (iw4 > 1 .and. iw5 > 1) then

      x34 = xew(iw5) - xew(iw4)
      y34 = yew(iw5) - yew(iw4)
      z34 = zew(iw5) - zew(iw4)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif

      itab_u(iu)%pgc12b = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc45b = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)

   endif

! U point PGF coefficients for (1-2-6-3) 4-point stencil

   if (iw3 > 1 .and. iw6 > 1) then

      x34 = xew(iw3) - xew(iw6)
      y34 = yew(iw3) - yew(iw6)
      z34 = zew(iw3) - zew(iw6)

      if (mdomain <= 1) then ! Spherical geometry case

         call matinv3x3(x12,y12,z12  &
                       ,x34,y34,z34  &
                       ,y12*z34-z12*y34,z12*x34-x12*z34,x12*y34-y12*x34   &
                       ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

      else                   ! Cartesian geometry case

         call matinv2x2(x12,y12,x34,y34,b11,b21,b12,b22)

         b13 = 0.
         b23 = 0.
      endif
                 
      itab_u(iu)%pgc12c = dnupgf * (unx(iu) * b11 + uny(iu) * b12 + unz(iu) * b13)
      itab_u(iu)%pgc63c = dnupgf * (unx(iu) * b21 + uny(iu) * b22 + unz(iu) * b23)
      
   endif

! U point PGF coefficients for (1-2) 2-point stencil

   if (mdomain <= 1) then ! Spherical geometry
      itab_u(iu)%pgc12d = dnupgf / (unx(iu) * x12 + uny(iu) * y12 + unz(iu) * z12)
   else                   ! Cartesian geometry
      itab_u(iu)%pgc12d = dnupgf / (unx(iu) * x12 + uny(iu) * y12)
   endif

enddo

! Loop over all W points

do iw = 2,mwa

! Indices of 3 M points, 9 U points, and 3 W points surrounding W point

   im1 = itab_w(iw)%im(1)                   
   im2 = itab_w(iw)%im(2)                   
   im3 = itab_w(iw)%im(3)                   

   iu1 = itab_w(iw)%iu(1)                   
   iu2 = itab_w(iw)%iu(2)                   
   iu3 = itab_w(iw)%iu(3)                   
   iu4 = itab_w(iw)%iu(4)                   
   iu5 = itab_w(iw)%iu(5)                   
   iu6 = itab_w(iw)%iu(6)                   
   iu7 = itab_w(iw)%iu(7)                   
   iu8 = itab_w(iw)%iu(8)                   
   iu9 = itab_w(iw)%iu(9)

   iw1 = itab_w(iw)%iw(1)                   
   iw2 = itab_w(iw)%iw(2)                   
   iw3 = itab_w(iw)%iw(3)                   

! Project neighbor U & W unit vectors onto W unit vector

   if (iw1 > 1)                                &
   call matrix_3x3(unx(iu4),unx(iu5),wnx(iw1)  &
                  ,uny(iu4),uny(iu5),wny(iw1)  &
                  ,unz(iu4),unz(iu5),wnz(iw1)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu(4)           &
                  ,itab_w(iw)%fwu(5)           &
                  ,itab_w(iw)%fww(1)           )

   if (iw2 > 1)                                &
   call matrix_3x3(unx(iu6),unx(iu7),wnx(iw2)  &
                  ,uny(iu6),uny(iu7),wny(iw2)  &
                  ,unz(iu6),unz(iu7),wnz(iw2)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu(6)           &
                  ,itab_w(iw)%fwu(7)           &
                  ,itab_w(iw)%fww(2)           )

   if (iw3 > 1)                                &
   call matrix_3x3(unx(iu8),unx(iu9),wnx(iw3)  &
                  ,uny(iu8),uny(iu9),wny(iw3)  &
                  ,unz(iu8),unz(iu9),wnz(iw3)  &
                  ,wnx(iw),wny(iw),wnz(iw)     &
                  ,itab_w(iw)%fwu(8)           &
                  ,itab_w(iw)%fwu(9)           &
                  ,itab_w(iw)%fww(3)           )

! Divide fwu4, fwu5, fwu6, fwu7, fwu8, and fwu9 by 2 for use with two W points
                  
   if (iw1 > 1) itab_w(iw)%fwu(4) = .5 * itab_w(iw)%fwu(4)
   if (iw1 > 1) itab_w(iw)%fwu(5) = .5 * itab_w(iw)%fwu(5)
   if (iw2 > 1) itab_w(iw)%fwu(6) = .5 * itab_w(iw)%fwu(6)
   if (iw2 > 1) itab_w(iw)%fwu(7) = .5 * itab_w(iw)%fwu(7)
   if (iw3 > 1) itab_w(iw)%fwu(8) = .5 * itab_w(iw)%fwu(8)
   if (iw3 > 1) itab_w(iw)%fwu(9) = .5 * itab_w(iw)%fwu(9)

! Find horizontal unit vector components for each pair of U's in this W

   call matinv3x3(unx(iu1),uny(iu1),unz(iu1)  &
                 ,unx(iu2),uny(iu2),unz(iu2)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,b11,b21,b31,b12,b22,b32,b13,b23,b33)
                 
   call matinv3x3(unx(iu2),uny(iu2),unz(iu2)  &
                 ,unx(iu3),uny(iu3),unz(iu3)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,c11,c21,c31,c12,c22,c32,c13,c23,c33)

   call matinv3x3(unx(iu3),uny(iu3),unz(iu3)  &
                 ,unx(iu1),uny(iu1),unz(iu1)  &
                 ,wnx(iw) ,wny(iw) ,wnz(iw)   &
                 ,d11,d21,d31,d12,d22,d32,d13,d23,d33)

! Coefficients for 3D wind (for HORIZONTAL wind, don't use tri*4)

   itab_w(iw)%vxu(1) = (b11 + d21) / 3.
   itab_w(iw)%vxu(2) = (b21 + c11) / 3.
   itab_w(iw)%vxu(3) = (c21 + d11) / 3.
   itab_w(iw)%vxw    = (b31 + c31 + d31) / 6.  ! Includes div by 2 for 2 W's

   itab_w(iw)%vyu(1) = (b12 + d22) / 3.
   itab_w(iw)%vyu(2) = (b22 + c12) / 3.
   itab_w(iw)%vyu(3) = (c22 + d12) / 3.
   itab_w(iw)%vyw    = (b32 + c32 + d32) / 6.

   itab_w(iw)%vzu(1) = (b13 + d23) / 3.
   itab_w(iw)%vzu(2) = (b23 + c13) / 3.
   itab_w(iw)%vzu(3) = (c23 + d13) / 3.
   itab_w(iw)%vzw    = (b33 + c33 + d33) / 6.

!!!!!!!!!!!! special (new form that omits umc(k,*) from its own fcor)
   itab_w(iw)%vxu_w(1) = (b11 + d21) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vxu_w(2) = (b21 + c11) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vxu_w(3) = (c21 + d11) / 3. ! Includes div by 3 for 3 tri vertices

   itab_w(iw)%vyu_w(1) = (b12 + d22) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vyu_w(2) = (b22 + c12) / 3. ! Includes div by 3 for 3 tri vertices
   itab_w(iw)%vyu_w(3) = (c22 + d12) / 3. ! Includes div by 3 for 3 tri vertices
!!!!!!!!!!!! end special


!!!!!!! SAMPLE ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! (xe,ye,ze) components of 3D wind

!   um = vxu1 * umc(k,iu1) + vxu2 * umc(k,iu2) + vxu3 * umc(k,iu3)  &
!      + vxw * (wmc(k,iw) + wmc(k+1,iw))
!   vm = vyu1 * umc(k,iu1) + vyu2 * umc(k,iu2) + vyu3 * umc(k,iu3)  &
!      + vyw * (wmc(k,iw) + wmc(k+1,iw))
!   wm = vzu1 * umc(k,iu1) + vzu2 * umc(k,iu2) + vzu3 * umc(k,iu3)  &
!      + vzw * (wmc(k,iw) + wmc(k+1,iw))

! (xe,ye,ze) components of HORIZONTAL wind

!   um = vxu1 * umc(k,iu1) + vxu2 * umc(k,iu2) + vxu3 * umc(k,iu3)
!   vm = vyu1 * umc(k,iu1) + vyu2 * umc(k,iu2) + vyu3 * umc(k,iu3)
!   wm = vzu1 * umc(k,iu1) + vzu2 * umc(k,iu2) + vzu3 * umc(k,iu3)

!!!!!!!!!!! END SAMPLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

return
end subroutine grid_geometry_tri

!===============================================================================

subroutine ctrlvols_tri()

use mem_ijtabs,  only: jtab_m, jtab_u, jtab_w, itab_u, itab_w
use misc_coms,   only: io6, mdomain, itopoflg
use consts_coms, only: erad
use mem_grid,    only: nsw_max, nza, nma, nua, nwa, lpu, lcu, lpw, lsw,  &
                       topm, topw, zm, dzt, zt, zfacm, zfact, &
                       dnv, arm0, arw0, aru, arw, volt, volti, volui, volwi, &
                       xem, yem, zem, glatm, glonm
use leaf_coms,   only: isfcl,nml
use sea_coms,    only: nms
use mem_sflux,   only: init_fluxcells,nseaflux,nlandflux,seaflux,landflux
use mem_leaf,    only: land
use mem_sea,     only: sea

implicit none

integer :: j,iw,iwp,iter,iu,iup,im1,im2,k,km,i,im  &
   ,im11,im21,im12,im22,iu1,iu2,iw1,iw2,kp  &
   ,iu1a,iu1b,iu2a,iu2b  &
   ,iuo1a,iuo1b,iuo2a,iuo2b  &
   ,im3,iu3,ka
integer :: isf,ilf,kw

integer, parameter :: npass = 2
integer :: ipass

real :: hmin,hmax,arwo4,sum1,sumk  &
   ,arwo3,w1,w2,t1,t2,t3,hm,dt13,dt32,t13,t32

!!!!!!!!!!!!! special quadrature parameters

integer, parameter :: np = 30
integer :: ip,jp,ipair,npoly,jm,jm1,jm2,ju

real :: top1,top2,top3,top4
real :: topp,topm13,topm32
real :: facj1,facj2,faci1,faci2
real :: del_arw,u_height

!!!!!!!!!!!!! end special

! Check whether LEAF and SEA models are being used

if (isfcl == 0) then

! LEAF and SEA models are NOT being used.  Initialize topography without
! reading a surface file.

   if (itopoflg == 1) then

! Read TOPM from dataset

      call topo_database_read(nma,xem,yem,zem,topm,arm0,glatm,glonm)

   else

! Initialize TOPM from default value (may customize here)

      call topo_init(nma,topm,glatm,glonm,xem,yem,zem)   

   endif

! Prevent TOPM lower than lowest model level zm(1)

   topm(2:nma) = max(topm(2:nma),zm(1))

! Fill TOPW by averaging TOPM values
! Fill ARW and VOLT from topography using quadrature method

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_w(1)%jend(1); iw = jtab_w(1)%iw(j)
   im1 = itab_w(iw)%im(1); im2 = itab_w(iw)%im(2); im3 = itab_w(iw)%im(3)
!----------------------------------------------------------------------
   call qsub('W',iw)

      t1 = topm(im1)
      t2 = topm(im2)
      t3 = topm(im3)
      
      topw(iw) = (t1 + t2 + t3) / 3.      
   
! Initialize arw, and volt to zero

      arw (1:nza,iw) = 0.
      volt(1:nza,iw) = 0.

! Loop over model levels

      do k = 2,nza-1

         if (t1 < zm(k-1) .and. t2 < zm(k-1) .and. t3 < zm(k-1)) then

! Bottom of this model layer is above all 3 topography points

            arw (k,iw) = arw0(iw)
            volt(k,iw) = arw0(iw) * dzt(k)

         else
 
            arw (k,iw) = 0.
            volt(k,iw) = 1.e-9


! Interpolate 3 topography heights to np*np sub-triangles in this sector

            dt13 = t3 - t1
            dt32 = t2 - t3

            del_arw = arw0(iw) / real(np*np)

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

         endif

      enddo  ! k
      
   enddo  ! j/iw
   call rsub('W',1)

else  ! isfcl = 1

! LEAF and SEA are being used.  

! Adjust topography information that was read from LANDFILE and SEAFILE,
! if necessary, to prevent values less than lowest model level zm(1)

   land%zm(1:nml) = max(land%zm(1:nml), zm(1))
   sea%zm (1:nms) = max( sea%zm(1:nms), zm(1))

!Fill TOPM and TOPW from surface file topography.
! Determine and initialize flux cells for entire model domain.
! Initialize ARW and VOLT from surface file topography.

   write(io6,'(/,a)') 'ctrvols calling init_fluxcells'
   call init_fluxcells()

endif

! ARU, LPU, and TOPM adjustment

do ipass = 1,npass

   write(io6,*) 'Defining contol volume areas: ipass = ',ipass

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(1)%jend(1); iu = jtab_u(1)%iu(j)
      im1 = itab_u(iu)%im(1); im2 = itab_u(iu)%im(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      aru(1,iu) = 0.
      lpu(iu) = nza

      if (dnv(iu) < 1.e-6) then

         do k = 2,nza
            aru(k,iu) = 0.
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

            if (zm(k) <= hmin) then

               aru(k,iu) = 0.

            elseif (zm(km) >= hmax) then

               aru(k,iu) = dnv(iu) * dzt(k)

            elseif(zm(k) < hmax .and. zm(km) < hmin) then

               aru(k,iu) = dnv(iu) * .5 * (zm(k) - hmin)**2 / (hmax - hmin)
                              
               if (ipass < npass .and. aru(k,iu) < .01 * dnv(iu) * dzt(k)) then
                  topm(im1) = max(topm(im1),zm(k) + .003)
                  topm(im2) = max(topm(im2),zm(k) + .003)
                  aru(k,iu) = 0.
               endif
               
            elseif(zm(k) <  hmax .and. zm(km) >=  hmin) then

               aru(k,iu) = dnv(iu) * dzt(k)  &
                  * (.5 * (zm(k) + zm(km)) - hmin) / (hmax - hmin)

            elseif(zm(k) >= hmax .and. zm(km) < hmin) then

               aru(k,iu) = dnv(iu) * (zm(k) - .5 * (hmax + hmin))
               
               if (ipass < npass .and. aru(k,iu) < .01 * dnv(iu) * dzt(k)) then
                  topm(im1) = max(topm(im1),zm(k))
                  topm(im2) = max(topm(im2),zm(k))
                  aru(k,iu) = 0.
               endif
               
            elseif(zm(k) >= hmax .and. zm(km) >=  hmin) then

               aru(k,iu) = dnv(iu)  &
                  * (dzt(k) - .5 * (hmax - zm(km)) ** 2 / (hmax - hmin))

            else

               write(io6,*) 'aru option not reached ',k,j,  &
                  zm(k),zm(km),hmax,hmin
               stop 'stop aru defn'   

            endif

            if (aru(k,iu) > .005 * dnv(iu) * dzt(k)) lpu(iu) = k

         enddo

      endif

! Expand ARU with height in case of spherical geometry

      do k = 2,nza
         aru(k,iu) = aru(k,iu) * zfact(k)
      enddo

   enddo
   call rsub('U',1)

enddo  ! end ipass loop  
      
! Set ARU to zero at non-topo walls
! [ONLY FOR ISFCL = 0]

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(3)%jend(1); iu = jtab_u(3)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = 2,nza
      aru(k,iu) = 0.
   enddo
   lpu(iu) = nza 
      
enddo
call rsub('U',3)

! Expand ARW and VOLT with height for spherical geometry, and compute VOLTI

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(3)%jend(1); iw = jtab_w(3)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Number of vertices of IW polygon

   npoly = itab_w(iw)%npoly

! Loop over vertical levels

   do k = 2,nza

      if (volt(k,iw) > 1.e-9) then
         arw (k,iw) = arw (k,iw) * zfacm(k)**2
         volt(k,iw) = volt(k,iw) * zfact(k)**2
      else
         arw (k,iw) = 0.
         volt(k,iw) = 1.e-9
      endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!go to 1
! Option for stability: expand volt if too small relative to any grid cell face
 
      volt(k,iw) = max(volt(k,iw),real(arw(k,iw),8) * dzt(k))

! Loop over faces of IW polygon

      do ju = 1,npoly
         iu = itab_w(iw)%iu(ju)
         
         u_height = aru(k,iu) / dnv(iu)
         
         volt(k,iw) = max(volt(k,iw), real(arw0(iw),8) * zfact(k) * u_height)

      enddo

1 continue

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      volti(k,iw) = 1. / volt(k,iw)

   enddo  ! k

! Set arw = 0 for bottom (k = 1), wall-on-top (k = nza-1), 
! and top (k = nza) levels

   arw(1,iw) = 0.   
   arw(nza-1,iw) = 0.   
   arw(nza,iw) = 0.   

enddo
call rsub('W',3)

! LPW and LSW, and possible ARW, VOLT, and VOLTI modifications

nsw_max = 1

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(6)%jend(1); iw = jtab_w(6)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Initialize LPW(IW) and LSW(IW)

   lpw(iw) = nza-1
   lsw(iw) = 1

! Number of vertices of IW polygon

   npoly = itab_w(iw)%npoly

! Loop over vertical levels from top to bottom

   do k = nza-1,2,-1

! Loop over faces of IW polygon

      do ju = 1,npoly
         iu = itab_w(iw)%iu(ju)
         
         if (lpu(iu) <= k) lpw(iw) = k
      enddo

      if (lpw(iw) > k) then
         arw  (k,iw) = 0.     ! close area if all surrounding U's are closed
         volt (k,iw) = 1.e-9  ! close volume if all surrounding U's are closed
         volti(k,iw) = 1.e9   ! close volume if all surrounding U's are closed
      endif

      if (k < nza-1 .and. lpw(iw) <= k .and. arw(k,iw) == 0.) then
         print*, 'Subroutine ctrlvols lpw-arw mismatch '
         print*, 'k,iw = ',k,iw
         stop 'stop_1 ctrlvols '
      endif

! Increase LSW if K-1 W level intersects topography in this cell

      if (arw(k,iw) > 0. .and. arw(k,iw) < .999 * arw0(iw) * zfacm(k)**2) then
         lsw(iw) = lsw(iw) + 1
      endif

   enddo  ! k
   
! Increase nsw_max if necessary

   if (lsw(iw) > nsw_max) nsw_max = lsw(iw)

enddo
call rsub('W',6)

! VOLUI from VOLT

lcu(1:nua) = nza

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_u(4)%jend(1); iu = jtab_u(4)%iu(j)
   iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = nza,2,-1
      volui(k,iu) = 1. / (volt(k,iw1) + volt(k,iw2))      

      if (volt(k,iw1) + volt(k,iw2) > 1.e-8) then
         lcu(iu) = k
      endif
   enddo

   volui(1,iu) = 1.e9

enddo
call rsub('U',4)

! VOLWI from VOLT

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(5)%jend(1); iw = jtab_w(5)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 2,nza
      kp = min(k+1,nza)
      volwi(k,iw) = 2. / (volt(k,iw) + volt(kp,iw)) 
   enddo

! modify volwi for lpw and lpw-1 levels

   ka = lpw(iw)
   volwi(ka,iw)   = 1. / (volt(ka,iw) + .5 * volt(ka+1,iw))
   volwi(ka-1,iw) = 1. / (.5 * volt(ka,iw))

enddo
call rsub('W',5)

! In case ARW has been reset to 0 anywhere (because it was nearly zero), 
! transfer the seaflux and landflux cell values to KW = LPW(IW).

if (isfcl == 1) then

   do isf = 2,nseaflux
      iw = seaflux(isf)%iw
      kw = seaflux(isf)%kw

      if (kw < lpw(iw)) seaflux(isf)%kw = lpw(iw)
      if (kw > lpw(iw) + lsw(iw) - 1) seaflux(isf)%kw = lpw(iw) + lsw(iw) - 1
   enddo

   do ilf = 2,nlandflux
      iw = landflux(ilf)%iw
      kw = landflux(ilf)%kw

      if (kw < lpw(iw)) landflux(ilf)%kw = lpw(iw)
      if (kw > lpw(iw) + lsw(iw) - 1) landflux(ilf)%kw = lpw(iw) + lsw(iw) - 1
   enddo

endif

return
end subroutine ctrlvols_tri


