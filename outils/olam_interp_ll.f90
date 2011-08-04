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
subroutine interp_hvn_ll(nlon,nlat,nlevin,nlevout,alon,alat,field,field_ll)

use mem_grid,   only: mza, mma, mva, mwa, zm, zt, lpu, lpv, lpw, &
                      xem, yem, zem, xeu, yeu, zeu, &
                      xev, yev, zev, xew, yew, zew, &
                      glatm, glonm, glatw, glonw
use mem_ijtabs, only: itab_m, itab_w, itab_u, itab_v
use misc_coms,  only: io6, meshtype
use consts_coms, only: erad,piu180,pio180

implicit none

integer, intent(in) :: nlon,nlat,nlevin,nlevout
real, intent(in) :: alon(nlon),alat(nlat)
real, intent(in) :: field(nlevin,mva)
real, intent(out) :: field_ll(nlon,nlat,nlevout)

integer :: k,im,iv,iw,npoly,j,jnext,ivnext,ilat,ilon,kll,lonflag
integer :: k1,k2,koff,lpuvmax,lpwmax

real :: x(3),y(3),z(3)
real :: xw(7),yw(7),xuv(7),yuv(7)

real :: a(nlevin),b(nlevin),c(nlevin)
real :: field_avg(nlevin)

real :: raxis,radmax,rad
real :: deglat,abslat,aminlat,amaxlat,aminlon,amaxlon,deglon
real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
real :: qlat,qlon,qx,qy

if (nlevin == 1) then
   k1 = 1
   k2 = 1
   koff = 0
else
   k2 = mza-1
   koff = 1
endif

x(1) = 0.
y(1) = 0.

!-----------------------------------------------------------
! First loop is over M points for interpolating U/V points
!-----------------------------------------------------------

do im = 2,mma

   npoly = itab_m(im)%npoly

! Initialize maximum radius, maximum lpu/lpv, and field average

   radmax = 0.
   lpuvmax = 2
   field_avg(1:nlevin) = 0.

! Loop over all U or V points that surround current M point

   do j = 1,npoly

! Current U/V point index   

      if (meshtype == 1) then
         iv = itab_m(im)%iu(j)
      else
         iv = itab_m(im)%iv(j)
      endif

! Skip current U/V point if index < 2

      if (iv < 2) go to 8

! Transform current U/V point to PS coordinates tangent at M point

      if (meshtype == 1) then
         call e_ps(xeu(iv),yeu(iv),zeu(iv),glatm(im),glonm(im),xuv(j),yuv(j))
         if (lpuvmax < lpu(iv)) lpuvmax = lpu(iv)
      else
         call e_ps(xev(iv),yev(iv),zev(iv),glatm(im),glonm(im),xuv(j),yuv(j))
         if (lpuvmax < lpv(iv)) lpuvmax = lpv(iv)
      endif

      rad = sqrt(xuv(j)**2 + yuv(j)**2)
      if (radmax < rad) radmax = rad

      if (nlevin > 1) then
         k1 = lpuvmax
      endif

      do k = k1,k2
         field_avg(k) = field_avg(k) + field(k,iv) / real(npoly)
      enddo

   enddo

! Find maximum range of latitude and longitude inside the polygon of U points

   deglat = piu180 * radmax / erad

   aminlat = max(-90.,glatm(im) - deglat)
   amaxlat = min( 90.,glatm(im) + deglat)

   abslat = max(abs(aminlat),abs(amaxlat))

   lonflag = 0

   if (abslat > 89.9) then

      aminlon = -180.01
      amaxlon =  180.01

   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonm(im) - deglon
      amaxlon = glonm(im) + deglon

      if (aminlon < -180.) lonflag = 1
      if (amaxlon >  180.) lonflag = 2

   endif

! Loop over all U/V points that surround current M point and fill field values

   do j = 1,npoly

      jnext = j + 1
      if (j == npoly) jnext = 1

! Current U/V point index   

      if (meshtype == 1) then
         iv = itab_m(im)%iu(j)
         ivnext = itab_m(im)%iu(jnext)
      else
         iv = itab_m(im)%iv(j)
         ivnext = itab_m(im)%iv(jnext)
      endif

      x(2) = xuv(j)
      y(2) = yuv(j)

      x(3) = xuv(jnext)
      y(3) = yuv(jnext)

! Loop over vertical levels

      do k = k1,k2
         z(1) = field_avg(k)
         z(2) = field(k,iv)
         z(3) = field(k,ivnext)

! Evaluate interpolation coefficients for current trio of points

         call matrix_3x3(1.  , x(1), y(1),  &
                         1.  , x(2), y(2),  &
                         1.  , x(3), y(3),  &
                         z(1), z(2), z(3),  &
                         a(k), b(k), c(k)   )
      enddo

! Set up some triangle-check coefficients

      v0x = x(2) - x(1)
      v0y = y(2) - y(1)

      v1x = x(3) - x(1)
      v1y = y(3) - y(1) 

      dot00 = v0x * v0x + v0y * v0y
      dot01 = v0x * v1x + v0y * v1y
      dot11 = v1x * v1x + v1y * v1y

      denomi = 1. / (dot00 * dot11 - dot01 * dot01)

! Loop over all possible lat-lon points in range

      do ilat = 1,nlat

         if (alat(ilat) < aminlat .or. alat(ilat) > amaxlat) cycle

         do ilon = 1,nlon  ! loop over longitude points

            if (alon(ilon) < aminlon .and. &
               (lonflag /= 2 .or. alon(ilon) > amaxlon)) cycle

            if (alon(ilon) > amaxlon .and. &
               (lonflag /= 1 .or. alon(ilon) < aminlon)) cycle

! Transform current lat-lon point to PS space

            call ll_xy(alat(ilat),alon(ilon),glatm(im),glonm(im),qx,qy)         

! Set up remaining triangle_check coefficients

            v2x = qx - x(1)
            v2y = qy - y(1)

            dot02 = v0x * v2x + v0y * v2y
            dot12 = v1x * v2x + v1y * v2y

            u = (dot11 * dot02 - dot01 * dot12) * denomi
            v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

            if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

               do k = k1,k2
                  kll = k - koff

! Interpolate to current field point

                  field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy

               enddo  ! k

            endif  ! q point inside triangle

         enddo ! ilon

      enddo  ! ilat

   enddo   ! j

8  continue

enddo   ! im

!-----------------------------------------------------------
! Second loop is over W points for interpolating U/V points
!-----------------------------------------------------------

do iw = 2,mwa

   npoly = itab_w(iw)%npoly

! Initialize maximum radius, maximum lpu/lpv, and field average

   radmax = 0.
   lpuvmax = 2
   field_avg(1:nlevin) = 0.

! Loop over all U or V points that surround current W point

   do j = 1,npoly

! Current U/V point index   

      if (meshtype == 1) then
         iv = itab_w(iw)%iu(j)
      else
         iv = itab_w(iw)%iv(j)
      endif

! Skip current U/V point if ndex < 2

      if (iv < 2) go to 9

! Transform current U/V point to PS coordinates tangent at W point

      if (meshtype == 1) then
         call e_ps(xeu(iv),yeu(iv),zeu(iv),glatw(iw),glonw(iw),xuv(j),yuv(j))
         if (lpuvmax < lpu(iv)) lpuvmax = lpu(iv)
      else
         call e_ps(xev(iv),yev(iv),zev(iv),glatw(iw),glonw(iw),xuv(j),yuv(j))
         if (lpuvmax < lpv(iv)) lpuvmax = lpv(iv)
      endif

      rad = sqrt(xuv(j)**2 + yuv(j)**2)
      if (radmax < rad) radmax = rad

      if (nlevin > 1) then
         k1 = lpuvmax
      endif

      do k = k1,k2
         field_avg(k) = field_avg(k) + field(k,iv) / real(npoly)
      enddo

    enddo

! Find maximum range of latitude and longitude inside the polygon of W points

   deglat = piu180 * radmax / erad

   aminlat = max(-90.,glatw(iw) - deglat)
   amaxlat = min( 90.,glatw(iw) + deglat)

   abslat = max(abs(aminlat),abs(amaxlat))

   lonflag = 0

   if (abslat > 89.9) then

      aminlon = -180.01
      amaxlon =  180.01

   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonw(iw) - deglon
      amaxlon = glonw(iw) + deglon

      if (aminlon < -180.) lonflag = 1
      if (amaxlon >  180.) lonflag = 2

   endif

! Loop over all U/V points that surround current W point and fill field values

   do j = 1,npoly

      jnext = j + 1
      if (j == npoly) jnext = 1

! Current U/V point index   

      if (meshtype == 1) then
         iv = itab_w(iw)%iu(j)
         ivnext = itab_w(iw)%iu(jnext)
      else
         iv = itab_w(iw)%iv(j)
         ivnext = itab_w(iw)%iv(jnext)
      endif

      x(2) = xuv(j)
      y(2) = yuv(j)

      x(3) = xuv(jnext)
      y(3) = yuv(jnext)

! Loop over vertical levels

      do k = k1,k2
         z(1) = field_avg(k)
         z(2) = field(k,iv)
         z(3) = field(k,ivnext)

! Evaluate interpolation coefficients for current trio of points

         call matrix_3x3(1.  , x(1), y(1),  &
                         1.  , x(2), y(2),  &
                         1.  , x(3), y(3),  &
                         z(1), z(2), z(3),  &
                         a(k), b(k), c(k)   )
      enddo

! Set up some triangle-check coefficients

      v0x = x(2) - x(1)
      v0y = y(2) - y(1)

      v1x = x(3) - x(1)
      v1y = y(3) - y(1) 

      dot00 = v0x * v0x + v0y * v0y
      dot01 = v0x * v1x + v0y * v1y
      dot11 = v1x * v1x + v1y * v1y

      denomi = 1. / (dot00 * dot11 - dot01 * dot01)

! Loop over all possible lat-lon points in range

      do ilat = 1,nlat

         if (alat(ilat) < aminlat .or. alat(ilat) > amaxlat) cycle

         do ilon = 1,nlon  ! loop over longitude points

            if (alon(ilon) < aminlon .and. &
               (lonflag /= 2 .or. alon(ilon) > amaxlon)) cycle

            if (alon(ilon) > amaxlon .and. &
               (lonflag /= 1 .or. alon(ilon) < aminlon)) cycle

! Transform current lat-lon point to PS space

            call ll_xy(alat(ilat),alon(ilon),glatw(iw),glonw(iw),qx,qy)     

! Set up remaining triangle_check coefficients

            v2x = qx - x(1)
            v2y = qy - y(1)

            dot02 = v0x * v2x + v0y * v2y
            dot12 = v1x * v2x + v1y * v2y

            u = (dot11 * dot02 - dot01 * dot12) * denomi
            v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

            if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

               do k = k1,k2
                  kll = k - koff

! Interpolate to current field point

                  field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy

               enddo  ! k

            endif  ! q point inside triangle

         enddo ! ilon

      enddo  ! ilat

   enddo  ! j

9  continue

enddo   ! iw

return
end subroutine interp_hvn_ll

!================================================================================

subroutine interp_htw_ll(nlon,nlat,nlevin,nlevout,alon,alat,field,field_ll)

use mem_grid,   only: mma, mza, mwa, zm, zt, lpw, &
                      xem, yem, zem, xew, yew, zew, glatm, glonm

use mem_ijtabs, only: itab_m
use misc_coms,  only: io6
use consts_coms, only: erad,piu180,pio180

implicit none

integer, intent(in) :: nlon,nlat,nlevin,nlevout
real, intent(in) :: alon(nlon),alat(nlat)
real, intent(in) :: field(nlevin,mwa)
real, intent(inout) :: field_ll(nlon,nlat,nlevout)

integer :: k,im,iw,npoly,lpwmax,j,jnext,iwnext,ilat,ilon,kll,lonflag
integer :: k1,k2,koff

real :: x(3),y(3),z(3)
real :: xw(8),yw(8)

real :: a(nlevin),b(nlevin),c(nlevin)
real :: field_avg(nlevin)

real :: raxis,radmax,rad
real :: deglat,abslat,aminlat,amaxlat,aminlon,amaxlon,deglon
real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
real :: beglon0,qlat,qlon,qx,qy

if (nlevin == 1) then
   k1 = 1
   k2 = 1
   koff = 0
else
   k2 = mza-1
   koff = 1
endif

x(1) = 0.
y(1) = 0.

! Loop over M points for interpolating W points

do im = 2,mma

   npoly = itab_m(im)%npoly

! Initialize maximum radius, maximum lpw, and field average

   radmax = 0.
   lpwmax = 2
   field_avg(1:nlevin) = 0.

! Loop over all W points that surround current M point

   do j = 1,npoly

! Current W point index

      iw = itab_m(im)%iw(j)

! Skip this IM point if iw < 2

      if (iw < 2) go to 9

! Transform current W point to PS coordinates tangent at M point

      call e_ps(xew(iw),yew(iw),zew(iw),glatm(im),glonm(im),xw(j),yw(j))

      rad = sqrt(xw(j)**2 + yw(j)**2)
      if (radmax < rad)     radmax = rad
      if (lpwmax < lpw(iw)) lpwmax = lpw(iw)

      if (nlevin > 1) then
         k1 = lpwmax
      endif

      do k = k1,k2
         field_avg(k) = field_avg(k) + field(k,iw) / real(npoly)
      enddo

   enddo

! Find maximum range of latitude and longitude inside the polygon of W points

   deglat = piu180 * radmax / erad

   aminlat = max(-90.,glatm(im) - deglat)
   amaxlat = min( 90.,glatm(im) + deglat)

   abslat = max(abs(aminlat),abs(amaxlat))

   lonflag = 0

   if (abslat > 89.9) then

      aminlon = -180.01
      amaxlon =  180.01

   else

      deglon = min(180.,deglat / cos(abslat * pio180))
      aminlon = glonm(im) - deglon
      amaxlon = glonm(im) + deglon
      
      if (aminlon < -180.) lonflag = 1
      if (amaxlon >  180.) lonflag = 2

   endif

! Loop over all W points that surround current M point and fill field values

   do j = 1,npoly
      jnext = j + 1
      if (j == npoly) jnext = 1

      iw     = itab_m(im)%iw(j)
      iwnext = itab_m(im)%iw(jnext)

      x(2) = xw(j)
      y(2) = yw(j)

      x(3) = xw(jnext)
      y(3) = yw(jnext)

! Loop over vertical levels

      do k = k1,k2
         z(1) = field_avg(k)
         z(2) = field(k,iw)
         z(3) = field(k,iwnext)

! Evaluate interpolation coefficients for current trio of points

         call matrix_3x3(1.  , x(1), y(1),  &
                         1.  , x(2), y(2),  &
                         1.  , x(3), y(3),  &
                         z(1), z(2), z(3),  &
                         a(k), b(k), c(k)   )
      enddo

! Set up some triangle-check coefficients

      v0x = x(2) - x(1)
      v0y = y(2) - y(1)

      v1x = x(3) - x(1)
      v1y = y(3) - y(1)

      dot00 = v0x * v0x + v0y * v0y
      dot01 = v0x * v1x + v0y * v1y
      dot11 = v1x * v1x + v1y * v1y

      denomi = 1. / (dot00 * dot11 - dot01 * dot01)

! Loop over all possible lat-lon points in range

      do ilat = 1,nlat

         if (alat(ilat) < aminlat .or. alat(ilat) > amaxlat) cycle

         do ilon = 1,nlon  ! loop over longitude points

            if (alon(ilon) < aminlon .and. &
               (lonflag /= 2 .or. alon(ilon) > amaxlon)) cycle

            if (alon(ilon) > amaxlon .and. &
               (lonflag /= 1 .or. alon(ilon) < aminlon)) cycle

! Transform current lat-lon point to PS space

            call ll_xy(alat(ilat),alon(ilon),glatm(im),glonm(im),qx,qy)

! Set up remaining triangle_check coefficients

            v2x = qx - x(1)
            v2y = qy - y(1)

            dot02 = v0x * v2x + v0y * v2y
            dot12 = v1x * v2x + v1y * v2y

            u = (dot11 * dot02 - dot01 * dot12) * denomi
            v = (dot00 * dot12 - dot01 * dot02) * denomi

! Check if current qx,qy point is inside or very near current triangle

            if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

! Point is inside or very near triangle; loop over vertical levels

               do k = k1,k2
                  kll = k - koff

! Interpolate to current field point

                  field_ll(ilon,ilat,kll) = a(k) + b(k) * qx + c(k) * qy

               enddo  ! k

            endif  ! q point inside triangle

         enddo ! ilon

      enddo  ! ilat

   enddo   ! j

9  continue

enddo   ! im

return
end subroutine interp_htw_ll


