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
subroutine icosahedron(nxp0)

use mem_ijtabs,  only: mrls, itab_md, itab_ud, itab_wd, alloc_itabsd
use mem_grid,    only: nza, nma, nua, nwa, xem, yem, zem, &
                       alloc_xyzem, impent
use misc_coms,   only: io6
use consts_coms, only: pi2, erad, erador5

implicit none

integer, intent(in) :: nxp0

real, parameter :: pwrd = 1.0  ! 0.9 is close to making uniform-sized triangles
                               ! 1.0 is original value

integer :: ibigd,i,j,idiamond,im_left,iu0,iu1,iu2,iu3,iu4,iw1,iw2,im &
   ,idiamond_top,im_top,im_right,im_bot,nn10,idiamond_right,idiamond_bot &
   ,iu,iw
integer :: id

real :: wts,wtn,wtw,wte,expansion,anglen,anglew,anglee,angles,wtw0,wte0,sumwt

integer, save, dimension(10) :: ibigd_ne = (/6,7,8,9,10,7,8,9,10,6/) &
                               ,ibigd_se = (/2,3,4,5,1,2,3,4,5,1/)

real, save, dimension(10) :: xed_s,xed_n,xed_w,xed_e &
                            ,yed_s,yed_n,yed_w,yed_e &
                            ,zed_s,zed_n,zed_w,zed_e

! Define triangles, edges, and vertices for icosahedral faces and subdivisions

mrls = 1  ! Default value

! For now, use nxp0 to divide each face

nn10 = nxp0 * nxp0 * 10

! ADD 1 to total number of points needed

nma =     nn10 + 2 + 1  ! ADDING 1 for reference point (index = 1)
nua = 3 * nn10     + 1  ! ADDING 1 for reference point (index = 1)
nwa = 2 * nn10     + 1  ! ADDING 1 for reference point (index = 1)

! Allocate memory for itabs and M earth coords
! Initialize all neighbor indices to zero

call alloc_itabsd(nma,nua,nwa)

call alloc_xyzem(nma)

do im = 2,nma
   itab_md(im)%itopm = im
   call mdloops('f',im,1,0,1,0)
enddo

do iu = 2,nua
   itab_ud(iu)%iup = iu
   call udloops('f',iu, 1, 4, 7, 8,11,12,13,14,16,20)
   call udloops('n',iu,21,22,23, 0, 0, 0, 0, 0, 0, 0)
enddo

do iw = 2,nwa
   itab_wd(iw)%iwp = iw
   call wdloops('f',iw, 1, 3, 5, 6, 7, 8,11,12,13,14)
   call wdloops('n',iw,15,16,17,18,19,20,21,23,25,26)
   call wdloops('n',iw,27,28,29,30,33,34, 0, 0, 0, 0)
enddo

! Fill big diamond corner coordinates

do id = 1,5

   anglen = .2 * (id-1) * pi2
   anglew = anglen - .1 * pi2
   anglee = anglen + .1 * pi2

   zed_s(id) = -erad
   xed_s(id) = 0.
   yed_s(id) = 0.

   zed_n(id) = erador5
   xed_n(id) = erador5 * 2. * cos(anglen)
   yed_n(id) = erador5 * 2. * sin(anglen)

   zed_w(id) = -erador5
   xed_w(id) = erador5 * 2. * cos(anglew)
   yed_w(id) = erador5 * 2. * sin(anglew)

   zed_e(id) = -erador5
   xed_e(id) = erador5 * 2. * cos(anglee)
   yed_e(id) = erador5 * 2. * sin(anglee)

enddo

do id = 6,10

   angles = .2 * (id-6) * pi2 + .1 * pi2
   anglew = angles - .1 * pi2
   anglee = angles + .1 * pi2

   zed_s(id) = -erador5
   xed_s(id) = erador5 * 2. * cos(angles)
   yed_s(id) = erador5 * 2. * sin(angles)

   zed_n(id) = erad
   xed_n(id) = 0.
   yed_n(id) = 0.

   zed_w(id) = erador5
   xed_w(id) = erador5 * 2. * cos(anglew)
   yed_w(id) = erador5 * 2. * sin(anglew)

   zed_e(id) = erador5
   xed_e(id) = erador5 * 2. * cos(anglee)
   yed_e(id) = erador5 * 2. * sin(anglee)

enddo

! Store IM index of south-pole and north-pole pentagonal points

impent(1) = 2
impent(12) = nma

do ibigd = 1,10

   do j = 1,nxp0
      do i = 1,nxp0

         idiamond = (ibigd - 1) * nxp0 * nxp0 &
                  + (j - 1)     * nxp0        &
                  +  i

! Indices that are "attached" to this diamond

         im_left  = idiamond + 2

         iu0 = 3 * idiamond
         iu1 = 3 * idiamond - 1
         iu3 = 3 * idiamond + 1

         iw1 = 2 * idiamond
         iw2 = 2 * idiamond + 1

! Store IM index of 10 out of 12 pentagonal points

         if (i == 1 .and. j == nxp0) impent(ibigd+1) = im_left

! Indices that are "attached" to another diamond

         if (ibigd < 6) then   ! Southern 5 diamonds

            ! Top diamond indices

            if (i < nxp0) then
               idiamond_top = idiamond + 1
            else
               idiamond_top = (ibigd_ne(ibigd) - 1) * nxp0 * nxp0 &
                            + (j - 1)               * nxp0        &
                            +  1
            endif

            im_top   = idiamond_top + 2
            iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)

            ! Right diamond indices

            if (j > 1 .and. i < nxp0) then
               idiamond_right = idiamond - nxp0 + 1
            elseif (j == 1) then
               idiamond_right = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                              + (i - 1)             * nxp0        &
                              +  1
               iu2 = 3 * idiamond_right - 1 ! (it's the iu1 for id_right)
            else            ! case for i = nxp0 and j > 1
               idiamond_right = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                              + (j - 2)             * nxp0        &
                              +  1
            endif

            im_right = idiamond_right + 2

            ! Bottom diamond indices
            
            if (j > 1) then
               idiamond_bot = idiamond - nxp0
               iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)
            else
               idiamond_bot = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                            + (i - 2)             * nxp0        &
                            +  1
            endif
            im_bot = idiamond_bot + 2
            
            if (i == 1 .and. j == 1) im_bot = 2

         else                  ! Northern 5 diamonds

           ! Top diamond indices

            if (i < nxp0) then
               idiamond_top = idiamond + 1
               iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)
            else
               idiamond_top = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                            + (nxp0 - 1)          * nxp0        &
                            +  j + 1
            endif

            im_top = idiamond_top + 2

            ! Right diamond indices

            if (j > 1 .and. i < nxp0) then
               idiamond_right = idiamond - nxp0 + 1
            elseif (j == 1 .and. i < nxp0) then
               idiamond_right = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                              + (nxp0 - 1)          * nxp0        &
                              + i + 1
            else            ! case for i = nxp0
               idiamond_right = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                              + (nxp0 - 1)          * nxp0        &
                              +  j
               iu4 = 3 * idiamond_right + 1 ! (it's the iu3 for id_right)
            endif

            im_right = idiamond_right + 2

            ! Bottom diamond indices
            
            if (j > 1) then
               idiamond_bot = idiamond - nxp0
            else
               idiamond_bot = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                            + (nxp0 - 1)          * nxp0        &
                            + i 
            endif

            im_bot = idiamond_bot + 2
            iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)
            
            if (i == nxp0 .and. j == nxp0) &
               im_top = 10 * nxp0 * nxp0 + 3

         endif

         call fill_diamond(im_left,im_right,im_top,im_bot &
            ,iu0,iu1,iu2,iu3,iu4,iw1,iw2)
            
! M point (xem,yem,zem) coordinates
         
         if (i + j <= nxp0) then
            wts  = max(0.,min(1.,real(nxp0 + 1 - i - j) / real(nxp0)))
            wtn  = 0.
            wtw0 = max(0.,min(1.,real(j) / real(i + j - 1)))
            wte0 = 1. - wtw0
         else
            wts  = 0.
            wtn  = max(0.,min(1.,real(i + j - nxp0 - 1) / real(nxp0)))
            wte0 = max(0.,min(1.,real(nxp0 - j) &
                 / real(2 * nxp0 + 1 - i - j)))
            wtw0 = 1. - wte0
         endif

! Experimental adjustment in spacing
! Compute sum of weights raised to pwrd

         wtw = (1. - wts - wtn) * wtw0
         wte = (1. - wts - wtn) * wte0
         sumwt = wts**pwrd + wtn**pwrd + wtw**pwrd + wte**pwrd

         wts = wts**pwrd / sumwt
         wtn = wtn**pwrd / sumwt
         wtw = wtw**pwrd / sumwt
         wte = wte**pwrd / sumwt

         xem(im_left) = wts * xed_s(ibigd) &
                      + wtn * xed_n(ibigd) &
                      + wtw * xed_w(ibigd) &
                      + wte * xed_e(ibigd)

         yem(im_left) = wts * yed_s(ibigd) &
                      + wtn * yed_n(ibigd) &
                      + wtw * yed_w(ibigd) &
                      + wte * yed_e(ibigd)

         zem(im_left) = wts * zed_s(ibigd) &
                      + wtn * zed_n(ibigd) &
                      + wtw * zed_w(ibigd) &
                      + wte * zed_e(ibigd)

! Push M point coordinates out to earth radius

         expansion = erad / sqrt(xem(im_left) ** 2 &
                               + yem(im_left) ** 2 &
                               + zem(im_left) ** 2 )


         xem(im_left) = xem(im_left) * expansion
         yem(im_left) = yem(im_left) * expansion
         zem(im_left) = zem(im_left) * expansion
         
      enddo  ! end i loop
   enddo     ! end j loop

enddo        ! end idbig loop

xem(2) = 0.
yem(2) = 0.
zem(2) = -erad

xem(nma) = 0.
yem(nma) = 0.
zem(nma) = erad

!call twist()

call tri_neighbors()

! This is the place to do spring dynamics

call spring_dynamics(1)

return
end subroutine icosahedron

!===============================================================================

subroutine fill_diamond(im_left,im_right,im_top,im_bot  &
   ,iu0,iu1,iu2,iu3,iu4,iw1,iw2)

use mem_ijtabs, only: itab_ud, itab_wd
use misc_coms,  only: io6

implicit none

integer, intent(in) :: im_left,im_right,im_top,im_bot
integer, intent(in) :: iu0,iu1,iu2,iu3,iu4,iw1,iw2     

itab_ud(iu0)%im(1) = im_left
itab_ud(iu0)%im(2) = im_right
itab_ud(iu0)%iw(1) = iw1
itab_ud(iu0)%iw(2) = iw2
itab_ud(iu0)%mrlu = 1
         
itab_ud(iu1)%im(1) = im_left
itab_ud(iu1)%im(2) = im_bot
itab_ud(iu1)%iw(2) = iw1
itab_ud(iu1)%mrlu = 1
         
itab_ud(iu2)%iw(1) = iw1

itab_ud(iu3)%im(1) = im_top
itab_ud(iu3)%im(2) = im_left
itab_ud(iu3)%iw(2) = iw2
itab_ud(iu3)%mrlu = 1

itab_ud(iu4)%iw(1) = iw2
         
itab_wd(iw1)%iu(1) = iu0 
itab_wd(iw1)%iu(2) = iu1 
itab_wd(iw1)%iu(3) = iu2 
itab_wd(iw1)%mrlw = 1
itab_wd(iw1)%mrlw_orig = 1
      
itab_wd(iw2)%iu(1) = iu0 
itab_wd(iw2)%iu(2) = iu4 
itab_wd(iw2)%iu(3) = iu3 
itab_wd(iw2)%mrlw = 1
itab_wd(iw2)%mrlw_orig = 1

return
end subroutine fill_diamond
   
!===============================================================================

subroutine twist()

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd
use mem_grid,    only: nma, nua, xem, yem, zem
use misc_coms,   only: io6, nxp
use consts_coms, only: pi1, pi2

implicit none

integer :: im,im1,im2,im1_east,im2_east
integer :: iu,iu1,iu2,iu3,iu_east,iueq,iueq_east,iueq_fill
integer :: iw,iw1

integer :: nxpo2  ! half of nxp (nxp must be even if twist is done)

integer :: iueq_tab(10*nxp)
integer :: iweq_tab(10*nxp)

real :: radm  ! Distance of M point from Earth axis [m]
real :: angm  ! longitude (radians) of M point before shift is applied

nxpo2 = nxp / 2

! Initialize iueq_tab and iweq_tab to zero

iueq_tab(:) = 0
iweq_tab(:) = 0

! Shift all M points in southern hemisphere (excluding equator) 
! 36 degrees to the east

do im = 2,nma
   if (zem(im) < -100.) then

      radm = sqrt(xem(im) ** 2 + yem(im) ** 2)
      angm = atan2(yem(im),xem(im))

      xem(im) = radm * cos(angm + .2 * pi1)
      yem(im) = radm * sin(angm + .2 * pi1)

   endif
enddo

! Loop over all U points and search for those that are on the equator

do iu = 1,nua
   im1 = itab_ud(iu)%im(1)
   im2 = itab_ud(iu)%im(2)

   if (abs(zem(im1)) < 100. .and. abs(zem(im2)) < 100.) then

! Compute special index for equatorial U points that increases with longitude.
! This uses knowledge of icosahedron subroutine that im2 is always east of iu.
! IUEQ skips approximately every other integer and therefore needs to be collapsed.

      iueq = int(10 * nxp * (atan2(yem(im2),xem(im2)) + pi1) / pi2) + 1   

! Fill table of iu and iw1 indices for each iueq

      iueq_tab(iueq) = iu
      iweq_tab(iueq) = itab_ud(iu)%iw(1)

   endif
enddo

! Collapse IUEQ tables in order to use consecutive integer indices

iueq_fill = 0

do iueq = 1,nxp * 10

   if (iueq_tab(iueq) > 1) then
      iueq_fill = iueq_fill + 1

      iueq_tab(iueq_fill) = iueq_tab(iueq)
      iweq_tab(iueq_fill) = iweq_tab(iueq)
   endif

enddo

! Loop over equatorial U points

do iueq = 1,nxp * 5

! Find equatorial index of U point that is 36 degrees to the east

   iueq_east = iueq + nxpo2
   if (iueq_east > nxp * 5) iueq_east = iueq_east - nxp * 5

! Get U point indices for iueq and iueq_east points

   iu      = iueq_tab(iueq)
   iu_east = iueq_tab(iueq_east)

! Get index for adjacent W point to the south of iu

   iw = iweq_tab(iueq)

! Get indices for M endpoints of iu and iu_east

   im1 = itab_ud(iu)%im(1)
   im2 = itab_ud(iu)%im(2)

   im1_east = itab_ud(iu_east)%im(1)
   im2_east = itab_ud(iu_east)%im(2)

! Get U neighbor indices for adjacent W point to the south

   iu1 = itab_wd(iw)%iu(1)
   iu2 = itab_wd(iw)%iu(2)
   iu3 = itab_wd(iw)%iu(3)

! For U points on equator: shift their south W point neighbor index 

   itab_ud(iu_east)%iw(1) = iw

! For W points bordering equator on south: shift their equatorial U point
! neighbor index

   if     (iu1 == iu) then
      itab_wd(iw)%iu(1) = iu_east
   elseif (iu2 == iu) then
      itab_wd(iw)%iu(2) = iu_east
   elseif (iu3 == iu) then
      itab_wd(iw)%iu(3) = iu_east
   endif

! For U points bordering equator on south: shift their equatorial M point 
! neighbor index

   if (iu1 /= iu) then
      if     (itab_ud(iu1)%im(1) == im1) then
              itab_ud(iu1)%im(1) =  im1_east
      elseif (itab_ud(iu1)%im(2) == im1) then
              itab_ud(iu1)%im(2) =  im1_east
      elseif (itab_ud(iu1)%im(1) == im2) then
              itab_ud(iu1)%im(1) =  im2_east
      elseif (itab_ud(iu1)%im(2) == im2) then
              itab_ud(iu1)%im(2) =  im2_east
      endif
   endif

   if (iu2 /= iu) then
      if     (itab_ud(iu2)%im(1) == im1) then
              itab_ud(iu2)%im(1) =  im1_east
      elseif (itab_ud(iu2)%im(2) == im1) then
              itab_ud(iu2)%im(2) =  im1_east
      elseif (itab_ud(iu2)%im(1) == im2) then
              itab_ud(iu2)%im(1) =  im2_east
      elseif (itab_ud(iu2)%im(2) == im2) then
              itab_ud(iu2)%im(2) =  im2_east
      endif
   endif

   if (iu3 /= iu) then
      if     (itab_ud(iu3)%im(1) == im1) then
              itab_ud(iu3)%im(1) =  im1_east
      elseif (itab_ud(iu3)%im(2) == im1) then
              itab_ud(iu3)%im(2) =  im1_east
      elseif (itab_ud(iu3)%im(1) == im2) then
              itab_ud(iu3)%im(1) =  im2_east
      elseif (itab_ud(iu3)%im(2) == im2) then
              itab_ud(iu3)%im(2) =  im2_east
      endif
   endif

enddo

return
end subroutine twist

