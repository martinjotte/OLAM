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
subroutine contslab_horiz_mp(iplt)

use oplot_coms, only: op
use mem_grid,   only: mza, mma, mwa, zm, zt, lpw, xem, yem, zem, xew, yew, zew
use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: npoly,j,kt,k,iw,jw,im,notavail
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real :: hcpn(7),vcpn(7),fldvals(7)

integer :: ktf(mwa),km(mma)
real :: wtbot(mma),wttop(mma)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

! Loop over W points for contouring M points

do jw = 1, jtab_w(jtw_prog)%jend(1)
   iw = jtab_w(jtw_prog)%iw(jw)

! Skip this W point if it is underground

   if (ktf(iw) /= 0) cycle

! Get plot coordinates of current W point.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

   npoly = itab_w(iw)%npoly

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over all M points that surround current W point

   do j = 1,npoly

! Current M point index   

      im = itab_w(iw)%im(j)

! Skip current M point if index < 2

      if (im < 2) go to 9

! Get plot coordinates of current M point

      call oplot_transform(iplt,xem(im),yem(im),zem(im),hcpn(j),vcpn(j))

! Skip this W point if current M point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(j)) > 1.e11) go to 9
      
! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(j))
         if (hcpn(j) < -180.001) iflag180 =  1
         if (hcpn(j) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(j) >= op%xmin) ipwx1 = 1 
      if (hcpn(j) <= op%xmax) ipwx2 = 1 
      if (vcpn(j) >= op%ymin) ipwy1 = 1 
      if (vcpn(j) <= op%ymax) ipwy2 = 1 

   enddo

! If any window flag is zero, all M points for this W point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! Loop over all M points that surround current W point and fill field values

   do j = 1,npoly
      im = itab_w(iw)%im(j)
      call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                     fldvals(j),notavail)
      if (notavail > 0) go to 9 
   enddo

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

   if (iflag180 /= 0) then

      do j = 1,npoly
         hcpn(j) = hcpn(j) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)
            
   endif

9  continue
            
enddo

! NOTE: For Hex grid, may in the future want underground to always be
!       triangle areas.  In that case, first do underground plot and then
!       call contpolyg selectively for trios of points in remaining 
!       (open) sectors, using U/V point at face of underground block.

return
end subroutine contslab_horiz_mp

!===============================================================================

subroutine contslab_horiz_vn(iplt)

use oplot_coms, only: op
use mem_grid,   only: mma, mva, mwa, zm, zt, lpv, lpw, &
                      xem, yem, zem, &
                      xev, yev, zev, xew, yew, zew
use mem_ijtabs, only: itab_m, itab_w, itab_v, jtab_m, jtm_vadj, &
                      jtab_w, jtw_prog
use misc_coms,  only: io6, isubdomain
use mem_para,   only: myrank

implicit none

integer, intent(in) :: iplt

integer :: npoly,j,jm,jn,jnn,kt,k,im,notavail,iw,iv,iw1,iw2,iv1,iv2,jw
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real :: hcpn(7),vcpn(7),fldvals(7)
real :: hcpn3(3),vcpn3(3),fldvals3(3)
real :: avail, avg

integer :: ktf(mwa),kv(mva)
real :: wtbot(mva),wttop(mva)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

!------------------------------------------------------
! FIRST LOOP is over M points for contouring V points
!------------------------------------------------------

do jm = 1, jtab_m(jtm_vadj)%jend(1)
   im = jtab_m(jtm_vadj)%im(jm)

! Get plot coordinates of current M point.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   npoly = itab_m(im)%npoly

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

   avail = 0.
   avg = 0.

! Loop over all U or V points that surround current M point

   do j = 1,npoly

! Current U/V point index   

      iv = itab_m(im)%iv(j)
      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

! TEMPORARY FIX TO AVOID ACCESSING UNDEFINED VALUES IN PARALLEL

      if (isubdomain == 1) then
         if (itab_v(iv)%irank /= myrank) goto 8
      endif

! Skip current V cell if index < 2

      if (iv < 2) go to 8

! Get plot coordinates of current V point

      call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hcpn(j),vcpn(j))

! Skip this M point if current V point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(j)) > 1.e11) go to 8
      
! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(j))
         if (hcpn(j) < -180.001) iflag180 =  1
         if (hcpn(j) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any V point is on window side of 
! respective boundary

      if (hcpn(j) >= op%xmin) ipwx1 = 1 
      if (hcpn(j) <= op%xmax) ipwx2 = 1 
      if (vcpn(j) >= op%ymin) ipwy1 = 1 
      if (vcpn(j) <= op%ymax) ipwy2 = 1 

      if (ktf(iw1) /= 0 .and. ktf(iw2) /= 0) cycle

      call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                     fldvals(j),notavail)

      if (ktf(iw1) == 0) then
         avail = avail + .5
         avg = avg + .5 * fldvals(j)
      endif

      if (ktf(iw2) == 0) then
         avail = avail + .5
         avg = avg + .5 * fldvals(j)
      endif

   enddo

   if (avail < .1) cycle
   
   avg = avg / avail
   
! If any window flag is zero, all V points for this M point are outside
! the same window boundary, so skip this M point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 8

! Loop over all V points that surround current M point and contour plot each
! available sector 

   do j = 1,npoly
      jn = j + 1
      if (jn > npoly) jn = 1
      jnn = jn + 1
      if (jnn > npoly) jnn = 1

      iv1 = itab_m(im)%iv(j)
      iv2 = itab_m(im)%iv(jn)

! Specific way to get IW since ordering of W and U/V neighbors of M is not
! identical for both grid systems

      iw = itab_m(im)%iw(jnn)

      if (ktf(iw) == 0) then
      
         hcpn3(1) = hpt
         vcpn3(1) = vpt
         fldvals3(1) = avg

         hcpn3(2) = hcpn(j)
         vcpn3(2) = vcpn(j)
         fldvals3(2) = fldvals(j)
         
         hcpn3(3) = hcpn(jn)
         vcpn3(3) = vcpn(jn)
         fldvals3(3) = fldvals(jn)

         call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)
         
! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
! again at other end of plot

         if (iflag180 /= 0) then

            hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180
      
            call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)
            
         endif

      endif

   enddo

8  continue

enddo

!------------------------------------------------------
! SECOND LOOP is over W points for contouring V points
!------------------------------------------------------

do jw = 1, jtab_w(jtw_prog)%jend(1)
   iw = jtab_w(jtw_prog)%iw(jw)

   if (ktf(iw) /= 0) cycle

! Get plot coordinates of current W point.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

   npoly = itab_w(iw)%npoly

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

! Loop over all V points that surround current W point

   do j = 1,npoly

! Current V point index   

      iv = itab_w(iw)%iv(j)

! TEMPORARY FIX TO AVOID ACCESSING UNDEFINED VALUES IN PARALLEL

      if (isubdomain == 1) then
         if (itab_v(iv)%irank /= myrank) goto 9
      endif

! Skip current V cell if index < 2

      if (iv < 2) go to 9

! Get plot coordinates of current V point

      call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hcpn(j),vcpn(j))

! Skip this M point if current V point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(j)) > 1.e11) go to 9
      
! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(j))
         if (hcpn(j) < -180.001) iflag180 =  1
         if (hcpn(j) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any M point is on window side of 
! respective boundary

      if (hcpn(j) >= op%xmin) ipwx1 = 1 
      if (hcpn(j) <= op%xmax) ipwx2 = 1 
      if (vcpn(j) >= op%ymin) ipwy1 = 1 
      if (vcpn(j) <= op%ymax) ipwy2 = 1 

   enddo
      
! If any window flag is zero, all M points for this W point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! Loop over all V points that surround current W point and fill field values

   do j = 1,npoly
      iv = itab_w(iw)%iv(j)
      
      call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                     fldvals(j),notavail)
      if (notavail > 0) cycle 
   enddo

! Contour plot cell of 2-D or 3-D field

   call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

! If this polygon crosses +/- 180 degrees longitude in lat/lon plot, re-plot
! at other end

   if (iflag180 /= 0) then

      do j = 1,npoly
         hcpn(j) = hcpn(j) + 360. * iflag180
      enddo
      
      call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)
            
   endif

9  continue
            
enddo

return
end subroutine contslab_horiz_vn

!===============================================================================

subroutine contslab_horiz_tw(iplt)

use oplot_coms, only: op
use mem_grid,   only: mva, mma, mwa, zm, zt, lpw, xem, yem, zem, &
                      xev, yev, zev, xew, yew, zew
use mem_ijtabs, only: itab_m, jtab_m, jtm_vadj
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt

integer :: npoly,j,jm,jn,jnn,kt,k,im,iv,iw,iw1,iw2,notavail
integer :: iflag180
integer :: ipwx1,ipwx2,ipwy1,ipwy2

real :: hpt,vpt
real :: hcpn(7),vcpn(7),fldvals(7)
real :: hcpn3(3),vcpn3(3),fldvals3(3)
real :: avail, avg

integer :: ktf(mwa),kw(mwa)
real :: wtbot(mwa),wttop(mwa)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

! Loop over M points for contouring W points

do jm = 1, jtab_m(jtm_vadj)%jend(1)
   im = jtab_m(jtm_vadj)%im(jm)

! Get plot coordinates of current M point.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   npoly = itab_m(im)%npoly

! Initialize iflag180 and plot window flags to zero

   iflag180 = 0

   ipwx1 = 0
   ipwx2 = 0
   ipwy1 = 0
   ipwy2 = 0

   avail = 0.
   avg = 0.

! Loop over all W points that surround current M point

   do j = 1,npoly

! Current W point index   

      iw = itab_m(im)%iw(j)

! Skip current W cell if index < 2

      if (iw < 2) go to 9

! Get plot coordinates of current W point

      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hcpn(j),vcpn(j))

! Skip this M point if current W point is far outside plot window 
! (which means that orthographic projection returned large value that
! indicates that point is on other side of Earth)

      if (abs(hcpn(j)) > 1.e11) go to 9

! Avoid wrap-around and set iflag180

      if (op%projectn(iplt)== 'L') then
         call ll_unwrap(hpt,hcpn(j))
         if (hcpn(j) < -180.001) iflag180 =  1
         if (hcpn(j) >  180.001) iflag180 = -1
      endif

! Set plot window flag to 1 if any W point is on window side of 
! respective boundary

      if (hcpn(j) >= op%xmin) ipwx1 = 1 
      if (hcpn(j) <= op%xmax) ipwx2 = 1 
      if (vcpn(j) >= op%ymin) ipwy1 = 1 
      if (vcpn(j) <= op%ymax) ipwy2 = 1 

      if (ktf(iw) /= 0) cycle

      call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                     fldvals(j),notavail)

      avail = avail + 1.
      avg = avg + fldvals(j)

   enddo

   if (avail < .1) cycle
   
   avg = avg / avail
   
! If any window flag is zero, all W points for this M point are outside
! the same window boundary, so skip this W point

   if (ipwx1 == 0 .or. ipwx2 == 0 .or. ipwy1 == 0 .or. ipwy2 == 0) go to 9

! If all W points around this M point are available, plot them together

   if (nint(avail) == npoly) then

      call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)

! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
! again at other end of plot

      if (iflag180 /= 0) then

         hcpn(1:npoly) = hcpn(1:npoly) + 360. * iflag180
      
         call contpolyg(op%icolortab(iplt),op%ifill,npoly,hcpn,vcpn,fldvals)
            
      endif

! If only some W points are available, plot them by sectors

   else
            
! Loop over all W points that surround current M point and contour plot each
! available sector 

      do j = 1,npoly
         jn = j + 1
         if (jn > npoly) jn = 1
         jnn = jn + 1
         if (jnn > npoly) jnn = 1
      
         iw1 = itab_m(im)%iw(j)
         iw2 = itab_m(im)%iw(jn)

         hcpn3(1) = hpt
         vcpn3(1) = vpt
         fldvals3(1) = avg

         if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then
      
            hcpn3(2) = hcpn(j)
            vcpn3(2) = vcpn(j)
            fldvals3(2) = fldvals(j)
         
            hcpn3(3) = hcpn(jn)
            vcpn3(3) = vcpn(jn)
            fldvals3(3) = fldvals(jn)

         elseif (ktf(iw1) == 0) then
            
! Specific way to get IV since ordering of W and U/V neighbors of M is not
! identical for both grid systems

            iv = itab_m(im)%iv(jnn)
            call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hcpn3(3),vcpn3(3))

            fldvals3(3) = fldvals(j)

            hcpn3(2) = hcpn(j)
            vcpn3(2) = vcpn(j)
            fldvals3(2) = fldvals(j)

         elseif (ktf(iw2) == 0) then
         
! Specific way to get IV since ordering of W and U/V neighbors of M is not
! identical for both grid systems

            iv = itab_m(im)%iv(jnn)
            call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hcpn3(2),vcpn3(2))

            fldvals3(2) = fldvals(jn)

            hcpn3(3) = hcpn(jn)
            vcpn3(3) = vcpn(jn)
            fldvals3(3) = fldvals(jn)

         else
         
            cycle

         endif

         call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)

! If lat/lon plot and this polygon crosses +/- 180 degrees longitude, plot
! again at other end of plot

         if (iflag180 /= 0) then

            hcpn3(1:3) = hcpn3(1:3) + 360. * iflag180
      
            call contpolyg(op%icolortab(iplt),op%ifill,3,hcpn3,vcpn3,fldvals3)
            
         endif

      enddo

   endif

9  continue

enddo

return
end subroutine contslab_horiz_tw

!===============================================================================

subroutine contslab_vert_v(iplt)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, lpw, zt
use misc_coms,  only: io6
use mem_ijtabs, only: jtab_w, jtw_prog
use consts_coms,only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: k,jw,iw,iv1,iv2,iok,notavail
real :: hpt
real :: hcpn(4),vcpn(4),fldvals(4)
real :: topo1,topo2
real :: wtbot = 1., wttop = 0.

! First plot underground T cells with underground color

call plot_underground_w(iplt,(/0/))

do jw = 1, jtab_w(jtw_prog)%jend(1)
   iw = jtab_w(jtw_prog)%iw(jw)

! Get horizontal plot coordinates for IW point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,hcpn)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,hcpn) ! Need to fix for hex??
   endif

! Skip current IW point if it does not intersect plot cone

   if (iok /= 1) cycle

! Avoid wrap_around

   if (op%projectn(iplt) == 'C') then
      call ll_unwrap(hcpn(1),hcpn(2))
   endif

   hpt = .5 * (hcpn(1) + hcpn(2))

! Skip current IV point if either cell side is outside plot window. 

   if (hcpn(1) < op%xmin .or. hcpn(1) > op%xmax .or.  &
       hcpn(2) < op%xmin .or. hcpn(2) > op%xmax) cycle
   
   hcpn(3) = hcpn(2)
   hcpn(4) = hcpn(1)
   
   do k = lpw(iw),mza-1   ! Loop is over W levels

! Skip this K point if either upper or lower cell center is outside plot window. 

   if (zt(k) < op%ymin .or. zt(k+1) > op%ymax) cycle
   
! Get T-cell vertical coordinates

      vcpn(1) = zt(k)
      vcpn(2) = vcpn(1)
      vcpn(3) = zt(k+1)
      vcpn(4) = vcpn(3)

! special for dudhia expts
!if (vcpn(4) > op%ymax) cycle
! end special

! Fill field values of 4 T points around current M point

      call oplot_lib(k,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldvals(1),notavail)
      if (notavail > 0) cycle 
      call oplot_lib(k,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldvals(2),notavail)
      if (notavail > 0) cycle 
      call oplot_lib(k+1,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldvals(3),notavail)
      if (notavail > 0) cycle 
      call oplot_lib(k+1,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldvals(4),notavail)
      if (notavail > 0) cycle 

! Contour plot cell around current M point

      call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)
      
   enddo

enddo
   
return
end subroutine contslab_vert_v

!===============================================================================

subroutine contslab_vert_tw(iplt)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, lpw, zm, zt, &
                      xem, yem, zem, xev, yev, zev, xew, yew, zew
use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
use misc_coms,  only: io6, mdomain
use consts_coms, only: erad, pio180

implicit none

integer, intent(in) :: iplt

integer :: k,k1,k2,iw,iw1,iw2,iw3,iv2,im1,im2,jv1,jv2,jv3,jw
integer :: npoly,iok,notavail,itri
real :: radcone
real :: wtbot = 1., wttop = 0.
real :: hcpn(4), vcpn(4), fldvals(4)
real :: valw(2), valw1(2), valw2(2), valw3(2), valm1(2), valm2(2), valv2(2)
real :: wta1, wta2, wta3, wtb1, wtb2, wtb3

! First plot underground T cells with underground color

call plot_underground_w(iplt,(/0/))

! Loop is over W for contouring W points
! Limit to primary W point or else we will go out of bounds for a parallel run

do jw = 1, jtab_w(jtw_prog)%jend(1)
   iw = jtab_w(jtw_prog)%iw(jw)

   npoly = itab_w(iw)%npoly

! Loop over U/V neighbors of W, and for each, get indices for
! M endpoints and opposite W neighbor

   do jv2 = 1,npoly

      jv1 = jv2 - 1
      if (jv2 == 1) jv1 = npoly

      jv3 = jv2 + 1
      if (jv2 == npoly) jv3 = 1

      iv2 = itab_w(iw)%iv(jv2)

      im1 = itab_w(iw)%im(jv1) ! check these
      im2 = itab_w(iw)%im(jv2) ! check these

      iw1 = itab_w(iw)%iw(jv1)
      iw2 = itab_w(iw)%iw(jv2)
      iw3 = itab_w(iw)%iw(jv3)

! Loop over both triangles in current JV sector

      do itri = 1,2

! Check for intersection of cone surface and current triangle

         if (itri == 1) then
            call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
               xem(im1),yem(im1),zem(im1),xev(iv2),yev(iv2),zev(iv2), &
               wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
         else
            call coneplot_tri(iplt,iw,xew(iw),yew(iw),zew(iw), &
               xev(iv2),yev(iv2),zev(iv2),xem(im2),yem(im2),zem(im2), &
               wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,hcpn)
         endif
                             
! Jump to end of loop if cone surface does not intersect triangle

         if (iok /= 1) cycle

! Skip triangle if it crosses +/- 180 degree point along cone circle
! (In future, maybe re-program to truncate cells)

         if (mdomain < 2 .and. op%projectn(iplt) == 'C') then
            radcone = erad * sin(op%coneang * pio180)

            if (abs(hcpn(2) - hcpn(1)) > 3. * radcone) cycle
         endif

! Skip this triangle if either cell side is outside plot window. 

         if (hcpn(1) < op%xmin .or. hcpn(1) > op%xmax .or.  &
             hcpn(2) < op%xmin .or. hcpn(2) > op%xmax) cycle         

         do k = lpw(iw),mza+1   ! Loop is over T levels

! Use k = mza+1 level only for contouring field values defined at T points
! (in which case the top half of the mza cell is contoured)

            if (k == mza+1 .and. op%stagpt == 'W') cycle

! Define vertical plot coordinates depending on whether field values
! are defined at T or W points

            if (op%stagpt == 'T') then
               vcpn(1:2) = zt(k-1)
               vcpn(3:4) = zt(k)
               if (k == lpw(iw)) vcpn(1:2) = zm(k-1)
               if (k == mza+1)     vcpn(3:4) = zm(k-1)
            else
               vcpn(1:2) = zm(k-1)
               vcpn(3:4) = zm(k)
            endif

! Fill field values at two adjacent vertical levels at each of three IW pts

            if (k < mza+1) then
               call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                              valw(2),notavail)
               if (notavail > 0) cycle 

               call oplot_lib(k,iw2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw2(2),notavail)
               if (notavail > 0) cycle
            endif

            if (k > lpw(iw) .or. op%stagpt == 'W') then
               call oplot_lib(k-1,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                              valw(1),notavail)
               if (notavail > 0) cycle 

               call oplot_lib(k-1,iw2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                             valw2(1),notavail)
               if (notavail > 0) cycle
            else
               valw (1) = valw (2)
               valw2(1) = valw2(2) 
            endif

! W values interpolated to V2

            valv2(1) = .5 * (valw(1) + valw2(1))
            valv2(2) = .5 * (valw(2) + valw2(2))

            if (itri == 1) then

               if (k < mza+1) then
                  call oplot_lib(k,iw1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                 valw1(2),notavail)
                  if (notavail > 0) cycle 
               endif

               if (k > lpw(iw) .or. op%stagpt == 'W') then
                  call oplot_lib(k-1,iw1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                 valw1(1),notavail)
                  if (notavail > 0) cycle 
               else
                  valw1(1) = valw1(2) 
               endif

               if (k == mza+1) then
                  valw (2) = valw (1)
                  valw1(2) = valw1(1)
                  valw2(2) = valw2(1)
               endif

! W values interpolated to M1: For now, take simple average

               valm1(1) = .3333333 * (valw(1) + valw1(1) + valw2(1))
               valm1(2) = .3333333 * (valw(2) + valw1(2) + valw2(2))

! W values interpolated to A and B points for plotting

               fldvals(1) = wta1 * valw(1) + wta2 * valm1(1) + wta3 * valv2(1)
               fldvals(2) = wtb1 * valw(1) + wtb2 * valm1(1) + wtb3 * valv2(1)
               fldvals(3) = wtb1 * valw(2) + wtb2 * valm1(2) + wtb3 * valv2(2)
               fldvals(4) = wta1 * valw(2) + wta2 * valm1(2) + wta3 * valv2(2)

           else

               if (k < mza+1) then
                  call oplot_lib(k,iw3,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                 valw3(2),notavail)
                  if (notavail > 0) cycle 
               endif

               if (k > lpw(iw) .or. op%stagpt == 'W') then
                  call oplot_lib(k-1,iw3,'VALUE',op%fldname(iplt),wtbot,wttop, &
                                 valw3(1),notavail)
                  if (notavail > 0) cycle
               else
                  valw3(1) = valw3(2) 
               endif

               if (k == mza+1) then
                  valw (2) = valw (1)
                  valw2(2) = valw2(1)
                  valw3(2) = valw3(1)
               endif

! W values interpolated to M2: For now, take simple average

               valm2(1) = .3333333 * (valw(1) + valw2(1) + valw3(1))
               valm2(2) = .3333333 * (valw(2) + valw2(2) + valw3(2))

! W values interpolated to A and B points for plotting

               fldvals(1) = wta1 * valw(1) + wta2 * valv2(1) + wta3 * valm2(1)
               fldvals(2) = wtb1 * valw(1) + wtb2 * valv2(1) + wtb3 * valm2(1)
               fldvals(3) = wtb1 * valw(2) + wtb2 * valv2(2) + wtb3 * valm2(2)
               fldvals(4) = wta1 * valw(2) + wta2 * valv2(2) + wta3 * valm2(2)

            endif ! itri

! Contour plot cell around current M point

            call contpolyg(op%icolortab(iplt),op%ifill,4,hcpn,vcpn,fldvals)

         enddo ! K

      enddo ! ITRI

   enddo ! JV

enddo ! IW

return
end subroutine contslab_vert_tw
