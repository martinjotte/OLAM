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
subroutine draw_cellbnd(ntpn,xtpn,ytpn)

implicit none

integer, intent(in) :: ntpn
real, intent(in) :: xtpn(ntpn),ytpn(ntpn)

integer :: itpn

call o_frstpt(xtpn(ntpn),ytpn(ntpn))
do itpn = 1,ntpn
   call o_vector(xtpn(itpn),ytpn(itpn))
enddo

return
end subroutine draw_cellbnd

!===============================================================================

subroutine ll_unwrap(xpt,xtpn) 

implicit none

real, intent(in) :: xpt
real, intent(inout) :: xtpn

if (xpt < -90. .and. xtpn > 0.) then
   xtpn = xtpn - 360.
elseif (xpt > 90. .and. xtpn < 0.) then
   xtpn = xtpn + 360.
endif

return
end subroutine ll_unwrap

!===============================================================================

subroutine oplot_prtvalue(value,xpt,ypt,vsprd,psiz,itab)

use plotcolors, only: clrtab
use misc_coms,  only: io6

implicit none

real, intent(in) :: value,xpt,ypt,vsprd,psiz
integer, intent(in) :: itab

integer :: i,ng,ip,iw

character(len=20) :: number,numbr
integer :: ln

call o_gsplci(10)
call o_gstxci(10)

if     (clrtab(itab)%ifmt(2) == 7) then
   write (number,'(f12.7)') value
elseif (clrtab(itab)%ifmt(2) == 6) then
   write (number,'(f11.6)') value
elseif (clrtab(itab)%ifmt(2) == 5) then
   write (number,'(f10.5)') value
elseif (clrtab(itab)%ifmt(2) == 4) then
   write (number,'(f9.4)') value
elseif (clrtab(itab)%ifmt(2) == 3) then
   write (number,'(f9.3)') value
elseif (clrtab(itab)%ifmt(2) == 2) then
   write (number,'(f9.2)') value
elseif (clrtab(itab)%ifmt(2) == 1) then
   write (number,'(f9.1)') value
elseif (clrtab(itab)%ifmt(2) == 0 .or. clrtab(itab)%ifmt(2) == 10) then
   write (number,'(f7.0)') value
elseif (clrtab(itab)%ifmt(2) == -1) then
   write (number,'(e9.1)') value
elseif (clrtab(itab)%ifmt(2) == -2) then
   write (number,'(e10.2)') value
elseif (clrtab(itab)%ifmt(2) == -3) then
   write (number,'(e11.3)') value
elseif (clrtab(itab)%ifmt(2) == -4) then
   write (number,'(e12.4)') value
elseif (clrtab(itab)%ifmt(2) == -5) then
   write (number,'(e13.5)') value
elseif (clrtab(itab)%ifmt(2) == -6) then
   write (number,'(e14.6)') value
else
   write(number,'(e15.7)') value
endif

numbr = adjustl(number)
if (numbr(len_trim(numbr):len_trim(numbr)) == '.') then
   ln = len_trim(numbr) - 1
else
   ln = len_trim(numbr)
endif

call o_plchlq(xpt,ypt+.4*vsprd,numbr(1:ln),psiz,0.,0.)

return
end subroutine oplot_prtvalue

!===============================================================================

subroutine oplot_locindex(i,x,y,psiz)

implicit none

integer, intent(in) :: i
real, intent(in) :: x,y,psiz

integer :: ln
character(len=8) :: number3

!call o_gsplci(10)
!call o_gstxci(10)

write(number3,'(i6)') i

call o_plchlq (x,y,trim(adjustl(number3)),psiz,0.,0.)

return
end subroutine oplot_locindex

!===============================================================================

subroutine contpolyg(itab,ifill,npts,xin,yin,zin)

use oplot_coms, only: op

implicit none

integer, intent(in) :: itab,ifill,npts
real, intent(in) :: xin(npts),yin(npts),zin(npts)

integer :: ntri,ntri2,ipt
real :: x1,x2,x3,y1,y2,y3,z1,z2,z3,c1
real :: fldval_min,fldval_max

if (npts == 3) then

   x1  = xin(1)
   x2  = xin(2)
   x3  = xin(3)
   y1  = yin(1)
   y2  = yin(2)
   y3  = yin(3)
   z1  = zin(1)
   z2  = zin(2)
   z3  = zin(3)
   call cont3(itab,ifill,z1,z2,z3,x1,x2,x3,y1,y2,y3)

elseif (npts > 3) then

! Define "3" point as average of all others

   x3 = 0.
   y3 = 0.
   z3 = 0.
   c1 = 1. / real(npts)

   do ipt = 1,npts
      x3 = x3 + xin(ipt)
      y3 = y3 + yin(ipt)
      z3 = z3 + zin(ipt)
   enddo

   x3 = x3 * c1
   y3 = y3 * c1
   z3 = z3 * c1

   do ntri = 1,npts
      ntri2 = mod(ntri,npts) + 1

      x1  = xin(ntri)
      y1  = yin(ntri)
      z1  = zin(ntri)

      x2  = xin(ntri2)
      y2  = yin(ntri2)
      z2  = zin(ntri2)

      call cont3(itab,ifill,z1,z2,z3,x1,x2,x3,y1,y2,y3)
   enddo

endif

fldval_min = minval(zin(1:npts))
fldval_max = maxval(zin(1:npts))

if (op%fldval_min > fldval_min) op%fldval_min = fldval_min
if (op%fldval_max < fldval_max) op%fldval_max = fldval_max

return
end subroutine contpolyg

!===============================================================================

subroutine cont3(itab,ifill,z1,z2,z3,x1,x2,x3,y1,y2,y3)

implicit none

integer, intent(in) :: itab,ifill
real, intent(in) :: z1,z2,z3,x1,x2,x3,y1,y2,y3

if     (z1 >= z2 .and. z2 >= z3) then
   call cont3lf(itab,ifill,z1,z2,z3,x1,x2,x3,y1,y2,y3)
elseif (z1 >= z3 .and. z3 >= z2) then
   call cont3lf(itab,ifill,z1,z3,z2,x1,x3,x2,y1,y3,y2)
elseif (z2 >= z1 .and. z1 >= z3) then
   call cont3lf(itab,ifill,z2,z1,z3,x2,x1,x3,y2,y1,y3)
elseif (z2 >= z3 .and. z3 >= z1) then
   call cont3lf(itab,ifill,z2,z3,z1,x2,x3,x1,y2,y3,y1)
elseif (z3 >= z1 .and. z1 >= z2) then
   call cont3lf(itab,ifill,z3,z1,z2,x3,x1,x2,y3,y1,y2)
elseif (z3 >= z2 .and. z2 >= z1) then
   call cont3lf(itab,ifill,z3,z2,z1,x3,x2,x1,y3,y2,y1)
endif

return
end subroutine cont3

!===============================================================================
 
subroutine cont3lf(itab,ifill,z1,z2,z3,x1,x2,x3,y1,y2,y3)

implicit none

integer, intent(in) :: itab,ifill
real, intent(in) :: z1,z2,z3,x1,x2,x3,y1,y2,y3

if (ifill == 0) then
   call cont3l(itab,z1,z2,z3,x1,x2,x3,y1,y2,y3)
else
   call cont3f(itab,z1,z2,z3,x1,x2,x3,y1,y2,y3)
endif
return
end subroutine cont3lf 

!===============================================================================

subroutine cont3l(itab,z1,z2,z3,x1,x2,x3,y1,y2,y3)

use plotcolors, only: clrtab
use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: itab
real, intent(in) :: z1,z2,z3,x1,x2,x3,y1,y2,y3

integer :: ival,iskip

real :: xp1,xp2,yp1,yp2,contlev
real :: xq1,xq2,yq1,yq2

real :: epsw,frac_remain

! Extract lowest contour value from (color) table that exceeds 
! the value of z3, which is the lowest data value

ival = 1
do while (z3 >= clrtab(itab)%vals(ival) .and. ival < clrtab(itab)%nvals)
   ival = ival + 1
enddo
contlev = clrtab(itab)%vals(ival)

! If all 3 points are in the same contour interval, return

if (contlev >= z1 .or. contlev <= z3) then
   return
endif

! Draw all contour lines of value less than z1

do while (contlev < z1 .and. ival <= clrtab(itab)%nvals)

   xp1 = x3 + (x1-x3) * (contlev-z3) / (z1-z3)
   yp1 = y3 + (y1-y3) * (contlev-z3) / (z1-z3)

   if (contlev < z2) then
      xp2 = x3 + (x2-x3) * (contlev-z3) / (z2-z3)
      yp2 = y3 + (y2-y3) * (contlev-z3) / (z2-z3)
   else
      xp2 = x2 + (x1-x2) * (contlev-z2) / (z1-z2)
      yp2 = y2 + (y1-y2) * (contlev-z2) / (z1-z2)
   endif

! Truncate segment if it crosses plot window boundary or skip if both endpoints
! are outside plot window

   call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

!hs SPECIAL FOR HS EXPERIMENT:  

!hs  Subroutine trunc_segment does not work correctly
!hs for plots that have pressure as vertical coordinate because pressure values
!hs increase downward.  Thus, comment out trunc_segment call above.

!hs Also, uncomment following 4 lines to bypass call to trunc_segment

!hs xq1 = xp1
!hs xq2 = xp2
!hs yq1 = yp1
!hs !hs yq2 = yp2

!hs END SPECIAL

   if (iskip == 1) go to 9

! Plot line segment

   call o_frstpt (xq1,yq1)
   call o_vector (xq2,yq2)

9 continue
   
   ival = ival + 1
   contlev = clrtab(itab)%vals(ival)

enddo

return
end subroutine cont3l

!===============================================================================
 
subroutine cont3f(itab,z1,z2,z3,x1,x2,x3,y1,y2,y3)

use plotcolors, only: clrtab
use misc_coms,  only: io6

implicit none

integer, intent(in) :: itab
real, intent(in) :: z1,z2,z3,x1,x2,x3,y1,y2,y3

integer :: ival,icolor,iflag
real :: xo1,xo2,yo1,yo2 ,contlev    ! "old" values of xcpn1,xcpn2,ycpn1,ycpn2
real, dimension(5) :: xcpn,ycpn

! Extract lowest contour value from color table that exceeds 
! the value of z3, which is the lowest data value

ival = 1
do while (z3 >= clrtab(itab)%vals(ival) .and. ival < clrtab(itab)%nvals)
   ival = ival + 1
enddo
contlev = clrtab(itab)%vals(ival)

! If all 3 points are in the same contour interval, tileplot the
! triangle with single color

if (contlev >= z1 .or. contlev <= z3) then

   xcpn(1) = x1
   ycpn(1) = y1
   xcpn(2) = x2
   ycpn(2) = y2
   xcpn(3) = x3
   ycpn(3) = y3
   icolor = clrtab(itab)%ipal(ival)

   call fillpolyg(3,xcpn,ycpn,icolor)
   return
endif

! Draw all contour lines of value less than z1

iflag = 0

do while (contlev < z1)

   if (iflag > 0) then
      xo1 = xcpn(1)
      yo1 = ycpn(1)
      xo2 = xcpn(2)
      yo2 = ycpn(2)
   endif
    
   xcpn(1) = x3 + (x1-x3) * (contlev-z3) / (z1-z3)
   ycpn(1) = y3 + (y1-y3) * (contlev-z3) / (z1-z3)

   if (contlev < z2) then

      xcpn(2) = x3 + (x2-x3) * (contlev-z3) / (z2-z3)
      ycpn(2) = y3 + (y2-y3) * (contlev-z3) / (z2-z3)

      if (iflag == 0) then  ! lowest contour color: z3 is a node
         xcpn(3) = x3   
         ycpn(3) = y3
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(3,xcpn,ycpn,icolor)
      else
         xcpn(3) = xo2
         ycpn(3) = yo2
         xcpn(4) = xo1
         ycpn(4) = yo1
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(4,xcpn,ycpn,icolor)
      endif

      if (ival == clrtab(itab)%nvals) then  ! highest table color:
         xcpn(3) = x2                       ! z1 is a node
         ycpn(3) = y2
         xcpn(4) = x1         
         ycpn(4) = y1
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(4,xcpn,ycpn,icolor)
         return
      elseif (clrtab(itab)%vals(ival+1) >= z1) then  ! highest table color:
         xcpn(3) = x2                                ! z1 is a node
         ycpn(3) = y2
         xcpn(4) = x1         
         ycpn(4) = y1
         icolor = clrtab(itab)%ipal(ival+1)
         call fillpolyg(4,xcpn,ycpn,icolor)
         return
      endif
      
      iflag = 1
         
   else   
      
      xcpn(2) = x2 + (x1-x2) * (contlev-z2) / (z1-z2)
      ycpn(2) = y2 + (y1-y2) * (contlev-z2) / (z1-z2)
      
      if (iflag == 0) then  ! lowest contour color: z3 is a node
         xcpn(3) = x2
         ycpn(3) = y2
         xcpn(4) = x3        
         ycpn(4) = y3         
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(4,xcpn,ycpn,icolor)
      elseif (iflag == 1) then ! switching from contlev < z2 to contlev > z2
         xcpn(3) = x2
         ycpn(3) = y2
         xcpn(4) = xo2
         ycpn(4) = yo2
         xcpn(5) = xo1         
         ycpn(5) = yo1
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(5,xcpn,ycpn,icolor)
      else                     ! keeping to contlev > z2
         xcpn(3) = xo2
         ycpn(3) = yo2
         xcpn(4) = xo1
         ycpn(4) = yo1
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(4,xcpn,ycpn,icolor)
      endif
      
      if (ival == clrtab(itab)%nvals) then  ! highest table color:
         xcpn(3) = x1                       ! z1 is a node
         ycpn(3) = y1
         icolor = clrtab(itab)%ipal(ival)
         call fillpolyg(3,xcpn,ycpn,icolor)
         return
      elseif (clrtab(itab)%vals(ival+1) >= z1) then  ! highest table color:
         xcpn(3) = x1                                ! z1 is a node
         ycpn(3) = y1
         icolor = clrtab(itab)%ipal(ival+1)
         call fillpolyg(3,xcpn,ycpn,icolor)
         return
      endif

      iflag = 2
         
   endif

   ival = ival + 1
   contlev = clrtab(itab)%vals(ival)
   
enddo

return
end subroutine cont3f

!===============================================================================

subroutine fillpolyg (np,xp,yp,icolor)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: np,icolor
real, intent(in) :: xp(np),yp(np)

integer :: ip,ipnext
integer :: iq,iqnext,nq
integer :: ir,irnext,nr

real :: xr(10),yr(10)

! Truncate polygon to plot window

call trunc_polyg(np,xp,yp,nr,xr,yr)

!hs SPECIAL FOR HS EXPERIMENT:  

!hs Subroutine trunc_polyg does not work correctly
!hs for plots that have pressure as vertical coordinate because pressure values
!hs increase downward.  Thus, comment out trunc_polyg call above.

!hs Also, uncomment following 3 lines to bypass call to trunc_polyg

!hs nr = np
!hs xr(:) = xp(:)
!hs yr(:) = yp(:)

! END SPECIAL

! Fill the polygon with single color

if (nr > 2) call o_sfsgfa (xr,yr,nr,icolor)
         
return
end subroutine fillpolyg

!===============================================================================

subroutine trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(out) :: iskip

real, intent(in)  :: xp1,xp2,yp1,yp2
real, intent(out) :: xq1,xq2,yq1,yq2

real :: frac_remain

! Initialize skip flag to zero

iskip = 0

xq1 = xp1
xq2 = xp2
yq1 = yp1
yq2 = yp2

! If either segment endpoint is more than 1.e11 m from plot center, it represents
! "missing" point on other side of earth in orthographic projection.
! Set iskip flag and return.

if (abs(xp1) > 1.e11 .or. abs(xp2) > 1.e11) then
   iskip = 1
   return
endif           
           
! If both P segment endpoints are outside same plot window boundary, 
! set iskip flag and return

if ((xp1 < op%xmin .and. xp2 < op%xmin) .or.  &
    (xp1 > op%xmax .and. xp2 > op%xmax) .or.  &
    (yp1 < op%ymin .and. yp2 < op%ymin) .or.  &
    (yp1 > op%ymax .and. yp2 > op%ymax)) then

   iskip = 1
   return
endif           
           
! Truncate line segment if it crosses x-boundary, removing any portion that is
! outside plot window

if (xp1 < op%xmin) then
   frac_remain = (op%xmin - xp2) / (xp1 - xp2)
   xq1 = op%xmin
   yq1 = yp2 + frac_remain * (yp1 - yp2)
elseif (xp1 > op%xmax) then
   frac_remain = (op%xmax - xp2) / (xp1 - xp2)
   xq1 = op%xmax
   yq1 = yp2 + frac_remain * (yp1 - yp2)
elseif (xp2 < op%xmin) then
   frac_remain = (op%xmin - xp1) / (xp2 - xp1)
   xq2 = op%xmin
   yq2 = yp1 + frac_remain * (yp2 - yp1)
elseif (xp2 > op%xmax) then
   frac_remain = (op%xmax - xp1) / (xp2 - xp1)
   xq2 = op%xmax
   yq2 = yp1 + frac_remain * (yp2 - yp1)
endif

! If both Q segment endpoints are outside same plot window boundary, 
! set iskip flag and return

if ((xq1 < op%xmin .and. xq2 < op%xmin) .or.  &
    (xq1 > op%xmax .and. xq2 > op%xmax) .or.  &
    (yq1 < op%ymin .and. yq2 < op%ymin) .or.  &
    (yq1 > op%ymax .and. yq2 > op%ymax)) then

   iskip = 1
   return
endif           
           
! Truncate line segment if it crosses y-boundary, removing any portion that is
! outside plot window

if (yq1 < op%ymin) then
   frac_remain = (op%ymin - yq2) / (yq1 - yq2)
   xq1 = xq2 + frac_remain * (xq1 - xq2)
   yq1 = op%ymin
elseif (yq1 > op%ymax) then
   frac_remain = (op%ymax - yq2) / (yq1 - yq2)
   xq1 = xq2 + frac_remain * (xq1 - xq2)
   yq1 = op%ymax
elseif (yq2 < op%ymin) then
   frac_remain = (op%ymin - yq1) / (yq2 - yq1)
   xq2 = xq1 + frac_remain * (xq2 - xq1)
   yq2 = op%ymin
elseif (yq2 > op%ymax) then
   frac_remain = (op%ymax - yq1) / (yq2 - yq1)
   xq2 = xq1 + frac_remain * (xq2 - xq1)
   yq2 = op%ymax
endif

return

end subroutine trunc_segment

!===============================================================================

subroutine trunc_polyg (np,xp,yp,nr,xr,yr)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: np
real,    intent(in) :: xp(np),yp(np)

integer, intent(out) :: nr
real,    intent(out) :: xr(10),yr(10)

integer :: ip,ipnext
integer :: iq,iqnext,nq

real :: xq(10),yq(10)

real :: frac_remain

! Assign xq,yq based on how xp(ip) are located inside or outside X boundaries

nq = 0

! Loop over all input ip points

do ip = 1,np

! Get index for next ip point.  We will consider segment between them.

   ipnext = ip + 1
   if (ip == np) ipnext = 1

! If current ip point is in bounds, copy to (xq,yq)

   if (xp(ip) >= op%xmin .and. xp(ip) <= op%xmax) then
      nq = nq + 1
      xq(nq) = xp(ip)
      yq(nq) = yp(ip)

! If next ip point is also in bounds, cycle

      if (xp(ipnext) >= op%xmin .and. xp(ipnext) <= op%xmax) then
         cycle

! If next ip point is out of bounds, check if it is out on left side

      elseif (xp(ipnext) < op%xmin) then

! Interpolate segment between current and next ip pts

         frac_remain = (op%xmin - xp(ip)) / (xp(ipnext) - xp(ip))
         nq = nq + 1
         xq(nq) = op%xmin
         yq(nq) = yp(ip) + frac_remain * (yp(ipnext) - yp(ip))
         cycle
         
      else

! Otherwise, ip point is out of bounds on right
! Interpolate segment between current and next ip pts

         frac_remain = (op%xmax - xp(ip)) / (xp(ipnext) - xp(ip))
         nq = nq + 1
         xq(nq) = op%xmax
         yq(nq) = yp(ip) + frac_remain * (yp(ipnext) - yp(ip))
         cycle

      endif

   endif

! If we got here, ip point is out of bounds.  Check if it is out on left side.

   if (xp(ip) <= op%xmin) then

! Ip point is outside on left side.  If next ip point is also, skip this point.

      if (xp(ipnext) < op%xmin) then

         cycle

      else

! If next ip point is inside, interpolate between current and next ip points

         frac_remain = (op%xmin - xp(ipnext)) / (xp(ip) - xp(ipnext))
         nq = nq + 1
         xq(nq) = op%xmin
         yq(nq) = yp(ipnext) + frac_remain * (yp(ip) - yp(ipnext))

      endif

   else

! If we got here, ip point is out of bounds on right side. If next point is also,
! skip this ip point.

      if (xp(ipnext) > op%xmax) then

         cycle

      else

! If next ip point is inside, interpolate between current and next ip points

         frac_remain = (op%xmax - xp(ipnext)) / (xp(ip) - xp(ipnext))
         nq = nq + 1
         xq(nq) = op%xmax
         yq(nq) = yp(ipnext) + frac_remain * (yp(ip) - yp(ipnext))

      endif

   endif     

enddo

! Assign xr,yr based on how xq(iq) are located inside or outside Y boundaries

nr = 0

! Loop over all iq points

do iq = 1,nq

! Get index for next iq point.  We will consider segment between them.

   iqnext = iq + 1
   if (iq == nq) iqnext = 1

! If current iq point is in bounds, copy to (xr,yr)

   if (yq(iq) >= op%ymin .and. yq(iq) <= op%ymax) then
      nr = nr + 1
      xr(nr) = xq(iq)
      yr(nr) = yq(iq)

! If next iq point is also in bounds, cycle

      if (yq(iqnext) >= op%ymin .and. yq(iqnext) <= op%ymax) then
         cycle

! If next iq point is out of bounds, check if it is out on bottom side

      elseif (yq(iqnext) < op%ymin) then

! Interpolate segment between current and next iq pts

         frac_remain = (op%ymin - yq(iq)) / (yq(iqnext) - yq(iq))
         nr = nr + 1
         xr(nr) = xq(iq) + frac_remain * (xq(iqnext) - xq(iq))
         yr(nr) = op%ymin
         cycle
         
      else

! Otherwise, iq point is out of bounds on top side
! Interpolate segment between current and next pts

         frac_remain = (op%ymax - yq(iq)) / (yq(iqnext) - yq(iq))
         nr = nr + 1
         xr(nr) = xq(iq) + frac_remain * (xq(iqnext) - xq(iq))
         yr(nr) = op%ymax
         cycle

      endif

   endif

! If we got here, iq point is out of bounds.  Check if it is out on bottom side.

   if (yq(iq) <= op%ymin) then

! Point is outside on bottom side.  If next point is also, skip this point.

      if (yq(iqnext) < op%ymin) then
         cycle
      else

! If next point is inside, interpolate between current and next

         frac_remain = (op%ymin - yq(iqnext)) / (yq(iq) - yq(iqnext))
         nr = nr + 1
         xr(nr) = xq(iqnext) + frac_remain * (xq(iq) - xq(iqnext))
         yr(nr) = op%ymin

      endif

   else

! If we got here, iq point is out of bounds on top side.
! If next point is also, skip this point.

      if (yq(iqnext) > op%ymax) then
         cycle
      else

! If next point is in, interpolate between current and next

         frac_remain = (op%ymax - yq(iqnext)) / (yq(iq) - yq(iqnext))
         nr = nr + 1
         xr(nr) = xq(iqnext) + frac_remain * (xq(iq) - xq(iqnext))
         yr(nr) = op%ymax

      endif

   endif     

enddo

return
end subroutine trunc_polyg

!===============================================================================

subroutine plotback()

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

real,    dimension(5)  :: rif, rjf
real,    dimension(8)  :: dst, ind
integer, dimension(13) :: iasf = 1
    
! This subroutine plots a white background for the entire plot.

call o_gsclip (0)

! Set all the GKS aspect source flags to "individual".

call o_gsasf (iasf)

! Force solid fill.

call o_gsfais (1)

call o_set (0.,1.,0.,1.,0.,1.,0.,1.,1)

call o_sfseti ('TYPE OF FILL',0)
     
rif(1) = 0.
rif(2) = 1.
rif(3) = 1.
rif(4) = 0.
rif(5) = 0.  ! not needed?

rjf(1) = 0.
rjf(2) = 0.
rjf(3) = 1.
rjf(4) = 1.
rjf(5) = 0.   ! not needed?

! The last parameter in the following call sets the color.

call o_sfsgfa (rif,rjf,4,7)

op%iplotback = 1

return
end subroutine plotback

!===============================================================================

subroutine plot_labelbar (iplt,fldname,units,projectn,slabloc  &
   ,slabmin,slabmax,time_istp8,ihour,idate,imonth,iyear)

use oplot_coms, only: op
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt,ihour,idate,imonth,iyear
character(len=*), intent(in) :: fldname,units,projectn
real, intent(in) :: slabloc,slabmin,slabmax
real(kind=8), intent(in) :: time_istp8

integer :: lc,lcdname,lcdunits
real :: bsize,xlabx
character(len=60) :: title,title1
character(len=5), dimension(12) :: plotmonth

data plotmonth(1:12) / ' JAN ',' FEB ',' MAR ',' APR ',' MAY ',' JUN '  &
                      ,' JUL ',' AUG ',' SEP ',' OCT ',' NOV ',' DEC '  /

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Specify font # and scale font size to designated plotter coordinates

call o_pcsetr('CL',1.)  ! set character line width to 1.
call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
bsize = .010 * (op%hp2 - op%hp1)

! Set color

call o_gsplci(10)
call o_gstxci(10)
call o_sflush()

! Plot field name and units

call o_plchhq(op%fx1,op%fnamey,trim(fldname)//trim(units),1.0*bsize,0.,-1.)

! Plot UTC time and date

write (title, '(i6.6,a4)') ihour,' UTC'
write (title1,'(3x,i2.2,a5,i4.4)') idate,plotmonth(imonth),iyear
call o_plchhq(op%fx2,op%fnamey,trim(title)//trim(title1),bsize,0.,1.)

! RETURN if not plotting information block

if (op%labelbar(iplt) /= 'i') return

! Plot elapsed time (sec) in simulation

write (title,'(F11.1)') time_istp8
call o_plchhq(op%timex,op%timsy,trim(adjustl(title))//' sec',bsize,0.,-1.)

! Plot elapsed time (day) in simulation

write (title,'(F10.2)') time_istp8/86400.
call o_plchhq(op%timex,op%timdy,trim(adjustl(title))//' days',bsize,0.,-1.)

! Plot slab location

write (title,'(F12.1)') slabloc

if (projectn == 'L' .or. projectn == 'P' .or.  &
    projectn == 'O' .or. projectn == 'Z') then
    
   if (op%pltlev(iplt) == 'p') then 
      call o_plchhq(op%slabx,op%slaby  &
         ,'P   '//trim(adjustl(title))//' mb',bsize,0.,-1.)
   elseif (op%pltlev(iplt) == 's') then 
      call o_plchhq(op%slabx,op%slaby  &
         ,'LOWEST LEVEL PLOT ',bsize,0.,-1.)
   else
      call o_plchhq(op%slabx,op%slaby  &
         ,'Z   '//trim(adjustl(title))//' m',bsize,0.,-1.)
   endif

elseif (projectn == 'X') then
   call o_plchhq(op%slabx,op%slaby  &
      ,'X   '//trim(adjustl(title))//' m',bsize,0.,-1.)
elseif (projectn == 'Y') then
   call o_plchhq(op%slabx,op%slaby  &
      ,'Y   '//trim(adjustl(title))//' m',bsize,0.,-1.)
elseif (projectn == 'V') then
   call o_plchhq(op%slabx,op%slaby  &
      ,'V   '//trim(adjustl(title))//' m',bsize,0.,-1.)
elseif (projectn == 'C') then
   write (title,'(F7.2)')  slabloc
   call o_plchhq(op%slabx,op%slaby  &
      ,'C   '//trim(adjustl(title))//' deg',bsize,0.,-1.)
endif

! Plot min value in plotted slab

write(title,'(g13.6)') slabmin
call o_plchhq (op%slabx,op%sminy,'MIN '//trim(adjustl(title)),bsize,0.,-1.)

! Plot max value in plotted slab

write(title,'(g13.6)') slabmax
call o_plchhq (op%slabx,op%smaxy,'MAX '//trim(adjustl(title)),bsize,0.,-1.)

return
end subroutine plot_labelbar

!===============================================================================

subroutine plot_colorbar(iplt,itab)

use oplot_coms, only: op
use plotcolors, only: clrtab
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt,itab

integer :: ibox,jbox,ln
real :: bsize,yinc
real, dimension(4) :: xbox,ybox
character(len=14)  :: number,numbr

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Specify font # and scale font size to designated plotter coordinates

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',1.)  ! set character line width to 1
bsize = .009 * (op%hp2 - op%hp1)

call o_gsplci(10)
call o_gstxci(10)

! Compute height increment of colorbar boxes

yinc = (op%fy2 - op%fy1) / float(clrtab(itab)%nvals)

! Loop over colorbar boxes

do ibox = 1,clrtab(itab)%nvals

! Skip bottom and top boxes if plotting integer values at center of box

   if (ibox == 1                  .and. clrtab(itab)%ifmt(1) >= 10) cycle
   if (ibox == clrtab(itab)%nvals .and. clrtab(itab)%ifmt(1) >= 10) cycle

! Copy colorbar left and right x coords

   xbox(1) = op%cbx1
   xbox(2) = op%cbx2
   xbox(3) = xbox(2)
   xbox(4) = xbox(1)

! Turn bottom and top boxes into triangles

   if (ibox == 1) then
      xbox(1) = op%cbx1 + .45 * (op%cbx2 - op%cbx1)
      xbox(2) = op%cbx1 + .55 * (op%cbx2 - op%cbx1)
   elseif (ibox == clrtab(itab)%nvals) then
      xbox(3) = op%cbx1 + .55 * (op%cbx2 - op%cbx1)
      xbox(4) = op%cbx1 + .45 * (op%cbx2 - op%cbx1)
   endif

! For current colorbar box, compute top and bottom y coords

   ybox(1) = op%fy1 + float(ibox-1) * yinc
   ybox(2) = ybox(1)
   ybox(3) = ybox(2) + yinc
   ybox(4) = ybox(3)         

! Select print format of colorbar boxes

   if     (clrtab(itab)%ifmt(1) == 5) then
      write (number,'(f8.5)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == 4) then
      write (number,'(f7.4)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == 3) then
      write (number,'(f6.3)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == 2) then
      write (number,'(f6.2)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == 1) then
      write (number,'(f6.1)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) ==  0 .or.  &
           clrtab(itab)%ifmt(1) == 10 .or.  &
           clrtab(itab)%ifmt(1) == 20) then
      write (number,'(f7.0)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -1) then
      write (number,'(e8.1)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -2) then
      write (number,'(e9.2)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -3) then
      write (number,'(e10.3)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -4) then
      write (number,'(e11.4)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -5) then
      write (number,'(e12.5)') clrtab(itab)%vals(ibox)
   elseif (clrtab(itab)%ifmt(1) == -6) then
      write (number,'(e13.6)') clrtab(itab)%vals(ibox)
   else
      write (number,'(e14.7)') clrtab(itab)%vals(ibox)
   endif

   numbr = adjustl(number)
   if(numbr(len_trim(numbr):len_trim(numbr)) == '.') then
      ln = len_trim(numbr) - 1
   else
      ln = len_trim(numbr)
   endif

! Print value of current colorbar box

   if (clrtab(itab)%ifmt(1) >= 10) then     ! integer value at box center
      call o_plchhq (op%cblx,.5*(ybox(1)+ybox(3)),numbr(1:ln),bsize,0.,-1.)
   elseif (ibox < clrtab(itab)%nvals) then  ! real value at box top
      call o_plchhq (op%cblx,ybox(3),numbr(1:ln),bsize,0.,-1.)
   endif

! Plot color in current colorbar box

   call fillpolyg (4,xbox,ybox,clrtab(itab)%ipal(ibox))

! Draw outline around current colorbar box

   call o_gsplci(10)
   call o_gstxci(10)

   call o_frstpt(xbox(4),ybox(4))
   do jbox = 1,4
      call o_vector(xbox(jbox),ybox(jbox))
   enddo

enddo

return
end subroutine plot_colorbar

!===============================================================================

subroutine oplot_xy2(panel,colorbar0,aspect,scalelab &
   ,n,xval,yval,xlab,ylab                            &
   ,xmin,xmax,xinc,labincx  ,ymin,ymax,yinc,labincy  )
 
use oplot_coms, only: op
use misc_coms,  only: io6

! This routine is a substitute for NCAR Graphics routine ezxy to allow 
! control over fonts, labels, axis labels, line width, scaling, etc.
! Pass in a value of 1 for n to not plot (only draw frame and ticks)

implicit none

character(len=1), intent(in) :: panel,colorbar0
real, intent(in) :: aspect,scalelab
integer, intent(in) :: n,labincx,labincy
real, intent(in) :: xval(n),yval(n)
character(len=*), intent(in) :: xlab,ylab
real, intent(in) :: xmin,xmax,xinc,ymin,ymax,yinc

integer :: i,logy,itickvalq
real :: xl,yl,dx,dy,tickval,xmargin,ymargin,sizelab,xlabx,ylaby,tickvalq
character(len=20)  :: numbr,numbr2

! Set plot color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_gstxci(10)
call o_sflush()

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

call oplot_panel(panel,colorbar0,aspect,'N')
call o_set(op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Draw frame

call o_frstpt(op%fx1,op%fy1)
call o_vector(op%fx2,op%fy1)
call o_vector(op%fx2,op%fy2)
call o_vector(op%fx1,op%fy2)
call o_vector(op%fx1,op%fy1)

! Specify font # and scale font size to designated plotter coordinates

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',1.)  ! set character line width to 1
sizelab = scalelab * (op%hp2 - op%hp1)
call o_sflush()

! Write x axis label

xlabx = .5 * (op%fx1 + op%fx2)
call o_plchhq(xlabx,op%xlaby,trim(xlab),sizelab, 0.,0.)

! Write y axis label

ylaby = .5 * (op%fy1 + op%fy2)
call o_plchhq(op%ylabx,ylaby,trim(ylab),sizelab,90.,0.)

! Scale local working window (xmin,xmax,0.,1.) 
! to plotter coordinates (op%h1,op%h2,op%vp1,op%vp2)

call o_set(op%h1,op%h2,op%vp1,op%vp2,xmin,xmax,0.,1.,1)

! Plot and label X-axis ticks

tickval = nint(xmin/xinc) * xinc
if (tickval < xmin - .001 * xinc) tickval = tickval + xinc
   
do while (tickval < xmax + .001 * xinc)
   
   if (mod(nint(tickval/xinc),labincx) == 0) then  ! Only for long ticks

      dy = .014

! Encode and plot current X tick label
 
      if (xinc * labincx >= .999) then
         write (numbr,'(i6)') nint(tickval)
      elseif (xinc * labincx >= .0999) then
         write (numbr,'(f5.1)') tickval
      elseif (xinc * labincx >= .00999) then
         write (numbr,'(f5.2)') tickval
      elseif (xinc * labincx >= .000999) then
         write (numbr,'(f6.3)') tickval
      else
         write (numbr,'(f7.4)') tickval
      endif
   
      call o_plchhq(tickval,op%xtlaby,trim(adjustl(numbr)),sizelab,0.,0.)

   else                                          ! Only for short ticks
      dy = .007
   endif
   
! Plot current X tick
   
   call o_frstpt(tickval,op%fy1)
   call o_vector(tickval,op%fy1 + dy)
   call o_frstpt(tickval,op%fy2)
   call o_vector(tickval,op%fy2 - dy)
   
   tickval = tickval + xinc
enddo

! Scale local working window (0.,1.,ymin,ymax) 
! to plotter coordinates (op%hp1,op%hp2,op%v1,op%v2)

call o_set(op%hp1,op%hp2,op%v1,op%v2,0.,1.,ymin,ymax,1)

! Plot and label Y-axis ticks

tickval = nint(ymin/yinc) * yinc
if ((tickval - ymax) / (ymin - ymax) > 1.001) tickval = tickval + yinc
   
do while ((tickval - ymin) / (ymax - ymin) < 1.001)
   
   if (mod(nint(tickval/yinc),labincy) == 0) then  ! Only for long ticks

      dx = .014
 
! Encode current Y tick label

      if (abs(yinc * labincy) >= .999) then
         write (numbr,'(i6)') nint(tickval)
      elseif (yinc * labincy >= .0999) then
         write (numbr,'(f5.1)') tickval
      elseif (yinc * labincy >= .00999) then
         write (numbr,'(f5.2)') tickval
      elseif (yinc * labincy >= .000999) then
         write (numbr,'(f6.3)') tickval
      else

         logy = int(log10(yinc * labincy)) - 2
         tickvalq = tickval * 10. ** (-logy)         
         itickvalq = nint(tickvalq)

! If significand is at or above 10, reduce

         do while (abs(itickvalq) >= 10)
            logy = logy + 1
            tickvalq = tickvalq * .1
            itickvalq = nint(tickvalq)
         enddo

! Use only one of the following 2 lines

!        write (numbr,'(f4.1)') tickvalq   ! If real significand is required
         write (numbr,'(i2)') itickvalq    ! If integer significand is ok

         write (numbr2,'(i3)') logy
         numbr = trim(adjustl(numbr))
         numbr2 = trim(adjustl(numbr2))

! Determine whether significand, power of 10, or both are to be plotted

         if (itickvalq == 1) then 
            numbr = '10:S3:'//trim(numbr2)//'        '
         elseif (itickvalq == -1) then 
            numbr = '-10:S3:'//trim(numbr2)//'        '
         elseif (itickvalq /= 0) then
            numbr = trim(adjustl(numbr))//'x'//'10:S3:'//trim(adjustl(numbr2))//'        '
         endif

      endif

! Plot Y tick label

!      call o_plchhq(op%ytlabx,tickval,numbr(1:len_trim(numbr)),sizelab,0.,1.)
      call o_plchhq(op%ytlabx,tickval,trim(adjustl(numbr)),sizelab,0.,1.)

   else                                          ! Only for short ticks
      dx = .007
   endif

! Plot current Y tick

   call o_frstpt(op%fx1     ,tickval)
   call o_vector(op%fx1 + dx,tickval)
   call o_frstpt(op%fx2     ,tickval)
   call o_vector(op%fx2 - dx,tickval)
   
   tickval = tickval + yinc
enddo

! Scale local working window (xmin,xmax,ymin,ymax)
!  to plotter coordinates (op%h1,op%h2,op%v1,op%v2)

call o_set(op%h1,op%h2,op%v1,op%v2,xmin,xmax,ymin,ymax,1)

! Plot values

call o_frstpt(xval(1),yval(1))
do i = 2,n
   call o_vector(xval(i),yval(i))
enddo

return
end subroutine oplot_xy2

!===============================================================================

subroutine oplot_xy2log10(panel,colorbar0,aspect,scalelab  &
   ,n,xval,yval,xlab,ylab                         &
   ,xmin,xmax,xinc,labincx  ,logymin,logymax      )
 
use oplot_coms, only: op
use misc_coms,  only: io6

! This routine is a substitute for NCAR Graphics routine ezxy to allow 
! control over fonts, labels, axis labels, line width, scaling, etc.
! Pass in a value of 1 for n to not plot (only draw frame and ticks)

implicit none

character(len=1), intent(in) :: panel,colorbar0
real, intent(in) :: aspect,scalelab
integer, intent(in) :: n,labincx
real, intent(in) :: xval(n),yval(n)
character(len=*), intent(in) :: xlab,ylab
real, intent(in) :: xmin,xmax,xinc
integer, intent(in) :: logymin,logymax

integer :: i,itick,jtick
real :: xl,yl,dx,dy,tickval,xmargin,ymargin,sizelab,xlabx,ylaby,yminlog,ymaxlog
character(len=20)  :: numbr

! Set plot color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_gstxci(10)
call o_sflush()

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

call oplot_panel(panel,colorbar0,aspect,'N')
call o_set(op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Draw frame

call o_frstpt(op%fx1,op%fy1)
call o_vector(op%fx2,op%fy1)
call o_vector(op%fx2,op%fy2)
call o_vector(op%fx1,op%fy2)
call o_vector(op%fx1,op%fy1)

! Specify font # and scale font size to designated plotter coordinates

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',1.)  ! set character line width to 1
sizelab = scalelab * (op%hp2 - op%hp1)
call o_sflush()

! Write x axis label

xlabx = .5 * (op%fx1 + op%fx2)
call o_plchhq(xlabx,op%xlaby,trim(xlab),sizelab, 0.,0.)

! Write y axis label

ylaby = .5 * (op%fy1 + op%fy2)
call o_plchhq(op%ylabx,ylaby,trim(ylab),sizelab,90.,0.)

! Scale local working window (xmin,xmax,0.,1.) 
! to plotter coordinates (op%h1,op%h2,op%vp1,op%vp2)

call o_set(op%h1,op%h2,op%vp1,op%vp2,xmin,xmax,0.,1.,1)

! Plot and label X-axis ticks

tickval = nint(xmin/xinc) * xinc
if (tickval < xmin - .001 * xinc) tickval = tickval + xinc
   
do while (tickval < xmax + .001 * xinc)
   
   if (mod(nint(tickval/xinc),labincx) == 0) then  ! Only for long ticks

      dy = .020

! Encode and plot current X tick label
 
      if (xinc * labincx >= .999) then
         write (numbr,'(i6)') nint(tickval)
      elseif (xinc * labincx >= .0999) then
         write (numbr,'(f5.1)') tickval
      elseif (xinc * labincx >= .00999) then
         write (numbr,'(f5.2)') tickval
      elseif (xinc * labincx >= .000999) then
         write (numbr,'(f6.3)') tickval
      else
         write (numbr,'(f7.4)') tickval
      endif
   
      call o_plchhq(tickval,op%xtlaby,trim(adjustl(numbr)),sizelab,0.,0.)

   else                                          ! Only for short ticks
      dy = .010
   endif
   
! Plot current X tick
   
   call o_frstpt(tickval,op%fy1)
   call o_vector(tickval,op%fy1 + dy)
   call o_frstpt(tickval,op%fy2)
   call o_vector(tickval,op%fy2 - dy)
   
   tickval = tickval + xinc
enddo

! Scale local working window (0.,1.,logymin,logymax) 
! to plotter coordinates (op%hp1,op%hp2,op%v1,op%v2)

call o_set(op%hp1,op%hp2,op%v1,op%v2,0.,1.,real(logymin),real(logymax),1)

! Plot and label Y-axis ticks

jtick = logymin
   
do while (jtick <= logymax)

   do itick = 1,9

      tickval = log10(real(itick) * 10. ** real(jtick))

      if (itick == 1) then  ! Only for long ticks

         dx = .020

! Encode and plot current Y tick label

         write (numbr,'(i3)') jtick
         numbr = adjustl(numbr)
         numbr = '10:S3:'//trim(numbr)//'        '
         call o_plchhq(op%ytlabx,tickval,trim(adjustl(numbr)),sizelab,0.,1.)
      else                                          ! Only for short ticks
         dx = .010
      endif
   
! Plot current Y tick
   
      call o_frstpt(op%fx1     ,tickval)
      call o_vector(op%fx1 + dx,tickval)
      call o_frstpt(op%fx2     ,tickval)
      call o_vector(op%fx2 - dx,tickval)
   
      if (itick == 1 .and. jtick >= logymax) exit

   enddo

jtick = jtick + 1

enddo

! Scale local working window (xmin,xmax,ymin,ymax)
!  to plotter coordinates (op%h1,op%h2,op%v1,op%v2)

call o_set(op%h1,op%h2,op%v1,op%v2,xmin,xmax,real(logymin),real(logymax),1)

! Plot values

yminlog = real(logymin)
ymaxlog = real(logymax)

call o_frstpt(xval(1),min(ymaxlog,max(yminlog,log10(yval(1)))))
do i = 2,n
   call o_vector(xval(i),min(ymaxlog,max(yminlog,log10(yval(i)))))
enddo

return
end subroutine oplot_xy2log10

!===============================================================================

subroutine oplot_xy2loglog10(panel,colorbar0,aspect,scalelab  &
   ,n,xval,yval,xlab,ylab                         &
   ,logxmin,logxmax  ,logymin,logymax      )
 
use oplot_coms, only: op
use misc_coms,  only: io6

! This routine is a substitute for NCAR Graphics routine ezxy to allow 
! control over fonts, labels, axis labels, line width, scaling, etc.
! Pass in a value of 1 for n to not plot (only draw frame and ticks)

implicit none

character(len=1), intent(in) :: panel,colorbar0
real, intent(in) :: aspect,scalelab
integer, intent(in) :: n
real, intent(in) :: xval(n),yval(n)
character(len=*), intent(in) :: xlab,ylab
integer, intent(in) :: logxmin,logxmax
integer, intent(in) :: logymin,logymax

integer :: i,itick,jtick
real :: xl,yl,dx,dy,tickval,xmargin,ymargin,sizelab,xlabx,ylaby  &
   ,xminlog,xmaxlog,yminlog,ymaxlog

character(len=20)  :: numbr

! Set plot color (black)

call o_gsplci(10)
call o_gsfaci(10)
call o_gstxci(10)
call o_sflush()

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

call oplot_panel(panel,colorbar0,aspect,'N')
call o_set(op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Draw frame

call o_frstpt(op%fx1,op%fy1)
call o_vector(op%fx2,op%fy1)
call o_vector(op%fx2,op%fy2)
call o_vector(op%fx1,op%fy2)
call o_vector(op%fx1,op%fy1)

! Specify font # and scale font size to designated plotter coordinates

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',1.)  ! set character line width to 1
sizelab = scalelab * (op%hp2 - op%hp1)
call o_sflush()

! Write x axis label

xlabx = .5 * (op%fx1 + op%fx2)
call o_plchhq(xlabx,op%xlaby,trim(xlab),sizelab, 0.,0.)

! Write y axis label

ylaby = .5 * (op%fy1 + op%fy2)
call o_plchhq(op%ylabx,ylaby,trim(ylab),sizelab,90.,0.)

! Scale local working window (xmin,xmax,0.,1.) 
! to plotter coordinates (op%h1,op%h2,op%vp1,op%vp2)

call o_set(op%h1,op%h2,op%vp1,op%vp2,real(logxmin),real(logxmax),0.,1.,1)

! Plot and label X-axis ticks

itick = logxmin
   
do while (itick <= logxmax)
   
   do jtick = 1,9
   
      tickval = log10(real(jtick) * 10. ** real(itick))

      if (jtick == 1) then  ! Only for long ticks

         dy = .020

! Encode and plot current X tick label

         write (numbr,'(i3)') itick
         numbr = adjustl(numbr)
         numbr = '10:S3:'//trim(numbr)//'        '
         call o_plchhq(tickval,op%xtlaby,trim(adjustl(numbr)),sizelab,0.,0.)
      else                                          ! Only for short ticks
         dy = .010
      endif

! Plot current X tick

      call o_frstpt(tickval,op%fy1)
      call o_vector(tickval,op%fy1 + dy)
      call o_frstpt(tickval,op%fy2)
      call o_vector(tickval,op%fy2 - dy)
   
      if (jtick == 1 .and. itick >= logxmax) exit

   enddo

   itick = itick + 1

enddo

! Scale local working window (0.,1.,logymin,logymax) 
! to plotter coordinates (op%hp1,op%hp2,op%v1,op%v2)

call o_set(op%hp1,op%hp2,op%v1,op%v2,0.,1.,real(logymin),real(logymax),1)

! Plot and label Y-axis ticks

jtick = logymin
   
do while (jtick <= logymax)
   
   do itick = 1,9
   
      tickval = log10(real(itick) * 10. ** real(jtick))

      if (itick == 1) then  ! Only for long ticks

         dx = .020
 
! Encode and plot current Y tick label

         write (numbr,'(i3)') jtick
         numbr = adjustl(numbr)
         numbr = '10:S3:'//trim(numbr)//'        '
         call o_plchhq(op%ytlabx,tickval,trim(adjustl(numbr)),sizelab,0.,1.)
      else                                          ! Only for short ticks
         dx = .010
      endif
   
! Plot current Y tick
   
      call o_frstpt(op%fx1     ,tickval)
      call o_vector(op%fx1 + dx,tickval)
      call o_frstpt(op%fx2     ,tickval)
      call o_vector(op%fx2 - dx,tickval)
   
      if (itick == 1 .and. jtick >= logymax) exit

   enddo

   jtick = jtick + 1

enddo

! Scale local working window (xmin,xmax,ymin,ymax)
!  to plotter coordinates (op%h1,op%h2,op%v1,op%v2)

call o_set(op%h1,op%h2,op%v1,op%v2,real(logxmin),real(logxmax)  &
           ,real(logymin),real(logymax),1)

! Plot values

xminlog = real(logxmin)
xmaxlog = real(logxmax)

yminlog = real(logymin)
ymaxlog = real(logymax)

xl = min(xmaxlog,max(xminlog,log10(xval(1))))
yl = min(ymaxlog,max(yminlog,log10(yval(1))))

call o_frstpt(xl,yl)
do i = 2,n
   xl = min(xmaxlog,max(xminlog,log10(xval(i))))
   yl = min(ymaxlog,max(yminlog,log10(yval(i))))

   call o_vector(xl,yl)
enddo

return
end subroutine oplot_xy2loglog10

!===============================================================================

subroutine niceinc20(x1,x2,xinc,labincx)

! Use from 20 to 50 small increments, and 4 to 10 large (labeled) increments

implicit none

real, intent(in) :: x1,x2
real, intent(out) :: xinc
integer, intent(out) :: labincx

integer :: j
real :: xdif,d

xdif = max(1.e-30,abs(x2-x1))
j = int(log10(xdif))

if (xdif < 1.) then
   d = xdif * 10. ** (-j+1)  ! between 1. and 10.
else
   d = xdif * 10. ** (-j)
endif

if (d < 2.0)then
   xinc = .04 * xdif / d
elseif (d < 4.0)then
   xinc = .1 * xdif / d
else
   xinc = .2 * xdif / d
endif

labincx = 5

return
end subroutine niceinc20




