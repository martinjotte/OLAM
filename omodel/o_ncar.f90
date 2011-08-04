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

! The subroutines in this file copy real quantities to real*4 quantities 
! and call NCAR Graphics subroutines with the latter.

subroutine o_opngks()

use oplot_coms, only: op
use misc_coms,  only: iparallel
use mem_para,   only: myrank

implicit none

integer, parameter :: ierrf=6, lunit=2, iwkid=1
integer            :: iwtype
character(len=10) :: number

! call opngks()

iwtype = 1     ! default to ncar graphics meta file

if (op%plttype == 0 .and. iparallel == 1) then

! For parallel run, append rank onto metacode filename

   write (number,'(i10)') myrank
   op%pltname = trim(op%pltname)//'_r'//trim(adjustl(number))

elseif (op%plttype == 1) then

! set ncar graphics postscript device

   if (op%pltorient == 0) then
      iwtype = 20
   else
      iwtype = 26
   endif
   
   op%pltname = trim(op%pltname)//".ps"

elseif (op%plttype == 2) then

! set ncar graphics pdf device

   if (op%pltorient == 0) then
      iwtype = 11
   else
      iwtype = 12
   endif
   
   op%pltname = trim(op%pltname)//".pdf"

endif

call gopks (ierrf, 0)
call ngsetc('ME', op%pltname)
call gopwk (iwkid, lunit, iwtype)
call gacwk (iwkid)

end subroutine o_opngks

!===============================================================================

subroutine o_clsgks()

implicit none

integer, parameter :: iwkid=1
integer            :: iopen
integer, external  :: ngckop

! call clsgks()

! Check if GKS workstation is active, and if so
! deactivate and close it.

iopen = ngckop(iwkid)
if (iopen .eq. 1) then
   call gdawk(iwkid)
   call gclwk(iwkid)
endif

! Close GKS if it is open

call gqops(iopen)
if (iopen > 0) call gclks()

end subroutine o_clsgks

!===============================================================================

subroutine o_frame()

use oplot_coms, only: op

implicit none

call frame()

op%iplotback = 0

return
end subroutine o_frame

!===============================================================================

subroutine o_frstpt(o_x,o_y)

implicit none

real, intent(in) :: o_x,o_y

real(kind=4) :: x,y

x = o_x
y = o_y

call frstpt(x,y)

return
end subroutine o_frstpt

!===============================================================================

subroutine o_vector(o_x,o_y)

implicit none

real, intent(in) :: o_x,o_y

real(kind=4) :: x,y

x = o_x
y = o_y

call vector(x,y)

return
end subroutine o_vector

!===============================================================================

subroutine o_mappos(o_x1,o_x2,o_y1,o_y2)

implicit none

real, intent(in) :: o_x1,o_x2,o_y1,o_y2

real(kind=4) :: x1,x2,y1,y2

x1 = o_x1
x2 = o_x2
y1 = o_y1
y2 = o_y2

call mappos(x1,x2,y1,y2)

return
end subroutine o_mappos

!===============================================================================

subroutine o_maproj(chr,o_x,o_y,o_r)

implicit none

character(len=*), intent(in) :: chr
real, intent(in) :: o_x,o_y,o_r

real(kind=4) :: x,y,r

x = o_x
y = o_y
r = o_r

call maproj(chr,x,y,r)

return
end subroutine o_maproj

!===============================================================================

subroutine o_mapset(chr,o_x1,o_x2,o_y1,o_y2)

implicit none

character(len=*), intent(in) :: chr
real, intent(in) :: o_x1,o_x2,o_y1,o_y2

real(kind=4) :: x1,x2,y1,y2

x1 = o_x1
x2 = o_x2
y1 = o_y1
y2 = o_y2

call mapset(chr,x1,x2,y1,y2)

return
end subroutine o_mapset

!===============================================================================

subroutine o_set(o_x1,o_x2,o_y1,o_y2,o_fx1,o_fx2,o_fy1,o_fy2,i)

implicit none

real, intent(in) :: o_x1,o_x2,o_y1,o_y2,o_fx1,o_fx2,o_fy1,o_fy2
integer, intent(in) :: i

real(kind=4) :: x1,x2,y1,y2,fx1,fx2,fy1,fy2

x1 = o_x1
x2 = o_x2
y1 = o_y1
y2 = o_y2
fx1 = o_fx1
fx2 = o_fx2
fy1 = o_fy1
fy2 = o_fy2

call set(x1,x2,y1,y2,fx1,fx2,fy1,fy2,i)

return
end subroutine o_set

!===============================================================================

subroutine o_plchhq(o_x,o_y,chr,o_r1,o_r2,o_r3)

implicit none

real, intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
character(len=*), intent(in) :: chr

real(kind=4) :: x,y,r1,r2,r3

x = o_x
y = o_y
r1 = o_r1
r2 = o_r2
r3 = o_r3

call plchhq(x,y,chr,r1,r2,r3)

return
end subroutine o_plchhq

!===============================================================================

subroutine o_plchmq(o_x,o_y,chr,o_r1,o_r2,o_r3)

implicit none

real, intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
character(len=*), intent(in) :: chr

real(kind=4) :: x,y,r1,r2,r3

x = o_x
y = o_y
r1 = o_r1
r2 = o_r2
r3 = o_r3

call plchmq(x,y,chr,r1,r2,r3)

return
end subroutine o_plchmq

!===============================================================================

subroutine o_plchlq(o_x,o_y,chr,o_r1,o_r2,o_r3)

implicit none

real, intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
character(len=*), intent(in) :: chr

real(kind=4) :: x,y,r1,r2,r3

x = o_x
y = o_y
r1 = o_r1
r2 = o_r2
r3 = o_r3

call plchlq(x,y,chr,r1,r2,r3)

return
end subroutine o_plchlq

!===============================================================================

subroutine o_sfsgfa(o_x,o_y,nr,icolor)

implicit none

real, intent(in) :: o_x(*),o_y(*)
integer, intent(in) :: nr,icolor

integer :: ind(18)
real(kind=4) :: dst(12)
real(kind=4), dimension(nr) :: x,y

x(1:nr) = o_x(1:nr)
y(1:nr) = o_y(1:nr)

call sfsgfa (x,y,nr,dst,12,ind,18,icolor)

return
end subroutine o_sfsgfa

!===============================================================================

subroutine o_hls(iwk,ic,o_h,o_l,o_s)

implicit none

integer, intent(in) :: iwk,ic
real, intent(in) :: o_h,o_l,o_s

real(kind=4) :: h,l,s
real(kind=4) :: r,g,b

h = o_h
l = o_l
s = o_s

call hlsrgb (h,l,s,r,g,b)
call gscr (iwk,ic,r,g,b)

return
end subroutine o_hls

!===============================================================================

subroutine o_gsplci(i)

implicit none

integer, intent(in) :: i

call gsplci(i)

end subroutine o_gsplci

!===============================================================================

subroutine o_gstxci(i)

implicit none

integer, intent(in) :: i

call gstxci(i)

end subroutine o_gstxci

!===============================================================================

subroutine o_gsfaci(i)

implicit none

integer, intent(in) :: i

call gsfaci(i)

end subroutine o_gsfaci

!===============================================================================

subroutine o_gsclip(i)

implicit none

integer, intent(in) :: i

call gsclip(i)

end subroutine o_gsclip

!===============================================================================

subroutine o_gsasf(i)

implicit none

integer, intent(in) :: i(13)

call gsasf(i)

end subroutine o_gsasf

!===============================================================================

subroutine o_gsfais(i)

implicit none

integer, intent(in) :: i

call gsfais(i)

end subroutine o_gsfais

!===============================================================================

subroutine o_sfseti(a,i)

implicit none

integer, intent(in) :: i
character(len=*) :: a

call sfseti(a,i)

end subroutine o_sfseti

!===============================================================================

subroutine o_sflush()

implicit none

call sflush()

end subroutine o_sflush

!===============================================================================

subroutine o_mapint()

implicit none

call mapint()

end subroutine o_mapint

!===============================================================================

subroutine o_mapsti(a,i)

implicit none

character(len=*), intent(in) :: a
integer, intent(in) :: i

call mapsti(a,i)

end subroutine o_mapsti

!===============================================================================

subroutine o_mapstc(a,b)

implicit none

character(len=*), intent(in) :: a,b

call mapstc(a,b)

end subroutine o_mapstc

!===============================================================================

subroutine o_maplot()

implicit none

call maplot()

end subroutine o_maplot

!===============================================================================

subroutine o_pcsetr(a,r)

implicit none

character(len=*), intent(in) :: a
real, intent(in) :: r

call pcsetr(a,r)

end subroutine o_pcsetr

!===============================================================================

subroutine o_pcseti(a,i)

implicit none

character(len=*), intent(in) :: a
integer, intent(in) :: i

call pcseti(a,i)

end subroutine o_pcseti

!===============================================================================

subroutine o_clswk()

use oplot_coms, only: op
implicit none

integer, parameter :: iwkid=1

call ngsrat(2, op%i_att, op%r_att)
call gdawk(iwkid)
call gclwk(iwkid)

end subroutine o_clswk

!===============================================================================

subroutine o_reopnwk()

use oplot_coms, only: op
use plotcolors, only: gks_colors
implicit none

integer, parameter :: iwkid=1
integer            :: iopen
integer, external  :: ngckop

iopen = ngckop(iwkid)

if (iopen .ne. 1) then

  call ngreop(iwkid, 2, 1, op%pltname, 1, op%i_att, op%r_att, 0, 0, 0)
  call gacwk(iwkid)
  call gks_colors(iwkid)

endif

end subroutine o_reopnwk
