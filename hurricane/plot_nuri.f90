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
subroutine plot_nuri(iplt)

use mem_basic,   only: press
use mem_grid,    only: mwa, lpw, glatw, glonw
use misc_coms,   only: time8
use oplot_coms,  only: op
use consts_coms, only: pio180, erad

implicit none

integer :: iplt

integer :: idat,num,nfhour,iwlp,iw
integer, save :: newcall=9
real, save, dimension(0:300) :: flat,flon
real :: zeh,reh,xeh,yeh,press_lowest,bsize
character(len=2) :: title

integer, save :: iplottime=0,iplottimex=0,jplottime
real, save, dimension(1000) :: xobs,yobs,xmodel,ymodel
 
if (iplt > 1) return

if (newcall /= 1) then
   newcall = 1
   open(32,file='nuri_location_19aug2008_00UTC',status='old',form='formatted')
   do idat = 0,96
      read(32,21) num,flat(idat),flon(idat)
   enddo
   close(32)
endif
21 format(i4,2f8.2)

iplottimex = iplottimex + 1

!-------------------------------------------------------
!if (mod(iplottimex,6) /= 1) RETURN  ! SPECIAL
!if (mod(iplottimex,6) /= 1) go to 5  ! SPECIAL
!-------------------------------------------------------

iplottime = iplottime + 1

nfhour = nint(time8/3600.)

! Find "earth" coordinates of observed hurricane center

zeh = erad * sin(flat(nfhour) * pio180)
reh = erad * cos(flat(nfhour) * pio180)  ! distance from earth center
xeh = reh  * cos(flon(nfhour) * pio180)
yeh = reh  * sin(flon(nfhour) * pio180)

! The following version uses offsets to shift from the 'best track'
! initial location to the effective GFS initial location

!n zeh = erad * sin((flat(nfhour)-.095) * pio180)
!n reh = erad * cos((flat(nfhour)-.095) * pio180)  ! distance from earth center
!n xeh = reh  * cos((flon(nfhour)-.212) * pio180)
!n yeh = reh  * sin((flon(nfhour)-.212) * pio180)

! Transform hurricane earth coords to whatever projection is in use

! write(6,'(a,3i5,6e13.3)') 'pltf ',iplottime,iplt,nfhour,xeh,yeh,zeh &
!                           ,xobs(iplottime),yobs(iplottime)

call oplot_transform(iplt,xeh,yeh,zeh,xobs(iplottime),yobs(iplottime))

! Search for hurricane lowest pressure

press_lowest = 2.e5
iwlp = 1

!----------------------------------------------------------------------
do iw = 2,mwa
!----------------------------------------------------------------------

! Distance of this IW point from search center

   if (glatw(iw) < 10.  .or. glatw(iw) > 30. .or. &
       glonw(iw) < 110. .or. glonw(iw) > 150.) cycle

   if (lpw(iw) <= 3 .and. real(press(3,iw)) < press_lowest) then
      press_lowest = real(press(3,iw))
      iwlp = iw
   endif
enddo

! Find "earth" coordinates of modeled hurricane center

zeh = erad * sin(glatw(iwlp) * pio180)
reh = erad * cos(glatw(iwlp) * pio180)  ! distance from earth center
xeh = reh  * cos(glonw(iwlp) * pio180)
yeh = reh  * sin(glonw(iwlp) * pio180)

! Transform hurricane earth coords to whatever projection is in use

call oplot_transform(iplt,xeh,yeh,zeh,xmodel(iplottime),ymodel(iplottime))

5 continue

! Set character font, line width, and plotted size

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',2.)
bsize = .008 * (op%hp2 - op%hp1) * 0.3

! Set color to black

call o_gsplci(10)
call o_gstxci(10)
call o_gsfaci(10)
call o_sflush()

! Plot initial modelled hurricane eye location with black color

  write(title,'(a1)') '+'
  call o_plchhq (xmodel(1),ymodel(1),trim(adjustl(title)),bsize,0.,0.)

! Plot modelled hurricane eye location with black color, excluding initial time

do jplottime = 2,iplottime
   write(title,'(i2)') jplottime-1
   call o_plchhq (xmodel(jplottime),ymodel(jplottime),trim(adjustl(title)) &
      ,bsize,0.,0.)
enddo

! Set color to red

call o_gsplci(1)
call o_gstxci(1)
call o_gsfaci(1)
call o_sflush()

! Plot initial observed hurricane eye location with red color

  write(title,'(a1)') '+'
  call o_plchhq (xobs(1),yobs(1),trim(adjustl(title)),bsize,0.,0.)

! Plot observed hurricane eye location with red color, excluding initial time

do jplottime = 2,iplottime
   write(title,'(i2)') jplottime-1
   call o_plchhq (xobs(jplottime),yobs(jplottime),trim(adjustl(title)),bsize,0.,0.)
enddo

end subroutine plot_nuri

