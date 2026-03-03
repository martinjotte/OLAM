subroutine plot_frances(iplt)

use mem_ijtabs
use mem_basic
use mem_cuparm
use mem_grid
use mem_micro
use mem_radiate
use mem_tend
use mem_turb
use mem_nudge

use var_tables

use misc_coms
use oplot_coms
use consts_coms

implicit none

integer :: iplt

integer :: idat,num,nfhour,iwlp,iw
integer, save :: newcall=9
real, save, dimension(0:300) :: flat,flon
real :: zef,ref,xef,yef,press_lowest,bsize
character(len=2) :: title

integer, save :: iplottime=0,iplottimex=0,jplottime
real, save, dimension(1000) :: xobs,yobs,xmodel,ymodel

if (iplt > 1) return

if (newcall /= 1) then
   newcall = 1
   open(32,file='../special/frances_location_1sept00z',status='old',form='formatted')
   do idat = 0,201
      read(32,21) num,flat(idat),flon(idat)
   enddo
   close(32)
endif
21 format(i4,2f8.2)

iplottimex = iplottimex + 1

!-------------------------------------------------------
if (mod(iplottimex,10) /= 1) RETURN  ! SPECIAL
!-------------------------------------------------------

iplottime = iplottime + 1

nfhour = nint(time8/3600.)

! Find "earth" coordinates of observed hurricane center

zef = erad * sin(flat(nfhour) * pio180)
ref = erad * cos(flat(nfhour) * pio180)  ! distance from earth center
xef = ref  * cos(flon(nfhour) * pio180)
yef = ref  * sin(flon(nfhour) * pio180)

! Transform hurricane earth coords to whatever projection is in use

! write(6,'(a,3i5,6e13.3)') 'pltf ',iplottime,iplt,nfhour,xef,yef,zef &
!                           ,xobs(iplottime),yobs(iplottime)

call oplot_transform(iplt,xef,yef,zef,flon(nfhour),flat(nfhour),xobs(iplottime),yobs(iplottime))

! Search for hurricane lowest pressure

press_lowest = 2.e5
iwlp = 1

!----------------------------------------------------------------------
do iw = 2,mwa
!----------------------------------------------------------------------

! Distance of this IW point from search center

   if (glatw(iw) < 20.  .or. glatw(iw) > 35. .or.  &
       glonw(iw) < -90. .or. glonw(iw) > -65.) cycle

   if (lpw(iw) <= 3 .and. real(press(3,iw)) < press_lowest) then
      press_lowest = real(press(3,iw))
      iwlp = iw
   endif
enddo

! Find "earth" coordinates of modeled hurricane center

   zef = erad * sin(glatw(iwlp) * pio180)
   ref = erad * cos(glatw(iwlp) * pio180)  ! distance from earth center
   xef = ref  * cos(glonw(iwlp) * pio180)
   yef = ref  * sin(glonw(iwlp) * pio180)

! Transform hurricane earth coords to whatever projection is in use

call oplot_transform(iplt,xef,yef,zef,glonw(iwlp),glatw(iwlp),xmodel(iplottime),ymodel(iplottime))

! Set character font, line width, and plotted size

call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
call o_pcsetr('CL',2.)
bsize = .008 * (op%hp2 - op%hp1)

! Set color to black

call o_sflush()
call o_gsplci(10)
call o_gstxci(10)
call o_gsfaci(10)

! Plot modelled hurricane eye location with black color, excluding initial time

do jplottime = 2,iplottime
   write(title,'(i2)') jplottime-1
   call o_plchhq (xmodel(jplottime),ymodel(jplottime),trim(adjustl(title))  &
      ,bsize,0.,0.)
enddo

! Set color to red

call o_sflush()
call o_gsplci(270)
call o_gstxci(270)
call o_gsfaci(270)

! Plot initial modelled hurricane eye location with red color

write(title,'(a1)') '+'
call o_plchhq (xmodel(1),ymodel(1),trim(adjustl(title)),bsize,0.,0.)

! Plot observed hurricane eye location with red color, excluding initial time

do jplottime = 2,iplottime
   write(title,'(i2)') jplottime-1
call o_plchhq (xobs(jplottime),yobs(jplottime),trim(adjustl(title)),bsize,0.,0.)
enddo

return
end subroutine plot_frances

