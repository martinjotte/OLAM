! The wrapper subroutines in this file call NCAR Graphics library routines.
! Use the "dummy" wrappers in o_ncar_dummy.f90 instead if you do not have
! the NCAR graphics library.

!===============================================================================

subroutine o_opngks()

  use oplot_coms, only: op
  use mem_para,   only: myrank

  implicit none

  integer, parameter :: ierrf=6, lunit=2, iwkid=1
  integer            :: iwtype, iopen

  if (myrank > 0) return

  ! Instead of just calling opngks, we use the NCAR GKS
  ! routines to setup the plotting types and devices
  ! call opngks()

  ! Default to ncar graphics meta file

  iwtype = 1

  if (op%plttype == 1) then

     ! Set ncar graphics postscript device

     if (op%pltorient == 0) then
        iwtype = 20
     else
        iwtype = 26
     endif

     op%pltname = trim(op%pltname)//".ps"

  elseif (op%plttype == 2) then

     ! Set ncar graphics pdf device

     if (op%pltorient == 0) then
        iwtype = 11
     else
        iwtype = 12
     endif

     op%pltname = trim(op%pltname)//".pdf"

  endif

  ! Open GKS if it is closed

  call gqops(iopen)
  if (iopen == 0) call gopks (ierrf, 0)

  ! Set output file and device

  call ngsetc('ME', op%pltname)
  call gopwk (iwkid, lunit, iwtype)
  call gacwk (iwkid)

end subroutine o_opngks

!===============================================================================

subroutine o_clsgks()

  use mem_para, only: myrank

  implicit none

  integer, parameter :: iwkid=1
  integer            :: iopen
  integer, external  :: ngckop

  if (myrank > 0) return

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

end subroutine o_frame

!===============================================================================

subroutine o_frstpt(x,y)

  implicit none

  real, intent(in) :: x,y

  call frstpt(x,y)

end subroutine o_frstpt

!===============================================================================

subroutine o_vector(x,y)

  implicit none

  real, intent(in) :: x,y

  call vector(x,y)

end subroutine o_vector

!===============================================================================

subroutine o_mappos(x1,x2,y1,y2)

  implicit none

  real, intent(in) :: x1,x2,y1,y2

  call mappos(x1,x2,y1,y2)

end subroutine o_mappos

!===============================================================================

subroutine o_maproj(chr,x,y,r)

  implicit none

  character(*), intent(in) :: chr
  real,         intent(in) :: x,y,r

  call maproj(chr,x,y,r)

end subroutine o_maproj

!===============================================================================

subroutine o_mapset(chr,x1,x2,y1,y2)

  implicit none

  character(*), intent(in) :: chr
  real,         intent(in) :: x1,x2,y1,y2

  call mapset(chr,x1,x2,y1,y2)

end subroutine o_mapset

!===============================================================================

subroutine o_set(x1,x2,y1,y2,fx1,fx2,fy1,fy2,i)

  implicit none

  real,    intent(in) :: x1,x2,y1,y2,fx1,fx2,fy1,fy2
  integer, intent(in) :: i

  call set(x1,x2,y1,y2,fx1,fx2,fy1,fy2,i)

end subroutine o_set

!===============================================================================

subroutine o_plchhq(x,y,chr,r1,r2,r3)

  implicit none

  real,         intent(in) :: x,y,r1,r2,r3
  character(*), intent(in) :: chr

  call plchhq(x,y,chr,r1,r2,r3)

end subroutine o_plchhq

!===============================================================================

subroutine o_plchmq(x,y,chr,r1,r2,r3)

  implicit none

  real,         intent(in) :: x,y,r1,r2,r3
  character(*), intent(in) :: chr

  call plchmq(x,y,chr,r1,r2,r3)

end subroutine o_plchmq

!===============================================================================

subroutine o_plchlq(x,y,chr,r1,r2,r3)

  implicit none

  real,         intent(in) :: x,y,r1,r2,r3
  character(*), intent(in) :: chr

  call plchlq(x,y,chr,r1,r2,r3)

end subroutine o_plchlq

!===============================================================================

subroutine o_sfsgfa(x,y,nr,icolor)

  implicit none

  real,    intent(in) :: x(nr), y(nr)
  integer, intent(in) :: nr, icolor
! integer             :: ind(21)
! real                :: dst(14)
!
!  call sfsgfa(x,y,nr,dst,14,ind,21,icolor)

  call gsfaci(icolor)
  call gfa(nr,x,y)

end subroutine o_sfsgfa

!===============================================================================

subroutine o_hls(iwk,ic,h,l,s)

  implicit none

  integer, intent(in) :: iwk,ic
  real,    intent(in) :: h,l,s
  real                :: r,g,b

  call hlsrgb (h,l,s,r,g,b)
  call gscr (iwk,ic,r,g,b)

end subroutine o_hls

!===============================================================================

subroutine o_gscr(iwk,ic,r,g,b)

  implicit none

  integer, intent(in) :: iwk, ic
  real,    intent(in) :: r, g, b

  call gscr (iwk, ic, r, g, b)

end subroutine o_gscr

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

subroutine o_gslwsc(x)

  implicit none

  real, intent(in) :: x

  call gslwsc(x)

end subroutine o_gslwsc

!===============================================================================

subroutine o_sfseti(a,i)

  implicit none

  integer,      intent(in) :: i
  character(*), intent(in) :: a

  call sfseti(a,i)

end subroutine o_sfseti

!===============================================================================

subroutine o_sflush()

  implicit none

! Flushes line buffers, updates worksatations, and flushes all I/O buffers:
! call sflush()

! Just flushing the line buffers is probably all that is necessary:
  call plotif(0.,0.,2)

end subroutine o_sflush

!===============================================================================

subroutine o_mapint()

  implicit none

  call mapint()

end subroutine o_mapint

!===============================================================================

subroutine o_mapgrd()

  implicit none

  call mapgrd()

end subroutine o_mapgrd

!===============================================================================

subroutine o_mapstr(a,r)

  implicit none

  character(*), intent(in) :: a
  real,         intent(in) :: r

  call mapstr(a,r)

end subroutine o_mapstr

!===============================================================================

subroutine o_mapsti(a,i)

  implicit none

  character(*), intent(in) :: a
  integer,      intent(in) :: i

  call mapsti(a,i)

end subroutine o_mapsti

!===============================================================================

subroutine o_mapstc(a,b)

  implicit none

  character(*), intent(in) :: a,b

  call mapstc(a,b)

end subroutine o_mapstc

!===============================================================================

subroutine o_maplot()

  implicit none

  call maplot()

end subroutine o_maplot

!===============================================================================

subroutine o_mplndr(a,i)

  implicit none

  character(*), intent(in) :: a
  integer,      intent(in) :: i

  call mplndr(a,i)

end subroutine o_mplndr

!===============================================================================

subroutine o_pcsetr(a,r)

  implicit none

  character(*), intent(in) :: a
  real,         intent(in) :: r

  call pcsetr(a,r)

end subroutine o_pcsetr

!===============================================================================

subroutine o_pcseti(a,i)

  implicit none

  character(*), intent(in) :: a
  integer,      intent(in) :: i

  call pcseti(a,i)

end subroutine o_pcseti

!===============================================================================

subroutine o_gngpat(a,b,i)

  implicit none

  character(*), intent( in) :: a, b
  integer,      intent(out) :: i

  call gngpat(a,b,i)

end subroutine o_gngpat

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
  use mem_para,   only: myrank

  implicit none

  integer, parameter :: iwkid=1
  integer            :: iopen
  integer, external  :: ngckop

  if (myrank > 0) return

  iopen = ngckop(iwkid)

  if (iopen .ne. 1) then

     call ngreop(iwkid, 2, 1, op%pltname, 1, op%i_att, op%r_att, 0, 0, 0)
     call gacwk(iwkid)
     call gks_colors(iwkid)

  endif

end subroutine o_reopnwk

!===============================================================================
