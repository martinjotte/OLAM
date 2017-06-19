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

! The subroutines in this file copy real quantities to real*4 quantities 
! and call NCAR Graphics subroutines with the latter.

subroutine o_opngks()
  use misc_coms, only: io6
  implicit none

  write(io6,*)
  write(io6,*) "You are not using the NCAR Graphics library."
  write(io6,*) "No plots will be created."
  write(io6,*)

end subroutine o_opngks

!===============================================================================

subroutine o_clsgks()
  implicit none
end subroutine o_clsgks

!===============================================================================

subroutine o_frame()
  implicit none
end subroutine o_frame

!===============================================================================

subroutine o_frstpt(o_x,o_y)
  implicit none
  real, intent(in) :: o_x,o_y
end subroutine o_frstpt

!===============================================================================

subroutine o_vector(o_x,o_y)
  implicit none
  real, intent(in) :: o_x,o_y
end subroutine o_vector

!===============================================================================

subroutine o_mappos(o_x1,o_x2,o_y1,o_y2)
  implicit none
  real, intent(in) :: o_x1,o_x2,o_y1,o_y2
end subroutine o_mappos

!===============================================================================

subroutine o_maproj(chr,o_x,o_y,o_r)
  implicit none
  character(*), intent(in) :: chr
  real,         intent(in) :: o_x,o_y,o_r
end subroutine o_maproj

!===============================================================================

subroutine o_mapset(chr,o_x1,o_x2,o_y1,o_y2)
  implicit none
  character(*), intent(in) :: chr
  real,         intent(in) :: o_x1,o_x2,o_y1,o_y2
end subroutine o_mapset

!===============================================================================

subroutine o_set(o_x1,o_x2,o_y1,o_y2,o_fx1,o_fx2,o_fy1,o_fy2,i)
  implicit none
  real,    intent(in) :: o_x1,o_x2,o_y1,o_y2,o_fx1,o_fx2,o_fy1,o_fy2
  integer, intent(in) :: i
end subroutine o_set

!===============================================================================

subroutine o_plchhq(o_x,o_y,chr,o_r1,o_r2,o_r3)
  implicit none
  real,         intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
  character(*), intent(in) :: chr
end subroutine o_plchhq

!===============================================================================

subroutine o_plchmq(o_x,o_y,chr,o_r1,o_r2,o_r3)
  implicit none
  real,         intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
  character(*), intent(in) :: chr
end subroutine o_plchmq

!===============================================================================

subroutine o_plchlq(o_x,o_y,chr,o_r1,o_r2,o_r3)
  implicit none
  real,         intent(in) :: o_x,o_y,o_r1,o_r2,o_r3
  character(*), intent(in) :: chr
end subroutine o_plchlq

!===============================================================================

subroutine o_sfsgfa(o_x,o_y,nr,icolor)
  implicit none
  real,    intent(in) :: o_x(*),o_y(*)
  integer, intent(in) :: nr,icolor
end subroutine o_sfsgfa

!===============================================================================

subroutine o_hls(iwk,ic,o_h,o_l,o_s)
  implicit none
  integer, intent(in) :: iwk,ic
  real,    intent(in) :: o_h,o_l,o_s
end subroutine o_hls

!===============================================================================

subroutine o_gscr(iwk,ic,o_r,o_g,o_b)
  implicit none
  integer, intent(in) :: iwk,ic
  real,    intent(in) :: o_r,o_g,o_b
end subroutine o_gscr

!===============================================================================

subroutine o_gsplci(i)
  implicit none
  integer, intent(in) :: i
end subroutine o_gsplci

!===============================================================================

subroutine o_gstxci(i)
  implicit none
  integer, intent(in) :: i
end subroutine o_gstxci

!===============================================================================

subroutine o_gsfaci(i)
  implicit none
  integer, intent(in) :: i
end subroutine o_gsfaci

!===============================================================================

subroutine o_gsclip(i)
  implicit none
  integer, intent(in) :: i
end subroutine o_gsclip

!===============================================================================

subroutine o_gsasf(i)
  implicit none
  integer, intent(in) :: i(13)
end subroutine o_gsasf

!===============================================================================

subroutine o_gsfais(i)
  implicit none
  integer, intent(in) :: i
end subroutine o_gsfais

!===============================================================================

subroutine o_sfseti(a,i)
  implicit none
  integer,      intent(in) :: i
  character(*), intent(in) :: a
end subroutine o_sfseti

!===============================================================================

subroutine o_sflush()
  implicit none
end subroutine o_sflush

!===============================================================================

subroutine o_mapint()
  implicit none
end subroutine o_mapint

!===============================================================================

subroutine o_mapsti(a,i)
  implicit none
  character(*), intent(in) :: a
  integer,      intent(in) :: i
end subroutine o_mapsti

!===============================================================================

subroutine o_mapstc(a,b)
  implicit none
  character(*), intent(in) :: a,b
end subroutine o_mapstc

!===============================================================================

subroutine o_maplot()
  implicit none
end subroutine o_maplot

!===============================================================================

subroutine o_pcsetr(a,r)
  implicit none
  character(*), intent(in) :: a
  real,         intent(in) :: r
end subroutine o_pcsetr

!===============================================================================

subroutine o_pcseti(a,i)
  implicit none
  character(*), intent(in) :: a
  integer,      intent(in) :: i
end subroutine o_pcseti

!===============================================================================

subroutine o_clswk()
  implicit none
end subroutine o_clswk

!===============================================================================

subroutine o_reopnwk()
  implicit none
end subroutine o_reopnwk

!===============================================================================
