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
subroutine mdloopf(init,im,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_md, nloops_m
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: im
integer, intent(in) :: j1,j2,j3,j4,j5,j6

if (init == 'f') itab_md(im)%loop(1:nloops_m) = .false.

if (j1 < 0) itab_md(im)%loop(abs(j1)) = .false.
if (j2 < 0) itab_md(im)%loop(abs(j2)) = .false.
if (j3 < 0) itab_md(im)%loop(abs(j3)) = .false.
if (j4 < 0) itab_md(im)%loop(abs(j4)) = .false.
if (j5 < 0) itab_md(im)%loop(abs(j5)) = .false.
if (j6 < 0) itab_md(im)%loop(abs(j6)) = .false.

if (j1 > 0) itab_md(im)%loop(j1) = .true.
if (j2 > 0) itab_md(im)%loop(j2) = .true.
if (j3 > 0) itab_md(im)%loop(j3) = .true.
if (j4 > 0) itab_md(im)%loop(j4) = .true.
if (j5 > 0) itab_md(im)%loop(j5) = .true.
if (j6 > 0) itab_md(im)%loop(j6) = .true.

return
end subroutine mdloopf

!===============================================================================

subroutine udloopf(init,iu,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_ud, nloops_u
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: iu
integer, intent(in) :: j1,j2,j3,j4,j5,j6
   
if (init == 'f') itab_ud(iu)%loop(1:nloops_u) = .false.

if (j1  < 0) itab_ud(iu)%loop(abs(j1))  = .false.
if (j2  < 0) itab_ud(iu)%loop(abs(j2))  = .false.
if (j3  < 0) itab_ud(iu)%loop(abs(j3))  = .false.
if (j4  < 0) itab_ud(iu)%loop(abs(j4))  = .false.
if (j5  < 0) itab_ud(iu)%loop(abs(j5))  = .false.
if (j6  < 0) itab_ud(iu)%loop(abs(j6))  = .false.

if (j1  > 0) itab_ud(iu)%loop(j1)  = .true.
if (j2  > 0) itab_ud(iu)%loop(j2)  = .true.
if (j3  > 0) itab_ud(iu)%loop(j3)  = .true.
if (j4  > 0) itab_ud(iu)%loop(j4)  = .true.
if (j5  > 0) itab_ud(iu)%loop(j5)  = .true.
if (j6  > 0) itab_ud(iu)%loop(j6)  = .true.

return
end subroutine udloopf

!===============================================================================

subroutine wdloopf(init,iw,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_wd, nloops_w
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: iw
integer, intent(in) :: j1,j2,j3,j4,j5,j6
                             
if (init == 'f') itab_wd(iw)%loop(1:nloops_w) = .false.

if (j1  < 0) itab_wd(iw)%loop(abs(j1))  = .false.
if (j2  < 0) itab_wd(iw)%loop(abs(j2))  = .false.
if (j3  < 0) itab_wd(iw)%loop(abs(j3))  = .false.
if (j4  < 0) itab_wd(iw)%loop(abs(j4))  = .false.
if (j5  < 0) itab_wd(iw)%loop(abs(j5))  = .false.
if (j6  < 0) itab_wd(iw)%loop(abs(j6))  = .false.

if (j1  > 0) itab_wd(iw)%loop(j1)  = .true.
if (j2  > 0) itab_wd(iw)%loop(j2)  = .true.
if (j3  > 0) itab_wd(iw)%loop(j3)  = .true.
if (j4  > 0) itab_wd(iw)%loop(j4)  = .true.
if (j5  > 0) itab_wd(iw)%loop(j5)  = .true.
if (j6  > 0) itab_wd(iw)%loop(j6)  = .true.

return
end subroutine wdloopf

!===============================================================================

subroutine mloopf(init,im,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_m, nloops_m
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: im
integer, intent(in) :: j1,j2,j3,j4,j5,j6

if (init == 'f') itab_m(im)%loop(1:nloops_m) = .false.

if (j1 < 0) itab_m(im)%loop(abs(j1)) = .false.
if (j2 < 0) itab_m(im)%loop(abs(j2)) = .false.
if (j3 < 0) itab_m(im)%loop(abs(j3)) = .false.
if (j4 < 0) itab_m(im)%loop(abs(j4)) = .false.
if (j5 < 0) itab_m(im)%loop(abs(j5)) = .false.
if (j6 < 0) itab_m(im)%loop(abs(j6)) = .false.

if (j1 > 0) itab_m(im)%loop(j1) = .true.
if (j2 > 0) itab_m(im)%loop(j2) = .true.
if (j3 > 0) itab_m(im)%loop(j3) = .true.
if (j4 > 0) itab_m(im)%loop(j4) = .true.
if (j5 > 0) itab_m(im)%loop(j5) = .true.
if (j6 > 0) itab_m(im)%loop(j6) = .true.

return
end subroutine mloopf

!===============================================================================

subroutine uloopf(init,iu,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_u, nloops_u
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: iu
integer, intent(in) :: j1,j2,j3,j4,j5,j6
   
if (init == 'f') itab_u(iu)%loop(1:nloops_u) = .false.

if (j1  < 0) itab_u(iu)%loop(abs(j1))  = .false.
if (j2  < 0) itab_u(iu)%loop(abs(j2))  = .false.
if (j3  < 0) itab_u(iu)%loop(abs(j3))  = .false.
if (j4  < 0) itab_u(iu)%loop(abs(j4))  = .false.
if (j5  < 0) itab_u(iu)%loop(abs(j5))  = .false.
if (j6  < 0) itab_u(iu)%loop(abs(j6))  = .false.

if (j1  > 0) itab_u(iu)%loop(j1)  = .true.
if (j2  > 0) itab_u(iu)%loop(j2)  = .true.
if (j3  > 0) itab_u(iu)%loop(j3)  = .true.
if (j4  > 0) itab_u(iu)%loop(j4)  = .true.
if (j5  > 0) itab_u(iu)%loop(j5)  = .true.
if (j6  > 0) itab_u(iu)%loop(j6)  = .true.

return
end subroutine uloopf

!===============================================================================

subroutine vloopf(init,iv,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_v, nloops_v
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: iv
integer, intent(in) :: j1,j2,j3,j4,j5,j6
   
if (init == 'f') itab_v(iv)%loop(1:nloops_v) = .false.

if (j1  < 0) itab_v(iv)%loop(abs(j1))  = .false.
if (j2  < 0) itab_v(iv)%loop(abs(j2))  = .false.
if (j3  < 0) itab_v(iv)%loop(abs(j3))  = .false.
if (j4  < 0) itab_v(iv)%loop(abs(j4))  = .false.
if (j5  < 0) itab_v(iv)%loop(abs(j5))  = .false.
if (j6  < 0) itab_v(iv)%loop(abs(j6))  = .false.

if (j1  > 0) itab_v(iv)%loop(j1)  = .true.
if (j2  > 0) itab_v(iv)%loop(j2)  = .true.
if (j3  > 0) itab_v(iv)%loop(j3)  = .true.
if (j4  > 0) itab_v(iv)%loop(j4)  = .true.
if (j5  > 0) itab_v(iv)%loop(j5)  = .true.
if (j6  > 0) itab_v(iv)%loop(j6)  = .true.

return
end subroutine vloopf

!===============================================================================

subroutine wloopf(init,iw,j1,j2,j3,j4,j5,j6)

use mem_ijtabs, only: itab_w, nloops_w
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: init

integer, intent(in) :: iw
integer, intent(in) :: j1,j2,j3,j4,j5,j6
                             
if (init == 'f') itab_w(iw)%loop(1:nloops_w) = .false.

if (j1  < 0) itab_w(iw)%loop(abs(j1))  = .false.
if (j2  < 0) itab_w(iw)%loop(abs(j2))  = .false.
if (j3  < 0) itab_w(iw)%loop(abs(j3))  = .false.
if (j4  < 0) itab_w(iw)%loop(abs(j4))  = .false.
if (j5  < 0) itab_w(iw)%loop(abs(j5))  = .false.
if (j6  < 0) itab_w(iw)%loop(abs(j6))  = .false.

if (j1  > 0) itab_w(iw)%loop(j1)  = .true.
if (j2  > 0) itab_w(iw)%loop(j2)  = .true.
if (j3  > 0) itab_w(iw)%loop(j3)  = .true.
if (j4  > 0) itab_w(iw)%loop(j4)  = .true.
if (j5  > 0) itab_w(iw)%loop(j5)  = .true.
if (j6  > 0) itab_w(iw)%loop(j6)  = .true.

return
end subroutine wloopf
