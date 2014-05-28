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
Module obnd

Contains

subroutine trsets()

use var_tables, only: num_scalar, scalar_tab
use mem_basic,  only: theta, sh_v
use micro_coms, only: level
use misc_coms,  only: io6

implicit none

integer :: n

! SCALARS - Lateral, top, and bottom boundary conditions.

! Prognostic scalars other than thil

do n = 2, num_scalar

   call latsett(scalar_tab(n)%var_p)
   call topbots(scalar_tab(n)%var_p)

enddo

! THETA and SH_V

call latsett(theta)
call topbots(theta)

if (level >= 1) then
   call latsett(sh_v)
   call topbots(sh_v)
endif

return
end subroutine trsets

!===============================================================================

subroutine latsett(sclr)

use mem_ijtabs, only: jtab_w, itab_w, istp, mrl_endl, jtw_lbcp
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: io6

implicit none

real, intent(inout) :: sclr(mza,mwa)

integer :: j,iw,iwp,k,mrl

! LBC for scalars (usually cyclic)

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
      sclr(k,iw) = sclr(k,iwp)
   enddo
enddo
endif
call rsub('W',31)

return
end subroutine latsett

!===============================================================================

subroutine topbots(sclr)

use mem_ijtabs, only: jtab_w, istp, mrl_endl, jtw_prog
use mem_grid,   only: mza, mwa, lpw

!$ use omp_lib

implicit none

real, intent(inout) :: sclr(mza,mwa)

integer :: iw,j,k,ka,mrl

! Top/bottom boundary condition for scalars: zero-gradient

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw,ka,k)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   sclr(mza,iw) = sclr(mza-1,iw)

   ka = lpw(iw)
   do k = ka-1,1,-1
      sclr(k,iw) = sclr(ka,iw)
   enddo
   
enddo
!$omp end parallel do
endif
call rsub('W',33)

return
end subroutine topbots

!===============================================================================

subroutine lbcopy_m(mrl, a1)

use mem_ijtabs, only: jtab_m, itab_m, jtm_lbcp
use mem_grid,   only: mza, mma
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1(mza,mma)

integer :: j,im,imp

! Lateral boundary copy (usually cyclic)

!----------------------------------------------------------------------
if (mrl > 0) then
do j = 1,jtab_m(jtm_lbcp)%jend(mrl); im = jtab_m(jtm_lbcp)%im(j)
   imp = itab_m(im)%imp
!----------------------------------------------------------------------

   if (present(a1)) a1(:,im) = a1(:,imp) 

enddo
endif

return
end subroutine lbcopy_m

!===============================================================================

subroutine lbcopy_u(mrl, a1, a2)

use mem_ijtabs, only: jtab_u, itab_u, jtu_lbcp
use mem_grid,   only: mza, mua
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1(mza,mua)
real, optional, intent(inout) :: a2(mza,mua)

integer :: j,iu,iup

! Lateral boundary copy (usually cyclic)

!----------------------------------------------------------------------
if (mrl > 0) then
do j = 1,jtab_u(jtu_lbcp)%jend(mrl); iu = jtab_u(jtu_lbcp)%iu(j)
   iup = itab_u(iu)%iup
!----------------------------------------------------------------------

   if (present(a1)) a1(:,iu) = a1(:,iup) 
   if (present(a2)) a2(:,iu) = a2(:,iup) 

enddo
endif

return
end subroutine lbcopy_u

!===============================================================================

subroutine lbcopy_v(mrl, vmc, vc)

use mem_ijtabs, only: jtab_v, itab_v, jtv_lbcp
use mem_grid,   only: mza, mva
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: vmc(mza,mva)
real, optional, intent(inout) :: vc (mza,mva)

integer :: j,iv,ivp

! Lateral boundary copy (usually cyclic)

!----------------------------------------------------------------------
if (mrl > 0) then
do j = 1,jtab_v(jtv_lbcp)%jend(mrl); iv = jtab_v(jtv_lbcp)%iv(j)
   ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------

   if (present(vmc)) vmc(:,iv) = vmc(:,ivp) 
   if (present(vc))  vc (:,iv) = vc (:,ivp) 

enddo
endif

return
end subroutine lbcopy_v

!===============================================================================

subroutine lbcopy_w(mrl, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, &
                         d1, d2)

use mem_ijtabs, only: jtab_w, itab_w, jtw_lbcp
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1 (mza,mwa)
real, optional, intent(inout) :: a2 (mza,mwa)
real, optional, intent(inout) :: a3 (mza,mwa)
real, optional, intent(inout) :: a4 (mza,mwa)
real, optional, intent(inout) :: a5 (mza,mwa)
real, optional, intent(inout) :: a6 (mza,mwa)
real, optional, intent(inout) :: a7 (mza,mwa)
real, optional, intent(inout) :: a8 (mza,mwa)
real, optional, intent(inout) :: a9 (mza,mwa)
real, optional, intent(inout) :: a10(mza,mwa)
real, optional, intent(inout) :: a11(mza,mwa)
real, optional, intent(inout) :: a12(mza,mwa)

real(8), optional, intent(inout) :: d1(mza,mwa)
real(8), optional, intent(inout) :: d2(mza,mwa)

integer :: j,iw,iwp

! Lateral boundary copy (usually cyclic)

!----------------------------------------------------------------------
if (mrl > 0) then
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

   if (present(a1 )) a1 (:,iw) = a1 (:,iwp) 
   if (present(a2 )) a2 (:,iw) = a2 (:,iwp) 
   if (present(a3 )) a3 (:,iw) = a3 (:,iwp) 
   if (present(a4 )) a4 (:,iw) = a4 (:,iwp) 
   if (present(a5 )) a5 (:,iw) = a5 (:,iwp) 
   if (present(a6 )) a6 (:,iw) = a6 (:,iwp) 
   if (present(a7 )) a7 (:,iw) = a7 (:,iwp) 
   if (present(a8 )) a8 (:,iw) = a8 (:,iwp) 
   if (present(a9 )) a9 (:,iw) = a9 (:,iwp) 
   if (present(a10)) a10(:,iw) = a10(:,iwp) 
   if (present(a11)) a11(:,iw) = a11(:,iwp) 
   if (present(a12)) a12(:,iw) = a12(:,iwp) 
   if (present(d1 )) d1 (:,iw) = d1 (:,iwp) 
   if (present(d2 )) d2 (:,iw) = d2 (:,iwp) 

enddo
endif

return
end subroutine lbcopy_w

!===============================================================================

subroutine lbcopy_w1d(mrl, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, d1)

use mem_ijtabs, only: jtab_w, itab_w, jtw_lbcp
use mem_grid,   only: mwa
use misc_coms,  only: io6
use consts_coms, only: r8

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1 (mwa)
real, optional, intent(inout) :: a2 (mwa)
real, optional, intent(inout) :: a3 (mwa)
real, optional, intent(inout) :: a4 (mwa)
real, optional, intent(inout) :: a5 (mwa)
real, optional, intent(inout) :: a6 (mwa)
real, optional, intent(inout) :: a7 (mwa)
real, optional, intent(inout) :: a8 (mwa)
real, optional, intent(inout) :: a9 (mwa)
real, optional, intent(inout) :: a10(mwa)
real, optional, intent(inout) :: a11(mwa)
real, optional, intent(inout) :: a12(mwa)

real(r8), optional, intent(inout) :: d1(mwa)

integer :: j,iw,iwp

! Lateral boundary copy (usually cyclic)

!----------------------------------------------------------------------
if (mrl > 0) then
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

   if (present(a1 )) a1 (iw) = a1 (iwp) 
   if (present(a2 )) a2 (iw) = a2 (iwp) 
   if (present(a3 )) a3 (iw) = a3 (iwp) 
   if (present(a4 )) a4 (iw) = a4 (iwp) 
   if (present(a5 )) a5 (iw) = a5 (iwp) 
   if (present(a6 )) a6 (iw) = a6 (iwp) 
   if (present(a7 )) a7 (iw) = a7 (iwp) 
   if (present(a8 )) a8 (iw) = a8 (iwp) 
   if (present(a9 )) a9 (iw) = a9 (iwp) 
   if (present(a10)) a10(iw) = a10(iwp) 
   if (present(a11)) a11(iw) = a11(iwp) 
   if (present(a12)) a12(iw) = a12(iwp) 
   if (present(d1 )) d1 (iw) = d1 (iwp) 

enddo
endif

return
end subroutine lbcopy_w1d

End Module obnd

