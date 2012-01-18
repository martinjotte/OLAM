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
Module massflux

Contains

subroutine zero_massflux(wmarwsc, rho_old, umarusc)

use mem_ijtabs, only: jtab_u, jtab_v, jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mua, mva, mwa, lpw
use mem_basic,  only: rho
use misc_coms,  only: io6

!$ use omp_lib

implicit none

real(kind=8), intent(out) :: wmarwsc(mza,mwa)
real(kind=8), intent(out) :: rho_old(mza,mwa)
real(kind=8), intent(out) :: umarusc(mza,mua)

integer :: j,k,iu,iv,iw,mrl

! Zero out long timestep mass flux components (used for scalar advective 
! transport) so they may be summed over small timesteps

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iu,k)
do j = 1,jtab_u(14)%jend(mrl); iu = jtab_u(14)%iu(j)
!----------------------------------------------------------------------
   call qsub('U',iu)
   do k = 1,mza-1          ! begin at level 2 even if below ground
      umarusc(k,iu) = 0.   ! initialize horiz mass flux for scalar advection
   enddo
enddo
!$omp end parallel do
endif
call rsub('U',14)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = 1,lpw(iw)-2
      rho_old(k,iw) = 0.0
   enddo
   do k = lpw(iw)-1,mza-1
      wmarwsc(k,iw) = 0.        ! initialize vert mass flux for scalar advection
      rho_old(k,iw) = rho(k,iw) ! Save DTL density for use with scalar updates
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',18)

return
end subroutine zero_massflux

!===========================================================================

subroutine timeavg_massflux(wmarwsc, umarusc)

use mem_ijtabs, only: jtab_u, jtab_v, itab_u, itab_v, jtab_w, itab_w, &
                      istp, mrl_endl
use mem_grid,   only: mza, mua, mva, mwa, lpu, lpv, lpw
use misc_coms,  only: io6, nacoust

!$ use omp_lib

implicit none

real(kind=8), intent(inout) :: wmarwsc(mza,mwa)
real(kind=8), intent(inout) :: umarusc(mza,mua)

integer :: j,k,iu,iv,iw,mrl,mrlu,mrlv,mrlw
real :: acoi,acoi2

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iu,mrlu,acoi,k)
do j = 1,jtab_u(20)%jend(mrl); iu = jtab_u(20)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   mrlu = itab_u(iu)%mrlu
   acoi = 1. / float(nacoust(mrlu))
   do k = lpu(iu),mza-1
      umarusc(k,iu) = umarusc(k,iu) * acoi  ! upsum
   enddo

enddo
!$omp end parallel do
endif
call rsub('U',20)

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,mrlw,acoi2,k)
do j = 1,jtab_w(25)%jend(mrl); iw = jtab_w(25)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   mrlw = itab_w(iw)%mrlw
   acoi2 = 1. / float(nacoust(mrlw))
   do k = lpw(iw),mza-2
      wmarwsc(k,iw) = wmarwsc(k,iw) * acoi2
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',25)

return
end subroutine timeavg_massflux

!===========================================================================

subroutine zero_momsc(vmsc,wmsc,rho_old)

use mem_ijtabs, only: jtab_v, jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mva, mwa, lpw
use mem_basic,  only: rho
use misc_coms,  only: io6

!$ use omp_lib

implicit none

real, intent(out) :: vmsc(mza,mva)
real, intent(out) :: wmsc(mza,mwa)
real(kind=8), intent(out) :: rho_old(mza,mwa)

integer :: j,k,iv,iw,mrl

! Zero out long timestep mass flux components (used for scalar advective 
! transport) so they may be summed over small timesteps

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iv,k)
do j = 1,jtab_v(14)%jend(mrl); iv = jtab_v(14)%iv(j)
!----------------------------------------------------------------------
call qsub('V',iv)
   do k = 1,mza-1
      vmsc(k,iv) = 0.
   enddo
enddo
!$omp end parallel do
endif
call rsub('V',14)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = 1,lpw(iw)-2
      rho_old(k,iw) = 0.0
   enddo
   do k = lpw(iw)-1,mza-1
      wmsc(k,iw) = 0.
      rho_old(k,iw) = rho(k,iw) ! Save DTL density for use with scalar updates
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',18)

return
end subroutine zero_momsc

!===========================================================================

subroutine timeavg_momsc(vmsc,wmsc)

use mem_ijtabs, only: jtab_v, itab_v, jtab_w, itab_w, istp, mrl_endl
use mem_grid,   only: mza, mva, mwa, lpv, lpw
use misc_coms,  only: io6, nacoust

!$ use omp_lib

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

integer :: j,k,iv,iw,mrl,mrlv,mrlw
real :: acoi,acoi2

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iv,mrlv,acoi,k)
do j = 1,jtab_v(20)%jend(mrl); iv = jtab_v(20)%iv(j)
!----------------------------------------------------------------------
call qsub('V',iv)
   mrlv = itab_v(iv)%mrlv
   acoi = 1. / float(nacoust(mrlv))
   do k = lpv(iv),mza-1
      vmsc(k,iv) = vmsc(k,iv) * acoi
   enddo
enddo
!$omp end parallel do
endif
call rsub('V',20)

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,mrlw,acoi2,k)
do j = 1,jtab_w(25)%jend(mrl); iw = jtab_w(25)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   mrlw = itab_w(iw)%mrlw
   acoi2 = 1. / float(nacoust(mrlw))
   do k = lpw(iw),mza-2
      wmsc(k,iw) = wmsc(k,iw) * acoi2
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',25)

return
end subroutine timeavg_momsc

!===============================================================================

subroutine diagnose_ucm()

use mem_basic,  only: uc,rho,vmc
use mem_ijtabs, only: mrl_ends,istp,jtab_v,itab_v
use mem_grid,   only: mza,lpv,aru,arv
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer :: iv1,iv2,iv3,iv4,iv5,iv8,iv9,iv12,iv13,iv14,iv15,iv16
integer :: mrl,kb,k,iv,j,iw1,iw2

call psub()
!------------------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iv,iv1,iv2,iv3,iv4,iv5,iv8,iv9, &
!$omp                     iv12,iv13,iv14,iv15,iv16,iw1,iw2,kb,k)
do j = 1,jtab_v(16)%jend(mrl); iv = jtab_v(16)%iv(j)
   iv1  = itab_v(iv)%iv(1) ; iv2  = itab_v(iv)%iv(2) ; iv3  = itab_v(iv)%iv(3)
   iv4  = itab_v(iv)%iv(4) ; iv5  = itab_v(iv)%iv(5) ; iv8  = itab_v(iv)%iv(8)
   iv9  = itab_v(iv)%iv(9) ; iv12 = itab_v(iv)%iv(12); iv13 = itab_v(iv)%iv(13)
   iv14 = itab_v(iv)%iv(14); iv15 = itab_v(iv)%iv(15); iv16 = itab_v(iv)%iv(16)
   iw1  = itab_v(iv)%iw(1) ; iw2  = itab_v(iv)%iw(2)
!------------------------------------------------------------------------------
call qsub('V',iv)

   kb = lpv(iv)

! Loop over T levels

   do k = kb,mza-1

      uc(k,iv) = (itab_v(iv)%fuv(1)  * vmc(k,iv1)  * arv(k,iv1)   &
               +  itab_v(iv)%fuv(2)  * vmc(k,iv2)  * arv(k,iv2)   &
               +  itab_v(iv)%fuv(3)  * vmc(k,iv3)  * arv(k,iv3)   &
               +  itab_v(iv)%fuv(4)  * vmc(k,iv4)  * arv(k,iv4)   &
               +  itab_v(iv)%fuv(5)  * vmc(k,iv5)  * arv(k,iv5)   &
               +  itab_v(iv)%fuv(8)  * vmc(k,iv8)  * arv(k,iv8)   &
               +  itab_v(iv)%fuv(9)  * vmc(k,iv9)  * arv(k,iv9)   &
               +  itab_v(iv)%fuv(12) * vmc(k,iv12) * arv(k,iv12)  &
               +  itab_v(iv)%fuv(13) * vmc(k,iv13) * arv(k,iv13)  &
               +  itab_v(iv)%fuv(14) * vmc(k,iv14) * arv(k,iv14)  &
               +  itab_v(iv)%fuv(15) * vmc(k,iv15) * arv(k,iv15)  &
               +  itab_v(iv)%fuv(16) * vmc(k,iv16) * arv(k,iv16)) &
               / (aru(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2)))

   enddo

enddo
!$omp end parallel do
endif
call rsub('Vb',16)

return
end subroutine diagnose_ucm

!===============================================================================

subroutine diagnose_uc()

use mem_basic,  only: uc,vc
use mem_ijtabs, only: mrl_ends,istp,jtab_v,itab_v
use mem_grid,   only: mza,lpv,dnu,dniv
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer :: iv1,iv2,iv3,iv4,iv5,iv8,iv9,iv12,iv13,iv14,iv15,iv16
integer :: mrl,kb,k,iv,j

call psub()
!------------------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iv,iv1,iv2,iv3,iv4,iv5,iv8,iv9, &
!$omp                     iv12,iv13,iv14,iv15,iv16,kb,k)
do j = 1,jtab_v(16)%jend(mrl); iv = jtab_v(16)%iv(j)
   iv1  = itab_v(iv)%iv(1) ; iv2  = itab_v(iv)%iv(2) ; iv3  = itab_v(iv)%iv(3)
   iv4  = itab_v(iv)%iv(4) ; iv5  = itab_v(iv)%iv(5) ; iv8  = itab_v(iv)%iv(8)
   iv9  = itab_v(iv)%iv(9) ; iv12 = itab_v(iv)%iv(12); iv13 = itab_v(iv)%iv(13)
   iv14 = itab_v(iv)%iv(14); iv15 = itab_v(iv)%iv(15); iv16 = itab_v(iv)%iv(16)
!------------------------------------------------------------------------------
call qsub('V',iv)

   kb = lpv(iv)

! Loop over T levels

   do k = kb,mza-1

      uc(k,iv) = (itab_v(iv)%fuv(1)  * vc(k,iv1) * dnu(iv1)   &
               +  itab_v(iv)%fuv(2)  * vc(k,iv2) * dnu(iv2)   &
               +  itab_v(iv)%fuv(3)  * vc(k,iv3) * dnu(iv3)   &
               +  itab_v(iv)%fuv(4)  * vc(k,iv4) * dnu(iv4)   &
               +  itab_v(iv)%fuv(5)  * vc(k,iv5) * dnu(iv5)   &
               +  itab_v(iv)%fuv(8)  * vc(k,iv8) * dnu(iv8)   &
               +  itab_v(iv)%fuv(9)  * vc(k,iv9) * dnu(iv9)   &
               +  itab_v(iv)%fuv(12) * vc(k,iv12) * dnu(iv12) &
               +  itab_v(iv)%fuv(13) * vc(k,iv13) * dnu(iv13) &
               +  itab_v(iv)%fuv(14) * vc(k,iv14) * dnu(iv14) &
               +  itab_v(iv)%fuv(15) * vc(k,iv15) * dnu(iv15) &
               +  itab_v(iv)%fuv(16) * vc(k,iv16) * dnu(iv16)) * dniv(iv)

   enddo

enddo
!$omp end parallel do
endif
call rsub('Vb',16)

return
end subroutine diagnose_uc

!===============================================================================

subroutine diagnose_vc()

use mem_basic,  only: uc,vc
use mem_ijtabs, only: mrl_ends,istp,jtab_u,itab_u
use mem_grid,   only: mza,lpu,   dnu,dniv
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer :: iu1,iu2,iu3,iu4,iu5,iu8,iu9,iu12,iu13,iu14,iu15,iu16
integer :: mrl,kb,k,iu,j

real :: tuu1, tuu2, tuu3, tuu4

call psub()
!------------------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu,iu1,iu2,iu3,iu4,tuu1,tuu2,tuu3,tuu4,kb,k)
do j = 1,jtab_u(22)%jend(mrl); iu = jtab_u(22)%iu(j)

   iu1 = itab_u(iu)%iu(1)
   iu2 = itab_u(iu)%iu(2)
   iu3 = itab_u(iu)%iu(3)
   iu4 = itab_u(iu)%iu(4)

   tuu1 = itab_u(iu)%tuu(1)
   tuu2 = itab_u(iu)%tuu(2)
   tuu3 = itab_u(iu)%tuu(3)
   tuu4 = itab_u(iu)%tuu(4)
!------------------------------------------------------------------------------
call qsub('U',iu)

   kb = lpu(iu)

! Loop over T levels

   do k = kb,mza-1

      vc(k,iu) = uc(k,iu1) * tuu1  &
               + uc(k,iu2) * tuu2  &
               + uc(k,iu3) * tuu3  &
               + uc(k,iu4) * tuu4

! Alternate form based on divergences will replace tuu(1) with fvu(1), etc?

   enddo

enddo
!$omp end parallel do
endif
call rsub('Ub',22)

return
end subroutine diagnose_vc

!===========================================================================

subroutine tridiffo(m1,ka,kz,cim1,ci,cip1,rhs,soln)

implicit none

integer, intent(in) :: m1,ka,kz
real, intent(in) :: cim1(m1),ci(m1),cip1(m1),rhs(m1)
real, intent(out) :: soln(m1)
real :: scr1(m1)  ! automatic array
real :: scr2(m1)  ! automatic array

integer :: k
real cji

scr1(ka) = cip1(ka) / ci(ka)
scr2(ka) = rhs(ka) / ci(ka)

do k = ka+1,kz
   soln(k) = ci(k) - cim1(k) * scr1(k-1)
   cji = 1. / soln(k)
   scr1(k) = cip1(k) * cji
   scr2(k) = (rhs(k) - cim1(k) * scr2(k-1)) * cji
enddo

soln(kz) = scr2(kz)

do k = kz-1,ka,-1
   soln(k) = scr2(k) - scr1(k) * soln(k+1)
enddo

return
end subroutine tridiffo

End Module massflux
