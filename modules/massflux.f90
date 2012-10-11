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

use mem_ijtabs, only: jtab_u, jtab_v, jtab_w, istp, mrl_begl, &
                      jtu_wstn, jtw_wstn
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
do j = 1,jtab_u(jtu_wstn)%jend(mrl); iu = jtab_u(jtu_wstn)%iu(j)
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
do j = 1,jtab_w(jtw_wstn)%jend(mrl); iw = jtab_w(jtw_wstn)%iw(j)
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
                      istp, mrl_endl, jtu_prog, jtw_prog
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
do j = 1,jtab_u(jtu_prog)%jend(mrl); iu = jtab_u(jtu_prog)%iu(j)
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
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
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

use mem_ijtabs, only: jtab_v, jtab_w, istp, mrl_begl, jtv_wstn, jtw_wstn
use mem_grid,   only: mza, mva, mwa, lpw
use mem_basic,  only: rho
use misc_coms,  only: io6

!$ use omp_lib

implicit none

real, intent(out) :: vmsc(mza,mva)
real, intent(out) :: wmsc(mza,mwa)
real(kind=8), intent(out) :: rho_old(mza,mwa)

integer :: j,k,iv,iw,mrl

mrl = mrl_begl(istp)

! Zero out long timestep mass flux components (used for scalar advective 
! transport) so they may be summed over small timesteps

if (mrl > 0) then

   call psub()
!----------------------------------------------------------------------
   !$omp parallel
   do iv = 2, mva
!----------------------------------------------------------------------
   call qsub('V',iv)

      vmsc(:,iv) = 0.0

   enddo
   !$omp end parallel do
   call rsub('V',14)

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(k)
   do iw = 2, mwa
!----------------------------------------------------------------------
   call qsub('W',iw)

      wmsc(:,iw) = 0.0

      ! Save DT_LONG density for use with scalar updates

      do k = 1, lpw(iw)-2
         rho_old(k,iw) = 0.0
      enddo
      do k = lpw(iw)-1, mza-1
         rho_old(k,iw) = rho(k,iw)
      enddo

   enddo
   !$omp end parallel do
   call rsub('W',18)

endif

return
end subroutine zero_momsc

!===========================================================================

subroutine timeavg_momsc(vmsc,wmsc)

use mem_ijtabs, only: jtab_v, itab_v, jtab_w, itab_w, istp, mrl_endl, &
                      jtv_prog, jtw_prog
use mem_grid,   only: mza, mva, mwa, lpv, lpw
use misc_coms,  only: io6, nacoust

!$ use omp_lib

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

integer :: j,k,iv,iw,mrl,mrlv,mrlw
real :: acoi

mrl = mrl_endl(istp)
if (mrl > 0) then

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(mrlv,acoi,k)
   do iv = 2, mva
!----------------------------------------------------------------------
   call qsub('V',iv)

      mrlv = itab_v(iv)%mrlv
      acoi = 1.0 / real(nacoust(mrlv))

      do k = lpv(iv), mza-1
         vmsc(k,iv) = vmsc(k,iv) * acoi
      enddo

   enddo
   !$omp end parallel do
   call rsub('V',20)

   call psub()
!----------------------------------------------------------------------
   !$omp parallel do private(mrlw,acoi,k)
   do iw = 2, mwa
!----------------------------------------------------------------------
   call qsub('W',iw)

      mrlw = itab_w(iw)%mrlw
      acoi = 1.0 / real(nacoust(mrlw))

      wmsc(lpw(iw)-1,iw) = 0.0

      do k = lpw(iw), mza-2
         wmsc(k,iw) = wmsc(k,iw) * acoi
      enddo

      wmsc(mza-1,iw) = 0.0

   enddo
   !$omp end parallel do
   call rsub('W',25)

endif

return
end subroutine timeavg_momsc

!===============================================================================
!
!subroutine diagnose_ucm()
!
!use mem_basic,  only: uc, rho, vmc
!use mem_ijtabs, only: mrl_ends, istp, jtab_v, itab_v, jtv_prog
!use mem_grid,   only: mza, lpv, aru, arv
!use misc_coms,  only: io6
!
!!$ use omp_lib
!
!implicit none
!
!! CURRENTLY NOT NEEDED WITH THE PEROT METHOD
!! TODO: If needed, compute uc from earth-cartesian velocities for hexagons
!
!return
!end subroutine diagnose_ucm
!
!===============================================================================
!
!subroutine diagnose_uc()
!
!use mem_basic,  only: uc, vc
!use mem_ijtabs, only: mrl_ends, istp, jtab_v, itab_v, jtv_prog
!use mem_grid,   only: mza, lpv, dnu, dniv
!use misc_coms,  only: io6
!
!!$ use omp_lib
!
!implicit none
!
!! CURRENTLY NOT NEEDED WITH THE PEROT METHOD
!! TODO: If needed, compute uc from earth-cartesian velocities for hexagons
!
!return
!end subroutine diagnose_uc
!
!===============================================================================

subroutine diagnose_vc()

use mem_basic,  only: uc, vc
use mem_ijtabs, only: mrl_ends, istp, jtab_u, itab_u, jtu_wadj
use mem_grid,   only: mza, lpu, dnu, dniv
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
do j = 1,jtab_u(jtu_wadj)%jend(mrl); iu = jtab_u(jtu_wadj)%iu(j)

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

End Module massflux
