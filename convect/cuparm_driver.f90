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
subroutine cuparm_driver(rhot)

use mem_grid,    only: mwa, mza, lpw
use grell_coms,  only: alloc_grell, iact_gr
use misc_coms,   only: io6, time_istp8, nqparm, nqparm_sh, confrq, dtlong,  &
                       initial, itime1
use mem_ijtabs,  only: itab_w, jtab_w, mrl_begl, istp, mrls
use mem_cuparm,  only: thsrc, rtsrc, thsrcsh, rtsrcsh, aconpr, conprr
use mem_tend,    only: thilt, sh_wt

use mem_basic,   only: wc, theta, press, rho, sh_v
use consts_coms, only: tkmin, r8
use mem_micro,   only: sh_c

!$ use omp_lib

implicit none

real, intent(inout) :: rhot(mza,mwa)
real, parameter :: cptime = 7200.0

integer :: j
integer :: iw
integer :: k
integer :: mrl, mrlw

! Memory allocation for Grell cumulus transport

if ( any(nqparm(1:mrls) == 2) ) then
   if (.not. allocated(iact_gr)) then
      call alloc_grell(mwa)
   endif
endif

! If model run has been initialized from observational data, avoid cumulus
! parameterization during an initial period to allow gravity waves to settle.

if ((initial == 2) .and. (time_istp8 < real(cptime,r8))) return

! Check whether it is time to update cumulus parameterization tendencies

if ((istp == 1) .and. (mod(time_istp8+0.001_r8,real(confrq,r8)) < dtlong)) then

! Print message that cumulus parameterization is being computed

   write(io6, '(A,f0.2,A)') &
        ' Convective tendencies updated at ', time_istp8/3600., &
        ' hrs into simulation'

! Initialize indices for search for maximum heating rate
! (commented out because incorrect in parallel operation)

!    dthmax = 0.
!    iwqmax = 0
!    kqmax = 0

! Loop over all IW grid cells where cumulus parameterization may be done

   call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,mrlw) 
   do j = 1,jtab_w(15)%jend(1); iw = jtab_w(15)%iw(j) ! jend(1) for mrl = 1
!----------------------------------------------------------------------
   call qsub('W',iw)

! MRL for current IW column

      mrlw = itab_w(iw)%mrlw

! Select cumulus convection scheme based on MRL of current IW column

      if (nqparm(mrlw) == 1) then
   
! Kuo deep convection

         call cuparm_kuo(iw)

      elseif (nqparm(mrlw) == 2) then
   
! Grell deep convection

         call cuparth(iw,dtlong)

      endif
   
      if (nqparm_sh(mrlw) == 1) then
   
! Grell shallow convection

         call cuparth_shal(iw,dtlong)
      
      endif

! Update indices for maximum heating rate
! (commented out because incorrect in parallel operation)

!   do k = lpw(iw),mza-1
!      ftcon = thsrc(k,iw) / rho(k,iw)
!      if (ftcon > dthmax) then
!         dthmax = ftcon
!         iwqmax = iw
!         kqmax = k
!      endif
!   enddo

   enddo
!$omp end parallel do
   call rsub('Wa',15)

! Print maximum heating rate
! (commented out because incorrect in parallel operation)

!   write(io6, '(A,I0,A,I0,A,F0.3,A)') " MAX CONVECTIVE HEATING RATE AT IW=",  &
!        iwqmax, " K=", kqmax, " IS ", dthmax*86400., " K/DAY"


endif

! Add current value of convective tendencies to thilt and sh_wt arrays
! every long timestep (whether cumulus parameterization is updated 
! this timestep or not)

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k) 
do j = 1,jtab_w(15)%jend(mrl); iw = jtab_w(15)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   aconpr(iw) = aconpr(iw) + conprr(iw) * dtlong

   do k = lpw(iw),mza-1
      thilt(k,iw) = thilt(k,iw)  &
                  + (thsrc(k,iw) + thsrcsh(k,iw)) * rho(k,iw)

      sh_wt(k,iw) = sh_wt(k,iw)  &
                  + (rtsrc(k,iw) + rtsrcsh(k,iw)) * rho(k,iw)

      rhot (k,iw) = rhot (k,iw)  &
                  + (rtsrc(k,iw) + rtsrcsh(k,iw)) * rho(k,iw)
   enddo

enddo
!$omp end parallel do
endif
call rsub('Wb',15)

return
end subroutine cuparm_driver
