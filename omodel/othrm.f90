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
subroutine thermo()

use mem_ijtabs, only: jtab_w, istp, mrl_endl, jtw_prog
use micro_coms, only: miclevel
use misc_coms,  only: io6

implicit none

integer iw,j,mrl

! Horizontal loop over W/T points

!-------------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then

!$omp parallel do private (iw)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!-------------------------------------------------------------------------

   if (miclevel <= 1) then
      call drythrm(iw)
   elseif (miclevel == 2) then
      call satadjst(iw)
   elseif (miclevel == 3) then
      call wetthrm3(iw)
   else
      stop 'Thermo option not supported...MICLEVEL out of bounds'
   endif

enddo
!$omp end parallel do
endif

end subroutine thermo

!===============================================================================

subroutine drythrm(iw)

! This routine calculates theta and rv for the case where no condensate is
! allowed.

use mem_basic,  only: theta, thil, tair, sh_v, sh_w, press
use micro_coms, only: miclevel
use mem_grid,   only: lpw, mza
use misc_coms,  only: io6
use consts_coms,only: p00i, rocp

implicit none

integer, intent(in) :: iw

integer k

do k = lpw(iw),mza
   theta(k,iw) = thil(k,iw)
   tair (k,iw) = thil(k,iw) * (real(press(k,iw)) * p00i) ** rocp
enddo

if (miclevel == 1) then
   do k = lpw(iw),mza
      sh_v(k,iw) = sh_w(k,iw)
   enddo
endif

end subroutine drythrm

!===============================================================================

subroutine satadjst(iw)

! This routine diagnoses theta, sh_v, and sh_c using a saturation adjustment
! for the case when sh_c is the only allowed condensate

use mem_basic,   only: thil, theta, tair, rho, sh_w, sh_v, press
use mem_micro,   only: sh_c
use consts_coms, only: p00i, cp, alvl, rocp
use mem_grid,    only: lpw, mza
use misc_coms,   only: io6
use therm_lib,   only: rhovsl
implicit none

integer, intent(in) :: iw

integer :: iterate,k
real :: temp,rlvs,rt,rc,t_il,rvls,exner,rhovs

do k = lpw(iw),mza
   exner = (real(press(k,iw)) * p00i) ** rocp  ! Defined WITHOUT CP factor
   t_il = thil(k,iw) * exner
   temp = t_il

   do iterate = 1,20
      rhovs = rhovsl(temp-273.15)
      sh_c(k,iw) = max(0.,sh_w(k,iw)-rhovs/real(rho(k,iw)))
      sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)
      temp = 0.7 * temp  &
           + 0.3 * t_il * (1. + alvl * sh_c(k,iw) / (cp * max(temp,253.)))
   enddo

   theta(k,iw) = temp / exner
   tair (k,iw) = temp
enddo

end subroutine satadjst

!===============================================================================

subroutine wetthrm3(iw)

! This routine calculates theta and sh_v for "miclevel 3 microphysics"
! given prognosed theta_il, cloud, drizzle, rain, pristine ice, snow, 
! aggregates, graupel, hail, q6, and q7.

use mem_basic,   only: press, theta, thil, tair, sh_v, sh_w
use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, q6, q7
use micro_coms,  only: jnmb, rxmin
use consts_coms, only: p00i, rocp, alvl, alvi, cpi4, cp253i
use mem_grid,    only: mza, lpw
use misc_coms,   only: io6
use therm_lib,   only: qtc

implicit none

integer, intent(in) :: iw

integer :: k
real :: tcoal,fracliq,tairstr

real :: exner   (mza)  ! automatic array
real :: airtempk(mza)  ! automatic array
real :: til     (mza)  ! automatic array
real :: totliq  (mza)  ! automatic array
real :: totice  (mza)  ! automatic array
real :: qhydm   (mza)  ! automatic array

do k = lpw(iw),mza
   exner(k) = (real(press(k,iw)) * p00i) ** rocp  ! exner WITHOUT CP factor
   til(k) = thil(k,iw) * exner(k)          ! ice-liquid temperature T_il
   totliq(k) = 0.                          ! total liquid spec. density
   totice(k) = 0.                          ! total ice mixing ratio
enddo

if (jnmb(1) >= 1) then
   do k = lpw(iw),mza
      totliq(k) = totliq(k) + sh_c(k,iw)
   enddo
endif

if (jnmb(2) >= 1) then
   do k = lpw(iw),mza
      totliq(k) = totliq(k) + sh_r(k,iw)
   enddo
endif

if (jnmb(3) >= 1) then
   do k = lpw(iw),mza
      totice(k) = totice(k) + sh_p(k,iw)
   enddo
endif

if (jnmb(4) >= 1) then
   do k = lpw(iw),mza
      totice(k) = totice(k) + sh_s(k,iw)
   enddo
endif

if (jnmb(5) >= 1) then
   do k = lpw(iw),mza
      totice(k) = totice(k) + sh_a(k,iw)
   enddo
endif

if (jnmb(6) >= 1) then
   do k = lpw(iw),mza
      if (sh_g(k,iw) > rxmin(6)) then
         call qtc(q6(k,iw)/sh_g(k,iw),tcoal,fracliq)
         totliq(k) = totliq(k) + sh_g(k,iw) * fracliq
         totice(k) = totice(k) + sh_g(k,iw) * (1. - fracliq)
      endif
   enddo
endif

if (jnmb(7) >= 1) then
   do k = lpw(iw),mza
      if (sh_h(k,iw) > rxmin(7)) then
         call qtc(q7(k,iw)/sh_g(k,iw),tcoal,fracliq)
         totliq(k) = totliq(k) + sh_h(k,iw) * fracliq
         totice(k) = totice(k) + sh_h(k,iw) * (1. - fracliq)
      endif
   enddo
endif

if (jnmb(8) >= 1) then
   do k = lpw(iw),mza
      totliq(k) = totliq(k) + sh_d(k,iw)
   enddo
endif

do k = lpw(iw),mza
   qhydm(k) = alvl * totliq(k) + alvi * totice(k)
   sh_v(k,iw) = sh_w(k,iw) - totliq(k) - totice(k)
enddo

do k = lpw(iw),mza
   if (tair(k,iw) > 253.) then
      tairstr = .5 * (til(k)  &
         + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
   else
      tairstr = til(k) * (1. + qhydm(k) * cp253i)
   endif
   theta(k,iw) = tairstr / exner(k)
   tair (k,iw) = tairstr
enddo

end subroutine wetthrm3


