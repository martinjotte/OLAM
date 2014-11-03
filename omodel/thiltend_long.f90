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
subroutine thiltend_long(mrl)

use mem_ijtabs, only: istp, jtab_w, mrl_begl, jtw_prog
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

integer :: j,iw

! Horizontal loop over W/T points

!----------------------------------------------------------------------
!$omp parallel do private (iw)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   call thiltend_long0(iw)

enddo
!$omp end parallel do

return
end subroutine thiltend_long

!===============================================================================

subroutine thiltend_long0(iw)

use mem_ijtabs,  only: itab_w
use misc_coms,   only: io6, dtlm
use mem_basic,   only: rho, thil
use mem_turb,    only: vkh, hkm, sxfer_tk, sxfer_rk
use mem_tend,    only: thilt
use mem_grid,    only: mza, mwa, arw, dzim, volt, volti, arv, dniv, lpv

implicit none

integer, intent(in) :: iw

integer :: k,ka,ks
integer :: j, iv, iwn, npoly

real :: dtl,dtli
real :: hdniv, hflux

! Automatic arrays:

real, dimension(mza) :: akodz,dtomass,vctr1,vctr5  &
                       ,vctr6,vctr7,vctr8,vctr9,del_thil

dtl  = dtlm(itab_w(iw)%mrlw)
dtli = 1. / dtl

! Sum horizontal diffusive fluxes 

! Number of edges of this IW polygon

npoly = itab_w(iw)%npoly
vctr1(1:mza) = 0.

! Loop over V neighbors of this W cell

do j = 1, npoly
   iv  = itab_w(iw)%iv(j)
   iwn = itab_w(iw)%iw(j)
   ka  = lpv(iv)

   hdniv = .5 * dniv(iv)  ! use this 1/dx form now - it seems better than A/V

   do k = ka, mza

! Horizontal turbulent flux across ARV

      hflux = arv(k,iv) * (hkm(k,iwn) + hkm(k,iw))  &
            * hdniv * (thil(k,iwn) - thil(k,iw))

      vctr1(k) = vctr1(k) + hflux
   enddo
enddo

! Vertical loop over T levels

do k = ka,mza

! Update thil tendency from horizontal turbulent fluxes

   thilt(k,iw) = thilt(k,iw) + volti(k,iw) * vctr1(k)

enddo

return
end subroutine thiltend_long0
