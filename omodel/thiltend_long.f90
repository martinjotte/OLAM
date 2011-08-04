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
subroutine thiltend_long(mrl,rhot)

use mem_ijtabs, only: istp, jtab_w, mrl_begl
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(inout) :: rhot(mza,mwa)

integer :: j,iw

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iw)
do j = 1,jtab_w(17)%jend(mrl); iw = jtab_w(17)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call thiltend_long0(iw,rhot)

enddo
!$omp end parallel do
call rsub('W',17)

return
end subroutine thiltend_long

!===============================================================================

subroutine thiltend_long0(iw,rhot)

use mem_ijtabs,  only: itab_w
use misc_coms,   only: io6, dtlm, meshtype
use mem_basic,   only: rho, thil
use mem_turb,    only: vkh, hkm, sxfer_tk, sxfer_rk
use mem_tend,    only: thilt
use mem_grid,    only: mza, mwa, lpw, lsw, arw, dzim, volt, volti, &
                       aru, arv, dniu, dniv
use massflux,    only: tridiffo

implicit none

integer, intent(in) :: iw

real, intent(inout) :: rhot(mza,mwa)

integer :: k,ka,ks
integer :: j, iu, iv, iwn, npoly

real :: dtl,dtli
real :: hdniu, hdniv, hflux

! Automatic arrays:

real, dimension(mza) :: akodz,dtomass,vctr1,vctr5  &
                       ,vctr6,vctr7,vctr8,vctr9,del_thil

dtl = dtlm(itab_w(iw)%mrlw)
dtli = 1. / dtl
ka = lpw(iw)

! Vertical loop over T levels 

do k = ka,mza-1
   dtomass(k) = dtl / (rho(k,iw) * volt(k,iw)) 
enddo

! Vertical loop over W levels: fill tri-diagonal matrix coefficients and r.h.s.

do k = ka,mza-2
   akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
   vctr5(k) = - akodz(k) * dtomass(k)
   vctr7(k) = - akodz(k) * dtomass(k+1)
   vctr6(k) = 1. - vctr5(k) - vctr7(k)
   vctr8(k) = akodz(k) * (thil(k,iw) - thil(k+1,iw))
enddo

! Vertical loop over T levels that are adjacent to surface

do ks = 1,lsw(iw)
   k = ka + ks - 1

! Apply surface heat xfer [kg_a K] directly to thilt [kg_a K / s]

   thilt(k,iw) = thilt(k,iw) + dtli * volti(k,iw) * sxfer_tk(ks,iw)

! Apply surface vapor xfer [kg_vap] directly to rhot [kg_air / s]

   rhot(k,iw) = rhot(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

! Change in thil from surface heat xfer

   del_thil(k) = sxfer_tk(ks,iw) / (rho(k,iw) * volt(k,iw))

! Zero out sxfer_tk(ks,iw) now that it has been transferred to the atm

   sxfer_tk(ks,iw) = 0.  
enddo

! Lowest T level that is not adjacent to surface

del_thil(ka+lsw(iw)) = 0.

! Vertical loop over W levels that are adjacent to surface

do ks = 1,lsw(iw)
   k = ka + ks - 1

! Change in vctr8 from surface heat xfer

   vctr8(k) = vctr8(k) + akodz(k) * (del_thil(k) - del_thil(k+1))
enddo

! Solve tri-diagonal matrix

if (ka < mza-2) then
   call tridiffo(mza,ka,mza-2,vctr5,vctr6,vctr7,vctr8,vctr9)
endif

! Set bottom and top internal turbulent fluxes to zero

vctr9(ka-1) = 0.
vctr9(mza-1) = 0.

! Sum horizontal diffusive fluxes 

! Number of edges of this IW polygon

npoly = itab_w(iw)%npoly
vctr1(1:mza) = 0.

if (meshtype == 1) then

! If mesh is triangular, loop over U neighbors of this W cell

   do j = 1,npoly
      iu  = itab_w(iw)%iu(j)
      iwn = itab_w(iw)%iw(j)

      hdniu = .5 * dniu(iu)  ! use this 1/dx form now - it seems better than A/V

      do k = ka,mza-1

! Horizontal turbulent flux across ARU

         hflux = aru(k,iu) * (hkm(k,iwn) + hkm(k,iw))  &
               * hdniu * (thil(k,iwn) - thil(k,iw))

         vctr1(k) = vctr1(k) + hflux
      enddo
   enddo

else

! If mesh is hexagonal, loop over V neighbors of this W cell

   do j = 1,npoly
      iv  = itab_w(iw)%iv(j)
      iwn = itab_w(iw)%iw(j)

      hdniv = .5 * dniv(iv)  ! use this 1/dx form now - it seems better than A/V

      do k = ka,mza-1

! Horizontal turbulent flux across ARV

         hflux = arv(k,iv) * (hkm(k,iwn) + hkm(k,iw))  &
               * hdniv * (thil(k,iwn) - thil(k,iw))

         vctr1(k) = vctr1(k) + hflux
      enddo
   enddo

endif

! Vertical loop over T levels

do k = ka,mza-1

! Update thil tendency from vertical and horizontal turbulent fluxes

   thilt(k,iw) = thilt(k,iw) &
               + volti(k,iw) * (vctr9(k-1) - vctr9(k) + vctr1(k))

enddo

return
end subroutine thiltend_long0
