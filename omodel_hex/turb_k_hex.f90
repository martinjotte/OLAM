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
subroutine turb_k_hex(mrl,vxe,vye,vze)

! This is Smagorinsky-Lilly-Hill turbulence parameterization

use mem_turb,   only: hkm, vkm, vkh
use mem_ijtabs, only: istp, jtab_w, itab_w, mrl_endl
use mem_grid,   only: mza, mwa, lpw, dzim
use misc_coms,  only: io6, idiffk, csx, csz, zkhkm, akmin, meshtype

use mem_basic,   only: rho, theta, sh_v
use consts_coms, only: vonk, alvl, grav2, cp, rvap
use mem_grid,    only: mza, lpw, arw0, dzm, dzim, zm, lpu, lpv,  &
                       unx, uny, unz, vnx, vny, vnz
use micro_coms,  only: level
!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in) :: vxe(mza,mwa)
real, intent(in) :: vye(mza,mwa)
real, intent(in) :: vze(mza,mwa)

integer :: j,iw,k,mrlw,ka

real :: richnum,ambda,ambda2,hill_term,richnum_term
real :: scalen_asympt,scalen_vert,scalen_horiz,bkmin

real :: thetav(mza)
real :: strain2(mza)
real :: bvfreq2(mza)
real :: vkz2(mza)
real :: dzim2 (mza)

real, parameter :: rchmax = 3.  ! Test with new asympt vert scale length
real, parameter :: rmin = -100.
real, parameter :: rmax = 1. / 3.  !  1. / zkhkm(mrl)

do k = 2,mza-2
   dzim2(k) = dzim(k) * dzim(k)
enddo

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iw,mrlw,ka,k,scalen_horiz,scalen_asympt,thetav, &
!$omp                      strain2,bvfreq2,richnum,richnum_term, &
!$omp                      scalen_vert,ambda,ambda2,vkz2,hill_term,bkmin)
do j = 1,jtab_w(34)%jend(mrl); iw = jtab_w(34)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   mrlw = itab_w(iw)%mrlw
   ka = lpw(iw)

! For IDIFFK not equal to 2, set diffusion coefficients to zero
   
   if (idiffk(mrlw) /= 2) then
      do k = ka-1,mza-1
         hkm(k,iw) = 0.
         vkm(k,iw) = 0.
         vkh(k,iw) = 0.
      enddo
      cycle
   endif
   
   scalen_horiz = csx(mrlw) * sqrt(arw0(iw))  ! change this later?
   scalen_asympt = csx(mrlw) * 300.

! Loop over T levels

   do k = ka,mza-1
      thetav(k) = theta(k,iw) * (1. + .61 * sh_v(k,iw))
   enddo

! Loop over W levels

   do k = ka,mza-2

! Vertical strain rate squared (based on dV/dz only)

      strain2(k) = dzim2(k) * ((vxe(k+1,iw) - vxe(k,iw))**2  &
                             + (vye(k+1,iw) - vye(k,iw))**2  &
                             + (vze(k+1,iw) - vze(k,iw))**2)

      bvfreq2(k) = grav2 * dzim(k)  &
         * (thetav(k+1) - thetav(k)) / (thetav(k+1) + thetav(k)) 

!  Compute Richardson number and Lilly Richardson-number term

      richnum = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-15)))
      richnum_term = min(rchmax,sqrt(max(0.,(1.-zkhkm(mrlw)*richnum))))

! Compute vertical and net scale lengths: scalen_vert & ambda

      scalen_vert = csz(mrlw) * dzm(k)
      ambda = max(scalen_vert,min(scalen_asympt,scalen_horiz))
      ambda2 = ambda ** 2

      vkz2(k) = (vonk * (zm(k) - zm(ka-1))) ** 2

      if (bvfreq2(k) < -1.e-12) then
         hill_term = sqrt(-bvfreq2(k))
      else
         hill_term = 0.
      endif

      vkm(k,iw) = .5 * (rho(k,iw) + rho(k+1,iw))        & ! density factor
                * vkz2(k) * ambda2 / (vkz2(k) + ambda2) & ! lengthscale^2 factor
                * (sqrt(strain2(k)) + hill_term)        & ! strain rate + Hill term
                * richnum_term                            ! Lilly Richnum term

      vkh(k,iw) = vkm(k,iw) * zkhkm(mrlw)  ! Need to change this factor later

   enddo
      
! Zero values for top and bottom boundaries   
   
   vkm(ka-1,iw) = 0.
   vkh(ka-1,iw) = 0.
   vkm(mza-1,iw) = 0.
   vkh(mza-1,iw) = 0.

! Horizontal diffusion coefficient (current version)

   bkmin = akmin(1) * .075 * arw0(iw) ** .666667

! akmin hardwired for "grid 1" value for now

   do k = ka,mza-1

!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      hkm(k,iw) = max(vkm(k-1,iw),vkm(k,iw),bkmin * real(rho(k,iw)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo

enddo
!$omp end parallel do
call rsub('W',34)

return
end subroutine turb_k_hex

