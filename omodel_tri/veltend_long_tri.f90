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
subroutine veltend_long()

use mem_ijtabs,   only: jtab_u, jtab_w, istp, mrl_begl, itab_u, &
                        jtu_prog, jtw_prog
use mem_grid,     only: mza, lpu, lpw, unx, uny, unz, vnx, vny, wnx, wny, wnz
use misc_coms,    only: io6, iparallel
use mem_tend,     only: vmxet, vmyet, vmzet, wmt, umt
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v
use obnd,         only: lbcopy_w

!$ use omp_lib

implicit none

integer :: j
integer :: k
integer :: iu
integer :: iw, iw1, iw2
integer :: mrl

mrl = mrl_begl(istp)
if (mrl == 0) return

 call lbcopy_w(mrl, a1=vmxet, a2=vmyet, a3=vmzet)

! MPI SEND of VMXET, VMYET, VMZET

if (iparallel == 1) call mpi_send_w('V',vmxet=vmxet,vmyet=vmyet,vmzet=vmzet)

! Horizontal loop over U points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iu)
do j = 1,jtab_u(jtu_prog)%jend(mrl); iu = jtab_u(jtu_prog)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call veltend_long_u(iu)

enddo
!$omp end parallel do
call rsub('U',12)

! Horizontal loop over W points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iw,k)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call veltend_long_w(iw)

   ! Update WM tendency from long-timestep earth-cartesian tendencies
   ! Can be computed before communication of tendencies is completed

   do k = lpw(iw), mza-2

      wmt(k,iw) = wmt(k,iw)                                        &
                + 0.5 * ( wnx(iw) * (vmxet(k,iw) + vmxet(k+1,iw) ) &
                        + wny(iw) * (vmyet(k,iw) + vmyet(k+1,iw) ) &
                        + wnz(iw) * (vmzet(k,iw) + vmzet(k+1,iw) ) )
   enddo

enddo
!$omp end parallel do
call rsub('W',16)

! MPI RECV of VMXET, VMYET, VMZET

if (iparallel == 1) call mpi_recv_w('V',vmxet=vmxet,vmyet=vmyet,vmzet=vmzet)

! Update UM tendency from long-timestep earth-cartesian tendencies

! Horizontal loop over U points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iu,iw1,iw2,k)
do j = 1,jtab_u(jtu_prog)%jend(mrl); iu = jtab_u(jtu_prog)%iu(j)
!----------------------------------------------------------------------
   call qsub('U',iu)

   iw1 = itab_u(iu)%iw(1)
   iw2 = itab_u(iu)%iw(2)

   do k = lpu(iu), mza-1
         
! Update UM tendency from turbulent fluxes

      umt(k,iu) = umt(k,iu)                                        &
                + 0.5 * ( unx(iu) * (vmxet(k,iw1) + vmxet(k,iw2) ) &
                        + uny(iu) * (vmyet(k,iw1) + vmyet(k,iw2) ) &
                        + unz(iu) * (vmzet(k,iw1) + vmzet(k,iw2) ) )
   enddo

enddo
!$omp end parallel do
call rsub('U',12)

return
end subroutine veltend_long

!===============================================================================

subroutine veltend_long_u(iu)

use mem_tend,    only: umt
use mem_ijtabs,  only: itab_u
use mem_basic,   only: rho, uc, wc
use misc_coms,   only: io6, dtlm, icorflg
use mem_turb,    only: hkm
use consts_coms, only: omega2
use mem_grid,    only: mza, aru, arw, volt, dniu, dzim, unx, uny, lpu, arw0
use massflux,    only: tridiffo

implicit none

integer, intent(in) :: iu

integer :: k,ka
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6

real :: fuu5,fuu6,fuu7,fuu8,fuu9,fuu10,fuu11,fuu12
real :: fuw3,fuw4,fuw5,fuw6
real :: vxu1_u,vyu1_u
real :: vxu2_u,vyu2_u
real :: vxu3_u,vyu3_u
real :: vxu4_u,vyu4_u
real :: vxw1_u,vyw1_u
real :: vxw2_u,vyw2_u

real :: dtl
real :: vproj1,vproj2,vproj3,vproj4
real :: vx1,vx2,vy1,vy2
real :: qdniu1,qdniu2,qdniu3,qdniu4
real :: umass
real :: wt1,wt2                     ! velocity weights for Coriolis force

! Automatic arrays:

real, dimension(mza) :: umt_cor

iu1  = itab_u(iu)%iu(1);  iu2  = itab_u(iu)%iu(2);  iu3  = itab_u(iu)%iu(3)
iu4  = itab_u(iu)%iu(4);  iu5  = itab_u(iu)%iu(5);  iu6  = itab_u(iu)%iu(6)
iu7  = itab_u(iu)%iu(7);  iu8  = itab_u(iu)%iu(8);  iu9  = itab_u(iu)%iu(9)
iu10 = itab_u(iu)%iu(10); iu11 = itab_u(iu)%iu(11); iu12 = itab_u(iu)%iu(12)

iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2); iw3 = itab_u(iu)%iw(3)
iw4 = itab_u(iu)%iw(4); iw5 = itab_u(iu)%iw(5); iw6 = itab_u(iu)%iw(6)

fuu5  = itab_u(iu)%fuu(5);  fuu6  = itab_u(iu)%fuu(6);  fuw3 = itab_u(iu)%fuw(3)
fuu7  = itab_u(iu)%fuu(7);  fuu8  = itab_u(iu)%fuu(8);  fuw4 = itab_u(iu)%fuw(4)
fuu9  = itab_u(iu)%fuu(9);  fuu10 = itab_u(iu)%fuu(10); fuw5 = itab_u(iu)%fuw(5)
fuu11 = itab_u(iu)%fuu(11); fuu12 = itab_u(iu)%fuu(12); fuw6 = itab_u(iu)%fuw(6)
   
vxu1_u = itab_u(iu)%vxu_u(1); vyu1_u = itab_u(iu)%vyu_u(1)
vxu2_u = itab_u(iu)%vxu_u(2); vyu2_u = itab_u(iu)%vyu_u(2)
vxu3_u = itab_u(iu)%vxu_u(3); vyu3_u = itab_u(iu)%vyu_u(3)
vxu4_u = itab_u(iu)%vxu_u(4); vyu4_u = itab_u(iu)%vyu_u(4)
vxw1_u = itab_u(iu)%vxw_u(1); vyw1_u = itab_u(iu)%vyw_u(1)
vxw2_u = itab_u(iu)%vxw_u(2); vyw2_u = itab_u(iu)%vyw_u(2)

wt1 = arw0(iw2) / (arw0(iw1) + arw0(iw2))
wt2 = 1. - wt1

dtl = dtlm(itab_u(iu)%mrlu)

ka = lpu(iu)

! Vertical loop over T levels

do k = ka,mza-1

! Mass in U control volume; and its inverse times dtl

   umass = rho(k,iw1) * volt(k,iw1) + rho(k,iw2) * volt(k,iw2)

! Coriolis force

   if (icorflg == 1) then

      vx1 = vxu1_u * uc(k,iu1)  &
          + vxu2_u * uc(k,iu2)  &
          + vxw1_u * (wc(k-1,iw1) + wc(k,iw1))

      vy1 = vyu1_u * uc(k,iu1)  &
          + vyu2_u * uc(k,iu2)  &
          + vyw1_u * (wc(k-1,iw1) + wc(k,iw1))

      vx2 = vxu3_u * uc(k,iu3)  &
          + vxu4_u * uc(k,iu4)  &
          + vxw2_u * (wc(k-1,iw2) + wc(k,iw2))

      vy2 = vyu3_u * uc(k,iu3)  &
          + vyu4_u * uc(k,iu4)  &
          + vyw2_u * (wc(k-1,iw2) + wc(k,iw2))

      umt_cor(k) = omega2 * umass * ((wt1 * vy1 + wt2 * vy2) * unx(iu)  &
                                   - (wt1 * vx1 + wt2 * vx2) * uny(iu))

   else
      umt_cor(k) = 0.
   endif

enddo

! Coefficients for horizontal turbulent fluxes: quarter dniu since summing over
! two K values and since horizontal gradient of u is effectively over 2 deltax

qdniu1 = .25 * dniu(iu1)
qdniu2 = .25 * dniu(iu2)
qdniu3 = .25 * dniu(iu3)
qdniu4 = .25 * dniu(iu4)

! Vertical loop over T levels

do k = ka,mza-1

! Projected velocities

   vproj1 = fuu5  * uc(k,iu5)               &
          + fuu6  * uc(k,iu6)               &
          + fuw3  * (wc(k,iw3) + wc(k-1,iw3))

   vproj2 = fuu7  * uc(k,iu7)               &
          + fuu8  * uc(k,iu8)               &
          + fuw4  * (wc(k,iw4) + wc(k-1,iw4))

   vproj3 = fuu9  * uc(k,iu9)               &
          + fuu10 * uc(k,iu10)              &
          + fuw5  * (wc(k,iw5) + wc(k-1,iw5))

   vproj4 = fuu11 * uc(k,iu11)              &
          + fuu12 * uc(k,iu12)              &
          + fuw6  * (wc(k,iw6) + wc(k-1,iw6))

! Update UM tendency from turbulent fluxes and Coriolis force

   umt(k,iu) = umt(k,iu)                        &

      + aru(k,iu1) * (hkm(k,iw1) + hkm(k,iw3))  &
          * qdniu1 * (vproj1 - uc(k,iu))        &

      + aru(k,iu2) * (hkm(k,iw1) + hkm(k,iw4))  &
          * qdniu2 * (vproj2 - uc(k,iu))        &

      + aru(k,iu3) * (hkm(k,iw2) + hkm(k,iw5))  &
          * qdniu3 * (vproj3 - uc(k,iu))        &
             
      + aru(k,iu4) * (hkm(k,iw2) + hkm(k,iw6))  &
          * qdniu4 * (vproj4 - uc(k,iu))        &
             
      + umt_cor(k)

enddo

return
end subroutine veltend_long_u

!===============================================================================

subroutine veltend_long_w(iw)

use mem_tend,    only: wmt
use mem_ijtabs,  only: itab_w
use mem_basic,   only: uc,wc,rho
use misc_coms,   only: io6,icorflg, dtlm
use mem_turb,    only: hkm
use consts_coms, only: omega2
use mem_grid,    only: mza, lpw, aru, arw, dniu, dzit, volt, volwi, wnx, wny
use massflux,    only: tridiffo

implicit none

integer, intent(in) :: iw

integer :: k,ka
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
integer :: iw1,iw2,iw3

real :: fwu4,fwu5,fwu6,fwu7,fwu8,fwu9
real :: fww1,fww2,fww3
real :: vxu1_w,vyu1_w
real :: vxu2_w,vyu2_w
real :: vxu3_w,vyu3_w

real :: term,dtl2,sflux
real :: vproj1,vproj2,vproj3
real :: vx_w,vy_w
real :: arukodx1,arukodx2,arukodx3

! Automatic arrays:

real, dimension(mza) :: wmt_cor

iu1 = itab_w(iw)%iu(1); iu2 = itab_w(iw)%iu(2); iu3 = itab_w(iw)%iu(3)
iu4 = itab_w(iw)%iu(4); iu5 = itab_w(iw)%iu(5); iu6 = itab_w(iw)%iu(6)
iu7 = itab_w(iw)%iu(7); iu8 = itab_w(iw)%iu(8); iu9 = itab_w(iw)%iu(9)

iw1 = itab_w(iw)%iw(1); iw2 = itab_w(iw)%iw(2); iw3 = itab_w(iw)%iw(3)

fwu4 = itab_w(iw)%fwu(4);  fwu5 = itab_w(iw)%fwu(5); fww1 = itab_w(iw)%fww(1)
fwu6 = itab_w(iw)%fwu(6);  fwu7 = itab_w(iw)%fwu(7); fww2 = itab_w(iw)%fww(2)
fwu8 = itab_w(iw)%fwu(8);  fwu9 = itab_w(iw)%fwu(9); fww3 = itab_w(iw)%fww(3)

vxu1_w = itab_w(iw)%vxu_w(1); vyu1_w = itab_w(iw)%vyu_w(1)
vxu2_w = itab_w(iw)%vxu_w(2); vyu2_w = itab_w(iw)%vyu_w(2)
vxu3_w = itab_w(iw)%vxu_w(3); vyu3_w = itab_w(iw)%vyu_w(3)

dtl2 = 2. * dtlm(itab_w(iw)%mrlw)
ka = lpw(iw)

! Vertical loop over T levels

do k = ka,mza-1

! Coriolis force

   if (icorflg == 1) then
      term = omega2 * rho(k,iw) * volt(k,iw)

      vx_w = vxu1_w * uc(k,iu1) + vxu2_w * uc(k,iu2) + vxu3_w * uc(k,iu3)
      vy_w = vyu1_w * uc(k,iu1) + vyu2_w * uc(k,iu2) + vyu3_w * uc(k,iu3)

      wmt_cor(k) = term * (vy_w * wnx(iw) - vx_w * wny(iw))
   else
      wmt_cor(k) = 0.
   endif

enddo

! Vertical loop over W levels

do k = ka,mza-2

! Diffusive flux coefficients

   arukodx1 = .25 * dniu(iu1)   &
      * (aru(k  ,iu1) * (hkm(k  ,iw1) + hkm(k  ,iw))  &
      +  aru(k+1,iu1) * (hkm(k+1,iw1) + hkm(k+1,iw)))

   arukodx2 = .25 * dniu(iu2)   &
      * (aru(k  ,iu2) * (hkm(k  ,iw2) + hkm(k  ,iw))  &
      +  aru(k+1,iu2) * (hkm(k+1,iw2) + hkm(k+1,iw)))

   arukodx3 = .25 * dniu(iu3)   &
      * (aru(k  ,iu3) * (hkm(k  ,iw3) + hkm(k  ,iw))  &
      +  aru(k+1,iu3) * (hkm(k+1,iw3) + hkm(k+1,iw)))

! Projected velocities

   vproj1 = fwu4 * (uc(k,iu4) + uc(k+1,iu4))  &
          + fwu5 * (uc(k,iu5) + uc(k+1,iu5))  &
          + fww1 * wc(k,iw1)

   vproj2 = fwu6 * (uc(k,iu6) + uc(k+1,iu6))  &
          + fwu7 * (uc(k,iu7) + uc(k+1,iu7))  &
          + fww2 * wc(k,iw2)

   vproj3 = fwu8 * (uc(k,iu8) + uc(k+1,iu8))  &
          + fwu9 * (uc(k,iu9) + uc(k+1,iu9))  &
          + fww3 * wc(k,iw3)

! Update WM tendency from turbulent fluxes and Coriolis force

   wmt(k,iw) = wmt(k,iw)                &

      + arukodx1 * (vproj1 - wc(k,iw))  &
      + arukodx2 * (vproj2 - wc(k,iw))  &
      + arukodx3 * (vproj3 - wc(k,iw))  &
         
      + .5 * (wmt_cor(k) + wmt_cor(k+1))

enddo
   
return
end subroutine veltend_long_w
