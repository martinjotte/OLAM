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
subroutine scalar_transport(umarusc,wmarwsc,rho_old)

use mem_ijtabs,   only: istp, itab_u, jtab_u, jtab_w, mrl_endl
use mem_grid,     only: mza, mua, mwa, zt, zm, dzim, lpu, volt
use misc_coms,    only: io6, dtlm, iparallel
use var_tables,   only: scalar_tab, num_scalar
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

!$ use omp_lib

implicit none

real(kind=8), intent(in) :: umarusc(mza,mua)
real(kind=8), intent(in) :: wmarwsc(mza,mwa)
real(kind=8), intent(in) :: rho_old(mza,mwa)

integer :: j,iw,iw1,iw2,k,mrl,n,kb,iu

real :: dt2

! Automatic arrays:

real :: zwt1(mza)
real :: zwt2(mza)

real :: hcnum_s(mza)
real :: hcnum1_s(mza)

real :: hflx_s(mza,mua,num_scalar)

real :: rpos(mza,mwa,num_scalar)
real :: rneg(mza,mwa,num_scalar)
real :: scp0(mza,mwa,num_scalar)

real, pointer :: scp(:,:)

! Compute vertical advective weights for a stretched grid

do k = 1,mza-1
   zwt1(k) = (zt(k+1) - zm(k)) * dzim(k)
   zwt2(k) =  (zm(k) - zt(k))  * dzim(k)
enddo

! Horizontal advective flux for scalars

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
! JTAB_U(22): ITAB_W(IW)%IU(1:3)
!$omp parallel do private (iu)
do j = 1,jtab_u(22)%jend(mrl); iu = jtab_u(22)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call hflux_s('S',num_scalar,iu,umarusc,wmarwsc,rho_old,hflx_s)

enddo
!$omp end parallel do
endif
call rsub('Uc',22)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw)
do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0('L','S',num_scalar,iw,umarusc,wmarwsc, &
                          zwt1,zwt2,rho_old,hflx_s,rpos,rneg,scp0)

enddo
!$omp end parallel do
endif
call rsub('Wa',26)

! Parallel send/recv of scalars after updates from low-order advective flux

if (iparallel == 1) then
   call mpi_send_w('L',scp0=scp0)
   call mpi_recv_w('L',scp0=scp0)
endif

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
! JTAB_U(21): ITAB_W(IW)%IU(1:9)
!$omp parallel do private (iu,iw1,iw2,dt2,kb,k,hcnum_s,hcnum1_s,n,scp)
do j = 1,jtab_u(21)%jend(mrl); iu = jtab_u(21)%iu(j)
iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
call qsub('U',iu)

   dt2 = 2. * dtlm(itab_u(iu)%mrlu)
   
   kb = lpu(iu)

! Loop over T levels

   do k = kb,mza-1

! Horizontal courant number for scalars (2 in numerator compensates for 
! double mass in denominator since flux used for T cell) 

      hcnum_s(k) = umarusc(k,iu) * dt2  &
         / (rho_old(k,iw1) * volt(k,iw1) + rho_old(k,iw2) * volt(k,iw2))

! Courant number magnitude 1 for low order flux

      hcnum1_s(k) = sign(1.,hcnum_s(k))
         
   enddo
      
! Loop over prognostic scalars

   do n = 1,num_scalar

! Point SCP to scalar array

      scp => scalar_tab(n)%var_p

! Loop over T levels

      hflx_s(1:kb-1,iu,n) = 0.

      do k = kb,mza-1

! Horizontal antidiffusive scalar flux

         hflx_s(k,iu,n) = umarusc(k,iu) * .5  &
!            * (hcnum_s(k) - hcnum1_s(k)) * (scp(k,iw1) - scp(k,iw2))
            * (.75 * hcnum_s(k) - hcnum1_s(k)) * (scp(k,iw1) - scp(k,iw2))

      enddo

   enddo

enddo
!$omp end parallel do
endif
call rsub('Ud',21)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
! JTAB_W(28): IW, ITAB_W(IW)%IW(1:3)
!$omp parallel do private (iw)
do j = 1,jtab_w(28)%jend(mrl); iw = jtab_w(28)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0('M','S',num_scalar,iw,umarusc,wmarwsc, &
                          zwt1,zwt2,rho_old,hflx_s,rpos,rneg,scp0)

enddo
!$omp end parallel do
endif
call rsub('Wb',28)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw)
do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0('H','S',num_scalar,iw,umarusc,wmarwsc, &
                          zwt1,zwt2,rho_old,hflx_s,rpos,rneg,scp0)

enddo
!$omp end parallel do
endif
call rsub('Wc',26)

return
end subroutine scalar_transport

!=========================================================================

subroutine hflux_s(sclr_type,num_sclr,iu,umarusc,wmarwsc,rho_old,hflx_s)

use mem_ijtabs, only: itab_u, itab_w
use mem_basic,  only: thil
use misc_coms,  only: io6, dtsm, dtlm
use mem_grid,   only: mza, mua, mwa, lpu, lpw, volt, arw0, dnv
use var_tables, only: scalar_tab
use mem_para,   only: myrank

implicit none

character(1), intent(in) :: sclr_type
integer, intent(in) :: num_sclr
integer, intent(in) :: iu

real(kind=8), intent(in) :: umarusc(mza,mua)
real(kind=8), intent(in) :: wmarwsc(mza,mwa)
real(kind=8), intent(in) :: rho_old(mza,mwa)

real, intent(inout) :: hflx_s(mza,mua,num_sclr)

integer :: iu1,iu2,iu3,iu4
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: iuc,iud
integer :: iwa,iwb,iwc,iwd

integer :: k,kb,n

real :: tuu1,tuu2,tuu3,tuu4
real :: guw1,guw2,guw3,guw4
real :: tuuc,tuud
real :: guwc,guwd

real :: scq,scq_updn,crossmm,crossww,fww

real :: dt,dt2,dto2
real :: cnum_w

! Automatic arrays

real :: vcnum_s(mza)
real :: htdisp(mza)
real :: c1(mza)
real :: c2(mza)
real :: htgrads(mza)
real :: htgrads_updn(mza)
real :: vterm(mza)
real :: scp_upwind(mza)

real, pointer :: scp(:,:)

iu1 = itab_u(iu)%iu(1); iu2 = itab_u(iu)%iu(2)
iu3 = itab_u(iu)%iu(3); iu4 = itab_u(iu)%iu(4)

iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
iw3 = itab_u(iu)%iw(3); iw4 = itab_u(iu)%iw(4)
iw5 = itab_u(iu)%iw(5); iw6 = itab_u(iu)%iw(6)

tuu1 = itab_u(iu)%tuu(1); tuu2 = itab_u(iu)%tuu(2)
tuu3 = itab_u(iu)%tuu(3); tuu4 = itab_u(iu)%tuu(4)

guw1 = itab_u(iu)%guw(1); guw2 = itab_u(iu)%guw(2)
guw3 = itab_u(iu)%guw(3); guw4 = itab_u(iu)%guw(4)

crossmm = itab_u(iu)%crossmm
crossww = itab_u(iu)%crossww

if (sclr_type == 'S') then
   dt = dtlm(itab_u(iu)%mrlu)
else
   dt = dtsm(itab_u(iu)%mrlu)
endif

dt2  = 2. * dt
dto2 = .5 * dt
   
kb = lpu(iu)
hflx_s(1:kb-1,iu,1:num_sclr) = 0.

! Loop over T levels

do k = kb,mza-1

   if (umarusc(k,iu) >= 0.) then

      iwa = iw2
      iwb = iw1
      iwc = iw3
      iwd = iw4

      iuc = iu1
      iud = iu2
            
      tuuc = tuu1
      tuud = tuu2

      guwc  = guw1
      guwd  = guw2
      
      fww = crossww - .5

   else

      iwa = iw1
      iwb = iw2
      iwc = iw5
      iwd = iw6

      iuc = iu3
      iud = iu4
            
      tuuc = tuu3
      tuud = tuu4

      guwc = guw3
      guwd = guw4

      fww = .5 - crossww

   endif

! Vertical Courant number at T level in IWB column

   vcnum_s(k) = dto2 * (wmarwsc(k-1,iwb) + wmarwsc(k,iwb)) &
              / (rho_old(k,iwb) * volt(k,iwb))

! Horizontal transverse displacement in single timestep at T level in IWB column

   htdisp(k) = dt2 * (tuuc * umarusc(k,iuc) / dnv(iuc)   &
                    + tuud * umarusc(k,iud) / dnv(iud))  &
             / (rho_old(k,iwb) * volt(k,iwb) / arw0(iwb))

   c1(k) = .3333333 * vcnum_s(k) * htdisp(k)
   c2(k) = c1(k) - .5 * htdisp(k)

! Loop over all scalars

   do n = 1,num_sclr

! Point SCP to scalar array

      if (sclr_type == 'S') then
         scp => scalar_tab(n)%var_p
      else
         scp => thil
      endif

! Zero out gradient/difference arrays

      vterm(k)        = 0.0
      htgrads(k)      = 0.0
      htgrads_updn(k) = 0.0
      scq             = scp(k,iwb)

! Horizontal transverse scalar gradient

      if (k >= lpw(iwc) .and. k >= lpw(iwd)) then

         htgrads(k) = guwc * (scp(k,iwc) - scp(k,iwb)) &
                    + guwd * (scp(k,iwd) - scp(k,iwb))

! Correction to upwind point from grid distortion

         scq = scp(k,iwb) + htgrads(k) * (.5 - crossmm) * dnv(iu) &
                          + (scp(k,iwa) - scp(k,iwb)) * fww

      endif

! Terms with vertical gradient component

      if (vcnum_s(k) >= 0.) then

! Case where mean IWB cell motion is upward

         if (k > lpw(iwb)) then

! Vertical difference term

            vterm(k) = .5 * vcnum_s(k) * (scp(k-1,iwb) - scp(k,iwb))

            if (k > lpw(iwc) .and. k > lpw(iwd)) then

! Lower-level horizontal transverse scalar gradient

               htgrads_updn(k) = guwc * (scp(k-1,iwc) - scp(k-1,iwb)) &
                               + guwd * (scp(k-1,iwd) - scp(k-1,iwb))

! Correction to upwind point from grid distortion

               scq_updn = scp(k-1,iwb) &
                        + htgrads_updn(k) * (.5 - crossmm) * dnv(iu) &
                        + (scp(k-1,iwa) - scp(k-1,iwb)) * fww

               vterm(k) = .5 * vcnum_s(k) * (scq_updn - scq)

            endif

         endif

      else

! Case where mean IWB cell motion is downward

         if (k < mza - 1) then

! Vertical difference term

            vterm(k) = .5 * vcnum_s(k) * (scp(k,iwb) - scp(k+1,iwb))

            if (k + 1 >= lpw(iwc) .and. k + 1 >= lpw(iwd)) then

! Upper-level horizontal transverse scalar gradient

               htgrads_updn(k) = guwc * (scp(k+1,iwc) - scp(k+1,iwb)) &
                               + guwd * (scp(k+1,iwd) - scp(k+1,iwb))

! Correction to upwind point from grid distortion

               scq_updn = scp(k+1,iwb) &
                        + htgrads_updn(k) * (.5 - crossmm) * dnv(iu) &
                        + (scp(k+1,iwa) - scp(k+1,iwb)) * fww

               vterm(k) = .5 * vcnum_s(k) * (scq - scq_updn)

            endif

         endif

      endif

! Upwinded upstream value

      scp_upwind(k) = scq + vterm(k) &
                    + c2(k) * htgrads(k) &
                    - c1(k) * htgrads_updn(k)

! SPECIAL---------------------------------------
!      scp_upwind(k) = scp(k,iwb)
!-----------------------------------------------

! Low order flux

      hflx_s(k,iu,n) = umarusc(k,iu) * scp_upwind(k)

   enddo ! n  

enddo ! k

return
end subroutine hflux_s

!===============================================================================

subroutine scalar_transport0(action,sclr_type,num_sclr,iw,umarusc,wmarwsc, &
                             zwt1,zwt2,rho_old,hflx_s,rpos,rneg,scp0)

use var_tables, only: scalar_tab
use mem_ijtabs, only: itab_w
use mem_basic,  only: rho, thil
use mem_tend,   only: thilt
use misc_coms,  only: io6, dtlm, dtsm
use mem_turb,   only: vkh, hkm, sxfer_rk
use mem_grid,   only: mza, mua, mwa, lpw, lsw,  &
                      dniu, volt, aru, arw, volwi, dzim, volti
use massflux,   only: tridiffo
use mem_para,   only: myrank

implicit none

character(1), intent(in) :: action
character(1), intent(in) :: sclr_type

integer, intent(in) :: num_sclr
integer, intent(in) :: iw

real(kind=8), intent(in) :: umarusc(mza,mua)
real(kind=8), intent(in) :: wmarwsc(mza,mwa)
real(kind=8), intent(in) :: rho_old(mza,mwa)

real, intent(in) :: zwt1(mza)
real, intent(in) :: zwt2(mza)

real, intent(inout) :: hflx_s(mza,mua,num_sclr)

real, intent(inout) :: rpos(mza,mwa,num_sclr)
real, intent(inout) :: rneg(mza,mwa,num_sclr)
real, intent(inout) :: scp0(mza,mwa,num_sclr)

real, pointer :: scp(:,:)
real, pointer :: sct(:,:)

integer :: n,k,ka,ks,iu1,iu2,iu3,iw1,iw2,iw3,km,kp

real :: dt,dt2,dti,hdti,hdniu1,hdniu2,hdniu3,diru1,diru2,diru3

! Automatic arrays:

real :: akodz   (mza)
real :: dtomass (mza)
real :: vctr5   (mza)
real :: vctr6   (mza)
real :: vctr7   (mza)
real :: vctr8   (mza)
real :: vctr9   (mza)
real :: vflx_s  (mza)
real :: tmass   (mza)
real :: akodn1  (mza)
real :: akodn2  (mza)
real :: akodn3  (mza)
real :: vcnum_s (mza)
real :: vcnum1_s(mza)
real :: wt1     (mza)
real :: scp_w   (mza)
real :: del_scp (mza)
real :: rhoi    (mza)

real :: cns1(mza)
real :: cns2(mza)
real :: cns3(mza)

real :: cb(mza)
real :: ct(mza)
real :: c1(mza)
real :: c2(mza)
real :: c3(mza)

real :: adv(mza)

real :: scp_upwind(mza)
real :: scp_min(mza)
real :: scp_max(mza)

real :: ppos(mza)
real :: pneg(mza)
real :: qpos(mza)
real :: qneg(mza)

real, parameter :: fourthirds = 4./3.

iu1 = itab_w(iw)%iu(1)
iu2 = itab_w(iw)%iu(2)
iu3 = itab_w(iw)%iu(3)

iw1 = itab_w(iw)%iw(1)
iw2 = itab_w(iw)%iw(2)
iw3 = itab_w(iw)%iw(3)

diru1 = itab_w(iw)%diru(1)
diru2 = itab_w(iw)%diru(2)
diru3 = itab_w(iw)%diru(3)

! Initial computations for this column that need not be repeated for
! each transported scalar field

if (sclr_type == 'S') then
   dt = dtlm(itab_w(iw)%mrlw)
else
   dt = dtsm(itab_w(iw)%mrlw)
endif

dt2 = dt * 2.
dti = 1. / dt
hdti = .5 * dti
ka = lpw(iw)
   
! Vertical loop over T levels - mass in T cells; inverse density

do k = ka,mza-1
   tmass(k) = rho_old(k,iw) * volt(k,iw)
   rhoi(k) = 1. / rho(k,iw)
enddo
   
! Vertical loop over W levels - vertical courant number for scalars (full and unit)

do k = ka,mza-2
   vcnum_s (k) = wmarwsc(k,iw) * dt2 / (tmass(k) + tmass(k+1))
   vcnum1_s(k) = sign(1.,vcnum_s(k))
enddo

vcnum_s(ka-1) = 0.

vcnum_s (mza-1) = 0.
vcnum1_s(mza-1) = 0.
   
if (action == 'L') then

! Vertical loop over T levels

   do k = ka,mza-1

! Set cnum coefficients (one-half courant number)

      cns1(k) = max(0.0_8,diru1 * dt * umarusc(k,iu1) &
              / (tmass(k) + rho_old(k,iw1) * volt(k,iw1)))

      cns2(k) = max(0.0_8,diru2 * dt * umarusc(k,iu2) &
              / (tmass(k) + rho_old(k,iw2) * volt(k,iw3)))

      cns3(k) = max(0.0_8,diru3 * dt * umarusc(k,iu3) &
              / (tmass(k) + rho_old(k,iw2) * volt(k,iw3)))

   enddo

   if (sclr_type == 'S') then

      hdniu1 = .5 * dniu(iu1)
      hdniu2 = .5 * dniu(iu2)
      hdniu3 = .5 * dniu(iu3)

! Vertical loop over T levels 

      do k = ka,mza-1

! delta_t over mass

         dtomass(k) = dt / tmass(k)

! Scalar horizontal turbulent flux coefficients

         akodn1(k) = aru(k,iu1) * (hkm(k,iw1) + hkm(k,iw)) * hdniu1
         akodn2(k) = aru(k,iu2) * (hkm(k,iw2) + hkm(k,iw)) * hdniu2
         akodn3(k) = aru(k,iu3) * (hkm(k,iw3) + hkm(k,iw)) * hdniu3

      enddo

! Vertical loop over W levels

      do k = ka,mza-2

! Prepare for vertical diffusion - Fill tri-diagonal matrix coefficients

         akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
         vctr5(k) = - akodz(k) * dtomass(k)
         vctr7(k) = - akodz(k) * dtomass(k+1)
         vctr6(k) = 1. - vctr5(k) - vctr7(k)
      enddo

   endif ! (sclr_type == 'S')
   
! Loop over scalars

   do n = 1,num_sclr

! Point SCP and SCT to scalar arrays

      if (sclr_type == 'S') then
         scp => scalar_tab(n)%var_p
         sct => scalar_tab(n)%var_t
      else
         scp => thil
         sct => thilt
      endif

      if (sclr_type == 'S') then

! Vertical loop over W levels: 
! Fill r.h.s. for tri-diagonal matrix eqn for vertical diffusion

         do k = ka,mza-2
            vctr8(k) = akodz(k) * (scp(k,iw) - scp(k+1,iw))
         enddo

! Special case for total water sh_w:  apply surface vapor flux

         if (scalar_tab(n)%name == 'SH_W') then

! Vertical loop over T levels that are adjacent to surface

            do ks = 1,lsw(iw)
               k = ka + ks - 1

! Apply surface vapor xfer [kg_vap] directly to SCT [kg_vap / (m^3 s)]

               sct(k,iw) = sct(k,iw) + dti * volti(k,iw) * sxfer_rk(ks,iw)

! Change in SCP from surface xfer

               del_scp(k) = sxfer_rk(ks,iw) / (rho_old(k,iw) * volt(k,iw))

! Zero out sxfer_rk(ks,iw) now that it has been transferred to the atm

               sxfer_rk(ks,iw) = 0.  

            enddo

! Lowest T level that is not adjacent to surface

            del_scp(ka+lsw(iw)) = 0.

! Vertical loop over W levels that are adjacent to surface

            do ks = 1,lsw(iw)
               k = ka + ks - 1

! Change in vctr8 from surface vapor xfer

               vctr8(k) = vctr8(k) + akodz(k) * (del_scp(k) - del_scp(k+1))
            enddo

         endif

! Solve tri-diagonal matrix equation

         if (ka <= mza-2) then
            call tridiffo(mza,ka,mza-2,vctr5,vctr6,vctr7,vctr8,vctr9)
         endif

! Set bottom and top vertical internal turbulent fluxes to zero

         vctr9(ka-1) = 0.
         vctr9(mza-1) = 0.
      
! Vertical loop over T levels

         do k = ka,mza-1

! Add contributions to scalar tendency from horizontal and vertical diffusion

            sct(k,iw) = sct(k,iw) + volti(k,iw) &

               * (akodn1(k) * (scp(k,iw1) - scp(k,iw)) &
               +  akodn2(k) * (scp(k,iw2) - scp(k,iw)) &
               +  akodn3(k) * (scp(k,iw3) - scp(k,iw)) &
               +  vctr9(k-1) - vctr9(k))

         enddo

      endif ! (sclr_type == 'S')

! Vertical loop over T levels

      do k = ka,mza-1

! Horizontally upwinded values for vertical advective fluxes

         scp_upwind(k) = scp(k,iw) &
            + cns1(k) * (scp(k,iw1) - scp(k,iw)) &
            + cns2(k) * (scp(k,iw2) - scp(k,iw)) &
            + cns3(k) * (scp(k,iw3) - scp(k,iw))

      enddo

! Set bottom and top vertical advective fluxes to zero

      vflx_s(ka-1) = 0.
      vflx_s(mza-1) = 0.

! Vertical loop over W levels

      do k = ka,mza-2

! Vertical low order scalar flux from horizontally upwinded values

         wt1(k) = .5 * (vcnum1_s(k) + 1.)
         scp_w(k) = wt1(k) * scp_upwind(k) + (1. - wt1(k)) * scp_upwind(k+1)
         vflx_s(k) = wmarwsc(k,iw) * scp_w(k)

! Alternate form...

!         if (wmarwsc(k,iw) > 0.) then
!            vflx_s(k) = wmarwsc(k,iw) * scp_upwind(k)
!         else
!            vflx_s(k) = wmarwsc(k,iw) * scp_upwind(k+1)
!         endif

      enddo

! Vertical loop over T levels

      do k = ka,mza-1

! Low order advective flux convergence

         adv(k) =         vflx_s(k-1)     &
                -         vflx_s(k  )     &
                + diru1 * hflx_s(k,iu1,n) &
                + diru2 * hflx_s(k,iu2,n) &
                + diru3 * hflx_s(k,iu3,n)

! Update of scp from low-order advection only

         scp0(k,iw,n) = rhoi(k) * (scp(k,iw) * rho_old(k,iw) &
                      + dt * volti(k,iw) * adv(k))

      enddo

      scp0(1:ka-1,iw,n) = scp0(ka,iw,n)
      
   enddo ! n
   
elseif (action == 'M') then

! Loop over scalars

   do n = 1,num_sclr

! Point SCP to scalar arrays

      if (sclr_type == 'S') then
         scp => scalar_tab(n)%var_p
      else
         scp => thil
      endif

! Set bottom and top vertical advective fluxes to zero

      vflx_s(ka-1) = 0.
      vflx_s(mza-1) = 0.

! Vertical loop over W levels

      do k = ka,mza-2

! Vertical antidiffusive scalar flux

!x         vflx_s(k) = wmarwsc(k,iw) * (zwt1(k) + vcnum_s(k) - vcnum1_s(k)) &
!x                   * (scp(k,iw) - scp(k+1,iw))


         wt1(k) = zwt1(k) + vcnum_s(k) * (.5 + vcnum1_s(k) * (.5 - zwt1(k))) &
                - .5 * (vcnum1_s(k) + 1.)
         
         scp_w(k) = wt1(k) * (scp(k,iw) - scp(k+1,iw))

         vflx_s(k) = wmarwsc(k,iw) * scp_w(k)

      enddo

! Vertical loop over T levels

      do k = ka,mza-1
         km = max(ka,k-1)
         kp = min(mza-1,k+1)

! Maximum and minimum of local and neighboring scalar values

         scp_max(k) = max(scp(k ,iw ),scp0(k ,iw ,n) &
                         ,scp(km,iw ),scp0(km,iw ,n) &
                         ,scp(kp,iw ),scp0(kp,iw ,n) &
                         ,scp(k ,iw1),scp0(k ,iw1,n) &
                         ,scp(k ,iw2),scp0(k ,iw2,n) &
                         ,scp(k ,iw3),scp0(k ,iw3,n) )

         scp_min(k) = min(scp(k ,iw ),scp0(k ,iw ,n) &
                         ,scp(km,iw ),scp0(km,iw ,n) &
                         ,scp(kp,iw ),scp0(kp,iw ,n) &
                         ,scp(k ,iw1),scp0(k ,iw1,n) &
                         ,scp(k ,iw2),scp0(k ,iw2,n) &
                         ,scp(k ,iw3),scp0(k ,iw3,n) )

! Compute 6 quantities for flux limiter (ppos,qpos,rpos,pneg,qneg,rneg)

         ppos(k) = max(0.,        vflx_s(k-1))       &
                 + max(0.,      - vflx_s(k  ))       &
                 + max(0.,diru1 * hflx_s(k  ,iu1,n)) &
                 + max(0.,diru2 * hflx_s(k  ,iu2,n)) &
                 + max(0.,diru3 * hflx_s(k  ,iu3,n))

         qpos(k) = (scp_max(k) - scp(k,iw)) * dti * volt(k,iw)

         if (ppos(k) > 1.e-12) then
            rpos(k,iw,n) = min(1.0, qpos(k) / ppos(k))
         else
            rpos(k,iw,n) = 0.
         endif

         pneg(k) = - min(0.,        vflx_s(k-1))       &
                   - min(0.,      - vflx_s(k  ))       &
                   - min(0.,diru1 * hflx_s(k  ,iu1,n)) &
                   - min(0.,diru2 * hflx_s(k  ,iu2,n)) &
                   - min(0.,diru3 * hflx_s(k  ,iu3,n))

         qneg(k) = (scp(k,iw) - scp_min(k)) * dti * volt(k,iw)

         if (pneg(k) > 1.e-12) then
            rneg(k,iw,n) = min(1.0, qneg(k) / pneg(k))
         else
            rneg(k,iw,n) = 0.
         endif

      enddo

      rpos(1:ka-1,iw,n) = 0.
      rneg(1:ka-1,iw,n) = 0.

      rpos(mza,iw,n) = 0.
      rneg(mza,iw,n) = 0.

   enddo ! n

elseif (action == 'H') then

! Loop over scalars

   do n = 1,num_sclr

! Point SCP and SCT to scalar arrays

      if (sclr_type == 'S') then
         scp => scalar_tab(n)%var_p
         sct => scalar_tab(n)%var_t
      else
         scp => thil
         sct => thilt
      endif

! Set bottom and top vertical antidiffusive advective fluxes to zero

      vflx_s(ka-1) = 0.
      vflx_s(mza-1) = 0.

! Vertical loop over W levels

      do k = ka,mza-2

! Vertical antidiffusive scalar flux

!         vflx_s(k) = wmarwsc(k,iw) * (scp(k,iw) - scp(k+1,iw)) &
!                   * (zwt1(k) + .5 * vcnum_s(k) - .5 - .5 * vcnum1_s(k)) &
                   
! (newer)

         wt1(k) = zwt1(k) + vcnum_s(k) * (.5 + vcnum1_s(k) * (.5 - zwt1(k))) &
                - .5 * (vcnum1_s(k) + 1.)
         
         scp_w(k) = wt1(k) * (scp(k,iw) - scp(k+1,iw))


         vflx_s(k) = wmarwsc(k,iw) * scp_w(k)

      enddo

! Vertical loop over T levels

      do k = ka,mza-1

! Limit the 5 fluxes

         if (vflx_s(k-1) >= 0) then
            cb(k) = min(rpos(k,iw,n),rneg(k-1,iw,n))
         else
            cb(k) = min(rpos(k-1,iw,n),rneg(k,iw,n))
         endif

         if (vflx_s(k) <= 0) then
            ct(k) = min(rpos(k,iw,n),rneg(k+1,iw,n))
         else
            ct(k) = min(rpos(k+1,iw,n),rneg(k,iw,n))
         endif

         if (diru1 * hflx_s(k,iu1,n) >= 0) then
            c1(k) = min(rpos(k,iw,n),rneg(k,iw1,n))
         else
            c1(k) = min(rpos(k,iw1,n),rneg(k,iw,n))
         endif

         if (diru2 * hflx_s(k,iu2,n) >= 0) then
            c2(k) = min(rpos(k,iw,n),rneg(k,iw2,n))
         else
            c2(k) = min(rpos(k,iw2,n),rneg(k,iw,n))
         endif

         if (diru3 * hflx_s(k,iu3,n) >= 0) then
            c3(k) = min(rpos(k,iw,n),rneg(k,iw3,n))
         else
            c3(k) = min(rpos(k,iw3,n),rneg(k,iw,n))
         endif

!=================================================================
! SPECIAL
!cb(k) = 1.0
!ct(k) = 1.0
!c1(k) = 1.0
!c2(k) = 1.0
!c3(k) = 1.0
!=================================================================

! Antidiffusive advective flux convergence

         adv(k) = cb(k) *         vflx_s(k-1)     &
                - ct(k) *         vflx_s(k)       &
                + c1(k) * diru1 * hflx_s(k,iu1,n) &
                + c2(k) * diru2 * hflx_s(k,iu2,n) &
                + c3(k) * diru3 * hflx_s(k,iu3,n)

! Update scp from tendency array and antidiffusive fluxes

         if (sclr_type == 'S') then

            scp(k,iw) = scp0(k,iw,n) &
                      + rhoi(k) * dt * (sct(k,iw) + adv(k) * volti(k,iw))

! Update scp0 (which is a provisional updated THIL from long timestep 
! tendencies plus horizontal and EXPLICIT vertical advection)

         else

            scp0(k,iw,n) = scp0(k,iw,n) &
                         + rhoi(k) * dt * (sct(k,iw) + adv(k) * volti(k,iw))
                                
         endif

      enddo

   enddo  ! end loop over scalars

endif  ! (action = 'H')

return
end subroutine scalar_transport0
