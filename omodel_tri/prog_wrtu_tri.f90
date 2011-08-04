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
subroutine prog_wrtu(umarusc,wmarwsc,alpha_press,rhot)

use mem_ijtabs, only: jtab_u, jtab_w, itab_u, istp, itab_w,  &
                      mrl_begs, mrl_ends, mrl_begl, mrl_endl
use mem_basic,  only: rho, thil, wc, press, wmc, ump, umc, uc, theta, sh_w, sh_v
use mem_grid,   only: zt, zm, dzim, lpw, mza, mua, mwa, aru, arw, lpu, lcu,  &
                      volui, volt, unx, volwi, dzm, dzt, glatw, glonw
use misc_coms,  only: io6, iparallel, rinit, rinit8
use consts_coms, only: cpocv, pc1, rdry, rvap

!$ use omp_lib

implicit none

real(kind=8), intent(inout) :: umarusc(mza,mua)
real(kind=8), intent(inout) :: wmarwsc(mza,mwa)

real, intent(inout) :: alpha_press(mza,mwa)
real, intent(in)    :: rhot       (mza,mwa)

integer :: j,iu,iw,k,ka,kb,iwp,mrl,iup

integer :: iw1,iw2

real :: dts,dts2,wmarw2,cnvel,pgf

! automatic arrays

real(kind=8) :: umaru(mza,mua) ! U mom density times U-face surface area [kg/s]
real(kind=8) :: wmarw(mza,mwa) ! W mom density times W-face surface area [kg/s]
real(kind=8) :: rho_s(mza,mwa)
real(kind=8) :: delex_rho(mza,mwa)

real :: hcnum_u(mza,mua)
real :: hcnum_w(mza,mua)
real :: delex_rhothil(mza,mwa,1)

real :: zwt1(mza)
real :: zwt2(mza)

real, parameter :: chi = .1
real, parameter :: chic = 1.5 + chi
real, parameter :: chip = -.5 - chi

! Initialize automatic arrays

umaru         = rinit8
wmarw         = rinit8
delex_rho     = rinit8
hcnum_u       = rinit
hcnum_w       = rinit
zwt1          = rinit
zwt2          = rinit

! Make copy of rho array

rho_s         = rho

! Compute vertical advective weights for stretched grid

do k = 1,mza-1
   zwt1(k) = (zt(k+1) - zm(k)) * dzim(k)
   zwt2(k) =  (zm(k) - zt(k))  * dzim(k)
enddo

! Horizontal mass fluxes for acoustic and long timesteps, and ump update

!$omp parallel do private(k)
do iu = 2,mua
   do k = 1,mza-1
      umaru  (k,iu) = (chic * umc(k,iu) + chip * ump(k,iu)) * aru(k,iu)
      umarusc(k,iu) = umarusc(k,iu) + umaru(k,iu)
      ump    (k,iu) = umc(k,iu)
   enddo
   umaru(mza,iu) = 0.
enddo
!$omp end parallel do

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_W(18): IW, ITAB_W(IW)%IW(1:9)
!$omp parallel do private(iw,kb)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   kb = lpw(iw)

! Explicit (current timestep) vertical mass flux

   wmarw(1:kb-1  ,iw) = 0.
   wmarw(kb:mza-2,iw) = wmc(kb:mza-2,iw) * arw(kb:mza-2,iw)
   wmarw(mza-1   ,iw) = 0.

enddo
!$omp end parallel do
endif
call rsub('Wa',18)

! Horizontal courant numbers for acoustic timestep quantities

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_U(21): IU, ITAB_W(IW)%IU(1:9)
!$omp parallel do private(iu)
do j = 1,jtab_u(21)%jend(mrl); iu = jtab_u(21)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call hcnum(iu,umaru,hcnum_u,hcnum_w,rho_s)

enddo
!$omp end parallel do
endif
call rsub('U',21)

! Compute the ALPHA term for pressure used in prdctw2

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,ka,k)
do j = 1,jtab_w(17)%jend(mrl); iw = jtab_w(17)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   ka = lpw(iw)

   do k = ka,mza-1

      alpha_press(k,iw) = pc1                               &
         * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)  &
         * theta(k,iw) / thil(k,iw)) ** cpocv

   enddo

enddo
!$omp end parallel do
endif

! Explicit rho & thil updates

call progex_rt(umaru,wmarw,rho_s,rhot,delex_rho,delex_rhothil,hcnum_u)

! Update RHO, THIL, WMC

call prog_wrt(umaru,wmarw,wmarwsc,hcnum_w,rho_s,alpha_press &
             ,delex_rho,delex_rhothil)

! Update UMC, UC

call prog_u(umaru,wmarw,hcnum_u,rho_s,zwt1,zwt2)

! WC diagnosis from WMC and RHO - Form from Wenneker    

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,kb,k)
! JTAB_W(18): IW, ITAB_W(IW)%IW(1:9)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   kb = lpw(iw)
   do k = kb,mza-2
      wc(k,iw) = 2. * wmc(k,iw) * dzm(k)  &
               / (dzt(k+1) * rho(k,iw) + dzt(k) * rho(k+1,iw))
   enddo

   wc(1:kb-1,iw) = wc(kb,iw)  ! WC topography BC
   wc(mza-1,iw) = 0.
enddo
!$omp end parallel do
endif
call rsub('Wg',18)

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for WMC, WC, THIL

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,iwp,k)
do j = 1,jtab_w(24)%jend(mrl); iw = jtab_w(24)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 1,mza
      wmc (k,iw) = wmc (k,iwp)
      wc  (k,iw) = wc  (k,iwp)
      thil(k,iw) = thil(k,iwp)
   enddo

enddo
!$omp end parallel do
endif
call rsub('W',24)

return
end subroutine prog_wrtu

!=========================================================================

subroutine hcnum(iu,umaru,hcnum_u,hcnum_w,rho_s)

use mem_ijtabs, only: itab_u, itab_w, itabg_w, itabg_u
use mem_basic,  only: thil
use misc_coms,  only: io6, dtsm
use mem_grid,   only: mza, mua, mwa, lpu, volt, glatu, glonu

implicit none

integer, intent(in) :: iu

real(kind=8), intent(in) :: umaru(mza,mua)

real, intent(out) :: hcnum_u(mza,mua)
real, intent(out) :: hcnum_w(mza,mua)

real(kind=8), intent(in) :: rho_s(mza,mwa)

integer :: iw1,iw2

integer :: k,ka,kb

real :: dts,dts2
real :: cnum_w

! Automatic array

real :: umass(mza)

iw1 = itab_u(iu)%iw(1)
iw2 = itab_u(iu)%iw(2)

kb = lpu(iu)
ka = max(kb-1,2)

dts = dtsm(itab_u(iu)%mrlu)
dts2 = dts * 2.

umass(kb-1) = 0.
hcnum_u(1:kb-1,iu) = 0.

! Loop over T levels

do k = ka,mza-1

! Mass in combined iw1 and iw2 cells

   umass(k) = rho_s(k,iw1) * volt(k,iw1) + rho_s(k,iw2) * volt(k,iw2)

! Horizontal courant number for T

   hcnum_u(k,iu) = umaru(k,iu) * dts / umass(k)

enddo

! Loop over W levels (of U column) for horizontal W Courant number

do k = ka,mza-2
   hcnum_w(k,iu) = (umaru(k,iu) + umaru(k+1,iu)) * dts2  &
                 / (umass(k) + umass(k+1))
enddo

if (kb == 2) then
   hcnum_w(1,iu) = umaru(2,iu) * dts2 / umass(2)
else
   hcnum_w(1:ka-1,iu) = 0.
endif

hcnum_w(mza-1,iu) = umaru(mza-1,iu) * dts2 / umass(mza-1)

return
end subroutine hcnum

!===============================================================================

subroutine progex_rt(umaru,wmarw,rho_s,rhot,delex_rho,thil0,hcnum_u)

use mem_ijtabs,   only: istp, itab_u, itab_w, jtab_u, jtab_w, mrl_begs, itabg_w
use mem_grid,     only: mza, mua, mwa, zt, zm, dzim, lpu, lpw, volt, volti
use mem_basic,    only: rho, thil, umc
use misc_coms,    only: io6, dtsm, iparallel, rinit, rinit8
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

!$ use omp_lib

implicit none

real(kind=8), intent(in)  :: umaru(mza,mua)
real(kind=8), intent(in)  :: wmarw(mza,mwa)
real,         intent(in)  :: hcnum_u(mza,mua)
real(kind=8), intent(in)  :: rho_s(mza,mwa)
real,         intent(in)  :: rhot(mza,mwa)
real(kind=8), intent(out) :: delex_rho(mza,mwa)
real,         intent(out) :: thil0(mza,mwa,1)

integer :: j,iu,iw,iw1,iw2,k,mrl,n,ka,kb
integer :: iu1,iu2,iu3

real :: diru1, diru2, diru3
real :: dts

! Automatic arrays:

real :: zwt1(mza)
real :: zwt2(mza)

real :: hcnum_s(mza)
real :: hcnum1_s(mza)

real :: hflx_s(mza,mua,1)
real :: rpos(mza,mwa,1)
real :: rneg(mza,mwa,1)

real, pointer :: scp(:,:)

! Default initializations

delex_rho = rinit8
thil0     = rinit
hflx_s    = rinit
rpos      = rinit
rneg      = rinit

! Compute vertical advective weights for a stretched grid

do k = 1,mza-1
   zwt1(k) = (zt(k+1) - zm(k)) * dzim(k)
   zwt2(k) =  (zm(k) - zt(k))  * dzim(k)
enddo

! Horizontal advective flux for thil

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_U(22): ITAB_W(IW)%IU(1:3)
!$omp parallel do private(iu)
do j = 1,jtab_u(22)%jend(mrl); iu = jtab_u(22)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call hflux_s('T',1,iu,umaru,wmarw,rho_s,hflx_s)

enddo
!$omp end parallel do
endif
call rsub('Ua',22)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,iu1,iu2,iu3,diru1,diru2,diru3,dts,kb,k)
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
iu1 = itab_w(iw)%iu(1); iu2 = itab_w(iw)%iu(2); iu3 = itab_w(iw)%iu(3)
diru1 = itab_w(iw)%diru(1); diru2 = itab_w(iw)%diru(2); diru3 = itab_w(iw)%diru(3)
!----------------------------------------------------------------------
call qsub('W',iw)

   dts = dtsm(itab_w(iw)%mrlw)

   kb = lpw(iw)

! Loop over T levels

   do k = kb,mza-1

! Change in RHO from horizontal and EXPLICIT vertical advection

      delex_rho(k,iw) = dts * (rhot(k,iw) + volti(k,iw) &
         * (diru1 * umaru(k,iu1)        &
         +  diru2 * umaru(k,iu2)        &
         +  diru3 * umaru(k,iu3)        &
         +  wmarw(k-1,iw) - wmarw(k,iw)))

      rho(k,iw) = rho(k,iw) + delex_rho(k,iw)   
      
   enddo

   call scalar_transport0('L','T',1,iw,umaru,wmarw, &
                          zwt1,zwt2,rho_s,hflx_s,rpos,rneg,thil0)

enddo
!$omp end parallel do
endif
call rsub('Wa',19)

! Parallel send/recv of thil after updates from low-order advective flux

if (iparallel == 1) then
   call mpi_send_w('L',thil0=thil0)
   call mpi_recv_w('L',thil0=thil0)
endif

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_U(21): ITAB_W(IW)%IU(1:9)
!$omp parallel do private(iu,iw1,iw2,kb,k,hcnum_s,hcnum1_s)
do j = 1,jtab_u(21)%jend(mrl); iu = jtab_u(21)%iu(j)
iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
call qsub('U',iu)

   kb = lpu(iu)
   hflx_s(1:kb-1,iu,1) = 0.0

! Loop over T levels

   do k = kb,mza-1

! Horizontal courant number for scalars (factor of 2 compensates
! for double mass in denominator since flux used for T cell) 

      hcnum_s(k) = 2.0 * hcnum_u(k,iu)

! Courant number magnitude 1 for low order flux

      hcnum1_s(k) = sign(1.,hcnum_s(k))

! Horizontal antidiffusive scalar flux

      hflx_s(k,iu,1) = umaru(k,iu) * .5  &
         * (.75 * hcnum_s(k) - hcnum1_s(k)) * (thil(k,iw1) - thil(k,iw2))

   enddo

enddo
!$omp end parallel do
endif
call rsub('Ub',21)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_W(28): IW, ITAB_W(IW)%IW(1:3)
!$omp parallel do private(iw)
do j = 1,jtab_w(28)%jend(mrl); iw = jtab_w(28)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0('M','T',1,iw,umaru,wmarw, &
                          zwt1,zwt2,rho_s,hflx_s,rpos,rneg,thil0)

enddo
!$omp end parallel do
endif
call rsub('Wb',28)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,kb,k)
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call scalar_transport0('H','T',1,iw,umaru,wmarw, &
                          zwt1,zwt2,rho_s,hflx_s,rpos,rneg,thil0)

! Compute DELEX_RHOTHIL and copy into THIL0

   kb = lpw(iw)

   do k = kb,mza-1
      thil0(k,iw,1) = thil0(k,iw,1) * rho(k,iw) - thil(k,iw) * rho_s(k,iw)
   enddo

enddo
!$omp end parallel do
endif
call rsub('Wc',19)

return
end subroutine progex_rt

!===============================================================================

subroutine prog_wrt(umaru,wmarw,wmarwsc,hcnum_w,rho_s,alpha_press &
   ,delex_rho,delex_rhothil)

use mem_ijtabs, only: jtab_u, jtab_w, itab_u, istp, itab_w,  &
                      mrl_begs, mrl_ends, mrl_endl
use mem_basic,  only: rho, thil, wc, press, wmc, ump, umc, uc
use mem_grid,   only: zt, zm, dzim, lpw, mza, mua, mwa, aru, lpu, lcu,  &
                      volui, volt, unx, volwi, dzm, dzt, arw
use misc_coms,  only: io6, iparallel, rinit

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

!$ use omp_lib

implicit none

real(kind=8), intent(in)    :: umaru        (mza,mua) ! U face mass flux [kg/s]
real(kind=8), intent(inout) :: wmarw        (mza,mwa) ! W face mass flux [kg/s]
real(kind=8), intent(inout) :: wmarwsc      (mza,mwa) ! scalar mass flux [kg/s]
real,         intent(in)    :: hcnum_w      (mza,mua)
real(kind=8), intent(in)    :: rho_s        (mza,mwa)
real,         intent(in)    :: alpha_press  (mza,mwa)
real(kind=8), intent(in)    :: delex_rho    (mza,mwa)
real,         intent(inout) :: delex_rhothil(mza,mwa,1)

integer :: j,iu,iw,iw1,iw2,k,ka,kb,iwp,mrl,iup

real :: dts,dts2,cnvel,pgf

! automatic arrays

real :: avflx(mza,mwa)  ! antidiffusive vertical WM flux

real :: ahflx1(mza,mwa) ! antidiffusive horizontal WM fluxes
real :: ahflx2(mza,mwa)
real :: ahflx3(mza,mwa)

real :: rpos(mza,mwa)
real :: rneg(mza,mwa)

real :: wmc0(mza,mwa)
real :: wc0(mza,mwa)

real :: zwt1(mza)
real :: zwt2(mza)

real, parameter :: chi = .1
real, parameter :: chic = 1.5 + chi
real, parameter :: chip = -.5 - chi

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55 ! wmc
real, parameter :: fr = .55 ! rho
real, parameter :: fp = .75 ! press

! default initializations

wc0  = rinit
rpos = rinit
rneg = rinit

! Horizontal loop over primary W/T points to compute wmc0

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw)
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call prog_wrt0('L',iw,umaru,wmarw,rho_s,delex_rho,alpha_press &
      ,delex_rhothil,hcnum_w,zwt1,zwt2,avflx,ahflx1,ahflx2,ahflx3,rpos,rneg &
      ,wmc0,wc0)

enddo
!$omp end parallel do
endif
call rsub('Wd',19)

! WC0 diagnosis from WMC0 and RHO - Form from Wenneker    

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,kb,k)
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   kb = lpw(iw)
   do k = kb,mza-2
      wc0(k,iw) = 2. * wmc0(k,iw) * dzm(k)  &
                / (dzt(k+1) * rho(k,iw) + dzt(k) * rho(k+1,iw))
   enddo

   wc0(1:kb-1,iw) = wc0(kb,iw)  ! WC topography BC
   wc0(mza-1,iw) = 0.
enddo
!$omp end parallel do
endif
call rsub('W',19)

! Parallel send/recv of WC0 after updates from low-order advective flux

if (iparallel == 1) then
   call mpi_send_w('L', wmc0=wc0)
   call mpi_recv_w('L', wmc0=wc0)
endif

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
! JTAB_W(28): IW, ITAB_W(IW)%IW(1:3)
!$omp parallel do private(iw)
do j = 1,jtab_w(28)%jend(mrl); iw = jtab_w(28)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call prog_wrt0('M',iw,umaru,wmarw,rho_s,delex_rho,alpha_press &
      ,delex_rhothil,hcnum_w,zwt1,zwt2,avflx,ahflx1,ahflx2,ahflx3,rpos,rneg &
      ,wmc0,wc0)
      
enddo
!$omp end parallel do
endif
call rsub('We',28)

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw)
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call prog_wrt0('H',iw,umaru,wmarw,rho_s,delex_rho,alpha_press &
      ,delex_rhothil,hcnum_w,zwt1,zwt2,avflx,ahflx1,ahflx2,ahflx3,rpos,rneg &
      ,wmc0,wc0)
      
enddo
!$omp end parallel do
endif
call rsub('Wf',19)

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for PRESS, RHO, WMARW

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,iwp,k)
do j = 1,jtab_w(22)%jend(mrl); iw = jtab_w(22)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = 1,mza
      press(k,iw) = press(k,iwp)
      rho(k,iw) = rho(k,iwp)
   enddo
   
enddo
!$omp end parallel do
endif
call rsub('W',22)

! Parallel send/recv of P group (wmc, press, rho)

if (iparallel == 1) then
   call mpi_send_w('P')
   call mpi_recv_w('P')
endif

! NOW THAT WE HAVE COMMUNICATED WMC, COMPUTE THE VERTICAL MASS FLUX
! AT (t+fw), AND SUM UP THE CONTRIBUTION FOR SCALAR TRANSPORT

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,kb,k)
do j = 1,jtab_w(18)%jend(mrl); iw = jtab_w(18)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   
   kb = lpw(iw)
   do k = kb, mza-2
      wmarw(k,iw)   = (1.0 - fw) * wmarw(k,iw) + fw * wmc(k,iw) * arw(k,iw)
      wmarwsc(k,iw) = wmarwsc(k,iw) + wmarw(k,iw)
   enddo

enddo
!$omp end parallel do
endif
call rsub('W',18)

return
end subroutine prog_wrt

!=========================================================================

subroutine prog_wrt0(action,iw,umaru,wmarw,rho_s,delex_rho &
                    ,alpha_press,delex_rhothil,hcnum_w &
                    ,zwt1,zwt2,avflx,ahflx1,ahflx2,ahflx3,rpos,rneg &
                    ,wmc0,wc0)

use mem_tend,    only: wmt
use mem_ijtabs,  only: itab_w, itab_u
use mem_basic,   only: wmc, rho, thil, wc, uc, theta, press
use misc_coms,   only: io6, initial, dn01d, th01d, &
                       deltax, nxp, mdomain, time8, dtsm
use consts_coms, only: cpocv, gravo2, grav
use mem_grid,    only: mza, mua, mwa, lpw, arw, volt, volti, volwi, dzt, &
                       dzim, xew, zm, glatw, glonw
use mem_rayf,    only: rayfw_distim, rayf_cofw, rayf_distim, rayf_cof
use massflux,    only: tridiffo

implicit none

character(1), intent(in) :: action

integer, intent(in) :: iw

real(kind=8), intent(in) :: umaru(mza,mua)
real(kind=8), intent(in) :: wmarw(mza,mwa)
real(kind=8), intent(in) :: rho_s(mza,mwa)
real(kind=8), intent(in) :: delex_rho(mza,mwa)

real, intent(in) :: alpha_press  (mza,mwa)
real, intent(inout) :: delex_rhothil(mza,mwa,1)
real, intent(in) :: hcnum_w      (mza,mua)

real, intent(in) :: zwt1(mza)
real, intent(in) :: zwt2(mza)

real, intent(inout) :: avflx (mza,mwa)
real, intent(inout) :: ahflx1(mza,mwa)
real, intent(inout) :: ahflx2(mza,mwa)
real, intent(inout) :: ahflx3(mza,mwa)

real, intent(inout) :: rpos(mza,mwa)
real, intent(inout) :: rneg(mza,mwa)
real, intent(inout) :: wmc0 (mza,mwa)
real, intent(inout) :: wc0 (mza,mwa)

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
integer :: iw1,iw2,iw3

real :: diru1,diru2,diru3
real :: fwu4,fwu5,fwu6,fwu7,fwu8,fwu9
real :: fww1,fww2,fww3

integer :: k,ka,kb,km,kp
integer :: k1,k2,k3

real :: dts,dtso2,dts2,dts8,dtsi,flux_rhothil
real :: vol1,vol2

real :: c6,c7,c8,c9,c10,pfrac,sinlat,coslat,theq
real :: fracx, rayfx

real, parameter :: ak = 1. / (40. * 86400.) ! HS expt only
real, parameter :: as = 1. / (4.  * 86400.) ! HS expt only

real :: cnum_w
real :: hcnsclr

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55 ! wmc
real, parameter :: fr = .55 ! rho
real, parameter :: fp = .75 ! press

real, parameter :: pc2 = fp * cpocv
real, parameter :: fourthirds = 4./3.

real :: sc7,sc8,sc9,sc11,sc12,sc13,sc14,sc21,sc22,sc23,sc24,sc25,sc26

! Automatic arrays

real :: del_rhothil(mza)
real :: delex_wm   (mza)
real :: del_wm     (mza)
real :: wmarw2     (mza)
real :: wmtharw    (mza)
real :: del_wmtharw(mza)
real :: wvertflx   (mza)
real :: wc_min     (mza)
real :: wc_max     (mza)
real :: vcnum_w    (mza)
real :: vcnum1_w   (mza)
real :: vflx       (mza)

real :: hcn1u1(mza)
real :: hcn1u2(mza)
real :: hcn1u3(mza)

real :: vproj1(mza)
real :: vproj2(mza)
real :: vproj3(mza)

real :: vproj01(mza)
real :: vproj02(mza)
real :: vproj03(mza)

real :: ppos(mza)
real :: pneg(mza)
real :: qpos(mza)
real :: qneg(mza)

real :: cb(mza)
real :: ct(mza)
real :: c1(mza)
real :: c2(mza)
real :: c3(mza)

real :: hmflx1(mza)
real :: hmflx2(mza)
real :: hmflx3(mza)

real :: hflx1(mza)
real :: hflx2(mza)
real :: hflx3(mza)

real(kind=8) :: del_wmarw(mza)
real(kind=8) :: rhothil  (mza)
real(kind=8) :: press_t  (mza)

real, dimension(mza) :: vctr1,vctr2,vctr5,vctr6,vctr10
real, dimension(mza) :: vctr31,vctr32,vctr33,vctr34,vctr36

iu1 = itab_w(iw)%iu(1); iu2 = itab_w(iw)%iu(2); iu3 = itab_w(iw)%iu(3)
iw1 = itab_w(iw)%iw(1); iw2 = itab_w(iw)%iw(2); iw3 = itab_w(iw)%iw(3)

diru1 = itab_w(iw)%diru(1); diru2 = itab_w(iw)%diru(2); diru3 = itab_w(iw)%diru(3)

kb = lpw(iw)
ka = max(kb-1,2)

dts  = dtsm(itab_w(iw)%mrlw)
dtso2 = .5 * dts
dts2 = 2. * dts
dts8 = 8. * dts
dtsi = 1. / dts

if (action == 'L' .or. action == 'M') then

   iu4 = itab_w(iw)%iu(4)
   iu5 = itab_w(iw)%iu(5)
   iu6 = itab_w(iw)%iu(6)
   iu7 = itab_w(iw)%iu(7)
   iu8 = itab_w(iw)%iu(8)
   iu9 = itab_w(iw)%iu(9)

   fwu4 = itab_w(iw)%fwu(4)
   fwu5 = itab_w(iw)%fwu(5)
   fwu6 = itab_w(iw)%fwu(6)
   fwu7 = itab_w(iw)%fwu(7)
   fwu8 = itab_w(iw)%fwu(8)
   fwu9 = itab_w(iw)%fwu(9)

   fww1 = itab_w(iw)%fww(1)
   fww2 = itab_w(iw)%fww(2)
   fww3 = itab_w(iw)%fww(3)

! Loop over W levels

   do k = kb-1,mza-2

      k1 = max(k,lpw(iw1)-1)
      k2 = max(k,lpw(iw2)-1)
      k3 = max(k,lpw(iw3)-1)

! Projected UC velocity in IW1, IW2, IW3

      vproj1(k) = fwu4 * (uc(k1,iu4) + uc(k1+1,iu4)) &
                + fwu5 * (uc(k1,iu5) + uc(k1+1,iu5)) &
                + fww1 * wc(k1,iw1)

      vproj2(k) = fwu6 * (uc(k2,iu6) + uc(k2+1,iu6)) &
                + fwu7 * (uc(k2,iu7) + uc(k2+1,iu7)) &
                + fww2 * wc(k2,iw2)

      vproj3(k) = fwu8 * (uc(k3,iu8) + uc(k3+1,iu8)) &
                + fwu9 * (uc(k3,iu9) + uc(k3+1,iu9)) &
                + fww3 * wc(k3,iw3)

      if (action == 'M') then

! Projected UC velocity in IW1, IW2, IW3 (updated WC0 replacing WC)

         vproj01(k) = vproj1(k) + fww1 * (wc0(k1,iw1) - wc(k1,iw1))
         vproj02(k) = vproj2(k) + fww2 * (wc0(k2,iw2) - wc(k2,iw2))
         vproj03(k) = vproj3(k) + fww3 * (wc0(k3,iw3) - wc(k3,iw3))
      
      endif

   enddo

endif

if (action == 'L') then

! Vertical loop over T levels

   do k = kb+1,mza-1

      wmarw2(k) = wmarw(k-1,iw) + wmarw(k,iw)

! vertical courant number & sign

      vcnum_w (k) = wmarw2(k) * dtso2 / (rho_s(k,iw) * volt(k,iw))
      vcnum1_w(k) = sign(1.,vcnum_w(k))

! vertical low order advective WC fluxes (use /= wts?)

      vflx(k) = .25 * wmarw2(k) * (wc(k-1,iw) + wc(k,iw)  &
                  + vcnum1_w(k) * (wc(k-1,iw) - wc(k,iw)))

   enddo

!===========================================================================
! MAYBE TRY THE FOLLOWING:

!!   vflx(mza-1) = 0.
!!   avflx(mza-1,iw) = 0.

! MAYBE ALSO TRY:

!!   VOLWI(MZA-2,IW) = 1.0 / (VOLT(MZA-1) + .5 * VOLT(MZA-2))

! AND:  hflux_w(mza-1) added to hflux_w(mza-2)

!===========================================================================

 ! Vertical fluxes at kb are 0 since wc(kb,iw) and wc(kb-1,iw) are set equal  

   vflx(kb) = 0.

! Loop over W levels (downward loop to handle kb-1 level)

   do k = mza-2,kb-1,-1

! Sign of horizontal courant numbers

      hcn1u1(k) = sign(1.,hcnum_w(k,iu1))
      hcn1u2(k) = sign(1.,hcnum_w(k,iu2))
      hcn1u3(k) = sign(1.,hcnum_w(k,iu3))

! Half horizontal mass fluxes at W levels

      hmflx1(k) = .25 * (umaru(k,iu1) + umaru(k+1,iu1))
      hmflx2(k) = .25 * (umaru(k,iu2) + umaru(k+1,iu2))
      hmflx3(k) = .25 * (umaru(k,iu3) + umaru(k+1,iu3))

! Horizontal low-order fluxes (inward directed to IW)

      hflx1(k) = hmflx1(k) * (diru1 * (vproj1(k) + wc(k,iw)) &
                        + hcn1u1(k) * (vproj1(k) - wc(k,iw)) )

      hflx2(k) = hmflx2(k) * (diru2 * (vproj2(k) + wc(k,iw)) &
                        + hcn1u2(k) * (vproj2(k) - wc(k,iw)) )

      hflx3(k) = hmflx3(k) * (diru3 * (vproj3(k) + wc(k,iw)) &
                        + hcn1u3(k) * (vproj3(k) - wc(k,iw)) )

      if (k >= kb) then

! Change in WM(k) from low order advective fluxes

         wmc0(k,iw) = wmc(k,iw) + dts * volwi(k,iw) &
            * (vflx(k) - vflx(k+1) + hflx1(k) + hflx2(k) + hflx3(k))

      else

!!!!!!!!!!!!!!!!!!!!!special   
         go to 54
!!!!!!!!!!!!!!!!!!!!!! end special

! Change in WM(kb) low order advective fluxes at kb-1

         wmc0(kb,iw) = wmc0(kb,iw) + dts * volwi(kb,iw) &
            * (hflx1(k) + hflx2(k) + hflx3(k))

!!!!!!!!!!!!!! special
         54 continue
!!!!!!!!!!!!!!!!end special
      endif

   enddo

   return

endif  ! (action == 'L')

if (action == 'M') then

! Loop over T levels to compute vertical antidiffusive WC fluxes

   do k = kb+1,mza-1

      wmarw2(k) = wmarw(k-1,iw) + wmarw(k,iw)

! vertical courant number & sign

      vcnum_w (k) = wmarw2(k) * dtso2 / (rho_s(k,iw) * volt(k,iw))
      vcnum1_w(k) = sign(1.,vcnum_w(k))

! Vertical antidiffusive advective WC fluxes

      avflx(k,iw) = .25 * wmarw2(k) &
                  * (vcnum_w(k) - vcnum1_w(k)) * (wc(k-1,iw) - wc(k,iw))

   enddo

! Vertical fluxes at kb are 0 since wc(kb,iw) and wc(kb-1,iw) are set equal

   avflx(kb,iw) = 0.

! Loop over W levels (downward loop to handle kb-1 level)

   do k = mza-2,kb-1,-1

! Sign of horizontal courant numbers

      hcn1u1(k) = sign(1.,hcnum_w(k,iu1))
      hcn1u2(k) = sign(1.,hcnum_w(k,iu2))
      hcn1u3(k) = sign(1.,hcnum_w(k,iu3))

! Half horizontal mass fluxes at W levels

      hmflx1(k) = .25 * (umaru(k,iu1) + umaru(k+1,iu1))
      hmflx2(k) = .25 * (umaru(k,iu2) + umaru(k+1,iu2))
      hmflx3(k) = .25 * (umaru(k,iu3) + umaru(k+1,iu3))

! Horizontal antidiffusive fluxes for WM

      ahflx1(k,iw) = hmflx1(k) &
                   * (hcnum_w(k,iu1) - hcn1u1(k)) * (vproj1(k) - wc(k,iw))

      ahflx2(k,iw) = hmflx2(k) &
                   * (hcnum_w(k,iu2) - hcn1u2(k)) * (vproj2(k) - wc(k,iw))

      ahflx3(k,iw) = hmflx3(k) &
                   * (hcnum_w(k,iu3) - hcn1u3(k)) * (vproj3(k) - wc(k,iw))

   enddo

! Loop over W levels

   do k = kb,mza-2
      km = max(ka,k-1)
      kp = min(mza-1,k+1)

! Maximum and minimum of local and projected velocities

      wc_max(k) = max(wc(k,iw) ,wc0(k,iw)  &
                     ,wc(km,iw),wc0(km,iw) &
                     ,wc(kp,iw),wc0(kp,iw) &
                     ,vproj1(k),vproj01(k) &
                     ,vproj2(k),vproj02(k) &
                     ,vproj3(k),vproj03(k) )

      wc_min(k) = min(wc(k,iw) ,wc0(k,iw)  &
                     ,wc(km,iw),wc0(km,iw) &
                     ,wc(kp,iw),wc0(kp,iw) &
                     ,vproj1(k),vproj01(k) &
                     ,vproj2(k),vproj02(k) &
                     ,vproj3(k),vproj03(k) )

! Compute 6 quantities for flux limiter (ppos,qpos,rpos,pneg,qneg,rneg)

      ppos(k) = max(0., avflx(k  ,iw)) &
              + max(0.,-avflx(k+1,iw)) &
              + max(0.,ahflx1(k  ,iw)) &
              + max(0.,ahflx2(k  ,iw)) &
              + max(0.,ahflx3(k  ,iw))

      qpos(k) = (wc_max(k) - wc(k,iw)) * dtsi * volt(k,iw)

      if (ppos(k) > 0.) then
         rpos(k,iw) = min(1.,qpos(k) / ppos(k))
      else
         rpos(k,iw) = 0.
      endif

      pneg(k) = - min(0., avflx(k  ,iw)) &
                - min(0.,-avflx(k+1,iw)) &
                - min(0.,ahflx1(k  ,iw)) &
                - min(0.,ahflx2(k  ,iw)) &
                - min(0.,ahflx3(k  ,iw))

      qneg(k) = (wc(k,iw) - wc_min(k)) * dtsi * volt(k,iw)

      if (pneg(k) > 0.) then
         rneg(k,iw) = min(1.,qneg(k) / pneg(k))
      else
         rneg(k,iw) = 0.
      endif

   enddo

   rpos(1:kb-1,iw) = 0.
   rneg(1:kb-1,iw) = 0.

   rpos(mza-1,iw) = 0.
   rneg(mza-1,iw) = 0.

   return

endif  ! (action == 'M')

!-----------------------------------
! This section done if action = 'H'
!-----------------------------------

! Vertical loop over T levels

do k = kb,mza-1

! RHOTHIL(t) and PRESS(t)

   rhothil(k) = rho_s(k,iw) * thil(k,iw)

   press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

enddo

! Loop over W levels for update of DELEX_WM (downward loop to handle kb-1 level)

do k = mza-2,kb-1,-1

! Limit the 5 fluxes for WMC

   if (ahflx1(k,iw) >= 0) then
      c1(k) = min(rpos(k,iw),rneg(k,iw1))
   else
      c1(k) = min(rpos(k,iw1),rneg(k,iw))
   endif

   if (ahflx2(k,iw) >= 0) then
      c2(k) = min(rpos(k,iw),rneg(k,iw2))
   else
      c2(k) = min(rpos(k,iw2),rneg(k,iw))
   endif

   if (ahflx3(k,iw) >= 0) then
      c3(k) = min(rpos(k,iw),rneg(k,iw3))
   else
      c3(k) = min(rpos(k,iw3),rneg(k,iw))
   endif

   if (k >= kb) then

      if (avflx(k,iw) >= 0) then
         cb(k) = min(rpos(k,iw),rneg(k-1,iw))
      else
         cb(k) = min(rpos(k-1,iw),rneg(k,iw))
      endif

      if (avflx(k+1,iw) <= 0) then
         ct(k) = min(rpos(k,iw),rneg(k+1,iw))
      else
         ct(k) = min(rpos(k+1,iw),rneg(k,iw))
      endif

! Change in WM from EXPLICIT terms (long timestep tendency, 3 horizontal
! advective fluxes, 2 vertical advective fluxes, vertical pgf, gravity)

      delex_wm(k) = wmc0(k,iw) - wmc(k,iw) &

         + dts * (volwi(k,iw) * (wmt(k,iw) &

         + cb(k) * avflx (k  ,iw)  &
         - ct(k) * avflx (k+1,iw)  &
         + c1(k) * ahflx1(k  ,iw)  &
         + c2(k) * ahflx2(k  ,iw)  &
         + c3(k) * ahflx3(k  ,iw)) &

         + dzim(k) * (press_t(k) - press_t(k+1)  &

         - gravo2 * (dzt(k) * rho_s(k,iw) + dzt(k+1) * rho_s(k+1,iw))))

   else
   
!!!!!!!!!!!!!!!!!!!!! special   
      go to 55
!!!!!!!!!!!!!!!!!!!!! end special

      delex_wm(kb) = delex_wm(kb) &

         + dts * volwi(k,iw) &

         * (c1(k) * ahflx1(k  ,iw) &
         +  c2(k) * ahflx2(k  ,iw) &
         +  c3(k) * ahflx3(k  ,iw) )

!!!!!!!!!!!!!! special
      55 continue
!!!!!!!!!!!!!! end special

   endif

enddo

c6  = dts * .5 * fw
c7  = dts * .25 * fw
c8  = dts * fp * cpocv
c9  = dts * (-.5) * fr * grav
c10 = dts * fw

! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN

fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax)

! END SPECIAL

! RAYLEIGH FRICTION ON WM

if (rayfw_distim > 1.e-6) then
   rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-2)
   rayfx = 0.   ! Default: no extra RAYF
   do k = kb,mza-2
      wmc(k,iw) = wmc(k,iw) - dts * max(rayf_cofw(k),rayfx) * wmc(k,iw)
   enddo
endif

! RAYLEIGH FRICTION ON THIL

if (rayf_distim > 1.e-6) then
   if (initial == 1) then   ! HHI case
      rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
      rayfx = 0.   ! Default: no extra RAYF
      do k = kb,mza-1
! Form based on theta alone
         delex_rhothil(k,iw,1) = delex_rhothil(k,iw,1) &
                               + max(rayf_cof(k),rayfx)  & 
                               * dn01d(k) * (th01d(k) - theta(k,iw))
      enddo
   else                     ! LHI/VARI case
      ! Need implementation for LHI/VARI (use vartp for merid. variation?)
   endif
endif

!!!!!!!!!!!!!!! Special - uncomment only for Held-Suarez experiment !!!!!!!

!hs   do k = kb,mza-1                                                     !hs
                     
!hs      pfrac = 1.e-5 * press(k,iw)                                      !hs
!hs      sinlat = sin(glatw(iw)*pio180)                                   !hs
!hs      coslat = cos(glatw(iw)*pio180)                                   !hs
!hs      theq = max(200.*pfrac**(-.285714)  &                             !hs
!hs             ,315. - 60. * sinlat**2 - 10. * log(pfrac) * coslat**2)   !hs

!hs      rayf_cof(k) = ak  &                                              !hs
!hs         + (as - ak) * max(0.,(press(k,iw)-7.e4)/3.e4) * coslat**4     !hs
            
!hs      delex_rhothil(k,iw,1) = delex_rhothil(k,iw,1) &                  !hs
!hs                            + dts * rayf_cof(k)  &                     !hs
!hs                            * rho_s(k,iw) * (theq - theta(k,iw))       !hs
               
!hs   enddo                                                               !hs
!!!!!!!!!!!!!!!!!!!!!! End special !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fill matrix coefficients for implicit update of WM

do k = kb,mza-1
   vctr1(k)  = wc(k,iw) + wc(k-1,iw)     ! T pts
   vctr2(k)  = thil(k,iw) + thil(k+1,iw) ! W pts
   vctr5(k)  = press_t(k) / rhothil(k)   ! T pts
   vctr6(k)  = c6 * volti(k,iw)          ! T pts
   vctr10(k) = c10 * volti(k,iw)         ! T pts
enddo
vctr2(kb-1) = vctr2(kb)

do k = kb,mza-2
         
   sc7  = c7   * volwi(k,iw) ! W pts
   sc8  = c8   * dzim(k)     ! W pts
   sc9  = c9   * dzim(k)     ! W pts
   sc11 = sc8  * vctr5(k)    ! W pts
   sc12 = sc8  * vctr5(k+1)  ! W pts
   sc13 = sc9  * dzt(k)      ! W pts
   sc14 = sc9  * dzt(k+1)    ! W pts
      
   sc21 = sc7  * vctr1(k)    ! W pts
   sc22 = sc7  * vctr1(k+1)  ! W pts
   sc23 = sc11 * vctr6(k)    ! W pts
   sc24 = sc12 * vctr6(k+1)  ! W pts
   sc25 = sc13 * vctr10(k)   ! W pts
   sc26 = sc14 * vctr10(k+1) ! W pts
      
   vctr32(k) = 1. + arw(k,iw) &
             * (sc22 - sc21 + vctr2(k) * (sc23 + sc24) + sc25 - sc26)
         
   vctr31(k) = - arw(k-1,iw) * (sc21 + sc23 * vctr2(k-1) + sc25)
      
   vctr33(k) =   arw(k+1,iw) * (sc22 - sc24 * vctr2(k+1) + sc26) 
     
   vctr34(k) = delex_wm(k) &
             + sc11 * delex_rhothil(k,iw,1) - sc12 * delex_rhothil(k+1,iw,1) &
             + sc13 * delex_rho(k,iw)       + sc14 * delex_rho(k+1,iw)

enddo

! Solve implicit tri-diagonal matrix equation for delta WM (del_wm)

if (kb <= mza-2) then
   call tridiffo(mza,kb,mza-2,vctr31,vctr32,vctr33,vctr34,del_wm)
endif
       
! Vertical loop over W points

do k = kb,mza-2

! Update WMC from (t) to (t+1) due to change in WM (del_wm)

   wmc(k,iw) = wmc(k,iw) + del_wm(k)

! Change in vertical mass and heat fluxes from (t) to (t+fw)

   del_wmarw(k)   = fw * del_wm(k) * arw(k,iw)
   del_wmtharw(k) = del_wmarw(k) * .5 * (thil(k,iw) + thil(k+1,iw))

enddo

! Set top and bottom values in mass-flux-change and heat-flux-change arrays to zero

del_wmarw(kb-1)  = 0.
del_wmarw(mza-1) = 0.

del_wmtharw(kb-1)  = 0.
del_wmtharw(mza-1) = 0.

! Vertical loop over T points

do k = kb,mza-1

! Change of rho from (*) to (t+1)

   rho(k,iw) = rho(k,iw) + dts * volti(k,iw) * (del_wmarw(k-1) - del_wmarw(k))

! Change of rhothil from (t) to (t+1)

   del_rhothil(k) = delex_rhothil(k,iw,1)  &
      + dts * volti(k,iw) * (del_wmtharw(k-1) - del_wmtharw(k))

! Update pressure from (t) to (t+tp)

   press(k,iw) = press_t(k) + pc2 * vctr5(k) * del_rhothil(k)

! Update thil from (t) to (t+1)

   thil(k,iw) = (rhothil(k) + del_rhothil(k)) / rho(k,iw)

enddo

return
end subroutine prog_wrt0

!===============================================================================

subroutine prog_u(umaru,wmarw,hcnum_u,rho_s,zwt1,zwt2)

use mem_ijtabs, only: jtab_u, jtab_w, itab_u, istp, itab_w,  &
                      mrl_begs, mrl_ends, mrl_endl
use mem_basic,  only: rho, wc, press, wmc, ump, umc, uc
use mem_grid,   only: zt, zm, dzim, lpw, mza, mua, mwa, aru, lpu, lcu,  &
                      volui, volt, unx, volwi, dzm, dzt
use misc_coms,  only: io6, iparallel, rinit

use olam_mpi_atm, only: mpi_send_u, mpi_recv_u

!$ use omp_lib

implicit none

real(kind=8), intent(in) :: umaru(mza,mua) ! U face mass flux [kg/s]
real(kind=8), intent(in) :: wmarw(mza,mwa) ! W face mass flux [kg/s]
real(kind=8), intent(in) :: rho_s(mza,mwa)

real, intent(in) :: hcnum_u(mza,mua)

real, intent(in) :: zwt1(mza)
real, intent(in) :: zwt2(mza)

integer :: j,iu,iw,k,ka,kb,iwp,mrl,iup

integer :: iw1,iw2

real :: dts,dts2,wmarw2,cnvel,pgf

! automatic arrays

real :: umc0(mza,mua)  ! update from low order advective flux
real :: uc0(mza,mua)  ! update from low order advective flux

real :: avflx(mza,mua)  ! antidiffusive vertical UM flux

real :: ahflx1(mza,mua) ! antidiffusive horizontal UM fluxes
real :: ahflx2(mza,mua)
real :: ahflx3(mza,mua)
real :: ahflx4(mza,mua)

real :: rpos(mza,mua)
real :: rneg(mza,mua)

! default initializations

uc0  = rinit
rpos = rinit
rneg = rinit

! Update UMC0 from low order flux and compute antidiffusive flux
! Diagnose UC0

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu)
do j = 1,jtab_u(16)%jend(mrl); iu = jtab_u(16)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call prog_u0('L',iu,umaru,wmarw,rho_s,hcnum_u, &
                avflx,ahflx1,ahflx2,ahflx3,ahflx4,rpos,rneg,umc0,uc0)

enddo
!$omp end parallel do
endif
call rsub('Ua',16)

! Parallel send/recv of UC0 after updates from low-order advective flux

if (iparallel == 1) then
   call mpi_send_u('L',uc0=uc0)
   call mpi_recv_u('L',uc0=uc0)
endif

! Update MAX & MIN UC values for flux limiter

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu)
do j = 1,jtab_u(16)%jend(mrl); iu = jtab_u(16)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call prog_u0('M',iu,umaru,wmarw,rho_s,hcnum_u, &
                avflx,ahflx1,ahflx2,ahflx3,ahflx4,rpos,rneg,umc0,uc0)

enddo
!$omp end parallel do
endif
call rsub('Ub',16)

if (iparallel == 1) then
   call mpi_send_u('R',rpos=rpos,rneg=rneg)  ! Send rpos, rneg (for flux limiter)
   call mpi_recv_u('R',rpos=rpos,rneg=rneg)  ! Receive rpos, rneg (for flux limiter)
endif

! Horizontal loop over U points to update UMC

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
!   call uwcomp(0,0,0,0.,0.,0.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(iu)
do j = 1,jtab_u(16)%jend(mrl); iu = jtab_u(16)%iu(j)
!----------------------------------------------------------------------
call qsub('U',iu)

   call prog_u0('H',iu,umaru,wmarw,rho_s,hcnum_u, &
                avflx,ahflx1,ahflx2,ahflx3,ahflx4,rpos,rneg,umc0,uc0)

enddo
!$omp end parallel do
endif
call rsub('Uc',16)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
!   call uwcomp(2,0,0,0.,0.,0.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! LATERAL BOUNDARY CONDITION ONLY FOR LIMITED-AREA-DOMAIN MODEL CONFIGURATION
! Copy LBC for UMC, UC

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iu,iup,k) 
do j = 1,jtab_u(18)%jend(mrl); iu = jtab_u(18)%iu(j)
   iup = itab_u(iu)%iup
!----------------------------------------------------------------------
call qsub('U',iu)

   do k = 1,mza-1
      umc(k,iu) = umc(k,iup)
      uc(k,iu)  = uc(k,iup)
   enddo

enddo
!$omp end parallel do 
endif
call rsub('U',18)

if (iparallel == 1) then
   call mpi_send_u('U')  ! Send U group
endif

return
end subroutine prog_u

!============================================================================

subroutine prog_u0(action,iu,umaru,wmarw,rho_s,hcnum_u, &
                   avflx,ahflx1,ahflx2,ahflx3,ahflx4,rpos,rneg,umc0,uc0)

use mem_tend,    only: umt
use mem_ijtabs,  only: itab_u
use mem_basic,   only: uc, wc, press, ump, umc, rho
use misc_coms,   only: io6, dtsm, initial, mdomain, u01d, v01d, dn01d, &
                       deltax, nxp
use consts_coms, only: erad
use mem_grid,    only: lpu, lcu, volt, aru, volui, xeu, yeu, zeu,  &
                       unx, uny, unz, mza, mua, mwa, dnu, dniu, dnv, arw0, zt
use mem_rayf,    only: rayf_distim, rayf_cof

implicit none

character(1), intent(in) :: action
integer, intent(in) :: iu

real(kind=8), intent(in) :: umaru(mza,mua)
real(kind=8), intent(in) :: wmarw(mza,mwa)
real(kind=8), intent(in) :: rho_s(mza,mwa)

real, intent(in) :: hcnum_u(mza,mua)

real, intent(inout) :: avflx(mza,mua)

real, intent(inout) :: ahflx1(mza,mua)
real, intent(inout) :: ahflx2(mza,mua)
real, intent(inout) :: ahflx3(mza,mua)
real, intent(inout) :: ahflx4(mza,mua)

real, intent(inout) :: rpos(mza,mua)
real, intent(inout) :: rneg(mza,mua)

real, intent(inout) :: umc0(mza,mua)
real, intent(inout) :: uc0(mza,mua)

integer :: k,ka,kb,km,kp

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12
integer :: iw1,iw2,iw3,iw4,iw5,iw6

real :: diru1,diru2,diru3,diru4

real :: fuu5,fuu6,fuu7,fuu8,fuu9,fuu10,fuu11,fuu12
real :: fuw3,fuw4,fuw5,fuw6
real :: pgc12,pgc45,pgc63,pgc12b,pgc45b,pgc12c,pgc63c,pgc12d

real :: dts,dts2,raxis,uv01dr,uv01dx,uv01dy,uv01dz,ucref,uc2
real :: fracx, rayfx, div1, div2, ump_eqdiv

real :: rho_us,qposs

!hs real :: pressloc

! Automatic arrays

real :: umass(mza)
real :: umt_rayf(mza)
real :: pgf(mza)

real :: wmarw2(mza)
real :: vcnum_u(mza)
real :: vcnum1_u(mza)
real :: vflx(mza)

real :: rho_u(mza)

real :: uc_max(mza)
real :: uc_min(mza)

real :: ppos(mza)
real :: pneg(mza)
real :: qpos(mza)
real :: qneg(mza)

real :: vproj1(mza)
real :: vproj2(mza)
real :: vproj3(mza)
real :: vproj4(mza)

real :: vproj01(mza)
real :: vproj02(mza)
real :: vproj03(mza)
real :: vproj04(mza)

real :: hflx1(mza)
real :: hflx2(mza)
real :: hflx3(mza)
real :: hflx4(mza)

real :: hcn1u1(mza)
real :: hcn1u2(mza)
real :: hcn1u3(mza)
real :: hcn1u4(mza)

real :: cb(mza)
real :: ct(mza)
real :: c1(mza)
real :: c2(mza)
real :: c3(mza)
real :: c4(mza)

! neighbor indices and coefficients for this point in the U stencil

iu1  = itab_u(iu)%iu(1)
iu2  = itab_u(iu)%iu(2)
iu3  = itab_u(iu)%iu(3)
iu4  = itab_u(iu)%iu(4) 
iu5  = itab_u(iu)%iu(5)
iu6  = itab_u(iu)%iu(6)
iu7  = itab_u(iu)%iu(7)
iu8  = itab_u(iu)%iu(8)
iu9  = itab_u(iu)%iu(9)
iu10 = itab_u(iu)%iu(10)
iu11 = itab_u(iu)%iu(11)
iu12 = itab_u(iu)%iu(12)

iw1 = itab_u(iu)%iw(1)
iw2 = itab_u(iu)%iw(2)
iw3 = itab_u(iu)%iw(3)
iw4 = itab_u(iu)%iw(4)
iw5 = itab_u(iu)%iw(5)
iw6 = itab_u(iu)%iw(6)

diru1 = itab_u(iu)%diru(1)
diru2 = itab_u(iu)%diru(2)
diru3 = itab_u(iu)%diru(3)
diru4 = itab_u(iu)%diru(4)

dts = dtsm(itab_u(iu)%mrlu)
dts2 = dts * 2.

ka = lcu(iu)
kb = lpu(iu)

!k if (iu == 25402) print*, 'ka,kb at iu=25402 ',ka,kb

if (action == 'L' .or. action == 'M') then

   fuu5  = itab_u(iu)%fuu(5)
   fuu6  = itab_u(iu)%fuu(6)
   fuu7  = itab_u(iu)%fuu(7)
   fuu8  = itab_u(iu)%fuu(8)
   fuu9  = itab_u(iu)%fuu(9)
   fuu10 = itab_u(iu)%fuu(10)
   fuu11 = itab_u(iu)%fuu(11)
   fuu12 = itab_u(iu)%fuu(12)

   fuw3 = itab_u(iu)%fuw(3)
   fuw4 = itab_u(iu)%fuw(4)
   fuw5 = itab_u(iu)%fuw(5)
   fuw6 = itab_u(iu)%fuw(6)

! Vertical loop over U points

   do k = ka,mza-1

! Projected velocities of horizontal neighbor points onto current IU point

      vproj1(k) = fuu5  * uc(k,iu5)               &
                + fuu6  * uc(k,iu6)               &
                + fuw3  * (wc(k,iw3) + wc(k-1,iw3))

      vproj2(k) = fuu7  * uc(k,iu7)               &
                + fuu8  * uc(k,iu8)               &
                + fuw4  * (wc(k,iw4) + wc(k-1,iw4))

      vproj3(k) = fuu9  * uc(k,iu9)               &
                + fuu10 * uc(k,iu10)              &
                + fuw5  * (wc(k,iw5) + wc(k-1,iw5))

      vproj4(k) = fuu11 * uc(k,iu11)              &
                + fuu12 * uc(k,iu12)              &
                + fuw6  * (wc(k,iw6) + wc(k-1,iw6))

   enddo
   
endif ! (action == 'L' .or. action == 'M')

if (action == 'L') then
 
! Vertical loop over T levels

   do k = ka,mza-1

! Mass in UM control volume

      umass(k) = rho_s(k,iw1) * volt(k,iw1) + rho_s(k,iw2) * volt(k,iw2)

   enddo

! Vertical loop over W levels

   do k = ka,mza-2

! Vertical mass flux at top of U control volume

      wmarw2(k) = wmarw(k,iw1) + wmarw(k,iw2)
      
! Vertical courant number (double time in numerator, double mass in denom)

      vcnum_u (k) = wmarw2(k) * dts2 / (umass(k) + umass(k+1))
      vcnum1_u(k) = sign(1.,vcnum_u(k))

! vertical low order advective UC flux (use /= wts?)

      vflx(k) = .5 * wmarw2(k) * (uc(k,iu) + uc(k+1,iu) &
                 + vcnum1_u(k) * (uc(k,iu) - uc(k+1,iu)))

! vertical antidiffusive advective UC flux

      avflx(k,iu) = .5 * wmarw2(k) &
                  * (vcnum_u(k) - vcnum1_u(k)) * (uc(k,iu) - uc(k+1,iu))

   enddo

! Zero out top and bottom vertical fluxes

   vflx(ka-1) = 0.
   avflx(ka-1,iu) = 0.

   vflx(mza-1) = 0.
   avflx(mza-1,iu) = 0.

! Vertical loop over U points

   do k = ka,mza-1

! Courant number magnitude 1 for low order fluxes

      hcn1u1(k) = sign(1.,hcnum_u(k,iu1))
      hcn1u2(k) = sign(1.,hcnum_u(k,iu2))
      hcn1u3(k) = sign(1.,hcnum_u(k,iu3))
      hcn1u4(k) = sign(1.,hcnum_u(k,iu4))

! Horizontal low-order fluxes (inward directed to IU)

      hflx1(k) = umaru(k,iu1) * .5 * (diru1 * (vproj1(k) + uc(k,iu))  &
                                + hcn1u1(k) * (vproj1(k) - uc(k,iu)))

      hflx2(k) = umaru(k,iu2) * .5 * (diru2 * (vproj2(k) + uc(k,iu))  &
                                + hcn1u2(k) * (vproj2(k) - uc(k,iu)))

      hflx3(k) = umaru(k,iu3) * .5 * (diru3 * (vproj3(k) + uc(k,iu))  &
                                + hcn1u3(k) * (vproj3(k) - uc(k,iu)))

      hflx4(k) = umaru(k,iu4) * .5 * (diru4 * (vproj4(k) + uc(k,iu))  &
                                + hcn1u4(k) * (vproj4(k) - uc(k,iu)))

! Horizontal antidiffusive fluxes

      ahflx1(k,iu) = umaru(k,iu1) * .5 * (hcnum_u(k,iu1) - hcn1u1(k)) &
                   * (vproj1(k) - uc(k,iu))
      ahflx2(k,iu) = umaru(k,iu2) * .5 * (hcnum_u(k,iu2) - hcn1u2(k)) &
                   * (vproj2(k) - uc(k,iu))
      ahflx3(k,iu) = umaru(k,iu3) * .5 * (hcnum_u(k,iu3) - hcn1u3(k)) &
                   * (vproj3(k) - uc(k,iu))
      ahflx4(k,iu) = umaru(k,iu4) * .5 * (hcnum_u(k,iu4) - hcn1u4(k)) &
                   * (vproj4(k) - uc(k,iu))

! Update UMC from low order advection

      umc0(k,iu) = umc(k,iu) + dts * volui(k,iu) &
        * (vflx(k-1) - vflx(k) + hflx1(k) + hflx2(k) + hflx3(k) + hflx4(k))

! UC diagnosis from UMC and RHO - Form from Wenneker      

      uc0(k,iu) = umc0(k,iu) / (volui(k,iu)  &
         * (volt(k,iw2) * rho(k,iw1) + volt(k,iw1) * rho(k,iw2)))

   enddo

   uc0(1:ka-1,iu) = uc0(ka,iu)

   return

endif ! (action == 'L')

if (action == 'M') then

! Vertical loop over U points

   do k = ka,mza-1
      km = max(ka,k-1)
      kp = min(mza-1,k+1)

! Projected velocities of horizontal neighbor points onto current IU point

      vproj01(k) = fuu5  * uc0(k,iu5)               &
                 + fuu6  * uc0(k,iu6)               &
                 + fuw3  * (wc(k,iw3) + wc(k-1,iw3))

      vproj02(k) = fuu7  * uc0(k,iu7)               &
                 + fuu8  * uc0(k,iu8)               &
                 + fuw4  * (wc(k,iw4) + wc(k-1,iw4))

      vproj03(k) = fuu9  * uc0(k,iu9)               &
                 + fuu10 * uc0(k,iu10)              &
                 + fuw5  * (wc(k,iw5) + wc(k-1,iw5))

      vproj04(k) = fuu11 * uc0(k,iu11)              &
                 + fuu12 * uc0(k,iu12)              &
                 + fuw6  * (wc(k,iw6) + wc(k-1,iw6))

! Maximum and minimum of local and projected velocities

      uc_max(k) = max(uc(k,iu) ,uc0(k,iu)  &
                     ,uc(km,iu),uc0(km,iu) &
                     ,uc(kp,iu),uc0(kp,iu) &
                     ,vproj1(k),vproj01(k) &
                     ,vproj2(k),vproj02(k) &
                     ,vproj3(k),vproj03(k) &
                     ,vproj4(k),vproj04(k) )

      uc_min(k) = min(uc(k,iu) ,uc0(k,iu)  &
                     ,uc(km,iu),uc0(km,iu) &
                     ,uc(kp,iu),uc0(kp,iu) &
                     ,vproj1(k),vproj01(k) &
                     ,vproj2(k),vproj02(k) &
                     ,vproj3(k),vproj03(k) &
                     ,vproj4(k),vproj04(k) )

! Compute 6 quantities for flux limiter

      rho_u(k) = volui(k,iu)  &
               * (volt(k,iw2) * rho(k,iw1) + volt(k,iw1) * rho(k,iw2))

      rho_us = volui(k,iu)  &
             * (volt(k,iw2) * rho_s(k,iw1) + volt(k,iw1) * rho_s(k,iw2))

      ppos(k) = max(0., avflx(k-1,iu)) &
              + max(0.,-avflx(k  ,iu)) &
              + max(0.,ahflx1(k  ,iu)) &
              + max(0.,ahflx2(k  ,iu)) &
              + max(0.,ahflx3(k  ,iu)) &
              + max(0.,ahflx4(k  ,iu))

      qpos(k) = (uc_max(k) * rho_u(k) - umc0(k,iu)) / (dts * volui(k,iu))

      qposs = (uc_max(k) * rho_us - umc0(k,iu)) / (dts * volui(k,iu))

      if (ppos(k) > 0.) then
         rpos(k,iu) = min(1.,qpos(k) / ppos(k))
      else
         rpos(k,iu) = 0.
      endif

      pneg(k) = - min(0., avflx(k-1,iu)) &
                - min(0.,-avflx(k  ,iu)) &
                - min(0.,ahflx1(k  ,iu)) &
                - min(0.,ahflx2(k  ,iu)) &
                - min(0.,ahflx3(k  ,iu)) &
                - min(0.,ahflx4(k  ,iu))

      qneg(k) = (umc0(k,iu) - uc_min(k) * rho_u(k)) / (dts * volui(k,iu))

      if (pneg(k) > 0.) then
         rneg(k,iu) = min(1.,qneg(k) / pneg(k))
      else
         rneg(k,iu) = 0.
      endif

   enddo

   rneg(1:ka-1,iu) = 0.
   rpos(1:ka-1,iu) = 0.

   rneg(mza,iu) = 0.
   rpos(mza,iu) = 0.

   return

endif ! (action == 'M')

!-----------------------------------
! This section done if action = 'H'
!-----------------------------------

! RAYLEIGH FRICTION ON UM

umt_rayf(:) = 0.

if (rayf_distim > 1.e-6) then

! Vertical loop over U points

   do k = kb,mza-1

      if (initial == 1) then      ! HHI case

! Must rotate reference wind to local UC orientation

         if (mdomain <= 1) then  ! Model uses "earth" coordinates
            raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

            if (raxis > 1.e3) then
               uv01dr = -v01d(k) * zeu(iu) / erad  ! radially outward from axis

               uv01dx = (-u01d(k) * yeu(iu) + uv01dr * xeu(iu)) / raxis 
               uv01dy = ( u01d(k) * xeu(iu) + uv01dr * yeu(iu)) / raxis 
               uv01dz =   v01d(k) * raxis / erad 

               ucref = uv01dx * unx(iu) + uv01dy * uny(iu) + uv01dz * unz(iu)
            else
               ucref = 0.
            endif
         else
            ucref = u01d(k) * unx(iu) + v01d(k) * uny(iu)
         endif

! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN

         fracx = abs(xeu(iu)) / (real(nxp-1) * .866 * deltax)
         rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
         rayfx = 0.   ! Default: no extra RAYF
         
! END SPECIAL

         umt_rayf(k) = max(rayf_cof(k),rayfx) * dn01d(k) * (ucref - uc(k,iu))

      else                     ! LHI/VARI case

! HORIZONTAL DIVERGENCE DAMPING

! UMP value that would equalize horizontal divergence in IW1 and IW2 cells
! (assuming no blockage by topography)

         ump_eqdiv = (arw0(iw2) * (ump(k,iu1) * dnv(iu1) * diru1   &
                                 + ump(k,iu2) * dnv(iu2) * diru2)  &
                   -  arw0(iw1) * (ump(k,iu3) * dnv(iu3) * diru3   &
                                 + ump(k,iu4) * dnv(iu4) * diru4)) &
                   / (dnv(iu) * (arw0(iw1) + arw0(iw2)))

         umt_rayf(k) = rayf_cof(k) * (ump_eqdiv - ump(k,iu))

         ! Need implementation for LHI/VARI (use vartp for merid. variation?)
         ! Ok to do this in veltend_long???
      endif

   enddo

endif ! (rayf_distim > 1.e-6)

! PRESSURE GRADIENT FORCE ON UMC

pgc12  = itab_u(iu)%pgc12
pgc45  = itab_u(iu)%pgc45
pgc63  = itab_u(iu)%pgc63
pgc12b = itab_u(iu)%pgc12b
pgc45b = itab_u(iu)%pgc45b
pgc12c = itab_u(iu)%pgc12c
pgc63c = itab_u(iu)%pgc63c
pgc12d = itab_u(iu)%pgc12d

pgf(ka:kb-1) = 0.

! Vertical loop over T levels

do k = kb,mza-1

   if ((lpu(iu1) > k .or. lpu(iu4) > k) .and.  &
       (lpu(iu2) > k .or. lpu(iu3) > k)) then

! 2-point stencil

      pgf(k) = pgc12d * (press(k,iw1) - press(k,iw2))

   elseif (lpu(iu2) > k .or. lpu(iu3) > k) then

! 4-point stencil with pts 6-3

      pgf(k) = (pgc12c * (press(k,iw1) - press(k,iw2))  &
             +  pgc63c * (press(k,iw6) - press(k,iw3))  )

   elseif (lpu(iu4) > k .or. lpu(iu1) > k) then

! 4-point stencil with pts 4-5

      pgf(k) = (pgc12b * (press(k,iw1) - press(k,iw2))  &
             +  pgc45b * (press(k,iw4) - press(k,iw5))  )

   else

! 6-point stencil

      pgf(k) = (pgc12 * (press(k,iw1) - press(k,iw2))  &
             +  pgc45 * (press(k,iw4) - press(k,iw5))  &
             +  pgc63 * (press(k,iw6) - press(k,iw3))  )

   endif

enddo

! Vertical loop over U points

do k = ka,mza-1

! Limit the 6 fluxes

   if (avflx(k-1,iu) >= 0) then
      cb(k) = min(rpos(k,iu),rneg(k-1,iu))
   else
      cb(k) = min(rpos(k-1,iu),rneg(k,iu))
   endif

   if (avflx(k,iu) <= 0) then
      ct(k) = min(rpos(k,iu),rneg(k+1,iu))
   else
      ct(k) = min(rpos(k+1,iu),rneg(k,iu))
   endif

   if (ahflx1(k,iu) >= 0) then
      c1(k) = min(rpos(k,iu),rneg(k,iu6))
   else
      c1(k) = min(rpos(k,iu6),rneg(k,iu))
   endif

   if (ahflx2(k,iu) >= 0) then
      c2(k) = min(rpos(k,iu),rneg(k,iu7))
   else
      c2(k) = min(rpos(k,iu7),rneg(k,iu))
   endif

   if (ahflx3(k,iu) >= 0) then
      c3(k) = min(rpos(k,iu),rneg(k,iu10))
   else
      c3(k) = min(rpos(k,iu10),rneg(k,iu))
   endif

   if (ahflx4(k,iu) >= 0) then
      c4(k) = min(rpos(k,iu),rneg(k,iu11))
   else
      c4(k) = min(rpos(k,iu11),rneg(k,iu))
   endif

! Update UM from long timestep tendencies, advection, and pressure gradient force

!-------SPECIAL--------------------
   if (k < kb) then
      c1(k) = .5 * c1(k)
      c2(k) = .5 * c2(k)
      c3(k) = .5 * c3(k)
      c4(k) = .5 * c4(k)
   endif
!----------------------------------

   umc(k,iu) = umc0(k,iu) &
      + dts * (pgf(k) + umt_rayf(k) + volui(k,iu) * (umt(k,iu) &

      + cb(k) * avflx (k-1,iu) &
      - ct(k) * avflx (k  ,iu) &
      + c1(k) * ahflx1(k  ,iu) &
      + c2(k) * ahflx2(k  ,iu) &
      + c3(k) * ahflx3(k  ,iu) &
      + c4(k) * ahflx4(k  ,iu)))

! UC diagnosis from UMC and RHO - Form from Wenneker      

      uc(k,iu) = umc(k,iu) / (volui(k,iu)  &
         * (volt(k,iw2) * rho(k,iw1) + volt(k,iw1) * rho(k,iw2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment for 10-meter mountain experiment only !!!!
! [MUST NOW ACCOUNT FOR LOW-ORDER + ANTIDIFFUSIVE FLUX]
!      uc2 = uc(k,iu) + uc(k+1,iu) + cnvel * (uc(k,iu) - uc(k+1,iu))
!      call uwcomp(1,iu,k,fadvz(k),wmarw2,uc2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

! Below-ground points

uc(1:ka-1,iu) = uc(ka,iu)

!!!!!!!!!!!! Special section (uncomment only for Held-Suarez experiment)
!hs   do k = lpu(iu),mza-1                                         !hs
!hs      pressloc = .5 * (press(k,iw1) + press(k,iw2))             !hs
!hs      rayf_cof(k) = max(0.,(pressloc-7.e4)/3.e4) / 86400.       !hs

!hs      umc(k,iu) = umc(k,iu) + dts * rayf_cof(k)  &              !hs
!hs         * .5 * (rho_s(k,iw1) + rho_s(k,iw2)) * (0. - uc(k,iu)) !hs
!hs   enddo                                                        !hs
!!!!!!!!!!!! End special section (uncomment only for Held-Suarez experiment)

return
end subroutine prog_u0


