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
subroutine timestep()

use misc_coms,  only: io6, time8, time_istp8, nqparm, initial, ilwrtyp,   &
                      iswrtyp, dtsm, nqparm_sh, dtlm, iparallel,   &
                      s1900_init, s1900_sim
use mem_ijtabs, only: nstp, istp, mrls, leafstep
use mem_nudge,  only: nudflag
use mem_grid,   only: mza, mva, mwa
use micro_coms, only: level
use leaf_coms,  only: isfcl
use mem_para,   only: myrank
use massflux,   only: zero_momsc, timeavg_momsc

use mem_basic  ! needed only when print statements below are uncommented
use mem_tend   ! needed only when print statements below are uncommented
use mem_leaf   ! needed only when print statements below are uncommented
use ed_options, only: frq_phenology

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_recv_v

use oplot_coms,  only: op

use mem_timeavg, only: accum_timeavg

implicit none

integer :: is,mvl,isl4,jstp,iw
real :: t1,w1

! automatic arrays

real :: vmsc(mza,mva) ! V face momentum for scalar advection
real :: wmsc(mza,mwa) ! W face momentum for scalar advection
real(kind=8) :: rho_old(mza,mwa) ! density at beginning of timestep [kg/m^3]

real :: alpha_press(mza,mwa) ! 
real :: rhot       (mza,mwa) ! grid-cell total mass tendency [kg/s]

! +----------------------------------------------------------------------------+
! |  Each call to subroutine timestep drives all steps in advancing by dtlong  |
! +----------------------------------------------------------------------------+

time_istp8 = time8

if (time8 < 1.e-3) then
!   call bubble()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wm2&5
   call diagn_global()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wm2&5

endif

do jstp = 1,nstp  ! nstp = no. of finest-grid-level aco steps in dtlm(1)
   istp = jstp

   call tend0(rhot)

! call check_nans(1)

   if (ilwrtyp + iswrtyp > 0) then
      call radiate()
   endif

! call check_nans(2)

   if (any(nqparm   (1:mrls) > 0) .or.  &
       any(nqparm_sh(1:mrls) > 0)) then

      call cuparm_driver(rhot)

      if (isfcl == 1) then
         call surface_cuparm_flux()
      endif
   endif

! call check_nans(6)

!nud   if (initial == 2 .and. nudflag == 1)  &
!nud   call obs_nudge(rhot)

! call check_nans(10)

!   call zero_momsc(vmsc,wmsc,rho_old)

! call check_nans(11)

   call prog_wrtv(vmsc,wmsc,alpha_press,rhot)

! call check_nans(12)

!   call timeavg_momsc(vmsc,wmsc)

! call check_nans(13)

!   call scalar_transport(vmsc,wmsc,rho_old)

! call check_nans(14)

!   call predtr(rho_old)

! call check_nans(15)

   if (iparallel == 1) then
      call mpi_recv_v('V')  ! Recv V group
   endif

! call check_nans(16)

   if (level /= 3) then
      call thermo()
   endif

! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
if (mod(real(time8),op%frqplt) < dtlm(1) .and. istp == nstp+1000) then

   allocate (op%extfld(mza,mwa))
   op%extfld(:,:) = thil(:,:)
   op%extfldname = 'THIL'
   call plot_fields(1)
   deallocate (op%extfld)

   allocate (op%extfld(mza,mwa))
   op%extfld(:,:) = theta(:,:)
   op%extfldname = 'THETA'
   call plot_fields(2)
   deallocate (op%extfld)

   allocate (op%extfld(mza,mwa))
   op%extfld(:,:) = (sh_w(:,:) - sh_v(:,:)) * 1.e3
   op%extfldname = 'SH_TOTCOND'
   call plot_fields(3)
   deallocate (op%extfld)

endif
! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

! call check_nans(17)

   if (level == 3) then
      call micro()  ! maybe later make freq. uniform

      if (isfcl == 1) then
!         call surface_precip_flux()
      endif
   endif

! call check_nans(18)

!   call trsets()  

! call check_nans(19)

   if (iparallel == 1) then
      call mpi_send_w('T')  ! Send W group
      call mpi_recv_w('T')  ! Recv W group
   endif

! call check_nans(20)

   if (iparallel == 1) then
      call mpi_send_w('S')  ! Send scalars
      call mpi_recv_w('S')  ! Recv scalars
   endif

! call check_nans(21)

   if (leafstep(istp) > 0) then
      call leaf3()

      if (iparallel == 1) call mpi_send_wl('T')
      if (iparallel == 1) call mpi_recv_wl('T')

      call seacells()

      if (iparallel == 1) call mpi_send_ws('T')
      if (iparallel == 1) call mpi_recv_ws('T')
   endif

! call check_nans(22)

   time_istp8 = time8 + float(istp) * dtsm(mrls)  ! Update precise time
   s1900_sim = s1900_init + time_istp8

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wm2&5
   call diagn_global()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wm2&5

! Call ED model if it is time to do vegetation dynamics

if (mod(real(time8)+dtlm(1),frq_phenology) < dtlm(1)) call ed_vegetation_dynamics()  

return
end subroutine timestep

!==========================================================================

!!subroutine check_nans(icall)
!!
!!use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,vmc,vc
!!use mem_micro,  only: sh_c,sh_r
!!use mem_grid,   only: mza,mwa,lpw,volti,mva
!!use mem_tend,   only: thilt, vmt
!!use micro_coms, only: level
!!use misc_coms,  only: io6, iparallel
!!use mem_ijtabs, only: itab_v, itab_w
!!
!!implicit none
!!  
!!integer :: i,k,iv
!!integer, intent(in) :: icall
!!
!!!do iv = 1,mva
!!!   if (itab_v(iv)%ivglobe == 1206) then
!!!      write(io6,'(a,i6,5e15.7)') 'cns ',icall,vmc(2,iv),vc(2,iv),vmt(2,iv)
!!!   endif
!!!enddo
!!
!!return
!!
!!do i = 2,mwa
!!   do k = lpw(i),mza
!!      if (isnan(sh_w (k,i)) .or.  &
!!          isnan(rho  (k,i)) .or.  &
!!          isnan(thil (k,i)) .or.  &
!!          isnan(thilt(k,i)) .or.  &
!!          isnan(press(k,i)) .or.  &
!!          isnan(wc   (k,i)) .or.  &
!!          isnan(wmc  (k,i)) .or.  &
!!          isnan(sh_v (k,i)) .or.  &
!!          thil(k,i) < 100.0) then
!!
!!         write(io6,*) ''
!!         write(io6,*) 'check_nans',k,i,icall
!!         write(io6,*) 'sh_w,rho,thil',sh_w(k,i),rho(k,i),thil(k,i)
!!         write(io6,*) 'thilt,press',thilt(k,i),press(k,i)
!!         write(io6,*) 'wc, wmc, sh_v',wc(k,i),wmc(k,i),sh_v(k,i)
!!         stop
!!      endif
!!
!!!      if (level >= 3) then
!!!         if (isnan(sh_c(k,i)) .or.  &
!!!             isnan(sh_r(k,i))) then
!!               
!!!            write(io6,*) 'nan',k,i,icall
!!!            write(io6,*) 'sh_c,sh_r',sh_c(k,i),sh_r(k,i)
!!!            stop
!!!         endif
!!!      endif
!!
!!   enddo
!!enddo
!!
!!return
!!end subroutine check_nans

!===============================================================================

subroutine modsched()

use mem_ijtabs, only: nstp, mrls,  &
                      mrl_endl, mrl_ends, mrl_begl, mrl_begs, leafstep

use misc_coms,  only: io6, nacoust, ndtrat, dtlm, dtlong, dtsm, nqparm
use leaf_coms,  only: dt_leaf, mrl_leaf, isfcl
use sea_coms,   only: dt_sea

implicit none

integer :: ndts
integer :: ng
integer :: j
integer :: i
integer :: ii
integer :: iw
integer :: jstp
integer :: mrl

integer :: nshort(mrls)  ! automatic array
integer :: nlong(mrls)   ! automatic array

! Find number of short timesteps of most refined mesh done per long and short
!    timestep of every refinement level.

nshort(mrls) = 1
nlong(mrls) = nacoust(mrls)

do mrl = mrls-1,1,-1
   nlong(mrl) = nlong(mrl+1) * ndtrat(mrl+1)
   nshort(mrl) = nlong(mrl) / nacoust(mrl)
   ndts = nshort(mrl) / nshort(mrl+1)
   if (nshort(mrl) /= ndts * nshort(mrl+1)) then
      write(io6,'(/,a)')   'Short timesteps not integer ratio between consec MRLs.'
      write(io6,'(a,2i5)') 'mrl, nshort(mrl) = ',mrl,nshort(mrl)
      write(io6,'(a,2i5)') 'mrl+1, nshort(mrl+1) = ',mrl+1,nshort(mrl+1)
      stop 'stop modsched_1'
   endif
enddo

nstp = nlong(1)

do mrl = 1,mrls
   write(io6,*) 'modsched-0 ',mrl,nlong(mrl),nshort(mrl)
enddo

write(io6,'(/,a)') '=== Timestep Schedule ===='
write(io6,'(a,/)') '              jstp    mrl_l  mrl_s'

! Long timestep for each mesh refinement level

dtlm(1) = dtlong
dtsm(1) = dtlm(1) / nacoust(1)

do mrl = 2,mrls
   dtlm(mrl) = dtlm(mrl-1) / float(ndtrat(mrl))
   dtsm(mrl) = dtlm(mrl) / nacoust(mrl)
enddo

! Default:  Fill dt_leaf, mrl_leaf with values for MRL = 1

dt_leaf = dtlm(1)
mrl_leaf = 1
   
! Allocate mrl-schedule arrays

allocate (mrl_begl(nstp))
allocate (mrl_begs(nstp))
allocate (mrl_endl(nstp))
allocate (mrl_ends(nstp))

allocate (leafstep(nstp))

leafstep(1:nstp) = 0
   
! Fill mrl-table values

do jstp = 1,nstp

   mrl_begl(jstp) = 0
   mrl_begs(jstp) = 0
   mrl_endl(jstp) = 0
   mrl_ends(jstp) = 0

! Fill lowest MRL value for which loop processes are to be carried out when
! jstp has its current value.  (However, a zero MRL value indicates process
! is inactive on current jstp.)

   do mrl = mrls,1,-1
      if (mod(jstp-1, nlong(mrl)) == 0) mrl_begl(jstp) = mrl
      if (mod(jstp-1,nshort(mrl)) == 0) mrl_begs(jstp) = mrl

      if (mod(jstp  , nlong(mrl)) == 0) mrl_endl(jstp) = mrl
      if (mod(jstp  ,nshort(mrl)) == 0) mrl_ends(jstp) = mrl
   enddo

   write(io6,333) jstp,mrl_begl(jstp),mrl_begs(jstp),mrl_endl(jstp),mrl_ends(jstp)
   333 format('modsched0 ',5i7)

   if (mrl_ends(jstp) == 0 .or. mrl_begs(jstp) == 0) then
      write(io6,*) 'mrl_s value = 0 is not allowed'
      stop 'stop - modsched_2'
   endif

! Set LEAFSTEP = 1 to do leaf timestep for selected jstp value(s):
!    1. Always do leaf at end of long timestep for mrl = 1
!    2. Also do leaf at end of long timestep for any mrl > 1 whose dtlm > 30 s

   mrl = mrl_endl(jstp)
   if (isfcl == 1 .and. mrl > 0) then
      if (mrl == 1 .or. dtlm(mrl) > 30.) then
         leafstep(jstp) = 1

! Set leaf mrl and timestep according to highest selected mrl

         mrl_leaf = max(mrl_leaf,mrl)
         dt_leaf = min(dt_leaf,dtlm(mrl))
      endif
   endif

enddo

! Set dt_sea equal to dt_leaf

dt_sea = dt_leaf

return
end subroutine modsched

!==========================================================================

subroutine tend0(rhot)

use mem_ijtabs, only: jtab_w, jtab_v, istp, mrl_begl
use var_tables, only: scalar_tab, num_scalar
use mem_grid,   only: mza, mwa, mva, lcv, lpw
use mem_tend,   only: wmt, vmt, thilt
use misc_coms,  only: io6

!$ use omp_lib

implicit none

real, intent(in) :: rhot

integer :: n,mrl,j,k,iw,iv

! SET SCALAR TENDENCIES TO ZERO

mrl = mrl_begl(istp)

if (mrl > 0) then
   do n = 1,num_scalar
      call tnd0(scalar_tab(n)%var_t)
   enddo
   call tnd0(thilt)
   call tnd0(rhot)
endif

! SET W MOMENTUM TENDENCY TO ZERO

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(14)%jend(mrl); iw = jtab_w(14)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
      wmt(k,iw) = 0.
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',14)

! SET V MOMENTUM TENDENCY TO ZERO

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iv,k)
do j = 1,jtab_v(11)%jend(mrl); iv = jtab_v(11)%iv(j)
!----------------------------------------------------------------------
call qsub('V',iv)
   do k = lcv(iv),mza-1
      vmt(k,iv) = 0.
   enddo
enddo
!$omp end parallel do
endif
call rsub('V',11)

return
end subroutine tend0

!==========================================================================

subroutine tnd0(vart)

use mem_ijtabs, only: jtab_w, istp, mrl_begl
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: io6
!$ use omp_lib

implicit none

real, intent(out) :: vart(mza,mwa)

integer :: j,iw,k,mrl

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(11)%jend(mrl); iw = jtab_w(11)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw)-1,mza
      vart(k,iw) = 0.
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',11)

return
end subroutine tnd0

!==========================================================================

subroutine predtr(rho_old)

use var_tables, only: num_scalar, scalar_tab
use mem_ijtabs, only: istp, jtab_w, mrl_endl
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6

implicit none

real(kind=8), intent(in) :: rho_old(mza,mwa)

integer :: n,mrl

!   -  Step thermodynamic variables from  t  to  t+1.
!   -  Set top, lateral and bottom boundary conditions on some variables
!        if needed.
!   -  Call adjustment to assure all positive definite quantities
!        remain positive.
!   -  Rediagnose some thermodynamic quantities for use on the small
!        timestep.

!     Update the scalars and apply lateral, top, and bottom boundary
!     conditions.

mrl = mrl_endl(istp)
if (mrl > 0) then
do n = 1,num_scalar
   call o_update(n,scalar_tab(n)%var_p,scalar_tab(n)%var_t,rho_old)
enddo
endif

return
end subroutine predtr

!==========================================================================

subroutine o_update(n,varp,vart,rho_old)

use mem_ijtabs, only: jtab_w, istp, itab_w, mrl_endl
use mem_basic,  only: rho
use misc_coms,  only: io6, dtlm
use mem_grid,   only: mza, mwa, lpw
!$ use omp_lib

implicit none

integer, intent(in) :: n

real, intent(out) :: varp(mza,mwa)
real, intent(in) :: vart(mza,mwa)

real(kind=8), intent(in) :: rho_old(mza,mwa)

integer :: iw,j,k,mrl
real :: dtl

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k,dtl)
do j = 1,jtab_w(27)%jend(mrl); iw = jtab_w(27)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)
   dtl = dtlm(itab_w(iw)%mrlw)
   do k = lpw(iw),mza-1
      varp(k,iw) = (varp(k,iw) * rho_old(k,iw) + dtl * vart(k,iw)) / rho(k,iw) 
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',27)!RRR

return
end subroutine o_update

!==========================================================================

subroutine bubble()

use mem_basic, only: thil, theta
use misc_coms, only: io6
!use mem_grid
!use mem_ijtabs

implicit none

integer :: iw,i,j,k


   do k = 2,2
      thil(k,17328) = thil(k,17328) + 5. 
      theta(k,17328) = theta(k,17328) + 5. 

      thil(k,17329) = thil(k,17329) + 5. 
      theta(k,17329) = theta(k,17329) + 5. 

      thil(k,17333) = thil(k,17333) + 5. 
      theta(k,17333) = theta(k,17333) + 5. 

      thil(k,17334) = thil(k,17334) + 5. 
      theta(k,17334) = theta(k,17334) + 5. 

      thil(k,17335) = thil(k,17335) + 5. 
      theta(k,17335) = theta(k,17335) + 5. 

      thil(k,17336) = thil(k,17336) + 5. 
      theta(k,17336) = theta(k,17336) + 5. 
   enddo


return
end subroutine bubble

!==========================================================================

subroutine diagn_global()

! subroutine for diagnosis of williamson et al expt 1
     
use mem_basic
use mem_micro
use micro_coms
use mem_ijtabs
use misc_coms
use consts_coms
use mem_grid
use oplot_coms,  only: op

!use mem_wm5_refsoln
use mem_wm5_refsoln_cubic

! The following NSLCON is used as a namelist flag to indicate shallow water 
! test case 1, 2, or 5

use leaf_coms, only: nslcon

implicit none

integer, save :: ncall = 0
integer :: j,iw,kk,iday
integer :: iu,iv,iw1,iw2
real :: raxis,u,v,height,vol,volu,volv

real :: alpha_wm1

real :: sum_abshdif,sum_abshtr,sum_hdif2 ,sum_htr2 ,abshdif_max,abshtr_max
real :: sum_vcdif  ,sum_vctr  ,sum_vcdif2,sum_vctr2,vcdif_max  ,vctr_max

real, save, dimension(18000) :: ge1,ge2,ge3,ve1,ve2,ve3,vn1,vn2,vn3,vctr18
real, save :: height_init(141000),uc_init(141000),vc_init(141000)

real, save :: aspect = .7
real, save :: scalelab = .014
real, save :: timebeg,timeend,timedif,timeinc

real :: uv01dz,uv01dy,uv01dx,uv01dr
real :: uc_tr, uc0_tr, vc_tr, vc0_tr

real :: zanal0_wm5,zanal_wm5

! On the first call to this subroutine, compute the plotting time increment

if (ncall == 0) then

   timebeg = time8   / 86400.
   timeend = timmax8 / 86400.
   timedif = timeend - timebeg

   if (timedif < .03) then
      timeinc = .001
   elseif (timedif < .06) then
      timeinc = .002
   elseif (timedif < .1) then
      timeinc = .004

   elseif (timedif < .3) then
      timeinc = .01
   elseif (timedif < .6) then
      timeinc = .02
   elseif (timedif < 1.) then
      timeinc = .04

   elseif (timedif < 3.) then
      timeinc = .1
   elseif (timedif < 6.) then
      timeinc = .2
   elseif (timedif < 10.) then
      timeinc = .4

   elseif (timedif < 30.) then
      timeinc = 1.
   elseif (timedif < 60.) then
      timeinc = 2.
   elseif (timedif < 100.) then
      timeinc = 4.

   elseif (timedif < 300.) then
      timeinc = 10.
   elseif (timedif < 600.) then
      timeinc = 20.
   elseif (timedif < 1000.) then
      timeinc = 40.
   endif
      
endif

ncall = ncall + 1

alpha_wm1 = 0.  ! Use four settings:  0., .05, pio2 - .05, pio2

! Initialize summation and max/min quantities to zero

sum_abshdif = 0.
sum_abshtr  = 0.

sum_hdif2   = 0.
sum_htr2    = 0.

abshdif_max = 0.
abshtr_max  = 0.

sum_vcdif  = 0.
sum_vctr   = 0.

sum_vcdif2 = 0.
sum_vctr2  = 0.

vcdif_max = 0.
vctr_max  = 0.

!------------------------------------------------------------

! Horizontal loop over all IW points

!----------------------------------------------------------------------
do j = 1,jtab_w(8)%jend(1); iw = jtab_w(8)%iw(j)
!---------------------------------------------------------------------

if (nslcon > 2) then

! Due to problem with reference solution (or its interpolation), 
! omit high latitude points

if (glatw(iw) > 85. .or. glatw(iw) < -85.) cycle

! Get WM5 reference solution (for days 0 and 15) for this point

!   call extract_wm5(glonw(iw),glatw(iw),zanal0_wm5,zanal_wm5)

   call npr_bicubics(zanal00_wm5,glatw(iw),glonw(iw),zanal0_wm5)
   call npr_bicubics(zanal15_wm5,glatw(iw),glonw(iw),zanal_wm5)

else

! Case for WM2

   zanal0_wm5 = 0.
   zanal_wm5 = 0.
   
endif

! Get height and volume for current IW point

   height = rho(2,iw)
   vol = volt(2,iw)

! first time in this routine, save initial values: uh(iw) = u; vh(iw) = v

   if (ncall == 1) then
      height_init(iw) = height
   endif

! Find sums and max values for height

   sum_abshdif = sum_abshdif + vol *  &
                 abs(height - height_init(iw) - zanal_wm5 + zanal0_wm5)

   sum_hdif2 = sum_hdif2 + vol *  &
               (height - height_init(iw) - zanal_wm5 + zanal0_wm5)**2

   abshdif_max = max(abshdif_max,  &
                 abs(height - height_init(iw) - zanal_wm5 + zanal0_wm5))

if (nslcon <= 2) then
   sum_abshtr = sum_abshtr + vol * abs(height_init(iw))
   sum_htr2   = sum_htr2   + vol * height_init(iw)**2
   abshtr_max = max(abshtr_max,    height_init(iw))
else
   sum_abshtr = sum_abshtr + vol * abs(zanal_wm5)
   sum_htr2   = sum_htr2   + vol * zanal_wm5**2
   abshtr_max = max(abshtr_max,    zanal_wm5)
endif

enddo

! First 3 normalized global errors (Williamson et al 1992 Eqs. 82-84)

ge1(ncall) = sum_abshdif / sum_abshtr
ge2(ncall) = sqrt(sum_hdif2) / sqrt(sum_htr2)
ge3(ncall) = abshdif_max / abshtr_max

! 3 normalized global normal-velocity errors

if (nslcon <= 2) then
   vn1(ncall) = sum_vcdif / sum_vctr
   vn2(ncall) = sum_vcdif2 / sum_vctr2
   vn3(ncall) = vcdif_max / vctr_max
endif

!write(io6,'(a,i6,9e10.3)') 'vns ',ncall,vn1(ncall),vn2(ncall),vn3(ncall),  &
!            sum_ucdif,sum_uctr,sum_ucdif2,sum_uctr2,ucdif_max,uctr_max


vctr18(ncall) = time8 / 86400.

if (nslcon <= 2 .and. time8 + .5 * dtlong <  432000.        ) return
if (nslcon >  2 .and. time8 + .5 * dtlong < 1296000. - 7200.) return

! Reopen the current graphics output workstation if it is closed

call o_reopnwk()

!----------------------------------------------------------------------
! Plot 3 height norms on single frame
!----------------------------------------------------------------------
call plotback()
call oplot_xy2log10('0','N',aspect,scalelab &
              ,ncall,  vctr18,ge1          &
              ,'TIME(DAYS)',' '             &
              ,timebeg,timeend,timeinc,5  ,-6,-1  )
call oplot_xy2log10('0','N',aspect,scalelab &
              ,ncall,  vctr18,ge2          &
              ,' ',' '                      &
              ,timebeg,timeend,timeinc,5  ,-6,-1  )
call oplot_xy2log10('0','N',aspect,scalelab &
              ,ncall,  vctr18,ge3          &
              ,' ',' '                      &
              ,timebeg,timeend,timeinc,5  ,-6,-1  )
call o_frame()

!----------------------------------------------------------------------
! Plot 3 velocity norms on single frame [using normal wind component]
!----------------------------------------------------------------------

if (nslcon <= 2) then

   call plotback()
   call oplot_xy2log10('0','N',aspect,scalelab &
                 ,ncall,  vctr18,vn1          &
                 ,'TIME(DAYS)',' '             &
                 ,timebeg,timeend,timeinc,5  ,-7,-0  )
   call oplot_xy2log10('0','N',aspect,scalelab &
                 ,ncall,  vctr18,vn2          &
                 ,' ',' '                      &
                 ,timebeg,timeend,timeinc,5  ,-7,-0  )
   call oplot_xy2log10('0','N',aspect,scalelab &
                 ,ncall,  vctr18,vn3          &
                 ,' ',' '                      &
                 ,timebeg,timeend,timeinc,5  ,-7,-0  )
   call o_frame()

endif

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

if ((trim(runtype) /= 'PLOTONLY') .and. (op%plttype == 0)) then
   call o_clswk()
endif

return
end subroutine diagn_global
