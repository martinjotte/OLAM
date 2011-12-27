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
use massflux,   only: zero_momsc, timeavg_momsc, diagnose_uc

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

   if (initial == 2 .and. nudflag == 1)  &
   call obs_nudge(rhot)

! call check_nans(10)

   call zero_momsc(vmsc,wmsc,rho_old)

! call check_nans(11)

   call prog_wrtv(vmsc,wmsc,alpha_press,rhot)

! call check_nans(12)

   call timeavg_momsc(vmsc,wmsc)

! call check_nans(13)

   call scalar_transport(vmsc,wmsc,rho_old)

! call check_nans(14)

   call predtr(rho_old)

! call check_nans(15)

   if (iparallel == 1) then
      call mpi_recv_v('V')  ! Recv V group
   endif

   call diagnose_uc()

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
         call surface_precip_flux()
      endif
   endif

! call check_nans(18)

   call trsets()  

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

   call accum_timeavg()

! call check_nans(22)

   time_istp8 = time8 + float(istp) * dtsm(mrls)  ! Update precise time
   s1900_sim = s1900_init + time_istp8

enddo

! Call ED model if it is time to do vegetation dynamics

if (mod(real(time8)+dtlm(1),frq_phenology) < dtlm(1)) call ed_vegetation_dynamics()  

return
end subroutine timestep

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

real, intent(inout) :: vart(mza,mwa)

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
