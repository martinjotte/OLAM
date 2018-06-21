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
subroutine timestep()

use misc_coms,   only: io6, time8, time_istp8, time_istp8p, time_bias, &
                       nqparm, initial, ilwrtyp, iswrtyp, dtsm, dtlm, &
                       iparallel, isubdomain, s1900_init, s1900_sim, do_chem
use mem_ijtabs,  only: nstp, istp, mrls, leafstep, mrl_begl, mrl_endl, mrl_ends
use mem_nudge,   only: nudflag, nudnxp, o3nudflag, io3
use mem_grid,    only: mza, mva, mwa
use micro_coms,  only: miclevel
use leaf_coms,   only: isfcl
use mem_basic,   only: vmc, vc, vxe, vye, vze, thil, rho, wmc, wc
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
use var_tables,  only: nvar_par, vtab_r, nptonv
use obnd,        only: trsets, lbcopy_v, lbcopy_w, botset, latsett
use oplot_coms,  only: op
use oname_coms,  only: nl
use mem_flux_accum, only: flux_accum
use consts_coms, only: r8
use vel_t3d,     only: diagvel_t3d
use var_tables,  only: num_scalar, scalar_tab
use mem_megan,   only: megan_avg_temp
use emis_defn,   only: get_emis
use depv_defn,   only: get_depv
use wrtv_mem,    only: prog_wrtv
use mem_turb,    only: vkm, vkh
use check_nan,   only: check_nans, compute_mass_sums

implicit none

integer :: jstp, mrl, n

real     :: vmsc       (mza,mva) ! V face momentum for scalar advection
real     :: wmsc       (mza,mwa) ! W face momentum for scalar advection
real(r8) :: rho_old    (mza,mwa) ! density at beginning of long timestep [kg/m^3]
real     :: alpha_press(mza,mwa) ! coefficient for computing pressure
real     :: rhot       (mza,mwa) ! grid-cell total mass tendency [kg/s]
real     :: vxesc      (mza,mwa)
real     :: vyesc      (mza,mwa)
real     :: vzesc      (mza,mwa)

real(r8) :: time0

! +----------------------------------------------------------------------------+
! |  Each call to subroutine timestep drives all steps in advancing by dtlong  |
! +----------------------------------------------------------------------------+

if (time_istp8 < 1.e-3_r8) then
   ! call bubble()

   ! For shallow water test cases, compute error norms at initial time

   if (nl%test_case == 2 .or. nl%test_case == 5) then
      call diagn_global_swtc()
   endif

   ! For ncar dcmip test cases, compute error norms at initial time

   if (nl%test_case >= 10 .and. nl%test_case < 900) then
      call diagn_global_dcmip()
   endif
endif

do jstp = 1,nstp  ! nstp = no. of finest-grid-level aco steps in dtlm(1)
   istp = jstp

   ! Bypass all processes except microphyics if running microphysics parcel simulation

   if (nl%test_case >= 901 .and. nl%test_case <= 950) go to 1311

   call tend0(rhot)

   mrl = mrl_begl(istp)
   if (mrl > 0) then
      call comp_alpha_press(mrl, alpha_press)
      call surface_turb_flux(mrl)
      call sea_spray(mrl)
      call dust_src(mrl)
   endif

   ! call check_nans(1,rvara1=rhot,rvara2=alpha_press)

   if (any( nqparm(1:mrls) > 0 )) then
      call cuparm_driver(rhot)
      if (isfcl == 1) then
         call surface_cuparm_flux()
      endif
   endif

   ! call check_nans(2,rvara1=rhot,rvara2=alpha_press)

   if (ilwrtyp + iswrtyp > 0) then
      call radiate()
   endif

   ! small-scale divergence/vorticity damping

   mrl = mrl_begl(istp)
   if (mrl > 0) then
      call vort_damp(mrl)
      call divh_damp(mrl)
   endif

   ! Compute nudging tendencies

   if (initial == 2 .and. nudflag == 1) then
      if (nudnxp == 0) then
         call  obs_nudge(rhot)
      else
         call spec_nudge(rhot)
      endif
   endif

   ! Nudging of ozone if it is in the scalar table

   if (initial == 2 .and. o3nudflag == 1 .and. io3 > 0) then
      call obs_nudge_o3()
   endif

   ! Long timestep PBL tendencies (computation of K's, scalar diffusion, 
   ! CMAQ emissions and deposition, and lateral friction with shaved cells)

   mrl = mrl_begl(istp)
   if (mrl > 0) then

      if (do_chem == 1) then
         call get_emis ( mrl )
         call get_depv ( mrl )
      endif

      call pbl_driver(mrl,rhot)

      ! MPI send of computed K's
      if (iparallel == 1) then
         call mpi_send_w(mrl, rvara1=vkm, rvara2=vkh)
      endif

      if (isfcl == 1) call lateral_friction(mrl)

      ! MPI recv of computed K's
      if (iparallel == 1) then
         call mpi_recv_w(mrl, rvara1=vkm, rvara2=vkh)
      endif
      call lbcopy_w(mrl, a1=vkm, a2=vkh)

      call comp_horiz_k(mrl)

   endif

   ! call check_nans(5,rvara1=rhot,rvara2=alpha_press)

   call zero_momsc(vmsc,wmsc,vxesc,vyesc,vzesc,rho_old)

   ! call check_nans(11,rvara1=rhot,rvara2=alpha_press)

   ! Call olam_dcmip_phys, which is the OLAM interface to DCMIP auxiliary
   ! physics subroutine that provides tendencies to some model fields

   mrl = mrl_begl(istp)
   if (mrl > 0) then
      !--------------------------------------
      if (nl%test_case ==  42 .or. &
          nl%test_case ==  43 .or. &
          nl%test_case == 110 .or. &
          nl%test_case == 111 .or. &
          nl%test_case == 112 .or. &
          nl%test_case == 113 .or. &
          nl%test_case == 114 .or. &
          nl%test_case == 121 .or. &
          nl%test_case == 122 .or. &
          nl%test_case == 131) then
      !--------------------------------------
         call olam_dcmip2016_phys(rhot)
      endif
   endif

   ! Call olam_dcmip_terminator, which is the OLAM interface to DCMIP auxiliary
   ! chemistry subroutine terminator that provides chemical tendencies.

   mrl = mrl_begl(istp)
   if (mrl > 0) then
      !--------------------------------------
      if (nl%test_case == 110 .or. &
          nl%test_case == 111 .or. &
          nl%test_case == 112 .or. &
          nl%test_case == 113 .or. &
          nl%test_case == 114) then
      !--------------------------------------
         call olam_dcmip_terminator()
      endif
   endif

   ! call check_nans(11,rvara1=rhot,rvara2=alpha_press)

   !    write(*,'(a)') ' calling mass_sums2 '
   !    call compute_mass_sums()

   ! Bypass call to thiltend_long if using prescribed flow for DCMIP tests
   !--------------------------------------
   if (nl%test_case == 11 .or. &
       nl%test_case == 12 .or. &
       nl%test_case == 13) go to 33
   !--------------------------------------

   call prog_wrtv(vmsc,wmsc,vxesc,vyesc,vzesc,alpha_press,rhot)

   33 continue

   !    write(*,'(a)') ' calling mass_sums3 '
   !    call compute_mass_sums()

   !--------------------------------------
   ! Get flow update if using DCMIP prescribed flow

   if (nl%test_case == 11 .or. &
       nl%test_case == 12) then

      ! Set time0 to half-forward time (of grid-1 short timestep)
      ! for time-centered advective tendencies

      time0 = time8 + .5 * dtsm(1)

      call olam_dcmip_prescribedflow(time0,vmsc,wmsc,rho_old,vxesc,vyesc,vzesc)

   endif
   !--------------------------------------

   ! call check_nans(12,rvara1=rhot,rvara2=alpha_press)

   mrl = mrl_endl(istp)
   if (nl%split_scalars > 0 .and. mrl > 0) then
      call predtr_split(mrl,rho_old)
      do n = 1, num_scalar
         call latsett(mrl,scalar_tab(n)%var_p)
         call botset (mrl,scalar_tab(n)%var_p)
      enddo
      if (iparallel == 1) call mpi_send_w(mrl, scalars='S')
      ! call check_pos(1)
   endif

   ! Bypass call to timeavg_momsc if using prescribed flow for DCMIP tests
   !--------------------------------------
   if (nl%test_case == 11 .or. &
       nl%test_case == 12 .or. &
       nl%test_case == 13) go to 34
   !--------------------------------------

   call timeavg_momsc(vmsc,wmsc,vxesc,vyesc,vzesc)

   34 continue

   ! call check_nans(13,rvara1=rhot,rvara2=alpha_press)

   mrl = mrl_endl(istp)
   if (nl%split_scalars > 0 .and. mrl > 0) then
      if (iparallel == 1) call mpi_recv_w(mrl, scalars='S')
   endif

   mrl = mrl_endl(istp)
   if (mrl > 0) then
      call scalar_transport(mrl,vmsc,wmsc,vxesc,vyesc,vzesc,rho_old)
      if (nl%split_scalars > 0) call scalar_hdiff_split(mrl, rho_old)
   endif

   ! call check_nans(14,rvara1=rhot,rvara2=alpha_press)

   call predtr(rho_old)

   ! call check_nans(15,rvara1=rhot,rvara2=alpha_press)

   ! Special diagnosis of water vapor for DCMIP tests; thil is dry theta in that case
   !---------------------------------------------------------------------------------
   if (nl%test_case == 110 .or. &
       nl%test_case == 111 .or. &
       nl%test_case == 112 .or. &
       nl%test_case == 113 .or. &
       nl%test_case == 114 .or. &
       nl%test_case == 121 .or. &
       nl%test_case == 122 .or. &
       nl%test_case == 131) then

      call thermo_dcmip()

      go to 35
   endif
   !--------------------------------------

   if (miclevel /= 3) then
      call thermo()
   endif

   35 continue

   ! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
   !if (mod(time8,op%frqplt) < dtlm(1) .and. istp == nstp+1000) then
   !
   !   allocate (op%extfld(mza,mwa))
   !   op%extfld(:,:) = thil(:,:)
   !   op%extfldname = 'THIL'
   !   call plot_fields(1)
   !   deallocate (op%extfld)
   !
   !   allocate (op%extfld(mza,mwa))
   !   op%extfld(:,:) = theta(:,:)
   !   op%extfldname = 'THETA'
   !   call plot_fields(2)
   !   deallocate (op%extfld)
   !
   !   allocate (op%extfld(mza,mwa))
   !   op%extfld(:,:) = (sh_w(:,:) - sh_v(:,:)) * 1.e3
   !   op%extfldname = 'SH_TOTCOND'
   !   call plot_fields(3)
   !   deallocate (op%extfld)
   !
   !endif
   ! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

   ! call check_nans(16,rvara1=rhot,rvara2=alpha_press)

   1311 continue

   mrl = mrl_endl(istp)
   if (miclevel == 3 .and. mrl > 0) then
      call micro()  ! maybe later make freq. uniform
   endif

   ! Bypass all processes except microphyics if running microphysics parcel simulation

   if (nl%test_case >= 901 .and. nl%test_case <= 950) go to 1312

   ! call check_nans(17,rvara1=rhot,rvara2=alpha_press)

   ! Call atmospheric chemistry here

   if (do_chem == 1) call cmaq_driver()

   ! Cyclic lateral boundaries and bottom boundary for scalars

   call trsets(mrl)

   ! call check_nans(18,rvara1=rhot,rvara2=alpha_press)

   mrl = mrl_ends(istp)

   if (iparallel == 1) then
      call mpi_send_w(mrl, dvara1=rho, rvara1=wmc, &
                           rvara2=wc,  rvara3=thil)

      call mpi_recv_w(mrl, dvara1=rho, rvara1=wmc, &
                           rvara2=wc,  rvara3=thil)
   endif

   if (mrl_ends(istp) > 0) then
      call lbcopy_w(mrl, a1=thil, a2=wmc, a3=wc, d1=rho)
   else
      call lbcopy_w(mrl, a1=thil, a2=wmc, a3=wc)
   endif

   ! call check_nans(19,rvara1=rhot,rvara2=alpha_press)

   mrl = mrl_endl(istp)
   if (mrl > 0) then

      if (iparallel == 1) call mpi_send_w(mrl, scalars='S')  ! Send scalars

      if (miclevel == 3) call omic_update_v_mom(mrl)

      if (iparallel == 1) call mpi_recv_w(mrl, scalars='S')  ! Recv scalars

      do n = 1, nvar_par
         call lbcopy_w(mrl, a1=vtab_r(nptonv(n))%rvar2_p)
      enddo

      ! call check_pos(3)
   endif

   ! call check_nans(20,rvara1=rhot,rvara2=alpha_press)

   ! MPI send/recv and lbc copy of VMC, VC

   mrl = mrl_ends(istp)

   if (iparallel == 1) then
      call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
      call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
   endif
   call lbcopy_v(1, vmc=vmc, vc=vc)

   ! Compute earth cartesian velocities

   if (mrl > 0) then
      call diagvel_t3d(mrl)
   endif

   ! call check_nans(21,rvara1=rhot,rvara2=alpha_press)

   if (leafstep(istp) > 0) then
      call leaf4()
      call seacells()
      if (do_chem == 1) call megan_avg_temp()
   endif

   ! call check_nans(22,rvara1=rhot,rvara2=alpha_press)

   1312 continue

   time_istp8  = time8 + istp * dtsm(mrls)  ! Update precise time
   time_istp8p = time_istp8 + time_bias
   s1900_sim   = s1900_init + time_istp8

   ! Add current fluxes to time integrals

   call flux_accum()

   ! if (jstp == nstp) then
   !    call compute_mass_sums()
   ! endif

enddo

! For ncar dcmip test cases, compute error norms 

if (nl%test_case >= 10 .and. nl%test_case < 900) then
   call diagn_global_dcmip()
endif

end subroutine timestep

!===============================================================================

subroutine modsched()

use mem_ijtabs, only: nstp, mrls,  &
                      mrl_endl, mrl_ends, mrl_begl, mrl_begs, leafstep

use misc_coms,  only: io6, nacoust, ndtrat, dtlm, dtlong, dtsm
use leaf_coms,  only: dt_leaf, mrl_leaf, isfcl
use sea_coms,   only: dt_sea
use consts_coms,only: r8

implicit none

integer :: ndts
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
   dtlm(mrl) = dtlm(mrl-1) / ndtrat(mrl)
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
      if (mrl == 1 .or. dtlm(mrl) > 30.0_r8) then
         leafstep(jstp) = 1

! Set leaf mrl and timestep according to highest selected mrl

         mrl_leaf = max(mrl_leaf, mrl)
         dt_leaf  = min(dt_leaf, real(dtlm(mrl)))
      endif
   endif

enddo

! Set dt_sea equal to dt_leaf

dt_sea = dt_leaf

end subroutine modsched

!==========================================================================

subroutine tend0(rhot)

use mem_ijtabs, only: jtab_w, istp, mrl_begl, jtw_wstn, jtab_v, jtv_prog
use var_tables, only: scalar_tab, num_scalar
use mem_grid,   only: mza, mwa, lpw, lpv
use mem_tend,   only: thilt, vmxet, vmyet, vmzet, vmt

implicit none

real, intent(inout) :: rhot(mza,mwa)

integer :: n,mrl,j,k,iw,iv

mrl = mrl_begl(istp)
if (mrl > 0) then

   ! SET TENDENCIES AT W POINTS TO ZERO

   !$omp parallel
   !$omp do private(iw,k,n)
   do j = 1,jtab_w(jtw_wstn)%jend(mrl); iw = jtab_w(jtw_wstn)%iw(j)

      do k = lpw(iw),mza
         vmxet(k,iw) = 0.0
         vmyet(k,iw) = 0.0
         vmzet(k,iw) = 0.0
         thilt(k,iw) = 0.0
         rhot (k,iw) = 0.0
      enddo

      do n = 1,num_scalar
         do k = lpw(iw),mza
            scalar_tab(n)%var_t(k,iw) = 0.0
         enddo
      enddo
   enddo
   !$omp end do

   ! SET TENDENCIES AT V PONTS TO ZERO

   !$omp do private(iv,k)
   do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

      do k = lpv(iv),mza
         vmt(k,iv) = 0.0
      enddo
   enddo
   !$omp end do
   !$omp end parallel

endif

end subroutine tend0

!==========================================================================

subroutine predtr(rho_old)

use var_tables, only: num_scalar, scalar_tab
use mem_ijtabs, only: istp, mrl_endl, jtab_w, itab_w, jtw_prog
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: dtlm
use consts_coms,only: r8
use mem_basic,  only: rho
use oname_coms, only: nl

implicit none

real(r8), intent(in) :: rho_old(mza,mwa)

integer  :: n,iw,j,k,mrl
real(r8) :: dtl

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

   !$omp parallel do collapse(2) private(n,j,iw,dtl,k)
   do n = 1, num_scalar
      do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

         dtl = dtlm(itab_w(iw)%mrlw)

         do k = lpw(iw),mza
            scalar_tab(n)%var_p(k,iw) = ( scalar_tab(n)%var_p(k,iw) * rho_old(k,iw)  &
                                         + dtl * scalar_tab(n)%var_t(k,iw) ) / rho(k,iw)

            ! With monotonic advection and splitting, positive-definite scalars
            ! should remain positive-definite to machine precision. Here we
            ! ensure that tiny negative concentrations aren't produeced.

            if ( scalar_tab(n)%pdef .and. nl%iscal_monot > 0 .and. &
                 nl%split_scalars > 0) then

               scalar_tab(n)%var_p(k,iw) = max( scalar_tab(n)%var_p(k,iw), 0.0 )

            endif

            ! Without monotonic advection and splitting, negative concentrations
            ! can be produced that are significant. If we set those to zero, we 
            ! will no longer conserve mass. 

         enddo

      enddo
   enddo
   !$omp end parallel do

endif

end subroutine predtr

!==========================================================================

subroutine comp_alpha_press(mrl, alpha_press)

  use mem_grid,    only: lpw, mza, mwa
  use mem_ijtabs,  only: jtab_w, jtw_prog
  use consts_coms, only: pc1, rdry, rvap, cpocv
  use mem_basic,   only: sh_w, sh_v, theta, thil

  implicit none

  integer, intent(in)  :: mrl
  real,    intent(out) :: alpha_press(mza,mwa)
  integer              :: j, iw, k

! Evaluate alpha coefficient for pressure

  !$omp parallel do private(iw,k) 
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        alpha_press(k,iw) = pc1 * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap) &
                          * theta(k,iw) / thil(k,iw)) ** cpocv
     enddo

  enddo
  !$omp end parallel do

end subroutine comp_alpha_press

!==========================================================================

subroutine omic_update_v_mom(mrl)
  
  use mem_grid,    only: lpv, mza
  use mem_ijtabs,  only: jtab_v, jtv_prog, itab_v
  use mem_basic,   only: vmc, vc, rho

  implicit none
  
  integer, intent(in) :: mrl
  integer             :: j, iv, iw1, iw2, k

! Update V momentum after microphysics update of rho

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1, jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vmc(k,iv) = 0.5 * vc(k,iv) * (rho(k,iw1) + rho(k,iw2))
     enddo
  enddo
  !$omp end parallel do

end subroutine omic_update_v_mom

!==========================================================================

subroutine bubble()

use mem_basic, only: thil, theta
!use mem_grid
!use mem_ijtabs

implicit none

integer :: k

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

end subroutine bubble

!==========================================================================

subroutine predtr_split(mrl,rho_old)

  use var_tables,  only: num_scalar, scalar_tab
  use mem_ijtabs,  only: jtab_w, itab_w, jtw_prog
  use mem_grid,    only: mza, mwa, lpw
  use misc_coms,   only: dtlm
  use consts_coms, only: r8

  implicit none

  integer,  intent(in) :: mrl
  real(r8), intent(in) :: rho_old(mza,mwa)
  integer              :: iw, j, n, k
  real(r8)             :: dtl

  ! Step scalars from t to t+1 in a time-split sub step

  if (mrl > 0) then

     !$omp parallel do collapse(2) private(n,j,iw,dtl,k)
     do n = 1, num_scalar
        do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

           dtl = dtlm(itab_w(iw)%mrlw)

           do k = lpw(iw), mza
              scalar_tab(n)%var_p(k,iw) = scalar_tab(n)%var_p(k,iw)  &
                                        + dtl * scalar_tab(n)%var_t(k,iw) / rho_old(k,iw)
              scalar_tab(n)%var_t(k,iw) = 0.0

              if (scalar_tab(n)%pdef) then
                 scalar_tab(n)%var_p(k,iw) = max( scalar_tab(n)%var_p(k,iw), 0.0)
              endif
           enddo

        enddo
     enddo
     !$omp end parallel do

  endif

  ! TODO: check if updated scalars are positive-definite??

end subroutine predtr_split
