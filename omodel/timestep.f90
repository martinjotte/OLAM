subroutine timestep()

use misc_coms,   only: time8, time_istp8, time_istp8p, time_bias, io6, &
                       nqparm, initial, ilwrtyp, iswrtyp, dtsm, i_o3, &
                       iparallel, s1900_init, s1900_sim, do_chem, &
                       nrk_wrtv, nrk_scal
use mem_ijtabs,  only: nstp, istp, mrls, mrl_begl, mrl_endl
use mem_nudge,   only: nudflag, nudnxp, o3nudflag
use mem_grid,    only: mza, mwa
use micro_coms,  only: miclevel
use leaf_coms,   only: isfcl
use mem_basic,   only: thil
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
use obnd,        only: set_scalars_lbc, set_scalars_bottom
use oname_coms,  only: nl
use mem_flux_accum, only: flux_accum
use consts_coms, only: r8
use vel_t3d,     only: diag_uzonal_umerid
use mem_megan,   only: megan_avg_temp
use emis_defn,   only: get_emis
use depv_defn,   only: get_depv
use wrtv_rk,     only: prog_wrtv_rk
use wrtv_orig,   only: prog_wrtv_orig
use check_nan,   only: check_nans
use pbl_drivers, only: pbl_driver, comp_horiz_k
use olam_mpi_sfc,only: mpi_send_wsfc, mpi_recv_wsfc
use hcane_rz,    only: ncycle_hurrinit, icycle_hurrinit, timmax_hurrinit, &
                       vortex_add_thetapert
use obs_nudge_mod,only: obs_nudge
!use oplot_coms,  only: op

implicit none

integer :: jstp

real(r8) :: rho_old(mza,mwa) ! density at beginning of long timestep [kg/m^3]
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

   ! call check_nans(15)

endif

do jstp = 1,nstp  ! nstp = no. of finest-grid-level aco steps in dtlm
   istp = jstp

   ! Bypass all processes except land/lake/sea if running groundwater spin-up simulation

   if (nl%igw_spinup == 1) go to 1400

   ! Bypass all processes except microphyics if running microphysics parcel simulation

   if (nl%test_case >= 901 .and. nl%test_case <= 950) go to 1311

   if (mrl_begl(istp) > 0) then
      call tend0(rho_old)
      call comp_alpha_press()

      call sea_spray()
      call dust_src()
   endif

   ! call check_nans(1)

   if (any( nqparm(1:mrls) > 0 )) then
      call cuparm_driver()

      if (isfcl == 1 .and. mrl_begl(istp) > 0) then
         call surface_cuparm_flux()
      endif
   endif

   ! update cloud fraction before radiation

   if (mrl_begl(istp) > 0) then
      call calc_3d_cloud_fraction()
   endif

   ! call check_nans(2)

   if (ilwrtyp + iswrtyp > 0) then
      call radiate()
   endif

   if (mrl_begl(istp) > 0) then

      if (isfcl > 0) then
         call surface_driver()
      else
         call surface_turb_flux()
      endif

      ! Add incremental axisymmetric potential temperature perturbation
      ! inside the hurricane core to increase vortex intensity

      if (icycle_hurrinit < ncycle_hurrinit .or. &
         (ncycle_hurrinit == 1 .and. time_istp8p < timmax_hurrinit)) then

         call vortex_add_thetapert()
      endif

      ! Nudging tendencies

      if (initial == 2 .and. nudflag == 1) then
         if (nudnxp == 0) then
            call  obs_nudge()
         else
            call spec_nudge()
         endif
      endif

      ! Nudging of ozone if it is in the scalar table

      if (initial == 2 .and. o3nudflag == 1 .and. i_o3 > 0) then
         call obs_nudge_o3()
      endif

      ! CMAQ emissions and deposition

      if (do_chem == 1) then
         call get_emis()
         call get_depv()
      endif

      ! Long timestep PBL tendencies

      call pbl_driver()

      ! and lateral friction with terrain with shaved cells

      if (isfcl == 1) call lateral_friction()

      ! Computation of horizontal K's

      call comp_horiz_k()

   endif  ! mrl_begl(istp) > 0

   ! call check_nans(11)

   ! Call olam_dcmip_phys, which is the OLAM interface to DCMIP auxiliary
   ! physics subroutine that provides tendencies to some model fields

   if (mrl_begl(istp) > 0) then
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
         call olam_dcmip2016_phys()
      endif
   endif

   ! Call olam_dcmip_terminator, which is the OLAM interface to DCMIP auxiliary
   ! chemistry subroutine terminator that provides chemical tendencies.

   if (mrl_begl(istp) > 0) then
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

   ! call check_nans(11,rvara1=alpha_press)

   ! Bypass call to thiltend_long if using prescribed flow for DCMIP tests
   !--------------------------------------
   if (nl%test_case == 11 .or. &
       nl%test_case == 12 .or. &
       nl%test_case == 13) go to 33
   !--------------------------------------

   if (nrk_wrtv == 1) then
      call prog_wrtv_orig()
   else
      call prog_wrtv_rk()
   endif

   33 continue  ! test_case == 11, 12, or 13

   if (mrl_endl(istp) > 0) then
      call diag_uzonal_umerid()
   endif

   ! call check_nans(12)

   !--------------------------------------
   ! Get flow update if using DCMIP prescribed flow

   if (nl%test_case == 11 .or. &
       nl%test_case == 12) then

      ! Set time0 to half-forward time (of grid-1 short timestep)
      ! for time-centered advective tendencies

      time0 = time8 + .5 * dtsm

      call olam_dcmip_prescribedflow(time0,rho_old)

   endif
   !--------------------------------------

   ! call check_nans(12,rvara1=alpha_press)

   ! Bypass call to timeavg_momsc if using prescribed flow for DCMIP tests
   !--------------------------------------
   if (nl%test_case == 11 .or. &
       nl%test_case == 12 .or. &
       nl%test_case == 13) go to 34
   !--------------------------------------

   call timeavg_momsc()

   34 continue

   ! call check_nans(13)

   if (mrl_endl(istp) > 0) then
      if (nrk_scal == 1) then
         call scalar_transport_orig(rho_old)
      else
         call scalar_transport_rk(rho_old)
      endif
   endif

   ! call check_pos(1)
   ! call check_nans(14)

   ! Special diagnosis of water vapor for DCMIP tests;
   ! thil is dry theta in these cases

   if (nl%test_case == 110 .or. &
       nl%test_case == 111 .or. &
       nl%test_case == 112 .or. &
       nl%test_case == 113 .or. &
       nl%test_case == 114 .or. &
       nl%test_case == 121 .or. &
       nl%test_case == 122 .or. &
       nl%test_case == 131) then

      if (mrl_endl(istp) > 0) then
         call thermo_dcmip()
         go to 35
      endif
   endif

   if (mrl_endl(istp) > 0 .and. miclevel /= 3) then
      call thermo()
   endif

   1311 continue

   if (mrl_endl(istp) > 0 .and. miclevel == 3) then
      call micro()
   endif

   ! call check_pos(2)
   ! call check_nans(16)

   ! Bypass all processes except microphyics if running microphysics parcel simulation

   if (nl%test_case >= 901 .and. nl%test_case <= 950) go to 1312

   ! Call atmospheric chemistry here

   if (do_chem == 1) then
      call cmaq_driver()
   endif

   35 continue

   ! call check_pos(3)
   ! call check_nans(20,rvara1=alpha_press)

   if (mrl_endl(istp) > 0) then

      ! Bottom boundary for scalars
      call set_scalars_bottom()

      ! Parallel send/recv of scalars
      if (iparallel == 1) call mpi_send_w(rvara1=thil, scalars='S')  ! Send scalars
      if (iparallel == 1) call mpi_recv_w(rvara1=thil, scalars='S')  ! Recv scalars

      ! Lateral boundary copy of scalars for limited-area runs
      call set_scalars_lbc()

   endif

   1400 continue  ! nl%igw_spinup == 1

   if (isfcl > 0 .and. mrl_endl(jstp) > 0) then

      if (nl%igw_spinup /= 1) then
         call sfcg_avgatm()
         if (do_chem == 1) call megan_avg_temp()
      endif

   endif

   ! call check_nans(21,rvara1=alpha_press)

   1312 continue

   time_istp8  = time8 + istp * dtsm  ! Update precise time
   time_istp8p = time_istp8 + time_bias
   s1900_sim   = s1900_init + time_istp8

   ! Add current fluxes to time integrals

   call flux_accum()

   ! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
   !if (mod(time8,op%frqplt) < dtlm .and. istp == nstp+1000) then
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
   !   op%extfld(:,:) = (rr_w(:,:) - rr_v(:,:)) * 1.e3
   !   op%extfldname = 'RR_TOTCOND'
   !   call plot_fields(3)
   !   deallocate (op%extfld)
   !
   !endif
   ! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

enddo

! For ncar dcmip test cases, compute error norms

if (nl%test_case >= 10 .and. nl%test_case < 900) then
   call diagn_global_dcmip()
endif

end subroutine timestep

!===============================================================================

subroutine modsched()

  use mem_ijtabs,  only: nstp, mrl_endl, mrl_begl
  use misc_coms,   only: io6, dtlong, nacoust, dtlm, dtsm
  use leaf_coms,   only: dt_leaf
  use lake_coms,   only: dt_lake
  use sea_coms,    only: dt_sea
  use sea_swm,     only: dt_swm, niter_swm

  implicit none

  integer :: jstp

  nstp = nacoust

  write(io6,'(/,a)') '=== Timestep Schedule ===='
  write(io6,'(a,/)') '              jstp    mrl_begl  mrl_endl'

  ! Set timestep lengths

  dtlm = dtlong
  dtsm = dtlm / real(nacoust)

  dt_leaf = dtlm
  dt_sea  = dtlm
  dt_lake = dtlm
  dt_swm  = dtlm / real(niter_swm)

  ! Allocate mrl-schedule arrays

  allocate (mrl_begl(nstp)) ; mrl_begl = 0
  allocate (mrl_endl(nstp)) ; mrl_endl = 0

  ! Fill acoustic timestep sub-cycling flags for processes
  ! to carry out for each jstp value

  do jstp = 1,nstp
     if (mod(jstp-1, nacoust) == 0) mrl_begl(jstp) = 1
     if (mod(jstp  , nacoust) == 0) mrl_endl(jstp) = 1

     write(io6,333) jstp,mrl_begl(jstp),mrl_endl(jstp)
     333 format('modsched0 ',3i10)
  enddo

end subroutine modsched

!==========================================================================

subroutine tend0(rho_old)

  use var_tables,  only: scalar_tab, num_scalar
  use mem_tend,    only: thilt, vmxet, vmyet, vmzet, vmt
  use misc_coms,   only: nrk_scal
  use mem_sfcg,    only: mwsfc, sfcg
  use mem_basic,   only: vmsc, wmsc, rho, vxesc, vyesc, vzesc
  use mem_grid,    only: mza, mwa, mva, lpw, lpv
  use consts_coms, only: r8

  implicit none

  real(r8), intent(out) :: rho_old(mza,mwa)
  integer               :: iw, iv, n, k, iwsfc

  !$omp parallel
  !$omp do private(k,n)
  do iw = 2, mwa

     do k = lpw(iw), mza
        vmxet  (k,iw) = 0.0
        vmyet  (k,iw) = 0.0
        vmzet  (k,iw) = 0.0
        thilt  (k,iw) = 0.0
        rho_old(k,iw) = rho(k,iw)
     enddo

     do k = lpw(iw)-1, mza
        wmsc(k,iw) = 0.0
     enddo

     do n = 1, num_scalar
        do k = lpw(iw), mza
           scalar_tab(n)%var_t(k,iw) = 0.0
        enddo
     enddo

     if (nrk_scal == 1) then
        do k = lpw(iw), mza
           vxesc(k,iw) = 0.0
           vyesc(k,iw) = 0.0
           vzesc(k,iw) = 0.0
        enddo
     endif

  enddo
  !$omp end do nowait

  !$omp do private(k)
  do iv = 2, mva
     do k = lpv(iv), mza
        vmt (k,iv) = 0.0
        vmsc(k,iv) = 0.0
     enddo
  enddo
  !$omp end do nowait

  if (allocated(sfcg%pcpg)) then
     !$omp do
     do iwsfc = 2, mwsfc
        sfcg%pcpg (iwsfc) = 0.
        sfcg%qpcpg(iwsfc) = 0.
        sfcg%dpcpg(iwsfc) = 0.
     enddo
     !$omp end do nowait
  endif

  !$omp end parallel

end subroutine tend0

!==========================================================================

subroutine comp_alpha_press()

  use mem_grid,    only: lpw, lpv, mza, dzim, dniv, zfacit
  use mem_ijtabs,  only: jtab_w, jtw_prog, jtab_v, jtv_prog, itab_v
  use consts_coms, only: pc1, rdry, rvap, cpocv, gravo2
  use mem_basic,   only: rr_v, rr_w, theta, thil, alpha_press, &
                         pwfac, pvfac

  implicit none

  integer :: j, iw, k, iv, iw1, iw2
  real    :: rw

  !$omp parallel
  !$omp do private(iw,k,rw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Evaluate alpha coefficient for pressure

     do k = lpw(iw), mza
        alpha_press(k,iw) = pc1 * ( (rdry + rvap * rr_v(k,iw)) &
                                  * theta(k,iw) / thil(k,iw) ) ** cpocv
     enddo

     ! Factor muliplying pressure gradient in vertical eq. of motion
     ! rho_dry / rho_total / dz

     do k = lpw(iw), mza-1
      ! pwfac(k,iw) = dzim(k) * ( zwgt_bot(k) / (1. + rr_w(k  ,iw)) &
      !                         + zwgt_top(k) / (1. + rr_w(k+1,iw)) )
        rw          = 0.5 * (rr_w(k+1,iw) + rr_w(k,iw))
        pwfac(k,iw) = dzim(k) * (1.0 - rw + rw * rw)
     enddo

  enddo
  !$omp end do nowait

  !$omp do private(iv,iw1,iw2,k,rw)
  do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! Factor muliplying pressure gradient in horizontal eq. of motion
     ! rho_dry / rho_total / dx

     do k = lpv(iv), mza
      ! pvfac(k,iv) = dnivo2(iv) * zfacit(k) * ( 1. / (1. + rr_w(k,iw1)) &
      !                                        + 1. / (1. + rr_w(k,iw2)) )
        rw          = 0.5 * (rr_w(k,iw1) + rr_w(k,iw2))
        pvfac(k,iv) = dniv(iv) * zfacit(k) * (1.0 - rw + rw * rw)
     enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

end subroutine comp_alpha_press

!===========================================================================

subroutine timeavg_momsc()

  use mem_ijtabs, only: istp, mrl_endl
  use mem_grid,   only: mza, mva, mwa, lpw, lpv
  use misc_coms,  only: nacoust
  use mem_basic,  only: vmsc, wmsc

  implicit none

  integer :: k, iv, iw
  real    :: acoi

  if (mrl_endl(istp) > 0) then
     acoi = 1.0 / nacoust

     !$omp parallel
     !$omp do private(k)
     do iv = 2, mva
        do k = lpv(iv), mza
           vmsc(k,iv) = vmsc(k,iv) * acoi
        enddo
     enddo
     !$omp end do nowait

     !$omp do private(k)
     do iw = 2, mwa
        do k = lpw(iw), mza-1
           wmsc(k,iw) = wmsc(k,iw) * acoi
        enddo
     enddo
     !$omp end do nowait
     !$omp end parallel

  endif

end subroutine timeavg_momsc

!==========================================================================

subroutine bubble()

  use mem_basic, only: thil, theta

  implicit none

  integer :: k

  do k = 2, 2
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
