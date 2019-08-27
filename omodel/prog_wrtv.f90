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

module wrtv_mem

  use consts_coms, only: r8
  implicit none

  real, allocatable :: vortp   (:,:)
  real, allocatable :: vortp_t (:,:)
  real, allocatable :: vortn   (:,:)

  real, allocatable :: thil_old(:,:)
  real, allocatable :: wmca    (:,:)

  logical, save     :: rotational = .false.
  integer, save     :: iflag      = 0

  real,     allocatable :: alpha_press(:,:)
  real(r8), allocatable :: rd_rt_w    (:,:)

  private
  public alloc_prog_wrtv_mem, prog_wrtv, comp_alpha_press, alpha_press

contains

!===============================================================================

subroutine alloc_prog_wrtv_mem()

  use mem_grid,   only: mma, mwa, mva, mza
  use misc_coms,  only: rinit, rinit8
  use oname_coms, only: nl

  implicit none

  if (nl%expnme(1:1) == 'R') then
     rotational = .true.
  else
     rotational = .false.
  endif

  if (rotational) then
     allocate(vortp  (mza,mma)) ; vortp   = rinit
     allocate(vortp_t(mza,mwa)) ; vortp_t = rinit
     allocate(vortn  (mza,mva)) ; vortn   = rinit
  endif

  allocate(thil_old   (mza,mwa))
  allocate(alpha_press(mza,mwa)) ; alpha_press = rinit
  allocate(rd_rt_w    (mza,mwa)) ; rd_rt_w     = rinit8

  if (nl%ithil_monot + nl%ivelc_monot > 0) then
     allocate(wmca(mza,mwa)) ; wmca(:,:) = 0.0
  endif

  if (nl%ithil_monot == 1 .or. nl%ivelc_monot == 1) then
     iflag = 1
  else
     iflag = 2
  endif


end subroutine alloc_prog_wrtv_mem

!===============================================================================

subroutine prog_wrtv(vmsc, wmsc)

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

  use mem_ijtabs,   only: jtab_v, jtab_w, istp, &
                          mrl_ends, jtw_wadj, jtm_vadj, jtv_prog, &
                          jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp, &
                          itab_w, itab_v
  use mem_basic,    only: rho, thil, wc, press, vmp, vmc, vc, wmc, &
                          vxe, vye, vze
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, arw, arv
  use mem_tend,     only: thilt, vmxet, vmyet, vmzet
  use misc_coms,    only: iparallel, dtsm
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m, &
                          mpi_send_v, mpi_recv_v
  use oplot_coms,   only: op
  use obnd,         only: lbcopy_m, lbcopy_w, lbcopy_v
  use vel_t3d,      only: vel_t3d_hex
  use oname_coms,   only: nl
  use mem_adv,      only: dxps_w, dyps_w, dzps_w, dzzps_w, &
                          dxps_v, dyps_v, dzps_v, &
                          dxyps_v, dxxps_v, dyyps_v, &
                          gxps_scp, gyps_scp, gzps_scp, &
                          gxyps_scp, gxxps_scp, gyyps_scp, gzzps_scp, &
                          gxps_vxe, gyps_vxe, gzps_vxe, &
                          gxyps_vxe, gxxps_vxe, gyyps_vxe, gzzps_vxe, &
                          gxps_vye, gyps_vye, gzps_vye, &
                          gxyps_vye, gxxps_vye, gyyps_vye, gzzps_vye, &
                          gxps_vze, gyps_vze, gzps_vze, &
                          gxyps_vze, gxxps_vze, gyyps_vze, gzzps_vze
  use mem_thuburn,  only: comp_cfls, comp_and_apply_monot_limits, &
                          comp_and_apply_pd_limits
  use grad_lib

  implicit none

  real, intent(inout) :: vmsc(mza,mva)
  real, intent(inout) :: wmsc(mza,mwa)

  integer :: j, iv, iw, k, mrl, kd, iwd, iw1, iw2

! automatic arrays

  integer :: iwdepv(mza,mva) ! donor cell IW index for V face
  integer :: kdepw(mza,mwa)  ! donor cell K index for W face

  real :: vmcf, dt
  real :: vmcfa(mza,mva) ! Time-extrapolated VMC
  real :: vcf(mza)

! Volume-weighted T-cell momentum component tendencies from advective
! transport; evaluated and applied on acoustic timestep

  real :: vmxet_volt(mza,mwa)
  real :: vmyet_volt(mza,mwa)
  real :: vmzet_volt(mza,mwa)

  real :: thil_upv(mza,mva) ! Upstreamed THIL at each V interface
  real :: vxe_upv (mza,mva) ! Upstreamed VXE  at each V interface
  real :: vye_upv (mza,mva) ! Upstreamed VYE  at each V interface
  real :: vze_upv (mza,mva) ! Upstreamed VZE  at each V interface

  real :: thil_upw(mza,mwa) ! Upstreamed THIL at each W level
  real :: vxe_upw (mza,mwa) ! Upstreamed VXE  at each W level
  real :: vye_upw (mza,mwa) ! Upstreamed VYE  at each W level
  real :: vze_upw (mza,mwa) ! Upstreamed VZE  at each W level

  real :: thilt_short(mza,mwa)
  real :: vmxet_short(mza,mwa)
  real :: vmyet_short(mza,mwa)
  real :: vmzet_short(mza,mwa)

! THIS ROUTINE IS PERFORMED AT THE END OF EACH SMALL TIMESTEP

  mrl = mrl_ends(istp)
  if (mrl == 0) return

! INCLUDE THE LONG TIMESTEP THIL AND MOMENTUM TENDENCIES IN EACH SHORT TIMESTEP

  !$omp parallel private(vcf)
  !$omp do private(iw,dt,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     dt = dtsm(itab_w(iw)%mrlw)

     do k = lpw(iw), mza
        thilt_short(k,iw) = thilt(k,iw)
        vmxet_short(k,iw) = vmxet(k,iw)
        vmyet_short(k,iw) = vmyet(k,iw)
        vmzet_short(k,iw) = vmzet(k,iw)
     enddo

     if (nl%ithil_monot + nl%ivelc_monot > 0) then
        do k = lpw(iw), mza-1
           wmca(k,iw) = wmc(k,iw) * arw(k,iw)
        enddo
     endif

     call donorpointw(iw, dt, wc(:,iw), vxe(:,iw), vye(:,iw), vze(:,iw), kdepw(:,iw))

  enddo
  !$omp end do nowait

  !$omp do private(iv,iw1,iw2,k,dt,vmcf)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     ! Extrapolate VM to time T + 1/2; update VMP

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vmcf        = 1.5 * vmc(k,iv) - 0.5 * vmp(k,iv)
        vmsc (k,iv) = vmsc(k,iv) + vmcf
        vmcfa(k,iv) = vmcf * arv(k,iv)
        vmp  (k,iv) = vmc (k,iv)
        vcf  (k)    = 2.0 * vmcf / real( rho(k,iw1) + rho(k,iw2) )
     enddo

     dt = dtsm(itab_v(iv)%mrlv)

     call donorpointv(iv, dt, vcf, vxe, vye, vze, iwdepv(:,iv))

  enddo
  !$omp end do nowait
  !$omp end parallel

  if (nl%ithil_monot + nl%ivelc_monot > 0) then

      ! Compute CFL numbers needed by Thuburn monotonic scheme
      call comp_cfls( mrl, dtsm, vmcfa, wmca, rho, iflag, do_check=.false. )

  endif

! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
! (Example of how to plot "external" field; one not available in module memory)
!
! if (mod(time8,op%frqplt) < dtlm(1) .and. istp == 900) then
!    allocate (op%extfld(mza,mwa))
!    op%extfld(:,:) = vxe(:,:)
!    op%extfldname = 'VXE'
!    call plot_fields(11)
!    deallocate (op%extfld)
! endif
! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

! EVALUATE HORIZONTAL GRADIENTS OF THIL, VXE, VYE, AND VZE FOR BEGS

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     if (nl%horiz_adv_order <= 2) then
        call grad_t2d  (iw, thil, gxps_scp(:,iw),  gyps_scp(:,iw))
        call grad_t2d  (iw, vxe,  gxps_vxe(:,iw),  gyps_vxe(:,iw))
        call grad_t2d  (iw, vye,  gxps_vye(:,iw),  gyps_vye(:,iw))
        call grad_t2d  (iw, vze,  gxps_vze(:,iw),  gyps_vze(:,iw))
     else
        call grad_t2d_quad(iw, thil, gxps_scp (:,iw), gyps_scp (:,iw), &
                                     gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
        call grad_t2d_quad(iw, vxe,  gxps_vxe (:,iw), gyps_vxe (:,iw), &
                                     gxxps_vxe(:,iw), gxyps_vxe(:,iw), gyyps_vxe(:,iw))
        call grad_t2d_quad(iw, vye,  gxps_vye (:,iw), gyps_vye (:,iw), &
                                     gxxps_vye(:,iw), gxyps_vye(:,iw), gyyps_vye(:,iw))
        call grad_t2d_quad(iw, vze,  gxps_vze (:,iw), gyps_vze (:,iw), &
                                      gxxps_vze(:,iw), gxyps_vze(:,iw), gyyps_vze(:,iw))
     endif
  enddo
  !$omp end parallel do

! MPI send of THIL, VXE, VYE, and VZE gradient components

  if (iparallel == 1) then

     if (nl%horiz_adv_order <= 2) then
        call mpi_send_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
                             rvara3=gxps_vxe, rvara4=gyps_vxe, &
                             rvara5=gxps_vye, rvara6=gyps_vye, &
                             rvara7=gxps_vze, rvara8=gyps_vze  )
     else
        call mpi_send_w(mrl, rvara1=gxps_scp,   rvara2=gyps_scp,                      &
                             rvara3=gxxps_scp,  rvara4=gxyps_scp,  rvara5=gyyps_scp,  &
                             rvara6=gxps_vxe,   rvara7=gyps_vxe,                      &
                             rvara8=gxxps_vxe,  rvara9=gxyps_vxe,  rvara10=gyyps_vxe, &
                             rvara11=gxps_vye,  rvara12=gyps_vye,                     &
                             rvara13=gxxps_vye, rvara14=gxyps_vye, rvara15=gyyps_vye, &
                             rvara16=gxps_vze,  rvara17=gyps_vze,                     &
                             rvara18=gxxps_vze, rvara19=gxyps_vze, rvara20=gyyps_vze  )
     endif

  endif ! iparallel == 1

! EVALUATE VERTICAL GRADIENT OF THIL, VXE, VYE, AND VZE FOR BEGS

  !$omp parallel
  !$omp do private(iw)
  do j = 1,jtab_w(jtw_wadj)%jend(mrl); iw = jtab_w(jtw_wadj)%iw(j)

     thil_old(2:mza,iw) = thil(2:mza,iw)

    call grad_z_quad(iw, thil(:,iw), gzps_scp(:,iw), gzzps_scp(:,iw))
    call grad_z_quad(iw, vxe (:,iw), gzps_vxe(:,iw), gzzps_vxe(:,iw))
    call grad_z_quad(iw, vye (:,iw), gzps_vye(:,iw), gzzps_vye(:,iw))
    call grad_z_quad(iw, vze (:,iw), gzps_vze(:,iw), gzzps_vze(:,iw))

  enddo
  !$omp end do

! EVALUATE UPWINDED THIL, VXE, VYE, AND VZE AT EACH W FACE FOR
! VERTICAL FLUX COMPUTATION

  !$omp do private(iw,k,kd)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza-1
        kd = kdepw(k,iw)

        thil_upw(k,iw) =                      thil(kd,iw) &
                       +  dxps_w(k,iw) *  gxps_scp(kd,iw) &
                       +  dyps_w(k,iw) *  gyps_scp(kd,iw) &
                       +  dzps_w(k,iw) *  gzps_scp(kd,iw) &
                       + dzzps_w(k,iw) * gzzps_scp(kd,iw)

        vxe_upw(k,iw)  =                       vxe(kd,iw) &
                       +  dxps_w(k,iw) *  gxps_vxe(kd,iw) &
                       +  dyps_w(k,iw) *  gyps_vxe(kd,iw) &
                       +  dzps_w(k,iw) *  gzps_vxe(kd,iw) &
                       + dzzps_w(k,iw) * gzzps_vxe(kd,iw)

        vye_upw(k,iw)  =                      vye(kd,iw) &
                       +  dxps_w(k,iw) * gxps_vye(kd,iw) &
                       +  dyps_w(k,iw) * gyps_vye(kd,iw) &
                       +  dzps_w(k,iw) * gzps_vye(kd,iw) &
                       + dzzps_w(k,iw) * gzzps_vye(kd,iw)

        vze_upw(k,iw)  =                       vze(kd,iw) &
                       +  dxps_w(k,iw) *  gxps_vze(kd,iw) &
                       +  dyps_w(k,iw) *  gyps_vze(kd,iw) &
                       +  dzps_w(k,iw) *  gzps_vze(kd,iw) &
                       + dzzps_w(k,iw) * gzzps_vze(kd,iw)
     enddo

  enddo
  !$omp end do nowait
  !$omp end parallel

! Finish MPI recv of THIL, VXE, VYE, and VZE gradient components

  if (iparallel == 1) then
     if (nl%horiz_adv_order <= 2) then
        call mpi_recv_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
                             rvara3=gxps_vxe, rvara4=gyps_vxe, &
                             rvara5=gxps_vye, rvara6=gyps_vye, &
                             rvara7=gxps_vze, rvara8=gyps_vze  )
     else
        call mpi_recv_w(mrl, rvara1=gxps_scp,   rvara2=gyps_scp,                      &
                             rvara3=gxxps_scp,  rvara4=gxyps_scp,  rvara5=gyyps_scp,  &
                             rvara6=gxps_vxe,   rvara7=gyps_vxe,                      &
                             rvara8=gxxps_vxe,  rvara9=gxyps_vxe,  rvara10=gyyps_vxe, &
                             rvara11=gxps_vye,  rvara12=gyps_vye,                     &
                             rvara13=gxxps_vye, rvara14=gxyps_vye, rvara15=gyyps_vye, &
                             rvara16=gxps_vze,  rvara17=gyps_vze,                     &
                             rvara18=gxxps_vze, rvara19=gxyps_vze, rvara20=gyyps_vze  )
     endif
  endif

! Lateral boundary copy of THIL, VXE, VYE, and VZE gradient components

  if (nl%horiz_adv_order <= 2) then
     call lbcopy_w(mrl, a1=gxps_scp, a2=gyps_scp, &
                        a3=gxps_vxe, a4=gyps_vxe, &
                        a5=gxps_vye, a6=gyps_vye, &
                        a7=gxps_vze, a8=gyps_vze  )
  else
     call lbcopy_w(mrl, a1=gxps_scp,   a2=gyps_scp,                  &
                        a3=gxxps_scp,  a4=gxyps_scp,  a5=gyyps_scp,  &
                        a6=gxps_vxe,   a7=gyps_vxe,                  &
                        a8=gxxps_vxe,  a9=gxyps_vxe,  a10=gyyps_vxe, &
                        a11=gxps_vye,  a12=gyps_vye,                 &
                        a13=gxxps_vye, a14=gxyps_vye, a15=gyyps_vye, &
                        a16=gxps_vze,  a17=gyps_vze,                 &
                        a18=gxxps_vze, a19=gxyps_vze, a20=gyyps_vze  )
  endif

! EVALUATE UPWINDED THIL, VXE, VYE, AND VZE AT EACH V FACE FOR
! HORIZONTAL FLUX COMPUTATION

  !$omp parallel do private(iv,k,iwd)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     if (nl%horiz_adv_order <= 2) then

        do k = lpv(iv), mza
           iwd = iwdepv(k,iv)

           thil_upv(k,iv) =                    thil(k,iwd) &
                          + dxps_v(k,iv) * gxps_scp(k,iwd) &
                          + dyps_v(k,iv) * gyps_scp(k,iwd) &
                          + dzps_v(k,iv) * gzps_scp(k,iwd)

           vxe_upv(k,iv)  =                     vxe(k,iwd) &
                          + dxps_v(k,iv) * gxps_vxe(k,iwd) &
                          + dyps_v(k,iv) * gyps_vxe(k,iwd) &
                          + dzps_v(k,iv) * gzps_vxe(k,iwd)

           vye_upv(k,iv)  =                     vye(k,iwd) &
                          + dxps_v(k,iv) * gxps_vye(k,iwd) &
                          + dyps_v(k,iv) * gyps_vye(k,iwd) &
                          + dzps_v(k,iv) * gzps_vye(k,iwd)

           vze_upv(k,iv)  =                     vze(k,iwd) &
                          + dxps_v(k,iv) * gxps_vze(k,iwd) &
                          + dyps_v(k,iv) * gyps_vze(k,iwd) &
                          + dzps_v(k,iv) * gzps_vze(k,iwd)
        enddo

     else

        do k = lpv(iv), mza
           iwd = iwdepv(k,iv)

           thil_upv(k,iv) =                      thil(k,iwd) &
                          +  dxps_v(k,iv) *  gxps_scp(k,iwd) &
                          + dxxps_v(k,iv) * gxxps_scp(k,iwd) &
                          + dxyps_v(k,iv) * gxyps_scp(k,iwd) &
                          +  dyps_v(k,iv) *  gyps_scp(k,iwd) &
                          + dyyps_v(k,iv) * gyyps_scp(k,iwd) &
                          +  dzps_v(k,iv) *  gzps_scp(k,iwd)

           vxe_upv(k,iv)  =                       vxe(k,iwd) &
                          +  dxps_v(k,iv) *  gxps_vxe(k,iwd) &
                          + dxxps_v(k,iv) * gxxps_vxe(k,iwd) &
                          + dxyps_v(k,iv) * gxyps_vxe(k,iwd) &
                          +  dyps_v(k,iv) *  gyps_vxe(k,iwd) &
                          + dyyps_v(k,iv) * gyyps_vxe(k,iwd) &
                          +  dzps_v(k,iv) *  gzps_vxe(k,iwd)

           vye_upv(k,iv)  =                       vye(k,iwd) &
                          +  dxps_v(k,iv) *  gxps_vye(k,iwd) &
                          + dxxps_v(k,iv) * gxxps_vye(k,iwd) &
                          + dxyps_v(k,iv) * gxyps_vye(k,iwd) &
                          +  dyps_v(k,iv) *  gyps_vye(k,iwd) &
                          + dyyps_v(k,iv) * gyyps_vye(k,iwd) &
                          +  dzps_v(k,iv) *  gzps_vye(k,iwd)

           vze_upv(k,iv)  =                       vze(k,iwd) &
                          +  dxps_v(k,iv) *  gxps_vze(k,iwd) &
                          + dxxps_v(k,iv) * gxxps_vze(k,iwd) &
                          + dxyps_v(k,iv) * gxyps_vze(k,iwd) &
                          +  dyps_v(k,iv) *  gyps_vze(k,iwd) &
                          + dyyps_v(k,iv) * gyyps_vze(k,iwd) &
                          +  dzps_v(k,iv) *  gzps_vze(k,iwd)
        enddo

     endif

  enddo
  !$omp end parallel do

! Compute the monotonic or positive-definite flux limiters and then apply them

  if (nl%ithil_monot == 1) then
     call comp_and_apply_monot_limits(mrl, thil, thil_upw, thil_upv, kdepw, iwdepv)
  elseif (nl%ithil_monot == 2) then
     call comp_and_apply_pd_limits(mrl, thil, thil_upw, thil_upv, kdepw, iwdepv)
  endif

  if (nl%ivelc_monot == 1) then
     call comp_and_apply_monot_limits(mrl, vxe, vxe_upw, vxe_upv, kdepw, iwdepv)
     call comp_and_apply_monot_limits(mrl, vye, vye_upw, vye_upv, kdepw, iwdepv)
     call comp_and_apply_monot_limits(mrl, vze, vze_upw, vze_upv, kdepw, iwdepv)
  endif

! MAIN LOOP OVER W COLUMNS FOR UPDATING WM, WC, RHO, THIL, AND PRESS

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Prognose vertical velocity, density, thil, and diagnose pressure
     call prog_wrt_begs( iw, vmcfa, wmsc,                      &
                         thil_upv, vxe_upv, vye_upv, vze_upv,  &
                         thil_upw, vxe_upw, vye_upw, vze_upw,  &
                         vmxet_volt, vmyet_volt, vmzet_volt,   &
                         thilt_short, vmxet_short, vmyet_short, vmzet_short )
  enddo
  !$omp end parallel do

! MPI SEND/RECV and LBC copy of quantities needed for prog_v:
! PRESS, RHO, VMXET_VOLT, VMYET_VOLT, and VMZET_VOLT

  if (iparallel == 1) then

     call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_volt,  rvara2=vmyet_volt,  rvara3=vmzet_volt, &
                     rvara4=vmxet_short, rvara5=vmyet_short, rvara6=vmzet_short )

     call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_volt,  rvara2=vmyet_volt,  rvara3=vmzet_volt, &
                     rvara4=vmxet_short, rvara5=vmyet_short, rvara6=vmzet_short )
  endif

  call lbcopy_w(mrl, a1=vmxet_volt,  a2=vmyet_volt,  a3=vmzet_volt,  &
                     a4=vmxet_short, a5=vmyet_short, a6=vmzet_short, &
                     d1=press,       d2=rho                          )

! COMPUTE TERMS FOR ROTATIONAL FORM OF HORIZONTAL MOMENTUM

  if (rotational) then
     call prep_rotational(mrl)
  endif

! MAIN LOOP OVER V POINTS TO UPDATE VMC AND VC

  !$omp parallel do private(iv)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     call prog_v_begs( iv, vmxet_volt, vmyet_volt, vmzet_volt,   &
                           vmxet_short, vmyet_short, vmzet_short )

  enddo

! MPI SEND/RECV and LBC of VC and VMC

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif
  call lbcopy_v(mrl, vmc=vmc, vc=vc)

end subroutine prog_wrtv

!=========================================================================

subroutine prog_wrt_begs( iw, vmcfa, wmsc,                      &
                          thil_upv, vxe_upv, vye_upv, vze_upv,  &
                          thil_upw, vxe_upw, vye_upw, vze_upw,  &
                          vmxet_volt, vmyet_volt, vmzet_volt,   &
                          thilt_short, vmxet_short, vmyet_short, vmzet_short )

  use mem_ijtabs,  only: itab_w
  use mem_basic,   only: wmc, rho, thil, wc, press, &
                         vxe, vye, vze, vxe2, vye2, vze2
  use misc_coms,   only: dtsm, deltax, nxp, initial, dn01d, th01d
  use consts_coms, only: cpocv, omega2, pi1, pio180, r8
  use mem_grid,    only: mza, mva, mwa, lpv, lpw, lve2, arw, dzim, gdzim, &
                         vnx, vny, vnz, wnx, wny, wnz, wnxo2, wnyo2, wnzo2, &
                         volti, volwi, glatw, glonw, dzt_top, dzt_bot, &
                         zwgt_top8, zwgt_bot8, gdz_abov8, gdz_belo8
  use tridiag,     only: tridiffo
  use oname_coms,  only: nl
  use mem_turb,    only: akmodx, akhodx
  use mem_rayf,    only: dorayfw, rayf_cofw, krayfw_bot, &
                         dorayf, rayf_cof, krayf_bot
  use mem_nudge,   only: rhot_nud, nudflag

  implicit none

  integer, intent(in) :: iw

  real, intent(in)    :: vmcfa(mza,mva)
  real, intent(inout) :: wmsc(mza,mwa)
  real, intent(in)    :: thil_upv(mza,mva)
  real, intent(in)    :: vxe_upv (mza,mva)
  real, intent(in)    :: vye_upv (mza,mva)
  real, intent(in)    :: vze_upv (mza,mva)

  real, intent(in)    :: thil_upw(mza,mwa)
  real, intent(in)    :: vxe_upw (mza,mwa)
  real, intent(in)    :: vye_upw (mza,mwa)
  real, intent(in)    :: vze_upw (mza,mwa)

  real, intent(inout) :: vmxet_volt(mza,mwa)
  real, intent(inout) :: vmyet_volt(mza,mwa)
  real, intent(inout) :: vmzet_volt(mza,mwa)

  real, intent(inout) :: thilt_short(mza,mwa)
  real, intent(inout) :: vmxet_short(mza,mwa)
  real, intent(inout) :: vmyet_short(mza,mwa)
  real, intent(inout) :: vmzet_short(mza,mwa)

  integer :: jv, iv, iwn
  integer :: k, ka, kbv, ksw
  integer :: npoly

  real :: dts
  real :: c6, c8, c9, c10
  real :: dirv, vmarv
  real :: del_rhothil, vmt1
  real :: rad0_swtc, rad_swtc, topo_swtc

  ! Vertical implicit scheme weighting parameters

  real, parameter :: fw = .55  ! wmc
  real, parameter :: fr = .80  ! rho
  real, parameter :: fp = .80  ! press

  real, parameter :: pc2 = fp * cpocv

  ! Automatic arrays

  real(r8) :: wmarw    (mza)
  real(r8) :: del_wmarw(mza)
  real(r8) :: hflux_rho(mza)
  real(r8) :: delex_rho(mza)
  real(r8) :: rhothil  (mza)
  real(r8) :: press_t  (mza)

  real :: delex_wm     (mza)
  real :: delex_rhothil(mza)
  real :: del_wm       (mza)
  real :: fwdel_wm     (mza)

  real :: wmt_other    (mza)

  real :: hflux_thil(mza)
  real :: hflux_vxe(mza)
  real :: hflux_vye(mza)
  real :: hflux_vze(mza)

  real :: hdiff_vxe(mza)
  real :: hdiff_vye(mza)
  real :: hdiff_vze(mza)

  real :: vflux_thil(mza)
  real :: vflux_vxe(mza)
  real :: vflux_vye(mza)
  real :: vflux_vze(mza)

  real(r8) :: b1(mza),b2(mza),b5(mza),b6(mza),b10(mza)
  real(r8) :: b7(mza),b8(mza),b9(mza),b11(mza),b12(mza),b13(mza),b14(mza)
  real(r8) :: b20(mza),b21(mza),b22(mza),b23(mza),b24(mza),b25(mza),b26(mza)
  real     :: b31(mza),b32(mza),b33(mza),b34(mza)

  real :: vxe1(mza)
  real :: vye1(mza)
  real :: vze1(mza)

  ka = lpw(iw)

  dts = dtsm(itab_w(iw)%mrlw)

  ! Set bottom & top vertical advective mass and heat fluxes to zero

  wmarw(1:ka-1) = 0.
  wmarw(mza)    = 0.

  vflux_thil(ka-1) = 0.
  vflux_thil(mza)  = 0.

  vflux_vxe(ka-1) = 0.
  vflux_vxe(mza)  = 0.

  vflux_vye(ka-1) = 0.
  vflux_vye(mza)  = 0.

  vflux_vze(ka-1) = 0.
  vflux_vze(mza)  = 0.

  ! Loop over W levels

  do k = ka,mza-1

     ! vertical mass flux for time level t

     wmarw(k) = wmc(k,iw) * arw(k,iw)

     ! vertical fluxes

     vflux_thil(k) = wmarw(k) * thil_upw(k,iw)
     vflux_vxe(k)  = wmarw(k) * vxe_upw (k,iw)
     vflux_vye(k)  = wmarw(k) * vye_upw (k,iw)
     vflux_vze(k)  = wmarw(k) * vze_upw (k,iw)

  enddo

  ! Initialize horizontal advection arrays to zero

  hflux_rho (1:mza) = 0.
  hflux_thil(1:mza) = 0.
  hflux_vxe (1:mza) = 0.
  hflux_vye (1:mza) = 0.
  hflux_vze (1:mza) = 0.

  hdiff_vxe (1:mza) = 0.
  hdiff_vye (1:mza) = 0.
  hdiff_vze (1:mza) = 0.

  ! Number of edges of this IW polygon

  npoly = itab_w(iw)%npoly

  ! Loop over V neighbors of this W cell

  do jv = 1,npoly

     iv   = itab_w(iw)%iv(jv)
     kbv  = lpv(iv)
     dirv = itab_w(iw)%dirv(jv)
     iwn  = itab_w(iw)%iw(jv)

     ! Loop over T levels

     do k = kbv, mza

        vmarv = dirv * vmcfa(k,iv)

        ! Sum horizontal advection fluxes over V faces

        hflux_rho (k) = hflux_rho(k)  + vmarv

        hflux_thil(k) = hflux_thil(k) + vmarv * thil_upv(k,iv) &
                      + akhodx(k,iv) * (thil_old(k,iwn) - thil_old(k,iw))

        hflux_vxe (k) = hflux_vxe(k) + vmarv * vxe_upv(k,iv)
        hflux_vye (k) = hflux_vye(k) + vmarv * vye_upv(k,iv)
        hflux_vze (k) = hflux_vze(k) + vmarv * vze_upv(k,iv)

        hdiff_vxe (k) = hdiff_vxe(k) + akmodx(k,iv) * (vxe(k,iwn) - vxe(k,iw))
        hdiff_vye (k) = hdiff_vye(k) + akmodx(k,iv) * (vye(k,iwn) - vye(k,iw))
        hdiff_vze (k) = hdiff_vze(k) + akmodx(k,iv) * (vze(k,iwn) - vze(k,iw))

     enddo

  enddo

  do k = ka, mza

     ! Add horizontal diffusive fluxes and Coriolis to momentum tendencies

     vmxet_short(k,iw) = vmxet_short(k,iw) + hdiff_vxe(k) * volti(k,iw) &
                       + omega2 * rho(k,iw) * vye(k,iw)
     vmyet_short(k,iw) = vmyet_short(k,iw) + hdiff_vye(k) * volti(k,iw) &
                       - omega2 * rho(k,iw) * vxe(k,iw)
     vmzet_short(k,iw) = vmzet_short(k,iw) + hdiff_vze(k) * volti(k,iw)

     ! Explicit momentum tendency from advective transport
     ! (weighted by T cell volume)

     vmxet_volt(k,iw) = hflux_vxe(k) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet_volt(k,iw) = hflux_vye(k) + vflux_vye(k-1) - vflux_vye(k)
     vmzet_volt(k,iw) = hflux_vze(k) + vflux_vze(k-1) - vflux_vze(k)

  enddo

  ! RAYLEIGH FRICTION ON THIL - only for horizontally homogeneous initialization

  if (dorayf .and. initial == 1) then

     do k = krayf_bot, mza
        thilt_short(k,iw) = thilt_short(k,iw) &
                          + rayf_cof(k) * dn01d(k) * (th01d(k) - thil(k,iw))
     enddo

     ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions
     ! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
     ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
     ! rayfx = 0.   ! Default: no extra RAYF
     ! do k = ka, mza
     !    thilt(k,iw) = thilt(k,iw) + max(rayf_cof(k),rayfx)  &
     !                * dn01d(k) * (th01d(k) - thil(k,iw))
     ! enddo

  endif ! (dorayf .and. initial == 1)

  ! Loop over T levels

  do k = ka, mza

     ! Explicit density tendency

     delex_rho(k) = dts * volti(k,iw) * (hflux_rho(k) + wmarw(k-1) - wmarw(k))

     ! Explicit density-thil tendency

     delex_rhothil(k) = dts * ( thilt_short(k,iw) &
                              + volti(k,iw) * (hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)) )

     ! RHOTHIL(t) and PRESS(t)

     rhothil(k) = rho(k,iw) * thil(k,iw)
     press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

  enddo

! NUDGING TENDENCY ON RHO

  if (nudflag == 1) then
     do k = ka, mza
        delex_rho(k) = delex_rho(k) + dts * rhot_nud(k,iw)
     enddo
  endif

  c6  = dts * .50 * fw
  c8  = dts * pc2
  c9  =-dts * fr
  c10 = dts * fw

  ! Rayleigh friction on W

  b20      (:) = 1.0
  wmt_other(:) = 0.0

  if (dorayfw) then

     do k = krayfw_bot, mza-1
        wmt_other(k) = - rayf_cofw(k) * wmc(k,iw)
        b20      (k) = b20(k) + c10 * rayf_cofw(k)
     enddo

     ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions

     ! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
     ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-1)
     ! rayfx = 0.   ! Default: no extra RAYF
     ! do k = ka, mza-1
     !    wmt_rayf(k) = - max(rayf_cofw(k),rayfx) * wmc(k,iw)
     ! enddo

  endif ! (dorayfw)

  ! Vertical loop over W levels

  do k = ka,mza-1

     ! Explicit vertical momentum tendency

     delex_wm(k) = dts * ( wmt_other(k) &

          + dzim(k) * real( (press_t(k) - press_t(k+1)) * rd_rt_w(k,iw) &

          - gdz_belo8(k) * rho(k,iw) - gdz_abov8(k) * rho(k+1,iw) ) &

          ! Average Coriolis force in vertical between T levels

          + wnxo2(iw) * (vmxet_short(k,iw) + vmxet_short(k+1,iw)) &
          + wnyo2(iw) * (vmyet_short(k,iw) + vmyet_short(k+1,iw)) &
          + wnzo2(iw) * (vmzet_short(k,iw) + vmzet_short(k+1,iw)) &

          ! ADD volume-weighted arrays VMXET_VOLT, VMYET_VOLT, VMZET_VOLT
          ! over both T levels and divide by combined VOLT of those T levels

          + ( wnx(iw) * (vmxet_volt(k,iw) + vmxet_volt(k+1,iw)) &
            + wny(iw) * (vmyet_volt(k,iw) + vmyet_volt(k+1,iw)) &
            + wnz(iw) * (vmzet_volt(k,iw) + vmzet_volt(k+1,iw)) &
            ) * volwi(k,iw) )

  enddo

  ! Fill matrix coefficients for implicit update of WM

  ! Loop over W pts
  do k = ka, mza-1
     b2(k)  = thil(k,iw) + thil(k+1,iw)
  enddo

  b2(ka-1)= 2. * thil(ka,iw)
  b2(mza) = 2. * thil(mza,iw)

  ! Loop over T pts
  do k = ka, mza
     b1(k)  = wc(k,iw) + wc(k-1,iw)
     b5(k)  = press_t(k) / rhothil(k)
     b6(k)  = c6  * volti(k,iw)
     b10(k) = c10 * volti(k,iw)
  enddo

  ! Loop over W pts
  do k = ka, mza-1

     b7(k)  = c6 * volwi(k,iw)
     b8(k)  = c8 * dzim (k) * rd_rt_w(k,iw)
     b9(k)  = c9 * gdzim(k)

     b11(k) = b8(k)  * b5(k)
     b12(k) = b8(k)  * b5(k+1)

     b13(k) = b9(k)  * dzt_top(k)
     b14(k) = b9(k)  * dzt_bot(k+1)

     b21(k) = b7(k)  * b1(k)
     b22(k) = b7(k)  * b1(k+1)

     b23(k) = b11(k) * b6(k)
     b24(k) = b12(k) * b6(k+1)

     b25(k) = b13(k) * b10(k)
     b26(k) = b14(k) * b10(k+1)

     b32(k) = b20(k) + arw(k,iw) &
            * (b22(k) - b21(k) + b2(k) * (b23(k) + b24(k)) + b25(k) - b26(k))

     b31(k) = - arw(k-1,iw) * (b21(k) + b23(k) * b2(k-1) + b25(k))

     b33(k) =   arw(k+1,iw) * (b22(k) - b24(k) * b2(k+1) + b26(k))

     b34(k) = delex_wm(k) &
            + b11(k) * delex_rhothil(k) - b12(k) * delex_rhothil(k+1) &
            + b13(k) * delex_rho(k)     + b14(k) * delex_rho(k+1)

  enddo

  ! Solve implicit tri-diagonal matrix equation for delta WM (del_wm)

  if (ka <= mza-1) then
     call tridiffo(mza,ka,mza-1,b31,b32,b33,b34,del_wm)
  endif

  ! Vertical loop over W points

  do k = ka,mza-1

     ! Change in vertical momentum from t to t+fw

     fwdel_wm(k) = fw * del_wm(k)

     ! Add vertical momentum at (t+fw) to array for long-timestep scalar transport

     wmsc(k,iw) = wmsc(k,iw) + wmc(k,iw) + fwdel_wm(k)

     ! Change in vertical fluxes from t to t + fw using del_wm value

     del_wmarw(k)  = fwdel_wm(k) * arw(k,iw)

     vflux_thil(k) = del_wmarw(k) * thil_upw(k,iw)
     vflux_vxe(k)  = del_wmarw(k) * vxe_upw(k,iw)
     vflux_vye(k)  = del_wmarw(k) * vye_upw(k,iw)
     vflux_vze(k)  = del_wmarw(k) * vze_upw(k,iw)

  enddo

  ! Set mass-flux-change array to 0 at top & bottom (other flux arrays already 0)

  del_wmarw(ka-1) = 0.
  del_wmarw(mza)  = 0.

  ! Vertical loop over T points

  do k = ka,mza

     ! Change of rho from (t) to (t+1)

     rho(k,iw) = rho(k,iw) + delex_rho(k) &
               + dts * volti(k,iw) * (del_wmarw(k-1) - del_wmarw(k))

     ! Change of rhothil from (t) to (t+1)

     del_rhothil = delex_rhothil(k) &
                 + dts * volti(k,iw) * (vflux_thil(k-1) - vflux_thil(k))

     ! Update pressure from (t) to (t+tp)

     press(k,iw) = press_t(k) + pc2 * b5(k) * del_rhothil

     ! Update thil from (t) to (t+1)

     thil(k,iw) = (rhothil(k) + del_rhothil) / rho(k,iw)

     ! Update volume-weighted momentum tendencies at T due to implicit flux correction

     vmxet_volt(k,iw) = vmxet_volt(k,iw) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet_volt(k,iw) = vmyet_volt(k,iw) + vflux_vye(k-1) - vflux_vye(k)
     vmzet_volt(k,iw) = vmzet_volt(k,iw) + vflux_vze(k-1) - vflux_vze(k)

  enddo

  if (lve2(iw) > 0) then

     do ksw = 1, lve2(iw)
        k = ksw + ka - 1

        ! Zero out vxe2, vye2, vze2 prior to new diagnosis
        vxe2(ksw,iw) = 0.
        vye2(ksw,iw) = 0.
        vze2(ksw,iw) = 0.

        ! Estimate velocity in T cells at (t+1) by prognostic method
        vxe1(k) = vxe(k,iw) + dts * (vmxet_short(k,iw) + vmxet_volt(k,iw) * volti(k,iw)) / real(rho(k,iw))
        vye1(k) = vye(k,iw) + dts * (vmyet_short(k,iw) + vmyet_volt(k,iw) * volti(k,iw)) / real(rho(k,iw))
        vze1(k) = vze(k,iw) + dts * (vmzet_short(k,iw) + vmzet_volt(k,iw) * volti(k,iw)) / real(rho(k,iw))
     enddo

     ! Loop over adjacent V faces

     do jv = 1, npoly
        iv = itab_w(iw)%iv(jv)

        ! Project full-forward-time vxe1, vye1, vze1 onto V faces that are
        ! below ground, and then project back to vxe2, vye2, vze2

        if (lpv(iv) > ka) then
           do k = ka, lpv(iv) - 1
              ksw = k - ka + 1
              vmt1 = vnx(iv) * vxe1(k) + vny(iv) * vye1(k) + vnz(iv) * vze1(k)

              vxe2(ksw,iw) = vxe2(ksw,iw) + itab_w(iw)%ecvec_vx(jv) * vmt1
              vye2(ksw,iw) = vye2(ksw,iw) + itab_w(iw)%ecvec_vy(jv) * vmt1
              vze2(ksw,iw) = vze2(ksw,iw) + itab_w(iw)%ecvec_vz(jv) * vmt1
           enddo
        endif

     enddo
  endif

  ! For shallow water test cases 2 & 5, rho & press are
  ! interpreted as water depth & height

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     topo_swtc = 0.

     if (nl%test_case == 5) then
        rad0_swtc = pi1 / 9.

        rad_swtc = sqrt((glonw(iw) * pio180 + 0.5 * pi1)**2 &
                 + (glatw(iw) * pio180 - pi1 / 6.) ** 2)

        topo_swtc = max(0., 2000. * (1. - rad_swtc / rad0_swtc))
     endif

     do k = ka,mza
        press(k,iw) = rho(k,iw) + topo_swtc + fp * delex_rho(k)
     enddo
  endif

  ! Vertical loop over W points

  do k = ka,mza-1

     ! Update WMC and WC to future value (at t+1) due to change in WM (del_wm)

     wmc(k,iw) = wmc(k,iw) + del_wm(k)
     wc (k,iw) = wmc(k,iw) / real(zwgt_bot8(k) * rho(k,iw) + zwgt_top8(k) * rho(k+1,iw))

  enddo

  ! Set top & bottom values of WC

  wc(ka-1,iw) = wc(ka,iw)

  wc(1:ka-2,iw) = 0.
  wc(mza   ,iw) = 0.

  ! Bottom boundary for THIL

  do k = 1, ka-1
     thil(k,iw) = thil(ka,iw)
  enddo

end subroutine prog_wrt_begs

!==========================================================================

subroutine comp_alpha_press(mrl)

  use mem_grid,    only: lpw, mza, zwgt_bot, zwgt_top
  use mem_ijtabs,  only: jtab_w, jtw_prog
  use consts_coms, only: pc1, rdry, rvap, cpocv
  use mem_basic,   only: rr_v, rr_w, theta, thil

  implicit none

  integer, intent(in)  :: mrl
  integer              :: j, iw, k

  !$omp parallel do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)


     ! Evaluate alpha coefficient for pressure
     do k = lpw(iw), mza
        alpha_press(k,iw) = pc1 * ( (rdry + rvap * rr_v(k,iw)) &
                                  * theta(k,iw) / thil(k,iw) ) ** cpocv
     enddo

     do k = lpw(iw), mza-1
        rd_rt_w(k,iw) = real( zwgt_bot(k) / (1. + rr_w(k  ,iw))     &
                            + zwgt_top(k) / (1. + rr_w(k+1,iw)), r8 )
     enddo

  enddo
  !$omp end parallel do

end subroutine comp_alpha_press

!============================================================================

subroutine prog_v_begs( iv, vmxet_volt, vmyet_volt, vmzet_volt,   &
                            vmxet_short, vmyet_short, vmzet_short )

  use mem_ijtabs,  only: itab_v
  use mem_basic,   only: vc, press, vmc, rho, vxe, vye, vze, rr_w
  use mem_tend,    only: vmt
  use misc_coms,   only: dtsm
  use consts_coms, only: eradi, gravo2
  use mem_grid,    only: mza, mwa, lpv, volvi, xev, yev, zev, &
                         unx, uny, unz, vnx, vny, vnz, vnxo2, vnyo2, vnzo2, &
                         dniu, dniv
  use oname_coms,  only: nl

  implicit none

  integer, intent(in) :: iv

  real, intent(in) :: vmxet_volt (mza,mwa)
  real, intent(in) :: vmyet_volt (mza,mwa)
  real, intent(in) :: vmzet_volt (mza,mwa)
  real, intent(in) :: vmxet_short(mza,mwa)
  real, intent(in) :: vmyet_short(mza,mwa)
  real, intent(in) :: vmzet_short(mza,mwa)

  integer :: k, kb
  integer :: iw1,iw2,im1,im2

  real :: dts
  real :: vx, vy, vz, uc, watv, tke1, tke2, vortp_v, dtso2dnu, dtso2dnv
  real :: pgf(mza)

  ! Extract neighbor indices and coefficients for this point in the U stencil

  iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
  im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)

  dts = dtsm(itab_v(iv)%mrlv)

  kb = lpv(iv)

  ! Vertical loop over V points to compute pressure gradient force

  do k = kb, mza
     pgf(k) = dniv(iv) * real(press(k,iw1) - press(k,iw2)) &
            * (1. - 0.5 * (rr_w(k,iw1) + rr_w(k,iw2)))
  enddo

  ! For shallow water test cases 2 & 5, rho & press are
  ! interpreted as water depth & height

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     do k = kb,mza
        pgf(k) = pgf(k) * gravo2 * (rho(k,iw1) + rho(k,iw2))
     enddo
  endif

  ! Update VMC

  if (.not. rotational) then

     do k = kb,mza
        vmc(k,iv) = vmc(k,iv) + dts * ( vmt(k,iv) + pgf(k)                   &

                  + vnxo2(iv) * (vmxet_short(k,iw1) + vmxet_short(k,iw2))    &
                  + vnyo2(iv) * (vmyet_short(k,iw1) + vmyet_short(k,iw2))    &
                  + vnzo2(iv) * (vmzet_short(k,iw1) + vmzet_short(k,iw2))    &

                  + ( vnx(iv) * (vmxet_volt(k,iw1) + vmxet_volt(k,iw2))      &
                    + vny(iv) * (vmyet_volt(k,iw1) + vmyet_volt(k,iw2))      &
                    + vnz(iv) * (vmzet_volt(k,iw1) + vmzet_volt(k,iw2)) )    &
                  * volvi(k,iv) )

        vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )
     enddo

  else  ! If using rotational form:

     dtso2dnu = 0.5 * dts * dniu(iv)
     dtso2dnv = 0.5 * dts * dniv(iv)

     ! Vertical loop over V levels

     do k = kb,mza
        vx = .5 * (vxe(k,iw1) + vxe(k,iw2))
        vy = .5 * (vye(k,iw1) + vye(k,iw2))
        vz = .5 * (vze(k,iw1) + vze(k,iw2))

        uc =    unx(iv) * vx + uny(iv) * vy + unz(iv) * vz
        watv = (xev(iv) * vx + yev(iv) * vy + zev(iv) * vz) * eradi

        tke1 = .5 * (vxe(k,iw1)**2 + vye(k,iw1)**2 + vze(k,iw1)**2)
        tke2 = .5 * (vxe(k,iw2)**2 + vye(k,iw2)**2 + vze(k,iw2)**2)

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@A1
        vortp_v = .5 * (vortp(k,im1) + vortp(k,im2)) &

             ! APVM method (Ringler et al. 2010) to obtain upwinded vortp at V

             + dtso2dnu * uc       * (vortp(k,im1) - vortp(k,im2)) &
             + dtso2dnv * vc(k,iv) * (vortp_t(k,iw1) - vortp_t(k,iw2))
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@A2

        vmc(k,iv) = vmc(k,iv) + dts * (vmt(k,iv) + pgf(k) &

                  + vnxo2(iv) * (vmxet_short(k,iw1) + vmxet_short(k,iw2))    &
                  + vnyo2(iv) * (vmyet_short(k,iw1) + vmyet_short(k,iw2))    &
                  + vnzo2(iv) * (vmzet_short(k,iw1) + vmzet_short(k,iw2))    &

                  + .5 * (rho(k,iw1) + rho(k,iw2))                           &
                  * (dniv(iv) * (tke1 - tke2)                                &
                  + uc * vortp_v                                             &
                  - watv * .5 * (vortn(k-1,iv) + vortn(k,iv))))

        vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )
     enddo

  endif ! rotational

end subroutine prog_v_begs

!===============================================================================

subroutine prep_rotational(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtab_m, itab_m, &
                          jtm_vadj, jtv_prog, jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp
  use mem_basic,    only: vc, wc
  use mem_grid,     only: mza, lpm, lpv, lpw, dzim, zfact, &
                          zfacit, zfacim, dnv, dniv, arm0
  use obnd,         only: lbcopy_m, lbcopy_w
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m

  implicit none

  integer, intent(in) :: mrl

  integer :: j, iv, k, kb
  integer :: iw, iw1, iw2
  integer :: im,npoly,jv
  integer :: jm
  real    :: arm0i

! IF USING ROTATIONAL METHOD, HORIZONTAL LOOP OVER M/P COLUMNS
! FOR COMPUTING VERTICAL VORTICITY.

  !$omp parallel do private(im,npoly,kb,jv,iv,k,arm0i)
  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)

     npoly = itab_m(im)%npoly
     kb    = lpm(im)

     vortp(:,im) = 0.

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1,npoly
        iv = itab_m(im)%iv(jv)

        if (itab_v(iv)%im(2) == im) then

           do k = kb,mza
              vortp(k,im) = vortp(k,im) + vc(k,iv) * dnv(iv)
           enddo

        else

           do k = kb,mza
              vortp(k,im) = vortp(k,im) - vc(k,iv) * dnv(iv)
           enddo

        endif
     enddo

     ! Convert circulation to relative vertical vorticity at M 
     ! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
     ! divide their quotient by zfact)

     arm0i = 1. / arm0(im)

     do k = kb,mza
        vortp(k,im) = vortp(k,im) * arm0i * zfacit(k)
     enddo

  enddo
  !$omp end parallel do 

  ! PARALLEL SEND/RECV OF VORTP

  if (iparallel == 1) then
     call mpi_send_m(mrl, rvara1=vortp)
     call mpi_recv_m(mrl, rvara1=vortp)
  endif
  call lbcopy_m(mrl, a1=vortp)

! IF USING ROTATIONAL FORM FOR PROGNOSING HORIZONTAL MOMENTUM,
! EVALUATE VORTP_T AND VORTN AND CORIOLIS FORCING TERMS

  ! Horizontal loop over W/T columns
  !$omp parallel do private(iw,npoly,kb,k,jm,im)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     npoly = itab_w(iw)%npoly
     kb = lpw(iw)

     vortp_t(:,iw) = 0.

     ! Loop over M neighbors of W
     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        ! Vertical loop over T levels; average vorticity from P points to T point.
        do k = kb,mza
           vortp_t(k,iw) = vortp_t(k,iw) + itab_w(iw)%farm(jm) * vortp(k,im)
        enddo
     enddo

  enddo
  !$omp end parallel do

  if (iparallel == 1) call mpi_send_w(mrl, rvara1=vortp_t)
  if (iparallel == 1) call mpi_recv_w(mrl, rvara1=vortp_t)
  call lbcopy_w(mrl, a1=vortp_t)

  ! Horizontal loop over V/N columns
  !$omp parallel do private(iv,iw1,iw2,k,kb) 
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
     kb = lpv(iv)

     ! Vertical loop over N levels; compute horizontal relative vorticity
     ! and vertical velocity at N points (at time T)
     do k = kb,mza-1
        vortn(k,iv) = zfacim(k) * ((wc(k,iw1) - wc(k,iw2)) * dniv(iv) &
                    + (vc(k+1,iv) * zfact(k+1) - vc(k,iv) * zfact(k)) * dzim(k))
     enddo

     vortn(1:kb-1,iv) = 0.
     vortn(mza,iv) = 0.

  enddo
  !$omp end parallel do

end subroutine prep_rotational


end module wrtv_mem
