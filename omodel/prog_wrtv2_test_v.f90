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

  real, allocatable :: alpha_press(:,:)
  real, allocatable :: pwfac      (:,:)
  real, allocatable :: pvfac      (:,:)

  real, allocatable :: b20(:)

! integer, parameter :: nrk = 2
  integer, parameter :: nrk = 3

  real, parameter :: xstg3(3) = (/ 1./3., 0.5, 1.0 /)
  real, parameter :: xstg2(2) =        (/ 0.5, 1.0 /)

  private
  public alloc_wrtv_mem, prog_wrtv, comp_alpha_press

contains

!===============================================================================

subroutine alloc_wrtv_mem()

  use mem_grid,    only: mza, mwa, mva
  use misc_coms,   only: rinit
  use consts_coms, only: pc1, rdry, cpocv

  implicit none

  real, parameter :: apd = pc1 * rdry ** cpocv

  allocate(alpha_press(mza,mwa)) ; alpha_press = apd
  allocate(pwfac      (mza,mwa)) ; pwfac       = rinit
  allocate(pvfac      (mza,mva)) ; pvfac       = rinit

  allocate(b20(mza)) ; b20 = 1.0

end subroutine alloc_wrtv_mem

!===============================================================================

subroutine prog_wrtv(vmsca, wmsca)

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

! The above description needs amending: The time differencing method is RK3.

  use mem_ijtabs,   only: jtab_v, jtab_w, istp, mrl_endl, &
                          mrl_ends, jtw_wadj, jtm_vadj, jtv_prog, &
                          jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp, &
                          itab_w, itab_v
  use mem_basic,    only: rho, thil, wc, press, vmc, vc, wmc, &
                          vxe, vye, vze, theta, tair
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, arv, dzto2, dztsqo6, &
                          volt, volti
  use mem_tend,     only: thilt, vmxet, vmyet, vmzet
  use misc_coms,    only: iparallel, dtsm, mstp, time8p, dtlm
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use oplot_coms,   only: op
  use obnd,         only: lbcopy_m, lbcopy_w, lbcopy_v
  use oname_coms,   only: nl
  use mem_adv,      only: gxps_scp, gyps_scp, gxyps_scp, gxxps_scp, gyyps_scp, &
                          xx0_v, yy0_v, xy0_v
  use grad_lib,     only: grad_z_quad, grad_t2d, grad_t2d_quad
  use vel_t3d,      only: diagvel_t3d
  use consts_coms,  only: p00i, rocp
  use mem_turb,     only: akmodx, akhodx, khtopv, kmtopv, khtop
  use pbl_drivers,  only: solve_eddy_diff_heat, solve_eddy_diff_vc, solve_eddy_diff_vxe
  use mem_rayf,     only: dorayfw, rayf_cofw, krayfw_bot
  use oname_coms,   only: nl

  implicit none

  real, intent(inout) :: vmsca(mza,mva)
  real, intent(inout) :: wmsca(mza,mwa)

  integer :: j, iv, iw, k, mrl, iw1, iw2, iwn, jv, istage
  real    :: dts, v4

! automatic arrays

  real :: vmca    (mza,mva)

  real :: thil_upv(mza,mva) ! Upstreamed THIL at each V interface
  real :: vxe_upv (mza,mva) ! Upstreamed VXE  at each V interface
  real :: vye_upv (mza,mva) ! Upstreamed VYE  at each V interface
  real :: vze_upv (mza,mva) ! Upstreamed VZE  at each V interface

  real :: thil_upw(mza,mwa) ! Upstreamed THIL at each W level
  real :: vxe_upw (mza,mwa) ! Upstreamed VXE  at each W level
  real :: vye_upw (mza,mwa) ! Upstreamed VYE  at each W level
  real :: vze_upw (mza,mwa) ! Upstreamed VZE  at each W level

  real :: vmxet_rk(mza,mwa)
  real :: vmyet_rk(mza,mwa)
  real :: vmzet_rk(mza,mwa)

  real :: rth0(mza,mwa)
  real :: wmc0(mza,mwa)
  real :: rho0(mza,mwa)
  real :: vmc0(mza,mva)

  real :: gzps_th (mza), gzps_vx (mza), gzps_vy (mza), gzps_vz (mza)
  real :: gzzps_th(mza), gzzps_vx(mza), gzzps_vy(mza), gzzps_vz(mza)

  real :: thilt_short(mza,mwa)
  real :: vmxet_short(mza,mwa)
  real :: vmyet_short(mza,mwa)
  real :: vmzet_short(mza,mwa)
  real :: vmt_short  (mza,mva)

  real :: hflux_thil(mza)
  real :: hflux_vxe (mza)
  real :: hflux_vye (mza)
  real :: hflux_vze (mza)

  real :: gxps_vxe(mza,mwa)
  real :: gyps_vxe(mza,mwa)
  real :: gxps_vye(mza,mwa)
  real :: gyps_vye(mza,mwa)
  real :: gxps_vze(mza,mwa)
  real :: gyps_vze(mza,mwa)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  real :: gxxps_vxe(mza,mwa)
!!  real :: gyyps_vxe(mza,mwa)
!!  real :: gxyps_vxe(mza,mwa)
!!
!!  real :: gxxps_vye(mza,mwa)
!!  real :: gyyps_vye(mza,mwa)
!!  real :: gxyps_vye(mza,mwa)
!!
!!  real :: gxxps_vze(mza,mwa)
!!  real :: gyyps_vze(mza,mwa)
!!  real :: gxyps_vze(mza,mwa)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS ROUTINE IS PERFORMED AT THE END OF EACH SMALL TIMESTEP

  mrl = mrl_ends(istp)
  if (mrl == 0) return

  !$omp parallel
  !$omp do private(iw,k,v4,jv,iv,iwn)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     thil_upw(      mza,iw) = 0.
     thil_upw(lpw(iw)-1,iw) = 0.

     do k = lpw(iw), mza
        v4 = real(volt(k,iw))

        ! Store starting values for R-K timestepping
        rth0(k,iw) = thil(k,iw) * rho(k,iw)
        rho0(k,iw) = rho (k,iw)
        wmc0(k,iw) = wmc (k,iw)

        ! Include long-timestep tendencies in each short timestep
        ! and weight by volume
        thilt_short(k,iw) = thilt(k,iw) * v4
        vmxet_short(k,iw) = vmxet(k,iw) * v4
        vmyet_short(k,iw) = vmyet(k,iw) * v4
        vmzet_short(k,iw) = vmzet(k,iw) * v4
     enddo

     if (khtop(iw) >= lpw(iw)) then
        call solve_eddy_diff_heat(iw, thilt_short(:,iw))
     endif
     call solve_eddy_diff_vxe (iw, vmxet_short(:,iw), vmyet_short(:,iw), vmzet_short(:,iw))

     do jv = 1, itab_w(iw)%npoly
        iv  = itab_w(iw)%iv(jv)
        iwn = itab_w(iw)%iw(jv)

        ! Loop over T levels
        do k = lpv(iv), khtopv(iv)
           thilt_short(k,iw) = thilt_short(k,iw) + akhodx(k,iv) * (thil(k,iwn) - thil(k,iw))
        enddo
        do k = lpv(iv), kmtopv(iv)
           vmxet_short(k,iw) = vmxet_short(k,iw) + akmodx(k,iv) * (vxe(k,iwn) - vxe(k,iw))
           vmyet_short(k,iw) = vmyet_short(k,iw) + akmodx(k,iv) * (vye(k,iwn) - vye(k,iw))
           vmzet_short(k,iw) = vmzet_short(k,iw) + akmodx(k,iv) * (vze(k,iwn) - vze(k,iw))
        enddo
     enddo

  enddo
  !$omp end do nowait

  !$omp do private(iv,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     vmc0(:,iv) = vmc(:,iv)

!     do k = lpv(iv), mza
!        vmt_short(k,iv) = vmt_short(k,iv) + vmt(k,iv)
!     enddo

     vmt_short(:,iv) = 0.

!    call solve_eddy_diff_vc(iv, vmt_short(:,iv))

     call extra_mix_top_vc(iv, vmt_short(:,iv))
  enddo
  !$omp end do
  !$omp end parallel

  call divh_damp(1, vmt_short)

  do istage = 1, nrk

     if (nrk == 2) then
        dts = xstg2(istage) * dtsm(1)
     else
        dts = xstg3(istage) * dtsm(1)
     endif

     ! Set coefficients for vertical velocity damping (Rayleigh friction)

     if     (mstp == 0 .and. time8p <         dtlm(1)) then
        b20 = 9.0
     elseif (mstp == 1 .and. time8p < 2._r8 * dtlm(1)) then
        b20 = 3.0
     elseif (mstp == 2 .and. time8p < 3._r8 * dtlm(1)) then
        b20 = 2.0
     elseif (mstp == 3 .and. time8p < 4._r8 * dtlm(1)) then
        b20 = 1.5
     else
        b20 = 1.0
     endif

     if (dorayfw) then
        do k = krayfw_bot, mza-1
           b20(k) = max(b20(k), 1.0 + dts * rayf_cofw(k))
        enddo
     endif

! INCLUDE THE LONG TIMESTEP THIL AND MOMENTUM TENDENCIES IN EACH SHORT TIMESTEP

!  ! omp parallel
!  ! omp do private(iw,dt,k)
!  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

!     dt = dtsm(itab_w(iw)%mrlw)

!!     do k = lpw(iw), mza
!!        thilt_short(k,iw) = thilt(k,iw)
!!        vmxet_short(k,iw) = vmxet(k,iw)
!!        vmyet_short(k,iw) = vmyet(k,iw)
!!        vmzet_short(k,iw) = vmzet(k,iw)
!!     enddo

!!! temp
!     do k = lpw(iw), mza-1
!        wmca(k,iw) = wmc(k,iw) * arw(k,iw)
!     enddo

!     call donorpointw(iw, 0.0, wc(:,iw), vxe(:,iw), vye(:,iw), vze(:,iw), kdepw(:,iw))

!  enddo
!  ! omp end do nowait

  ! Compute CFL numbers needed by Thuburn monotonic scheme

!!!!!! temp
!!
!  if (nl%ithil_monot > 0 .and. istage == 3) then
!     call comp_cfls_short( mrl, vmca, wmca )
!  endif

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


  !$omp parallel
  !$omp do private(iv,k)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     do k = lpv(iv), mza
        vmca(k,iv) = vmc(k,iv) * arv(k,iv)
        if (istage == nrk) vmsca(k,iv) = vmsca(k,iv) + vmca(k,iv)
     enddo

!   call wind_vec_at_v(iv, vxe_upv(:,iv), vye_upv(:,iv), vze_upv(:,iv))

  enddo
  !$omp end do

  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        vmxet_rk(k,iw) = vmxet_short(k,iw)
        vmyet_rk(k,iw) = vmyet_short(k,iw)
        vmzet_rk(k,iw) = vmzet_short(k,iw)
     enddo

     if (nl%horiz_adv_order <= 2) then

        call grad_t2d(iw, thil, gxps_scp(:,iw),  gyps_scp(:,iw))

        if (nl%iscal_monot > 0) then
           call limit_h2(iw, thil, gxps_scp(:,iw), gyps_scp(:,iw))
        endif

     else

        call grad_t2d_quad(iw, thil, gxps_scp (:,iw), gyps_scp (:,iw), &
                                     gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))

        if (nl%iscal_monot > 0) then
           call limit_h3(iw, thil, gxps_scp(:,iw), gyps_scp(:,iw),&
                                   gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
        endif

     endif

     call grad_t2d(iw, vxe, gxps_vxe(:,iw), gyps_vxe(:,iw))
     call grad_t2d(iw, vye, gxps_vye(:,iw), gyps_vye(:,iw))
     call grad_t2d(iw, vze, gxps_vze(:,iw), gyps_vze(:,iw))

!!     call grad_t2d_quad(iw, vxe, gxps_vxe (:,iw), gyps_vxe (:,iw), &
!!                                 gxxps_vxe(:,iw), gxyps_vxe(:,iw), gyyps_vxe(:,iw))
!!
!!     call grad_t2d_quad(iw, vye, gxps_vye (:,iw), gyps_vye(:,iw), &
!!                                 gxxps_vye(:,iw), gxyps_vye(:,iw), gyyps_vye(:,iw))
!!
!!     call grad_t2d_quad(iw, vze, gxps_vze (:,iw), gyps_vze(:,iw), &
!!                                 gxxps_vze(:,iw), gxyps_vze(:,iw), gyyps_vze(:,iw))


!!     call grad_t2d_v(iw, vxe_upv, vye_upv, vze_upv,  &
!!                     gxps_vxe(:,iw), gyps_vxe(:,iw), &
!!                     gxps_vye(:,iw), gyps_vye(:,iw), &
!!                     gxps_vze(:,iw), gyps_vze(:,iw)  )

  enddo
  !$omp end do nowait
  !$omp end parallel

! MPI send of THIL, VXE, VYE, and VZE gradient components

  if (iparallel == 1) then

     if (nl%horiz_adv_order <= 2) then
        call mpi_send_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
                             rvara3=gxps_vxe, rvara4=gyps_vxe, &
                             rvara5=gxps_vye, rvara6=gyps_vye, &
                             rvara7=gxps_vze, rvara8=gyps_vze  )
     else
        call mpi_send_w(mrl, rvara1=gxps_scp,   rvara2=gyps_scp, &
                             rvara3=gxxps_scp,  rvara4=gxyps_scp, rvara5=gyyps_scp,  &
                             rvara6=gxps_vxe,   rvara7=gyps_vxe, &
                             rvara8=gxps_vye,   rvara9=gyps_vye, &
                             rvara10=gxps_vze,  rvara11=gyps_vze )
     endif
!!
!!        call mpi_send_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
!!                             rvara3=gxps_vxe, rvara4=gyps_vxe, &
!!                             rvara5=gxps_vye, rvara6=gyps_vye, &
!!                             rvara7=gxps_vze, rvara8=gyps_vze, &
!!                             rvara9=gxxps_vxe, rvara10=gyyps_vxe, rvara11=gxyps_vxe, &
!!                             rvara12=gxxps_vye, rvara13=gyyps_vye,  rvara14=gxyps_vye,&
!!                             rvara15=gxxps_vze, rvara16=gyyps_vze,  rvara17=gxyps_vze )
  endif

  !$omp parallel private(gzps_th,gzps_vx,gzps_vy,gzps_vz,   &
  !$                     gzzps_th,gzzps_vx,gzzps_vy,gzzps_vz)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Evaluate vertical gradient of thil, vxe, vye, and vze for begs

     call grad_z_quad(iw, thil(:,iw), gzps_th, gzzps_th)
     call grad_z_quad(iw, vxe (:,iw), gzps_vx, gzzps_vx)
     call grad_z_quad(iw, vye (:,iw), gzps_vy, gzzps_vy)
     call grad_z_quad(iw, vze (:,iw), gzps_vz, gzzps_vz)

     if (nl%iscal_monot > 0) then
        call limit_z(iw, thil(:,iw), gzps_th, gzzps_th)
     endif

    ! Evaluate upwinded thil, vxe, vye, and vze at each w face

    do k = lpw(iw), mza-1

        if (wmc(k,iw) > 0.) then

           thil_upw(k,iw) = thil(k,iw)  &
                          +   dzto2(k) *  gzps_th(k) &
                          + dztsqo6(k) * gzzps_th(k)

           vxe_upw(k,iw)  = vxe(k,iw)                &
                          +   dzto2(k) *  gzps_vx(k) &
                          + dztsqo6(k) * gzzps_vx(k)

           vye_upw(k,iw)  = vye(k,iw)                &
                          +   dzto2(k) *  gzps_vy(k) &
                          + dztsqo6(k) * gzzps_vy(k)

           vze_upw(k,iw)  = vze(k,iw)                &
                          +   dzto2(k) *  gzps_vz(k) &
                          + dztsqo6(k) * gzzps_vz(k)
        else

           thil_upw(k,iw) = thil(k+1,iw)                 &
                          -   dzto2(k+1) *  gzps_th(k+1) &
                          + dztsqo6(k+1) * gzzps_th(k+1)

           vxe_upw(k,iw) = vxe(k+1,iw)                  &
                         -   dzto2(k+1) *  gzps_vx(k+1) &
                         + dztsqo6(k+1) * gzzps_vx(k+1)

           vye_upw(k,iw) = vye(k+1,iw)                  &
                         -   dzto2(k+1) *  gzps_vy(k+1) &
                         + dztsqo6(k+1) * gzzps_vy(k+1)

           vze_upw(k,iw) = vze(k+1,iw)                  &
                         -   dzto2(k+1) *  gzps_vz(k+1) &
                         + dztsqo6(k+1) * gzzps_vz(k+1)
        endif

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
        call mpi_recv_w(mrl, rvara1=gxps_scp,   rvara2=gyps_scp, &
                             rvara3=gxxps_scp,  rvara4=gxyps_scp, rvara5=gyyps_scp, &
                             rvara6=gxps_vxe,   rvara7=gyps_vxe, &
                             rvara8=gxps_vye,   rvara9=gyps_vye, &
                             rvara10=gxps_vze,  rvara11=gyps_vze )
     endif
!!
!!        call mpi_recv_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
!!                             rvara3=gxps_vxe, rvara4=gyps_vxe, &
!!                             rvara5=gxps_vye, rvara6=gyps_vye, &
!!                             rvara7=gxps_vze, rvara8=gyps_vze, &
!!                             rvara9=gxxps_vxe, rvara10=gyyps_vxe, rvara11=gxyps_vxe, &
!!                             rvara12=gxxps_vye, rvara13=gyyps_vye,  rvara14=gxyps_vye,&
!!                             rvara15=gxxps_vze, rvara16=gyyps_vze,  rvara17=gxyps_vze )


  endif

! Lateral boundary copy of THIL, VXE, VYE, and VZE gradient components

  if (nl%horiz_adv_order <= 2) then
     call lbcopy_w(mrl, a1=gxps_scp, a2=gyps_scp, &
                        a3=gxps_vxe, a4=gyps_vxe, &
                        a5=gxps_vye, a6=gyps_vye, &
                        a7=gxps_vze, a8=gyps_vze  )
  else
     call lbcopy_w(mrl, a1=gxps_scp,  a2=gyps_scp, &
                        a3=gxxps_scp, a4=gxyps_scp,  a5=gyyps_scp,  &
                        a6=gxps_vxe,  a7=gyps_vxe, &
                        a8=gxps_vye,  a9=gyps_vye, &
                        a10=gxps_vze, a11=gyps_vze )
  endif

! EVALUATE UPWINDED THIL, VXE, VYE, AND VZE AT EACH V FACE FOR
! HORIZONTAL FLUX COMPUTATION

  !$omp parallel
  !$omp do private(iv,k,iw1,iw2)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza

        if (vmca(k,iv) > 0.) then

           if (nl%horiz_adv_order <= 2) then
              thil_upv(k,iv) = thil(k,iw1) &
                             + itab_v(iv)%dxps(1) * gxps_scp(k,iw1) &
                             + itab_v(iv)%dyps(1) * gyps_scp(k,iw1)
           else
              thil_upv(k,iv) = thil(k,iw1) &
                             + itab_v(iv)%dxps(1) *  gxps_scp(k,iw1) &
                             + itab_v(iv)%dyps(1) *  gyps_scp(k,iw1) &
                             + xx0_v(1,iv)        * gxxps_scp(k,iw1) &
                             + yy0_v(1,iv)        * gyyps_scp(k,iw1) &
                             + xy0_v(1,iv)        * gxyps_scp(k,iw1)
           endif

           vxe_upv(k,iv)  = vxe(k,iw1) &
                          + itab_v(iv)%dxps(1) * gxps_vxe(k,iw1) &
                          + itab_v(iv)%dyps(1) * gyps_vxe(k,iw1)

           vye_upv(k,iv)  = vye(k,iw1) &
                          + itab_v(iv)%dxps(1) * gxps_vye(k,iw1) &
                          + itab_v(iv)%dyps(1) * gyps_vye(k,iw1)

           vze_upv(k,iv)  = vze(k,iw1) &
                          + itab_v(iv)%dxps(1) * gxps_vze(k,iw1) &
                          + itab_v(iv)%dyps(1) * gyps_vze(k,iw1)
!!
!!           vxe_upv(k,iv)  = vxe(k,iw1) &
!!                          + itab_v(iv)%dxps(1) * gxps_vxe(k,iw1) &
!!                          + itab_v(iv)%dyps(1) * gyps_vxe(k,iw1) &
!!                          + xx0_v(1,iv)        * gxxps_vxe(k,iw1) &
!!                          + yy0_v(1,iv)        * gyyps_vxe(k,iw1) &
!!                          + xy0_v(1,iv)        * gxyps_vxe(k,iw1)
!!
!!           vye_upv(k,iv)  = vye(k,iw1) &
!!                          + itab_v(iv)%dxps(1) * gxps_vye(k,iw1) &
!!                          + itab_v(iv)%dyps(1) * gyps_vye(k,iw1) &
!!                          + xx0_v(1,iv)        * gxxps_vye(k,iw1) &
!!                          + yy0_v(1,iv)        * gyyps_vye(k,iw1) &
!!                          + xy0_v(1,iv)        * gxyps_vye(k,iw1)
!!
!!           vze_upv(k,iv)  = vze(k,iw1) &
!!                          + itab_v(iv)%dxps(1) * gxps_vze(k,iw1) &
!!                          + itab_v(iv)%dyps(1) * gyps_vze(k,iw1) &
!!                          + xx0_v(1,iv)        * gxxps_vze(k,iw1) &
!!                          + yy0_v(1,iv)        * gyyps_vze(k,iw1) &
!!                          + xy0_v(1,iv)        * gxyps_vze(k,iw1)

        else

           if (nl%horiz_adv_order <= 2) then
              thil_upv(k,iv) = thil(k,iw2) &
                             + itab_v(iv)%dxps(2) * gxps_scp(k,iw2) &
                             + itab_v(iv)%dyps(2) * gyps_scp(k,iw2)
           else
              thil_upv(k,iv) = thil(k,iw2) &
                             + itab_v(iv)%dxps(2) *  gxps_scp(k,iw2) &
                             + itab_v(iv)%dyps(2) *  gyps_scp(k,iw2) &
                             + xx0_v(2,iv)        * gxxps_scp(k,iw2) &
                             + yy0_v(2,iv)        * gyyps_scp(k,iw2) &
                             + xy0_v(2,iv)        * gxyps_scp(k,iw2)
           endif

           vxe_upv(k,iv)  = vxe(k,iw2) &
                          + itab_v(iv)%dxps(2) * gxps_vxe(k,iw2) &
                          + itab_v(iv)%dyps(2) * gyps_vxe(k,iw2)

           vye_upv(k,iv)  = vye(k,iw2) &
                          + itab_v(iv)%dxps(2) * gxps_vye(k,iw2) &
                          + itab_v(iv)%dyps(2) * gyps_vye(k,iw2)

           vze_upv(k,iv)  = vze(k,iw2) &
                          + itab_v(iv)%dxps(2) * gxps_vze(k,iw2) &
                          + itab_v(iv)%dyps(2) * gyps_vze(k,iw2)
!!
!!           vxe_upv(k,iv)  = vxe(k,iw2) &
!!                          + itab_v(iv)%dxps(2) * gxps_vxe(k,iw2) &
!!                          + itab_v(iv)%dyps(2) * gyps_vxe(k,iw2) &
!!                          + xx0_v(2,iv)        * gxxps_vxe(k,iw1) &
!!                          + yy0_v(2,iv)        * gyyps_vxe(k,iw1) &
!!                          + xy0_v(2,iv)        * gxyps_vxe(k,iw1)
!!
!!           vye_upv(k,iv)  = vye(k,iw2) &
!!                          + itab_v(iv)%dxps(2) * gxps_vye(k,iw2) &
!!                          + itab_v(iv)%dyps(2) * gyps_vye(k,iw2) &
!!                          + xx0_v(2,iv)        * gxxps_vye(k,iw1) &
!!                          + yy0_v(2,iv)        * gyyps_vye(k,iw1) &
!!                          + xy0_v(2,iv)        * gxyps_vye(k,iw1)
!!
!!           vze_upv(k,iv)  = vze(k,iw2) &
!!                          + itab_v(iv)%dxps(2) * gxps_vze(k,iw2) &
!!                          + itab_v(iv)%dyps(2) * gyps_vze(k,iw2) &
!!                          + xx0_v(2,iv)        * gxxps_vze(k,iw1) &
!!                          + yy0_v(2,iv)        * gyyps_vze(k,iw1) &
!!                          + xy0_v(2,iv)        * gxyps_vze(k,iw1)

        endif

     enddo

  enddo
  !$omp end do

! Compute the monotonic or positive-definite flux limiters and then apply them

!  if (nl%ithil_monot == 1 .and. istage==3) then
!     call comp_and_apply_monot_limits(mrl, thil, thil_upw, thil_upv, wmca, vmca)
!  endif

! MAIN LOOP OVER W COLUMNS FOR UPDATING WM, WC, RHO, THIL, AND PRESS

  !$omp do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Prognose vertical velocity, density, thil, and diagnose pressure
     call prog_wrt_begs( iw, istage, vmca, wmsca, dts,         &
                         thil_upv, vxe_upv, vye_upv, vze_upv,  &
                         thil_upw, vxe_upw, vye_upw, vze_upw,  &
                         thilt_short, vmxet_rk, vmyet_rk, vmzet_rk, &
                         rho0, rth0, wmc0 )

  enddo
  !$omp end do nowait
  !$omp end parallel


! MPI SEND/RECV and LBC copy of quantities needed for prog_v:
! PRESS, RHO, VMXET_SHORT, VMYET_SHORT, and VMZET_SHORT

  if (iparallel == 1) then

     call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_rk, rvara2=vmyet_rk, rvara3=vmzet_rk, &
                     rvara4=thil, rvara5=wmc, rvara6=wc)

     call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_rk, rvara2=vmyet_rk, rvara3=vmzet_rk, &
                     rvara4=thil, rvara5=wmc, rvara6=wc)
  endif

  call lbcopy_w(mrl, a1=vmxet_rk,  a2=vmyet_rk,  a3=vmzet_rk, &
                     a4=thil,      a5=wmc,       a6=wc,       &
                     d1=rho,       d2=press )

! enddo  ! istage temp!!!!!!!!!!!!
!!!!!!!!!!!!!!
! stop!!!!!!!!!!!!!!!!!!!!
!

  !!!!!!!!!!! temp !!!!!
  !!! theta = thil

! COMPUTE TERMS FOR ROTATIONAL FORM OF HORIZONTAL MOMENTUM

!!  if (rotational) then
!!     call prep_rotational(mrl)
!!  endif

!!  if (dorayfdiv) then
!!
!!     !$omp parallel do private(iw,jv,iv,k)
!!     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!!
!!        div2d(:,iw) = 0.0
!!        do jv = 1, itab_w(iw)%npoly
!!           iv = itab_w(iw)%iv(jv)
!!
!!           do k = krayfdiv_bot, mza
!!              div2d(k,iw) = div2d(k,iw) - itab_w(iw)%dirv(jv) * vmc(k,iv) * dnu(iv)
!!           enddo
!!        enddo
!!     enddo
!!     !$omp end parallel do
!!
!!     ! MPI send of div2d
!!     if (iparallel == 1) then
!!        call mpi_send_w(mrl, svara1=div2d)
!!        call mpi_recv_w(mrl, svara1=div2d)
!!     endif
!!     call lbcopy_w(mrl, s1=div2d)
!!
!!  endif ! dorayfdiv

! MAIN LOOP OVER V POINTS TO UPDATE VMC AND VC

  !$omp parallel do private(iv)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     call prog_v_begs( iv, dts, vmc0, vmxet_rk, vmyet_rk, vmzet_rk, vmt_short )

  enddo
  !$omp end parallel do

! MPI SEND/RECV and LBC of VC and VMC

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif
  call lbcopy_v(mrl, vmc=vmc, vc=vc)

  call diagvel_t3d(mrl)

  enddo ! istage loop

!!  mrl = mrl_endl(istp)
!!  if (mrl > 0) then
!!
!!     !$omp parallel do private(k)
!!     do iw = 2, mwa
!!        do k = lpw(iw), mza
!!!        theta(k,iw) = thil(k,iw)
!!           tair(k,iw) = theta(k,iw) * (real(press(k,iw)) * p00i) ** rocp
!!        enddo
!!        tair(1:lpw(iw)-1,iw) = tair(lpw(iw),iw)
!!
!!        thilt(:,iw) = 0.0
!!     enddo
!!     !$omp end parallel do
!!
!!  endif

end subroutine prog_wrtv


!=========================================================================


subroutine prog_wrt_begs( iw, istage, vmca, wmsca, dts,         &
                          thil_upv, vxe_upv, vye_upv, vze_upv,  &
                          thil_upw, vxe_upw, vye_upw, vze_upw,  &
                          thilt_short, vmxet_rk, vmyet_rk, vmzet_rk, &
                          rho0, rth0, wmc0 )

  use mem_ijtabs,  only: itab_w
  use mem_basic,   only: wmc, rho, thil, wc, press, vxe, vye
  use misc_coms,   only: deltax, nxp, initial, dn01d, th01d
  use consts_coms, only: cpocv, rocv, omega2, pi1, pio180, r8
  use mem_grid,    only: mza, mva, mwa, lpv, lpw, arw, wnx, wny, wnz, volt, &
                         gravm, volti, volwi, glatw, glonw, &
                         zwgt_top8, zwgt_bot8, gdz_wgtm8, gdz_wgtp8, &
                         zwgt_top, zwgt_bot
  use tridiag,     only: tridiffo
  use oname_coms,  only: nl
! use mem_rayf,    only: dorayfw, rayf_cofw, krayfw_bot, &
!                        dorayf, rayf_cof, krayf_bot
  use mem_nudge,   only: rhot_nud, nudflag

  implicit none

  integer, intent(in) :: iw, istage
  real,    intent(in) :: dts

  real, intent(in)    :: vmca(mza,mva)
  real, intent(inout) :: wmsca(mza,mwa)
  real, intent(in)    :: thil_upv(mza,mva)
  real, intent(in)    :: vxe_upv (mza,mva)
  real, intent(in)    :: vye_upv (mza,mva)
  real, intent(in)    :: vze_upv (mza,mva)

  real, intent(in)    :: thil_upw(mza,mwa)
  real, intent(in)    :: vxe_upw (mza,mwa)
  real, intent(in)    :: vye_upw (mza,mwa)
  real, intent(in)    :: vze_upw (mza,mwa)

  real, intent(in)    :: thilt_short(mza,mwa)
  real, intent(inout) :: vmxet_rk(mza,mwa)
  real, intent(inout) :: vmyet_rk(mza,mwa)
  real, intent(inout) :: vmzet_rk(mza,mwa)

  real,    intent(in) :: rth0(mza,mva)
  real,    intent(in) :: wmc0(mza,mva)
  real,    intent(in) :: rho0(mza,mva)
! real(r8), intent(in) :: rho0(mza,mva)

  integer :: jv, iv, iwn
  integer :: k, ka, kbv
  integer :: npoly

  real :: c6, c8, c9, c10, c11
  real :: dirv, vmarv
  real :: del_rhothil
  real :: rad0_swtc, rad_swtc, topo_swtc

  ! Vertical implicit scheme weighting parameters

  real, parameter :: fr  = .55 ! rho
  real, parameter :: fp2 = .55 ! press in W
  real, parameter :: fp3 = .70 ! press in V

  real, parameter :: pc2 = fp2 * cpocv
  real, parameter :: pc3 = fp3 * cpocv

  ! Automatic arrays

  real     :: wmarw    (mza)
  real(r8) :: hflux_rho(mza)
  real     :: delex_rho(mza)

  real :: rhothil  (mza)
! real(r8) :: rhothil  (mza)
! real(r8) :: press_t  (mza)
  real :: press_t  (mza)

  real(r8) :: press_ex (mza)
  real :: mass_ex  (mza)
! real(r8) :: rd_rt_w  (mza)

  real :: delex_wm     (mza)
  real :: delex_rhothil(mza)
  real :: mass

! real :: hflux_thil(mza)
! real :: hflux_vxe (mza)
! real :: hflux_vye (mza)
! real :: hflux_vze (mza)

  real :: thilt_rk(mza)

  real :: vflux_thil(mza)
  real :: vflux_vxe (mza)
  real :: vflux_vye (mza)
  real :: vflux_vze (mza)

  real :: b1(mza), b2(mza)
  real :: b7

  real :: b4(mza),b5(mza),b6(mza),b10(mza),b15(mza)
  real :: b21, b22, b23, b24, b25, b26
  real :: b31(mza),b32(mza),b33(mza),b34(mza)

  real :: b8, b9, b11, b12, b13, b14

  ka = lpw(iw)

  thilt_rk(ka:mza) = thilt_short(ka:mza,iw)

  ! Loop over W levels

!!  do k = ka, mza-1
!!
!!     ! vertical mass flux for time level t
!!
!!     wmarw(k) = wmc(k,iw) * arw(k,iw)
!!
!!     ! vertical fluxes
!!
!!     vflux_thil(k) = wmarw(k) * thil_upw(k,iw)
!!     vflux_vxe(k)  = wmarw(k) * vxe_upw (k,iw)
!!     vflux_vye(k)  = wmarw(k) * vye_upw (k,iw)
!!     vflux_vze(k)  = wmarw(k) * vze_upw (k,iw)
!!
!!  enddo

  ! Initialize horizontal advection arrays to zero

  hflux_rho  = 0._r8
! hflux_thil = 0.
! hflux_vxe  = 0.
! hflux_vye  = 0.
! hflux_vze  = 0.

  ! Number of edges of this IW polygon

  npoly = itab_w(iw)%npoly

  ! Sum advective tendencies over V neighbors of this cell

  !dir$ loop count max=7
  do jv = 1,npoly

     iv   = itab_w(iw)%iv(jv)
     kbv  = lpv(iv)
     dirv = itab_w(iw)%dirv(jv)

     ! Loop over T levels
     do k = kbv, mza

       vmarv = dirv * vmca(k,iv)

        ! Sum horizontal advection fluxes over V faces

        hflux_rho(k)   = hflux_rho(k)   + vmarv

        thilt_rk(k)    = thilt_rk(k)    + vmarv * thil_upv(k,iv)
        vmxet_rk(k,iw) = vmxet_rk(k,iw) + vmarv * vxe_upv (k,iv)
        vmyet_rk(k,iw) = vmyet_rk(k,iw) + vmarv * vye_upv (k,iv)
        vmzet_rk(k,iw) = vmzet_rk(k,iw) + vmarv * vze_upv (k,iv)
     enddo

  enddo

  ! Add Coriolis terms to momentum tendencies

  do k = ka, mza
     mass = rho(k,iw) * real(volt(k,iw))

     vmxet_rk(k,iw) = vmxet_rk(k,iw) + mass * omega2 * vye(k,iw)
     vmyet_rk(k,iw) = vmyet_rk(k,iw) - mass * omega2 * vxe(k,iw)
  enddo

  ! RAYLEIGH FRICTION ON THIL - only for horizontally homogeneous initialization
!  if (dorayf .and. initial == 1) then

!     do k = krayf_bot, mza
!        thilt_short(k,iw) = thilt_short(k,iw) &
!                          + rayf_cof(k) * dn01d(k) * (th01d(k) - thil(k,iw))
!     enddo

     ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions
     ! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
     ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
     ! rayfx = 0.   ! Default: no extra RAYF
     ! do k = ka, mza
     !    thilt(k,iw) = thilt(k,iw) + max(rayf_cof(k),rayfx)  &
     !                * dn01d(k) * (th01d(k) - thil(k,iw))
     ! enddo

!  endif ! (dorayf .and. initial == 1)

  c8  = dts * pc2
  c9  =-dts * fr

  do k = ka, mza

    b10(k) = dts * volti(k,iw)

     rhothil(k) = rho(k,iw) * thil(k,iw)

     ! Explicit density tendency

     delex_rho(k) = b10(k) * real(hflux_rho(k))

     if (nudflag == 1) delex_rho(k) = delex_rho(k) + dts * rhot_nud(k,iw)

     ! Explicit density-thil tendency

     delex_rhothil(k) = b10(k) * thilt_rk(k)

     if (istage > 1) then
        delex_rho    (k) = delex_rho    (k) + rho0(k,iw) - rho(k,iw)
        delex_rhothil(k) = delex_rhothil(k) + rth0(k,iw) - rhothil(k)
     endif

     b5 (k) = alpha_press(k,iw) * rhothil(k) ** rocv
     b15(k) = b10(k) * b5(k)

     press_ex(k) = b5(k) * (rhothil(k) + pc2 * delex_rhothil(k))

     mass_ex(k) = real(volt(k,iw)) * (rho(k,iw) + fr * delex_rho(k))

  enddo

  ! Vertical loop over W levels

  do k = ka, mza-1

     b4(k) = gravm(k) * volwi(k,iw)

     ! Explicit vertical momentum tendency

     delex_wm(k) = wmc0(k,iw) + dts * ( & !wmt_other(k) &

          ! Pressure gradient at W level

          + real( press_ex(k) - press_ex(k+1) ) * pwfac(k,iw) &

          ! Volume-weighted buoyancy forcing in vertical between T levels

          - (mass_ex(k) + mass_ex(k+1)) * b4(k)  &

          ! Volume-weighted momentum forcings in vertical between T levels

          + ( wnx(iw) * (vmxet_rk(k,iw) + vmxet_rk(k+1,iw)) &
            + wny(iw) * (vmyet_rk(k,iw) + vmyet_rk(k+1,iw)) &
            + wnz(iw) * (vmzet_rk(k,iw) + vmzet_rk(k+1,iw)) ) * volwi(k,iw) )

  enddo

  ! Fill matrix coefficients for implicit update of WM

  ! Loop over W pts
!  do k = ka, mza-1
!     b1(k)  = wnx(iw) * vxe_upw(k,iw) &
!            + wny(iw) * vye_upw(k,iw) &
!            + wnz(iw) * vze_upw(k,iw)
!  enddo

! b2(ka-1) = thil(ka,iw)
! b2(mza)  = thil(mza,iw)
! b2(ka-1)= 2. * thil(ka,iw)
! b2(mza) = 2. * thil(mza,iw)

! b1(ka-1) = 0.
! b1(mza)  = 0.

  ! Loop over T pts
! do k = ka, mza
!    b1(k)  = wc(k,iw) + wc(k-1,iw)
!    b5(k)  = press_t(k) / rhothil(k)
!    b6(k)  = c6  * volti(k,iw)
!     b10(k) = c10 * volti(k,iw)
!  enddo

  ! Loop over W pts
  do k = ka, mza-1

!    b7  = c11 * volwi(k,iw)
!    b7(k)  = c6 * volwi(k,iw)

!    b9(k)  = c9 * dzim(k)
!    b9(k)  = c9 * grav * volwi(k,iw) * vfac(k)

!    b11 = b8 * b5(k)
!    b12 = b8 * b5(k+1)

!     b13(k) = b9(k)  * dzt_top(k)
!     b14(k) = b9(k)  * dzt_bot(k+1)

!    b9 = c9 * b4(k)
     b9 = c9 * dts * b4(k)

!    b13 = b9 * volt(k,iw)
!    b14 = b9 * volt(k+1,iw)

!     b13 = c9 * gdz_wgtm8(k)
!     b14 = c9 * gdz_wgtp8(k)

!    b13(k) = b9(k)  * volt(k  ,iw)
!    b14(k) = b9(k)  * volt(k+1,iw)

!    b21 = b7 * b1(k-1)
!    b22 = b7 * b1(k+1)

     b7  = dts * volwi(k,iw)

     b21 = b7 * wc(k-1,iw)
     b22 = b7 * wc(k+1,iw)
!    b21(k) = b7(k)  * wc(k-1,iw)
!    b22(k) = b7(k)  * wc(k+1,iw)

!    b23(k) = b11(k) * b6(k)
!    b24(k) = b12(k) * b6(k+1)

     b8  = c8 * pwfac(k,iw)

     b23 = b8 * b15(k)
     b24 = b8 * b15(k+1)

!    b25 = b13 * b10(k)
!    b26 = b14 * b10(k+1)

!    b25(k) = b9(k) * c10
!    b26(k) = b9(k) * c10

!            * (b2(k) * (b23(k) + b24(k)) + b25(k) - b26(k))
!           * (b22(k) - b21(k) + b2(k) * (b23(k) + b24(k)) + b25(k) - b26(k))

     b31(k) =        - arw(k-1,iw) * (b21 + b23 * thil_upw(k-1,iw) + b9)

     b32(k) = b20(k) + arw(k  ,iw) * thil_upw(k,iw) * (b23 + b24)

     b33(k) =          arw(k+1,iw) * (b22 - b24 * thil_upw(k+1,iw) + b9)

  enddo

  ! Solve implicit tri-diagonal matrix equation for delta WM (del_wm)

  if (ka <= mza-1) then
     call tridiffo(mza,ka,mza-1,b31,b32,b33,delex_wm,wmc(:,iw))
  endif

  ! Vertical loop over W points

  do k = ka, mza-1

     wmarw(k)  = wmc(k,iw) * arw(k,iw)

     vflux_thil(k) = wmarw(k) * thil_upw(k,iw)
     vflux_vxe(k)  = wmarw(k) * vxe_upw (k,iw)
     vflux_vye(k)  = wmarw(k) * vye_upw (k,iw)
     vflux_vze(k)  = wmarw(k) * vze_upw (k,iw)

     ! Add vertical momentum at (t+fw) to array for long-timestep scalar transport

     if (istage == nrk) then
        wmsca(k,iw) = wmsca(k,iw) + wmarw(k)
     endif

  enddo

  ! Set bottom & top vertical advective mass and heat fluxes to zero

  wmarw(ka-1) = 0.
  wmarw(mza)  = 0.

  vflux_thil(ka-1) = 0.
  vflux_thil(mza)  = 0.

  vflux_vxe(ka-1) = 0.
  vflux_vxe(mza)  = 0.

  vflux_vye(ka-1) = 0.
  vflux_vye(mza)  = 0.

  vflux_vze(ka-1) = 0.
  vflux_vze(mza)  = 0.

  ! Vertical loop over T points

  do k = ka, mza

     ! Change of rho from (t) to (t+1)

     rho(k,iw) = rho(k,iw) + delex_rho(k) + b10(k) * (wmarw(k-1) - wmarw(k))

!     rho(k,iw) = rho0(k,iw) &
!               + b10(k) * real(hflux_rho(k) + wmarw(k-1) - wmarw(k))

     ! Change of rhothil from (t) to (t+1)

     del_rhothil = delex_rhothil(k) &
                 + b10(k) * (vflux_thil(k-1) - vflux_thil(k))

     ! Update pressure from (t) to (t+tp)

     press(k,iw) = b5(k) * (rhothil(k) + pc3 * del_rhothil)

     ! Update thil from (t) to (t+1)

!    thil(k,iw) = (rhothil(k) + del_rhothil) / real(rho(k,iw))
     thil(k,iw) = (rhothil(k) + del_rhothil) / rho(k,iw)

     ! Update volume-weighted momentum tendencies at T due to implicit flux correction

     vmxet_rk(k,iw) = vmxet_rk(k,iw) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet_rk(k,iw) = vmyet_rk(k,iw) + vflux_vye(k-1) - vflux_vye(k)
     vmzet_rk(k,iw) = vmzet_rk(k,iw) + vflux_vze(k-1) - vflux_vze(k)

  enddo

!!  if (lve2(iw) > 0) then
!!
!!     do ksw = 1, lve2(iw)
!!        k = ksw + ka - 1
!!
!!        ! Zero out vxe2, vye2, vze2 prior to new diagnosis
!!        vxe2(ksw,iw) = 0.
!!        vye2(ksw,iw) = 0.
!!        vze2(ksw,iw) = 0.
!!
!!        ! Estimate velocity in T cells at (t+1) by prognostic method
!!        vxe1(k) = vxe(k,iw) + dts * vmxet_rk(k,iw) * volti(k,iw) / real(rho(k,iw))
!!        vye1(k) = vye(k,iw) + dts * vmyet_rk(k,iw) * volti(k,iw) / real(rho(k,iw))
!!        vze1(k) = vze(k,iw) + dts * vmzet_rk(k,iw) * volti(k,iw) / real(rho(k,iw))
!!     enddo
!!
!!     ! Loop over adjacent V faces
!!
!!     do jv = 1, npoly
!!        iv = itab_w(iw)%iv(jv)
!!
!!        ! Project full-forward-time vxe1, vye1, vze1 onto V faces that are
!!        ! below ground, and then project back to vxe2, vye2, vze2
!!
!!        if (lpv(iv) > ka) then
!!           do k = ka, lpv(iv) - 1
!!              ksw = k - ka + 1
!!              vmt1 = vnx(iv) * vxe1(k) + vny(iv) * vye1(k) + vnz(iv) * vze1(k)
!!
!!              vxe2(ksw,iw) = vxe2(ksw,iw) + itab_w(iw)%ecvec_vx(jv) * vmt1
!!              vye2(ksw,iw) = vye2(ksw,iw) + itab_w(iw)%ecvec_vy(jv) * vmt1
!!              vze2(ksw,iw) = vze2(ksw,iw) + itab_w(iw)%ecvec_vz(jv) * vmt1
!!           enddo
!!        endif
!!
!!     enddo
!!  endif

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
        press(k,iw) = rho(k,iw) + topo_swtc + fp2 * delex_rho(k)
     enddo
  endif

  ! Vertical loop over W points

  do k = ka, mza-1
     wc (k,iw) = wmc(k,iw) / (zwgt_bot(k) * rho(k,iw) + zwgt_top(k) * rho(k+1,iw))
  enddo

  ! Set top & bottom values of WC

!  wc(ka-1,iw) = wc(ka,iw)

!  wc(1:ka-2,iw) = 0.
!  wc(mza   ,iw) = 0.

  ! Bottom boundary for THIL

  do k = 1, ka-1
     thil(k,iw) = thil(ka,iw)
  enddo

end subroutine prog_wrt_begs

!==========================================================================

subroutine comp_alpha_press(mrl)

  use mem_grid,    only: lpw, lpv, mza, dzim, dniv, zfacit
  use mem_ijtabs,  only: jtab_w, jtw_prog, jtab_v, jtv_wadj, itab_v
  use consts_coms, only: pc1, rdry, rvap, cpocv, r8
  use mem_basic,   only: rr_v, rr_w, theta, thil

  implicit none

  integer, intent(in)  :: mrl
  integer              :: j, iw, k, iv, iw1, iw2

  !$omp parallel
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Evaluate alpha coefficient for pressure
     do k = lpw(iw), mza
!       alpha_press(k,iw) = pc1 * (rdry + rvap * rr_v(k,iw)) ** cpocv
        alpha_press(k,iw) = pc1 * ( (rdry + rvap * rr_v(k,iw)) &
                                  * theta(k,iw) / thil(k,iw) ) ** cpocv
     enddo

     do k = lpw(iw), mza-1
!       pwfac(k,iw) = dzim(k) * ( zwgt_bot(k) / (1. + rr_w(k  ,iw)) &
!                               + zwgt_top(k) / (1. + rr_w(k+1,iw)) )
        pwfac(k,iw) = dzim(k) * (1.0 - 0.5 *(rr_w(k,iw) + rr_w(k+1,iw)))
     enddo
  enddo
  !$omp end do nowait

  !$omp do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
!       pvfac(k,iv) = dnivo2(iv) * zfacit(k) * ( 1. / (1. + rr_w(k,iw1)) + 1. / (1. + rr_w(k,iw2)) )
        pvfac(k,iv) = dniv(iv) * zfacit(k) * (1.0 - 0.5 * (rr_w(k,iw1) + rr_w(k,iw2)))
     enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

end subroutine comp_alpha_press

!============================================================================

subroutine prog_v_begs( iv, dts, vmc0, vmxet_rk, vmyet_rk, vmzet_rk, vmt_short )

  use mem_ijtabs,  only: itab_v
  use mem_basic,   only: vc, press, vmc, rho
  use consts_coms, only: eradi, gravo2
  use mem_grid,    only: mza, mva, mwa, lpv, dniu, dniv, &
                         zfacit, vnx, vny, vnz, dnivo2, volvi, lpw
  use oname_coms,  only: nl
  use mem_rayf

  implicit none

  integer, intent(in) :: iv
  real,    intent(in) :: dts
  real,    intent(in) :: vmxet_rk(mza,mwa)
  real,    intent(in) :: vmyet_rk(mza,mwa)
  real,    intent(in) :: vmzet_rk(mza,mwa)

  real,    intent(in) :: vmt_short(mza,mva)
  real,    intent(in) :: vmc0     (mza,mva)

  integer :: k, kb, kbot
  integer :: iw1, iw2
  real    :: pgf(mza)

  ! Extract neighbor indices and coefficients for this point in the V stencil

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  kb = lpv(iv)

  ! Vertical loop over V points to compute pressure gradient force

  do k = kb, mza
     pgf(k) = pvfac(k,iv) * real( press(k,iw1) - press(k,iw2) )
  enddo

  ! For shallow water test cases 2 & 5, rho & press are
  ! interpreted as water depth & height

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     do k = kb, mza
        pgf(k) = pgf(k) * gravo2 * (rho(k,iw1) + rho(k,iw2))
     enddo
  endif

!!  if (dorayfdiv) then
!!
!!     ! Vertical loop over V points
!!     do k = krayfdiv_bot, mza
!!
!!        ! Divergence in IW1 excluding IV
!!        sum1 = div2d(k,iw1) - vmc(k,iv) * dnu(iv)
!!
!!        ! Divergence in IW2 excluding IV
!!        sum2 = div2d(k,iw2) + vmc(k,iv) * dnu(iv)
!!
!!        ! VMP value that would equalize horizontal divergence in IW1 and IW2
!!        ! cells (assuming no blockage by topography)
!!
!!        vmp_eqdiv = ( arw0(iw1) * sum2 - arw0(iw2) * sum1 ) * rayfdivfac(iv)
!!
!!        pgf(k) = pgf(k) + rayf_cofdiv(k) * (vmp_eqdiv - vmc(k,iv))
!!
!!     enddo
!!
!!  endif ! (dorayfdiv)




  ! Update VMC

!!  if (.not. rotational) then

!     kbot = max(lpw(iw1),lpw(iw2))

!!     if (lpw(iw1) < kbot) then
!!
!!        do k = lpw(iw1), kbot-1
!!           vmc(k,iv) = vmp(k,iv) + dts * ( vmt_short(k,iv) &
!!                     + ( vnx(iv) * vmxet_rk(k,iw1) &
!!                       + vny(iv) * vmyet_rk(k,iw1) &
!!                       + vnz(iv) * vmzet_rk(k,iw1) ) * volvi(k,iv) )
!!
!!           vc(k,iv) = vmc(k,iv) / real( rho(k,iw1) )
!!
!!
!!        enddo
!!
!!     elseif (lpw(iw2) < kbot) then
!!
!!        do k = lpw(iw2), kbot-1
!!           vmc(k,iv) = vmp(k,iv) + dts * ( vmt_short(k,iv) &
!!                     + ( vnx(iv) * vmxet_rk(k,iw2) &
!!                       + vny(iv) * vmyet_rk(k,iw2) &
!!                       + vnz(iv) * vmzet_rk(k,iw2) ) * volvi(k,iv) )
!!
!!           vc(k,iv) = vmc(k,iv) / real( rho(k,iw2) )
!!      enddo
!!
!!     endif

     do k = min(lpw(iw1),lpw(iw2)), kb-1
        vmc(k,iv) = vmc0(k,iv) + dts * vmt_short(k,iv)
     enddo

     do k = kb, mza
!     do k = kbot, mza

        vmc(k,iv) = vmc0(k,iv) + dts * ( vmt_short(k,iv) + pgf(k)    &
                  + ( vnx(iv) * (vmxet_rk(k,iw1) + vmxet_rk(k,iw2)) &
                    + vny(iv) * (vmyet_rk(k,iw1) + vmyet_rk(k,iw2)) &
                    + vnz(iv) * (vmzet_rk(k,iw1) + vmzet_rk(k,iw2)) ) * volvi(k,iv) )

!!!        vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )

     enddo

     kbot = max(lpw(iw1),lpw(iw2))

     do k = kbot, mza
        vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )
     enddo

     if (lpw(iw1) < kbot) then
        do k = lpw(iw1), kbot-1
           vc(k,iv) = vmc(k,iv) / real( rho(k,iw1) )
        enddo
     elseif (lpw(iw2) < kbot) then
        do k = lpw(iw2), kbot-1
           vc(k,iv) = vmc(k,iv) / real( rho(k,iw2) )
        enddo
     endif

!!     if (any(itab_v(iv)%ivglobe == [340103,340105])) then
!!        write(*,*) itab_v(iv)%ivglobe
!!        do k = mza, min(lpw(iw1),lpw(iw2)), -1
!!           write(*,'(I4,10g13.5)') k, vmc(k,iv), vc(k,iv), pgf(k), vmt_short(k,iv), ( vnx(iv) * (vmxet_rk(k,iw1) + vmxet_rk(k,iw2)) &
!!                    + vny(iv) * (vmyet_rk(k,iw1) + vmyet_rk(k,iw2)) &
!!                    + vnz(iv) * (vmzet_rk(k,iw1) + vmzet_rk(k,iw2)) ) * volvi(k,iv)
!!        enddo
!!        write(*,*)
!!     endif
!!


!!  else  ! If using rotational form:
!!
!!     dtso2dnu = 0.5 * dts * dniu(iv)
!!     dtso2dnv = 0.5 * dts * dniv(iv)
!!
!!     ! Vertical loop over V levels
!!
!!     do k = kb,mza
!!        vx = .5 * (vxe(k,iw1) + vxe(k,iw2))
!!        vy = .5 * (vye(k,iw1) + vye(k,iw2))
!!        vz = .5 * (vze(k,iw1) + vze(k,iw2))
!!
!!        uc =    unx(iv) * vx + uny(iv) * vy + unz(iv) * vz
!!        watv = (xev(iv) * vx + yev(iv) * vy + zev(iv) * vz) * eradi
!!
!!        tke1 = .5 * (vxe(k,iw1)**2 + vye(k,iw1)**2 + vze(k,iw1)**2)
!!        tke2 = .5 * (vxe(k,iw2)**2 + vye(k,iw2)**2 + vze(k,iw2)**2)
!!
!!        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@A1
!!        vortp_v = .5 * (vortp(k,im1) + vortp(k,im2)) &
!!
!!             ! APVM method (Ringler et al. 2010) to obtain upwinded vortp at V
!!
!!             + dtso2dnu * uc       * (vortp(k,im1) - vortp(k,im2)) &
!!             + dtso2dnv * vc(k,iv) * (vortp_t(k,iw1) - vortp_t(k,iw2))
!!        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@A2
!!
!!        vmc(k,iv) = vmc(k,iv) + dts * (vmt_short(k,iv) + pgf(k) &
!!
!!                  + ( vnx(iv) * (vmxet_rk(k,iw1) + vmxet_rk(k,iw2))    &
!!                    + vny(iv) * (vmyet_rk(k,iw1) + vmyet_rk(k,iw2))    &
!!                    + vnz(iv) * (vmzet_rk(k,iw1) + vmzet_rk(k,iw2)) ) * volvi(k,iv) &
!!
!!                  + .5 * (rho(k,iw1) + rho(k,iw2))                           &
!!                  * (dniv(iv) * (tke1 - tke2)                                &
!!                  + uc * vortp_v                                             &
!!                  - watv * .5 * (vortn(k-1,iv) + vortn(k,iv))))
!!
!!        vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )
!!     enddo
!!
!!  endif ! rotational

end subroutine prog_v_begs

!===============================================================================

subroutine extra_mix_top_vc( iv, vmt )

  use mem_ijtabs,  only: itab_v
  use mem_basic,   only: vc, rho
  use mem_grid,    only: mza, arw, volvi
  use mem_rayf

  implicit none

  integer, intent(in   ) :: iv
  real,    intent(inout) :: vmt(mza)

  real    :: vflux(mza)
  integer :: k, iw1, iw2

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  vflux (mza)           = 0.0
  vflux(krayfdif_bot-1) = 0.0

  do k = krayfdif_bot, mza-1
    vflux(k) = 0.25 * (arw(k,iw1) + arw(k,iw2)) * rayf_cofdif(k) * (vc(k,iv) - vc(k+1,iv)) &
             * (rho(k+1,iw1) + rho(k+1,iw2) + rho(k,iw1) + rho(k,iw2))
!     vflux(k) = (arw(k,iw1) + arw(k,iw2)) * rayf_cofdif(k) * (vmc(k,iv) - vmc(k+1,iv))
  enddo

  do k = krayfdif_bot, mza
     vmt(k) = vmt(k) + (vflux(k-1) - vflux(k)) * volvi(k,iv)
  enddo

end subroutine extra_mix_top_vc

!===============================================================================

!!subroutine prep_rotational(mrl)
!!
!!  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtab_m, itab_m, &
!!                          jtm_vadj, jtv_prog, jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp
!!  use mem_basic,    only: vc, wc
!!  use mem_grid,     only: mza, lpm, lpv, lpw, dzim, zfact, &
!!                          zfacit, zfacim, dnv, dniv, arm0
!!  use obnd,         only: lbcopy_m, lbcopy_w
!!  use misc_coms,    only: iparallel
!!  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
!!
!!  implicit none
!!
!!  integer, intent(in) :: mrl
!!
!!  integer :: j, iv, k, kb
!!  integer :: iw, iw1, iw2
!!  integer :: im,npoly,jv
!!  integer :: jm
!!  real    :: arm0i
!!
!!! IF USING ROTATIONAL METHOD, HORIZONTAL LOOP OVER M/P COLUMNS
!!! FOR COMPUTING VERTICAL VORTICITY.
!!
!!  !$omp parallel do private(im,npoly,kb,jv,iv,k,arm0i)
!!  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)
!!
!!     npoly = itab_m(im)%npoly
!!     kb    = lpm(im)
!!
!!     vortp(:,im) = 0.
!!
!!     ! Loop over V neighbors to evaluate circulation around M (at time T)
!!     do jv = 1,npoly
!!        iv = itab_m(im)%iv(jv)
!!
!!        if (itab_v(iv)%im(2) == im) then
!!
!!           do k = kb,mza
!!              vortp(k,im) = vortp(k,im) + vc(k,iv) * dnv(iv)
!!           enddo
!!
!!        else
!!
!!           do k = kb,mza
!!              vortp(k,im) = vortp(k,im) - vc(k,iv) * dnv(iv)
!!           enddo
!!
!!        endif
!!     enddo
!!
!!     ! Convert circulation to relative vertical vorticity at M
!!     ! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
!!     ! divide their quotient by zfact)
!!
!!     arm0i = 1. / arm0(im)
!!
!!     do k = kb,mza
!!        vortp(k,im) = vortp(k,im) * arm0i * zfacit(k)
!!     enddo
!!
!!  enddo
!!  !$omp end parallel do 
!!
!!  ! PARALLEL SEND/RECV OF VORTP
!!
!!  if (iparallel == 1) then
!!     call mpi_send_m(mrl, rvara1=vortp)
!!     call mpi_recv_m(mrl, rvara1=vortp)
!!  endif
!!  call lbcopy_m(mrl, a1=vortp)
!!
!!! IF USING ROTATIONAL FORM FOR PROGNOSING HORIZONTAL MOMENTUM,
!!! EVALUATE VORTP_T AND VORTN AND CORIOLIS FORCING TERMS
!!
!!  ! Horizontal loop over W/T columns
!!  !$omp parallel do private(iw,npoly,kb,k,jm,im)
!!  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!!
!!     npoly = itab_w(iw)%npoly
!!     kb = lpw(iw)
!!
!!     vortp_t(:,iw) = 0.
!!
!!     ! Loop over M neighbors of W
!!     do jm = 1,npoly
!!        im = itab_w(iw)%im(jm)
!!
!!        ! Vertical loop over T levels; average vorticity from P points to T point.
!!        do k = kb,mza
!!           vortp_t(k,iw) = vortp_t(k,iw) + itab_w(iw)%farm(jm) * vortp(k,im)
!!        enddo
!!     enddo
!!
!!  enddo
!!  !$omp end parallel do
!!
!!  if (iparallel == 1) call mpi_send_w(mrl, rvara1=vortp_t)
!!  if (iparallel == 1) call mpi_recv_w(mrl, rvara1=vortp_t)
!!  call lbcopy_w(mrl, a1=vortp_t)
!!
!!  ! Horizontal loop over V/N columns
!!  !$omp parallel do private(iv,iw1,iw2,k,kb) 
!!  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!!
!!     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!!     kb = lpv(iv)
!!
!!     ! Vertical loop over N levels; compute horizontal relative vorticity
!!     ! and vertical velocity at N points (at time T)
!!     do k = kb,mza-1
!!        vortn(k,iv) = zfacim(k) * ((wc(k,iw1) - wc(k,iw2)) * dniv(iv) &
!!                    + (vc(k+1,iv) * zfact(k+1) - vc(k,iv) * zfact(k)) * dzim(k))
!!     enddo
!!
!!     vortn(1:kb-1,iv) = 0.
!!     vortn(mza,iv) = 0.
!!
!!  enddo
!!  !$omp end parallel do
!!
!!end subroutine prep_rotational


subroutine wind_vec_at_v(iv, vxe_v, vye_v, vze_v)

  use mem_grid,   only: mza, lpv, vnx, vny, vnz
  use mem_ijtabs, only: itab_v
  use mem_basic,  only: vc, vxe, vye, vze

  implicit none

  integer, intent( in) :: iv
  real,    intent(out) :: vxe_v(mza), vye_v(mza), vze_v(mza)
  integer              :: iw1, iw2, k
  real                 :: vx, vy, vz
  real                 :: vcc, vdiff

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  do k = lpv(iv), mza
     vx = 0.5 * (vxe(k,iw1) + vxe(k,iw2))
     vy = 0.5 * (vye(k,iw1) + vye(k,iw2))
     vz = 0.5 * (vze(k,iw1) + vze(k,iw2))

     vcc   = vnx(iv) * vx + vny(iv) * vy + vnz(iv) * vz
     vdiff = vc(k,iv) - vcc

     vxe_v(k) = vx + vdiff * vnx(iv)
     vye_v(k) = vy + vdiff * vny(iv)
     vze_v(k) = vz + vdiff * vnz(iv)
  enddo

end subroutine wind_vec_at_v


subroutine grad_t2d_v(iw, vxe_v, vye_v, vze_v, &
                          gxps_vxe, gyps_vxe, gxps_vye, gyps_vye, &
                          gxps_vze, gyps_vze)

  use mem_basic,  only: vxe, vye, vze
  use mem_grid,   only: mza, mva, lpv, gxps_coef, gyps_coef
  use mem_ijtabs, only: itab_w

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: vxe_v(mza,mva), vye_v(mza,mva), vze_v(mza,mva)

  real,    intent(out) :: gxps_vxe(mza), gyps_vxe(mza)
  real,    intent(out) :: gxps_vye(mza), gyps_vye(mza)
  real,    intent(out) :: gxps_vze(mza), gyps_vze(mza)

  integer :: jv, iv, k
  real    :: dvx, dvy, dvz, gx, gy

  gxps_vxe = 0.
  gyps_vxe = 0.

  gxps_vye = 0.
  gyps_vye = 0.

  gxps_vze = 0.
  gyps_vze = 0.

! Loop over V faces of this T cell

  do jv = 1, itab_w(iw)%npoly
     iv = itab_w(iw)%iv(jv)

! Vertical loop over T levels
! Zero-gradient lateral B.C. below lpv(iv)

     gx = 2.0 * gxps_coef(iw,jv)
     gy = 2.0 * gyps_coef(iw,jv)

     do k = lpv(iv), mza
        dvx = vxe_v(k,iv) - vxe(k,iw)
        dvy = vye_v(k,iv) - vye(k,iw)
        dvz = vze_v(k,iv) - vze(k,iw)

        gxps_vxe(k) = gxps_vxe(k) + gx * dvx
        gyps_vxe(k) = gyps_vxe(k) + gy * dvx

        gxps_vye(k) = gxps_vye(k) + gx * dvy
        gyps_vye(k) = gyps_vye(k) + gy * dvy

        gxps_vze(k) = gxps_vze(k) + gx * dvz
        gyps_vze(k) = gyps_vze(k) + gy * dvz

!!        gxps_vxe(k) = gxps_vxe(k) + gxps_coefv(iw,jv) * dvx
!!        gyps_vxe(k) = gyps_vxe(k) + gyps_coefv(iw,jv) * dvx
!!
!!        gxps_vye(k) = gxps_vye(k) + gxps_coefv(iw,jv) * dvy
!!        gyps_vye(k) = gyps_vye(k) + gyps_coefv(iw,jv) * dvy
!!
!!        gxps_vze(k) = gxps_vze(k) + gxps_coefv(iw,jv) * dvz
!!        gyps_vze(k) = gyps_vze(k) + gyps_coefv(iw,jv) * dvz
     enddo

  enddo

end subroutine grad_t2d_v


end module wrtv_mem
