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

module wrtv_orig
  implicit none

  real, allocatable :: vortp   (:,:)
  real, allocatable :: vortp_t (:,:)
  real, allocatable :: vortn   (:,:)

  logical, save     :: rotational = .false.

  private
  public prog_wrtv_orig, init_wrtv_orig

contains

!===============================================================================

subroutine init_wrtv_orig()

  use mem_grid,   only: mma, mwa, mva, mza
  use misc_coms,  only: rinit
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

end subroutine init_wrtv_orig

!===============================================================================

subroutine prog_wrtv_orig()

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, istp, mrl_begl, &
                          mrl_begs, mrl_ends, jtm_vadj, jtv_prog, jtw_wadj, &
                          jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp
  use mem_basic,    only: rho, thil, press, wmc, wc, vmc, vc, vxe, vye, vze, vmp, &
                          vmsc, wmsc, vxesc, vyesc, vzesc
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, arw, arv, volt
  use mem_tend,     only: thilt, vmxet, vmyet, vmzet, vmt
  use misc_coms,    only: iparallel, dtsm
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use oplot_coms,   only: op
  use obnd,         only: lbcopy_v, lbcopy_w
  use vel_t3d,      only: vel_t3d_hex
  use oname_coms,   only: nl
  use mem_turb,     only: akmodx, akhodx, khtopv, kmtopv, khtop
  use mem_rayf,     only: dorayfmix, rayf_mix_top_vxe
  use grad_lib,     only: grad_t2d, grad_z
  use pbl_drivers,  only: solve_eddy_diff_heat, solve_eddy_diff_vxe
  use vel_t3d,      only: diagvel_t3d

!!  use mem_adv,      only: dxps_w, dyps_w, dzps_w, &
!!                          dxyps_w, dxxps_w, dyyps_w, dzzps_w, &
!!                          dxps_v, dyps_v, dzps_v, &
!!                          dxyps_v, dxxps_v, dyyps_v, dzzps_v, &
!!                          gxps_scp, gyps_scp, gzps_scp, &
!!                          gxyps_scp, gxxps_scp, gyyps_scp, gzzps_scp, &
!!                          gxps_vxe, gyps_vxe, gzps_vxe, &
!!                          gxyps_vxe, gxxps_vxe, gyyps_vxe, gzzps_vxe, &
!!                          gxps_vye, gyps_vye, gzps_vye, &
!!                          gxyps_vye, gxxps_vye, gyyps_vye, gzzps_vye, &
!!                          gxps_vze, gyps_vze, gzps_vze, &
!!                          gxyps_vze, gxxps_vze, gyyps_vze, gzzps_vze

  implicit none

  integer :: j, iv, iw, k, mrl, kd, iwd, jv, iwn, iw1, iw2

! automatic arrays

  integer :: iwdepv(mza,mva) ! donor cell IW index for V face
  integer :: kdepw(mza,mwa)  ! donor cell K index for W face

  real :: v4, dt
  real :: vmcf(mza,mva) ! Time-extrapolated VMC

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

  real :: vmt_short(mza,mva)

  real :: gxps_the(mza,mwa), gyps_the(mza,mwa), gzps_the(mza,mwa)
  real :: gxps_vxe(mza,mwa), gyps_vxe(mza,mwa), gzps_vxe(mza,mwa)
  real :: gxps_vye(mza,mwa), gyps_vye(mza,mwa), gzps_vye(mza,mwa)
  real :: gxps_vze(mza,mwa), gyps_vze(mza,mwa), gzps_vze(mza,mwa)

  real :: vcf(mza)
  real :: dxps_w(mza), dyps_w(mza), dzps_w(mza)
  real :: dxps_v(mza), dyps_v(mza), dzps_v(mza)

! THIS ROUTINE IS PERFORMED AT THE END OF EACH SMALL TIMESTEP

  mrl = mrl_ends(istp)
  if (mrl == 0) return

  dt = dtsm(1)

! INCLUDE THE LONG TIMESTEP THIL AND MOMENTUM TENDENCIES IN EACH SHORT TIMESTEP

  !$omp parallel do private(iw,k,v4,jv,iv,iwn)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        v4 = real(volt(k,iw))

        thilt_short(k,iw) = thilt(k,iw) * v4
        vmxet_short(k,iw) = vmxet(k,iw) * v4
        vmyet_short(k,iw) = vmyet(k,iw) * v4
        vmzet_short(k,iw) = vmzet(k,iw) * v4
     enddo

     call solve_eddy_diff_heat(iw, thilt_short(:,iw))

     call solve_eddy_diff_vxe (iw, vmxet_short(:,iw), vmyet_short(:,iw), &
                                   vmzet_short(:,iw))

     call grad_t2d(iw, thil, gxps_the(:,iw), gyps_the(:,iw))
     call grad_t2d(iw, vxe,  gxps_vxe(:,iw), gyps_vxe(:,iw))
     call grad_t2d(iw, vye,  gxps_vye(:,iw), gyps_vye(:,iw))
     call grad_t2d(iw, vze,  gxps_vze(:,iw), gyps_vze(:,iw))

     do jv = 1, itab_w(iw)%npoly
        iv  = itab_w(iw)%iv(jv)
        iwn = itab_w(iw)%iw(jv)

        ! Loop over T levels
        do k = lpv(iv), khtopv(iv)
           thilt_short(k,iw) = thilt_short(k,iw) &
                             + akhodx(k,iv) * (thil(k,iwn) - thil(k,iw))
        enddo

        ! Loop over T levels
        do k = lpv(iv), kmtopv(iv)
           vmxet_short(k,iw) = vmxet_short(k,iw) &
                             + akmodx(k,iv) * (vxe(k,iwn) - vxe(k,iw))
           vmyet_short(k,iw) = vmyet_short(k,iw) &
                             + akmodx(k,iv) * (vye(k,iwn) - vye(k,iw))
           vmzet_short(k,iw) = vmzet_short(k,iw) &
                             + akmodx(k,iv) * (vze(k,iwn) - vze(k,iw))
        enddo
     enddo

     ! Rayleigh friction on velocity gradient

     if (dorayfmix) then
        call rayf_mix_top_vxe (iw, vmxet_short(:,iw), vmyet_short(:,iw), vmzet_short(:,iw))
     endif

  enddo
  !$omp end parallel do

  ! MPI send of THIL, VXE, VYE, and VZE gradient components

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=gxps_the, rvara2=gyps_the, &
                          rvara3=gxps_vxe, rvara4=gyps_vxe, &
                          rvara5=gxps_vye, rvara6=gyps_vye, &
                          rvara7=gxps_vze, rvara8=gyps_vze  )
  endif

  !$omp parallel private(dxps_w,dyps_w,dzps_w)
  !$omp do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     ! Extrapolate VM to time T + 1/2; update VMP

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vmcf (k,iv) = 1.5 * vmc(k,iv) - 0.5 * vmp(k,iv)
        vmsc (k,iv) = vmsc(k,iv) + vmcf(k,iv)
        vmp  (k,iv) = vmc (k,iv)
     enddo

  enddo
  !$omp end do nowait

  ! EVALUATE VERTICAL GRADIENT OF THIL, VXE, VYE, AND VZE FOR BEGS

  !$omp do private(iw)
  do j = 1, jtab_w(jtw_wadj)%jend(mrl); iw = jtab_w(jtw_wadj)%iw(j)
     call grad_z(iw, thil(:,iw), gzps_the(:,iw))
     call grad_z(iw, vxe (:,iw), gzps_vxe(:,iw))
     call grad_z(iw, vye (:,iw), gzps_vye(:,iw))
     call grad_z(iw, vze (:,iw), gzps_vze(:,iw))
  enddo
  !$omp end do nowait

  !$omp do private(iw,k,kd)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     call donorpointw(iw, dt, wc(:,iw), vxe(:,iw), vye(:,iw), vze(:,iw), &
                      dxps_w, dyps_w, dzps_w)

     do k = lpw(iw), mza-1

        if (wmc(k,iw) > 0.) then
           kd = k
        else
           kd = k+1
        endif

        thil_upw(k,iw) =                  thil(kd,iw) &
                       +  dxps_w(k) * gxps_the(kd,iw) &
                       +  dyps_w(k) * gyps_the(kd,iw) &
                       +  dzps_w(k) * gzps_the(kd,iw)

        vxe_upw(k,iw)  =                   vxe(kd,iw) &
                       +  dxps_w(k) * gxps_vxe(kd,iw) &
                       +  dyps_w(k) * gyps_vxe(kd,iw) &
                       +  dzps_w(k) * gzps_vxe(kd,iw)

        vye_upw(k,iw)  =                   vye(kd,iw) &
                       +  dxps_w(k) * gxps_vye(kd,iw) &
                       +  dyps_w(k) * gyps_vye(kd,iw) &
                       +  dzps_w(k) * gzps_vye(kd,iw)

        vze_upw(k,iw)  =                   vze(kd,iw) &
                       +  dxps_w(k) * gxps_vze(kd,iw) &
                       +  dyps_w(k) * gyps_vze(kd,iw) &
                       +  dzps_w(k) * gzps_vze(kd,iw)
     enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Finish MPI recv of THIL, VXE, VYE, and VZE gradient components

  if (iparallel == 1) then
     call mpi_recv_w(mrl, rvara1=gxps_the, rvara2=gyps_the, &
                          rvara3=gxps_vxe, rvara4=gyps_vxe, &
                          rvara5=gxps_vye, rvara6=gyps_vye, &
                          rvara7=gxps_vze, rvara8=gyps_vze  )
  endif

  ! Lateral boundary copy of THIL, VXE, VYE, and VZE gradient components

  call lbcopy_w(mrl, a1=gxps_the, a2=gyps_the, &
                     a3=gxps_vxe, a4=gyps_vxe, &
                     a5=gxps_vye, a6=gyps_vye, &
                     a7=gxps_vze, a8=gyps_vze  )

! EVALUATE UPWINDED THIL, VXE, VYE, AND VZE AT EACH V FACE FOR
! HORIZONTAL FLUX COMPUTATION.

  !$omp parallel private(vcf,dxps_v,dyps_v,dzps_v)
  !$omp do private(iv,iw1,iw2,k,iwd)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vcf(k) = 2.0 * vmcf(k,iv) / real( rho(k,iw1) + rho(k,iw2) )
     enddo

     call donorpointv(iv, dt, vcf, vxe, vye, vze, dxps_v, dyps_v, dzps_v)

     do k = lpv(iv), mza

        if (vcf(k) > 0.) then
           iwd = iw1
        else
           iwd = iw2
        endif

        thil_upv(k,iv) =                 thil(k,iwd) &
                       + dxps_v(k) * gxps_the(k,iwd) &
                       + dyps_v(k) * gyps_the(k,iwd) &
                       + dzps_v(k) * gzps_the(k,iwd)

        vxe_upv(k,iv)  =                  vxe(k,iwd) &
                       + dxps_v(k) * gxps_vxe(k,iwd) &
                       + dyps_v(k) * gyps_vxe(k,iwd) &
                       + dzps_v(k) * gzps_vxe(k,iwd)

        vye_upv(k,iv)  =                  vye(k,iwd) &
                       + dxps_v(k) * gxps_vye(k,iwd) &
                       + dyps_v(k) * gyps_vye(k,iwd) &
                       + dzps_v(k) * gzps_vye(k,iwd)

        vze_upv(k,iv)  =                  vze(k,iwd) &
                       + dxps_v(k) * gxps_vze(k,iwd) &
                       + dyps_v(k) * gyps_vze(k,iwd) &
                       + dzps_v(k) * gzps_vze(k,iwd)
     enddo
  enddo
  !$omp end do

! MAIN LOOP OVER W COLUMNS FOR UPDATING WM, WC, RHO, THIL, AND PRESS

  !$omp do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Prognose vertical velocity, density, thil, and diagnose pressure
     call prog_wrt_begs( iw, vmcf, wmsc,                       &
                         thil_upv, vxe_upv, vye_upv, vze_upv,  &
                         thil_upw, vxe_upw, vye_upw, vze_upw,  &
                         thilt_short, vmxet_short, vmyet_short, vmzet_short, &
                         vxesc, vyesc, vzesc )

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! MPI send of quantities needed for prog_v:
  ! PRESS, RHO, VMXET_VOLT, VMYET_VOLT, and VMZET_VOLT

  if (iparallel == 1) then
     call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_short, rvara2=vmyet_short, rvara3=vmzet_short )
  endif

  ! MPI recv and LBC copy of quantities needed for prog_v:
  ! PRESS, RHO, VMXET_VOLT, VMYET_VOLT, and VMZET_VOLT

  if (iparallel == 1) then
     call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=vmxet_short, rvara2=vmyet_short, rvara3=vmzet_short )
  endif

  call lbcopy_w(mrl, a1=vmxet_short, a2=vmyet_short, a3=vmzet_short,  &
                     d1=press,       d2=rho                           )

  ! Compute terms for rotational form of horizontal momentum

  if (rotational) then
     call prep_rotational(mrl)
  endif

  ! A good place to do divergence damping

  vmt_short = vmt

  call divh_damp(1, vmt_short)

  ! Main loop over v points to update vmc and vc

  !$omp parallel do private(iv)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     call prog_v_begs( iv, vmt_short(:,iv), vmxet_short, vmyet_short, vmzet_short )

  enddo
  !$omp end parallel do

  ! MPI SEND/RECV and LBC of VC and VMC

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif
  call lbcopy_v(mrl, vmc=vmc, vc=vc)

  ! UPDATE EARTH CARTESIAN VELOCITIES

  call diagvel_t3d(mrl)

end subroutine prog_wrtv_orig

!=========================================================================

subroutine prog_wrt_begs( iw, vmcf, wmsc,                       &
                          thil_upv, vxe_upv, vye_upv, vze_upv,  &
                          thil_upw, vxe_upw, vye_upw, vze_upw,  &
                          thilt_short, vmxet_short, vmyet_short, vmzet_short,   &
                          vxesc, vyesc, vzesc )

  use mem_ijtabs,  only: itab_w
  use mem_basic,   only: wmc, rho, thil, wc, press, vxe, vye, vze, &
                         alpha_press, pwfac
  use vel_t3d,     only: vxe1, vye1, vze1
  use misc_coms,   only: dtsm, deltax, nxp, initial, dn01d, th01d, nrk_scal
  use consts_coms, only: cpocv, omega2, pi1, pio180, r8
  use mem_grid,    only: mza, mva, mwa, lpv, lpw, lve2, arw, arv, &
                         vnx, vny, vnz, wnx, wny, wnz, wnxo2, wnyo2, wnzo2, &
                         dzim, volt, volti, volwi, glatw, glonw, gdzim, &
                         dzt_top, dzt_bot, zwgt_bot, gravm, &
                         zwgt_top8, zwgt_bot8, gdz_wgtm8, gdz_wgtp8
  use tridiag,     only: tridiffo
  use oname_coms,  only: nl
  use mem_turb,    only: akmodx, akhodx
  use mem_rayf,    only: dorayfw, rayf_cofw, krayfw_bot, &
                         dorayf, rayf_cof, krayf_bot

  implicit none

  integer, intent(in) :: iw

  real, intent(in)    :: vmcf(mza,mva)
  real, intent(inout) :: wmsc(mza,mwa)
  real, intent(in)    :: thil_upv(mza,mva)
  real, intent(in)    :: vxe_upv (mza,mva)
  real, intent(in)    :: vye_upv (mza,mva)
  real, intent(in)    :: vze_upv (mza,mva)

  real, intent(in)    :: thil_upw(mza,mwa)
  real, intent(in)    :: vxe_upw (mza,mwa)
  real, intent(in)    :: vye_upw (mza,mwa)
  real, intent(in)    :: vze_upw (mza,mwa)

  real, intent(inout) :: thilt_short(mza,mwa)
  real, intent(inout) :: vmxet_short(mza,mwa)
  real, intent(inout) :: vmyet_short(mza,mwa)
  real, intent(inout) :: vmzet_short(mza,mwa)

  real, intent(inout) :: vxesc(mza,mwa)
  real, intent(inout) :: vyesc(mza,mwa)
  real, intent(inout) :: vzesc(mza,mwa)

  integer :: jv, iv, iwn
  integer :: k, ka, kbv, ksw
  integer :: npoly, ktop

  real     :: dts, r4, dtov
  real(r8) :: dt8, v4, v8i

  real :: c6, c7, c8, c9, c10
  real :: dirv, vmarv
  real :: vmt1
  real :: rad0_swtc, rad_swtc, topo_swtc

  ! Vertical implicit scheme weighting parameters

  real, parameter :: fw = .55  ! wmc
  real, parameter :: fr = .80  ! rho
  real, parameter :: fp = .80  ! press

  real, parameter :: pc2 = fp * cpocv

  ! Automatic arrays

! real(r8) :: wmarw    (mza)
  real     :: del_wmarw
  real(r8) :: del_wmar8(mza)

  real(r8) :: hflux_rho(mza)
  real(r8) :: delex_rho(mza)
  real(r8) :: rhothil  (mza)
  real(r8) :: press_t  (mza)
  real(r8) :: delex_rhothil(mza)
  real(r8) :: hflux_thil(mza)
  real(r8) :: vflux_thil(mza)
  real(r8) :: thil_tend (mza)
  real(r8) :: rhoi, del_rhothil

  real :: delex_wm     (mza)
  real :: del_wm       (mza)
  real :: fwdel_wm     (mza)
  real :: r4i

  real :: hflux_vxe(mza)
  real :: hflux_vye(mza)
  real :: hflux_vze(mza)

  real :: hdiff_vxe(mza)
  real :: hdiff_vye(mza)
  real :: hdiff_vze(mza)

  real :: vflux_vxe(mza)
  real :: vflux_vye(mza)
  real :: vflux_vze(mza)

  real :: b1(mza),b2(mza),b3(mza),b5(mza),b6(mza),b10(mza)
  real :: b7(mza),b8(mza),b9(mza),b11(mza),b12(mza),b13(mza),b14(mza)
  real :: b20(mza),b21(mza),b22(mza),b23(mza),b24(mza),b25(mza),b26(mza)
  real :: b31(mza),b32(mza),b33(mza),b34(mza)

  real :: vmxet(mza), vmyet(mza), vmzet(mza)
  real :: vmxe (mza), vmye (mza), vmze (mza)
  real :: vxef (mza), vyef (mza), vzef (mza)

  real     :: wmarw
  real(r8) :: wmarw8(mza)

  ka = lpw(iw)

  dt8 = dtsm(itab_w(iw)%mrlw)  ! double precision DT
  dts = dt8                    ! single precision DT

  ! Set bottom & top vertical advective mass and heat fluxes to zero

  wmarw8(1:ka-1) = 0._r8
  wmarw8(mza)    = 0._r8

  vflux_thil(ka-1) = 0._r8
  vflux_thil(mza)  = 0._r8

  vflux_vxe(ka-1) = 0.
  vflux_vxe(mza)  = 0.

  vflux_vye(ka-1) = 0.
  vflux_vye(mza)  = 0.

  vflux_vze(ka-1) = 0.
  vflux_vze(mza)  = 0.

  ! Loop over W levels

  do k = ka,mza-1

     ! vertical mass flux for time level t

     wmarw8(k) = wmc(k,iw) * arw(k,iw)
     wmarw     = wmarw8(k)

     ! vertical fluxes

     vflux_thil(k) = wmarw * thil_upw(k,iw)
     vflux_vxe(k)  = wmarw * vxe_upw (k,iw)
     vflux_vye(k)  = wmarw * vye_upw (k,iw)
     vflux_vze(k)  = wmarw * vze_upw (k,iw)

  enddo

  ! Initialize horizontal advection arrays to zero

  hflux_rho (1:mza) = 0._r8
  hflux_thil(1:mza) = 0._r8
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

        vmarv = dirv * vmcf(k,iv) * arv(k,iv)

        ! Sum horizontal advection fluxes over V faces

        hflux_rho (k) = hflux_rho (k) + vmarv
        hflux_thil(k) = hflux_thil(k) + vmarv * thil_upv(k,iv)
        hflux_vxe (k) = hflux_vxe (k) + vmarv * vxe_upv (k,iv)
        hflux_vye (k) = hflux_vye (k) + vmarv * vye_upv (k,iv)
        hflux_vze (k) = hflux_vze (k) + vmarv * vze_upv (k,iv)

     enddo

  enddo

  do k = ka, mza
     r4 = real(rho (k,iw))
     v4 = real(volt(k,iw))

     ! Compute current T cell momentum and store in temp array

     vmxe(k) = vxe(k,iw) * r4
     vmye(k) = vye(k,iw) * r4
     vmze(k) = vze(k,iw) * r4

     ! Coriolis tendencies

     vmxet_short(k,iw) = vmxet_short(k,iw) + omega2 * v4 * vmye(k)
     vmyet_short(k,iw) = vmyet_short(k,iw) - omega2 * v4 * vmxe(k)

     ! Explicit momentum tendencies from advective transport
     ! (weighted by T cell volume)

     vmxet(k) = vmxet_short(k,iw) + hflux_vxe(k) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet(k) = vmyet_short(k,iw) + hflux_vye(k) + vflux_vye(k-1) - vflux_vye(k)
     vmzet(k) = vmzet_short(k,iw) + hflux_vze(k) + vflux_vze(k-1) - vflux_vze(k)

     ! Explicit THIL tendencies from advective transport
     ! (weighted by T cell volume)

     thil_tend(k) = thilt_short(k,iw) + hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)

  enddo

  ! RAYLEIGH FRICTION ON THIL - only for horizontally homogeneous initialization

  if (dorayf .and. initial == 1) then

     do k = krayf_bot, mza
        v4 = real(volt(k,iw))

        thil_tend(k) = thil_tend(k) &
                        + real( rayf_cof(k) * dn01d(k) * (th01d(k) - thil(k,iw)) * v4, r8)
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
     v8i = real(volti(k,iw),r8)

     ! Explicit density tendency

     delex_rho(k) = dt8 * v8i * (hflux_rho(k) + wmarw8(k-1) - wmarw8(k))

     ! Explicit density-thil tendency

     delex_rhothil(k) = dt8 * v8i * thil_tend(k)

     ! RHOTHIL(t) and PRESS(t)

     rhothil(k) = rho(k,iw) * thil(k,iw)
     press_t(k) = alpha_press(k,iw) * real(rhothil(k)) ** cpocv

  enddo

  c6  = dts * .50 * fw
  c7  = dts * .25 * fw
  c8  = dts * pc2
  c9  =-dts * fr
  c10 = dts * fw

  ! Vertical loop over W levels

  do k = ka,mza-1

     ! Explicit vertical momentum tendency

     delex_wm(k) = dts * ( real( press_t(k) - press_t(k+1) ) * pwfac(k,iw) &

          - real( gdz_wgtm8(k) * rho(k,iw) + gdz_wgtp8(k) * rho(k+1,iw) ) &

          ! ADD volume-weighted arrays VMXET_VOLT, VMYET_VOLT, VMZET_VOLT
          ! over both T levels and divide by combined VOLT of those T levels

          + ( wnx(iw) * (vmxet(k) + vmxet(k+1)) &
            + wny(iw) * (vmyet(k) + vmyet(k+1)) &
            + wnz(iw) * (vmzet(k) + vmzet(k+1)) ) * volwi(k,iw) )

  enddo

  ! Rayleigh friction on W

  b20(:) = 1.0

  if (dorayfw) then
     do k = krayfw_bot, mza-1
        delex_wm(k) = delex_wm(k) - dts * rayf_cofw(k) * wmc(k,iw)
        b20     (k) = b20(k) + c10 * rayf_cofw(k)
     enddo
  endif

  ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions

  ! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
  ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-1)
  ! rayfx = 0.   ! Default: no extra RAYF
  ! do k = ka, mza-1
  !    wmt_rayf(k) = - max(rayf_cofw(k),rayfx) * wmc(k,iw)
  ! enddo

  ! Fill matrix coefficients for implicit update of WM

  do k = ka, mza-1
     b2(k)  = thil(k,iw) + thil(k+1,iw)       ! W pts
  enddo

  b2(ka-1)= 2. * thil(ka,iw)
  b2(mza) = 2. * thil(mza,iw)

  do k = ka, mza
     b1(k)  = wc(k,iw) + wc(k-1,iw)    ! T pts
     b5(k)  = real(press_t(k)) / real(rhothil(k))  ! T pts
     b6(k)  = c6  * volti(k,iw)        ! T pts
     b10(k) = c10 * volti(k,iw)        ! T pts
  enddo

  do k = ka,mza-1

     b7(k)  = c6     * volwi(k,iw)  ! W pts
     b8(k)  = c8     * pwfac(k,iw)  ! W pts
     b9(k)  = c9     * gdzim(k)     ! W pts

     b11(k) = b8(k)  * b5(k)        ! W pts
     b12(k) = b8(k)  * b5(k+1)      ! W pts

     b13(k) = b9(k)  * dzt_top(k)   ! W pts
     b14(k) = b9(k)  * dzt_bot(k+1) ! W pts

     b21(k) = b7(k)  * b1(k)    ! W pts
     b22(k) = b7(k)  * b1(k+1)  ! W pts

     b23(k) = b11(k) * b6(k)    ! W pts
     b24(k) = b12(k) * b6(k+1)  ! W pts

     b25(k) = b13(k) * b10(k)   ! W pts
     b26(k) = b14(k) * b10(k+1) ! W pts

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

     del_wmar8(k)  = fwdel_wm(k) * arw(k,iw)
     del_wmarw     = del_wmar8(k)

     vflux_thil(k) = del_wmarw * thil_upw(k,iw)
     vflux_vxe (k) = del_wmarw * vxe_upw (k,iw)
     vflux_vye (k) = del_wmarw * vye_upw (k,iw)
     vflux_vze (k) = del_wmarw * vze_upw (k,iw)

  enddo

  ! Set mass-flux-change array to 0 at top & bottom (other flux arrays already 0)

  del_wmar8(ka-1) = 0._r8
  del_wmar8(mza)  = 0._r8

  ! Vertical loop over T points

  do k = ka,mza

     ! Change of rho from (t) to (t+1)

     rho(k,iw) = rho(k,iw) + delex_rho(k) &
               + dt8 * real(volti(k,iw),r8) * (del_wmar8(k-1) - del_wmar8(k))

     rhoi = 1._r8 / rho(k,iw)

     ! Change of rhothil from (t) to (t+1)

     del_rhothil = delex_rhothil(k) &
                 + dt8 * real(volti(k,iw),r8) * (vflux_thil(k-1) - vflux_thil(k))

     ! Update pressure from (t) to (t+tp)

     press(k,iw) = press_t(k) + pc2 * b5(k) * del_rhothil

     ! Update thil from (t) to (t+1)

     thil(k,iw) = (rhothil(k) + del_rhothil) * rhoi

     ! Update volume-weighted momentum tendencies at T due to implicit flux correction

     vmxet(k) = vmxet(k) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet(k) = vmyet(k) + vflux_vye(k-1) - vflux_vye(k)
     vmzet(k) = vmzet(k) + vflux_vze(k-1) - vflux_vze(k)

  enddo

  if (nrk_scal == 1) then
     ktop = mza
  elseif (nl%icut_vel == 1) then
     ktop = ka + lve2(iw) - 1
  else
     ktop = ka - 1
  endif

  ! Estimate velocity in T cells at (t+1) by prognostic method

  do k = ka, ktop
     dtov = dts * volti(k,iw)
     r4i  = 1.0 / real(rho(k,iw))

     vxef(k) = (vmxe(k) + dtov * vmxet(k)) * r4i
     vyef(k) = (vmye(k) + dtov * vmyet(k)) * r4i
     vzef(k) = (vmze(k) + dtov * vmzet(k)) * r4i
  enddo

  ! Add half-forward T cell velocity to scalar earth-cartesian velocity arrays

  if (nrk_scal == 1) then
     do k = ka, mza
        vxesc(k,iw) = vxesc(k,iw) + .5 * (vxe(k,iw) + vxef(k))
        vyesc(k,iw) = vyesc(k,iw) + .5 * (vye(k,iw) + vyef(k))
        vzesc(k,iw) = vzesc(k,iw) + .5 * (vze(k,iw) + vzef(k))
     enddo
  endif

  ! Save velocity in T cells at (t+1) for prognostic method

  if (nl%icut_vel == 1) then
     do ksw = 1, lve2(iw)
        k = ksw + ka - 1

        vxe1(ksw,iw) = vxef(k)
        vye1(ksw,iw) = vyef(k)
        vze1(ksw,iw) = vzef(k)
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

  do k = ka, mza-1

     ! Update WMC and WC to future value (at t+1) due to change in WM (del_wm)

     wmc(k,iw) = wmc(k,iw) + del_wm(k)
     wc (k,iw) = wmc(k,iw) / real(zwgt_bot8(k) * rho(k,iw) + zwgt_top8(k) * rho(k+1,iw))

  enddo

  ! Set top & bottom values of WC

!  wc(ka-1,iw) = wc(ka,iw)

!  wc(1:ka-2,iw) = 0.
!  wc(mza   ,iw) = 0.

  ! Bottom boundary for THIL

  do k = 1, ka-1
     thil(k,iw) = thil(ka,iw)
  enddo

  if (.not. rotational) then
     do k = ka, mza
        vmxet_short(k,iw) = vmxet(k)
        vmyet_short(k,iw) = vmyet(k)
        vmzet_short(k,iw) = vmzet(k)
     enddo
  endif

end subroutine prog_wrt_begs

!============================================================================

subroutine prog_v_begs( iv, vmt_short, vmxet_short, vmyet_short, vmzet_short )

  use mem_ijtabs,  only: itab_v
  use mem_basic,   only: vc, press, vmc, rho, vxe, vye, vze, pvfac
  use mem_tend,    only: vmt
  use misc_coms,   only: dtsm, initial, mdomain, deltax, nxp
  use consts_coms, only: eradi, gravo2
  use mem_grid,    only: mza, mwa, lpv, volt, xev, yev, zev, &
                         unx, uny, unz, vnx, vny, vnz, vnxo2, vnyo2, vnzo2, &
                         dniu, dniv, volvi
  use mem_rayf!,    only: dorayf, rayf_cof, vc03d, dn03d, krayf_bot
  use oname_coms,  only: nl

  implicit none

  integer, intent(in) :: iv

  real, intent(in) :: vmt_short(mza)
  real, intent(in) :: vmxet_short(mza,mwa)
  real, intent(in) :: vmyet_short(mza,mwa)
  real, intent(in) :: vmzet_short(mza,mwa)

  integer :: k, kb
  integer :: iw1,iw2,im1,im2

  real :: dts
  real :: fracx, rayfx
  real :: vx, vy, vz, uc, watv, tke1, tke2, vortp_v, dtso2dnu, dtso2dnv
  real :: pgf(mza)

  ! Extract neighbor indices and coefficients for this point in the U stencil

  iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
  im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)

  dts = dtsm(itab_v(iv)%mrlv)

  kb = lpv(iv)

  ! Vertical loop over V points to compute pressure gradient force.

  do k = kb, mza
     pgf(k) = real( press(k,iw1) - press(k,iw2) ) * pvfac(k,iv)
  enddo

  ! For shallow water test cases 2 & 5, rho & press are interpreted
  ! as water depth & height.

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     do k = kb, mza
        pgf(k) = pgf(k) * real( rho(k,iw1) + rho(k,iw2) ) * gravo2
     enddo
  endif

  ! Update VMC

  if (.not. rotational) then

     do k = kb, mza
        vmc(k,iv) = vmc(k,iv) + dts * ( vmt_short(k) + pgf(k)             &
                  + ( vnx(iv) * (vmxet_short(k,iw1) + vmxet_short(k,iw2)) &
                    + vny(iv) * (vmyet_short(k,iw1) + vmyet_short(k,iw2)) &
                    + vnz(iv) * (vmzet_short(k,iw1) + vmzet_short(k,iw2)) ) * volvi(k,iv) )

        vc(k,iv) = 2.0 * vmc(k,iv) / real(rho(k,iw1) + rho(k,iw2))
     enddo

  else  ! If using rotational form:

     dtso2dnu = 0.5 * dts * dniu(iv)
     dtso2dnv = 0.5 * dts * dniv(iv)

     ! Vertical loop over V levels

     do k = kb, mza
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

        vmc(k,iv) = vmc(k,iv) + dts * ( vmt_short(k) + pgf(k)             &

                  + ( vnx(iv) * (vmxet_short(k,iw1) + vmxet_short(k,iw2)) &
                    + vny(iv) * (vmyet_short(k,iw1) + vmyet_short(k,iw2)) &
                    + vnz(iv) * (vmzet_short(k,iw1) + vmzet_short(k,iw2)) ) * volvi(k,iv) &

                  + .5 * real(rho(k,iw1) + rho(k,iw2))                    &
                    * (dniv(iv) * (tke1 - tke2)                           &
                  + uc * vortp_v                                          &
                  - watv * .5 * (vortn(k-1,iv) + vortn(k,iv))))

        vc(k,iv) = 2.0 * vmc(k,iv) / real(rho(k,iw1) + rho(k,iw2))
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


end module wrtv_orig
