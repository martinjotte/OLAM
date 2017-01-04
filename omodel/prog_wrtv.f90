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
subroutine prog_wrtv(vmsc,wmsc,vxesc,vyesc,vzesc,alpha_press,rhot)

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, istp, itab_w, jtab_m, itab_m, &
                        mrl_begl, mrl_begs, mrl_ends, &
                        jtm_vadj, jtv_prog, jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp
use mem_basic,    only: rho, thil, theta, wc, wmc, press, vmp, vmc, vp, vc, &
                        vxe, vye, vze, vxe2, vye2, vze2, &
                        strict_wvt_donorpoint
use mem_grid,     only: mza, mma, mva, mwa, lpm, lpv, lpw, &
                        dzim, zfact, zfacit, zfacim, dnv, dniv, dnu, &
                        arm0, vnx, vny, vnz, wnxo2, wnyo2, wnzo2
use mem_tend,     only: thilt, vmxet, vmyet, vmzet, sh_wt
use misc_coms,    only: iparallel, time8, dtlm, rinit, initial, dn01d, th01d, &
                        deltax, nxp, mdomain
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
use oplot_coms,   only: op
use obnd,         only: lbcopy_m, lbcopy_w
use mem_rayf,     only: dorayf, dorayfw, dorayfdiv, krayfdiv_bot, &
                        rayf_cof, rayf_cofw, krayf_bot, krayfw_bot

use vel_t3d,      only: vel_t3d_hex
use oname_coms,   only: nl
use mem_adv,      only: dxps_w, dyps_w, dzps_w, &
                        dxyps_w, dxxps_w, dyyps_w, dzzps_w, &
                        dxps_v, dyps_v, dzps_v, &
                        dxyps_v, dxxps_v, dyyps_v, dzzps_v, &
                        gxps_scp, gyps_scp, gzps_scp, &
                        gxyps_scp, gxxps_scp, gyyps_scp, gzzps_scp, &
                        gxps_vxe, gyps_vxe, gzps_vxe, &
                        gxyps_vxe, gxxps_vxe, gyyps_vxe, gzzps_vxe, &
                        gxps_vye, gyps_vye, gzps_vye, &
                        gxyps_vye, gxxps_vye, gyyps_vye, gzzps_vye, &
                        gxps_vze, gyps_vze, gzps_vze, &
                        gxyps_vze, gxxps_vze, gyyps_vze, gzzps_vze

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

real, intent(in)    :: alpha_press(mza,mwa)
real, intent(inout) :: rhot       (mza,mwa)

integer :: j, iv, k, ka, kb, mrl, kbv, kd
integer :: iw, iw1, iw2

! automatic arrays

integer :: iwdepv(mza,mva) ! donor cell IW index for V face
integer :: kdepw(mza,mwa)  ! donor cell K index for W face

real :: vmcf(mza,mva) ! Time-extrapolated VMC

real, allocatable :: vcf (:,:) ! Time-extrapolated VC
real, allocatable :: vxef(:,:) ! Time-extrapolated XE velocity component at T point
real, allocatable :: vyef(:,:) ! Time-extrapolated YE velocity component at T point
real, allocatable :: vzef(:,:) ! Time-extrapolated ZE velocity component at T point

real, allocatable, save :: vortp(:,:)
real :: vortn  (mza,mva)
real :: vortp_t(mza,mwa)
real :: div2d  (mza,mwa)

! Volume-weighted T-cell momentum component tendencies from advective and
! turbulent transport; evaluated and applied on acoustic timestep

real :: vmxet_volt(mza,mwa)
real :: vmyet_volt(mza,mwa)
real :: vmzet_volt(mza,mwa)

real :: vmx_cor(mza,mwa)
real :: vmy_cor(mza,mwa)

real :: thil_upv(mza,mva) ! Upstreamed THIL at each V interface
real :: vxe_upv (mza,mva) ! Upstreamed VXE  at each V interface
real :: vye_upv (mza,mva) ! Upstreamed VYE  at each V interface
real :: vze_upv (mza,mva) ! Upstreamed VZE  at each V interface

real :: thil_upw(mza,mwa) ! Upstreamed THIL at each W level
real :: vxe_upw (mza,mwa) ! Upstreamed VXE  at each W level
real :: vye_upw (mza,mwa) ! Upstreamed VYE  at each W level
real :: vze_upw (mza,mwa) ! Upstreamed VZE  at each W level

integer :: im,npoly,jv,iwd
integer :: jm
real :: arm0i
real    :: fracx, rayfx

logical, save :: firstime = .true.

logical :: rotational

if (firstime) then
   firstime = .false.
   allocate(vortp(mza,mma))
   vortp = rinit
endif

rotational = .false.
!if (nl%expnme(1:1) == 'R') rotational = .true.

! Half-forward velocities for computing donor points:

if (strict_wvt_donorpoint) then
   allocate(vcf (mza,mva)) ; vcf  = 0.0
   allocate(vxef(mza,mwa)) ; vxef = 0.0
   allocate(vyef(mza,mwa)) ; vyef = 0.0
   allocate(vzef(mza,mwa)) ; vzef = 0.0
endif

vmcf(:,1) = 0.

! Add remaining contributions to long timestep tendencies

mrl = mrl_begl(istp)
if (mrl > 0) then

   ! Horizontal loop over W columns for BEGL

   !----------------------------------------------------------------------
   !$omp parallel do private(iw,ka,k) 
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   !----------------------------------------------------------------------

      ka = lpw(iw)

      ! Vertical loop over T levels

      do k = ka, mza

         ! Include moisture changes in total density tendency

         rhot(k,iw) = rhot(k,iw) + sh_wt(k,iw)

         ! Apply density tendency to momentum tendencies to conserve momentum

         vmxet(k,iw) = vmxet(k,iw) + vxe(k,iw) * rhot(k,iw)
         vmyet(k,iw) = vmyet(k,iw) + vye(k,iw) * rhot(k,iw)
         vmzet(k,iw) = vmzet(k,iw) + vze(k,iw) * rhot(k,iw)

      enddo

      ! RAYLEIGH FRICTION ON THIL - only for horizontally homogeneous initialization

      if (dorayf .and. initial == 1) then

         do k = krayf_bot, mza
            thilt(k,iw) = thilt(k,iw) &
                        + rayf_cof(k) * dn01d(k) * (th01d(k) - theta(k,iw))
         enddo

         ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions
         ! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
         ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
         ! rayfx = 0.   ! Default: no extra RAYF
         ! do k = ka, mza
         !    thilt(k,iw) = thilt(k,iw) + max(rayf_cof(k),rayfx)  & 
         !                * dn01d(k) * (th01d(k) - theta(k,iw))
         ! enddo

      endif ! (dorayf .and. initial == 1)

   enddo
   !$omp end parallel do

! MPI SEND/RECV of VMXET, VMYET, VMZET

   if (iparallel == 1) then
      call mpi_send_w(mrl, rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)
      call mpi_recv_w(mrl, rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)
   endif

   call lbcopy_w(mrl, a1=vmxet, a2=vmyet, a3=vmzet)

endif ! mrl = mrl_begl(istp) > 0

! Horizontal loop over M/P columns for computing vertical vorticity.  Always do this
! on BEGL timestep; if using rotational method, also do this on BEGS timestep

if (rotational) then
   mrl = mrl_begs(istp)
else
   mrl = mrl_begl(istp)
endif

if (mrl > 0) then
!----------------------------------------------------------------------
   !$omp parallel do private(im,npoly,kb,jv,iv,k,arm0i)
   do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)
!----------------------------------------------------------------------

      npoly = itab_m(im)%npoly

      kb = lpm(im)

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

   ! Parallel send/recv of vortp

   if (iparallel == 1) then
      call mpi_send_m(mrl, rvara1=vortp)
      call mpi_recv_m(mrl, rvara1=vortp)
   endif
   call lbcopy_m(mrl, a1=vortp)

endif ! mrl > 0

! If using rotational form for prognosing horizontal momentum,
! evaluate vortp_t and vortn

if (rotational) then
   mrl = mrl_begs(istp)
   if (mrl > 0) then

! Horizontal loop over W/T columns

   !----------------------------------------------------------------------
   !$omp parallel do private(iw,npoly,kb,k,jm,im)
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   !----------------------------------------------------------------------

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

!----------------------------------------------------------------------
   !$omp parallel do private(iv,iw1,iw2,k,kb) 
   do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
   !----------------------------------------------------------------------

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
   endif ! mrl > 0

endif ! rotational

! Horizontal loop over V/N columns

!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(k) 
do iv = 2,mva
!----------------------------------------------------------------------

! Extrapolate VM to time T + 1/2; update VMP

   do k = lpv(iv),mza
      vmcf(k,iv) = (1.5 * vmc(k,iv) - 0.5 * vmp(k,iv))
      vmsc(k,iv) = vmsc(k,iv) + vmcf(k,iv)
      vmp (k,iv) = vmc (k,iv)
   enddo

! Extrapolate V to time T + 1/2; update VP

   if (strict_wvt_donorpoint) then
      do k = 1,mza
         vcf (k,iv) = (1.5 * vc (k,iv) - 0.5 * vp (k,iv))
         vp  (k,iv) = vc  (k,iv)
      enddo
   endif

enddo
!$omp end parallel do 
endif

if (strict_wvt_donorpoint) then

! Compute donor point locations using half-forward velocities.

   mrl = mrl_begs(istp)

   if (mrl > 0) then

      ! Diagnose 3D velocity at T points for BEGS using half-future vcf
      ! and current wc (assume not necessary to extrapolate vxe2,vye2,vze2)

      call vel_t3d_hex(mrl, vcf, wc, vxef, vyef, vzef, vxe2, vye2, vze2)

      ! MPI send of VXE, VYE, VZE

      if (iparallel == 1) then
         call mpi_send_w(mrl, rvara1=vxef, rvara2=vyef, rvara3=vzef)
      endif

      ! Diagnose advective donor point locations for all primary W faces
      ! No parallel communication is necessary to compute this

      if (nl%adv_order <= 2) then
         call donorpointw  (1, mrl, wc, vxef, vyef, vzef, kdepw)
      else
         call donorpointw_3(1, mrl, wc, vxef, vyef, vzef, kdepw)
      endif

      ! Finish MPI recv of VXE, VYE, VZE and do a LBC copy

      if (iparallel == 1) then
         call mpi_recv_w(mrl, rvara1=vxef, rvara2=vyef, rvara3=vzef)
      endif

      call lbcopy_w(mrl, a1=vxef, a2=vyef, a3=vzef)

      ! Diagnose advective donor point locations for the V faces surrounding all
      ! primary W points. Communication of velocities must have been completed

      if (nl%adv_order <= 2) then
         call donorpointv  (1, mrl, vcf, vxef, vyef, vzef, iwdepv)
      else
         call donorpointv_3(1, mrl, vcf, vxef, vyef, vzef, iwdepv)
      endif

   endif

else

! Compute donor point locations using current velocities.

   mrl = mrl_begs(istp)
   if (mrl > 0) then

      if (nl%adv_order <= 2) then
         call donorpointw  (1, mrl, wc, vxe, vye, vze, kdepw)
         call donorpointv  (1, mrl, vc, vxe, vye, vze, iwdepv)
      else
         call donorpointw_3(1, mrl, wc, vxe, vye, vze, kdepw)
         call donorpointv_3(1, mrl, vc, vxe, vye, vze, iwdepv)
      endif

   endif

endif  ! strict_wvt_donorpoint

! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
! (Example of how to plot "external" field; one not available in module memory)
!
!if (mod(time8,op%frqplt) < dtlm(1) .and. istp == 900) then
!
!   allocate (op%extfld(mza,mwa))
!   op%extfld(:,:) = vxe(:,:)
!   op%extfldname = 'VXE'
!   call plot_fields(11)
!   deallocate (op%extfld)
!
!endif
! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

mrl = mrl_begs(istp)
if (mrl > 0) then

! Compute horizontal divergence if we are damping it at the model top

   if (dorayfdiv) then
      
      !$omp parallel do private(iw,jv,iv,k) 
      do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
         
         div2d(:,iw) = 0.0
         do jv = 1, itab_w(iw)%npoly
            iv = itab_w(iw)%iv(jv)

            do k = krayfdiv_bot, mza
               div2d(k,iw) = div2d(k,iw) - itab_w(iw)%dirv(jv) * vmp(k,iv) * dnu(iv)
            enddo
         enddo
      enddo
      !$omp end parallel do

      if (iparallel == 1) then
         call mpi_send_w(mrl, rvara1=div2d)
      endif
      call lbcopy_w(mrl, a1=div2d)

   endif ! dorayfdiv

! Evaluate horizontal gradients of THIL, VXE, VYE, and VZE for BEGS

  if (nl%adv_order <= 2) then
     call grad_t2d  (mrl, thil, gxps_scp,  gyps_scp)
     call grad_t2d  (mrl, vxe,  gxps_vxe,  gyps_vxe)
     call grad_t2d  (mrl, vye,  gxps_vye,  gyps_vye)
     call grad_t2d  (mrl, vze,  gxps_vze,  gyps_vze)
  else
     call grad_t2d_3(mrl, thil, gxps_scp,  gyps_scp,  &
                                gxxps_scp, gxyps_scp, gyyps_scp)
     call grad_t2d_3(mrl, vxe,  gxps_vxe,  gyps_vxe,  &
                                gxxps_vxe, gxyps_vxe, gyyps_vxe)
     call grad_t2d_3(mrl, vye,  gxps_vye,  gyps_vye,  &
                                gxxps_vye, gxyps_vye, gyyps_vye)
     call grad_t2d_3(mrl, vze,  gxps_vze,  gyps_vze,  &
                                gxxps_vze, gxyps_vze, gyyps_vze)
  endif

! Finish MPI RECV of DIV2D, and MPI SEND of THIL, VXE, VYE, and VZE 
! gradient components (12 in all)

   if (iparallel == 1) then

      if (dorayfdiv) then
         call mpi_recv_w(mrl, rvara1=div2d)
      endif

      if (nl%adv_order <= 2) then
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

   endif

! Evaluate vertical gradient of THIL, VXE, VYE, and VZE for BEGS

   if (nl%adv_order <= 2) then
      call grad_z  (mrl, thil, gzps_scp)
      call grad_z  (mrl, vxe,  gzps_vxe)
      call grad_z  (mrl, vye,  gzps_vye)
      call grad_z  (mrl, vze,  gzps_vze)
   else
      call grad_z_3(mrl, thil, gzps_scp, gzzps_scp)
      call grad_z_3(mrl, vxe,  gzps_vxe, gzzps_vxe)
      call grad_z_3(mrl, vye,  gzps_vye, gzzps_vye)
      call grad_z_3(mrl, vze,  gzps_vze, gzzps_vze)
   endif

!  Horizontal loop over all primary W columns

   !$omp parallel do private(iw,k,kd) 
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
      
      do k = lpw(iw), mza

         kd = kdepw(k,iw)

         if (nl%adv_order <= 2) then

            thil_upw(k,iw) = thil(kd,iw)                    &
                           + dxps_w(k,iw) * gxps_scp(kd,iw) &
                           + dyps_w(k,iw) * gyps_scp(kd,iw) &
                           + dzps_w(k,iw) * gzps_scp(kd,iw)

            vxe_upw(k,iw)  = vxe(kd,iw)                     &
                           + dxps_w(k,iw) * gxps_vxe(kd,iw) &
                           + dyps_w(k,iw) * gyps_vxe(kd,iw) &
                           + dzps_w(k,iw) * gzps_vxe(kd,iw)

            vye_upw(k,iw)  = vye(kd,iw)                     &
                           + dxps_w(k,iw) * gxps_vye(kd,iw) &
                           + dyps_w(k,iw) * gyps_vye(kd,iw) &
                           + dzps_w(k,iw) * gzps_vye(kd,iw)

            vze_upw(k,iw)  = vze(kd,iw)                     &
                           + dxps_w(k,iw) * gxps_vze(kd,iw) &
                           + dyps_w(k,iw) * gyps_vze(kd,iw) &
                           + dzps_w(k,iw) * gzps_vze(kd,iw)
         else

            thil_upw(k,iw) = thil(kd,iw)                                                       &
                           + dxps_w(k,iw) * gxps_scp(kd,iw) + dxxps_w(k,iw) * gxxps_scp(kd,iw) &
                                                            + dxyps_w(k,iw) * gxyps_scp(kd,iw) &
                           + dyps_w(k,iw) * gyps_scp(kd,iw) + dyyps_w(k,iw) * gyyps_scp(kd,iw) &
                           + dzps_w(k,iw) * gzps_scp(kd,iw) + dzzps_w(k,iw) * gzzps_scp(kd,iw)

            vxe_upw(k,iw)  = vxe(kd,iw)                                                        &
                           + dxps_w(k,iw) * gxps_vxe(kd,iw) + dxxps_w(k,iw) * gxxps_vxe(kd,iw) &
                           + dyps_w(k,iw) * gyps_vxe(kd,iw) + dxyps_w(k,iw) * gxyps_vxe(kd,iw) &
                           + dzps_w(k,iw) * gzps_vxe(kd,iw) + dzzps_w(k,iw) * gzzps_vxe(kd,iw)

            vye_upw(k,iw)  = vye(kd,iw)                                                        &
                           + dxps_w(k,iw) * gxps_vye(kd,iw) + dxxps_w(k,iw) * gxxps_vye(kd,iw) &
                           + dyps_w(k,iw) * gyps_vye(kd,iw) + dxyps_w(k,iw) * gxyps_vye(kd,iw) &
                           + dzps_w(k,iw) * gzps_vye(kd,iw) + dzzps_w(k,iw) * gzzps_vye(kd,iw)

            vze_upw(k,iw)  = vze(kd,iw)                                                        &
                           + dxps_w(k,iw) * gxps_vze(kd,iw) + dxxps_w(k,iw) * gxxps_vze(kd,iw) &
                           + dyps_w(k,iw) * gyps_vze(kd,iw) + dxyps_w(k,iw) * gxyps_vze(kd,iw) &
                           + dzps_w(k,iw) * gzps_vze(kd,iw) + dzzps_w(k,iw) * gzzps_vze(kd,iw)
         endif

      enddo

   enddo
   !$omp end parallel do

   if (iparallel == 1) then
      if (nl%adv_order <= 2) then
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

   if (nl%adv_order <= 2) then
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

!  Horizontal loop over V

!----------------------------------------------------------------------
   !$omp parallel do private(iv,iwd,kbv,k) 
   do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)
!----------------------------------------------------------------------

      kbv = lpv(iv)
      do  k = kbv, mza

         iwd = iwdepv(k,iv)

         if (nl%adv_order <= 2) then

            thil_upv(k,iv) = thil(k,iwd)                    &
                           + dxps_v(k,iv) * gxps_scp(k,iwd) &
                           + dyps_v(k,iv) * gyps_scp(k,iwd) &
                           + dzps_v(k,iv) * gzps_scp(k,iwd)

            vxe_upv(k,iv)  = vxe(k,iwd)                     &
                           + dxps_v(k,iv) * gxps_vxe(k,iwd) &
                           + dyps_v(k,iv) * gyps_vxe(k,iwd) &
                           + dzps_v(k,iv) * gzps_vxe(k,iwd)

            vye_upv(k,iv)  = vye(k,iwd)                     &
                           + dxps_v(k,iv) * gxps_vye(k,iwd) &
                           + dyps_v(k,iv) * gyps_vye(k,iwd) &
                           + dzps_v(k,iv) * gzps_vye(k,iwd)

            vze_upv(k,iv)  = vze(k,iwd)                     &
                           + dxps_v(k,iv) * gxps_vze(k,iwd) &
                           + dyps_v(k,iv) * gyps_vze(k,iwd) &
                           + dzps_v(k,iv) * gzps_vze(k,iwd)
         else

            thil_upv(k,iv) = thil(k,iwd)                                                       &
                           + dxps_v(k,iv) * gxps_scp(k,iwd) + dxxps_v(k,iv) * gxxps_scp(k,iwd) &
                                                            + dxyps_v(k,iv) * gxyps_scp(k,iwd) &
                           + dyps_v(k,iv) * gyps_scp(k,iwd) + dyyps_v(k,iv) * gyyps_scp(k,iwd) &
                           + dzps_v(k,iv) * gzps_scp(k,iwd) + dzzps_v(k,iv) * gzzps_scp(k,iwd)

            vxe_upv(k,iv) = vxe(k,iwd)                                                         &
                           + dxps_v(k,iv) * gxps_vxe(k,iwd) + dxxps_v(k,iv) * gxxps_vxe(k,iwd) &
                                                            + dxyps_v(k,iv) * gxyps_vxe(k,iwd) &
                           + dyps_v(k,iv) * gyps_vxe(k,iwd) + dyyps_v(k,iv) * gyyps_vxe(k,iwd) &
                           + dzps_v(k,iv) * gzps_vxe(k,iwd) + dzzps_v(k,iv) * gzzps_vxe(k,iwd)

            vye_upv(k,iv) = vye(k,iwd)                                                         &
                           + dxps_v(k,iv) * gxps_vye(k,iwd) + dxxps_v(k,iv) * gxxps_vye(k,iwd) &
                                                            + dxyps_v(k,iv) * gxyps_vye(k,iwd) &
                           + dyps_v(k,iv) * gyps_vye(k,iwd) + dyyps_v(k,iv) * gyyps_vye(k,iwd) &
                           + dzps_v(k,iv) * gzps_vye(k,iwd) + dzzps_v(k,iv) * gzzps_vye(k,iwd)

            vze_upv(k,iv) = vze(k,iwd)                                                         &
                           + dxps_v(k,iv) * gxps_vze(k,iwd) + dxxps_v(k,iv) * gxxps_vze(k,iwd) &
                                                            + dxyps_v(k,iv) * gxyps_vze(k,iwd) &
                           + dyps_v(k,iv) * gyps_vze(k,iwd) + dyyps_v(k,iv) * gyyps_vze(k,iwd) &
                           + dzps_v(k,iv) * gzps_vze(k,iwd) + dzzps_v(k,iv) * gzzps_vze(k,iwd)
         endif

      enddo

   enddo
   !$omp end parallel do

endif

! Main loop over W columns for updating WM, WC, RHO, THIL, and PRESS

!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Prognose vertical velocity, density, thil, and diagnose pressure

   call prog_wrt_begs( iw, vmcf, wmsc, alpha_press, rhot,    &
                       thil_upv, vxe_upv, vye_upv, vze_upv,  &
                       thil_upw, vxe_upw, vye_upw, vze_upw,  &
                       vmxet_volt, vmyet_volt, vmzet_volt,   &
                       vxesc, vyesc, vzesc, vmx_cor, vmy_cor )

enddo
!$omp end parallel do 
endif

! MPI SEND/RECV and LBC copy of quantities needed for prog_v: 
! PRESS, RHO, VMXET, VMYET, and VMZET

if (iparallel == 1) then

   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=vmxet, rvara2=vmyet, rvara3=vmzet, &
                   rvara4=vmx_cor, rvara5=vmy_cor)

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=vmxet, rvara2=vmyet, rvara3=vmzet, &
                   rvara4=vmx_cor, rvara5=vmy_cor)

endif

call lbcopy_w(mrl, a1=vmxet, a2=vmyet, a3=vmzet, &
              a4=vmx_cor, a5=vmy_cor, d1=press, d2=rho)

! Horizontal loop over V points to update VMC

!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iv) 
do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!----------------------------------------------------------------------

   call prog_v_begs(iv,vmxet_volt,vmyet_volt,vmzet_volt,div2d,vmx_cor,vmy_cor, &
                    vortp,vortn,vortp_t,rotational)

enddo
!$omp end parallel do 
endif

if (strict_wvt_donorpoint) then
   deallocate( vcf)
   deallocate(vxef)
   deallocate(vyef)
   deallocate(vzef)
endif

end subroutine prog_wrtv

!=========================================================================

subroutine prog_wrt_begs( iw, vmcf, wmsc, alpha_press, rhot,    &
                          thil_upv, vxe_upv, vye_upv, vze_upv,  &
                          thil_upw, vxe_upw, vye_upw, vze_upw,  &  
                          vmxet_volt, vmyet_volt, vmzet_volt,   &
                          vxesc, vyesc, vzesc, vmx_cor, vmy_cor )

use mem_tend,    only: thilt, vmxet, vmyet, vmzet
use mem_ijtabs,  only: itab_w
use mem_basic,   only: wmc, rho, thil, wc, press, &
                       vxe, vye, vze, vxe2, vye2, vze2
use misc_coms,   only: dtsm, initial, dn01d, th01d, deltax, nxp, &
                       mdomain, icorflg
use consts_coms, only: cpocv, omega2, pi1, pio180, r8
use mem_grid,    only: mza, mva, mwa, lpv, lpw, lve2, arv, arw, &
                       vnx, vny, vnz, wnx, wny, wnz, wnxo2, wnyo2, wnzo2, &
                       dzim, volt, volti, glatw, glonw, &
                       dzt_top, dzt_bot, zwgt_top, zwgt_bot, gravm
use tridiag,     only: tridiffo
use oname_coms,  only: nl
use mem_turb,    only: akmodx, akhodx
use mem_rayf,    only: dorayfw, rayf_cofw, krayfw_bot

implicit none

integer, intent(in) :: iw

real, intent(in) :: vmcf(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)
real, intent(in) :: alpha_press(mza,mwa)
real, intent(in) :: rhot(mza,mwa)

real, intent(in) :: thil_upv(mza,mva)
real, intent(in) :: vxe_upv (mza,mva)
real, intent(in) :: vye_upv (mza,mva)
real, intent(in) :: vze_upv (mza,mva)

real, intent(in) :: thil_upw(mza,mwa)
real, intent(in) :: vxe_upw (mza,mwa)
real, intent(in) :: vye_upw (mza,mwa)
real, intent(in) :: vze_upw (mza,mwa)

real, intent(inout) :: vmxet_volt(mza,mwa)
real, intent(inout) :: vmyet_volt(mza,mwa)
real, intent(inout) :: vmzet_volt(mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

real, intent(inout) :: vmx_cor(mza,mwa)
real, intent(inout) :: vmy_cor(mza,mwa)

integer :: jv, iv, iwn
integer :: k, ka, kbv, kp, ksw
integer :: npoly

real :: dts, omeg2
real :: c6, c7, c8, c9, c10
real :: dirv, vmarv
real :: del_rhothil, vmt1
real :: rad0_swtc, rad_swtc, topo_swtc

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55  ! wmc
real, parameter :: fr = .75  ! rho
real, parameter :: fp = .75  ! press

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

real :: wmt_rayf     (mza)

real :: hflux_thil(mza)
real :: hflux_vxe(mza)
real :: hflux_vye(mza)
real :: hflux_vze(mza)

real :: vflux_thil(mza)
real :: vflux_vxe(mza)
real :: vflux_vye(mza)
real :: vflux_vze(mza)

real :: b1(mza),b2(mza),b3(mza),b5(mza),b6(mza),b10(mza)
real :: b7(mza),b8(mza),b9(mza),b11(mza),b12(mza),b13(mza),b14(mza)
real :: b21(mza),b22(mza),b23(mza),b24(mza),b25(mza),b26(mza)
real :: b31(mza),b32(mza),b33(mza),b34(mza)

real :: vmxe1(mza)
real :: vmye1(mza)
real :: vmze1(mza)

real :: vxe1(mza)
real :: vye1(mza)
real :: vze1(mza)

ka = lpw(iw)

dts = dtsm(itab_w(iw)%mrlw)

! Set bottom & top vertical advective mass and heat fluxes to zero

wmarw(1:ka-1) = 0.
wmarw(mza) = 0.

vflux_thil(ka-1)  = 0.
vflux_thil(mza) = 0.

vflux_vxe(ka-1)  = 0.
vflux_vxe(mza) = 0.

vflux_vye(ka-1)  = 0.
vflux_vye(mza) = 0.

vflux_vze(ka-1)  = 0.
vflux_vze(mza) = 0.

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

! Number of edges of this IW polygon

npoly = itab_w(iw)%npoly

! Loop over V neighbors of this W cell

do jv = 1,npoly

   iv   = itab_w(iw)%iv(jv)
   kbv  = lpv(iv)
   dirv = itab_w(iw)%dirv(jv)
   iwn  = itab_w(iw)%iw(jv)

! Loop over T levels

   do k = kbv,mza

      vmarv = dirv * vmcf(k,iv) * arv(k,iv)

! Sum horizontal advection fluxes over V faces      

      hflux_rho(k)  = hflux_rho(k)  + vmarv

      hflux_thil(k) = hflux_thil(k) + vmarv * thil_upv(k,iv) &
                                    + akhodx(k,iv) * (thil(k,iwn) - thil(k,iw))

      hflux_vxe(k)  = hflux_vxe(k)  + vmarv * vxe_upv(k,iv)  &
                                    + akmodx(k,iv) * (vxe(k,iwn) - vxe(k,iw))

      hflux_vye(k)  = hflux_vye(k)  + vmarv * vye_upv(k,iv)  &
                                    + akmodx(k,iv) * (vye(k,iwn) - vye(k,iw))

      hflux_vze(k)  = hflux_vze(k)  + vmarv * vze_upv(k,iv)  &
                                    + akmodx(k,iv) * (vze(k,iwn) - vze(k,iw))

   enddo

enddo

! Coriolis force

if (icorflg == 1) then
   omeg2 = omega2
else
   omeg2 = 0.
endif

! Loop over T levels

do k = ka,mza

   ! Coriolis force for T cell

   vmx_cor(k,iw) =  omeg2 * vye(k,iw) * rho(k,iw)
   vmy_cor(k,iw) = -omeg2 * vxe(k,iw) * rho(k,iw)

   ! Explicit density tendency

   delex_rho(k) = dts * (rhot(k,iw) &
      + volti(k,iw) * (hflux_rho(k) + wmarw(k-1) - wmarw(k)))

   ! Explicit density-thil tendency

   if (nl%split_scalars == 1) then
      delex_rhothil(k) = dts * (thil(k,iw) * rhot(k,iw) &
           + volti(k,iw) * (hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)))
   else
      delex_rhothil(k) = dts * (thilt(k,iw) + thil(k,iw) * rhot(k,iw) &
           + volti(k,iw) * (hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)))
   endif

   ! Explicit momentum tendency from advective and turbulent transport
   ! (weighted by T cell volume)

   vmxet_volt(k,iw) = hflux_vxe(k) + vflux_vxe(k-1) - vflux_vxe(k)
   vmyet_volt(k,iw) = hflux_vye(k) + vflux_vye(k-1) - vflux_vye(k)
   vmzet_volt(k,iw) = hflux_vze(k) + vflux_vze(k-1) - vflux_vze(k)

   ! RHOTHIL(t) and PRESS(t)

   rhothil(k) = rho(k,iw) * thil(k,iw)
   press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

   ! Compute current T cell momentum and store in temp array

   vmxe1(k) = vxe(k,iw) * rho(k,iw)
   vmye1(k) = vye(k,iw) * rho(k,iw)
   vmze1(k) = vze(k,iw) * rho(k,iw)
enddo

! Rayleigh friction on W

wmt_rayf(:) = 0.

if (dorayfw) then
   do k = krayfw_bot, mza-1
      wmt_rayf(k) = - rayf_cofw(k) * wmc(k,iw)
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

   delex_wm(k) = dts * (wmt_rayf(k) &

      + dzim(k) * (press_t(k) - press_t(k+1) &

      - gravm(k) * (dzt_top(k) * rho(k,iw) + dzt_bot(k+1) * rho(k+1,iw))) &

! Average long-timestep tendencies and Coriolis force in vertical between T levels

      + wnxo2(iw) * (vmxet(k,iw) + vmxet(k+1,iw) + vmx_cor(k,iw) + vmx_cor(k+1,iw)) &
      + wnyo2(iw) * (vmyet(k,iw) + vmyet(k+1,iw) + vmy_cor(k,iw) + vmy_cor(k+1,iw)) &
      + wnzo2(iw) * (vmzet(k,iw) + vmzet(k+1,iw)) &

! ADD volume-weighted VMXET, VMYET, VMZET over both T levels and divide by
! combined VOLT of those T levels

      + (wnx(iw) * (vmxet_volt(k,iw) + vmxet_volt(k+1,iw)) &
      +  wny(iw) * (vmyet_volt(k,iw) + vmyet_volt(k+1,iw)) &
      +  wnz(iw) * (vmzet_volt(k,iw) + vmzet_volt(k+1,iw))) &
      / (volt(k,iw) + volt(k+1,iw)) )

enddo

c6  = dts * .5 * fw
c7  = dts * .25 * fw
c8  = dts * pc2
c9  =-dts * fr
c10 = dts * fw

! Fill matrix coefficients for implicit update of WM

do k = ka,mza
   kp = min(k+1,mza)
   b1(k)  = wc(k,iw) + wc(k-1,iw)           ! T pts
   b2(k)  = thil(k,iw) + thil(kp,iw)        ! W pts
   b3(k)  = 2. / (volt(k,iw) + volt(kp,iw)) ! W pts [b3 replaces volwi]
   b5(k)  = press_t(k) / rhothil(k)         ! T pts
   b6(k)  = c6 * volti(k,iw)                ! T pts
   b10(k) = c10 * volti(k,iw)               ! T pts
enddo
b2(ka-1) = b2(ka)
b3(ka)   = 1. / (volt(ka,iw) + .5 * volt(ka+1,iw))

do k = ka,mza-1

   b7(k)  = c7     * b3(k)        ! W pts
   b8(k)  = c8     * dzim(k)      ! W pts
   b9(k)  = c9     * dzim(k) * gravm(k) ! W pts
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

   b32(k) = 1. + arw(k,iw) &
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

del_wmarw(ka-1)  = 0.
del_wmarw(mza) = 0.

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

! Estimate velocity in T cells at (t+1) by prognostic method

   vxe1(k) = (vmxe1(k) + dts * (vmxet(k,iw) + vmxet_volt(k,iw) * volti(k,iw) + vmx_cor(k,iw))) / rho(k,iw)
   vye1(k) = (vmye1(k) + dts * (vmyet(k,iw) + vmyet_volt(k,iw) * volti(k,iw) + vmy_cor(k,iw))) / rho(k,iw)
   vze1(k) = (vmze1(k) + dts * (vmzet(k,iw) + vmzet_volt(k,iw) * volti(k,iw)                )) / rho(k,iw)

! Add half-forward T cell velocity to scalar arrays

   vxesc(k,iw) = vxesc(k,iw) + .5 * (vxe(k,iw) + vxe1(k))
   vyesc(k,iw) = vyesc(k,iw) + .5 * (vye(k,iw) + vye1(k))
   vzesc(k,iw) = vzesc(k,iw) + .5 * (vze(k,iw) + vze1(k))

enddo

if (lve2(iw) > 0) then

! Zero out vxe2, vye2, vze2 prior to new diagnosis

   do ksw = 1, lve2(iw)
      vxe2(ksw,iw) = 0.
      vye2(ksw,iw) = 0.
      vze2(ksw,iw) = 0.
   enddo

! Loop over adjacent V faces

   do jv = 1, npoly
      iv = itab_w(iw)%iv(jv)

! Project full-forward-time vxe1, vye1, vze1 onto V faces that are below ground,
! and then project back to vxe2, vye2, vze2

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
   wc (k,iw) = wmc(k,iw) / (zwgt_bot(k) * rho(k,iw) + zwgt_top(k) * rho(k+1,iw))

enddo

! Set top & bottom values of WC

!wc(ka-1,iw) = wnx(iw) * vxe1(ka) + wny(iw) * vye1(ka) + wnz(iw) * vze1(ka)
wc(ka-1,iw) = wc(ka,iw)

wc(1:ka-2,iw) = 0.
wc(mza   ,iw) = 0.

end subroutine prog_wrt_begs

!============================================================================

subroutine prog_v_begs(iv,vmxet_volt,vmyet_volt,vmzet_volt,div2d,vmx_cor,vmy_cor, &
                       vortp,vortn,vortp_t,rotational)

use mem_ijtabs,  only: itab_v
use mem_basic,   only: vc, press, vmp, vmc, rho, vxe, vye, vze
use misc_coms,   only: dtsm, initial, mdomain, u01d, v01d, dn01d, &
                       deltax, nxp
use consts_coms, only: eradi, gravo2
use mem_grid,    only: mza, mma, mva, mwa, lpv, volt, xev, yev, zev, &
                       unx, uny, unz, vnx, vny, vnz, vnxo2, vnyo2, vnzo2, &
                       dniu, dniv, arw0, dnu, c1, c2 
use mem_tend,    only: vmxet, vmyet, vmzet
use mem_rayf,    only: dorayf, rayf_cof, vc03d, dn03d, krayf_bot, &
                       dorayfdiv, krayfdiv_bot, rayf_cofdiv
use oname_coms,  only: nl

implicit none

integer, intent(in) :: iv

real, intent(in) :: vmxet_volt(mza,mwa)
real, intent(in) :: vmyet_volt(mza,mwa)
real, intent(in) :: vmzet_volt(mza,mwa)
real, intent(in) :: div2d(mza,mwa)
real, intent(in) :: vmx_cor(mza,mwa)
real, intent(in) :: vmy_cor(mza,mwa)
real, intent(in) :: vortp  (mza,mma)
real, intent(in) :: vortn  (mza,mva)
real, intent(in) :: vortp_t(mza,mwa)

logical, intent(in) :: rotational

integer :: k, kb
integer :: iw1,iw2,im1,im2,im3,im4,im5,im6
real :: sum1,sum2,vmp_eqdiv
real :: dts
real :: fracx, rayfx
real :: vx, vy, vz, uc, watv, tke1, tke2, vortp_v, dtso2dnu, dtso2dnv
real :: vort_big1, vort_big2

real :: vmt_rayf(mza)
real :: vmt_vorf(mza)
real :: pgf     (mza)

real, parameter :: onethird = 1./3.

! Extract neighbor indices and coefficients for this point in the U stencil

iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2); im3 = itab_v(iv)%im(3)
im4 = itab_v(iv)%im(4); im5 = itab_v(iv)%im(5); im6 = itab_v(iv)%im(6)

dts = dtsm(itab_v(iv)%mrlv)

kb = lpv(iv)

vmt_rayf(:) = 0.

! Horizontal divergence damping

if (dorayfdiv) then

   ! Vertical loop over V points

   do k = krayfdiv_bot, mza

      ! Divergence in IW1 excluding IV

      sum1 = div2d(k,iw1) - vmp(k,iv) * dnu(iv)

      ! Divergence in IW2 excluding IV

      sum2 = div2d(k,iw2) + vmp(k,iv) * dnu(iv)

      ! VMP value that would equalize horizontal divergence in IW1 and IW2
      ! cells (assuming no blockage by topography)

      vmp_eqdiv = (arw0(iw1) * sum2 - arw0(iw2) * sum1) &
           / (dnu(iv) * (arw0(iw1) + arw0(iw2)))

      vmt_rayf(k) = vmt_rayf(k) + rayf_cofdiv(k) * (vmp_eqdiv - vmp(k,iv))

   enddo

endif ! (dorayfdiv)

! Rayleigh friction on vc: Only for horizontally homogeneous initialization

if (dorayf .and. initial == 1) then

   ! Vertical loop over V points
      
   do k = krayf_bot, mza
      vmt_rayf(k) = vmt_rayf(k) + rayf_cof(k) * dn03d(k,iv) * (vc03d(k,iv) - vc(k,iv))
   enddo

   ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions

   ! fracx = abs(xev(iv)) / (real(nxp-1) * .866 * deltax)
   ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
   ! rayfx = 0.   ! Default: no extra RAYF
   ! do k = kb, mza
   !    vmt_rayf(k) = vmt_rayf(k) &
   !                + max(rayf_cof(k),rayfx) * dn03d(k) *  (vc03d(k,iv) - vc(k,iv))
   ! enddo

endif ! (dorayf and initial == 1)
   
! Vertical loop over V points

do k = kb,mza

   ! Horizontal filter for vertical vorticity

   vort_big1 = onethird * (vortp(k,im2) + vortp(k,im3) + vortp(k,im4))
   vort_big2 = onethird * (vortp(k,im1) + vortp(k,im5) + vortp(k,im6))
   
   vmt_vorf(k) = (rho(k,iw1) + rho(k,iw2)) &
               * ( (vort_big1 - vortp(k,im1)) * c1(iv) &
                 + (vort_big2 - vortp(k,im2)) * c2(iv) )

   ! Pressure gradient force

   pgf(k) = dniv(iv) * (press(k,iw1) - press(k,iw2))
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
      vmc(k,iv) = vmc(k,iv) + dts * (vmt_rayf(k) + vmt_vorf(k) + pgf(k) &
                + vnxo2(iv) * (vmxet(k,iw1) + vmxet(k,iw2) + vmx_cor(k,iw1) + vmx_cor(k,iw2)) &
                + vnyo2(iv) * (vmyet(k,iw1) + vmyet(k,iw2) + vmy_cor(k,iw1) + vmy_cor(k,iw2)) &
                + vnzo2(iv) * (vmzet(k,iw1) + vmzet(k,iw2))  &
                +  (vnx(iv) * (vmxet_volt(k,iw1) + vmxet_volt(k,iw2))  &
                +   vny(iv) * (vmyet_volt(k,iw1) + vmyet_volt(k,iw2))  &
                +   vnz(iv) * (vmzet_volt(k,iw1) + vmzet_volt(k,iw2))) &
                / (volt(k,iw1) + volt(k,iw2)) )

      vc(k,iv) = 2.0 * vmc(k,iv) / (rho(k,iw1) + rho(k,iw2))
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

      vmc(k,iv) = vmc(k,iv) + dts * (vmt_rayf(k) + vmt_vorf(k) + pgf(k) &
                + vnxo2(iv) * (vmxet(k,iw1) + vmxet(k,iw2) + vmx_cor(k,iw1) + vmx_cor(k,iw2)) &
                + vnyo2(iv) * (vmyet(k,iw1) + vmyet(k,iw2) + vmy_cor(k,iw1) + vmy_cor(k,iw2)) &
                + vnzo2(iv) * (vmzet(k,iw1) + vmzet(k,iw2)) &
                + .5 * (rho(k,iw1) + rho(k,iw2)) &
                * (dniv(iv) * (tke1 - tke2) &
                + uc * vortp_v &
                - watv * .5 * (vortn(k-1,iv) + vortn(k,iv))))

      vc(k,iv) = 2.0 * vmc(k,iv) / (rho(k,iw1) + rho(k,iw2))
   enddo

endif

end subroutine prog_v_begs
