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
                        mrl_begl, mrl_begs, mrl_ends, mrl_endl, &
                        jtm_vadj, jtv_prog, jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp
use mem_basic,    only: rho, thil, theta, wc, press, wmc, vmp, vmc, vp, vc, &
                        vxe2, vye2, vze2, &
                        sh_w, sh_v, vxe, vye, vze, strict_wvt_donorpoint
use mem_grid,     only: mza, mma, mva, mwa, nsw_max, lpm, lpv, lpw, lsw, &
                        zt, zm, dzim, zfacit, dzm, dzt, dnv, dnu, arm0, arv, &
                        vnx, vny, vnz, volt, glatw, glonw
use mem_tend,     only: vmt, vmxet, vmyet, vmzet, sh_wt
use misc_coms,    only: io6, iparallel, time8, dtlm, rinit
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
use oplot_coms,   only: op
use obnd,         only: lbcopy_m, lbcopy_w
use mem_rayf,     only: dorayfdiv, krayfdiv_bot

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

real, intent(in)    :: alpha_press(mza,mwa)
real, intent(inout) :: rhot       (mza,mwa)

integer :: j, iv, k, ka, kb, mrl, kbv, kd, ksw
integer :: iw, iw1, iw2, iwp, ivp

! automatic arrays

integer :: iwdepv(mza,mva) ! donor cell IW index for V face
integer :: iwrecv(mza,mva) ! recvr cell IW index for V face
integer :: kdepw(mza,mwa)  ! donor cell K index for W face
integer :: krecw(mza,mwa)  ! recvr cell K index for W face

real :: dxps_v(mza,mva) ! X component in PS projection of displacement for V face
real :: dyps_v(mza,mva) ! Y component in PS projection of displacement for V face
real :: dzps_v(mza,mva) ! Z component in PS projection of displacement for V face

real :: dxps_w(mza,mwa) ! X component in PS projection of displacement for W face
real :: dyps_w(mza,mwa) ! Y component in PS projection of displacement for W face
real :: dzps_w(mza,mwa) ! Z component in PS projection of displacement for W face

real :: vmcf(mza,mva) ! Time-extrapolated VMC

real, allocatable :: vcf (:,:) ! Time-extrapolated VC
real, allocatable :: vxef(:,:) ! Time-extrapolated XE velocity component at T point
real, allocatable :: vyef(:,:) ! Time-extrapolated YE velocity component at T point
real, allocatable :: vzef(:,:) ! Time-extrapolated ZE velocity component at T point

real :: gxps_thil(mza,mwa) ! X component in PS projection of THIL gradient 
real :: gyps_thil(mza,mwa) ! Y component in PS projection of THIL gradient 
real :: gzps_thil(mza,mwa) ! Z component in PS projection of THIL gradient 

real :: gxps_vxe(mza,mwa) ! X component in PS projection of VXE gradient
real :: gyps_vxe(mza,mwa) ! Y component in PS projection of VXE gradient
real :: gzps_vxe(mza,mwa) ! Z component in PS projection of VXE gradient

real :: gxps_vye(mza,mwa) ! X component in PS projection of VYE gradient
real :: gyps_vye(mza,mwa) ! Y component in PS projection of VYE gradient
real :: gzps_vye(mza,mwa) ! Z component in PS projection of VYE gradient

real :: gxps_vze(mza,mwa) ! X component in PS projection of VZE gradient
real :: gyps_vze(mza,mwa) ! Y component in PS projection of VZE gradient
real :: gzps_vze(mza,mwa) ! Z component in PS projection of VZE gradient

real :: thil_s(mza,mwa)

real :: vortp(mza,mma)
real :: div2d(mza,mwa)

real :: thil_upv(mza,mva) ! Upstreamed THIL at each V interface
real :: vxe_upv (mza,mva) ! Upstreamed VXE  at each V interface
real :: vye_upv (mza,mva) ! Upstreamed VYE  at each V interface
real :: vze_upv (mza,mva) ! Upstreamed VZE  at each V interface

real :: thil_upw(mza,mwa) ! Upstreamed THIL at each W level
real :: vxe_upw (mza,mwa) ! Upstreamed VXE  at each W level
real :: vye_upw (mza,mwa) ! Upstreamed VYE  at each W level
real :: vze_upw (mza,mwa) ! Upstreamed VZE  at each W level

integer :: iv1,iv2,iv3,iv4,im,npoly,jv,im1,im2,im3,im4,im5,im6,iwd
real :: c1,c2,vort_big1,vort_big2
real :: arm0i

! Parameters for vorticity diffusion

real, parameter :: tvort = 3600.
real, parameter :: c0 = .125 / tvort

vortp = rinit

! Half-forward velocities for computing donor points:

if (strict_wvt_donorpoint) then
   allocate(vcf (mza,mva)) ; vcf  = 0.0
   allocate(vxef(mza,mwa)) ; vxef = 0.0
   allocate(vyef(mza,mwa)) ; vyef = 0.0
   allocate(vzef(mza,mwa)) ; vzef = 0.0
endif

! Save copy of thil

thil_s(:,:) = thil(:,:)

vmcf(:,1) = 0.

! First compute long timestep tendencies
! Maybe move to a separate subroutine veltend_long_hex?

mrl = mrl_begl(istp)
if (mrl > 0) then

! Horizontal loop over W columns for BEGL

!----------------------------------------------------------------------
   !$omp parallel do private(iw,k) 
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Include moisture changes in total density tendency

      do k = lpw(iw), mza
         rhot(k,iw) = rhot(k,iw) + sh_wt(k,iw)
      enddo

      call prog_wrt_begl(iw,rhot)
   
   enddo
   !$omp end parallel do

! MPI SEND of VMXET, VMYET, VMZET

   if (iparallel == 1) then
      call mpi_send_w(mrl, rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)
   endif

! Horizontal loop over M/P columns for BEGL; Diagnose vertical vorticity
! in preparation for horizontal filter

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

! Parallel send of vortp

   if (iparallel == 1) call mpi_send_m(mrl, rvara1=vortp)

! MPI RECV of VMXET, VMYET, VMZET

   if (iparallel == 1) then
      call mpi_recv_w(mrl, rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)
   endif

   call lbcopy_w(mrl, a1=vmxet, a2=vmyet, a3=vmzet)

endif ! mrl = mrl_begl(istp) > 0

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

      call donorpointw(0, mrl, wc, vxef, vyef, vzef, kdepw, krecw, &
                       dxps_w, dyps_w, dzps_w)

      ! Finish MPI recv of VXE, VYE, VZE and do a LBC copy

      if (iparallel == 1) then
         call mpi_recv_w(mrl, rvara1=vxef, rvara2=vyef, rvara3=vzef)
      endif

      call lbcopy_w(mrl, a1=vxef, a2=vyef, a3=vzef)

      ! Diagnose advective donor point locations for the V faces surrounding all
      ! primary W points. Communication of velocities must have been completed

      call donorpointv(0, mrl, vcf, vxef, vyef, vzef, iwdepv, iwrecv, &
                       dxps_v, dyps_v, dzps_v)
   endif

else

! Compute donor point locations using current velocities.
   
   mrl = mrl_begs(istp)

   if (mrl > 0) then
      call donorpointw(0, mrl, wc, vxe, vye, vze, kdepw, krecw, &
                       dxps_w, dyps_w, dzps_w)

      call donorpointv(0, mrl, vc, vxe, vye, vze, iwdepv, iwrecv, &
                       dxps_v, dyps_v, dzps_v)
   endif

endif  ! strict_wvt_donorpoint

mrl = mrl_begl(istp)
if (mrl > 0) then

! MPI RECV and LBC of VORTP

   if (iparallel == 1) call mpi_recv_m(mrl, rvara1=vortp)
   call lbcopy_m(mrl, a1=vortp)

! Horizontal loop over V columns for PROG_V_BEGL

!----------------------------------------------------------------------
   !$omp parallel do private(iv,iw1,iw2,im1,im2,im3,im4,im5,im6,iv1,iv2,iv3,iv4, &
   !$omp                     kb,k,c1,c2,vort_big1,vort_big2) 
   do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

      kb = lpv(iv)

! Vertical loop over T levels

      do k = kb,mza

! Update VM tendency from turbulent fluxes

         vmt(k,iv) = vmt(k,iv) + .5 &
                   * (vnx(iv) * (vmxet(k,iw1) + vmxet(k,iw2)) &
                   +  vny(iv) * (vmyet(k,iw1) + vmyet(k,iw2)) &
                   +  vnz(iv) * (vmzet(k,iw1) + vmzet(k,iw2)))

      enddo

!------------------------------------------------------
! SPECIAL - HORIZONTAL FILTER FOR VERTICAL VORTICITY

      im1  = itab_v(iv)%im(1)
      im2  = itab_v(iv)%im(2)
      im3  = itab_v(iv)%im(3)
      im4  = itab_v(iv)%im(4)
      im5  = itab_v(iv)%im(5)
      im6  = itab_v(iv)%im(6)

      iv1  = itab_v(iv)%iv(1)
      iv2  = itab_v(iv)%iv(2)
      iv3  = itab_v(iv)%iv(3)
      iv4  = itab_v(iv)%iv(4)
   
      c1 = -c0 * arm0(im1) * dnu(iv1) * dnu(iv2) / (dnu(iv1) * dnu(iv2) * dnv(iv) &
         + dnu(iv) * dnu(iv2) * dnv(iv1) + dnu(iv) * dnu(iv1) * dnv(iv2))

      c2 = c0 * arm0(im2) * dnu(iv3) * dnu(iv4) / (dnu(iv3) * dnu(iv4) * dnv(iv) &
         + dnu(iv) * dnu(iv4) * dnv(iv3) + dnu(iv) * dnu(iv3) * dnv(iv4))

! Vertical loop over V levels (for now, don't check for k >= lpm, etc.)

      do k = kb,mza
   
         vort_big1 = .3333333 * (vortp(k,im2) + vortp(k,im3) + vortp(k,im4))
         vort_big2 = .3333333 * (vortp(k,im1) + vortp(k,im5) + vortp(k,im6))
   
         vmt(k,iv) = vmt(k,iv) + (rho(k,iw1) + rho(k,iw2)) &
            * ((vort_big1 - vortp(k,im1)) * c1 + (vort_big2 - vortp(k,im2)) * c2)

      enddo

! END - HORIZONTAL FILTER FOR VERTICAL VORTICITY
!------------------------------------------------------
              
   enddo
   !$omp end parallel do 

endif ! mrl = mrl_begl(istp) > 0

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

! [Now, we can reuse vmxet, vmyet, vmzet arrays]

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

! Evaluate T3D gradients of THIL, VXE, VYE, and VZE for BEGS

   call grad_t3d(mrl,thil,gxps_thil,gyps_thil,gzps_thil)
   call grad_t3d(mrl,vxe,gxps_vxe,gyps_vxe,gzps_vxe)
   call grad_t3d(mrl,vye,gxps_vye,gyps_vye,gzps_vye)
   call grad_t3d(mrl,vze,gxps_vze,gyps_vze,gzps_vze)

! Finish MPI RECV of DIV2D, and MPI SEND of THIL, VXE, VYE, and VZE 
! gradient components (12 in all)

   if (iparallel == 1) then

      if (dorayfdiv) then
         call mpi_recv_w(mrl, rvara1=div2d)
      endif

      call mpi_send_w(mrl, &
         rvara1 =gxps_thil, rvara2 =gyps_thil, rvara3 =gzps_thil, &
         rvara4 =gxps_vxe,  rvara5 =gyps_vxe,  rvara6 =gzps_vxe,  &
         rvara7 =gxps_vye,  rvara8 =gyps_vye,  rvara9 =gzps_vye,  &
         rvara10=gxps_vze,  rvara11=gyps_vze,  rvara12=gzps_vze   )

   endif

!  Horizontal loop over all primary W columns

   !$omp parallel do private(iw,k,kd) 
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
      
      do k = lpw(iw), mza

         kd = kdepw(k,iw)

         thil_upw(k,iw) = thil_s(kd,iw)                   &
                        + dxps_w(k,iw) * gxps_thil(kd,iw) &
                        + dyps_w(k,iw) * gyps_thil(kd,iw) &
                        + dzps_w(k,iw) * gzps_thil(kd,iw)

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

      enddo

   enddo
   !$omp end parallel do

   if (iparallel == 1) then

      call mpi_recv_w(mrl, &
         rvara1 =gxps_thil, rvara2 =gyps_thil, rvara3 =gzps_thil, &
         rvara4 =gxps_vxe,  rvara5 =gyps_vxe,  rvara6 =gzps_vxe,  &
         rvara7 =gxps_vye,  rvara8 =gyps_vye,  rvara9 =gzps_vye,  &
         rvara10=gxps_vze,  rvara11=gyps_vze,  rvara12=gzps_vze   )

   endif

   call lbcopy_w(mrl, a1 =gxps_thil, a2 =gyps_thil, a3 =gzps_thil, &
                      a4 =gxps_vxe,  a5 =gyps_vxe,  a6 =gzps_vxe,  &
                      a7 =gxps_vye,  a8 =gyps_vye,  a9 =gzps_vye,  &
                      a10=gxps_vze,  a11=gyps_vze,  a12=gzps_vze   )

!  Horizontal loop over V

!----------------------------------------------------------------------
   !$omp parallel do private(iv,iwd,kbv,k) 
   do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)
!----------------------------------------------------------------------

      kbv = lpv(iv)
      do  k = kbv, mza

         iwd = iwdepv(k,iv)

         thil_upv(k,iv) = thil(k,iwd)                     &
                        + dxps_v(k,iv) * gxps_thil(k,iwd) &
                        + dyps_v(k,iv) * gyps_thil(k,iwd) &
                        + dzps_v(k,iv) * gzps_thil(k,iwd)

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
         
      enddo

   enddo
   !$omp end parallel do

endif

! Main loop over W columns for updating WM, WC, RHO, THIL, and PRESS

!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw,ksw) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Prognose vertical velocity, density, thil, and diagnose pressure

   call prog_wrt_begs( iw, vmcf, wmsc, alpha_press, rhot, thil_s, &
                       thil_upv, vxe_upv, vye_upv, vze_upv,       &
                       thil_upw, vxe_upw, vye_upw, vze_upw,       &
                       vxesc, vyesc, vzesc                        )

enddo
!$omp end parallel do 
endif

! MPI SEND/RECV and LBC copy of quantities needed for prog_v: 
! PRESS, RHO, VMXET, VMYET, and VMZET

if (iparallel == 1) then
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=vmxet, rvara2=vmyet, rvara3=vmzet)
endif

call lbcopy_w(mrl, a1=vmxet, a2=vmyet, a3=vmzet, a4=wmc, &
                   a5=thil,  a6=wc,    d1=press, d2=rho)

! Horizontal loop over V points to update VMC

!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iv) 
do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!----------------------------------------------------------------------

   call prog_v_begs(iv,vmxet,vmyet,vmzet,div2d)

enddo
!$omp end parallel do 
endif

if (strict_wvt_donorpoint) then
   deallocate( vcf)
   deallocate(vxef)
   deallocate(vyef)
   deallocate(vzef)
endif

return
end subroutine prog_wrtv

!=========================================================================

subroutine prog_wrt_begl(iw,rhot)

! This version includes turbulent fluxes of VXE, VYE, VZE through V faces
! Vertical mixing is computed in subroutine pbl_driver

! All diffusive tendencies are evaluated at T points

use mem_tend,    only: thilt, wmt, vmxet, vmyet, vmzet
use mem_ijtabs,  only: istp, itab_w, mrl_begl, mrl_begr, mrl_begs, mrl_endr
use mem_basic,   only: wmc, rho, thil, wc, vc, theta, vxe, vye, vze
use misc_coms,   only: io6, initial, dn01d, th01d, &
                       deltax, nxp, mdomain, time8, dtlm
use mem_grid,    only: mza, mva, mwa, lpv, lpw, arv, dniv, volt, volti, &
                       xew, vnx, vny, vnz, wnxo2, wnyo2, wnzo2
use mem_turb,    only: hkm
use mem_rayf,    only: rayf_cof, rayf_cofw, dorayf, dorayfw, krayf_bot, krayfw_bot

implicit none

integer, intent(in) :: iw
real,    intent(in) :: rhot(mza,mwa)

integer :: iv, iwn, k, ka, kbv, npoly, jv, ksw
real    :: fracx, rayfx
real    :: arvkodx, hdniv, vmt1

! Automatic arrays:

real :: hdiff_vxe(mza)
real :: hdiff_vye(mza)
real :: hdiff_vze(mza)

ka = lpw(iw)

! Number of edges of this IW polygon

npoly = itab_w(iw)%npoly

hdiff_vxe(:) = 0.
hdiff_vye(:) = 0.
hdiff_vze(:) = 0.

! Loop over V neighbors of this W cell

do jv = 1,npoly
   iv  = itab_w(iw)%iv(jv)
   iwn = itab_w(iw)%iw(jv)
   kbv  = lpv(iv)
   hdniv = .5 * dniv(iv)

! Vertical loop over T levels

   do k = kbv, mza

! Horizontal diffusive flux coefficient

      arvkodx = hdniv * arv(k,iv) * (hkm(k,iwn) + hkm(k,iw))

! Compute and sum horizontal diffusive flux across this V neighbor

      hdiff_vxe(k) = hdiff_vxe(k) + arvkodx * (vxe(k,iwn) - vxe(k,iw))
      hdiff_vye(k) = hdiff_vye(k) + arvkodx * (vye(k,iwn) - vye(k,iw))
      hdiff_vze(k) = hdiff_vze(k) + arvkodx * (vze(k,iwn) - vze(k,iw))

   enddo

enddo

! Vertical loop over T levels

do k = ka,mza

! Evaluate momentum tendency in T cell from horizontal turbulent transport

   vmxet(k,iw) = vmxet(k,iw) + volti(k,iw) * hdiff_vxe(k) &
               + vxe(k,iw) * rhot(k,iw)
   vmyet(k,iw) = vmyet(k,iw) + volti(k,iw) * hdiff_vye(k) &
               + vye(k,iw) * rhot(k,iw)
   vmzet(k,iw) = vmzet(k,iw) + volti(k,iw) * hdiff_vze(k) &
               + vze(k,iw) * rhot(k,iw)

enddo

! Vertical loop over W levels

do k = ka,mza-1

! Update WM tendency from turbulent fluxes

   wmt(k,iw) = wmt(k,iw)                                   &
             + ( wnxo2(iw) * (vmxet(k,iw) + vmxet(k+1,iw)) &
               + wnyo2(iw) * (vmyet(k,iw) + vmyet(k+1,iw)) &
               + wnzo2(iw) * (vmzet(k,iw) + vmzet(k+1,iw)) )
enddo

! RAYLEIGH FRICTION ON WM

if (dorayfw) then

!! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN
!! fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
!! rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-1)
!! rayfx = 0.   ! Default: no extra RAYF
!! do k = ka, mza-1
!!    wmt(k,iw) = wmt(k,iw) - max(rayf_cofw(k),rayfx) * wmc(k,iw)
!! enddo
!! END SPECIAL

   do k = krayfw_bot, mza-1
      wmt(k,iw) = wmt(k,iw) - rayf_cofw(k) * wmc(k,iw)
   enddo

endif ! (dorayfw)

! RAYLEIGH FRICTION ON THIL

if (dorayf) then

   if (initial == 1) then   ! HHI case

!! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN
!!    fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
!!    rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
!!    rayfx = 0.   ! Default: no extra RAYF
!!    do k = ka, mza
!!       thilt(k,iw) = thilt(k,iw) + max(rayf_cof(k),rayfx)  & 
!!                   * dn01d(k) * (th01d(k) - theta(k,iw))
!!    enddo
!! END SPECIAL

      do k = krayf_bot, mza
         thilt(k,iw) = thilt(k,iw) + &
                       rayf_cof(k) * dn01d(k) * (th01d(k) - theta(k,iw))
      enddo

   else                     ! LHI/VARI case
! Need implementation for LHI/VARI (use vartp for merid. variation?)
   endif

endif ! (dorayf)

return
end subroutine prog_wrt_begl

!=========================================================================

subroutine prog_wrt_begs( iw, vmcf, wmsc, alpha_press, rhot, thil_s, &
                          thil_upv, vxe_upv, vye_upv, vze_upv,       &
                          thil_upw, vxe_upw, vye_upw, vze_upw,       &  
                          vxesc, vyesc, vzesc                        )

use mem_tend,    only: thilt, wmt, vmxet, vmyet, vmzet
use mem_ijtabs,  only: itab_w
use mem_basic,   only: wmc, rho, thil, wc, vc, theta, press, &
                       vxe, vye, vze, vxe2, vye2, vze2
use misc_coms,   only: io6, dtsm, initial, dn01d, th01d, deltax, nxp, &
                       mdomain, time8, icorflg
use consts_coms, only: cpocv, grav, omega2, pi1, pio180, r8
use mem_grid,    only: mza, mva, mwa, nsw_max, lpv, lpw, lve2, arv, arw, &
                       vnx, vny, vnz, wnx, wny, wnz, wnxo2, wnyo2, wnzo2, &
                       xew, dzim, dzt, volt, volti, glatw, glonw, &
                       dzt_top, dzt_bot, zwgt_top, zwgt_bot
use tridiag,     only: tridiffo
use oname_coms,  only: nl

implicit none

integer, intent(in) :: iw

real, intent(in) :: vmcf(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)
real, intent(in) :: alpha_press(mza,mwa)
real, intent(in) :: rhot(mza,mwa)

real, intent(in) :: thil_s(mza,mwa)

real, intent(in) :: thil_upv(mza,mva)
real, intent(in) :: vxe_upv (mza,mva)
real, intent(in) :: vye_upv (mza,mva)
real, intent(in) :: vze_upv (mza,mva)

real, intent(in) :: thil_upw(mza,mwa)
real, intent(in) :: vxe_upw (mza,mwa)
real, intent(in) :: vye_upw (mza,mwa)
real, intent(in) :: vze_upw (mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

integer :: jv, iv, iwn
integer :: k, ka, kbv, kp, ksw
integer :: npoly

real :: dts
real :: c6, c7, c8, c9, c10
real :: dirv, vmarv
real :: del_rhothil, vmt1, wmt1
real :: rad0_swtc, rad_swtc, topo_swtc

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55  ! wmc
real, parameter :: fr = .75  ! rho
real, parameter :: fp = .75  ! press

real, parameter :: pc2 = fp * cpocv

! Automatic arrays

real :: wmarw        (mza)
real :: del_wmarw    (mza)
real :: delex_wm     (mza)
real :: delex_rhothil(mza)
real :: del_wm       (mza)
real :: fwdel_wm     (mza)

real :: vmx_cor(mza)
real :: vmy_cor(mza)

real :: hflux_rho(mza)
real :: hflux_thil(mza)
real :: hflux_vxe(mza)
real :: hflux_vye(mza)
real :: hflux_vze(mza)

real :: vflux_thil(mza)
real :: vflux_vxe(mza)
real :: vflux_vye(mza)
real :: vflux_vze(mza)

real(r8) :: delex_rho(mza)
real(r8) :: rhothil  (mza)
real(r8) :: press_t  (mza)

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

! Loop over T levels

   do k = kbv,mza

      vmarv = dirv * vmcf(k,iv) * arv(k,iv)

! Sum horizontal advection fluxes over V faces      

      hflux_rho(k)  = hflux_rho(k)  + vmarv

      hflux_thil(k) = hflux_thil(k) + vmarv * thil_upv(k,iv)

      hflux_vxe(k)  = hflux_vxe(k)  + vmarv * vxe_upv(k,iv)

      hflux_vye(k)  = hflux_vye(k)  + vmarv * vye_upv(k,iv)

      hflux_vze(k)  = hflux_vze(k)  + vmarv * vze_upv(k,iv)

   enddo

enddo

! Coriolis force

if (icorflg == 1) then
   do k = ka, mza
      vmx_cor(k) =  omega2 * vye(k,iw) * rho(k,iw)
      vmy_cor(k) = -omega2 * vxe(k,iw) * rho(k,iw)
   enddo
else
   do k = ka, mza
      vmx_cor(k) = 0.0
      vmy_cor(k) = 0.0
   enddo
endif

! Loop over T levels

do k = ka,mza

! Prognostic changes from long-timestep tendencies, explicit advective fluxes,
! and Coriolis force

   delex_rho(k) = dts * (rhot(k,iw) &
      + volti(k,iw) * (hflux_rho(k) + wmarw(k-1) - wmarw(k)))

   delex_rhothil(k) = dts * (thilt(k,iw) + thil(k,iw) * rhot(k,iw) &
      + volti(k,iw) * (hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)))

   vmxet(k,iw) = volti(k,iw) * (hflux_vxe(k) + vflux_vxe(k-1) - vflux_vxe(k)) &
               + vmx_cor(k)
   vmyet(k,iw) = volti(k,iw) * (hflux_vye(k) + vflux_vye(k-1) - vflux_vye(k)) &
               + vmy_cor(k)
   vmzet(k,iw) = volti(k,iw) * (hflux_vze(k) + vflux_vze(k-1) - vflux_vze(k))

! RHOTHIL(t) and PRESS(t)

   rhothil(k) = rho(k,iw) * thil_s(k,iw)
   press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

! Compute current T cell momentum and store in temp array

   vmxe1(k) = vxe(k,iw) * rho(k,iw)
   vmye1(k) = vye(k,iw) * rho(k,iw)
   vmze1(k) = vze(k,iw) * rho(k,iw)
enddo

! Loop over W levels for update of DELEX_WM

do k = ka,mza-1

! Change in WM from EXPLICIT terms (long timestep tendency, 3 horizontal
! advective fluxes, 2 vertical advective fluxes, vertical pgf, gravity)

   delex_wm(k) = dts * (wmt(k,iw) &

      + dzim(k) * (press_t(k) - press_t(k+1) &

      - grav * (dzt_top(k) * rho(k,iw) + dzt_bot(k+1) * rho(k+1,iw))) &

! In the following, do we want to use unequal weights between k and k+1 T levels?

      + wnxo2(iw) * (vmxet(k,iw) + vmxet(k+1,iw)) &
      + wnyo2(iw) * (vmyet(k,iw) + vmyet(k+1,iw)) &
      + wnzo2(iw) * (vmzet(k,iw) + vmzet(k+1,iw)))

enddo

c6  = dts * .5 * fw
c7  = dts * .25 * fw
c8  = dts * pc2
c9  =-dts * fr * grav
c10 = dts * fw

! Fill matrix coefficients for implicit update of WM

do k = ka,mza
   kp = min(k+1,mza)
   b1(k)  = wc(k,iw) + wc(k-1,iw)            ! T pts
   b2(k)  = thil_s(k,iw) + thil_s(kp,iw)    ! W pts
   b3(k)  = 2. / (volt(k,iw) + volt(kp,iw)) ! W pts [b3 replaces volwi]
   b5(k)  = press_t(k) / rhothil(k)          ! T pts
   b6(k)  = c6 * volti(k,iw)                 ! T pts
   b10(k) = c10 * volti(k,iw)                ! T pts
enddo
b2(ka-1) = b2(ka)
b3(ka)   = 1. / (volt(ka,iw) + .5 * volt(ka+1,iw))

do k = ka,mza-1

   b7(k)  = c7     * b3(k)        ! W pts
   b8(k)  = c8     * dzim(k)      ! W pts
   b9(k)  = c9     * dzim(k)      ! W pts
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

! Update velocity tendencies at T due to change in fluxes

   vmxet(k,iw) = vmxet(k,iw) + volti(k,iw) * (vflux_vxe(k-1) - vflux_vxe(k))
   vmyet(k,iw) = vmyet(k,iw) + volti(k,iw) * (vflux_vye(k-1) - vflux_vye(k))
   vmzet(k,iw) = vmzet(k,iw) + volti(k,iw) * (vflux_vze(k-1) - vflux_vze(k))

! Update velocity in T cells, and add half-forward value to scalar arrays

   vxe1(k) = (vmxe1(k) + dts * vmxet(k,iw)) / rho(k,iw)
   vye1(k) = (vmye1(k) + dts * vmyet(k,iw)) / rho(k,iw)
   vze1(k) = (vmze1(k) + dts * vmzet(k,iw)) / rho(k,iw)

   vxesc(k,iw) = vxesc(k,iw) + .5 * (vxe(k,iw) + vxe1(k))
   vyesc(k,iw) = vyesc(k,iw) + .5 * (vye(k,iw) + vye1(k))
   vzesc(k,iw) = vzesc(k,iw) + .5 * (vze(k,iw) + vze1(k))
enddo

! Contribution of surface W face to v[xyz]e2

wmt1 = wnx(iw) * vxe1(ka) + wny(iw) * vye1(ka) + wnz(iw) * vze1(ka)

vxe2(1,iw) = wmt1 * wnxo2(iw)
vye2(1,iw) = wmt1 * wnyo2(iw)
vze2(1,iw) = wmt1 * wnzo2(iw)

! Zero out v[xyz]e2 values about lpw

if (lve2(iw) > 1) then
   do ksw = 2, lve2(iw)
      vxe2(ksw,iw) = 0.
      vye2(ksw,iw) = 0.
      vze2(ksw,iw) = 0.
   enddo
endif

! Loop over adjacent V faces

do jv = 1, npoly
   iv = itab_w(iw)%iv(jv)

! Project vxe1, vye1, vze1 onto V faces that are below ground, and then
! project back to vxe2, vye2, vze2

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

wc(ka-1,iw) = 0.
wc(mza ,iw) = 0.

return
end subroutine prog_wrt_begs

!============================================================================

subroutine prog_v_begs(iv,vmxet,vmyet,vmzet,div2d)

use mem_tend,    only: vmt
use mem_ijtabs,  only: itab_v, itab_w
use mem_basic,   only: vc, press, vmp, vmc, rho
use misc_coms,   only: io6, dtsm, initial, mdomain, u01d, v01d, dn01d, &
                       deltax, nxp
use consts_coms, only: erad, eradi, gravo2
use mem_grid,    only: lpv, lpw, volt, xev, yev, zev, &
                       vnxo2, vnyo2, vnzo2, mza, mva, mwa, dniv, arw0, dnu
use mem_rayf,    only: dorayf, rayf_cof, vc03d, dn03d, krayf_bot, &
                       dorayfdiv, krayfdiv_bot, rayf_cofdiv
use oname_coms,  only: nl

implicit none

integer, intent(in) :: iv

real, intent(in) :: vmxet(mza,mwa)
real, intent(in) :: vmyet(mza,mwa)
real, intent(in) :: vmzet(mza,mwa)
real, intent(in) :: div2d(mza,mwa)

integer :: jv,ivn,k,kb,npoly

integer :: iw1,iw2

real :: sum1,sum2,vmp_eqdiv

real :: dts,raxis,uv01dr,uv01dx,uv01dy,uv01dz,vcref
real :: fracx, rayfx

! Automatic arrays

real :: vmt_rayf(mza)
real :: pgf     (mza)

! Extract neighbor indices and coefficients for this point in the U stencil

iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

dts = dtsm(itab_v(iv)%mrlv)

kb = lpv(iv)

vmt_rayf(:) = 0.

if (dorayf) then

! FOR HORIZONTAL HOMOGENEOUS CASE, APPLY
! RAYLEIGH FRICTION DIRECTLY TO VMC

   if (initial == 1) then      ! HHI case

!! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN
!!    fracx = abs(xev(iv)) / (real(nxp-1) * .866 * deltax)
!!    rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
!!    rayfx = 0.   ! Default: no extra RAYF
!!    do k = kb, mza
!!        vmt_rayf(k) = max(rayf_cof(k),rayfx) * dn03d(k) *  (vc03d(k,iv) - vc(k,iv))
!!     enddo
!! END SPECIAL

! Vertical loop over V points
      
      do k = krayf_bot, mza
         vmt_rayf(k) = rayf_cof(k) * dn03d(k,iv) * (vc03d(k,iv) - vc(k,iv))
      enddo

   else                     ! LHI/VARI case
   endif

endif ! (dorayf)
   
if (dorayfdiv) then

! FOR RUNS VARYING LATITUDINALLY AND/OR LONGITUDINALLY,
! PERFORM HORIZONTAL DIVERGENCE DAMPING

! Vertical loop over V points

   do k = krayfdiv_bot, mza

! Divergence in IW1 excluding IV

      sum1 = div2d(k,iw1) - vmp(k,iv) * dnu(iv)

! Divergence in IW2 excluding IV

      sum2 = div2d(k,iw2) + vmp(k,iv) * dnu(iv)

! VMP value that would equalize horizontal divergence in IW1 and IW2 cells
! (assuming no blockage by topography)

      vmp_eqdiv = (arw0(iw1) * sum2 - arw0(iw2) * sum1) &
           / (dnu(iv) * (arw0(iw1) + arw0(iw2)))

      vmt_rayf(k) = vmt_rayf(k) + rayf_cofdiv(k) * (vmp_eqdiv - vmp(k,iv))

   enddo

endif ! (dorayfdiv)

! Vertical loop over V points

do k = kb,mza
   pgf(k) = dniv(iv) * (press(k,iw1) - press(k,iw2))
enddo

! For shallow water test cases 2 & 5, rho & press are
! interpreted as water depth & height

if (nl%test_case == 2 .or. nl%test_case == 5) then
   do k = kb,mza
      pgf(k) = pgf(k) * gravo2 * (rho(k,iw1) + rho(k,iw2))
   enddo
endif

! Vertical loop over V points

do k = kb,mza

! Update VM from long timestep tendencies, advection, and pressure gradient force

   vmc(k,iv) = vmc(k,iv) + dts * (vmt(k,iv) + vmt_rayf(k) + pgf(k) &
             + vnxo2(iv) * (vmxet(k,iw1) + vmxet(k,iw2)) &
             + vnyo2(iv) * (vmyet(k,iw1) + vmyet(k,iw2)) &
             + vnzo2(iv) * (vmzet(k,iw1) + vmzet(k,iw2)) )

   vc(k,iv) = 2.0 * vmc(k,iv) / (rho(k,iw1) + rho(k,iw2))

enddo

return
end subroutine prog_v_begs
