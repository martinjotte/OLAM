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
subroutine prog_wrtv(vmsc,wmsc,alpha_press,rhot)

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

use mem_ijtabs, only: jtab_v, jtab_w, itab_v, istp, itab_w, jtab_m, itab_m, &
                      mrl_begl, mrl_begs, mrl_ends, mrl_endl
use mem_basic,  only: rho, thil, theta, wc, press, wmc, vmp, vmc, vp, vc, &
                      sh_w, sh_v
use mem_grid,   only: mza, mma, mva, mwa, lpm, lpv, lcv, lpw, &
                      zt, zm, dzim, zfacit, dzm, dzt, dnv, dnu, arm0, arv, &
                      vnx, vny, vnz, volt, volwi
use mem_tend,   only: vmt
use mem_turb,   only: vels
use misc_coms,  only: io6, iparallel, time8, dtlm
use consts_coms, only: cpocv, pc1, rdry, rvap
use massflux,   only: diagnose_uc

use olam_mpi_atm, only: mpi_send_vf, mpi_recv_vf, &
                        mpi_send_w, mpi_recv_w, &
                        mpi_send_v

use oplot_coms,  only: op

!$ use omp_lib

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

real, intent(inout) :: alpha_press(mza,mwa)
real, intent(in)    :: rhot       (mza,mwa)

integer :: j,iv,iw,k,ka,kb,iwp,mrl,ivp

integer :: iw1,iw2

real :: dts,dts2

! automatic arrays

integer :: iwdepv(mza,mva)
integer :: kdepw(mza,mwa)

real :: dxps_v(mza,mva)
real :: dyps_v(mza,mva)
real :: dzps_v(mza,mva)

real :: dxps_w(mza,mwa)
real :: dyps_w(mza,mwa)
real :: dzps_w(mza,mwa)

real :: vmcf(mza,mva) ! Time-extrapolated vmc
real :: vcf(mza,mva)  ! Time-extrapolated vc
real :: wmcf(mza,mwa) ! Time-interpolated wmc

real :: vxe(mza,mwa) ! XE velocity component at T point
real :: vye(mza,mwa) ! YE velocity component at T point
real :: vze(mza,mwa) ! ZE velocity component at T point

real :: vmxet(mza,mwa) ! XE velocity tendency component at T point
real :: vmyet(mza,mwa) ! YE velocity tendency component at T point
real :: vmzet(mza,mwa) ! ZE velocity tendency component at T point

real :: gxps_thil(mza,mwa)
real :: gyps_thil(mza,mwa)
real :: gzps_thil(mza,mwa)

real :: gxps_vxe(mza,mwa)
real :: gyps_vxe(mza,mwa)
real :: gzps_vxe(mza,mwa)

real :: gxps_vye(mza,mwa)
real :: gyps_vye(mza,mwa)
real :: gzps_vye(mza,mwa)

real :: gxps_vze(mza,mwa)
real :: gyps_vze(mza,mwa)
real :: gzps_vze(mza,mwa)

real :: thil_s(mza,mwa)

real :: vortp(mza,mma)

integer :: iv1,iv2,iv3,iv4,im,npoly,jv,im1,im2,im3,im4,im5,im6
real :: c0,c1,c2,vort_big1,vort_big2
real :: arm0i,tvort

! Save copy of thil

thil_s(:,:) = thil(:,:)

vmcf(:,1) = 0.
vcf(:,1) = 0.
wmcf(:,1) = 0.

! Horizontal loop over V/N columns

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(k) 
do iv = 2,mva
!----------------------------------------------------------------------
call qsub('V',iv)

! Extrapolate VM and V to time T + 1/2; update VMP

   do k = 1,mza-1
      vmcf(k,iv) = (1.5 * vmc(k,iv) - 0.5 * vmp(k,iv))
      vcf (k,iv) = (1.5 * vc (k,iv) - 0.5 * vp (k,iv))

      vmsc(k,iv) = vmsc(k,iv) + vmcf(k,iv)
      vmp (k,iv) = vmc(k,iv)
   enddo

enddo
!$omp end parallel do 
endif
call rsub('Va',16)

! Diagnose 3D velocity at T points for BEGS using half-future vcf and current wc

mrl = mrl_begs(istp)
if (mrl > 0) then
   call vel_t3d(mrl,vcf,wc,vxe,vye,vze)
endif

! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
! (Example of how to plot "external" field; one not available in module memory)

if (mod(real(time8),op%frqplt) < dtlm(1) .and. istp == 900) then

   allocate (op%extfld(mza,mwa))
   op%extfld(:,:) = vxe(:,:)
   op%extfldname = 'VXE'
   call plot_fields(11)
   deallocate (op%extfld)

endif
! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

!---------------------------------------------------------
! Parallel send/recv of VXE, VYE, VZE
!---------------------------------------------------------

! Diagnose advective donor point location for all V and W faces

mrl = mrl_begs(istp)
if (mrl > 0) then
   call donorpoint3d(0,mrl,vcf,wc,vxe,vye,vze, &
      iwdepv,kdepw,dxps_v,dyps_v,dzps_v,dxps_w,dyps_w,dzps_w)
endif

! Diagnose 3D velocity at T points for BEGS using current vc and wc

mrl = mrl_begs(istp)
if (mrl > 0) then
   call vel_t3d(mrl,vc,wc,vxe,vye,vze)
endif
! Compute turbulent mixing coefficients

mrl = mrl_begl(istp)
if (mrl > 0) then
   call turb_k_hex(mrl,vxe,vye,vze)
endif

!---------------------------------------------------------
! Parallel send/recv of VXE, VYE, VZE
!---------------------------------------------------------

if (mrl > 0 .and. iparallel == 1) then
   call mpi_send_w('K')  ! Send K's
endif

!----------------------------------------------------------------------
if (mrl > 0) then
!$omp parallel do private(iw,ka) 
do j = 1,jtab_w(20)%jend(mrl); iw = jtab_w(20)%iw(j)
!----------------------------------------------------------------------
   ka = lpw(iw)

   vels(iw) = sqrt(vxe(ka,iw) ** 2 + vye(ka,iw) ** 2 + vze(ka,iw) ** 2)

enddo

call surface_turb_flux(mrl)
endif

if (mrl > 0 .and. iparallel == 1) then
   call mpi_recv_w('K')  ! Recv K's
endif

if (mrl > 0) then
   call thiltend_long(mrl,rhot)
endif

! Horizontal loop over W columns for BEGL

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw) 
do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Evaluate alpha coefficient for pressure

   alpha_press(:,iw) = pc1 * (((1. - sh_w(:,iw)) * rdry + sh_v(:,iw) * rvap) &
                     * theta(:,iw) / thil(:,iw)) ** cpocv

! Long timestep tendencies for WM (turbulent mixing), 
! and preparation for horizontal turbulent fluxes of VM

   vmxet(:,iw) = 0.
   vmyet(:,iw) = 0.
   vmzet(:,iw) = 0.

   call prog_wrt_begl(iw,vxe,vye,vze,vmxet,vmyet,vmzet)
   
enddo
!$omp end parallel do
endif
call rsub('Wa',16)

!---------------------------------------------------------
! Parallel send/recv of vmxet, vmyet, vmzet
!---------------------------------------------------------

! Horizontal loop over M/P columns for BEGL; Diagnose vertical vorticity
! in preparation for horizontal filter

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(im,npoly,kb,jv,iv,k,arm0i)
do j = 1,jtab_m(3)%jend(mrl); im = jtab_m(3)%im(j)
!----------------------------------------------------------------------
call qsub('M',im)

   npoly = itab_m(im)%npoly

   kb = lpm(im)

   vortp(:,im) = 0.

! Loop over V neighbors to evaluate circulation around M (at time T)

   do jv = 1,npoly
      iv = itab_m(im)%iv(jv)

      if (itab_v(iv)%im(2) == im) then

         do k = kb,mza-1
            vortp(k,im) = vortp(k,im) + vc(k,iv) * dnv(iv)
         enddo

      else

         do k = kb,mza-1
            vortp(k,im) = vortp(k,im) - vc(k,iv) * dnv(iv)
         enddo

      endif
   enddo

! Convert circulation to relative vertical vorticity at M 
! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
! divide their quotient by zfact)

   arm0i = 1. / arm0(im)

   do k = kb,mza-1
      vortp(k,im) = vortp(k,im) * arm0i * zfacit(k)
   enddo

enddo
!$omp end parallel do 
endif
call rsub('M',3)

tvort = 3600.
c0 = .125 / tvort
!c0 = 0.

! Horizontal loop over V columns for PROG_V_BEGL

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iv,iw1,iw2,im1,im2,im3,im4,im5,im6,iv1,iv2,iv3,iv4, &
!$omp                     kb,k,c1,c2,vort_big1,vort_big2) 
do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)
iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
call qsub('V',iv)

   kb = lpv(iv)

! Vertical loop over T levels

   do k = kb,mza-1

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

   do k = kb,mza-1
   
      vort_big1 = .3333333 * (vortp(k,im2) + vortp(k,im3) + vortp(k,im4))
      vort_big2 = .3333333 * (vortp(k,im1) + vortp(k,im5) + vortp(k,im6))
   
      vmt(k,iv) = vmt(k,iv) + (rho(k,iw1) + rho(k,iw2)) &
         * ((vort_big1 - vortp(k,im1)) * c1 + (vort_big2 - vortp(k,im2)) * c2)

   enddo
              
enddo
endif
call rsub('V',12)

! [Now, we can reuse vmxet, vmyet, vmzet arrays]

! Evaluate T3D gradients of THIL, VXE, VYE, and VZE for BEGS

mrl = mrl_begs(istp)
if (mrl > 0) then
   call grad_t3d(mrl,thil,gxps_thil,gyps_thil,gzps_thil)
   call grad_t3d(mrl,vxe,gxps_vxe,gyps_vxe,gzps_vxe)
   call grad_t3d(mrl,vye,gxps_vye,gyps_vye,gzps_vye)
   call grad_t3d(mrl,vze,gxps_vze,gyps_vze,gzps_vze)
endif

!---------------------------------------------------------
! Parallel send/recv of T3D gradients (12 of them)
!---------------------------------------------------------

! Main loop over W columns for updating WM, WC, RHO, THIL, and PRESS

call psub()
!----------------------------------------------------------------------
mrl = mrl_begs(istp)
if (mrl > 0) then
!$omp parallel do private(iw) 
do j = 1,jtab_w(19)%jend(mrl); iw = jtab_w(19)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call prog_wrt_begs(iw,vmcf,wmsc,alpha_press,rhot, &
                      vxe,vye,vze,vmxet,vmyet,vmzet,kdepw,iwdepv, &
                      dxps_v,dyps_v,dzps_v,dxps_w,dyps_w,dzps_w, &
                      gxps_thil,gyps_thil,gzps_thil, &
                      gxps_vxe,gyps_vxe,gzps_vxe, &
                      gxps_vye,gyps_vye,gzps_vye, &
                      gxps_vze,gyps_vze,gzps_vze,thil_s)

enddo
endif
call rsub('Wa',19)

!---------------------------------------------------------
! Parallel send/recv of vmxet, vmyet, vmzet
!---------------------------------------------------------

! Horizontal loop over V points to update VMC

call psub()
!----------------------------------------------------------------------
mrl = mrl_ends(istp)
if (mrl > 0) then
!$omp parallel do private(iv) 
do j = 1,jtab_v(16)%jend(mrl); iv = jtab_v(16)%iv(j)
!----------------------------------------------------------------------
call qsub('V',iv)

   call prog_v_begs(iv,vmxet,vmyet,vmzet)

enddo
endif
call rsub('Va',16)

call diagnose_uc()

if (iparallel == 1) then
   call mpi_send_v('V')  ! Send V group
endif

return
end subroutine prog_wrtv

!=========================================================================

subroutine prog_wrt_begl(iw,vxe,vye,vze,vmxet,vmyet,vmzet)

! This version includes turbulent fluxes of VXE, VYE, VZE 
! through both V and W faces

! All diffusive tendencies are evaluated at T points

use mem_tend,    only: thilt, wmt
use mem_ijtabs,  only: istp, itab_w, mrl_begl, mrl_begr, mrl_begs, mrl_endr
use mem_basic,   only: wmc, rho, thil, wc, vc, theta, press
use misc_coms,   only: io6, initial, dn01d, th01d, &
                       deltax, nxp, mdomain, time8, dtlm
use consts_coms, only: gravo2, grav
use mem_grid,    only: mza, mva, mwa, lpv, lpw, lsw, arv, arw, arw0, &
                       dnu, dniv, dzim, dzt, dzit, volt, volti, volwi, &
                       xew, zm, unx, uny, vnx, vny, wnx, wny, wnz, &
                       glatw, glonw, topw
use mem_turb,    only: hkm, vkm, sflux_w, vkm_sfc
use mem_rayf,    only: rayfw_distim, rayf_cofw, rayf_distim, rayf_cof
use massflux,    only: tridiffo
use mem_rayf,    only: rayfw_distim, rayf_cofw, rayf_distim, rayf_cof

implicit none

integer, intent(in) :: iw

real, intent(in) :: vxe(mza,mwa)
real, intent(in) :: vye(mza,mwa)
real, intent(in) :: vze(mza,mwa)

real, intent(inout) :: vmxet(mza,mwa)
real, intent(inout) :: vmyet(mza,mwa)
real, intent(inout) :: vmzet(mza,mwa)

integer :: iv, iwn, k, kb, kbv, npoly, jv, kn, ivn

real :: dirv,arw0i

real :: fracx, rayfx

real :: fwv,fww,vproj,arvkodx,qdniv

real :: dtl,dtl2,dtli,hdtli,sflux

! Automatic arrays:

real :: akodz(mza),tmass(mza),dtomass(mza)
real :: vctr2a(mza),vctr2b(mza),vctr2c(mza)
real :: vctr3(mza),vctr5(mza),vctr6(mza),vctr7(mza)
real :: vctr8a(mza),vctr8b(mza),vctr8c(mza)
real :: vctr9a(mza),vctr9b(mza),vctr9c(mza)

real :: hdiff_w(mza)
real :: hdiff_vxe(mza)
real :: hdiff_vye(mza)
real :: hdiff_vze(mza)

kb = lpw(iw)

arw0i = 1. / arw0(iw)

! Initial computations for this column

dtl = dtlm(itab_w(iw)%mrlw)
dtl2 = dtl * 2.
dtli = 1. / dtl
hdtli = .5 * dtli
   
! Vertical loop over W levels

do k = kb,mza-2
   akodz(k) = arw(k,iw) * vkm(k,iw) * dzim(k)
enddo

akodz(kb-1) = 0.
akodz(mza-1) = 0.

vctr3(1:mza) = 0.

! Vertical loop over T levels 

do k = kb,mza-1

! Mass in t control volume; and its inverse times dtl

   tmass(k)  = rho(k,iw) * volt(k,iw)
   dtomass(k) = dtl / tmass(k)

! Distribution of surface flux over multiple levels in steep topography

   if (k <= lpw(iw) + lsw(iw) - 1) then
      vctr3(k) = (arw(k,iw) - arw(k-1,iw)) * vkm_sfc(iw) * dzim(k-1) * 2.
   endif

! Fill tri-diagonal matrix coefficients

   vctr5(k) = -dtomass(k) * akodz(k-1)
   vctr7(k) = -dtomass(k) * akodz(k)
   vctr6(k) = 1. - vctr5(k) - vctr7(k) + dtomass(k) * vctr3(k)

! Fill r.h.s. vectors

   vctr8a(k) = vxe(k,iw)
   vctr8b(k) = vye(k,iw)
   vctr8c(k) = vze(k,iw)

enddo

! Solve tri-diagonal matrix for each component

if (kb <= mza-1) then
   call tridiffo(mza,kb,mza-1,vctr5,vctr6,vctr7,vctr8a,vctr9a)
   call tridiffo(mza,kb,mza-1,vctr5,vctr6,vctr7,vctr8b,vctr9b)
   call tridiffo(mza,kb,mza-1,vctr5,vctr6,vctr7,vctr8c,vctr9c)
endif

! Now, vctr9 contains velocity(t+1) values

! Vertical loop over W levels

do k = kb,mza-2

! Compute internal vertical turbulent fluxes

   vctr2a(k) = akodz(k) * (vctr9a(k) - vctr9a(k+1))
   vctr2b(k) = akodz(k) * (vctr9b(k) - vctr9b(k+1))
   vctr2c(k) = akodz(k) * (vctr9c(k) - vctr9c(k+1))
   
enddo

! Set bottom and top internal fluxes to zero

vctr2a(kb-1) = 0.
vctr2b(kb-1) = 0.
vctr2c(kb-1) = 0.

vctr2a(mza-1) = 0.
vctr2b(mza-1) = 0.
vctr2c(mza-1) = 0.

! Number of edges of this IW polygon

npoly = itab_w(iw)%npoly

hdiff_vxe(1:mza) = 0.
hdiff_vye(1:mza) = 0.
hdiff_vze(1:mza) = 0.

! Loop over V neighbors of this W cell

do jv = 1,npoly
   iv  = itab_w(iw)%iv(jv)
   iwn = itab_w(iw)%iw(jv)

   qdniv = .25 * dniv(iv)

! Vertical loop over T levels

   do k = kb,mza-1

! Horizontal diffusive flux coefficient

      arvkodx = qdniv * arv(k,iv) * (hkm(k,iwn) + hkm(k,iw))

! Compute and sum horizontal diffusive flux across this V neighbor

      hdiff_vxe(k) = hdiff_vxe(k) + arvkodx * (vxe(k,iwn) - vxe(k,iw))
      hdiff_vye(k) = hdiff_vye(k) + arvkodx * (vye(k,iwn) - vye(k,iw))
      hdiff_vze(k) = hdiff_vze(k) + arvkodx * (vze(k,iwn) - vze(k,iw))

   enddo

enddo

! Vertical loop over T levels

do k = kb,mza-1

   vmxet(k,iw) = volti(k,iw) &
               * (vctr2a(k-1) - vctr2a(k) - vctr3(k) * vctr9a(k) + hdiff_vxe(k))
   vmyet(k,iw) = volti(k,iw) &
               * (vctr2b(k-1) - vctr2b(k) - vctr3(k) * vctr9b(k) + hdiff_vye(k))
   vmzet(k,iw) = volti(k,iw) &
               * (vctr2c(k-1) - vctr2c(k) - vctr3(k) * vctr9c(k) + hdiff_vze(k))

enddo

! Vertical loop over W levels

do k = kb,mza-2

! Update WM tendency from turbulent fluxes

   wmt(k,iw) = wmt(k,iw) + .5 &
             * (wnx(iw) * (vmxet(k,iw) + vmxet(k+1,iw)) &
             *  wny(iw) * (vmyet(k,iw) + vmyet(k+1,iw)) &
             *  wnz(iw) * (vmzet(k,iw) + vmzet(k+1,iw)))

enddo

! RAYLEIGH FRICTION ON WM

if (rayfw_distim > 1.e-6) then
   fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
   rayfx = .2 * (-2. + 3. * fracx) * rayf_cofw(mza-2)
   rayfx = 0.   ! Default: no extra RAYF
   do k = kb,mza-2
      wmt(k,iw) = wmt(k,iw) - max(rayf_cofw(k),rayfx) * wmc(k,iw)
   enddo
endif

! RAYLEIGH FRICTION ON THIL

if (rayf_distim > 1.e-6) then
   if (initial == 1) then   ! HHI case
      fracx = abs(xew(iw)) / (real(nxp-1) * .866 * deltax) ! ENDS OF CYC DOMAIN
      rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
      rayfx = 0.   ! Default: no extra RAYF
      do k = kb,mza-1
! Form based on theta alone
         thilt(k,iw) = thilt(k,iw) + max(rayf_cof(k),rayfx)  & 
            * dn01d(k) * (th01d(k) - theta(k,iw))
      enddo
   else                     ! LHI/VARI case
! Need implementation for LHI/VARI (use vartp for merid. variation?)
   endif

endif ! (rayf_distim > 1.e-6)

return
end subroutine prog_wrt_begl

!=========================================================================

subroutine prog_wrt_begs(iw,vmcf,wmsc,alpha_press,rhot, &
                         vxe,vye,vze,vmxet,vmyet,vmzet,kdepw,iwdepv, &
                         dxps_v,dyps_v,dzps_v,dxps_w,dyps_w,dzps_w, &
                         gxps_thil,gyps_thil,gzps_thil, &
                         gxps_vxe,gyps_vxe,gzps_vxe, &
                         gxps_vye,gyps_vye,gzps_vye, &
                         gxps_vze,gyps_vze,gzps_vze,thil_s)

use mem_tend,    only: thilt, wmt
use mem_ijtabs,  only: itab_w
use mem_basic,   only: wmc, rho, thil, wc, vc, theta, press
use misc_coms,   only: io6, dtsm, initial, dn01d, th01d, &
                       deltax, nxp, mdomain, time8
use consts_coms, only: cpocv, gravo2, grav, omega2
use mem_grid,    only: mza, mva, mwa, lpv, lpw, &
                       arv, arw, volt, volti, volwi, dzm, dzim, dzt, xew, zm, &
                       unx, uny, vnx, vny, wnx, wny, wnz, glatw, glonw
use massflux,    only: tridiffo

implicit none

integer, intent(in) :: iw

real, intent(in) :: vmcf(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)
real, intent(in) :: alpha_press(mza,mwa)
real, intent(in) :: rhot(mza,mwa)

integer, intent(in) :: iwdepv(mza,mva)
integer, intent(in) :: kdepw(mza,mwa)

real, intent(in) :: vxe(mza,mwa)
real, intent(in) :: vye(mza,mwa)
real, intent(in) :: vze(mza,mwa)

real, intent(out) :: vmxet(mza,mwa)
real, intent(out) :: vmyet(mza,mwa)
real, intent(out) :: vmzet(mza,mwa)

real, intent(in) :: dxps_v(mza,mva)
real, intent(in) :: dyps_v(mza,mva)
real, intent(in) :: dzps_v(mza,mva)

real, intent(in) :: dxps_w(mza,mwa)
real, intent(in) :: dyps_w(mza,mwa)
real, intent(in) :: dzps_w(mza,mwa)

real, intent(in) :: gxps_thil(mza,mwa)
real, intent(in) :: gyps_thil(mza,mwa)
real, intent(in) :: gzps_thil(mza,mwa)

real, intent(in) :: gxps_vxe(mza,mwa)
real, intent(in) :: gyps_vxe(mza,mwa)
real, intent(in) :: gzps_vxe(mza,mwa)

real, intent(in) :: gxps_vye(mza,mwa)
real, intent(in) :: gyps_vye(mza,mwa)
real, intent(in) :: gzps_vye(mza,mwa)

real, intent(in) :: gxps_vze(mza,mwa)
real, intent(in) :: gyps_vze(mza,mwa)
real, intent(in) :: gzps_vze(mza,mwa)

real, intent(in) :: thil_s(mza,mwa)

integer :: jv, iv, iwn, iwd

real :: dirv
real :: fwv,fww

integer :: k,ka,kd,kbv
integer :: k1,k2,k3
integer :: npoly,kn,ivn

real :: dts,dtso2,dts2,dts8,flux_rhothil

real :: c6,c7,c8,c9,c10
real :: fracx, rayfx

real :: cnum_w
real :: vproj
real :: del_rhothil
real :: hcnsclr

! Vertical implicit scheme weighting parameters

real, parameter :: fw = .55  ! wmc
real, parameter :: fr = .55  ! rho
real, parameter :: fp = .75  ! press

real, parameter :: pc2 = fp * cpocv

! Automatic arrays

real :: vmarv        (mza)
real :: wmarw        (mza)
real :: del_wmarw    (mza)
real :: delex_wm     (mza)
real :: delex_rhothil(mza)
real :: del_wm       (mza)
real :: fwdel_wm     (mza)
real :: wmtharw      (mza)
real :: del_wmtharw  (mza)
real :: wvertflx     (mza)
real :: hadv_rho     (mza)
real :: hadv_rhothil (mza)
real :: hadv_wm      (mza)
real :: thilw        (mza)
real :: vxew         (mza)
real :: vyew         (mza)
real :: vzew         (mza)

real :: hflux_rho(mza)
real :: hflux_thil(mza)
real :: hflux_vxe(mza)
real :: hflux_vye(mza)
real :: hflux_vze(mza)

real :: vflux_thil(mza)
real :: vflux_vxe(mza)
real :: vflux_vye(mza)
real :: vflux_vze(mza)

real(kind=8) :: delex_rho(mza)
real(kind=8) :: rhothil  (mza)
real(kind=8) :: press_t  (mza)

real :: b1(mza),b2(mza),b5(mza),b6(mza),b10(mza)
real :: b7(mza),b8(mza),b9(mza),b11(mza),b12(mza),b13(mza),b14(mza)
real :: b21(mza),b22(mza),b23(mza),b24(mza),b25(mza),b26(mza)
real :: b31(mza),b32(mza),b33(mza),b34(mza)

ka = lpw(iw)

dts = dtsm(itab_w(iw)%mrlw)
dtso2 = .5 * dts
dts2 = 2. * dts
dts8 = 8. * dts

! Set bottom & top vertical advective mass and heat fluxes to zero

wmarw(1:ka-1) = 0.
wmarw(mza-1) = 0.

vflux_thil(ka-1)  = 0.
vflux_thil(mza-1) = 0.

vflux_vxe(ka-1)  = 0.
vflux_vxe(mza-1) = 0.

vflux_vye(ka-1)  = 0.
vflux_vye(mza-1) = 0.

vflux_vze(ka-1)  = 0.
vflux_vze(mza-1) = 0.

! Loop over W levels

do k = ka,mza-2

! vertical mass flux for time level t

   wmarw(k) = wmc(k,iw) * arw(k,iw)

! Advective donor cell k level

   kd = kdepw(k,iw)

! mean scalar values in vertical fluxes

   thilw(k) = thil_s(kd,iw) &
            + dxps_w(k,iw) * gxps_thil(kd,iw) &
            + dyps_w(k,iw) * gyps_thil(kd,iw) &
            + dzps_w(k,iw) * gzps_thil(kd,iw)

   vxew(k) = vxe(kd,iw) &
           + dxps_w(k,iw) * gxps_vxe(kd,iw) &
           + dyps_w(k,iw) * gyps_vxe(kd,iw) &
           + dzps_w(k,iw) * gzps_vxe(kd,iw)

   vyew(k) = vye(kd,iw) &
           + dxps_w(k,iw) * gxps_vye(kd,iw) &
           + dyps_w(k,iw) * gyps_vye(kd,iw) &
           + dzps_w(k,iw) * gzps_vye(kd,iw)

   vzew(k) = vze(kd,iw) &
           + dxps_w(k,iw) * gxps_vze(kd,iw) &
           + dyps_w(k,iw) * gyps_vze(kd,iw) &
           + dzps_w(k,iw) * gzps_vze(kd,iw)

! vertical fluxes

   vflux_thil(k) = wmarw(k) * thilw(k)
   vflux_vxe(k)  = wmarw(k) * vxew(k)
   vflux_vye(k)  = wmarw(k) * vyew(k)
   vflux_vze(k)  = wmarw(k) * vzew(k)

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
   iv  = itab_w(iw)%iv(jv)
   iwn = itab_w(iw)%iw(jv)
   
   kbv = lpv(iv)

   dirv = itab_w(iw)%dirv(jv)

! Loop over T levels

   do k = kbv,mza-1

      iwd = iwdepv(k,iv)

      vmarv(k) = dirv * vmcf(k,iv) * arv(k,iv)

! Sum horizontal advection fluxes over V faces      

      hflux_rho(k) = hflux_rho(k) + vmarv(k)

      hflux_thil(k) = hflux_thil(k) &
                    + vmarv(k) * (thil_s(k,iwd) &
                    + dxps_v(k,iv) * gxps_thil(k,iwd) &
                    + dyps_v(k,iv) * gyps_thil(k,iwd) &
                    + dzps_v(k,iv) * gzps_thil(k,iwd))

      hflux_vxe(k) = hflux_vxe(k) &
                   + vmarv(k) * (vxe(k,iwd) &
                   + dxps_v(k,iv) * gxps_vxe(k,iwd) &
                   + dyps_v(k,iv) * gyps_vxe(k,iwd) &
                   + dzps_v(k,iv) * gzps_vxe(k,iwd))

      hflux_vye(k) = hflux_vye(k) &
                   + vmarv(k) * (vye(k,iwd) &
                   + dxps_v(k,iv) * gxps_vye(k,iwd) &
                   + dyps_v(k,iv) * gyps_vye(k,iwd) &
                   + dzps_v(k,iv) * gzps_vye(k,iwd))

      hflux_vze(k) = hflux_vze(k) &
                   + vmarv(k) * (vze(k,iwd) &
                   + dxps_v(k,iv) * gxps_vze(k,iwd) &
                   + dyps_v(k,iv) * gyps_vze(k,iwd) &
                   + dzps_v(k,iv) * gzps_vze(k,iwd))

   enddo

enddo

! Loop over T levels

do k = ka,mza-1

! Prognostic changes from long-timestep tendencies, explicit advective fluxes,
! and Coriolis force

   delex_rho(k) = dts * (rhot(k,iw) &
      + volti(k,iw) * (hflux_rho(k) + wmarw(k-1) - wmarw(k)))

   delex_rhothil(k) = dts * (thilt(k,iw) &
      + volti(k,iw) * (hflux_thil(k) + vflux_thil(k-1) - vflux_thil(k)))

   vmxet(k,iw) = volti(k,iw) * (hflux_vxe(k) + vflux_vxe(k-1) - vflux_vxe(k)) &
               + omega2 * vye(k,iw) * rho(k,iw)
   vmyet(k,iw) = volti(k,iw) * (hflux_vye(k) + vflux_vye(k-1) - vflux_vye(k)) &
               - omega2 * vxe(k,iw) * rho(k,iw)
   vmzet(k,iw) = volti(k,iw) * (hflux_vze(k) + vflux_vze(k-1) - vflux_vze(k))

! RHOTHIL(t) and PRESS(t)

   rhothil(k) = rho(k,iw) * thil_s(k,iw)
   press_t(k) = alpha_press(k,iw) * rhothil(k) ** cpocv

enddo

! Loop over W levels for update of DELEX_WM

do k = ka,mza-2

! Change in WM from EXPLICIT terms (long timestep tendency, 3 horizontal
! advective fluxes, 2 vertical advective fluxes, vertical pgf, gravity)

   delex_wm(k) = dts * (wmt(k,iw) &

      + dzim(k) * (press_t(k) - press_t(k+1) &

      - gravo2 * (dzt(k) * rho(k,iw) + dzt(k+1) * rho(k+1,iw))) &

! In the following, do we want to use unequal weights between k and k+1 T levels?

      + wnx(iw) * .5 * (vmxet(k,iw) + vmxet(k+1,iw)) &
      + wny(iw) * .5 * (vmyet(k,iw) + vmyet(k+1,iw)) &
      + wnz(iw) * .5 * (vmzet(k,iw) + vmzet(k+1,iw)))

enddo

!!!!!!!!!!!!!!!!!!!!!special
   go to 55
!!!!!!!!!!!!!!!!!!!!!! end special

! Change in WM(ka) from EXPLICIT terms (3 horizontal advective fluxes at ka-1)

delex_wm(ka) = delex_wm(ka) + dts * volwi(ka,iw) * .25 * hadv_wm(ka-1)

!!!!!!!!!!!!!! special
55 continue
!!!!!!!!!!!!!!!!end special

c6  = dts * .5 * fw
c7  = dts * .25 * fw
c8  = dts * fp * cpocv
c9  = dts * (-.5) * fr * grav
c10 = dts * fw

! Fill matrix coefficients for implicit update of WM

do k = ka,mza-1
   b1(k)  = wc(k,iw) + wc(k-1,iw)     ! T pts
   b2(k)  = thil_s(k,iw) + thil_s(k+1,iw) ! W pts
   b5(k)  = press_t(k) / rhothil(k)   ! T pts
   b6(k)  = c6 * volti(k,iw)          ! T pts
   b10(k) = c10 * volti(k,iw)         ! T pts
enddo
b2(ka-1) = b2(ka)

do k = ka,mza-2

   b7(k)  = c7     * volwi(k,iw) ! W pts
   b8(k)  = c8     * dzim(k)     ! W pts
   b9(k)  = c9     * dzim(k)     ! W pts
   b11(k) = b8(k)  * b5(k)       ! W pts
   b12(k) = b8(k)  * b5(k+1)     ! W pts
   b13(k) = b9(k)  * dzt(k)      ! W pts
   b14(k) = b9(k)  * dzt(k+1)    ! W pts

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

if (ka <= mza-2) then
   call tridiffo(mza,ka,mza-2,b31,b32,b33,b34,del_wm)
endif

! Vertical loop over W points

do k = ka,mza-2

! Change in vertical momentum from t to t+fw

   fwdel_wm(k) = fw * del_wm(k)

! Add vertical momentum at (t+fw) to array for long-timestep scalar transport

   wmsc(k,iw) = wmsc(k,iw) + wmc(k,iw) + fwdel_wm(k)

! Change in vertical fluxes from t to t + fw using del_wm value

   del_wmarw(k)  = fwdel_wm(k) * arw(k,iw)
   vflux_thil(k) = del_wmarw(k) * thilw(k)
   vflux_vxe(k)  = del_wmarw(k) * vxew(k)
   vflux_vye(k)  = del_wmarw(k) * vyew(k)
   vflux_vze(k)  = del_wmarw(k) * vzew(k)

enddo

! Set mass-flux-change array to 0 at top & bottom (other flux arrays already 0)

del_wmarw(ka-1)  = 0.
del_wmarw(mza-1) = 0.

! Vertical loop over T points

do k = ka,mza-1

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

enddo

! Vertical loop over W points

do k = ka,mza-2

! Update WMC and WC to future value (at t+1) due to change in WM (del_wm)

   wmc(k,iw) = wmc(k,iw) + del_wm(k)
   wc(k,iw) = wmc(k,iw) * 2. * dzm(k) &
            / (dzt(k+1) * rho(k,iw) + dzt(k) * rho(k+1,iw))
enddo

! Set top & bottom values of WC

wc(1:ka-1,iw) = wc(ka,iw)
wc(mza-1,iw) = 0.

return
end subroutine prog_wrt_begs

!============================================================================

subroutine prog_v_begs(iv,vmxet,vmyet,vmzet)

use mem_tend,    only: vmt
use mem_ijtabs,  only: itab_v, itab_w
use mem_basic,   only: vp, vc, press, vmp, vmc, rho
use misc_coms,   only: io6, dtsm, initial, mdomain, u01d, v01d, dn01d, &
                       deltax, nxp
use consts_coms, only: erad, eradi
use mem_grid,    only: lpv, lcv, lpw, volt, aru, volvi, xev, yev, zev, &
                       unx, uny, vnx, vny, vnz, mza, mva, mwa, dniv, arw0, dnu
use mem_rayf,    only: rayf_distim, rayf_cof

implicit none

integer, intent(in) :: iv

real, intent(in) :: vmxet(mza,mwa)
real, intent(in) :: vmyet(mza,mwa)
real, intent(in) :: vmzet(mza,mwa)

integer :: jv,ivn,k,kb,npoly

integer :: iw1,iw2

real :: sum1,sum2,vmp_eqdiv

real :: dts,raxis,uv01dr,uv01dx,uv01dy,uv01dz,vcref,vc2
real :: fracx, rayfx, div1, div2

! Automatic array

real :: vmt_rayf(mza)

! Extract neighbor indices and coefficients for this point in the U stencil

iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

dts = dtsm(itab_v(iv)%mrlv)

kb = lpv(iv)

! RAYLEIGH FRICTION ON VMC

vmt_rayf(:) = 0.

if (rayf_distim > 1.e-6) then

! Vertical loop over V points

   do k = kb,mza-1

      if (initial == 1) then      ! HHI case

! Must rotate reference wind to local VC orientation

         if (mdomain <= 1) then  ! Model uses "earth" coordinates
            raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis
            
            if (raxis > 1.e3) then
               uv01dr = -v01d(k) * zev(iv) / erad  ! radially outward from axis

               uv01dx = (-u01d(k) * yev(iv) + uv01dr * xev(iv)) / raxis 
               uv01dy = ( u01d(k) * xev(iv) + uv01dr * yev(iv)) / raxis 
               uv01dz =   v01d(k) * raxis / erad 

               vcref = uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv)
            else
               vcref = 0.
            endif
         else
            vcref = u01d(k) * vnx(iv) + v01d(k) * vny(iv)
         endif

! SPECIAL - EXTRA RAYF AT ENDS OF CYCLIC DOMAIN

         fracx = abs(xev(iv)) / (real(nxp-1) * .866 * deltax)
         rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza-1)
         rayfx = 0.   ! Default: no extra RAYF
         
! END SPECIAL

         vmt_rayf(k) = max(rayf_cof(k),rayfx) * dn01d(k) * (vcref - vc(k,iv))

      else                     ! LHI/VARI case

! HORIZONTAL DIVERGENCE DAMPING

! Divergence in IW1 excluding IV

         npoly = itab_w(iw1)%npoly
         sum1 = 0.
   
         do jv = 1,npoly
            ivn = itab_w(iw1)%iv(jv)
      
            if (ivn /= iv) &
               sum1 = sum1 - itab_w(iw1)%dirv(jv) * vmp(k,ivn) * dnu(ivn)
         enddo

! Divergence in IW2 excluding IV

         npoly = itab_w(iw2)%npoly
         sum2 = 0.
   
         do jv = 1,npoly
            ivn = itab_w(iw2)%iv(jv)
      
            if (ivn /= iv) &
               sum2 = sum2 - itab_w(iw2)%dirv(jv) * vmp(k,ivn) * dnu(ivn)
         enddo

! VMP value that would equalize horizontal divergence in IW1 and IW2 cells
! (assuming no blockage by topography)

         vmp_eqdiv = (arw0(iw1) * sum2 - arw0(iw2) * sum1) &
                   / (dnu(iv) * (arw0(iw1) + arw0(iw2)))

         vmt_rayf(k) = rayf_cof(k) * (vmp_eqdiv - vmp(k,iv))

         ! Need implementation for LHI/VARI (use vartp for merid. variation?)
         ! Ok to do this in veltend_long???

      endif

   enddo

endif ! (rayf_distim > 1.e-6)

! Vertical loop over V points

do k = kb,mza-1

! Update VM from long timestep tendencies, advection, and pressure gradient force

   vmc(k,iv) = vmc(k,iv) + dts * (vmt(k,iv) + vmt_rayf(k) &
             + dniv(iv) * (press(k,iw1) - press(k,iw2)) &
             + .5 * (vnx(iv) * (vmxet(k,iw1) + vmxet(k,iw2)) &
                   + vny(iv) * (vmyet(k,iw1) + vmyet(k,iw2)) &
                   + vnz(iv) * (vmzet(k,iw1) + vmzet(k,iw2))))

   vc(k,iv) = vmc(k,iv) / (.5 * (rho(k,iw1) + rho(k,iw2)))

enddo

vc(1:kb-1,iv) = vc(kb,iv)

return
end subroutine prog_v_begs

!===============================================================================

subroutine vel_t3d(mrl,vs,ws,vxe,vye,vze)

use mem_ijtabs, only: jtab_w, itab_v, itab_w
use mem_grid,   only: mza, mva, mwa, lpw, vnx, vny, vnz, wnx, wny, wnz
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in) :: vs(mza,mva)
real, intent(in) :: ws(mza,mwa)

real, intent(out) :: vxe(mza,mwa)
real, intent(out) :: vye(mza,mwa)
real, intent(out) :: vze(mza,mwa)

integer :: j,iw,npoly,kb,k,jv,iv
real :: farv2

real :: wst(mza)

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,kb,k,jv,iv,wst,farv2) 
do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

! Vertical loop over T levels

   do k = kb,mza-1

! Diagnose 3D earth-velocity vector at T points; W contribution first

      wst(k) = .5 * (ws(k-1,iw) + ws(k,iw))

      vxe(k,iw) = wst(k) * wnx(iw)
      vye(k,iw) = wst(k) * wny(iw)
      vze(k,iw) = wst(k) * wnz(iw)

   enddo

! Loop over V neighbors of this W cell

   do jv = 1,npoly
      iv  = itab_w(iw)%iv(jv)

      farv2 = 2. * itab_w(iw)%farv(jv)

! Vertical loop over T levels

      do k = kb,mza-1

! Diagnose 3D earth-velocity vector at T points; VC contribution

         vxe(k,iw) = vxe(k,iw) + farv2 * vs(k,iv) * vnx(iv)
         vye(k,iw) = vye(k,iw) + farv2 * vs(k,iv) * vny(iv)
         vze(k,iw) = vze(k,iw) + farv2 * vs(k,iv) * vnz(iv)

      enddo

   enddo

enddo
!$omp end parallel do
call rsub('Wa',16)

return
end subroutine vel_t3d
