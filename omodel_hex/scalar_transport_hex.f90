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
subroutine scalar_transport(vmsc,wmsc,rho_old)

use mem_ijtabs, only: istp, jtab_v, jtab_w, mrl_endl, itab_v, itab_w
use mem_grid,   only: mza, mva, mwa, lpv, lpw, lsw, zt, zm, dzim, &
                      dniv, volt, arv, arw, dzim, volti
use misc_coms,  only: io6, dtlm
use var_tables, only: num_scalar, scalar_tab
use mem_turb,   only: vkh, hkm, sxfer_rk
use massflux,   only: tridiffo

!$ use omp_lib

implicit none

real, intent(in) :: vmsc(mza,mva)
real, intent(in) :: wmsc(mza,mwa)
real(kind=8), intent(in) :: rho_old(mza,mwa)

integer :: j,iw,iw1,iw2,iwd,mrl
integer :: n,k,kb,kbv,ks,kd,iv,iwn,jv,npoly
real :: dtl,dtli
real :: dirv,hdniv

! Automatic arrays:

real :: vsc(mza,mva)
real :: wsc(mza,mwa)

real :: vxe(mza,mwa) ! XE velocity component at T point
real :: vye(mza,mwa) ! YE velocity component at T point
real :: vze(mza,mwa) ! ZE velocity component at T point

real :: gxps_scp(mza,mwa)
real :: gyps_scp(mza,mwa)
real :: gzps_scp(mza,mwa)

integer :: iwdepv(mza,mva)
integer :: kdepw(mza,mwa)

real :: dxps_v(mza,mva)
real :: dyps_v(mza,mva)
real :: dzps_v(mza,mva)

real :: dxps_w(mza,mwa)
real :: dyps_w(mza,mwa)
real :: dzps_w(mza,mwa)

real :: akodz(mza),dtomass(mza)
real :: vctr5(mza),vctr6(mza),vctr7(mza),vctr8(mza)
real :: hfluxadv(mza),hfluxdif(mza),vfluxadv(mza),vfluxdif(mza)

real :: del_scp(mza)

real, pointer :: scp(:,:)
real, pointer :: sct(:,:)

! Horizontal loop over W columns for ENDL

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,kb,k) 
do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   kb = lpw(iw)

! Diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

! Vertical loop over W levels

   do k = kb,mza-2
   
!write(6,'(a,2i7,3e15.3)') 'sct1 ',k,iw,wmsc(k,iw),rho_old(k,iw),rho_old(k+1,iw)
   
      wsc(k,iw) = wmsc(k,iw) / (.5 * (rho_old(k,iw) + rho_old(k+1,iw)))
   enddo
   
   wsc(kb-1,iw) = wsc(kb,iw)
   wsc(mza-1,iw) = 0.

enddo
!$omp end parallel do
endif
call rsub('Wa',16)

! Horizontal loop over V columns for ENDL

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private(iv,iw1,iw2,kb,k) 
do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)
iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
call qsub('V',iv)

   kb = lpv(iv)

! Diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

! Vertical loop over T levels

   do k = kb,mza-1
      vsc(k,iv) = vmsc(k,iv) / (.5 * (rho_old(k,iw1) + rho_old(k,iw2)))
   enddo

enddo
endif
call rsub('V',12)

! Diagnose 3D velocity at T points for ENDL using velocities for scalar advection

mrl = mrl_endl(istp)
if (mrl > 0) then
   call vel_t3d(mrl,vsc,wsc,vxe,vye,vze)
endif

!---------------------------------------------------------
! Parallel send/recv of VXE, VYE, VZE
!---------------------------------------------------------

! Diagnose advective donor point location for all V and W faces

mrl = mrl_endl(istp)
if (mrl > 0) then
   call donorpoint3d(1,mrl,vsc,wsc,vxe,vye,vze, &
      iwdepv,kdepw,dxps_v,dyps_v,dzps_v,dxps_w,dyps_w,dzps_w)
endif

! LOOP OVER SCALARS HERE

do n = 1,num_scalar

! Point SCP and SCT to scalar table arrays

   scp => scalar_tab(n)%var_p
   sct => scalar_tab(n)%var_t

! Evaluate T3D gradient of scalar field

   mrl = mrl_endl(istp)
   if (mrl > 0) then
      call grad_t3d(mrl,scp,gxps_scp,gyps_scp,gzps_scp)
   endif

!---------------------------------------------------------
! Parallel send/recv of T3D gradients (3 of them)
!---------------------------------------------------------

! Horizontal loop over W/T points

   call psub()
!----------------------------------------------------------------------
   mrl = mrl_endl(istp)
   if (mrl > 0) then
   !$omp parallel do private (iw,kb,k,npoly,jv,iv,iwn,kbv,iwd,kd,ks, &
   !$omp                      dtl,dtli,dtomass,hfluxadv,hfluxdif,dirv, &
   !$omp                      vfluxadv,akodz,vctr5,vctr7,vctr6,vctr8, &
   !$omp                      hdniv,del_scp,vfluxdif)
   do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)
!----------------------------------------------------------------------
   call qsub('W',iw)

      kb = lpw(iw)

      dtl = dtlm(itab_w(iw)%mrlw)
      dtli = 1. / dtl
   
! Vertical loop over T levels 

      do k = kb,mza-1
         dtomass(k) = dtl / (rho_old(k,iw)  * volt(k,iw))
      enddo

! Loop over neighbor V points of this W cell

      npoly = itab_w(iw)%npoly

      hfluxadv(1:mza) = 0.
      hfluxdif(1:mza) = 0.

      do jv = 1,npoly
         iv  = itab_w(iw)%iv(jv)
         iwn = itab_w(iw)%iw(jv)
         
         kbv = lpv(iv)

         dirv = itab_w(iw)%dirv(jv)

         hdniv = .5 * dniv(iv)
   
! Vertical loop over T levels 

         do k = kbv,mza-1

! Horizontal advective donor cell

            iwd = iwdepv(k,iv)

! Horizontal advective and diffusive scalar fluxes

            hfluxadv(k) = hfluxadv(k) &
                        + dirv * vmsc(k,iv) * arv(k,iv) * (scp(k,iwd) &
                        + dxps_v(k,iv) * gxps_scp(k,iwd) &
                        + dyps_v(k,iv) * gyps_scp(k,iwd) &
                        + dzps_v(k,iv) * gzps_scp(k,iwd))

            hfluxdif(k) = hfluxdif(k) &
                        + hdniv * arv(k,iv) * (hkm(k,iwn) + hkm(k,iw)) &
                        * (scp(k,iwn) - scp(k,iw))

         enddo
      enddo

! Vertical loop over W levels

      do k = kb,mza-2

! Vertical advective donor cell k level

         kd = kdepw(k,iw)

! vertical scalar advective flux

         vfluxadv(k) = wmsc(k,iw) * arw(k,iw) &
                     * (scp(kd,iw) &
                      + dxps_w(k,iw) * gxps_scp(kd,iw) &
                      + dyps_w(k,iw) * gyps_scp(kd,iw) &
                      + dzps_w(k,iw) * gzps_scp(kd,iw))

! Prepare for vertical diffusion - Fill tri-diagonal matrix coefficients

         akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
         vctr5(k) = - akodz(k) * dtomass(k)
         vctr7(k) = - akodz(k) * dtomass(k+1)
         vctr6(k) = 1. - vctr5(k) - vctr7(k)
         vctr8(k) = akodz(k) * (scp(k,iw) - scp(k+1,iw))
      enddo

! Special case for total water SH_W: Apply surface vapor flux

      if (scalar_tab(n)%name == 'SH_W') then

! Vertical loop over T levels that are adjacent to surface

         do ks = 1,lsw(iw)
            k = kb + ks - 1

! Apply surface vapor xfer [kg_vap] directly to SCT [kg_vap / (m^3 s)]

            sct(k,iw) = sct(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

! Change in SCP from surface xfer

            del_scp(k) = sxfer_rk(ks,iw) / (rho_old(k,iw) * volt(k,iw))

! Zero out sxfer_rk(ks,iw) now that it has been transferred to the atm

            sxfer_rk(ks,iw) = 0.  

         enddo

! Lowest T level that is not adjacent to surface

         del_scp(kb+lsw(iw)) = 0.

! Vertical loop over W levels that are adjacent to surface

         do ks = 1,lsw(iw)
            k = kb + ks - 1

! Change in vctr8 from surface vapor xfer

            vctr8(k) = vctr8(k) + akodz(k) * (del_scp(k) - del_scp(k+1))
         enddo

      endif

! Solve tri-diagonal matrix equation

      if (kb < mza-2) then
         call tridiffo(mza,kb,mza-2,vctr5,vctr6,vctr7,vctr8,vfluxdif)
      endif

! Set bottom and top vertical internal turbulent fluxes to zero

      vfluxdif(kb-1)  = 0.
      vfluxdif(mza-1) = 0.

! Set bottom & top vertical advective fluxes to zero

      vfluxadv(kb-1)  = 0.
      vfluxadv(mza-1) = 0.

! Vertical loop over T levels

      do k = kb,mza-1

! Add contributions to scalar tendency from horizontal and vertical 
! advection and diffusion

         sct(k,iw) = sct(k,iw) + volti(k,iw) &
            * (hfluxadv(k) + vfluxadv(k-1) - vfluxadv(k) &
             + hfluxdif(k) + vfluxdif(k-1) - vfluxdif(k))

      enddo

   enddo
   !$omp end parallel do
   endif
   call rsub('W',26)

enddo ! n
return
end subroutine scalar_transport

!===============================================================================

subroutine donorpoint3d(ldt,mrl,vs,ws,vxe,vye,vze, &
   iwdepv,kdepw,dxps_v,dyps_v,dzps_v,dxps_w,dyps_w,dzps_w)

use mem_ijtabs, only: jtab_v, jtab_w, itab_v, itab_w
use mem_grid,   only: mza, mva, mwa, lpv, lpw, &
                      zt, zm, unx, uny, unz, xev, yev, zev
use misc_coms,  only: io6, dtlm, dtsm
use max_dims,   only: maxgrds
use consts_coms, only: eradi

!$ use omp_lib

implicit none

integer, intent(in) :: ldt,mrl

real, intent(in) :: vs(mza,mva)
real, intent(in) :: ws(mza,mwa)

real, intent(in) :: vxe(mza,mwa)
real, intent(in) :: vye(mza,mwa)
real, intent(in) :: vze(mza,mwa)

integer, intent(out) :: iwdepv(mza,mva)
integer, intent(out) :: kdepw(mza,mwa)

real, intent(out) :: dxps_v(mza,mva)
real, intent(out) :: dyps_v(mza,mva)
real, intent(out) :: dzps_v(mza,mva)

real, intent(out) :: dxps_w(mza,mwa)
real, intent(out) :: dyps_w(mza,mwa)
real, intent(out) :: dzps_w(mza,mwa)

integer :: j,kb,k,iv,iw,iw1,iw2

real :: dto2
real :: dxps1,dyps1,dxps2,dyps2,cosv1,sinv1,cosv2,sinv2
real :: unx_w, uny_w, vnx_w, vny_w, vnz_w, wnx_v, wny_v, wnz_v
real :: vxeface,vyeface,vzeface,uface,vface,wface

real :: dtm(maxgrds)

if (ldt == 1) then
   dtm(:) = dtlm(:)
else
   dtm(:) = dtsm(:)
endif

! Horizontal loop over V points

call psub()
!----------------------------------------------------------------------------
!$omp parallel do private(iv,iw1,iw2,kb,k,dto2,dxps1,dyps1,dxps2,dyps2, &
!$omp                     cosv1,sinv1,cosv2,sinv2,wnx_v,wny_v,wnz_v, &
!$omp                     vxeface,vyeface,vzeface,uface,wface) 
do j = 1,jtab_v(15)%jend(mrl); iv = jtab_v(15)%iv(j)
iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------------
call qsub('V',iv)

   kb = lpv(iv)

   dto2 = .5 * dtm(itab_v(iv)%mrlv)

   dxps1 = itab_v(iv)%dxps(1)
   dyps1 = itab_v(iv)%dyps(1)

   dxps2 = itab_v(iv)%dxps(2)
   dyps2 = itab_v(iv)%dyps(2)

   cosv1 = itab_v(iv)%cosv(1)
   sinv1 = itab_v(iv)%sinv(1)

   cosv2 = itab_v(iv)%cosv(2)
   sinv2 = itab_v(iv)%sinv(2)
   
   wnx_v = xev(iv) * eradi
   wny_v = yev(iv) * eradi
   wnz_v = zev(iv) * eradi

! Vertical loop over T/V levels

   do k = kb,mza-1

! Average 3 earth velocity components from T points to V face

      vxeface = .5 * (vxe(k,iw1) + vxe(k,iw2))
      vyeface = .5 * (vye(k,iw1) + vye(k,iw2))
      vzeface = .5 * (vze(k,iw1) + vze(k,iw2))

! Project earth velocity components at V face onto U and W directions

      uface = unx(iv) * vxeface + uny(iv) * vyeface + unz(iv) * vzeface
      wface = wnx_v   * vxeface + wny_v   * vyeface + wnz_v   * vzeface

! Compute displacement components for V face relative to T point

      if (vs(k,iv) > 0.) then

         iwdepv(k,iv) = iw1

         dxps_v(k,iv) = -dto2 * (vs(k,iv) * cosv1 - uface * sinv1) + dxps1
         dyps_v(k,iv) = -dto2 * (vs(k,iv) * sinv1 + uface * cosv1) + dyps1
         dzps_v(k,iv) = -dto2 * wface

      else

         iwdepv(k,iv) = iw2

         dxps_v(k,iv) = -dto2 * (vs(k,iv) * cosv2 - uface * sinv2) + dxps2
         dyps_v(k,iv) = -dto2 * (vs(k,iv) * sinv2 + uface * cosv2) + dyps2
         dzps_v(k,iv) = -dto2 * wface

      endif         

   enddo

enddo
call rsub('V',15)

! Horizontal loop over W points

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,kb,k,dto2,unx_w,uny_w,vnx_w,vny_w,vnz_w, &
!$omp                     vxeface,vyeface,vzeface,uface,vface) 
do j = 1,jtab_w(15)%jend(mrl); iw = jtab_w(15)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   kb = lpw(iw)

   dto2 = dtm(itab_w(iw)%mrlw)

   unx_w = itab_w(iw)%unx_w
   uny_w = itab_w(iw)%uny_w

   vnx_w = itab_w(iw)%vnx_w
   vny_w = itab_w(iw)%vny_w
   vnz_w = itab_w(iw)%vnz_w

! Vertical loop over W levels

   do k = kb,mza-1

! Average 3 velocity components from T points to W point

      vxeface = .5 * (vxe(k,iw) + vxe(k+1,iw))
      vyeface = .5 * (vye(k,iw) + vye(k+1,iw))
      vzeface = .5 * (vze(k,iw) + vze(k+1,iw))

! Project earth velocity components at W face onto U and V directions

      uface = unx_w * vxeface + uny_w * vyeface
      vface = vnx_w * vxeface + vny_w * vyeface + vnz_w * vzeface

! Compute displacement components for W face relative to T point

      if (ws(k,iw) > 0.) then

         kdepw(k,iw) = k

         dxps_w(k,iw) = -dto2 * uface
         dyps_w(k,iw) = -dto2 * vface
         dzps_w(k,iw) = -dto2 * ws(k,iw) + zm(k) - zt(k)

      else

         kdepw(k,iw) = k + 1

         dxps_w(k,iw) = -dto2 * uface
         dyps_w(k,iw) = -dto2 * vface
         dzps_w(k,iw) = -dto2 * ws(k,iw) + zm(k) - zt(k+1)

      endif         

   enddo
   
enddo
call rsub('V',15)

return
end subroutine donorpoint3d

!=========================================================================

subroutine grad_t3d(mrl,scp,gxps,gyps,gzps)

use mem_ijtabs, only: jtab_w, itab_w
use mem_grid,   only: mza, mwa, lpw, zt, zm, dzim, vnx, vny, vnz, wnx, wny, wnz
use misc_coms,  only: io6
use max_dims,   only: maxgrds

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in) :: scp(mza,mwa)

real, intent(out) :: gxps(mza,mwa)
real, intent(out) :: gyps(mza,mwa)
real, intent(out) :: gzps(mza,mwa)

integer :: j,iw,npoly,kb,jw1,jw2,iw1,iw2,kb1,kb2,k

real :: farm,wti
real :: gxps1,gyps1,gxps2,gyps2

real :: wt(mza),gwz(mza),gz(mza)

! Horizontal loop over W columns for BEGS

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,kb,jw1,jw2,iw1,iw2,kb1,kb2,k, &
!$omp                     wt,farm,gxps1,gyps1,gxps2,gyps2,gwz,wti)
do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)
   
   gxps(:,iw) = 0.
   gyps(:,iw) = 0.
   
   wt(:) = 0.

! Loop over W neighbors of this W cell

   do jw1 = 1,npoly
      jw2 = jw1 + 1
      if (jw1 == npoly) jw2 = 1

      iw1 = itab_w(iw)%iw(jw1)
      iw2 = itab_w(iw)%iw(jw2)
      
      kb1 = lpw(iw1)
      kb2 = lpw(iw2)
      
      farm = itab_w(iw)%farm(jw1)

      gxps1 = itab_w(iw)%gxps1(jw1)
      gyps1 = itab_w(iw)%gyps1(jw1)

      gxps2 = itab_w(iw)%gxps2(jw1)
      gyps2 = itab_w(iw)%gyps2(jw1)

! Vertical loop over T levels

      do k = kb,mza-1
      
! Gradient calculation when all 3 W points are above ground

         if (k >= kb1 .and. k >= kb2) then

            gxps(k,iw) = gxps(k,iw) + farm * (gxps1 * (scp(k,iw1) - scp(k,iw)) &
                                            + gxps2 * (scp(k,iw2) - scp(k,iw)))

            gyps(k,iw) = gyps(k,iw) + farm * (gyps1 * (scp(k,iw1) - scp(k,iw)) &
                                            + gyps2 * (scp(k,iw2) - scp(k,iw)))

            wt(k) = wt(k) + farm
            
         elseif (k >= kb1) then

! Gradient calculation when IW2 point is below ground

         elseif (k >= kb2) then

! Gradient calculation when IW1 point is below ground

         endif
         
      enddo
      
   enddo
   
! Vertical loop over W levels

   do k = kb,mza-2
      gwz(k) = dzim(k) * (scp(k+1,iw) - scp(k,iw))
   enddo
   
   gwz(kb-1) = gwz(kb)
   gwz(mza-1) = gwz(mza-2)

! Vertical loop over T levels

   do k = kb,mza-1

! Diagnose 3D gradient at T points

      wti = 1. / wt(k)

      gxps(k,iw) = gxps(k,iw) * wti
      gyps(k,iw) = gyps(k,iw) * wti
      gzps(k,iw) = .5 * (gwz(k-1) + gwz(k))

   enddo

enddo
!$omp end parallel do
call rsub('Wa',16)

return
end subroutine grad_t3d

