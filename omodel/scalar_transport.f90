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

subroutine scalar_transport(mrl,vmsc, wmsc, vxesc, vyesc, vzesc, rho_old)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, &
                          jtv_wadj, jtw_prog, jtw_wadj
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, arv, arw, volti, &
                          zwgt_bot8, zwgt_top8
  use misc_coms,    only: dtlm, iparallel
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: akhodx
  use oname_coms,   only: nl
  use mem_thuburn,  only: comp_cfls_long, comp_and_apply_monot_limits, &
                          comp_and_apply_pd_limits
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use consts_coms,  only: r8
  use grad_lib,     only: grad_t2d, grad_t2d_quad, grad_z_quad
  use mem_adv,      only: dxps_w, dyps_w, dzps_w, dzzps_w, &
                          dxps_v, dyps_v, dzps_v, dxyps_v, dxxps_v, dyyps_v, &
                          gxps_scp, gyps_scp, gzps_scp, &
                          gxyps_scp, gxxps_scp, gyyps_scp, gzzps_scp
  implicit none

  integer, intent(in) :: mrl

  real, intent(in) :: vmsc(mza,mva)
  real, intent(in) :: wmsc(mza,mwa)

  real, intent(inout) :: vxesc(mza,mwa)
  real, intent(inout) :: vyesc(mza,mwa)
  real, intent(inout) :: vzesc(mza,mwa)

  real(r8), intent(in) :: rho_old(mza,mwa)

  integer  :: j,iw,iw1,iw2,iwd
  integer  :: n,k,kb,kd,iv,iwn,jv
  real     :: dirv, dt

! Automatic arrays:

  real :: vsc(mza)
  real :: wsc(mza)

  real :: vmsca(mza,mva)
  real :: wmsca(mza,mwa)

  integer :: iwdepv(mza,mva)

  integer :: kdepw(mza,mwa)

  real :: scp_upv(mza,mva)
  real :: scp_upw(mza,mwa)

  real :: hflux(mza), vfluxadv(mza)

  real, pointer, contiguous :: scp(:,:)
  real, pointer, contiguous :: sct(:,:)

! Return if this is not the end of the long timestep on any MRL

  if (mrl == 0) return

! Horizontal loop over all primary W columns to diagnose
! face-normal vertical velocity at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp parallel private(wsc,vsc)
  !$omp do private(iw,kb,dt,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     kb = lpw(iw)
     dt = dtlm(itab_w(iw)%mrlw)

     do k = kb, mza-1
        wsc  (k)    = wmsc(k,iw) &
                    / real(zwgt_bot8(k) * rho_old(k,iw) + zwgt_top8(k) * rho_old(k+1,iw))
        wmsca(k,iw) = wmsc(k,iw) * arw(k,iw)
     enddo

     wsc  (kb-1)    = 0.0
     wmsca(kb-1,iw) = 0.0

     wsc  (mza)    = 0.0
     wmsca(mza,iw) = 0.0

     call donorpointw(iw, dt, wsc, vxesc(:,iw), vyesc(:,iw), vzesc(:,iw), kdepw(:,iw))
  enddo
  !$omp end do nowait

! Loop over V columns (the immediate neighbors of all primary W points) to
! diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp do private(iv,iw1,iw2,kb,dt,k)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)
     kb  = lpv(iv)
     dt  = dtlm(itab_v(iv)%mrlv)

     do k = kb, mza
        vsc  (k)    = 2.0 * vmsc(k,iv) / real( rho_old(k,iw1) + rho_old(k,iw2) )
        vmsca(k,iv) = vmsc(k,iv) * arv(k,iv)
     enddo

     call donorpointv(iv, dt, vsc, vxesc, vyesc, vzesc, iwdepv(:,iv))

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Compute outflow CFL numbers for long timestep stability check;
  ! also needed by Thuburn monotonic scheme

  call comp_cfls_long( mrl, vmsca, wmsca, rho_old, nl%iscal_monot, do_check=.true. )

! LOOP OVER SCALARS HERE
  do n = 1, num_scalar

! Point SCP and SCT to scalar table arrays

     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

! Evaluate horizontal gradients of scalar field on a local polar
! stereographic projection

     !$omp parallel do private(iw)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        if (nl%horiz_adv_order <= 2) then
           call grad_t2d  (iw, scp, gxps_scp(:,iw), gyps_scp(:,iw))
        else
           call grad_t2d_quad(iw, scp, gxps_scp (:,iw), gyps_scp (:,iw), &
                                       gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
        endif

     enddo
     !$omp end parallel do

! MPI send of SCP horizontal gradient components

     if (iparallel == 1) then
        if (nl%horiz_adv_order <= 2) then
           call mpi_send_w(mrl, rvara1=gxps_scp,  rvara2=gyps_scp)
        else
           call mpi_send_w(mrl, rvara1=gxps_scp,  rvara2=gyps_scp,  &
                                rvara3=gxxps_scp, rvara4=gxyps_scp, rvara5=gyyps_scp)
        endif
     endif

! Evaluate vertical gradient of scalar field. This can overlap with
! communication of horizontal gradients

     !$omp parallel
     !$omp do private(iw)
     do j = 1,jtab_w(jtw_wadj)%jend(mrl); iw = jtab_w(jtw_wadj)%iw(j)

        call grad_z_quad(iw, scp(:,iw), gzps_scp(:,iw), gzzps_scp(:,iw))

     enddo
     !$omp end do

! Horizontal loop over all primary W columns to compute the
! upwinded scalar value at each W interface. This can be
! computed before scalar gradients are received at the borders

     !$omp do private(iw,k,kd)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        do k = lpw(iw), mza-1
           kd = kdepw(k,iw)

           scp_upw(k,iw) =                       scp(kd,iw) &
                         +  dxps_w(k,iw) *  gxps_scp(kd,iw) &
                         +  dyps_w(k,iw) *  gyps_scp(kd,iw) &
                         +  dzps_w(k,iw) *  gzps_scp(kd,iw) &
                         + dzzps_w(k,iw) * gzzps_scp(kd,iw)
        enddo
     enddo
     !$omp end do nowait
     !$omp end parallel


! MPI recv of SCP horizontal gradient components

     if (iparallel == 1) then
        if (nl%horiz_adv_order <= 2) then
           call mpi_recv_w(mrl, rvara1=gxps_scp,  rvara2=gyps_scp)
        else
           call mpi_recv_w(mrl, rvara1=gxps_scp,  rvara2=gyps_scp,  &
                                rvara3=gxxps_scp, rvara4=gxyps_scp, rvara5=gyyps_scp)
        endif
     endif

     if (nl%horiz_adv_order <= 2) then
        call lbcopy_w(mrl, a1=gxps_scp,  a2=gyps_scp)
     else
        call lbcopy_w(mrl, a1=gxps_scp,  a2=gyps_scp,  &
                           a3=gxxps_scp, a4=gxyps_scp, a5=gyyps_scp)
     endif

! Horizontal loop over all V points surrounding primary W/T columns
! to compute upwinded scalar value at the V points

     !$omp parallel do private(iv,k,iwd)
     do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

        if (nl%horiz_adv_order <= 2) then

           do k = lpv(iv), mza
              iwd = iwdepv(k,iv)
              scp_upv(k,iv) =                     scp(k,iwd) &
                            + dxps_v(k,iv) * gxps_scp(k,iwd) &
                            + dyps_v(k,iv) * gyps_scp(k,iwd) &
                            + dzps_v(k,iv) * gzps_scp(k,iwd)
           enddo

        else

           do k = lpv(iv), mza
              iwd = iwdepv(k,iv)
              scp_upv(k,iv) =                       scp(k,iwd) &
                            +  dxps_v(k,iv) *  gxps_scp(k,iwd) &
                            + dxxps_v(k,iv) * gxxps_scp(k,iwd) &
                            +  dyps_v(k,iv) *  gyps_scp(k,iwd) &
                            + dyyps_v(k,iv) * gyyps_scp(k,iwd) &
                            + dxyps_v(k,iv) * gxyps_scp(k,iwd) &
                            +  dzps_v(k,iv) *  gzps_scp(k,iwd)
           enddo

        endif

     enddo
     !$omp end parallel do

! Compute and apply the monotonic or positive-definite flux limiters.

     if (nl%iscal_monot == 1) then
        call comp_and_apply_monot_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv, &
             wmsca, vmsca)
     elseif (nl%iscal_monot == 2 .and. scalar_tab(n)%pdef) then
        call comp_and_apply_pd_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv, &
             wmsca, vmsca)
     endif

! Horizontal loop over W/T points

     !$omp parallel private(hflux,vfluxadv)
     !$omp do private (iw,kb,jv,iv,iwn,dirv,k)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        kb = lpw(iw)

! Loop over neighbor V points of this W cell

        hflux(kb:mza) = 0.

        do jv = 1, itab_w(iw)%npoly

           iv   = itab_w(iw)%iv(jv)
           iwn  = itab_w(iw)%iw(jv)
           dirv = itab_w(iw)%dirv(jv)

           if (nl%split_scalars < 1) then

              ! Horizontal advective and diffusive scalar fluxes

              ! Vertical loop over T levels
              do k = lpv(iv), mza
                 hflux(k) = hflux(k) + dirv * vmsca(k,iv) * scp_upv(k,iv)    &
                                     + akhodx(k,iv) * (scp(k,iwn) - scp(k,iw))
              enddo

           else

              ! Horizontal advective fluxes

              ! Vertical loop over T levels
              do k = lpv(iv), mza
                 hflux(k) = hflux(k) + dirv * vmsca(k,iv) * scp_upv(k,iv)
              enddo

           endif

        enddo

! Vertical loop over W levels

        do k = kb, mza-1
           vfluxadv(k) = wmsca(k,iw) * scp_upw(k,iw)
        enddo

! Set bottom & top vertical advective fluxes to zero

        vfluxadv(kb-1) = 0.
        vfluxadv(mza)  = 0.

! Vertical loop over T levels

        do k = kb,mza

! Add contributions to scalar tendency from horizontal and vertical
! advection and diffusion

           sct(k,iw) = sct(k,iw) &
                     + volti(k,iw) * (hflux(k) + vfluxadv(k-1) - vfluxadv(k))

        enddo

     enddo
     !$omp end do nowait
     !$omp end parallel

  enddo ! n

end subroutine scalar_transport


subroutine donorpointv(iv, dt, vs, vxe, vye, vze, iwdepv)

  use mem_ijtabs,  only: itab_v, jtv_wadj
  use mem_grid,    only: mza, mwa, lpv, unx, uny, unz, xev, yev, zev, dzto4
  use misc_coms,   only: mdomain
  use consts_coms, only: eradi
  use oname_coms,  only: nl
  use mem_adv,     only: dxps_v, dyps_v, dzps_v, dxyps_v, dxxps_v, dyyps_v, &
                         xx0_v, xy0_v, yy0_v, dxps_m1, dyps_m1, dxps_m2, dyps_m2

  implicit none

  integer, intent(in)  :: iv
  real,    intent(in)  :: dt
  real,    intent(in)  :: vs (mza)
  real,    intent(in)  :: vxe(mza,mwa)
  real,    intent(in)  :: vye(mza,mwa)
  real,    intent(in)  :: vze(mza,mwa)
  integer, intent(out) :: iwdepv(mza)

  real :: dto2
  real :: dxps1, dyps1, dxps2, dyps2
  real :: cosv1, sinv1, cosv2, sinv2
  real :: wnx_v, wny_v, wnz_v
  real :: vxeface, vyeface, vzeface

  real :: dx, dxo2, dy, dyo2
  real :: dx1, dy1, dx2, dy2

  real :: ufacev(mza), wfacev

  real, parameter :: wt1 = 8.0 / 12.
  real, parameter :: wt2 = 1.0 / 12.

  integer :: kb, k, iw1, iw2

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  kb = lpv(iv)

  dto2 = .5 * dt

  dxps1 = itab_v(iv)%dxps(1)
  dyps1 = itab_v(iv)%dyps(1)

  dxps2 = itab_v(iv)%dxps(2)
  dyps2 = itab_v(iv)%dyps(2)

  cosv1 = itab_v(iv)%cosv(1)
  sinv1 = itab_v(iv)%sinv(1)

  cosv2 = itab_v(iv)%cosv(2)
  sinv2 = itab_v(iv)%sinv(2)

  if (mdomain <= 1) then
     wnx_v = xev(iv) * eradi
     wny_v = yev(iv) * eradi
     wnz_v = zev(iv) * eradi
  else
     wnx_v = 0.
     wny_v = 0.
     wnz_v = 1.
  endif

! Vertical loop over T/V levels

  do k = kb, mza

! Average 3 earth velocity components from T points to V face

     vxeface = .5 * (vxe(k,iw1) + vxe(k,iw2))
     vyeface = .5 * (vye(k,iw1) + vye(k,iw2))
     vzeface = .5 * (vze(k,iw1) + vze(k,iw2))

! Project earth velocity components at V face onto U and W directions

     ufacev(k) = unx(iv) * vxeface + uny(iv) * vyeface + unz(iv) * vzeface
     wfacev    = wnx_v   * vxeface + wny_v   * vyeface + wnz_v   * vzeface

! Compute z displacement component for V face relative to T point

     dzps_v(k,iv)  = min( max( -dto2 * wfacev, -dzto4(k) ), dzto4(k) )
  enddo

! Compute X, Y displacement components for V face relative to T point
! Vertical loop over T/V levels

  if (nl%horiz_adv_order <= 2) then

     do k = kb, mza

        if (vs(k) > 0.0) then

           iwdepv(k) = iw1
           dxps_v(k,iv) = -dto2 * (vs(k) * cosv1 - ufacev(k) * sinv1) + dxps1
           dyps_v(k,iv) = -dto2 * (vs(k) * sinv1 + ufacev(k) * cosv1) + dyps1

        else

           iwdepv(k) = iw2
           dxps_v(k,iv) = -dto2 * (vs(k) * cosv2 - ufacev(k) * sinv2) + dxps2
           dyps_v(k,iv) = -dto2 * (vs(k) * sinv2 + ufacev(k) * cosv2) + dyps2

        endif

     enddo

  else

     do k = kb, mza

        if (vs(k) > 0.0) then

           iwdepv(k) = iw1

           dxo2 = -dto2 * (vs(k) * cosv1 - ufacev(k) * sinv1)
           dx   = 2.0 * dxo2

           dyo2 = -dto2 * (vs(k) * sinv1 + ufacev(k) * cosv1)
           dy   = 2.0 * dyo2

           dxps_v(k,iv) = dxo2 + dxps1
           dyps_v(k,iv) = dyo2 + dyps1

           dx1 = dxps_m1(1,iv) + dx
           dy1 = dyps_m1(1,iv) + dy

           dx2 = dxps_m2(1,iv) + dx
           dy2 = dyps_m2(1,iv) + dy

           dxxps_v(k,iv) = xx0_v(1,iv) + wt1 * dxps_v(k,iv) * dxps_v(k,iv) &
                                       + wt2 * (dx1 * dx1 + dx2 * dx2)

           dxyps_v(k,iv) = xy0_v(1,iv) + wt1 * dxps_v(k,iv) * dyps_v(k,iv) &
                                       + wt2 * (dx1 * dy1 + dx2 * dy2)

           dyyps_v(k,iv) = yy0_v(1,iv) + wt1 * dyps_v(k,iv) * dyps_v(k,iv) &
                                       + wt2 * (dy1 * dy1 + dy2 * dy2)

        else

           iwdepv(k) = iw2

           dxo2 = -dto2 * (vs(k) * cosv2 - ufacev(k) * sinv2)
           dx   = 2.0 * dxo2

           dyo2 = -dto2 * (vs(k) * sinv2 + ufacev(k) * cosv2)
           dy   = 2.0 * dyo2

           dxps_v(k,iv) = dxo2 + dxps2
           dyps_v(k,iv) = dyo2 + dyps2

           dx1 = dxps_m1(2,iv) + dx
           dy1 = dyps_m1(2,iv) + dy

           dx2 = dxps_m2(2,iv) + dx
           dy2 = dyps_m2(2,iv) + dy

           dxxps_v(k,iv) = xx0_v(2,iv) + wt1 * dxps_v(k,iv) * dxps_v(k,iv) &
                                       + wt2 * (dx1 * dx1 + dx2 * dx2)

           dxyps_v(k,iv) = xy0_v(2,iv) + wt1 * dxps_v(k,iv) * dyps_v(k,iv) &
                                       + wt2 * (dx1 * dy1 + dx2 * dy2)

           dyyps_v(k,iv) = yy0_v(2,iv) + wt1 * dyps_v(k,iv) * dyps_v(k,iv) &
                                       + wt2 * (dy1 * dy1 + dy2 * dy2)

        endif

     enddo

  endif

end subroutine donorpointv


!===============================================================================


subroutine donorpointw(iw, dt, ws, vxe, vye, vze, kdepw)

  use mem_ijtabs, only: itab_w, jtw_prog
  use mem_grid,   only: mza, lpw, dzto2, dzt, dztsqo6
  use mem_adv,    only: dxps_w, dyps_w, dzps_w, dzzps_w

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: dt
  real,    intent(in)  :: ws (mza)
  real,    intent(in)  :: vxe(mza)
  real,    intent(in)  :: vye(mza)
  real,    intent(in)  :: vze(mza)

  integer, intent(out) :: kdepw(mza)

  real :: dto2, zp
  real :: unx_w, uny_w, vnx_w, vny_w, vnz_w
  real :: vxeface, vyeface, vzeface
  real :: uface,vface

  real, parameter :: onethird = 1.0 / 3.0

  integer :: kb, k

  kb = lpw(iw)

  dto2 = 0.5 * dt

  unx_w = itab_w(iw)%unx_w
  uny_w = itab_w(iw)%uny_w

  vnx_w = itab_w(iw)%vnx_w
  vny_w = itab_w(iw)%vny_w
  vnz_w = itab_w(iw)%vnz_w

! Vertical loop over W levels

  do k = kb, mza-1

! Average 3 velocity components from T points to W point

     vxeface = .5 * (vxe(k) + vxe(k+1))
     vyeface = .5 * (vye(k) + vye(k+1))
     vzeface = .5 * (vze(k) + vze(k+1))

! Project earth velocity components at W face onto U and V directions

     uface = unx_w * vxeface + uny_w * vyeface
     vface = vnx_w * vxeface + vny_w * vyeface + vnz_w * vzeface

! Compute x and y displacement components for W face relative to T point

     dxps_w(k,iw) = -dto2 * uface
     dyps_w(k,iw) = -dto2 * vface

  enddo

! Compute z displacement components for W face relative to T point,
! and store the vertical donor and recv cell indices.
! Vertical loop over W levels

  do k = kb, mza-1

     if (ws(k) > 0.0) then

        kdepw(k) = k
        dzps_w(k,iw) = max( -dto2 * ws(k) + dzto2(k), 0.0)

        zp = min( dt * abs(ws(k)), dzt(k) )
        dzzps_w(k,iw) = dztsqo6(k) + zp * (zp * onethird - dzto2(k))

     else

        kdepw(k) = k + 1
        dzps_w(k,iw) = min( -dto2 * ws(k) - dzto2(k+1), 0.0)

        zp = min( dt * abs(ws(k)), dzt(k+1) )
        dzzps_w(k,iw) = dztsqo6(k+1) + zp * (zp * onethird - dzto2(k+1))

     endif

  enddo

  kdepw(mza) = mza

end subroutine donorpointw


!===========================================================================


subroutine zero_momsc(vmsc,wmsc,vxesc,vyesc,vzesc,rho_old)

use mem_ijtabs, only: istp, mrl_begl, jtv_wstn, jtw_wstn
use mem_grid,   only: mza, mva, mwa, lpw
use mem_basic,  only: rho
use consts_coms,only: r8

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

real(r8), intent(out) :: rho_old(mza,mwa)

integer :: k,iv,iw,mrl

mrl = mrl_begl(istp)
if (mrl > 0) then

! Zero out long timestep mass flux components (used for scalar advective
! transport) so they may be summed over small timesteps

   !$omp parallel
   !$omp do
   do iv = 2, mva
      vmsc(:,iv) = 0.0
   enddo
   !$omp end do nowait

   !$omp do private(k)
   do iw = 2, mwa
      wmsc (:,iw) = 0.0

      vxesc(:,iw) = 0.0
      vyesc(:,iw) = 0.0
      vzesc(:,iw) = 0.0

      ! Save DT_LONG density for use with scalar updates
      do k = 1, lpw(iw)-2
         rho_old(k,iw) = 0.0
      enddo
      do k = lpw(iw)-1, mza
         rho_old(k,iw) = rho(k,iw)
      enddo

   enddo
   !$omp end do nowait
   !$omp end parallel

endif

end subroutine zero_momsc

!===========================================================================

subroutine timeavg_momsc(vmsc,wmsc,vxesc,vyesc,vzesc)

  use mem_ijtabs, only: itab_v, itab_w, istp, mrl_endl, jtv_prog, jtw_prog
  use mem_grid,   only: mza, mva, mwa, lpw
  use misc_coms,  only: nacoust

  implicit none

  real, intent(inout) :: vmsc(mza,mva)
  real, intent(inout) :: wmsc(mza,mwa)

  real, intent(inout) :: vxesc(mza,mwa)
  real, intent(inout) :: vyesc(mza,mwa)
  real, intent(inout) :: vzesc(mza,mwa)

  integer :: k,iv,iw,mrl,mrlv,mrlw
  real    :: acoi

  mrl = mrl_endl(istp)
  if (mrl > 0) then

     !$omp parallel
     !$omp do private(mrlv,acoi)
     do iv = 2, mva
        mrlv = itab_v(iv)%mrlv
        acoi = 1.0 / real(nacoust(mrlv))

        vmsc(:,iv) = vmsc(:,iv) * acoi
     enddo
     !$omp end do nowait

     !$omp do private(mrlw,acoi,k)
     do iw = 2, mwa
        mrlw = itab_w(iw)%mrlw
        acoi = 1.0 / real(nacoust(mrlw))

        wmsc(lpw(iw)-1,iw) = 0.0

        do k = lpw(iw), mza-1
           wmsc(k,iw) = wmsc(k,iw) * acoi
        enddo

        vxesc(:,iw) = vxesc(:,iw) * acoi
        vyesc(:,iw) = vyesc(:,iw) * acoi
        vzesc(:,iw) = vzesc(:,iw) * acoi

        wmsc(mza,iw) = 0.0
     enddo
     !$omp end do
     !$omp end parallel

endif

end subroutine timeavg_momsc

!===========================================================================

subroutine scalar_hdiff_split(mrl)

  use mem_grid,    only: mza, volti, lpw, lpv
  use mem_turb,    only: akhodx
  use mem_ijtabs,  only: jtab_w, jtw_prog, jtw_wadj, itab_w
  use var_tables,  only: num_scalar, scalar_tab
  use olam_mpi_atm,only: mpi_recv_w, mpi_send_w
  use consts_coms, only: r8

  implicit none

  integer,  intent(in) :: mrl

  integer              :: n, j, iw, k, iwn, iv, jv
  real                 :: hflux(mza)

  if (mrl == 0) return

  !$omp parallel private(hflux)
  !$omp do collapse(2) private(n,j,iw,jv,iv,iwn,k)
  do n = 1, num_scalar
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     hflux = 0.0

     do jv = 1, itab_w(iw)%npoly

        iv  = itab_w(iw)%iv(jv)
        iwn = itab_w(iw)%iw(jv)

        do k = lpv(iv), mza
           hflux(k) = hflux(k) + akhodx(k,iv) * ( scalar_tab(n)%var_p(k,iwn) &
                                                - scalar_tab(n)%var_p(k,iw ) )
        enddo

     enddo

     do k = lpw(iw), mza
        scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                  + volti(k,iw) * hflux(k)
     enddo

  enddo
  enddo
  !$omp end do
  !$omp end parallel

end subroutine scalar_hdiff_split
