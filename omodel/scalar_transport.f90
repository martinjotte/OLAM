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
subroutine scalar_transport(vmsc, wmsc, vxesc, vyesc, vzesc, rho_old)

  use mem_ijtabs,   only: istp, jtab_v, jtab_w, mrl_endl, itab_v, itab_w, &
                          jtv_wadj, jtw_prog
  use mem_grid,     only: mza, mva, mwa, nsw_max, lpv, lpw, lsw, zt, zm, dzim, &
                          volt, arv, arw, dzim, volti, zwgt_bot, zwgt_top
  use misc_coms,    only: io6, dtlm, iparallel, time8p
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: akhodx
  use mem_basic,    only: rho
  use oname_coms,   only: nl
  use mem_thuburn,  only: comp_cfl1, comp_cfl2, comp_vert_limits, &
                          comp_horiz_limits, apply_flux_limiters
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use consts_coms,  only: r8

  implicit none

  real, intent(in) :: vmsc(mza,mva)
  real, intent(in) :: wmsc(mza,mwa)

  real, intent(inout) :: vxesc(mza,mwa)
  real, intent(inout) :: vyesc(mza,mwa)
  real, intent(inout) :: vzesc(mza,mwa)

  real(r8), intent(in) :: rho_old(mza,mwa)

  integer  :: j,iw,iw1,iw2,iw3,iw4,iwd,iwr,mrl
  integer  :: n,k,kb,kbv,kd,kr,iv,iwn,jv,npoly
  real     :: dirv

! Automatic arrays:

  real :: vsc(mza,mva)
  real :: wsc(mza,mwa)

  real :: vmsca(mza,mva)
  real :: wmsca(mza,mwa)

  real :: gxps_scp(mza,mwa)
  real :: gyps_scp(mza,mwa)
  real :: gzps_scp(mza,mwa)

  integer :: iwdepv(mza,mva)
  integer :: iwrecv(mza,mva)

  integer :: kdepw(mza,mwa)
  integer :: krecw(mza,mwa)

  real :: dxps_v(mza,mva)
  real :: dyps_v(mza,mva)
  real :: dzps_v(mza,mva)

  real :: dxps_w(mza,mwa)
  real :: dyps_w(mza,mwa)
  real :: dzps_w(mza,mwa)

  real :: scp_upv(mza,mva)
  real :: scp_upw(mza,mwa)

  real :: hfluxadv(mza), hfluxdif(mza), vfluxadv(mza)

  real, pointer :: scp(:,:)
  real, pointer :: sct(:,:)

! Return if this is not the end of the long timestep on any MRL

  mrl = mrl_endl(istp)
  if (mrl == 0) return

! MPI send of VXESC, VYESC, VZESC

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxesc, rvara2=vyesc, rvara3=vzesc)
  endif

! Horizontal loop over all primary W columns to diagnose
! face-normal vertical velocity at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp parallel 
  !$omp do private(iw,kb,k) 
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     kb = lpw(iw)

     do k = kb, mza-1
        wsc  (k,iw) = wmsc(k,iw) &
                    / (zwgt_bot(k) * rho_old(k,iw) + zwgt_top(k) * rho_old(k+1,iw))
        wmsca(k,iw) = wmsc(k,iw) * arw(k,iw)
     enddo

     wsc  (kb-1,iw) = 0.0
     wmsca(kb-1,iw) = 0.0

     wsc  (mza,iw) = 0.0
     wmsca(mza,iw) = 0.0

  enddo
  !$omp end do

! Loop over V columns (the immediate neighbors of all primary W points) to
! diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp do private(iv,iw1,iw2,kb,k) 
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)
     kb  = lpv(iv)
     
     do k = kb, mza
        vsc  (k,iv) = 2.0 * vmsc(k,iv) / ( rho_old(k,iw1) + rho_old(k,iw2) )
        vmsca(k,iv) = vmsc(k,iv) * arv(k,iv)
     enddo

     vsc  (1:kb-1,iv) = 0.0
     vmsca(1:kb-1,iv) = 0.0

  enddo
  !$omp end do
  !$omp end parallel

! Diagnose advective donor point locations for all primary W faces
! No parallel communication is necessary to compute this

  call donorpointw(1, mrl, wsc, vxesc, vyesc, vzesc, kdepw, krecw, &
                   dxps_w, dyps_w, dzps_w)

! Complete MPI recv of VXESC, VYESC, VZESC and do LBC copy

  if (iparallel == 1) then
     call mpi_recv_w(mrl, rvara1=vxesc, rvara2=vyesc, rvara3=vzesc)
  endif

  call lbcopy_w(mrl, a1=vxesc, a2=vyesc, a3=vzesc)

! Diagnose advective donor point locations for the V faces surrounding all
! primary W points. Communication of velocities must have been completed

  call donorpointv(1, mrl, vsc, vxesc, vyesc, vzesc, iwdepv, iwrecv, &
                   dxps_v, dyps_v, dzps_v)

  ! Compute outflow CFL numbers for long timestep stability check;
  ! also needed by Thuburn monotonic scheme

  call comp_cfl1( mrl, dtlm, vmsca, wmsca, vsc, wsc, rho_old, &
                  iwdepv, kdepw, do_check=.true. )

  ! Diagnose inlfow CFL numbers at V and W interfaces for Thuburn monotonic scheme

  if (nl%iscal_monot == 1) then
     call comp_cfl2( mrl, dtlm, vmsca, wmsca, rho_old, rho, dzps_v, &
                     iwrecv, krecw )
  endif

! LOOP OVER SCALARS HERE
! (skip n=1 which is THIL and computed elsewhere)

  do n = 2, num_scalar

! Point SCP and SCT to scalar table arrays

     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

! Evaluate T3D gradient of scalar field

     call grad_t3d(mrl, scp, gxps_scp, gyps_scp, gzps_scp)

! MPI send of SCP gradient components

     if (iparallel == 1) then
        call mpi_send_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
                             rvara3=gzps_scp)
     endif

! Horizontal loop over all primary W columns to compute the
! upwinded scalar value at each W interface. This can be 
! computed before scalar gradients are received at the borders

! Loop over all primary W/T columns:

     !$omp parallel do private(iw,k,kd)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        ! Loop over W levels
        do k = lpw(iw), mza-1
           kd = kdepw(k,iw)

           ! First-order upstream:
           ! scp_upw(k,iw) = scp(kd,iw)

           scp_upw(k,iw) = scp(kd,iw)                     &
                         + dxps_w(k,iw) * gxps_scp(kd,iw) &
                         + dyps_w(k,iw) * gyps_scp(kd,iw) &
                         + dzps_w(k,iw) * gzps_scp(kd,iw)
        enddo
     enddo
     !$omp end parallel do

! Compute bounds on the scalar values at each W level for the Thuburn
! flux limiter (can be done before gradient communication is finished)
     
     if (nl%iscal_monot == 1) then
        call comp_vert_limits(mrl, scp, scp_upw, iwdepv, iwrecv, kdepw, krecw)
     endif

! MPI recv of SCP gradient components

     if (iparallel == 1) then
        call mpi_recv_w(mrl, rvara1=gxps_scp, rvara2=gyps_scp, &
                             rvara3=gzps_scp)
     endif

     call lbcopy_w(mrl, a1=gxps_scp, a2=gyps_scp, a3=gzps_scp)

! Horizontal loop over all V points surrounding primary W/T columns
! to compute upwinded scalar value at the V points

     !$omp parallel do private(iv,k,iwd)
     do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

        do k = lpv(iv), mza
           iwd = iwdepv(k,iv)

           ! 1st-ORDER UPWIND:         
           ! scp_upv(k,iv) = scp(k,iwd)

           ! Miura (2007, JAS):
           scp_upv(k,iv) = scp(k,iwd)                       &
                           + dxps_v(k,iv) * gxps_scp(k,iwd) &
                           + dyps_v(k,iv) * gyps_scp(k,iwd) &
                           + dzps_v(k,iv) * gzps_scp(k,iwd)

        enddo
     enddo

! Compute bounds on the scalar values at each V interface 
! for the Thuburn flux limiter

     if (nl%iscal_monot == 1) then
        call comp_horiz_limits(mrl, scp, scp_upv, iwdepv, iwrecv)
     endif


! Limit the advective fluxes if using monotonic advection

     if (nl%iscal_monot == 1) then
        call apply_flux_limiters(mrl, kdepw, iwdepv, scp_upw, scp_upv)
     endif

! Horizontal loop over W/T points

     !$omp parallel do private (iw, kb, hfluxadv, hfluxdif, jv, iv, iwn, &
     !$omp                      dirv, k, vfluxadv)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        kb = lpw(iw)

! Loop over neighbor V points of this W cell

        hfluxadv(1:mza) = 0.
        hfluxdif(1:mza) = 0.

        do jv = 1, itab_w(iw)%npoly

           iv   = itab_w(iw)%iv(jv)
           iwn  = itab_w(iw)%iw(jv)
           dirv = itab_w(iw)%dirv(jv)

! Vertical loop over T levels 

           do k = lpv(iv), mza

! Horizontal advective and diffusive scalar fluxes

              hfluxadv(k) = hfluxadv(k) + dirv * vmsca(k,iv) * scp_upv(k,iv)

              hfluxdif(k) = hfluxdif(k) + akhodx(k,iv) * (scp(k,iwn) - scp(k,iw))
           enddo
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

           sct(k,iw) = sct(k,iw) + volti(k,iw) &
                     * ( hfluxadv(k) + vfluxadv(k-1) - vfluxadv(k) &
                       + hfluxdif(k) )
        enddo

     enddo
     !$omp end parallel do

  enddo ! n

end subroutine scalar_transport


!=========================================================================


subroutine grad_t3d(mrl, scp, gxps, gyps, gzps)

  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
  use mem_grid,   only: mza, mwa, lpw, lpv, zt, zm, dzim, gxps_coef, gyps_coef
  use misc_coms,  only: io6
  use max_dims,   only: maxgrds

  implicit none

  integer, intent(in) :: mrl
  real,    intent(in) :: scp(mza,mwa)
  real,    intent(out) :: gxps(mza,mwa)
  real,    intent(out) :: gyps(mza,mwa)
  real,    intent(out) :: gzps(mza,mwa)

  integer :: j, iw, npoly, kb, jw1, jw2, iw1, iw2, iv1, iv2, k
  real    :: gxps1, gyps1, gxps2, gyps2
  real    :: gwz(mza)

! Horizontal loop over W columns for BEGS

!----------------------------------------------------------------------
  !$omp parallel do private(iw,npoly,kb,jw1,jw2,iw1,iw2,iv1, &
  !$omp                     iv2,gxps1,gyps1,gxps2,gyps2,k,gwz)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

     npoly = itab_w(iw)%npoly
     kb = lpw(iw)
   
     gxps(:,iw) = 0.
     gyps(:,iw) = 0.

! Loop over W neighbors of this W cell

     do jw1 = 1, npoly

        iw1 = itab_w(iw)%iw(jw1)
        iv1 = itab_w(iw)%iv(jw1)

! Vertical loop over T levels
! Zero-gradient lateral B.C. below lpv(iv1)

        do k = lpv(iv1), mza
           gxps(k,iw) = gxps(k,iw) + gxps_coef(iw,jw1) * (scp(k,iw1) - scp(k,iw))
           gyps(k,iw) = gyps(k,iw) + gyps_coef(iw,jw1) * (scp(k,iw1) - scp(k,iw))
        enddo

     enddo

! Vertical loop over W levels

     do k = kb, mza-1
        gwz(k) = dzim(k) * (scp(k+1,iw) - scp(k,iw))
     enddo

! Constant gradient top and bottom:
     gwz(kb-1)  = gwz(kb)
     gwz(mza) = gwz(mza-1)

! Zero-gradient top and bottom:
!    gwz(kb-1)  = 0.0
!    gwz(mza) = 0.0

! Vertical loop over T levels
 
     do k = kb, mza
        gzps(k,iw) = 0.5 * (gwz(k-1) + gwz(k))
     enddo

  enddo
  !$omp end parallel do

end subroutine grad_t3d


!===============================================================================


subroutine donorpointv(ldt, mrl, vs, vxe, vye, vze, iwdepv, iwrecv, &
                       dxps_v, dyps_v, dzps_v)

  use mem_ijtabs,  only: jtab_v, itab_v, jtv_wadj
  use mem_grid,    only: mza, mva, mwa, lpv, zt, zm, &
                         unx, uny, unz, xev, yev, zev
  use misc_coms,   only: io6, dtlm, dtsm, mdomain
  use max_dims,    only: maxgrds
  use consts_coms, only: eradi

  implicit none

  integer, intent(in)  :: ldt, mrl
  real,    intent(in)  :: vs (mza,mva)
  real,    intent(in)  :: vxe(mza,mwa)
  real,    intent(in)  :: vye(mza,mwa)
  real,    intent(in)  :: vze(mza,mwa)

  integer, intent(out) :: iwdepv(mza,mva)
  integer, intent(out) :: iwrecv(mza,mva)
  real,    intent(out) :: dxps_v(mza,mva)
  real,    intent(out) :: dyps_v(mza,mva)
  real,    intent(out) :: dzps_v(mza,mva)

  real :: dto2
  real :: dxps1, dyps1, dxps2, dyps2
  real :: cosv1, sinv1, cosv2, sinv2
  real :: wnx_v, wny_v, wnz_v
  real :: vxeface, vyeface, vzeface

  real :: dtm(maxgrds)
  real :: ufacev(mza), wfacev

  integer :: j, kb, k, iv, iw1, iw2

  if (ldt == 1) then
     dtm(:) = dtlm(:)
  else
     dtm(:) = dtsm(:)
  endif

! Horizontal loop over V points

!----------------------------------------------------------------------------
  !$omp parallel do private(iv,iw1,iw2,kb,k,dto2,dxps1,dyps1,dxps2,dyps2, &
  !$omp                     cosv1,sinv1,cosv2,sinv2,wnx_v,wny_v,wnz_v, &
  !$omp                     vxeface,vyeface,vzeface,ufacev,wfacev)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)
  iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------------

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

        dzps_v(k,iv) = -dto2 * wfacev

     enddo

! Compute X, Y displacement components for V face relative to T point
! Vertical loop over T/V levels

     do k = kb, mza

        if (vs(k,iv) > 0.0) then

           iwdepv(k,iv) = iw1
           iwrecv(k,iv) = iw2

           dxps_v(k,iv) = -dto2 * (vs(k,iv) * cosv1 - ufacev(k) * sinv1) + dxps1
           dyps_v(k,iv) = -dto2 * (vs(k,iv) * sinv1 + ufacev(k) * cosv1) + dyps1

        else

           iwdepv(k,iv) = iw2
           iwrecv(k,iv) = iw1

           dxps_v(k,iv) = -dto2 * (vs(k,iv) * cosv2 - ufacev(k) * sinv2) + dxps2
           dyps_v(k,iv) = -dto2 * (vs(k,iv) * sinv2 + ufacev(k) * cosv2) + dyps2

        endif

     enddo

  enddo
  !$omp end parallel do

  return
end subroutine donorpointv


!===============================================================================


subroutine donorpointw(ldt, mrl, ws, vxe, vye, vze, kdepw, krecw, &
                       dxps_w, dyps_w, dzps_w)

  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
  use mem_grid,   only: mza, mwa, lpw, dzt_top, dzt_bot
  use misc_coms,  only: io6, dtlm, dtsm
  use max_dims,   only: maxgrds

  implicit none

  integer, intent(in)  :: ldt, mrl
  real,    intent(in)  :: ws (mza,mwa)
  real,    intent(in)  :: vxe(mza,mwa)
  real,    intent(in)  :: vye(mza,mwa)
  real,    intent(in)  :: vze(mza,mwa)

  integer, intent(out) :: kdepw(mza,mwa)
  integer, intent(out) :: krecw(mza,mwa)
  real,    intent(out) :: dxps_w(mza,mwa)
  real,    intent(out) :: dyps_w(mza,mwa)
  real,    intent(out) :: dzps_w(mza,mwa)

  real :: dto2
  real :: unx_w, uny_w, vnx_w, vny_w, vnz_w
  real :: vxeface, vyeface, vzeface
  real :: uface,vface
  real :: dtm(maxgrds)

  integer :: j, kb, k, kp, iw

  if (ldt == 1) then
     dtm(:) = dtlm(:)
  else
     dtm(:) = dtsm(:)
  endif

! Horizontal loop over all primary W points

!------------------------------------------------------------------------
  !$omp parallel do private(iw,kb,k,kp,dto2,unx_w,uny_w,vnx_w,vny_w,vnz_w, &
  !$omp                     vxeface,vyeface,vzeface,uface,vface) 
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!------------------------------------------------------------------------

     kb = lpw(iw)

     dto2 = 0.5 * dtm(itab_w(iw)%mrlw)

     unx_w = itab_w(iw)%unx_w
     uny_w = itab_w(iw)%uny_w
   
     vnx_w = itab_w(iw)%vnx_w
     vny_w = itab_w(iw)%vny_w
     vnz_w = itab_w(iw)%vnz_w

! Vertical loop over W levels

     do k = kb, mza
        kp = min(k+1,mza)

! Average 3 velocity components from T points to W point

        vxeface = .5 * (vxe(k,iw) + vxe(kp,iw))
        vyeface = .5 * (vye(k,iw) + vye(kp,iw))
        vzeface = .5 * (vze(k,iw) + vze(kp,iw))

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

        if (ws(k,iw) > 0.0) then

           kdepw(k,iw) = k
           krecw(k,iw) = k + 1

           dzps_w(k,iw) = max( -dto2 * ws(k,iw) + dzt_top(k), 0.0)

        else

           kdepw(k,iw) = k + 1
           krecw(k,iw) = k

           dzps_w(k,iw) = min( -dto2 * ws(k,iw) - dzt_bot(k), 0.0)

        endif

     enddo
   
     kdepw (mza,iw) = mza
     krecw (mza,iw) = mza
     dzps_w(mza,iw) = dzt_top(mza)

  enddo
  !$omp end parallel do

  return
end subroutine donorpointw

!===========================================================================

subroutine zero_momsc(vmsc,wmsc,vxesc,vyesc,vzesc,rho_old)

use mem_ijtabs, only: jtab_v, jtab_w, istp, mrl_begl, jtv_wstn, jtw_wstn
use mem_grid,   only: mza, mva, mwa, lpw
use mem_basic,  only: rho
use misc_coms,  only: io6
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

! Zero out long timestep mass flux components (used for scalar advective 
! transport) so they may be summed over small timesteps

if (mrl > 0) then

!----------------------------------------------------------------------
   !$omp parallel do
   do iv = 2, mva
!----------------------------------------------------------------------

      vmsc(:,iv) = 0.0

   enddo
   !$omp end parallel do

!----------------------------------------------------------------------
   !$omp parallel do private(k)
   do iw = 2, mwa
!----------------------------------------------------------------------

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
   !$omp end parallel do

endif

return
end subroutine zero_momsc

!===========================================================================

subroutine timeavg_momsc(vmsc,wmsc,vxesc,vyesc,vzesc)

use mem_ijtabs, only: jtab_v, itab_v, jtab_w, itab_w, istp, mrl_endl, &
                      jtv_prog, jtw_prog
use mem_grid,   only: mza, mva, mwa, nsw_max, lpv, lpw
use misc_coms,  only: io6, nacoust

implicit none

real, intent(inout) :: vmsc(mza,mva)
real, intent(inout) :: wmsc(mza,mwa)

real, intent(inout) :: vxesc(mza,mwa)
real, intent(inout) :: vyesc(mza,mwa)
real, intent(inout) :: vzesc(mza,mwa)

integer :: k,iv,iw,mrl,mrlv,mrlw
real :: acoi

mrl = mrl_endl(istp)
if (mrl > 0) then

!----------------------------------------------------------------------
   !$omp parallel do private(mrlv,acoi)
   do iv = 2, mva
!----------------------------------------------------------------------

      mrlv = itab_v(iv)%mrlv
      acoi = 1.0 / real(nacoust(mrlv))

      vmsc(:,iv) = vmsc(:,iv) * acoi

   enddo
   !$omp end parallel do

!----------------------------------------------------------------------
   !$omp parallel do private(mrlw,acoi,k)
   do iw = 2, mwa
!----------------------------------------------------------------------

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
   !$omp end parallel do

endif

return
end subroutine timeavg_momsc
