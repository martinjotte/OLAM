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
subroutine scalar_transport(vmsc, wmsc, rho_old)

  use mem_ijtabs,   only: istp, jtab_v, jtab_w, mrl_endl, itab_v, itab_w
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, lsw, zt, zm, dzim, &
                          dniv, volt, arv, arw, dzim, volti
  use misc_coms,    only: io6, dtlm, iparallel
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: vkh, hkm, sxfer_rk, fqtpbl
  use massflux,     only: tridiffo
  use mem_basic,    only: rho
  use mem_para,     only: myrank
  use oname_coms,   only: nl
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_thuburn

  !$ use omp_lib

  implicit none

  real,    intent(in) :: vmsc(mza,mva)
  real,    intent(in) :: wmsc(mza,mwa)
  real(8), intent(in) :: rho_old(mza,mwa)

  integer :: j,iw,iw1,iw2,iw3,iw4,iwd,iwr,mrl
  integer :: n,k,kb,kbv,ks,kd,kr,iv,iwn,jv,npoly
  real    :: dtl,dtli
  real    :: dirv,hdniv

! Automatic arrays:

  real :: vsc(mza,mva)
  real :: wsc(mza,mwa)

  real :: vxesc(mza,mwa) ! XE velocity component at T point
  real :: vyesc(mza,mwa) ! YE velocity component at T point
  real :: vzesc(mza,mwa) ! ZE velocity component at T point

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

  real :: akodz(mza), dtomass(mza)
  real :: vctr5(mza), vctr6(mza), vctr7(mza), vctr8(mza)
  real :: hfluxadv(mza), hfluxdif(mza), vfluxadv(mza), vfluxdif(mza)
  real :: del_scp(mza)

  real, pointer :: scp(:,:)
  real, pointer :: sct(:,:)

! Extra variables for Thuburn scheme
  
  real                 :: scp_vin_min(mza), scp_vin_max(mza)
  real                 :: cfl_vout, cfl_wout

! Return if this is not the end of the long timestep on any MRL

  mrl = mrl_endl(istp)
  if (mrl == 0) return

! Horizontal loop over all primary W columns to diagnose
! face-normal vertical velocity at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp parallel 
  !$omp do private(iw,kb,k) 
  do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)

     kb = lpw(iw)

     do k = kb, mza-2
        wsc(k,iw) = 2.0 * wmsc(k,iw) / (rho_old(k,iw) + rho_old(k+1,iw))
     enddo
   
     wsc(kb-1 ,iw) = wsc(kb,iw)
     wsc(mza-1,iw) = 0.
   
  enddo
  !$omp end do

! Loop over V columns (the immediate neighbors of all primary W points) to
! diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp do private(iv,iw1,iw2,kb,k) 
  do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     kb = lpv(iv)

     vsc(2:kb-1,iv) = 0.0

     do k = kb, mza-1
        vsc(k,iv) = 2.0 * vmsc(k,iv) / ( rho_old(k,iw1) + rho_old(k,iw2) )
     enddo

  enddo
  !$omp end do
  !$omp end parallel

! Diagnose 3D velocity at T points using velocities for scalar advection

  call vel_t3d_hex(mrl, vsc, wsc, vxesc, vyesc, vzesc)

! MPI send of VXESC, VYESC, VZESC

  if (iparallel == 1) then
     call mpi_send_w('V', vxe=vxesc, vye=vyesc, vze=vzesc)
  endif

! Diagnose advective donor point locations for all primary W faces
! No parallel communication is necessary to compute this

  call donorpointw(1, mrl, wsc, vxesc, vyesc, vzesc, kdepw, krecw, &
                   dxps_w, dyps_w, dzps_w)

! Complete MPI recv of VXESC, VYESC, VZESC

  if (iparallel == 1) then
     call mpi_recv_w('V', vxe=vxesc, vye=vyesc, vze=vzesc)
  endif

! Diagnose advective donor point locations for the V faces surrounding all
! primary W points. Communication of velocities must have been completed

  call donorpointv(1, mrl, vsc, vxesc, vyesc, vzesc, iwdepv, iwrecv, &
                   dxps_v, dyps_v, dzps_v)

! Diagnose CFL numbers at V and W interfaces for Thuburn flux limiter

  if (nl%iscal_monot == 1) then

     ! Begin OpenMP parallel block
     !$omp parallel 

     !$omp workshare
     cfl_out_sum(:,:) = 0.0
   ! cfl_in_sum (:,:) = 0.0
     !$omp end workshare 

     !$omp do private(iv,k,kbv,iwr,iwd,cfl_vout)
     do j = 1, jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)
        kbv = lpv(iv)
      
        do k = kbv, mza-1
           iwr = iwrecv(k,iv)
           iwd = iwdepv(k,iv)

           kdepv(k,iv) = merge( min(k+1,mza-1), max(k-1,kbv), dzps_v(k,iv) >= 0.0)

           cfl_vin(k,iv) = abs(vmsc(k,iv)) * arv(k,iv) * dtlm(itab_w(iwr)%mrlw) &
                           * volti(k,iwr) / rho_old(k,iwr) 

           cfl_vout      = abs(vmsc(k,iv)) * arv(k,iv) * dtlm(itab_w(iwd)%mrlw) &
                           * volti(k,iwd) / rho_old(k,iwd) 

         ! Not needed for scalars  
         ! cfl_in_sum (k,iwr) = cfl_in_sum (k,iwr) + cfl_vin(k,iv)

           cfl_out_sum(k,iwd) = cfl_out_sum(k,iwd) + cfl_vout
        enddo
     enddo
     !$omp end do

     !$omp do private(iw,kb,dtl,k,kd,kr,cfl_wout)
     do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)
        kb = lpw(iw)

        ! note: use dtsm for short time step
        dtl = dtlm(itab_w(iw)%mrlw)
      
        do k = kb, mza-2
           kd = kdepw(k,iw)
           kr = krecw(k,iw)

           cfl_win(k,iw) = abs(wmsc(k,iw)) * arw(k,iw) * dtl * volti(kr,iw) &
                           / rho_old(kr,iw)

           cfl_wout      = abs(wmsc(k,iw)) * arw(k,iw) * dtl * volti(kd,iw) &
                           / rho_old(kd,iw)

         ! Not needed for scalars
         ! cfl_in_sum (kr,iw) = cfl_in_sum (kr,iw) + cfl_win(k,iw)

           cfl_out_sum(kd,iw) = cfl_out_sum(kd,iw) + cfl_wout
        enddo

        cfl_win(mza-1,iw) = 0.0

        do k = kb, mza-1
           tfact(k,iw) = rho(k,iw) / rho_old(k,iw)

         ! if we don't have future rho (for thil, vxe, vye, vze):
         ! tfact(k,iw) = 1.0 + cfl_in_sum(k,iw) - cfl_out_sum(k,iw)

           ! cfl_out_sum is reused as 1 / cfl_out_sum
           cfl_out_sum(k,iw) = 1.0 / max( cfl_out_sum(k,iw), 1.e-6)
        enddo

     enddo
     !$omp end do

     ! End OpenMP parallel block
     !$omp end parallel

  endif ! monotonic

! LOOP OVER SCALARS HERE

  do n = 1, num_scalar

! Point SCP and SCT to scalar table arrays

     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

! Evaluate T3D gradient of scalar field

     call grad_t3d(mrl, scp, gxps_scp, gyps_scp, gzps_scp)

! MPI send of SCP gradient components

     if (iparallel == 1) then
        call mpi_send_w('G', gxps_scp=gxps_scp, gyps_scp=gyps_scp, &
                             gzps_scp=gzps_scp)
     endif

! Compute bounds on the scalar values at each W level for the Thuburn
! flux limiter (can be done before gradient communication is finished)

     ! Begin OpenMP parallel block
     !$omp parallel

     if (nl%iscal_monot == 1) then

        !$omp workshare
        scp_local_min(:,:) = scp(:,:)
        scp_local_max(:,:) = scp(:,:)

        scp_in_max(:,:) = scp(:,:)
        scp_in_min(:,:) = scp(:,:)

        c_scp_in_max_sum(:,:) = 0.0
        c_scp_in_min_sum(:,:) = 0.0
        !$omp end workshare

! Expand inflow bounds at each W level based on the transverse cells 
! in the upstream neighborhood of each W level. 
! Loop over all immediate V neighbors of each primary W/T columns:

        !$omp do private(iv,k,iwd,iwr)
        do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)
           do k = lpv(iv), mza-1
              iwd = iwdepv(k,iv)
              iwr = iwrecv(k,iv)
              scp_in_max(k,iwr) = max(scp_in_max(k,iwr), scp(k,iwd))
              scp_in_min(k,iwr) = min(scp_in_min(k,iwr), scp(k,iwd))
           enddo
        enddo
        !$omp end do

     endif ! monotonic

! Horizontal loop over all primary W columns to compute the
! upwinded scalar value at each W interface. This can be 
! computed before scalar gradients are received at the borders

! Loop over all primary W/T columns:

     !$omp do private(iw,k,kr,kd)
     do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)

        do k = lpw(iw), mza-1
           kd = kdepw(k,iw)

           ! First-order upstream:
           ! scp_upw(k,iw) = scp(kd,iw)

           scp_upw(k,iw) = scp(kd,iw)                       &
                           + dxps_w(k,iw) * gxps_scp(kd,iw) &
                           + dyps_w(k,iw) * gyps_scp(kd,iw) &
                           + dzps_w(k,iw) * gzps_scp(kd,iw)
        enddo

        if (nl%iscal_monot == 1) then

! Make sure the upwinded scalar value at each W level is properly bounded

           do k = lpw(iw), mza-1
              kr = krecw(k,iw)
              kd = kdepw(k,iw)
              scp_upw(k,iw) = min( scp_upw(k,iw), max( scp_in_max(kd,iw), scp(kr,iw)))
              scp_upw(k,iw) = max( scp_upw(k,iw), min( scp_in_min(kd,iw), scp(kr,iw)))
           enddo

! Compute the W contribution to the max/min allowed values in each grid box, and
! the sum of the max and min inputs to each grid box

           do k = lpw(iw), mza-1
              kr = krecw(k,iw)
              kd = kdepw(k,iw)

              scp_local_min(kr,iw) = min( scp_local_min(kr,iw), scp_in_min(kd,iw))
              scp_local_max(kr,iw) = max( scp_local_max(kr,iw), scp_in_max(kd,iw))

              c_scp_in_max_sum(kr,iw) = c_scp_in_max_sum(kr,iw) + cfl_win(k,iw) * &
                                        max( scp_in_max(kd,iw), scp_upw(k,iw))

              c_scp_in_min_sum(kr,iw) = c_scp_in_min_sum(kr,iw) + cfl_win(k,iw) * &
                                        min( scp_in_min(kd,iw), scp_upw(k,iw))
           enddo

        endif ! monotonic

     enddo
     !$omp end do

     ! End OpenMP parallel block
     !$omp end parallel

! MPI recv of SCP gradient components

     if (iparallel == 1) then
        call mpi_recv_w('G', gxps_scp=gxps_scp, gyps_scp=gyps_scp, &
                             gzps_scp=gzps_scp)
     endif

! Horizontal loop over all V points surrounding primary W/T columns
! to compute upwinded scalar value at the V points

     !$omp parallel do private(iv,k,kd,kbv,iwd,iwr,iw3,iw4,scp_vin_min,scp_vin_max)
     do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)

        kbv = lpv(iv)
        do k = kbv, mza-1
           iwd = iwdepv(k,iv)

           ! 1st-ORDER UPWIND:         
           ! scp_upv(k,iv) = scp(k,iwd)

           ! Miura (2007, JAS):
           scp_upv(k,iv) = scp(k,iwd)                       &
                           + dxps_v(k,iv) * gxps_scp(k,iwd) &
                           + dyps_v(k,iv) * gyps_scp(k,iwd) &
                           + dzps_v(k,iv) * gzps_scp(k,iwd)
        enddo

! Compute bounds on the scalar values at each V interface 
! for the Thuburn flux limiter

        if (nl%iscal_monot == 1) then

           iw3 = itab_v(iv)%iw(3)
           iw4 = itab_v(iv)%iw(4)

           ! Compute inflow bounds at each V interface. Here we include the 
           ! immediate upwind cell and vertically upstream-diagonal neighbor

           do k = kbv, mza-1
              iwd = iwdepv(k,iv)
              kd  =  kdepv(k,iv)
              scp_vin_min(k) = min( scp(k,iwd), scp(kd,iwd))
              scp_vin_max(k) = max( scp(k,iwd), scp(kd,iwd))
           enddo

           ! Include contribution of lateral neighbor iw3 to inflow bounds
      
           do k = max(lpw(iw3),kbv), mza-1
              scp_vin_min(k) = min( scp_vin_min(k), scp(k,iw3))
              scp_vin_max(k) = max( scp_vin_max(k), scp(k,iw3))
           enddo

           ! Include contribution of lateral neighbor iw4 to inflow bounds
      
           do k = max(lpw(iw4),kbv), mza-1
              scp_vin_min(k) = min( scp_vin_min(k), scp(k,iw4))
              scp_vin_max(k) = max( scp_vin_max(k), scp(k,iw4))
           enddo

           ! Make sure the upwinded scalar value at each V interface is properly bounded

           do k = kbv, mza-1
              iwr = iwrecv(k,iv)
              scp_upv(k,iv) = min(scp_upv(k,iv), max(scp_vin_max(k), scp(k,iwr)))
              scp_upv(k,iv) = max(scp_upv(k,iv), min(scp_vin_min(k), scp(k,iwr)))
           enddo

           ! Compute the V contribution to the max/min allowed values in each grid box
           ! and the sum of the max and min inputs to each grid box

           do k = kbv, mza-1
              iwr = iwrecv(k,iv)

              scp_local_min(k,iwr) = min(scp_local_min(k,iwr), scp_vin_min(k))
              scp_local_max(k,iwr) = max(scp_local_max(k,iwr), scp_vin_max(k))
         
              c_scp_in_max_sum(k,iwr) = c_scp_in_max_sum(k,iwr) + cfl_vin(k,iv) * &
                                        max( scp_vin_max(k), scp_upv(k,iv))

              c_scp_in_min_sum(k,iwr) = c_scp_in_min_sum(k,iwr) + cfl_vin(k,iv) * &
                                        min( scp_vin_min(k), scp_upv(k,iv))
           enddo

        endif ! monotonic

     enddo
     !$omp end parallel do

!  Thuburn limiter: compute scalar out min,max for each cell
!  Horizontal loop over all primary W/T columns

     if (nl%iscal_monot == 1) then

        !$omp parallel do private(iw,k)
        do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)

           ! Vertical loop over T levels
           do k = lpw(iw), mza-1

              scp_out_min(k,iw) = (scp(k,iw) + c_scp_in_max_sum(k,iw) - &
                   scp_local_max(k,iw) * tfact(k,iw)) * cfl_out_sum(k,iw)

              scp_out_max(k,iw) = (scp(k,iw) + c_scp_in_min_sum(k,iw) - &
                   scp_local_min(k,iw) * tfact(k,iw)) * cfl_out_sum(k,iw)

           enddo

        enddo
        !$omp end parallel do

        ! MPI send of scalar max/min outflow values

        if (iparallel == 1) then
           call mpi_send_w('G', gxps_scp=scp_out_min, gyps_scp=scp_out_max)
        endif

! Limit the vertical fluxes based on the computed scalar outgoing max/min values
! Can be done before outgoing max/mins are received at the boundary cells
! Horizontal loop over all primary W/T columns
      
        !$omp parallel do private(iw,kd,k)
        do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)

           ! Vertical loop over W levels
           do k = lpw(iw), mza-2
              kd = kdepw(k,iw)
              scp_upw(k,iw) = min( scp_upw(k,iw), scp_out_max(kd,iw))
              scp_upw(k,iw) = max( scp_upw(k,iw), scp_out_min(kd,iw))
           enddo

        enddo
        !$omp end parallel do
      
! MPI receive of scalar max/min outflow values

        if (iparallel == 1) then
           call mpi_recv_w('G', gxps_scp=scp_out_min, gyps_scp=scp_out_max)
        endif

! Limit the horizontal fluxes based on the computed scalar outgoing max/min
! Horizontal loop over all V points bordering primary W/T columns

        !$omp parallel do private(iv,iwd,k) 
        do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)

           ! Vertical loop over T levels
           do  k = lpv(iv), mza-1
              iwd = iwdepv(k,iv)
              scp_upv(k,iv) = min( scp_upv(k,iv), scp_out_max(k,iwd) )
              scp_upv(k,iv) = max( scp_upv(k,iv), scp_out_min(k,iwd) )
           enddo

        enddo
        !$omp end parallel do

     endif ! monotonic

! Horizontal loop over W/T points

     !$omp parallel do private (iw,kb,k,npoly,jv,iv,iwn,kbv,iwd,kd,ks,   &
     !$omp                      dtl,dtli,dtomass,hfluxadv,hfluxdif,dirv, &
     !$omp                      vfluxadv,akodz,vctr5,vctr7,vctr6,vctr8,  &
     !$omp                      hdniv,del_scp,vfluxdif)
     do j = 1,jtab_w(26)%jend(mrl); iw = jtab_w(26)%iw(j)

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

        do jv = 1, npoly
           iv  = itab_w(iw)%iv(jv)
           iwn = itab_w(iw)%iw(jv)

           kbv = lpv(iv)

           dirv = itab_w(iw)%dirv(jv)

           hdniv = .5 * dniv(iv)
   
! Vertical loop over T levels 

           do k = kbv, mza-1

! Horizontal advective and diffusive scalar fluxes

              hfluxadv(k) = hfluxadv(k) &
                          + dirv * vmsc(k,iv) * arv(k,iv) * scp_upv(k,iv)

              hfluxdif(k) = hfluxdif(k) &
                          + hdniv * arv(k,iv) * (hkm(k,iwn) + hkm(k,iw)) &
                          * (scp(k,iwn) - scp(k,iw))
              
           enddo
        enddo

! Vertical loop over W levels

        do k = kb, mza-2

           vfluxadv(k) = wmsc(k,iw) * arw(k,iw) * scp_upw(k,iw)

! Prepare for vertical diffusion - Fill tri-diagonal matrix coefficients

           akodz(k) = arw(k,iw) * vkh(k,iw) * dzim(k)
           vctr5(k) = - akodz(k) * dtomass(k)
           vctr7(k) = - akodz(k) * dtomass(k+1)
           vctr6(k) = 1. - vctr5(k) - vctr7(k)
           vctr8(k) = akodz(k) * (scp(k,iw) - scp(k+1,iw))
        enddo

! Special case for total water SH_W: Apply surface vapor flux

        if (scalar_tab(n)%name == 'SH_W') then

           if (allocated(fqtpbl)) fqtpbl(:,iw) = 0.0

! Vertical loop over T levels that are adjacent to surface

           do ks = 1,lsw(iw)
              k = kb + ks - 1

! Apply surface vapor xfer [kg_vap] directly to SCT [kg_vap / (m^3 s)]

              sct(k,iw) = sct(k,iw) + dtli * volti(k,iw) * sxfer_rk(ks,iw)

              if (allocated(fqtpbl)) then
                 fqtpbl(k,iw)  = dtli * volti(k,iw) * sxfer_rk(ks,iw) / rho(k,iw)
              endif

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

! Grell scheme needs PBL moisture tendency stored

        if (allocated(fqtpbl) .and. scalar_tab(n)%name == 'SH_W') then
           
           do k = kb,mza-1
              fqtpbl(k,iw) = fqtpbl(k,iw) + volti(k,iw) * (vfluxdif(k-1) - vfluxdif(k)) / rho(k,iw)
           enddo

        endif

     enddo
     !$omp end parallel do

  enddo ! n

end subroutine scalar_transport


!=========================================================================


subroutine grad_t3d(mrl, scp, gxps, gyps, gzps)

  use mem_ijtabs, only: jtab_w, itab_w
  use mem_grid,   only: mza, mwa, lpw, zt, zm, dzim, vnx, vny, vnz, wnx, wny, wnz
  use misc_coms,  only: io6
  use max_dims,   only: maxgrds

  !$ use omp_lib

  implicit none

  integer, intent(in) :: mrl
  real,    intent(in) :: scp(mza,mwa)
  real,    intent(out) :: gxps(mza,mwa)
  real,    intent(out) :: gyps(mza,mwa)
  real,    intent(out) :: gzps(mza,mwa)

  integer :: j,iw,npoly,kb,jw1,jw2,iw1,iw2,kb1,kb2,k,kbot,kbv

  real :: farm,wti
  real :: gxps1,gyps1,gxps2,gyps2

  real :: wt(mza), gwz(mza)

! Horizontal loop over W columns for BEGS

  call psub()
!----------------------------------------------------------------------
  !$omp parallel do private(iw,npoly,kb,jw1,jw2,iw1,iw2,kb1,kb2,k,kbot,&
  !$omp                     kbv,wt,farm,gxps1,gyps1,gxps2,gyps2,gwz,wti)
  do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
  call qsub('W',iw)

     npoly = itab_w(iw)%npoly
     kb = lpw(iw)
   
     gxps(:,iw) = 0.
     gyps(:,iw) = 0.
   
     wt(:) = 0.

     kbot = mza-1

! Loop over W neighbors of this W cell

     do jw1 = 1, npoly

        jw2 = mod(jw1,npoly) + 1

        iw1 = itab_w(iw)%iw(jw1)
        iw2 = itab_w(iw)%iw(jw2)
      
        kb1 = lpw(iw1)
        kb2 = lpw(iw2)
      
        farm = itab_w(iw)%farm(jw1)

        gxps1 = itab_w(iw)%gxps1(jw1)
        gyps1 = itab_w(iw)%gyps1(jw1)

        gxps2 = itab_w(iw)%gxps2(jw1)
        gyps2 = itab_w(iw)%gyps2(jw1)

        kbv  = max(kb,kb1,kb2)
        kbot = min(kbv,kbot)

! Vertical loop over T levels

        do k = kbv, mza-1

!!      do k = kb,mza-1
!!      
!! Gradient calculation when all 3 W points are above ground
!!
!!         if (k >= kb1 .and. k >= kb2) then

           gxps(k,iw) = gxps(k,iw) + farm * (gxps1 * (scp(k,iw1) - scp(k,iw)) &
                                           + gxps2 * (scp(k,iw2) - scp(k,iw)))

           gyps(k,iw) = gyps(k,iw) + farm * (gyps1 * (scp(k,iw1) - scp(k,iw)) &
                                           + gyps2 * (scp(k,iw2) - scp(k,iw)))

           wt(k) = wt(k) + farm
            
!!         elseif (k >= kb1) then
!!
!! Gradient calculation when IW2 point is below ground
!!
!!         elseif (k >= kb2) then
!!
!! Gradient calculation when IW1 point is below ground
!!
!!         endif
         
        enddo
      
     enddo

! Vertical loop over T levels to normalize horizontal gradients.
! Below kbot, flow is dominated by terrain and 1st-order upstream
! is probably sufficient (gxps=0 and gyps=0 below kbot)

     do k = kbot, mza-1
        wti = 1.0 / wt(k)
        gxps(k,iw) = gxps(k,iw) * wti
        gyps(k,iw) = gyps(k,iw) * wti
     enddo

! Vertical loop over W levels

     do k = kb, mza-2
        gwz(k) = dzim(k) * (scp(k+1,iw) - scp(k,iw))
     enddo
   
     gwz(kb-1)  = gwz(kb)
     gwz(mza-1) = gwz(mza-2)

! Vertical loop over T levels
     
     do k = kb, mza-1
        gzps(k,iw) = 0.5 * (gwz(k-1) + gwz(k))
     enddo
   
  enddo
  !$omp end parallel do
  call rsub('Wa',16)

  return
end subroutine grad_t3d


!===============================================================================


subroutine donorpointv(ldt, mrl, vs, vxe, vye, vze, iwdepv, iwrecv, &
                       dxps_v, dyps_v, dzps_v)

  use mem_ijtabs,  only: jtab_v, itab_v
  use mem_grid,    only: mza, mva, mwa, lpv, zt, zm, &
                         unx, uny, unz, xev, yev, zev
  use misc_coms,   only: io6, dtlm, dtsm
  use max_dims,    only: maxgrds
  use consts_coms, only: eradi

  !$ use omp_lib

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

  call psub()
!----------------------------------------------------------------------------
  !$omp parallel do private(iv,iw1,iw2,kb,k,dto2,dxps1,dyps1,dxps2,dyps2, &
  !$omp                     cosv1,sinv1,cosv2,sinv2,wnx_v,wny_v,wnz_v, &
  !$omp                     vxeface,vyeface,vzeface,ufacev,wfacev)
  do j = 1,jtab_v(12)%jend(mrl); iv = jtab_v(12)%iv(j)
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

     do k = kb, mza-1

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

     do k = kb, mza-1

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
  call rsub('V',15)

  return
end subroutine donorpointv


!===============================================================================


subroutine donorpointw(ldt, mrl, ws, vxe, vye, vze, kdepw, krecw, &
                       dxps_w, dyps_w, dzps_w)

  use mem_ijtabs, only: jtab_w, itab_w
  use mem_grid,   only: mza, mwa, lpw, zt, zm
  use misc_coms,  only: io6, dtlm, dtsm
  use max_dims,   only: maxgrds

  !$ use omp_lib

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

  integer :: j, kb, k, iw

  if (ldt == 1) then
     dtm(:) = dtlm(:)
  else
     dtm(:) = dtsm(:)
  endif

! Horizontal loop over all primary W points

  call psub()
!------------------------------------------------------------------------
  !$omp parallel do private(iw,kb,k,dto2,unx_w,uny_w,vnx_w,vny_w,vnz_w, &
  !$omp                     vxeface,vyeface,vzeface,uface,vface) 
  do j = 1,jtab_w(15)%jend(mrl); iw = jtab_w(15)%iw(j)
!------------------------------------------------------------------------
  call qsub('W',iw)

     kb = lpw(iw)

     dto2 = dtm(itab_w(iw)%mrlw)

     unx_w = itab_w(iw)%unx_w
     uny_w = itab_w(iw)%uny_w
   
     vnx_w = itab_w(iw)%vnx_w
     vny_w = itab_w(iw)%vny_w
     vnz_w = itab_w(iw)%vnz_w

! Vertical loop over W levels

     do k = kb, mza-1

! Average 3 velocity components from T points to W point

        vxeface = .5 * (vxe(k,iw) + vxe(k+1,iw))
        vyeface = .5 * (vye(k,iw) + vye(k+1,iw))
        vzeface = .5 * (vze(k,iw) + vze(k+1,iw))

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

     do k = kb, mza-2

        if (ws(k,iw) > 0.0) then

           kdepw(k,iw) = k
           krecw(k,iw) = k + 1

           dzps_w(k,iw) = -dto2 * ws(k,iw) + zm(k) - zt(k)

        else

           kdepw(k,iw) = k + 1
           krecw(k,iw) = k

           dzps_w(k,iw) = -dto2 * ws(k,iw) + zm(k) - zt(k+1)

        endif

     enddo
   
     kdepw (mza-1,iw) = mza-1
     krecw (mza-1,iw) = mza
     dzps_w(mza-1,iw) = zm(mza-1) - zt(mza-1)

  enddo
  call rsub('V',15)

  return
end subroutine donorpointw
