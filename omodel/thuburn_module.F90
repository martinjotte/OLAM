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

Module mem_thuburn
  implicit none

  real,    allocatable, private :: scp_local_min(:,:)
  real,    allocatable, private :: scp_local_max(:,:)
  real,    allocatable, private :: cfl_out_sum(:,:)
  real,    allocatable, private :: cfl_vin(:,:)
  real,    allocatable, private :: cfl_win(:,:)
  real,    allocatable, private :: c_scp_in_max_sum(:,:)
  real,    allocatable, private :: c_scp_in_min_sum(:,:)
  real,    allocatable, private :: scp_out_min(:,:)
  real,    allocatable, private :: scp_out_max(:,:)
  real,    allocatable, private :: scp_in_min(:,:)
  real,    allocatable, private :: scp_in_max(:,:)
  real,    allocatable, private :: tfact(:,:)
  integer, allocatable, private :: kdepv(:,:)

Contains

!===============================================================================

  subroutine alloc_thuburn(imonot, mza, mva, mwa)
    use misc_coms, only: rinit
    implicit none
    
    integer, intent(in) :: imonot, mza, mva, mwa

    allocate( cfl_out_sum(mza,mwa)) ; cfl_out_sum = rinit

    if (imonot == 1) then
       
       allocate( scp_local_min   (mza,mwa)) ; scp_local_min    = rinit
       allocate( scp_local_max   (mza,mwa)) ; scp_local_max    = rinit
       allocate( cfl_vin         (mza,mva)) ; cfl_vin          = rinit
       allocate( c_scp_in_max_sum(mza,mwa)) ; c_scp_in_max_sum = rinit
       allocate( c_scp_in_min_sum(mza,mwa)) ; c_scp_in_min_sum = rinit
       allocate( scp_out_min     (mza,mwa)) ; scp_out_min      = rinit
       allocate( scp_out_max     (mza,mwa)) ; scp_out_max      = rinit
       allocate( cfl_win         (mza,mwa)) ; cfl_win          = rinit
       allocate( scp_in_min      (mza,mwa)) ; scp_in_min       = rinit
       allocate( scp_in_max      (mza,mwa)) ; scp_in_max       = rinit
       allocate( tfact           (mza,mwa)) ; tfact            = rinit
       allocate( kdepv           (mza,mva)) ; kdepv            = 0

    endif

  end subroutine alloc_thuburn

!===============================================================================

  subroutine dealloc_thuburn()
    implicit none

    if (allocated( scp_local_min    )) deallocate( scp_local_min )
    if (allocated( scp_local_max    )) deallocate( scp_local_max )
    if (allocated( cfl_out_sum      )) deallocate( cfl_out_sum )
    if (allocated( cfl_vin          )) deallocate( cfl_vin )
    if (allocated( c_scp_in_max_sum )) deallocate( c_scp_in_max_sum )
    if (allocated( c_scp_in_min_sum )) deallocate( c_scp_in_min_sum )
    if (allocated( scp_out_min      )) deallocate( scp_out_min )
    if (allocated( scp_out_max      )) deallocate( scp_out_max )
    if (allocated( cfl_win          )) deallocate( cfl_win )
    if (allocated( scp_in_min       )) deallocate( scp_in_min )
    if (allocated( scp_in_max       )) deallocate( scp_in_max )
    if (allocated( tfact            )) deallocate( tfact )
    if (allocated( kdepv            )) deallocate( kdepv )

  end subroutine dealloc_thuburn

!===============================================================================

  subroutine comp_cfl1(mrl, dtm, vmca, wmca, vs, ws, rho_old, iwdepv, kdepw, do_check)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, glatw, glonw, mza, mva, mwa
    use consts_coms, only: r8
    use oname_coms,  only: nl
    use misc_coms,   only: io6, iparallel, time8p
    use mem_para,    only: myrank, mgroupsize
    use max_dims,    only: maxgrds

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer,           intent(in) :: mrl
    real,              intent(in) :: vmca   (mza,mva)
    real,              intent(in) :: wmca   (mza,mwa)
    real,              intent(in) :: vs     (mza,mva)
    real,              intent(in) :: ws     (mza,mwa)
    real(r8),          intent(in) :: dtm    (maxgrds)
    real(r8),          intent(in) :: rho_old(mza,mwa)
    integer,           intent(in) :: iwdepv (mza,mva)
    integer,           intent(in) :: kdepw  (mza,mwa)
    logical, optional, intent(in) :: do_check

    integer :: j, iv, k, kd, iw, iwd, npoly, n
    real    :: cfl_max, cfl_maxs(mgroupsize)
    integer :: ier, inode, imax(3), imaxs(3,mgroupsize)

    !$omp parallel
    !$omp workshare
    cfl_out_sum(:,:) = 0.0
    !$omp end workshare 

    !$omp do private(iv,k,iwd)
    do j = 1, jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Loop over T/V level
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          cfl_out_sum(k,iwd) = cfl_out_sum(k,iwd) + abs(vmca(k,iv))
       enddo

    enddo
    !$omp end do

    !$omp do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
  
       ! Loop over W/M level
       do k = lpw(iw), mza-1
          kd = kdepw(k,iw)
          cfl_out_sum(kd,iw) = cfl_out_sum(kd,iw) + abs(wmca(k,iw))
       enddo
     
       ! Loop over T/V level
       do k = lpw(iw), mza
          cfl_out_sum(k,iw) = cfl_out_sum(k,iw) * dtm(itab_w(iw)%mrlw) &
                            * volti(k,iw) / rho_old(k,iw)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    if (.not. present(do_check)) return
    if (.not.         do_check)  return

    ! Find the max CFL number on each node, and print an error message
    ! if the CFL number is greater than 1

    cfl_max = -1.0

    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
       do k = lpw(iw), mza

          if ( cfl_out_sum(k,iw) > cfl_max  .or.  &
               cfl_out_sum(k,iw) /= cfl_out_sum(k,iw) ) then
             cfl_max = cfl_out_sum(k,iw)
             imax    = (/ k, iw, itab_w(iw)%iwglobe /)
          endif

          if ( cfl_out_sum(k,iw) > 1.0  .or.  &
               cfl_out_sum(k,iw) /= cfl_out_sum(k,iw) ) then
             npoly = itab_w(iw)%npoly
             write(*,*)
             write(*,'(4(A,I0),2(A,f0.3),/,2(A,f0.3),A,7(f0.3,1x))')           &
                  "!!! CFL VIOLATION at node ", myrank,                        &
                  ", iw=", iw, ", iwglobe=", itab_w(iw)%iwglobe, ", k=", k,    &
                  " lat=", glatw(iw), " lon=", glonw(iw),                      &
                  "!!! CFL = ", cfl_out_sum(k,iw),                             &
                  ", W = ", ws(k,iw), ", VC = ", vs(k,itab_w(iw)%iv(1:npoly))
          endif

       enddo
    enddo

    ! Print the global max CFL number

    if (nl%cfl_prtfrq > 1.d-12) then
       if (mod(time8p,nl%cfl_prtfrq) < dtm(1)) then

#ifdef OLAM_MPI
          if (iparallel == 1) then
             call MPI_Gather(cfl_max, 1, MPI_REAL, cfl_maxs, 1, MPI_REAL, &
                             0, MPI_COMM_WORLD, ier)
             call MPI_Gather(imax, 3, MPI_INTEGER, imaxs, 3, MPI_INTEGER, &
                             0, MPI_COMM_WORLD, ier)
          endif
#endif

          if (myrank == 0) then
             inode = 0
             if (iparallel == 1) then
                if (any( cfl_maxs(:) /= cfl_maxs(:) )) then
                   do n = 1, mgroupsize
                      if (cfl_maxs(n) /= cfl_maxs(n)) then
                         inode = n
                         cfl_max = cfl_maxs(inode)
                         imax(:) = imaxs(:,inode)
                         exit
                      endif
                   enddo
                else
                   inode    = maxloc(cfl_maxs, dim=1)
                   cfl_max  = cfl_maxs(inode)
                   imax(:)  = imaxs(:,inode)
                endif
             endif
             write(*,'(5x,A,f0.3,3(A,I0))') "Max CFL# = ", cfl_max,  &
                  " at node ", inode-1, ", iwglobe=", imax(3), ", k=", imax(1)
          endif
       endif
    endif

  end subroutine comp_cfl1

!===============================================================================

  subroutine comp_cfl2(mrl, dtm, vmca, wmca, rho_old, rho, dzps_v, iwrecv, krecw)

    ! Diagnose inflow CFL numbers

    use mem_ijtabs,  only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use consts_coms, only: r8
    use oname_coms,  only: nl
    use misc_coms,   only: io6, iparallel, time8p
    use mem_para,    only: myrank, mgroupsize
    use max_dims,    only: maxgrds

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer,  intent(in) :: mrl
    real,     intent(in) :: vmca   (mza,mva)
    real,     intent(in) :: wmca   (mza,mwa)
    real(r8), intent(in) :: dtm    (maxgrds)
    real(r8), intent(in) :: rho_old(mza,mwa)
    real(r8), intent(in) :: rho    (mza,mwa)
    real,     intent(in) :: dzps_v (mza,mva)
    integer,  intent(in) :: iwrecv (mza,mva)
    integer,  intent(in) :: krecw  (mza,mwa)

    integer :: j, iv, k, kr, iw, iwr

    !$omp parallel 
    !$omp do private(iv,k,iwr)
    do j = 1, jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Loop over T/V levels
       do k = lpv(iv), mza
          iwr = iwrecv(k,iv)

          kdepv(k,iv) = merge( min(k+1,mza), max(k-1,lpv(iv)), dzps_v(k,iv) >= 0.0)

          cfl_vin(k,iv) = abs(vmca(k,iv)) * dtm(itab_w(iwr)%mrlw) &
                        * volti(k,iwr) / rho_old(k,iwr) 
       enddo
    enddo
   !$omp end do

    !$omp do private(iw,k,kr)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W/M levels
       do k = lpw(iw), mza-1
          kr = krecw(k,iw)
          cfl_win(k,iw) = abs(wmca(k,iw)) * dtm(itab_w(iw)%mrlw) &
                        * volti(kr,iw) / rho_old(kr,iw)
       enddo

       cfl_win(mza,iw) = 0.0

       ! Loop over T levels
       do k = lpw(iw), mza
          tfact(k,iw) = rho(k,iw) / rho_old(k,iw)

          ! cfl_out_sum is reused as 1 / cfl_out_sum
          cfl_out_sum(k,iw) = 0.999999 / max( cfl_out_sum(k,iw), 1.e-6)
       enddo

    enddo
    !$omp end do
    !$omp end parallel

  end subroutine comp_cfl2

!===============================================================================

  subroutine comp_vert_limits(mrl, scp, scp_upw, iwdepv, iwrecv, kdepw, krecw)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, itab_v, jtab_w, jtw_prog
    use mem_grid,     only: lpv, lpw, mza, mva, mwa

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upw(mza,mwa)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: iwrecv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)
    integer, intent(in)    :: krecw  (mza,mwa)

    integer :: j, iv, k, kr, kd, iw, iwd, iwr
    real    :: scp_win_min, scp_win_max

    !$omp parallel
    !$omp workshare
    scp_local_min(:,:) = scp(:,:)
    scp_local_max(:,:) = scp(:,:)

    scp_in_max(:,:) = scp(:,:)
    scp_in_min(:,:) = scp(:,:)

    c_scp_in_max_sum(:,:) = 0.0
    c_scp_in_min_sum(:,:) = 0.0
    !$omp end workshare

! Expand inflow bounds of each cell based on the upstream neighbors
! Loop over all immediate V neighbors of each primary W/T column:

    !$omp do private(iv,k,iwd,iwr)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Loop over T levels
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          iwr = iwrecv(k,iv)
          scp_in_max(k,iwr) = max(scp_in_max(k,iwr), scp(k,iwd))
          scp_in_min(k,iwr) = min(scp_in_min(k,iwr), scp(k,iwd))
       enddo
    enddo
    !$omp end do

    !$omp do private(iw,k,kr,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W levels
       do k = lpw(iw), mza-1
          kr = krecw(k,iw)
          kd = kdepw(k,iw)

          scp_win_max = max( scp_in_max(kd,iw), scp(kr,iw) )
          scp_win_min = min( scp_in_min(kd,iw), scp(kr,iw) )

          scp_upw(k,iw) = max( scp_win_min, min(scp_upw(k,iw), scp_win_max) )

! Compute the W contribution to the max/min allowed values in each grid box, and
! the sum of the max and min inputs to each grid box

          scp_local_min(kr,iw) = min( scp_local_min(kr,iw), scp_win_min )
          scp_local_max(kr,iw) = max( scp_local_max(kr,iw), scp_win_max )

          c_scp_in_max_sum(kr,iw) = c_scp_in_max_sum(kr,iw) + cfl_win(k,iw) * &
                                    max( scp(kd,iw), scp_upw(k,iw) )

          c_scp_in_min_sum(kr,iw) = c_scp_in_min_sum(kr,iw) + cfl_win(k,iw) * &
                                    min( scp(kd,iw), scp_upw(k,iw) )
       enddo

    enddo
    !$omp end do
    !$omp end parallel

  end subroutine comp_vert_limits

!===============================================================================

  subroutine comp_horiz_limits(mrl, scp, scp_upv, iwdepv, iwrecv)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, itab_v, jtab_w, jtw_prog
    use olam_mpi_atm, only: mpi_send_w
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use misc_coms,    only: iparallel

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: iwrecv (mza,mva)

    integer :: j, iv, k, kr, kd, iw, iwd, iwr
    real    :: scp_vin_min, scp_vin_max
    integer :: iw1, iw2, iw3, iw4

    !$omp parallel
    !$omp do private(iv,k,iwd,iwr,kd,iw3,iw4,scp_vin_max,scp_vin_min)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Compute inflow bounds at each V interface. Here we include the 
       ! immediate upwind cell and vertically upstream-diagonal neighbor

       iw1 = itab_v(iv)%iw(1)
       iw2 = itab_v(iv)%iw(2)

       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          iwr = iwrecv(k,iv)
          kd  = kdepv(k,iv)

          iw3 = merge( itab_v(iv)%iw(3), iw1, k >= lpw(itab_v(iv)%iw(3)) )
          iw4 = merge( itab_v(iv)%iw(4), iw1, k >= lpw(itab_v(iv)%iw(4)) )

          scp_vin_max = max( scp(k,iw1), scp(k,iw2), scp(k,iw3), scp(k,iw4), scp(kd,iwd) )
          scp_vin_min = min( scp(k,iw1), scp(k,iw2), scp(k,iw3), scp(k,iw4), scp(kd,iwd) )

          scp_upv(k,iv) = max( scp_vin_min, min(scp_upv(k,iv), scp_vin_max) )

          scp_local_max(k,iwr) = max(scp_local_max(k,iwr), scp_vin_max)
          scp_local_min(k,iwr) = min(scp_local_min(k,iwr), scp_vin_min)

          c_scp_in_max_sum(k,iwr) = c_scp_in_max_sum(k,iwr) + cfl_vin(k,iv) * &
                                    max( scp(k,iwd), scp_upv(k,iv) )

          c_scp_in_min_sum(k,iwr) = c_scp_in_min_sum(k,iwr) + cfl_vin(k,iv) * &
                                    min( scp(k,iwd), scp_upv(k,iv) )
       enddo

    enddo
    !$omp end do

    ! Compute scalar out min,max for each cell
    ! Horizontal loop over all primary W/T columns

    !$omp do private(iw,k)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Vertical loop over T levels
       do k = lpw(iw), mza

          scp_out_min(k,iw) = (scp(k,iw) + c_scp_in_max_sum(k,iw) - &
               scp_local_max(k,iw) * tfact(k,iw)) * cfl_out_sum(k,iw)

          scp_out_max(k,iw) = (scp(k,iw) + c_scp_in_min_sum(k,iw) - &
               scp_local_min(k,iw) * tfact(k,iw)) * cfl_out_sum(k,iw)
       enddo

    enddo
    !$omp end do
    !$omp end parallel

    ! MPI send of scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=scp_out_min, rvara2=scp_out_max)
    endif

  end subroutine comp_horiz_limits

!===============================================================================

  subroutine apply_flux_limiters(mrl, kdepw, iwdepv, scp_upw, scp_upv)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, itab_v, jtab_w, jtw_prog
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use olam_mpi_atm, only: mpi_recv_w
    use misc_coms,    only: iparallel
    use obnd,         only: lbcopy_w

    implicit none

    integer, intent(in)    :: mrl
    integer, intent(in)    :: kdepw  (mza,mwa)
    integer, intent(in)    :: iwdepv (mza,mva)
    real,    intent(inout) :: scp_upw(mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)

  integer :: j, iv, k, kd, iw, iwd
    
    ! Limit the vertical fluxes based on the computed scalar outgoing max/min values
    ! Can be done before outgoing max/mins are received at the boundary cells

    ! Horizontal loop over all primary W/T columns
    !$omp parallel do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Vertical loop over W levels
       do k = lpw(iw), mza-1
          kd = kdepw(k,iw)
          scp_upw(k,iw) = max( scp_out_min(kd,iw), min(scp_upw(k,iw), scp_out_max(kd,iw)) )
       enddo

    enddo
    !$omp end parallel do

    ! MPI receive of scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=scp_out_min, rvara2=scp_out_max)
    endif
    call lbcopy_w(mrl, a1=scp_out_min, a2=scp_out_max)

    ! Limit the horizontal fluxes based on the computed scalar outgoing max/min

    ! Horizontal loop over all V points bordering primary W/T columns
    !$omp parallel do private(iv,k,iwd) 
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Vertical loop over T levels
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          scp_upv(k,iv) = max( scp_out_min(k,iwd), min(scp_upv(k,iv), scp_out_max(k,iwd)) )
       enddo

    enddo
    !$omp end parallel do

  end subroutine apply_flux_limiters

End Module mem_thuburn
