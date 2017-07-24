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

  real,    allocatable          :: scp_min(:,:)
  real,    allocatable          :: scp_max(:,:)

  real,    allocatable          :: cfl_out_sum(:,:)
  real,    allocatable, private :: cfl_vin(:,:,:)

  real,    allocatable, private :: scp_out_min(:,:)
  real,    allocatable, private :: scp_out_max(:,:)

  real,    allocatable, private :: tfact(:,:)

  real,    allocatable, private :: cfl_win_t(:,:)
  real,    allocatable, private :: cfl_win_b(:,:)

Contains

!===============================================================================

  subroutine alloc_thuburn(imonot, mza, mwa)
    use misc_coms, only: rinit
    implicit none
    
    integer, intent(in) :: imonot, mza, mwa

    allocate( cfl_out_sum(mza,mwa)) ; cfl_out_sum = rinit
    
    if (imonot == 1) then
       allocate( scp_min    (mza,mwa)) ; scp_min     = rinit
       allocate( scp_max    (mza,mwa)) ; scp_max     = rinit
       allocate( scp_out_min(mza,mwa)) ; scp_out_min = rinit
       allocate( tfact      (mza,mwa)) ; tfact       = rinit
    endif

    if (imonot > 0) then
       allocate( scp_out_max(mza,mwa)) ; scp_out_max = rinit
       allocate( cfl_win_t  (mza,mwa)) ; cfl_win_t   = rinit
       allocate( cfl_win_b  (mza,mwa)) ; cfl_win_b   = rinit
       allocate( cfl_vin    (mza,7,mwa)) ; cfl_vin   = rinit
    endif

  end subroutine alloc_thuburn

!===============================================================================

  subroutine dealloc_thuburn()
    implicit none

    if (allocated( cfl_out_sum)) deallocate( cfl_out_sum )
    if (allocated( scp_min    )) deallocate( scp_min )
    if (allocated( scp_max    )) deallocate( scp_max )
    if (allocated( scp_out_min)) deallocate( scp_out_min )
    if (allocated( scp_out_max)) deallocate( scp_out_max )
    if (allocated( cfl_vin    )) deallocate( cfl_vin )
    if (allocated( cfl_win_t  )) deallocate( cfl_win_t )
    if (allocated( cfl_win_b  )) deallocate( cfl_win_b )
    if (allocated( tfact      )) deallocate( tfact )

  end subroutine dealloc_thuburn

!===============================================================================

  subroutine comp_cfl1(mrl, dtm, vmsca, wmsca, vs, ws, rho_old, imonot, do_check)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, glatw, glonw, mza, mva, mwa
    use consts_coms, only: r8
    use misc_coms,   only: iparallel
    use mem_para,    only: myrank, mgroupsize
    use max_dims,    only: maxgrds
    use mem_basic,   only: wc

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer,           intent(in) :: mrl
    real,              intent(in) :: vmsca   (mza,mva)
    real,              intent(in) :: wmsca   (mza,mwa)
    real,              intent(in) :: vs     (mza,mva)
    real,              intent(in) :: ws     (mza,mwa)
    real(r8),          intent(in) :: dtm    (maxgrds)
    real(r8),          intent(in) :: rho_old(mza,mwa)
    integer,           intent(in) :: imonot
    logical, optional, intent(in) :: do_check

    integer :: j, jv, iv, k, iw, npoly, n
    real    :: cfl_max, cfl_maxs(mgroupsize)
    real    :: w_max, w_maxs(mgroupsize)
    integer :: ier, inode, imax(3), imaxs(3,mgroupsize), iwmax(3), iwmaxs(3,mgroupsize), iwnode
    logical :: check

    check = .false. 
    if (present(do_check)) check = do_check

    !$omp parallel do private(iw,k,jv,iv)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       do k = lpw(iw), mza
          cfl_out_sum(k,iw) = max(wmsca(k,iw),0.0) - min(wmsca(k-1,iw),0.0)
       enddo
       
       do jv = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(jv)

          do k = lpv(iv), mza
             cfl_out_sum(k,iw) = cfl_out_sum(k,iw) - min(itab_w(iw)%dirv(jv) * vmsca(k,iv), 0.0)
          enddo
       enddo

       do k = lpw(iw), mza
          cfl_out_sum(k,iw) = cfl_out_sum(k,iw) &
                            * dtm(itab_w(iw)%mrlw) * volti(k,iw) / rho_old(k,iw)
       enddo

    enddo
    !$omp end parallel do

    if (check) then

       ! Find the max CFL number on each node, and print an error message
       ! if the CFL number is greater than 1

       cfl_max = -1.0
       w_max   =  0.0

       do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
          do k = lpw(iw), mza

             if ( cfl_out_sum(k,iw) > cfl_max  .or.  &
                  cfl_out_sum(k,iw) /= cfl_out_sum(k,iw) ) then
                cfl_max = cfl_out_sum(k,iw)
                imax    = (/ k, iw, itab_w(iw)%iwglobe /)
             endif

             if ( abs(wc(k,iw)) > abs(w_max)  .or.  &
                  wc(k,iw) /= wc(k,iw) ) then
                w_max = wc(k,iw)
                iwmax    = (/ k, iw, itab_w(iw)%iwglobe /)
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

#ifdef OLAM_MPI
       if (iparallel == 1) then
          call MPI_Gather(cfl_max, 1, MPI_REAL, cfl_maxs, 1, MPI_REAL, &
                          0, MPI_COMM_WORLD, ier)
          call MPI_Gather(imax, 3, MPI_INTEGER, imaxs, 3, MPI_INTEGER, &
                          0, MPI_COMM_WORLD, ier)

          call MPI_Gather(w_max, 1, MPI_REAL, w_maxs, 1, MPI_REAL, &
                          0, MPI_COMM_WORLD, ier)
          call MPI_Gather(iwmax, 3, MPI_INTEGER, iwmaxs, 3, MPI_INTEGER, &
                          0, MPI_COMM_WORLD, ier)
       endif
#endif

       if (myrank == 0) then
          inode  = 1
          iwnode = 1
          if (iparallel == 1) then

             if (any( cfl_maxs(:) /= cfl_maxs(:) )) then
                do n = 1, mgroupsize
                   if (cfl_maxs(n) /= cfl_maxs(n)) then
                      inode   = n
                      cfl_max = cfl_maxs(inode)
                      imax(:) = imaxs(:,inode)
                      exit
                   endif
                enddo
             else
                inode   = maxloc(cfl_maxs, dim=1)
                cfl_max = cfl_maxs(inode)
                imax(:) = imaxs(:,inode)
             endif

             if (any( w_maxs(:) /= w_maxs(:) )) then
                do n = 1, mgroupsize
                   if (w_maxs(n) /= w_maxs(n)) then
                      iwnode   = n
                      w_max    = w_maxs(iwnode)
                      iwmax(:) = iwmaxs(:,iwnode)
                      exit
                   endif
                enddo
             else
                iwnode   = maxloc(abs(w_maxs), dim=1)
                w_max    = w_maxs(iwnode)
                iwmax(:) = iwmaxs(:,iwnode)
             endif
          endif

          write(*,'(5x,A,f0.3,3(A,I0))') "Max CFL = ", cfl_max,  &
               " at node ", inode-1, ", iwglobe=", imax(3), ", k=", imax(1)

          write(*,'(5x,A,f0.3,3(A,I0))') "Max  W  = ", w_max,  &
               " at node ", iwnode-1, ", iwglobe=", iwmax(3), ", k=", iwmax(1)
       endif
    endif

    if (imonot > 0) then

       ! cfl_out_sum is reused as 1 / cfl_out_sum

       !$omp parallel do private(iw,k)
       do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
          do k = lpw(iw), mza
             cfl_out_sum(k,iw) = 0.999998 / max( cfl_out_sum(k,iw), 1.e-7)
          enddo
       enddo
       !$omp end parallel do
    
    endif

  end subroutine comp_cfl1

!===============================================================================

  subroutine comp_cfl2(mrl, dtm, vmca, wmca, rho_old, rho, imonot)

    ! Diagnose inflow CFL numbers

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use consts_coms, only: r8
    use max_dims,    only: maxgrds
    use misc_coms,   only: iparallel
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w
    use obnd,        only: lbcopy_w

    implicit none

    integer,  intent(in) :: mrl
    real,     intent(in) :: vmca   (mza,mva)
    real,     intent(in) :: wmca   (mza,mwa)
    real(r8), intent(in) :: dtm    (maxgrds)
    real(r8), intent(in) :: rho_old(mza,mwa)
    real(r8), intent(in) :: rho    (mza,mwa)
    integer,  intent(in) :: imonot

    integer              :: j, iw, k, jv, iv
    real                 :: fact(mza)

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=cfl_out_sum)
    endif

    !$omp parallel do private(iw,k,fact,jv,iv)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W/M levels

!      do k = lpw(iw)+1, mza-1
       do k = lpw(iw), mza
          fact(k)         = dtm(itab_w(iw)%mrlw) * volti(k,iw) / rho_old(k,iw)

          cfl_win_b(k,iw) =  max(wmca(k-1,iw),0.0) * fact(k)
          cfl_win_t(k,iw) = -min(wmca(k,  iw),0.0) * fact(k)
       enddo

       if (imonot == 1) then
          do k = lpw(iw), mza
             tfact(k,iw) = rho(k,iw) / rho_old(k,iw)
          enddo
       endif

       do jv = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(jv)
          do k = lpv(iv), mza
             cfl_vin(k,jv,iw) = max(itab_w(iw)%dirv(jv) * vmca(k,iv), 0.0) * fact(k)
          enddo
       enddo

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=cfl_out_sum)
    endif
    call lbcopy_w(mrl, a1=cfl_out_sum)

  end subroutine comp_cfl2

!===============================================================================

  subroutine comp_and_apply_monot_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w, itab_v
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
    use misc_coms,    only: iparallel
    use obnd,         only: lbcopy_w

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upw(mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)

    integer :: j, iw, ka, k, jv, iv, iwn, kd, iwd, iw1, iw2, iw3, iw4
    real    :: c_scp_in_max_sum(mza), c_scp_in_min_sum(mza)
    real    :: scp_int, scp_inb, scpup
    real    :: scp_upwmax(mza), scp_upwmin(mza)
    real    :: scpmin(mza), scpmax(mza), smin, smax

    !$omp parallel 
    !$omp do private(iv,iw1,iw2,iw3,iw4,k,scpmin,scpmax)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       iw1 = itab_v(iv)%iw(1)
       iw2 = itab_v(iv)%iw(2)
       iw3 = itab_v(iv)%iw(3)
       iw4 = itab_v(iv)%iw(4)

       ! Vertical loop over T levels
       do k = lpv(iv), mza-1
          scpmin(k) = min( scp(k-1,iw1), scp(k,iw1), scp(k+1,iw1),&
                           scp(k-1,iw2), scp(k,iw2), scp(k+1,iw2) )

          scpmax(k) = max( scp(k-1,iw1), scp(k,iw1), scp(k+1,iw1),&
                           scp(k-1,iw2), scp(k,iw2), scp(k+1,iw2) )
       enddo

       do k = max(lpv(iv), lpw(iw3)), mza-1
          scpmin(k) = min(scpmin(k), scp(k-1,iw3), scp(k,iw3), scp(k+1,iw3) )
          scpmax(k) = max(scpmax(k), scp(k-1,iw3), scp(k,iw3), scp(k+1,iw3) )
       enddo

       do k = max(lpv(iv), lpw(iw4)), mza-1
          scpmin(k) = min(scpmin(k), scp(k-1,iw4), scp(k,iw4), scp(k+1,iw4) )
          scpmax(k) = max(scpmax(k), scp(k-1,iw4), scp(k,iw4), scp(k+1,iw4) )
       enddo

       scpmin(mza) = min( scp(mza-1,iw1), scp(mza,iw1), scp(mza-1,iw2), scp(mza,iw2),&
                          scp(mza-1,iw3), scp(mza,iw3), scp(mza-1,iw4), scp(mza,iw4) )
       scpmax(mza) = max( scp(mza-1,iw1), scp(mza,iw1), scp(mza-1,iw2), scp(mza,iw2),&
                          scp(mza-1,iw3), scp(mza,iw3), scp(mza-1,iw4), scp(mza,iw4) )

       do k = lpv(iv), mza
          scp_upv(k,iv) = max( min(scp_upv(k,iv), scpmax(k)), scpmin(k) )
       enddo

    enddo
    !$omp end do

    !$omp do private(iw,ka,k,scp_upwmax,scp_upwmin,jv,iv,iwn,&
    !$omp            scp_int,scp_inb,c_scp_in_max_sum,c_scp_in_min_sum,&
    !$omp            scpup,smax,smin)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka, mza-1
          scp_upwmax(k) = max( scp(k,iw), scp(k+1,iw) )
          scp_upwmin(k) = min( scp(k,iw), scp(k+1,iw) )
       enddo

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          do k = lpv(iv), mza-1
             scp_upwmax(k) = max(scp_upwmax(k), scp(k,iwn), scp(k+1,iwn))
             scp_upwmin(k) = min(scp_upwmin(k), scp(k,iwn), scp(k+1,iwn))
          enddo
       enddo

       scp_upwmax(ka-1) = scp_upwmax(ka)
       scp_upwmin(ka-1) = scp_upwmin(ka)

       scp_upwmax(mza) = scp_upwmax(mza-1)
       scp_upwmin(mza) = scp_upwmin(mza-1)

       do k = ka, mza-1
          scp_upw(k,iw) = max( min(scp_upw(k,iw), scp_upwmax(k)), scp_upwmin(k) )
       enddo

       scp_int = scp(ka+1,iw) * min(cfl_out_sum(ka,iw),cfl_out_sum(ka+1,iw),1.0)

       c_scp_in_max_sum(ka) = cfl_win_t(ka,iw) * max(scp_int,scp_upw(ka,iw))
       c_scp_in_min_sum(ka) = cfl_win_t(ka,iw) * min(scp_int,scp_upw(ka,iw))

       ! Loop over T levels
       do k = ka+1, mza-1

          scp_int = scp(k+1,iw) * min(cfl_out_sum(k,iw),cfl_out_sum(k+1,iw),1.0)
          scp_inb = scp(k-1,iw) * min(cfl_out_sum(k,iw),cfl_out_sum(k-1,iw),1.0)

          c_scp_in_max_sum(k) = cfl_win_t(k,iw) * max(scp_int,scp_upw(k  ,iw)) &
                              + cfl_win_b(k,iw) * max(scp_inb,scp_upw(k-1,iw))

          c_scp_in_min_sum(k) = cfl_win_t(k,iw) * min(scp_int,scp_upw(k  ,iw)) &
                              + cfl_win_b(k,iw) * min(scp_inb,scp_upw(k-1,iw))
       enddo

       scp_inb = scp(mza-1,iw) * min(cfl_out_sum(mza,iw),cfl_out_sum(mza-1,iw),1.0)

       c_scp_in_max_sum(mza) = cfl_win_b(mza,iw) * max(scp_inb,scp_upw(mza-1,iw))
       c_scp_in_min_sum(mza) = cfl_win_b(mza,iw) * min(scp_inb,scp_upw(mza-1,iw))

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          ! Loop over T levels
          do k = lpv(iv), mza

             scpup = scp(k,iwn) * min(cfl_out_sum(k,iwn),cfl_out_sum(k,iw),1.0)

             c_scp_in_max_sum(k) = c_scp_in_max_sum(k) + cfl_vin(k,jv,iw) * &
                                   max( scpup, scp_upv(k,iv) )

             c_scp_in_min_sum(k) = c_scp_in_min_sum(k) + cfl_vin(k,jv,iw) * &
                                   min( scpup, scp_upv(k,iv) )
          enddo
       enddo

       ! Loop over T levels
       do k = ka, mza

          smax = max( scp_upwmax(k-1), scp_upwmax(k) )
          smin = min( scp_upwmin(k-1), scp_upwmin(k) )

          scp_out_min(k,iw) = (scp(k,iw) + c_scp_in_max_sum(k) - &
               smax * tfact(k,iw)) * cfl_out_sum(k,iw)

          scp_out_max(k,iw) = (scp(k,iw) + c_scp_in_min_sum(k) - &
               smin * tfact(k,iw)) * cfl_out_sum(k,iw)

       enddo

    enddo
    !$omp end do
    !$omp end parallel

    ! MPI send of allowed scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=scp_out_min, rvara2=scp_out_max)
    endif
    
    ! Limit the vertical fluxes based on the computed scalar outgoing max/min values
    ! Can be done before outgoing max/mins are received at the boundary cells

    !$omp parallel do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Vertical loop over W levels
       do k = lpw(iw), mza-1
          kd = kdepw(k,iw)
          scp_upw(k,iw) = min( scp_out_max(kd,iw), max(scp_upw(k,iw), scp_out_min(kd,iw)) )
       enddo

    enddo
    !$omp end parallel do

    ! MPI receive of allowed scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=scp_out_min, rvara2=scp_out_max)
    endif
    call lbcopy_w(mrl, a1=scp_out_min, a2=scp_out_max)

    ! Limit the horizontal fluxes based on the computed scalar outgoing max/min

    !$omp parallel do private(iv,k,iwd) 
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Vertical loop over T levels
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          scp_upv(k,iv) = min( scp_out_max(k,iwd), max(scp_upv(k,iv), scp_out_min(k,iwd)) )
       enddo

    enddo
    !$omp end parallel do

  end subroutine comp_and_apply_monot_limits

!===============================================================================

  subroutine comp_and_apply_pd_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
    use misc_coms,    only: iparallel
    use obnd,         only: lbcopy_w
    use var_tables
    use mem_para

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upw(mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)

    integer :: j, iw, ka, k, jv, iv, iwn, kd, iwd
    real    :: c_scp_in_min_sum(mza)
    real    :: scpup, scp_int, scp_inb

    !$omp parallel
    !$omp do private(iv,k)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       do k = lpv(iv), mza
          scp_upv(k,iv) = max( scp_upv(k,iv), 0.0 )
       enddo

    enddo
    !$omp end do

    !$omp do private(iw,ka,k,c_scp_in_min_sum,scp_int,scp_inb,&
    !$omp            jv,iv,iwn,scpup)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka, mza-1
          scp_upw(k,iw) = max(scp_upw(k,iw), 0.0)
       enddo

       scp_int = scp(ka+1,iw) * min(cfl_out_sum(ka,iw),cfl_out_sum(ka+1,iw),1.0)
       c_scp_in_min_sum(ka) = cfl_win_t(ka,iw) * min(scp_int,scp_upw(ka,iw))

       ! Loop over T levels
       do k = ka+1, mza-1

          scp_int = scp(k+1,iw) * min(cfl_out_sum(k,iw),cfl_out_sum(k+1,iw),1.0)
          scp_inb = scp(k-1,iw) * min(cfl_out_sum(k,iw),cfl_out_sum(k-1,iw),1.0)

          c_scp_in_min_sum(k) = cfl_win_t(k,iw) * min(scp_int,scp_upw(k  ,iw)) &
                              + cfl_win_b(k,iw) * min(scp_inb,scp_upw(k-1,iw))
       enddo

       scp_inb = scp(mza-1,iw) * min(cfl_out_sum(mza,iw),cfl_out_sum(mza-1,iw),1.0)
       c_scp_in_min_sum(mza) = cfl_win_b(mza,iw) * min(scp_inb,scp_upw(mza-1,iw))

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          ! Loop over T levels
          do k = lpv(iv), mza

             scpup = scp(k,iwn) * min(cfl_out_sum(k,iwn), cfl_out_sum(k,iw))

             c_scp_in_min_sum(k) = c_scp_in_min_sum(k) + cfl_vin(k,jv,iw) * &
                  min( scpup, scp_upv(k,iv) )
          enddo
       enddo

       ! Loop over T levels
       do k = ka, mza
          scp_out_max(k,iw) = (scp(k,iw) + c_scp_in_min_sum(k)) * cfl_out_sum(k,iw)
       enddo

    enddo
    !$omp end do
    !$omp end parallel

    ! MPI send of allowed scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=scp_out_max)
    endif
    
    ! Limit the vertical fluxes based on the computed scalar outgoing max/min values
    ! Can be done before outgoing max/mins are received at the boundary cells

    !$omp parallel do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Vertical loop over W levels
       do k = lpw(iw), mza-1
          kd = kdepw(k,iw)
          scp_upw(k,iw) = min(scp_upw(k,iw), scp_out_max(kd,iw))
       enddo

    enddo
    !$omp end parallel do

    ! MPI receive of allowed scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=scp_out_max)
    endif
    call lbcopy_w(mrl, a1=scp_out_max)

    ! Limit the horizontal fluxes based on the computed scalar outgoing max/min

    !$omp parallel do private(iv,k,iwd) 
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Vertical loop over T levels
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          scp_upv(k,iv) = min(scp_upv(k,iv), scp_out_max(k,iwd))
       enddo

    enddo
    !$omp end parallel do

  end subroutine comp_and_apply_pd_limits

!===============================================================================

  subroutine comp_cfl2_t(mrl, dtm, vmca, wmca, rho, imonot)

    ! Diagnose inflow CFL numbers

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use consts_coms, only: r8
    use max_dims,    only: maxgrds
    use misc_coms,   only: iparallel
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w
    use obnd,        only: lbcopy_w

    implicit none

    integer,  intent(in) :: mrl
    real,     intent(in) :: vmca(mza,mva)
    real,     intent(in) :: wmca(mza,mwa)
    real(r8), intent(in) :: dtm (maxgrds)
    real(r8), intent(in) :: rho (mza,mwa)
    integer,  intent(in) :: imonot

    integer              :: j, iw, k, jv, iv
    real                 :: fact(mza)
    real(r8)             :: flux(mza)

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=cfl_out_sum)
    endif

    !$omp parallel do private(iw,k,fact,jv,iv,flux)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W/M levels
       do k = lpw(iw), mza
          fact(k)         = dtm(itab_w(iw)%mrlw) * volti(k,iw) / rho(k,iw)

          cfl_win_t(k,iw) = -min(wmca(k,  iw),0.0) * fact(k)
          cfl_win_b(k,iw) =  max(wmca(k-1,iw),0.0) * fact(k)
       enddo

       if (imonot == 1) then
          flux(lpw(iw):mza) = 0.0
       endif

       do jv = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(jv)

          do k = lpv(iv), mza
             cfl_vin(k,jv,iw) = max(itab_w(iw)%dirv(jv) * vmca(k,iv), 0.0) * fact(k)
             if (imonot == 1) flux(k) = flux(k) + itab_w(iw)%dirv(jv) * vmca(k,iv)
          enddo
       enddo

       if (imonot == 1) then
          do k = lpw(iw), mza
             tfact(k,iw) = 1.0 + fact(k) * (flux(k) + wmca(k-1,iw) - wmca(k,iw))
          enddo
       endif

    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=cfl_out_sum)
    endif
    call lbcopy_w(mrl, a1=cfl_out_sum)

  end subroutine comp_cfl2_t

End Module mem_thuburn
