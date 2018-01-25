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

  real,    allocatable          :: cfl_out_sum(:,:)
  real,    allocatable, private :: cfl_vin(:,:,:)

  real,    allocatable, private :: scp_out_min(:,:)
  real,    allocatable, private :: scp_out_max(:,:)

  real,    allocatable, private :: tfact(:,:)

  real,    allocatable, private :: cfl_win_t(:,:)
  real,    allocatable, private :: cfl_win_b(:,:)

  logical, allocatable, private :: isstab(:)

Contains

!===============================================================================

  subroutine alloc_thuburn(imonot1, imonot2, mza, mwa)
    use misc_coms, only: rinit
    implicit none
    

    integer, intent(in) :: imonot1, imonot2, mza, mwa

    allocate( cfl_out_sum(mza,mwa)) ; cfl_out_sum = rinit
    
    if (imonot1 == 1 .or. imonot2 == 1) then
       allocate( scp_out_min(mza,mwa)) ; scp_out_min = rinit
       allocate( tfact      (mza,mwa)) ; tfact       = rinit
       allocate( isstab         (mwa)) ; isstab      = .true.
    endif

    if (imonot1 > 0 .or. imonot2 > 0) then
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
    if (allocated( scp_out_min)) deallocate( scp_out_min )
    if (allocated( scp_out_max)) deallocate( scp_out_max )
    if (allocated( cfl_vin    )) deallocate( cfl_vin )
    if (allocated( cfl_win_t  )) deallocate( cfl_win_t )
    if (allocated( cfl_win_b  )) deallocate( cfl_win_b )
    if (allocated( tfact      )) deallocate( tfact )

  end subroutine dealloc_thuburn

!===============================================================================

  subroutine comp_cfl1(mrl, dtm, vmsca, wmsca, vs, ws, rho, imonot, do_check)

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
    real,              intent(in) :: vmsca(mza,mva)
    real,              intent(in) :: wmsca(mza,mwa)
    real,              intent(in) :: vs   (mza,mva)
    real,              intent(in) :: ws   (mza,mwa)
    real(r8),          intent(in) :: dtm  (maxgrds)
    real(r8),          intent(in) :: rho  (mza,mwa)
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
                            * dtm(itab_w(iw)%mrlw) * volti(k,iw) / rho(k,iw)
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
             if (imonot == 1) tfact(k,iw) = cfl_out_sum(k,iw)
             cfl_out_sum(k,iw) = 1.0 / max( cfl_out_sum(k,iw), 1.e-7)
          enddo
       enddo
       !$omp end parallel do

    endif

  end subroutine comp_cfl1

!===============================================================================

  subroutine comp_cfl2(mrl, dtm, vmca, wmca, rho, imonot)

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
    real                 :: flux(mza)

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=cfl_out_sum)
    endif

    !$omp parallel private(fact,flux)
    !$omp do private(iw,k,jv,iv)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W/M levels
       do k = lpw(iw), mza
          fact     (k)    = dtm(itab_w(iw)%mrlw) * volti(k,iw) / rho(k,iw)
          cfl_win_t(k,iw) = -min(wmca(k,  iw),0.0) * fact(k)
          cfl_win_b(k,iw) =  max(wmca(k-1,iw),0.0) * fact(k)

          if (imonot == 1) flux(k) = cfl_win_t(k,iw) + cfl_win_b(k,iw)
       enddo

       do jv = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(jv)

          do k = lpv(iv), mza
             cfl_vin(k,jv,iw) = max(itab_w(iw)%dirv(jv) * vmca(k,iv), 0.0) * fact(k)
             if (imonot == 1) flux(k) = flux(k) + cfl_vin(k,jv,iw)
          enddo
       enddo

       if (imonot == 1) then
          do k = lpw(iw), mza
             tfact(k,iw) = 1.0 + flux(k) - tfact(k,iw)
          enddo
       endif

    enddo
    !$omp end do
    !$omp end parallel

    if (iparallel == 1) then
       call mpi_recv_w(mrl, rvara1=cfl_out_sum)
    endif
    call lbcopy_w(mrl, a1=cfl_out_sum)

    if (imonot == 1) then
       do iw = 2, mwa
          if (all(cfl_out_sum(lpw(iw):mza,iw) >= 1.0)) then
             isstab(iw) = .true.
          else
             isstab(iw) = .false.
          endif
       enddo
    endif

  end subroutine comp_cfl2

!===============================================================================

  subroutine comp_and_apply_monot_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w, itab_v
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
    use misc_coms,    only: iparallel
    use obnd,         only: lbcopy_w

    use var_tables

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upw(mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)

    integer :: j, iw, ka, k, jv, iv, iwn, kd, iwd, iw1, iw2, iw3, iw4
    integer :: iv1, iv2, iv3, iv4
    real    :: c_scp_in_max_sum(mza), c_scp_in_min_sum(mza)
    real    :: scp_int, scp_inb, scpup
    real    :: smin, smax
    real    :: scpmin(mza+1), scpmax(mza+1)

    real    :: scpminv(mza,mva), scpmaxv(mza,mva)

    !$omp parallel private(c_scp_in_max_sum,c_scp_in_min_sum,scpmax,scpmin)
    !$omp do private(iv,iw1,iw2,iw3,iw4,iv1,iv2,iv3,iv4,ka,k,smax,smin)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       iw1 = itab_v(iv)%iw(1)
       iw2 = itab_v(iv)%iw(2)
       iw3 = itab_v(iv)%iw(3)
       iw4 = itab_v(iv)%iw(4)

       iv1 = itab_v(iv)%iv(1)
       iv2 = itab_v(iv)%iv(2)
       iv3 = itab_v(iv)%iv(3)
       iv4 = itab_v(iv)%iv(4)

       ka = lpv(iv)

       ! Vertical loop over T levels
       do k = 1, mza
          scpminv(k,iv) = min(scp(k,iw1),scp(k,iw2))
          scpmaxv(k,iv) = max(scp(k,iw1),scp(k,iw2))
       enddo

       if (lpv(iv1) <= ka .or. lpv(iv2) <= ka) then
          scpminv(ka-1,iv) = min(scpminv(ka-1,iv), scp(ka-1,iw3))
          scpmaxv(ka-1,iv) = max(scpmaxv(ka-1,iv), scp(ka-1,iw3))
       endif

       do k = max(lpv(iv), min(lpv(iv1),lpv(iv2))), mza
          scpminv(k,iv) = min(scpminv(k,iv), scp(k,iw3))
          scpmaxv(k,iv) = max(scpmaxv(k,iv), scp(k,iw3))
       enddo

       if (lpv(iv3) <= ka .or. lpv(iv4) <= ka) then
          scpminv(ka-1,iv) = min(scpminv(ka-1,iv), scp(ka-1,iw4))
          scpmaxv(ka-1,iv) = max(scpmaxv(ka-1,iv), scp(ka-1,iw4))
       endif

       do k = max(lpv(iv), min(lpv(iv3),lpv(iv4))), mza
          scpminv(k,iv) = min(scpminv(k,iv), scp(k,iw4))
          scpmaxv(k,iv) = max(scpmaxv(k,iv), scp(k,iw4))
       enddo

       do k = lpv(iv), mza-1
          smax = max(scpmaxv(k-1,iv), scpmaxv(k,iv), scpmaxv(k+1,iv))
          smin = min(scpminv(k-1,iv), scpminv(k,iv), scpminv(k+1,iv))
          scp_upv(k,iv) = max( min(scp_upv(k,iv), smax), smin )
       enddo

       smax = max(scpmaxv(mza-1,iv), scpmaxv(mza,iv))
       smin = min(scpminv(mza-1,iv), scpminv(mza,iv))
       scp_upv(k,iv) = max( min(scp_upv(k,iv), smax), smin )

    enddo
    !$omp end do

    !$omp do private(iw,ka,k,jv,iv,iwn,scp_int,scp_inb,scpup,smax,smin)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka-1, mza
          scpmax(k) = maxval( scpmaxv( k, itab_w(iw)%iv( 1:itab_w(iw)%npoly ) ) )
          scpmin(k) = minval( scpminv( k, itab_w(iw)%iv( 1:itab_w(iw)%npoly ) ) )
       enddo

       do k = ka, mza-1
          smax = max( scpmax(k), scpmax(k+1) )
          smin = min( scpmin(k), scpmin(k+1) )
          scp_upw(k,iw) = max( min(scp_upw(k,iw), smax), smin )
       enddo

       if (isstab(iw)) then

          c_scp_in_max_sum(ka) = cfl_win_t(ka,iw) * max(scp(ka+1,iw),scp_upw(ka,iw))
          c_scp_in_min_sum(ka) = cfl_win_t(ka,iw) * min(scp(ka+1,iw),scp_upw(ka,iw))

          ! Loop over T levels
          do k = ka+1, mza-1
             c_scp_in_max_sum(k) = cfl_win_t(k,iw) * max(scp(k+1,iw),scp_upw(k  ,iw)) &
                                 + cfl_win_b(k,iw) * max(scp(k-1,iw),scp_upw(k-1,iw))

             c_scp_in_min_sum(k) = cfl_win_t(k,iw) * min(scp(k+1,iw),scp_upw(k  ,iw)) &
                                 + cfl_win_b(k,iw) * min(scp(k-1,iw),scp_upw(k-1,iw))
          enddo

          c_scp_in_max_sum(mza) = cfl_win_b(mza,iw) * max(scp(mza-1,iw),scp_upw(mza-1,iw))
          c_scp_in_min_sum(mza) = cfl_win_b(mza,iw) * min(scp(mza-1,iw),scp_upw(mza-1,iw))

       else

          scp_int = max(scp(ka+1,iw) * min(cfl_out_sum(ka+1,iw),1.0), scpmin(ka+1))

          c_scp_in_max_sum(ka) = cfl_win_t(ka,iw) * max(scp(ka+1,iw),scp_upw(ka,iw))
          c_scp_in_min_sum(ka) = cfl_win_t(ka,iw) * min(scp_int,scp_upw(ka,iw))

          ! Loop over T levels
          do k = ka+1, mza-1
             scp_int = max(scp(k+1,iw) * min(cfl_out_sum(k+1,iw),1.0), scpmin(k+1))
             scp_inb = max(scp(k-1,iw) * min(cfl_out_sum(k-1,iw),1.0), scpmin(k-1))

             c_scp_in_max_sum(k) = cfl_win_t(k,iw) * max(scp(k+1,iw),scp_upw(k  ,iw)) &
                                 + cfl_win_b(k,iw) * max(scp(k-1,iw),scp_upw(k-1,iw))

             c_scp_in_min_sum(k) = cfl_win_t(k,iw) * min(scp_int,scp_upw(k  ,iw)) &
                                 + cfl_win_b(k,iw) * min(scp_inb,scp_upw(k-1,iw))
          enddo

          scp_inb = max(scp(mza-1,iw) * min(cfl_out_sum(mza-1,iw),1.0), scpmin(mza-1))

          c_scp_in_max_sum(mza) = cfl_win_b(mza,iw) * max(scp(mza-1,iw),scp_upw(mza-1,iw))
          c_scp_in_min_sum(mza) = cfl_win_b(mza,iw) * min(scp_inb,scp_upw(mza-1,iw))

       endif

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          if (isstab(iwn)) then

             ! Loop over T levels
             do k = lpv(iv), mza
                c_scp_in_max_sum(k) = c_scp_in_max_sum(k) + cfl_vin(k,jv,iw) * &
                                      max( scp(k,iwn), scp_upv(k,iv) )

                c_scp_in_min_sum(k) = c_scp_in_min_sum(k) + cfl_vin(k,jv,iw) * &
                                      min( scp(k,iwn), scp_upv(k,iv) )
             enddo

          else

             ! Loop over T levels
             do k = lpv(iv), mza
                scpup = max(scp(k,iwn) * min(cfl_out_sum(k,iwn),1.0), scpminv(k,iv))

                c_scp_in_max_sum(k) = c_scp_in_max_sum(k) + cfl_vin(k,jv,iw) * &
                                      max( scp(k,iwn), scp_upv(k,iv) )

                c_scp_in_min_sum(k) = c_scp_in_min_sum(k) + cfl_vin(k,jv,iw) * &
                                      min( scpup, scp_upv(k,iv) )
             enddo

          endif

       enddo

       scpmax(mza+1) = scpmax(mza)
       scpmin(mza+1) = scpmin(mza)

       ! Loop over T levels
       do k = ka, mza
          smax = max( scpmax(k-1), scpmax(k), scpmax(k+1) )
          smin = min( scpmin(k-1), scpmin(k), scpmin(k+1) )

          scp_out_min(k,iw) = (scp(k,iw) + c_scp_in_max_sum(k) - &
               smax * tfact(k,iw)) * cfl_out_sum(k,iw)

          scp_out_max(k,iw) = 0.9999*(scp(k,iw) + c_scp_in_min_sum(k) - &
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

    !$omp parallel private(c_scp_in_min_sum)
    !$omp do private(iw,ka,k,scp_int,scp_inb,jv,iv,iwn,scpup)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka, mza-1
          scp_upw(k,iw) = max(scp_upw(k,iw), 0.0)
       enddo

       scp_int              = scp(ka+1,iw) * min(cfl_out_sum(ka+1,iw),1.0)
       c_scp_in_min_sum(ka) = cfl_win_t(ka,iw) * min(scp_int,scp_upw(ka,iw))

       ! Loop over T levels
       do k = ka+1, mza-1

          scp_int = scp(k+1,iw) * min(cfl_out_sum(k+1,iw),1.0)
          scp_inb = scp(k-1,iw) * min(cfl_out_sum(k-1,iw),1.0)

          c_scp_in_min_sum(k) = cfl_win_t(k,iw) * min(scp_int,scp_upw(k  ,iw)) &
                              + cfl_win_b(k,iw) * min(scp_inb,scp_upw(k-1,iw))
       enddo

       scp_inb               = scp(mza-1,iw) * min(cfl_out_sum(mza-1,iw),1.0)
       c_scp_in_min_sum(mza) = cfl_win_b(mza,iw) * min(scp_inb,scp_upw(mza-1,iw))

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          ! Loop over T levels
          do k = lpv(iv), mza
             scpup = scp(k,iwn) * min(cfl_out_sum(k,iwn), 1.0)

             c_scp_in_min_sum(k) = c_scp_in_min_sum(k) + cfl_vin(k,jv,iw) * &
                  max(0., min( scpup, scp_upv(k,iv) ) )
          enddo
       enddo

       ! Loop over T levels
       do k = ka, mza
          scp_out_max(k,iw) = 0.9999 * (scp(k,iw) + c_scp_in_min_sum(k)) * cfl_out_sum(k,iw)
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
          scp_upv(k,iv) = max(0., min(scp_upv(k,iv), scp_out_max(k,iwd)))
       enddo

    enddo
    !$omp end parallel do

  end subroutine comp_and_apply_pd_limits

End Module mem_thuburn
