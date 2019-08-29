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

  use consts_coms, only: r8

  real(r8), allocatable, private :: cfl_out_sum(:,:)
  real,     allocatable, private :: scp_out_min(:,:)
  real,     allocatable, private :: scp_out_max(:,:)
  real(r8), allocatable, private :: tfact      (:,:)
  real,     allocatable, private :: dtom       (:,:)

  integer,  allocatable, private :: kbot12(:)
  integer,  allocatable, private :: kbot34(:)

  private :: r8

Contains

!===============================================================================

  subroutine alloc_thuburn(imonot1, imonot2, imonot3, mza, mwa, mva)

    use misc_coms, only: rinit, rinit8
    implicit none

    integer, intent(in) :: imonot1, imonot2, imonot3, mza, mwa, mva

    if (imonot1 == 1 .or. imonot2 == 1 .or. imonot3 == 1) then
       allocate( scp_out_min(mza,mwa)) ; scp_out_min = rinit
       allocate( tfact      (mza,mwa)) ; tfact       = rinit8
       allocate( kbot12         (mva)) ; kbot12      = 1
       allocate( kbot34         (mva)) ; kbot34      = 1
    endif

    if (imonot1 > 0 .or. imonot2 > 0 .or. imonot3 > 0) then
       allocate( cfl_out_sum(mza,mwa)) ; cfl_out_sum = rinit8
       allocate( scp_out_max(mza,mwa)) ; scp_out_max = rinit
       allocate( dtom       (mza,mwa)) ; dtom        = rinit
    endif

  end subroutine alloc_thuburn

!===============================================================================

  subroutine dealloc_thuburn()
    implicit none

    if (allocated( cfl_out_sum )) deallocate( cfl_out_sum )
    if (allocated( scp_out_min )) deallocate( scp_out_min )
    if (allocated( scp_out_max )) deallocate( scp_out_max )
    if (allocated( tfact       )) deallocate( tfact )
    if (allocated( dtom        )) deallocate( dtom )

  end subroutine dealloc_thuburn

!===============================================================================

  subroutine comp_cfls_long(mrl, vmsca, wmsca, rho, imonot, do_check)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, glatw, glonw, mza, mva, mwa
    use consts_coms, only: r8
    use misc_coms,   only: iparallel, dtlm
    use mem_para,    only: myrank, mgroupsize
    use mem_basic,   only: wc, vc, rho2=>rho
!$  use omp_lib,     only: omp_get_max_threads, omp_get_thread_num
    use mem_nudge,   only: nudflag, rhot_nud
    use oname_coms,  only: nl

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer,           intent(in) :: mrl
    real,              intent(in) :: vmsca(mza,mva)
    real,              intent(in) :: wmsca(mza,mwa)
    real(r8),          intent(in) :: rho  (mza,mwa)
    integer,           intent(in) :: imonot
    logical, optional, intent(in) :: do_check

    integer :: j, iv, k, iw, n, jv
    integer :: ier, inode, iwnode
    logical :: check
    integer :: mlocc, mlocw, nthreads, myid
    real    :: dt
    real(r8):: flux(mza), cfl8
    real    :: cfl(mza), dovr(mza)

    integer, allocatable :: imax(:,:), iwmax(:,:)
    real,    allocatable :: cfl_max(:), w_max(:)

#ifdef OLAM_MPI
    real :: sbuf(6), rbuf(6,mgroupsize)
#endif

    check = .false.
    if (present(do_check)) check = do_check

    if (.not. check .and. imonot < 1) return

    myid     = 1
    nthreads = 1
 !$ nthreads = omp_get_max_threads()

    if (check) then
       allocate(cfl_max(nthreads))
       allocate(imax (2,nthreads))

       allocate(w_max  (nthreads))
       allocate(iwmax(2,nthreads))

       cfl_max = -1.0
       w_max   =  0.0
    endif

    !$omp parallel private(myid,cfl,flux,dovr) shared(cfl_max,w_max,imax,iwmax)
    !$ myid = omp_get_thread_num() + 1

    !$omp do private(iw,dt,k,n,iv,jv,cfl8)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       dt = dtlm(itab_w(iw)%mrlw)

       if (imonot == 2 .or. (imonot == 0 .and. check)) then

          do k = lpw(iw), mza
             cfl (k) = max(wmsca(k,iw),0.0) - min(wmsca(k-1,iw),0.0)
             dovr(k) = dt * volti(k,iw) / real(rho(k,iw))
          enddo

          do n = 1, itab_w(iw)%npoly
             iv = itab_w(iw)%iv(n)
             do k = lpv(iv), mza
                cfl(k) = cfl(k) - min(itab_w(iw)%dirv(n) * vmsca(k,iv), 0.0)
             enddo
          enddo

          do k = lpw(iw), mza
             cfl(k) = cfl(k) * dovr(k)
          enddo

          if (imonot == 2) then
             do k = lpw(iw), mza
                dtom       (k,iw) = dovr(k)
                cfl_out_sum(k,iw) = 1._r8 / real( max( cfl(k), 1.e-7 ), r8)
             enddo
          endif

       elseif (imonot == 1) then

          ! Loop over W/M levels
          do k = lpw(iw), mza
             flux(k)    = max(wmsca(k-1,iw),0.0) - min(wmsca(k,iw),0.0)
             dtom(k,iw) = dt * volti(k,iw) / real(rho(k,iw))
          enddo

          do jv = 1, itab_w(iw)%npoly
             iv = itab_w(iw)%iv(jv)
             do k = lpv(iv), mza
                flux(k) = flux(k) + max(itab_w(iw)%dirv(jv) * vmsca(k,iv), 0.0)
             enddo
          enddo

          if ( nudflag > 0 .and. nl%nud_preserve_mix_ratio .and. &
               itab_w(iw)%mrlw <= nl%max_nud_mrl ) then
             do k = lpw(iw), mza
                tfact(k,iw) = (rho2(k,iw) - dtlm(itab_w(iw)%mrlw) * rhot_nud(k,iw)) / rho(k,iw)
             enddo
          else
             do k = lpw(iw), mza
                tfact(k,iw) = rho2(k,iw) / rho(k,iw)
             enddo
          endif

          do k = lpw(iw), mza
             cfl8   =  max(1._r8 + flux(k) * dtom(k,iw) - tfact(k,iw), 1.e-9_r8)
             if (check) cfl(k) = cfl8
             cfl_out_sum(k,iw) = 1._r8 / cfl8
          enddo

       endif

       if (check) then
          do k = lpw(iw), mza

             if (cfl(k) > cfl_max(myid) .or. cfl(k) /= cfl(k)) then
                cfl_max(myid) = cfl(k)
                imax (:,myid) = (/ k, itab_w(iw)%iwglobe /)
             endif

             if (abs(wc(k,iw)) > abs(w_max(myid)) .or. wc(k,iw) /= wc(k,iw)) then
                w_max  (myid) = wc(k,iw)
                iwmax(:,myid) = (/ k, itab_w(iw)%iwglobe /)
             endif

             if ( cfl(k) > 1.0 .or. cfl(k) /= cfl(k) ) then
                n = itab_w(iw)%npoly
                write(*,*)
                write(*,'(4(A,I0),2(A,f0.3),/,2(A,f0.3),A,7(f0.3,1x))')           &
                     "!!! CFL VIOLATION at node ", myrank,                        &
                     ", iw=", iw, ", iwglobe=", itab_w(iw)%iwglobe, ", k=", k,    &
                     " lat=", glatw(iw), " lon=", glonw(iw),                      &
                     "!!! CFL = ", cfl(k),                             &
                     ", W = ", wc(k,iw), ", VC = ", vc(k,itab_w(iw)%iv(1:n))
             endif

          enddo
       endif

    enddo
    !$omp end do nowait
    !$omp end parallel

    if (check) then

       mlocc = 1
       mlocw = 1

       !$ if (nthreads > 1) then
       !$
       !$    mlocw = maxloc(abs(w_max), dim=1)
       !$    mlocc = maxloc(cfl_max, dim=1)
       !$
       !$    do n = 1, nthreads
       !$       if (cfl_max(n) /= cfl_max(n)) mlocc = n
       !$       if (  w_max(n) /=   w_max(n)) mlocw = n
       !$    enddo
       !$ endif

#ifdef OLAM_MPI
       if (iparallel == 1) then
          sbuf(1)   = cfl_max(mlocc)
          sbuf(2:3) = real(imax(:,mlocc))
          sbuf(4)   = w_max(mlocw)
          sbuf(5:6) = real(iwmax(:,mlocw))

          call MPI_Gather(sbuf, 6, MPI_REAL, rbuf, 6, MPI_REAL, &
                          0, MPI_COMM_WORLD, ier)
       endif
#endif

       if (myrank == 0) then
          inode  = 1
          iwnode = 1

#ifdef OLAM_MPI
          if (iparallel == 1) then

             if (any( rbuf(1,:) /= rbuf(1,:) )) then
                do n = 1, mgroupsize
                   if (rbuf(1,n) /= rbuf(1,n)) then
                      inode          = n
                      cfl_max(mlocc) = rbuf(1,inode)
                      imax (:,mlocc) = nint(rbuf(2:3,inode))
                      exit
                   endif
                enddo
             else
                inode          = maxloc(rbuf(1,:), dim=1)
                cfl_max(mlocc) = rbuf(1,inode)
                imax (:,mlocc) = nint(rbuf(2:3,inode))
             endif

             if (any( rbuf(4,:) /= rbuf(4,:) )) then
                do n = 1, mgroupsize
                   if (rbuf(4,n) /= rbuf(4,n)) then
                      iwnode         = n
                      w_max  (mlocw) = rbuf(4,iwnode)
                      iwmax(:,mlocw) = nint(rbuf(5:6,iwnode))
                      exit
                   endif
                enddo
             else
                iwnode         = maxloc(abs(rbuf(4,:)), dim=1)
                w_max  (mlocw) = rbuf(4,iwnode)
                iwmax(:,mlocw) = nint(rbuf(5:6,iwnode))
             endif
          endif
#endif

          write(*,'(5x,A,f0.3,3(A,I0))') "Max CFL = ", cfl_max(mlocc),  &
               " at node ", inode-1, ", iwglobe=", imax(2,mlocc), ", k=", imax(1,mlocc)

          write(*,'(5x,A,f0.3,3(A,I0))') "Max  W  = ", w_max(mlocw),  &
               " at node ", iwnode-1, ", iwglobe=", iwmax(2,mlocw), ", k=", iwmax(1,mlocw)
       endif
    endif

  end subroutine comp_cfls_long

!===============================================================================

  subroutine comp_cfls_short(mrl, vmsca, wmsca)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use consts_coms, only: r8
    use misc_coms,   only: dtsm
    use mem_basic,   only: rho

    implicit none

    integer,           intent(in) :: mrl
    real,              intent(in) :: vmsca(mza,mva)
    real,              intent(in) :: wmsca(mza,mwa)

    integer :: j, iv, k, iw, jv
    real    :: dt, vmca
    real(r8):: flux(mza), cfl(mza)

    !$omp parallel private(cfl,flux)
    !$omp do private(iw,dt,k,jv,iv,vmca)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       dt = dtsm(itab_w(iw)%mrlw)

       ! Loop over T levels
       do k = lpw(iw), mza
          cfl (k)    = max(wmsca(k,  iw),0.0) - min(wmsca(k-1,iw),0.0)
          flux(k)    = max(wmsca(k-1,iw),0.0) - min(wmsca(k,  iw),0.0)
          dtom(k,iw) = dt * volti(k,iw) / real(rho(k,iw))
       enddo

       ! Loop over adjacent V faces
       do jv = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(jv)

          ! Loop over T levels
          do k = lpv(iv), mza
             vmca    = itab_w(iw)%dirv(jv) * vmsca(k,iv)
             cfl (k) = cfl (k) - min(vmca, 0.0)
             flux(k) = flux(k) + max(vmca, 0.0)
          enddo
       enddo

       ! Loop over T levels
       do k = lpw(iw), mza
          flux       (k)    = flux(k) * dtom(k,iw)
          tfact      (k,iw) = 1._r8 + flux(k) - cfl(k) * dtom(k,iw)
          cfl_out_sum(k,iw) = 1._r8 / max(1._8 + flux(k) - tfact(k,iw), 1.e-12_r8)
       enddo

    enddo
    !$omp end do nowait
    !$omp end parallel

  end subroutine comp_cfls_short

!===============================================================================

  subroutine comp_and_apply_monot_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv, &
                                         wmsca, vmsca)

    use mem_ijtabs,   only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w, &
                            itab_v, jtv_prog
    use mem_grid,     only: lpv, lpw, mza, mva, mwa
    use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
    use misc_coms,    only: iparallel
    use obnd,         only: lbcopy_w, lbcopy_v

    implicit none

    integer, intent(in)    :: mrl
    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upw(mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)

    real, intent(in) :: wmsca(mza,mwa)
    real, intent(in) :: vmsca(mza,mva)

    integer :: j, iw, ka, k, jv, iv, iwn, kd, iwd, iw1, iw2, iw3, iw4
    integer :: iv1, iv2, iv3, iv4
    real(r8):: c_scp_in_max_sum(mza), c_scp_in_min_sum(mza)
    real    :: scp_int, scp_inb, scpup
    real    :: smin, smax, vin, wint, winb
    real    :: scpmin(mza), scpmax(mza)
    real    :: scpmint(mza), scpmaxt(mza)

    logical, save :: firstime = .true.

    if (firstime) then

       !$omp parallel do private(iv, iv1, iv2, iv3, iv4)
       do j = 1, jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)

          iv1 = itab_v(iv)%iv(1)
          iv2 = itab_v(iv)%iv(2)
          iv3 = itab_v(iv)%iv(3)
          iv4 = itab_v(iv)%iv(4)

          kbot12(iv) = max(lpv(iv), min(lpv(iv1), lpv(iv2))) - 1
          kbot34(iv) = max(lpv(iv), min(lpv(iv3), lpv(iv4))) - 1
       enddo
       !$omp end parallel do

       if (iparallel == 1) then
          call mpi_send_v(1, i1dvara1=kbot12, i1dvara2=kbot34)
          call mpi_recv_v(1, i1dvara1=kbot12, i1dvara2=kbot34)
       endif
       call lbcopy_v(1, iv1=kbot12, iv2=kbot34)

    endif

    !$omp parallel private(c_scp_in_max_sum,c_scp_in_min_sum,scpmax,scpmin,&
    !$                     scpmaxt,scpmint)
    !$omp do private(iv,iw1,iw2,iw3,iw4,k,smax,smin)
    do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       iw1 = itab_v(iv)%iw(1)
       iw2 = itab_v(iv)%iw(2)
       iw3 = itab_v(iv)%iw(3)
       iw4 = itab_v(iv)%iw(4)

       ! Vertical loop over T levels
       do k = lpv(iv)-1, mza
          scpmin(k) = min(scp(k,iw1),scp(k,iw2))
          scpmax(k) = max(scp(k,iw1),scp(k,iw2))
       enddo

       do k = kbot12(iv), mza
          scpmin(k) = min(scpmin(k), scp(k,iw3))
          scpmax(k) = max(scpmax(k), scp(k,iw3))
       enddo

       do k = kbot34(iv), mza
          scpmin(k) = min(scpmin(k), scp(k,iw4))
          scpmax(k) = max(scpmax(k), scp(k,iw4))
       enddo

       do k = lpv(iv), mza-1
          smax = max(scpmax(k-1), scpmax(k), scpmax(k+1))
          smin = min(scpmin(k-1), scpmin(k), scpmin(k+1))
          scp_upv(k,iv) = max( min(scp_upv(k,iv), smax), smin )
       enddo

       smax = max(scpmax(mza-1), scpmax(mza))
       smin = min(scpmin(mza-1), scpmin(mza))
       scp_upv(mza,iv) = max( min(scp_upv(mza,iv), smax), smin )

    enddo
    !$omp end do

    !$omp do private(iw,ka,k,jv,iv,iwn,scp_int,scp_inb,scpup,smax,smin,vin,wint,winb)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka, mza
          scpmaxt(k) = scp(k,iw)
          scpmint(k) = scp(k,iw)
       enddo

       c_scp_in_max_sum = 0._r8
       c_scp_in_min_sum = 0._r8

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          do k = lpv(iv), mza
             scpmaxt(k) = max(scpmaxt(k), scp(k,iwn))
             scpmint(k) = min(scpmint(k), scp(k,iwn))

             vin = max(itab_w(iw)%dirv(jv) * vmsca(k,iv), 0.0)

             c_scp_in_max_sum(k) = c_scp_in_max_sum(k)  &
                                 + vin * max( scp(k,iwn), scp_upv(k,iv) )

             c_scp_in_min_sum(k) = c_scp_in_min_sum(k)  &
                                 + vin * min( scp(k,iwn), scp_upv(k,iv) )
          enddo
       enddo

       ! Loop over W levels
       do k = ka, mza-1
          scpmax(k) = max( scpmaxt(k), scpmaxt(k+1) )
          scpmin(k) = min( scpmint(k), scpmint(k+1) )
          scp_upw(k,iw) = max( min(scp_upw(k,iw), scpmax(k)), scpmin(k) )
       enddo

       scpmax(mza) = scpmaxt(mza)
       scpmin(mza) = scpmint(mza)

       scpmax(ka-1) = scpmaxt(ka)
       scpmin(ka-1) = scpmint(ka)

       wint = -min(wmsca(ka,iw),0.0)
       c_scp_in_max_sum(ka) = c_scp_in_max_sum(ka) + wint * max(scp(ka+1,iw),scp_upw(ka,iw))
       c_scp_in_min_sum(ka) = c_scp_in_min_sum(ka) + wint * min(scp(ka+1,iw),scp_upw(ka,iw))

       ! Loop over T levels
       do k = ka+1, mza-1
          wint = min(wmsca(k  ,iw),0.0)
          winb = max(wmsca(k-1,iw),0.0)

          c_scp_in_max_sum(k) = c_scp_in_max_sum(k) &
                              + winb * max(scp(k-1,iw),scp_upw(k-1,iw)) &
                              - wint * max(scp(k+1,iw),scp_upw(k  ,iw))

          c_scp_in_min_sum(k) = c_scp_in_min_sum(k) &
                              + winb * min(scp(k-1,iw),scp_upw(k-1,iw)) &
                              - wint * min(scp(k+1,iw),scp_upw(k  ,iw))
       enddo

       winb = max(wmsca(mza-1,iw),0.0)
       c_scp_in_max_sum(mza) = c_scp_in_max_sum(mza) + winb * max(scp(mza-1,iw),scp_upw(mza-1,iw))
       c_scp_in_min_sum(mza) = c_scp_in_min_sum(mza) + winb * min(scp(mza-1,iw),scp_upw(mza-1,iw))

       ! Loop over T levels
       do k = ka, mza
          smax = max( scpmax(k-1), scpmax(k) )
          smin = min( scpmin(k-1), scpmin(k) )

          scp_out_min(k,iw) = ( scp(k,iw) + c_scp_in_max_sum(k) * dtom(k,iw) - &
                                smax * tfact(k,iw) ) * cfl_out_sum(k,iw)

          scp_out_max(k,iw) = ( scp(k,iw) + c_scp_in_min_sum(k) * dtom(k,iw) - &
                                smin * tfact(k,iw) ) * cfl_out_sum(k,iw)
       enddo

       do k = ka, mza-1
          kd = kdepw(k,iw)
          scp_upw(k,iw) = min( scp_out_max(kd,iw), max(scp_upw(k,iw), scp_out_min(kd,iw)) )
       enddo

    enddo
    !$omp end do
    !$omp end parallel

    ! MPI send of allowed scalar max/min outflow values

    if (iparallel == 1) then
       call mpi_send_w(mrl, rvara1=scp_out_min, rvara2=scp_out_max)
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

  subroutine comp_and_apply_pd_limits(mrl, scp, scp_upw, scp_upv, kdepw, iwdepv, &
                                      wmsca, vmsca)

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

    real, intent(in) :: wmsca(mza,mwa)
    real, intent(in) :: vmsca(mza,mva)

    integer :: j, iw, ka, k, jv, iv, iwn, kd, iwd
    real(r8):: c_scp_in_min_sum(mza)
    real    :: scpup, scp_int, scp_inb

    !$omp parallel private(c_scp_in_min_sum)
    !$omp do private(iw,ka,k,scp_int,scp_inb,jv,iv,iwn,scpup)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ka = lpw(iw)

       do k = ka, mza-1
          scp_upw(k,iw) = max(scp_upw(k,iw), 0.0)
       enddo

       c_scp_in_min_sum(ka) = -min(wmsca(ka,iw),0.0) * min(scp(ka+1,iw),scp_upw(ka,iw))

       ! Loop over T levels
       do k = ka+1, mza-1
          c_scp_in_min_sum(k) = max(wmsca(k-1,iw),0.0) * min(scp(k-1,iw),scp_upw(k-1,iw)) &
                              - min(wmsca(k,  iw),0.0) * min(scp(k+1,iw),scp_upw(k  ,iw))
       enddo

       c_scp_in_min_sum(mza) = max(wmsca(mza-1,iw),0.0) * min(scp(mza-1,iw),scp_upw(mza-1,iw))

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)
          iwn = itab_w(iw)%iw(jv)

          ! Loop over T levels
          do k = lpv(iv), mza
             c_scp_in_min_sum(k) = c_scp_in_min_sum(k) &
                                 + max(itab_w(iw)%dirv(jv) * vmsca(k,iv), 0.0) * &
                                   max(0., min( scp(k,iwn), scp_upv(k,iv) ) )
          enddo
       enddo

       ! Loop over T levels
       do k = ka, mza
          scp_out_max(k,iw) = (scp(k,iw) + c_scp_in_min_sum(k) * dtom(k,iw)) * cfl_out_sum(k,iw)
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
