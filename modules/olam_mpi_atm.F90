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

Module olam_mpi_atm

Contains

  subroutine olam_mpi_init()

#ifdef OLAM_MPI
    use mpi
#endif

    use mem_para, only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8
    implicit none

#ifdef OLAM_MPI

    integer :: ierr, iomp
!$  integer :: iprov

! Initialize MPI and determine process groupsize and myrank

    iomp = 0

!$  iomp = 1
!$  call MPI_Init_thread(MPI_THREAD_MULTIPLE, iprov, ierr)
!$  if (iprov < MPI_THREAD_MULTIPLE) then
!$     write(*,*) "MPI_THREAD_MULTIPLE is not provided by MPI"
!$     call olam_stop("Stopping model run")
!$  endif

    if (iomp == 0) then
       call MPI_Init(ierr)
    endif

    call MPI_Comm_size(MPI_COMM_WORLD,mgroupsize,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

    call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, nbytes_int  , ierr)
    call MPI_Pack_size(1, MPI_REAL   , MPI_COMM_WORLD, nbytes_real , ierr)
    call MPI_Pack_size(1, MPI_REAL8  , MPI_COMM_WORLD, nbytes_real8, ierr)

#else

    mgroupsize = 1
    myrank     = 0

    nbytes_int   = 4
    nbytes_real  = 4
    nbytes_real8 = 8

#endif

end subroutine olam_mpi_init

!===============================================================================

subroutine alloc_mpi_sndrcv_bufs()

  use mem_para,  only: mgroupsize,           &
                       send_v,    recv_v,    &
                       send_w,    recv_w,    &
                       send_m,    recv_m,    &
                       send_wnud, recv_wnud, &
                       send_wsfc, recv_wsfc

  use misc_coms, only: iparallel
  use mem_nudge, only: nudflag, nudnxp

  implicit none

  if (iparallel == 0) return

! Allocate send and recv tables

  allocate(send_v(mgroupsize))
  allocate(recv_v(mgroupsize))

  allocate(send_m(mgroupsize))
  allocate(recv_m(mgroupsize))

  allocate(send_w(mgroupsize))
  allocate(recv_w(mgroupsize))

  if (nudflag > 0 .and. nudnxp > 0) then
     allocate(send_wnud(mgroupsize))
     allocate(recv_wnud(mgroupsize))
  endif

  allocate(send_wsfc(mgroupsize))
  allocate(recv_wsfc(mgroupsize))

end subroutine alloc_mpi_sndrcv_bufs

!===============================================================================

subroutine olam_mpi_finalize()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,  only: send_v,    recv_v,    ireqr_v,    inext_v,    nrecvs_v,    &
                       send_w,    recv_w,    ireqr_w,    inext_w,    nrecvs_m,    &
                       send_m,    recv_m,    ireqr_m,    inext_m,    nrecvs_w,    &
                       send_wnud, recv_wnud, ireqr_wnud, inext_wnud, nrecvs_wnud, &
                       send_wsfc, recv_wsfc, ireqr_wsfc, inext_wsfc, nrecvs_wsfc
  use mem_nudge, only: nudflag, nudnxp
  use misc_coms, only: iparallel

  implicit none

#ifdef OLAM_MPI

  integer :: ierr, jrecv

  if (iparallel == 1) then

     ! First make sure all processes get here so that all sends are finished

     call olam_mpi_barrier()

     ! Cancel pending receives

     do jrecv = 1, nrecvs_v(1)
        if (ireqr_v(jrecv,inext_v) /= MPI_REQUEST_NULL) then
           call MPI_Cancel(ireqr_v(jrecv,inext_v), ierr)
        endif
     enddo

     do jrecv = 1, nrecvs_m(1)
        if (ireqr_m(jrecv,inext_m) /= MPI_REQUEST_NULL) then
           call MPI_Cancel(ireqr_m(jrecv,inext_m), ierr)
        endif
     enddo

     do jrecv = 1, nrecvs_w(1)
        if (ireqr_w(jrecv,inext_w) /= MPI_REQUEST_NULL) then
           call MPI_Cancel(ireqr_w(jrecv,inext_w), ierr)
        endif
     enddo

     if (nudflag > 0 .and. nudnxp > 0) then
        do jrecv = 1, nrecvs_wnud
           if (ireqr_wnud(jrecv,inext_wnud) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_wnud(jrecv,inext_wnud), ierr)
           endif
        enddo
     endif

     do jrecv = 1, nrecvs_wsfc
        if (ireqr_wsfc(jrecv,inext_wsfc) /= MPI_REQUEST_NULL) then
           call MPI_Cancel(ireqr_wsfc(jrecv,inext_wsfc), ierr)
        endif
     enddo

  endif

  ! Now terminate MPI

  call MPI_Finalize(ierr)

  ! Deallocate unused arrays

  if (allocated(send_v )) deallocate(send_v)
  if (allocated(recv_v )) deallocate(recv_v)

  if (allocated(send_m )) deallocate(send_m)
  if (allocated(recv_m )) deallocate(recv_m)

  if (allocated(send_w )) deallocate(send_w)
  if (allocated(recv_w )) deallocate(recv_w)

  if (allocated(send_wnud )) deallocate(send_wnud)
  if (allocated(recv_wnud )) deallocate(recv_wnud)

  if (allocated(send_wsfc )) deallocate(send_wsfc)
  if (allocated(recv_wsfc )) deallocate(recv_wsfc)

#endif

end subroutine olam_mpi_finalize

!==================================================================

subroutine olam_alloc_mpi(mza, mrls)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_ijtabs, only: itab_v, itabg_v, jtab_v, itab_w, itabg_w, jtab_w, &
                        itab_m, itabg_m, jtab_m, mloops
  use mem_nudge,  only: itab_wnud, itabg_wnud, jtab_wnud, nudflag, nudnxp
  use mem_para,   only: nbytes_int, nbytes_real, nbytes_real8,     &
                        nrecvs_v, nrecvs_w, nrecvs_m, nrecvs_wnud, &
                        nsends_v, nsends_w, nsends_m, nsends_wnud, &
                        recv_v,   recv_w,   recv_m,   recv_wnud,   &
                        send_v,   send_w,   send_m,   send_wnud,   &
                        ireqr_v,  ireqr_w,  ireqr_m,  ireqr_wnud,  &
                        ireqs_v,  ireqs_w,  ireqs_m,  ireqs_wnud,  &
                        icurr_v,  icurr_w,  icurr_m,  icurr_wnud,  &
                        inext_v,  inext_w,  inext_m,  inext_wnud,  &
                        itagv,    itagw,    itagm,    itagwnud

  use var_tables, only: nvar_par

  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: mrls

#ifdef OLAM_MPI

  integer :: nbytes_per_iv
  integer :: nbytes_per_iw
  integer :: nbytes_per_im
  integer :: nbytes_per_iwnud

  integer :: itag10  =  10
  integer :: itag110 = 110
  integer :: itag210 = 210
  integer :: itag310 = 310

  integer :: itag11  =  11
  integer :: itag111 = 111
  integer :: itag211 = 211
  integer :: itag311 = 311

  integer :: ierr, ipos
  integer :: jsend, jrecv, jtmp, j

  integer :: nv
  integer :: mrl

  integer, allocatable :: iglobe(:)

! Post V receives

  allocate( ireqr_v(nrecvs_v(1),2) )
  allocate( ireqs_v(nsends_v(1),2) )

  do jrecv = 1,nrecvs_v(1)
     allocate( recv_v(jrecv)%npts(mrls) )
     call MPI_Irecv(recv_v(jrecv)%npts, mrls, MPI_INTEGER, recv_v(jrecv)%iremote, &
          itag10, MPI_COMM_WORLD, ireqr_v(jrecv,icurr_v), ierr)
  enddo

! Post M receives

  allocate( ireqr_m(nrecvs_m(1),2) )
  allocate( ireqs_m(nsends_m(1),2) )

  do jrecv = 1,nrecvs_m(1)
     allocate( recv_m(jrecv)%npts(mrls) )
     call MPI_Irecv(recv_m(jrecv)%npts, mrls, MPI_INTEGER, recv_m(jrecv)%iremote, &
          itag110, MPI_COMM_WORLD, ireqr_m(jrecv,icurr_m), ierr)
  enddo

! Post W receives

  allocate( ireqr_w(nrecvs_w(1),2) )
  allocate( ireqs_w(nsends_w(1),2) )

  do jrecv = 1,nrecvs_w(1)
     allocate( recv_w(jrecv)%npts(mrls) )
     call MPI_Irecv(recv_w(jrecv)%npts, mrls, MPI_INTEGER, recv_w(jrecv)%iremote, &
          itag210, MPI_COMM_WORLD, ireqr_w(jrecv,icurr_w), ierr)
  enddo

! Post WNUD receives

  if (nudflag > 0 .and. nudnxp > 0) then

     allocate( ireqr_w(nrecvs_wnud,2) )
     allocate( ireqs_w(nsends_wnud,2) )

     do jrecv = 1,nrecvs_wnud
        allocate( recv_wnud(jrecv)%npts(1) )
        call MPI_Irecv(recv_wnud(jrecv)%npts, 1, MPI_INTEGER, recv_wnud(jrecv)%iremote, &
             itag310, MPI_COMM_WORLD, ireqr_wnud(jrecv,icurr_wnud), ierr)
     enddo
  endif

! Determine number of bytes to send per IV column

  nbytes_per_iv =       3 * nbytes_int  &
                + mza * 4 * nbytes_real

! Loop over all V sends for mrl = 1

  do jsend = 1, nsends_v(1)

     ! Send buffer sizes to receive ranks

     call MPI_Isend(jtab_v(mloops+jsend)%jend, mrls, MPI_INTEGER, &
                    send_v(jsend)%iremote, itag10, MPI_COMM_WORLD, &
                    ireqs_v(jsend,icurr_v), ierr)

     ! Determine size of send_v buffer for mrl = 1

     send_v(jsend)%nbytes = nbytes_per_iv * jtab_v(mloops+jsend)%jend(1)

     ! Allocate buffer

     allocate(send_v(jsend)%buff( send_v(jsend)%nbytes) )

     ! If at least 1 V point needs to be sent to current remote rank for
     ! current mrl, increase nsends_v(mrl) by 1.

     do mrl = 2, mrls
        if (jtab_v(mloops+jsend)%jend(mrl) > 0) then
           nsends_v(mrl) = nsends_v(mrl) + 1
        endif
     enddo

  enddo

! Determine number of bytes to send per IM column

  nbytes_per_im = mza * 2 * nbytes_real

! Loop over all M sends for mrl = 1

  do jsend = 1, nsends_m(1)

     ! Send buffer sizes to receive ranks

     call MPI_Isend(jtab_m(mloops+jsend)%jend, mrls, MPI_INTEGER, &
                    send_m(jsend)%iremote, itag110, MPI_COMM_WORLD, &
                    ireqs_m(jsend,icurr_m), ierr)

     ! Determine size of send_m buffer for mrl = 1

     send_m(jsend)%nbytes = nbytes_per_im * jtab_m(mloops+jsend)%jend(1)

     ! Allocate buffer

     allocate(send_m(jsend)%buff( send_m(jsend)%nbytes) )

     ! If at least 1 M point needs to be sent to current remote rank for
     ! current mrl, increase nsends_m(mrl) by 1.

     do mrl = 2, mrls
        if (jtab_m(mloops+jsend)%jend(mrl) > 0) then
           nsends_m(mrl) = nsends_m(mrl) + 1
        endif
     enddo

  enddo

! Determine number of bytes to send per IW column

  ! Extra room for the 'G' communication group
  nv = max(nvar_par, 20)

  nbytes_per_iw = mza *  2 * nbytes_real8 &
                + mza * nv * nbytes_real

! Loop over all W sends for mrl = 1

  do jsend = 1, nsends_w(1)

     ! Send buffer sizes to receive ranks

     call MPI_Isend(jtab_w(mloops+jsend)%jend, mrls, MPI_INTEGER, &
                    send_w(jsend)%iremote, itag210, MPI_COMM_WORLD, &
                    ireqs_w(jsend,icurr_w), ierr)

     ! Determine size of send_w buffer for mrl = 1

     send_w(jsend)%nbytes = nbytes_per_iw * jtab_w(mloops+jsend)%jend(1)

     ! Allocate buffer

     allocate(send_w(jsend)%buff(send_w(jsend)%nbytes))

     ! If at least 1 W point needs to be sent to current remote rank for
     ! current mrl, increase nsends_w(mrl) by 1.

     do mrl = 2, mrls
        if (jtab_w(mloops+jsend)%jend(mrl) > 0) then
           nsends_w(mrl) = nsends_w(mrl) + 1
        endif
     enddo

  enddo

  if (nudflag > 0 .and. nudnxp > 0) then

     ! Determine number of bytes to send per IWNUD column

     nbytes_per_iwnud = mza * 6 * nbytes_real8

     ! Loop over all WNUD sends

     do jsend = 1, nsends_wnud

        ! Send buffer sizes to receive ranks

        call MPI_Isend(jtab_wnud(jsend)%jend, 1, MPI_INTEGER, &
                       send_wnud(jsend)%iremote, itag310, MPI_COMM_WORLD, &
                       ireqs_wnud(jsend,icurr_wnud), ierr)

        ! Determine size of send_wnud buffer

        send_wnud(jsend)%nbytes = nbytes_per_iwnud * jtab_wnud(jsend)%jend

        ! Allocate buffer

        allocate(send_wnud (jsend)%buff(send_wnud(jsend)%nbytes))

     enddo
  endif

! Loop over all V receives for mrl = 1

  do jtmp = 1, nrecvs_v(1)

     ! Get recv_v buffer sizes from any node

     call MPI_Waitany(nrecvs_v(1), ireqr_v(:,icurr_v), jrecv, MPI_STATUS_IGNORE, ierr)

     recv_v(jrecv)%nbytes = nbytes_per_iv * recv_v(jrecv)%npts(1)

     ! Allocate recv_v buffers

     allocate(recv_v(jrecv)%buff( recv_v(jrecv)%nbytes ) )
     allocate(recv_v(jrecv)%ipts( recv_v(jrecv)%npts(1) ) )

     ! Loop over all mrl values greater than 1. If at least 1 V point needs to be
     ! received from current remote rank for each mrl, increase nrecvs_v(mrl) by 1.

     do mrl = 2, mrls
        if (recv_v(jrecv)%npts(mrl) > 0) then
           nrecvs_v(mrl) = nrecvs_v(mrl) + 1
        endif
     enddo

     ! Post receive to get the list of iv points we are getting

     call MPI_Irecv(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, MPI_PACKED, &
                    recv_v(jrecv)%iremote, itag11, MPI_COMM_WORLD,        &
                    ireqr_v(jrecv,inext_v), ierr)
  enddo

! Loop over all M receives for mrl = 1

  do jtmp = 1, nrecvs_m(1)

     ! Get recv_m buffer sizes from any node

     call MPI_Waitany(nrecvs_m(1), ireqr_m(:,icurr_m), jrecv, MPI_STATUS_IGNORE, ierr)

     recv_m(jrecv)%nbytes = nbytes_per_im * recv_m(jrecv)%npts(1)

     ! Allocate recv_m buffers

     allocate(recv_m(jrecv)%buff( recv_m(jrecv)%nbytes) )
     allocate(recv_m(jrecv)%ipts( recv_m(jrecv)%npts(1) ) )

     ! Loop over all mrl values greater than 1. If at least 1 M point needs to be
     ! received from current remote rank for each mrl, increase nrecvs_m(mrl) by 1.

     do mrl = 2, mrls
        if (recv_m(jrecv)%npts(mrl) > 0) then
           nrecvs_m(mrl) = nrecvs_m(mrl) + 1
        endif
     enddo

     ! Post receive to get the list of im points we are getting

     call MPI_Irecv(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, MPI_PACKED, &
                    recv_m(jrecv)%iremote, itag111, MPI_COMM_WORLD,       &
                    ireqr_m(jrecv,inext_m), ierr)
  enddo

! Loop over all W receives for mrl = 1

  do jtmp = 1, nrecvs_w(1)

     ! Get recv_w buffer sizes from any node

     call MPI_Waitany(nrecvs_w(1), ireqr_w(:,icurr_w), jrecv, MPI_STATUS_IGNORE, ierr)

     recv_w(jrecv)%nbytes = nbytes_per_iw * recv_w(jrecv)%npts(1)

     ! Allocate recv_w buffers

     allocate(recv_w(jrecv)%buff( recv_w(jrecv)%nbytes))
     allocate(recv_w(jrecv)%ipts( recv_w(jrecv)%npts(1) ) )

     ! Loop over all mrl values greater than 1. If at least 1 W point needs to be
     ! received from current remote rank for each mrl, increase nrecvs_w(mrl) by 1.

     do mrl = 2, mrls
        if (recv_w(jrecv)%npts(mrl) > 0) then
           nrecvs_w(mrl) = nrecvs_w(mrl) + 1
        endif
     enddo

     ! Post receive to get the list of iw points we are getting

     call MPI_Irecv(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, MPI_PACKED, &
                    recv_w(jrecv)%iremote, itag211, MPI_COMM_WORLD,       &
                    ireqr_w(jrecv,inext_w), ierr)
  enddo

! Loop over all WNUD receives

  if (nudflag > 0 .and. nudnxp > 0) then
     do jtmp = 1, nrecvs_wnud

        ! Get recv_wnud buffer sizes from any node

        call MPI_Waitany(nrecvs_wnud, ireqr_wnud(:,icurr_wnud), jrecv, MPI_STATUS_IGNORE, ierr)

        recv_wnud(jrecv)%nbytes = nbytes_per_iwnud * recv_wnud(jrecv)%npts(1)

        ! Allocate recv_wnud buffers

        allocate(recv_wnud(jrecv)%buff( recv_wnud(jrecv)%nbytes) )
        allocate(recv_wnud(jrecv)%ipts( recv_wnud(jrecv)%npts(1) ) )

        ! Post receive to get the list of iwnud points we are getting

        call MPI_Irecv(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, MPI_PACKED, &
                       recv_wnud(jrecv)%iremote, itag311, MPI_COMM_WORLD, &
                       ireqr_wnud(jrecv,inext_w), ierr)
     enddo
  endif

  ! Communicate the list of iv points we are sending to adjacent nodes

  do jtmp = 1, nsends_v(1)

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_v(1), ireqs_v(:,icurr_v), jsend, MPI_STATUS_IGNORE, ierr)

     allocate( iglobe( jtab_v(mloops+jsend)%jend(1) ) )

     do j = 1, jtab_v(mloops+jsend)%jend(1)
        iglobe(j) = itab_v( jtab_v(mloops+jsend)%iv(j) )%ivglobe
     enddo

     ipos = 0

     call MPI_Pack(iglobe, jtab_v(mloops+jsend)%jend(1), MPI_INTEGER, &
          send_v(jsend)%buff, send_v(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_v(jsend)%buff, ipos, MPI_PACKED, &
                    send_v(jsend)%iremote, itag11, MPI_COMM_WORLD, &
                    ireqs_v(jsend,inext_v), ierr)

     deallocate(iglobe)
  enddo

  ! Communicate the list of im points we are sending to adjacent nodes

  do jtmp = 1, nsends_m(1)

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_m(1), ireqs_m(:,icurr_m), jsend, MPI_STATUS_IGNORE, ierr)

     allocate( iglobe( jtab_m(mloops+jsend)%jend(1) ) )

     do j = 1, jtab_m(mloops+jsend)%jend(1)
        iglobe(j) = itab_m( jtab_m(mloops+jsend)%im(j) )%imglobe
     enddo

     ipos = 0

     call MPI_Pack(iglobe, jtab_m(mloops+jsend)%jend(1), MPI_INTEGER, &
          send_m(jsend)%buff, send_m(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_m(jsend)%buff, ipos, MPI_PACKED,           &
                    send_m(jsend)%iremote, itag111, MPI_COMM_WORLD, &
                    ireqs_m(jsend,inext_m), ierr)

     deallocate(iglobe)
  enddo

  ! Communicate the list of iw points we are sending to adjacent nodes

  do jtmp = 1, nsends_w(1)

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_w(1), ireqs_w(:,icurr_w), jsend, MPI_STATUS_IGNORE, ierr)

     allocate( iglobe( jtab_w(mloops+jsend)%jend(1) ) )

     do j = 1, jtab_w(mloops+jsend)%jend(1)
        iglobe(j) = itab_w( jtab_w(mloops+jsend)%iw(j) )%iwglobe
     enddo

     ipos = 0

     call MPI_Pack(iglobe, jtab_w(mloops+jsend)%jend(1), MPI_INTEGER, &
          send_w(jsend)%buff, send_w(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_w(jsend)%buff, ipos, MPI_PACKED,           &
                    send_w(jsend)%iremote, itag211, MPI_COMM_WORLD, &
                    ireqs_w(jsend,inext_w), ierr)

     deallocate(iglobe)
  enddo

  ! Communicate the list of iwnud points we are sending to adjacent nodes

  if (nudflag > 0 .and. nudnxp > 0) then
     do jtmp = 1, nsends_wnud

        ! Make sure previous send is finished
        call MPI_Waitany(nsends_wnud, ireqs_wnud(:,icurr_wnud), jsend, MPI_STATUS_IGNORE, ierr)

        allocate( iglobe( jtab_wnud(jsend)%jend ) )

        do j = 1, jtab_wnud(jsend)%jend
           iglobe(j) = itab_wnud( jtab_wnud(jsend)%iwnud(j) )%iwnudglobe
        enddo

        ipos = 0

        call MPI_Pack(iglobe, jtab_wnud(jsend)%jend, MPI_INTEGER, &
             send_wnud(jsend)%buff, send_wnud(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

        call MPI_Isend(send_wnud(jsend)%buff, ipos, MPI_PACKED,           &
                       send_wnud(jsend)%iremote, itag311, MPI_COMM_WORLD, &
                       ireqs_wnud(jsend,inext_w), ierr)

        deallocate(iglobe)
     enddo
  endif

  ! Unpack and store the list of iv points we are getting from adjacent nodes

  icurr_v = mod(icurr_v,2) + 1
  inext_v = mod(inext_v,2) + 1

  do jtmp = 1, nrecvs_v(1)
     call MPI_Waitany(nrecvs_v(1), ireqr_v(:,icurr_v), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iglobe( recv_v(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, ipos, &
                     iglobe, recv_v(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_v(jrecv)%npts(1)
        recv_v(jrecv)%ipts(j) = itabg_v( iglobe(j) )%iv_myrank
        if (recv_v(jrecv)%ipts(j) < 2) then
           recv_v(jrecv)%ipts(j) = itabg_v( iglobe(j) )%iv_myrank_ivp
        endif
     enddo

     deallocate(iglobe)

     call MPI_Irecv(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, MPI_PACKED, &
                    recv_v(jrecv)%iremote, itagv, MPI_COMM_WORLD,         &
                    ireqr_v(jrecv,inext_v), ierr)
  enddo

  ! Unpack and store the list of im points we are getting from adjacent nodes

  icurr_m = mod(icurr_m,2) + 1
  inext_m = mod(inext_m,2) + 1

  do jtmp = 1, nrecvs_m(1)
     call MPI_Waitany(nrecvs_m(1), ireqr_m(:,icurr_m), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iglobe( recv_m(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, ipos, &
                     iglobe, recv_m(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_m(jrecv)%npts(1)
        recv_m(jrecv)%ipts(j) = itabg_m( iglobe(j) )%im_myrank
        if (recv_m(jrecv)%ipts(j) < 2) then
           recv_m(jrecv)%ipts(j) = itabg_m( iglobe(j) )%im_myrank_imp
        endif
     enddo

     deallocate(iglobe)

     call MPI_Irecv(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, MPI_PACKED, &
                    recv_m(jrecv)%iremote, itagm, MPI_COMM_WORLD,         &
                    ireqr_m(jrecv,inext_m), ierr)
  enddo

  ! Unpack and store the list of iw points we are getting from adjacent nodes

  icurr_w = mod(icurr_w,2) + 1
  inext_w = mod(inext_w,2) + 1

  do jtmp = 1, nrecvs_w(1)
     call MPI_Waitany(nrecvs_w(1), ireqr_w(:,icurr_w), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iglobe( recv_w(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, ipos, &
                     iglobe, recv_w(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_w(jrecv)%npts(1)
        recv_w(jrecv)%ipts(j) = itabg_w( iglobe(j) )%iw_myrank
        if (recv_w(jrecv)%ipts(j) < 2) then
           recv_w(jrecv)%ipts(j) = itabg_w( iglobe(j) )%iw_myrank_iwp
        endif
     enddo

     deallocate(iglobe)

     call MPI_Irecv(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, MPI_PACKED, &
                    recv_w(jrecv)%iremote, itagw, MPI_COMM_WORLD,         &
                    ireqr_w(jrecv,inext_w), ierr)
  enddo

  ! Unpack and store the list of iwnud points we are getting from adjacent nodes

  if (nudflag > 0 .and. nudnxp > 0) then

     icurr_wnud = mod(icurr_wnud,2) + 1
     inext_wnud = mod(inext_wnud,2) + 1

     do jtmp = 1, nrecvs_wnud
        call MPI_Waitany(nrecvs_wnud, ireqr_wnud(:,icurr_wnud), jrecv, MPI_STATUS_IGNORE, ierr)

        ipos = 0

        allocate( iglobe( recv_wnud(jrecv)%npts(1) ) )

        call MPI_Unpack(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, ipos, &
                        iglobe, recv_wnud(jrecv)%npts(1), MPI_INTEGER, &
                        MPI_COMM_WORLD, ierr)

        do j = 1, recv_wnud(jrecv)%npts(1)
           recv_wnud(jrecv)%ipts(j) = itabg_wnud( iglobe(j) )%iwnud_myrank
        enddo

        deallocate(iglobe)

        call MPI_Irecv(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, MPI_PACKED, &
                       recv_wnud(jrecv)%iremote, itagwnud, MPI_COMM_WORLD,         &
                       ireqr_wnud(jrecv,inext_wnud), ierr)
     enddo

  endif

#endif

end subroutine olam_alloc_mpi

!===============================================================================

subroutine mpi_send_v(mrl, rvara1, rvara2, rvara3, rvara4, &
                           i1dvara1, i1dvara2, i1dvara3)

! Subroutine to perform a parallel MPI send of a "V group" of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: send_v, nsends_v, itagv, ireqs_v, icurr_v, inext_v
  use mem_ijtabs, only: jtab_v, mloops
  use mem_grid,   only: mza, mva

  implicit none

  integer, intent(in) :: mrl

  real, optional, intent(in) :: rvara1(mza,mva)
  real, optional, intent(in) :: rvara2(mza,mva)
  real, optional, intent(in) :: rvara3(mza,mva)
  real, optional, intent(in) :: rvara4(mza,mva)

  integer, optional, intent(in) :: i1dvara1(mva)
  integer, optional, intent(in) :: i1dvara2(mva)
  integer, optional, intent(in) :: i1dvara3(mva)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jtmp, jsend
  integer :: j
  integer :: iv

  if (mrl < 1) return

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_v(mrl)

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_v(mrl), ireqs_v(:,icurr_v), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,ierr,j,iv) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1,jtab_v(mloops+jsend)%jend(mrl)
        iv = jtab_v(mloops+jsend)%iv(j)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Pack(rvara3(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Pack(rvara4(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Pack(i1dvara1(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Pack(i1dvara2(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Pack(i1dvara3(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_v(jsend)%buff, ipos, MPI_PACKED, &
                    send_v(jsend)%iremote, itagv, MPI_COMM_WORLD, &
                    ireqs_v(jsend,inext_v), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_v = mod(icurr_v,2) + 1
  inext_v = mod(inext_v,2) + 1

#endif

end subroutine mpi_send_v

!=============================================================================

subroutine mpi_send_m(mrl, rvara1, rvara2)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: nsends_m, send_m, itagm, ireqs_m, icurr_m, inext_m
  use mem_ijtabs, only: jtab_m, mloops
  use mem_grid,   only: mza, mma

  implicit none

  integer,        intent(in) :: mrl
  real, optional, intent(in) :: rvara1(mza,mma)
  real, optional, intent(in) :: rvara2(mza,mma)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jtmp, jsend
  integer :: j, im

  if (mrl < 1) return

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_m(mrl)

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_m(mrl), ireqs_m(:,icurr_m), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,ierr,j,im) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1,jtab_m(mloops+jsend)%jend(mrl)
        im = jtab_m(mloops+jsend)%im(j)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,im),mza,MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,im),mza,MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_m(jsend)%buff, ipos, MPI_PACKED, &
                    send_m(jsend)%iremote, itagm, MPI_COMM_WORLD, &
                    ireqs_m(jsend,inext_m), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_m = mod(icurr_m,2) + 1
  inext_m = mod(inext_m,2) + 1

#endif

end subroutine mpi_send_m

!=============================================================================

subroutine mpi_send_w(mrl, scalars,                                    &
                      rvara1,   rvara2,   rvara3,   rvara4,   rvara5,  &
                      rvara6,   rvara7,   rvara8,   rvara9,   rvara10, &
                      rvara11,  rvara12,  rvara13,  rvara14,  rvara15, &
                      rvara16,  rvara17,  rvara18,  rvara19,  rvara20, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5,&
                      dvara1,   dvara2,   i1dvara1, i1dvara2, i1dvara3,&
                      svara1,   svara2,   svara3,   svara4,   svara5,  &
                      swvar1,   swvar2,   swvar3,   swvar4,   swvar5   )

! Subroutine to perform a parallel MPI send of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use var_tables, only: nvar_par, vtab_r, nptonv
  use mem_ijtabs, only: jtab_w, mloops
  use mem_grid,   only: mza, mwa, lsw, nsw_max
  use consts_coms,only: r8
  use mem_para,   only: send_w, nsends_w, itagw, ireqs_w, icurr_w, inext_w

  implicit none

  integer, intent(in) :: mrl

  character(1), optional, intent(in) :: scalars

  real(r8), optional, intent(in) :: dvara1 (mza,mwa)
  real(r8), optional, intent(in) :: dvara2 (mza,mwa)

  real, optional, contiguous, intent(in) :: svara1(:,:)
  real, optional, contiguous, intent(in) :: svara2(:,:)
  real, optional, contiguous, intent(in) :: svara3(:,:)
  real, optional, contiguous, intent(in) :: svara4(:,:)
  real, optional, contiguous, intent(in) :: svara5(:,:)

  real, optional, intent(in) :: rvara1 (mza,mwa)
  real, optional, intent(in) :: rvara2 (mza,mwa)
  real, optional, intent(in) :: rvara3 (mza,mwa)
  real, optional, intent(in) :: rvara4 (mza,mwa)
  real, optional, intent(in) :: rvara5 (mza,mwa)
  real, optional, intent(in) :: rvara6 (mza,mwa)
  real, optional, intent(in) :: rvara7 (mza,mwa)
  real, optional, intent(in) :: rvara8 (mza,mwa)
  real, optional, intent(in) :: rvara9 (mza,mwa)
  real, optional, intent(in) :: rvara10(mza,mwa)
  real, optional, intent(in) :: rvara11(mza,mwa)
  real, optional, intent(in) :: rvara12(mza,mwa)
  real, optional, intent(in) :: rvara13(mza,mwa)
  real, optional, intent(in) :: rvara14(mza,mwa)
  real, optional, intent(in) :: rvara15(mza,mwa)
  real, optional, intent(in) :: rvara16(mza,mwa)
  real, optional, intent(in) :: rvara17(mza,mwa)
  real, optional, intent(in) :: rvara18(mza,mwa)
  real, optional, intent(in) :: rvara19(mza,mwa)
  real, optional, intent(in) :: rvara20(mza,mwa)

  real, optional, intent(in) :: r1dvara1(mwa)
  real, optional, intent(in) :: r1dvara2(mwa)
  real, optional, intent(in) :: r1dvara3(mwa)
  real, optional, intent(in) :: r1dvara4(mwa)
  real, optional, intent(in) :: r1dvara5(mwa)

  integer, optional, intent(in) :: i1dvara1(mwa)
  integer, optional, intent(in) :: i1dvara2(mwa)
  integer, optional, intent(in) :: i1dvara3(mwa)

  real, optional, intent(in) :: swvar1(nsw_max,mwa)
  real, optional, intent(in) :: swvar2(nsw_max,mwa)
  real, optional, intent(in) :: swvar3(nsw_max,mwa)
  real, optional, intent(in) :: swvar4(nsw_max,mwa)
  real, optional, intent(in) :: swvar5(nsw_max,mwa)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: i, ivar
  integer :: j, iw

  if (mrl < 1) return

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_w(mrl)

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_w(mrl), ireqs_w(:,icurr_w), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,ierr,j,iw,i,ivar) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, jtab_w(mloops+jsend)%jend(mrl)
        iw = jtab_w(mloops+jsend)%iw(j)

        if (present(dvara1)) then
           call MPI_Pack(dvara1(1,iw),mza,MPI_REAL8, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Pack(dvara2(1,iw),mza,MPI_REAL8, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara1)) then
           call MPI_Pack(svara1(:,iw),size(svara1,1),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara2)) then
           call MPI_Pack(svara2(:,iw),size(svara2,1),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara3)) then
           call MPI_Pack(svara3(:,iw),size(svara3,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara4)) then
           call MPI_Pack(svara4(:,iw),size(svara4,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara5)) then
           call MPI_Pack(svara5(:,iw),size(svara5,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Pack(rvara3(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Pack(rvara4(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara5)) then
           call MPI_Pack(rvara5(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara6)) then
           call MPI_Pack(rvara6(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara7)) then
           call MPI_Pack(rvara7(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara8)) then
           call MPI_Pack(rvara8(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara9)) then
           call MPI_Pack(rvara9(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara10)) then
           call MPI_Pack(rvara10(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara11)) then
           call MPI_Pack(rvara11(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara12)) then
           call MPI_Pack(rvara12(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara13)) then
           call MPI_Pack(rvara13(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara14)) then
           call MPI_Pack(rvara14(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara15)) then
           call MPI_Pack(rvara15(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara16)) then
           call MPI_Pack(rvara16(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara17)) then
           call MPI_Pack(rvara17(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara18)) then
           call MPI_Pack(rvara18(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara19)) then
           call MPI_Pack(rvara19(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara20)) then
           call MPI_Pack(rvara20(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara1)) then
           call MPI_Pack(r1dvara1(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara2)) then
           call MPI_Pack(r1dvara2(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara3)) then
           call MPI_Pack(r1dvara3(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara4)) then
           call MPI_Pack(r1dvara4(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara5)) then
           call MPI_Pack(r1dvara5(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Pack(i1dvara1(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Pack(i1dvara2(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Pack(i1dvara3(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar1)) then
           call MPI_Pack(swvar1(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar2)) then
           call MPI_Pack(swvar2(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar3)) then
           call MPI_Pack(swvar3(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar4)) then
           call MPI_Pack(swvar4(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar5)) then
           call MPI_Pack(swvar5(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(scalars)) then
           do i = 1, nvar_par
              ivar = nptonv(i)

              call MPI_Pack(vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL, &
                   send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           enddo
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_w(jsend)%buff, ipos, MPI_PACKED, &
                    send_w(jsend)%iremote, itagw, MPI_COMM_WORLD, &
                    ireqs_w(jsend,inext_w), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_w = mod(icurr_w,2) + 1
  inext_w = mod(inext_w,2) + 1

#endif

end subroutine mpi_send_w

!=============================================================================

subroutine mpi_send_wnud(dvara1, dvara2, dvara3, dvara4, dvara5, dvara6)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: nsends_wnud, send_wnud, itagwnud, ireqs_wnud, &
                        icurr_wnud, inext_wnud
  use mem_grid,   only: mza
  use mem_nudge,  only: mwnud, jtab_wnud
  use consts_coms,only: r8

  implicit none

  real(r8), optional, intent(in) :: dvara1(mza,mwnud)
  real(r8), optional, intent(in) :: dvara2(mza,mwnud)
  real(r8), optional, intent(in) :: dvara3(mza,mwnud)
  real(r8), optional, intent(in) :: dvara4(mza,mwnud)
  real(r8), optional, intent(in) :: dvara5(mza,mwnud)
  real(r8), optional, intent(in) :: dvara6(mza,mwnud)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: j, iwnud

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_wnud

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_wnud, ireqs_wnud(:,icurr_wnud), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,ierr,j,iwnud) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, jtab_wnud(jsend)%jend
        iwnud = jtab_wnud(jsend)%iwnud(j)

        if (present(dvara1)) then
           call MPI_Pack(dvara1(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Pack(dvara2(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara3)) then
           call MPI_Pack(dvara3(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara4)) then
           call MPI_Pack(dvara4(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara5)) then
           call MPI_Pack(dvara5(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara6)) then
           call MPI_Pack(dvara6(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_wnud(jsend)%buff, ipos, MPI_PACKED, &
                    send_wnud(jsend)%iremote, itagwnud, MPI_COMM_WORLD, &
                    ireqs_wnud(jsend, inext_wnud), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_wnud = mod(icurr_wnud,2) + 1
  inext_wnud = mod(inext_wnud,2) + 1

#endif

end subroutine mpi_send_wnud

!=============================================================================

subroutine mpi_recv_v(mrl, rvara1, rvara2, rvara3, rvara4, &
                           i1dvara1, i1dvara2, i1dvara3)

! Subroutine to perform a parallel MPI receive of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: recv_v, nrecvs_v, ireqr_v, itagv, icurr_v, inext_v
  use mem_grid,   only: mza, mva

  implicit none

  integer, intent(in) :: mrl

  real, optional, intent(inout) :: rvara1(mza,mva)
  real, optional, intent(inout) :: rvara2(mza,mva)
  real, optional, intent(inout) :: rvara3(mza,mva)
  real, optional, intent(inout) :: rvara4(mza,mva)

  integer, optional, intent(in) :: i1dvara1(mva)
  integer, optional, intent(in) :: i1dvara2(mva)
  integer, optional, intent(in) :: i1dvara3(mva)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: iv

  if (mrl < 1) return

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_v(mrl)

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_v(mrl), ireqr_v(:,icurr_v), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,iv,ierr) default(shared)

     ipos = 0

     do j = 1, recv_v(jrecv)%npts(mrl)
        iv = recv_v(jrecv)%ipts(j)

        if (present(rvara1)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara1(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara2(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara3(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara4(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara1(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara2(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara3(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Irecv(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, MPI_PACKED, &
                    recv_v(jrecv)%iremote, itagv, MPI_COMM_WORLD,         &
                    ireqr_v(jrecv,inext_v), ierr                          )

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_v

!=============================================================================

subroutine mpi_recv_m(mrl, rvara1, rvara2)

! Subroutine to perform a parallel MPI receive of a "M group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: recv_m, nrecvs_m, ireqr_m, itagm, icurr_m, inext_m
  use mem_grid,   only: mza, mma

  implicit none

  integer,        intent(in)    :: mrl
  real, optional, intent(inout) :: rvara1(mza,mma)
  real, optional, intent(inout) :: rvara2(mza,mma)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: im

  if (mrl < 1) return

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_m(mrl)

     !  Now, let's wait on our receives

     call MPI_Waitany(nrecvs_m(mrl), ireqr_m(:,icurr_m), jrecv, MPI_STATUS_IGNORE, ierr)

     !  We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,im,ierr) default(shared)

     ipos = 0

     do j = 1, recv_m(jrecv)%npts(mrl)
        im = recv_m(jrecv)%ipts(j)

        if (present(rvara1)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                rvara1(1,im),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                rvara2(1,im),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Irecv(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, MPI_PACKED, &
                    recv_m(jrecv)%iremote, itagm, MPI_COMM_WORLD,         &
                    ireqr_m(jrecv,inext_m), ierr                          )

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_m

!=============================================================================

subroutine mpi_recv_w(mrl, scalars,                                    &
                      rvara1,   rvara2,   rvara3,   rvara4,   rvara5,  &
                      rvara6,   rvara7,   rvara8,   rvara9,   rvara10, &
                      rvara11,  rvara12,  rvara13,  rvara14,  rvara15, &
                      rvara16,  rvara17,  rvara18,  rvara19,  rvara20, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5,&
                      dvara1,   dvara2,   i1dvara1, i1dvara2, i1dvara3,&
                      svara1,   svara2,   svara3,   svara4,   svara5,  &
                      swvar1,   swvar2,   swvar3,   swvar4,   swvar5   )

! Subroutine to perform a parallel MPI receive of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use var_tables, only: vtab_r, nvar_par, nptonv
  use mem_para,   only: recv_w, nrecvs_w, ireqr_w, itagw, icurr_w, inext_w
  use mem_grid,   only: mza, mwa, lsw, nsw_max
  use consts_coms,only: r8

  implicit none

  integer, intent(in) :: mrl

  character(1), optional, intent(in) :: scalars

  real(r8), optional, intent(inout) :: dvara1 (mza,mwa)
  real(r8), optional, intent(inout) :: dvara2 (mza,mwa)

  real, optional, contiguous, intent(inout) :: svara1(:,:)
  real, optional, contiguous, intent(inout) :: svara2(:,:)
  real, optional, contiguous, intent(inout) :: svara3(:,:)
  real, optional, contiguous, intent(inout) :: svara4(:,:)
  real, optional, contiguous, intent(inout) :: svara5(:,:)

  real, optional, intent(inout) :: rvara1 (mza,mwa)
  real, optional, intent(inout) :: rvara2 (mza,mwa)
  real, optional, intent(inout) :: rvara3 (mza,mwa)
  real, optional, intent(inout) :: rvara4 (mza,mwa)
  real, optional, intent(inout) :: rvara5 (mza,mwa)
  real, optional, intent(inout) :: rvara6 (mza,mwa)
  real, optional, intent(inout) :: rvara7 (mza,mwa)
  real, optional, intent(inout) :: rvara8 (mza,mwa)
  real, optional, intent(inout) :: rvara9 (mza,mwa)
  real, optional, intent(inout) :: rvara10(mza,mwa)
  real, optional, intent(inout) :: rvara11(mza,mwa)
  real, optional, intent(inout) :: rvara12(mza,mwa)
  real, optional, intent(inout) :: rvara13(mza,mwa)
  real, optional, intent(inout) :: rvara14(mza,mwa)
  real, optional, intent(inout) :: rvara15(mza,mwa)
  real, optional, intent(inout) :: rvara16(mza,mwa)
  real, optional, intent(inout) :: rvara17(mza,mwa)
  real, optional, intent(inout) :: rvara18(mza,mwa)
  real, optional, intent(inout) :: rvara19(mza,mwa)
  real, optional, intent(inout) :: rvara20(mza,mwa)

  real, optional, intent(inout) :: r1dvara1(mwa)
  real, optional, intent(inout) :: r1dvara2(mwa)
  real, optional, intent(inout) :: r1dvara3(mwa)
  real, optional, intent(inout) :: r1dvara4(mwa)
  real, optional, intent(inout) :: r1dvara5(mwa)

  integer, optional, intent(inout) :: i1dvara1(mwa)
  integer, optional, intent(inout) :: i1dvara2(mwa)
  integer, optional, intent(inout) :: i1dvara3(mwa)

  real, optional, intent(inout) :: swvar1(nsw_max,mwa)
  real, optional, intent(inout) :: swvar2(nsw_max,mwa)
  real, optional, intent(inout) :: swvar3(nsw_max,mwa)
  real, optional, intent(inout) :: swvar4(nsw_max,mwa)
  real, optional, intent(inout) :: swvar5(nsw_max,mwa)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, ivar, jtmp
  integer :: i, j
  integer :: iw
!  integer :: iwglobe

  if (mrl < 1) return

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_w(mrl)

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_w(mrl), ireqr_w(:,icurr_w), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,iw,ierr,i,ivar) default(shared)

     ipos = 0

     do j = 1, recv_w(jrecv)%npts(mrl)
        iw = recv_w(jrecv)%ipts(j)

        if (present(dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                dvara1(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                dvara2(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara1(:,iw),size(svara1,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara2(:,iw),size(svara2,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara3(:,iw),size(svara3,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara4(:,iw),size(svara4,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara5(:,iw),size(svara5,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara1(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara2(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara3(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara4(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara5(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara6)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara6(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara7)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara7(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara8)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara8(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara9)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara9(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara10)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara10(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara11)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara11(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara12)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara12(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara13)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara13(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara14)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara14(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara15)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara15(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara16)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara16(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara17)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara17(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara18)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara18(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara19)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara19(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara20)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara20(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara1(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara2(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara3(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara4(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara5(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara1(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara2(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara3(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar1(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar2(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar3(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar4(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar5(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(scalars)) then
           do i = 1, nvar_par
              ivar = nptonv(i)

              call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                   vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           enddo
        endif

     enddo

     call MPI_Irecv(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, MPI_PACKED, &
                    recv_w(jrecv)%iremote, itagw, MPI_COMM_WORLD,         &
                    ireqr_w(jrecv,inext_w), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_w

!=============================================================================

subroutine mpi_recv_wnud(dvara1, dvara2, dvara3, dvara4, dvara5, dvara6)

! Subroutine to perform a parallel MPI receive of a "WNUD group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,    only: recv_wnud, nrecvs_wnud, ireqr_wnud, itagwnud, &
                         icurr_wnud, inext_wnud
  use mem_grid,    only: mza
  use mem_nudge,   only: mwnud
  use consts_coms, only: r8

  implicit none

  real(r8), optional, intent(inout) :: dvara1(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara2(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara3(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara4(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara5(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara6(mza,mwnud)

#ifdef OLAM_MPI

  integer  :: ierr, ipos
  integer  :: jrecv, jtmp
  integer  :: j, iwnud
  real(r8) :: vctr1(mza)

  do jtmp = 1, nrecvs_wnud

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_wnud, ireqr_wnud(:,icurr_wnud), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     ! For now we don't use OpenMP to unpack these buffers because more then one receive
     ! may be added to an individual iwnud point.

     ipos = 0

     do j = 1, recv_wnud(jrecv)%npts(1)
        iwnud = recv_wnud(jrecv)%ipts(j)

        if (present(dvara1)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara1(2:mza,iwnud) = dvara1(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara2)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara2(2:mza,iwnud) = dvara2(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara3)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara3(2:mza,iwnud) = dvara3(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara4)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara4(2:mza,iwnud) = dvara4(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara5)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara5(2:mza,iwnud) = dvara5(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara6)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara6(2:mza,iwnud) = dvara6(2:mza,iwnud) + vctr1(2:mza)
        endif

     enddo

     call MPI_Irecv(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, MPI_PACKED, &
                    recv_wnud(jrecv)%iremote, itagwnud, MPI_COMM_WORLD,         &
                    ireqr_wnud(jrecv,inext_wnud), ierr)
  enddo

#endif

end subroutine mpi_recv_wnud

!=============================================================================

subroutine olam_stop(message)

#ifdef OLAM_MPI
  use mem_para, only: myrank
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI
  integer :: ierr
#endif

  character(*), intent(in) :: message

#ifdef OLAM_MPI
  write(*,'(A,I0,A)') "Node ", myrank, ":"
  write(*,'(A)') "STOP "//message
  call mpi_abort(MPI_COMM_WORLD,1,ierr)
  stop
#else
  write(*,*) "STOPPING: "//message
  stop
#endif

end subroutine olam_stop

!================================================================================

subroutine olam_mpi_barrier()

#ifdef OLAM_MPI
  use mpi
#endif
  implicit none

#ifdef OLAM_MPI
  integer :: ierr
  call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

end subroutine olam_mpi_barrier

!=============================================================================

End Module olam_mpi_atm
