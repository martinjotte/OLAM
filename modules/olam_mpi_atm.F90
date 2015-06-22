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

  use mem_para,  only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8
  implicit none

#ifdef OLAM_MPI

  integer :: ierr

! Initialize MPI and determine process groupsize and myrank
  
  call MPI_Init(ierr)
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

  use mem_para,  only: mgroupsize,     &
                       send_v, recv_v, &
                       send_w, recv_w, &
                       send_m, recv_m, &
                       send_wnud, recv_wnud

  use misc_coms, only: iparallel
  implicit none

  if (iparallel == 0) return

! Allocate send and recv tables

  allocate(send_v (mgroupsize))
  allocate(recv_v (mgroupsize))

  allocate(send_m (mgroupsize))
  allocate(recv_m (mgroupsize))

  allocate(send_w (mgroupsize))
  allocate(recv_w (mgroupsize))

  allocate(send_wnud(mgroupsize))
  allocate(recv_wnud(mgroupsize))

end subroutine alloc_mpi_sndrcv_bufs

!===============================================================================

subroutine olam_mpi_finalize()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para, only: send_v, recv_v, &
                      send_w, recv_w, &
                      send_m, recv_m, &
                      send_wnud, recv_wnud

  implicit none

#ifdef OLAM_MPI

  integer :: ierr

  call MPI_Finalize(ierr)
  
  if (allocated(send_v))  deallocate(send_v)
  if (allocated(recv_v))  deallocate(recv_v)

  if (allocated(send_m))  deallocate(send_m)
  if (allocated(recv_m))  deallocate(recv_m)

  if (allocated(send_w))  deallocate(send_w)
  if (allocated(recv_w))  deallocate(recv_w)

  if (allocated(send_wnud))  deallocate(send_wnud)
  if (allocated(recv_wnud))  deallocate(recv_wnud)

#endif

end subroutine olam_mpi_finalize

!==================================================================

subroutine olam_alloc_mpi(mza, mrls)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_ijtabs, only: jtab_v, jtab_w, jtab_m, mloops
  use mem_nudge,  only: jtab_wnud
  use mem_para,   only: myrank,   nbytes_int, nbytes_real, nbytes_real8,&
                        nrecvs_v, nrecvs_w,   nrecvs_m,    nrecvs_wnud, &
                        nsends_v, nsends_w,   nsends_m,    nsends_wnud, &
                        recv_v,   recv_w,     recv_m,      recv_wnud,   &
                        send_v,   send_w,     send_m,      send_wnud
  use misc_coms,  only: io6
  use var_tables, only: nvar_par

  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: mrls

#ifdef OLAM_MPI

  integer :: nbytes_per_iu
  integer :: nbytes_per_iv
  integer :: nbytes_per_iw
  integer :: nbytes_per_im
  integer :: nbytes_per_iwnud

  integer :: itag10  =  10
  integer :: itag110 = 110
  integer :: itag210 = 210
  integer :: itag310 = 310

  integer :: ierr
  integer :: jsend
  integer :: jrecv
  integer :: jtmp

  integer :: nupts, nvpts, nwpts, nwnudpts
  integer :: nv
  integer :: mrl

  integer              :: sbuf(mrls)
  integer, allocatable :: rbuf(:,:)

  integer              :: wsbuf(mrls)
  integer, allocatable :: wrbuf(:,:)

  integer              :: msbuf(mrls)
  integer, allocatable :: mrbuf(:,:)

  integer              :: wnudsbuf
  integer, allocatable :: wnudrbuf(:)

! Post V receives

  allocate(rbuf(mrls,nrecvs_v(1)))

  do jrecv = 1,nrecvs_v(1)
     call MPI_IRecv(rbuf(1,jrecv), mrls, MPI_INTEGER, recv_v(jrecv)%iremote, &
          itag10, MPI_COMM_WORLD, recv_v(jrecv)%irequest, ierr)
  enddo

! Post M receives

  allocate(mrbuf(mrls,nrecvs_m(1)))

  do jrecv = 1,nrecvs_m(1)
     call MPI_IRecv(mrbuf(1,jrecv), mrls, MPI_INTEGER, recv_m(jrecv)%iremote, &
          itag210, MPI_COMM_WORLD, recv_m(jrecv)%irequest, ierr)
  enddo

! Post W receives

  allocate(wrbuf(mrls,nrecvs_w(1)))

  do jrecv = 1,nrecvs_w(1)

     call MPI_IRecv(wrbuf(1,jrecv), mrls, MPI_INTEGER, recv_w(jrecv)%iremote, &
          itag110, MPI_COMM_WORLD, recv_w(jrecv)%irequest, ierr)

  enddo

! Post WNUD receives

  allocate(wnudrbuf(nrecvs_wnud))

  do jrecv = 1,nrecvs_wnud

     call MPI_IRecv(wnudrbuf(jrecv), 1, MPI_INTEGER, recv_wnud(jrecv)%iremote, &
          itag310, MPI_COMM_WORLD, recv_wnud(jrecv)%irequest, ierr)

  enddo

! Determine number of bytes to send per IV column

  nbytes_per_iv = nbytes_int &
                + mza * 4 * nbytes_real

! Loop over all V sends for mrl = 1

  do jsend = 1,nsends_v(1)

! Determine size of send_v buffer for mrl = 1

     send_v(jsend)%nbytes = nbytes_int &
                          + nbytes_per_iv * jtab_v(mloops+jsend)%jend(1)

! Allocate buffer

     allocate(send_v(jsend)%buff(send_v(jsend)%nbytes))

! Send buffer sizes to receive ranks

     sbuf(1) = send_v(jsend)%nbytes

     do mrl = 2,mrls
        sbuf(mrl) = jtab_v(mloops+jsend)%jend(mrl)

! If at least 1 V point needs to be sent to current remote rank for 
! current mrl, increase nsends_v(mrl) by 1.

        if (jtab_v(mloops+jsend)%jend(mrl) > 0) &
             nsends_v(mrl) = nsends_v(mrl) + 1

     enddo

     call MPI_Send(sbuf, mrls, MPI_INTEGER, send_v(jsend)%iremote, itag10, &
          MPI_COMM_WORLD,ierr)

  enddo

 ! Determine number of bytes to send per IM column

  nbytes_per_im = nbytes_int &
                + mza * 2 * nbytes_real

! Loop over all V sends for mrl = 1

  do jsend = 1, nsends_m(1)

! Determine size of send_m buffer for mrl = 1

     send_m(jsend)%nbytes = nbytes_int &
                          + nbytes_per_im * jtab_m(mloops+jsend)%jend(1)

! Allocate buffer

     allocate(send_m(jsend)%buff(send_m(jsend)%nbytes))

! Send buffer sizes to receive ranks

     msbuf(1) = send_m(jsend)%nbytes

     do mrl = 2,mrls
        msbuf(mrl) = jtab_v(mloops+jsend)%jend(mrl)

! If at least 1 M point needs to be sent to current remote rank for 
! current mrl, increase nsends_m(mrl) by 1.

        if (jtab_m(mloops+jsend)%jend(mrl) > 0) &
             nsends_m(mrl) = nsends_m(mrl) + 1

     enddo

     call MPI_Send(msbuf, mrls, MPI_INTEGER, send_m(jsend)%iremote, itag210, &
          MPI_COMM_WORLD,ierr)

  enddo

! Determine number of bytes to send per IW column

  ! Extra room for the 'G' communication group
  nv = max(nvar_par, 15)

  nbytes_per_iw = nbytes_int &
                + mza * max( 3*nbytes_real8 + 5*nbytes_real, nv*nbytes_real)

! Loop over all W sends for mrl = 1
  
  do jsend = 1, nsends_w(1)

! Determine size of send_w buffer for mrl = 1

     send_w(jsend)%nbytes = nbytes_int &
                          + nbytes_per_iw * jtab_w(mloops+jsend)%jend(1)

! Allocate buffer

     allocate(send_w(jsend)%buff(send_w(jsend)%nbytes))

! Send buffer sizes to receive ranks

     wsbuf(1) = send_w(jsend)%nbytes

     do mrl = 2,mrls
        wsbuf(mrl) = jtab_w(mloops+jsend)%jend(mrl)

! If at least 1 W point needs to be sent to current remote rank for 
! current mrl, increase nsends_w(mrl) by 1.

        if (jtab_w(mloops+jsend)%jend(mrl) > 0) nsends_w(mrl) = nsends_w(mrl) + 1

     enddo

     call MPI_Send(wsbuf, mrls, MPI_INTEGER, send_w(jsend)%iremote, itag110, &
          MPI_COMM_WORLD, ierr)

  enddo

! Determine number of bytes to send per IWNUD column

  nbytes_per_iwnud = nbytes_int + mza * 6 * nbytes_real

! Loop over all WNUD sends
  
  do jsend = 1, nsends_wnud

! Determine size of send_wnud buffer

     send_wnud(jsend)%nbytes = nbytes_int &
                             + nbytes_per_iwnud * jtab_wnud(jsend)%jend

! Allocate buffer

     allocate(send_wnud(jsend)%buff(send_wnud(jsend)%nbytes))

! Send buffer sizes to receive ranks

     wnudsbuf = send_wnud(jsend)%nbytes

     call MPI_Send(wnudsbuf, 1, MPI_INTEGER, send_wnud(jsend)%iremote, &
                   itag310, MPI_COMM_WORLD, ierr)

  enddo

! Loop over all V receives for mrl = 1

  do jtmp = 1, nrecvs_v(1)

! Get recv_v buffer sizes from any node

     call MPI_Waitany(nrecvs_v(1), recv_v(1:nrecvs_v(1))%irequest, jrecv, &
          MPI_STATUS_IGNORE, ierr)

     recv_v(jrecv)%nbytes  = rbuf(1,jrecv)

! Allocate recv_v buffers

     allocate(recv_v(jrecv)%buff(recv_v(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 V point needs to be 
! received from current remote rank for each mrl, increase nrecvs_v(mrl) by 1.

     do mrl = 2, mrls
        if (rbuf(mrl,jrecv) > 0) nrecvs_v(mrl) = nrecvs_v(mrl) + 1
     enddo

  enddo

! Loop over all M receives for mrl = 1

  do jtmp = 1, nrecvs_m(1)

! Get recv_m buffer sizes from any node

     call MPI_Waitany(nrecvs_m(1), recv_m(1:nrecvs_m(1))%irequest, jrecv, &
          MPI_STATUS_IGNORE, ierr)

     recv_m(jrecv)%nbytes = mrbuf(1,jrecv)

! Allocate recv_m buffers

     allocate(recv_m(jrecv)%buff(recv_m(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 M point needs to be 
! received from current remote rank for each mrl, increase nrecvs_m(mrl) by 1.

     do mrl = 2, mrls
        if (mrbuf(mrl,jrecv) > 0) nrecvs_m(mrl) = nrecvs_m(mrl) + 1
     enddo

  enddo

! Loop over all W receives for mrl = 1

  do jtmp = 1,nrecvs_w(1)

! Get recv_w buffer sizes

     call MPI_Waitany(nrecvs_w(1), recv_w(1:nrecvs_w(1))%irequest, jrecv, &
          MPI_STATUS_IGNORE, ierr)

! Allocate recv_w buffers

     recv_w(jrecv)%nbytes = wrbuf(1,jrecv)
     allocate(recv_w(jrecv)%buff(recv_w(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 W point needs to be 
! received from current remote rank for each mrl, increase nrecvs_w(mrl) by 1.

     do mrl = 2,mrls
        if (wrbuf(mrl,jrecv) > 0) nrecvs_w(mrl) = nrecvs_w(mrl) + 1
     enddo

  enddo

! Loop over all WNUD receives

  do jtmp = 1,nrecvs_wnud

! Get recv_wnud buffer sizes

     call MPI_Waitany(nrecvs_wnud, recv_wnud(1:nrecvs_wnud)%irequest, jrecv, &
          MPI_STATUS_IGNORE, ierr)

! Allocate recv_wnud buffers

     recv_wnud(jrecv)%nbytes = wnudrbuf(jrecv)
     allocate(recv_wnud(jrecv)%buff(recv_wnud(jrecv)%nbytes))

  enddo

#endif

  return
end subroutine olam_alloc_mpi

!===============================================================================

subroutine mpi_send_v(mrl, rvara1, rvara2, rvara3, rvara4)

! Subroutine to perform a parallel MPI send of a "V group" of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_para,   only: send_v, recv_v, nsends_v, nrecvs_v, itagv

use mem_ijtabs, only: itab_v, jtab_v, mrl_begs, istp, mloops
use mem_grid,   only: mza, mva
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(in) :: rvara1(mza,mva)
real, optional, intent(in) :: rvara2(mza,mva)
real, optional, intent(in) :: rvara3(mza,mva)
real, optional, intent(in) :: rvara4(mza,mva)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: j
integer :: iv
integer :: ivglobe

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_v(mrl)

   call MPI_Irecv(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,MPI_PACKED, &
                  recv_v(jrecv)%iremote,itagv,MPI_COMM_WORLD,         &
                  recv_v(jrecv)%irequest,ierr                         )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_v(mrl)

   ipos = 0

   call MPI_Pack(jtab_v(mloops+jsend)%jend(mrl),1,MPI_INTEGER, &
      send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
   do j = 1,jtab_v(mloops+jsend)%jend(mrl)
      iv = jtab_v(mloops+jsend)%iv(j)
      ivglobe = itab_v(iv)%ivglobe
!----------------------------------------------------------------

      call MPI_Pack(ivglobe,1,MPI_INTEGER, &
         send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

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

   enddo

   call MPI_Isend(send_v(jsend)%buff,ipos,MPI_PACKED,         &
                  send_v(jsend)%iremote,itagv,MPI_COMM_WORLD, &
                  send_v(jsend)%irequest,ierr                 )

enddo

#endif

return
end subroutine mpi_send_v

!=============================================================================

subroutine mpi_send_m(mrl, rvara1, rvara2)

#ifdef OLAM_MPI
  use mpi
#endif
  
  use misc_coms,  only: io6
  use mem_para,   only: nrecvs_m, nsends_m, recv_m, send_m, itagm
  use mem_ijtabs, only: jtab_m, itab_m, mloops
  use mem_grid,   only: mza, mma

  implicit none

  integer,        intent(in) :: mrl
  real, optional, intent(in) :: rvara1(mza,mma)
  real, optional, intent(in) :: rvara2(mza,mma)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jsend, ivar
  integer :: j
  integer :: im
  integer :: imglobe

  if (mrl < 1) return

! Before we send anything, post the receives

  do jrecv = 1, nrecvs_m(mrl)

     call MPI_Irecv(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,MPI_PACKED, &
                    recv_m(jrecv)%iremote,itagm,MPI_COMM_WORLD,         &
                    recv_m(jrecv)%irequest,ierr                         )
  enddo

! Now we can actually go on to sending the stuff

  do jsend = 1, nsends_m(mrl)

     ipos = 0

     call MPI_Pack(jtab_m(mloops+jsend)%jend(mrl),1,MPI_INTEGER, &
          send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
     do j = 1,jtab_m(mloops+jsend)%jend(mrl)
        im = jtab_m(mloops+jsend)%im(j)
        imglobe = itab_m(im)%imglobe
!----------------------------------------------------------------
        
        call MPI_Pack(imglobe,1,MPI_INTEGER, &
             send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,im),mza,MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,im),mza,MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Isend(send_m(jsend)%buff,ipos,MPI_PACKED,         &
                    send_m(jsend)%iremote,itagm,MPI_COMM_WORLD, &
                    send_m(jsend)%irequest,ierr                 )
     
  enddo

#endif

end subroutine mpi_send_m

!=============================================================================

subroutine mpi_send_w(mrl, scalars, dvara1, dvara2,                      &
                      rvara1, rvara2, rvara3, rvara4,  rvara5,  rvara6,  &
                      rvara7, rvara8, rvara9, rvara10, rvara11, rvara12, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5   )

! Subroutine to perform a parallel MPI send of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use var_tables, only: nvar_par, vtab_r, num_scalar, nptonv
use mem_ijtabs, only: jtab_w, itab_w, mrl_begs, mrl_begl, mrl_endl, istp, mloops
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use consts_coms,only: r8
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w, itagw

implicit none

integer, intent(in) :: mrl

character(1), optional, intent(in) :: scalars

real(r8), optional, intent(in) :: dvara1 (mza,mwa)
real(r8), optional, intent(in) :: dvara2 (mza,mwa)

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

real, optional, intent(in) :: r1dvara1(mwa)
real, optional, intent(in) :: r1dvara2(mwa)
real, optional, intent(in) :: r1dvara3(mwa)
real, optional, intent(in) :: r1dvara4(mwa)
real, optional, intent(in) :: r1dvara5(mwa)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: i, j
integer :: iw
integer :: iwglobe

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_w(mrl)

   call MPI_Irecv(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,MPI_PACKED, &
                  recv_w(jrecv)%iremote,itagw,MPI_COMM_WORLD,         &
                  recv_w(jrecv)%irequest,ierr                         )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_w(mrl)

   ipos = 0

   call MPI_Pack(jtab_w(mloops+jsend)%jend(mrl),1,MPI_INTEGER, &
      send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
   do j = 1,jtab_w(mloops+jsend)%jend(mrl)
      iw = jtab_w(mloops+jsend)%iw(j)
      iwglobe = itab_w(iw)%iwglobe
!----------------------------------------------------------------

      call MPI_Pack(iwglobe,1,MPI_INTEGER, &
         send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (present(dvara1)) then
         call MPI_Pack(dvara1(1,iw),mza,MPI_REAL8, &
         send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
      endif

      if (present(dvara2)) then
         call MPI_Pack(dvara2(1,iw),mza,MPI_REAL8, &
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

      if (present(scalars)) then
         do i = 1, nvar_par
            ivar = nptonv(i)

            call MPI_Pack(vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL, &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         enddo
      endif

   enddo

   call MPI_Isend(send_w(jsend)%buff,ipos,MPI_PACKED,         &
                  send_w(jsend)%iremote,itagw,MPI_COMM_WORLD, &
                  send_w(jsend)%irequest,ierr                 )

enddo

#endif

return
end subroutine mpi_send_w

!=============================================================================

subroutine mpi_send_wnud(rvara1, rvara2, rvara3, rvara4, rvara5, rvara6)

#ifdef OLAM_MPI
  use mpi
#endif
  
  use misc_coms,  only: io6
  use mem_para,   only: nrecvs_wnud, nsends_wnud, recv_wnud, send_wnud, itagwnud
  use mem_grid,   only: mza
  use mem_nudge,  only: mwnud, jtab_wnud, itab_wnud

  implicit none

  real, optional, intent(in) :: rvara1(mza,mwnud)
  real, optional, intent(in) :: rvara2(mza,mwnud)
  real, optional, intent(in) :: rvara3(mza,mwnud)
  real, optional, intent(in) :: rvara4(mza,mwnud)
  real, optional, intent(in) :: rvara5(mza,mwnud)
  real, optional, intent(in) :: rvara6(mza,mwnud)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jsend, ivar
  integer :: j
  integer :: iwnud
  integer :: iwnudglobe

! Before we send anything, post the receives

  do jrecv = 1, nrecvs_wnud

     call MPI_Irecv(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,MPI_PACKED, &
                    recv_wnud(jrecv)%iremote,itagwnud,MPI_COMM_WORLD,         &
                    recv_wnud(jrecv)%irequest,ierr                            )
  enddo

! Now we can actually go on to sending the stuff

  do jsend = 1, nsends_wnud

     ipos = 0

     call MPI_Pack(jtab_wnud(jsend)%jend,1,MPI_INTEGER, &
          send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
     do j = 1,jtab_wnud(jsend)%jend
        iwnud = jtab_wnud(jsend)%iwnud(j)
        iwnudglobe = itab_wnud(iwnud)%iwnudglobe
!----------------------------------------------------------------
        
        call MPI_Pack(iwnudglobe,1,MPI_INTEGER, &
             send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Pack(rvara3(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Pack(rvara4(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara5)) then
           call MPI_Pack(rvara5(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara6)) then
           call MPI_Pack(rvara6(1,iwnud),mza,MPI_REAL, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Isend(send_wnud(jsend)%buff,ipos,MPI_PACKED, &
                    send_wnud(jsend)%iremote,itagwnud,MPI_COMM_WORLD, &
                    send_wnud(jsend)%irequest,ierr                    )
     
  enddo

#endif

end subroutine mpi_send_wnud

!=============================================================================

subroutine mpi_recv_v(mrl, rvara1, rvara2, rvara3, rvara4)

! Subroutine to perform a parallel MPI receive of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_para,   only: send_v, recv_v, nsends_v, nrecvs_v, mgroupsize
use mem_ijtabs, only: itabg_v, mrl_begs, istp, mloops
use mem_grid,   only: mza, mva
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: rvara1(mza,mva)
real, optional, intent(inout) :: rvara2(mza,mva)
real, optional, intent(inout) :: rvara3(mza,mva)
real, optional, intent(inout) :: rvara4(mza,mva)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nvpts
integer :: j
integer :: iv
integer :: ivglobe

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_v(mrl)

   call MPI_Waitany(nrecvs_v(mrl), recv_v(1:nrecvs_v(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
      nvpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
   do j = 1,nvpts
      call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
         ivglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iv = itabg_v(ivglobe)%iv_myrank
      if (iv < 2) iv = itabg_v(ivglobe)%iv_myrank_ivp
!----------------------------------------------------------------

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

   enddo

enddo

! Make sure all of our sends are finished and de-allocated

call MPI_Waitall(nsends_v(mrl), send_v(1:nsends_v(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_v

!=============================================================================

subroutine mpi_recv_m(mrl, rvara1, rvara2)

! Subroutine to perform a parallel MPI receive of a "M group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,   only: send_m, recv_m, nsends_m, nrecvs_m
  use mem_ijtabs, only: itabg_m
  use mem_grid,   only: mza, mma
  use misc_coms,  only: io6

  implicit none

  integer,        intent(in)    :: mrl
  real, optional, intent(inout) :: rvara1(mza,mma)
  real, optional, intent(inout) :: rvara2(mza,mma)

#ifdef OLAM_MPI

  integer :: ierr,ipos
  integer :: jrecv,jsend,ivar,jtmp
  integer :: nmpts
  integer :: j
  integer :: im
  integer :: imglobe

  if (mrl < 1) return

!  Now, let's wait on our receives

  do jtmp = 1, nrecvs_m(mrl)

     call MPI_Waitany(nrecvs_m(mrl), recv_m(1:nrecvs_m(mrl))%irequest, jrecv, &
                      MPI_STATUS_IGNORE, ierr)

     !  We got all our stuff.  Now unpack it into appropriate space.
     
     ipos = 0

     call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                     nmpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
     do j = 1, nmpts
        call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                        imglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        im = itabg_m(imglobe)%im_myrank
        if (im < 2) im = itabg_m(imglobe)%im_myrank_imp
!----------------------------------------------------------------

        if (present(rvara1)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
           rvara1(1,im),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

     enddo

  enddo

! Make sure all of our sends are finished and de-allocated

  call MPI_Waitall(nsends_m(mrl), send_m(1:nsends_m(mrl))%irequest, &
       MPI_STATUSES_IGNORE, ierr)

#endif

end subroutine mpi_recv_m

!=============================================================================

subroutine mpi_recv_w(mrl, scalars, dvara1, dvara2,                      &
                      rvara1, rvara2, rvara3, rvara4,  rvara5,  rvara6,  &
                      rvara7, rvara8, rvara9, rvara10, rvara11, rvara12, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5   )

! Subroutine to perform a parallel MPI receive of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use var_tables, only: vtab_r, nvar_par, num_scalar, nptonv
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w, mgroupsize
use mem_ijtabs, only: itabg_w, mrl_begs, mrl_begl, mrl_endl, istp, mloops
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use consts_coms,only: r8

implicit none

integer, intent(in) :: mrl

character(1), optional, intent(in) :: scalars

real(r8), optional, intent(inout) :: dvara1 (mza,mwa)
real(r8), optional, intent(inout) :: dvara2 (mza,mwa)

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

real, optional, intent(inout) :: r1dvara1(mwa)
real, optional, intent(inout) :: r1dvara2(mwa)
real, optional, intent(inout) :: r1dvara3(mwa)
real, optional, intent(inout) :: r1dvara4(mwa)
real, optional, intent(inout) :: r1dvara5(mwa)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwpts
integer :: i, j
integer :: iw
integer :: iwglobe

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_w(mrl)

   call MPI_Waitany(nrecvs_w(mrl), recv_w(1:nrecvs_w(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
      nwpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
   do j = 1,nwpts
      call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
         iwglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iw = itabg_w(iwglobe)%iw_myrank
      if (iw < 2) iw = itabg_w(iwglobe)%iw_myrank_iwp
!----------------------------------------------------------------

      if (present(dvara1)) then
         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
            dvara1(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
      endif

      if (present(dvara2)) then
         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
            dvara2(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
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

      if (present(scalars)) then
         do i = 1, nvar_par
            ivar = nptonv(i)

            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
               vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         enddo
      endif

   enddo

enddo

! Make sure all of our sends are finished and de-allocated

 call MPI_Waitall(nsends_w(mrl), send_w(1:nsends_w(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_w

!=============================================================================

subroutine mpi_recv_wnud(rvara1, rvara2, rvara3, rvara4, rvara5, rvara6)

! Subroutine to perform a parallel MPI receive of a "WNUD group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,  only: send_wnud, recv_wnud, nsends_wnud, nrecvs_wnud
  use mem_grid,  only: mza
  use mem_nudge, only: mwnud, itabg_wnud
  use misc_coms, only: io6

  implicit none

  real, optional, intent(inout) :: rvara1(mza,mwnud)
  real, optional, intent(inout) :: rvara2(mza,mwnud)
  real, optional, intent(inout) :: rvara3(mza,mwnud)
  real, optional, intent(inout) :: rvara4(mza,mwnud)
  real, optional, intent(inout) :: rvara5(mza,mwnud)
  real, optional, intent(inout) :: rvara6(mza,mwnud)

  real :: vctr1(mza)

#ifdef OLAM_MPI

  integer :: ierr,ipos
  integer :: jrecv,jsend,ivar,jtmp
  integer :: nwnudpts
  integer :: j
  integer :: iwnud
  integer :: iwnudglobe

!  Now, let's wait on our receives

  do jtmp = 1, nrecvs_wnud

     call MPI_Waitany(nrecvs_wnud, recv_wnud(1:nrecvs_wnud)%irequest, jrecv, &
                      MPI_STATUS_IGNORE, ierr)

     !  We got all our stuff.  Now unpack it into appropriate space.
     
     ipos = 0

     call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                     nwnudpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------
     do j = 1, nwnudpts
        call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                        iwnudglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        iwnud = itabg_wnud(iwnudglobe)%iwnud_myrank
!----------------------------------------------------------------

        if (present(rvara1)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara1(2:mza,iwnud) = rvara1(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara2(2:mza,iwnud) = rvara2(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(rvara3)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara3(2:mza,iwnud) = rvara3(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(rvara4)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara4(2:mza,iwnud) = rvara4(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(rvara5)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara5(2:mza,iwnud) = rvara5(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(rvara6)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           rvara6(2:mza,iwnud) = rvara6(2:mza,iwnud) + vctr1(2:mza)
        endif

     enddo

  enddo

! Make sure all of our sends are finished and de-allocated

  call MPI_Waitall(nsends_wnud, send_wnud(1:nsends_wnud)%irequest, &
       MPI_STATUSES_IGNORE, ierr)

#endif

end subroutine mpi_recv_wnud

!=============================================================================

subroutine olam_stop(message)

  use mem_para,  only: myrank

#ifdef OLAM_MPI
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

End Module olam_mpi_atm
