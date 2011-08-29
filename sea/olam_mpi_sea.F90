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
subroutine olam_alloc_mpi_sea(mrls)

#ifdef OLAM_MPI
  use mpi
#endif

use mem_sea,    only: jtab_ws_mpi
use mem_para,   only: nsends_ws, nsends_wsf, nrecvs_ws, nrecvs_wsf,  &
                      send_ws, send_wsf, recv_ws, recv_wsf
use mem_sflux,  only: seaflux, jseaflux
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrls

#ifdef OLAM_MPI

integer :: nbytes_int
integer :: nbytes_real

integer :: nbytes_per_iws
integer :: nbytes_per_iwsf

integer :: itag310 = 310
integer :: itag350 = 350

integer :: ierr
integer :: jsend
integer :: jrecv
integer :: jtmp

integer :: nwspts, nwsfpts
integer :: mrl

! FOR SEA CELLS, THIS ROUTINE ASSUMES MRL=1. IF MPI DATA TRANSFERS OF SEA
! CELLS IS EVER NEEDED, MODIFY THE SEA PORTION OF THE CODE BY FOLLOWING WHAT
! IS DONE FOR THE SEAFLUX COMMUNICATION

integer :: sfsbuf(mrls)
integer :: sfrbuf(mrls, nrecvs_wsf(1))

! allocate send buffers

call MPI_Pack_size(1,MPI_INTEGER,MPI_COMM_WORLD,nbytes_int  ,ierr)
call MPI_Pack_size(1,MPI_REAL   ,MPI_COMM_WORLD,nbytes_real ,ierr)

! Post sea cell receives

do jrecv = 1,nrecvs_ws(1)

   ! Hardwired for 1 mrl only
   call MPI_IRecv(recv_ws(jrecv)%nbytes, 1, MPI_INTEGER, recv_ws(jrecv)%iremote, &
        itag310, MPI_COMM_WORLD, recv_ws(jrecv)%irequest, ierr)

enddo

! Post seaflux cell receives

do jrecv = 1,nrecvs_wsf(1)

   call MPI_IRecv(sfrbuf(1,jrecv), mrls, MPI_INTEGER, recv_wsf(jrecv)%iremote, &
        itag350, MPI_COMM_WORLD, recv_wsf(jrecv)%irequest, ierr)

enddo

! Determine number of bytes to send per IWS column

nbytes_per_iws = 1 * nbytes_int   &
               + 4 * nbytes_real

! Loop over all WS sends

do jsend = 1,nsends_ws(1)

! Determine size of send_ws buffer for mrl = 1

   send_ws(jsend)%nbytes = nbytes_int  &
                         + nbytes_per_iws * jtab_ws_mpi(jsend)%jend(1)

! Allocate buffer

   allocate(send_ws(jsend)%buff(send_ws(jsend)%nbytes))

! Send buffer sizes to receive ranks

   ! Hardwired for mrl=1. If other mrls become necessary follow what is done
   ! in the following landflux section
   call MPI_Send(send_ws(jsend)%nbytes, 1, MPI_INTEGER, &
                 send_ws(jsend)%iremote, itag310, MPI_COMM_WORLD, ierr)

enddo

! Determine number of bytes to send per SEAFLUX cell

nbytes_per_iwsf = 1 * nbytes_int   &
                + 4 * nbytes_real

! Loop over all WSF sends for mrl = 1

do jsend = 1, nsends_wsf(1)

! Determine size of send_wsf buffer for mrl = 1

   send_wsf(jsend)%nbytes = nbytes_int  &
                          + nbytes_per_iwsf * jseaflux(2+jsend)%jend(1)

! Allocate buffer

   allocate(send_wsf(jsend)%buff(send_wsf(jsend)%nbytes))

! Send buffer sizes to receive ranks

   sfsbuf(1) = send_wsf(jsend)%nbytes

   do mrl = 2, mrls
      sfsbuf(mrl) = jseaflux(2+jsend)%jend(mrl)
      
! If at least 1 WSF point needs to be sent to current remote rank for 
! current mrl, increase nsends_wsf(mrl) by 1.

      if (jseaflux(2+jsend)%jend(mrl) > 0) nsends_wsf(mrl) = nsends_wsf(mrl) + 1
   enddo

   call MPI_Send(sfsbuf, mrls, MPI_INTEGER, send_wsf(jsend)%iremote, itag350, &
        MPI_COMM_WORLD, ierr)

enddo

! Loop over all WS receives

do jtmp = 1, nrecvs_ws(1)

! Get completed recv_ws buffer sizes

   ! Hardwired for mrl=1. If other mrls become necessary follow what is done
   ! in the following seaflux section
   call MPI_Waitany(nrecvs_ws(1), recv_ws(1:nrecvs_ws(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

! Allocate recv_ws buffers for completed transfer

   allocate(recv_ws(jrecv)%buff(recv_ws(jrecv)%nbytes))

enddo

! Loop over all WSF receives

do jtmp = 1,nrecvs_wsf(1)

! Get completed recv_wsf buffer sizes

   call MPI_Waitany(nrecvs_wsf(1), recv_wsf(1:nrecvs_wsf(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

! Allocate recv_wsf buffers for completed transfer

   recv_wsf(jrecv)%nbytes = sfrbuf(1,jrecv)
   allocate(recv_wsf(jrecv)%buff(recv_wsf(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 WSF point needs to be 
! received from current remote rank for each mrl, increase nrecvs_wsf(mrl) by 1
   
   do mrl = 2,mrls
      if (sfrbuf(mrl,jrecv) > 0) nrecvs_wsf(mrl) = nrecvs_wsf(mrl) + 1
   enddo

enddo

#endif

return
end subroutine olam_alloc_mpi_sea

!===============================================================================

subroutine mpi_send_ws(sendgroup)

! Subroutine to perform a parallel MPI send of a "WS group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_sea,    only: sea, itab_ws, jtab_ws_mpi
use mem_para,   only: nrecvs_ws, nsends_ws, send_ws, recv_ws

implicit none

character(1), intent(in) :: sendgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag6 = 6
integer :: j
integer :: iws
integer :: iwsglobe

! Before we send anything, post the receives

do jrecv = 1,nrecvs_ws(1)

   call MPI_Irecv(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,MPI_PACKED,  &
                  recv_ws(jrecv)%iremote,itag6,MPI_COMM_WORLD,          &
                  recv_ws(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_ws(1)

   ipos = 0

   call MPI_Pack(jtab_ws_mpi(jsend)%jend(1),1,MPI_INTEGER,  &
      send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_ws_mpi(jsend)%jend(1)
      iws = jtab_ws_mpi(jsend)%iws(j)
      iwsglobe = itab_ws(iws)%iwglobe
!----------------------------------------------------------------
      call qsub('WS',iws)

      call MPI_Pack(iwsglobe,1,MPI_INTEGER,  &
         send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'A' .or. sendgroup == 'T') then

         call MPI_Pack(sea%rough(iws),1,MPI_INTEGER,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(sea%can_temp(iws),1,MPI_REAL,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(sea%can_shv(iws),1,MPI_REAL,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'R') then

         call MPI_Pack(sea%rlongup(iws),1,MPI_INTEGER,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(sea%rlong_albedo(iws),1,MPI_INTEGER,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(sea%albedo_beam(iws),1,MPI_INTEGER,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(sea%albedo_diffuse(iws),1,MPI_INTEGER,  &
            send_ws(jsend)%buff,send_ws(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WSsend',jsend)

   call MPI_Isend(send_ws(jsend)%buff,ipos,MPI_PACKED,        &
                  send_ws(jsend)%iremote,itag6,MPI_COMM_WORLD,  &
                  send_ws(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_ws

!===============================================================================

subroutine mpi_send_wsf(sendgroup,mrl)

! Subroutine to perform a parallel MPI send of a "WSF group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_para,   only: nrecvs_wsf, nsends_wsf, send_wsf, recv_wsf
use mem_sflux,  only: seaflux, jseaflux

implicit none

character(1), intent(in) :: sendgroup
integer,      intent(in) :: mrl

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag7 = 7
integer :: j
integer :: isf
integer :: isfglobe

real :: rscr(4)

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_wsf(mrl)

   call MPI_Irecv(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,MPI_PACKED,  &
                  recv_wsf(jrecv)%iremote,itag7,MPI_COMM_WORLD,          &
                  recv_wsf(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_wsf(mrl)

   ipos = 0

   call MPI_Pack(jseaflux(2+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jseaflux(2+jsend)%jend(mrl)
      isf = jseaflux(2+jsend)%iseaflux(j)
      isfglobe = seaflux(isf)%ifglobe
!----------------------------------------------------------------
      call qsub('WSF',isf)

      call MPI_Pack(isfglobe,1,MPI_INTEGER,  &
         send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'A') then ! for initialization

         rscr(1) = seaflux(isf)%rhos
         rscr(2) = seaflux(isf)%airtemp
         rscr(3) = seaflux(isf)%airshv

         call MPI_Pack(rscr,3,MPI_REAL,  &
         send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'T') then ! for turbulent fluxes

         rscr(1) = seaflux(isf)%rhos
         rscr(2) = seaflux(isf)%sxfer_t
         rscr(3) = seaflux(isf)%sxfer_r
         rscr(4) = seaflux(isf)%ustar

         call MPI_Pack(rscr,4,MPI_REAL,  &
            send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'C') then ! for cuparm fluxes

         rscr(1) = seaflux(isf)%pcpg
         rscr(2) = seaflux(isf)%qpcpg
         rscr(3) = seaflux(isf)%dpcpg

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'M') then ! for microphysics precip fluxes

         rscr(1) = seaflux(isf)%pcpg
         rscr(2) = seaflux(isf)%qpcpg
         rscr(3) = seaflux(isf)%dpcpg

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'R') then ! for radiative fluxes

         rscr(1) = seaflux(isf)%rlong
         rscr(2) = seaflux(isf)%rshort
         rscr(3) = seaflux(isf)%rshort_diffuse

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wsf(jsend)%buff,send_wsf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WSFsend',jsend)

   call MPI_Isend(send_wsf(jsend)%buff,ipos,MPI_PACKED,        &
                  send_wsf(jsend)%iremote,itag7,MPI_COMM_WORLD,  &
                  send_wsf(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_wsf

!=============================================================================

subroutine mpi_recv_ws(recvgroup)

! Subroutine to perform a parallel MPI receive of a "WS group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_sea,    only: sea, itabg_ws
use mem_para,   only: nsends_ws, nrecvs_ws, send_ws, recv_ws, myrank

implicit none

character(1), intent(in) :: recvgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: nwspts
integer :: j, jtmp
integer :: iws
integer :: iwsglobe

! Now, let's wait on our receives

do jtmp = 1,nrecvs_ws(1)

   call MPI_Waitany(nrecvs_ws(1), recv_ws(1:nrecvs_ws(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
      nwspts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwspts
      call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
         iwsglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iws = itabg_ws(iwsglobe)%iws_myrank
!----------------------------------------------------------------
      call qsub('WS',iws)

      if (recvgroup == 'A' .or. recvgroup == 'T') then

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%rough(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%can_temp(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%can_shv(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'R') then

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%rlongup(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%rlong_albedo(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%albedo_beam(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_ws(jrecv)%buff,recv_ws(jrecv)%nbytes,ipos,  &
            sea%albedo_diffuse(iws),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WSrecv',jrecv)

enddo

! Make sure sends are all finished and de-allocated
call MPI_Waitall(nsends_ws(1), send_ws(1:nsends_ws(1))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_ws

!=============================================================================

subroutine mpi_recv_wsf(recvgroup,mrl)

! Subroutine to perform a parallel MPI receive of a "WSF group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_sea,    only: sea, itabg_ws
use mem_para,   only: nsends_wsf, nrecvs_wsf, send_wsf, recv_wsf
use mem_sflux,  only: seaflux, seafluxg

implicit none

character(1), intent(in) :: recvgroup
integer,      intent(in) :: mrl

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwsfpts
integer :: j
integer :: isf
integer :: isfglobe
real    :: rscr(4)

if (mrl < 1) return

! Now, let's wait on our receives

do jtmp = 1,nrecvs_wsf(mrl)

   call MPI_Waitany(nrecvs_wsf(mrl), recv_wsf(1:nrecvs_wsf(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
      nwsfpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwsfpts
      call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
         isfglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      isf = seafluxg(isfglobe)%isf_myrank
!----------------------------------------------------------------
      call qsub('WSF',isf)
      
      if (recvgroup == 'A') then ! for initialization

         call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)
         
         seaflux(isf)%rhos    = rscr(1)
         seaflux(isf)%airtemp = rscr(2)
         seaflux(isf)%airshv  = rscr(3)

      elseif (recvgroup == 'T') then ! for turbulent fluxes

         call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
            rscr,4,MPI_REAL,MPI_COMM_WORLD,ierr)

         seaflux(isf)%rhos    = rscr(1)
         seaflux(isf)%sxfer_t = rscr(2)
         seaflux(isf)%sxfer_r = rscr(3)
         seaflux(isf)%ustar   = rscr(4)

      elseif (recvgroup == 'C') then ! for cuparm fluxes

         call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         seaflux(isf)%pcpg  = rscr(1)
         seaflux(isf)%qpcpg = rscr(2)
         seaflux(isf)%dpcpg = rscr(3)

      elseif (recvgroup == 'M') then ! for microphysics precip fluxes

         call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         seaflux(isf)%pcpg  = rscr(1)
         seaflux(isf)%qpcpg = rscr(2)
         seaflux(isf)%dpcpg = rscr(3)

      elseif (recvgroup == 'R') then ! for radiative fluxes

         call MPI_Unpack(recv_wsf(jrecv)%buff,recv_wsf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         seaflux(isf)%rlong          = rscr(1)
         seaflux(isf)%rshort         = rscr(2)
         seaflux(isf)%rshort_diffuse = rscr(3)

      endif

   enddo
   call rsub('WSFrecv',jrecv)

enddo

! Make sure sends are all finished and de-allocated
call MPI_Waitall(nsends_wsf(mrl), send_wsf(1:nsends_wsf(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_wsf
