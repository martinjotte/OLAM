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
subroutine olam_alloc_mpi_land(mrls)

#ifdef OLAM_MPI
  use mpi
#endif

use mem_leaf,  only: jtab_wl_mpi
use mem_para,  only: nsends_wl, nsends_wlf, nrecvs_wl, nrecvs_wlf,  &
                     send_wl, send_wlf, recv_wl, recv_wlf
use mem_sflux, only: landflux, jlandflux
use misc_coms, only: io6

implicit none

integer, intent(in) :: mrls

#ifdef OLAM_MPI

integer :: nbytes_int
integer :: nbytes_real

integer :: nbytes_per_iwl
integer :: nbytes_per_iwlf

integer :: itag210 = 210
integer :: itag250 = 250

integer :: ierr
integer :: jsend
integer :: jrecv
integer :: jtmp

integer :: nwlpts, nwlfpts
integer :: mrl

! FOR LAND CELLS, THIS ROUTINE ASSUMES MRL=1. IF MPI DATA TRANSFERS OF LAND
! CELLS IS EVER NEEDED, MODIFY THE LAND PORTION OF THE CODE BY FOLLOWING WHAT
! IS DONE FOR THE LANDFLUX COMMUNICATION

integer :: lfsbuf(mrls)
integer :: lfrbuf(mrls, nrecvs_wlf(1))

! Allocate send buffers

call MPI_Pack_size(1,MPI_INTEGER,MPI_COMM_WORLD,nbytes_int  ,ierr)
call MPI_Pack_size(1,MPI_REAL   ,MPI_COMM_WORLD,nbytes_real ,ierr)

! Post land cell receives

do jrecv = 1,nrecvs_wl(1)

   ! Hardwired for 1 mrl only
   call MPI_IRecv(recv_wl(jrecv)%nbytes, 1, MPI_INTEGER, recv_wl(jrecv)%iremote, &
        itag210, MPI_COMM_WORLD, recv_wl(jrecv)%irequest, ierr)

enddo

! Post landflux cell receives

do jrecv = 1,nrecvs_wlf(1)

   call MPI_IRecv(lfrbuf(1,jrecv), mrls, MPI_INTEGER, recv_wlf(jrecv)%iremote, &
        itag250, MPI_COMM_WORLD, recv_wlf(jrecv)%irequest, ierr)

enddo

! Determine number of bytes to send per IWL column

nbytes_per_iwl = 1 * nbytes_int   &
               + 4 * nbytes_real

! Loop over all WL sends

do jsend = 1,nsends_wl(1)

! Determine size of send_wl buffer for mrl = 1

   send_wl(jsend)%nbytes = nbytes_int  &
                         + nbytes_per_iwl * jtab_wl_mpi(jsend)%jend(1)

! Allocate buffer

   allocate(send_wl(jsend)%buff(send_wl(jsend)%nbytes))

! Send buffer sizes to receive ranks

   ! Hardwired for mrl=1. If other mrls become necessary follow what is done
   ! in the following landflux section
   call MPI_Send(send_wl(jsend)%nbytes, 1, MPI_INTEGER, &
                 send_wl(jsend)%iremote, itag210, MPI_COMM_WORLD, ierr)

enddo

! Determine number of bytes to send per LANDFLUX cell

nbytes_per_iwlf = 1 * nbytes_int   &
                + 6 * nbytes_real

! Loop over all WLF sends for mrl = 1

do jsend = 1, nsends_wlf(1)

! Determine size of send_wlf buffer for mrl = 1

   send_wlf(jsend)%nbytes = nbytes_int  &
                          + nbytes_per_iwlf * jlandflux(2+jsend)%jend(1)

! Allocate buffer

   allocate(send_wlf(jsend)%buff(send_wlf(jsend)%nbytes))

! Send buffer sizes to receive ranks

   lfsbuf(1) = send_wlf(jsend)%nbytes

   do mrl = 2, mrls
      lfsbuf(mrl) = jlandflux(2+jsend)%jend(mrl)
      
! If at least 1 WLF point needs to be sent to current remote rank for 
! current mrl, increase nsends_wlf(mrl) by 1.

      if (jlandflux(2+jsend)%jend(mrl) > 0) nsends_wlf(mrl) = nsends_wlf(mrl) + 1
   enddo

   call MPI_Send(lfsbuf, mrls, MPI_INTEGER, send_wlf(jsend)%iremote, itag250, &
        MPI_COMM_WORLD, ierr)

enddo

! Loop over all WL receives

do jtmp = 1, nrecvs_wl(1)

! Get completed recv_wl buffer sizes

   ! Hardwired for mrl=1. If other mrls become necessary follow what is done
   ! in the following landflux section
   call MPI_Waitany(nrecvs_wl(1), recv_wl(1:nrecvs_wl(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

! Allocate recv_wl buffers for completed transfer

   allocate(recv_wl(jrecv)%buff(recv_wl(jrecv)%nbytes))

enddo

! Loop over all WLF receives

do jtmp = 1,nrecvs_wlf(1)

! Get completed recv_wlf buffer sizes

   call MPI_Waitany(nrecvs_wlf(1), recv_wlf(1:nrecvs_wlf(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

! Allocate recv_wlf buffers for completed transfer

   recv_wlf(jrecv)%nbytes = lfrbuf(1,jrecv)
   allocate(recv_wlf(jrecv)%buff(recv_wlf(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 WLF point needs to be 
! received from current remote rank for each mrl, increase nrecvs_wlf(mrl) by 1

   do mrl = 2,mrls
      if (lfrbuf(mrl,jrecv) > 0) nrecvs_wlf(mrl) = nrecvs_wlf(mrl) + 1
   enddo

enddo

#endif

return
end subroutine olam_alloc_mpi_land

!===============================================================================

subroutine mpi_send_wl(sendgroup)

! Subroutine to perform a parallel MPI send of a "WL group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_leaf,   only: land, itab_wl, jtab_wl_mpi

use mem_para,   only: nrecvs_wl, nsends_wl, send_wl, recv_wl

implicit none

character(1), intent(in) :: sendgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag4 = 4
integer :: j
integer :: iwl
integer :: iwlglobe

! Before we send anything, post the receives

do jrecv = 1,nrecvs_wl(1)

   call MPI_Irecv(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,MPI_PACKED,  &
                  recv_wl(jrecv)%iremote,itag4,MPI_COMM_WORLD,          &
                  recv_wl(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_wl(1)

   ipos = 0

   call MPI_Pack(jtab_wl_mpi(jsend)%jend(1),1,MPI_INTEGER,  &
      send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_wl_mpi(jsend)%jend(1)
      iwl = jtab_wl_mpi(jsend)%iwl(j)
      iwlglobe = itab_wl(iwl)%iwglobe
!----------------------------------------------------------------
      call qsub('WL',iwl)

      call MPI_Pack(iwlglobe,1,MPI_INTEGER,  &
         send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'A' .or. sendgroup == 'T') then

         call MPI_Pack(land%rough(iwl),1,MPI_INTEGER,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(land%can_temp(iwl),1,MPI_REAL,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(land%can_shv(iwl),1,MPI_REAL,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'R') then

         call MPI_Pack(land%rlongup(iwl),1,MPI_INTEGER,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(land%rlong_albedo(iwl),1,MPI_INTEGER,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(land%albedo_beam(iwl),1,MPI_INTEGER,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(land%albedo_diffuse(iwl),1,MPI_INTEGER,  &
            send_wl(jsend)%buff,send_wl(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WLsend',jsend)

   call MPI_Isend(send_wl(jsend)%buff,ipos,MPI_PACKED,        &
                  send_wl(jsend)%iremote,itag4,MPI_COMM_WORLD,  &
                  send_wl(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_wl

!===============================================================================

subroutine mpi_send_wlf(sendgroup,mrl)

! Subroutine to perform a parallel MPI send of a "WLF group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms,  only: io6
use mem_para,   only: nrecvs_wlf, nsends_wlf, send_wlf, recv_wlf
use mem_sflux,  only: landflux, jlandflux

implicit none

character(1), intent(in) :: sendgroup
integer,      intent(in) :: mrl

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag5 = 5
integer :: j
integer :: ilf
integer :: ilfglobe

real :: rscr(6)

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_wlf(mrl)

   call MPI_Irecv(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,MPI_PACKED,  &
                  recv_wlf(jrecv)%iremote,itag5,MPI_COMM_WORLD,          &
                  recv_wlf(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_wlf(mrl)

   ipos = 0

   call MPI_Pack(jlandflux(2+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jlandflux(2+jsend)%jend(mrl)
      ilf = jlandflux(2+jsend)%ilandflux(j)
      ilfglobe = landflux(ilf)%ifglobe
!----------------------------------------------------------------
      call qsub('WLF',ilf)

      call MPI_Pack(ilfglobe,1,MPI_INTEGER,  &
         send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'A') then ! for initialization

         rscr(1) = landflux(ilf)%rhos
         rscr(2) = landflux(ilf)%airtemp
         rscr(3) = landflux(ilf)%airshv

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'T') then ! for turbulent fluxes
      
         rscr(1) = landflux(ilf)%vels
         rscr(2) = landflux(ilf)%prss
         rscr(3) = landflux(ilf)%rhos
         rscr(4) = landflux(ilf)%sxfer_t
         rscr(5) = landflux(ilf)%sxfer_r
         rscr(6) = landflux(ilf)%ustar

         call MPI_Pack(rscr,6,MPI_REAL,  &
            send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'C') then ! for cuparm fluxes

         rscr(1) = landflux(ilf)%pcpg
         rscr(2) = landflux(ilf)%qpcpg
         rscr(3) = landflux(ilf)%dpcpg

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'M') then ! for microphysics precip fluxes

         rscr(1) = landflux(ilf)%pcpg
         rscr(2) = landflux(ilf)%qpcpg
         rscr(3) = landflux(ilf)%dpcpg

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'R') then ! for radiative fluxes

         rscr(1) = landflux(ilf)%rlong
         rscr(2) = landflux(ilf)%rshort
         rscr(3) = landflux(ilf)%rshort_diffuse

         call MPI_Pack(rscr,3,MPI_REAL,  &
            send_wlf(jsend)%buff,send_wlf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WLFsend',jsend)

   call MPI_Isend(send_wlf(jsend)%buff,ipos,MPI_PACKED,        &
                  send_wlf(jsend)%iremote,itag5,MPI_COMM_WORLD,  &
                  send_wlf(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_wlf

!=============================================================================

subroutine mpi_recv_wl(recvgroup)

! Subroutine to perform a parallel MPI receive of a "WL group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms, only: io6
use mem_leaf,  only: land, itabg_wl
use mem_para,  only: nsends_wl, nrecvs_wl, send_wl, recv_wl

implicit none

character(1), intent(in) :: recvgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwlpts
integer :: j
integer :: iwl
integer :: iwlglobe

! Now, let's wait on our receives

do jtmp = 1,nrecvs_wl(1)

   call MPI_Waitany(nrecvs_wl(1), recv_wl(1:nrecvs_wl(1))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
      nwlpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwlpts
      call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
         iwlglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iwl = itabg_wl(iwlglobe)%iwl_myrank
!----------------------------------------------------------------
      call qsub('WL',iwl)

      if (recvgroup == 'A' .or. recvgroup == 'T') then

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%rough(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%can_temp(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%can_shv(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'R') then

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%rlongup(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%rlong_albedo(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%albedo_beam(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_wl(jrecv)%buff,recv_wl(jrecv)%nbytes,ipos,  &
            land%albedo_diffuse(iwl),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('WLrecv',jrecv)

enddo

! Make sure sends are all finished and de-allocated
call MPI_Waitall(nsends_wl(1), send_wl(1:nsends_wl(1))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_wl

!=============================================================================

subroutine mpi_recv_wlf(recvgroup,mrl)

! Subroutine to perform a parallel MPI receive of a "WLF group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use misc_coms, only: io6
use mem_para,  only: nsends_wlf, nrecvs_wlf, send_wlf, recv_wlf
use mem_sflux, only: landflux, landfluxg

implicit none

character(1), intent(in) :: recvgroup
integer,      intent(in) :: mrl

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwlfpts
integer :: j
integer :: ilf
integer :: ilfglobe
real    :: rscr(6)

if (mrl < 1) return

! Now, let's wait on our receives

do jtmp = 1,nrecvs_wlf(mrl)

   call MPI_Waitany(nrecvs_wlf(mrl), recv_wlf(1:nrecvs_wlf(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
      nwlfpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwlfpts
      call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
         ilfglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      ilf = landfluxg(ilfglobe)%ilf_myrank
!----------------------------------------------------------------
      call qsub('WLF',ilf)

      if (recvgroup == 'A') then ! for initialization

         call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)
            
         landflux(ilf)%rhos    = rscr(1)
         landflux(ilf)%airtemp = rscr(2)
         landflux(ilf)%airshv  = rscr(3)

      elseif (recvgroup == 'T') then ! for turbulent fluxes

         call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
            rscr,6,MPI_REAL,MPI_COMM_WORLD,ierr)
            
         landflux(ilf)%vels    = rscr(1)
         landflux(ilf)%prss    = rscr(2)
         landflux(ilf)%rhos    = rscr(3)
         landflux(ilf)%sxfer_t = rscr(4)
         landflux(ilf)%sxfer_r = rscr(5)
         landflux(ilf)%ustar   = rscr(6)

      elseif (recvgroup == 'C') then ! for cuparm fluxes

         call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         landflux(ilf)%pcpg  = rscr(1)
         landflux(ilf)%qpcpg = rscr(2)
         landflux(ilf)%dpcpg = rscr(3)

      elseif (recvgroup == 'M') then ! for microphysics precip fluxes

         call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         landflux(ilf)%pcpg  = rscr(1)
         landflux(ilf)%qpcpg = rscr(2)
         landflux(ilf)%dpcpg = rscr(3)

      elseif (recvgroup == 'R') then ! for radiative fluxes

         call MPI_Unpack(recv_wlf(jrecv)%buff,recv_wlf(jrecv)%nbytes,ipos,  &
            rscr,3,MPI_REAL,MPI_COMM_WORLD,ierr)

         landflux(ilf)%rlong          = rscr(1)
         landflux(ilf)%rshort         = rscr(2)
         landflux(ilf)%rshort_diffuse = rscr(3)

      endif

   enddo
   call rsub('WLFrecv',jrecv)

enddo

! Make sure sends are all finished and de-allocated
call MPI_Waitall(nsends_wlf(mrl), send_wlf(1:nsends_wlf(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_wlf
