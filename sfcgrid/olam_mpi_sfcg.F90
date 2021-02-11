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
Module olam_mpi_sfcg

Contains

subroutine olam_alloc_mpi_sfcg(nzg, nzs)

  ! FOR SFC GRID CELLS, THIS ROUTINE ASSUMES MRL=1.

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,  only: jtab_wsfc_mpi
  use mem_para,  only: nbytes_int, nbytes_real, nbytes_real8,           &
                       nsends_wsfc, nrecvs_wsfc, send_wsfc, recv_wsfc
  use misc_coms, only: io6

  implicit none

  integer, intent(in) :: nzg
  integer, intent(in) :: nzs

#ifdef OLAM_MPI

  integer :: nbytes_per_iwsfc

  integer :: itag410 = 410

  integer :: ierr
  integer :: jsend
  integer :: jrecv
  integer :: jtmp

  integer, allocatable :: ireqs(:)

  integer              :: wsbuf(1)
  integer, allocatable :: wrbuf(:,:)

  ! Post SFC grid cell receives

  allocate(wrbuf(1,nrecvs_wsfc(1)))

  do jrecv = 1,nrecvs_wsfc(1)

     ! Hardwired for 1 mrl only
     call MPI_IRecv(wrbuf(1,jrecv), 1, MPI_INTEGER, recv_wsfc(jrecv)%iremote, &
          itag410, MPI_COMM_WORLD, recv_wsfc(jrecv)%irequest, ierr)

  enddo

  ! Determine number of bytes to send per IWSFC column

  nbytes_per_iwsfc = 3       * nbytes_int  &
                   + 18      * nbytes_real &
                   + 2 * nzg * nbytes_real &
                   + 3 * nzs * nbytes_real

  ! Loop over all WSFC sends

  do jsend = 1,nsends_wsfc(1)

  ! Determine size of send_wsfc buffer for mrl = 1

     send_wsfc(jsend)%nbytes = nbytes_int  &
                             + nbytes_per_iwsfc * jtab_wsfc_mpi(jsend)%jend(1)

  ! Allocate buffer

     allocate(send_wsfc(jsend)%buff(send_wsfc(jsend)%nbytes))

  ! Send buffer sizes to receive ranks

     wsbuf(1) = send_wsfc(jsend)%nbytes

     ! Hardwired for mrl=1.
     call MPI_Send(wsbuf, 1, MPI_INTEGER, send_wsfc(jsend)%iremote, itag410, &
          MPI_COMM_WORLD, ierr)

  enddo

  allocate( ireqs(nrecvs_wsfc(1)))

  ireqs(1:nrecvs_wsfc(1)) = recv_wsfc(1:nrecvs_wsfc(1))%irequest

  ! Loop over all WSFC receives

  do jtmp = 1, nrecvs_wsfc(1)

     ! Get completed recv_wsfc buffer sizes

     ! Hardwired for mrl=1.
     call MPI_Waitany(nrecvs_wsfc(1), ireqs, jrecv, MPI_STATUS_IGNORE, ierr)

     recv_wsfc(jrecv)%nbytes = wrbuf(1,jrecv)

  ! Allocate recv_wsfc buffers for completed transfer

     allocate(recv_wsfc(jrecv)%buff(recv_wsfc(jrecv)%nbytes))

  enddo

#endif

end subroutine olam_alloc_mpi_sfcg

!===============================================================================

subroutine mpi_send_wsfc(pcpg, qpcpg, dpcpg, head, soil_watfrac)

  ! Subroutine to perform a parallel MPI send of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use misc_coms,   only: io6
  use consts_coms, only: r8
  use mem_sfcg,    only: sfcg, itab_wsfc, jtab_wsfc_mpi, mwsfc
  use mem_para,    only: nrecvs_wsfc, nsends_wsfc, send_wsfc, recv_wsfc, itagwsfc, nbytes_int
  use leaf_coms,   only: nzs
  use sea_coms,    only: nzi
  use mem_land,    only: nzg, land, omland
  use mem_lake,    only: lake, omlake
  use mem_sea,     only: sea, omsea

  implicit none

  real, optional, intent(inout) :: pcpg (8,mwsfc)
  real, optional, intent(inout) :: qpcpg(8,mwsfc)
  real, optional, intent(inout) :: dpcpg(8,mwsfc)
  real, optional, intent(inout) :: head        (nzg,mwsfc)
  real, optional, intent(inout) :: soil_watfrac(nzg,mwsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos, ipos0
  integer :: jrecv,jsend
  integer :: j, jv
  integer :: iwsfc, iland, ilake, isea
  integer :: iwglobe
  integer :: iposs(nsends_wsfc(1))

  ! Before we send anything, post the receives

  do jrecv = 1,nrecvs_wsfc(1)

     call MPI_Irecv(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,MPI_PACKED,  &
                    recv_wsfc(jrecv)%iremote,itagwsfc,MPI_COMM_WORLD,          &
                    recv_wsfc(jrecv)%irequest,ierr                          )

  enddo

  ! Make sure previous sends are finished

  do jsend = 1, nsends_wsfc(1)

     if (send_wsfc(jsend)%nbytes > 0) then
        call MPI_Wait( send_wsfc(jsend)%irequest, MPI_STATUS_IGNORE, ierr)
    endif
  enddo

  ! Pack the messages into send buffers

!s  !$omp parallel do private(ipos,ierr,j,iwsfc,iland,ilake,isea,iwglobe,ipos0) shared(iposs)
  do jsend = 1,nsends_wsfc(1)

     ipos = 0

     ! Pack number of columns for this jsend

     call MPI_Pack(jtab_wsfc_mpi(jsend)%jend(1),1,MPI_INTEGER,  &
        send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

     ! reserve space for starting positions of all columns for this jsend

     ipos = ipos + nbytes_int * jtab_wsfc_mpi(jsend)%jend(1)

     ! Loop over number of columns for this jsend 

     do j = 1,jtab_wsfc_mpi(jsend)%jend(1)
        iwsfc = jtab_wsfc_mpi(jsend)%iwsfc(j)
        iwglobe = itab_wsfc(iwsfc)%iwglobe

        ! Pack starting position (data value = ipos) for this column, placing
        ! it at reserved position ipos0 in the buffer

        ipos0 = nbytes_int * j

        call MPI_Pack(ipos,1,MPI_INTEGER, &
             send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos0,MPI_COMM_WORLD,ierr)

        ! Pack sfcg global W index for this column, resuming usage of ipos value

        call MPI_Pack(iwglobe,1,MPI_INTEGER,  &
           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

        if (present(pcpg)) then

           ! Pack PCPG, QPCPG, DPCPG

           if (.not. present(qpcpg) .or. .not. present(dpcpg)) then
              write(io6,'(a)') 'In mpi_send_wsfc, pcpg is present but qpcpg and/or dpcpg are not'
              stop 'stop mpi_send_wsfc '
           endif

           call MPI_Pack(pcpg (1,iwsfc),8,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(qpcpg(1,iwsfc),8,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(dpcpg(1,iwsfc),8,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
 
        elseif (present(head)) then

           ! Pack HEAD, HYDRESIST

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(head        (1,iland),nzg,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(soil_watfrac(1,iland),nzg,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           endif

        else

           ! Pack SFCG quantities

           call MPI_Pack(sfcg%vels    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%prss    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rhos    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtemp (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtheta(iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airrrv  (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%cantemp (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%canrrv  (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rough   (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%head1   (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%wthv    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              ! Pack LAND quantities

              call MPI_Pack(land%nlev_sfcwater    (iland),1,MPI_INTEGER,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd             (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd_diffuse     (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%cosz             (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_fracarea     (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_albedo       (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_rough        (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_temp         (iland),1,   MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_mass  (1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_energy(1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_depth (1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_water     (1,iland),nzg, MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_energy    (1,iland),nzg, MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              ! Pack LAKE quantities

              call MPI_Pack(lake%lake_energy(ilake),1,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              ! Pack SEA quantities

              call MPI_Pack(sea%nlev_seaice (isea),1,MPI_INTEGER,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_cantemp (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_cantemp (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_canrrv  (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_canrrv  (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_wthv    (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_wthv    (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_rough   (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_rough   (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaicec     (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seatc       (isea),1,MPI_REAL,   send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaice_tempk(1,isea),nzi,MPI_REAL,send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

          endif

        endif

     enddo

     ! Save total number of packed bytes for this jsend

     iposs(jsend) = ipos

  enddo
!s  !$omp end parallel do

  ! Now we can actually go on to sending the stuff

  do jsend = 1, nsends_wsfc(1)

     call MPI_Isend(send_wsfc(jsend)%buff,iposs(jsend),MPI_PACKED,        &
                    send_wsfc(jsend)%iremote,itagwsfc,MPI_COMM_WORLD,  &
                    send_wsfc(jsend)%irequest,ierr                  )

  enddo

#endif

end subroutine mpi_send_wsfc

!=============================================================================

subroutine mpi_recv_wsfc(pcpg, qpcpg, dpcpg, head, soil_watfrac)

  ! Subroutine to perform a parallel MPI receive of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use misc_coms,   only: io6
  use consts_coms, only: r8
  use mem_sfcg,    only: sfcg, itabg_wsfc, mwsfc, itab_wsfc
  use leaf_coms,   only: nzs
  use sea_coms,    only: nzi
  use mem_land,    only: nzg, land, omland
  use mem_lake,    only: lake, omlake
  use mem_sea,     only: sea, omsea
  use mem_para,    only: nsends_wsfc, nrecvs_wsfc, send_wsfc, recv_wsfc, nbytes_int

  implicit none

  real, optional, intent(inout) :: pcpg (8,mwsfc)
  real, optional, intent(inout) :: qpcpg(8,mwsfc)
  real, optional, intent(inout) :: dpcpg(8,mwsfc)
  real, optional, intent(inout) :: head        (nzg,mwsfc)
  real, optional, intent(inout) :: soil_watfrac(nzg,mwsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos, ipos0
  integer :: jrecv,jtmp
  integer :: jend
  integer :: j, jv
  integer :: iwsfc, iland, ilake, isea
  integer :: iwglobe
  integer :: ireqs(nrecvs_wsfc(1))

  ! Now, let's wait on our receives

  ireqs(1:nrecvs_wsfc(1)) = recv_wsfc(1:nrecvs_wsfc(1))%irequest

  do jtmp = 1,nrecvs_wsfc(1)

     call MPI_Waitany(nrecvs_wsfc(1), ireqs, jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got all our stuff.  Now unpack it into appropriate space.

     ipos = 0

     ! Unpack number of columns for this jtmp/jrecv

     call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,  &
        jend,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

     ! Loop over number of columns for this jtmp/jrecv

!s   !$omp parallel do private(ipos,ipos0,iwglobe,iwsfc,iland,ilake,isea,ierr)
     do j = 1,jend

        ! Unpack starting position (data value = ipos) for this j column, accessing
        ! it at reserved position ipos0 in the buffer

        ipos0 = nbytes_int * j

        call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos0, &
           ipos,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        ! Unpack sfcg global W index for this column, using just-unpacked ipos value

        call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,  &
           iwglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        iwsfc = itabg_wsfc(iwglobe)%iwsfc_myrank

        if (present(pcpg)) then

           ! Unpack PCPG, QPCPG, DPCPG

           if (.not. present(qpcpg) .or. .not. present(dpcpg)) then
              write(io6,'(a)') 'In mpi_recv_wsfc, pcpg is present but qpcpg and/or dpcpg are not'
              stop 'stop mpi_recv_wsfc '
           endif

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,pcpg (1,iwsfc),8,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,qpcpg(1,iwsfc),8,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,dpcpg(1,iwsfc),8,MPI_REAL,MPI_COMM_WORLD,ierr)

        elseif (present(head)) then

           ! Pack HEAD, SOIL_WATFRAC

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,head        (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,soil_watfrac(1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
           endif

        else
 
           ! Unpack SFCG quantities

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%vels    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%prss    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%rhos    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%airtemp (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%airtheta(iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%airrrv  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%cantemp (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%canrrv  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%rough   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%head1   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sfcg%wthv    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              ! Unpack LAND quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%nlev_sfcwater    (iland),1,  MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%ppfd             (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%ppfd_diffuse     (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%cosz             (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%veg_fracarea     (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%veg_albedo       (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%veg_rough        (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%veg_temp         (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%sfcwater_mass  (1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%sfcwater_energy(1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%sfcwater_depth (1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%soil_water     (1,iland),nzg,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,land%soil_energy    (1,iland),nzg,MPI_REAL,   MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              ! Unpack LAKE quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,lake%lake_energy(ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              ! Unpack SEA quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%nlev_seaice (isea),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%sea_cantemp (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%ice_cantemp (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%sea_canrrv  (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%ice_canrrv  (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%sea_wthv    (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%ice_wthv    (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%sea_rough   (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%ice_rough   (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%seaicec     (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%seatc       (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos,sea%seaice_tempk(1,isea),nzi,MPI_REAL,   MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo
!s   !$omp end parallel do

  enddo

#endif

end subroutine mpi_recv_wsfc

End Module olam_mpi_sfcg

