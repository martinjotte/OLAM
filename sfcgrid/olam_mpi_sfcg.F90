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

  use mem_sfcg,  only: itab_wsfc, itabg_wsfc, jtab_wsfc_mpi
  use mem_para,  only: nbytes_int, nbytes_real, ireqr_wsfc, ireqs_wsfc, nsends_wsfc, &
                       nrecvs_wsfc, send_wsfc, recv_wsfc, icurr_wsfc, inext_wsfc, itagwsfc
  implicit none

  integer, intent(in) :: nzg
  integer, intent(in) :: nzs

#ifdef OLAM_MPI

  integer, parameter :: itag410 = 410
  integer, parameter :: itag411 = 411

  integer :: nbytes_per_iwsfc
  integer :: ipos, ierr
  integer :: jsend, jrecv, jtmp, j

  integer, allocatable :: iwglobe(:)
! integer, allocatable :: ireqr_wsfc2(:)

  allocate( ireqr_wsfc(nrecvs_wsfc,2) )
  allocate( ireqs_wsfc(nsends_wsfc,2) )

! allocate( ireqr_wsfc2(nrecvs_wsfc) )

  ! Post SFC grid cell receives

  do jrecv = 1, nrecvs_wsfc

     ! Hardwired for 1 mrl only
     allocate( recv_wsfc(jrecv)%npts(1) )

     call MPI_IRecv(recv_wsfc(jrecv)%npts, 1, MPI_INTEGER, recv_wsfc(jrecv)%iremote, &
                    itag410, MPI_COMM_WORLD, ireqr_wsfc(jrecv,icurr_wsfc), ierr)
  enddo

  ! Determine number of bytes to send per IWSFC column

  nbytes_per_iwsfc = 3       * nbytes_int  &
                   + 18      * nbytes_real &
                   + 2 * nzg * nbytes_real &
                   + 3 * nzs * nbytes_real

  ! Loop over all WSFC sends


  do jsend = 1, nsends_wsfc

     ! Send buffer sizes to receive ranks
     ! Hardwired for mrl=1.

     call MPI_Isend(jtab_wsfc_mpi(jsend)%jend(1), 1, MPI_INTEGER, &
                    send_wsfc(jsend)%iremote, itag410, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,icurr_wsfc), ierr)

     ! Determine size of send_wsfc buffer for mrl = 1

     send_wsfc(jsend)%nbytes = nbytes_per_iwsfc * jtab_wsfc_mpi(jsend)%jend(1)

     ! Allocate buffer

     allocate(send_wsfc(jsend)%buff(send_wsfc(jsend)%nbytes))

  enddo

  ! Loop over all WSFC receives

  do jtmp = 1, nrecvs_wsfc

     ! Get completed recv_wsfc buffer sizes
     ! Hardwired for mrl=1.

     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     recv_wsfc(jrecv)%nbytes = nbytes_per_iwsfc * recv_wsfc(jrecv)%npts(1)

     ! Allocate recv_wsfc buffers for completed transfer

     allocate( recv_wsfc(jrecv)%buff( recv_wsfc(jrecv)%nbytes ) )
     allocate( recv_wsfc(jrecv)%ipts( recv_wsfc(jrecv)%npts(1) ) )

     ! Post next receive to get the list of iwsfc points we are getting

     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itag411, MPI_COMM_WORLD, &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr)
  enddo

  ! Communicate the list of iwsfc points we are sending to adjacent nodes

  do jtmp = 1, nsends_wsfc

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_wsfc, ireqs_wsfc(:,icurr_wsfc), jsend, MPI_STATUS_IGNORE, ierr)

     allocate(iwglobe( jtab_wsfc_mpi(jsend)%jend(1) ) )

     do j = 1, jtab_wsfc_mpi(jsend)%jend(1)
        iwglobe(j) = itab_wsfc( jtab_wsfc_mpi(jsend)%iwsfc(j) )%iwglobe
     enddo

     ipos = 0

     call MPI_Pack(iwglobe, jtab_wsfc_mpi(jsend)%jend(1), MPI_INTEGER, &
          send_wsfc(jsend)%buff, send_wsfc(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_wsfc(jsend)%buff, ipos, MPI_PACKED, &
                    send_wsfc(jsend)%iremote, itag411, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,inext_wsfc), ierr)

     deallocate(iwglobe)
  enddo

  ! increment MPI request pointers

  icurr_wsfc = mod(icurr_wsfc,2) + 1
  inext_wsfc = mod(inext_wsfc,2) + 1

  ! Unpack and store the list of iwsfc points we are getting from adjacent nodes

  do jtmp = 1, nrecvs_wsfc

     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iwglobe( recv_wsfc(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, ipos, &
                     iwglobe, recv_wsfc(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_wsfc(jrecv)%npts(1)
        recv_wsfc(jrecv)%ipts(j) = itabg_wsfc( iwglobe(j) )%iwsfc_myrank
     enddo

     deallocate(iwglobe)

     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr                          )
  enddo

!!deallocate(ireqr_wsfc2)
#endif

end subroutine olam_alloc_mpi_sfcg

!===============================================================================

subroutine mpi_send_wsfc(head, soil_watfrac)

  ! Subroutine to perform a parallel MPI send of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,    only: sfcg, itab_wsfc, jtab_wsfc_mpi, mwsfc
  use mem_para,    only: nsends_wsfc, send_wsfc, itagwsfc, ireqs_wsfc, &
                         icurr_wsfc, inext_wsfc
  use leaf_coms,   only: nzs
  use sea_coms,    only: nzi
  use mem_land,    only: nzg, land, omland
  use mem_lake,    only: lake, omlake
  use mem_sea,     only: sea, omsea

  implicit none

  real, optional, intent(inout) :: head        (nzg,mwsfc)
  real, optional, intent(inout) :: soil_watfrac(nzg,mwsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: j
  integer :: iwsfc, iland, ilake, isea
! integer :: iwglobe

  ! Before we send anything, post the receives

!  !$omp do private(ierr) schedule(static,1)
!  do jrecv = 1, nrecvs_wsfc
!
!     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
!                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
!                    ireqr_wsfc(jrecv), ierr                                     )
!
!  enddo
!  !$omp end do nowait

!!!$omp do private(ipos,ierr,j,iwsfc,iwglobe,iland,ilake,isea) schedule(dynamic,1)

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_wsfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_wsfc, ireqs_wsfc(:,icurr_wsfc), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,j,iwsfc,iland,ilake,isea,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0

     ! Loop over number of columns for this jsend

!----------------------------------------------------------------
     do j = 1,jtab_wsfc_mpi(jsend)%jend(1)
        iwsfc = jtab_wsfc_mpi(jsend)%iwsfc(j)
!       iwglobe = itab_wsfc(iwsfc)%iwglobe
!----------------------------------------------------------------

        ! Pack the messages into send buffers

        if (present(head)) then

           ! Pack HEAD, HYDRESIST

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(head        (1,iland),nzg,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(soil_watfrac(1,iland),nzg,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           endif

        else

           ! Pack SFCG quantities

           call MPI_Pack(sfcg%vels    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%prss    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rhos    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtemp (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtheta(iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airrrv  (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%cantemp (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%canrrv  (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rough   (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%head1   (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%wthv    (iwsfc),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              ! Pack LAND quantities

              call MPI_Pack(land%nlev_sfcwater    (iland),1,MPI_INTEGER,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd             (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd_diffuse     (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%cosz             (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_fracarea     (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_albedo       (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_rough        (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_temp         (iland),1,   MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_mass  (1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_energy(1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_depth (1,iland),nzs, MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_water     (1,iland),nzg, MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_energy    (1,iland),nzg, MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              ! Pack LAKE quantities

              call MPI_Pack(lake%lake_energy(ilake),1,MPI_REAL,send_wsfc(jsend)%buff, &
                            send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              ! Pack SEA quantities

              call MPI_Pack(sea%nlev_seaice (isea),1,MPI_INTEGER, &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_cantemp (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_cantemp (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_canrrv  (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_canrrv  (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_wthv    (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_wthv    (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_rough   (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_rough   (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaicec     (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seatc       (isea),1,MPI_REAL,    &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaice_tempk(1,isea),nzi,MPI_REAL, &
                           send_wsfc(jsend)%buff,send_wsfc(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

          endif

        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_wsfc(jsend)%buff, ipos, MPI_PACKED,            &
                    send_wsfc(jsend)%iremote, itagwsfc, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,inext_wsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_wsfc = mod(icurr_wsfc,2) + 1
  inext_wsfc = mod(inext_wsfc,2) + 1

#endif

end subroutine mpi_send_wsfc

!=============================================================================

subroutine mpi_recv_wsfc(head, soil_watfrac)

  ! Subroutine to perform a parallel MPI receive of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,    only: sfcg, mwsfc
  use leaf_coms,   only: nzs
  use sea_coms,    only: nzi
  use mem_land,    only: nzg, land, omland
  use mem_lake,    only: lake, omlake
  use mem_sea,     only: sea, omsea
  use mem_para,    only: nrecvs_wsfc, itagwsfc, recv_wsfc, ireqr_wsfc, icurr_wsfc, inext_wsfc

  implicit none

  real, optional, intent(inout) :: head        (nzg,mwsfc)
  real, optional, intent(inout) :: soil_watfrac(nzg,mwsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: iwsfc, iland, ilake, isea
! integer :: iwglobe

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_wsfc

     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(j,ipos,iwsfc,iland,ilake,isea,ierr) &
     !$omp      firstprivate(jrecv) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

!----------------------------------------------------------------------
     do j = 1, recv_wsfc(jrecv)%npts(1)
        iwsfc = recv_wsfc(jrecv)%ipts(j)
!----------------------------------------------------------------------

        if (present(head)) then

           ! Unpack HEAD, SOIL_WATFRAC

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              head        (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              soil_watfrac(1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
           endif

        else

           ! Unpack SFCG quantities

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%vels    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%prss    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%rhos    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%airtemp (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%airtheta(iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%airrrv  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%cantemp (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%canrrv  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%rough   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%head1   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%wthv    (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              ! Unpack LAND quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%nlev_sfcwater    (iland),1,  MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd             (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd_diffuse     (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%cosz             (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_fracarea     (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_albedo       (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_rough        (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_temp         (iland),1,  MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_mass  (1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_energy(1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_depth (1,iland),nzs,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_water     (1,iland),nzg,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_energy    (1,iland),nzg,MPI_REAL,   MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              ! Unpack LAKE quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              lake%lake_energy(ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              ! Unpack SEA quantities

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%nlev_seaice (isea),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_cantemp (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_cantemp (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_canrrv  (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_canrrv  (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_wthv    (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_wthv    (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_rough   (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_rough   (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaicec     (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seatc       (isea),1,MPI_REAL,   MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaice_tempk(1,isea),nzi,MPI_REAL,   MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo

     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_wsfc

End Module olam_mpi_sfcg
