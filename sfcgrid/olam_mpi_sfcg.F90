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

  use mem_sfcg,  only: itab_wsfc, itabg_wsfc, jtab_wsfc_mpi, &
                       itab_vsfc, itabg_vsfc, jtab_vsfc_mpi
  use mem_para,  only: nbytes_int, nbytes_real, &
                        ireqr_wsfc,  ireqs_wsfc,  ireqr_vsfc,  ireqs_vsfc, &
                       nsends_wsfc, nrecvs_wsfc, nsends_vsfc, nrecvs_vsfc, &
                         send_wsfc,   recv_wsfc,   send_vsfc,   recv_vsfc, &
                        icurr_wsfc,  inext_wsfc,  icurr_vsfc,  inext_vsfc, &
                          itagwsfc,                 itagvsfc
  implicit none

  integer, intent(in) :: nzg
  integer, intent(in) :: nzs

#ifdef OLAM_MPI

  integer :: nbytes_per_iwsfc
  integer :: nbytes_per_ivsfc

  integer, parameter :: itag410 = 410
  integer, parameter :: itag411 = 411
  integer, parameter :: itag510 = 510
  integer, parameter :: itag511 = 511

  integer :: ipos, ierr
  integer :: jsend, jrecv, jtmp, j

  integer, allocatable :: iglobe(:)

  ! Post WSFC receives

  allocate( ireqr_wsfc(nrecvs_wsfc,2) )
  allocate( ireqs_wsfc(nsends_wsfc,2) )

  do jrecv = 1, nrecvs_wsfc
     allocate( recv_wsfc(jrecv)%npts(1) )
     call MPI_IRecv(recv_wsfc(jrecv)%npts, 1, MPI_INTEGER, recv_wsfc(jrecv)%iremote, &
                    itag410, MPI_COMM_WORLD, ireqr_wsfc(jrecv,icurr_wsfc), ierr)
  enddo

  ! Post VSFC receives

  allocate( ireqr_vsfc(nrecvs_vsfc,2) )
  allocate( ireqs_vsfc(nsends_vsfc,2) )

  do jrecv = 1, nrecvs_vsfc
     allocate( recv_vsfc(jrecv)%npts(1) )
     call MPI_IRecv(recv_vsfc(jrecv)%npts, 1, MPI_INTEGER, recv_vsfc(jrecv)%iremote, &
                    itag510, MPI_COMM_WORLD, ireqr_vsfc(jrecv,icurr_vsfc), ierr)
  enddo

  ! Determine number of bytes to send per IWSFC column

  nbytes_per_iwsfc = 3       * nbytes_int  &
                   + 14      * nbytes_real &
                   + 2 * nzg * nbytes_real &
                   + 3 * nzs * nbytes_real

  ! Loop over all WSFC sends

  do jsend = 1, nsends_wsfc

     ! Send buffer sizes to receive ranks

     call MPI_Isend(jtab_wsfc_mpi(jsend)%jend, 1, MPI_INTEGER, &
                    send_wsfc(jsend)%iremote, itag410, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,icurr_wsfc), ierr)

     ! Determine size of send_wsfc buffer

     send_wsfc(jsend)%nbytes = nbytes_per_iwsfc * jtab_wsfc_mpi(jsend)%jend

     ! Allocate buffer

     allocate(send_wsfc(jsend)%buff(send_wsfc(jsend)%nbytes))

  enddo

  ! Determine number of bytes to send per IVSFC column

  nbytes_per_ivsfc = 8       * nbytes_real &
                   + 2 * nzg * nbytes_real

  ! Loop over all VSFC sends

  do jsend = 1, nsends_vsfc

     ! Send buffer sizes to receive ranks

     call MPI_Isend(jtab_vsfc_mpi(jsend)%jend, 1, MPI_INTEGER, &
                    send_vsfc(jsend)%iremote, itag510, MPI_COMM_WORLD, &
                    ireqs_vsfc(jsend,icurr_vsfc), ierr)

     ! Determine size of send_vsfc buffer

     send_vsfc(jsend)%nbytes = nbytes_per_ivsfc * jtab_vsfc_mpi(jsend)%jend

     ! Allocate buffer

     allocate(send_vsfc(jsend)%buff(send_vsfc(jsend)%nbytes))

  enddo

  ! Loop over all WSFC receives

  do jtmp = 1, nrecvs_wsfc

     ! Get completed recv_wsfc buffer sizes

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

  ! Loop over all VSFC receives

  do jtmp = 1, nrecvs_vsfc

     ! Get completed recv_vsfc buffer sizes

     call MPI_Waitany(nrecvs_vsfc, ireqr_vsfc(:,icurr_vsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     recv_vsfc(jrecv)%nbytes = nbytes_per_ivsfc * recv_vsfc(jrecv)%npts(1)

     ! Allocate recv_vsfc buffers for completed transfer

     allocate( recv_vsfc(jrecv)%buff( recv_vsfc(jrecv)%nbytes ) )
     allocate( recv_vsfc(jrecv)%ipts( recv_vsfc(jrecv)%npts(1) ) )

     ! Post next receive to get the list of ivsfc points we are getting

     call MPI_Irecv(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_vsfc(jrecv)%iremote, itag511, MPI_COMM_WORLD, &
                    ireqr_vsfc(jrecv,inext_vsfc), ierr)
  enddo

  ! Communicate the list of iwsfc points we are sending to adjacent nodes

  do jtmp = 1, nsends_wsfc

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_wsfc, ireqs_wsfc(:,icurr_wsfc), jsend, MPI_STATUS_IGNORE, ierr)

     allocate(iglobe( jtab_wsfc_mpi(jsend)%jend ) )

     do j = 1, jtab_wsfc_mpi(jsend)%jend
        iglobe(j) = itab_wsfc( jtab_wsfc_mpi(jsend)%iwsfc(j) )%iwglobe
     enddo

     ipos = 0

     call MPI_Pack(iglobe, jtab_wsfc_mpi(jsend)%jend, MPI_INTEGER, &
          send_wsfc(jsend)%buff, send_wsfc(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_wsfc(jsend)%buff, ipos, MPI_PACKED, &
                    send_wsfc(jsend)%iremote, itag411, MPI_COMM_WORLD, &
                    ireqs_wsfc(jsend,inext_wsfc), ierr)

     deallocate(iglobe)
  enddo

  ! Communicate the list of ivsfc points we are sending to adjacent nodes

  do jtmp = 1, nsends_vsfc

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_vsfc, ireqs_vsfc(:,icurr_vsfc), jsend, MPI_STATUS_IGNORE, ierr)

     allocate(iglobe( jtab_vsfc_mpi(jsend)%jend ) )

     do j = 1, jtab_vsfc_mpi(jsend)%jend
        iglobe(j) = itab_vsfc( jtab_vsfc_mpi(jsend)%ivsfc(j) )%ivglobe
     enddo

     ipos = 0

     call MPI_Pack(iglobe, jtab_vsfc_mpi(jsend)%jend, MPI_INTEGER, &
          send_vsfc(jsend)%buff, send_vsfc(jsend)%nbytes, ipos, MPI_COMM_WORLD, ierr)

     call MPI_Isend(send_vsfc(jsend)%buff, ipos, MPI_PACKED, &
                    send_vsfc(jsend)%iremote, itag511, MPI_COMM_WORLD, &
                    ireqs_vsfc(jsend,inext_vsfc), ierr)

     deallocate(iglobe)
  enddo

  ! Unpack and store the list of iwsfc points we are getting from adjacent nodes

  icurr_wsfc = mod(icurr_wsfc,2) + 1
  inext_wsfc = mod(inext_wsfc,2) + 1

  do jtmp = 1, nrecvs_wsfc
     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iglobe( recv_wsfc(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, ipos, &
                     iglobe, recv_wsfc(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_wsfc(jrecv)%npts(1)
        recv_wsfc(jrecv)%ipts(j) = itabg_wsfc( iglobe(j) )%iwsfc_myrank
     enddo

     deallocate(iglobe)

     call MPI_Irecv(recv_wsfc(jrecv)%buff, recv_wsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_wsfc(jrecv)%iremote, itagwsfc, MPI_COMM_WORLD,         &
                    ireqr_wsfc(jrecv,inext_wsfc), ierr                          )
  enddo

  ! Unpack and store the list of ivsfc points we are getting from adjacent nodes

  icurr_vsfc = mod(icurr_vsfc,2) + 1
  inext_vsfc = mod(inext_vsfc,2) + 1

  do jtmp = 1, nrecvs_vsfc
     call MPI_Waitany(nrecvs_vsfc, ireqr_vsfc(:,icurr_vsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ipos = 0

     allocate( iglobe( recv_vsfc(jrecv)%npts(1) ) )

     call MPI_Unpack(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, ipos, &
                     iglobe, recv_vsfc(jrecv)%npts(1), MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

     do j = 1, recv_vsfc(jrecv)%npts(1)
        recv_vsfc(jrecv)%ipts(j) = itabg_vsfc( iglobe(j) )%ivsfc_myrank
     enddo

     deallocate(iglobe)

     call MPI_Irecv(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_vsfc(jrecv)%iremote, itagvsfc, MPI_COMM_WORLD,         &
                    ireqr_vsfc(jrecv,inext_vsfc), ierr                          )
  enddo

#endif

end subroutine olam_alloc_mpi_sfcg

!===============================================================================

subroutine mpi_send_wsfc(set, soil_watfrac, div2d_ex)

  ! Subroutine to perform a parallel MPI send of a "WSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,    only: sfcg, itab_wsfc, jtab_wsfc_mpi, mwsfc
  use mem_land,    only: nzg, land, mland, omland
  use mem_lake,    only: lake, omlake
  use mem_sea,     only: sea, msea, omsea
  use leaf_coms,   only: nzs
  use sea_coms,    only: nzi
  use mem_para,    only: nsends_wsfc, send_wsfc, itagwsfc, ireqs_wsfc, &
                         icurr_wsfc, inext_wsfc

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: soil_watfrac(nzg,mland)
  real, optional, intent(inout) :: div2d_ex(msea)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: j, nb
  integer :: iwsfc, iland, ilake, isea

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_wsfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_wsfc, ireqs_wsfc(:,icurr_wsfc), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,j,iwsfc,iland,ilake,isea,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0

     ! Loop over number of columns for this jsend

     do j = 1,jtab_wsfc_mpi(jsend)%jend
        iwsfc = jtab_wsfc_mpi(jsend)%iwsfc(j)
        nb = send_wsfc(jsend)%nbytes

        ! Pack the messages into send buffers

        if (set == 'head_swm_grad') then

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(land%head   (1,iland),nzg,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(soil_watfrac(1,iland),nzg,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%gxps_vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gxps_vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gxps_vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%gyps_vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_progw') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%vmxet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmyet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmzet     (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmxet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmyet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vmzet_area(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%wdepth    (isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%head1   (iwsfc),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_div2d_ex') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(div2d_ex(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           endif
  
        elseif (set == 'swm_diagvel') then  ! [W call from surface_driver after swm_diagvel]

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then

              isea = iwsfc - omsea

              call MPI_Pack(sea%vxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%vze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           endif
           
        elseif (set == 'sfc_driv_end') then

           call MPI_Pack(sfcg%cantemp(iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%canrrv (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%wthv   (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rough  (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%nlev_seaice   (isea),  1,MPI_INTEGER, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_cantemp   (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_cantemp   (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_canrrv    (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_canrrv    (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_wthv      (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_wthv      (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%sea_rough     (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%ice_rough     (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaicec       (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seatc         (isea),  1,MPI_REAL,    &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%seaice_tempk(1,isea),nzi,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              call MPI_Pack(lake%lake_energy(ilake),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(lake%depth      (ilake),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Pack(land%sfcwater_mass  (1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_energy(1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%sfcwater_depth (1,iland),nzs, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%nlev_sfcwater    (iland),1,MPI_INTEGER, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_fracarea     (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_albedo       (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_rough        (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%veg_temp         (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_water     (1,iland),nzg, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%soil_energy    (1,iland),nzg, MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd             (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%ppfd_diffuse     (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(land%cosz             (iland),1,   MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'avgatm') then

           call MPI_Pack(sfcg%vels    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%prss    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%rhos    (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtemp (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airtheta(iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(sfcg%airrrv  (iwsfc),1,MPI_REAL, &
                         send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Pack(sea%windxe(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%windye(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sea%windze(isea),1,MPI_REAL, &
                            send_wsfc(jsend)%buff,nb,ipos,MPI_COMM_WORLD,ierr)

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

!===============================================================================

subroutine mpi_send_vsfc(watflux, energyflux, vc_ex)

  ! Subroutine to perform a parallel MPI send of a "VSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,    only: sfcg, itab_vsfc, jtab_vsfc_mpi, mvsfc
  use mem_land,    only: nzg
  use mem_para,    only: nsends_vsfc, send_vsfc, itagvsfc, ireqs_vsfc, &
                         icurr_vsfc, inext_vsfc
  use oname_coms,  only: nl

  implicit none

  real, optional, intent(inout) :: watflux   (nzg,mvsfc)
  real, optional, intent(inout) :: energyflux(nzg,mvsfc)
  real, optional, intent(inout) :: vc_ex         (mvsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: j, nb
  integer :: ivsfc, iw1, iw2

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_vsfc

     ! Make sure the previous sends are finished
     call MPI_Waitany(nsends_vsfc, ireqs_vsfc(:,icurr_vsfc), jsend, MPI_STATUS_IGNORE, ierr)

     !$omp task private(ipos,j,ivsfc,ierr) &
     !$omp      firstprivate(jsend) default(shared)

     ipos = 0

     ! Loop over number of VSFC points for this jsend

     do j = 1,jtab_vsfc_mpi(jsend)%jend
        ivsfc = jtab_vsfc_mpi(jsend)%ivsfc(j)
        nb = send_vsfc(jsend)%nbytes

        iw1 = itab_vsfc(ivsfc)%iwn(1)
        iw2 = itab_vsfc(ivsfc)%iwn(2)

        if (present(watflux)) then  ! follows ivsfc horiz flux loop in surface_driver

           call MPI_Pack(watflux   (1,ivsfc),nzg,MPI_REAL,send_vsfc(jsend)%buff, &
                                                 nb,ipos,MPI_COMM_WORLD,ierr)
           call MPI_Pack(energyflux(1,ivsfc),nzg,MPI_REAL,send_vsfc(jsend)%buff, &
                                                 nb,ipos,MPI_COMM_WORLD,ierr)

           ! Updates from swm_hflux; limited to swm_active points

           if ((sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) .and. &
                nl%igw_spinup /= 1) then

              call MPI_Pack(sfcg%hflux_wat(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_enr(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vxe(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vye(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%hflux_vze(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vmp      (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vmc      (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vc       (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                                    nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        elseif (present(vc_ex)) then ! follows vc_ex computation in swm_progv

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Pack(vc_ex(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                           nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        else  ! follows swm_progv call in surface_driver

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Pack(sfcg%vmc(ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                              nb,ipos,MPI_COMM_WORLD,ierr)
              call MPI_Pack(sfcg%vc (ivsfc),1,MPI_REAL,send_vsfc(jsend)%buff, &
                                              nb,ipos,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_vsfc(jsend)%buff, ipos, MPI_PACKED,            &
                    send_vsfc(jsend)%iremote, itagvsfc, MPI_COMM_WORLD, &
                    ireqs_vsfc(jsend,inext_vsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_vsfc = mod(icurr_vsfc,2) + 1
  inext_vsfc = mod(inext_vsfc,2) + 1

#endif

end subroutine mpi_send_vsfc

!=============================================================================

subroutine mpi_recv_wsfc(set, soil_watfrac, div2d_ex)

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
  use mem_sea,     only: sea, msea, omsea
  use mem_para,    only: nrecvs_wsfc, itagwsfc, recv_wsfc, ireqr_wsfc, icurr_wsfc, inext_wsfc

  implicit none

  character(*), optional, intent(in) :: set
  real, optional, intent(inout) :: soil_watfrac(nzg,mwsfc)
  real, optional, intent(inout) :: div2d_ex(msea)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: iwsfc, iland, ilake, isea
! integer :: iwglobe

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_wsfc

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_wsfc, ireqr_wsfc(:,icurr_wsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(j,ipos,iwsfc,iland,ilake,isea,ierr) &
     !$omp      firstprivate(jrecv) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

     do j = 1, recv_wsfc(jrecv)%npts(1)
        iwsfc = recv_wsfc(jrecv)%ipts(j)

        if (set == 'head_swm_grad') then

           if (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%head   (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              soil_watfrac(1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gxps_vze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%gyps_vze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_progw') then

           if (sfcg%leaf_class(iwsfc) == 0 .and. sfcg%swm_active(iwsfc)) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmxet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmyet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmzet     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmxet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmyet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vmzet_area(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%wdepth    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sfcg%head1   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_div2d_ex') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              div2d_ex     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'swm_diagvel') then

           if (sfcg%swm_active(iwsfc) .and. sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vxe    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vye    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%vze    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif
           
        elseif (set == 'sfc_driv_end') then

           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%cantemp  (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%canrrv   (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%wthv      (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                           sfcg%rough     (iwsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%nlev_seaice (isea),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_cantemp (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_cantemp (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_canrrv  (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_canrrv  (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_wthv    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_wthv    (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%sea_rough   (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%ice_rough   (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaicec     (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seatc       (isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%seaice_tempk(1,isea),nzi,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) == 1) then
              ilake = iwsfc - omlake

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              lake%lake_energy(ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              lake%depth      (ilake),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           elseif (sfcg%leaf_class(iwsfc) >= 2) then
              iland = iwsfc - omland

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_mass  (1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_energy(1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%sfcwater_depth (1,iland),nzs,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%nlev_sfcwater    (iland),  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_fracarea     (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_albedo       (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_rough        (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%veg_temp         (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_water     (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%soil_energy    (1,iland),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd             (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%ppfd_diffuse     (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              land%cosz             (iland),  1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (set == 'avgatm') then

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

           if (sfcg%leaf_class(iwsfc) == 0) then
              isea = iwsfc - omsea

              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windxe(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windye(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_wsfc(jrecv)%buff,recv_wsfc(jrecv)%nbytes,ipos, &
                              sea%windze(isea),1,MPI_REAL,MPI_COMM_WORLD,ierr)

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

!=============================================================================

subroutine mpi_recv_vsfc(watflux, energyflux, vc_ex)

  ! Subroutine to perform a parallel MPI receive of a "VSFC group"
  ! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_sfcg,    only: sfcg, mvsfc, itab_vsfc
  use mem_land,    only: nzg, land, omland
  use mem_para,    only: nrecvs_vsfc, itagvsfc, recv_vsfc, ireqr_vsfc, icurr_vsfc, inext_vsfc
  use oname_coms,  only: nl

  implicit none

  real, optional, intent(inout) :: watflux   (nzg,mvsfc)
  real, optional, intent(inout) :: energyflux(nzg,mvsfc)
  real, optional, intent(inout) :: vc_ex         (mvsfc)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: ivsfc, iw1, iw2

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_vsfc

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_vsfc, ireqr_vsfc(:,icurr_vsfc), jrecv, MPI_STATUS_IGNORE, ierr)

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task private(j,ipos,ivsfc,ierr) &
     !$omp      firstprivate(jrecv) default(shared)

     ipos = 0

     ! Loop over number of columns for this jtmp/jrecv

     do j = 1, recv_vsfc(jrecv)%npts(1)
        ivsfc = recv_vsfc(jrecv)%ipts(j)

        iw1 = itab_vsfc(ivsfc)%iwn(1)
        iw2 = itab_vsfc(ivsfc)%iwn(2)

        if (present(watflux)) then  ! follows ivsfc horiz flux loop in surface_driver

           call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                           watflux   (1,ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)
           call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                           energyflux(1,ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

           ! Updates from swm_hflux; limited to swm_active points

           if ((sfcg%swm_active(iw1) .or. sfcg%swm_active(iw2)) .and. &
                nl%igw_spinup /= 1) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_wat(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_enr(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vxe(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vye(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%hflux_vze(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmp      (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmc      (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vc       (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        elseif (present(vc_ex)) then ! follows vc_ex computation in swm_progv

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              vc_ex(ivsfc),nzg,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        else  ! follows swm_progv call in surface_driver

           if ((sfcg%swm_active(iw1) .and. sfcg%leaf_class(iw1) == 0) .and. &
               (sfcg%swm_active(iw2) .and. sfcg%leaf_class(iw2) == 0)) then

              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vmc(ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)
              call MPI_Unpack(recv_vsfc(jrecv)%buff,recv_vsfc(jrecv)%nbytes,ipos, &
                              sfcg%vc (ivsfc),1,MPI_REAL,MPI_COMM_WORLD,ierr)

           endif

        endif

     enddo

     call MPI_Irecv(recv_vsfc(jrecv)%buff, recv_vsfc(jrecv)%nbytes, MPI_PACKED, &
                    recv_vsfc(jrecv)%iremote, itagvsfc, MPI_COMM_WORLD,         &
                    ireqr_vsfc(jrecv,inext_vsfc), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_vsfc

End Module olam_mpi_sfcg
