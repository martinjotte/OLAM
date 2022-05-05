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

subroutine init_sfcgrid()

  use misc_coms, only: iparallel
  implicit none

  if (iparallel == 1) then

     call para_init_sfc()

  else

     call serial_init_sfc()

  endif

end subroutine init_sfcgrid

!===============================================================================

subroutine para_init_sfc()

  use misc_coms,    only: iparallel
  use mem_ijtabs,   only: itabg_w
  use mem_sfcg,     only: mmsfc, mvsfc, mwsfc, nmsfc, nvsfc, nwsfc, &
                          itab_msfc, itabg_msfc, itab_msfc_pd, &
                          itab_vsfc, itabg_vsfc, itab_vsfc_pd, &
                          itab_wsfc, itabg_wsfc, itab_wsfc_pd, &
                          alloc_sfcgrid1, sfcg
  use mem_para,     only: myrank, mgroupsize, nbytes_int, nbytes_real
  use olam_mpi_sfc
  use max_dims,     only: maxremote
  use mem_land,     only: mland, omland, onland, itab_land, nzg
  use mem_lake,     only: mlake, omlake, onlake, itab_lake
  use mem_sea,      only: msea, omsea, onsea, itab_sea
  use pom2k1d,      only: alloc_pomgrid
  use leaf_coms,    only: nzs

  implicit none

  integer :: iw, npoly, j, iwn, imn, ivn, imsfc, ivsfc, iwsfc, iwatm
  integer :: jsend, jrecv, jend, ips

  integer :: imsfc_myrank ! Counter for SFC grid M points to be included on this rank
  integer :: ivsfc_myrank ! Counter for SFC grid V points to be included on this rank
  integer :: iwsfc_myrank ! Counter for SFC grid W points to be included on this rank

  integer :: iland_myrank ! Counter for LAND points to be included on this rank
  integer :: ilake_myrank ! Counter for LAKE points to be included on this rank
  integer :: isea_myrank  ! Counter for SEA points to be included on this rank

  integer :: nbytes_per_iwsfc
  integer :: nbytes_per_ivsfc
  integer :: nbytes_per_imsfc

  ! Automatic arrays

  logical :: myrankflag_msfc(nmsfc) ! Flag for SFC grid M points existing on this rank
  logical :: myrankflag_vsfc(nvsfc) ! Flag for SFC grid V points existing on this rank
  logical :: myrankflag_wsfc(nwsfc) ! Flag for SFC grid W points existing on this rank

  integer :: send_vsfc_iremotes(maxremote) ! holds ranks we will send to
  integer :: send_wsfc_iremotes(maxremote) ! holds ranks we will send to
  integer :: send_msfc_iremotes(maxremote) ! holds ranks we will send to

  integer :: recv_vsfc_iremotes(maxremote) ! holds ranks we will recv from
  integer :: recv_wsfc_iremotes(maxremote) ! holds ranks we will recv from
  integer :: recv_msfc_iremotes(maxremote) ! holds ranks we will recv from

  ! Allocatable arrays

  logical, allocatable :: ivsfc_sends(:,:) ! temporary send table
  logical, allocatable :: iwsfc_sends(:,:) ! temporary send table
  logical, allocatable :: imsfc_sends(:,:) ! temporary send table

  integer, allocatable :: ivsfc_recvs(:) ! temporary recv table
  integer, allocatable :: iwsfc_recvs(:) ! temporary recv table
  integer, allocatable :: imsfc_recvs(:) ! temporary recv table

  ! Initialize send & recv counters to zero

  nsends_vsfc = 0
  nsends_msfc = 0
  nsends_wsfc = 0

  nrecvs_vsfc = 0
  nrecvs_wsfc = 0
  nrecvs_msfc = 0

  ! Initialize myrank flag arrays to .false.

  myrankflag_msfc(:) = .false.
  myrankflag_vsfc(:) = .false.
  myrankflag_wsfc(:) = .false.

  ! Loop over all VSFC points, and for each whose assigned irank is equal to
  ! myrank, flag all MSFC, VSFC, and WSFC points in its computational stencil
  ! for inclusion on this rank

  do ivsfc = 2,nvsfc
     if (itabg_vsfc(ivsfc)%irank == myrank) then

        myrankflag_vsfc( ivsfc ) = .true.
        myrankflag_msfc( itab_vsfc_pd(ivsfc)%imn(1:2) ) = .true.
        myrankflag_wsfc( itab_vsfc_pd(ivsfc)%iwn(1:2) ) = .true.

     endif
  enddo

  ! Loop over all WSFC points, and for each whose assigned irank is equal to
  ! myrank, flag all MSFC, VSFC, and WSFC points in its computational stencil
  ! for inclusion on this rank

  do iwsfc = 2,nwsfc
     if (itabg_wsfc(iwsfc)%irank == myrank) then

        npoly = itab_wsfc_pd(iwsfc)%npoly

        myrankflag_wsfc( iwsfc ) = .true.
        myrankflag_msfc( itab_wsfc_pd(iwsfc)%imn(1:npoly) ) = .true.
        myrankflag_vsfc( itab_wsfc_pd(iwsfc)%ivn(1:npoly) ) = .true.
        myrankflag_wsfc( itab_wsfc_pd(iwsfc)%iwn(1:npoly) ) = .true.

     endif

     ! If the SFC grid W point is coupled to an ATM cell whose rank is myrank,
     ! also flag the SFC grid point for inclusion on this rank.

     do j = 1,itab_wsfc_pd(iwsfc)%nwatm
        iw = itab_wsfc_pd(iwsfc)%iwatm(j) ! global index
        if (itabg_w(iw)%irank == myrank) myrankflag_wsfc(iwsfc) = .true.
     enddo

  enddo

  ! Loop over all MSFC points, and for each whose assigned irank is equal to
  ! myrank, flag all MSFC, VSFC, and WSFC points in its computational stencil
  ! for inclusion on this rank

  do imsfc = 2,nmsfc
     if (itabg_msfc(imsfc)%irank == myrank) then

        myrankflag_msfc( imsfc ) = .true.
        myrankflag_msfc( itab_msfc_pd(imsfc)%imn(1:3) ) = .true.
        myrankflag_vsfc( itab_msfc_pd(imsfc)%ivn(1:3) ) = .true.
        myrankflag_wsfc( itab_msfc_pd(imsfc)%iwn(1:3) ) = .true.

     endif
  enddo

  ! Ignore the "dummy" 1st point in case it was activated

  myrankflag_msfc(1) = .false.
  myrankflag_vsfc(1) = .false.
  myrankflag_wsfc(1) = .false.

  ! Loop over all MSFC, VSFC, and WSFC points and count the ones that
  ! have been flagged for inclusion on this rank.

  imsfc_myrank = 1
  ivsfc_myrank = 1
  iwsfc_myrank = 1

  iland_myrank = 1
  ilake_myrank = 1
  isea_myrank  = 1

  do imsfc = 2,nmsfc
     if (myrankflag_msfc(imsfc)) imsfc_myrank = imsfc_myrank + 1
  enddo

  do ivsfc = 2,nvsfc
     if (myrankflag_vsfc(ivsfc)) ivsfc_myrank = ivsfc_myrank + 1
  enddo

  do iwsfc = 2,nwsfc
     if (myrankflag_wsfc(iwsfc)) then
        iwsfc_myrank = iwsfc_myrank + 1

        if (itab_wsfc_pd(iwsfc)%leaf_class >= 2) then
           iland_myrank = iland_myrank + 1
        elseif (itab_wsfc_pd(iwsfc)%leaf_class == 1) then
           ilake_myrank = ilake_myrank + 1
        else
           isea_myrank = isea_myrank + 1
        endif

     endif
  enddo

  ! Set mmsfc, mvsfc, mwsfc values for this rank

  mmsfc = imsfc_myrank
  mvsfc = ivsfc_myrank
  mwsfc = iwsfc_myrank

  mland = iland_myrank
  mlake = ilake_myrank
  msea  = isea_myrank

  omland = 0
  omlake = mland - 1
  omsea  = mland + mlake - 2

  ! Allocate SFC grid itab and lsg data structures and arrays

  call alloc_sfcgrid1(mmsfc, mvsfc, mwsfc)

  allocate (itab_land(mland))
  allocate (itab_lake(mlake))
  allocate (itab_sea (msea ))

  call alloc_pomgrid(msea)

  ! Reset myrank point counters to 1

  imsfc_myrank = 1
  ivsfc_myrank = 1
  iwsfc_myrank = 1

  iland_myrank = 1
  ilake_myrank = 1
  isea_myrank  = 1

  ! Store new local rank MSFC, VSFC, and WSFC indices in itabg data structure
  ! (Land, lake, sea W cells automatically in correct order on this rank)

  do imsfc = 2, nmsfc
     if (myrankflag_msfc(imsfc)) then
        imsfc_myrank = imsfc_myrank + 1

        itabg_msfc(imsfc)%imsfc_myrank = imsfc_myrank
        itab_msfc(imsfc_myrank)%imglobe = imsfc
        itab_msfc(imsfc_myrank)%irank = itabg_msfc(imsfc)%irank
     endif
  enddo

  do ivsfc = 2, nvsfc
     if (myrankflag_vsfc(ivsfc)) then
        ivsfc_myrank = ivsfc_myrank + 1

        itabg_vsfc(ivsfc)%ivsfc_myrank = ivsfc_myrank
        itab_vsfc(ivsfc_myrank)%ivglobe = ivsfc
        itab_vsfc(ivsfc_myrank)%irank = itabg_vsfc(ivsfc)%irank
     endif
  enddo

  do iwsfc = 2, nwsfc
     if (myrankflag_wsfc(iwsfc)) then
        iwsfc_myrank = iwsfc_myrank + 1

        itabg_wsfc(iwsfc)%iwsfc_myrank = iwsfc_myrank
        itab_wsfc(iwsfc_myrank)%iwglobe = iwsfc
        itab_wsfc(iwsfc_myrank)%irank = itabg_wsfc(iwsfc)%irank

        ! Fill itab_land, itab_lake, and itab_sea global indices

        if (itab_wsfc_pd(iwsfc)%leaf_class >= 2) then
           iland_myrank = iwsfc_myrank - omland
           itab_land(iland_myrank)%iwglobe = iwsfc - onland
        elseif (itab_wsfc_pd(iwsfc)%leaf_class == 1) then
           ilake_myrank = iwsfc_myrank - omlake
           itab_lake(ilake_myrank)%iwglobe = iwsfc - onlake
        else
           isea_myrank = iwsfc_myrank - omsea
           itab_sea(isea_myrank)%iwglobe = iwsfc - onsea
        endif
     endif
  enddo

  ! Read surface grid information (that does not contain neighbor or coupling
  ! indices) for all points in local parallel subdomain for this rank (or for
  ! all points in domain if run is sequential)

  call sfcgfile_read()

  ! Copy and convert global itab_msfc_pd index information to local rank itab_msfc

  do imsfc = 2, nmsfc ! global index

     if (myrankflag_msfc(imsfc)) then

        imsfc_myrank = itabg_msfc(imsfc)%imsfc_myrank

        ! Set indices of IMSFC neighbors that are present on this rank

        do j = 1, 3
           ivn = itab_msfc_pd(imsfc)%ivn(j) ! global index
           iwn = itab_msfc_pd(imsfc)%iwn(j) ! global index
           imn = itab_msfc_pd(imsfc)%imn(j) ! global index

           if (myrankflag_vsfc(ivn)) &
                itab_msfc(imsfc_myrank)%ivn(j) = itabg_vsfc(ivn)%ivsfc_myrank

           if (myrankflag_wsfc(iwn)) &
                itab_msfc(imsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank

           if (myrankflag_msfc(imn)) &
                itab_msfc(imsfc_myrank)%imn(j) = itabg_msfc(imn)%imsfc_myrank
        enddo

     endif
  enddo

  ! Copy and convert global itab_vsfc_pd index information to local rank itab_vsfc

  do ivsfc = 2, nvsfc ! global index

     if (myrankflag_vsfc(ivsfc)) then

        ivsfc_myrank = itabg_vsfc(ivsfc)%ivsfc_myrank

        ! Set indices of IVSFC neighbors that are present on this rank

        do j = 1,2
           imn = itab_vsfc_pd(ivsfc)%imn(j) ! global index
           iwn = itab_vsfc_pd(ivsfc)%iwn(j) ! global index

           if (myrankflag_msfc(imn)) &
                itab_vsfc(ivsfc_myrank)%imn(j) = itabg_msfc(imn)%imsfc_myrank

           if (myrankflag_wsfc(iwn)) &
                itab_vsfc(ivsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank
        enddo

     endif
  enddo

  ! Copy and convert global itab_wsfc_pd index information to local rank itab_wsfc

  do iwsfc = 2, nwsfc ! global index

     if (myrankflag_wsfc(iwsfc)) then

        iwsfc_myrank = itabg_wsfc(iwsfc)%iwsfc_myrank

        itab_wsfc(iwsfc_myrank)%npoly = itab_wsfc_pd(iwsfc)%npoly
        itab_wsfc(iwsfc_myrank)%nwatm = itab_wsfc_pd(iwsfc)%nwatm
        sfcg%leaf_class(iwsfc_myrank) = itab_wsfc_pd(iwsfc)%leaf_class

        do j = 1,itab_wsfc_pd(iwsfc)%npoly
           imn = itab_wsfc_pd(iwsfc)%imn(j) ! global index
           ivn = itab_wsfc_pd(iwsfc)%ivn(j) ! global index
           iwn = itab_wsfc_pd(iwsfc)%iwn(j) ! global index

           if (myrankflag_msfc(imn)) &
                itab_wsfc(iwsfc_myrank)%imn(j) = itabg_msfc(imn)%imsfc_myrank

           if (myrankflag_vsfc(ivn)) &
                itab_wsfc(iwsfc_myrank)%ivn(j) = itabg_vsfc(ivn)%ivsfc_myrank

           if (myrankflag_wsfc(iwn)) &
                itab_wsfc(iwsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank
        enddo

        do j = 1,itab_wsfc_pd(iwsfc)%nwatm
           iwatm = itab_wsfc_pd(iwsfc)%iwatm(j)  ! global index
           iw    = itabg_w(iwatm)%iw_myrank      ! local index (or -1 if not on this rank)

           itab_wsfc(iwsfc_myrank)%iwatm(j) = iw ! local index (or -1 if not on this rank)
        enddo

     endif
  enddo

  ! Allocate/initialize temporary send/recv tables

  ips = min(maxremote, mgroupsize-1)

  allocate(ivsfc_sends(mvsfc,ips)) ; ivsfc_sends = .false.
  allocate(imsfc_sends(mmsfc,ips)) ; imsfc_sends = .false.
  allocate(iwsfc_sends(mwsfc,ips)) ; iwsfc_sends = .false.

  allocate(ivsfc_recvs(mvsfc)) ; ivsfc_recvs = -1
  allocate(imsfc_recvs(mmsfc)) ; imsfc_recvs = -1
  allocate(iwsfc_recvs(mwsfc)) ; iwsfc_recvs = -1

  ! Fill temporary MPI send/receive tables

  do ivsfc = 2, nvsfc ! global index
     if (itabg_vsfc(ivsfc)%irank /= myrank) then

        ! IVSFC point is primary on a remote rank. Add to send table any V
        ! or W point that is primary on myrank and is in the stencil of IVSFC.

        do j = 1,2
           imn = itab_vsfc_pd(ivsfc)%imn(j)
           if (itabg_msfc(imn)%irank == myrank) &
                call send_table_msfc(imn,itabg_vsfc(ivsfc)%irank)

           iwn = itab_vsfc_pd(ivsfc)%iwn(j)
           if (itabg_wsfc(iwn)%irank == myrank) &
                call send_table_wsfc(iwn,itabg_vsfc(ivsfc)%irank)
        enddo

        ! If IVSFC is in our memory and primary on a remote rank,
        ! add to receive table

        if (myrankflag_vsfc(ivsfc)) &
             call recv_table_vsfc(ivsfc,itabg_vsfc(ivsfc)%irank)

     endif
  enddo

  do iwsfc = 2, nwsfc ! global index
     if (itabg_wsfc(iwsfc)%irank /= myrank) then

        ! IWSFC point is primary on a remote rank. Add to send table any V, W,
        ! or M point that is primary on myrank and is in the stencil of IWSFC.

        do j = 1, itab_wsfc_pd(iwsfc)%npoly

           ivn = itab_wsfc_pd(iwsfc)%ivn(j)
           if (itabg_vsfc(ivn)%irank == myrank) &
                call send_table_vsfc(ivn,itabg_wsfc(iwsfc)%irank)

           iwn = itab_wsfc_pd(iwsfc)%iwn(j)
           if (itabg_wsfc(iwn)%irank == myrank) &
                call send_table_wsfc(iwn,itabg_wsfc(iwsfc)%irank)

           imn = itab_wsfc_pd(iwsfc)%imn(j)
           if (itabg_msfc(imn)%irank == myrank) &
                call send_table_msfc(imn,itabg_wsfc(iwsfc)%irank)

        enddo

        ! If IWSFC is in our memory and primary on a remote rank,
        ! add to receive table

        if (myrankflag_wsfc(iwsfc)) &
             call recv_table_wsfc(iwsfc,itabg_wsfc(iwsfc)%irank)

     else

        ! If this SFC grid W point is primary on myrank and is coupled to an ATM
        ! grid column on a remote rank, add this SFC grid W point global index
        ! and that remote rank to MPI send table.

        do j = 1,itab_wsfc_pd(iwsfc)%nwatm
           iwatm = itab_wsfc_pd(iwsfc)%iwatm(j) ! Global index
           if (itabg_w(iwatm)%irank /= myrank) &
                call send_table_wsfc(iwsfc,itabg_w(iwatm)%irank)
        enddo

     endif
  enddo

  do imsfc = 2,nmsfc ! global index
     if (itabg_msfc(imsfc)%irank /= myrank) then

        ! IMSFC point is primary on a remote rank. Add to send table any V, W,
        ! or M point that is primary on myrank and is in the stencil of IMSFC.

        do j = 1, 3
           ivn = itab_msfc_pd(imsfc)%ivn(j)
           if (itabg_vsfc(ivn)%irank == myrank) &
                call send_table_vsfc(ivn,itabg_msfc(imsfc)%irank)

           iwn = itab_msfc_pd(imsfc)%iwn(j)
           if (itabg_wsfc(iwn)%irank == myrank) &
                call send_table_wsfc(iwn,itabg_msfc(imsfc)%irank)

           imn = itab_msfc_pd(imsfc)%imn(j)
           if (itabg_msfc(imn)%irank == myrank) &
                call send_table_msfc(imn,itabg_msfc(imsfc)%irank)
        enddo

        ! If IMSFC is in our memory and primary on a remote rank,
        ! add to receive table

        if (myrankflag_msfc(imsfc)) &
             call recv_table_msfc(imsfc,itabg_msfc(imsfc)%irank)

     endif
  enddo

  ! Allocate main send and receive tables

  if (iparallel == 1) then

  ! Determine number of bytes to send per IWSFC column

  nbytes_per_iwsfc =  3       * nbytes_int  &
                   + 20       * nbytes_real &
                   +  2 * nzg * nbytes_real &
                   +  3 * nzs * nbytes_real

  ! Copy data from temp arrays to wsfc send table

  allocate( send_wsfc(nsends_wsfc) )

  do jsend = 1, nsends_wsfc
     jend = count( iwsfc_sends(2:mwsfc,jsend) )

     send_wsfc(jsend)%jend    = jend
     send_wsfc(jsend)%iremote = send_wsfc_iremotes(jsend)
     send_wsfc(jsend)%nbytes  = jend * nbytes_per_iwsfc

     allocate( send_wsfc(jsend)%ipts( send_wsfc(jsend)%jend   ) )
     allocate( send_wsfc(jsend)%buff( send_wsfc(jsend)%nbytes ) )

     j = 0
     do iwsfc = 2, mwsfc
        if (iwsfc_sends(iwsfc,jsend)) then
           j = j + 1
           send_wsfc(jsend)%ipts(j) = iwsfc
        endif
     enddo
  enddo

  ! Determine number of bytes to send per IVSFC column

  nbytes_per_ivsfc = 8       * nbytes_real &
                   + 2 * nzg * nbytes_real

  ! Copy data from temp arrays to vsfc send table

  allocate( send_vsfc(nsends_vsfc) )

  do jsend = 1, nsends_vsfc
     jend = count( ivsfc_sends(2:mvsfc,jsend) )

     send_vsfc(jsend)%jend    = jend
     send_vsfc(jsend)%iremote = send_vsfc_iremotes(jsend)
     send_vsfc(jsend)%nbytes  = jend * nbytes_per_ivsfc

     allocate( send_vsfc(jsend)%ipts( send_vsfc(jsend)%jend   ) )
     allocate( send_vsfc(jsend)%buff( send_vsfc(jsend)%nbytes ) )

     j = 0
     do ivsfc = 2, mvsfc
        if (ivsfc_sends(ivsfc,jsend)) then
           j = j + 1
           send_vsfc(jsend)%ipts(j) = ivsfc
        endif
     enddo
  enddo

  ! Determine number of bytes to send per IMSFC column

  nbytes_per_imsfc = 1 * nbytes_real

  ! Copy data from temp arrays to msfc send table

  allocate( send_msfc(nsends_msfc) )

  do jsend = 1, nsends_msfc
     jend = count( imsfc_sends(2:mmsfc,jsend) )

     send_msfc(jsend)%jend    = jend
     send_msfc(jsend)%iremote = send_msfc_iremotes(jsend)
     send_msfc(jsend)%nbytes  = jend * nbytes_per_imsfc

     allocate( send_msfc(jsend)%ipts( send_msfc(jsend)%jend   ) )
     allocate( send_msfc(jsend)%buff( send_msfc(jsend)%nbytes ) )

     j = 0
     do imsfc = 2, mmsfc
        if (imsfc_sends(imsfc,jsend)) then
           j = j + 1
           send_msfc(jsend)%ipts(j) = imsfc
        endif
     enddo
  enddo

  ! Copy data from temp arrays to wsfc recv table

  allocate( recv_wsfc(nrecvs_wsfc) )

  do jrecv = 1, nrecvs_wsfc
     jend = count( iwsfc_recvs(2:mwsfc) == jrecv )

     recv_wsfc(jrecv)%jend    = jend
     recv_wsfc(jrecv)%iremote = recv_wsfc_iremotes(jrecv)
     recv_wsfc(jrecv)%nbytes  = jend * nbytes_per_iwsfc

     allocate( recv_wsfc(jrecv)%ipts( recv_wsfc(jrecv)%jend   ) )
     allocate( recv_wsfc(jrecv)%buff( recv_wsfc(jrecv)%nbytes ) )

     j = 0
     do iwsfc = 2, mwsfc
        if (iwsfc_recvs(iwsfc) == jrecv) then
           j = j + 1
           recv_wsfc(jrecv)%ipts(j) = iwsfc
        endif
     enddo
  enddo

  ! Copy data from temp arrays to vsfc recv table

  allocate( recv_vsfc(nrecvs_vsfc) )

  do jrecv = 1, nrecvs_vsfc
     jend = count( ivsfc_recvs(2:mvsfc) == jrecv )

     recv_vsfc(jrecv)%jend    = jend
     recv_vsfc(jrecv)%iremote = recv_vsfc_iremotes(jrecv)
     recv_vsfc(jrecv)%nbytes  = jend * nbytes_per_ivsfc

     allocate( recv_vsfc(jrecv)%ipts( recv_vsfc(jrecv)%jend   ) )
     allocate( recv_vsfc(jrecv)%buff( recv_vsfc(jrecv)%nbytes ) )

     j = 0
     do ivsfc = 2, mvsfc
        if (ivsfc_recvs(ivsfc) == jrecv) then
           j = j + 1
           recv_vsfc(jrecv)%ipts(j) = ivsfc
        endif
     enddo
  enddo

  ! Copy data from temp arrays to msfc recv table

  allocate( recv_msfc(nrecvs_msfc) )

  do jrecv = 1, nrecvs_msfc
     jend = count( imsfc_recvs(2:mmsfc) == jrecv )

     recv_msfc(jrecv)%jend    = jend
     recv_msfc(jrecv)%iremote = recv_msfc_iremotes(jrecv)
     recv_msfc(jrecv)%nbytes  = jend * nbytes_per_imsfc

     allocate( recv_msfc(jrecv)%ipts( recv_msfc(jrecv)%jend   ) )
     allocate( recv_msfc(jrecv)%buff( recv_msfc(jrecv)%nbytes ) )

     j = 0
     do imsfc = 2, mmsfc
        if (imsfc_recvs(imsfc) == jrecv) then
           j = j + 1
           recv_msfc(jrecv)%ipts(j) = imsfc
        endif
     enddo
  enddo

  endif

  ! Deallocate temporary arrays

  deallocate(itab_wsfc_pd, itab_vsfc_pd, itab_msfc_pd)


Contains


  subroutine recv_table_wsfc(iwsfc,iremote)

    use mem_sfcg,     only: itabg_wsfc
    use olam_mpi_sfc, only: nrecvs_wsfc

    implicit none

    integer, intent(in)  :: iwsfc    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jrecv

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_wsfc
       if (recv_wsfc_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_wsfc, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_wsfc.

    if (jrecv > nrecvs_wsfc) then
       nrecvs_wsfc = jrecv
       recv_wsfc_iremotes(jrecv) = iremote
    endif

    ! Enter remote rank in recv-remote-rank table.

    iwsfc_recvs( itabg_wsfc(iwsfc)%iwsfc_myrank ) = jrecv

  end subroutine recv_table_wsfc


  subroutine recv_table_vsfc(ivsfc,iremote)

    use mem_sfcg,     only: itabg_vsfc
    use olam_mpi_sfc, only: nrecvs_vsfc

    implicit none

    integer, intent(in)  :: ivsfc    ! global index
    integer, intent(in)  :: iremote  ! remote rank
    integer              :: jrecv

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_vsfc
       if (recv_vsfc_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_vsfc, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_vsfc.

    if (jrecv > nrecvs_vsfc) then
       nrecvs_vsfc = jrecv
       recv_vsfc_iremotes(jrecv) = iremote
    endif

    ! Enter remote rank in recv-remote-rank table.

    ivsfc_recvs( itabg_vsfc(ivsfc)%ivsfc_myrank ) = jrecv

  end subroutine recv_table_vsfc


  subroutine recv_table_msfc(imsfc,iremote)

    use mem_sfcg,     only: itabg_msfc
    use olam_mpi_sfc, only: nrecvs_msfc

    implicit none

    integer, intent(in)  :: imsfc    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jrecv

    ! Check whether iremote is already in table of ranks to receive from

    do jrecv = 1, nrecvs_msfc
       if (recv_msfc_iremotes(jrecv) == iremote) exit
    enddo

    ! If jrecv exceeds nrecvs_msfc, jrecv represents a rank not yet entered in
    ! the table, so increase nrecvs_msfc.

    if (jrecv > nrecvs_msfc) then
       nrecvs_msfc = jrecv
       recv_msfc_iremotes(jrecv) = iremote
    endif

    ! Enter remote rank in recv-remote-rank table.

    imsfc_recvs( itabg_msfc(imsfc)%imsfc_myrank ) = jrecv

  end subroutine recv_table_msfc


  subroutine send_table_wsfc(iwsfc,iremote)

    use mem_sfcg,     only: itabg_wsfc
    use olam_mpi_sfc, only: nsends_wsfc

    implicit none

    integer,  intent(in) :: iwsfc    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_wsfc is already in table of ranks to send to

    do jsend = 1, nsends_wsfc
       if (send_wsfc_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_wsfc, jsend represents a rank not yet entered in
    ! the table, so increase nsends_wsfc.

    if (jsend > nsends_wsfc) then
       nsends_wsfc = jsend
       send_wsfc_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table.

    iwsfc_sends( itabg_wsfc(iwsfc)%iwsfc_myrank, jsend ) = .true.

  end subroutine send_table_wsfc


  subroutine send_table_vsfc(ivsfc,iremote)

    use mem_sfcg,     only: itabg_vsfc
    use olam_mpi_sfc, only: nsends_vsfc

    implicit none

    integer,  intent(in) :: ivsfc    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_vsfc is already in table of ranks to send to

    do jsend=1,nsends_vsfc
       if (send_vsfc_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_vsfc, jsend represents a rank not yet entered in
    ! the table, so increase nsends_vsfc.

    if (jsend > nsends_vsfc) then
       nsends_vsfc = jsend
       send_vsfc_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table.

    ivsfc_sends( itabg_vsfc(ivsfc)%ivsfc_myrank, jsend ) = .true.

  end subroutine send_table_vsfc


  subroutine send_table_msfc(imsfc,iremote)

    use mem_sfcg,     only: itabg_msfc
    use olam_mpi_sfc, only: nsends_msfc

    implicit none

    integer,  intent(in) :: imsfc    ! global index
    integer,  intent(in) :: iremote  ! remote rank
    integer              :: jsend

    ! Check whether iremote_msfc is already in table of ranks to send to

    do jsend=1,nsends_msfc
       if (send_msfc_iremotes(jsend) == iremote) exit
    enddo

    ! If jsend exceeds nsends_msfc, jsend represents a rank not yet entered in the
    ! the table, so increase nsends_msfc.

    if (jsend > nsends_msfc) then
       nsends_msfc = jsend
       send_msfc_iremotes(jsend) = iremote
    endif

    ! Enter point in send-point table.

    imsfc_sends( itabg_msfc(imsfc)%imsfc_myrank, jsend ) = .true.

  end subroutine send_table_msfc


end subroutine para_init_sfc

!===============================================================================

subroutine serial_init_sfc()

  use mem_sfcg,     only: mmsfc, mvsfc, mwsfc, nmsfc, nvsfc, nwsfc, &
                          itab_msfc, itabg_msfc, itab_msfc_pd, &
                          itab_vsfc, itabg_vsfc, itab_vsfc_pd, &
                          itab_wsfc, itabg_wsfc, itab_wsfc_pd, &
                          alloc_sfcgrid1, sfcg
  use mem_land,     only: mland, nland, itab_land
  use mem_lake,     only: mlake, nlake, itab_lake
  use mem_sea,      only: msea, nsea, itab_sea
  use pom2k1d,      only: alloc_pomgrid

  implicit none

  integer :: imsfc, ivsfc, iwsfc
  integer :: iland, ilake, isea

  ! Allocate SFC grid itab and lsg data structures and arrays

  call alloc_sfcgrid1(mmsfc, mvsfc, mwsfc)

  allocate (itab_land(mland))
  allocate (itab_lake(mlake))
  allocate (itab_sea (msea ))

  call alloc_pomgrid(msea)

  do imsfc = 1, nmsfc
     itabg_msfc(imsfc)%imsfc_myrank = imsfc
     itab_msfc (imsfc)%imglobe      = imsfc
     itab_msfc (imsfc)%irank        = 0
  enddo

  do ivsfc = 1, nvsfc
     itabg_vsfc(ivsfc)%ivsfc_myrank = ivsfc
     itab_vsfc (ivsfc)%ivglobe      = ivsfc
     itab_vsfc (ivsfc)%irank        = 0
  enddo

  do iwsfc = 1, nwsfc
     itabg_wsfc(iwsfc)%iwsfc_myrank = iwsfc
     itab_wsfc (iwsfc)%iwglobe      = iwsfc
     itab_wsfc (iwsfc)%irank        = 0
  enddo

  ! Fill itab_land, itab_lake, and itab_sea global indices

  do iland = 1, nland
     itab_land(iland)%iwglobe = iland
  enddo

  do ilake = 1, nlake
     itab_lake(ilake)%iwglobe = ilake
  enddo

  do isea = 1, nsea
     itab_sea(isea)%iwglobe = isea
  enddo

  ! Read surface grid information (that does not contain neighbor or coupling
  ! indices) for all points in local parallel subdomain for this rank (or for
  ! all points in domain if run is sequential)

  call sfcgfile_read()

  ! Copy global itab_msfc_pd index information to itab_msfc

  do imsfc = 2, mmsfc
     itab_msfc(imsfc)%ivn(1:3) = itab_msfc_pd(imsfc)%ivn(1:3)
     itab_msfc(imsfc)%iwn(1:3) = itab_msfc_pd(imsfc)%iwn(1:3)
     itab_msfc(imsfc)%imn(1:3) = itab_msfc_pd(imsfc)%imn(1:3)
  enddo

  ! Copy global itab_vsfc_pd index information to itab_vsfc

  do ivsfc = 2, mvsfc
     itab_vsfc(ivsfc)%imn(1:2) = itab_vsfc_pd(ivsfc)%imn(1:2)
     itab_vsfc(ivsfc)%iwn(1:2) = itab_vsfc_pd(ivsfc)%iwn(1:2)
  enddo

  ! Copy global itab_wsfc_pd index information to itab_wsfc

  do iwsfc = 2, mwsfc
     itab_wsfc(iwsfc)%npoly = itab_wsfc_pd(iwsfc)%npoly
     itab_wsfc(iwsfc)%nwatm = itab_wsfc_pd(iwsfc)%nwatm

     sfcg%leaf_class(iwsfc) = itab_wsfc_pd(iwsfc)%leaf_class

     itab_wsfc(iwsfc)%imn(1:7) = itab_wsfc_pd(iwsfc)%imn(1:7)
     itab_wsfc(iwsfc)%ivn(1:7) = itab_wsfc_pd(iwsfc)%ivn(1:7)
     itab_wsfc(iwsfc)%iwn(1:7) = itab_wsfc_pd(iwsfc)%iwn(1:7)

     itab_wsfc(iwsfc)%iwatm(1:8) = itab_wsfc_pd(iwsfc)%iwatm(1:8)
  enddo

end subroutine serial_init_sfc
