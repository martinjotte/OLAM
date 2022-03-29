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
subroutine para_init_sfcg()

  use misc_coms,  only: io6
  use mem_ijtabs, only: itabg_w, itab_w, mrls
  use mem_sfcg,   only: mmsfc, mvsfc, mwsfc, nmsfc, nvsfc, nwsfc, &
                        itab_msfc, itabg_msfc, itab_msfc_pd, &
                        itab_vsfc, itabg_vsfc, itab_vsfc_pd, &
                        itab_wsfc, itabg_wsfc, itab_wsfc_pd, &
                        alloc_sfcgrid1, sfcg
  use mem_para,   only: myrank, nsends_wsfc, nrecvs_wsfc, &
                        nsends_vsfc, nrecvs_vsfc
  use mem_land,   only: mland, omland, onland, itab_land
  use mem_lake,   only: mlake, omlake, onlake, itab_lake
  use mem_sea,    only: msea, omsea, onsea, itab_sea
  use pom2k1d,    only: alloc_pomgrid

  implicit none

  integer :: iw, npoly, j, jj, iwn, imn, ivn, ivnn, imsfc, ivsfc, iwsfc, iwatm

  integer :: imsfc_myrank = 1 ! Counter for SFC grid M points to be included on this rank
  integer :: ivsfc_myrank = 1 ! Counter for SFC grid V points to be included on this rank
  integer :: iwsfc_myrank = 1 ! Counter for SFC grid W points to be included on this rank

  integer :: iland_myrank = 1 ! Counter for LAND points to be included on this rank
  integer :: ilake_myrank = 1 ! Counter for LAKE points to be included on this rank
  integer :: isea_myrank  = 1 ! Counter for SEA points to be included on this rank

  ! Automatic arrays

  logical :: myrankflag_msfc(nmsfc) ! Flag for SFC grid M points existing on this rank
  logical :: myrankflag_vsfc(nvsfc) ! Flag for SFC grid V points existing on this rank
  logical :: myrankflag_wsfc(nwsfc) ! Flag for SFC grid W points existing on this rank

  ! Allocate send & recv counter arrays and initialize to zero

  nsends_vsfc = 0
  nsends_wsfc = 0

  nrecvs_vsfc = 0
  nrecvs_wsfc = 0

  myrankflag_msfc(:) = .false.
  myrankflag_vsfc(:) = .false.
  myrankflag_wsfc(:) = .false.

  allocate(itabg_msfc(nmsfc))
  allocate(itabg_vsfc(nvsfc))
  allocate(itabg_wsfc(nwsfc))

  ! Loop over all global SFC grid W points

  do iwsfc = 2,nwsfc

     ! Set the rank of the WSFC point to that of the ATM cell to which
     ! it is coupled with the largest coupling area

     iw = itab_wsfc_pd(iwsfc)%iwatm(1) ! global index
     itabg_wsfc(iwsfc)%irank = itabg_w(iw)%irank

     ! If the rank of the WSFC cell is the same as myrank, then
     !    (1) flag the WSFC point for inclusion on this rank
     !    (2) flag all IWN and IMN neighbors of WSFC cell for inclusion on this rank.
     !    (3) flag all IVN neighbors of flagged IMN neighbors for inclusion on this rank

     if (itabg_wsfc(iwsfc)%irank == myrank) then
        myrankflag_wsfc(iwsfc) = .true.

        npoly = itab_wsfc_pd(iwsfc)%npoly

        do j = 1,npoly
           iwn = itab_wsfc_pd(iwsfc)%iwn(j)
           myrankflag_wsfc(iwn) = .true.

           imn = itab_wsfc_pd(iwsfc)%imn(j)
           myrankflag_msfc(imn) = .true.

           do jj = 1,3
              ivnn = itab_msfc_pd(imn)%ivn(jj)
              myrankflag_vsfc(ivnn) = .true.
           enddo

        enddo

     endif

     ! If the SFC grid W point is coupled to an ATM cell whose rank is myrank,
     ! flag the SFC grid point for inclusion on this rank

     do j = 1,itab_wsfc_pd(iwsfc)%nwatm
        iw = itab_wsfc_pd(iwsfc)%iwatm(j) ! global index

        if (itabg_w(iw)%irank == myrank) then
           myrankflag_wsfc(iwsfc) = .true.
        endif
     enddo

  enddo

  ! Loop over all global SFC grid W points.  For those that have been flagged for
  ! inclusion on this rank, increment counters of WSFC, LAND, LAKE, and SEA cells.

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

  ! Count SFC grid M and V points that have been flagged for inclusion,
  ! and set their rank to that of their first WSFC neighbor

  do imsfc = 2,nmsfc
     if (myrankflag_msfc(imsfc)) imsfc_myrank = imsfc_myrank + 1
     iwn = itab_msfc_pd(imsfc)%iwn(1)
     itabg_msfc(imsfc)%irank = itabg_wsfc(iwn)%irank
  enddo

  do ivsfc = 2,nvsfc
     if (myrankflag_vsfc(ivsfc)) ivsfc_myrank = ivsfc_myrank + 1
     iwn = itab_vsfc_pd(ivsfc)%iwn(1)
     itabg_vsfc(ivsfc)%irank = itabg_wsfc(iwn)%irank
  enddo

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

  do imsfc = 2,nmsfc
     if (myrankflag_msfc(imsfc)) then
        imsfc_myrank = imsfc_myrank + 1
        itabg_msfc(imsfc)%imsfc_myrank = imsfc_myrank
     endif
  enddo

  do ivsfc = 2,nvsfc
     if (myrankflag_vsfc(ivsfc)) then
        ivsfc_myrank = ivsfc_myrank + 1
        itabg_vsfc(ivsfc)%ivsfc_myrank = ivsfc_myrank
     endif
  enddo

  do iwsfc = 2,nwsfc
     if (myrankflag_wsfc(iwsfc)) then
        iwsfc_myrank = iwsfc_myrank + 1
        itabg_wsfc(iwsfc)%iwsfc_myrank = iwsfc_myrank
     endif
  enddo

  ! Copy and convert global itab_wsfc_pd index information to local rank itab_wsfc

  do iwsfc = 2,nwsfc ! global index

     if (myrankflag_wsfc(iwsfc)) then

        iwsfc_myrank = itabg_wsfc(iwsfc)%iwsfc_myrank

        itab_wsfc(iwsfc_myrank)%irank    = itabg_wsfc(iwsfc)%irank
        itab_wsfc(iwsfc_myrank)%iwglobe  = iwsfc
        itab_wsfc(iwsfc_myrank)%npoly    = itab_wsfc_pd(iwsfc)%npoly
        itab_wsfc(iwsfc_myrank)%nwatm    = itab_wsfc_pd(iwsfc)%nwatm

        do j = 1,7
           imn = itab_wsfc_pd(iwsfc)%imn(j) ! global index
           ivn = itab_wsfc_pd(iwsfc)%ivn(j) ! global index
           iwn = itab_wsfc_pd(iwsfc)%iwn(j) ! global index

           if (myrankflag_msfc(imn)) then
              itab_wsfc(iwsfc_myrank)%imn(j) = itabg_msfc(imn)%imsfc_myrank  ! local index
           else
              itab_wsfc(iwsfc_myrank)%imn(j) = 1
           endif

           if (myrankflag_vsfc(ivn)) then
              itab_wsfc(iwsfc_myrank)%ivn(j) = itabg_vsfc(ivn)%ivsfc_myrank  ! local index
           else
              itab_wsfc(iwsfc_myrank)%ivn(j) = 1
           endif

           if (myrankflag_wsfc(iwn)) then
              itab_wsfc(iwsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank  ! local index
           else
              itab_wsfc(iwsfc_myrank)%iwn(j) = 1
           endif
        enddo

        do j = 1,itab_wsfc_pd(iwsfc)%nwatm 
           iwatm = itab_wsfc_pd(iwsfc)%iwatm(j)  ! global index
           iw    = itabg_w(iwatm)%iw_myrank      ! local index (or -1 if not on this rank)

           itab_wsfc(iwsfc_myrank)%iwatm(j) = iw ! local index (or -1 if not on this rank)
        enddo

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

  ! Fill itab_msfc table

  do imsfc = 2,nmsfc
     if (myrankflag_msfc(imsfc)) then
        imsfc_myrank = itabg_msfc(imsfc)%imsfc_myrank
        itab_msfc(imsfc_myrank)%irank   = itabg_msfc(imsfc)%irank
        itab_msfc(imsfc_myrank)%imglobe = imsfc

        ! The following is apparently needed only for contslab.f90

        do j = 1,3
           ivn = itab_msfc_pd(imsfc)%ivn(j) ! global index
           iwn = itab_msfc_pd(imsfc)%iwn(j) ! global index
           itab_msfc(imsfc_myrank)%ivn(j) = itabg_vsfc(ivn)%ivsfc_myrank  ! local index
           itab_msfc(imsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank  ! local index
        enddo
     endif
  enddo

  ! Fill itab_vsfc table

  do ivsfc = 2,nvsfc
     if (myrankflag_vsfc(ivsfc)) then
        ivsfc_myrank = itabg_vsfc(ivsfc)%ivsfc_myrank
        itab_vsfc(ivsfc_myrank)%irank   = itabg_vsfc(ivsfc)%irank
        itab_vsfc(ivsfc_myrank)%ivglobe = ivsfc

        do j = 1,2
           imn = itab_vsfc_pd(ivsfc)%imn(j) ! global index
           iwn = itab_vsfc_pd(ivsfc)%iwn(j) ! global index
           itab_vsfc(ivsfc_myrank)%imn(j) = itabg_msfc(imn)%imsfc_myrank  ! local index
           itab_vsfc(ivsfc_myrank)%iwn(j) = itabg_wsfc(iwn)%iwsfc_myrank  ! local index
        enddo
     endif
  enddo

  ! Read surface grid information (that does not contain neighbor or coupling
  ! indices) for all points in local parallel subdomain for this rank (or for
  ! all points in domain if run is sequential)

  call sfcgfile_read()

  ! Fill MPI send/recv tables

  ! Loop over all SFC grid W points in global domain

  do iwsfc = 2,nwsfc

     ! If this SFC grid W point is primary on a remote rank and is also in the
     ! memory of myrank, add that remote irank to MPI receive table.

     if (itabg_wsfc(iwsfc)%irank /= myrank .and. myrankflag_wsfc(iwsfc)) then
        call recv_table_wsfc(itabg_wsfc(iwsfc)%irank)
     endif

     ! If this SFC grid W point is primary on a remote rank and is adjacent to
     ! a point that is primary on myrank, add the myrank-primary point global
     ! index and that remote rank to MPI send table.

     if (itabg_wsfc(iwsfc)%irank /= myrank) then
        do j = 1,itab_wsfc_pd(iwsfc)%npoly
           iwn = itab_wsfc_pd(iwsfc)%iwn(j) ! Global index

           if (itabg_wsfc(iwn)%irank == myrank) then
              call send_table_wsfc(iwn,itabg_wsfc(iwsfc)%irank)
           endif
        enddo
     endif

     ! If this SFC grid W point is primary on myrank and is coupled to an ATM
     ! grid column on a remote rank, add this SFC grid W point global index
     ! and that remote rank to MPI send table.

     if (itabg_wsfc(iwsfc)%irank == myrank) then
        do j = 1,itab_wsfc_pd(iwsfc)%nwatm
           iwatm = itab_wsfc_pd(iwsfc)%iwatm(j) ! Global index

           if (itabg_w(iwatm)%irank /= myrank) then
              call send_table_wsfc(iwsfc,itabg_w(iwatm)%irank)
           endif
        enddo
     endif

  enddo

  ! Loop over all VSFC points in global domain

  do ivsfc = 2,nvsfc

     ! If this VSFC grid point is primary on a remote rank and is also in the
     ! memory of myrank, add that remote irank to MPI receive table.

     if (itabg_vsfc(ivsfc)%irank /= myrank .and. myrankflag_vsfc(ivsfc)) then
        call recv_table_vsfc(itabg_vsfc(ivsfc)%irank)
     endif

     ! If this VSFC grid point is primary on this rank and an adjacent WSFC point
     ! (of which there are 4) is primary on a remote rank, add this VSFC point
     ! index and that remote rank to send table.

     if (itabg_vsfc(ivsfc)%irank == myrank) then
        do j = 1,2
           imn = itab_vsfc_pd(ivsfc)%iwn(j)

           do jj = 1,3
              iwn = itab_msfc_pd(imn)%iwn(jj)

              if (itabg_wsfc(iwn)%irank /= myrank) &
                 call send_table_vsfc(ivsfc,itabg_wsfc(iwn)%irank)
           enddo
        enddo
     endif

  enddo

  ! Deallocate temporary arrays

  deallocate(itab_wsfc_pd, itab_vsfc_pd, itab_msfc_pd)
  
end subroutine para_init_sfcg 

!===============================================================================

subroutine recv_table_wsfc(iremote)

use mem_para,  only: nrecvs_wsfc, recv_wsfc
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote_wsfc is already in table of ranks to receive from

do jrecv = 1, nrecvs_wsfc
   if (recv_wsfc(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_wsfc, jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wsfc.

if (jrecv > nrecvs_wsfc) nrecvs_wsfc = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wsfc(jrecv)%iremote = iremote

end subroutine recv_table_wsfc

!===============================================================================

subroutine recv_table_vsfc(iremote)

use mem_para,  only: nrecvs_vsfc, recv_vsfc
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote_vsfc is already in table of ranks to receive from

do jrecv=1,nrecvs_vsfc
   if (recv_vsfc(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_vsfc, jrecv represents a rank not yet entered in the
! table, so increase nrecvs_vsfc.

if (jrecv > nrecvs_vsfc) nrecvs_vsfc = jrecv

! Enter remote rank in recv-remote-rank table.

recv_vsfc(jrecv)%iremote = iremote

end subroutine recv_table_vsfc

!===============================================================================

subroutine send_table_wsfc(iwsfc,iremote)

use mem_sfcg, only: itab_wsfc, itabg_wsfc
use mem_para,   only: nsends_wsfc, send_wsfc
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iwsfc
integer, intent(in) :: iremote

integer :: jsend
integer :: iwsfc_myrank

! Check whether iremote_wsfc is already in table of ranks to send to

do jsend=1, nsends_wsfc
   if (send_wsfc(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_wsfc, jsend represents a rank not yet entered in the
! table, so increase nsends_wsfc.

if (jsend > nsends_wsfc) nsends_wsfc = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iwsfc_myrank = itabg_wsfc(iwsfc)%iwsfc_myrank

itab_wsfc(iwsfc_myrank)%send(jsend) = .true.
send_wsfc(jsend)%iremote = iremote

end subroutine send_table_wsfc

!===============================================================================

subroutine send_table_vsfc(ivsfc,iremote)

use mem_sfcg,  only: itab_vsfc, itabg_vsfc
use mem_para,  only: nsends_vsfc, send_vsfc
use misc_coms, only: io6

implicit none

integer, intent(in) :: ivsfc
integer, intent(in) :: iremote

integer :: jsend
integer :: ivsfc_myrank

! Check whether iremote_vsfc is already in table of ranks to send to

do jsend=1,nsends_vsfc
   if (send_vsfc(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_vsfc, jsend represents a rank not yet entered in the
! table, so increase nsends_vsfc.

if (jsend > nsends_vsfc) nsends_vsfc = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

ivsfc_myrank = itabg_vsfc(ivsfc)%ivsfc_myrank

itab_vsfc(ivsfc_myrank)%send(jsend) = .true.
send_vsfc(jsend)%iremote = iremote

end subroutine send_table_vsfc

