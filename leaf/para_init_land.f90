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
subroutine para_init_land(landflag, landflux_temp)

use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank,  &
!                     send_ul, recv_ul,  &
                      send_wl, recv_wl,  &
                      send_wlf, recv_wlf,  &
!                     nsends_ul, nrecvs_ul,  &
                      nsends_wl, nrecvs_wl,  &
                      nsends_wlf, nrecvs_wlf

use leaf_coms,  only: mml, mul, mwl, nml, nul, nwl, nzg

use mem_leaf,   only: itab_ml, itab_ul, itab_wl,       &
                      itabg_ml, itabg_ul, itabg_wl,    &
                      land_vars, land, alloc_land_grid

use mem_mksfc,  only: itab_mls_vars, itab_uls_vars, itab_wls_vars
                        
use mem_sflux,  only: landflux, nlandflux, jlandflux, flux_vars

implicit none

logical,         intent(in) :: landflag(nwl)
type(flux_vars), intent(in) :: landflux_temp(nlandflux)

integer :: jml
integer :: iw,iml,iul,iwl
integer :: ilf
integer :: im1l,im2l,iw1l,iw2l

integer :: iml_myrank = 1 ! Counter for ML points to be included on this rank
integer :: iul_myrank = 1 ! Counter for UL points to be included on this rank
integer :: iwl_myrank = 1 ! Counter for WL points to be included on this rank

! Automatic arrays

logical :: myrankflag_ml(nml) ! Flag for ML points existing on this rank
logical :: myrankflag_ul(nul) ! Flag for UL points existing on this rank
logical :: myrankflag_wl(nwl) ! Flag for WL points existing on this rank

! Temporary datatypes

type(itab_mls_vars), allocatable :: ltab_ml(:)
type(itab_uls_vars), allocatable :: ltab_ul(:)
type(itab_wls_vars), allocatable :: ltab_wl(:)
type(land_vars)                  :: land_t

! Move data to temporary data structures, nullifying the old datatype

call move_alloc (itab_ml, ltab_ml)
call move_alloc (itab_ul, ltab_ul)
call move_alloc (itab_wl, ltab_wl)

call move_alloc (land%leaf_class, land_t%leaf_class)
call move_alloc (land%ntext_soil, land_t%ntext_soil)

call move_alloc (land%area , land_t%area )
call move_alloc (land%xew  , land_t%xew  )
call move_alloc (land%yew  , land_t%yew  )
call move_alloc (land%zew  , land_t%zew  )
call move_alloc (land%glatw, land_t%glatw)
call move_alloc (land%glonw, land_t%glonw)
call move_alloc (land%wnx  , land_t%wnx  )
call move_alloc (land%wny  , land_t%wny  )
call move_alloc (land%wnz  , land_t%wnz  )

call move_alloc (land%xem  , land_t%xem  )
call move_alloc (land%yem  , land_t%yem  )
call move_alloc (land%zem  , land_t%zem  )
call move_alloc (land%zm   , land_t%zm   )
call move_alloc (land%glatm, land_t%glatm)
call move_alloc (land%glonm, land_t%glonm)

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_wl (mrls)) ; nsends_wl  = 0
allocate (nsends_wlf(mrls)) ; nsends_wlf = 0
allocate (nrecvs_wl (mrls)) ; nrecvs_wl  = 0
allocate (nrecvs_wlf(mrls)) ; nrecvs_wlf = 0

! Initialize myrank flag arrays to .false.

myrankflag_ml(1:nml) = .false.
myrankflag_ul(1:nul) = .false.
myrankflag_wl(1:nwl) = .false.

!------------------------------------------------------------------------------
! The following section will be required for landcell-to-landcell transport of 
! water via rivers and aquifers, which is not yet implemented in LEAF.  The MPI
! send/recv tables set up by this section will be distinct from those set up
! for landcell-atmosphere fluxes in the next section, and will exchange all
! prognostic LEAF quantities.)
!------------------------------------------------------------------------------

!t Loop over all UL points, and for each whose assigned irank is equal to myrank, 
!t flag its WL neighbors for inclusion on this rank.

!t do iul = 2,nul
!t    if (itabg_ul(iul)%irank == myrank) then

!t       iw1l = ltab_ul(iul)%iw(1)
!t       iw2l = ltab_ul(iul)%iw(2)

!t       myrankflag_wl(iw1l) = .true.
!t       myrankflag_wl(iw2l) = .true.

!t    endif
!t enddo

! Loop over all WL points, and for each whose landflag value is .true., flag 
! the WL point for inclusion on this rank.  Then, for any WL point flagged, 
! flag its UL and ML neighbors for inclusion on this rank.

do iwl = 2,nwl
   if (landflag(iwl)) then
      myrankflag_wl(iwl) = .true.
   endif

   if (myrankflag_wl(iwl)) then

      do jml = 1,ltab_wl(iwl)%npoly
         iml = ltab_wl(iwl)%im(jml)
         iul = ltab_wl(iwl)%iu(jml)

         myrankflag_ml(iml) = .true.
         myrankflag_ul(iul) = .true.
      enddo

   endif
enddo

! Loop over all ML, UL, and WL points and count the ones that have been flagged
! for inclusion on this rank.

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = iml_myrank + 1
   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = iul_myrank + 1
   endif
enddo

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1
   endif
enddo

! Set mml, mul, mwl values for this rank

mml = iml_myrank
mul = iul_myrank
mwl = iwl_myrank

! Re-allocate itab data structures and main grid coordinate arrays

call alloc_land_grid(mml,mul,mwl,nzg)

! Reset point counts to 1

iml_myrank = 1
iul_myrank = 1
iwl_myrank = 1

! Store new myrank ML, UL, WL indices in itabg data structures

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = iml_myrank + 1

      itabg_ml(iml)%iml_myrank = iml_myrank      
   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = iul_myrank + 1

      itabg_ul(iul)%iul_myrank = iul_myrank
   endif
enddo

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1

      itabg_wl(iwl)%iwl_myrank = iwl_myrank      
   endif
enddo

! Memory copy to main tables

do iml = 2,nml
   if (myrankflag_ml(iml)) then
      iml_myrank = itabg_ml(iml)%iml_myrank

      itab_ml(iml_myrank)%imglobe = iml

      land%xem  (iml_myrank) = land_t%xem  (iml)
      land%yem  (iml_myrank) = land_t%yem  (iml)
      land%zem  (iml_myrank) = land_t%zem  (iml)
      land%zm   (iml_myrank) = land_t%zm   (iml)
      land%glatm(iml_myrank) = land_t%glatm(iml)
      land%glonm(iml_myrank) = land_t%glonm(iml)

   endif
enddo

do iul = 2,nul
   if (myrankflag_ul(iul)) then
      iul_myrank = itabg_ul(iul)%iul_myrank

      itab_ul(iul_myrank)%irank = itabg_ul(iul)%irank
      itab_ul(iul_myrank)%iuglobe = iul

      im1l = ltab_ul(iul)%im(1)
      im2l = ltab_ul(iul)%im(2)
      iw1l = ltab_ul(iul)%iw(1)
      iw2l = ltab_ul(iul)%iw(2)

      if (myrankflag_ml(im1l)) itab_ul(iul_myrank)%im(1) = itabg_ml(im1l)%iml_myrank
      if (myrankflag_ml(im2l)) itab_ul(iul_myrank)%im(2) = itabg_ml(im2l)%iml_myrank
      if (myrankflag_wl(iw1l)) itab_ul(iul_myrank)%iw(1) = itabg_wl(iw1l)%iwl_myrank
      if (myrankflag_wl(iw2l)) itab_ul(iul_myrank)%iw(2) = itabg_wl(iw2l)%iwl_myrank
   endif
enddo

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = itabg_wl(iwl)%iwl_myrank

      itab_wl(iwl_myrank)%irank   = itabg_wl(iwl)%irank
      itab_wl(iwl_myrank)%iwglobe = iwl
      itab_wl(iwl_myrank)%npoly   = ltab_wl(iwl)%npoly

      land%leaf_class  (iwl_myrank) = land_t%leaf_class  (iwl)
      land%ntext_soil(:,iwl_myrank) = land_t%ntext_soil(:,iwl)

      land%area (iwl_myrank) = land_t%area (iwl)
      land%xew  (iwl_myrank) = land_t%xew  (iwl)
      land%yew  (iwl_myrank) = land_t%yew  (iwl)
      land%zew  (iwl_myrank) = land_t%zew  (iwl)
      land%glatw(iwl_myrank) = land_t%glatw(iwl)
      land%glonw(iwl_myrank) = land_t%glonw(iwl)
      land%wnx  (iwl_myrank) = land_t%wnx  (iwl)
      land%wny  (iwl_myrank) = land_t%wny  (iwl)
      land%wnz  (iwl_myrank) = land_t%wnz  (iwl)

      do jml = 1,ltab_wl(iwl)%npoly
         iml = ltab_wl(iwl)%im(jml)
         iul = ltab_wl(iwl)%iu(jml)
         
         if (myrankflag_ml(iml)) itab_wl(iwl_myrank)%im(jml) =  &
                                 itabg_ml(iml)%iml_myrank

         if (myrankflag_ul(iul)) itab_wl(iwl_myrank)%iu(jml) =  &
                                 itabg_ul(iul)%iul_myrank
      enddo

   endif
enddo

!------------------------------------------------------------------------------
! The following section will be required for landcell-to-landcell transport of 
! water via rivers and aquifers, which is not yet implemented in LEAF.  The MPI
! send/recv tables set up by this section will be distinct from those set up
! for landcell-atmosphere fluxes in the next section, and will exchange all
! prognostic LEAF quantities.)
!------------------------------------------------------------------------------

! Loop over all UL points and get the 2 neighbor WL indices of each

!t do iul = 2,nul
!t    iw1l = ltab_ul(iul)%iw(1)
!t    iw2l = ltab_ul(iul)%iw(2)

!t    if (itabg_ul(iul)%irank == myrank) then

! UL point is on myrank.  If either WL neighbor is on remote rank, myrank must
! receive that WL cell.

!t       if (iw1l > 1 .and. itabg_wl(iw1l)%irank /= myrank) then
!t          call recv_table_wl(itabg_wl(iw1l)%irank)
!t       endif

!t       if (iw2l > 1 .and. itabg_wl(iw2l)%irank /= myrank) then
!t          call recv_table_wl(itabg_wl(iw2l)%irank)
!t       endif

!t    elseif (itabg_ul(iul)%irank /= myrank) then

! UL point is on remote rank.  If either WL neighbor is on myrank, myrank must
! send that WL cell.

!t       if (iw1l > 1 .and. itabg_wl(iw1l)%irank == myrank) then
!t          call send_table_wl(iw1l,itabg_ul(iul)%irank)
!t       endif

!t       if (iw2l > 1 .and. itabg_wl(iw2l)%irank == myrank) then
!t          call send_table_wl(iw2l,itabg_ul(iul)%irank)
!t       endif

!t    endif

!t enddo

! Loop over all landflux cells
! (This section is required for computing landcell-atmosphere fluxes and 
! involves only surface and canopy properties.)

do ilf = 2,nlandflux

   iw  = landflux_temp(ilf)%iw   ! full-domain index
   iwl = landflux_temp(ilf)%iwls ! full-domain index

   if (itabg_w(iw)%irank /= myrank .and. itabg_wl(iwl)%irank == myrank) then  

! ATM cell is on remote rank while LAND cell is on myrank, so landflux is computed
! on remote rank and must be received by myrank, and LAND cell fields must be
! sent to remote rank.

      call send_table_wl(iwl,itabg_w(iw)%irank)
      call recv_table_wlf(itabg_w(iw)%irank)

   endif

   if (itabg_w(iw)%irank == myrank .and. itabg_wl(iwl)%irank /= myrank) then  

! ATM cell is on myrank while LAND cell is on remote rank, so landflux is computed
! on myrank and must be sent to remote rank, and LAND cell fields must be
! received from remote rank.

      call send_table_wlf(ilf,itabg_wl(iwl)%irank)
      call recv_table_wl(itabg_wl(iwl)%irank)

   endif

enddo

! Deallocate temporary data structures and arrays

deallocate (ltab_ml, ltab_ul, ltab_wl)

deallocate (land_t%leaf_class, land_t%ntext_soil)

deallocate (land_t%area, land_t%xew, land_t%yew, land_t%zew)
deallocate (land_t%glatw, land_t%glonw)
deallocate (land_t%wnx, land_t%wny, land_t%wnz)

deallocate (land_t%xem, land_t%yem, land_t%zem, land_t%zm)
deallocate (land_t%glatm, land_t%glonm)

end subroutine para_init_land

!===============================================================================

subroutine recv_table_wl(iremote)

use mem_para,  only: nrecvs_wl, recv_wl
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wl(1) .and. recv_wl(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wl(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wl(1).

if (jrecv > nrecvs_wl(1)) nrecvs_wl(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wl(jrecv)%iremote = iremote

return
end subroutine recv_table_wl

!===============================================================================

subroutine recv_table_wlf(iremote)

use mem_para,  only: nrecvs_wlf, recv_wlf, myrank
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wlf(1) .and. recv_wlf(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wlf(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wlf(1).

if (jrecv > nrecvs_wlf(1)) nrecvs_wlf(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wlf(jrecv)%iremote = iremote

return
end subroutine recv_table_wlf

!===============================================================================

subroutine send_table_wl(iwl,iremote)

use mem_leaf,  only: itab_wl, itabg_wl
use mem_para,  only: nsends_wl, send_wl
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwl
integer, intent(in) :: iremote

integer :: jsend
integer :: iwl_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_wl(1) .and. send_wl(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_wl(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wl(1).

if (jsend > nsends_wl(1)) nsends_wl(1) = jsend

! If nsends_wl(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wl, nsends_wl(1) exceeds maxremote'
   write(io6,*) 'nsends_wl(1) = ',nsends_wl(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wl '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iwl_myrank = itabg_wl(iwl)%iwl_myrank

itab_wl(iwl_myrank)%send(jsend) = .true.

send_wl(jsend)%iremote = iremote

return
end subroutine send_table_wl

!===============================================================================

subroutine send_table_wlf(ilf,iremote)

use mem_sflux, only: landflux, landfluxg
use mem_para,  only: nsends_wlf, send_wlf
use max_dims,   only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: ilf
integer, intent(in) :: iremote

integer :: jsend
integer :: ilf_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_wlf(1) .and. send_wlf(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_wlf(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wlf(1).

if (jsend > nsends_wlf(1)) nsends_wlf(1) = jsend

! If nsends_wlf(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wlf, nsends_wlf(1) exceeds maxremote'
   write(io6,*) 'nsends_wlf(1) = ',nsends_wlf(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wlf '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

ilf_myrank = landfluxg(ilf)%ilf_myrank

landflux(ilf_myrank)%sendf(jsend) = .true.

send_wlf(jsend)%iremote = iremote

return
end subroutine send_table_wlf

