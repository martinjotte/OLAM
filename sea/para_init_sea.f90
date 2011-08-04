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
subroutine para_init_sea(seaflag, seaflux_temp)

use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank,  &
!                     send_us, recv_us,  &
                      send_ws, recv_ws,  &
                      send_wsf, recv_wsf,  &
!                     nsends_us, nrecvs_us,  &
                      nsends_ws, nrecvs_ws,  &
                      nsends_wsf, nrecvs_wsf

use sea_coms,   only: mms, mus, mws, nms, nus, nws

use mem_sea,    only: itab_ms, itab_us, itab_ws,             &
                      itabg_ms, itabg_us, itabg_ws,          &
                      sea_vars, sea, alloc_sea_grid

use mem_mksfc,  only: itab_mls_vars, itab_uls_vars, itab_wls_vars

use mem_sflux,  only: seaflux, nseaflux, jseaflux, flux_vars

implicit none

logical,         intent(in) :: seaflag(nws)
type(flux_vars), intent(in) :: seaflux_temp(nseaflux)

integer :: j,jms,jsend
integer :: iw,ims,ius,iws
integer :: isf
integer :: im1s,im2s, iw1s,iw2s

integer :: ims_myrank = 1 ! Counter for MS points to be included on this rank
integer :: ius_myrank = 1 ! Counter for US points to be included on this rank
integer :: iws_myrank = 1 ! Counter for WS points to be included on this rank

! Automatic arrays

logical :: myrankflag_ms(nms) ! Flag for MS points existing on this rank
logical :: myrankflag_us(nus) ! Flag for US points existing on this rank
logical :: myrankflag_ws(nws) ! Flag for WS points existing on this rank

! Temporary datatypes

type(itab_mls_vars), allocatable :: ltab_ms(:)
type(itab_uls_vars), allocatable :: ltab_us(:)
type(itab_wls_vars), allocatable :: ltab_ws(:)
type(sea_vars)                   :: sea_t

! Move data to temporary data structures, nullifying the old datatype

call move_alloc(itab_ms, ltab_ms)
call move_alloc(itab_us, ltab_us)
call move_alloc(itab_ws, ltab_ws)

call move_alloc (sea%leaf_class, sea_t%leaf_class)

call move_alloc (sea%area , sea_t%area )
call move_alloc (sea%xew  , sea_t%xew  )
call move_alloc (sea%yew  , sea_t%yew  )
call move_alloc (sea%zew  , sea_t%zew  )
call move_alloc (sea%glatw, sea_t%glatw)
call move_alloc (sea%glonw, sea_t%glonw)

call move_alloc (sea%xem  , sea_t%xem  )
call move_alloc (sea%yem  , sea_t%yem  )
call move_alloc (sea%zem  , sea_t%zem  )
call move_alloc (sea%zm   , sea_t%zm   )
call move_alloc (sea%glatm, sea_t%glatm)
call move_alloc (sea%glonm, sea_t%glonm)

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_ws (mrls)) ; nsends_ws  = 0
allocate (nsends_wsf(mrls)) ; nsends_wsf = 0
allocate (nrecvs_ws (mrls)) ; nrecvs_ws  = 0
allocate (nrecvs_wsf(mrls)) ; nrecvs_wsf = 0

! Initialize myrank flag arrays to .false.

myrankflag_ms(:) = .false.
myrankflag_us(:) = .false.
myrankflag_ws(:) = .false.

!------------------------------------------------------------------------------
! The following section (which would correspond to landcell-to-landcell exchange
! of information between parallel processors) is likely not required for sea 
! cells.  The prognostic ocean model will perform its own parallel exchange of
! information.
!------------------------------------------------------------------------------

!t Loop over all US points, and for each whose assigned irank is equal to myrank, 
!t flag its WS neighbors for inclusion on this rank.

!t do ius = 2,nus
!t    if (itabg_us(ius)%irank == myrank) then

!t       iw1s = ltab_us(ius)%iw(1)
!t       iw2s = ltab_us(ius)%iw(2)

!t       myrankflag_ws(iw1s) = .true.
!t       myrankflag_ws(iw2s) = .true.

!t    endif
!t enddo

! Loop over all WS points, and for each whose seaflag value is .true., flag
! the WS point for inclusion on this rank.  Then, for any WS point flagged, 
! flag its US and MS neighbors for inclusion on this rank.

do iws = 2,nws
   if (seaflag(iws)) then
      myrankflag_ws(iws) = .true.
   endif

   if (myrankflag_ws(iws)) then

      do jms = 1,ltab_ws(iws)%npoly
         ims = ltab_ws(iws)%im(jms)
         ius = ltab_ws(iws)%iu(jms)

         myrankflag_ms(ims) = .true.
         myrankflag_us(ius) = .true.
      enddo

   endif
enddo

! Loop over all MS, US, and WS points and count the ones that have been flagged
! for inclusion on this rank.

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = ims_myrank + 1
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = ius_myrank + 1
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1
   endif
enddo

! Set mms, mus, mws values for this rank

mms = ims_myrank
mus = ius_myrank
mws = iws_myrank

! Re-allocate itab data structures and main grid coordinate arrays

call alloc_sea_grid(mms,mus,mws)

! Reset point counts to 1

ims_myrank = 1
ius_myrank = 1
iws_myrank = 1

! Store new myrank MS, US, WS indices in itabg data structures

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = ims_myrank + 1

      itabg_ms(ims)%ims_myrank = ims_myrank      
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = ius_myrank + 1

      itabg_us(ius)%ius_myrank = ius_myrank
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1

      itabg_ws(iws)%iws_myrank = iws_myrank
   endif
enddo

! Memory copy to main tables

do ims = 2,nms
   if (myrankflag_ms(ims)) then
      ims_myrank = itabg_ms(ims)%ims_myrank

      itab_ms(ims_myrank)%imglobe = ims

      sea%xem  (ims_myrank) = sea_t%xem  (ims)
      sea%yem  (ims_myrank) = sea_t%yem  (ims)
      sea%zem  (ims_myrank) = sea_t%zem  (ims)
      sea%zm   (ims_myrank) = sea_t%zm   (ims)
      sea%glatm(ims_myrank) = sea_t%glatm(ims)
      sea%glonm(ims_myrank) = sea_t%glonm(ims)
   endif
enddo

do ius = 2,nus
   if (myrankflag_us(ius)) then
      ius_myrank = itabg_us(ius)%ius_myrank

      itab_us(ius_myrank)%irank = itabg_us(ius)%irank
      itab_us(ius_myrank)%iuglobe = ius

      im1s = ltab_us(ius)%im(1)
      im2s = ltab_us(ius)%im(2)
      iw1s = ltab_us(ius)%iw(1)
      iw2s = ltab_us(ius)%iw(2)

      if (myrankflag_ms(im1s)) itab_us(ius_myrank)%im(1) = itabg_ms(im1s)%ims_myrank
      if (myrankflag_ms(im2s)) itab_us(ius_myrank)%im(2) = itabg_ms(im2s)%ims_myrank
      if (myrankflag_ws(iw1s)) itab_us(ius_myrank)%iw(1) = itabg_ws(iw1s)%iws_myrank
      if (myrankflag_ws(iw2s)) itab_us(ius_myrank)%iw(2) = itabg_ws(iw2s)%iws_myrank
   endif
enddo

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = itabg_ws(iws)%iws_myrank

      itab_ws(iws_myrank)%irank   = itabg_ws(iws)%irank
      itab_ws(iws_myrank)%iwglobe = iws
      itab_ws(iws_myrank)%npoly   = ltab_ws(iws)%npoly

      sea%leaf_class  (iws_myrank) = sea_t%leaf_class(iws)

      sea%area (iws_myrank) = sea_t%area (iws)
      sea%xew  (iws_myrank) = sea_t%xew  (iws)
      sea%yew  (iws_myrank) = sea_t%yew  (iws)
      sea%zew  (iws_myrank) = sea_t%zew  (iws)
      sea%glatw(iws_myrank) = sea_t%glatw(iws)
      sea%glonw(iws_myrank) = sea_t%glonw(iws)

       do jms = 1,ltab_ws(iws)%npoly
         ims = ltab_ws(iws)%im(jms)
         ius = ltab_ws(iws)%iu(jms)

         if (myrankflag_ms(ims)) itab_ws(iws_myrank)%im(jms) =  &
                                 itabg_ms(ims)%ims_myrank

         if (myrankflag_us(ius)) itab_ws(iws_myrank)%iu(jms) =  &
                                 itabg_us(ius)%ius_myrank
      enddo

   endif
enddo

!------------------------------------------------------------------------------
! The following section (which would correspond to landcell-to-landcell exchange
! of information between parallel processors) is likely not required for sea 
! cells.  The prognostic ocean model will perform its own parallel exchange of
! information.
!------------------------------------------------------------------------------

!t Loop over all US points and get the 2 neighbor WS indices of each

!t do ius = 2,nus
!t    iw1s = ltab_us(ius)%iw(1)
!t    iw2s = ltab_us(ius)%iw(2)

!t    if (itabg_us(ius)%irank == myrank) then

! US point is on myrank.  If either WS neighbor is on remote rank, myrank must
! receive that WS cell.

!t       if (iw1s > 1 .and. itabg_ws(iw1s)%irank /= myrank) then
!t          call recv_table_ws(itabg_ws(iw1s)%irank)
!t       endif

!t       if (iw2s > 1 .and. itabg_ws(iw2s)%irank /= myrank) then
!t          call recv_table_ws(itabg_ws(iw2s)%irank)
!t       endif

!t    elseif (itabg_us(ius)%irank /= myrank) then

! US point is on remote rank.  If either WS neighbor is on myrank, myrank must
! send that WS cell.

!t       if (iw1s > 1 .and. itabg_ws(iw1s)%irank == myrank) then
!t          call send_table_ws(iw1s,itabg_us(ius)%irank)
!t       endif

!t       if (iw2s > 1 .and. itabg_ws(iw2s)%irank == myrank) then
!t          call send_table_ws(iw2s,itabg_us(ius)%irank)
!t       endif

!t    endif

!t enddo

! Loop over all seaflux cells
! (This section is required for computing seacell-atmosphere fluxes and 
! involves only surface and canopy properties.)

do isf = 2,nseaflux

   iw  = seaflux_temp(isf)%iw   ! full-domain index
   iws = seaflux_temp(isf)%iwls ! full-domain index

   if (itabg_w(iw)%irank /= myrank .and. itabg_ws(iws)%irank == myrank) then  

! ATM cell is on remote rank while SEA cell is on myrank, so seaflux is computed
! on remote rank and must be received by myrank, and SEA cell fields must be
! sent to remote rank.

      call send_table_ws(iws,itabg_w(iw)%irank)
      call recv_table_wsf(itabg_w(iw)%irank)

   endif

   if (itabg_w(iw)%irank == myrank .and. itabg_ws(iws)%irank /= myrank) then  

! ATM cell is on myrank while SEA cell is on remote rank, so seaflux is computed
! on myrank and must be sent to remote rank, and SEA cell fields must be
! received from remote rank.

      call send_table_wsf(isf,itabg_ws(iws)%irank)
      call recv_table_ws(itabg_ws(iws)%irank)

   endif

enddo

! Deallocate temporary data structures and arrays

deallocate (ltab_ms, ltab_us, ltab_ws)
deallocate (sea_t%leaf_class)

deallocate (sea_t%area, sea_t%xew, sea_t%yew, sea_t%zew)
deallocate (sea_t%glatw, sea_t%glonw)

deallocate (sea_t%xem, sea_t%yem, sea_t%zem, sea_t%zm)
deallocate (sea_t%glatm, sea_t%glonm)

end subroutine para_init_sea

!===============================================================================

subroutine recv_table_ws(iremote)

use mem_para,  only: nrecvs_ws, recv_ws
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_ws(1) .and. recv_ws(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_ws(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_ws(1).

if (jrecv > nrecvs_ws(1)) nrecvs_ws(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_ws(jrecv)%iremote = iremote

return
end subroutine recv_table_ws

!===============================================================================

subroutine recv_table_wsf(iremote)

use mem_para,  only: nrecvs_wsf, recv_wsf
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

jrecv = 1
do while (jrecv <= nrecvs_wsf(1) .and. recv_wsf(jrecv)%iremote /= iremote)
   jrecv = jrecv + 1
enddo

! If jrecv exceeds nrecvs_wsf(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wsf(1).

if (jrecv > nrecvs_wsf(1)) nrecvs_wsf(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wsf(jrecv)%iremote = iremote

return
end subroutine recv_table_wsf

!===============================================================================

subroutine send_table_ws(iws,iremote)

use mem_sea,  only: itab_ws, itabg_ws
use mem_para,  only: nsends_ws, send_ws
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: iws
integer, intent(in) :: iremote

integer :: jsend
integer :: iws_myrank

! Check whether iremote is already in table of ranks to send to

jsend = 1
do while (jsend <= nsends_ws(1) .and. send_ws(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_ws(1), jsend represents a rank not yet entered in the
! table, so increase nsends_ws(1).

if (jsend > nsends_ws(1)) nsends_ws(1) = jsend

! If nsends_ws(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_ws, nsends_ws(1) exceeds maxremote'
   write(io6,*) 'nsends_ws(1) = ',nsends_ws(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_ws '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iws_myrank = itabg_ws(iws)%iws_myrank

itab_ws(iws_myrank)%send(jsend) = .true.

send_ws(jsend)%iremote = iremote

return
end subroutine send_table_ws

!===============================================================================

subroutine send_table_wsf(isf,iremote)

use mem_sflux, only: seaflux, seafluxg
use mem_para,  only: nsends_wsf, send_wsf
use max_dims,  only: maxremote
use misc_coms, only: io6

implicit none

integer, intent(in) :: isf
integer, intent(in) :: iremote

integer :: jsend
integer :: isf_myrank

! Check whether iremote is already in table of ranks to send to

!write(io6,*) 'pbs1 ',isf,iremote

jsend = 1
do while (jsend <= nsends_wsf(1) .and. send_wsf(jsend)%iremote /= iremote)
   jsend = jsend + 1
enddo

! If jsend exceeds nsends_wsf(1), jsend represents a rank not yet entered in the
! table, so increase nsends_wsf(1).

if (jsend > nsends_wsf(1)) nsends_wsf(1) = jsend

! If nsends_wsf(1) exceeds maxremote, print error message and stop

if (jsend > maxremote) then
   write(io6,*) 'In subroutine send_table_wsf, nsends_wsf(1) exceeds maxremote'
   write(io6,*) 'nsends_wsf(1) = ',nsends_wsf(1)
   write(io6,*) 'maxremote = ',maxremote
   write(io6,*) 'Stopping model '
   stop 'stopping in send_table_wsf '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

isf_myrank = seafluxg(isf)%isf_myrank

seaflux(isf_myrank)%sendf(jsend) = .true.

send_wsf(jsend)%iremote = iremote

return
end subroutine send_table_wsf

