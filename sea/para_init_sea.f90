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
subroutine para_init_sea(seaflag)

use max_dims,   only: maxnlspoly
use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank

use sea_coms,   only: mws, nws, iseagrid

use mem_sea,    only: itab_ws, itabg_ws, sea_vars, sea, alloc_sea_grid

use mem_mksfc,  only: itab_wls_vars

implicit none

logical, intent(in) :: seaflag(nws)

integer :: j,jsend
integer :: iw,iws
integer :: isf

integer :: iws_myrank = 1 ! Counter for WS points to be included on this rank

! Automatic arrays

logical :: myrankflag_ws(nws) ! Flag for WS points existing on this rank

! Temporary datatypes

type(itab_wls_vars), allocatable :: ltab_ws(:)
type(sea_vars)                   :: sea_t

! Move data to temporary data structures, nullifying the old datatype

call move_alloc(itab_ws, ltab_ws)

call move_alloc (sea%leaf_class, sea_t%leaf_class)

call move_alloc (sea%area , sea_t%area )
call move_alloc (sea%glatw, sea_t%glatw)
call move_alloc (sea%glonw, sea_t%glonw)
call move_alloc (sea%xew  , sea_t%xew  )
call move_alloc (sea%yew  , sea_t%yew  )
call move_alloc (sea%zew  , sea_t%zew  )
call move_alloc (sea%topw , sea_t%topw )

! Initialize myrank flag arrays to .false.

myrankflag_ws(:) = .false.

! Loop over all WS points, and for each whose seaflag value is .true., flag
! the WS point for inclusion on this rank.

do iws = 2,nws
   if (seaflag(iws)) then
      myrankflag_ws(iws) = .true.
   endif
enddo

! Loop over all WS points and count the ones that have been flagged
! for inclusion on this rank.

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1
   endif
enddo

mws = iws_myrank

! Re-allocate itab data structures and main grid coordinate arrays

print*, 'in para_init_sea1 ',mws

call alloc_sea_grid(mws)

! Reset point counts to 1

iws_myrank = 1

! Store new myrank WS indices in itabg data structures

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = iws_myrank + 1

      itabg_ws(iws)%iws_myrank = iws_myrank
   endif
enddo

! Memory copy to main tables

do iws = 2,nws
   if (myrankflag_ws(iws)) then
      iws_myrank = itabg_ws(iws)%iws_myrank

      itab_ws(iws_myrank)%irank   = itabg_ws(iws)%irank
      itab_ws(iws_myrank)%iwglobe = iws
      itab_ws(iws_myrank)%iw      = ltab_ws(iws)%iw
      itab_ws(iws_myrank)%kw      = ltab_ws(iws)%kw
      itab_ws(iws_myrank)%npoly   = ltab_ws(iws)%npoly
      itab_ws(iws_myrank)%arf_iw  = ltab_ws(iws)%arf_iw
      itab_ws(iws_myrank)%arf_kw  = ltab_ws(iws)%arf_kw

      itab_ws(iws_myrank)%xem(1:maxnlspoly)  = ltab_ws(iws)%xem(1:maxnlspoly)
      itab_ws(iws_myrank)%yem(1:maxnlspoly)  = ltab_ws(iws)%yem(1:maxnlspoly)
      itab_ws(iws_myrank)%zem(1:maxnlspoly)  = ltab_ws(iws)%zem(1:maxnlspoly)

      sea%area (iws_myrank) = sea_t%area (iws)
      sea%glatw(iws_myrank) = sea_t%glatw(iws)
      sea%glonw(iws_myrank) = sea_t%glonw(iws)
      sea%xew  (iws_myrank) = sea_t%xew  (iws)
      sea%yew  (iws_myrank) = sea_t%yew  (iws)
      sea%zew  (iws_myrank) = sea_t%zew  (iws)
      sea%topw (iws_myrank) = sea_t%topw (iws)

      sea%leaf_class(iws_myrank) = sea_t%leaf_class(iws)
   endif
enddo

! Deallocate temporary data structures and arrays

deallocate (ltab_ws)

deallocate (sea_t%leaf_class)

deallocate (sea_t%area, sea_t%xew, sea_t%yew, sea_t%zew)
deallocate (sea_t%glatw, sea_t%glonw)

end subroutine para_init_sea
