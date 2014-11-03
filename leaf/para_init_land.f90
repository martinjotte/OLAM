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
subroutine para_init_land(landflag)

use max_dims,   only: maxnlspoly
use misc_coms,  only: io6

use mem_ijtabs, only: itabg_w, mrls

use mem_para,   only: mgroupsize, myrank

use leaf_coms,  only: mwl, nwl, nzg

use mem_leaf,   only: itab_wl, itabg_wl, land_vars, land, alloc_land_grid

use mem_mksfc,  only: itab_wls_vars
                        
implicit none

logical, intent(in) :: landflag(nwl)

integer :: iw,iwl

integer :: iwl_myrank = 1 ! Counter for WL points to be included on this rank

! Automatic arrays

logical :: myrankflag_wl(nwl) ! Flag for WL points existing on this rank

! Temporary datatypes

type(itab_wls_vars), allocatable :: ltab_wl(:)
type(land_vars)                  :: land_t

! Move data to temporary data structures, nullifying the old datatype

call move_alloc (itab_wl, ltab_wl)

call move_alloc (land%leaf_class, land_t%leaf_class)
call move_alloc (land%ntext_soil, land_t%ntext_soil)

call move_alloc (land%area , land_t%area )
call move_alloc (land%glatw, land_t%glatw)
call move_alloc (land%glonw, land_t%glonw)
call move_alloc (land%xew  , land_t%xew  )
call move_alloc (land%yew  , land_t%yew  )
call move_alloc (land%zew  , land_t%zew  )
call move_alloc (land%topw , land_t%topw )
call move_alloc (land%wnx  , land_t%wnx  )
call move_alloc (land%wny  , land_t%wny  )
call move_alloc (land%wnz  , land_t%wnz  )

! Initialize myrank flag arrays to .false.

myrankflag_wl(1:nwl) = .false.

! Loop over all WL points, and for each whose landflag value is .true., flag 
! the WL point for inclusion on this rank.

do iwl = 2,nwl
   if (landflag(iwl)) then
      myrankflag_wl(iwl) = .true.
   endif
enddo

! Loop over all WL points and count the ones that have been flagged
! for inclusion on this rank.

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1
   endif
enddo

mwl = iwl_myrank

! Re-allocate itab data structures and main grid coordinate arrays

call alloc_land_grid(mwl,nzg)

! Reset point counts to 1

iwl_myrank = 1

! Store new myrank WL indices in itabg data structures

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = iwl_myrank + 1

      itabg_wl(iwl)%iwl_myrank = iwl_myrank      
   endif
enddo

! Memory copy to main tables

do iwl = 2,nwl
   if (myrankflag_wl(iwl)) then
      iwl_myrank = itabg_wl(iwl)%iwl_myrank

      itab_wl(iwl_myrank)%irank   = itabg_wl(iwl)%irank
      itab_wl(iwl_myrank)%iwglobe = iwl
      itab_wl(iwl_myrank)%iw      = ltab_wl(iwl)%iw
      itab_wl(iwl_myrank)%kw      = ltab_wl(iwl)%kw
      itab_wl(iwl_myrank)%npoly   = ltab_wl(iwl)%npoly
      itab_wl(iwl_myrank)%arf_iw  = ltab_wl(iwl)%arf_iw
      itab_wl(iwl_myrank)%arf_kw  = ltab_wl(iwl)%arf_kw

      itab_wl(iwl_myrank)%xem(1:maxnlspoly) = ltab_wl(iwl)%xem(1:maxnlspoly)
      itab_wl(iwl_myrank)%yem(1:maxnlspoly) = ltab_wl(iwl)%yem(1:maxnlspoly)
      itab_wl(iwl_myrank)%zem(1:maxnlspoly) = ltab_wl(iwl)%zem(1:maxnlspoly)

      land%area (iwl_myrank) = land_t%area (iwl)
      land%glatw(iwl_myrank) = land_t%glatw(iwl)
      land%glonw(iwl_myrank) = land_t%glonw(iwl)
      land%xew  (iwl_myrank) = land_t%xew  (iwl)
      land%yew  (iwl_myrank) = land_t%yew  (iwl)
      land%zew  (iwl_myrank) = land_t%zew  (iwl)
      land%topw (iwl_myrank) = land_t%topw (iwl)
      land%wnx  (iwl_myrank) = land_t%wnx  (iwl)
      land%wny  (iwl_myrank) = land_t%wny  (iwl)
      land%wnz  (iwl_myrank) = land_t%wnz  (iwl)

      land%leaf_class  (iwl_myrank) = land_t%leaf_class  (iwl)
      land%ntext_soil(:,iwl_myrank) = land_t%ntext_soil(:,iwl)
   endif
enddo

! Deallocate temporary data structures and arrays

deallocate(ltab_wl)

deallocate (land_t%leaf_class, land_t%ntext_soil)

deallocate (land_t%area, land_t%xew, land_t%yew, land_t%zew)
deallocate (land_t%glatw, land_t%glonw)
deallocate (land_t%wnx, land_t%wny, land_t%wnz)

end subroutine para_init_land
