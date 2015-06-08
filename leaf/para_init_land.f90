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
subroutine para_init_land()

  use misc_coms,  only: io6
  use mem_ijtabs, only: itabg_w
  use mem_para,   only: myrank
  use leaf_coms,  only: mwl, nwl, nzg
  use mem_leaf,   only: itab_wl, alloc_land_grid, itab_wl_pd_iw
                        
  implicit none

  integer :: iw, iwl
  integer :: iwl_myrank ! Counter for WL points to be included on this rank

! Loop over all WL points and count the ones that have been flagged
! for inclusion on this rank.

! Land and sea cells are included on this node if the corresponding
! atmospheric cell is primary on this node

  iwl_myrank = 1

  do iwl = 2,nwl
     iw = itab_wl_pd_iw(iwl)
     if (itabg_w(iw)%irank == myrank) then
        iwl_myrank = iwl_myrank + 1
     endif
  enddo

  mwl = iwl_myrank

  ! Allocate itab data structures and main grid coordinate arrays

  call alloc_land_grid(mwl,nzg)

  iwl_myrank = 1
  itab_wl(1)%iwglobe = 1
  itab_wl(1)%iw      = 1

  do iwl = 2,nwl
     iw = itab_wl_pd_iw(iwl)
     if (itabg_w(iw)%irank == myrank) then
        iwl_myrank = iwl_myrank + 1
        itab_wl(iwl_myrank)%iwglobe = iwl
        itab_wl(iwl_myrank)%iw      = iw
     endif
  enddo

  ! Read land grid information

  call landfile_read()

  ! Deallocate temporary arrays

  deallocate(itab_wl_pd_iw)
  
end subroutine para_init_land
