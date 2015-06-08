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
subroutine para_init_sea()

  use misc_coms,  only: io6
  use mem_ijtabs, only: itabg_w
  use mem_para,   only: myrank
  use sea_coms,   only: mws, nws
  use mem_sea,    only: itab_ws, alloc_sea_grid, itab_ws_pd_iw

  implicit none

  integer :: iw, iws
  integer :: iws_myrank ! Counter for WS points to be included on this rank

! Loop over all WS points and count the ones that have been flagged
! for inclusion on this rank.

! Land and sea cells are included on this node if the corresponding
! atmospheric cell is primary on this node

  iws_myrank = 1

  do iws = 2,nws
     iw = itab_ws_pd_iw(iws)
     if (itabg_w(iw)%irank == myrank) then
        iws_myrank = iws_myrank + 1
     endif
  enddo

  mws = iws_myrank

  ! Allocate itab data structures and sea grid coordinate arrays

  call alloc_sea_grid(mws)

  iws_myrank = 1
  itab_ws(1)%iwglobe = 1
  itab_ws(1)%iw      = 1

  do iws = 2,nws
     iw = itab_ws_pd_iw(iws)
     if (itabg_w(iw)%irank == myrank) then
        iws_myrank = iws_myrank + 1
        itab_ws(iws_myrank)%iwglobe = iws
        itab_ws(iws_myrank)%iw      = iw
     endif
  enddo

  ! Read sea grid information

  call seafile_read()

  ! Deallocate temporary arrays

  deallocate(itab_ws_pd_iw)
  
end subroutine para_init_sea
