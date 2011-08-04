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
Module mem_para

  integer :: mgroupsize
  integer :: myrank

  integer, allocatable :: nsends_u(:)  ! dimensioned to mrls
  integer, allocatable :: nsends_v(:)  ! dimensioned to mrls
  integer, allocatable :: nsends_w(:)  ! dimensioned to mrls

  integer, allocatable :: nrecvs_u(:)  ! dimensioned to mrls
  integer, allocatable :: nrecvs_v(:)  ! dimensioned to mrls
  integer, allocatable :: nrecvs_w(:)  ! dimensioned to mrls

! LAND and SEA sends/receives are scalars

! integer, allocatable :: nsends_ul(:)   ! dimensioned to mrls
  integer, allocatable :: nsends_wl(:)   ! dimensioned to mrls
  integer, allocatable :: nsends_wlf(:)  ! dimensioned to mrls

! integer, allocatable :: nrecvs_ul(:)   ! dimensioned to mrls
  integer, allocatable :: nrecvs_wl(:)   ! dimensioned to mrls
  integer, allocatable :: nrecvs_wlf(:)  ! dimensioned to mrls

! integer, allocatable :: nsends_us(:)   ! dimensioned to mrls
  integer, allocatable :: nsends_ws(:)   ! dimensioned to mrls
  integer, allocatable :: nsends_wsf(:)  ! dimensioned to mrls

! integer, allocatable :: nrecvs_us(:)   ! dimensioned to mrls
  integer, allocatable :: nrecvs_ws(:)   ! dimensioned to mrls
  integer, allocatable :: nrecvs_wsf(:)  ! dimensioned to mrls

  Type nodebuffs
     character, allocatable :: buff(:)
     integer :: nbytes  =  0
     integer :: iremote = -1
     integer :: irequest
  End Type nodebuffs

  type(nodebuffs), allocatable :: send_u(:)
  type(nodebuffs), allocatable :: send_v(:)
  type(nodebuffs), allocatable :: send_w(:)
  type(nodebuffs), allocatable :: send_uf(:)
  type(nodebuffs), allocatable :: send_vf(:)

  type(nodebuffs), allocatable :: recv_u(:)
  type(nodebuffs), allocatable :: recv_v(:)
  type(nodebuffs), allocatable :: recv_w(:)
  type(nodebuffs), allocatable :: recv_uf(:)
  type(nodebuffs), allocatable :: recv_vf(:)

! type(nodebuffs), allocatable :: send_ul(:)
  type(nodebuffs), allocatable :: send_wl(:)
  type(nodebuffs), allocatable :: send_wlf(:)

! type(nodebuffs), allocatable :: recv_ul(:)
  type(nodebuffs), allocatable :: recv_wl(:)
  type(nodebuffs), allocatable :: recv_wlf(:)

! type(nodebuffs), allocatable :: send_us(:)
  type(nodebuffs), allocatable :: send_ws(:)
  type(nodebuffs), allocatable :: send_wsf(:)

! type(nodebuffs), allocatable :: recv_us(:)
  type(nodebuffs), allocatable :: recv_ws(:)
  type(nodebuffs), allocatable :: recv_wsf(:)

End Module mem_para

