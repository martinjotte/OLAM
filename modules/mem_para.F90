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
Module mem_para

#ifdef OLAM_MPI
   use mpi, only: MPI_REQUEST_NULL
   implicit none
   private     :: MPI_REQUEST_NULL
#else
   implicit none
   integer, parameter, private :: MPI_REQUEST_NULL = 0
#endif

  integer :: mgroupsize
  integer :: myrank

  integer, allocatable :: nsends_v(:)  ! dimensioned to mrls
  integer, allocatable :: nsends_w(:)  ! dimensioned to mrls
  integer, allocatable :: nsends_m(:)  ! dimensioned to mrls

  integer :: nsends_wnud

  integer, allocatable :: nrecvs_v(:)  ! dimensioned to mrls
  integer, allocatable :: nrecvs_w(:)  ! dimensioned to mrls
  integer, allocatable :: nrecvs_m(:)  ! dimensioned to mrls

  integer :: nrecvs_wnud

  Type nodebuffs
     character, allocatable :: buff(:)
     integer :: nbytes   =  0
     integer :: iremote  = -1
     integer :: irequest = MPI_REQUEST_NULL
  End Type nodebuffs

  type(nodebuffs), allocatable :: send_v(:)
  type(nodebuffs), allocatable :: send_w(:)
  type(nodebuffs), allocatable :: send_m(:)
  type(nodebuffs), allocatable :: send_wnud(:)

  type(nodebuffs), allocatable :: recv_v(:)
  type(nodebuffs), allocatable :: recv_w(:)
  type(nodebuffs), allocatable :: recv_m(:)
  type(nodebuffs), allocatable :: recv_wnud(:)

  integer                      :: mva_primary = 0
  integer, target, allocatable :: iva_globe_primary(:)
  integer, target, allocatable :: iva_local_primary(:)

  integer                      :: mwa_primary = 0
  integer, target, allocatable :: iwa_globe_primary(:)
  integer, target, allocatable :: iwa_local_primary(:)

  integer                      :: mma_primary = 0
  integer, target, allocatable :: ima_globe_primary(:)
  integer, target, allocatable :: ima_local_primary(:)

  integer                      :: mwl_primary = 0
  integer, target, allocatable :: iwl_globe_primary(:)
  integer, target, allocatable :: iwl_local_primary(:)

  integer                      :: mws_primary = 0
  integer, target, allocatable :: iws_globe_primary(:)
  integer, target, allocatable :: iws_local_primary(:)

  integer                      :: mwnud_primary = 0
  integer, target, allocatable :: iwnud_globe_primary(:)
  integer, target, allocatable :: iwnud_local_primary(:)

  integer, parameter :: itagv   = 1
  integer, parameter :: itagm   = 2
  integer, parameter :: itagw   = 3
  integer, parameter :: itagwnud= 4

  integer :: nbytes_int
  integer :: nbytes_real
  integer :: nbytes_real8

End Module mem_para
