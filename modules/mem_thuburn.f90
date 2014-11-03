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

Module mem_thuburn
  implicit none

  real,    allocatable :: scp_local_min(:,:)
  real,    allocatable :: scp_local_max(:,:)
  real,    allocatable :: cfl_out_sum(:,:)
  real,    allocatable :: cfl_vin(:,:)
  real,    allocatable :: cfl_win(:,:)
  real,    allocatable :: c_scp_in_max_sum(:,:)
  real,    allocatable :: c_scp_in_min_sum(:,:)
  real,    allocatable :: scp_out_min(:,:)
  real,    allocatable :: scp_out_max(:,:)
  real,    allocatable :: scp_in_min(:,:)
  real,    allocatable :: scp_in_max(:,:)
  real,    allocatable :: tfact(:,:)
  integer, allocatable :: kdepv(:,:)

Contains

!===============================================================================

  subroutine alloc_thuburn(imonot, mza, mva, mwa)
    use misc_coms, only: rinit
    implicit none
    
    integer, intent(in) :: imonot, mza, mva, mwa

    if (imonot == 1) then
       
       allocate( scp_local_min   (mza,mwa)) ; scp_local_min    = rinit
       allocate( scp_local_max   (mza,mwa)) ; scp_local_max    = rinit
       allocate( cfl_out_sum     (mza,mwa)) ; cfl_out_sum      = rinit
       allocate( cfl_vin         (mza,mva)) ; cfl_vin          = rinit
       allocate( c_scp_in_max_sum(mza,mwa)) ; c_scp_in_max_sum = rinit
       allocate( c_scp_in_min_sum(mza,mwa)) ; c_scp_in_min_sum = rinit
       allocate( scp_out_min     (mza,mwa)) ; scp_out_min      = rinit
       allocate( scp_out_max     (mza,mwa)) ; scp_out_max      = rinit
       allocate( cfl_win         (mza,mwa)) ; cfl_win          = rinit
       allocate( scp_in_min      (mza,mwa)) ; scp_in_min       = rinit
       allocate( scp_in_max      (mza,mwa)) ; scp_in_max       = rinit
       allocate( tfact           (mza,mwa)) ; tfact            = rinit
       allocate( kdepv           (mza,mva)) ; kdepv            = 0

    endif

  end subroutine alloc_thuburn

!===============================================================================

  subroutine dealloc_thuburn()
    implicit none

    if (allocated( scp_local_min    )) deallocate( scp_local_min )
    if (allocated( scp_local_max    )) deallocate( scp_local_max )
    if (allocated( cfl_out_sum      )) deallocate( cfl_out_sum )
    if (allocated( cfl_vin          )) deallocate( cfl_vin )
    if (allocated( c_scp_in_max_sum )) deallocate( c_scp_in_max_sum )
    if (allocated( c_scp_in_min_sum )) deallocate( c_scp_in_min_sum )
    if (allocated( scp_out_min      )) deallocate( scp_out_min )
    if (allocated( scp_out_max      )) deallocate( scp_out_max )
    if (allocated( cfl_win          )) deallocate( cfl_win )
    if (allocated( scp_in_min       )) deallocate( scp_in_min )
    if (allocated( scp_in_max       )) deallocate( scp_in_max )
    if (allocated( tfact            )) deallocate( tfact )
    if (allocated( kdepv            )) deallocate( kdepv )

  end subroutine dealloc_thuburn

End Module mem_thuburn
