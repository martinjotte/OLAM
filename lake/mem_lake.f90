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
   ! Colorado State University Relakerch Foundation ; ATMET, LLC 

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
Module mem_lake

  implicit none

  integer :: nzlake    ! Max number of vertical levels in lake model 

  integer :: nlake     ! Total # of lake cell W pts in model domain
  integer :: mlake     ! Total # of lake cell W pts in model parallel sub-domain

  integer :: onlake    ! Offset # of lake cell W pts in model global domain
  integer :: omlake    ! Offset # of lake cell W pts in model parallel sub-domain
                        ! (ilake + omlake = iwsfc)

  real, parameter :: timescale_lake = 2. * 86400. ! Lake cell-to-cell equilibration time scale [s]

! LAKE GRID TABLES

  Type itab_lake_vars
     integer :: iwglobe = 1
  End type

  type (itab_lake_vars), allocatable, target :: itab_lake(:)

! lake MODEL VARIABLES

  Type lake_vars

     real, allocatable :: depth       (:) ! lake water mean depth [m]
     real, allocatable :: lake_energy (:) ! lake energy [J/kg]
     real, allocatable :: surface_srrv(:) ! lake surface sat vapor mixing ratio [kg_vap/kg_dryair]

  End Type lake_vars

  type (lake_vars), target :: lake

Contains

!=========================================================================

  subroutine alloc_lake(mlake)

     use misc_coms, only: rinit

     implicit none

     integer, intent(in) :: mlake

!    Allocate and initialize lake arrays

     allocate (lake%depth       (mlake)) ; lake%depth        = rinit
     allocate (lake%lake_energy (mlake)) ; lake%lake_energy  = rinit
     allocate (lake%surface_srrv(mlake)) ; lake%surface_srrv = rinit

  end subroutine alloc_lake

!=========================================================================

  subroutine filltab_lake()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(lake%depth))        call increment_vtable('LAKE%DEPTH'       , 'RW', rvar1=lake%depth)
     if (allocated(lake%lake_energy))  call increment_vtable('LAKE%LAKE_ENERGY' , 'RW', rvar1=lake%lake_energy)
     if (allocated(lake%surface_srrv)) call increment_vtable('LAKE%SURFACE_SRRV', 'RW', rvar1=lake%surface_srrv)
    
  end subroutine filltab_lake

End Module mem_lake
