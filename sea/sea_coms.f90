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
Module sea_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

   use max_dims,    only: pathlen
   use consts_coms, only: r8

   implicit none

   ! Don't re-export symbols from other modules
   private :: pathlen, r8

   ! SST VARIABLES
   integer            :: nsstfiles, isstfile, isstflg
   integer            :: iupdsst, isstcyclic
   character(pathlen) :: sst_database

   character(pathlen), allocatable :: fnames_sst  (:)
   character(14),      allocatable :: ctotdate_sst(:)
   real(r8),           allocatable :: s1900_sst   (:)

   ! SEAICE VARIABLES
   integer            :: nseaicefiles, iseaicefile, iseaiceflg
   integer            :: iupdseaice, iseaicecyclic
   character(pathlen) :: seaice_database

   character(pathlen), allocatable :: fnames_seaice  (:)
   character(14),      allocatable :: ctotdate_seaice(:)
   real(r8),           allocatable :: s1900_seaice   (:)

   real :: seatmp     ! default sea sfc temperature [K]
   real :: seaice     ! default sea ice fraction [0-1]
   real :: dt_sea     ! sea timestep [s]

   integer, parameter :: nzi    = 3      ! max number of seaice layers
   real,    parameter :: emi    = 1.0    ! emissivity of ice
   real,    parameter :: emw    = 1.0    ! emissivity of water
   real,    parameter :: t00sea = 271.38 ! Freezing temperature of sea water [K]

End Module sea_coms

