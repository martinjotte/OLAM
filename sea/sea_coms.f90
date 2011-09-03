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
Module sea_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

   use max_dims,    only: maxsstfiles, pathlen
   use consts_coms, only: r8
   implicit none

   private :: maxsstfiles, pathlen, r8

   character(pathlen) :: seafile

   ! SST VARIABLES
   integer            :: nsstfiles, isstfile, isstflg
   integer            :: iupdsst, isstcyclic
   character(pathlen) :: sst_database
   character(pathlen) :: fnames_sst(maxsstfiles)
   character(14)      :: ctotdate_sst(maxsstfiles)
   real(r8)           :: s1900_sst(maxsstfiles)

   ! SEAICE VARIABLES
   integer            :: nseaicefiles, iseaicefile, iseaiceflg
   integer            :: iupdseaice, iseaicecyclic
   character(pathlen) :: seaice_database
   character(pathlen) :: fnames_seaice  (maxsstfiles)
   character(14)      :: ctotdate_seaice(maxsstfiles)
   real(r8)           :: s1900_seaice(maxsstfiles)

   integer :: nms  ! Total # of sea cell M pts in model domain
   integer :: nus  ! Total # of sea cell U pts in model domain
   integer :: nws  ! Total # of sea cell W pts in model domain

   integer :: mms  ! Total # of sea cell M pts in model parallel sub_domain
   integer :: mus  ! Total # of sea cell U pts in model parallel sub_domain
   integer :: mws  ! Total # of sea cell W pts in model parallel sub-domain

   real :: seatmp     ! default sea sfc temp [K]
   real :: dt_sea     ! sea timestep [s]
   
End Module sea_coms

