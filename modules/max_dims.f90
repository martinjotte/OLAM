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

Module max_dims

  implicit none

  integer, parameter :: maxgrds      = 20   ! Max # of grids
  integer, parameter :: maxremote    = 30   ! Max # of remote send/recv processes
  integer, parameter :: maxsndg      = 200  ! Max # of vertical levels for the input sounding
  integer, parameter :: maxndvifiles = 2000 ! Max # of input NDVI files
  integer, parameter :: maxisdirs    = 30   ! Max # of directories that contain data files
  integer, parameter :: maxpr        = 100  ! Max # of press levels in the input data files
  integer, parameter :: maxnplt      = 150  ! Max # of fields to plot
  integer, parameter :: maxpltfiles  = 2000 ! Max # of input files for a plotonly run
  integer, parameter :: maxngrdll    = 20   ! Max # of geog. pts for each grid refinement
  integer, parameter :: pathlen      = 128  ! Max length of character strings for file paths
  integer, parameter :: maxnlspoly   = 7    ! Max # of M pts for a single land/sea cell
  integer, parameter :: maxlite      = 150  ! Max # of output "lite" variables

End Module max_dims

