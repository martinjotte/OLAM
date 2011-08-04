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
Module max_dims

   integer, parameter ::    &

       maxgrds      = 20    & ! Max # of grids
      ,maxremote    = 30    & ! Max # of remote send/recv processes
      ,nzgmax       = 25    & ! Max # of soil levels
      ,maxsndg      = 200   & ! Max # of vertical levels for the input sounding
      ,maxsstfiles  = 2000  & ! Max # of input SST files
      ,maxndvifiles = 2000  & ! Max # of input NDVI files
      ,maxisfiles   = 2000  & ! Max # of input data files for the isan stage
      ,maxisdirs    = 30    & ! Max # of directories that contain data files
      ,maxpr        = 100   & ! Max # of press levels in the input data files
      ,maxnplt      = 150   & ! Max # of fields to plot
      ,maxpltfiles  = 2000  & ! Max # of input files for a plotonly run
      ,maxngrdll    = 20    & ! Max # of geog. pts for each grid refinement
      ,pathlen      = 128   & ! Max length of character strings for file paths
      ,maxnlspoly   = 7       ! max # of M pts for a single land/sea cell

End Module max_dims

