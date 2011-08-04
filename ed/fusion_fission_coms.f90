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
Module fusion_fission_coms

! Minimum biomass density (kgC/m2) required to form a new recruit.
real, parameter :: min_recruit_size  = 1.0e-3

! Minimum DBH class used in patch profiling
real, parameter :: min_dbh_class = 0.0  

! Maximum DBH (cm) used in patch profiling
real, parameter :: maxdbh = 200.0 

! Number of DBH bins in patch profiling
integer, parameter :: ff_ndbh = 20

! Minimum height class in patch profiling
real, parameter :: min_hgt_class = 0.0

! Cohort fusion tolerance on DBH (dimensionless)
real, parameter :: fusetol = 0.4

! Cohort fusion tolerance on height (m)
real, parameter :: fusetol_h = 0.5

! Cohort fusion tolerance on LAI (m2 leaf/m2 ground)
real, parameter :: lai_fuse_tol = 0.8

! Cohort splitting tolerance on LAI (m2 leaf/ m2 ground)
real, parameter :: lai_tol = 1.0

! Min plant density for height bin to be used in height profile comparisons (plants/m2)
real, parameter :: ntol = 0.001

! Fractional tolerance for patch pft height comparisons (dimensionless)
real, parameter :: profile_tol = 0.2

! Maximum patch age for fusion (years)
real, parameter :: max_patch_age = 500.0

end Module fusion_fission_coms
