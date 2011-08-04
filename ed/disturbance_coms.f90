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
Module disturbance_coms

! GENERAL PARAMETERS
!--------------------------

integer, parameter :: patch_dynamics = 1  !  Set to 1 to incorporate the effects of disturbance, and to do patch fusion.

integer, parameter :: n_dist_types = 3 ! Number of disturbance types.  Disturbance type 1 corresponds to agriculture, disturbance type 2 corresponds to secondary forest, disturbance type 3 corresponds to primary forest.

real, parameter :: min_new_patch_area = 0.005  !  minimum fractional area required to form a new patch.

integer, parameter :: num_lu_trans = 19 ! number of different types of land use transitions in George Hurtt's GLU data set.

! TREEFALL DISTURBANCE
!--------------------------

real :: treefall_disturbance_rate ! Rate (1/years) at which treefall gaps form.  Read from OLAMIN.

real, parameter :: treefall_hite_threshold = 10.0  !  Only trees above this height create a gap when they fall.

real, parameter :: treefall_age_threshold = 0.0  !  Minimum patch age for treefall disturbance.

! FORESTRY
!--------------------------

integer, parameter :: forestry_on = 0  ! Set to 1 if to do forest harvesting.
integer, parameter :: agriculture_on = 0  ! Set to 1 if to do agriculture.
integer, parameter :: plantation_year = 1960 ! Earliest year at which plantations occur
real, parameter :: plantation_rotation = 25.0 ! Number of years that a plantation requires to reach maturity
real, parameter :: mature_harvest_age = 50.0 ! Years that a non-plantation patch requires to reach maturity

! FIRE
!--------------------------

integer, parameter :: fire_on = 1  ! Set to 1 if to do fire
real, parameter :: fire_dryness_threshold = 0.2  !  (meters) Fire can occur if total soil water falls below this threshold.
real, parameter :: fire_parameter = 1.0  ! Dimensionless parameter controlling speed of fire spread.



type lutime
   integer :: landuse_year  ! the year
   real, dimension(num_lu_trans) :: landuse  ! the landuse information (e.g., area harvested, biomass harvested, area abandoned, area converted to agriculture, etc.)
   type(lutime), pointer :: next_lutime ! pointer to the next year
end type lutime

end Module disturbance_coms
