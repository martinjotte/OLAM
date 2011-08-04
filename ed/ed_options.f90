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
Module ed_options

use max_dims, only: pathlen

integer :: ied_init_mode ! 0 signifies a start from near-bare ground.   
                         ! 1 allows you to start from version 1 database files.

integer :: istoma_scheme ! 0 means do the full stomatal conductance, photosynthesis calculation each time.  1 means do the abbreviated approximate calculation.

integer :: ianth_disturb ! 0 means no anthropogenic disturbances such as agriculture and forestry.

integer :: iphen_scheme ! 0 means the old predictive scheme.  Other schemes will be coming soon.

character(pathlen) :: ed_inputs_dir

character(pathlen) :: ed_offline_db ! File specifying where to obtain the 
! meteorological driver data for an offline run.

integer :: n_plant_lim   ! Determines whether (1) or not (0) plants can
                         ! be limited by nitrogen
integer :: n_decomp_lim  ! Determines whether (1) or not (0) decomposition
                         ! can be limited by nitrogen

integer :: include_fire  ! Determines whether (1) or not (0) fire the simulation will include fire.


integer :: ied_offline   ! Run ED offline (1) or coupled to OLAM (0)

real :: frq_rad_ol  ! For an offline run, frequency at which short 
! wave radiation is updated [s].
real :: frq_met_ol  ! For an offline run, frequency at which all met 
! excluding short wave radiation is updated [s].

integer :: metcyc1 ! For an offline run, first year of meteorological forcing
integer :: metcyc2 ! For an offline run, last year of meteorological forcing

real, parameter :: frq_phenology = 86400.0  ! time scale on which phenology
!and daily dynamics are updated.  You could change this, but note that 
!structural dynamics occurs on the first of every month.  So, phenology
!and structural dynamics could potentially get out of sync, which would be bad.

real :: runoff_time          ! In offline runs (and maybe online, too),
                             ! the forest could get saturated and 
                             ! develop a sfcwater pool of meters or 
                             ! more.  In many cases, this is not realistic.
                             ! This is the time scale in seconds for 
                             ! this water to disappear.         

character(pathlen) :: ed_hfilin

end Module ed_options
