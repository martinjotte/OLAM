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
Module pft_coms

use mem_ed, only: n_pft

!--------------------
! GENERAL
!--------------------
integer, dimension(n_pft) :: include_pft ! Set to 1 if you want to include this PFT; 0 otherwise.

! This flag specifies what PFTs can grow on agriculture patches.  Set 
! to 1 if you want to include this PFT on agriculture patches
integer, dimension(n_pft), parameter :: include_pft_ag = (/  &
     1,  & ! C4 grass
     0,  & ! early successional broadleaf evergreen
     0,  & ! mid successional broadleaf evergreen
     0,  & ! late successional broadleaf evergreen
     1,  & ! C3 grass
     0, & ! northern pines
     0, & ! southern pines
     0, & ! late successional conifers
     0, & ! early successional broadleaf deciduous
     0, & ! early successional broadleaf deciduous
     0  /) ! early successional broadleaf deciduous


!--------------------
! PHOTOSYNTHESIS AND STOMATAL CONDUCTANCE
!--------------------
real, dimension(n_pft) :: D0 ! (mol H2O/mol air).   Stomata begin to rapidly close once the difference between intercellular and boundary layer H2O mixing ratios exceed this value.

real, dimension(n_pft) :: Vm_low_temp ! Temperature (C) below which leaf metabolic activity begins to rapidly decline.

real, dimension(n_pft) :: Vm0  ! maximum photosynthetic capacity at a reference temperature. (umol/m2/s)

real, dimension(n_pft) :: stomatal_slope !  Slope of the Ball/Berry stomatal conductance-photosynthesis relationship.

real, dimension(n_pft) :: cuticular_cond ! Intercept of the Ball/Berry stomatal conductance relationship.  (umol/m2/s)

real, dimension(n_pft) :: quantum_efficiency  ! efficiency of using PAR to fix CO2 (mol CO2 / Einstein)

integer, dimension(n_pft) :: photosyn_pathway ! specifies photosynthetic pathway.  3 corresponds to C3, 4 corresponds to C4.

!--------------------
! RESPIRATION AND TURNOVER
!--------------------
real, dimension(n_pft) :: growth_resp_factor  !  Determines level of growth respiration.  Starting with accumulated photosynthesis (P), leaf (Rl) and root respiration (Rr) are first subtracted.  Then, the growth respiration = (growth_resp_factor) * (P - Rl - Rr).

real, dimension(n_pft) :: leaf_turnover_rate ! This is 1/(leaf life span).  Units are 1/year.

real, dimension(n_pft) :: root_turnover_rate ! This is 1/(fine root life span).  Units are 1/year.

real, dimension(n_pft) :: dark_respiration_factor ! Sets the rate of dark (i.e., leaf) respiration.  (dimensionless; it is relative to Vm0.)

real, dimension(n_pft) :: storage_turnover_rate ! Turnover rate of plant storage pools (1/year).

real, dimension(n_pft) :: root_respiration_factor ! (umol CO2)/(kg fine roots)/second

!--------------------
! MORTALITY AND SURVIVORSHIP
!--------------------
real, dimension(n_pft) :: mort1 ! Sets the time scale at which plants out of carbon balance suffer mortality (1/year)

real, dimension(n_pft) :: mort2 ! Determines how poor the carbon balance needs to be before plants suffer large mortality rates.

real, dimension(n_pft) :: mort3 ! Density-independent mortality rate (1/years)

real, parameter :: frost_mort = 3.0 ! Determines how rapidly trees die if it is too cold for them (1/year)

real, dimension(n_pft) :: seedling_mortality ! Fraction of seedlings that suffer mortality without becoming a recruit. 

real, dimension(n_pft) :: treefall_s_gtht ! Survivorship fraction for trees with heights greater than treefall_hite_threshold (see disturbance_coms.f90)

real, dimension(n_pft) :: treefall_s_ltht ! Survivorship fraction for trees with heights less than treefall_hite_threshold (see disturbance_coms.f90)

real, dimension(n_pft) :: plant_min_temp ! Below this temperature, mortality rapidly increases.

!--------------------
! NITROGEN AND WATER REQUIREMENTS
!--------------------
real, parameter :: c2n_slow = 10.0 ! Carbon to Nitrogen ratio, slow pool.
real, parameter :: c2n_structural = 150.0 ! Carbon to Nitrogen ratio, structural pool.
real, parameter :: c2n_storage = 150.0 ! Carbon to Nitrogen ratio, storage pool.
real, parameter :: c2n_stem = 150.0 ! Carbon to Nitrogen ratio, structural stem.
real, parameter :: l2n_stem = 150.0 ! Carbon to Nitrogen ratio, structural stem.
real, dimension(n_pft) :: c2n_leaf ! Leaf carbon to nitrogen ratio
real, dimension(n_pft) :: c2n_recruit ! Recruit carbon to nitrogen ratio
real, parameter :: C2B = 2.0  !  Carbon-to-biomass ratio of plant tissues.
real, parameter :: agf_bs = 0.7  ! fraction of structural stem that is assumed to be above ground.
real, parameter :: plant_N_supply_scale = 0.5  ! Supply coefficient for plant nitrogen uptake (m2/(kgC fine root)/day)
real, dimension(n_pft) :: water_conductance  ! Supply coefficient for plant water uptake:  (kg H2O) (m2 ground) / (m3 H2O) / (kgC root) / (seconds)

!--------------------
! ALLOCATION AND ALLOMETRY
!--------------------
real, dimension(n_pft) :: rho  ! wood density.  Used only for tropical PFTs (kg/m3).
real, dimension(n_pft) :: SLA ! specific leaf area (m2 leaf / kg C)
real, dimension(n_pft) :: q ! Ratio of (kg fine roots) / (kg leaves)
real, dimension(n_pft) :: qsw ! Ratio of (kg sapwood) / (kg leaves)
real, dimension(n_pft) :: hgt_min ! minimum height of an individual (m)
real, dimension(n_pft) :: b1Ht  !  DBH-height allometry intercept (m).  Temperate PFTs only.
real, dimension(n_pft) :: b2Ht  !  DBH-height allometry slope (1/cm).  Temperate PFTs only.
real, dimension(n_pft) :: b1Bs  !  DBH-stem allometry intercept (kg stem biomass / plant * cm^{-b2Bs}).  Temperate PFTs only.
real, dimension(n_pft) :: b2Bs  !  DBH-stem allometry slope (dimensionless).  Temperate PFTs only.
real, dimension(n_pft) :: b1Bl  !  DBH-leaf allometry intercept (kg leaf biomass / plant * cm^{-b2Bl}).  Temperate PFTs only.
real, dimension(n_pft) :: b2Bl  !  DBH-leaf allometry slope (dimensionless).  Temperate PFTs only.
real, dimension(n_pft) :: max_dbh !  Maximum DBH attainable by this PFT (cm)

!--------------------
! LEAF HABIT AND PHYSICAL PROPERTIES
!--------------------
integer, dimension(n_pft) :: phenology ! Indicates leaf habit.  0 - evergreen coniferous; 1 - drought deciduous; 2 - cold deciduous.

real, dimension(n_pft) :: clumping_factor ! 0-1 factor indicating degree of clumpiness of leaves and shoots.

real, dimension(n_pft) :: leaf_width  ! leaf width used to compute the aerodynamic resistance (m).

!--------------------
! REPRODUCTION AND RECRUITMENT
!--------------------
real, dimension(n_pft) :: r_fract ! fraction of (positive) carbon balance devoted to reproduction.

real, dimension(n_pft) :: seed_rain ! external input of seeds (kgC/m2/year)

real, dimension(n_pft) :: nonlocal_dispersal !  Fraction of seed dispersal that is gridcell-wide.

real, dimension(n_pft) :: repro_min_h ! Minimum height plants need to attain before allocating to reproduction.

End Module pft_coms
