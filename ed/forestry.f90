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
subroutine apply_forestry(cs, year)

  use ed_structure_defs
  use disturbance_coms, only: lutime, min_new_patch_area, plantation_year
  use lphys_interface,  only: plant_patch
  implicit none

  type(site), target :: cs
  integer, intent(in) :: year
  real :: primary_harvest_target
  real :: secondary_harvest_target
  real :: total_harvest_target
  type(lutime), pointer :: clutime
  real :: total_site_biomass
  type(patch), pointer :: np
  real :: area_mature_primary
  real :: agb_mature_primary
  real :: area_mature_secondary
  real :: agb_mature_secondary
  real :: area_mature_plantation
  real :: agb_mature_plantation
  real :: total_harvested_area
  real :: lambda_mature_primary
  real :: lambda_mature_secondary
  real :: lambda_mature_plantation
  real :: harvest_deficit

  ! The site is harvested up to a target biomass.  Patches are harvested 
  ! first from those above the harvest_age, with equal rates.  If the 
  ! biomass target is not met, the remaining biomass is harvested from the 
  ! patches below the minimum age, starting with the oldest.
  ! Harvest rates are taken from George Hurtt's GLU, Global landuse files. 
  ! Elements 12 and 16 are secondary harvesting and 14 and 18 are primary.

  ! First, find the right year on the linked list.
  clutime => cs%first_lutime
  find_lutime: do while(clutime%landuse_year /= year)
     clutime => clutime%next_lutime
     if(.not.associated(clutime))then
        clutime => cs%last_lutime
        exit find_lutime
     endif
  enddo find_lutime

  ! Set primary and secondary targets based on current rates and unapplied
  ! harvest from previous years (memory)
  primary_harvest_target = clutime%landuse(14) + clutime%landuse(18) +  &
       cs%primary_harvest_memory
  secondary_harvest_target = clutime%landuse(12) + clutime%landuse(16) +  &
       cs%secondary_harvest_memory
  total_harvest_target = primary_harvest_target + secondary_harvest_target

  ! Decide whether or not to create a harvest patch: 
  ! (a) must have site agb > 0 and (b) harvest must exceed some minimum 
  ! threshold.  Note: this is not necessarily the total harvest area.  
  total_site_biomass = sum(cs%agb(1:n_pft, 1:n_dbh))

  if(total_site_biomass == 0.0 .or.  &
       total_harvest_target <= total_site_biomass * min_new_patch_area)then
     ! Update memory and return
     cs%primary_harvest_memory = primary_harvest_target
     cs%secondary_harvest_memory = secondary_harvest_target
     return
  endif

  ! Allocate and initialize the new patch.
  nullify(np)
  allocate(np)
  np%dist_type = 2
  call initialize_disturbed_patch(np)
  nullify(np%tallest)
  nullify(np%shortest)
  np%siteptr => cs

  ! Compute current stocks of agb in mature forests.
  call inventory_mature_forests(cs,   &
       area_mature_primary, agb_mature_primary, area_mature_secondary,  &
       agb_mature_secondary, area_mature_plantation, agb_mature_plantation)

  ! Compute the mature-forest harvest rates
  call mature_forest_harvest_rates(agb_mature_primary,  &
     agb_mature_secondary, agb_mature_plantation, primary_harvest_target,   &
     secondary_harvest_target, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation, harvest_deficit)

  ! Apply harvesting to the mature stands
  call harvest_mature_patches(cs, np, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation)

  ! Compute harvested area from mature patches.  This is also updated
  ! in harvest_immature_patches().
  total_harvested_area = lambda_mature_primary * area_mature_primary +  &
       lambda_mature_secondary * area_mature_secondary +  &
       lambda_mature_plantation * area_mature_plantation

  call harvest_immature_patches(cs, np, harvest_deficit, total_harvested_area)

  ! Now we know the area of the new patch, and can normalize the averaged
  ! patch quantities.

  np%area = total_harvested_area
  call normalize_harvest_patch(np, np%area)

  ! If any patches now have zero area, terminate them.
  call terminate_patches(cs)

  ! Set the pointers for this patch.
  cs%youngest_patch%younger => np
  np%older => cs%youngest_patch
  nullify(np%younger)
  cs%youngest_patch => np

  ! Plant the patch if it is a plantation.
  if(cs%plantation == 1 .and. year > plantation_year)then
     call plant_patch(np, cs%plantation_stocking_pft,   &
          cs%plantation_stocking_density, 2.0)
     np%plantation = 1
  endif

  call update_patch_derived_props(np)
        
  call new_patch_sfc_props(np)
  
  ! Clear out the primary harvest memory.
  cs%primary_harvest_memory = 0.0
  ! There still may be a deficit if we have harvested all of the patch agb.
  cs%secondary_harvest_memory = harvest_deficit

  return
end subroutine apply_forestry

!==============================================================================

subroutine inventory_mature_forests(cs, area_mature_primary,   &
     agb_mature_primary, area_mature_secondary, agb_mature_secondary,  &
     area_mature_plantation, agb_mature_plantation)
  
  use ed_structure_defs
  use disturbance_coms, only: plantation_rotation, mature_harvest_age

  implicit none

  type(site)        :: cs
  real, intent(out) :: area_mature_primary
  real, intent(out) :: agb_mature_primary
  real, intent(out) :: area_mature_secondary
  real, intent(out) :: agb_mature_secondary
  real, intent(out) :: area_mature_plantation
  real, intent(out) :: agb_mature_plantation
  type(patch),  pointer :: cp
  type(cohort), pointer :: cc
  real, external :: ed_agb

  ! Initialize inventory
  area_mature_primary = 0.0
  area_mature_secondary = 0.0
  area_mature_plantation = 0.0
  agb_mature_primary = 0.0
  agb_mature_secondary = 0.0
  agb_mature_plantation = 0.0

  ! Loop over patches
  cp => cs%oldest_patch
  do while(associated(cp))
     
     ! Compute the patch agb
     cp%plant_ag_biomass = 0.0
     cc => cp%tallest
     do while(associated(cc))
        cp%plant_ag_biomass = cp%plant_ag_biomass + ed_agb(cc%bdead,   &
             cc%balive, cc%bleaf, cc%pft, cc%hite, cc%bstorage) * cc%nplant
        cc => cc%shorter
     enddo

     ! Increment appropriate counter
     if(cp%plantation == 1 .and. cp%age > plantation_rotation)then
        
        ! Mature plantation
        area_mature_plantation = area_mature_plantation + cp%area
        agb_mature_plantation = agb_mature_plantation +   &
             cp%plant_ag_biomass * cp%area

     elseif(cp%dist_type == 2 .and. cp%plantation /= 1 .and.   &
          cp%age > mature_harvest_age)then

        ! Mature secondary
        area_mature_secondary = area_mature_secondary + cp%area
        agb_mature_secondary = agb_mature_secondary +   &
             cp%plant_ag_biomass * cp%area

     elseif(cp%dist_type == 3 .and. cp%age > mature_harvest_age)then

        ! Mature primary
        area_mature_primary = area_mature_primary + cp%area
        agb_mature_primary = agb_mature_primary +   &
             cp%plant_ag_biomass * cp%area

     endif

     cp => cp%younger
  enddo

  return
end subroutine inventory_mature_forests

!======================================================================

subroutine mature_forest_harvest_rates(agb_mature_primary,  &
     agb_mature_secondary, agb_mature_plantation, primary_harvest_target,   &
     secondary_harvest_target, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation, harvest_deficit)

  use ed_structure_defs

  implicit none

  real, intent(in) :: agb_mature_primary
  real, intent(in) :: agb_mature_secondary
  real, intent(in) :: agb_mature_plantation
  real, intent(in) :: primary_harvest_target
  real, intent(inout) :: secondary_harvest_target
  real, intent(out) :: lambda_mature_primary
  real, intent(out) :: lambda_mature_plantation
  real, intent(out) :: lambda_mature_secondary
  real, intent(out) :: harvest_deficit

  ! Compute harvesting rate in mature primary forest.  If there is 
  ! not enough biomass to harvest, harvest what is possible from primary
  ! and attempt to harvest the remainder of the target from secondary.
  if(agb_mature_primary > primary_harvest_target)then
     lambda_mature_primary = primary_harvest_target / agb_mature_primary
  else
     lambda_mature_primary = 1.0
     harvest_deficit = primary_harvest_target - agb_mature_primary
     secondary_harvest_target = secondary_harvest_target + harvest_deficit
  endif

  ! Compute harvesting rate in mature plantations and mature secondary
  ! forests.  First try to remove all biomass from plantations.  If this
  ! is not possible, remove what you could and try to remove the remainder
  ! from mature secondary forests.  If this is again not possible, store
  ! what is leaf in harvest_deficit.
  if(agb_mature_plantation > secondary_harvest_target)then
     lambda_mature_plantation = secondary_harvest_target /   &
          agb_mature_plantation
     lambda_mature_secondary = 0.0
     harvest_deficit = 0.0
  else
     lambda_mature_plantation = 1.0
     harvest_deficit = secondary_harvest_target - agb_mature_plantation
     if(agb_mature_secondary > harvest_deficit)then
        lambda_mature_secondary = harvest_deficit / agb_mature_secondary
        harvest_deficit = 0.0
     else
        lambda_mature_secondary = 1.0
        harvest_deficit = harvest_deficit - agb_mature_secondary
     endif
  endif

  return
end subroutine mature_forest_harvest_rates

!======================================================================

subroutine harvest_mature_patches(cs, np, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation)
  
  use ed_structure_defs
  use disturbance_coms, only: mature_harvest_age, plantation_rotation

  implicit none

  type(site)           :: cs
  type(patch)          :: np
  type(patch), pointer :: cp
  real :: dA
  real, intent(in) :: lambda_mature_plantation
  real, intent(in) :: lambda_mature_secondary
  real, intent(in) :: lambda_mature_primary

  ! Loop over patches
  cp => cs%oldest_patch
  do while(associated(cp))
     
     if(cp%plantation == 1 .and. cp%age > plantation_rotation)then
        ! Harvest mature plantations
        dA = cp%area * lambda_mature_plantation
     elseif(cp%dist_type == 2 .and. cp%plantation /= 1 .and.   &
          cp%age > mature_harvest_age)then
        ! Harvest mature secondary
        dA = cp%area * lambda_mature_secondary
     elseif(cp%dist_type == 3 .and. cp%age > mature_harvest_age)then
        ! Harvest mature primary
        dA = cp%area * lambda_mature_primary
     else
        dA = 0.0  ! Immature patches not harvested here.
     endif
     
     if(dA > 0.0)then
        cp%area = cp%area - dA
        call increment_patch_vars(np, cp, dA)
        call accumulate_disturbance_litter(np, cp, 1, dA)
     endif

     cp => cp%younger
  enddo


  return
end subroutine harvest_mature_patches

!====================================================================

subroutine harvest_immature_patches(cs, np, harvest_deficit,   &
     total_harvest_area)
  
  use disturbance_coms, only: plantation_rotation, mature_harvest_age
  use ed_structure_defs

  implicit none

  real, intent(inout)  :: harvest_deficit
  type(site)           :: cs
  type(patch)          :: np
  type(patch), pointer :: cp
  real :: lambda
  real :: dA
  real, intent(inout) :: total_harvest_area

  ! Loop over patches
  cp => cs%oldest_patch
  do while(associated(cp))

     ! First harvest the immature secondary
     if(harvest_deficit > 0.0 .and.  &  ! There is still a deficit
          cp%dist_type == 2 .and.  &  ! Secondary forest
          ( (cp%plantation == 1 .and. cp%age < plantation_rotation) .or.  &
          ! either immature plantation or immature secondary
          (cp%plantation /= 1 .and. cp%age < mature_harvest_age) ) ) then

        if( (cp%area * cp%plant_ag_biomass) > harvest_deficit)then

           ! Patch is not totally harvested
           lambda = harvest_deficit / (cp%area * cp%plant_ag_biomass)
           dA = cp%area * lambda
           harvest_deficit = 0.0

        else

           ! Patch is totally harvested
           dA = cp%area
           harvest_deficit = harvest_deficit - cp%area * cp%plant_ag_biomass

        endif

        total_harvest_area = total_harvest_area + dA
        cp%area = cp%area - dA

        call increment_patch_vars(np, cp, dA)
        call accumulate_disturbance_litter(np, cp, 2, dA)

     endif

     cp => cp%younger
  enddo

  ! Return if we have reached our harvest target.
  if(harvest_deficit <= 0.0)return

  ! If we did not reach our target, loop again through patches, this time 
  ! harvesting from immature primary.
  cp => cs%oldest_patch
  do while(associated(cp))

     ! If necessary, harvest the immature primary
     if(harvest_deficit > 0.0 .and.  &  ! There is still a deficit
          cp%dist_type == 3 .and.    &  ! and this is primary forest
          cp%age < mature_harvest_age)then     ! and it is immature

        if( (cp%area * cp%plant_ag_biomass) > harvest_deficit)then
           
           ! Patch is not totally harvested
           lambda = harvest_deficit / (cp%area * cp%plant_ag_biomass)
           dA = cp%area * lambda
           harvest_deficit = 0.0

        else

           ! Patch is totally harvested
           dA = cp%area
           harvest_deficit = harvest_deficit - cp%area * cp%plant_ag_biomass
           
        endif

        total_harvest_area = total_harvest_area + dA
        cp%area = cp%area - dA

        call increment_patch_vars(np, cp, dA)
        call accumulate_disturbance_litter(np, cp, 2, dA)

     endif

     cp => cp%younger
  enddo

  return
end subroutine harvest_immature_patches

!==================================================================

subroutine normalize_harvest_patch(np, area)

  use ed_structure_defs
  use mem_ed, only: n_pft
  use leaf_coms, only: nzg, nzs

  implicit none

  type(patch)      :: np
  real, intent(in) :: area
  real :: area_fac
  integer :: k

  area_fac = 1.0 / area

  np%fast_soil_C = np%fast_soil_C * area_fac

  np%slow_soil_C = np%slow_soil_C * area_fac

  np%structural_soil_C = np%structural_soil_C * area_fac

  np%structural_soil_L = np%structural_soil_L * area_fac

  np%mineralized_soil_N = np%mineralized_soil_N * area_fac

  np%fast_soil_N = np%fast_soil_N * area_fac

  np%sum_dgd = np%sum_dgd * area_fac

  np%sum_chd = np%sum_chd * area_fac

  np%can_temp = np%can_temp * area_fac

  np%can_shv = np%can_shv * area_fac

  np%can_depth = np%can_depth * area_fac

  do k = 1, nzs
     np%sfcwater_mass(k) = np%sfcwater_mass(k) * area_fac
     np%sfcwater_energy(k) = np%sfcwater_energy(k) * area_fac
     np%sfcwater_depth(k) = np%sfcwater_depth(k) * area_fac
  enddo

  do k = 1, nzg
     np%soil_energy(k) = np%soil_energy(k) * area_fac
     np%soil_water(k) = np%soil_water(k) * area_fac
  enddo

  np%rough = np%rough * area_fac

  np%omean_rh = np%omean_rh * area_fac

  np%dmean_A_decomp = np%dmean_A_decomp * area_fac

  np%dmean_Af_decomp = np%dmean_Af_decomp * area_fac

  np%repro(1:n_pft) = np%repro(1:n_pft) * area_fac

  np%fsc_in = np%fsc_in * area_fac

  np%ssc_in = np%ssc_in * area_fac

  np%ssl_in = np%ssl_in * area_fac

  np%fsn_in = np%fsn_in * area_fac

  np%total_plant_nitrogen_uptake = np%total_plant_nitrogen_uptake * area_fac

  return
end subroutine normalize_harvest_patch

