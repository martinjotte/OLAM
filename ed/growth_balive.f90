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
subroutine dbalive_dt(cs, tfact)

  use ed_structure_defs
  use pft_coms, only: q, qsw, plant_N_supply_scale, c2n_storage
  use ed_options, only: N_plant_lim

  implicit none

  type(site)            :: cs
  type(patch),  pointer :: ed_patch
  type(cohort), pointer :: cc
  real :: salloc
  real :: salloci
  real :: bl
  real :: br
  real :: tfact
  real :: daily_C_gain
  real :: carbon_balance
  real :: carbon_balance_pot
  real :: carbon_balance_max
  real :: balive_in
  real :: nitrogen_supply
  real :: mortality_rates
  real :: dndt
  real :: nitrogen_uptake
  real :: N_uptake_pot

  ed_patch => cs%oldest_patch
  do while(associated(ed_patch))
     
     ! Reset averaged variables
     ed_patch%total_plant_nitrogen_uptake = 0.0
     
     ! Loop over cohorts
     cc => ed_patch%tallest
     do while(associated(cc))

        ! Initialize cohort nitrogen uptake
        nitrogen_uptake = 0.0
        N_uptake_pot = 0.0

        ! Set allocation factors
        salloc = 1.0 + qsw(cc%pft) * cc%hite + q(cc%pft)
        salloci = 1.0 / salloc
        
        ! Calculate leaf, fine root biomass

        if(cc%phenology_status /= 2)then
           bl = cs%green_leaf_factor(cc%pft) * cc%balive * salloci
        else
           bl = 0.0
        endif
        br = q(cc%pft) * cc%balive * salloci 

        call transfer_C_from_storage(cc, salloc, nitrogen_uptake, N_uptake_pot)

        call plant_maintenance_and_resp(cc, br, bl, tfact, daily_C_gain)

        ! Update the Balive
        cc%balive = cc%balive - cc%maintenance_costs

        ! Update the storage
        cc%bstorage = cc%bstorage - cc%storage_respiration

        ! When you lose storage carbon, allow the associated nitrogen 
        ! to go to litter in order to maintain prescribed C2N ratio.
        ed_patch%fsn_in = ed_patch%fsn_in + cc%storage_respiration /   &
             c2n_storage * cc%nplant

        ! Calculate actual, potential and maximum carbon balances
        call plant_carbon_balances(cc, daily_C_gain, carbon_balance,  &
             carbon_balance_pot, carbon_balance_max)
      
        ! Allocate plant carbon balance to balive and bstorage
        balive_in = cc%balive
        call alloc_plant_c_balance(cc, salloc, salloci, carbon_balance,  &
             nitrogen_uptake)

        ! Do a shadow calculation to see what would have happened if stomata 
        ! were open.  This is used to calculate potential nitrogen uptake, 
        ! N_uptake_pot.
        if(N_plant_lim == 1)call potential_N_uptake(cc, salloc, salloci,  &
             balive_in, carbon_balance_pot, N_uptake_pot)

        !  Increment the [kgN/m2] taken up during previous frq_phenology.
        ed_patch%total_plant_nitrogen_uptake =   &
             ed_patch%total_plant_nitrogen_uptake   &
             + nitrogen_uptake * cc%nplant

        ! Calculate plant N limitation factor
        if(n_plant_lim == 0 .or. N_uptake_pot <= 0.0)then
           cc%fsn = 1.0
        else
           br = q(cc%pft) * cc%balive * salloci 
           nitrogen_supply = plant_N_supply_scale * br   &
                * ed_patch%mineralized_soil_N
           cc%fsn = nitrogen_supply / (nitrogen_supply + N_uptake_pot)
        endif
        
        ! Do mortality --- note that only frost mortality changes daily.
        dndt = - mortality_rates(cc) * cc%nplant * tfact

        ! Update monthly mortality rate [plants/m2/month]
        cc%monthly_dndt = cc%monthly_dndt + dndt

        cc => cc%shorter
     enddo
     
     ! Update litter
     call Litter(ed_patch)

     ! Recompute patch LAI
     call update_patch_derived_props(ed_patch)
     
     !reset average daily temperature
     ed_patch%avg_daily_temp = 0.0 
     ed_patch => ed_patch%younger

  enddo

  return
end subroutine dbalive_dt

!====================================================================

subroutine transfer_C_from_storage(cc, salloc, nitrogen_uptake, N_uptake_pot)

  use ed_structure_defs
  use pft_coms, only: c2n_leaf, c2n_storage, c2n_stem
  use decomposition_coms, only: f_labile

  implicit none

  type(cohort)        :: cc
  real, intent(in)    :: salloc
  real, intent(inout) :: nitrogen_uptake
  real, intent(inout) :: N_uptake_pot
  real :: off_allometry_cb
  real :: dbh2bl
  real :: increment

  ! Only do the transfer there are supposed to be leaves

  if(cc%phenology_status < 2)then
     
     ! If plants have storage, transfer it to balive
     off_allometry_cb = dbh2bl(cc%dbh,cc%pft) * salloc - cc%balive
     increment = min( max(0.0, off_allometry_cb), cc%bstorage)
     cc%balive = cc%balive + increment
     cc%bstorage = cc%bstorage - increment

     ! N uptake is required since c2n_leaf < c2n_storage.
     ! Units are kgN/plant/frq_phenology.
     nitrogen_uptake = increment * ( f_labile (cc%pft) / c2n_leaf(cc%pft) +   &
          (1.0 - f_labile(cc%pft)) / c2n_stem - 1.0 / c2n_storage)
     N_uptake_pot = nitrogen_uptake
     
  endif

  return
end subroutine transfer_C_from_storage

!====================================================================

subroutine plant_maintenance_and_resp(cc, br, bl, tfact, daily_C_gain)

  use ed_structure_defs
  use pft_coms, only: phenology, root_turnover_rate, leaf_turnover_rate,  &
       growth_resp_factor, storage_turnover_rate, q, qsw
  use leaf_coms, only: nzg
  use ed_options, only: frq_phenology

  implicit none

  type(cohort)     :: cc
  real, intent(in) :: br
  real, intent(in) :: bl
  real, intent(in) :: tfact
  real, intent(out) :: daily_C_gain
  real :: maintenance_temp_dep

  ! Get the temperature dependence
  if(phenology(cc%pft) /= 1)then
     maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 -   &
          cc%patchptr%soil_tempk(nzg))))
  else
     maintenance_temp_dep = 1.0
  endif

  ! Calculate maintenance demand (kgC/plant/year)
  
  cc%maintenance_costs = (root_turnover_rate(cc%pft) * br +  &
       leaf_turnover_rate(cc%pft) * bl) * maintenance_temp_dep
  
  ! Convert units of maintenance to [kgC/plant/frq_phenology]
  
  cc%maintenance_costs = cc%maintenance_costs * tfact
        
  ! Compute daily C uptake [kgC/plant/frq_phenology]

  daily_C_gain = 1.2e-8 * frq_phenology * (cc%dmean_gpp -   &
       cc%dmean_leaf_resp - cc%dmean_root_resp) / cc%nplant

  ! Compute respiration rates [kgC/plant/frq_phenology]

  cc%growth_respiration = max(0.0, daily_C_gain * growth_resp_factor(cc%pft))
  cc%storage_respiration = cc%bstorage * storage_turnover_rate(cc%pft) * tfact
  cc%vleaf_respiration = (1.0 - cc%siteptr%green_leaf_factor(cc%pft))  &
       / (1.0 + q(cc%pft) + qsw(cc%pft) * cc%hite)  &
       * cc%balive * storage_turnover_rate(cc%pft) * tfact
  
  return
end subroutine plant_maintenance_and_resp

!===================================================================

subroutine plant_carbon_balances(cc, daily_C_gain, carbon_balance,  &
     carbon_balance_pot, carbon_balance_max)

  use ed_structure_defs
  use pft_coms, only: growth_resp_factor
  use ed_options, only: frq_phenology

  implicit none

  type(cohort)      :: cc
  real, intent(in)  :: daily_C_gain
  real, intent(out) :: carbon_balance
  real, intent(out) :: carbon_balance_pot
  real, intent(out) :: carbon_balance_max
  real :: daily_C_gain_pot
  real :: daily_C_gain_max
  real :: growth_respiration_pot
  real :: growth_respiration_max

  ! Calculate actual daily carbon balance: kgC/plant/frq_phenology.

  carbon_balance = daily_C_gain - cc%growth_respiration - cc%vleaf_respiration

  ! Calculate potential carbon balance (used for nitrogen 
  ! demand function).  [kgC/plant/frq_phenology]
  daily_C_gain_pot = 1.2e-8 * frq_phenology * (cc%dmean_gpp_pot -   &
       cc%dmean_leaf_resp - cc%dmean_root_resp) / cc%nplant
  growth_respiration_pot = max(0.0, daily_C_gain_pot *   &
       growth_resp_factor(cc%pft))
  carbon_balance_pot = daily_C_gain_pot - growth_respiration_pot -   &
       cc%vleaf_respiration
        
  ! Calculate maximum carbon balance (used for mortality) 
  daily_C_gain_max = 1.2e-8 * frq_phenology * (cc%dmean_gpp_max -   &
       cc%dmean_leaf_resp - cc%dmean_root_resp) / cc%nplant
  growth_respiration_max = max(0.0, daily_C_gain_max *   &
       growth_resp_factor(cc%pft))
  carbon_balance_max = daily_C_gain_max - growth_respiration_max -   &
       cc%vleaf_respiration
        
  ! Carbon balances for mortality
  cc%cb(13) = cc%cb(13) + carbon_balance - cc%maintenance_costs
  cc%cb_max(13) = cc%cb_max(13) + carbon_balance_max -   &
       cc%maintenance_costs

  return
end subroutine plant_carbon_balances

!====================================================================

subroutine alloc_plant_c_balance(cc, salloc, salloci, carbon_balance,   &
     nitrogen_uptake)

  use ed_structure_defs
  use pft_coms, only: c2n_storage, c2n_leaf, sla, c2n_stem
  use decomposition_coms, only: f_labile

  implicit none

  type(cohort)     :: cc
  real, intent(in) :: salloc
  real, intent(in) :: salloci
  real, intent(in) :: carbon_balance
  real, intent(inout) :: nitrogen_uptake
  real :: bl_max
  real :: bl_pot
  real :: dbh2bl
  real :: increment

  if(cc%phenology_status == 0 .and. carbon_balance > 0.0 )then

     ! Simply update monthly carbon gain.  This will be 
     ! used for structural growth at the end of the month.
     cc%bstorage = cc%bstorage + carbon_balance
     nitrogen_uptake = nitrogen_uptake +  &
          carbon_balance / c2n_storage
     cc%bleaf = cc%balive * salloci * cc%siteptr%green_leaf_factor(cc%pft)

  else

     ! are there leaves?
     if(cc%phenology_status < 2)then

        bl_max = dbh2bl(cc%dbh,cc%pft) * cc%siteptr%green_leaf_factor(cc%pft)
        bl_pot = cc%siteptr%green_leaf_factor(cc%pft)  &
             * (cc%balive + carbon_balance) * salloci

        ! will this increment take us over the limit?
        if(bl_pot > bl_max)then

           ! if so, put remainder in storage
           increment = carbon_balance   &
                - (dbh2bl(cc%dbh,cc%pft) * salloc - cc%balive)
           cc%bstorage = cc%bstorage + increment
           nitrogen_uptake = nitrogen_uptake +  &
                increment / c2n_storage
           increment = dbh2bl(cc%dbh,cc%pft) * salloc - cc%balive
           cc%balive = cc%balive + increment
           nitrogen_uptake = nitrogen_uptake +  &
                increment * (f_labile(cc%pft) / c2n_leaf(cc%pft) +  &
                (1.0 - f_labile(cc%pft)) / c2n_stem)
           cc%bleaf = bl_max
           cc%phenology_status = 0

        else

           ! it will not exceed limit, so just add to balive
           cc%balive = cc%balive + carbon_balance
           cc%phenology_status = 1
           cc%bleaf = cc%balive * salloci *   &
                cc%siteptr%green_leaf_factor(cc%pft)
           if(carbon_balance < 0.0)then
              cc%patchptr%fsn_in = cc%patchptr%fsn_in - carbon_balance *  &
                   (f_labile(cc%pft) / c2n_leaf(cc%pft) + (1.0 -  &
                   f_labile(cc%pft)) / c2n_stem ) * cc%nplant
           else
              nitrogen_uptake = nitrogen_uptake +  &
                   carbon_balance * (f_labile(cc%pft) / c2n_leaf(cc%pft) +  &
                   (1.0 - f_labile(cc%pft)) / c2n_stem)
           endif

        endif

     else

        ! in this case, carbon balance in negative so just subtract
        cc%balive = cc%balive + carbon_balance
        cc%bleaf = 0.0
        cc%lai = 0.0
        cc%patchptr%fsn_in = cc%patchptr%fsn_in - carbon_balance *  &
             (f_labile(cc%pft) / c2n_leaf(cc%pft) + (1.0 -   &
             f_labile(cc%pft)) / c2n_stem) * cc%nplant

     endif
  endif

  cc%lai = cc%bleaf * cc%nplant * sla(cc%pft)

  return
end subroutine alloc_plant_c_balance

!====================================================================

subroutine potential_N_uptake(cc, salloc, salloci, balive_in,   &
     carbon_balance_pot, N_uptake_pot)

  use ed_structure_defs
  use pft_coms, only: c2n_storage, c2n_leaf, c2n_stem
  use decomposition_coms, only: f_labile

  implicit none

  type(cohort)     :: cc
  real, intent(in) :: salloc
  real, intent(in) :: salloci
  real, intent(in) :: balive_in
  real, intent(in) :: carbon_balance_pot
  real, intent(inout) :: N_uptake_pot
  real :: bl_max
  real :: bl_pot
  real :: dbh2bl
  real :: increment

  if(cc%phenology_status == 0 .and. carbon_balance_pot > 0.0 )then

     ! Positive carbon balance with plants fully flushed
     N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage

  else

     if(cc%phenology_status < 2)then

        ! There are at least some leaves
        bl_max = dbh2bl(cc%dbh,cc%pft) * cc%siteptr%green_leaf_factor(cc%pft)
        bl_pot = cc%siteptr%green_leaf_factor(cc%pft)  &
             * (balive_in + carbon_balance_pot) * salloci

        if(bl_pot > bl_max)then

           ! this increment took us over the limit, so remainder is 
           ! put in storage
           increment = carbon_balance_pot   &
                - (dbh2bl(cc%dbh,cc%pft) * salloc - balive_in)
           N_uptake_pot = N_uptake_pot + increment / c2n_storage
           increment = dbh2bl(cc%dbh,cc%pft) * salloc - balive_in
           N_uptake_pot = N_uptake_pot + increment * (f_labile(cc%pft) /   &
                c2n_leaf(cc%pft) + (1.0 - f_labile(cc%pft)) / c2n_stem)
        else

           ! this increment did not exceed the limit.

           if(carbon_balance_pot > 0.0)then

              ! There is uptake if carbon_balance_pot is positive

              N_uptake_pot = N_uptake_pot + carbon_balance_pot *   &
                   (f_labile(cc%pft) / c2n_leaf(cc%pft) + (1.0 -   &
                   f_labile(cc%pft)) / c2n_stem)
           endif

        endif

     endif

  endif

  return
end subroutine potential_N_uptake

!====================================================================

subroutine Litter(cp)

  use ed_structure_defs
  use pft_coms, only: c2n_leaf, c2n_stem, l2n_stem
  use decomposition_coms, only: f_labile

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: cc
  real :: plant_litter
  real :: plant_litter_f
  real :: plant_litter_s

  ! Add fine root and leaf turnover to the litter

  ! Loop over cohorts
  cc => cp%shortest
  do while(associated(cc))

     plant_litter = cc%maintenance_costs * cc%nplant
     plant_litter_f = plant_litter * f_labile(cc%pft)
     plant_litter_s = plant_litter - plant_litter_f

     cp%fsc_in = cp%fsc_in + plant_litter_f
     cp%fsn_in = cp%fsn_in + plant_litter_f / c2n_leaf(cc%pft)

     cp%ssc_in = cp%ssc_in + plant_litter_s
     cp%ssl_in = cp%ssl_in + plant_litter_s * l2n_stem / c2n_stem

     cc => cc%taller
  enddo

  return
end subroutine Litter
