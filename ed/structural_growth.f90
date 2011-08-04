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
subroutine structural_growth(cs, month)

  use ed_structure_defs
  use pft_coms, only: q, qsw, seedling_mortality, c2n_leaf, c2n_storage,   &
       c2n_recruit, c2n_stem, l2n_stem
  use decomposition_coms, only: f_labile

  implicit none

  type(site) :: cs
  integer, intent(in) :: month
  type(patch), pointer :: ed_patch
  type(cohort), pointer :: cc
  real :: salloc
  real :: salloci
  real :: balive_in
  real :: bdead_in
  real :: hite_in
  real :: dbh_in
  real :: nplant_in
  real :: bstorage_in
  real :: f_bseeds
  real :: f_bdead
  real :: balive_mort_litter
  real :: bstorage_mort_litter
  real :: struct_litter
  real :: mort_litter
  real :: seed_litter
  real :: net_seed_N_uptake
  real :: net_stem_N_uptake
  integer :: update_month
  real :: cb_act
  real :: cb_max
  integer :: imonth

  ! Initialization
  cs%basal_area(1:n_pft, 1:n_dbh) = 0.0
  cs%agb(1:n_pft, 1:n_dbh) = 0.0

  ! Loop over patches
  ed_patch => cs%oldest_patch
  do while(associated(ed_patch))

     ! Loop over cohorts
     cc => ed_patch%tallest
     do while(associated(cc))

        salloc = 1.0 + q(cc%pft) + qsw(cc%pft) * cc%hite
        salloci = 1.0 / salloc

        ! Remember inputs in order to calculate increments later on
        balive_in = cc%balive
        bdead_in = cc%bdead
        hite_in = cc%hite
        dbh_in = cc%dbh
        nplant_in = cc%nplant
        bstorage_in = cc%bstorage

        call plant_structural_allocation(cc, month, f_bseeds, f_bdead)

        ! Grow plants
        cc%bdead = cc%bdead + f_bdead * cc%bstorage

        ! Rebalance the plant nitrogen uptake considering the actual 
        ! allocation to structural growth.
        net_stem_N_uptake = (cc%bdead - bdead_in) * ( 1.0 / c2n_stem &
             - 1.0 / c2n_storage) * cc%nplant

        ! Calculate total seed production and seed litter
        cc%bseeds = f_bseeds * cc%bstorage
        seed_litter = cc%bseeds * cc%nplant * seedling_mortality(cc%pft)

        ! Rebalance the plant nitrogen uptake considering the actual 
        ! allocation to seeds.
        net_seed_N_uptake = cc%bseeds * (1.0 / c2n_recruit(cc%pft) &
             - 1.0 / c2n_storage) * cc%nplant

        ! Adjust storage pool
        cc%bstorage = cc%bstorage * (1.0 - f_bdead - f_bseeds)

        ! Apply mortality and reset tracker
        if(cc%nplant + cc%monthly_dndt < 3.0e-8)cc%monthly_dndt = 3.0e-8 -   &
             cc%nplant

        cc%nplant = cc%nplant + cc%monthly_dndt
        ! Plants suffering mortality do not disperse seeds.
        seed_litter = seed_litter - cc%bseeds * cc%monthly_dndt * (1.0 -   &
             seedling_mortality(cc%pft))
        cc%monthly_dndt = 0.0

        ! Calculate contribution to litter pools.
        balive_mort_litter = cc%balive * (nplant_in - cc%nplant)
        bstorage_mort_litter = cc%bstorage * (nplant_in - cc%nplant)
        struct_litter = cc%bdead * (nplant_in - cc%nplant)
        mort_litter = balive_mort_litter + bstorage_mort_litter + struct_litter

        ! Finalize litter inputs
        ed_patch%fsc_in = ed_patch%fsc_in + f_labile(cc%pft) *   &
             balive_mort_litter + bstorage_mort_litter + seed_litter
        ed_patch%fsn_in = ed_patch%fsn_in + f_labile(cc%pft) *   &
             balive_mort_litter / c2n_leaf(cc%pft) +   &
             bstorage_mort_litter/ c2n_storage + seed_litter /  &
             c2n_recruit(cc%pft)
        ed_patch%ssc_in = ed_patch%ssc_in + (1.0 - f_labile(cc%pft)) *  &
             balive_mort_litter + struct_litter
        ed_patch%ssl_in = ed_patch%ssl_in + ( (1.0 - f_labile(cc%pft)) *  &
             balive_mort_litter + struct_litter ) * l2n_stem / c2n_stem
        cc%patchptr%total_plant_nitrogen_uptake = &
             cc%patchptr%total_plant_nitrogen_uptake   &
             + net_seed_N_uptake + net_stem_N_uptake

        ! Calculate the derived cohort properties
        call update_derived_cohort_props(cc)

        ! Update annual average carbon balances for mortality
        update_month = month - 1
        if(update_month == 0)update_month = 12
        cc%cb(update_month) = cc%cb(13)
        cc%cb_max(update_month) = cc%cb_max(13)
        cc%cb(13) = 0.0
        cc%cb_max(13) = 0.0
        cb_act = 0.0
        cb_max = 0.0
        do imonth = 1,12
           cb_act = cb_act + cc%cb(imonth)
           cb_max = cb_max + cc%cb_max(imonth)
        enddo
        if(cb_max > 0.0)then
           cc%cbr_bar = cb_act / cb_max
        else
           cc%cbr_bar = 0.0
        endif

        ! Update interesting output quantities
        call update_vital_rates(cc%siteptr, cc, dbh_in, bdead_in,   &
             balive_in, hite_in, bstorage_in, nplant_in, mort_litter)

        cc => cc%shorter
     enddo
           
     ! Age the patch
     if(ed_patch%dist_type /= 1)ed_patch%age = ed_patch%age + 1.0/12.0
     
     ed_patch => ed_patch%younger
  enddo

  return
end subroutine structural_growth

!===============================================================

subroutine plant_structural_allocation(cc, month, f_bseeds, f_bdead)

  use pft_coms, only: phenology, repro_min_h, r_fract
  use ed_structure_defs

  implicit none

  type(cohort) :: cc
  integer, intent(in) :: month
  real, intent(out) :: f_bseeds
  real, intent(out) :: f_bdead

  ! Calculate fraction of bstorage going to bdead and reproduction
  if(phenology(cc%pft) /= 2   .or.  &  ! for NOT broad leaf deciduous
       (cc%siteptr%lat >= 0.0 .and. month == 6) .or.  &  ! or Jun in north
       (cc%siteptr%lat < 0.0 .and. month == 12) )then    ! or Dec in south
     
     ! For all PFTs except broadleaf deciduous
     
     if(cc%hite <= repro_min_h(cc%pft))then
        f_bseeds = 0.0
     else
        f_bseeds = r_fract(cc%pft)
     endif
     f_bdead = 1.0 - f_bseeds
     
  else
     
     f_bdead = 0.0
     f_bseeds = 0.0
     
  endif
        
  return
end subroutine plant_structural_allocation

!===================================================================

subroutine update_derived_cohort_props(cc)

  use ed_structure_defs
  use pft_coms, only: phenology, sla, q, qsw
  use mem_leaf, only: land

  implicit none

  type(cohort) :: cc
  real :: bd2dbh
  real :: dbh2h
  real :: bl
  real :: bl_max
  real :: dbh2bl
  real :: rootdepth
  real :: calc_root_depth
  integer, external :: assign_root_depth

  cc%dbh = bd2dbh(cc%pft, cc%bdead) 
  cc%hite = dbh2h(cc%pft, cc%dbh)
     
  if(cc%phenology_status /= 2)then

     ! Update status
     bl = cc%balive * cc%siteptr%green_leaf_factor(cc%pft) / (1.0 +   &
          q(cc%pft) + qsw(cc%pft) * cc%hite)
     bl_max = dbh2bl(cc%dbh,cc%pft) * cc%siteptr%green_leaf_factor(cc%pft)
     if(bl.lt.bl_max)then
        cc%phenology_status = 1
     else
        cc%phenology_status = 0
     endif
     
     ! Update LAI
     cc%lai = bl * cc%nplant * sla(cc%pft)
     cc%bleaf = bl

  endif

  ! Update rooting depth
  rootdepth = calc_root_depth(cc%hite,cc%dbh,cc%pft)
  
  ! See which discrete soil level this corresponds to
  cc%krdepth = assign_root_depth(rootdepth, cc%patchptr,   &
       land%lsl(cc%siteptr%iland))
  

  return
end subroutine update_derived_cohort_props

!=====================================================================

subroutine update_vital_rates(cs, cc, dbh_in, bdead_in, balive_in, hite_in,  &
     bstorage_in, nplant_in, mort_litter)
  
  use ed_structure_defs

  implicit none

  real, intent(in) :: dbh_in
  real, intent(in) :: bdead_in
  real, intent(in) :: balive_in
  real, intent(in) :: hite_in
  real, intent(in) :: bstorage_in
  real, intent(in) :: nplant_in
  real, intent(in) :: mort_litter
  type(site) :: cs
  type(cohort) :: cc
  integer :: bdbh
  real, external :: ed_agb

  ! Only update for cohorts on the first census.
  if(cc%first_census /= 1)return

  bdbh = min(int(dbh_in*0.1),10)+1

  ! Computed for plants alive both at past census and current census
  cs%basal_area_growth(cc%pft,bdbh) = cs%basal_area_growth(cc%pft,bdbh) +  &
       cc%patchptr%area * cc%nplant * 3.1415 * 0.25 * (cc%dbh**2 - dbh_in**2)
  cs%agb_growth(cc%pft,bdbh) = cs%agb_growth(cc%pft,bdbh) +  &
       cc%patchptr%area * cc%nplant * (ed_agb(cc%bdead, cc%balive,   &
       cc%bleaf, cc%pft, cc%hite, cc%bstorage) -  &
       ed_agb(bdead_in, balive_in, cc%bleaf, cc%pft,  &
       hite_in, bstorage_in))
  ! note bleaf is unchanged.
  
  ! Computed for plants alive at past census but dead at current census
  cs%basal_area_mort(cc%pft,bdbh) = cs%basal_area_mort(cc%pft,bdbh) +  &
       cc%patchptr%area * (nplant_in - cc%nplant) * 3.1415 * 0.25 * dbh_in**2
        
  cs%agb_mort(cc%pft,bdbh) = cs%agb_mort(cc%pft,bdbh) +  &
       cc%patchptr%area * mort_litter

  return
end subroutine update_vital_rates
        
!====================================================================

subroutine update_patch_derived_props(cp)
  
  use ed_structure_defs

  implicit none

  type(patch) :: cp
  real :: norm_fac
  type(cohort), pointer :: cc
  real :: ba

  ! call derived patch-level structural quantities.  These depend
  ! on the results from reproduction, which in turn depends on 
  ! structural growth results from all patches.

  ! Reset height
  norm_fac = 0.0
  cp%veg_height = 0.0
  cp%lai = 0.0
  
  ! Loop over cohorts
  cc => cp%tallest
  do while(associated(cc))
     
     ! Compute contribution to height
     ba = cc%nplant * cc%dbh**2
     norm_fac = norm_fac + ba
     cp%veg_height = cp%veg_height + cc%hite * ba
     
     ! Update LAI
     cp%lai = cp%lai + cc%lai
     
     cc => cc%shorter
  enddo
  
  ! Update vegetation height
  if(norm_fac > 0.0)then
     cp%veg_height = cp%veg_height / norm_fac
  else
     ! this branch if there aren't any cohorts
     cp%veg_height = 0.2
  endif
  cp%veg_rough = 0.13 * cp%veg_height

  return
end subroutine update_patch_derived_props

!====================================================================

subroutine update_site_derived_props(cs, census_flag)
  
  use ed_structure_defs

  implicit none

  type(site)            :: cs
  type(patch),  pointer :: cp
  type(cohort), pointer :: cc
  real :: ba
  integer :: bdbh
  real, external :: ed_agb
  integer, intent(in) :: census_flag


  ! Loop over patches
  cp => cs%oldest_patch
  do while(associated(cp))
     
     ! Loop over cohorts
     cc => cp%tallest
     do while(associated(cc))

        ! Update basal area, agb
        if(census_flag == 0 .or. cc%first_census == 1)then
           bdbh = min( int(cc%dbh * 0.1), 10) + 1
           ba = cc%nplant * cc%dbh**2
           cs%basal_area(cc%pft, bdbh) = cs%basal_area(cc%pft, bdbh) +   &
                cp%area * ba * 3.1415 * 0.25
           cs%agb(cc%pft, bdbh) = cs%agb(cc%pft, bdbh) +  &
                ed_agb(cc%bdead, cc%balive, cc%bleaf, cc%pft, cc%hite,   &
                cc%bstorage) * cc%nplant
        endif

        cc => cc%shorter
     enddo
     
     cp => cp%younger
  enddo

  return
end subroutine update_site_derived_props
