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
subroutine apply_disturbances(cs)
  
  use ed_structure_defs
  use misc_coms, only: current_time
  use disturbance_coms, only: n_dist_types, treefall_age_threshold,  &
       min_new_patch_area
  use lphys_interface, only: insert_survivors, plant_patch, split_cohorts, &
                             apply_forestry

  implicit none

  type(site), target :: cs
  integer :: q
  real :: area
  type(patch), pointer :: np
  type(patch), pointer :: cp
  real :: dA
  real :: area_fac
  real, dimension(n_pft, n_dbh) :: initial_agb
  real, dimension(n_pft, n_dbh) :: initial_basal_area

  ! Store AGB, basal area profiles in memory.
  call update_site_derived_props(cs, 1)
  initial_agb(1:n_pft, 1:n_dbh) = cs%agb(1:n_pft, 1:n_dbh)
  initial_basal_area(1:n_pft, 1:n_dbh) = cs%basal_area(1:n_pft, 1:n_dbh)
  
  ! First take care of harvesting: secondary -> secondary and 
  ! primary -> secondary.
  call apply_forestry(cs, current_time%year)
  
  ! Update the cut output variables
  call update_site_derived_props(cs, 1)
  cs%agb_cut(1:n_pft, 1:n_dbh) = cs%agb_cut(1:n_pft, 1:n_dbh) +  &
       initial_agb(1:n_pft, 1:n_dbh) - cs%agb(1:n_pft, 1:n_dbh)
  cs%basal_area_cut(1:n_pft, 1:n_dbh) = cs%basal_area_cut(1:n_pft, 1:n_dbh) + &
       initial_basal_area(1:n_pft, 1:n_dbh) - cs%basal_area(1:n_pft, 1:n_dbh)

  ! Loop over q, the *destination* landuse type
  do q = 1, n_dist_types
     
     ! First, decide if enough area is disturbed to warrant creating
     ! a new patch.
     area = 0.0
     cp => cs%oldest_patch
     do while(associated(cp))
        
        if( (q == 1 .and. cp%dist_type > 1) .or.  & ! conversion to ag
             (q == 2 .and. cp%dist_type == 1) .or.  & ! abandonment
             (q == 3 .and. cp%dist_type > 1 .and.   & ! natural.....
             (cp%age > treefall_age_threshold   &  ! old enough for treefall
             .or. cs%nat_dist_type == 1)) )then  ! or it's a fire
           
           area = area + cp%area * (1.0 - exp(   &
                - (cs%disturbance_rates(q,cp%dist_type)  &
                + cs%disturbance_memory(q,cp%dist_type))))
           
        endif
        cp => cp%younger
     enddo

     if(area > min_new_patch_area)then
!        print*,'making new patch',area,q

        nullify(np)
        allocate(np)
        np%dist_type = q
        np%plantation = 0
        np%area = area
        np%siteptr => cs
        call initialize_disturbed_patch(np)
        
        ! Now go through patches, adding its contribution to the new patch.
        cp => cs%oldest_patch
        do while(associated(cp))
           if( (q == 1 .and. cp%dist_type > 1) .or.  & ! conversion to ag
                (q == 2 .and. cp%dist_type == 1) .or.  & ! abandonment
                (q == 3 .and. cp%dist_type > 1 .and.   & ! natural.....
                (cp%age > treefall_age_threshold   &  ! old enough for treefall
                .or. cs%nat_dist_type == 1)) )then  ! or it's a fire
              
              dA = cp%area * (1.0 - exp(   &
                   - (cs%disturbance_rates(q,cp%dist_type)  &
                   + cs%disturbance_memory(q,cp%dist_type))))
              
              area_fac = dA / np%area
              
              call increment_patch_vars(np, cp, area_fac)
              call insert_survivors(np, cp, q, area_fac)
              call accumulate_disturbance_litter(np, cp, q, area_fac)
              
              ! update patch area
              cp%area = cp%area - dA
              
           endif
           
           cp => cp%younger
        enddo
        
        !  insert the new patch
        nullify(np%younger)
        np%older => cs%youngest_patch
        cs%youngest_patch%younger => np
        cs%youngest_patch => np
        
        ! if the new patch is agriculture, plant it with grasses.
        if(q == 1)call plant_patch(np, cs%agri_stocking_pft,   &
             cs%agri_stocking_density, 1.0)
        
        !  fuse then terminate cohorts
        if(associated(np%tallest))then
           call fuse_cohorts(np)
           call terminate_cohorts(np)
           call split_cohorts(np)
        endif
        
        ! Store AGB, basal area profiles in memory.
        initial_agb(1:n_pft, 1:n_dbh) = cs%agb(1:n_pft, 1:n_dbh)
        initial_basal_area(1:n_pft, 1:n_dbh) = cs%basal_area(1:n_pft, 1:n_dbh)

        ! Update the derived properties including veg_height, lai
        call update_patch_derived_props(np)

        ! Update soil temp, fracliq, etc.
        call new_patch_sfc_props(np)

        ! Update AGB, basal area.
        call update_site_derived_props(cs,1)

        ! Update either cut or mortality
        if(q /= 3)then
           cs%agb_cut(1:n_pft, 1:n_dbh) = cs%agb_cut(1:n_pft, 1:n_dbh) +  &
                initial_agb(1:n_pft, 1:n_dbh) - cs%agb(1:n_pft, 1:n_dbh)
           cs%basal_area_cut(1:n_pft, 1:n_dbh) =   &
                cs%basal_area_cut(1:n_pft, 1:n_dbh) +  &
                initial_basal_area(1:n_pft, 1:n_dbh) -   &
                cs%basal_area(1:n_pft, 1:n_dbh)
        else
           cs%agb_mort(1:n_pft, 1:n_dbh) = cs%agb_mort(1:n_pft, 1:n_dbh) +  &
                initial_agb(1:n_pft, 1:n_dbh) - cs%agb(1:n_pft, 1:n_dbh)
           cs%basal_area_mort(1:n_pft, 1:n_dbh) =   &
                cs%basal_area_mort(1:n_pft, 1:n_dbh) +  &
                initial_basal_area(1:n_pft, 1:n_dbh) -   &
                cs%basal_area(1:n_pft, 1:n_dbh)
        endif
        
        !  clear the disturbance memory for this disturbance type
        cs%disturbance_memory(q,1:n_dist_types) = 0.0
        
     else
        if(area > 0.0)then
           ! the patch creation has been skipped because 
           !    the area was too small 
           ! put the current disturbance rates in memory to be 
           !    added at the next timestep
           cs%disturbance_memory(q,1:n_dist_types) =   &
                cs%disturbance_memory(q,1:n_dist_types) +   &
                cs%disturbance_rates(q,1:n_dist_types)
        endif
     endif
     
  enddo

  return
end subroutine apply_disturbances

!====================================================================

subroutine site_disturbance_rates(month, year, cs)

  use ed_structure_defs
  use ed_options, only: include_fire
  use disturbance_coms, only: treefall_disturbance_rate, lutime
  use pft_coms, only: agf_bs

  implicit none

  integer, intent(in) :: month
  integer, intent(in) :: year
  type(site)          :: cs
  integer :: im
  type(lutime), pointer :: clutime
  real :: fire_dist_rate

  !  calculate fire disturbance rates only if fire is on.
  if(include_fire == 1)then
     fire_dist_rate = sum(cs%lambda_fire(1:12)) / 12.0
  else
     fire_dist_rate = 0.0
  endif
  
  !  treefall disturbance is currently spatiotemporally constant, from OLAMIN
  
  ! For natural disturbance use largest disturbance mode
  if(fire_dist_rate > treefall_disturbance_rate)then
     cs%nat_disturbance_rate = fire_dist_rate
     cs%nat_dist_type = 1
  else
     cs%nat_disturbance_rate = treefall_disturbance_rate
     cs%nat_dist_type = 0
  endif
  
  ! Set disturbance rates assuming only natural disturbance
  cs%disturbance_rates(1:2,1:3)= 0.0
  cs%disturbance_rates(3,1)= 0.0
  cs%disturbance_rates(3,2:3)= cs%nat_disturbance_rate
        
  ! Now it is time for anthropogenic disturbance rates.  First make
  ! sure landuse files exist.
  if(associated(cs%first_lutime))then
     
     ! Only consider years after the start of the dataset.
     ! For years after the end of the dataset, repeat the last year.
     if(year >= cs%first_lutime%landuse_year)then
        
        clutime => cs%first_lutime
        find_lu_year: do while(clutime%landuse_year /= year)
           clutime => clutime%next_lutime
           
           if(.not.associated(clutime))then
              ! If we are at the end of the linked list, use the disturbance
              ! rates of the last year.
              clutime => cs%last_lutime
              exit find_lu_year
           endif
           
        enddo find_lu_year
        
        !    update land-use transition matrix
        
        ! AA2AA, agriculture to agriculture
        cs%disturbance_rates(1,1)= 0.0
        
        ! NA2AA, secondary forest to agriculture
        cs%disturbance_rates(1,2)= clutime%landuse(7) + clutime%landuse(9)
        
        ! NN2AA, primary forest to agriculture
        cs%disturbance_rates(1,3)= clutime%landuse(4) + clutime%landuse(5)
        
        ! AA2NA, agriculture to secondary forest
        cs%disturbance_rates(2,1)= clutime%landuse(8) + clutime%landuse(10)
        
        ! NA2NA, secondary forest to secondary forest (this is taken care
        ! of in the harvesting.)
        cs%disturbance_rates(2,2)= 0.0
        
        ! NN2NA, primary forest to secondary forest (zero here because we
        ! meet a biomass target instead.)
        cs%disturbance_rates(2,3)= 0.0 ! clutime%landuse(11)

        ! AA2NN, agriculture to primary forest (should be zero)
        cs%disturbance_rates(3,1)= 0.0 !clutime%landuse(3) + clutime%landuse(6)
        
        ! NA2NN, secondary forest to primary forest
        cs%disturbance_rates(3,2)= cs%nat_disturbance_rate
        
        ! NN2NN, primary forest to primary forest
        cs%disturbance_rates(3,3)= cs%nat_disturbance_rate
        
     endif ! after first year
  endif ! landuse exists

  ! fraction of above ground litter from disturbance that is 
  ! removed from patch
  cs%loss_fraction(1) = agf_bs
  cs%loss_fraction(2) = agf_bs
  cs%loss_fraction(3) = 0.0

  return
end subroutine site_disturbance_rates

!=====================================================================

subroutine initialize_disturbed_patch(np)
  
  use ed_structure_defs
  use leaf_coms, only: nzs, nzg

  implicit none
  
  type(patch) :: np

  ! dist_type is not set here.

  np%age = 0.0

  ! area is not set here.

  np%fast_soil_C = 0.0

  np%slow_soil_C = 0.0

  np%structural_soil_C = 0.0

  np%structural_soil_L = 0.0

  np%mineralized_soil_N = 0.0

  np%fast_soil_N = 0.0

  np%sum_dgd = 0.0

  np%sum_chd = 0.0

  ! plantation is not set here.

  np%can_temp = 0.0

  np%can_shv = 0.0

  np%can_depth = 0.0

  allocate(np%sfcwater_mass(nzs))
  np%sfcwater_mass(1:nzs) = 0.0

  allocate(np%sfcwater_energy(nzs))
  np%sfcwater_energy(1:nzs) = 0.0

  allocate(np%sfcwater_depth(nzs))
  np%sfcwater_depth(1:nzs) = 0.0

  allocate(np%rshort_s(nzs))

  allocate(np%rshort_s_beam(nzs))

  allocate(np%rshort_s_diffuse(nzs))

  allocate(np%sfcwater_tempk(nzs))
  ! sfcwater_tempk is not set here.

  allocate(np%sfcwater_fracliq(nzs))
  ! sfcwater_fracliq is not set here.

  ! nlev_sfcwater is not set here.

  allocate(np%ntext_soil(nzg))
  ! ntext_soil is not set here.

  allocate(np%soil_energy(nzg))
  np%soil_energy(1:nzg) = 0.0

  allocate(np%soil_water(nzg))
  np%soil_water(1:nzg) = 0.0

  allocate(np%soil_tempk(nzg))
  ! soil_tempk is not set here.

  allocate(np%soil_fracliq(nzg))
  ! soil_fracliq is not set here.

  ! ground_shv is not set here.

  ! surface_ssh is not set here.

  np%rough = 0.0

  call init_ed_patch_vars(np)

  np%fsc_in = 0.0

  np%ssc_in = 0.0

  np%ssl_in = 0.0

  np%fsn_in = 0.0

  np%total_plant_nitrogen_uptake = 0.0

  nullify(np%tallest, np%shortest)

  return
end subroutine initialize_disturbed_patch

!======================================================================

subroutine increment_patch_vars(np, cp, area_fac)

  use ed_structure_defs
  use mem_ed, only: n_pft
  use leaf_coms, only: nzg

  implicit none

  type(patch) :: np
  type(patch) :: cp
  real :: area_fac
  integer :: k

  np%fast_soil_C = np%fast_soil_C + cp%fast_soil_C * area_fac

  np%slow_soil_C = np%slow_soil_C + cp%slow_soil_C * area_fac

  np%structural_soil_C = np%structural_soil_C + cp%structural_soil_C * area_fac

  np%structural_soil_L = np%structural_soil_L + cp%structural_soil_L * area_fac

  np%mineralized_soil_N = np%mineralized_soil_N + cp%mineralized_soil_N *   &
       area_fac

  np%fast_soil_N = np%fast_soil_N + cp%fast_soil_N * area_fac

  np%sum_dgd = np%sum_dgd + cp%sum_dgd * area_fac

  np%sum_chd = np%sum_chd + cp%sum_chd * area_fac

  np%can_temp = np%can_temp + cp%can_temp * area_fac

  np%can_shv = np%can_shv + cp%can_shv * area_fac

  np%can_depth = np%can_depth + cp%can_depth * area_fac

  do k = 1, cp%nlev_sfcwater
     np%sfcwater_mass(k) = np%sfcwater_mass(k) + cp%sfcwater_mass(k) * area_fac
     np%sfcwater_energy(k) = np%sfcwater_energy(k) +   &
          cp%sfcwater_energy(k) * cp%sfcwater_mass(k) * area_fac
     np%sfcwater_depth(k) = np%sfcwater_depth(k) + cp%sfcwater_depth(k) *   &
          area_fac
  enddo

  do k = 1, nzg
     np%ntext_soil(k) = cp%ntext_soil(k)
     np%soil_energy(k) = np%soil_energy(k) +cp%soil_energy(k) * area_fac
     np%soil_water(k) = np%soil_water(k) + cp%soil_water(k) * area_fac
  enddo

  np%rough = np%rough + cp%rough * area_fac

  np%omean_rh = np%omean_rh + cp%omean_rh * area_fac

  np%dmean_A_decomp = np%dmean_A_decomp + cp%dmean_A_decomp * area_fac

  np%dmean_Af_decomp = np%dmean_Af_decomp + cp%dmean_Af_decomp * area_fac

  np%repro(1:n_pft) = np%repro(1:n_pft) + cp%repro(1:n_pft) * area_fac

  np%fsc_in = np%fsc_in + cp%fsc_in * area_fac

  np%ssc_in = np%ssc_in + cp%ssc_in * area_fac

  np%ssl_in = np%ssl_in + cp%ssl_in * area_fac

  np%fsn_in = np%fsn_in + cp%fsn_in * area_fac

  np%total_plant_nitrogen_uptake = np%total_plant_nitrogen_uptake +   &
       cp%total_plant_nitrogen_uptake * area_fac

  return
end subroutine increment_patch_vars

!=======================================================================

subroutine insert_survivors(np, cp, q, area_fac)

  use ed_structure_defs
  use lphys_interface, only: insert_cohort
  implicit none

  type(patch), target :: np
  type(patch)         :: cp
  integer, intent(in) :: q
  real, intent(in) :: area_fac
  type(cohort), pointer :: cc
  real :: n_survivors
  real :: survivorship
  type(cohort), pointer :: nc
  real :: cohort_area_fac

  cc => cp%shortest
  do while(associated(cc))

     cohort_area_fac = survivorship(q, cc) * area_fac
     n_survivors = cc%nplant * cohort_area_fac

     ! If something survived, make a new cohort
     if(n_survivors > 0.0)then

        ! Allocate
        nullify(nc)
        allocate(nc)

        ! Copy
        call copy_cohort(cc, nc)

        ! Adjust area-based variables
        nc%nplant = nc%nplant * cohort_area_fac
        nc%lai = nc%lai * cohort_area_fac
        ! nc%hcapveg = nc%hcapveg * cohort_area_fac
        nc%veg_water = nc%veg_water * cohort_area_fac
        nc%omean_gpp = nc%omean_gpp * cohort_area_fac
        nc%omean_leaf_resp = nc%omean_leaf_resp * cohort_area_fac
        nc%omean_root_resp = nc%omean_root_resp * cohort_area_fac
        nc%growth_respiration = nc%growth_respiration * cohort_area_fac
        nc%storage_respiration = nc%storage_respiration * cohort_area_fac
        nc%vleaf_respiration = nc%vleaf_respiration * cohort_area_fac
        nc%Psi_open = nc%Psi_open * cohort_area_fac

        ! Set pointers
        nc%siteptr => cp%siteptr
        nc%patchptr => np
        call insert_cohort(nc,np)
     endif
     cc => cc%taller
  enddo  ! end loop over cohorts
  
  return
end subroutine insert_survivors

!======================================================================
subroutine accumulate_disturbance_litter(np, cp, q, area_fac)
  
  use ed_structure_defs
  use decomposition_coms, only: f_labile
  use mem_ed, only: n_pft
  use pft_coms, only: c2n_storage, c2n_leaf, c2n_recruit, c2n_stem, l2n_stem

  implicit none

  type(patch) :: np
  type(patch) :: cp
  integer, intent(in) :: q
  real, intent(in) :: area_fac
  real :: fast_litter
  real :: struct_litter
  real :: fast_litter_n
  type(cohort), pointer :: cc
  real :: survivorship

  fast_litter = 0.0
  struct_litter = 0.0
  fast_litter_n = 0.0

  cc => cp%shortest
  do while(associated(cc))

     fast_litter = fast_litter + (1.0 - survivorship(q,cc)) * (  &
          f_labile(cc%pft) * cc%balive + cc%bstorage ) * cc%nplant

     fast_litter_n = fast_litter_n + (1.0 - survivorship(q,cc)) * (  &
          f_labile(cc%pft) * cc%balive / c2n_leaf(cc%pft) + cc%bstorage /  &
          c2n_storage ) * cc%nplant

     struct_litter = struct_litter + cc%nplant *   &
          (1.0 - survivorship(q, cc)) * ( (1.0 -   &
          cp%siteptr%loss_fraction(q)) * cc%bdead +   &
          (1.0 - f_labile(cc%pft)) * cc%balive)
     
     cc => cc%taller
  enddo

  !  Load disturbance litter directly into carbon and N pools
  np%fast_soil_C = np%fast_soil_C + fast_litter * area_fac

  np%structural_soil_C = np%structural_soil_C + struct_litter * area_fac

  np%structural_soil_L = np%structural_soil_L + l2n_stem / c2n_stem *   &
       struct_litter * area_fac

  np%fast_soil_N = np%fast_soil_N + fast_litter_n * area_fac

  return
end subroutine accumulate_disturbance_litter

!=================================================================

subroutine new_patch_sfc_props(np)

  use ed_structure_defs
  use leaf_coms, only: nzg, slcpd, nzs
  use mem_leaf, only: land

  implicit none

  type(patch) :: np
  integer :: k

  do k = 1, nzg
     call qwtk(np%soil_energy(k), np%soil_water(k)*1.e3,  &
          slcpd(np%ntext_soil(k)), np%soil_tempk(k), np%soil_fracliq(k))
  enddo

  np%nlev_sfcwater = 0
  k = 1
  do while(np%sfcwater_mass(min(k,nzs)) > 1.0e-6 .and. k <= nzs)
     np%nlev_sfcwater = k
     np%sfcwater_energy(k) = np%sfcwater_energy(k) / np%sfcwater_mass(k)
     call qtk(np%sfcwater_energy(k), np%sfcwater_tempk(k),   &
          np%sfcwater_fracliq(k))
     k = k+1
  enddo

  call grndvap(np%siteptr%iland,                              &
             np%nlev_sfcwater,          np%ntext_soil  (nzg), &
             np%soil_water     (nzg),   np%soil_energy (nzg), &
             np%sfcwater_energy(max(1,np%nlev_sfcwater)),     &
             land%rhos(np%siteptr%iland),                     &
             np%can_shv,                np%ground_shv,        &
             np%surface_ssh                                   )
             
  return
end subroutine new_patch_sfc_props

!============================================================================

real function survivorship(dest_type, cc)

  use ed_structure_defs
  use disturbance_coms, only: treefall_hite_threshold
  use pft_coms, only: treefall_s_ltht, treefall_s_gtht

  implicit none

  integer, intent(in) :: dest_type
  type(cohort)        :: cc

  if(dest_type == 1)then
     survivorship = 0.0  !  agric
  else
     if(dest_type == 2)then
        survivorship = 0.0   ! secondary land
     else  !   natural land 
        if(dest_type == 3)then
           if(cc%siteptr%nat_dist_type == 1)then
              survivorship = 0.0 ! fire
           elseif(cc%siteptr%nat_dist_type == 0)then
              if(cc%hite < treefall_hite_threshold)then
                 survivorship =  treefall_s_ltht(cc%pft) 
              else 
                 survivorship = treefall_s_gtht(cc%pft)
              endif
           endif
        endif
     endif
  endif

  return
end function survivorship

!=============================================================

subroutine plant_patch(np, pft, density, height_factor)

  use ed_structure_defs
  use pft_coms,        only: q, qsw, sla, hgt_min
  use leaf_coms,       only: dt_leaf
  use lphys_interface, only: insert_cohort
  implicit none

  type(patch), target :: np
  integer, intent(in) :: pft
  real, intent(in) :: density
  real, intent(in) :: height_factor
  type(cohort), pointer :: nc
  real :: h2dbh
  real :: dbh2bd
  real :: dbh2bl

  ! Just make one cohort.  It will soon be split at the splitting call of
  ! apply_disturbances().
  nullify(nc)
  allocate(nc)
  nc%pft = pft
  nc%nplant = density
  nc%hite = hgt_min(nc%pft) * height_factor
  nc%dbh = h2dbh(nc%hite,nc%pft)
  nc%bdead = dbh2bd(nc%dbh,nc%hite,nc%pft)
  nc%bleaf = dbh2bl(nc%dbh,nc%pft)
  nc%phenology_status = 0
  nc%balive = nc%bleaf * (1.0 + q(nc%pft) + qsw(nc%pft) * nc%hite)
  nc%lai = nc%bleaf * nc%nplant * sla(nc%pft)
  nc%bstorage = 0.0
  nc%hcapveg =  3.e4 * max(1.,.025 * dt_leaf)
  nc%veg_temp = np%can_temp
  nc%veg_water = 0.0

  call init_ed_cohort_vars(nc, np)

  nc%new_recruit_flag = 1 ! should plantations be considered recruits?

  nc%siteptr => np%siteptr
  nc%patchptr => np
  call insert_cohort(nc, np)

  return
end subroutine plant_patch
