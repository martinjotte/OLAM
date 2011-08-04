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
subroutine ed_init_atm()
use ed_options, only: ied_init_mode, ied_offline
use ed_structure_defs
use mem_leaf, only: land, first_site
use leaf_coms, only: nzg, nzs, slcpd, soil_rough
use consts_coms, only: alli1000, cliq1000

implicit none
type(site), pointer :: cs
type(patch), pointer :: cp
type(cohort), pointer :: cc
integer :: k
real :: soil_tempc
integer :: nls
integer :: nlsw1

if(.not.associated(first_site))return

if(ied_init_mode == 0 .or. ied_init_mode == 2)then

   !  If we are here, then we don't have any restart information 
   !  available for the surface biophysical properties.

   ! Loop over sites
   cs => first_site
   do while(associated(cs))

      ! Loop over patches
      cp => cs%oldest_patch
      do while(associated(cp))

         ! Canopy air properties
         if(ied_offline == 0)then
            cp%can_temp = land%can_temp(cs%iland)
            cp%can_shv = land%can_shv(cs%iland)
         else
            cp%can_temp = cs%metinput%atm_tmp
            cp%can_shv = cs%metinput%atm_shv
         endif
         cp%can_depth = land%can_depth(cs%iland)

         ! Allocate sfcwater state variables
         allocate(cp%sfcwater_mass(nzs))
         allocate(cp%sfcwater_energy(nzs))
         allocate(cp%sfcwater_depth(nzs))
         allocate(cp%rshort_s(nzs))
         allocate(cp%rshort_s_beam(nzs))
         allocate(cp%rshort_s_diffuse(nzs))
         allocate(cp%sfcwater_tempk(nzs))
         allocate(cp%sfcwater_fracliq(nzs))

         ! Initialize sfcwater state variables
         cp%sfcwater_mass(1) = land%sfcwater_mass(1,cs%iland)
         cp%sfcwater_energy(1) = land%sfcwater_energy(1,cs%iland)
         cp%sfcwater_depth(1) = land%sfcwater_depth(1,cs%iland)
         cp%nlev_sfcwater = land%nlev_sfcwater(cs%iland)
         do k = 1,cp%nlev_sfcwater
            call qtk(cp%sfcwater_energy(k), cp%sfcwater_tempk(k),  &
                 cp%sfcwater_fracliq(k))
         enddo

         ! Allocate soil state variables
         allocate(cp%ntext_soil(nzg))
         allocate(cp%soil_energy(nzg))
         allocate(cp%soil_water(nzg))
         allocate(cp%soil_tempk(nzg))
         allocate(cp%soil_fracliq(nzg))

         ! Initialize soil properties
         do k = 1,nzg
            cp%ntext_soil(k) = land%ntext_soil(k,cs%iland)
            cp%soil_water(k) = land%soil_water(k,cs%iland)
            if(ied_offline == 0)then
               cp%soil_energy(k) = land%soil_energy(k,cs%iland)
            else
               soil_tempc = cp%can_temp - 273.15
               cp%soil_energy(k) = soil_tempc * (slcpd(cp%ntext_soil(k)) +   &
                    cp%soil_water(k) * cliq1000) + cp%soil_water(k) * alli1000
            endif
            call qwtk(cp%soil_energy(k), cp%soil_water(k)*1.0e3,  &
                 slcpd(cp%ntext_soil(k)), cp%soil_tempk(k), cp%soil_fracliq(k))
         enddo

         cp%rough = soil_rough
         if(ied_offline == 0) then
            cp%ground_shv = land%ground_shv(cs%iland)
            cp%surface_ssh = land%surface_ssh(cs%iland)
         else

            nls   = cp%nlev_sfcwater
            nlsw1 = max(nls,1)
   
            call grndvap(cs%iland,                                    &
                 nls,                                      &
                 cp%ntext_soil       (nzg),  &
                 cp%soil_water       (nzg),  &
                 cp%soil_energy      (nzg),  &
                 cp%sfcwater_energy(nlsw1),  &
                 land%rhos(cs%iland)                 ,  &
                 cp%can_shv              ,  &
                 cp%ground_shv           ,  &
                 cp%surface_ssh          )

         endif

         ! Initialize vegetation properties
         cc => cp%tallest
         do while(associated(cc))

            cc%hcapveg = land%hcapveg(cs%iland)
            if(ied_offline == 0)then
               cc%veg_temp = land%can_temp(cs%iland)
            else
               cc%veg_temp = cs%metinput%atm_tmp
            endif
            cc%veg_water = 0.0

            cc => cc%shorter
         enddo

         call update_patch_derived_props(cp)

         cp => cp%younger
      enddo

      call update_site_derived_props(cs,0)

      cs => cs%next_site
   enddo
endif

return
end subroutine ed_init_atm

!======================================================================

subroutine init_ed_site_vars(cs)

  use ed_structure_defs
  use mem_ed, only: n_pft, n_dbh
  use ed_options, only: ied_offline, frq_rad_ol, frq_met_ol
  
  implicit none

  type(site) :: cs

  cs%basal_area_growth(1:n_pft, 1:n_dbh) = 0.0
  cs%basal_area_mort(1:n_pft, 1:n_dbh) = 0.0
  cs%basal_area_cut(1:n_pft, 1:n_dbh) = 0.0
  cs%basal_area_recruit(1:n_pft, 1:n_dbh) = 0.0

  cs%agb_growth(1:n_pft, 1:n_dbh) = 0.0
  cs%agb_mort(1:n_pft, 1:n_dbh) = 0.0
  cs%agb_cut(1:n_pft, 1:n_dbh) = 0.0
  cs%agb_recruit(1:n_pft, 1:n_dbh) = 0.0

  allocate(cs%green_leaf_factor(n_pft))
  cs%green_leaf_factor(1:n_pft) = 1.0

  allocate(cs%leaf_aging_factor(n_pft))
  cs%leaf_aging_factor(1:n_pft) = 1.0

  cs%min_monthly_temp = 0.0

  cs%mmean_gpp = 0.0
  cs%mmean_plresp = 0.0
  cs%mmean_rh = 0.0
  cs%mmean_nep = 0.0

  cs%omean_precip = 0.0
  cs%omean_qprecip = 0.0
  cs%omean_netrad = 0.0

  cs%lambda_fire(1:12) = 0.0

  cs%disturbance_memory(1:n_dist_types, 1:n_dist_types) = 0.0

  cs%agri_stocking_density = 10.0

  if(cs%lat > 40.0)then
     cs%agri_stocking_pft = 5
     cs%plantation_stocking_pft = 6
  else
     cs%agri_stocking_pft = 1
     cs%plantation_stocking_pft = 7
  endif
  cs%plantation_stocking_density = 4.0

  cs%primary_harvest_memory = 0.0
  cs%secondary_harvest_memory = 0.0

  return
end subroutine init_ed_site_vars

!======================================================================

subroutine init_ed_patch_vars(cp)

  use ed_structure_defs
  use mem_ed, only: n_pft
  
  implicit none

  type(patch) :: cp
  integer :: count_cohorts
  integer :: ipft

  allocate(cp%A_o_max(n_pft))
  allocate(cp%A_c_max(n_pft))
  allocate(cp%old_stoma_data_max(n_pft))

  do ipft = 1,n_pft
     cp%old_stoma_data_max(ipft)%recalc = 1
  enddo

  cp%avg_daily_temp = 0.0

  cp%omean_rh = 0.0
  cp%omean_nep = 0.0

  cp%omean_runoff = 0.0
  cp%omean_wflux = 0.0
  cp%omean_latflux = 0.0
  cp%omean_qrunoff = 0.0
  cp%omean_hflux = 0.0

  cp%dmean_A_decomp = 0.0
  cp%dmean_Af_decomp = 0.0

  cp%cohort_count = count_cohorts(cp)

  cp%repro(1:n_pft) = 0.0

  return
end subroutine init_ed_patch_vars

!======================================================================

subroutine init_ed_cohort_vars(cc, cp)

  use mem_leaf, only: land
  use ed_structure_defs
  
  implicit none

  type(patch)  :: cp
  type(cohort) :: cc
  real :: root_depth
  real :: calc_root_depth
  integer, external :: assign_root_depth

  cc%omean_gpp = 0.0
  cc%omean_leaf_resp = 0.0
  cc%omean_root_resp = 0.0
  
  cc%dmean_leaf_resp = 0.0
  cc%dmean_root_resp = 0.0
  cc%dmean_gpp = 0.0
  cc%dmean_gpp_pot = 0.0
  cc%dmean_gpp_max = 0.0
  
  cc%growth_respiration = 0.0
  cc%storage_respiration = 0.0
  cc%vleaf_respiration = 0.0
  
  cc%fsn = 1.0
  cc%monthly_dndt = 0.0
  
  nullify(cc%old_stoma_data)
  allocate(cc%old_stoma_data)
  cc%old_stoma_data%recalc = 1

  cc%Psi_open = 0.0

  cc%cb(1:13) = 1.0
  cc%cb_max(1:13) = 1.0
  cc%cbr_bar = 1.0

  root_depth = calc_root_depth(cc%hite,cc%dbh, cc%pft)
  cc%krdepth = assign_root_depth(root_depth, cp, land%lsl(cc%siteptr%iland))
  cc%first_census = 0
  cc%new_recruit_flag = 0

  return
end subroutine init_ed_cohort_vars

