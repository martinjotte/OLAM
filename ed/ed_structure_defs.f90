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
Module ed_structure_defs

  use mem_ed, only: n_pft, n_dbh
  use disturbance_coms, only: lutime, n_dist_types
  use fusion_fission_coms, only: ff_ndbh
  use c34constants, only: stoma_data
  use offline_coms

   Type site

      ! (I) INITIALIZING A SITE
      !     When you add a variable to this list, you need to make sure that
      !     it is initialized properly.  init_ed_site_vars() is probably
      !     a good place to do this for most variables.
      ! (II) OUTPUT QUANTITIES
      !     These also are often average variables, and so need to be
      !     normalized before output and re-initialized after output.
      !============================================================

      ! GEOGRAPHICAL INFORMATION
      !  All of the land% (see mem_leaf) arrays are indexed on iland.  
      !  So, each ED site corresponds to exactly one of these land cells.  
      !  Here, iland tells you which one.      
      integer :: iland ! initialized in spawn_sites()

      ! Latitude of the ED site
      real :: lat   ! initialized in spawn_sites()

      ! Longitude of the ED site
      real :: lon  ! initialized in spawn_sites()

      ! Total basal area (m2/ha or cm2/m2)
      real, dimension(n_pft, n_dbh) :: basal_area  ! no init. necessary

      ! Total above ground biomass (tC/ha)
      real, dimension(n_pft, n_dbh) :: agb  ! no init. necessary

      ! ALL GROWTH, MORTALITY, CUT and RECUIT RATES ARE INITIALIZED
      ! IN init_ed_site_vars() and re-initialized in zero_ed_yearly_vars()
      !------------------------------------------------
      ! Growth, basal area metric [m2/ha/y]
      real, dimension(n_pft, n_dbh) :: basal_area_growth

      ! Mortality, basal area metric [m2/ha/y]
      real, dimension(n_pft, n_dbh) :: basal_area_mort

      ! Cut, basal area metric [m2/ha/y]
      real, dimension(n_pft, n_dbh) :: basal_area_cut

      ! Recruitment, basal area metric [m2/ha/y]
      real, dimension(n_pft, n_dbh) :: basal_area_recruit

      ! Growth, AGB metric [tC/ha/y]
      real, dimension(n_pft, n_dbh) :: agb_growth

      ! Mortality, AGB metric [tC/ha/y]
      real, dimension(n_pft, n_dbh) :: agb_mort

      ! Cut, AGB metric [tC/ha/y]
      real, dimension(n_pft, n_dbh) :: agb_cut

      ! Recruitment, AGB metric [tC/ha/y]
      real, dimension(n_pft, n_dbh) :: agb_recruit

      ! In the spring, this is an elongation factor and gives the fraction 
      ! of elongation relative to full leaf size.  In the fall, this is a 
      ! coloration factor and gives the proportion of leaves that are still 
      ! green.  This and leaf_aging_factor are initialized in 
      ! init_ed_site_vars()
      real, allocatable, dimension(:) :: green_leaf_factor

      ! This is a dimensionless aging factor that accounts for declines in 
      ! photosynthetic capacity with leaf age.
      real, allocatable, dimension(:) :: leaf_aging_factor

      ! Minimum daily-average temperature in the current month (K).  Used
      ! to determine whether recruitment can occur.  Initialized in 
      ! init_ed_site_vars().
      real :: min_monthly_temp

      real :: omean_precip
      real :: omean_qprecip
      real :: omean_netrad

      ! All mmean quantities are initialized in init_ed_site_vars() and 
      ! re-initialized in zero_ed_monthly_vars.
      !---------------------------------------------------
      ! Monthly-mean GPP [tC/ha/month]
      real :: mmean_gpp

      ! Monthly-mean plant respiration [tC/ha/month]
      real :: mmean_plresp

      ! Monthly-mean heterotrophic respiration [tC/ha/month]
      real :: mmean_rh

      ! Monthly-mean NEP [tC/ha/month]
      real :: mmean_nep

      ! Monthly-mean fire disturbance rate [1/month].  Initialized in 
      ! init_ed_state_vars().
      real, dimension(12) :: lambda_fire

      ! Effective natural disturbance rate considering treefall and fire [1/y]
      real :: nat_disturbance_rate

      ! Type of natural disturbance to be applied [0=treefall, 1=fire]
      integer :: nat_dist_type

      ! Matrix of disturbance rates [1/y]
      real, dimension(n_dist_types, n_dist_types) :: disturbance_rates

      ! fraction of above ground litter from disturbance that is 
      ! removed from patch
      real, dimension(n_dist_types) :: loss_fraction

      ! number of years of landuse data for this site
      integer :: num_landuse_years

      ! If the area disturbed is very small, but it into memory here instead
      ! of creating a new patch.  Use this to increment the disturbance_rates
      ! next time around. [1/y].  Initialized in init_ed_site_vars().
      real, dimension(n_dist_types, n_dist_types) :: disturbance_memory

      ! Upon creating an agriculture patch in this site, stock it with this 
      ! PFT.  Set, along with the other stocking parameters, in 
      ! init_ed_site_vars().
      integer :: agri_stocking_pft

      ! Upon creating an agriculture patch in this site, stock it with 
      ! this density of plants [plants/m2]
      real :: agri_stocking_density

      ! Upon creating a plantation patch in this site, stock it with this PFT
      integer :: plantation_stocking_pft

      ! Upon creating an plantation patch in this site, stock it with 
      ! this density of plants [plants/m2]
      real :: plantation_stocking_density

      ! Unapplied primary forest harvest from previous years (save until 
      ! harvest is above minimum threshold.) [kgC/m2].  Initialized 
      ! together with secondary memory in init_ed_site_vars().
      real :: primary_harvest_memory

      ! Unapplied secondary forest harvest from previous years (save until 
      ! harvest is above minimum threshold.) [kgC/m2]
      real :: secondary_harvest_memory

      ! In an offline run, store all meteorological driver info here.
      type(met_driv_data) :: metinput

      ! Flag specifying whether (1) or not (0) this site is a plantation
      integer :: plantation

      ! POINTERS
      ! pointer to the oldest patch in the site
      type(patch), pointer :: oldest_patch
      ! pointer to the youngest patch in the site
      type(patch), pointer :: youngest_patch
      ! pointer to the next site
      type(site), pointer :: next_site  
      ! pointer to the first year of landuse information
      type(lutime), pointer :: first_lutime
      ! pointer to the last year of landuse information
      type(lutime), pointer :: last_lutime

   End Type

   Type patch

      ! (I) INITIALIZING A PATCH
      !     When you add a variable to this list, you need to make sure that
      !     it is initialized properly.  
      !    (A) Model restarts
      !       (1) Upon initializing a model run, the restart data are read
      !           in bare_ground_init() or in history restart.
      !       (2) State variables depending on the atmospheric state are
      !           set in ed_init_atm().
      !       (3) Other fundamental variables are set in init_ed_patch_vars()
      !       (4) Derived variables are set in update_patch_derived_props()
      !    (B) New patch from disturbance
      !       (1) When you create a new patch, all quantities are allocated
      !           and initialized to zero in initialize_disturbed_patch().
      !       (2) The contribution from each pre-existing patch to the new
      !           disturbed patch is set in increment_patch_vars().
      !       (3) For forest harvest, you also need to normalize the 
      !           patch variables; this is done in normalize_harvest_patch
      ! (II) DAILY AVERAGES USED IN VEGETATION_DYNAMICS().
      !    (A) These must be initialized according to (I) above.
      !    (B) These must be normalized before use.  This is done in 
      !        normalize_ed_daily_vars().
      !    (C) These must be reinitialized after use.  This is done 
      !        in reinitialize_ed_daily_vars().
      ! (III) OUTPUT QUANTITIES
      !    (A) These must be initialized according to (I) above.
      !    (B) Before write-out, they must be normalized.
      !       (1) For variables written on FRQSTATE, this is done 
      !           in normalized_ed_output_vars().
      !    (C) After write-out, they must be reset to zero.
      !       (1) For variables written on FRQSTATE, this is done 
      !           in zero_ed_output_vars().
      ! (IV) FUSE PATCH
      !      When you fuse patches, in fuse_2_patches(), you may need to 
      !      average your new variable.
      !============================================================

      integer :: dist_type         ! Patch type index:
                                   !       Agriculture = 1
                                   !       Secondary Forest = 2
                                   !       Primary Forest = 3

      ! Time since last disturbance (years)
      real :: age

      ! Fractional area of the site
      real :: area

      ! Soil carbon concentration, fast pool (kg/m2)
      real :: fast_soil_C

      ! Soil carbon concentration, slow pool (kg/m2)
      real :: slow_soil_C

      ! Soil carbon concentration, structural pool (kg/m2)
      real :: structural_soil_C

      ! Soil lignin concentration, structural pool (kg/m2)
      real :: structural_soil_L

      ! Soil nitrogen concentration, mineralized pool (kg/m2)
      real :: mineralized_soil_N

      ! Soil nitrogen concentration, fast pool (kg/m2)
      real :: fast_soil_N

      ! Number of degree days
      ! Degree days --- sum of daily average temperatures above 278.15 K 
      real :: sum_dgd

      ! Number of chilling days
      ! Chill days --- number of days with average temperatures below 278.15 K
      real :: sum_chd

      ! Flag specifying whether (1) or not (0) this patch is a plantation
      integer :: plantation

      ! Temperature (K) of canopy air
      real :: can_temp

      ! Water vapor specific humidity (kg/kg) of canopy air
      real :: can_shv

      ! Canopy depth (m)
      real :: can_depth

      ! Surface water mass (kg/m2)
      real, allocatable, dimension(:) :: sfcwater_mass

      ! Surface water internal energy (J/kg)
      real, allocatable, dimension(:) :: sfcwater_energy

      ! Depth of surface water (m)
      real, allocatable, dimension(:) :: sfcwater_depth

      ! Short wave radiation absorbed by the surface water (W/m2)
      real, allocatable, dimension(:) :: rshort_s

      ! Short wave radiation absorbed by the surface water, 
      ! beam component (W/m2)
      real, allocatable, dimension(:) :: rshort_s_beam

      ! Short wave radiation absorbed by the surface water, 
      ! diffuse component (W/m2)
      real, allocatable, dimension(:) :: rshort_s_diffuse

      ! Temperature of surface water (K)
      real, allocatable, dimension(:) :: sfcwater_tempk

      ! Liquid fraction of surface water
      real, allocatable, dimension(:) :: sfcwater_fracliq

      ! Number of surface water layers
      integer :: nlev_sfcwater

      ! Soil textural class index
      integer, allocatable, dimension(:) :: ntext_soil

      ! Soil internal energy (J/m3)
      real, allocatable, dimension(:) :: soil_energy

      ! Soil water (m3 water / m3 soil)
      real, allocatable, dimension(:) :: soil_water

      ! Temperature of soil (K)
      real, allocatable, dimension(:) :: soil_tempk

      ! Liquid fraction of soil
      real, allocatable, dimension(:) :: soil_fracliq

      ! Effective specific humidity (kg/kg) just above soil
      real :: ground_shv

      ! Surface saturation specific humidity (kg/kg)
      real :: surface_ssh

      ! Net roughness length (m)
      real :: rough

      ! Photosynthetic rates for different PFTs, if they were at the top 
      ! of the canopy (umol/m2 leaf/s).  Used by mortality function.
      real, allocatable, dimension(:) :: A_o_max  ! open stomata
      real, allocatable, dimension(:) :: A_c_max  ! closed stomata

      ! This will hold the stomatal conductance data from the previous 
      ! time step corresponding to A_o_max
      type(stoma_data), allocatable, dimension(:) :: old_stoma_data_max

      ! Average daily temperature [K]
      real :: avg_daily_temp

      ! average of rh [umol/m2/s] over FRQSTATE
      real :: omean_rh

      ! average of net ecosystem productivity (NEP) [umol/m2/s] over FRQSTATE
      real :: omean_nep

      ! Mean water flux from the canopy air to the atmosphere (kg_H2O/m2/s),
      ! averaged over FRQSTATE
      real :: omean_wflux
      real :: omean_latflux
      real :: omean_hflux

      ! Mean runoff (kg_H2O/m2/s), averaged over FRQSTATE
      real :: omean_runoff
      real :: omean_qrunoff

      ! Daily average of A_decomp, the temperature and moisture dependence
      ! of heterotrophic respiration.
      real :: dmean_A_decomp

      ! Daily average of the product A_decomp * f_decomp, which incorporates
      ! temperature, moisture, and N dependence of decomposition.
      real :: dmean_Af_decomp

      ! number of cohorts in the patch
      integer :: cohort_count

      ! Carbon available to establish recruits [kgC/m2]
      real, dimension(n_pft) :: repro

      ! Vegetation roughness length (m)
      real :: veg_rough

      ! Vegetation height (m)
      real :: veg_height

      ! Leaf area index (m2 leaf / m2 ground)
      real :: lai

      ! Input to fast soil carbon pool [kgC/m2/frq_phenology]
      real :: fsc_in

      ! Input to structural soil carbon pool [kgC/m2/frq_phenology]
      real :: ssc_in

      ! Input to soil lignin pool [kgC/m2/frq_phenology]
      real :: ssl_in

      ! Input to fast soil nitrogen pool [kgN/m2/frq_phenology]
      real :: fsn_in

      ! Plant nitrogen update summed over all cohorts [kgN/m2/frq_phenology]
      real :: total_plant_nitrogen_uptake

      ! Short wave radiation absorbed by the ground (W/m2)
      real :: rshort_g

      ! Short wave radiation absorbed by the ground, beam component (W/m2)
      real :: rshort_g_beam

      ! Short wave radiation absorbed by the ground, diffuse component (W/m2)
      real :: rshort_g_diffuse

      ! Long wave radiation absorbed by the ground (W/m2)
      real :: rlong_g

      ! Long wave radiation absorbed by the ground (W/m2), due to the 
      ! surface and vegetation alone
      real :: rlong_g_surf

      ! Long wave radiation absorbed by the ground (W/m2), due to the 
      ! incident long wave alone
      real :: rlong_g_incid

      ! Long wave radiation absorbed by the surface water (W/m2)
      real :: rlong_s

      ! Long wave radiation absorbed by the surface water (W/m2), due to 
      ! the surface and vegetation alone
      real :: rlong_s_surf

      ! Long wave radiation absorbed by the surface water (W/m2), due 
      ! to the incident atmospheric long wave alone
      real :: rlong_s_incid

      ! Patch albedo
      real :: albedo

      ! Patch albedo, beam component
      real :: albedo_beam

      ! Patch albedo, diffuse component
      real :: albedo_diffuse

      ! Upward long wave radiation at the top of the canopy (W/m2)
      real :: rlongup

      ! Albedo for long wave radiation
      real :: rlong_albedo

      ! Total snow depth as calculated in the radiation scheme.  Used for 
      ! checking if cohorts are buried.
      real :: total_snow_depth

      ! Fraction of vegetation covered with snow.  Used for computing 
      ! surface roughness.
      real :: snowfac

      ! can-to-atm heat xfer this step [kg_air K/m^2]
      real :: sxfer_t

      ! can-to-atm vapor xfer this step [kg_vap/m^2]
      real :: sxfer_r

      ! friction velocity [m/s]
      real :: ustar

      ! limitation of heterotrophic respiration due to physical 
      ! environmental factors (0-1 coefficient)
      real :: A_decomp

      ! damping of decomposition due to nitrogen immobilization 
      ! (0-1 coefficient)
      real :: f_decomp

      ! total heterotrophic respiration (umol/m2/s)
      real :: rh

      ! coarse woody debris contribution to rh (umol/m2/s)
      real :: cwd_rh

      ! Integer flag specifying whether this patch is to be fused
      integer :: fuse_flag

      ! Plant density broken down into size and PFT bins.  Used in patch fusion
      real, dimension(n_pft, ff_ndbh) :: pft_density_profile

      ! Above ground biomass in this patch [kgC/m2]
      real :: plant_ag_biomass

      ! Patch name (only used when restarting from ED1)
      character(len=64) :: pname

      ! Pointers
      type(cohort), pointer :: tallest  ! tallest cohort in the patch
      type(cohort), pointer :: shortest ! shortest cohort in the patch
      type(patch), pointer :: older     ! next-oldest patch
      type(patch), pointer :: younger   ! next-youngest patch
      type(site), pointer :: siteptr    ! pointer to the site

   End Type

   Type cohort
      
      ! (I) INITIALIZING A COHORT
      !     When you add a variable to this list, you need to make sure that
      !     it is initialized properly.  
      !    (A) Model restarts
      !       (1) Upon initializing a model run, the restart data are read
      !           in bare_ground_init().
      !       (2) State variables depending on the atmospheric state are
      !           set in ed_init_atm().
      !       (3) Other fundamental variables are set in init_ed_cohort_vars()
      !    (B) Recruitment: Basic state variables are set in 
      !        reproduction(); others set by a call to init_ed_cohort_vars().
      !    (C) Survivors from disturbance:  All area-based quantities
      !        need to be re-weighted.
      !    (D) Plantation: similar to reproduction.
      ! (II) DAILY AVERAGES USED IN VEGETATION_DYNAMICS().
      !    (A) These must be initialized according to (I) above.
      !    (B) These must be normalized before use.  This is done in 
      !        normalize_ed_daily_vars().
      !    (C) These must be reinitialized after use.  This is done 
      !        in reinitialize_ed_daily_vars().
      ! (III) OUTPUT QUANTITIES
      !    (A) These must be initialized according to (I) above.
      !    (B) Before write-out, they must be normalized.
      !       (1) For variables written on FRQSTATE, this is done 
      !           in normalized_ed_output_vars().
      !    (C) After write-out, they must be reset to zero.
      !       (1) For variables written on FRQSTATE, this is done 
      !           in zero_ed_output_vars().
      ! (IV) FUSE PATCH
      !      When you fuse patches, in fuse_2_patches(), the densities
      !      need to be re-weighted.
      ! (IV) FUSE COHORTS
      ! (V) SPLITTING COHORTS
      ! (VI) COPYING COHORTS
      !============================================================

      !================================================================
      ! Plant functional type:
      ! 1 - C4 grass
      ! 2 - tropical early successional tree
      ! 3 - tropical mid successional tree
      ! 4 - tropical late successional tree
      ! 5 - C3 grass
      ! 6 - northern pines (temperate)
      ! 7 - southern pines (temperate)
      ! 8 - late successional conifers
      ! 9 - early successional cold-deciduous hardwood
      ! 10 - mid successional cold-deciduous hardwood
      ! 11 - late successional cold-deciduous hardwood
      integer :: pft

      ! Density of plants (number/m2)
      real :: nplant

      ! Plant height (m)
      real :: hite

      ! Plant diameter at breast height (cm)
      real :: dbh

      ! Biomass of the structural stem (kgC/plant)
      real :: bdead

      ! Biomass of leaves (kgC/plant)
      real :: bleaf

      ! phenology_status codes:
      ! 0 - plant has the maximum LAI, given its size
      ! 1 - plant has an LAI between 0 and its maximum
      ! 2 - plant has no leaves
      integer :: phenology_status  

      ! Biomass of live tissue (kgC/plant)
      real :: balive

      ! Leaf area index (m2 leaf / m2 ground)
      real :: lai

      ! Plant storage pool of carbon [kgC/plant]
      real :: bstorage

      ! Heat capacity of vegetation (J/m2/K): 
      ! NOTE: AFTER DISTURBANCE, THIS IS NOT CURRENTLY
      ! RE-WEIGHTED BY SURVIVORSHIP AND DISTURBED AREA.  
      ! IT SHOULD BE WHEN I START USING REALISTIC HEAT CAPACITIES.
      ! IT SHOULD ALSO BE DONE IN FUSE_PATCHES.  ITS CALCULATION SHOULD
      ! ALSO BE MODIFIED IN COHORT FUSION.
      real :: hcapveg

      ! Vegetation temperature (K)
      real :: veg_temp

      ! Vegetation surface water (kg/m2)
      real :: veg_water

      ! Gross primary productivity (GPP) [umol/m2/s], averaged over the 
      ! output frequency (FRQSTATE)
      real :: omean_gpp

      ! Mean leaf respiration rate (umol/m2 ground/s), averaged over FRQSTATE
      real :: omean_leaf_resp

      ! Mean root respiration rate (umol/m2 ground/s), averaged over FRQSTATE
      real :: omean_root_resp

      ! Mean leaf respiration rate (umol/m2 ground/s), averaged over 
      ! frq_phenology (phenology update time; of the order of 1 day)
      real :: dmean_leaf_resp

      ! Mean root respiration rate (umol/m2 ground/s), averaged over 
      ! frq_phenology (phenology update time; of the order of 1 day)
      real :: dmean_root_resp

      ! Gross primary productivity (GPP) [umol/m2 ground/s], averaged over 
      ! frq_phenology (phenology update time; of the order of 1 day)
      real :: dmean_gpp

      ! Potential GPP in the absence of N limitation [umol/m2 ground/s], 
      ! averaged over frq_phenology (phenology update time; of the order 
      ! of 1 day)
      real :: dmean_gpp_pot

      ! Maximum GPP if cohort were at the top of the canopy 
      ! [umol/m2 ground/s], averaged over frq_phenology (phenology update 
      ! time; of the order of 1 day)
      real :: dmean_gpp_max

      ! Plant growth respiration (kgC/plant/frq_phenology)
      real :: growth_respiration

      ! Plant storage respiration (kgC/plant/frq_phenology)
      real :: storage_respiration

      ! Plant virtual leaf respiration (kgC/plant/frq_phenology)
      real :: vleaf_respiration

      ! Weighting factor for open and closed stomata due to N limitation
      real :: fsn

      ! Plant mortality rate [plants/m2/month]
      real :: monthly_dndt

      ! This is where you keep the derivatives of the 
      ! stomatal conductance residual and the old met conditions.
      type(stoma_data), pointer :: old_stoma_data

      ! Transpiration rate, open stomata (mm/s)
      real :: Psi_open

      ! Monthly carbon balance for past 12 months and the current month 
      ! (kgC/plant)
      real, dimension(13) :: cb

      ! Maximum monthly carbon balance for past 12 months and the current 
      ! month  if cohort were at the top of the canopy (kgC/plant)
      real, dimension(13) :: cb_max

      ! Annual average ratio of cb/cb_max
      real :: cbr_bar

      ! This specifies the index of the deepest soil layer of which the 
      ! cohort can access water.
      integer :: krdepth

      ! The model reports annual growth, mortality, cut and recruitment.  
      ! These rates are calculated with respect to two censuses.  This 
      ! is the flag specifying if a cohort was present at the first
      ! census.
      integer :: first_census

      ! Flag specifying if this cohort is a new recruit with respect
      ! to the first census.
      integer :: new_recruit_flag

      ! Photosynthetically active radiation (PAR) absorbed by the 
      ! cohort (units are Einsteins/m2/s)
      real :: par_v

      ! Photosynthetically active radiation (PAR) absorbed by the 
      ! cohort (units are Einsteins/m2/s), beam component
      real :: par_v_beam

      ! Photosynthetically active radiation (PAR) absorbed by the 
      ! cohort (units are Einsteins/m2/s), diffuse component
      real :: par_v_diffuse

      ! Total short wave radiation absorbed by the cohort, W/m2
      real :: rshort_v

      ! Total short wave radiation absorbed by the cohort, W/m2, beam component
      real :: rshort_v_beam

      ! Total short wave radiation absorbed by the cohort, W/m2, diffuse 
      ! component
      real :: rshort_v_diffuse

      ! Total long wave radiation absorbed by the cohort (W/m2)
      real :: rlong_v

      ! Total long wave radiation absorbed by the cohort (W/m2), due to 
      ! the temperature of the vegetation and surface alone
      real :: rlong_v_surf

      ! Total long wave radiation absorbed by the cohort (W/m2), due to 
      ! the incident atmospheric long wave alone
      real :: rlong_v_incid

      ! Leaf aerodynamic resistance (s/m)
      real :: rb

      ! Photosynthesis rate, open stomata (umol/m2 leaf/s)
      real :: A_open

      ! Photosynthesis rate, closed stomata (umol/m2 leaf/s)
      real :: A_closed

      ! Transpiration rate, closed stomata (mm/s)
      real :: Psi_closed

      ! Stomatal resistance for water, open stomata (s/m)
      real :: rsw_open

      ! Stomatal resistance for water, closed stomata (s/m)
      real :: rsw_closed

      ! Weighting factor for open and closed stomata (fsw=1 => fully open)
      real :: fsw

      ! Product of fsw and fsn
      real :: fs_open
      
      ! Net stomatal resistance [s/m]
      real :: stomatal_resistance

      ! Plant maintenance costs due to turnover of leaves and fine 
      ! roots [kgC/plant/frq_phenology]
      real :: maintenance_costs

      ! Amount of seeds produced for dispersal [kgC/plant]
      real :: bseeds

      ! POINTERS 
      type(cohort), pointer :: shorter ! next-shortest cohort
      type(cohort), pointer :: taller  ! next-tallest cohort
      type(patch), pointer :: patchptr ! pointer to its patch
      type(site), pointer :: siteptr   ! pointer to its site

   End Type

end Module ed_structure_defs
