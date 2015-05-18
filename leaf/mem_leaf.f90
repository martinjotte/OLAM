!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

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

!===============================================================================
Module mem_leaf

   use mem_mksfc, only: itab_wls_vars

   implicit none
   
   private :: itab_wls_vars

! LAND SURFACE GRID TABLES

   type (itab_wls_vars), target, allocatable :: itab_wl(:)

! DERIVED TYPES TO HOLD GLOBAL GRID INDICES FOR A PARALLEL RUN

   Type itabg_wl_vars
      integer :: iwl_myrank = -1
      integer :: irank      = -1
   End Type itabg_wl_vars

   type (itabg_wl_vars), allocatable, target :: itabg_wl(:)

! LAND SURFACE MODEL VARIABLES

   Type land_vars

! W-point grid variables stored on LANDFILE:

      real, allocatable :: area (:) ! cell surface area [m^2]
      real, allocatable :: glatw(:) ! latitude of land cell W points
      real, allocatable :: glonw(:) ! longitude of land cell W points
      real, allocatable :: xew  (:) ! earth x coord of land W points
      real, allocatable :: yew  (:) ! earth y coord of land W points
      real, allocatable :: zew  (:) ! earth z coord of land W points
      real, allocatable :: topw (:) ! topographic height of land W points
      real, allocatable :: wnx  (:) ! norm unit vector x comp of land cells
      real, allocatable :: wny  (:) ! norm unit vector y comp of land cells
      real, allocatable :: wnz  (:) ! norm unit vector z comp of land cells

      integer, allocatable :: leaf_class  (:)  ! leaf ("vegetation") class
      integer, allocatable :: ntext_soil(:,:)  ! soil textural class

! Atmospheric near-surface properties

      real, allocatable :: rhos   (:) ! air density [kg_air/m^3]
      real, allocatable :: vels   (:) ! wind speed [m/s]
      real, allocatable :: prss   (:) ! air pressure [Pa]
      real, allocatable :: airtemp(:) ! air temperature [K]
      real, allocatable :: airshv (:) ! air specific humidity[kg_vap/kg_air]

! Canopy to atmosphere turbulent flux quantities

      real, allocatable :: ustar (:) ! friction velocity [m/s]
      real, allocatable :: vkmsfc(:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: sfluxt(:)
      real, allocatable :: sfluxr(:)
      real, allocatable :: sfluxc(:)

      real, allocatable :: sxfer_t   (:) ! canair-to-atm heat xfer this step [kg_air K/m^2]
      real, allocatable :: sxfer_r   (:) ! canair-to-atm vapor xfer this step [kg_vap/m^2]
      real, allocatable :: sxfer_c   (:) ! canair-to-atm CO2 xfer this step [ppm/m^2]

      real, allocatable :: ggaer  (:) ! canopy-atmosphere aerodynamic conductance [m/s]
      real, allocatable :: ed_zeta(:)
      real, allocatable :: ed_rib (:)

! Radiative flux quantities

      real, allocatable :: albedo_beam   (:) ! net canopy/soil s/w beam albedo [0-1]
      real, allocatable :: albedo_diffuse(:) ! net canopy/soil s/w diffuse albedo [0-1]
      real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
      real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
      real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
      real, allocatable :: rlong_albedo  (:) ! net canopy/soil l/w albedo [0-1]
      real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]

      real, allocatable :: rshort_g  (:) ! s/w net rad flux to soil [W/m^2]
      real, allocatable :: rshort_s(:,:) ! s/w net rad flux to sfc water [W/m^2]
      real, allocatable :: rshort_v  (:) ! s/w net rad flux to veg [W/m^2]
      real, allocatable :: rlong_g   (:) ! l/w net rad flux to soil [W/m^2]
      real, allocatable :: rlong_s   (:) ! l/w net rad flux to sfc water [W/m^2]
      real, allocatable :: rlong_v   (:) ! l/w net rad flux to veg [W/m^2]

      real, allocatable :: cosz      (:) ! cosine of the solar zenith angle of the land cell

! Precipitation fluxes

      real, allocatable :: pcpg (:) ! new pcp amount this leaf timestep [kg/m^2]
      real, allocatable :: qpcpg(:) ! new pcp energy this leaf timestep [J/m^2]
      real, allocatable :: dpcpg(:) ! new pcp depth this leaf timestep [m]

! Canopy and surface quantities:

      real, allocatable :: can_depth(:) ! canopy depth for heat & vap capacity [m]
      real, allocatable :: cantemp (:) ! canopy air temp [K]
      real, allocatable :: canshv  (:) ! canopy air vapor spec hum [kg_vap/kg_air]

      integer, allocatable :: nlev_sfcwater (:) ! # of active surface water levels
      real, allocatable :: soil_water     (:,:) ! soil water content [vol_water/vol_tot]
      real, allocatable :: soil_energy    (:,:) ! soil energy [J/m^3]
      real, allocatable :: head0            (:) ! LBC total hydraulic head [m]
      real, allocatable :: head1            (:) ! UBC total hydraulic head [m]
      real, allocatable :: sfcwater_mass  (:,:) ! surface water mass [kg/m^2]
      real, allocatable :: sfcwater_energy(:,:) ! surface water energy [J/kg]
      real, allocatable :: sfcwater_depth (:,:) ! surface water depth [m]

      real, allocatable :: hcapveg     (:) ! veg heat capacity [J/(m^2 K)]
      real, allocatable :: surface_ssh (:) ! surface saturation spec hum [kg_vap/kg_air]
      real, allocatable :: ground_shv  (:) ! soil vapor spec hum [kg_vap/kg_air]
      real, allocatable :: rough       (:) ! land cell roughness height [m]
      real, allocatable :: veg_fracarea(:) ! veg fractional area
      real, allocatable :: veg_lai     (:) ! veg leaf area index
      real, allocatable :: veg_rough   (:) ! veg roughness height [m]
      real, allocatable :: veg_height  (:) ! veg height [m]
      real, allocatable :: veg_albedo  (:) ! veg albedo
      real, allocatable :: veg_tai     (:) ! veg total area index
      real, allocatable :: veg_water   (:) ! veg sfc water content [kg/m^2]
      real, allocatable :: veg_temp    (:) ! veg temp [K]
      real, allocatable :: veg_ndvip   (:) ! veg past ndvi (obs time)
      real, allocatable :: veg_ndvif   (:) ! veg future ndvi (obs time)
      real, allocatable :: veg_ndvic   (:) ! veg current ndvi
      real, allocatable :: stom_resist (:) ! veg stomatal resistance [s/m]
      real, allocatable :: snowfac     (:) ! frac veg burial by snowcover
      real, allocatable :: vf          (:) ! frac coverage of non-buried part of veg

      ! ED model variables:

      integer, allocatable :: ed_flag (:) ! if ED is being run in this cell (0 or 1)
      integer, allocatable :: ed_ifm  (:) ! if ED is being run in this cell (0 or 1)
      integer, allocatable :: ed_ipy  (:) ! if ED is being run in this cell (0 or 1)
      real, allocatable :: gpp        (:) ! Gross primary productivity (umol/m2/s)
      real, allocatable :: rh         (:) ! Heterotrophic respiration (umol/m2/s)
      real, allocatable :: nep        (:) ! Net ecosystem productivity (umol/m2/s)
      real, allocatable :: agb        (:) ! Site above ground biomass (kgC/m2)
      real, allocatable :: basal_area (:) ! Basal area (m2/ha)
      real, allocatable :: agb_growth (:) ! Net plant growth rate [kgC/m2/y]
      real, allocatable :: agb_mort   (:) ! Net plant mortality rate [kgC/m2/y]
      real, allocatable :: agb_cut    (:) ! Net plant cut rate [kgC/m2/y]
      real, allocatable :: agb_recruit(:) ! Net plant recruitment rate [kgC/m2/y]

   End Type land_vars

   type (land_vars), target :: land

Contains

!=========================================================================

   subroutine alloc_land_grid(mwl, nzg)

     use misc_coms, only: rinit, runtype

     implicit none

     integer, intent(in) :: mwl, nzg

!    Allocate land grid tables

     allocate (itab_wl(mwl))

!    Allocate and initialize land grid information from mksfc

     allocate (land%area (mwl)) ; land%area  = rinit
     allocate (land%glatw(mwl)) ; land%glatw = rinit
     allocate (land%glonw(mwl)) ; land%glonw = rinit
     allocate (land%xew  (mwl)) ; land%xew   = rinit
     allocate (land%yew  (mwl)) ; land%yew   = rinit
     allocate (land%zew  (mwl)) ; land%zew   = rinit
     allocate (land%topw (mwl)) ; land%topw  = rinit
     allocate (land%wnx  (mwl)) ; land%wnx   = rinit
     allocate (land%wny  (mwl)) ; land%wny   = rinit
     allocate (land%wnz  (mwl)) ; land%wnz   = rinit

     allocate (land%leaf_class    (mwl)) ; land%leaf_class = 0
     allocate (land%ntext_soil(nzg,mwl)) ; land%ntext_soil = 0

   end subroutine alloc_land_grid

!=========================================================================

   subroutine alloc_leaf(mwl, nzg, nzs)
     use misc_coms, only: rinit
     implicit none

     integer, intent(in) :: mwl, nzg, nzs

!    Allocate and initialize land arrays

     allocate (land%rhos               (mwl)) ; land%rhos            = rinit
     allocate (land%vels               (mwl)) ; land%vels            = rinit
     allocate (land%prss               (mwl)) ; land%prss            = rinit
     allocate (land%airtemp            (mwl)) ; land%airtemp         = rinit
     allocate (land%airshv             (mwl)) ; land%airshv          = rinit

     allocate (land%ustar              (mwl)) ; land%ustar           = rinit
     allocate (land%vkmsfc             (mwl)) ; land%vkmsfc          = rinit
     allocate (land%sfluxt             (mwl)) ; land%sfluxt          = rinit
     allocate (land%sfluxr             (mwl)) ; land%sfluxr          = rinit
     allocate (land%sfluxc             (mwl)) ; land%sfluxc          = rinit

     allocate (land%sxfer_t            (mwl)) ; land%sxfer_t         = 0.0
     allocate (land%sxfer_r            (mwl)) ; land%sxfer_r         = 0.0
     allocate (land%sxfer_c            (mwl)) ; land%sxfer_c         = 0.0

     allocate (land%ggaer              (mwl)) ; land%ggaer           = 0.0
     allocate (land%ed_zeta            (mwl)) ; land%ed_zeta         = 0.0
     allocate (land%ed_rib             (mwl)) ; land%ed_rib          = 0.0

     allocate (land%albedo_beam        (mwl)) ; land%albedo_beam     = rinit
     allocate (land%albedo_diffuse     (mwl)) ; land%albedo_diffuse  = rinit
     allocate (land%rshort             (mwl)) ; land%rshort          = rinit
     allocate (land%rshort_diffuse     (mwl)) ; land%rshort_diffuse  = rinit
     allocate (land%rlong              (mwl)) ; land%rlong           = rinit
     allocate (land%rlong_albedo       (mwl)) ; land%rlong_albedo    = rinit
     allocate (land%rlongup            (mwl)) ; land%rlongup         = rinit
     allocate (land%rshort_g           (mwl)) ; land%rshort_g        = rinit
     allocate (land%rshort_s       (nzs,mwl)) ; land%rshort_s        = rinit
     allocate (land%rshort_v           (mwl)) ; land%rshort_v        = rinit
     allocate (land%rlong_g            (mwl)) ; land%rlong_g         = rinit
     allocate (land%rlong_s            (mwl)) ; land%rlong_s         = rinit
     allocate (land%rlong_v            (mwl)) ; land%rlong_v         = rinit

     allocate (land%cosz               (mwl)) ; land%cosz            = rinit

     allocate (land%pcpg               (mwl)) ; land%pcpg            = 0.0
     allocate (land%qpcpg              (mwl)) ; land%qpcpg           = 0.0
     allocate (land%dpcpg              (mwl)) ; land%dpcpg           = 0.0

     allocate (land%can_depth          (mwl)) ; land%can_depth       = rinit
     allocate (land%cantemp            (mwl)) ; land%cantemp         = rinit
     allocate (land%canshv             (mwl)) ; land%canshv          = rinit

     allocate (land%nlev_sfcwater      (mwl)) ; land%nlev_sfcwater   = 0
     allocate (land%soil_water     (nzg,mwl)) ; land%soil_water      = rinit
     allocate (land%soil_energy    (nzg,mwl)) ; land%soil_energy     = rinit
     allocate (land%head0              (mwl)) ; land%head0           = rinit
     allocate (land%head1              (mwl)) ; land%head1           = rinit
     allocate (land%sfcwater_mass  (nzs,mwl)) ; land%sfcwater_mass   = rinit
     allocate (land%sfcwater_energy(nzs,mwl)) ; land%sfcwater_energy = rinit
     allocate (land%sfcwater_depth (nzs,mwl)) ; land%sfcwater_depth  = rinit

     allocate (land%hcapveg            (mwl)) ; land%hcapveg         = rinit
     allocate (land%surface_ssh        (mwl)) ; land%surface_ssh     = rinit
     allocate (land%ground_shv         (mwl)) ; land%ground_shv      = rinit
     allocate (land%rough              (mwl)) ; land%rough           = rinit
     allocate (land%veg_fracarea       (mwl)) ; land%veg_fracarea    = rinit
     allocate (land%veg_lai            (mwl)) ; land%veg_lai         = rinit
     allocate (land%veg_rough          (mwl)) ; land%veg_rough       = rinit
     allocate (land%veg_height         (mwl)) ; land%veg_height      = rinit
     allocate (land%veg_albedo         (mwl)) ; land%veg_albedo      = rinit
     allocate (land%veg_tai            (mwl)) ; land%veg_tai         = rinit
     allocate (land%veg_water          (mwl)) ; land%veg_water       = rinit
     allocate (land%veg_temp           (mwl)) ; land%veg_temp        = rinit
     allocate (land%veg_ndvip          (mwl)) ; land%veg_ndvip       = rinit
     allocate (land%veg_ndvif          (mwl)) ; land%veg_ndvif       = rinit
     allocate (land%veg_ndvic          (mwl)) ; land%veg_ndvic       = rinit
     allocate (land%stom_resist        (mwl)) ; land%stom_resist     = rinit
     allocate (land%snowfac            (mwl)) ; land%snowfac         = rinit
     allocate (land%vf                 (mwl)) ; land%vf              = rinit

     allocate (land%ed_flag            (mwl)) ; land%ed_flag         = 0
     allocate (land%ed_ifm             (mwl)) ; land%ed_ifm          = 0
     allocate (land%ed_ipy             (mwl)) ; land%ed_ipy          = 0
     allocate (land%gpp                (mwl)) ; land%gpp             = rinit
     allocate (land%rh                 (mwl)) ; land%rh              = rinit
     allocate (land%nep                (mwl)) ; land%nep             = rinit
     allocate (land%agb                (mwl)) ; land%agb             = rinit
     allocate (land%basal_area         (mwl)) ; land%basal_area      = rinit
     allocate (land%agb_growth         (mwl)) ; land%agb_growth      = rinit
     allocate (land%agb_mort           (mwl)) ; land%agb_mort        = rinit
     allocate (land%agb_cut            (mwl)) ; land%agb_cut         = rinit
     allocate (land%agb_recruit        (mwl)) ; land%agb_recruit     = rinit

   end subroutine alloc_leaf

!=========================================================================

   subroutine filltab_leaf()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(land%rhos)) call increment_vtable('LAND%RHOS', 'LW', rvar1=land%rhos)

     if (allocated(land%vels)) call increment_vtable('LAND%VELS', 'LW', rvar1=land%vels)

     if (allocated(land%prss)) call increment_vtable('LAND%PRSS', 'LW', rvar1=land%prss)

     if (allocated(land%airtemp)) call increment_vtable('LAND%AIRTEMP', 'LW', rvar1=land%airtemp)

     if (allocated(land%airshv)) call increment_vtable('LAND%AIRSHV', 'LW', rvar1=land%airshv)

     if (allocated(land%ustar)) call increment_vtable('LAND%USTAR', 'LW', rvar1=land%ustar)

     if (allocated(land%vkmsfc)) call increment_vtable('LAND%VKMSFC', 'LW', rvar1=land%vkmsfc)

     if (allocated(land%sfluxt)) call increment_vtable('LAND%SFLUXT', 'LW', rvar1=land%sfluxt)

     if (allocated(land%sfluxr)) call increment_vtable('LAND%SFLUXR', 'LW', rvar1=land%sfluxr)

     if (allocated(land%sfluxc)) call increment_vtable('LAND%SFLUXC', 'LW', rvar1=land%sfluxc)

     if (allocated(land%sxfer_t)) call increment_vtable('LAND%SXFER_T', 'LW', rvar1=land%sxfer_t)

     if (allocated(land%sxfer_r)) call increment_vtable('LAND%SXFER_R', 'LW', rvar1=land%sxfer_r)

     if (allocated(land%sxfer_c)) call increment_vtable('LAND%SXFER_C', 'LW', rvar1=land%sxfer_c)

     if (allocated(land%ggaer)) call increment_vtable('LAND%ED_GGAER', 'LW', rvar1=land%ggaer)

     if (allocated(land%ed_zeta)) call increment_vtable('LAND%ED_ZETA', 'LW', rvar1=land%ed_zeta)

     if (allocated(land%ed_rib)) call increment_vtable('LAND%ED_RIB', 'LW', rvar1=land%ed_rib)

     if (allocated(land%albedo_beam)) call increment_vtable('LAND%ALBEDO_BEAM', 'LW', rvar1=land%albedo_beam)

     if (allocated(land%albedo_diffuse)) call increment_vtable('LAND%ALBEDO_DIFFUSE', 'LW', rvar1=land%albedo_diffuse)

     if (allocated(land%rshort)) call increment_vtable('LAND%RSHORT', 'LW', rvar1=land%rshort)

     if (allocated(land%rshort_diffuse)) call increment_vtable('LAND%RSHORT_DIFFUSE', 'LW', rvar1=land%rshort_diffuse)

     if (allocated(land%rlong)) call increment_vtable('LAND%RLONG', 'LW', rvar1=land%rlong)

     if (allocated(land%rlong_albedo)) call increment_vtable('LAND%RLONG_ALBEDO', 'LW', rvar1=land%rlong_albedo)

     if (allocated(land%rlongup)) call increment_vtable('LAND%RLONGUP', 'LW', rvar1=land%rlongup)

     if (allocated(land%rshort_g)) call increment_vtable('LAND%RSHORT_G', 'LW', rvar1=land%rshort_g)

     if (allocated(land%rshort_s)) call increment_vtable('LAND%RSHORT_S', 'LW', rvar2=land%rshort_s)

     if (allocated(land%rshort_v)) call increment_vtable('LAND%RSHORT_V', 'LW', rvar1=land%rshort_v)

     if (allocated(land%rlong_g)) call increment_vtable('LAND%RLONG_G', 'LW', rvar1=land%rlong_g)

     if (allocated(land%rlong_s)) call increment_vtable('LAND%RLONG_S', 'LW', rvar1=land%rlong_s)

     if (allocated(land%rlong_v)) call increment_vtable('LAND%RLONG_V', 'LW', rvar1=land%rlong_v)

     if (allocated(land%cosz)) call increment_vtable('LAND%COSZ', 'LW', rvar1=land%cosz)

     if (allocated(land%pcpg)) call increment_vtable('LAND%PCPG', 'LW', rvar1=land%pcpg)

     if (allocated(land%qpcpg)) call increment_vtable('LAND%QPCPG', 'LW', rvar1=land%qpcpg)

     if (allocated(land%dpcpg)) call increment_vtable('LAND%DPCPG', 'LW', rvar1=land%dpcpg)

     if (allocated(land%can_depth)) call increment_vtable('LAND%CAN_DEPTH', 'LW', rvar1=land%can_depth)

     if (allocated(land%cantemp)) call increment_vtable('LAND%CANTEMP', 'LW', rvar1=land%cantemp)

     if (allocated(land%canshv)) call increment_vtable('LAND%CANSHV', 'LW', rvar1=land%canshv)

     if (allocated(land%nlev_sfcwater)) call increment_vtable('LAND%NLEV_SFCWATER', 'LW', ivar1=land%nlev_sfcwater)

     if (allocated(land%soil_water)) call increment_vtable('LAND%SOIL_WATER', 'LW', rvar2=land%soil_water)

     if (allocated(land%soil_energy)) call increment_vtable('LAND%SOIL_ENERGY', 'LW', rvar2=land%soil_energy)
     
     if (allocated(land%head0)) call increment_vtable('LAND%HEAD0', 'LW', rvar1=land%head0)

     if (allocated(land%head1)) call increment_vtable('LAND%HEAD1', 'LW', rvar1=land%head1)

     if (allocated(land%sfcwater_mass)) call increment_vtable('LAND%SFCWATER_MASS', 'LW', rvar2=land%sfcwater_mass)

     if (allocated(land%sfcwater_energy)) call increment_vtable('LAND%SFCWATER_ENERGY', 'LW', rvar2=land%sfcwater_energy)

     if (allocated(land%sfcwater_depth)) call increment_vtable('LAND%SFCWATER_DEPTH', 'LW', rvar2=land%sfcwater_depth)

     if (allocated(land%hcapveg)) call increment_vtable('LAND%HCAPVEG', 'LW', rvar1=land%hcapveg)

     if (allocated(land%surface_ssh)) call increment_vtable('LAND%SURFACE_SSH', 'LW', rvar1=land%surface_ssh)

     if (allocated(land%ground_shv)) call increment_vtable('LAND%GROUND_SHV', 'LW', rvar1=land%ground_shv)

     if (allocated(land%rough)) call increment_vtable('LAND%ROUGH', 'LW', rvar1=land%rough)

     if (allocated(land%veg_fracarea)) call increment_vtable('LAND%VEG_FRACAREA', 'LW', rvar1=land%veg_fracarea)

     if (allocated(land%veg_lai)) call increment_vtable('LAND%VEG_LAI', 'LW', rvar1=land%veg_lai)

     if (allocated(land%veg_rough)) call increment_vtable('LAND%VEG_ROUGH', 'LW', rvar1=land%veg_rough)

     if (allocated(land%veg_height)) call increment_vtable('LAND%VEG_HEIGHT', 'LW', rvar1=land%veg_height)

     if (allocated(land%veg_albedo)) call increment_vtable('LAND%VEG_ALBEDO', 'LW', rvar1=land%veg_albedo)

     if (allocated(land%veg_tai)) call increment_vtable('LAND%VEG_TAI', 'LW', rvar1=land%veg_tai)

     if (allocated(land%veg_water)) call increment_vtable('LAND%VEG_WATER', 'LW', rvar1=land%veg_water)

     if (allocated(land%veg_temp)) call increment_vtable('LAND%VEG_TEMP', 'LW', rvar1=land%veg_temp)

     if (allocated(land%veg_ndvic)) call increment_vtable('LAND%VEG_NDVIC', 'LW', rvar1=land%veg_ndvic)

     if (allocated(land%stom_resist)) call increment_vtable('LAND%STOM_RESIST', 'LW', rvar1=land%stom_resist)

     if (allocated(land%snowfac)) call increment_vtable('LAND%SNOWFAC', 'LW', rvar1=land%snowfac)

     if (allocated(land%vf)) call increment_vtable('LAND%VF', 'LW', rvar1=land%vf)

   end subroutine filltab_leaf

!=========================================================================

   subroutine filltab_ED()
     use var_tables, only: increment_EDtab, num_ED, vtab_ED
     implicit none

     if (allocated(land%gpp)) call increment_EDtab('LAND%GPP', hist=.true., mavg=.true., yavg=.true., rvar1=land%gpp)

     if (allocated(land%rh)) call increment_EDtab('LAND%RH', hist=.true., mavg=.true., yavg=.true., rvar1=land%rh)

     if (allocated(land%nep)) call increment_EDtab('LAND%NEP', hist=.true., mavg=.true., yavg=.true., rvar1=land%nep)

     if (allocated(land%agb)) call increment_EDtab('LAND%AGB', yavg=.true., rvar1=land%agb)

     if (allocated(land%basal_area)) call increment_EDtab('LAND%BASAL_AREA', yavg=.true., rvar1=land%basal_area)

     if (allocated(land%agb_growth)) call increment_EDtab('LAND%AGB_GROWTH', yavg=.true., rvar1=land%agb_growth)

     if (allocated(land%agb_mort)) call increment_EDtab('LAND%AGB_MORT', yavg=.true., rvar1=land%agb_mort)

     if (allocated(land%agb_cut)) call increment_EDtab('LAND%AGB_CUT', yavg=.true., rvar1=land%agb_cut)

     if (allocated(land%agb_recruit)) call increment_EDtab('LAND%AGB_RECRUIT', yavg=.true., rvar1=land%agb_recruit)
   end subroutine filltab_ED

End Module mem_leaf
