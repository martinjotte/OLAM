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
Module mem_leaf

   use ed_structure_defs, only: site
   use max_dims,          only: maxremote
   use mem_mksfc,         only: itab_mls_vars, itab_uls_vars, itab_wls_vars
   implicit none
   
   private :: site, maxremote, itab_mls_vars, itab_uls_vars, itab_wls_vars

! LAND SURFACE GRID TABLES

   type (itab_mls_vars), target, allocatable :: itab_ml(:)
   type (itab_uls_vars), target, allocatable :: itab_ul(:)
   type (itab_wls_vars), target, allocatable :: itab_wl(:)

! DERIVED TYPES TO HOLD GLOBAL GRID INDICES FOR A PARALLEL RUN

   Type itabg_ml_vars
      integer :: iml_myrank = -1
      integer :: irank      = -1
   End Type itabg_ml_vars

   Type itabg_ul_vars
      integer :: iul_myrank = -1
      integer :: irank      = -1
   End Type itabg_ul_vars

   Type itabg_wl_vars
      integer :: iwl_myrank = -1
      integer :: irank      = -1
   End Type itabg_wl_vars

   type (itabg_ml_vars), allocatable         :: itabg_ml(:)
   type (itabg_ul_vars), allocatable         :: itabg_ul(:)
   type (itabg_wl_vars), allocatable, target :: itabg_wl(:)

! DERIVED TYPE TO HOLD MPI TABLES FOR A PARALLEL RUN

   Type jtab_wl_mpi_vars
      integer, allocatable :: iwl(:)
      integer, allocatable :: jend(:)
   End Type jtab_wl_mpi_vars

   type (jtab_wl_mpi_vars) :: jtab_wl_mpi(maxremote) 

! LAND SURFACE MODEL VARIABLES

   Type land_vars

      ! Surface and soil types (currently read in from mksfc):

      integer, allocatable :: leaf_class   (:)  ! leaf ("vegetation") class
      integer, allocatable :: ntext_soil (:,:)  ! soil textural class

      ! M-point grid variables read in from mksfc:

      real, allocatable :: xem         (:) ! earth x coord of land M points
      real, allocatable :: yem         (:) ! earth y coord of land M points
      real, allocatable :: zem         (:) ! earth z coord of land M points
      real, allocatable :: zm          (:) ! topography height at land M points
      real, allocatable :: glatm       (:) ! latitude of land cell M points
      real, allocatable :: glonm       (:) ! longitude of land cell M points

      ! W-point grid variables read in from mksfc:

      real, allocatable :: area        (:) ! land cell area [m^2]
      real, allocatable :: xew         (:) ! earth x coord of land W points
      real, allocatable :: yew         (:) ! earth y coord of land W points
      real, allocatable :: zew         (:) ! earth z coord of land W points
      real, allocatable :: glatw       (:) ! latitude of land cell W points
      real, allocatable :: glonw       (:) ! longitude of land cell W points
      real, allocatable :: wnx         (:) ! norm unit vector x comp of land cells
      real, allocatable :: wny         (:) ! norm unit vector y comp of land cells
      real, allocatable :: wnz         (:) ! norm unit vector z comp of land cells

      ! LEAF3 model variables:

      integer, allocatable :: nlev_sfcwater (:) ! # of active surface water levels
      real, allocatable :: soil_water     (:,:) ! soil water content [vol_water/vol_tot]
      real, allocatable :: soil_energy    (:,:) ! soil energy [J/m^3]
      real, allocatable :: head0            (:) ! LBC total hydraulic head [m]
      real, allocatable :: head1            (:) ! UBC total hydraulic head [m]
      real, allocatable :: sfcwater_mass  (:,:) ! surface water mass [kg/m^2]
      real, allocatable :: sfcwater_energy(:,:) ! surface water energy [J/kg]
      real, allocatable :: sfcwater_depth (:,:) ! surface water depth [m]
      real, allocatable :: rshort_s       (:,:) ! s/w net rad flux to sfc water [W/m^2]

      real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
      real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
      real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
      real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: rlong_albedo  (:) ! net canopy/soil l/w albedo [0-1]
      real, allocatable :: albedo_beam   (:) ! net canopy/soil s/w beam albedo [0-1]
      real, allocatable :: albedo_diffuse(:) ! net canopy/soil s/w diffuse albedo [0-1]

      real, allocatable :: rshort_g    (:) ! s/w net rad flux to soil [W/m^2]
      real, allocatable :: rshort_v    (:) ! s/w net rad flux to veg [W/m^2]
      real, allocatable :: rlong_g     (:) ! l/w net rad flux to soil [W/m^2]
      real, allocatable :: rlong_s     (:) ! l/w net rad flux to sfc water [W/m^2]
      real, allocatable :: rlong_v     (:) ! l/w net rad flux to veg [W/m^2]
      real, allocatable :: rhos        (:) ! atmosphere air density [kg_air/m^3]
      real, allocatable :: vels        (:) ! surface wind speed [m/s]
      real, allocatable :: prss        (:) ! air pressure [Pa]
      real, allocatable :: ustar       (:) ! friction velocity [m/s]
      real, allocatable :: sxfer_t     (:) ! canair-to-atm heat xfer this step [kg_air K/m^2]
      real, allocatable :: sxfer_r     (:) ! canair-to-atm vapor xfer this step [kg_vap/m^2]
      real, allocatable :: sxfer_tsav  (:) ! saved previous value of sxfer_t
      real, allocatable :: sxfer_rsav  (:) ! saved previous value of sxter_r
      real, allocatable :: can_depth   (:) ! canopy depth for heat & vap capacity [m]
      real, allocatable :: hcapveg     (:) ! veg heat capacity [J/(m^2 K)]
      real, allocatable :: can_temp    (:) ! canopy air temp [K]
      real, allocatable :: can_shv     (:) ! canopy air vapor spec hum [kg_vap/kg_air]
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
      real, allocatable :: pcpg        (:) ! new pcp amount this leaf timestep [kg/m^2]
      real, allocatable :: qpcpg       (:) ! new pcp energy this leaf timestep [J/m^2]
      real, allocatable :: dpcpg       (:) ! new pcp depth this leaf timestep [m]
      real, allocatable :: cosz        (:) ! cosine of the solar zenith angle of the land cell
      integer, allocatable :: lsl      (:) ! Lowest soil layer for a given iland cell.

      ! ED model variables:

      integer, allocatable :: ed_flag  (:) ! if ED is being run in this cell (0 or 1)
      real, allocatable :: gpp         (:) ! Gross primary productivity (umol/m2/s)
      real, allocatable :: rh          (:) ! Heterotrophic respiration (umol/m2/s)
      real, allocatable :: nep         (:) ! Net ecosystem productivity (umol/m2/s)
      real, allocatable :: agb         (:) ! Site above ground biomass (kgC/m2)
      real, allocatable :: basal_area  (:) ! Basal area (m2/ha)
      real, allocatable :: agb_growth  (:) ! Net plant growth rate [kgC/m2/y]
      real, allocatable :: agb_mort    (:) ! Net plant mortality rate [kgC/m2/y]
      real, allocatable :: agb_cut     (:) ! Net plant cut rate [kgC/m2/y]
      real, allocatable :: agb_recruit (:) ! Net plant recruitment rate [kgC/m2/y]

   End Type land_vars

   type (land_vars), target :: land

! The basic unit of the ED model is the site, which loosely corresponds to a
! grid cell.  Now, the entire ED memory structure resides on linked lists.  
! This is a pointer to the first site of the linked list. 

   type (site), pointer :: first_site => null()

Contains

!=========================================================================

   subroutine alloc_land_grid(mml, mul, mwl, nzg)
     use misc_coms, only: rinit
     implicit none

     integer, intent(in) :: mml, mul, mwl, nzg

!    Allocate sea grid tables

     allocate (itab_ml(mml))
     allocate (itab_ul(mul))
     allocate (itab_wl(mwl))

!    Allocate and initialize land grid information from mksfc

     allocate (land%leaf_class    (mwl)) ; land%leaf_class = 0
     allocate (land%ntext_soil(nzg,mwl)) ; land%ntext_soil = 0

     allocate (land%area (mwl)) ; land%area  = rinit
     allocate (land%xew  (mwl)) ; land%xew   = rinit
     allocate (land%yew  (mwl)) ; land%yew   = rinit
     allocate (land%zew  (mwl)) ; land%zew   = rinit
     allocate (land%glatw(mwl)) ; land%glatw = rinit
     allocate (land%glonw(mwl)) ; land%glonw = rinit
     allocate (land%wnx  (mwl)) ; land%wnx   = rinit
     allocate (land%wny  (mwl)) ; land%wny   = rinit
     allocate (land%wnz  (mwl)) ; land%wnz   = rinit

     allocate (land%xem  (mml)) ; land%xem   = rinit
     allocate (land%yem  (mml)) ; land%yem   = rinit
     allocate (land%zem  (mml)) ; land%zem   = rinit
     allocate (land%zm   (mml)) ; land%zm    = rinit
     allocate (land%glatm(mml)) ; land%glatm = rinit
     allocate (land%glonm(mml)) ; land%glonm = rinit

   end subroutine alloc_land_grid

!=========================================================================

   subroutine alloc_leaf(mwl, nzg, nzs)
     use misc_coms, only: rinit
     implicit none

     integer, intent(in) :: mwl, nzg, nzs

!    Allocate and initialize land arrays

     allocate (land%nlev_sfcwater      (mwl)) ; land%nlev_sfcwater   = 0
     allocate (land%rhos               (mwl)) ; land%rhos            = rinit
     allocate (land%vels               (mwl)) ; land%vels            = rinit
     allocate (land%prss               (mwl)) ; land%prss            = rinit
     allocate (land%ustar              (mwl)) ; land%ustar           = rinit
     allocate (land%sxfer_t            (mwl)) ; land%sxfer_t         = rinit
     allocate (land%sxfer_r            (mwl)) ; land%sxfer_r         = rinit
     allocate (land%sxfer_tsav         (mwl)) ; land%sxfer_tsav      = rinit
     allocate (land%sxfer_rsav         (mwl)) ; land%sxfer_rsav      = rinit
     allocate (land%can_depth          (mwl)) ; land%can_depth       = rinit
     allocate (land%hcapveg            (mwl)) ; land%hcapveg         = rinit
     allocate (land%can_temp           (mwl)) ; land%can_temp        = rinit
     allocate (land%can_shv            (mwl)) ; land%can_shv         = rinit
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
     allocate (land%veg_ndvic          (mwl)) ; land%veg_ndvic       = rinit
     allocate (land%veg_ndvif          (mwl)) ; land%veg_ndvif       = rinit
     allocate (land%stom_resist        (mwl)) ; land%stom_resist     = rinit
     allocate (land%rshort             (mwl)) ; land%rshort          = rinit
     allocate (land%rshort_diffuse     (mwl)) ; land%rshort_diffuse  = rinit
     allocate (land%rlong              (mwl)) ; land%rlong           = rinit
     allocate (land%rlongup            (mwl)) ; land%rlongup         = rinit
     allocate (land%rlong_albedo       (mwl)) ; land%rlong_albedo    = rinit
     allocate (land%albedo_beam        (mwl)) ; land%albedo_beam     = rinit
     allocate (land%albedo_diffuse     (mwl)) ; land%albedo_diffuse  = rinit
     allocate (land%rshort_g           (mwl)) ; land%rshort_g        = rinit
     allocate (land%rshort_v           (mwl)) ; land%rshort_v        = rinit
     allocate (land%rlong_g            (mwl)) ; land%rlong_g         = rinit
     allocate (land%rlong_s            (mwl)) ; land%rlong_s         = rinit
     allocate (land%rlong_v            (mwl)) ; land%rlong_v         = rinit
     allocate (land%snowfac            (mwl)) ; land%snowfac         = rinit
     allocate (land%vf                 (mwl)) ; land%vf              = rinit
     allocate (land%pcpg               (mwl)) ; land%pcpg            = rinit
     allocate (land%qpcpg              (mwl)) ; land%qpcpg           = rinit
     allocate (land%dpcpg              (mwl)) ; land%dpcpg           = rinit
     allocate (land%ed_flag            (mwl)) ; land%ed_flag         = 0
     allocate (land%cosz               (mwl)) ; land%cosz            = rinit
     allocate (land%gpp                (mwl)) ; land%gpp             = rinit
     allocate (land%agb                (mwl)) ; land%agb             = rinit
     allocate (land%basal_area         (mwl)) ; land%basal_area      = rinit
     allocate (land%agb_growth         (mwl)) ; land%agb_growth      = rinit
     allocate (land%agb_mort           (mwl)) ; land%agb_mort        = rinit
     allocate (land%agb_cut            (mwl)) ; land%agb_cut         = rinit
     allocate (land%agb_recruit        (mwl)) ; land%agb_recruit     = rinit
     allocate (land%rh                 (mwl)) ; land%rh              = rinit
     allocate (land%nep                (mwl)) ; land%nep             = rinit
     allocate (land%soil_water     (nzg,mwl)) ; land%soil_water      = rinit
     allocate (land%soil_energy    (nzg,mwl)) ; land%soil_energy     = rinit
     allocate (land%head0              (mwl)) ; land%head0           = rinit
     allocate (land%head1              (mwl)) ; land%head1           = rinit
     allocate (land%lsl                (mwl)) ; land%lsl             = 0
     allocate (land%sfcwater_mass  (nzs,mwl)) ; land%sfcwater_mass   = rinit
     allocate (land%sfcwater_energy(nzs,mwl)) ; land%sfcwater_energy = rinit
     allocate (land%sfcwater_depth (nzs,mwl)) ; land%sfcwater_depth  = rinit
     allocate (land%rshort_s       (nzs,mwl)) ; land%rshort_s        = rinit

   end subroutine alloc_leaf

!=========================================================================

   subroutine filltab_leaf()

     use var_tables, only: vtab_r, num_var, increment_vtable
     use misc_coms,  only: iparallel, runtype

     implicit none

     if (iparallel==1 .or. runtype=='PLOTONLY' .or. runtype=='PARCOMBINE') then

        ! THESE ONLY NEED TO BE WRITTEN TO HISTORY FILE FOR PARALLEL RUNS,
        ! AND READ FOR PLOTONLY OR PARCOMBINE RUNS.

        if (allocated(itab_ml)) then

           call increment_vtable('IMGLOBE_L', 'LM', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ml%imglobe

        endif

        if (allocated(itab_ul)) then

           call increment_vtable('IUGLOBE_L', 'LU', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ul%iuglobe

           call increment_vtable('IRANKU_L',  'LU', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ul%irank

        endif

        if (allocated(itab_wl)) then

           call increment_vtable('IWGLOBE_L', 'LW', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_wl%iwglobe

           call increment_vtable('IRANKW_L',  'LW', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_wl%irank

        endif

     endif

     if (allocated(land%nlev_sfcwater)) then
        call increment_vtable('LAND%NLEV_SFCWATER', 'LW')
        vtab_r(num_var)%ivar1_p => land%nlev_sfcwater
     endif

     if (allocated(land%sxfer_t)) then
        call increment_vtable('LAND%SXFER_T', 'LW')
        vtab_r(num_var)%rvar1_p => land%sxfer_t
     endif

     if (allocated(land%sxfer_r)) then
        call increment_vtable('LAND%SXFER_R', 'LW')
        vtab_r(num_var)%rvar1_p => land%sxfer_r
     endif

     if (allocated(land%can_depth)) then
        call increment_vtable('LAND%CAN_DEPTH', 'LW')
        vtab_r(num_var)%rvar1_p => land%can_depth
     endif

     if (allocated(land%hcapveg)) then
        call increment_vtable('LAND%HCAPVEG', 'LW')
        vtab_r(num_var)%rvar1_p => land%hcapveg
     endif

     if (allocated(land%can_temp)) then
        call increment_vtable('LAND%CAN_TEMP', 'LW')
        vtab_r(num_var)%rvar1_p => land%can_temp
     endif

     if (allocated(land%can_shv)) then
        call increment_vtable('LAND%CAN_SHV', 'LW')
        vtab_r(num_var)%rvar1_p => land%can_shv
     endif

     if (allocated(land%surface_ssh)) then
        call increment_vtable('LAND%SURFACE_SSH', 'LW')
        vtab_r(num_var)%rvar1_p => land%surface_ssh
     endif

     if (allocated(land%ground_shv)) then
        call increment_vtable('LAND%GROUND_SHV', 'LW')
        vtab_r(num_var)%rvar1_p => land%ground_shv
     endif

     if (allocated(land%rough)) then
        call increment_vtable('LAND%ROUGH', 'LW')
        vtab_r(num_var)%rvar1_p => land%rough
     endif

     if (allocated(land%veg_fracarea)) then
        call increment_vtable('LAND%VEG_FRACAREA', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_fracarea
     endif

     if (allocated(land%veg_lai)) then
        call increment_vtable('LAND%VEG_LAI', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_lai
     endif

     if (allocated(land%veg_rough)) then
        call increment_vtable('LAND%VEG_ROUGH', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_rough
     endif

     if (allocated(land%veg_height)) then
        call increment_vtable('LAND%VEG_HEIGHT', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_height
     endif

     if (allocated(land%veg_albedo)) then
        call increment_vtable('LAND%VEG_ALBEDO', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_albedo
     endif

     if (allocated(land%veg_tai)) then
        call increment_vtable('LAND%VEG_TAI', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_tai
     endif

     if (allocated(land%veg_water)) then
        call increment_vtable('LAND%VEG_WATER', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_water
     endif

     if (allocated(land%veg_temp)) then
        call increment_vtable('LAND%VEG_TEMP', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_temp
     endif

     if (allocated(land%veg_ndvic)) then
        call increment_vtable('LAND%VEG_NDVIC', 'LW')
        vtab_r(num_var)%rvar1_p => land%veg_ndvic
     endif

     if (allocated(land%stom_resist)) then
        call increment_vtable('LAND%STOM_RESIST', 'LW')
        vtab_r(num_var)%rvar1_p => land%stom_resist
     endif

     if (allocated(land%rshort)) then
        call increment_vtable('LAND%RSHORT', 'LW')
        vtab_r(num_var)%rvar1_p => land%rshort
     endif

     if (allocated(land%rshort_diffuse)) then
        call increment_vtable('LAND%RSHORT_DIFFUSE', 'LW')
        vtab_r(num_var)%rvar1_p => land%rshort_diffuse
     endif

     if (allocated(land%rlong)) then
        call increment_vtable('LAND%RLONG', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlong
     endif

     if (allocated(land%rlongup)) then
        call increment_vtable('LAND%RLONGUP', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlongup
     endif

     if (allocated(land%rlong_albedo)) then
        call increment_vtable('LAND%RLONG_ALBEDO', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlong_albedo
     endif

     if (allocated(land%albedo_beam)) then
        call increment_vtable('LAND%ALBEDO_BEAM', 'LW')
        vtab_r(num_var)%rvar1_p => land%albedo_beam
     endif

     if (allocated(land%albedo_diffuse)) then
        call increment_vtable('LAND%ALBEDO_DIFFUSE', 'LW')
        vtab_r(num_var)%rvar1_p => land%albedo_diffuse
     endif

     if (allocated(land%rshort_g)) then
        call increment_vtable('LAND%RSHORT_G', 'LW')
        vtab_r(num_var)%rvar1_p => land%rshort_g
     endif

     if (allocated(land%rshort_v)) then
        call increment_vtable('LAND%RSHORT_V', 'LW')
        vtab_r(num_var)%rvar1_p => land%rshort_v
     endif

     if (allocated(land%rlong_g)) then
        call increment_vtable('LAND%RLONG_G', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlong_g
     endif

     if (allocated(land%rlong_s)) then
        call increment_vtable('LAND%RLONG_S', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlong_s
     endif

     if (allocated(land%rlong_v)) then
        call increment_vtable('LAND%RLONG_V', 'LW')
        vtab_r(num_var)%rvar1_p => land%rlong_v
     endif

     if (allocated(land%snowfac)) then
        call increment_vtable('LAND%SNOWFAC', 'LW')
        vtab_r(num_var)%rvar1_p => land%snowfac
     endif

     if (allocated(land%vf)) then
        call increment_vtable('LAND%VF', 'LW')
        vtab_r(num_var)%rvar1_p => land%vf
     endif

     if (allocated(land%pcpg)) then
        call increment_vtable('LAND%PCPG', 'LW')
        vtab_r(num_var)%rvar1_p => land%pcpg
     endif

     if (allocated(land%qpcpg)) then
        call increment_vtable('LAND%QPCPG', 'LW')
        vtab_r(num_var)%rvar1_p => land%qpcpg
     endif

     if (allocated(land%dpcpg)) then
        call increment_vtable('LAND%DPCPG', 'LW')
        vtab_r(num_var)%rvar1_p => land%dpcpg
     endif

     if (allocated(land%lsl)) then
        call increment_vtable('LAND%LSL', 'LW')
        vtab_r(num_var)%ivar1_p => land%lsl
     endif

     if (allocated(land%head0)) then
        call increment_vtable('LAND%HEAD0', 'LW')
        vtab_r(num_var)%rvar1_p => land%head0
     endif

     if (allocated(land%head1)) then
        call increment_vtable('LAND%HEAD1', 'LW')
        vtab_r(num_var)%rvar1_p => land%head1
     endif

     if (allocated(land%ntext_soil)) then
        call increment_vtable('LAND%NTEXT_SOIL', 'LW')
        vtab_r(num_var)%ivar2_p => land%ntext_soil
     endif

     if (allocated(land%soil_water)) then
        call increment_vtable('LAND%SOIL_WATER', 'LW')
        vtab_r(num_var)%rvar2_p => land%soil_water
     endif

     if (allocated(land%soil_energy)) then
        call increment_vtable('LAND%SOIL_ENERGY', 'LW')
        vtab_r(num_var)%rvar2_p => land%soil_energy
     endif
     
     if (allocated(land%sfcwater_mass)) then
        call increment_vtable('LAND%SFCWATER_MASS', 'LW')
        vtab_r(num_var)%rvar2_p => land%sfcwater_mass
     endif

     if (allocated(land%sfcwater_energy)) then
        call increment_vtable('LAND%SFCWATER_ENERGY', 'LW')
        vtab_r(num_var)%rvar2_p => land%sfcwater_energy
     endif

     if (allocated(land%sfcwater_depth)) then
        call increment_vtable('LAND%SFCWATER_DEPTH', 'LW')
        vtab_r(num_var)%rvar2_p => land%sfcwater_depth
     endif

     if (allocated(land%rshort_s)) then
        call increment_vtable('LAND%RSHORT_S', 'LW')
        vtab_r(num_var)%rvar2_p => land%rshort_s
     endif

   end subroutine filltab_leaf

!=========================================================================

   subroutine filltab_ED()
     use var_tables, only: increment_EDtab, num_ED, vtab_ED
     implicit none

     if (allocated(land%gpp)) then
        call increment_EDtab('LAND%GPP', hist=.true., mavg=.true., yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%gpp
     endif

     if (allocated(land%rh)) then
        call increment_EDtab('LAND%RH', hist=.true., mavg=.true., yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%rh
     endif

     if (allocated(land%nep)) then
        call increment_EDtab('LAND%NEP', hist=.true., mavg=.true., yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%nep
     endif

     if (allocated(land%agb)) then
        call increment_EDtab('LAND%AGB', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb
     endif

     if (allocated(land%basal_area)) then
        call increment_EDtab('LAND%BASAL_AREA', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb
     endif

     if (allocated(land%agb_growth)) then
        call increment_EDtab('LAND%AGB_GROWTH', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb_growth
     endif

     if (allocated(land%agb_mort)) then
        call increment_EDtab('LAND%AGB_MORT', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb_mort
     endif

     if (allocated(land%agb_cut)) then
        call increment_EDtab('LAND%AGB_CUT', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb_cut
     endif

     if (allocated(land%agb_recruit)) then
        call increment_EDtab('LAND%AGB_RECRUIT', yavg=.true.)
        vtab_ED(num_ED)%rvar1_p => land%agb_recruit
     endif

   end subroutine filltab_ED

!===============================================================================

   subroutine fill_jland()

   use mem_ijtabs, only: mrls

   use misc_coms,  only: io6, iparallel
   
   use mem_para,   only: mgroupsize, myrank,  &
!                        send_ul, recv_ul,  &
                         send_wl, recv_wl,  &
                         send_wlf, recv_wlf,  &
!                        nsends_ul, nrecvs_ul,  &
                         nsends_wl, nrecvs_wl,  &
                         nsends_wlf, nrecvs_wlf

   use leaf_coms,   only: mml, mul, mwl, nml, nul, nwl
   
   implicit none
   
   integer :: jsend, iwl, jend

! Allocate and zero-fill JTAB_WL_MPI%JEND

   do jsend = 1,maxremote
      allocate (jtab_wl_mpi(jsend)%jend(mrls))
                jtab_wl_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Return if run is not parallel (jtab not needed)

   if (iparallel == 0) return

! Compute and store JTAB_WL_MPI%JEND(1)

   do jsend = 1,nsends_wl(1)
      jtab_wl_mpi(jsend)%jend(1) = 0
      do iwl = 2,mwl
         if (itab_wl(iwl)%send(jsend)) then
            jtab_wl_mpi(jsend)%jend(1) = jtab_wl_mpi(jsend)%jend(1) + 1
         endif
      enddo
      jtab_wl_mpi(jsend)%jend(1) = max(1,jtab_wl_mpi(jsend)%jend(1))
   enddo

! Allocate and zero-fill JTAB_WL_MPI%IWL

   do jsend = 1,nsends_wl(1)
      jend = jtab_wl_mpi(jsend)%jend(1)
      allocate (jtab_wl_mpi(jsend)%iwl(jend))
                jtab_wl_mpi(jsend)%iwl(1:jend) = 0
   enddo

! Initialize JTAB_WL_MPI%JEND counters to zero

   do jsend = 1,nsends_wl(1)
      jtab_wl_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Compute JTAB_WL_MPI%IWL (only fill for MRL = 1)

   do iwl = 2,mwl
      do jsend = 1,nsends_wl(1)

         if (itab_wl(iwl)%send(jsend)) then
            jtab_wl_mpi(jsend)%jend(1) = jtab_wl_mpi(jsend)%jend(1) + 1
            jtab_wl_mpi(jsend)%iwl(jtab_wl_mpi(jsend)%jend(1)) = iwl
         endif

      enddo
   enddo

   return
   end subroutine fill_jland

End Module mem_leaf
