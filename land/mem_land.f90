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
Module mem_land

  implicit none

  integer :: nzg         ! number of soil grid levels

  integer :: nland       ! Total # of land cell W pts in model domain
  integer :: mland       ! Total # of land cell W pts in model parallel sub-domain

  integer :: onland = 0  ! Always zero, so need not use in code
  integer :: omland = 0  ! Always zero, so need not use in code

  integer :: kperc = 2   ! Vertical index of percolation level; perc depth is slz(kperc)

  real :: landgrid_dztop ! Thickness (m) of top (shallowest) soil grid level 
  real :: landgrid_depth ! Depth (m) of soil grid lower boundary

  real, save, allocatable :: slz   (:) ! height (negative) of soil layer M pt [m]
  real, save, allocatable :: dslz  (:) ! soil layer thickness at T pt [m]
  real, save, allocatable :: dslzo2(:) ! HALF soil layer thickness at T pt [m]
  real, save, allocatable :: dslzi (:) ! inverse soil layer thickness at T pt [1/m]
  real, save, allocatable :: slzt  (:) ! height (negative) of soil T pt [m]

  ! LAND GRID TABLES

  Type itab_land_vars
     integer :: iwglobe = 1
  End type

  type (itab_land_vars), allocatable, target :: itab_land(:)

  ! LAND MODEL VARIABLES

  Type land_vars

     real, allocatable :: slope_fact(:) ! orographic roughness parameter

     real, allocatable :: rshort_g  (:) ! s/w net rad flux to soil [W/m^2]
     real, allocatable :: rshort_s(:,:) ! s/w net rad flux to sfc water [W/m^2]
     real, allocatable :: rshort_v  (:) ! s/w net rad flux to veg [W/m^2]
     real, allocatable :: rlong_g   (:) ! l/w net rad flux to soil [W/m^2]
     real, allocatable :: rlong_s   (:) ! l/w net rad flux to sfc water [W/m^2]
     real, allocatable :: rlong_v   (:) ! l/w net rad flux to veg [W/m^2]

     real, allocatable :: cosz      (:) ! cosine of the solar zenith angle of the land cell

!     real, allocatable :: par         (:) ! total photosynthetic active radiation (W/m^2)
!     real, allocatable :: par_diffuse (:) ! diffuse photosynthetic active radiation (W/m^2)
      real, allocatable :: ppfd        (:) ! total photosynthetic photon flux density (uMol/m^2/s)
      real, allocatable :: ppfd_diffuse(:) ! diffuse photosynthetic photon flux density (uMol/m^2/s)
      
     ! Canopy and surface quantities:

     integer, allocatable :: nlev_sfcwater (:) ! # of active surface water levels
     real, allocatable :: sfcwater_mass  (:,:) ! surface water mass [kg/m^2]
     real, allocatable :: sfcwater_energy(:,:) ! surface water energy [J/kg]
     real, allocatable :: sfcwater_depth (:,:) ! surface water depth [m]

     real, allocatable :: hcapveg     (:) ! veg heat capacity [J/(m^2 K)]
     real, allocatable :: surface_srrv(:) ! surface saturation vapor mixing ratio [kg_vap/kg_dryair]
     real, allocatable :: ground_rrv  (:) ! soil vapor mixing ratio [kg_vap/kg_dryair]
     real, allocatable :: veg_fracarea(:) ! veg fractional area
     real, allocatable :: veg_lai     (:) ! veg leaf area index
     real, allocatable :: veg_rough   (:) ! veg roughness height [m]
     real, allocatable :: veg_height  (:) ! veg height [m]
     real, allocatable :: veg_albedo  (:) ! veg albedo
     real, allocatable :: veg_tai     (:) ! veg total area index
     real, allocatable :: veg_water   (:) ! veg sfc water content [kg/m^2]
     real, allocatable :: veg_energy  (:) ! (veg + veg_water) energy [J/m^2]
     real, allocatable :: veg_temp    (:) ! veg temp [K]
     real, allocatable :: veg_ndvip   (:) ! veg past ndvi (obs time)
     real, allocatable :: veg_ndvif   (:) ! veg future ndvi (obs time)
     real, allocatable :: veg_ndvic   (:) ! veg current ndvi
     real, allocatable :: stom_resist (:) ! veg stomatal resistance [s/m]
     real, allocatable :: snowfac     (:) ! frac veg burial by snowcover
     real, allocatable :: vf          (:) ! frac coverage of non-buried part of veg

     ! Soil quantities (time-dependent)

     real, allocatable :: soil_water (:,:) ! soil water content [vol_water/vol_tot]
     real, allocatable :: soil_energy(:,:) ! soil energy [J/m^3]
     real, allocatable :: head       (:,:) ! hydraulic head [m]
     real, allocatable :: head0        (:) ! LBC total hydraulic head [m]

     ! Soil quantities read from SoilGrids datasets (constant in time)

     ! In future, may also add "srsthick" fields

     integer, allocatable :: usdatext       (:) ! USDA soil textural class converted from FAO

     real, allocatable :: z_bedrock         (:) ! height (non-positive) of top of bedrock [m]
     real, allocatable :: gpp               (:) ! gross primary production (of carbon) [gC/(m^2 yr)]
     real, allocatable :: glhymps_ksat      (:) ! glhymps ksat [m/s]
     real, allocatable :: glhymps_ksat_pfr  (:) ! glhymps ksat with permafrost [m/s]
     real, allocatable :: glhymps_poros     (:) ! glhymps porosity []

     real, allocatable :: sand            (:,:) ! soil sand fraction []
     real, allocatable :: clay            (:,:) ! soil clay fraction []
     real, allocatable :: silt            (:,:) ! soil silt fraction []
     real, allocatable :: organ           (:,:) ! soil organic specific density [kg/kg]
     real, allocatable :: bulkdens_drysoil(:,:) ! soil bulk density [kg/m^3]
     real, allocatable :: pH_soil         (:,:) ! soil pH []
     real, allocatable :: cec_soil        (:,:) ! soil cation exchange capacity [cmol+/kg]

     ! Soil quantities computed from SoilGrids-read parameters (constant in time)

     real, allocatable :: wresid_vg          (:,:) ! residual water content (vG) []
     real, allocatable :: wsat_vg            (:,:) ! saturation water content (porosity) []
     real, allocatable :: ksat_vg            (:,:) ! saturation hydraulic conductivity [m/s]
     real, allocatable :: alpha_vg           (:,:) ! van Genuchten alpha term [1/m]
     real, allocatable :: en_vg              (:,:) ! van Genuchten N term []
     real, allocatable :: lambda_vg          (:,:) ! van Genuchten lambda term []
     real, allocatable :: specifheat_drysoil (:,:) ! Specific heat of dry soil [J/(m^3 K)]

  End Type land_vars

  type (land_vars) :: land

Contains

!=========================================================================

  subroutine alloc_landcol()

  implicit none

  ! Allocate leaf column arrays

  allocate (slz   (nzg+1))
  allocate (dslz  (nzg))
  allocate (dslzo2(nzg))
  allocate (dslzi (nzg))
  allocate (slzt  (nzg))

  end subroutine alloc_landcol

!=========================================================================

  subroutine alloc_land(mland, nzg, nzs)

  use misc_coms,  only: rinit
  use oname_coms, only: nl

  implicit none

  integer, intent(in) :: mland, nzg, nzs

  ! Surface/vegetation quantities

  allocate (land%slope_fact         (mland)) ; land%slope_fact      = 1.0

  allocate (land%rshort_g           (mland)) ; land%rshort_g        = rinit
  allocate (land%rshort_s       (nzs,mland)) ; land%rshort_s        = rinit
  allocate (land%rshort_v           (mland)) ; land%rshort_v        = rinit
  allocate (land%rlong_g            (mland)) ; land%rlong_g         = rinit
  allocate (land%rlong_s            (mland)) ; land%rlong_s         = rinit
  allocate (land%rlong_v            (mland)) ; land%rlong_v         = rinit

  allocate (land%cosz               (mland)) ; land%cosz            = rinit

! allocate (land%par                (mland)) ; land%par             = rinit
! allocate (land%par_diffuse        (mland)) ; land%par_diffuse     = rinit
  allocate (land%ppfd               (mland)) ; land%ppfd            = rinit
  allocate (land%ppfd_diffuse       (mland)) ; land%ppfd_diffuse    = rinit

  allocate (land%nlev_sfcwater      (mland)) ; land%nlev_sfcwater   = 0
  allocate (land%sfcwater_mass  (nzs,mland)) ; land%sfcwater_mass   = rinit
  allocate (land%sfcwater_energy(nzs,mland)) ; land%sfcwater_energy = rinit
  allocate (land%sfcwater_depth (nzs,mland)) ; land%sfcwater_depth  = rinit

  allocate (land%hcapveg            (mland)) ; land%hcapveg         = rinit
  allocate (land%surface_srrv       (mland)) ; land%surface_srrv    = rinit
  allocate (land%ground_rrv         (mland)) ; land%ground_rrv      = rinit
  allocate (land%veg_fracarea       (mland)) ; land%veg_fracarea    = rinit
  allocate (land%veg_lai            (mland)) ; land%veg_lai         = rinit
  allocate (land%veg_rough          (mland)) ; land%veg_rough       = rinit
  allocate (land%veg_height         (mland)) ; land%veg_height      = rinit
  allocate (land%veg_albedo         (mland)) ; land%veg_albedo      = rinit
  allocate (land%veg_tai            (mland)) ; land%veg_tai         = rinit
  allocate (land%veg_water          (mland)) ; land%veg_water       = rinit
  allocate (land%veg_energy         (mland)) ; land%veg_energy      = rinit
  allocate (land%veg_temp           (mland)) ; land%veg_temp        = rinit
  allocate (land%veg_ndvip          (mland)) ; land%veg_ndvip       = rinit
  allocate (land%veg_ndvif          (mland)) ; land%veg_ndvif       = rinit
  allocate (land%veg_ndvic          (mland)) ; land%veg_ndvic       = rinit
  allocate (land%stom_resist        (mland)) ; land%stom_resist     = rinit
  allocate (land%snowfac            (mland)) ; land%snowfac         = rinit
  allocate (land%vf                 (mland)) ; land%vf              = rinit

  ! Soil quantities

  allocate (land%soil_water     (nzg,mland)) ; land%soil_water  = rinit
  allocate (land%soil_energy    (nzg,mland)) ; land%soil_energy = rinit
  allocate (land%head           (nzg,mland)) ; land%head        = rinit
  allocate (land%head0              (mland)) ; land%head0       = rinit

 ! The following are allocated in makesfc3.f90 because they are on the SFCGRIDFILE

 ! allocate (land%usdatext            (mland)) ; land%usdatext         = 0
 ! allocate (land%z_bedrock           (mland)) ; land%z_bedrock        = rinit
 ! allocate (land%gpp                 (mland)) ; land%gpp              = rinit
 ! allocate (land%sand            (nzg,mland)) ; land%sand             = rinit
 ! allocate (land%clay            (nzg,mland)) ; land%clay             = rinit
 ! allocate (land%silt            (nzg,mland)) ; land%silt             = rinit
 ! allocate (land%organ           (nzg,mland)) ; land%organ            = rinit
 ! allocate (land%bulkdens_drysoil(nzg,mland)) ; land%bulkdens_drysoil = rinit
 ! allocate (land%pH_soil         (nzg,mland)) ; land%pH_soil          = rinit
 ! allocate (land%cec_soil        (nzg,mland)) ; land%cec_soil         = rinit

  allocate (land%wresid_vg         (nzg,mland)) ; land%wresid_vg          = rinit
  allocate (land%wsat_vg           (nzg,mland)) ; land%wsat_vg            = rinit
  allocate (land%ksat_vg           (nzg,mland)) ; land%ksat_vg            = rinit
  allocate (land%alpha_vg          (nzg,mland)) ; land%alpha_vg           = rinit
  allocate (land%en_vg             (nzg,mland)) ; land%en_vg              = rinit
  allocate (land%lambda_vg         (nzg,mland)) ; land%lambda_vg          = rinit
  allocate (land%specifheat_drysoil(nzg,mland)) ; land%specifheat_drysoil = rinit

  end subroutine alloc_land

!=========================================================================

  subroutine filltab_land()

  use var_tables, only: increment_vtable

  implicit none

  ! Surface/vegetation quantities

  if (allocated(land%rshort_g))        call increment_vtable('LAND%RSHORT_G',        'LW', rvar1=land%rshort_g)
  if (allocated(land%rshort_s))        call increment_vtable('LAND%RSHORT_S',        'LW', rvar2=land%rshort_s)
  if (allocated(land%rshort_v))        call increment_vtable('LAND%RSHORT_V',        'LW', rvar1=land%rshort_v)
  if (allocated(land%rlong_g))         call increment_vtable('LAND%RLONG_G',         'LW', rvar1=land%rlong_g)
  if (allocated(land%rlong_s))         call increment_vtable('LAND%RLONG_S',         'LW', rvar1=land%rlong_s)
  if (allocated(land%rlong_v))         call increment_vtable('LAND%RLONG_V',         'LW', rvar1=land%rlong_v)
  if (allocated(land%cosz))            call increment_vtable('LAND%COSZ',            'LW', rvar1=land%cosz)
  if (allocated(land%nlev_sfcwater))   call increment_vtable('LAND%NLEV_SFCWATER',   'LW', ivar1=land%nlev_sfcwater)
  if (allocated(land%sfcwater_mass))   call increment_vtable('LAND%SFCWATER_MASS',   'LW', rvar2=land%sfcwater_mass)
  if (allocated(land%sfcwater_energy)) call increment_vtable('LAND%SFCWATER_ENERGY', 'LW', rvar2=land%sfcwater_energy)
  if (allocated(land%sfcwater_depth))  call increment_vtable('LAND%SFCWATER_DEPTH',  'LW', rvar2=land%sfcwater_depth)
  if (allocated(land%hcapveg))         call increment_vtable('LAND%HCAPVEG',         'LW', rvar1=land%hcapveg)
  if (allocated(land%surface_srrv))    call increment_vtable('LAND%SURFACE_SRRV',    'LW', rvar1=land%surface_srrv)
  if (allocated(land%ground_rrv))      call increment_vtable('LAND%GROUND_RRV',      'LW', rvar1=land%ground_rrv)
  if (allocated(land%veg_fracarea))    call increment_vtable('LAND%VEG_FRACAREA',    'LW', rvar1=land%veg_fracarea)
  if (allocated(land%veg_lai))         call increment_vtable('LAND%VEG_LAI',         'LW', rvar1=land%veg_lai)
  if (allocated(land%veg_rough))       call increment_vtable('LAND%VEG_ROUGH',       'LW', rvar1=land%veg_rough)
  if (allocated(land%veg_height))      call increment_vtable('LAND%VEG_HEIGHT',      'LW', rvar1=land%veg_height)
  if (allocated(land%veg_albedo))      call increment_vtable('LAND%VEG_ALBEDO',      'LW', rvar1=land%veg_albedo)
  if (allocated(land%veg_tai))         call increment_vtable('LAND%VEG_TAI',         'LW', rvar1=land%veg_tai)
  if (allocated(land%veg_water))       call increment_vtable('LAND%VEG_WATER',       'LW', rvar1=land%veg_water)
  if (allocated(land%veg_energy))      call increment_vtable('LAND%VEG_ENERGY',      'LW', rvar1=land%veg_energy)
  if (allocated(land%veg_temp))        call increment_vtable('LAND%VEG_TEMP',        'LW', rvar1=land%veg_temp)
  if (allocated(land%veg_ndvic))       call increment_vtable('LAND%VEG_NDVIC',       'LW', rvar1=land%veg_ndvic)
  if (allocated(land%stom_resist))     call increment_vtable('LAND%STOM_RESIST',     'LW', rvar1=land%stom_resist)
  if (allocated(land%snowfac))         call increment_vtable('LAND%SNOWFAC',         'LW', rvar1=land%snowfac)
  if (allocated(land%vf))              call increment_vtable('LAND%VF',              'LW', rvar1=land%vf)

! if (allocated(land%par))          call increment_vtable('LAND%PAR',          'LW', rvar1=land%par)
! if (allocated(land%par_diffuse))  call increment_vtable('LAND%PAR_DIFFUSE',  'LW', rvar1=land%par_diffuse)
  if (allocated(land%ppfd))         call increment_vtable('LAND%PPFD',         'LW', rvar1=land%ppfd)
  if (allocated(land%ppfd_diffuse)) call increment_vtable('LAND%PPFD_DIFFUSE', 'LW', rvar1=land%ppfd_diffuse)
     
  ! Soil quantities

  if (allocated(land%soil_water))  call increment_vtable('LAND%SOIL_WATER',  'LW', rvar2=land%soil_water)
  if (allocated(land%soil_energy)) call increment_vtable('LAND%SOIL_ENERGY', 'LW', rvar2=land%soil_energy)
  if (allocated(land%head))        call increment_vtable('LAND%HEAD',        'LW', rvar2=land%head)
  if (allocated(land%head0))       call increment_vtable('LAND%HEAD0',       'LW', rvar1=land%head0)

  end subroutine filltab_land

End Module mem_land
