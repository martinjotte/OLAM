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

  real, allocatable :: slz   (:) ! height (negative) of soil layer M pt [m]
  real, allocatable :: dslz  (:) ! soil layer thickness at T pt [m]
  real, allocatable :: dslzo2(:) ! HALF soil layer thickness at T pt [m]
  real, allocatable :: dslzi (:) ! inverse soil layer thickness at T pt [1/m]
  real, allocatable :: slzt  (:) ! height (negative) of soil T pt [m]

  ! LAND GRID TABLES

  Type itab_land_vars
     integer :: iwglobe
  End type itab_land_vars

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

      real, allocatable :: par         (:) ! total photosynthetic active radiation (W/m^2)
      real, allocatable :: par_diffuse (:) ! diffuse photosynthetic active radiation (W/m^2)
      real, allocatable :: ppfd        (:) ! total photosynthetic photon flux density (uMol/m^2/s)
      real, allocatable :: ppfd_diffuse(:) ! diffuse photosynthetic photon flux density (uMol/m^2/s)

     ! Canopy and surface quantities:

     integer, allocatable :: nlev_sfcwater (:) ! # of active surface water levels
     real, allocatable :: sfcwater_mass  (:,:) ! surface water mass [kg/m^2]
     real, allocatable :: sfcwater_energy(:,:) ! surface water energy [J/kg]
     real, allocatable :: sfcwater_depth (:,:) ! surface water depth [m]

     real, allocatable :: hcapveg     (:) ! veg heat capacity [J/(m^2 K)]
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
!    real, allocatable :: glhymps_ksat_pfr  (:) ! glhymps ksat with permafrost [m/s]
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
     real, allocatable :: wfrac_low          (:,:)

     ! Other constant derived soil quantities

     integer, allocatable :: k_bedrock (:) ! level of the heighest bedrock/transition layer
     real,    allocatable :: soilfldcap(:) ! field capacity [Vol/Vol] of top soil layer
     real,    allocatable :: soilwilt  (:) ! wilting point [Vol/Vol] of top soil layer

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

  subroutine alloc_land1(mland0)

    use misc_coms,  only: rinit
    implicit none

    integer, intent(in) :: mland0
    integer             :: iland

    ! This routine allocates sea arrays necessary for the MAKEGRID stage
    ! or for reading sfcgrid information

    allocate( itab_land                (mland0) )
    allocate( land%usdatext            (mland0) )
    allocate( land%z_bedrock           (mland0) )
    allocate( land%gpp                 (mland0) )

    allocate( land%glhymps_ksat        (mland0) )
!   allocate( land%glhymps_ksat_pfr    (mland0) )
    allocate( land%glhymps_poros       (mland0) )

    allocate( land%sand            (nzg,mland0) )
    allocate( land%clay            (nzg,mland0) )
    allocate( land%silt            (nzg,mland0) )
    allocate( land%organ           (nzg,mland0) )
    allocate( land%bulkdens_drysoil(nzg,mland0) )
    allocate( land%pH_soil         (nzg,mland0) )
    allocate( land%cec_soil        (nzg,mland0) )

    !$omp parallel do
    do iland = 1, mland0
       itab_land              (iland) = itab_land_vars( iwglobe = 1 )
       land%usdatext          (iland) = 0
       land%z_bedrock         (iland) = rinit
       land%gpp               (iland) = rinit
       land%glhymps_ksat      (iland) = rinit
!      land%glhymps_ksat_pfr  (iland) = rinit
       land%glhymps_poros     (iland) = rinit
       land%sand            (:,iland) = rinit
       land%clay            (:,iland) = rinit
       land%silt            (:,iland) = rinit
       land%organ           (:,iland) = rinit
       land%bulkdens_drysoil(:,iland) = rinit
       land%pH_soil         (:,iland) = rinit
       land%cec_soil        (:,iland) = rinit
    enddo
    !$omp end parallel do

  end subroutine alloc_land1

!=========================================================================

  subroutine alloc_land2()

  use misc_coms, only: rinit, do_chem
  use leaf_coms, only: nzs

  implicit none

  integer :: iland

  ! Surface/vegetation quantities

  allocate( land%slope_fact            (mland) )

  allocate( land%rshort_g              (mland) )
  allocate( land%rshort_s          (nzs,mland) )
  allocate( land%rshort_v              (mland) )
  allocate( land%rlong_g               (mland) )
  allocate( land%rlong_s               (mland) )
  allocate( land%rlong_v               (mland) )
  allocate( land%cosz                  (mland) )

  ! Photosynthetically active radiation (PAR) currently unused
! allocate( land%par                   (mland) )
! allocate( land%par_diffuse           (mland) )

  ! Photosynthetic photon flux density (PPFD) only needed with full chemistry
  if (do_chem == 1) then
     allocate( land%ppfd               (mland) )
     allocate( land%ppfd_diffuse       (mland) )
  endif

  allocate( land%nlev_sfcwater         (mland) )
  allocate( land%sfcwater_mass     (nzs,mland) )
  allocate( land%sfcwater_energy   (nzs,mland) )
  allocate( land%sfcwater_depth    (nzs,mland) )

  allocate( land%hcapveg               (mland) )
  allocate( land%veg_fracarea          (mland) )
  allocate( land%veg_lai               (mland) )
  allocate( land%veg_rough             (mland) )
  allocate( land%veg_height            (mland) )
  allocate( land%veg_albedo            (mland) )
  allocate( land%veg_tai               (mland) )
  allocate( land%veg_water             (mland) )
  allocate( land%veg_energy            (mland) )
  allocate( land%veg_temp              (mland) )
  allocate( land%veg_ndvip             (mland) )
  allocate( land%veg_ndvif             (mland) )
  allocate( land%veg_ndvic             (mland) )
  allocate( land%stom_resist           (mland) )
  allocate( land%snowfac               (mland) )
  allocate( land%vf                    (mland) )

  ! Soil quantities

  allocate( land%soil_water        (nzg,mland) )
  allocate( land%soil_energy       (nzg,mland) )
  allocate( land%head              (nzg,mland) )
  allocate( land%head0                 (mland) )
  allocate( land%wresid_vg         (nzg,mland) )
  allocate( land%wsat_vg           (nzg,mland) )
  allocate( land%ksat_vg           (nzg,mland) )
  allocate( land%alpha_vg          (nzg,mland) )
  allocate( land%en_vg             (nzg,mland) )
  allocate( land%lambda_vg         (nzg,mland) )
  allocate( land%specifheat_drysoil(nzg,mland) )
  allocate( land%wfrac_low         (nzg,mland) )
  allocate( land%k_bedrock             (mland) )
  allocate( land%soilfldcap            (mland) )
  allocate( land%soilwilt              (mland) )

  !$omp parallel do
  do iland = 1, mland

     if ( allocated( land%slope_fact         ) ) land%slope_fact          (iland) = 1.0

     if ( allocated( land%rshort_g           ) ) land%rshort_g            (iland) = rinit
     if ( allocated( land%rshort_s           ) ) land%rshort_s          (:,iland) = rinit
     if ( allocated( land%rshort_v           ) ) land%rshort_v            (iland) = rinit
     if ( allocated( land%rlong_g            ) ) land%rlong_g             (iland) = rinit
     if ( allocated( land%rlong_s            ) ) land%rlong_s             (iland) = rinit
     if ( allocated( land%rlong_v            ) ) land%rlong_v             (iland) = rinit
     if ( allocated( land%cosz               ) ) land%cosz                (iland) = rinit

     if ( allocated( land%par                ) ) land%par                 (iland) = rinit
     if ( allocated( land%par_diffuse        ) ) land%par_diffuse         (iland) = rinit
     if ( allocated( land%ppfd               ) ) land%ppfd                (iland) = rinit
     if ( allocated( land%ppfd_diffuse       ) ) land%ppfd_diffuse        (iland) = rinit

     if ( allocated( land%nlev_sfcwater      ) ) land%nlev_sfcwater       (iland) = 0
     if ( allocated( land%sfcwater_mass      ) ) land%sfcwater_mass     (:,iland) = rinit
     if ( allocated( land%sfcwater_energy    ) ) land%sfcwater_energy   (:,iland) = rinit
     if ( allocated( land%sfcwater_depth     ) ) land%sfcwater_depth    (:,iland) = rinit

     if ( allocated( land%hcapveg            ) ) land%hcapveg             (iland) = rinit
     if ( allocated( land%veg_fracarea       ) ) land%veg_fracarea        (iland) = rinit
     if ( allocated( land%veg_lai            ) ) land%veg_lai             (iland) = rinit
     if ( allocated( land%veg_rough          ) ) land%veg_rough           (iland) = rinit
     if ( allocated( land%veg_height         ) ) land%veg_height          (iland) = rinit
     if ( allocated( land%veg_albedo         ) ) land%veg_albedo          (iland) = rinit
     if ( allocated( land%veg_tai            ) ) land%veg_tai             (iland) = rinit
     if ( allocated( land%veg_water          ) ) land%veg_water           (iland) = rinit
     if ( allocated( land%veg_energy         ) ) land%veg_energy          (iland) = rinit
     if ( allocated( land%veg_temp           ) ) land%veg_temp            (iland) = rinit
     if ( allocated( land%veg_ndvip          ) ) land%veg_ndvip           (iland) = rinit
     if ( allocated( land%veg_ndvif          ) ) land%veg_ndvif           (iland) = rinit
     if ( allocated( land%veg_ndvic          ) ) land%veg_ndvic           (iland) = rinit
     if ( allocated( land%stom_resist        ) ) land%stom_resist         (iland) = rinit
     if ( allocated( land%snowfac            ) ) land%snowfac             (iland) = rinit
     if ( allocated( land%vf                 ) ) land%vf                  (iland) = rinit

     ! Soil quantities

     if ( allocated( land%soil_water         ) ) land%soil_water        (:,iland) = rinit
     if ( allocated( land%soil_energy        ) ) land%soil_energy       (:,iland) = rinit
     if ( allocated( land%head               ) ) land%head              (:,iland) = rinit
     if ( allocated( land%head0              ) ) land%head0               (iland) = rinit

     if ( allocated( land%wresid_vg          ) ) land%wresid_vg         (:,iland) = rinit
     if ( allocated( land%wsat_vg            ) ) land%wsat_vg           (:,iland) = rinit
     if ( allocated( land%ksat_vg            ) ) land%ksat_vg           (:,iland) = rinit
     if ( allocated( land%alpha_vg           ) ) land%alpha_vg          (:,iland) = rinit
     if ( allocated( land%en_vg              ) ) land%en_vg             (:,iland) = rinit
     if ( allocated( land%lambda_vg          ) ) land%lambda_vg         (:,iland) = rinit
     if ( allocated( land%specifheat_drysoil ) ) land%specifheat_drysoil(:,iland) = rinit
     if ( allocated( land%wfrac_low          ) ) land%wfrac_low         (:,iland) = rinit
     if ( allocated( land%k_bedrock          ) ) land%k_bedrock           (iland) = 0
     if ( allocated( land%soilfldcap         ) ) land%soilfldcap          (iland) = rinit
     if ( allocated( land%soilwilt           ) ) land%soilwilt            (iland) = rinit

  enddo
  !$omp end parallel do

  end subroutine alloc_land2

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
