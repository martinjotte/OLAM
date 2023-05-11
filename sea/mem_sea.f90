Module mem_sea

   use max_dims, only: maxgrds, maxngrdll

   implicit none

   integer :: nzsea = 0 ! Max number of vertical levels in sea/ocean model

   integer :: nsea      ! Total # of sea cell W pts in model global domain
   integer :: msea      ! Total # of sea cell W pts in model parallel sub-domain

   integer :: onsea     ! Offset # of sea cell W pts in model global domain
   integer :: omsea     ! Offset # of sea cell W pts in model parallel sub-domain
                        ! (isea + omsea = iwsfc)

! SEA GRID TABLES

   Type itab_sea_vars
      integer :: iwglobe = 1
   End type

   type (itab_sea_vars), allocatable, target :: itab_sea(:)

! SEA MODEL VARIABLES

   Type sea_vars

! Canopy to atmosphere turbulent flux quantities

      real, allocatable :: sea_ustar  (:) ! friction velocity [m/s]
      real, allocatable :: ice_ustar  (:) ! friction velocity [m/s]

      real, allocatable :: sea_vkmsfc (:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: ice_vkmsfc (:) ! surface drag coefficient [kg/(m s)]

      real, allocatable :: sea_sfluxt (:)
      real, allocatable :: ice_sfluxt (:)

      real, allocatable :: sea_sfluxr (:)
      real, allocatable :: ice_sfluxr (:)

      real, allocatable :: sea_sfluxc (:)
      real, allocatable :: ice_sfluxc (:)

      real, allocatable :: sea_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: ice_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]

      real, allocatable :: sea_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: ice_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]

      real, allocatable :: sea_sxfer_c(:) ! can_air-to-atm CO2 xfer [ppm/m^2]
      real, allocatable :: ice_sxfer_c(:) ! can_air-to-atm CO2 xfer [ppm/m^2]

      real, allocatable :: sea_ggaer  (:) ! surface aerodynamic conductance over water [m/s]
      real, allocatable :: ice_ggaer  (:) ! surface aerodynamic conductance over ice [m/s]

      real, allocatable :: sea_wthv   (:) ! surface buoyancy flux over water [K m/s]
      real, allocatable :: ice_wthv   (:) ! surface buoyancy flux over ice [K m/s]

      real, allocatable :: windxe     (:) ! XE wind component [m/s] for sfc stress
      real, allocatable :: windye     (:) ! YE wind component [m/s] for sfc stress
      real, allocatable :: windze     (:) ! ZE wind component [m/s] for sfc stress

! Radiative flux quantities

      real, allocatable :: sea_albedo    (:) ! water s/w albedo [0-1]
      real, allocatable :: ice_albedo    (:) ! seaice s/w albedo [0-1]

      real, allocatable :: sea_rlongup   (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: ice_rlongup   (:) ! upward can-top l/w flux [W/m^2]

      real, allocatable :: ice_net_rlong (:)
      real, allocatable :: ice_net_rshort(:)

! Canopy and surface quantities:

      real, allocatable :: sea_cantemp   (:) ! "canopy" air temperature [K]
      real, allocatable :: ice_cantemp   (:) ! "canopy" air temperature [K]

      real, allocatable :: sea_canrrv    (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: ice_canrrv    (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]

      integer, allocatable :: nlev_seaice(:) ! number of seaice layers in current cell

      real, allocatable :: seatp         (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf         (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc         (:) ! current sea temperature [K]
      real, allocatable :: seaicep       (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef       (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec       (:) ! seaice fraction [0-1]

      real, allocatable :: surface_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: sea_sfc_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: ice_sfc_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]

      real, allocatable :: sea_rough     (:) ! water surface roughness height [m]
      real, allocatable :: ice_rough     (:) ! water surface roughness height [m]

      real, allocatable :: seaice_energy(:,:)
      real, allocatable :: seaice_tempk (:,:)

      ! Shallow Water Model (SWM) quantities:

      real, allocatable :: swmdepth(:) ! Depth of water at T point [m]
      real, allocatable :: vxe(:)
      real, allocatable :: vye(:)
      real, allocatable :: vze(:)

      real, allocatable :: vmxet(:) ! Tendency of (velocity * depth) [m^2/s^2]
      real, allocatable :: vmyet(:) ! Tendency of (velocity * depth) [m^2/s^2]
      real, allocatable :: vmzet(:) ! Tendency of (velocity * depth) [m^2/s^2]

      real, allocatable :: vmxet_area(:) ! Advective tendency of (velocity * depth * area) [m^4/s^2]
      real, allocatable :: vmyet_area(:) ! Advective tendency of (velocity * depth * area) [m^4/s^2]
      real, allocatable :: vmzet_area(:) ! Advective tendency of (velocity * depth * area) [m^4/s^2]

      real, allocatable :: vxe1(:)
      real, allocatable :: vye1(:)
      real, allocatable :: vze1(:)

      real, allocatable :: gxps_vxe(:)
      real, allocatable :: gyps_vxe(:)

      real, allocatable :: gxps_vye(:)
      real, allocatable :: gyps_vye(:)

      real, allocatable :: gxps_vze(:)
      real, allocatable :: gyps_vze(:)

      ! 1-D Princeton Ocean Model (POM1D) active flag:

      logical, allocatable :: pom_active(:) ! POM1D model active in these sfcg cells

   End Type sea_vars

   type (sea_vars) :: sea

! POM1D ZONE INFORMATION

  integer :: npomzons

  integer, target :: npomzonll(maxgrds)
  real   , target :: pomzonrad(maxgrds,maxngrdll)
  real   , target :: pomzonlat(maxgrds,maxngrdll)
  real   , target :: pomzonlon(maxgrds,maxngrdll)

Contains

!=========================================================================

   subroutine alloc_sea(msea)

     use misc_coms, only: rinit
     use sea_coms,  only: nzi
     use mem_co2,   only: co2flag
     use mem_sfcg,  only: nswmzons

     implicit none

     integer, intent(in) :: msea

!    Allocate and initialize sea arrays

     allocate (sea%sea_ustar     (msea)) ; sea%sea_ustar      = rinit
     allocate (sea%ice_ustar     (msea)) ; sea%ice_ustar      = rinit

     allocate (sea%sea_vkmsfc    (msea)) ; sea%sea_vkmsfc     = rinit
     allocate (sea%ice_vkmsfc    (msea)) ; sea%ice_vkmsfc     = rinit

     allocate (sea%sea_sfluxt    (msea)) ; sea%sea_sfluxt     = rinit
     allocate (sea%ice_sfluxt    (msea)) ; sea%ice_sfluxt     = rinit

     allocate (sea%sea_sfluxr    (msea)) ; sea%sea_sfluxr     = rinit
     allocate (sea%ice_sfluxr    (msea)) ; sea%ice_sfluxr     = rinit

     allocate (sea%sea_sxfer_t   (msea)) ; sea%sea_sxfer_t    = 0.0
     allocate (sea%ice_sxfer_t   (msea)) ; sea%ice_sxfer_t    = 0.0

     allocate (sea%sea_sxfer_r   (msea)) ; sea%sea_sxfer_r    = 0.0
     allocate (sea%ice_sxfer_r   (msea)) ; sea%ice_sxfer_r    = 0.0

     if (co2flag /= 0) then
        allocate (sea%sea_sfluxc (msea)) ; sea%sea_sfluxc     = rinit
        allocate (sea%ice_sfluxc (msea)) ; sea%ice_sfluxc     = rinit

        allocate (sea%sea_sxfer_c(msea)) ; sea%sea_sxfer_c    = 0.0
        allocate (sea%ice_sxfer_c(msea)) ; sea%ice_sxfer_c    = 0.0
     endif

     allocate (sea%sea_ggaer     (msea)) ; sea%sea_ggaer      = rinit
     allocate (sea%ice_ggaer     (msea)) ; sea%ice_ggaer      = rinit

     allocate (sea%sea_wthv      (msea)) ; sea%sea_wthv       = rinit
     allocate (sea%ice_wthv      (msea)) ; sea%ice_wthv       = rinit

     allocate (sea%windxe        (msea)) ; sea%windxe         = rinit
     allocate (sea%windye        (msea)) ; sea%windye         = rinit
     allocate (sea%windze        (msea)) ; sea%windze         = rinit

     allocate (sea%sea_albedo    (msea)) ; sea%sea_albedo     = 0.0
     allocate (sea%ice_albedo    (msea)) ; sea%ice_albedo     = 0.0

     allocate (sea%sea_rlongup   (msea)) ; sea%sea_rlongup    = 0.0
     allocate (sea%ice_rlongup   (msea)) ; sea%ice_rlongup    = 0.0

     allocate (sea%ice_net_rlong (msea)) ; sea%ice_net_rlong  = 0.0
     allocate (sea%ice_net_rshort(msea)) ; sea%ice_net_rshort = 0.0

     allocate (sea%sea_cantemp   (msea)) ; sea%sea_cantemp    = rinit
     allocate (sea%ice_cantemp   (msea)) ; sea%ice_cantemp    = rinit

     allocate (sea%sea_canrrv    (msea)) ; sea%sea_canrrv     = rinit
     allocate (sea%ice_canrrv    (msea)) ; sea%ice_canrrv     = rinit

     allocate (sea%nlev_seaice   (msea)) ; sea%nlev_seaice    = 0

     allocate (sea%seatp         (msea)) ; sea%seatp          = rinit
     allocate (sea%seatf         (msea)) ; sea%seatf          = rinit
     allocate (sea%seatc         (msea)) ; sea%seatc          = rinit
     allocate (sea%seaicep       (msea)) ; sea%seaicep        = rinit
     allocate (sea%seaicef       (msea)) ; sea%seaicef        = rinit
     allocate (sea%seaicec       (msea)) ; sea%seaicec        = rinit

     allocate (sea%surface_srrv  (msea)) ; sea%surface_srrv   = rinit
     allocate (sea%sea_sfc_srrv  (msea)) ; sea%sea_sfc_srrv   = rinit
     allocate (sea%ice_sfc_srrv  (msea)) ; sea%ice_sfc_srrv   = rinit

     allocate (sea%sea_rough     (msea)) ; sea%sea_rough      = rinit
     allocate (sea%ice_rough     (msea)) ; sea%ice_rough      = rinit

     allocate (sea%seaice_energy(nzi,msea)) ; sea%seaice_energy = rinit
     allocate (sea%seaice_tempk (nzi,msea)) ; sea%seaice_tempk  = rinit

     if (nswmzons > 0) then

        allocate (sea%swmdepth      (msea)) ; sea%swmdepth    = 0.0

        allocate (sea%vxe           (msea)) ; sea%vxe         = 0.0
        allocate (sea%vye           (msea)) ; sea%vye         = 0.0
        allocate (sea%vze           (msea)) ; sea%vze         = 0.0

        allocate (sea%vmxet         (msea)) ; sea%vmxet       = 0.0
        allocate (sea%vmyet         (msea)) ; sea%vmyet       = 0.0
        allocate (sea%vmzet         (msea)) ; sea%vmzet       = 0.0

        allocate (sea%vmxet_area    (msea)) ; sea%vmxet_area  = 0.0
        allocate (sea%vmyet_area    (msea)) ; sea%vmyet_area  = 0.0
        allocate (sea%vmzet_area    (msea)) ; sea%vmzet_area  = 0.0

        allocate (sea%vxe1          (msea)) ; sea%vxe1        = 0.0
        allocate (sea%vye1          (msea)) ; sea%vye1        = 0.0
        allocate (sea%vze1          (msea)) ; sea%vze1        = 0.0

        allocate (sea%gxps_vxe      (msea)) ; sea%gxps_vxe    = 0.0
        allocate (sea%gyps_vxe      (msea)) ; sea%gyps_vxe    = 0.0

        allocate (sea%gxps_vye      (msea)) ; sea%gxps_vye    = 0.0
        allocate (sea%gyps_vye      (msea)) ; sea%gyps_vye    = 0.0

        allocate (sea%gxps_vze      (msea)) ; sea%gxps_vze    = 0.0
        allocate (sea%gyps_vze      (msea)) ; sea%gyps_vze    = 0.0

     endif

   end subroutine alloc_sea

!=========================================================================

   subroutine filltab_sea()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(sea%sea_ustar))      call increment_vtable('SEA%SEA_USTAR',      'SW', rvar1=sea%sea_ustar)
     if (allocated(sea%ice_ustar))      call increment_vtable('SEA%ICE_USTAR',      'SW', rvar1=sea%ice_ustar)
     if (allocated(sea%sea_vkmsfc))     call increment_vtable('SEA%SEA_VKMSFC',     'SW', rvar1=sea%sea_vkmsfc)
     if (allocated(sea%ice_vkmsfc))     call increment_vtable('SEA%ICE_VKMSFC',     'SW', rvar1=sea%ice_vkmsfc)
     if (allocated(sea%sea_sfluxt))     call increment_vtable('SEA%SEA_SFLUXT',     'SW', rvar1=sea%sea_sfluxt)
     if (allocated(sea%ice_sfluxt))     call increment_vtable('SEA%ICE_SFLUXT',     'SW', rvar1=sea%ice_sfluxt)
     if (allocated(sea%sea_sfluxr))     call increment_vtable('SEA%SEA_SFLUXR',     'SW', rvar1=sea%sea_sfluxr)
     if (allocated(sea%ice_sfluxr))     call increment_vtable('SEA%ICE_SFLUXR',     'SW', rvar1=sea%ice_sfluxr)
     if (allocated(sea%sea_sfluxc))     call increment_vtable('SEA%SEA_SFLUXC',     'SW', rvar1=sea%sea_sfluxc)
     if (allocated(sea%ice_sfluxc))     call increment_vtable('SEA%ICE_SFLUXC',     'SW', rvar1=sea%ice_sfluxc)
     if (allocated(sea%sea_sxfer_t))    call increment_vtable('SEA%SEA_SXFER_T',    'SW', rvar1=sea%sea_sxfer_t)
     if (allocated(sea%ice_sxfer_t))    call increment_vtable('SEA%ICE_SXFER_T',    'SW', rvar1=sea%ice_sxfer_t)
     if (allocated(sea%sea_sxfer_r))    call increment_vtable('SEA%SEA_SXFER_R',    'SW', rvar1=sea%sea_sxfer_r)
     if (allocated(sea%ice_sxfer_r))    call increment_vtable('SEA%ICE_SXFER_R',    'SW', rvar1=sea%ice_sxfer_r)
     if (allocated(sea%sea_sxfer_c))    call increment_vtable('SEA%SEA_SXFER_C',    'SW', rvar1=sea%sea_sxfer_c)
     if (allocated(sea%ice_sxfer_c))    call increment_vtable('SEA%ICE_SXFER_C',    'SW', rvar1=sea%ice_sxfer_c)
     if (allocated(sea%sea_wthv))       call increment_vtable('SEA%SEA_WTHV',       'SW', rvar1=sea%sea_wthv)
     if (allocated(sea%ice_wthv))       call increment_vtable('SEA%ICE_WTHV',       'SW', rvar1=sea%ice_wthv)
     if (allocated(sea%windxe))         call increment_vtable('SEA%WINDXE',         'SW', rvar1=sea%windxe)
     if (allocated(sea%windye))         call increment_vtable('SEA%WINDYE',         'SW', rvar1=sea%windye)
     if (allocated(sea%windze))         call increment_vtable('SEA%WINDZE',         'SW', rvar1=sea%windze)
     if (allocated(sea%sea_albedo))     call increment_vtable('SEA%SEA_ALBEDO',     'SW', rvar1=sea%sea_albedo)
     if (allocated(sea%ice_albedo))     call increment_vtable('SEA%ICE_ALBEDO',     'SW', rvar1=sea%ice_albedo)
     if (allocated(sea%sea_rlongup))    call increment_vtable('SEA%SEA_RLONGUP',    'SW', rvar1=sea%sea_rlongup)
     if (allocated(sea%ice_rlongup))    call increment_vtable('SEA%ICE_RLONGUP',    'SW', rvar1=sea%ice_rlongup)
     if (allocated(sea%ice_net_rlong))  call increment_vtable('SEA%ICE_NET_RLONG',  'SW', rvar1=sea%ice_net_rlong)
     if (allocated(sea%ice_net_rshort)) call increment_vtable('SEA%ICE_NET_RSHORT', 'SW', rvar1=sea%ice_net_rshort)
     if (allocated(sea%seatc))          call increment_vtable('SEA%SEATC',          'SW', rvar1=sea%seatc)
     if (allocated(sea%sea_cantemp))    call increment_vtable('SEA%SEA_CANTEMP',    'SW', rvar1=sea%sea_cantemp)
     if (allocated(sea%ice_cantemp))    call increment_vtable('SEA%ICE_CANTEMP',    'SW', rvar1=sea%ice_cantemp)
     if (allocated(sea%sea_canrrv))     call increment_vtable('SEA%SEA_CANRRV',     'SW', rvar1=sea%sea_canrrv)
     if (allocated(sea%ice_canrrv))     call increment_vtable('SEA%ICE_CANRRV',     'SW', rvar1=sea%ice_canrrv)
     if (allocated(sea%nlev_seaice))    call increment_vtable('SEA%NLEV_SEAICE',    'SW', ivar1=sea%nlev_seaice)
     if (allocated(sea%seaicec))        call increment_vtable('SEA%SEAICEC',        'SW', rvar1=sea%seaicec)
     if (allocated(sea%surface_srrv))   call increment_vtable('SEA%SURFACE_SRRV',   'SW', rvar1=sea%surface_srrv)
     if (allocated(sea%sea_sfc_srrv))   call increment_vtable('SEA%SEA_SFC_SRRV',   'SW', rvar1=sea%sea_sfc_srrv)
     if (allocated(sea%ice_sfc_srrv))   call increment_vtable('SEA%ICE_SFC_SRRV',   'SW', rvar1=sea%ice_sfc_srrv)
     if (allocated(sea%sea_rough))      call increment_vtable('SEA%SEA_ROUGH',      'SW', rvar1=sea%sea_rough)
     if (allocated(sea%ice_rough))      call increment_vtable('SEA%ICE_ROUGH',      'SW', rvar1=sea%ice_rough)
     if (allocated(sea%seaice_energy))  call increment_vtable('SEA%SEAICE_ENERGY',  'SW', rvar2=sea%seaice_energy)
     if (allocated(sea%seaice_tempk))   call increment_vtable('SEA%SEAICE_TEMPK',   'SW', rvar2=sea%seaice_tempk)

     if (allocated(sea%swmdepth))       call increment_vtable('SEA%SWMDEPTH',       'SW', rvar1=sea%swmdepth)
     if (allocated(sea%vxe))            call increment_vtable('SEA%VXE',            'SW', rvar1=sea%vxe)
     if (allocated(sea%vye))            call increment_vtable('SEA%VYE',            'SW', rvar1=sea%vye)
     if (allocated(sea%vze))            call increment_vtable('SEA%VZE',            'SW', rvar1=sea%vze)

   end subroutine filltab_sea

End Module mem_sea
