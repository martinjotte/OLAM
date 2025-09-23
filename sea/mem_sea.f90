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
      integer :: iwglobe
   End type itab_sea_vars

   type (itab_sea_vars), allocatable, target :: itab_sea(:)

! SEA MODEL VARIABLES

   Type sea_vars

! Canopy to atmosphere turbulent flux quantities

      real, allocatable :: sea_ustar  (:) ! friction velocity [m/s]
      real, allocatable :: ice_ustar  (:) ! friction velocity [m/s]

      real, allocatable :: sea_vkmsfc (:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: ice_vkmsfc (:) ! surface drag coefficient [kg/(m s)]

      real, allocatable :: sea_vkhsfc (:) ! surface heat and vapor transfer coefficient [kg/(m s)]
      real, allocatable :: ice_vkhsfc (:) ! surface heat and vapor transfer coefficient [kg/(m s)]

      real, allocatable :: sea_sfluxt (:) ! can_air-to-atm sensible heat flux over sea water[W m^-2]
      real, allocatable :: ice_sfluxt (:) ! can_air-to-atm sensible heat flux over sea ice [W m^-2]

      real, allocatable :: sea_sfluxr (:) ! can_air-to-atm water vapor flux over sea water[kg_vap m^-2 s^-1]
      real, allocatable :: ice_sfluxr (:) ! can_air-to-atm water vapor flux over sea ice [kg_vap m^-2 s^-1]

      real, allocatable :: sea_sfluxc (:)
      real, allocatable :: ice_sfluxc (:)

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

      real, allocatable :: spraytemp     (:) ! seaspray temperature [K]
      real, allocatable :: spray2temp    (:) ! seaspray2 temperature [K]

      integer, allocatable :: nlev_seaice(:) ! number of seaice layers in current cell

      real, allocatable :: seatp         (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf         (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc         (:) ! current sea temperature [K]
      real, allocatable :: seaicep       (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef       (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec       (:) ! seaice fraction [0-1]

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
      integer, allocatable :: pom_kba   (:) ! number of POM1D levels

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

  subroutine alloc_sea1(msea0)

    use misc_coms, only: rinit, runtype
    implicit none

    integer, intent(in) :: msea0
    integer             :: isea

    ! This routine allocates sea arrays necessary for the MAKEGRID stage
    ! or for reading sfcgrid information

    allocate( itab_sea      (msea0) )
    allocate( sea%pom_active(msea0) )

    if (npomzons > 0 .or. runtype == 'MAKEGRID') then
       allocate( sea%pom_kba(msea0) )
    endif

    !$omp parallel do
    do isea = 1, msea0

       if ( allocated( itab_sea       ) ) itab_sea      (isea) = itab_sea_vars( iwglobe = 1 )
       if ( allocated( sea%pom_active ) ) sea%pom_active(isea) = .false.
       if ( allocated( sea%pom_kba    ) ) sea%pom_kba   (isea) = 0

    enddo
    !$omp end parallel do

  end subroutine alloc_sea1

!=========================================================================

  subroutine alloc_sea2()

     use misc_coms,  only: rinit
     use sea_coms,   only: nzi
     use mem_co2,    only: co2flag
     use mem_sfcg,   only: nswmzons
     use oname_coms, only: nl

     implicit none

     integer :: isea

     ! This routine allocates and initializes sea arrays needed for a full
     ! model integration that weren't allocated in alloc_sea1

     allocate (sea%sea_ustar     (msea))
     allocate (sea%ice_ustar     (msea))

     allocate (sea%sea_vkmsfc    (msea))
     allocate (sea%ice_vkmsfc    (msea))

     allocate (sea%sea_vkhsfc    (msea))
     allocate (sea%ice_vkhsfc    (msea))

     allocate (sea%sea_sfluxt    (msea))
     allocate (sea%ice_sfluxt    (msea))

     allocate (sea%sea_sfluxr    (msea))
     allocate (sea%ice_sfluxr    (msea))

     if (co2flag /= 0) then
        allocate (sea%sea_sfluxc (msea))
        allocate (sea%ice_sfluxc (msea))
     endif

     allocate (sea%sea_ggaer     (msea))
     allocate (sea%ice_ggaer     (msea))

     allocate (sea%sea_wthv      (msea))
     allocate (sea%ice_wthv      (msea))

     allocate (sea%windxe        (msea))
     allocate (sea%windye        (msea))
     allocate (sea%windze        (msea))

     allocate (sea%sea_albedo    (msea))
     allocate (sea%ice_albedo    (msea))

     allocate (sea%sea_rlongup   (msea))
     allocate (sea%ice_rlongup   (msea))

     allocate (sea%ice_net_rlong (msea))
     allocate (sea%ice_net_rshort(msea))

     allocate (sea%sea_cantemp   (msea))
     allocate (sea%ice_cantemp   (msea))

     allocate (sea%sea_canrrv    (msea))
     allocate (sea%ice_canrrv    (msea))

     if (nl%iseasprayflg > 0) &
          allocate (sea%spraytemp(msea))

     if (nl%iseasprayflg > 1) &
          allocate (sea%spray2temp(msea))

     allocate (sea%nlev_seaice   (msea))

     allocate (sea%seatp         (msea))
     allocate (sea%seatf         (msea))
     allocate (sea%seatc         (msea))
     allocate (sea%seaicep       (msea))
     allocate (sea%seaicef       (msea))
     allocate (sea%seaicec       (msea))

     allocate (sea%sea_rough     (msea))
     allocate (sea%ice_rough     (msea))

     allocate (sea%seaice_energy(nzi,msea))
     allocate (sea%seaice_tempk (nzi,msea))

     if (nswmzons > 0) then

        allocate (sea%swmdepth   (msea))

        allocate (sea%vxe        (msea))
        allocate (sea%vye        (msea))
        allocate (sea%vze        (msea))

        allocate (sea%vmxet      (msea))
        allocate (sea%vmyet      (msea))
        allocate (sea%vmzet      (msea))

        allocate (sea%vmxet_area (msea))
        allocate (sea%vmyet_area (msea))
        allocate (sea%vmzet_area (msea))

        allocate (sea%vxe1       (msea))
        allocate (sea%vye1       (msea))
        allocate (sea%vze1       (msea))

        allocate (sea%gxps_vxe   (msea))
        allocate (sea%gyps_vxe   (msea))

        allocate (sea%gxps_vye   (msea))
        allocate (sea%gyps_vye   (msea))

        allocate (sea%gxps_vze   (msea))
        allocate (sea%gyps_vze   (msea))

     endif

     !$omp parallel do
     do isea = 1, msea

        if ( allocated( sea%sea_ustar     ) ) sea%sea_ustar        (isea) = rinit
        if ( allocated( sea%ice_ustar     ) ) sea%ice_ustar        (isea) = rinit
        if ( allocated( sea%sea_vkmsfc    ) ) sea%sea_vkmsfc       (isea) = rinit
        if ( allocated( sea%ice_vkmsfc    ) ) sea%ice_vkmsfc       (isea) = rinit
        if ( allocated( sea%sea_vkhsfc    ) ) sea%sea_vkhsfc       (isea) = rinit
        if ( allocated( sea%ice_vkhsfc    ) ) sea%ice_vkhsfc       (isea) = rinit
        if ( allocated( sea%sea_sfluxt    ) ) sea%sea_sfluxt       (isea) = rinit
        if ( allocated( sea%ice_sfluxt    ) ) sea%ice_sfluxt       (isea) = rinit
        if ( allocated( sea%sea_sfluxr    ) ) sea%sea_sfluxr       (isea) = rinit
        if ( allocated( sea%ice_sfluxr    ) ) sea%ice_sfluxr       (isea) = rinit
        if ( allocated( sea%sea_sfluxc    ) ) sea%sea_sfluxc       (isea) = rinit
        if ( allocated( sea%ice_sfluxc    ) ) sea%ice_sfluxc       (isea) = rinit
        if ( allocated( sea%sea_ggaer     ) ) sea%sea_ggaer        (isea) = rinit
        if ( allocated( sea%ice_ggaer     ) ) sea%ice_ggaer        (isea) = rinit
        if ( allocated( sea%sea_wthv      ) ) sea%sea_wthv         (isea) = rinit
        if ( allocated( sea%ice_wthv      ) ) sea%ice_wthv         (isea) = rinit
        if ( allocated( sea%windxe        ) ) sea%windxe           (isea) = rinit
        if ( allocated( sea%windye        ) ) sea%windye           (isea) = rinit
        if ( allocated( sea%windze        ) ) sea%windze           (isea) = rinit
        if ( allocated( sea%sea_albedo    ) ) sea%sea_albedo       (isea) = 0.0
        if ( allocated( sea%ice_albedo    ) ) sea%ice_albedo       (isea) = 0.0
        if ( allocated( sea%sea_rlongup   ) ) sea%sea_rlongup      (isea) = 0.0
        if ( allocated( sea%ice_rlongup   ) ) sea%ice_rlongup      (isea) = 0.0
        if ( allocated( sea%ice_net_rlong ) ) sea%ice_net_rlong    (isea) = 0.0
        if ( allocated( sea%ice_net_rshort) ) sea%ice_net_rshort   (isea) = 0.0
        if ( allocated( sea%sea_cantemp   ) ) sea%sea_cantemp      (isea) = rinit
        if ( allocated( sea%ice_cantemp   ) ) sea%ice_cantemp      (isea) = rinit
        if ( allocated( sea%sea_canrrv    ) ) sea%sea_canrrv       (isea) = rinit
        if ( allocated( sea%ice_canrrv    ) ) sea%ice_canrrv       (isea) = rinit
        if ( allocated( sea%spraytemp     ) ) sea%spraytemp        (isea) = rinit
        if ( allocated( sea%spray2temp    ) ) sea%spray2temp       (isea) = rinit
        if ( allocated( sea%nlev_seaice   ) ) sea%nlev_seaice      (isea) = 0
        if ( allocated( sea%seatp         ) ) sea%seatp            (isea) = rinit
        if ( allocated( sea%seatf         ) ) sea%seatf            (isea) = rinit
        if ( allocated( sea%seatc         ) ) sea%seatc            (isea) = rinit
        if ( allocated( sea%seaicep       ) ) sea%seaicep          (isea) = rinit
        if ( allocated( sea%seaicef       ) ) sea%seaicef          (isea) = rinit
        if ( allocated( sea%seaicec       ) ) sea%seaicec          (isea) = rinit
        if ( allocated( sea%sea_rough     ) ) sea%sea_rough        (isea) = rinit
        if ( allocated( sea%ice_rough     ) ) sea%ice_rough        (isea) = rinit
        if ( allocated( sea%seaice_energy ) ) sea%seaice_energy  (:,isea) = rinit
        if ( allocated( sea%seaice_tempk  ) ) sea%seaice_tempk   (:,isea) = rinit
        if ( allocated( sea%swmdepth      ) ) sea%swmdepth         (isea) = 0.0
        if ( allocated( sea%vxe           ) ) sea%vxe              (isea) = 0.0
        if ( allocated( sea%vye           ) ) sea%vye              (isea) = 0.0
        if ( allocated( sea%vze           ) ) sea%vze              (isea) = 0.0
        if ( allocated( sea%vmxet         ) ) sea%vmxet            (isea) = 0.0
        if ( allocated( sea%vmyet         ) ) sea%vmyet            (isea) = 0.0
        if ( allocated( sea%vmzet         ) ) sea%vmzet            (isea) = 0.0
        if ( allocated( sea%vmxet_area    ) ) sea%vmxet_area       (isea) = 0.0
        if ( allocated( sea%vmyet_area    ) ) sea%vmyet_area       (isea) = 0.0
        if ( allocated( sea%vmzet_area    ) ) sea%vmzet_area       (isea) = 0.0
        if ( allocated( sea%vxe1          ) ) sea%vxe1             (isea) = 0.0
        if ( allocated( sea%vye1          ) ) sea%vye1             (isea) = 0.0
        if ( allocated( sea%vze1          ) ) sea%vze1             (isea) = 0.0
        if ( allocated( sea%gxps_vxe      ) ) sea%gxps_vxe         (isea) = 0.0
        if ( allocated( sea%gyps_vxe      ) ) sea%gyps_vxe         (isea) = 0.0
        if ( allocated( sea%gxps_vye      ) ) sea%gxps_vye         (isea) = 0.0
        if ( allocated( sea%gyps_vye      ) ) sea%gyps_vye         (isea) = 0.0
        if ( allocated( sea%gxps_vze      ) ) sea%gxps_vze         (isea) = 0.0
        if ( allocated( sea%gyps_vze      ) ) sea%gyps_vze         (isea) = 0.0

     enddo
     !$omp end parallel do

  end subroutine alloc_sea2

!=========================================================================

  subroutine filltab_sea()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(sea%sea_ustar))      call increment_vtable('SEA%SEA_USTAR',      'SW', rvar1=sea%sea_ustar)
     if (allocated(sea%ice_ustar))      call increment_vtable('SEA%ICE_USTAR',      'SW', rvar1=sea%ice_ustar)
     if (allocated(sea%sea_vkmsfc))     call increment_vtable('SEA%SEA_VKMSFC',     'SW', rvar1=sea%sea_vkmsfc)
     if (allocated(sea%ice_vkmsfc))     call increment_vtable('SEA%ICE_VKMSFC',     'SW', rvar1=sea%ice_vkmsfc)
     if (allocated(sea%sea_vkhsfc))     call increment_vtable('SEA%SEA_VKHSFC',     'SW', rvar1=sea%sea_vkhsfc)
     if (allocated(sea%ice_vkhsfc))     call increment_vtable('SEA%ICE_VKHSFC',     'SW', rvar1=sea%ice_vkhsfc)
     if (allocated(sea%sea_sfluxt))     call increment_vtable('SEA%SEA_SFLUXT',     'SW', rvar1=sea%sea_sfluxt)
     if (allocated(sea%ice_sfluxt))     call increment_vtable('SEA%ICE_SFLUXT',     'SW', rvar1=sea%ice_sfluxt)
     if (allocated(sea%sea_sfluxr))     call increment_vtable('SEA%SEA_SFLUXR',     'SW', rvar1=sea%sea_sfluxr)
     if (allocated(sea%ice_sfluxr))     call increment_vtable('SEA%ICE_SFLUXR',     'SW', rvar1=sea%ice_sfluxr)
     if (allocated(sea%sea_sfluxc))     call increment_vtable('SEA%SEA_SFLUXC',     'SW', rvar1=sea%sea_sfluxc)
     if (allocated(sea%ice_sfluxc))     call increment_vtable('SEA%ICE_SFLUXC',     'SW', rvar1=sea%ice_sfluxc)
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
     if (allocated(sea%spraytemp))      call increment_vtable('SEA%SPRAYTEMP',      'SW', rvar1=sea%spraytemp)
     if (allocated(sea%spray2temp))     call increment_vtable('SEA%SPRAY2TEMP',     'SW', rvar1=sea%spray2temp)
     if (allocated(sea%nlev_seaice))    call increment_vtable('SEA%NLEV_SEAICE',    'SW', ivar1=sea%nlev_seaice)
     if (allocated(sea%seaicec))        call increment_vtable('SEA%SEAICEC',        'SW', rvar1=sea%seaicec)
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
