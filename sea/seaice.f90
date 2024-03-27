subroutine prep_seaice( sst, seaice, sea_cantemp, ice_cantemp,        &
                        seaice_energy, seaice_tempk, nlev_seaice,     &
                        ice_albedo, ice_rlongup, rshort_i, rlong_i,   &
                        rshort, rlong, rough, sea_canrrv, ice_canrrv, &
                        sea_ustar, ice_ustar, sea_ggaer, ice_ggaer,   &
                        sea_wthv, ice_wthv                            )

  use consts_coms, only: cice, t00
  use sea_coms,    only: nzi, t00sea

  implicit none

  real,    intent(in)    :: sst
  real,    intent(in)    :: seaice
  real,    intent(in)    :: sea_cantemp
  real,    intent(out)   :: ice_cantemp
  real,    intent(out)   :: seaice_energy(nzi) ! seaice layer energy [J/kg]
  real,    intent(out)   :: seaice_tempk (nzi) ! seaice layer temperature [K]
  integer, intent(inout) :: nlev_seaice
  real,    intent(out)   :: ice_albedo
  real,    intent(out)   :: ice_rlongup
  real,    intent(out)   :: rshort_i
  real,    intent(out)   :: rlong_i
  real,    intent(in)    :: rshort
  real,    intent(in)    :: rlong
  real,    intent(out)   :: rough  ! sea cell roughess height [m]
  real,    intent(in)    :: sea_canrrv
  real,    intent(out)   :: ice_canrrv
  real,    intent(in)    :: sea_ustar
  real,    intent(out)   :: ice_ustar
  real,    intent(in)    :: sea_ggaer
  real,    intent(out)   :: ice_ggaer
  real,    intent(in)    :: sea_wthv
  real,    intent(out)   :: ice_wthv
  integer                :: k

! First check whether sea ice is present

  if (seaice < 0.01) then

     ! If no sea ice present, nullify quantities and return

     nlev_seaice      = 0
     seaice_tempk (:) = 0.0
     seaice_energy(:) = 0.0
     ice_cantemp      = 0.0
     ice_canrrv       = 0.0
     ice_albedo       = 0.0
     ice_rlongup      = 0.0
     rshort_i         = 0.0
     rlong_i          = 0.0
     rough            = 0.0
     ice_ustar        = 0.0
     ice_ggaer        = 0.0
     ice_wthv         = 0.0

     return

  else

     ! Always set bottom layer temperature equal to the SST

     seaice_tempk (1) = min(sst, t00sea)
     seaice_energy(1) = cice * (seaice_tempk(1) - t00sea)

     ! Ice roughness height

     rough = 5.0e-4 ! [m]; taken from Los Alamos sea ice model

  endif

! Seaice is present. If there was no seaice the previous timestep,
! initialize the ice temperature and energy based on the SST and air
! temperature, and initialize the ice canopy temperature, vapor, and
! ustar to the sea canopy values

  if (nlev_seaice == 0) then

     nlev_seaice = nzi
     ice_cantemp = sea_cantemp
     ice_canrrv  = sea_canrrv
     ice_ustar   = sea_ustar
     ice_ggaer   = sea_ggaer
     ice_wthv    = sea_wthv

     ! Initialize top layer to the canopy temperature

     seaice_tempk (nlev_seaice) = min(ice_cantemp, t00sea)
     seaice_energy(nlev_seaice) = cice * (seaice_tempk(nlev_seaice) - t00sea)

     ! Interpolate other layers between the top and bottom temperature

     do k = 2, nlev_seaice-1

        seaice_tempk(k) = real(k-1)/real(nlev_seaice-1) * seaice_tempk(nlev_seaice) &
                        + real(nlev_seaice-k)/real(nlev_seaice-1) * seaice_tempk(1)

        seaice_energy(k) = cice * (seaice_tempk(k) - t00sea)

     enddo

     ! Estimate ice albedo and outgoing longwave in case we have
     ! created seaice between calls to the radiation scheme

     call sfcrad_seaice_1( ice_rlongup,        &
                           ice_albedo,         &
                           nlev_seaice,        &
                           ice_cantemp,        &
                           seaice_tempk(1:nzi) )

     ! Estimate net seaice longwave and shortwave radiative fluxes in case
     ! seaice has been created between calls to the radiation scheme

     call sfcrad_seaice_2( rshort_i,    &
                           rlong_i,     &
                           nlev_seaice, &
                           rshort,      &
                           rlong,       &
                           ice_rlongup, &
                           ice_albedo   )

  endif

end subroutine prep_seaice

!===============================================================================

subroutine seaice( isea, iwsfc, nlev_seaice, rhos, ustar, vkhsfc, can_depth, &
                   rshort_i, rlong_i, glatw, glonw, airtheta, airrrv, canexner, &
                   seaice_energy, seaice_tempk, cantemp, canrrv, sfluxt, sfluxr, &
                   surface_srrv )

  use consts_coms, only: alvi, cice, t00, cp, alli, r8
  use sea_coms,    only: dt_sea, t00sea, nzi
  use therm_lib,   only: rhovsi, qtk_sea
  use mem_sfcg,    only: sfcg
  use matrix,      only: matrix8_NxN
  use leaf4_canopy,only: sing_print

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
  integer, intent(in)    :: iwsfc        ! current sfc grid cell index
  integer, intent(in)    :: nlev_seaice  ! number of seaice levels
  real,    intent(in)    :: rhos         ! air density [kg/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: vkhsfc       ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: rshort_i     ! s/w net rad flux to seaice [W/m^2]
  real,    intent(in)    :: rlong_i      ! l/w net rad flux to seaice [W/m^2]
  real,    intent(in)    :: glatw        ! Latitude of lake cell 'center' [deg]
  real,    intent(in)    :: glonw        ! Longitude of lake cell 'center' [deg]
  real,    intent(in)    :: airtheta     ! atm potential temp [K]
  real,    intent(in)    :: airrrv       ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: canexner     ! canopy Exner function []
  real,    intent(inout) :: seaice_energy(nzi) ! seaice layer energy [J/kg]
  real,    intent(inout) :: seaice_tempk (nzi) ! seaice layer temperature [K]
  real,    intent(inout) :: cantemp      ! ice "canopy" air temp [K]
  real,    intent(inout) :: canrrv       ! ice "canopy" air vapor mix ratio [kg_vap/kg_dryair]
  real,    intent(out)   :: sfluxt       ! can_air to atm heat flux [W m^-2]
  real,    intent(out)   :: sfluxr       ! can_air to atm vapor flux [kg_vap m^-2 s^-1]
  real,    intent(out)   :: surface_srrv ! ice surface sat mix ratio [kg_vap/kg_dryair]

! Local parameters

  real, parameter :: fcn = 0.75         ! Crank-Nicolson future time weight
  real, parameter :: iceden   = 900.0   ! Density of seaice layers [kg/m^3]
  real, parameter :: icethick = 0.5     ! Thickness of seaice layers [m]
  real, parameter :: icek     = 2.0     ! Thermal conductivity of seaice [W/(K m)]

  real, parameter :: icemass  = icethick * iceden ! mass of unit area of ice [kg/m^2]
  real, parameter :: icer     = icethick / icek   ! ice thermal resistivity [K m^2/W]

! Local variables

  real :: energy_per_m2(nzi) ! seaice energy [J/m^2]
  real :: hxfers(nzi+1)      ! seaice heat xfer [J/m2]

  real :: fracliq   ! fraction of water in liquid phase
  real :: rdi       ! canopy conductance [m/s]
  real :: hxferic   ! seaice-to-can_air heat xfer this step [J/m^2]
  real :: wxferic   ! seaice-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxferca   ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca   ! vapor xfer from can_air to atm this step [kg/m^2]
  real :: sfc_rhovs ! sat vapor density at seaice surface temp [kg_vap/m^3]
  real :: can_rhov  ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair    ! Canopy air mass [kg/m^2]
  real :: hcapcan   ! Canopy air heat capacity [J/(m^2 K)]
  real :: canairi   ! Inverse of canair
  real :: hcapcani  ! Inverse of hcapcan

  integer :: k

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
  real(r8) :: h1, h2, h4, h5, h6, h7, h8
  real(r8) :: y1, y2, y3, y4, y5, y6, y7, y8, y9, y10

  real(r8) :: aa4(4,4), xx4(4), yy4(4)         ! 4x4 matrix equation terms

  logical :: sing

  if (nlev_seaice == 0) return

  ! Compute sfcwater energy per m^2

  do k = 1, nlev_seaice
     energy_per_m2(k) = seaice_energy(k) * icemass
  enddo

! Evaluate surface saturation vapor density and mixing ratio of sea surface

  sfc_rhovs    = rhovsi(seaice_tempk(nlev_seaice)-t00)
  surface_srrv = sfc_rhovs / rhos

  ! rdi = ustar/5 is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  ! Canopy air quantities

  can_rhov = canrrv * rhos
  canair   = rhos * can_depth
  canairi  = 1. / canair
  hcapcan  = cp * canair
  hcapcani = 1. / hcapcan

  ! Set up and solve a 4x4 linear system of equations that use trapezoidal-implicit
  ! differencing.  The solution of the system consists of turbulent heat and
  ! water vapor fluxes between canopy air, the seaice surface, and the free
  ! atmosphere, and the consequent changes to water and temperature of canopy
  ! air and the seaice surface.

  ! It is assumed that the heat capacity of the (upper level of) seaice
  ! is sufficiently high that seaice temperature change is negligible in the
  ! context of the matrix solver.

  a5  = dt_sea * rdi   ! sfc sublayer vap xfer coef
  a6  = cp * rhos * a5 ! sfc sublayer heat xfer coef
  a9  = dt_sea * vkhsfc / sfcg%dzt_bot(iwsfc)
  a10 = cp * a9

  h4 = fcn * rhos * canairi  ! = fcn / can_depth
  h7 = fcn * hcapcani
  h8 = fcn * canairi

  y2  = sfc_rhovs - can_rhov
  y5  = seaice_tempk(nlev_seaice) - cantemp
  y9  = canrrv  - airrrv
  y10 = cantemp - canexner * airtheta

  aa4(1,1) = 1._r8 + a5 * h4
  aa4(1,2) = 0._r8
  aa4(1,3) =       - a5 * h4
  aa4(1,4) = 0._r8
  yy4(1)   =         a5 * y2   ! WSC row

  aa4(2,1) = 0._r8
  aa4(2,2) = 1._r8 + a6 * h7
  aa4(2,3) = 0._r8
  aa4(2,4) =       - a6 * h7
  yy4(2)   =         a6 * y5   ! HSC row

  aa4(3,1) =       - a9 * h8
  aa4(3,2) = 0._r8
  aa4(3,3) = 1._r8 + a9 * h8
  aa4(3,4) = 0._r8
  yy4(3)   =         a9 * y9   ! WCA row

  aa4(4,1) = 0._r8
  aa4(4,2) =       - a10 * h7
  aa4(4,3) = 0._r8
  aa4(4,4) = 1._r8 + a10 * h7
  yy4(4)   =         a10 * y10 ! HCA row

  call matrix8_NxN(4,aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'seaice1',4,aa4,yy4,glatw,glonw)

  wxferic = xx4(1)
  hxferic = xx4(2)
  wxferca = xx4(3)
  hxferca = xx4(4)

  cantemp = cantemp + (hxferic - hxferca) * hcapcani
  canrrv  = canrrv  + (wxferic - wxferca) * canairi

  sfluxt = hxferca / dt_sea
  sfluxr = wxferca / dt_sea

! Zero out sfcwater internal heat transfer array at top and bottom surfaces.
! Energy transfer at top is applied separately.

  hxfers(1)             = 0.0
  hxfers(nlev_seaice+1) = 0.0

! Compute internal seaice energy xfer if at least two layers exist [J/m2]

  do k = 2, nlev_seaice
     hxfers(k) = dt_sea * (seaice_tempk(k-1) - seaice_tempk(k)) / icer
  enddo

! Add contributions to seaice energy from internal transfers of heat

  do k = 2, nlev_seaice
     energy_per_m2(k) = energy_per_m2(k) + hxfers(k) - hxfers(k+1)
  enddo

! Apply long- and short-wave radiative transfer and seaice-to-canopy
! sensible and latent heat transfer to top seaice layer

  energy_per_m2(nlev_seaice) = energy_per_m2(nlev_seaice)  &
       + dt_sea * (rshort_i + rlong_i) - hxferic - wxferic * alvi

! Compute new seaice_energy

  do k = 2, nlev_seaice

     seaice_energy(k) = energy_per_m2(k) / icemass

     ! Bound seaice_energy so that the fraction of liquid water in each layer
     ! is never larger than 50%. For now, we trust that ice should exist when
     ! the input seaice data indicates that ice is present, and this keeps the
     ! seaice temperature at or below its freezing point

     seaice_energy(k) = min( seaice_energy(k), 0.5 * alli )

  enddo

! Update seaice temperature

  do k = 2, nlev_seaice
     call qtk_sea( seaice_energy(k), seaice_tempk(k), fracliq )
  enddo

end subroutine seaice

!===============================================================================

subroutine sfcrad_seaice_1( rlongup, albedo, nlev_seaice, cantemp, seaice_tempk )

  use sea_coms,    only: nzi, emi
  use consts_coms, only: stefan
  implicit none

  real,    intent(out) :: rlongup   ! seaice outgoing L/W radiation[W/m^2]
  real,    intent(out) :: albedo    ! seaice albedo                [0-1]

  integer, intent(in)  :: nlev_seaice       ! number of seaice levels
  real,    intent(in)  :: cantemp           ! seaice canopy temperature [K]
  real,    intent(in)  :: seaice_tempk(nzi) ! seaice layer temperatures [K]

  if (nlev_seaice > 0) then

     rlongup = emi * stefan * seaice_tempk(nlev_seaice) ** 4

     ! In the Los Alamos sea ice model, visible albedo is 0.78 and
     ! NIR albedo is 0.36 for temperatures below -1C.  Here, we assume
     ! equal parts visible and NIR, yielding an albedo of 0.57.  This
     ! albedo decreases to 0.5 as the temperature increases to 0C.

     if (cantemp < 272.15) then
        albedo = 0.72  ! Dave modification
     !  albedo = 0.57
     else
        albedo = 0.72 - 0.22 * ( min( 273.15, cantemp) - 272.15 )
     endif

  else

     rlongup = 0.0
     albedo  = 0.0

  end if

end subroutine sfcrad_seaice_1

!===============================================================================

subroutine sfcrad_seaice_2( rshort_i, rlong_i, nlev_seaice, &
                            rshort, rlong, rlongup, albedo  )

  use sea_coms, only: nzi, emi
  implicit none

  real,    intent(out) :: rshort_i    ! net S/W radiation at ice surface [W/m^2]
  real,    intent(out) :: rlong_i     ! net L/W radiation at ice surface [W/m^2]

  integer, intent(in)  :: nlev_seaice ! number of seaice levls
  real,    intent(in)  :: rshort      ! incoming S/W radiation       [W/m^2]
  real,    intent(in)  :: rlong       ! incoming L/W radiation       [W/m^2]
  real,    intent(in)  :: rlongup     ! seaice outgoing L/W radiation[W/m^2]
  real,    intent(in)  :: albedo      ! seaice albedo                [0-1]

  if (nlev_seaice > 0) then

     rshort_i = (1.0 - albedo) * rshort
     rlong_i  = rlong * emi - rlongup

  else

     rshort_i = 0.0
     rlong_i  = 0.0

  end if

end subroutine sfcrad_seaice_2
