Module leaf4_canopy

Contains

  subroutine canopy(iland, iwsfc, nlsw1, icomb, ktrans, transp,                  &
                    leaf_class,    can_depth,       rhos,            vels,       &
                    ustar,         vkhsfc,          sfluxt,          sfluxr,     &
                    pcpg,          qpcpg,           dpcpg,           rshort,     &
                    cantemp,       canrrv,          glatw,           glonw,      &
                    airtheta,      airrrv,          canexner,                    &
                    snowfac,       vf,              stom_resist,     veg_height, &
                    veg_rough,     veg_tai,         veg_lai,         hcapveg,    &
                    rshort_v,      rshort_g,        rlong_v,         rlong_s,    &
                    rlong_g,       veg_water,       veg_energy,      veg_temp,   &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,              &
                    energy_per_m2, sfcwater_tempk,  sfcwater_fracliq,soil_wfrac, &
                    soil_water,    soil_energy,     wsat_vg,         wresid_vg,  &
                    soilfldcap,    ksat_vg,         specifheat_drysoil,          &
                    head,          head_slope,      soil_tempk,      soil_fracliq)

  use leaf_coms,     only: soil_rough, dt_leaf, kroot, rcmin, snowmin_expl, &
                           wcap_vmin
  use mem_sfcg,      only: sfcg
  use mem_land,      only: nzg, dslz, dslzi
  use leaf4_surface, only: sfcwater_soil_comb, grndvap_ab
  use consts_coms,   only: cp, vonk, alvi, alvl, alli, cliq, cice, rvap, t00, r8
  use matrix,        only: matrix8_2x2, matrix8_3x3, matrix8_4x4, matrix8_NxN
  use therm_lib,     only: qwtk, rhovsil, eslf

  implicit none

  integer, intent(in)    :: iland      ! index of current land cell
  integer, intent(in)    :: iwsfc      ! index of current SFC grid cell
  integer, intent(in)    :: nlsw1      ! k index of top active sfcwater layer
  integer, intent(inout) :: icomb      ! implicit heat balance flag [0=no, 1=yes]
  integer, intent(out)   :: ktrans     ! k index of soil layer supplying transp
  real,    intent(out)   :: transp     ! transpiration xfer this LEAF timestep [kg_vap/m^2]
  integer, intent(in)    :: leaf_class ! leaf class (vegetation class)

  real, intent(in)    :: can_depth        ! canopy depth for heat and vap capacity [m]
  real, intent(in)    :: rhos             ! atmospheric air density [kg_dryair/m^3]
  real, intent(in)    :: vels             ! surface wind speed [m/s]
  real, intent(in)    :: ustar            ! friction velocity [m/s]
  real, intent(in)    :: vkhsfc           ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real, intent(out)   :: sfluxt           ! can_air to atm heat flux [W m^-2]
  real, intent(out)   :: sfluxr           ! can_air to atm vapor flux [kg_vap m^-2 s^-1]
  real, intent(in)    :: pcpg             ! added precip mass this leaf timestep [kg/m^2]
  real, intent(in)    :: qpcpg            ! added precip energy this leaf timestep [J/m^2]
  real, intent(in)    :: dpcpg            ! added precip depth this leaf timestep [m]
  real, intent(in)    :: rshort           ! downward sfc s/w rad flux [W/m^2]
  real, intent(inout) :: cantemp          ! canopy air temp [K]
  real, intent(inout) :: canrrv           ! canopy air vapor mixing ratio [kg_vap/kg_dryair]
  real, intent(in)    :: glatw            ! Latitude of land cell 'center' [deg]
  real, intent(in)    :: glonw            ! Longitude of land cell 'center' [deg]
  real, intent(in)    :: airtheta         ! atm potential temp [K]
  real, intent(in)    :: airrrv           ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real, intent(in)    :: canexner         ! canopy Exner function []
  real, intent(in)    :: snowfac          ! fractional veg burial by snowcover
  real, intent(in)    :: vf               ! fractional coverage of non-buried part of veg
  real, intent(inout) :: stom_resist      ! veg stomatal resistance [s/m]
  real, intent(in)    :: veg_height       ! veg height [m]
  real, intent(in)    :: veg_rough        ! veg roughess height [m]
  real, intent(in)    :: veg_tai          ! veg total area index
  real, intent(in)    :: veg_lai          ! veg leaf area index
  real, intent(in)    :: hcapveg          ! veg heat capacity [J/(m^2 K)]
  real, intent(in)    :: rshort_v         ! s/w rad flux absorbed by veg [W/m^2]
  real, intent(in)    :: rshort_g         ! s/w rad flux absorbed by gnd [W/m^2]
  real, intent(in)    :: rlong_v          ! l/w rad flux absorbed by veg [W/m^2]
  real, intent(in)    :: rlong_s          ! l/w rad flux absorbed by sfc [W/m^2]
  real, intent(in)    :: rlong_g          ! l/w rad flux absorbed by gnd [W/m^2]
  real, intent(inout) :: veg_water        ! veg sfc water content [kg/m^2]
  real, intent(inout) :: veg_energy       ! (veg + veg_water) energy [J/m^2]
  real, intent(inout) :: veg_temp         ! veg temp [K]
  real, intent(inout) :: sfcwater_mass    ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_energy  ! surface water energy [J/kg]
  real, intent(inout) :: sfcwater_depth   ! surface water depth [m]
  real, intent(inout) :: energy_per_m2    ! surface water energy [J/m^2]
  real, intent(inout) :: sfcwater_tempk   ! surface water temperature [K]
  real, intent(inout) :: sfcwater_fracliq ! fraction of sfc water in liquid phase
  real, intent(inout) :: soil_wfrac        (nzg) ! soil water fractionn []
  real, intent(inout) :: soil_water        (nzg) ! soil water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy       (nzg) ! soil energy [J/m^3]
  real, intent(in)    :: wsat_vg           (nzg) ! saturation water content (porosity) []
  real, intent(in)    :: wresid_vg         (nzg) ! residual water content []
  real, intent(in)    :: soilfldcap              ! top-layer soil field capacity []
  real, intent(in)    :: ksat_vg           (nzg) ! saturation hydraulic conductivity [m/s]
  real, intent(in)    :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]
  real, intent(inout) :: head              (nzg) ! hydraulic head [m] (relative to local topo datum)
  real, intent(in)    :: head_slope        (nzg) ! d(head) / d(soil_water) [m]
  real, intent(in)    :: soil_tempk        (nzg) ! soil temp [K]
  real, intent(in)    :: soil_fracliq      (nzg) ! fraction of soil moisture in liquid phase

  ! Local parameters

  real, parameter :: exar = 2.5  ! for computing rasveg
  real, parameter :: covr = 2.16 ! scaling tai value for computing wtveg
  real, parameter :: c1 = 116.6  ! for computing rb

  ! Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar))
  ! from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
  ! 100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
  ! The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
  ! total expression inside the radical in Eq. 3.37.
  !bob: parameter(exar=3.5,covr=2.16,c1=98.8)

  ! Intercept, slope parameters for stomatal conductance factors

  real, parameter :: brad = 196.   , srad = .047    ! for s/w radiative flux
  real, parameter :: btlo = 281.5  , stlo = .26     ! for low canopy temperature
  real, parameter :: bthi = 310.1  , sthi = -.124   ! for high canopy temperature
  real, parameter :: bvpd = 4850.  , svpd = -.0051  ! for vapor pressure deficit
  real, parameter :: bswp = -1.07e6, sswp = 7.42e-6 ! for soil water potential

  real, parameter :: fcn = 0.5 ! Crank-Nicolson future time weight for canopy
                               ! turbulent flux balance

  ! Local variables

  integer :: k        ! loop index over soil layers
  integer :: iveg     ! flag for exposed vegetation (0=no, 1=yes)
  logical :: iwetveg  ! flag for wet vegetation

  real :: alpha       ! "alpha" term in soil surface humidity formulation
  real :: beta        ! "beta" term in soil surface humidity formulation
  real :: hxfersc     ! sfc-to-can_air heat xfer this step [J/m^2]
  real :: wxfersc     ! sfc-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxferca     ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca     ! vapor xfer from can_air to atm this step [kg_vap/m^2]
  real :: wxfergc     ! gnd-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxfervc     ! veg-to-can_air heat xfer this step [J/m^2]
  real :: wxfervc     ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
  real :: vegw_energy ! energy of veg_water [J/m^2]
  real :: wshed       ! water shed from veg this LEAF timestep [kg/m^2]
  real :: qshed       ! water energy shed from veg this LEAF timestep [J/m^2]
  real :: dshed       ! depth of water shed from veg [m]
  real :: factv       ! for computing rasveg
  real :: aux         ! for computing rasveg
  real :: rasveg      ! full-veg value of rd [s/m]
  real :: wtveg       ! weighting of rasveg in computing rd
! real :: rasgnd      ! not used
  real :: fracliqv    ! fraction of veg surface water in liquid phase
  real :: fthi        ! high-temp environ factor for stomatal resist
  real :: ftlo        ! low-temp environ factor for stomatal resist
  real :: frad        ! s/w radiative environ factor for stomatal resist
  real :: fswp        ! soil water potential environ factor for stomatal resist
  real :: fvpd        ! vap press deficit environ factor for stomatal resist
  real :: esat_veg    ! saturation vapor pressure at veg temp [Pa]
  real :: veg_rhovs   ! saturation vapor density at veg temp [kg_vap/m^3]
  real :: veg_rhovsp  ! saturation vapor density gradient at veg temp [kg_vap/(K m^3)]
  real :: e_can       ! vapor pressure of canopy air [Pa]
  real :: e_leaf      ! vapor pressure at leaf surface [Pa]
  real :: vpd         ! vapor pressure deficit across stomata [Pa]
  real :: rc_inf      ! asymptotic stomatal resistance at current veg environment [s/m]
  real :: sigmaw      ! fractional coverage of leaf surface by veg_water
  real :: slai        ! effective veg lai uncovered by surface water (snowcover)
  real :: stai        ! effective veg tai uncovered by surface water (snowcover)
  real :: swp         ! soil water potential factor for stomatal control
  real :: tvegc       ! vegetation temperature (deg C)
  real :: zognd       ! soil roughness height [m]
  real :: zoveg       ! vegetation roughness height [m]
  real :: zdisp       ! vegetation displacement height remaining after snowcover [m]
  real :: zveg        ! vegetation height remaining after snowcover [m]
  real :: transp_test ! test value of transpiration flux [kg_vap/(m^2 s)]
  real :: rdi         ! (soil or surface water)-to-can_air conductance [m/s]
  real :: rb          ! veg-to-can_air resistance [s/m]
  real :: rc          ! stomatal resistance [s/m]

  real :: can_rhov    ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair      ! Canopy air mass [kg/m^2]
  real :: hcapcan     ! Canopy air heat capacity [J/(m^2 K)]
  real :: hcapsfc     ! Heat capacity of top soil + any sfcwater, or only sfcwater [J/(m^2 K)]
  real :: hcapvegw    ! Heat capacity of vegetation + any veg_water [J/(m^2 K)]

  real :: canairi     ! Inverse of canair
  real :: hcapcani    ! Inverse of hcapcan
  real :: hcapsfci    ! Inverse of hcapsfc
  real :: hcapvegwi   ! Inverse of hcapvegw

  real :: alvveg      ! Latent heat of vaporization at veg surface, adjusted for phase of veg_water [J/(kg)]
  real :: alvsfc      ! Latent heat of vaporization at ground surface, adjusted for phase of sfcwater [J/(kg)]

  real :: radsfc      ! Radiation absorbed by surface [J/m^2]
  real :: radveg      ! Radiation absorbed by vegetation [J/m^2]

  real :: sfc_energy  ! temporary variable for surface energy [J/m^2]
  real :: sfc_wmass   ! temporary variable for surface water [kg/m^2]
  real :: sfc_soilhc  ! temporary variable for soil heat capacity [J/(m^2 K)]
  real :: sfc_tempk   ! temporary variable for surface temperature [K]
  real :: sfc_tempk1  ! temporary variable for sfc_tempk + 1 K

  real :: sfc_rhovs   ! surface saturation water vapor density [kg_vap/m^3]
  real :: sfc_rhovs1  ! sfc_rhovs at 1 K warmer [kg_vap/m^3]
  real :: sfc_rhovsp  ! sfc_rhovs derivative with respect to temperature

  real :: gnd_rhov    ! ground sfc evaporative water vapor density [kg_vap/m^3]
  real :: gnd_rhovp   ! gnd_rhov derivative with respect to temperature

  real :: eps_wxfer   ! very small water transfer to accomodate truncation error
  real :: sfcwater_fracarea
  real :: evap        ! amount evaporated from veg sfc [kg/m^2]

  real :: specvol, vw_max, fact
  real :: tf, wadd, wshed2, dshed2, delw

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
  real(r8) :: h1, h2, h3, h4, h5, h6, h7, h8
  real(r8) :: y1, y2, y3, y4, y5, y9, y10

  real(r8) :: aa3(3,3), xx3(3), yy3(3)  ! 3x3 matrix equation terms
  real(r8) :: aa4(4,4), xx4(4), yy4(4)  ! 4x4 matrix equation terms
  real(r8) :: aa5(5,5), xx5(5), yy5(5)  ! 5x5 matrix equation terms
  real(r8) :: aa6(6,6), xx6(6), yy6(6)  ! 6x6 matrix equation terms
  real(r8) :: aa7(7,7), xx7(7), yy7(7)  ! 7x7 matrix equation terms

  logical :: sing, skiptest(15)

  integer :: itest  ! test number

  integer, parameter :: zero8 = 0.0_r8

  eps_wxfer = dt_leaf * 1.e-10 ! 1.e-10 is 1 mm pcp/dew in 316 years  [kg/m^2]

  ! Set IVEG flag = 1 if there is sufficient exposed vegetation; set to 0 otherwise

  iveg = 0
  if (vf > .001) iveg = 1

  wshed = 0.
  qshed = 0.
  dshed = 0.

  ! Canopy air quantities

  can_rhov = canrrv * rhos
  canair   = rhos * can_depth
  canairi  = 1. / canair
  hcapcan  = cp * canair
  hcapcani = 1. / hcapcan

  if (iveg == 1) then  ! Vegetation case

     ! Vegetation is sufficiently abundant and not covered by snow.

     ! TAI and LAI reduced by ground snowcover (for low vegetation)

     stai = veg_tai * (1. - snowfac)
     slai = veg_lai * (1. - snowfac)

     ! Maximum allowed vegetation surface water

     vw_max = 0.2 * stai

     ! COMPUTE RDI WITH VEGETATION INFLUENCE

     ! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
     ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
     ! the values (zognd,0) when veg covered.

     zognd = soil_rough
     zoveg = veg_rough * (1. - snowfac) + zognd * snowfac
     zveg  = veg_height * (1. - snowfac)
     zdisp = zveg * .63
     !bob  rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd) &
     !bob         / (vonk * vonk * vels)

     ! Aerodynamic resistance (rd) between surface and canopy air are weighted
     ! between areas without and with vegetation.

     !bob   factv  = log(zts / zoveg) / (vonk * vonk * vels)
     factv  = 1. / (vonk * ustar)
     aux    = exp(exar * (1. - (zdisp + zoveg) / zveg))
     rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
     wtveg  = max(0.,min(1., 1.1 * veg_tai / covr))
     rdi    = ustar / (5. * (1. - wtveg) + ustar * rasveg * wtveg)

     ! Add any precipitation intercepted mass and energy to vegetation surface

     veg_water  = veg_water  + pcpg  * vf
     veg_energy = veg_energy + qpcpg * vf

     ! Check for presence of veg_water and set iwetveg flag

     iwetveg = (veg_water > wcap_vmin)

     ! If veg water is below a tiny threshold, remove from surface
     ! (conserving water mass and energy)

     if (.not. iwetveg) then
        canrrv     = canrrv + veg_water * canairi
        veg_energy = veg_energy - alvi * veg_water
        veg_water  = 0.
     endif

     ! Recompute vegetation temperature and liquid water fraction

     call qwtk(veg_energy, veg_water, hcapveg, veg_temp, fracliqv)

     tvegc = veg_temp - t00

     ! Remove any existing veg water over shedding threshold (which would have
     ! come from any dew or frost deposition on the previous timestep or precip
     ! this timestep). This will NOT change temperature or liquid fraction.

     if (veg_water > vw_max) then
        vegw_energy = veg_energy - hcapveg * tvegc

        wshed = veg_water - vw_max
        qshed = wshed * vegw_energy / veg_water

        veg_water  = vw_max
        veg_energy = veg_energy - qshed

        ! accumulated precip shed
        wshed2 = min(pcpg * vf, wshed)

        ! depth of accumulated precip shed
        dshed2 = wshed2 * max(1.e-3, dpcpg / max(pcpg,1.e-30))

        ! include depth of existing water shed (assume ice or liquid)
        dshed = dshed2 + (wshed - wshed2) * 1.e-3
     endif

     ! Vegetation saturation vapor density and derivative

     veg_rhovs  = rhovsil(tvegc)
     veg_rhovsp = rhovsil(tvegc+1.) - veg_rhovs

     ! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
     ! WITHOUT the LAI factor; the LAI factor is included later in computing
     ! xfers that involve rb and/or rc.

     rb = (1. + .5 * veg_tai) / (.01 * sqrt(ustar * c1))

     ! Set soil water potential to -600 m prior to search for wettest soil layer
     ! in root zone.  This value is low enough to effectively prevent transpiration.

     swp = -600.

     ! Initialize ktrans to zero prior to search for wettest soil layer in root zone

     ktrans = 0

     ! Loop over soil layers in the root zone

     do k = kroot(leaf_class), nzg

        ! Find layer in root zone with highest head (unless mostly ice)

        if (swp < head(k) .and. soil_fracliq(k) > .5) then
           swp = head(k)
           ktrans = k
        endif

     enddo

     ! Bypass stomatal resistance computation if root zone is too dry

     if (ktrans < 1) then

        rc = 1.e18

     else

        ! Compute saturation vapor pressure at veg temp

        esat_veg  = eslf(tvegc)

        ! Compute vapor pressure in canopy air from equation of state

        e_can = can_rhov * rvap * cantemp

        ! Compute vapor pressure at leaf surface using rc from previous timestep

        rc = stom_resist
        e_leaf = (rb * esat_veg + rc * e_can) / (rb + rc)
        vpd = max(0.,esat_veg - e_leaf)

        ! Evaluate 5 environmental factors and new rc
        ! (Swp multiplier converts from meters to Pascals (hydrostatic eqn for water)

        ftlo = 1. + exp(-stlo * (veg_temp - btlo))
        fthi = 1. + exp(-sthi * (veg_temp - bthi))
        frad = 1. + exp(-srad * (rshort   - brad))
        fswp = 1. + exp(-sswp * (swp*9810.- bswp))
        fvpd = 1. + exp(-svpd * (vpd      - bvpd))

        ! Compute asymptotoc value of stomatal resistance based on environmental factors

        rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(leaf_class)

        ! Update stomatal conductance assuming 15-minute response time

        tf = min(0.0011 * dt_leaf, 0.9)
        rc = rc / (1. + tf * (rc / rc_inf - 1.))

        ! Limit maximum transpiration to be <= 500 W/m^2 by increasing rc if necessary.

        transp_test = alvl * veg_lai * (veg_rhovs - can_rhov) / (rb + rc)
        if (transp_test > 500.) then
           rc = (rb + rc) * transp_test * .002 - rb
        endif

     endif

     stom_resist = rc

  else  ! No-vegetation case

     ! Aerodynamic conductance for bare soil or snow based on Garratt.

     rdi = .2 * ustar

     wshed  = 0.
     qshed  = 0.
     ktrans = 0

     stom_resist = 1.e18

  endif

  ! Apply contributions to sfcwater mass, depth, and energy_per_m2 from
  ! precipitation and shedding of moisture from vegetation.  Sfcwater
  ! quantities represent the NLSW1 layer in the sfcwater arrays in case
  ! nlev_sfcwater = 0.

  wadd = pcpg * (1. - vf) + wshed

  if (wadd > 1.e-30) then
     sfcwater_mass   = sfcwater_mass   + wadd
     sfcwater_depth  = sfcwater_depth  + dpcpg * (1. - vf) + dshed
     energy_per_m2   = energy_per_m2   + qpcpg * (1. - vf) + qshed
  endif

  ! If nlsw1 = 1 (lowest sfcwater layer), and sfcwater mass is present but
  ! below limiting stability value for explicit heat transfer or if sfcwater
  ! is mostly liquid, perform implicit thermal balance between sfcwater and
  ! the top soil layer and set icomb flag to 1.

  if (nlsw1 == 1 .and. &
     (sfcwater_mass < snowmin_expl .or. sfcwater_fracliq > .5)) then

     call sfcwater_soil_comb(iland, iwsfc,                  &
                             soil_water(nzg),soil_energy(nzg), &
                             specifheat_drysoil(nzg), sfcwater_mass,energy_per_m2, &
                             sfcwater_tempk,sfcwater_fracliq)

     icomb = 1

     ! Diagnose combined sfc heat capacity (hcapsfc) for sfcwater and top soil layer.
     ! See description of hcapsfc for sfcwater alone in the ELSE part of the IF block.

     sfc_energy = energy_per_m2 + soil_energy(nzg) * dslz(nzg)
     sfc_wmass  = sfcwater_mass + soil_water(nzg)  * dslz(nzg) * 1.e3
     sfc_soilhc = specifheat_drysoil(nzg) * dslz(nzg)

     if     (sfc_energy < 0.) then
        hcapsfc = sfc_wmass * cice + sfc_soilhc
     elseif (sfc_energy > sfc_wmass * alli) then
        hcapsfc = sfc_wmass * cliq + sfc_soilhc
     else
        hcapsfc = sfc_wmass * (alli / 1.) + sfc_soilhc ! Assume small positive dT/dE
     endif

  else

     ! If sfcwater and soil(nzg) temperatures are NOT implicitly coupled, set
     ! icomb flag to 0 and diagnose sfc heat capacity (hcapsfc) for sfcwater alone.

     icomb = 0

     ! First, diagnose sfcwater temperature and liquid fraction.  Use qwtk
     ! instead of qtk because energy_per_m2 is used instead of sfcwater_energy.
     ! ("dryhcap" = 100 is very small value)

     call qwtk(energy_per_m2,sfcwater_mass,100.,sfcwater_tempk,sfcwater_fracliq)

     if     (energy_per_m2 < 0.) then
        hcapsfc = sfcwater_mass * cice
     elseif (energy_per_m2 > sfcwater_mass * alli) then
        hcapsfc = sfcwater_mass * cliq
     else
        hcapsfc = sfcwater_mass * (alli / 1.)  ! Assume small but nonzero dT/dE
     endif

  endif

  hcapsfci = 1. / hcapsfc

  ! Use the following latent heat for vapor flux at the surface when multiplying
  ! inverse surface heat capacity to get implicit temperature change of the
  ! surface.  (When surface water is a mixture of liquid and ice, inverse surface
  ! heat capacity is zero, so alvsfc has no effect.)

  alvsfc = alvi - sfcwater_fracliq * alli

  ! Sfcwater saturation vapor density and derivative

  sfc_tempk = sfcwater_tempk

  sfc_tempk1 = sfc_tempk + 1.

  sfc_rhovs  = rhovsil(sfc_tempk -273.15)
  sfc_rhovs1 = rhovsil(sfc_tempk1-273.15)
  sfc_rhovsp = sfc_rhovs1 - sfc_rhovs

  ! Reduced saturation vapor density and gradient of soil during evaporation

  call grndvap_ab(iland, sfc_tempk, soil_water(nzg), wresid_vg(nzg), soilfldcap, &
                  head(nzg), alpha, beta)

  gnd_rhov  = alpha * sfc_rhovs
  gnd_rhovp = alpha * sfc_rhovsp

  ! The remainder of subroutine canopy sets up and solves a linear system of
  ! equations using the trapezoidal-implicit method.  The solution of the system
  ! consists of turbulent heat and water vapor fluxes between canopy air, the
  ! surface (sfcwater or soil), vegetation (if present), and the free atmosphere,
  ! plus the consequent changes to water and energy content of canopy air, the surface,
  ! and vegetation.  The solution obtains values for the following quantities:

  !  1. Surface to canopy air heat transfer
  !  2. Surface to canopy air water vapor transfer
  !  3. Vegetation to canopy air heat transfer
  !  4. Vegetation to canopy air water vapor transfer from transpiration
  !  5. Vegetation to canopy air water vapor transfer from veg surface water
  !  6. Canopy air to free atmosphere heat transfer
  !  7. Canopy air to free atmosphere vapor transfer
  !  8. Canopy air temperature update
  !  9. Canopy specific humidity update
  ! 10. Vegetation temperature update
  ! 11. Vegetation surface water update

  ! External forcing terms for this system are:

  ! 1. Radiative fluxes to vegetation and to the surface

  ! (Precipitation/shedding fluxes to vegetation and to the surface were already
  ! applied above.  Updates to surface water and energy are applied later.)

  ! Quantities that apply both WITH and WITHOUT VEGETATION
  ! [radsfc is net shortwave + longwave input to surface].  In LEAF4, all shortwave
  ! radiation absorbed by the surface is put into rshort_g (in subroutine
  ! sfcrad_land), whether sfcwater is present or not.  Here in subroutine
  ! leaf4_canopy, the decision is made, based on abundance of sfcwater, whether
  ! rshort_g goes to the top sfcwater layer or to the combined sfcwater-soil surface.

  ! Define sfcwater_fracarea to be 0.0 with zero sfcwater_mass, to be 1.0 when
  ! sfcwater_mass is at or above a threshold value (currently set to 2 kg/m^2),
  ! and to increase as the sqrt of sfcwater_mass below the threshold.
  ! Sfcwater_fracarea partitions surface evaporation into sfcwater and bare
  ! ground portions.  The sfcwater_mass threshold should be set large enough
  ! that sfcwater_mass cannot all evaporate in a single model time step.  The
  ! sqrt functional form (in which sfcwater_fracarea increases very rapidly
  ! at low sfcwater_mass) ensures that small amounts of sfcwater_mass will be
  ! able to evaporate in a timely manner.

  if (sfcwater_mass > 0.) then
     sfcwater_fracarea = sqrt(min(1.,0.5 * sfcwater_mass))
  else
     sfcwater_fracarea = 0.
  endif

  radsfc = dt_leaf * (rlong_s + rshort_g + rlong_g)

  a5  = dt_leaf * rdi                          ! sfc vap xfer coef  (all area)
  a6  = cp * rhos * a5                         ! sfc heat xfer coef (all area)
  a7  = a5 * sfcwater_fracarea                 ! sfc vap xfer coef  (wet areas only)
  a8  = a5 * beta * (1. - sfcwater_fracarea)   ! sfc vap xfer coef  (soil areas only)
  a9  = dt_leaf * vkhsfc / sfcg%dzt_bot(iwsfc) ! can-atm vap xfer coef
  a10 = cp * a9                                ! can-atm heat xfer coef 

  h2 = fcn * sfc_rhovsp * hcapsfci
  h3 = fcn * gnd_rhovp  * hcapsfci
  h4 = fcn * rhos * canairi
  h6 = fcn * hcapsfci
  h7 = fcn * hcapcani
  h8 = fcn * canairi

  y2  = sfc_rhovs - can_rhov + h2 * radsfc
  y5  = sfc_tempk - cantemp  + h6 * radsfc
  y3  = gnd_rhov  - can_rhov + h3 * radsfc
  y9  = canrrv    - airrrv
  y10 = cantemp   - canexner * airtheta

  if (iveg == 0) then

     ! Case with NO VEGETATION

     ! Set all skiptest values to .false.

     skiptest(:) = .false.

     ! Fill arrays for matrix solution.

! Test 101:  Solve 4x4 system to test for either condensation onto surface or
!            partial evaporation from sfcwater that completely covers the surface.
!            (When sfcwater completely covers the surface, it is too abundant to
!            completely evaporate in a single time step.)

     aa4(1,1) = 1._r8 + a5 * (h2 * alvsfc + h4)
     aa4(1,2) =         a5 * h2
     aa4(1,3) =       - a5 * h4
     aa4(1,4) = 0._r8
     yy4(1)   =         a5 * y2         ! WSC row

     aa4(2,1) =         a6 * h6 * alvsfc
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) = 0._r8
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * y5         ! HSC row

     aa4(3,1) =       - a9 * h8
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a9 * h8
     aa4(3,4) = 0._r8
     yy4(3)   =         a9 * y9         ! WCA row

     aa4(4,1) = 0._r8
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * h7
     yy4(4)   =         a10 * y10       ! HCA row

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land101',4,aa4,yy4,glatw,glonw)

     if (xx4(1) < eps_wxfer .or.  &
        (xx4(1) <= sfcwater_mass .and. sfcwater_fracarea > 0.999)) then ! test was successful

        wxfersc = xx4(1)
        hxfersc = xx4(2)
        wxferca = xx4(3)
        hxferca = xx4(4)
        wxfergc = 0.
        itest = 101
        go to 120
     endif

! All remaining situations do NOT have complete coverage of surface by sfcwater
! and do NOT have condensation onto surface.

     if (sfcwater_fracarea > 0.999) then
        write(*,'(a,i10,2f9.2)') 'NOVEG CASE sfcwater_fracarea should be less than 1', &
            iwsfc,glatw,glonw
        print*, 'sf1 ',sfcwater_fracarea,eps_wxfer,sfcwater_mass
        print*, 'sf2 ',xx4(1),xx4(2),xx4(3),xx4(4)
        print*, 'sf3 ',a5,a6,a9,a10
        print*, 'sf4 ',h2,h4,h6,h7,h8
        print*, 'sf5 ',y2,y5,y9,y10
        stop 'sfcwater_fracarea = 1, NOVEG CASE'
     endif

! Test 102:  If some sfcwater is present, solve 5x5 system to test for
!            partial evaporation from sfcwater and evaporation from soil

     if (sfcwater_mass > 0.) then

        aa5(1,1) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa5(1,2) =         a7 * h2
        aa5(1,3) =         a7 * h4
        aa5(1,4) =       - a7 * h4
        aa5(1,5) = 0._r8
        yy5(1)   =         a7 * y2    ! WSC row

        aa5(2,1) =         a6 * h6 * alvsfc
        aa5(2,2) = 1._r8 + a6 * (h6 + h7)
        aa5(2,3) =         a6 * h6 * alvsfc
        aa5(2,4) = 0._r8
        aa5(2,5) =       - a6 * h7
        yy5(2)   =         a6 * y5    ! HSC row

        aa5(3,1) =         a8 * h4
        aa5(3,2) =         a8 * h3
        aa5(3,3) = 1._r8 + a8 * (h3 * alvsfc + h4)
        aa5(3,4) =       - a8 * h4
        aa5(3,5) = 0._r8
        yy5(3)   =         a8 * y3    ! WGC row

        aa5(4,1) =       - a9 * h8
        aa5(4,2) = 0._r8
        aa5(4,3) =       - a9 * h8
        aa5(4,4) = 1._r8 + a9 * h8
        aa5(4,5) = 0._r8
        yy5(4)   =         a9 * y9         ! WCA row

        aa5(5,1) = 0._r8
        aa5(5,2) =       - a10 * h7
        aa5(5,3) = 0._r8
        aa5(5,4) = 0._r8
        aa5(5,5) = 1._r8 + a10 * h7
        yy5(5)   =         a10 * y10       ! HCA row

        call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land102',5,aa5,yy5,glatw,glonw)

        if (xx5(1) > -eps_wxfer .and. xx5(1) <= sfcwater_mass .and. &
            xx5(3) > -eps_wxfer) then ! test was successful

           wxfersc = xx5(1)
           hxfersc = xx5(2)
           wxfergc = xx5(3)
           wxferca = xx5(4)
           hxferca = xx5(5)
           itest = 102
           go to 120
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 103, 104, and 105

        if (xx5(1) < -eps_wxfer) then
           skiptest(3) = .true.
           skiptest(4) = .true.
           skiptest(5) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 103

        if (xx5(1) <= sfcwater_mass) then
           skiptest(3) = .true.
        endif

     endif

     103 continue
     if (skiptest(3)) go to 104

! Test 103: Given complete evaporation of any sfcwater, solve 4x4 system to
!           test for positive evaporation from soil

     aa4(1,1) = 1._r8 + a8 * (h3 * alvsfc + h4)
     aa4(1,2) =         a8 * h3
     aa4(1,3) =       - a8 * h4
     aa4(1,4) = 0._r8
     yy4(1)   =         a8 * (y3 - (h3 * alvsfc + h4) * sfcwater_mass)  ! WGC row

     aa4(2,1) =         a6 * h6 * alvsfc
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) = 0._r8
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)         ! HSC row

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

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land103',4,aa4,yy4,glatw,glonw)

     if (xx4(1) > -eps_wxfer) then ! test was successful

        wxfersc = sfcwater_mass
        wxfergc = xx4(1)
        hxfersc = xx4(2)
        wxferca = xx4(3)
        hxferca = xx4(4)
        itest = 103
        go to 120
     endif

     104 continue
     if (skiptest(4)) go to 105

! Test 104: If some sfcwater is present, and given zero vapor flux with soil,
!           solve 4x4 system to test for partial evaporation of sfcwater

     if (sfcwater_mass > 0.) then

        aa4(1,1) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa4(1,2) =         a7 * h2
        aa4(1,3) =       - a7 * h4
        aa4(1,4) = 0._r8
        yy4(1)   =         a7 * y2        ! WSC row

        aa4(2,1) =         a6 * h6 * alvsfc
        aa4(2,2) = 1._r8 + a6 * (h6 + h7)
        aa4(2,3) = 0._r8
        aa4(2,4) =       - a6 * h7
        yy4(2)   =         a6 * y5        ! HSC row

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

        call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land104',4,aa4,yy4,glatw,glonw)

        if (xx4(1) > -eps_wxfer .and. xx4(1) < sfcwater_mass) then ! test was successful

           wxfersc = xx4(1)
           hxfersc = xx4(2)
           wxferca = xx4(3)
           hxferca = xx4(4)
           wxfergc = 0.
           itest = 104
           go to 120
        endif

     endif

     105 continue
     if (skiptest(5)) go to 106

! Test 105: Given zero vapor flux with soil and given complete evaporation of
!           any sfcwater, solve 3x3 system (a reduction of the 4x4 system in
!           Test 104) to get surface heat flux

     aa3(1,1) = 1._r8 + a6 * (h6 + h7)
     aa3(1,2) = 0._r8
     aa3(1,3) =       - a6 * h7
     yy3(1)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)  ! HSC row

     aa3(2,1) = 0._r8
     aa3(2,2) = 1._r8 + a9 * h8
     aa3(2,3) = 0._r8
     yy3(2)   =         a9 * y9   ! WCA row

     aa3(3,1) =       - a10 * h7
     aa3(3,2) = 0._r8
     aa3(3,3) = 1._r8 + a10 * h7
     yy3(3)   =         a10 * y10 ! HCA row

     call matrix8_3x3(aa3,yy3,xx3,sing); if (sing) call sing_print(iwsfc,'land105',3,aa3,yy3,glatw,glonw)

     wxfersc = sfcwater_mass
     wxfergc = 0.
     hxfersc = xx3(1)
     wxferca = xx3(2)
     hxferca = xx3(3)

     itest = 105
     go to 120

     106 continue

     ! If we got here, none of the above tests was satisfied, an outcome that should not happen

     itest = 92

     print*, 'itest = 92 result obtained in leaf4_canopy ',iwsfc,glatw,glonw
     stop 'stop test 92 '

     120 continue

     ! Update components

     cantemp = cantemp + (hxfersc - hxferca) * hcapcani

     ! The latent heat of sublimation, alvi, must be used in computing the following gain/loss of
     ! energy_per_m2 due to vapor flux.  This is because extensive energy and mass carried by the
     ! vapor are both explicitly added to or subtracted from surface water, and the energy of vapor,
     ! defined relative to a zero-energy reference of ice at 0 deg C, is its mass times alvi
     ! (neglecting a small contribution from vapor temperature differing from 0 deg C).  Following
     ! the gain or loss of mass and energy via vapor flux, the resultant intensive energy of surface
     ! water is diagnosed as its extensive energy divided by its extensive mass.  If surface water
     ! is in the liquid phase, its intensive energy includes the latent heat of fusion.  Vapor
     ! carries this energy by virtue of its mass flux alone, so the total energy alvi carried by
     ! vapor is only alvl above this value, which correctly accounts for evaporation/condensation
     ! occurring in this case rather than sublimation/deposition.

     energy_per_m2 = energy_per_m2 + radsfc - hxfersc - (wxfersc + wxfergc) * alvi

     canrrv = canrrv + (wxfersc + wxfergc - wxferca) * canairi

     sfluxt = hxferca / dt_leaf
     sfluxr = wxferca / dt_leaf

     delw            = dslzi(nzg) * wxfergc * .001
     fact            = wsat_vg(nzg) - wresid_vg(nzg)

     soil_water(nzg) = soil_water(nzg) - delw
     head      (nzg) = head      (nzg) - delw * head_slope(nzg)
     soil_wfrac(nzg) = soil_wfrac(nzg) - delw / fact

     soil_wfrac(nzg) = max(0., min(1., soil_wfrac(nzg)))

     sfcwater_mass = sfcwater_mass - wxfersc

     specvol = .001
     if (wxfersc > 0. .and. sfcwater_mass > 1.e-3) specvol = sfcwater_depth / sfcwater_mass

     sfcwater_depth = sfcwater_depth - wxfersc * specvol

     transp = 0.

  else

     ! Case WITH VEGETATION

     if     (veg_energy < 0.) then
        hcapvegw = veg_water * cice + hcapveg
     elseif (veg_energy > veg_water * alli) then
        hcapvegw = veg_water * cliq + hcapveg
     else
        hcapvegw = veg_water * (alli / 1.) + hcapveg ! Assume small positive dT/dE
     endif

     hcapvegwi = 1. / hcapvegw

     ! Use the following latent heat for vapor flux at the vegetation surface
     ! when multiplying inverse vegw heat capacity to get implicit temperature
     ! change of the vegetation.  (When veg_water is a mixture of liquid and ice,
     ! inverse veg heat capacity is zero, so alvveg has no effect.)

     alvveg = alvi - fracliqv * alli

     radveg = dt_leaf * (rshort_v + rlong_v)

     ! veg_water fractional coverage

     sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66666666)

     a1 = dt_leaf * 2. * stai / rb                  ! veg vap xfer coef  (all area)
     a2 = cp * rhos * a1                            ! veg heat xfer coef (all area)
     a3 = a1 * sigmaw                               ! veg vap xfer coef  (veg_water evap)
     a4 = dt_leaf * slai * (1.-sigmaw) / (rb + rc)  ! veg vap xfer coef  (transp)

     ! Auxiliary quantities

     h1 = fcn * veg_rhovsp * hcapvegwi
     h5 = fcn * hcapvegwi

     y1 = veg_rhovs - can_rhov + h1 * radveg
     y4 = veg_temp  - cantemp  + h5 * radveg

     ! Fill arrays for matrix solution (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, vegetation, and the surface.
     ! For vegetation, there are 3 distinct situations to consider:

     ! 1. Condensation onto vegetation (Vc)
     ! 2. Evaporation from vegetation with veg water not depleted (Vew)
     ! 3. Evaporation from vegetation with veg water depleted or zero (Ve)

     ! For the surface, there are 5 distinct situations to consider:

     ! 1. Condensation onto surface (Sc) or full sfcwater coverage with partial evaporation (Sew)
     ! 2. Evaporation from both sfcwater and soil with sfcwater not depleting (GeSew)
     ! 3. Evaporation from both sfcwater and soil with sfcwater depleting or zero (GeSe)
     ! 4. Evaporation from sfcwater only with sfcwater not depleting (Sew)
     ! 5. Evaporation from sfcwater only with sfcwater depleting or zero (Se)

     ! (Zero vapor flux from exposed soil occurs when soil is too warm for condensation
     ! and too dry for evaporation.)

     ! There are 15 possible pairs containing one item from each list,
     ! only one of which is the correct description for a particular land
     ! grid cell and timestep.

     ! 6x6  Test  1  Vc  (Sc or Sew)               (or sfcwater_fracarea = 1)
     ! 6x6  Test  2  Vew (Sc or Sew)  (iwetveg=t)  (or sfcwater_fracarea = 1)
     ! 6x6  Test  3  Ve  (Sc or Sew)               (or sfcwater_fracarea = 1)

     ! 7x7  Test  4  Vc  GeSew                     (sfcwater_mass > 0.)
     ! 7x7  Test  5  Vew GeSew        (iwetveg=t)  (sfcwater_mass > 0.)
     ! 7x7  Test  6  Ve  GeSew                     (sfcwater_mass > 0.)

     ! 6x6  Test  7  Vc  GeSe
     ! 6x6  Test  8  Vew GeSe         (iwetveg=t)
     ! 6x6  Test  9  Ve  GeSe

     ! 6x6  Test 10  Vc  Sew                       (sfcwater_mass > 0.)
     ! 6x6  Test 11  Vew Sew          (iwetveg=t)  (sfcwater_mass > 0.)
     ! 6x6  Test 12  Ve  Sew                       (sfcwater_mass > 0.)

     ! 5x5  Test 13  Vc  Se
     ! 5x5  Test 14  Vew Se           (iwetveg=t)
     ! 5x5  Test 15  Ve  Se

     ! We proceed to solve each equation in the above order until finding one
     ! that produces fluxes of the correct sign for the particular equation.

     ! Set all skiptest values to .false.

     skiptest(:) = .false.

! Test 1: Solve 6x6 system to test for condensation onto vegetation and either
!         condensation onto surface or partial evaporation from sfcwater that
!         completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     aa6(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa6(1,2) =         a1 * h4
     aa6(1,3) =         a1 * h1
     aa6(1,4) = 0._r8
     aa6(1,5) =       - a1 * h4
     aa6(1,6) = 0._r8
     yy6(1)   =         a1 * y1   ! WVC row

     aa6(2,1) =         a5 * h4
     aa6(2,2) = 1._r8 + a5 * (h2 * alvsfc + h4)
     aa6(2,3) = 0._r8
     aa6(2,4) =         a5 * h2
     aa6(2,5) =       - a5 * h4
     aa6(2,6) = 0._r8
     yy6(2)   =         a5 * y2   ! WSC row

     aa6(3,1) =         a2 * h5 * alvveg
     aa6(3,2) = 0._r8
     aa6(3,3) = 1._r8 + a2 * (h5 + h7)
     aa6(3,4) =         a2 * h7
     aa6(3,5) = 0._r8
     aa6(3,6) =       - a2 * h7
     yy6(3)   =         a2 * y4   ! HVC row

     aa6(4,1) = 0._r8
     aa6(4,2) =         a6 * h6 * alvsfc
     aa6(4,3) =         a6 * h7
     aa6(4,4) = 1._r8 + a6 * (h6 + h7)
     aa6(4,5) = 0._r8
     aa6(4,6) =       - a6 * h7
     yy6(4)   =         a6 * y5   ! HSC row

     aa6(5,1) =       - a9 * h8
     aa6(5,2) =       - a9 * h8
     aa6(5,3) = 0._r8
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a9 * h8
     aa6(5,6) = 0._r8
     yy6(5)   =         a9 * y9   ! WCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a10 * h7
     aa6(6,4) =       - a10 * h7
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a10 * h7
     yy6(6)   =         a10 * y10 ! HCA row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land1',6,aa6,yy6,glatw,glonw)

     if (xx6(1) < eps_wxfer .and. &
        (xx6(2) < eps_wxfer .or. (xx6(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

        wxfervc = xx6(1)
        wxfersc = xx6(2)
        hxfervc = xx6(3)
        hxfersc = xx6(4)
        wxferca = xx6(5)
        hxferca = xx6(6)
        wxfergc = 0.
        transp  = 0.
        itest = 1
        go to 20
     endif

     ! If this test resulted in upward flux from the surface with
     ! sfcwater_fracarea less than 1, skip past tests 2 and 3

     if (xx6(2) > eps_wxfer .and. sfcwater_fracarea <= 0.999) then
        skiptest(2) = .true.
        skiptest(3) = .true.
     endif

     2 continue
     if (skiptest(2)) go to 3

! Test 2: If vegetation is wet, solve 6x6 system to test for partial evaporation from
!         vegetation and either condensation onto surface or partial evaporation from
!         sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     if (iwetveg) then

        aa6(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
        aa6(1,2) =         (a3 + a4) * h4
        aa6(1,3) =         (a3 + a4) * h1
        aa6(1,4) = 0._r8
        aa6(1,5) =       - (a3 + a4) * h4
        aa6(1,6) = 0._r8
        yy6(1)   =         (a3 + a4) * y1  ! WVC row

        aa6(2,1) =         a5 * h4
        aa6(2,2) = 1._r8 + a5 * (h2 * alvsfc + h4)
        aa6(2,3) = 0._r8
        aa6(2,4) =         a5 * h2
        aa6(2,5) =       - a5 * h4
        aa6(2,6) = 0._r8
        yy6(2)   =         a5 * y2         ! WSC row

        aa6(3,1) =         a2 * h5 * alvveg
        aa6(3,2) = 0._r8
        aa6(3,3) = 1._r8 + a2 * (h5 + h7)
        aa6(3,4) =         a2 * h7
        aa6(3,5) = 0._r8
        aa6(3,6) =       - a2 * h7
        yy6(3)   =         a2 * y4         ! HVC row

        aa6(4,1) = 0._r8
        aa6(4,2) =         a6 * h6 * alvsfc
        aa6(4,3) =         a6 * h7
        aa6(4,4) = 1._r8 + a6 * (h6 + h7)
        aa6(4,5) = 0._r8
        aa6(4,6) =       - a6 * h7
        yy6(4)   =         a6 * y5         ! HSC row

        aa6(5,1) =       - a9 * h8
        aa6(5,2) =       - a9 * h8
        aa6(5,3) = 0._r8
        aa6(5,4) = 0._r8
        aa6(5,5) = 1._r8 + a9 * h8
        aa6(5,6) = 0._r8
        yy6(5)   =         a9 * y9   ! WCA row

        aa6(6,1) = 0._r8
        aa6(6,2) = 0._r8
        aa6(6,3) =       - a10 * h7
        aa6(6,4) =       - a10 * h7
        aa6(6,5) = 0._r8
        aa6(6,6) = 1._r8 + a10 * h7
        yy6(6)   =         a10 * y10 ! HCA row

        call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land2',6,aa6,yy6,glatw,glonw)

        evap = xx6(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx6(1) > -eps_wxfer .and. evap <= veg_water .and. &
           (xx6(2) <  eps_wxfer .or. (xx6(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

           wxfervc = evap
           wxfersc = xx6(2)
           hxfervc = xx6(3)
           hxfersc = xx6(4)
           wxferca = xx6(5)
           hxferca = xx6(6)
           wxfergc = 0.
           transp  = xx6(1) - evap
           itest = 2
           go to 20
        endif

        ! If this test resulted in upward flux from the surface with
        ! sfcwater_fracarea less than 1, skip past test 3

        ! If this test resulted in less than complete evaporation of veg_water
        ! skip past test 3

        if (xx6(2) > eps_wxfer .and. sfcwater_fracarea <= 0.999) then
           skiptest(3) = .true.
        elseif (evap <= veg_water) then
           skiptest(3) = .true.
        endif

     endif

     3 continue
     if (skiptest(3)) go to 4

! Test 3: Given complete evaporation of any veg_water, solve 6x6 system to test
!         for transpiration and either condensation onto surface or partial
!         evaporation from sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     aa6(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa6(1,2) =         a4 * h4
     aa6(1,3) =         a4 * h1
     aa6(1,4) = 0._r8
     aa6(1,5) =       - a4 * h4
     aa6(1,6) = 0._r8
     yy6(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water)  ! WVC row

     aa6(2,1) =         a5 * h4
     aa6(2,2) = 1._r8 + a5 * (h2 * alvsfc + h4)
     aa6(2,3) = 0._r8
     aa6(2,4) =         a5 * h2
     aa6(2,5) =       - a5 * h4
     aa6(2,6) = 0._r8
     yy6(2)   =         a5 * (y2 - h4 * veg_water)                ! WSC row

     aa6(3,1) =         a2 * h5 * alvveg
     aa6(3,2) = 0._r8
     aa6(3,3) = 1._r8 + a2 * (h5 + h7)
     aa6(3,4) =         a2 * h7
     aa6(3,5) = 0._r8
     aa6(3,6) =       - a2 * h7
     yy6(3)   =         a2 * (y4 - h5 * alvveg * veg_water)         ! HVC row

     aa6(4,1) = 0._r8
     aa6(4,2) =         a6 * h6 * alvsfc
     aa6(4,3) =         a6 * h7
     aa6(4,4) = 1._r8 + a6 * (h6 + h7)
     aa6(4,5) = 0._r8
     aa6(4,6) =       - a6 * h7
     yy6(4)   =         a6 * y5                                   ! HSC row

     aa6(5,1) =       - a9 * h8
     aa6(5,2) =       - a9 * h8
     aa6(5,3) = 0._r8
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a9 * h8
     aa6(5,6) = 0._r8
     yy6(5)   =         a9 * y9   ! WCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a10 * h7
     aa6(6,4) =       - a10 * h7
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a10 * h7
     yy6(6)   =         a10 * y10 ! HCA row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land3',6,aa6,yy6,glatw,glonw)

     if (xx6(1) > -eps_wxfer .and. &
        (xx6(2) <  eps_wxfer .or. (xx6(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

        wxfervc = veg_water
        wxfersc = xx6(2)
        hxfervc = xx6(3)
        hxfersc = xx6(4)
        wxferca = xx6(5)
        hxferca = xx6(6)
        wxfergc = 0.
        transp  = xx6(1)
        itest = 3
        go to 20
     endif

     4 continue

! All remaining situations do NOT have complete coverage of surface by sfcwater
! and do NOT have condensation onto surface.

     if (sfcwater_fracarea > 0.999) then
        write(*,'(a,i10,2f9.2)') 'VEG CASE sfcwater_fracarea should be less than 1', &
            iwsfc,glatw,glonw
        print*, 'sf11 ',sfcwater_fracarea,eps_wxfer,sfcwater_mass
        print*, 'sf12 ',xx4(1),xx4(2),xx4(3),xx4(4)
        print*, 'sf13 ',veg_water
        stop 'sfcwater_fracarea = 1, VEG CASE'
     endif

     if (skiptest(4)) go to 5

! Test 4: If some sfcwater is present, solve 7x7 system to test for
!         condensation onto vegetation, partial evaporation of sfcwater,
!         and evaporation from soil

     if (sfcwater_mass > 0.) then

        aa7(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
        aa7(1,2) =         a1 * h4
        aa7(1,3) =         a1 * h1
        aa7(1,4) = 0._r8
        aa7(1,5) =         a1 * h4
        aa7(1,6) =       - a1 * h4
        aa7(1,7) = 0._r8
        yy7(1)   =         a1 * y1  ! WVC row

        aa7(2,1) =         a7 * h4
        aa7(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa7(2,3) = 0._r8
        aa7(2,4) =         a7 * h2
        aa7(2,5) =         a7 * h4
        aa7(2,6) =       - a7 * h4
        aa7(2,7) = 0._r8
        yy7(2)   =         a7 * y2  ! WSC row

        aa7(3,1) =         a2 * h5 * alvveg
        aa7(3,2) = 0._r8
        aa7(3,3) = 1._r8 + a2 * (h5 + h7)
        aa7(3,4) =         a2 * h7
        aa7(3,5) =         a2 * h5 * alvveg
        aa7(3,6) = 0._r8
        aa7(3,7) =       - a2 * h7
        yy7(3)   =         a2 * y4  ! HVC row

        aa7(4,1) = 0._r8
        aa7(4,2) =         a6 * h6 * alvsfc
        aa7(4,3) =         a6 * h7
        aa7(4,4) = 1._r8 + a6 * (h6 + h7)
        aa7(4,5) =         a6 * h6 * alvsfc
        aa7(4,6) = 0._r8
        aa7(4,7) =       - a6 * h7
        yy7(4)   =         a6 * y5  ! HSC row

        aa7(5,1) =         a8 * h4
        aa7(5,2) =         a8 * h4
        aa7(5,3) = 0._r8
        aa7(5,4) =         a8 * h3
        aa7(5,5) = 1._r8 + a8 * (h3 * alvsfc + h4)
        aa7(5,6) =       - a8 * h4
        aa7(5,7) = 0._r8
        yy7(5)   =         a8 * y3  ! WGC row

        aa7(6,1) =       - a9 * h8
        aa7(6,2) =       - a9 * h8
        aa7(6,3) = 0._r8
        aa7(6,4) = 0._r8
        aa7(6,5) =       - a9 * h8
        aa7(6,6) = 1._r8 + a9 * h8
        aa7(6,7) = 0._r8
        yy7(6)   =         a9 * y9   ! WCA row

        aa7(7,1) = 0._r8
        aa7(7,2) = 0._r8
        aa7(7,3) =       - a10 * h7
        aa7(7,4) =       - a10 * h7
        aa7(7,5) = 0._r8
        aa7(7,6) = 0._r8
        aa7(7,7) = 1._r8 + a10 * h7
        yy7(7)   =         a10 * y10 ! HCA row

        call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land4',7,aa7,yy7,glatw,glonw)

        if (xx7(1) <  eps_wxfer .and. &
            xx7(2) > -eps_wxfer .and. xx7(2) <= sfcwater_mass .and. &
            xx7(5) > -eps_wxfer) then ! test was successful

           wxfervc = xx7(1)
           wxfersc = xx7(2)
           hxfervc = xx7(3)
           hxfersc = xx7(4)
           wxfergc = xx7(5)
           wxferca = xx7(6)
           hxferca = xx7(7)
           transp  = 0.
           itest = 4
           go to 20
        endif

        ! If this test resulted in condensation onto vegetation,
        ! skip past tests 5 and 6

        if (xx7(1) < -eps_wxfer) then
           skiptest(5) = .true.
           skiptest(6) = .true.
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 7, 10, and 13

        if (xx7(2) < -eps_wxfer) then
           skiptest( 7) = .true.
           skiptest(10) = .true.
           skiptest(13) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 7

        if (xx7(2) <= sfcwater_mass) then
           skiptest(7) = .true.
        endif

     endif

     5 continue
     if (skiptest(5)) go to 6

! Test 5: If some sfcwater is present, and also if vegetation is wet, solve 7x7
!         system to test for partial evaporation from vegetation, partial
!         evaporation of sfcwater, and zero or positive evaporation from soil

     if (sfcwater_mass > 0. .and. iwetveg) then

        aa7(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
        aa7(1,2) =         (a3 + a4) * h4
        aa7(1,3) =         (a3 + a4) * h1
        aa7(1,4) = 0._r8
        aa7(1,5) =         (a3 + a4) * h4
        aa7(1,6) =       - (a3 + a4) * h4
        aa7(1,7) = 0._r8
        yy7(1)   =         (a3 + a4) * y1  ! WVC row

        aa7(2,1) =         a7 * h4
        aa7(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa7(2,3) = 0._r8
        aa7(2,4) =         a7 * h2
        aa7(2,5) =         a7 * h4
        aa7(2,6) =       - a7 * h4
        aa7(2,7) = 0._r8
        yy7(2)   =         a7 * y2         ! WSC row

        aa7(3,1) =         a2 * h5 * alvveg
        aa7(3,2) = 0._r8
        aa7(3,3) = 1._r8 + a2 * (h5 + h7)
        aa7(3,4) =         a2 * h7
        aa7(3,5) =         a2 * h5 * alvveg
        aa7(3,6) = 0._r8
        aa7(3,7) =       - a2 * h7
        yy7(3)   =         a2 * y4         ! HVC row

        aa7(4,1) = 0._r8
        aa7(4,2) =         a6 * h6 * alvsfc
        aa7(4,3) =         a6 * h7
        aa7(4,4) = 1._r8 + a6 * (h6 + h7)
        aa7(4,5) =         a6 * h6 * alvsfc
        aa7(4,6) = 0._r8
        aa7(4,7) =       - a6 * h7
        yy7(4)   =         a6 * y5         ! HSC row

        aa7(5,1) =         a8 * h4
        aa7(5,2) =         a8 * h4
        aa7(5,3) = 0._r8
        aa7(5,4) =         a8 * h3
        aa7(5,5) = 1._r8 + a8 * (h3 * alvsfc + h4)
        aa7(5,6) =       - a8 * h4
        aa7(5,7) = 0._r8
        yy7(5)   =         a8 * y3         ! WGC row

        aa7(6,1) =       - a9 * h8
        aa7(6,2) =       - a9 * h8
        aa7(6,3) = 0._r8
        aa7(6,4) = 0._r8
        aa7(6,5) =       - a9 * h8
        aa7(6,6) = 1._r8 + a9 * h8
        aa7(6,7) = 0._r8
        yy7(6)   =         a9 * y9   ! WCA row

        aa7(7,1) = 0._r8
        aa7(7,2) = 0._r8
        aa7(7,3) =       - a10 * h7
        aa7(7,4) =       - a10 * h7
        aa7(7,5) = 0._r8
        aa7(7,6) = 0._r8
        aa7(7,7) = 1._r8 + a10 * h7
        yy7(7)   =         a10 * y10 ! HCA row

        call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land5',7,aa7,yy7,glatw,glonw)

        evap = xx7(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx7(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx7(2) > -eps_wxfer .and. xx7(2) <= sfcwater_mass .and. &
            xx7(5) > -eps_wxfer) then ! test was successful

           wxfervc = evap
           wxfersc = xx7(2)
           hxfervc = xx7(3)
           hxfersc = xx7(4)
           wxfergc = xx7(5)
           wxferca = xx7(6)
           hxferca = xx7(7)
           transp  = xx7(1) - evap
           itest = 5
           go to 20
        endif

        ! If this test resulted in less than complete evaporation of veg_water,
        ! skip past test 6

        if (evap <= veg_water) then
           skiptest(6) = .true.
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 8, 11, and 14

        if (xx7(2) < -eps_wxfer) then
           skiptest( 8) = .true.
           skiptest(11) = .true.
           skiptest(14) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 8

        if (xx7(2) <= sfcwater_mass) then
           skiptest(8) = .true.
        endif

     endif

     6 continue
     if (skiptest(6)) go to 7

! Test 6:  If some sfcwater is present, and given complete evaporation of any
!         veg_water, solve 7x7 system to test for zero or positive transpiration,
!         partial evaporation of sfcwater, and zero or positive evaporation from soil

     if (sfcwater_mass > 0.) then

        aa7(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
        aa7(1,2) =         a4 * h4
        aa7(1,3) =         a4 * h1
        aa7(1,4) = 0._r8
        aa7(1,5) =         a4 * h4
        aa7(1,6) =       - a4 * h4
        aa7(1,7) = 0._r8
        yy7(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water)  ! WVC row

        aa7(2,1) =         a7 * h4
        aa7(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa7(2,3) = 0._r8
        aa7(2,4) =         a7 * h2
        aa7(2,5) =         a7 * h4
        aa7(2,6) =       - a7 * h4
        aa7(2,7) = 0._r8
        yy7(2)   =         a7 * (y2 - h4 * veg_water)                ! WSC row

        aa7(3,1) =         a2 * h5 * alvveg
        aa7(3,2) = 0._r8
        aa7(3,3) = 1._r8 + a2 * (h5 + h7)
        aa7(3,4) =         a2 * h7
        aa7(3,5) =         a2 * h5 * alvveg
        aa7(3,6) = 0._r8
        aa7(3,7) =       - a2 * h7
        yy7(3)   =         a2 * (y4 - h5 * alvveg * veg_water)         ! HVC row

        aa7(4,1) = 0._r8
        aa7(4,2) =         a6 * h6 * alvsfc
        aa7(4,3) =         a6 * h7
        aa7(4,4) = 1._r8 + a6 * (h6 + h7)
        aa7(4,5) =         a6 * h6 * alvsfc
        aa7(4,6) = 0._r8
        aa7(4,7) =       - a6 * h7
        yy7(4)   =         a6 * y5                                   ! HSC row

        aa7(5,1) =         a8 * h4
        aa7(5,2) =         a8 * h4
        aa7(5,3) = 0._r8
        aa7(5,4) =         a8 * h3
        aa7(5,5) = 1._r8 + a8 * (h3 * alvsfc + h4)
        aa7(5,6) =       - a8 * h4
        aa7(5,7) = 0._r8
        yy7(5)   =         a8 * (y3 - h4 * veg_water)                ! WGC row

        aa7(6,1) =       - a9 * h8
        aa7(6,2) =       - a9 * h8
        aa7(6,3) = 0._r8
        aa7(6,4) = 0._r8
        aa7(6,5) =       - a9 * h8
        aa7(6,6) = 1._r8 + a9 * h8
        aa7(6,7) = 0._r8
        yy7(6)   =         a9 * y9   ! WCA row

        aa7(7,1) = 0._r8
        aa7(7,2) = 0._r8
        aa7(7,3) =       - a10 * h7
        aa7(7,4) =       - a10 * h7
        aa7(7,5) = 0._r8
        aa7(7,6) = 0._r8
        aa7(7,7) = 1._r8 + a10 * h7
        yy7(7)   =         a10 * y10 ! HCA row

        call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land6',7,aa7,yy7,glatw,glonw)

        if (xx7(1) > -eps_wxfer .and. &
            xx7(2) > -eps_wxfer .and. xx7(2) <= sfcwater_mass .and. &
            xx7(5) > -eps_wxfer) then ! test was successful

           wxfervc = veg_water
           wxfersc = xx7(2)
           hxfervc = xx7(3)
           hxfersc = xx7(4)
           wxfergc = xx7(5)
           wxferca = xx7(6)
           hxferca = xx7(7)
           transp  = xx7(1)
           itest = 6
           go to 20
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 9, 12, and 15

        if (xx7(2) < -eps_wxfer) then
           skiptest( 9) = .true.
           skiptest(12) = .true.
           skiptest(15) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 9

        if (xx7(2) <= sfcwater_mass) then
           skiptest(9) = .true.
        endif

     endif

     7 continue
     if (skiptest(7)) go to 8

! Test 7:  Given complete evaporation of any sfcwater, solve 6x6 system to test
!         for condensation onto vegetation and zero or positive evaporation from soil

     aa6(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa6(1,2) =         a1 * h4
     aa6(1,3) =         a1 * h1
     aa6(1,4) = 0._r8
     aa6(1,5) =       - a1 * h4
     aa6(1,6) = 0._r8
     yy6(1)   =         a1 * (y1 - h4 * sfcwater_mass)                ! WVC row

     aa6(2,1) =         a8 * h4
     aa6(2,2) = 1._r8 + a8 * (h3 * alvsfc + h4)
     aa6(2,3) = 0._r8
     aa6(2,4) =         a8 * h3
     aa6(5,5) =       - a8 * h4
     aa6(5,6) = 0._r8
     yy6(2)   =         a8 * (y3 - (h3 * alvsfc + h4) * sfcwater_mass)  ! WGC row

     aa6(3,1) =         a2 * h5 * alvveg
     aa6(3,2) = 0._r8
     aa6(3,3) = 1._r8 + a2 * (h5 + h7)
     aa6(3,4) =         a2 * h7
     aa6(3,5) = 0._r8
     aa6(3,6) =       - a2 * h7
     yy6(3)   =         a2 * y4                                       ! HVC row

     aa6(4,1) = 0._r8
     aa6(4,2) =         a6 * h6 * alvsfc
     aa6(4,3) =         a6 * h7
     aa6(4,4) = 1._r8 + a6 * (h6 + h7)
     aa6(4,5) = 0._r8
     aa6(4,6) =       - a6 * h7
     yy6(4)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)         ! HSC row

     aa6(5,1) =       - a9 * h8
     aa6(5,2) =       - a9 * h8
     aa6(5,3) = 0._r8
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a9 * h8
     aa6(5,6) = 0._r8
     yy6(5)   =         a9 * y9   ! WCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a10 * h7
     aa6(6,4) =       - a10 * h7
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a10 * h7
     yy6(6)   =         a10 * y10 ! HCA row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land7',6,aa6,yy6,glatw,glonw)

     if (xx6(1) <  eps_wxfer .and. &
         xx6(2) > -eps_wxfer) then ! test was successful

        wxfervc = xx6(1)
        wxfersc = sfcwater_mass
        hxfervc = xx6(3)
        hxfersc = xx6(4)
        wxfergc = xx6(2)
        wxferca = xx6(5)
        hxferca = xx6(6)
        transp  = 0.
        itest = 7
        go to 20
     endif

     ! If this test resulted in condensation onto vegetation,
     ! skip past tests 8 and 9

     if (xx6(1) < -eps_wxfer) then
        skiptest(8) = .true.
        skiptest(9) = .true.
     endif

     8 continue
     if (skiptest(8)) go to 9

! Test 8:  If vegetation is wet, and given complete evaporation of any sfcwater,
!         solve 4x4 system to test for partial evaporation from vegetation
!         and zero or positive evaporation from soil

     if (iwetveg) then

        aa6(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
        aa6(1,2) =         (a3 + a4) * h4
        aa6(1,3) =         (a3 + a4) * h1
        aa6(1,4) = 0._r8
        aa6(1,5) =       - (a3 + a4) * h4
        aa6(1,6) = 0._r8
        yy6(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass)         ! WVC row

        aa6(2,1) =         a8 * h4
        aa6(2,2) = 1._r8 + a8 * (h3 * alvsfc + h4)
        aa6(2,3) = 0._r8
        aa6(2,4) =         a8 * h3
        aa6(5,5) =       - a8 * h4
        aa6(5,6) = 0._r8
        yy6(2)   =         a8 * (y3 - (h3 * alvsfc + h4) * sfcwater_mass)  ! WGC row

        aa6(3,1) =         a2 * h5 * alvveg
        aa6(3,2) = 0._r8
        aa6(3,3) = 1._r8 + a2 * (h5 + h7)
        aa6(3,4) =         a2 * h7
        aa6(3,5) = 0._r8
        aa6(3,6) =       - a2 * h7
        yy6(3)   =         a2 * y4                                       ! HVC row

        aa6(4,1) = 0._r8
        aa6(4,2) =         a6 * h6 * alvsfc
        aa6(4,3) =         a6 * h7
        aa6(4,4) = 1._r8 + a6 * (h6 + h7)
        aa6(4,5) = 0._r8
        aa6(4,6) =       - a6 * h7
        yy6(4)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)         ! HSC row

        aa6(5,1) =       - a9 * h8
        aa6(5,2) =       - a9 * h8
        aa6(5,3) = 0._r8
        aa6(5,4) = 0._r8
        aa6(5,5) = 1._r8 + a9 * h8
        aa6(5,6) = 0._r8
        yy6(5)   =         a9 * y9   ! WCA row

        aa6(6,1) = 0._r8
        aa6(6,2) = 0._r8
        aa6(6,3) =       - a10 * h7
        aa6(6,4) =       - a10 * h7
        aa6(6,5) = 0._r8
        aa6(6,6) = 1._r8 + a10 * h7
        yy6(6)   =         a10 * y10 ! HCA row

        call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land8',6,aa6,yy6,glatw,glonw)

        evap = xx6(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx6(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx6(2) > -eps_wxfer) then ! test was successful

           wxfervc = evap
           wxfersc = sfcwater_mass
           hxfervc = xx6(3)
           hxfersc = xx6(4)
           wxfergc = xx6(2)
           wxferca = xx6(5)
           hxferca = xx6(6)
           transp  = xx6(1) - evap
           itest = 8
           go to 20
        endif

        ! If this test resulted in less than complete evaporation of veg_water,
        ! skip past test 9

        if (evap <= veg_water) then
           skiptest(9) = .true.
        endif

     endif

     9 continue
     if (skiptest(9)) go to 10

! Test 9:  Given complete evaporation of any sfcwater, and complete evaporation
!         of any veg_water, solve 6x6 system to test for zero or positive
!         transpiration and zero or positive evaporation from soil

     aa6(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa6(1,2) =         a4 * h4
     aa6(1,3) =         a4 * h1
     aa6(1,4) = 0._r8
     aa6(1,5) =       - a1 * h4
     aa6(1,6) = 0._r8
     yy6(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water - h4 * sfcwater_mass)  ! WVC row

     aa6(2,1) =         a8 * h4
     aa6(2,2) = 1._r8 + a8 * (h3 * alvsfc + h4)
     aa6(2,3) = 0._r8
     aa6(2,4) =         a8 * h3
     aa6(5,5) =       - a8 * h4
     aa6(5,6) = 0._r8
     yy6(2)   =         a8 * (y3 - h4 * veg_water - (h3 * alvsfc + h4) * sfcwater_mass)  ! WGC row

     aa6(3,1) =         a2 * h5 * alvveg
     aa6(3,2) = 0._r8
     aa6(3,3) = 1._r8 + a2 * (h5 + h7)
     aa6(3,4) =         a2 * h7
     aa6(3,5) = 0._r8
     aa6(3,6) =       - a2 * h7
     yy6(3)   =         a2 * (y4 - h5 * alvveg * veg_water)                              ! HVC row

     aa6(4,1) = 0._r8
     aa6(4,2) =         a6 * h6 * alvsfc
     aa6(4,3) =         a6 * h7
     aa6(4,4) = 1._r8 + a6 * (h6 + h7)
     aa6(4,5) = 0._r8
     aa6(4,6) =       - a6 * h7
     yy6(4)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)                          ! HSC row

     aa6(5,1) =       - a9 * h8
     aa6(5,2) =       - a9 * h8
     aa6(5,3) = 0._r8
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a9 * h8
     aa6(5,6) = 0._r8
     yy6(5)   =         a9 * y9   ! WCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a10 * h7
     aa6(6,4) =       - a10 * h7
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a10 * h7
     yy6(6)   =         a10 * y10 ! HCA row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land9',6,aa6,yy6,glatw,glonw)

     if (xx6(1) > -eps_wxfer .and. &
         xx6(2) > -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass
        hxfervc = xx6(3)
        hxfersc = xx6(4)
        wxfergc = xx6(2)
        transp  = xx6(1)
        wxferca = xx6(5)
        hxferca = xx6(6)
        itest = 9
        go to 20
     endif

     10 continue
     if (skiptest(10)) go to 11

! Test 10:  If some sfcwater is present, and given zero vapor flux with soil,
!          solve 6x6 system to test for condensation onto vegetation and
!          partial evaporation of sfcwater

     if (sfcwater_mass > 0.) then

        aa6(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
        aa6(1,2) =         a1 * h4
        aa6(1,3) =         a1 * h1
        aa6(1,4) = 0._r8
        aa6(1,5) =       - a1 * h4
        aa6(1,6) = 0._r8
        yy6(1)   =         a1 * y1  ! WVC row

        aa6(2,1) =         a7 * h4
        aa6(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa6(2,3) = 0._r8
        aa6(2,4) =         a7 * h2
        aa6(2,5) =       - a7 * h4
        aa6(2,6) = 0._r8
        yy6(2)   =         a7 * y2  ! WSC row

        aa6(3,1) =         a2 * h5 * alvveg
        aa6(3,2) = 0._r8
        aa6(3,3) = 1._r8 + a2 * (h5 + h7)
        aa6(3,4) =         a2 * h7
        aa6(3,5) = 0._r8
        aa6(3,6) =       - a2 * h7
        yy6(3)   =         a2 * y4  ! HVC row

        aa6(4,1) = 0._r8
        aa6(4,2) =         a6 * h6 * alvsfc
        aa6(4,3) =         a6 * h7
        aa6(4,4) = 1._r8 + a6 * (h6 + h7)
        aa6(4,5) = 0._r8
        aa6(4,6) =       - a6 * h7
        yy6(4)   =         a6 * y5  ! HSC row

        aa6(5,1) =       - a9 * h8
        aa6(5,2) =       - a9 * h8
        aa6(5,3) = 0._r8
        aa6(5,4) = 0._r8
        aa6(5,5) = 1._r8 + a9 * h8
        aa6(5,6) = 0._r8
        yy6(5)   =         a9 * y9   ! WCA row

        aa6(6,1) = 0._r8
        aa6(6,2) = 0._r8
        aa6(6,3) =       - a10 * h7
        aa6(6,4) =       - a10 * h7
        aa6(6,5) = 0._r8
        aa6(6,6) = 1._r8 + a10 * h7
        yy6(6)   =         a10 * y10 ! HCA row

        call matrix8_NxN(6,AA6,YY6,XX6,sing); if (sing) call sing_print(iwsfc,'land10',6,aa6,yy6,glatw,glonw)

        if (xx6(1) <  eps_wxfer .and. &
            xx6(2) > -eps_wxfer .and. xx6(2) < sfcwater_mass) then ! test was successful

           wxfervc = xx6(1)
           wxfersc = xx6(2)
           hxfervc = xx6(3)
           hxfersc = xx6(4)
           wxferca = xx6(5)
           hxferca = xx6(6)
           wxfergc = 0.
           transp  = 0.
           itest = 10
           go to 20
        endif

        ! If this test resulted in condensation onto vegetation,
        ! skip past tests 11 and 12

        if (xx6(1) < -eps_wxfer) then
           skiptest(11) = .true.
           skiptest(12) = .true.
        endif

     endif

     11 continue
     if (skiptest(11)) go to 12

! Test 11:  If some sfcwater is present and vegetation is wet, and given zero
!          vapor flux with soil, solve 6x6 system to test for partial
!          evaporation from vegetation and partial evaporation of sfcwater

     if (sfcwater_mass > 0. .and. iwetveg) then

        aa6(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
        aa6(1,2) =         (a3 + a4) * h4
        aa6(1,3) =         (a3 + a4) * h1
        aa6(1,4) = 0._r8
        aa6(1,5) =       - (a3 + a4) * h4
        aa6(1,6) = 0._r8
        yy6(1)   =         (a3 + a4) * y1  ! WVC row

        aa6(2,1) =         a7 * h4
        aa6(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa6(2,3) = 0._r8
        aa6(2,4) =         a7 * h2
        aa6(2,5) =       - a7 * h4
        aa6(2,6) = 0._r8
        yy6(2)   =         a7 * y2         ! WSC row

        aa6(3,1) =         a2 * h5 * alvveg
        aa6(3,2) = 0._r8
        aa6(3,3) = 1._r8 + a2 * (h5 + h7)
        aa6(3,4) =         a2 * h7
        aa6(3,5) = 0._r8
        aa6(3,6) =       - a2 * h7
        yy6(3)   =         a2 * y4         ! HVC row

        aa6(4,1) = 0._r8
        aa6(4,2) =         a6 * h6 * alvsfc
        aa6(4,3) =         a6 * h7
        aa6(4,4) = 1._r8 + a6 * (h6 + h7)
        aa6(4,5) = 0._r8
        aa6(4,6) =       - a6 * h7
        yy6(4)   =         a6 * y5         ! HSC row

        aa6(5,1) =       - a9 * h8
        aa6(5,2) =       - a9 * h8
        aa6(5,3) = 0._r8
        aa6(5,4) = 0._r8
        aa6(5,5) = 1._r8 + a9 * h8
        aa6(5,6) = 0._r8
        yy6(5)   =         a9 * y9   ! WCA row

        aa6(6,1) = 0._r8
        aa6(6,2) = 0._r8
        aa6(6,3) =       - a10 * h7
        aa6(6,4) =       - a10 * h7
        aa6(6,5) = 0._r8
        aa6(6,6) = 1._r8 + a10 * h7
        yy6(6)   =         a10 * y10 ! HCA row

        call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land11',6,aa6,yy6,glatw,glonw)

        evap = xx6(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx6(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx6(2) > -eps_wxfer .and. xx6(2) <= sfcwater_mass) then ! test was successful

           wxfervc = evap
           wxfersc = xx6(2)
           hxfervc = xx6(3)
           hxfersc = xx6(4)
           wxferca = xx6(5)
           hxferca = xx6(6)
           wxfergc = 0.
           transp  = xx6(1) - evap
           itest = 11
           go to 20
        endif

        ! If this test resulted in less than complete evaporation of veg_water,
        ! skip past test 12

        if (evap <= veg_water) then
           skiptest(12) = .true.
        endif

     endif

     12 continue
     if (skiptest(12)) go to 13

! Test 12:  If some sfcwater is present, given complete evaporation of any
!          veg_water, and given zero vapor flux with soil, solve 6x6 system
!          to test for zero or positive transpiration and partial evaporation of
!          any surface water

     if (sfcwater_mass > 0.) then

        aa6(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
        aa6(1,2) =         a4 * h4
        aa6(1,3) =         a4 * h1
        aa6(1,4) = 0._r8
        aa6(1,5) =       - a4 * h4
        aa6(1,6) = 0._r8
        yy6(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water)  ! WVC row

        aa6(2,1) =         a7 * h4
        aa6(2,2) = 1._r8 + a7 * (h2 * alvsfc + h4)
        aa6(2,3) = 0._r8
        aa6(2,4) =         a7 * h2
        aa6(2,5) =       - a7 * h4
        aa6(2,6) = 0._r8
        yy6(2)   =         a7 * (y2 - h4               * veg_water)  ! WSC row

        aa6(3,1) =         a2 * h5 * alvveg
        aa6(3,2) = 0._r8
        aa6(3,3) = 1._r8 + a2 * (h5 + h7)
        aa6(3,4) =         a2 * h7
        aa6(3,5) = 0._r8
        aa6(3,6) =       - a2 * h7
        yy6(3)   =         a2 * (y4 - h5 * alvveg        * veg_water)  ! HVC row

        aa6(4,1) = 0._r8
        aa6(4,2) =         a6 * h6 * alvsfc
        aa6(4,3) =         a6 * h7
        aa6(4,4) = 1._r8 + a6 * (h6 + h7)
        aa6(4,5) = 0._r8
        aa6(4,6) =       - a6 * h7
        yy6(4)   =         a6 * y5                                   ! HSC row

        aa6(5,1) =       - a9 * h8
        aa6(5,2) =       - a9 * h8
        aa6(5,3) = 0._r8
        aa6(5,4) = 0._r8
        aa6(5,5) = 1._r8 + a9 * h8
        aa6(5,6) = 0._r8
        yy6(5)   =         a9 * y9   ! WCA row

        aa6(6,1) = 0._r8
        aa6(6,2) = 0._r8
        aa6(6,3) =       - a10 * h7
        aa6(6,4) =       - a10 * h7
        aa6(6,5) = 0._r8
        aa6(6,6) = 1._r8 + a10 * h7
        yy6(6)   =         a10 * y10 ! HCA row

        call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land12',6,aa6,yy6,glatw,glonw)

        if (xx6(1) > -eps_wxfer .and. &
            xx6(2) > -eps_wxfer .and. xx6(2) <= sfcwater_mass) then ! test was successful

           wxfervc = veg_water
           wxfersc = xx6(2)
           hxfervc = xx6(3)
           hxfersc = xx6(4)
           wxferca = xx6(5)
           hxferca = xx6(6)
           wxfergc = 0.
           transp  = xx6(1)
           itest = 12
           go to 20
        endif

     endif

     13 continue
     if (skiptest(13)) go to 14

! Test 13:  Given zero vapor flux with soil and given complete evaporation of
!          any sfcwater, solve 5x5 system to test for condensation onto vegetation

     aa5(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa5(1,2) =         a1 * h1
     aa5(1,3) = 0._r8
     aa5(1,4) =       - a1 * h4
     aa5(1,5) = 0._r8
     yy5(1)   =         a1 * (y1 - h4 * sfcwater_mass)         ! WVC row

     aa5(2,1) =         a2 * h5 * alvveg
     aa5(2,2) = 1._r8 + a2 * (h5 + h7)
     aa5(2,3) =         a2 * h7
     aa5(2,4) = 0._r8
     aa5(2,5) =       - a2 * h7
     yy5(2)   =         a2 * y4                                ! HVC row

     aa5(3,1) = 0._r8
     aa5(3,2) =         a6 * h7
     aa5(3,3) = 1._r8 + a6 * (h6 + h7)
     aa5(3,4) = 0._r8
     aa5(3,5) =       - a6 * h7
     yy5(3)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)  ! HSC row

     aa5(4,1) =       - a9 * h8
     aa5(4,2) = 0._r8
     aa5(4,3) = 0._r8
     aa5(4,4) = 1._r8 + a9 * h8
     aa5(4,5) = 0._r8
     yy5(4)   =         a9 * y9   ! WCA row

     aa5(5,1) = 0._r8
     aa5(5,2) =       - a10 * h7
     aa5(5,3) =       - a10 * h7
     aa5(5,4) = 0._r8
     aa5(5,5) = 1._r8 + a10 * h7
     yy5(5)   =         a10 * y10 ! HCA row

     call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land13',5,aa5,yy5,glatw,glonw)

     if (xx5(1) < eps_wxfer) then ! test was successful

        wxfervc = xx5(1)
        wxfersc = sfcwater_mass
        hxfervc = xx5(2)
        hxfersc = xx5(3)
        wxferca = xx5(4)
        hxferca = xx5(5)
        wxfergc = 0.
        transp  = 0.
        itest = 13
        go to 20
     endif

     ! If this test resulted in condensation onto vegetation,
     ! skip past tests 14 and 15

     if (xx5(1) < -eps_wxfer) then
        skiptest(14) = .true.
        skiptest(15) = .true.
     endif

     14 continue

     if (skiptest(14)) go to 15

! Test 14:  If vegetation is wet, and given zero vapor flux with soil and
!          complete evaporation of any sfcwater, solve 5x5 system to test for
!          partial evaporation from vegetation

     if (iwetveg) then

        aa5(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
        aa5(1,2) =         (a3 + a4) * h1
        aa5(1,3) = 0._r8
        aa5(1,4) =       - (a3 + a4) * h4
        aa5(1,5) = 0._r8
        yy5(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass)  ! WVC row

        aa5(2,1) =         a2 * h5 * alvveg
        aa5(2,2) = 1._r8 + a2 * (h5 + h7)
        aa5(2,3) =         a2 * h7
        aa5(2,4) = 0._r8
        aa5(2,5) =       - a2 * h7
        yy5(2)   =         a2 * y4                                ! HVC row

        aa5(3,1) = 0._r8
        aa5(3,2) =         a6 * h7
        aa5(3,3) = 1._r8 + a6 * (h6 + h7)
        aa5(3,4) = 0._r8
        aa5(3,5) =       - a6 * h7
        yy5(3)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)  ! HSC row

        aa5(4,1) =       - a9 * h8
        aa5(4,2) = 0._r8
        aa5(4,3) = 0._r8
        aa5(4,4) = 1._r8 + a9 * h8
        aa5(4,5) = 0._r8
        yy5(4)   =         a9 * y9   ! WCA row

        aa5(5,1) = 0._r8
        aa5(5,2) =       - a10 * h7
        aa5(5,3) =       - a10 * h7
        aa5(5,4) = 0._r8
        aa5(5,5) = 1._r8 + a10 * h7
        yy5(5)   =         a10 * y10 ! HCA row

        call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land14',5,aa5,yy5,glatw,glonw)

        evap = xx5(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx5(1) > -eps_wxfer .and. evap <= veg_water) then ! test was successful

           wxfervc = evap
           wxfersc = sfcwater_mass
           hxfervc = xx5(2)
           hxfersc = xx5(3)
           wxferca = xx5(4)
           hxferca = xx5(5)
           wxfergc = 0.
           transp  = xx5(1) - evap
           itest = 14
           go to 20
        endif

        ! If this test resulted in less than complete evaporation of veg_water,
        ! skip past test 15

        if (evap <= veg_water) then
           skiptest(15) = .true.
        endif

     endif

     15 continue
     if (skiptest(15)) go to 16


! Test 15:  Given complete evaporation of any veg_water, complete evaporation
!          of any sfcwater, and zero vapor flux with soil, solve 3x3 system
!          to test for zero or positive transpiration

     aa5(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa5(1,2) =         a4 * h1
     aa5(1,3) = 0._r8
     aa5(1,4) =       - a4 * h4
     aa5(1,5) = 0._r8
     yy5(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water - h4 * sfcwater_mass)  ! WVC row

     aa5(2,1) =         a2 * h5 * alvveg
     aa5(2,2) = 1._r8 + a2 * (h5 + h7)
     aa5(2,3) =         a2 * h7
     aa5(2,4) = 0._r8
     aa5(2,5) =       - a2 * h7
     yy5(2)   =         a2 * (y4 - h5 * alvveg        * veg_water)              ! HVC row

     aa5(3,1) = 0._r8
     aa5(3,2) =         a6 * h7
     aa5(3,3) = 1._r8 + a6 * (h6 + h7)
     aa5(3,4) = 0._r8
     aa5(3,5) =       - a6 * h7
     yy5(3)   =         a6 * (y5 - h6 * alvsfc * sfcwater_mass)                 ! HSC row

     aa5(4,1) =       - a9 * h8
     aa5(4,2) = 0._r8
     aa5(4,3) = 0._r8
     aa5(4,4) = 1._r8 + a9 * h8
     aa5(4,5) = 0._r8
     yy5(4)   =         a9 * y9   ! WCA row

     aa5(5,1) = 0._r8
     aa5(5,2) =       - a10 * h7
     aa5(5,3) =       - a10 * h7
     aa5(5,4) = 0._r8
     aa5(5,5) = 1._r8 + a10 * h7
     yy5(5)   =         a10 * y10 ! HCA row

     call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land15',5,aa5,yy5,glatw,glonw)

     if (xx5(1) > -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass
        hxfervc = xx5(2)
        hxfersc = xx5(3)
        wxferca = xx5(4)
        hxferca = xx5(5)
        wxfergc = 0.
        transp  = xx5(1)
        itest = 15
        go to 20
     endif

     16 continue

     ! If we got here, none of the above tests was satisfied, an outcome that should not happen

     itest = 91

     print*, 'itest = 91 result obtained in leaf4_canopy ',iwsfc,glatw,glonw
     stop 'stop test 91 '

     20 continue

     ! Update components

     cantemp = cantemp + (hxfervc + hxfersc - hxferca) * hcapcani

     veg_water = max(0.,veg_water - wxfervc)

     ! The latent heat of sublimation, alvi, must be used in computing the following gain/loss of
     ! veg_energy and energy_per_m2 due to vapor flux.  (See explanation above in similar context.)

     veg_energy = veg_energy + radveg - hxfervc - wxfervc * alvi - transp * alvl

     call qwtk(veg_energy, veg_water, hcapveg, veg_temp, fracliqv)

     ! Energy_per_m2 may pertain to sfcwater only, soil only, or both combined.  In cases
     ! where it pertains to sfcwater only, wxfergc will be zero.

     energy_per_m2 = energy_per_m2 + radsfc - hxfersc - (wxfersc + wxfergc) * alvi

     canrrv = canrrv + (wxfervc + transp + wxfersc + wxfergc - wxferca) * canairi

     sfluxt = hxferca / dt_leaf
     sfluxr = wxferca / dt_leaf

     delw            = dslzi(nzg) * wxfergc * .001
     soil_water(nzg) = soil_water(nzg) - delw
     head      (nzg) = head      (nzg) - delw * head_slope(nzg)

     sfcwater_mass = sfcwater_mass - wxfersc

     specvol = .001
     if (wxfersc > 0. .and. sfcwater_mass > 1.e-3) specvol = sfcwater_depth / sfcwater_mass

     sfcwater_depth = sfcwater_depth - wxfersc * specvol

  endif  ! End of vegetation section

  ! If icomb flag is set to 1, meaning that any sfcwater mass, if present, was
  ! brought into thermal equilibrium with the top soil layer, perform a second
  ! thermal equilibration to correctly distribute changes in energy_per_m2 that
  ! were added in this subroutine.

  if (icomb == 1) then

     call sfcwater_soil_comb(iland, iwsfc,                                         &
                             soil_water(nzg),soil_energy(nzg),                     &
                             specifheat_drysoil(nzg), sfcwater_mass,energy_per_m2, &
                             sfcwater_tempk,sfcwater_fracliq)

  endif

  end subroutine canopy

!===============================================================================

  subroutine vegndvi(iland, iwsfc, timefac_ndvi, leaf_class, veg_height, &
                     veg_rough, veg_ndvip, veg_ndvif, veg_ndvic,         &
                     veg_tai, veg_lai, veg_fracarea, veg_albedo          )

  use leaf_coms, only: veg_frac, albv_brown, albv_green, sai, dead_frac,  &
                       fpar_max, veg_clump, glai_max, dfpardsr, fpar_min, &
                       sr_max, sr_min, tai_max

  implicit none

  integer, intent(in) :: iland        ! index of current land cell
  integer, intent(in) :: iwsfc        ! index of current SFC grid cell
  integer, intent(in) :: leaf_class ! leaf class

  real, intent(in)  :: timefac_ndvi ! frac of time from past to future NDVI obs
  real, intent(in)  :: veg_height   ! veg height [m]
  real, intent(in)  :: veg_ndvip    ! veg past ndvi (obs time)
  real, intent(in)  :: veg_ndvif    ! veg future ndvi (obs time)
  real, intent(out) :: veg_ndvic    ! veg current ndvi
  real, intent(out) :: veg_tai      ! veg total area index
  real, intent(out) :: veg_lai      ! veg leaf area index
  real, intent(out) :: veg_fracarea ! veg fractional area
  real, intent(out) :: veg_albedo   ! veg albedo
  real, intent(out) :: veg_rough    ! veg roughness height [m]

  ! Local parameters

  real, parameter :: bz         = .91       ! for computing veg roughness height
  real, parameter :: hz         = .0075     ! for computing veg roughness height
  real, parameter :: extinc_veg = .5        ! veg extinction coefficient for rad flux
  real, parameter :: fpcon      = -.3338082 ! for computing veg leaf area index

  ! Local variables

  real :: sr         ! simple ratio
  real :: fpar       ! fraction of photosynthetically-active radiation
  real :: dead_lai   ! dead-matter leaf area index
  real :: green_frac ! lai fraction of tai

  ! Update vegetation TAI, LAI, fractional area, albedo, and roughness

  ! Compute LAI, vegetation roughness, albedo, vegfrac from time-dependent NDVI

  if (tai_max(leaf_class) < .1) then

     veg_lai = 0.
     veg_tai = 0.
     veg_rough = 0.
     veg_albedo = 0.
     veg_fracarea = 0.
     veg_ndvic = 0.

  else

     ! Time-interpolate ndvi to get current value veg_ndvic5 for this land cell

     veg_ndvic = veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi

     ! Limit ndvi to prevent values > .99 to prevent division by zero.

     if (veg_ndvic > .99) veg_ndvic = .99

     ! Compute "simple ratio" and limit between sr_min and sr_max(leaf_class).

     sr = (1. + veg_ndvic) / (1. - veg_ndvic)

     if (sr < sr_min) sr = sr_min
     if (sr > sr_max(leaf_class)) sr = sr_max(leaf_class)

     ! Compute fpar

     fpar = fpar_min + (sr - sr_min) * dfpardsr(leaf_class)

     ! Compute green leaf area index (veg_lai), dead leaf area index (dead_lai),
     ! total area index (tai), and green fraction

     veg_lai = glai_max(leaf_class) * (veg_clump(leaf_class) * fpar / fpar_max  &
        + (1. - veg_clump(leaf_class)) * alog(1. - fpar) * fpcon)

     dead_lai = (glai_max(leaf_class) - veg_lai) * dead_frac(leaf_class)

     veg_tai = sai(leaf_class) + veg_lai + dead_lai
     green_frac = veg_lai / veg_tai

     ! Compute vegetation roughness height, albedo, and fractional area

     veg_rough = veg_height * (1. - bz * exp(-hz * veg_tai))  * 0.3  ! Bob's test of reduced veg_rough
     veg_albedo = albv_green(leaf_class) * green_frac  &
                + albv_brown(leaf_class) * (1. - green_frac)
     veg_fracarea = veg_frac(leaf_class) * (1. - exp(-extinc_veg * veg_tai))

  endif

  end subroutine vegndvi

!===============================================================================

  subroutine sing_print(iwsfc,ieqn,nsize,aa,yy,glatw,glonw)

  use consts_coms, only: r8
  use mem_sfcg,    only: itab_wsfc
  use mem_ijtabs,  only: itab_w

  implicit none

  integer,          intent(in) :: iwsfc, nsize
  character(len=*), intent(in) :: ieqn
  real(r8),         intent(in) :: aa(nsize,nsize), yy(nsize)
  real,             intent(in) :: glatw, glonw

  integer :: irow

  print*, 'singular matrix '
  print*, 'iwsfc, ieqn, nsize: ', iwsfc, ieqn, nsize
  print*, 'glatw, glonw: ', glatw, glonw
  print*, 'iwsfc_globe, nwatm ', itab_wsfc(iwsfc)%iwglobe, itab_wsfc(iwsfc)%nwatm
  print*, 'iwatm_globe ', itab_w( itab_wsfc(iwsfc)%iwatm( 1:itab_wsfc(iwsfc)%nwatm ) )%iwglobe
  print*, 'kwatm ', itab_wsfc(iwsfc)%kwatm( 1:itab_wsfc(iwsfc)%nwatm )

  do irow = 1,nsize
     write(*,'(a,i5,12e15.5)') 'row, aa, yy: ',irow,aa(irow,1:nsize),yy(irow)
  enddo

  stop 'stopping singular matrix '

  end subroutine sing_print

!===============================================================================

  subroutine fast_canopy(iland, iwsfc, nlsw1, icomb,                         &
                         sfcwat_nud,    sfctemp_nud,     fracliq_nud,        &
                         sfcwater_mass, sfcwater_energy, sfcwater_depth,     &
                         energy_per_m2, sfcwater_tempk,  sfcwater_fracliq,   &
                         soil_water,    soil_energy,     specifheat_drysoil, &
                         soil_tempk,    soil_fracliq                         )

  use leaf_coms,      only: dt_leaf
  use mem_land,       only: nzg
  use therm_lib,      only: qwtk
  use leaf4_surface,  only: sfcwater_soil_comb

  implicit none

  integer, intent(in)    :: iland      ! index of current land cell
  integer, intent(in)    :: iwsfc      ! index of current SFC grid cell
  integer, intent(in)    :: nlsw1      ! k index of top active sfcwater layer
  integer, intent(inout) :: icomb      ! implicit heat balance flag [0=no, 1=yes]

  real, intent(in) :: sfcwat_nud   ! clim net water mass input (P-ET) [kg/(m^2 s)]
  real, intent(in) :: sfctemp_nud  ! clim sfc temp [K]
  real, intent(in) :: fracliq_nud  ! clim sfc fracliq []

  real, intent(inout) :: sfcwater_mass    ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_energy  ! surface water energy [J/kg]
  real, intent(inout) :: sfcwater_depth   ! surface water depth [m]
  real, intent(inout) :: energy_per_m2    ! surface water energy [J/m^2]
  real, intent(inout) :: sfcwater_tempk   ! surface water temperature [K]
  real, intent(inout) :: sfcwater_fracliq ! fraction of sfc water in liquid phase

  real, intent(inout) :: soil_water        (nzg) ! soil water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy       (nzg) ! soil energy [J/m^3]
  real, intent(in)    :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]
  real, intent(in)    :: soil_tempk        (nzg) ! soil temp [K]
  real, intent(in)    :: soil_fracliq      (nzg) ! fraction of soil moisture in liquid phase

  ! Local variables

  real :: flux, ediff

! PLAN: nudge temp/energy of sfcwater and soil(nzg); add net water (PCP - E - T)
! to sfcwater and soil.  Do runoff ONLY in locations/times (months?) when there
! is net water loss, and do runoff only when there is surface water present.
! Do not let sfcwater exceed a prescribed upper limit (of 2 m depth water equiv).

  ! Apply net change (P - ET) to sfcwater

  sfcwater_mass = max(0.,min(2000.,sfcwater_mass + sfcwat_nud * dt_leaf))

  ! In FAST_CANOPY, always set icomb flag to 1.

  call sfcwater_soil_comb(iland, iwsfc,                  &
                          soil_water(nzg),soil_energy(nzg), &
                          specifheat_drysoil(nzg), sfcwater_mass,energy_per_m2, &
                          sfcwater_tempk, sfcwater_fracliq)

  icomb = 1

  ! Compute surface energy flux based on difference between model and nudging
  ! values of sfcwater_tempk and sfcwater_fracliq

  ediff = sfctemp_nud - sfcwater_tempk + 80. * (fracliq_nud - sfcwater_fracliq)

  flux = max(-500.,min(500.,20. * ediff))

  energy_per_m2 = energy_per_m2 + flux * dt_leaf

  ! With icomb flag set to 1 in FAST_CANOPY, meaning that any sfcwater mass was
  ! brought into thermal equilibrium with the top soil layer, perform a second
  ! thermal equilibration to correctly distribute changes in energy_per_m2 that
  ! were added in this subroutine.

  call sfcwater_soil_comb(iland, iwsfc,                                         &
                          soil_water(nzg),soil_energy(nzg),                     &
                          specifheat_drysoil(nzg), sfcwater_mass,energy_per_m2, &
                          sfcwater_tempk,sfcwater_fracliq)

  end subroutine fast_canopy

! Remaining issues:

! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Build lookup tables, especially for things with exponentials?

End module leaf4_canopy
