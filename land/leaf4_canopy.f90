Module leaf4_canopy

  use consts_coms, only: cliq

  ! Parameters from Bonan et al. (2018) for leaf heat capacity
  ! (see section 2.4 and table 1 in paper):

  real, parameter :: fwat    = 0.70  ! fraction of fresh biomass that is water
  real, parameter :: cvegdry = 1500. ! specific heat dry biomass (J/kg/K)
  real, parameter :: cvegwet = fwat / (1. - fwat) * cliq ! Bonan eq. 29
  real, parameter :: madry   = 0.25  ! leaf dry mass per area (range .067 - .25)

  private :: cliq

Contains

  subroutine canopy(iland,         iwsfc,         ktrans,             transp,         &
                    leaf_class,    can_depth,     rhos,               vels,           &
                    ustar,         vkhsfc,        sfluxt,             sfluxr,         &
                    pcpg,          qpcpg,         dpcpg,              rshort,         &
                    cantemp,       canrrv,        glatw,              glonw,          &
                    airtheta,      airrrv,        canexner,           pbl_height,     &
                    snowfac,       vf,            stom_resist,        veg_height,     &
                    veg_rough,     veg_tai,       veg_lai,            hcapveg,        &
                    rshort_s,      rshort_v,      rlong_s,            rlong_v,        &
                    veg_water,     veg_energy,    veg_temp,           skncomp,        &
                    sfcwater_mass, sfcwater_epm2, sfcwater_depth,     soil_wfrac,     &
                    soil_water,    soil_energy,   wsat_vg,            wresid_vg,      &
                    soilfldcap,    ksat_vg,       specifheat_drysoil, head,           &
                    head_slope,    soil_tempk,    soil_fracliq,       thermcond_soil, &
                    sfcwater_tcond)

  use leaf_coms,     only: soil_rough, dt_leaf, kroot, rcmin, wcap_min, wcap_vmin, &
                           emisv, snow_rough
  use mem_sfcg,      only: sfcg, itab_wsfc
  use mem_land,      only: nzg, nzs_max, dslz, dslzi, dslzo2
  use leaf4_surface, only: sfcwater_soil_comb, grndvap_ab
  use consts_coms,   only: cp, vonk, alvi, alvl, alli, cliq, cice, rvap, t00, &
                           r8, stefan, cice1000, cliq1000, alli1000
  use matrix,        only: matrix8_2x2, matrix8_3x3, matrix8_4x4, matrix8_NxN
  use therm_lib,     only: qwtk, rhovsil, eslf

  implicit none

  integer, intent(in)  :: iland      ! index of current land cell
  integer, intent(in)  :: iwsfc      ! index of current SFC grid cell
  integer, intent(out) :: ktrans     ! k index of soil layer supplying transp
  real,    intent(out) :: transp     ! transpiration xfer this LEAF timestep [kg_vap/m^2]
  integer, intent(in)  :: skncomp    ! surface skinlayer composition type []
  integer, intent(in)  :: leaf_class ! leaf class (vegetation class)

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
  real, intent(in)    :: pbl_height       ! boundary layer height [m]
  real, intent(in)    :: snowfac          ! fractional veg burial by snowcover
  real, intent(in)    :: vf               ! fractional coverage of non-buried part of veg
  real, intent(inout) :: stom_resist      ! veg stomatal resistance [s/m]
  real, intent(in)    :: veg_height       ! veg height [m]
  real, intent(in)    :: veg_rough        ! veg roughess height [m]
  real, intent(in)    :: veg_tai          ! veg total area index
  real, intent(in)    :: veg_lai          ! veg leaf area index
  real, intent(inout) :: hcapveg          ! veg heat capacity [J/(m^2 K)]
  real, intent(in)    :: rshort_s         ! s/w rad flux absorbed by surface [W/m^2]
  real, intent(in)    :: rshort_v         ! s/w rad flux absorbed by veg [W/m^2]
  real, intent(in)    :: rlong_s          ! l/w rad flux absorbed by surface [W/m^2]
  real, intent(in)    :: rlong_v          ! l/w rad flux absorbed by veg [W/m^2]
  real, intent(inout) :: veg_water        ! veg sfc water content [kg/m^2]
  real, intent(inout) :: veg_energy       ! (veg + veg_water) energy [J/m^2]
  real, intent(inout) :: veg_temp         ! veg temp [K]

  real, intent(inout) :: sfcwater_mass (nzs_max)  ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_epm2 (nzs_max)  ! surface water energy [J/m^2]
  real, intent(inout) :: sfcwater_depth(nzs_max)  ! surface water depth [m]
  real, intent(in)    :: sfcwater_tcond(nzs_max)  ! surface water thermal conductivity [W/m/K]

  real, intent(inout) :: soil_wfrac        (nzg)  ! soil water fraction []
  real, intent(inout) :: soil_water        (nzg)  ! soil water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy       (nzg)  ! soil energy [J/m^3]
  real, intent(in)    :: wsat_vg           (nzg)  ! saturation water content (porosity) []
  real, intent(in)    :: wresid_vg         (nzg)  ! residual water content []
  real, intent(in)    :: soilfldcap               ! top-layer soil field capacity []
  real, intent(in)    :: ksat_vg           (nzg)  ! saturation hydraulic conductivity [m/s]
  real, intent(in)    :: specifheat_drysoil(nzg)  ! specific heat of dry soil [J/(m^3 K)]
  real, intent(inout) :: head              (nzg)  ! hydraulic head [m] (relative to local topo datum)
  real, intent(in)    :: head_slope        (nzg)  ! d(head) / d(soil_water) [m]
  real, intent(inout) :: soil_tempk        (nzg)  ! soil temp [K]
  real, intent(inout) :: soil_fracliq      (nzg)  ! fraction of soil moisture in liquid phase
  real, intent(in)    :: thermcond_soil    (nzg)  ! soil thermal conductivity [W/m/K]

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
  real, parameter :: bswp = -6, sswp = 7.42e-6 ! for soil water potential

  real, parameter :: fcn = 0.75 ! Crank-Nicolson future time weight for canopy
                                ! turbulent flux balance
  ! Local variables

  integer :: k        ! vertical index over soil layers
  integer :: iveg     ! flag for exposed vegetation (0=no, 1=yes)
  logical :: iwetveg  ! flag for wet vegetation
  logical :: iwetsfc  ! flag for wet surface
  logical :: drygnd   ! flag for dry soil surface

  real :: alpha       ! "alpha" term in soil surface humidity formulation
  real :: beta        ! "beta" term in soil surface humidity formulation
  real :: hxfersc     ! sfc-to-can_air heat xfer this step [J/m^2]
  real :: wxfersc     ! sfc-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxferca     ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca     ! vapor xfer from can_air to atm this step [kg_vap/m^2]
  real :: wxfergc     ! gnd-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxfervc     ! veg-to-can_air heat xfer this step [J/m^2]
  real :: wxfervc     ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
  real :: hxferus     ! head xfer between skin layer and layer belowe [J/m^2]
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
  real :: hcapskn     ! Heat capacity of surface skinlayer (sfcwater and/or soil(nzg)  ) [J/(m^2 K)]
  real :: hcapund     ! Heat capacity of layer (snow or soil) under skinlayer [J/(m^2 K)]
  real :: hcapveg_new ! Heat capacity of vegetation [J/m^2/K]
  real :: hcapvegw    ! Heat capacity of vegetation + any veg_water [J/(m^2 K)]
  real :: atmair      ! Boundary layer air mass [kg/m^2]
  real :: hcapatm     ! Heat capacity of atm boundary layer

  real :: canairi     ! Inverse of canair
  real :: hcapcani    ! Inverse of hcapcan
  real :: hcapskni    ! Inverse of hcapskn
  real :: hcapundi    ! Inverse of hcapund
  real :: hcapvegwi   ! Inverse of hcapvegw
  real :: atmairi     ! Inverse of atmair
  real :: hcapatmi    ! Inverse of hcapatm

  real :: alvveg      ! Latent heat of vaporization at veg surface, adjusted for phase of veg_water [J/(kg)]
  real :: alvskn      ! Latent heat of vaporization at ground surface, adjusted for phase of sfcwater [J/(kg)]

  real :: radsfc      ! Radiation absorbed by surface [J/m^2]
  real :: radveg      ! Radiation absorbed by vegetation [J/m^2]

  real :: skn_epm2    ! skinlayer energy per m^2 [J/m^2]
  real :: skn_wmass   ! skinlayer water [kg/m^2]
  real :: skn_fracliq ! skinlayer liquid fraction []
  real :: skn_tempk   ! skinlayer temperature [K]
  real :: skn_rhovs   ! skinlayer saturation water vapor density [kg_vap/m^3]
  real :: skn_rhovsp  ! skn_rhovs derivative with respect to temperature
  real :: und_tempk   ! temperature of layer under the skinlayer [K]

  real :: gnd_rhov    ! ground sfc evaporative water vapor density [kg_vap/m^3]
  real :: gnd_rhovp   ! gnd_rhov derivative with respect to temperature

  real :: eps_wxfer   ! very small water transfer to accomodate truncation error
  real :: evap        ! amount evaporated from veg sfc [kg/m^2]

  integer :: nzs           ! top surface water layer
  real :: sfcwater_tempk   ! temperature of top surface water layer [K]
  real :: sfcwater_fracliq ! fraction of sfc water in liquid phase

  real :: sfcwater_fracarea ! fractional area coverage of surface water []
  real :: sfcwater_tot      ! total surface water mass            [kg/m^2]
  real :: gndwater_tot      ! soil water mass in top (skin) layer [kg/m^2]

  real :: specvol, vw_max, fact, rfactor(2)
  real :: tf, wshed2, dshed2, delw

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
  real(r8) :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11
  real(r8) :: y1, y2, y3, y4, y5, y9, y10, y11

  real(r8) :: aa3(3,3), xx3(3), yy3(3)  ! 3x3 matrix equation terms
  real(r8) :: aa4(4,4), xx4(4), yy4(4)  ! 4x4 matrix equation terms
  real(r8) :: aa5(5,5), xx5(5), yy5(5)  ! 5x5 matrix equation terms
  real(r8) :: aa6(6,6), xx6(6), yy6(6)  ! 6x6 matrix equation terms
  real(r8) :: aa7(7,7), xx7(7), yy7(7)  ! 7x7 matrix equation terms
  real(r8) :: aa8(8,8), xx8(8), yy8(8)  ! 8x8 matrix equation terms

  logical :: sing, skiptest(22), didtest(22)
  integer :: itest  ! test number

  integer, parameter :: zero8 = 0.0_r8

  didtest(:) = .false.

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

  ! Atmosphere air quantities

  atmair   = rhos * max( pbl_height, 3.0 * sfcg%dzt_bot(iwsfc) )
  atmairi  = 1. / atmair
  hcapatm  = cp * atmair
  hcapatmi = 1. / hcapatm

  if (iveg == 1) then  ! Vegetation case

     ! Vegetation is sufficiently abundant and not covered by snow.

     ! TAI and LAI reduced by ground snowcover (for low vegetation)

     stai = veg_tai * (1. - snowfac)
     slai = veg_lai * (1. - snowfac)

     ! Vegetation heat capacity

     hcapveg_new = madry * (cvegdry * stai + cvegwet * (slai + .08))
     veg_energy  = veg_energy + (hcapveg_new - hcapveg) * (veg_temp - 273.15)
     hcapveg     = hcapveg_new

     ! Maximum allowed vegetation surface water

     vw_max = 0.2 * stai

     ! COMPUTE RDI WITH VEGETATION INFLUENCE

     ! Compute ground-canopy resistance rd.
     ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
     ! the values (zognd,0) when veg covered.

     sfcwater_fracarea = min(1., sqrt( 0.5 * sfcwater_mass(1) ))
     zognd = soil_rough * (1.-sfcwater_fracarea) + snow_rough * sfcwater_fracarea
     zoveg = max( veg_rough * (1. - snowfac), zognd )
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
!    wtveg  = max(0.,min(1., 1.1 * veg_tai / covr))

     ! This is for just the vegetated part (vf)

     wtveg  = max(0.,min(1., 1.1 * stai / (vf*covr)))
     rdi    = ustar / (5. * (1. - wtveg) + ustar * rasveg * wtveg)

     ! Now averaged conductances for vegetated and non-vegetated fraction

     rdi = 0.2 * ustar * (1. - vf) + vf * rdi

     ! Add any precipitation intercepted mass and energy to vegetation surface

     veg_water  = veg_water  + pcpg  * vf
     veg_energy = veg_energy + qpcpg * vf

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

        rc     = stom_resist
        e_leaf = (rb * esat_veg + rc * e_can) / (rb + rc)
        vpd    = max(0., esat_veg - e_leaf)

        ! Evaluate 5 environmental factors and new rc
        ! (Swp multiplier converts from meters to Pascals (hydrostatic eqn for water)

        swp = max(swp, -570.)  ! limit fswp to about 1.e18
        vpd = min(vpd, 13.e3)  ! limit fvpd to about 1.e18

        ftlo = 1. + exp(-stlo * (veg_temp - btlo))
        fthi = 1. + exp(-sthi * (veg_temp - bthi))
        frad = 1. + exp(-srad * (rshort   - brad))
        fswp = 1. + exp(-sswp * (swp*9810.- bswp))
        fvpd = 1. + exp(-svpd * (vpd      - bvpd))

        ! Compute asymptotic value of stomatal resistance based on environmental factors

        rc_inf = ftlo * fthi * frad * min(fvpd * fswp, 1.e29) * rcmin(leaf_class)

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

     ktrans      = 0
     stom_resist = 1.e18

  endif

  ! For skncomp = 1, thermally combine sfcwater(1) and soil(nzg) layers,
  ! which also provides sfcwater_tempk and sfcwater_fracliq.

  ! For skncomp = 2 or 3, just diagnose sfcwater_tempk and sfcwater_fracliq.
  ! Use qwtk instead of qtk because sfcwater_epm2(1) is used
  ! ("dryhcap" = 10 is very small value)

  ! Diagnose surface skinlayer heat capacity (hcapskn) based on skncomp value.

  nzs = 1
  if (skncomp == 3) nzs = 2

  sfcwater_mass (nzs) = sfcwater_mass (nzs) +  pcpg * (1. - vf) + wshed
  sfcwater_depth(nzs) = sfcwater_depth(nzs) + dpcpg * (1. - vf) + dshed
  sfcwater_epm2 (nzs) = sfcwater_epm2 (nzs) + qpcpg * (1. - vf) + qshed

  ! Use the following latent heat for vapor flux at the skinlayer when multiplying
  ! inverse skinlayer heat capacity to get implicit temperature change of the
  ! skinlayer.  (When skinlayer water is a mixture of liquid and ice, inverse skinlayer
  ! heat capacity is zero, so alvskn has no effect.)

  if (skncomp == 1) then

     call sfcwater_soil_comb( iland, iwsfc, soil_water(nzg), soil_energy(nzg), &
                              specifheat_drysoil(nzg), sfcwater_mass(1), sfcwater_epm2(1), &
                              sfcwater_tempk, sfcwater_fracliq )

     soil_tempk  (nzg) = sfcwater_tempk
     soil_fracliq(nzg) = sfcwater_fracliq

     skn_tempk   = sfcwater_tempk
     skn_fracliq = sfcwater_fracliq
     skn_epm2    = sfcwater_epm2(1) + soil_energy(nzg) * dslz(nzg)
     skn_wmass   = sfcwater_mass(1) + soil_water (nzg) * dslz(nzg) * 1.e3

     if     (skn_epm2 < 0.) then
        hcapskn = skn_wmass * cice + specifheat_drysoil(nzg) * dslz(nzg)
        alvskn  = alvi - skn_fracliq * alli - (skn_tempk-273.15) * cice
     elseif (skn_epm2 > skn_wmass * alli) then
        hcapskn = skn_wmass * cliq + specifheat_drysoil(nzg) * dslz(nzg)
        alvskn  = alvi - skn_fracliq * alli - (skn_tempk-273.15) * cliq
     else
        hcapskn = skn_wmass * (alli / 1.) + specifheat_drysoil(nzg) * dslz(nzg)   ! Assume small positive dT/dE
        alvskn  = alvi - skn_fracliq * alli
     endif

     sfcwater_tot = sfcwater_mass(1)

     rfactor(2) = dslzo2(nzg)   / thermcond_soil(nzg)
     rfactor(1) = dslzo2(nzg-1) / thermcond_soil(nzg-1)

     und_tempk  = soil_tempk(nzg-1)

     if (soil_energy(nzg-1) < 0.) then
        hcapund = (soil_water(nzg-1) * cice1000 + specifheat_drysoil(nzg-1)) * dslz(nzg-1)
     elseif (soil_energy(nzg-1) > soil_water(nzg-1) * alli1000) then
        hcapund = (soil_water(nzg-1) * cliq1000 + specifheat_drysoil(nzg-1)) * dslz(nzg-1)
     else
        hcapund = (soil_water(nzg-1) * alli1000 + specifheat_drysoil(nzg-1)) * dslz(nzg-1)
     endif

  else ! skncomp = 2 or 3

     call qwtk( sfcwater_epm2(nzs), sfcwater_mass(nzs), 0., sfcwater_tempk, sfcwater_fracliq )

     skn_tempk   = sfcwater_tempk
     skn_fracliq = sfcwater_fracliq

     if     (sfcwater_epm2(nzs) < 0.) then
        hcapskn = sfcwater_mass(nzs) * cice
        alvskn  = alvi - skn_fracliq * alli - (skn_tempk-273.15) * cice
     elseif (sfcwater_epm2(nzs) > sfcwater_mass(nzs) * alli) then
        hcapskn = sfcwater_mass(nzs) * cliq
        alvskn  = alvi - skn_fracliq * alli - (skn_tempk-273.15) * cliq
     else
        hcapskn = sfcwater_mass(nzs) * (alli / 1.)  ! Assume small but nonzero dT/dE
        alvskn  = alvi - skn_fracliq * alli
     endif

     sfcwater_tot = sum( sfcwater_mass(1:nzs) )

     rfactor(2) = 0.5 * sfcwater_depth(nzs) / sfcwater_tcond(nzs)

     if (nzs > 1) then

        rfactor(1) = 0.5 * sfcwater_depth(nzs-1) / sfcwater_tcond(nzs-1)
        call qwtk( sfcwater_epm2(nzs-1), sfcwater_mass(nzs-1), 0., und_tempk, fact )

        if     (sfcwater_epm2(nzs-1) < 0.) then
           hcapund = max(10., sfcwater_mass(nzs-1)) * cice
        elseif (sfcwater_epm2(nzs-1) > sfcwater_mass(nzs-1) * alli) then
           hcapund = max(10., sfcwater_mass(nzs-1)) * cliq
        else
           hcapund = max(10., sfcwater_mass(nzs-1)) * (alli / 1.)  ! Assume small but nonzero dT/dE
        endif

     else

        rfactor(1) = dslzo2(nzg) / thermcond_soil(nzg)
        und_tempk = soil_tempk(nzg)

        if (soil_energy(nzg) < 0.) then       ! dslz(nzg-1) used to give increased storage capacity
           hcapund = soil_water(nzg) * dslz(nzg) * cice1000 + specifheat_drysoil(nzg) * dslz(nzg-1)
        elseif (soil_energy(nzg) > soil_water(nzg) * alli1000) then
           hcapund = soil_water(nzg) * dslz(nzg) * cliq1000 + specifheat_drysoil(nzg) * dslz(nzg-1)
        else
           hcapund = soil_water(nzg) * dslz(nzg) * alli1000 + specifheat_drysoil(nzg) * dslz(nzg-1)
        endif

     endif

  endif

  hcapskni = 1. / hcapskn
  hcapundi = 1. / hcapund

  ! Sfcwater saturation vapor density and derivative

  skn_rhovs  = rhovsil(skn_tempk - 273.15)
  skn_rhovsp = rhovsil(skn_tempk - 272.15) - skn_rhovs

  ! Reduced saturation vapor density and gradient of soil during evaporation

  drygnd = .false.

  if (skncomp == 1) then

     call grndvap_ab(iland, soil_tempk(nzg), soil_water(nzg), wresid_vg(nzg), soilfldcap, &
                     head(nzg), alpha, beta)

     gnd_rhov  = alpha * skn_rhovs
     gnd_rhovp = alpha * skn_rhovsp

     ! If soil is very dry, skip evaporation tests
     if (beta < .0001 .or. gnd_rhov < can_rhov) drygnd = .true.

  else

     ! In this case, sfcwater_fracarea will equal 1, and the gnd_rhov,
     ! gnd_rhovp, and beta values will just get multiplied by 0.

     beta      = 1.0
     gnd_rhov  = skn_rhovs
     gnd_rhovp = skn_rhovsp

  endif

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
  ! [radsfc is net shortwave + longwave absorbed by surface].

  ! Define sfcwater_fracarea to be 0.0 with zero sfcwater_mass, to be 1.0 when
  ! sfcwater_mass is at or above a threshold value (currently set to 2 kg/m^2),
  ! and to increase as the sqrt of sfcwater_mass below the threshold.
  ! Sfcwater_fracarea partitions surface evaporation into sfcwater and bare
  ! ground portions.  The sfcwater_mass threshold should be set large enough
  ! that sfcwater_mass cannot all evaporate in a single model time step.  The
  ! sqrt functional form (in which sfcwater_fracarea increases very rapidly
  ! at low sfcwater_mass) ensures that small amounts of sfcwater_mass will be
  ! able to evaporate in a timely manner.

  ! Sfcwater_mass(1) is used in the following tests as a proxy for surface water
  ! abundance, even though sfcwater_mass(2) may be the actual layer that exchanges
  ! vapor with canopy air.  This is permitted because if sfcwater(2) is active,
  ! then there is always enough water in it that total evaporation cannot occur in
  ! one timestep.  When there is little enough sfcwater_mass that all may evaporate
  ! in one timestep or that some ground will be exposed, then sfcwater(1) will be
  ! the layer that exchanges vapor with canopy air.

  if (sfcwater_mass(1) > wcap_min) then
     sfcwater_fracarea = sqrt( 0.5 * sfcwater_mass(1) )
     if (sfcwater_fracarea >= .999) sfcwater_fracarea = 1.0
     iwetsfc = .true.
  else
     sfcwater_fracarea = 0.
     iwetsfc = .false.
  endif

  radsfc = dt_leaf * (rshort_s + rlong_s)

  a5  = dt_leaf * rdi                          ! sfc vap xfer coef  (all area)
  a6  = cp * rhos * a5                         ! sfc heat xfer coef (all area)
  a7  = a5 * sfcwater_fracarea                 ! sfc vap xfer coef  (wet areas only)
  a8  = a5 * beta * (1. - sfcwater_fracarea)   ! sfc vap xfer coef  (soil areas only)
  a9  = dt_leaf * vkhsfc / sfcg%dzt_bot(iwsfc) ! can-atm vap xfer coef
  a10 = cp * a9                                ! can-atm heat xfer coef
  a11 = dt_leaf / (rfactor(1) + rfactor(2))    ! snow/soil to sfc heat xfer coef

  h2 = fcn * skn_rhovsp * hcapskni
  h3 = fcn * gnd_rhovp  * hcapskni
  h4 = fcn * rhos * canairi
  h6 = fcn * hcapskni
  h7 = fcn * hcapcani
  h8 = fcn * canairi

  ! Give soil/snow a finite (but large) heat capacity in implicit solver
  h9 = fcn * hcapundi

  ! Constant soil/snow temperature boundary condition (infinity heat capacity)
  ! h9 = 0.

  ! Give atmosphere heat/vapor a finite storage capacity proportional to the PBL depth
  h10 = fcn * atmairi * rhos  ! atm inverse water vapor storage
  h11 = fcn * hcapatmi        ! atm inverse heat capacity

  ! Keep atmosphere temperature/vapor constant in implicit solver (infinite storage capacity)
  ! h10 = 0.
  ! h11 = 0.

  y2  = skn_rhovs - can_rhov + h2 * radsfc
  y5  = skn_tempk - cantemp  + h6 * radsfc
  y3  = gnd_rhov  - can_rhov + h3 * radsfc

  y9  = canrrv    - airrrv
  y10 = cantemp   - canexner * airtheta
  y11 = und_tempk - skn_tempk - h6 * radsfc

  if (iveg == 0) then

     ! Case with NO VEGETATION

     ! Set all skiptest values to .true.

     skiptest(:) = .true.
     skiptest(1) = .false.

     ! Fill arrays for matrix solution.

! Test 101:  Solve 4x4 system to test for either condensation onto surface or
!            partial evaporation from sfcwater that completely covers the surface.
!            (When sfcwater completely covers the surface, it is too abundant to
!            completely evaporate in a single time step.)

     if (skiptest(1)) goto 102

     didtest(1) = .true.

     aa5(1,1) = 1._r8 + a5 * (h2 * alvskn + h4)
     aa5(1,2) =         a5 * h2
     aa5(1,3) =       - a5 * h4
     aa5(1,4) = 0._r8
     aa5(1,5) =       - a5 * h2
     yy5(1)   =         a5 * y2         ! WSC row

     aa5(2,1) =         a6 * h6 * alvskn
     aa5(2,2) = 1._r8 + a6 * (h6 + h7)
     aa5(2,3) = 0._r8
     aa5(2,4) =       - a6 * h7
     aa5(2,5) =       - a6 * h6
     yy5(2)   =         a6 * y5         ! HSC row

     aa5(3,1) =       - a9 * h8
     aa5(3,2) = 0._r8
     aa5(3,3) = 1._r8 + a9 * (h8 + h10)
     aa5(3,4) = 0._r8
     aa5(3,5) = 0._r8
     yy5(3)   =         a9 * y9         ! WCA row

     aa5(4,1) = 0._r8
     aa5(4,2) =       - a10 * h7
     aa5(4,3) = 0._r8
     aa5(4,4) = 1._r8 + a10 * (h7 + h11)
     aa5(4,5) = 0._r8
     yy5(4)   =         a10 * y10       ! HCA row

     aa5(5,1) =       - a11 * h6 * alvskn
     aa5(5,2) =       - a11 * h6
     aa5(5,3) = 0._r8
     aa5(5,4) = 0._r8
     aa5(5,5) = 1._r8 + a11 * (h6 + h9)
     yy5(5)   =         a11 * y11       ! HUS row

     call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land101',5,aa5,yy5,glatw,glonw)

     if (xx5(1) <= eps_wxfer .or. sfcwater_fracarea > 0.999) then ! test was successful

        wxfersc = min( real(xx5(1)), sfcwater_tot )
        hxfersc = xx5(2)
        wxferca = xx5(3)
        hxferca = xx5(4)
        wxfergc = 0.
        hxferus = xx5(5)
        itest   = 101

        go to 120
     endif

     ! Condensation/sfcwater evap test failed; continue to test 102

     skiptest(2) = .false.

     102 continue

! All remaining situations do NOT have complete coverage of surface by sfcwater

     if (sfcwater_fracarea > 0.999) then
        write(*,'(a,i10,2f9.2)') 'NOVEG CASE sfcwater_fracarea should be less than 1', &
            iwsfc,glatw,glonw
        print*, 'sf1 ',sfcwater_fracarea,eps_wxfer,sfcwater_mass(1)
        print*, 'sf2 ',xx4(1),xx4(2),xx4(3),xx4(4)
        print*, 'sf3 ',a5,a6,a9,a10
        print*, 'sf4 ',h2,h4,h6,h7,h8
        print*, 'sf5 ',y2,y5,y9,y10
        stop 'sfcwater_fracarea = 1, NOVEG CASE'
     endif

! All remaining tests also do NOT have condensation onto surface.

! Test 102:  If some sfcwater is present, solve 5x5 system to test for
!            partial evaporation from sfcwater and evaporation from soil

     ! If we are here and surface is not wet, do test 103 instead
     if ((.not. skiptest(2)) .and. (.not. iwetsfc)) then
        skiptest(2) = .true.
        skiptest(3) = .false.
     endif

     ! If we are here and soil is dry, do test 104 instead
     if ((.not. skiptest(2)) .and. drygnd) then
        skiptest(2) = .true.
        skiptest(4) = .false.
     endif

     if (skiptest(2)) goto 103

     didtest(2) = .true.

     aa6(1,1) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa6(1,2) =         a7 * h2
     aa6(1,3) =         a7 * h4
     aa6(1,4) =       - a7 * h4
     aa6(1,5) = 0._r8
     aa6(1,6) =       - a7 * h2
     yy6(1)   =         a7 * y2    ! WSC row

     aa6(2,1) =         a6 * h6 * alvskn
     aa6(2,2) = 1._r8 + a6 * (h6 + h7)
     aa6(2,3) =         a6 * h6 * alvskn
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a6 * h7
     aa6(2,6) =       - a6 * h6
     yy6(2)   =         a6 * y5    ! HSC row

     aa6(3,1) =         a8 * h4
     aa6(3,2) =         a8 * h3
     aa6(3,3) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa6(3,4) =       - a8 * h4
     aa6(3,5) = 0._r8
     aa6(3,6) =       - a8 * h3
     yy6(3)   =         a8 * y3    ! WGC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) =       - a9 * h8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * y9    ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) = 0._r8
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10  ! HCA row

     aa6(6,1) =       - a11 * h6 * alvskn
     aa6(6,2) =       - a11 * h6
     aa6(6,3) =       - a11 * h6 * alvskn
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * y11       ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land102',6,aa6,yy6,glatw,glonw)

     if ( xx6(1) >= -eps_wxfer .and. xx6(1) <= sfcwater_mass(1) .and. &
          xx6(3) >= -eps_wxfer ) then ! test was successful

        wxfersc = xx6(1)
        hxfersc = xx6(2)
        wxfergc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        itest   = 102

        go to 120
     endif

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 3

     if (xx6(1) > sfcwater_mass(1)) then

        skiptest(3) = .false.

     ! If this test resulted in condensation onto just the ground,
     ! go to test 104

     elseif (xx6(1) >= 0. .and. xx6(3) <= 0.) then

        skiptest(4) = .false.

     ! If this test resulted in condensation onto the ground and surface water,
     ! go to test 106

     else

        skiptest(6) = .false.

     endif

     103 continue

! Test 103: Given complete evaporation of any sfcwater, solve 4x4 system to
!           test for positive evaporation from soil

     ! If we are here and soil is dry, do test 105 instead
     if ((.not. skiptest(3)) .and. drygnd) then
        skiptest(3) = .true.
        skiptest(5) = .false.
     endif

     if (skiptest(3)) go to 104

     didtest(3) = .true.

     aa5(1,1) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa5(1,2) =         a8 * h3
     aa5(1,3) =       - a8 * h4
     aa5(1,4) = 0._r8
     aa5(1,5) =       - a8 * h3
     yy5(1)   =         a8 * (y3 - (h3 * alvskn + h4) * sfcwater_mass(1))  ! WGC row

     aa5(2,1) =         a6 * h6 * alvskn
     aa5(2,2) = 1._r8 + a6 * (h6 + h7)
     aa5(2,3) = 0._r8
     aa5(2,4) =       - a6 * h7
     aa5(2,5) =       - a6 * h6
     yy5(2)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))         ! HSC row

     aa5(3,1) =       - a9 * h8
     aa5(3,2) = 0._r8
     aa5(3,3) = 1._r8 + a9 * (h8 + h10)
     aa5(3,4) = 0._r8
     aa5(3,5) = 0._r8
     yy5(3)   =         a9 * (y9 + h8 * sfcwater_mass(1))                  ! WCA row

     aa5(4,1) = 0._r8
     aa5(4,2) =       - a10 * h7
     aa5(4,3) = 0._r8
     aa5(4,4) = 1._r8 + a10 * (h7 + h11)
     aa5(4,5) = 0._r8
     yy5(4)   =         a10 * y10                                          ! HCA row

     aa5(5,1) =       - a11 * h6 * alvskn
     aa5(5,2) =       - a11 * h6
     aa5(5,3) = 0._r8
     aa5(5,4) = 0._r8
     aa5(5,5) = 1._r8 + a11 * (h6 + h9)
     yy5(5)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))       ! HUS row

     call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land103',5,aa5,yy5,glatw,glonw)

     if (xx5(1) > -eps_wxfer) then ! test was successful

        wxfersc = sfcwater_mass(1)
        wxfergc = xx5(1)
        hxfersc = xx5(2)
        wxferca = xx5(3)
        hxferca = xx5(4)
        hxferus = xx5(5)
        itest   = 103

        go to 120
     endif

     ! If test failed, go to test 105

     skiptest(5) = .false.

     104 continue

! Test 104: If some sfcwater is present, and given zero vapor flux with soil,
!           solve 4x4 system to test for partial evaporation of sfcwater

     ! If we are here and surface is not wet, do test 105 instead
     if ((.not. skiptest(4)) .and. (.not. iwetsfc)) then
        skiptest(4) = .true.
        skiptest(5) = .false.
     endif

     if (skiptest(4)) go to 105

     didtest(4) = .true.

     aa5(1,1) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa5(1,2) =         a7 * h2
     aa5(1,3) =       - a7 * h4
     aa5(1,4) = 0._r8
     aa5(1,5) =       - a7 * h2
     yy5(1)   =         a7 * y2        ! WSC row

     aa5(2,1) =         a6 * h6 * alvskn
     aa5(2,2) = 1._r8 + a6 * (h6 + h7)
     aa5(2,3) = 0._r8
     aa5(2,4) =       - a6 * h7
     aa5(2,5) =       - a6 * h6
     yy5(2)   =         a6 * y5        ! HSC row

     aa5(3,1) =       - a9 * h8
     aa5(3,2) = 0._r8
     aa5(3,3) = 1._r8 + a9 * (h8 + h10)
     aa5(3,4) = 0._r8
     aa5(3,5) = 0._r8
     yy5(3)   =         a9 * y9        ! WCA row

     aa5(4,1) = 0._r8
     aa5(4,2) =       - a10 * h7
     aa5(4,3) = 0._r8
     aa5(4,4) = 1._r8 + a10 * (h7 + h11)
     aa5(4,5) = 0._r8
     yy5(4)   =         a10 * y10      ! HCA row

     aa5(5,1) =       - a11 * h6 * alvskn
     aa5(5,2) =       - a11 * h6
     aa5(5,3) = 0._r8
     aa5(5,4) = 0._r8
     aa5(5,5) = 1._r8 + a11 * (h6 + h9)
     yy5(5)   =         a11 * y11       ! HUS row

     call matrix8_NxN(5,aa5,yy5,xx5,sing); if (sing) call sing_print(iwsfc,'land104',5,aa5,yy5,glatw,glonw)

     if (xx5(1) >= -eps_wxfer .and. xx5(1) <= sfcwater_mass(1)) then ! test was successful

        wxfersc = xx5(1)
        hxfersc = xx5(2)
        wxferca = xx5(3)
        hxferca = xx5(4)
        hxferus = xx5(5)
        wxfergc = 0.
        itest   = 104

        go to 120
     endif

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 105

     if (xx5(1) > sfcwater_mass(1)) then

        skiptest(5) = .false.

     ! If this test resulted in condensation onto surface,
     ! go to test 106

     else

        skiptest(6) = .false.

     endif

     105 continue

! Test 105: Given zero vapor flux with soil and given complete evaporation of
!           any sfcwater, solve 3x3 system (a reduction of the 4x4 system in
!           Test 104) to get surface heat fluxes

     if (skiptest(5)) go to 106

     didtest(5) = .true.

     aa3(1,1) = 1._r8 + a6 * (h6 + h7)
     aa3(1,2) =       - a6 * h7
     aa3(1,3) =       - a6 * h6
     yy3(1)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))    ! HSC row

     aa3(2,1) =       - a10 * h7
     aa3(2,2) = 1._r8 + a10 * (h7 + h11)
     aa3(2,3) = 0._r8
     yy3(2)   =         a10 * y10                                     ! HCA row

     aa3(3,1) =       - a11 * h6
     aa3(3,2) = 0._r8
     aa3(3,3) = 1._r8 + a11 * (h6 + h9)
     yy3(3)   =         a11 * (y11 - h6 * alvskn * sfcwater_mass(1))  ! HUS row

     call matrix8_3x3(aa3,yy3,xx3,sing); if (sing) call sing_print(iwsfc,'land105',3,aa3,yy3,glatw,glonw)

     wxfersc = sfcwater_mass(1)
     wxfergc = 0.
     hxfersc = xx3(1)
     wxferca = a9 * (y9 + h8 * sfcwater_mass(1)) / (1._r8 + a9 * (h8 + h10))
     hxferca = xx3(2)
     hxferus = xx3(3)
     itest   = 105

     go to 120

     106 continue

! Test 106: Given zero vapor flux with soil and no evaporation of
!           sfcwater, solve 3x3 system (a reduction of the 4x4 system in
!           Test 104) to get surface heat flux

     if (skiptest(6)) go to 107

     didtest(6) = .true.

     aa3(1,1) = 1._r8 + a6 * (h6 + h7)
     aa3(1,2) =       - a6 * h7
     aa3(1,3) =       - a6 * h6
     yy3(1)   =         a6 * y5         ! HSC row

     aa3(2,1) =       - a10 * h7
     aa3(2,2) = 1._r8 + a10 * (h7 + h11)
     aa3(2,3) = 0._r8
     yy3(2)   =         a10 * y10       ! HCA row

     aa3(3,1) =       - a11 * h6
     aa3(3,2) = 0._r8
     aa3(3,3) = 1._r8 + a11 * (h6 + h9)
     yy3(3)   =         a11 * y11       ! HUS row

     call matrix8_3x3(aa3,yy3,xx3,sing); if (sing) call sing_print(iwsfc,'land106',3,aa3,yy3,glatw,glonw)

     wxfersc = 0.
     wxfergc = 0.
     hxfersc = xx3(1)
     wxferca = a9 * y9 / (1._r8 + a9 * (h8 + h10))
     hxferca = xx3(2)
     hxferus = xx3(3)
     itest   = 106

     go to 120

     107 continue

     ! If we got here, none of the above tests was satisfied, an outcome that should not happen

     itest = 92

     print*, 'itest = 92 result obtained in leaf4_canopy ',iwsfc,glatw,glonw
     stop 'stop test 92 '

     120 continue

     ! Update components

     ! The latent heat of sublimation, alvi, must be used in computing the gain/loss of
     ! sfcwater_epm2 due to vapor flux.  This is because extensive energy and mass carried by the
     ! vapor are both explicitly added to or subtracted from surface water, and the energy of vapor,
     ! defined relative to a zero-energy reference of ice at 0 deg C, is its mass times alvi
     ! (neglecting a small contribution from vapor temperature differing from 0 deg C).  Following
     ! the gain or loss of mass and energy via vapor flux, the resultant intensive energy of surface
     ! water is diagnosed as its extensive energy divided by its extensive mass.  If surface water
     ! is in the liquid phase, its intensive energy includes the latent heat of fusion.  Vapor
     ! carries this energy by virtue of its mass flux alone, so the total energy alvi carried by
     ! vapor is only alvl above this value, which correctly accounts for evaporation/condensation
     ! occurring in this case rather than sublimation/deposition.

     sfluxt = hxferca / dt_leaf
     sfluxr = wxferca / dt_leaf

     cantemp = cantemp + (hxfersc - hxferca) * hcapcani
     canrrv  = canrrv  + (wxfersc + wxfergc - wxferca) * canairi

     sfcwater_mass(nzs) = sfcwater_mass(nzs) - wxfersc
     sfcwater_epm2(nzs) = sfcwater_epm2(nzs) + radsfc - hxfersc - wxfersc * alvi + hxferus

     specvol = .001
     if (wxfersc > 0. .and. sfcwater_mass(nzs) > 1.e-3) specvol = sfcwater_depth(nzs) / sfcwater_mass(nzs)

     sfcwater_depth(nzs) = sfcwater_depth(nzs) - wxfersc * specvol

     if (skncomp == 1) then

        ! Sensible and latent heat contributions to the ground (soil) are added to
        ! sfcwater_epm2(1), but the subsequent call to subroutine sfcwater_soil_comb
        ! combines the energies of sfcwater(1) and soil(nzg)  .

        gndwater_tot     = max(0., soil_water(nzg) - wresid_vg(nzg)) * dslz(nzg) * 1.e3

        sfcwater_epm2(1) = sfcwater_epm2(1) - min(wxfergc, gndwater_tot) * alvi &
                                       + max(0., wxfergc - gndwater_tot) * alvl

        delw             = dslzi(nzg) * min(wxfergc, gndwater_tot) * .001
        fact             = wsat_vg(nzg) - wresid_vg(nzg)

        soil_water(nzg)  = soil_water(nzg) - delw
        head      (nzg)  = head      (nzg) - delw * head_slope(nzg)
        soil_wfrac(nzg)  = soil_wfrac(nzg) - delw / fact
        soil_wfrac(nzg)  = max(0., min(1., soil_wfrac(nzg)))

        ! With skncomp = 1, any sfcwater mass present was brought into thermal equilibrium
        ! with nzg soil layer.  Perform a second thermal equilibration to correctly
        ! distribute changes in sfcwater_epm2 that were added in this subroutine.

        call sfcwater_soil_comb( iland, iwsfc, soil_water(nzg), soil_energy(nzg), &
                                 specifheat_drysoil(nzg), sfcwater_mass(1), sfcwater_epm2(1), &
                                 sfcwater_tempk, sfcwater_fracliq )

        soil_tempk  (nzg) = sfcwater_tempk
        soil_fracliq(nzg) = sfcwater_fracliq

        soil_energy(nzg-1) = soil_energy(nzg-1) - hxferus * dslzi(nzg-1)

        call qwtk(soil_energy(nzg-1), soil_water(nzg-1)*1.e3, &
                  specifheat_drysoil(nzg-1), soil_tempk(nzg-1), soil_fracliq(nzg-1))

        ! In case too much water evaporated from the skin layer, try to remove it from layer nzg-1
        ! This should not be necessary but just in case...

        if (wxfergc - gndwater_tot > eps_wxfer) then
           delw              = dslzi(nzg-1) * (wxfergc-gndwater_tot) * .001
           delw              = min(delw, soil_water(nzg-1) - wresid_vg(nzg-1))
           fact              = wsat_vg(nzg) - wresid_vg(nzg-1)

           soil_energy(nzg-1) = soil_energy(nzg-1) + (specifheat_drysoil(nzg-1) &
                              * (soil_tempk(nzg-1) - t00) - soil_energy(nzg-1)) &
                              * delw / max(soil_water(nzg-1), wresid_vg(nzg-1))

           soil_water (nzg-1) = soil_water(nzg-1) - delw
           head       (nzg-1) = head      (nzg-1) - delw * head_slope(nzg-1)
           soil_wfrac (nzg-1) = soil_wfrac(nzg-1) - delw / fact
           soil_wfrac (nzg-1) = max(0., min(1., soil_wfrac(nzg-1)))
        endif

     else

        if (nzs > 1) then
           sfcwater_epm2(nzs-1) = sfcwater_epm2(nzs-1) - hxferus
        else
           soil_energy(nzg) = soil_energy(nzg) - hxferus * dslzi(nzg)
           call qwtk(soil_energy(nzg), soil_water(nzg)*1.e3, &
                     specifheat_drysoil(nzg), soil_tempk(nzg), soil_fracliq(nzg))
        endif

     endif

     transp = 0.

     ! set vegetation temperature to canopy temperature even if there is no
     ! vegetation or if it is snow covered

     veg_temp   = cantemp
     veg_water  = 0.0
     veg_energy = hcapveg * (veg_temp - 273.15)

  else

     ! Case WITH VEGETATION

     ! Use the following latent heat for vapor flux at the vegetation surface
     ! when multiplying inverse vegw heat capacity to get implicit temperature
     ! change of the vegetation.  (When veg_water is a mixture of liquid and ice,
     ! inverse veg heat capacity is zero, so alvveg has no effect.)

     alvveg = alvi - fracliqv * alli

     if     (veg_energy < 0.) then
        hcapvegw = veg_water * cice + hcapveg
        alvveg   = alvi - fracliqv * alli - tvegc * cice
     elseif (veg_energy > veg_water * alli) then
        hcapvegw = veg_water * cliq + hcapveg
        alvveg   = alvi - fracliqv * alli - tvegc * cliq
     else
        hcapvegw = veg_water * (alli / 1.) + hcapveg ! Assume small positive dT/dE
        alvveg   = alvi - fracliqv * alli
     endif

     hcapvegwi = 1. / hcapvegw

     radveg = dt_leaf * (rshort_v + rlong_v)

     ! veg_water fractional coverage.
     ! Check for presence of veg_water and set iwetveg flag

     if (veg_water > wcap_vmin) then
        sigmaw = ( veg_water / (.2 * stai) )**(2./3.)
        if (sigmaw >= .999) sigmaw = 1.0
        iwetveg = .true.
     else
        sigmaw  = 0.
        iwetveg = .false.
     endif

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

     skiptest(:) = .true.
     skiptest(1) = .false.

! Test 1: Solve 6x6 system to test for condensation onto vegetation and either
!         condensation onto surface or partial evaporation from sfcwater that
!         completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     if (skiptest(1)) goto 2

     didtest(1) = .true.

     aa7(1,1) = 1._r8 + a1 * (h1  * alvveg + h4)
     aa7(1,2) =         a1 * h4
     aa7(1,3) =         a1 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a1 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a1 * y1        ! WVC row

     aa7(2,1) =         a5 * h4
     aa7(2,2) = 1._r8 + a5 * (h2  * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a5 * h2
     aa7(2,5) =       - a5 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a5 * h2
     yy7(2)   =         a5 * y2        ! WSC row

     aa7(3,1) =         a2 * h5  * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5  + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4        ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6  * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6  + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5        ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * y9        ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10      ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11      ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land1',7,aa7,yy7,glatw,glonw)

     if ( ( (xx7(1) <= eps_wxfer) .or. (sigmaw            > 0.999 .and. xx7(1) <= veg_water) ) .and.  &
          ( (xx7(2) <= eps_wxfer) .or. (sfcwater_fracarea > 0.999) ) ) then ! test was sucessfull

        wxfervc = xx7(1)
        wxfersc = min( real(xx7(2)), sfcwater_tot )
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = 0.
        itest   = 1

        go to 91
     endif

     ! If surface is covered with water or there is evaporation from vegetation
     ! and condensation onto ground, go to test 2

     if ((sfcwater_fracarea > 0.999) .or. (xx7(1) >= 0. .and. xx7(2) <= 0.)) then

        skiptest(2) = .false.

     ! If there is evaporation from vegetation and ground, go to test 9

     elseif (xx7(1) >= 0. .and. xx7(2) >= 0.) then

        skiptest(9) = .false.

     ! If there is condensation onto vegetation and evaporation from surface,
     ! go to test 4

     else  ! xx7(1) < 0. .and. xx7(2) > 0.

        skiptest(4) = .false.

     endif

     2 continue

! Test 2: If vegetation is wet, solve 6x6 system to test for partial evaporation from
!         vegetation and either condensation onto surface or partial evaporation from
!         sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)


     ! If surface is not wet, go to test 3 instead
     if ((.not. skiptest(2)) .and. (.not. iwetveg)) then
        skiptest(2) = .true.
        skiptest(3) = .false.
     endif

     if (skiptest(2)) go to 3

     didtest(2) = .true.

     aa7(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa7(1,2) =         (a3 + a4) * h4
     aa7(1,3) =         (a3 + a4) * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - (a3 + a4) * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         (a3 + a4) * y1  ! WVC row

     aa7(2,1) =         a5 * h4
     aa7(2,2) = 1._r8 + a5 * (h2 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a5 * h2
     aa7(2,5) =       - a5 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a5 * h2
     yy7(2)   =         a5 * y2         ! WSC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4         ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5         ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * y9         ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10       ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11       ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land2',7,aa7,yy7,glatw,glonw)

     evap = xx7(1) * a3 / max(a3 + a4, 1.e-30_r8) ! evaporation of veg_water

     if ( xx7(1) >= -eps_wxfer .and. evap <= veg_water .and. &
         (xx7(2) <=  eps_wxfer .or. sfcwater_fracarea > 0.999) ) then ! test was successful

        wxfervc = evap
        wxfersc = min( real(xx7(2)), sfcwater_tot )
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = xx7(1) - evap
        itest   = 2

        go to 91
     endif

     ! If this test evaporates too much veg_water, go to test 3 and recheck

     if ( evap > veg_water ) then

        skiptest(3) = .false.

     ! If this test results in condensation onto vegetation and evaporation from sfc,
     ! go to test 19

     elseif (xx7(1) <= 0. .and. xx7(2) >= 0.) then

        skiptest(19) = .false.

     ! Go to test 13 if there is evaporation from ground and vegetation

     elseif (xx7(1) >= 0. .and. xx7(2) >= 0.) then

        skiptest(13) = .false.

     ! If this test results in condensation onto vegetation and ground,
     ! repeat test with no veg water transfer

     else  ! (xx7(1) < 0. .and. xx7(2) < 0.)

        a3 = 0._r8
        a4 = 0._r8
        goto 2

     endif

     3 continue

! Test 3: Given complete evaporation of any veg_water, solve 6x6 system to test
!         for transpiration and either condensation onto surface or partial
!         evaporation from sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     if (skiptest(3)) go to 4

     didtest(3) = .true.

     aa7(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa7(1,2) =         a4 * h4
     aa7(1,3) =         a4 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a4 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water) ! WVC row

     aa7(2,1) =         a5 * h4
     aa7(2,2) = 1._r8 + a5 * (h2 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a5 * h2
     aa7(2,5) =       - a5 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a5 * h2
     yy7(2)   =         a5 * (y2 - h4 * veg_water)                 ! WSC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * (y4 - h5 * alvveg * veg_water)        ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5                                    ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * (y9 + h8 * veg_water)                 ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10                                  ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11                                  ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land3',7,aa7,yy7,glatw,glonw)

     if (xx7(1) >= -eps_wxfer .and. (xx7(2) <= eps_wxfer .or. sfcwater_fracarea > 0.999)) then ! test was successful

        wxfervc = veg_water
        wxfersc = min( real(xx7(2)), sfcwater_tot )
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = xx7(1)
        itest   = 3

        go to 91
     endif

     ! If this test results in condensation onto vegetation and evaporation from sfc,
     ! go to test 20

     if (xx7(1) <= 0. .and. xx7(2) >= 0.) then

        skiptest(20) = .false.

     ! Go to test 14 if there is evaporation from ground and vegetation

     elseif (xx7(1) >= 0. .and. xx7(2) >= 0.) then

        skiptest(14) = .false.

     ! If this test results in condensation onto vegetation and ground,
     ! REPEAT test with NO vegetation water transfer

     else  ! (xx7(1) < 0. .and. xx7(2) < 0.)

        a4 = 0._r8
        goto 3

     endif

     4 continue

! All remaining situations (except the fallback tests 19-22) do NOT have
! complete coverage of surface by sfcwater

     if (sfcwater_fracarea > 0.999 .and. any(.not. skiptest(4:18))) then
        write(*,'(a,i10,2f9.2)') 'VEG CASE sfcwater_fracarea should be less than 1', &
            itab_wsfc(iwsfc)%iwglobe,glatw,glonw
        print*, 'sf11 ',sfcwater_fracarea,eps_wxfer,sfcwater_mass(1)
        print*, 'sf12 ',real( [xx6(1),xx6(2),xx6(3),xx6(4)] )
        print*, 'sf13 ',veg_water

        write(*,*) didtest
        write(*,*) skiptest
        write(*,*) cantemp, veg_temp, skn_tempk

        stop 'sfcwater_fracarea = 1, VEG CASE'
     endif

! All remaining tests also do NOT have condensation onto surface.

! Test 4: If some sfcwater is present, solve 7x7 system to test for
!         condensation onto vegetation, partial evaporation of sfcwater,
!         and evaporation from soil

     ! If we are here and surface is not wet, go to test 5
     if ((.not. skiptest(4)) .and. (.not. iwetsfc)) then
        skiptest(4) = .true.
        skiptest(5) = .false.
     endif

     ! If we are here and soil is dry, do test 6 instead
     if ((.not. skiptest(4)) .and. drygnd) then
        skiptest(4) = .true.
        skiptest(6) = .false.
     endif

     if (skiptest(4)) go to 5

     didtest(4) = .true.

     aa8(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa8(1,2) =         a1 * h4
     aa8(1,3) =         a1 * h1
     aa8(1,4) = 0._r8
     aa8(1,5) =         a1 * h4
     aa8(1,6) =       - a1 * h4
     aa8(1,7) = 0._r8
     aa8(1,8) = 0._r8
     yy8(1)   =         a1 * y1  ! WVC row

     aa8(2,1) =         a7 * h4
     aa8(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa8(2,3) = 0._r8
     aa8(2,4) =         a7 * h2
     aa8(2,5) =         a7 * h4
     aa8(2,6) =       - a7 * h4
     aa8(2,7) = 0._r8
     aa8(2,8) =       - a7 * h2
     yy8(2)   =         a7 * y2  ! WSC row

     aa8(3,1) =         a2 * h5 * alvveg
     aa8(3,2) = 0._r8
     aa8(3,3) = 1._r8 + a2 * (h5 + h7)
     aa8(3,4) =         a2 * h7
     aa8(3,5) = 0._r8
     aa8(3,6) = 0._r8
     aa8(3,7) =       - a2 * h7
     aa8(3,8) = 0._r8
     yy8(3)   =         a2 * y4  ! HVC row

     aa8(4,1) = 0._r8
     aa8(4,2) =         a6 * h6 * alvskn
     aa8(4,3) =         a6 * h7
     aa8(4,4) = 1._r8 + a6 * (h6 + h7)
     aa8(4,5) =         a6 * h6 * alvskn
     aa8(4,6) = 0._r8
     aa8(4,7) =       - a6 * h7
     aa8(4,8) =       - a6 * h6
     yy8(4)   =         a6 * y5  ! HSC row

     aa8(5,1) =         a8 * h4
     aa8(5,2) =         a8 * h4
     aa8(5,3) = 0._r8
     aa8(5,4) =         a8 * h3
     aa8(5,5) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa8(5,6) =       - a8 * h4
     aa8(5,7) = 0._r8
     aa8(5,8) =       - a8 * h3
     yy8(5)   =         a8 * y3  ! WGC row

     aa8(6,1) =       - a9 * h8
     aa8(6,2) =       - a9 * h8
     aa8(6,3) = 0._r8
     aa8(6,4) = 0._r8
     aa8(6,5) =       - a9 * h8
     aa8(6,6) = 1._r8 + a9 * (h8 + h10)
     aa8(6,7) = 0._r8
     aa8(6,8) = 0._r8
     yy8(6)   =         a9 * y9   ! WCA row

     aa8(7,1) = 0._r8
     aa8(7,2) = 0._r8
     aa8(7,3) =       - a10 * h7
     aa8(7,4) =       - a10 * h7
     aa8(7,5) = 0._r8
     aa8(7,6) = 0._r8
     aa8(7,7) = 1._r8 + a10 * (h7 + h11)
     aa8(7,8) = 0._r8
     yy8(7)   =         a10 * y10 ! HCA row

     aa8(8,1) = 0._r8
     aa8(8,2) =       - a11 * h6 * alvskn
     aa8(8,3) = 0._r8
     aa8(8,4) =       - a11 * h6
     aa8(8,5) =       - a11 * h6 * alvskn
     aa8(8,6) = 0._r8
     aa8(8,7) = 0._r8
     aa8(8,8) = 1._r8 + a11 * (h6 + h9)
     yy8(8)   =         a11 * y11      ! HUS row

     call matrix8_NxN(8,aa8,yy8,xx8,sing); if (sing) call sing_print(iwsfc,'land4',8,aa8,yy8,glatw,glonw)

     if (xx8(1) <=  eps_wxfer .and. &
         xx8(2) >= -eps_wxfer .and. xx8(2) <= sfcwater_mass(1) .and. &
         xx8(5) >= -eps_wxfer) then ! test was successful

        wxfervc = xx8(1)
        wxfersc = xx8(2)
        hxfervc = xx8(3)
        hxfersc = xx8(4)
        wxfergc = xx8(5)
        wxferca = xx8(6)
        hxferca = xx8(7)
        hxferus = xx8(8)
        transp  = 0.
        itest   = 4

        go to 91
     endif

     ! If this test resulted in too much evaporation of surface water,
     ! go to test 5 (evaporate remaining surface water) and recheck

     if ( xx8(2) > sfcwater_mass(1))  then

        skiptest(5) = .false.

     ! If this test resulted in evaporation from vegetation, ground, and sfc water
     ! go to test 9

     elseif ( xx8(1) > eps_wxfer .and. xx8(2) > eps_wxfer .and. xx8(5) > eps_wxfer ) then

        skiptest(9) = .false.

     ! If this test resulted in evaporation from vegetation and sfc water but
     ! condensation from ground, go to test 13

     elseif ( xx8(1) > eps_wxfer .and. xx8(2) > eps_wxfer ) then

        skiptest(13) = .false.

     ! If this test resulted in evaporation from vegetation and condensation onto
     ! surface water, go to test 17

     elseif ( xx8(1) > eps_wxfer ) then

        skiptest(17) = .false.

     ! If this test resulted in condensation onto vegetation, evaporation from
     ! surface water, and condensation onto the ground, go to test 6

     elseif ( xx8(1) <= eps_wxfer .and. xx8(2) >= -eps_wxfer .and. xx8(5) < -eps_wxfer) then

        skiptest(6) = .false.

     ! If this test resulted in condensation onto vegetation and sfc,
     ! go to test 8

     else

        skiptest(8) = .false.

     endif

     5 continue

! Test 5: Given complete evaporation of any sfcwater, solve 6x6 system to test
!         for condensation onto vegetation and zero or positive evaporation from soil

     ! If we are here and soil is dry, do test 7 instead
     if ((.not. skiptest(5)) .and. drygnd) then
        skiptest(5) = .true.
        skiptest(7) = .false.
     endif

     if (skiptest(5)) go to 6

     didtest(5) = .true.

     aa7(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa7(1,2) =         a1 * h4
     aa7(1,3) =         a1 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a1 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a1 * (y1 - h4 * sfcwater_mass(1))                 ! WVC row

     aa7(2,1) =         a8 * h4
     aa7(2,2) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a8 * h3
     aa7(2,5) =       - a8 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a8 * h3
     yy7(2)   =         a8 * (y3 - (h3 * alvskn + h4) * sfcwater_mass(1)) ! WGC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4                                           ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))        ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * (y9 + h8 * sfcwater_mass(1))                 ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10                                         ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))      ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land5',7,aa7,yy7,glatw,glonw)

     if (xx7(1) <=  eps_wxfer .and. &
         xx7(2) >= -eps_wxfer ) then ! test was successful

        wxfervc = xx7(1)
        wxfersc = sfcwater_mass(1)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxfergc = xx7(2)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        transp  = 0.
        itest   = 5

        go to 91
     endif

     ! If this test resulted in evaporation from vegetation and ground,
     ! go to test 11

     if ( xx7(1) > eps_wxfer .and. xx7(2) > eps_wxfer ) then

        skiptest(11) = .false.

     ! If this test resulted in evaporation from vegetation but
     ! condensation onto ground, go to test 15

     elseif ( xx7(1) > eps_wxfer .and. xx7(2) < eps_wxfer ) then

        skiptest(15) = .false.

     ! If this test resulted in condensation onto vegetation and ground,
     ! go to test 7

     else

        skiptest(7) = .false.

     endif

     6 continue

! Test 6:  If some sfcwater is present, and given zero vapor flux with soil,
!          solve 6x6 system to test for condensation onto vegetation and
!          partial evaporation of sfcwater

     ! If we are here and surface is not wet, do test 7 instead
     if ((.not. skiptest(6)) .and. (.not. iwetsfc)) then
        skiptest(6) = .true.
        skiptest(7) = .false.
     endif

     if (skiptest(6)) go to 7

     didtest(6) = .true.

     aa7(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa7(1,2) =         a1 * h4
     aa7(1,3) =         a1 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a1 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a1 * y1  ! WVC row

     aa7(2,1) =         a7 * h4
     aa7(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a7 * h2
     aa7(2,5) =       - a7 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a7 * h2
     yy7(2)   =         a7 * y2  ! WSC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4  ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5  ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * y9   ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10 ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11 ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land6',7,aa7,yy7,glatw,glonw)

     if ( xx7(1) <=  eps_wxfer .and. &
          xx7(2) >= -eps_wxfer .and. xx7(2) <= sfcwater_mass(1) ) then ! test was successful

        wxfervc = xx7(1)
        wxfersc = xx7(2)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = 0.
        itest   = 6

        go to 91
     endif

     ! If this test evaporated too much surface water, go to test 7
     ! and redo test with the available surface water evaporated

     if (xx7(2) > sfcwater_mass(1)) then

        skiptest(7) = .false.

     ! If this test resulted in evaporation from vegetation and sfc water
     ! go to test 13

     elseif ( xx7(1) > eps_wxfer .and. xx7(2) > eps_wxfer ) then

        skiptest(13) = .false.

     ! If this test resulted in evaporation from vegetation and condensation
     ! onto sfc water go to test 17

     elseif ( xx7(1) > eps_wxfer ) then

        skiptest(17) = .false.

     ! If this test resulted in condensation onto vegetation and sfc water
     ! go to test 8

     else

        skiptest(8) = .false.

     endif

     7 continue

! Test 7:  Given zero vapor flux with soil and given complete evaporation of
!          any sfcwater, solve 5x5 system to test for condensation onto vegetation

     if (skiptest(7)) go to 8

     didtest(7) = .true.

     aa6(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa6(1,2) =         a1 * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - a1 * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         a1 * (y1 - h4 * sfcwater_mass(1))            ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * y4                                      ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))   ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * (y9 + h8 * sfcwater_mass(1))            ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10                                    ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1)) ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land7',6,aa6,yy6,glatw,glonw)

     if (xx6(1) <= eps_wxfer) then ! test was successful

        wxfervc = xx6(1)
        wxfersc = sfcwater_mass(1)
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = 0.
        itest = 7

        go to 91
     endif

     ! This test resulted in evaporation from vegetation, go to test 15/16

     if (iwetveg) then
        skiptest(15) = .false.
     else
        skiptest(16) = .false.
     endif

     8 continue

! Test 8:  Given zero vapor flux with soil and no evaporation of
!          sfcwater, solve 5x5 system to test for condensation onto vegetation

     if (skiptest(8)) go to 9

     didtest(8) = .true.

     aa6(1,1) = 1._r8 + a1 * (h1 * alvveg + h4)
     aa6(1,2) =         a1 * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - a1 * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         a1 * y1          ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * y4          ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * y5          ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * y9          ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10        ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * y11        ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land8',6,aa6,yy6,glatw,glonw)

     if (xx6(1) <= eps_wxfer) then ! test was successful

        wxfervc = xx6(1)
        wxfersc = 0.
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = 0.
        itest = 8

        go to 91
     endif

     ! This test resulted in evaporation from vegetation, go to test 17

     skiptest(17) = .false.

     9 continue

! All remaining tests also do NOT have condensation onto vegetation.

! Test 9: If some sfcwater is present, and also if vegetation is wet, solve 7x7
!         system to test for partial evaporation from vegetation, partial
!         evaporation of sfcwater, and zero or positive evaporation from soil

     ! If we are here and soil is dry, do test 13 instead
     if ((.not. skiptest(9)) .and. drygnd) then
        skiptest(9)  = .true.
        skiptest(13) = .false.
     endif

     if (.not. skiptest(9)) then

        ! If vegetation is not wet, go to test 10
        if (.not. iwetveg) then
           skiptest(9)  = .true.
           skiptest(10) = .false.

        ! If surface is not wet, go to test 11
        elseif (.not. iwetsfc) then
           skiptest(9)  = .true.
           skiptest(11) = .false.
        endif

     endif

     if (skiptest(9)) go to 10

     didtest(9) = .true.

     aa8(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa8(1,2) =         (a3 + a4) * h4
     aa8(1,3) =         (a3 + a4) * h1
     aa8(1,4) = 0._r8
     aa8(1,5) =         (a3 + a4) * h4
     aa8(1,6) =       - (a3 + a4) * h4
     aa8(1,7) = 0._r8
     aa8(1,8) = 0._r8
     yy8(1)   =         (a3 + a4) * y1  ! WVC row

     aa8(2,1) =         a7 * h4
     aa8(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa8(2,3) = 0._r8
     aa8(2,4) =         a7 * h2
     aa8(2,5) =         a7 * h4
     aa8(2,6) =       - a7 * h4
     aa8(2,7) = 0._r8
     aa8(2,8) =       - a7 * h2
     yy8(2)   =         a7 * y2         ! WSC row

     aa8(3,1) =         a2 * h5 * alvveg
     aa8(3,2) = 0._r8
     aa8(3,3) = 1._r8 + a2 * (h5 + h7)
     aa8(3,4) =         a2 * h7
     aa8(3,5) = 0._r8
     aa8(3,6) = 0._r8
     aa8(3,7) =       - a2 * h7
     aa8(3,8) = 0._r8
     yy8(3)   =         a2 * y4         ! HVC row

     aa8(4,1) = 0._r8
     aa8(4,2) =         a6 * h6 * alvskn
     aa8(4,3) =         a6 * h7
     aa8(4,4) = 1._r8 + a6 * (h6 + h7)
     aa8(4,5) =         a6 * h6 * alvskn
     aa8(4,6) = 0._r8
     aa8(4,7) =       - a6 * h7
     aa8(4,8) =       - a6 * h6
     yy8(4)   =         a6 * y5         ! HSC row

     aa8(5,1) =         a8 * h4
     aa8(5,2) =         a8 * h4
     aa8(5,3) = 0._r8
     aa8(5,4) =         a8 * h3
     aa8(5,5) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa8(5,6) =       - a8 * h4
     aa8(5,7) = 0._r8
     aa8(5,8) =       - a8 * h3
     yy8(5)   =         a8 * y3         ! WGC row

     aa8(6,1) =       - a9 * h8
     aa8(6,2) =       - a9 * h8
     aa8(6,3) = 0._r8
     aa8(6,4) = 0._r8
     aa8(6,5) =       - a9 * h8
     aa8(6,6) = 1._r8 + a9 * (h8 + h10)
     aa8(6,7) = 0._r8
     aa8(6,8) = 0._r8
     yy8(6)   =         a9 * y9         ! WCA row

     aa8(7,1) = 0._r8
     aa8(7,2) = 0._r8
     aa8(7,3) =       - a10 * h7
     aa8(7,4) =       - a10 * h7
     aa8(7,5) = 0._r8
     aa8(7,6) = 0._r8
     aa8(7,7) = 1._r8 + a10 * (h7 + h11)
     aa8(7,8) = 0._r8
     yy8(7)   =         a10 * y10      ! HCA row

     aa8(8,1) = 0._r8
     aa8(8,2) =       - a11 * h6 * alvskn
     aa8(8,3) = 0._r8
     aa8(8,4) =       - a11 * h6
     aa8(8,5) =       - a11 * h6 * alvskn
     aa8(8,6) = 0._r8
     aa8(8,7) = 0._r8
     aa8(8,8) = 1._r8 + a11 * (h6 + h9)
     yy8(8)   =         a11 * y11      ! HUS row

     call matrix8_NxN(8,aa8,yy8,xx8,sing); if (sing) call sing_print(iwsfc,'land9',8,aa8,yy8,glatw,glonw)

     evap = xx8(1) * a3 / max(a3 + a4, 1.e-30_r8) ! evaporation of veg_water

     if ( xx8(1) >= -eps_wxfer .and. evap   <= veg_water        .and. &
          xx8(2) >= -eps_wxfer .and. xx8(2) <= sfcwater_mass(1) .and. &
          xx8(5) >= -eps_wxfer ) then  ! test was successful

        wxfervc = evap
        wxfersc = xx8(2)
        hxfervc = xx8(3)
        hxfersc = xx8(4)
        wxfergc = xx8(5)
        wxferca = xx8(6)
        hxferca = xx8(7)
        hxferus = xx8(8)
        transp  = xx8(1) - evap
        itest   = 9

        go to 91
     endif

     ! If this test resulted in more than complete evaporation of veg_water
     ! and more than complete evaporation of sfc water, go to test 12

     if (evap > veg_water .and. xx8(2) > sfcwater_mass(1)) then

        skiptest(12) = .false.

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 11

     elseif (xx8(2) > sfcwater_mass(1)) then

        skiptest(11) = .false.

     ! If this test resulted in more than complete evaporation of veg_water,
     ! go to test 10

     elseif (evap > veg_water) then

        skiptest(10) = .false.

     ! If this test resulted in condensation onto the ground sfc with transpiration
     ! and sfcwater evaporation, go to test 13

     elseif (xx8(1) >= -eps_wxfer .and. xx8(2) >= -eps_wxfer .and. xx8(5) < -eps_wxfer) then

        skiptest(13) = .false.

     ! If this test resulted in condensation onto the ground sfc and surface water
     ! with transpiration, go to test 17

     elseif (xx8(1) >= -eps_wxfer) then

        skiptest(17) = .false.

     ! If this test resulted in condensation onto vegetation and evaporation from
     ! surface water and ground, REPEAT test with no veg evaporation/transpiration

     elseif ( xx8(2) > eps_wxfer .and. xx8(5) > eps_wxfer) then

        a3 = 0._r8
        a4 = 0._r8
        go to 9

     ! If this test resulted in condensation onto vegetation and ground surface
     ! and evaporation from surface water, go to test 13 with no veg
     ! evaporation/transpiration

     elseif ( xx8(2) > eps_wxfer) then

        a3 = 0._r8
        a4 = 0._r8
        skiptest(13) = .false.

     ! If this test resulted in condensation onto vegetation, surface water, and
     ! and ground surface then go to test 19

     else

        skiptest(19) = .false.

     endif

     10 continue

! Test 10: If some sfcwater is present, and given complete evaporation of any
!          veg_water, solve 7x7 system to test for zero or positive transpiration,
!          partial evaporation of sfcwater, and zero or positive evaporation from soil

     ! If we are here and soil is dry, do test 14 instead
     if ((.not. skiptest(10)) .and. drygnd) then
        skiptest(10)  = .true.
        skiptest(14) = .false.
     endif

     ! If we are here and surface is not wet, do test 12 instead
     if ((.not. skiptest(10)) .and. (.not. iwetsfc)) then
        skiptest(10) = .true.
        skiptest(12) = .false.
     endif

     if (skiptest(10)) go to 11

     didtest(10) = .true.

     aa8(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa8(1,2) =         a4 * h4
     aa8(1,3) =         a4 * h1
     aa8(1,4) = 0._r8
     aa8(1,5) =         a4 * h4
     aa8(1,6) =       - a4 * h4
     aa8(1,7) = 0._r8
     aa8(1,8) = 0._r8
     yy8(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water)  ! WVC row

     aa8(2,1) =         a7 * h4
     aa8(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa8(2,3) = 0._r8
     aa8(2,4) =         a7 * h2
     aa8(2,5) =         a7 * h4
     aa8(2,6) =       - a7 * h4
     aa8(2,7) = 0._r8
     aa8(2,8) =       - a7 * h2
     yy8(2)   =         a7 * (y2 - h4 * veg_water)                  ! WSC row

     aa8(3,1) =         a2 * h5 * alvveg
     aa8(3,2) = 0._r8
     aa8(3,3) = 1._r8 + a2 * (h5 + h7)
     aa8(3,4) =         a2 * h7
     aa8(3,5) = 0._r8
     aa8(3,6) = 0._r8
     aa8(3,7) =       - a2 * h7
     aa8(3,8) = 0._r8
     yy8(3)   =         a2 * (y4 - h5 * alvveg * veg_water)         ! HVC row

     aa8(4,1) = 0._r8
     aa8(4,2) =         a6 * h6 * alvskn
     aa8(4,3) =         a6 * h7
     aa8(4,4) = 1._r8 + a6 * (h6 + h7)
     aa8(4,5) =         a6 * h6 * alvskn
     aa8(4,6) = 0._r8
     aa8(4,7) =       - a6 * h7
     aa8(4,8) =       - a6 * h6
     yy8(4)   =         a6 * y5                                     ! HSC row

     aa8(5,1) =         a8 * h4
     aa8(5,2) =         a8 * h4
     aa8(5,3) = 0._r8
     aa8(5,4) =         a8 * h3
     aa8(5,5) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa8(5,6) =       - a8 * h4
     aa8(5,7) = 0._r8
     aa8(5,8) =       - a8 * h3
     yy8(5)   =         a8 * (y3 - h4 * veg_water)                  ! WGC row

     aa8(6,1) =       - a9 * h8
     aa8(6,2) =       - a9 * h8
     aa8(6,3) = 0._r8
     aa8(6,4) = 0._r8
     aa8(6,5) =       - a9 * h8
     aa8(6,6) = 1._r8 + a9 * (h8 + h10)
     aa8(6,7) = 0._r8
     aa8(6,8) = 0._r8
     yy8(6)   =         a9 * (y9 + h8 * veg_water)                  ! WCA row

     aa8(7,1) = 0._r8
     aa8(7,2) = 0._r8
     aa8(7,3) =       - a10 * h7
     aa8(7,4) =       - a10 * h7
     aa8(7,5) = 0._r8
     aa8(7,6) = 0._r8
     aa8(7,7) = 1._r8 + a10 * (h7 + h11)
     aa8(7,8) = 0._r8
     yy8(7)   =         a10 * y10                                   ! HCA row

     aa8(8,1) = 0._r8
     aa8(8,2) =       - a11 * h6 * alvskn
     aa8(8,3) = 0._r8
     aa8(8,4) =       - a11 * h6
     aa8(8,5) =       - a11 * h6 * alvskn
     aa8(8,6) = 0._r8
     aa8(8,7) = 0._r8
     aa8(8,8) = 1._r8 + a11 * (h6 + h9)
     yy8(8)   =         a11 * y11                                   ! HUS row

     call matrix8_NxN(8,aa8,yy8,xx8,sing); if (sing) call sing_print(iwsfc,'land10',8,aa8,yy8,glatw,glonw)

     if ( xx8(1) >= -eps_wxfer .and. &
          xx8(2) >= -eps_wxfer .and. xx8(2) <= sfcwater_mass(1) .and. &
          xx8(5) >= -eps_wxfer ) then ! test was successful

        wxfervc = veg_water
        wxfersc = xx8(2)
        hxfervc = xx8(3)
        hxfersc = xx8(4)
        wxfergc = xx8(5)
        wxferca = xx8(6)
        hxferca = xx8(7)
        hxferus = xx8(8)
        transp  = xx8(1)
        itest   = 10

        go to 91
     endif

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 12

     if (xx8(2) > sfcwater_mass(1)) then

        skiptest(12) = .false.

     ! If this test resulted in condensation onto the ground sfc with transpiration
     ! and sfcwater evaporation, go to test 14

     elseif (xx8(1) >= -eps_wxfer .and. xx8(2) >= -eps_wxfer .and. xx8(5) < -eps_wxfer) then

        skiptest(14) = .false.

     ! If this test resulted in condensation onto the ground and surface water
     ! with transpiration, go to test 18

     elseif (xx8(1) >= -eps_wxfer) then

        skiptest(18) = .false.

     ! If this test resulted in condensation onto vegetation and evaporation from
     ! surface water and ground, REPEAT test with no transpiration

     elseif ( xx8(2) > eps_wxfer .and. xx8(5) > eps_wxfer) then

        a4 = 0._r8
        go to 10

     ! If this test resulted in condensation onto vegetation and ground surface
     ! and evaporation from surface water, go to test 14 with NO transpiration

     elseif ( xx8(2) > eps_wxfer) then

        a4 = 0._r8
        skiptest(14) = .false.

     ! If this test resulted in condensation onto vegetation, surface water, and
     ! and ground surface then go to test 20

     else

        skiptest(20) = .false.

     endif

     11 continue

! Test 11: If vegetation is wet, and given complete evaporation of any sfcwater,
!          solve 4x4 system to test for partial evaporation from vegetation
!          and zero or positive evaporation from soil

     ! If we are here and soil is dry, do test 15 instead
     if ((.not. skiptest(11)) .and. drygnd) then
        skiptest(11)  = .true.
        skiptest(15) = .false.
     endif

     ! If we are here and vegetation is not wet, do test 12 instead
     if ((.not. skiptest(11)) .and. (.not. iwetveg)) then
        skiptest(11) = .true.
        skiptest(12) = .false.
     endif

     if (skiptest(11)) go to 12

     didtest(11) = .true.

     aa7(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa7(1,2) =         (a3 + a4) * h4
     aa7(1,3) =         (a3 + a4) * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - (a3 + a4) * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass(1))          ! WVC row

     aa7(2,1) =         a8 * h4
     aa7(2,2) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a8 * h3
     aa7(2,5) =       - a8 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a8 * h3
     yy7(2)   =         a8 * (y3 - (h3 * alvskn + h4) * sfcwater_mass(1)) ! WGC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4                                           ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))        ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * (y9 + h8 * sfcwater_mass(1))                 ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10                                         ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))      ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land11',7,aa7,yy7,glatw,glonw)

     evap = xx7(1) * a3 / max(a3 + a4, 1.e-30_r8)  ! evaporation of veg_water

     if ( xx7(1) >= -eps_wxfer .and. evap <= veg_water .and. &
          xx7(2) >= -eps_wxfer ) then ! test was successful

        wxfervc = evap
        wxfersc = sfcwater_mass(1)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxfergc = xx7(2)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        transp  = xx7(1) - evap
        itest   = 11

        go to 91
     endif

     ! If this test resulted in more than complete evaporation of veg water,
     ! go to test 12

     if (evap > veg_water) then

        skiptest(12) = .false.

     ! If this test resulted in condensation onto the ground sfc with transpiration,
     ! go to test 15

     elseif (xx7(1) >= -eps_wxfer) then

        skiptest(15) = .false.

     ! If this test resulted in condensation onto vegetation and evaporation
     ! from the ground, REPEAT test with no transpiration

     elseif (xx7(2) > eps_wxfer) then

        a3 = 0._r8
        a4 = 0._r8
        go to 11

     ! If this test resulted in condensation onto vegetation and ground,
     ! go to test 21

     else

        skiptest(21) = .false.

     endif

     12 continue

! Test 12: Given complete evaporation of any sfcwater, and complete evaporation
!          of any veg_water, solve 6x6 system to test for zero or positive
!          transpiration and zero or positive evaporation from soil

     ! If we are here and soil is dry, do test 16 instead
     if ((.not. skiptest(12)) .and. drygnd) then
        skiptest(12)  = .true.
        skiptest(16) = .false.
     endif

     if (skiptest(12)) go to 13

     didtest(12) = .true.

     aa7(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa7(1,2) =         a4 * h4
     aa7(1,3) =         a4 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a4 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water - h4 * sfcwater_mass(1))  ! WVC row

     aa7(2,1) =         a8 * h4
     aa7(2,2) = 1._r8 + a8 * (h3 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a8 * h3
     aa7(2,5) =       - a8 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a8 * h3
     yy7(2)   =         a8 * (y3 - h4 * veg_water - (h3 * alvskn + h4) * sfcwater_mass(1))  ! WGC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * (y4 - h5 * alvveg * veg_water)                                 ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))                          ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * (y9 + h8 * (veg_water + sfcwater_mass(1)))                     ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10                                                           ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))                        ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land12',7,aa7,yy7,glatw,glonw)

     if ( xx7(1) >= -eps_wxfer .and. &
          xx7(2) >= -eps_wxfer ) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass(1)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxfergc = xx7(2)
        transp  = xx7(1)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        itest   = 12

        go to 91
     endif

     ! If this test resulted in condensation onto the surface with transpiration,
     ! go to test 16

     if (xx7(1) >= -eps_wxfer) then

        skiptest(16) = .false.

     ! If this test resulted in condensation onto vegetation and evaporation from
     ! the ground, REPEAT test with no transpiration

     elseif (xx7(2) > eps_wxfer) then

        a4 = 0._r8
        go to 12

     ! If this test resulted in condensation onto vegetation and the ground,
     ! go to test 22

     else

        skiptest(22) = .false.

     endif

     13 continue

! Test 13: If some sfcwater is present and vegetation is wet, and given zero
!          vapor flux with soil, solve 6x6 system to test for partial
!          evaporation from vegetation and partial evaporation of sfcwater

     if (.not. skiptest(13)) then

        ! If vegetation is not wet, go to test 14
        if (.not. iwetveg) then
           skiptest(13) = .true.
           skiptest(14) = .false.

        ! If surface is not wet, go to test 15
        elseif (.not. iwetsfc) then
           skiptest(13) = .true.
           skiptest(15) = .false.
        endif

     endif

     if (skiptest(13)) go to 14

     didtest(13) = .true.

     aa7(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa7(1,2) =         (a3 + a4) * h4
     aa7(1,3) =         (a3 + a4) * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - (a3 + a4) * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         (a3 + a4) * y1  ! WVC row

     aa7(2,1) =         a7 * h4
     aa7(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a7 * h2
     aa7(2,5) =       - a7 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a7 * h2
     yy7(2)   =         a7 * y2         ! WSC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * y4         ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5         ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * y9         ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10       ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11       ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land13',7,aa7,yy7,glatw,glonw)

     evap = xx7(1) * a3 / max(a3 + a4, 1.e-30_r8) ! evaporation of veg_water

     if ( xx7(1) >= -eps_wxfer .and. evap   <= veg_water  .and. &
          xx7(2) >= -eps_wxfer .and. xx7(2) <= sfcwater_mass(1) ) then ! test was successful

        wxfervc = evap
        wxfersc = xx7(2)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = xx7(1) - evap
        itest   = 13

        go to 91
     endif

     ! If this test resulted in more than complete evaporation of veg_water
     ! and more than complete evaporation of sfc water, go to test 16

     if (evap > veg_water .and. xx7(2) > sfcwater_mass(1)) then

        skiptest(16) = .false.

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 15

     elseif (xx7(2) > sfcwater_mass(1)) then

        skiptest(15) = .false.

     ! If this test resulted in more than complete evaporation of veg_water,
     ! go to test 14

     elseif (evap > veg_water) then

        skiptest(14) = .false.

     ! If this test resulted in condensation onto the sfc water with transpiration,
     ! go to test 17

     elseif (xx7(1) >= -eps_wxfer) then

        skiptest(17) = .false.

     ! If this test resultd in condensation onto vegetation with evaportation
     ! from the surface water, REPEAT test with no veg water transfer

     elseif (xx7(2) > eps_wxfer) then

        a3 = 0._r8
        a4 = 0._r8
        go to 13

     ! If this test resulted in condensation onto vegetation and surface water,
     ! skip to test 19

     else

        skiptest(19) = .false.

     endif

     14 continue

! Test 14: If some sfcwater is present, given complete evaporation of any
!          veg_water, and given zero vapor flux with soil, solve 6x6 system
!          to test for zero or positive transpiration and partial evaporation of
!          any surface water

     ! If we are here and surface is not wet, do test 16 instead
     if ((.not. skiptest(14)) .and. (.not. iwetsfc)) then
        skiptest(14) = .true.
        skiptest(16) = .false.
     endif

     if (skiptest(14)) go to 15

     didtest(14) = .true.

     aa7(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa7(1,2) =         a4 * h4
     aa7(1,3) =         a4 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a4 * h4
     aa7(1,6) = 0._r8
     aa7(1,7) = 0._r8
     yy7(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water) ! WVC row

     aa7(2,1) =         a7 * h4
     aa7(2,2) = 1._r8 + a7 * (h2 * alvskn + h4)
     aa7(2,3) = 0._r8
     aa7(2,4) =         a7 * h2
     aa7(2,5) =       - a7 * h4
     aa7(2,6) = 0._r8
     aa7(2,7) =       - a7 * h2
     yy7(2)   =         a7 * (y2 - h4 * veg_water)                 ! WSC row

     aa7(3,1) =         a2 * h5 * alvveg
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) = 0._r8
     aa7(3,6) =       - a2 * h7
     aa7(3,7) = 0._r8
     yy7(3)   =         a2 * (y4 - h5 * alvveg * veg_water)        ! HVC row

     aa7(4,1) = 0._r8
     aa7(4,2) =         a6 * h6 * alvskn
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * (h6 + h7)
     aa7(4,5) = 0._r8
     aa7(4,6) =       - a6 * h7
     aa7(4,7) =       - a6 * h6
     yy7(4)   =         a6 * y5                                    ! HSC row

     aa7(5,1) =       - a9 * h8
     aa7(5,2) =       - a9 * h8
     aa7(5,3) = 0._r8
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a9 * (h8 + h10)
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a9 * (y9 + h8 * veg_water)                 ! WCA row

     aa7(6,1) = 0._r8
     aa7(6,2) = 0._r8
     aa7(6,3) =       - a10 * h7
     aa7(6,4) =       - a10 * h7
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a10 * (h7 + h11)
     aa7(6,7) = 0._r8
     yy7(6)   =         a10 * y10                                  ! HCA row

     aa7(7,1) = 0._r8
     aa7(7,2) =       - a11 * h6 * alvskn
     aa7(7,3) = 0._r8
     aa7(7,4) =       - a11 * h6
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a11 * (h6 + h9)
     yy7(7)   =         a11 * y11                                  ! HUS row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'land14',7,aa7,yy7,glatw,glonw)

     if ( xx7(1) >= -eps_wxfer .and. &
          xx7(2) >= -eps_wxfer .and. xx7(2) <= sfcwater_mass(1) ) then ! test was successful

        wxfervc = veg_water
        wxfersc = xx7(2)
        hxfervc = xx7(3)
        hxfersc = xx7(4)
        wxferca = xx7(5)
        hxferca = xx7(6)
        hxferus = xx7(7)
        wxfergc = 0.
        transp  = xx7(1)
        itest = 14

        go to 91
     endif

     ! If this test resulted in more than complete evaporation of sfc water,
     ! go to test 16

     if (xx7(2) > sfcwater_mass(1)) then

        skiptest(16) = .false.

     ! If this test resulted in condensation onto the sfc water with transpiration,
     ! go to test 18

     elseif (xx7(1) >= -eps_wxfer) then

        skiptest(18) = .false.

     ! If this test resulted in condensation onto vegetation with evaporation
     ! from the surfacw water, REPEAT test with no veg water transfer

     elseif (xx7(2) > eps_wxfer) then

        a4 = 0._r8
        go to 14

     ! If this test resulted in condensation onto vegetation and surface water,
     ! skip to test 19

     else

        skiptest(20) = .false.

     endif

     15 continue

! Test 15: If vegetation is wet, and given zero vapor flux with soil and
!          complete evaporation of any sfcwater, solve 5x5 system to test for
!          partial evaporation from vegetation

     ! If we are here and vegetation is not wet, do test 16 instead
     if ((.not. skiptest(15)) .and. (.not. iwetveg)) then
        skiptest(15) = .true.
        skiptest(16) = .false.
     endif

     if (skiptest(15)) go to 16

     didtest(15) = .true.

     aa6(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa6(1,2) =         (a3 + a4) * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - (a3 + a4) * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass(1))   ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * y4                                    ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1)) ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * (y9 + h8 * sfcwater_mass(1))          ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10                                  ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * y11                                  ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land15',6,aa6,yy6,glatw,glonw)

     evap = xx6(1) * a3 / max(a3 + a4, 1.e-30_r8) ! evaporation of veg_water

     if (xx6(1) > -eps_wxfer .and. evap <= veg_water) then ! test was successful

        wxfervc = evap
        wxfersc = sfcwater_mass(1)
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = xx6(1) - evap
        itest = 15
        go to 91
     endif

     ! If this test resulted in more than complete evaporation of veg_water,
     ! go to test 16

     if (evap > veg_water) then

        skiptest(16) = .false.

     ! If this test resulted in condensation onto vegetation,
     ! go to test 21

     else

        skiptest(21) = .false.

     endif

    16 continue

! Test 16: Given complete evaporation of any veg_water, complete evaporation
!          of any sfcwater, and zero vapor flux with soil, solve 3x3 system
!          to test for zero or positive transpiration

     if (skiptest(16)) go to 17

     didtest(16) = .true.

     aa6(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa6(1,2) =         a4 * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - a4 * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water - h4 * sfcwater_mass(1)) ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * (y4 - h5 * alvveg * veg_water)                                ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))                         ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * (y9 + h8 * (veg_water + sfcwater_mass(1)))                    ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10                                                          ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))                       ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land16',6,aa6,yy6,glatw,glonw)

     if (xx6(1) >= -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass(1)
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = xx6(1)
        itest   = 16

        go to 91
     endif

     ! This test failed due to condensation onto vegetation,
     ! go to test 22

     skiptest(22) = .false.

     17 continue

! Test 17: If vegetation is wet and given zero vapor flux with soil and
!          no evaporation of sfcwater, solve 5x5 system to test for
!          partial evaporation from vegetation

     ! If we are here and vegetation is not wet, do test 18 instead
     if ((.not. skiptest(17)) .and. (.not. iwetveg)) then
        skiptest(17) = .true.
        skiptest(18) = .false.
     endif

     if (skiptest(17)) go to 18

     didtest(17) = .true.

     aa6(1,1) = 1._r8 + (a3 + a4) * (h1 * alvveg + h4)
     aa6(1,2) =         (a3 + a4) * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - (a3 + a4) * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         (a3 + a4) * y1     ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * y4            ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * y5            ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * y9            ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10          ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * y11          ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land17',6,aa6,yy6,glatw,glonw)

     evap = xx6(1) * a3 / max(a3 + a4, 1.e-30_r8) ! evaporation of veg_water

     if (xx6(1) >= -eps_wxfer .and. evap <= veg_water) then ! test was successful

        wxfervc = evap
        wxfersc = 0.0
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = xx6(1) - evap
        itest   = 17
        go to 91
     endif

     ! If this test resulted in more than complete evaporation of veg_water,
     ! go to test 18

     if (evap > veg_water) then

        skiptest(18) = .false.

     ! If this test resulted in condensation onto vegetation,
     ! go to test 19

     else

        skiptest(19) = .false.

     endif

    18 continue

! Test 18: Given complete evaporation of any veg_water, no evaporation
!          of any sfcwater, and zero vapor flux with soil, solve 3x3 system
!          to test for zero or positive transpiration

     if (skiptest(18)) go to 19

     didtest(18) = .true.

     aa6(1,1) = 1._r8 + a4 * (h1 * alvveg + h4)
     aa6(1,2) =         a4 * h1
     aa6(1,3) = 0._r8
     aa6(1,4) =       - a4 * h4
     aa6(1,5) = 0._r8
     aa6(1,6) = 0._r8
     yy6(1)   =         a4 * (y1 - (h1 * alvveg + h4) * veg_water)  ! WVC row

     aa6(2,1) =         a2 * h5 * alvveg
     aa6(2,2) = 1._r8 + a2 * (h5 + h7)
     aa6(2,3) =         a2 * h7
     aa6(2,4) = 0._r8
     aa6(2,5) =       - a2 * h7
     aa6(2,6) = 0._r8
     yy6(2)   =         a2 * (y4 - h5 * alvveg * veg_water)         ! HVC row

     aa6(3,1) = 0._r8
     aa6(3,2) =         a6 * h7
     aa6(3,3) = 1._r8 + a6 * (h6 + h7)
     aa6(3,4) = 0._r8
     aa6(3,5) =       - a6 * h7
     aa6(3,6) =       - a6 * h6
     yy6(3)   =         a6 * y5                                     ! HSC row

     aa6(4,1) =       - a9 * h8
     aa6(4,2) = 0._r8
     aa6(4,3) = 0._r8
     aa6(4,4) = 1._r8 + a9 * (h8 + h10)
     aa6(4,5) = 0._r8
     aa6(4,6) = 0._r8
     yy6(4)   =         a9 * (y9 + h8 * veg_water)                  ! WCA row

     aa6(5,1) = 0._r8
     aa6(5,2) =       - a10 * h7
     aa6(5,3) =       - a10 * h7
     aa6(5,4) = 0._r8
     aa6(5,5) = 1._r8 + a10 * (h7 + h11)
     aa6(5,6) = 0._r8
     yy6(5)   =         a10 * y10                                   ! HCA row

     aa6(6,1) = 0._r8
     aa6(6,2) = 0._r8
     aa6(6,3) =       - a11 * h6
     aa6(6,4) = 0._r8
     aa6(6,5) = 0._r8
     aa6(6,6) = 1._r8 + a11 * (h6 + h9)
     yy6(6)   =         a11 * y11                                   ! HUS row

     call matrix8_NxN(6,aa6,yy6,xx6,sing); if (sing) call sing_print(iwsfc,'land18',6,aa6,yy6,glatw,glonw)

     if (xx6(1) >= -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = 0.
        hxfervc = xx6(2)
        hxfersc = xx6(3)
        wxferca = xx6(4)
        hxferca = xx6(5)
        hxferus = xx6(6)
        wxfergc = 0.
        transp  = xx6(1)
        itest   = 18
        go to 91
     endif

     ! This test failed due to condensation onto vegetation,
     ! go to test 19

     skiptest(20) = .false.

     19 continue

! Test 19: Given NO evaporation of veg water, NO evaporation of the
!          sfc water, NO tranpiration, and NO soil vapor flux, solve
!          3x3 system to compute the heat fluxes.

     if (.not. skiptest(19)) then

        ! If vegetation is not wet, go to test 20
        if (.not. iwetveg) then
           skiptest(19) = .true.
           skiptest(20) = .false.

        ! If surface is not wet, go to test 21
        elseif (.not. iwetsfc) then
           skiptest(19)  = .true.
           skiptest(21) = .false.
        endif

     endif

     if (skiptest(19)) go to 20

     didtest(19) = .true.

     aa4(1,1) = 1._r8 + a2 * (h5 + h7)
     aa4(1,2) =         a2 * h7
     aa4(1,3) = 0._r8
     aa4(1,4) =       - a2 * h7
     yy4(1)   =         a2 * y4            ! HVC row

     aa4(2,1) =         a6 * h7
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) =       - a6 * h6
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * y5            ! HSC row

     aa4(3,1) = 0._r8
     aa4(3,2) =       - a11 * h6
     aa4(3,3) = 1._r8 + a11 * (h6 + h9)
     aa4(3,4) = 0._r8
     yy4(3)   =         a11 * y11          ! HUS row

     aa4(4,1) =       - a10 * h7
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * (h7 + h11)
     yy4(4)   =         a10 * y10          ! HCA row

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land19',4,aa4,yy4,glatw,glonw)

     wxfervc = 0.
     wxfersc = 0.
     hxfervc = xx4(1)
     hxfersc = xx4(2)
     hxferus = xx4(3)
     wxferca = a9 * y9 / (1._r8 + a9 * (h8 + h10))
     hxferca = xx4(4)
     wxfergc = 0.
     transp  = 0.
     itest   = 19
     goto 91

     20 continue

! Test 20: Given COMPLETE evaporation of veg water, NO evaporation
!          of sfc water, and NO tranpiration and NO soil vapor flux,
!          solve 3x3 system to compute the heat fluxes

     ! If we are here and surface is not wet, do test 22 instead
     if ((.not. skiptest(20)) .and. (.not. iwetsfc)) then
        skiptest(20) = .true.
        skiptest(22) = .false.
     endif

     if (skiptest(20)) go to 21

     didtest(20) = .true.

     aa4(1,1) = 1._r8 + a2 * (h5 + h7)
     aa4(1,2) =         a2 * h7
     aa4(1,3) = 0._r8
     aa4(1,4) =       - a2 * h7
     yy4(1)   =         a2 * (y4 - h5 * alvveg * veg_water)   ! HVC row

     aa4(2,1) =         a6 * h7
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) =       - a6 * h6
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * y5                               ! HSC row

     aa4(3,1) = 0._r8
     aa4(3,2) =       - a11 * h6
     aa4(3,3) = 1._r8 + a11 * (h6 + h9)
     aa4(3,4) = 0._r8
     yy4(3)   =         a11 * y11                             ! HUS row

     aa4(4,1) =       - a10 * h7
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * (h7 + h11)
     yy4(4)   =         a10 * y10                             ! HCA row

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land20',4,aa4,yy4,glatw,glonw)

     wxfervc = veg_water
     wxfersc = 0.
     hxfervc = xx4(1)
     hxfersc = xx4(2)
     hxferus = xx4(3)
     wxferca = a9 * (y9 + h8 * veg_water) / (1._r8 + a9 * (h8 + h10))
     hxferca = xx4(4)
     wxfergc = 0.
     transp  = 0.
     itest   = 20
     go to 91

     21 continue

! Test 21: Given no evaporation of veg water, COMPLETE evaporation
!          of sfc water, and NO tranpiration and NO soil vapor flux,
!          solve 3x3 system to compute the heat fluxes

     ! If we are here and vegetation is not wet, do test 22 instead
     if ((.not. skiptest(21)) .and. (.not. iwetveg)) then
        skiptest(21) = .true.
        skiptest(22) = .false.
     endif

     if (skiptest(21)) go to 22

     didtest(21) = .true.

     aa4(1,1) = 1._r8 + a2 * (h5 + h7)
     aa4(1,2) =         a2 * h7
     aa4(1,3) = 0._r8
     aa4(1,4) =       - a2 * h7
     yy4(1)   =         a2 * y4                                       ! HVC row

     aa4(2,1) =         a6 * h7
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) =       - a6 * h6
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))    ! HSC row

     aa4(3,1) = 0._r8
     aa4(3,2) =       - a11 * h6
     aa4(3,3) = 1._r8 + a11 * (h6 + h9)
     aa4(3,4) = 0._r8
     yy4(3)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))  ! HUS row

     aa4(4,1) =       - a10 * h7
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * (h7 + h11)
     yy4(4)   =         a10 * y10                                     ! HCA row

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land21',4,aa4,yy4,glatw,glonw)

     wxfervc = 0.
     wxfersc = sfcwater_mass(1)
     hxfervc = xx4(1)
     hxfersc = xx4(2)
     hxferus = xx4(3)
     wxferca = a9 * (y9 + h8 * sfcwater_mass(1)) / (1._r8 + a9 * (h8 + h10))
     hxferca = xx4(4)
     wxfergc = 0.
     transp  = 0.
     itest   = 21
     go to 91

     22 continue

! Test 22: Given complete evaporation of veg water, complete evaporation
!          of sfc water, and NO tranpiration and NO soil vapor flux,
!          solve 3x3 system to compute the heat fluxes

     if (skiptest(22)) go to 23

     didtest(22) = .true.

     aa4(1,1) = 1._r8 + a2 * (h5 + h7)
     aa4(1,2) =         a2 * h7
     aa4(1,3) = 0._r8
     aa4(1,4) =       - a2 * h7
     yy4(1)   =         a2 * (y4 - h5 * alvveg * veg_water)           ! HVC row

     aa4(2,1) =         a6 * h7
     aa4(2,2) = 1._r8 + a6 * (h6 + h7)
     aa4(2,3) =       - a6 * h6
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * (y5 - h6 * alvskn * sfcwater_mass(1))    ! HSC row

     aa4(3,1) = 0._r8
     aa4(3,2) =       - a11 * h6
     aa4(3,3) = 1._r8 + a11 * (h6 + h9)
     aa4(3,4) = 0._r8
     yy4(3)   =         a11 * (y11 + h6 * alvskn * sfcwater_mass(1))  ! HUS row

     aa4(4,1) =       - a10 * h7
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * (h7 + h11)
     yy4(4)   =         a10 * y10                                     ! HCA row

     call matrix8_4x4(aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'land22',4,aa4,yy4,glatw,glonw)

     wxfervc = veg_water
     wxfersc = sfcwater_mass(1)
     hxfervc = xx4(1)
     hxfersc = xx4(2)
     hxferus = xx4(3)
     wxferca = a9 * (y9 + h8 * (veg_water + sfcwater_mass(1))) / (1._r8 + a9 * (h8 + h10))
     hxferca = xx4(4)
     wxfergc = 0.
     transp  = 0.
     itest   = 22
     go to 91

     23 continue

     ! If we got here, none of the above tests was satisfied, an outcome that should not happen

     itest = 91

     print*, 'itest = 91 result obtained in leaf4_canopy ',itab_wsfc(iwsfc)%iwglobe,glatw,glonw
     print*, cantemp, veg_temp, skn_tempk
     print*, radsfc, radveg
     stop 'stop test 91 '

     91 continue

     ! Update components

     sfluxt = hxferca / dt_leaf
     sfluxr = wxferca / dt_leaf

     cantemp = cantemp + (hxfervc + hxfersc - hxferca) * hcapcani
     canrrv  = canrrv  + (wxfervc + transp + wxfersc + wxfergc - wxferca) * canairi

     ! The latent heat of sublimation, alvi, must be used in computing the gain/loss of
     ! veg_energy and sfcwater_epm2 due to vapor flux.  (See explanation above in similar context.)

     veg_water  = max(0.,veg_water - wxfervc)
     veg_energy = veg_energy + radveg - hxfervc - wxfervc * alvi - transp * alvl

     call qwtk(veg_energy, veg_water, hcapveg, veg_temp, fracliqv)

     sfcwater_mass(nzs) = sfcwater_mass(nzs) - wxfersc
     sfcwater_epm2(nzs) = sfcwater_epm2(nzs) + radsfc - hxfersc - wxfersc * alvi + hxferus

     specvol = .001
     if (wxfersc > 0. .and. sfcwater_mass(nzs) > 1.e-3) specvol = sfcwater_depth(nzs) / sfcwater_mass(nzs)

     sfcwater_depth(nzs) = sfcwater_depth(nzs) - wxfersc * specvol

     if (skncomp == 1) then

        ! Sensible and latent heat contributions to the ground (soil) are added to
        ! sfcwater_epm2(1), but the subsequent call to subroutine sfcwater_soil_comb
        ! combines the energies of sfcwater(1) and soil(nzg)  .

        gndwater_tot     = max(0., soil_water(nzg) - wresid_vg(nzg)) * dslz(nzg) * 1.e3

        sfcwater_epm2(1) = sfcwater_epm2(1) - min(wxfergc, gndwater_tot) * alvi &
                                       + max(0., wxfergc - gndwater_tot) * alvl

        delw             = dslzi(nzg) * min(wxfergc, gndwater_tot) * .001
        fact             = wsat_vg(nzg) - wresid_vg(nzg)

        soil_water(nzg)  = soil_water(nzg) - delw
        head      (nzg)  = head      (nzg) - delw * head_slope(nzg)
        soil_wfrac(nzg)  = soil_wfrac(nzg) - delw / fact
        soil_wfrac(nzg)  = max(0., min(1., soil_wfrac(nzg)))

        ! With skncomp = 1, any sfcwater mass present was brought into thermal equilibrium
        ! with nzg soil layer.  Perform a second thermal equilibration to correctly
        ! distribute changes in sfcwater_epm2 that were added in this subroutine.

        call sfcwater_soil_comb( iland, iwsfc, soil_water(nzg), soil_energy(nzg), &
                                 specifheat_drysoil(nzg), sfcwater_mass(1), sfcwater_epm2(1), &
                                 sfcwater_tempk, sfcwater_fracliq )

        soil_tempk  (nzg) = sfcwater_tempk
        soil_fracliq(nzg) = sfcwater_fracliq

        soil_energy(nzg-1) = soil_energy(nzg-1) - hxferus * dslzi(nzg-1)

        call qwtk(soil_energy(nzg-1), soil_water(nzg-1)*1.e3, &
                  specifheat_drysoil(nzg-1), soil_tempk(nzg-1), soil_fracliq(nzg-1))

        ! In case too much water evaporated from the skin layer, try to remove it from layer nzg-1
        ! This should not be necessary but just in case...

        if (wxfergc - gndwater_tot > eps_wxfer) then
           delw              = dslzi(nzg-1) * (wxfergc-gndwater_tot) * .001
           delw              = min(delw, soil_water(nzg-1) - wresid_vg(nzg-1))
           fact              = wsat_vg(nzg) - wresid_vg(nzg-1)

           soil_energy(nzg-1) = soil_energy(nzg-1) + (specifheat_drysoil(nzg-1) &
                              * (soil_tempk(nzg-1) - t00) - soil_energy(nzg-1)) &
                              * delw / max(soil_water(nzg-1), wresid_vg(nzg-1))

           soil_water (nzg-1) = soil_water(nzg-1) - delw
           head       (nzg-1) = head      (nzg-1) - delw * head_slope(nzg-1)
           soil_wfrac (nzg-1) = soil_wfrac(nzg-1) - delw / fact
           soil_wfrac (nzg-1) = max(0., min(1., soil_wfrac(nzg-1)))
        endif

     else

        if (nzs > 1) then
           sfcwater_epm2(nzs-1) = sfcwater_epm2(nzs-1) - hxferus
        else
           soil_energy(nzg) = soil_energy(nzg) - hxferus * dslzi(nzg)
           call qwtk(soil_energy(nzg), soil_water(nzg)*1.e3, &
                     specifheat_drysoil(nzg), soil_tempk(nzg), soil_fracliq(nzg))
        endif

     endif

  endif  ! End of vegetation section

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
        + (1. - veg_clump(leaf_class)) * log(1. - fpar) * fpcon)

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

  subroutine fast_canopy(iland, iwsfc,                                       &
                         sfcwat_nud,     sfctemp_nud,    fracliq_nud,        &
                         sfcwater_mass,  sfcwater_epm2,  sfcwater_depth,     &
                         soil_water,     soil_energy,    specifheat_drysoil, &
                         soil_tempk,     soil_fracliq                        )

  use leaf_coms,      only: dt_leaf
  use therm_lib,      only: qwtk
  use leaf4_surface,  only: sfcwater_soil_comb

  implicit none

  integer, intent(in)    :: iland      ! index of current land cell
  integer, intent(in)    :: iwsfc      ! index of current SFC grid cell

  real, intent(in) :: sfcwat_nud   ! clim net water mass input (P-ET) [kg/(m^2 s)]
  real, intent(in) :: sfctemp_nud  ! clim sfc temp [K]
  real, intent(in) :: fracliq_nud  ! clim sfc fracliq []

  real, intent(inout) :: sfcwater_mass    ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_epm2    ! surface water energy [J/m^2]
  real, intent(inout) :: sfcwater_depth   ! surface water depth [m]

  real, intent(inout) :: soil_water         ! top soil layer water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy        ! top soil layer energy [J/m^3]
  real, intent(in)    :: specifheat_drysoil ! specific heat of dry soil [J/(m^3 K)]
  real, intent(inout) :: soil_tempk         ! top soil layer temperature
  real, intent(inout) :: soil_fracliq       ! top soil layer liquid fraction

  ! Local variables

  real :: flux, ediff
  real :: sfcwater_tempk   ! surface water temperature [K]
  real :: sfcwater_fracliq ! fraction of sfc water in liquid phase

! PLAN: nudge temp/energy of sfcwater and soil(nzg); add net water (PCP - E - T)
! to sfcwater and soil.  Do runoff ONLY in locations/times (months?) when there
! is net water loss, and do runoff only when there is surface water present.
! Do not let sfcwater exceed a prescribed upper limit (of 2 m depth water equiv).

  ! Apply net change (P - ET) to sfcwater

  sfcwater_mass = max(0.,min(2000.,sfcwater_mass + sfcwat_nud * dt_leaf))

  ! In FAST_CANOPY, always thermally combine sfcwater and soil(nzg).

  call sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy, &
                          specifheat_drysoil, sfcwater_mass, sfcwater_epm2, &
                          sfcwater_tempk, sfcwater_fracliq)

  ! Compute surface energy flux based on difference between model and nudging
  ! values of sfcwater_tempk and sfcwater_fracliq

  ediff = sfctemp_nud - sfcwater_tempk + 80. * (fracliq_nud - sfcwater_fracliq)

  flux = max(-500.,min(500.,20. * ediff))

  sfcwater_epm2 = sfcwater_epm2 + flux * dt_leaf

  ! Perform a second thermal equilibration to correctly distribute changes in
  ! sfcwater_epm2 that were added in this subroutine.

  call sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy, &
                          specifheat_drysoil, sfcwater_mass, sfcwater_epm2, &
                          sfcwater_tempk, sfcwater_fracliq)

  soil_tempk   = sfcwater_tempk
  soil_fracliq = sfcwater_fracliq

end subroutine fast_canopy

! Remaining issues:

! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Build lookup tables, especially for things with exponentials?

End module leaf4_canopy
