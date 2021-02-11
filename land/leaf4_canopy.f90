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
Module leaf4_canopy

Contains

  subroutine canopy(iland, iwsfc, nlsw1, icomb, ktrans, transp,                  &
                    leaf_class,    can_depth,       rhos,            vels,       &
                    ustar,         sxfer_t,         sxfer_r,         pcpg,       &
                    qpcpg,         dpcpg,           rshort,          cantemp,    &
                    canrrv,        glatw,           glonw,                       &
                    snowfac,       vf,              stom_resist,     veg_height, &
                    veg_rough,     veg_tai,         veg_lai,         hcapveg,    &
                    rshort_v,      rshort_g,        rlong_v,         rlong_s,    &
                    rlong_g,       veg_water,       vegwater_energy, veg_temp,   &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,              &
                    energy_per_m2, sfcwater_tempk,  sfcwater_fracliq,            &
                    soil_water,    soil_energy,                                  &
                    wsat_vg,       ksat_vg,         specifheat_drysoil,          &
                    head,          soil_tempk,      soil_fracliq                 )

  use leaf_coms,     only: soil_rough, dt_leaf, kroot, rcmin, snowmin_expl, &
                           wcap_min, wcap_vmin
  use mem_sfcg,      only: sfcg
  use mem_land,      only: nzg, dslz, dslzi, slzt
  use leaf4_surface, only: sfcwater_soil_comb, grndvap_ab
  use consts_coms,   only: cp, vonk, alvl, cliq, cice, alli, rvap, t00, r8
  use misc_coms,     only: io6
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
  real, intent(in)    :: sxfer_t          ! canopy-to-atm sens heat xfer this step [kg_dryair K/m^2]
  real, intent(in)    :: sxfer_r          ! canopy-to-atm vapor xfer this step [kg_vap/m^2]
  real, intent(in)    :: pcpg             ! added precip mass this leaf timestep [kg/m^2]
  real, intent(in)    :: qpcpg            ! added precip energy this leaf timestep [J/m^2]
  real, intent(in)    :: dpcpg            ! added precip depth this leaf timestep [m]
  real, intent(in)    :: rshort           ! downward sfc s/w rad flux [W/m^2]
  real, intent(inout) :: cantemp          ! canopy air temp [K]
  real, intent(inout) :: canrrv           ! canopy air vapor mixing ratio [kg_vap/kg_dryair]
  real, intent(in)    :: glatw            ! Latitude of land cell 'center' [deg]
  real, intent(in)    :: glonw            ! Longitude of land cell 'center' [deg]
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
  real, intent(inout) :: vegwater_energy  ! veg sfc water energy [J/m^2]
  real, intent(inout) :: veg_temp         ! veg temp [K]
  real, intent(inout) :: sfcwater_mass    ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_energy  ! surface water energy [J/kg]
  real, intent(inout) :: sfcwater_depth   ! surface water depth [m]
  real, intent(inout) :: energy_per_m2    ! surface water energy [J/m^2]
  real, intent(inout) :: sfcwater_tempk   ! surface water temperature [K]
  real, intent(inout) :: sfcwater_fracliq ! fraction of sfc water in liquid phase
  real, intent(inout) :: soil_water        (nzg) ! soil water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy       (nzg) ! soil energy [J/m^3]
  real, intent(in)    :: wsat_vg           (nzg) ! saturation water content (porosity) []
  real, intent(in)    :: ksat_vg           (nzg) ! saturation hydraulic conductivity [m/s]
  real, intent(in)    :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]
  real, intent(in)    :: head              (nzg) ! hydraulic head [m] (relative to local topo datum)
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

  real :: hxfersc     ! sfc-to-can_air heat xfer this step [J/m^2]
  real :: wxfersc     ! sfc-to-can_air vap xfer this step [kg_vap/m^2]
  real :: wxfergc     ! gnd-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxfervc     ! veg-to-can_air heat xfer this step [J/m^2]
  real :: wxfervc     ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
  real :: wshed       ! water shed from veg this LEAF timestep [kg/m^2]
  real :: qshed       ! water energy shed from veg this LEAF timestep [J/m^2]
  real :: dshed       ! depth of water shed from veg [m]
  real :: veg_energy  ! combined veg and veg_water energy [J/m^2]
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
  real :: hcapsfc     ! Surface heat capacity [J/(m^2 K)]

  real :: canairi     ! Inverse of canair
  real :: hcapcani    ! Inverse of hcapcan
  real :: hcapsfci    ! Inverse of hcapsfc
  real :: hcapvegi    ! Inverse of hcapveg

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
  real :: gnd_rhov1   ! gnd_rhov at 1 K warmer [kg_vap/m^3]
  real :: gnd_rhovp   ! gnd_rhov derivative with respect to temperature

  real :: eps_wxfer   ! very small water transfer to accomodate truncation error
  real :: sfcwater_fracarea
  real :: evap        ! amount evaporated from veg sfc [kg/m^2]

  real :: specvol, vw_max
  real :: tf, wadd, wshed2, qshed2

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8
  real(r8) :: h1, h2, h3, h4, h5, h6, h7
  real(r8) :: y1, y2, y3, y4, y5

  real(r8) :: aa2(2,2), xx2(2), yy2(2)  ! 2x2 matrix equation terms
  real(r8) :: aa3(3,3), xx3(3), yy3(3)  ! 3x3 matrix equation terms
  real(r8) :: aa4(4,4), xx4(4), yy4(4)  ! 4x4 matrix equation terms
  real(r8) :: aa5(5,5), xx5(5), yy5(5)  ! 5x5 matrix equation terms

  logical :: sing, skiptest(15)

  integer :: itest  ! test number

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

     ! Evaluate vegetation Celsius temperature

     tvegc      = veg_temp - t00
     veg_energy = hcapveg * tvegc + vegwater_energy

     ! Remove any existing veg water over threshold

     if (veg_water > vw_max) then
        wshed = veg_water - vw_max
        qshed = wshed * vegwater_energy / veg_water
        dshed = wshed * .001

        veg_water       = vw_max
        vegwater_energy = vegwater_energy - qshed
        veg_energy      = veg_energy      - qshed
     endif

     ! Add any precipitation intercepted mass and energy to vegetation surface

     if (pcpg > 1.e-30) then
        veg_water       = veg_water       + pcpg  * vf
        vegwater_energy = vegwater_energy + qpcpg * vf

        if (veg_water > vw_max) then
           wshed2 = veg_water - vw_max
           qshed2 = wshed2 * vegwater_energy / veg_water

           veg_water       = vw_max
           vegwater_energy = vegwater_energy - qshed2
           veg_energy      = veg_energy - qshed2

           wshed = wshed + wshed2
           qshed = qshed + qshed2
           dshed = dshed + wshed2 * dpcpg / pcpg
        endif

        call qwtk(veg_energy, veg_water, hcapveg, veg_temp, fracliqv)
        tvegc = veg_temp - t00
     endif

     ! Check for presence of veg_water and set iwetveg flag

     iwetveg = (veg_water > wcap_vmin)
     if (.not. iwetveg) then
        if (veg_water > 1.e-20) canrrv = canrrv + veg_water * canairi
        veg_water       = 0.
        vegwater_energy = 0.
        veg_energy      = tvegc * hcapveg
     endif

     ! Vegetation saturation vapor density and gradient

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
     sfcwater_energy = sfcwater_energy + qpcpg * (1. - vf) + qshed
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

     ! Diagnose combined sfc heat capacity for sfcwater and top soil layer

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
     ! icomb flag to 0 and diagnose sfc heat capacity for sfcwater alone.
     ! First, diagnose sfcwater temperature and liquid fraction.  Use qwtk
     ! instead of qtk because energy_per_m2 is used instead of sfcwater_energy.
     ! ("dryhcap" = 100 is very small value)

     icomb = 0

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

  ! Sfcwater saturation vapor density and derivative

  sfc_tempk = sfcwater_tempk

  sfc_tempk1 = sfc_tempk + 1.

  sfc_rhovs  = rhovsil(sfc_tempk -273.15)
  sfc_rhovs1 = rhovsil(sfc_tempk1-273.15)
  sfc_rhovsp = sfc_rhovs1 - sfc_rhovs

  ! Reduced saturation vapor density and gradient of soil during evaporation

  call grndvap_ab(iland, sfc_tempk, soil_water(nzg), head(nzg), can_rhov, sfc_rhovs, gnd_rhov)

  gnd_rhovp = sfc_rhovsp * gnd_rhov / sfc_rhovs

  ! The remainder of subroutine canopy sets up and solves a linear system of
  ! equations using the trapezoidal-implicit method.  The solution of the system 
  ! consists of turbulent heat and water vapor fluxes between canopy air, the
  ! surface (sfcwater or soil), and vegetation (if present), and the consequent
  ! changes to water and energy content of canopy air, the surface, and vegetation.
  ! The solution obtains values for the following quantities:

  ! 1. Surface to canopy air heat transfer
  ! 2. Surface to canopy air water vapor transfer
  ! 3. Vegetation to canopy air heat transfer
  ! 4. Vegetation to canopy air water vapor transfer from transpiration
  ! 5. Vegetation to canopy air water vapor transfer from veg surface water
  ! 6. Canopy air temperature update
  ! 7. Canopy specific humidity update
  ! 8. Vegetation temperature update
  ! 9. Vegetation surface water update

  ! External forcing terms for this system are:

  ! 1. Canopy air to free atmosphere fluxes of sensible heat and vapor
  ! 2. Radiative fluxes to vegetation and to the surface

  ! (Precipitation/shedding fluxes to vegetation and to the surface were already
  ! applied above.  Updates to surface water and energy are applied later.)

  ! It may be desirable at some point to incorporate the sensible heat and vapor
  ! fluxes between canopy air and the free atmosphere within the implicit system.
  ! This would require that subroutine stars not compute these fluxes itself, but
  ! that it instead compute surface exchange coefficients for sensible heat and
  ! vapor (as it currently does for momentum).  It would also require that most 
  ! of the code below in this subroutine be transferred to a separate subroutine
  ! so that it can be called for every flux cell.  In this structural change,
  ! turbulent transfers between canopy air, vegetation, and the surface would be
  ! attached to the flux cell data structure during computation and later summed
  ! over each land cell as canopy air to free atmosphere fluxes already are.
  ! This change to the code would permit the canopy storage capacity to be
  ! arbitrarily small, but would then place the stability burden on the lowest
  ! predicted atmospheric layer.

  ! Quantities that apply both WITH and WITHOUT VEGETATION
  ! [radsfc is net shortwave + longwave input to surface.  In LEAF4, all shortwave
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

  a5 = dt_leaf * rdi                  ! sfc vap xfer coef  (all area)
  a6 = cp * rhos * a5                 ! sfc heat xfer coef (all area)
  a7 = a5 * sfcwater_fracarea         ! sfc vap xfer coef  (wet areas only)
  a8 = a5 * (1. - sfcwater_fracarea)  ! sfc vap xfer coef  (wet areas only)

  h2 = fcn * sfc_rhovsp * hcapsfci
  h3 = fcn * gnd_rhovp  * hcapsfci
  h4 = fcn * rhos * canairi
  h6 = fcn * hcapsfci
  h7 = fcn * hcapcani
                                                                    !!       CAN VAP            CAN VAP           RHOVS & TEMP
  y2 = sfc_rhovs - can_rhov + h2 * radsfc + h4 *      sxfer_r  !!  -h4 * veg_water - h4 * sfcwater_mass - h2 * alvl * sfcwater_mass
  y5 = sfc_tempk - cantemp  + h6 * radsfc + h7 * cp * sxfer_t  !!                                       - h6 * alvl * sfcwater_mass
  y3 = gnd_rhov  - can_rhov + h3 * radsfc + h4 *      sxfer_r  !!  -h4 * veg_water - h4 * sfcwater_mass - h3 * alvl * sfcwater_mass

  if (iveg == 0) then

     ! Case with NO VEGETATION

     ! Set all skiptest values to .false.

     skiptest(:) = .false.

     ! Fill arrays for matrix solution.

! Test 101:  Solve 2x2 system to test for either condensation onto surface or
!           partial evaporation from sfcwater that completely covers the surface.
!           (When sfcwater completely covers the surface, it is too abundant to
!           completely evaporate in a single time step.)
 
     aa2(1,1) = 1._r8 + a5 * (h2 * alvl + h4)
     aa2(1,2) =         a5 * h2
     yy2(1)   =         a5 * y2         ! WSC row

     aa2(2,1) =         a6 * h6 * alvl
     aa2(2,2) = 1._r8 + a6 * (h6 + h7)
     yy2(2)   =         a6 * y5         ! HSC row

     call matrix8_2x2(AA2,YY2,XX2,sing); if (sing) call sing_print(iland,101,2,aa2,yy2,glatw,glonw)

     if (xx2(1) < eps_wxfer .or.  &
        (xx2(1) <= sfcwater_mass .and. sfcwater_fracarea > 0.999)) then ! test was successful

        wxfersc = xx2(1)
        hxfersc = xx2(2)
        wxfergc = 0.
        itest = 101
        go to 120
     endif

! All remaining situations do NOT have complete coverage of surface by sfcwater
! and do NOT have condensation onto surface.

     if (sfcwater_fracarea > 0.999) then
        write(io6,'(a,i10,2f9.2)') 'NOVEG CASE sfcwater_fracarea should be less than 1', &
            iwsfc,glatw,glonw
        stop 'sfcwater_fracarea = 1, NOVEG CASE'
     endif

! Test 102:  If some sfcwater is present, solve 3x3 system to test for
!            partial evaporation from sfcwater and evaporation from soil

     if (sfcwater_mass > 0.) then

        aa3(1,1) = 1._r8 + a7 * (h2 * alvl + h4)
        aa3(1,2) =         a7 * h2
        aa3(1,3) =         a7 * h4
        yy3(1)   =         a7 * y2    ! WSC row

        aa3(2,1) =         a6 * h6 * alvl
        aa3(2,2) = 1._r8 + a6 * (h6 + h7)
        aa3(2,3) =         a6 * h6 * alvl
        yy3(2)   =         a6 * y5    ! HSC row

        aa3(3,1) =         a8 * h4
        aa3(3,2) =         a8 * h3
        aa3(3,3) = 1._r8 + a8 * (h3 * alvl + h4)
        yy3(3)   =         a8 * y3    ! WGC row

        call matrix8_3x3(AA3,YY3,XX3,sing); if (sing) call sing_print(iland,102,3,aa3,yy3,glatw,glonw)

        if (xx3(1) > -eps_wxfer .and. xx3(1) <= sfcwater_mass .and. &
            xx3(3) > -eps_wxfer) then ! test was successful
 
           wxfersc = xx3(1)
           hxfersc = xx3(2)
           wxfergc = xx3(3)
           itest = 102
           go to 120
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 103, 104, and 105

        if (xx3(1) < -eps_wxfer) then
           skiptest(3) = .true.
           skiptest(4) = .true.
           skiptest(5) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 103

        if (xx3(1) <= sfcwater_mass) then
           skiptest(3) = .true.
        endif
 
     endif

     103 continue
     if (skiptest(3)) go to 104

! Test 103: Given complete evaporation of any sfcwater, solve 2x2 system to
!           test for positive evaporation from soil

     aa2(1,1) = 1._r8 + a8 * (h3 * alvl + h4)
     aa2(1,2) =         a8 * h3
     yy2(1)   =         a8 * (y3 - (h3 * alvl + h4) * sfcwater_mass)  ! WGC row

     aa2(2,1) =         a6 * h6 * alvl
     aa2(2,2) = 1._r8 + a6 * (h6 + h7)
     yy2(2)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)         ! HSC row

     call matrix8_2x2(AA2,YY2,XX2,sing); if (sing) call sing_print(iland,103,2,aa2,yy2,glatw,glonw)

     if (xx2(1) > -eps_wxfer) then ! test was successful

        wxfersc = sfcwater_mass
        hxfersc = xx2(2)
        wxfergc = xx2(1)
        itest = 103
        go to 120
     endif

     104 continue
     if (skiptest(4)) go to 105

! Test 104: If some sfcwater is present, and given zero vapor flux with soil,
!           solve 2x2 system to test for partial evaporation of sfcwater

     if (sfcwater_mass > 0.) then

        aa2(1,1) = 1._r8 + a7 * (h2 * alvl + h4)
        aa2(1,2) =         a7 * h2
        yy2(1)   =         a7 * y2        ! WSC row

        aa2(2,1) =         a6 * h6 * alvl
        aa2(2,2) = 1._r8 + a6 * (h6 + h7)
        yy2(2)   =         a6 * y5        ! HSC row

        call matrix8_2x2(AA2,YY2,XX2,sing); if (sing) call sing_print(iland,104,2,aa2,yy2,glatw,glonw)

        if (xx2(1) > -eps_wxfer .and. xx2(1) < sfcwater_mass) then ! test was successful

           wxfersc = xx2(1)
           hxfersc = xx2(2)
           wxfergc = 0.
           itest = 104
           go to 120
        endif

     endif

     105 continue
     if (skiptest(5)) go to 106

! Test 105: Given zero vapor flux with soil and given complete evaporation of
!           any sfcwater, solve 1x1 system (a reduction of the 2x2 system in
!           Test 104) to get surface heat flux

     aa2(2,2) = 1._r8 + a6 * (h6 + h7)
     yy2(2)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)  ! HSC row

     wxfersc = sfcwater_mass
     wxfergc = 0.
     hxfersc = yy2(2) / aa2(2,2)

     itest = 105
     go to 120

     106 continue

     ! If we got here, none of the above tests was satisfied, an outcome that should not happen

     itest = 92

     print*, 'itest = 92 result obtained in leaf4_canopy ',iwsfc,glatw,glonw
     stop 'stop test 92 '

     120 continue

     ! Update components

     cantemp = cantemp + (hxfersc - cp * sxfer_t) * hcapcani

     energy_per_m2 = energy_per_m2 + radsfc - hxfersc - (wxfersc + wxfergc) * alvl

     canrrv = canrrv + (wxfersc + wxfergc - sxfer_r) * canairi



     soil_water(nzg) = soil_water(nzg) - dslzi(nzg) * wxfergc * .001

     sfcwater_mass = sfcwater_mass - wxfersc

     specvol = .001
     if (wxfersc > 0. .and. sfcwater_mass > 1.e-3) specvol = sfcwater_depth / sfcwater_mass

     sfcwater_depth = sfcwater_depth - wxfersc * specvol

     transp = 0.

  else

     ! Case WITH VEGETATION

     hcapvegi = 1. / hcapveg

     radveg = dt_leaf * (rshort_v + rlong_v)

     ! veg_water fractional coverage

     sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66666666)

     a1 = dt_leaf * 2. * stai / rb                  ! veg vap xfer coef  (all area)
     a2 = cp * rhos * a1                            ! veg heat xfer coef (all area)
     a3 = a1 * sigmaw                               ! veg vap xfer coef  (veg_water evap)
     a4 = dt_leaf * slai * (1.-sigmaw) / (rb + rc)  ! veg vap xfer coef  (transp)

     ! Auxiliary quantities

     h1 = fcn * veg_rhovsp * hcapvegi
     h5 = fcn * hcapvegi
                                                                       !!       CAN VAP            CAN VAP           RHOVS & TEMP
     y1 = veg_rhovs - can_rhov + h1 * radveg + h4 *      sxfer_r  !!  -h4 * veg_water - h4 * sfcwater_mass - h1 * alvl * veg_water
     y4 = veg_temp  - cantemp  + h5 * radveg + h7 * cp * sxfer_t  !!                                       - h5 * alvl * veg_water

     ! Fill arrays for matrix solution (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, vegetation, and the surface.
     ! For vegetation, there are 3 distinct situations to consider:

     ! 1. Condensation onto vegetation (Vc)
     ! 2. Evaporation from vegetation with veg water not depleted (Vew)
     ! 3. Evaporation from vegetation with veg water depleted or zero (Ve)

     ! For the surface, there are 4 distinct situations to consider:

     ! 1. Condensation onto surface (Gc) or full sfcwater coverage with partial evaporation (Sew)
     ! 2. Evaporation from both sfcwater and soil with sfcwater not depleting (GeSew)
     ! 3. Evaporation from both sfcwater and soil with sfcwater depleting or zero (GeSe)
     ! 4. Evaporation from sfcwater only with sfcwater not depleting (Sew)
     ! 5. Evaporation from sfcwater only with sfcwater depleting or zero (Se)

     ! (Zero vapor flux from exposed soil occurs when soil is too warm for condensation
     ! and too dry for evaporation.)

     ! There are 15 possible pairs containing one item from each list,
     ! only one of which is the correct description for a particular land
     ! grid cell and timestep.

     ! 4x4  Test  1  Vc  (Sc or Sew)               (or sfcwater_fracarea = 1)
     ! 4x4  Test  2  Vew (Sc or Sew)  (iwetveg=t)  (or sfcwater_fracarea = 1)
     ! 4x4  Test  3  Ve  (Sc or Sew)               (or sfcwater_fracarea = 1)

     ! 5x5  Test  4  Vc  GeSew                     (sfcwater_mass > 0.)
     ! 5x5  Test  5  Vew GeSew        (iwetveg=t)  (sfcwater_mass > 0.)
     ! 5x5  Test  6  Ve  GeSew                     (sfcwater_mass > 0.)

     ! 4x4  Test  7  Vc  GeSe                 
     ! 4x4  Test  8  Vew GeSe         (iwetveg=t) 
     ! 4x4  Test  9  Ve  GeSe                  

     ! 4x4  Test 10  Vc  Sew                       (sfcwater_mass > 0.)
     ! 4x4  Test 11  Vew Sew          (iwetveg=t)  (sfcwater_mass > 0.)
     ! 4x4  Test 12  Ve  Sew                       (sfcwater_mass > 0.)

     ! 3x3  Test 13  Vc  Se                    
     ! 3x3  Test 14  Vew Se           (iwetveg=t)  
     ! 3x3  Test 15  Ve  Se                    

     ! We proceed to solve each equation in the above order until finding one
     ! that produces fluxes of the correct sign for the particular equation.

     ! Set all skiptest values to .false.

     skiptest(:) = .false.

! Test 1: Solve 4x4 system to test for condensation onto vegetation and either
!         condensation onto surface or partial evaporation from sfcwater that
!         completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)
 
     aa4(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
     aa4(1,2) =         a1 * h4
     aa4(1,3) =         a1 * h1
     aa4(1,4) = 0._r8
     yy4(1)   =         a1 * y1  ! WVC row

     aa4(2,1) =         a5 * h4
     aa4(2,2) = 1._r8 + a5 * (h2 * alvl + h4)
     aa4(2,3) = 0._r8
     aa4(2,4) =         a5 * h2
     yy4(2)   =         a5 * y2  ! WSC row

     aa4(3,1) =         a2 * h5 * alvl
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a2 * (h5 + h7)
     aa4(3,4) =         a2 * h7
     yy4(3)   =         a2 * y4  ! HVC row

     aa4(4,1) = 0._r8
     aa4(4,2) =         a6 * h6 * alvl
     aa4(4,3) =         a6 * h7
     aa4(4,4) = 1._r8 + a6 * (h6 + h7)
     yy4(4)   =         a6 * y5  ! HSC row

     call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,1,4,aa4,yy4,glatw,glonw)

     if (xx4(1) < eps_wxfer .and. &
        (xx4(2) < eps_wxfer .or. (xx4(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

        wxfervc = xx4(1)
        wxfersc = xx4(2)
        hxfervc = xx4(3)
        hxfersc = xx4(4)
        wxfergc = 0.
        transp  = 0.
        itest = 1
        go to 20
     endif

     ! If this test resulted in upward flux from the surface with
     ! sfcwater_fracarea less than 1, skip past tests 2 and 3

     if (xx4(2) > eps_wxfer .and. sfcwater_fracarea <= 0.999) then
        skiptest(2) = .true.
        skiptest(3) = .true.
     endif
 
     2 continue
     if (skiptest(2)) go to 3

! Test 2: If vegetation is wet, solve 4x4 system to test for partial evaporation from
!         vegetation and either condensation onto surface or partial evaporation from
!         sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     if (iwetveg) then

        aa4(1,1) = 1._r8 + (a3 + a4) * (h1 * alvl + h4)
        aa4(1,2) =         (a3 + a4) * h4
        aa4(1,3) =         (a3 + a4) * h1
        aa4(1,4) = 0._r8
        yy4(1)   =         (a3 + a4) * y1  ! WVC row

        aa4(2,1) =         a5 * h4
        aa4(2,2) = 1._r8 + a5 * (h2 * alvl + h4)
        aa4(2,3) = 0._r8
        aa4(2,4) =         a5 * h2
        yy4(2)   =         a5 * y2         ! WSC row

        aa4(3,1) =         a2 * h5 * alvl
        aa4(3,2) = 0._r8
        aa4(3,3) = 1._r8 + a2 * (h5 + h7)
        aa4(3,4) =         a2 * h7
        yy4(3)   =         a2 * y4         ! HVC row

        aa4(4,1) = 0._r8
        aa4(4,2) =         a6 * h6 * alvl
        aa4(4,3) =         a6 * h7
        aa4(4,4) = 1._r8 + a6 * (h6 + h7)
        yy4(4)   =         a6 * y5         ! HSC row

        call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,2,4,aa4,yy4,glatw,glonw)

        evap = xx4(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx4(1) > -eps_wxfer .and. evap <= veg_water .and. &
           (xx4(2) <  eps_wxfer .or. (xx4(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

           wxfervc = evap
           wxfersc = xx4(2)
           hxfervc = xx4(3)
           hxfersc = xx4(4)
           wxfergc = 0.
           transp  = xx4(1) - evap
           itest = 2
           go to 20
        endif

        ! If this test resulted in upward flux from the surface with
        ! sfcwater_fracarea less than 1, skip past test 3

        ! If this test resulted in less than complete evaporation of veg_water
        ! skip past test 3

        if (xx4(2) > eps_wxfer .and. sfcwater_fracarea <= 0.999) then
           skiptest(3) = .true.
        elseif (evap <= veg_water) then
           skiptest(3) = .true.
        endif
 
     endif

     3 continue
     if (skiptest(3)) go to 4

! Test 3: Given complete evaporation of any veg_water, solve 4x4 system to test
!         for transpiration and either condensation onto surface or partial
!         evaporation from sfcwater that completely covers the surface.
!         (When sfcwater completely covers the surface, it is too abundant to
!         completely evaporate in a single time step.)

     aa4(1,1) = 1._r8 + a4 * (h1 * alvl + h4)
     aa4(1,2) =         a4 * h4
     aa4(1,3) =         a4 * h1
     aa4(1,4) = 0._r8
     yy4(1)   =         a4 * (y1 - (h1 * alvl + h4) * veg_water)  ! WVC row

     aa4(2,1) =         a5 * h4
     aa4(2,2) = 1._r8 + a5 * (h2 * alvl + h4)
     aa4(2,3) = 0._r8
     aa4(2,4) =         a5 * h2
     yy4(2)   =         a5 * (y2 - h4 * veg_water)                ! WSC row

     aa4(3,1) =         a2 * h5 * alvl
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a2 * (h5 + h7)
     aa4(3,4) =         a2 * h7
     yy4(3)   =         a2 * (y4 - h5 * alvl * veg_water)         ! HVC row

     aa4(4,1) = 0._r8
     aa4(4,2) =         a6 * h6 * alvl
     aa4(4,3) =         a6 * h7
     aa4(4,4) = 1._r8 + a6 * (h6 + h7)
     yy4(4)   =         a6 * y5                                   ! HSC row

     call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,3,4,aa4,yy4,glatw,glonw)

     if (xx4(1) > -eps_wxfer .and. &
        (xx4(2) <  eps_wxfer .or. (xx4(2) <= sfcwater_mass .and. sfcwater_fracarea > 0.999))) then ! test was successful

        wxfervc = veg_water
        wxfersc = xx4(2)
        hxfervc = xx4(3)
        hxfersc = xx4(4)
        wxfergc = 0.
        transp  = xx4(1)
        itest = 3
        go to 20
     endif

     4 continue

! All remaining situations do NOT have complete coverage of surface by sfcwater
! and do NOT have condensation onto surface.

     if (sfcwater_fracarea > 0.999) then
        write(io6,'(a,i10,2f9.2)') 'VEG CASE sfcwater_fracarea should be less than 1', &
            iwsfc,glatw,glonw
        stop 'sfcwater_fracarea = 1, VEG CASE'
     endif

     if (skiptest(4)) go to 5

! Test 4: If some sfcwater is present, solve 5x5 system to test for
!         condensation onto vegetation, partial evaporation of sfcwater,
!         and evaporation from soil

     if (sfcwater_mass > 0.) then

        aa5(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
        aa5(1,2) =         a1 * h4
        aa5(1,3) =         a1 * h1
        aa5(1,4) = 0._r8
        aa5(1,5) =         a1 * h4
        yy5(1)   =         a1 * y1  ! WVC row

        aa5(2,1) =         a7 * h4
        aa5(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa5(2,3) = 0._r8
        aa5(2,4) =         a7 * h2
        aa5(2,5) =         a7 * h4
        yy5(2)   =         a7 * y2  ! WSC row

        aa5(3,1) =         a2 * h5 * alvl
        aa5(3,2) = 0._r8
        aa5(3,3) = 1._r8 + a2 * (h5 + h7)
        aa5(3,4) =         a2 * h7
        aa5(3,5) =         a2 * h5 * alvl
        yy5(3)   =         a2 * y4  ! HVC row

        aa5(4,1) = 0._r8
        aa5(4,2) =         a6 * h6 * alvl
        aa5(4,3) =         a6 * h7
        aa5(4,4) = 1._r8 + a6 * (h6 + h7)
        aa5(4,5) =         a6 * h6 * alvl
        yy5(4)   =         a6 * y5  ! HSC row

        aa5(5,1) =         a8 * h4
        aa5(5,2) =         a8 * h4
        aa5(5,3) = 0._r8
        aa5(5,4) =         a8 * h3
        aa5(5,5) = 1._r8 + a8 * (h3 * alvl + h4)
        yy5(5)   =         a8 * y3  ! WGC row

        call matrix8_NxN(5,AA5,YY5,XX5,sing); if (sing) call sing_print(iland,4,5,aa5,yy5,glatw,glonw)

        if (xx5(1) <  eps_wxfer .and. &
            xx5(2) > -eps_wxfer .and. xx5(2) <= sfcwater_mass .and. &
            xx5(5) > -eps_wxfer) then ! test was successful

           wxfervc = xx5(1)
           wxfersc = xx5(2)
           hxfervc = xx5(3)
           hxfersc = xx5(4)
           wxfergc = xx5(5)
           transp  = 0.
           itest = 4
           go to 20
        endif

        ! If this test resulted in condensation onto vegetation,
        ! skip past tests 5 and 6

        if (xx5(1) < -eps_wxfer) then
           skiptest(5) = .true.
           skiptest(6) = .true.
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 7, 10, and 13

        if (xx5(2) < -eps_wxfer) then
           skiptest( 7) = .true.
           skiptest(10) = .true.
           skiptest(13) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 7

        if (xx5(2) <= sfcwater_mass) then
           skiptest(7) = .true.
        endif
 
     endif

     5 continue
     if (skiptest(5)) go to 6

! Test 5: If some sfcwater is present, and also if vegetation is wet, solve 5x5
!         system to test for partial evaporation from vegetation, partial
!         evaporation of sfcwater, and zero or positive evaporation from soil

     if (sfcwater_mass > 0. .and. iwetveg) then

        aa5(1,1) = 1._r8 + (a3 + a4) * (h1 * alvl + h4)
        aa5(1,2) =         (a3 + a4) * h4
        aa5(1,3) =         (a3 + a4) * h1
        aa5(1,4) = 0._r8
        aa5(1,5) =         (a3 + a4) * h4
        yy5(1)   =         (a3 + a4) * y1  ! WVC row

        aa5(2,1) =         a7 * h4
        aa5(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa5(2,3) = 0._r8
        aa5(2,4) =         a7 * h2
        aa5(2,5) =         a7 * h4
        yy5(2)   =         a7 * y2         ! WSC row

        aa5(3,1) =         a2 * h5 * alvl
        aa5(3,2) = 0._r8
        aa5(3,3) = 1._r8 + a2 * (h5 + h7)
        aa5(3,4) =         a2 * h7
        aa5(3,5) =         a2 * h5 * alvl
        yy5(3)   =         a2 * y4         ! HVC row

        aa5(4,1) = 0._r8
        aa5(4,2) =         a6 * h6 * alvl
        aa5(4,3) =         a6 * h7
        aa5(4,4) = 1._r8 + a6 * (h6 + h7)
        aa5(4,5) =         a6 * h6 * alvl
        yy5(4)   =         a6 * y5         ! HSC row

        aa5(5,1) =         a8 * h4
        aa5(5,2) =         a8 * h4
        aa5(5,3) = 0._r8
        aa5(5,4) =         a8 * h3
        aa5(5,5) = 1._r8 + a8 * (h3 * alvl + h4)
        yy5(5)   =         a8 * y3         ! WGC row

        call matrix8_NxN(5,AA5,YY5,XX5,sing); if (sing) call sing_print(iland,5,5,aa5,yy5,glatw,glonw)

        evap = xx5(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx5(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx5(2) > -eps_wxfer .and. xx5(2) <= sfcwater_mass .and. &
            xx5(5) > -eps_wxfer) then ! test was successful

           wxfervc = evap
           wxfersc = xx5(2)
           hxfervc = xx5(3)
           hxfersc = xx5(4)
           wxfergc = xx5(5)
           transp  = xx5(1) - evap
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

        if (xx5(2) < -eps_wxfer) then
           skiptest( 8) = .true.
           skiptest(11) = .true.
           skiptest(14) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 8

        if (xx5(2) <= sfcwater_mass) then
           skiptest(8) = .true.
        endif
 
     endif

     6 continue
     if (skiptest(6)) go to 7

! Test 6:  If some sfcwater is present, and given complete evaporation of any
!         veg_water, solve 5x5 system to test for zero or positive transpiration,
!         partial evaporation of sfcwater, and zero or positive evaporation from soil

     if (sfcwater_mass > 0.) then

        aa5(1,1) = 1._r8 + a4 * (h1 * alvl + h4)
        aa5(1,2) =         a4 * h4
        aa5(1,3) =         a4 * h1
        aa5(1,4) = 0._r8
        aa5(1,5) =         a4 * h4
        yy5(1)   =         a4 * (y1 - (h1 * alvl + h4) * veg_water)  ! WVC row

        aa5(2,1) =         a7 * h4
        aa5(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa5(2,3) = 0._r8
        aa5(2,4) =         a7 * h2
        aa5(2,5) =         a7 * h4
        yy5(2)   =         a7 * (y2 - h4 * veg_water)                ! WSC row

        aa5(3,1) =         a2 * h5 * alvl
        aa5(3,2) = 0._r8
        aa5(3,3) = 1._r8 + a2 * (h5 + h7)
        aa5(3,4) =         a2 * h7
        aa5(3,5) =         a2 * h5 * alvl
        yy5(3)   =         a2 * (y4 - h5 * alvl * veg_water)         ! HVC row

        aa5(4,1) = 0._r8
        aa5(4,2) =         a6 * h6 * alvl
        aa5(4,3) =         a6 * h7
        aa5(4,4) = 1._r8 + a6 * (h6 + h7)
        aa5(4,5) =         a6 * h6 * alvl
        yy5(4)   =         a6 * y5                                   ! HSC row

        aa5(5,1) =         a8 * h4
        aa5(5,2) =         a8 * h4
        aa5(5,3) = 0._r8
        aa5(5,4) =         a8 * h3
        aa5(5,5) = 1._r8 + a8 * (h3 * alvl + h4)
        yy5(5)   =         a8 * (y3 - h4 * veg_water)                ! WGC row

        call matrix8_NxN(5,AA5,YY5,XX5,sing); if (sing) call sing_print(iland,6,5,aa5,yy5,glatw,glonw)

        if (xx5(1) > -eps_wxfer .and. &
            xx5(2) > -eps_wxfer .and. xx5(2) <= sfcwater_mass .and. &
            xx5(5) > -eps_wxfer) then ! test was successful

           wxfervc = veg_water
           wxfersc = xx5(2)
           hxfervc = xx5(3)
           hxfersc = xx5(4)
           wxfergc = xx5(5)
           transp  = xx5(1)
           itest = 6
           go to 20
        endif

        ! If this test resulted in condensation onto the surface,
        ! skip past tests 9, 12, and 15

        if (xx5(2) < -eps_wxfer) then
           skiptest( 9) = .true.
           skiptest(12) = .true.
           skiptest(15) = .true.
        endif

        ! If this test resulted in less than complete evaporation of sfcwater,
        ! skip past test 9

        if (xx5(2) <= sfcwater_mass) then
           skiptest(9) = .true.
        endif
 
     endif


     7 continue
     if (skiptest(7)) go to 8

! Test 7:  Given complete evaporation of any sfcwater, solve 4x4 system to test
!         for condensation onto vegetation and zero or positive evaporation from soil

     aa4(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
     aa4(1,2) =         a1 * h4
     aa4(1,3) =         a1 * h1
     aa4(1,4) = 0._r8
     yy4(1)   =         a1 * (y1 - h4 * sfcwater_mass)                ! WVC row

     aa4(2,1) =         a8 * h4
     aa4(2,2) = 1._r8 + a8 * (h3 * alvl + h4)
     aa4(2,3) = 0._r8
     aa4(2,4) =         a8 * h3
     yy4(2)   =         a8 * (y3 - (h3 * alvl + h4) * sfcwater_mass)  ! WGC row

     aa4(3,1) =         a2 * h5 * alvl
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a2 * (h5 + h7)
     aa4(3,4) =         a2 * h7
     yy4(3)   =         a2 * y4                                       ! HVC row

     aa4(4,1) = 0._r8
     aa4(4,2) =         a6 * h6 * alvl
     aa4(4,3) =         a6 * h7
     aa4(4,4) = 1._r8 + a6 * (h6 + h7)
     yy4(4)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)         ! HSC row

     call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,7,4,aa4,yy4,glatw,glonw)

     if (xx4(1) <  eps_wxfer .and. &
         xx4(2) > -eps_wxfer) then ! test was successful

        wxfervc = xx4(1)
        wxfersc = sfcwater_mass
        hxfervc = xx4(3)
        hxfersc = xx4(4)
        wxfergc = xx4(2)
        transp  = 0.
        itest = 7
        go to 20
     endif

     ! If this test resulted in condensation onto vegetation,
     ! skip past tests 8 and 9

     if (xx4(1) < -eps_wxfer) then
        skiptest(8) = .true.
        skiptest(9) = .true.
     endif

     8 continue
     if (skiptest(8)) go to 9

! Test 8:  If vegetation is wet, and given complete evaporation of any sfcwater,
!         solve 4x4 system to test for partial evaporation from vegetation
!         and zero or positive evaporation from soil

     if (iwetveg) then

        aa4(1,1) = 1._r8 + (a3 + a4) * (h1 * alvl + h4)
        aa4(1,2) =         (a3 + a4) * h4
        aa4(1,3) =         (a3 + a4) * h1
        aa4(1,4) = 0._r8
        yy4(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass)         ! WVC row

        aa4(2,1) =         a8 * h4
        aa4(2,2) = 1._r8 + a8 * (h3 * alvl + h4)
        aa4(2,3) = 0._r8
        aa4(2,4) =         a8 * h3
        yy4(2)   =         a8 * (y3 - (h3 * alvl + h4) * sfcwater_mass)  ! WGC row

        aa4(3,1) =         a2 * h5 * alvl
        aa4(3,2) = 0._r8
        aa4(3,3) = 1._r8 + a2 * (h5 + h7)
        aa4(3,4) =         a2 * h7
        yy4(3)   =         a2 * y4                                       ! HVC row

        aa4(4,1) = 0._r8
        aa4(4,2) =         a6 * h6 * alvl
        aa4(4,3) =         a6 * h7
        aa4(4,4) = 1._r8 + a6 * (h6 + h7)
        yy4(4)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)         ! HSC row

        call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,8,4,aa4,yy4,glatw,glonw)

        evap = xx4(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx4(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx4(2) > -eps_wxfer) then ! test was successful

           wxfervc = evap
           wxfersc = sfcwater_mass
           hxfervc = xx4(3)
           hxfersc = xx4(4)
           wxfergc = xx4(2)
           transp  = xx4(1) - evap
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
!         of any veg_water, solve 4x4 system to test for zero or positive
!         transpiration and zero or positive evaporation from soil
 
     aa4(1,1) = 1._r8 + a4 * (h1 * alvl + h4)
     aa4(1,2) =         a4 * h4
     aa4(1,3) =         a4 * h1
     aa4(1,4) = 0._r8
     yy4(1)   =         a4 * (y1 - (h1 * alvl + h4) * veg_water - h4 * sfcwater_mass)  ! WVC row

     aa4(2,1) =         a8 * h4
     aa4(2,2) = 1._r8 + a8 * (h3 * alvl + h4)
     aa4(2,3) = 0._r8
     aa4(2,4) =         a8 * h3
     yy4(2)   =         a8 * (y3 - h4 * veg_water - (h3 * alvl + h4) * sfcwater_mass)  ! WGC row

     aa4(3,1) =         a2 * h5 * alvl
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a2 * (h5 + h7)
     aa4(3,4) =         a2 * h7
     yy4(3)   =         a2 * (y4 - h5 * alvl * veg_water)                              ! HVC row

     aa4(4,1) = 0._r8
     aa4(4,2) =         a6 * h6 * alvl
     aa4(4,3) =         a6 * h7
     aa4(4,4) = 1._r8 + a6 * (h6 + h7)
     yy4(4)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)                          ! HSC row

     call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,9,4,aa4,yy4,glatw,glonw)

     if (xx4(1) > -eps_wxfer .and. &
         xx4(2) > -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass
        hxfervc = xx4(3)
        hxfersc = xx4(4)
        wxfergc = xx4(2)
        transp  = xx4(1)
        itest = 9
        go to 20
     endif

     10 continue
     if (skiptest(10)) go to 11

! Test 10:  If some sfcwater is present, and given zero vapor flux with soil,
!          solve 4x4 system to test for condensation onto vegetation and
!          partial evaporation of sfcwater

     if (sfcwater_mass > 0.) then

        aa4(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
        aa4(1,2) =         a1 * h4
        aa4(1,3) =         a1 * h1
        aa4(1,4) = 0._r8
        yy4(1)   =         a1 * y1  ! WVC row

        aa4(2,1) =         a7 * h4
        aa4(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa4(2,3) = 0._r8
        aa4(2,4) =         a7 * h2
        yy4(2)   =         a7 * y2  ! WSC row

        aa4(3,1) =         a2 * h5 * alvl
        aa4(3,2) = 0._r8
        aa4(3,3) = 1._r8 + a2 * (h5 + h7)
        aa4(3,4) =         a2 * h7
        yy4(3)   =         a2 * y4  ! HVC row

        aa4(4,1) = 0._r8
        aa4(4,2) =         a6 * h6 * alvl
        aa4(4,3) =         a6 * h7
        aa4(4,4) = 1._r8 + a6 * (h6 + h7)
        yy4(4)   =         a6 * y5  ! HSC row

        call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,10,4,aa4,yy4,glatw,glonw)

        if (xx4(1) <  eps_wxfer .and. &
            xx4(2) > -eps_wxfer .and. xx4(2) < sfcwater_mass) then ! test was successful

           wxfervc = xx4(1)
           wxfersc = xx4(2)
           hxfervc = xx4(3)
           hxfersc = xx4(4)
           wxfergc = 0.
           transp  = 0.
           itest = 10
           go to 20
        endif

        ! If this test resulted in condensation onto vegetation,
        ! skip past tests 11 and 12

        if (xx4(1) < -eps_wxfer) then
           skiptest(11) = .true.
           skiptest(12) = .true.
        endif

     endif

     11 continue
     if (skiptest(11)) go to 12

! Test 11:  If some sfcwater is present and vegetation is wet, and given zero
!          vapor flux with soil, solve 4x4 system to test for partial
!          evaporation from vegetation and partial evaporation of sfcwater

     if (sfcwater_mass > 0. .and. iwetveg) then

        aa4(1,1) = 1._r8 + (a3 + a4) * (h1 * alvl + h4)
        aa4(1,2) =         (a3 + a4) * h4
        aa4(1,3) =         (a3 + a4) * h1
        aa4(1,4) = 0._r8
        yy4(1)   =         (a3 + a4) * y1  ! WVC row

        aa4(2,1) =         a7 * h4
        aa4(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa4(2,3) = 0._r8
        aa4(2,4) =         a7 * h2
        yy4(2)   =         a7 * y2         ! WSC row

        aa4(3,1) =         a2 * h5 * alvl
        aa4(3,2) = 0._r8
        aa4(3,3) = 1._r8 + a2 * (h5 + h7)
        aa4(3,4) =         a2 * h7
        yy4(3)   =         a2 * y4         ! HVC row

        aa4(4,1) = 0._r8
        aa4(4,2) =         a6 * h6 * alvl
        aa4(4,3) =         a6 * h7
        aa4(4,4) = 1._r8 + a6 * (h6 + h7)
        yy4(4)   =         a6 * y5         ! HSC row

        call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,11,4,aa4,yy4,glatw,glonw)

        evap = xx4(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx4(1) > -eps_wxfer .and. evap <= veg_water .and. &
            xx4(2) > -eps_wxfer .and. xx4(2) <= sfcwater_mass) then ! test was successful

           wxfervc = evap
           wxfersc = xx4(2)
           hxfervc = xx4(3)
           hxfersc = xx4(4)
           wxfergc = 0.
           transp  = xx4(1) - evap
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
!          veg_water, and given zero vapor flux with soil, solve 4x4 system
!          to test for zero or positive transpiration and partial evaporation of
!          any surface water
 
     if (sfcwater_mass > 0.) then
 
        aa4(1,1) = 1._r8 + a4 * (h1 * alvl + h4)
        aa4(1,2) =         a4 * h4
        aa4(1,3) =         a4 * h1
        aa4(1,4) = 0._r8
        yy4(1)   =         a4 * (y1 - (h1 * alvl + h4) * veg_water)  ! WVC row

        aa4(2,1) =         a7 * h4
        aa4(2,2) = 1._r8 + a7 * (h2 * alvl + h4)
        aa4(2,3) = 0._r8
        aa4(2,4) =         a7 * h2
        yy4(2)   =         a7 * (y2 - h4               * veg_water)  ! WSC row

        aa4(3,1) =         a2 * h5 * alvl
        aa4(3,2) = 0._r8
        aa4(3,3) = 1._r8 + a2 * (h5 + h7)
        aa4(3,4) =         a2 * h7
        yy4(3)   =         a2 * (y4 - h5 * alvl        * veg_water)  ! HVC row

        aa4(4,1) = 0._r8
        aa4(4,2) =         a6 * h6 * alvl
        aa4(4,3) =         a6 * h7
        aa4(4,4) = 1._r8 + a6 * (h6 + h7)
        yy4(4)   =         a6 * y5                                   ! HSC row

        call matrix8_4x4(AA4,YY4,XX4,sing); if (sing) call sing_print(iland,12,4,aa4,yy4,glatw,glonw)

        if (xx4(1) > -eps_wxfer .and. &
            xx4(2) > -eps_wxfer .and. xx4(2) <= sfcwater_mass) then ! test was successful

           wxfervc = veg_water
           wxfersc = xx4(2)
           hxfervc = xx4(3)
           hxfersc = xx4(4)
           wxfergc = 0.
           transp  = xx4(1)
           itest = 12
           go to 20
        endif

     endif

     13 continue
     if (skiptest(13)) go to 14

! Test 13:  Given zero vapor flux with soil and given complete evaporation of
!          any sfcwater, solve 3x3 system to test for condensation onto vegetation

     aa3(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
     aa3(1,2) =         a1 * h1
     aa3(1,3) = 0._r8
     yy3(1)   =         a1 * (y1 - h4 * sfcwater_mass)         ! WVC row

     aa3(2,1) =         a2 * h5 * alvl
     aa3(2,2) = 1._r8 + a2 * (h5 + h7)
     aa3(2,3) =         a2 * h7
     yy3(2)   =         a2 * y4                                ! HVC row

     aa3(3,1) = 0._r8
     aa3(3,2) =         a6 * h7
     aa3(3,3) = 1._r8 + a6 * (h6 + h7)
     yy3(3)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)  ! HSC row 

     call matrix8_3x3(AA3,YY3,XX3,sing); if (sing) call sing_print(iland,13,3,aa3,yy3,glatw,glonw)

     if (xx3(1) < eps_wxfer) then ! test was successful

        wxfervc = xx3(1)
        wxfersc = sfcwater_mass
        hxfervc = xx3(2)
        hxfersc = xx3(3)
        wxfergc = 0.
        transp  = 0.
        itest = 13
        go to 20
     endif

     ! If this test resulted in condensation onto vegetation,
     ! skip past tests 14 and 15

     if (xx3(1) < -eps_wxfer) then
        skiptest(14) = .true.
        skiptest(15) = .true.
     endif

     14 continue

     if (skiptest(14)) go to 15

! Test 14:  If vegetation is wet, and given zero vapor flux with soil and
!          complete evaporation of any sfcwater, solve 3x3 system to test for
!          partial evaporation from vegetation
 
     if (iwetveg) then

        aa3(1,1) = 1._r8 + (a3 + a4) * (h1 * alvl + h4)
        aa3(1,2) =         (a3 + a4) * h1
        aa3(1,3) = 0._r8
        yy3(1)   =         (a3 + a4) * (y1 - h4 * sfcwater_mass)  ! WVC row

        aa3(2,1) =         a2 * h5 * alvl
        aa3(2,2) = 1._r8 + a2 * (h5 + h7)
        aa3(2,3) =         a2 * h7
        yy3(2)   =         a2 * y4                                ! HVC row

        aa3(3,1) = 0._r8
        aa3(3,2) =         a6 * h7
        aa3(3,3) = 1._r8 + a6 * (h6 + h7)
        yy3(3)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)  ! HSC row 

        call matrix8_3x3(AA3,YY3,XX3,sing); if (sing) call sing_print(iland,14,3,aa3,yy3,glatw,glonw)

        evap = xx3(1) * a3 / (a3 + a4) ! evaporation of veg_water

        if (xx3(1) > -eps_wxfer .and. evap <= veg_water) then ! test was successful

           wxfervc = evap
           wxfersc = sfcwater_mass
           hxfervc = xx3(2)
           hxfersc = xx3(3)
           wxfergc = 0.
           transp  = xx3(1) - evap
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
 
     aa3(1,1) = 1._r8 + a4 * (h1 * alvl + h4)
     aa3(1,2) =         a4 * h1
     aa3(1,3) = 0._r8
     yy3(1)   =         a4 * (y1 - (h1 * alvl + h4) * veg_water - h4 * sfcwater_mass)  ! WVC row

     aa3(2,1) =         a2 * h5 * alvl
     aa3(2,2) = 1._r8 + a2 * (h5 + h7)
     aa3(2,3) =         a2 * h7
     yy3(2)   =         a2 * (y4 - h5 * alvl        * veg_water)                        ! HVC row

     aa3(3,1) = 0._r8
     aa3(3,2) =         a6 * h7
     aa3(3,3) = 1._r8 + a6 * (h6 + h7)
     yy3(3)   =         a6 * (y5 - h6 * alvl * sfcwater_mass)                           ! HSC row 

     call matrix8_3x3(AA3,YY3,XX3,sing); if (sing) call sing_print(iland,15,3,aa3,yy3,glatw,glonw)

     if (xx3(1) > -eps_wxfer) then ! test was successful

        wxfervc = veg_water
        wxfersc = sfcwater_mass
        hxfervc = xx3(2)
        hxfersc = xx3(3)
        wxfergc = 0.
        transp  = xx3(1)
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

     cantemp = cantemp + (hxfervc + hxfersc - cp * sxfer_t) * hcapcani

     veg_temp = veg_temp + (radveg - hxfervc - (wxfervc + transp) * alvl) * hcapvegi

     ! Energy_per_m2 may pertain to sfcwater only, soil only, or both combined.  In cases
     ! where it pertains to sfcwater only, wxfergc will be zero.

     energy_per_m2 = energy_per_m2 + radsfc - hxfersc - (wxfersc + wxfergc) * alvl

     canrrv = canrrv + (wxfervc + transp + wxfersc + wxfergc - sxfer_r) * canairi

     veg_water = max(0.,veg_water - wxfervc)

     soil_water(nzg) = soil_water(nzg) - dslzi(nzg) * wxfergc * .001

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
  use misc_coms, only: io6

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

     veg_rough = veg_height * (1. - bz * exp(-hz * veg_tai))
     veg_albedo = albv_green(leaf_class) * green_frac  &
                + albv_brown(leaf_class) * (1. - green_frac)
     veg_fracarea = veg_frac(leaf_class) * (1. - exp(-extinc_veg * veg_tai))

  endif         

  end subroutine vegndvi

!===============================================================================

  subroutine sing_print(iland,ieqn,nsize,aa,yy,glatw,glonw)

  use misc_coms,   only: io6
  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: iland, ieqn, nsize

  real(r8), intent(in) :: aa(nsize,nsize), yy(nsize)
  real,     intent(in) :: glatw, glonw

  integer :: irow

  print*, 'singular matrix '
  print*, 'iland, ieqn, nsize ', iland, ieqn, nsize
  print*, 'glatw, glonw ', glatw, glonw

  do irow = 1,nsize
     write(io6,'(a,i5,5e15.5)') 'row, aa, yy ',irow,aa(irow,1:nsize),yy(irow)
  enddo

  stop 'singular matrix '

  end subroutine sing_print

!===============================================================================

  subroutine fast_canopy(iland, iwsfc, nlsw1, icomb,                         &
                         sfcwat_nud,    sfctemp_nud,     fracliq_nud,        &
                         sfcwater_mass, sfcwater_energy, sfcwater_depth,     &
                         energy_per_m2, sfcwater_tempk,  sfcwater_fracliq,   &
                         soil_water,    soil_energy,     specifheat_drysoil, &
                         soil_tempk,    soil_fracliq                         )

  use leaf_coms,      only: dt_leaf, snowmin_expl
  use mem_sfcg,       only: sfcg
  use mem_land,       only: nzg
  use misc_coms,      only: io6
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

  integer :: k        ! loop index over soil layers
  integer :: iveg     ! flag for exposed vegetation (0=no, 1=yes)

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
