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

 subroutine canopy(iwl, nlsw1, icomb, leaf_class, ktrans, ntext_soil, &
                   soil_water, soil_energy, soil_tempk, soil_fracliq, &
                   sfcwater_mass, sfcwater_energy, sfcwater_depth,    &
                   energy_per_m2, sfcwater_tempk, sfcwater_fracliq,   &
                   veg_height, veg_rough, veg_tai, veg_lai,           &
                   hcapveg, can_depth,                                &
                   rhos, vels, prss, ustar, pcpg, qpcpg, dpcpg,       &
                   rshort, rshort_v, rshort_g,                        &
                   rlong_v, rlong_s, rlong_g, sxfer_t, sxfer_r,       &
                   snowfac, vf, surface_ssh, ground_shv,              &
                   veg_water, veg_temp, cantemp, canshv,              &
                   transp, stom_resist, snowmin_expl,                 &
                   glatw, glonw, flag_vg                              ) 

 use leaf_coms, only: nzg, soil_rough, dt_leaf, dslz, dslzi, slzt, &
                      slmstsh0_ch, slmstsh0_vg, kroot, rcmin,      &
                      slcpd, slcons_ch, slcons_vg

 use consts_coms,    only: cp, vonk, alvl, cliq, cice, alli, rvap, r8
 use misc_coms,      only: io6
 use mem_leaf,       only: itab_wl
 use leaf4_sfcwater, only: sfcwater_soil_comb
 use leaf4_soil,     only: soil_wat2pot
 use matrix,         only: matrix8_2x2, matrix8_3x3, matrix8_4x4

 implicit none

 integer, intent(in)    :: iwl        ! index of current land cell
 integer, intent(in)    :: nlsw1      ! k index of top active sfcwater layer
 integer, intent(inout) :: icomb      ! implicit heat balance flag [0=no, 1=yes]
 integer, intent(in)    :: leaf_class ! leaf class (vegetation class)
 integer, intent(out)   :: ktrans     ! k index of soil layer supplying transp

 integer, intent(in) :: ntext_soil  (nzg) ! soil textural class
 real, intent(inout) :: soil_water  (nzg) ! soil water content [vol_water/vol_tot]
 real, intent(inout) :: soil_energy (nzg) ! soil energy [J/m^3]
 real, intent(in)    :: soil_tempk  (nzg) ! soil temp [K]
 real, intent(in)    :: soil_fracliq(nzg) ! fraction of soil moisture in liquid phase

 real, intent(inout) :: sfcwater_mass    ! surface water mass [kg/m^2]
 real, intent(inout) :: sfcwater_energy  ! surface water energy [J/kg]
 real, intent(inout) :: sfcwater_depth   ! surface water depth [m]
 real, intent(inout) :: energy_per_m2    ! surface water energy [J/m^2]
 real, intent(inout) :: sfcwater_tempk   ! surface water temperature [K]
 real, intent(inout) :: sfcwater_fracliq ! fraction of sfc water in liquid phase

 real, intent(in) :: veg_height  ! veg height [m]
 real, intent(in) :: veg_rough   ! veg roughess height [m]
 real, intent(in) :: veg_tai     ! veg total area index
 real, intent(in) :: veg_lai     ! veg leaf area index
 real, intent(in) :: hcapveg     ! veg heat capacity [J/(m^2 K)]
 real, intent(in) :: can_depth   ! canopy depth for heat and vap capacity [m]
 real, intent(in) :: rhos        ! atmospheric air density [kg/m^3]
 real, intent(in) :: vels        ! surface wind speed [m/s]
 real, intent(in) :: prss        ! air pressure [Pa]
 real, intent(in) :: ustar       ! friction velocity [m/s]
 real, intent(in) :: pcpg        ! added precip mass this leaf timestep [kg/m^2]
 real, intent(in) :: qpcpg       ! added precip energy this leaf timestep [J/m^2]
 real, intent(in) :: dpcpg       ! added precip depth this leaf timestep [m]
 real, intent(in) :: rshort      ! downward sfc s/w rad flux [W/m^2]
 real, intent(in) :: rshort_v    ! s/w rad flux absorbed by veg [W/m^2]
 real, intent(in) :: rshort_g    ! s/w rad flux absorbed by gnd [W/m^2]
 real, intent(in) :: rlong_v     ! l/w rad flux absorbed by veg [W/m^2]
 real, intent(in) :: rlong_s     ! l/w rad flux absorbed by sfc [W/m^2]
 real, intent(in) :: rlong_g     ! l/w rad flux absorbed by gnd [W/m^2]
 real, intent(in) :: sxfer_t     ! canopy-to-atm sens heat xfer this step [kg_air K/m^2]
 real, intent(in) :: sxfer_r     ! canopy-to-atm vapor xfer this step [kg_vap/m^2]
 real, intent(in) :: snowfac     ! fractional veg burial by snowcover
 real, intent(in) :: vf          ! fractional coverage of non-buried part of veg
 real, intent(in) :: surface_ssh ! surface sat spec hum [kg_vap/kg_air]
 real, intent(in) :: ground_shv  ! soil vapor spec hum [kg_vap/kg_air]

 real, intent(inout) :: veg_water   ! veg sfc water content [kg/m^2]
 real, intent(inout) :: veg_temp    ! veg temp [K]
 real, intent(inout) :: cantemp     ! canopy air temp [K]
 real, intent(inout) :: canshv      ! canopy air vapor spec hum [kg_vap/kg_air]
 real, intent(inout) :: stom_resist ! veg stomatal resistance [s/m]

 real, intent(out) :: transp       ! transpiration xfer this LEAF timestep [kg_vap/m^2]
 real, intent(in)  :: snowmin_expl ! min sfcwater mass for explicit heat xfer [kg/m^2]
 real, intent(in)  :: glatw        ! Latitude of land cell 'center' [deg]
 real, intent(in)  :: glonw        ! Longitude of land cell 'center' [deg]

 logical, intent(in) :: flag_vg

! Local parameters

 real, parameter :: exar = 2.5  ! for computing rasveg
 real, parameter :: covr = 2.16 ! scaling tai value for computing wtveg
 real, parameter :: c1 = 116.6  ! for computing rb

!     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar)) 
!     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
!     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
!     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
!     total expression inside the radical in Eq. 3.37.
!bob      parameter(exar=3.5,covr=2.16,c1=98.8)

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
 integer :: nts      ! soil textural class for current soil layer
 integer :: iveg     ! flag for exposed vegetation (0=no, 1=yes)
 integer :: iwetsfc  ! flag for wet surface (0=no, 1=yes)
 integer :: iwetveg  ! flag for wet vegetation (0=no, 1=yes)

 real :: hxfersc  ! sfc-to-can_air heat xfer this step [J/m^2]
 real :: wxfersc  ! sfc-to-can_air vap xfer this step [kg_vap/m^2]
 real :: hxfervc  ! veg-to-can_air heat xfer this step [J/m^2]
 real :: wxfervc  ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
 real :: wshed    ! water shed from veg this LEAF timestep [kg/m^2]
 real :: qwshed   ! water energy shed from veg this LEAF timestep [J/m^2]

 real :: factv       ! for computing rasveg
 real :: aux         ! for computing rasveg
 real :: rasveg      ! full-veg value of rd [s/m]
 real :: wtveg       ! weighting of rasveg in computing rd
 real :: rasgnd      ! not used
 real :: c3          ! veg_sfc-to-can_air vapor density difference [kg_vap/m^3]
 real :: fracliqv    ! fraction of veg surface water in liquid phase
 real :: fthi        ! high-temp environ factor for stomatal resist
 real :: ftlo        ! low-temp environ factor for stomatal resist
 real :: frad        ! s/w radiative environ factor for stomatal resist
 real :: fswp        ! soil water potential environ factor for stomatal resist
 real :: fvpd        ! vap press deficit environ factor for stomatal resist
 real :: qwtot       ! total internal energy of veg plus new precip on veg [J/m^2]
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
 real :: tvegk       ! vegetation temperature (K)
 real :: wtroot      ! not used
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
 real :: sfc_fracliq ! temporary variable for surface water liquid fraction

 real :: sfc_rhovs   ! surface saturation water vapor density [kg_vap/m^3]
 real :: sfc_rhovs1  ! sfc_rhovs at 1 K warmer [kg_vap/m^3]
 real :: sfc_rhovsp  ! sfc_rhovs derivative with respect to temperature

 real :: gnd_rhov    ! ground sfc evaporative water vapor density [kg_vap/m^3]
 real :: gnd_rhov1   ! gnd_rhov at 1 K warmer [kg_vap/m^3]
 real :: gnd_rhovp   ! gnd_rhov derivative with respect to temperature

 real :: wcap_soil   ! unfilled soil water capacity [kg/m^2]
 real :: wcap_rate   ! potential infiltration at sat hyd conductivity [kg/m^2]
 real :: wcap_both   ! minimum of (wcap_soil, wcap_rate) [kg/m^2]
 real :: wcap_min    ! minimum surface water water [kg/m^2]

 real :: waterfrac   ! Fraction of water saturation of top soil layer
 real :: vw          ! same as veg_water (shorter word)
 real :: vwm         ! minimum value of vw

 real :: water_frac_ul
 real :: water_frac
 real :: head        ! soil water potential [m]

 real :: psi, headp  ! unused but required subroutine arguments

 real :: a1, a2, a3, a4, a5, a6, a134
 real :: y1, y2, y3, y4, y5, y23
 real :: h1, h2, h3, h4, h5, h6, h7, h23
 real :: specvol

 real(r8) :: aa2(2,2), xx2(2),   yy2(2)  ! 2x2 matrix equation terms
 real(r8) :: aa3(3,3), xx3(3,3), yy3(3)  ! 3x3 matrix equation terms
 real(r8) :: aa4(4,4), xx4(4,6), yy4(4)  ! 4x4 matrix equation terms

 logical :: sing

 real, external :: rhovsil ! function to compute sat vapor density
 real, external :: eslf    ! function to compute sat vapor pressure

 integer :: itest  ! test number
 integer :: ieqn   ! matrix equation number; loop index

! Set IVEG flag = 1 if there is sufficient exposed vegetation; set to 0 otherwise

 iveg = 0
 if (vf > .001) iveg = 1

 if (iveg == 1) then  ! Vegetation case

! Vegetation is sufficiently abundant and not covered by snow.
! COMPUTE RDI WITH VEGETATION INFLUENCE

! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
! the values (zognd,0) when veg covered.

    zognd = soil_rough
    zoveg = veg_rough * (1. - snowfac) + zognd * snowfac
    zveg  = veg_height * (1. - snowfac)
    zdisp = zveg * .63
!bob      rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd)
!bob     +      / (vonk * vonk * vels)

! Aerodynamic resistance (rd) between surface and canopy air are weighted
! between areas without and with vegetation.

!bob   factv  = log(zts / zoveg) / (vonk * vonk * vels)
    factv  = 1. / (vonk * ustar)
    aux    = exp(exar * (1. - (zdisp + zoveg) / zveg))
    rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
    wtveg  = max(0.,min(1., 1.1 * veg_tai / covr))
    rdi    = ustar / (5. * (1. - wtveg) + ustar * rasveg * wtveg)
   
! TAI and LAI reduced by ground snowcover (for low vegetation)

    stai = veg_tai * (1. - snowfac)
    slai = veg_lai * (1. - snowfac)

! Evaluate vegetation Celsius temperature

    tvegc = veg_temp - 273.15

! Add any precipitation intercepted mass and energy to vegetation surface

    veg_water = veg_water + pcpg * vf
    qwtot = hcapveg * tvegc + qpcpg * vf

! Compute equilbrium temperature of veg + precipitation
! (Ignore mass and energy of water already on vegetation before precipitation)

    call qwtk(qwtot,pcpg * vf,hcapveg,veg_temp,fracliqv)
    tvegc = veg_temp - 273.15

! Shed any excess intercepted precipitation and its energy

    if (veg_water > .2 * stai) then
       wshed = veg_water - .2 * stai

       if (fracliqv <= .0001) then
          qwshed = cice * tvegc * wshed
       else
          qwshed = (cliq * tvegc + fracliqv * alli) * wshed
       endif

       veg_water = veg_water - wshed
    else
       wshed  = 0.
       qwshed = 0.
    endif

! Check for presence of veg_water and set iwetveg flag

    wcap_min = dt_leaf * 1.e-7 ! 1.e-7 is 1 mm pcp/dew in 115 days   [kg/m^2]

    if (veg_water > wcap_min) then
       iwetveg = 1
    else
       iwetveg = 0
       veg_water = 0.
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

    do k = kroot(leaf_class),nzg
       nts = ntext_soil(k)

! Soil water potential

       call soil_wat2pot(iwl,nts,flag_vg,soil_water(k), slzt(k), &
                         water_frac_ul, water_frac, psi, head, headp)

! Find layer in root zone with highest head (unless mostly ice)

       if (swp < head .and. soil_fracliq(k) > .5) then 
          swp = head
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

       e_can = canshv * rhos * rvap * cantemp

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
! (requires dt_leaf <= 900.) 

       rc = 1. / (1. / rc + .0011 * dt_leaf * (1. / rc_inf - 1. / rc))

! Limit maximum transpiration to be <= 500 W/m^2 by increasing rc if necessary.

       transp_test = alvl * veg_lai * (veg_rhovs - rhos * canshv) / (rb + rc)
       if (transp_test > 500.) then
          rc = (rb + rc) * transp_test * .002 - rb
       endif      
      
    endif
      
    stom_resist = rc

 else  ! No-vegetation case

! Aerodynamic conductance for bare soil or snow based on Garratt.

    rdi = .2 * ustar

    wshed  = 0.
    qwshed = 0.
    ktrans = 0

 endif

! Canopy air quantities

 can_rhov = canshv * rhos
 canair = rhos * can_depth
 hcapcan = cp * canair

! Apply contributions to sfcwater mass, depth, and energy_per_m2 from
! precipitation and shedding of moisture from vegetation.  Sfcwater
! quantities represent the NLSW1 layer in the sfcwater arrays in case
! nlev_sfcwater = 0.

 sfcwater_mass  = sfcwater_mass  + pcpg  * (1. - vf) + wshed
 sfcwater_depth = sfcwater_depth + dpcpg * (1. - vf) + wshed * .001
 energy_per_m2  = energy_per_m2  + qpcpg * (1. - vf) + qwshed

! Check whether nlsw1 = 1 (lowest sfcwater layer), in which case sfcwater
! content may be very small which requires adjustments to be made.
! First, set implicit energy exchange flag (icomb) to default of 0
! and surface wetness flag (iwetsfc) to default of 1.

 icomb = 0
 iwetsfc = 1  ! Any surface evaporation will be from sfcwater
 nts = ntext_soil(nzg)

 if (nlsw1 == 1) then

! If sfcwater mass is very small or zero, or if it is mostly liquid and and
! scarce enough to infiltrate this timestep, transfer its mass and energy to
! top soil layer.  First compute water capacity available in soil for new
! infiltration, water capacity that could infiltrate based on saturated
! hydraulic conductivity, and minimum water capacity below which infiltration
! is forced even if sfcwater is frozen.  This helps to avoid situations of
! over-evaporation of sfcwater; for simplicity, evaporation is not strictly
! limited to available sfcwater mass in matrix solution below.

    if (flag_vg) then
       wcap_soil = dslz(nzg) * (slmstsh0_vg(nts) - soil_water(nzg)) * 1000. ! [kg/m^2]
       wcap_rate = dt_leaf * slcons_vg(nts) * 1000.                         ! [kg/m^2]
    else
       wcap_soil = dslz(nzg) * (slmstsh0_ch(nts) - soil_water(nzg)) * 1000. ! [kg/m^2]
       wcap_rate = dt_leaf * slcons_ch(nts) * 1000.                         ! [kg/m^2]
    endif

    wcap_both = min(wcap_soil, wcap_rate)                           ! [kg/m^2]

    wcap_min  = dt_leaf * 1.e-6 ! 1.e-6 is 1 mm pcp/dew in 11 days    [kg/m^2]

    if (sfcwater_mass < wcap_min .or. &
       (sfcwater_mass < wcap_both .and. sfcwater_fracliq > .5)) then

       soil_water (nzg) = soil_water (nzg) + dslzi(nzg) * sfcwater_mass * .001
       soil_energy(nzg) = soil_energy(nzg) + dslzi(nzg) * energy_per_m2

       sfcwater_mass   = 0.
       sfcwater_energy = 0.
       sfcwater_depth  = 0.
       energy_per_m2   = 0.

       sfcwater_tempk   = soil_tempk(nzg)
       sfcwater_fracliq = 0.

       icomb = 1
       iwetsfc = 0   ! Any surface evaporation will be from soil

    elseif (sfcwater_mass < snowmin_expl .or. sfcwater_fracliq > .5) then

! If nlsw1 = 1 (lowest sfcwater layer), and sfcwater mass is present but
! below limiting stability value for explicit heat transfer, or if sfcwater
! is mostly liquid, perform implicit thermal balance between sfcwater and
! the top soil layer.  

       call sfcwater_soil_comb(iwl,ntext_soil(nzg), &
                            soil_water(nzg),soil_energy(nzg), &
                            sfcwater_mass,energy_per_m2, &
                            sfcwater_tempk,sfcwater_fracliq)

       icomb = 1

    endif

 endif ! (nlsw1 == 1)

 if (icomb == 0) then

! If sfcwater and soil(nzg) temperatures are NOT implicitly coupled, diagnose
! sfc heat capacity for sfcwater alone.  First, diagnose sfcwater temperature
! and liquid fraction.  Use qwtk instead of qtk because energy_per_m2 is used
! instead of sfcwater_energy. ("dryhcap" = 100 is very small value)

    call qwtk(energy_per_m2,sfcwater_mass,100.,sfcwater_tempk,sfcwater_fracliq)

    if     (energy_per_m2 < 0.) then
       hcapsfc = sfcwater_mass * cice
    elseif (energy_per_m2 > sfcwater_mass * alli) then
       hcapsfc = sfcwater_mass * cliq
    else
       hcapsfc = sfcwater_mass * (alli / 1.)  ! Assume small but nonzero dT/dE
    endif

 else

! If sfcwater and soil(nzg) temperatures ARE implicitly coupled, diagnose
! sfc heat capacity for both together.

    sfc_energy = energy_per_m2 + soil_energy(nzg) * dslz(nzg)
    sfc_wmass  = sfcwater_mass + soil_water(nzg)  * dslz(nzg) * 1.e3
    sfc_soilhc = slcpd(nts) * dslz(nzg)

    if     (sfc_energy < 0.) then
       hcapsfc = sfc_wmass * cice + sfc_soilhc
    elseif (sfc_energy > sfc_wmass * alli) then
       hcapsfc = sfc_wmass * cliq + sfc_soilhc
    else
       hcapsfc = sfc_wmass * (alli / 1.) + sfc_soilhc ! Assume small positive dT/dE
    endif

 endif

! Sfcwater saturation vapor density and derivative

 sfc_tempk = sfcwater_tempk

 sfc_tempk1 = sfc_tempk + 1.

 sfc_rhovs  = rhovsil(sfc_tempk -273.15)
 sfc_rhovs1 = rhovsil(sfc_tempk1-273.15)
 sfc_rhovsp = sfc_rhovs1 - sfc_rhovs

! Reduced saturation vapor density and gradient of soil during evaporation

 call grndvap_ab(iwl,nts,sfc_tempk, soil_water(nzg),can_rhov,sfc_rhovs, gnd_rhov, flag_vg)
 call grndvap_ab(iwl,nts,sfc_tempk1,soil_water(nzg),can_rhov,sfc_rhovs1,gnd_rhov1,flag_vg)

 gnd_rhovp = gnd_rhov1 - gnd_rhov

! The remainder of subroutine canopy sets up and solves a linear system of
! equations using the trapezoidal-implicit method.  The set of equations 
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

 canairi  = 1. / canair
 hcapcani = 1. / hcapcan
 hcapsfci = 1. / hcapsfc

 radsfc = dt_leaf * (rlong_s + rshort_g + rlong_g)

 a5 = dt_leaf * rdi                             ! sfc vap xfer coef
 a6 = cp * rhos * a5                            ! sfc heat xfer coef

 y2 = sfc_rhovs - rhos * canshv &
    + fcn * (sfc_rhovsp * radsfc * hcapsfci + rhos * sxfer_r * canairi)

 y3 = gnd_rhov - rhos * canshv &
    + fcn * (gnd_rhovp * radsfc * hcapsfci + rhos * sxfer_r * canairi)

 y5 = sfc_tempk - cantemp &
    + fcn * (radsfc * hcapsfci + cp * sxfer_t * hcapcani)

 h2 = fcn * sfc_rhovsp * hcapsfci
 h3 = fcn * gnd_rhovp  * hcapsfci
 h4 = fcn * rhos * canairi

 h6 = fcn * hcapsfci
 h7 = fcn * hcapcani

 if (iveg == 0) then

! Case with NO VEGETATION

! Fill arrays for matrix solution.

! Test 1 for WET SURFACE CONDITION and/or CONDENSATION

! WSC row

    aa2(1,1) = 1._r8 + a5 * (h2 * alvl + h4)
    aa2(1,2) =         a5 * h2
    yy2(1)   =         a5 * y2

! HSC row

    aa2(2,1) =         a6 * h6 * alvl
    aa2(2,2) = 1._r8 + a6 * (h6 + h7)
    yy2(2)   =         a6 * y5

    call matrix8_2x2(AA2,YY2,XX2,sing)

    if (sing) call sing_print(iwl,1,2,aa2,yy2,glatw,glonw)

    if (iwetsfc == 1 .or. xx2(1) <= 0._r8) then
       wxfersc = xx2(1)
       hxfersc = xx2(2)

       itest = 51
       go to 10
    endif

! Test 2 for DRY SURFACE EVAPORATION

    aa2(1,1) = 1._r8 + a5 * (h3 * alvl + h4)
    aa2(1,2) =         a5 * h3
    yy2(1)   =         a5 * y3

    call matrix8_2x2(AA2,YY2,XX2,sing)

    if (sing) call sing_print(iwl,2,2,aa2(:,2),yy2,glatw,glonw)

    if (xx2(1) >= 0._r8) then
       wxfersc = xx2(1)
       hxfersc = xx2(2)

       itest = 52
       go to 10
    endif

! If we got here, surface is dry and surface vapor flux is zero

    wxfersc = 0.
    hxfersc = yy2(2) / aa2(2,2)

    itest = 53

    10 continue

! Update components

    cantemp = cantemp + (hxfersc - cp * sxfer_t) * hcapcani

    energy_per_m2 = energy_per_m2 + radsfc - hxfersc - wxfersc * alvl

    canshv = canshv + (wxfersc - sxfer_r) * canairi

    specvol = .001
    if (wxfersc > 0. .and. sfcwater_mass > 1.e-3) specvol = sfcwater_depth / sfcwater_mass

    sfcwater_mass = sfcwater_mass - wxfersc

    sfcwater_depth = sfcwater_depth - wxfersc * specvol

    transp = 0.

 else

! Case WITH VEGETATION

    hcapvegi = 1. / hcapveg

    radveg = dt_leaf * (rshort_v + rlong_v)

! veg_water fractional coverage

    sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66666666)

    a1 = dt_leaf * 2. * stai / rb                  ! veg vap xfer coef (dew formation)
    a2 = cp * rhos * a1                            ! veg heat xfer coef
    a3 = a1 * sigmaw                               ! veg vap xfer coef (veg_water evap)
    a4 = dt_leaf * slai * (1.-sigmaw) / (rb + rc)  ! veg vap xfer coef (transp)

! Auxiliary quantities

    y1 = veg_rhovs - rhos * canshv &
       + fcn * (veg_rhovsp * radveg * hcapvegi + rhos * sxfer_r * canairi)

    y4 = veg_temp - cantemp &
       + fcn * (radveg * hcapvegi + cp * sxfer_t * hcapcani)

    h1 = fcn * veg_rhovsp * hcapvegi

    h5 = fcn * hcapvegi

! Fill arrays for matrix solution (trapezoidal implicit method) to balance
! vapor and heat fluxes between canopy air, vegetation, and the surface.
! For vegetation, there are 3 distinct situations to consider:

! 1. Condensation onto vegetation (Vc)
! 2. Evaporation from vegetation with veg water not depleted (Vew)
! 3. Evaporation from vegetation with veg water depleted or zero (Ve)

! For the surface, there are 4 distinct situations to consider:

! 1. Wet surface (iwet = 1) (Gw)
! 2. Dry surface (iwet = 0) with condensation (Gc)
! 3. Dry surface (iwet = 0) with evaporation (Ge)
! 4. Dry surface (iwet = 0) with zero vapor flux (Go)

! (Zero vapor flux from the soil occurs when soil is too warm for condensation
! and too dry for evaporation.)

! There are 12 possible pairs containing one item from each list.  However, 
! equations for a wet soil surface (either condensation or evaporation) are
! identical to those with condensation onto a dry soil surface, so there are
! 9 distinct matrix equations, only one of which is the correct description
! for a particular land grid cell and timestep.

! We choose to perform all matrix computations first and to determine which is
! the correct choice after they are complete.  This leads to extra computation
! in most cases, but it gives code that is easier to understand and debug.
! This approach also provides all solutions together which facilitates
! understanding of the system behavior.

! Rather than always solving all 9 matrix equations, however, we first check
! if iwetsfc = 1, in which case we can eliminate 3 equations.  We also check if
! iwetveg = 0, in which case we can eliminate 1 equation for iwetsfc = 0 and
! one equation for iwetsfc = 1.

    do ieqn = 1,6

       if (iwetsfc == 1 .and. ieqn >= 4) cycle

       if (iwetveg == 0 .and. (ieqn == 2 .or. ieqn == 5)) cycle

       ! EVT row

       if (ieqn == 1 .or. ieqn == 4) a134 = a1
       if (ieqn == 2 .or. ieqn == 5) a134 = a3 + a4
       if (ieqn == 3 .or. ieqn == 6) a134 = a4

       aa4(1,1) = 1._r8 + a134 * (h1 * alvl + h4)
       aa4(1,2) =         a134 * h4
       aa4(1,3) =         a134 * h1
       aa4(1,4) = 0._r8
       yy4(1)   =         a134 * y1

       ! WSC row

       if (ieqn == 1 .or. ieqn == 4) then

          h23 = h2
          y23 = y2
          if (ieqn == 4) then
             h23 = h3
             y23 = y3
          endif
        
          aa4(2,1) =         a5 * h4
          aa4(2,2) = 1._r8 + a5 * (h23 * alvl + h4)
          aa4(2,3) = 0._r8
          aa4(2,4) =         a5 * h23
          yy4(2)   =         a5 * y23

       endif

       if (ieqn == 1) then

       ! HVC row

          aa4(3,1) =         a2 * h5 * alvl
          aa4(3,2) = 0._r8
          aa4(3,3) = 1._r8 + a2 * (h5 + h7)
          aa4(3,4) =         a2 * h7
          yy4(3)   =         a2 * y4

       ! HSC row

          aa4(4,1) = 0._r8
          aa4(4,2) =         a6 * h6 * alvl
          aa4(4,3) =         a6 * h7
          aa4(4,4) = 1._r8 + a6 * (h6 + h7)
          yy4(4)   =         a6 * y5

       endif

! Solve 4x4 matrix equation

       call matrix8_4x4(AA4,YY4,XX4(:,ieqn),sing)

       if (sing) call sing_print(iwl,ieqn,4,aa4,yy4,glatw,glonw)

! For dry soil condition, also solve 3x3 matrix equation (zero vapor flux)

       if (ieqn >= 4) then

          aa3(1,1) = aa4(1,1)
          aa3(1,2) = aa4(1,3)
          aa3(1,3) = aa4(1,4)
          aa3(2,1) = aa4(3,1)
          aa3(2,2) = aa4(3,3)
          aa3(2,3) = aa4(3,4)
          aa3(3,1) = aa4(4,1)
          aa3(3,2) = aa4(4,3)
          aa3(3,3) = aa4(4,4)

          yy3(1) = yy4(1)
          yy3(2) = yy4(3)
          yy3(3) = yy4(4)
                
          call matrix8_3x3(AA3,YY3,XX3(:,ieqn-3),sing)

          if (sing) call sing_print(iwl,ieqn,3,aa3,yy3,glatw,glonw)

       endif

    enddo

! Proceed to test results; first, consider cases for iwetsfc = 1

    vw = veg_water
    vwm = -1.e-7 * dt_leaf  ! minimum veg water; 1.e-7 is 1 mm pcp/dew in 115 days

    itest = 0

    if (iwetsfc == 1) then

       if     (xx4(1,1) <= 0.) then
          itest = 1     ! VcGw
       elseif (iwetveg  == 1   .and. &
               xx4(1,2) >= 0.  .and. &
               xx4(1,2) <= vw) then
          itest = 2     ! VewGw
       elseif (xx4(1,3) > vwm) then       ! was (xx4(1,3) >= 0.) 
          itest = 3     ! VeGw
       else
          itest = 91
       endif

    else ! iwetsfc = 0

       if     (xx4(1,1) <= 0.  .and. &
               xx4(2,1) <= 0.) then
          itest = 4     ! VcGc
       elseif (xx4(1,4) <= 0.  .and. &
               xx4(2,4) >= 0.) then
          itest = 5     ! VcGe
       elseif (xx3(1,1) <= 0.) then
          itest = 6     ! VcGo
       elseif (iwetveg  == 1   .and. &
               xx4(1,2) >= 0.  .and. &
               xx4(1,2) <= vw  .and. &
               xx4(2,2) <= 0.) then
          itest = 7     ! VewGc
       elseif (iwetveg  == 1   .and. &
               xx4(1,5) >= 0.  .and. &
               xx4(1,5) <= vw  .and. &
               xx4(2,5) >= 0.) then
          itest = 8     ! VeGe
       elseif (iwetveg  == 1   .and. & 
               xx3(1,2) >= 0.  .and. &
               xx3(1,2) <= vw) then
             itest = 9  ! VewGo
       elseif (xx4(1,3) >= 0.  .and. &
               xx4(2,3) <= 0.) then
          itest = 10    ! VeGc
       elseif (xx4(1,6) >= 0.  .and. &
               xx4(2,6) >= 0.) then
          itest = 11    ! VeGe
       elseif (xx3(1,3) > vwm) then       ! was (xx4(1,3) >= 0.) 
          itest = 12    ! VeGo
       else
          itest = 92
       endif

    endif  ! iwetsfc

    20 continue

    if (itest > 90) then
       print*, '[itest > 90] ',itest,iwetsfc,iwetveg,iwl

       print*, 'veg_water ',veg_water

       print*, '! xx4(1,1), xx4(2,1) ',xx4(1,1), xx4(2,1)
       print*, '! xx4(1,2), xx4(2,2) ',xx4(1,2), xx4(2,2)
       print*, '! xx4(1,3), xx4(2,3) ',xx4(1,3), xx4(2,3)
       print*, '! xx4(1,4), xx4(2,4) ',xx4(1,4), xx4(2,4)
       print*, '! xx4(1,5), xx4(2,5) ',xx4(1,5), xx4(2,5)
       print*, '! xx4(1,6), xx4(2,6) ',xx4(1,6), xx4(2,6)
       print*, '! xx3(1,1)           ',xx3(1,1)
       print*, '! xx3(1,2)           ',xx3(1,2)
       print*, '! xx3(1,3)           ',xx3(1,3)

       stop 'itest > 90 '
    endif 

! Apply results of successful test

    if     (itest == 1 .or. itest == 4) then
       wxfervc = xx4(1,1)
       wxfersc = xx4(2,1)
       hxfervc = xx4(3,1)
       hxfersc = xx4(4,1)

       transp  = 0.
    elseif (itest == 5) then
       wxfervc = xx4(1,4)
       wxfersc = xx4(2,4)
       hxfervc = xx4(3,4)
       hxfersc = xx4(4,4)
       transp  = 0.
    elseif (itest == 6) then
       wxfervc = xx3(1,1)
       wxfersc = 0.
       hxfervc = xx3(2,1)
       hxfersc = xx3(3,1)
       transp  = 0.
    elseif (itest == 2 .or. itest == 7) then
       wxfervc = xx4(1,2) * a3 / (a3 + a4)
       wxfersc = xx4(2,2)
       hxfervc = xx4(3,2)
       hxfersc = xx4(4,2)
       transp  = xx4(1,2) - wxfervc
    elseif (itest == 8) then
       wxfervc = xx4(1,5) * a3 / (a3 + a4)
       wxfersc = xx4(2,5)
       hxfervc = xx4(3,5)
       hxfersc = xx4(4,5)
       transp  = xx4(1,5) - wxfervc
    elseif (itest == 9) then
       wxfervc = xx3(1,2) * a3 / (a3 + a4)
       wxfersc = 0.
       hxfervc = xx3(2,2)
       hxfersc = xx3(3,2)
       transp  = xx3(1,2) - wxfervc
    elseif (itest == 3 .or. itest == 10) then
       wxfervc = veg_water
       wxfersc = xx4(2,3)
       hxfervc = xx4(3,3)
       hxfersc = xx4(4,3)
       transp  = xx4(1,3) - wxfervc
       if (a4 == 0.0) transp = 0.0
    elseif (itest == 11) then
       wxfervc = veg_water
       wxfersc = xx4(2,6)
       hxfervc = xx4(3,6)
       hxfersc = xx4(4,6)
       transp  = xx4(1,6) - wxfervc
       if (a4 == 0.0) transp = 0.0
    elseif (itest == 12) then
       wxfervc = veg_water
       wxfersc = 0.
       hxfervc = xx3(2,3)
       hxfersc = xx3(3,3)
       transp  = xx3(1,3) - wxfervc
       if (a4 == 0.0) transp = 0.0
    endif

! Update components

    cantemp = cantemp + (hxfervc + hxfersc - cp * sxfer_t) * hcapcani

    veg_temp = veg_temp + (radveg - hxfervc - (wxfervc + transp) * alvl) * hcapvegi

    energy_per_m2 = energy_per_m2 + radsfc - hxfersc - wxfersc * alvl

    canshv = canshv + (wxfervc + transp + wxfersc - sxfer_r) * canairi

    veg_water = max(0.,veg_water - wxfervc)

    specvol = .001
    if (wxfersc > 0. .and. sfcwater_mass > 1.e-3) specvol = sfcwater_depth / sfcwater_mass

    sfcwater_mass = sfcwater_mass - wxfersc

    sfcwater_depth = sfcwater_depth - wxfersc * specvol

 endif  ! End of vegetation section

return
end subroutine canopy

!===============================================================================

subroutine vegndvi(iwl, leaf_class, timefac_ndvi, veg_height,           &
                   veg_ndvip, veg_ndvif, veg_ndvic,                     &
                   veg_tai, veg_lai, veg_fracarea, veg_albedo, veg_rough)

use leaf_coms, only: veg_frac, albv_brown, albv_green, sai, dead_frac,  &
                     fpar_max, veg_clump, glai_max, dfpardsr, fpar_min, &
                     sr_max, sr_min, tai_max
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwl        ! index of current land cell
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

return
end subroutine vegndvi

!===============================================================================

 subroutine sing_print(iwl,ieqn,nsize,aa,yy,glatw,glonw)

 use misc_coms,   only: io6
 use consts_coms, only: r8

 implicit none

 integer, intent(in) :: iwl, ieqn, nsize

 real(r8), intent(in) :: aa(nsize,nsize), yy(nsize)
 real,     intent(in) :: glatw, glonw

 integer :: irow

 print*, 'singular matrix '
 print*, 'iwl, ieqn, nsize ', iwl, ieqn, nsize
 print*, 'glatw, glonw ', glatw, glonw

 do irow = 1,nsize
    write(io6,'(a,i5,5e15.5)') 'row, aa, yy ',irow,aa(irow,1:nsize),yy(irow)
 enddo

 stop 'singular matrix '

 end subroutine sing_print

! Remaining issues:
! 
! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Build lookup tables, especially for things with exponentials?

End module leaf4_canopy
