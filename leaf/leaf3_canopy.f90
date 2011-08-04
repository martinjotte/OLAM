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
Module leaf3_canopy

Contains

subroutine canopy(iwl, nlev_sfcwater, ntext_soil, leaf_class, ktrans,   &
                  soil_water, soil_fracliq, soil_tempk,                 &
                  sfcwater_mass, sfcwater_tempk,                        &
                  veg_height, veg_rough, veg_tai, veg_lai,              &
                  hcapveg, can_depth,                                   &
                  rhos, vels, prss, pcpg, qpcpg,                        &
                  rshort, rshort_v, rlong_v, sxfer_t, sxfer_r, ustar,   &
                  snowfac, vf, surface_ssh, ground_shv,                 &
                  veg_water, veg_temp, can_temp, can_shv,               &
                  wshed, qwshed, transp, stom_resist,                   &
                  hxfergc, wxfergc, hxfersc, wxfersc, hxferca,          &
                  hxfervc, wxfervc, rdi, rb, time8,                     &
                  lsl, ed_transp, ed_patch)

use leaf_coms,   only: nzg, soil_rough, dt_leaf,  &
                       slpots, slmsts, slbs, kroot, rcmin, soilcp, dslz

use consts_coms, only: cp, vonk, eps_vap, alvl, cliq, cice, alli, rvap

use ed_structure_defs
use misc_coms, only: io6
use mem_leaf,  only: itab_wl

implicit none

integer, intent(in)  :: iwl           ! index of current land cell
integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
integer, intent(in)  :: ntext_soil(nzg) ! soil textural class
integer, intent(in)  :: leaf_class      ! leaf class (vegetation class)
integer, intent(in)  :: lsl
integer, intent(out) :: ktrans          ! k index of soil layer supplying transp

real, intent(in) :: soil_water(nzg)   ! soil water content [vol_water/vol_tot]
real, intent(in) :: soil_fracliq(nzg) ! fraction of soil moisture in liquid phase
real, intent(in) :: soil_tempk        ! soil temp [K]
real, intent(in) :: sfcwater_mass     ! surface water mass [kg/m^2]
real, intent(in) :: sfcwater_tempk    ! surface water temp [K]
real, intent(in) :: veg_height  ! veg height [m]
real, intent(in) :: veg_rough   ! veg roughess height [m]
real, intent(in) :: veg_tai     ! veg total area index
real, intent(in) :: veg_lai     ! veg leaf area index
real, intent(in) :: hcapveg     ! veg heat capacity [J/(m^2 K)]
real, intent(in) :: can_depth   ! canopy depth for heat and vap capacity [m]
real, intent(in) :: rhos        ! atmospheric air density [kg/m^3]
real, intent(in) :: vels        ! surface wind speed [m/s]
real, intent(in) :: prss        ! air pressure [Pa]
real, intent(in) :: pcpg        ! new precip amount this leaf timestep [kg/m^2]
real, intent(in) :: qpcpg       ! new precip energy this leaf timestep [J/m^2]
real, intent(in) :: rshort      ! downward sfc s/w rad flux [W/m^2]
real, intent(in) :: rshort_v    ! s/w rad flux absorbed by veg [W/m^2]
real, intent(in) :: rlong_v     ! l/w rad flux absorbed by veg [W/m^2]
real, intent(in) :: sxfer_t     ! surface heat xfer this step [kg_air K/m^2]
real, intent(in) :: sxfer_r     ! surface vapor xfer this step [kg_vap/m^2]
real, intent(in) :: ustar       ! friction velocity [m/s]
real, intent(in) :: snowfac     ! fractional veg burial by snowcover
real, intent(in) :: vf          ! fractional coverage of non-buried part of veg
real, intent(in) :: surface_ssh ! surface sat spec hum [kg_vap/kg_air]
real, intent(in) :: ground_shv  ! soil vapor spec hum [kg_vap/kg_air]

real, intent(inout) :: veg_water   ! veg sfc water content [kg/m^2]
real, intent(inout) :: veg_temp    ! veg temp [K]
real, intent(inout) :: can_temp    ! canopy air temp [K]
real, intent(inout) :: can_shv     ! canopy air vapor spec hum [kg_vap/kg_air]
real, intent(inout) :: stom_resist ! veg stomatal resistance [s/m]

real, intent(out) :: hxfergc ! soil-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfergc ! soil-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: wshed   ! water shed from veg this LEAF timestep [kg/m^2]
real, intent(out) :: qwshed  ! water energy shed from veg this LEAF timestep [J/m^2]
real, intent(out) :: transp  ! transpiration xfer this LEAF timestep [kg_vap/m^2]
real, intent(out) :: hxferca ! can_air-to-atm heat xfer this step [J/m^2]
real, intent(out) :: hxfervc ! veg-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfervc ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
real, intent(out) :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
real, intent(out) :: rb      ! veg-to-can_air resistance [s/m]
real(kind=8), intent(in) :: time8   ! model time [s]

real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]
type(patch), target, optional :: ed_patch

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

! Local variables

integer :: k        ! loop index over soil layers
integer :: nts      ! soil textural class for current soil layer

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
real :: slpotv      ! soil water potential [m]
real :: swp         ! soil water potential factor for stomatal control
real :: tvegc       ! vegetation temperature (deg C)
real :: tvegk       ! vegetation temperature (K)
real :: wtroot      ! not used
real :: zognd       ! soil roughness height [m]
real :: zoveg       ! vegetation roughness height [m]
real :: zdisp       ! vegetation displacement height remaining after snowcover [m]
real :: zveg        ! vegetation height remaining after snowcover [m]
real :: wxfer       ! (saturated soil or sfc water)-to-can_air vap xfer this step [kg_vap/m^2]
real :: transp_test ! test value of transpiration flux [kg_vap/(m^2 s)]
real :: rc          ! stomatal resistance [s/m]

real :: canair
real :: canhcap
real :: f1
real :: f2
real :: f3
real :: a1
real :: a2
real :: a3
real :: a4
real :: b1
real :: b2
real :: b3
real :: b4
real :: b5
real :: b6
real :: b7
real :: b8
real :: b9
real :: evapotransp

real, external :: rhovsil ! function to compute sat vapor density
real, external :: eslf    ! function to compute sat vapor pressure

! Can_air-to-atm heat xfer in [J/m^2]

hxferca = cp * sxfer_t

! Canopy air (vapor) and heat capacity

canair = rhos * can_depth
canhcap = cp * canair

! Initialize wshed = qwshed = 0

wshed  = 0.
qwshed = 0.

!print*, ' '
!write(io6,*) 'can005 ',iwl,can_shv,can_temp,canair
!write(io6,*) 'can006 ',rhos,can_depth

! Check whether this land cell has exposed vegetation

if ((.not.present(ed_patch)) .and. (vf < .001)) then

! If the TAI is very small or if snow mostly covers the vegetation, 
! then BYPASS THE VEGETATION COMPUTATIONS.

! Aerodynamic conductance for bare soil or snow based on Garratt.

   rdi = .2 * ustar

! Set transpiration to zero

   transp = 0.
   ktrans = 0

! Check if any surface water layers currently exist

   if (nlev_sfcwater >= 1) then

! If surface water is present, compute heat and vapor xfer between 
! surface water and canopy air.

      hxfergc = 0.
      hxfersc = cp * rhos * dt_leaf * rdi * (sfcwater_tempk - can_temp)

      wxfergc = 0.
      wxfersc = min(sfcwater_mass,                                 &
                    rhos * dt_leaf * rdi * (surface_ssh - can_shv) )

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw01 ',wxfersc,sfcwater_mass,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw02 ',surface_ssh, can_shv
!endif

! Update canopy temperature and vapor specific humidity, and set vegetation
! temperature to canopy air value

!if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!     iwl == 821 .or. iwl == 822) then
!   write(io6,*) ' '
!   write(io6,*) 'can0 ',iwl,can_temp,dble(can_temp)
!    write(io6,*) 'can01 ',iwl,can_shv,wxfersc,sxfer_r,canair
!endif

      can_temp = can_temp + (hxfersc - hxferca) / canhcap
      can_shv  = can_shv  + (wxfersc - sxfer_r) / canair
      veg_temp = can_temp

!if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!     iwl == 821 .or. iwl == 822) then
!   write(io6,*) ' '
!   write(io6,*) 'can1 ',can_temp,dble(can_temp)
!   write(io6,*) 'can2 ',hxfersc,dble(hxfersc)
!   write(io6,*) 'can3 ',hxferca,dble(hxferca)
!   write(io6,*) 'can4 ',canhcap,dble(canhcap)
!   write(io6,*) 'can5 ',rdi,dble(rdi)
!   write(io6,*) 'can6 ',dt_leaf,dble(dt_leaf)
!   write(io6,*) 'can7 ',rhos,dble(rhos)
!   write(io6,*) 'can8 ',cp,dble(cp)
!   write(io6,*) 'can9 ',sfcwater_tempk,dble(sfcwater_tempk)
!   write(io6,*) ' '
!endif

      return

   else

! If surface water is absent, compute heat xfer between soil and canopy air

      hxfergc = cp * rhos * dt_leaf * rdi * (soil_tempk - can_temp)
      hxfersc = 0.

! Compare saturation vapor specific humidity at soil temperature against 
! canopy air vapor specific humidity

      if (surface_ssh < can_shv) then

! If saturation vapor specific humidity at soil temperature is less than 
! canopy air vapor specific humidity, compute dew formation contribution 
! to surface water

         wxfergc = 0.
         wxfersc = rhos * dt_leaf * rdi * (surface_ssh - can_shv)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw03 ',wxfersc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw04 ',surface_ssh, can_shv
!endif

      elseif (ground_shv > can_shv) then

! Else, if equilibrium soil vapor specific humidity is greater than 
! canopy air vapor specific humidity, compute evaporation from soil

!print*, 'can099 ',iwl,rhos,dt_leaf,rdi,ground_shv,can_shv

         wxfergc = max(0.,rhos * dt_leaf * rdi * (ground_shv - can_shv))
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw05 ',wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw06 ',ground_shv, can_shv
!endif

      else

! If neither of the above is true, both vapor fluxes are zero.

         wxfergc = 0.
         wxfersc = 0.

      endif

! Update canopy temperature and vapor specific humidity, and set vegetation
! temperature to canopy air value


! if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!      iwl == 821 .or. iwl == 822) then
!    write(io6,*) ' '
!    write(io6,*) 'can10 ',iwl,can_temp,dble(can_temp)
!    write(io6,*) 'can101 ',hxfergc,hxfersc,hxferca,canhcap
! endif

      can_temp = can_temp + (hxfergc + hxfersc - hxferca) / canhcap
      can_shv  = can_shv  + (wxfergc + wxfersc - sxfer_r) / canair
      veg_temp = can_temp


! if (itab_wl(iwl)%iwglobe == 821 .or. itab_wl(iwl)%iwglobe == 822 .or.  &
!      iwl == 821 .or. iwl == 822) then
!    write(io6,*) ' '
!    write(io6,*) 'can11 ',can_temp,dble(can_temp)
!    write(io6,*) 'can12 ',hxfergc,dble(hxfergc)
!    write(io6,*) 'can13 ',hxfersc,dble(hxfersc)
!    write(io6,*) 'can14 ',hxferca,dble(hxferca)
!    write(io6,*) 'can15 ',canhcap,dble(canhcap)
!    write(io6,*) 'can16 ',rdi,dble(rdi)
!    write(io6,*) 'can17 ',dt_leaf,dble(dt_leaf)
!    write(io6,*) 'can18 ',rhos,dble(rhos)
!    write(io6,*) 'can19 ',cp,dble(cp)
!    write(io6,*) ' '
! endif

      return

   endif

else

! If vegetation is sufficiently abundant and not covered by snow,
! COMPUTE CANOPY XFER WITH VEGETATION INFLUENCE

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
   
! Check if any surface water layers currently exist

   if (nlev_sfcwater >= 1) then

! If surface water is present, compute heat and vapor xfer between 
! surface water and canopy air.

      hxfergc = 0.
      hxfersc = cp * rhos * dt_leaf * rdi * (sfcwater_tempk - can_temp)

      wxfergc = 0.
      wxfersc = min(sfcwater_mass,                                 &
                    rhos * dt_leaf * rdi * (surface_ssh - can_shv) )

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw07 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw08 ',surface_ssh, can_shv
!endif

   else

! If surface water is absent, compute heat xfer between soil and canopy air

      hxfergc = cp * rhos * dt_leaf * rdi * (soil_tempk - can_temp)
      hxfersc = 0.

! Compare saturation vapor specific humidity at soil temperature against 
! canopy air vapor specific humidity

      if (surface_ssh < can_shv) then

! If saturation vapor specific humidity at soil temperature is less than 
! canopy air vapor specific humidity, compute dew formation contribution 
! to surface water

         wxfergc = 0.
         wxfersc = rhos * dt_leaf * rdi * (surface_ssh - can_shv)

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw09 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw010 ',surface_ssh, can_shv
!endif

      elseif (ground_shv > can_shv) then

! Else, if equilibrium soil vapor specific humidity is greater than 
! canopy air vapor specific humidity, compute evaporation from soil

         wxfergc = max(0.,rhos * dt_leaf * rdi * (ground_shv - can_shv))
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw011 ',wxfersc,wxfergc,rhos,dt_leaf,rdi
!   write(io6,*) 'sfcw012 ',ground_shv, can_shv
!endif

      else

! If neither of the above is true, both vapor fluxes are zero.

         wxfergc = 0.
         wxfersc = 0.

         
!if (iwl == 102) then
!   write(io6,*) 'sfcw013 ',wxfersc,wxfergc
!endif

      endif

   endif

! Here, LEAF3 and ED2 branch off.

   if(.not.present(ed_patch))then

! TAI and LAI reduced by ground snowcover (for low vegetation)

   stai = veg_tai * (1. - snowfac)
   slai = veg_lai * (1. - snowfac)

! Evaluate vegetation Celsius temperature

   tvegc = veg_temp - 273.15

! Check if precipitation is occurring

   if (pcpg > 1.e-12) then
   
! If precipitation, add intercepted mass and energy to vegetation surface
!   (Ignore energy of water already on vegetation surface)

      veg_water = veg_water + pcpg * vf
      qwtot = hcapveg * tvegc + qpcpg * vf

! Compute equilbrium temperature of veg + precipitation

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
      endif
      
   endif

! Vegetation saturation vapor density and gradient

   veg_rhovs = rhovsil(tvegc)
   veg_rhovsp = rhovsil(tvegc+1.) - veg_rhovs

! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
! WITHOUT the LAI factor; the LAI factor is included later in computing
! xfers that involve rb and/or rc.

   rb = (1. + .5 * veg_tai) / (.01 * sqrt(ustar * c1))

! Initialize soil water potential to -600 m (later to be converted to Pascals)
! prior to search for wettest soil layer in root zone.  This value is low 
! enough to effectively prevent transpiration.

   swp = -600.

! Initialize ktrans to zero prior to search for wettest soil layer in root zone
   
   ktrans = 0

! Loop over soil layers in the root zone

   do k = kroot(leaf_class),nzg
      nts = ntext_soil(k)

! Soil water potential

      slpotv = slpots(nts) * (slmsts(nts) / soil_water(k)) ** slbs(nts)

! Multiply by liquid fraction (ice is unavailable for transpiration)

      slpotv = slpotv * soil_fracliq(k)

! Find layer in root zone with highest slpotv AND soil_water above minimum soilcp
! Set ktrans to this layer

      if (slpotv > swp .and. soil_water(k) > soilcp(nts)) then
         swp = slpotv
         ktrans = k
      endif
   enddo

  swp = swp * 9810. ! convert from meters to Pascals (hydrostatic eqn for water)

! Bypass stomatal resistance computation if root zone is too dry
   
   if (ktrans < 1) then

      rc = 1.e18

   else

! Compute saturation vapor pressure at veg temp

      esat_veg  = eslf(tvegc)

! Compute vapor pressure in canopy air from equation of state

      e_can = can_shv * rhos * rvap * can_temp

! Compute vapor pressure at leaf surface using rc from previous timestep

      rc = stom_resist
      e_leaf = (rb * esat_veg + rc * e_can) / (rb + rc)
      vpd = max(0.,esat_veg - e_leaf)

! Evaluate 5 environmental factors and new rc

      ftlo = 1. + exp(-stlo * (veg_temp - btlo))
      fthi = 1. + exp(-sthi * (veg_temp - bthi))
      frad = 1. + exp(-srad * (rshort   - brad))
      fswp = 1. + exp(-sswp * (swp      - bswp))
      fvpd = 1. + exp(-svpd * (vpd      - bvpd))

! Compute asymptotoc value of stomatal resistance based on environmental factors

      rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(leaf_class)

! Update stomatal conductance assuming 15-minute response time 
! (requires dt_leaf <= 900.) 

      rc = 1. / (1. / rc + .0011 * dt_leaf * (1. / rc_inf - 1. / rc))

! Limit maximum transpiration to be <= 500 W/m^2 by increasing rc if necessary.

      transp_test = alvl * veg_lai * (veg_rhovs - rhos * can_shv) / (rb + rc)
      if (transp_test > 500.) then
         rc = (rb + rc) * transp_test * .002 - rb
      endif      
      
   endif
      
   stom_resist = rc

!-------------------------------------------------------------------------------
! Begin implicit exchange of heat and moisture between vegetation and canopy air
!-------------------------------------------------------------------------------

! veg_water fractional coverage
   
   sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66667)

! auxiliary quantities

   f1 = (hxfergc + hxfersc - hxferca) / canhcap
   f2 = (wxfergc + wxfersc - sxfer_r) / canair
   f3 = dt_leaf * (rshort_v + rlong_v) / hcapveg

   a1 = dt_leaf * 2. * stai / rb                  ! vap xfer coef (dew formation)
   a2 = cp * rhos * a1                            ! heat xfer coef
   a3 = a1 * sigmaw                               ! vap xfer coef (veg_water evap)
   a4 = dt_leaf * slai * (1.-sigmaw) / (rb + rc)  ! vap xfer coef (transp)

   b1 = canhcap + a2
   b2 = canhcap * a2
   b3 = hcapveg * b1 + b2
   b4 = b2 * (can_temp + f1)
   b5 = b1 * hcapveg * (veg_temp + f3)
   b6 = veg_rhovs + veg_rhovsp * ((b5 + b4) / b3 - veg_temp)
   b7 = veg_rhovsp * b1 * alvl / b3
   b8 = b7 + 1. / can_depth
   b9 = b6 - rhos * (can_shv + f2)

! First, assume case of condensation/deposition of vapor onto veg

   evapotransp = a1 * b9 / (1. + a1 * b8)
   
! Now test this assumption

   if (evapotransp <= 0.) then

! Assumption was correct.  Compute individual transpiration and evaporation

      wxfervc = evapotransp
      transp = 0.

   elseif (veg_water > 1.e-12) then

! Assumption was incorrect.  Now, consider the case where veg_water is present.

! Assume that veg_water does not all evaporate

      evapotransp = (a3 + a4) * b9 / (1. + (a3 + a4) * b8)
      wxfervc = evapotransp * a3 / (a3 + a4)

! Now test this assumption

      if (wxfervc <= veg_water) then

! Assumption was correct.  Compute transpiration

         transp = evapotransp - wxfervc

      else

! Assumption was incorrect.  All of veg_water evaporates this leaf step
!  Compute combined and individual transpiration and evaporation

         evapotransp = (veg_water + a4 * b9) / (1. + a4 * b8)
         wxfervc = veg_water
         transp = evapotransp - wxfervc
      
      endif

      if(a4 == 0.0)transp = 0.0

   else

! Original condensation/deposition assumption was incorrect, and there is 
!   no veg_water.

      veg_water = 0.
      evapotransp = a4 * b9 / (1. + a4 * b8)
      wxfervc = 0.
      transp = evapotransp

   endif

! Compute remaining quantities

   can_shv = can_shv + f2 + evapotransp / canair
   veg_temp = (b5 + b4 - b1 * alvl * evapotransp) / b3
   hxfervc = (b2 * veg_temp - b4) / b1
   can_temp = can_temp + f1 + hxfervc / canhcap
   veg_water = veg_water - wxfervc

   else

      ktrans = 0 
      transp = 0.
      ed_transp(:) = 0.

      call ed_canopy_update(ed_patch, vels, rhos, prss, pcpg, qpcpg,   &
           wshed, qwshed, canair, canhcap, dt_leaf, hxfergc, hxferca,  &
           wxfergc, hxfersc, wxfersc, sxfer_r, ed_transp, ntext_soil,  &
           soil_water, soil_fracliq, lsl)

   endif

endif

return
end subroutine canopy

!===============================================================================

subroutine vegndvi(iwl, leaf_class, timefac_ndvi, veg_height,         &
                   veg_ndvip, veg_ndvif, veg_ndvic,                     &
                   veg_tai, veg_lai, veg_fracarea, veg_albedo, veg_rough)

use leaf_coms, only: veg_frac, albv_brown, albv_green, sai, dead_frac,   &
                     fpar_max, veg_clump, glai_max, dfpardsr, fpar_min,  &
                     sr_max, sr_min, tai_max
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwl      ! index of current land cell
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

! Remaining issues:
! 
! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Build lookup tables, especially for things with exponentials?

End module leaf3_canopy
