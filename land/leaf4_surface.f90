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
Module leaf4_surface

Contains

  subroutine sfcwater(iland, iwsfc, icomb, wfree1, qwfree1, dwfree1,       &
                      head1, nlev_sfcwater,                                &
                      sfcwater_mass, sfcwater_energy, sfcwater_depth,      &
                      sfcwater_tempk, sfcwater_fracliq, energy_per_m2,     &
                      soil_water, soil_energy, specifheat_drysoil,         &
                      soil_watfrac, soil_tempk, soil_fracliq , thermcond_soil )

  use leaf_coms,   only: nzs, dt_leaf, snowmin_expl
  use mem_land,    only: nzg, dslz, dslzi, dslzo2
  use consts_coms, only: alvi, cice, cliq, alli
  use misc_coms,   only: io6
  use leaf4_plot,  only: leaf_plot
  use therm_lib,   only: qwtk

  implicit none

  integer, intent(in)    :: iland         ! current land cell index
  integer, intent(in)    :: iwsfc         ! current sfcg cell index
  integer, intent(in)    :: icomb         ! implicit heat balance flag [0=no, 1=yes]
  real,    intent(out)   :: wfree1        ! free liquid in lowest sfcwater layer [kg/m^2]
  real,    intent(out)   :: qwfree1       ! energy carried by wfree1 [J/m^2]
  real,    intent(out)   :: dwfree1       ! depth carried by wfree1 [m]
  real,    intent(inout) :: head1         ! Top boundary hydraulic head for soil model [m]
  integer, intent(inout) :: nlev_sfcwater ! # of active sfc water levels

  real,    intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
  real,    intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
  real,    intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
  real,    intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
  real,    intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
  real,    intent(inout) :: energy_per_m2   (nzs) ! sfcwater energy [J/m^2]

  real,    intent(inout) :: soil_water        (nzg) ! soil water [water_vol/total_vol]
  real,    intent(inout) :: soil_energy       (nzg) ! soil internal energy [J/m^3]
  real,    intent(in)    :: specifheat_drysoil(nzg) ! specific heat of dry soil [J/(m^3 K)]
  real,    intent(in)    :: soil_tempk        (nzg) ! soil temperature [K]
  real,    intent(in)    :: soil_fracliq      (nzg) ! fraction of soil water in liq phase
  real,    intent(in)    :: soil_watfrac      (nzg) ! water fraction in soil layer [vol_water/vol_total]
  real,    intent(in)    :: thermcond_soil    (nzg) ! soil thermal conductivity [W/(K m)]

  ! Local variables

  real :: hxfers        (nzs+1) ! sfcwater heat xfer [J/m2] 
  real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]

  integer :: k         ! vertical index over sfcwater layers
  integer :: kold      ! vertical index of adjacent lower sfcwater layer

  real :: soil_rfactor ! soil thermal resistance [K m^2/W]
  real :: hxfergs   ! energy transfer from soil to sfcwater this step [J/m^2]
  real :: snden     ! sfcwater density [kg/m^3]
  real :: vegfracc  ! 1 minus veg fractional area
  real :: wfree     ! free liquid in sfcwater layer that can percolate out [kg/m^2]
  real :: qwfree    ! energy carried by wfree [J/m^2]
  real :: dwfree    ! depth carried by wfree [m]
  real :: fracstep  ! ratio of leaf timestep to snow density exponential decay time
  real :: totsnow   ! sum of mass over sfcwater layers [kg/m^2]
  real :: hcapsoil  ! soil(nzg) heat capacity [J/(m^2 K)]
  real :: sndenmin  ! minimum sfcwater density [kg/m^3]
  real :: tempk     ! Kelvin temperature of snowcover, limited to max of 273.15 [K]
  real :: wcap_min  ! minimum surface water water [kg/m^2]

  ! Local parameters

  integer, parameter :: iland_print = 0

  ! Initialize free water values to zero

  wfree1   = 0.
  qwfree1  = 0.
  dwfree1  = 0.

  head1    = 0.

  wcap_min = dt_leaf * 1.e-9 ! 1.e-9 is 1 mm pcp/dew in 31 years    [kg/m^2]

  ! If there are multiple water layers and the top layer mass is below
  ! threshold, iteratively transfer quantities to the layer below until
  ! the threshold is reached

  if (nlev_sfcwater > 1) then
     do while (nlev_sfcwater > 1 .and. sfcwater_mass(nlev_sfcwater) < 10.) ! 10 kg/m^2 threshold 
        k = nlev_sfcwater

        sfcwater_mass(k-1)  = sfcwater_mass(k-1)  + sfcwater_mass(k)
        energy_per_m2(k-1)  = energy_per_m2(k-1)  + energy_per_m2(k)
        sfcwater_depth(k-1) = sfcwater_depth(k-1) + sfcwater_depth(k)

        sfcwater_mass   (k) = 0.
        sfcwater_energy (k) = 0.
        sfcwater_depth  (k) = 0.
        energy_per_m2   (k) = 0.
        sfcwater_fracliq(k) = 0.

        nlev_sfcwater = nlev_sfcwater - 1
     enddo
  endif

  ! Take inventory of sfcwater_mass(1) now that exchange with canopy is
  ! complete.  If sfcwater_mass(1) is less than threshold, set it to zero and return.

  if (sfcwater_mass(1) < wcap_min) then

 ! if (iland == 10110) print*, 'swat15.1 ',soil_water(nzg),sfcwater_mass(1),dslzi(nzg) * sfcwater_mass(1) * .001

     sfcwater_mass  (1) = 0.
     sfcwater_energy(1) = 0.
     sfcwater_depth (1) = 0.
     energy_per_m2  (1) = 0.

     sfcwater_tempk  (1) = soil_tempk(nzg)
     sfcwater_fracliq(1) = 0.

     nlev_sfcwater = 0

 ! if (iland == 10110) print*, 'swat15.2 ',soil_water(nzg),sfcwater_mass(1),dslzi(nzg) * sfcwater_mass(1) * .001

     return
  endif

  ! If surface water was present at the beginning of this leaf step
  ! AND implicit heat balance between sfcwater(1) and soil(nzg) was NOT done,
  ! compute heat transfer between snowcover layers and between soil and snowcover

  if (nlev_sfcwater > 0 .and. icomb == 0) then

     ! Surface water was present.

     ! Loop over existing sfcwater layers

     do k = 1,nlev_sfcwater

        ! Compute snow heat resistance times HALF layer depth (sfcwater_rfactor).
        ! Sfcwater_tempk(k) should be correctly balanced value at this point, so 
        ! sfcwater_rfactor(k) should have correct value.  Formula applies to snow,
        ! so limit temperature to no greater than 273.15.

        snden = sfcwater_mass(k) / sfcwater_depth(k)
        tempk = min(273.15,sfcwater_tempk(k))

        sfcwater_rfactor(k) = .5 * sfcwater_depth(k) &
           / (1.093e-3 * exp(.028 * tempk) *         &
           (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

     enddo

     ! Internal heat transfer if multiple sfcwater layers are present

     if (nlev_sfcwater >= 2) then

        ! Zero out sfcwater internal heat transfer array at bottom and top surfaces.
        ! Energy transfer at bottom and top are applied separately.

        hxfers(1) = 0.
        hxfers(nlev_sfcwater+1) = 0. 

        ! Compute internal sfcwater energy xfer if at least two layers exist [J/m2]

        do k = 2,nlev_sfcwater
           hxfers(k) = dt_leaf * (sfcwater_tempk(k-1) - sfcwater_tempk(k)) &
                     / (sfcwater_rfactor(k-1) + sfcwater_rfactor(k))      
        enddo

        ! Add contributions to sfcwater energy_per_m2 from internal transfers of 
        ! sensible heat (and latent heat, implicitly included in sfcwater_rfactor).
        ! This excludes energy transfer from internal gravitational draining of 
        ! free water mass.

        do k = 1,nlev_sfcwater
           energy_per_m2(k) = energy_per_m2(k) + hxfers(k) - hxfers(k+1)
        enddo

     endif

     ! Compute heat resistance of top HALF of top soil layer (soil_rfactor).

     soil_rfactor = dslzo2(nzg) / thermcond_soil(nzg)

     ! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
     ! layer, and apply to both.

     hxfergs = dt_leaf * (soil_tempk(nzg) - sfcwater_tempk(1))   &
             / (soil_rfactor + sfcwater_rfactor(1))

     energy_per_m2(1) = energy_per_m2(1) + hxfergs

     soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)

  endif

  ! If we are here sufficient sfcwater mass remains,
  ! and if nlev_sfcwater is less than 1, set it to 1.

  if (nlev_sfcwater < 1) then
     nlev_sfcwater = 1
  endif

  ! Now that new mass and energy contributions have been applied, loop downward
  ! through all existing sfcwater layers beginning with top layer

  do k = nlev_sfcwater,1,-1

     ! Diagnose sfcwater temperature and liquid fraction.  Use qwtk instead
     ! of qtk because energy_per_m2 is used instead of sfcwater_energy.
     ! ("dryhcap" = 100 is very small value)

     call qwtk(energy_per_m2(k),sfcwater_mass(k),100., &
               sfcwater_tempk(k),sfcwater_fracliq(k))

     ! If k = 1 (lowest sfcwater layer) and sfcwater mass is below limiting value,
     ! or sfcwater is mostly liquid, assume that it is in thermal equilibrium with
     ! the top soil layer.  Diagnose the sfcwater and soil temperatures and liquid
     ! water fractions together.

     if (k == 1 .and. &
        (sfcwater_mass(1) < snowmin_expl .or. sfcwater_fracliq(1) > .5)) then

        call sfcwater_soil_comb(iland, iwsfc,                      &
                                soil_water(nzg),soil_energy(nzg),  &
                                specifheat_drysoil(nzg),           &
                                sfcwater_mass(1),energy_per_m2(1), &
                                sfcwater_tempk(1),sfcwater_fracliq(1))
     endif

     ! Diagnose sfcwater density.  Make sure sfcwater_depth is not too near zero.

     sfcwater_depth(k) = max(sfcwater_depth(k), .001 * sfcwater_mass(k))

     snden = sfcwater_mass(k) / sfcwater_depth(k)

     ! Assume that as snow ages on ground, its density difference with a limiting
     ! maximum value (currently set to 400 kg/m^3) decays exponentially (with a
     ! decay time currently set to about 3 weeks).  If sfcwater density is less
     ! than this limiting value, apply the density increase for this timestep.

     ! This formulation and decay constants are very crude approximations to a few
     ! widely variable observations of snowpack density and are intended only as a 
     ! rough representation of the tendency for snowcover to compress with time.  
     ! A better formulation that accounts for several environmental factors would
     ! be desirable here.

     if (snden < 400.) then
        fracstep = .5e-6 * dt_leaf  ! .5e-6 is inverse decay time scale
        snden = snden * (1. - fracstep) + 400. * fracstep
        sfcwater_depth(k) = sfcwater_mass(k) / snden
     endif

     ! If liquid exists in current sfcwater layer, any low-density ice structure
     ! tends to collapse.  Increase density accordingly using simple linear relation.

     if (snden < 1.e3 * sfcwater_fracliq(k)) then
        snden = 1.e3 * sfcwater_fracliq(k)
        sfcwater_depth(k) = sfcwater_mass(k) / snden
     endif

     ! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

     if (sfcwater_fracliq(k) > .10) then
        wfree = sfcwater_mass(k) * (sfcwater_fracliq(k) - .10) / .90

        ! Evaluate energy and depth contained in wfree (which is in liquid phase)

        qwfree = wfree * (cliq * (sfcwater_tempk(k) - 273.15) + alli)
        dwfree = wfree * .001

        ! If this is NOT lowest sfcwater layer, drain mass, energy, and depth from this layer
        ! and add to layer below.

        if (k >= 2) then

           ! Check if essentially all of sfcwater_mass(k) will drain from layer

           if (wfree > .999 * sfcwater_mass(k)) then

              ! All sfcwater_mass(k) drains from layer.  Set layer quantities to zero to
              ! avoid truncation error.

              sfcwater_mass(k)   = 0.
              sfcwater_energy(k) = 0.
              energy_per_m2(k)   = 0.
              sfcwater_depth(k)  = 0.

           else

              ! Not all sfcwater_mass(k) drains from layer.  Drain mass, energy, and depth 
              ! of free water out of current layer

              sfcwater_mass(k)  = sfcwater_mass(k)  - wfree
              energy_per_m2(k)  = energy_per_m2(k)  - qwfree
              sfcwater_depth(k) = sfcwater_depth(k) - dwfree

           endif

           sfcwater_mass(k-1)  = sfcwater_mass(k-1)  + wfree
           energy_per_m2(k-1)  = energy_per_m2(k-1)  + qwfree
           sfcwater_depth(k-1) = sfcwater_depth(k-1) + dwfree

        else

           ! If this IS lowest sfcwater layer, copy free water to output arrays

           wfree1  = wfree
           qwfree1 = qwfree
           dwfree1 = dwfree

        endif

     endif

  enddo

  ! 11 Feb 2015: Define head1 from total sfcwater_mass (but leaf4_soil will
  ! still allow no more than wfree1 to actually enter soil on current timestep

  ! July 2017: With LEAF as a 3D groundwater model, limit on head1 has been removed.

  head1 = .001 * sum(sfcwater_mass(1:nlev_sfcwater))

  if (iland == iland_print) then

     call leaf_plot(iland,                             &
                    nlev_sfcwater,                     &
                    linit           = 1,               &
                    lframe          = 1,               &
                    sfcwater_mass   = sfcwater_mass,   & 
                    sfcwater_energy = sfcwater_energy, & 
                    energy_per_m2   = energy_per_m2,   & 
                    sfcwater_depth  = sfcwater_depth   ) 

  endif

  end subroutine sfcwater

!===============================================================================

  subroutine sfcwater_soil_comb(iland, iwsfc, soil_water, soil_energy,       &
                                specifheat_drysoil, sfcwater_mass, energy_per_m2, &
                                sfcwater_tempk, sfcwater_fracliq                  )

  use mem_land,    only: nzg, dslz, dslzi
  use consts_coms, only: cice, cliq, alli
  use misc_coms,   only: io6
  use therm_lib,   only: qwtk

  implicit none

  integer, intent(in) :: iland              ! current land cell index
  integer, intent(in) :: iwsfc              ! current sfcg cell index

  real, intent(in)    :: soil_water         ! soil water [water_vol/total_vol]
  real, intent(inout) :: soil_energy        ! soil internal energy [J/m^3]
  real, intent(in)    :: specifheat_drysoil ! specific heat of dry soil [J/(m^3 K)]
  real, intent(in)    :: sfcwater_mass      ! surface water mass [kg/m^2]
  real, intent(inout) :: energy_per_m2      ! sfcwater energy [J/m^2]
  real, intent(inout) :: sfcwater_tempk     ! surface water temp [K]
  real, intent(inout) :: sfcwater_fracliq   ! fraction of sfc water in liq phase

  ! Local variables

  real :: w_comb       ! (sfcwater + soil) water mass [kg/m^2]
  real :: qw_comb      ! (sfcwater + soil) energy [J/m^2]
  real :: tempk_comb   ! (sfcwater + soil) Kelvin temp [K]
  real :: fracliq_comb ! (sfcwater + soil) frac of water in liq phase
  real :: hcapsoil     ! soil heat capacity [J/(m^2 K)]
  real :: flmin        ! lower bound on sfcwater_fracliq in balance with soil
  real :: flmax        ! upper bound on sfcwater_fracliq in balance with soil

  ! Combined sfcwater and soil water mass per square meter

  w_comb = sfcwater_mass + soil_water * 1.e3 * dslz(nzg)

  ! Combined sfcwater and soil energy per square meter

  qw_comb = energy_per_m2 + soil_energy * dslz(nzg)

  ! Soil heat capacity per square meter

  hcapsoil = specifheat_drysoil * dslz(nzg)

  ! Diagnose equilibrium temperature and fractional liquid/ice water phases

  call qwtk(qw_comb,w_comb,hcapsoil,tempk_comb,fracliq_comb)

  ! Diagnose new energy value for sfcwater based on qw_comb value.

  if (qw_comb < 0.) then

     ! Case of equilibrium temperature below 0 deg C.  Sfcwater_fracliq = 0.

     sfcwater_fracliq = 0.
     sfcwater_tempk   = tempk_comb
     energy_per_m2    = sfcwater_mass * cice * (tempk_comb - 273.15)

  elseif (qw_comb > w_comb * alli) then

     ! Case of equilibrium temperature above 0 deg C.  Sfcwater_fracliq = 1.

     sfcwater_fracliq = 1.
     sfcwater_tempk   = tempk_comb
     energy_per_m2    = sfcwater_mass * (cliq * (tempk_comb - 273.15) + alli)

  else

     ! Equilibrium temperature is 0 deg C.  If sfcwater_mass is near zero,
     ! assume identical values for sfcwater_fracliq and soil_fracliq.

     if (sfcwater_mass < 1.e-6) then

        sfcwater_fracliq = fracliq_comb
        sfcwater_tempk   = 273.15
        energy_per_m2    = sfcwater_mass * sfcwater_fracliq * alli

     else

        ! If sfcwater_mass is larger, determine separate values for
        ! sfcwater_fracliq and soil_fracliq using constraint that the sum
        ! of (mass * fracliq) over both components is (w_comb * fracliq_comb).

        ! Lower bound on sfcwater_fracliq: case with soil_water all liquid

        flmin = (fracliq_comb * w_comb - soil_water * 1.e3 * dslz(nzg)) &
              / sfcwater_mass

        ! Upper bound on sfcwater_fracliq: case with soil_water all ice

        flmax = fracliq_comb * w_comb / sfcwater_mass         

        ! New sfcwater_fracliq value becomes closest value within bounds to old value

        sfcwater_fracliq = max(0.,flmin,min(1.0,flmax,sfcwater_fracliq))
        sfcwater_tempk   = 273.15
        energy_per_m2    = sfcwater_mass * sfcwater_fracliq * alli

     endif

  endif

  ! New energy value for soil is combined energy minus new sfcwater energy

  soil_energy = (qw_comb - energy_per_m2) * dslzi(nzg)

  end subroutine sfcwater_soil_comb

!===============================================================================

  subroutine sfcwater_adjust(iland, iwsfc, glatw, glonw, nlev_sfcwater, sfcwater_mass, &
                             sfcwater_energy, sfcwater_depth, energy_per_m2)

  use leaf_coms, only: nzs, dt_leaf

  use consts_coms, only: alli
  use misc_coms, only: io6

  implicit none

  integer, intent(in)    :: iland         ! current land cell index
  integer, intent(in)    :: iwsfc         ! current sfcg cell index
  real,    intent(in)    :: glatw         ! Latitude of land cell 'center' [deg]
  real,    intent(in)    :: glonw         ! Longitude of land cell 'center' [deg]
  integer, intent(inout) :: nlev_sfcwater ! # of active sfc water levels

  real, intent(inout) :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
  real, intent(inout) :: sfcwater_energy(nzs) ! surface water energy [J/kg]
  real, intent(inout) :: sfcwater_depth (nzs) ! surface water depth [m]
  real, intent(inout) :: energy_per_m2  (nzs) ! sfcwater energy [J/m^2]

  ! Local variables

  real :: mass_new   (nzs) ! mass of new set of sfcwater layers [kg/m^2]
  real :: energy_new (nzs) ! energy of new set of sfcwater layers [J/m^2]
  real :: depth_new  (nzs) ! depth of new set of sfcwater layers [m]

  integer :: k         ! vertical index over sfcwater layers
  integer :: kold      ! vertical index of adjacent lower sfcwater layer
  integer :: icelayers ! # of sfcwater layers that contain some ice
  integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)
  integer :: nlev_new  ! new number of sfcwater layers after adjustment

  real :: totsnow  ! sum of mass over sfcwater layers [kg/m^2]
  real :: wtnew    ! weight for new sfcwater layer when adjusting layer thickness
  real :: wtold    ! weight for old sfcwater layer when adjusting layer thickness
  real :: dwtold   ! change in wtold for partial mass transfer from old layer
  real :: wdiff    ! change in sfcwater mass when adjusting layer thickness [kg/m^2]

  ! Local parameters

  real, parameter :: snowmin = 11.      ! min sfcwater layer mass with multiple layers [kg/m^2] 

  real, save, dimension(10,10) :: thick  ! snowlayer thickness scaling factor

  integer, parameter :: iland_print = 0

  real :: sfcwat_tempk, sfcwat_fracliq

  data thick(1:10, 1)/  1., .00, .00, .00, .00, .00, .00, .00, .00, .00/
  data thick(1:10, 2)/ .50, .50, .00, .00, .00, .00, .00, .00, .00, .00/
  data thick(1:10, 3)/ .25, .50, .25, .00, .00, .00, .00, .00, .00, .00/
  data thick(1:10, 4)/ .17, .33, .33, .17, .00, .00, .00, .00, .00, .00/
  data thick(1:10, 5)/ .10, .20, .40, .20, .10, .00, .00, .00, .00, .00/
  data thick(1:10, 6)/ .07, .14, .29, .29, .14, .07, .00, .00, .00, .00/
  data thick(1:10, 7)/ .05, .09, .18, .36, .18, .09, .05, .00, .00, .00/
  data thick(1:10, 8)/ .03, .07, .13, .27, .27, .13, .07, .03, .00, .00/
  data thick(1:10, 9)/ .02, .04, .09, .18, .34, .18, .09, .04, .02, .00/
  data thick(1:10,10)/ .02, .03, .06, .13, .26, .26, .13, .06, .03, .02/

  ! Re-distribute mass, energy, and depth among layers to maintain prescribed 
  ! distribution of mass

  ! Count up all existing snow layers that are not totally liquid

  icelayers = 0

  totsnow = 0.

  do k = 1,nlev_sfcwater
     totsnow = totsnow + sfcwater_mass(k)

     if (sfcwater_mass(k) > 1.e-9 .and.  &
         energy_per_m2(k) < sfcwater_mass(k) * alli) then

        icelayers = icelayers + 1
     endif
  enddo

  ! If total sfcwater mass is very small or zero, make it zero and return

  if (totsnow < 1.e-9) then
     nlev_sfcwater = 0
  endif

  if (nlev_sfcwater == 0) then
     sfcwater_mass  (1:nzs) = 0.
     sfcwater_energy(1:nzs) = 0.
     sfcwater_depth (1:nzs) = 0.
     return
  endif

  ! If more than one snowcover layer is allowed, perform adjustment of layers if needed

  ! Find maximum number of layers for which thinnest layer (top and bottom) would 
  ! be thicker than snowmin.

  maxlayers = 1

  do while (maxlayers < nzs .and. maxlayers < 10 .and.  &
            thick(1,maxlayers+1) * totsnow > snowmin)

     maxlayers = maxlayers + 1
  enddo

  ! Determine new number of layers.  This number may be at most one greater than
  ! nlev_sfcwater and one greater than icelayers, but no greater than maxlayers.

  nlev_new = min(nlev_sfcwater+1, icelayers+1, maxlayers)

  if (nlev_new == 1 .and. nlev_sfcwater == 1) return

  ! Set index for first old layer; set transfer weights for first old layer

  kold = 1
  wtold = 1. ! fraction of "remaining" mass in old layer

  ! Loop over new set of layers

  do k = 1,nlev_new

     ! To begin current new layer, set fraction of "unfilled" mass in new layer to 1.

     wtnew = 1.

     ! Set mass of current new layer (already determined)

     mass_new(k) = totsnow * thick(k,nlev_new)

     ! Initialize energy, depth, and s/w flux of current new layer to zero

     energy_new(k) = 0.
     depth_new (k) = 0.

     10 continue

     ! Compute "unfilled" mass in current new layer minus remaining mass in
     ! current old layer

     wdiff = wtnew * mass_new(k) - wtold * sfcwater_mass(kold)

     ! Check sign of wdiff

     if (wdiff > 0.) then

        ! If "unfilled" mass in current new layer exceeds remaining mass in current old
        ! layer, transfer all of remaining energy, depth, and s/w flux from old layer

        energy_new(k) = energy_new(k) + wtold * energy_per_m2 (kold)
        depth_new (k) = depth_new (k) + wtold * sfcwater_depth(kold)

        ! Reduce fraction of "unfilled" mass in current new layer, jump to next old
        ! layer, and return wtold to 1.0 to indicate no mass in old layer yet removed.

        wtnew = wtnew - wtold * sfcwater_mass(kold) / mass_new(k)
        kold = kold + 1
        wtold = 1.

        ! If old-layer counter does not exceed top old layer, repeat transfer operation 
        ! for current old layer

        if (kold <= nlev_sfcwater) go to 10

     else

        ! If "unfilled" mass in current new layer is less than remaining mass in current
        ! old layer, transfer only the portion of remaining energy, depth, and s/w flux
        ! from old layer that fit into new layer.

        ! Change in wtold

        dwtold = wtnew * mass_new(k) / sfcwater_mass(kold)

        wtold = wtold - dwtold

        ! Energy, depth, and s/w flux transfer to current new layer

        energy_new(k) = energy_new(k) + dwtold * energy_per_m2 (kold)
        depth_new (k) = depth_new (k) + dwtold * sfcwater_depth(kold)

     endif

  enddo

  ! Now that mass, energy, depth, and s/w flux have been transferred to new layers,
  ! copy back to original arrays

  do k = 1,nlev_new
     sfcwater_mass(k)   = mass_new(k)
     sfcwater_energy(k) = energy_new(k) / sfcwater_mass(k)
     sfcwater_depth(k)  = depth_new(k)
  enddo

  if (nlev_new < nlev_sfcwater) then
     do k = nlev_new+1, nlev_sfcwater
        sfcwater_mass(k)   = 0.0
        sfcwater_energy(k) = 0.0
        sfcwater_depth(k)  = 0.0
     enddo
  endif

  nlev_sfcwater = nlev_new

  ! Replace sfcwater energy limiter from earlier code versions with energy
  ! check and message.  If message ever gets printed, investigate reasons.

  do k = 1,nlev_sfcwater
     if (sfcwater_energy(k) > 6.e5 .or. sfcwater_energy(k) < -2.5e5) then

        write(*,*) ' '
        write(*,*) 'Sfcwater energy is outside allowable range. '
        write(*,*) 'iwsfc,iland,k,lat,lon = ',iwsfc,iland,k,glatw,glonw
        write(*,*) 'sfcwater_energy = ',sfcwater_energy(k)
        write(*,*) 'p1',energy_new(k),sfcwater_mass(k),nlev_sfcwater,nlev_new
        write(*,*) 'p2',kold, energy_per_m2(kold)
        stop 'stop sfcwater energy '

     endif
  enddo

  end subroutine sfcwater_adjust

!===============================================================================

  subroutine remove_runoff(iland, iwsfc, leaf_class, sfcwater_mass, sfcwater_energy, &
                           sfcwater_depth, runoff)

  use leaf_coms,   only: nzs, dt_leaf
  use consts_coms, only: alli, cliq
  use misc_coms,   only: io6
  use therm_lib,   only: qtk

  implicit none

  integer, intent(in) :: iland
  integer, intent(in) :: iwsfc
  integer, intent(in) :: leaf_class

  real, intent(inout) :: sfcwater_mass
  real, intent(inout) :: sfcwater_energy
  real, intent(inout) :: sfcwater_depth
  real, intent(inout) :: runoff

  ! Local variables

  real :: sfcwater_tempk
  real :: sfcwater_fracliq
  real :: energy_per_m2
  real :: sfcwater_mass_thresh
  real :: wfree
  real :: qrunoff, drunoff

  real, parameter :: runoff_time = 86400.  ! time scale for surface-water runoff [s]

  ! Set a threshold value of sfcwater_mass for runoff to occur that is based on
  ! leaf_class

  if (leaf_class == 17 .or. leaf_class == 20) then
     sfcwater_mass_thresh = 1000.  ! 1000 kg/m^2 equivalent to 1.0 m depth
  else
     sfcwater_mass_thresh = 1.     ! 1 kg/m^2 equivalent to 1.0 mm depth
  endif

  if (sfcwater_mass <= sfcwater_mass_thresh) then
     runoff = 0.0
     return
  endif

  call qtk(sfcwater_energy, sfcwater_tempk, sfcwater_fracliq)

  ! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

  if (sfcwater_fracliq > .10) then
     wfree = sfcwater_mass * (sfcwater_fracliq - .10) / .90

     ! Fraction of wfree removed as runoff this time step
     ! (In future development, runoff_time should be function of topography)

     runoff = wfree * dt_leaf / runoff_time

     ! Evaluate energy and depth contained in runoff (which is in liquid phase)

     qrunoff = runoff * (cliq * (sfcwater_tempk - 273.15) +  alli)
     drunoff = runoff * 0.001

     ! Subtract runoff from sfcwater mass, energy, and depth

     energy_per_m2 = sfcwater_mass * sfcwater_energy
     sfcwater_mass = sfcwater_mass - runoff
     sfcwater_energy = (energy_per_m2 - qrunoff) / sfcwater_mass
     sfcwater_depth = sfcwater_depth - drunoff
  endif
  
  end subroutine remove_runoff

!===============================================================================

  subroutine grndvap(iland, rhos, canshv, nlev_sfcwater, surface_ssh, &
                     ground_shv, sfcwater_energy, soil_water, soil_energy, &
                     head, specifheat_drysoil)

  use mem_land,    only: nzg
  use consts_coms, only: grav, rvap
  use misc_coms,   only: io6
  use therm_lib,   only: qtk, rhovsil, qwtk

  implicit none

  integer, intent(in)  :: iland           ! current land cell number 
  real,    intent(in)  :: rhos            ! air density [kg/m^3]
  real,    intent(in)  :: canshv          ! canopy vapor spec hum [kg_vap/kg_air]
  integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
  real,    intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]
  real,    intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
  real,    intent(in)  :: sfcwater_energy ! [J/kg]
  real,    intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
  real,    intent(in)  :: soil_energy     ! [J/m^3]
  real,    intent(in)  :: head            ! hydraulic head of top soil layer [m]
  real,    intent(in)  :: specifheat_drysoil

  ! Local variables

  real :: tempk     ! surface water temp [K]
  real :: fracliq   ! fraction of surface water in liquid phase
  real :: can_rhov  ! canopy water vapor density [kg_vap/m^3]
  real :: sfc_rhovs ! ground sfc saturation vapor density [kg_vap/m^3]
  real :: gnd_rhov  ! ground sfc evaporative vapor density [kg_vap/m^3]

  if (nlev_sfcwater > 0) then

     ! With surface water (or snowcover) present, sfc_rhovs is the saturation
     ! vapor density of the water/snow surface and is used for computing
     ! dew/frost formation and sfcwater evaporation.

     call qtk(sfcwater_energy,tempk,fracliq)
     sfc_rhovs = rhovsil(tempk-273.15)
     surface_ssh = sfc_rhovs / rhos

  else

     ! Without snowcover, sfc_rhovs is the saturation vapor density at the
     ! temperature of the soil surface, and is used for computing dew/frost
     ! formation only.  

     call qwtk(soil_energy,soil_water*1.e3,specifheat_drysoil,tempk,fracliq)
     sfc_rhovs = rhovsil(tempk-273.15)

     ! gnd_rhov is the effective saturation vapor density of soil and is used
     ! for computing soil evaporation.

     can_rhov = canshv * rhos

     call grndvap_ab(iland,tempk,soil_water,head,can_rhov,sfc_rhovs,gnd_rhov)

     surface_ssh = sfc_rhovs / rhos
     ground_shv  = gnd_rhov  / rhos

  endif

  end subroutine grndvap

  !===============================================================================

  subroutine grndvap_ab(iland,tempk,soil_water,head,can_rhov,sfc_rhovs,gnd_rhov)

  use mem_land,  only: nzg
  use consts_coms, only: grav, rvap
  use misc_coms,   only: io6

  implicit none

  integer, intent(in) :: iland       ! current land cell number 

  real, intent(in)  :: tempk       ! soil temperature [K]
  real, intent(in)  :: soil_water  ! soil water content [vol_water/vol_tot]
  real, intent(in)  :: head        ! hydraulic head [m]
  real, intent(in)  :: can_rhov    ! canopy vapor density [kg_vap/m^3]
  real, intent(in)  :: sfc_rhovs   ! surface saturation vapor density [kg_vap/m^3]
  real, intent(out) :: gnd_rhov    ! ground equilibrium vapor density [kg_vap/m^3]

 ! Local parameter

 real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

 ! Local variables

  real :: alpha   ! "alpha" term in Lee and Pielke (1993)
  real :: beta    ! "beta" term in Lee and Pielke (1993)
  real, save :: sfldcap(12)  ! soil water field capacity [vol_water/vol_tot]

  !HARDWIRE NTS NOW, BUT LATER COMPUTE SFLDCAP
  integer, parameter :: nts = 5
  data sfldcap/.135,.150,.195,.255,.240,.255,.322,.325,.310,.370,.367,.535/

  ! Without snowcover, gnd_rhov is the effective saturation mixing
  ! ratio of soil and is used for soil evaporation.  First, compute the
  ! "alpha" term or soil "relative humidity" and the "beta" term.
  ! Use limiter on head so that positive values are not used.

  alpha = exp(gorvap * min(0.,head) / tempk)
  beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nts)) * 3.14159)) ** 2

  gnd_rhov = sfc_rhovs * alpha * beta + (1. - beta) * can_rhov

  ! Require that gnd_rhov be at least slightly less than sfc_rhovs

  gnd_rhov = min(gnd_rhov, 0.999 * sfc_rhovs)

  end subroutine grndvap_ab

  !===============================================================================

  subroutine sfcrad_land(iland, leaf_class, xewl, yewl, zewl, wnxl, wnyl, wnzl, &
                       rshort, rlong, rlongup, rlong_albedo, albedo_beam,     &
                       sfcwater_energy, sfcwater_depth, rshort_s,             &
                       rshort_g, rshort_v, rlong_g, rlong_s, rlong_v,         &
                       nlev_sfcwater, veg_temp, veg_fracarea, veg_height,     &
                       veg_albedo, snowfac, vf, cosz,                         &
                       soil_energy, soil_water, wresid_vg, wsat_vg,           &
                       specifheat_drysoil, sand                               )

  use leaf_coms,   only: nzs, emisv, emisg, emisw
  use consts_coms, only: stefan, eradi
  use misc_coms,   only: io6
  use mem_radiate, only: sunx, suny, sunz
  use therm_lib,   only: qwtk, qtk

  implicit none

  integer, intent(in)  :: iland         ! index of current land cell
  integer, intent(in)  :: leaf_class    ! leaf class
  real,    intent(in)  :: xewl          ! land cell earth x coordinate [m]
  real,    intent(in)  :: yewl          ! land cell earth y coordinate [m]
  real,    intent(in)  :: zewl          ! land cell earth z coordinate [m]
  real,    intent(in)  :: wnxl          ! land cell norm unit vec x comp [m]
  real,    intent(in)  :: wnyl          ! land cell norm unit vec y comp [m]
  real,    intent(in)  :: wnzl          ! land cell norm unit vec z comp [m]
  real,    intent(in)  :: rshort        ! downward surface incident s/w rad flux [W/m^2]
  real,    intent(in)  :: rlong         ! downward surface incident l/w rad flux [W/m^2]
  real,    intent(out) :: rlongup       ! upward sfc l/w rad flux [W/m^2]
  real,    intent(out) :: rlong_albedo  ! albedo for upward sfc l/w
  real,    intent(out) :: albedo_beam   ! land cell albedo (beam)
  real,    intent(in)  :: sfcwater_energy(nzs) ! surface water internal energy [J/kg]
  real,    intent(in)  :: sfcwater_depth(nzs)  ! surface water depth [m]
  real,    intent(out) :: rshort_s(nzs) ! s/w net rad flux to sfc water [W/m^2]
  real,    intent(out) :: rshort_g      ! s/w net rad flux to soil [W/m^2]
  real,    intent(out) :: rshort_v      ! s/w net rad flux to veg [W/m^2]
  real,    intent(out) :: rlong_g       ! l/w net rad flux to soil [W/m^2]
  real,    intent(out) :: rlong_s       ! l/w net rad flux to sfc water [W/m^2]
  real,    intent(out) :: rlong_v       ! l/w net rad flux to veg [W/m^2]
  integer, intent(in)  :: nlev_sfcwater ! # of active surface water layers
  real,    intent(in)  :: veg_temp      ! veg temp [K]
  real,    intent(in)  :: veg_fracarea  ! veg fractional area coverage
  real,    intent(in)  :: veg_height    ! veg height [m]
  real,    intent(in)  :: veg_albedo    ! veg albedo
  real,    intent(out) :: snowfac       ! fractional veg burial by snowcover
  real,    intent(out) :: vf            ! fractional coverage of non-buried part of veg
  real,    intent(out) :: cosz          ! solar zenith angle for land cell
  real,    intent(in)  :: soil_energy   ! soil internal energy [J/m^3]
  real,    intent(in)  :: soil_water    ! soil water content [vol_water/vol_tot]
  real,    intent(in)  :: wresid_vg
  real,    intent(in)  :: wsat_vg
  real,    intent(in)  :: specifheat_drysoil
  real,    intent(in)  :: sand

  ! Local variables

  integer :: k             ! vertical index over surface water layers

  real :: soil_tempk       ! soil temp [K]
  real :: soil_fracliq     ! fraction of soil water in liquid phase
  real :: sfcwater_tempk   ! surface water temp [K]
  real :: sfcwater_fracliq ! fraction of surface water in liquid phase
  real :: rlonga_v         ! longwave radiative flux from atm  to veg  [W/m^2]
  real :: rlonga_s         ! longwave radiative flux from atm  to snow [W/m^2]
  real :: rlonga_g         ! longwave radiative flux from atm  to soil [W/m^2]
  real :: rlongv_a         ! longwave radiative flux from veg  to atm  [W/m^2]
  real :: rlongv_s         ! longwave radiative flux from veg  to snow [W/m^2]
  real :: rlongv_g         ! longwave radiative flux from veg  to soil [W/m^2]
  real :: rlongs_a         ! longwave radiative flux from snow to atm  [W/m^2]
  real :: rlongs_v         ! longwave radiative flux from snow to veg  [W/m^2]
  real :: rlongg_a         ! longwave radiative flux from soil to atm  [W/m^2]
  real :: rlongg_v         ! longwave radiative flux from soil to veg  [W/m^2]
  real :: fracabs          ! fraction of rshort that is absorbed by snowcover
  real :: soil_watfrac_ul  ! soil water fraction (unlimited)
  real :: soil_watfrac     ! soil water fraction (limited)
  real :: vfc       ! 1 - vf
  real :: fc50      ! minimum of soil water fraction and 50%
  real :: alg_dry   ! ground (soil) albedo for dry soil
  real :: alg       ! ground (soil) albedo
  real :: als       ! snowcover albedo (needs better formula based on age of snow)
  real :: alv       ! veg albedo
  real :: abss      ! fraction s/w rad absorbed into top sfcwater layer
  real :: fractrans ! fraction of s/w rad flux transmitted through snowcover layer
  real :: absg      ! fraction of rshort that is absorbed by ground (soil)
  real :: emv       ! veg emissivity
  real :: glong     ! soil l/w rad emission [W/m^2]
  real :: slong     ! sfc water l/w rad emission [W/m^2]
  real :: vlong     ! veg l/w rad emission [W/m^2]

  ! This routine is called twice (for each land cell) by the radiation
  ! parameterization driver, performing exactly the same operations on both calls.
  ! All that is used from the first call are net surface albedo and upward
  ! longwave radiative flux for each land cell.  The radiation parameterization
  ! carries out atmospheric radiative transfer computations following this first 
  ! call using the albedo and upward longwave flux.  The second call to this 
  ! subroutine is made after the atmospheric radiative fluxes are computed.  
  ! This call provides net radiative fluxes to vegetation, snowcover, and soil in
  ! each land cell, plus functions of snowcover, all of which are used in leaf4.

  !-------------------------------------------------------------------------------
  ! Two forms - but still need to consider shadowing and exact energy conservation.
  ! A correct formulation must exchange radiative energy between model columns.

  ! Compute solar zenith angle for land cells

  cosz = (xewl * sunx + yewl * suny + zewl * sunz) * eradi

  ! Compute solar incidence angle for land cells (accounts for local topo slope)

  cosz = wnxl * sunx + wnyl * suny + wnzl * sunz

  ! COSZ IS NOT CURRENTLY USED IN THIS SUBROUTINE
  !-------------------------------------------------------------------------------
  ! Change for LEAF4: rshort_s(1:nzs) is always set to zero here, and rshort_g
  ! always gets all shortwave that is absorbed by the surface, whether it is
  ! sfcwater, soil, or a combination.  Subroutine leaf4_canopy further sorts out
  ! where rshort_g should go.

  rshort_s(1:nzs) = 0.

  if (nlev_sfcwater == 0) then

     ! Case with no surface water

     ! Diagnose soil temperature and liquid fraction

     call qwtk(soil_energy,soil_water*1.e3,  &
               specifheat_drysoil,soil_tempk,soil_fracliq)

     ! Shortwave radiation calculations

     soil_watfrac_ul = (soil_water - wresid_vg) / (wsat_vg - wresid_vg)

     soil_watfrac = min(1.0,max(0.001,soil_watfrac_ul))

     fc50 = min(.50, soil_watfrac)

     if (leaf_class == 2) then

        alg = .80  ! Firn/glacier albedo

     else

        if (sand > 0.7) then
           alg_dry = .41  ! Higher albedo for sandy soil
        else
           alg_dry = .31  ! Soils other than sand and loamy sand
        endif

        alg = alg_dry - .34 * fc50  ! ground (soil) albedo

     endif

     snowfac = 0.

     alv = veg_albedo

     vf = veg_fracarea
     vfc = 1. - vf

     albedo_beam = vf * alv + vfc * vfc * alg

     rshort_g = rshort * vfc * (1. - alg)
     rshort_v = rshort * vf * (1. - alv + vfc * alg)
    !rshort_a = rshort * albedo_beam

     ! Longwave radiation calculations

     emv   = emisv(leaf_class)
     glong = emisg * stefan * soil_tempk ** 4
     vlong = emv * stefan * veg_temp ** 4

     rlonga_v = rlong * vf * (emv + vfc * (1. - emisg))
     rlonga_g = rlong * vfc * emisg
     rlongv_g = vlong * vf * emisg
     rlongv_a = vlong * vf * (2. - emisg - vf + emisg * vf)
     rlongg_v = glong * vf * emv
     rlongg_a = glong * vfc

     rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - emisg))

     rlongup = rlongv_a + rlongg_a

     rlong_g = rlonga_g - rlongg_a + rlongv_g - rlongg_v
     rlong_s = 0.
     rlong_v = rlonga_v - rlongv_a + rlongg_v - rlongv_g

  else

     ! Case with surface water

     ! Diagnose surface water temperature and liquid fraction

     call qtk(sfcwater_energy(nlev_sfcwater),sfcwater_tempk,sfcwater_fracliq)

     ! Shortwave radiation calculations

     if (leaf_class == 2) then

        ! Sfcwater albedo ALS is set to .80 over firn/glacier

        als = .80  ! Firn/glacier albedo

     else

        ! Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid
        ! to .5 for all-ice

        als = .5 - .36 * sfcwater_fracliq

     endif

     ! Sum over sfcwater layers to get total depth

     snowfac = 0.

     do k = 1,nlev_sfcwater
        snowfac = snowfac + sfcwater_depth(k)
     enddo

     snowfac = snowfac / max(.001,veg_height)

     ! If sfcwater (snowcover) is deep enough to cover most vegetation, assume that 
     ! all is covered

     if (snowfac > .9) snowfac = 1.   

     alv = veg_albedo

     vf = veg_fracarea * (1. - snowfac)
     vfc = 1. - vf

     albedo_beam = vf * alv + vfc * vfc * als

     rshort_g = rshort * vfc * (1. - als)
     rshort_v = rshort * vf * (1. - alv + vfc * als)
    !rshort_a = rshort * albedo_beam

     ! Longwave radiation calculations

     emv   = emisv(leaf_class)
     slong = emisw * stefan * sfcwater_tempk ** 4
     vlong = emv * stefan * veg_temp ** 4

     rlonga_v = rlong * vf * (emv + vfc * (1. - emisw))
     rlonga_s = rlong * vfc * emisw
     rlongv_s = vlong * vf * emisw
     rlongv_a = vlong * vf * (2. - emisw - vf + emisw * vf)
     rlongs_v = slong * vf * emv
     rlongs_a = slong * vfc

     rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - emisw))

     rlongup = rlongv_a + rlongs_a

     rlong_g = 0.
     rlong_s = rlonga_s - rlongs_a + rlongv_s - rlongs_v
     rlong_v = rlonga_v - rlongv_a + rlongs_v - rlongv_s

  endif

  end subroutine sfcrad_land

End Module leaf4_surface
