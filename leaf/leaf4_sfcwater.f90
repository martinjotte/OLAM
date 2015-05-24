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
Module leaf4_sfcwater

Contains

subroutine sfcwater(iwl, icomb, nlev_sfcwater, ntext_soil,             &
                    soil_water, soil_energy, soil_tempk, soil_fracliq, &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,    &
                    sfcwater_tempk, sfcwater_fracliq, energy_per_m2,   &
                    head1, vf, wfree1, qwfree1, dwfree1,               &
                    snowmin_expl                                       )

use leaf_coms,   only: nzg, nzs, dt_leaf, slz, dslz, dslzi, dslzo2,  &
                       slmsts, soilcond0, soilcond1, soilcond2, slcpd
                       
use consts_coms, only: alvi, cice, cliq, alli
use misc_coms,   only: io6
use leaf4_plot,  only: leaf_plot

implicit none

integer, intent(in)    :: iwl                ! current land cell index
integer, intent(in)    :: icomb              ! implicit heat balance flag [0=no, 1=yes]
integer, intent(inout) :: nlev_sfcwater      ! # of active sfc water levels
integer, intent(in)    :: ntext_soil   (nzg) ! soil textural class

real, intent(inout) :: soil_water      (nzg) ! soil water [water_vol/total_vol]
real, intent(inout) :: soil_energy     (nzg) ! soil internal energy [J/m^3]
real, intent(in)    :: soil_tempk      (nzg) ! soil temperature [K]
real, intent(in)    :: soil_fracliq    (nzg) ! fraction of soil water in liq phase

real, intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
real, intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
real, intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
real, intent(inout) :: energy_per_m2   (nzs) ! sfcwater energy [J/m^2]

real, intent(inout) :: head1   ! UBC hydraulic head for soil model [m]
real, intent(in)    :: vf      ! fractional coverage of non-buried part of veg
real, intent(out)   :: wfree1  ! free liquid in lowest sfcwater layer [kg/m^2]
real, intent(out)   :: qwfree1 ! energy carried by wfree1 [J/m^2]
real, intent(out)   :: dwfree1 ! depth carried by wfree1 [m]
real, intent(in)    :: snowmin_expl ! min sfcwater mass for explicit heat xfer [kg/m^2]

! Local variables

real :: hxfers        (nzs+1) ! sfcwater heat xfer [J/m2] 
real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]

integer :: k         ! vertical index over sfcwater layers
integer :: kold      ! vertical index of adjacent lower sfcwater layer
integer :: icelayers ! # of sfcwater layers that contain some ice
integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)

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
real :: wtnew     ! weight for new sfcwater layer when adjusting layer thickness
real :: wtold     ! weight for old sfcwater layer when adjusting layer thickness
real :: dwtold    ! change in wtold for partial mass transfer from old layer
real :: wdiff     ! change in sfcwater mass when adjusting layer thickness [kg/m^2]
real :: soilcond  ! soil thermal conductivity [W/(K m)]
real :: waterfrac ! soil water fraction in soil layer [vol_water/vol_total]
real :: tempk     ! Kelvin temperature of snowcover, limited to max of 273.15 [K]
real :: flmin     ! lower bound on sfcwater_fracliq(1) in balance with soil
real :: flmax     ! upper bound on sfcwater_fracliq(1) in balance with soil
real :: specvol   ! specific volume of sfcwater involved in vapor xfer [m^3/kg]
real :: wcap_min  ! minimum surface water water [kg/m^2]

! Local parameters

integer, parameter :: iwl_print = 0
integer :: linit, lframe

! Initialize free water values to zero

wfree1  = 0.
qwfree1 = 0.
dwfree1 = 0.

head1 = 0.

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

   waterfrac = soil_water(nzg) / slmsts(ntext_soil(nzg))
   soilcond =        soilcond0(ntext_soil(nzg))  &
      + waterfrac * (soilcond1(ntext_soil(nzg))  &
      + waterfrac *  soilcond2(ntext_soil(nzg))  )
   soil_rfactor = dslzo2(nzg) / soilcond

! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
! layer, and apply to both.

   hxfergs = dt_leaf * (soil_tempk(nzg) - sfcwater_tempk(1))   &
           / (soil_rfactor + sfcwater_rfactor(1))

   energy_per_m2(1) = energy_per_m2(1) + hxfergs

   soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)

endif

! Take inventory of sfcwater_mass(1) now that exchange with canopy is
! complete.  It is possible for sfcwater_mass(1) to be slightly negative
! at this point.  If sfcwater_mass(1) is below threshold, transfer sfcwater
! mass and energy to soil, set sfcwater quantities to zero and return.

wcap_min = dt_leaf * 1.e-6 ! 1.e-6 is 1 mm pcp/dew in 11 days    [kg/m^2]

if (sfcwater_mass(1) < wcap_min) then

   soil_water (nzg) = soil_water (nzg) + dslzi(nzg) * sfcwater_mass(1) * .001
   soil_energy(nzg) = soil_energy(nzg) + dslzi(nzg) * energy_per_m2(1)

   sfcwater_mass  (1) = 0.
   sfcwater_energy(1) = 0.
   sfcwater_depth (1) = 0.
   energy_per_m2  (1) = 0.

   sfcwater_tempk  (1) = soil_tempk(nzg)
   sfcwater_fracliq(1) = 0.

   nlev_sfcwater = 0

   return

elseif (nlev_sfcwater < 1) then

! If sufficient sfcwater mass remains and nlev_sfcwater is less than 1,
! set it to 1.

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

      call sfcwater_soil_comb(iwl,ntext_soil(nzg), &
                              soil_water(nzg),soil_energy(nzg), &
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

head1 = .001 * sum(sfcwater_mass(1:nlev_sfcwater))

if (iwl == iwl_print) then

   call leaf_plot(iwl,             &
                  nlev_sfcwater,   &
                  linit           = 1,               &
                  lframe          = 1,               &
                  sfcwater_mass   = sfcwater_mass,   & 
                  sfcwater_energy = sfcwater_energy, & 
                  energy_per_m2   = energy_per_m2,   & 
                  sfcwater_depth  = sfcwater_depth   ) 

endif

return
end subroutine sfcwater

!===============================================================================

subroutine sfcwater_soil_comb(iwl, ntext_soil, soil_water, soil_energy, &
                              sfcwater_mass, energy_per_m2, &
                              sfcwater_tempk, sfcwater_fracliq)

use leaf_coms,   only: nzg, dslz, dslzi, slcpd
                       
use consts_coms, only: cice, cliq, alli
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iwl              ! current land cell index
integer, intent(in) :: ntext_soil       ! soil textural class

real, intent(in)    :: soil_water       ! soil water [water_vol/total_vol]
real, intent(inout) :: soil_energy      ! soil internal energy [J/m^3]
real, intent(in)    :: sfcwater_mass    ! surface water mass [kg/m^2]
real, intent(inout) :: energy_per_m2    ! sfcwater energy [J/m^2]
real, intent(inout) :: sfcwater_tempk   ! surface water temp [K]
real, intent(inout) :: sfcwater_fracliq ! fraction of sfc water in liq phase

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

  hcapsoil = slcpd(ntext_soil) * dslz(nzg)

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

  return
  end subroutine sfcwater_soil_comb

!===============================================================================

subroutine sfcwater_adjust(iwl, nlev_sfcwater, sfcwater_mass, sfcwater_energy, &
                           sfcwater_depth, energy_per_m2, glatw, glonw)

use leaf_coms, only: nzs, dt_leaf
                       
use consts_coms, only: alli
use misc_coms, only: io6

implicit none

integer, intent(in)    :: iwl              ! current land cell index
integer, intent(inout) :: nlev_sfcwater    ! # of active sfc water levels

real, intent(inout) :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy(nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth (nzs) ! surface water depth [m]
real, intent(inout) :: energy_per_m2  (nzs) ! sfcwater energy [J/m^2]
real, intent(in) :: glatw         ! Latitude of land cell 'center' [deg]
real, intent(in) :: glonw         ! Longitude of land cell 'center' [deg]

! Local variables

real :: mass_new   (nzs) ! mass of new set of sfcwater layers [kg/m^2]
real :: energy_new (nzs) ! energy of new set of sfcwater layers [J/m^2]
real :: depth_new  (nzs) ! depth of new set of sfcwater layers [m]

integer :: k         ! vertical index over sfcwater layers
integer :: kold      ! vertical index of adjacent lower sfcwater layer
integer :: icelayers ! # of sfcwater layers that contain some ice
integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)
integer :: nlev_new  ! new number of sfcwater layers after adjustment

real :: snden    ! sfcwater density [kg/m^3]
real :: totsnow  ! sum of mass over sfcwater layers [kg/m^2]
real :: sndenmin ! minimum sfcwater density [kg/m^3]
real :: wtnew    ! weight for new sfcwater layer when adjusting layer thickness
real :: wtold    ! weight for old sfcwater layer when adjusting layer thickness
real :: dwtold   ! change in wtold for partial mass transfer from old layer
real :: wdiff    ! change in sfcwater mass when adjusting layer thickness [kg/m^2]
real :: glatwl   ! latitude of land cell 'center'
real :: glonwl   ! longitude of land cell 'center'

! Local parameters

real, parameter :: snowmin = 11.      ! min sfcwater layer mass with multiple layers [kg/m^2] 

real, save, dimension(10,10) :: thick  ! snowlayer thickness scaling factor

integer, parameter :: iwl_print = 0

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

   energy_new (k) = 0.
   depth_new  (k) = 0.

10 continue

! Compute "unfilled" mass in current new layer minus remaining mass in
! current old layer

   wdiff = wtnew * mass_new(k) - wtold * sfcwater_mass(kold)

! Check sign of wdiff

   if (wdiff > 0.) then

! If "unfilled" mass in current new layer exceeds remaining mass in current old
! layer, transfer all of remaining energy, depth, and s/w flux from old layer

      energy_new (k) = energy_new (k) + wtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + wtold * sfcwater_depth(kold)

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

      energy_new (k) = energy_new (k) + dwtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + dwtold * sfcwater_depth(kold)

   endif

enddo

! Now that mass, energy, depth, and s/w flux have been transferred to new layers,
! copy back to original arrays

do k = 1,nlev_new
   sfcwater_mass(k)   = mass_new(k)
   sfcwater_energy(k) = energy_new(k) / sfcwater_mass(k)
   sfcwater_depth(k)  = depth_new(k)
enddo

nlev_sfcwater = nlev_new

! Replace sfcwater energy limiter from earlier code versions with energy
! check and message.  If message ever gets printed, investigate reasons.

do k = 1,nlev_sfcwater
   if (sfcwater_energy(k) > 6.e5 .or. sfcwater_energy(k) < -2.5e5) then

      write(io6,*) ' '
      write(io6,*) 'Sfcwater energy is outside allowable range. '
      write(io6,*) 'iwl,k,lat,lon = ',iwl,k,glatw,glonw
      write(io6,*) 'sfcwater_energy = ',sfcwater_energy(k)
      write(io6,*) 'p1',energy_new(k),sfcwater_mass(k),nlev_sfcwater,nlev_new
      write(io6,*) 'p2',kold, energy_per_m2(kold)
      stop 'stop sfcwater energy'
   endif

enddo

return
end subroutine sfcwater_adjust

!===============================================================================

subroutine remove_runoff(iwl, ksn, sfcwater_fracliq, sfcwater_mass,   &
     sfcwater_tempk, sfcwater_energy, sfcwater_depth, runoff, qrunoff)

  use leaf_coms, only: nzs, dt_leaf
  use consts_coms, only: alli, cliq
  use misc_coms, only: io6

  implicit none

  real, parameter :: runoff_time = 86400.

  integer, intent(in) :: iwl
  integer, intent(in) :: ksn

  real, intent(inout) :: sfcwater_fracliq(nzs)
  real, intent(inout) :: sfcwater_tempk  (nzs)
  real, intent(inout) :: sfcwater_mass   (nzs)
  real, intent(inout) :: sfcwater_energy (nzs)
  real, intent(inout) :: sfcwater_depth  (nzs)
  real, intent(out)   :: runoff
  real, intent(out)   :: qrunoff

  ! Get rid of runoff

  runoff = 0.0
  qrunoff = 0.0

  if(ksn >= 1)then
     if(sfcwater_fracliq(ksn) > 0.1)then
        call qtk(sfcwater_energy(ksn), sfcwater_tempk(ksn),   &
             sfcwater_fracliq(ksn))
        runoff = sfcwater_mass(ksn) * (sfcwater_fracliq(ksn) - 0.1) /   &
             (0.9 * runoff_time) * dt_leaf
        qrunoff = runoff * (cliq * (sfcwater_tempk(ksn) - 273.15) +  alli)
        
        sfcwater_energy(ksn) = (sfcwater_energy(ksn) *   &
             sfcwater_mass(ksn) - qrunoff ) / (sfcwater_mass(ksn) - runoff)
        sfcwater_mass(ksn) = sfcwater_mass(ksn) - runoff
        sfcwater_depth(ksn) = sfcwater_depth(ksn) - 0.001 * runoff
     endif
  endif
  
  return
end subroutine remove_runoff

End Module leaf4_sfcwater
