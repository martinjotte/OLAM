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
Module leaf3_sfcwater

Contains

subroutine sfcwater(iwl, nlev_sfcwater, ntext_soil,                  &
                    soil_rfactor, soil_water, soil_energy,           &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,  &
                    soil_tempk, soil_fracliq, sfcwater_tempk,        &
                    sfcwater_fracliq, energy_per_m2, rshort_s,       &
                    hxfersc, wxfersc, rlong_g, rlong_s,              &
                    pcpg, qpcpg, dpcpg, wshed, qwshed, head1,        &
                    vf, wfree1, qwfree1, dwfree1, time8              )

use leaf_coms,   only: nzg, nzs, dt_leaf, slz, dslz, dslzi, dslzo2,  &
                       slmsts, soilcond0, soilcond1, soilcond2, slcpd
                       
use consts_coms, only: alvi, cice, cliq, alli
use misc_coms,   only: io6
use leaf3_plot,  only: leaf_plot

implicit none

integer, intent(in)    :: iwl              ! current land cell index
integer, intent(inout) :: nlev_sfcwater      ! # of active sfc water levels
integer, intent(in)    :: ntext_soil   (nzg) ! soil textural class

real, intent(out)   :: soil_rfactor    (nzg) ! soil thermal resistance [K m^2/W]
real, intent(inout) :: soil_water      (nzg) ! soil water [water_vol/total_vol]
real, intent(inout) :: soil_energy     (nzg) ! soil internal energy [J/m^3]
real, intent(inout) :: soil_tempk      (nzg) ! soil temperature [K]
real, intent(inout) :: soil_fracliq    (nzg) ! fraction of soil water in liq phase

real, intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
real, intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
real, intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
real, intent(inout) :: energy_per_m2   (nzs) ! sfcwater energy [J/m^2]
real, intent(inout) :: rshort_s        (nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(in)    :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(in)    :: rlong_g ! l/w net rad flux to soil [W/m^2]
real, intent(in)    :: rlong_s ! l/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: pcpg    ! new pcp amount this leaf timestep [kg/m^2]
real, intent(in)    :: qpcpg   ! new pcp energy this leaf timestep [J/m^2]
real, intent(in)    :: dpcpg   ! new pcp depth this leaf timestep [m]
real, intent(in)    :: wshed   ! water shed from veg this timestep [kg/m^2]
real, intent(in)    :: qwshed  ! water energy shed from veg this timestep [J/m^2]
real, intent(inout) :: head1   ! UBC hydraulic head for soil model [m]
real, intent(in)    :: vf      ! fractional coverage of non-buried part of veg

real, intent(out) :: wfree1   ! free liquid in lowest sfcwater layer [kg/m^2]
real, intent(out) :: qwfree1  ! energy carried by wfree1 [J/m^2]
real, intent(out) :: dwfree1  ! depth carried by wfree1 [m]

! Local variables

real :: hxfers        (nzs+1) ! sfcwater heat xfer [J/m2] 
real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]

integer :: k         ! vertical index over sfcwater layers
integer :: kold      ! vertical index of adjacent lower sfcwater layer
integer :: icelayers ! # of sfcwater layers that contain some ice
integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)

real :: hxfergs   ! energy transfer from soil to sfcwater this step [J/m^2]
real :: rfac      ! bounded sfcwater thermal resistivity at k=1 [K m^2/W]
real :: snden     ! sfcwater density [kg/m^3]
real :: vegfracc  ! 1 minus veg fractional area
real :: wfree     ! free liquid in sfcwater layer that can percolate out [kg/m^2]
real :: qwfree    ! energy carried by wfree [J/m^2]
real :: dwfree    ! depth carried by wfree [m]
real :: fracstep  ! ratio of leaf timestep to snow density exponential decay time
real :: totsnow   ! sum of mass over sfcwater layers [kg/m^2]
real :: wt        ! sfcwater(1) + soil(nzg) water masses (impl balance) [kg/m^2]
real :: qwt       ! sfcwater(1) + soil(nzg) energies (impl balance) [J/m^2]
real :: soilhcap  ! soil(nzg) heat capacity [J/(m^2 K)]
real :: sndenmin  ! minimum sfcwater density [kg/m^3]
real :: wtnew     ! weight for new sfcwater layer when adjusting layer thickness
real :: wtold     ! weight for old sfcwater layer when adjusting layer thickness
real :: dwtold    ! change in wtold for partial mass transfer from old layer
real :: wdiff     ! change in sfcwater mass when adjusting layer thickness [kg/m^2]
real :: soilcond  ! soil thermal conductivity [W/(K m)]
real :: waterfrac ! soil water fraction in soil layer [vol_water/vol_total]
real :: tempk     ! Kelvin temperature [K]
real :: fracliq   ! fraction of water in liquid phase returned from qwtk
real :: flmin     ! lower bound on sfcwater_fracliq(1) in balance with soil
real :: flmax     ! upper bound on sfcwater_fracliq(1) in balance with soil
real :: specvol   ! specific volume of sfcwater involved in vapor xfer [m^3/kg]
real(kind=8), intent(in) :: time8   ! model time [s]

! Local parameters

real, parameter :: sndenmax = 1000.   ! max sfcwater density [kg/m^3]
real, parameter :: snowmin = 11.      ! min sfcwater layer mass with multiple layers [kg/m^2] 
real, parameter :: snowmin_expl = 10. ! min sfcwater mass for explicit heat xfer [kg/m^2]
real, parameter :: rfac_snowmin = .01 ! min sfcwater rfactor [K m^2/W]

integer, parameter :: iwl_print = 0
integer :: linit, lframe

! Initialize free water values to zero

wfree1  = 0.
qwfree1 = 0.
dwfree1 = 0.

head1 = 0.

! Check whether surface water was present at the beginning of this leaf step

if (nlev_sfcwater > 0) then

! Surface water was present.

! Loop over existing sfcwater layers

   do k = 1,nlev_sfcwater
   
! Sfcwater energy per m^2   

      energy_per_m2(k) = sfcwater_energy(k) * sfcwater_mass(k)

! Compute snow heat resistance times HALF layer depth (sfcwater_rfactor).
! Sfcwater_tempk(k) should be correctly balanced value at this point, so 
! sfcwater_rfactor(k) should have correct value.  Formula applies to snow,
! so limit temperature to no greater than 273.15.

      snden = sfcwater_mass(k) / sfcwater_depth(k)
      tempk = min(273.15,sfcwater_tempk(k))

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

   enddo


! Zero out sfcwater internal heat transfer array at bottom and top surfaces.
! Energy transfer at bottom and top are applied separately.

   hxfers(1) = 0.
   hxfers(nlev_sfcwater+1) = 0. 

! Compute internal sfcwater energy xfer if at least two layers exist [J/m2]

   if (nlev_sfcwater >= 2) then
      do k = 2,nlev_sfcwater
         hxfers(k) = dt_leaf * (sfcwater_tempk(k-1) - sfcwater_tempk(k))  &
                   / (sfcwater_rfactor(k-1) + sfcwater_rfactor(k))      
      enddo
   endif

! Add contributions to sfcwater energy_per_m2 from internal transfers of 
! sensible heat and vapor (latent heat) and from shortwave radiative transfer.
! This excludes energy transfer from internal gravitational draining of 
! free water mass.

   do k = 1,nlev_sfcwater
      energy_per_m2(k) = energy_per_m2(k)  &
         + hxfers(k) - hxfers(k+1) + dt_leaf * rshort_s(k)
   enddo

! Compute heat resistance of top HALF of top soil layer (soil_rfactor).

   waterfrac = soil_water(nzg) / slmsts(ntext_soil(nzg))
   soilcond =        soilcond0(ntext_soil(nzg))  &
      + waterfrac * (soilcond1(ntext_soil(nzg))  &
      + waterfrac *  soilcond2(ntext_soil(nzg))  )
   soil_rfactor(nzg) = dslzo2(nzg) / soilcond

! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
! layer.  Impose minimum value on sfcwater_rfactor(1).  If minimum applies, 
! energy will be implicitly re-balanced later.  Apply heat transfer to bottom
! sfcwater layer and top soil layer.

   rfac = max(rfac_snowmin,sfcwater_rfactor(1))

   hxfergs = dt_leaf * (soil_tempk(nzg) - sfcwater_tempk(1))   &
           / (soil_rfactor(nzg) + rfac)

   energy_per_m2(1) = energy_per_m2(1) + hxfergs

   soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)

! Apply longwave radiative transfer and sfcwater-to-can_air sensible heat
! transfer to top sfcwater layer

   energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
      + dt_leaf * rlong_s - hxfersc

endif

! If no sfcwater layers exist and there are net positive contributions to 
! sfcwater from precipitation, shedding of water from vegetation, and 
! vapor flux with canopy air, create a new sfcwater layer and initialize 
! prognostic sfcwater fields to zero.

if (nlev_sfcwater == 0 .and.  &
    pcpg * (1. - vf) + wshed - wxfersc > 1.e-9) then

   nlev_sfcwater = 1

   sfcwater_mass  (1) = 0.
   sfcwater_energy(1) = 0.
   energy_per_m2  (1) = 0.
   sfcwater_depth (1) = 0.
endif

! Return if no sfcwater layers now exist

if (nlev_sfcwater < 1) return

! Sfcwater layers do exist

! Apply mass, energy, and depth contributions to top sfcwater layer from
! precipitation, shedding of water from vegetation, and vapor flux with 
! canopy air.  Get value for specific volume of sfcwater involved in vapor xfer.

specvol = .001
if (wxfersc > 0.) specvol =  &
   sfcwater_depth(nlev_sfcwater) / sfcwater_mass(nlev_sfcwater)

sfcwater_mass(nlev_sfcwater) = sfcwater_mass(nlev_sfcwater)  &
   + pcpg * (1. - vf)                                        &
   + wshed                                                   &
   - wxfersc

energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
   + qpcpg * (1. - vf)                                       &
   + qwshed                                                  &
   - wxfersc * alvi 

sfcwater_depth(nlev_sfcwater) = sfcwater_depth(nlev_sfcwater)  &
   + dpcpg * (1. - vf)                                         &
   + wshed * .001                                              &
   - wxfersc * specvol

! If nlev_sfcwater = 1, check whether any sfcwater mass still remains 
! (after possibly having all evaporated into canopy).  If mass is below
! threshold value, set sfcwater quantities to zero and return.

if (nlev_sfcwater == 1 .and. sfcwater_mass(nlev_sfcwater) < 1.e-9) then

   nlev_sfcwater = 0

   sfcwater_mass(1:nzs)   = 0.
   sfcwater_energy(1:nzs) = 0.
   sfcwater_depth(1:nzs)  = 0.

   return
   
endif

! Prepare to transfer water downward through snow layers by percolation.
! Fracliq is the fraction of liquid in the snowcover or surface water.
! wfree [kg/m^2] is the quantity of that liquid that is free (not attached to
! snowcover) and therefore available to drain into the layer below.

! First, prepare to sum sfcwater mass over all existing layers

totsnow = 0

! Loop downward through all existing sfcwater layers beginning with top layer

do k = nlev_sfcwater,1,-1

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

! Diagnose sfcwater temperature and liquid fraction now that new mass and
! energy contributions have been applied.  Use qwtk instead of qtk in case
! sfcwater_mass(k) is too small; "dryhcap" = 100 is very small value    

   call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
             sfcwater_tempk(k),sfcwater_fracliq(k))

! If this is bottom layer, diagnose sfcwater_rfactor.  Since sfcwater_tempk(k)
! may not be a stable computation at this point, assume that tempk = 0, which
! gives the minimum rfactor for a given density.

   if (k == 1) then
      tempk = 273.15

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

   endif

! If this is bottom layer and sfcwater rfactor is less than minimum stable 
! value, bring bottom sfcwater and top soil layer into thermal equilibrium 
! by exchanging heat between them.

   if (k == 1 .and.   &
       (sfcwater_mass(1) < snowmin_expl .or.  &
        sfcwater_rfactor(1) < rfac_snowmin)) then

! Combined sfcwater and soil water mass per square meter

      wt = sfcwater_mass(1) + soil_water(nzg) * 1.e3 * dslz(nzg)

! Combined sfcwater and soil energy per square meter.

      qwt = energy_per_m2(1) + soil_energy(nzg) * dslz(nzg)

! Soil heat capacity per square meter

      soilhcap = slcpd(ntext_soil(nzg)) * dslz(nzg)

! Diagnose equilibrium temperature and fractional liquid/ice water phases

      call qwtk(qwt,wt,soilhcap,tempk,fracliq)

! Diagnose new energy value for sfcwater based on qwt value.

      if (qwt < 0.) then

! Case of equilibrium temperature below 0 deg C.  Sfcwater_fracliq = 0.

         sfcwater_fracliq(1) = 0.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * cice * (tempk - 273.15)

      elseif (qwt > wt * alli) then

! Case of equilibrium temperature above 0 deg C.  Sfcwater fracliq = 1.

         sfcwater_fracliq(1) = 1.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * (cliq * (tempk - 273.15) + alli)

      else

! Equilibrium temperature is 0 deg C.  Determine separate values for
! sfcwater_fracliq(1) and soil_fracliq(nzg) using constraint that the sum
! of (mass * fracliq) over both components is (wt * fracliq).

! Lower bound on sfcwater_fracliq(1): case with soil_water(nzg) all liquid

         flmin = (fracliq * wt - soil_water(nzg) * 1.e3 * dslz(nzg))  &
               / sfcwater_mass(1)         

! Upper bound on sfcwater_fracliq(1): case with soil_water(nzg) all ice

         flmax = fracliq * wt / sfcwater_mass(1)         

! New sfcwater_fracliq(1) value becomes closest value within bounds to old value.

         sfcwater_fracliq(1) = max(0.,flmin,min(1.0,flmax,sfcwater_fracliq(1)))
         sfcwater_tempk(1) = 273.15
         energy_per_m2(1) = sfcwater_mass(1) * sfcwater_fracliq(1) * alli

      endif

! New energy value for soil is combined energy minus new sfcwater energy

      soil_energy(nzg) = (qwt - energy_per_m2(1)) * dslzi(nzg)

   else

! Current sfcwater layer is either not the bottom one or is thick enough
! to not require implicit thermal balance with top soil layer.

! Diagnose sfcwater temperature and liquid fraction.  Use qwtk instead of qtk
! in case sfcwater_mass(k) is too small; "dryhcap" = 100 is neglible value    
  
      call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
                sfcwater_tempk(k),sfcwater_fracliq(k))

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

            sfcwater_mass(k)  = 0.
            energy_per_m2(k)  = 0.
            sfcwater_depth(k) = 0.
          
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

! Add remaining sfcwater mass in current layer to mass sum

   totsnow = totsnow + sfcwater_mass(k)

enddo

head1 = .001 * totsnow

! Check whether any sfcwater mass remains (after possibly having all drained into soil)

if (totsnow < 1.e-9) then

! Total sfcwater mass is very close to zero.  Set sfcwater layer count to zero,
! set sfcwater arrays to exactly zero, and return

   nlev_sfcwater = 0

   sfcwater_mass  (1:nzs) = 0.
   sfcwater_energy(1:nzs) = 0.
   sfcwater_depth (1:nzs) = 0.

   head1 = 0.
   
endif

if (iwl == iwl_print) then

   call leaf_plot(iwl,             &
                  nlev_sfcwater,   &
                  time8,           &
                  linit           = 1,               &
                  lframe          = 1,               &
                  sfcwater_mass   = sfcwater_mass,   & 
                  sfcwater_energy = sfcwater_energy, & 
                  energy_per_m2   = energy_per_m2,   & 
                  sfcwater_depth  = sfcwater_depth,  & 
                  pcpg            = pcpg,            &
                  qpcpg           = qpcpg,           &
                  dpcpg           = dpcpg,           &
                  wshed           = wshed,           &
                  qwshed          = qwshed,          &
                  hxfersc         = hxfersc,         &
                  wxfersc         = wxfersc          )

endif

return
end subroutine sfcwater

!===============================================================================

subroutine sfcwater_adjust(iwl, nlev_sfcwater, sfcwater_mass, sfcwater_energy, &
                           sfcwater_depth, energy_per_m2, rshort_s, &
                           xewl, yewl, zewl, time8)

use leaf_coms, only: nzs, dt_leaf
                       
use consts_coms, only: alvi, cice, cliq, alli, piu180
use misc_coms, only: io6

implicit none

integer, intent(in)    :: iwl              ! current land cell index
integer, intent(inout) :: nlev_sfcwater    ! # of active sfc water levels

real, intent(inout) :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy(nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth (nzs) ! surface water depth [m]
real, intent(inout) :: energy_per_m2  (nzs) ! sfcwater energy [J/m^2]
real, intent(inout) :: rshort_s       (nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(in) :: xewl         ! Earth X-coordinate of land cell 'center' [m]
real, intent(in) :: yewl         ! Earth Y-coordinate of land cell 'center' [m]
real, intent(in) :: zewl         ! Earth Z-coordinate of land cell 'center' [m]

real(kind=8), intent(in) :: time8   ! model time [s]

! Local variables

real :: mass_new   (nzs) ! mass of new set of sfcwater layers [kg/m^2]
real :: energy_new (nzs) ! energy of new set of sfcwater layers [J/kg]
real :: depth_new  (nzs) ! depth of new set of sfcwater layers [m]
real :: rshort_snew(nzs) ! s/w rad flux to new set of sfcwater layers [W/m^2]

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
real :: raxis    ! distance of land cell 'center' from Earth axis [m]

! Local parameters

real, parameter :: sndenmax = 1000.   ! max sfcwater density [kg/m^3]
real, parameter :: snowmin = 11.      ! min sfcwater layer mass with multiple layers [kg/m^2] 
real, parameter :: snowmin_expl = 10. ! min sfcwater mass for explicit heat xfer [kg/m^2]
real, parameter :: rfac_snowmin = .01 ! min sfcwater rfactor [K m^2/W]

real, save, dimension(10,10) :: thick  ! snowlayer thickness scaling factor

integer, parameter :: iwl_print = 0

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
   rshort_snew(k) = 0.

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
      rshort_snew(k) = rshort_snew(k) + wtold * rshort_s      (kold)

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
      rshort_snew(k) = rshort_snew(k) + dwtold * rshort_s      (kold)

   endif

enddo

! Now that mass, energy, depth, and s/w flux have been transferred to new layers,
! copy back to original arrays

do k = 1,nlev_new
   sfcwater_mass(k)   = mass_new(k)
   sfcwater_energy(k) = energy_new(k) / sfcwater_mass(k)
   sfcwater_depth(k)  = depth_new(k)
   rshort_s(k)        = rshort_snew(k)

! Replace sfcwater energy limiter from earlier code versions with energy
! check and message.  If message ever gets printed, investigate reasons.

   if (sfcwater_energy(k) > 6.e5 .or. sfcwater_energy(k) < -2.e5) then
      raxis = sqrt(xewl ** 2 + yewl ** 2)

      glatwl = atan2(zewl,raxis)   * piu180
      glonwl = atan2(yewl,xewl) * piu180

      write(io6,*) ' '
      write(io6,*) 'Sfcwater energy is outside allowable range. '
      write(io6,*) 'iwl,k,lat,lon = ',iwl,k,glatwl,glonwl
      write(io6,*) 'sfcwater_energy = ',sfcwater_energy(k)
      write(io6,*) 'p1',energy_new(k),sfcwater_mass(k),nlev_sfcwater,nlev_new
      write(io6,*) 'p2',kold, energy_per_m2(kold)
      stop 'stop sfcwater energy'
   endif

enddo

nlev_sfcwater = nlev_new

return
end subroutine sfcwater_adjust

!===============================================================================

subroutine remove_runoff(iwl, ksn, sfcwater_fracliq, sfcwater_mass,   &
     sfcwater_tempk, sfcwater_energy, sfcwater_depth, runoff, qrunoff)

  use leaf_coms, only: nzs, dt_leaf
  use ed_options, only: runoff_time
  use consts_coms, only: alli, cliq
  use misc_coms, only: io6

  implicit none

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

End Module leaf3_sfcwater
