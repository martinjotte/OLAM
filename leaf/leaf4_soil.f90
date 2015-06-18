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
Module leaf4_soil

Contains

subroutine soil(iwl, nlev_sfcwater, ktrans, ntext_soil,            &
                soil_water, soil_energy, soil_tempk, soil_fracliq, &
                sfcwater_mass, sfcwater_energy, sfcwater_depth,    &
                energy_per_m2, psi, head, head0, head1,            &
                wfree1, qwfree1, transp, glatw, glonw              )

use leaf_coms, only: nzg, nzs, dslz, dslzi, slzt, dslzo2, dt_leaf, &
                     slreso2, soilcp, slbs, slpots, slmsts, kroot, slcpd, &
                     soilcond0, soilcond1, soilcond2, slmstsi, headp_ph, &
                     water_frac_ph0, water_def_ph0, slpott

use consts_coms, only: cliq1000, cice1000, alli1000, alvi
use misc_coms,   only: io6
use tridiag,     only: tridiffo
use leaf4_plot,  only: leaf_plot
use mem_flux_accum, only: wxfer1_l_accum

implicit none

integer, intent(in) :: iwl               ! current land cell number 
integer, intent(inout) :: nlev_sfcwater  ! # active levels of surface water
integer, intent(in) :: ktrans  ! k index of soil layer supplying transpiration

integer, intent(in) :: ntext_soil  (nzg) ! soil textural class

real, intent(inout) :: soil_water  (nzg) ! soil water content [vol_water/vol_tot]
real, intent(inout) :: soil_energy (nzg) ! soil energy [J/m^3]
real, intent(in)    :: soil_tempk  (nzg) ! soil temperature (K)
real, intent(in)    :: soil_fracliq(nzg) ! fraction of soil water that is liquid

real, intent(inout) :: sfcwater_mass   ! surface water mass(1) [kg/m^2]
real, intent(inout) :: sfcwater_energy ! surface water energy(1) [J/kg]
real, intent(inout) :: sfcwater_depth  ! surface water depth(1) [m]
real, intent(inout) :: energy_per_m2   ! sfcwater energy(1) [J/m^2]

real, intent(out) :: psi (nzg) ! soil water potential [m]
real, intent(out) :: head(nzg) ! total hydraulic head [m]

real, intent(inout) :: head0   ! LBC total hydraulic head [m]
real, intent(inout) :: head1   ! UBC total hydraulic head [m]
real, intent(inout) :: wfree1  ! UBC free water mass [kg/m^2]
real, intent(inout) :: qwfree1 ! UBC free water energy [J/m^2]
real, intent(in)    :: transp  ! transpiration loss [kg/m^2]
real, intent(in)    :: glatw   ! Latitude of land cell 'center' [deg]
real, intent(in)    :: glonw   ! Longitude of land cell 'center' [deg]

! Local variables

real :: hxferg(nzg+1) ! soil internal heat xfer (J/m^2]
real :: wxfer (nzg+1) ! soil water xfer [m]
real :: qwxfer(nzg+1) ! soil energy xfer from water xfer [J/m^2] 

real :: soil_wat_avail(nzg) ! liq water in soil layer available for transport [m]
real :: water_frac    (nzg) ! Fractional water content in layer
real :: water_fraci   (nzg) ! Inverse fractional water content in layer
real :: hydresist_bot (nzg)
real :: hydresist_top (nzg)
real :: soil_rfactor  (nzg) ! soil thermal resistance
real :: hrterm_top    (nzg)
real :: hrterm_bot    (nzg)

! Defined at bottom face of soil layers & of sfcwater layer 1

real :: water_frac_bnd(nzg+1)
real :: xfercoef      (nzg+1)

! Defined at soil layer centers

real :: headp   (nzg)  ! Derivative of soil water tension head wrt soil water 

real :: vctr5 (nzg+1)
real :: vctr6 (nzg+1)
real :: vctr7 (nzg+1)
real :: vctr8 (nzg+1)

integer :: k     ! vertical index over soil layers
integer :: nts   ! soil textural class

real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
real :: soilcond  ! soil thermal conductivity [W/(K m)]

integer, parameter :: iwl_print = 0

! Remove transpiration water from ktrans soil layer
! Units of wloss are [vol_water/vol_tot], of transp are [kg/m^2].

if (ktrans > 0) then

   wloss = transp * dslzi(ktrans) * 1.e-3

   soil_water(ktrans) = soil_water(ktrans) - wloss

   qwloss = wloss * (cliq1000 * (soil_tempk(ktrans) - 273.15) + alli1000)

   soil_energy(ktrans) = soil_energy(ktrans) - qwloss

endif

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Fractional water content and its inverse in middle of soil layer

   water_frac(k)  = soil_water(k) * slmstsi(nts)
   water_fraci(k) = 1.0 / water_frac(k)

! Fractional water content at W levels

   if (k == 1) then
      water_frac_bnd(k) = water_frac(k)
   else
      water_frac_bnd(k) = .5 * (water_frac(k-1) + water_frac(k))
   endif

! Soil heat resistance times HALF layer depth (soil_rfactor).

   soilcond =            soilcond0(ntext_soil(k)) &
      + water_frac(k) * (soilcond1(ntext_soil(k)) &
      + water_frac(k) *  soilcond2(ntext_soil(k)) )

   soil_rfactor(k) = dslzo2(k) / soilcond

enddo

! Fractional water content at top (surface) W level
! (Assume mean between cell center value and saturation)

water_frac_bnd(nzg+1) = .5 * (1.0 + water_frac(nzg))

! Soil bottom, top, and internal sensible heat xfers [J/m2]

hxferg(1) = 0.
hxferg(nzg+1) = 0.

do k = 2,nzg
   hxferg(k) = dt_leaf * (soil_tempk(k-1) - soil_tempk(k)) &
             / (soil_rfactor(k-1) + soil_rfactor(k))      
enddo

! Update soil_energy values [J/m3] at all levels from internal heat conduction

do k = 1,nzg
   soil_energy(k) = soil_energy(k) + dslzi(k) * (hxferg(k) - hxferg(k+1))
enddo

! [12/07/04] Revisit the computation of water xfer between soil layers in
! relation to the extreme nonlinearity of hydraulic conductivity with respect
! to soil moisture.  What is the best value for hydraulic conductivity (or
! resistivity) at the interface between the two layers?  The answer is
! definitely not the average of the resistivity values of the layers
! because a very dry layer would shut down xfer of water into it.  The average 
! conductivity would, on the other hand, over-estimate water xfer between layers 
! when one is wet and the other dry.  A good compromise seems to be to average
! the fractional moisture content between the two layers and to apply this
! average value in computing hydraulic resistance for the bottom half of the
! upper layer and the top half of the lower layer.  Then, add these resistances
! to obtain the total hydraulic resistance between the two layers.

! [2/22/2013] Numerical simulations with an offline 100-layer model (soiltest.f90)
! show that the simple average of soil moisture between two layers (as suggested
! above) is within a few per cent of the moisture that gives the mean
! (vertically-integrated) hydraulic conductivity between the centers of those
! two layers in an equilibrium (constant flux) situation.  The actual steady-state
! moisture midway between the layers is higher, and would therefore bias the fluxes
! toward excessively high values if used for computing hydraulic conductivity.

do k = 1, nzg
   nts = ntext_soil(k)
   hrterm_bot(k) = water_frac_bnd(k) ** (-2. * slbs(nts) - 3.)
enddo

do k = 1, nzg-1
   nts = ntext_soil(k)
   if (nts == ntext_soil(k+1)) then
      hrterm_top(k) = hrterm_bot(k+1)
   else
      hrterm_top(k) = water_frac_bnd(k+1) ** (-2. * slbs(nts) - 3.)
   endif
enddo

nts = ntext_soil(nzg)
hrterm_top(nzg) = water_frac_bnd(nzg+1) ** (-2. * slbs(nts) - 3.)

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Hydraulic resistance of top and bottom half of each soil layer

   hydresist_bot(k) = slreso2(k,nts) * hrterm_bot(k)

   hydresist_top(k) = slreso2(k,nts) * hrterm_top(k)

! Molecular water potential in middle of soil layer [m]

   psi(k) = slpots(nts) * water_fraci(k) ** slbs(nts)

! [5/14/09] Compute total hydraulic head and its derivative with respect to 
! water content in each soil layer.  This is used for (1) solving water fluxes
! implicitly between layers and (2) utilizing utilize new soil bottom boundary
! condition of specified total head.  

   head(k) = psi(k) + slzt(k)
   headp(k) = -psi(k) * slbs(nts) * water_fraci(k) * slmstsi(nts)
   
! If soil layer is nearly full, compute additional head that simulates 
! pressure contribution

   if (water_frac(k) > water_frac_ph0) then

     head(k)  = head(k)  &
              + headp_ph(nts) * (water_frac(k) - water_frac_ph0) * slmsts(nts)

! headp_ph addition to headp

      headp(k) = headp(k) + headp_ph(nts)

   endif

! Soil water available for transport in middle of soil layer (lesser of
! excess soil water above minimum value and unfrozen portion of water)

   soil_wat_avail(k) = dslz(k) &
      * min(soil_water(k) - soilcp(nts) , soil_water(k) * soil_fracliq(k))

enddo

! Set up quantities for implicit solution

do k = 2,nzg

   xfercoef(k) = dt_leaf                                    &
               * .5 * (soil_fracliq(k) + soil_fracliq(k-1)) &
               / (hydresist_top(k-1) + hydresist_bot(k))

   vctr5(k) = -xfercoef(k) * headp(k-1) * dslzi(k-1)
   vctr7(k) = -xfercoef(k) * headp(k  ) * dslzi(k  )
   vctr6(k) = 1. - vctr5(k) - vctr7(k)
   vctr8(k) = xfercoef(k) * (head(k-1) - head(k))

enddo

! Hydraulic head (head0) is specified as constant in time at bottom of 
! soil level k = 1.

xfercoef(1) = dt_leaf         &
            * soil_fracliq(1) &
            / hydresist_bot(1)

vctr5(1) = 0.
vctr7(1) = -xfercoef(1) * headp(1) * dslzi(1)
vctr6(1) = 1. - vctr7(1)
vctr8(1) = xfercoef(1) * (head0 - head(1))

! For first trial, assume that hydraulic head (head1) at top of top soil layer
! will change due to implicit water flux (could be up or down).

xfercoef(nzg+1) = dt_leaf           &
                * soil_fracliq(nzg) &
                / hydresist_top(nzg)

vctr5(nzg+1) = -xfercoef(nzg+1) * headp(nzg) * dslzi(nzg)
vctr7(nzg+1) = 0.
vctr6(nzg+1) = 1. - vctr5(nzg+1) + xfercoef(nzg+1)
vctr8(nzg+1) = xfercoef(nzg+1) * (head(nzg) - head1)

! Get first trial implicit fluxes (loop from 1 to nzg+1)

call tridiffo(nzg+1,1,nzg+1,vctr5,vctr6,vctr7,vctr8,wxfer)

! Check sign and magnitude of wxfer(nzg+1)

if (wxfer(nzg+1) >= 0.) then

! If wxfer(nzg+1) is positive or zero, compute qwxfer(nzg+1) based on
! soil energy at nzg

   qwxfer(nzg+1) = wxfer(nzg+1) * (cliq1000 * (soil_tempk(nzg) - 273.15) + alli1000)

else

! wxfer(nzg+1) was computed as negative in tridiffo

   if (wxfer(nzg+1) > -1.e-3 * wfree1) then

! If wxfer(nzg+1) is negative but does not deplete all of free 
! sfcwater in level 1 (wfree1), compute qwxfer(nzg+1) based on sfcwater
! energy at k = 1.  Apply full magnitude of wxfer(nzg+1).  
! (Factor of 1.e-3 converts from kg/m^2 to m depth.)

      qwxfer(nzg+1) = wxfer(nzg+1) * (1000. * qwfree1 / wfree1)

   else

! If wxfer(nzg+1) is negative and can deplete more than available free water
! in level 1 (wfree1), set wxfer(nzg+1) and qwxfer(nzg+1) based on the free
! water limits. (Factor of 1.e-3 converts from kg/m^2 to m depth.)

      wxfer (nzg+1) = -.001 * wfree1
      qwxfer(nzg+1) = -qwfree1

! vctr5(nzg) and vctr6(nzg) are correct for this as set in the k loop above

      vctr7(nzg) = 0.
      vctr8(nzg) = xfercoef(nzg) &
         * (head(nzg-1) - head(nzg) - headp(nzg) * .001 * wfree1 * dslzi(nzg))

! Get second-time implicit fluxes (loop from 1 to nzg)

      call tridiffo(nzg,1,nzg,vctr5,vctr6,vctr7,vctr8,wxfer)

   endif

endif

! Loop over soil W levels

do k = 1,nzg

! Limit wxfer fluxes so that none exceeds half of the water content of
! the current donor cell

   if (wxfer(k) > 0. .and. k > 1) then
      wxfer(k) = min(wxfer(k),soil_water(k-1) * dslzo2(k-1))
   elseif (wxfer(k) < 0.) then
      wxfer(k) = max(wxfer(k),-soil_water(k) * dslzo2(k))
   endif

! Compute q transfers between soil layers (qwxfer) [J/m2]

   if (wxfer(k) < 0.) then
      qwxfer(k) = wxfer(k) * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000)
   elseif (k > 1) then
      qwxfer(k) = wxfer(k) * (cliq1000 * (soil_tempk(k-1) - 273.15) + alli1000)
   else       ! Upward flow through bottom boundary: Assume liquid at 0 C
      qwxfer(k) = wxfer(k) * alli1000
   endif

enddo

if (iwl == iwl_print)             &
   call leaf_plot(iwl,            &
                  nlev_sfcwater,  &
                  linit          = 1,            &
                  lframe         = 1,            &
                  ntext_soil     = ntext_soil,   &
                  ktrans         = ktrans,       &
                  soil_water     = soil_water,   &
                  soil_energy    = soil_energy,  &
                  soil_rfactor   = soil_rfactor, &
                  soil_tempk     = soil_tempk,   &
                  soil_fracliq   = soil_fracliq, &
                  hxferg         = hxferg,       &
                  wxfer          = wxfer,        &
                  qwxfer         = qwxfer,       &
                  psi            = psi,          &
                  transp         = transp,       &
                  head0          = head0,        &
                  head1          = head1,        &
                  head           = head          )

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Update soil water (impose minimum value of soilcp) and energy.

   soil_water(k) = max(soilcp(nts),soil_water(k) &
      + dslzi(k) * (wxfer(k) - wxfer(k+1)))
   soil_energy(k) = soil_energy(k) + dslzi(k) * (qwxfer(k) - qwxfer(k+1))

enddo

! Add bottom water flux to accumulation array

wxfer1_l_accum(iwl) = wxfer1_l_accum(iwl) + wxfer(1)

! Apply fluxes at nzg + 1 to sfcwater at k = 1

sfcwater_mass  = sfcwater_mass  + wxfer (nzg+1) * 1000.
energy_per_m2  = energy_per_m2  + qwxfer(nzg+1)
sfcwater_depth = sfcwater_depth + wxfer (nzg+1)

! Check for minimum sfcwater_mass

if (sfcwater_mass > 1.e-6) then

! Sfcwater mass is above minimum

   nlev_sfcwater = max(1,nlev_sfcwater)

   sfcwater_energy = energy_per_m2 / sfcwater_mass

! Check sfcwater_energy.  If message ever gets printed, investigate reasons.

   if (sfcwater_energy > 6.e5 .or. sfcwater_energy < -2.5e5) then

      write(io6,*) ' '
      write(io6,*) 'Sfcwater energy is outside allowable range.'
      write(io6,*) 'iwl,lat,lon = ',iwl,glatw,glonw
      write(io6,*) 'sfcwater energy & mass = ',sfcwater_energy,sfcwater_mass
      write(io6,*) 'sfcwater depth & energy_per_m2 = ',sfcwater_depth,energy_per_m2
      write(io6,*) 'wxfer,qwxfer = ',wxfer(nzg+1),qwxfer(nzg+1)
      write(io6,*) 'head,head1 = ',head(nzg),head1
      write(io6,*) 'wfree1,qwfree1 = ',wfree1,qwfree1
      write(io6,*) 'soil_tempk = ',soil_tempk(nzg)
      stop 'stop sfcwater energy (subroutine soil)'
   endif

else

! If insufficient sfcwater mass remains, transfer residual mass and energy to soil

   soil_water (nzg) = soil_water (nzg) + dslzi(nzg) * sfcwater_mass * .001
   soil_energy(nzg) = soil_energy(nzg) + dslzi(nzg) * energy_per_m2

   sfcwater_mass   = 0. 
   energy_per_m2   = 0.
   sfcwater_energy = 0.
   sfcwater_depth  = 0.

   if (nzs == 1) nlev_sfcwater = 0

endif

return
end subroutine soil

End Module leaf4_soil
