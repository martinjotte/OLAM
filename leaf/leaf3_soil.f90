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
Module leaf3_soil

Contains

subroutine soil(iwl, leaf_class, nlev_sfcwater, ntext_soil, ktrans,      &
                soil_tempk, soil_fracliq, soil_rfactor,                  &
                hxfergc, wxfergc, rshort_g, rlong_g, transp,             &
                soil_water, soil_energy,                                 &
                hxferg, wxfer, qwxfer,                                   &
                psi, lsl, head, head0, head1, wfree1, qwfree1, dwfree1,  &
                sfcwater_mass, energy_per_m2, sfcwater_depth,            &
                ed_transp, ed_patch                                      )

use leaf_coms, only: nzg, nzs, dslz, dslzi, slzt, dslzo2, dt_leaf,  &
                     slcons1, soilcp, slbs, slpots, slmsts, kroot, slcpd,  &
                     soilcond0, soilcond1, soilcond2

use consts_coms, only: cliq1000, cice1000, alli1000, alvi
use misc_coms,   only: io6, time8
use massflux,    only: tridiffo
use leaf3_plot,  only: leaf_plot

use ed_structure_defs

implicit none

integer, intent(in) :: iwl              ! current land cell number 
integer, intent(in) :: leaf_class       ! leaf class
integer, intent(in) :: nlev_sfcwater    ! # active levels of surface water
integer, intent(in) :: ntext_soil (nzg) ! soil textural class

integer, intent(in) :: ktrans  ! k index of soil layer supplying transpiration

real, intent(in) :: soil_tempk   (nzg) ! soil temperature (K)
real, intent(in) :: soil_fracliq (nzg) ! fraction of soil water that is liquid
real, intent(out) :: soil_rfactor (nzg) ! soil thermal resistance
real, intent(in) :: hxfergc            ! heat xfer from ground to canopy [kg/m^2]
real, intent(in) :: wxfergc            ! water xfer from ground to canopy [kg/m^2]
real, intent(in) :: rshort_g           ! s/w radiative flux abs by ground [W/m^2]
real, intent(in) :: rlong_g            ! l/w radiative flux abs by ground [W/m^2]
real, intent(in) :: transp             ! transpiration loss [kg/m^2]
real, intent(inout) :: head0           ! LBC total hydraulic head [m]
real, intent(inout) :: head1           ! UBC total hydraulic head [m]
real, intent(inout) :: wfree1          ! UBC free water mass [kg/m^2]
real, intent(inout) :: qwfree1         ! UBC free water energy [J/m^2]
real, intent(inout) :: dwfree1         ! UBC free water depth [m]

real, intent(inout) :: sfcwater_mass (nzs)
real, intent(inout) :: energy_per_m2 (nzs)
real, intent(inout) :: sfcwater_depth(nzs)

real, intent(inout) :: soil_water  (nzg) ! soil water content [vol_water/vol_tot]
real, intent(inout) :: soil_energy (nzg) ! soil energy [J/m^3]

real, intent(out) :: hxferg (nzg+1) ! soil internal heat xfer (J/m^2]
real, intent(out) :: wxfer  (nzg+1) ! soil water xfer [m]
real, intent(out) :: qwxfer (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
real, intent(out) :: psi    (nzg)   ! soil water potential [m]
real, intent(out) :: head   (nzg)   ! total hydraulic head [m]

type(patch), target, optional :: ed_patch

integer, intent(in) :: lsl

! Local variables

real :: soil_wat_avail(nzg) ! liq water in soil layer available for transport [m]
real :: water_frac    (nzg) ! Fractional water content in layer
real :: water_fraci   (nzg) ! Inverse fractional water content in layer
real :: hydresist_bot (nzg)
real :: hydresist_top (nzg)

! Defined at bottom face of soil layers & sfcwater(1)

real :: water_frac_bnd(nzg+1)
real :: xfercoef      (nzg+1)

! Defined at soil layer centers

real :: headp   (nzg)
real :: headp_ph(nzg)

real :: vctr5   (nzg+1)
real :: vctr6   (nzg+1)
real :: vctr7   (nzg+1)
real :: vctr8   (nzg+1)

integer :: k     ! vertical index over soil layers
integer :: nts   ! soil textural class
integer :: linit, lframe

real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
real :: runoff   ! runoff loss [kg/m^2]
real :: total_water
real :: total_energy
real :: soil_tempc
real :: soilcond  ! soil thermal conductivity [W/(K m)]

real, parameter :: water_frac_ph0 = .99
real, parameter :: water_def_ph0 = 1. - water_frac_ph0
   
real, dimension(nzg) :: ed_transp

integer, parameter :: iwl_print = 0

! Remove transpiration water from ktrans soil layer
! Units of wloss are [vol_water/vol_tot], of transp are [kg/m^2].

if ((.not. present(ed_patch)) .and. ktrans > 0) then

   wloss = transp * dslzi(ktrans) * 1.e-3

   soil_water(ktrans) = soil_water(ktrans) - wloss

   qwloss = wloss * (cliq1000 * (soil_tempk(ktrans) - 273.15) + alli1000)

   soil_energy(ktrans) = soil_energy(ktrans) - qwloss

elseif (present(ed_patch)) then

   do k = lsl, nzg
      wloss = ed_transp(k) * dslzi(k) * 1.e-3
      qwloss = wloss * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000)
      soil_water(k) = soil_water(k) - wloss
      soil_energy(k) = soil_energy(k) - qwloss
      ed_patch%omean_latflux =   &
           ed_patch%omean_latflux + qwloss * dslz(k) / dt_leaf
   enddo

endif

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Fractional water content and its inverse in middle of soil layer

   water_frac(k)  = soil_water(k) / slmsts(nts)
   water_fraci(k) = 1.0 / water_frac(k)

! Fractional water content and its inverse at W levels

   if (k == 1) then
      water_frac_bnd(k) = water_frac(k)
   else
      water_frac_bnd(k) = .5 * (water_frac(k-1) + water_frac(k))
   endif

! Soil heat resistance times HALF layer depth (soil_rfactor).

   soilcond =        soilcond0(ntext_soil(k))  &
      + water_frac(k) * (soilcond1(ntext_soil(k))  &
      + water_frac(k) *  soilcond2(ntext_soil(k))  )

   soil_rfactor(k) = dslzo2(k) / soilcond

enddo

water_frac_bnd(nzg+1) = water_frac(nzg)

! Soil bottom, top, and internal sensible heat xfers [J/m2]

hxferg(1) = 0.
hxferg(nzg+1) = 0.

do k = 2,nzg
   hxferg(k) = dt_leaf * (soil_tempk(k-1) - soil_tempk(k))   &
             / (soil_rfactor(k-1) + soil_rfactor(k))      
enddo

! Update soil Q values [J/m3] from internal sensible heat and upward water 
! vapor (latent heat) xfers, and from longwave and shortwave fluxes at the
! top of the soil.  This excludes effects of dew/frost formation, 
! precipitation, shedding, and percolation, which were already applied 
! to the top soil layer in subroutine sfcwater.  Update top soil moisture 
! from evaporation only if sfcwater was absent.

do k = 1,nzg
   soil_energy(k) = soil_energy(k) + dslzi(k) * (hxferg(k) - hxferg(k+1))
enddo

soil_energy(nzg) = soil_energy(nzg)  &
   + dslzi(nzg) * (dt_leaf * (rshort_g + rlong_g) - hxfergc - wxfergc * alvi)

if (nlev_sfcwater == 0) then
   soil_water(nzg) = soil_water(nzg) - 1.e-3 * wxfergc * dslzi(nzg)
endif

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

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Hydraulic resistance of top and bottom half of each soil layer

   hydresist_bot(k) = dslzo2(k) / &
      (slcons1(k,nts) * water_frac_bnd(k)   ** (2. * slbs(nts) + 3.))
      
   hydresist_top(k) = dslzo2(k) / &
      (slcons1(k,nts) * water_frac_bnd(k+1) ** (2. * slbs(nts) + 3.))
      
! Molecular water potential in middle of soil layer [m]

   psi(k) = slpots(nts) * water_fraci(k) ** slbs(nts)

! [5/14/09] Compute total hydraulic head and its derivative with respect to 
! water content in each soil layer.  This is used for (1) solving water fluxes
! implicitly between layers and (2) utilizing utilize new soil bottom boundary
! condition of specified total head.  

   head(k) = psi(k) + slzt(k)
   headp(k) = -psi(k) * slbs(nts) * water_fraci(k)
   
! If soil layer is nearly full, compute additional head that simulates 
! pressure contribution

   if (water_frac(k) > water_frac_ph0) then

! linear form      

     headp_ph(k) = (10. - slpots(nts)) / (water_def_ph0 * slmsts(nts))
     head(k)  = head(k)  &
              + headp_ph(k) * (water_frac(k) - water_frac_ph0) * slmsts(nts)

! quadratic form      

!      headp_ph(k) = (10. - slpots(nts))  &
!                  * 2. * (water_frac(k) - water_frac_ph0) / water_def_ph0 ** 2

!      head(k) = head(k) + headp_ph(k) * .5 * (water_frac(k) - water_frac_ph0)

! headp_ph addition to headp

      headp(k) = headp(k) + headp_ph(k)
   endif

! Soil water available for transport in middle of soil layer (lesser of
! excess soil water above minimum value and unfrozen portion of water)

   soil_wat_avail(k) = dslz(k)  &
      * min(soil_water(k) - soilcp(nts) , soil_water(k) * soil_fracliq(k))

enddo

! Set up quantities for implicit solution

do k = 2,nzg

   xfercoef(k) = dt_leaf                                     &
               * .5 * (soil_fracliq(k) + soil_fracliq(k-1))  &
               / (hydresist_top(k-1) + hydresist_bot(k))

   vctr5(k) = -xfercoef(k) * headp(k-1) * dslzi(k-1)
   vctr7(k) = -xfercoef(k) * headp(k  ) * dslzi(k  )
   vctr6(k) = 1. - vctr5(k) - vctr7(k)
   vctr8(k) = xfercoef(k) * (head(k-1) - head(k))

enddo

! Hydraulic head (head0) is specified as constant in time at bottom of 
! soil level k = 1.

xfercoef(1) = dt_leaf          &
            * soil_fracliq(1)  &
            / hydresist_bot(1)

vctr5(1) = 0.
vctr7(1) = -xfercoef(1) * headp(1) * dslzi(1)
vctr6(1) = 1. - vctr7(1)
vctr8(1) = xfercoef(1) * (head0 - head(1))

! For first trial, assume that hydraulic head (head1) at top of top soil layer
! will change due to implicit water flux (could be up or down).

xfercoef(nzg+1) = dt_leaf          &
                * soil_fracliq(nzg)  &
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
! soil energy at nzg, and apply to sfcwater

   qwxfer(nzg+1) = wxfer(nzg+1) * (cliq1000 * (soil_tempk(nzg) - 273.15) + alli1000)

else

! wxfer(nzg+1) was computed as negative in tridiffo

   if (wxfer(nzg+1) > -.99e-3 * wfree1) then

! If wxfer(nzg+1) is negative but does not deplete more than 99% of free 
! sfcwater in level 1 (wfree1), compute qwxfer(nzg+1) based on sfcwater
! energy at k = 1.  Apply full magnitude of wxfer(nzg+1).  
! (Factor of 1.e-3 converts from kg/m^2 to m depth.)

      qwxfer(nzg+1) = wxfer(nzg+1) * (1000. * qwfree1 / wfree1)

   else

! If wxfer(nzg+1) is negative and depletes more than 99% of free 
! sfcwater in level 1 (wfree1), compute qwxfer(nzg+1) based on sfcwater
! energy at k = 1.  Limit wxfer(nzg+1) to free water content (wfree1).  
! (Factor of 1.e-3 converts from kg/m^2 to m depth.)

      wxfer(nzg+1) = -.001 * wfree1
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

! Compute q transfers between soil layers (qwxfer) [J/m2]

   if (wxfer(k) > 0. .and. k > 1) then
      qwxfer(k) = wxfer(k) * (cliq1000 * (soil_tempk(k-1) - 273.15) + alli1000)
   else
      qwxfer(k) = wxfer(k) * (cliq1000 * (soil_tempk(k) - 273.15) + alli1000)
   endif

enddo

! Check sign of wxfer(nzg+1) and make transfers with sfcwater

if (iwl == iwl_print)                  &
   call leaf_plot(iwl,            &
                  nlev_sfcwater,  &
                  time8,          &
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
                  hxfergc        = hxfergc,      &
                  wxfergc        = wxfergc,      &
                  head0          = head0,        &
                  head1          = head1,        &
                  head           = head          )

! Apply fluxes at nzg + 1 to sfcwater at k = 1

if (wxfer(nzg+1) > -.00099 * sfcwater_mass(1)) then

! If less than 99% of sfcwater mass is removed by flux

   sfcwater_mass(1)  = sfcwater_mass(1)  + 1000. * wxfer(nzg+1) 
   energy_per_m2(1)  = energy_per_m2(1)  + qwxfer(nzg+1)
   sfcwater_depth(1) = sfcwater_depth(1) + wxfer(nzg+1)

else

! If all sfcwater mass is removed by flux

   sfcwater_mass(1)  = 0. 
   energy_per_m2(1)  = 0.
   sfcwater_depth(1) = 0.

endif

! Loop over soil T levels

do k = 1,nzg
   nts = ntext_soil(k)

! Update soil water (impose minimum value of soilcp) and energy.

   soil_water(k) = max(soilcp(nts),soil_water(k)  &
      + dslzi(k) * (wxfer(k) - wxfer(k+1)))
   soil_energy(k) = soil_energy(k) + dslzi(k) * (qwxfer(k) - qwxfer(k+1))

enddo

! Compute soil respiration if ED is being run

if (present(ed_patch)) then
   call soil_respiration(ed_patch)
endif

return
end subroutine soil

End Module leaf3_soil
