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
                wfree1, qwfree1, transp, glatw, glonw, flag_vg     )

use leaf_coms, only: nzg, nzs, slz, dslz, dslzi, slzt, dslzo2, dt_leaf,     &
                     slreso2_ch, slreso2_vg, soilcp_ch, soilcp_vg, slbs_ch, &
                     kroot, slcpd, soilcond0, slcons_ch, slcons_vg,         &
                     soilcond1, soilcond2, wfrac_high1, em_vg, emi_vg

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
logical, intent(in) :: flag_vg ! flag for van Genuchten model
                               ! (instead of Clapp & Hornberger)

! Local variables

real ::   soilcpk(nzg) ! copy of soilcp_ch or soilcp_vg to soil column
real ::   slconsk(nzg) ! copy of slcons_ch or slcons_vg to soil column
real ::  slreso2k(nzg) ! copy of slreso2_ch or slreso2_vg to soil column
real ::       emk(nzg) ! copy of em_vg to soil column
real ::      emik(nzg) ! copy of emi_vg to soil column
real ::     slbsk(nzg) ! copy of slbs to soil column

real :: hxferg(nzg+1) ! soil internal heat xfer (J/m^2]
real :: wxfer (nzg+1) ! soil water xfer [m]
real :: qwxfer(nzg+1) ! soil energy xfer from water xfer [J/m^2] 

real :: soil_wat_avail(nzg) ! liq water in soil layer available for transport [m]
real :: water_frac_ul (nzg) ! Fractional water content in layer [unlimited]
real :: water_frac    (nzg) ! Fractional water content in layer [limited to (0,1)]
real :: water_fraci   (nzg) ! Inverse fractional water content in layer
real :: hydresist_bot (nzg) ! hydraulic resistance of bottom half of layer [s]
real :: hydresist_top (nzg) ! hydraulic resistance of top half of layer [s]
real :: soil_rfactor  (nzg) ! soil thermal resistance

! Defined at bottom face of soil layers & of sfcwater layer 1

real :: water_frac_bnd(nzg+1)
real :: xfercoef      (nzg+1)

! Defined at soil layer centers

real :: headp(nzg)  ! Derivative of soil water tension head wrt soil water 

real :: vctr5(nzg+1)
real :: vctr6(nzg+1)
real :: vctr7(nzg+1)
real :: vctr8(nzg+1)

integer :: k, kk, km ! vertical index over soil layers
integer :: nts       ! soil textural class

real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
real :: soilcond ! soil thermal conductivity [W/(K m)]
real :: flxlim   ! water flux limiter (prior to implicit solution)

real :: aspect, scalelab, xmin, xmax, xinc, ymin, ymax, yinc
real :: wfrac, wfraci

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

   ! Get water fraction; soil water potential & derivative

   call soil_wat2pot(iwl, nts, flag_vg, soil_water(k), slzt(k), &
                     water_frac_ul(k), water_frac(k), psi(k), head(k), headp(k))

   if (flag_vg) then ! van Genuchten form

       slreso2k(k) =  slreso2_vg(k,nts)
        soilcpk(k) =   soilcp_vg  (nts)
        slconsk(k) =   slcons_vg  (nts)
            emk(k) =       em_vg  (nts)
           emik(k) =      emi_vg  (nts)

   else              ! Clapp & Hornberger form

       slreso2k(k) =  slreso2_ch(k,nts)
        soilcpk(k) =   soilcp_ch  (nts)
        slconsk(k) =   slcons_ch  (nts)
          slbsk(k) =     slbs_ch  (nts)

   endif

   ! Soil heat resistance times HALF layer depth (soil_rfactor).

   soilcond =            soilcond0(nts) &
      + water_frac(k) * (soilcond1(nts) &
      + water_frac(k) *  soilcond2(nts) )

   soil_rfactor(k) = dslzo2(k) / soilcond
enddo

! Loop over soil W levels - fractional water content and heat transfer

do k = 2,nzg
   water_frac_bnd(k) = .5 * (water_frac(k-1) + water_frac(k))
   hxferg(k) = dt_leaf * (soil_tempk(k-1) - soil_tempk(k)) &
             / (soil_rfactor(k-1) + soil_rfactor(k))      
enddo

! At top, water_frac_bnd assumes mean between nzg cell value and saturation

water_frac_bnd(1)     = water_frac(1)
water_frac_bnd(nzg+1) = .5 * (1.0 + water_frac(nzg))

! Soil bottom, top, and internal sensible heat xfers [J/m2]

hxferg(1) = 0.
hxferg(nzg+1) = 0.

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

if (flag_vg) then ! van Genuchten form

   ! Loop over soil T levels

   do k = 1, nzg

      ! Hydraulic resistance of top and bottom half of each soil layer

      hydresist_bot(k) = slreso2k(k) / (sqrt(water_frac_bnd(k)) &
                       * (1. - min( (1. - water_frac_bnd(k  )**emik(k))**emk(k), 0.9999) )**2)

      hydresist_top(k) = slreso2k(k) / (sqrt(water_frac_bnd(k+1)) &
                       * (1. - min( (1. - water_frac_bnd(k+1)**emik(k))**emk(k), 0.9999) )**2)

      ! Soil water available for transport in middle of soil layer (lesser of
      ! excess soil water above minimum value and unfrozen portion of water)

      soil_wat_avail(k) = dslz(k) &
         * min(soil_water(k) - soilcpk(k) , soil_water(k) * soil_fracliq(k))

   enddo

else          ! Clapp & Hornberger form

   ! Loop over soil T levels

   do k = 1, nzg

      ! Hydraulic resistance of top and bottom half of each soil layer

      hydresist_bot(k) = slreso2k(k) * water_frac_bnd(k  ) ** (-2. * slbsk(k) - 3.)
      hydresist_top(k) = slreso2k(k) * water_frac_bnd(k+1) ** (-2. * slbsk(k) - 3.)

      ! Soil water available for transport in middle of soil layer (lesser of
      ! excess soil water above minimum value and unfrozen portion of water)

      soil_wat_avail(k) = dslz(k) &
         * min(soil_water(k) - soilcpk(k) , soil_water(k) * soil_fracliq(k))

   enddo

endif  ! van Genuchten vs. Clapp & Hornberger

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

! Compute bottom transfer coefficient assuming that hydraulic head (head0)
! applies at bottom of soil level k = 1.

xfercoef(1) = dt_leaf         &
            * soil_fracliq(1) &
            / hydresist_bot(1)

! If water table is below bottom of soil model, reduce bottom transfer
! coefficient by re-assigning the effective depth at which head0 applies
! to be that of the water table itself.

if (head0 < slz(1)) then
   xfercoef(1) = xfercoef(1) * dslzo2(1) / (slzt(1) - head0)
endif

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

! Limit estimated fluxes in case vertical gradient of water potential
! temporarily becomes large

do k = 1,nzg+1
   km = max(1,k-1)
   kk = min(nzg,k)
   flxlim = .5 * (slconsk(km) + slconsk(kk)) * dt_leaf
   if     (vctr8(k) >  flxlim) then
           vctr8(k) =  flxlim + 0.1 * (vctr8(k) - flxlim)
   elseif (vctr8(k) < -flxlim * 2.) then
           vctr8(k) = -flxlim * 2. + 0.1 * (vctr8(k) + flxlim * 2.)
   endif
enddo

! Get first trial implicit fluxes (loop from 1 to nzg+1)

call tridiffo(nzg+1,1,nzg+1,vctr5,vctr6,vctr7,vctr8,wxfer)

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

      ! If wxfer(nzg+1) is negative and can deplete more than available free
      ! water (wfree1), set wxfer(nzg+1) and qwxfer(nzg+1) based on the free
      ! water limits. (Factor of 1.e-3 converts from kg/m^2 to m depth.)

      wxfer (nzg+1) = -.001 * wfree1
      qwxfer(nzg+1) = -qwfree1

      ! vctr5(nzg) and vctr6(nzg) are correct here as set in the k loop above

      vctr7(nzg) = 0.
      vctr8(nzg) = xfercoef(nzg) &
         * (head(nzg-1) - head(nzg) - headp(nzg) * .001 * wfree1 * dslzi(nzg))

      ! Limit estimated fluxes in case vertical gradient of water potential
      ! temporarily becomes large

      flxlim = slconsk(nzg) * dt_leaf
      if     (vctr8(nzg) >  flxlim) then
              vctr8(nzg) =  flxlim + 0.1 * (vctr8(nzg) - flxlim)
      elseif (vctr8(nzg) < -flxlim * 2.) then
              vctr8(nzg) = -flxlim * 2. + 0.1 * (vctr8(nzg) + flxlim * 2.)
      endif

      ! Get second-time implicit fluxes (loop from 1 to nzg)

      call tridiffo(nzg,1,nzg,vctr5,vctr6,vctr7,vctr8,wxfer)

   endif

endif

! Loop over soil W levels

do k = 1,nzg

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

! Update soil water (impose minimum value of soilcp) and energy.

do k = 1,nzg
   soil_water(k) = max(soilcpk(k),soil_water(k) &
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

end subroutine soil

!===============================================================================

subroutine soil_wat2pot( iwl, nts, flag_vg, soil_water, slzt,        &
                         water_frac_ul, water_frac, psi, head, headp )

  use leaf_coms, only: soilcp_ch, soilcp_vg,             &
                       slmsts_ch, slmsts_vg,             &
                       slpots_ch, slbs_ch,               &
                       slpott_high1_ch, slpott_high1_vg, &
                       headp_high_ch, headp_high_vg,     &
                       alphai_vg, emi_vg, eni_vg,        &
                       slpott_low2_vg, headp_low_vg,     &
                       wfrac_high1, wfrac_low2,          &
                       slmsts_mscp_vg, slmsts_mscpi_vg,  &
                       slmsts_ch, slmstsi_ch
         

  implicit none

  integer, intent(in) :: iwl
  integer, intent(in) :: nts
  logical, intent(in) :: flag_vg

  real, intent(in) :: soil_water
  real, intent(in) :: slzt

  real, intent(inout) :: water_frac_ul
  real, intent(inout) :: water_frac
  real, intent(inout) :: psi
  real, intent(inout) :: head
  real, intent(inout) :: headp

  real :: wfiemi

  if (flag_vg) then ! van Genuchten soil model

     ! Compute "unlimited" water fractional content in soil, which allows small
     ! exceedence of upper and lower bounds

     water_frac_ul = (soil_water - soilcp_vg(nts)) & ! "unlimited"
                   * slmsts_mscpi_vg(nts)

     ! Apply upper and lower limits to water fractional content

     water_frac = max(.001,min(.999,water_frac_ul))

     ! Compute soil water potential based on limited water fraction.  This
     ! value will be used to evaluate hydraulic head only if unlimited water
     ! fraction is within these limits, and in fact within even stricter ones.

     wfiemi = water_frac**(-emi_vg(nts))
     psi = alphai_vg(nts) * (wfiemi - 1.)**eni_vg(nts)

     ! [5/14/09] Compute total hydraulic head and its derivative with 
     ! respect towater content in each soil layer.  This is used for (1)
     ! solving water fluxes implicitly between layers and (2) utilizing
     ! utilize new soil bottom boundary condition of specified total head.

     if (water_frac_ul >= wfrac_high1) then

        ! For water fraction greater than or equal to wfrac_high1, use
        ! linear head model to represent confinement pressure

        headp = headp_high_vg(nts)

        head = slpott_high1_vg(nts) &
             + headp * slmsts_mscp_vg(nts) * (water_frac_ul - wfrac_high1) &
             + slzt

     elseif (water_frac_ul <= wfrac_low2) then

        ! For water fraction less than or equal to wfrac_low2, use
        ! linear head model to represent confinement pressure

        headp = headp_low_vg(nts)

        head = slpott_low2_vg(nts) &
             + headp * slmsts_mscp_vg(nts) * (water_frac_ul - wfrac_low2) &
             + slzt

     else

        ! For water fraction other than extremely high or low, use standard
        ! water potential formula

        headp = psi * eni_vg(nts) * emi_vg(nts) * wfiemi &
              / ((wfiemi - 1.) * (soilcp_vg(nts) - soil_water))

        head = psi + slzt

     endif

  else              ! Clapp & Hornberger model

     ! Compute "unlimited" water fractional content in soil, which allows small
     ! exceedence of upper bound

     water_frac_ul = soil_water * slmstsi_ch(nts) ! "unlimited"

     ! Apply upper and lower limits to water fractional content

     water_frac = max(.001,min(.999,water_frac_ul))

     ! Compute soil water potential based on limited water fraction.  This
     ! value will be used to evaluate hydraulic head only if unlimited water
     ! fraction is within upper limit, and in fact within even a stricter one.

     psi = slpots_ch(nts) * water_frac ** (-slbs_ch(nts))

     ! [5/14/09] Compute total hydraulic head and its derivative with 
     ! respect towater content in each soil layer.  This is used for (1)
     ! solving water fluxes implicitly between layers and (2) utilizing
     ! utilize new soil bottom boundary condition of specified total head.

     if (water_frac_ul >= wfrac_high1) then

        ! For water fraction greater than or equal to wfrac_high1, use
        ! linear head model to represent confinement pressure

        headp = headp_high_ch(nts)

        head = slpott_high1_ch(nts) &
             + headp * slmsts_ch(nts) * (water_frac_ul - wfrac_high1) &
             + slzt

     else

        ! For water fraction other than extremely high or low, use standard
        ! water potential formula

        headp = -psi * slbs_ch(nts) / (water_frac * slmsts_ch(nts))

        head = psi + slzt

     endif

  endif  ! van Genuchten model vs Clapp & Hornberger model

end subroutine soil_wat2pot

!===============================================================================

subroutine soil_pot2wat( iwl, nts, flag_vg, head, slzt, &
                         water_frac_ul, soil_water      )

  use leaf_coms, only: soilcp_ch, soilcp_vg,             &
                       slmsts_ch, slmsts_vg,             &
                       slpots_ch, slbsi_ch,              &
                       slpott_high1_ch, slpott_high1_vg, &
                       headp_highi_ch, headp_highi_vg,   &
                       alpha_vg, em_vg, en_vg,           &
                       slpott_low2_vg, headp_low_vg,     &
                       wfrac_high1, wfrac_low2,          &
                       slmsts_mscp_vg, slmsts_mscpi_vg,  &
                       slmsts_ch, slmstsi_ch,            &
                       headp_lowi_vg

  implicit none

  integer, intent(in) :: iwl
  integer, intent(in) :: nts
  logical, intent(in) :: flag_vg

  real, intent(in) :: head
  real, intent(in) :: slzt

  real, intent(inout) :: water_frac_ul
  real, intent(inout) :: soil_water

  real :: head_z

  head_z = head - slzt

  if (flag_vg) then

     ! Using van Genuchten soil water model for this land cell

     if (head_z >= slpott_high1_vg(nts)) then
      
        ! If initial head exceeds slpott_high1_vg, estimate soil_water from
        ! high-end linear water potential equation       

        water_frac_ul = wfrac_high1 &
                      + (head_z - slpott_high1_vg(nts)) &
                      * headp_highi_vg(nts) * slmsts_mscpi_vg(nts)

     elseif (head_z <= slpott_low2_vg(nts)) then
      
        ! If initial head is less than slpott_low2_vg, estimate soil_water from
        ! low-end linear water potential equation       

        water_frac_ul = wfrac_low2 &
                      + (head_z - slpott_low2_vg(nts)) &
                      * headp_lowi_vg(nts) * slmsts_mscpi_vg(nts)

     else

        ! If initial head does not exceed slpott_high1_vg, estimate soil_water
        ! from nonlinear water potential equation       

        water_frac_ul = ((head_z * alpha_vg(nts))**en_vg(nts) + 1.)**(-em_vg(nts))

     endif

     soil_water = soilcp_vg(nts) + water_frac_ul * slmsts_mscp_vg(nts)

  else

     ! Using Clapp & Hornberger soil water model for this land cell

     if (head_z >= slpott_high1_ch(nts)) then
      
        ! If initial head exceeds slpott, estimate soil_water from
        ! linear water potential equation       

        water_frac_ul = wfrac_high1 &
                      + (head_z - slpott_high1_ch(nts)) &
                      * headp_highi_ch(nts) * slmstsi_ch(nts)

     else

        ! If initial head does not exceed slpott, estimate soil_water
        ! from nonlinear water potential equation       

        water_frac_ul = (slpots_ch(nts) / head_z) ** slbsi_ch(nts)

     endif

     soil_water = water_frac_ul * slmsts_ch(nts)

  endif

end subroutine soil_pot2wat

End Module leaf4_soil
