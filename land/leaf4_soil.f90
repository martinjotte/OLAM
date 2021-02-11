!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the ame team working at other institutions (University of
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

  subroutine soil(iland, iwsfc, ktrans, transp, wfree1, qwfree1,         &
                  leaf_class, glatw, glonw, head1, head0, nlev_sfcwater, &
                  gpp, sfcwater_mass, sfcwater_energy, sfcwater_depth,   &
                  energy_per_m2, soil_water, soil_energy,                &
                  wresid_vg, wsat_vg, ksat_vg, lambda_vg, en_vg,alpha_vg,&
                  head_press, soil_watfrac, head_slope, head,            &
                  soil_tempk, soil_fracliq, thermcond_soil               )

  use leaf_coms, only: nzs, dt_leaf, z_root, kroot
  use mem_sfcg,  only: itab_wsfc
  use mem_land,  only: nzg, slz, dslz, dslzi, slzt, dslzo2, kperc

  use consts_coms, only: cliq1000, cice1000, alli1000, alvi
  use misc_coms,   only: io6
  use tridiag,     only: tridiffo
  use leaf4_plot,  only: leaf_plot
  use mem_flux_accum, only: wxferi_accum, wxferp_accum, wxfer1_accum
  use oname_coms,  only: nl

  implicit none

  integer, intent(in)    :: iland           ! current land cell number 
  integer, intent(in)    :: iwsfc           ! index of current SFC grid cell
  integer, intent(in)    :: ktrans          ! k index of soil layer supplying transpiration
  real,    intent(in)    :: transp          ! transpiration loss [kg/m^2]
  real,    intent(inout) :: wfree1          ! Soil top boundary free water mass [kg/m^2]
  real,    intent(inout) :: qwfree1         ! Soil top boundary free water energy [J/m^2]
  integer, intent(in)    :: leaf_class      ! leaf class (vegetation class)
  real,    intent(in)    :: glatw           ! Latitude of land cell 'center' [deg]
  real,    intent(in)    :: glonw           ! Longitude of land cell 'center' [deg]
  real,    intent(inout) :: head1           ! Upper boundary hydraulic head [m]
  real,    intent(inout) :: head0           ! Lower boundary hydraulic head [m]
  integer, intent(inout) :: nlev_sfcwater   ! # active levels of surface water
  real,    intent(in)    :: gpp             ! Carbon gross primary production [?]
  real,    intent(inout) :: sfcwater_mass   ! surface water mass(1) [kg/m^2]
  real,    intent(inout) :: sfcwater_energy ! surface water energy(1) [J/kg]
  real,    intent(inout) :: sfcwater_depth  ! surface water depth(1) [m]
  real,    intent(inout) :: energy_per_m2   ! sfcwater energy(1) [J/m^2]

  real, intent(inout) :: soil_water    (nzg) ! soil water content [vol_water/vol_tot]
  real, intent(inout) :: soil_energy   (nzg) ! soil energy [J/m^3]
  real, intent(in   ) :: wresid_vg     (nzg) ! residual water content (vG) []
  real, intent(in   ) :: wsat_vg       (nzg) ! saturation water content (porosity) []
  real, intent(in   ) :: ksat_vg       (nzg) ! saturation hydraulic conductivity [m/s]
  real, intent(in   ) :: lambda_vg     (nzg) ! van Genuchten lambda parameter [ ]
  real, intent(in   ) :: en_vg         (nzg) ! van Genuchten n parameter [ ]
  real, intent(in   ) :: alpha_vg      (nzg) ! van Genuchten alpha parameter [ ]
  real, intent(in   ) :: head_press    (nzg) ! head from positive (confinement) pressure [m]
  real, intent(in   ) :: soil_watfrac  (nzg) ! fractional water content within vG range [limited to (0,1)]
  real, intent(in   ) :: head_slope    (nzg) ! derivative of head with respect to soil_water [m]
  real, intent(inout) :: head          (nzg) ! hydraulic head [m] (relative to local topo datum)
  real, intent(in   ) :: soil_tempk    (nzg) ! soil temperature (K)
  real, intent(in   ) :: soil_fracliq  (nzg) ! fraction of soil water that is liquid
  real, intent(in   ) :: thermcond_soil(nzg) ! soil thermal conductivity [W/(K m)]

  ! Local variables

  real :: hxferg(nzg+1) ! soil internal heat xfer (J/m^2]
  real :: wxfer (nzg+1) ! soil water xfer [m]
  real :: qwxfer(nzg+1) ! soil energy xfer from water xfer [J/m^2] 

  real :: hydresist_bot(nzg) ! hydraulic resistance of bottom half of layer [s]
  real :: hydresist_top(nzg) ! hydraulic resistance of top half of layer [s]
  real :: soil_rfactor (nzg) ! soil thermal resistance across W level

  ! Defined at bottom face of soil layers & of sfcwater layer 1

  real :: soil_watfracw(nzg+1) ! soil_watfrac interpolated to W levels [ ]
  real :: xfercoef     (nzg+1) ! hydraulic transfer coefficient

  real :: qow(0:nzg+1)
  real :: wx

  ! Defined at soil layer centers

  real :: vctr5(nzg+1)
  real :: vctr6(nzg+1)
  real :: vctr7(nzg+1)
  real :: vctr8(nzg+1)

  integer :: k, kk, km ! vertical index over soil layers

  real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
  real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
  real :: flxlim   ! water flux limiter (prior to implicit solution)
  real :: khyd_top, khyd_bot

  integer, parameter :: iland_print = 0

  ! Loop over soil T levels

  do k = 1,nzg

     ! Soil heat resistance times HALF layer depth (soil_rfactor).

     soil_rfactor(k) = dslzo2(k) / thermcond_soil(k)

     qow(k) = cliq1000 * (soil_tempk(k) - 273.15) + alli1000

  enddo
  qow(0) = alli1000
  if (wfree1 > 0.) then
     qow(nzg+1) = 1000. * qwfree1 / wfree1
  else
     qow(nzg+1) = alli1000
  endif

  ! Remove transpiration water from ktrans soil layer
  ! Units of wloss are [vol_water/vol_tot], of transp are [kg/m^2].

  if (ktrans > 0) then
     wloss  = transp * dslzi(ktrans) * 1.e-3
     qwloss = wloss * qow(ktrans)

     soil_water (ktrans) = soil_water (ktrans) - wloss
     soil_energy(ktrans) = soil_energy(ktrans) - qwloss
     head       (ktrans) = head       (ktrans) - wloss * head_slope(ktrans)
  endif

  ! Loop over soil W levels - fractional water content and heat transfer

  do k = 2,nzg
     soil_watfracw(k) = .5 * (soil_watfrac(k-1) + soil_watfrac(k))
     hxferg(k) = dt_leaf * (soil_tempk(k-1) - soil_tempk(k)) &
               / (soil_rfactor(k-1) + soil_rfactor(k))      
  enddo

  ! At top, soil_watfracw assumes mean between nzg cell value and saturation

  soil_watfracw(1)     = soil_watfrac(1)
  soil_watfracw(nzg+1) = .5 * (1.0 + soil_watfrac(nzg))

  ! Soil bottom, top, and internal sensible heat xfers [J/m2]

  hxferg(1) = 0.
  hxferg(nzg+1) = 0.

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

  ! [2/3/2017] ParFlow takes arithmetic mean of saturation hydraulic resistance
  ! (which is the same as the harmonic mean of saturation hydraulic conductivity),
  ! but it uses upwinding of (soil_water/wsat_wg).  Maybe leaf should now do the same?

  ! [12/12/2018] Tested ParFlow upwinding method in long-timestep spin-up simulation
  ! and found it to be too unstable.

  ! Loop over soil T levels

  do k = 1, nzg

     ! Update soil_energy values [J/m3] at all levels from internal heat conduction

     soil_energy(k) = soil_energy(k) + dslzi(k) * (hxferg(k) - hxferg(k+1))

     ! Compute hydraulic conductivity for relative saturations at bottom and
     ! top of cell using central hydraulic properties of cell

     call soil_wat2khyd(soil_watfracw(k),   ksat_vg(k), lambda_vg(k), en_vg(k), &
                        khyd_bot, gpp, slzt(k), z_root(leaf_class), head(k))

     call soil_wat2khyd(soil_watfracw(k+1), ksat_vg(k), lambda_vg(k), en_vg(k), &
                        khyd_top, gpp, slzt(k), z_root(leaf_class), head(k))

     ! Compute hydraulic resistance between cells

     hydresist_bot(k) = dslzo2(k) / khyd_bot
     hydresist_top(k) = dslzo2(k) / khyd_top

  enddo

  ! Set up tridiagonal matrix quantities for implicit solution.  Soil_fracliq
  ! factor in xfercoef represents immobility of frozen water.

  do k = 2,nzg
     xfercoef(k) = dt_leaf                                    &
                 * .5 * (soil_fracliq(k) + soil_fracliq(k-1)) &
                 / (hydresist_top(k-1) + hydresist_bot(k))

     vctr5(k) = -xfercoef(k) * head_slope(k-1) * dslzi(k-1)
     vctr7(k) = -xfercoef(k) * head_slope(k  ) * dslzi(k  )
     vctr6(k) = 1. - vctr5(k) - vctr7(k)
     vctr8(k) = xfercoef(k) * (head(k-1) - head(k))
  enddo

  ! If soil model is more than 50 meters deep AND this column participates in 
  ! lateral groundwater transport, assume that there is no drainage through
  ! the lower boundary.  Set the boundary transfer coefficient to zero.

  if (slz(1) < -50. .and. itab_wsfc(iwsfc)%ivoronoi == 3) then

     xfercoef(1) = 0.

  else

     ! Otherwise, compute bottom transfer coefficient assuming that hydraulic
     ! head (head0) applies at bottom of soil level k = 1.


     xfercoef(1) = dt_leaf         &
                 * soil_fracliq(1) &
                 / hydresist_bot(1)

     ! If water table is below bottom of soil model, reduce bottom transfer
     ! coefficient by re-assigning the effective depth at which head0 applies
     ! to be that of the water table itself.

     if (head0 < slz(1)) then
        xfercoef(1) = xfercoef(1) * dslzo2(1) / (slzt(1) - head0)
     endif

  endif

  vctr5(1) = 0.
  vctr7(1) = -xfercoef(1) * head_slope(1) * dslzi(1)
  vctr6(1) = 1. - vctr7(1)
  vctr8(1) = xfercoef(1) * (head0 - head(1))

  ! For first trial, assume that hydraulic head (head1) at top of top soil layer
  ! will change due to implicit water flux (could be up or down).

  xfercoef(nzg+1) = dt_leaf           &
                  * soil_fracliq(nzg) &
                  / hydresist_top(nzg)

  vctr5(nzg+1) = -xfercoef(nzg+1) * head_slope(nzg) * dslzi(nzg)
  vctr7(nzg+1) = 0.
  vctr6(nzg+1) = 1. - vctr5(nzg+1) + xfercoef(nzg+1)
  vctr8(nzg+1) = xfercoef(nzg+1) * (head(nzg) - head1)

  ! Get first trial implicit fluxes (loop from 1 to nzg+1)

  call tridiffo(nzg+1,1,nzg+1,vctr5,vctr6,vctr7,vctr8,wxfer)

  if (wxfer(nzg+1) >= 0.) then

     ! If wxfer(nzg+1) is positive or zero, compute qwxfer(nzg+1) based on
     ! soil energy at nzg

     wx = soil_water(nzg) * soil_fracliq(nzg) * dslz(nzg)
     if (wxfer(nzg+1) < wx .or. wxfer(nzg) < 0.) then
        qwxfer(nzg+1) = wxfer(nzg+1) * qow(nzg)
     else
        ! For Courant number > 1, include qow at next upstream point
        qwxfer(nzg+1) = wx * qow(nzg) + (wxfer(nzg+1) - wx) * qow(nzg-1)
     endif

  else

     ! wxfer(nzg+1) was computed as negative in tridiffo

     if (wxfer(nzg+1) > -1.e-3 * wfree1) then

        ! If wxfer(nzg+1) is negative but does not deplete all of free 
        ! sfcwater in level 1 (wfree1), compute qwxfer(nzg+1) based on sfcwater
        ! energy at k = 1.  Apply full magnitude of wxfer(nzg+1).  
        ! (Factor of 1.e-3 converts from kg/m^2 to m depth.)

        qwxfer(nzg+1) = wxfer(nzg+1) * qow(nzg+1)

     else

        ! If wxfer(nzg+1) is negative and can deplete more than available free
        ! water (wfree1), set wxfer(nzg+1) and qwxfer(nzg+1) based on the free
        ! water limits. (Factor of 1.e-3 converts from kg/m^2 to m depth.)

        wxfer (nzg+1) = -.001 * wfree1
        qwxfer(nzg+1) = -qwfree1

        ! vctr5(nzg) and vctr6(nzg) are correct here as set in the k loop above

        vctr7(nzg) = 0.
        vctr8(nzg) = xfercoef(nzg) &
           * (head(nzg-1) - head(nzg) - head_slope(nzg) * .001 * wfree1 * dslzi(nzg))

        ! Get second-time implicit fluxes (loop from 1 to nzg)

        call tridiffo(nzg,1,nzg,vctr5,vctr6,vctr7,vctr8,wxfer)

     endif

  endif

  ! Loop over soil W levels

  do k = 1,nzg

     ! Compute q transfers between soil layers (qwxfer) [J/m2]

     if (wxfer(k) < 0.) then

        wx = soil_water(k) * soil_fracliq(k) * dslz(k)
        if (wxfer(k) > -wx .or. wxfer(k+1) > 0.) then
           qwxfer(k) = wxfer(k) * qow(k)
        else
           ! For Courant number > 1, include qow at next upstream point
           qwxfer(k) = -wx * qow(k) + (wxfer(k) + wx) * qow(k+1)
        endif

     elseif (k > 1) then

        wx = soil_water(k-1) * soil_fracliq(k-1) * dslz(k-1)
        if (wxfer(k) < wx .or. wxfer(k-1) < 0.) then
           qwxfer(k) = wxfer(k) * qow(k-1)
        else
           ! For Courant number > 1, include qow at next upstream point
           qwxfer(k) = wx * qow(k-1) + (wxfer(k) - wx) * qow(k-2)
        endif

     else       ! Upward flow through bottom boundary: Assume liquid at 0 C
        qwxfer(k) = wxfer(k) * alli1000
     endif

  enddo

  if (iland == iland_print) &
     call leaf_plot(iland,                         &
                    nlev_sfcwater,                 &
                    linit          = 1,            &
                    lframe         = 1,            &
                    ktrans         = ktrans,       &
                    soil_water     = soil_water,   &
                    soil_energy    = soil_energy,  &
                    soil_rfactor   = soil_rfactor, &
                    soil_tempk     = soil_tempk,   &
                    soil_fracliq   = soil_fracliq, &
                    hxferg         = hxferg,       &
                    wxfer          = wxfer,        &
                    qwxfer         = qwxfer,       &
                    transp         = transp,       &
                    head0          = head0,        &
                    head1          = head1,        &
                    head           = head          )

  ! Loop over soil T levels

  ! Update soil water and energy from water fluxes

  do k = 1,nzg
     soil_water (k) = soil_water (k) + dslzi(k) * (wxfer (k) - wxfer (k+1))
     soil_energy(k) = soil_energy(k) + dslzi(k) * (qwxfer(k) - qwxfer(k+1))
  enddo

  ! Add infiltration, percolation, and bottom water fluxes to accumulation arrays

  wxferi_accum(iland) = wxferi_accum(iland) + wxfer(nzg)
  wxferp_accum(iland) = wxferp_accum(iland) + wxfer(kperc)
  wxfer1_accum(iland) = wxfer1_accum(iland) + wxfer(1)

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
        write(io6,*) 'iwsfc,iland,lat,lon = ',iwsfc,iland,glatw,glonw
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

! ===============================================================================

  subroutine soil_wat2pot( k, iland, soil_water, wresid_vg, wsat_vg, alpha_vg, &
                           en_vg, psi, psi_slope          )

  implicit none

  integer, intent(in) :: k,iland
  real, intent(in   ) :: soil_water
  real, intent(in   ) :: wresid_vg
  real, intent(in   ) :: wsat_vg
  real, intent(in   ) :: alpha_vg
  real, intent(in   ) :: en_vg
  real, intent(inout) :: psi
  real, intent(inout) :: psi_slope

  real :: eni_vg
  real :: em_vg
  real :: emi_vg

  real, parameter :: psi_slope_trans = 1000.  ! Inverse specific storage (for transients)
                                               ! of single grid cell at supersaturation
  real, parameter :: psi_slope_static = 10000. ! Inverse specific storage (asymptotic)
                                               ! of single grid cell at supersaturation

  real, parameter :: wfrac_inc = 0.0001 ! Small increment of wfrac for computing
                                        ! slope within normal range of wfrac
  
  real, parameter :: psi_low2 =  -50000.
  real, parameter :: psi_low1 = -100000.
  real, parameter :: psi_inc12 = psi_low2 - psi_low1

  real :: wfrac_low2, wfrac_low1, wfrac_inc12, soil_watfrac_ul

  real :: psi0

  ! van Genuchten auxiliary variables. Option for em_vg should be identical
  ! between subroutines soil_wat2pot, soil_pot2wat, and soil_wat2khyd.

  eni_vg = 1.0 / en_vg
  em_vg  = 1.0 - eni_vg  ! Option 1
 !em_vg  = 1.0           ! Option 2
  emi_vg = 1.0 / em_vg

  ! Compute soil water potential and its derivative with respect to water content
  ! based on water content in the soil layer

  if (soil_water >= wsat_vg) then

     ! For soil water exceeding saturation value, use linear relationship
     ! assuming a constant specific storage parameter.  Use more elastic
     ! specific storage value for water potential slope

     psi = (soil_water - wsat_vg) * psi_slope_static
     psi_slope = psi_slope_trans

  else

     ! Evaluate two water fraction values near the dry limit for this soil
     ! layer and compare with actual water fraction

     wfrac_low2 = ((psi_low2 * alpha_vg)**en_vg + 1.)**(-em_vg)
     wfrac_low1 = ((psi_low1 * alpha_vg)**en_vg + 1.)**(-em_vg)

     soil_watfrac_ul = (soil_water - wresid_vg) / (wsat_vg - wresid_vg)

     if (soil_watfrac_ul < wfrac_low2) then

        ! At the lower end of the range of soil water content in the van Genuchten
        ! model, the matric potential becomes singular.  To avoid computational
        ! problems associated with the singularity and nearby steepness of the
        ! matric potential curve, the matric potential is approximated by a linear
        ! function of water content.  The constant slope is defined using two
        ! values of matric potential, psi_low2 and psi_low1, along with a
        ! corresponding value of water fraction for each.

        wfrac_inc12  = wfrac_low2 - wfrac_low1

        psi_slope = psi_inc12 / (wfrac_inc12 * (wsat_vg - wresid_vg))

        psi = psi_low1 + (soil_watfrac_ul - wfrac_low1) &
            * psi_slope * (wsat_vg - wresid_vg)

     else

        ! For water fraction within the range (wfrac_low2, wsat_vg), use standard
        ! van Genuchten potential formula

        psi = ((soil_watfrac_ul)**(-emi_vg) - 1.)**eni_vg / alpha_vg

        if (soil_watfrac_ul - wfrac_inc > wfrac_low1) then

           psi0 = ((soil_watfrac_ul - wfrac_inc)**(-emi_vg) - 1.)**eni_vg / alpha_vg
           psi_slope = (psi - psi0) / (wfrac_inc * (wsat_vg - wresid_vg))

        else

           psi_slope = (psi - psi_low1) &
                     / ((soil_watfrac_ul - wfrac_low1) * (wsat_vg - wresid_vg))
        endif

     endif

  endif

  end subroutine soil_wat2pot

!===============================================================================

  subroutine soil_pot2wat( psi, wresid_vg, wsat_vg, alpha_vg, en_vg, soil_water)

  implicit none

  real, intent(in) :: psi
  real, intent(in) :: wresid_vg
  real, intent(in) :: wsat_vg
  real, intent(in) :: alpha_vg
  real, intent(in) :: en_vg

  real, intent(inout) :: soil_water

  real :: eni_vg
  real :: em_vg
  real :: emi_vg

  real, parameter :: head_slope_static = 10000. ! Inverse specific storage (asymptotic)
                                               ! of single grid cell at supersaturation

  real, parameter :: psi_low2 =  -50000.
  real, parameter :: psi_low1 = -100000.
  real, parameter :: psi_inc12 = psi_low2 - psi_low1

  real :: wfrac_low2, wfrac_low1, wfrac_inc12, soil_watfrac_ul

  ! van Genuchten auxiliary variables. Option for em_vg should be identical
  ! between subroutines soil_wat2pot, soil_pot2wat, and soil_wat2khyd.

  eni_vg = 1.0 / en_vg
  em_vg  = 1.0 - eni_vg  ! Option 1
 !em_vg  = 1.0           ! Option 2
  emi_vg = 1.0 / em_vg

  ! Compute soil water content based on soil water potential in the soil layer

  if (psi >= 0.) then
      
     ! For water potential zero or higher, use linear relationship
     ! assuming a constant specific storage parameter.

     soil_water = wsat_vg + psi / head_slope_static

  elseif (psi <= psi_low2) then

     ! At the lower end of the range of soil water content in the van Genuchten
     ! model, the matric potential becomes singular.  To avoid computational
     ! problems associated with the singularity and nearby steepness of the
     ! matric potential curve, the matric potential is approximated by a linear
     ! function of water content.  The constant slope is defined using two
     ! values of matric potential, psi_low2 and psi_low1, along with a
     ! corresponding value of water fraction for each.

     wfrac_low2 = ((psi_low2 * alpha_vg)**en_vg + 1.)**(-em_vg)
     wfrac_low1 = ((psi_low1 * alpha_vg)**en_vg + 1.)**(-em_vg)
     wfrac_inc12  = wfrac_low2 - wfrac_low1

     ! Estimate soil_water from low-end linear water potential equation       

     soil_watfrac_ul = wfrac_low1 + (psi - psi_low1) * wfrac_inc12 / psi_inc12

     soil_water = wresid_vg + soil_watfrac_ul * (wsat_vg - wresid_vg)

  else

     ! For water potential within the range (psi_low2, 0.0), use standard
     ! van Genuchten inverse potential formula

     soil_watfrac_ul = ((psi * alpha_vg)**en_vg + 1.)**(-em_vg)

     soil_water = wresid_vg + soil_watfrac_ul * (wsat_vg - wresid_vg)

  endif  

  end subroutine soil_pot2wat

!===============================================================================

  subroutine soil_wat2khydX( soil_watfrac, ksat_vg, lambda_vg, en_vg, khyd)

  implicit none

  real, intent(in   ) :: soil_watfrac
  real, intent(in   ) :: ksat_vg
  real, intent(in   ) :: lambda_vg
  real, intent(in   ) :: en_vg
  real, intent(inout) :: khyd

  real :: eni_vg
  real :: em_vg
  real :: emi_vg

  ! van Genuchten auxiliary variables. Option for em_vg should be identical
  ! between subroutines soil_wat2pot, soil_pot2wat, and soil_wat2khyd.

  eni_vg = 1.0 / en_vg
  em_vg  = 1.0 - eni_vg  ! Option 1
 !em_vg  = 1.0           ! Option 2
  emi_vg = 1.0 / em_vg

  ! Compute hydraulic conductivity

  khyd = ksat_vg * ((soil_watfrac**lambda_vg) * (1. - (1. - soil_watfrac**emi_vg)**em_vg)**2)

  end subroutine soil_wat2khydX

!===============================================================================

  subroutine soil_wat2khyd( soil_watfrac, ksat_vg, lambda_vg, en_vg, khyd, &
                            gpp, slzt, z_root, head)

  implicit none

  real, intent(in   ) :: soil_watfrac
  real, intent(in   ) :: ksat_vg
  real, intent(in   ) :: lambda_vg
  real, intent(in   ) :: en_vg
  real, intent(inout) :: khyd

  ! The following input quantities are for soil structure effect

  real, intent(in) :: gpp
  real, intent(in) :: slzt
  real, intent(in) :: z_root
  real, intent(in) :: head

  real :: ksokt_gpp, ksokt_root, ksokt_head, ksokt

  real :: eni_vg
  real :: em_vg
  real :: emi_vg

  real, parameter :: khyd_max = 0.0001  ! Maximum allowed hydraulic conductivity [m/s]

  ! van Genuchten auxiliary variables. Option for em_vg should be identical
  ! between subroutines soil_wat2pot, soil_pot2wat, and soil_wat2khyd.

  eni_vg = 1.0 / en_vg
  em_vg  = 1.0 - eni_vg  ! Option 1
 !em_vg  = 1.0           ! Option 2
  emi_vg = 1.0 / em_vg

  ! Compute hydraulic conductivity

  khyd = ksat_vg * ((soil_watfrac**lambda_vg) * (1. - (1. - soil_watfrac**emi_vg)**em_vg)**2)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SOIL STRUCTURE EXPERIMENT: BYPASS SOIL STRUCTURE EFFECT
!  GO TO 10
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Compute and apply soil structure factor if GPP > 0, soil layer is
  ! inside root zone, and head > -0.1 m

  if (gpp > 0. .and. slzt > z_root .and. head > -0.1) then

     ksokt_gpp = min(1., .333333e-3 * gpp)  ! gpp / 3000

     ! Soil structure hydraulic conductivity factor from root depth
     ! (slzt is negative; z_root is negative or zero)

     ksokt_root = (z_root - slzt) / z_root

     ! Soil structure hydraulic conductivity factor from head

     ksokt_head = min(1., 10. * (head + 0.1))

     ! Scaling factor of 999.0 suggested by Fatichi, but it could instead
     ! range from 9.0 for sand to 99999.0 for clay

     ksokt = 1.0 + 999.0 * ksokt_gpp * ksokt_root * ksokt_head

     khyd = min(ksokt * khyd, khyd_max)

  endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END OF CODE FOR SOIL STRUCTURE EXPERIMENT: BYPASS SOIL STRUCTURE EFFECT
  10 CONTINUE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  khyd = max(khyd,1.e-15) ! to prevent division by zero outside this subroutine

  end subroutine soil_wat2khyd

!===============================================================================

  subroutine head_column(nzg, iland, slzt, &
                         soil_water, wresid_vg, wsat_vg, alpha_vg, en_vg, &
                         head_press, head, head_slope, soil_watfrac)

  ! Diagnose hydraulic head (relative to local topographic datum) and
  ! its derivative with respect to soil water, and perform any necessary
  ! temporal relaxation of head_press.

  use mem_land,   only: hptimi
  use leaf_coms,  only: dt_leaf

  implicit none

  integer, intent(in) :: nzg, iland

  real, intent(in)    :: slzt        (nzg)
  real, intent(in)    :: soil_water  (nzg)
  real, intent(in)    :: wresid_vg   (nzg)
  real, intent(in)    :: wsat_vg     (nzg)
  real, intent(in)    :: alpha_vg    (nzg)
  real, intent(in)    :: en_vg       (nzg)
  real, intent(inout) :: head_press  (nzg)
  real, intent(inout) :: head        (nzg)
  real, intent(inout) :: head_slope  (nzg)
  real, intent(inout) :: soil_watfrac(nzg)

  real, parameter :: head_slope_trans = 1000.  ! Inverse specific storage (for transients)
                                               ! of single grid cell at supersaturation

  integer :: k

  real :: psi, psi_slope
  real :: head_press0
  real :: soil_watfrac_ul
  real :: head_press_change

  do k = 1,nzg

     soil_watfrac_ul = (soil_water(k) - wresid_vg(k)) / (wsat_vg(k) - wresid_vg(k))

     soil_watfrac(k) = min(1.0,max(0.001,soil_watfrac_ul))

     call soil_wat2pot(k, iland, soil_water(k), wresid_vg(k), wsat_vg(k), &
                       alpha_vg(k), en_vg(k), psi, psi_slope)

     ! Check if pressure head is required and/or active

     if (psi > 1.e-2 .or. head_press(k) > 1.e-2) then

        ! Soil water potential exceeds zero and/or it did recently resulting in
        ! a nonzero value of land%head_press.  Adjust head_press(k) toward psi
        ! at rate defined by relaxation time scale.

        head_press_change = dt_leaf * hptimi * (psi - head_press(k))

        head_press(k) = max(0.,head_press(k) + head_press_change)
           
        ! Compute total head relative to local topographic datum

        head(k) = head_press(k) + slzt(k)
        head_slope(k) = head_slope_trans

     else

        head_press(k) = 0.
        head(k) = psi + slzt(k)
        head_slope(k) = psi_slope

     endif

  enddo

  end subroutine head_column

End Module leaf4_soil
