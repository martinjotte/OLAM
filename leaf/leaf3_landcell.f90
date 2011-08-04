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
Module leaf3_landcell

Contains

subroutine landcell(iwl, nlev_sfcwater, leaf_class, ntext_soil,             &
                    soil_water, soil_energy,                                &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,         &
                    rshort_s,  rshort_v, rshort_g, rshort,                  &
                    rlong_v, rlong_s, rlong_g,                              &
                    veg_height, veg_rough, veg_tai, veg_lai, veg_fracarea,  &
                    hcapveg, can_depth,                                     &
                    rhos, vels, prss, pcpg, qpcpg, dpcpg,                   &
                    sxfer_t, sxfer_r, ustar, snowfac, vf,                   &
                    surface_ssh, ground_shv, veg_water, veg_temp,           &
                    can_temp, can_shv, stom_resist, veg_ndvip, veg_ndvif,   &
                    veg_ndvic, veg_albedo, rough, lsl, head0, head1,        &
                    xewl, yewl, zewl, timefac_ndvi, time8,                  &
                    ed_patch)
                              
use leaf_coms,         only: nzg, nzs, soil_rough, snow_rough, slcpd, dt_leaf
use misc_coms,         only: io6
use ed_structure_defs, only: patch
use leaf3_canopy,      only: canopy, vegndvi
use leaf3_sfcwater,    only: sfcwater, sfcwater_adjust, remove_runoff
use leaf3_soil,        only: soil
use leaf3_plot,        only: leaf_plot

implicit none

integer, intent(in   ) :: iwl               ! current land cell index
integer, intent(inout) :: nlev_sfcwater     ! # active levels of surface water
integer, intent(in   ) :: leaf_class        ! leaf ("vegetation") class
integer, intent(in   ) :: ntext_soil  (nzg) ! soil textural class

real, intent(inout) :: soil_water     (nzg) ! soil water [vol_water/vol_tot]
real, intent(inout) :: soil_energy    (nzg) ! soil energy [J/m^3]
real, intent(inout) :: sfcwater_mass  (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy(nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth (nzs) ! surface water depth [m]
real, intent(inout) :: rshort_s       (nzs) ! s/w net rad flux to sfc water [W/m^2]

real, intent(in   ) :: rshort_v     ! s/w net rad flux to veg [W/m^2]
real, intent(in   ) :: rshort_g     ! s/w net rad flux to soil [W/m^2]
real, intent(in   ) :: rshort       ! s/w incident sfc rad flux [W/m^2]
real, intent(in   ) :: rlong_v      ! l/w net rad flux to veg [W/m^2]
real, intent(in   ) :: rlong_s      ! l/w net rad flux to sfc water [W/m^2]
real, intent(in   ) :: rlong_g      ! l/w net rad flux to soil [W/m^2]
real, intent(in   ) :: veg_height   ! veg height [m]
real, intent(inout) :: veg_rough    ! veg roughness height [m]
real, intent(inout) :: veg_tai      ! veg total area index
real, intent(inout) :: veg_lai      ! veg leaf area index
real, intent(inout) :: veg_fracarea ! veg fractional area
real, intent(in   ) :: hcapveg      ! veg heat capacity [J/(m^2 K)]
real, intent(in   ) :: can_depth    ! canopy depth for heat and vap capacity [m]
real, intent(in   ) :: timefac_ndvi ! frac of time from past to future NDVI obs
real, intent(in   ) :: rhos         ! atmosphere air density [kg/m^3]
real, intent(in   ) :: vels         ! surface wind speed [m/s]
real, intent(in   ) :: prss         ! pressure [Pa]
real, intent(inout) :: pcpg         ! new precip amount this leaf timestep [kg/m^2]
real, intent(inout) :: qpcpg        ! new precip energy this leaf timestep [J/m^2]
real, intent(inout) :: dpcpg        ! new precip depth this leaf timestep [m]
real, intent(in   ) :: sxfer_t      ! can-to-atm heat xfer this step [kg_air K/m^2]
real, intent(in   ) :: sxfer_r      ! can-to-atm vapor xfer this step [kg_vap/m^2]
real, intent(in   ) :: ustar        ! friction velocity [m/s]
real, intent(in   ) :: snowfac      ! fractional veg burial by snowcover
real, intent(in   ) :: vf           ! fractional coverage of non-buried part of veg
real, intent(inout) :: surface_ssh  ! surface saturation spec hum [kg_vap/kg_air]
real, intent(inout) :: ground_shv   ! soil vapor spec hum [kg_vap/kg_air]
real, intent(inout) :: veg_water    ! veg sfc water content [kg/m^2]
real, intent(inout) :: veg_temp     ! veg temperature [K]
real, intent(inout) :: can_temp     ! canopy air temperature [K]
real, intent(inout) :: can_shv      ! canopy air vapor spec hum [kg/kg]
real, intent(inout) :: stom_resist  ! veg stomatal resistance [s/m]
real, intent(inout) :: veg_ndvip    ! past veg ndvi (obs time)
real, intent(inout) :: veg_ndvif    ! past veg ndvi (obs time)
real, intent(inout) :: veg_ndvic    ! current veg ndvi
real, intent(inout) :: veg_albedo   ! veg albedo
real, intent(inout) :: rough        ! net land cell roughness height [m]
real, intent(inout) :: head0        ! LBC total hydraulic head [m]
real, intent(inout) :: head1        ! UBC total hydraulic head [m]
real, intent(in   ) :: xewl         ! Earth X-coordinate of land cell 'center' [m]
real, intent(in   ) :: yewl         ! Earth Y-coordinate of land cell 'center' [m]
real, intent(in   ) :: zewl         ! Earth Z-coordinate of land cell 'center' [m]

integer, intent(in) :: lsl          ! Lowest soil layer

real(kind=8), intent(in) :: time8   ! model time [s]

type(patch), target, optional :: ed_patch

! Local arrays

real :: soil_tempk      (nzg) ! soil temperature [K]
real :: soil_fracliq    (nzg) ! fraction of soil moisture in liquid phase
real :: soil_rfactor    (nzg) ! soil thermal resistivity [K m^2/W]
real :: sfcwater_tempk  (nzs) ! surface water temperature [K]
real :: sfcwater_fracliq(nzs) ! fraction of sfc water in liquid phase
real :: energy_per_m2   (nzs) ! sfcwater energy [J/m^2]

real :: hxferg        (nzg+1) ! heat xfer between soil layers [J/m^2]

real :: wxfer       (nzg+1) ! soil water xfer [m]
real :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
real :: psi         (nzg)   ! soil water potential [m]
real :: head        (nzg)   ! total hydraulic head (including pressure) [m]

! Local variables

integer :: linit, lframe
integer :: k     ! vertical index over soil layers
integer :: nlsw1 ! maximum of (1,nlev_sfcwater)

integer :: ktrans ! vertical index of soil layer supplying transpiration

real :: transp  ! transpiration xfer this LEAF timestep [kg/m^2]
real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]
real :: hxfergc ! heat xfer from ground (soil) to can_air this step [J/m^2]
real :: wxfergc ! vapor xfer from ground (soil) to can_air this step [kg_vap/m^2]
real :: hxfersc ! heat xfer from sfcwater to can_air this step [J/m^2]
real :: wxfersc ! vapor xfer from sfcwater to can_air this step [kg_vap/m^2]
real :: wshed   ! water shed from veg this timestep [kg/m^2]
real :: qwshed  ! water energy shed from veg this timestep [J/m^2]
real :: hxferca ! can_air-to-atm heat xfer this step [J/m^2]
real :: hxfervc ! veg-to-can_air heat xfer this step [J/m^2]
real :: wxfervc ! veg-to-can_air vapor xfer this step [kg_vap/m^2]
real :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
real :: rb      ! veg-to-can_air resistance [s/m]

real :: wfree1  ! sfcwater bottom layer free water mass [kg/m^2]
real :: qwfree1 ! sfcwater bottom layer free water energy [J/m^2]
real :: dwfree1 ! sfcwater bottom layer free water depth [m]

integer, parameter :: iwl_print = 0

real, external :: rhovsil ! function to compute sat spec hum over liquid or ice

!=================================================================
! TEMPORARY RUNOFF VARIABLES
real :: runoff
real :: qrunoff
!=================================================================

if(.not.present(ed_patch))then

! Diagnose soil temperature and liquid fraction

do k = 1,nzg
   call qwtk(soil_energy(k),soil_water(k)*1.e3,  &
      slcpd(ntext_soil(k)),soil_tempk(k),soil_fracliq(k))
enddo

! Diagnose surface water temperature and liquid fraction

do k = 1,nlev_sfcwater
   call qtk(sfcwater_energy(k),sfcwater_tempk(k),sfcwater_fracliq(k))
enddo

! Update vegetation TAI, LAI, fractional area, albedo, and roughness

call vegndvi(iwl,                      &
             leaf_class,   timefac_ndvi, &
             veg_height,   veg_ndvip,    &
             veg_ndvif,    veg_ndvic,    &
             veg_tai,      veg_lai,      &
             veg_fracarea, veg_albedo,   &
             veg_rough                   ) 

else
   soil_tempk(lsl:nzg) = ed_patch%soil_tempk(lsl:nzg)
   soil_fracliq(lsl:nzg) = ed_patch%soil_fracliq(lsl:nzg)
   sfcwater_tempk(lsl:nlev_sfcwater) =   &
        ed_patch%sfcwater_tempk(lsl:nlev_sfcwater)
   sfcwater_fracliq(lsl:nlev_sfcwater) =   &
        ed_patch%sfcwater_fracliq(lsl:nlev_sfcwater)
   veg_fracarea = 1.0
endif

! Compute roughness length based on vegetation and snow.

rough = max(soil_rough,veg_rough) * (1. - snowfac) + snow_rough


! Evaluate turbulent exchanges of heat and moisture between vegetation and canopy air
! and also between soil or snow surface and canopy air.  Evaluate transfer of
! precipitation moisture and heat to vegetation and shed from vegetation to surface.
! Update vegetation and canopy air temperatures resulting from these 
! plus radiative fluxes.

nlsw1 = max(nlev_sfcwater,1)

call canopy(iwl,                                          &
            nlev_sfcwater,         ntext_soil,            &
            leaf_class,            ktrans,                &
            soil_water,            soil_fracliq,          &
            soil_tempk    (nzg),   sfcwater_mass (nlsw1), &
            sfcwater_tempk(nlsw1), veg_height,            &
            veg_rough,             veg_tai,               &
            veg_lai,               hcapveg,               &
            can_depth,             rhos,                  &
            vels,                  prss,                  &
            pcpg,                  qpcpg,                 &
            rshort,                rshort_v,              &
            rlong_v,               sxfer_t,               &
            sxfer_r,               ustar,                 &
            snowfac,               vf,                    &
            surface_ssh,           ground_shv,            &
            veg_water,             veg_temp,              &
            can_temp,              can_shv,               &
            wshed,                 qwshed,                &
            transp,                stom_resist,           &
            hxfergc,               wxfergc,               &
            hxfersc,               wxfersc,               &
            hxferca,               hxfervc,               &
            wxfervc,               rdi,                   &
            rb,                    time8,                 &
            lsl,                   ed_transp,             &
            ed_patch               )

! CALL SFCWATER:
!  1. Compute soil and sfcwater heat conductivities
!  2. Compute surface heat flux for top layer of soil or sfcwater
!  3. Evaluate internal and bottom fluxes for sfcwater
!  4. Update sfcwater layer energies due to heat flux and solar radiation
!  5. Evaluate melting and percolation of liquid through sfcwater layers

call sfcwater(iwl,                                &
              nlev_sfcwater,    ntext_soil,       &
              soil_rfactor,     soil_water,       &
              soil_energy,      sfcwater_mass,    &
              sfcwater_energy,  sfcwater_depth,   &
              soil_tempk,       soil_fracliq ,    &
              sfcwater_tempk,   sfcwater_fracliq, &
              energy_per_m2,                      &
              rshort_s,         hxfersc,          &
              wxfersc,          rlong_g,          &
              rlong_s,          pcpg,             &
              qpcpg,            dpcpg,            &
              wshed,            qwshed,           &
              head1,                              &
              vf,               wfree1,           &
              qwfree1,          dwfree1,          &
              time8                               )

call soil(iwl,                            &
          leaf_class,     nlev_sfcwater,  &
          ntext_soil,     ktrans,         &
          soil_tempk,     soil_fracliq,   &
          soil_rfactor,   hxfergc,        &
          wxfergc,        rshort_g,       &
          rlong_g,        transp,         &
          soil_water,     soil_energy,    &
          hxferg,         wxfer,          &
          qwxfer,         psi,            &
          lsl,                            &
          head,           head0,          &
          head1,          wfree1,         &
          qwfree1,        dwfree1,        &
          sfcwater_mass,  energy_per_m2,  &
          sfcwater_depth,                 &
          ed_transp,       ed_patch       )

if (iwl == iwl_print)  then

   call leaf_plot(iwl,                                 &
                  nlev_sfcwater,                       &
                  time8,                               &
                  linit            = 1,                &
                  lframe           = 1,                &
                  ntext_soil       = ntext_soil,       &
                  leaf_class       = leaf_class,       &   
                  ktrans           = ktrans,           &
                  soil_water       = soil_water,       &   
                  soil_energy      = soil_energy,      &
                  soil_rfactor     = soil_rfactor,     &   
                  soil_tempk       = soil_tempk,       &
                  soil_fracliq     = soil_fracliq,     &   
                  hxferg           = hxferg,           &
                  wxfer            = wxfer,            &  
                  qwxfer           = qwxfer,           &
                  psi              = psi,              &   
                  sfcwater_mass    = sfcwater_mass,    & 
                  sfcwater_energy  = sfcwater_energy,  &
                  sfcwater_depth   = sfcwater_depth,   & 
                  sfcwater_tempk   = sfcwater_tempk,   &
                  sfcwater_fracliq = sfcwater_fracliq, & 
                  rshort_s         = rshort_s,         &
                  rshort_v         = rshort_v,         & 
                  rshort_g         = rshort_g,         &
                  rshort           = rshort,           & 
                  rlong_v          = rlong_v,          &
                  rlong_s          = rlong_s,          &
                  rlong_g          = rlong_g,          &
                  veg_height       = veg_height,       &  
                  veg_rough        = veg_rough,        &
                  veg_tai          = veg_tai,          & 
                  veg_lai          = veg_lai,          &
                  hcapveg          = hcapveg,          & 
                  can_depth        = can_depth,        &
                  rhos             = rhos,             &  
                  vels             = vels,             &
                  prss             = prss,             &
                  pcpg             = pcpg,             &
                  qpcpg            = qpcpg,            &
                  dpcpg            = dpcpg,            &
                  sxfer_t          = sxfer_t,          &
                  sxfer_r          = sxfer_r,          &
                  ustar            = ustar,            &
                  snowfac          = snowfac,          &
                  vf               = vf,               &
                  surface_ssh      = surface_ssh,      &
                  ground_shv       = ground_shv,       &
                  veg_water        = veg_water,        &
                  veg_temp         = veg_temp,         &
                  can_temp         = can_temp,         &
                  can_shv          = can_shv,          &
                  wshed            = wshed,            &
                  qwshed           = qwshed,           &
                  transp           = transp,           &
                  stom_resist      = stom_resist,      &
                  hxfergc          = hxfergc,          &
                  wxfergc          = wxfergc,          &
                  hxfersc          = hxfersc,          &
                  wxfersc          = wxfersc,          &
                  veg_fracarea     = veg_fracarea,     &
                  hxferca          = hxferca,          &
                  hxfervc          = hxfervc,          &
                  wxfervc          = wxfervc,          & 
                  rdi              = rdi,              &
                  rb               = rb,               &
                  head0            = head0,            &
                  head1            = head1,            &
                  head             = head              )
                  
endif

! Inventory sfcwater layer(s) and adjust thicknesses if more than one

call sfcwater_adjust(iwl,                                &
                     nlev_sfcwater,    sfcwater_mass,    &
                     sfcwater_energy,  sfcwater_depth,   &
                     energy_per_m2,    rshort_s,         &
                     xewl, yewl, zewl,                   &
                     time8                               )

! Compute surface and ground vap mxrat for next timestep; put into surface_ssh 
! and ground_ssv.

nlsw1 = max(nlev_sfcwater,1)

call grndvap(iwl,                                       &
             nlev_sfcwater,          ntext_soil  (nzg), &
             soil_water     (nzg),   soil_energy (nzg), &
             sfcwater_energy(nlsw1), rhos,              &
             can_shv,                ground_shv,        &
             surface_ssh                                )
             
if (present(ed_patch)) then

   do k = lsl,nzg
      call qwtk(soil_energy(k),soil_water(k)*1.e3,  &
           slcpd(ntext_soil(k)),soil_tempk(k),soil_fracliq(k))
   enddo
   
   ! Diagnose surface water temperature and liquid fraction
   
   do k = 1,nlev_sfcwater
      call qtk(sfcwater_energy(k),sfcwater_tempk(k),sfcwater_fracliq(k))
   enddo
   
   do k = lsl,nzg
      ed_patch%soil_tempk(k) = soil_tempk(k)
      ed_patch%soil_fracliq(k) = soil_fracliq(k)
   enddo
   
   do k = 1,nlev_sfcwater
      ed_patch%sfcwater_tempk(k) = sfcwater_tempk(k)
      ed_patch%sfcwater_fracliq(k) = sfcwater_fracliq(k)
   enddo

endif

!-----------------------------------------------------------------------------
! TEMPORARY UNTIL FULL LEAF-HYDRO MODEL IS COMPLETED WITH STREAM/RIVER RUNOFF:
! Simple representation of runoff

call remove_runoff(iwl,            nlev_sfcwater,    &
                   sfcwater_fracliq, sfcwater_mass,    &
                   sfcwater_tempk,   sfcwater_energy,  &
                   sfcwater_depth,   runoff,           &
                   qrunoff                             )
!-----------------------------------------------------------------------------

if(present(ed_patch))then
   ed_patch%omean_runoff = ed_patch%omean_runoff + runoff / dt_leaf
   ed_patch%omean_qrunoff = ed_patch%omean_qrunoff + qrunoff / dt_leaf
endif

! Call land patch plot routine for selected iwl values.
! Use '(iwl == 0)' to turn off this plot call.

if (iwl == 0)  &
   call leaf_plot(iwl,                                 &
                  nlev_sfcwater,                       &
                  time8,                               &
                  linit            = 1,                &
                  lframe           = 1,                &
                  ntext_soil       = ntext_soil,       &
                  leaf_class       = leaf_class,       &   
                  ktrans           = ktrans,           &
                  soil_water       = soil_water,       &   
                  soil_energy      = soil_energy,      &
                  soil_rfactor     = soil_rfactor,     &   
                  soil_tempk       = soil_tempk,       &
                  soil_fracliq     = soil_fracliq,     &   
                  hxferg           = hxferg,           &
                  wxfer            = wxfer,            &  
                  qwxfer           = qwxfer,           &
                  psi              = psi,              &   
                  sfcwater_mass    = sfcwater_mass,    & 
                  sfcwater_energy  = sfcwater_energy,  &
                  sfcwater_depth   = sfcwater_depth,   & 
                  sfcwater_tempk   = sfcwater_tempk,   &
                  sfcwater_fracliq = sfcwater_fracliq, & 
                  rshort_s         = rshort_s,         &
                  rshort_v         = rshort_v,         & 
                  rshort_g         = rshort_g,         &
                  rshort           = rshort,           & 
                  rlong_v          = rlong_v,          &
                  rlong_s          = rlong_s,          &
                  rlong_g          = rlong_g,          &
                  veg_height       = veg_height,       &  
                  veg_rough        = veg_rough,        &
                  veg_tai          = veg_tai,          & 
                  veg_lai          = veg_lai,          &
                  hcapveg          = hcapveg,          & 
                  can_depth        = can_depth,        &
                  rhos             = rhos,             &  
                  vels             = vels,             &
                  prss             = prss,             &
                  pcpg             = pcpg,             &
                  qpcpg            = qpcpg,            &
                  dpcpg            = dpcpg,            &
                  sxfer_t          = sxfer_t,          &
                  sxfer_r          = sxfer_r,          &
                  ustar            = ustar,            &
                  snowfac          = snowfac,          &
                  vf               = vf,               &
                  surface_ssh      = surface_ssh,      &
                  ground_shv       = ground_shv,       &
                  veg_water        = veg_water,        &
                  veg_temp         = veg_temp,         &
                  can_temp         = can_temp,         &
                  can_shv          = can_shv,          &
                  wshed            = wshed,            &
                  qwshed           = qwshed,           &
                  transp           = transp,           &
                  stom_resist      = stom_resist,      &
                  hxfergc          = hxfergc,          &
                  wxfergc          = wxfergc,          &
                  hxfersc          = hxfersc,          &
                  wxfersc          = wxfersc,          &
                  veg_fracarea     = veg_fracarea,     &
                  hxferca          = hxferca,          &
                  hxfervc          = hxfervc,          &
                  wxfervc          = wxfervc,          & 
                  rdi              = rdi,              &
                  rb               = rb,               &
                  head0            = head0,            &
                  head1            = head1,            &
                  head             = head              )

! At this point, any precipitation has been added to vegetation and surface
! water, so zero out precipitation variables prior to new contribution by
! surface flux subroutines

pcpg  = 0.
qpcpg = 0.
dpcpg = 0.

if(present(ed_patch))then
   if(associated(ed_patch%younger))return
endif

return
end subroutine landcell

End Module leaf3_landcell
