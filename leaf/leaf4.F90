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

subroutine leaf4()

use leaf_coms, only: nzg, nzs, mwl, iupdndvi, s1900_ndvi, indvifile, nndvifiles

use mem_leaf,  only: land, itab_wl
use misc_coms, only: io6, s1900_sim, isubdomain

use leaf4_landcell,    only: landcell
use mem_para,          only: myrank

implicit none

! Local variables

integer :: iwl
real :: timefac_ndvi  ! fraction of elapsed time from past to future NDVI obs

! Time interpolation factor for updating NDVI

timefac_ndvi = 0.

if (iupdndvi == 1 .and. nndvifiles > 1) then
   timefac_ndvi = (s1900_sim             - s1900_ndvi(indvifile-1))  &
                / (s1900_ndvi(indvifile) - s1900_ndvi(indvifile-1))
endif

! Loop over ALL LAND CELLS

!$omp parallel do
do iwl = 2,mwl

   if (land%ed_flag(iwl) == 0) then 

! Update LAND CELL

   call landcell(iwl                 ,                                 &
      land%nlev_sfcwater        (iwl), land%leaf_class          (iwl), &
      land%ntext_soil     (1:nzg,iwl), land%soil_water    (1:nzg,iwl), &
      land%soil_energy    (1:nzg,iwl), land%sfcwater_mass (1:nzs,iwl), &
      land%sfcwater_energy(1:nzs,iwl), land%sfcwater_depth(1:nzs,iwl), &
      land%rshort_v             (iwl), land%rshort_g            (iwl), &
      land%rshort               (iwl), land%rlong_v             (iwl), &
      land%rlong_s              (iwl), land%rlong_g             (iwl), &
      land%veg_height           (iwl),                                 &
      land%veg_rough            (iwl), land%veg_tai             (iwl), &
      land%veg_lai              (iwl), land%veg_fracarea        (iwl), &
      land%hcapveg              (iwl), land%can_depth           (iwl), &
      land%rhos                 (iwl), land%vels                (iwl), &
      land%prss                 (iwl), land%pcpg                (iwl), &
      land%qpcpg                (iwl), land%dpcpg               (iwl), &
      land%sxfer_t              (iwl), land%sxfer_r             (iwl), &
      land%ustar                (iwl), land%snowfac             (iwl), &
      land%vf                   (iwl), land%surface_ssh         (iwl), &
      land%ground_shv           (iwl), land%veg_water           (iwl), &
      land%veg_temp             (iwl), land%cantemp             (iwl), &
      land%canshv               (iwl), land%stom_resist         (iwl), &
      land%veg_ndvip            (iwl), land%veg_ndvif           (iwl), &
      land%veg_ndvic            (iwl), land%veg_albedo          (iwl), &
      land%rough                (iwl), land%ggaer               (iwl), &
      land%head0                (iwl), land%head1               (iwl), &
      land%glatw                (iwl), land%glonw               (iwl), &
      land%flag_vg              (iwl), timefac_ndvi                    )

   elseif(land%ed_flag(iwl) == 1)then

#ifdef USE_ED2
      ! Time to call ED.
      call ed_biophys_wrapper()
#endif

   endif

! Zero out LAND%SXFER_T(iwl) and LAND%SXFER_R(iwl) now that they have 
! been applied to the canopy

   land%sxfer_t(iwl) = 0.
   land%sxfer_r(iwl) = 0.
   land%sxfer_c(iwl) = 0.

enddo
!$omp end parallel do

return
end subroutine leaf4

!===============================================================================

subroutine grndvap(iwl, nlev_sfcwater, nts, soil_water, soil_energy, &
                   sfcwater_energy, rhos, canshv, surface_ssh, ground_shv, flag_vg)

use leaf_coms,   only: nstyp, slcpd, nzg
use consts_coms, only: grav, rvap
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iwl         ! current land cell number 
integer, intent(in) :: nlev_sfcwater ! # active levels of surface water
integer, intent(in) :: nts           ! soil textural class (local name)

real, intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
real, intent(in)  :: soil_energy     ! [J/m^3]
real, intent(in)  :: sfcwater_energy ! [J/kg]
real, intent(in)  :: rhos            ! air density [kg/m^3]
real, intent(in)  :: canshv          ! canopy vapor spec hum [kg_vap/kg_air]
real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]
real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
logical, intent(in) :: flag_vg

! Local variables

real :: tempk     ! surface water temp [K]
real :: fracliq   ! fraction of surface water in liquid phase
real :: can_rhov  ! canopy water vapor density [kg_vap/m^3]
real :: sfc_rhovs ! ground sfc saturation vapor density [kg_vap/m^3]
real :: gnd_rhov  ! ground sfc evaporative vapor density [kg_vap/m^3]

real, external :: rhovsil  ! function to compute sat vapor density (over ice or liq)

  can_rhov = canshv * rhos

  if (nlev_sfcwater > 0) then

  ! sfc_rhovs is the saturation vapor density of the top soil or snow surface
  ! and is used for dew formation and snow evaporation.

     call qtk(sfcwater_energy,tempk,fracliq)
     sfc_rhovs = rhovsil(tempk-273.15)
     surface_ssh = sfc_rhovs / rhos

  else

  ! Without snowcover, gnd_rhov is the effective saturation vapor density
  ! of soil and is used for soil evaporation.

     call qwtk(soil_energy,soil_water*1.e3,slcpd(nts),tempk,fracliq)
     sfc_rhovs = rhovsil(tempk-273.15)

     call grndvap_ab(iwl,nts,tempk,soil_water,can_rhov,sfc_rhovs,gnd_rhov,flag_vg)

     surface_ssh = sfc_rhovs / rhos
     ground_shv  = gnd_rhov  / rhos

  endif

end subroutine grndvap

!===============================================================================

subroutine grndvap_ab(iwl,nts,tempk,soil_water,can_rhov,sfc_rhovs,gnd_rhov,flag_vg)

use leaf_coms, only: nstyp, nzg, slzt

use consts_coms, only: grav, rvap
use misc_coms,   only: io6
 use leaf4_soil, only: soil_wat2pot

implicit none

integer, intent(in) :: iwl       ! current land cell number 
integer, intent(in) :: nts       ! soil textural class (local name)

real, intent(in)  :: tempk       ! soil temperature [K]
real, intent(in)  :: soil_water  ! soil water content [vol_water/vol_tot]
real, intent(in)  :: can_rhov    ! canopy vapor density [kg_vap/m^3]
real, intent(in)  :: sfc_rhovs   ! surface saturation vapor density [kg_vap/m^3]
real, intent(out) :: gnd_rhov    ! ground equilibrium vapor density [kg_vap/m^3]
logical, intent(in) :: flag_vg

! Local parameter

real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

! Local variables

real :: slpotvn ! soil water potential [m]
real :: alpha   ! "alpha" term in Lee and Pielke (1993)
real :: beta    ! "beta" term in Lee and Pielke (1993)
real, save :: sfldcap(nstyp)  ! soil water field capacity [vol_water/vol_tot]
real :: water_frac_ul
real :: water_frac
real :: psi
real :: head ! soil water potential [m]
real :: headp

data sfldcap/.135,.150,.195,.255,.240,.255,.322,.325,.310,.370,.367,.535/

  ! Without snowcover, gnd_rhov is the effective saturation mixing
  ! ratio of soil and is used for soil evaporation.  First, compute the
  ! "alpha" term or soil "relative humidity" and the "beta" term.

  call soil_wat2pot(iwl,nts,flag_vg,soil_water,slzt(nzg), &
                  water_frac_ul, water_frac, psi, head, headp)

  alpha = exp(gorvap * head / tempk)
  beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nts)) * 3.14159)) ** 2

  gnd_rhov = sfc_rhovs * alpha * beta + (1. - beta) * can_rhov

end subroutine grndvap_ab

!===============================================================================

subroutine sfcrad_land(iwl, leaf_class, ntext_soil, nlev_sfcwater,      &
                       sfcwater_energy, sfcwater_depth,                 &
                       soil_energy, soil_water,                         &
                       veg_temp, veg_fracarea, veg_height,              &
                       veg_albedo, rshort, rlong,                       &
                       rshort_s, rshort_g, rshort_v,                    &
                       rlong_g, rlong_s, rlong_v,                       &
                       rlongup, rlong_albedo, albedo_beam, snowfac, vf, &
                       cosz, xewl, yewl, zewl, wnxl, wnyl, wnzl, flag_vg)

use leaf_coms,   only: nzs, slcpd, slmstsi_ch, emisv, emisg, slmstsi_vg
use consts_coms, only: stefan, eradi
use misc_coms,   only: io6
use mem_radiate, only: sunx, suny, sunz

implicit none

integer, intent(in) :: iwl         ! index of current land cell
integer, intent(in) :: leaf_class    ! leaf class
integer, intent(in) :: ntext_soil    ! soil textural class
integer, intent(in) :: nlev_sfcwater ! # of active surface water layers

real, intent(in) :: sfcwater_energy(nzs) ! surface water internal energy [J/kg]
real, intent(in) :: sfcwater_depth(nzs)  ! surface water depth [m]
real, intent(in) :: soil_energy  ! soil internal energy [J/m^3]
real, intent(in) :: soil_water   ! soil water content [vol_water/vol_tot]
real, intent(in) :: veg_temp     ! veg temp [K]
real, intent(in) :: veg_fracarea ! veg fractional area coverage
real, intent(in) :: veg_height   ! veg height [m]
real, intent(in) :: veg_albedo   ! veg albedo
real, intent(in) :: rshort       ! downward surface incident s/w rad flux [W/m^2]
real, intent(in) :: rlong        ! downward surface incident l/w rad flux [W/m^2]
real, intent(in) :: xewl         ! land cell earth x coordinate [m]
real, intent(in) :: yewl         ! land cell earth y coordinate [m]
real, intent(in) :: zewl         ! land cell earth z coordinate [m]
real, intent(in) :: wnxl         ! land cell norm unit vec x comp [m]
real, intent(in) :: wnyl         ! land cell norm unit vec y comp [m]
real, intent(in) :: wnzl         ! land cell norm unit vec z comp [m]

logical, intent(in) :: flag_vg   ! van Genuchten flag

real, intent(out) :: rshort_s(nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(out) :: rshort_g ! s/w net rad flux to soil [W/m^2]
real, intent(out) :: rshort_v ! s/w net rad flux to veg [W/m^2]
real, intent(out) :: rlong_g  ! l/w net rad flux to soil [W/m^2]
real, intent(out) :: rlong_s  ! l/w net rad flux to sfc water [W/m^2]
real, intent(out) :: rlong_v  ! l/w net rad flux to veg [W/m^2]
real, intent(out) :: rlongup  ! upward sfc l/w rad flux [W/m^2]
real, intent(out) :: rlong_albedo  ! albedo for upward sfc l/w
real, intent(out) :: albedo_beam   ! land cell albedo (beam)
real, intent(out) :: snowfac  ! fractional veg burial by snowcover
real, intent(out) :: vf       ! fractional coverage of non-buried part of veg
real, intent(out) :: cosz     ! solar zenith angle for land cell

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
real :: emg       ! soil emissivity
real :: ems       ! surface water emissivity
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
             slcpd(ntext_soil),soil_tempk,soil_fracliq)

! Shortwave radiation calculations

   if (flag_vg) then
      fc50 = min(.50, soil_water * slmstsi_vg(ntext_soil))
   else
      fc50 = min(.50, soil_water * slmstsi_ch(ntext_soil))
   endif

   if (leaf_class == 2) then

      alg = .80  ! Firn/glacier albedo

   else

      if (ntext_soil == 1 .or. ntext_soil == 2) then
         alg_dry = .41  ! Higher albedo for sand or loamy sand
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
!  rshort_a = rshort * albedo_beam

! Longwave radiation calculations

   emv    = emisv(leaf_class)
   emg    = emisg(ntext_soil)
   glong  = emg * stefan * soil_tempk ** 4
   vlong  = emv * stefan * veg_temp ** 4

   rlonga_v = rlong * vf * (emv + vfc * (1. - emg))
   rlonga_g = rlong * vfc * emg
   rlongv_g = vlong * vf * emg
   rlongv_a = vlong * vf * (2. - emg - vf + emg * vf)
   rlongg_v = glong * vf * emv
   rlongg_a = glong * vfc

   rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - emg))

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
!  rshort_a = rshort * albedo_beam

! Longwave radiation calculations

   emv   = emisv(leaf_class)
   ems   = 1.0
   slong = ems * stefan * sfcwater_tempk ** 4
   vlong = emv * stefan * veg_temp ** 4

   rlonga_v = rlong * vf * (emv + vfc * (1. - ems))
   rlonga_s = rlong * vfc * ems
   rlongv_s = vlong * vf * ems
   rlongv_a = vlong * vf * (2. - ems - vf + ems * vf)
   rlongs_v = slong * vf * emv
   rlongs_a = slong * vfc

   rlong_albedo = (vf * (1. - emv) + vfc * vfc * (1. - ems))

   rlongup = rlongv_a + rlongs_a

   rlong_g = 0.
   rlong_s = rlonga_s - rlongs_a + rlongv_s - rlongs_v
   rlong_v = rlonga_v - rlongv_a + rlongs_v - rlongv_s

endif

end subroutine sfcrad_land

