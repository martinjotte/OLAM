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
subroutine leaf3()

use leaf_coms, only: nzg, nzs, mwl, iupdndvi, s1900_ndvi, indvifile, nndvifiles

use mem_leaf,  only: land, first_site, itab_wl
use misc_coms, only: io6, time8, s1900_sim, iparallel

use ed_structure_defs, only: patch,site
use ed_options,        only: ied_offline
use leaf3_landcell,    only: landcell
use mem_para,          only: myrank

!$ use omp_lib

implicit none

! Local variables

integer :: iwl
real :: timefac_ndvi  ! fraction of elapsed time from past to future NDVI obs

type(patch), pointer :: ed_patch
type(site),  pointer :: ed_site

! Time interpolation factor for updating NDVI

timefac_ndvi = 0.

if (iupdndvi == 1 .and. nndvifiles > 1) then
   timefac_ndvi = (s1900_sim             - s1900_ndvi(indvifile-1))  &
                / (s1900_ndvi(indvifile) - s1900_ndvi(indvifile-1))
endif

! Loop over ALL LAND CELLS

ed_site => first_site

!$omp parallel do
do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   if (land%ed_flag(iwl) == 0 .and. ied_offline == 0) then 

! Update LAND CELL

   call landcell(iwl                 ,                                  &
      land%nlev_sfcwater        (iwl), land%leaf_class           (iwl), &
      land%ntext_soil     (1:nzg,iwl), land%soil_water     (1:nzg,iwl), &
      land%soil_energy    (1:nzg,iwl), land%sfcwater_mass  (1:nzs,iwl), &
      land%sfcwater_energy(1:nzs,iwl), land%sfcwater_depth (1:nzs,iwl), &
      land%rshort_s       (1:nzs,iwl), land%rshort_v             (iwl), &
      land%rshort_g             (iwl), land%rshort               (iwl), &
      land%rlong_v              (iwl), land%rlong_s              (iwl), &
      land%rlong_g              (iwl), land%veg_height           (iwl), &
      land%veg_rough            (iwl), land%veg_tai              (iwl), &
      land%veg_lai              (iwl), land%veg_fracarea         (iwl), &
      land%hcapveg              (iwl), land%can_depth            (iwl), &
      land%rhos                 (iwl), land%vels                 (iwl), &
      land%prss                 (iwl), land%pcpg                 (iwl), &
      land%qpcpg                (iwl), land%dpcpg                (iwl), &
      land%sxfer_t              (iwl), land%sxfer_r              (iwl), &
      land%ustar                (iwl), land%snowfac              (iwl), &
      land%vf                   (iwl), land%surface_ssh          (iwl), &
      land%ground_shv           (iwl), land%veg_water            (iwl), &
      land%veg_temp             (iwl), land%can_temp             (iwl), &
      land%can_shv              (iwl), land%stom_resist          (iwl), &
      land%veg_ndvip            (iwl), land%veg_ndvif            (iwl), &
      land%veg_ndvic            (iwl), land%veg_albedo           (iwl), &
      land%rough                (iwl), land%lsl                  (iwl), &
      land%head0                (iwl), land%head1                (iwl), &
      land%xew                  (iwl), land%yew                  (iwl), &
      land%zew                  (iwl)                                 , &
      timefac_ndvi                   , time8                            )

   elseif(land%ed_flag(iwl) == 1)then

! in this case, ED is being run.  Therefore, we want to loop over patches 
! and send the patch variables.

      ed_site%omean_precip = ed_site%omean_precip + land%pcpg(ed_site%iland)
      ed_site%omean_qprecip = ed_site%omean_qprecip + land%qpcpg(ed_site%iland)

      ed_patch => ed_site%oldest_patch
      
      do while(associated(ed_patch))

         ed_site%omean_netrad = ed_site%omean_netrad +   &
              ed_patch%area * (  &
              (1.0 - ed_patch%albedo_beam) *   &
              (land%rshort(ed_site%iland) -   &
              land%rshort_diffuse(ed_site%iland)) +  &
              (1.0 - ed_patch%albedo_diffuse) *   &
              land%rshort_diffuse(ed_site%iland) ) +   &
              land%rlong(ed_site%iland) *  &
              (1.0 - ed_patch%rlong_albedo) -   &
              ed_patch%rlongup

         call landcell(iwl                                                      , &
            ed_patch%nlev_sfcwater             , land%leaf_class           (iwl), &
            ed_patch%ntext_soil     (1:nzg)    , ed_patch%soil_water     (1:nzg), &
            ed_patch%soil_energy    (1:nzg), &
            ed_patch%sfcwater_mass  (1:nzs), &
            ed_patch%sfcwater_energy(1:nzs)    , ed_patch%sfcwater_depth (1:nzs), &
            ed_patch%rshort_s       (1:nzs)    , land%rshort_v           (iwl)  , &
            ed_patch%rshort_g                  , land%rshort             (iwl)  , &
            land%rlong_v            (iwl)      , ed_patch%rlong_s               , &
            ed_patch%rlong_g                   , ed_patch%veg_height            , &
            ed_patch%veg_rough                 , ed_patch%lai                   , &
            land%veg_lai            (iwl)      , land%veg_fracarea       (iwl)  , &
            land%hcapveg            (iwl)      , ed_patch%can_depth             , &
            land%rhos               (iwl)      , land%vels               (iwl)  , &
            land%prss               (iwl)      , land%pcpg               (iwl)  , &
            land%qpcpg              (iwl)      , land%dpcpg              (iwl)  , &
            ed_patch%sxfer_t                   , ed_patch%sxfer_r               , &
            ed_patch%ustar                     , ed_patch%snowfac               , &
            land%vf                 (iwl)      , ed_patch%surface_ssh           , &
            ed_patch%ground_shv                , land%veg_water          (iwl)  , &
            land%veg_temp           (iwl)      , ed_patch%can_temp              , &
            ed_patch%can_shv                   , land%stom_resist        (iwl)  , &
            land%veg_ndvip          (iwl)      , land%veg_ndvif          (iwl)  , &
            land%veg_ndvic          (iwl)      , land%veg_albedo         (iwl)  , &
            ed_patch%rough                     , land%lsl                (iwl)  , &
            land%head0              (iwl)      , land%head1              (iwl)  , &
            land%xew                (iwl)      , land%yew                  (iwl), &
            land%zew                (iwl)                                       , &
            timefac_ndvi                       , time8                          , &
            ed_patch                                                              )

         ed_patch => ed_patch%younger
         
      enddo
      
   endif

! Zero out LAND%SXFER_T(iwl) and LAND%SXFER_R(iwl) now that they have 
! been applied to the canopy (but save previous values for plotting)

   land%sxfer_tsav(iwl) = land%sxfer_t(iwl)
   land%sxfer_rsav(iwl) = land%sxfer_r(iwl)

   land%sxfer_t(iwl) = 0.
   land%sxfer_r(iwl) = 0.

   if (land%ed_flag(iwl) == 1) then
      ed_patch => ed_site%oldest_patch
      do while(associated(ed_patch))
         ed_patch%sxfer_t = 0.0
         ed_patch%sxfer_r = 0.0
         ed_patch => ed_patch%younger
      enddo
      ed_site => ed_site%next_site
   endif

enddo
!$omp end parallel do

return
end subroutine leaf3

!===============================================================================

subroutine grndvap(iwl, nlev_sfcwater, nts, soil_water, soil_energy,    &
                   sfcwater_energy, rhos, can_shv, ground_shv, surface_ssh)

use leaf_coms,   only: nstyp, slbs, slmsts, slpots, slcpd, nzg
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
real, intent(in)  :: can_shv         ! canopy vapor spec hum [kg_vap/kg_air]
real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]

! Local parameter

real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant

! Local variables

real :: slpotvn ! soil water potential [m]
real :: alpha   ! "alpha" term in Lee and Pielke (1993)
real :: beta    ! "beta" term in Lee and Pielke (1993)
real :: tempk   ! surface water temp [K]
real :: fracliq ! fraction of surface water in liquid phase
real, save, dimension(nstyp) :: sfldcap  ! soil water field capacity [vol_water/vol_tot]

data sfldcap/.135,.150,.195,.255,.240,.255,.322,.325,.310,.370,.367,.535/

real, external :: rhovsil  ! function to compute sat vapor density (over ice or liq)

! surface_ssh is the saturation mixing ratio of the top soil or snow surface
! and is used for dew formation and snow evaporation.

if (nlev_sfcwater > 0) then
   call qtk(sfcwater_energy,tempk,fracliq)
   surface_ssh = rhovsil(tempk-273.15) / rhos
else

! Without snowcover, ground_shv is the effective saturation mixing
! ratio of soil and is used for soil evaporation.  First, compute the
! "alpha" term or soil "relative humidity" and the "beta" term.

   call qwtk(soil_energy,soil_water*1.e3,slcpd(nts),tempk,fracliq)
   surface_ssh = rhovsil(tempk-273.15) / rhos

   slpotvn = slpots(nts) * (slmsts(nts) / soil_water) ** slbs(nts)
   alpha = exp(gorvap * slpotvn / tempk)
   beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nts)) * 3.14159)) ** 2
   ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv

endif

return
end subroutine grndvap

!===============================================================================

subroutine sfcrad_land(iwl, leaf_class, ntext_soil, nlev_sfcwater,      &
                       sfcwater_energy, sfcwater_depth,                 &
                       soil_energy, soil_water,                         &
                       veg_temp, veg_fracarea, veg_height,              &
                       veg_albedo, rshort, rlong,                       &
                       rshort_s, rshort_g, rshort_v,                    &
                       rlong_g, rlong_s, rlong_v,                       &
                       rlongup, rlong_albedo, albedo_beam, snowfac, vf, &
                       cosz, xewl, yewl, zewl, wnxl, wnyl, wnzl         )

use leaf_coms,   only: nzs, slcpd, slmsts, emisv, emisg
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
real :: fracabs (nzs)    ! fraction of rshort that is absorbed by snowcover
real :: vfc       ! 1 - vf
real :: fcpct     ! soil water fraction
real :: alg       ! soil albedo
real :: als       ! snowcover albedo (needs better formula based on age of snow)
real :: alv       ! veg albedo
real :: rad       ! fraction s/w rad absorbed into sfcwater + soil
real :: fractrans ! fraction of s/w rad flux transmitted through snowcover layer
real :: absg      ! fraction of rshort that is absorbed by ground
real :: algs      ! albedo from snow plus ground
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
! each land cell, plus functions of snowcover, all of which are used in leaf3.

!-------------------------------------------------------------------------------
! Two forms - but still need to consider shadowing and exact energy conservation.
! A correct formulation must exchange radiative energy between model columns.

! Compute solar zenith angle for land cells

cosz = (xewl * sunx + yewl * suny + zewl * sunz) * eradi

! Compute solar incidence angle for land cells (accounts for local topo slope)

cosz = wnxl * sunx + wnyl * suny + wnzl * sunz

! COSZ IS NOT CURRENTLY USED IN THIS SUBROUTINE
!-------------------------------------------------------------------------------

if (nlev_sfcwater == 0) then

! Case with no surface water

! Shortwave radiation calculations

   snowfac = 0.
   alv = veg_albedo

   fcpct = soil_water / slmsts(ntext_soil)  ! soil water fraction

   if (leaf_class == 2) then

      alg = .50  ! Firn/glacier albedo

   else

      if (fcpct > .5) then
         alg = .14                ! ground albedo
      else
         alg = .31 - .34 * fcpct  ! ground albedo
      endif

   endif

   vf = veg_fracarea
   vfc = 1. - vf

   albedo_beam = vf * alv + vfc * vfc * alg

   rshort_g = rshort * vfc * (1. - alg)
   rshort_s(1:nzs) = 0.
   rshort_v = rshort * vf * (1. - alv + vfc * alg)
!  rshort_a = rshort * albedo_beam

! Longwave radiation calculations

! Diagnose soil temperature and liquid fraction

   call qwtk(soil_energy,soil_water*1.e3,  &
             slcpd(ntext_soil),soil_tempk,soil_fracliq)

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

   alv = veg_albedo

! Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid
! to .5 for all-ice

   als = .5 - .36 * sfcwater_fracliq
   rad = 1. - als    ! fraction shortwave absorbed into sfcwater + soil

   snowfac = 0.
   do k = nlev_sfcwater,1,-1

      snowfac = snowfac + sfcwater_depth(k)

! fractrans is fraction of shortwave entering each sfcwater layer that
! gets transmitted through that layer

      fractrans = exp(-20. * sfcwater_depth(k))

! fracabs(k) is fraction of total incident shortwave (at top of top sfcwater
! layer) that is absorbed in each sfcwater layer

      fracabs(k) = rad * (1. - fractrans)

! rad is fraction of total incident shortwave (at top of top sfcwater layer)
! that remains at bottom of current sfcwater layer

      rad = rad * fractrans
   enddo

   snowfac = snowfac / max(.001,veg_height)
   
! If sfcwater (snowcover) is deep enough to cover most vegetation, assume that 
! all is covered

   if (snowfac > .9) snowfac = 1.   

   vf = veg_fracarea * (1. - snowfac)
   vfc = 1. - vf

   fcpct = soil_water / slmsts(ntext_soil)

   if (leaf_class == 2) then

      alg = .50  ! Firn/glacier albedo

   else

      if (fcpct > .5) then
         alg = .14
      else
         alg = .31 - .34 * fcpct
      endif

   endif

   absg = (1. - alg) * rad
   algs = 1. - absg
   do k = nlev_sfcwater,1,-1
      algs = algs - fracabs(k)
      rshort_s(k) = rshort * vfc * fracabs(k)
   enddo

   albedo_beam = vf * alv + vfc * vfc * algs

   rshort_g = rshort * vfc * absg
   rshort_s(nlev_sfcwater+1:nzs) = 0.
   rshort_v = rshort * vf * (1. - alv + vfc * algs)
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

return
end subroutine sfcrad_land

