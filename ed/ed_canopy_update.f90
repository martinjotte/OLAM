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
subroutine ed_canopy_update(ed_patch, vels, rhos, prss, pcpg, qpcpg,  &
     wshed_canopy, qwshed_canopy, canair, canhcap, dt_leaf, hxfergc,  &
     hxferca, wxfergc, hxfersc, wxfersc, sxfer_r, ed_transp,          &
     ntext_soil, soil_water, soil_fracliq, lsl)

  use ed_structure_defs
  use leaf_coms, only: nzg
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only: alvi

  implicit none

  type(patch)      :: ed_patch
  real, intent(in) :: vels
  real, intent(in) :: rhos
  real, intent(in) :: prss
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(inout) :: wshed_canopy
  real, intent(inout) :: qwshed_canopy
  real, intent(in) :: canair
  real, intent(in) :: canhcap
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: hxferca
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_r
  real, dimension(nzg) :: ed_transp
  type(cohort), pointer :: cc
  
  ! Local
  real :: sum_lai_rbi
  integer :: ndims
  integer, dimension(nzg) :: ed_ktrans
  integer, dimension(nzg), intent(in) :: ntext_soil
  real, dimension(nzg), intent(in) :: soil_water
  real, dimension(nzg), intent(in) :: soil_fracliq
  integer, intent(in) :: lsl

  call canopy_photosynthesis(ed_patch, vels, rhos, prss, sum_lai_rbi,  &
       ed_ktrans, ntext_soil, soil_water, soil_fracliq, lsl)

  call canopy_precip_interception(ed_patch, pcpg, qpcpg, wshed_canopy,  &
       qwshed_canopy)

  ndims = 2
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        ndims = ndims + 2
     endif
     cc => cc%shorter
  enddo

!  call canopy_implicit_driver(ed_patch, ndims, rhos, canhcap, canair,   &
!       sum_lai_rbi, dt_leaf, hxfergc, hxferca, wxfergc, hxfersc, wxfersc,  &
!       sxfer_r, ed_transp, ed_ktrans)
  call canopy_explicit_driver(ed_patch, ndims, rhos, canhcap, canair,   &
       sum_lai_rbi, dt_leaf, hxfergc, hxferca, wxfergc, hxfersc, wxfersc,  &
       sxfer_r, ed_transp, ed_ktrans)

  ed_patch%omean_wflux = ed_patch%omean_wflux + (sxfer_r) / dt_leaf
  ed_patch%omean_latflux = ed_patch%omean_latflux + (wxfergc + wxfersc) /   &
       dt_leaf * alvi
  ed_patch%omean_hflux = ed_patch%omean_hflux + hxferca / dt_leaf

  return
end subroutine ed_canopy_update

!==============================================================

subroutine canopy_photosynthesis(ed_patch, vels, rhos, prss, sum_lai_rbi,  &
     ed_ktrans, ntext_soil, soil_water, soil_fracliq, lsl)

  use ed_structure_defs
  use pft_coms, only: n_pft, leaf_width, water_conductance, q, qsw, include_pft
  use canopy_radiation_coms, only: lai_min
  use leaf_coms, only: nzg, slpots, slmsts, slbs, soilcp, dslz

  implicit none

  type(patch)      :: ed_patch
  real, intent(in) :: vels
  real, intent(in) :: rhos
  real, intent(in) :: prss
  real, intent(out) :: sum_lai_rbi
  integer, dimension(nzg), intent(out) :: ed_ktrans
  integer, dimension(nzg), intent(in) :: ntext_soil
  real, dimension(nzg), intent(in) :: soil_water
  real, dimension(nzg), intent(in) :: soil_fracliq

  integer, dimension(nzg) :: root_depth_indices
  type(cohort), pointer :: cc
  integer :: ipft
  real :: P_op
  real :: P_cl
  real :: leaf_resp
  real :: cumulative_lai
  real :: water_demand
  real :: water_supply
  real :: gpp
  integer :: k1
  integer :: k2
  real :: swp
  integer :: nts
  real :: slpotv
  real, dimension(nzg) :: available_liquid_water
  integer, intent(in) :: lsl

  ! calculate liquid water available for transpiration
  available_liquid_water(nzg) = max(0.0, 1.0e3 * dslz(nzg) *   &
       soil_fracliq(nzg) * (soil_water(nzg) - soilcp(ntext_soil(nzg))))
  do k1 = nzg-1, lsl, -1
     available_liquid_water(k1) = available_liquid_water(k1+1) +  &
          dslz(k1) * 1.0e3 * soil_fracliq(k1) * max(0.0, soil_water(k1) -  &
          soilcp(ntext_soil(k1)))
  enddo

  ! Initialize the array of maximum photosynthesis rates used in the 
  ! mortality function.
  ed_patch%A_o_max(1:n_pft) = 0.0
  ed_patch%A_c_max(1:n_pft) = 0.0

  ! Find the first cohort with leaves above snow
  cc => ed_patch%tallest
  cohort_with_leaves: do while(associated(cc))
     if(cc%lai > lai_min .and.   &
          cc%hite > ed_patch%total_snow_depth) exit cohort_with_leaves
     cc => cc%shorter
  enddo cohort_with_leaves

  if(associated(cc))then

     ! Compute maximum photosynthetic rates,ie, the rate the cohort would have 
     ! if it were at the top of the canopy (used for the mortality function)
     do ipft = 1, n_pft

        ! Compute aerodynamic resistance between leaf and canopy air space

        cc%rb = 1.0/(  &
             0.003*sqrt(vels/leaf_width(cc%pft))  &
             + 0.5 * (2.06e-5)   &
             * ( (1.6e8)*abs(cc%veg_temp-ed_patch%can_temp)  &
             *leaf_width(cc%pft)**3 )**(0.25) &
             /leaf_width(cc%pft))

        if(include_pft(ipft) == 1) call lphysiol_full(  &
             cc%veg_temp-273.15  &  ! Vegetation temperature (C)
             ,ed_patch%can_shv*29.0/18.0  &  ! canopy specific humidity (mol/mol)
             !,initp%canopy_co2*1.0e-6  &
             ,2.8e-4  &                      ! effective CO2 mixing ratio (mol/mol)
             !,3.7e-4  &
             ,cc%par_v/cc%lai  &             ! absorbed PAR (Ein/m2 leaf/s)
             ,cc%rb  &                       ! aerodynamic resistance (s/m)
             ,rhos  &                        ! air density (kg/m3)
             ,ed_patch%A_o_max(ipft)   &     ! maximum open photosynthetic rate (umol/m2 leaf/s)
             ,ed_patch%A_c_max(ipft)  &      ! maximum closed photosynthetic rate (umol/m2 leaf/s)
             ,P_op  & ! open stomata resistance for water [s/m]
             ,P_cl  & ! closed stomata resistance for water [s/m]
             ,ipft  & ! PFT
             ,prss  & ! pressure (kg/m/s2)
             ,leaf_resp  &  ! leaf respiration rate (umol/m2 leaf/s)
             ,ed_patch%siteptr%green_leaf_factor(ipft)  &  ! fraction of actual green leaves relative to on-allometry value.
             ,ed_patch%siteptr%leaf_aging_factor(ipft)  &
             ,ed_patch%old_stoma_data_max(ipft)) ! type containing the exact stomatal derivatives and met info
     enddo

  else

     ed_patch%A_o_max(1:n_pft) = 0.0
     ed_patch%A_c_max(1:n_pft) = 0.0

  endif

  ! cumulative LAI
  cumulative_lai = 0.0
  
  ! LAI/rb, summed over all cohorts. Used in the implicit scheme.
  sum_lai_rbi = 0.0
  
  ! Initialize variables for transpiration calculation
  root_depth_indices(:) = 0

  ! Loop over all cohorts
  
  cc => ed_patch%tallest
  do while(associated(cc))
     
     ! Only need to worry about photosyn if radiative transfer has been 
     ! done for this cohort
     
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        
        ! Aerodynamic resistance [s/m]
        cc%rb = 1.0/(  &
             0.003*sqrt(vels*exp(-0.5*cumulative_lai)  &
             /leaf_width(cc%pft))  &
             + 0.5 * (2.06e-5)   &
             * ( (1.6e8)*abs(cc%veg_temp-ed_patch%can_temp)  &
             *leaf_width(cc%pft)**3 )**(0.25) &
             /leaf_width(cc%pft))

        call lphysiol_full(cc%veg_temp-273.15  &  ! Vegetation temperature (C)
             ,ed_patch%can_shv*29.0/18.0  &  ! canopy specific humidity (mol/mol)
             !,initp%canopy_co2*1.0e-6  &
             ,2.8e-4  &                      ! effective CO2 mixing ratio (mol/mol)
             !,3.7e-4  &
             ,cc%par_v/cc%lai  &             ! absorbed PAR (Ein/m2 leaf/s)
             ,cc%rb  &                       ! aerodynamic resistance (s/m)
             ,rhos  &                        ! air density (kg/m3)
             ,cc%A_open   &     ! maximum open photosynthetic rate (umol/m2 leaf/s)
             ,cc%A_closed  &      ! maximum closed photosynthetic rate (umol/m2 leaf/s)
             ,cc%rsw_open  & ! open stomata resistance for water [s/m]
             ,cc%rsw_closed  & ! closed stomata resistance for water [s/m]
             ,cc%pft  & ! PFT
             ,prss  & ! pressure (kg/m/s2)
             ,leaf_resp  &  ! leaf respiration rate (umol/m2 leaf/s)
             ,ed_patch%siteptr%green_leaf_factor(cc%pft)  &  ! fraction of actual green leaves relative to on-allometry value.
             ,ed_patch%siteptr%leaf_aging_factor(cc%pft)  &
             ,cc%old_stoma_data) ! type containing the exact stomatal derivatives and met info
        
        ! Term for the implicit scheme
        sum_lai_rbi = sum_lai_rbi + cc%lai / cc%rb

        ! Leaf respiration
        cc%omean_leaf_resp = cc%omean_leaf_resp + leaf_resp * cc%lai
        cc%dmean_leaf_resp = cc%dmean_leaf_resp + leaf_resp * cc%lai

        ! Demand for water
        water_demand = cc%lai * cc%Psi_open ! kg/m2/s; Psi_open is from last time step    

        ! Supply of water
        water_supply = water_conductance(cc%pft) *   &
             available_liquid_water(cc%krdepth) * 1.0e-3 *   &
             q(cc%pft) * cc%balive / (1.0 + q(cc%pft) + cc%hite *   &
             qsw(cc%pft)) * cc%nplant

        root_depth_indices(cc%krdepth) = 1

        ! Weighting between open/closed stomata
        cc%fsw = water_supply / (water_supply + water_demand)

        ! Account for nitrogen limitation
        cc%fs_open = cc%fsw * cc%fsn

        ! Photorespiration can become important at high temperatures.  If so,
        ! close down the stomata.
        if(cc%A_open < cc%A_closed)cc%fs_open = 0.0

        ! Net stomatal resistance
        cc%stomatal_resistance = 1.0 / (cc%fs_open / cc%rsw_open +   &
             (1.0 - cc%fs_open) / cc%rsw_closed)

        ! GPP, averaged over frqstate
        gpp = cc%lai * (cc%fs_open * cc%A_open + (1.0 - cc%fs_open) *   &
             cc%A_closed + leaf_resp)
        cc%omean_gpp = cc%omean_gpp + gpp

        ! GPP, averaged over frq_phenology
        cc%dmean_gpp = cc%dmean_gpp + gpp

        ! Potential GPP if no N limitation
        cc%dmean_gpp_pot = cc%dmean_gpp_pot + cc%lai * (cc%fsw * cc%A_open +  &
             (1.0 - cc%fsw) * cc%A_closed + leaf_resp)

        ! Maximum GPP if at the top of the canopy
        cc%dmean_gpp_max = cc%dmean_gpp_max + cc%lai * (cc%fs_open *   &
             ed_patch%A_o_max(cc%pft) + (1.0 - cc%fs_open) *   &
             ed_patch%A_c_max(cc%pft) + leaf_resp)

        ! Update cumulative LAI
        cumulative_lai = cumulative_lai + cc%lai

     else

        ! This is the case where a cohort does not have leaves or is 
        ! snow-covered
        cc%A_open = 0.0
        cc%A_closed = 0.0
        cc%Psi_open = 0.0
        cc%Psi_closed = 0.0
        cc%rsw_open = 0.0
        cc%rsw_closed = 0.0
        cc%rb = 0.0
     endif

     cc => cc%shorter
  enddo

  ! For plants of a given rooting depth, determine soil level from which 
  ! transpired water is to be extracted.
  ed_ktrans(:) = 0
  do k1 = lsl, nzg
     swp = -1.e10
     if(root_depth_indices(k1) == 1)then
        do k2 = k1, nzg
           nts = ntext_soil(k2)
           slpotv = slpots(nts) * (slmsts(nts) / soil_water(k2)) ** slbs(nts)
           ! Multiply by liquid fraction (ice is unavailable for transpiration)
           slpotv = slpotv * soil_fracliq(k2)
           ! Find layer in root zone with highest slpotv AND soil_water 
           ! above minimum soilcp.  Set ktrans to this layer.
           if (slpotv > swp .and. soil_water(k2) > soilcp(nts)) then
              swp = slpotv
              ed_ktrans(k1) = k2
           endif
        enddo
     endif
  enddo

  return
end subroutine canopy_photosynthesis

!==============================================================

subroutine canopy_precip_interception(ed_patch, pcpg, qpcpg, wshed_canopy,  &
     qwshed_canopy)

  use ed_structure_defs
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only: cice, cliq, alli

  implicit none

  type(patch)      :: ed_patch
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(inout) :: wshed_canopy
  real, intent(inout) :: qwshed_canopy

  real :: laii
  type(cohort), pointer :: cc
  real :: tvegc
  real :: qwtot
  real :: lai_fraction
  real :: fracliqv
  real :: wshed_layer
  real :: hcapveg_factor

  laii = 0.0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%hite > ed_patch%total_snow_depth)then
        laii = laii + cc%lai
     endif
     cc => cc%shorter
  enddo

  if (laii > lai_min) then
     
     laii = 1.0 / laii
     
     cc => ed_patch%tallest
     do while(associated(cc))
        if(cc%hite > ed_patch%total_snow_depth .and. cc%lai > lai_min)then
           tvegc = cc%veg_temp - 273.15

           hcapveg_factor = cc%lai / ed_patch%lai

! If precipitation, add intercepted mass and energy to vegetation surface
!   (Ignore energy of water already on vegetation surface)

! Vegetation layers intercept precipitation in proportion to their LAI
           lai_fraction = cc%lai * laii
           if(tvegc > 0.0)then
              qwtot = cc%hcapveg * hcapveg_factor * tvegc + qpcpg * lai_fraction +   &
                   cc%veg_water * (cliq * tvegc + alli)
           else
              qwtot = cc%hcapveg * hcapveg_factor * tvegc + qpcpg * lai_fraction +   &
                   cc%veg_water * cice * tvegc
           endif
           cc%veg_water = cc%veg_water + pcpg * lai_fraction

! Compute equilbrium temperature of veg + precipitation

           call qwtk(qwtot, cc%veg_water, cc%hcapveg * hcapveg_factor, cc%veg_temp,  &
                fracliqv)
           tvegc = cc%veg_temp - 273.15
      
! Shed any excess intercepted precipitation and its energy

           if (cc%veg_water > .22 * cc%lai) then
              wshed_layer = cc%veg_water - .22 * cc%lai
              wshed_canopy = wshed_canopy + wshed_layer

              if (fracliqv <= .0001) then
                 qwshed_canopy = qwshed_canopy + cice * tvegc * wshed_layer
              else
                 qwshed_canopy = qwshed_canopy +   &
                      (cliq * tvegc + fracliqv * alli) * wshed_layer
              endif
              
              cc%veg_water = cc%veg_water - wshed_layer
           endif

        endif
        cc => cc%shorter
     enddo
  else
     wshed_canopy = pcpg
     qwshed_canopy = qpcpg
  endif

  return
end subroutine canopy_precip_interception

!==============================================================

subroutine canopy_implicit_driver(ed_patch, ndims, rhos, canhcap, canair,  &
     sum_lai_rbi, dt_leaf, hxfergc, hxferca, wxfergc, hxfersc, wxfersc,    &
     sxfer_r, ed_transp, ed_ktrans)

!-----------------------------------------------------------------------------
! Execute the implicit exchange of heat and moisture between vegetation and 
! canopy air
!-----------------------------------------------------------------------------

  use ed_structure_defs
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only: cp, alvl, alvi
  use leaf_coms, only: nzg

  implicit none

  type(patch)         :: ed_patch
  integer, intent(in) :: ndims
  real, intent(in) :: rhos
  real, intent(in) :: canhcap
  real, intent(in) :: canair
  real, intent(in) :: sum_lai_rbi
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: hxferca
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_r
  real, dimension(nzg), intent(out) :: ed_transp
  integer, dimension(nzg), intent(in) :: ed_ktrans

  ! Locals

  real, dimension(ndims, ndims) :: t_evolve_matrix
  integer :: ic
  type(cohort), pointer :: cc
  real :: tvegc
  real :: veg_rhovs
  real :: rhovsil
  real :: veg_rhovsp
  real :: vp_gradient
  real :: a1
  real :: sigmaw
  real :: dsigmaw_dW
  real :: a4
  real :: et_conductance
  real, dimension(ndims) :: explicit_deriv_portion
  real, dimension(ndims) :: original_state
  integer :: ind1
  integer :: ind2
  integer, dimension(ndims) :: indx
  real :: d
  real, dimension(ndims) :: implicit_new_state
  integer :: k
  real :: mult
  real :: a10
  integer :: idim

  ! Initialize time evolution matrix and the explicit contribution
  t_evolve_matrix(1:ndims, 1:ndims) = 0.0
  explicit_deriv_portion(1:ndims) = 0.0
  ed_transp(:) = 0.0

  ! explicitly integrated contribution to the canopy air temperature
  explicit_deriv_portion(1) = (hxfergc + hxfersc - hxferca) /   &
       (dt_leaf * canhcap)

  ! Canopy air temperature  (d can_temp, Tcan)
  t_evolve_matrix(1,1) = - 2.2 * cp * rhos / canhcap * sum_lai_rbi

  ! explicitly integrated contribution to the canopy specific humidity
  explicit_deriv_portion(2) = (wxfergc + wxfersc - sxfer_r) /   &
       (dt_leaf * canair)

  ic = 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then

        ! Set indices
        ic = ic + 1
        ind1 = 1 + 2 * ic
        ind2 = 2 + 2 * ic

        ! Compute heat transfer coefficient
        a10 = 2.2 * cp * rhos * cc%lai / cc%rb

        ! compute ET using variables at time n
        tvegc = cc%veg_temp - 273.15
        veg_rhovs= rhovsil(tvegc)
        veg_rhovsp = rhovsil(tvegc+1.0) - veg_rhovs
        vp_gradient = veg_rhovs - rhos * ed_patch%can_shv
        a1 = 2.2 * cc%lai / cc%rb
        sigmaw = min(1.,(cc%veg_water / (.22 * cc%lai))**.66667)
        if(cc%veg_water > 0.0)then
           dsigmaw_dW = 1.0 / (0.33 * cc%lai * sqrt(sigmaw))
        else
           dsigmaw_dW = 0.0
        endif
        a4 = cc%lai / (cc%rb + cc%stomatal_resistance)

        if(vp_gradient <= 0.0)then

           et_conductance = a1

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * a1

           ! d veg_water, can_shv
           mult = rhos * a1
           t_evolve_matrix(ind2,2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * ed_patch%can_shv

           ! d veg_water, veg_temp
           mult = - veg_rhovsp * a1
           t_evolve_matrix(ind2,ind1) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cc%veg_temp

        else

           et_conductance = a1 * sigmaw + a4 * (1.0 - sigmaw)

           ! d can_shv, veg_water
           mult = vp_gradient * (a1 - a4) * dsigmaw_dW / canair
           t_evolve_matrix(2,ind2) = t_evolve_matrix(2,ind2) + mult
           explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
                cc%veg_water

           ! d veg_temp, veg_water
           mult = -alvl / cc%hcapveg * vp_gradient * (a1 - a4) * dsigmaw_dW
           t_evolve_matrix(ind1,ind2) = mult
           explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
                cc%veg_water * mult

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * sigmaw * a1

           ! d veg_water, can_shv
           mult = rhos * sigmaw * a1
           t_evolve_matrix(ind2,2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * ed_patch%can_shv

           ! d veg_water, veg_temp
           mult = - veg_rhovsp * sigmaw * a1
           t_evolve_matrix(ind2,ind1) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cc%veg_temp

           ! d veg_water, veg_water
           mult = - vp_gradient * a1 * dsigmaw_dW
           t_evolve_matrix(ind2,ind2) = mult
           explicit_deriv_portion(ind2) = explicit_deriv_portion(ind2) -   &
                mult * cc%veg_water

        endif

        ! d can_temp, veg_temp
        t_evolve_matrix(1,ind1) = a10 / canhcap

        ! dcan_shv, explicit
        explicit_deriv_portion(2) = explicit_deriv_portion(2) +   &
             vp_gradient * et_conductance / canair

        ! dcan_shv, can_shv
        mult = - rhos / canair * et_conductance
        t_evolve_matrix(2,2) = t_evolve_matrix(2,2) + mult
        explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
             ed_patch%can_shv

        ! dcan_shv, veg_temp
        mult = veg_rhovsp * et_conductance / canair
        t_evolve_matrix(2,ind1) = t_evolve_matrix(2,ind1) + mult
        explicit_deriv_portion(2) = explicit_deriv_portion(2) - mult *   &
             cc%veg_temp

        ! dveg_temp, explicit
        explicit_deriv_portion(ind1) = (cc%rshort_v + cc%rlong_v) /   &
             cc%hcapveg - a10 / cc%hcapveg * (cc%veg_temp -   &
             ed_patch%can_temp) - alvl * et_conductance * vp_gradient /   &
             cc%hcapveg

        ! d veg_temp, can_temp
        mult = a10 / cc%hcapveg
        t_evolve_matrix(ind1,1) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             ed_patch%can_temp * mult

        ! d veg_temp, can_shv
        mult = rhos * alvl / cc%hcapveg * et_conductance
        t_evolve_matrix(ind1,2) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             ed_patch%can_shv * mult

        ! d veg_temp, veg_temp
        mult = - a10 / cc%hcapveg - alvl / cc%hcapveg * veg_rhovsp *   &
             et_conductance
        t_evolve_matrix(ind1,ind1) = mult
        explicit_deriv_portion(ind1) = explicit_deriv_portion(ind1) -   &
             cc%veg_temp * mult

        ! derivative matrix is done.  Now load initial state of the cohorts:

        original_state(ind1) = cc%veg_temp
        original_state(ind2) = cc%veg_water

     endif
     cc => cc%shorter
  enddo

  ! original state of canopy
  original_state(1) = ed_patch%can_temp
  original_state(2) = ed_patch%can_shv
  
  if(ic > 0)then
     
     ! Compute (Identity matrix) - dt_leaf * (Derivs matrix)
     do ic=1,ndims
        do k=1,ndims
           if(ic == k)then
              t_evolve_matrix(ic,k) = 1.0 - dt_leaf * t_evolve_matrix(ic,k)
           else
              t_evolve_matrix(ic,k) = - dt_leaf * t_evolve_matrix(ic,k)
           endif
        enddo
     enddo
     
     ! Do the matrix inversion
     call ludcmp_dble(t_evolve_matrix,ndims,ndims,indx,d)
     
  endif
     
  ! This is the vector that the inverse matrix needs to multiply
  implicit_new_state(1:ndims) = original_state(1:ndims) +   &
       explicit_deriv_portion(1:ndims) * dt_leaf

  ! Do the multiplication
  if(ic > 0)call lubksb_dble(t_evolve_matrix,ndims,ndims,indx,  &
       implicit_new_state)
     
  ! Load the new state into the patch and cohort structures
  ed_patch%can_temp = implicit_new_state(1)
  ed_patch%can_shv = implicit_new_state(2)
  idim = 1
  ic = 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        ic = ic + 1
        idim = idim + 2
        cc%veg_temp = implicit_new_state(idim)
        cc%veg_water = implicit_new_state(idim+1)
        
        veg_rhovs= rhovsil(original_state(idim)-273.15)
        veg_rhovsp = rhovsil(original_state(idim)-272.15) - veg_rhovs
        vp_gradient = veg_rhovs - rhos * original_state(2)
        if(vp_gradient > 0.0)then
          ! Calculate the resulting transpiration
           sigmaw = min(1.,(original_state(idim+1) / (.22 * cc%lai))**.66667)
           if(original_state(idim+1) > 0.0)then
              dsigmaw_dW = 1.0 / (0.33 * cc%lai * sqrt(sigmaw))
           else
              dsigmaw_dW = 0.0
           endif
           a4 = cc%lai / (cc%rb + cc%stomatal_resistance)
           ed_transp(ed_ktrans(cc%krdepth)) =   &
                ed_transp(ed_ktrans(cc%krdepth)) +  &
                vp_gradient * (1.0-sigmaw) * a4 - &
                rhos * (1.0 - sigmaw) * a4 * (ed_patch%can_shv -   &
                original_state(2)) + veg_rhovsp * (1.0 - sigmaw) * a4 *   &
                (cc%veg_temp - original_state(idim))  - vp_gradient *   &
                dsigmaw_dW * a4 * (cc%veg_water - original_state(idim+1))
        endif

     endif
     cc => cc%shorter
  enddo

  ed_transp(:) = ed_transp(:) * dt_leaf

  ed_patch%omean_latflux = ed_patch%omean_latflux +   &
       ((ed_patch%can_shv - original_state(2)) * canair -  & 
       (wxfergc + wxfersc - sxfer_r)) / dt_leaf * alvl

  ! Require veg_water >= 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        if(cc%veg_water < 0.0)then
           
           ! Take away from the canopy humidity...
           ed_patch%can_shv = ed_patch%can_shv + cc%veg_water / canair
           
           ! ... and set the veg_water to zero
           cc%veg_water = 0.0
           
        endif
     endif
     cc => cc%shorter
  enddo

  if(ed_patch%can_temp + 1.0 == 0.0 .or. ed_patch%can_temp > 400.0 .or.   &
       ed_patch%can_temp < 100.0 .or. ed_patch%can_shv < 0.0001)then
     print*,'can temp is bad'
     print*,ed_patch%can_temp,original_state(1),explicit_deriv_portion(1)
     print*,ed_patch%can_shv,original_state(2),explicit_deriv_portion(2)
     print*,original_state(1:ndims)
     print*,hxfergc/canhcap,hxfersc/canhcap,hxferca/canhcap,ed_patch%ustar
     stop
  endif

  ed_patch%avg_daily_temp = ed_patch%avg_daily_temp + ed_patch%can_temp

  return
end subroutine canopy_implicit_driver
!==============================================================

subroutine canopy_explicit_driver(ed_patch, ndims, rhos, canhcap, canair,  &
     sum_lai_rbi, dt_leaf, hxfergc, hxferca, wxfergc, hxfersc, wxfersc,    &
     sxfer_r, ed_transp, ed_ktrans)

!-----------------------------------------------------------------------------
! Execute the explicit exchange of heat and moisture between vegetation and 
! canopy air
!-----------------------------------------------------------------------------

  use ed_structure_defs
  use canopy_radiation_coms, only: lai_min
  use consts_coms, only: cp, alvl, alvi, cliq, cice, alli
  use leaf_coms, only: nzg

  implicit none

  type(patch)         :: ed_patch
  integer, intent(in) :: ndims
  real, intent(in) :: rhos
  real, intent(in) :: canhcap
  real, intent(in) :: canair
  real, intent(in) :: sum_lai_rbi
  real, intent(in) :: dt_leaf
  real, intent(in) :: hxfergc
  real, intent(in) :: hxferca
  real, intent(in) :: wxfergc
  real, intent(in) :: hxfersc
  real, intent(in) :: wxfersc
  real, intent(in) :: sxfer_r
  real, dimension(nzg), intent(out) :: ed_transp
  integer, dimension(nzg), intent(in) :: ed_ktrans

  ! Locals

  integer :: ic
  type(cohort), pointer :: cc
  real :: tvegc
  real :: veg_rhovs
  real :: rhovsil
  real :: veg_rhovsp
  real :: vp_gradient
  real :: a1
  real :: sigmaw
  real :: dsigmaw_dW
  real :: a4
  real :: et_conductance
  real, dimension(ndims) :: explicit_deriv_portion
  real, dimension(ndims) :: original_state
  integer :: ind1
  integer :: ind2
  real :: a10
  real, dimension(ndims) :: explicit_new_state
  integer :: idim
  real :: qwtot
  real :: fracliqv
  real :: dQdt

  ! Initialize time evolution matrix and the explicit contribution
  explicit_deriv_portion(1:ndims) = 0.0
  ed_transp(:) = 0.0

  ! explicitly integrated contribution to the canopy air temperature
  explicit_deriv_portion(1) = (hxfergc + hxfersc - hxferca) /   &
       (dt_leaf * canhcap)

  ! explicitly integrated contribution to the canopy specific humidity
  explicit_deriv_portion(2) = (wxfergc + wxfersc - sxfer_r) /   &
       (dt_leaf * canair)

  ic = 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then

        ! Set indices
        ic = ic + 1
        ind1 = 1 + 2 * ic
        ind2 = 2 + 2 * ic

        ! Compute heat transfer coefficient
        a10 = 2.2 * cp * rhos * cc%lai / cc%rb

        ! compute ET using variables at time n
        tvegc = cc%veg_temp - 273.15
        veg_rhovs= rhovsil(tvegc)
        veg_rhovsp = rhovsil(tvegc+1.0) - veg_rhovs
        vp_gradient = veg_rhovs - rhos * ed_patch%can_shv
        a1 = 2.2 * cc%lai / cc%rb
        sigmaw = min(1.,(cc%veg_water / (.22 * cc%lai))**.66667)
        if(cc%veg_water > 0.0)then
           dsigmaw_dW = 1.0 / (0.33 * cc%lai * sqrt(sigmaw))
        else
           dsigmaw_dW = 0.0
        endif
        a4 = cc%lai / (cc%rb + cc%stomatal_resistance)

        if(vp_gradient <= 0.0)then

           et_conductance = a1

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * a1

        else

           et_conductance = a1 * sigmaw + a4 * (1.0 - sigmaw)

           ! d veg_water, explicit
           explicit_deriv_portion(ind2) = - vp_gradient * sigmaw * a1
           ed_transp(ed_ktrans(cc%krdepth)) =   &
                ed_transp(ed_ktrans(cc%krdepth)) +  &
                vp_gradient * (1.0 - sigmaw) * a4

        endif

        ! contribution to dcan_temp/dt, from this vegetation layer's veg_temp
        explicit_deriv_portion(1) = explicit_deriv_portion(1) + a10 /  &
             canhcap * (cc%veg_temp - ed_patch%can_temp)

        ! contribution to dcan_shv/dt, from this vegetation layer's veg_water
        explicit_deriv_portion(2) = explicit_deriv_portion(2) +   &
             vp_gradient * et_conductance / canair

        ! dveg_temp, explicit
        dQdt = cc%rshort_v + cc%rlong_v -   &
             a10 * (cc%veg_temp - ed_patch%can_temp) -   &
             alvl * et_conductance * vp_gradient
        tvegc = cc%veg_temp - 273.15
        if(tvegc > 0.0)then
           ed_patch%omean_latflux = ed_patch%omean_latflux +   &
                vp_gradient * et_conductance * alvl -   &
                (cliq * tvegc + alli) * explicit_deriv_portion(ind2)
           explicit_deriv_portion(ind1) = dQdt / (cc%hcapveg * cc%lai /   &
                ed_patch%lai + cliq *   &
                cc%veg_water)
        else
           ed_patch%omean_latflux = ed_patch%omean_latflux +   &
                vp_gradient * et_conductance * alvl -   &
                cice * tvegc * explicit_deriv_portion(ind2)
           explicit_deriv_portion(ind1) = dQdt / (cc%hcapveg * cc%lai /   &
                ed_patch%lai + cice * cc%veg_water)
        endif

        ! Load initial state of the cohorts:
        original_state(ind1) = cc%veg_temp
        original_state(ind2) = cc%veg_water

     endif
     cc => cc%shorter
  enddo

  ! original state of canopy
  original_state(1) = ed_patch%can_temp
  original_state(2) = ed_patch%can_shv
  
  ! new state of canopy
  do ic=1,ndims
     explicit_new_state(ic) = original_state(ic) +   &
          explicit_deriv_portion(ic) * dt_leaf
  enddo

  ! Load the new state into the patch and cohort structures
  ed_patch%can_temp = explicit_new_state(1)
  ed_patch%can_shv = explicit_new_state(2)
  idim = 1
  ic = 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        ic = ic + 1
        idim = idim + 2
        cc%veg_temp = explicit_new_state(idim)
        cc%veg_water = explicit_new_state(idim+1)
     endif
     cc => cc%shorter
  enddo

  ed_transp(:) = ed_transp(:) * dt_leaf

  ! Require 0.22 LAI >= veg_water >= 0
  cc => ed_patch%tallest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > ed_patch%total_snow_depth)then
        if(cc%veg_water < 0.0)then
           
           ! Take away from the canopy humidity...
           ed_patch%can_shv = ed_patch%can_shv + cc%veg_water / canair
           
           ! ... and set the veg_water to zero
           cc%veg_water = 0.0
           
        endif
     endif
     cc => cc%shorter
  enddo

  ed_patch%avg_daily_temp = ed_patch%avg_daily_temp + ed_patch%can_temp

  return
end subroutine canopy_explicit_driver

!======================================================================

subroutine ludcmp(a,n,np,indx,d)

  implicit none
  
  real, parameter :: tiny=1.0e-20
  integer :: n
  integer :: np
  real :: d
  real, dimension(np,np) :: a
  integer, dimension(n) :: indx
  real, dimension(np) :: vv
  integer :: i
  integer :: j
  integer :: k
  integer :: imax
  real :: sum
  real :: aamax
  real :: dum

  d = 1.0
  do i=1,n
     aamax = 0.0
     do j=1,n
        if(abs(a(i,j)) > aamax)aamax = abs(a(i,j))
     enddo
     if(aamax == 0.0)then
        print*,'singular matrix in ludcmp'
        do j=1,n
           print*,i,j,abs(a(i,j)),aamax
        enddo
        stop
     endif
     vv(i) = 1.0 / aamax
  enddo

  do j=1,n
     if(j.gt.1)then
        do i=1,j-1
           sum = a(i,j)
           if(i > 1)then
              do k=1,i-1
                 sum = sum - a(i,k) * a(k,j)
              enddo
              a(i,j) = sum
           endif
        enddo
     endif
     aamax = 0.0
     do i=j,n
        sum = a(i,j)
        if (j > 1)then
           do k=1,j-1
              sum = sum - a(i,k) * a(k,j)
           enddo
           a(i,j) = sum
        endif
        dum = vv(i) * abs(sum)
        if(dum >= aamax)then
           imax = i
           aamax = dum
        endif
     enddo
     if(j /= imax)then
        do k=1,n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        enddo
        d = -d
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     if(j /= n)then
        if(a(j,j) == 0.0) a(j,j) = tiny
        dum = 1.0 / a(j,j)
        do i=j+1,n
           a(i,j) = a(i,j)*dum
        enddo
     endif
  enddo
  if(a(n,n) == 0.0)a(n,n) = tiny

  return
end subroutine ludcmp

!======================================================================

subroutine ludcmp_dble(a,n,np,indx,d)

  implicit none
  
  real, parameter :: tiny = 1.0d-20
  integer, intent(in) :: n
  integer, intent(in) :: np
  real, intent(out) :: d
  real, dimension(np,np), intent(inout) :: a
  real(kind=8), dimension(np,np) :: ad
  integer, dimension(n), intent(out) :: indx
  real(kind=8), dimension(np) :: vv
  integer :: i
  integer :: j
  integer :: k
  integer :: imax
  real(kind=8) :: sum
  real(kind=8) :: aamax
  real(kind=8) :: dum

  ad = dble(a)

  d = 1.0

  do i = 1, n
     aamax = 0.0d0
     do j = 1, n
        if(abs(ad(i,j)) > aamax)aamax = abs(ad(i,j))
     enddo
     if(aamax == 0.0d0)then
        print*,'singular matrix in ludcmp'
        do j=1,n
           print*,i,j,abs(a(i,j)),aamax
        enddo
        stop
     endif
     vv(i) = 1.0d0 / aamax
  enddo

  do j = 1, n
     if(j.gt.1)then
        do i = 1, j - 1
           sum = ad(i,j)
           if(i > 1)then
              do k = 1, i - 1
                 sum = sum - ad(i,k) * ad(k,j)
              enddo
              ad(i,j) = sum
           endif
        enddo
     endif
     aamax = 0.0d0
     do i=j,n
        sum = ad(i,j)
        if (j > 1)then
           do k = 1, j - 1
              sum = sum - ad(i,k) * ad(k,j)
           enddo
           ad(i,j) = sum
        endif
        dum = vv(i) * abs(sum)
        if(dum >= aamax)then
           imax = i
           aamax = dum
        endif
     enddo
     if(j /= imax)then
        do k = 1, n
           dum = ad(imax,k)
           ad(imax,k) = ad(j,k)
           ad(j,k) = dum
        enddo
        d = -d
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     if(j /= n)then
        if(ad(j,j) == 0.0d0) ad(j,j) = tiny
        dum = 1.0d0 / ad(j,j)
        do i = j + 1, n
           ad(i,j) = ad(i,j) * dum
        enddo
     endif
  enddo
  if(ad(n,n) == 0.0d0)ad(n,n) = tiny

  a = real(ad)

  return
end subroutine ludcmp_dble

!==================================================================

subroutine lubksb(a,n,np,indx,b)
  implicit none
  integer :: n
  integer :: np
  integer, dimension(n) :: indx
  real, dimension(n) :: b
  real, dimension(np,np) :: a
  integer :: ii
  integer :: i
  integer :: ll
  real :: sum
  integer :: j

  ii = 0
  do i=1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if(ii /= 0)then
        do j=ii,i-1
           sum = sum - a(i,j) * b(j)
        enddo
     elseif(sum /= 0.0)then
        ii = i
     endif
     b(i) = sum
  enddo
  do i=n,1,-1
     sum = b(i)
     if ( i < n )then
        do j=i+1,n
           sum = sum - a(i,j) * b(j)
        enddo
     endif
     b(i) = sum / a(i,i)
  enddo

  return
end subroutine lubksb

!==================================================================

subroutine lubksb_dble(a,n,np,indx,b)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: np
  integer, dimension(n), intent(in) :: indx
  real, dimension(n), intent(inout) :: b
  real(kind=8), dimension(n) :: bd
  real, dimension(np,np), intent(in) :: a
  real(kind=8), dimension(np,np) :: ad
  integer :: ii
  integer :: i
  integer :: ll
  real(kind=8) :: sum
  integer :: j

  ad = dble(a)
  bd = dble(b)

  ii = 0
  do i=1,n
     ll = indx(i)
     sum = bd(ll)
     bd(ll) = bd(i)
     if(ii /= 0)then
        do j=ii,i-1
           sum = sum - ad(i,j) * bd(j)
        enddo
     elseif(sum /= 0.0d0)then
        ii = i
     endif
     bd(i) = sum
  enddo
  do i=n,1,-1
     sum = bd(i)
     if ( i < n )then
        do j=i+1,n
           sum = sum - ad(i,j) * bd(j)
        enddo
     endif
     bd(i) = sum / ad(i,i)
  enddo

  b = real(bd)

  return
end subroutine lubksb_dble

