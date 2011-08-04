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
subroutine phenology_driver(cs, doy, month, tfact)

  use ed_structure_defs
  use ed_options, only: iphen_scheme

  implicit none

  type(site) :: cs
  integer :: doy
  integer :: month
  real :: tfact
  type(patch), pointer :: ed_patch

  ! Get the patch-level average daily temperature, which is needed for 
  ! mortality, recruitment and some phenology schemes

  ed_patch => cs%oldest_patch
  do while(associated(ed_patch))
     ed_patch%avg_daily_temp = ed_patch%avg_daily_temp * tfact
     ed_patch => ed_patch%younger
  enddo
     
  if(iphen_scheme == 0)then
     
     !  Default predictive scheme

     call update_thermal_sums(month, cs)
     call update_phenology(doy, cs)

  elseif(iphen_scheme == 1)then

     !  Scheme driven by remote sensing

  elseif(iphen_scheme == 2)then

     !  My new predictive scheme

  endif

  return
end subroutine phenology_driver

!======================================================================
subroutine update_thermal_sums(month, cs)
  
  use ed_structure_defs

  implicit none

  type(site) :: cs
  integer, intent(in) :: month

  type(patch), pointer :: ed_patch

  type(cohort), pointer :: cc
  integer :: i
  integer :: j
  integer :: ip

  ! Chill days --- number of days with average temperatures below 278.15 K
  ! Degree days --- sum of daily average temperatures above 278.15 K 

  ! loop over patches

  ed_patch => cs%oldest_patch
  do while(associated(ed_patch))
     
     ! Minimum monthly temperature of the site
     cs%min_monthly_temp = min(cs%min_monthly_temp, ed_patch%avg_daily_temp)
     
     if(ed_patch%avg_daily_temp > 278.15)then  
        
        ! update dgd
        
        if(cs%lat >= 0.0)then  
           
           ! northern hemisphere
           
           if(month <= 8)then 
              
              ! update only for Jan-Aug.
              
              ed_patch%sum_dgd = ed_patch%sum_dgd   &
                   + (ed_patch%avg_daily_temp-278.15)
              
           else 
              
              ! not accumulating degree days in Sep-Dec
              
              ed_patch%sum_dgd = 0.0
              
           endif
           
        else 
           
           ! in southern hemisphere
           
           if(month <= 2 .or. month >= 7)then 
              
              ! Accumulating only Jul-Feb
              
              ed_patch%sum_dgd = ed_patch%sum_dgd +   &
                   ed_patch%avg_daily_temp - 278.15
              
           else 
              
              ! not accumulating degree days Mar-Jun
              
              ed_patch%sum_dgd = 0.0
              
           endif
        endif
     else 
        
        ! update chilling days
        
        if(cs%lat >= 0.0)then  
           
           ! northern hemisphere
           
           if(month >= 11 .or. month <= 6)then 
              
              ! accumulate only Nov-Jun
              
              ed_patch%sum_chd = ed_patch%sum_chd + 1.0
              
           else 
              
              ! not accumulating chilling days, Jul-Oct
              
              ed_patch%sum_chd = 0.0
              
           endif
        else 
           
           ! in southern hemisphere
           
           if(month >= 5)then 
              
              ! accumulate only in May-Dec
              
              ed_patch%sum_chd = ed_patch%sum_chd + 1.0
              
           else 
              
              ! not accumulating chilling days Jan-Apr
              
              ed_patch%sum_chd = 0.0
              
           endif
        endif
     endif
     
     ed_patch => ed_patch%younger
  enddo

  return
end subroutine update_thermal_sums

!==================================================================

subroutine update_phenology(day, cs)

  use ed_structure_defs
  use leaf_coms, only: nzg
  use pft_coms, only: phenology, sla, c2n_leaf, q, qsw, l2n_stem, c2n_stem
  use decomposition_coms, only: f_labile
  use phenology_coms, only: retained_carbon_fraction, theta_crit
  use mem_leaf, only: land

  implicit none
  
  integer :: day
  type(site) :: cs
  integer :: isoil_lev
  real :: daylight
  real :: daylength
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  integer :: drop_cold
  integer :: leaf_out_cold
  real, dimension(nzg) :: theta
  real :: leaf_litter
  real :: bl

  ! Level to evaluate the soil temperature

  isoil_lev = nzg 

  ! Calculate daylength for this gridcell

  daylight = daylength(cs%lat,day)  ! day (= the Julian day) is input

  ! Loop over patches

  cp => cs%oldest_patch
  do while(associated(cp))

     ! Re-initialize litter inputs
     cp%fsc_in = 0.0
     cp%fsn_in = 0.0
     cp%ssc_in = 0.0
     cp%ssl_in = 0.0

     ! Determine what phenology thresholds have been crossed
     call phenology_thresholds(daylight, cp%soil_tempk(isoil_lev),   &
          cp%soil_water, cp%ntext_soil, cp%sum_chd, cp%sum_dgd, drop_cold,   &
          leaf_out_cold, theta, land%lsl(cs%iland))

     ! Loop over cohorts
     cc => cp%tallest
     do while(associated(cc))

        if(cc%phenology_status < 2 .and. phenology(cc%pft) == 2)then

           ! Is this a cold deciduous with leaves?

           if(drop_cold == 1)then

              ! drop leaves

              cc%phenology_status = 2 

              ! compute litter inputs

              leaf_litter = (1.0 - retained_carbon_fraction)  &
                   * cc%lai / sla(cc%pft)
              cp%fsc_in = cp%fsc_in + leaf_litter * f_labile(cc%pft)
              cp%fsn_in = cp%fsn_in + leaf_litter * f_labile(cc%pft) /   &
                   c2n_leaf(cc%pft)
              cp%ssc_in = cp%ssc_in + leaf_litter * (1.0 - f_labile(cc%pft))
              cp%ssl_in = cp%ssl_in + leaf_litter *   &
                   (1.0 - f_labile(cc%pft)) * l2n_stem / c2n_stem

              ! adjust plant carbon pools

              cc%lai = 0.0
              cc%bleaf = 0.0
              cc%balive = cc%balive - leaf_litter / cc%nplant
              cc%cb(13) = cc%cb(13) - leaf_litter / cc%nplant
              cc%cb_max(13) = cc%cb_max(13) - leaf_litter / cc%nplant

           endif

        elseif(cc%phenology_status == 2 .and. phenology(cc%pft) == 2)then 
           
           ! Cold deciduous? 

           if(leaf_out_cold == 1)then

              ! Update plant carbon pools
              cc%phenology_status = 1 ! 1 indicates leaves are growing
              bl = cc%balive / (1.0 + qsw(cc%pft) * cc%hite + q(cc%pft))
              cc%lai = cc%nplant * bl * sla(cc%pft)
              cc%veg_temp = cp%can_temp
              cc%veg_water = 0.0
              cc%bleaf = bl

           endif

        elseif(phenology(cc%pft) == 1)then 

           ! Drought deciduous?

           if(theta(cc%krdepth) < theta_crit)then

              !  it is time to drop leaves

              if(cc%phenology_status < 2)then

                 ! update litter pools
                 leaf_litter = (1.0 - retained_carbon_fraction) * cc%lai /   &
                      sla(cc%pft)
                 cp%fsc_in = cp%fsc_in + leaf_litter * f_labile(cc%pft)
                 cp%fsn_in = cp%fsn_in + leaf_litter * f_labile(cc%pft) /   &
                      c2n_leaf(cc%pft)
                 cp%ssc_in = cp%ssc_in + leaf_litter * (1.0 - f_labile(cc%pft))
                 cp%ssl_in = cp%ssl_in + leaf_litter *   &
                      (1.0 - f_labile(cc%pft)) * l2n_stem / c2n_stem

                 ! update plant carbon pools
                 cc%phenology_status = 2
                 cc%lai = 0.0
                 cc%bleaf = 0.0
                 cc%balive = cc%balive - leaf_litter / cc%nplant
                 cc%cb(13) = cc%cb(13) - leaf_litter / cc%nplant
                 cc%cb_max(13) = cc%cb_max(13) - leaf_litter/cc%nplant

              endif

           elseif(theta(cc%krdepth) > theta_crit)then

              ! it is time to flush

              if(cc%phenology_status == 2)then

                 ! update carbon pools
                 cc%phenology_status = 1
                 bl = cc%balive / (1.0 + qsw(cc%pft) * cc%hite + q(cc%pft))
                 cc%lai = cc%nplant * bl * sla(cc%pft)
                 cc%veg_temp = cp%can_temp
                 cc%veg_water = 0.0
                 cc%bleaf = bl

              endif
           endif  ! critical moisture

        endif  ! phenology type

        cc => cc%shorter
     enddo  ! cohorts

     cp => cp%younger
  enddo  ! patches

  return
end subroutine update_phenology

!=====================================================================

real function daylength(lat,day)

  use consts_coms, only: pio180

  implicit none
  
  real :: lat
  real :: arg
  integer :: day
  
  arg = -tan(lat*pio180)*tan(-23.5*pio180*cos(6.283/365.0*(float(day)+9.0)))
  if( abs(arg) < 1.0 )then
     daylength = 120.0 * acos(arg)/(15.0*pio180)
  elseif(arg >= 1.0)then
     daylength = 0.0
  elseif(arg <= -1.0)then
     daylength = 1440.0
  endif

  return
end function daylength

!=======================================================================

subroutine phenology_thresholds(daylight, soil_temp, soil_water, soil_class,  &
     sum_chd, sum_dgd, drop_cold, leaf_out_cold, theta, lsl)

  use leaf_coms, only: nzg, slmsts, slz
  use phenology_coms, only: dl_tr, st_tr1, st_tr2, phen_a, phen_b, phen_c

  implicit none

  real, intent(in) :: daylight
  real, intent(in) :: soil_temp
  real, dimension(nzg), intent(in) :: soil_water
  integer, dimension(nzg), intent(in) :: soil_class
  real, intent(inout) :: sum_dgd
  real, intent(inout) :: sum_chd
  integer, intent(out) :: drop_cold
  integer, intent(out) :: leaf_out_cold
  real, dimension(nzg), intent(out) :: theta
  real :: gdd_threshold
  integer :: k1
  integer :: k2
  integer, intent(in) :: lsl

  ! initialize
  drop_cold = 0
  leaf_out_cold = 0
  theta(1:nzg) = 0.0

  !  determine whether or not we have cold deciduous leaf drop

  if( (daylight <= dl_tr .and. soil_temp < st_tr1) .or.   &
       soil_temp < st_tr2 )then 
     
     ! there is leaf drop
     
     drop_cold = 1
     
     ! reset degree day and chill day counters
     
     sum_dgd = 0.0
     sum_chd = 0.0
     
  endif

  ! do we have cold deciduous leaf-out?

  gdd_threshold = phen_a + phen_b * exp(phen_c * sum_chd)
  if(sum_dgd >= gdd_threshold) leaf_out_cold = 1

  ! calculate average theta for drought deciduous PFTs.  The different
  ! k1's represent different rooting depths.
  
  theta(1:nzg) = 0.0
  do k1 = lsl, nzg
     do k2 = k1,nzg-1
        theta(k1) = theta(k1) + soil_water(k2) *   &
             (slz(k2+1)-slz(k2)) / slmsts(soil_class(k2))
     enddo
     theta(k1) = theta(k1) - soil_water(nzg) * slz(nzg) /  &
          slmsts(soil_class(nzg))
     theta(k1) = - theta(k1) / slz(k1)
  enddo

  return
end subroutine phenology_thresholds

