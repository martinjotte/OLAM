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
subroutine soil_respiration(ed_patch)

  use ed_structure_defs
  use leaf_coms, only: nzg, slmsts
  use pft_coms, only: root_respiration_factor, q, qsw

  implicit none

  type(patch) :: ed_patch
  
  real :: r_resp_temp_fac
  type(cohort), pointer :: cc
  real :: Lc
  real :: r_resp

  ! This is the temperature dependence of root respiration.  Same for all 
  ! cohorts.
  r_resp_temp_fac = 1.0 / (1.0 + exp(0.4   &
       * ( 278.15 - ed_patch%soil_tempk(nzg) ) ) ) &
       / (1.0 + exp(0.4 * ( ed_patch%soil_tempk(nzg) - 318.15 ) ) )  &
       * exp( 10.41 - 3000.0/ed_patch%soil_tempk(nzg) )
  
  cc => ed_patch%tallest
  do while(associated(cc))

     r_resp = root_respiration_factor(cc%pft) * r_resp_temp_fac *   &
          cc%balive * q(cc%pft) / (1.0 + q(cc%pft) + qsw(cc%pft) *   &
          cc%hite) * cc%nplant
     cc%omean_root_resp = cc%omean_root_resp + r_resp
     cc%dmean_root_resp = cc%dmean_root_resp + r_resp
          
     cc => cc%shorter
  enddo

  ! Compute soil/temperature modulation of heterotrophic respiration
  call resp_index(ed_patch%ntext_soil(nzg),ed_patch%soil_tempk(nzg),  &
       ed_patch%soil_water(nzg),slmsts(ed_patch%ntext_soil(nzg)),  &
       ed_patch%A_decomp)

  ! Compute nitrogen immobilization factor
  call resp_f_decomp(ed_patch,Lc)

  ! Compute heterotrophic respiration
  call resp_rh(ed_patch,Lc)

  ! Update averaged variables
  ed_patch%dmean_A_decomp = ed_patch%dmean_A_decomp + ed_patch%A_decomp
  ed_patch%dmean_Af_decomp = ed_patch%dmean_Af_decomp +   &
       ed_patch%A_decomp * ed_patch%f_decomp
  ed_patch%omean_rh = ed_patch%omean_rh + ed_patch%rh

  return
end subroutine soil_respiration

!=====================================================================

subroutine resp_index(nsoil,tempk,theta,slmsts,resp_weight)

  use decomposition_coms, only: resp_temperature_increase, resp_opt_water,  &
       resp_water_below_opt, resp_water_above_opt

  implicit none

  integer :: nsoil
  real :: tempk
  real :: theta
  real :: slmsts
  real :: resp_weight
  real :: temperature_limitation
  real :: water_limitation
  real :: Ws

  ! temperature dependence
  temperature_limitation = min(1.0,exp(resp_temperature_increase *   &
       (tempk-318.15)))

  ! moisture dependence
  Ws = theta/slmsts
  if(Ws.le.resp_opt_water)then
     water_limitation = exp( (Ws - resp_opt_water) * resp_water_below_opt)
  else
     water_limitation = exp((resp_opt_water-Ws) * resp_water_above_opt)
  endif
  
  ! compute the weight
  resp_weight = temperature_limitation * water_limitation
     
  return
end subroutine resp_index

!===================================================================
subroutine resp_f_decomp(cp,Lc)

  use ed_structure_defs
  use decomposition_coms, only: r_stsc, N_immobil_supply_scale, K1
  use pft_coms, only: c2n_structural, c2n_slow
  use ed_options, only: n_decomp_lim

  implicit none

  type(patch) :: cp

  real :: N_immobilization_demand
  real, intent(out) :: Lc

  if(cp%structural_soil_C > 0.0)then
     if(cp%structural_soil_L == cp%structural_soil_C)then
        Lc = 0.049787 ! = exp(-3.0)
     else
        Lc = exp(-3.0*cp%structural_soil_L/cp%structural_soil_C)
     endif
  else
     Lc=0.0
  endif
  
  if(n_decomp_lim == 1)then
     N_immobilization_demand = cp%A_decomp * Lc * K1   &
          * cp%structural_soil_C   &
          *((1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural) 
     
     cp%f_decomp = N_immobil_supply_scale * cp%mineralized_soil_N   &
          / (N_immobilization_demand + N_immobil_supply_scale *   &
          cp%mineralized_soil_N)
  else
     ! Option for no plant N limitation
     cp%f_decomp = 1.0
  endif

  return
end subroutine resp_f_decomp

!=============================================================

subroutine resp_rh(cp,Lc)

  use decomposition_coms, only: K1, K2, K3, r_fsc, r_ssc, r_stsc, cwd_frac
  use ed_structure_defs

  implicit none

  type(patch) :: cp
  real, intent(in) :: Lc
  real :: fast_C_loss
  real :: structural_C_loss
  real :: slow_C_loss

  ! These have units [kgC/m2/day]
  fast_C_loss = cp%A_decomp * K2 * cp%fast_soil_C
  structural_C_loss = cp%A_decomp * Lc * K1   &
       * cp%structural_soil_C * cp%f_decomp
  slow_C_loss = cp%A_decomp * K3 * cp%slow_soil_C

  ! Unit conversion is (kg C)/day  to (umol CO2)/s
  cp%rh = 964.5062 * (r_fsc*fast_C_loss + r_stsc*structural_C_loss    &
       + r_ssc*slow_C_loss)
  cp%cwd_rh = 964.5062 * (r_stsc*structural_C_loss + r_ssc*slow_C_loss) *   &
       cwd_frac

  return
end subroutine resp_rh

!==========================================================================

subroutine update_C_and_N_pools(cs)
  
  use ed_structure_defs
  use decomposition_coms, only: K1, K2, K3, r_stsc
  use pft_coms, only: c2n_slow, c2n_structural

  implicit none

  type(site) :: cs
  type(patch), pointer :: ed_patch
  real :: Lc
  real :: fast_C_loss
  real :: fast_N_loss
  real :: structural_C_loss
  real :: structural_L_loss
  real :: slow_C_input
  real :: slow_C_loss
  real :: mineralized_N_input
  real :: mineralized_N_loss

  ed_patch => cs%oldest_patch
  do while(associated(ed_patch))

     if(ed_patch%structural_soil_C > 0.0)then
        if(ed_patch%structural_soil_L == ed_patch%structural_soil_C)then
           Lc = 0.049787 ! = exp(-3.0)
        else
           Lc = exp(-3.0*ed_patch%structural_soil_L/ed_patch%structural_soil_C)
        endif
     else
        Lc=0.0
     endif
        
     ! fast pools
     fast_C_loss = ed_patch%dmean_A_decomp * K2 * ed_patch%fast_soil_C
     fast_N_loss = ed_patch%dmean_A_decomp * K2 * ed_patch%fast_soil_N
        
     ! structural pools
     structural_C_loss = ed_patch%dmean_Af_decomp * Lc * K1 *   &
          ed_patch%structural_soil_C
     structural_L_loss = ed_patch%dmean_Af_decomp * Lc * K1 *   &
          ed_patch%structural_soil_L

     ! slow pools
     slow_C_input = (1.0 - r_stsc) * structural_C_loss
     slow_C_loss = ed_patch%dmean_A_decomp * K3 * ed_patch%slow_soil_C

     ! mineralized pool
     mineralized_N_input = fast_N_loss + slow_C_loss / c2n_slow
     mineralized_N_loss = ed_patch%total_plant_nitrogen_uptake +   &
          ed_patch%dmean_Af_decomp * Lc * K1 * ed_patch%structural_soil_C   &
          * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
        
     ! all C fluxes have units kgC/m2/frq_phenology, and we are updating on 
     ! the frq_phenology time step.
     ed_patch%fast_soil_C = ed_patch%fast_soil_C + ed_patch%fsc_in -   &
          fast_C_loss
     ed_patch%structural_soil_C = ed_patch%structural_soil_C +   &
          ed_patch%ssc_in - structural_C_loss
     ed_patch%structural_soil_L = ed_patch%structural_soil_L +   &
          ed_patch%ssl_in - structural_L_loss
     ed_patch%slow_soil_C = ed_patch%slow_soil_C + slow_C_input - slow_C_loss

     ! all N fluxes have units kgN/m2/frq_phenology, and we are updating on
     ! the frq_phenology time step
     ed_patch%fast_soil_N = ed_patch%fast_soil_N + ed_patch%fsn_in -   &
          fast_N_loss
     ed_patch%mineralized_soil_N = ed_patch%mineralized_soil_N +   &
          mineralized_N_input - mineralized_N_loss
        
     ! reset average variables
     ed_patch%dmean_A_decomp = 0.0
     ed_patch%dmean_Af_decomp = 0.0

     ! require pools to be >= 0.0
     ed_patch%fast_soil_C = max(0.0,ed_patch%fast_soil_C)
     ed_patch%structural_soil_C = max(0.0,ed_patch%structural_soil_C)
     ed_patch%structural_soil_L = max(0.0,ed_patch%structural_soil_L)
     ed_patch%slow_soil_C = max(0.0,ed_patch%slow_soil_C)
     ed_patch%fast_soil_N = max(0.0,ed_patch%fast_soil_N)
     ed_patch%mineralized_soil_N = max(0.0,ed_patch%mineralized_soil_N)
     
     ed_patch => ed_patch%younger
  enddo
  
  return
end subroutine update_C_and_N_pools
