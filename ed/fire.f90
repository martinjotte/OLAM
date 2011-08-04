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
subroutine fire_frequency(month, cs)

  use ed_structure_defs
  use pft_coms, only: agf_bs, qsw, q
  use leaf_coms, only: nzg, slz
  use disturbance_coms, only: fire_dryness_threshold, fire_parameter
  use mem_leaf, only: land

  implicit none

  integer, intent(in) :: month
  type(site)          :: cs

  real :: ignition_rate
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  real :: babove
  real :: patch_water_depth
  integer :: k
  real, external :: ed_agb
  real :: fuel

  ! Initialize
  ignition_rate = 0.0

  ! loop over patches
  cp => cs%youngest_patch
  do while(associated(cp))  
     
     ! Initialize patches
     fuel = 0.0
     
     ! add up fuel from all the cohorts
     cc => cp%shortest
     do while(associated(cc))
        babove = ed_agb(cc%bdead, cc%balive, cc%bleaf, cc%pft, cc%hite,   &
             cc%bstorage) * cc%nplant
        fuel = fuel + babove
        cc => cc%taller
     enddo
     
     ! calculate patch water in meters
     patch_water_depth = sum(cp%sfcwater_depth(1:cp%nlev_sfcwater)) *  &
          0.001 - cp%soil_water(nzg) * slz(nzg)
     do k = land%lsl(cs%iland), nzg-1
        patch_water_depth = patch_water_depth + cp%soil_water(k)  &
             *(slz(k+1) - slz(k))
     enddo

     ! calculate patch contribution to the ignition rate
     if(patch_water_depth < fire_dryness_threshold)then
        ignition_rate = ignition_rate + fuel * cp%area
     endif

     ! end loop over patches
     cp => cp%older
  enddo
  
  ! calculate fire dist rate [1/month]
  cs%lambda_fire(month) = fire_parameter * ignition_rate

  return
end subroutine fire_frequency
