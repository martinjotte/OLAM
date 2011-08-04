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
Module decomposition_coms
use mem_ed, only: n_pft

real, parameter :: resp_opt_water = 0.94  !  Optimal soil porosity, as a fraction of total porosity, for heterotrophic respiration (dimensionless).
real, parameter :: resp_water_below_opt = 4.96  ! Determines rate at which heterotrophic respiration declines for relative porosities below resp_opt_water (dimensionless).
real, parameter :: resp_water_above_opt = 4.88  ! Determines rate at which heterotrophic respiration declines for relative porosities above resp_opt_water (dimensionless).
real, parameter :: resp_temperature_increase = 0.0741  ! Determines how rapidly heterotrophic respiration increases with increasing temperature (1/K).

real, parameter :: N_immobil_supply_scale = 40.0 / 365.25! Supply coefficient for nitrogen immobilization (1/day)

real, parameter :: cwd_frac = 0.2  ! Fraction of structural material that goes to coarse woody debris upon mortality.  Note that currently CWD decomposed at a rate identical to structural soil C.

real, parameter :: r_fsc=1.0   ! Fraction of structural pool decomposition going to heterotrophic respiration

real, parameter :: r_stsc=0.3  ! Fraction of structural pool decomposition going to heterotrophic respiration instead of the slow pool.

real, parameter :: r_ssc=1.0   ! Fraction of structural pool decomposition going to heterotrophic respiration

real, parameter :: K1=4.5 / 365.25 ! Intrinsic decay rate of fast pool soil carbon (1/days); this is modulated by Lc

real, parameter :: K2=11.0 / 365.25 ! Intrinsic decay rate of fast pool soil carbon (1/days)

real, parameter :: K3=100.2 / 365.25 ! Intrinsic decay rate of slow pool soil carbon (1/days).  This pool has already decayed from the structural pool.

real, dimension(n_pft), parameter :: f_labile =  (/  &
     1.0,  & ! C4 grass
     1.0,  & ! early successional broadleaf evergreen
     1.0,  & ! mid successional broadleaf evergreen
     1.0,  & ! late successional broadleaf evergreen
     1.0,  & ! C3 grass
     0.79, & ! northern pines
     0.79, & ! southern pines
     0.79, & ! late successional conifers
     0.79, & ! early successional broadleaf deciduous
     0.79, & ! early successional broadleaf deciduous
     0.79  /) ! early successional broadleaf deciduous

end Module decomposition_coms
