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
subroutine ed_vegetation_dynamics()

  use ed_structure_defs
  use mem_leaf, only: first_site
  use misc_coms, only: simtime, current_time, dtlm
  use leaf_coms, only: dt_leaf
  use ed_options, only: include_fire, frq_phenology
  use lphys_interface, only: apply_disturbances, reproduction

  implicit none

  type(site), pointer :: cs
  real :: tfact1
  real :: tfact2
  integer :: doy
  integer :: julday
  
  ! find the day of year
  doy = julday(current_time%month, current_time%date, current_time%year)
  
  ! Time factor for averaging FRQSTATE
  tfact1 = dt_leaf / frq_phenology

  ! Time factor for averaging dailies 
  tfact2 = frq_phenology / (365.25 * 86400.0)

  cs => first_site
  do while(associated(cs))

     call normalize_ed_daily_vars(cs)

     call phenology_driver(cs, doy, current_time%month, tfact1)

     call dbalive_dt(cs, tfact2)

     if(current_time%date == 1)then

        call structural_growth(cs, current_time%month)
        call reproduction(cs, current_time%month)
        cs%min_monthly_temp = 500.0

        if(include_fire == 1)call fire_frequency(current_time%month, cs)
        call site_disturbance_rates(current_time%month, current_time%year, cs)

        if(current_time%month == 1)then

           call fuse_patches(cs)
           call apply_disturbances(cs)

        endif

     endif

     call update_C_and_N_pools(cs)

     call reinitialize_ed_daily_vars(cs)

     cs => cs%next_site
  enddo

  return
end subroutine ed_vegetation_dynamics
