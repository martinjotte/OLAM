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
subroutine reproduction(cs, month)

  use ed_structure_defs
  use pft_coms, only: nonlocal_dispersal, seedling_mortality, phenology,  &
       sla, c2n_stem, l2n_stem, c2n_recruit, seed_rain, include_pft,  &
       include_pft_ag, qsw, q, hgt_min, plant_min_temp
  use decomposition_coms, only: f_labile
  use leaf_coms, only: dt_leaf
  use fusion_fission_coms, only: min_recruit_size
  use lphys_interface, only: split_cohorts, insert_cohort

  implicit none

  type(site), target :: cs
  integer, intent(in) :: month
  type(patch), pointer :: ptarget
  type(patch), pointer :: psource
  type(cohort), pointer :: cc
  integer :: count_cohorts
  type(cohort), pointer :: dc
  real :: nplant
  real :: balive
  real :: bleaf
  real :: dbh2bd
  real :: dbh2bl
  real :: bdead
  real :: dbh
  real :: h2dbh
  real :: hite
  integer :: pft

  !--------------------------
  ! First, sort the cohorts after the growth
  !----------------------------
  ptarget => cs%oldest_patch
  do while(associated(ptarget))
     call sort_cohorts(ptarget)
     ptarget => ptarget%younger
  enddo

  !----------------------------
  ! Next, do the dispersal
  !------------------------------
  ptarget => cs%oldest_patch
  do while(associated(ptarget))
     psource => cs%oldest_patch
     do while(associated(psource))

        ! Non-local, gridcell-wide dispersal
        cc => psource%tallest
        do while(associated(cc))

           if(phenology(cc%pft) /= 2   .or.  &  ! for NOT broad leaf deciduous
                (cs%lat >= 0.0 .and. month == 6) .or.   &  ! or Jun in north
                (cs%lat < 0.0 .and. month == 12) )then     ! or Dec in south
              ptarget%repro(cc%pft) = ptarget%repro(cc%pft) +   &
                   nonlocal_dispersal(cc%pft) * cc%nplant *   &
                   ( 1.0 - seedling_mortality(cc%pft) ) &
                   * cc%bseeds * psource%area
           endif
           cc => cc%shorter
        enddo

        ! Local dispersal (seeds stay in this patch)
        if(associated(psource,target=ptarget))then
           cc => psource%tallest
           do while(associated(cc))
              if(phenology(cc%pft) /= 2 .or.  &  ! for NOT broad leaf deciduous
                   (cs%lat >= 0.0 .and. month == 6) .or.   &  ! or Jun in north
                   (cs%lat < 0.0 .and. month == 12) )then     ! or Dec in south
                 ptarget%repro(cc%pft) = ptarget%repro(cc%pft) +   &
                      cc%nplant * (1.0 - nonlocal_dispersal(cc%pft))   &
                      * ( 1.0 - seedling_mortality(cc%pft) ) * cc%bseeds
              endif
              cc => cc%shorter
           enddo
        endif

        psource => psource%younger
     enddo

     ptarget => ptarget%younger
  enddo

  !------------------------------
  ! Finally, do the recruitment 
  !------------------------------

  ptarget => cs%oldest_patch
  do while(associated(ptarget))
     do pft = 1, n_pft
        if(include_pft(pft) == 1 .and.   &
             cs%min_monthly_temp >= plant_min_temp(pft) - 5.0)then
           if(ptarget%dist_type /= 1 .or. include_pft_ag(pft) == 1)then
              hite = hgt_min(pft)
              dbh = h2dbh(hite, pft)
              bdead = dbh2bd(dbh, hite, pft)
              bleaf = dbh2bl(dbh, pft)
              balive = bleaf * (1.0 + q(pft) + qsw(pft) * hite)
              nplant = ptarget%repro(pft) / (balive + bdead)
              if(include_pft(pft) == 1) nplant = nplant + seed_rain(pft)
              if( (nplant * (balive + bdead)) > min_recruit_size)then
                 nullify(dc)
                 allocate(dc)
                 dc%pft = pft
                 dc%nplant = nplant
                 dc%hite = hite
                 dc%dbh = dbh
                 dc%bdead = bdead
                 dc%bleaf = bleaf
                 dc%phenology_status = 0
                 dc%balive = balive
                 dc%lai = dc%bleaf * dc%nplant * sla(pft)
                 dc%bstorage = 0.0
                 dc%hcapveg =  3.e4 * max(1.,.025 * dt_leaf)
                 dc%veg_temp = ptarget%can_temp
                 dc%veg_water = 0.0
                 dc%siteptr => cs
                 dc%patchptr => ptarget
                 call init_ed_cohort_vars(dc, ptarget)
                 dc%new_recruit_flag = 1

                 call insert_cohort(dc, ptarget)
                 ptarget%repro(pft) = 0.0
              endif
           else

              ! seed litter from an agricultural patch.
              ptarget%fast_soil_N = ptarget%fast_soil_N   &
                   + ptarget%repro(pft) / c2n_recruit(pft)
              ptarget%fast_soil_C = ptarget%fast_soil_C + ptarget%repro(pft)
              ptarget%repro(pft) = 0.0

           endif
        endif
     enddo
     ptarget => ptarget%younger
  enddo

  !----------------------------------------------
  !  Now terminate, fuse, and split cohorts
  !----------------------------------------------
  ptarget => cs%oldest_patch
  do while(associated(ptarget))
     if(associated(ptarget%tallest))then
        call terminate_cohorts(ptarget)
        call fuse_cohorts(ptarget)
        call split_cohorts(ptarget)
     endif

     ptarget%cohort_count = count_cohorts(ptarget)

     call update_patch_derived_props(ptarget)

     ptarget => ptarget%younger
  enddo

  call update_site_derived_props(cs, 0)

  return
end subroutine reproduction
