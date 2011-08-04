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
subroutine sort_cohorts(cp)

  use ed_structure_defs

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: tallestc
  type(cohort), pointer :: shortestc
  type(cohort), pointer :: currentc
  type(cohort), pointer :: next_c
  type(cohort), pointer :: shortptr
  type(cohort), pointer :: tallptr
  type(cohort), pointer :: shortestc_save
  real :: tsp

  nullify(tallestc)
  nullify(shortestc)
  currentc => cp%tallest
  do while(associated(currentc))
     next_c => currentc%shorter

     tsp = currentc%hite

     tallptr => shortestc
     shortestc_save => shortestc
     next_taller: do while(associated(shortestc))
        if(shortestc%hite.ge.tsp)exit next_taller
        tallptr => tallptr%taller
        shortestc => shortestc%taller
     enddo next_taller
     shortestc => shortestc_save

     if(.not.associated(tallptr))then
        shortptr => tallestc
        tallestc => currentc
     else
        shortptr => tallptr%shorter
        tallptr%shorter => currentc
     endif

     if(.not.associated(shortptr))then
        shortestc => currentc
     else
        shortptr%taller => currentc
     endif

     currentc%taller => tallptr
     currentc%shorter => shortptr

     currentc => next_c
  enddo

  cp%shortest => shortestc
  cp%tallest => tallestc

  return

end subroutine sort_cohorts

!---------------------------------------------------------

subroutine terminate_cohorts(cp)

  use ed_structure_defs
  use fusion_fission_coms, only: min_recruit_size
  use pft_coms, only: l2n_stem, c2n_stem, c2n_storage, c2n_leaf
  use decomposition_coms, only: f_labile

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: cc
  type(cohort), pointer :: nextc
  integer :: count_cohorts
  real :: csize

  cc => cp%tallest
  do while(associated(cc))
     nextc => cc%shorter
     csize = cc%nplant * (cc%balive + cc%bdead + cc%bstorage)
     if(csize < (0.1 * min_recruit_size) )then

        if(.not.associated(cc%taller))then
           cp%tallest => cc%shorter
        else
           cc%taller%shorter => cc%shorter
        endif

        if(.not.associated(cc%shorter))then
           cp%shortest => cc%taller
        else
           cc%shorter%taller => cc%taller
        endif

        ! Update litter pools

        cp%fsc_in = cp%fsc_in + cc%nplant * (f_labile(cc%pft) *   &
             cc%balive + cc%bstorage)
        cp%fsn_in = cp%fsn_in + cc%nplant * (f_labile(cc%pft) * cc%balive /   &
             c2n_leaf(cc%pft) + cc%bstorage / c2n_storage)
        cp%ssc_in = cp%ssc_in + cc%nplant * ((1.0 - f_labile(cc%pft)) *  &
             cc%balive + cc%bdead)
        cp%ssl_in = cp%ssl_in + cc%nplant * ( (1.0 - f_labile(cc%pft)) *  &
             cc%balive + cc%bdead ) * l2n_stem / c2n_stem

        deallocate(cc)

     endif
     cc => nextc
  enddo

  cp%cohort_count = count_cohorts(cp)

  return
end subroutine terminate_cohorts

!===================================================================

integer function count_cohorts(cp)

  use ed_structure_defs

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: cc
  
  count_cohorts = 0
  cc => cp%tallest
  do while(associated(cc))
     count_cohorts = count_cohorts + 1
     cc => cc%shorter
  enddo

  return
end function count_cohorts

!----------------------------------------------------------

subroutine fuse_cohorts(cp)

  use ed_structure_defs
  use pft_coms, only: rho, b1Ht, max_dbh, sla
  use fusion_fission_coms, only: fusetol_h, fusetol, lai_fuse_tol

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: cc
  integer :: fusion_took_place
  type(cohort), pointer :: nextc
  type(cohort), pointer :: nextnextc
  real :: hite_threshold
  logical :: fusion_test
  real :: newn
  real :: total_lai
  integer :: count_cohorts
  real :: dbh2h
  real :: dbh2bl

  cc => cp%tallest
  if(.not.associated(cc))return  ! return if there aren't any cohorts

  fusion_took_place = 0

  do while(.not.associated(cc,target=cp%shortest))
     nextc => cc%shorter
     nextnextc => nextc%shorter

     do while(.not.associated(cc,target=nextc))

        ! get fusion height threshold
        if(rho(nextc%pft) == 0.0)then

           hite_threshold = b1Ht(nextc%pft)

        else

           hite_threshold = dbh2h(nextc%pft, max_dbh(nextc%pft))

        endif

        ! test for similarity
        if(nextc%hite < (0.95 * hite_threshold ))then

           fusion_test = (abs(cc%hite - nextc%hite) < fusetol_h)

        else

           fusion_test = (abs(cc%dbh - nextc%dbh)/(0.5*(cc%dbh + nextc%dbh))  &
                < fusetol)

        endif

        if(fusion_test)then
        ! Cohorts have a similar size

           newn = cc%nplant + nextc%nplant

           total_lai = (nextc%nplant * dbh2bl(nextc%dbh,nextc%pft)  &
                + cc%nplant * dbh2bl(cc%dbh,cc%pft)) * sla(cc%pft)

           if( cc%pft == nextc%pft .and.   &  ! cohorts are the same PFT and
                total_lai < lai_fuse_tol .and.  &  ! LAI won't be too big and
                cc%first_census == nextc%first_census .and. &  ! won't mess
                cc%new_recruit_flag == nextc%new_recruit_flag)then ! up output

              fusion_took_place = 1

              call fuse_2_cohorts(cc, nextc, newn)

              nextc%taller%shorter => nextnextc  
              if (.not.associated(nextc%shorter))then
                 cp%shortest => nextc%taller
              else
                 nextnextc%taller => nextc%taller
              endif

              deallocate(nextc)
              nullify(nextc)

           endif
        endif

        if(associated(nextnextc))then
           nextc => nextnextc
           nextnextc => nextc%shorter
        else 
           nextc => cc
        endif
     enddo
     if(.not.associated(cc,target=cp%shortest)) cc => cc%shorter
  enddo

  if(fusion_took_place == 1)call sort_cohorts(cp)

  cp%cohort_count = count_cohorts(cp)

  return
end subroutine fuse_cohorts

!=============================================================

subroutine split_cohorts(cp)

  use ed_structure_defs
  use pft_coms, only: q, qsw, sla
  use fusion_fission_coms, only: lai_tol
  use lphys_interface, only: insert_cohort
  implicit none
  
  type(patch),  target  :: cp
  type(cohort), pointer :: cc
  type(cohort), pointer :: copyc
  real, parameter :: epsilon=0.0001
  integer :: i
  integer :: count_cohorts
  real :: dbh2h
  real :: slai
  real :: bd2dbh

  cc => cp%tallest
  do while(associated(cc))

     slai =  cc%nplant * cc%balive * cp%siteptr%green_leaf_factor(cc%pft)   &
          /(1.0 + q(cc%pft) + qsw(cc%pft) * cc%hite) * sla(cc%pft)

     if(slai > lai_tol)then

        ! Half the densities
        cc%nplant = cc%nplant * 0.5
        cc%lai = cc%lai * 0.5
        ! cc%hcapveg = cc%hcapveg * 0.5
        cc%veg_water = cc%veg_water * 0.5
        cc%omean_gpp = cc%omean_gpp * 0.5
        cc%omean_leaf_resp = cc%omean_leaf_resp * 0.5
        cc%omean_root_resp = cc%omean_root_resp * 0.5
        cc%growth_respiration = cc%growth_respiration * 0.5
        cc%storage_respiration = cc%storage_respiration * 0.5
        cc%vleaf_respiration = cc%vleaf_respiration * 0.5
        cc%Psi_open = cc%Psi_open * 0.5

        ! Allocate memory and copy cohort
        nullify(copyc)
        allocate(copyc)
        call copy_cohort(cc, copyc)

        ! Tweak the heights and DBHs
        cc%bdead = cc%bdead - epsilon
        cc%dbh = bd2dbh(cc%pft, cc%bdead)
        cc%hite = dbh2h(cc%pft, cc%dbh)

        copyc%bdead = copyc%bdead + epsilon
        copyc%dbh = bd2dbh(copyc%pft, copyc%bdead)
        copyc%hite = dbh2h(copyc%pft, copyc%dbh)

        ! Set pointers
        copyc%patchptr => cp
        copyc%siteptr => cp%siteptr
        call insert_cohort(copyc,cp)
     endif
     cc => cc%shorter
  enddo

  cp%cohort_count = count_cohorts(cp)

  return
end subroutine split_cohorts

!===============================================

subroutine insert_cohort(pcc,cp)

  use ed_structure_defs

  implicit none

  type(patch)           :: cp
  type(cohort), target  :: pcc
  type(cohort), pointer :: icohort
  type(cohort), pointer :: ptallest
  type(cohort), pointer :: pshortest
  type(cohort), pointer :: pshortest_save
  type(cohort), pointer :: tallptr
  type(cohort), pointer :: shortptr

  real :: tsp

  icohort => pcc
  tsp = icohort%hite

  pshortest => cp%shortest
  ptallest => cp%tallest

  pshortest_save => cp%shortest
 
  if(associated(cp%shortest))then
     insert_shortest: do while(associated(pshortest))
        if(pshortest%hite .lt. tsp)then
           pshortest => pshortest%taller
        else
           exit insert_shortest
        endif
     enddo insert_shortest
  endif

  tallptr => pshortest
  pshortest => pshortest_save

  if(.not.associated(tallptr))then
     shortptr => ptallest
     ptallest => icohort
  else
     shortptr => tallptr%shorter
     tallptr%shorter => icohort
  endif

  cp%tallest => ptallest

  if(.not.associated(shortptr))then
     pshortest => icohort
  else
     shortptr%taller => icohort
  endif

  cp%shortest => pshortest

  pcc%taller => tallptr
  pcc%shorter => shortptr

  return
end subroutine insert_cohort

!===================================================================

subroutine copy_cohort(src,dest)

  use ed_structure_defs

  implicit none

  type(cohort) :: src
  type(cohort) :: dest

  dest%pft = src%pft

  dest%nplant = src%nplant

  dest%hite = src%hite

  dest%dbh = src%dbh

  dest%bdead = src%bdead

  dest%bleaf = src%bleaf

  dest%phenology_status = src%phenology_status

  dest%balive = src%balive

  dest%lai = src%lai

  dest%bstorage = src%bstorage

  dest%hcapveg = src%hcapveg

  dest%veg_temp = src%veg_temp

  dest%veg_water = src%veg_water

  dest%omean_gpp = src%omean_gpp

  dest%omean_leaf_resp = src%omean_leaf_resp

  dest%omean_root_resp = src%omean_root_resp

  dest%growth_respiration = src%growth_respiration

  dest%storage_respiration = src%storage_respiration

  dest%vleaf_respiration = src%vleaf_respiration

  dest%fsn = src%fsn

  allocate(dest%old_stoma_data)
  dest%old_stoma_data = src%old_stoma_data

  dest%Psi_open = src%Psi_open

  dest%cb(1:13) = src%cb(1:13)

  dest%cb_max(1:13) = src%cb_max(1:13)

  dest%cbr_bar = src%cbr_bar

  dest%krdepth = src%krdepth

  dest%first_census = src%first_census

  dest%new_recruit_flag = src%new_recruit_flag

  return
end subroutine copy_cohort

!=======================================================================

subroutine fuse_patches(cs)

  use ed_structure_defs
  use fusion_fission_coms, only: ff_ndbh, ntol, profile_tol
  use mem_ed, only: n_pft
  use lphys_interface, only: fuse_2_patches

  implicit none

  type(site) :: cs
  real :: tolerance_mult
  type(patch), pointer :: cp
  type(patch), pointer :: nextp
  integer :: istop
  type(patch), pointer :: tpp
  type(patch), pointer :: tmpptr
  integer :: i
  integer :: j
  real :: norm

  integer :: count_cohorts,npatches
  integer, parameter :: maxpatches = 15

  ! ALGORITHM
  ! set all fusion flags to true
  ! create size profiles
  ! goto every patch
  ! find next older patch with same dist_type
  ! check fusion criterion
  ! if within criterion, fuse, otherwise, skip

  tolerance_mult = 1.0
  max_patch: do
     
     !  loop over patches and create pft size profiles
     cp => cs%youngest_patch
     do while(associated(cp))
        call patch_pft_size_profile(cp,ff_ndbh)
        cp%fuse_flag = 1
        cp => cp%older
     enddo
     
     cp => cs%youngest_patch
     do while(associated(cp) .and. associated(cp%older))

        ! find fusion candidate for a given patch
        nextp => cp%older
        istop=0
        nullify(tpp)
        do while((associated(nextp)) .and. (istop .eq. 0))
           if(cp%dist_type == nextp%dist_type .and.   &
                cp%plantation == nextp%plantation)then
              istop = 1
              tpp => nextp
           endif
           nextp => nextp%older
        enddo
        
        if(associated(tpp))then ! found a fusion candidate?

           do i=1, n_pft          ! loop over pft 

              do j=1, ff_ndbh !      loop over hgt bins

                 if( (cp%pft_density_profile(i,j) >   &
                      (tolerance_mult * ntol)) .or.  &
                      (tpp%pft_density_profile(i,j) >   &
                      (tolerance_mult * ntol)))then

                    norm = abs(cp%pft_density_profile(i,j) -   &
                         tpp%pft_density_profile(i,j)) / (0.5 *   &
                         (cp%pft_density_profile(i,j) +   &
                         tpp%pft_density_profile(i,j)))

                    if(norm > profile_tol) cp%fuse_flag = 0   ! reject

                 endif
              enddo
           enddo
           
           ! fusion
           if(cp%fuse_flag == 1)then
              tmpptr => cp%older   !  tmpptr is needed bc f2p frees cp
              call fuse_2_patches(cp,tpp)
              cp => tmpptr
           else
              cp => cp%older
           endif
        else
           cp => cp%older
        endif
     enddo
     
     cp => cs%oldest_patch
     npatches = 0
     do while(associated(cp))
        npatches = npatches + 1
        cp%cohort_count = count_cohorts(cp)
        cp => cp%younger
     enddo
     
     if(npatches <= maxpatches)exit max_patch
     tolerance_mult = tolerance_mult * 1.1
  enddo max_patch
  
  return

end subroutine fuse_patches

!=====================================================================

subroutine patch_pft_size_profile(cp, nbins)

  use ed_structure_defs
  use fusion_fission_coms, only: maxdbh
  use mem_ed, only: n_pft

  implicit none

  type(patch)         :: cp
  integer, intent(in) :: nbins
  type(cohort), pointer :: cc
  real :: dh
  real :: rmin
  real :: rmax
  integer :: j

  dh = maxdbh / nbins

  ! initialize bins
  cp%pft_density_profile(1:n_pft,1:nbins)=0.0

  ! Loop over cohorts
  cc => cp%shortest
  do while(associated(cc)) 

     ! Loop over bins
     do j=1,nbins
        rmin = (j - 1) * dh
        rmax = j * dh
        if((cc%dbh >= rmin) .and. (cc%dbh < rmax))then
           cp%pft_density_profile(cc%pft,j) =   &
                cp%pft_density_profile(cc%pft,j) + cc%nplant
        endif
     enddo

     ! deal with largest dbh bin
     j = nbins
     rmin = j * dh
     if(cc%dbh > rmin)then
        cp%pft_density_profile(cc%pft,j) = cp%pft_density_profile(cc%pft,j)  &
             + cc%nplant
     endif

     cc => cc%taller
  enddo
  
  return

end subroutine patch_pft_size_profile

!====================================================================
subroutine fuse_2_patches(p1, p2)
  
  use ed_structure_defs
  use leaf_coms, only: nzg, slcpd, nzs
  use fusion_fission_coms, only: ff_ndbh
  use lphys_interface, only: insert_cohort
  implicit none

  type(patch),  target  :: p1
  type(patch),  target  :: p2
  type(patch),  pointer :: dp
  type(patch),  pointer :: rp
  type(cohort), pointer :: cc
  real :: newarea
  real :: newareai
  integer :: k
  integer :: count_cohorts
  type(cohort), pointer :: nextc

  dp => p1
  rp => p2

  newarea = rp%area + dp%area
  newareai = 1.0 / newarea

  rp%age = (rp%age * rp%area + dp%age * dp%area) * newareai

  rp%fast_soil_C = (rp%fast_soil_C * rp%area + dp%fast_soil_C * dp%area) *   &
       newareai

  rp%slow_soil_C = (rp%slow_soil_C * rp%area + dp%slow_soil_C * dp%area) *   &
       newareai

  rp%structural_soil_C = (rp%structural_soil_C * rp%area +   &
       dp%structural_soil_C * dp%area) * newareai

  rp%structural_soil_L = (rp%structural_soil_L * rp%area +   &
       dp%structural_soil_L * dp%area) * newareai

  rp%mineralized_soil_N = (rp%mineralized_soil_N * rp%area +   &
       dp%mineralized_soil_N * dp%area) * newareai

  rp%fast_soil_N = (rp%fast_soil_N * rp%area + dp%fast_soil_N * dp%area) *   &
       newareai

  rp%sum_dgd = (rp%sum_dgd * rp%area + dp%sum_dgd * dp%area) * newareai

  rp%sum_chd = (rp%sum_chd * rp%area + dp%sum_chd * dp%area) * newareai

  rp%can_temp = (rp%can_temp * rp%area + dp%can_temp * dp%area) * newareai

  rp%can_shv = (rp%can_shv * rp%area + dp%can_shv * dp%area) * newareai

  rp%can_depth = (rp%can_depth * rp%area + dp%can_depth * dp%area) * newareai

  rp%sfcwater_energy(1:nzs) = (rp%sfcwater_energy(1:nzs) *   &
       rp%sfcwater_mass(1:nzs) * rp%area + dp%sfcwater_energy(1:nzs) *   &
       dp%sfcwater_mass(1:nzs) * dp%area) * newareai

  rp%sfcwater_mass(1:nzs) = (rp%sfcwater_mass(1:nzs) * rp%area +  &
       dp%sfcwater_mass(1:nzs) * dp%area) * newareai

  rp%sfcwater_depth(1:nzs) = (rp%sfcwater_depth(1:nzs) * rp%area +  &
       dp%sfcwater_depth(1:nzs) * dp%area) * newareai

  rp%soil_energy(1:nzg) = (rp%soil_energy(1:nzg) * rp%area +   &
       dp%soil_energy(1:nzg) * dp%area) * newareai
  
  rp%soil_water(1:nzg) = (rp%soil_water(1:nzg) * rp%area +   &
       dp%soil_water(1:nzg) * dp%area) * newareai

  call new_patch_sfc_props(rp)

  rp%rough = (rp%rough * rp%area + dp%rough * dp%area) * newareai

  rp%omean_rh = (rp%omean_rh * rp%area + dp%omean_rh * dp%area) * newareai
  rp%omean_runoff = (rp%omean_runoff * rp%area + dp%omean_runoff *   &
       dp%area) * newareai
  rp%omean_wflux = (rp%omean_wflux * rp%area + dp%omean_wflux *   &
       dp%area) * newareai
  rp%omean_latflux = (rp%omean_latflux * rp%area + dp%omean_latflux *   &
       dp%area) * newareai
  rp%omean_qrunoff = (rp%omean_qrunoff * rp%area + dp%omean_qrunoff *   &
       dp%area) * newareai
  rp%omean_hflux = (rp%omean_hflux * rp%area + dp%omean_hflux *   &
       dp%area) * newareai

  rp%dmean_A_decomp = (rp%dmean_A_decomp * rp%area + dp%dmean_A_decomp *   &
       dp%area) * newareai

  rp%dmean_Af_decomp = (rp%dmean_Af_decomp * rp%area + dp%dmean_Af_decomp *   &
       dp%area) * newareai

  rp%repro(1:n_pft) = (rp%repro(1:n_pft) * rp%area + dp%repro(1:n_pft) *   &
       dp%area) * newareai

  rp%fsc_in = (rp%fsc_in * rp%area + dp%fsc_in * dp%area) * newareai

  rp%ssc_in = (rp%ssc_in * rp%area + dp%ssc_in * dp%area) * newareai

  rp%ssl_in = (rp%ssl_in * rp%area + dp%ssl_in * dp%area) * newareai

  rp%fsn_in = (rp%fsn_in * rp%area + dp%fsn_in * dp%area) * newareai

  rp%total_plant_nitrogen_uptake = (rp%total_plant_nitrogen_uptake *   &
       rp%area + dp%total_plant_nitrogen_uptake * dp%area) * newareai

  ! adjust densities of cohorts in recipient patch
  cc => rp%shortest
  do while(associated(cc))
     cc%nplant = cc%nplant * rp%area * newareai
     cc%lai = cc%lai * rp%area * newareai
     !     cc%hcapveg = cc%hcapveg * rp%area * newareai
     cc%veg_water = cc%veg_water * rp%area * newareai
     cc%omean_gpp = cc%omean_gpp * rp%area * newareai
     cc%omean_leaf_resp = cc%omean_leaf_resp * rp%area * newareai
     cc%omean_root_resp = cc%omean_root_resp * rp%area * newareai
     cc%growth_respiration = cc%growth_respiration * rp%area * newareai
     cc%storage_respiration = cc%storage_respiration * rp%area * newareai
     cc%vleaf_respiration = cc%vleaf_respiration * rp%area * newareai
     cc%Psi_open = cc%Psi_open * rp%area * newareai
     cc => cc%taller
  enddo
  
  ! adjust densities of cohorts in recipient patch
  cc => dp%shortest
  do while(associated(cc))
     cc%nplant = cc%nplant * dp%area * newareai
     cc%lai = cc%lai * dp%area * newareai
     !     cc%hcapveg = cc%hcapveg * dp%area * newareai
     cc%veg_water = cc%veg_water * dp%area * newareai
     cc%omean_gpp = cc%omean_gpp * dp%area * newareai
     cc%omean_leaf_resp = cc%omean_leaf_resp * dp%area * newareai
     cc%omean_root_resp = cc%omean_root_resp * dp%area * newareai
     cc%growth_respiration = cc%growth_respiration * dp%area * newareai
     cc%storage_respiration = cc%storage_respiration * dp%area * newareai
     cc%vleaf_respiration = cc%vleaf_respiration * dp%area * newareai
     cc%Psi_open = cc%Psi_open * dp%area * newareai
     cc => cc%taller
  enddo
  
  ! insert donor cohorts into recipient patch
  nullify(nextc)
  cc => dp%shortest
  if(associated(cc)) nextc => cc%taller
  do while(associated(dp%shortest))
     call insert_cohort(cc,rp)
     cc%patchptr => rp
     cc => nextc
     dp%shortest => cc
     if(associated(cc)) nextc => cc%taller
  enddo

  ! Reset veg_height, LAI, etc.
  call update_patch_derived_props(rp)

  ! update size profile within patch
  call patch_pft_size_profile(rp, ff_ndbh)

  rp%cohort_count = count_cohorts(rp)
  rp%area = newarea

  ! update pointers
  if(.not.associated(dp,target=dp%siteptr%youngest_patch))then
     ! donor not youngest
     dp%younger%older => dp%older
  else 
     dp%siteptr%youngest_patch => dp%older
  endif
  if(.not.associated(dp,target=dp%siteptr%oldest_patch))then
     ! donor not oldest
     dp%older%younger => dp%younger
  else
     dp%siteptr%oldest_patch => dp%younger
  endif

  ! free donor patch
  deallocate(dp)

  return
end subroutine fuse_2_patches

!=====================================================================

subroutine terminate_patches(cs)
  
  use ed_structure_defs
  implicit none

  type(site)            :: cs
  type(patch),  pointer :: cp
  type(patch),  pointer :: nextp
  type(cohort), pointer :: cc
  type(cohort), pointer :: nextc

  cp => cs%oldest_patch
  do while(associated(cp))
     nextp => cp%younger

     if(cp%area <= 0.0)then
        ! terminate this patch
     
        if(.not.associated(cp%older))then
           ! it is the oldest.  There is necessarily a younger patch because
           ! areas must sum to 1.
           cs%oldest_patch => cp%younger
           nullify(cs%oldest_patch%older)
        elseif(.not.associated(cp%younger))then
           cs%youngest_patch => cp%older
           nullify(cs%youngest_patch%younger)
        else
           cp%older%younger => cp%younger
           cp%younger%older => cp%older
        endif

        cc => cp%tallest
        do while(associated(cc))
           nextc => cc%shorter
           deallocate(cc)
           cc => nextc
        enddo

        deallocate(cp)

     endif

     cp => nextp
  enddo

  return
end subroutine terminate_patches

!===================================================================

subroutine fuse_2_cohorts(cc, nextc, newn)
  
  use ed_structure_defs
  use pft_coms, only: q, qsw, sla
  use mem_leaf, only: land

  implicit none

  type(cohort) :: cc
  type(cohort) :: nextc
  real, intent(in) :: newn
  real :: newni
  real :: bd2dbh
  real :: dbh2h
  real :: cb_max
  real :: root_depth
  real :: calc_root_depth
  integer, external :: assign_root_depth

  newni = 1.0 / newn
  
  ! Conserve carbon by calculating bdead first.
  cc%bdead = (cc%nplant * cc%bdead + nextc%nplant * nextc%bdead) * newni

  ! Then get dbh and hite from bdead.
  cc%dbh = bd2dbh(cc%pft, cc%bdead)
  cc%hite = dbh2h(cc%pft, cc%dbh)

  ! keep the phenology_status of cc and conserve carbon to get balive.
  cc%balive = (cc%nplant * cc%balive + nextc%nplant * nextc%balive) *newni

  ! Update bleaf and lai.
  if(cc%phenology_status < 2)then
     cc%bleaf = cc%siteptr%green_leaf_factor(cc%pft) * cc%balive /   &
          (1.0 + q(cc%pft) + cc%hite * qsw(cc%pft))
     cc%lai = cc%bleaf * sla(cc%pft) * newn
  else
     cc%lai = 0.0
  endif

  cc%bstorage = (cc%nplant * cc%bstorage + nextc%nplant *  &
       nextc%bstorage) * newni

  ! Do our best to conserve energy.  May be problematic
  ! if cohorts have different phenology_status codes.
  if(cc%phenology_status < 2 .and. nextc%phenology_status < 2)then

     cc%veg_temp = (cc%veg_temp * cc%nplant * cc%hcapveg +  &
          nextc%veg_temp * nextc%nplant * nextc%hcapveg) /   &
          (cc%hcapveg * cc%nplant + nextc%hcapveg * nextc%nplant)

     cc%veg_water = (cc%veg_water * cc%nplant +   &
          nextc%veg_water * nextc%nplant) * newni

     ! Remember, phenology_status is determined solely by cc.
     ! So if cc didn't have leaves before, it doesn't have
     ! any now.

  elseif(cc%phenology_status < 2)then

     cc%veg_water = cc%veg_water * cc%nplant * newni

  endif

  ! No need to modify hcapveg until the new implementation.
  
  cc%omean_gpp = (cc%omean_gpp * cc%nplant + nextc%omean_gpp *   &
       nextc%nplant) * newni

  cc%omean_leaf_resp = (cc%omean_leaf_resp * cc%nplant +   &
       nextc%omean_leaf_resp * nextc%nplant) * newni

  cc%omean_root_resp = (cc%omean_root_resp * cc%nplant +   &
       nextc%omean_root_resp * nextc%nplant) * newni
  
  cc%growth_respiration = (cc%growth_respiration * cc%nplant +  &
       nextc%growth_respiration * nextc%nplant) * newni
  
  cc%storage_respiration = (cc%storage_respiration * cc%nplant +  &
       nextc%storage_respiration * nextc%nplant) * newni
  
  cc%vleaf_respiration = (cc%vleaf_respiration * cc%nplant +  &
       nextc%vleaf_respiration * nextc%nplant) * newni
  
  cc%fsn = (cc%fsn * cc%nplant + nextc%fsn * nextc%nplant) * newni

  cc%Psi_open = (cc%Psi_open * cc%nplant + nextc%Psi_open *   &
       nextc%nplant) * newni

  cc%cb(1:13) = (cc%cb(1:13) * cc%nplant + nextc%nplant *   &
       nextc%cb(1:13)) * newni

  cc%cb_max(1:13) = (cc%cb_max(1:13) * cc%nplant +   &
       nextc%nplant * nextc%cb_max(1:13)) * newni

  cb_max = sum(cc%cb_max(1:12))
  if(cb_max > 0.0)then
     cc%cbr_bar = sum(cc%cb(1:12)) / cb_max
  else
     cc%cbr_bar = 0.0
  endif
  
  root_depth = calc_root_depth(cc%hite, cc%dbh, cc%pft)
  cc%krdepth = assign_root_depth(root_depth, cc%patchptr,   &
       land%lsl(cc%siteptr%iland))

  cc%nplant = newn

  return
end subroutine fuse_2_cohorts
