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
subroutine ed_history_write()

  use misc_coms, only: current_time, hfilepref, time8
  use ed_structure_defs
  use mem_leaf,  only: first_site
  use leaf_coms, only: nzg
  use max_dims,  only: pathlen

  implicit none

  character(pathlen+40) :: fname
  integer :: nsites
  type(site), pointer :: cs
  type(patch), pointer :: cp
  type(cohort), pointer :: cc
  integer :: k
  integer :: idist
  integer :: jdist
  integer :: npatches
  integer :: idbh
  integer :: ipft
  integer :: int_time
  integer :: nhours
  integer :: nminutes
  integer :: nseconds

  nhours = int(current_time%time / 3600.0)
  nminutes = int((current_time%time - 3600.0 * nhours) / 60.0)
  nseconds = current_time%time - 3600.0 * nhours - 60.0 * nminutes
  int_time = 100 * (100 * nhours + nminutes) + nseconds
  write(fname,'(a,i4.4,a,i2.2,a,i2.2,a,i6.6,a)')trim(hfilepref)//'-',  &
       current_time%year,'-',current_time%month,'-',current_time%date,  &
       '-',int_time,'-ED-RESTART.dat'

  open(12,file=trim(fname),form='unformatted',status='replace')

  write(12)current_time%year, current_time%month, current_time%date,  &
       current_time%time

  write(12)time8

  nsites = 0
  cs => first_site
  do while(associated(cs))
     nsites = nsites + 1
     cs => cs%next_site
  enddo

  write(12)nsites
  cs => first_site
  do while(associated(cs))
     write(12)cs%iland, cs%lat, cs%lon

     write(12)cs%agri_stocking_density, cs%plantation_stocking_density,  &
          cs%agri_stocking_pft, cs%plantation_stocking_pft

     do ipft = 1,n_pft
        do idbh = 1,n_dbh
           write(12)  &
                cs%basal_area_growth(ipft,idbh),   &
                cs%basal_area_mort(ipft,idbh),     &
                cs%basal_area_cut(ipft,idbh),      &
                cs%basal_area_recruit(ipft,idbh),  &
                cs%agb_growth(ipft,idbh),          &
                cs%agb_mort(ipft,idbh),            &
                cs%agb_cut(ipft,idbh),             &
                cs%agb_recruit(ipft,idbh)
        enddo
     enddo

     write(12)cs%min_monthly_temp, cs%mmean_gpp, cs%mmean_plresp,  &
          cs%mmean_rh, cs%mmean_nep

     write(12)cs%lambda_fire(1:12), cs%nat_disturbance_rate,  &
          cs%nat_dist_type

     do idist = 1,n_dist_types
        do jdist = 1,n_dist_types
           write(12)cs%disturbance_rates(idist,jdist),  &
                cs%disturbance_memory(idist,jdist)
        enddo
     enddo

     write(12)cs%primary_harvest_memory,   &
          cs%secondary_harvest_memory

     npatches = 0
     cp => cs%oldest_patch
     do while(associated(cp))
        npatches = npatches + 1
        cp => cp%younger
     enddo
     write(12)npatches

     cp => cs%oldest_patch
     do while(associated(cp))
        write(12)  &
             cp%dist_type,  &
             cp%age,   &
             cp%area,   &
             cp%fast_soil_C,  &
             cp%slow_soil_C,  &
             cp%structural_soil_C,  &
             cp%structural_soil_L,  &
             cp%mineralized_soil_N,  &
             cp%fast_soil_N,  &
             cp%sum_dgd,  &
             cp%sum_chd,  &
             cp%plantation,  &
             cp%can_temp,  &
             cp%can_shv,  &
             cp%can_depth,  &
             cp%nlev_sfcwater, &
             cp%dmean_A_decomp, &
             cp%dmean_Af_decomp, &
             cp%avg_daily_temp

        do k = 1,cp%nlev_sfcwater
           write(12)  &
                cp%sfcwater_mass(k),  &
                cp%sfcwater_energy(k),  &
                cp%sfcwater_depth(k)
        enddo

        do k = 1,nzg
           write(12)  &
                cp%ntext_soil(k),  &
                cp%soil_energy(k),  &
                cp%soil_water(k)
        enddo

        write(12)cp%rough, cp%cohort_count

        do ipft = 1,n_pft
           write(12)cp%repro(ipft)
        enddo

        write(12)cp%ground_shv, cp%surface_ssh

        cc => cp%tallest
        do while(associated(cc))
           write(12)  &
                cc%pft,  &
                cc%nplant,  &
                cc%hite,  &
                cc%dbh,  &
                cc%bdead,  &
                cc%bleaf,  &
                cc%phenology_status,  &
                cc%balive,  &
                cc%bstorage,  &
                cc%hcapveg,  &
                cc%veg_temp,  &
                cc%veg_water,  &
                cc%fsn,  &
                cc%monthly_dndt,  &
                cc%Psi_open,  &
                cc%cb(1:13),  &
                cc%cb_max(1:13),  &
                cc%first_census,  &
                cc%new_recruit_flag,  &
                cc%lai,  &
                cc%dmean_gpp, &
                cc%dmean_gpp_pot, &
                cc%dmean_gpp_max, &
                cc%dmean_leaf_resp, &
                cc%dmean_root_resp

           cc => cc%shorter
        enddo

        cp => cp%younger
     enddo

     cs => cs%next_site
  enddo


  close(12)
  return
end subroutine ed_history_write

!====================================================================

subroutine ed_history_read()

  use ed_structure_defs
  use misc_coms, only: hfilepref, current_time, time8
  use mem_leaf, only: land, first_site
  use leaf_coms, only: nzg, nzs, slcpd
  use pft_coms, only: sla
  use ed_options, only: istoma_scheme, ed_hfilin
  use max_dims, only: pathlen

  implicit none

  character(pathlen) :: fname
  type(site), pointer :: cs
  type(patch), pointer :: new_patch
  type(cohort), pointer :: new_cohort
  integer :: ipft
  integer :: idbh
  integer :: npatches
  integer :: k
  integer :: idist
  integer :: jdist
  integer :: ipatch
  integer :: icohort
  real :: root_depth
  real, external :: calc_root_depth
  integer, external :: assign_root_depth
  real :: cb_act
  real :: cb_max
  real :: lat
  real :: lon
  integer :: iland
  type(site), pointer :: tcs
  logical :: l1

  fname = trim(ed_hfilin)

  inquire(file=trim(fname),exist=l1)
  if(.not.l1)then
     print*,'File:'
     print*,trim(fname)
     print*,'does not exist.  Specify ED_HFILIN properly in OLAMIN'
     print*
     print*,'Message brought to you by subroutine ed_history_read in ed_history.f90'
     stop
  endif

  open(12,file=trim(fname),form='unformatted',status='old')

  ! Read the time stamp
  read(12)current_time%year, current_time%month, current_time%date,  &
       current_time%time

  read(12)time8

  ! Skip over the number of sites
  read(12)

  cs => first_site

  cs%omean_precip = 0.0
  cs%omean_qprecip = 0.0
  cs%omean_netrad = 0.0

  do while(associated(cs))

     ! Skip iland, lat, lon
     read(12)iland,lat,lon
     do while(cs%iland /= iland)
        if(associated(cs,target=first_site))then
           first_site => cs%next_site
        else
           tcs => first_site
           find_last_site: do while(associated(tcs))
              if(associated(tcs%next_site,target=cs))then
                 tcs%next_site => cs%next_site
                 exit find_last_site
              endif
              tcs => tcs%next_site
           enddo find_last_site
        endif
        land%ed_flag(cs%iland) = 0
        tcs => cs%next_site
        deallocate(cs)
        cs => tcs
     enddo

     read(12)cs%agri_stocking_density, cs%plantation_stocking_density,  &
          cs%agri_stocking_pft, cs%plantation_stocking_pft

     do ipft = 1,n_pft
        do idbh = 1,n_dbh
           read(12)  &
                cs%basal_area_growth(ipft,idbh),   &
                cs%basal_area_mort(ipft,idbh),     &
                cs%basal_area_cut(ipft,idbh),      &
                cs%basal_area_recruit(ipft,idbh),  &
                cs%agb_growth(ipft,idbh),          &
                cs%agb_mort(ipft,idbh),            &
                cs%agb_cut(ipft,idbh),             &
                cs%agb_recruit(ipft,idbh)
        enddo
     enddo

     read(12)cs%min_monthly_temp, cs%mmean_gpp, cs%mmean_plresp,  &
          cs%mmean_rh, cs%mmean_nep

     read(12)cs%lambda_fire(1:12), cs%nat_disturbance_rate,  &
          cs%nat_dist_type

     do idist = 1,n_dist_types
        do jdist = 1,n_dist_types
           read(12)cs%disturbance_rates(idist,jdist),  &
                cs%disturbance_memory(idist,jdist)
        enddo
     enddo

     read(12)cs%primary_harvest_memory,   &
          cs%secondary_harvest_memory

     nullify(cs%oldest_patch)
     nullify(cs%youngest_patch)

     read(12)npatches

     do ipatch = 1,npatches

        nullify(new_patch)
        allocate(new_patch)
        if(.not.associated(cs%oldest_patch))then
           cs%oldest_patch => new_patch
           nullify(new_patch%older)
        else
           new_patch%older => cs%youngest_patch
           cs%youngest_patch%younger => new_patch
        endif
        cs%youngest_patch => new_patch
        nullify(new_patch%younger)
        new_patch%siteptr => cs
        nullify(new_patch%tallest)
        nullify(new_patch%shortest)

        read(12)  &
             new_patch%dist_type,  &
             new_patch%age,   &
             new_patch%area,   &
             new_patch%fast_soil_C,  &
             new_patch%slow_soil_C,  &
             new_patch%structural_soil_C,  &
             new_patch%structural_soil_L,  &
             new_patch%mineralized_soil_N,  &
             new_patch%fast_soil_N,  &
             new_patch%sum_dgd,  &
             new_patch%sum_chd,  &
             new_patch%plantation,  &
             new_patch%can_temp,  &
             new_patch%can_shv,  &
             new_patch%can_depth,  &
             new_patch%nlev_sfcwater, &
             new_patch%dmean_A_decomp, &
             new_patch%dmean_Af_decomp, &
             new_patch%avg_daily_temp

        allocate(new_patch%sfcwater_mass(nzs))
        allocate(new_patch%sfcwater_energy(nzs))
        allocate(new_patch%sfcwater_depth(nzs))
        allocate(new_patch%rshort_s(nzs))
        allocate(new_patch%rshort_s_beam(nzs))
        allocate(new_patch%rshort_s_diffuse(nzs))
        allocate(new_patch%sfcwater_tempk(nzs))
        allocate(new_patch%sfcwater_fracliq(nzs))

        do k = 1,new_patch%nlev_sfcwater
           read(12)  &
                new_patch%sfcwater_mass(k),  &
                new_patch%sfcwater_energy(k),  &
                new_patch%sfcwater_depth(k)
           call qtk(new_patch%sfcwater_energy(k),new_patch%sfcwater_tempk(k), &
                new_patch%sfcwater_fracliq(k))
        enddo
        do k = new_patch%nlev_sfcwater+1,nzs
           new_patch%sfcwater_mass(k) = 0.0
           new_patch%sfcwater_energy(k) = 0.0
           new_patch%sfcwater_depth(k) = 0.0
        enddo

        allocate(new_patch%ntext_soil(nzg))
        allocate(new_patch%soil_energy(nzg))
        allocate(new_patch%soil_water(nzg))
        allocate(new_patch%soil_tempk(nzg))
        allocate(new_patch%soil_fracliq(nzg))

        do k = 1,nzg
           read(12)  &
                new_patch%ntext_soil(k),  &
                new_patch%soil_energy(k),  &
                new_patch%soil_water(k)

           call qwtk(new_patch%soil_energy(k),new_patch%soil_water(k)*1.0e3,  &
                slcpd(new_patch%ntext_soil(k)), new_patch%soil_tempk(k),  &
                new_patch%soil_fracliq(k))
        enddo

        read(12)new_patch%rough, new_patch%cohort_count
        
!        if(istoma_scheme == 1)then
           allocate(new_patch%old_stoma_data_max(n_pft))
           new_patch%old_stoma_data_max(1:n_pft)%recalc = 1
!        endif

        new_patch%omean_rh = 0.0
        new_patch%omean_nep = 0.0
        new_patch%omean_runoff = 0.0
        new_patch%omean_qrunoff = 0.0
        new_patch%omean_wflux = 0.0
        new_patch%omean_latflux = 0.0
        new_patch%omean_hflux = 0.0
        
        allocate(new_patch%A_o_max(n_pft))
        allocate(new_patch%A_c_max(n_pft))

        do ipft = 1,n_pft
           read(12)new_patch%repro(ipft)
        enddo

        read(12)new_patch%ground_shv, new_patch%surface_ssh

        do icohort = 1,new_patch%cohort_count
           nullify(new_cohort)
           allocate(new_cohort)
           new_cohort%siteptr => cs
           new_cohort%patchptr => new_patch
           if(.not.associated(new_patch%tallest))then
              new_patch%tallest => new_cohort
              nullify(new_cohort%taller)
           else
              new_cohort%taller => new_patch%shortest
              new_patch%shortest%shorter => new_cohort
           endif
           new_patch%shortest => new_cohort
           nullify(new_cohort%shorter)

           read(12)  &
                new_cohort%pft,  &
                new_cohort%nplant,  &
                new_cohort%hite,  &
                new_cohort%dbh,  &
                new_cohort%bdead,  &
                new_cohort%bleaf,  &
                new_cohort%phenology_status,  &
                new_cohort%balive,  &
                new_cohort%bstorage,  &
                new_cohort%hcapveg,  &
                new_cohort%veg_temp,  &
                new_cohort%veg_water,  &
                new_cohort%fsn,  &
                new_cohort%monthly_dndt,  &
                new_cohort%Psi_open,  &
                new_cohort%cb(1:13),  &
                new_cohort%cb_max(1:13),  &
                new_cohort%first_census,  &
                new_cohort%new_recruit_flag,  &
                new_cohort%lai,  &
                new_cohort%dmean_gpp, &
                new_cohort%dmean_gpp_pot, &
                new_cohort%dmean_gpp_max, &
                new_cohort%dmean_leaf_resp, &
                new_cohort%dmean_root_resp

           new_cohort%omean_gpp = 0.0
           new_cohort%omean_leaf_resp = 0.0
           new_cohort%omean_root_resp = 0.0

!           if(istoma_scheme == 1)then
              nullify(new_cohort%old_stoma_data)
              allocate(new_cohort%old_stoma_data)
              new_cohort%old_stoma_data%recalc = 1
!           endif

           root_depth = calc_root_depth(new_cohort%hite,  &
                new_cohort%dbh, new_cohort%pft)
           new_cohort%krdepth = assign_root_depth(root_depth, new_patch,  &
                land%lsl(cs%iland))

           cb_act = sum(new_cohort%cb(1:12))
           cb_max = sum(new_cohort%cb_max(1:12))
           if(cb_max > 0.0)then
              new_cohort%cbr_bar = cb_act / cb_max
           else
              new_cohort%cbr_bar = 0.0
           endif

        enddo

        call update_patch_derived_props(new_patch)

     enddo

     call update_site_derived_props(cs, 0)

     allocate(cs%green_leaf_factor(n_pft))
     allocate(cs%leaf_aging_factor(n_pft))
     cs%green_leaf_factor(1:n_pft) = 1.0
     cs%leaf_aging_factor(1:n_pft) = 1.0

     cs => cs%next_site
  enddo

  close(12)

  return

end subroutine ed_history_read
