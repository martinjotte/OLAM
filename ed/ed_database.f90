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
subroutine read_ed_database()
  use ed_options,   only: ied_init_mode, iphen_scheme
  implicit none

  if(ied_init_mode == 0)then

     ! Initialize everything with near-bare ground
     call bare_ground_init()

     ! Initialize the thermal sums
     if(iphen_scheme == 0)call initialize_thermal_sums()

  elseif(ied_init_mode == 1)then

     ! Initialize with complete restart information
     call ed_history_read()

  elseif(ied_init_mode == 2)then

     ! Initialize with ED1-type restart info
     call read_ed1_history_file()

     ! Initialize the thermal sums
     if(iphen_scheme == 0)call initialize_thermal_sums()

  endif

  call landuse_init()

  return
end subroutine read_ed_database


!======================================================================
subroutine bare_ground_init()

  use ed_structure_defs
  use mem_leaf,  only: first_site
  use pft_coms,  only: SLA, q, qsw, hgt_min, include_pft, n_pft

  implicit none

  type(site),   pointer :: cs
  type(patch),  pointer :: new_patch
  type(cohort), pointer :: new_cohort
  real :: h2dbh
  real :: dbh2bd
  real :: dbh2bl
  integer :: ipft

  ! Loop over all sites
  cs => first_site
  do while(associated(cs))

     ! Allocate the only patch
     nullify(new_patch)
     allocate(new_patch)

     ! Define the bare-ground setup
     new_patch%dist_type = 3
     new_patch%age = 0.0
     new_patch%area = 1.0
     new_patch%fast_soil_C = 0.2
     new_patch%slow_soil_C = 0.01
     new_patch%structural_soil_C = 10.0
     new_patch%structural_soil_L = new_patch%structural_soil_C
     new_patch%mineralized_soil_N = 1.0
     new_patch%fast_soil_N = 1.0
     new_patch%sum_dgd = 0.0
     new_patch%sum_chd = 0.0
     new_patch%plantation = 0

     ! Add patch to the linked list
     nullify(new_patch%older)
     nullify(new_patch%younger)
     cs%oldest_patch => new_patch
     cs%youngest_patch => new_patch
     new_patch%siteptr => cs
     nullify(new_patch%tallest)
     nullify(new_patch%shortest)

     ! Add the cohorts
     do ipft = 1,n_pft
        if(include_pft(ipft) == 1)then

           ! Allocate memory
           nullify(new_cohort)
           allocate(new_cohort)

           ! Define the near-bare ground
           new_cohort%pft = ipft
           new_cohort%hite = hgt_min(ipft)
           new_cohort%dbh = h2dbh(new_cohort%hite,ipft)
           new_cohort%bdead = dbh2bd(new_cohort%dbh,new_cohort%hite,ipft)
           new_cohort%bleaf = dbh2bl(new_cohort%dbh,ipft)
           new_cohort%nplant = 0.1
           new_cohort%phenology_status = 0
           new_cohort%balive = new_cohort%bleaf * (1.0 + q(ipft) +  &
                qsw(ipft) * new_cohort%hite)
           new_cohort%lai = new_cohort%bleaf * new_cohort%nplant * SLA(ipft)
           new_cohort%bstorage = 0.0

           ! Set pointers
           new_cohort%siteptr => cs
           new_cohort%patchptr => new_patch
           nullify(new_cohort%shorter)
           if(.not.associated(new_patch%tallest))then
              new_patch%tallest => new_cohort
              nullify(new_cohort%taller)
           else
              new_cohort%taller => new_patch%shortest
              new_patch%shortest%shorter => new_cohort
           endif
           new_patch%shortest => new_cohort

           call init_ed_cohort_vars(new_cohort, new_patch)

        endif
     enddo

     call init_ed_patch_vars(new_patch)

     call init_ed_site_vars(cs)

     cs => cs%next_site
  enddo

  return
end subroutine bare_ground_init

!======================================================================
subroutine single_pft_init()

  use ed_structure_defs
  use mem_leaf, only: first_site
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, n_pft

  implicit none

  type(site),   pointer :: cs
  type(patch),  pointer :: new_patch
  type(cohort), pointer :: new_cohort
  real :: h2dbh
  real :: dbh2bd
  real :: dbh2bl
  integer :: ipft
  integer :: ii

  ! Loop over all sites
  cs => first_site
  do while(associated(cs))

     ! Allocate the only patch
     nullify(new_patch)
     allocate(new_patch)

     ! Define the bare-ground setup
     new_patch%dist_type = 3
     new_patch%age = 0.0
     new_patch%area = 1.0
     new_patch%fast_soil_C = 0.2
     new_patch%slow_soil_C = 0.01
     new_patch%structural_soil_C = 10.0
     new_patch%structural_soil_L = new_patch%structural_soil_C
     new_patch%mineralized_soil_N = 1.0
     new_patch%fast_soil_N = 1.0
     new_patch%sum_dgd = 0.0
     new_patch%sum_chd = 0.0
     new_patch%plantation = 0

     ! Add patch to the linked list
     nullify(new_patch%older)
     nullify(new_patch%younger)
     cs%oldest_patch => new_patch
     cs%youngest_patch => new_patch
     new_patch%siteptr => cs
     nullify(new_patch%tallest)
     nullify(new_patch%shortest)

     ! Add the cohorts
     do ipft = 1,n_pft
        if(include_pft(ipft) == 1)then

           do ii = 1,8
              ! Allocate memory
              nullify(new_cohort)
              allocate(new_cohort)
              
              ! Define the near-bare ground
              new_cohort%pft = ipft
              if(ipft == 1)then
                 new_cohort%hite = hgt_min(ipft)
              else
                 new_cohort%hite = hgt_min(ipft) + 1.0 * ii
              endif
              new_cohort%dbh = h2dbh(new_cohort%hite,ipft)
              new_cohort%bdead = dbh2bd(new_cohort%dbh,new_cohort%hite,ipft)
              new_cohort%bleaf = dbh2bl(new_cohort%dbh,ipft)
              new_cohort%nplant = 0.5 / (new_cohort%bleaf * SLA(ipft))
              new_cohort%phenology_status = 0
              new_cohort%balive = new_cohort%bleaf * (1.0 + q(ipft) +  &
                   qsw(ipft) * new_cohort%hite)
              new_cohort%lai = new_cohort%bleaf * new_cohort%nplant * SLA(ipft)
              new_cohort%bstorage = 0.0
              
              ! Set pointers
              new_cohort%siteptr => cs
              new_cohort%patchptr => new_patch
              nullify(new_cohort%shorter)
              if(.not.associated(new_patch%tallest))then
                 new_patch%tallest => new_cohort
                 nullify(new_cohort%taller)
              else
                 new_cohort%taller => new_patch%shortest
                 new_patch%shortest%shorter => new_cohort
              endif
              new_patch%shortest => new_cohort
              
              call init_ed_cohort_vars(new_cohort, new_patch)
           enddo
        endif
     enddo
     
     call init_ed_patch_vars(new_patch)

     call init_ed_site_vars(cs)

     cs => cs%next_site
  enddo

  return
end subroutine single_pft_init

!======================================================================
subroutine read_ed1_history_file()

  use ed_structure_defs
  use mem_leaf,   only: land, first_site
  use pft_coms,   only: SLA, q, qsw, hgt_min, include_pft, n_pft
  use ed_options, only: ed_hfilin
  use max_dims,   only: pathlen

  implicit none

  type(site), pointer :: cs
  type(patch), pointer :: new_patch
  type(cohort), pointer :: new_cohort
  real :: h2dbh
  real :: dbh2bd
  real :: dbh2bl
  integer :: ipft
  integer :: ii
  real, parameter :: edres=1.0
!  character(len=256), parameter :: ed_fprefix = '/cluster/home/dmm31/ed_inputs/amazon_restart/lu.sa.05.restart.01.'
  character(pathlen+40) :: ed_fname
  real :: flat
  real :: flon
  character(pathlen+40) :: pss_name
  character(pathlen+40) :: css_name
  logical :: restart_exist
  type(site), pointer :: tcs
  integer :: nsites

  character(len=256) :: cdum
  integer :: nwater
  real, dimension(100) :: depth
  real :: time
  character(len=256) :: pname
  integer :: trk
  real :: age
  real :: area
  real :: fsc
  real :: stsc
  real :: stsl
  real :: ssc
  real :: psc
  real :: msn
  real :: fsn
  real, dimension(100) :: water
  integer :: ierr
  real, dimension(12) :: cb
  real, dimension(12) :: cb_max
  integer :: leaves_on
  real :: balive
  real :: avgRg
  real :: bdead
  real :: nplant
  real :: hite
  real :: dbh
  character(len=64) :: cname
  type(patch), pointer :: cp
  real :: area_tot

  real :: area_sum
  real :: patch_lai
  type(cohort), pointer :: cc
  real :: site_lai
  integer :: ncohorts

  ! Loop over all sites
  cs => first_site
  site_loop: do while(associated(cs))
     ! Make file name
     if(cs%lat >= 0.0)then
        flat = edres * int(cs%lat / edres) + 0.5 * edres 
     else
        flat = - edres * int(-cs%lat / edres) - 0.5 * edres
     endif

     if(cs%lon >= 0.0)then
        flon = edres * int(cs%lon / edres) + 0.5 * edres 
     else
        flon = - edres * int(-cs%lon / edres) - 0.5 * edres
     endif

     if(edres == 1.0)then
        if(flat <= -10.0)then
           write(ed_fname,'(a,f6.1)')trim(ed_hfilin)//'lat',flat
        elseif(flat < 0.0 .or. flat >= 10.0)then
           write(ed_fname,'(a,f5.1)')trim(ed_hfilin)//'lat',flat
        else
           write(ed_fname,'(a,f4.1)')trim(ed_hfilin)//'lat',flat
        endif
        if(flon <= -100.0)then
           write(ed_fname,'(a,f7.1)')trim(ed_fname)//'lon',flon
        elseif(flon <= -10.0 .or. flon >= 100.0)then
           write(ed_fname,'(a,f6.1)')trim(ed_fname)//'lon',flon
        elseif(flon < 0.0)then
           write(ed_fname,'(a,f5.1)')trim(ed_fname)//'lon',flon
        elseif(flon < 10.0)then
           write(ed_fname,'(a,f4.1)')trim(ed_fname)//'lon',flon
        else
           write(ed_fname,'(a,f5.1)')trim(ed_fname)//'lon',flon
        endif
     else
        print*,'bad edres; need to set up YOUR edres in ed_database.'
        stop
     endif

     pss_name = trim(ed_fname)//'.pss'
     css_name = trim(ed_fname)//'.css'
     
     ! If file does not exist, deallocate and have LEAF do this site.
     inquire(file=trim(pss_name),exist=restart_exist)
     if(.not.restart_exist)then
        ! Remove site from linked list
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
        ! Put LEAF in charge
        land%ed_flag(cs%iland) = 0
        ! Save pointer to next site
        tcs => cs%next_site
        ! Deallocate
        deallocate(cs)
        ! Set pointer to next site and cycle
        cs => tcs
        cycle site_loop
     endif

     nullify(cs%oldest_patch)
     nullify(cs%youngest_patch)

     ! Open file
     open(12,file=trim(pss_name),form='formatted',status='old')
     read(12,*)  ! skip header
     read(12,*)cdum,nwater
     read(12,*)cdum,depth(1:nwater)
     read(12,*)
     area_tot = 0.0
     read_patches: do
        read(12,*,iostat=ierr)time,pname,trk,age,area,fsc,stsc,stsl,ssc, &
             psc,msn,fsn,water(1:nwater)
        if(ierr /= 0)exit read_patches

        ! Allocate new patch
        nullify(new_patch)
        allocate(new_patch)

        new_patch%dist_type = trk + 1
        new_patch%age = age
        new_patch%area = area
        area_tot = area_tot + area
        new_patch%fast_soil_C = fsc
        new_patch%slow_soil_C = ssc
        new_patch%structural_soil_C = stsc
        new_patch%structural_soil_L = stsl
        new_patch%mineralized_soil_N = msn
        new_patch%fast_soil_N = fsn

        new_patch%pname = trim(pname)

        new_patch%sum_dgd = 0.0
        new_patch%sum_chd = 0.0
        new_patch%plantation = 0

        ! Add patch to the linked list
        if(associated(cs%oldest_patch))then
           cs%youngest_patch%younger => new_patch
           new_patch%older => cs%youngest_patch
        else
           cs%oldest_patch => new_patch
           nullify(new_patch%older)
        endif
        cs%youngest_patch => new_patch
        nullify(new_patch%younger)
        new_patch%siteptr => cs
        nullify(new_patch%tallest)
        nullify(new_patch%shortest)

     enddo read_patches
     close(12)

     ! Renormalize areas
     cp => cs%oldest_patch
     do while(associated(cp))
        cp%area = cp%area / area_tot
        cp => cp%younger
     enddo

     ! Add the cohorts
     open(12,file=trim(css_name),form='formatted',status='old')
     read(12,*)  ! skip header
     read(12,*)  ! skip header
     read_cohorts: do

        read(12,*,iostat=ierr)time,pname,cname,dbh,hite,ipft,nplant,  &
             bdead,balive,avgRg,leaves_on,cb(1:12),cb_max(1:12)
        if(ierr /= 0)exit read_cohorts

        ! Find patch
        cp => cs%oldest_patch
        put_cohort: do while(associated(cp))
           if(trim(cp%pname) == trim(pname))exit put_cohort
           cp => cp%younger
        enddo put_cohort

        ! Allocate memory
        nullify(new_cohort)
        allocate(new_cohort)
        ipft = ipft + 1
        if(ipft >= 5)ipft = ipft - 3
        new_cohort%pft = ipft
        new_cohort%dbh = dbh
        new_cohort%hite = hite
        new_cohort%bdead = bdead

        ! Setting balive to default instead of file value
        new_cohort%bleaf = dbh2bl(new_cohort%dbh,ipft)
        new_cohort%balive = new_cohort%bleaf * (1.0 + q(ipft) +  &
             qsw(ipft) * new_cohort%hite)
        new_cohort%phenology_status = 0
        new_cohort%nplant = nplant / (cp%area * area_tot)
        new_cohort%lai = new_cohort%bleaf * new_cohort%nplant * SLA(ipft)
        new_cohort%bstorage = 0.0
        new_cohort%cb(1:12) = cb(1:12)
        new_cohort%cb(13) = 0.0
        new_cohort%cb_max(1:12) = cb_max(1:12)
        new_cohort%cb_max(13) = 0.0

        ! Set pointers
        new_cohort%siteptr => cs
        new_cohort%patchptr => cp
        
        nullify(new_cohort%shorter)
        if(.not.associated(cp%tallest))then
           cp%tallest => new_cohort
           nullify(new_cohort%taller)
        else
           new_cohort%taller => cp%shortest
           cp%shortest%shorter => new_cohort
        endif
        cp%shortest => new_cohort

     enddo read_cohorts
     close(12)

     cp => cs%oldest_patch
     area_sum = 0.0
     site_lai = 0.0
     ncohorts = 0
     do while(associated(cp))
        area_sum = area_sum + cp%area
        patch_lai = 0.0
        cc => cp%tallest
        do while(associated(cc))
           patch_lai = patch_lai + cc%lai
           ncohorts =ncohorts + 1
           cc => cc%shorter
        enddo
        site_lai = site_lai + cp%area * patch_lai
        cp => cp%younger
     enddo
     print*,cs%lat,cs%lon,site_lai,ncohorts

     ! If there are no cohorts, set some up
     ! Add the cohorts
     if(ncohorts == 0)then
        cp => cs%youngest_patch
        do ipft = 1,n_pft
           if(include_pft(ipft) == 1)then

              ! Allocate memory
              nullify(new_cohort)
              allocate(new_cohort)
              
              ! Define the near-bare ground
              new_cohort%pft = ipft
              new_cohort%hite = hgt_min(ipft)
              new_cohort%dbh = h2dbh(new_cohort%hite,ipft)
              new_cohort%bdead = dbh2bd(new_cohort%dbh,new_cohort%hite,ipft)
              new_cohort%bleaf = dbh2bl(new_cohort%dbh,ipft)
              new_cohort%nplant = 0.1
              new_cohort%phenology_status = 0
              new_cohort%balive = new_cohort%bleaf * (1.0 + q(ipft) +  &
                   qsw(ipft) * new_cohort%hite)
              new_cohort%lai = new_cohort%bleaf * new_cohort%nplant * SLA(ipft)
              new_cohort%bstorage = 0.0
              
              ! Set pointers
              new_cohort%siteptr => cs
              new_cohort%patchptr => cp
              nullify(new_cohort%shorter)
              if(.not.associated(cp%tallest))then
                 cp%tallest => new_cohort
                 nullify(new_cohort%taller)
              else
                 new_cohort%taller => cp%shortest
                 cp%shortest%shorter => new_cohort
              endif
              cp%shortest => new_cohort

           endif
        enddo
     endif
     
     cp => cs%oldest_patch
     do while(associated(cp))
        cc => cp%tallest
        do while(associated(cc))
           call init_ed_cohort_vars(cc,cp)
           cc => cc%shorter
        enddo
        call init_ed_patch_vars(cp)
        cp => cp%younger
     enddo
     call init_ed_site_vars(cs)

     cs => cs%next_site
  enddo site_loop

  return
end subroutine read_ed1_history_file
!====================================================================
subroutine initialize_thermal_sums()

  use ed_structure_defs
  use mem_leaf,  only: first_site
  use misc_coms, only: imonth1

  implicit none

  type(site), pointer :: cs
  real, dimension(12) :: chd_m
  real, dimension(12) :: dgd_m
  real :: chd
  real :: dgd
  type(patch), pointer :: cp

  ! Loop over all sites

  cs => first_site
  do while(associated(cs))

     call read_thermal_sum('chd', cs%lat, cs%lon, chd_m)
     call read_thermal_sum('dgd', cs%lat, cs%lon, dgd_m)

     if(cs%lat >= 0.0)then

        if(imonth1 <= 8)then
           dgd = sum(dgd_m(1:imonth1))
        else
           dgd = 0.0
        endif

        if(imonth1 >= 11)then
           chd = sum(chd_m(imonth1:12))
        elseif(imonth1 <= 6)then
           chd = sum(chd_m(11:12)) + sum(chd_m(1:imonth1))
        else
           chd = 0.0
        endif


     else

        if(imonth1 <= 2)then
           dgd = sum(dgd_m(7:12)) + sum(dgd_m(1:imonth1))
        elseif(imonth1 >= 7)then
           dgd = sum(dgd_m(7:imonth1))
        else
           dgd = 0.0
        endif

        if(imonth1 >= 5)then
           chd = sum(chd_m(5:imonth1))
        else
           chd = 0.0
        endif

     endif

     ! Loop over patches

     cp => cs%oldest_patch
     do while(associated(cp))
        cp%sum_chd = chd
        cp%sum_dgd = dgd
        cp => cp%younger
     enddo

     cs => cs%next_site
  enddo

  return
end subroutine initialize_thermal_sums

!=======================================================================

subroutine read_thermal_sum(type, lat, lon, var_out)

  use misc_coms, only: iyear1, imonth1, idate1, itime1
  use ed_options, only: ed_inputs_dir

  implicit none

  character(len=3), intent(in) :: type
  real, intent(in) :: lat
  real, intent(in) :: lon
  real, dimension(12), intent(out) :: var_out
  character(len=256) :: fname
  logical :: exans
  real :: tlat
  real :: tlon
  real, dimension(12) :: var_current_year
  real, dimension(12) :: var_past_year
  real :: partial_month_fraction
  integer :: ierr
  logical :: l1

!  write(fname,'(a,i4.4,a)')trim(ed_inputs_dir)//'temp.'//type//'.y',  &
!       iyear1, '.dat'
!  inquire(file=trim(fname),exist=exans)
!  if(.not.exans)fname = trim(ed_inputs_dir)//'temp.'//type//'.avg.dat'
  fname = trim(ed_inputs_dir)//type//'/temp.'//type//'.avg.dat'
  inquire(file=trim(fname),exist=l1)
  if(.not.l1)then
     print*,'File ',trim(fname),' does not exist.'
     print*,'Stopping in ed_database.f90'
     stop
  endif
  open(12,file=trim(fname),form='formatted',status='old')
  find_thermal_sum_c: do
     read(12,*,iostat=ierr)tlat,tlon,var_current_year(1:12)
     if(ierr /= 0)then
        ! Fill missing data with zeros.
        var_current_year(1:12) = 0.0
        exit find_thermal_sum_c
     endif
     if(abs(tlat-lat) < 1.0 .and. abs(tlon-lon) < 1.0)exit find_thermal_sum_c
  enddo find_thermal_sum_c
  close(12)

  ! Contribution from partial month
  if(imonth1 /= 2 .and. imonth1 /= 4 .and. imonth1 /= 6 .and.   &
       imonth1 /= 9 .and. imonth1 /= 11)then
     partial_month_fraction = (float(idate1-1) + itime1 * 0.01 / 24.0) / 31.0
  elseif(imonth1 /= 2)then
     partial_month_fraction = (float(idate1-1) + itime1 * 0.01 / 24.0) / 30.0
  elseif(mod(iyear1,4) /= 0)then
     partial_month_fraction = (float(idate1-1) + itime1 * 0.01 / 24.0) / 28.0
  else
     partial_month_fraction = (float(idate1-1) + itime1 * 0.01 / 24.0) / 29.0
  endif

  var_out(imonth1) = partial_month_fraction * var_current_year(imonth1)
  var_out(1:(imonth1-1)) = var_current_year(1:(imonth1-1))

!  write(fname,'(a,i4.4,a)')trim(ed_inputs_dir)//type//'/temp.'//type//'.y',  &
!       iyear1-1, '.dat'
!  inquire(file=trim(fname),exist=exans)
!  if(.not.exans)fname = trim(ed_inputs_dir)//'temp.'//type//'.avg.dat'
  fname = trim(ed_inputs_dir)//type//'/temp.'//type//'.avg.dat'
  inquire(file=trim(fname),exist=l1)
  if(.not.l1)then
     print*,'File ',trim(fname),' does not exist.'
     print*,'Stopping in ed_database.f90'
     stop
  endif
  open(12,file=trim(fname),form='formatted',status='old')
  find_thermal_sum_p: do
     read(12,*,iostat=ierr)tlat,tlon,var_past_year(1:12)
     if(ierr /= 0)then
        var_past_year(1:12) = 0.0
        exit find_thermal_sum_p
     endif
     if(abs(tlat-lat) < 1.0 .and. abs(tlon-lon) < 1.0)exit find_thermal_sum_p
  enddo find_thermal_sum_p
  close(12)

  var_out((imonth1+1):12) = var_past_year((imonth1+1):12)

  return
end subroutine read_thermal_sum
