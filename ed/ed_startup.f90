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
subroutine initialize_ed_sites()
  implicit none
  
  ! This function allocates the memory, sets the latitude and longitude, and finds the corresponding iland index for all grid cells at which ED will be run.
  call spawn_sites()

  call read_ed_database()

  return
end subroutine initialize_ed_sites

!====================================================================

subroutine spawn_sites()
  use mem_leaf,    only: land, first_site, itab_wl
  use mem_ed,      only: n_soi, soi_lat, soi_lon, n_ed_region,  &
       ed_reg_latmin, ed_reg_latmax, ed_reg_lonmin, ed_reg_lonmax
  use ed_structure_defs
  use consts_coms, only: erad, pio180, piu180
  use leaf_coms,   only: mwl

  implicit none

  integer :: isoi
  type(site), pointer :: new_site
  real :: xe_soi
  real :: ye_soi
  real :: ze_soi
  real :: shortest_distance2
  integer :: iwl
  real :: xe_mean
  real :: ye_mean
  real :: ze_mean
  integer :: ipt
  real :: distance2
  integer :: ireg
  real :: lat_mean
  real :: lon_mean
  type(site), pointer :: current_site
  type(site), pointer :: cs
  type(site), pointer :: pbelow
  type(site), pointer :: pabove
  integer :: min_iland_diff
  integer :: iland_diff

  ! First, allocate memory for the SOIs (sites of interest)
  do isoi = 1,n_soi
     nullify(new_site)
     allocate(new_site)
     ! Loop over OLAM land surface locations and find nearest to this SOI.
     xe_soi = erad * cos(pio180*soi_lat(isoi)) * cos(pio180*soi_lon(isoi))
     ye_soi = erad * cos(pio180*soi_lat(isoi)) * sin(pio180*soi_lon(isoi))
     ze_soi = erad * sin(pio180*soi_lat(isoi))
     shortest_distance2 = 1.0e30

     do iwl = 2,mwl

        xe_mean = land%xew(iwl)
        ye_mean = land%yew(iwl)
        ze_mean = land%zew(iwl)

        distance2 = (xe_soi - xe_mean)**2 + (ye_soi - ye_mean)**2 + &
             (ze_soi - ze_mean)**2

        if(distance2.lt.shortest_distance2)then
           shortest_distance2 = distance2
           new_site%iland = iwl
           new_site%lat = land%glatw(iwl)
           new_site%lon = land%glonw(iwl)
        endif
     enddo

     ! add site to the linked list
     if(associated(first_site))then
        ! find the site that is just below
        min_iland_diff = 1e8
        nullify(pbelow)
        cs => first_site
        do while(associated(cs))
           iland_diff = new_site%iland - cs%iland
           if(iland_diff > 0 .and. iland_diff < min_iland_diff)then
              min_iland_diff = iland_diff
              pbelow => cs
           endif
           cs => cs%next_site
        enddo
        if(associated(pbelow))then
           if(associated(pbelow%next_site))then
              pabove => pbelow%next_site
              new_site%next_site => pabove
              pbelow%next_site => new_site
           else
              nullify(new_site%next_site)
              pbelow%next_site => new_site
           endif
        else
           pabove => first_site
           new_site%next_site => pabove
           first_site => new_site
        endif
     else
        nullify(new_site%next_site)
        first_site => new_site
     endif

  enddo

  ! Now, allocate memory for the rectangular areas
  do ireg = 1,n_ed_region
     do iwl = 2,mwl

        lat_mean = land%glatw(iwl)
        lon_mean = land%glonw(iwl)

        if(lat_mean >= ed_reg_latmin(ireg)     .and. &
             lat_mean <= ed_reg_latmax(ireg)   .and. &
             lon_mean >= ed_reg_lonmin(ireg)   .and. &
             lon_mean <= ed_reg_lonmax(ireg))then
           nullify(new_site)
           allocate(new_site)
           new_site%iland = iwl
           new_site%lat = lat_mean
           new_site%lon = lon_mean

           ! add site to the linked list
           if(associated(first_site))then
              ! find the site that is just below
              min_iland_diff = 1e8
              nullify(pbelow)
              cs => first_site
              do while(associated(cs))
                 iland_diff = new_site%iland - cs%iland
                 if(iland_diff > 0 .and. iland_diff < min_iland_diff)then
                    min_iland_diff = iland_diff
                    pbelow => cs
                 endif
                 cs => cs%next_site
              enddo
              if(associated(pbelow))then
                 if(associated(pbelow%next_site))then
                    pabove => pbelow%next_site
                    new_site%next_site => pabove
                    pbelow%next_site => new_site
                 else
                    nullify(new_site%next_site)
                    pbelow%next_site => new_site
                 endif
              else
                 pabove => first_site
                 new_site%next_site => pabove
                 first_site => new_site
              endif
           else
              nullify(new_site%next_site)
              first_site => new_site
           endif
        endif
     enddo
  enddo

  current_site => first_site
  do while(associated(current_site))
     land%ed_flag(current_site%iland) = 1
     current_site => current_site%next_site
  enddo

  return
end subroutine spawn_sites

