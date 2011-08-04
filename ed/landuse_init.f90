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
subroutine landuse_init()

  use consts_coms, only: erad, pio180
  use mem_leaf, only: first_site
  use ed_structure_defs
  use disturbance_coms, only: lutime, num_lu_trans
  use misc_coms, only: iyear1
  use ed_options, only: ianth_disturb

  implicit none

  type(site), pointer :: cs
  real :: file_lat
  real :: file_lon
  character(len=256) :: fname
  real :: lu_area
  logical :: exans
  integer :: iyear
  type(lutime), pointer :: newlutime

  cs => first_site
  do while(associated(cs))

     nullify(cs%first_lutime)
     nullify(cs%last_lutime)

     ! Generate the landuse file name
     call landuse_file_name(cs%lat, cs%lon, file_lat, file_lon, fname)

     ! Use file_lat to compute the physical area sampled by the file
     lu_area = (erad * pio180)**2 * abs(cos(pio180 * file_lat))

     inquire(file=trim(fname),exist=exans)

     if(exans .and. ianth_disturb==1)then

        ! Land use file exists
        open(12,file=trim(fname),form='formatted',status='old')
        ! Skip header
        read(12,*)
        ! Each GLU file has 300 years, 1700--1999.
        cs%num_landuse_years = 300
        ! Loop over years
        do iyear = 1,cs%num_landuse_years
           ! Allocate the memory
           nullify(newlutime)
           allocate(newlutime)
           ! Read the file
           read(12,*)newlutime%landuse_year, newlutime%landuse(1:19)
           ! Set the pointers
           nullify(newlutime%next_lutime)
           if(associated(cs%first_lutime))then
              cs%last_lutime%next_lutime => newlutime
           else
              cs%first_lutime => newlutime
           endif
           cs%last_lutime => newlutime
           ! Normalize by the area
           newlutime%landuse(12) = newlutime%landuse(12) / lu_area
           newlutime%landuse(14) = newlutime%landuse(14) / lu_area
           newlutime%landuse(16) = newlutime%landuse(16) / lu_area
           newlutime%landuse(18) = newlutime%landuse(18) / lu_area

        enddo
        close(12)

     else
        ! no GLU data for this site.  probably water, or disturbance is off.
        print*,'file: ',trim(fname),' not found.  assigning 0s for landuse'
        cs%num_landuse_years = 1
        nullify(newlutime)
        allocate(newlutime)
        newlutime%landuse_year = iyear1
        newlutime%landuse(1:num_lu_trans) = 0.0
        nullify(newlutime%next_lutime)
        cs%first_lutime => newlutime
        cs%last_lutime => newlutime
     endif

     call read_plantation_fractions(cs, file_lat, file_lon)

     cs => cs%next_site
  enddo

  return
end subroutine landuse_init

!===================================================================

subroutine landuse_file_name(lat, lon, file_lat, file_lon, fname)
  
  use ed_options, only: ed_inputs_dir

  implicit none

  real, intent(in) :: lat
  real, intent(in) :: lon
  real, intent(out) :: file_lat
  character(len=256), intent(out) :: fname
  character(len=5) :: file_lat_string
  real, intent(out) :: file_lon
  character(len=6) :: file_lon_string
  
  if(lat > 0.0)then
     file_lat = 0.5 + 1.0 * int(lat)
     if(file_lat < 10.0)then
        write(file_lat_string,'(f4.1)')file_lat
     else
        write(file_lat_string,'(f5.1)')file_lat
     endif
  else
     file_lat = - (0.5 + 1.0 * int(-lat))
     if(file_lat > -10.0)then
        write(file_lat_string,'(f4.1)')file_lat
     else
        write(file_lat_string,'(f5.1)')file_lat
     endif
  endif

  if(lon > 0.0)then
     file_lon = 0.5 + 1.0 * int(lon)
     if(file_lon < 10.0)then
        write(file_lon_string,'(f4.1)')file_lon
     elseif(lon < 100.0)then
        write(file_lon_string,'(f5.1)')file_lon
     else
        write(file_lon_string,'(f6.1)')file_lon
     endif
  else
     file_lon = - (0.5 + 1.0 * int(-lon))
     if(lon > -10.0)then
        write(file_lon_string,'(f4.1)')file_lon
     elseif(lon > -100.0)then
        write(file_lon_string,'(f5.1)')file_lon
     else
        write(file_lon_string,'(f6.1)')file_lon
     endif
  endif
  
  fname = trim(ed_inputs_dir)//'glu/lat'//trim(file_lat_string)
  fname = trim(fname)//'lon'//trim(file_lon_string)//'.lu'

  return
end subroutine landuse_file_name

!===========================================================================

subroutine read_plantation_fractions(cs, file_lat, file_lon)

  use ed_structure_defs
  use ed_options, only: ed_inputs_dir

  implicit none

  type(site)       :: cs
  real, intent(in) :: file_lat
  real, intent(in) :: file_lon
  character(len=256) :: fname
  logical :: exans
  integer :: irec
  real :: lat
  real :: lon
  real :: fracplant

  ! Set plantation fraction
  fname = trim(ed_inputs_dir)//'fraction.plantation'
  inquire(file=trim(fname),exist=exans)
  if(.not.exans)then
     print*,'There is no plantation file.  Exiting.'
     print*,trim(fname)
     stop
  endif

  open(12, file=trim(fname), form='formatted', status='old')
  read_plantation: do
     read(12,*,iostat=irec)lat, lon, fracplant
     if(irec /= 0)exit read_plantation
     if(lat == file_lat .and. lon == file_lon)then
        if(fracplant > 0.125)cs%plantation = 1
        exit read_plantation
     endif
  enddo read_plantation
  close(12)

  return
end subroutine read_plantation_fractions
