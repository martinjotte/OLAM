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
subroutine init_offline_met()

  use mem_leaf, only: first_site
  use ed_options, only: ed_offline_db
  use offline_coms, only: nformats, ed_ol_names, ed_ol_nlon, ed_ol_nlat,  &
       ed_ol_dx, ed_ol_dy, ed_ol_xmin, ed_ol_ymin, ed_ol_nv, ed_ol_vars,  &
       ed_ol_frq, ed_ol_interp, max_ol_vars
  use ed_structure_defs

  implicit none

  integer :: iformat
  integer :: iv
  type(site), pointer :: cs
  integer :: mem_size
  logical :: l1

  ! Open the config file
  inquire(file=trim(ed_offline_db),exist=l1)
  if(.not.l1)then
     print*,'File:'
     print*,trim(ed_offline_db)
     print*,'does not exist.  Specify ED_OFFLINE_DB properly in OLAMIN'
     print*
     print*,'Message brought to you by subroutine init_offline_met in ed_offline_met.f90.'
     stop
  endif

  open(12,file=trim(ed_offline_db),form='formatted',status='old')

  ! Read the number of different file formats
  read(12,*)nformats

  ! Allocate the header information for each format
  allocate(ed_ol_names(nformats))
  allocate(ed_ol_nlon(nformats))
  allocate(ed_ol_nlat(nformats))
  allocate(ed_ol_dx(nformats))
  allocate(ed_ol_dy(nformats))
  allocate(ed_ol_xmin(nformats))
  allocate(ed_ol_ymin(nformats))
  allocate(ed_ol_nv(nformats))
  allocate(ed_ol_vars(nformats, max_ol_vars))
  allocate(ed_ol_frq(nformats, max_ol_vars))
  allocate(ed_ol_interp(nformats, max_ol_vars))

  ! Read the information for each format
  do iformat = 1,nformats
     read(12,'(a)')ed_ol_names(iformat)
     read(12,*)ed_ol_nlon(iformat), ed_ol_nlat(iformat),   &
          ed_ol_dx(iformat), ed_ol_dy(iformat), ed_ol_xmin(iformat),   &
          ed_ol_ymin(iformat)
     read(12,*)ed_ol_nv(iformat)
     read(12,*)ed_ol_vars(iformat,1:ed_ol_nv(iformat))
     read(12,*)ed_ol_frq(iformat,1:ed_ol_nv(iformat))
     read(12,*)ed_ol_interp(iformat,1:ed_ol_nv(iformat))

     !  Loop over sites
     cs => first_site
     do while(associated(cs))

        ! Make sure site falls within file domain
        if(cs%lon < (ed_ol_xmin(iformat) - 0.5 * ed_ol_dx(iformat)) .or.  &
             cs%lat < (ed_ol_ymin(iformat) - 0.5 * ed_ol_dy(iformat)) .or.  &
             cs%lon > (ed_ol_xmin(iformat) + (ed_ol_nlon(iformat)-1) *  &
             ed_ol_dx(iformat) + 0.5 * ed_ol_dx(iformat)) .or.  &
             cs%lat > (ed_ol_ymin(iformat) + (ed_ol_nlat(iformat)-1) *  &
             ed_ol_dy(iformat) + 0.5 * ed_ol_dy(iformat)) )then
           print*
           print*,'========================================================'
           print*,'Site is not within the domain of the meteorological drivers'
           print*,'========================================================'
           print*,iformat,trim(ed_ol_names(iformat))
           print*,cs%lon, ed_ol_xmin(iformat) - 0.5 * ed_ol_dx(iformat)
           print*,cs%lat, ed_ol_ymin(iformat) - 0.5 * ed_ol_dy(iformat)
           print*,cs%lon, ed_ol_xmin(iformat) + (ed_ol_nlon(iformat)-1) *  &
                ed_ol_dx(iformat) + 0.5 * ed_ol_dx(iformat)
           print*,cs%lat, ed_ol_ymin(iformat) + (ed_ol_nlat(iformat)-1) *  &
                ed_ol_dy(iformat) + 0.5 * ed_ol_dy(iformat)
           stop
        endif

        ! Allocate memory
        do iv = 1,ed_ol_nv(iformat)

           ! Calculate number of points in a month
           mem_size = nint(86400.0 / ed_ol_frq(iformat,iv)) * 31
           
           ! If this is an interpolated variable, it is best
           ! to have two months read in at a time.
           if(ed_ol_interp(iformat,iv) == 1)mem_size = 2 * mem_size
           
           if(trim(ed_ol_vars(iformat,iv)) == 'nbdsf')then
              allocate(cs%metinput%nbdsf(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'nddsf')then
              allocate(cs%metinput%nddsf(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vbdsf')then
              allocate(cs%metinput%vbdsf(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vddsf')then
              allocate(cs%metinput%vddsf(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'prate')then
              allocate(cs%metinput%prate(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'dlwrf')then
              allocate(cs%metinput%dlwrf(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'pres')then
              allocate(cs%metinput%pres(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'hgt')then
              allocate(cs%metinput%hgt(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'ugrd')then
              allocate(cs%metinput%ugrd(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vgrd')then
              allocate(cs%metinput%vgrd(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'sh')then
              allocate(cs%metinput%sh(mem_size))
           elseif(trim(ed_ol_vars(iformat,iv)) == 'tmp')then
              allocate(cs%metinput%tmp(mem_size))
           endif

        enddo

        cs => cs%next_site
     enddo
     
  enddo

  ! Close the header file and return
  close(12)

  return
end subroutine init_offline_met

!======================================================================

subroutine read_offline_met_init()

  use offline_coms, only: nformats, ed_ol_names, ed_ol_nv, ed_ol_interp,  &
       ed_ol_frq
  use misc_coms,  only: current_time
  use ed_options, only: metcyc1, metcyc2
  use hdf5_utils, only: shdf5_open, shdf5_close
  use max_dims,   only: pathlen

  implicit none

  integer :: year_use
  integer :: ncyc
  integer :: iformat
  character(pathlen) :: infile
  character(3), dimension(12), parameter :: mname = (/'JAN', 'FEB',   &
       'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/) 
  logical :: exans
  integer :: iv
  integer :: offset
  integer :: m2
  integer :: y2
  integer :: year_use_2

  ! If we need to recycle over years, find the appropriate year to apply.
  year_use = current_time%year
  ncyc = metcyc2 - metcyc1 + 1

  ! If we are after the last year...
  do while(year_use > metcyc2)
     year_use = year_use - ncyc
  enddo

  ! If we are before the first year...
  do while(year_use < metcyc1)
     year_use = year_use + ncyc
  enddo

  ! Loop over the different file formats
  do iformat = 1, nformats

     ! Open the file
     write(infile,'(a,i4.4,a,a)')trim(ed_ol_names(iformat)), year_use,   &
          mname(current_time%month),'.h5'
     inquire(file=trim(infile),exist=exans)
     if(exans)then
        call shdf5_open(trim(infile),'R')
     else
        print*,'Cannot open met driver input file',trim(infile)
        stop
     endif
     
     ! Loop over variables
     do iv = 1, ed_ol_nv(iformat)
        
        offset = 0
        call read_ol_file(iformat, iv, year_use, mname(current_time%month),  &
             current_time%year, offset)

     enddo

     ! Close the HDF5 file.
     call shdf5_close()

     ! For all interpolated variables, we also need the next time.

     ! Find next month and year
     m2 = current_time%month + 1
     y2 = current_time%year
     year_use_2 = year_use
     
     ! If this takes us into the next year, increment year and 
     ! reset month to January.
     if(m2 == 13)then
        m2 = 1
        y2 = current_time%year + 1
        year_use_2 = y2
        
        ! If we are now after the last year...
        do while(year_use_2 > metcyc2)
           year_use_2 = year_use_2 - ncyc
        enddo
        
        ! If we are now before the first year...
        do while(year_use_2 < metcyc1)
           year_use_2 = year_use_2 + ncyc
        enddo
     endif

     ! Now, open the file once.
     write(infile,'(a,i4.4,a,a)')trim(ed_ol_names(iformat)), year_use_2,   &
          mname(m2),'.h5'
     inquire(file=trim(infile),exist=exans)
     if(exans)then
        call shdf5_open(trim(infile),'R')
     else
        print*,'Cannot open met driver input file',trim(infile)
        stop
     endif
     
     ! Loop over variables
     do iv = 1, ed_ol_nv(iformat)
        
        if(ed_ol_interp(iformat,iv) == 1)then
           
           ! Read the file.
           offset = nint(86400.0 / ed_ol_frq(iformat,iv)) * 31
           call read_ol_file(iformat, iv, year_use_2, mname(m2),  &
                y2, offset)

        endif
     enddo

     ! Close the HDF5 file.
     call shdf5_close()

  enddo

  return
end subroutine read_offline_met_init

!======================================================================

subroutine read_offline_met()

  use offline_coms, only: nformats, ed_ol_names, ed_ol_nv, ed_ol_interp,  &
       ed_ol_frq, ed_ol_vars
  use misc_coms,  only: current_time
  use ed_options, only: metcyc1, metcyc2
  use hdf5_utils, only: shdf5_open, shdf5_close
  use max_dims,   only: pathlen

  implicit none

  integer :: year_use
  integer :: ncyc
  integer :: iformat
  character(pathlen) :: infile
  character(len=3), dimension(12), parameter :: mname = (/'JAN', 'FEB',   &
       'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/) 
  logical :: exans
  integer :: iv
  integer :: offset
  integer :: m2
  integer :: y2
  integer :: year_use_2

  ! If we need to recycle over years, find the appropriate year to apply.
  year_use = current_time%year
  ncyc = metcyc2 - metcyc1 + 1

  ! If we are after the last year...
  do while(year_use > metcyc2)
     year_use = year_use - ncyc
  enddo

  ! If we are before the first year...
  do while(year_use < metcyc1)
     year_use = year_use + ncyc
  enddo

  ! Loop over the different file formats
  do iformat = 1, nformats

     ! Open the file
     write(infile,'(a,i4.4,a,a)')trim(ed_ol_names(iformat)), year_use,   &
          mname(current_time%month),'.h5'
     inquire(file=trim(infile),exist=exans)
     if(exans)then
        call shdf5_open(trim(infile),'R')
     else
        print*,'Cannot open met driver input file',trim(infile)
        stop
     endif
     
     ! Loop over variables
     do iv = 1, ed_ol_nv(iformat)
        
        ! See if it is an interpolation variable.
        if(ed_ol_interp(iformat,iv) == 0)then
           
           ! If not, things are simple.  Just read in the month.
           offset = 0
           call read_ol_file(iformat, iv, year_use,   &
                mname(current_time%month), current_time%year, offset)

        else

           ! Here, just transfer future to current month.
           call transfer_ol_month(trim(ed_ol_vars(iformat,iv)),   &
                ed_ol_frq(iformat,iv))

        endif

     enddo

     ! Close the HDF5 file.
     call shdf5_close()

     ! For all interpolated variables, get the future month.

     ! Find next month and year
     m2 = current_time%month + 1
     y2 = current_time%year
     year_use_2 = year_use
     
     ! If this takes us into the next year, increment year and 
     ! reset month to January.
     if(m2 == 13)then
        m2 = 1
        y2 = current_time%year + 1
        year_use_2 = y2
        
        ! If we are now after the last year...
        do while(year_use_2 > metcyc2)
           year_use_2 = year_use_2 - ncyc
        enddo
        
        ! If we are now before the first year...
        do while(year_use_2 < metcyc1)
           year_use_2 = year_use_2 + ncyc
        enddo
     endif

     ! Now, open the file once.
     write(infile,'(a,i4.4,a,a)')trim(ed_ol_names(iformat)), year_use_2,   &
          mname(m2),'.h5'
     inquire(file=trim(infile),exist=exans)
     if(exans)then
        call shdf5_open(trim(infile),'R')
     else
        print*,'Cannot open met driver input file',trim(infile)
        stop
     endif
     
     ! Loop over variables
     do iv = 1, ed_ol_nv(iformat)
        
        if(ed_ol_interp(iformat,iv) == 1)then
           
           ! Read the file.
           offset = nint(86400.0 / ed_ol_frq(iformat,iv)) * 31
           call read_ol_file(iformat, iv, year_use_2, mname(m2),  &
                y2, offset)

        endif
     enddo

     ! Close the HDF5 file.
     call shdf5_close()

  enddo

  return
end subroutine read_offline_met

!===============================================================

subroutine update_offline_met()
  
  use offline_coms, only: ed_ol_frq
  use misc_coms, only: current_time
  use ed_options, only: frq_met_ol
  use ed_structure_defs
  use mem_leaf, only: land, first_site
  use consts_coms, only: rdry, cice, cliq, alli
  use misc_coms, only: dtlm

  implicit none

  integer :: points_per_day
  integer :: nday
  integer :: np
  integer :: ndays_elapsed
  integer :: nseconds_elapsed
  integer :: mlo
  integer :: mhi
  real :: t1
  real :: t2
  type(site), pointer :: cs
  integer :: iformat
  integer :: iv

  ! Initialize vels
  cs => first_site
  do while(associated(cs))
     land%vels(cs%iland) = 0.0
     cs => cs%next_site
  enddo

  ! Loop over the different file formats
  do iformat = 1, nformats

     ! Loop over variables
     do iv = 1, ed_ol_nv(iformat)

        if(ed_ol_interp(iformat,iv) == 0)then

           ! If this is not an interpolation variable, just find the time
           ! point.
           ndays_elapsed = current_time%date - 1
           nseconds_elapsed = nint(current_time%time) + ndays_elapsed * 86400
           mlo = int(float(nseconds_elapsed) / ed_ol_frq(iformat,iv)) + 1

           ! Find which variable it is, and then fill the sites.
           if(trim(ed_ol_vars(iformat,iv)) == 'nbdsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%nir_beam = cs%metinput%nbdsf(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'nddsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%nir_diffuse = cs%metinput%nddsf(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vbdsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%par_beam = cs%metinput%vbdsf(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vddsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%par_diffuse = cs%metinput%vddsf(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'prate')then
              cs => first_site
              do while(associated(cs))
                 land%pcpg(cs%iland) = dtlm(1) * cs%metinput%prate(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'dlwrf')then
              cs => first_site
              do while(associated(cs))
                 land%rlong(cs%iland) = cs%metinput%dlwrf(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'pres')then
              cs => first_site
              do while(associated(cs))
                 land%prss(cs%iland) = cs%metinput%pres(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'hgt')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%geoht = cs%metinput%hgt(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'ugrd')then
              cs => first_site
              do while(associated(cs))
                 land%vels(cs%iland) = land%vels(cs%iland) +   &
                      cs%metinput%ugrd(mlo)**2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vgrd')then
              cs => first_site
              do while(associated(cs))
                 land%vels(cs%iland) = land%vels(cs%iland) +   &
                      cs%metinput%vgrd(mlo)**2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'sh')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%atm_shv = cs%metinput%sh(mlo)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'tmp')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%atm_tmp = cs%metinput%tmp(mlo)
                 cs => cs%next_site
              enddo
           endif

        else

           ! In this case, we need to interpolate.
           
           ! First, get the number of points per day and per month.
           points_per_day = nint(86400.0/ed_ol_frq(iformat,iv))
           if(current_time%month /= 2 .and. current_time%month /= 4 .and.   &
                current_time%month == 6 .and. current_time%month == 9  .and.  &
                current_time%month == 11)then
              nday = 31
           elseif(current_time%month /= 2)then
              nday = 30
           elseif(mod(current_time%year,4) /= 0)then
              nday = 28
           else
              nday = 29
           endif
           np = nday * points_per_day

           ! Get indices
           mlo = int(float(nseconds_elapsed)/ed_ol_frq(iformat,iv)) + 1
           mhi = mlo + 1
           if(mhi > np)then
              mhi = 1 + nint(86400.0 / ed_ol_frq(iformat,iv)) * 31
           endif

           ! Get interpolation factors
           t1 = mod(float(nseconds_elapsed), ed_ol_frq(iformat,iv)) /   &
                ed_ol_frq(iformat,iv)
           t2 = 1.0 - t1

           ! Find the variable and fill the sites.
           if(trim(ed_ol_vars(iformat,iv)) == 'prate')then
              cs => first_site
              do while(associated(cs))
                 land%pcpg(cs%iland) = dtlm(1) * (cs%metinput%prate(mhi) *   &
                      t1 + cs%metinput%prate(mlo) * t2)
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'dlwrf')then
              cs => first_site
              do while(associated(cs))
                 land%rlong(cs%iland) = cs%metinput%dlwrf(mhi) * t1 +  &
                      cs%metinput%dlwrf(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'pres')then
              cs => first_site
              do while(associated(cs))
                 land%prss(cs%iland) = cs%metinput%pres(mhi) * t1 +   &
                      cs%metinput%pres(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'hgt')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%geoht = cs%metinput%hgt(mhi) * t1 +  &
                      cs%metinput%hgt(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'ugrd')then
              cs => first_site
              do while(associated(cs))
                 land%vels(cs%iland) = land%vels(cs%iland) +   &
                      (cs%metinput%ugrd(mhi) * t1 +   &
                      cs%metinput%ugrd(mlo) * t2)**2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vgrd')then
              cs => first_site
              do while(associated(cs))
                 land%vels(cs%iland) = land%vels(cs%iland) +   &
                      (cs%metinput%vgrd(mhi) * t1 +  &
                      cs%metinput%vgrd(mlo) * t2)**2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'sh')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%atm_shv = cs%metinput%sh(mhi) * t1 +   &
                      cs%metinput%sh(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'tmp')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%atm_tmp = cs%metinput%tmp(mhi) * t1 +  &
                      cs%metinput%tmp(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'nbdsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%nir_beam = cs%metinput%nbdsf(mhi) * t1 +  &
                      cs%metinput%nbdsf(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'nddsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%nir_diffuse = cs%metinput%nddsf(mhi) * t1 +  &
                      cs%metinput%nddsf(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vbdsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%par_beam = cs%metinput%vbdsf(mhi) * t1 +  &
                      cs%metinput%vbdsf(mlo) * t2
                 cs => cs%next_site
              enddo
           elseif(trim(ed_ol_vars(iformat,iv)) == 'vddsf')then
              cs => first_site
              do while(associated(cs))
                 cs%metinput%par_diffuse = cs%metinput%vddsf(mhi) * t1 +  &
                      cs%metinput%vddsf(mlo) * t2
                 cs => cs%next_site
              enddo
           endif
           
        endif
        
     enddo
  enddo
  
  ! Change from velocity squared to velocity, get the rhos, and compute
  ! qpcpg and dpcpg.
  cs => first_site
  do while(associated(cs))

     ! vels
     land%vels(cs%iland) = sqrt(max(0.0,land%vels(cs%iland)))

     ! rho
     land%rhos(cs%iland) = land%prss(cs%iland) / (rdry *   &
          cs%metinput%atm_tmp * (1.0 + 0.61 * cs%metinput%atm_shv))

     ! qpcpg, dpcpg
     if(cs%metinput%atm_tmp > 273.15)then
        land%qpcpg(cs%iland) = (cliq * (cs%metinput%atm_tmp - 273.15) +   &
             alli) * land%pcpg(cs%iland)
     else
        land%qpcpg(cs%iland) = cice * (cs%metinput%atm_tmp - 273.15) *  &
             land%pcpg(cs%iland)
     endif
     land%dpcpg(cs%iland) = max(0.0, land%pcpg(cs%iland) * 0.001)

     cs => cs%next_site
  enddo

  return
end subroutine update_offline_met

!========================================================================

subroutine read_ol_file(iformat, iv, year_use, mname, year, offset)

  use offline_coms, only: ed_ol_frq, ed_ol_nlon, ed_ol_nlat, ed_ol_vars,   &
       ed_ol_xmin, ed_ol_dx, ed_ol_ymin, ed_ol_dy
  use hdf5_utils, only: shdf5_irec
  use ed_structure_defs
  use mem_leaf, only: first_site

  implicit none

  integer, intent(in) :: iformat
  integer, intent(in) :: iv
  integer, intent(in) :: year_use
  character(len=3), intent(in) :: mname
  integer, intent(in) :: year
  integer, intent(in) :: offset

  integer :: points_per_day
  integer :: nday
  integer :: np
  real, allocatable, dimension(:,:,:) :: metvar
  integer :: ndims
  integer, dimension(3) :: idims
  type(site), pointer :: cs
  integer :: ilon
  integer :: ilat
  
  ! Determine how many times to read
  points_per_day = nint(86400.0/ed_ol_frq(iformat,iv))
  if(mname /= 'FEB' .and. mname /= 'APR' .and. mname /= 'JUN' .and.   &
       mname /= 'SEP' .and. mname /= 'NOV')then
     nday = 31
  elseif(mname /= 'FEB')then
     nday = 30
  elseif(mod(year_use,4) /= 0)then
     nday = 28
  else
     nday = 29
  endif
  np = nday * points_per_day
  
  ! Allocate the buffer space
  allocate(metvar(np,ed_ol_nlon(iformat),ed_ol_nlat(iformat)))

  ! Get the variable
  ndims = 3
  idims(1) = np 
  idims(2) = ed_ol_nlon(iformat)
  idims(3) = ed_ol_nlat(iformat)
  call shdf5_irec(ndims, idims, trim(ed_ol_vars(iformat,iv)),  &
       rvara = metvar)
  
  ! Assign the data to sites
  cs => first_site
  do while(associated(cs))
     
     ! Get the indices.  Remember, latitude is flipped.
     ilon = min( max(1, 1 + nint( (cs%lon - ed_ol_xmin(iformat)) /   &
          ed_ol_dx(iformat) ) ), ed_ol_nlon(iformat))
     ilat = ed_ol_nlat(iformat) -  &
          min( max(1, 1 + nint( (cs%lat - ed_ol_ymin(iformat)) /   &
          ed_ol_dy(iformat) ) ), ed_ol_nlat(iformat)) + 1

     ! Deal with possible leap-year mismatch
     if(mname == 'FEB' .and. mod(year,4) == 0 .and.   &
          mod(year_use,4) /= 0)then
        metvar((np+1):(np+points_per_day), ilon, ilat) =   &
             metvar((np-points_per_day+1):np, ilon, ilat)
        np = np + points_per_day
     endif
     
     ! Get the time series.
     if(trim(ed_ol_vars(iformat,iv)) == 'nbdsf')then
        cs%metinput%nbdsf((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'nddsf')then
        cs%metinput%nddsf((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'vbdsf')then
        cs%metinput%vbdsf((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'vddsf')then
        cs%metinput%vddsf((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'prate')then
        cs%metinput%prate((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'dlwrf')then
        cs%metinput%dlwrf((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'pres')then
        cs%metinput%pres((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'hgt')then
        cs%metinput%hgt((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'ugrd')then
        cs%metinput%ugrd((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'vgrd')then
        cs%metinput%vgrd((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'sh')then
        cs%metinput%sh((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     elseif(trim(ed_ol_vars(iformat,iv)) == 'tmp')then
        cs%metinput%tmp((offset+1):(offset+np)) = metvar(1:np,ilon,ilat)
     endif
     
     cs => cs%next_site
  enddo
  
  ! Deallocate the buffer
  deallocate(metvar)

  return
end subroutine read_ol_file

!======================================================================================

subroutine transfer_ol_month(vname, frq)

  use ed_structure_defs
  use mem_leaf, only: first_site

  implicit none

  character(*) :: vname
  type(site), pointer :: cs
  real, intent(in) :: frq
  integer :: mem_size

  mem_size = nint(86400.0 / frq) * 31

  cs => first_site
  do while(associated(cs))
  
     ! Get the time series.
     if(trim(vname) == 'nbdsf')then
        cs%metinput%nbdsf(1:mem_size) = cs%metinput%nbdsf((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'nddsf')then
        cs%metinput%nddsf(1:mem_size) = cs%metinput%nddsf((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'vbdsf')then
        cs%metinput%vbdsf(1:mem_size) = cs%metinput%vbdsf((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'vddsf')then
        cs%metinput%vddsf(1:mem_size) = cs%metinput%vddsf((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'prate')then
        cs%metinput%prate(1:mem_size) = cs%metinput%prate((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'dlwrf')then
        cs%metinput%dlwrf(1:mem_size) = cs%metinput%dlwrf((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'pres')then
        cs%metinput%pres(1:mem_size) = cs%metinput%pres((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'hgt')then
        cs%metinput%hgt(1:mem_size) = cs%metinput%hgt((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'ugrd')then
        cs%metinput%ugrd(1:mem_size) = cs%metinput%ugrd((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'vgrd')then
        cs%metinput%vgrd(1:mem_size) = cs%metinput%vgrd((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'sh')then
        cs%metinput%sh(1:mem_size) = cs%metinput%sh((mem_size+1):(2*mem_size))
     elseif(trim(vname) == 'tmp')then
        cs%metinput%tmp(1:mem_size) = cs%metinput%tmp((mem_size+1):(2*mem_size))
     endif
     
     cs => cs%next_site
  enddo
  
  

  return
end subroutine transfer_ol_month
