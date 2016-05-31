!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

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

program grib_to_gdf

  use grib_get_mod
  use hdf5_utils

  implicit none

  integer, parameter   :: max_hours=24
  character(len= 10)   :: c_date, c_hr(max_hours), out_ext
  character(len=128)   :: filein,fileout,out_prefix,dstring,namein
  character(len= 20)   :: cargv
  integer, allocatable :: itlevs(:)
  real,    allocatable :: a3(:,:,:,:), a2(:,:,:)

  integer :: ia,iargc,istr,iend,inc,nsea
  integer :: i,j,idprep=0,ndate,nhr,ii,nh,nb,n
  integer :: nlev,inproj,ivertcoord,miss,ircode,num_hours
  integer :: nsoilw, nsoilt, ngnd
  
  integer :: iyyyy,imm,idd,ihh,nuyyyy,numm,nudd,nuhh
  real :: tinc,xnelat,xnelon,aver
  logical :: exists, dousage

  integer, parameter :: maxlev=200, maxvar=100, maxgnd=9
  integer :: iplevs(maxlev), iswlevs(maxgnd), istlevs(maxgnd), islevs(maxgnd)

  integer :: nvar3d, nvarsfc, nvargnd, gdf_file_ver, ndims, idims(3)
  character(len=10) :: var3d(maxvar), varsfc(maxvar) !, vargnd(maxvar)
  character(len=10) :: out3d(maxvar), outsfc(maxvar) !, outgnd(maxvar)

  character(len=10) :: sstvar, landmask, icevar
  logical           :: have_sst, have_ice

  character(len=1)  :: zsoil

  real, parameter   :: t00sea  = 271.38
  real, parameter   :: icefmin = 0.001

  real, allocatable :: landf(:,:), sst(:,:), ice(:,:), snow(:,:)
  integer, allocatable :: land(:,:)

  real, allocatable :: soilw(:,:,:), soilt(:,:,:)

  character(1), parameter :: anum(maxgnd) = ['1','2','3','4','5','6','7','8','9']

  character(10) :: soilw_var
  character(10) :: soilw_vars(maxgnd)
  logical       :: have_soilw

  character(10) :: soilt_var
  character(10) :: soilt_vars(maxgnd)
  logical       :: have_soilt

  character(10) :: snow_var
  logical       :: have_snow

  namelist /dgrib_in/ nvar3d, var3d, nvarsfc, varsfc, c_date, num_hours, c_hr, &
                     filein, out3d, outsfc, wgrib1_exe, wgrib2_exe, out_prefix, &
!                    nvargnd, vargnd, outgnd, &
                     sstvar, landmask, icevar, &
                     soilw_var, soilt_var, snow_var

  dousage = .false.
  out_ext = '.h5'
  namein  = 'DGRIB_IN'

! Initialize the namelist variables

  wgrib1_exe = './wgrib'
  wgrib2_exe = './wgrib2'
  num_hours  = 1
  c_date     = '99999999'
  c_hr(:)    = '99999999'
  out_prefix = ''
  var3d(:)   = ''
  varsfc(:)  = ''
!  vargnd(:)  = ''
  out3d(:)   = ''
  outsfc(:)  = ''
!  outgnd(:)  = ''
  nvar3d     = 0
  nvarsfc    = 0
!  nvargnd    = 0
  filein     = ''

  sstvar     = ''
  icevar     = ''
  landmask   = ''
  soilw_var  = ''
  soilt_var  = ''
  snow_var   = ''

  ia = command_argument_count()
! write(*,*) 'num args: ', ia

  do i = 1, ia
     call get_command_argument(i,cargv)
     if (cargv(1:2) == '-f') then
        namein = ""
        call get_command_argument(i+1,namein)
        exit
     endif
  enddo

  inquire( file=namein, exist=exists )
  if (.not. exists) then
     write(*,*) "Error: namelist file " // trim(namein) // " does not exist."
     dousage = .true.
  else
     open (10, file=namein, status='old')
     read (10, dgrib_in)
     close(10)
  endif

  i = 1
  do while( i <= ia )
     
     call get_command_argument(i,cargv)
     
     if (cargv(1:2) == '-t') then
        c_date = ""
        call get_command_argument(i+1,c_date)
        if (len_trim(c_date)<1) dousage = .true.
        i = i + 2
     elseif (cargv(1:2) == '-h') then
        c_hr(1) = ""
        call get_command_argument(i+1,c_hr(1))
        if (len_trim(c_hr(1))<1) dousage = .true.
        num_hours = 1   ! Only 1 time allowed if argument
        i = i + 2
     elseif (cargv(1:2) == '-g') then
        filein = ""
        call get_command_argument(i+1,filein)
        if (len_trim(filein)<1) dousage = .true.
        i = i + 2
     elseif (cargv(1:3) == '-e1') then
        wgrib1_exe = ""
        call get_command_argument(i+1,wgrib1_exe)
        if (len_trim(wgrib1_exe)<1) dousage = .true.
        i = i + 2
     elseif (cargv(1:3) == '-e2') then
        wgrib2_exe = ""
        call get_command_argument(i+1,wgrib2_exe)
        if (len_trim(wgrib2_exe)<1) dousage = .true.
        i = i + 2
     elseif (cargv(1:2) == '-f') then
        ! we already processed the namelist
        i = i + 2
     else
        dousage = .true.
        exit
     endif
  enddo

  if (dousage) then
     print*
     print*,'Usage: grib2olam [-f namelist_file] [-g gribfile] [-t yyyymmddhh] [-h hh]'
     print*,'      namelist file name defaults to DGRIB_IN'
     print*,'      -t : data time to process - yyyymmddhh'
     print*,'      -h : forecast hour - hh'
     print*,'      -e1: wgrib  executable path'
     print*,'      -e2: wgrib2 executable path'
     print*
     print*,'           If c_date or c_hr set to ''99999999'', '
     print*,'           the date/hour of the first record'
     print*,'           in the grib file will be used.'
     print*
     print*,'      Command line arguments will override the namelist!'
     stop 'grib2olam usage'
  endif

  print*
  print*,'---------------------------------------------------------------'
  print*,'GRIB to GDF converter'
  print*,'---------------------------------------------------------------'
  print*,'Namelist file name: ',trim(namein)
  print*,'GRIB file name    : ',trim(filein)
  print*,'Data date/time    : ', c_date
  print*,'Forecast hour     : ', c_hr(1:num_hours)
  print*,'3D  variables: ', var3d(1:nvar3d)
  print*,'SFC variables: ', varsfc(1:nvarsfc)
  print*

  do i = 1, nvar3d
     if (len_trim(var3d(i)) < 1) then
        write(*,*) " No input name given for 3d variable ", i
        stop
     endif
     if (len_trim(out3d(i)) < 1) then
        write(*,*) " No output name given for 3d variable ", i
        stop
     endif
  enddo

  do i = 1, nvarsfc
     if (len_trim(varsfc(i)) < 1) then
        write(*,*) " No input name given for sfc variable ", i
        stop
     endif
     if (len_trim(outsfc(i)) < 1) then
        write(*,*) " No output name given for sfc variable ", i
        stop
     endif
  enddo

!!  do i = 1, nvargnd
!!     if (len_trim(vargnd(i)) < 1) then
!!        write(*,*) " No input name given for gnd variable ", i
!!        stop
!!     endif
!!     if (len_trim(outgnd(i)) < 1) then
!!        write(*,*) " No output name given for gnd variable ", i
!!        stop
!!     endif
!!  enddo

  ! Does the grib file exist?

  inquire( file=filein, exist=exists )
  if (.not. exists) then
     write(*,*)
     write(*,*) "Error: grib file " // trim(filein) // " does not exist."
     write(*,*)
     stop
  endif

  call grib_get_vers(filein)

  call grib_get_date(filein)

  call grib_get_proj(filein)

  call grib_get_fields(filein)

  print*
  print '(a,a,i4)',' Found projection: ',trim(projection), iscan
  print '(a)', '  params:  lat1     lon1     nx     ny       dx          dy'
  print '(a,2f8.2,2i7,2f12.4)','        ',alat1,alon1, nx,ny,dx,dy
  print*

  ! Special for OLAM: We need global a dataset

  if ( (projection /= "latlon" .and. projection /= "gaussian") .or. &
       (nx*dx < 360.-dx) .or. (ny*dy < 180.-dy) ) then
     write(*,*)
     write(*,*) "WARNING: OLAM requires a global lat/lon dataset!"
     write(*,*) "This dataset is not sufficient to initialize OLAM."
     write(*,*)
  endif

  ! Some bookeeping

  xnelat=0.
  xnelon=0.
  ivertcoord=1
  if (projection(1:6)=="latlon") then
     inproj=1
     if (iscan == 0 ) alat1=alat1 - dy * (ny-1)
  elseif (projection(1:8)=="gaussian") then
     inproj=2
  elseif (projection(1:7)=="Lambert") then
     inproj=3
     dx=dx*1000.
     dy=dy*1000.
     aorient=alov
  elseif (projection(1:5)=="polar") then
     inproj=4
     ! Assume the polar stereo is true polar stereo
     reflat1=90.
     reflat2=amiss
  endif

  ndate=0

  hour_loop: do nh = 1, num_hours
   
     nhr=0
   
     print*,'Processing date/time: ', c_date, c_hr(nh)
     print*
   
     read(c_hr(nh),'(i10)') nhr
     if(c_hr(nh)   == '99999999') nhr=ifhrs(1)
     if(c_date == '99999999') then
        ndate=longdate
        write(c_date,'(i10)')ndate
     endif

     print*,'dates1: ', c_hr(nh), c_date, nhr, ndate

     !     wgrib only prints out the 2 digit year in short inventory
     !     wgrib2 has 4
     if (grib_ver == 1) read(c_date(3:),'(i8)')ndate

     !     Quick check to see if there are any fields matching time wanted
     print*,'dates2:',nrec
     j=0
     do i=1,nrec
        !if (trim(fields(i)) == trim(var3d(1)) .and.   &
        !    idates(i) == ndate .and. ifhrs(i) == nhr) then
        if (idates(i) == ndate .and. ifhrs(i) == nhr) then
           j=1
           exit
        endif
     enddo
     if ( j == 0) then
        print*,'Nothing found of the proper date'
        stop 'no stuff'
     endif

     ! Determine the number of pressure levels we are going to use

     nlev=0
     do i=1,nrec
        if (levs(i)<=0) cycle  !don't count the surface fields here
        if ( trim(fields(i)) == trim(var3d(1)) .and.  &
             idates(i) == ndate .and. ifhrs(i) == nhr ) then
           nlev=nlev+1
           iplevs(nlev)=levs(i)
        endif
     enddo
   
     ! Make sure pressure levels monotonically decrease

     allocate (itlevs(nlev))
     itlevs=-1
     do i=1,nlev
        do j=1,nlev
           if(iplevs(j) > itlevs(i))then
              itlevs(i)=iplevs(j)
              ia=j
           endif
        enddo
        !  Effectively remove current highest pressure level
        iplevs(ia)=-1
     enddo
     iplevs(1:nlev)=itlevs(1:nlev)
     deallocate (itlevs)

     print*
     print '(1x,I0,A)', nlev, " Pressure levels found:"
     print*, iplevs(1:nlev)
     print*

!!     ! Determine the number of soil layers we are going to use
!!
!!     if (nvargnd > 0) then
!!
!!        ngnd = 0
!!        do i=1, nrec
!!           if (levs(i)>=0) cycle  ! don't count the surface and atm fields here
!!
!!           if (ngnd < 9) write(zsoil,'(I1)') ngnd+1
!!
!!           if ( (trim(fields(i)) == trim(vargnd(1)) .or. trim(fields(i)) == trim(vargnd(1))//'L'//zsoil) .and. &
!!                idates(i) == ndate .and. ifhrs(i) == nhr ) then
!!              ngnd=ngnd+1
!!              islevs(ngnd)=levs(i)
!!           endif
!!        enddo
!!
!!        if (ngnd > 0) then
!!
!!           ! Make sure soil depths monotonically decrease
!!
!!           allocate (itlevs(ngnd))
!!           itlevs=-99999
!!           do i=1,ngnd
!!              do j=1,ngnd
!!                 if(islevs(j) > itlevs(i))then
!!                    itlevs(i)=islevs(j)
!!                    ia=j
!!                 endif
!!              enddo
!!              !  Effectively remove current highest level
!!              islevs(ia)=-99999
!!           enddo
!!           islevs(1:ngnd)=itlevs(1:ngnd)
!!           deallocate (itlevs)
!!
!!        else
!!
!!           ! Check for ECMWF style names
!!
!!        endif
!!
!!     endif
!!
!!     print*
!!     print '(1x,I0,A)', ngnd, " Soil layers found:"
!!     print*, islevs(1:ngnd)
!!     print*

     read(c_date(1:4) ,'(i4)')iyyyy
     read(c_date(5:6) ,'(i2)')imm
     read(c_date(7:8) ,'(i2)')idd
     read(c_date(9:10),'(i2)')ihh

     tinc=real(nhr)
     call date_add_to(iyyyy,imm,idd,ihh*10000,tinc,funit,nuyyyy,numm,nudd,nuhh)
     nuhh=nuhh/100

!    allocate (a4(nx,ny,ngnd,nvargnd))
     allocate (a3(nx,ny,nlev,nvar3d))
     allocate (a2(nx,ny,nvarsfc))
   
     ! Fill the 3 and 2-d variable arrays
   
     do i = 1, nvar3d
        do j = 1, nlev
           ! print*,'nnn:',var3d(i),nx,ny,j,i,a3(1,1,j,i)
           a3(1:nx,1:ny,j,i) = amiss
           call getfield(filein,nx,ny,a3(:,:,j,i),var3d(i),iplevs(j),ndate,nhr,ircode)
           ! print*,'nnn:',var3d(i),nx,ny,j,i,a3(1,1,j,i)
        enddo
     enddo

     do i = 1, nvarsfc
        a2(1:nx,1:ny,i) = amiss
        j=0 !specify j=0 for the surface
        call getfield(filein,nx,ny,a2(:,:,i),varsfc(i),j,ndate,nhr,ircode)
     enddo
     
!!     do i = 1, nvargnd
!!        do j = 1, ngnd
!!           ! print*,'nnn:',vargnd(i),nx,ny,j,i,a4(1,1,j,i)
!!           a4(1:nx,1:ny,j,i) = amiss
!!           call getfield(filein,nx,ny,a4(:,:,j,i),vargnd(i),islevs(j),ndate,nhr,ircode)
!!           ! print*,'nnn:',vargnd(i),nx,ny,j,i,a4(1,1,j,i)
!!        enddo
!!     enddo

     ! Check for SST variable

     have_sst = .false.
     if (len_trim(sstvar) > 0) then
        allocate(sst(nx,ny))
        call getfield(filein, nx, ny, sst, sstvar, 0, ndate, nhr, ircode)
        if (ircode == 1) then
           have_sst = .true.
           call prepfield(sst, nx, ny, projection, 'SST')
           where (sst >= -100.0) sst = max(sst, t00sea)
        else
           deallocate(sst)
        endif
     endif

     ! Check for SEAICE variable

     have_ice = .false.
     if (len_trim(icevar) > 0) then
        allocate(ice(nx,ny))
        call getfield(filein, nx,ny, ice, icevar, 0, ndate, nhr, ircode)
        if (ircode == 1) then
           have_ice = .true.
           call prepfield(ice, nx, ny, projection, 'ICEC')
           where (ice >= -100.0) ice = min( max(ice, 0.0), 1.0 )
        else
           deallocate(ice)
        endif
     endif

     ! Check for LANDMASK variable

     if (have_sst .or. have_ice) then
        allocate(landf(nx,ny))
        landf = amiss
        if (len_trim(landmask) > 0) then
           call getfield(filein, nx,ny, landf, landmask, 0, ndate, nhr, ircode)
           if (ircode == 1) then
              call prepfield(landf, nx, ny, projection, 'LAND')
           endif
        endif
     endif

     ! Check for SOIL WATER variable

     have_soilw = .false.
     nb = len_trim(soilw_var)
     nsoilw = 0

     if (nb > 0) then

        if (soilw_var(nb:nb) == '#') then
           do n = 1, maxgnd
              soilw_vars(n) = soilw_var(1:nb-1) // anum(n)
           enddo
        else
           do n = 1, maxgnd
              soilw_vars(n) = soilw_var
           enddo
        endif

        do n = 1, maxgnd
           do i = 1, nrec
              if (levs(i) >= 0) cycle  ! don't count the surface and atm fields here

              if (fields(i) == soilw_vars(n) .and. idates(i) == ndate .and. ifhrs(i) == nhr) then

                 ! skip levels we already found
                 if (n > 1) then
                    if (any(iswlevs(1:n-1) == levs(i))) cycle
                 endif

                 nsoilw = nsoilw + 1
                 iswlevs(nsoilw) = levs(i)
                 exit
              endif

           enddo
           
           if (n > nsoilw) exit
        enddo

        if (nsoilw > 0) then
           have_soilw = .true.

           ! Make sure soil depths monotonically decrease if we don't specify the levels

           if (soilw_var(nb:nb) /= '#') then
              allocate (itlevs(nsoilw))
              itlevs=-99999
              do i=1,nsoilw
                 do j=1,nsoilw
                    if(iswlevs(j) > itlevs(i))then
                       itlevs(i)=iswlevs(j)
                       ia=j
                    endif
                 enddo
                 !  Effectively remove current highest level
                 iswlevs(ia)=-99999
              enddo
              iswlevs(1:nsoilw)=itlevs(1:nsoilw)
              deallocate (itlevs)
           endif

           allocate(soilw(nx,ny,nsoilw))
           soilw = amiss

           do n = 1, nsoilw
              call getfield(filein,nx,ny,soilw(:,:,n),soilw_vars(n),iswlevs(n),ndate,nhr,ircode)
              if (ircode == 1) then
                 call prepfield(soilw(:,:,n), nx, ny, projection, 'SOILW')
              else
                 soilw(:,:,n) = amiss
              endif
           enddo
        endif

     endif

     ! Check for SOIL TEMPERATURE variable

     have_soilt = .false.
     nb = len_trim(soilt_var)
     nsoilt = 0

     if (nb > 0) then

        if (soilt_var(nb:nb) == '#') then
           do n = 1, maxgnd
              soilt_vars(n) = soilt_var(1:nb-1) // anum(n)
           enddo
        else
           do n = 1, maxgnd
              soilt_vars(n) = soilt_var
           enddo
        endif

        do n = 1, maxgnd
           do i = 1, nrec
              if (levs(i) >= 0) cycle  ! don't count the surface and atm fields here

              if (fields(i) == soilt_vars(n) .and. idates(i) == ndate .and. ifhrs(i) == nhr) then

                 ! skip levels we already found
                 if (n > 1) then
                    if (any(istlevs(1:n-1) == levs(i))) cycle
                 endif

                 nsoilt = nsoilt + 1
                 istlevs(nsoilt) = levs(i)
                 exit
              endif

           enddo
           
           if (n > nsoilt) exit
        enddo

        if (nsoilt > 0) then
           have_soilt = .true.

           ! Make sure soil depths monotonically decrease if we don't specify the levels

           if (soilt_var(nb:nb) /= '#') then
              allocate (itlevs(nsoilt))
              itlevs=-99999
              do i=1,nsoilt
                 do j=1,nsoilt
                    if(istlevs(j) > itlevs(i))then
                       itlevs(i)=istlevs(j)
                       ia=j
                    endif
                 enddo
                 !  Effectively remove current highest level
                 istlevs(ia)=-99999
              enddo
              istlevs(1:nsoilt)=itlevs(1:nsoilt)
              deallocate (itlevs)
           endif

           allocate(soilt(nx,ny,nsoilt))
           soilt = amiss

           do n = 1, nsoilt
              call getfield(filein,nx,ny,soilt(:,:,n),soilt_vars(n),istlevs(n),ndate,nhr,ircode)
              if (ircode == 1) then
                 call prepfield(soilt(:,:,n), nx, ny, projection, 'SOILT')
              else
                 soilt(:,:,n) = amiss
              endif
           enddo
        endif

        if (have_soilt .and. allocated(soilw) .and. nsoilw > 0) then
           do n = 1, nsoilt
              where (soilw(:,:,1) < amiss0) soilt(:,:,n) = amiss
           enddo
        endif

     endif

     ! Check for SNOW variable

     have_snow = .false.
     if (len_trim(snow_var) > 0) then
        allocate(snow(nx,ny))
        call getfield(filein, nx, ny, snow, snow_var, 0, ndate, nhr, ircode)
        if (ircode == 1) then
           have_snow = .true.
           call prepfield(snow, nx, ny, projection, snow_var)
           where (snow >= amiss0) snow = max(snow, 0.0)
        else
           deallocate(snow)
        endif
     endif

     ! Relax SST using a poison solver over land areas to get a continuous field

     if (have_sst) then

        allocate(land(nx,ny))

        where ( (landf > -100.0 .and. landf >= 0.5) .or. sst < -100.0 )
           land = 1
        elsewhere
           land = 0
        endwhere

        ! Initialize SST over land for relaxation

        do j = 1, ny

           nsea = count( land(:,j) == 0 )

           if (nsea == 0) then
              ! If all longitudes are land, we are probably over Antartica
              ! Initialize SST to minimum value
              sst(:,j) = t00sea
           elseif (nsea < nx) then
              ! Initialize SST over land to the average SST at that latitude
              aver = sum(sst(:,j), mask=(land(:,j)==0)) / real(nsea)
              where (land(:,j) /= 0) sst(:,j) = aver
           endif

        enddo

        call relax(nx,ny,sst,land)
        sst = max(sst, t00sea)

     endif

     ! Relax SEAICE using a poison solver over land areas to get a continuous field

     if (have_ice) then

        if (.not. allocated(land)) allocate(land(nx,ny))

        where ( (landf > -100.0 .and. landf >= 0.5) .or. ice < -100.0 )
           land = 1
        elsewhere
           land = 0
        endwhere

        ! Initialize SEAICE over land for relaxation

        do j = 1, ny

           nsea = count( land(:,j) == 0 )

           if (nsea == 0) then
              ! If all longitudes are land, we are probably over Antartica
              ! Initialize seaice concentration to 1
              ice(:,j) = 1.0
           elseif (nsea < nx) then
              ! Initialize seacie over land to the average at that latitude
              aver = sum(ice(:,j), mask=(land(:,j)==0)) / real(nsea)
              where (land(:,j) /= 0) ice(:,j) = aver
           endif

        enddo

        call relax(nx,ny,ice,land)

        ! Set ice fraction values less than icefmin to 0
        where (ice < icefmin) ice = 0.0

        ! Make sure maximum ice fraction is 1.0 after relaxation
        ice = min(ice, 1.0)

     endif

     if (allocated(land)) deallocate(land)

     !call plotfield(nx,ny,a2(:,:,iice)*10.)
     !call plotfield(nx,ny,sst-273.15)

     ! Run the arrays through prepfield for possible conversion and flipping
     do i=1,nvar3d
        do j=1,nlev
           call prepfield(a3(:,:,j,i),nx,ny,projection,var3d(i))
        enddo
     enddo
   
     do i=1,nvarsfc
        call prepfield(a2(:,:,i),nx,ny,projection,varsfc(i))
     enddo

!!     do i=1,nvargnd
!!        do j=1,ngnd
!!           call prepfield(a4(:,:,j,i),nx,ny,projection,vargnd(i))
!!        enddo
!!     enddo
     
     ngnd = 0
     if (nsoilt > 0 .and. nsoilw == nsoilt) then
        ngnd = nsoilt
        islevs(1:ngnd) = istlevs(1:ngnd)
     endif

     !-----------------------------------------------------
     ! Output file name
     !-----------------------------------------------------
   
     write(dstring,'(i4.4,a1,i2.2,a1,i2.2,a1,i4.4)') &
          nuyyyy,'-',numm,'-',nudd,'-',nuhh
     fileout = trim(out_prefix)//trim(dstring)//trim(out_ext)
   
     ! Write HDF5

     print*
     print*, "Writing HDF5 file: " // trim(fileout)

     call shdf5_open(fileout,'W',1)
   
     gdf_file_ver=3
     ndims=1 ; idims(1)=1
     call shdf5_orec(ndims,idims,'version',ivars=gdf_file_ver)
     call shdf5_orec(ndims,idims,'year'   ,ivars=iyyyy)     
     call shdf5_orec(ndims,idims,'month'  ,ivars=imm)    
     call shdf5_orec(ndims,idims,'day'    ,ivars=idd)     
     call shdf5_orec(ndims,idims,'hour'   ,ivars=ihh*100)     
     call shdf5_orec(ndims,idims,'ftime'  ,ivars=nhr*100)     
     call shdf5_orec(ndims,idims,'nx'     ,ivars=nx)  
     call shdf5_orec(ndims,idims,'ny'     ,ivars=ny)  
     call shdf5_orec(ndims,idims,'nlev'   ,ivars=nlev)  
     call shdf5_orec(ndims,idims,'iproj'  ,ivars=inproj)  
     call shdf5_orec(ndims,idims,'vcoord' ,ivars=ivertcoord)  
     call shdf5_orec(ndims,idims,'swlat'  ,rvars=alat1)        
     call shdf5_orec(ndims,idims,'swlon'  ,rvars=alon1)       
     call shdf5_orec(ndims,idims,'nelat'  ,rvars=xnelat)        
     call shdf5_orec(ndims,idims,'nelon'  ,rvars=xnelon)       
     call shdf5_orec(ndims,idims,'dx'     ,rvars=dx)    
     call shdf5_orec(ndims,idims,'dy'     ,rvars=dx)     
     call shdf5_orec(ndims,idims,'reflat1',rvars=reflat1)        
     call shdf5_orec(ndims,idims,'reflat2',rvars=reflat2)        
     call shdf5_orec(ndims,idims,'ngnd'   ,ivars=ngnd)

     if (allocated(alats)) then
        ndims=1 ; idims(1)=ny
        call shdf5_orec(ndims, idims, 'glat', rvara=alats)
     endif

     if (nlev > 0) then
        ndims=1 ; idims(1)=nlev
        call shdf5_orec(ndims, idims, 'levels', ivara=iplevs)
     endif

     if (ngnd > 0) then
        ndims=1 ; idims(1)=ngnd
        call shdf5_orec(ndims, idims, 'sdepths', ivara=islevs)
     endif

     do i = 1, nvar3d
        write(*,*) out3d(i)
        ndims=3 ; idims(1)=nx; idims(2)=ny; idims(3)=nlev
        call shdf5_orec(ndims, idims, out3d(i), rvara=a3(:,:,:,i))
     enddo

     do i = 1, nvarsfc
        write(*,*) outsfc(i)
        ndims=2 ; idims(1)=nx; idims(2)=ny
        call shdf5_orec(ndims, idims, outsfc(i), rvara=a2(:,:,i))
     enddo

     if (have_sst .and. allocated(sst)) then
        write(*,*) 'SST'
        ndims=2 ; idims(1)=nx; idims(2)=ny
        call shdf5_orec(ndims, idims, 'SST', rvara=sst)
        deallocate(sst)
     endif

     if (have_ice .and. allocated(ice)) then
        write(*,*) 'ICEC'
        ndims=2 ; idims(1)=nx; idims(2)=ny
        call shdf5_orec(ndims, idims, 'ICEC', rvara=ice)
        deallocate(ice)
     endif

     if (have_soilw .and. allocated(soilw) .and. ngnd > 0) then
        write(*,*) 'SOILW'
        ndims=3 ; idims(1)=nx; idims(2)=ny; idims(3)=ngnd
        call shdf5_orec(ndims, idims, 'SOILW', rvara=soilw)
     endif

     if (have_soilt .and. allocated(soilt) .and. ngnd > 0) then
        write(*,*) 'SOILT'
        ndims=3 ; idims(1)=nx; idims(2)=ny; idims(3)=ngnd
        call shdf5_orec(ndims, idims, 'SOILT', rvara=soilt)
     endif

     if (have_snow .and. allocated(snow)) then
        write(*,*) 'SNOWMASS'
        ndims=2 ; idims(1)=nx; idims(2)=ny
        call shdf5_orec(ndims, idims, 'SNOWMASS', rvara=snow)
     endif


     call shdf5_close()

  enddo hour_loop

  !call clsgks

contains

!************************************************************************

subroutine getfield(filein,nx,ny,a,var,iplev,ndate,nhr,ircode)

  use grib_get_mod
  implicit none

  integer,      intent(in)    :: nx, ny, iplev, ndate, nhr
  real,         intent(inout) :: a(nx,ny)
  integer,      intent(out)   :: ircode
  character(*), intent(in)    :: filein, var
  integer                     :: i

  ircode = 0

  do i = 1, nrec

     if ( var   == fields(i) .and. iplev ==  levs(i) .and. &
          ndate == idates(i) .and. nhr   == ifhrs(i) ) then

        call grib_get_recf(filein, irecs(i), a, nx, ny)
        
        write(*, '(a,a,i8,i14,i8)') ' Found record: ', &
             trim(fields(i)), levs(i), idates(i), ifhrs(i)

        ircode = 1
        exit
     endif

  enddo
  
  if (ircode == 0) then
     write(*,*)
     write(*,*) "Level ", iplev
     write(*,*) "WARNING: Could not find " // trim(var) // " in grib file."
     write(*,*) "It will be set to missing values in the output file."
     write(*,*)
  endif

end subroutine getfield

!************************************************************************

subroutine prepfield(a,nx,ny,projection,type)

  use grib_get_mod, only: iscan, amiss, amiss0
  implicit none

  integer,      intent(in)    :: nx, ny
  real,         intent(inout) :: a(nx,ny)
  character(*), intent(in)    :: type, projection

  real,         allocatable   :: b(:,:)
  real,         parameter     :: g = 9.80616
  integer                     :: i, j

  if ( (projection=='latlon' .or. projection=='gaussian') .and. iscan==0) then
     write(*,*) 'flipping the projection'
     allocate (b(nx,ny))
     call flip(nx,ny,a,b)
     a=b
     deallocate (b)
  endif

  ! WGRIB sets missing/undefined values to 9.999e20
  where(a >= 1.e20) a = amiss

  if (type=='RH' .or. type=='R') then

     ! Assuming RH in %
     write(*,*) 'Converting RH from percent to fraction'
     where (a > amiss0) a = a / 100.

  elseif (type=='SPFH' .or. type=='Q') then

     write(*,*) 'Converting spec hum to g/kg'
     where (a > amiss0) a = a * 1000.

  elseif (type=='PRMSL' .or. type=='PRES') then

     write(*,*) 'Converting pressure to mb'
     where (a > amiss0) a = a / 100.

  elseif (type=='Z') then

     write(*,*) 'Factoring geopotential to height using g=', g
     where (a > amiss0) a = a / g

  elseif (type=='SOILW') then

     write(*,*) 'Setting negative soil water to missing values'
     where (a <= 0.0) a = amiss

  elseif (type=='SD') then

     ! ECMWF snow depth is in m of water equivalent
     write(*,*) 'Converting snow water depth to kg/m^2'
     where (a > amiss0) a = min(a * 1000., 175.0)

  endif

  !call ezcntr(a,nx,ny)

end subroutine prepfield

!************************************************************************

subroutine flip(nx,ny,vin,vout)
  implicit none

  integer, intent( in) :: nx,ny
  real,    intent( in) :: vin (nx,ny)
  real,    intent(out) :: vout(nx,ny)
  integer              :: i,j,jj

  do j=1,ny
     jj=ny-j+1
     do i=1,nx
        vout(i,j)=vin(i,jj)
     enddo
  enddo

end subroutine flip

!************************************************************************

subroutine date_add_to (inyear,inmonth,indate,inhour,  &
                        tinc,tunits,outyear,outmonth,outdate,outhour)

  implicit none

  ! add a time increment to a date and output new date
  ! -> uses hhmmss for hours, 4 digit year

  integer,      intent(in)  :: inyear, inmonth, indate, inhour
  real,         intent(in)  :: tinc
  character(1), intent(in)  :: tunits
  integer,      intent(out) :: outyear,outmonth,outdate,outhour

  integer,      parameter   :: mondays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
  real                      :: ttinc, strtim, xhourin, xminin, xsecin, xround
  integer                   :: izhours, izmin, izsec, iadddays, idays
  logical                   :: isleap

  ! convert input time to seconds

  ttinc = tinc
  if     (tunits=='m') then
     ttinc = tinc * 60.
  elseif (tunits=='h') then
     ttinc = tinc * 3600.
  elseif (tunits=='d') then
     ttinc=tinc * 86400.
  endif

  xhourin = inhour/10000
  xminin = mod(inhour,10000)/100
  xsecin = mod(inhour,100)
  strtim = xhourin + xminin/60. + xsecin/3600.

  izhours = int(mod(strtim+ttinc/3600.,24.)+.001)
  izmin   = int(mod(strtim+ttinc/3600.,1.)*60+.001)
  izsec   = int(mod(strtim*3600.+ttinc,60.)+.001)

  outhour = izhours*10000 + izmin*100 + izsec

  if (ttinc < 0.) then
     xround =-.001
  else
     xround = .001
  endif
     
  iadddays = int((strtim+ttinc/3600.)/24.+xround)

  outyear  = inyear
  outdate  = indate + iadddays
  outmonth = inmonth

  idays = mondays(outmonth)
  isleap = ( (mod(outyear,400) == 0) .or. &
             (mod(outyear,100) /= 0 .and. mod(outyear,4) == 0) )
  if (outmonth==2 .and. isleap) idays=29

  do while (outdate > idays .or. outdate < 1)

     if (outdate > idays) then

        outdate  = outdate - idays
        outmonth = outmonth + 1
        if (outmonth > 12) then
           outyear  = outyear + 1
           outmonth = 1
        endif

     elseif (outdate < 1) then

        if (outmonth==1) outmonth = 13
        idays = mondays(outmonth-1)
        isleap = ( (mod(outyear,400) == 0) .or. &
                   (mod(outyear,100) /= 0 .and. mod(outyear,4) == 0) )
        if (outmonth-1==2 .and. isleap) idays = 29
        outdate = idays+outdate
        outmonth = outmonth-1
        if (outmonth==12) outyear = outyear-1

     endif

     idays = mondays(outmonth)
     isleap = ( (mod(outyear,400) == 0) .or. &
                (mod(outyear,100) /= 0 .and. mod(outyear,4) == 0) )
     if (outmonth==2 .and. isleap) idays = 29

  enddo

end subroutine date_add_to

end program grib_to_gdf
