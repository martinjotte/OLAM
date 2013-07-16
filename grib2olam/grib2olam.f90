!
! Copyright (C) 1991-2003  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

program grib_to_gdf

  use grib_get_mod
  use hdf5_utils

  implicit none

  integer, parameter   :: max_hours=24
  character(len= 10)   :: c_date, c_hr(max_hours), out_ext
  character(len=128)   :: filein,fileout,out_prefix,dstring,namein
  character(len= 20)   :: cargv
  integer, allocatable :: itlevs(:)
  real,    allocatable :: a4(:,:,:,:),a3(:,:,:,:),a2(:,:,:)

  integer :: ia,iargc,istr,iend,inc
  integer :: i,j,idprep=0,ndate,nhr,ii,nh,ngnd
  integer :: nlev,inproj,ivertcoord,miss,ircode,num_hours
  
  integer :: iyyyy,imm,idd,ihh,nuyyyy,numm,nudd,nuhh
  real :: tinc, amiss=-999.,xnelat,xnelon
  logical :: exists, dousage

  integer, parameter :: maxlev=200, maxvar=100
  integer :: iplevs(maxlev), islevs(maxlev)
  integer :: nvar3d, nvarsfc, nvargnd, gdf_file_ver, ndims, idims(3)
  character(len=10) :: var3d(maxvar), varsfc(maxvar), vargnd(maxvar)
  character(len=10) :: out3d(maxvar), outsfc(maxvar), outgnd(maxvar)

  real, parameter   :: t00sea  = 271.38
  real, parameter   :: icefmin = 0.001
  real, allocatable :: landf(:,:), sst(:,:)
  integer, allocatable :: land(:,:)
  integer  :: iice

  namelist /dgrib_in/ nvar3d, var3d, nvarsfc, varsfc, c_date, num_hours, c_hr, &
                     filein, out3d, outsfc, wgrib1_exe, wgrib2_exe, out_prefix, &
                     nvargnd, vargnd, outgnd

  dousage = .false.
  out_ext = ".h5"
  namein  = "DGRIB_IN"

! Initialize the namelist variables

  wgrib1_exe        = "./wgrib"
  wgrib2_exe        = "./wgrib2"
  num_hours         = 1
  c_date            = "99999999"
  c_hr(1:max_hours) = "99999999"
  out_prefix        = ""
  var3d(:)          = ""
  varsfc(:)         = ""
  vargnd(:)         = ""
  out3d(:)          = ""
  outsfc(:)         = ""
  outgnd(:)         = ""
  nvar3d            = 0
  nvarsfc           = 0
  nvargnd           = 0
  filein            = ""

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

do i = 1, nvargnd
   if (len_trim(vargnd(i)) < 1) then
      write(*,*) " No input name given for gnd variable ", i
      stop
   endif
   if (len_trim(outgnd(i)) < 1) then
      write(*,*) " No output name given for gnd variable ", i
      stop
   endif
enddo

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
print '(a,a,i4)','Found projection: ',trim(projection), iscan
print '(a)', '  params:  lat1     lon1     nx     ny       dx          dy'
print '(a,2f8.2,2i7,2f12.4)','        ',alat1,alon1, nx,ny,dx,dy
print*

! Special for OLAM: We need global a dataset

if ((projection /= "latlon") .or. (nx*dx < 360.-dx) .or. (ny*dy < 180.-dy)) then
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
elseif (projection(1:7)=="Lambert") then
   inproj=2
   dx=dx*1000.
   dy=dy*1000.
   aorient=alov
elseif (projection(1:5)=="polar") then
   inproj=3
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
      if( trim(fields(i)) == trim(var3d(1)) .and.  &
               idates(i) == ndate .and. ifhrs(i) == nhr)then
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
      if(levs(i)<=0)cycle  !don't count the surface fields here
      if(trim(fields(i)) == trim(var3d(1)) .and.  &
              idates(i) == ndate .and. ifhrs(i) == nhr)then
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

   ! Determine the number of soil layers we are going to use

   ngnd=0
   do i=1, nrec
      if(levs(i)>=0)cycle  !don't count the surface and atm fields here
      if(trim(fields(i)) == trim(vargnd(1)) .and.  &
              idates(i) == ndate .and. ifhrs(i) == nhr) then
         ngnd=ngnd+1
         islevs(ngnd)=levs(i)
      endif
   enddo

   ! Make sure soil depths monotonically decrease

   allocate (itlevs(ngnd))
   itlevs=-99999
   do i=1,ngnd
      do j=1,ngnd
         if(islevs(j) > itlevs(i))then
            itlevs(i)=islevs(j)
            ia=j
         endif
      enddo
      !  Effectively remove current highest level
      islevs(ia)=-99999
   enddo
   islevs(1:ngnd)=itlevs(1:ngnd)
   deallocate (itlevs)

   print*
   print '(1x,I0,A)', ngnd, " Soil layers found:"
   print*, islevs(1:ngnd)
   print*

   read(c_date(1:4) ,'(i4)')iyyyy
   read(c_date(5:6) ,'(i2)')imm
   read(c_date(7:8) ,'(i2)')idd
   read(c_date(9:10),'(i2)')ihh

   tinc=real(nhr)
   call date_add_to(iyyyy,imm,idd,ihh*10000,tinc,funit,nuyyyy,numm,nudd,nuhh)
   nuhh=nuhh/100

   allocate (a4(nx,ny,nlev,nvargnd))
   allocate (a3(nx,ny,nlev,nvar3d))
   allocate (a2(nx,ny,nvarsfc))
   
   ! Fill the 3 and 2-d variable arrays
   
   do i=1,nvar3d
      do j=1,nlev
         ! print*,'nnn:',var3d(i),nx,ny,j,i,a3(1,1,j,i)
         a3(1:nx,1:ny,j,i) = amiss
         call getfield(filein,  &
                       nx*ny,a3(1,1,j,i),var3d(i),iplevs(j),ndate,nhr,ircode)
         ! print*,'nnn:',var3d(i),nx,ny,j,i,a3(1,1,j,i)
      enddo
   enddo

   !---------------------------------
   ! The Iberinco ECMWF exception
   !if (nvar3d == 5 .and. var3d(6) == 'Z') then
   !   do j=1,nlev
   !      ! a3(1:nx,1:ny,j,i)=amiss ! Don't fill in missing data
   !      call getfield(trim(filein)  &
   !                   ,nx*ny,a3(1,1,j,4),var3d(6),iplevs(j),ndate,nhr,ircode)
   !      if (ircode == 1) then
   !         ! Convert Z to GH
   !         a3(1:nx,1:ny,j,4) = a3(1:nx,1:ny,j,4)/9.8
   !      else
   !         !a3(1:nx,1:ny,j,4)=amiss
   !      endif
   !      print*,'zzz:',var3d(6),ircode,j,a3(1,1,j,4)
   !   enddo
   !endif
   !---------------------------------
   
   do i=1,nvarsfc
      a2(1:nx,1:ny,i)=amiss
      j=0 !specify j=0 for the surface
      call getfield(filein,  &
                    nx*ny,a2(1,1,i),varsfc(i),j,ndate,nhr,ircode)
   enddo

   do i = 1, nvargnd
      do j = 1, ngnd
         ! print*,'nnn:',vargnd(i),nx,ny,j,i,a4(1,1,j,i)
         a4(1:nx,1:ny,j,i) = amiss
         call getfield(filein,  &
                       nx*ny,a4(1,1,j,i),vargnd(i),islevs(j),ndate,nhr,ircode)
         ! print*,'nnn:',vargnd(i),nx,ny,j,i,a4(1,1,j,i)
      enddo
   enddo

   allocate(landf(nx,ny))
   allocate(land(nx,ny))
   call getfield(filein, nx*ny, landf, 'LAND', 0, ndate, nhr, ircode)
   write(*,*) ircode
   land = nint(landf)
   deallocate(landf)

   allocate(sst(nx,ny))
   do i = 1, nvarsfc
      if (varsfc(i) == "TMP") then
         sst = max(a2(:,:,i),t00sea)
      elseif (varsfc(i) == "ICEC") then
         iice = i
      endif
   enddo

   ! Make sure ice fraction is between 0 and 1
   a2(:,:,iice) = min( max( a2(:,:,iice), 0.0 ), 1.0 )

   ! Initialize ice fraction over land to 0
   where (land(:,:) == 1) a2(:,:,iice) = 0.0

   ! The relaxation for SST/SEAICE assumes constant top and bottom rows and
   ! periodic east and west boundaries. Set SST/SEAICE along bottom row 
   ! to 100% seaice for relaxation, and assume the top row is defined 

   sst(:,1)      = t00sea
   a2 (:,1,iice) = 1.0

   call relax(nx,ny,sst,land)
   call relax(nx,ny,a2(:,:,iice),land)

   ! Set ice fractions less than icefmin to 0
   where (a2(:,:,iice) < icefmin) a2(:,:,iice) = 0.0

   ! Make sure max ice fraction is 1.0 after relaxation
   a2(:,:,iice) = min( a2(:,:,iice), 1.0 )

!   call plotfield(nx,ny,a2(:,:,iice)*10.)
!   call plotfield(nx,ny,sst-273.15)

   ! Run the arrays through prepfield for possible conversion and flipping
   do i=1,nvar3d
      do j=1,nlev
         call prepfield(a3(1,1,j,i),nx,ny,projection,var3d(i),iscan,amiss)
      enddo
   enddo
   
   do i=1,nvarsfc
      call prepfield(a2(1,1,i),nx,ny,projection,varsfc(i),iscan,amiss)
   enddo

   do i=1,nvargnd
      do j=1,ngnd
         call prepfield(a4(1,1,j,i),nx,ny,projection,vargnd(i),iscan,amiss)
      enddo
   enddo
   
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

   ndims=1 ; idims(1)=nlev
   call shdf5_orec(ndims,idims,'levels'  ,ivara=iplevs)
      
   ndims=1 ; idims(1)=ngnd
   call shdf5_orec(ndims,idims,'sdepths' ,ivara=islevs)

   do i=1,nvar3d
      write(*,*) i, trim(out3d(i))
      ndims=3 ; idims(1)=nx; idims(2)=ny; idims(3)=nlev
      call shdf5_orec(ndims,idims,out3d(i),rvara=a3(:,:,:,i))
   enddo

   do i=1,nvarsfc
      write(*,*) nvar3d+i, trim(outsfc(i))
      ndims=2 ; idims(1)=nx; idims(2)=ny
      call shdf5_orec(ndims,idims,outsfc(i),rvara=a2(:,:,i))
   enddo

   do i=1,nvargnd
      write(*,*) i+nvar3d+nvarsfc, trim(outgnd(i))
      ndims=3 ; idims(1)=nx; idims(2)=ny; idims(3)=ngnd
      call shdf5_orec(ndims,idims,outgnd(i),rvara=a4(:,:,:,i))
   enddo

   write(*,*) 1+nvar3d+nvarsfc+nvargnd, "SST"
   ndims=2 ; idims(1)=nx; idims(2)=ny
   call shdf5_orec(ndims,idims,"SST",rvara=sst)
      
   call shdf5_close()

enddo hour_loop

!!call clsgks

end program grib_to_gdf

!************************************************************************

subroutine getfield(filein,nxy,a,var,iplev,ndate,nhr,ircode)

  use grib_get_mod
  implicit none

  integer      :: iplev,ndate,nhr,nxy,ircode
  real         :: a(nxy)
  character(*) :: filein,var

  integer :: i

  ircode = 0

  do i=1,nrec

     if ( var == fields(i) .and. iplev == levs(i) .and.  &
          ndate == idates(i) .and. nhr == ifhrs(i) ) then
        ! print*,'call grib_get_rec:',trim(filein),nxy,a(1),trim(irecs(i))

        call grib_get_recf(filein, irecs(i), a, nxy)
        
        print '(a,a,i8,i14,i8)','Found record: ' &
             ,trim(fields(i)),levs(i),idates(i),ifhrs(i)

        ircode = 1
        exit
     endif

  enddo
  
  if (ircode == 0) then
     write(*,*)
     write(*,*) "WARNING: Could not find " // trim(var) // " in grib file."
     write(*,*) "It will be set to missing values in the output file."
     write(*,*)
  endif

end subroutine getfield

!************************************************************************

subroutine prepfield(a,nx,ny,projection,type,iscan,amiss)
implicit none

integer nx,ny
real, dimension(nx*ny) :: a
character(len=*) :: type, projection

real, allocatable :: b(:)
real :: amiss
real :: g=9.80616
integer :: iscan

integer i

if(projection(1:6) == 'latlon' .and. iscan == 0)then
   print*,'   flipping the projection'
   allocate (b(nx*ny))
   call flip(nx,ny,a,b)
   a=b
   deallocate (b)
endif

do i=1,nx*ny
   if(a(i) >= 1.e10)a(i)=amiss
enddo

if(type(1:1)=='R')then
   ! Assuming RH in %
   print*,'   factoring R'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/100.
   enddo
elseif(type(1:4)=='SPFH')then
   print*,'   converting spec hum to g/kg'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)*1000.
   enddo
elseif(type=='PRMSL'.or.type=='PRES')then
   print*,'   factoring to mb'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/100.
   enddo
elseif(type=='Z')then
   print*,'   factoring geopotential to height using g=',g
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/g
   enddo

!!!!       TO DO:      !!!!!!!!
!
!  If Z is used and geo height less than about -100, this gets 
!     confused with the amiss. A suggested change:
!
! if(a(i) > amiss*g) a(i)=a(i)/g
!
!  I don't like operating on missing values. The other option is to change
!     amiss, especially now that hdf5 files are used.


endif

!!call ezcntr(a,nx,ny)


return
end subroutine prepfield

!************************************************************************

subroutine flip(nx,ny,vin,vout)
  implicit none

  integer, intent( in) :: nx,ny
  real,    intent( in) :: vin (nx,ny)
  real,    intent(out) :: vout(nx,ny)
  integer              :: i,j,ii,jj

  do i=1,nx
   ! ii=nx/2+i
   ! if(ii.gt.nx)ii=i-nx/2
     ii=i
     do j=1,ny
        jj=ny-j+1
        vout(i,j)=vin(ii,jj)
     enddo
  enddo

end subroutine flip

!************************************************************************

subroutine date_add_to (inyear,inmonth,indate,inhour  &
                        ,tinc,tunits,outyear,outmonth,outdate,outhour)

! add a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

integer inyear,inmonth,indate,inhour  &
       ,outyear,outmonth,outdate,outhour
real tinc
character(len=1) :: tunits
dimension mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/

! convert input time to seconds

ttinc=tinc
if(tunits.eq.'m') ttinc=tinc*60.
if(tunits.eq.'h') ttinc=tinc*3600.
if(tunits.eq.'d') ttinc=tinc*86400.
!print*,'inc:',tinc,tunits,ttinc

xhourin=inhour/10000
xminin=mod(inhour,10000)/100
xsecin=mod(inhour,100)
strtim=xhourin+xminin/60.+xsecin/3600.
!print*,'strtim=',strtim

izhours=int(mod(strtim+ttinc/3600.,24.)+.001)
izmin  =int(mod(strtim+ttinc/3600.,1.)*60+.001)
izsec  =int(mod(strtim*3600.+ttinc,60.)+.001)
!print*,'izs=',izhours,izmin,izsec

outhour= izhours*10000+izmin*100+izsec

iround=.001
if(ttinc<0.) iround=-.001
iadddays=int((strtim+ttinc/3600.)/24.+iround)

outyear=inyear
outdate=indate+iadddays
outmonth=inmonth

20 continue
   idays=mondays(outmonth)
   if(outmonth==2.and.mod(outyear,4)==0)  idays=29

   if(outdate.gt.idays) then
      outdate=outdate-idays
      outmonth=outmonth+1
      if(outmonth.gt.12) then
         outyear=outyear+1
         outmonth=1
      endif
   elseif(outdate.lt.1) then
      if(outmonth.eq.1)outmonth=13
      idays=mondays(outmonth-1)
      if(outmonth-1.eq.2.and.mod(outyear,4).eq.0)  idays=29
      outdate=idays+outdate
      outmonth=outmonth-1
      if(outmonth.eq.12)outyear=outyear-1
   else
      goto 21
   endif

   goto 20

21 continue

!print*,'out stuff:',outyear,outmonth,outdate,outhour

end subroutine date_add_to
