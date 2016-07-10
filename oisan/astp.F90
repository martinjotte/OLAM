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

!===============================================================================
subroutine read_press_header(fform)

use isan_coms,  only: iyear, nprz, levpr, ivertcoord, secondlat, cntlon, &
                      cntlat, xnelon, xnelat, itinc, inproj, gdatdy, gdatdx, &
                      xswlat, xswlon, npry, nprx, ihh, idd, imm, &
                      iyy, isversion, marker, innpr, imonth, idate, ihour, &
                      ipoffset, glat
use misc_coms,  only: io6, iparallel
use hdf5_utils, only: shdf5_open, shdf5_irec
use mem_para,   only: myrank, nbytes_int, nbytes_real

#ifdef OLAM_MPI
use mpi
#endif

implicit none

character(len=3), intent(inout) :: fform

character(len=16) :: ext

logical :: exists
integer :: lv, n, ier, igloberr, ipry
integer :: ndims, idims(2)
integer :: bytes, isize

integer, allocatable :: buffer(:)

! Read the header of input pressure file.

write(io6,'(/,a,/)') 'Reading pressure gridded data header '//trim(innpr)

! Set flag for type of file:
!     "GDF" - text Ralph II file (none or .vfm extension)
!     "HD5"  - HDF5 GDF file ("h5" or "hdf5")
!                     2 or 3D field records (x,y) or (x,y,z/p)
!                       (1,1)   - (SW lon,lat) --
!                       (1,1,1) - (SW lon,lat,bottom z/p level)
!                     range of header info

! Find file name extension

ext = trim( innpr( index(innpr,'.',back=.true.)+1:) )

fform = "GDF"
if (ext == 'h5' .or. ext == 'hdf5' .or. ext == 'hdf' .or. &
    ext == 'H5' .or. ext == 'HDF5' .or. ext == 'HDF' )    &
    fform = 'HD5'

inquire(file=innpr, exist=exists)
if (.not. exists) then
   write(*,*) "read_press_header: Error opening analysis file:"
   write(*,*) trim(innpr)
   stop       " File does not exist."
endif

if (fform == 'GDF') then

   open(11, file=innpr, status="OLD", action="READ")
   read(11,*) marker,isversion
   if(marker.ne.999999) isversion=1

   if (isversion == 1) then
      rewind 11
      read(11,*) iyy,imm,idd,ihh,nprz,nprx,npry,xswlon,xswlat,gdatdx,gdatdy
      read(11,*) (levpr(n),n=1,nprz)
      inproj = 1
      ihh = ihh * 100
   elseif (isversion == 2) then
      write(io6,*) 'doing RALPH 2 format'
      read(11,*) iyy,imm,idd,ihh,itinc,nprz,nprx,npry
      read(11,*) inproj,gdatdx,gdatdy,xswlat,xswlon &
                ,xnelat,xnelon,cntlat,cntlon,secondlat
      read(11,*) ivertcoord,(levpr(lv),lv=1,nprz)
   endif

elseif (fform == 'HD5') then

#ifdef OLAM_MPI
   if (iparallel == 1) then
      bytes = 0
      isize = nbytes_int*12 + nbytes_real*9
      allocate( buffer( isize ) )
   endif
#endif

   if (myrank == 0) then

      call shdf5_open (innpr, 'R')

      ndims    = 1
      idims(1) = 1
      idims(2) = 1

      call shdf5_irec(ndims, idims, 'version',ivars=isversion)
      call shdf5_irec(ndims, idims, 'year'   ,ivars=iyy)
      call shdf5_irec(ndims, idims, 'month'  ,ivars=imm)
      call shdf5_irec(ndims, idims, 'day'    ,ivars=idd)
      call shdf5_irec(ndims, idims, 'hour'   ,ivars=ihh)
      call shdf5_irec(ndims, idims, 'ftime'  ,ivars=itinc)
      call shdf5_irec(ndims, idims, 'nx'     ,ivars=nprx)
      call shdf5_irec(ndims, idims, 'ny'     ,ivars=npry)
      call shdf5_irec(ndims, idims, 'nlev'   ,ivars=nprz)
      call shdf5_irec(ndims, idims, 'iproj'  ,ivars=inproj)
      call shdf5_irec(ndims, idims, 'vcoord' ,ivars=ivertcoord)
      call shdf5_irec(ndims, idims, 'swlat'  ,rvars=xswlat)
      call shdf5_irec(ndims, idims, 'swlon'  ,rvars=xswlon)
      call shdf5_irec(ndims, idims, 'nelat'  ,rvars=xnelat)
      call shdf5_irec(ndims, idims, 'nelon'  ,rvars=xnelon)
      call shdf5_irec(ndims, idims, 'dx'     ,rvars=gdatdx)
      call shdf5_irec(ndims, idims, 'dy'     ,rvars=gdatdy)
      call shdf5_irec(ndims, idims, 'reflat1',rvars=cntlat)
      call shdf5_irec(ndims, idims, 'reflat2',rvars=secondlat)
      
      idims(1) = nprz
      call shdf5_irec(ndims, idims, 'levels' ,ivara=levpr)

      if (inproj == 2) then
         if (allocated(glat)) deallocate(glat)
         allocate(glat(npry))

         idims(1) = npry
         call shdf5_irec(ndims, idims, 'glat' ,rvara=glat)
      endif

#ifdef OLAM_MPI
      if (iparallel == 1) then
         call MPI_Pack(isversion , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(iyy       , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(imm       , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(idd       , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(ihh       , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(itinc     , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(nprx      , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(npry      , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(nprz      , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(inproj    , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(ivertcoord, 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(xswlat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(xswlon    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(xnelat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(xnelon    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(gdatdx    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(gdatdy    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(cntlat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
         call MPI_Pack(secondlat , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
      endif
#endif

   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then

      call MPI_Bcast(buffer, isize, MPI_PACKED, 0, MPI_COMM_WORLD, ier)

      if (myrank /= 0) then

         bytes = 0
         call MPI_Unpack(buffer, isize, bytes, isversion , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, iyy       , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, imm       , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, idd       , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, ihh       , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, itinc     , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, nprx      , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, npry      , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, nprz      , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, inproj    , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, ivertcoord, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, xswlat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, xswlon    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, xnelat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, xnelon    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, gdatdx    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, gdatdy    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, cntlat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         call MPI_Unpack(buffer, isize, bytes, secondlat , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
         
      endif

      call MPI_Bcast(levpr, nprz, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

      if (inproj == 2) then
         if (myrank /= 0) then
            if (allocated(glat)) deallocate(glat)
            allocate(glat(npry))
         endif
         call MPI_Bcast(glat, npry, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
      endif

      deallocate(buffer)

   endif
#endif

endif

write(io6,*) 'nprz1 ',nprz,nprx,npry

! Check for consistency between file parameters and namelist parameters

if (iyy /= iyear .or. imm /= imonth .or. idd /= idate .or. ihh /= ihour) then

   write(io6,*) 'Pressure file dates not the same as namelist!'
   write(io6,*) 'Year :',iyy,iyear
   write(io6,*) 'Month:',imm,imonth
   write(io6,*) 'Day  :',idd,idate
   write(io6,*) 'Hour :',ihh,ihour
   stop 'pr_dates'
endif

! Check pressure data domain size, location, and type

if (inproj == 1) then

   ! INPROJ = 1 denotes an input gridded atmospheric dataset defined on a
   ! latitude-longitude grid, with uniformly-spaced latitude and longitude,
   ! defined by parameters (nprx, npry, nprz, gdatdx, gdatdy, xswlat, xswlon)

   ! We make the requirement that a full global domain of pressure-level data
   ! be read in.  Check this here.  Following the convention for the NCEP/DOE
   ! Reanalysis2 data, it is assumed that nprx * gdatdx should equal 360 degrees
   ! (of longitude).  If this is not the case, this check will stop execution.

   igloberr = 0

   if (abs(nprx * gdatdx - 360.) > .1) igloberr = 1

   ! Data points may be defined at latitudinal coordinates that include both
   ! geographic poles (-90. and 90. degrees), in which (npry-1) * gdatdy should
   ! equal 180 degrees, or data points may be offset by 1/2 gdatdy from polar
   ! locations, in which case npry * gdatdy should equal 180 degrees.  Both
   ! possibilities are checked here, and if neither is satisfied, this check
   ! will stop execution.  For either case, the beginning latitude of the dataset
   ! is checked for consistency.

   if (abs((npry-1) * gdatdy - 180.) < .1) then
      if (abs(xswlat + 90.) > .1) igloberr = 1
   elseif (abs(npry * gdatdy - 180.) < .1) then
      if (abs(xswlat - 0.5 * gdatdy + 90.) > .1) igloberr = 1
   else
      igloberr = 1
   endif

   if (igloberr == 1) then
      write(io6,*) 'INPROJ = ',inproj
      write(io6,*) 'Gridded pressure level data must have global coverage'
      write(io6,*) 'nprx,npry = ',nprx,npry
      write(io6,*) 'gdatdx,gdatdy = ',gdatdx,gdatdy
      write(io6,*) 'xswlat,xswlon= ',xswlat,xswlon
      stop 'astp stop1 - non-global domain in input pressure data'
   endif

   ! Compute longitudinal offset index for copying input data to the expanded
   ! isan pressure arrays.

   ipoffset = int((xswlon + 180.) / gdatdx) + 2

elseif (inproj == 2) then

   ! INPROJ = 2 denotes an input gridded atmospheric dataset defined on a
   ! latitude-longitude grid, with uniformly-spaced longitude and variable
   ! latitudinal spacing, and defined by parameters (nprx, npry, nprz, gdatdx,
   ! xswlon) and glat, an array of specified latitudes.

   ! We make the requirement that a full global domain of pressure-level data
   ! be read in.  Check this here.  Following the convention for the NCEP/DOE
   ! Reanalysis2 data, it is assumed that nprx * gdatdx should equal 360 degrees
   ! (of longitude).  If this is not the case, this check will stop execution.

   igloberr = 0

   if (abs(nprx * gdatdx - 360.) > .1) igloberr = 1

   ! Data points may be defined at latitudinal coordinates that include both
   ! geographic poles (-90. and 90. degrees) or data points may be offset by
   ! (approximately) 1/2 gdatdy from polar locations.  Here, we check that at
   ! least one of these possibilities is satisfied.  If not, this check will
   ! stop execution.

   if (glat(1) - (glat(2) - glat(1)) > -90.) igloberr = 1
   if (glat(npry) + (glat(npry) - glat(npry-1)) < 90.) igloberr = 1

   if (igloberr == 1) then
      write(io6,*) 'INPROJ = ',inproj
      write(io6,*) 'Gridded pressure level data must have global coverage'
      write(io6,*) 'nprx,npry = ',nprx,npry
      write(io6,*) 'gdatdx = ',gdatdx
      write(io6,*) 'xswlat,xswlon= ',xswlat,xswlon
      do ipry = 1,npry
         write(io6,*) 'ipry, glat = ',ipry,glat(ipry)
      enddo
      stop 'astp stop2 - non-global domain in input pressure data'
   endif

   ! Compute longitudinal offset index for copying input data to the expanded
   ! isan pressure arrays.

   ipoffset = int((xswlon + 180.) / gdatdx) + 2

else

   write(io6,*) 'The input gridded atmospheric dataset does not conform '
   write(io6,*) 'to currently-implemented formats, which are: '
   write(io6,*) '(iproj=1) latitude-longitude with uniform spacing '
   write(io6,*) '(iproj=2) latitude-longitude with specified variable '
   write(io6,*) '          latitude and uniformly-spaced longitude '
   write(io6,*) ' '
   write(io6,*) 'Other formats will require additional coding.'

   stop 'astp stop - input atmospheric dataset format'

endif

end subroutine read_press_header

!===============================================================================

subroutine pressure_stage(fform, p_u, p_v, p_t, p_z, p_r, p_o, &
                          p_topo, p_prsfc, p_tsfc, p_shsfc)

use isan_coms,   only: pnpr, levpr, nprx, npry, nprz, nprz_rh
use consts_coms, only: rocp, p00, eps_vap
use misc_coms,   only: io6
use hdf5_utils,  only: shdf5_close
use mem_para,    only: myrank

implicit none

character(len=*), intent(in) :: fform

real, intent(inout) :: p_u(nprx+4,npry+4,nprz)
real, intent(inout) :: p_v(nprx+4,npry+4,nprz)
real, intent(inout) :: p_t(nprx+4,npry+4,nprz)
real, intent(inout) :: p_z(nprx+4,npry+4,nprz)
real, intent(inout) :: p_r(nprx+4,npry+4,nprz)
real, intent(inout) :: p_o(nprx+4,npry+4,nprz)

real, intent(inout) :: p_topo (nprx+4,npry+4)
real, intent(inout) :: p_prsfc(nprx+4,npry+4)
real, intent(inout) :: p_tsfc (nprx+4,npry+4)
real, intent(inout) :: p_shsfc(nprx+4,npry+4)

real :: thmax,thmin,vapor_press,qmax
integer :: i,j,k,iunit
logical :: isrh, doconvert

real, external :: eslf

iunit = 11

write(io6, *) ''
write(io6, *) '*****************************************************'
write(io6, *) '     Access pressure level data'
write(io6, *) '*****************************************************'

do k = 1,nprz
   pnpr(k) = levpr(k) * 100.
enddo

! Call routine to fill pressure arrays from the chosen dataset.

call get_press (fform, iunit, p_u, p_v, p_t, p_z, p_r, p_o, &
                p_topo, p_prsfc, p_tsfc, p_shsfc, isrh)

!!!!!!!! Be careful !!!!!!!!!
!  Check input humidity variable p_r.  Assume that if the maximum of the field
!  is greater than 1.1 (which allows for some machine roundoff),
!  it is specific humidity (in g/kg) which needs to be converted to kg/kg,
!  else it is R.H. which needs to be converted to specific humidity

if (.not. isrh) then

   ! Old RALPH files should have specific humidity in g/kg
   doconvert = .true.
   
   ! New HDF files should also be in g/kg, but if the grib file had a 
   ! nonstandard name for humidity, grib2olam might not have converted
   ! the units properly. Check the magnitude of humdity to guess the units.
   ! This will be changed when we write the units to the degribbed file!

   if (fform == 'HD5') then

      ! Level 1 should always the lowest level!
      qmax = maxval(p_r(1:nprx+4, 1:npry+4, 1))

      ! If all specific humidities at lowest level are less than 1,
      ! assume it is alread kg/kg
      if (qmax < 1.0) doconvert = .false.

   endif

   if (doconvert) then

      ! Convert specific humidity to kg/kg units

      write(io6, *) ''
      write(io6, *) 'Converting specific humidity to kg/kg'

      do k = 1,nprz
         do j = 1,npry+4
            do i = 1,nprx+4
               p_r(i,j,k) = .001 * p_r(i,j,k)
            enddo
         enddo
      enddo

   endif

else

   ! Convert R.H. to specific humidity

   write(io6, *) ''
   write(io6, *) 'Converting relative humidity to specific humidity'

   do k = 1,nprz
      do j = 1,npry+4
         do i = 1,nprx+4

            ! Compute ambient vapor pressure based on R.H.
            ! and saturation vapor pressure (eslf)

            vapor_press = p_r(i,j,k) * eslf(p_t(i,j,k)-273.15)

            ! Do not allow vapor pressure to exceed ambient pressure

            vapor_press = min(pnpr(k),vapor_press)

            ! Compute specific humidity from vapor press and ambient press

            p_r(i,j,k) = eps_vap * vapor_press &
                 / (pnpr(k) + vapor_press * (eps_vap - 1.))

         enddo
      enddo
   enddo

endif

! Print max-min theta at bottom and top levels

thmax = maxval(p_t(:,:,nprz)) * (p00/pnpr(nprz))**rocp
thmin = minval(p_t(:,:,1   )) * (p00/pnpr(1))**rocp

write(io6, *) ''
write(io6, "(' Minimum THETA at ',I4,' mb: ',F9.3)") levpr(1), thmin
write(io6, "(' Maximum THETA at ',I4,' mb: ',F9.3)") levpr(nprz), thmax

if (fform == 'GDF') then
   close(iunit)
elseif (fform == 'HD5' .and. myrank == 0) then
   call shdf5_close()
endif

write(io6,*) ''
write(io6,*) 'nprz2 ',nprz

end subroutine pressure_stage

!===============================================================================

subroutine get_press (fform, iunit, p_u, p_v, p_t, p_z, p_r, p_o, &
                      p_topo, p_prsfc, p_tsfc, p_shsfc, isrh)

use max_dims,   only: maxpr
use isan_coms,  only: nprz, npry, nprx, pnpr, iyear, imonth, idate, &
                      ihour, levpr, nprz_rh, ihydsfc, nbot_o3, haso3
use misc_coms,  only: io6, iparallel
use hdf5_utils, only: shdf5_irec, shdf5_info
use mem_para,   only: myrank

#ifdef OLAM_MPI
use mpi
#endif

implicit none

character(len=*), intent(in) :: fform
integer, intent(in)  :: iunit
logical, intent(inout) :: isrh

real, intent(inout) :: p_u(nprx+4,npry+4,nprz)
real, intent(inout) :: p_v(nprx+4,npry+4,nprz)
real, intent(inout) :: p_t(nprx+4,npry+4,nprz)
real, intent(inout) :: p_z(nprx+4,npry+4,nprz)
real, intent(inout) :: p_r(nprx+4,npry+4,nprz)
real, intent(inout) :: p_o(nprx+4,npry+4,nprz)

real, intent(inout) :: p_topo (nprx+4,npry+4)
real, intent(inout) :: p_prsfc(nprx+4,npry+4)
real, intent(inout) :: p_tsfc (nprx+4,npry+4)
real, intent(inout) :: p_shsfc(nprx+4,npry+4)

real :: as(nprx,npry)
real :: as3(nprx,npry,nprz)

integer :: i,j,k,lv,n,ier
integer :: ithere(maxpr,5),isfthere(5)
character(len=1) :: idat(5) = (/ 'T','U','V','H','R' /)
integer :: ndims, idims(3), njdims, jdims(3)
character(10) :: varname

! Set ihydsfc = 1 to read surface variables (topography, pressure, temperature, 
! moisture) and to compute geopotential heights of pressure levels by
! hydrostatic integration from surface; set ihydint = 0 otherwise.

ihydsfc = 0

! Initialize with missing data flag

ithere   = -999
isfthere = -999
haso3    = .false.

!  Read upper air fields

isrh = .true.
write(io6,*) ' '

if (fform == 'GDF') then

   do lv = 1,nprz

! Zonal wind component

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,as,p_u(1,1,lv))
      write(io6, '('' ==  Read UE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Meridional wind component

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,as,p_v(1,1,lv))
      write(io6, '('' ==  Read VE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Temperature

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,as,p_t(1,1,lv))
      write(io6, '('' ==  Read T on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear
      write(io6,*) 'prread1 ',as(5,5),nprx,npry,p_t(5,5,lv)

! Geopotential height

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,as,p_z(1,1,lv))
      write(io6, '('' ==  Read Z on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Relative humidity (or specific humidity)

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,as,p_r(1,1,lv))
      write(io6, '('' ==  Read RH on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

   enddo

   write(io6,*) ''

   goto 71

   70 CONTINUE
   write(io6,*) 'Premature end of file or error in pressure input file!'
   write(io6,*) 'We''ll close our eyes and pretend it didn''t happen!'
   71 continue

elseif (fform == 'HD5') then

! Read 2D surface vars if ihydsfc = 1

   if (ihydsfc == 1) then

      ndims = 2
      idims(1) = nprx
      idims(2) = npry
      idims(3) = 1

   ! Read topography

      if (myrank == 0) then
         call shdf5_irec(ndims, idims,'TOPO',rvara = as)
         call prfill(nprx,npry,as,p_topo)
      endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_topo, (nprx+4)*(npry+4), MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read surface pressure

      if (myrank == 0) then
         call shdf5_irec(ndims, idims,'PRSFC',rvara = as)
         call prfill(nprx,npry,as,p_prsfc)
      endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_prsfc, (nprx+4)*(npry+4), MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read surface temperature

      if (myrank == 0) then
         call shdf5_irec(ndims, idims,'TSFC',rvara = as)
         call prfill(nprx,npry,as,p_tsfc)
      endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_tsfc, (nprx+4)*(npry+4), MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read surface specific humidity

      if (myrank == 0) then
         call shdf5_irec(ndims, idims,'SHSFC',rvara = as)
         call prfill(nprx,npry,as,p_shsfc)
      endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_shsfc, (nprx+4)*(npry+4), MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   endif  ! ihydsfc

! Read 3D met vars

   ndims = 3
   idims(1) = nprx
   idims(2) = npry
   idims(3) = nprz

   ! Read east-west velocity U. It may be called UP, UE, or U in the file

   if (myrank == 0) then

      varname = 'U'
      call shdf5_info(varname, njdims, jdims)

      if (njdims <= 0) then
         varname = 'UP'
         call shdf5_info(varname, njdims, jdims)
      endif

      if (njdims <= 0) then
         varname = 'UE'
         call shdf5_info(varname, njdims, jdims)
      endif

      call shdf5_irec(ndims, idims, varname, rvara = as3)
      call prfill3(nprx,npry,nprz,as3,p_u)

   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_u, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read north-south velocity V. It may be called VP, VE, or V in the file

   if (myrank == 0) then

      varname = 'V'
      call shdf5_info(varname, njdims, jdims)

      if (njdims <= 0) then
         varname = 'VP'
         call shdf5_info(varname, njdims, jdims)
      endif

      if (njdims <= 0) then
         varname = 'VE'
         call shdf5_info(varname, njdims, jdims)
      endif

      call shdf5_irec(ndims, idims, varname, rvara = as3)
      call prfill3(nprx,npry,nprz,as3,p_v)

   endif
   
#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_v, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read temperature 

   if (myrank == 0) then

      call shdf5_irec(ndims, idims,'TEMP',rvara = as3)
      call prfill3(nprx,npry,nprz,as3,p_t)

   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_t, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read geopotential height

   if (myrank == 0) then
      call shdf5_irec(ndims, idims,'GEO',rvara = as3)
      call prfill3(nprx,npry,nprz,as3,p_z)
   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_z, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read water vapor specific humidity if it exists in the file

   if (myrank == 0) then

      njdims = 0
      call shdf5_info("SHV", njdims, jdims)

      if (njdims > 0) then

         isrh = .false.
         call shdf5_irec(ndims, idims,'SHV',rvara = as3)
         call prfill3(nprx,npry,nprz,as3,p_r)
     
      else
      
   ! Read relative humidity if we don't have specific humidity, 
   ! It may be called RELHUM or RH

         varname = 'RH'
         call shdf5_info(varname, njdims, jdims)

         if (njdims <= 0) then
            varname = 'RELHUM'
         endif

         isrh = .true.
         call shdf5_irec(ndims, idims, varname, rvara = as3)
         call prfill3(nprx,npry,nprz,as3,p_r)

      endif
   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( p_r, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
      call MPI_Bcast( isrh, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier )
   endif
#endif

   ! Read ozone if it exists in the file

   if (myrank == 0) then

      njdims = 0
      call shdf5_info("O3MR", njdims, jdims)

      if (njdims > 0) then

         haso3 = .true.
         call shdf5_irec(ndims, idims, 'O3MR', rvara = as3)
         call prfill3(nprx,npry,nprz,as3,p_o)

      endif

   endif

#ifdef OLAM_MPI
   if (iparallel == 1) then
      call MPI_Bcast( haso3, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier )
         
      if (haso3) then
         call MPI_Bcast( p_o, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, ier )
      endif

   endif
#endif

!  if (ivertcoord == 3) call shdf5_irec('PRESS',rvara = p_p)

!  print*,'uu max-min:', maxval(p_u), minval(p_u)
!  print*,'vv max-min:', maxval(p_v), minval(p_v)
!  print*,'tt max-min:', maxval(p_t), minval(p_t)
!  print*,'zz max-min:', maxval(p_z), minval(p_z)
!  print*,'rr max-min:', maxval(p_r), minval(p_r)
!  print*,'pp max-min:', maxval(p_p), minval(p_p)

   ! 3D scalar vars

!   if (num_scvars > 0) then
!      do nv = 1, num_scvars
!         call shdf5_irec(trim(in_scvars(nv)),rvara = p_scalar(1,1,1,nv))

!         print*,trim(in_scvars(nv)),' max-min:', maxval(p_scalar(:,:,:,nv)), &
!                      minval(p_scalar(:,:,:,nv))

!      enddo
!   endif

endif

! Special for RH:
! Reanalysis typically only reports RH up to 100mb, but still reports the other
! fields up to higher levels. Check for the highest level that reports RH:

do k = nprz,1,-1
   if (all(p_r(:,:,k) > -998.)) then
      nprz_rh = k
      exit
   endif
enddo

! Special for GDF text files:
! The old ralph format does not differentiate between specific humidity
! and relative humidity. We assume that if the maximum of the humidity 
! variable is greater than 1.1 (which allows for some machine roundoff), 
! it is specific humidity (in g/kg), else it is relative humidity

if (fform == 'GDF') then
   if (maxval(p_r(1:nprx+4,1:npry+4,1:nprz_rh)) > 1.1) then
      isrh = .false.
   else
      isrh = .true.
   endif
endif

! Special for OZONE:
! GFS only reports ozone ABOVE 100 mb, whereas the CFSR reanalysis reports the
! whole column. Check the lowest level at which ozone is reported if it is in
! the analysis file:

if (haso3) then
   nbot_o3 = 0

   do k = 1, nprz
      if (all(p_o(:,:,k) > -998.)) then
         nbot_o3 = k
         exit
      endif
   enddo
   
   if (nbot_o3 == 0) haso3 = .false.
endif

write(io6, *) '----------------------------------------------------'
write(io6, "(A,I0,A,I0,A)") ' Pressure-level data has ', nprz, &
     ' levels, goes up to ', levpr(nprz), ' mb.'
write(io6, "(A,I0,A,I0,A)") ' Water vapor has ', nprz_rh, &
     ' levels, is reported up to ', levpr(nprz_rh), ' mb.'
write(io6, *) ''

! Check for missing data

do k = 1,nprz
   ithere(k,1) = count(p_t(:,:,k) < -998.0)
   ithere(k,2) = count(p_u(:,:,k) < -998.0)
   ithere(k,3) = count(p_v(:,:,k) < -998.0)
   ithere(k,4) = count(p_z(:,:,k) < -998.0)
   ithere(k,5) = count(p_r(:,:,k) < -998.0)
enddo

where (ithere(:,:) > nprx*npry) ithere(:,:) = -1

write(io6,*) '---------------------------------------------------'
write(io6,*) ' # of missing values per level (-1 = all missing): '
write(io6,*) '---------------------------------------------------'
write(io6, '(a,5(6x,a1))') '   P (mb)', (idat(n),n=1,5)

do k = 1,nprz
   write(io6, '(f10.2,t10,5(i7))') pnpr(k)/100., (ithere(k,n),n=1,5)
enddo

if (any(ithere(1:nprz,(/1,2,3,4/)) /= 0) .or. any(ithere(1:nprz_rh,5) /= 0)) then
   write(io6,*) ''
   write(io6,*) '-------------------------------------------------------'
   write(io6,*) 'WARNING - There appears to be missing data in the input'
   write(io6,*) '          pressure-level data files'
   write(io6,*) '-------------------------------------------------------'
endif

end subroutine get_press

!===============================================================================

subroutine prfill (nprx,npry,xx,dn)

use isan_coms,  only: gdatdy, xswlat, ipoffset, inproj

implicit none

integer, intent(in) :: nprx,npry

real, intent(in)    :: xx(nprx,npry)
real, intent(inout) :: dn(nprx+4,npry+4)

integer :: i,j,nprxo2,iin

nprxo2 = nprx / 2

! Copy pressure-level or surface data from xx to dn array, shifting by ipoffset
! in x direction and by 2 in y direction

do j = 1,npry
   do i = 1,nprx+4
      iin = mod(i+nprx-ipoffset-1,nprx)+1
      dn(i,j+2) = xx(iin,j)
   enddo
enddo

! Fill in N/S boundary rows, based on whether xswlat = -90. or -90. + gdatdy/2.

if (inproj==1 .and. abs(xswlat + 90.) < .1) then

   do i = 1,nprxo2
      dn(i,1)      = dn(i+nprxo2,5)
      dn(i,2)      = dn(i+nprxo2,4)
      dn(i,npry+3) = dn(i+nprxo2,npry+1)
      dn(i,npry+4) = dn(i+nprxo2,npry)
   enddo

   do i = nprxo2+1,nprx+4
      dn(i,1)      = dn(i-nprxo2,5)
      dn(i,2)      = dn(i-nprxo2,4)
      dn(i,npry+3) = dn(i-nprxo2,npry+1)
      dn(i,npry+4) = dn(i-nprxo2,npry)
   enddo

elseif (inproj==2 .or. (inproj==1 .and. abs(xswlat - .5 * gdatdy + 90.) < .1)) then

   do i = 1,nprxo2
      dn(i,1)      = dn(i+nprxo2,4)
      dn(i,2)      = dn(i+nprxo2,3)
      dn(i,npry+3) = dn(i+nprxo2,npry+2)
      dn(i,npry+4) = dn(i+nprxo2,npry+1)
   enddo

   do i = nprxo2+1,nprx+4
      dn(i,1)      = dn(i-nprxo2,4)
      dn(i,2)      = dn(i-nprxo2,3)
      dn(i,npry+3) = dn(i-nprxo2,npry+2)
      dn(i,npry+4) = dn(i-nprxo2,npry+1)
   enddo

else
   print*, 'xswlat & gdatdy combo not found ',xswlat,gdatdy
endif

end subroutine prfill

!===============================================================================

subroutine prfill3 (nprx,npry,nprz,xx,dn)

use isan_coms,  only: gdatdy, xswlat, ipoffset, inproj

implicit none

integer, intent(in) :: nprx,npry,nprz

real, intent(in)    :: xx(nprx,npry,nprz)
real, intent(inout) :: dn(nprx+4,npry+4,nprz)

integer :: i,j,k,nprxo2,iin

nprxo2 = nprx / 2

! Copy pressure-level or surface data from xx to dn array, shifting by ipoffset
! in x direction and by 2 in y direction

do k = 1,nprz
   do j = 1,npry
      do i = 1,nprx+4
         iin = mod(i+nprx-ipoffset-1,nprx)+1
         dn(i,j+2,k) = xx(iin,j,k)
      enddo
   enddo

! Fill in N/S boundary rows, based on whether xswlat = -90. or -90. + gdatdy/2.

   if (inproj==1 .and. abs(xswlat + 90.) < .1) then

      do i = 1,nprxo2
         dn(i,1,k)      = dn(i+nprxo2,5,k)
         dn(i,2,k)      = dn(i+nprxo2,4,k)
         dn(i,npry+3,k) = dn(i+nprxo2,npry+1,k)
         dn(i,npry+4,k) = dn(i+nprxo2,npry,k)
      enddo

      do i = nprxo2+1,nprx+4
         dn(i,1,k)      = dn(i-nprxo2,5,k)
         dn(i,2,k)      = dn(i-nprxo2,4,k)
         dn(i,npry+3,k) = dn(i-nprxo2,npry+1,k)
         dn(i,npry+4,k) = dn(i-nprxo2,npry,k)
      enddo

   elseif (inproj==2 .or. (inproj==1 .and. abs(xswlat - .5 * gdatdy + 90.) < .1)) then

      do i = 1,nprxo2
         dn(i,1,k)      = dn(i+nprxo2,4,k)
         dn(i,2,k)      = dn(i+nprxo2,3,k)
         dn(i,npry+3,k) = dn(i+nprxo2,npry+2,k)
         dn(i,npry+4,k) = dn(i+nprxo2,npry+1,k)
      enddo

      do i = nprxo2+1,nprx+4
         dn(i,1,k)      = dn(i-nprxo2,4,k)
         dn(i,2,k)      = dn(i-nprxo2,3,k)
         dn(i,npry+3,k) = dn(i-nprxo2,npry+2,k)
         dn(i,npry+4,k) = dn(i-nprxo2,npry+1,k)
      enddo

   else
      print*, 'xswlat & gdatdy combo not found ',xswlat,gdatdy
   endif
   
enddo

end subroutine prfill3
