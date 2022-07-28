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

subroutine read_press_header()

  use isan_coms,  only: iyear, nprz, levpr, ivertcoord, secondlat, cntlat, &
                        xnelon, xnelat, itinc, inproj, gdatdy, gdatdx, &
                        xswlat, xswlon, npry, nprx, ihh, idd, imm, iyy, &
                        isversion, innpr, imonth, idate, ihour, ipoffset, &
                        glat, pnpr, o3name, haso3
  use misc_coms,  only: io6, iparallel
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info
  use mem_para,   only: myrank, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  character(len=16) :: ext

  logical :: exists
  integer :: ier, igloberr, ipry, k
  integer :: ndims, idims(3)
  integer :: bytes, isize

  integer, allocatable :: buffer(:)

  ndims = 0
  idims = 0

  ! Read the header of input pressure file.

  write(io6,'(/,a,/)') 'Reading pressure gridded data header '//trim(innpr)

  ! Find file name extension

  ext = trim( innpr( index(innpr,'.',back=.true.)+1:) )

  inquire(file=innpr, exist=exists)
  if (.not. exists) then
     write(*,*) "read_press_header: Error opening analysis file:"
     write(*,*) trim(innpr)
     stop       " File does not exist."
  endif

  ! Read header only on processor rank 0

  if (myrank == 0) then
     call shdf5_open (innpr, 'R')

     ndims    = 1
     idims(1) = 1

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

     if (allocated(levpr)) deallocate(levpr)
     allocate(levpr(nprz))

     ndims    = 1
     idims(1) = nprz
     call shdf5_irec(ndims, idims, 'levels' ,ivar1=levpr)

     if (inproj == 2) then
        if (allocated(glat)) deallocate(glat)
        allocate(glat(npry))

        ndims    = 1
        idims(1) = npry
        call shdf5_irec(ndims, idims, 'glat' ,rvar1=glat)
     endif

     ! Check for ozone variable in file

     ndims = 0
     idims = 0
     haso3 = .false.

     call shdf5_info('O3MR', ndims, idims)

     if (ndims == 3 .and. all( idims == [nprx,npry,nprz] ) ) then

        haso3  = .true.
        o3name = 'O3MR'

     else

        call shdf5_info('OZONE', ndims, idims)

        if (ndims == 3 .and. all( idims == [nprx,npry,nprz] ) ) then
           haso3  = .true.
           o3name = 'OZONE'
        endif

     endif

  endif  ! myrank == 0

  ! If parallel run, communicate header data to other ranks

#ifdef OLAM_MPI
  if (iparallel == 1) then

     bytes = 0
     isize = nbytes_int*13 + nbytes_real*9
     allocate( buffer( isize ) )

     if (myrank == 0) then
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
        call MPI_Pack(haso3     , 1, MPI_LOGICAL, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlon    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xnelat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xnelon    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdx    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdy    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(cntlat    , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(secondlat , 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
     endif

     call MPI_Bcast(buffer, isize, MPI_PACKED, 0, MPI_COMM_WORLD, ier)

     if (myrank /= 0) then
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
        call MPI_Unpack(buffer, isize, bytes, haso3     , 1, MPI_LOGICAL, MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xswlat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xswlon    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xnelat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xnelon    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, gdatdx    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, gdatdy    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, cntlat    , 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, secondlat , 1, MPI_REAL   , MPI_COMM_WORLD, ier)

        if (allocated(levpr)) deallocate(levpr)
        allocate(levpr(nprz))

        if (inproj == 2) then
           if (allocated(glat)) deallocate(glat)
           allocate(glat(npry))
        endif
     endif

     call MPI_Bcast(levpr, nprz, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (inproj == 2) then
        call MPI_Bcast(glat, npry, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
     endif

     deallocate(buffer)

  endif
#endif

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

  ! Convert input pressure levels to Pascals

  if (allocated(pnpr)) deallocate(pnpr)
  allocate(pnpr(nprz))

  if (levpr(1) < 2000) then

     write(io6,*)
     write(io6,*) "Assuming input pressure levels are in units of mb"
     write(io6,*)
     do k = 1,nprz
        pnpr(k) = real( levpr(k) ) * 100.
     enddo

  else

     write(io6,*)
     write(io6,*) "Assuming input pressure levels are in units of Pa"
     write(io6,*)
     do k = 1,nprz
        pnpr(k) = real( levpr(k) )
     enddo

  endif

end subroutine read_press_header



subroutine pressure_stage()

  use isan_coms,   only: npd, nprx, npry, nprz, gdatdy, gdatdx, pnpr, &
                         pcol_p, pcol_z, ipoffset, kzonoff, lzon_bot, &
                         inproj, plat, glat, xswlat, o_press, o_theta, o_rho, &
                         o_rrw, o_uzonal, o_umerid, o_ozone, haso3, o3name
  use hdf5_utils,  only: shdf5_info, shdf5_irec, shdf5_close
  use mem_para,    only: myrank
  use misc_coms,   only: rinit, io6, iparallel, i_o3
  use mem_ijtabs,  only: jtab_w, jtw_init
  use mem_zonavg,  only: zonz, zont, zonr, zonu, zono, zonp_vect
  use consts_coms, only: p00, rocp, eps_vap, cvocp, p00kord, eps_vapi
  use therm_lib,   only: eslf
  use mem_grid,    only: mza, mwa, zt
  use prfill_mod,  only: prfill, prfill3

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real, parameter :: mwair  = 28.9628             ! molecular weight of air
  real, parameter :: mwo3   = 48.0                ! molecular weight of ozone
  real, parameter :: cnvto3 = mwair / mwo3 * 1.e6 ! ozone mixing ratio to ppmV

  real :: yoffset, xoffset, gdyi, gdxi
  real :: dprat1, dprat2, dg, vapor_press

  integer :: k, j, i, iw
  integer :: levp, ilat, nprz_rh, nbot_o3
  integer :: ndims, idims(3)
  logical :: isrh

  real :: pcol_thet  (npd)
  real :: pcol_exneri(npd)
  real :: pvect      (npd)

  character(10) :: varname

  ! Array to store input lat/lon/pressure data
  real, allocatable :: as3(:,:,:), at3(:,:,:)

  ! Define expanded array with 2 added rows/colums at N, S, E, and W
  ! boundaries so that overlapping quadratic interpolation (in subroutine
  ! gdtost) has sufficient points to work from.  Input data may be offset by
  ! 1/2 grid cell from nodal latitudes and longitudes (e.g., (-180.,-90.)),
  ! in which case all added points would be required.
  real, allocatable :: p3d(:,:,:)

  ! Arrays to store data on analysis pressure levels interpolated to each
  ! model horizontal (iw) location
  real, allocatable :: field(:,:)

#ifdef OLAM_MPI
  integer :: ier, ireq
#endif

  yoffset = 3.0
  xoffset = 1.0 + real(ipoffset)

  gdyi = 1.0 / gdatdy
  gdxi = 1.0 / gdatdx

  ! Fill column array of pressure level values in MKS

  do k = 1, nprz
     pcol_p(k+2) = pnpr(k)
  enddo
  pcol_p(2) = 110000.   ! Phony underground level
  pcol_p(1) = 120000.   ! Phony underground level

  ! If necessary, fill upper pressure levels from ZONAVG data

  k = nprz + 2                ! highest ppd level filled so far
  kzonoff = k + 1 - lzon_bot  ! k offset for copying pzon levels to ppd

  do levp = lzon_bot, 22
     k = levp + kzonoff
     pcol_p(k) = zonp_vect(levp)
  enddo

  ! Inverse exner function (defined without cp factor)

  do k = 1, npd
     pcol_exneri(k) = ( p00 / pcol_p(k) )**rocp
  enddo

  ! Allocate space for 3d arrays

  if (myrank == 0) then
     allocate( as3(nprx,npry,nprz) )
  endif
  allocate(p3d(nprx+4,npry+4,nprz))

  allocate(field (npd,mwa)) ; field  = rinit

! For INPROJ = 2, where latitudes are specified, fill expanded latitude array

  if (inproj == 2) then

     if (allocated(plat)) deallocate(plat)
     allocate(plat(npry+4))

     do ilat = 1,npry
        plat(ilat+2) = glat(ilat)
     enddo

     plat(2) = plat(3) - (plat(4) - plat(3))
     plat(1) = plat(2) - (plat(3) - plat(2))

     plat(npry+3) = plat(npry+2) + (plat(npry+2) - plat(npry+1))
     plat(npry+4) = plat(npry+3) + (plat(npry+3) - plat(npry+2))

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read geopotential height
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) then

     ndims = 0
     idims = 0

     varname = 'GEO'
     call shdf5_info(varname, ndims, idims)

     if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
        write(*,*) "Geopotential height (GEO) not found in analysis file."
        stop
     endif

     write(io6,*) "Reading geopotential height " // trim(varname)
     call shdf5_irec(ndims, idims,varname,rvar3 = as3)
     call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                     ireq, ier )
     if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
  endif
#endif

  ! Horizontally interpolate geopotential of each pressure level
  ! to the model column lat/lon locations

  call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, pcol_z )

  ! Extra levels above analysis (from zonavg arrays)

  if (npd > nprz+2) then
     call hinterp_zonavg(pcol_z, zonz)
  endif

  ! Special for extrapolated levels at bottom

  dprat2 = 10000. / ((pcol_p(3) - pcol_p(4)) * 1.05)
  dprat1 = 10000. / ((pcol_p(3) - pcol_p(4)) * 1.13)

  ! Set phony underground height field and vertically interpolate pressure
  ! to OLAM vertical coordinate levels

  !$omp parallel do private(iw,dg)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
     dg           = pcol_z(4,iw) - pcol_z(3,iw)
     pcol_z(2,iw) = pcol_z(3,iw) - dg * dprat2
     pcol_z(1,iw) = pcol_z(2,iw) - dg * dprat1

     call hintrp_cc(npd, pcol_p, pcol_z(:,iw), mza, o_press(:,iw), zt)
  enddo
  !$omp end parallel do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) then
     ndims = 0
     idims = 0

     varname = 'TEMP'
     call shdf5_info(varname, ndims, idims)

     if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
        write(*,*) "Air temperature (TEMP) not found in analysis file."
        stop
     endif

     write(io6,*) "Reading air temperature " // trim(varname)
     call shdf5_irec(ndims, idims,varname,rvar3 = as3)

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
#endif
     call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                      ireq, ier )
     if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
  endif
#endif

  ! Horizontally interpolate geopotential of each pressure level
  ! to the model column lat/lon locations

  call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )

  ! Extra levels above analysis (from zonavg arrays)

  if (npd > nprz+2) then
     call hinterp_zonavg(field, zont)
  endif

  ! Set phony underground height field and vertically interpolate pressure
  ! levels to OLAM pressure field

  !$omp parallel private(pcol_thet)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     do k = 3, npd
        pcol_thet(k) = field(k,iw) * pcol_exneri(k)
     enddo

     pcol_thet(2) = pcol_thet(3)
     pcol_thet(1) = pcol_thet(2)

     call hintrp_cc(npd, pcol_thet, pcol_z(:,iw), mza, o_theta(:,iw), zt)
  enddo
  !$omp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read water vapor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) then

     ndims = 0
     idims = 0

     ! First check if specific humidity is in analysis file

     varname = 'SHV'
     isrh    = .false.
     call shdf5_info(varname, ndims, idims)

     if (ndims /= 3) then

        ! If we don't have specific humidity, check for relative humidity.
        ! It may be called RELHUM or RH

        varname = 'RH'
        isrh   = .true.
        call shdf5_info(varname, ndims, idims)

        if (ndims /= 3) then

           varname = 'RELHUM'
           isrh   = .true.
           call shdf5_info(varname, ndims, idims)

           if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
              write(*,*) "Water vapor not found in analysis file."
              stop
           endif

        endif
     endif

     ! If we need temperature to convert RH to mixing ratio, store it in
     ! p3d. Note: as3 is still the air temperature read in from the reanalysis

     if (isrh) then
        allocate( at3(nprx,npry,nprz) )
        at3 = as3
     endif

     write(io6,*) "Reading water vapor " // trim(varname)
     call shdf5_irec(ndims, idims, varname, rvar3=as3)

     if (isrh) then

        ! RH should be stored as a ratio (0-1). If there are larger values,
        ! assume RH is in percent and convert to decimal RH. This is probably
        ! not necessary with recent versions of grib2olam but we will keep
        ! this check anyway. TODO: Read and store units from grib2olam.
        if ( any(as3(:,:,1) > 2.0) ) then
           write(io6, *) '    Converting relative humidity (%) to ratio (0-1)'
           as3 = 0.01 * as3
        endif

        ! Convert relative humidity to mixing ratio

        write(io6,*) '    Converting relative humidity to mixing ratio'

        !$omp parallel do collapse(2) private(k,j,i,vapor_press)
        do k = 1, nprz
           do j = 1, npry
              do i = 1, nprx

                 ! Compute ambient vapor pressure based on R.H.
                 ! and saturation vapor pressure (eslf) (at3 temporarily
                 ! stores air temperature).

                 vapor_press = as3(i,j,k) * eslf(at3(i,j,k)-273.15)

                 ! Do not allow vapor pressure to exceed ambient pressure

                 vapor_press = min( 0.9*pnpr(k), vapor_press )

                 ! Compute mixing ratio from vapor press and ambient press

                 as3(i,j,k) = eps_vap * vapor_press &
                            / ( pnpr(k) - vapor_press )
              enddo
           enddo
        enddo
        !$omp end parallel do

        deallocate(at3)

     else

        ! If any specific humidities at lowest level are greater than 1,
        ! assume it is g/kg and convert to kg/kg. This is probably not
        ! necessary with recent versions of grib2olam but we will keep
        ! this check anyway. TODO: Read and store units from grib2olam.
        if ( any(as3(:,:,1) > 1.0) ) then
           write(io6, *) ' '
           write(io6, *) 'Converting g/kg specific humidity to kg/kg'
           as3 = 0.001 * as3
        endif

        ! Convert specific humidity to mixing ratio

        write(io6, *) ' '
        write(io6, *) 'Converting specific humidity to mixing ratio'

        !$omp parallel do
        do k = 1, nprz
           as3(:,:,k) = as3(:,:,k) / (1.0 - as3(:,:,k))
        enddo
        !$omp end parallel do

     endif

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
#endif
     call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                      ireq, ier )
     if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
  endif
#endif

! Special for RH:
! Reanalysis typically only reports RH up to 100mb, but still reports the other
! fields up to higher levels. Check for the highest level that reports RH:

  nprz_rh = 0

  do k = nprz,1,-1
     if (all(p3d(:,:,k) > -998.)) then
        nprz_rh = k
        exit
     endif
  enddo

  ! Horizontally interpolate humidity at each pressure level
  ! to the model column lat/lon locations

  call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )

  ! Extra levels above analysis (from zonavg arrays)

  if (npd > nprz+2) then
     call hinterp_zonavg(field, zonr)
  endif

  !$omp parallel private(pvect)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     ! Copy humidity pressure levels to temprorary vector
     pvect(3:npd) = field(3:npd,iw)

     ! Phony underground levels
     pvect(1:2) = pvect(3)

     ! Set any missing humidity levels from Mclatchy soundings
     ! (only for reanalyses that don't report humdidity at all levels)
     if (nprz_rh < nprz) then
        call fill_wvap_levs_mclat(iw, nprz_rh+3, nprz+2, pvect)
     endif

     ! Vertically interpolate humidity from pressure levels to model levels
     call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_rrw(:,iw), zt)

     do k = 1, mza
        o_rrw(k,iw) = max( 1.e-8, o_rrw(k,iw) )
     enddo

  enddo
  !$iomp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read zonal wind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) then

     ndims = 0
     idims = 0

     ! East-wst velocity may be called UP, UE, or U

     varname = 'U'
     call shdf5_info(varname, ndims, idims)

     if (ndims /= 3) then

        varname = 'UP'
        call shdf5_info(varname, ndims, idims)

        if (ndims /= 3) then

           varname = 'UE'
           call shdf5_info(varname, ndims, idims)

           if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
              write(*,*) "Zonal wind (U) not found in analysis file."
              stop
           endif

        endif
     endif

     write(io6,*) "Reading zonal wind " // trim(varname)
     call shdf5_irec(ndims, idims, varname, rvar3 = as3)

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
#endif
     call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                     ireq, ier )
     if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
  endif
#endif

  ! Horizontally interpolate zonal wind at each pressure level
  ! to the model column lat/lon locations

  call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )

  ! Extra levels above analysis (from zonavg arrays)

  if (npd > nprz+2) then
     call hinterp_zonavg(field, zonu)
  endif

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
     call hintrp_cc(npd, field(:,iw), pcol_z(:,iw), mza, o_uzonal(:,iw), zt)
  enddo
  !$omp end parallel do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read meridional wind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myrank == 0) then

     ndims = 0
     idims = 0

     ! East-wst velocity may be called VP, VE, or V

     varname = 'V'
     call shdf5_info(varname, ndims, idims)

     if (ndims /= 3) then

        varname = 'VP'
        call shdf5_info(varname, ndims, idims)

        if (ndims /= 3) then

           varname = 'VE'
           call shdf5_info(varname, ndims, idims)

           if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
              write(*,*) "Meridional wind (V) not found in analysis file."
              stop
           endif

        endif
     endif

     write(io6,*) "Reading meridional wind " // trim(varname)
     call shdf5_irec(ndims, idims, varname, rvar3 = as3)

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
#endif
     call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                      ireq, ier )
     if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
  endif
#endif

  ! Horizontally interpolate meridional wind at each pressure level
  ! to the model column lat/lon locations

  call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     if (npd > nprz+2) then
        field(nprz+3:npd,iw) = 0.0
     endif

     call hintrp_cc(npd, field(:,iw), pcol_z(:,iw), mza, o_umerid(:,iw), zt)
  enddo
  !$omp end parallel do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read ozone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (i_o3 > 0) then

     nbot_o3 = nprz + 1

     if (haso3) then

        if (myrank == 0) then

           ndims = 0
           idims = 0

           call shdf5_info(o3name, ndims, idims)

           if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
              write(*,*) "Ozone mixing ration not found in analysis file."
              stop
           endif

           write(io6,*) "Reading ozone mixing ratio " // trim(varname)
           call shdf5_irec(ndims, idims, varname, rvar3 = as3)

#ifdef OLAM_MPI
           if (iparallel == 1) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
#endif
           call prfill3(nprx,npry,nprz,as3,p3d,gdatdy,xswlat,ipoffset,inproj)

        endif

#ifdef OLAM_MPI
        if (iparallel == 1) then
           call MPI_Ibcast( p3d, (nprx+4)*(npry+4)*nprz, MPI_REAL, 0, MPI_COMM_WORLD, &
                            ireq, ier )
           if (myrank > 0) call MPI_Wait( ireq, MPI_STATUS_IGNORE, ier )
        endif
#endif

        ! Special for OZONE:
        ! Some analysis only reports ozone ABOVE 100 mb. Check for the lowest
        ! level at which ozone is reported in the analysis file

        do k = 1, nprz
           if (all(p3d(:,:,k) > -998.)) then
              nbot_o3 = k
              exit
           endif
        enddo

        ! Horizontally interpolate ozone at each pressure level
        ! to the model column lat/lon locations

        if (nbot_o3 < nprz+1) then
           call hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )
        endif

     endif

     ! Extra levels above analysis (from zonavg arrays)

     if (npd > nprz+2) then
        call hinterp_zonavg(field, zono)
     endif

     !$omp parallel private(pvect)
     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

        ! Copy ozone pressure levels to temprorary vector
        pvect(nbot_o3+2:npd) = field(nbot_o3:npd,iw)

        ! Set any missing ozone levels from Mclatchy soundings
        ! (only for reanalyses that don't report ozone at all levels)
        if ( nbot_o3 > 1 ) then
           call fill_ozone_levs_mclat(iw, 3, nbot_o3+1, pvect)
        endif

        ! Phony underground levels
        pvect(1:2) = pvect(3)

        ! Vertically interpolate ozone from pressure levels to model levels
        call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_ozone(:,iw), zt)

        ! Convert to ppmV
        do k = 1, mza
           o_ozone(k,iw) = max( 1.e-30, cnvto3 * o_ozone(k,iw) )
        enddo

     enddo
     !$omp end do
     !$omp end parallel

  endif

#ifdef OLAM_MPI
  if (iparallel == 1 .and. myrank == 0) then
     call MPI_Wait( ireq , MPI_STATUS_IGNORE, ier )
  endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Deallocate arrays and close HDF5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (allocated(as3)) deallocate(as3)
  if (allocated(p3d)) deallocate(p3d)

  if (myrank == 0) then
     call shdf5_close()
  endif

end subroutine pressure_stage



subroutine hinterp_plevs( gdxi, gdyi, xoffset, yoffset, p3d, field )

  use mem_grid,   only: mwa, glatw, glonw
  use mem_ijtabs, only: jtab_w, jtw_init
  use isan_coms

  implicit none

  ! This routine horizontally interpolates the input lon/lat/press
  ! analysis (p3d array) to each model column at the standard pressure
  ! levels.

  real, intent(in ) :: gdxi, gdyi, xoffset, yoffset
  real, intent(in ) :: p3d(nprx+4,npry+4,nprz)
  real, intent(out) :: field(npd,mwa)

  integer :: j, iw, ilat, k
  real    :: grx, gry

  !$omp parallel do private(iw,grx,gry,ilat,k)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     if (inproj == 1) then

        grx = (glonw(iw) - xswlon) * gdxi + xoffset
        gry = (glatw(iw) - xswlat) * gdyi + yoffset

     elseif (inproj == 2) then

        ! estimate latitude index assuming uniform spacing of plat

        ilat = 2 + npry * int((glatw(iw) - plat(2)) / 180.)

        ! find correct latitude index

        if (plat(ilat) > glatw(iw)) then
           do while(plat(ilat) > glatw(iw))
              ilat = ilat - 1
           enddo
        elseif (plat(ilat+1) < glatw(iw)) then
           do while(plat(ilat+1) < glatw(iw))
              ilat = ilat + 1
           enddo
        endif

        grx = (glonw(iw) - xswlon) * gdxi + xoffset
        gry = (glatw(iw) - plat(ilat)) / (plat(ilat+1) - plat(ilat)) + real(ilat)

     endif

     ! Horizontally interpolate gridded pressure-level data to column
     ! at location of current W point

     do k = 1, nprz
        call gdtost(p3d(:,:,k), nprx+4, npry+4, grx, gry, field(k+2,iw))
     enddo

     ! Temporarily set underground levels to value at lowest analysis level
     field(1:2,iw) = field(3,iw)

     ! Temporarily set added levels above analysis to value at highest analysis level
     if (npd > nprz+2) then
        field(nprz+3:npd,iw) = field(nprz+2,iw)
     endif

  enddo
  !$omp end parallel do

end subroutine hinterp_plevs




subroutine hinterp_zonavg(field, zona)

  use mem_grid,   only: mwa, glatw
  use mem_ijtabs, only: jtab_w, jtw_init
  use mem_zonavg, only: nlata, nplev
  use isan_coms

  implicit none

  real, intent(in   ) :: zona(nlata,nplev)
  real, intent(inout) :: field(npd,mwa)

  integer :: j, iw, ilat, k, levp
  real    :: rlat, wt2

  if (npd <= nprz+2) return

  !$omp parallel do private(iw,rlat,ilat,wt2,levp,k)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     rlat = .4 * (glatw(iw) + 93.75)
     ilat = int(rlat)
     wt2  = rlat - float(ilat)

     do levp = lzon_bot, 22
        k = levp + kzonoff
        field(k,iw) = (1. - wt2) * zona(ilat,levp) + wt2 * zona(ilat+1,levp)
     enddo

  enddo
  !$omp end parallel do

end subroutine hinterp_zonavg




subroutine fill_wvap_levs_mclat(iw, nbot, ntop, qvect)

  use mem_mclat, only: sslat, mclat, ypp_mclat, mclat_spline
  use mem_grid,  only: glatw
  use isan_coms, only: npd, pcol_p

  implicit none

  integer, intent(in)  :: iw, nbot, ntop
  real,    intent(out) :: qvect(npd)

  real    :: mcol(33,6)
  integer :: lv

  ! Fills pressure levels from nbot to ntop with water vapor values
  ! from the Mclatchy soundings

  if (nbot > ntop) return

  ! This assumes that mclat_spline has already been called
  ! with the current date

  call spline2_vec(13, 33*6, sslat, mclat, ypp_mclat, glatw(iw), mcol)

  ! Compute water vapor mixing ratios by dividing by dry density

  do lv = 1,33
     mcol(lv,4) = mcol(lv,4) / (mcol(lv,6) - mcol(lv,4))
  enddo

  ! Set any missing water vapor data

  lv = ntop - nbot + 1

  call pintrp_ee(33, mcol(1,5), mcol(1,2), lv, qvect(nbot:ntop), pcol_p(nbot:ntop))

end subroutine fill_wvap_levs_mclat




subroutine fill_ozone_levs_mclat(iw, nbot, ntop, ovect)

  use mem_mclat, only: sslat, mclat, ypp_mclat, mclat_spline
  use mem_grid,  only: glatw
  use isan_coms, only: npd, pcol_p

  implicit none

  integer, intent(in)  :: iw, nbot, ntop
  real,    intent(out) :: ovect(npd)

  real    :: mcol(33,6)
  integer :: lv

  ! Fills pressure levels from nbot to ntop with ozone values
  ! from the Mclatchy soundings

  if (nbot > ntop) return

  ! This assumes that mclat_spline has already been called
  ! with the current date

  call spline2_vec(13, 33*6, sslat, mclat, ypp_mclat, glatw(iw), mcol)

  ! Compute ozone mixing ratios by dividing by dry density

  do lv = 1,33
     mcol(lv,5) = mcol(lv,5) / (mcol(lv,6) - mcol(lv,4))
  enddo

  ! Vertically interpolate Mclatchy water vapor BY PRESSURE to analysis levels

  lv = ntop - nbot + 1

  call pintrp_ee(33, mcol(1,4), mcol(1,2), lv, ovect(nbot:ntop), pcol_p(nbot:ntop))

end subroutine fill_ozone_levs_mclat
