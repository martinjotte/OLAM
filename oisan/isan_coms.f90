Module isan_coms

  use max_dims,    only: maxisdirs, pathlen
  use consts_coms, only: r8

  implicit none

  integer :: iyear, imonth, idate, ihour

  integer :: nfgfiles = -1
  integer :: ifgfile  =  0

  character(pathlen) :: iapr(maxisdirs)

  character(pathlen), allocatable :: fnames_fg  (:)
  character(14),      allocatable :: ctotdate_fg(:)
  real(r8),           allocatable :: s1900_fg   (:)

! Pressure header variables:
  integer :: isversion, iyy, imm, idd, ihh, itinc, inproj, ivertcoord
  integer :: nprx, npry, nprz
  real    :: xswlon, xswlat, gdatdx, gdatdy

  logical :: irev_ns = .false.

  real, allocatable :: pnpr (:)
  real, allocatable :: plats(:)

  ! input variables interpolated to height at each OLAM column

  real, allocatable :: o_press (:,:)
  real, allocatable :: o_rho   (:,:)
  real, allocatable :: o_theta (:,:)
  real, allocatable :: o_rrw   (:,:)
  real, allocatable :: o_uzonal(:,:)
  real, allocatable :: o_umerid(:,:)
  real, allocatable :: o_ozone (:,:)

  real              :: pbc
  real, allocatable :: z_pbc(:)

  private :: maxisdirs, pathlen, r8


Contains


subroutine read_analysis_header(noplevs, nosoil)

  use misc_coms,  only: io6
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info

! Fortran 2018
! import, only: nprx, npry, nprz, pnpr, plats, irev_ns, inproj, &
!               xswlat, xswlon, gdatdx, gdatdy

  implicit none


  logical, optional, intent(in) :: noplevs
  logical, optional, intent(in) :: nosoil

  logical :: doplevs, dosoil
  integer :: igloberr, ipry, k
  integer :: ndims, idims(3)

  integer, allocatable :: levpr(:)

  doplevs = .true.
  if (present(noplevs)) then
     doplevs = (.not. noplevs)
  endif

  dosoil = .true.
  if (present(nosoil)) then
     dosoil = (.not. nosoil)
  endif

  nprz = 0
  if (allocated(pnpr)) deallocate(pnpr)

  if (allocated(plats)) deallocate(plats)

  ndims    = 1
  idims(1) = 1

  ! Only for netCDF files with reversed latitudes N->S
  irev_ns = .false.

! call shdf5_irec(ndims, idims, 'version',ivars=isversion)
! call shdf5_irec(ndims, idims, 'year'   ,ivars=iyy)
! call shdf5_irec(ndims, idims, 'month'  ,ivars=imm)
! call shdf5_irec(ndims, idims, 'day'    ,ivars=idd)
! call shdf5_irec(ndims, idims, 'hour'   ,ivars=ihh)
! call shdf5_irec(ndims, idims, 'ftime'  ,ivars=itinc)
  call shdf5_irec(ndims, idims, 'nx'     ,ivars=nprx)
  call shdf5_irec(ndims, idims, 'ny'     ,ivars=npry)
  call shdf5_irec(ndims, idims, 'iproj'  ,ivars=inproj)
! call shdf5_irec(ndims, idims, 'vcoord' ,ivars=ivertcoord)
  call shdf5_irec(ndims, idims, 'swlat'  ,rvars=xswlat)
  call shdf5_irec(ndims, idims, 'swlon'  ,rvars=xswlon)
! call shdf5_irec(ndims, idims, 'nelat'  ,rvars=xnelat)
! call shdf5_irec(ndims, idims, 'nelon'  ,rvars=xnelon)
  call shdf5_irec(ndims, idims, 'dx'     ,rvars=gdatdx)
  call shdf5_irec(ndims, idims, 'dy'     ,rvars=gdatdy)
! call shdf5_irec(ndims, idims, 'reflat1',rvars=cntlat)
! call shdf5_irec(ndims, idims, 'reflat2',rvars=secondlat)

  if (.not. allocated(plats)) allocate( plats(npry) )
  if (inproj == 2) then
     ! allocate( plats(npry) )

     ndims    = 1
     idims(1) = npry
     call shdf5_irec(ndims, idims, 'glat', rvar1=plats)
  endif

!!  write(io6,*) 'nprz1 ',nprz,nprx,npry
!!
!!  ! Check for consistency between file parameters and namelist parameters
!!
!!  if (iyy /= iyear .or. imm /= imonth .or. idd /= idate .or. ihh /= ihour) then
!!
!!     write(io6,*) 'Pressure file dates not the same as namelist!'
!!     write(io6,*) 'Year :',iyy,iyear
!!     write(io6,*) 'Month:',imm,imonth
!!     write(io6,*) 'Day  :',idd,idate
!!     write(io6,*) 'Hour :',ihh,ihour
!!     stop 'pr_dates'
!!  endif

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

     if (plats(1) - (plats(2) - plats(1)) > -90.) igloberr = 1
     if (plats(npry) + (plats(npry) - plats(npry-1)) < 90.) igloberr = 1

     if (igloberr == 1) then
        write(io6,*) 'INPROJ = ',inproj
        write(io6,*) 'Gridded pressure level data must have global coverage'
        write(io6,*) 'nprx,npry = ',nprx,npry
        write(io6,*) 'gdatdx = ',gdatdx
        write(io6,*) 'xswlat,xswlon= ',xswlat,xswlon
        do ipry = 1,npry
           write(io6,*) 'ipry, glat = ',ipry,plats(ipry)
        enddo
        stop 'astp stop2 - non-global domain in input pressure data'
     endif

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

  if (doplevs) then

     ndims    = 1
     idims(1) = 1
     call shdf5_irec(ndims, idims, 'nlev', ivars=nprz)

     allocate(levpr(nprz))

     ndims    = 1
     idims(1) = nprz
     call shdf5_irec(ndims, idims, 'levels', ivar1=levpr)

     ! Convert input pressure levels to Pascals

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

  endif

end subroutine read_analysis_header


End module isan_coms
