module geia_emis

  real, allocatable ::  cl_emis(:)
  real, allocatable :: hcl_emis(:)

contains

subroutine geia_init()

  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use misc_coms,  only: io6
  use mem_ijtabs, only: jtab_w, jtw_prog
  use mem_grid,   only: mwa, arw0, glatw, glonw
  use oname_coms, only: nl
  use prfill_mod, only: prfill

  implicit none

  integer, parameter :: nlon = 360
  integer, parameter :: nlat = 180

  integer, parameter :: nio = nlon + 4
  integer, parameter :: njo = nlat + 4

  real, parameter :: year2sec = 365. * 24. * 3600.
  real, parameter :: cl_mwght = 35.5

  logical :: exists
  integer :: j, iw, ier
  integer :: ndims, idims(2)
  real    :: gdatdx, gdatdy, xswlat, xswlon
  integer :: ipoffset, inproj
  real    :: grx, gry

  real :: buffer(nlon,nlat)

  real ::  cl(nio, njo)
  real :: hcl(nio, njo)

  allocate( cl_emis(mwa))
  allocate(hcl_emis(mwa))

   cl_emis = 0.0
  hcl_emis = 0.0

  inquire(file=nl%geia_emis_file, exist=exists)
  if (.not. exists) then
     write(io6,*) "Error opening emissions file " // trim(nl%geia_emis_file)
     write(io6,*) "CL and HCL emissions will be set to 0"
     return
  endif

  call shdf5_open(nl%geia_emis_file, 'R', trypario=.true.)

   cl = 0.0
  hcl = 0.0

  ndims    = 2
  idims(1) = nlon
  idims(2) = nlat

  gdatdx = 1.0
  gdatdy = 1.0

  xswlat =  -89.5
  xswlon = -179.5

  ipoffset = int((xswlon + 180.) / gdatdx) + 2
  inproj = 1

  call shdf5_irec(ndims, idims,  'CL', rvar2=buffer)
  call prfill(nlon, nlat, buffer,  cl, gdatdy, xswlat, ipoffset, inproj)

  call shdf5_irec(ndims, idims, 'HCL', rvar2=buffer)
  call prfill(nlon, nlat, buffer, hcl, gdatdy, xswlat, ipoffset, inproj)

  call shdf5_close()

  ! Fill emissions arrays by interpolation

  !$omp parallel do private(j,iw,gry,grx)
  do j = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(j)

     ! fractional x/y indices in data arrays at current iw point location

     gry = (glatw(iw) - xswlat) / gdatdy + 3.
     grx = (glonw(iw) - xswlon) / gdatdx + 1. + real(ipoffset)

     call gdtost( cl, nio, njo, grx, gry,  cl_emis(iw))
     call gdtost(hcl, nio, njo, grx, gry, hcl_emis(iw))

     ! convert cl emissions to gm/sec emissions for particulates
     cl_emis(iw) = cl_emis(iw) * arw0(iw) / year2sec

     ! convert hcl to mol/sec emissions for gases
     hcl_emis(iw) = hcl_emis(iw) * arw0(iw) / year2sec / cl_mwght

  enddo
  !$omp end parallel do

end subroutine geia_init

end module geia_emis
