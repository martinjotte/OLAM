module geia_emis

  real, allocatable ::  cl_emis(:)
  real, allocatable :: hcl_emis(:)

contains

subroutine geia_init()

  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use mem_para,   only: myrank
  use isan_coms,  only: gdatdx, gdatdy, xswlat, xswlon, ipoffset
  use misc_coms,  only: io6, iparallel
  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
  use mem_grid,   only: mwa, arw0, glatw, glonw
  use oname_coms, only: nl

#ifdef OLAM_MPI
  use mpi
#endif

  integer, parameter :: nlon = 360
  integer, parameter :: nlat = 180

  integer, parameter :: nio = nlon + 4
  integer, parameter :: njo = nlat + 4

  real, parameter :: year2sec = 365. * 24. * 3600.
  real, parameter :: cl_mwhgt = 35.5

  logical :: exists
  integer :: j, iw
  integer :: ndims, idims(2)

  real :: buffer(nlon, nlat, 2)

  real ::  cl(nio, njo)
  real :: hcl(nio, njo)

  allocate( cl_emis(mwa))
  allocate(hcl_emis(mwa))

   cl = 0.0
  hcl = 0.0

   cl_emis = 0.0
  hcl_emis = 0.0

  inquire(file=nl%geia_emis_file, exist=exists)
  if (.not. exists) then
     write(io6,*) "Error opening emissions file " // trim(nl%geia_emis_file)
     write(io6,*) "CL and HCL emissions will be set to 0"
     return
  endif

  if (myrank == 0) then
     call shdf5_open(nl%geia_emis_file, 'R')

     ndims    = 2
     idims(1) = nlon
     idims(2) = nlat

     call shdf5_irec(ndims, idims,  'CL', rvara=buffer(:,:,1))
     call shdf5_irec(ndims, idims, 'HCL', rvara=buffer(:,:,2))

     call shdf5_close()
  endif
        
#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(buffer, nlon*nlat*2, MPI_REAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  gdatdx = 1.0
  gdatdy = 1.0
  
  xswlat =  -89.5
  xswlon = -179.5
  ipoffset = int((xswlon + 180.) / gdatdx) + 2

  call prfill(nlon, nlat, buffer(:,:,1),  cl)
  call prfill(nlon, nlat, buffer(:,:,2), hcl)

  ! Fill emissions arrays by interpolation

  do j = 1, jtab_w(jtw_prog)%jend(1)
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

end subroutine geia_init

end module geia_emis
