module geia_emis

  real, allocatable ::  cl_emis(:)
  real, allocatable :: hcl_emis(:)

contains

subroutine geia_init()

  use hdf5_utils,  only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use misc_coms,   only: io6
  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: mwa, arw0, glatw, glonw
  use oname_coms,  only: nl
  use analysis_lib,only: gdtost_ll

  implicit none

  integer, parameter :: nlon = 360
  integer, parameter :: nlat = 180

  real, parameter :: year2sec = 365. * 24. * 3600.
  real, parameter :: cl_mwght = 35.5

  logical :: exists
  integer :: j, iw, ier, inproj
  integer :: ndims, idims(2)
  real    :: dlon, dlat, swlon, swlat

  real, allocatable :: buffer(:,:)

  allocate( cl_emis(mwa))
  allocate(hcl_emis(mwa))

  !$omp parallel do
  do iw = 1, mwa
      cl_emis(iw) = 0.0
     hcl_emis(iw) = 0.0
  enddo
  !$omp end parallel do

  inquire(file=nl%geia_emis_file, exist=exists)
  if (.not. exists) then
     write(io6,*) "Error opening emissions file " // trim(nl%geia_emis_file)
     write(io6,*) "CL and HCL emissions will be set to 0"
     return
  endif

  call shdf5_open(nl%geia_emis_file, 'R', trypario=.true.)

  ndims    = 2
  idims(1) = nlon
  idims(2) = nlat

  dlon = 1.0
  dlat = 1.0

  swlon = -179.5
  swlat =  -89.5

  inproj = 1

  allocate(buffer(nlon,nlat))

  ! Read and interpolate CL emissions

  call shdf5_irec(ndims, idims, 'CL', rvar2=buffer)

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(j)

     call gdtost_ll(nlon, nlat, glonw(iw), glatw(iw), cl_emis(iw), &
                    swlon, swlat, dlon, dlat, inproj, r2d=buffer)

     ! convert cl emissions to gm/sec for particulates
     cl_emis(iw) = cl_emis(iw) * arw0(iw) / year2sec

  enddo
  !$omp end parallel do

  ! Read and interpolate HCL emissions

  call shdf5_irec(ndims, idims, 'HCL', rvar2=buffer)

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(j)

     call gdtost_ll(nlon, nlat, glonw(iw), glatw(iw), hcl_emis(iw), &
                    swlon, swlat, dlon, dlat, inproj, r2d=buffer)

     ! convert hcl emissions to mol/sec for gases
     hcl_emis(iw) = hcl_emis(iw) * arw0(iw) / year2sec / cl_mwght

  enddo
  !$omp end parallel do

  call shdf5_close()

end subroutine geia_init

end module geia_emis
