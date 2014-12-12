subroutine read_seaice_analysis(iaction)

  use mem_sea,    only: sea, itab_ws
  use sea_coms,   only: mws, iseaicefile, iseaiceflg
  use misc_coms,  only: io6, s1900_sim, s1900_init, iparallel
  use max_dims,   only: pathlen
  use isan_coms,  only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, &
                        inproj, xswlat, xswlon, gdatdx, gdatdy, ipoffset
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use mem_para,   only: myrank

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: iaction

  integer            :: nf
  character(pathlen) :: fname
  character(16)      :: ext
  integer            :: ndims, idims(3)
  real               :: grx, gry
  integer            :: nio, njo, iws
  logical            :: exists, has_seaice
  integer            :: bytes, nbytes_int, nbytes_real, isize, ier
  
  real,    allocatable :: ice(:,:)  ! sea ice concentration [0 - 1]
  real,    allocatable :: a2d(:,:)
  integer, allocatable :: buffer(:)

! Nothing to do here if iseaiceflg is not 2

  if (iseaiceflg /= 2) return

  has_seaice  = .false.

  if (iaction == 0) then

     ! Loop over analysis files and search for the one that corresponds
     ! to the model start time

     iseaicefile = 0
     do nf = 1, nfgfiles
        if (s1900_fg(nf) <= s1900_sim) then
           iseaicefile = nf
        endif
     enddo
  
     if (iseaicefile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find analysis file for sst/seaice'
        write(io6,*) 'Stopping model.'
        stop 'stop: no current analysis file for sst/seaice file'
     endif

  elseif (iaction == 1) then

     ! Processing next seaice file (only called with iaction = 1 if iupdseaice = 1)

     iseaicefile = iseaicefile + 1

     if (iseaicefile > nfgfiles) then
        write(io6,*) ' '
        write(io6,*) 'No future analysis file is available for sst/seaice'
        write(io6,*) 'Stopping model.'
        stop 'stop: no future analysis file for sst/seaice file'
     endif

     sea%seaicep(:) = sea%seaicef(:)   

  endif

  fname = fnames_fg(iseaicefile)

! Make sure we have a HDF5 analysis file

  ext = trim( fname( index(fname,'.',back=.true.)+1:) )

  if ( .not. ( &
       ext == 'h5' .or. ext == 'hdf5' .or. ext == 'hdf' .or. &
       ext == 'H5' .or. ext == 'HDF5' .or. ext == 'HDF' )) then
     write(io6,*) "read_seaice: Analysis file is not hdf5 format"
     write(io6,*) "Cannot read sst/seaice data."
     stop
  endif

! Process selected analysis file

  inquire(file=fname, exist=exists)
  if (.not. exists) then
     write(io6,'(A)') " read_seaice: Error opening analysis file " // trim(fname)
     write(io6,'(A)') " Using default soil initialization instead."
     return
  endif

  write(io6,'(A)') ' read_seaice: opening ' // trim(fname)

#ifdef OLAM_MPI
  if (iparallel == 1) then

     call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, nbytes_int , ier)
     call MPI_Pack_size(1, MPI_REAL   , MPI_COMM_WORLD, nbytes_real, ier)

     bytes = 0
     isize = nbytes_int*3 + nbytes_real*4
     allocate( buffer( isize ) )
  endif
#endif

  if (myrank == 0) then

     call shdf5_open (fname, 'R')

     ndims    = 1
     idims(1) = 1
     idims(2) = 1
     idims(3) = 1

     call shdf5_irec(ndims, idims, 'nx'   , ivars=nprx)
     call shdf5_irec(ndims, idims, 'ny'   , ivars=npry)
     call shdf5_irec(ndims, idims, 'iproj', ivars=inproj)
     call shdf5_irec(ndims, idims, 'swlat', rvars=xswlat)
     call shdf5_irec(ndims, idims, 'swlon', rvars=xswlon)
     call shdf5_irec(ndims, idims, 'dx'   , rvars=gdatdx)
     call shdf5_irec(ndims, idims, 'dy'   , rvars=gdatdy)

#ifdef OLAM_MPI
     if (iparallel == 1) then
        call MPI_Pack(nprx  , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(npry  , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(inproj, 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlat, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlon, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdx, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdy, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
     endif
#endif

  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then

     call MPI_Bcast(buffer, isize, MPI_PACKED, 0, MPI_COMM_WORLD, ier)

     if (myrank /= 0) then
        call MPI_Unpack(buffer, isize, bytes, nprx  , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, npry  , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, inproj, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xswlat, 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, xswlon, 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, gdatdx, 1, MPI_REAL   , MPI_COMM_WORLD, ier)
        call MPI_Unpack(buffer, isize, bytes, gdatdy, 1, MPI_REAL   , MPI_COMM_WORLD, ier)
     endif

     deallocate(buffer)

  endif
#endif

  ! Check data domain size and location

  if (inproj /= 1) then
     write(io6,*) 'You must input a lat-lon grid for sst/seaice initialization.'
     write(io6,*) 'Stopping run.'
     if (myrank == 0) call shdf5_close()
     stop
  endif

  ! We make the requirement that a full global domain of data be
  ! read in.  Check this here.  Following the convention for the NCEP/DOE
  ! Reanalysis2 data, assume that data exists at both latitudinal boundaries
  ! (-90. and 90. degrees) but that the longitudinal boundary is not repeated.
  ! If either is not the case, this check will stop execution. 

  if (abs(nprx      * gdatdx - 360.) > .1 .or. &
      abs ((npry-1) * gdatdy - 180.) > .1) then
     write(io6,*) 'Gridded seaice data does not have global coverage.'
     write(io6,*) 'Stopping model run.'
     if (myrank == 0) call shdf5_close()
     stop
  endif

  ! Compute longitudinal offset index, which is the index in the expanded
  ! arrays where the first input data point (at xswlon) is located.

  ipoffset = int((xswlon + 180.) / gdatdx) + 2
  nio = nprx + 4
  njo = npry + 4

  allocate(ice(nio,njo))

  ! Check if seaice is in the analysis file, and read it

  if (myrank == 0) then
     call shdf5_info('ICEC', ndims, idims)

     if (ndims > 0) then
        allocate(a2d(nprx,npry))

        call shdf5_irec(ndims, idims, 'ICEC', rvara = a2d)
        call prfill(nprx, npry, a2d, ice)

        has_seaice = .true.
        deallocate(a2d)
     endif

     call shdf5_close()
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(has_seaice, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  if (.not. has_seaice) then
     write(io6,*) "read_seaice: Analysis file does not contain seaice."
     write(io6,*) "Stopping run."
     stop
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(ice, nio*njo, MPI_REAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  ! Fill seaice array

  do iws = 2, mws

! fractional x/y indices in pressure data arrays at current iw point location 

     gry = (sea%glatw(iws) - xswlat) / gdatdy + 3.
     grx = (sea%glonw(iws) - xswlon) / gdatdx + 1. + real(ipoffset) 

     call gdtost(ice, nprx+4, npry+4, grx, gry, sea%seaicef(iws))

     sea%seaicef(iws) = max( min( sea%seaicef(iws), 1.0 ), 0.0 )

  enddo

  deallocate( ice )

  if (iaction == 0) then
     sea%seaicep(:) = sea%seaicef(:)
  endif

end subroutine read_seaice_analysis
