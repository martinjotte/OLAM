subroutine read_sst_analysis(iaction)

  use mem_sea,    only: sea, itab_ws
  use sea_coms,   only: mws, isstfile, isstflg
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
  real               :: xoffpix, yoffpix, xperdeg, yperdeg
  real               :: glat, glon, rio, rjo
  real               :: wio1, wio2, wjo1, wjo2
  integer            :: nio, njo, iws
  integer            :: io1, io2, jo1, jo2
  logical            :: exists, has_sst
  integer            :: bytes, nbytes_int, nbytes_real, isize, ier

  real,    allocatable :: sst(:,:)  ! sea surface temperature [K]
  real,    allocatable :: a2d(:,:)
  integer, allocatable :: buffer(:)

! Nothing to do here if isstflg is not 2

  if (isstflg /= 2) return

  has_sst  = .false.

  if (iaction == 0) then

     ! Loop over analysis files and search for the one that corresponds
     ! to the model start time

     isstfile = 0
     do nf = 1, nfgfiles
        if (s1900_fg(nf) <= s1900_sim) then
           isstfile = nf
        endif
     enddo

     if (isstfile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find analysis file for sst/seaice'
        write(io6,*) 'Stopping model.'
        stop 'stop: no current analysis file for sst/seaice file'
     endif

  elseif (iaction == 1) then

     ! Processing next sst file (only called with iaction = 1 if iupdsst = 1)

     isstfile = isstfile + 1

     if (isstfile > nfgfiles) then
        write(io6,*) ' '
        write(io6,*) 'No future analysis file is available for sst/seaice'
        write(io6,*) 'Stopping model.'
        stop 'stop: no future analysis file for sst/seaice file'
     endif

     sea%seatp(:) = sea%seatf(:)   

  endif

  fname = fnames_fg(isstfile)

! Make sure we have a HDF5 analysis file

  ext = trim( fname( index(fname,'.',back=.true.)+1:) )

  if ( .not. ( &
       ext == 'h5' .or. ext == 'hdf5' .or. ext == 'hdf' .or. &
       ext == 'H5' .or. ext == 'HDF5' .or. ext == 'HDF' )) then
     write(io6,*) "read_sst: Analysis file is not hdf5 format"
     write(io6,*) "Cannot read sst/seaice data."
     stop
  endif

! Process selected analysis file

  inquire(file=fname, exist=exists)
  if (.not. exists) then
     write(io6,'(A)') " read_sst: Error opening analysis file " // trim(fname)
     write(io6,'(A)') " Using default soil initialization instead."
     return
  endif

  write(io6,'(A)') ' read_sst: opening ' // trim(fname)

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
     write(io6,*) 'Gridded sst data does not have global coverage.'
     write(io6,*) 'Stoppig model run.'
     if (myrank == 0) call shdf5_close()
     stop
  endif

  ! Compute longitudinal offset index, which is the index in the expanded
  ! arrays where the first input data point (at xswlon) is located. We
  ! follow the GRIB standard and keep the SW corner at (0,-90)

  ipoffset = xswlon / gdatdx + 1
  nio = nprx + 3
  njo = npry + 2

  allocate(sst(nio,njo))

  ! Check if sst is in the analysis file, and read it

  if (myrank == 0) then
     call shdf5_info('SST', ndims, idims)

     if (ndims > 0) then
        allocate(a2d(nprx,npry))

        call shdf5_irec(ndims, idims, 'SST', rvara = a2d)
        call prfill(nprx, npry, a2d, sst)

        has_sst = .true.
        deallocate(a2d)
     endif

     call shdf5_close()
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(has_sst, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  if (.not. has_sst) then
     write(io6,*) "read_sst: Analysis file does not contain sst."
     write(io6,*) "Stopping run."
     stop
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(sst, nio*njo, MPI_REAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  ! GRIB data is unstaggered, so there should be no offsets, but we 
  ! offset the grid by 1 row and column so we take that into account here

  xoffpix = 1.0
  yoffpix = 1.0

  xperdeg = 1.0 / gdatdx
  yperdeg = 1.0 / gdatdy

  ! Fill sst array

  do iws = 2, mws

     glat = sea%glatw(iws)
     glon = sea%glonw(iws)

     ! Convert OLAM's -180,180 longitude range to GRIB's 0,360 

     if (glon < 0.0) glon = 360.0 + glon 
     glon = max(0.001,min(359.999,glon))

     ! Find the nearest analysis point to this land cell

     rio = 1. + (glon      ) * xperdeg + xoffpix
     rjo = 1. + (glat + 90.) * yperdeg + yoffpix

     io1 = int(rio)
     jo1 = int(rjo)
         
     wio2 = rio - real(io1)
     wjo2 = rjo - real(jo1)
           
     wio1 = 1. - wio2
     wjo1 = 1. - wjo2

     io2 = min(nio, io1 + 1)
     jo2 = min(njo, jo1 + 1)

     sea%seatf(iws) =  &
          + wio1 * (wjo1 * sst(io1,jo1) + wjo2 * sst(io1,jo2))  &
          + wio2 * (wjo1 * sst(io2,jo1) + wjo2 * sst(io2,jo2))
     
  enddo

  deallocate( sst )

  if (iaction == 0) then
     sea%seatp(:) = sea%seatf(:)
  endif

end subroutine read_sst_analysis
