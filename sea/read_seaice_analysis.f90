subroutine read_seaice_analysis(iaction)

  use mem_sea,    only: sea, itab_ws
  use sea_coms,   only: mws, iseaicefile, iseaiceflg
  use misc_coms,  only: io6, s1900_sim, s1900_init
  use max_dims,   only: pathlen
  use isan_coms,  only: nfgfiles, s1900_fg, fnames_fg
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close

  implicit none

  integer, intent(in) :: iaction

  integer            :: nf, ipoffset
  character(pathlen) :: fname
  character(16)      :: ext
  integer            :: ndims, idims(3)
  integer            :: nx, ny, inproj
  real               :: xswlat, xswlon, gdatdx, gdatdy
  real               :: xoffpix, yoffpix, xperdeg, yperdeg
  real               :: glat, glon, rio, rjo
  real               :: wio1, wio2, wjo1, wjo2
  integer            :: nio, njo, iws
  integer            :: io1, io2, jo1, jo2
  logical            :: exists, has_seaice

  real, allocatable  :: ice(:,:)  ! sea ice concentration [0 - 1]
  real, allocatable  :: a2d(:,:)

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
     write(io6,*) "read_soil: Error opening analysis file " // trim(fname)
     write(io6,*) "Using default soil initialization instead."
     return
  endif

  write(io6,*) 'read_seaice_analysis: opening ' // trim(fname)

  call shdf5_open (fname, 'R')

  ndims    = 1
  idims(1) = 1
  idims(2) = 1
  idims(3) = 1

  call shdf5_irec(ndims, idims, 'nx'   , ivars=nx)
  call shdf5_irec(ndims, idims, 'ny'   , ivars=ny)
  call shdf5_irec(ndims, idims, 'iproj', ivars=inproj)
  call shdf5_irec(ndims, idims, 'swlat', rvars=xswlat)
  call shdf5_irec(ndims, idims, 'swlon', rvars=xswlon)
  call shdf5_irec(ndims, idims, 'dx'   , rvars=gdatdx)
  call shdf5_irec(ndims, idims, 'dy'   , rvars=gdatdy)

  ! Check data domain size and location

  if (inproj /= 1) then
     write(io6,*) 'You must input a lat-lon grid for sst/seaice initialization.'
     write(io6,*) 'Stopping run.'
     call shdf5_close()
     stop
  endif

  ! We make the requirement that a full global domain of data be
  ! read in.  Check this here.  Following the convention for the NCEP/DOE
  ! Reanalysis2 data, assume that data exists at both latitudinal boundaries
  ! (-90. and 90. degrees) but that the longitudinal boundary is not repeated.
  ! If either is not the case, this check will stop execution. 

  if (abs(nx      * gdatdx - 360.) > .1 .or. &
      abs ((ny-1) * gdatdy - 180.) > .1) then
     write(io6,*) 'Gridded seaice data does not have global coverage.'
     write(io6,*) 'Stopping model run.'
     call shdf5_close()
     stop
  endif

  ! Compute longitudinal offset index, which is the index in the expanded
  ! arrays where the first input data point (at xswlon) is located. We
  ! follow the GRIB standard and keep the SW corner at (0,-90)

  ipoffset = xswlon / gdatdx + 1

  ! Check if seaice is in the analysis file, and read it

  call shdf5_info('ICEC', ndims, idims)

  if (ndims > 0) then

     nio = nx + 3
     njo = ny + 2

     allocate(a2d(nx,ny))
     allocate(ice(nio,njo))

     call shdf5_irec(ndims, idims, 'ICEC', rvara = a2d)
     call prfill(nx, ny, ipoffset, a2d, ice)

     has_seaice = .true.
     deallocate(a2d)
  endif

  call shdf5_close()

  if (.not. has_seaice) then
     write(io6,*) "read_seaice: Analysis file does not contain seaice."
     write(io6,*) "Stopping run."
     stop
  endif

  ! GRIB data is unstaggered, so there should be no offsets, but we 
  ! offset the grid by 1 row and column so we take that into account here

  xoffpix = 1.0
  yoffpix = 1.0

  xperdeg = 1.0 / gdatdx
  yperdeg = 1.0 / gdatdy

  ! Fill seaice array

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

     sea%seaicef(iws) =  &
          + wio1 * (wjo1 * ice(io1,jo1) + wjo2 * ice(io1,jo2))  &
          + wio2 * (wjo1 * ice(io2,jo1) + wjo2 * ice(io2,jo2))
     
     sea%seaicef(iws) = max( min( sea%seaicef(iws), 1.0 ), 0.0 )

  enddo

  deallocate( ice )

end subroutine read_seaice_analysis
