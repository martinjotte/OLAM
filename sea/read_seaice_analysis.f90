subroutine read_seaice_analysis(iaction)

  use mem_sea,    only: sea, msea, omsea
  use mem_sfcg,   only: sfcg
  use sea_coms,   only: iseaicefile, iseaiceflg
  use misc_coms,  only: io6, s1900_sim
  use max_dims,   only: pathlen
  use isan_coms,  only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, glat, &
                        inproj, xswlat, xswlon, gdatdx, gdatdy, ipoffset
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use prfill_mod, only: prfill

  implicit none

  integer, intent(in) :: iaction

  integer            :: nf
  character(pathlen) :: fname
  character(16)      :: ext
  integer            :: ndims, idims(3)
  real               :: grx, gry
  integer            :: nio, njo, isea, iwsfc
  logical            :: exists, has_seaice
  integer            :: igloberr, ilat, ipry

  real,    allocatable :: ice(:,:)  ! sea ice concentration [0 - 1]
  real,    allocatable :: a2d(:,:)
  real,    allocatable :: plat(:)

! Nothing to do here if iseaiceflg is not 2

  if (iseaiceflg /= 2) return

  has_seaice  = .false.

  if (iaction == 0) then

     ! Get the listing of analysis files if the inventory was
     ! never called
     if (nfgfiles < 0) call isan_file_inv()

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
     write(io6,'(A)') " Stopping model."
     stop 'stop: no future analysis file for sst/seaice file'
  endif

  write(io6,'(A)') ' read_seaice: opening ' // trim(fname)

  call shdf5_open (fname, 'R', trypario=.true.)

  ndims    = 1
  idims(1) = 1

  call shdf5_irec(ndims, idims, 'nx'   , ivars=nprx)
  call shdf5_irec(ndims, idims, 'ny'   , ivars=npry)
  call shdf5_irec(ndims, idims, 'iproj', ivars=inproj)
  call shdf5_irec(ndims, idims, 'swlat', rvars=xswlat)
  call shdf5_irec(ndims, idims, 'swlon', rvars=xswlon)
  call shdf5_irec(ndims, idims, 'dx'   , rvars=gdatdx)
  call shdf5_irec(ndims, idims, 'dy'   , rvars=gdatdy)

  if (inproj == 2) then
     if (allocated(glat)) deallocate(glat)
     allocate(glat(npry))

     idims(1) = npry
     call shdf5_irec(ndims, idims, 'glat' ,rvar1=glat)
  endif

  ! Check data domain size, location, and type

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
        call shdf5_close()
        stop 'astp stop1 - non-global domain in input sst/seaice data'
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
        stop 'astp stop2 - non-global domain in input sst/seaice data'
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

     stop 'astp stop - input sst/seaice dataset format'

  endif

  nio = nprx + 4
  njo = npry + 4

  allocate(ice(nio,njo))

  ! Check if seaice is in the analysis file, and read it

  call shdf5_info('ICEC', ndims, idims)

  if (ndims > 0) then
     allocate(a2d(nprx,npry))

     call shdf5_irec(ndims, idims, 'ICEC', rvar2 = a2d)
     call prfill(nprx, npry, a2d, ice, gdatdy, xswlat, ipoffset, inproj)

     has_seaice = .true.
     deallocate(a2d)
  endif

  call shdf5_close()

  if (.not. has_seaice) then
     write(io6,*) "read_seaice: Analysis file does not contain seaice."
     write(io6,*) "Stopping run."
     stop
  endif

  if (inproj == 2) then
     allocate(plat(npry+4))

     do ilat = 1,npry
        plat(ilat+2) = glat(ilat)
     enddo

     plat(2) = plat(3) - (plat(4) - plat(3))
     plat(1) = plat(2) - (plat(3) - plat(2))

     plat(npry+3) = plat(npry+2) + (plat(npry+2) - plat(npry+1))
     plat(npry+4) = plat(npry+3) + (plat(npry+3) - plat(npry+2))
  endif

  ! Fill seaice array

  do isea = 2, msea
     iwsfc = isea + omsea

! fractional x/y indices in pressure data arrays at current iw point location

     if (inproj == 1) then

        gry = (sfcg%glatw(iwsfc) - xswlat) / gdatdy + 3.
        grx = (sfcg%glonw(iwsfc) - xswlon) / gdatdx + 1. + real(ipoffset)

     elseif (inproj == 2) then

        ! estimate latitude index assuming uniform spacing of plat

        ilat = 2 + npry * int((sfcg%glatw(iwsfc) - plat(2)) / 180.)

        ! find correct latitude index

        if (plat(ilat) > sfcg%glatw(iwsfc)) then
           do while(plat(ilat) > sfcg%glatw(iwsfc))
              ilat = ilat - 1
           enddo
        elseif (plat(ilat+1) < sfcg%glatw(iwsfc)) then
           do while(plat(ilat+1) < sfcg%glatw(iwsfc))
              ilat = ilat + 1
           enddo
        endif

        gry = (sfcg%glatw(iwsfc) - plat(ilat)) / (plat(ilat+1) - plat(ilat)) + real(ilat)
        grx = (sfcg%glonw(iwsfc) - xswlon) / gdatdx + 1. + real(ipoffset)

     endif

     call gdtost(ice, nprx+4, npry+4, grx, gry, sea%seaicef(isea))

     sea%seaicef(isea) = max( min( sea%seaicef(isea), 1.0 ), 0.0 )

  enddo

  deallocate( ice )

  if (iaction == 0) then
     sea%seaicep(:) = sea%seaicef(:)
  endif

end subroutine read_seaice_analysis
