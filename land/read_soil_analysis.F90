subroutine read_soil_analysis(soil_tempc)

  use misc_coms,  only: io6, s1900_sim, s1900_init, isubdomain, iparallel
  use leaf_coms,  only: dt_leaf, wcap_min
  use mem_land,   only: land, mland, omland, nzg, slzt
  use mem_sfcg,   only: sfcg, itab_wsfc
  use consts_coms,only: pio180, piu180, erad, cliq1000, alli1000, cice, &
                         cice1000, r8
  use max_dims,   only: pathlen
  use isan_coms,  only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, glat, &
                        inproj, xswlat, xswlon, gdatdx, gdatdy, ipoffset
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use mem_para,   only: myrank, nbytes_int, nbytes_real
  use prfill_mod, only: prfill, prfill3
  use mem_ijtabs, only: itab_w

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real, intent(inout) :: soil_tempc(nzg,mland)

  integer            :: ifgfile, nf
  character(pathlen) :: fname
  character(16)      :: ext
  logical            :: exists
  integer            :: ndims, idims(3)
  real               :: grx, gry
  real               :: snowdens, mass, tempc, tempk
  integer            :: nio, njo, ngnd
  integer            :: iland, iwsfc, i, j, k, kk, ntext
  logical            :: has_snow, has_soilt, has_soilw
  integer            :: bytes, isize, ier, igloberr, ilat, ipry

  real, allocatable  :: snow (:,:)   ! snow mass  [kg/m2]
  real, allocatable  :: soilt(:,:,:) ! soil temp  [K]
  real, allocatable  :: soilw(:,:,:) ! soil water [Vol. fraction]

  real, allocatable  :: a2d(:,:), a3d(:,:,:)
  real, allocatable  :: tcol(:), wcol(:), tcol2(:), wcol2(:)
  real, allocatable  :: zcol(:), ztmp(:)
  real               :: wprof(nzg)

  integer, allocatable :: buffer(:)
  real,    allocatable :: plat(:)

  has_snow  = .false.
  has_soilt = .false.
  has_soilw = .false.
  ngnd      =  0

! Loop over analysis files and search for the one that corresponds
! to the model start time

  ifgfile = 0
  do nf = 1, nfgfiles
     if (s1900_fg(nf) <= s1900_sim) then
        ifgfile = nf
     endif
  enddo

  if (ifgfile < 1) then
     write(io6,*) ' '
     write(io6,*) 'Unable to find analysis file for soil initialization'
     write(io6,*) 'Using default initialization'
     return
  elseif (s1900_fg(ifgfile) < s1900_init - 1800.0_r8) then
     write(io6,*) ' '
     write(io6,*) 'Available analysis file is more than 1800 seconds'
     write(io6,*) 'prior to simulation initialization time.'
     write(io6,*) 'Ising default soil initialization'
     return
  endif

  fname = fnames_fg(ifgfile)

! Make sure we have a HDF5 analysis file

  ext = trim( fname( index(fname,'.',back=.true.)+1:) )

  if ( .not. ( &
       ext == 'h5' .or. ext == 'hdf5' .or. ext == 'hdf' .or. &
       ext == 'H5' .or. ext == 'HDF5' .or. ext == 'HDF' )) then
     write(io6,*) "read_soil: Analysis file is not hdf5 format"
     write(io6,*) "Using default soil initialization instead."
     return
  endif

! Process selected analysis file

  inquire(file=fname, exist=exists)
  if (.not. exists) then
     write(io6,*) "read_soil: Error opening analysis file " // trim(fname)
     write(io6,*) "Using default soil initialization instead."
     return
  endif

  write(io6,'(A)') ' read_soil: opening ' // trim(fname)

#ifdef OLAM_MPI
  if (iparallel == 1) then
     bytes = 0
     isize = nbytes_int*4 + nbytes_real*4
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

     if (inproj == 2) then
        if (allocated(glat)) deallocate(glat)
        allocate(glat(npry))

        idims(1) = npry
        call shdf5_irec(ndims, idims, 'glat' ,rvar1=glat)
     endif

     ! Check if ngnd, the # of soil levels, is in the analysis file and read it

     ngnd = 0
     call shdf5_info('ngnd', ndims, idims)
     if (ndims > 0) call shdf5_irec(ndims, idims, 'ngnd' , ivars=ngnd)

     ! Check if sdepths, the soil depth array, is in the analysis file and read it

     call shdf5_info('sdepths', ndims, idims)
     if (ndims > 0) then
        if (ngnd == 0) ngnd = idims(1)
        allocate( ztmp(ngnd), zcol(ngnd) )
        call shdf5_irec(ndims, idims, 'sdepths', rvar1=ztmp)
     endif

     ! OLAM stores the soil arrays from bottom to top, so we need to reverse
     ! the input soil depth array, and convert to m

     do k = 1, ngnd
        kk = ngnd - k + 1
        zcol(kk) = ztmp(k) * 0.01
     enddo

     deallocate(ztmp)

#ifdef OLAM_MPI
     if (iparallel == 1) then
        call MPI_Pack(nprx  , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(npry  , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(inproj, 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlat, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(xswlon, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdx, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(gdatdy, 1, MPI_REAL   , buffer, isize, bytes, MPI_COMM_WORLD, ier)
        call MPI_Pack(ngnd  , 1, MPI_INTEGER, buffer, isize, bytes, MPI_COMM_WORLD, ier)
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
        call MPI_Unpack(buffer, isize, bytes, ngnd  , 1, MPI_INTEGER, MPI_COMM_WORLD, ier)
        if (ngnd > 0) allocate(zcol(ngnd))
     endif

     if (ngnd > 0) call MPI_Bcast(zcol, ngnd, MPI_REAL, 0, MPI_COMM_WORLD, ier)

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

  if (ngnd == 0) then
     write(io6,*) "Gridded analysis dataset does not contain soil depth information."
     write(io6,*) "Soil initialization from analysis file will be skipped."
     write(io6,*)
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
        if (myrank == 0) call shdf5_close()
        stop 'astp stop1 - non-global domain in input soil data'
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
        write(io6,*) 'Gridded soil data must have global coverage'
        write(io6,*) 'nprx,npry = ',nprx,npry
        write(io6,*) 'gdatdx = ',gdatdx
        write(io6,*) 'xswlat,xswlon= ',xswlat,xswlon
        do ipry = 1,npry
           write(io6,*) 'ipry, glat = ',ipry,glat(ipry)
        enddo
        stop 'astp stop2 - non-global domain in input soil data'
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

     stop 'astp stop - input soil dataset format'

  endif

  nio = nprx + 4
  njo = npry + 4

  if (myrank == 0) then

     ! Check if snow mass is in the analysis file, and read it

     call shdf5_info('SNOWMASS', ndims, idims)

     if (ndims > 0) then

        allocate(a2d(nprx,npry))
        allocate(snow(nio,njo))

        call shdf5_irec(ndims, idims, 'SNOWMASS', rvar2 = a2d)
        call prfill(nprx, npry, a2d, snow, gdatdy, xswlat, ipoffset, inproj)

        has_snow = .true.
        deallocate(a2d)
     endif

     ! Check if soil temperature and/or soil moisture are in the analysis file,
     ! and read them

     if (ngnd > 0 .and. allocated(zcol)) then
        allocate(a3d(nprx,npry,ngnd))

        call shdf5_info('SOILT', ndims, idims)

        if (ndims > 0) then
           allocate(soilt(nio,njo,ngnd))

           call shdf5_irec(ndims, idims, 'SOILT', rvar3 = a3d)
           call prfill3(nprx, npry, ngnd, a3d, soilt, gdatdy, xswlat, ipoffset, inproj)
           has_soilt = .true.
        endif

        call shdf5_info('SOILW', ndims, idims)

        if (ndims > 0) then
           allocate(soilw(nio,njo,ngnd))

           call shdf5_irec(ndims, idims, 'SOILW', rvar3 = a3d)
           call prfill3(nprx, npry, ngnd, a3d, soilw, gdatdy, xswlat, ipoffset, inproj)
           has_soilw = .true.
        endif

        deallocate(a3d)
     endif

     call shdf5_close()
  endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(has_snow , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Bcast(has_soilt, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Bcast(has_soilw, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  if (myrank > 0 .and. has_snow ) allocate(snow (nio,njo))
  if (myrank > 0 .and. has_soilt) allocate(soilt(nio,njo,ngnd))
  if (myrank > 0 .and. has_soilw) allocate(soilw(nio,njo,ngnd))

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Bcast(snow,  nio*njo     , MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Bcast(soilt, nio*njo*ngnd, MPI_REAL, 0, MPI_COMM_WORLD, ier)
     call MPI_Bcast(soilw, nio*njo*ngnd, MPI_REAL, 0, MPI_COMM_WORLD, ier)
  endif
#endif

  if (.not. has_snow) then
     write(io6,*) "read_soil: Analysis file does not contain snow mass."
     write(io6,*) "Skipping snow depth initialization."
     write(io6,*)
  endif

  if (.not. has_soilt) then
     write(io6,*) "read_soil: Analysis file does not contain soil temperature."
     write(io6,*) "Using default soil temperatures instead."
     write(io6,*)
  endif

  if (.not. has_soilw) then
     write(io6,*) "read_soil: Analysis file does not contain soil moisture."
     write(io6,*) "Using default soil water instead."
     write(io6,*)
  endif

  ! Make sure that the input fields don't contain all missing data. This
  ! can happen if we specified to save the soil fields during degribbing,
  ! but the grib file didn't have any soil or snow fields

  if (has_snow) then
     if ( count( snow(:,:) < -1.0 .or. snow(:,:) > 1.e30) >= nprx*npry ) then
        write(io6,*) "read_soil: Snow mass from analysis contains mostly missing data."
        write(io6,*) "Skipping snow depth initialization."
        has_snow = .false.
     endif
  endif

  if (has_soilt) then
     if ( count( soilt(:,:,:) < -1.0  .or. soilt(:,:,:) > 1.e30) >= nprx*npry*ngnd ) then
        write(io6,*) "read_soil: Soil temperature from analysis contains mostly missing data."
        write(io6,*) "Skipping soil temperature initialization."
        has_soilt = .false.
     endif
  endif

  if (has_soilw) then
     if ( count( soilw(:,:,:) < -1.0 .or. soilw(:,:,:) > 1.e30) >= nprx*npry*ngnd ) then
        write(io6,*) "read_soil: Soil moisture from analysis contains mostly missing data."
        write(io6,*) "Skipping soil moisture initialization."
        has_soilw = .false.
     endif
  endif

  ! If all data missing, return

  if ((.not. has_snow) .and. (.not. has_soilt) .and. (.not. has_soilw)) return

  ! Set masks to indicate missing data

  if (has_snow) then
     do i = 1, nio
        do j = 1, njo
           if (snow(i,j) < -1.0 .or. snow(i,j) > 1.e30) then
              snow(i,j) = 2.e30
           else
              snow(i,j) = max( snow(i,j), 0.0 )
           endif
        enddo
     enddo
  endif

  if (has_soilt) then
     do i = 1, nio
        do j = 1, njo
           do k = 1, ngnd
              if ( soilt(i,j,k) < 100.0 .or. soilt(i,j,k) > 1.e30 ) then
                 soilt(i,j,k) = 2.e30
              else
                 soilt(i,j,k) = max( min(soilt(i,j,k), 340.0), 200.0)
              endif
           enddo
        enddo
     enddo
  endif

  if (has_soilw) then
     do i = 1, nio
        do j = 1, njo
           do k = 1, ngnd
              if ( soilw(i,j,k) < -1.0 .or. soilw(i,j,k) > 1.e30 ) then
                 soilw(i,j,k) = 2.e30
              else
                 soilw(i,j,k) = max(soilw(i,j,k), 1.e-10)
              endif
           enddo
        enddo
     enddo
  endif

  allocate(tcol (ngnd))
  allocate(tcol2(ngnd))
  allocate(wcol (ngnd))
  allocate(wcol2(ngnd))

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

  ! Fill soil arrays

  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     !if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle
     if (isubdomain == 1) then
        if (.not. any( itab_w( itab_wsfc(iwsfc)%iwatm( 1:itab_wsfc(iwsfc)%nwatm ) )%irank == myrank ) ) cycle
     endif

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

     ! Interpolate soil temperature to this land cell if any of the 4 closest
     ! analysis points have non-missing temperature data

     if (has_soilt) then

        tcol(1:ngnd) = -999.

        do k = 1, ngnd
           ! where soilt is masked out, tcol will be returned as missing
           call gdtost(soilt(:,:,k), nprx+4, npry+4, grx, gry, tcol(k))
        enddo

        if ( all(tcol(1:ngnd) < 1.e20 .and. tcol(1:ngnd) > 0.) ) then

           tcol(1:ngnd) = tcol(1:ngnd) - 273.15

           if (ngnd == 1) then

              soil_tempc(1:nzg,iland) = tcol(1)

           else

              ! OLAM stores the soil arrays from bottom to top, so
              ! the input soil array needs to be reversed
              do k = 1, ngnd
                 kk = ngnd - k + 1
                 tcol2(kk)  = tcol(k)
              enddo

              call hintrp_cc( ngnd, tcol2, zcol, nzg, soil_tempc(:,iland), slzt )
           endif

        endif
     endif

     ! Interpolate soil moisture to this land cell if any of the 4 closest
     ! analysis points have non-missing soil moisture data

     if (has_soilw) then

        wcol(1:ngnd) = -999.

        do k = 1, ngnd
           ! where soilw is masked out, wcol will be returned as missing
           call gdtost(soilw(:,:,k), nprx+4, npry+4, grx, gry, wcol(k))
        enddo

        if ( all(wcol(1:ngnd) < 0.5 .and. wcol(1:ngnd) >= 0.) ) then

           if (ngnd == 1) then

              wprof(1:nzg)= wcol(1)

           else

              ! OLAM stores the soil arrays from bottom to top, so
              ! the input soil array needs to be reversed
              do k = 1, ngnd
                 kk = ngnd - k + 1
                 wcol2(kk)  = wcol(k)
              enddo

              call hintrp_cc( ngnd, wcol2, zcol, nzg, wprof, slzt )
           endif

           ! Bound soil moisture between residual and saturation values
           ! Soil moisture is only modified above the water table depth

           do k = 1, nzg
              if (land%head0(iland) - slzt(k) < 0.) then
                 land%soil_water(k,iland) = max( 1.03*land%wresid_vg(k,iland), &
                                            min( wprof(k), .97*land%wsat_vg(k,iland) ))
              endif
           enddo

        endif
     endif

     ! Interpolate snow depth to this land cell if any of the 4 closest
     ! analysis points have non-missing snow data

     if (has_snow) then

        ! where snow is masked out, mass will be returned as missing
        call gdtost(snow, nprx+4, npry+4, grx, gry, mass)

        if (mass > wcap_min .and. mass < 1.e20) then

           land%sfcwater_mass(1,iland) = mass

           tempc = sfcg%cantemp(iwsfc) - 273.15
           if (has_soilt) tempc = min(tempc, soil_tempc(nzg,iland))

           land%sfcwater_energy(1,iland) = min(0., tempc * cice)

           ! snow density calculation comes from CLM3.0 documentation
           ! which is based on Anderson 1975 NWS Technical Doc # 19

           snowdens = 50.0
           tempk    = tempc + 273.15
           if (tempk > 258.15) snowdens = 50.0 + 1.5 * (tempk - 258.15)**1.5

           land%sfcwater_depth(1,iland) = land%sfcwater_mass(1,iland) / snowdens

        else

           land%sfcwater_mass  (1,iland) = 0.
           land%sfcwater_energy(1,iland) = 0.
           land%sfcwater_depth (1,iland) = 0.

        endif

     endif

  enddo

end subroutine read_soil_analysis
