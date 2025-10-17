subroutine read_soil_analysis(soil_tempc)

  use misc_coms,   only: io6, s1900_sim, s1900_init, iparallel
  use leaf_coms,   only: wcap_min
  use mem_land,    only: land, mland, omland, nzg, slzt
  use mem_sfcg,    only: sfcg, itab_wsfc
  use consts_coms, only: pio180, piu180, cliq1000, alli1000, cice, &
                         cice1000, r8
  use max_dims,    only: pathlen
  use isan_coms,   only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, plats, &
                         inproj, xswlat, xswlon, gdatdx, gdatdy, irev_ns, &
                         read_analysis_header
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use mem_para,    only: myrank
  use analysis_lib,only: gdtost_ll

  implicit none

  real, intent(inout) :: soil_tempc(nzg,mland)

  integer            :: ifgfile, nf
  character(pathlen) :: fname
  logical            :: exists
  integer            :: ndims, idims(3)
  real               :: snowdens, tempc, tempk
  integer            :: ngnd, iland, iwsfc, k, kk

  real, allocatable  :: snow (:)   ! snow mass  [kg/m2]
  real, allocatable  :: soilt(:,:) ! soil temp  [K]
  real, allocatable  :: soilw(:,:) ! soil water [Vol. fraction]

  real, allocatable  :: a2d(:,:)
  real, allocatable  :: tcol(:), wcol(:)
  real, allocatable  :: zcol(:), ztmp(:)
  real, allocatable  :: wprof(:)

  ngnd = 0

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

! Process selected analysis file

  call shdf5_exists(fname, exists)

  if (.not. exists) then
     write(io6,*) "read_soil: Error opening analysis file " // trim(fname)
     write(io6,*) "Using default soil initialization instead."
     return
  endif

  write(io6,'(A)') ' read_soil: opening ' // trim(fname)

  call shdf5_open(fname, 'R', mpio_collective_read=.false.)

  call read_analysis_header(noplevs=.true.)

  ! Check if ngnd, the # of soil levels, is in the analysis file and read it

  ngnd = 0
  call shdf5_info('ngnd', ndims, idims)
  if (ndims > 0) call shdf5_irec(ndims, idims, 'ngnd' , ivars=ngnd)

  ! Check if sdepths, the soil depth array, is in the analysis file and read it

  call shdf5_info('sdepths', ndims, idims)
  if (ndims > 0) then

     if (ngnd /= idims(1)) then
        ngnd = idims(1)
     endif

     allocate( ztmp(ngnd), zcol(ngnd) )
     call shdf5_irec(ndims, idims, 'sdepths', rvar1=ztmp)

     ! OLAM stores the soil arrays from bottom to top, so we need to reverse
     ! the input soil depth array, and convert to m

     do k = 1, ngnd
        kk = ngnd - k + 1
        zcol(kk) = ztmp(k) * 0.01
     enddo

     deallocate(ztmp)

  else

     ngnd = 0

  endif

  allocate(a2d(nprx,npry))

  ! Check if snow mass is in the analysis file, and read it

  call shdf5_info('SNOWMASS', ndims, idims)

  if (ndims > 0) then
     allocate(snow(mland))

     call shdf5_irec(ndims, idims, 'SNOWMASS', rvar2=a2d)

     where(a2d < -1.0 .or. a2d > 1.e30)
        a2d = -1.e30
     elsewhere
        a2d = max( a2d, 0.0 )
     endwhere

     !$omp parallel do private(iwsfc)
     do iland = 2, mland
        iwsfc = iland + omland

        call gdtost_ll( nprx, npry, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                        snow(iland), xswlon, xswlat, gdatdx, gdatdy, inproj, &
                        r2d=a2d, plats=plats, irev_ns=irev_ns )
     enddo
     !$omp end parallel do

  endif

  if (ngnd > 0 .and. allocated(zcol)) then

     ! Check if soil temperature is in the analysis file, and read it

     call shdf5_info('SOILT', ndims, idims)

     if (ndims > 0) then
        allocate(soilt(ngnd,mland))

        do k = 1, ngnd

           call shdf5_irec( ndims, idims, 'SOILT', rvar2=a2d, &
                            start=[1,1,k], counts=[nprx,npry,1] )

           where(a2d < 100. .or. a2d > 1.e30)
              a2d = -1.e30
           elsewhere
              a2d = max( min(a2d, 340.0), 200.0)
           endwhere

           !$omp parallel do private(iwsfc)
           do iland = 2, mland
              iwsfc = iland + omland

              call gdtost_ll( nprx, npry, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                              soilt(k,iland), xswlon, xswlat, gdatdx, gdatdy, inproj, &
                              r2d=a2d, plats=plats, irev_ns=irev_ns )
           enddo
           !$omp end parallel do

        enddo

     endif

     ! Check if soil moisture is in the analysis file, and read it

     call shdf5_info('SOILW', ndims, idims)

     if (ndims > 0) then
        allocate(soilw(ngnd,mland))

        do k = 1, ngnd

           call shdf5_irec( ndims, idims, 'SOILW', rvar2=a2d, &
                            start=[1,1,k], counts=[nprx,npry,1] )

           where(a2d < -1.0 .or. a2d > 1.e30)
              a2d = -1.e30
           elsewhere
              a2d = max( min(a2d, 1.0), 0.0 )
           endwhere

           !$omp parallel do private(iwsfc)
           do iland = 2, mland
              iwsfc = iland + omland

              call gdtost_ll( nprx, npry, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                              soilw(k,iland), xswlon, xswlat, gdatdx, gdatdy, inproj, &
                              r2d=a2d, plats=plats, irev_ns=irev_ns )
           enddo
           !$omp end parallel do

        enddo

     endif

  endif

  call shdf5_close()

  ! Make sure that the input fields don't contain all missing data. This
  ! can happen if we specified to save the soil fields during degribbing,
  ! but the grib file didn't have any soil or snow fields

  if (allocated(snow)) then
     write(io6,*) "read_soil: Initializing snow depth from analysis file."
     write(io6,*)
  else
     write(io6,*) "read_soil: Analysis file does not contain snow mass."
     write(io6,*) "Skipping snow depth initialization."
     write(io6,*)
  endif

  if (allocated(soilt)) then
     write(io6,*) "read_soil: Initializing soil temperature from analysis file."
     write(io6,*)
  else
     write(io6,*) "read_soil: Analysis file does not contain soil temperature."
     write(io6,*) "Using default soil temperatures instead."
     write(io6,*)
  endif

  if (allocated(soilw)) then
     write(io6,*) "read_soil: Initializing soil water from analysis file."
     write(io6,*)
  else
     write(io6,*) "read_soil: Analysis file does not contain soil moisture."
     write(io6,*) "Using default soil water instead."
     write(io6,*)
  endif

  ! If all data missing, return

  if ( (.not. allocated(snow) ) .and. &
       (.not. allocated(soilt)) .and. &
       (.not. allocated(soilw)) ) return

  ! Fill soil arrays

  allocate(tcol(ngnd))
  allocate(wcol(ngnd))
  allocate(wprof(nzg))

  !$omp parallel private(tcol,wcol,wprof)
  !$omp do private(iwsfc,k,kk,tempc,snowdens,tempk)
  do iland = 2, mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     if (allocated(soilt)) then

        if ( all( soilt(:,iland) < 350. .and. soilt(:,iland) > 190. ) ) then

           if (ngnd == 1) then

              soil_tempc(1:nzg,iland) = soilt(1,iland) - 273.15

           else

              ! OLAM stores the soil arrays from bottom to top, so
              ! the input soil array needs to be reversed
              do k = 1, ngnd
                 kk = ngnd - k + 1
                 tcol(kk) = soilt(k,iland) - 273.15
              enddo

              call hintrp_cc( ngnd, tcol, zcol, nzg, soil_tempc(:,iland), slzt )
           endif

        endif
     endif

     ! Interpolate soil moisture to this land cell if any of the 4 closest
     ! analysis points have non-missing soil moisture data

     if (allocated(soilw)) then

        if ( all( soilw(:,iland) <= 1.0 .and. soilw(:,iland) >= 0. ) ) then

           if (ngnd == 1) then

              wprof(1:nzg) = soilw(1,iland)

           else

              ! OLAM stores the soil arrays from bottom to top, so
              ! the input soil array needs to be reversed
              do k = 1, ngnd
                 kk = ngnd - k + 1
                 wcol(kk) = soilw(k,iland)
              enddo

              call hintrp_cc( ngnd, wcol, zcol, nzg, wprof, slzt )
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

     if (allocated(snow)) then

        if (snow(iland) > wcap_min .and. snow(iland) < 1.e20) then

           land%sfcwater_mass(1,iland) = snow(iland)

           tempc = sfcg%cantemp(iwsfc) - 273.15
           tempc = min(tempc, soil_tempc(nzg,iland))

           land%sfcwater_epm2(1,iland) = land%sfcwater_mass(1,iland) * min(0., tempc * cice)

           ! snow density calculation comes from CLM3.0 documentation
           ! which is based on Anderson 1975 NWS Technical Doc # 19

           ! NOTE: this is only appropriate for newly fallen snow, NOT existing snowpack
           !snowdens = 50.0
           !tempk    = tempc + 273.15
           !if (tempk > 258.15) snowdens = 50.0 + 1.5 * (tempk - 258.15)**1.5

           ! go from 300 to 650 kg/m^3
           snowdens = 250. + 5.0 * max(0., min(17., tempc+15.))**1.5

           land%sfcwater_depth(1,iland) = land%sfcwater_mass(1,iland) / snowdens

        else

           land%sfcwater_mass (1,iland) = 0.
           land%sfcwater_epm2 (1,iland) = 0.
           land%sfcwater_depth(1,iland) = 0.

        endif

     endif

  enddo
  !$omp end do
  !$omp end parallel

end subroutine read_soil_analysis
