module point_io

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  use                   minterp_lib, only: minterp_wghts, get_weights_lonlat
  implicit none

  real,                 parameter  :: rmiss = -1.e34

  integer                          :: nobs, nobs_loc
  integer                          :: iobs_current
  character(8)                     :: date
  character(30)                    :: time_units
  character(30)                    :: time_longname

  integer,             allocatable :: itimes     (:) ! obspack uses integer secs since 1970
  real(r8),            allocatable :: s1900_times(:)
  real,                allocatable :: latitudes  (:)
  real,                allocatable :: longitudes (:)
  real,                allocatable :: altitudes  (:)
  integer,             allocatable :: iobs_loc   (:)
  type(minterp_wghts), allocatable :: wgts       (:)

  ! variables to be written each observation point
  real,                allocatable :: obs      (:,:) ! dimmed to (naddsc,nobs)
! real,                allocatable :: airtemp    (:)
! real,                allocatable :: airprss    (:)

contains


  subroutine read_point_file()

    use oname_coms,   only: nl
    use max_dims,     only: pathlen
    use misc_coms,    only: current_time, io6, iparallel, s1900_sim, naddsc
    use mem_para,     only: olam_stop, myrank
    use mem_ijtabs,   only: itabg_m
    use netcdf_utils, only: ncf_times_to_s1900
    use hdf5_utils,   only: shdf5_exists, shdf5_open, shdf5_close, shdf5_info, &
                            shdf5_irec
    implicit none

    character(pathlen)               :: filename
    integer                          :: iyear, imonth, idate, ihour
    logical                          :: exists, error
    integer                          :: ndims, idims(3), i, il
    integer,             allocatable :: itim_tmp(:)
    real(r8),            allocatable :: dtim_tmp(:)
    real,                allocatable :: alts_tmp(:)
    real,                allocatable :: lons_tmp(:)
    real,                allocatable :: lats_tmp(:)
    type(minterp_wghts), allocatable :: wgts_tmp(:)

    ! character string YYYYMMDD
    write(date,'(I4,I2.2,I2.2)') current_time%year, current_time%month, current_time%date

    filename = trim(nl%point_obsin_header) // date // trim(nl%point_obsin_suffix)

    call shdf5_exists(filename, exists)

    if (.not. exists) then
       write(io6,*) "point_io: Error opening data file " // trim(filename)
       write(io6,*) "File to read point output locations cannot be found."
       call olam_stop("Stopping model")
    endif

    call shdf5_open(filename, 'R')

    ndims = 0
    idims = 0
    time_units    = ' '
    time_longname = ' '

    call shdf5_info('time', ndims, idims)

    nobs = 0
    if (ndims == 1) nobs = idims(1)

    if (nobs > 0) then

       ! ObsPack files use integer seconds since 1970
       if (allocated(itimes)) deallocate(itimes)
       allocate(itimes(nobs))
       call shdf5_irec(ndims, idims, 'time', ivar1=itimes, units=time_units, long_name=time_longname)

       ! OLAM uses time as double precision real seconds since 1900
       if (allocated(s1900_times)) deallocate(s1900_times)
       allocate(s1900_times(nobs))
       s1900_times(:) = real( itimes(:), r8)

       if (allocated(latitudes)) deallocate(latitudes)
       allocate(latitudes(nobs))
       call shdf5_irec(ndims, idims, 'latitude', rvar1=latitudes)

       if (allocated(longitudes)) deallocate(longitudes)
       allocate(longitudes(nobs))
       call shdf5_irec(ndims, idims, 'longitude', rvar1=longitudes)

       if (allocated(altitudes)) deallocate(altitudes)
       allocate(altitudes(nobs))
       call shdf5_irec(ndims, idims, 'altitude', rvar1=altitudes)

    endif

    call shdf5_close()

    ! Convert time to seconds since 1900 to match OLAM

    call ncf_times_to_s1900(time_units, nobs, s1900_times, error)

    ! Get interpolation weights and M location that corresponds to each
    ! lat-lon observation

    if (allocated(wgts)) deallocate(wgts)
    allocate(wgts(nobs))
    call get_weights_lonlat(longitudes, latitudes, wgts, nobs)

    ! For a parallel run, only store the points that are primary on this rank

    nobs_loc = nobs

    if (iparallel == 1) then

       ! number of obs primary to this MPI rank
       nobs_loc = count( itabg_m( wgts(:)%imglobe )%irank == myrank )

       if (allocated(iobs_loc)) deallocate(iobs_loc)
       allocate(iobs_loc(nobs_loc))

       ! allocate temporary arrays to store information local to this MPI rank
       allocate(lons_tmp(nobs_loc))
       allocate(lats_tmp(nobs_loc))
       allocate(alts_tmp(nobs_loc))
       allocate(itim_tmp(nobs_loc))
       allocate(dtim_tmp(nobs_loc))
       allocate(wgts_tmp(nobs_loc))

       ! copy data that is local to this MPI rank to the temporary arrays
       if (nobs_loc > 0) then
          il = 0
          do i = 1, nobs
             if ( itabg_m( wgts(i)%imglobe )%irank == myrank ) then
                il = il + 1
                lons_tmp(il) = longitudes (i)
                lats_tmp(il) = latitudes  (i)
                alts_tmp(il) = altitudes  (i)
                wgts_tmp(il) = wgts       (i)
                dtim_tmp(il) = s1900_times(i)
                itim_tmp(il) = itimes     (i)
                iobs_loc(il) = i
             endif
          enddo
       endif

       ! replace the global arrays with the information only needed on this rank
       call move_alloc(lons_tmp, longitudes)
       call move_alloc(lats_tmp, latitudes)
       call move_alloc(alts_tmp, altitudes)
       call move_alloc(dtim_tmp, s1900_times)
       call move_alloc(itim_tmp, itimes)
       call move_alloc(wgts_tmp, wgts)

    endif

    ! reset obs counter and skip to current time if the model start does not
    ! correspond to the first observation in the input file

    iobs_current = 1

    if (nobs_loc > 0) then

       ! if observations are more then 10 minutes before simulation time
       if (s1900_times(1) < s1900_sim - 600.) then
          do i = 1, nobs_loc
             if (s1900_times(i) > s1900_sim - 600.) exit
          enddo
          iobs_current = i
       endif

    endif

    ! Allocate local space for point output observations. Fortran allows zero-size
    ! arrays if no points exist on a rank

    if (allocated(obs)) deallocate(obs)
    allocate(obs(naddsc,nobs_loc))
    obs(:,:) = rmiss

  end subroutine read_point_file



  subroutine fill_point_obs()

    use mem_grid!,   only: mza, lpw, zt, dzim
    use mem_ijtabs, only: itab_m, itabg_m
    use mem_addsc,  only: addsc
    use misc_coms,  only: s1900_sim, time_bias, naddsc

    implicit none

    integer :: i, inext
    integer :: k, kobs
    real    :: zwt_top, zwt_bot
    integer :: im, iw, iaddsc, n
    real    :: field(3)

    if (iobs_current <= nobs_loc) then

       if ( s1900_times(iobs_current) < s1900_sim + time_bias ) then

          do i = iobs_current+1, nobs_loc
             if (s1900_times(i) > s1900_sim + time_bias) exit
          enddo
          inext = i - 1

          do i = iobs_current, inext

             do k = 2, mza
                if (altitudes(i) < zt(k)) exit
             enddo

             ! point is between zt(k-1) and zt(k)
             kobs = k - 1
             zwt_top = (altitudes(i) - zt(kobs)) * dzim(kobs)
             zwt_bot = 1.0 - zwt_top

             ! M point corresponding to this obs
             im = itabg_m( wgts(i)%imglobe )%im_myrank

             ! interpolate from 3 surrounding W cells
             ! do n = 1, 3
             !    iw = itab_m(im)%iw(i)
             !
             !    if (kobs >= mza) then
             !       field(n) = tair(mza,iw)
             !    elseif (kobs <= lpw(iw)) then
             !       field(n) = tair(lpw(iw),iw)
             !    else
             !       field(n) = zwt_top * tair(kobs+1,iw) + zwt_bot * tair(kobs,iw)
             !    endif
             ! enddo
             !
             ! airtemp(i) = sum( field(1:3) * wgts(i)%wt(1:3) )

             do iaddsc = 1, naddsc

                ! interpolate from 3 surrounding W cells
                do n = 1, 3
                   iw = itab_m(im)%iw(n)

                   if (kobs >= mza) then
                      field(n) = addsc(iaddsc)%sclp(mza,iw)
                   elseif (kobs < lpw(iw)) then
                      field(n) = addsc(iaddsc)%sclp(lpw(iw),iw)
                   else
                      field(n) = zwt_top * addsc(iaddsc)%sclp(kobs+1,iw) &
                               + zwt_bot * addsc(iaddsc)%sclp(kobs  ,iw)
                   endif
                enddo

                obs(iaddsc,i) = sum( field(1:3) * wgts(i)%wt(1:3) )

             enddo

          enddo

          iobs_current = inext + 1

       end if

    end if

  end subroutine fill_point_obs



  subroutine write_point_file

      use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
      use misc_coms,  only: iclobber, naddsc
      use max_dims,   only: pathlen
      use oname_coms, only: nl

      implicit none

      character(pathlen)   :: filename
      integer              :: ndims, idims(2), i
      integer, allocatable :: iobs_number(:)
      integer, allocatable :: itrc_number(:)

      filename = trim(nl%point_obsout_header) // date // trim(nl%point_obsout_suffix)

      call shdf5_open(filename, 'W', iclobber)

      ! Write obs # as a coordinate variable

      ndims    = 1
      idims(1) = nobs

      allocate(iobs_number(nobs))
      iobs_number = [ (i, i=1, nobs) ]

      call shdf5_orec(ndims, idims, "obs", ivar1=iobs_number, isdim=.true., long_name = "obs number")

      ! Write tracer # as a coordinate variable (for A. Schuh)

      if (naddsc > 0) then

         ndims    = 1
         idims(1) = naddsc

         allocate(itrc_number(naddsc))
         iobs_number = [ (i, i=1, naddsc) ]

         call shdf5_orec(ndims, idims, "tracer", ivar1=itrc_number, isdim=.true., long_name = "tracer number")

      endif

      ! Write lat/lon/time. These are not coordinate variable for point output with NetCDF!

      ndims = 1
      idims(1) = nobs_loc

      call shdf5_orec( ndims, idims, "time", ivar1=itimes, &
                       gpoints   = iobs_loc,      &
                       nglobe    = nobs,          &
                       dimnames  = [ "obs" ],     &
                       long_name = time_longname, &
                       units     = time_units     )

      call shdf5_orec( ndims, idims, "longitude", rvar1=longitudes, &
                       gpoints   = iobs_loc,      &
                       nglobe    = nobs,          &
                       dimnames  = [ "obs" ],     &
                       long_name = "longitude",   &
                       units     = "degrees_east" )

      call shdf5_orec( ndims, idims, "latitude", rvar1=latitudes, &
                       gpoints   = iobs_loc,       &
                       nglobe    = nobs,           &
                       dimnames  = [ "obs" ],      &
                       long_name = "latitude",     &
                       units     = "degrees_north" )

      call shdf5_orec( ndims, idims, "altitude", rvar1=altitudes, &
                       gpoints   = iobs_loc,       &
                       nglobe    = nobs,           &
                       dimnames  = [ "obs" ],      &
                       long_name = "altitude in meters above sea level", &
                       units     = "meters" )

      ! Write point observations. For now these are just the tracers for A. Schuh

      if (naddsc > 0) then

         ndims = 2
         idims = [ naddsc, nobs_loc ]

         call shdf5_orec( ndims, idims, "CO2_Tracers", rvar2=obs, &
                          gpoints = iobs_loc,                     &
                          nglobe = nobs,                          &
                          dimnames = ["tracer", "obs   "],        &
                          rmissing = rmiss,                       &
                          long_name = "simulated mole fraction of carbon dioxide in dry air" )

      endif

      call shdf5_close()

    end subroutine write_point_file


end module point_io
