module point_io

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  use                   minterp_lib, only: minterp_wghts, get_weights_lonlat
  implicit none

  real, parameter  :: rmiss  = -1.e34
  integer          :: nfiles = 0
  character(8)     :: date   = ' '

  Type point_data
     character(1)                     :: type
     integer                          :: nobs
     integer                          :: nobs_loc
     integer                          :: iobs_current
     character(30)                    :: time_units
     character(30)                    :: time_longname
     character(30)                    :: time_dimname
     real(r8),            allocatable :: dtimes     (:)
     real(r8),            allocatable :: s1900_times(:)
     real,                allocatable :: latitude   (:)
     real,                allocatable :: longitude  (:)
     real,                allocatable :: altitude   (:)
     real,                allocatable :: pressure   (:)
     integer,             allocatable :: iobs_loc   (:)
     type(minterp_wghts), allocatable :: wgts       (:)
     real,                allocatable :: obs(:,:)  ! dimmed to (naddsc,nobs_loc)
     character(30)                    :: obs_dimname
     character(30)                    :: obs_longname
     real(r8),            allocatable :: obs_dimvals(:)
  End type point_data

  type(point_data), allocatable, target :: pfiles(:)


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
    use string_lib,   only: lowercase

    implicit none

    character(pathlen)               :: filename
    integer                          :: nfiles, nf
    logical                          :: exists, error
    integer                          :: ndims, idims(3), i, il
    real(r8),            allocatable :: tims_tmp(:)
    real(r8),            allocatable :: dtim_tmp(:)
    real,                allocatable :: alts_tmp(:)
    real(r8),            allocatable :: dims_tmp(:)
    real,                allocatable :: prss_tmp(:)
    real,                allocatable :: lons_tmp(:)
    real,                allocatable :: lats_tmp(:)
    type(minterp_wghts), allocatable :: wgts_tmp(:)
    logical                          :: firstime = .true.
    type(point_data),        pointer :: pf
    integer,                 pointer :: nobs, nobs_loc

    character(30) :: attached_dim(1)
    character(30) :: units

    nfiles = nl%point_files

    if (firstime .and. nfiles > 0) then
       allocate(pfiles(nfiles))
       firstime = .false.
    endif

    ! character string YYYYMMDD
    write(date,'(I4,I2.2,I2.2)') current_time%year, current_time%month, current_time%date

    do nf = 1, nfiles
       pf       => pfiles(nf)
       nobs     => pfiles(nf)%nobs
       nobs_loc => pfiles(nf)%nobs_loc

       nobs     = 0
       nobs_loc = 0

       pf%time_units    = ' '
       pf%time_longname = ' '
       pf%type          = ' '

       pf%type = nl%pfile_types(nf)

       if (.not. any( pf%type == ['Z', 'P', 'S', 'C'] ) ) then
          write(*,'(A,I0)') " Illegal point file type for point_file # ", nf
          write(*,*)        "Skipping this file"
          cycle
       endif

       filename = trim(nl%point_obsin_headers(nf)) // date // trim(nl%point_obsin_suffixes(nf))

       call shdf5_exists(filename, exists)

       if (.not. exists) then
          write(io6,'(A)') " point_io: Error opening data file " // trim(filename)
          write(io6,'(A)') " Skipping this date."
          cycle
       endif

       call shdf5_open(filename, 'R')

       ndims = 0
       idims = 0
       call shdf5_info('time', ndims, idims)

       if (ndims == 1) nobs = idims(1)

       if (nobs > 0) then

          pf%time_dimname = ' '
          pf%obs_dimname = ' '
          pf%time_longname = ' '
          pf%obs_longname = ' '

          ! ObsPack files use integer seconds since 1970
          if (allocated(pf%dtimes)) deallocate(pf%dtimes)
          allocate(pf%dtimes(nobs))
          call shdf5_irec(ndims, idims, 'time', dvar1=pf%dtimes, units=pf%time_units, long_name=pf%time_longname, &
               dimname=pf%time_dimname, attached_dimnames=attached_dim)

          if (len_trim(pf%time_dimname) == 0 .and. len_trim(attached_dim(1)) > 0) then
             pf%obs_dimname = attached_dim(1)
          endif

          ! OLAM uses time as double precision real seconds since 1900
          if (allocated(pf%s1900_times)) deallocate(pf%s1900_times)
          allocate(pf%s1900_times(nobs))
          pf%s1900_times(:) = pf%dtimes(:)

          if (allocated(pf%latitude)) deallocate(pf%latitude)
          allocate(pf%latitude(nobs))
          call shdf5_irec(ndims, idims, 'latitude', rvar1=pf%latitude)

          if (allocated(pf%longitude)) deallocate(pf%longitude)
          allocate(pf%longitude(nobs))
          call shdf5_irec(ndims, idims, 'longitude', rvar1=pf%longitude)

          if (pf%type == 'Z') then

             call shdf5_info('altitude', ndims, idims)
             if (ndims == 0 .and. idims(1) == 0) then
                write(*,*) "Altitude cannot be found in NetCDF file for a Z (height) point file"
                write(*,*) "when reading ", trim(filename)
                write(*,*) "Skipping this file"
                cycle
             endif

             if (allocated(pf%altitude)) deallocate(pf%altitude)
             allocate(pf%altitude(nobs))
             call shdf5_irec(ndims, idims, 'altitude', rvar1=pf%altitude)

          elseif (pf%type == 'P') then

             call shdf5_info('pressure', ndims, idims)
             if (ndims == 0 .and. idims(1) == 0) then
                write(*,*) "Pressure cannot be found in NetCDF file for a P (pressure) point file"
                write(*,*) "when reading ", trim(filename)
                write(*,*) "Skipping this file"
                cycle
             endif

             if (allocated(pf%pressure)) deallocate(pf%pressure)
             allocate(pf%pressure(nobs))
             call shdf5_irec(ndims, idims, 'pressure', rvar1=pf%pressure, units=units)
             call lowercase(units)
             if (units == 'mb' .or. units == 'millibars') pf%pressure = pf%pressure * 100.

          endif

          if (len_trim(pf%obs_dimname) > 0) then
             call shdf5_info(pf%obs_dimname, ndims, idims)
             if (ndims == 0 .and. idims(1) /= nobs) then
                write(*,*) trim(pf%obs_dimname) // " cannot be found in NetCDF file."
                write(*,*) "Assuming no obs coordinate variable."
                pf%obs_dimname = " "
             else
                if (allocated(pf%obs_dimvals)) deallocate(pf%obs_dimvals)
                allocate(pf%obs_dimvals(nobs))
                call shdf5_irec(ndims, idims, pf%obs_dimname, dvar1=pf%obs_dimvals, long_name=pf%obs_longname)
             endif
          endif

       endif

       call shdf5_close()

       ! Convert time to seconds since 1900 to match OLAM

       call ncf_times_to_s1900(pf%time_units, nobs, pf%s1900_times, error)

       ! Get interpolation weights and M location that corresponds to each
       ! lat-lon observation

       if (allocated(pf%wgts)) deallocate(pf%wgts)
       allocate(pf%wgts(nobs))
       call get_weights_lonlat(pf%longitude, pf%latitude, pf%wgts, nobs)

       ! For a parallel run, only store the points that are primary on this rank

       nobs_loc = nobs

       if (iparallel == 1) then

          ! number of obs primary to this MPI rank
          nobs_loc = count( itabg_m( pf%wgts(:)%imglobe )%irank == myrank )

          if (allocated(pf%iobs_loc)) deallocate(pf%iobs_loc)
          allocate(pf%iobs_loc(nobs_loc))

          ! allocate temporary arrays to store information local to this MPI rank
          allocate(lons_tmp(nobs_loc))
          allocate(lats_tmp(nobs_loc))
          allocate(tims_tmp(nobs_loc))
          allocate(dtim_tmp(nobs_loc))
          allocate(wgts_tmp(nobs_loc))

          if (allocated(pf%altitude)   ) allocate(alts_tmp(nobs_loc))
          if (allocated(pf%pressure)   ) allocate(prss_tmp(nobs_loc))
          if (allocated(pf%obs_dimvals)) allocate(dims_tmp(nobs_loc))

          ! copy data that is local to this MPI rank to the temporary arrays
          if (nobs_loc > 0) then
             il = 0
             do i = 1, nobs
                if ( itabg_m( pf%wgts(i)%imglobe )%irank == myrank ) then
                   il = il + 1

                   lons_tmp(il) = pf%longitude  (i)
                   lats_tmp(il) = pf%latitude   (i)
                   wgts_tmp(il) = pf%wgts       (i)
                   dtim_tmp(il) = pf%s1900_times(i)
                   tims_tmp(il) = pf%dtimes     (i)

                   if (allocated(pf%altitude)   ) alts_tmp(il) = pf%altitude   (i)
                   if (allocated(pf%pressure)   ) prss_tmp(il) = pf%pressure   (i)
                   if (allocated(pf%obs_dimvals)) dims_tmp(il) = pf%obs_dimvals(i)

                   pf%iobs_loc(il) = i
                endif
             enddo
          endif

          ! replace the global arrays with the information only needed on this rank
          call move_alloc(lons_tmp, pf%longitude)
          call move_alloc(lats_tmp, pf%latitude)
          call move_alloc(dtim_tmp, pf%s1900_times)
          call move_alloc(tims_tmp, pf%dtimes)
          call move_alloc(wgts_tmp, pf%wgts)

          if (allocated(alts_tmp)) call move_alloc(alts_tmp, pf%altitude)
          if (allocated(prss_tmp)) call move_alloc(prss_tmp, pf%pressure)
          if (allocated(dims_tmp)) call move_alloc(dims_tmp, pf%obs_dimvals)

       endif

       ! reset obs counter and skip to current time if the model start does not
       ! correspond to the first observation in the input file

       pf%iobs_current = 1

       if (nobs_loc > 0) then

          ! if observations are more then 10 minutes before simulation time
          if (pf%s1900_times(1) < s1900_sim - 600.) then
             do i = 1, nobs_loc
                if (pf%s1900_times(i) > s1900_sim - 600.) exit
             enddo
             pf%iobs_current = i
          endif

       endif

       ! Allocate local space for point output observations. Fortran allows zero-size
       ! arrays if no points exist on a rank

       if (allocated(pf%obs)) deallocate(pf%obs)
       allocate(pf%obs(naddsc,nobs_loc))
       pf%obs(:,:) = rmiss

    enddo

  end subroutine read_point_file



  subroutine fill_point_obs()

    use mem_grid,   only: mza, lpw, zt, dzt, dzim
    use mem_ijtabs, only: itab_m, itabg_m
    use mem_basic,  only: press, rho
    use mem_addsc,  only: addsc
    use misc_coms,  only: s1900_sim, time_bias, naddsc
    use oname_coms, only: nl

    implicit none

    integer                   :: i, inext
    integer                   :: k, kobs, ka
    real                      :: zwt_top, zwt_bot
    integer                   :: im, iw, iaddsc, n
    real                      :: field(3)
    integer                   :: nfiles, nf
    type(point_data), pointer :: pf
    integer,          pointer :: iobs_current, nobs_loc, nobs

    nfiles = nl%point_files

    do nf = 1, nfiles
       pf           => pfiles(nf)
       iobs_current => pfiles(nf)%iobs_current
       nobs_loc     => pfiles(nf)%nobs_loc
       nobs         => pfiles(nf)%nobs

       ! skip if no obs available for current date/time
       if (nobs == 0) cycle

       if (iobs_current <= nobs_loc) then

          if ( pf%s1900_times(iobs_current) < s1900_sim + time_bias ) then

             do i = iobs_current+1, nobs_loc
                if (pf%s1900_times(i) > s1900_sim + time_bias) exit
             enddo
             inext = i - 1

             do i = iobs_current, inext

                ! M point corresponding to this obs
                im = itabg_m( pf%wgts(i)%imglobe )%im_myrank

                if (pf%type == 'Z') then

                   do k = 2, mza
                      if (pf%altitude(i) < zt(k)) exit
                   enddo

                   ! point is between zt(k-1) and zt(k)
                   kobs = k - 1
                   zwt_top = (pf%altitude(i) - zt(kobs)) * dzim(kobs)
                   zwt_bot = 1.0 - zwt_top

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

                      pf%obs(iaddsc,i) = sum( field(1:3) * pf%wgts(i)%wt(1:3) )

                   enddo

                elseif (pf%type == 'S') then

                   do iaddsc = 1, naddsc

                      ! interpolate from 3 surrounding W cells
                      do n = 1, 3
                         iw = itab_m(im)%iw(n)
                         field(n) = addsc(iaddsc)%sclp(lpw(iw),iw)
                      enddo

                      pf%obs(iaddsc,i) = sum( field(1:3) * pf%wgts(i)%wt(1:3) )

                   enddo

                elseif (pf%type == 'P') then

                   do iaddsc = 1, naddsc

                      ! interpolate from 3 surrounding W cells
                      do n = 1, 3
                         iw = itab_m(im)%iw(n)
                         if (pf%pressure(i) >= press(lpw(iw),iw)) then
                            field(n) = addsc(iaddsc)%sclp(lpw(iw),iw)
                         elseif (pf%pressure(i) <= press(mza,iw)) then
                            field(n) = addsc(iaddsc)%sclp(mza,iw)
                         else
                            do k = lpw(iw)+1, mza
                               if (pf%pressure(i) >= press(k,iw)) exit
                            enddo
                            ! point is between press(k-1) and press(k)
                            kobs = k - 1
                            zwt_top = (real(press(kobs,iw)) - pf%pressure(i)) / real(press(kobs,iw) - press(kobs+1,iw))
                            zwt_bot = 1.0 - zwt_top
                            field(n) = zwt_top * addsc(iaddsc)%sclp(kobs+1,iw) &
                                     + zwt_bot * addsc(iaddsc)%sclp(kobs  ,iw)
                         endif
                      enddo

                      pf%obs(iaddsc,i) = sum( field(1:3) * pf%wgts(i)%wt(1:3) )

                   enddo

                elseif (pf%type == 'C') then

                   do iaddsc = 1, naddsc

                      ! interpolate from 3 surrounding W cells
                      do n = 1, 3
                         iw = itab_m(im)%iw(n)
                         ka = lpw(iw)

                         field(n) = sum( addsc(iaddsc)%sclp(ka:mza,iw) * real(rho(ka:mza,iw)) * dzt(ka:mza) ) &
                                  / sum( real(rho(ka:mza,iw)) * dzt(ka:mza) )
                      enddo

                      pf%obs(iaddsc,i) = sum( field(1:3) * pf%wgts(i)%wt(1:3) )

                   enddo

                else

                   write(*,*) "Illegal point file type in fill_point_obs, skipping."
                   exit

                endif

             enddo

             iobs_current = inext + 1

          endif

       endif

    enddo

  end subroutine fill_point_obs



  subroutine write_point_file

      use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
      use misc_coms,  only: iclobber, naddsc
      use max_dims,   only: pathlen
      use oname_coms, only: nl

      implicit none

      character(pathlen)        :: filename
      integer                   :: ndims, idims(2), i
      integer,      allocatable :: itrc_number(:)
      integer                   :: nfiles, nf
      type(point_data), pointer :: pf
      integer,          pointer :: iobs_current, iobs_loc(:), nobs_loc, nobs
      character(30)             :: name = 'tracer'
      character(30)             :: dimname

      nfiles = nl%point_files

      do nf = 1, nfiles
         pf           => pfiles(nf)
         iobs_current => pfiles(nf)%iobs_current
         iobs_loc     => pfiles(nf)%iobs_loc
         nobs_loc     => pfiles(nf)%nobs_loc
         nobs         => pfiles(nf)%nobs

         ! Skip if no obs for current date/time
         if (nobs == 0) cycle

         filename = trim(nl%point_obsout_headers(nf)) // date // trim(nl%point_obsout_suffixes(nf))

         call shdf5_open(filename, 'W', iclobber)

         ! If there is a non-time coordinate variable, write it here

         ndims    = 1
         idims(1) = nobs_loc

         if (len_trim(pf%obs_dimname) > 0 .and. allocated(pf%obs_dimvals)) then

            call shdf5_orec( ndims, idims, pf%obs_dimname, dvar1=pf%obs_dimvals, &
                             isdim     = .true.,          &
                             long_name = pf%obs_longname, &
                             gpoints   = iobs_loc,        &
                             nglobe    = nobs             )
         endif

         ! Write tracer # as a coordinate variable (for A. Schuh)

         if (naddsc > 0) then

            ndims    = 1
            idims(1) = naddsc

            if (allocated(itrc_number)) deallocate(itrc_number)
            allocate(itrc_number(naddsc))
            itrc_number = [ (i, i=1, naddsc) ]

            call shdf5_orec( ndims, idims, "tracer", ivar1=itrc_number, &
                             isdim=.true.,               &
                             long_name = "tracer number" )

         endif

         ! Write lat/lon/time. These are not coordinate variable for point output with NetCDF!

         ndims = 1
         idims(1) = nobs_loc

         call shdf5_orec( ndims, idims, "time", dvar1=pf%dtimes,     &
                          gpoints   = iobs_loc,                      &
                          nglobe    = nobs,                          &
                          isdim     = len_trim(pf%time_dimname) > 0, &
                          dimnames  = [ pf%obs_dimname ],            &
                          long_name = pf%time_longname,              &
                          units     = pf%time_units                  )

         if (len_trim(pf%time_dimname) > 0) then
            dimname = pf%time_dimname
         else
            dimname = pf%obs_dimname
         endif

         call shdf5_orec( ndims, idims, "longitude", rvar1=pf%longitude, &
                          gpoints   = iobs_loc,      &
                          nglobe    = nobs,          &
                          dimnames  = [ dimname ],   &
                          long_name = "longitude",   &
                          units     = "degrees_east" )

         call shdf5_orec( ndims, idims, "latitude", rvar1=pf%latitude, &
                          gpoints   = iobs_loc,       &
                          nglobe    = nobs,           &
                          dimnames  = [ dimname ],    &
                          long_name = "latitude",     &
                          units     = "degrees_north" )

         if (allocated(pf%altitude)) then

            call shdf5_orec( ndims, idims, "altitude", rvar1=pf%altitude, &
                             gpoints   = iobs_loc,       &
                             nglobe    = nobs,           &
                             dimnames  = [ dimname ],    &
                             long_name = "altitude in meters above sea level", &
                             units     = "meters" )

         elseif (allocated(pf%pressure)) then

               call shdf5_orec( ndims, idims, "pressure", rvar1=pf%pressure, &
                                gpoints   = iobs_loc,       &
                                nglobe    = nobs,           &
                                dimnames  = [ dimname ],    &
                                long_name = "atmospheric pressure", &
                                units     = "Pa" )
         endif


         ! Write point observations. For now these are just the tracers for A. Schuh

         if (naddsc > 0) then

            ndims = 2
            idims = [ naddsc, nobs_loc ]

            call shdf5_orec( ndims, idims, "co2_tracers", rvar2=pf%obs, &
                             gpoints = iobs_loc,                        &
                             nglobe = nobs,                             &
                             dimnames = [ name, dimname ],              &
                             rmissing = rmiss,                          &
                             units = "ppm",                             &
                             long_name = "simulated mole fraction of carbon dioxide in dry air" )

         endif

         call shdf5_close()

      enddo

    end subroutine write_point_file


end module point_io
