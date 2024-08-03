subroutine isan_driver(iaction)

  use isan_coms,  only: ihour, idate, imonth, iyear, nfgfiles, ifgfile, &
                        ctotdate_fg, fnames_fg, s1900_fg, pnpr, plats, &
                        o_rho, o_press, o_theta, o_rrw, o_uzonal, o_umerid, &
                        o_ozone, read_analysis_header, z_pbc
  use misc_coms,  only: io6, runtype, s1900_init, s1900_sim, i_o3
  use mem_zonavg, only: zonavg_init
  use mem_grid,   only: mza, mwa
  use mem_nudge,  only: nudflag, nudnxp, o3nudflag
  use hdf5_utils, only: shdf5_exists, shdf5_open, shdf5_close
  use mem_para,   only: olam_mpi_finalize
  use obs_nudge_mod, only: nudge_prep_obs

  implicit none

  integer, intent(in) :: iaction

  integer :: nf
  logical :: exists

  ! Check type of call to isan_driver

  if (iaction == 0) then

! Model run is being started, either for initialization or history restart

! Do inventory of isan file names and times

     call ISAN_file_inv()

! Loop over isan files and search for the one that corresponds to current
! or most recent model time.

     ifgfile = 0
     do nf = 1,nfgfiles
        if (s1900_fg(nf) <= s1900_sim) then
           ifgfile = nf
        endif
     enddo

     if (ifgfile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find isan file for current model time '
        write(io6,*) 'Stopping model '
        stop 'stop: no current isan file'
     elseif (runtype == 'INITIAL' .and. &
             s1900_fg(ifgfile) < s1900_init - 1800.d0) then
        write(io6,*) ' '
        write(io6,*) 'Available isan file is more than 1800 seconds '
        write(io6,*) 'prior to simulation initialization time. '
        write(io6,*) 'Stopping model '
        stop 'stop: initial isan file'
     endif

  elseif (iaction == 1) then

! Processing next isan file (only called with iaction = 1 if nudflag = 1)

     ifgfile = ifgfile + 1

     if (ifgfile > nfgfiles) then
        write(io6,*) ' '
        write(io6,*) 'No future isan file is available for nudging '
        write(io6,*) 'Stopping model '
        stop 'stop: no future isan file'
     endif

  endif

! Get current time

  call date_unmake_big(iyear,imonth,idate,ihour,ctotdate_fg(ifgfile))
  ihour = ihour / 100

! Fill zonavg arrays for current time

  call zonavg_init(idate,imonth,iyear)

! Read header information from gridded pressure-level files for this file time.
! This information includes input data array dimensions.

  write(io6,*)
  write(io6,'(1x,a,i0)')  'Reading ISAN file ifgfile = ', ifgfile
  write(io6,'(1x,a)')     fnames_fg(ifgfile)
  write(io6,'(1x,a,4i6)') ctotdate_fg(ifgfile),iyear,imonth,idate,ihour

  call shdf5_exists(fnames_fg(ifgfile), exists)

  if (.not. exists) then
     write(io6,*) "Error in isan_driver:"
     write(io6,*) "input analysis file cannot be found."
     write(io6,*) "Stopping model run."
     call olam_mpi_finalize()
     stop
  endif

  call shdf5_open(fnames_fg(ifgfile), 'R')

  call read_analysis_header(nosoil=.true.)

! Allocate memory for ISAN processing

  allocate( o_rho   (mza,mwa) )
  allocate( o_press (mza,mwa) )
  allocate( o_theta (mza,mwa) )
  allocate( o_rrw   (mza,mwa) )
  allocate( o_uzonal(mza,mwa) )
  allocate( o_umerid(mza,mwa) )

  if (i_o3 > 0) then
     allocate( o_ozone (mza,mwa) )
  endif

  allocate( z_pbc(mwa) )

! Read gridded pressure-level data, add any ZONAVG fields as necessary, and
! interpolate analysis fields to OLAM grid points

  call pressure_stage()

! Perform iterative hydrostatic balancing on analysis fields

  call vterpp()

  deallocate( z_pbc )

! Copy interpolated and adjusted analysis fields to model arrays if initializing

  if (iaction == 0 .and. runtype == 'INITIAL') then
     call isnstage()
  endif

! Pressure no longer needed

  deallocate( o_press )

! If nudging, prepare observational nudging fields

  if (nudflag > 0) then
     if (nudnxp == 0) then
        call nudge_prep_obs ( iaction )
     else
        call nudge_prep_spec( iaction )
     endif
  endif

  if (o3nudflag > 0) then
     call nudge_prep_o3( iaction )
  endif

! Deallocate ISAN arrays if necessary

  if (allocated( o_rho    )) deallocate( o_rho    )
  if (allocated( o_theta  )) deallocate( o_theta  )
  if (allocated( o_rrw    )) deallocate( o_rrw    )
  if (allocated( o_uzonal )) deallocate( o_uzonal )
  if (allocated( o_umerid )) deallocate( o_umerid )
  if (allocated( o_ozone  )) deallocate( o_ozone  )

! Allocated in read_press_header

  if (allocated( pnpr  )) deallocate( pnpr  )
  if (allocated( plats )) deallocate( plats )

! Close HDF5 file

  call shdf5_close()

end subroutine isan_driver
