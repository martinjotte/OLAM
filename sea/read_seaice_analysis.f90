subroutine read_seaice_analysis(iaction)

  use mem_sea,     only: sea, msea, omsea
  use mem_sfcg,    only: sfcg
  use sea_coms,    only: iseaicefile, iseaiceflg
  use misc_coms,   only: io6, s1900_sim
  use max_dims,    only: pathlen
  use isan_coms,   only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, plats, &
                         inproj, xswlat, xswlon, gdatdx, gdatdy, irev_ns, &
                         read_analysis_header
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use analysis_lib,only: gdtost_ll
  use mem_para,    only: olam_mpi_finalize

  implicit none

  integer, intent(in) :: iaction

  integer            :: nf
  character(pathlen) :: fname
  integer            :: ndims, idims(3)
  integer            :: isea, iwsfc
  real,  allocatable :: ice(:,:)  ! sea ice concentration [0 - 1]
  logical            :: exists

! Nothing to do here if iseaiceflg is not 2

  if (iseaiceflg /= 2) return

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

! Process selected analysis file

  write(io6,'(A)') ' read_seaice: opening ' // trim(fname)

  call shdf5_exists(fname, exists)

  if (.not. exists) then
     write(io6,*) "read_seaice_analysis: Error opening analysis file " // trim(fname)
     write(io6,*) "Stopping run."
     call olam_mpi_finalize()
     stop
  endif

  call shdf5_open(fname, 'R')

  call read_analysis_header(noplevs=.true., nosoil=.true.)

  ! Check if seaice is in the analysis file, and read it

  call shdf5_info('ICEC', ndims, idims)

  if (ndims>0 .and. idims(1)==nprx .and. idims(2)==npry) then

     allocate(ice(nprx,npry))
     call shdf5_irec(ndims, idims, 'ICEC', rvar2=ice)

  else

     write(io6,*) "read_seaice: Analysis file does not contain seaice."
     write(io6,*) "Stopping run."
     call olam_mpi_finalize()
     stop

  endif

  call shdf5_close()

  ! Fill seaice array

  !$omp parallel do private(iwsfc)
  do isea = 2, msea
     iwsfc = isea + omsea

     call gdtost_ll( nprx, npry, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                     sea%seaicef(isea), xswlon, xswlat, gdatdx, gdatdy, inproj, &
                     r2d=ice, plats=plats, irev_ns=irev_ns )

     sea%seaicef(isea) = max( min( sea%seaicef(isea), 1.0 ), 0.0 )
  enddo
  !$omp end parallel do

  if (iaction == 0) then
     sea%seaicep(:) = sea%seaicef(:)
  endif

end subroutine read_seaice_analysis
