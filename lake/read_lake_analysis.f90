subroutine read_lake_analysis()

  use mem_lake,    only: lake, mlake, omlake
  use mem_sfcg,    only: sfcg
  use sea_coms,    only: isstflg
  use misc_coms,   only: io6, s1900_sim
  use max_dims,    only: pathlen
  use consts_coms, only: t00, cliq, alli
  use isan_coms,   only: nfgfiles, s1900_fg, fnames_fg, nprx, npry, plats, &
                         inproj, xswlat, xswlon, gdatdx, gdatdy, irev_ns, &
                         read_analysis_header
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use analysis_lib,only: gdtost_ll
  use therm_lib,   only: rhovsl

  implicit none

  integer            :: nf, ilakefile
  character(pathlen) :: fname
  character(16)      :: ext
  integer            :: ndims, idims(3)
  real               :: laketemp
  integer            :: nio, njo, ilake, iwsfc
  real, allocatable  :: wtemp(:,:)  ! lake temperature [K]
  logical            :: exists

! Nothing to do here if isstflg is not 2

  if (isstflg /= 2) return

  ! Get the listing of analysis files if the inventory was
  ! never called
  if (nfgfiles < 0) call isan_file_inv()

  ! Loop over analysis files and search for the one that corresponds
  ! to the model start time

  ilakefile = 0
  do nf = 1, nfgfiles
     if (s1900_fg(nf) <= s1900_sim) then
        ilakefile = nf
     endif
  enddo

  if (ilakefile < 1) then
     write(io6,*) ' '
     write(io6,*) 'Unable to find analysis for lake temperature'
     write(io6,*) 'Stopping model.'
     stop 'stop: no current analysis file for sst/seaice file'
  endif

  fname = fnames_fg(ilakefile)

! Process selected analysis file

  write(io6,'(A)') ' read_lake: opening ' // trim(fname)

  call shdf5_exists(fname, exists)

  if (.not. exists) then
     write(io6,*) "read_lake_analysis: Error opening analysis file " // trim(fname)
     write(io6,*) "Using default lake initialization instead."
     return
  endif

  call shdf5_open(fname, 'R')

  call read_analysis_header(noplevs=.true., nosoil=.true.)

  ! Check if lake temperature is in the analysis file, and read it

  call shdf5_info('SST', ndims, idims)

  if (ndims>0 .and. idims(1)==nprx .and. idims(2)==npry) then

     allocate(wtemp(nprx,npry))
     call shdf5_irec(ndims, idims, 'SST', rvar2=wtemp)

  else

     write(io6,*) "read_lake: Analysis file does not contain lake temperature."
     write(io6,*) "Stopping run."
     stop

  endif

  call shdf5_close()

  ! Fill lake energy

  !$omp parallel do private(iwsfc,laketemp)
  do ilake = 2, mlake
     iwsfc = ilake + omlake

     call gdtost_ll( nprx, npry, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                     laketemp, xswlon, xswlat, gdatdx, gdatdy, inproj, &
                     r2d=wtemp, plats=plats, irev_ns=irev_ns )

     lake%lake_energy(ilake) = (laketemp - t00) * cliq + alli

  enddo
  !$omp end parallel do

end subroutine read_lake_analysis
