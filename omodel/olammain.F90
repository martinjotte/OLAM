program main

  use misc_coms,    only: io6, iparallel
  use mem_para,     only: myrank, mgroupsize, olam_mpi_init, &
                          olam_mpi_finalize, olam_mpi_barrier
  use misc_coms,    only: tmpdir
  use max_dims,     only: pathlen
  use hdf5,         only: h5open_f, h5close_f
  use hdf5_f2f,     only: fh5_close_caches
  use oname_coms,   only: cmdlne_runtype, cmdlne_fields, numcf, nl
  use sys_utils,    only: set_environment_variable
  !$ use omp_lib,   only: omp_get_max_threads

  implicit none

  character(pathlen) :: name_name = 'OLAMIN'
  character(pathlen) :: cargv
  character(30)      :: io6file
  integer :: numarg
  integer :: i
  integer :: nthreads
  integer :: bad
  integer :: hdferr

#if defined(__xlC__) || defined(__ibmxl__)
  character(10)      :: dowrite
#endif

  ! Some versions of parallel HDF5 hang without this set:
  call set_environment_variable("HDF5_USE_FILE_LOCKING", "FALSE")

  ! OpenMPI up to version 4.1.4 has a bad bug with parallel IO
  ! https://www.hdfgroup.org/2022/11/workarounds-for-openmpi-bug-exposed-by-make-check-in-hdf5-1-13-3
  call set_environment_variable("OMPI_MCA_io_ompio_priority", "1")
  call set_environment_variable("OMPI_MCA_fbtl_posix_read_datasieving", "0")

! If run is sequential, default choice is to set io6 to standard output unit 6.

  io6 = 6

! Determine if this run is parallel, and determine myrank and mgroupsize

  call olam_mpi_init()

  if (mgroupsize > 1) then
     iparallel = 1
  else
     iparallel = 0
  endif

! Initialize HDF5 library

  call h5open_f(hdferr)

! Parse the command line arguments

  numarg = command_argument_count()
  if (myrank == 0) write(*,*) 'numarg:', numarg

  allocate(cmdlne_fields(numarg/2 + 2))

  i = 1
  bad = 0
  do while (i <= numarg)
     call get_command_argument(i,cargv)
     ! write(io6,*) 'args: ',i,cargv

     if (cargv(1:1) == '-') then

        if (cargv(2:2) == 'f') then

           call get_command_argument(i+1,name_name)
           if (len_trim(name_name) < 1) bad = bad + 1
           i = i + 2

        elseif (cargv(2:2) == 'r') then

           call get_command_argument(i+1,cmdlne_runtype)
           if (len_trim(cmdlne_runtype) < 1) bad = bad + 1
           i = i + 2

        elseif (cargv(2:2) == 'z') then

           if (numcf < size(cmdlne_fields)) then

              numcf = numcf + 1

              call get_command_argument(i+1,cmdlne_fields(numcf))
              if ( len_trim( cmdlne_fields(numcf) ) < 1  .or. &
                   scan( cmdlne_fields(numcf), '=') == 0 ) then
                 bad = bad + 1
                 numcf = numcf - 1
              endif

           else

              if (myrank == 0) write(*,*) "OLAM too many 'z' arguments"

           endif
           i = i + 2

        else
           i = i + 1
        endif
     else
        i = i + 1
     endif
  enddo

  if (bad > 0) then
     if (myrank == 0) write(*,*) 'OLAM usage: ''exec name'' '
     if (myrank == 0) write(*,*) '  [-f ''Namelist file''] '
     if (iparallel == 1) call olam_mpi_barrier()
     if (myrank == 0) then
        stop 'Bad command line arguments'
     else
        stop
     endif
  endif

! Read Fortran namelist

  if (myrank == 0) write(*,*) 'OLAM input namelist file: ',trim(name_name)
  if (myrank == 0) write(*,'(/,a)') 'olammain calling namelist read'
  call read_nl(name_name)

! If run is parallel, default choice is to attach output unit io6 to separate files

  if (iparallel == 1 .and. myrank > 0) then

     io6 = 20

     if (nl%save_node_logs) then

        write (io6file,'(i10)') myrank
        io6file = 'o.io6_r'//trim(adjustl(io6file))

        ! Output file is replaced if it exists
        open(unit=io6, file=io6file, status='replace', form='formatted', action='write')

#ifdef __PGI
        ! Prevent Portland Group compiler from buffering io6
        call setvbuf3f(io6,1,1024)
#endif

     else

        io6file = '/dev/null'

#if defined(__xlC__) || defined(__ibmxl__)
        ! Work around an issue with IBM compiler and MPI
        inquire(file=io6file, write=dowrite, number=i)
        if (dowrite == 'YES' .and. i > 0) close(i)
#endif

        open(unit=io6, file=io6file, action='write')

     endif
  endif

  nthreads = 1
  !$ nthreads = omp_get_max_threads()

  write(io6,'(/,a,i6  )') ' myrank     = ', myrank
  write(io6,'(  a,i6  )') ' mgroupsize = ', mgroupsize
  write(io6,'(  a,i6  )') ' iparallel  = ', iparallel
  write(io6,'(  a,i6,/)') ' nthreads   = ', nthreads

! Get unix tmp directory

  call get_environment_variable("TMPDIR", cargv, length=i)
  if ( i > 0 ) then
     if (cargv(i:i) == '/') then
        tmpdir = cargv(1:i-1)
     else
        tmpdir = cargv
     endif
  else
     tmpdir = "/tmp"
  endif

  write(io6,*) "Using ", trim(tmpdir), " as the temporary directory."

! Initialize, execute, and end olam run

  call olam_run(name_name)

! stop HDF5 library

  call fh5_close_caches(hdferr)
  call h5close_f       (hdferr)

! If this run is parallel, finalize MPI and close io6 file

  call olam_mpi_finalize()

  write(io6,'(A)')
  write(io6,'(A)') "olam_end"

  if (iparallel == 1) then
     close(io6)
  endif

end program main
