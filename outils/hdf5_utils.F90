Module hdf5_utils

!!!!!!!!!!!!!!!! temp
  use hdf5_f2f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use max_dims, only: pathlen

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

  use hdf5, only: FORTRAN_REAL4_TYPE => H5T_IEEE_F32LE, &
                  FORTRAN_REAL8_TYPE => H5T_IEEE_F64LE, &
                  FORTRAN_INT1_TYPE  => H5T_STD_I8LE,   &
                  FORTRAN_INT2_TYPE  => H5T_STD_I16LE,  &
                  FORTRAN_INT4_TYPE  => H5T_STD_I32LE,  &
                  HID_T

  implicit none (external, type)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
  logical, parameter :: has_phdf5 = .true.
#else
  logical, parameter :: has_phdf5 = .false.
#endif

  logical :: do_phdf5
  logical :: do_hdf5
  logical :: do_comm
  logical :: indepio

  character(pathlen) :: fname = ' '
  character(1)       :: fmode = ' '

  private
  public :: shdf5_open, shdf5_info, shdf5_orec, shdf5_irec, shdf5_close, &
            shdf5_io, shdf5_orec_ll, shdf5_write_global_attribute, &
            shdf5_exists, FORTRAN_REAL4_TYPE, FORTRAN_REAL8_TYPE, &
            FORTRAN_INT1_TYPE, FORTRAN_INT2_TYPE, FORTRAN_INT4_TYPE

Contains

!===============================================================================

subroutine shdf5_exists(locfn, exists, serial)

#ifdef OLAM_MPI
  use mpi_f08,    only: MPI_Bcast, MPI_LOGICAL, MPI_Allreduce, MPI_IN_PLACE, &
                        MPI_LOGICAL, MPI_LAND
  use mem_para,   only: MPI_COMM_OLAM
#endif
  use misc_coms,  only: iparallel
  use mem_para,   only: myrank
  use oname_coms, only: nl
  use hdf5_f2f,   only: fh5f_exists
  import,         only: has_phdf5

  implicit none (external, type)

  character(*),      intent(in)  :: locfn
  logical,           intent(out) :: exists
  logical, optional, intent(in)  :: serial
  logical                        :: dophdf5, dohdf5
  logical                        :: docomm, doser

  ! Are we using parallel HDF5
  dophdf5 = (iparallel==1) .and. has_phdf5 .and. (.not. nl%disable_phdf5_reads )

  ! Does this rank directly call HDF5 routines
  dohdf5 = (myrank==0) .or. (iparallel==1 .and. nl%allranks_read_hdf5) .or. dophdf5

  ! Do we need MPI routines to send/receive data to/from node 0
  docomm = (iparallel==1) .and. (.not. dophdf5) .and. (.not. nl%allranks_read_hdf5 )

  ! Is each process doing I/O independent from the other processes
  doser = iparallel /= 1

  if ( (iparallel == 1) .and. present(serial) ) then
     dophdf5 = dophdf5 .and. (.not. serial)
     docomm  = docomm  .and. (.not. serial)
     dohdf5  = dohdf5  .or.         serial
     doser   = serial
  endif

  if (dohdf5) call fh5f_exists(locfn, exists, pario=dophdf5)

#ifdef OLAM_MPI
  if (docomm) then
     call MPI_Bcast(exists, 1, MPI_LOGICAL, 0, MPI_COMM_OLAM)
  elseif (.not. doser) then
     call MPI_Allreduce(MPI_IN_PLACE, exists, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_OLAM)
  endif
#endif

end subroutine shdf5_exists

!===============================================================================

subroutine shdf5_open(locfn, access, idelete, serial)

#ifdef OLAM_MPI
  use mpi_f08,    only: MPI_Bcast, MPI_LOGICAL, MPI_Allreduce, MPI_IN_PLACE, MPI_LOR
  use mem_para,   only: MPI_COMM_OLAM
#endif
  use misc_coms,  only: iparallel, iclobber, io6
  use mem_para,   only: olam_mpi_finalize, myrank
  use oname_coms, only: nl
  use hdf5_f2f,   only: fh5f_open, fh5f_create
  import,         only: has_phdf5, do_phdf5, do_hdf5, do_comm, indepio, fname, fmode

  implicit none (external, type)

  character(*),      intent(in) :: locfn   ! file name
  character(1),      intent(in) :: access  ! File access ('R' or 'W')
                                           ! read, write, or append
  integer, optional, intent(in) :: idelete ! If W, overwrite file if exists?
                                           !    1=yes, 0=no
  logical, optional, intent(in) :: serial  ! Disable collective/parallel features
                                           ! (all nodes independent or each other)

  integer                       :: hdferr  ! Error flag
  integer                       :: iaccess ! int access flag
  logical                       :: exists  ! File existence
  logical                       :: ierror  ! Error opening file
  logical                       :: idel    ! Delete existing file flag

  if (access /= 'R' .and. access /= 'W') then
     write(io6,*) 'Error in shdf5_open:'
     write(io6,*) '   Access must be either R or W.'
     write(io6,*) 'Returning to model without opening ' // trim(locfn)
     return
  endif

  if (access == 'R') then  ! When reading, all ranks may access the file.

     ! Are we using parallel HDF5
     do_phdf5 = (iparallel==1) .and. has_phdf5 .and. (.not. nl%disable_phdf5_reads )

     ! Does this rank directly call HDF5 routines
     do_hdf5 = (myrank==0) .or. nl%allranks_read_hdf5 .or. do_phdf5

     ! Do we need MPI routines to send/receive data to/from node 0
     do_comm = (iparallel==1) .and. (.not. do_phdf5) .and. (.not. nl%allranks_read_hdf5 )

     ! Are all nodes doing serial reads/writes on same file
     indepio = (iparallel==1) .and. (.not. do_phdf5) .and. nl%allranks_read_hdf5

  else  ! When writing, only one rank may access the file at a time if not using PHDF5

     ! Are we using parallel HDF5
     do_phdf5 = (iparallel==1) .and. has_phdf5 .and. (.not. nl%disable_phdf5_writes )

     ! Does this rank directly call HDF5 routines
     do_hdf5 = (myrank==0) .or. do_phdf5

     ! Do we need MPI routines to send/receive data to/from node 0
     do_comm = (iparallel==1) .and. (.not. do_phdf5)

     ! Are all nodes doing serial reads/writes on same file
     indepio = .false.

  endif

  if (present(serial)) then
     ! If nodes read different files, or only some nodes open file
     do_phdf5 = do_phdf5 .and. (.not. serial)
     do_hdf5  = do_hdf5  .or.         serial
     do_comm  = do_comm  .and. (.not. serial)
     indepio  = indepio  .and. (.not. serial)
  endif

  ! Check if a HDF5 file is already open. Currently we can only have one file
  ! open at a time.

  ierror = .false.

  if (do_hdf5 .and. len_trim(fname) > 0) then
     write(io6,*) 'Error in shdf5_open:'
     write(io6,*) '   File ', trim(locfn)
     write(io6,*) '   is being opened while another file:'
     write(io6,*) '   ', trim(fname)
     write(io6,*) '   is still open. Currently only one HDF5 file can be open at a time.'
     ierror = .true.
  endif

  ! Open/create file if this process rank is performing I/O

  if (do_hdf5 .and. .not. ierror) then

     hdferr = 0

     ! Create a new file or open an existing RAMS file.

     if (access == 'R' .or. access == 'A') then

        if (access == 'R') iaccess = 1
        if (access == 'A') iaccess = 2

        call fh5f_open(locfn, iaccess, hdferr, pario=do_phdf5)

        if (hdferr < 0) then
           write(io6,*) 'shdf5_open:'
           write(io6,*) '   Error opening hdf5 file - error -', hdferr
           write(io6,*) '   Filename: ',trim(locfn)
           write(io6,*) 'shdf5_open: open error'
           ierror = .true.
        endif

     elseif (access == 'W') then

        if (present(idelete)) then
           idel = (idelete  == 1)
        else
           idel = (iclobber == 1)
        endif

        if (idel) then
           iaccess = 1
        else
           iaccess = 2
        endif

        ! Check for existence of file.
        inquire(file=locfn, exist=exists)

        if (.not. exists) then

           call fh5f_create(locfn, iaccess, hdferr, pario=do_phdf5)

        else

           if (idel) then
              write(io6,*) 'Existing file ' // trim(locfn) // ' is being deleted on write'
              call fh5f_create(locfn, iaccess, hdferr, pario=do_phdf5)
              if (hdferr < 0) then
                 write(io6,*) 'shdf5_open:'
                 write(io6,*) '   Error creating hdf5 file - error -', hdferr
                 write(io6,*) '   Filename: ',trim(locfn)
                 write(io6,*) 'shdf5_open: create error'
                 ierror = .true.
              endif
           else
              write(io6,*) 'In shdf5_open:'
              write(io6,*) '   Attempt to open an existing file for writing,'
              write(io6,*) '   but overwrite is disabled.'
              write(io6,*) '   Filename: ', trim(locfn)
              write(io6,*) 'Delete existing file or change ICLOBBER in OLAMIN namelist'
              ierror = .true.
           endif

        endif
     endif
  endif

#ifdef OLAM_MPI
  if (do_comm) then
     call MPI_Bcast(ierror, 1, MPI_LOGICAL, 0, MPI_COMM_OLAM, hdferr)
  elseif (do_phdf5 .or. indepio) then
     call MPI_Allreduce(MPI_IN_PLACE, ierror, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_OLAM, hdferr)
  endif
#endif

  if (ierror) then
     write(io6,*)
     write(io6,*) 'Stopping model in shdf5_open due to error opening ' // trim(locfn)
     call olam_mpi_finalize()
     stop
  endif

  if (do_hdf5) then
     fname = locfn
     fmode = access
  endif

end subroutine shdf5_open

!===============================================================================

subroutine shdf5_info(dsetname, ndims, dims, dimname, attached_dimnames, &
                      units, chunk_dims)

#ifdef OLAM_MPI
  use mpi_f08,  only: MPI_Bcast, MPI_INTEGER, MPI_CHARACTER
  use mem_para, only: MPI_COMM_OLAM, olam_mpi_barrier
#endif
  use hdf5_f2f, only: fh5_open_dataset, fh5_get_info, fh5_close_dataset, &
                      fh5f_query_dimname, fh5_get_attached_scales, &
                      fh5f_read_attribute, fh5_get_chunk_dims
  import,       only: do_hdf5, do_comm, indepio

  implicit none (external, type)

  character(*),        intent(in)    :: dsetname ! Dataset name
  integer, contiguous, intent(inout) :: dims(:)
  integer,             intent(inout) :: ndims    ! Dataset rank (in file)

  ! Optional arrays to read common NetCDF convention attributes or dimension information
  character(*),        intent(inout), optional :: units
  character(*),        intent(inout), optional :: dimname
  character(*),        intent(inout), optional :: attached_dimnames(:)
  integer, contiguous, intent(inout), optional :: chunk_dims(:)

  integer :: hdferr

  ndims = -1
  dims  =  0

  if (present(attached_dimnames)) attached_dimnames(:) = ' '
  if (present(dimname))           dimname              = ' '
  if (present(units))             units                = ' '
  if (present(chunk_dims))        chunk_dims(:)        = 0

  ! Open the dataset.

  if (do_hdf5) call fh5_open_dataset(dsetname, hdferr)
#ifdef OLAM_MPI
  if (do_comm) call MPI_Bcast(hdferr, 1, MPI_INTEGER, 0, MPI_COMM_OLAM)
#endif

  ! Return if there was an error (dataset not in file)

  if (hdferr < 0) return

  ! Get dimension information for the dataset

  if (do_hdf5) call fh5_get_info(dsetname, ndims, dims)
#ifdef OLAM_MPI
  if (do_comm) call MPI_Bcast(ndims, 1,         MPI_INTEGER, 0, MPI_COMM_OLAM)
  if (do_comm) call MPI_Bcast(dims, size(dims), MPI_INTEGER, 0, MPI_COMM_OLAM)
#endif

  ! Return if there was an error

  if (ndims < 0) then
     if (do_hdf5) call fh5_close_dataset(hdferr)
     return
  endif

  ! Is this variable a dimension?

  if (present(dimname)) then
     if (do_hdf5) call fh5f_query_dimname(dimname)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(dimname, len(dimname), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  ! Does the variable have dimension scales?

  if (present(attached_dimnames)) then
     if (do_hdf5) call fh5_get_attached_scales(attached_dimnames)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(attached_dimnames, len(attached_dimnames)*size(attached_dimnames), &
                                  MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  ! Read any attributes

  if (present(units)) then
     if (do_hdf5) call fh5f_read_attribute("units", cvalue=units)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(units, len(units), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(chunk_dims)) then
     if (do_hdf5) call fh5_get_chunk_dims(ndims,chunk_dims)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(chunk_dims, size(chunk_dims), MPI_INTEGER, 0, &
                                 MPI_COMM_OLAM)
#endif
  endif

  ! Close the dataset

  if (do_hdf5) call fh5_close_dataset(hdferr)

  ! If all nodes are independently reading the same file outside of parallel HDF5,
  ! wait here for all nodes to finish

#ifdef OLAM_MPI
  if (indepio) call olam_mpi_barrier()
#endif

end subroutine shdf5_info

!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars,    &
                                          bvar1,ivar1,rvar1,dvar1,lvar1,    &
                                          bvar2,ivar2,rvar2,dvar2,lvar2,    &
                                          bvar3,ivar3,rvar3,dvar3,lvar3,    &
                                          bvar4,ivar4,rvar4,dvar4,          &
                                          nglobe, lpoints, gpoints,         &
                                          units, long_name, positive,       &
                                          imissing, rmissing, dmissing,     &
                                          isdim, dimnames, standard_name,   &
                                          cell_methods, dims_chunk,         &
                                          storage_type                      )

  use oname_coms, only: nl
  use misc_coms,  only: io6
  use hdf5_f2f,   only: fh5_prepare_write, fh5_write, fh5_close_write, &
                        fh5f_attach_dims, fh5f_create_dim, fh5f_write_attribute, &
                        fh5_close_dataset
  import,         only: i1, r8, do_hdf5, do_comm, fname, HID_T, FORTRAN_INT1_TYPE, &
                        FORTRAN_INT4_TYPE, FORTRAN_REAL4_TYPE, FORTRAN_REAL8_TYPE

  implicit none (external, type)

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional :: bvars
  integer,      intent(in), optional :: ivars
  real,         intent(in), optional :: rvars
  real(r8),     intent(in), optional :: dvars
  logical,      intent(in), optional :: lvars

  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional             :: nglobe
  integer,      intent(in), optional, contiguous :: lpoints(:), gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

! Compression/chunking options
  integer,      intent(in), optional :: dims_chunk(ndims) ! Compression dimensions.

! Type of variable
  integer(HID_T), intent(in), optional :: storage_type

! Local variables
  integer                  :: hdferr, n
  integer(HID_T)           :: stype
  integer(i1)              :: lbufs
  integer(i1), allocatable :: lbuf1(:), lbuf2(:,:), lbuf3(:,:,:)

  ! Maybe we can get around specifying all possible ranks and types by using
  ! assumed-rank and/or assumed-type arrays and C interoperability once
  ! most compilers support this

  n = count( [present(bvars), present(bvar1), present(bvar2), present(bvar3), present(bvar4), &
              present(ivars), present(ivar1), present(ivar2), present(ivar3), present(ivar4), &
              present(rvars), present(rvar1), present(rvar2), present(rvar3), present(rvar4), &
              present(dvars), present(dvar1), present(dvar2), present(dvar3), present(dvar4), &
              present(lvars), present(lvar1), present(lvar2), present(lvar3)] )

  if (n /= 1) then
     write(io6,*) 'Error in shdf5_orec writing ' // trim(dsetname)
     write(io6,*) 'Only one variable type should be specified for HDF5 output.'
     return
  endif

  ! check arguments for parallel output

  if (present(gpoints)) then

     if (.not. present(nglobe)) then
        write(io6,*) 'Error in shdf5_orec reading ' // trim(dsetname)
        write(io6,*) 'Both gpoints and nglobe must be specified for parallel output.'
        return
     endif

     if (present(lpoints)) then
        if (size(gpoints) /= size(lpoints)) then
           write(io6,*) 'Error in shdf5_orec reading ' // trim(dsetname)
           write(io6,*) 'gpoints and lpoints must have same number of points for parallel output.'
           return
        endif
     endif

  endif

  ! Check dimensions

! if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
!    print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims), trim(dsetname)
!    stop    'shdf5_orec: bad dims'
! endif

  ! Do we want to change the storage type

  if (present(storage_type)) then
     stype = storage_type
  elseif (present(bvars) .or. present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
     stype = FORTRAN_INT1_TYPE
  elseif (present(ivars) .or. present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
     stype = FORTRAN_INT4_TYPE
  elseif (present(rvars) .or. present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
     stype = FORTRAN_REAL4_TYPE
  elseif (present(dvars) .or. present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
     stype = FORTRAN_REAL8_TYPE
  elseif (present(lvars) .or. present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
     stype = FORTRAN_INT1_TYPE
  endif

  ! Convert logicals to 1-byte integers, since HDF5 does not handle Fortran logical/bools

  if (present(lvars)) then
     lbufs = merge(-1_i1, 0_i1, lvars)
  elseif (present(lvar1)) then
     allocate(lbuf1, source=merge(-1_i1, 0_i1, lvar1))
  elseif (present(lvar2)) then
     allocate(lbuf2, source=merge(-1_i1, 0_i1, lvar2))
  elseif (present(lvar3)) then
     allocate(lbuf3, source=merge(-1_i1, 0_i1, lvar3))
  endif

! write(io6,'(A,8(1x,I0))') " Writing: "//trim(dsetname), dims(1:ndims)

  if (do_hdf5) then

     ! Prepare memory and options for the write

     call fh5_prepare_write(ndims, dims, dsetname, stype, hdferr, &
                            icompress=nl%icompress, &
                            mcoords=lpoints, fcoords=gpoints, ifsize=nglobe, &
                            dims_chunk=dims_chunk)

     if (hdferr < 0) then
        print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
        stop
     endif

     ! Write the dataset.

     if     (present(ivars)) then ; call fh5_write(ivars, hdferr)
     elseif (present(rvars)) then ; call fh5_write(rvars, hdferr)
     elseif (present(dvars)) then ; call fh5_write(dvars, hdferr)
     elseif (present(lvars)) then ; call fh5_write(lbufs, hdferr)
     elseif (present(bvars)) then ; call fh5_write(bvars, hdferr)

     elseif (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
     elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
     elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
     elseif (present(lvar1)) then ; call fh5_write(lbuf1, hdferr)
     elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

     elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
     elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
     elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
     elseif (present(lvar2)) then ; call fh5_write(lbuf2, hdferr)
     elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

     elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
     elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
     elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
     elseif (present(lvar3)) then ; call fh5_write(lbuf3, hdferr)
     elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)

     elseif (present(ivar4)) then ; call fh5_write(ivar4, hdferr)
     elseif (present(rvar4)) then ; call fh5_write(rvar4, hdferr)
     elseif (present(dvar4)) then ; call fh5_write(dvar4, hdferr)
     elseif (present(bvar4)) then ; call fh5_write(bvar4, hdferr)

     else
        print*, 'Incorrect or missing data field argument in shdf5_orec'
        stop    'shdf5_orec: bad data field'
     endif

     if (hdferr /= 0) then
        print*, 'In shdf5_orec: hdf5 write error =', hdferr
        stop    'shdf5_orec: hdf5 write error'
     endif

     ! Close the write buffers but keep the data set open

     call fh5_close_write(hdferr, no_dsetid=.true.)

     ! Indicate if this variable is a dimension

     if (present(isdim)) then
        if (isdim) call fh5f_create_dim(dsetname, hdferr)
     endif

     ! Link each variable to its dimensions

     if ((present(dimnames)) .and. (.not. present(isdim))) then
        if (size(dimnames) >= ndims) call fh5f_attach_dims(ndims, dimnames, hdferr)
     endif

     ! Write any dataset attributes to match common NetCDF conventions

     if (present(units)) then
        if (len_trim(units) > 0) call fh5f_write_attribute("units", cvalue=units)
     endif

     if (present(long_name)) then
        if (len_trim(long_name) > 0) call fh5f_write_attribute("long_name", cvalue=long_name)
     endif

     if (present(standard_name)) then
        if (len_trim(standard_name) > 0) call fh5f_write_attribute("standard_name", cvalue=standard_name)
     endif

     if (present(cell_methods)) then
        if (len_trim(cell_methods) > 0) call fh5f_write_attribute("cell_methods", cvalue=cell_methods)
     endif

     if (present(positive)) then
        if (len_trim(positive) > 0) call fh5f_write_attribute("positive", cvalue=positive)
     endif

     if (present(imissing)) then
        call fh5f_write_attribute("missing_value", ivalue=imissing)
     endif

     if (present(rmissing)) then
        call fh5f_write_attribute("missing_value", rvalue=rmissing)
     endif

     if (present(dmissing)) then
        call fh5f_write_attribute("missing_value", dvalue=dmissing)
     endif

     call fh5_close_dataset(hdferr)

  endif

#ifdef OLAM_MPI
  if ( do_comm .and. present(gpoints) ) then

     if (nl%allranks_write_hdf5) then
        call shdf5_orec_no_phdf5_allranks()
     else
        call shdf5_orec_no_phdf5_rank0()
     endif

  endif

  Contains

  subroutine shdf5_orec_no_phdf5_rank0()

    use mem_para, only: mgroupsize, myrank, MPI_COMM_OLAM

    use mpi_f08,  only: MPI_Igather, MPI_INTEGER, MPI_Send, MPI_Isend, MPI_REAL, &
                        MPI_REAL8, MPI_INTEGER1, MPI_Waitall, MPI_STATUSES_IGNORE, &
                        MPI_Recv, MPI_STATUS_IGNORE, MPI_Request, MPI_Wait

    use hdf5_f2f, only: fh5_prepare_write, fh5_write, fh5_prepare_write, fh5_close_write

    import,       only: dsetname, ndims, dims, nglobe, lpoints, gpoints, stype, &
                        bvar1, bvar2, bvar3, bvar4, ivar1, ivar2, ivar3, ivar4, &
                        rvar1, rvar2, rvar3, rvar4, dvar1, dvar2, dvar3, dvar4, &
                        lvar1, lvar2, lvar3, r8, i1, fh5_prepare_write, &
                        lbuf1, lbuf2, lbuf3

    implicit none (external, type)

    integer,     allocatable :: ibuff(:), points(:), nus(:)
    real,        allocatable :: rbuff(:)
    real(r8),    allocatable :: dbuff(:)
    integer(i1), allocatable :: bbuff(:)
    integer                  :: hdferr, maxbuff, locbuff
    integer                  :: nu, ier, base, n, i, is
    integer,       parameter :: itag1 = 2098
    integer,       parameter :: itag2 = 2099
    type(MPI_Request)        :: ireqs(2)
    integer                  :: dimsn(5)

    ! Collect on rank 0 the number of cells output by each process

    nu = size(gpoints)

    if (myrank == 0) then
       allocate(nus(mgroupsize))
    else
       allocate(nus(1))
    endif

    call MPI_Igather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_OLAM, ireqs(1), ier)

    ! Send to rank 0 the list of points and data to be written to disk

    if (myrank > 0 .and. nu > 0) then

       locbuff = nu
       is      = 1
       do n = ndims-1, 1, -1
          locbuff = locbuff * dims(n)
          is      = is      * dims(n)
       enddo

       call MPI_Isend(gpoints, nu, MPI_INTEGER, 0, itag1, MPI_COMM_OLAM, ireqs(2), ier)

       if (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then

          allocate(ibuff(locbuff))
          do n = 1, nu
             if (present(lpoints)) then
                i = lpoints(n)
             else
                i = n
             endif

             if (present(ivar1)) ibuff( (n-1)*is+1 : n*is) =   ivar1(i)
             if (present(ivar2)) ibuff( (n-1)*is+1 : n*is) =   ivar2(:,i)
             if (present(ivar3)) ibuff( (n-1)*is+1 : n*is) = [ ivar3(:,:,i)   ]
             if (present(ivar4)) ibuff( (n-1)*is+1 : n*is) = [ ivar4(:,:,:,i) ]

          enddo
          call MPI_Send(ibuff, locbuff, MPI_INTEGER, 0, itag2, MPI_COMM_OLAM, ier)

       elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then

          allocate(rbuff(locbuff))
          do n = 1, nu

             if (present(lpoints)) then
                i = lpoints(n)
             else
                i = n
             endif

             if (present(rvar1)) rbuff( (n-1)*is+1 : n*is) =   rvar1(i)
             if (present(rvar2)) rbuff( (n-1)*is+1 : n*is) =   rvar2(:,i)
             if (present(rvar3)) rbuff( (n-1)*is+1 : n*is) = [ rvar3(:,:,i)   ]
             if (present(rvar4)) rbuff( (n-1)*is+1 : n*is) = [ rvar4(:,:,:,i) ]
          enddo

          call MPI_Send(rbuff, locbuff, MPI_REAL, 0, itag2, MPI_COMM_OLAM, ier)

       elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then

          allocate(dbuff(locbuff))
          do n = 1, nu
             if (present(lpoints)) then
                i = lpoints(n)
             else
                i = n
             endif

             if (present(dvar1)) dbuff( (n-1)*is+1 : n*is) =   dvar1(i)
             if (present(dvar2)) dbuff( (n-1)*is+1 : n*is) =   dvar2(:,i)
             if (present(dvar3)) dbuff( (n-1)*is+1 : n*is) = [ dvar3(:,:,i)   ]
             if (present(dvar4)) dbuff( (n-1)*is+1 : n*is) = [ dvar4(:,:,:,i) ]

          enddo
          call MPI_Send(dbuff, locbuff, MPI_REAL8, 0, itag2, MPI_COMM_OLAM, ier)

       elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then

          allocate(bbuff(locbuff))
          do n = 1, nu
             if (present(lpoints)) then
                i = lpoints(n)
             else
                i = n
             endif

             if (present(bvar1)) bbuff( (n-1)*is+1 : n*is) =   bvar1(i)
             if (present(bvar2)) bbuff( (n-1)*is+1 : n*is) =   bvar2(:,i)
             if (present(bvar3)) bbuff( (n-1)*is+1 : n*is) = [ bvar3(:,:,i)   ]
             if (present(bvar4)) bbuff( (n-1)*is+1 : n*is) = [ bvar4(:,:,:,i) ]

          enddo
          call MPI_Send(bbuff, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, ier)

       elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then

          allocate(bbuff(locbuff))
          do n = 1, nu
             if (present(lpoints)) then
                i = lpoints(n)
             else
                i = n
             endif

             if (present(lvar1)) bbuff( (n-1)*is+1 : n*is) =   lbuf1(i)
             if (present(lvar2)) bbuff( (n-1)*is+1 : n*is) =   lbuf2(:,i)
             if (present(lvar2)) bbuff( (n-1)*is+1 : n*is) = [ lbuf3(:,:,i) ]

          enddo
          call MPI_Send(bbuff, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, ier)

       endif

       call MPI_Waitall(2, ireqs, MPI_STATUSES_IGNORE, ier)

    elseif (myrank == 0) then

       ! Rank 0 only: collect data from other nodes and write to disk

       call MPI_Wait(ireqs(1), MPI_STATUS_IGNORE, ier)

       base = maxval(nus(2:mgroupsize))
       maxbuff = base
       do n = ndims-1, 1, -1
          maxbuff = maxbuff * dims(n)
       enddo

       allocate(points(base))

       if     (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
          allocate(ibuff(maxbuff))
       elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
          allocate(rbuff(maxbuff))
       elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
          allocate(dbuff(maxbuff))
       elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
          allocate(bbuff(maxbuff))
       elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)                    ) then
          allocate(bbuff(maxbuff))
       endif

       do n = 2, mgroupsize

          if (nus(n) < 1) cycle

          ! Prepare memory and options for the write

          call MPI_Recv(points, base, MPI_INTEGER, n-1, itag1, MPI_COMM_OLAM, &
                        MPI_STATUS_IGNORE, ier)

          dimsn(1:ndims) = [ dims(1:ndims-1), nus(n) ]

          call fh5_prepare_write(ndims, dimsn, dsetname, stype, hdferr, &
                                 fcoords=points(1:nus(n)), ifsize=nglobe)
          if (hdferr /= 0) then
             print*, "shdf5_orec_no_phdf5_rank0: can't prepare requested field:", trim(dsetname)
             call fh5_close_write(hdferr, no_dsetid=.true.)
             return
          endif

          ! Write the dataset.

          if (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then

             call MPI_Recv( ibuff, maxbuff, MPI_INTEGER, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE, ier )
             call fh5_write( ibuff, hdferr )

          elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then

             call MPI_Recv( rbuff, maxbuff, MPI_REAL, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE, ier )
             call fh5_write( rbuff, hdferr )

          elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then

             call MPI_Recv( dbuff, maxbuff, MPI_REAL8,   n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE, ier )
             call fh5_write( dbuff, hdferr )

          elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then

             call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE, ier )
             call fh5_write( bbuff, hdferr )

          elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then

             call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE, ier )
             call fh5_write( bbuff, hdferr )

          else
             print*, 'Incorrect or missing data field argument in shdf5_orec2'
             stop    'shdf5_orec2: bad data field'
          endif

          if (hdferr /= 0) then
             print*, 'In shdf5_orec2: hdf5 write error =', hdferr
             stop    'shdf5_orec2: hdf5 write error'
          endif

          call fh5_close_write(hdferr)

       enddo

    endif

  end subroutine shdf5_orec_no_phdf5_rank0


  subroutine shdf5_orec_no_phdf5_allranks()

    use mem_para, only: mgroupsize, myrank, olam_mpi_barrier

    use hdf5_f2f, only: fh5f_open, fh5_prepare_write, fh5_write, fh5_close_write, &
                        fh5f_close

    import, only: dsetname, ndims, dims, nglobe, lpoints, gpoints, stype, &
                  bvar1, bvar2, bvar3, bvar4, ivar1, ivar2, ivar3, ivar4, &
                  rvar1, rvar2, rvar3, rvar4, dvar1, dvar2, dvar3, dvar4, &
                  lvar1, lvar2, lvar3, r8, i1, fname, lbuf1, lbuf2, lbuf3

    implicit none (external, type)

    integer :: hdferr, maxbuff, locbuff
    integer :: nu, ier, base, n, i, is

    if (myrank == 0) call fh5f_close(hdferr)

    call olam_mpi_barrier()

    do n = 1, mgroupsize-1

       if (myrank == n) then

          call fh5f_open(fname, 2, hdferr)

          call fh5_prepare_write(ndims, dims, dsetname, stype, hdferr, &
                                 mcoords=lpoints, fcoords=gpoints, ifsize=nglobe)

          if     (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
          elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
          elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
          elseif (present(lvar1)) then ; call fh5_write(lbuf1, hdferr)
          elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

          elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
          elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
          elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
          elseif (present(lvar2)) then ; call fh5_write(lbuf2, hdferr)
          elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

          elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
          elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
          elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
          elseif (present(lvar3)) then ; call fh5_write(lbuf3, hdferr)
          elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)

          elseif (present(ivar4)) then ; call fh5_write(ivar4, hdferr)
          elseif (present(rvar4)) then ; call fh5_write(rvar4, hdferr)
          elseif (present(dvar4)) then ; call fh5_write(dvar4, hdferr)
          elseif (present(bvar4)) then ; call fh5_write(bvar4, hdferr)

          else
             print*, 'Incorrect or missing data field argument in shdf5_orec'
             stop    'shdf5_orec: bad data field'
          endif

          if (hdferr /= 0) then
             print*, 'In shdf5_orec: hdf5 write error =', hdferr
             stop    'shdf5_orec: hdf5 write error'
          endif

          call fh5_close_write(hdferr)
          call fh5f_close     (hdferr)

          !write(*,'(A,I0,A)') "Rank ", myrank, " finished writing " // trim(dsetname)
       endif

       call olam_mpi_barrier()
    enddo

    if (myrank == 0) call fh5f_open(fname, 2, hdferr)

  end subroutine shdf5_orec_no_phdf5_allranks
#endif

end subroutine shdf5_orec

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars,     &
                                          bvar1,ivar1,rvar1,dvar1,lvar1,     &
                                          bvar2,ivar2,rvar2,dvar2,lvar2,     &
                                          bvar3,ivar3,rvar3,dvar3,lvar3,     &
                                          bvar4,ivar4,rvar4,dvar4,           &
                                          points, start, counts,             &
                                          imissing, rmissing, dmissing,      &
                                          dimname, attached_dimnames, units, &
                                          standard_name, long_name,          &
                                          rscale, roffset, dscale, doffset)

#ifdef OLAM_MPI
  use mpi_f08,   only: MPI_Bcast, MPI_INTEGER, MPI_REAL, MPI_REAL8, MPI_LOGICAL, &
                       MPI_INTEGER1, MPI_CHARACTER
  use mem_para,  only: olam_mpi_barrier, MPI_COMM_OLAM
#endif
  use misc_coms, only: io6
  use hdf5_f2f,  only: fh5_prepare_read, fh5_read, fh5_close_read, fh5f_query_dimname, &
                       fh5_get_attached_scales, fh5f_read_attribute, fh5_close_dataset
  import,        only: i1, r8, do_comm, do_hdf5, indepio

  implicit none (external, type)

  character(*), intent(IN)             :: dsetname ! Dataset name
  integer,      intent(IN)             :: ndims    ! Number of dimensions or rank
  integer,      intent(IN), contiguous :: dims(:)  ! Dataset dimensions

! Array and scalar arguments for different types. Only specify one in each call.
  integer(i1), intent(inout), optional :: bvars
  integer,     intent(inout), optional :: ivars
  real,        intent(inout), optional :: rvars
  real(r8),    intent(inout), optional :: dvars
  logical,     intent(inout), optional :: lvars

  integer(i1), intent(inout), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,     intent(inout), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,        intent(inout), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),    intent(inout), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,     intent(inout), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(IN), optional, contiguous :: points(:)
  integer,      intent(IN), optional, contiguous :: start (:)
  integer,      intent(IN), optional, contiguous :: counts(:)

! Optional arrays to read common NetCDF convention attributes or dimension information
  character(*), intent(inout), optional :: units
  character(*), intent(inout), optional :: dimname
  character(*), intent(inout), optional :: attached_dimnames(:)
  character(*), intent(inout), optional :: standard_name, long_name
  integer,      intent(inout), optional :: imissing
  real,         intent(inout), optional :: rmissing, rscale, roffset
  real(r8),     intent(inout), optional :: dmissing, dscale, doffset

! Local variables
  integer                  :: hdferr  ! Error flag
  logical                  :: exists, do_bcst
  integer                  :: n
  integer(i1)              :: lbufs
  integer(i1), allocatable :: lbuf1(:), lbuf2(:,:), lbuf3(:,:,:)

  ! Maybe we can get around specifying all possible ranks and types by using
  ! assumed-rank and/or assumed-type arrays and C interoperability once
  ! most compilers support this

  n = count( [present(bvars), present(bvar1), present(bvar2), present(bvar3), present(bvar4), &
              present(ivars), present(ivar1), present(ivar2), present(ivar3), present(ivar4), &
              present(rvars), present(rvar1), present(rvar2), present(rvar3), present(rvar4), &
              present(dvars), present(dvar1), present(dvar2), present(dvar3), present(dvar4), &
              present(lvars), present(lvar1), present(lvar2), present(lvar3)] )

  if (n /= 1) then
     write(io6,*) 'Error in shdf5_irec reading ' // trim(dsetname)
     write(io6,*) 'Only one variable type should be specified for HDF5 input.'
     return
  endif

! Check dimensions

  if (ndims <= 0 .or. minval(dims(1:ndims)) <= 0) then
     write(io6,*) 'Dimension error in shdf5_irec:', ndims, dims(1:ndims)
     write(io6,*) 'Dataset name: ' // dsetname
     stop    'shdf5_irec: bad dims'
  endif

! write(io6,'(A,8(1x,I0))') " Reading: "//trim(dsetname), dims(1:ndims)

  do_bcst = do_comm .and. (.not. present(points))

  ! Prepare file and memory space for the read

  if (do_hdf5) call fh5_prepare_read(dsetname, ndims, dims, hdferr, coords=points, &
                                     start=start, counts=counts)
#ifdef OLAM_MPI
  if (do_comm) call MPI_Bcast(hdferr, 1, MPI_INTEGER, 0, MPI_COMM_OLAM)
#endif

  if (hdferr < 0) then
     print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
     return
  endif

  ! Read data, and broadcast to other nodes for a parallel run without
  ! parallel HDF5 or partial I/O

  if (present(ivars)) then

     if (do_hdf5) call fh5_read(ivars, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(ivars, 1, MPI_INTEGER,  0, MPI_COMM_OLAM)
#endif

  elseif (present(rvars)) then

     if (do_hdf5) call fh5_read(rvars, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(rvars, 1, MPI_REAL,     0, MPI_COMM_OLAM)
#endif

  elseif (present(dvars)) then

     if (do_hdf5) call fh5_read(dvars, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(dvars, 1, MPI_REAL8,    0, MPI_COMM_OLAM)
#endif

  elseif (present(lvars)) then

     if (do_hdf5) then
        call fh5_read(lbufs, hdferr)
        lvars = ( lbufs /= 0_i1 )
     endif
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(lvars, 1, MPI_LOGICAL,  0, MPI_COMM_OLAM)
#endif

  elseif (present(bvars)) then

     if (do_hdf5) call fh5_read(bvars, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(bvars, 1, MPI_INTEGER1, 0, MPI_COMM_OLAM)
#endif

  elseif (present(ivar1)) then

     if (do_hdf5) call fh5_read(ivar1, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(ivar1, size(ivar1), MPI_INTEGER,  0, MPI_COMM_OLAM)
#endif

  elseif (present(rvar1)) then

     if (do_hdf5) call fh5_read(rvar1, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(rvar1, size(rvar1), MPI_REAL,     0, MPI_COMM_OLAM)
#endif

  elseif (present(dvar1)) then

     if (do_hdf5) call fh5_read(dvar1, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(dvar1, size(dvar1), MPI_REAL8,    0, MPI_COMM_OLAM)
#endif
  elseif (present(lvar1)) then

     if (do_hdf5) then
        allocate( lbuf1( size(lvar1) ) )
        call fh5_read(lbuf1, hdferr)
        lvar1 = ( lbuf1 /= 0_i1 )
        deallocate(lbuf1)
     endif
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(lvar1, size(lvar1), MPI_LOGICAL,  0, MPI_COMM_OLAM)
#endif

  elseif (present(bvar1)) then

     if (do_hdf5) call fh5_read(bvar1, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(bvar1, size(bvar1), MPI_INTEGER1, 0, MPI_COMM_OLAM)
#endif

  elseif (present(ivar2)) then

     if (do_hdf5) call fh5_read(ivar2, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(ivar2, size(ivar2), MPI_INTEGER,  0, MPI_COMM_OLAM)
#endif

  elseif (present(rvar2)) then

     if (do_hdf5) call fh5_read(rvar2, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(rvar2, size(rvar2), MPI_REAL,     0, MPI_COMM_OLAM)
#endif

  elseif (present(dvar2)) then

     if (do_hdf5) call fh5_read(dvar2, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(dvar2, size(dvar2), MPI_REAL8,    0, MPI_COMM_OLAM)
#endif

  elseif (present(lvar2)) then

     if (do_hdf5) then
        allocate( lbuf2( size(lvar2,dim=1), size(lvar2,dim=2) ) )
        call fh5_read(lbuf2, hdferr)
        lvar2 = ( lbuf2 /= 0_i1 )
        deallocate(lbuf2)
     endif
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(lvar2, size(lvar2), MPI_LOGICAL,  0, MPI_COMM_OLAM)
#endif

  elseif (present(bvar2)) then

     if (do_hdf5) call fh5_read(bvar2, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(bvar2, size(bvar2), MPI_INTEGER1, 0, MPI_COMM_OLAM)
#endif

  elseif (present(ivar3)) then

     if (do_hdf5) call fh5_read(ivar3, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(ivar3, size(ivar3), MPI_INTEGER,  0, MPI_COMM_OLAM)
#endif

  elseif (present(rvar3)) then

     if (do_hdf5) call fh5_read(rvar3, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(rvar3, size(rvar3), MPI_REAL,     0, MPI_COMM_OLAM)
#endif

  elseif (present(dvar3)) then

     if (do_hdf5) call fh5_read(dvar3, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(dvar3, size(dvar3), MPI_REAL8,    0, MPI_COMM_OLAM)
#endif

  elseif (present(lvar3)) then

     if (do_hdf5) then
        allocate( lbuf3( size(lvar3,dim=1), size(lvar3,dim=2), size(lvar3,dim=3) ) )
        call fh5_read(lbuf3, hdferr)
        lvar3 = ( lbuf3 /= 0_i1 )
        deallocate(lbuf3)
     endif
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(lvar3, size(lvar3), MPI_LOGICAL,  0, MPI_COMM_OLAM)
#endif

  elseif (present(bvar3)) then

     if (do_hdf5) call fh5_read(bvar3, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(bvar3, size(bvar3), MPI_INTEGER1, 0, MPI_COMM_OLAM)
#endif

  elseif (present(ivar4)) then

     if (do_hdf5) call fh5_read(ivar4, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(ivar4, size(ivar4), MPI_INTEGER,  0, MPI_COMM_OLAM)
#endif

  elseif (present(rvar4)) then

     if (do_hdf5) call fh5_read(rvar4, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(rvar4, size(rvar4), MPI_REAL,     0, MPI_COMM_OLAM)
#endif

  elseif (present(dvar4)) then

     if (do_hdf5) call fh5_read(dvar4, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(dvar4, size(dvar4), MPI_REAL8,    0, MPI_COMM_OLAM)
#endif

  elseif (present(bvar4)) then

     if (do_hdf5) call fh5_read(bvar4, hdferr)
#ifdef OLAM_MPI
     if (do_bcst) call MPI_Bcast(bvar4, size(bvar4), MPI_INTEGER1, 0, MPI_COMM_OLAM)
#endif

  else
     print*,'Incorrect or missing data field argument in shdf5_irec'
     print*, 'field = ', dsetname
     stop    'shdf5_irec: bad data field'
  endif

#ifdef OLAM_MPI
  if (do_comm) call MPI_Bcast(hdferr, 1, MPI_INTEGER, 0, MPI_COMM_OLAM)
#endif

  if (hdferr /= 0) then
     print*, 'shdf5_irec: call fh5d_read: hdf5 error =', hdferr
     print*, 'Error reading ', trim(dsetname)
     print*, 'ndims = ', ndims
     print*, 'dims  = ', dims(1:ndims)
     stop
  endif

  ! Close the read buffers but leave dataset open

  if (do_hdf5) call fh5_close_read(hdferr, no_dsetid=.true.)

  ! Is this variable a dimension?

  if (present(dimname)) then
     dimname = ' '
     if (do_hdf5) call fh5f_query_dimname(dimname)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(dimname, len(dimname), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  ! Does the variable have dimension scales?

  if (present(attached_dimnames)) then
     attached_dimnames(:) = ' '
     if (do_hdf5) call fh5_get_attached_scales(attached_dimnames)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(attached_dimnames, len(attached_dimnames)*size(attached_dimnames), &
                                 MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  ! Read any attributes

  if (present(units)) then
     units = ' '
     if (do_hdf5) call fh5f_read_attribute("units", cvalue=units)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(units, len(units), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(imissing)) then
     imissing = 1
     if (do_hdf5) call fh5f_read_attribute("missing_value", ivalue=imissing, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(imissing, 1, MPI_INTEGER, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(rmissing)) then
     rmissing = 1.0
     if (do_hdf5) call fh5f_read_attribute("missing_value", rvalue=rmissing, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(rmissing, 1, MPI_REAL, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(dmissing)) then
     dmissing = 1.0_r8
     if (do_hdf5) call fh5f_read_attribute("missing_value", dvalue=dmissing, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(dmissing, 1, MPI_REAL8, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(rscale)) then
     rscale = 1.0
     if (do_hdf5) call fh5f_read_attribute("scale_factor", rvalue=rscale, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(rscale, 1, MPI_REAL, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(dscale)) then
     dscale = 1._r8
     if (do_hdf5) call fh5f_read_attribute("scale_factor", dvalue=dscale, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(dscale, 1, MPI_REAL8, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(roffset)) then
     roffset = 0.0
     if (do_hdf5) call fh5f_read_attribute("add_offset", rvalue=roffset, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(roffset, 1, MPI_REAL, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(doffset)) then
     doffset = 0.0_r8
     if (do_hdf5) call fh5f_read_attribute("add_offset", dvalue=doffset, exists=exists)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(doffset, 1, MPI_REAL8, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(long_name)) then
     long_name = ' '
     if (do_hdf5) call fh5f_read_attribute("long_name", cvalue=long_name)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(dimname, len(dimname), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  if (present(standard_name)) then
     standard_name = ' '
     if (do_hdf5) call fh5f_read_attribute("standard_name", cvalue=standard_name)
#ifdef OLAM_MPI
     if (do_comm) call MPI_Bcast(standard_name, len(standard_name), MPI_CHARACTER, 0, MPI_COMM_OLAM)
#endif
  endif

  if (do_hdf5) call fh5_close_dataset(hdferr)

#ifdef OLAM_MPI

  ! If all nodes are independently reading the same file outside of parallel HDF5,
  ! wait here for all nodes to finish

  if (indepio) call olam_mpi_barrier()

  if (do_comm .and. present(points)) then
     call shdf5_irec_no_phdf5_rank0()
  endif

Contains

  subroutine shdf5_irec_no_phdf5_rank0()

    use mem_para, only: mgroupsize, myrank, MPI_COMM_OLAM

    use mpi_f08,  only: MPI_Igather, MPI_INTEGER, MPI_Isend, MPI_Recv, MPI_STATUS_IGNORE, &
                        MPI_REAL, MPI_REAL8, MPI_INTEGER1, MPI_LOGICAL, MPI_Waitall, &
                        MPI_STATUSES_IGNORE, MPI_Wait, MPI_Send, MPI_Request

    use hdf5_f2f, only: fh5_prepare_read, fh5_read, fh5_close_read

    import,       only: bvar1, bvar2, bvar3, bvar4, ivar1, ivar2, ivar3, ivar4, &
                        rvar1, rvar2, rvar3, rvar4, dvar1, dvar2, dvar3, dvar4, &
                        lvar1, lvar2, lvar3, r8, i1, dsetname, ndims, dims, points

    implicit none (external, type)

    integer,     allocatable :: ibuff(:), rpoints(:), nus(:)
    real,        allocatable :: rbuff(:)
    real(r8),    allocatable :: dbuff(:)
    logical,     allocatable :: lbuff(:)
    integer(i1), allocatable :: bbuff(:)
    integer                  :: nu, base, n, i, is, isize
    integer                  :: maxbuff, locbuff, hdferr
    integer,       parameter :: itag1 = 2096
    integer,       parameter :: itag2 = 2097
    type(MPI_Request)        :: ireqs(2)

    ! Collect on rank 0 the number of cells to be read by each process

    nu = size(points)

    if (myrank == 0) then
       allocate(nus(mgroupsize))
    else
       allocate(nus(1))
    endif

    call MPI_Igather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_OLAM, ireqs(1))

    ! Send to rank 0 the list of points and receive the data read in from disk

    if (myrank > 0 .and. nu > 0) then

       call MPI_Isend(points, nu, MPI_INTEGER, 0, itag1, MPI_COMM_OLAM, ireqs(2))

       if     (present(ivar1)) then
          call MPI_Recv(ivar1, size(ivar1), MPI_INTEGER, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(ivar2)) then
          call MPI_Recv(ivar2, size(ivar2), MPI_INTEGER, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(ivar3)) then
          call MPI_Recv(ivar3, size(ivar3), MPI_INTEGER, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(ivar4)) then
          call MPI_Recv(ivar4, size(ivar4), MPI_INTEGER, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)

       elseif (present(rvar1)) then
          call MPI_Recv(rvar1, size(rvar1), MPI_REAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(rvar2)) then
          call MPI_Recv(rvar2, size(rvar2), MPI_REAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(rvar3)) then
          call MPI_Recv(rvar3, size(rvar3), MPI_REAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(rvar4)) then
          call MPI_Recv(rvar4, size(rvar4), MPI_REAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)

       elseif (present(dvar1)) then
          call MPI_Recv(dvar1, size(dvar1), MPI_REAL8, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(dvar2)) then
          call MPI_Recv(dvar2, size(dvar2), MPI_REAL8, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(dvar3)) then
          call MPI_Recv(dvar3, size(dvar3), MPI_REAL8, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(dvar4)) then
          call MPI_Recv(dvar4, size(dvar4), MPI_REAL8, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)

       elseif (present(bvar1)) then
          call MPI_Recv(bvar1, size(bvar1), MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(bvar2)) then
          call MPI_Recv(bvar2, size(bvar2), MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(bvar3)) then
          call MPI_Recv(bvar3, size(bvar3), MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(bvar4)) then
          call MPI_Recv(bvar4, size(bvar4), MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)

       elseif (present(lvar1)) then
          call MPI_Recv(lvar1, size(lvar1), MPI_LOGICAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(lvar2)) then
          call MPI_Recv(lvar2, size(lvar2), MPI_LOGICAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       elseif (present(lvar3)) then
          call MPI_Recv(lvar3, size(lvar3), MPI_LOGICAL, 0, itag2, MPI_COMM_OLAM, MPI_STATUS_IGNORE)
       endif

       call MPI_Waitall(2, ireqs, MPI_STATUSES_IGNORE)

    elseif (myrank == 0) then

       ! Rank 0 only: read data from disk and send to other nodes

       call MPI_Wait(ireqs(1), MPI_STATUS_IGNORE)

       base = maxval(nus(2:mgroupsize))
       maxbuff = base
       locbuff = nu
       is      = 1
       do n = ndims-1, 1, -1
          maxbuff = maxbuff * dims(n)
          locbuff = locbuff * dims(n)
          is      = is      * dims(n)
       enddo

       allocate(rpoints(base))

       if     (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
          allocate(ibuff(maxbuff))
       elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
          allocate(rbuff(maxbuff))
       elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
          allocate(dbuff(maxbuff))
       elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
          allocate(bbuff(maxbuff))
       elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)                    ) then
          allocate(bbuff(maxbuff))
          allocate(lbuff(maxbuff))
       endif

       do n = 2, mgroupsize

          if (nus(n) < 1) cycle

          call MPI_Recv(rpoints, base, MPI_INTEGER, n-1, itag1, MPI_COMM_OLAM, &
                        MPI_STATUS_IGNORE)

          call fh5_prepare_read(dsetname, ndims, dims, hdferr, coords=rpoints(1:nus(n)))

          if (hdferr /= 0) then
             print*, "shdf5_orec2: can't prepare requested field:", trim(dsetname)
             return
          endif

          ! Write the dataset.

          isize = nus(n) * is

          if (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then

             call fh5_read(ibuff, hdferr)
             call MPI_Send(ibuff, isize, MPI_INTEGER, n-1, itag2, MPI_COMM_OLAM)

          elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then

             call fh5_read(rbuff, hdferr)
             call MPI_Send(rbuff, isize, MPI_REAL, n-1, itag2, MPI_COMM_OLAM)

          elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then

             call fh5_read(dbuff, hdferr)
             call MPI_Send(dbuff, isize, MPI_REAL8, n-1, itag2, MPI_COMM_OLAM)

          elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then

             call fh5_read(bbuff, hdferr)
             call MPI_Send(bbuff, isize, MPI_INTEGER1, n-1, itag2, MPI_COMM_OLAM)

          elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then

             call fh5_read(bbuff, hdferr)
             lbuff(1:isize) = ( bbuff(1:isize) /= 0_i1 )
             call MPI_Send(lbuff, isize, MPI_LOGICAL, n-1, itag2, MPI_COMM_OLAM)

          else

             print*, 'Incorrect or missing data field argument in shdf5_orec2'
             stop    'shdf5_orec2: bad data field'

          endif

          if (hdferr /= 0) then
             print*, 'In shdf5_orec2: hdf5 write error =', hdferr
             stop    'shdf5_orec2: hdf5 write error'
          endif

          call fh5_close_read(hdferr)

       enddo

    endif

  end subroutine shdf5_irec_no_phdf5_rank0
#endif

end subroutine shdf5_irec

!===============================================================================

subroutine shdf5_close()

  use hdf5_f2f, only: fh5f_close
  import,       only: do_hdf5, fname, fmode

  implicit none (external, type)

  integer :: hdferr

  ! Close hdf file.

  if (do_hdf5) then
     call fh5f_close(hdferr)
     fname = ' '
     fmode = ' '
  endif

end subroutine shdf5_close

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                               bvar1,ivar1,rvar1,dvar1,lvar1, &
                                               bvar2,ivar2,rvar2,dvar2,lvar2, &
                                               bvar3,ivar3,rvar3,dvar3,lvar3, &
                                               bvar4,ivar4,rvar4,dvar4        )

  import, only: i1, r8, shdf5_irec, shdf5_orec

  implicit none (external, type)

  character(*), intent(in)             :: dsetname, action
  integer,      intent(in)             :: ndims
  integer,      intent(in), contiguous :: dims(:)

  integer(i1),  intent(inout), optional :: bvars
  integer,      intent(inout), optional :: ivars
  real,         intent(inout), optional :: rvars
  real(r8),     intent(inout), optional :: dvars
  logical,      intent(inout), optional :: lvars

  integer(i1),  intent(inout), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(inout), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(inout), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),     intent(inout), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(inout), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (action == 'READ') then

     call shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,dvar4        )
  elseif (action == 'WRITE') then

     call shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,dvar4        )
  else

     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif

end subroutine shdf5_io

!===============================================================================

subroutine shdf5_orec_ll(ndims,dims,dsetname,bvar1,ivar1,rvar1,dvar1,lvar1,  &
                                             bvar2,ivar2,rvar2,dvar2,lvar2,  &
                                             bvar3,ivar3,rvar3,dvar3,lvar3,  &
                                             gpoints, storage_type,          &
                                             units, long_name, positive,     &
                                             imissing, rmissing, dmissing,   &
                                             isdim, dimnames, standard_name, &
                                             cell_methods, dims_chunk        )
  use oname_coms,  only: nl
  use misc_coms,   only: io6
  use hdf5_f2f,    only: fh5_prepare_write_ll, fh5_write, fh5_close_write, &
                         fh5f_create_dim, fh5f_attach_dims, fh5f_write_attribute, &
                         fh5_close_dataset
  import,          only: i1, r8, do_comm, HID_t, FORTRAN_INT1_TYPE, FORTRAN_INT4_TYPE, &
                         FORTRAN_REAL4_TYPE, FORTRAN_REAL8_TYPE, do_hdf5, fname

  implicit none (external, type)

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

  ! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

  ! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional, contiguous :: gpoints(:)

  ! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

  ! Type of variable
  integer(HID_T), intent(in), optional :: storage_type

  ! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

  ! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

  ! Compression/chunking options
  integer,      intent(in), optional :: dims_chunk(ndims) ! Compression dimensions.

  ! Local variables
  integer                  :: hdferr, n
  integer(HID_T)           :: stype
  integer(i1), allocatable :: lbuf1(:), lbuf2(:,:), lbuf3(:,:,:)

  ! Check subroutine arguments

  n = count( [present(bvar1), present(bvar2), present(bvar3), &
              present(ivar1), present(ivar2), present(ivar3), &
              present(rvar1), present(rvar2), present(rvar3), &
              present(dvar1), present(dvar2), present(dvar3), &
              present(lvar1), present(lvar2), present(lvar3)] )

  if (n /= 1) then
     write(io6,*) 'Error in shdf5_orec_ll writing ' // trim(dsetname)
     write(io6,*) 'Only one variable type should be specified for HDF5 output.'
     return
  endif

  ! Check dimensions

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec_ll:', ndims, dims(1:ndims)
     stop    'shdf5_orec_ll: bad dims'
  endif

  ! Prepare the storage type

  if (present(storage_type)) then
     stype = storage_type
  elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then
     stype = FORTRAN_INT1_TYPE
  elseif (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then
     stype = FORTRAN_INT4_TYPE
  elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then
     stype = FORTRAN_REAL4_TYPE
  elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then
     stype = FORTRAN_REAL8_TYPE
  elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
     stype = FORTRAN_INT1_TYPE
  endif

  ! Convert logicals to 1-byte integers, since HDF5 does not handle Fortran logical/bools

  if (present(lvar1)) then
     allocate(lbuf1, source=merge(-1_i1, 0_i1, lvar1))
  elseif (present(lvar2)) then
     allocate(lbuf2, source=merge(-1_i1, 0_i1, lvar2))
  elseif (present(lvar3)) then
     allocate(lbuf3, source=merge(-1_i1, 0_i1, lvar3))
  endif

! write(io6,'(A,8(1x,I0))') " Writing: "//trim(dsetname), dims(1:ndims)

  if (do_hdf5) then

     ! Prepare memory and options for the write

     call fh5_prepare_write_ll( ndims, dims, dsetname, stype, hdferr, &
                                icompress=nl%icompress, &
                                fcoords=gpoints, dims_chunk=dims_chunk )

     if (hdferr < 0) then
        print*, "shdf5_orec_ll: can't prepare requested field:", trim(dsetname)
        return
     endif

     ! Write the dataset

     if     (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
     elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
     elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
     elseif (present(lvar1)) then ; call fh5_write(lbuf1, hdferr)
     elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

     elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
     elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
     elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
     elseif (present(lvar2)) then ; call fh5_write(lbuf2, hdferr)
     elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

     elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
     elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
     elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
     elseif (present(lvar3)) then ; call fh5_write(lbuf3, hdferr)
     elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)
     else
        print*, 'Incorrect or missing data field argument in shdf5_orec_ll'
        stop    'shdf5_orec_ll: bad data field'
     endif

     if (hdferr /= 0) then
        print*, 'In shdf5_orec_ll: hdf5 write error =', hdferr
        stop    'shdf5_orec_ll: hdf5 write error'
     endif

     ! Close the write buffers but keep the data set open

     call fh5_close_write(hdferr, no_dsetid=.true.)

     ! Indicate if this variable is a dimension

     if (present(isdim)) then
        if (isdim) call fh5f_create_dim(dsetname, hdferr)
     endif

     ! Link each variable to its dimensions

     if ((present(dimnames)) .and. (.not. present(isdim))) then
        if (size(dimnames) >= ndims) call fh5f_attach_dims(ndims, dimnames, hdferr)
     endif

     ! Write any dataset attributes to match common NetCDF conventions

     if (present(units)) then
        if (len_trim(units) > 0) call fh5f_write_attribute("units", cvalue=units)
     endif

     if (present(long_name)) then
        if (len_trim(long_name) > 0) call fh5f_write_attribute("long_name", cvalue=long_name)
     endif

     if (present(standard_name)) then
        if (len_trim(standard_name) > 0) call fh5f_write_attribute("standard_name", cvalue=standard_name)
     endif

     if (present(cell_methods)) then
        if (len_trim(cell_methods) > 0) call fh5f_write_attribute("cell_methods", cvalue=cell_methods)
     endif

     if (present(positive)) then
        if (len_trim(positive) > 0) call fh5f_write_attribute("positive", cvalue=positive)
     endif

     if (present(imissing)) then
        call fh5f_write_attribute("missing_value", ivalue=imissing)
     endif

     if (present(rmissing)) then
        call fh5f_write_attribute("missing_value", rvalue=rmissing)
     endif

     if (present(dmissing)) then
        call fh5f_write_attribute("missing_value", dvalue=dmissing)
     endif

     call fh5_close_dataset(hdferr)

  endif

#ifdef OLAM_MPI

  if ( do_comm .and. present(gpoints) ) then

     if (nl%allranks_write_hdf5) then
        call shdf5_orec_ll_no_phdf5_allranks()
     else
        call shdf5_orec_ll_no_phdf5_rank0()
     endif

  endif

Contains

  subroutine shdf5_orec_ll_no_phdf5_rank0()

    use mem_para, only: mgroupsize, myrank, MPI_COMM_OLAM

    use mpi_f08,  only: MPI_Request, MPI_Igather, MPI_INTEGER, MPI_Isend, &
                        MPI_REAL, MPI_REAL8, MPI_INTEGER1, MPI_Waitall, &
                        MPI_STATUSES_IGNORE, MPI_Wait, MPI_STATUS_IGNORE, &
                        MPI_Recv, MPI_Send

    use hdf5_f2f, only: fh5_prepare_write_ll, fh5_write, fh5_close_write

    import,       only: dsetname, ndims, dims, gpoints, stype, r8, i1, &
                        bvar1, bvar2, bvar3, ivar1, ivar2, ivar3, &
                        rvar1, rvar2, rvar3, dvar1, dvar2, dvar3, &
                        lvar1, lvar2, lvar3, lbuf1, lbuf2, lbuf3

    implicit none (external, type)

    integer,     allocatable :: ibuff(:), points(:), nus(:)
    real,        allocatable :: rbuff(:)
    real(r8),    allocatable :: dbuff(:)
    integer(i1), allocatable :: bbuff(:)
    integer                  :: nu, base, n
    integer                  :: hdferr, maxbuff, locbuff
    integer,       parameter :: itag1 = 3040
    integer,       parameter :: itag2 = 3041
    type(MPI_Request)        :: ireqs(2)

    nu = size(gpoints)

    if (myrank == 0) then
       allocate(nus(mgroupsize))
    else
       allocate(nus(1))
    endif

    call MPI_Igather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_OLAM, ireqs(1))

    if (myrank > 0 .and. nu > 0) then

       ! Send to rank 0 the list of points and data to be written to disk

       call MPI_Isend(gpoints, nu, MPI_INTEGER, 0, itag1, MPI_COMM_OLAM, ireqs(2))

       locbuff = nu
       do n = 3, ndims
          locbuff = locbuff * dims(n)
       enddo

       if     (present(ivar1)) then
          call MPI_Send(ivar1, locbuff, MPI_INTEGER,  0, itag2, MPI_COMM_OLAM)
       elseif (present(rvar1)) then
          call MPI_Send(rvar1, locbuff, MPI_REAL,     0, itag2, MPI_COMM_OLAM)
       elseif (present(dvar1)) then
          call MPI_Send(dvar1, locbuff, MPI_REAL8,    0, itag2, MPI_COMM_OLAM)
       elseif (present(lvar1)) then
          call MPI_Send(lbuf1, locbuff, MPI_INTEGER1,  0, itag2, MPI_COMM_OLAM)
       elseif (present(bvar1)) then
          call MPI_Send(bvar1, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM)

       elseif (present(ivar2)) then
          call MPI_Send(ivar2, locbuff, MPI_INTEGER,  0, itag2, MPI_COMM_OLAM)
       elseif (present(rvar2)) then
          call MPI_Send(rvar2, locbuff, MPI_REAL,     0, itag2, MPI_COMM_OLAM)
       elseif (present(dvar2)) then
          call MPI_Send(dvar2, locbuff, MPI_REAL8,    0, itag2, MPI_COMM_OLAM)
       elseif (present(lvar2)) then
          call MPI_Send(lbuf2, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM)
       elseif (present(bvar2)) then
          call MPI_Send(bvar2, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM)

       elseif (present(ivar3)) then
          call MPI_Send(ivar3, locbuff, MPI_INTEGER,  0, itag2, MPI_COMM_OLAM)
       elseif (present(rvar3)) then
          call MPI_Send(rvar3, locbuff, MPI_REAL,     0, itag2, MPI_COMM_OLAM)
       elseif (present(dvar3)) then
          call MPI_Send(dvar3, locbuff, MPI_REAL8,    0, itag2, MPI_COMM_OLAM)
       elseif (present(lvar3)) then
          call MPI_Send(lbuf3, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM)
       elseif (present(bvar3)) then
          call MPI_Send(bvar3, locbuff, MPI_INTEGER1, 0, itag2, MPI_COMM_OLAM)
       endif

       call MPI_Waitall(2, ireqs, MPI_STATUSES_IGNORE)

    elseif (myrank == 0) then

       ! Rank 0 only: collect data from other nodes and write to disk

       call MPI_Wait(ireqs(1), MPI_STATUS_IGNORE)

       base = maxval(nus(2:mgroupsize))
       maxbuff = base
       do n = 3, ndims
          maxbuff = maxbuff * dims(n)
       enddo

       allocate(points(base))

       if     (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then
          allocate(ibuff(maxbuff))
       elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then
          allocate(rbuff(maxbuff))
       elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then
          allocate(dbuff(maxbuff))
       elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
          allocate(bbuff(maxbuff))
       elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then
          allocate(bbuff(maxbuff))
       endif

       do n = 2, mgroupsize

          if (nus(n) < 1) cycle

          call MPI_Recv(points, base, MPI_INTEGER, n-1, itag1, MPI_COMM_OLAM, &
                        MPI_STATUS_IGNORE)

          call fh5_prepare_write_ll(ndims, dims, dsetname, stype, hdferr, &
                                    fcoords=points(1:nus(n)))

          if (hdferr /= 0) then
             print*, "shdf5_orec_ll2: can't prepare requested field:", trim(dsetname)
             return
          endif

          ! Write the dataset.

          if (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then

             call MPI_Recv( ibuff, maxbuff, MPI_INTEGER, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE )
             call fh5_write( ibuff, hdferr )

          elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then

             call MPI_Recv( rbuff, maxbuff, MPI_REAL, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE )
             call fh5_write( rbuff, hdferr )

          elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then

             call MPI_Recv( dbuff, maxbuff, MPI_REAL8, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE )
             call fh5_write( dbuff, hdferr )

          elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then

             call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE )
             call fh5_write( bbuff, hdferr )

          elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then

             call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag2, &
                            MPI_COMM_OLAM, MPI_STATUS_IGNORE )
             call fh5_write( bbuff, hdferr )

          else
             print*, 'Incorrect or missing data field argument in shdf5_orec_ll2'
             stop    'shdf5_orec_ll2: bad data field'
          endif

          if (hdferr /= 0) then
             print*, 'In shdf5_orec_ll2: hdf5 write error =', hdferr, myrank
             stop    'shdf5_orec_ll2: hdf5 write error'
          endif

          call fh5_close_write(hdferr)

       enddo

    endif

  end subroutine shdf5_orec_ll_no_phdf5_rank0


  subroutine shdf5_orec_ll_no_phdf5_allranks()

    use mem_para, only: mgroupsize, myrank, olam_mpi_barrier

    use hdf5_f2f, only: fh5_prepare_write_ll, fh5_write, fh5_close_write, &
                        fh5f_open, fh5f_close

    import,       only: dsetname, ndims, dims, gpoints, stype, &
                        bvar1, bvar2, bvar3, ivar1, ivar2, ivar3, &
                        rvar1, rvar2, rvar3, dvar1, dvar2, dvar3, &
                        lvar1, lvar2, lvar3, r8, i1, fname, lbuf1, lbuf2, lbuf3

    implicit none (external, type)

    integer,     allocatable :: ibuff(:), points(:), nus(:)
    real,        allocatable :: rbuff(:)
    real(r8),    allocatable :: dbuff(:)
    integer(i1), allocatable :: bbuff(:)
    integer                  :: nu, base, n
    integer                  :: hdferr, maxbuff, locbuff

    if (myrank == 0) call fh5f_close(hdferr)

    call olam_mpi_barrier()

    do n = 1, mgroupsize-1

       if (myrank == n) then

          call fh5f_open(fname, 2, hdferr)

          call fh5_prepare_write_ll(ndims, dims, dsetname, stype, hdferr, &
                                    fcoords=gpoints)

          if     (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
          elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
          elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
          elseif (present(lvar1)) then ; call fh5_write(lbuf1, hdferr)
          elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

          elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
          elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
          elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
          elseif (present(lvar2)) then ; call fh5_write(lbuf2, hdferr)
          elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

          elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
          elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
          elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
          elseif (present(lvar3)) then ; call fh5_write(lbuf3, hdferr)
          elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)

          else
             print*, 'Incorrect or missing data field argument in shdf5_orec_ll'
             stop    'shdf5_orec_ll: bad data field'
          endif

          if (hdferr /= 0) then
             print*, 'In shdf5_orec_ll: hdf5 write error =', hdferr
             stop    'shdf5_orec_ll: hdf5 write error'
          endif

          call fh5_close_write(hdferr)
          call fh5f_close     (hdferr)

          !write(*,'(A,I0,A)') "Rank ", myrank, " finished writing " // trim(dsetname)
       endif

       call olam_mpi_barrier()

    enddo

    if (myrank == 0) call fh5f_open(fname, 2, hdferr)

  end subroutine shdf5_orec_ll_no_phdf5_allranks
#endif

end subroutine shdf5_orec_ll

!===============================================================================

subroutine shdf5_write_global_attribute(name, ivalue, rvalue, dvalue, cvalue)

  use misc_coms, only: io6
  use hdf5_f2f,  only: fh5f_write_global_attribute
  import,        only: do_hdf5, r8

  implicit none (external, type)

  character(*),           intent(in) :: name
  integer,      optional, intent(in) :: ivalue
  real,         optional, intent(in) :: rvalue
  real(r8),     optional, intent(in) :: dvalue
  character(*), optional, intent(in) :: cvalue
  integer                            :: n

  n = count( [present(ivalue), present(rvalue), present(dvalue), present(cvalue)] )

  if (n /= 1) then
     write(io6,*) "Error in shdf5_write_global_attribute:"
     write(io6,*) "One (and only one) variable type must be present"
     write(io6,*) "Not reading global attribute " // trim(name)
     return
  endif

  if (do_hdf5) call fh5f_write_global_attribute(name, ivalue=ivalue, rvalue=rvalue, &
                                                      dvalue=dvalue, cvalue=cvalue)

end subroutine shdf5_write_global_attribute

!===============================================================================

subroutine shdf5_read_global_attribute(name, exists, ivalue, rvalue, dvalue, cvalue)

#ifdef OLAM_MPI
  use mpi_f08,   only: MPI_Bcast, MPI_INTEGER, MPI_REAL, MPI_REAL8, MPI_CHARACTER
  use mem_para,  only: MPI_COMM_OLAM
#endif
  use misc_coms, only: io6
  use hdf5_f2f,  only: fh5f_read_global_attribute
  import,        only: r8, do_hdf5, do_comm

  implicit none (external, type)

  character(*),           intent(in)  :: name
  integer,      optional, intent(out) :: ivalue
  real,         optional, intent(out) :: rvalue
  real(r8),     optional, intent(out) :: dvalue
  character(*), optional, intent(out) :: cvalue
  logical,      optional, intent(out) :: exists
  integer                             :: n

  if (present(exists)) exists = .false.

  n = count( [present(ivalue), present(rvalue), present(dvalue), present(cvalue)] )
  if (n /= 1) then
     write(io6,*) "Error in shdf5_read_global_attribute:"
     write(io6,*) "One (and only one) variable type must be present"
     write(io6,*) "Not reading global attribute " // trim(name)
     return
  endif

  if (do_hdf5) call fh5f_read_global_attribute(name, exists=exists, &
                    ivalue=ivalue, rvalue=rvalue, dvalue=dvalue, cvalue=cvalue)

#ifdef OLAM_MPI
  if (do_comm) then
     if (present(ivalue)) then
        call MPI_Bcast(ivalue, 1, MPI_INTEGER, 0, MPI_COMM_OLAM)
     elseif (present(rvalue)) then
        call MPI_Bcast(rvalue, 1, MPI_REAL, 0, MPI_COMM_OLAM)
     elseif (present(dvalue)) then
        call MPI_Bcast(dvalue, 1, MPI_REAL8, 0, MPI_COMM_OLAM)
     elseif (present(cvalue)) then
        call MPI_Bcast(cvalue, len(cvalue), MPI_CHARACTER, 0, MPI_COMM_OLAM)
     endif
  endif
#endif

end subroutine shdf5_read_global_attribute

!===============================================================================

end module hdf5_utils
