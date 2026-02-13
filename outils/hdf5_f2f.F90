module hdf5_f2f

  use                          hdf5, only: HID_T, HSIZE_T
  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8, i2=>int16

  implicit none (external, type)

  logical, parameter :: bigendian = ichar(transfer(1,'a')) == 0

  integer(HID_T)   :: fileid

  integer(HID_T)   :: xferid
  integer(HID_T)   :: xferid_col
  integer(HID_T)   :: xferid_ind

  integer(HID_T)   :: dsetid = -1
  integer(HID_T)   :: mspcid
  integer(HID_T)   :: dspcid

  integer(HSIZE_T) :: dimsf(7)
  integer          :: ndimsf

  logical :: dopario = .false.
  integer :: is_linux = -1

  interface fh5_write
     module procedure         &
          fh5_write_int1,     &
          fh5_write_int2,     &
          fh5_write_integer,  &
          fh5_write_real,     &
          fh5_write_real8
  end interface fh5_write

  interface fh5_read
     module procedure        &
          fh5_read_int1,     &
          fh5_read_int2,     &
          fh5_read_integer,  &
          fh5_read_real,     &
          fh5_read_real8
  end interface fh5_read

  private
  public :: fh5f_open, fh5f_close, fh5f_create, fh5_prepare_write, fh5_write, &
            fh5_close_write, fh5_open_dataset, fh5_close_dataset, fh5_get_info, &
            fh5_prepare_read, fh5_read, fh5_prepare_write_ll, &
            fh5f_write_attribute, fh5f_create_dim, fh5f_attach_dims, &
            fh5f_read_global_attribute, fh5f_write_global_attribute, &
            fh5f_read_attribute, fh5f_query_dimname, fh5_get_attached_scales, &
            fh5_get_chunk_dims, fh5f_exists, fh5_close_read

contains

!===============================================================================

  subroutine fh5f_exists(locfn, exists, pario)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    use mpi_f08,   only: MPI_INFO_NULL
    use mem_para,  only: MPI_COMM_OLAM
    use misc_coms, only: iparallel
    use hdf5,      only: h5pset_fapl_mpio_f, h5pset_all_coll_metadata_ops_f
#endif
#ifdef OLD_HDF5
    use hdf5,      only: h5fis_hdf5_f
#else
    use hdf5,      only: h5fis_accessible_f
#endif
    use hdf5,      only: HID_T, h5pcreate_f, H5P_FILE_ACCESS_F, h5eset_auto_f, &
                         h5pclose_f
    import,        only: dopario

    implicit none (external, type)

    character(*),      intent(IN)  :: locfn
    logical,           intent(OUT) :: exists
    logical, optional, intent(IN)  :: pario

    integer(HID_T)                 :: access_prp
    integer                        :: hdferr

    dopario = .false.

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (iparallel == 1 .and. present(pario)) dopario = pario
#endif

    call h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (dopario) then
       call h5pset_fapl_mpio_f(access_prp, MPI_COMM_OLAM, MPI_INFO_NULL, hdferr)
       call h5pset_all_coll_metadata_ops_f(access_prp, .true., hdferr)
    endif
#endif

    call h5eset_auto_f(0, hdferr)

#ifdef OLD_HDF5
    call h5fis_hdf5_f      (locfn, exists, hdferr)
#else
    call h5fis_accessible_f(locfn, exists, hdferr, access_prp)
#endif

    if (hdferr /= 0) exists = .false.
    call h5eset_auto_f(1, hdferr)

    call h5pclose_f(access_prp, hdferr)

  end subroutine fh5f_exists

!===============================================================================

  subroutine fh5f_open(locfn, iaccess, hdferr, pario)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    use mpi_f08,    only: MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_MIN, &
                          MPI_INFO_NULL, MPI_Bcast
    use mem_para,   only: MPI_COMM_OLAM, myrank
    use misc_coms,  only: iparallel
    use string_lib, only: lowercase, strip_char
    use sys_utils,  only: get_command_as_string
    use hdf5,       only: h5pcopy_f, h5pset_fapl_mpio_f, h5pset_all_coll_metadata_ops_f, &
                          H5Pset_alignment_f, H5Pset_dxpl_mpio_f, h5fd_mpio_collective_f
#endif
#ifndef OLD_HDF5
    use hdf5,       only: H5F_LIBVER_V110_F, h5fis_accessible_f
#else
    use hdf5,       only: h5fis_hdf5_f
#endif
    use hdf5,       only: HID_T, HSIZE_T, SIZE_t, H5F_ACC_RDONLY_F, H5F_ACC_RDWR_F, &
                          h5pcreate_f, H5P_FILE_ACCESS_F, h5pget_cache_f, h5pset_cache_f, &
                          h5pset_sieve_buf_size_f, h5pset_libver_bounds_f, H5F_LIBVER_LATEST_F, &
                          H5Pset_meta_block_size_f, H5Pset_small_data_block_size_f, &
                          h5eset_auto_f, h5pclose_f, h5fopen_f, h5pcreate_f, h5pset_buffer_f, &
                          h5pset_hyper_vector_size_f, h5p_dataset_xfer_f
    use misc_coms,  only: io6
    import,         only: dopario, is_linux, fileid, xferid, xferid_col, xferid_ind

    implicit none (external, type)

    character(*),      intent(IN)  :: locfn
    integer,           intent(IN)  :: iaccess
    integer,           intent(OUT) :: hdferr
    logical, optional, intent(IN)  :: pario

    integer(HID_T)                 :: access_prp
    integer                        :: flags
    integer(HSIZE_T),    parameter :: msize = 1024 * 256
    integer(SIZE_t)                :: nelmts, nbytes
    real                           :: w0
    integer                        :: mdc
    logical                        :: ish5
#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    integer(HSIZE_T)               :: aligns1, aligns2
    integer,             parameter :: alignsmax = 1024 * 1024 * 128
    character(12)                  :: string
    integer                        :: ierr
    integer,                target :: finfo(2), nstripes
    integer,               pointer :: blksiz, fstype
#endif

    dopario = .false.

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)

    if (is_linux < 0) then
       is_linux = 0
       string   = ' '
       call get_environment_variable("OSTYPE", string)
       call lowercase(string)
       if (string(1:5) == 'linux') is_linux = 1

       if (iparallel == 1) then
          call MPI_Allreduce(MPI_IN_PLACE, is_linux, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_OLAM, hdferr)
       endif
    endif

    if (iparallel == 1 .and. present(pario)) dopario = pario

#endif

    ! Turn off default error handling?
    ! call h5eset_auto_f(0, hdferr)

    if(iaccess == 1) flags = H5F_ACC_RDONLY_F
    if(iaccess == 2) flags = H5F_ACC_RDWR_F

    call h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)

    call h5pget_cache_f(access_prp, mdc, nelmts, nbytes, w0, hdferr)
    nbytes = max(nbytes, int(1024 * 1024 * 32, SIZE_t) )
    call h5pset_cache_f(access_prp, mdc, nelmts, nbytes, w0, hdferr)

    nbytes = 1024 * 1024 * 32
    call h5pset_sieve_buf_size_f(access_prp, nbytes, hdferr)

#ifdef OLD_HDF5
    call h5pset_libver_bounds_f(access_prp, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, hdferr)
#else
    call h5pset_libver_bounds_f(access_prp, H5F_LIBVER_V110_F, H5F_LIBVER_LATEST_F, hdferr)
#endif

    ! Increase some buffer sizes
    call H5Pset_meta_block_size_f(access_prp, msize, hdferr)
    call H5Pset_small_data_block_size_f(access_prp, msize, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (dopario) then

       write(io6,*)
       write(io6,'(A)') "Enabling parallel HDF5 input for file " // trim(locfn)

       call h5pset_fapl_mpio_f(access_prp, MPI_COMM_OLAM, MPI_INFO_NULL, hdferr)
       call h5pset_all_coll_metadata_ops_f(access_prp, .true., hdferr)
    endif
#endif

    call h5eset_auto_f(0, hdferr)

#ifdef OLD_HDF5
    call h5fis_hdf5_f      (locfn, ish5, hdferr)
#else
    call h5fis_accessible_f(locfn, ish5, hdferr, access_prp)
#endif

    if (hdferr /= 0) ish5 = .false.
    call h5eset_auto_f(1, hdferr)

    if (.not. ish5) then
       write(io6,*) "File " // trim(locfn) // " does not exist or is not a hdf5/nc4 file."
       write(io6,*) "If this file is an older NetCDF3 format file convert it to NetCDF4 using:"
       write(io6,*) "nccopy -k nc4 INFILE OUTFILE"
       write(io6,*)
       call h5pclose_f(access_prp, hdferr)
       hdferr = -1
       return
    endif

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if ( dopario .and. iaccess == 2 .and. is_linux > 0) then

       finfo = 0
       blksiz => finfo(1)
       fstype => finfo(2)

       if (myrank == 0) then

          ! Get file system blocksize if we are using gpfs, or
          ! stripe size if we are using lustre and the file is striped

          call get_command_as_string( "df --output=fstype | tail -n +2"//trim(locfn), string)
          call lowercase(string)

          if (string(1:4) == "gpfs") then
             fstype = 1

             call get_command_as_string( "stat -f --printf='%s' "//trim(locfn), string)
             if (len_trim(string) > 0) then
                if ( verify(trim(string),"0123456789") == 0 ) read(string,*,iostat=ierr) blksiz
             endif

          elseif (string(1:6) == "lustre") then
             fstype = 2

             nstripes = 1
             call get_command_as_string( "lfs getstripe -c "//trim(locfn), string)
             if (len_trim(string) > 0) then
                call strip_char(string, new_line(' '))
                if ( verify(trim(string),"0123456789") == 0 ) read(string,*,iostat=ierr) nstripes
             endif

             if (nstripes > 1) then
                call get_command_as_string( "lfs getstripe -S "//trim(locfn), string)
                if (len_trim(string) > 0) then
                   call strip_char(string, new_line(' '))
                   if ( verify(trim(string),"0123456789") == 0 ) read(string,*,iostat=ierr) blksiz
                endif
             endif

          endif
       endif

       call MPI_Bcast(finfo, 2, MPI_INTEGER,  0, MPI_COMM_OLAM, hdferr)

       ! Align writes to file-system block size if using gpfs or
       ! stripe size if we are using lustre

       if (blksiz >= 4096) then
          aligns1 = min(blksiz,alignsmax)
          aligns2 = aligns1 / 4
          call H5Pset_alignment_f( access_prp, aligns2, aligns1, hdferr )
       endif

    endif
#endif

    call h5fopen_f(locfn, flags, fileid, hdferr, access_prp)

    call h5pclose_f(access_prp, hdferr)

    call h5pcreate_f(h5p_dataset_xfer_f, xferid, hdferr)

    nbytes = 1024 * 1024 * 128
    call h5pset_buffer_f(xferid, nbytes, hdferr)

    nelmts = 10000
    call h5pset_hyper_vector_size_f(xferid, nelmts, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (dopario) then
       call h5pcopy_f(xferid, xferid_col, hdferr)
       call h5pcopy_f(xferid, xferid_ind, hdferr)
       call H5Pset_dxpl_mpio_f(xferid_col, h5fd_mpio_collective_f, hdferr)
       call h5pclose_f(xferid, hdferr)
    endif
#endif

  end subroutine fh5f_open

!=============================================================================

  subroutine fh5f_close(hdferr)

    use hdf5, only: h5fclose_f, h5pclose_f
    import,   only: fileid, xferid_ind, xferid_col, xferid, dopario

    implicit none (external, type)

    integer, intent(OUT) :: hdferr

    call h5fclose_f(fileid, hdferr)

    if (dopario) then
       call h5pclose_f(xferid_ind, hdferr)
       call h5pclose_f(xferid_col, hdferr)
    else
       call h5pclose_f(xferid, hdferr)
    endif

    dopario = .false.
  end subroutine fh5f_close

!===============================================================================

  subroutine fh5f_create(locfn, iaccess, hdferr, pario)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    use mpi_f08,    only: MPI_Info, MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, MPI_MIN, &
                          MPI_INFO_NULL, MPI_Bcast, MPI_Info_create, MPI_Info_set, MPI_Info_free
    use mem_para,   only: myrank, MPI_COMM_OLAM
    use string_lib, only: lowercase
    use sys_utils,  only: get_command_as_string
    use misc_coms,  only: iparallel
    use hdf5,       only: H5Pset_alignment_f, h5pset_fapl_mpio_f, h5pset_coll_metadata_write_f, &
                          h5pcopy_f, H5Pset_dxpl_mpio_f, h5fd_mpio_collective_f
#endif
#ifndef OLD_HDF5
    use hdf5,       only: H5F_LIBVER_V110_F
#endif
    use misc_coms,  only: io6
    use hdf5,       only: HID_T, HSIZE_T, SIZE_t, H5F_ACC_TRUNC_F, H5F_ACC_EXCL_F, &
                          h5pcreate_f, H5P_FILE_ACCESS_F, h5pget_cache_f, h5pset_cache_f, &
                          h5pset_sieve_buf_size_f, h5pset_libver_bounds_f, H5F_LIBVER_LATEST_F, &
                          H5Pset_meta_block_size_f, H5Pset_small_data_block_size_f, &
                          h5fcreate_f, h5pclose_f, h5pset_buffer_f, h5pset_hyper_vector_size_f, &
                          h5p_dataset_xfer_f
    import,         only: dopario, is_linux, fileid, xferid, xferid_col, xferid_ind

    implicit none (external, type)

    character(*),      intent(IN)  :: locfn
    integer,           intent(IN)  :: iaccess
    integer,           intent(OUT) :: hdferr
    logical, optional, intent(IN)  :: pario

    integer(HID_T)                 :: access_prp
    integer                        :: flags, mdc
    integer(HSIZE_T),    parameter :: msize = 1024 * 256
    integer(SIZE_t)                :: nelmts, nbytes
    real                           :: w0
#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    integer                        :: is, ierr
    character(12)                  :: string
    integer,               target  :: finfo(2)
    integer,               pointer :: blksiz, fstype
    integer(HSIZE_T)               :: aligns1, aligns2
    integer,             parameter :: alignsmax = 1024 * 1024 * 128
    integer,             parameter :: lustre_striping_factor = -1
    integer,             parameter :: lustre_striping_unit   = 1024 * 1024 * 4
    type(MPI_Info)                 :: info
    logical                        :: doinfo
#endif

    dopario = .false.

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (is_linux < 0) then
       is_linux = 0
       string   = ' '
       call get_environment_variable("OSTYPE", string)
       if (string == 'linux') is_linux = 1

       if (iparallel == 1) then
          call MPI_Allreduce(MPI_IN_PLACE, is_linux, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_OLAM, hdferr)
       endif
    endif

    if (iparallel == 1 .and. present(pario)) dopario = pario
#endif

    ! Turn off default error handling?
    ! call h5eset_auto_f(0, hdferr)

    if (iaccess == 1) flags = H5F_ACC_TRUNC_F
    if (iaccess == 2) flags = H5F_ACC_EXCL_F

    call h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)

    call h5pget_cache_f(access_prp, mdc, nelmts, nbytes, w0, hdferr)
    nbytes = max(nbytes, int(1024 * 1024 * 32, SIZE_t) )
    call h5pset_cache_f(access_prp, mdc, nelmts, nbytes, w0, hdferr)

    nbytes = 1024 * 1024 * 32
    call h5pset_sieve_buf_size_f(access_prp, nbytes, hdferr)

#ifdef OLD_HDF5
    call h5pset_libver_bounds_f(access_prp, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, hdferr)
#else
    call h5pset_libver_bounds_f(access_prp, H5F_LIBVER_V110_F, H5F_LIBVER_LATEST_F, hdferr)
#endif

    ! Increase some buffer sizes
    call H5Pset_meta_block_size_f(access_prp, msize, hdferr)
    call H5Pset_small_data_block_size_f(access_prp, msize, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (dopario) then

       doinfo = .false.
       finfo  = 0

       if (is_linux > 0) then

          blksiz => finfo(1)
          fstype => finfo(2)

          if (myrank == 0) then

             ! Get file system blocksize if we are using gpfs, or
             ! stripe size if we are using lustre and the file is striped

             is = scan(locfn, "/", .true.)

             if (is == 0) then
                call get_command_as_string( "df --output=fstype ./ | tail -n +2", string)
             else
                call get_command_as_string( "df --output=fstype "//locfn(1:is)//" | tail -n +2", string)
             endif
             call lowercase(string)

             if (string(1:4) == "gpfs") then
                fstype = 1

                if (is == 0) then
                   call get_command_as_string( "stat -f --printf='%s' ./", string)
                else
                   call get_command_as_string( "stat -f --printf='%s' "//locfn(1:is), string)
                endif

                if (len_trim(string) > 0) then
                   if ( verify(trim(string),"0123456789") == 0 ) read(string,*,iostat=ierr) blksiz
                endif

             elseif (string(1:6) == "lustre") then

                fstype = 2
                blksiz = lustre_striping_unit

             endif

          endif

          call MPI_Bcast(finfo, 2, MPI_INTEGER,  0, MPI_COMM_OLAM, hdferr)

          ! Align writes to file-system block size if using gpfs or
          ! stripe size if we are using lustre

          if (blksiz >= 4096) then
             aligns1 = min(blksiz,alignsmax)
             aligns2 = aligns1 / 4
             call H5Pset_alignment_f( access_prp, aligns2, aligns1, hdferr )
          endif

          if (fstype == 2) then
             doinfo = .true.
             call MPI_Info_create(info, hdferr)
             write(string,'(I0)') lustre_striping_unit
             call MPI_Info_set(info, "striping_unit", string, hdferr)
             call MPI_Info_set(info, "striping_factor", "-1", hdferr)
          endif

       endif

       write(io6,*)
       write(io6,'(A)') "Enabling parallel HDF5 output for file " // trim(locfn)

       if (doinfo) then
          call h5pset_fapl_mpio_f(access_prp, MPI_COMM_OLAM, info, hdferr)
       else
          call h5pset_fapl_mpio_f(access_prp, MPI_COMM_OLAM, MPI_INFO_NULL, hdferr)
       endif
       call h5pset_coll_metadata_write_f(access_prp, .true., hdferr)
       if (doinfo) call MPI_Info_free(info, hdferr)
    endif
#endif

    call h5fcreate_f(locfn, flags, fileid, hdferr, access_prp=access_prp)

    if (hdferr < 0) then
       write(io6,*) "File " // trim(locfn) // " could not be created by the HDF5 library."
       write(io6,*)
       call h5pclose_f(access_prp, hdferr)
       return
    endif

    call H5Pclose_f(access_prp, hdferr)

    call h5pcreate_f(h5p_dataset_xfer_f, xferid, hdferr)

    nbytes = 1024 * 1024 * 128
    call h5pset_buffer_f(xferid, nbytes, hdferr)

    nelmts = 10000
    call h5pset_hyper_vector_size_f(xferid, nelmts, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (dopario) then
       call h5pcopy_f(xferid, xferid_col, hdferr)
       call h5pcopy_f(xferid, xferid_ind, hdferr)
       call H5Pset_dxpl_mpio_f(xferid_col, h5fd_mpio_collective_f, hdferr)
       call h5pclose_f(xferid, hdferr)
    endif
#endif

  end subroutine fh5f_create

!===============================================================================

  subroutine fh5_prepare_write(ndims, dims, dname, stype, hdferr, icompress, &
                               mcoords, fcoords, ifsize, dims_chunk)

    use misc_coms, only: io6
    use hdf5,      only: HID_T, HSIZE_T, H5P_DEFAULT_F, h5iis_valid_f, H5Lexists_f, &
                         H5Dopen_f, H5Dget_space_status_f, H5D_SPACE_STS_ALLOCATED_F, &
                         H5D_SPACE_STS_PART_ALLOCATED_F, h5pcreate_f, h5pset_chunk_f, &
                         H5P_DATASET_CREATE_F, h5pset_shuffle_f, h5pset_deflate_f, &
                         h5screate_simple_f, h5dcreate_f, h5sclose_f, h5pclose_f, &
                         h5sselect_none_f, H5S_SELECT_SET_F, h5sselect_hyperslab_f, &
                         H5S_SELECT_OR_F, H5Sselect_all_f, h5pcopy_f
    import,        only: dimsf, ndimsf, fileid, dsetid, mspcid, dspcid, dopario, &
                         xferid, xferid_col

    implicit none (external, type)

    integer,        intent(IN)                       :: ndims
    integer,        intent(IN)                       :: dims(ndims)
    character(*),   intent(IN)                       :: dname
    integer(HID_T), intent(IN)                       :: stype
    integer,        intent(OUT)                      :: hdferr
    integer,        intent(IN), optional             :: icompress, ifsize
    integer,        intent(IN), optional, contiguous :: mcoords(:), fcoords(:)
    integer,        intent(IN), optional             :: dims_chunk(ndims)

    logical          :: docompress, valid
    integer          :: j, iop, flag
    integer(HSIZE_T) :: dims_file(ndims), dims_compress(ndims)
    integer(HSIZE_T) :: offset(ndims), countf(ndims)
    integer(HID_T)   :: propid
    integer(HID_T)   :: fspcid

    ! Output dimensions

    ndimsf         = ndims
    dimsf(1:ndims) = dims
    dims_file      = dims

    ! For distributed/partial output, global data (file) size will be different
    ! then the local array size

    if (present(fcoords) .and. present(ifsize)) then
       dims_file(ndims) = ifsize
    endif

    hdferr     = 0
    propid     = H5P_DEFAULT_F
    docompress = .false.

    ! Is the HDF5 file open?

    call h5iis_valid_f(fileid, valid, hdferr)
    if (hdferr < 0 .or. .not. valid) then
       hdferr = -1
       write(io6,*) "fh5_prepare_write: Trying to write " // trim(dname) // &
                    " but no HDF5 file is open."
       return
    endif

    ! Make sure dsetid is not already open

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (valid) then
       write(io6,*) "fh5_prepare_write: dsetid is already open"
       hdferr = -1
       return
    endif

    call H5Lexists_f(fileid, dname, valid, hdferr)
    if (valid) then

       call H5Dopen_f(fileid, dname, dsetid, hdferr)
       call H5Dget_space_status_f(dsetid, flag, hdferr)

       if ( flag /= H5D_SPACE_STS_ALLOCATED_F .and. &
            flag /= H5D_SPACE_STS_PART_ALLOCATED_F ) then
          write(io6,*) "HDF5 error: a data set has already been created for variable " // trim(dname)
          write(io6,*) "in the HDF5 file but has no space allocated to it. Don't know what to do here."
          hdferr = -1
          return
       endif

    else

       ! Setup HFS compression if we are not doing partial writes

       if ( present(icompress) .and. (.not. present(fcoords)) ) then

          if (present(dims_chunk)) then
             dims_compress = dims_chunk
          else
             dims_compress = dims
          endif

          if (icompress > 0 .and. icompress < 10 .and. product(dims_compress) > 1) then
             docompress = .true.
          endif

          ! Create a property list for compression/chunking/filters

          if (docompress) then
             call h5pcreate_f(H5P_DATASET_CREATE_F, propid, hdferr)
             call h5pset_chunk_f(propid, ndims, dims_compress, hdferr)
             call h5pset_shuffle_f(propid, hdferr)
             call h5pset_deflate_f(propid, icompress, hdferr)
          endif

       endif

       ! Create the global file space for the data if not previously created

       call h5screate_simple_f(ndims, dims_file, fspcid, hdferr)
       call h5dcreate_f(fileid, dname, stype, fspcid, dsetid, hdferr, dcpl_id=propid)
       call h5sclose_f(fspcid, hdferr)

       if (docompress) call h5pclose_f(propid, hdferr)

    endif

    ! Create the local memory space for the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! If we are only writing a subset of the local data,
    ! select the points in memory that we will writing from

    if (present(mcoords)) then

       if (size(mcoords) == 0) then

          call h5sselect_none_f(mspcid, hdferr)

       elseif (size(mcoords) < dims(ndims)) then

          if (ndims > 1) then
             offset(1:ndims-1) = 0
             countf(1:ndims-1) = dims(1:ndims-1)
          endif

          offset(ndims) = mcoords(1) - 1
          countf(ndims) = 1
          iop           = H5S_SELECT_SET_F

          ! Select a contiguous group of points as one slab for efficiency

          do j = 2, size(mcoords)
             if (mcoords(j) /= mcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(mspcid, iop, offset, countf, hdferr)
                offset(ndims) = mcoords(j) - 1
                countf(ndims) = 1
                iop           = H5S_SELECT_OR_F
             else
                countf(ndims) = countf(ndims) + 1
             endif
          enddo

          call h5sselect_hyperslab_f(mspcid, iop, offset, countf, hdferr)

       else

          call H5Sselect_all_f(mspcid, hdferr)

       endif

    else

       call H5Sselect_all_f(mspcid, hdferr)

    endif

    ! Create the dataspace in the file for the data

    call h5screate_simple_f(ndims, dims_file, dspcid, hdferr)

    ! For partial output, select the points in the output file that
    ! we will be writing to

    if (present(fcoords) .and. present(ifsize)) then

       if (size(fcoords) == 0) then

          call h5sselect_none_f(dspcid, hdferr)

       elseif (size(fcoords) < ifsize) then

          if (ndims > 1) then
             offset(1:ndims-1) = 0
             countf(1:ndims-1) = dims(1:ndims-1)
          endif

          offset(ndims) = fcoords(1) - 1
          countf(ndims) = 1
          iop           = H5S_SELECT_SET_F

          ! Select a contiguous group of points as one slab for efficiency

          do j = 2, size(fcoords)
             if (fcoords(j) /= fcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)
                offset(ndims) = fcoords(j) - 1
                countf(ndims) = 1
                iop           = H5S_SELECT_OR_F
             else
                countf(ndims) = countf(ndims) + 1
             endif
          enddo

          call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)

       else

          call H5Sselect_all_f(dspcid, hdferr)

       endif

    else

       call H5Sselect_all_f(dspcid, hdferr)

    endif

    ! For now, all parallel WRITES are done collectively
    if (dopario) then
       call h5pcopy_f(xferid_col, xferid, hdferr)
    endif

  end subroutine fh5_prepare_write

!===============================================================================

  subroutine fh5_write_int1(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dwrite_f, H5T_NATIVE_INTEGER_KIND
    import,                          only: i1, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer(i1), contiguous, target, intent(IN)  :: buf_integer1(..)  ! assumed-rank scalar or array
    integer,                         intent(OUT) :: hdferr
    type(c_ptr)                                  :: cptr

    cptr = c_loc(buf_integer1)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER_KIND(1), cptr, hdferr, &
                    mspcid, dspcid, xferid)

  end subroutine fh5_write_int1

!===============================================================================

  subroutine fh5_write_int2(buf_integer2, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dwrite_f, H5T_NATIVE_INTEGER_KIND
    import,                          only: i2, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer(i2), contiguous, target, intent(IN)  :: buf_integer2(..)  ! assumed-rank scalar or array
    integer,                         intent(OUT) :: hdferr
    type(c_ptr)                                  :: cptr

    cptr = c_loc(buf_integer2)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER_KIND(2), cptr, hdferr, &
                    mspcid, dspcid, xferid)

  end subroutine fh5_write_int2

!===============================================================================

  subroutine fh5_write_integer(buf_integer, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dwrite_f, H5T_NATIVE_INTEGER
    import,                          only: dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer, contiguous, target, intent(IN)  :: buf_integer(..)  ! assumed-rank scalar or array
    integer,                     intent(OUT) :: hdferr
    type(c_ptr)                              :: cptr

    cptr = c_loc(buf_integer)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, cptr, hdferr, &
                    mspcid, dspcid, xferid)

  end subroutine fh5_write_integer

!===============================================================================

  subroutine fh5_write_real(buf_real, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dwrite_f, H5T_NATIVE_REAL
    import,                          only: dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    real, contiguous, target, intent(IN)  :: buf_real(..)  ! assumed-rank scalar or array
    integer,                  intent(OUT) :: hdferr
    type(c_ptr)                           :: cptr

    cptr = c_loc(buf_real)

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, cptr, hdferr, &
                    mspcid, dspcid, xferid)

  end subroutine fh5_write_real

!===============================================================================

  subroutine fh5_write_real8(buf_real8, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dwrite_f, H5T_NATIVE_DOUBLE
    import,                          only: r8, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    real(r8), contiguous, target, intent(IN)  :: buf_real8(..)  ! assumed-rank scalar or array
    integer,                      intent(OUT) :: hdferr
    type(c_ptr)                               :: cptr

    cptr = c_loc(buf_real8)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, cptr, hdferr, &
                    mspcid, dspcid, xferid)

  end subroutine fh5_write_real8

!===============================================================================

  subroutine fh5_close_write(hdferr, no_dsetid)

    import, only: fh5_close_read

    implicit none (external, type)

    integer,           intent(OUT) :: hdferr
    logical, optional, intent(IN)  :: no_dsetid

    call fh5_close_read(hdferr, no_dsetid)

  end subroutine fh5_close_write

!===============================================================================

  subroutine fh5_open_dataset(dname, hdferr)

    use hdf5, only: h5iis_valid_f, H5Lexists_f, H5Dopen_f
    import,   only: dsetid, fileid

    implicit none (external, type)

    character(*), intent(IN ) :: dname
    integer,      intent(OUT) :: hdferr
    logical                   :: exists

    call h5iis_valid_f(dsetid, exists, hdferr)

    if (exists) then
       hdferr = -1
       write(*,*) "Error in fh5_open_dataset:"
       write(*,*) "Dataset id is already in use by HDF"
       return
    endif

    call H5Lexists_f(fileid, dname, exists, hdferr)

    if (.not. exists) hdferr = -1
    if (hdferr < 0) return

    call H5Dopen_f(fileid, dname, dsetid, hdferr)

  end subroutine fh5_open_dataset

!===============================================================================

  subroutine fh5_close_dataset(hdferr)

    use hdf5, only: h5iis_valid_f, h5dclose_f
    import,   only: dsetid

    implicit none (external, type)

    integer, intent(OUT) :: hdferr
    logical              :: valid

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (valid) call h5dclose_f(dsetid, hdferr)

  end subroutine fh5_close_dataset

!===============================================================================

  subroutine fh5_get_info(dname, ndims, dims)

    use hdf5, only: HSIZE_T, HID_T, h5iis_valid_f, H5Dget_space_f, &
                    H5Dclose_f, H5Sclose_f, H5Sget_simple_extent_dims_f
    import,   only: dsetid, fh5_open_dataset

    implicit none (external, type)

    character(*), intent(IN )             :: dname
    integer,      intent(OUT)             :: ndims
    integer,      intent(OUT), contiguous :: dims(:)

    logical                               :: valid
    integer(HSIZE_T)                      :: maxdimsc(7)
    integer(HSIZE_T)                      :: dimsc(7)
    integer(HID_T)                        :: dspcid
    integer                               :: hdferr

    ndims = -1
    dims  =  0

    call h5iis_valid_f(dsetid, valid, hdferr)

    if (.not. valid) then
       call fh5_open_dataset(dname, hdferr)
       if (hdferr < 0) return
    endif

    call H5Dget_space_f(dsetid, dspcid, hdferr)

    if (hdferr < 0) then
       if (.not. valid) call H5Dclose_f(dsetid, hdferr)
       return
    endif

    call H5Sget_simple_extent_dims_f(dspcid, dimsc, maxdimsc, ndims)

    if (ndims > 0) then
       dims(1:ndims)  = dimsc(1:ndims)
    endif

    call H5Sclose_f(dspcid, hdferr)
    if (.not. valid) call H5Dclose_f(dsetid, hdferr)

  end subroutine fh5_get_info

!===============================================================================

  subroutine fh5_get_chunk_dims(ndims,dims)

    use hdf5, only: HSIZE_T, HID_T, h5iis_valid_f, H5Dget_create_plist_f, &
                    H5Pget_layout_f, H5D_CHUNKED_F, H5Pget_chunk_f, H5Pclose_f
    import,   only: dsetid

    implicit none (external, type)

    integer, intent(IN)  :: ndims
    integer, intent(OUT) :: dims( ndims )

    integer(HSIZE_T)     :: dimsc( ndims )
    integer(HID_T)       :: creation_id
    integer              :: layout, rank, hdferr
    logical              :: valid

    dims   =  0
    layout = -1

    if (ndims <= 0) return

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (.not. valid) return

    call H5Dget_create_plist_f(dsetid, creation_id, hdferr)
    if (hdferr < 0) return

    call H5Pget_layout_f(creation_id, layout, hdferr)

    if (layout == H5D_CHUNKED_F) then
       call H5Pget_chunk_f(creation_id, ndims, dimsc, rank)
       rank = min(rank, ndims)
       if (rank > 0) dims(1:rank) = dimsc(1:rank)
    endif

    call H5Pclose_f(creation_id, hdferr)

  end subroutine fh5_get_chunk_dims

!===============================================================================

  subroutine fh5_read_int1(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dread_f, H5T_NATIVE_INTEGER_KIND
    import,                          only: i1, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer(i1), contiguous, target, intent(INOUT) :: buf_integer1(..)
    integer,                         intent(  OUT) :: hdferr
    type(c_ptr)                                    :: cptr

    cptr = c_loc(buf_integer1)

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER_KIND(1), cptr, hdferr, &
                   mspcid, dspcid, xferid)

  end subroutine fh5_read_int1


!===============================================================================

  subroutine fh5_read_int2(buf_integer2, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dread_f, H5T_NATIVE_INTEGER_KIND
    import,                          only: i2, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer(i2), contiguous, target, intent(INOUT) :: buf_integer2(..)
    integer,                         intent(  OUT) :: hdferr
    type(c_ptr)                                    :: cptr

    cptr = c_loc(buf_integer2)

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER_KIND(2), cptr, hdferr, &
                   mspcid, dspcid, xferid)

  end subroutine fh5_read_int2


!===============================================================================

  subroutine fh5_read_integer(buf_integer, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dread_f, H5T_NATIVE_INTEGER
    import,                          only: dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    integer, contiguous, target, intent(INOUT) :: buf_integer(..)
    integer,                     intent(  OUT) :: hdferr
    type(c_ptr)                                :: cptr

    cptr = c_loc(buf_integer)

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, cptr, hdferr, &
                   mspcid, dspcid, xferid)

  end subroutine fh5_read_integer


!===============================================================================

  subroutine fh5_read_real(buf_real, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dread_f, H5T_NATIVE_REAL
    import,                          only: dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    real, contiguous, target, intent(INOUT) :: buf_real(..)
    integer,                  intent(  OUT) :: hdferr
    type(c_ptr)                             :: cptr

    cptr = c_loc(buf_real)

    call h5dread_f(dsetid, H5T_NATIVE_REAL, cptr, hdferr, &
                   mspcid, dspcid, xferid)

  end subroutine fh5_read_real

!===============================================================================

  subroutine fh5_read_real8(buf_real8, hdferr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: h5dread_f, H5T_NATIVE_DOUBLE
    import,                          only: r8, dsetid, mspcid, dspcid, xferid

    implicit none (external, type)

    real(r8), contiguous, target, intent(INOUT) :: buf_real8(..)
    integer,                      intent(  OUT) :: hdferr
    type(c_ptr)                                 :: cptr

    cptr = c_loc(buf_real8)

    call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, cptr, hdferr, &
                   mspcid, dspcid, xferid)

  end subroutine fh5_read_real8

!===============================================================================

  subroutine fh5_prepare_read(dname, ndims, dims, hdferr, coords, start, counts)

    use misc_coms, only: io6
    use hdf5,      only: hsize_t, size_t, h5iis_valid_f, H5Lexists_f, h5fget_name_f, &
                         h5dopen_f, h5dget_space_f, h5sget_simple_extent_dims_f, &
                         h5screate_simple_f, h5pcopy_f, h5sselect_none_f, H5S_SELECT_SET_F, &
                         h5sselect_hyperslab_f, H5S_SELECT_OR_F
    import,        only: dimsf, ndimsf, fileid, dsetid, dspcid, mspcid, dopario, xferid, &
                         xferid_col, xferid_ind

    implicit none (external, type)

    character(*),      intent(IN)             :: dname
    integer,           intent(IN)             :: ndims
    integer,           intent(IN)             :: dims(ndims)
    integer,           intent(OUT)            :: hdferr

    integer, optional, intent(IN), contiguous :: coords(:)
    integer, optional, intent(IN), contiguous :: start(:)
    integer, optional, intent(IN), contiguous :: counts(:)

    integer          :: i, j, iop, ndims_file
    integer(HSIZE_T) :: dims_file(ndims), maxdims_file(ndims)
    integer(HSIZE_T) :: offset(ndims), countf(ndims)
    logical          :: exists
    character(30)    :: fname
    integer (SIZE_T) :: len

    dimsf  = 1
    ndimsf = ndims

    do i = 1, ndims
       dimsf(i) = dims(i)
    enddo

    if (present(coords)) then
       dimsf(ndims) = size(coords)
    else
       if (present(start) .and. present(counts)) then
          if (size(counts) >= ndims) dimsf(1:ndims) = counts(1:ndims)
       endif
    endif

    ! Is the HDF5 file open?

    call h5iis_valid_f(fileid, exists, hdferr)
    if (hdferr < 0 .or. .not. exists) then
       hdferr = -1
       write(io6,*) "fh5_prepare_read: Trying to read " // trim(dname) // &
                    " but no HDF5 file is open."
       return
    endif

    ! First make sure data exists in file

    call H5Lexists_f(fileid, dname, exists, hdferr)
    if (hdferr < 0 .or. .not. exists) then
       call h5fget_name_f(fileid, fname, len, hdferr)
       write(io6,*) "fh5_prepare_read: " // trim(dname) // " is not in file " // trim(fname)
       hdferr = -1
       return
    endif

    ! Make sure dsetid is not already open

    call h5iis_valid_f(dsetid, exists, hdferr)
    if (exists) then
       write(io6,*) "fh5_prepare_read: dsetid is already open"
       hdferr = -1
       return
    endif

    ! Open the dataset from file

    call h5dopen_f(fileid, dname, dsetid, hdferr)
    if (hdferr < 0) return

    ! Get the dataspace of the dataset (dimension, size, element type etc)

    call h5dget_space_f(dsetid, dspcid, hdferr)
    if (hdferr < 0) return

    ! get the dataspace dimensions in the file

    call h5sget_simple_extent_dims_f(dspcid, dims_file, maxdims_file, ndims_file)

    ! prepare the memspace to receive the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! If coords is present as an argument, use it to select the points
    ! in the file "dataspace" that we want to read

    if (present(coords)) then

       ! For now, parallel READS are only done collectively only if coords is present

       if (dopario) then
          call h5pcopy_f(xferid_col, xferid, hdferr)
       endif

       if (size(coords) == 0) then

          call h5sselect_none_f(dspcid, hdferr)

       else if (dimsf(ndims) < dims_file(ndims)) then

          ! If the number of points we want is less than that in the file,
          ! select the points we want (else we read the entire dataspace)

          if (ndims > 1) then
             offset(1:ndims-1) = 0
             countf(1:ndims-1) = dims(1:ndims-1)
          endif

          offset(ndims) = coords(1) - 1
          countf(ndims) = 1
          iop           = H5S_SELECT_SET_F

          ! Select a contiguous group of points as one slab for efficiency

          do j = 2, size(coords)
             if (coords(j) /= coords(j-1) + 1) then
                call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)
                offset(ndims) = coords(j) - 1
                countf(ndims) = 1
                iop           = H5S_SELECT_OR_F
             else
                countf(ndims) = countf(ndims) + 1
             endif
          enddo

          call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)

       endif

    else if (present(start) .and. present(counts)) then

       if (dopario) then
          call h5pcopy_f(xferid_ind, xferid, hdferr)
       endif

       ! If start and count are present, select the slab to read
       ! based on the start point (offset in file) and
       ! the count (number to read in each dimension)

       if (size(start) >= ndims .and. size(counts) >= ndims) then

          offset(1:ndims) = start (1:ndims) - 1
          countf(1:ndims) = counts(1:ndims)

          if ( all(offset >= 0) .and. all(countf >= 1) ) then
             call h5sselect_hyperslab_f(dspcid, H5S_SELECT_SET_F, offset, countf, hdferr)
          else
             write(*,*) "Error in fh5_prepare_read:"
             write(*,*) "Invalid start and count values"
          endif

       else
          write(*,*) "Error in fh5_prepare_read:"
          write(*,*) "Invalid start and count dimensions"
       endif

    else if (dopario) then

       call h5pcopy_f(xferid_ind, xferid, hdferr)

    endif

  end subroutine fh5_prepare_read

!===============================================================================

  subroutine fh5_close_read(hdferr, no_dsetid)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    use mem_para, only: myrank
    use hdf5,     only: H5Pget_mpio_actual_io_mode_f, H5D_MPIO_NO_COLLECTIVE_F, &
                        H5D_MPIO_CHUNK_INDEPENDENT_F, H5D_MPIO_CHUNK_COLLECTIVE_F, &
                        H5D_MPIO_CHUNK_MIXED_F, H5D_MPIO_CONTIG_COLLECTIVE_F
#endif
    use hdf5,     only: h5pclose_f, h5sclose_f, h5dclose_f
    import,       only: xferid, dsetid, mspcid, dspcid, dopario

    implicit none (external, type)

    integer,           intent(OUT) :: hdferr
    logical, optional, intent(IN)  :: no_dsetid
    logical                        :: close_dset

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    integer              :: mode

    if (.false. .and. dopario .and. myrank == 1) then
       call H5Pget_mpio_actual_io_mode_f(xferid, mode, hdferr)

       if ( mode == H5D_MPIO_NO_COLLECTIVE_F ) then
          write(*,*) " H5D_MPIO_NO_COLLECTIVE"
       elseif ( mode == H5D_MPIO_CHUNK_INDEPENDENT_F ) then
          write(*,*) " H5D_MPIO_CHUNK_INDEPENDENT"
       elseif ( mode == H5D_MPIO_CHUNK_COLLECTIVE_F ) then
          write(*,*) " H5D_MPIO_CHUNK_COLLECTIVE"
       elseif ( mode == H5D_MPIO_CHUNK_MIXED_F ) then
          write(*,*) " H5D_MPIO_CHUNK_MIXED"
       elseif ( mode == H5D_MPIO_CONTIG_COLLECTIVE_F ) then
          write(*,*) " H5D_MPIO_CONTIG_COLLECTIVE"
       else
          write(*,*) " Unknown H5D MPIO mode"
       endif
    endif
#endif

    if (dopario) then
       call h5pclose_f(xferid, hdferr)
    endif

    call h5sclose_f(mspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)

    close_dset = .true.
    if (present(no_dsetid)) close_dset = (.not. no_dsetid)

    if (close_dset) call h5dclose_f(dsetid, hdferr)

  end subroutine fh5_close_read

!===============================================================================

  subroutine fh5_prepare_write_ll(ndims, dims, dname, stype, hdferr, icompress, &
                                  fcoords, dims_chunk)

    use misc_coms, only: io6
    use hdf5,      only: hsize_t, hid_t, H5P_DEFAULT_F, h5iis_valid_f, H5Lexists_f, &
                         H5Dopen_f, H5Dget_space_status_f, H5D_SPACE_STS_ALLOCATED_F, &
                         H5D_SPACE_STS_PART_ALLOCATED_F, h5pcreate_f, h5pset_chunk_f, &
                         H5P_DATASET_CREATE_F, h5pset_shuffle_f, h5pset_deflate_f, &
                         h5screate_simple_f, h5dcreate_f, h5sclose_f, h5pclose_f, &
                         h5sselect_none_f, h5sselect_all_f, h5sselect_hyperslab_f, &
                         H5S_SELECT_SET_F, H5S_SELECT_OR_F, h5pcopy_f
    import,        only: fileid, dsetid, dspcid, mspcid, xferid, xferid_col, dimsf, &
                         ndimsf, dopario

    implicit none (external, type)

    integer,        intent(IN)                       :: ndims
    integer,        intent(IN)                       :: dims(ndims)
    character(*),   intent(IN)                       :: dname
    integer(HID_T), intent(IN)                       :: stype
    integer,        intent(OUT)                      :: hdferr
    integer,        intent(IN), optional             :: icompress
    integer,        intent(IN), optional, contiguous :: fcoords(:)
    integer,        intent(IN), optional             :: dims_chunk(ndims)

    logical          :: docompress, valid
    integer          :: i, j, iop, ilat, ilon, ilatp, ilonp
    integer(HSIZE_T) :: dims_file(ndims), dims_compress(ndims)
    integer(HSIZE_T) :: offset(ndims), countf(ndims)
    integer          :: ndims_file, flag
    integer(HID_T)   :: propid
    integer(HID_T)   :: fspcid

    ! Output dimensions, global file size

    ndims_file = ndims
    dims_file  = dims

    ! Local array memory size

    if (present(fcoords)) then

       ndimsf   = ndims - 1
       dimsf(1) = size(fcoords)
       do i = 3, ndims
          dimsf(i-1) = dims(i)
       end do

    else

       ndimsf = ndims - 1
       dimsf(1) = dims(1) * dims(2)
       do i = 3, ndims
          dimsf(i-1) = dims(i)
       end do

    endif

    ndimsf = max(ndimsf,1)

    hdferr     = 0
    propid     = H5P_DEFAULT_F
    docompress = .false.

   ! Is the HDF5 file open?

    call h5iis_valid_f(fileid, valid, hdferr)
    if (hdferr < 0 .or. .not. valid) then
       hdferr = -1
       write(io6,*) "fh5_prepare_write_ll: Trying to write " // trim(dname) // &
                    " but no HDF5 file is open."
       return
    endif

    ! Make sure dsetid is not already open

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (valid) then
       write(io6,*) "fh5_prepare_write_ll: dsetid is already open"
       hdferr = -1
       return
    endif

    call H5Lexists_f(fileid, dname, valid, hdferr)
    if (valid) then

       call H5Dopen_f(fileid, dname, dsetid, hdferr)
       call H5Dget_space_status_f(dsetid, flag, hdferr)

       if ( flag /= H5D_SPACE_STS_ALLOCATED_F .and. &
            flag /= H5D_SPACE_STS_PART_ALLOCATED_F ) then
          write(*,*) "HDF5 error: a data set has already been created for variable " // trim(dname)
          write(*,*) "in the HDF5 file but has no space allocated to it. Don't know what to do here."
          stop "HDF5 write error."
       endif

    else

       ! Setup HFS compression if we are not doing partial writes

       if ( present(icompress) .and. (.not. present(fcoords)) ) then

          if (present(dims_chunk)) then
             dims_compress = dims_chunk
          else
             dims_compress = dims
          endif

          if (icompress > 0 .and. icompress < 10 .and. product(dims_compress) > 1) then
             docompress = .true.
          endif

          ! Create a property list for compression/chunking/filters

          if (docompress) then
             call h5pcreate_f(H5P_DATASET_CREATE_F, propid, hdferr)
             call h5pset_chunk_f(propid, ndimsf, dimsf, hdferr)
             call h5pset_shuffle_f(propid, hdferr)
             call h5pset_deflate_f(propid, icompress, hdferr)
          endif

       endif

       ! Create the global file space for the data if not previously created

       call h5screate_simple_f(ndims, dims_file, fspcid, hdferr)
       call h5dcreate_f(fileid, dname, stype, fspcid, dsetid, hdferr, dcpl_id=propid)
       call h5sclose_f(fspcid, hdferr)

       if (docompress) call h5pclose_f(propid, hdferr)

    endif

    ! Create the local memory space for the data

    call h5screate_simple_f(ndimsf, dimsf, mspcid, hdferr)

    ! Select the points that we will output (either all or none for lat/lon output)

    if (present(fcoords)) then

       if (size(fcoords) == 0) then
          call h5sselect_none_f(mspcid, hdferr)
       else
          call h5sselect_all_f(mspcid, hdferr)
       endif

    endif

    ! Create the dataspace for the data

    call h5screate_simple_f(ndims_file, dims_file, dspcid, hdferr)

    ! For parallel output, select the points in the output file that
    ! we will be writing to

    if (present(fcoords)) then

       if (size(fcoords) == 0) then

          call h5sselect_none_f(dspcid, hdferr)

       else

          if (ndims > 2) then
             offset(3:ndims) = 0
             countf(3:ndims) = dims(3:ndims)
          endif

          ilon = mod(fcoords(1)-1,dims(1)) + 1
          ilat = (fcoords(1)-1)/dims(1) + 1

          offset(1) = ilon - 1
          offset(2) = ilat - 1
          countf(1) = 1
          countf(2) = 1
          iop       = H5S_SELECT_SET_F

          ! Select a contiguous group of points along a longitude as one slab for efficiency

          do j = 2, size(fcoords)

             ilonp = ilon
             ilatp = ilat

             ilon = mod(fcoords(j)-1,dims(1)) + 1
             ilat = (fcoords(j)-1)/dims(1) + 1

             if ((ilon /= ilonp + 1) .or. (ilat /= ilatp)) then
                call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)
                offset(1) = ilon - 1
                offset(2) = ilat - 1
                countf(1) = 1
                countf(2) = 1
                iop          = H5S_SELECT_OR_F
             else
                countf(1) = countf(1) + 1
             endif
          enddo

          call h5sselect_hyperslab_f(dspcid, iop, offset, countf, hdferr)

       endif

    endif

    ! For now, all parallel WRITES are done collectively
    if (dopario) then
       call h5pcopy_f(xferid_col, xferid, hdferr)
    endif

  end subroutine fh5_prepare_write_ll

!===============================================================================

  subroutine fh5f_write_attribute(name, ivalue, rvalue, dvalue, cvalue)

    use hdf5, only: hid_t, size_t, hsize_t, h5screate_f, H5S_SCALAR_F, H5Acreate_f, &
                    H5Awrite_f, H5T_STD_I32LE, H5T_NATIVE_INTEGER, H5T_IEEE_F32LE, &
                    H5T_NATIVE_REAL, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, h5tcopy_f, &
                    H5T_NATIVE_CHARACTER, h5tset_size_f, H5Tset_strpad_f, H5Tclose_f, &
                    H5Sclose_f, H5Aclose_f, H5T_STR_NULLTERM_F
    import,   only: r8, dsetid

    implicit none (external, type)

    character(*),           intent(in) :: name
    integer,      optional, intent(in) :: ivalue
    real,         optional, intent(in) :: rvalue
    real(r8),     optional, intent(in) :: dvalue
    character(*), optional, intent(in) :: cvalue

    integer(hid_t)   :: space_id, attr_id, atype_id
    integer(size_t)  :: attrlen
    integer(hsize_t) :: dims(1)
    integer          :: hdferr

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ! Create a scalar data space for the attribute
    CALL h5screate_f(H5S_SCALAR_F, space_id, hdferr)

    if     (present(ivalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, H5T_STD_I32LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)

    elseif (present(rvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, H5T_IEEE_F32LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)

    elseif (present(dvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, H5T_IEEE_F64LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, dvalue, dims, hdferr)

    elseif (present(cvalue)) then

       ! Create a type for fortran strings by expanding a native character
       attrlen = len_trim(cvalue) + 1

       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
       CALL h5tset_size_f(atype_id, attrlen, hdferr)

       ! Write null-terminated strings ("C-style")
       call H5Tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, atype_id, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, atype_id, trim(cvalue)//char(0), dims, hdferr)

       ! close character type
       CALL h5tclose_f(atype_id, hdferr)
    endif

    ! Close the attribute

    call H5Sclose_f(space_id, hdferr)
    call H5Aclose_f(attr_id,  hdferr)

  end subroutine fh5f_write_attribute

!===============================================================================

  subroutine fh5f_read_attribute(name, exists, ivalue, rvalue, dvalue, cvalue)

    use hdf5, only: hid_t, hsize_t, size_t, H5Aexists_f, H5Aopen_f, H5Aget_type_f, &
                    H5Tget_class_f, H5T_STRING_F, H5Aread_f, H5T_NATIVE_INTEGER, &
                    H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, H5T_NATIVE_CHARACTER, &
                    H5Tcopy_f, H5Tset_size_f, H5Aread_f, H5Tclose_f, H5Aclose_f
    import,   only: r8, dsetid

    implicit none (external, type)

    character(*),           intent(in ) :: name
    integer,      optional, intent(out) :: ivalue
    real,         optional, intent(out) :: rvalue
    real(r8),     optional, intent(out) :: dvalue
    character(*), optional, intent(out) :: cvalue
    logical,      optional, intent(out) :: exists

    integer(hid_t)   :: attr_id, atype_id, stype_id
    integer(hsize_t) :: dims(1)
    integer(size_t)  :: attrlen
    integer          :: class
    integer          :: ndims
    integer          :: hdferr
    logical          :: aexists

    if (present(exists)) exists = .false.

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ndims = 1
    dims  = 1

    call H5Aexists_f(dsetid, name, aexists, hdferr)

    if (.not. aexists) return

    call H5Aopen_f(dsetid, name, attr_id, hdferr)

    call H5Aget_type_f(attr_id, atype_id, hdferr)

    call H5Tget_class_f(atype_id, class, hdferr)

    if (present(cvalue) .and. (class /= H5T_STRING_F)) then

       write(*,*) "Input error reading attribute " // trim(name)
       write(*,*) "Character attribute expected but datatype is not a string in the file."
       write(*,*) "Not reading attribute"

    elseif ((.not. present(cvalue)) .and. (class == H5T_STRING_F)) then

       write(*,*) "Input error reading attribute " // trim(name)
       write(*,*) "Integer or float attribute expected but datatype is a string in the file."
       write(*,*) "Not reading attribute"

    else

       if (present(ivalue)) then

          ! read attribute from file
          ivalue = 0
          call H5Aread_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(rvalue)) then

          ! read attribute from file
          rvalue = 0.0
          call H5Aread_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(dvalue)) then

          ! read attribute from file
          dvalue = 0.0_r8
          call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, dvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(cvalue)) then

          attrlen = len(cvalue)
          call H5Tcopy_f(H5T_NATIVE_CHARACTER, stype_id, hdferr)
          call H5Tset_size_f(stype_id, attrlen, hdferr)

          ! read attribute from file
          cvalue = ' '
          call H5Aread_f(attr_id, stype_id, cvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

          ! Close the string type
          CALL H5Tclose_f(stype_id, hdferr)

       endif

    endif

    ! Close the inquired type
    CALL H5Tclose_f(atype_id, hdferr)

    ! Close the attribute
    call H5Aclose_f(attr_id,  hdferr)

  end subroutine fh5f_read_attribute

!===============================================================================

  subroutine fh5f_read_global_attribute(name, exists, ivalue, rvalue, dvalue, cvalue)

    use hdf5, only: hid_t, hsize_t, size_t, H5Gopen_f, H5Aexists_f, H5Gclose_f, &
                    H5Aopen_f, H5Aget_type_f, H5Tget_class_f, H5Aread_f, &
                    H5T_NATIVE_INTEGER, H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, &
                    H5T_NATIVE_CHARACTER, H5T_STRING_F, H5Tclose_f, H5Aclose_f, &
                    H5Gclose_f, H5Tcopy_f, H5TSET_SIZE_F
    import,   only: r8, fileid

    implicit none (external, type)

    character(*),           intent(in ) :: name
    integer,      optional, intent(out) :: ivalue
    real,         optional, intent(out) :: rvalue
    real(r8),     optional, intent(out) :: dvalue
    character(*), optional, intent(out) :: cvalue
    logical,      optional, intent(out) :: exists

    integer(hid_t)   :: attr_id, atype_id, stype_id, grp_id
    integer(hsize_t) :: dims(1)
    integer(size_t)  :: attrlen
    integer          :: class
    integer          :: ndims
    integer          :: hdferr
    logical          :: aexists

    if (present(exists)) exists = .false.

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ndims = 1
    dims  = 1

    ! Open the root group of the HDF5 file
    call H5Gopen_f(fileid, '/', grp_id, hdferr)

    call H5Aexists_f(grp_id, name, aexists, hdferr)

    if (.not. aexists) then
       call H5Gclose_f(grp_id, hdferr)
       return
    endif

    call H5Aopen_f(grp_id, name, attr_id, hdferr)

    call H5Aget_type_f(attr_id, atype_id, hdferr)

    call H5Tget_class_f(atype_id, class, hdferr)

    if (present(cvalue) .and. (class /= H5T_STRING_F)) then

       write(*,*) "Input error reading attribute " // trim(name)
       write(*,*) "Character attribute expected but datatype is not a string in the file."
       write(*,*) "Not reading attribute"

    elseif ((.not. present(cvalue)) .and. (class == H5T_STRING_F)) then

       write(*,*) "Input error reading attribute " // trim(name)
       write(*,*) "Integer or float attribute expected but datatype is a string in the file."
       write(*,*) "Not reading attribute"

    else

       if (present(ivalue)) then

          ! read attribute from file
          ivalue = 0
          call H5Aread_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(rvalue)) then

          ! read attribute from file
          rvalue = 0.0
          call H5Aread_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(dvalue)) then

          ! read attribute from file
          dvalue = 0.0_r8
          call H5Aread_f(attr_id, H5T_NATIVE_DOUBLE, dvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

       elseif (present(cvalue)) then

          attrlen = len(cvalue)
          call H5Tcopy_f(H5T_NATIVE_CHARACTER, stype_id, hdferr)
          call H5Tset_size_f(stype_id, attrlen, hdferr)

          ! read attribute from file
          cvalue = ' '
          call H5Aread_f(attr_id, stype_id, cvalue, dims, hdferr)
          if (present(exists)) exists = (hdferr >= 0)

          ! Close the string type
          CALL H5Tclose_f(stype_id, hdferr)

       endif

    endif

    ! Close the inquired type
    CALL H5Tclose_f(atype_id, hdferr)

    ! Close the attribute
    call H5Aclose_f(attr_id,  hdferr)

    ! Close the root group
    call H5Gclose_f(grp_id, hdferr)

  end subroutine fh5f_read_global_attribute

!===============================================================================

  subroutine fh5f_create_dim(dname, hdferr)

    use h5ds, only: H5DSset_scale_f
    import,   only: dsetid

    implicit none (external, type)

    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call H5DSset_scale_f(dsetid, hdferr, dimname=dname)

  end subroutine fh5f_create_dim

!===============================================================================

  subroutine fh5f_query_dimname(dname)

    use h5ds, only: H5DSis_scale_f, H5Dsget_scale_name_f
    use hdf5, only: size_t, h5iis_valid_f
    import,   only: dsetid

    implicit none (external, type)

    character(*), intent(OUT) :: dname

    integer(size_t)           :: dlen
    logical                   :: isscale, valid
    integer                   :: hdferr

    dname = ' '

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (.not. valid) then
       write(*,*) "No data set open in fh5f_query_dimname"
       return
    endif

    call H5DSis_scale_f(dsetid, isscale, hdferr)
    if (.not. isscale) return

    dlen = len(dname)
    call H5Dsget_scale_name_f(dsetid, dname, dlen, hdferr)

  end subroutine fh5f_query_dimname

!===============================================================================

  subroutine fh5f_attach_dims(ndims, dnames, hdferr)

    use h5ds, only: H5DSis_scale_f, H5Dsget_scale_name_f, H5DSattach_scale_f
    use hdf5, only: HID_T, H5Dopen_f, H5Dclose_f
    import,   only: fileid, dsetid

    implicit none (external, type)

    integer,      intent(IN)  :: ndims
    character(*), intent(IN)  :: dnames(ndims)
    integer,      intent(OUT) :: hdferr
    integer                   :: i, indx
    integer(HID_T)            :: dimid

    do i = 1, ndims
       indx = ndims - i + 1
       call H5Dopen_f(fileid, dnames(i), dimid, hdferr)
       call H5DSattach_scale_f(dsetid, dimid, indx, hdferr)
       call H5Dclose_f(dimid, hdferr)
    enddo

  end subroutine fh5f_attach_dims

!===============================================================================

  subroutine fh5f_write_global_attribute(name, ivalue, rvalue, dvalue, cvalue)

    use hdf5, only: hid_t, hsize_t, size_t, H5Gopen_f, h5screate_f, h5acreate_f, &
                    h5awrite_f, h5acreate_f, h5tcopy_f, h5tset_size_f, H5Tset_strpad_f, &
                    h5tclose_f, H5Sclose_f, H5Aclose_f, H5Gclose_f, h5s_scalar_f, &
                    h5t_std_i32le, h5t_native_integer, h5t_ieee_f32le, h5t_native_real, &
                    h5t_ieee_f64le, h5t_native_double, h5t_native_character, h5t_str_nullterm_f
    import,   only: r8, fileid

    implicit none (external, type)

    character(*),           intent(in) :: name
    integer,      optional, intent(in) :: ivalue
    real,         optional, intent(in) :: rvalue
    real(r8),     optional, intent(in) :: dvalue
    character(*), optional, intent(in) :: cvalue

    integer(hid_t)   :: space_id, attr_id, atype_id, grp_id
    integer(hsize_t) :: dims(1)
    integer(size_t)  :: attrlen
    integer          :: hdferr

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ! Open the root group of the HDF5 file
    call H5Gopen_f(fileid, '/', grp_id, hdferr)

    ! Create the data space for the attribute
    call h5screate_f(H5S_SCALAR_F, space_id, hdferr)

    if     (present(ivalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, H5T_STD_I32LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)

    elseif (present(rvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, H5T_IEEE_F32LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)

    elseif (present(dvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, H5T_IEEE_F64LE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, dvalue, dims, hdferr)

    elseif (present(cvalue)) then

       ! Create a type for fortran strings by expanding a native character
       attrlen = len_trim(cvalue) + 1

       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
       CALL h5tset_size_f(atype_id, attrlen, hdferr)

       ! Write null-terminated strings ("C-style")
       call H5Tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr)

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, atype_id, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, atype_id, trim(cvalue)//char(0), dims, hdferr)

       ! close character type
       CALL h5tclose_f(atype_id, hdferr)
    endif

    ! Close the attribute

    call H5Sclose_f(space_id, hdferr)
    call H5Aclose_f(attr_id,  hdferr)

    ! Close the root group
    call H5Gclose_f(grp_id, hdferr)

  end subroutine fh5f_write_global_attribute

!===============================================================================

  subroutine fh5_get_attached_scales(dnames)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
    use               hdf5,          only: HID_T, HSIZE_T, HVL_T, h5iis_valid_f, h5dget_space_f, &
                                           h5sget_simple_extent_ndims_f, H5Sclose_f, H5Aopen_f, &
                                           H5Aget_space_f, h5sget_simple_extent_dims_f, h5sclose_f, &
                                           h5aclose_f, H5Aget_type_f, h5aread_f, h5rget_name_f, &
                                           H5R_OBJECT_F
    use               h5ds,          only: H5DSget_num_scales_f
    import,                          only: dsetid

    implicit none (external, type)

    character(*), contiguous, intent(INOUT) :: dnames(:)
    integer                                 :: hdferr, nscal, rank, arank
    integer                                 :: j, i, d1, d2, is
    integer(HID_T)                          :: aid, asid, tid
    integer(HSIZE_T)                        :: adims(2), maxdims(2)
    character(30)                           :: dstring
    type(HVL_T),        allocatable, target :: ref_out(:)
    type(c_ptr)                             :: cptr
    integer(HID_T)                          :: lspcid
    logical                                 :: valid

    is = size(dnames)
    do i = 1, is
       dnames(i) = ' '
    enddo

    call h5iis_valid_f(dsetid, valid, hdferr)
    if (.not. valid) then
       write(*,*) "No data set open in fh5_get_attached_scales"
       return
    endif

    call H5Dget_space_f(dsetid, lspcid, hdferr)
    call h5sget_simple_extent_ndims_f(lspcid, rank, hdferr)
    call H5Sclose_f(lspcid, hdferr)

    if (size(dnames) < rank) then
       write(*,*) "Character array is not large enough to hold the dimension list."
       return
    endif

    do j = 1, rank
       call H5DSget_num_scales_f(dsetid, j, nscal, hdferr)
       if (nscal /= 1) return
    enddo

    call H5Aopen_f(dsetid, "DIMENSION_LIST", aid, hdferr)

    call H5Aget_space_f(aid, asid, hdferr)
    call h5sget_simple_extent_dims_f(asid, adims, maxdims, arank)
    call h5sclose_f(asid, hdferr)

    if (arank /= 1 .or. adims(1) /= rank) then
       write(*,*) "Dimension list does not match the dimensions of the variable in the file."
       call h5aclose_f(aid, hdferr)
       return
    endif

    call H5Aget_type_f(aid, tid, hdferr)

    ! This is a workaround to get the fortran interface to read the variable-length
    ! dimension scale references in the file

    allocate( ref_out( rank ) )
    cptr = c_loc(ref_out(1))
    call h5aread_f(aid, tid, cptr, hdferr)

    ! Now get the names of the dimensions referenced; we also need to reverse the order
    ! since they are stored in C order.

    do j = 1, rank
       i = rank - j + 1

       if (i <= is) then
          call h5rget_name_f(aid, H5R_OBJECT_F, ref_out(j)%p, dstring, hdferr)

          d1 = len(dnames(i))
          d2 = len_trim(dstring)

          if (dstring(1:1) == '/') then
             dnames(i) = dstring(2:)
          else
             dnames(i) = dstring
          endif
       endif
    enddo

    call h5aclose_f(aid, hdferr)

  end subroutine fh5_get_attached_scales

!===============================================================================

end module hdf5_f2f
