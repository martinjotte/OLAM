module hdf5_f2f

  use hdf5
  use h5ds

  implicit none

  integer, parameter :: r8 = selected_real_kind(13,300)
  integer, parameter :: i1 = selected_int_kind(2)
  logical, parameter :: bigendian = ichar(transfer(1,'a')) == 0

  integer(HID_T)   :: fileid
  integer(HID_T)   :: propid
  integer(HID_T)   :: dsetid
  integer(HID_T)   :: fspcid
  integer(HID_T)   :: mspcid
  integer(HID_T)   :: dspcid
  integer(HSIZE_T) :: dimsf(7)
  integer          :: ndimsf

  interface fh5_write
     module procedure                 &
          fh5_write_int1_scalar,      &
          fh5_write_integer_scalar,   &
          fh5_write_real_scalar,      &
          fh5_write_character_scalar, &
          fh5_write_real8_scalar,     &
          fh5_write_logical_scalar,   &
          fh5_write_int1_array1,      &
          fh5_write_integer_array1,   &
          fh5_write_real_array1,      &
          fh5_write_character_array1, &
          fh5_write_real8_array1,     &
          fh5_write_logical_array1,   &
          fh5_write_int1_array2,      &
          fh5_write_integer_array2,   &
          fh5_write_real_array2,      &
          fh5_write_character_array2, &
          fh5_write_real8_array2,     &
          fh5_write_logical_array2,   &
          fh5_write_int1_array3,      &
          fh5_write_integer_array3,   &
          fh5_write_real_array3,      &
          fh5_write_character_array3, &
          fh5_write_real8_array3,     &
          fh5_write_logical_array3,   &
          fh5_write_int1_array4,      &
          fh5_write_integer_array4,   &
          fh5_write_real_array4,      &
          fh5_write_real8_array4
  end interface fh5_write

  interface fh5_read
     module procedure                 &
          fh5_read_int1_scalar,       &
          fh5_read_integer_scalar,    &
          fh5_read_real_scalar,       &
          fh5_read_character_scalar,  &
          fh5_read_real8_scalar,      &
          fh5_read_logical_scalar,    &
          fh5_read_int1_array1,       &
          fh5_read_integer_array1,    &
          fh5_read_real_array1,       &
          fh5_read_character_array1,  &
          fh5_read_real8_array1,      &
          fh5_read_logical_array1,    &
          fh5_read_int1_array2,       &
          fh5_read_integer_array2,    &
          fh5_read_real_array2,       &
          fh5_read_character_array2,  &
          fh5_read_real8_array2,      &
          fh5_read_logical_array2,    &
          fh5_read_int1_array3,       &
          fh5_read_integer_array3,    &
          fh5_read_real_array3,       &
          fh5_read_character_array3,  &
          fh5_read_real8_array3,      &
          fh5_read_logical_array3,    &
          fh5_read_int1_array4,       &
          fh5_read_integer_array4,    &
          fh5_read_real_array4,       &
          fh5_read_real8_array4
  end interface fh5_read

  private
  public :: fh5f_open, fh5f_close, fh5f_create, fh5_prepare_write, &
            fh5_write, fh5_close_write, fh5d_open, fh5s_get_ndims, &
            fh5s_get_dims, fh5d_close, fh5_prepare_read, &
            fh5_read, fh5_close_read, fh5f_write_attribute, &
            fh5f_create_dim, fh5f_attach_dims, fh5f_write_global_attribute

contains

!===============================================================================

  subroutine fh5f_open(locfn, iaccess, hdferr)
    implicit none

    character(*), intent(IN)  :: locfn
    integer,      intent(IN)  :: iaccess
    integer,      intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer        :: flags

    ! turn off default error handling
!    call h5eset_auto_f(0, hdferr)

    access_id = H5P_DEFAULT_F
    if(iaccess == 1) flags = H5F_ACC_RDONLY_F
    if(iaccess == 2) flags = H5F_ACC_RDWR_F

    call h5fopen_f(locfn, flags, fileid, hdferr, access_id)
  end subroutine fh5f_open

!=============================================================================

  subroutine fh5f_close(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5fclose_f(fileid, hdferr)
  end subroutine fh5f_close

!===============================================================================

  subroutine fh5f_create(locfn, iaccess, hdferr)
    implicit none

    character(*), intent(IN)  :: locfn
    integer,      intent(IN)  :: iaccess
    integer,      intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer(HID_T) :: create_id
    integer        :: flags
    integer        :: info

    create_id = H5P_DEFAULT_F
    access_id = H5P_DEFAULT_F
    if (iaccess == 1) flags = H5F_ACC_TRUNC_F
    if (iaccess == 2) flags = H5F_ACC_EXCL_F

    call h5fcreate_f(locfn, flags, fileid, hdferr, create_id, access_id)

  end subroutine fh5f_create

!===============================================================================

  subroutine fh5_prepare_write(ndims, dims, hdferr, icompress, &
                               mcoords, fcoords, ifsize, dims_chunk)
    implicit none

    integer, intent(IN)                       :: ndims
    integer, intent(IN)                       :: dims(ndims)
    integer, intent(OUT)                      :: hdferr
    integer, intent(IN),             optional :: icompress, ifsize
    integer, intent(IN), contiguous, optional :: mcoords(:), fcoords(:)
    integer, intent(IN),             optional :: dims_chunk(ndims)

    logical          :: docompress
    integer          :: i, j, iop
    integer(HSIZE_T) :: dims_file(ndims), dims_compress(ndims)
    integer(HSIZE_T) :: offset(ndims), countf(ndims)

    ! Output dimensions

    ndimsf         = ndims
    dimsf(1:ndims) = dims
    dims_file      = dims

    ! For distributed output, global data (file) size will be different then
    ! the local array size

    if (present(fcoords) .and. present(ifsize)) then
       dims_file(ndims) = ifsize
    endif

    hdferr = 0
    propid = H5P_DEFAULT_F

    ! Setup HFS compression

    if (present(icompress)) then

       if (icompress > 0 .and. icompress < 10) then
          docompress = .true.
       else
          docompress = .false.
       endif

       if (docompress) then

          if (present(dims_chunk)) then
             dims_compress = dims_chunk
          else
             dims_compress = dims
          endif

          if (product(dims_compress) > 1) then

             ! Create a property list for compression/chunking/filters

             call h5pcreate_f(H5P_DATASET_CREATE_F, propid, hdferr)
             call h5pset_chunk_f(propid, ndims, dims_compress, hdferr)
             call h5pset_shuffle_f(propid, hdferr)
             call h5pset_deflate_f(propid, icompress, hdferr)

          endif

       endif

    endif

    ! Create the global file space for the data

    call h5screate_simple_f(ndims, dims_file, fspcid, hdferr)

    ! Create the local memory space for the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! If we are only writing a subset of the local data,
    ! select the points that we will output

    if (present(mcoords)) then

       if (size(mcoords) == 0) then

          call h5sselect_none_f(mspcid, hdferr)

       else

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

       endif

    endif

    ! Create the dataspace for the data

    call h5screate_simple_f(ndims, dims_file, dspcid, hdferr)

    ! For parallel output, select the points in the output file that
    ! we will be writing to

    if (present(fcoords) .and. present(ifsize)) then

       if (size(fcoords) == 0) then

          call h5sselect_none_f(dspcid, hdferr)

       else

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

       endif

    endif

  end subroutine fh5_prepare_write

!===============================================================================

  subroutine fh5_write_int1_array1(buf_integer1, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(IN)  :: buf_integer1(:)
    character(*),                    intent(IN)  :: dname
    integer,                         intent(OUT) :: hdferr
    integer,               optional, intent(IN)  :: n

    logical                      :: create
    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1, dimsf, &
!                   hdferr, mspcid, dspcid

    cptr = c_loc(buf_integer1(1))
    call c_f_pointer(cptr, fptr, (/1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_int1_array1

!===============================================================================

  subroutine fh5_write_int1_array2(buf_integer1, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(IN)  :: buf_integer1(:,:)
    character(*),                    intent(IN)  :: dname
    integer,                         intent(OUT) :: hdferr
    integer,               optional, intent(IN)  :: n

    logical                      :: create
    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1(1,1), dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer1(1,1))
    call c_f_pointer(cptr, fptr, (/1,1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_int1_array2

!===============================================================================

  subroutine fh5_write_int1_array3(buf_integer1, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(IN)  :: buf_integer1(:,:,:)
    character(*),                    intent(IN)  :: dname
    integer,                         intent(OUT) :: hdferr
    integer,               optional, intent(IN)  :: n

    logical                      :: create
    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1(1,1,1), dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer1(1,1,1))
    call c_f_pointer(cptr, fptr, (/1,1,1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_int1_array3

!===============================================================================

  subroutine fh5_write_int1_array4(buf_integer1, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(IN)  :: buf_integer1(:,:,:,:)
    character(*),                    intent(IN)  :: dname
    integer,                         intent(OUT) :: hdferr
    integer,               optional, intent(IN)  :: n

    logical                      :: create
    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:,:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1(1,1,1,1), dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer1(1,1,1,1))
    call c_f_pointer(cptr, fptr, (/1,1,1,1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_int1_array4

!===============================================================================

  subroutine fh5_write_integer_array1(buf_integer, dname, hdferr, n)
    implicit none

    integer, contiguous, intent(IN)  :: buf_integer(:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_integer_array1

!===============================================================================

  subroutine fh5_write_integer_array2(buf_integer, dname, hdferr, n)
    implicit none

    integer, contiguous, intent(IN)  :: buf_integer(:,:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_integer_array2

!===============================================================================

  subroutine fh5_write_integer_array3(buf_integer, dname, hdferr, n)
    implicit none

    integer, contiguous, intent(IN)  :: buf_integer(:,:,:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_integer_array3

!===============================================================================

  subroutine fh5_write_integer_array4(buf_integer, dname, hdferr, n)
    implicit none

    integer, contiguous, intent(IN)  :: buf_integer(:,:,:,:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_integer_array4

!===============================================================================

  subroutine fh5_write_real_array1(buf_real, dname, hdferr, n)
    implicit none

    real,  contiguous, intent(IN)  :: buf_real(:)
    character(*),      intent(IN)  :: dname
    integer,           intent(OUT) :: hdferr
    integer, optional, intent(IN)  :: n
    logical                        :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real_array1

!===============================================================================

  subroutine fh5_write_real_array2(buf_real, dname, hdferr, n)
    implicit none

    real,  contiguous, intent(IN)  :: buf_real(:,:)
    character(*),      intent(IN)  :: dname
    integer,           intent(OUT) :: hdferr
    integer, optional, intent(IN)  :: n
    logical                        :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real_array2

!===============================================================================

  subroutine fh5_write_real_array3(buf_real, dname, hdferr, n)
    implicit none

    real,  contiguous, intent(IN)  :: buf_real(:,:,:)
    character(*),      intent(IN)  :: dname
    integer,           intent(OUT) :: hdferr
    integer, optional, intent(IN)  :: n
    logical                        :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real_array3

!===============================================================================

  subroutine fh5_write_real_array4(buf_real, dname, hdferr, n)
    implicit none

    real,  contiguous, intent(IN)  :: buf_real(:,:,:,:)
    character(*),      intent(IN)  :: dname
    integer,           intent(OUT) :: hdferr
    integer, optional, intent(IN)  :: n
    logical                        :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real_array4

!===============================================================================

  subroutine fh5_write_character_array1(buf_character, dname, hdferr, n)
    implicit none

    character(*), contiguous, intent(IN)  :: buf_character(:)
    character(*),             intent(IN)  :: dname
    integer,                  intent(OUT) :: hdferr
    integer,        optional, intent(IN)  :: n
    logical                               :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_character_array1

!===============================================================================

  subroutine fh5_write_character_array2(buf_character, dname, hdferr, n)
    implicit none

    character(*), contiguous, intent(IN)  :: buf_character(:,:)
    character(*),             intent(IN)  :: dname
    integer,                  intent(OUT) :: hdferr
    integer,        optional, intent(IN)  :: n
    logical                               :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_character_array2

!===============================================================================

  subroutine fh5_write_character_array3(buf_character, dname, hdferr, n)
    implicit none

    character(*), contiguous, intent(IN)  :: buf_character(:,:,:)
    character(*),             intent(IN)  :: dname
    integer,                  intent(OUT) :: hdferr
    integer,        optional, intent(IN)  :: n
    logical                               :: create

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_character_array3

!===============================================================================

  subroutine fh5_write_real8_array1(buf_real8, dname, hdferr, n)
    implicit none

    real(r8), contiguous, intent(IN)  :: buf_real8(:)
    character(*),         intent(IN)  :: dname
    integer,              intent(OUT) :: hdferr
    integer,    optional, intent(IN)  :: n
    logical                           :: create

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_REAL8_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real8_array1

!===============================================================================

  subroutine fh5_write_real8_array2(buf_real8, dname, hdferr, n)
    implicit none

    real(r8), contiguous, intent(IN)  :: buf_real8(:,:)
    character(*),         intent(IN)  :: dname
    integer,              intent(OUT) :: hdferr
    integer,    optional, intent(IN)  :: n
    logical                           :: create

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_REAL8_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real8_array2

!===============================================================================

  subroutine fh5_write_real8_array3(buf_real8, dname, hdferr, n)
    implicit none

    real(r8), contiguous, intent(IN)  :: buf_real8(:,:,:)
    character(*),         intent(IN)  :: dname
    integer,              intent(OUT) :: hdferr
    integer,    optional, intent(IN)  :: n
    logical                           :: create

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_REAL8_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real8_array3

!===============================================================================

  subroutine fh5_write_real8_array4(buf_real8, dname, hdferr, n)
    implicit none

    real(r8), contiguous, intent(IN)  :: buf_real8(:,:,:,:)
    character(*),         intent(IN)  :: dname
    integer,              intent(OUT) :: hdferr
    integer,    optional, intent(IN)  :: n
    logical                           :: create

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_REAL8_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    call h5dwrite_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real8_array4

!===============================================================================

  subroutine fh5_write_logical_array1(buf_logical, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    logical, contiguous, intent(IN)  :: buf_logical(:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    type(c_ptr)                      :: cptr
    integer, pointer, contiguous     :: fptr(:)

    integer(i1), allocatable, target :: buf_integer(:)
    integer                          :: i, size

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! converting logical to integer
    size = 1
    do i = 1, ndimsf
       size = size * dimsf(i)
    enddo

    allocate (buf_integer(size))
    buf_integer = 0

    do i = 1, size
       if(buf_logical(i)) buf_integer(i) = -1
    enddo

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer, dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer(1))
    call c_f_pointer(cptr, fptr, (/1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

    deallocate(buf_integer)

  end subroutine fh5_write_logical_array1

!===============================================================================

  subroutine fh5_write_logical_array2(buf_logical, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    logical, contiguous, intent(IN)  :: buf_logical(:,:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    type(c_ptr)                      :: cptr
    integer, pointer, contiguous     :: fptr(:)

    integer(i1), allocatable, target :: buf_integer(:,:)
    integer                          :: i, j, is, js

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! converting logical to integer
    is = size(buf_logical, dim=1)
    js = size(buf_logical, dim=2)

    allocate(buf_integer(is,js))
    buf_integer = 0

    do j = 1, js
       do i = 1, is
          if(buf_logical(i,j)) then
             buf_integer(i,j) = -1
          else
             buf_integer(i,j) = 0
          endif
       enddo
    enddo

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer, dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer(1,1))
    call c_f_pointer(cptr, fptr, (/1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

    deallocate(buf_integer)

  end subroutine fh5_write_logical_array2

!===============================================================================

  subroutine fh5_write_logical_array3(buf_logical, dname, hdferr, n)

    use, intrinsic :: iso_c_binding
    implicit none

    logical, contiguous, intent(IN)  :: buf_logical(:,:,:)
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr
    integer,   optional, intent(IN)  :: n
    logical                          :: create

    type(c_ptr)                      :: cptr
    integer, pointer, contiguous     :: fptr(:)

    integer(i1), allocatable, target :: buf_integer(:,:,:)
    integer                          :: i, j, k, is, js, ks

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    create = .true.
    if (present(n)) then
       if (n > 1) create = .false.
    endif

    if (create) then
       call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                        hdferr, propid)
    endif

    ! converting logical to integer
    is = size(buf_logical, dim=1)
    js = size(buf_logical, dim=2)
    ks = size(buf_logical, dim=3)

    allocate(buf_integer(is,js,ks))
    buf_integer = 0

    do k = 1, ks
       do j = 1, js
          do i = 1, is
             if(buf_logical(i,j,k)) then
                buf_integer(i,j,k) = -1
             else
                buf_integer(i,j,k) = 0
             endif
          enddo
       enddo
    enddo

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer, dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer(1,1,1))
    call c_f_pointer(cptr, fptr, (/1/))

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

    deallocate(buf_integer)

  end subroutine fh5_write_logical_array3

!===============================================================================

  subroutine fh5_write_int1_scalar(buf_integer1, dname, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, intent(IN)  :: buf_integer1
    character(*),        intent(IN)  :: dname
    integer,             intent(OUT) :: hdferr

    type(c_ptr)      :: cptr
    integer, pointer :: fptr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                     hdferr, propid)

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1, dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer1)
    call c_f_pointer(cptr, fptr)

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)


  end subroutine fh5_write_int1_scalar

!===============================================================================

  subroutine fh5_write_integer_scalar(buf_integer, dname, hdferr)
    implicit none

    integer,      intent(IN)  :: buf_integer
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, &
                     hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_integer_scalar

!===============================================================================

  subroutine fh5_write_real_scalar(buf_real, dname, hdferr)
    implicit none

    real,         intent(IN)  :: buf_real
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, &
                     hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real_scalar

!===============================================================================

  subroutine fh5_write_character_scalar(buf_character, dname, hdferr)
    implicit none

    character(*), intent(IN)  :: buf_character
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, &
                     hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_character_scalar

!===============================================================================

  subroutine fh5_write_real8_scalar(buf_real8, dname, hdferr)
    implicit none

    real(r8),     intent(IN)  :: buf_real8
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8, fortran real*8 is H5T_NATIVE_REAL_8, but
    ! in 1.10 and older is called H5T_NATIVE_DOUBLE
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dcreate_f(fileid, dname, FORTRAN_REAL8_TYPE, fspcid, dsetid, &
                     hdferr, propid)

    call h5dwrite_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_real8_scalar

!===============================================================================

  subroutine fh5_write_logical_scalar(buf_logical, dname, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    logical,         intent(IN)  :: buf_logical
    character(*),    intent(IN)  :: dname
    integer,         intent(OUT) :: hdferr

    integer(i1), target          :: buf_integer
    integer(HID_T)               :: FORTRAN_INT1_TYPE

    type(c_ptr)      :: cptr
    integer, pointer :: fptr

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    call h5dcreate_f(fileid, dname, FORTRAN_INT1_TYPE, fspcid, dsetid, &
                     hdferr, propid)

    buf_integer = 0
    if(buf_logical) buf_integer = -1

    ! HDF5 1.8 does not like passing integer*1
!   call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, buf_integer, dimsf, &
!                   hdferr, mspcid, dspcid)

    cptr = c_loc(buf_integer)
    call c_f_pointer(cptr, fptr)

    call h5dwrite_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                    hdferr, mspcid, dspcid)

  end subroutine fh5_write_logical_scalar

!===============================================================================

  subroutine fh5_close_write(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(fspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)
    call h5sclose_f(mspcid, hdferr)

    call h5dclose_f(dsetid, hdferr)
    call h5pclose_f(propid, hdferr)

  end subroutine fh5_close_write

!===============================================================================

  subroutine fh5d_open(dname, hdferr)
    implicit none

    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr
    logical                   :: exists

    call H5Lexists_f(fileid, dname, exists, hdferr)

    if (hdferr < 0 .or. .not. exists) then
       hdferr = -1
       return
    endif

    call h5dopen_f(fileid, dname, dsetid, hdferr)

    if(dsetid < 0 .or. hdferr < 0) then
       hdferr = -1
       return
    end if

    call h5dget_space_f(dsetid, dspcid, hdferr)

  end subroutine fh5d_open

!===============================================================================

  subroutine fh5s_get_ndims(ndims)
    implicit none

    integer, intent(OUT) :: ndims

    integer :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

  end subroutine fh5s_get_ndims

!===============================================================================

  subroutine fh5s_get_dims(dims)
    implicit none

    integer, intent(OUT), contiguous :: dims(:)
    integer(HSIZE_T)                 :: maxdimsc(7), dimsc(7)
    integer                          :: ndims, i
    integer                          :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

    call h5sget_simple_extent_dims_f(dspcid, dimsc, maxdimsc, hdferr)

    do i = 1, ndims
       dims(i) = int(dimsc(i))
    end do

  end subroutine fh5s_get_dims

!===============================================================================

  subroutine fh5d_close(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

  end subroutine fh5d_close

!===============================================================================

  subroutine fh5_read_int1_array1(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(INOUT) :: buf_integer1(:)
    integer,                         intent(  OUT) :: hdferr

    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dread_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1, dimsf, &
!                  hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    cptr = c_loc(buf_integer1(1))
    call c_f_pointer(cptr, fptr, (/1/))

    call h5dread_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_int1_array1

!===============================================================================

  subroutine fh5_read_int1_array2(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(INOUT) :: buf_integer1(:,:)
    integer,                         intent(  OUT) :: hdferr

    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dread_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1, dimsf, &
!                  hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    cptr = c_loc(buf_integer1(1,1))
    call c_f_pointer(cptr, fptr, (/1,1/))

    call h5dread_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_int1_array2

!===============================================================================

  subroutine fh5_read_int1_array3(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(INOUT) :: buf_integer1(:,:,:)
    integer,                         intent(  OUT) :: hdferr

    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dread_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1(1,1,1), dimsf, &
!                  hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    cptr = c_loc(buf_integer1(1,1,1))
    call c_f_pointer(cptr, fptr, (/1,1,1/))

    call h5dread_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_int1_array3

!===============================================================================

  subroutine fh5_read_int1_array4(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(i1), target, contiguous, intent(INOUT) :: buf_integer1(:,:,:,:)
    integer,                         intent(  OUT) :: hdferr

    type(c_ptr)                  :: cptr
    integer, pointer, contiguous :: fptr(:,:,:,:)

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dread_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1(1,1,1,1), dimsf, &
!                  hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    cptr = c_loc(buf_integer1(1,1,1,1))
    call c_f_pointer(cptr, fptr, (/1,1,1,1/))

    call h5dread_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_int1_array4

!===============================================================================

  subroutine fh5_read_integer_array1(buf_integer, hdferr)
    implicit none

    integer, contiguous, intent(INOUT) :: buf_integer(:)
    integer,             intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_integer_array1

!===============================================================================

  subroutine fh5_read_integer_array2(buf_integer, hdferr)
    implicit none

    integer, contiguous, intent(INOUT) :: buf_integer(:,:)
    integer,             intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_integer_array2

!===============================================================================

  subroutine fh5_read_integer_array3(buf_integer, hdferr)
    implicit none

    integer, contiguous, intent(INOUT) :: buf_integer(:,:,:)
    integer,             intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_integer_array3

!===============================================================================

  subroutine fh5_read_integer_array4(buf_integer, hdferr)
    implicit none

    integer, contiguous, intent(INOUT) :: buf_integer(:,:,:,:)
    integer,             intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_integer_array4

!===============================================================================

  subroutine fh5_read_real_array1(buf_real, hdferr)
    implicit none

    real, contiguous, intent(INOUT) :: buf_real(:)
    integer,          intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real_array1

!===============================================================================

  subroutine fh5_read_real_array2(buf_real, hdferr)
    implicit none

    real, contiguous, intent(INOUT) :: buf_real(:,:)
    integer,          intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real_array2

!===============================================================================

  subroutine fh5_read_real_array3(buf_real, hdferr)
    implicit none

    real, contiguous, intent(INOUT) :: buf_real(:,:,:)
    integer,          intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real_array3

!===============================================================================

  subroutine fh5_read_real_array4(buf_real, hdferr)
    implicit none

    real, contiguous, intent(INOUT) :: buf_real(:,:,:,:)
    integer,          intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real_array4

!===============================================================================

  subroutine fh5_read_character_array1(buf_character, hdferr)
    implicit none

    character, contiguous, intent(INOUT) :: buf_character(:)
    integer,               intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_character_array1

!===============================================================================

  subroutine fh5_read_character_array2(buf_character, hdferr)
    implicit none

    character, contiguous, intent(INOUT) :: buf_character(:,:)
    integer,               intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_character_array2

!===============================================================================

  subroutine fh5_read_character_array3(buf_character, hdferr)
    implicit none

    character, contiguous, intent(INOUT) :: buf_character(:,:,:)
    integer,               intent(  OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_character_array3

!===============================================================================

  subroutine fh5_read_real8_array1(buf_real8, hdferr)
    implicit none

    real(r8), contiguous, intent(INOUT) :: buf_real8(:)
    integer,              intent(OUT)   :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dread_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real8_array1

!===============================================================================

  subroutine fh5_read_real8_array2(buf_real8, hdferr)
    implicit none

    real(r8), contiguous, intent(INOUT) :: buf_real8(:,:)
    integer,              intent(  OUT) :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dread_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real8_array2

!===============================================================================

  subroutine fh5_read_real8_array3(buf_real8, hdferr)
    implicit none

    real(r8), contiguous, intent(INOUT) :: buf_real8(:,:,:)
    integer,              intent(  OUT) :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dread_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real8_array3

!===============================================================================

  subroutine fh5_read_real8_array4(buf_real8, hdferr)
    implicit none

    real(r8), contiguous, intent(INOUT) :: buf_real8(:,:,:,:)
    integer,              intent(  OUT) :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dread_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

  end subroutine fh5_read_real8_array4

!===============================================================================

  subroutine fh5_read_logical_array1(buf_logical, hdferr)
    implicit none

    logical, contiguous, intent(INOUT) :: buf_logical(:)
    integer,             intent(  OUT) :: hdferr

    integer, allocatable               :: buf_integer(:)
    integer                            :: i, size

    hdferr = 1

    size = 1
    do i = 1, ndimsf
       size = size * dimsf(i)
    enddo
    allocate (buf_integer(size))

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    ! converting integer to logical
    do i = 1, size
       if(buf_integer(i) == -1) then
          buf_logical(i) = .true.
       else
          buf_logical(i) = .false.
       endif
    enddo

    deallocate(buf_integer)

  end subroutine fh5_read_logical_array1

!===============================================================================

  subroutine fh5_read_logical_array2(buf_logical, hdferr)
    implicit none

    logical, contiguous, intent(INOUT) :: buf_logical(:,:)
    integer,             intent(  OUT) :: hdferr

    integer, allocatable               :: buf_integer(:,:)
    integer                            :: i, j, is, js

    hdferr = 1

    is = size(buf_logical,dim=1)
    js = size(buf_logical,dim=2)

    allocate(buf_integer(is,js))

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    ! converting integer to logical
    do j = 1, js
       do i = 1, is
          if(buf_integer(i,j) == -1) then
             buf_logical(i,j) = .true.
          else
             buf_logical(i,j) = .false.
          endif
       enddo
    enddo

    deallocate(buf_integer)

  end subroutine fh5_read_logical_array2

!===============================================================================

  subroutine fh5_read_logical_array3(buf_logical, hdferr)
    implicit none

    logical, contiguous, intent(INOUT) :: buf_logical(:,:,:)
    integer,             intent(  OUT) :: hdferr

    integer, allocatable               :: buf_integer(:,:,:)
    integer                            :: i, j, k, is, js, ks

    hdferr = 1

    is = size(buf_logical,dim=1)
    js = size(buf_logical,dim=2)
    ks = size(buf_logical,dim=3)

    allocate(buf_integer(is,js,ks))

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
                   hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    ! converting integer to logical
    do k = 1, ks
       do j = 1, js
          do i = 1, is
             if(buf_integer(i,j,k) == -1) then
                buf_logical(i,j,k) = .true.
             else
                buf_logical(i,j,k) = .false.
             endif
          enddo
       enddo
    enddo

    deallocate(buf_integer)

  end subroutine fh5_read_logical_array3

!===============================================================================

  subroutine fh5_read_int1_scalar(buf_integer1, hdferr)

    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr)      :: cptr
    integer, pointer :: fptr

    integer(i1), target, intent(INOUT) :: buf_integer1
    integer,             intent(  OUT) :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran integer*1 is H5T_NATIVE_INTEGER_1,
    ! but in hdf5 1.10 and newer it is H5T_NATIVE_INTEGER_KIND(1)
    integer(HID_T) :: FORTRAN_INT1_TYPE

    if (bigendian) then
       FORTRAN_INT1_TYPE = H5T_STD_I8BE
    else
       FORTRAN_INT1_TYPE = H5T_STD_I8LE
    endif

    ! HDF5 1.8 does not like passing integer*1
!   call h5dread_f(dsetid, FORTRAN_INT1_TYPE, buf_integer1, dimsf, hdferr)

    cptr = c_loc(buf_integer1)
    call c_f_pointer(cptr, fptr)

    call h5dread_f(dsetid, FORTRAN_INT1_TYPE, fptr, dimsf, hdferr)

  end subroutine fh5_read_int1_scalar

!===============================================================================

  subroutine fh5_read_integer_scalar(buf_integer, hdferr)
    implicit none

    integer, intent(INOUT) :: buf_integer
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

  end subroutine fh5_read_integer_scalar

!===============================================================================

  subroutine fh5_read_real_scalar(buf_real, hdferr)
    implicit none

    real,    intent(INOUT) :: buf_real
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr)

  end subroutine fh5_read_real_scalar

!===============================================================================

  subroutine fh5_read_character_scalar(buf_character, hdferr)
    implicit none

    character, intent(INOUT) :: buf_character
    integer,   intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr)

  end subroutine fh5_read_character_scalar

!===============================================================================

  subroutine fh5_read_real8_scalar(buf_real8, hdferr)
    implicit none

    real(r8), intent(INOUT) :: buf_real8
    integer,  intent(OUT)   :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    call h5dread_f(dsetid, FORTRAN_REAL8_TYPE, buf_real8, dimsf, hdferr)

  end subroutine fh5_read_real8_scalar

!===============================================================================

  subroutine fh5_read_logical_scalar(buf_logical, hdferr)
    implicit none

    logical, intent(INOUT) :: buf_logical
    integer, intent(OUT)   :: hdferr
    integer                :: buf_integer

    hdferr = 1

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

    buf_logical = .false.
    if(buf_integer == -1) buf_logical = .true.

  end subroutine fh5_read_logical_scalar

!===============================================================================

  subroutine fh5_prepare_read(dname, ndims, dims, hdferr, coords, &
                              start, counts)
    implicit none

    character(*),      intent(IN)             :: dname
    integer,           intent(IN)             :: ndims
    integer,           intent(IN), contiguous :: dims(:)
    integer,           intent(OUT)            :: hdferr

    integer, optional, intent(IN), contiguous :: coords(:)
    integer, optional, intent(IN), contiguous :: start(:)
    integer, optional, intent(IN), contiguous :: counts(:)

    integer          :: i, j, iop, ndims_file
    integer(HSIZE_T) :: dims_file(7), maxdims_file(7)
    integer(HSIZE_T) :: offset(ndims), countf(ndims)
    logical          :: exists

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

    ! First make sure data exists in file

    call H5Lexists_f(fileid, dname, exists, hdferr)
    if (hdferr < 0 .or. .not. exists) then
       hdferr = -1
       return
    endif

    ! Open the dataset from file

    call h5dopen_f(fileid, dname, dsetid, hdferr)
    if (hdferr < 0) return

    ! getting the dataspace of the dataset (dimension, size, element type etc)

    call h5dget_space_f(dsetid, dspcid, hdferr)
    if (hdferr < 0) return

    ! prepare the memspace to receive the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! If coords is present as an argument, use it to select the points
    ! in the file "dataspace" that we want to read

    if (present(coords)) then

       ! check the dataspace dimensions in the file

       call h5sget_simple_extent_ndims_f(dspcid, ndims_file, hdferr)
       call h5sget_simple_extent_dims_f(dspcid, dims_file, maxdims_file, hdferr)

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

    endif

  end subroutine fh5_prepare_read

!===============================================================================

  subroutine fh5_close_read(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(mspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

  end subroutine fh5_close_read

!===============================================================================

  subroutine fh5f_write_attribute(name, ivalue, rvalue, dvalue, cvalue)
    implicit none

    character(*),           intent(in) :: name
    integer,      optional, intent(in) :: ivalue
    real,         optional, intent(in) :: rvalue
    real(r8),     optional, intent(in) :: dvalue
    character(*), optional, intent(in) :: cvalue

    integer(hid_t)   :: space_id, attr_id, atype_id
    integer(hsize_t) :: dims(1)
    integer(size_t)  :: attrlen

    integer          :: ndims
    integer          :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ndims = 1
    dims  = 1

    ! Create the data space for the attribute
    call H5Screate_simple_f(ndims, dims, space_id, hdferr)

    if     (present(ivalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)

    elseif (present(rvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, H5T_NATIVE_REAL, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)

    elseif (present(dvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(dsetid, name, FORTRAN_REAL8_TYPE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, FORTRAN_REAL8_TYPE, dvalue, dims, hdferr)

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

  subroutine fh5f_create_dim(dname, hdferr)
    implicit none

    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call H5DSset_scale_f(dsetid, hdferr, dimname=dname)

  end subroutine fh5f_create_dim

!===============================================================================

 subroutine fh5f_attach_dims(ndims, dnames, hdferr)
    implicit none

    integer,      intent(IN)  :: ndims
    character(*), intent(IN)  :: dnames(ndims)
    integer,      intent(OUT) :: hdferr
    integer                   :: i, indx
    integer(HID_T)            :: dimid

    do i = 1, ndims
       indx = ndims - i + 1
       call h5dopen_f(fileid, dnames(i), dimid, hdferr)
       call H5DSattach_scale_f(dsetid, dimid, indx, hdferr)
       call h5dclose_f(dimid, hdferr)
    enddo

  end subroutine fh5f_attach_dims

!===============================================================================

  subroutine fh5f_write_global_attribute(name, ivalue, rvalue, dvalue, cvalue)
    implicit none

    character(*),           intent(in) :: name
    integer,      optional, intent(in) :: ivalue
    real,         optional, intent(in) :: rvalue
    real(r8),     optional, intent(in) :: dvalue
    character(*), optional, intent(in) :: cvalue

    integer(hid_t)   :: space_id, attr_id, atype_id, grp_id
    integer(hsize_t) :: dims(1)
    integer(size_t)  :: attrlen

    integer          :: ndims
    integer          :: hdferr

    ! This is to remain compatible with both HDF5 1.8 and 1.10.
    ! In hdf5 1.8 and older, fortran real*8 is H5T_NATIVE_REAL_8
    ! but this is removed in 1.10 and later
    integer(HID_T) :: FORTRAN_REAL8_TYPE

    if (bigendian) then
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64BE
    else
       FORTRAN_REAL8_TYPE = H5T_IEEE_F64LE
    endif

    if ( (.not. present(ivalue))  .and. &
         (.not. present(rvalue))  .and. &
         (.not. present(dvalue))  .and. &
         (.not. present(cvalue)) ) return

    ! Open the root group of the HDF5 file
    call H5Gopen_f(fileid, '/', grp_id, hdferr)

    ndims = 1
    dims  = 1

    ! Create the data space for the attribute
    call H5Screate_simple_f(ndims, dims, space_id, hdferr)

    if     (present(ivalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, ivalue, dims, hdferr)

    elseif (present(rvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, H5T_NATIVE_REAL, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, H5T_NATIVE_REAL, rvalue, dims, hdferr)

    elseif (present(dvalue)) then

       ! Create an attribute for this dataset
       call H5Acreate_f(grp_id, name, FORTRAN_REAL8_TYPE, space_id, attr_id, hdferr)

       ! write attribute to file
       call H5Awrite_f(attr_id, FORTRAN_REAL8_TYPE, dvalue, dims, hdferr)

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

end module hdf5_f2f
