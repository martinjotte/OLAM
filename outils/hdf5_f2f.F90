module hdf5_f2f

#ifdef OLAM_HDF5_FORTRAN

  use hdf5

  implicit none

  integer(HID_T) :: fileid
  integer(HID_T) :: mspcid
  integer(HID_T) :: dsetid
  integer(HID_T) :: dspcid
  integer(HSIZE_T), dimension(7) :: dimsf
  integer(HSIZE_T) :: ndimsf

  interface fh5_write
     module procedure fh5_write_integer_scalar,    &
          fh5_write_real_scalar,       &
          fh5_write_character_scalar,  &
          fh5_write_real8_scalar,      &
          fh5_write_logical_scalar,    &
          fh5_write_integer_array,    &
          fh5_write_real_array,       &
          fh5_write_character_array,  &
          fh5_write_real8_array,      &
          fh5_write_logical_array
  end interface fh5_write
  interface fh5d_read
     module procedure fh5d_read_integer_scalar,    &
          fh5d_read_real_scalar,       &
          fh5d_read_character_scalar,  &
          fh5d_read_real8_scalar,      &
          fh5d_read_logical_scalar,    &
          fh5d_read_integer_array,     &
          fh5d_read_real_array,        &
          fh5d_read_character_array,   &
          fh5d_read_real8_array,       &
          fh5d_read_logical_array
  end interface fh5d_read

contains

  subroutine fh5f_open(locfn, iaccess, hdferr)

    implicit none

    character(len=*), intent(IN) :: locfn
    integer, intent(IN) :: iaccess
    integer, intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer :: flags
    integer :: nulo

    ! turn off default error handling
    !call h5eset_auto_f(1,nulo)

    access_id = H5P_DEFAULT_F
    if(iaccess == 1) flags = H5F_ACC_RDONLY_F
    if(iaccess == 2) flags = H5F_ACC_RDWR_F

    call h5fopen_f(locfn, flags, fileid, hdferr, access_id)

    hdferr = int(fileid) + hdferr

    return

  end subroutine fh5f_open

  subroutine fh5f_close(hdferr)

    implicit none

    integer, intent(OUT) :: hdferr

    call h5fclose_f(fileid, hdferr)

    return

  end subroutine fh5f_close

  subroutine fh5f_create(locfn, iaccess, hdferr)

    implicit none

    character(len=*), intent(IN) :: locfn
    integer, intent(IN) :: iaccess
    integer, intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer(HID_T) :: create_id
    integer :: flags

    integer :: nulo

    create_id = H5P_DEFAULT_F
    access_id = H5P_DEFAULT_F
    if (iaccess == 1) flags = H5F_ACC_TRUNC_F
    if (iaccess == 2) flags = H5F_ACC_EXCL_F

    call h5fcreate_f(locfn, flags, fileid, hdferr, create_id, access_id)

    hdferr = int(fileid) + hdferr

    return

  end subroutine fh5f_create

  subroutine fh5_prepare_write(ndims, dims, hdferr)

    implicit none

    integer, intent(IN) :: ndims
    integer, dimension(*), intent(IN) :: dims
    integer, intent(OUT) :: hdferr

    integer :: i

    dimsf = 1
    ndimsf = ndims

    do i = 1, ndims
       dimsf(i) = dims(i)
    end do

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    return

  end subroutine fh5_prepare_write

  subroutine fh5_write_integer_array(h5type, buf_integer, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    integer, intent(IN) :: buf_integer(*)
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_integer_array

  subroutine fh5_write_real_array(h5type, buf_real, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real, intent(IN) :: buf_real(*)
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_real_array

  subroutine fh5_write_character_array(h5type, buf_character, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    character(len=*), intent(IN) :: buf_character(*)
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_character_array

  subroutine fh5_write_real8_array(h5type, buf_real8, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real(kind=8), intent(IN) :: buf_real8(*)
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_DOUBLE, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_real8_array

  subroutine fh5_write_logical_array(h5type, buf_logical, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    logical, intent(IN) :: buf_logical(*)
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    integer, allocatable :: buf_integer(:)
    integer :: i, size

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, mspcid, dsetid, hdferr)

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

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid)

    deallocate(buf_integer)

    return

  end subroutine fh5_write_logical_array

  subroutine fh5_write_integer_scalar(h5type, buf_integer, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    integer, intent(IN) :: buf_integer
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_integer_scalar

  subroutine fh5_write_real_scalar(h5type, buf_real, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real, intent(IN) :: buf_real
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_real_scalar

  subroutine fh5_write_character_scalar(h5type, buf_character, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    character(len=*), intent(IN) :: buf_character
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_character_scalar

  subroutine fh5_write_real8_scalar(h5type, buf_real8, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real(kind=8), intent(IN) :: buf_real8
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_DOUBLE, mspcid, dsetid, hdferr)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_real8_scalar

  subroutine fh5_write_logical_scalar(h5type, buf_logical, dname, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    logical, intent(IN) :: buf_logical
    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    integer :: buf_integer

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, mspcid, dsetid, hdferr)

    buf_integer = 0
    if(buf_logical) buf_integer = -1

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid)

    return

  end subroutine fh5_write_logical_scalar

  subroutine fh5_close_write(hdferr)

    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(mspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

    return

  end subroutine fh5_close_write

  subroutine fh5d_open(dname, hdferr)

    implicit none

    character(len=*), intent(IN) :: dname
    integer, intent(OUT) :: hdferr

    call h5dopen_f(fileid, dname, dsetid, hdferr)

    if(dsetid < 0) then
       hdferr = int(dsetid)
       return
    end if

    call h5dget_space_f(dsetid, dspcid, hdferr)

    hdferr = int(dspcid) + hdferr

    return

  end subroutine fh5d_open

  subroutine fh5s_get_ndims(ndims)

    implicit none

    integer, intent(OUT) :: ndims

    integer :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

    return

  end subroutine fh5s_get_ndims

  subroutine fh5s_get_dims(dims)

    implicit none

    integer, intent(OUT) :: dims(*)
    integer(HSIZE_T), dimension(*) :: maxdimsc(7), dimsc(7)

    integer :: ndims, i
    integer :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

    call h5sget_simple_extent_dims_f(dspcid, dimsc, maxdimsc, hdferr)

    do i = 1, ndims
       dims(i) = int(dimsc(i))
    end do

    return

  end subroutine fh5s_get_dims

  subroutine fh5d_close(hdferr)

    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

    return

  end subroutine fh5d_close

  subroutine fh5d_read_integer_array(h5type, buf_integer, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    integer, intent(INOUT) :: buf_integer(*)
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return

  end subroutine fh5d_read_integer_array

  subroutine fh5d_read_real_array(h5type, buf_real, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real, intent(INOUT) :: buf_real(*)
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return

  end subroutine fh5d_read_real_array

  subroutine fh5d_read_character_array(h5type, buf_character, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    character, intent(INOUT) :: buf_character(*)
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return

  end subroutine fh5d_read_character_array

  subroutine fh5d_read_real8_array(h5type, buf_real8, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real(KIND=8), intent(INOUT) :: buf_real8(*)
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return

  end subroutine fh5d_read_real8_array

  subroutine fh5d_read_logical_array(h5type, buf_logical, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    logical, intent(INOUT) :: buf_logical(*)
    integer, intent(OUT) :: hdferr

    integer, allocatable :: buf_integer(:)
    integer :: i, size

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

    return

  end subroutine fh5d_read_logical_array

  subroutine fh5d_read_integer_scalar(h5type, buf_integer, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    integer, intent(INOUT) :: buf_integer
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

    return

  end subroutine fh5d_read_integer_scalar

  subroutine fh5d_read_real_scalar(h5type, buf_real, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real, intent(INOUT) :: buf_real
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr)

    return

  end subroutine fh5d_read_real_scalar

  subroutine fh5d_read_character_scalar(h5type, buf_character, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    character, intent(INOUT) :: buf_character
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr)

    return

  end subroutine fh5d_read_character_scalar

  subroutine fh5d_read_real8_scalar(h5type, buf_real8, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    real(KIND=8), intent(INOUT) :: buf_real8
    integer, intent(OUT) :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr)

    return

  end subroutine fh5d_read_real8_scalar

  subroutine fh5d_read_logical_scalar(h5type, buf_logical, hdferr)

    implicit none

    integer, intent(IN) :: h5type
    logical, intent(INOUT) :: buf_logical
    integer, intent(OUT) :: hdferr

    integer :: buf_integer

    hdferr = 1

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

    buf_logical = .false.
    if(buf_integer == -1) buf_logical = .true.

    return

  end subroutine fh5d_read_logical_scalar

  subroutine fh5_prepare_read(dname, ndims, dims, hdferr, coords)

    implicit none

    character(len=*), intent(IN) :: dname
    integer, intent(IN) :: ndims
    integer, dimension(*), intent(IN) :: dims
    integer, intent(OUT) :: hdferr
    integer(HSIZE_T), allocatable, dimension(:,:) :: coordsf
    integer(HSIZE_T), allocatable, dimension(:,:) :: coordsfmem

    integer(HSIZE_T) :: num_elements, npoints1, npoints2

    integer :: i, j, k, ndims1, ndims2

    integer, optional :: coords(*)

    integer(HSIZE_T) :: dims1(7), dims2(7), maxdims1(7), maxdims2(7)

    dimsf = 1
    ndimsf = ndims

    do i = 1, ndims
       dimsf(i) = dims(i)          
    enddo

    ! opening the dataset from file

    call h5dopen_f(fileid, dname, dsetid, hdferr)
    if(dsetid < 0) then
       hdferr = int(dsetid)
       return
    end if

    ! getting the dataspace of the dataset (dimension, size, element type etc)

    call h5dget_space_f(dsetid, dspcid, hdferr)
    if(dspcid < 0) then
       hdferr = int(dspcid)
       return
    end if

    ! copying the dataspace to the memspace

    call h5scopy_f(dspcid, mspcid, hdferr)
    if(dspcid < 0) then
       hdferr = int(dspcid)
       return
    end if

    ! if coords is present select the points on dataspace and memspace

    if (present(coords)) then

       ! calculating the number of elements of the coordsf and coordsfmem
       ! arrays

       num_elements = 1
       do i = 1, ndims
          num_elements = num_elements * dimsf(i)
       enddo

       ! allocating the arrays, one for dataspace and one for memspace

       allocate(coordsf(1:ndims,1:num_elements))
       allocate(coordsfmem(1:ndims,1:num_elements))

       ! ensuring that all elements are declared
       ! (a sort of fall-back)

       coordsf = -1
       coordsfmem = -1

       if (ndims == 1) then

          ! the first element points to the right position

          coordsf(1,1) = 1
          coordsfmem(1,1) = 1

          ! with one dimension is simple

          do j = 2, dims(1)

             coordsf(1,j) = coords(j)
             coordsfmem(1,j) = j

          enddo

       else if (ndims == 2) then

          ! the first element points to the right position

          coordsf(1:ndims,1:dims(1)) = 1
          coordsfmem(1:ndims,1:dims(1)) = 1

          ! with two dimension is different: we have to do some aritmethic work

          do j = 2, dims(2)
             do k = 1, dims(1)

                coordsf(2,(j-1) * dims(1) + k) = coords(j)
                coordsf(1,(j-1) * dims(1) + k) = k

                coordsfmem(2,(j-1) * dims(1) + k) = j
                coordsfmem(1,(j-1) * dims(1) + k) = k

             end do
          end do

       else

          stop 'ndims > 2 using partial I/O is not implemented'

       endif

       ! selecting elements of the dataspace, according to coords

       call h5sselect_elements_f(dspcid, H5S_SELECT_SET_F, ndims, &
            num_elements, coordsf, hdferr)

       ! selecting elements of memspace (sequential)

       call h5sselect_elements_f(mspcid, H5S_SELECT_SET_F, ndims, &
            num_elements, coordsfmem, hdferr)

       deallocate(coordsf, coordsfmem)

    endif

    return

  end subroutine fh5_prepare_read

  subroutine fh5_close_read(hdferr)

    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(mspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

    return

  end subroutine fh5_close_read

#endif

end module hdf5_f2f

