! These routines are not in the HDF5 module before version 1.10.0
! These are only dummy routines so that the compiler won't complain
! about missing functions if you are compiling with an older HDF5
!===============================================================================

subroutine h5pset_coll_metadata_write_f(plist_id, is_collective, hdferr)

  use hdf5, only: HID_T
  implicit none

  integer(HID_T), intent(in)  :: plist_id
  logical,        intent(in)  :: is_collective
  integer,        intent(out) :: hdferr

  hdferr = 0
end subroutine h5pset_coll_metadata_write_f

!===============================================================================

subroutine h5pset_all_coll_metadata_ops_f(plist_id, is_collective, hdferr)

  use hdf5, only: HID_T
  implicit none

  integer(HID_T), intent(in)  :: plist_id
  logical,        intent(in)  :: is_collective
  integer,        intent(out) :: hdferr

  hdferr = 0
end subroutine h5pset_all_coll_metadata_ops_f

!===============================================================================

subroutine h5aread_cptr(attr_id, mem_type_id, buf, hdferr)

  use, intrinsic :: iso_c_binding, only: c_ptr
  use hdf5,                        only: HID_T

  implicit none

  interface
     integer function h5aread_f_c(attr_id, mem_type_id, buf) &
          bind(c, name='h5aread_f_c')
       import                     :: HID_T, c_ptr
       integer(HID_T), intent(in) :: attr_id
       integer(HID_T), intent(in) :: mem_type_id
       type(c_ptr),    value      :: buf
     end function h5aread_f_c
  end interface

  integer(HID_T),      intent(in)    :: attr_id     ! attribute identifier
  integer(HID_T),      intent(in)    :: mem_type_id ! memory datatype identifier
  type(c_ptr), target, intent(inout) :: buf
  integer,             intent(out)   :: hdferr      ! error code

  hdferr = h5aread_f_c(attr_id, mem_type_id, buf)

end subroutine h5aread_cptr

!===============================================================================

subroutine h5dread_cptr(dset_id, mem_type_id, buf, hdferr, &
                        mem_space_id, file_space_id, xfer_prp)

  use, intrinsic :: iso_c_binding, only: c_ptr
  use hdf5,                        only: HID_T, H5P_DEFAULT_F, H5S_ALL_F

  implicit none

  interface
     integer function h5dread_f_c(dset_id, mem_type_id, mem_space_id_default, &
          file_space_id_default, xfer_prp_default, buf) bind(c, name='h5dread_f_c')
       import                     :: HID_T, c_ptr
       integer(HID_T), intent(in) :: dset_id
       integer(HID_T), intent(in) :: mem_type_id
       integer(HID_T)             :: mem_space_id_default
       integer(HID_T)             :: file_space_id_default
       integer(HID_T)             :: xfer_prp_default
       type(c_ptr), value         :: buf
     end function h5dread_f_c
  end interface

  integer(HID_T),           intent(in)    :: dset_id       ! dataset identifier
  integer(HID_T),           intent(in)    :: mem_type_id   ! memory datatype identifier
  type(c_ptr),              intent(inout) :: buf
  integer,                  intent(out)   :: hdferr        ! error code
  integer(HID_T), optional, intent(in)    :: mem_space_id  ! memory dataspace identfier
  integer(HID_T), optional, intent(in)    :: file_space_id ! file dataspace identfier
  integer(HID_T), optional, intent(in)    :: xfer_prp      ! transfer property list identifier

  integer(HID_T) :: xfer_prp_default
  integer(HID_T) :: mem_space_id_default
  integer(HID_T) :: file_space_id_default

  xfer_prp_default      = H5P_DEFAULT_F
  mem_space_id_default  = H5S_ALL_F
  file_space_id_default = H5S_ALL_F

  if (present(xfer_prp))      xfer_prp_default      = xfer_prp
  if (present(mem_space_id))  mem_space_id_default  = mem_space_id
  if (present(file_space_id)) file_space_id_default = file_space_id

  hdferr = h5dread_f_c(dset_id, mem_type_id, mem_space_id_default, &
       file_space_id_default, xfer_prp_default, buf)

end subroutine h5dread_cptr

!===============================================================================

subroutine h5dwrite_cptr(dset_id, mem_type_id, buf, hdferr, mem_space_id, &
                         file_space_id, xfer_prp)

  use, intrinsic :: iso_c_binding, only: c_ptr
  use hdf5,                        only: HID_T, H5P_DEFAULT_F, H5S_ALL_F

  implicit none

  interface
     integer function h5dwrite_f_c(dset_id, mem_type_id, mem_space_id_default, &
          file_space_id_default, xfer_prp_default, buf ) bind(c, name='h5dwrite_f_c')
       import                     :: HID_T, c_ptr
       integer(HID_T), intent(in) :: dset_id
       integer(HID_T), intent(in) :: mem_type_id
       integer(HID_T)             :: mem_space_id_default
       integer(HID_T)             :: file_space_id_default
       integer(HID_T)             :: xfer_prp_default
       type(c_ptr), value         :: buf
     end function h5dwrite_f_c
  end interface

  integer(HID_T),           intent(in)  :: dset_id       ! dataset identifier
  integer(HID_T),           intent(in)  :: mem_type_id   ! memory datatype identifier
  type(c_ptr),              intent(in)  :: buf
  integer,                  intent(out) :: hdferr        ! error code
  integer(HID_T), optional, intent(in)  :: mem_space_id  ! memory dataspace identfier
  integer(HID_T), optional, intent(in)  :: file_space_id ! file dataspace identfier
  integer(HID_T), optional, intent(in)  :: xfer_prp      ! transfer property list identifier

  integer(HID_T) :: xfer_prp_default
  integer(HID_T) :: mem_space_id_default
  integer(HID_T) :: file_space_id_default

  xfer_prp_default      = H5P_DEFAULT_F
  mem_space_id_default  = H5S_ALL_F
  file_space_id_default = H5S_ALL_F

  if (present(xfer_prp))      xfer_prp_default      = xfer_prp
  if (present(mem_space_id))  mem_space_id_default  = mem_space_id
  if (present(file_space_id)) file_space_id_default = file_space_id

  hdferr = h5dwrite_f_c(dset_id, mem_type_id, mem_space_id_default, &
       file_space_id_default, xfer_prp_default, buf)

end subroutine h5dwrite_cptr

!===============================================================================

subroutine h5rget_name_cptr(loc_id, ref_type, ref, name, hdferr, size)

  use, intrinsic :: iso_c_binding, only: c_ptr, c_char
  use hdf5,                        only: HID_T, SIZE_T

  implicit none

  interface
     integer function h5rget_name_ptr_c(loc_id, ref_type, ref, name, name_len, size_default) &
          bind(c, name='h5rget_name_ptr_c')
       import :: c_char, c_ptr, HID_T, SIZE_T
       integer(HID_T), intent(in) :: loc_id
       integer, intent(in) :: ref_type
       type(c_ptr), value, intent(in) :: ref
       character(kind=c_char), dimension(*), intent(in) :: name
       integer(SIZE_T) :: name_len
       integer(SIZE_T) :: size_default
     end function h5rget_name_ptr_c
  end interface

  integer(HID_T), intent(in) :: loc_id
  integer, intent(in) :: ref_type
  type(c_ptr), intent(in) :: ref
  character(*), intent(inout) :: name
  integer, intent(out) :: hdferr
  integer(SIZE_T), optional, intent(out) :: size

  integer(SIZE_T) :: size_default
  integer(SIZE_T) :: name_len

  name_len = len(name)
  hdferr = h5rget_name_ptr_c(loc_id, ref_type, ref, name, name_len, size_default)
  if (present(size)) size = size_default

end subroutine h5rget_name_cptr
