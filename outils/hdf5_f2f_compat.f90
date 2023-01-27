! These routines are not in the HDF5 module before version 1.10.0
! These are only dummy routines so that the compiler won't complain
! about missing functions if you are compiling with an older HDF5

subroutine h5pset_coll_metadata_write_f(plist_id, is_collective, hdferr)

  use hdf5, only: hid_t
  implicit none

  integer(hid_t), intent(in)  :: plist_id
  logical,        intent(in)  :: is_collective
  integer,        intent(out) :: hdferr

  hdferr = 0
end subroutine h5pset_coll_metadata_write_f



subroutine h5pset_all_coll_metadata_ops_f(plist_id, is_collective, hdferr)

  use hdf5, only: hid_t
  implicit none

  integer(hid_t), intent(in)  :: plist_id
  logical,        intent(in)  :: is_collective
  integer,        intent(out) :: hdferr

  hdferr = 0
end subroutine h5pset_all_coll_metadata_ops_f
