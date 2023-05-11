program test
  use hdf5_utils
  implicit none

  character(128) :: fileout="test.h5"
  integer        :: ndims, idims(3)
  integer        :: gdf_file_ver

  ndims = 0
  idims = 0

  print*,'Writing HDF5 file:',fileout

  call shdf5_open(fileout,'W',1)

  gdf_file_ver=3
  ndims=1 ; idims(1)=1
  call shdf5_orec(ndims,idims,'version',ivars=gdf_file_ver)

  call shdf5_close()

end program test
