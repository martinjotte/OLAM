Module max_dims

  implicit none

  integer, parameter :: maxgrds      = 20   ! Max # of grids
  integer, parameter :: maxremote    = 30   ! Max # of remote send/recv processes
  integer, parameter :: maxsndg      = 200  ! Max # of vertical levels for the input sounding
  integer, parameter :: maxndvifiles = 2000 ! Max # of input NDVI files
  integer, parameter :: maxisdirs    = 30   ! Max # of directories that contain data files
  integer, parameter :: maxnplt      = 150  ! Max # of fields to plot
  integer, parameter :: maxpltfiles  = 2000 ! Max # of input files for a plotonly run
  integer, parameter :: maxngrdll    = 20   ! Max # of geog. pts for each grid refinement
  integer, parameter :: pathlen      = 256  ! Max length of character strings for file paths
  integer, parameter :: maxnlspoly   = 7    ! Max # of M pts for a single land/sea cell
  integer, parameter :: maxlite      = 150  ! Max # of output "lite" variables
  integer, parameter :: maxlatlon    = 150  ! Max # of output "latlon" variables

End Module max_dims

