Module isan_coms

  use max_dims,    only: maxisdirs, pathlen
  use consts_coms, only: r8

  implicit none

  integer :: ioflgisz, ioflgvar, iszstage, ivrstage
  integer :: iyear, imonth, idate, ihour, ipoffset
  integer :: npd,lzon_bot,kzonoff

  integer :: nfgfiles = -1
  integer :: ifgfile  =  0

  character(pathlen) :: iapr(maxisdirs)
  character(pathlen) :: innpr

  character(pathlen), allocatable :: fnames_fg  (:)
  character(14),      allocatable :: ctotdate_fg(:)
  real(r8),           allocatable :: s1900_fg   (:)

  logical :: haso3
  integer :: nbot_o3

  ! Only used on MPI node 0:
  character(6) :: o3name = " "

! Pressure header variables:
  integer :: isversion, iyy, imm, idd, ihh, itinc, inproj, ivertcoord
  real    :: xnelat, xnelon, cntlat, secondlat

  integer :: nprx, npry, nprz
  integer, allocatable :: levpr(:)

  real    :: xswlon, xswlat, gdatdx, gdatdy
  real, allocatable :: pnpr(:)

  real, allocatable :: glat(:), plat(:)

  real, allocatable :: pcol_p(:)
  real, allocatable :: pcol_z(:,:)

  ! input variables interpolated to height at each OLAM column
  real, allocatable :: o_press (:,:)
  real, allocatable :: o_rho   (:,:)

  real, allocatable :: o_theta (:,:)
  real, allocatable :: o_rrw   (:,:)
  real, allocatable :: o_uzonal(:,:)
  real, allocatable :: o_umerid(:,:)
  real, allocatable :: o_ozone (:,:)

End module isan_coms
