Module leaf_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

  use max_dims,    only: maxndvifiles, pathlen
  use consts_coms, only: r8

  implicit none

  private :: maxndvifiles, pathlen, r8

  !-----------------------------------------------------------------------
  ! VEGETATION, CANOPY, SURFACE WATER, AND NAMELIST VARIABLE SECTION
  !-----------------------------------------------------------------------

  ! Parameters for subroutine vegndvi

  real, parameter :: sr_min     = 1.081  ! Minimum simple ratio
  real, parameter :: fpar_min   =  .001  ! Minimum fpar
  real, parameter :: fpar_max   =  .950  ! Maximum fpar

  real, parameter :: soil_rough =  .050  ! soil roughness height
  real, parameter :: snow_rough =  .010  ! snowcover roughness height

  character(pathlen) :: veg_database
  character(pathlen) :: soil_database
  character(pathlen) :: soilgrids_database
  character(pathlen) :: glhymps_database
  character(pathlen) :: ndvi_database
  character(pathlen) :: watertab_db

  integer            :: nndvifiles, indvifile
  character(pathlen) :: fnames_ndvi  (maxndvifiles)
  character(14)      :: ctotdate_ndvi(maxndvifiles)
  real(r8)           :: s1900_ndvi   (maxndvifiles)

  integer :: isoilflg
  integer :: isoilptf
  integer :: ivegflg
  integer :: ndviflg
  integer :: isfcl
  integer :: iupdndvi

  integer :: isoilstateinit
  integer :: iwatertabflg

  integer :: nzs      ! maximum # of snowcover layers

  integer :: nvgcon   ! leaf class (optionally used OLAMIN parameter)

  real :: dt_leaf      ! leaf timestep [s]
  real :: snowmin_expl ! minimum surface water mass for explicit computation [kg/m^2]
  real :: wcap_min     ! minimum surface water [kg/m^2]
  real :: wcap_vmin    ! minimum vegetation surface water [kg/m^2]

! Leaf vegetation parameters

  integer, parameter :: nvtyp = 21    ! total # of leaf ("veg") classes

  integer :: kroot  (nvtyp)  ! k index of soil layer of deepest roots

  real :: albv_green(nvtyp)  ! green vegetation albedo
  real :: albv_brown(nvtyp)  ! brown vegetation albedo
  real :: emisv     (nvtyp)  ! vegetation infrared emissivity
  real :: sr_max    (nvtyp)  ! maximum simple ratio (for using NDVI)
  real :: tai_max   (nvtyp)  ! maximum vegetation total area index
  real :: sai       (nvtyp)  ! vegetation stem area index
  real :: veg_clump (nvtyp)  ! vegetation clumping factor
  real :: veg_frac  (nvtyp)  ! vegetation fractional coverage
  real :: veg_ht    (nvtyp)  ! vegetation height [m]
  real :: z_root    (nvtyp)  ! height (negative) of deepest roots [m]
  real :: glai_max  (nvtyp)  ! vegation maximum green leaf area index
  real :: dead_frac (nvtyp)  ! vegetation dead-material fraction
  real :: rcmin     (nvtyp)  ! vegetation minimum stomatal resistance [s/m]
  real :: dfpardsr  (nvtyp)  ! rate of change of fpar with simple ratio (for using NDVI)

  !-----------------------------------------------------------------------
  ! SOIL SECTION
  !-----------------------------------------------------------------------

  real, parameter :: emisg = 0.98      ! Ground surface longwave emissivity []
  real, parameter :: emisw = 1.0       ! Surface water longwave emissivity []

End Module leaf_coms

