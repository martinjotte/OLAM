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

  real, parameter :: soil_rough =  .02   ! soil roughness height
  real, parameter :: snow_rough =  .003  ! snowcover roughness height

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

  integer :: nvgcon   ! leaf class (optionally used OLAMIN parameter)

  real :: dt_leaf      ! leaf timestep [s]
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

  real, parameter :: thermcond_dry_organic = 0.05  ! [W/(m K)]
  real, parameter :: thermcond_sat_organic = 0.25  ! [W/(m K)]
  real, parameter :: thermcond_bedrock     = 3.0   ! [W/(m K)]
  real, parameter :: thermcond_liq         = 0.6   ! [W/(m K)]
  real, parameter :: thermcond_ice         = 2.29  ! [W/(m K)]
  real, parameter :: thermcond_firn        = 0.8   ! [W/(m K)] Based on firn density of 600 kg/m^3
                                                   !           and the following reference:
                             ! Oster, S.E. and M.R. Albert (2002): Thermal conductivity of polar firn,
                             ! Journal of Glaciology, Vol. 68, Issue 272.

  real, parameter :: specifheat_bedrock = 2.2e6   ! Volumetric heat capacity of deep ground [J/(m^3 K)]
  real, parameter :: specifheat_coarse  = 2.128e6 ! Volumetric heat capacity of sand [J/(m^3 K)]
  real, parameter :: specifheat_fine    = 2.385e6 ! Volumetric heat capacity of silt/clay [J/(m^3 K)]
  real, parameter :: specifheat_organic = 2.5e6   ! Volumetric heat capacity of organic matter [J/(m^3 K)]

End Module leaf_coms
