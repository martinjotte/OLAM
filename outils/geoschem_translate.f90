Module gchem_aeros

  use consts_coms, only: r8

  logical, parameter :: anthOff = .true.  ! If true, use anthOff set of files
!!logical, parameter :: anthOff = .false. ! If false, use anthOn set of files

  integer :: ngchemfiles, igchemfile

  real(r8), allocatable :: s1900_gchem(:)

  Type gchem_vars

     ! Dimensions for GeosCHEM aerosol fields at a single time

     integer :: nlon
     integer :: nlat
     integer :: nlev

     ! Geoschem arrays (aerosols on lat-lon-press grid)

     real, allocatable :: dust1(:,:,:) ! dust1
     real, allocatable :: dust2(:,:,:) ! dust2
     real, allocatable :: dust3(:,:,:) ! dust3
     real, allocatable :: dust4(:,:,:) ! dust4
     real, allocatable :: asalt(:,:,:) ! accumulation mode salt
     real, allocatable :: csalt(:,:,:) ! coarse mode salt
     real, allocatable :: so4  (:,:,:) ! so4
     real, allocatable :: no3  (:,:,:) ! no3
     real, allocatable :: nh4  (:,:,:) ! nh4
     real, allocatable :: bcpi (:,:,:) ! hydrophillic black carbon
     real, allocatable :: bcpo (:,:,:) ! hydrophobic black carbon
     real, allocatable :: ocpi (:,:,:) ! hydrophillic organic carbon
     real, allocatable :: ocpo (:,:,:) ! hydrophobic organic carbon
     real, allocatable :: presg(:,:,:) ! pressure
     real, allocatable :: tempg(:,:,:) ! temperature

     real, allocatable :: lat(:) ! latitude
     real, allocatable :: lon(:) ! longitude

     real, allocatable :: readvar(:,:,:,:) ! 4D array for reading Kodros files

     character(len=80), allocatable :: fnames(:)

  End type gchem_vars

  type (gchem_vars) :: gchem(2)

  ! Nudging array (containing geoschem fields interpolated to OLAM grid)

  real, allocatable :: geoschem_nud(:,:,:) ! Dimensions (mza,mwa,nccntyp) 

Contains

  subroutine geoschem_init()

  use misc_coms, only: io6, s1900_sim

  implicit none

  ! Assign array dimensions to global (GL) and North America (NA) geoschem files

  gchem(1)%nlon = 144
  gchem(1)%nlat = 91
  gchem(1)%nlev = 30

  gchem(2)%nlon = 225
  gchem(2)%nlat = 202
  gchem(2)%nlev = 30

  ! Assign number of geoschem files for this simulation (hardwired value)

  ngchemfiles = 36

  allocate (s1900_gchem (ngchemfiles))

  allocate (gchem(1)%fnames(ngchemfiles))
  allocate (gchem(2)%fnames(ngchemfiles))

  ! Fill fnames_gchem array with geoschem file names - including relative or absolute path

  IF (anthOff) THEN  ! Use this for anthOff set of files

  gchem(1)%fnames(1) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_23_0000.h5'
  gchem(2)%fnames(1) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_23_0000.h5'

  gchem(1)%fnames(2) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_23_0600.h5'
  gchem(2)%fnames(2) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_23_0600.h5'

  gchem(1)%fnames(3) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_23_1200.h5'
  gchem(2)%fnames(3) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_23_1200.h5'

  gchem(1)%fnames(4) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_23_1800.h5'
  gchem(2)%fnames(4) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_23_1800.h5'

  gchem(1)%fnames(5) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_24_0000.h5'
  gchem(2)%fnames(5) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_24_0000.h5'

  gchem(1)%fnames(6) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_24_0600.h5'
  gchem(2)%fnames(6) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_24_0600.h5'

  gchem(1)%fnames(7) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_24_1200.h5'
  gchem(2)%fnames(7) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_24_1200.h5'

  gchem(1)%fnames(8) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_24_1800.h5'
  gchem(2)%fnames(8) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_24_1800.h5'

  gchem(1)%fnames(9) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_25_0000.h5'
  gchem(2)%fnames(9) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_25_0000.h5'

  gchem(1)%fnames(10) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_25_0600.h5'
  gchem(2)%fnames(10) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_25_0600.h5'

  gchem(1)%fnames(11) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_25_1200.h5'
  gchem(2)%fnames(11) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_25_1200.h5'

  gchem(1)%fnames(12) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_25_1800.h5'
  gchem(2)%fnames(12) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_25_1800.h5'

  gchem(1)%fnames(13) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_26_0000.h5'
  gchem(2)%fnames(13) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_26_0000.h5'

  gchem(1)%fnames(14) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_26_0600.h5'
  gchem(2)%fnames(14) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_26_0600.h5'

  gchem(1)%fnames(15) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_26_1200.h5'
  gchem(2)%fnames(15) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_26_1200.h5'

  gchem(1)%fnames(16) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_26_1800.h5'
  gchem(2)%fnames(16) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_26_1800.h5'

  gchem(1)%fnames(17) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_27_0000.h5'
  gchem(2)%fnames(17) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_27_0000.h5'

  gchem(1)%fnames(18) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_27_0600.h5'
  gchem(2)%fnames(18) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_27_0600.h5'

  gchem(1)%fnames(19) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_27_1200.h5'
  gchem(2)%fnames(19) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_27_1200.h5'

  gchem(1)%fnames(20) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_27_1800.h5'
  gchem(2)%fnames(20) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_27_1800.h5'

  gchem(1)%fnames(21) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_28_0000.h5'
  gchem(2)%fnames(21) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_28_0000.h5'

  gchem(1)%fnames(22) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_28_0600.h5'
  gchem(2)%fnames(22) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_28_0600.h5'

  gchem(1)%fnames(23) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_28_1200.h5'
  gchem(2)%fnames(23) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_28_1200.h5'

  gchem(1)%fnames(24) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_28_1800.h5'
  gchem(2)%fnames(24) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_28_1800.h5'

  gchem(1)%fnames(25) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_29_0000.h5'
  gchem(2)%fnames(25) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_29_0000.h5'

  gchem(1)%fnames(26) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_29_0600.h5'
  gchem(2)%fnames(26) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_29_0600.h5'

  gchem(1)%fnames(27) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_29_1200.h5'
  gchem(2)%fnames(27) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_29_1200.h5'

  gchem(1)%fnames(28) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_29_1800.h5'
  gchem(2)%fnames(28) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_29_1800.h5'

  gchem(1)%fnames(29) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_30_0000.h5'
  gchem(2)%fnames(29) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_30_0000.h5'

  gchem(1)%fnames(30) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_30_0600.h5'
  gchem(2)%fnames(30) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_30_0600.h5'

  gchem(1)%fnames(31) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_30_1200.h5'
  gchem(2)%fnames(31) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_30_1200.h5'

  gchem(1)%fnames(32) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_30_1800.h5'
  gchem(2)%fnames(32) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_30_1800.h5'

  gchem(1)%fnames(33) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_31_0000.h5'
  gchem(2)%fnames(33) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_31_0000.h5'

  gchem(1)%fnames(34) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_31_0600.h5'
  gchem(2)%fnames(34) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_31_0600.h5'

  gchem(1)%fnames(35) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_31_1200.h5'
  gchem(2)%fnames(35) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_31_1200.h5'

  gchem(1)%fnames(36) = '../olam_geoschem_files/anthOff_GL_Harvey_2017_08_31_1800.h5'
  gchem(2)%fnames(36) = '../olam_geoschem_files/anthOff_NA_Harvey_2017_08_31_1800.h5'

  ELSE  ! Use this for anthOn set of files

  gchem(1)%fnames(1) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_23_0000.h5'
  gchem(2)%fnames(1) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_23_0000.h5'

  gchem(1)%fnames(2) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_23_0600.h5'
  gchem(2)%fnames(2) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_23_0600.h5'

  gchem(1)%fnames(3) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_23_1200.h5'
  gchem(2)%fnames(3) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_23_1200.h5'

  gchem(1)%fnames(4) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_23_1800.h5'
  gchem(2)%fnames(4) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_23_1800.h5'

  gchem(1)%fnames(5) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_24_0000.h5'
  gchem(2)%fnames(5) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_24_0000.h5'

  gchem(1)%fnames(6) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_24_0600.h5'
  gchem(2)%fnames(6) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_24_0600.h5'

  gchem(1)%fnames(7) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_24_1200.h5'
  gchem(2)%fnames(7) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_24_1200.h5'

  gchem(1)%fnames(8) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_24_1800.h5'
  gchem(2)%fnames(8) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_24_1800.h5'

  gchem(1)%fnames(9) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_25_0000.h5'
  gchem(2)%fnames(9) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_25_0000.h5'

  gchem(1)%fnames(10) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_25_0600.h5'
  gchem(2)%fnames(10) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_25_0600.h5'

  gchem(1)%fnames(11) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_25_1200.h5'
  gchem(2)%fnames(11) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_25_1200.h5'

  gchem(1)%fnames(12) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_25_1800.h5'
  gchem(2)%fnames(12) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_25_1800.h5'

  gchem(1)%fnames(13) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_26_0000.h5'
  gchem(2)%fnames(13) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_26_0000.h5'

  gchem(1)%fnames(14) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_26_0600.h5'
  gchem(2)%fnames(14) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_26_0600.h5'

  gchem(1)%fnames(15) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_26_1200.h5'
  gchem(2)%fnames(15) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_26_1200.h5'

  gchem(1)%fnames(16) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_26_1800.h5'
  gchem(2)%fnames(16) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_26_1800.h5'

  gchem(1)%fnames(17) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_27_0000.h5'
  gchem(2)%fnames(17) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_27_0000.h5'

  gchem(1)%fnames(18) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_27_0600.h5'
  gchem(2)%fnames(18) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_27_0600.h5'

  gchem(1)%fnames(19) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_27_1200.h5'
  gchem(2)%fnames(19) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_27_1200.h5'

  gchem(1)%fnames(20) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_27_1800.h5'
  gchem(2)%fnames(20) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_27_1800.h5'

  gchem(1)%fnames(21) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_28_0000.h5'
  gchem(2)%fnames(21) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_28_0000.h5'

  gchem(1)%fnames(22) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_28_0600.h5'
  gchem(2)%fnames(22) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_28_0600.h5'

  gchem(1)%fnames(23) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_28_1200.h5'
  gchem(2)%fnames(23) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_28_1200.h5'

  gchem(1)%fnames(24) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_28_1800.h5'
  gchem(2)%fnames(24) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_28_1800.h5'

  gchem(1)%fnames(25) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_29_0000.h5'
  gchem(2)%fnames(25) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_29_0000.h5'

  gchem(1)%fnames(26) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_29_0600.h5'
  gchem(2)%fnames(26) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_29_0600.h5'

  gchem(1)%fnames(27) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_29_1200.h5'
  gchem(2)%fnames(27) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_29_1200.h5'

  gchem(1)%fnames(28) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_29_1800.h5'
  gchem(2)%fnames(28) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_29_1800.h5'

  gchem(1)%fnames(29) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_30_0000.h5'
  gchem(2)%fnames(29) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_30_0000.h5'

  gchem(1)%fnames(30) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_30_0600.h5'
  gchem(2)%fnames(30) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_30_0600.h5'

  gchem(1)%fnames(31) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_30_1200.h5'
  gchem(2)%fnames(31) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_30_1200.h5'

  gchem(1)%fnames(32) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_30_1800.h5'
  gchem(2)%fnames(32) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_30_1800.h5'

  gchem(1)%fnames(33) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_31_0000.h5'
  gchem(2)%fnames(33) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_31_0000.h5'

  gchem(1)%fnames(34) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_31_0600.h5'
  gchem(2)%fnames(34) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_31_0600.h5'

  gchem(1)%fnames(35) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_31_1200.h5'
  gchem(2)%fnames(35) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_31_1200.h5'

  gchem(1)%fnames(36) = '../olam_geoschem_files/anthOn_GL_Harvey_2017_08_31_1800.h5'
  gchem(2)%fnames(36) = '../olam_geoschem_files/anthOn_NA_Harvey_2017_08_31_1800.h5'

  ENDIF

  ! Fill s1900_gchem array with time for each geoschem file

  call date_abs_secs2(2017,08,23,000000,s1900_gchem(1))
  call date_abs_secs2(2017,08,23,060000,s1900_gchem(2))
  call date_abs_secs2(2017,08,23,120000,s1900_gchem(3))
  call date_abs_secs2(2017,08,23,180000,s1900_gchem(4))

  call date_abs_secs2(2017,08,24,000000,s1900_gchem(5))
  call date_abs_secs2(2017,08,24,060000,s1900_gchem(6))
  call date_abs_secs2(2017,08,24,120000,s1900_gchem(7))
  call date_abs_secs2(2017,08,24,180000,s1900_gchem(8))

  call date_abs_secs2(2017,08,25,000000,s1900_gchem(9))
  call date_abs_secs2(2017,08,25,060000,s1900_gchem(10))
  call date_abs_secs2(2017,08,25,120000,s1900_gchem(11))
  call date_abs_secs2(2017,08,25,180000,s1900_gchem(12))

  call date_abs_secs2(2017,08,26,000000,s1900_gchem(13))
  call date_abs_secs2(2017,08,26,060000,s1900_gchem(14))
  call date_abs_secs2(2017,08,26,120000,s1900_gchem(15))
  call date_abs_secs2(2017,08,26,180000,s1900_gchem(16))

  call date_abs_secs2(2017,08,27,000000,s1900_gchem(17))
  call date_abs_secs2(2017,08,27,060000,s1900_gchem(18))
  call date_abs_secs2(2017,08,27,120000,s1900_gchem(19))
  call date_abs_secs2(2017,08,27,180000,s1900_gchem(20))

  call date_abs_secs2(2017,08,28,000000,s1900_gchem(21))
  call date_abs_secs2(2017,08,28,060000,s1900_gchem(22))
  call date_abs_secs2(2017,08,28,120000,s1900_gchem(23))
  call date_abs_secs2(2017,08,28,180000,s1900_gchem(24))

  call date_abs_secs2(2017,08,29,000000,s1900_gchem(25))
  call date_abs_secs2(2017,08,29,060000,s1900_gchem(26))
  call date_abs_secs2(2017,08,29,120000,s1900_gchem(27))
  call date_abs_secs2(2017,08,29,180000,s1900_gchem(28))

  call date_abs_secs2(2017,08,30,000000,s1900_gchem(29))
  call date_abs_secs2(2017,08,30,060000,s1900_gchem(30))
  call date_abs_secs2(2017,08,30,120000,s1900_gchem(31))
  call date_abs_secs2(2017,08,30,180000,s1900_gchem(32))

  call date_abs_secs2(2017,08,31,000000,s1900_gchem(33))
  call date_abs_secs2(2017,08,31,060000,s1900_gchem(34))
  call date_abs_secs2(2017,08,31,120000,s1900_gchem(35))
  call date_abs_secs2(2017,08,31,180000,s1900_gchem(36))

  end subroutine geoschem_init

  !============================================================================

  subroutine geoschem_read0()

  use misc_coms,   only: io6, s1900_sim
  use mem_grid,    only: mza, mwa
  use ccnbin_coms, only: nccntyp 

  implicit none

  integer :: nf

  ! Loop over number of gchem_DATABASE file times and search for the one that
  ! corresponds to current or most recent model time.

  igchemfile = 0
  do nf = 1, ngchemfiles
     write(io6,*) 'ngchemf0 ',nf,gchem(1)%fnames(nf),' ',s1900_gchem(nf),' ',s1900_sim

     if (s1900_gchem(nf) <= s1900_sim) then
        igchemfile = nf
     endif
  enddo

  write(io6,*) 'geoschem_read0: Starting previous igchemfile = ',igchemfile

  if (igchemfile < 1) then
     write(io6,*) ' '
     write(io6,*) 'Unable to find previous or current gchem file for current'
     write(io6,*) 'model time.  Stopping model.'
     stop 'stop: no current gchem file'
  endif

  ! Allocate nudging arrays (containing geoschem fields interpolated to OLAM grid)

  allocate(geoschem_nud(mza,mwa,nccntyp)); geoschem_nud = 0.

  ! Call geoschem_read to read in nearest-in-time future geoschem aerosols
  ! and interpolate them to OLAM grid

  write(io6,'(a,i5)') 'geoschem_read0 calling geoschem_read ',igchemfile

  call geoschem_read()

  end subroutine geoschem_read0

  !============================================================================

  subroutine geoschem_read()

  use max_dims,   only: pathlen
  use misc_coms,  only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                        time8, runtype, s1900_init, s1900_sim
  use mem_grid,   only: mza, glatw, glonw, lpw
  use mem_basic,  only: press, rho
  use mem_ijtabs, only: jtab_w, jtw_prog
  use hdf5_utils, only: shdf5_irec, shdf5_orec, shdf5_open, shdf5_close

  implicit none

  integer :: io1,io2,jo1,jo2
  integer :: iw, j, kb
  integer :: iaer
  integer :: nlat, nlon, nlev, ids
  integer :: ndims, idims(3)

  real :: xperdeg,yperdeg
  real :: glat,glon
  real :: rio,rjo,wio1,wio2,wjo1,wjo2

  ! Column vectors of geos-chem fields at geos-chem vertical levels horizontally
  ! interpolated to IW location of OLAM grid

  real, allocatable :: presscol(:)   ! geos-chem pressure (hPa)
  real, allocatable :: aeroscol(:,:) ! geos-chem aerosol concentration (micrograms per m^3 at STP)

  ! Column vectors of OLAM fields at OLAM vertical levels at IW location

  real :: press4(mza) ! pressure [hPa]

  real :: vctr(mza,13) ! interpolated geos-chem aerosol concentration

  logical :: exans

  character(pathlen) :: flnm

  ! Processing next gchem file (only called with iaction = 1 if iupdgchem = 1)

  write(io6,'(a,i5)') 'geoschem_read1 ',igchemfile

  igchemfile = igchemfile + 1
   
  if (igchemfile > ngchemfiles) then
     write(io6,*) ' '
     write(io6,*) 'No future gchem file is available for nudging '
     write(io6,*) 'Stopping model '
     stop 'stop: no future gchem file'
  endif

  ! Loop over both (GL and NA) geoschem dataset regions.

  do ids = 1,2
     nlat = gchem(ids)%nlat
     nlon = gchem(ids)%nlon
     nlev = gchem(ids)%nlev

     write(io6,'(a,5i5)') 'gr11 ',ids,nlat,nlon,nlev

     ! Allocate gchem arrays for current dataset region

     allocate(gchem(ids)%dust1(nlon,nlat,nlev)) ! dust1
     allocate(gchem(ids)%dust2(nlon,nlat,nlev)) ! dust2
     allocate(gchem(ids)%dust3(nlon,nlat,nlev)) ! dust3
     allocate(gchem(ids)%dust4(nlon,nlat,nlev)) ! dust4
     allocate(gchem(ids)%asalt(nlon,nlat,nlev)) ! accumulation mode salt
     allocate(gchem(ids)%csalt(nlon,nlat,nlev)) ! coarse mode salt
     allocate(gchem(ids)%so4  (nlon,nlat,nlev)) ! so4
     allocate(gchem(ids)%no3  (nlon,nlat,nlev)) ! no3
     allocate(gchem(ids)%nh4  (nlon,nlat,nlev)) ! nh4
     allocate(gchem(ids)%bcpi (nlon,nlat,nlev)) ! hydrophillic black carbon
     allocate(gchem(ids)%bcpo (nlon,nlat,nlev)) ! hydrophobic black carbon
     allocate(gchem(ids)%ocpi (nlon,nlat,nlev)) ! hydrophillic organic carbon
     allocate(gchem(ids)%ocpo (nlon,nlat,nlev)) ! hydrophobic organic carbon
     allocate(gchem(ids)%presg(nlon,nlat,nlev)) ! pressure
     allocate(gchem(ids)%tempg(nlon,nlat,nlev)) ! temperature

     allocate(gchem(ids)%lat(nlat)) ! latitude
     allocate(gchem(ids)%lon(nlon)) ! longitude

     if (ids == 1) then
        allocate(presscol(nlev))
        allocate(aeroscol(nlev,13))
     endif

     ! Check if gchemfile exists

     flnm = trim(gchem(ids)%fnames(igchemfile))

     inquire(file=flnm, exist=exans)

     if (.not. exans) then

        ! gchemfile does not exist.

        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(io6,*) '!!!  GCHEMfile does not exist:'
        write(io6,*) '!!!  '//flnm
        write(io6,*) '!!!  Stopping run'
        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

        stop 'stop - no gchemfile'

     endif

     ! gchemfile exists; open it

     write(io6,'(a,5i5)') 'geoschem_read2 ',igchemfile,ids,flnm

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening gchemfile ', trim(flnm)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(flnm,'R')

     ! Open and read gchem_database file for current dataset region

     ndims    = 3
     idims(1) = nlon
     idims(2) = nlat
     idims(3) = nlev

  write(io6,'(a,5i5)') 'gr12 ',idims(1:3)

     call shdf5_irec(ndims,idims,'DUST1' ,rvar3=gchem(ids)%dust1)
     call shdf5_irec(ndims,idims,'DUST2' ,rvar3=gchem(ids)%dust2)
     call shdf5_irec(ndims,idims,'DUST3' ,rvar3=gchem(ids)%dust3)
     call shdf5_irec(ndims,idims,'DUST4' ,rvar3=gchem(ids)%dust4)
     call shdf5_irec(ndims,idims,'ASALT' ,rvar3=gchem(ids)%asalt)
     call shdf5_irec(ndims,idims,'CSALT' ,rvar3=gchem(ids)%csalt)
     call shdf5_irec(ndims,idims,'SO4'   ,rvar3=gchem(ids)%so4)
     call shdf5_irec(ndims,idims,'NO3'   ,rvar3=gchem(ids)%no3)
     call shdf5_irec(ndims,idims,'NH4'   ,rvar3=gchem(ids)%nh4)
     call shdf5_irec(ndims,idims,'BCPI'  ,rvar3=gchem(ids)%bcpi)
     call shdf5_irec(ndims,idims,'BCPO'  ,rvar3=gchem(ids)%bcpo)
     call shdf5_irec(ndims,idims,'OCPI'  ,rvar3=gchem(ids)%ocpi)
     call shdf5_irec(ndims,idims,'OCPO'  ,rvar3=gchem(ids)%ocpo)
     call shdf5_irec(ndims,idims,'PRESG' ,rvar3=gchem(ids)%presg)
     call shdf5_irec(ndims,idims,'TEMPG' ,rvar3=gchem(ids)%tempg)

     ndims    = 1
     idims(1) = nlat

     call shdf5_irec(ndims,idims,'LAT'  ,rvar1=gchem(ids)%lat)

     idims(1) = nlon

     call shdf5_irec(ndims,idims,'LON'  ,rvar1=gchem(ids)%lon)

     call shdf5_close()

  enddo ! ids

  ! Loop over OLAM grid columns for horizontal and vertical interpolation

  !$omp parallel private(press4,presscol,aeroscol,vctr)
  !$omp do private(iw,kb,glat,glon,ids,xperdeg,yperdeg,rio,rjo,io1,jo1, &
  !$omp            wio2,wjo2,wio1,wjo1,nlon,io2,jo2,iaer)
  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     kb = lpw(iw)
   
     glat = glatw(iw)
     glon = glonw(iw)

     glon = max(-179.999,min(179.999,glon))

     ! If iw point is inside NA dataset, use ids = 2; otherwise, use ids = 1 for GL dataset

     if (glat >   9.75 .and. &
         glat <   60.0 .and. &
         glon > -130.0 .and. &
         glon <  -60.0) then

        ids = 2

        xperdeg = 3.2
        yperdeg = 4.0

        rio = 1. + (glon + 130.) * xperdeg
        rjo = 1. + (glat - 9.75) * yperdeg
     else
        ids = 1

        xperdeg = 0.4
        yperdeg = 0.5

        rio = 1. + (glon + 180.) * xperdeg
        rjo = 1. + (glat +  90.) * yperdeg
     endif

     io1 = int(rio)
     jo1 = int(rjo)
         
     wio2 = rio - real(io1)
     wjo2 = rjo - real(jo1)
           
     wio1 = 1. - wio2
     wjo1 = 1. - wjo2

     nlon = gchem(ids)%nlon

     io2 = io1 + 1
     if (io2 > nlon) io2 = io2 - nlon
     jo2 = min(nlat, jo1 + 1)

     ! Fill column vector with OLAM model pressures at horizontal location iw,
     ! converting from Pa to hPa

     press4(:) = real(press(:,iw)) * .01

     ! Horizontally interpolate geos-chem pressure (in hPa) to location of
     ! OLAM grid column iw to fill a column vector at geos-chem vertical levels

     presscol(:) = wio1 * (wjo1 * gchem(ids)%presg(io1,jo1,:) + wjo2 * gchem(ids)%presg(io1,jo2,:)) &
                 + wio2 * (wjo1 * gchem(ids)%presg(io2,jo1,:) + wjo2 * gchem(ids)%presg(io2,jo2,:))

     ! Horizontally interpolate each geos-chem aerosol field to location of
     ! OLAM grid column iw to fill a column vector of the aeroscol array at
     ! geos-chem vertical levels

     aeroscol(:, 1) = wio1 * (wjo1 * gchem(ids)%dust1(io1,jo1,:) + wjo2 * gchem(ids)%dust1(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%dust1(io2,jo1,:) + wjo2 * gchem(ids)%dust1(io2,jo2,:))

     aeroscol(:, 2) = wio1 * (wjo1 * gchem(ids)%dust2(io1,jo1,:) + wjo2 * gchem(ids)%dust2(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%dust2(io2,jo1,:) + wjo2 * gchem(ids)%dust2(io2,jo2,:))

     aeroscol(:, 3) = wio1 * (wjo1 * gchem(ids)%dust3(io1,jo1,:) + wjo2 * gchem(ids)%dust3(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%dust3(io2,jo1,:) + wjo2 * gchem(ids)%dust3(io2,jo2,:))

     aeroscol(:, 4) = wio1 * (wjo1 * gchem(ids)%dust4(io1,jo1,:) + wjo2 * gchem(ids)%dust4(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%dust4(io2,jo1,:) + wjo2 * gchem(ids)%dust4(io2,jo2,:))

     aeroscol(:, 5) = wio1 * (wjo1 * gchem(ids)%asalt(io1,jo1,:) + wjo2 * gchem(ids)%asalt(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%asalt(io2,jo1,:) + wjo2 * gchem(ids)%asalt(io2,jo2,:))

     aeroscol(:, 6) = wio1 * (wjo1 * gchem(ids)%csalt(io1,jo1,:) + wjo2 * gchem(ids)%csalt(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%csalt(io2,jo1,:) + wjo2 * gchem(ids)%csalt(io2,jo2,:))

     aeroscol(:, 7) = wio1 * (wjo1 * gchem(ids)%so4(io1,jo1,:) + wjo2 * gchem(ids)%so4(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%so4(io2,jo1,:) + wjo2 * gchem(ids)%so4(io2,jo2,:))

     aeroscol(:, 8) = wio1 * (wjo1 * gchem(ids)%no3(io1,jo1,:) + wjo2 * gchem(ids)%no3(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%no3(io2,jo1,:) + wjo2 * gchem(ids)%no3(io2,jo2,:))

     aeroscol(:, 9) = wio1 * (wjo1 * gchem(ids)%nh4(io1,jo1,:) + wjo2 * gchem(ids)%nh4(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%nh4(io2,jo1,:) + wjo2 * gchem(ids)%nh4(io2,jo2,:))

     aeroscol(:,10) = wio1 * (wjo1 * gchem(ids)%bcpi(io1,jo1,:) + wjo2 * gchem(ids)%bcpi(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%bcpi(io2,jo1,:) + wjo2 * gchem(ids)%bcpi(io2,jo2,:))

     aeroscol(:,11) = wio1 * (wjo1 * gchem(ids)%bcpo(io1,jo1,:) + wjo2 * gchem(ids)%bcpo(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%bcpo(io2,jo1,:) + wjo2 * gchem(ids)%bcpo(io2,jo2,:))

     aeroscol(:,12) = wio1 * (wjo1 * gchem(ids)%ocpi(io1,jo1,:) + wjo2 * gchem(ids)%ocpi(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%ocpi(io2,jo1,:) + wjo2 * gchem(ids)%ocpi(io2,jo2,:))

     aeroscol(:,13) = wio1 * (wjo1 * gchem(ids)%ocpo(io1,jo1,:) + wjo2 * gchem(ids)%ocpo(io1,jo2,:)) &
                    + wio2 * (wjo1 * gchem(ids)%ocpo(io2,jo1,:) + wjo2 * gchem(ids)%ocpo(io2,jo2,:))

     ! Convert geos-chem aerosol concentrations in horizontally interpolated
     ! column at geos-chem vertical levels to number per kilogram of air

     ! Jack Kodros reported that aerosol abundances are given in micrograms per cubic meter at STP
     ! (email to Bill Cotton 23 May 2018).  Thus, to convert to micrograms per kg of air, we divide
     ! by the STP air density value, which is 1.276 kg/m^3.  To then convert to number (of aerosol
     ! particles) per kg of air, we further divide by the mean particle mass (in micrograms) for
     ! each aerosol category.  The mean particle mass is the "ccn_mass_ug" quantity printed out
     ! from the OLAM run during initialization.

     aeroscol(:, 1) = aeroscol(:, 1) / (1.276 * 5.9e-7 ) ! ug/m^3 @STP to #/kg_air dust1
     aeroscol(:, 2) = aeroscol(:, 2) / (1.276 * 3.4e-5 ) ! ug/m^3 @STP to #/kg_air dust2
     aeroscol(:, 3) = aeroscol(:, 3) / (1.276 * 1.9e-4 ) ! ug/m^3 @STP to #/kg_air dust3
     aeroscol(:, 4) = aeroscol(:, 4) / (1.276 * 9.8e-4 ) ! ug/m^3 @STP to #/kg_air dust4
     aeroscol(:, 5) = aeroscol(:, 5) / (1.276 * 6.4e-8 ) ! ug/m^3 @STP to #/kg_air asalt
     aeroscol(:, 6) = aeroscol(:, 6) / (1.276 * 1.3e-5 ) ! ug/m^3 @STP to #/kg_air csalt
     aeroscol(:, 7) = aeroscol(:, 7) / (1.276 * 2.9e-8 ) ! ug/m^3 @STP to #/kg_air so4
     aeroscol(:, 8) = aeroscol(:, 8) / (1.276 * 2.9e-8 ) ! ug/m^3 @STP to #/kg_air no3
     aeroscol(:, 9) = aeroscol(:, 9) / (1.276 * 2.9e-8 ) ! ug/m^3 @STP to #/kg_air nh4
     aeroscol(:,10) = aeroscol(:,10) / (1.276 * 7.0e-11) ! ug/m^3 @STP to #/kg_air bcpi
     aeroscol(:,11) = aeroscol(:,11) / (1.276 * 7.0e-11) ! ug/m^3 @STP to #/kg_air bcpo
     aeroscol(:,12) = aeroscol(:,12) / (1.276 * 3.2e-9 ) ! ug/m^3 @STP to #/kg_air ocpi
     aeroscol(:,13) = aeroscol(:,13) / (1.276 * 3.2e-9 ) ! ug/m^3 @STP to #/kg_air ocpo

     ! Vertically interpolate aerosol concentrations from geos-chem levels to
     ! model levels by comparing pressures.  Results are put into columns of
     ! the vctr array.

     do iaer = 1,13
        call pintrp_aa(nlev, aeroscol(:,iaer), presscol, mza, kb, vctr(:,iaer), press4)
     enddo

     ! Copy data from columns of vctr array to respective locations in
     ! geoschem_nud array, whose third index corresponds to index in ccntyp
     ! array.  Some geos-chem aerosols get combined into a single OLAM ccn 
     ! nudging array.

     geoschem_nud(kb:mza,iw,1) = vctr(kb:mza,1)                                     ! Dust1
     geoschem_nud(kb:mza,iw,2) = vctr(kb:mza,2)                                     ! Dust2
     geoschem_nud(kb:mza,iw,3) = vctr(kb:mza,3)                                     ! Dust3
     geoschem_nud(kb:mza,iw,4) = vctr(kb:mza,4)                                     ! Dust4
     geoschem_nud(kb:mza,iw,5) = vctr(kb:mza,5)  + vctr(kb:mza,6)                   ! Sea Salt
     geoschem_nud(kb:mza,iw,6) = vctr(kb:mza,7)  + vctr(kb:mza,8)  + vctr(kb:mza,9) ! SO4_NO3_NH4
     geoschem_nud(kb:mza,iw,7) = vctr(kb:mza,10) + vctr(kb:mza,11)                  ! Black Carbon
     geoschem_nud(kb:mza,iw,8) = vctr(kb:mza,12)                                    ! OCPI
     geoschem_nud(kb:mza,iw,9) = vctr(kb:mza,13)                                    ! OCPO

  enddo
  !$omp end do
  !$omp end parallel

  do ids = 1,2
     deallocate (gchem(ids)%dust1)
     deallocate (gchem(ids)%dust2)
     deallocate (gchem(ids)%dust3)
     deallocate (gchem(ids)%dust4)
     deallocate (gchem(ids)%asalt)
     deallocate (gchem(ids)%csalt)
     deallocate (gchem(ids)%so4)
     deallocate (gchem(ids)%no3)
     deallocate (gchem(ids)%nh4)
     deallocate (gchem(ids)%bcpi)
     deallocate (gchem(ids)%bcpo)
     deallocate (gchem(ids)%ocpi)
     deallocate (gchem(ids)%ocpo)
     deallocate (gchem(ids)%presg)
     deallocate (gchem(ids)%tempg)
     deallocate (gchem(ids)%lat)
     deallocate (gchem(ids)%lon)
  enddo

  end subroutine geoschem_read

!===============================================================================

  subroutine pintrp_aa(na,vctra,pressa,nb,kb0,vctrb,pressb)

  implicit none

  integer, intent(in)  :: na
  real,    intent(in)  :: vctra (na)
  real,    intent(in)  :: pressa(na)

  integer, intent(in)  :: nb
  integer, intent(in)  :: kb0
  real,    intent(out) :: vctrb (nb)
  real,    intent(in)  :: pressb(nb)

  integer :: ka
  integer :: kb
  real    :: grada

  ka = 1
  do kb = kb0,nb
     do while(ka < na-1 .and. pressb(kb) < pressa(ka+1))
        ka = ka + 1
     enddo
     
     if (pressb(kb) >= pressa(1)) then
        vctrb(kb) = vctra(1)
     elseif (pressb(kb) <= pressa(na)) then
        vctrb(kb) = vctra(na)
     else
        grada = (vctra(ka+1) - vctra(ka)) / (pressa(ka+1) - pressa(ka))
        vctrb(kb) = vctra(ka) + grada * (pressb(kb) - pressa(ka))
     endif
  enddo

  end subroutine pintrp_aa

  !===========================================================================

  subroutine nudge_geoschem()

  use mem_grid,    only: mza, lpw
  use mem_ijtabs,  only: istp, jtab_w, mrl_begl, jtw_prog, itab_w
  use mem_nudge,   only: tnudcent
  use ccnbin_coms, only: nccntyp
  use mem_micro,   only: ccntyp

  implicit none

  integer :: k, iw, mrl, ic, j
  real    :: tnudi

  ! Check whether it is time to nudge

  mrl = mrl_begl(istp)
  if (mrl < 1) return

  tnudi = 1. / 86400. ! Use 24-hour nudging timescale - as proposed

  ! Compute ozone nudging tendency using point-by-point (non-spectral) information

  !$omp parallel do private(iw,ic,k)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do ic = 1,nccntyp

        do k = lpw(iw), mza
           ccntyp(ic)%con_ccnt(k,iw) = ccntyp(ic)%con_ccnt(k,iw) &
              + (geoschem_nud(k,iw,ic) - ccntyp(ic)%con_ccn(k,iw)) * tnudi
        enddo

     enddo

  enddo
  !$omp end parallel do

  end subroutine nudge_geoschem

  !===========================================================================

  subroutine nudge_geoschem_startup()

  use mem_grid,    only: mza, lpw
  use mem_ijtabs,  only: istp, jtab_w, jtw_prog, itab_w
  use ccnbin_coms, only: nccntyp
  use mem_micro,   only: ccntyp

  implicit none

  integer :: k, iw, ic, j

  ! For start of new simulation, set CCN fields equal to geoschem nudging fields

  !$omp parallel do private(iw,ic,k)
  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     do ic = 1,nccntyp

        do k = lpw(iw), mza
           ccntyp(ic)%con_ccn(k,iw) = geoschem_nud(k,iw,ic)
        enddo

     enddo

  enddo
  !$omp end parallel do

  end subroutine nudge_geoschem_startup

  !============================================================================

  subroutine geoschem_translate()

  use misc_coms,  only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                        time8, runtype, s1900_init, s1900_sim
  use mem_grid,   only: mza, glatw, glonw, lpw
  use mem_basic,  only: press, rho
  use mem_ijtabs, only: jtab_w, jtw_prog
  use hdf5_utils, only: shdf5_irec, shdf5_orec, shdf5_open, shdf5_close

  implicit none

  integer :: iclobber = 1
  integer :: ndims, idims(4)
  integer :: ntim, nlat, nlon, nlev, itim, ids

  character(80) :: flnm

  ntim = 124

  ! Loop over both (GL and NA) geoschem dataset regions.

  do ids = 1,2
     nlat = gchem(ids)%nlat
     nlon = gchem(ids)%nlon
     nlev = gchem(ids)%nlev

     ! Allocate gchem arrays for current dataset region

     allocate(gchem(ids)%dust1(nlon,nlat,nlev)) ! dust1
     allocate(gchem(ids)%dust2(nlon,nlat,nlev)) ! dust2
     allocate(gchem(ids)%dust3(nlon,nlat,nlev)) ! dust3
     allocate(gchem(ids)%dust4(nlon,nlat,nlev)) ! dust4
     allocate(gchem(ids)%asalt(nlon,nlat,nlev)) ! accumulation mode salt
     allocate(gchem(ids)%csalt(nlon,nlat,nlev)) ! coarse mode salt
     allocate(gchem(ids)%so4  (nlon,nlat,nlev)) ! so4
     allocate(gchem(ids)%no3  (nlon,nlat,nlev)) ! no3
     allocate(gchem(ids)%nh4  (nlon,nlat,nlev)) ! nh4
     allocate(gchem(ids)%bcpi (nlon,nlat,nlev)) ! hydrophillic black carbon
     allocate(gchem(ids)%bcpo (nlon,nlat,nlev)) ! hydrophobic black carbon
     allocate(gchem(ids)%ocpi (nlon,nlat,nlev)) ! hydrophillic organic carbon
     allocate(gchem(ids)%ocpo (nlon,nlat,nlev)) ! hydrophobic organic carbon
     allocate(gchem(ids)%presg(nlon,nlat,nlev)) ! pressure
     allocate(gchem(ids)%tempg(nlon,nlat,nlev)) ! temperature

     allocate(gchem(ids)%lat(nlat)) ! latitude
     allocate(gchem(ids)%lon(nlon)) ! longitude

     allocate(gchem(ids)%readvar(nlon,nlat,nlev,ntim)) ! 4D array for reading Kodros files

     ! Loop over times to be processed

     do igchemfile = 1,ngchemfiles  ! "production" loop limits
   !  do igchemfile = 1,2           ! "testing" loop limits

        itim = igchemfile + 88  ! 88 skips the first 22 days (with 4 data times each day) in the kodros file

        write(6,'(a,3i10)') 'ids,igchemfile,itim ',ids,igchemfile,itim

        ! Specify Kodros file for current ids region (NA or GL)

        IF (anthOff) THEN  ! Use this for anthOff set of files

           if (ids == 1) then
              flnm = '../kodros_files/anthOff_2x25_Harvey2017.h5'
           else
              flnm = '../kodros_files/anthOff_NA_Harvey2017.h5'
           endif

        ELSE  ! Use this for anthOn set of files

           if (ids == 1) then
              flnm = '../kodros_files/anthOn_2x25_Harvey2017.h5'
           else
              flnm = '../kodros_files/anthOn_NA_Harvey2017.h5'
           endif

        ENDIF

        ! Open Kodros file

        write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6,*) 'Opening kodros file ', trim(flnm)
        write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

        call shdf5_open(flnm,'R')

        ! Read field values for current time from Kodros file

        ndims    = 4
        idims(1) = nlon
        idims(2) = nlat
        idims(3) = nlev
        idims(4) = ntim

        call shdf5_irec(ndims,idims,'DST1' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%dust1(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'DST2' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%dust2(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'DST3' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%dust3(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'DST4' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%dust4(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'SALA' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%asalt(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'SALC' ,rvar4=gchem(ids)%readvar)
        gchem(ids)%csalt(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'SO4'   ,rvar4=gchem(ids)%readvar)
        gchem(ids)%so4(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'NO3'   ,rvar4=gchem(ids)%readvar)
        gchem(ids)%no3(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'NH4'   ,rvar4=gchem(ids)%readvar)
        gchem(ids)%nh4(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'BCPI'  ,rvar4=gchem(ids)%readvar)
        gchem(ids)%bcpi(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'BCPO'  ,rvar4=gchem(ids)%readvar)
        gchem(ids)%bcpo(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'OCPI'  ,rvar4=gchem(ids)%readvar)
        gchem(ids)%ocpi(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'OCPO'  ,rvar4=gchem(ids)%readvar)
        gchem(ids)%ocpo(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'pressure',rvar4=gchem(ids)%readvar)
        gchem(ids)%presg(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        call shdf5_irec(ndims,idims,'temperature',rvar4=gchem(ids)%readvar)
        gchem(ids)%tempg(:,:,:) = gchem(ids)%readvar(:,:,:,itim)

        ndims    = 1
        idims(1) = nlat

        call shdf5_irec(ndims,idims,'lat'  ,rvar1=gchem(ids)%lat)

        idims(1) = nlon

        call shdf5_irec(ndims,idims,'lon'  ,rvar1=gchem(ids)%lon)

        ! Close Kodros file

        call shdf5_close()

        ! Open OLAM geoschem file

        flnm = trim(gchem(ids)%fnames(igchemfile))

        write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(io6,*) 'Opening gchemfile ', trim(flnm)
        write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

        call shdf5_open(flnm,'W',iclobber)

        ! Write field values for current time to OLAM geoschem file

        ndims    = 3
        idims(1) = nlon
        idims(2) = nlat
        idims(3) = nlev

        call shdf5_orec(ndims,idims,'DUST1'  ,rvar3=gchem(ids)%dust1)
        call shdf5_orec(ndims,idims,'DUST2'  ,rvar3=gchem(ids)%dust2)
        call shdf5_orec(ndims,idims,'DUST3'  ,rvar3=gchem(ids)%dust3)
        call shdf5_orec(ndims,idims,'DUST4'  ,rvar3=gchem(ids)%dust4)
        call shdf5_orec(ndims,idims,'ASALT'  ,rvar3=gchem(ids)%asalt)
        call shdf5_orec(ndims,idims,'CSALT'  ,rvar3=gchem(ids)%csalt)
        call shdf5_orec(ndims,idims,'SO4'    ,rvar3=gchem(ids)%so4)
        call shdf5_orec(ndims,idims,'NO3'    ,rvar3=gchem(ids)%no3)
        call shdf5_orec(ndims,idims,'NH4'    ,rvar3=gchem(ids)%nh4)
        call shdf5_orec(ndims,idims,'BCPI'   ,rvar3=gchem(ids)%bcpi)
        call shdf5_orec(ndims,idims,'BCPO'   ,rvar3=gchem(ids)%bcpo)
        call shdf5_orec(ndims,idims,'OCPI'   ,rvar3=gchem(ids)%ocpi)
        call shdf5_orec(ndims,idims,'OCPO'   ,rvar3=gchem(ids)%ocpo)
        call shdf5_orec(ndims,idims,'PRESG'  ,rvar3=gchem(ids)%presg)
        call shdf5_orec(ndims,idims,'TEMPG'  ,rvar3=gchem(ids)%tempg)

        ndims    = 1
        idims(1) = nlat

        call shdf5_orec(ndims,idims,'LAT'  ,rvar1=gchem(ids)%lat)

        idims(1) = nlon

        call shdf5_orec(ndims,idims,'LON'  ,rvar1=gchem(ids)%lon)

        ! Close OLAM goschem file

        call shdf5_close()

     enddo ! itim

     ! Deallocate arrays for current ids region 

     deallocate (gchem(ids)%dust1)
     deallocate (gchem(ids)%dust2)
     deallocate (gchem(ids)%dust3)
     deallocate (gchem(ids)%dust4)
     deallocate (gchem(ids)%asalt)
     deallocate (gchem(ids)%csalt)
     deallocate (gchem(ids)%so4)
     deallocate (gchem(ids)%no3)
     deallocate (gchem(ids)%nh4)
     deallocate (gchem(ids)%bcpi)
     deallocate (gchem(ids)%bcpo)
     deallocate (gchem(ids)%ocpi)
     deallocate (gchem(ids)%ocpo)
     deallocate (gchem(ids)%presg)
     deallocate (gchem(ids)%tempg)
     deallocate (gchem(ids)%lat)
     deallocate (gchem(ids)%lon)

  enddo ! ids

  end subroutine geoschem_translate

End module gchem_aeros

