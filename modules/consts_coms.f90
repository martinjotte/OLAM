Module consts_coms

  use iso_fortran_env, only: r4 => real32, & ! single precision real
                             r8 => real64, & ! double precision real
                             i1 => int8,   & ! one byte integer (char)
                             i2 => int16,  & ! two byte integer (short)
                             i4 => int32,  & ! four byte integer (int)
                             i8 => int64     ! eight byte integer (long)

  implicit none (external, type)

  ! OLAM double-precision constants

  real(r8), parameter ::     grav_r8 =      9.81_r8  ! acceleration of gravity [m/s^2]
  real(r8), parameter ::       cp_r8 =     1004._r8  ! dry air spec heat at cnst P [J/(kg K)]
  real(r8), parameter ::       cv_r8 =      717._r8  ! dry air spec heat at cnst vol [J/(kg K)]
  real(r8), parameter ::     rvap_r8 =      461._r8  ! water vapor gas constant [J/(kg K)]
  real(r8), parameter ::      p00_r8 =      1.e5_r8  ! reference pressure [Pa]
  real(r8), parameter ::      t00_r8 =    273.15_r8  ! 0 Celsius temperature [K]
  real(r8), parameter ::     alvl_r8 =    2.50e6_r8  ! latent heat of evaporation [J/kg]
  real(r8), parameter ::     alvi_r8 =   2.834e6_r8  ! latent heat of sublimation [J/kg]
  real(r8), parameter ::     alli_r8 =   0.334e6_r8  ! latent heat of fusion [J/kg]
  real(r8), parameter ::     cliq_r8 =     4186._r8  ! spec heat of liquid [J/(kg K)]
  real(r8), parameter ::     cice_r8 =     2106._r8  ! spec heat of ice [J/(kg K)]
  real(r8), parameter ::     vonk_r8 =      0.40_r8  ! von Karman constant
  real(r8), parameter ::       nu_r8 =  1.568e-5_r8  ! kinematic viscosity of air [m^2/s]
  real(r8), parameter ::    tkmin_r8 =     5.e-4_r8  ! minimum allowed TKE
  real(r8), parameter ::    solar_r8 =   1368.22_r8  ! solar constant [W/m^2]
  real(r8), parameter ::   stefan_r8 = 5.6696e-8_r8  ! Stefan-Boltzmann constant

  ! OLAM double-precision derived constants

  real(r8), parameter ::     rdry_r8 =    cp_r8 -    cv_r8
  real(r8), parameter ::     rocp_r8 =  rdry_r8 /    cp_r8
  real(r8), parameter ::    mrocp_r8 = -rdry_r8 /    cp_r8
  real(r8), parameter ::     cpor_r8 =    cp_r8 /  rdry_r8
  real(r8), parameter ::     rocv_r8 =  rdry_r8 /    cv_r8
  real(r8), parameter ::     cvor_r8 =    cv_r8 /  rdry_r8
  real(r8), parameter ::    cpocv_r8 =    cp_r8 /    cv_r8
  real(r8), parameter ::    cvocp_r8 =    cv_r8 /    cp_r8
  real(r8), parameter ::      cpi_r8 =   1.0_r8 /    cp_r8
  real(r8), parameter ::    cpio2_r8 =   0.5_r8 /    cp_r8
  real(r8), parameter ::     cpi4_r8 =   4.0_r8 *   cpi_r8
  real(r8), parameter ::     p00i_r8 =    1._r8 /   p00_r8
  real(r8), parameter ::     p00k_r8 =   p00_r8 ** rocp_r8
  real(r8), parameter ::    p00ki_r8 =   1.0_r8 /  p00k_r8
  real(r8), parameter ::  p00kord_r8 =  p00k_r8 /  rdry_r8
  real(r8), parameter ::      pc1_r8 =  p00i_r8 ** rocv_r8
  real(r8), parameter ::    grav2_r8 =   2.0_r8 *  grav_r8
  real(r8), parameter ::   gravo2_r8 =   0.5_r8 *  grav_r8
  real(r8), parameter ::    gravi_r8 =   1.0_r8 /  grav_r8
  real(r8), parameter ::     cpog_r8 =    cp_r8 /  grav_r8
  real(r8), parameter ::     gocp_r8 =  grav_r8 /    cp_r8
  real(r8), parameter ::   gordry_r8 =  grav_r8 /  rdry_r8
  real(r8), parameter ::   rdryog_r8 =  rdry_r8 /  grav_r8
  real(r8), parameter ::    alvl2_r8 =  alvl_r8 *  alvl_r8  ! alvl * alvl
  real(r8), parameter ::    alvi2_r8 =  alvi_r8 *  alvi_r8  ! alvi * alvi
  real(r8), parameter :: alli1000_r8 =  alli_r8 * 1000._r8  ! latent heat of fusion [J/(1000 kg)]
  real(r8), parameter :: cliq1000_r8 =  cliq_r8 * 1000._r8  ! spec heat of liquid [J/(1000 kg K)]
  real(r8), parameter :: cice1000_r8 =  cice_r8 * 1000._r8  ! spec heat of ice [J/(1000 kg K)]
  real(r8), parameter ::    cliqi_r8 =    1._r8 /  cliq_r8
  real(r8), parameter ::    cicei_r8 =    1._r8 /  cice_r8
  real(r8), parameter ::    allii_r8 =    1._r8 /  alli_r8
  real(r8), parameter ::   alvlor_r8 =  alvl_r8 /  rdry_r8
  real(r8), parameter ::  alvlocp_r8 =  alvl_r8 /    cp_r8
  real(r8), parameter ::  alviocp_r8 =  alvi_r8 /    cp_r8
  real(r8), parameter ::  eps_vap_r8 =  rdry_r8 /  rvap_r8
  real(r8), parameter :: eps_vapi_r8 =  rvap_r8 /  rdry_r8
  real(r8), parameter :: eps_virt_r8 = (rvap_r8 -  rdry_r8) / rdry_r8
  real(r8), parameter ::   cp253i_r8 =   cpi_r8 /  253._r8
  real(r8), parameter ::    vonki_r8 =    1._r8 /  vonk_r8

  ! double-precision pi

  real(r8), parameter ::   pio4_r8 = atan(1._r8)
  real(r8), parameter ::   pio2_r8 =   2._r8 * pio4_r8
  real(r8), parameter ::    pi1_r8 =   2._r8 * pio2_r8
  real(r8), parameter ::    pi2_r8 =   2._r8 *  pi1_r8
  real(r8), parameter ::    pi4_r8 =   2._r8 *  pi2_r8
  real(r8), parameter :: pio180_r8 =  pi1_r8 / 180._r8   ! radians per degree
  real(r8), parameter :: piu180_r8 = 180._r8 /  pi1_r8   ! degrees per radian

  ! OLAM single-precision constants

  real, parameter :: grav   = real(   grav_r8 )   ! acceleration of gravity [m/s^2]
  real, parameter :: cp     = real(     cp_r8 )   ! dry air spec heat at const P [J/(kg K)]
  real, parameter :: cv     = real(     cv_r8 )   ! dry air spec heat at const vol [J/(kg K)]
  real, parameter :: rvap   = real(   rvap_r8 )   ! water vapor gas constant [J/(kg K)]
  real, parameter :: p00    = real(    p00_r8 )   ! reference pressure [Pa]
  real, parameter :: t00    = real(    t00_r8 )   ! 0 Celsius temperature [K]
  real, parameter :: alvl   = real(   alvl_r8 )   ! latent heat of evaporation [J/kg]
  real, parameter :: alvi   = real(   alvi_r8 )   ! latent heat of sublimation [J/kg]
  real, parameter :: alli   = real(   alli_r8 )   ! latent heat of fusion [J/kg]
  real, parameter :: cliq   = real(   cliq_r8 )   ! spec heat of liquid [J/(kg K)]
  real, parameter :: cice   = real(   cice_r8 )   ! spec heat of ice [J/(kg K)]
  real, parameter :: vonk   = real(   vonk_r8 )   ! von Karman constant
  real, parameter :: nu     = real(     nu_r8 )   ! kinematic viscosity of air [m^2/s]
  real, parameter :: tkmin  = real(  tkmin_r8 )   ! minimum allowed TKE
  real, parameter :: solar  = real(  solar_r8 )   ! solar constant [W/m^2]
  real, parameter :: stefan = real( stefan_r8 )   ! Stefan-Boltzmann constant

  ! OLAM single-precision derived constants

  real, parameter :: rdry     =    cp - cv
  real, parameter :: rocp     =  rdry / cp
  real, parameter :: mrocp    = -rdry / cp
  real, parameter :: cpor     =    cp / rdry
  real, parameter :: rocv     =  rdry / cv
  real, parameter :: cvor     =    cv / rdry
  real, parameter :: cpocv    =    cp / cv
  real, parameter :: cvocp    =    cv / cp
  real, parameter :: cpi      =   1.0 / cp
  real, parameter :: cpio2    =   0.5 / cp
  real, parameter :: cpi4     =   4.0 * cpi
  real, parameter :: p00i     =    1. / p00
  real, parameter :: p00k     =   p00 ** rocp
  real, parameter :: p00ki    =   1.0 / p00k
  real, parameter :: p00kord  =  p00k / rdry
  real, parameter :: pc1      =  p00i ** rocv
  real, parameter :: grav2    =   2.0 * grav
  real, parameter :: gravo2   =   0.5 * grav
  real, parameter :: gravi    =   1.0 / grav
  real, parameter :: cpog     =    cp / grav
  real, parameter :: gocp     =  grav / cp
  real, parameter :: gordry   =  grav / rdry
  real, parameter :: rdryog   =  rdry / grav
  real, parameter :: alvl2    =  alvl * alvl   ! alvl * alvl
  real, parameter :: alvi2    =  alvi * alvi   ! alvi * alvi
  real, parameter :: alli1000 =  alli * 1000.  ! latent heat of fusion [J/(1000 kg)]
  real, parameter :: cliq1000 =  cliq * 1000.  ! spec heat of liquid [J/(1000 kg K)]
  real, parameter :: cice1000 =  cice * 1000.  ! spec heat of ice [J/(1000 kg K)]
  real, parameter :: cliqi    =    1. / cliq
  real, parameter :: cicei    =    1. / cice
  real, parameter :: allii    =    1. / alli
  real, parameter :: alvlor   =  alvl / rdry
  real, parameter :: alvlocp  =  alvl / cp
  real, parameter :: alviocp  =  alvi / cp
  real, parameter :: eps_vap  =  rdry / rvap
  real, parameter :: eps_vapi =  rvap / rdry
  real, parameter :: eps_virt = (rvap - rdry) / rdry
  real, parameter :: cp253i   =   cpi / 253.
  real, parameter :: vonki    =    1. / vonk

  ! single-precision pi

  real, parameter :: pio4   = atan(1.)
  real, parameter :: pio2   =   2. * pio4
  real, parameter :: pi1    =   2. * pio2
  real, parameter :: pi2    =   2. * pi1
  real, parameter :: pi4    =   2. * pi2
  real, parameter :: pio180 =  pi1 / 180.  ! radians per degree
  real, parameter :: piu180 = 180. / pi1   ! degrees per radian

  ! single-precision constants defined in init_consts

  real            :: dlat     ! Earth 1-latitude-degree dist [m]
  real            :: omega    ! Earth rotational angular velocity
  real            :: omega2   ! omega * 2
  real            :: xscale   ! erad scaling factor (for NCAR DCMIP test cases)
  real            :: erad     ! Earth radius [m]
  real            :: erad2
  real            :: erad4
  real            :: eradsq
  real            :: erador5
  real            :: eradi
  real            :: erad2sq
  real            :: fcoriol  ! omega2 for global run,
                              ! omega2 * sin(lat) for Cartesian run

Contains

!===============================================================================

  subroutine init_consts(mdomain,test_case,rlat0)

    import, only: dlat,  omega,   omega2,  xscale,  erad, erad2, erad4, eradsq, &
                  eradi, erador5, erad2sq, fcoriol, pio180

    implicit none (external, type)

    integer, intent(in) :: mdomain
    integer, intent(in) :: test_case
    real,    intent(in) :: rlat0

    ! Standard (Earth) values

    erad   = 6371.22e3          ! Earth radius [m]
    omega  = 7.29212e-5         ! Earth rotational angular velocity
    xscale = 1.0

    ! Alternate values depending on test case

    if     (test_case == 200) then
       omega = 0.
    elseif (test_case == 201) then
       xscale = 500.
       erad = erad / xscale
       omega = 0.
    elseif (test_case == 21 ) then
       xscale = 500.
       erad = erad / xscale
       omega = 0.
    elseif (test_case == 22 ) then
       xscale = 500.
       erad = erad / xscale
       omega = 0.
    elseif (test_case == 31 ) then
       xscale = 125.
       erad = erad / xscale
       omega = 0.
    elseif (test_case == 411) then
       xscale = 10.
       erad = erad / xscale
       omega = omega * xscale
    elseif (test_case == 412) then
       xscale = 100.
       erad = erad / xscale
       omega = omega * xscale
    elseif (test_case == 413) then
       xscale = 1000.
       erad = erad / xscale
       omega = omega * xscale
    elseif (test_case == 131) then
       xscale = 120.
       erad = erad / xscale
       omega = 0.
    endif

    ! Secondary values

    erad2    = erad * 2.
    erad4    = erad * 4.
    eradsq   = erad * erad
    erador5  = erad / sqrt(5.)
    eradi    = 1.0  / erad
    erad2sq  = erad2 * erad2
    dlat     = erad * pio180
    omega2   = omega * 2.

    if (mdomain <= 1) then

       ! global run
       fcoriol = omega2

    else

       ! limited area run; vxe is defined E-W, vye is defined N-S
       fcoriol = omega2 * sin(pio180 * rlat0)

    endif

  end subroutine init_consts

End Module consts_coms
