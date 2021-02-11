MODULE test41alt

  !=======================================================================
  !
  !  Functions for setting up idealized initial conditions for the
  !  Staniforth, Ullrich, Melvin and Jablonowski baroclinic instability.
  !
  !  Options:
  !    deep    deep atmosphere (1 = yes or 0 = no)
  !
  !  Given a point specified by: 
  !     lon    longitude (radians) 
  !     lat    latitude (radians) 
  !     p/z    pressure (Pa) / height (m)
  ! zcoords    1 if z is specified, 0 if p is specified
  !
  !  the functions will return:
  !       p    pressure if z is specified and zcoords = 1 (Pa)
  !       u    zonal wind (m s^-1)
  !       v    meridional wind (m s^-1)
  !       w    vertical velocity (m s^-1)
  !       t    temperature (K)
  !    phis    surface geopotential (m^2 s^-2)
  !      ps    surface pressure (Pa)
  !     rho    density (kj m^-3)
  !
  !
  !  Authors: Paul Ullrich
  !           (University of Michigan, dcmip@ucar.edu)
  !
  !=======================================================================

  use dcmip_initial_conditions_test_4, only: epv, theta

  IMPLICIT NONE

!=======================================================================
! use physical constants
!=======================================================================

  REAL(8), PARAMETER ::               &
       a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
       Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616d0,             & ! Gravity (m s^2)
       cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
       pi    = 3.14159265358979d0,    & ! pi
       p0    = 100000.0d0,            & ! surface pressure (Pa)
       kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
       omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

!-----------------------------------------------------------------------
! parameters
!----------------------------------------------------------------------- 
  REAL(8), PARAMETER ::               &
       T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
       T0P        = 240.d0     ,      & ! temperature at polar surface (K)
       B          = 2.d0       ,      & ! jet half-width parameter
       K          = 3.d0       ,      & ! jet width parameter
       lapse      = 0.005d0             ! lapse rate parameter

  REAL(8), PARAMETER ::               &
       pertu0     = 0.5d0      ,      & ! Perturbation wind velocity (m/s)
       pertr      = 1.d0/6.d0  ,      & ! Perturbation radius (Earth radii)
       pertlon    = pi/9.d0    ,      & ! Perturbation longitude
       pertlat    = 2.d0*pi/9.d0,     & ! Perturbation latitude
       pertz      = 15000.d0   ,      & ! Perturbation height cap
       dxepsilon  = 1.d-5               ! Small value for numerical derivatives
       

CONTAINS

  SUBROUTINE baroclinic_instability_alt(deep,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q1,q2)
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    REAL(8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q1,         & ! Tracer q1 - Potential temperature
                q2            ! Tracer q2 - Ertel's potential vorticity

    REAL(8) :: aref, omegaref, eta
    REAL(8) :: T0, constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT, inttermU, bigU, rcoslat, omegarcoslat

!-----------------------------------------------------------------------
!    constants / parameters
!-----------------------------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    ! Won't work when input is pressure
    if (zcoords .eq. 0) then
      print *, "Pressure based vertical coordinate not supported for initialization"
      return
    end if

    scaledZ = z / (B * constH)

!-----------------------------------------------------------------------
!    tau values
!-----------------------------------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

!-----------------------------------------------------------------------
!    radius ratio
!-----------------------------------------------------------------------
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

!-----------------------------------------------------------------------
!    interior term on temperature expression
!-----------------------------------------------------------------------
    inttermT = (rratio * cos(lat))**K &
             - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

!-----------------------------------------------------------------------
!    temperature
!-----------------------------------------------------------------------
    t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))

!-----------------------------------------------------------------------
!    hydrostatic pressure
!-----------------------------------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

!-----------------------------------------------------------------------
!    density (via ideal gas law)
!-----------------------------------------------------------------------
    rho = p / (Rd * t)

!-----------------------------------------------------------------------
!    calculate pointwise surface pressure
!-----------------------------------------------------------------------
    ps = p0

!-----------------------------------------------------------------------
!    velocity field
!-----------------------------------------------------------------------
    inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    bigU = g / aref * K * inttau2 * inttermU * t
    if (deep .eq. 0) then
      rcoslat = aref * cos(lat)
    else
      rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0
    w = 0.d0

!-----------------------------------------------------------------------
!    perturbation on velocity field (streamfunction type)
!-----------------------------------------------------------------------
    u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
        ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
        - evaluate_streamfunction(lon, lat - dxepsilon, z))

    v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
        ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
        - evaluate_streamfunction(lon - dxepsilon, lat, z))

!-----------------------------------------------------------------------
!    surface geopotential
!-----------------------------------------------------------------------
    phis = 0.d0

!-----------------------------------------------------------------------
!    initialize density from ideal gas law
!-----------------------------------------------------------------------
    rho = p/(t*Rd)

!-----------------------------------------------------------------------
!     tracer q1, potential temperature
!-----------------------------------------------------------------------
    eta = p / 1.0d5
    q1 = theta(lon,lat,eta)

!-----------------------------------------------------------------------
!     tracer q2, absolute value of EPV
!-----------------------------------------------------------------------
    q2 = abs(epv(lon,lat,eta)) * X ! the value of |EPV| scales with X
!
!   The absolute value of Ertel's potential vorticity
!   is selected to avoid the negative EPV values in the
!   Southern Hemisphere. Such negative values might interfere with positive-definite 
!   constraints in the tracer advection algorithm.

  END SUBROUTINE baroclinic_instability_alt

!-----------------------------------------------------------------------
!    Stream function perturbation function
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_streamfunction(lon, lat, z)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(8) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.d0 / pertr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
      cospert = cos(0.5d0 * pi * greatcircler)
    else
      cospert = 0.d0
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

END MODULE test41alt
