MODULE test2Xdeep

  !=======================================================================
  !
  !  Functions for setting up idealized initial conditions for the
  !  mountain-generated waves over a Schar-type mountain on the sphere
  !  in a deep atmosphere.
  !
  !  Options:
  ! subtest    sub-test case (1 or 2)
  !       X    scale factor of the Earth
  !
  !  Given a point specified by: 
  !     lon    longitude (radians) 
  !     lat    latitude (radians) 
  !     p/z    pressure (Pa)/height (m)
  ! rcoords    1 if r is specified, 0 if p is specified
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

  IMPLICIT NONE

!=======================================================================
! use physical constants
!=======================================================================

  REAL(8), PARAMETER ::                               &
       a     = 6371220.0d0,                           & ! Reference Earth's Radius (m)
       Rd    = 287.0d0,                               & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616d0,                             & ! Gravity (m s^2)
       cp    = 1004.5d0,                              & ! Specific heat capacity (J kg^-1 K^1)
       pi    = 3.14159265358979d0,                    & ! pi
       p0    = 100000.0d0,                            & ! surface pressure (Pa)
       kappa = 2.d0/7.d0,                             & ! Ratio of Rd to cp
       omega = 7.29212d-5,                            & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180.d0                             ! Conversion factor of degrees to radians

!-----------------------------------------------------------------------
! parameters
!----------------------------------------------------------------------- 
  REAL(8), PARAMETER ::                            &
       ueq        = 20.d0     ,                    & ! horizontal wind at equator surface (m/s)
       utop1      = 20.d0     ,                    & ! horizontal wind at equator top of atmosphere (m/s, subtest 1)
       utop2      = 80.d0     ,                    & ! horizontal wind at equator top of atmosphere (m/s, subtest 2)
       ztop       = 30000.d0  ,                    & ! top of the atmosphere (m)
       peq        = 100000.d0 ,                    & ! pressure at equator surface
       Teq        = 300.d0                           ! temperature at equator surface

CONTAINS

  SUBROUTINE deep_mountain_wave(subtest,X,lon,lat,p,r,rcoords,u,v,w,t,phis,ps,rho)
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: subtest      ! Unsheared (1) or sheared (0) test case

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Scale factor

    REAL(8), INTENT(INOUT) :: &
                p,            & ! Pressure  (Pa)
                r               ! Radius (m)

    INTEGER, INTENT(IN) :: rcoords     ! 0 or 1 see below

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho           ! density (kg m^-3)

    REAL(8) :: aref, utop, K, coeffA, coeffB, coeffC

!-----------------------------------------------------------------------
!    subtest
!-----------------------------------------------------------------------
     if (subtest .eq. 1) then
       utop = utop1
     elseif (subtest .eq. 2) then
       utop = utop2
     else
       print *,"Invalid subtest parameter (expected 1 or 2)"
       return
     end if

!-----------------------------------------------------------------------
!    scaled Earth radius
!-----------------------------------------------------------------------
     aref = a / X

!-----------------------------------------------------------------------
!    K coefficient
!-----------------------------------------------------------------------
     K = ((aref + ztop) * ueq)**2 - (aref * utop)**2
     if (K .eq. 0.d0) then
       K = 0.d0
     else
       K = 1.d0 / (0.5d0 - aref * g * ztop * (aref + ztop) / K)
     end if

!-----------------------------------------------------------------------
!    calculate radius associated with pointwise pressure
!-----------------------------------------------------------------------
    if (rcoords .eq. 0) then
      coeffA = (ueq * cos(lat))**2 / (2.d0 * Rd * Teq * (K - 2.d0))
      coeffB = 2.d0 * g / (Rd * Teq * (K - 2)) + g * K * cos(lat)**2 / (Rd * Teq * (K - 2.d0))
      coeffC = - log(p / p0)

      r = (- coeffB + sqrt(coeffB**2 - 4.d0 * coeffA * coeffC)) / (2.d0 * coeffA)
    end if

!-----------------------------------------------------------------------
!    calculate pointwise pressure associated with radius
!-----------------------------------------------------------------------
    if (rcoords .eq. 1) then
      p = (0.5d0 * (ueq * r / aref)**2 + r * g * K / (K - 2.d0) * (1.d0 - r / aref))
      p = 1.d0 / (Rd * Teq) * (2.d0 * r * g / (K - 2.d0) + p * cos(lat)**2)
      p = p0 * exp(- 1.d0 / (Rd * Teq) * (2.d0 * aref * g / (K - 2.d0) + 0.5d0 * ueq**2)) * exp(p)
    end if

!-----------------------------------------------------------------------
!    calculate pointwise surface pressure
!-----------------------------------------------------------------------
    ps = p0 * exp(- 1.d0 / (Rd * Teq) * 0.5d0 * ueq**2 * (1.d0 - cos(lat)**2))

!-----------------------------------------------------------------------
!    temperature field
!-----------------------------------------------------------------------
    t = Teq * (K - 2.d0) / (K * cos(lat)**2 - 2.d0)

!-----------------------------------------------------------------------
!    velocity field
!-----------------------------------------------------------------------
    u = 1.d0 / (K * cos(lat)**2 - 2.d0) * (               &
            + 2.d0 * r * g * K                            &
            - 2.d0 * r**2 * g * K / aref                  &
            + r**2 * ueq**2 * K / aref**2                 &
            - 2.d0 * r**2 * ueq**2 / aref**2)

    u = sqrt(u) * cos(lat)

    v = 0.d0
    w = 0.d0

!-----------------------------------------------------------------------
!    surface geopotential
!-----------------------------------------------------------------------
    phis = 0.d0

!-----------------------------------------------------------------------
!    initialize density from ideal gas law
!-----------------------------------------------------------------------
    rho = p/(t*Rd)

  END SUBROUTINE deep_mountain_wave

END MODULE test2Xdeep
