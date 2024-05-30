module umwm_sheltering
  ! Module that provides wave sheltering functions.
  use consts_coms, only: pi1

  implicit none

  private
  public :: sheltering_coare35, sheltering_reynolds

contains

!===============================================================================

  pure real function polyval(p, x)
    ! Evaluates a polynomial of order p over values of x:
    !   y = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)
    real, intent(in) :: p(:), x
    integer :: n
    polyval = 0
    do n = 1, size(p)
      polyval = polyval + p(n) * x**(size(p) - n)
    enddo
  end function polyval

!===============================================================================

  pure elemental real function reynolds_number(wspd, wind_height, swh, ust,&
    cp_peak, k_peak, air_viscosity) result(res)
    ! Computes the wave Reynolds number.
    real, intent(in) :: wspd, wind_height, swh, ust, cp_peak, k_peak,&
                            air_viscosity
    real, parameter :: re_max = 4e7, re_min = 6.62e2
    res = (wspd + 2.5 * ust * log(pi1 / (wind_height * k_peak)) - cp_peak)&
        * swh / air_viscosity
    res = max(min(res, re_max), re_min)
  end function reynolds_number

!===============================================================================

  pure elemental real function sheltering_reynolds(wspd, wind_height, swh,&
    ust, cp_peak, k_peak, air_viscosity) result(res)
    real, intent(in) :: wspd, wind_height, swh, ust, cp_peak, k_peak,&
                            air_viscosity
    real :: reynolds
    real, parameter :: pp(8) = [&
      8.08282075e-07, -7.22452970e-05, 2.69676878e-03, -5.45054760e-02,&
      6.43417148e-01, -4.42303614e+00, 1.63178789e+01, -2.47019533e+01]
    reynolds = reynolds_number(wspd, wind_height, swh, ust, cp_peak, k_peak, air_viscosity)
    res = polyval(pp, log(reynolds))
  end function sheltering_reynolds

!===============================================================================

  pure elemental real function sheltering_coare35(wspd) result(res)
    ! Returns the sheltering coefficient given input wind speed
    ! to approximately match the momentum flux of COARE 3.5 algorithm.
    real,    intent(in) :: wspd ! wind speed [m/s]
    real, parameter :: x1 = 15, x2 = 33, x3 = 60
    real, parameter :: y1 = 0.10, y2 = 0.09, y3 = 0.06
    real, parameter :: intercept = 0.04, curvature = 0.65, decay = 1.6

    real :: a, b, c, s1, slope1, slope2

    ! low range linear growth
    slope1 = (y1 - intercept) / x1
    slope2 = (y3 - y2) / (x3 - x2)

    ! medium range fitting (quadratic)
    c = curvature * (slope2 - slope1) / (x2 - x1)
    b = slope1 - 2 * c * x1
    a = y1 - b * x1 - c * x1**2

    ! high end decay
    s1 = a + b * x2 + c * x2**2

    if (wspd <= x1) then
      res = intercept + slope1 * wspd
    elseif (wspd > x1 .and. wspd <= x2) then
      res = a + b * wspd + c * wspd**2
    else
      res = s1 * exp(- (wspd - x2) / (decay * wspd))
    endif

  end function sheltering_coare35

end module umwm_sheltering
