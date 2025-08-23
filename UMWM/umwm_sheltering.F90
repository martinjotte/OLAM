module umwm_sheltering

  ! Module that provides wave sheltering functions.
  implicit none

  private
  public :: sheltering_coare35, sheltering_reynolds

contains

!===============================================================================

pure elemental real function reynolds_number(wspd, wind_height, swh, ust, &
                    cp_peak, k_peak) result(res)

  ! Computes the wave Reynolds number.

  use consts_coms, only: pi1, vonki
  use umwm_module, only: nu_air

  implicit none

  real, intent(in) :: wspd, wind_height, swh, ust, cp_peak, k_peak
  real, parameter  :: re_max = 4.00e7
  real, parameter  :: re_min = 6.62e2

  res = (wspd + vonki * ust * log(pi1 / (wind_height * k_peak)) - cp_peak) * swh / nu_air

  res = max(min(res, re_max), re_min)

end function reynolds_number

!===============================================================================

pure elemental real function sheltering_reynolds(wspdw, swh) result(res)

  use umwm_module, only: nu_air

  implicit none

  real, intent(in) :: wspdw, swh
  real             :: reynolds, x

  real,  parameter :: re_max = 4.00e7
  real,  parameter :: re_min = 6.62e2

  real,  parameter :: p(8) = [ &
       8.08282075e-07, -7.22452970e-05, 2.69676878e-03, -5.45054760e-02, &
       6.43417148e-01, -4.42303614e+00, 1.63178789e+01, -2.47019533e+01  ]

  reynolds = max( re_min, min(re_max, abs(wspdw) * swh / nu_air) )

  x = log(reynolds)

  res = p(8) + x * (p(7) + x * (p(6) + x * (p(5) + x * (p(4) + x * (p(3) + x * (p(2) + x * p(1)))))))

end function sheltering_reynolds

!===============================================================================

pure elemental real function sheltering_coare35(wspd) result(res)

  ! Returns the sheltering coefficient given input wind speed
  ! to approximately match the momentum flux of COARE 3.5 algorithm.

  implicit none

  real, intent(in) :: wspd ! wind speed [m/s]

  real, parameter  :: x1 = 15.0, x2 = 33.0, x3 = 60.0
  real, parameter  :: y1 = 0.10, y2 = 0.09, y3 = 0.06
  real, parameter  :: intercept = 0.04, curvature = 0.65, decay = 1.6

  ! low range linear growth
  real, parameter  :: slope1 = (y1 - intercept) / x1
  real, parameter  :: slope2 = (y3 - y2) / (x3 - x2)

  ! medium range fitting (quadratic)
  real, parameter  :: c = curvature * (slope2 - slope1) / (x2 - x1)
  real, parameter  :: b = slope1 - 2. * c * x1
  real, parameter  :: a = y1 - x1 * (b + c * x1)

  ! high end decay
  real, parameter  :: s1 = a + x2 * (b + c * x2)

  if (wspd <= x1) then
     res = intercept + slope1 * wspd
  elseif (wspd > x1 .and. wspd <= x2) then
     res = a + wspd * (b + c * wspd)
  else
     res = s1 * exp( -(wspd - x2) / (decay * wspd) )
  endif

end function sheltering_coare35

!===============================================================================

end module umwm_sheltering
