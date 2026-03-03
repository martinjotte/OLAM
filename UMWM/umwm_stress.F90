module umwm_stress

  !! Module with functions and subroutines to evaluate stresses.

  implicit none

  private
  public :: stress_atm, friction_velocity_skin

contains

!===============================================================================

subroutine stress_atm(i,iwsfc)

  use umwm_module,   only: om, pm, umwm, k, evs, ssin, invcp0, cth, sth, rhosw, dthg, stokes, &
                           kdk, ucurr, vcurr, taux_form, tauy_form, taux_skin, tauy_skin, nu_air, th
  use umwm_stokes,   only: ustoksfc, vstoksfc
  use mem_sfcg,      only: sfcg
  use consts_coms,   only: vonk
  use umwm_oforcing, only: zo_andreas

  implicit none

  ! Compute form drag from atmosphere: (positive into waves)

  integer, intent(in) :: i, iwsfc
  integer             :: o, p
  real                :: taux_util(om), tauy_util(om)
  real                :: strs, wspdrel, wdirrel, uc, vc
  real                :: cd_form, cd_skin, tau, rho_dthg
  real                :: tail, tailatmx, tailatmy, taux_sum, tauy_sum
  integer,  parameter :: nswmzons = 0

  if (.not. umwm%iactive(i)) then

     call zo_andreas( umwm%wspd10m(i), umwm%alogzo(i) )
     umwm%ustar(i) = umwm%wspd(i) * vonk / (umwm%alogzs(i) - umwm%alogzo(i))
     return

  endif

  rho_dthg = rhosw * dthg

  ! Evaluate wind speed dependent tail

  tail = rho_dthg * stress_tail( umwm%wspd10m(i), k(om,i) )

  ! Stress vector integral

  taux_util(1:om) = 0.
  tauy_util(1:om) = 0.

  do o = 1, om
     do p = 1, pm
        strs = evs(p,o,i) * ssin(p,o,i)
        taux_util(o) = taux_util(o) + strs * cth(p)
        tauy_util(o) = tauy_util(o) + strs * sth(p)
     enddo
  enddo

  taux_util(1:om) = taux_util(1:om) * invcp0(1:om,i)
  tauy_util(1:om) = tauy_util(1:om) * invcp0(1:om,i)

  taux_sum = sum( taux_util(1:om) * kdk(1:om,i) )
  tauy_sum = sum( tauy_util(1:om) * kdk(1:om,i) )

  ! Compute the tail
  tailatmx = taux_util(om) * k(om,i) * tail
  tailatmy = tauy_util(om) * k(om,i) * tail

  taux_form(i) = taux_sum * rho_dthg + tailatmx
  tauy_form(i) = tauy_sum * rho_dthg + tailatmy

  if (stokes .or. nswmzons > 0) then

     uc = 0.
     vc = 0.

     if (nswmzons > 0) then
        uc = uc + ucurr(i)
        vc = vc + vcurr(i)
     endif

     if (stokes) then
        uc = uc + ustoksfc(i)
        vc = vc + vstoksfc(i)
     endif

     ! wind speed and direction relative to surface velocity
     call wind_relative(umwm%uwind(i), umwm%vwind(i), &
                        uc, vc, wspdrel, wdirrel)

  else

     wspdrel = umwm%wspd(i)
     wdirrel = umwm%wdir(i)

  endif

  ! form-induced drag coefficient
  cd_form = drag_coefficient(taux_form(i), tauy_form(i), sfcg%rhos(iwsfc), umwm%wspd(i))

  ! skin-induced drag coefficient
  cd_skin = drag_coefficient_skin(cd_form, umwm%wspd(i), wspdrel, sfcg%dzt_bot(iwsfc))

  tau          = sfcg%rhos(iwsfc) * cd_skin * wspdrel**2
  taux_skin(i) = tau * cos(wdirrel)
  tauy_skin(i) = tau * sin(wdirrel)

  umwm%taux(i) = taux_form(i) + taux_skin(i)
  umwm%tauy(i) = tauy_form(i) + tauy_skin(i)

  ! total (form + skin) drag coefficient
  ! umwm%cd(i) = drag_coefficient(umwm%taux(i), umwm%tauy(i), sfcg%rhos(iwsfc), umwm%wspd(i))

  ! Update friction velocity based on total (form + skin) stress
  umwm%ustar(i) = sqrt( sqrt(umwm%taux(i)**2 + umwm%tauy(i)**2) / sfcg%rhos(iwsfc) )

  ! Update roughness length
  umwm%alogzo(i) = umwm%alogzs(i) - umwm%wspd(i) * vonk / umwm%ustar(i)

end subroutine stress_atm

!===============================================================================

real pure elemental function drag_coefficient(taux, tauy, rhoa, wspd) result(cd)
  implicit none

  ! Computes drag coefficient from stress vector (taux, tauy, N/m^2),
  ! air density (rhoa, kg/m^3), and wind speed (wspd, m/s).

  real, intent(in) :: taux, tauy, rhoa, wspd

  cd = sqrt(taux**2 + tauy**2) / (rhoa * wspd**2)

end function drag_coefficient

!===============================================================================

real pure elemental function drag_coefficient_skin(cd_form, wspd, wspdrel, z) result(cd)
  implicit none

  ! Computes the skin drag coefficient attenuated by from drag,
  ! given input absolute (wspd) and relative (wspdrel) wind speed (m/s),
  ! height z (m), air viscosity (m^2/s), and Von Karman constant.

  real, intent(in) :: wspd, z
  real, intent(in) :: cd_form, wspdrel

  cd = friction_velocity_skin(wspdrel, z)**2 / wspd**2
  cd = cd * (1. + 2. * cd / (cd + cd_form + tiny(cd))) / 3.
  if (cd > 1.e-2) cd = 1.e-2

end function drag_coefficient_skin

!===============================================================================

real pure elemental function friction_velocity_skin(wspd, z) result(ustar)

  use consts_coms, only: vonk
  use umwm_module, only: nu_air

  implicit none

  ! Computes the skin friction velocity (flow over smooth surface) in m/s
  ! given input wind speed (m/s), height z (m), air viscosity (m^2/s),
  ! and Von Karman constant. Input wind speed should relative to the
  ! water surface.

  real, intent(in) :: wspd, z
  real             :: z0
  integer          :: n

  z0 = 1.e-3
  do n = 1, 6
     ustar = vonk * wspd / log(z / z0)
     z0 = 0.135 * nu_air / ustar
  enddo

end function friction_velocity_skin

!===============================================================================

pure elemental real function stress_tail(wspd, kmax) result(tail)
  implicit none

  ! Computes the wind speed-dependent stress tail
  ! from kmax (highest wavenumber bin) to capillary.

  real, intent(in) :: wspd, kmax
  real             :: kmax_pow_tail
  real,  parameter :: kcap = 1.e3 ! capillary wavenumber
  real,  parameter :: a = 0.000112, b = -0.01451, c = -1.0186

! tail = a * wspd**2 + b * wspd + c
  tail = c + wspd * (b + wspd * a)

  kmax_pow_tail = kmax**tail
  tail = (kcap**(tail + 1.) - kmax * kmax_pow_tail) &
       / (kmax_pow_tail * (tail + 1.))

end function stress_tail

!===============================================================================

pure elemental subroutine wind_relative(ua, va, us, vs, wspdrel, wdirrel)
  implicit none

  ! Computes the wind speed and direction relative to the water surface,
  ! given input absolute wind speed and direction, and u- and v- components
  ! of the surface velocity.

  real, intent(in)  :: ua, va, us, vs
  real, intent(out) :: wspdrel, wdirrel
  real              :: urel, vrel

  urel = ua - us
  vrel = va - vs
  wspdrel = sqrt(urel**2 + vrel**2)
  wdirrel = atan2(vrel, urel)

end subroutine wind_relative

!===============================================================================

end module umwm_stress
