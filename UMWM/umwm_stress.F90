module umwm_stress

  !! Module with functions and subroutines to evaluate stresses.

  use umwm_module
  use consts_coms, only: vonk
  use umwm_stokes, only: u_stokes => ustok, v_stokes => vstok
  use mem_sea, only: msea, omsea
  use mem_sfcg, only: sfcg

  implicit none

  private
  public :: stress

contains

!===============================================================================

subroutine stress(option)

  character(3), intent(in) :: option

  integer :: i, o, p, iwsfc

  real :: taux_util(om,msea), tauy_util(om,msea)
  real :: tail(msea)

  ! wind speed and direction relative to surface velocity
  real :: wspdrel(msea), wdirrel(msea)

  real :: cd_form(msea), cd_skin(msea)

  if (option == 'atm') then

     ! compute form drag from atmosphere: (positive into waves)

     !$omp parallel do private(iwsfc,o,p)
     do i = 2,msea
        iwsfc = i + omsea
        tail(i) = stress_tail(umwm%wspd(i), k(om,i))  ! evaluate wind speed dependent tail

        taux_util(:,i) = 0
        tauy_util(:,i) = 0

        ! Stress vector integral
        do p = 1,pm
           do o = 1,om
              taux_util(o,i) = taux_util(o,i) + evs(o,p,i) * ssin(o,p,i) * cth(p) / cp0(o,i)
              tauy_util(o,i) = tauy_util(o,i) + evs(o,p,i) * ssin(o,p,i) * sth(p) / cp0(o,i)
           enddo
        enddo

        ! compute the tail
        tailatmx(i) = taux_util(om,i) * k(om,i) * tail(i) * rhosw * dthg
        tailatmy(i) = tauy_util(om,i) * k(om,i) * tail(i) * rhosw * dthg

        taux_form(i) = sum(taux_util(:,i) * kdk(:,i)) * rhosw * dthg + tailatmx(i)
        tauy_form(i) = sum(tauy_util(:,i) * kdk(:,i)) * rhosw * dthg + tailatmy(i)

        taux_diag(i) = sum(taux_util(oc(i):om,i) * kdk(oc(i):om,i)) * rhosw * dthg
        tauy_diag(i) = sum(tauy_util(oc(i):om,i) * kdk(oc(i):om,i)) * rhosw * dthg

        ! wind speed and direction relative to surface velocity
        call wind_relative(umwm%wspd(i), umwm%wdir(i), &
                           ucurr(i) + u_stokes(1,i),      &
                           vcurr(i) + v_stokes(1,i),      &
                           wspdrel(i), wdirrel(i))

        ! form-induced drag coefficient
        cd_form(i) = drag_coefficient(taux_form(i), tauy_form(i), sfcg%rhos(iwsfc), umwm%wspd(i))

        ! skin-induced drag coefficient
        cd_skin(i) = drag_coefficient_skin(cd_form(i),umwm% wspd(i), wspdrel(i), sfcg%dzt_bot(iwsfc), nu_air, vonk)

        taux_skin(i) = sfcg%rhos(iwsfc) * cd_skin(i) * wspdrel(i)**2 * cos(wdirrel(i))
        tauy_skin(i) = sfcg%rhos(iwsfc) * cd_skin(i) * wspdrel(i)**2 * sin(wdirrel(i))

        umwm%taux(i) = taux_form(i) + taux_skin(i)
        umwm%tauy(i) = tauy_form(i) + tauy_skin(i)

        ! total (form + skin) drag coefficient
        umwm%cd(i) = drag_coefficient(umwm%taux(i), umwm%tauy(i), sfcg%rhos(iwsfc), umwm%wspd(i))

        ! Update friction velocity based on total (form + skin) stress
        umwm%ustar(i) = sqrt(sqrt(umwm%taux(i)**2 + umwm%tauy(i)**2) / sfcg%rhos(iwsfc))

     enddo
     !$omp end parallel do
  endif ! if(option=='atm')

  if (option == 'ocn') then

     ! compute stress into ocean top (positive into ocean)
     !$omp parallel do private(o,p)
     do i = 2,msea
        iwsfc = i + omsea
        tail(i) = stress_tail(umwm%wspd(i), k(om,i))  ! evaluate wind speed dependent tail

        taux_util(:,i) = 0.
        tauy_util(:,i) = 0.
        do p = 1, pm

           ! dissipation into currents
           do o = 1, om
              taux_util(o,i) = taux_util(o,i) + evs(o,p,i) * cth(p) * invcp0(o,i) &
                           * (sds(o,p,i) + sdt(o,i) + sdv(o,i))
                              
              tauy_util(o,i) = tauy_util(o,i) + evs(o,p,i) * sth(p) * invcp0(o,i) &
                           * (sds(o,p,i) + sdt(o,i) + sdv(o,i))
           enddo

           ! correction for the Snl in diagnostic part
           do o = 3, om
              taux_util(o,i) = taux_util(o,i) - snl(o,p,i) * cth(p) * invcp0(o,i) &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))

              tauy_util(o,i) = tauy_util(o,i) - snl(o,p,i) * sth(p) * invcp0(o,i) &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))
           enddo
        enddo

        ! compute the tail
        tailocnx(i) = taux_util(om,i) * k(om,i) * tail(i) * rhosw * dthg
        tailocny(i) = tauy_util(om,i) * k(om,i) * tail(i) * rhosw * dthg

        ! integrate over frequencies
        taux_ocntop(i) = sum(taux_util(:,i) * kdk(:,i)) * rhosw * dthg &
                       + tailocnx(i) + taux_skin(i)
        tauy_ocntop(i) = sum(tauy_util(:,i) * kdk(:,i)) * rhosw * dthg &
                       + tailocny(i) + tauy_skin(i)

        ! compute stress into ocean bottom (positive into ocean)
        taux_util(:,i) = 0.
        tauy_util(:,i) = 0.
        do p = 1, pm
           do o = 1, om
              taux_util(o,i) = taux_util(o,i) + evs(o,p,i) * sbf(o,i) * cth(p) * invcp0(o,i)
              tauy_util(o,i) = tauy_util(o,i) + evs(o,p,i) * sbf(o,i) * sth(p) * invcp0(o,i)
           enddo
        enddo

        taux_ocnbot(i) = sum(taux_util(:,i) * kdk(:,i)) * rhosw * dthg
        tauy_ocnbot(i) = sum(tauy_util(:,i) * kdk(:,i)) * rhosw * dthg

        ! Snl conserves energy, but not momentum. the following term is
        ! the momentum loss due to non-linear downshifting of energy.
        taux_util(:,i) = 0.
        tauy_util(:,i) = 0.
        do p = 1, pm
           do o = 3, oc(i)
              taux_util(o,i) = taux_util(o,i) + snl(o,p,i)              &
                           * (bf1 * (1 - cp0(o,i) * invcp0(o-1,i))  &
                           +  bf2 * (1 - cp0(o,i) * invcp0(o-2,i))) &
                           * cth(p) * invcp0(o,i)
              tauy_util(o,i) = tauy_util(o,i) + snl(o,p,i)              &
                           * (bf1 * (1 - cp0(o,i) * invcp0(o-1,i))  &
                           +  bf2 * (1 - cp0(o,i) * invcp0(o-2,i))) &
                           * sth(p) * invcp0(o,i)
           enddo
        enddo

        taux_snl(i) = sum(taux_util(:,i) * kdk(:,i)) * rhosw * dthg
        tauy_snl(i) = sum(tauy_util(:,i) * kdk(:,i)) * rhosw * dthg

        ! compute energy flux into ocean:
        taux_util(:,i) = 0.
        tauy_util(:,i) = 0.
        do p = 1, pm

           do o = 1, om
              taux_util(o,i) = taux_util(o,i) + evs(o,p,i) * sds(o,p,i) * cth(p)
              tauy_util(o,i) = tauy_util(o,i) + evs(o,p,i) * sds(o,p,i) * sth(p)
           enddo

           do o = 3, om
              taux_util(o,i) = taux_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * cth(p)
              tauy_util(o,i) = tauy_util(o,i) - snl(o,p,i)     &
                           * (bf1 * cp0(o,i) * invcp0(o-1,i) &
                           +  bf2 * cp0(o,i) * invcp0(o-2,i))&
                           * sth(p)
           enddo

        enddo

        epsx_ocn(i) = sum(taux_util(:,i) * kdk(:,i)) * rhosw * dthg
        epsy_ocn(i) = sum(tauy_util(:,i) * kdk(:,i)) * rhosw * dthg

        ! compute energy flux from air:
        taux_util(:,i) = 0.
        tauy_util(:,i) = 0.
        do p = 1, pm
           do o = 1, om
              taux_util(o,i) = taux_util(o,i) + evs(o,p,i) * ssin(o,p,i) * cth(p)
              tauy_util(o,i) = tauy_util(o,i) + evs(o,p,i) * ssin(o,p,i) * sth(p)
           enddo
        enddo

        epsx_atm(i) = sum(taux_util(:,i) * kdk(:,i), dim=1) * rhosw * dthg
        epsy_atm(i) = sum(tauy_util(:,i) * kdk(:,i), dim=1) * rhosw * dthg

        ! This part calculates the components of form drag.  It has nothing to do with
        ! momentum fluxes into ocean, but we do it here because it is needed only for output.
        taux1(i) = 0.
        tauy1(i) = 0.
        taux2(i) = 0.
        tauy2(i) = 0.
        taux3(i) = 0.
        tauy3(i) = 0.

        do p = 1, pm
           do o = 1, om

              dummy(o,p,i) = evs(o,p,i) * ssin(o,p,i) * invcp0(o,i) * kdk(o,i)

              if (ssin(o,p,i) > 0.) then ! positive stress
                 ! wind pushing waves
                 taux1(i) = taux1(i) + dummy(o,p,i) * cth(p)
                 tauy1(i) = tauy1(i) + dummy(o,p,i) * sth(p)
              else ! negative stress, two cases
                 if (cos(umwm%wdir(i)-th(p)) < 0.) then
                    ! waves against wind
                    taux2(i) = taux2(i) + dummy(o,p,i) * cth(p)
                    tauy2(i) = tauy2(i) + dummy(o,p,i) * sth(p)
                 else
                    ! waves overrunning wind
                    taux3(i) = taux3(i) + dummy(o,p,i) * cth(p)
                    tauy3(i) = tauy3(i) + dummy(o,p,i) * sth(p)
                 endif
              endif

           enddo
        enddo

        taux1(i) = taux1(i) * rhosw * dthg
        tauy1(i) = tauy1(i) * rhosw * dthg
        taux2(i) = taux2(i) * rhosw * dthg
        tauy2(i) = tauy2(i) * rhosw * dthg
        taux3(i) = taux3(i) * rhosw * dthg
        tauy3(i) = tauy3(i) * rhosw * dthg

     enddo
     !$omp end parallel do

  endif ! if(option=='ocn')

end subroutine stress

!===============================================================================

  real pure elemental function drag_coefficient(taux, tauy, rhoa, wspd) result(cd)
    !! Computes drag coefficient from stress vector (taux, tauy, N/m^2),
    !! air density (rhoa, kg/m^3), and wind speed (wspd, m/s).
    real, intent(in) :: taux, tauy, rhoa, wspd
    cd = sqrt(taux**2 + tauy**2) / (rhoa * wspd**2)
  end function drag_coefficient

!===============================================================================

  real pure elemental function drag_coefficient_skin(cd_form, wspd, wspdrel, z, air_viscosity, von_karman) result(cd)
    !! Computes the skin drag coefficient attenuated by from drag,
    !! given input absolute (wspd) and relative (wspdrel) wind speed (m/s),
    !! height z (m), air viscosity (m^2/s), and Von Karman constant.
    real, intent(in) :: wspd, z, air_viscosity, von_karman
    real, intent(in) :: cd_form, wspdrel

    cd = friction_velocity_skin(wspdrel, z, air_viscosity, von_karman)**2 / wspd**2
    cd = cd * (1. + 2. * cd / (cd + cd_form + tiny(cd))) / 3.
    if (cd > 1.e-2) cd = 1.e-2
  end function drag_coefficient_skin

!===============================================================================

  real pure elemental function friction_velocity_skin(wspd, z, air_viscosity, von_karman) result(ustar)
    !! Computes the skin friction velocity (flow over smooth surface) in m/s
    !! given input wind speed (m/s), height z (m), air viscosity (m^2/s),
    !! and Von Karman constant. Input wind speed should relative to the
    !! water surface.
    real, intent(in) :: z, air_viscosity, von_karman
    real, intent(in) :: wspd
    real :: z0
    integer :: n
    z0 = 1.e-3
    do n = 1, 6
      ustar = von_karman * wspd / log(z / z0)
      z0 = 0.132 * air_viscosity / ustar
    end do
  end function friction_velocity_skin

!===============================================================================

  pure elemental real function stress_tail(wspd, kmax) result(tail)
    !! Computes the wind speed-dependent stress tail
    !! from kmax (highest wavenumber bin) to capillary.
    real, intent(in) :: wspd, kmax
    real :: kmax_pow_tail
    real, parameter :: kcap = 1.e3 ! capillary wavenumber
    real, parameter :: a = 0.000112, b = -0.01451, c = -1.0186
    tail = a * wspd**2 + b * wspd + c
    kmax_pow_tail = kmax**tail
    tail = (kcap**(tail + 1.) - kmax * kmax_pow_tail)&
         / (kmax_pow_tail * (tail + 1.))
  end function stress_tail

!===============================================================================

  pure subroutine stress_vector_integral(evs, src, cp, taux, tauy)
    !! Integrates the input source function src over the spectrum with
    !! wave variance e over the directions.
    real, intent(in) :: evs(:,:,:), src(:,:,:), cp(:,:)
    real, intent(out) :: taux(:,:), tauy(:,:)
    integer :: i, p, o
    integer :: dim(3)
    dim = shape(evs)
    taux = 0.
    tauy = 0.
    do i = 1, dim(3)
      do p = 1, dim(2)
        do o = 1, dim(1)
          taux(o,i) = taux(o,i) + evs(o,p,i) * src(o,p,i) * cth(p) / cp(o,i)
          tauy(o,i) = tauy(o,i) + evs(o,p,i) * src(o,p,i) * sth(p) / cp(o,i)
        end do
      end do
    end do
  end subroutine stress_vector_integral

!===============================================================================

  pure elemental subroutine wind_relative(wspd, wdir, us, vs, wspdrel, wdirrel)
    !! Computes the wind speed and direction relative to the water surface,
    !! given input absolute wind speed and direction, and u- and v- components
    !! of the surface velocity.
    real, intent(in) :: wspd, wdir, us, vs
    real, intent(out) :: wspdrel, wdirrel
    real :: urel, vrel
    urel = wspd * cos(wdir) - us
    vrel = wspd * sin(wdir) - vs
    wspdrel = sqrt(urel**2 + vrel**2)
    wdirrel = atan2(vrel, urel)
  end subroutine wind_relative

end module umwm_stress
