module umwm_source_functions

contains

!===============================================================================

subroutine sin_d12( i, iwsfc )

  use consts_coms,     only: pi2, grav, vonki, r8
  use umwm_sheltering, only: sheltering_coare35, sheltering_reynolds
  use mem_sea,         only: sea
  use mem_sfcg,        only: sfcg
  use umwm_module,     only: pm, om, oc, rhoswi, umwm, fprog, ucurr, vcurr, cp0, &
                             cth, sth, logl2overz, ssin, fkovg, sin_diss1, sin_diss2, &
                             opeak, ppeak, swh, dcp0, nu_air, freq
  implicit none

  ! Wind input function based on Jeffreys's sheltering hypothesis
  ! and described by Donelan et al. (2012).

  integer, intent(in) :: i, iwsfc

  integer :: o, p
  real    :: c7, s, fcutoff, shelt, a(2), sheltr, f, wspd10m, sd1
  real    :: wsl2(om), coswdir(pm), wcurr(pm), fs(pm)

  !!!!!!!!!!!!! temp
  integer, parameter :: nswmzons = 0
  !!!!!!!!!!!!!
  real,    parameter :: c6 = pi2 * rhoswi

  ! cut-off frequency (4*pierson-moskowitz peak frequency)
  fcutoff = min(fprog, 0.53 * grav / umwm%wspd10m(i) )

  ! compute variable sheltering coefficient
  shelt   = sheltering_coare35( umwm%wspd10m(i) )

  c7 = c6 * (1. - sea%seaicec(i)) * sfcg%rhos(iwsfc)

  ! Special for light winds
  if (umwm%wspd10m(i) <= 3.0) then
     sd1 = sin_diss2
  else
     f = min(umwm%wspd10m(i)-3., 7.0) / 7.
     sd1 = min(sin_diss1, shelt) * f + sin_diss2 * (1. - f)
  endif

  a = [umwm%uwind(i), umwm%vwind(i)] / umwm%wspd(i)

  do p = 1, pm
     coswdir(p) = a(1) * cth(p) + a(2) * sth(p)

     f     = 0.5*coswdir(p) + 0.5
     fs(p) = sin_diss2 * f + sd1 * (1. - f)
  enddo

  if (nswmzons > 0) then
     do p = 1, pm
        wcurr(p) = ucurr(i) * cth(p) + vcurr(i) * sth(p)
     enddo
  endif

  ! search for the cut-off frequency bin
  do o = om-2, 2, -1
     if (fcutoff > freq(o)) exit
  enddo
  oc(i) = o

  ! wind speed at half wavelength (todo: add M.O. stability functions psim)
  do o = 1, om
     wsl2(o) = umwm%wspd(i) + vonki * umwm%ustar(i) * logl2overz(o,i)
  enddo

  if (nswmzons > 0) then
     do o = 1, om
        do p = 1, pm
           s = wsl2(o) * coswdir(p) - cp0(o,i) - wcurr(p)
           ssin(p,o,i) = c7 * abs(s) * s * fkovg(o,i)
        enddo
     enddo
  else
     do o = 1, om
        do p = 1, pm
           s = wsl2(o) * coswdir(p) - cp0(o,i)
           ssin(p,o,i) = c7 * abs(s) * s * fkovg(o,i)
        enddo
     enddo
  endif

  ! prevent negative ssin for diagnostic tail
  do o = oc(i)+1, om
     do p = 1, pm
        ssin(p,o,i) = max(ssin(p,o,i), 0.) * shelt
     enddo
  enddo

  ! apply variable sheltering coefficient for prognostic part
  do o = 1, oc(i)
     do p = 1, pm
        ssin(p,o,i) = max( ssin(p,o,i), 0. ) * shelt &
                    + min( ssin(p,o,i), 0. ) * fs(p)
     enddo
  enddo

end subroutine sin_d12

!===============================================================================

subroutine sds_d12( i, dummy, sds )

  use umwm_module, only: om, pm, mss_fac, evs, cth2pp, k3dk, fsdsf, k4, &
                         sds_power, inv_sds_power, oc
  implicit none

  ! Wave dissipation function described by Donelan et al. (2012).

  integer, intent(in)  :: i
  real,    intent(out) :: dummy(pm,om), sds(pm,om)
  integer              :: o, p, pp
  real                 :: vect(pm,om)
  real,      parameter :: eps = 1.e3 * tiny(1.)**inv_sds_power

  do o = 1, om-1
     vect(1:pm,o) = 0.0

     do pp = 1, pm
        do p = 1, pm
           vect(p,o) = vect(p,o) + cth2pp(p,pp) * evs(pp,o,i)
        enddo
     enddo
  enddo

  dummy(1:pm,2) = vect(1:pm,1) * k3dk(1,i)

  do o = 3, om
     do p = 1, pm
        dummy(p,o) = dummy(p,o-1) + vect(p,o-1) * k3dk(o-1,i)
     enddo
  enddo

  dummy(1:pm,1) = 1.0

  do o = 2, om
     do p = 1, pm
        dummy(p,o) = (1. + mss_fac * dummy(p,o))**2
     enddo
  enddo

  do o = 1, oc(i)+2
     do p = 1, pm
        sds(p,o) = fsdsf(o) * dummy(p,o) * max(eps, evs(p,o,i) * k4(o,i))**sds_power
     enddo
  enddo

end subroutine sds_d12

!===============================================================================

subroutine snl_d12( i, iwsfc, snl, sdt, sds )

  use mem_sfcg,    only: sfcg
  use umwm_module, only: pm, om, oc, sdt_fac, rhoswi, sds, evs, k, sdt, &
                         bf1_renorm, bf2_renorm, snl_fac, cothkd, umwm, dta
  implicit none

  integer, intent(in   ) :: i, iwsfc
  real,    intent(out  ) :: snl(pm,om)
  real,    intent(out  ) :: sdt(om)
  real,    intent(inout) :: sds(pm,om)

  integer             :: o, p
  real                :: f, snfdt, sdsdt, eloss(pm,om)
  real,     parameter :: c6 = sdt_fac * sqrt(rhoswi)

  ! spread wave energy to 2 next longer wavenumbers exponentially decaying
  ! as distance from donating wavenumber, and remove the energy from
  ! donating wavenumbers:

  snfdt = snl_fac * dta

  do o = 1, oc(i)+2
     do p = 1, pm
        sdsdt = snfdt * sds(p,o)
!!      eloss(p,o) = evs(p,o,i) * min(sdsdt, 0.95)     ! explicit differencing
!!      eloss(p,o) = evs(p,o,o) * (1. - exp(-sdsdt))   ! exponential solution
        eloss(p,o) = evs(p,o,i) * sdsdt / (1. + sdsdt) ! implicit differencing
     enddo
  enddo

  do o = 1, oc(i)
     do p = 1, pm

        ! Apply spreading loss/gain tendencies
        snl(p,o) = bf2_renorm(o,i) * eloss(p,o+2)  &
                 + bf1_renorm(o,i) * eloss(p,o+1)  &
                 -                   eloss(p,o)

        ! account for plunging breakers
        sds(p,o) = sds(p,o) * cothkd(o,i)

     enddo
  enddo

  ! compute dissipation due to turbulence

  f = c6 * sqrt(sfcg%rhos(iwsfc)) * umwm%ustar(i)

  do o = 1, om
     sdt(o) = f * k(o,i)
  enddo

end subroutine snl_d12

!===============================================================================

subroutine s_ice( i, sice )

  use mem_sea,     only: sea
  use umwm_module, only: om, cg0, fice_lth, swh

  implicit none

  ! Wave attenuation by sea ice, following Kohout et al. (2014).

  integer, intent(in)  :: i
  real,    intent(out) :: sice(om)
  integer              :: o
  real                 :: wvatt

  ! parameters from Kohout et al. 2014
  real,     parameter :: h_th =  3.0      ! [m]
! real,     parameter :: c1   = -5.35e-6  ! [m^-1]
  real,     parameter :: c1   = -1.0e-5   ! [m^-1]
  real,     parameter :: s1   =  2.0 * c1

  if (sea%seaicec(i) > fice_lth) then

     wvatt = sea%seaicec(i) * s1 * min(h_th, max(.01, swh(i)))

     do o = 1, om
        sice(o) = cg0(o,i) * wvatt
     enddo

  else

     sice(1:om) = 0.0

  endif

end subroutine s_ice

!===============================================================================

subroutine source( i, dummy, sds, snl, sdt, sice )

  use umwm_module, only: oc, om, pm, dta, ssin, sdv, evsf, evs, oneoverk4, &
                         fsdsf, cothkd, inv_sds_power
  implicit none

  integer, intent(in) :: i
  real,    intent(in) :: dummy(pm,om), sds(pm,om), snl(pm,om)
  real,    intent(in) :: sdt(om), sice(om)
  integer             :: o, p
  real                :: vec(om), src
  real,     parameter :: eps = 10.*tiny(1.)

  do o = 1, om
     vec(o) = sice(o) - sdt(o) - sdv(o,i)
  enddo

  ! Integrate source terms for the prognostic range (o <= ol)

  do o = 1, oc(i)
     do p = 1, pm

        src         = dta * (ssin(p,o,i) - sds(p,o) + vec(o))
        evsf(p,o,i) = (evs(p,o,i) + snl(p,o)) * (1. + max(src,0.)) / (1. - min(src,0.))

     enddo
  enddo

  ! Integrate source terms for the diagnostic range (o > ol)

  do o = oc(i)+1, om
     do p = 1, pm
        src = ssin(p,o,i) + vec(o)

        evsf(p,o,i) = 0.0
        if (src > eps) evsf(p,o,i) = oneoverk4(o,i) * &
             (src / (fsdsf(o) * dummy(p,o) * cothkd(o,i)))**inv_sds_power
     enddo
  enddo

  do o = 1, om
     do p = 1, pm
        evs(p,o,i) = 0.5 * (evs(p,o,i) + evsf(p,o,i))
     enddo
  enddo

end subroutine source

!===============================================================================

end module umwm_source_functions
