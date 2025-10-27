subroutine surface_turb_flux()

  use leaf_coms,   only: isfcl
  use mem_land,    only: omland, land
  use mem_sea,     only: omsea, sea
  use mem_basic,   only: theta, tair, rr_v
  use mem_ijtabs,  only: itab_w, jtab_w, jtw_prog
  use mem_sfcg,    only: itab_wsfc, sfcg, mwsfc
  use misc_coms,   only: iparallel, dtlm
  use mem_grid,    only: lsw, lpw, volti
  use mem_turb,    only: akm_sfc, vkm_sfc, ustar, sfluxt, sfluxr, arw_sfc, &
                         wstar, wtv0, pblh, moli, ustar_k, wtv0_k
  use consts_coms, only: grav, p00, rocp, cp, cpi, alvl, eps_virt, vonk, p00i
  use oname_coms,  only: nl
  use mem_para,    only: myrank
  use mem_tend,    only: thilt, rr_wt
! use mem_co2,     only: rr_co2, rr_co2t

  implicit none

  integer :: j,iw,isea
  integer :: ks,kw,ka
  integer :: jsfc, jasfc, iwsfc

  real :: exneri, dtl
  real :: airthetav

  real, parameter :: onethird = 1./3.

  dtl = dtlm

  !$omp parallel do private(iw,ks,kw,exneri)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

  ! This subroutine is now called only if ISFCL = 0, which is the no-LEAF option.
  ! Assign atmospheric surface fluxes here.
  ! Example with sensible flux = 250 W/m^2 and latent flux of 150 W/m^2:

  !  sfluxt = 250.
  !  sfluxr = 150. / alvl

     sfluxt(iw) = 0.
     sfluxr(iw) = 0. / alvl

     ustar (iw) = 0.1  ! Minimum value

     wstar (iw) = 0.
     wtv0  (iw) = 0.

     akm_sfc (:,iw) = 0.

     do ks = 1, lsw(iw)
        kw = ks + lpw(iw) - 1

        exneri = theta(kw,iw) / tair(kw,iw)

        ustar_k(ks,iw) = 0.1
        wtv0_k (ks,iw) = 0.0

        thilt(kw,iw) = thilt(kw,iw) + sfluxt(iw) * cpi * exneri * arw_sfc(ks,iw) * volti(kw,iw)
        rr_wt(kw,iw) = rr_wt(kw,iw) + sfluxr(iw)                * arw_sfc(ks,iw) * volti(kw,iw)
     enddo

     ! albedt (iw) = albedo
     ! rlongup(iw) = stefan * 288.15 ** 4  ! std msl temp; make decision to not
                                              ! run radiation if not running leaf?
  enddo
  !$omp end parallel do

end subroutine surface_turb_flux

!===============================================================================

subroutine stars( zts, rough, vels, rhos, wstar, air_thetav, can_thetav, &
                  vkmsfc, vkhsfc, ustar, ggaero )

  ! Subroutine stars computes surface heat and vapor fluxes and momentum drag
  ! coefficient from Louis (1981) equations

  use consts_coms, only: grav

  implicit none

  ! Input variables

  real, intent(in) :: zts        ! height above surface of {vels, ths, rrv} [m]
  real, intent(in) :: rough      ! surface roughness height [m]
  real, intent(in) :: vels       ! atmos near-surface wind speed [m/s]
  real, intent(in) :: rhos       ! atmos near-surface density [kg/m^3]
  real, intent(in) :: wstar      ! PBL free-convective velocity [m/s]
  real, intent(in) :: air_thetav ! atmos near-surface virt. pot. temp [K]
  real, intent(in) :: can_thetav ! canopy air virt. pot. temp [K]

  ! Output variables

  real, intent(out) :: vkmsfc    ! surface drag coefficient
  real, intent(out) :: vkhsfc    ! surface heat and vapor transfer coefficient
  real, intent(out) :: ustar     ! surface friction velocity
  real, intent(out) :: ggaero    ! bare ground conductance m/s

  ! Local parameters

  real, parameter :: ustmin = .05 ! lower bound on ustar (friction velocity)
  real, parameter :: ubmin  = .1  ! lower bound on wind speed

  ! Local variables

  real :: vels0  ! wind speed with minimum imposed [m/s]
  real :: ri     ! bulk richardson number, eq. 3.45 in Garratt
  real :: a2fm   ! drag coefficient for momentum
  real :: a2fh   ! drag coefficient for heat/tracers
  real :: vtscr

  ! Routine to compute Louis (1981) surface layer parameterization.

  ! Louis surface layer profiles:
  ! u*^2 = U(zobs)^2                 * a2 * Fm(Ri)
  ! u*t* = U(zobs) * (T(zobs)-T(z0)) * a2 * Fh(Ri)

  vels0 = max(ubmin, sqrt(vels*vels + wstar*wstar))

  ri = 2.0 * grav * zts * (air_thetav - can_thetav)  &
     / ( (air_thetav + can_thetav) * vels0 * vels0 )

  ! Get the Louis drag coefficients
  call a2fmfh(ri, zts, rough, a2fm, a2fh)

  ustar = max(ustmin, vels0 * sqrt(a2fm))

  vtscr  = rhos * zts * vels0
  vkmsfc = vtscr * a2fm   ! convert drag coefficient to diffusivity
  vkhsfc = vtscr * a2fh   ! convert drag coefficient to diffusivity

  ggaero = vels0 * a2fh   ! aerodynamic conductance, currently unused

end subroutine stars

!===============================================================================

subroutine a2fmfh(ri, zobs, zrough, a2fm, a2fh)

  use consts_coms, only: vonk

  implicit none

  real, intent(in)  :: ri     ! surface layer bulk Richardson number
  real, intent(in)  :: zobs   ! reference height
  real, intent(in)  :: zrough ! roughness height (level where wind speed -> 0)
  real, intent(out) :: a2fm   ! a2 * Fm, drag coefficient for momentum
  real, intent(out) :: a2fh   ! a2 * Fh, drag coefficient for heat/tracers
  real              :: a2, c2

  ! This routine computes the product a2*Fm and a2*Fh, which are the surface
  ! drag coefficients from the Louis (1979) model with some updates from
  ! ECMWF given in "A Short History of the Operational PBL Parameterization
  ! at ECMWF", by Louis, Tiedtke, and Geleyn.

  ! Louis profiles:
  ! u*^2 = U(zobs)^2                 * a2 * Fm(Ri)
  ! u*t* = U(zobs) * (T(zobs)-T(z0)) * a2 * Fh(Ri)

  a2 = ( vonk / log(zobs/zrough) )**2

  if (ri >= 0.) then              ! STABLE CASE

     c2 = sqrt(1. + 5. * ri)
     a2fm = a2 / (1. + 10. * ri / c2)
     a2fh = a2 / (1. + 15. * ri * c2)

  else                            ! UNSTABLE CASE

     c2 = ri / (1. + 75. * a2 * sqrt(-zobs*ri/zrough))
     a2fm = a2 * (1. - 10. * c2)
     a2fh = a2 * (1. - 15. * c2)

  endif

end subroutine a2fmfh

!===============================================================================

subroutine surface_cuparm_flux()

  use mem_cuparm,  only: conprr, iactcu
  use mem_ijtabs,  only: istp, mrl_begl
  use mem_sfcg,    only: itab_wsfc
  use consts_coms, only: cliq, cice, alli, t00
  use mem_basic,   only: tair, rho, rr_v
  use mem_sfcg,    only: mwsfc, itab_wsfc, sfcg
  use misc_coms,   only: iparallel, dtlm
  use therm_lib,   only: rhovsl
  use mem_para,    only: myrank

  implicit none

  integer :: iw
  integer :: iwsfc, j
  integer :: kw

  real :: dtl, airtempc, tempc, qpcp

  ! Subroutine to transfer atmospheric cumulus parameterization
  ! precipitation FLUX to surface cells

  if (mrl_begl(istp) > 0) then

     dtl = dtlm

     ! Transfer precipitation FLUX to SFC grid cells

     !$omp parallel do private(j,iw,kw,airtempc,tempc,qpcp) schedule(guided)
     do iwsfc = 2, mwsfc

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        do j = 1,itab_wsfc(iwsfc)%nwatm
           iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
           kw = itab_wsfc(iwsfc)%kwatm(j)

           if (iactcu(iw) > 0) then

              ! Compute air temperature in C

              airtempc = tair(kw,iw) - t00

              ! Estimate wet bulb temp using computation from subroutine
              ! each_column in micphys. Assume that convective precip reaches
              ! surface at this wet bulb temp.

              tempc = airtempc - min(25., max(0., &
                   700. * (rhovsl(airtempc) / real(rho(kw,iw)) - rr_v(kw,iw))))

              if (tempc < 0.) then
                 qpcp = cice * tempc
              else
                 qpcp = cliq * tempc + alli
              endif

              sfcg%pcpg (iwsfc) = sfcg%pcpg (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * dtl * conprr(iw)
              sfcg%qpcpg(iwsfc) = sfcg%qpcpg(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * dtl * conprr(iw) * qpcp
              sfcg%dpcpg(iwsfc) = sfcg%dpcpg(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * dtl * conprr(iw) * .001

           endif
        enddo

     enddo
     !$omp end parallel do

  endif

end subroutine surface_cuparm_flux

!==============================================================================

subroutine sfclyr_profile (vels, rhos, canexner, ustar, sfluxt, sfluxr, dzt_bot, zrough, wstar, &
                           cantheta, canthetav, canrrv, airthetav, &
                           zobs, wind_zobs, theta_zobs, rrv_zobs)

  ! This subroutine diagnoses wind speed, temperature, and vapor mixing ratio
  ! over a sfc grid cell at height zobs that is within the surface layer profile

  use consts_coms, only: grav, cp, vonk

  implicit none

  real, intent(in) :: vels
  real, intent(in) :: rhos
  real, intent(in) :: canexner
  real, intent(in) :: ustar
  real, intent(in) :: sfluxt
  real, intent(in) :: sfluxr
  real, intent(in) :: dzt_bot
  real, intent(in) :: zrough
  real, intent(in) :: wstar
  real, intent(in) :: cantheta
  real, intent(in) :: canthetav
  real, intent(in) :: canrrv
  real, intent(in) :: airthetav
  real, intent(in) :: zobs

  real, intent(inout) :: theta_zobs
  real, intent(inout) ::  wind_zobs
  real, intent(inout) ::   rrv_zobs

  real :: vels0, richnum, a2fm, a2fh, fact

  ! Louis profiles:
  ! u*^2 = U(zobs)^2                 * a2 * Fm(Ri)
  ! u*t* = U(zobs) * (T(zobs)-T(z0)) * a2 * Fh(Ri)

  real, parameter :: ubmin  = .1  ! lower bound on wind speed

  vels0 = max(ubmin, sqrt(vels*vels + wstar*wstar))

  richnum = 2.0 * grav * dzt_bot * (airthetav - canthetav)  &
          / ( (airthetav + canthetav) * vels0 * vels0 )

  call a2fmfh(richnum, zobs, zrough, a2fm, a2fh)

  wind_zobs = ustar / sqrt(a2fm)

  fact = 1. / (rhos * wind_zobs * a2fh)

  theta_zobs = cantheta - sfluxt * fact / (cp * canexner)
  rrv_zobs   = canrrv   - sfluxr * fact

end subroutine sfclyr_profile

