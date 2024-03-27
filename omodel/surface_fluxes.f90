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
  real :: airthetav, ufree

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

subroutine stars( zts, rough, vels, rhos, ufree, air_thetav, can_thetav, &
                  vkmsfc, vkhsfc, ustar, ggaero )

  ! Subroutine stars computes surface heat and vapor fluxes and momentum drag
  ! coefficient from Louis (1981) equations

  use consts_coms, only: vonk, grav

  implicit none

  ! Input variables

  real, intent(in) :: zts        ! height above surface of {vels, ths, rrv} [m]
  real, intent(in) :: rough      ! surface roughness height [m]
  real, intent(in) :: vels       ! atmos near-surface wind speed [m/s]
  real, intent(in) :: rhos       ! atmos near-surface density [kg/m^3]
  real, intent(in) :: ufree      ! surface layer free-convective velocity [m/s]
  real, intent(in) :: air_thetav ! atmos near-surface virt. pot. temp [K]
  real, intent(in) :: can_thetav ! canopy air virt. pot. temp [K]

  ! Output variables

  real, intent(out) :: vkmsfc    ! surface drag coefficient
  real, intent(out) :: vkhsfc    ! surface heat and vapor transfer coefficient
  real, intent(out) :: ustar     ! surface friction velocity
  real, intent(out) :: ggaero    ! bare ground conductance m/s

  ! Local parameters

  real, parameter :: b = 5.
  real, parameter :: csm = 7.5
  real, parameter :: csh = 5.
  real, parameter :: d = 5.
  real, parameter :: ustmin = .05 ! lower bound on ustar (friction velocity)
  real, parameter :: ubmin  = .1  ! lower bound on wind speed

  ! Local variables

  real :: vels0  ! wind speed with minimum imposed [m/s]
  real :: a2     ! drag coefficient in neutral conditions, here same for h/m
  real :: ri     ! bulk richardson number, eq. 3.45 in Garratt
  real :: c1
  real :: c2
  real :: c3
  real :: cm
  real :: ch
  real :: fm
  real :: fh
  real :: vtscr  ! ustar times density

  ! Routine to compute Louis (1981) surface layer parameterization.

  vels0 = max(vels,ubmin,ufree)

  a2 = ( vonk / log(zts / rough) ) ** 2
  c1 = a2 * vels0

  ri = 2.0 * grav * zts * (air_thetav - can_thetav)  &
     / ( (air_thetav + can_thetav) * vels0 * vels0 )

  if (ri > 0.) then

     fm = 1. / (1. + (2. * b * ri / sqrt(1. + d * ri)))
     fh = 1. / (1. + (3. * b * ri * sqrt(1. + d * ri)))

  else                            ! UNSTABLE CASE

     c2 = b * a2 * sqrt(zts / rough * (abs(ri)))
     cm = csm * c2
     ch = csh * c2
     fm = (1. - 2. * b * ri / (1. + 2. * cm))
     fh = (1. - 3. * b * ri / (1. + 3. * ch))

  endif

  ustar = max(ustmin,sqrt(c1 * vels0 * fm))
  c3 = c1 * fh / ustar

  vtscr = ustar * rhos

  vkmsfc =   vtscr * ustar * zts / vels0
  vkhsfc =   vtscr * c3 * zts

  ! Store the aerodynamic conductance between the surface canopy and
  ! the lowest model level

  ggaero  = c3 * ustar

end subroutine stars

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

subroutine sfclyr_profile (vels, rhos, canexner, ustar, sfluxt, sfluxr, dzt_bot, zrough, ufree, &
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
  real, intent(in) :: ufree
  real, intent(in) :: cantheta
  real, intent(in) :: canthetav
  real, intent(in) :: canrrv
  real, intent(in) :: airthetav
  real, intent(in) :: zobs

  real, intent(inout) :: theta_zobs 
  real, intent(inout) ::  wind_zobs
  real, intent(inout) ::   rrv_zobs

  real :: vels0, a2, richnum

  real, parameter :: ubmin  = .1  ! lower bound on wind speed 

  vels0 = max(vels,ubmin,ufree)

  a2 = (vonk / log(zobs / zrough)) ** 2

  richnum = 2.0 * grav * dzt_bot * (airthetav - canthetav)  &
          / ( (airthetav + canthetav) * vels0 * vels0 )

  if (airthetav >= canthetav) then

     wind_zobs = sqrt((ustar**2 / a2) &
               * (1. + 10. * richnum / sqrt(1. + 5. * richnum)) )

     theta_zobs = cantheta - (sfluxt / (cp * canexner * a2 * vels0 * rhos)) &
                * (1. + 15. * richnum / sqrt(1. + 5. * richnum))

     rrv_zobs = canrrv - (sfluxr / (a2 * vels0 * rhos)) &
                * (1. + 15. * richnum / sqrt(1. + 5. * richnum))

  else

     wind_zobs = sqrt((ustar**2 / a2) &
               / (1. - 10. * richnum / (1. + 75. * a2 * sqrt(-zobs * richnum/zrough))))

     theta_zobs = cantheta - (sfluxt / (cp * canexner * a2 * vels0 * rhos)) &
                / (1. - 15. * richnum / (1. + 75. * a2 * sqrt(-zobs * richnum / zrough)))

     rrv_zobs = canrrv - (sfluxr / (a2 * vels0 * rhos)) &
                / (1. - 15. * richnum / (1. + 75. * a2 * sqrt(-zobs * richnum / zrough)))

  endif

end subroutine sfclyr_profile

