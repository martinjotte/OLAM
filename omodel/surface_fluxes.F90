!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================

! TURB and CUPARM flux scheduling:

! 1. Cuparm gives precip RATES
! 2. Stars gives heat and vapor flux RATES
! 3. Fluxes computed once per interval of dtlong = dtlm(1) = dt_leaf = dt_sea,
!    which (as of July 2011) are hardwired to be all the same.
! 4. Fluxes converted to AMOUNTS TRANSFERED based on rate * dtlm(1) 
! 5. Each time leaf runs, it uses at once all it has gotten and zeroes xfer arrays
! 6. When fluxes are done, they are done for all MRL = 1 atm points

! RADIATIVE fluxes only:

! 1. Radiative fluxes transfer RATES only, with no timestep information

!----------------------------------------------------------------------------
! AS OF JULY 2011, SUBROUTINE SURFACE_TURB_FLUXP IS HARDWIRED FOR 
! DT_LEAF = DTLM(1) AND MRL_LEAF = 1, WHICH FOR A LONG TIME HAS BEEN 
! ASSIGNED IN SUBROUTINE MODSCHED AND USED SUCCESSFULLY.
! IF DT_LEAF EVER GETS LARGE ENOUGH TO CAUSE INSTABILITY, THE INSTABILITY
! SHOULD BE CONTROLLED BY INCREASING CAPACITANCE OF THE LEAF COMPONENTS, 
! NOT BY REDUCING TIMESTEP; I.E., WE WANT LEAF TO BE CAPABLE OF RUNNING ON
! AS LONG A TIMESTEP AS THE ATMOSPHERIC MODEL.
!----------------------------------------------------------------------------

subroutine surface_turb_flux(mrl)

  use leaf_coms,   only: isfcl
  use mem_land,    only: omland, land
  use mem_sea,     only: omsea, sea
  use mem_ijtabs,  only: itab_w, itabg_w, jtab_w, jtw_prog, jtw_wstn
  use mem_sfcg,    only: itab_wsfc, sfcg, mwsfc
  use misc_coms,   only: isubdomain, dtlm
  use mem_grid,    only: lsw, lpw, arw
  use mem_turb,    only: akm_sfc, vkm_sfc, ustar, sfluxt, sfluxr, &
                         sxfer_tk, sxfer_rk, wstar, wtv0, pblh, moli, &
                         akh_dzi, akhth_dzi, akhrv_dzi
  use mem_basic,   only: press, rho, theta, tair, rr_v, vxe, vye, vze
  use mem_micro,   only: rr_c
  use consts_coms, only: grav, p00, rocp, cp, alvl, eps_virt, vonk, p00i
  use oname_coms,  only: nl
  use mem_para,    only: myrank

  implicit none

  integer, intent(in) :: mrl

  integer :: j,iw,isea,iland
  integer :: ks,kw,ka
  integer :: jsfc, jasfc, iwsfc

  real :: exneri, dtl
  real :: canexner, canexneri, cantheta, canthetav
  real :: airthetav, ufree
  real :: shflx  ! Specified surface sensible heat flux for ISFCL = 0 case [W/m^2]
  real :: srflx  ! Specified surface latent heat flux for ISFCL = 0 case [W/m^2]

  real :: sfc_cantheta(mwsfc)

  real, parameter :: onethird = 1./3.

  dtl = dtlm(1)

  if (isfcl == 0) then

     ! ISFCL = 0 is the no-LEAF option.  Assign surface fluxes here, noting the
     ! following examples.  SHFLX has units of [K m/s kg/m^3].

     ! Default surface sensible and latent heat fluxes = 0 W/m^2

     shflx = 0. / cp
     srflx = 0. / alvl

     ! Example with sensible flux = 250 W/m^2 and latent flux of 150 W/m^2:

     !  shflx = 250. / cp
     !  shflx = 150. / alvl

     !$omp parallel do private(iw,ks,kw,exneri)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        sfluxt(iw) = 0.
        sfluxr(iw) = 0.
        ustar (iw) = 0.1  ! Minimum value

        wstar (iw) = 0.
        wtv0  (iw) = 0.

        akm_sfc (:,iw) = 0.

        do ks = 1,lsw(iw)
           kw = ks + lpw(iw) - 1

           exneri = theta(kw,iw) / tair(kw,iw)

           sxfer_tk(ks,iw) = dtl * shflx * exneri * (arw(kw,iw) - arw(kw-1,iw))

           sxfer_rk(ks,iw) = dtl * srflx          * (arw(kw,iw) - arw(kw-1,iw))
        enddo

        ! albedt (iw) = albedo
        ! rlongup(iw) = stefan * 288.15 ** 4  ! std msl temp; make decision to not
                                              ! run radiation if not running leaf?
     enddo
     !$omp end parallel do

     return
  endif

  ! ISFCL = 1 is the LEAF option...

  ! Reset to zero the atm values of VKM_SFC, USTAR, SFLUXT, and SFLUXR
  ! for same mrl's that fluxes will be evaluated for this time.
  ! THUS, VKM_SFC, USTAR, SFLUXT, AND SFLUXR ARE ONLY SUMMED OVER
  ! SPACE, BUT NOT OVER TIME. (ON THE OTHER HAND, SXFER_TK AND SXFER_RK ARE SUMMED
  ! OVER BOTH SPACE AND TIME; THEY ARE RESET TO ZERO IN THILTEND_LONG AND
  ! SCALAR_TRANSPORT AFTER THEY ARE TRANSFERRED TO THE ATMOSPHERE.)

  ! Set sea and land fluxes to be done for SURFACE SIMILARITY:
  !    Do fluxes at beginning of long timestep at given mrl

  ustar   = 0.
  wtv0    = 0.
  sfluxt  = 0.
  sfluxr  = 0.
  vkm_sfc = 0.
  akm_sfc = 0.

  sxfer_tk = 0.
  sxfer_rk = 0.
! sxfer_ck = 0.  ! placeholder for CO2

  ! Loop over all SFC grid cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
  ! so that all surface fluxes are computed beneath all ATM columns that are primary

  !$omp parallel
  !$omp do private(airthetav,canexner,canexneri,cantheta,canthetav,ufree,isea)
  do iwsfc = 2,mwsfc

     if (isubdomain == 1) then
        if (.not. any( itab_w( max(1,itab_wsfc(iwsfc)%iwatm( 1:itab_wsfc(iwsfc)%nwatm )) )%irank == myrank ) ) cycle
     endif

     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))
     canexner  = (sfcg%prss(iwsfc) * p00i) ** rocp
     canexneri = 1. / canexner

     ! Compute turbulent fluxes based on whether SFC grid cell is land, lake, or sea

     if (sfcg%leaf_class(iwsfc) >= 1) then

        ! This is a land or lake cell

        cantheta  = sfcg%cantemp(iwsfc) * canexneri
        canthetav = cantheta * (1.0 + eps_virt * sfcg%canrrv(iwsfc))

        ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sfcg%wthv(iwsfc),0.0) / airthetav) ** onethird

        call stars(sfcg%dzt_bot (iwsfc), &
                   sfcg%rough   (iwsfc), &
                   sfcg%vels    (iwsfc), &
                   sfcg%rhos    (iwsfc), &
                   ufree               , &
                   sfcg%airtheta(iwsfc), &
                   airthetav           , &
                   sfcg%airrrv  (iwsfc), &
                   cantheta            , &
                   canthetav           , &
                   sfcg%canrrv  (iwsfc), &
                   sfcg%vkmsfc  (iwsfc), &
                   sfcg%sfluxt  (iwsfc), &
                   sfcg%sfluxr  (iwsfc), &
                   sfcg%ustar   (iwsfc), &
                   sfcg%ggaer   (iwsfc)  )

     else

        ! This is sea cell.  First, compute turbulent fluxes over open water areas.

        isea = iwsfc - omsea

        cantheta  = sea%sea_cantemp(isea) * canexneri
        canthetav = cantheta * (1.0 + eps_virt * sea%sea_canrrv(isea))

        ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sea%sea_wthv(isea),0.0) / airthetav) ** onethird

        call stars(sfcg%dzt_bot  (iwsfc), &
                   sea%sea_rough  (isea), &
                   sfcg%vels     (iwsfc), &
                   sfcg%rhos     (iwsfc), &
                   ufree                , &
                   sfcg%airtheta (iwsfc), &
                   airthetav            , &
                   sfcg%airrrv   (iwsfc), &
                   cantheta             , &
                   canthetav            , &
                   sea%sea_canrrv (isea), &
                   sea%sea_vkmsfc (isea), &
                   sea%sea_sfluxt (isea), &
                   sea%sea_sfluxr (isea), &
                   sea%sea_ustar  (isea), &
                   sea%sea_ggaer  (isea)  )

        sea%sea_wthv(isea) = ( sea%sea_sfluxt(isea) * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
           + sea%sea_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

        ! When we have CO2:
!       sea%sea_sfluxc(isea) = sfcg%rhos(iwsfc) * sea%sea_ggaer(isea) &
!                            * (sea%sea_co2(isea) - air_co2)

        ! Flux contributions to water
        sea%sea_sxfer_t(isea) = dtl * sea%sea_sfluxt(isea) * canexner
        sea%sea_sxfer_r(isea) = dtl * sea%sea_sfluxr(isea)
!       sea%sea_sxfer_c(isea) = dtl * sea%sea_sfluxc(isea)

        ! Check if sea ice is present

        if (sea%nlev_seaice(isea) == 0) then

           ! If no sea ice is present in this cell, zero out fluxes over ice
           ! and set cell flux to open water part

           sfcg%vkmsfc(iwsfc) = sea%sea_vkmsfc(isea)
           sfcg%ustar (iwsfc) = sea%sea_ustar (isea)
           sfcg%ggaer (iwsfc) = sea%sea_ggaer (isea)
           sfcg%sfluxt(iwsfc) = sea%sea_sfluxt(isea)
           sfcg%sfluxr(iwsfc) = sea%sea_sfluxr(isea)
!          sfcg%sfluxc(iwsfc) = sea%sea_sfluxc(isea)
           sfcg%wthv  (iwsfc) = sea%sea_wthv  (isea)

           sea%ice_vkmsfc(isea) = 0.0
           sea%ice_ustar (isea) = 0.0
           sea%ice_ggaer (isea) = 0.0
           sea%ice_sfluxt(isea) = 0.0
           sea%ice_sfluxr(isea) = 0.0
!          sea%ice_sfluxc(isea) = 0.0

           sea%ice_sxfer_t(isea) = 0.0
           sea%ice_sxfer_r(isea) = 0.0
!          sea%ice_sxfer_c(isea) = 0.0

        else

           ! If sea ice is present in this cell, compute turbulent fluxes over ice

           cantheta  = sea%ice_cantemp(isea) * canexneri
           canthetav = cantheta * (1.0 + eps_virt * sea%ice_canrrv(isea))

           ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sea%ice_wthv(isea),0.0) / airthetav) ** onethird

           call stars(sfcg%dzt_bot  (iwsfc), &
                      sea%ice_rough  (isea), &
                      sfcg%vels     (iwsfc), &
                      sfcg%rhos     (iwsfc), &
                      ufree                , &
                      sfcg%airtheta (iwsfc), &
                      airthetav            , &
                      sfcg%airrrv   (iwsfc), &
                      cantheta             , &
                      canthetav            , &
                      sea%ice_canrrv (isea), &
                      sea%ice_vkmsfc (isea), &
                      sea%ice_sfluxt (isea), &
                      sea%ice_sfluxr (isea), &
                      sea%ice_ustar  (isea), &
                      sea%ice_ggaer  (isea)  )

           ! When we have CO2:
!          sea%ice_sfluxc(isea) = sfcg%rhos(iwsfc) * sea%ice_ggaer(isea) &
!                               * (sea%ice_co2(isea) - airco2)

           sea%ice_wthv(isea) = ( sea%ice_sfluxt(isea) * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
              + sea%ice_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

           ! Flux contributions to seaice
           sea%ice_sxfer_t(isea) = dtl * sea%ice_sfluxt(isea) * canexner
           sea%ice_sxfer_r(isea) = dtl * sea%ice_sfluxr(isea)
!          sea%sea_sxfer_c(isea) = dtl * sea%sea_sfluxc(isea)

           ! Combine sea and ice values based on ice fraction:

           sfcg%vkmsfc(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_vkmsfc(isea) &
                                     + sea%seaicec(isea)  * sea%ice_vkmsfc(isea)

           sfcg%ustar(iwsfc)  = (1.0 - sea%seaicec(isea)) * sea%sea_ustar(isea) &
                                     + sea%seaicec(isea)  * sea%ice_ustar(isea)

           sfcg%ggaer(iwsfc)  = (1.0 - sea%seaicec(isea)) * sea%sea_ggaer(isea) &
                                     + sea%seaicec(isea)  * sea%ice_ggaer(isea)

           sfcg%sfluxt(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxt(isea) &
                                     + sea%seaicec(isea)  * sea%ice_sfluxt(isea)

           sfcg%sfluxr(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxr(isea) &
                                    + sea%seaicec(isea)  * sea%ice_sfluxr(isea)

           sfcg%wthv(iwsfc)   = (1.0 - sea%seaicec(isea)) * sea%sea_wthv(isea) &
                                     + sea%seaicec(isea)  * sea%ice_wthv(isea)

!          sfcg%sfluxc(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxc(isea) &
!                                    + sea%seaicec(isea)  * sea%ice_sfluxc(isea)

        endif

     endif  ! if this is sea cell

     sfcg%sxfer_t(iwsfc) = dtl * sfcg%sfluxt(iwsfc) * canexner
     sfcg%sxfer_r(iwsfc) = dtl * sfcg%sfluxr(iwsfc)
!    sfcg%sxfer_c(iwsfc) = dtl * sfcg%sfluxc(iwsfc)

  enddo
  !$omp end do

  ! Loop over ATM grid columns that are primary in this subdomain

  !$omp do private(iw, jsfc, iwsfc, jasfc, kw, ka, ks)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ! Loop over all SFC grid cells that couple to this ATM grid column

     do jsfc = 1,itab_w(iw)%jsfc2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)
        ka = lpw(iw)
        ks = kw - ka + 1

        ! Add flux contributions to IW atmospheric column

        ustar  (iw) = ustar  (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%ustar (iwsfc)
        wtv0   (iw) = wtv0   (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%wthv  (iwsfc)
        sfluxt (iw) = sfluxt (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%sfluxt(iwsfc)
        sfluxr (iw) = sfluxr (iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%sfluxr(iwsfc)
        vkm_sfc(iw) = vkm_sfc(iw) + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfcg%vkmsfc(iwsfc)
                    ! * land%slope_fact(iwsfc-omland)

        akm_sfc(ks,iw) = akm_sfc(ks,iw) + itab_wsfc(iwsfc)%arc(jasfc) * sfcg%vkmsfc(iwsfc)
                       ! * land%slope_fact(iwsfc-omland)

        sxfer_tk(ks,iw) = sxfer_tk(ks,iw) &
                        + itab_wsfc(iwsfc)%arc(jasfc) * dtl * sfcg%sfluxt(iwsfc)

        sxfer_rk(ks,iw) = sxfer_rk(ks,iw) &
                        + itab_wsfc(iwsfc)%arc(jasfc) * dtl * sfcg%sfluxr(iwsfc)

!       sxfer_ck(ks,iw) = sxfer_ck(ks,iw) &
!                       + itab_wsfc(iwsfc)%arc(jasfc) * dtl * sfcg%sfluxc(iwsfc)
     enddo
  enddo
  !$omp end do

  ! Compute some derived surface quantities

  !$omp do private(iw,ka)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ka = lpw(iw)

     moli(iw) = - grav * vonk * wtv0(iw) /  &
                ( ustar(iw)**3 * theta(ka,iw) * (1.0 + eps_virt * rr_v(ka,iw)) )

     if (wtv0(iw) > 0.0) then
        wstar(iw) = (grav * pblh(iw) * wtv0(iw) / theta(ka,iw)) ** onethird
     else
        wstar(iw) = 0.0
     endif

  enddo
  !$omp end do
  !$omp end parallel

  ! No MPI send/recv communication of ATM turbulent fluxes are required for
  ! averaging to SFC grid cells because the fluxes are directly computed on
  ! SFC grid cells.

end subroutine surface_turb_flux

!===============================================================================

subroutine stars( zts, rough, vels, rhos, ufree,  &
                  air_theta, air_thetav, air_rrv, &
                  can_theta, can_thetav, can_rrv, &
                  vkmsfc, sfluxt, sfluxr, ustar,  &
                  ggaero                          )

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
  real, intent(in) :: air_theta  ! atmos near-surface pot. temp [K]
  real, intent(in) :: air_thetav ! atmos near-surface virt. pot. temp [K]
  real, intent(in) :: air_rrv    ! atmos near-surface vapor spec hum [kg_vap/m^3]
  real, intent(in) :: can_theta  ! canopy air pot. temp [K]
  real, intent(in) :: can_thetav ! canopy air virt. pot. temp [K]
  real, intent(in) :: can_rrv    ! canopy air vapor spec hum [kg_vap/m^3]

  ! Output variables

  real, intent(out) :: vkmsfc    ! surface drag coefficient for this flux cell
  real, intent(out) :: sfluxt    ! surface sensible heat flux for this flux cell
  real, intent(out) :: sfluxr    ! surface vapor flux for this flux cell
  real, intent(out) :: ustar     ! surface friction velocity for this flux cell
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
  real :: tstar  !
  real :: rstar  !
  real :: vtscr  ! ustar times density

  ! Routine to compute Louis (1981) surface layer parameterization.

  vels0 = max(vels,ubmin,ufree)

  a2 = ( vonk / log(zts / rough) ) ** 2
  c1 = a2 * vels0

  ri = 2.0 * grav * zts * (air_thetav - can_thetav)  &
     / ( (air_thetav + can_thetav) * vels0 * vels0 )

  if (air_thetav >= can_thetav) then

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
  tstar = c3 * (air_theta - can_theta)
  rstar = c3 * (air_rrv   - can_rrv)

  vtscr = ustar * rhos

  vkmsfc =   vtscr * ustar * zts / vels0
  sfluxt = - vtscr * tstar
  sfluxr = - vtscr * rstar

  ! Store the aerodynamic conductance between the surface canopy and
  ! the lowest model level

  ggaero  = c3 * ustar

end subroutine stars

!===============================================================================

subroutine surface_cuparm_flux()

  use mem_cuparm,  only: conprr, iactcu
  use mem_ijtabs,  only: istp, mrl_begl, itabg_w
  use mem_sfcg,    only: itab_wsfc
  use consts_coms, only: cliq, cice, alli, t00
  use mem_basic,   only: tair, rho, rr_v
  use leaf_coms,   only: mrl_leaf
  use mem_sfcg,    only: mwsfc, itab_wsfc, sfcg
  use misc_coms,   only: io6, isubdomain, dtlm
  use therm_lib,   only: rhovsl
  use mem_para,    only: myrank

  implicit none

  integer :: mrl
  integer :: iw
  integer :: iwsfc, j
  integer :: kw

  real :: dtl, airtempc, tempc, qpcp

  ! Subroutine to transfer atmospheric cumulus parameterization 
  ! precipitation FLUX to surface cells

  ! Set land fluxes to be done for CUPARM
  !    1. for mrl = 0, no fluxes 
  !    2. For mrl > mrl_leaf, no fluxes
  !    3. For mrl <= mrl_leaf, do fluxes at beginning of long timestep 
  !                            for all points in mrl = 1

  mrl = mrl_begl(istp)
  if (mrl > 0) then

     dtl = dtlm(1)

     ! Transfer precipitation FLUX to SFC grid cells

     !$omp parallel do private(j, iw, kw, airtempc, tempc, qpcp) 
     do iwsfc = 2, mwsfc

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        do j = 1,itab_wsfc(iwsfc)%nwatm
           iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
           kw = itab_wsfc(iwsfc)%kwatm(j)

           ! Compute air temperature in C

           airtempc = tair(kw,iw) - t00

           ! Estimate wet bulb temp using computation from subroutine each_column in micphys.
           ! Assume that convective precip reaches surface at this wet bulb temp.

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
        enddo

     enddo
     !$omp end parallel do

  endif

end subroutine surface_cuparm_flux

