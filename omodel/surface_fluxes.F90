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

! TURB, CUPARM, and PRECIP flux scheduling:

! 1. Cuparm and micphys give precip RATES
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

use leaf_coms,   only: mwl, isfcl
use sea_coms,    only: mws
use mem_ijtabs,  only: itab_w, itabg_w, jtab_w, jtw_prog, jtw_wstn
use misc_coms,   only: isubdomain, dtlm
use mem_grid,    only: lsw, lpw, dzt_bot, arw
use mem_sea,     only: sea, itab_ws
use mem_leaf,    only: land, itab_wl
use mem_turb,    only: vkm_sfc, sfluxt, sfluxr, sxfer_tk, sxfer_rk, &
                       ustar, wstar, wtv0, pblh
use mem_basic,   only: press, rho, theta, tair, sh_v, vxe, vye, vze
use mem_micro,   only: sh_c
use consts_coms, only: grav, p00, rocp, cpi, eps_virt, vonk

implicit none

integer, intent(in) :: mrl

integer :: j,iw,iws,iwl
integer :: ks,kw,ka
integer :: nsea,jws,nland,jwl


real :: area_dt
real :: arf_kw
real :: arf_iw
real :: exneri
real :: my_co2
real :: canpress, canexneri, cantheta, canthetav
real :: airthetav, airshv, ufree, exner, vels, rhos
real :: shflx

real, parameter :: onethird = 1./3.

if (isfcl == 0) then

! ISFCL = 0 is the no-LEAF option.  Assign surface fluxes here, noting the
! following examples.  SHFLX has units of [K m/s kg/m^3].

! Default surface flux = 0 W/m^2

   shflx = 0. * cpi

! Example with surface flux = 250 W/m^2

!  shflx = 250. * cpi

!-----------------------------------------------------------------------------
   !$omp parallel do private(iw,ks,kw,exneri)
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!-----------------------------------------------------------------------------

      sfluxt(iw) = 0.
      sfluxr(iw) = 0.
      ustar (iw) = .1  ! Minimum value

      wstar (iw) = 0.
      wtv0  (iw) = 0.

      do ks = 1,lsw(iw)
         kw = ks + lpw(iw) - 1

         vkm_sfc (ks,iw) = 0.
         sxfer_tk(ks,iw) = 0.
         sxfer_rk(ks,iw) = 0.

         exneri = theta(kw,iw) / tair(kw,iw)

         sxfer_tk(ks,iw) = dtlm(1) * shflx * exneri &
                         * (arw(kw,iw) - arw(kw-1,iw))

      enddo

!      albedt (iw) = albedo
!      rlongup(iw) = stefan * 288.15 ** 4  ! std msl temp; make decision to not
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

!----------------------------------------------------------------------
!$omp parallel do private(iw)
do j = 1,jtab_w(jtw_wstn)%jend(mrl); iw = jtab_w(jtw_wstn)%iw(j)
!----------------------------------------------------------------------

   vkm_sfc(:,iw) = 0.
   ustar    (iw) = 0.
   wtv0     (iw) = 0.
   sfluxt   (iw) = 0.
   sfluxr   (iw) = 0.

enddo
!$omp end parallel do

! Fluxes with SEA cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and SEA cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in SEA cells

!----------------------------------------------------------------------
!dir$ novector
!$omp parallel do private(iw,kw,airshv,airthetav,vels,rhos,canpress,&
!$omp                     canexneri,cantheta,canthetav,ufree,exner)
do iws = 2,mws

   iw = itab_ws(iws)%iw
   if (isubdomain == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif
   kw = itab_ws(iws)%kw
!----------------------------------------------------------------------

! Calculate turbulent fluxes between atmosphere and "water canopy"

   if (allocated(sh_c)) then
      airshv = sh_v(kw,iw) + sh_c(kw,iw)
   else
      airshv = sh_v(kw,iw)
   endif

   airthetav = theta(kw,iw) * (1.0 + eps_virt * sh_v(kw,iw))
   vels      = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   rhos      = rho(kw,iw)

   canpress  = press(kw,iw) + dzt_bot(kw) * rho(kw,iw) * grav
   canexneri = (p00 / canpress) ** rocp
   cantheta  = sea%sea_cantemp(iws) * canexneri
   canthetav = cantheta * (1.0 + eps_virt * sea%sea_canshv(iws))

   ufree = (grav / airthetav * max(sea%sea_wthv(iws),0.0) * dzt_bot(kw)) ** onethird

! Sea (open water) component always computed

   call stars(dzt_bot         (kw), &
              sea%sea_rough  (iws), &
              vels                , &
              rhos                , &
              ufree               , &
              theta        (kw,iw), &
              airthetav           , &
              airshv              , &
              cantheta            , &
              canthetav           , &
              sea%sea_canshv (iws), &
              sea%sea_vkmsfc (iws), &
              sea%sea_sfluxt (iws), &
              sea%sea_sfluxr (iws), & 
              sea%sea_ustar  (iws), &
              sea%sea_rib    (iws), &
              sea%sea_ggaer  (iws)  )

   sea%sea_wthv(iws) = ( sea%sea_sfluxt(iws) * (1.0 + eps_virt * sh_v(kw,iw)) &
                       + sea%sea_sfluxr(iws) * eps_virt * theta(kw,iw) ) / rhos

   sea%sea_zeta(iws) = dzt_bot(kw) * grav * vonk * sea%sea_wthv(iws)  &
                     / (airthetav * sea%sea_ustar(iws) ** 3)

 ! When we have CO2:
 ! sea%sea_sfluxc(iws) = rhos * sea_ggaer(iws) * (sea%sea_co2 - air_co2)
 ! sea%sea_sfluxc(iws) = 0.

! Include fractional seaice component if seaice layers exist

   if (sea%nlev_seaice(iws) > 0) then

      cantheta  = sea%ice_cantemp(iws) * canexneri
      canthetav = cantheta * (1.0 + eps_virt * sea%ice_canshv(iws))

      ufree = (grav / airthetav * max(sea%ice_wthv(iws),0.0) * dzt_bot(kw)) ** onethird

      call stars(dzt_bot         (kw), &
                 sea%ice_rough  (iws), &
                 vels                , &
                 rhos                , &
                 ufree               , &
                 theta        (kw,iw), &
                 airthetav           , &
                 airshv              , &
                 cantheta            , &
                 canthetav           , &
                 sea%ice_canshv (iws), &
                 sea%ice_vkmsfc (iws), &
                 sea%ice_sfluxt (iws), &
                 sea%ice_sfluxr (iws), &
                 sea%ice_ustar  (iws), &
                 sea%ice_rib    (iws), &
                 sea%ice_ggaer  (iws)  )

      sea%ice_wthv(iws) = ( sea%ice_sfluxt(iws) * (1.0 + eps_virt * sh_v(kw,iw)) &
                          + sea%ice_sfluxr(iws) * eps_virt * theta(kw,iw) ) / rhos

      sea%ice_zeta(iws) = dzt_bot(kw) * grav * vonk * sea%ice_wthv(iws)  &
                        / (airthetav * sea%ice_ustar(iws) ** 3)

    ! When we have CO2:
    ! sea%ice_sfluxc(iws) = rhos * sea_ggaer(iws) * (sea%ice_co2 - air_co2)
    ! sea%ice_sfluxc(iws) = 0.

      ! Combine sea and ice values based on ice fraction:

      sea%vkmsfc(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_vkmsfc(iws) &
                             + sea%seaicec(iws)  * sea%ice_vkmsfc(iws)

      sea%ustar(iws)  = (1.0 - sea%seaicec(iws)) * sea%sea_ustar(iws) &
                             + sea%seaicec(iws)  * sea%ice_ustar(iws) 

      sea%ggaer(iws)  = (1.0 - sea%seaicec(iws)) * sea%sea_ggaer(iws) &
                             + sea%seaicec(iws)  * sea%ice_ggaer(iws) 

      sea%sfluxt(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfluxt(iws) &
                             + sea%seaicec(iws)  * sea%ice_sfluxt(iws)

      sea%sfluxr(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfluxr(iws) &
                             + sea%seaicec(iws)  * sea%ice_sfluxr(iws)

      sea%wthv(iws)   = (1.0 - sea%seaicec(iws)) * sea%sea_wthv(iws) &
                             + sea%seaicec(iws)  * sea%ice_wthv(iws)

      sea%zeta(iws)   = (1.0 - sea%seaicec(iws)) * sea%sea_zeta(iws) &
                             + sea%seaicec(iws)  * sea%ice_zeta(iws)

      sea%rib(iws)    = (1.0 - sea%seaicec(iws)) * sea%sea_rib(iws) &
                             + sea%seaicec(iws)  * sea%ice_rib(iws)

    ! sea%sfluxc(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfluxc(iws) &
    !                        + sea%seaicec(iws)  * sea%ice_sfluxc(iws)

   else

      sea%vkmsfc(iws) = sea%sea_vkmsfc(iws)
      sea%ustar (iws) = sea%sea_ustar (iws)
      sea%ggaer (iws) = sea%sea_ggaer (iws)
      sea%sfluxt(iws) = sea%sea_sfluxt(iws)
      sea%sfluxr(iws) = sea%sea_sfluxr(iws)
      sea%wthv  (iws) = sea%sea_wthv  (iws)
      sea%zeta  (iws) = sea%sea_zeta  (iws)
      sea%rib   (iws) = sea%sea_rib   (iws)
    ! sea%sfluxc(iws) = sea%sea_sfluxc(iws)

      sea%ice_vkmsfc(iws) = 0.0
      sea%ice_ustar (iws) = 0.0
      sea%ice_ggaer (iws) = 0.0
      sea%ice_sfluxt(iws) = 0.0
      sea%ice_sfluxr(iws) = 0.0
    ! sea%ice_sfluxc(iws) = 0.0
      
   endif

! Store atmospheric properties and flux contributions in SEA cell

   exner                = tair(kw,iw) / theta(kw,iw)

   sea%sxfer_t    (iws) = dtlm(1) * sea%sfluxt    (iws) * exner
   sea%sea_sxfer_t(iws) = dtlm(1) * sea%sea_sfluxt(iws) * exner
   sea%ice_sxfer_t(iws) = dtlm(1) * sea%ice_sfluxt(iws) * exner

   sea%sxfer_r    (iws) = dtlm(1) * sea%sfluxr    (iws)
   sea%sea_sxfer_r(iws) = dtlm(1) * sea%sea_sfluxr(iws)
   sea%ice_sxfer_r(iws) = dtlm(1) * sea%ice_sfluxr(iws)

 ! sea%sxfer_c    (iws) = dtlm(1) * sea%sfluxc    (iws)
 ! sea%sea_sxfer_c(iws) = dtlm(1) * sea%sea_sfluxc(iws)
 ! sea%ice_sxfer_c(iws) = dtlm(1) * sea%ice_sfluxc(iws)

enddo
!$omp end parallel do

!$omp parallel do private(iw, nsea, ka, jws, iws, kw, &
!$omp                     arf_iw, arf_kw, area_dt, ks)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   nsea = itab_w(iw)%nsea
   ka = lpw(iw)

   do jws = 1,nsea
      iws = itab_w(iw)%isea(jws)
      kw  = itab_ws(iws)%kw
      arf_iw  = itab_ws(iws)%arf_iw     ! sea cell area frac of atm IW col
      arf_kw  = itab_ws(iws)%arf_kw     ! sea cell area frac of atm (KW,IW)
                                        !    cell contact with surface
      area_dt = sea%area(iws) * dtlm(1) ! sea cell area * timestep [m^2 * s]

      ks = kw - ka + 1

! Add flux contributions to IW atmospheric column

      ustar      (iw) = ustar      (iw) + arf_iw  * sea%ustar (iws) 
      wtv0       (iw) = wtv0       (iw) + arf_iw  * sea%wthv  (iws)
      sfluxt     (iw) = sfluxt     (iw) + arf_iw  * sea%sfluxt(iws)
      sfluxr     (iw) = sfluxr     (iw) + arf_iw  * sea%sfluxr(iws)
      vkm_sfc (ks,iw) = vkm_sfc (ks,iw) + arf_kw  * sea%vkmsfc(iws)
      sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dt * sea%sfluxt(iws)
      sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dt * sea%sfluxr(iws)
   enddo
enddo
!$omp end parallel do

! Fluxes with LAND cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and LAND cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in LAND cells

!----------------------------------------------------------------------

!dir$ novector
!$omp parallel do private(iw,kw,my_co2,airshv,airthetav,vels,rhos,&
!$omp                     canpress,canexneri,cantheta,canthetav,ufree)
do iwl = 2,mwl

   iw = itab_wl(iwl)%iw
   if (isubdomain == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif
   kw = itab_wl(iwl)%kw

!----------------------------------------------------------------------

! Calculate turbulent fluxes between atmosphere and "land canopy"

   if (allocated(sh_c)) then
      airshv = sh_v(kw,iw) + sh_c(kw,iw)
   else
      airshv = sh_v(kw,iw)
   endif

   airthetav = theta(kw,iw) * (1.0 + eps_virt * sh_v(kw,iw))
   vels      = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   rhos      = rho(kw,iw)

   canpress  = press(kw,iw) + dzt_bot(kw) * rho(kw,iw) * grav
   canexneri = (p00 / canpress) ** rocp
   cantheta  = land%cantemp(iwl) * canexneri
   canthetav = cantheta * (1.0 + eps_virt * land%canshv(iwl))
   ufree     = (grav / airthetav * max(land%wthv(iwl),0.0) * dzt_bot(kw)) ** onethird

! Calculate turbulent fluxes between atmosphere and land canopy

   if (land%ed_flag(iwl) == 0) then

      call stars(dzt_bot      (kw), &
                 land%rough  (iwl), &
                 vels             , &
                 rhos             , &
                 ufree            , &
                 theta     (kw,iw), &
                 airthetav        , &
                 airshv           , &
                 cantheta         , &
                 canthetav        , &
                 land%canshv (iwl), &
                 land%vkmsfc (iwl), &
                 land%sfluxt (iwl), &
                 land%sfluxr (iwl), & 
                 land%ustar  (iwl), &
                 land%rib    (iwl), &
                 land%ggaer  (iwl)  )

      land%wthv(iwl) = ( land%sfluxt(iwl) * (1.0 + eps_virt * sh_v(kw,iw)) &
                       + land%sfluxr(iwl) * eps_virt * theta(kw,iw) ) / rhos

      land%zeta(iwl) = dzt_bot(kw) * grav * vonk * land%wthv(iwl)  &
                     / (airthetav * land%ustar(iwl) ** 3)

    ! If we have co2 or other scalars:
    ! land%sfluxc(iwl) = rhos * land%ggaer(iwl) * (land%canco2(iwl) - atm_co2(kw,iw))
      land%sfluxc(iwl) = 0.

   else

#ifdef USE_ED2

      ! Someday we may track CO2 in OLAM...
      call get_ed2_atm_co2(iwl,my_co2)

      call ed_stars_wrapper( iwl, &
               dzt_bot      (kw), &
               vels             , &
               rhos             , &
               ufree            , &
               theta     (kw,iw), &
               airthetav        , &
               airshv           , &
               my_co2,          , &
               land%vkmsfc (iwl), &
               land%sfluxt (iwl), &
               land%sfluxr (iwl), &
               land%sfluxc (iwl), &
               land%ustar  (iwl), &
               land%rib    (iwl), &
               land%ggaer  (iwl), &
               land%wthv   (iwl), &
               land%zeta   (iwl)  )

#endif

   endif

! Store flux contributions in LAND cell

   land%sxfer_t(iwl) = dtlm(1) * land%sfluxt(iwl) * tair(kw,iw) / theta(kw,iw)
   land%sxfer_r(iwl) = dtlm(1) * land%sfluxr(iwl)
   land%sxfer_c(iwl) = dtlm(1) * land%sfluxc(iwl)

enddo
!$omp end parallel do

!$omp parallel do private(iw, nland, ka, jwl, iwl, kw, &
!$omp                     arf_iw, arf_kw, area_dt, ks) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
   nland = itab_w(iw)%nland
   ka = lpw(iw)

   do jwl = 1,nland
      iwl = itab_w(iw)%iland(jwl)
      kw      = itab_wl(iwl)%kw
      arf_iw  = itab_wl(iwl)%arf_iw      ! land cell area frac of atm IW col
      arf_kw  = itab_wl(iwl)%arf_kw      ! land cell area frac of atm (KW,IW)
                                         !    cell contact with surface
      area_dt = land%area(iwl) * dtlm(1) ! land cell area * timestep [m^2 * s]

      ks = kw - ka + 1

! Add flux contributions to IW atmospheric column

      ustar      (iw) = ustar      (iw) + arf_iw  * land%ustar (iwl) 
      wtv0       (iw) = wtv0       (iw) + arf_iw  * land%wthv  (iwl)
      sfluxt     (iw) = sfluxt     (iw) + arf_iw  * land%sfluxt(iwl)
      sfluxr     (iw) = sfluxr     (iw) + arf_iw  * land%sfluxr(iwl)
      vkm_sfc (ks,iw) = vkm_sfc (ks,iw) + arf_kw  * land%vkmsfc(iwl)
      sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dt * land%sfluxt(iwl)
      sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dt * land%sfluxr(iwl)

   enddo
enddo
!$omp end parallel do

! Compute some derived surface quantities

!$omp parallel do private(iw,ka)
do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

   ka = lpw(iw)

   if (wtv0(iw) > 0.0) then
      wstar(iw) = (grav * pblh(iw) * wtv0(iw) / theta(ka,iw)) ** onethird
   else
      wstar(iw) = 0.0
   endif

enddo
!$omp end parallel do

end subroutine surface_turb_flux

!===============================================================================

subroutine stars( zts, rough, vels, rhos, ufree,  &
                  air_theta, air_thetav, air_shv, &
                  can_theta, can_thetav, can_shv, &
                  vkmsfc, sfluxt, sfluxr, ustar,  &
                  ri, ggaero                      )

! Subroutine stars computes surface heat and vapor fluxes and momentum drag
! coefficient from Louis (1981) equations

use consts_coms, only: vonk, grav

implicit none

! Input variables

real, intent(in) :: zts        ! height above surface of {vels, ths, shv} [m]
real, intent(in) :: rough      ! surface roughness height [m]
real, intent(in) :: vels       ! atmos near-surface wind speed [m/s]
real, intent(in) :: rhos       ! atmos near-surface density [kg/m^3]
real, intent(in) :: ufree      ! surface layer free-convective velocity [m/s]
real, intent(in) :: air_theta  ! atmos near-surface pot. temp [K]
real, intent(in) :: air_thetav ! atmos near-surface virt. pot. temp [K]
real, intent(in) :: air_shv    ! atmos near-surface vapor spec hum [kg_vap/m^3]
real, intent(in) :: can_theta  ! canopy air pot. temp [K]
real, intent(in) :: can_thetav ! canopy air virt. pot. temp [K]
real, intent(in) :: can_shv    ! canopy air vapor spec hum [kg_vap/m^3]

! Output variables

real, intent(out) :: vkmsfc    ! surface drag coefficient for this flux cell
real, intent(out) :: sfluxt    ! surface sensible heat flux for this flux cell
real, intent(out) :: sfluxr    ! surface vapor flux for this flux cell
real, intent(out) :: ustar     ! surface friction velocity for this flux cell
real, intent(out) :: ri        ! bulk richardson number, eq. 3.45 in Garratt
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
rstar = c3 * (air_shv   - can_shv)

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

use mem_cuparm,  only: conprr
use mem_leaf,    only: land, itab_wl
use mem_sea,     only: sea, itab_ws
use mem_ijtabs,  only: istp, mrl_begl, itabg_w
use consts_coms, only: cliq, alli, t00
use mem_basic,   only: tair, rho, sh_v
use leaf_coms,   only: mwl, mrl_leaf
use sea_coms,    only: mws
use misc_coms,   only: io6, isubdomain, dtlm
use ed_misc_coms,only: ed2_active

implicit none

integer :: mrl
integer :: iw
integer :: iwl, iws
integer :: kw

real :: airtempc
real :: tempc

real, external :: rhovsl

! Subroutine to transfer atmospheric cumulus parameterization 
! precipitation FLUX to surface cells

! Set land fluxes to be done for CUPARM
!    1. for mrl = 0, no fluxes 
!    2. For mrl > mrl_leaf, no fluxes
!    3. For mrl <= mrl_leaf, do fluxes at beginning of long timestep 
!                            for all points in mrl = 1

mrl = mrl_begl(istp)
if (mrl > 0) then

! Transfer precipitation FLUX to leaf land cells

   !$omp parallel do private(iw,kw,airtempc,tempc)
   do iwl = 2, mwl

      iw = itab_wl(iwl)%iw
      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      kw = itab_wl(iwl)%kw

! Compute air temperature in C

      airtempc = tair(kw,iw) - t00

! Estimate wet bulb temp using computation from subroutine each_column in micphys
! Assume that convective precip reaches surface at this wet bulb temp

      tempc = airtempc - min(25.,  &
           700. * (rhovsl(airtempc) / real(rho(kw,iw)) - sh_v(kw,iw)))

      land%pcpg (iwl) = land%pcpg (iwl) + dtlm(1) * conprr(iw)
      land%qpcpg(iwl) = land%qpcpg(iwl) + dtlm(1) * conprr(iw) * (cliq * tempc + alli)
      land%dpcpg(iwl) = land%dpcpg(iwl) + dtlm(1) * conprr(iw) * .001

   enddo
   !$omp end parallel do

#ifdef USE_ED2
   if (ed2_active == 1) then
      call copy_cuparm_to_ed()
   endif
#endif

! Transfer precipitation FLUX to sea cells

   !$omp parallel do private(iw,kw,airtempc,tempc)
   do iws = 2, mws

      iw = itab_ws(iws)%iw
      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      kw = itab_ws(iws)%kw

! Compute air temperature in C

      airtempc = tair(kw,iw) - t00

! Estimate wet bulb temp using computation from subroutine each_column in micphys
! Assume that convective precip reaches surface at this wet bulb temp

      tempc = airtempc - min(25.,  &
           700. * (rhovsl(airtempc) / real(rho(kw,iw)) - sh_v(kw,iw)))

      sea%pcpg (iws) = sea%pcpg (iws) + dtlm(1) * conprr(iw)
      sea%qpcpg(iws) = sea%qpcpg(iws) + dtlm(1) * conprr(iw) * (cliq * tempc + alli)
      sea%dpcpg(iws) = sea%dpcpg(iws) + dtlm(1) * conprr(iw) * .001

   enddo
   !$omp end parallel do

endif

end subroutine surface_cuparm_flux

