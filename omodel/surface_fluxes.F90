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
use misc_coms,   only: io6, iparallel, isubdomain, dtlm, mdomain
use mem_grid,    only: mza, mwa, lsw, lpw, zt, zm, arw
use mem_sea,     only: sea, itab_ws
use mem_leaf,    only: land, itab_wl
use mem_turb,    only: vkm_sfc, sfluxt, sfluxr, sxfer_tk, sxfer_rk, &
                       ustar, wstar, wtv0, pblh
use mem_basic,   only: press, rho, theta, tair, sh_v, vxe, vye, vze
use mem_micro,   only: sh_c
use consts_coms, only: grav, p00i, rocp, cpi

implicit none

integer, intent(in) :: mrl

integer :: j,jv,iw,iv,iws,iwl
integer :: ks,kw,ka,k
integer :: nsea,jws,nland,jwl


real :: area_dt
real :: arf_kw
real :: arf_iw
real :: exneri
real :: vkmsfc, vkmsfcs, vkmsfci
real :: vels
real :: my_co2

real :: shflx, sh_vc

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

         exneri = 1. / ((p00i * press(kw,iw)) ** rocp)

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
!$omp parallel do private(iw,ka)
do j = 1,jtab_w(jtw_wstn)%jend(mrl); iw = jtab_w(jtw_wstn)%iw(j)
!----------------------------------------------------------------------

   ka = lpw(iw)

   vkm_sfc(:,iw) = 0.
   ustar(iw)  = 0.
   sfluxt(iw) = 0.
   sfluxr(iw) = 0.

enddo
!$omp end parallel do

! Fluxes with SEA cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and SEA cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in SEA cells

!----------------------------------------------------------------------

!$omp parallel do private(iw,kw)
do iws = 2,mws

   iw = itab_ws(iws)%iw
   if (isubdomain == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif
   kw = itab_ws(iws)%kw
!----------------------------------------------------------------------

! Calculate turbulent fluxes between atmosphere and "water canopy"

   sea%rhos   (iws) = real(rho(kw,iw))
   sea%vels   (iws) = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   sea%prss   (iws) = real(press(kw,iw))
   sea%airtemp(iws) = tair(kw,iw)

   if (allocated(sh_c)) then
      sea%airshv(iws) = sh_v(kw,iw) + sh_c(kw,iw)
   else
      sea%airshv(iws) = sh_v(kw,iw)
   endif

! Sea (open water) component always computed

   call stars(zt(kw)-zm(kw-1),      &
              sea%sea_rough  (iws), &
              sea%vels       (iws), &
              sea%rhos       (iws), &
              sea%airtemp    (iws), &
              sea%airshv     (iws), &
              sea%sea_cantemp(iws), &
              sea%sea_canshv (iws), &
              sea%sea_vkmsfc (iws), &
              sea%sea_sfluxt (iws), &
              sea%sea_sfluxr (iws), & 
              sea%sea_ustar  (iws), &
              sea%sea_ggaer  (iws)  )

    sea%sfluxc(iws) = 0.

! Include fractional seaice component if seaice layers exist

   if (sea%nlev_seaice(iws) > 0) then

      call stars(zt(kw)-zm(kw-1),     &
                sea%ice_rough  (iws), &
                sea%vels       (iws), &
                sea%rhos       (iws), &
                sea%airtemp    (iws), &
                sea%airshv     (iws), &
                sea%ice_cantemp(iws), &
                sea%ice_canshv (iws), &
                sea%ice_vkmsfc (iws), &
                sea%ice_sfluxt (iws), &
                sea%ice_sfluxr (iws), &
                sea%ice_ustar  (iws), &
                sea%ice_ggaer  (iws)  )

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

   else

      sea%vkmsfc(iws) = sea%sea_vkmsfc(iws)
      sea%ustar (iws) = sea%sea_ustar (iws)
      sea%ggaer (iws) = sea%sea_ggaer (iws)
      sea%sfluxt(iws) = sea%sea_sfluxt(iws)
      sea%sfluxr(iws) = sea%sea_sfluxr(iws)

      sea%ice_vkmsfc(iws) = 0.0
      sea%ice_ustar (iws) = 0.0
      sea%ice_ggaer (iws) = 0.0
      sea%ice_sfluxt(iws) = 0.0
      sea%ice_sfluxr(iws) = 0.0
      
   endif

! Store atmospheric properties and flux contributions in SEA cell

   sea%sxfer_t    (iws) = dtlm(1) * sea%sfluxt    (iws)
   sea%sea_sxfer_t(iws) = dtlm(1) * sea%sea_sfluxt(iws)
   sea%ice_sxfer_t(iws) = dtlm(1) * sea%ice_sfluxt(iws)

   sea%sxfer_r    (iws) = dtlm(1) * sea%sfluxr    (iws)
   sea%sea_sxfer_r(iws) = dtlm(1) * sea%sea_sfluxr(iws)
   sea%ice_sxfer_r(iws) = dtlm(1) * sea%ice_sfluxr(iws)

   sea%sxfer_c    (iws) = dtlm(1) * sea%sfluxc    (iws)

! Reset surface precipitation flux arrays

   sea%pcpg (iws) = 0.
   sea%qpcpg(iws) = 0.
   sea%dpcpg(iws) = 0.

enddo
!$omp end parallel do

!$omp parallel do private(iw, nsea, ka, jws, iws, kw, &
!$omp                     arf_iw, arf_kw, area_dt, ks, exneri)
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

      exneri = theta(kw,iw) / tair(kw,iw)

      ustar      (iw) = ustar      (iw) + arf_iw  * sea%ustar (iws) 
       sfluxt    (iw) = sfluxt     (iw) + arf_iw  * sea%sfluxt(iws) * exneri
       sfluxr    (iw) = sfluxr     (iw) + arf_iw  * sea%sfluxr(iws)
      vkm_sfc (ks,iw) = vkm_sfc (ks,iw) + arf_kw  * sea%vkmsfc(iws)
      sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dt * sea%sfluxt(iws) * exneri
      sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dt * sea%sfluxr(iws)
   enddo
enddo
!$omp end parallel do

! Fluxes with LAND cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and LAND cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in LAND cells

!----------------------------------------------------------------------

!$omp parallel do private(iw,kw,my_co2)
do iwl = 2,mwl

   iw = itab_wl(iwl)%iw
   if (isubdomain == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif
   kw = itab_wl(iwl)%kw

!----------------------------------------------------------------------

! Calculate turbulent fluxes between atmosphere and "land canopy"

   land%prss   (iwl) = real(press(kw,iw))
   land%vels   (iwl) = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   land%rhos   (iwl) = real(rho(kw,iw))
   land%airtemp(iwl) = tair(kw,iw)

   if (allocated(sh_c)) then
      land%airshv(iwl) = sh_v(kw,iw) + sh_c(kw,iw)
   else
      land%airshv(iwl) = sh_v(kw,iw)
   endif

! Calculate turbulent fluxes between atmosphere and land canopy

   if (land%ed_flag(iwl) == 0) then

      call stars(zt(kw)-zm(kw-1), &
               land%rough  (iwl), &
               land%vels   (iwl), &
               land%rhos   (iwl), &
               land%airtemp(iwl), &
               land%airshv (iwl), &
               land%cantemp(iwl), &
               land%canshv (iwl), &
               land%vkmsfc (iwl), &
               land%sfluxt (iwl), &
               land%sfluxr (iwl), & 
               land%ustar  (iwl), &
               land%ggaer  (iwl)  )

      land%sfluxc (iwl) = 0.
      land%ed_zeta(iwl) = 0.
      land%ed_rib (iwl) = 0.

   else

#ifdef USE_ED2

      ! Someday we may track CO2 in OLAM...
      call get_ed2_atm_co2(iwl,my_co2)

      call ed_stars_wrapper(iwl, zt(kw)-zm(kw-1), &
           land%vels   (iwl), &
           land%rhos   (iwl), &
           land%airtemp(iwl), &
           land%airshv (iwl), &
           my_co2,            &
           land%vkmsfc (iwl), &
           land%sfluxt (iwl), &
           land%sfluxr (iwl), & 
           land%sfluxc (iwl), &
           land%ustar  (iwl), &
           land%ed_zeta(iwl), &
           land%ed_rib (iwl), &
           land%ggaer  (iwl)  )

#endif

   endif

! Store flux contributions in LAND cell

   land%sxfer_t(iwl) = dtlm(1) * land%sfluxt(iwl)
   land%sxfer_r(iwl) = dtlm(1) * land%sfluxr(iwl)
   land%sxfer_c(iwl) = dtlm(1) * land%sfluxc(iwl)

! Reset surface precipitation flux arrays

   land%pcpg (iwl) = 0.
   land%qpcpg(iwl) = 0.
   land%dpcpg(iwl) = 0.

enddo
!$omp end parallel do

!$omp parallel do private(iw, nland, ka, jwl, iwl, kw, &
!$omp                     arf_iw, arf_kw, area_dt, ks, exneri) 
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

      exneri = theta(kw,iw) / tair(kw,iw)

      ustar      (iw) = ustar      (iw) + arf_iw  * land%ustar (iwl) 
      sfluxt     (iw) = sfluxt     (iw) + arf_iw  * land%sfluxt(iwl) * exneri
      sfluxr     (iw) = sfluxr     (iw) + arf_iw  * land%sfluxr(iwl)
      vkm_sfc (ks,iw) = vkm_sfc (ks,iw) + arf_kw  * land%vkmsfc(iwl)
      sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dt * land%sfluxt(iwl) * exneri
      sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dt * land%sfluxr(iwl)

   enddo
enddo
!$omp end parallel do

! Compute some derived surface quantities

!$omp parallel do private(iw,ka)
do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

   ka = lpw(iw)
   
   wtv0(iw) = sfluxt(iw) * (1. + .61 * sh_v(ka,iw)) &
            + sfluxr(iw) * .61 * theta(ka,iw)

   if (wtv0(iw) > 0.0) then
      wstar(iw) = (grav * pblh(iw) * wtv0(iw) / theta(ka,iw)) ** 0.33333333
   else
      wstar(iw) = 0.0
   endif

enddo
!$omp end parallel do

return
end subroutine surface_turb_flux

!===============================================================================

subroutine stars(zts, rough, vels, rhos, airtemp, sh_vs, cantemp, &
                 canshv, vkmsfc, sfluxt, sfluxr, ustar0, ggaero  )

! Subroutine stars computes surface heat and vapor fluxes and momentum drag
! coefficient from Louis (1981) equations

use consts_coms, only: vonk, grav
use misc_coms,   only: io6

implicit none

! Input variables

real, intent(in) :: zts       ! height above surface of {vels, ths, sh_vs} [m]
real, intent(in) :: rough     ! surface roughness height [m]
real, intent(in) :: vels      ! atmos near-surface wind speed [m/s]
real, intent(in) :: airtemp   ! atmos near-surface temp [K]
real, intent(in) :: sh_vs     ! atmos near-surface vapor spec hum [kg_vap/m^3]
real, intent(in) :: cantemp   ! canopy air temp [K]
real, intent(in) :: canshv   ! canopy air vapor spec hum [kg_vap/m^3]

real, intent(in) :: rhos  ! atmos near-surface density [kg/m^3]

! Output variables

real, intent(out) :: vkmsfc   ! surface drag coefficient for this flux cell
real, intent(out) :: sfluxt   ! surface sensible heat flux for this flux cell
real, intent(out) :: sfluxr   ! surface vapor flux for this flux cell
real, intent(out) :: ustar0   ! surface friction velocity for this flux cell
real, intent(out) :: ggaero   ! bare ground conductance m/s

! Local parameters

real, parameter :: b = 5.
real, parameter :: csm = 7.5
real, parameter :: csh = 5.
real, parameter :: d = 5.
real, parameter :: ustmin = .1  ! lower bound on ustar (friction velocity)
real, parameter :: ubmin = .25  ! lower bound on wind speed (should use 1.0 for
                                !   convec case and 0.1 for stable case)
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
real :: ri     ! bulk richardson number, eq. 3.45 in Garratt
real :: tstar  ! 
real :: rstar  !
real :: vtscr  ! ustar0 times density

! Routine to compute Louis (1981) surface layer parameterization.

vels0 = max(vels,ubmin)

a2 = (vonk / log(zts / rough)) ** 2
c1 = a2 * vels0
ri = grav * zts * (airtemp - cantemp)  &
   / (.5 * (airtemp + cantemp) * vels0 * vels0)

if (airtemp - cantemp > 0.) then   ! STABLE CASE

   fm = 1. / (1. + (2. * b * ri / sqrt(1. + d * ri)))
   fh = 1. / (1. + (3. * b * ri * sqrt(1. + d * ri)))

else                            ! UNSTABLE CASE

   c2 = b * a2 * sqrt(zts / rough * (abs(ri)))
   cm = csm * c2
   ch = csh * c2
   fm = (1. - 2. * b * ri / (1. + 2. * cm))
   fh = (1. - 3. * b * ri / (1. + 3. * ch))
   
endif

ustar0 = max(ustmin,sqrt(c1 * vels0 * fm))
c3 = c1 * fh / ustar0
tstar = c3 * (airtemp - cantemp)
rstar = c3 * (sh_vs - canshv)

vtscr = ustar0 * rhos

vkmsfc = vtscr * ustar0 * zts / vels0
sfluxt = - vtscr * tstar
sfluxr = - vtscr * rstar

! Store the aerodynamic conductance between the surface canopy and
! the lowest model level

ggaero  = c3 * ustar0

return
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

!===============================================================================

subroutine surface_precip_flux()

use mem_micro,   only: pcpgr, qpcpgr, dpcpgr
use mem_leaf,    only: land, itab_wl
use mem_sea,     only: sea, itab_ws
use mem_ijtabs,  only: istp, mrl_endl, itabg_w
use leaf_coms,   only: mwl, mrl_leaf
use sea_coms,    only: mws
use misc_coms,   only: io6, isubdomain, dtlm
use ed_misc_coms,only: ed2_active

implicit none

integer :: iw
integer :: iwl, iws
integer :: mrl

! Subroutine to transfer atmospheric microphysics parameterization 
! precipitation flux to surface cells.

! Land fluxes to be done for PRECIP 
!    1. for mrl = 0, no fluxes 
!    2. For mrl > mrl_leaf, do fluxes at end of long timestep at given mrl
!    3. For mrl <= mrl_leaf, do fluxes at end of long timestep
!                            for all points in mrl = 1

mrl = mrl_endl(istp)
if (mrl > 0) then

! Transfer precipitation FLUX to leaf land cells

   !$omp parallel do private(iw)
   do iwl = 2, mwl

      iw = itab_wl(iwl)%iw
      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      land%pcpg (iwl) = land%pcpg (iwl) + dtlm(1) * pcpgr (iw)
      land%qpcpg(iwl) = land%qpcpg(iwl) + dtlm(1) * qpcpgr(iw)
      land%dpcpg(iwl) = land%dpcpg(iwl) + dtlm(1) * dpcpgr(iw)

   enddo
   !$omp end parallel do

#ifdef USE_ED2
   if (ed2_active == 1) then
      call copy_micro_to_ed()
   endif
#endif

! Transfer precipitation FLUX to sea cells

   !$omp parallel do private(iw)
   do iws = 2, mws

      iw = itab_ws(iws)%iw
      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      sea%pcpg (iws) = sea%pcpg (iws) + dtlm(1) * pcpgr (iw)
      sea%qpcpg(iws) = sea%qpcpg(iws) + dtlm(1) * qpcpgr(iw)
      sea%dpcpg(iws) = sea%dpcpg(iws) + dtlm(1) * dpcpgr(iw)

   enddo
   !$omp end parallel do

endif

! pcpgr, qpcpgr, dpcpgr have now been "transferred" to leaf cells.
! No need to zero pcpgr, qpcpgr, dpcpgr since microphysics will replace them.

end subroutine surface_precip_flux
