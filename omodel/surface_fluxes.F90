!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

! TURB, CUPARM, and PRECIP flux scheduling:

! 1. Cuparm and micphys give precip RATES
! 2. Stars gives heat and vapor flux RATES
! 3. Fluxes computed as often as the MORE FREQUENT of atm and leaf
! 4. Fluxes converted to AMOUNTS TRANSFERED based on rate * dtf 
!    (DTF SAME FOR MICRO & TURB?)
! 5. Each time leaf runs, it uses at once all it has gotten and zeroes xfer arrays
! 6. For MRL with dtlm < dt_leaf, only do fluxes for MRL atm points
! 7. For MRL with dtlm >= dt_leaf, do fluxes for MRL = 1 atm points

! RADIATIVE fluxes only:

! 1. Radiative fluxes transfer RATES only, with no timestep information

!----------------------------------------------------------------------------
! AS OF JULY 2011, SUBROUTINE SURFACE_TURB_FLUXP IS HARDWIRED FOR 
! DT_LEAF = DTLONG(1) AND MRL_LEAF = 1, WHICH FOR A LONG TIME HAS BEEN 
! ASSIGNED IN SUBROUTINE MODSCHED AND USED SUCCESSFULLY.
! IF DT_LEAF EVER GETS LARGE ENOUGH TO CAUSE INSTABILITY, THE INSTABILITY
! SHOULD BE CONTROLLED BY INCREASING CAPACITANCE OF THE LEAF COMPONENTS, 
! NOT BY REDUCING TIMESTEP; I.E., WE WANT LEAF TO BE CAPABLE OF RUNNING ON
! AS LONG A TIMESTEP AS THE ATMOSPHERIC MODEL.
!----------------------------------------------------------------------------

subroutine surface_turb_flux(mrl)

use leaf_coms,   only: mwl, dt_leaf, isfcl
use sea_coms,    only: mws, dt_sea
use mem_ijtabs,  only: itab_w, itabg_w, jtab_w, jtw_prog, jtw_wstn
use misc_coms,   only: io6, iparallel
use mem_para,    only: myrank
use mem_grid,    only: mza, mwa, lsw, lpw, zt, zm
use mem_sea,     only: sea, itabg_ws, itab_ws
use mem_leaf,    only: land, itabg_wl, itab_wl
use mem_sflux,   only: seaflux, landflux, jseaflux, jlandflux
use mem_turb,    only: vkm_sfc, sflux_t, sflux_r, sxfer_tk, sxfer_rk, &
                       ustar, wstar, wtv0, pblh
use mem_basic,   only: press, rho, theta, tair, sh_v, vxe, vye, vze
use consts_coms, only: grav

implicit none

integer, intent(in) :: mrl

integer :: j,jv,iw,iv,iws,iwl
integer :: ks,kw,ka,k
integer :: isf,ilf

real :: area_dtf
real :: arf_kw
real :: arf_atm
real :: arf_sea
real :: arf_land
real :: arf_sea_dtf
real :: arf_land_dtf
real :: exneri
real :: vkmsfc, vkmsfcs, vkmsfci
real :: vels
real :: my_co2, ed_zeta, ed_rib

if (isfcl == 0) then

! ISFCL = 0 is the no-LEAF option:  Set all surface fluxes (to zero by default).
! These flux values can be altered if desired.

   call psub()
!-----------------------------------------------------------------------------
   do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!-----------------------------------------------------------------------------
   call qsub('W',iw)

      sflux_t(iw) = 0.
      sflux_r(iw) = 0.
      ustar  (iw) = .1  ! Minimum value

      wstar  (iw) = 0.
      wtv0   (iw) = 0.

      do ks = 1,lsw(iw)
         vkm_sfc (ks,iw) = 0.
         sxfer_tk(ks,iw) = 0.
         sxfer_rk(ks,iw) = 0.
      enddo

!      albedt (iw) = albedo
!      rlongup(iw) = stefan * 288.15 ** 4  ! std msl temp; make decision to not
                                           ! run radiation if not running leaf?
   enddo
   call rsub('W_noleaf',13)

   return
endif

! ISFCL = 1 is the LEAF option...

! Reset to zero the atm values of VKM_SFC, USTAR, SFLUX_T, and SFLUX_R
! for same mrl's that fluxes will be evaluated for this time.  
! THUS, VKM_SFC, USTAR, SFLUX_T, AND SFLUX_R ARE ONLY SUMMED OVER
! SPACE, BUT NOT OVER TIME. (ON THE OTHER HAND, SXFER_TK AND SXFER_RK ARE SUMMED
! OVER BOTH SPACE AND TIME; THEY ARE RESET TO ZERO IN THILTEND_LONG AND
! SCALAR_TRANSPORT AFTER THEY ARE TRANSFERRED TO THE ATMOSPHERE.)

! Set sea and land fluxes to be done for SURFACE SIMILARITY: 
!    Do fluxes at beginning of long timestep at given mrl

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_wstn)%jend(mrl); iw = jtab_w(jtw_wstn)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   ka = lpw(iw)

   vkm_sfc(:,iw) = 0.
   ustar(iw)   = 0.
   sflux_t(iw) = 0.
   sflux_r(iw) = 0.

enddo
call rsub('W',20)

! Fluxes with SEA cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and SEA cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in SEAFLUX cells

call psub()
!----------------------------------------------------------------------
do j = 1,jseaflux(1)%jend(mrl)
   isf = jseaflux(1)%iseaflux(j)
   iw  = seaflux(isf)%iw         ! global index
   iws = seaflux(isf)%iwls       ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iw  = itabg_w (iw )%iw_myrank
      iws = itabg_ws(iws)%iws_myrank
   endif

   kw          = seaflux(isf)%kw
   arf_kw      = seaflux(isf)%arf_kw
   arf_atm     = seaflux(isf)%arf_atm    ! area ratio of flux cell to atm cell
   arf_sea     = seaflux(isf)%arf_sfc    ! area ratio of flux cell to sea cell
   arf_sea_dtf = seaflux(isf)%arf_sfc  & ! area ratio of flux cell to sea cell
               * seaflux(isf)%dtf        !   timestep of flux cell [m^2 * s]
   area_dtf    = seaflux(isf)%area     & ! area of flux cell
               * seaflux(isf)%dtf        !   timestep of flux cell [m^2 * s]

   ka = lpw(iw)
   ks = kw - ka + 1
!----------------------------------------------------------------------
   call qsub('SF1',iw)

! Calculate turbulent fluxes between atmosphere and "water canopy"

   vels   = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   exneri = theta(kw,iw) / tair(kw,iw)

! Sea (open water) component always computed

   call stars(zt(kw)-zm(kw-1),         &
              sea%sea_rough(iws),      &
              vels,                    &
              rho      (kw,iw),        &
              tair     (kw,iw),        &
              sh_v     (kw,iw),        &
              sea%seacan_temp(iws),    &
              sea%seacan_shv (iws),    &
              vkmsfcs,                 &
              seaflux(isf)%sea_sfluxt, &
              seaflux(isf)%sea_sfluxr, & 
              seaflux(isf)%sea_ustar0, &
              seaflux(isf)%sea_ggaer0  )

! Include fractional seaice component if seaice layers exist

   if (sea%nlev_seaice(iws) > 0) then

      call stars(zt(kw)-zm(kw-1),        &
                sea%ice_rough(iws),      &
                vels,                    &
                rho      (kw,iw),        &
                tair     (kw,iw),        &
                sh_v     (kw,iw),        &
                sea%icecan_temp(iws),    &
                sea%icecan_shv (iws),    &
                vkmsfci,                 &
                seaflux(isf)%ice_sfluxt, &
                seaflux(isf)%ice_sfluxr, &
                seaflux(isf)%ice_ustar0, &
                seaflux(isf)%ice_ggaer0  )

      vkmsfc              = (1.0 - sea%seaicec(iws)) * vkmsfcs + &
                                   sea%seaicec(iws)  * vkmsfci

      seaflux(isf)%ustar0 = (1.0 - sea%seaicec(iws)) * seaflux(isf)%sea_ustar0 + &
                                   sea%seaicec(iws)  * seaflux(isf)%ice_ustar0

      seaflux(isf)%ggaer0 = (1.0 - sea%seaicec(iws)) * seaflux(isf)%sea_ggaer0 + &
                                   sea%seaicec(iws)  * seaflux(isf)%ice_ggaer0

      seaflux(isf)%sfluxt = (1.0 - sea%seaicec(iws)) * seaflux(isf)%sea_sfluxt + &
                                   sea%seaicec(iws)  * seaflux(isf)%ice_sfluxt

      seaflux(isf)%sfluxr = (1.0 - sea%seaicec(iws)) * seaflux(isf)%sea_sfluxr + &
                                   sea%seaicec(iws)  * seaflux(isf)%ice_sfluxr

   else
      
      vkmsfc              = vkmsfcs
      seaflux(isf)%ustar0 = seaflux(isf)%sea_ustar0
      seaflux(isf)%ggaer0 = seaflux(isf)%sea_ggaer0
      seaflux(isf)%sfluxt = seaflux(isf)%sea_sfluxt
      seaflux(isf)%sfluxr = seaflux(isf)%sea_sfluxr

      vkmsfci                 = 0.0
      seaflux(isf)%ice_ustar0 = 0.0
      seaflux(isf)%ice_ggaer0 = 0.0
      seaflux(isf)%ice_sfluxt = 0.0
      seaflux(isf)%ice_sfluxr = 0.0
      
   endif

! Add flux contributions to IW atmospheric column

   vkm_sfc (ks,iw) = vkm_sfc(ks,iw)  + arf_kw   * vkmsfc
   ustar      (iw) = ustar  (iw)     + arf_atm  * seaflux(isf)%ustar0
   sflux_t    (iw) = sflux_t(iw)     + arf_atm  * seaflux(isf)%sfluxt * exneri
   sflux_r    (iw) = sflux_r(iw)     + arf_atm  * seaflux(isf)%sfluxr
   sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dtf * seaflux(isf)%sfluxt * exneri
   sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dtf * seaflux(isf)%sfluxr

! Store flux and ATM density contributions to SEA cell in SEAFLUX cell

   seaflux(isf)%rhos        = arf_sea     * rho(kw,iw)
   seaflux(isf)%sxfer_t     = arf_sea_dtf * seaflux(isf)%sfluxt
   seaflux(isf)%sxfer_r     = arf_sea_dtf * seaflux(isf)%sfluxr
   seaflux(isf)%sxfer_c     = 0.
   seaflux(isf)%ustar       = arf_sea     * seaflux(isf)%ustar0
   seaflux(isf)%ggaer       = arf_sea     * seaflux(isf)%ggaer0

   seaflux(isf)%sea_sxfer_t = arf_sea_dtf * seaflux(isf)%sea_sfluxt
   seaflux(isf)%sea_sxfer_r = arf_sea_dtf * seaflux(isf)%sea_sfluxr
   seaflux(isf)%sea_ustar   = arf_sea     * seaflux(isf)%sea_ustar0
   seaflux(isf)%sea_ggaer   = arf_sea     * seaflux(isf)%sea_ggaer0

   seaflux(isf)%ice_sxfer_t = arf_sea_dtf * seaflux(isf)%ice_sfluxt
   seaflux(isf)%ice_sxfer_r = arf_sea_dtf * seaflux(isf)%ice_sfluxr
   seaflux(isf)%ice_ustar   = arf_sea     * seaflux(isf)%ice_ustar0
   seaflux(isf)%ice_ggaer   = arf_sea     * seaflux(isf)%ice_ggaer0

enddo

! Do parallel send of TURBULENT fluxes and ATM properties for SEA cells

if (iparallel == 1) call mpi_send_wsf('T',mrl)

call rsub('JSEAFLUX',1)

! Fluxes with LAND cells

! 1. Evaluate turbulent surface fluxes between atmosphere columns and LAND cells
! 2. Add flux contribution from each flux cell to atmosphere columns
! 3. Store fluxes and atmospheric properties in LANDFLUX cells

call psub()
!----------------------------------------------------------------------
do j = 1,jlandflux(1)%jend(mrl)
   ilf = jlandflux(1)%ilandflux(j)
   iw  = landflux(ilf)%iw         ! global index
   iwl = landflux(ilf)%iwls       ! global index

! If run is parallel, get local rank indices

    if (iparallel == 1) then
      iw  = itabg_w (iw )%iw_myrank
      iwl = itabg_wl(iwl)%iwl_myrank
   endif

   kw           = landflux(ilf)%kw
   arf_kw       = landflux(ilf)%arf_kw
   arf_atm      = landflux(ilf)%arf_atm    ! flux cell to atm cell area ratio
   arf_land     = landflux(ilf)%arf_sfc    ! flux cell to land cell area ratio
   arf_land_dtf = landflux(ilf)%arf_sfc  & ! flux cell to land cell area ratio
                * landflux(ilf)%dtf        !   timestep of flux cell [m^2 * s]
   area_dtf     = landflux(ilf)%area     & ! flux cell area 
                * landflux(ilf)%dtf        !   timestep of flux cell [m^2 * s]

   ka = lpw(iw)
   ks = kw - ka + 1

!----------------------------------------------------------------------
   call qsub('LF1',iw)

   vels   = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
   exneri = theta(kw,iw) / tair(kw,iw)

! Store atmospheric properties in landflux cell

   landflux(ilf)%vels = arf_land * vels
   landflux(ilf)%prss = arf_land * press(kw,iw)
   landflux(ilf)%rhos = arf_land * rho(kw,iw)

! Calculate turbulent fluxes between atmosphere and land canopy

   if (land%ed_flag(iwl) == 0) then

      call stars(zt(kw)-zm(kw-1),      &
                 land%rough    (iwl),  &
                 vels,                 &
                 rho         (kw,iw),  &
                 tair        (kw,iw),  &
                 sh_v        (kw,iw),  &
                 land%can_temp (iwl),  &
                 land%can_shv  (iwl),  &
                 vkmsfc,               &
                 landflux(ilf)%sfluxt, &
                 landflux(ilf)%sfluxr, &
                 landflux(ilf)%ustar0, &
                 landflux(ilf)%ggaer0  )

      landflux(ilf)%sfluxc = 0.
      ed_zeta = 0.
      ed_rib  = 0.

   else

#ifdef USE_ED2

      ! Someday we may track CO2 in OLAM...
      call get_ed2_atm_co2(iwl,my_co2)

      call ed_stars_wrapper(iwl, zt(kw)-zm(kw-1),      &
           vels,                 &
           rho         (kw,iw),  &
           tair        (kw,iw),  &
           sh_v        (kw,iw),  &
           my_co2,               &
           vkmsfc,               &
           landflux(ilf)%sfluxt, &
           landflux(ilf)%sfluxr, &
           landflux(ilf)%sfluxc, &
           landflux(ilf)%ustar0, &
           ed_zeta, ed_rib,      &
           landflux(ilf)%ggaer0  )
#endif

   endif

! Add flux contributions to IW atmospheric column

   vkm_sfc(ks,iw)  = vkm_sfc(ks,iw)  + arf_kw   * vkmsfc
   ustar  (iw)     = ustar  (iw)     + arf_atm  * landflux(ilf)%ustar0
   sflux_t(iw)     = sflux_t(iw)     + arf_atm  * landflux(ilf)%sfluxt * exneri
   sflux_r(iw)     = sflux_r(iw)     + arf_atm  * landflux(ilf)%sfluxr
   sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dtf * landflux(ilf)%sfluxt * exneri
   sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dtf * landflux(ilf)%sfluxr

! Store flux contributions to land cell in landflux cell

   landflux(ilf)%sxfer_t = arf_land_dtf * landflux(ilf)%sfluxt
   landflux(ilf)%sxfer_r = arf_land_dtf * landflux(ilf)%sfluxr
   landflux(ilf)%sxfer_c = arf_land_dtf * landflux(ilf)%sfluxc
   landflux(ilf)%ustar   = arf_land     * landflux(ilf)%ustar0
   landflux(ilf)%ggaer   = arf_land     * landflux(ilf)%ggaer0
   landflux(ilf)%ed_zeta = arf_land     * ed_zeta
   landflux(ilf)%ed_rib  = arf_land     * ed_rib

enddo

! Do parallel send of TURBULENT fluxes and ATM properties

if (iparallel == 1) call mpi_send_wlf('T',mrl)

call rsub('JLANDFLUX',1)

! Zero out sea-cell and land-cell copies of atmospheric properties prior 
! to summation over flux cells.
! THUS, RHOS, VELS, PRSS, AND USTAR ARE ONLY SUMMED OVER SPACE, BUT NOT OVER TIME.

! (ON THE OTHER HAND, LAND%SXFER_T, LAND%SXFER_R, SEA%SXFER_T, AND SEA%SXFER_R
! ARE SUMMED OVER BOTH SPACE AND TIME; THEY ARE NOT RESET TO ZERO HERE BUT INSTEAD
! IN LEAF OR SEA AFTER THEY ARE TRANSFERRED TO THE LEAF OR SEA CANOPY.)

! Loop over all SEA cells

do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

   if (iparallel == 1 .and. itab_ws(iws)%irank /= myrank) cycle

   sea%rhos     (iws) = 0.
   sea%ustar    (iws) = 0.
   sea%sea_ustar(iws) = 0.
   sea%ice_ustar(iws) = 0.
   sea%ggaer    (iws) = 0.
   sea%sea_ggaer(iws) = 0.
   sea%ice_ggaer(iws) = 0.

enddo

! Loop over ALL LAND cells

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   land%ustar  (iwl) = 0.
   land%rhos   (iwl) = 0.
   land%prss   (iwl) = 0.
   land%vels   (iwl) = 0.
   land%ggaer  (iwl) = 0.
   land%ed_zeta(iwl) = 0.
   land%ed_rib (iwl) = 0.
enddo

! Sum fluxes for SEA cells

call psub()
!----------------------------------------------------------------------
! Do parallel recv of TURBULENT fluxes and ATM properties for SEA cells

if (iparallel == 1) call mpi_recv_wsf('T',mrl)

do j = 1,jseaflux(2)%jend(mrl)
   isf = jseaflux(2)%iseaflux(j)
   iws = seaflux(isf)%iwls        ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iws = itabg_ws(iws)%iws_myrank
   endif
!----------------------------------------------------------------------
   call qsub('SF',isf)

   sea%rhos(iws)        = sea%rhos(iws)        + seaflux(isf)%rhos

   sea%ustar    (iws)   = sea%ustar    (iws)   + seaflux(isf)%ustar
   sea%sea_ustar(iws)   = sea%sea_ustar(iws)   + seaflux(isf)%sea_ustar
   sea%ice_ustar(iws)   = sea%ice_ustar(iws)   + seaflux(isf)%ice_ustar

   sea%sxfer_t    (iws) = sea%sxfer_t    (iws) + seaflux(isf)%sxfer_t
   sea%sea_sxfer_t(iws) = sea%sea_sxfer_t(iws) + seaflux(isf)%sea_sxfer_t
   sea%ice_sxfer_t(iws) = sea%ice_sxfer_t(iws) + seaflux(isf)%ice_sxfer_t

   sea%sxfer_r    (iws) = sea%sxfer_r    (iws) + seaflux(isf)%sxfer_r
   sea%sea_sxfer_r(iws) = sea%sea_sxfer_r(iws) + seaflux(isf)%sea_sxfer_r
   sea%ice_sxfer_r(iws) = sea%ice_sxfer_r(iws) + seaflux(isf)%ice_sxfer_r

   sea%sxfer_c(iws)     = sea%sxfer_c(iws)     + seaflux(isf)%sxfer_c

   sea%ggaer    (iws)   = sea%ggaer    (iws)   + seaflux(isf)%ggaer
   sea%sea_ggaer(iws)   = sea%sea_ggaer(iws)   + seaflux(isf)%sea_ggaer
   sea%ice_ggaer(iws)   = sea%ice_ggaer(iws)   + seaflux(isf)%ice_ggaer

enddo
call rsub('JSEAFLUX',2)

! Sum atmospheric properties and fluxes for land cells

call psub()
!----------------------------------------------------------------------
! Do parallel recv of landflux TURBULENT fluxes and ATM properties for land cells

if (iparallel == 1) call mpi_recv_wlf('T',mrl)

do j = 1,jlandflux(2)%jend(mrl)
   ilf = jlandflux(2)%ilandflux(j)
   iwl = landflux(ilf)%iwls        ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iwl = itabg_wl(iwl)%iwl_myrank
   endif

!----------------------------------------------------------------------
   call qsub('LF',iw)

   land%rhos     (iwl) = land%rhos     (iwl) + landflux(ilf)%rhos
   land%vels     (iwl) = land%vels     (iwl) + landflux(ilf)%vels
   land%prss     (iwl) = land%prss     (iwl) + landflux(ilf)%prss
   land%sxfer_t  (iwl) = land%sxfer_t  (iwl) + landflux(ilf)%sxfer_t
   land%sxfer_r  (iwl) = land%sxfer_r  (iwl) + landflux(ilf)%sxfer_r
   land%sxfer_c  (iwl) = land%sxfer_c  (iwl) + landflux(ilf)%sxfer_c
   land%ustar    (iwl) = land%ustar    (iwl) + landflux(ilf)%ustar
   land%ggaer    (iwl) = land%ggaer    (iwl) + landflux(ilf)%ggaer
   land%ed_zeta  (iwl) = land%ed_zeta  (iwl) + landflux(ilf)%ed_zeta
   land%ed_rib   (iwl) = land%ed_rib   (iwl) + landflux(ilf)%ed_rib

enddo
call rsub('JLANDFLUX',2)

! Compute some derived surface quantities

do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

   ka = lpw(iw)
   
   wtv0(iw) = sflux_t(iw) * (1. + .61 * sh_v(ka,iw)) &
            + sflux_r(iw) * .61 * theta(ka,iw)

   if (wtv0(iw) > 0.0) then
      wstar(iw) = (grav * pblh(iw) * wtv0(iw) / theta(ka,iw)) ** 0.33333333
   else
      wstar(iw) = 0.0
   endif

enddo

return
end subroutine surface_turb_flux

!===============================================================================

subroutine stars(zts, rough, vels, rhos, airtemp, sh_vs, cantemp, &
                 can_shv, vkmsfc, sfluxt, sfluxr, ustar0, ggaero  )

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
real, intent(in) :: can_shv   ! canopy air vapor spec hum [kg_vap/m^3]

real(kind=8), intent(in) :: rhos  ! atmos near-surface density [kg/m^3]

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
rstar = c3 * (sh_vs - can_shv)

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

use mem_sflux,   only: jlandflux, landflux
use mem_cuparm,  only: conprr
use mem_leaf,    only: land, itabg_wl, itab_wl
use mem_ijtabs,  only: istp, mrl_begl, itabg_w
use consts_coms, only: cliq, alli, t00
use mem_basic,   only: tair, rho, sh_v
use leaf_coms,   only: mwl, dt_leaf, mrl_leaf
use misc_coms,   only: io6, iparallel
use mem_para,    only: myrank
use ed_misc_coms,only: ed2_active

implicit none

integer :: ilf
integer :: mrl
integer :: j
integer :: iw
integer :: iwl
integer :: kw

real :: arf_land_dtf
real :: airtempc
real :: tempc

real, external :: rhovsl

! Subroutine to transfer atmospheric cumulus parameterization 
! precipitation FLUX to leaf land cells

! Set land fluxes to be done for CUPARM
!    1. for mrl = 0, no fluxes 
!    2. For mrl > mrl_leaf, no fluxes
!    3. For mrl <= mrl_leaf, do fluxes at beginning of long timestep 
!                            for all points in mrl = 1

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0 .and. mrl <= mrl_leaf) then

   mrl = 1  ! special mrl set for leaf

   do j = 1,jlandflux(1)%jend(mrl)
      ilf = jlandflux(1)%ilandflux(j)
!----------------------------------------------------------------------
   call qsub('LF',iw)

   iw = landflux(ilf)%iw  ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif

   kw = landflux(ilf)%kw

! dt_leaf is correct timestep in next line since this subroutine is called at
! same frequency as leaf timestep   

   arf_land_dtf = landflux(ilf)%arf_sfc * dt_leaf

! Compute air temperature in C

   airtempc = tair(kw,iw) - t00

! Estimate wet bulb temp using computation from subroutine each_column in micphys
! Assume that convective precip reaches surface at this wet bulb temp

   tempc = airtempc - min(25.,  &
       700. * (rhovsl(airtempc) / real(rho(kw,iw)) - sh_v(kw,iw)))

   landflux(ilf)%pcpg  = arf_land_dtf * conprr(iw)
   landflux(ilf)%qpcpg = arf_land_dtf * conprr(iw) * (cliq * tempc + alli)
   landflux(ilf)%dpcpg = arf_land_dtf * conprr(iw) * .001

enddo

! Do parallel send/recv of CUPARM FLUXES

if (iparallel == 1) then
  call mpi_send_wlf('C',mrl)
  call mpi_recv_wlf('C',mrl)
endif

endif
call rsub('JLANDFLUX_cuparm',1)

! Sum landflux cell precipitation to get land cell precipitation

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0 .and. mrl <= mrl_leaf) then

mrl = 1  ! special mrl set for leaf

do j = 1,jlandflux(2)%jend(mrl)
   ilf = jlandflux(2)%ilandflux(j)
!----------------------------------------------------------------------
   call qsub('LF',iw)

   iwl = landflux(ilf)%iwls ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iwl = itabg_wl(iwl)%iwl_myrank
   endif

   land%pcpg (iwl) = land%pcpg (iwl) + landflux(ilf)%pcpg
   land%qpcpg(iwl) = land%qpcpg(iwl) + landflux(ilf)%qpcpg
   land%dpcpg(iwl) = land%dpcpg(iwl) + landflux(ilf)%dpcpg

enddo

#ifdef USE_ED2
if (ed2_active == 1) then
   call copy_cuparm_to_ed()
endif
#endif

endif
call rsub('JLANDFLUX_cuparm',2)

! WILL NEED TO ADD SEA LOOP HERE WHEN OCEAN MODEL IS COUPLED WITH OLAM

return
end subroutine surface_cuparm_flux

!===============================================================================

subroutine surface_precip_flux()

use mem_sflux,   only: jlandflux, landflux
use mem_micro,   only: pcpgr, qpcpgr, dpcpgr
use mem_leaf,    only: land, itabg_wl, itab_wl
use mem_ijtabs,  only: istp, mrl_endl, itabg_w
use leaf_coms,   only: mwl, mrl_leaf
use misc_coms,   only: io6, iparallel
use mem_para,    only: myrank
use ed_misc_coms,only: ed2_active

implicit none

integer :: j
integer :: ilf
integer :: iw
integer :: iwl
integer :: mrl

real :: arf_land
real :: arf_land_dtf

! Subroutine to transfer atmospheric microphysics parameterization 
! precipitation flux to leaf land cells.

! Land fluxes to be done for PRECIP 
!    1. for mrl = 0, no fluxes 
!    2. For mrl > mrl_leaf, do fluxes at end of long timestep at given mrl
!    3. For mrl <= mrl_leaf, do fluxes at end of long timestep
!                            for all points in mrl = 1

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then

if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

do j = 1,jlandflux(1)%jend(mrl)
   ilf = jlandflux(1)%ilandflux(j)
!----------------------------------------------------------------------
call qsub('LF',iw)

   iw = landflux(ilf)%iw  ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iw = itabg_w(iw)%iw_myrank
   endif

   arf_land_dtf = landflux(ilf)%arf_sfc * landflux(ilf)%dtf

   landflux(ilf)%pcpg  = arf_land_dtf * pcpgr(iw)
   landflux(ilf)%qpcpg = arf_land_dtf * qpcpgr(iw)
   landflux(ilf)%dpcpg = arf_land_dtf * dpcpgr(iw)

enddo

! Do parallel send/recv of MICPHYS PRECIP FLUXES

if (iparallel == 1) then
   call mpi_send_wlf('M',mrl)
   call mpi_recv_wlf('M',mrl)
endif

endif
call rsub('JLANDFLUX_pcp',1)

! Sum landflux cell precipitation to get land cell precipitation

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then

if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

do j = 1,jlandflux(2)%jend(mrl)
   ilf = jlandflux(2)%ilandflux(j)
!----------------------------------------------------------------------
call qsub('LF',iw)

   iwl = landflux(ilf)%iwls ! global index

! If run is parallel, get local rank indices

   if (iparallel == 1) then
      iwl = itabg_wl(iwl)%iwl_myrank
   endif

   land%pcpg (iwl) = land%pcpg (iwl) + landflux(ilf)%pcpg
   land%qpcpg(iwl) = land%qpcpg(iwl) + landflux(ilf)%qpcpg
   land%dpcpg(iwl) = land%dpcpg(iwl) + landflux(ilf)%dpcpg

enddo

#ifdef USE_ED2
if (ed2_active == 1) then
   call copy_micro_to_ed()
endif
#endif

endif
call rsub('JLANDFLUX_pcp',2)

! WILL NEED TO ADD SEA LOOP HERE WHEN OCEAN MODEL IS COUPLED WITH OLAM

! pcpgr, qpcpgr, dpcpgr have now been "transferred" to leaf cells.
! No need to zero pcpgr, qpcpgr, dpcpgr since microphysics will replace them.

return
end subroutine surface_precip_flux
