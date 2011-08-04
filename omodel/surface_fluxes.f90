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
use mem_ijtabs,  only: itab_w, itabg_w, jtab_w
use misc_coms,   only: io6, iparallel
use mem_para,    only: myrank
use mem_grid,    only: mza, mwa, lsw, lpw, zt, zm
use mem_sea,     only: sea, itabg_ws, itab_ws
use mem_leaf,    only: land, itabg_wl, itab_wl, first_site
use consts_coms, only: p00i, rocp
use mem_sflux,   only: seaflux, landflux, jseaflux, jlandflux
use mem_turb,    only: vkm_sfc, sflux_t, sflux_r, sxfer_tk, sxfer_rk, &
                       vels, ustar
use mem_basic,   only: press, rho, theta, sh_v

use ed_structure_defs

implicit none

integer, intent(in) :: mrl

integer :: j,jv,iw,iv,iws,iwl
integer :: ks,kw,ka,k
integer :: isf,ilf

real :: area_dtf
real :: arf_atm
real :: arf_sea
real :: arf_land
real :: arf_sea_dtf
real :: arf_land_dtf
real :: exner
real :: airtemp
real :: vkmsfc
real :: ustar0

type(site), pointer :: ed_site
type(patch), pointer :: ed_patch

if (isfcl == 0) then

! ISFCL = 0 is the no-LEAF option:  Set all surface fluxes (to zero by default).
! These flux values can be altered if desired.

   call psub()
!-----------------------------------------------------------------------------
   do j = 1,jtab_w(13)%jend(mrl); iw = jtab_w(13)%iw(j)
!-----------------------------------------------------------------------------
   call qsub('W',iw)

      vkm_sfc(iw) = 0.
      sflux_t(iw) = 0.
      sflux_r(iw) = 0.
      ustar  (iw) = .1  ! Minimum value

      do ks = 1,lsw(iw)
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
do j = 1,jtab_w(20)%jend(mrl); iw = jtab_w(20)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   ka = lpw(iw)

   vkm_sfc(iw) = 0.
   ustar(iw) = 0.
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

   exner = (p00i * press(kw,iw)) ** rocp
   airtemp = theta(kw,iw) * exner

   call stars(zt(kw)-zm(kw-1),     &
              sea%rough(iws),      &
              vels        (iw),    &
              rho      (kw,iw),    &
              airtemp,             &
              sh_v     (kw,iw),    &
              sea%can_temp(iws),   &
              sea%can_shv(iws),    &
              vkmsfc,              &
              seaflux(isf)%sfluxt, &
              seaflux(isf)%sfluxr, &
              ustar0               )

! Add flux contributions to IW atmospheric column

   vkm_sfc(iw) = vkm_sfc(iw) + arf_atm * vkmsfc
   ustar  (iw) = ustar  (iw) + arf_atm * ustar0
   sflux_t(iw) = sflux_t(iw) + arf_atm * seaflux(isf)%sfluxt / exner ! for Taylor PBL
   sflux_r(iw) = sflux_r(iw) + arf_atm * seaflux(isf)%sfluxr         ! for Taylor PBL
   sxfer_tk(ks,iw) = sxfer_tk(ks,iw) + area_dtf * seaflux(isf)%sfluxt / exner
   sxfer_rk(ks,iw) = sxfer_rk(ks,iw) + area_dtf * seaflux(isf)%sfluxr

! Store flux and ATM density contributions to SEA cell in SEAFLUX cell

   seaflux(isf)%sxfer_t = arf_sea_dtf * seaflux(isf)%sfluxt
   seaflux(isf)%sxfer_r = arf_sea_dtf * seaflux(isf)%sfluxr
   seaflux(isf)%ustar   = arf_sea     * ustar0
   seaflux(isf)%rhos    = arf_sea     * rho(kw,iw)
      
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

! Store atmospheric properties in landflux cell

   landflux(ilf)%vels = arf_land * vels(iw)
   landflux(ilf)%prss = arf_land * press(kw,iw)
   landflux(ilf)%rhos = arf_land * rho(kw,iw)

! Calculate turbulent fluxes between atmosphere and land canopy

   exner = (p00i * press(kw,iw)) ** rocp
   airtemp = theta(kw,iw) * exner

   if (land%ed_flag(iwl) == 0) then

      call stars(zt(kw)-zm(kw-1),      &
                 land%rough    (iwl),  &
                 vels           (iw),  &
                 rho         (kw,iw),  &
                 airtemp,              &
                 sh_v        (kw,iw),  &
                 land%can_temp (iwl),  &
                 land%can_shv  (iwl),  &
                 vkmsfc,               &
                 landflux(ilf)%sfluxt, &
                 landflux(ilf)%sfluxr, &
                 ustar0                )

! Add flux contributions to IW atmospheric column

      vkm_sfc(iw) = vkm_sfc(iw)  &
         + arf_atm * vkmsfc

      ustar  (iw) = ustar  (iw)  &
         + arf_atm * ustar0

      sflux_t(iw) = sflux_t(iw)  &
         + arf_atm * landflux(ilf)%sfluxt / exner ! for Taylor PBL

      sflux_r(iw) = sflux_r(iw)  &
         + arf_atm * landflux(ilf)%sfluxr         ! for Taylor PBL

      sxfer_tk(ks,iw) = sxfer_tk(ks,iw)  &
         + area_dtf * landflux(ilf)%sfluxt / exner

      sxfer_rk(ks,iw) = sxfer_rk(ks,iw)  &
         + area_dtf * landflux(ilf)%sfluxr

! Store flux contributions to land cell in landflux cell

      landflux(ilf)%sxfer_t = arf_land_dtf * landflux(ilf)%sfluxt
      landflux(ilf)%sxfer_r = arf_land_dtf * landflux(ilf)%sfluxr
      landflux(ilf)%ustar   = arf_land     * ustar0

   else

! Zero landflux values prior to summation over ED patches

      landflux(ilf)%sxfer_t = 0.
      landflux(ilf)%sxfer_r = 0.
      landflux(ilf)%ustar   = 0.

! Find the site

      ed_site => first_site
      find_site: do while(associated(ed_site))
         if(ed_site%iland == iwl)exit find_site
         ed_site => ed_site%next_site
      enddo find_site
      ed_patch => ed_site%oldest_patch
      do while(associated(ed_patch))

         call stars(zt(kw)-zm(kw-1), &
              ed_patch%rough  ,      &
              vels           (iw),   &
              rho         (kw,iw),   &
              airtemp,               &
              sh_v        (kw,iw),   &
              ed_patch%can_temp,     &
              ed_patch%can_shv ,     &
              vkmsfc,                &
              landflux(ilf)%sfluxt,  &
              landflux(ilf)%sfluxr,  &
              ustar0                 )

! Add flux contributions to IW atmospheric column

         vkm_sfc(iw) = vkm_sfc(iw) + arf_atm * vkmsfc * ed_patch%area
         ustar  (iw) = ustar  (iw) + arf_atm * ustar0 * ed_patch%area

! sflux_t(iw) and sflux_r(iw) are for Taylor PBL

         sflux_t(iw) = sflux_t(iw)  &
            + arf_atm * landflux(ilf)%sfluxt * ed_patch%area / exner

         sflux_r(iw) = sflux_r(iw)  &
            + arf_atm * landflux(ilf)%sfluxr * ed_patch%area

         sxfer_tk(ks,iw) = sxfer_tk(ks,iw)  &
            + area_dtf * landflux(ilf)%sfluxt * ed_patch%area / exner

         sxfer_rk(ks,iw) = sxfer_rk(ks,iw)  &
            + area_dtf * landflux(ilf)%sfluxr * ed_patch%area

!----------------------------------------------------------------------------------
! Conceptually, we need the following here: 3 arrays over all patches in current 
! IWL land cell for each flux cell.

     !    landflux(ilf)%sxfer_t(ed_patch) = arf_land_dtf * landflux(ilf)%sfluxt
     !    landflux(ilf)%sxfer_r(ed_patch) = arf_land_dtf * landflux(ilf)%sfluxr
     !    landflux(ilf)%ustar(ed_patch)   = arf_land     * ustar0
!----------------------------------------------------------------------------------

! Add flux contributions to landflux cell

         landflux(ilf)%sxfer_t = landflux(ilf)%sxfer_t  &
            + arf_land_dtf * landflux(ilf)%sfluxt * ed_patch%area

         landflux(ilf)%sxfer_r = landflux(ilf)%sxfer_r  &
            + arf_land_dtf * landflux(ilf)%sfluxr * ed_patch%area

         landflux(ilf)%ustar   = landflux(ilf)%ustar    &
            + arf_land * ustar0 * ed_patch%area

         ed_patch => ed_patch%younger
      enddo

   endif

enddo

! Do parallel send of TURBULENT fluxes and ATM properties

if (iparallel == 1) call mpi_send_wlf('T',mrl)

call rsub('JLANDFLUX',1)

! Zero out sea-cell and land-cell copies of atmospheric properties prior 
! to summation over flux cells.
! THUS, RHOS, VELS, PRSS, AND USTAR ARE ONLY SUMMED OVER SPACE, BUT NOT OVER TIME.

! (ON THE OTHER HAND, LAND%SXFER_T, LAND%SXFER_R, SEA%SXFER_T, AND SEA%SXFER_R
! ARE SUMMED OVER BOTH SPACE AND TIME; THEY ARE NOT RESET TO ZERO HERE BUT INSTEAD
! IN LEAF3 OR SEA AFTER THEY ARE TRANSFERRED TO THE LEAF3 OR SEA CANOPY.)

! Loop over all SEA cells

do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

   if (iparallel == 1 .and. itab_ws(iws)%irank /= myrank) cycle

   sea%ustar(iws) = 0.
   sea%rhos(iws)  = 0.
enddo

! Loop over ALL LAND cells

do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

   if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   land%ustar(iwl) = 0.
   land%rhos (iwl) = 0.
   land%prss (iwl) = 0.
   land%vels (iwl) = 0.
enddo

ed_site => first_site
do while(associated(ed_site))
   ed_patch => ed_site%oldest_patch
   do while(associated(ed_patch))
      ed_patch%ustar = 0.
      ed_patch => ed_patch%younger
   enddo
   ed_site => ed_site%next_site
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

   sea%sxfer_t(iws) = sea%sxfer_t(iws) + seaflux(isf)%sxfer_t
   sea%sxfer_r(iws) = sea%sxfer_r(iws) + seaflux(isf)%sxfer_r
   sea%rhos(iws)    = sea%rhos(iws)    + seaflux(isf)%rhos
   sea%ustar(iws)   = sea%ustar(iws)   + seaflux(isf)%ustar

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

   land%rhos   (iwl) = land%rhos   (iwl) + landflux(ilf)%rhos
   land%vels   (iwl) = land%vels   (iwl) + landflux(ilf)%vels
   land%prss   (iwl) = land%prss   (iwl) + landflux(ilf)%prss
   land%sxfer_t(iwl) = land%sxfer_t(iwl) + landflux(ilf)%sxfer_t
   land%sxfer_r(iwl) = land%sxfer_r(iwl) + landflux(ilf)%sxfer_r
   land%ustar  (iwl) = land%ustar  (iwl) + landflux(ilf)%ustar

   if (land%ed_flag(iwl) /= 0) then

! Find the site

      ed_site => first_site
      find_site2: do while(associated(ed_site))
         if(ed_site%iland == iwl)exit find_site2
         ed_site => ed_site%next_site
      enddo find_site2
      ed_patch => ed_site%oldest_patch
      do while(associated(ed_patch))

!----------------------------------------------------------------------------------
! Conceptually, we need the following here: 3 arrays over all patches in current 
! IWL land cell for each flux cell.

  !       ed_patch%sxfer_t = ed_patch%sxfer_t + landflux(ilf)%sxfer_t(ed_patch) 
  !       ed_patch%sxfer_r = ed_patch%sxfer_r + landflux(ilf)%sxfer_r(ed_patch)
  !       ed_patch%ustar   = ed_patch%ustar   + landflux(ilf)%ustar(ed_patch)
!----------------------------------------------------------------------------------

      enddo

   endif

enddo
call rsub('JLANDFLUX',2)

return
end subroutine surface_turb_flux

!===============================================================================

subroutine stars(zts, rough, vels, rhos, airtemp, sh_vs,        &
                 cantemp, can_shv, vkmsfc, sfluxt, sfluxr, ustar0)

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
real :: ri     ! bulk richardson numer, eq. 3.45 in Garratt
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

return
end subroutine stars

!===============================================================================

subroutine surface_cuparm_flux()

use mem_sflux,   only: jlandflux, landflux
use mem_cuparm,  only: conprr
use mem_leaf,    only: land, itabg_wl, itab_wl
use mem_ijtabs,  only: istp, mrl_begl, itabg_w
use consts_coms, only: p00i, rocp, cliq, alli
use mem_basic,   only: press, theta, rho, sh_v
use leaf_coms,   only: mwl, dt_leaf, mrl_leaf
use misc_coms,   only: io6, iparallel
use mem_para,    only: myrank

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
! precipitation FLUX to leaf3 land cells

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

   airtempc = theta(kw,iw) * (p00i * press(kw,iw)) ** rocp - 273.15

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
endif
call rsub('JLANDFLUX_cuparm',2)

! WILL NEED TO ADD SEA LOOP HERE WHEN OCEAN MODEL IS COUPLED WITH OLAM

return
end subroutine surface_cuparm_flux

!===============================================================================

subroutine surface_precip_flux()

use mem_sflux,  only: jlandflux, landflux
use mem_micro,  only: pcpgr, qpcpgr, dpcpgr
use mem_leaf,   only: land, itabg_wl, itab_wl
use mem_ijtabs, only: istp, mrl_endl, itabg_w
use leaf_coms,  only: mwl, mrl_leaf
use misc_coms,  only: io6, iparallel
use mem_para,   only: myrank

implicit none

integer :: j
integer :: ilf
integer :: iw
integer :: iwl
integer :: mrl

real :: arf_land
real :: arf_land_dtf

! Subroutine to transfer atmospheric microphysics parameterization 
! precipitation flux to leaf3 land cells.

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
endif
call rsub('JLANDFLUX_pcp',2)

! WILL NEED TO ADD SEA LOOP HERE WHEN OCEAN MODEL IS COUPLED WITH OLAM

! pcpgr, qpcpgr, dpcpgr have now been "transferred" to leaf cells.
! No need to zero pcpgr, qpcpgr, dpcpgr since microphysics will replace them.

return
end subroutine surface_precip_flux



