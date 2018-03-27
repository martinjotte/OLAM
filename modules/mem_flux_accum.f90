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

Module mem_flux_accum

  use consts_coms, only: r8

! Atmosphere grid arrays

  real(r8), allocatable, target :: rshort_accum      (:)
  real(r8), allocatable, target :: rshortup_accum    (:)
  real(r8), allocatable, target :: rlong_accum       (:)
  real(r8), allocatable, target :: rlongup_accum     (:)
  real(r8), allocatable, target :: rshort_top_accum  (:)
  real(r8), allocatable, target :: rshortup_top_accum(:)
  real(r8), allocatable, target :: rlongup_top_accum (:)

  real(r8), allocatable, target :: rshort_clr_accum      (:)
  real(r8), allocatable, target :: rshortup_clr_accum    (:)
  real(r8), allocatable, target :: rlong_clr_accum       (:)
  real(r8), allocatable, target :: rlongup_clr_accum     (:)
  real(r8), allocatable, target :: rshort_top_clr_accum  (:)
  real(r8), allocatable, target :: rshortup_top_clr_accum(:)
  real(r8), allocatable, target :: rlongup_top_clr_accum (:)

  real(r8), allocatable, target :: sfluxt_accum(:)
  real(r8), allocatable, target :: sfluxr_accum(:)

  real(r8), allocatable, target :: vc_accum   (:,:)
  real(r8), allocatable, target :: wc_accum   (:,:)
  real(r8), allocatable, target :: press_accum(:,:)
  real(r8), allocatable, target :: tair_accum (:,:)
  real(r8), allocatable, target :: sh_v_accum (:,:)

! Land arrays

  real(r8), allocatable, target :: vels_l_accum    (:)
  real(r8), allocatable, target :: airtemp_l_accum (:)
  real(r8), allocatable, target :: airshv_l_accum  (:)
  real(r8), allocatable, target :: cantemp_l_accum (:)
  real(r8), allocatable, target :: canshv_l_accum  (:)
  real(r8), allocatable, target :: skintemp_l_accum(:)
  real(r8), allocatable, target :: sfluxt_l_accum  (:)
  real(r8), allocatable, target :: sfluxr_l_accum  (:)
  real(r8), allocatable, target :: wxfer1_l_accum  (:) ! filled in leaf4_soil

! Sea arrays

  real(r8), allocatable, target :: vels_s_accum    (:)
  real(r8), allocatable, target :: airtemp_s_accum (:)
  real(r8), allocatable, target :: airshv_s_accum  (:)
  real(r8), allocatable, target :: cantemp_s_accum (:)
  real(r8), allocatable, target :: canshv_s_accum  (:)
  real(r8), allocatable, target :: skintemp_s_accum(:)
  real(r8), allocatable, target :: sfluxt_s_accum  (:)
  real(r8), allocatable, target :: sfluxr_s_accum  (:)


Contains

!===============================================================================

subroutine alloc_flux_accum(mza,mva,mwa,mwl,mws)

  use misc_coms, only: ilwrtyp, iswrtyp
  use leaf_coms, only: isfcl
  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: mza, mva, mwa, mwl, mws

! Allocate arrays for accumulated flux quantities
! Initialize arrays to zero

  if (iswrtyp > 0) then
     allocate (rshort_accum      (mwa)) ; rshort_accum       = 0._r8
     allocate (rshortup_accum    (mwa)) ; rshortup_accum     = 0._r8
     allocate (rshort_top_accum  (mwa)) ; rshort_top_accum   = 0._r8
     allocate (rshortup_top_accum(mwa)) ; rshortup_top_accum = 0._r8
  endif

  if (ilwrtyp > 0) then
     allocate (rlong_accum       (mwa)) ; rlong_accum        = 0._r8
     allocate (rlongup_accum     (mwa)) ; rlongup_accum      = 0._r8
     allocate (rlongup_top_accum (mwa)) ; rlongup_top_accum  = 0._r8
  endif

  if (.true.) then
     allocate (rshort_clr_accum      (mwa)) ; rshort_clr_accum       = 0._r8
     allocate (rshortup_clr_accum    (mwa)) ; rshortup_clr_accum     = 0._r8
     allocate (rlong_clr_accum       (mwa)) ; rlong_clr_accum        = 0._r8
     allocate (rlongup_clr_accum     (mwa)) ; rlongup_clr_accum      = 0._r8
     allocate (rshort_top_clr_accum  (mwa)) ; rshort_top_clr_accum   = 0._r8
     allocate (rshortup_top_clr_accum(mwa)) ; rshortup_top_clr_accum = 0._r8
     allocate (rlongup_top_clr_accum (mwa)) ; rlongup_top_clr_accum  = 0._r8
  endif

  allocate (vc_accum   (mza,mva)) ; vc_accum    = 0._r8
  allocate (wc_accum   (mza,mwa)) ; wc_accum    = 0._r8
  allocate (press_accum(mza,mwa)) ; press_accum = 0._r8
  allocate (tair_accum (mza,mwa)) ; tair_accum  = 0._r8
  allocate (sh_v_accum (mza,mwa)) ; sh_v_accum  = 0._r8

  if (isfcl > 0) then
     allocate (sfluxt_accum    (mwa)) ; sfluxt_accum     = 0._r8
     allocate (sfluxr_accum    (mwa)) ; sfluxr_accum     = 0._r8

     allocate (vels_l_accum    (mwl)) ; vels_l_accum     = 0._r8
     allocate (airtemp_l_accum (mwl)) ; airtemp_l_accum  = 0._r8
     allocate (airshv_l_accum  (mwl)) ; airshv_l_accum   = 0._r8
     allocate (cantemp_l_accum (mwl)) ; cantemp_l_accum  = 0._r8
     allocate (canshv_l_accum  (mwl)) ; canshv_l_accum   = 0._r8
     allocate (skintemp_l_accum(mwl)) ; skintemp_l_accum = 0._r8
     allocate (sfluxt_l_accum  (mwl)) ; sfluxt_l_accum   = 0._r8
     allocate (sfluxr_l_accum  (mwl)) ; sfluxr_l_accum   = 0._r8
     allocate (wxfer1_l_accum  (mwl)) ; wxfer1_l_accum   = 0._r8

     allocate (vels_s_accum    (mws)) ; vels_s_accum     = 0._r8
     allocate (airtemp_s_accum (mws)) ; airtemp_s_accum  = 0._r8
     allocate (airshv_s_accum  (mws)) ; airshv_s_accum   = 0._r8
     allocate (cantemp_s_accum (mws)) ; cantemp_s_accum  = 0._r8
     allocate (canshv_s_accum  (mws)) ; canshv_s_accum   = 0._r8
     allocate (skintemp_s_accum(mws)) ; skintemp_s_accum = 0._r8
     allocate (sfluxt_s_accum  (mws)) ; sfluxt_s_accum   = 0._r8
     allocate (sfluxr_s_accum  (mws)) ; sfluxr_s_accum   = 0._r8
  endif
  
end subroutine alloc_flux_accum

!===============================================================================

subroutine filltab_flux_accum()

  use var_tables, only: increment_vtable
  implicit none

  if (allocated (rshort_accum)) call increment_vtable('RSHORT_ACCUM', 'AW', dvar1=rshort_accum)

  if (allocated(rshortup_accum)) call increment_vtable('RSHORTUP_ACCUM', 'AW', dvar1=rshortup_accum)

  if (allocated(rlong_accum)) call increment_vtable('RLONG_ACCUM', 'AW', dvar1=rlong_accum)

  if (allocated(rlongup_accum)) call increment_vtable('RLONGUP_ACCUM', 'AW', dvar1=rlongup_accum)

  if (allocated(rshort_top_accum)) call increment_vtable('RSHORT_TOP_ACCUM', 'AW', dvar1=rshort_top_accum)

  if (allocated(rshortup_top_accum)) call increment_vtable('RSHORTUP_TOP_ACCUM','AW', dvar1=rshortup_top_accum)

  if (allocated(rlongup_top_accum)) call increment_vtable('RLONGUP_TOP_ACCUM', 'AW', dvar1=rlongup_top_accum)

  if (allocated (rshort_accum)) call increment_vtable('RSHORT_CLR_ACCUM', 'AW', dvar1=rshort_clr_accum)

  if (allocated(rshortup_accum)) call increment_vtable('RSHORTUP_CLR_ACCUM', 'AW', dvar1=rshortup_clr_accum)

  if (allocated(rlong_accum)) call increment_vtable('RLONG_CLR_ACCUM', 'AW', dvar1=rlong_clr_accum)

  if (allocated(rlongup_accum)) call increment_vtable('RLONGUP_CLR_ACCUM', 'AW', dvar1=rlongup_clr_accum)

  if (allocated(rshort_top_accum)) call increment_vtable('RSHORT_TOP_CLR_ACCUM', 'AW', dvar1=rshort_top_clr_accum)

  if (allocated(rshortup_top_accum)) call increment_vtable('RSHORTUP_TOP_CLR_ACCUM','AW', dvar1=rshortup_top_clr_accum)

  if (allocated(rlongup_top_accum)) call increment_vtable('RLONGUP_TOP_CLR_ACCUM', 'AW', dvar1=rlongup_top_clr_accum)

  if (allocated(sfluxt_accum)) call increment_vtable('SFLUXT_ACCUM', 'AW', dvar1=sfluxt_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('SFLUXR_ACCUM', 'AW', dvar1=sfluxr_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('VC_ACCUM', 'AV', dvar2=vc_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('WC_ACCUM', 'AW', dvar2=wc_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('PRESS_ACCUM', 'AW', dvar2=press_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('TAIR_ACCUM', 'AW', dvar2=tair_accum)

  if (allocated(sfluxr_accum)) call increment_vtable('SH_V_ACCUM', 'AW', dvar2=sh_v_accum)

  if (allocated(vels_l_accum)) call increment_vtable('VELS_L_ACCUM', 'LW', dvar1=vels_l_accum)

  if (allocated(airtemp_l_accum)) call increment_vtable('AIRTEMP_L_ACCUM', 'LW', dvar1=airtemp_l_accum)

  if (allocated(airshv_l_accum)) call increment_vtable('AIRSHV_L_ACCUM', 'LW', dvar1=airshv_l_accum)

  if (allocated(cantemp_l_accum)) call increment_vtable('CANTEMP_L_ACCUM', 'LW', dvar1=cantemp_l_accum)

  if (allocated(canshv_l_accum)) call increment_vtable('CANSHV_L_ACCUM', 'LW', dvar1=canshv_l_accum)

  if (allocated(skintemp_l_accum)) call increment_vtable('SKINTEMP_L_ACCUM', 'LW', dvar1=skintemp_l_accum)

  if (allocated(sfluxt_l_accum)) call increment_vtable('SFLUXT_L_ACCUM', 'LW', dvar1=sfluxt_l_accum)

  if (allocated(sfluxr_l_accum)) call increment_vtable('SFLUXR_L_ACCUM', 'LW', dvar1=sfluxr_l_accum)

  if (allocated(wxfer1_l_accum)) call increment_vtable('WXFER1_L_ACCUM', 'LW', dvar1=wxfer1_l_accum)

  if (allocated(vels_s_accum)) call increment_vtable('VELS_S_ACCUM', 'SW', dvar1=vels_s_accum)

  if (allocated(airtemp_s_accum)) call increment_vtable('AIRTEMP_S_ACCUM', 'SW', dvar1=airtemp_s_accum)

  if (allocated(airshv_s_accum)) call increment_vtable('AIRSHV_S_ACCUM', 'SW', dvar1=airshv_s_accum)

  if (allocated(cantemp_s_accum)) call increment_vtable('CANTEMP_S_ACCUM', 'SW', dvar1=cantemp_s_accum)

  if (allocated(canshv_s_accum)) call increment_vtable('CANSHV_S_ACCUM', 'SW', dvar1=canshv_s_accum)

  if (allocated(skintemp_s_accum)) call increment_vtable('SKINTEMP_S_ACCUM', 'SW', dvar1=skintemp_s_accum)

  if (allocated(sfluxt_s_accum)) call increment_vtable('SFLUXT_S_ACCUM', 'SW', dvar1=sfluxt_s_accum)

  if (allocated(sfluxr_s_accum)) call increment_vtable('SFLUXR_S_ACCUM', 'SW', dvar1=sfluxr_s_accum)

end subroutine filltab_flux_accum

!===============================================================================

subroutine flux_accum()

  use misc_coms,   only: io6, time_istp8, dtlm, dtsm, time_prevhist, &
                         ilwrtyp, iswrtyp, isubdomain

  use mem_ijtabs,  only: istp, itab_v, itab_w, itabg_w, jtab_v, jtab_w, &
                         jtv_prog, jtw_prog, mrl_begl, mrl_begs, mrl_endl

  use mem_basic,   only: vc, wc, press, tair, sh_v, vxe, vye, vze

  use mem_radiate, only: albedt, rshort, rlong, rlongup, &
                         rshort_top, rshortup_top, rlongup_top, &
                         rshort_clr, rshortup_clr, rshort_top_clr, &
                         rshortup_top_clr, rlong_clr, rlongup_clr, &
                         rlongup_top_clr

  use leaf_coms,   only: nzg, mwl, mrl_leaf, dt_leaf, isfcl, slcpd
  use mem_turb,    only: sfluxt, sfluxr
  use consts_coms, only: r8
  use sea_coms,    only: mws
  use mem_sea,     only: sea, itab_ws
  use mem_leaf,    only: land, itab_wl
  use mem_grid,    only: mza, lpv, lpw
  use therm_lib,   only: qtk, qwtk

  implicit none

  integer :: mrl, j, iv, iw, iwl, iws, k, nls, kw

  real :: tempk, fracliq, fldval, wind

  real(r8) :: dta

! Update accumulations of ATM velocity and pressure

!----------------------------------------------------------------------
  mrl = mrl_begs(istp)
  if (mrl > 0) then

     !$omp parallel do private(iv,dta,k)
     do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!----------------------------------------------------------------------

! Timestep for accumulating velocity and pressure, DTA, is the small
! (acoustic) timestep. The frequency that each IV cell is processed in this
! loop should be consistent with (inversely proportional to) that DTA.

        dta = dtsm(itab_v(iv)%mrlv)

        do k = lpv(iv),mza
           vc_accum(k,iv) = vc_accum(k,iv) + dta * real(vc(k,iv),r8)
        enddo

     enddo
     !$omp end parallel do

  endif

!----------------------------------------------------------------------
  mrl = mrl_begs(istp)
  if (mrl > 0) then

     !$omp parallel do private(iw,dta,k)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Timestep for accumulating velocity and pressure, DTA, is the small
! (acoustic) timestep. The frequency that each IW cell is processed in this
! loop should be consistent with (inversely proportional to) that DTA.

        dta = dtsm(itab_w(iw)%mrlw)

        do k = lpw(iw),mza
              wc_accum(k,iw) =    wc_accum(k,iw) + dta * real(wc(k,iw),r8)
           press_accum(k,iw) = press_accum(k,iw) + dta *   press(k,iw)
        enddo

     enddo
     !$omp end parallel do

  endif

! Update accumulations of ATM temperature and specific humidity

!----------------------------------------------------------------------
  mrl = mrl_endl(istp)
  if (mrl > 0) then

     !$omp parallel do private(iw,dta,k)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Timestep for accumulating velocity and pressure, DTA, is the long
! timestep. The frequency that each IW cell is processed in this
! loop should be consistent with (inversely proportional to) that DTA.

        dta = dtlm(itab_w(iw)%mrlw)

        do k = lpw(iw),mza
           tair_accum(k,iw) = tair_accum(k,iw) + dta * real(tair(k,iw),r8)
           sh_v_accum(k,iw) = sh_v_accum(k,iw) + dta * real(sh_v(k,iw),r8)
        enddo

     enddo
     !$omp end parallel do

  endif

! Update accumulations of ATM SFC turbulent fluxes

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

     !$omp parallel do private(iw,dta)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Timestep for accumulating fluxes, DTA, is the lesser of dtlm for the given
! IW cell and dt_leaf; the frequency that each IW cell is processed in this
! loop should be consistent with (inversely proportional to) that DTA.

        dta = min( dtlm(itab_w(iw)%mrlw), real(dt_leaf,r8) )

        sfluxt_accum(iw) = sfluxt_accum(iw) + dta * real(sfluxt(iw),r8)
        sfluxr_accum(iw) = sfluxr_accum(iw) + dta * real(sfluxr(iw),r8)
     enddo
     !$omp end parallel do

  endif

! Update accumulations of ATM SFC and TOP radiative fluxes

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. ilwrtyp + iswrtyp > 0) then

     !$omp parallel do private(iw,dta)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! Timestep for accumulating fluxes, DTA is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DTA.

        dta = dtlm(itab_w(iw)%mrlw)

        if (iswrtyp > 0) then
           rshort_accum      (iw) = rshort_accum      (iw) + dta * real(rshort(iw),r8)
           rshortup_accum    (iw) = rshortup_accum    (iw) + dta * real(rshort(iw) * albedt(iw),r8)
           rshort_top_accum  (iw) = rshort_top_accum  (iw) + dta * real(rshort_top(iw),r8)
           rshortup_top_accum(iw) = rshortup_top_accum(iw) + dta * real(rshortup_top(iw),r8)

           rshort_clr_accum      (iw) = rshort_clr_accum      (iw) + dta * real(rshort_clr(iw),r8)
           rshortup_clr_accum    (iw) = rshortup_clr_accum    (iw) + dta * real(rshortup_clr(iw),r8)
           rshort_top_clr_accum  (iw) = rshort_top_clr_accum  (iw) + dta * real(rshort_top_clr(iw),r8)
           rshortup_top_clr_accum(iw) = rshortup_top_clr_accum(iw) + dta * real(rshortup_top_clr(iw),r8)
        endif

        if (ilwrtyp > 0) then
           rlong_accum      (iw) = rlong_accum      (iw) + dta * real(rlong(iw),r8)
           rlongup_accum    (iw) = rlongup_accum    (iw) + dta * real(rlongup(iw),r8)
           rlongup_top_accum(iw) = rlongup_top_accum(iw) + dta * real(rlongup_top(iw),r8)

           rlong_clr_accum      (iw) = rlong_clr_accum      (iw) + dta * real(rlong_clr(iw),r8)
           rlongup_clr_accum    (iw) = rlongup_clr_accum    (iw) + dta * real(rlongup_clr(iw),r8)
           rlongup_top_clr_accum(iw) = rlongup_top_clr_accum(iw) + dta * real(rlongup_top_clr(iw),r8)
        endif

     enddo
     !$omp end parallel do
  endif

! Update accumulations of LAND cells

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

     !$omp parallel do private(iw,kw,dta,wind,nls,tempk,fracliq,fldval)
     do iwl = 2,mwl
        iw = itab_wl(iwl)%iw         ! global index
        kw = itab_wl(iwl)%kw

! If run is parallel, get local rank indices

        if (isubdomain == 1) then
           iw = itabg_w(iw)%iw_myrank
        endif

! Timestep for accumulating fluxes, DTA is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DTA.

        dta  = dtlm(itab_w(iw)%mrlw)
        wind = sqrt(vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2)

           vels_l_accum(iwl) =    vels_l_accum(iwl) + dta * real(wind,r8)
        airtemp_l_accum(iwl) = airtemp_l_accum(iwl) + dta * real(tair(kw,iw),r8)
         airshv_l_accum(iwl) =  airshv_l_accum(iwl) + dta * real(sh_v(kw,iw),r8)

        cantemp_l_accum(iwl) = cantemp_l_accum(iwl) + dta * real(land%cantemp(iwl),r8)
         canshv_l_accum(iwl) =  canshv_l_accum(iwl) + dta * real(land%canshv (iwl),r8)
         sfluxt_l_accum(iwl) =  sfluxt_l_accum(iwl) + dta * real(land%sfluxt (iwl),r8)
         sfluxr_l_accum(iwl) =  sfluxr_l_accum(iwl) + dta * real(land%sfluxr (iwl),r8)

        nls = land%nlev_sfcwater(iwl)

        if (nls > 0) then
           call qtk(land%sfcwater_energy(nls,iwl),tempk,fracliq)
        else
           call qwtk(land%soil_energy(nzg,iwl)       &
                    ,land%soil_water(nzg,iwl)*1.e3   &
                    ,slcpd(land%ntext_soil(nzg,iwl)) &
                    ,tempk, fracliq)
        endif

        fldval = (1. - land%vf(iwl)) * tempk &
                     + land%vf(iwl)  * land%veg_temp(iwl)

        skintemp_l_accum(iwl) = skintemp_l_accum(iwl) + dta * real(fldval,r8)

     enddo
     !$omp end parallel do 

  endif

! Update accumulations of SEA cells

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

     !$omp parallel do private(iw,kw,dta,wind,nls,fldval)
     do iws = 2,mws
        iw = itab_ws(iws)%iw         ! global index
        kw = itab_ws(iws)%kw

! If run is parallel, get local rank indices

        if (isubdomain == 1) then
           iw = itabg_w (iw )%iw_myrank
       endif

! Timestep for accumulating fluxes, DTA is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DTA.

        dta  = dtlm(itab_w(iw)%mrlw)
        wind = sqrt(vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2)

           vels_s_accum(iws) =    vels_s_accum(iws) + dta * real(wind,r8)
        airtemp_s_accum(iws) = airtemp_s_accum(iws) + dta * real(tair(kw,iw),r8)
         airshv_s_accum(iws) =  airshv_s_accum(iws) + dta * real(sh_v(kw,iw),r8)

        cantemp_s_accum(iws) = cantemp_s_accum(iws) + dta * real(sea%cantemp(iws),r8)
         canshv_s_accum(iws) =  canshv_s_accum(iws) + dta * real(sea%canshv (iws),r8)
         sfluxt_s_accum(iws) =  sfluxt_s_accum(iws) + dta * real(sea%sfluxt (iws),r8)
         sfluxr_s_accum(iws) =  sfluxr_s_accum(iws) + dta * real(sea%sfluxr (iws),r8)

        nls = sea%nlev_seaice(iws)

        if (nls > 0) then
           fldval = (1.0 - sea%seaicec(iws)) * sea%seatc(iws) &
                  +        sea%seaicec(iws)  * sea%seaice_tempk(nls,iws)
        else
           fldval = sea%seatc(iws)
        endif

        skintemp_s_accum(iws) =  skintemp_s_accum(iws) + dta * real(fldval,r8)
     enddo
     !$omp end parallel do 

  endif

end subroutine flux_accum

!===============================================================================

End Module mem_flux_accum
