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
  real(r8), allocatable, target :: sfluxt_accum      (:)
  real(r8), allocatable, target :: sfluxr_accum      (:)

! Land arrays

  real(r8), allocatable, target :: sfluxt_l_accum (:)
  real(r8), allocatable, target :: sfluxr_l_accum (:)
  real(r8), allocatable, target :: airtemp_l_accum(:)
  real(r8), allocatable, target :: airshv_l_accum (:)
  real(r8), allocatable, target :: cantemp_l_accum(:)
  real(r8), allocatable, target :: canshv_l_accum (:)

! Sea arrays

  real(r8), allocatable, target :: sfluxt_s_accum (:)
  real(r8), allocatable, target :: sfluxr_s_accum (:)
  real(r8), allocatable, target :: airtemp_s_accum(:)
  real(r8), allocatable, target :: airshv_s_accum (:)
  real(r8), allocatable, target :: cantemp_s_accum(:)
  real(r8), allocatable, target :: canshv_s_accum (:)


Contains

!===============================================================================

subroutine alloc_flux_accum(mza,mwa,mwl,mws)

  use misc_coms, only: ilwrtyp, iswrtyp
  use leaf_coms, only: isfcl
  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: mza, mwa, mwl, mws

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

  if (isfcl > 0) then
     allocate (sfluxt_accum   (mwa)) ; sfluxt_accum    = 0._r8
     allocate (sfluxr_accum   (mwa)) ; sfluxr_accum    = 0._r8

     allocate (sfluxt_l_accum (mwl)) ; sfluxt_l_accum  = 0._r8
     allocate (sfluxr_l_accum (mwl)) ; sfluxr_l_accum  = 0._r8
     allocate (airtemp_l_accum(mwl)) ; airtemp_l_accum = 0._r8
     allocate (airshv_l_accum (mwl)) ; airshv_l_accum  = 0._r8
     allocate (cantemp_l_accum(mwl)) ; cantemp_l_accum = 0._r8
     allocate (canshv_l_accum (mwl)) ; canshv_l_accum  = 0._r8

     allocate (sfluxt_s_accum (mws)) ; sfluxt_s_accum  = 0._r8
     allocate (sfluxr_s_accum (mws)) ; sfluxr_s_accum  = 0._r8
     allocate (airtemp_s_accum(mws)) ; airtemp_s_accum = 0._r8
     allocate (airshv_s_accum (mws)) ; airshv_s_accum  = 0._r8
     allocate (cantemp_s_accum(mws)) ; cantemp_s_accum = 0._r8
     allocate (canshv_s_accum (mws)) ; canshv_s_accum  = 0._r8
  endif
  
end subroutine alloc_flux_accum

!===============================================================================

subroutine filltab_flux_accum()

  use var_tables, only: vtab_r, num_var, increment_vtable
  implicit none

  if (allocated (rshort_accum)) then
     call increment_vtable('RSHORT_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rshort_accum
  endif

  if (allocated(rshortup_accum)) then
     call increment_vtable('RSHORTUP_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rshortup_accum
  endif

  if (allocated(rlong_accum)) then
     call increment_vtable('RLONG_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rlong_accum
  endif

  if (allocated(rlongup_accum)) then
     call increment_vtable('RLONGUP_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rlongup_accum
  endif

  if (allocated(rshort_top_accum)) then
     call increment_vtable('RSHORT_TOP_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rshort_top_accum
  endif

  if (allocated(rshortup_top_accum)) then
     call increment_vtable('RSHORTUP_TOP_ACCUM','AW')
     vtab_r(num_var)%dvar1_p => rshortup_top_accum
  endif

  if (allocated(rlongup_top_accum)) then
     call increment_vtable('RLONGUP_TOP_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => rlongup_top_accum
  endif

  if (allocated(sfluxt_accum)) then
     call increment_vtable('SFLUXT_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => sfluxt_accum
  endif

  if (allocated(sfluxr_accum)) then
     call increment_vtable('SFLUXR_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => sfluxr_accum
  endif

  if (allocated(sfluxt_l_accum)) then
     call increment_vtable('SFLUXT_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => sfluxt_l_accum
  endif

  if (allocated(sfluxr_l_accum)) then
     call increment_vtable('SFLUXR_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => sfluxr_l_accum
  endif

  if (allocated(airtemp_l_accum)) then
     call increment_vtable('AIRTEMP_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => airtemp_l_accum
  endif

  if (allocated(airshv_l_accum)) then
     call increment_vtable('AIRSHV_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => airshv_l_accum
  endif

  if (allocated(cantemp_l_accum)) then
     call increment_vtable('CANTEMP_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => cantemp_l_accum
  endif

  if (allocated(canshv_l_accum)) then
     call increment_vtable('CANSHV_L_ACCUM', 'LW')
     vtab_r(num_var)%dvar1_p => canshv_l_accum
  endif

  if (allocated(sfluxt_s_accum)) then
     call increment_vtable('SFLUXT_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => sfluxt_s_accum
  endif

  if (allocated(sfluxr_s_accum)) then
     call increment_vtable('SFLUXR_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => sfluxr_s_accum
  endif

  if (allocated(airtemp_s_accum)) then
     call increment_vtable('AIRTEMP_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => airtemp_s_accum
  endif

  if (allocated(airshv_s_accum)) then
     call increment_vtable('AIRSHV_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => airshv_s_accum
  endif

  if (allocated(cantemp_s_accum)) then
     call increment_vtable('CANTEMP_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => cantemp_s_accum
  endif

  if (allocated(canshv_s_accum)) then
     call increment_vtable('CANSHV_S_ACCUM', 'SW')
     vtab_r(num_var)%dvar1_p => canshv_s_accum
  endif

end subroutine filltab_flux_accum

!===============================================================================

subroutine flux_accum()

  use misc_coms,   only: io6, time_istp8, dtlm, time_prevhist, &
                         ilwrtyp, iswrtyp, isubdomain

  use mem_ijtabs,  only: istp, itab_w, itabg_w, jtab_w, mrl_begl, jtw_prog

  use mem_radiate, only: albedt, rshort, rlong, rlongup, &
                         rshort_top, rshortup_top, rlongup_top

  use leaf_coms,   only: mwl, mrl_leaf, dt_leaf, isfcl

  use mem_turb,    only: sfluxt, sfluxr

  use consts_coms, only: r8

  use sea_coms,    only: mws
  use mem_sea,     only: sea, itabg_ws, itab_ws

  use mem_leaf,    only: land, itabg_wl, itab_wl

  implicit none

  integer :: mrl, j, iw, iwl, iws

  real(r8) :: dta

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

        dta = min(dtlm(itab_w(iw)%mrlw),dt_leaf)

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
        endif

        if (ilwrtyp > 0) then
           rlong_accum      (iw) = rlong_accum      (iw) + dta * real(rlong(iw),r8)
           rlongup_accum    (iw) = rlongup_accum    (iw) + dta * real(rlongup(iw),r8)
           rlongup_top_accum(iw) = rlongup_top_accum(iw) + dta * real(rlongup_top(iw),r8)
        endif

     enddo
!$omp end parallel do
  endif

! Update accumulations of LAND cells

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

!$omp parallel do private(iw,dta)
     do iwl = 2,mwl
        iw = itab_wl(iwl)%iw         ! global index

! If run is parallel, get local rank indices

        if (isubdomain == 1) then
           iw = itabg_w (iw )%iw_myrank
        endif

! Timestep for accumulating fluxes, DTA is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DTA.

        dta = dtlm(itab_w(iw)%mrlw)

         sfluxt_l_accum(iwl) =  sfluxt_l_accum(iwl) + dta * real(land%sfluxt  (iwl),r8)
         sfluxr_l_accum(iwl) =  sfluxr_l_accum(iwl) + dta * real(land%sfluxr  (iwl),r8)
        airtemp_l_accum(iwl) = airtemp_l_accum(iwl) + dta * real(land%airtemp (iwl),r8)
         airshv_l_accum(iwl) =  airshv_l_accum(iwl) + dta * real(land%airshv  (iwl),r8)
        cantemp_l_accum(iwl) = cantemp_l_accum(iwl) + dta * real(land%cantemp(iwl),r8)
         canshv_l_accum(iwl) =  canshv_l_accum(iwl) + dta * real(land%canshv (iwl),r8)

     enddo
  endif

! Update accumulations of SEA cells

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

!$omp parallel do private(iw,dta)
     do iws = 2,mws
        iw = itab_ws(iws)%iw         ! global index

! If run is parallel, get local rank indices

        if (isubdomain == 1) then
           iw = itabg_w (iw )%iw_myrank
       endif

! Timestep for accumulating fluxes, DTA is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DTA.

        dta = dtlm(itab_w(iw)%mrlw)

         sfluxt_s_accum(iws) =  sfluxt_s_accum(iws) + dta * real(sea%sfluxt  (iws),r8)
         sfluxr_s_accum(iws) =  sfluxr_s_accum(iws) + dta * real(sea%sfluxr  (iws),r8)
        airtemp_s_accum(iws) = airtemp_s_accum(iws) + dta * real(sea%airtemp (iws),r8)
         airshv_s_accum(iws) =  airshv_s_accum(iws) + dta * real(sea%airshv  (iws),r8)
        cantemp_s_accum(iws) = cantemp_s_accum(iws) + dta * real(sea%cantemp(iws),r8)
         canshv_s_accum(iws) =  canshv_s_accum(iws) + dta * real(sea%canshv (iws),r8)

     enddo
  endif

end subroutine flux_accum

!===============================================================================

End Module mem_flux_accum
