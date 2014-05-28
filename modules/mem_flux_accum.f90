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

Module mem_flux_accum

  use consts_coms, only: r8

  real(r8), allocatable, target :: rshort_accum      (:)
  real(r8), allocatable, target :: rshortup_accum    (:)
  real(r8), allocatable, target :: rlong_accum       (:)
  real(r8), allocatable, target :: rlongup_accum     (:)
  real(r8), allocatable, target :: rshort_top_accum  (:)
  real(r8), allocatable, target :: rshortup_top_accum(:)
  real(r8), allocatable, target :: rlongup_top_accum (:)
  real(r8), allocatable, target :: sflux_t_accum     (:)
  real(r8), allocatable, target :: sflux_r_accum     (:)

Contains

!===============================================================================

subroutine alloc_flux_accum(mza,mwa)

  use misc_coms, only: ilwrtyp, iswrtyp
  use leaf_coms, only: isfcl
  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: mwa

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
     allocate (sflux_t_accum     (mwa)) ; sflux_t_accum      = 0._r8
     allocate (sflux_r_accum     (mwa)) ; sflux_r_accum      = 0._r8
  endif
  
end subroutine alloc_flux_accum

!===============================================================================

subroutine dealloc_flux_accum()

  implicit none

  if (allocated (rshort_accum      )) deallocate (rshort_accum      )
  if (allocated (rshortup_accum    )) deallocate (rshortup_accum    )
  if (allocated (rlong_accum       )) deallocate (rlong_accum       )
  if (allocated (rlongup_accum     )) deallocate (rlongup_accum     )
  if (allocated (rshort_top_accum  )) deallocate (rshort_top_accum  )
  if (allocated (rshortup_top_accum)) deallocate (rshortup_top_accum)
  if (allocated (rlongup_top_accum )) deallocate (rlongup_top_accum )
  if (allocated (sflux_t_accum     )) deallocate (sflux_t_accum     )
  if (allocated (sflux_r_accum     )) deallocate (sflux_r_accum     )

end subroutine dealloc_flux_accum

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

  if (allocated(sflux_t_accum)) then
     call increment_vtable('SFLUX_T_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => sflux_t_accum
  endif

  if (allocated(sflux_r_accum)) then
     call increment_vtable('SFLUX_R_ACCUM', 'AW')
     vtab_r(num_var)%dvar1_p => sflux_r_accum
  endif

end subroutine filltab_flux_accum

!===============================================================================

subroutine flux_accum()

  use misc_coms,   only: io6, time_istp8, dtlm, time_prevhist, &
                         ilwrtyp, iswrtyp

  use mem_ijtabs,  only: istp, itab_w, jtab_w, mrl_begl, jtw_prog

  use mem_radiate, only: albedt, rshort, rlong, rlongup, &
                         rshort_top, rshortup_top, rlongup_top

  use leaf_coms,   only: mrl_leaf, dt_leaf, isfcl

  use mem_turb,    only: sflux_t, sflux_r

  use consts_coms, only: r8

!$ use omp_lib

  implicit none

  integer :: mrl, j, iw

  real(r8) :: dta

! Update accumulations of ATM SFC turbulent fluxes

  call psub()
!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. isfcl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

!$omp parallel do private(iw,dta)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
        call qsub('W',iw)

! Timestep for accumulating fluxes, DTA, is the lesser of dtlm for the given
! IW cell and dt_leaf; the frequency that each IW cell is processed in this
! loop should be consistent with (inversely proportional to) that DTA.

        dta = min(dtlm(itab_w(iw)%mrlw),dt_leaf)

        sflux_t_accum(iw) = sflux_t_accum(iw) + dta * real(sflux_t(iw),r8)
        sflux_r_accum(iw) = sflux_r_accum(iw) + dta * real(sflux_r(iw),r8)

     enddo
!$omp end parallel do

  endif
  call rsub('Waf',20)

! Update accumulations of ATM SFC and TOP radiative fluxes

  call psub()
!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0 .and. ilwrtyp + iswrtyp > 0) then
!$omp parallel do private(iw,dta)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
        call qsub('W',iw)

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
  call rsub('Waf',12)

end subroutine flux_accum

!===============================================================================

End Module mem_flux_accum
