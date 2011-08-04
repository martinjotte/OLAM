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

Module mem_timeavg

  real, allocatable, target :: rshort_avg      (:)
  real, allocatable, target :: rshortup_avg    (:)
  real, allocatable, target :: rlong_avg       (:)
  real, allocatable, target :: rlongup_avg     (:)
  real, allocatable, target :: rshort_top_avg  (:)
  real, allocatable, target :: rshortup_top_avg(:)
  real, allocatable, target :: rlongup_top_avg (:)
  real, allocatable, target :: sflux_t_avg     (:)
  real, allocatable, target :: sflux_r_avg     (:)

Contains

!===============================================================================

subroutine alloc_timeavg(mza,mwa)

  use misc_coms, only: io6
   
  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: mwa

! Allocate arrays for time-averaged quantities
! Initialize arrays to zero

  allocate (rshort_avg      (mwa)) ; rshort_avg       = 0.
  allocate (rshortup_avg    (mwa)) ; rshortup_avg     = 0.
  allocate (rlong_avg       (mwa)) ; rlong_avg        = 0.
  allocate (rlongup_avg     (mwa)) ; rlongup_avg      = 0.
  allocate (rshort_top_avg  (mwa)) ; rshort_top_avg   = 0.
  allocate (rshortup_top_avg(mwa)) ; rshortup_top_avg = 0.
  allocate (rlongup_top_avg (mwa)) ; rlongup_top_avg  = 0.
  allocate (sflux_t_avg     (mwa)) ; sflux_t_avg      = 0.
  allocate (sflux_r_avg     (mwa)) ; sflux_r_avg      = 0.
  
end subroutine alloc_timeavg

!===============================================================================

subroutine dealloc_timeavg()

  implicit none

  if (allocated (rshort_avg      )) deallocate (rshort_avg      )
  if (allocated (rshortup_avg    )) deallocate (rshortup_avg    )
  if (allocated (rlong_avg       )) deallocate (rlong_avg       )
  if (allocated (rlongup_avg     )) deallocate (rlongup_avg     )
  if (allocated (rshort_top_avg  )) deallocate (rshort_top_avg  )
  if (allocated (rshortup_top_avg)) deallocate (rshortup_top_avg)
  if (allocated (rlongup_top_avg )) deallocate (rlongup_top_avg )
  if (allocated (sflux_t_avg     )) deallocate (sflux_t_avg     )
  if (allocated (sflux_r_avg     )) deallocate (sflux_r_avg     )

end subroutine dealloc_timeavg

!===============================================================================

subroutine filltab_timeavg()

  use var_tables, only: vtab_r, num_var, increment_vtable
  implicit none

  if (allocated (rshort_avg)) then
     call increment_vtable('RSHORT_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rshort_avg
  endif

  if (allocated(rshortup_avg)) then
     call increment_vtable('RSHORTUP_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rshortup_avg
  endif

  if (allocated(rlong_avg)) then
     call increment_vtable('RLONG_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rlong_avg
  endif

  if (allocated(rlongup_avg)) then
     call increment_vtable('RLONGUP_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rlongup_avg
  endif

  if (allocated(rshort_top_avg)) then
     call increment_vtable('RSHORT_TOP_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rshort_top_avg
  endif

  if (allocated(rshortup_top_avg)) then
     call increment_vtable('RSHORTUP_TOP_AVG','AW')
     vtab_r(num_var)%rvar1_p => rshortup_top_avg
  endif

  if (allocated(rlongup_top_avg)) then
     call increment_vtable('RLONGUP_TOP_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => rlongup_top_avg
  endif

  if (allocated(sflux_t_avg)) then
     call increment_vtable('SFLUX_T_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => sflux_t_avg
  endif

  if (allocated(sflux_r_avg)) then
     call increment_vtable('SFLUX_R_AVG', 'AW')
     vtab_r(num_var)%rvar1_p => sflux_r_avg
  endif

end subroutine filltab_timeavg

!===============================================================================

subroutine accum_timeavg()

  use misc_coms,   only: io6, time_istp8, dtlm, time_prevhist

  use mem_ijtabs,  only: istp, itab_w, jtab_w, mrl_begl

  use mem_radiate, only: albedt, rshort, rlong, rlongup, &
                         rshort_top, rshortup_top, rlongup_top

  use leaf_coms,   only: mrl_leaf, dt_leaf

  use mem_turb,    only: sflux_t, sflux_r

!$ use omp_lib

  implicit none

  integer :: mrl, j, iw
  real    :: dt_avg, wt_new, wt_old

! Update averages of ATM SFC turbulent fluxes

  call psub()
!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0) then

     if (mrl <= mrl_leaf) mrl = 1  ! special mrl set for leaf

!$omp parallel do private(iw,dt_avg,wt_new,wt_old)
     do j = 1,jtab_w(20)%jend(mrl); iw = jtab_w(20)%iw(j)
!----------------------------------------------------------------------
        call qsub('W',iw)

! Timestep for computing averages, DT_AVG, is the lesser of dtlm for the given
! IW cell and dt_leaf; the frequency that each IW cell is processed in this
! loop should be consistent with (inversely proportional to) that DT_AVG.

        dt_avg = min(dtlm(itab_w(iw)%mrlw),dt_leaf)

! Compute averaging weights based on DT_AVG and elapsed time since the
! previous history write.  This keeps the values in the time-averaged arrays
! equal to true averages (over time so far elapsed), allowing them to be
! correctly plotted at any time.  It also avoids the need to re-zero each
! average after each history write.  The following formulation assumes that
! this subroutine is called from subroutine timestep inside the sub-cycle loop
! AFTER any new flux update is peformed and BEFORE time_istp8 is updated.

        wt_new = dt_avg / (dt_avg + time_istp8 - time_prevhist)
        wt_old = 1. - wt_new

! Modify averages 

        sflux_t_avg(iw) = wt_old * sflux_t_avg(iw) &
             + wt_new * sflux_t(iw)

        sflux_r_avg(iw) = wt_old * sflux_r_avg(iw) &
             + wt_new * sflux_r(iw)

     enddo
!$omp end parallel do

  endif
  call rsub('Waf',20)

! Update averages of ATM SFC and TOP radiative fluxes

  call psub()
!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0) then
!$omp parallel do private(iw,dt_avg,wt_new,wt_old)
     do j = 1,jtab_w(12)%jend(mrl); iw = jtab_w(12)%iw(j)
!----------------------------------------------------------------------
        call qsub('W',iw)

! Timestep for computing averages, DT_AVG is dtlm for the given IW cell; 
! the frequency that each IW cell is processed in this loop should be 
! consistent with (inversely proportional to) that DT_AVG.

        dt_avg = dtlm(itab_w(iw)%mrlw)

! Compute averaging weights based on DT_AVG and elapsed time since the
! previous history write.  This keeps the values in the time-averaged arrays
! equal to true averages (over time so far elapsed), allowing them to be
! correctly plotted at any time.  It also avoids the need to re-zero each
! average after each history write.  The following formulation assumes that
! this subroutine is called from subroutine timestep inside the sub-cycle loop
! AFTER any new flux update is peformed and BEFORE time_istp8 is updated.

        wt_new = dt_avg / (dt_avg + time_istp8 - time_prevhist)
        wt_old = 1. - wt_new

! Modify averages 

        rshort_avg(iw) = wt_old * rshort_avg(iw) &
             + wt_new * rshort(iw)

        rshortup_avg(iw) = wt_old * rshortup_avg(iw) &
             + wt_new * rshort(iw) * albedt(iw)

        rlong_avg(iw) = wt_old * rlong_avg(iw) &
             + wt_new * rlong(iw)

        rlongup_avg(iw) = wt_old * rlongup_avg(iw) &
             + wt_new * rlongup(iw)

        rshort_top_avg(iw) = wt_old * rshort_top_avg(iw) &
             + wt_new * rshort_top(iw)

!       if (mod(iw,1000) == 1) then
!          print*, ' '
!          print*, 'mta1 ',iw,dt_avg,wt_old,wt_new
!          print*, 'mta2 ',rshort_top(iw),rshort_top_avg(iw)
!       endif

        rshortup_top_avg(iw) = wt_old * rshortup_top_avg(iw) &
             + wt_new * rshortup_top(iw)

        rlongup_top_avg(iw) = wt_old * rlongup_top_avg(iw) &
             + wt_new * rlongup_top(iw)

     enddo
!$omp end parallel do
  endif
  call rsub('Waf',12)

end subroutine accum_timeavg

!===============================================================================

End Module mem_timeavg
