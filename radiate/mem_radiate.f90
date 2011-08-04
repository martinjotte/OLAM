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

Module mem_radiate

  integer :: jday ! Julian day

  real :: solfac  ! solar-constant coefficient for variable Earth-Sun dist
  real :: sunx    ! x-component of unit vector pointing to sun [m]
  real :: suny    ! y-component of unit vector pointing to sun [m]
  real :: sunz    ! z-component of unit vector pointing to sun [m]

  real, allocatable, target :: fthrd       (:,:)
  real, allocatable, target :: fthrd_lw    (:,:)

  real, allocatable, target :: rshort        (:)
  real, allocatable, target :: rlong         (:)
  real, allocatable, target :: rlongup       (:)
  real, allocatable, target :: rshort_top    (:)
  real, allocatable, target :: rshortup_top  (:)
  real, allocatable, target :: rlongup_top   (:)
  real, allocatable, target :: albedt        (:)
  real, allocatable         :: albedt_beam   (:)
  real, allocatable         :: albedt_diffuse(:)
  real, allocatable, target :: cosz          (:)
  real, allocatable         :: rlong_albedo  (:)
  real, allocatable         :: rshort_diffuse(:)

  integer, allocatable      :: rad_region    (:)

! These are used for adding extra levels at the top with the Mclatchy soundings

  integer, parameter :: maxadd_rad = 10 ! max allowed # of added rad levels
  integer            :: nadd_rad        ! actual # of added radiation levels
  real               :: zmrad = 30.e3   ! top of radiation grid

Contains

!===============================================================================

  subroutine alloc_radiate(mza,mwa,ilwrtyp,iswrtyp)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza
    integer, intent(in) :: mwa
    integer, intent(in) :: ilwrtyp
    integer, intent(in) :: iswrtyp

!   Allocate arrays based on options (if necessary)

    if (ilwrtyp + iswrtyp > 0)  then

       write(io6,*) 'allocating rad ', mwa, mza

       allocate (fthrd     (mza,mwa)) ; fthrd          = rinit
       allocate (fthrd_lw  (mza,mwa)) ; fthrd_lw       = rinit
       allocate (rshort        (mwa)) ; rshort         = rinit
       allocate (rlong         (mwa)) ; rlong          = rinit
       allocate (rlongup       (mwa)) ; rlongup        = rinit
       allocate (rshort_top    (mwa)) ; rshort_top     = rinit
       allocate (rshortup_top  (mwa)) ; rshortup_top   = rinit
       allocate (rlongup_top   (mwa)) ; rlongup_top    = rinit
       allocate (rshort_diffuse(mwa)) ; rshort_diffuse = rinit
       allocate (rlong_albedo  (mwa)) ; rlong_albedo   = rinit
       allocate (albedt        (mwa)) ; albedt         = rinit
       allocate (albedt_beam   (mwa)) ; albedt_beam    = rinit
       allocate (albedt_diffuse(mwa)) ; albedt_diffuse = rinit
       allocate (cosz          (mwa)) ; cosz           = rinit

       if (iswrtyp == 3 .or. ilwrtyp == 3) then
          allocate (rad_region(mwa))  ; rad_region = 0
       endif

    endif

  end subroutine alloc_radiate

!===============================================================================

  subroutine dealloc_radiate()

    implicit none

    if (allocated(fthrd))          deallocate (fthrd)
    if (allocated(fthrd_lw))       deallocate (fthrd_lw)
    if (allocated(rshort))         deallocate (rshort)
    if (allocated(rlong))          deallocate (rlong)
    if (allocated(rlongup))        deallocate (rlongup)
    if (allocated(rshort_top))     deallocate (rshort_top)
    if (allocated(rshortup_top))   deallocate (rshortup_top)
    if (allocated(rlongup_top))    deallocate (rlongup_top)
    if (allocated(rad_region))     deallocate (rad_region)
    if (allocated(rshort_diffuse)) deallocate (rshort_diffuse)
    if (allocated(rlong_albedo))   deallocate (rlong_albedo)
    if (allocated(albedt))         deallocate (albedt)
    if (allocated(albedt_beam))    deallocate (albedt_beam)
    if (allocated(albedt_diffuse)) deallocate (albedt_diffuse)
    if (allocated(cosz))           deallocate (cosz)

  end subroutine dealloc_radiate

!===============================================================================

  subroutine filltab_radiate()

    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(fthrd)) then
       call increment_vtable('FTHRD', 'AW')
       vtab_r(num_var)%rvar2_p => fthrd
    endif

    if (allocated(fthrd_lw)) then
       call increment_vtable('FTHRD_LW', 'AW')
       vtab_r(num_var)%rvar2_p => fthrd_lw
    endif

    if (allocated(rshort)) then
       call increment_vtable('RSHORT', 'AW')
       vtab_r(num_var)%rvar1_p => rshort
    endif

    if (allocated(rlong)) then
       call increment_vtable('RLONG', 'AW')
       vtab_r(num_var)%rvar1_p => rlong
    endif

    if (allocated(rlongup)) then
       call increment_vtable('RLONGUP', 'AW')
       vtab_r(num_var)%rvar1_p => rlongup
    endif

    if (allocated(rshort_top)) then
       call increment_vtable('RSHORT_TOP', 'AW')
        vtab_r(num_var)%rvar1_p => rshort_top
    endif

    if (allocated(rshortup_top)) then
       call increment_vtable('RSHORTUP_TOP', 'AW')
        vtab_r(num_var)%rvar1_p => rshortup_top
    endif

    if (allocated(rlongup_top)) then
       call increment_vtable('RLONGUP_TOP', 'AW')
        vtab_r(num_var)%rvar1_p => rlongup_top
    endif

    if (allocated(albedt)) then
       call increment_vtable('ALBEDT', 'AW')
        vtab_r(num_var)%rvar1_p => albedt
    endif

    if (allocated(cosz)) then
       call increment_vtable('COSZ', 'AW')
        vtab_r(num_var)%rvar1_p => cosz
    endif

   end subroutine filltab_radiate

End Module mem_radiate
