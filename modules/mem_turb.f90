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
Module mem_turb

  real, allocatable, target :: tkep    (:,:)
  real, allocatable, target :: epsp    (:,:)
  real, allocatable, target :: hkm     (:,:)
  real, allocatable, target :: vkm     (:,:)
  real, allocatable, target :: vkh     (:,:)
  real, allocatable, target :: sxfer_tk(:,:)
  real, allocatable, target :: sxfer_rk(:,:)

  real, allocatable, target :: vkm_sfc (:)
  real, allocatable, target :: sflux_w (:)
  real, allocatable, target :: sflux_t (:)
  real, allocatable, target :: sflux_r (:)
  real, allocatable, target :: vels    (:)
  real, allocatable, target :: ustar   (:)

Contains

!===============================================================================

  subroutine alloc_turb(mza,mwa,nsw_max,idiffk)

    use misc_coms, only: rinit
    implicit none

    integer, intent(in) :: mza,mwa,nsw_max,idiffk

!   Allocate arrays based on options (if necessary)
!   Initialize arrays to zero

!   CHECK FOR IDIFFK BASED ON GRID 1 FOR NOW

    if (idiffk == 1 .or. idiffk == 4 .or. idiffk == 5 .or. idiffk == 6) then
       allocate (tkep(mza,mwa)) ; tkep = rinit
    endif

    if (idiffk == 6) then
       allocate (epsp(mza,mwa)) ; epsp = rinit
    endif

    allocate (hkm(mza,mwa)) ; hkm = rinit
    allocate (vkm(mza,mwa)) ; vkm = rinit
    allocate (vkh(mza,mwa)) ; vkh = rinit

    allocate (sxfer_tk(nsw_max,mwa)) ; sxfer_tk = rinit
    allocate (sxfer_rk(nsw_max,mwa)) ; sxfer_rk = rinit

    allocate (vkm_sfc(mwa)) ; vkm_sfc = rinit
    allocate (sflux_t(mwa)) ; sflux_t = rinit
    allocate (sflux_r(mwa)) ; sflux_r = rinit
    allocate (ustar  (mwa)) ; ustar   = rinit
    allocate (vels   (mwa)) ; vels    = rinit
    allocate (sflux_w(mwa)) ; sflux_w = rinit

  end subroutine alloc_turb
  
  !===============================================================================

  subroutine dealloc_turb()

    implicit none

    if (allocated(tkep))    deallocate (tkep)
    if (allocated(epsp))    deallocate (epsp)
    if (allocated(hkm))     deallocate (hkm)
    if (allocated(vkm))     deallocate (vkm)
    if (allocated(vkh))     deallocate (vkh)
    if (allocated(vkm_sfc)) deallocate (vkm_sfc)
    if (allocated(sflux_w)) deallocate (sflux_w)
    if (allocated(sflux_t)) deallocate (sflux_t)
    if (allocated(sflux_r)) deallocate (sflux_r)
    if (allocated(vels))    deallocate (vels)
    if (allocated(ustar))   deallocate (ustar)

  end subroutine dealloc_turb

  !===============================================================================

  subroutine filltab_turb()

    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(tkep)) then
       call increment_vtable('TKEP', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => tkep
    endif

    if (allocated(epsp)) then
       call increment_vtable('EPSP', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => epsp
    endif

    if (allocated(hkm)) then
       call increment_vtable('HKM', 'AW')
       vtab_r(num_var)%rvar2_p => hkm
    endif

    if (allocated(vkm)) then
       call increment_vtable('VKM', 'AW')
       vtab_r(num_var)%rvar2_p => vkm
    endif

    if (allocated(vkh)) then
       call increment_vtable('VKH', 'AW')
       vtab_r(num_var)%rvar2_p => vkh
    endif

    if (allocated(sxfer_tk)) then
       call increment_vtable('SXFER_TK', 'AW')
       vtab_r(num_var)%rvar2_p => sxfer_tk
    endif

    if (allocated(sxfer_rk)) then
       call increment_vtable('SXFER_RK', 'AW')
       vtab_r(num_var)%rvar2_p => sxfer_rk
    endif

    if (allocated(vkm_sfc)) then
       call increment_vtable('VKM_SFC', 'AW')
       vtab_r(num_var)%rvar1_p => vkm_sfc
    endif

    if (allocated(sflux_w)) then
       call increment_vtable('SFLUX_W', 'AW')
       vtab_r(num_var)%rvar1_p => sflux_w
    endif

    if (allocated(sflux_t)) then
       call increment_vtable('SFLUX_T', 'AW')
       vtab_r(num_var)%rvar1_p => sflux_t
    endif

    if (allocated(sflux_r)) then
       call increment_vtable('SFLUX_R', 'AW')
       vtab_r(num_var)%rvar1_p => sflux_r
    endif

    if (allocated(vels)) then
       call increment_vtable('VELS', 'AW')
       vtab_r(num_var)%rvar1_p => vels
    endif

    if (allocated(ustar)) then
       call increment_vtable('USTAR', 'AW')
       vtab_r(num_var)%rvar1_p => ustar
    endif

  end subroutine filltab_turb

End Module mem_turb
