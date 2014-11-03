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
Module mem_turb

  real, allocatable, target :: tkep    (:,:)
  real, allocatable, target :: epsp    (:,:)
  real, allocatable, target :: hkm     (:,:)
  real, allocatable, target :: hkh     (:,:)
  real, allocatable, target :: vkm     (:,:)
  real, allocatable, target :: vkh     (:,:)
  real, allocatable, target :: sxfer_tk(:,:)
  real, allocatable, target :: sxfer_rk(:,:)

  real, allocatable, target :: vkm_sfc (:,:)
  real, allocatable, target :: sfluxw  (:)
  real, allocatable, target :: sfluxt  (:)
  real, allocatable, target :: sfluxr  (:)
  real, allocatable, target :: ustar   (:)
  real, allocatable, target :: wstar   (:)
  real, allocatable, target :: wtv0    (:)
  real, allocatable, target :: pblh    (:)

  integer, allocatable, target :: kpblh(:)

  real,    allocatable      :: frac_land (:)
  real,    allocatable      :: frac_urb  (:)
  real,    allocatable      :: frac_sfc(:,:)

  real, allocatable, target :: fthpbl(:,:)
  real, allocatable, target :: fqtpbl(:,:)

Contains

!===============================================================================

  subroutine alloc_turb(mza, mwa, nsw_max, idiffk, mrls)

    use misc_coms, only: rinit, nqparm
    implicit none

    integer, intent(in) :: mza, mwa, nsw_max, idiffk, mrls

!   Allocate arrays based on options (if necessary)
!   Initialize arrays to zero

    allocate (hkm(mza,mwa)) ; hkm = rinit
!   allocate (hkh(mza,mwa)) ; hkh = rinit
!   allocate (vkm(mza,mwa)) ; vkm = rinit
!   allocate (vkh(mza,mwa)) ; vkh = rinit

    allocate (sxfer_tk(nsw_max,mwa)) ; sxfer_tk = 0.0
    allocate (sxfer_rk(nsw_max,mwa)) ; sxfer_rk = 0.0

    allocate (vkm_sfc (nsw_max,mwa)) ; vkm_sfc  = rinit
    allocate (frac_sfc(nsw_max,mwa)) ; frac_sfc = rinit

    allocate (sfluxt  (mwa)) ; sfluxt   = rinit
    allocate (sfluxr  (mwa)) ; sfluxr   = rinit
!   allocate (sfluxw  (mwa)) ; sfluxw   = rinit

    allocate (ustar    (mwa)) ; ustar     = rinit
    allocate (wstar    (mwa)) ; wstar     = rinit
    allocate (wtv0     (mwa)) ; wtv0      = rinit
    allocate (pblh     (mwa)) ; pblh      = rinit
    allocate (kpblh    (mwa)) ; kpblh     = 1
    allocate (frac_urb (mwa)) ; frac_urb  = 0.0 
    allocate (frac_land(mwa)) ; frac_land = 0.0

    if ( any(nqparm(1:mrls) == 1) .or. any(nqparm(1:mrls) == 2) .or. &
         any(nqparm(1:mrls) == 5) ) then
       allocate (fthpbl(mza,mwa)) ; fthpbl = 0.0
       allocate (fqtpbl(mza,mwa)) ; fqtpbl = 0.0
    endif

  end subroutine alloc_turb
  
  !===============================================================================

  subroutine dealloc_turb()

    implicit none

    if (allocated(tkep))    deallocate (tkep)
    if (allocated(epsp))    deallocate (epsp)
    if (allocated(hkm))     deallocate (hkm)
    if (allocated(hkh))     deallocate (hkh)
    if (allocated(vkm))     deallocate (vkm)
    if (allocated(vkh))     deallocate (vkh)
    if (allocated(vkm_sfc)) deallocate (vkm_sfc)
    if (allocated(sfluxw))  deallocate (sfluxw)
    if (allocated(sfluxt))  deallocate (sfluxt)
    if (allocated(sfluxr))  deallocate (sfluxr)
    if (allocated(ustar))   deallocate (ustar)
    if (allocated(wstar))   deallocate (wstar)
    if (allocated(wtv0))    deallocate (wtv0)
    if (allocated(pblh))    deallocate (pblh)
    if (allocated(kpblh))   deallocate (kpblh)
    if (allocated(fthpbl))  deallocate (fthpbl)
    if (allocated(fqtpbl))  deallocate (fqtpbl)

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

    if (allocated(hkh)) then
       call increment_vtable('HKH', 'AW')
       vtab_r(num_var)%rvar2_p => hkh
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
       vtab_r(num_var)%rvar2_p => vkm_sfc
    endif

    if (allocated(sfluxw)) then
       call increment_vtable('SFLUXW', 'AW')
       vtab_r(num_var)%rvar1_p => sfluxw
    endif

    if (allocated(sfluxt)) then
       call increment_vtable('SFLUXT', 'AW')
       vtab_r(num_var)%rvar1_p => sfluxt
    endif

    if (allocated(sfluxr)) then
       call increment_vtable('SFLUXR', 'AW')
       vtab_r(num_var)%rvar1_p => sfluxr
    endif

    if (allocated(ustar)) then
       call increment_vtable('USTAR', 'AW')
       vtab_r(num_var)%rvar1_p => ustar
    endif

    if (allocated(wstar)) then
       call increment_vtable('WSTAR', 'AW')
       vtab_r(num_var)%rvar1_p => wstar
    endif

    if (allocated(wtv0)) then
       call increment_vtable('WTV0', 'AW')
       vtab_r(num_var)%rvar1_p => wtv0
    endif

    if (allocated(pblh)) then
       call increment_vtable('PBLH', 'AW')
       vtab_r(num_var)%rvar1_p => pblh
    endif
    
    if (allocated(kpblh)) then
       call increment_vtable('KPBLH', 'AW')
       vtab_r(num_var)%ivar1_p => kpblh
    endif

    if (allocated(fthpbl)) then
       call increment_vtable('FTHPBL', 'AW')
       vtab_r(num_var)%rvar2_p => fthpbl
    endif

    if (allocated(fqtpbl)) then
       call increment_vtable('FQTPBL', 'AW')
       vtab_r(num_var)%rvar2_p => fqtpbl
    endif

  end subroutine filltab_turb

End Module mem_turb
