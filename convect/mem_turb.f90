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

  real,    allocatable, target :: tkep    (:,:)
  real,    allocatable, target :: epsp    (:,:)
  real,    allocatable, target :: vkm     (:,:)
  real,    allocatable, target :: vkh     (:,:)
  real,    allocatable, target :: sxfer_tk(:,:)
  real,    allocatable, target :: sxfer_rk(:,:)
  real,    allocatable, target :: vkm_sfc (:,:)

  real,    allocatable, target :: sfluxt(:)
  real,    allocatable, target :: sfluxr(:)
  real,    allocatable, target :: ustar (:)
  real,    allocatable, target :: wstar (:)
  real,    allocatable, target :: moli  (:)
  real,    allocatable, target :: wtv0  (:)
  real,    allocatable, target :: pblh  (:)
  integer, allocatable, target :: kpblh (:)

  real,    allocatable, target :: fthpbl(:,:)
  real,    allocatable, target :: fqtpbl(:,:)

  real,    allocatable         :: frac_land (:)
  real,    allocatable         :: frac_sea  (:)
  real,    allocatable         :: frac_lake (:)
  real,    allocatable         :: frac_urb  (:)
  real,    allocatable         :: frac_sfc(:,:)

  real,    allocatable         :: akmodx(:,:)
  real,    allocatable         :: akhodx(:,:)

Contains

!===============================================================================

  subroutine alloc_turb(mza, mwa, mva, nsw_max, idiffk, mrls)

    use misc_coms, only: rinit
    implicit none

    integer, intent(in) :: mza, mwa, mva, nsw_max, idiffk(:), mrls

!   Allocate arrays based on options (if necessary)
!   Initialize arrays to zero

    allocate (sxfer_tk(nsw_max,mwa)) ; sxfer_tk = 0.0
    allocate (sxfer_rk(nsw_max,mwa)) ; sxfer_rk = 0.0

    allocate (vkm_sfc (nsw_max,mwa)) ; vkm_sfc  = 0.0
    allocate (frac_sfc(nsw_max,mwa)) ; frac_sfc = 0.0

    allocate (vkm   (mza,mwa)) ; vkm       = rinit
    allocate (vkh   (mza,mwa)) ; vkh       = rinit
    allocate (fthpbl(mza,mwa)) ; fthpbl    = 0.0
    allocate (fqtpbl(mza,mwa)) ; fqtpbl    = 0.0

    allocate (sfluxt    (mwa)) ; sfluxt    = rinit
    allocate (sfluxr    (mwa)) ; sfluxr    = rinit
    allocate (ustar     (mwa)) ; ustar     = rinit
    allocate (wstar     (mwa)) ; wstar     = rinit
    allocate (moli      (mwa)) ; moli      = rinit
    allocate (wtv0      (mwa)) ; wtv0      = rinit
    allocate (pblh      (mwa)) ; pblh      = rinit
    allocate (kpblh     (mwa)) ; kpblh     = 1
    allocate (frac_urb  (mwa)) ; frac_urb  = 0.0
    allocate (frac_land (mwa)) ; frac_land = 0.0
    allocate (frac_lake (mwa)) ; frac_lake = 0.0
    allocate (frac_sea  (mwa)) ; frac_sea  = 0.0

    allocate (akmodx(mza,mva)) ; akmodx    = 0.0
    allocate (akhodx(mza,mva)) ; akhodx    = 0.0

  end subroutine alloc_turb
  
  !===============================================================================

  subroutine dealloc_turb()

    implicit none

    if (allocated(tkep))    deallocate (tkep)
    if (allocated(epsp))    deallocate (epsp)
    if (allocated(vkm))     deallocate (vkm)
    if (allocated(vkh))     deallocate (vkh)
    if (allocated(vkm_sfc)) deallocate (vkm_sfc)
    if (allocated(sfluxt))  deallocate (sfluxt)
    if (allocated(sfluxr))  deallocate (sfluxr)
    if (allocated(ustar))   deallocate (ustar)
    if (allocated(wstar))   deallocate (wstar)
    if (allocated(moli))    deallocate (moli)
    if (allocated(wtv0))    deallocate (wtv0)
    if (allocated(pblh))    deallocate (pblh)
    if (allocated(kpblh))   deallocate (kpblh)
    if (allocated(fthpbl))  deallocate (fthpbl)
    if (allocated(fqtpbl))  deallocate (fqtpbl)
    if (allocated(akmodx))  deallocate (akmodx)
    if (allocated(akhodx))  deallocate (akhodx)

  end subroutine dealloc_turb

!===============================================================================

  subroutine filltab_turb()

    use var_tables, only: increment_vtable
    implicit none

    if (allocated(tkep))     call increment_vtable('TKEP',    'AW', rvar2=tkep, mpt1=.true.)

    if (allocated(epsp))     call increment_vtable('EPSP',    'AW', rvar2=epsp, mpt1=.true.)

    if (allocated(vkm))      call increment_vtable('VKM',     'AW', rvar2=vkm)

    if (allocated(vkh))      call increment_vtable('VKH',     'AW', rvar2=vkh)

    if (allocated(sxfer_tk)) call increment_vtable('SXFER_TK','AW', rvar2=sxfer_tk)

    if (allocated(sxfer_rk)) call increment_vtable('SXFER_RK','AW', rvar2=sxfer_rk)

    if (allocated(vkm_sfc))  call increment_vtable('VKM_SFC', 'AW', rvar2=vkm_sfc)

    if (allocated(sfluxt))   call increment_vtable('SFLUXT',  'AW', rvar1=sfluxt)

    if (allocated(sfluxr))   call increment_vtable('SFLUXR',  'AW', rvar1=sfluxr)

    if (allocated(ustar))    call increment_vtable('USTAR',   'AW', rvar1=ustar)

    if (allocated(wstar))    call increment_vtable('WSTAR',   'AW', rvar1=wstar)

    if (allocated(wtv0))     call increment_vtable('WTV0',    'AW', rvar1=wtv0)

    if (allocated(pblh))     call increment_vtable('PBLH',    'AW', rvar1=pblh)
    
    if (allocated(kpblh))    call increment_vtable('KPBLH',   'AW', ivar1=kpblh)

    if (allocated(fthpbl))   call increment_vtable('FTHPBL',  'AW', rvar2=fthpbl)

    if (allocated(fqtpbl))   call increment_vtable('FQTPBL',  'AW', rvar2=fqtpbl)

  end subroutine filltab_turb

End Module mem_turb
