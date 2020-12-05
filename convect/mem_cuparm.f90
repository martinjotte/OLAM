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

Module mem_cuparm

   use consts_coms, only: r8

   real,    allocatable :: thsrc (:,:) ! heat / cp  tend ( rho X temp )
   real,    allocatable :: rtsrc (:,:) ! water mass tend ( rho X rr_w )
   real,    allocatable :: umsrc (:,:) ! rho * ue tend
   real,    allocatable :: vmsrc (:,:) ! rho * ve tend
   real,    allocatable :: qwcon (:,:) ! convective cloud water
   real(r8),allocatable :: aconpr  (:)
   real,    allocatable :: conprr  (:)
   real,    allocatable :: cbmf    (:) ! updraft mass flux
   integer, allocatable :: kcutop  (:)
   integer, allocatable :: kcubot  (:)
   integer, allocatable :: iactcu  (:)

   ! There are only needed for scalar mixing or subgrid CMAQ aqueous
   ! chemistry with the Grell scheme:
   real,    allocatable :: cddf    (:)
   integer, allocatable :: kudbot  (:)
   integer, allocatable :: kddtop  (:)
   integer, allocatable :: kddmax  (:)
   integer, allocatable :: kstabi  (:)

   ! Convective precipitation flux and evaporation rate; only
   ! needed with Grell deep convection and CMAQ chemistry
   real,    allocatable :: cu_pwa(:,:)
   real,    allocatable :: cu_pev(:,:)

   private :: r8

Contains

!===============================================================================

  subroutine alloc_cuparm(mza, mwa, mrls, nqparm)

    use consts_coms, only: r8
    use oname_coms,  only: nl

    implicit none

    integer, intent(in) :: mza, mwa, mrls
    integer, intent(in) :: nqparm(:)

    if ( any(nqparm(1:mrls) > 0) ) then

       ! Base tendency arrays for all deep convective schemes

       allocate (thsrc(mza,mwa)) ; thsrc  = 0.0
       allocate (rtsrc(mza,mwa)) ; rtsrc  = 0.0
       allocate (aconpr   (mwa)) ; aconpr = 0.0_r8
       allocate (conprr   (mwa)) ; conprr = 0.0
       allocate (qwcon(mza,mwa)) ; qwcon  = 0.0
       allocate (iactcu   (mwa)) ; iactcu = 0
       allocate (cbmf     (mwa)) ; cbmf   = 0.0
       allocate (kcutop   (mwa)) ; kcutop = -1
       allocate (kcubot   (mwa)) ; kcubot = -1

       ! Momentum tendencies for schemes that mix wind
       ! (currently Tiedtke, Grell deep, and Emanuel)

       if (nl%conv_uv_mix > 0) then
          if ( any(nqparm(1:mrls) == 1) .or. &
               any(nqparm(1:mrls) == 2) .or. &
               any(nqparm(1:mrls) == 4) ) then

             allocate (umsrc(mza,mwa)) ; umsrc = 0.0
             allocate (vmsrc(mza,mwa)) ; vmsrc = 0.0
          endif
       endif

       ! Do we need to save the convective precipitation flux?
       ! Only needed with CMAQ and Grell deep convection.

       if (nl%do_chem > 0) then
          if ( any(nqparm(1:mrls) == 2) ) then
             allocate (cu_pwa(mza,mwa)) ; cu_pwa = 0.0
             allocate (cu_pev(mza,mwa)) ; cu_pev = 0.0
          endif
       endif

       ! Arrays for tracer mixing and/or aqueous chemistry with
       ! the Grell scheme (deep and shallow)

       if ( any(nqparm(1:mrls) == 2) .or. &
            any(nqparm(1:mrls) == 5) ) then

          allocate(cddf  (mwa)) ; cddf   = 0.0
          allocate(kudbot(mwa)) ; kudbot = -1
          allocate(kddtop(mwa)) ; kddtop = -1
          allocate(kddmax(mwa)) ; kddmax = -1
          allocate(kstabi(mwa)) ; kstabi = -1
       endif

  end subroutine alloc_cuparm

!===============================================================================

  subroutine dealloc_cuparm()
    implicit none

    if (allocated(thsrc))    deallocate (thsrc)
    if (allocated(rtsrc))    deallocate (rtsrc)
    if (allocated(umsrc))    deallocate (umsrc)
    if (allocated(vmsrc))    deallocate (vmsrc)
    if (allocated(aconpr))   deallocate (aconpr)
    if (allocated(conprr))   deallocate (conprr)
    if (allocated(qwcon))    deallocate (qwcon)
    if (allocated(cbmf))     deallocate (cbmf)
    if (allocated(cddf))     deallocate (cddf)
    if (allocated(kcutop))   deallocate (kcutop)
    if (allocated(kcubot))   deallocate (kcubot)
    if (allocated(kudbot))   deallocate (kudbot)
    if (allocated(kddtop))   deallocate (kddtop)
    if (allocated(kddmax))   deallocate (kddmax)
    if (allocated(iactcu))   deallocate (iactcu)
    if (allocated(cu_pwa))   deallocate (cu_pwa)
    if (allocated(cu_pev))   deallocate (cu_pev)

  end subroutine dealloc_cuparm

!===============================================================================

  subroutine filltab_cuparm()

    use var_tables, only: increment_vtable
    implicit none

     if (allocated(thsrc))  call increment_vtable('THSRC', 'AW', rvar2=thsrc)

     if (allocated(rtsrc))  call increment_vtable('RTSRC', 'AW', rvar2=rtsrc)

     if (allocated(umsrc))  call increment_vtable('UMSRC', 'AW', rvar2=umsrc)

     if (allocated(vmsrc))  call increment_vtable('VMSRC', 'AW', rvar2=vmsrc)

     if (allocated(aconpr)) call increment_vtable('ACONPR','AW', dvar1=aconpr)

     if (allocated(conprr)) call increment_vtable('CONPRR','AW', rvar1=conprr)

     if (allocated(qwcon))  call increment_vtable('QWCON', 'AW', rvar2=qwcon)

     if (allocated(cbmf))   call increment_vtable('CBMF',  'AW', rvar1=cbmf)

     if (allocated(cddf))   call increment_vtable('CDDF',  'AW', rvar1=cddf)

     if (allocated(kcutop)) call increment_vtable('KCUTOP','AW', ivar1=kcutop)

     if (allocated(kcubot)) call increment_vtable('KCUBOT','AW', ivar1=kcubot)

     if (allocated(kudbot)) call increment_vtable('KUDBOT','AW', ivar1=kudbot)

     if (allocated(kddtop)) call increment_vtable('KDDTOP','AW', ivar1=kddtop)

     if (allocated(kddmax)) call increment_vtable('KDDMAX','AW', ivar1=kddmax)

     if (allocated(kstabi)) call increment_vtable('KSTABI','AW', ivar1=kstabi)

     if (allocated(iactcu)) call increment_vtable('IACTCU','AW', ivar1=iactcu)

     if (allocated(cu_pwa)) call increment_vtable('CU_PWA','AW', rvar2=cu_pwa)

     if (allocated(cu_pev)) call increment_vtable('CU_PEV','AW', rvar2=cu_pev)

   end subroutine filltab_cuparm

End Module mem_cuparm
