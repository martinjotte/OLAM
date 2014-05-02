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

Module mem_cuparm

   real,    allocatable, target :: thsrc (:,:)
   real,    allocatable, target :: rtsrc (:,:)
   real,    allocatable, target :: aconpr  (:)
   real,    allocatable, target :: conprr  (:)
   real,    allocatable, target :: vxsrc (:,:)
   real,    allocatable, target :: vysrc (:,:)
   real,    allocatable, target :: vzsrc (:,:)

   real,    allocatable, target :: cbmf  (:)
   integer, allocatable, target :: kcutop(:)
   integer, allocatable, target :: kcubot(:)

Contains

!===============================================================================

  subroutine alloc_cuparm(mza, mwa, mrls, nqparm)

    use misc_coms,  only: rinit
    use oname_coms, only: nl
    implicit none

    integer, intent(in) :: mza, mwa, mrls
    integer, intent(in) :: nqparm(:)
   
    if ( any(nqparm(1:mrls) > 0) ) then      
       
       ! Base tendency arrays for all deep convective schemes
       
       allocate (thsrc(mza,mwa)) ; thsrc  = 0.0
       allocate (rtsrc(mza,mwa)) ; rtsrc  = 0.0
       allocate (aconpr   (mwa)) ; aconpr = 0.0
       allocate (conprr   (mwa)) ; conprr = 0.0

       ! Extra arrays for momentum mixing

       if (nl%conv_uv_mix > 0) then
          allocate (vxsrc(mza,mwa)) ; vxsrc = 0.0
          allocate (vysrc(mza,mwa)) ; vysrc = 0.0
          allocate (vzsrc(mza,mwa)) ; vzsrc = 0.0
       endif

       ! Diagnostic arrays for clouds/radiation/tracer mixing

       allocate(cbmf  (mwa)) ; cbmf   = 0.0
       allocate(kcutop(mwa)) ; kcutop = -1
       allocate(kcubot(mwa)) ; kcubot = -1

    endif
       
  end subroutine alloc_cuparm

!===============================================================================

  subroutine dealloc_cuparm()
    implicit none

    if (allocated(thsrc))   deallocate (thsrc)
    if (allocated(rtsrc))   deallocate (rtsrc)
    if (allocated(aconpr))  deallocate (aconpr)
    if (allocated(conprr))  deallocate (conprr)
    if (allocated(cbmf))    deallocate (cbmf)
    if (allocated(kcutop))  deallocate (kcutop)
    if (allocated(kcubot))  deallocate (kcubot)
    if (allocated(vxsrc))   deallocate (vxsrc)
    if (allocated(vysrc))   deallocate (vysrc)
    if (allocated(vzsrc))   deallocate (vzsrc)

  end subroutine dealloc_cuparm

!===============================================================================

  subroutine filltab_cuparm()

    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

     if (allocated(thsrc)) then
        call increment_vtable('THSRC', 'AW')
        vtab_r(num_var)%rvar2_p => thsrc
     endif

     if (allocated(rtsrc)) then
        call increment_vtable('RTSRC', 'AW')
        vtab_r(num_var)%rvar2_p => rtsrc
     endif

     if (allocated(aconpr)) then
        call increment_vtable('ACONPR', 'AW')
        vtab_r(num_var)%rvar1_p => aconpr
     endif

     if (allocated(conprr)) then
        call increment_vtable('CONPRR', 'AW')
        vtab_r(num_var)%rvar1_p => conprr
     endif

     if (allocated(cbmf)) then
        call increment_vtable('CBMF', 'AW')
        vtab_r(num_var)%rvar1_p => cbmf
     endif

     if (allocated(kcutop)) then
        call increment_vtable('KCUTOP', 'AW')
        vtab_r(num_var)%ivar1_p => kcutop
     endif

     if (allocated(kcubot)) then
        call increment_vtable('KCUBOT', 'AW')
        vtab_r(num_var)%ivar1_p => kcubot
     endif

     if (allocated(vxsrc)) then
        call increment_vtable('VXSRC', 'AW')
        vtab_r(num_var)%rvar2_p => vxsrc
     endif

     if (allocated(vysrc)) then
        call increment_vtable('VYSRC', 'AW')
        vtab_r(num_var)%rvar2_p => vysrc
     endif

     if (allocated(vzsrc)) then
        call increment_vtable('VZSRC', 'AW')
        vtab_r(num_var)%rvar2_p => vzsrc
     endif

   end subroutine filltab_cuparm

End Module mem_cuparm
