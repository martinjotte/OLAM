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

   real, allocatable, target :: thsrc  (:,:)
   real, allocatable, target :: rtsrc  (:,:)
   real, allocatable, target :: thsrcsh(:,:)
   real, allocatable, target :: rtsrcsh(:,:)
   real, allocatable, target :: aconpr   (:)
   real, allocatable, target :: conprr   (:)

Contains

!===============================================================================

  subroutine alloc_cuparm(mza,mwa,nqparm,nqparm_sh)
    use misc_coms, only: rinit
    implicit none

    integer, intent(in) :: mza,mwa
    integer, intent(in) :: nqparm,nqparm_sh ! these are max values over all MRLs
   
!   Allocate all cuparm arrays if either deep or shallow cuparm is activated
!   on any MRL
      
    if (nqparm > 0 .or. nqparm_sh > 0 ) then
       allocate (thsrc(mza,mwa)) ; thsrc = rinit
       allocate (rtsrc(mza,mwa)) ; rtsrc = rinit

       allocate (aconpr(mwa)) ; aconpr = rinit
       allocate (conprr(mwa)) ; conprr = rinit

       allocate (thsrcsh(mza,mwa)) ; thsrcsh = rinit
       allocate (rtsrcsh(mza,mwa)) ; rtsrcsh = rinit
    endif

  end subroutine alloc_cuparm

!===============================================================================

  subroutine dealloc_cuparm()
    implicit none

    if (allocated(thsrc))   deallocate (thsrc)
    if (allocated(rtsrc))   deallocate (rtsrc)
    if (allocated(thsrcsh)) deallocate (thsrcsh)
    if (allocated(rtsrcsh)) deallocate (rtsrcsh)
    if (allocated(aconpr))  deallocate (aconpr)
    if (allocated(conprr))  deallocate (conprr)

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

     if (allocated(thsrcsh)) then
        call increment_vtable('THSRCSH', 'AW')
        vtab_r(num_var)%rvar2_p => thsrcsh
     endif

     if (allocated(rtsrcsh)) then
        call increment_vtable('RTSRCSH', 'AW')
        vtab_r(num_var)%rvar2_p => rtsrcsh
     endif

     if (allocated(aconpr)) then
        call increment_vtable('ACONPR', 'AW')
        vtab_r(num_var)%rvar1_p => aconpr
     endif

     if (allocated(conprr)) then
        call increment_vtable('CONPRR', 'AW')
        vtab_r(num_var)%rvar1_p => conprr
     endif

   end subroutine filltab_cuparm

End Module mem_cuparm
