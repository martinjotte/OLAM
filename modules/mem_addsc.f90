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
Module mem_addsc
   
   ! Added scalar variables and tendencies

   Type addsc_vars   
      real, allocatable :: sclp(:,:) ! somethings per kg_air
      real, allocatable :: sclt(:,:) ! somethings per kg_air per sec
      real, allocatable :: drydep(:) ! somethings per m^2 per sec (not currently implemented)
   End Type addsc_vars
   
   type (addsc_vars), allocatable, target :: addsc(:)
  
Contains

!===============================================================================

   subroutine alloc_addsc(mza,mwa,naddsc)

     use misc_coms, only: io6, rinit
     implicit none

     integer, intent(in) :: mza,mwa,naddsc
     integer             :: iaddsc

!    Allocate arrays based on options (if necessary).

     allocate (addsc(naddsc))

     do iaddsc = 1,naddsc
   
        write(io6,*) 'alloc_addsc ', iaddsc, mza, mwa, naddsc
   
        allocate (addsc(iaddsc)%sclp(mza,mwa)) ; addsc(iaddsc)%sclp   = rinit
        allocate (addsc(iaddsc)%drydep  (mwa)) ; addsc(iaddsc)%drydep = rinit

     enddo
   end subroutine alloc_addsc

!===============================================================================

   subroutine dealloc_addsc(naddsc)

   implicit none

   integer, intent(in) :: naddsc
   integer             :: iaddsc

!  Deallocate arrays

   if (allocated(addsc)) then
      
      do iaddsc = 1,naddsc
         if (allocated(addsc(iaddsc)%sclp  )) deallocate (addsc(iaddsc)%sclp  )
         if (allocated(addsc(iaddsc)%drydep)) deallocate (addsc(iaddsc)%drydep)
      enddo

      deallocate(addsc)

   endif
   end subroutine dealloc_addsc

!===============================================================================
               
   subroutine filltab_addsc(naddsc)
     use var_tables, only: vtab_r, num_var, increment_vtable
     implicit none

     integer, intent(in) :: naddsc

     integer      :: iaddsc
     character(7) :: sname

     do iaddsc = 1,naddsc

        if (allocated (addsc(iaddsc)%sclp)) then

           write(sname,'(a4,i3.3)') 'SCLP', iaddsc
           call increment_vtable(sname, 'AW', mpt1=.true.)
           vtab_r(num_var)%rvar2_p => addsc(iaddsc)%sclp

        endif
        
        if (allocated (addsc(iaddsc)%drydep)) then

           write(sname,'(a4,i3.3)') 'SCDD', iaddsc
           call increment_vtable(sname, 'AW')
           vtab_r(num_var)%rvar1_p => addsc(iaddsc)%drydep
           
        endif

     enddo
   end subroutine filltab_addsc

End Module mem_addsc
