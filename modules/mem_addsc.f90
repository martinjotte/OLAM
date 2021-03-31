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
Module mem_addsc

  implicit none

  ! Added scalar variables and tendencies

  Type addsc_vars
     real, allocatable :: sclp(:,:) ! somethings per kg_air
     real, allocatable :: sclt(:,:) ! somethings per kg_air per sec
     real, allocatable :: drydep(:) ! somethings per m^2 per sec (not currently implemented)
  End Type addsc_vars

  type (addsc_vars), allocatable :: addsc(:)

Contains

!===============================================================================

  subroutine alloc_addsc(mza,mwa,naddsc)

    use misc_coms, only: io6, rinit
    implicit none

    integer, intent(in) :: mza,mwa,naddsc
    integer             :: iaddsc

    ! Allocate arrays based on options (if necessary).

    allocate (addsc(naddsc))

    do iaddsc = 1,naddsc

       write(io6,*) 'alloc_addsc ', iaddsc, mza, mwa, naddsc

       allocate (addsc(iaddsc)%sclp(mza,mwa)) ; call random_number(addsc(iaddsc)%sclp)
       allocate (addsc(iaddsc)%drydep  (mwa)) ; addsc(iaddsc)%drydep = rinit

       addsc(iaddsc)%sclp = addsc(iaddsc)%sclp + 2.

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

    use var_tables, only: increment_vtable
    implicit none

    integer, intent(in) :: naddsc

    integer      :: iaddsc
    character(7) :: sname

    do iaddsc = 1,naddsc

       if (allocated (addsc(iaddsc)%sclp)) then

          write(sname,'(a4,i3.3)') 'SCLP', iaddsc
          call increment_vtable(sname, 'AW', mpt1=.true., rvar2=addsc(iaddsc)%sclp)

       endif

       if (allocated (addsc(iaddsc)%drydep)) then

          write(sname,'(a4,i3.3)') 'SCDD', iaddsc
          call increment_vtable(sname, 'AW', rvar1=addsc(iaddsc)%drydep)

       endif

    enddo

  end subroutine filltab_addsc

End Module mem_addsc
