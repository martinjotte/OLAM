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

Module mem_co2

  integer :: co2flag  ! CO2 flag:  0 = CO2 is not prognosed; 1 = CO2 is prognosed

  real, allocatable, target :: sh_co2(:,:) ! cloud water spec dens [kg_cld/kg_air]

Contains

!===============================================================================

  subroutine alloc_co2(mza,mwa)

    use misc_coms, only: rinit
    use oname_coms, only: nl

    implicit none

    integer, intent(in) :: mza,mwa
    integer :: ic

    ! Allocate arrays based on options
    ! Initialize arrays to zero

    if (co2flag > 0) then
       allocate (sh_co2(mza,mwa)) ; sh_co2 = rinit
    endif

  end subroutine alloc_co2

!===============================================================================

  subroutine dealloc_co2()

    implicit none

    if (allocated(sh_co2)) deallocate (sh_co2)
  end subroutine dealloc_co2

!===============================================================================

  subroutine filltab_co2()

    use var_tables, only: increment_vtable

    implicit none

    if (allocated(sh_co2)) call increment_vtable('SH_CO2', 'AW', mpt1=.true., rvar2=sh_co2)

  end subroutine filltab_co2

!===============================================================================

  subroutine co2init()

!  use misc_coms,   only: io6, runtype, co2flag
  use misc_coms,   only: io6, runtype
  use mem_ijtabs,  only: jtab_w, jtw_init
  use mem_grid,    only: mza
  use consts_coms, only: r8
  use mem_basic,   only: sh_w

  implicit none

  integer :: j, iw

  ! Initialize 3D and 2D microphysics fields

  if (runtype == 'INITIAL' .and. allocated(sh_co2)) then

     do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
        sh_co2(1:mza,iw) = 400.e-6 * (1.0 - sh_w(1:mza,iw))
     enddo

  endif ! runtype == 'INITIAL'

  end subroutine co2init

End Module mem_co2
