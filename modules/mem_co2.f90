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

  use const_data, only: mwair

  ! CO2 flag:  0 = CO2 is not prognosed; 1 = CO2 is prognosed
  integer           :: co2flag     = 0
  real              :: co2_initppm = 400.0  ! ppmv = mol / 10^6 mol

  real, parameter   :: co2_mw      = 44.0  ! g / mol
  real, parameter   :: co2_ppm2sh  = co2_mw / mwair  * 1.e-6
  real, parameter   :: co2_sh2ppm  = mwair  / co2_mw * 1.e+6

  real, allocatable :: rr_co2(:,:) ! CO2 specific density [kg_CO2/kg_air]

  integer           :: i_co2 = 0  ! index of co2 in scalar table

  private           :: mwair

Contains

!===============================================================================

  subroutine alloc_co2(mza,mwa)

    use misc_coms, only: rinit

    implicit none

    integer, intent(in) :: mza,mwa

    ! Allocate arrays based on options
    ! Initialize arrays to zero

    if (co2flag > 0) then
       allocate (rr_co2(mza,mwa)) ; rr_co2 = rinit
    endif

  end subroutine alloc_co2

!===============================================================================

  subroutine dealloc_co2()

    implicit none

    if (allocated(rr_co2)) deallocate (rr_co2)
  end subroutine dealloc_co2

!===============================================================================

  subroutine filltab_co2()

    use var_tables, only: increment_vtable

    implicit none

    if (allocated(rr_co2)) call increment_vtable('RR_CO2', 'AW', mpt1=.true., rvar2=rr_co2)

  end subroutine filltab_co2

!===============================================================================

  subroutine co2init()

  use misc_coms,   only: runtype
  use mem_ijtabs,  only: jtab_w, jtw_init
  use mem_grid,    only: mza
  use consts_coms, only: r8

  implicit none

  integer :: j, iw
  real    :: co2_initsh

  ! Initialize 3D and 2D microphysics fields

  if (runtype == 'INITIAL' .and. allocated(rr_co2) .and. co2flag /= 0) then

     co2_initsh = co2_initppm * co2_ppm2sh

     !$omp parallel do private(iw)
     do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
        rr_co2(2:mza,iw) = co2_initsh
     enddo
     !$omp end parallel do

  endif ! runtype == 'INITIAL'

  end subroutine co2init

End Module mem_co2
