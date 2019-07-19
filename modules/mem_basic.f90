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

Module mem_basic

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable :: vmp  (:,:) ! past V horiz momentum [kg/(m^2 s)]
  real, allocatable :: vmc  (:,:) ! current V horiz momentum [kg/(m^2 s)]
  real, allocatable :: vp   (:,:) ! past V horiz velocity [m/s]
  real, allocatable :: vc   (:,:) ! current V horiz velocity [m/s]

  real, allocatable :: wmc  (:,:) ! current vert momentum [kg/(m^2 s)]
  real, allocatable :: wc   (:,:) ! current vert velocity [m/s]
  real, allocatable :: sh_w (:,:) ! tot water spec dens [kg_wat/kg_air]
  real, allocatable :: sh_v (:,:) ! spec hum [kg_vap/kg_air]
  real, allocatable :: thil (:,:) ! ice-liquid pot temp [K]
  real, allocatable :: theta(:,:) ! pot temp [K]
  real, allocatable :: tair (:,:)  ! temperature [K]

  real, allocatable :: vxe  (:,:) ! earth-relative x velocity at T point [m/s]
  real, allocatable :: vye  (:,:) ! earth-relative y velocity at T point [m/s]
  real, allocatable :: vze  (:,:) ! earth-relative z velocity at T point [m/s]

  real, allocatable :: vxe2 (:,:) ! earth-relative x velocity at T point [m/s]
  real, allocatable :: vye2 (:,:) ! earth-relative y velocity at T point [m/s]
  real, allocatable :: vze2 (:,:) ! earth-relative z velocity at T point [m/s]

  real, allocatable :: ue   (:,:) ! easterly wind
  real, allocatable :: ve   (:,:) ! northerly wind

  real(r8), allocatable :: press(:,:) ! air pressure [Pa]
  real(r8), allocatable :: rho  (:,:) ! total air density [kg/m^3]

Contains

!===============================================================================

  subroutine alloc_basic(mza,mva,mwa,nve2_max)

    use misc_coms, only: rinit, rinit8

    implicit none

    integer, intent(in) :: mza,mva,mwa,nve2_max

!   Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
!   and initialize allocated arrays to zero

    allocate (vmp(mza,mva)) ; vmp = rinit
    allocate (vmc(mza,mva)) ; vmc = rinit
    allocate (vc (mza,mva)) ; vc  = rinit
!   allocate (vp (mza,mva)) ; vp  = rinit

    allocate (rho  (mza,mwa)) ; rho   = rinit8
    allocate (press(mza,mwa)) ; press = rinit8
    allocate (wmc  (mza,mwa)) ; wmc   = rinit
    allocate (wc   (mza,mwa)) ; wc    = rinit
    allocate (thil (mza,mwa)) ; thil  = rinit
    allocate (theta(mza,mwa)) ; theta = rinit
    allocate (tair (mza,mwa)) ; tair  = rinit
    allocate (sh_w (mza,mwa)) ; sh_w  = rinit
    allocate (sh_v (mza,mwa)) ; sh_v  = rinit

    allocate (vxe  (mza,mwa)) ; vxe   = rinit
    allocate (vye  (mza,mwa)) ; vye   = rinit
    allocate (vze  (mza,mwa)) ; vze   = rinit

    allocate (vxe2 (nve2_max,mwa)) ; vxe2 = rinit
    allocate (vye2 (nve2_max,mwa)) ; vye2 = rinit
    allocate (vze2 (nve2_max,mwa)) ; vze2 = rinit

!   allocate(ue(mza,mwa)) ; ue = rinit
!   allocate(ve(mza,mwa)) ; ve = rinit

  end subroutine alloc_basic

!===============================================================================

  subroutine dealloc_basic()
    implicit none

    if (allocated(vmp))   deallocate (vmp)
    if (allocated(vmc))   deallocate (vmc)
    if (allocated(vp))    deallocate (vp)
    if (allocated(vc))    deallocate (vc)
    if (allocated(wmc))   deallocate (wmc)
    if (allocated(wc))    deallocate (wc)
    if (allocated(rho))   deallocate (rho)
    if (allocated(sh_w))  deallocate (sh_w)
    if (allocated(sh_v))  deallocate (sh_v)
    if (allocated(press)) deallocate (press)
    if (allocated(thil))  deallocate (thil)
    if (allocated(theta)) deallocate (theta)
    if (allocated(tair))  deallocate (tair)
    if (allocated(vxe))   deallocate (vxe)
    if (allocated(vye))   deallocate (vye)
    if (allocated(vze))   deallocate (vze)
    if (allocated(vxe2))  deallocate (vxe2)
    if (allocated(vye2))  deallocate (vye2)
    if (allocated(vze2))  deallocate (vze2)
    if (allocated(ue))    deallocate (ue)
    if (allocated(ve))    deallocate (ve)

  end subroutine dealloc_basic

!===============================================================================

  subroutine filltab_basic()

    use var_tables, only: increment_vtable
    implicit none

    if (allocated(vmp))   call increment_vtable('VMP',  'AV', rvar2=vmp)

    if (allocated(vmc))   call increment_vtable('VMC',  'AV', rvar2=vmc)

    if (allocated(vp))    call increment_vtable('VP',   'AV', rvar2=vp)

    if (allocated(vc))    call increment_vtable('VC',   'AV', rvar2=vc)

    if (allocated(wmc))   call increment_vtable('WMC',  'AW', rvar2=wmc)

    if (allocated(wc))    call increment_vtable('WC',   'AW', rvar2=wc)

    if (allocated(sh_w))  call increment_vtable('SH_W', 'AW', rvar2=sh_w,  mpt1=.true.)

    if (allocated(sh_v))  call increment_vtable('SH_V', 'AW', rvar2=sh_v,  mpt1=.true.)

    if (allocated(thil))  call increment_vtable('THIL', 'AW', rvar2=thil)

    if (allocated(theta)) call increment_vtable('THETA','AW', rvar2=theta, mpt1=.true.)

    if (allocated(tair))  call increment_vtable('TAIR', 'AW', rvar2=tair,  mpt1=.true.)

    if (allocated(rho))   call increment_vtable('RHO',  'AW', dvar2=rho)

    if (allocated(press)) call increment_vtable('PRESS','AW', dvar2=press)

    if (allocated(vxe2))  call increment_vtable('VXE2', 'AW', rvar2=vxe2)

    if (allocated(vye2))  call increment_vtable('VYE2', 'AW', rvar2=vye2)

    if (allocated(vze2))  call increment_vtable('VZE2', 'AW', rvar2=vze2)

    if (allocated(ue))    call increment_vtable('UE',   'AW', rvar2=ue)

    if (allocated(ve))    call increment_vtable('VE',   'AW', rvar2=ve)

  end subroutine filltab_basic

End Module mem_basic
