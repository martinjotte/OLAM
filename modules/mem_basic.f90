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

  real, allocatable :: vmc  (:,:) ! current V horiz momentum [kg/(m^2 s)]
  real, allocatable :: vc   (:,:) ! current V horiz velocity [m/s]

  real, allocatable :: wmc  (:,:) ! current vert momentum [kg/(m^2 s)]
  real, allocatable :: wc   (:,:) ! current vert velocity [m/s]
  real, allocatable :: rr_w (:,:) ! tot water mixing ratio [kg_wat/kg_dryair]
  real, allocatable :: rr_v (:,:) ! water vapor mixing ratio [kg_vap/kg_dryair]
  real, allocatable :: thil (:,:) ! ice-liquid pot temp [K]
  real, allocatable :: theta(:,:) ! pot temp [K]
  real, allocatable :: tair (:,:) ! temperature [K]

  real, allocatable :: vxe  (:,:) ! earth-relative x velocity at T point [m/s]
  real, allocatable :: vye  (:,:) ! earth-relative y velocity at T point [m/s]
  real, allocatable :: vze  (:,:) ! earth-relative z velocity at T point [m/s]

  real, allocatable :: ue   (:,:) ! easterly wind
  real, allocatable :: ve   (:,:) ! northerly wind

  real(r8), allocatable :: press(:,:) ! air pressure [Pa]
  real(r8), allocatable :: rho  (:,:) ! dry air density [kg/m^3]

Contains

!===============================================================================

  subroutine alloc_basic(mza,mva,mwa)

    use misc_coms, only: rinit, rinit8

    implicit none

    integer, intent(in) :: mza,mva,mwa

!   Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
!   and initialize allocated arrays to zero

    allocate (vmc  (mza,mva)) ; vmc = rinit
    allocate (vc   (mza,mva)) ; vc  = rinit

    allocate (rho  (mza,mwa)) ; rho   = 0.0_r8 !rinit
    allocate (press(mza,mwa)) ; press = rinit8
    allocate (wmc  (mza,mwa)) ; wmc   = rinit
    allocate (wc   (mza,mwa)) ; wc    = 0.0 !rinit
    allocate (thil (mza,mwa)) ; thil  = rinit
    allocate (theta(mza,mwa)) ; theta = rinit
    allocate (tair (mza,mwa)) ; tair  = rinit
    allocate (rr_w (mza,mwa)) ; rr_w  = rinit
    allocate (rr_v (mza,mwa)) ; rr_v  = rinit

    allocate (vxe  (mza,mwa)) ; vxe   = rinit
    allocate (vye  (mza,mwa)) ; vye   = rinit
    allocate (vze  (mza,mwa)) ; vze   = rinit

    allocate(ue    (mza,mwa)) ; ue    = rinit
    allocate(ve    (mza,mwa)) ; ve    = rinit

  end subroutine alloc_basic

!===============================================================================

  subroutine dealloc_basic()
    implicit none

    if (allocated(vmc))   deallocate (vmc)
    if (allocated(vc))    deallocate (vc)

    if (allocated(wmc))   deallocate (wmc)
    if (allocated(wc))    deallocate (wc)
    if (allocated(rho))   deallocate (rho)
    if (allocated(rr_w))  deallocate (rr_w)
    if (allocated(rr_v))  deallocate (rr_v)
    if (allocated(press)) deallocate (press)
    if (allocated(thil))  deallocate (thil)
    if (allocated(theta)) deallocate (theta)
    if (allocated(tair))  deallocate (tair)

    if (allocated(vxe))   deallocate (vxe)
    if (allocated(vye))   deallocate (vye)
    if (allocated(vze))   deallocate (vze)

    if (allocated(ue))    deallocate (ue)
    if (allocated(ve))    deallocate (ve)

  end subroutine dealloc_basic

!===============================================================================

  subroutine filltab_basic()

    use var_tables, only: increment_vtable
    implicit none

    if (allocated(vmc))   call increment_vtable('VMC',  'AV', rvar2=vmc)

    if (allocated(vc))    call increment_vtable('VC',   'AV', rvar2=vc)

    if (allocated(wmc))   call increment_vtable('WMC',  'AW', rvar2=wmc)

    if (allocated(wc))    call increment_vtable('WC',   'AW', rvar2=wc)

    if (allocated(rr_w))  call increment_vtable('RR_W', 'AW', rvar2=rr_w,  mpt1=.true.)

    if (allocated(rr_v))  call increment_vtable('RR_V', 'AW', rvar2=rr_v,  mpt1=.true.)

    if (allocated(thil))  call increment_vtable('THIL', 'AW', rvar2=thil)

    if (allocated(theta)) call increment_vtable('THETA','AW', rvar2=theta, mpt1=.true.)

    if (allocated(tair))  call increment_vtable('TAIR', 'AW', rvar2=tair,  mpt1=.true.)

    if (allocated(rho))   call increment_vtable('RHO',  'AW', dvar2=rho)

    if (allocated(press)) call increment_vtable('PRESS','AW', dvar2=press)

    if (allocated(vxe))   call increment_vtable('VXE',  'AW', rvar2=vxe)

    if (allocated(vye))   call increment_vtable('VYE',  'AW', rvar2=vye)

    if (allocated(vze))   call increment_vtable('VZE',  'AW', rvar2=vze)

    if (allocated(ue))    call increment_vtable('UE',   'AW', rvar2=ue)

    if (allocated(ve))    call increment_vtable('VE',   'AW', rvar2=ve)

  end subroutine filltab_basic

End Module mem_basic
