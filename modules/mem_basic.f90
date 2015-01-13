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

  real, allocatable, target :: vmp  (:,:) ! past V horiz momentum [kg/(m^2 s)]
  real, allocatable, target :: vmc  (:,:) ! current V horiz momentum [kg/(m^2 s)]
  real, allocatable, target :: vp   (:,:) ! past V horiz velocity [m/s]
  real, allocatable, target :: vc   (:,:) ! current V horiz velocity [m/s]

  real, allocatable, target :: wmc  (:,:) ! current vert momentum [kg/(m^2 s)]
  real, allocatable, target :: wc   (:,:) ! current vert velocity [m/s]
  real, allocatable, target :: sh_w (:,:) ! tot water spec dens [kg_wat/kg_air]
  real, allocatable, target :: sh_v (:,:) ! spec hum [kg_vap/kg_air]
  real, allocatable, target :: thil (:,:) ! ice-liquid pot temp [K]
  real, allocatable, target :: theta(:,:) ! pot temp [K]
  real, allocatable, target :: tair (:,:)  ! temperature [K]

  real, allocatable, target :: vxe  (:,:) ! earth-relative x velocity at T point [m/s]
  real, allocatable, target :: vye  (:,:) ! earth-relative y velocity at T point [m/s]
  real, allocatable, target :: vze  (:,:) ! earth-relative z velocity at T point [m/s]

  real(r8), allocatable, target :: press(:,:) ! air pressure [Pa]
  real(r8), allocatable, target :: rho  (:,:) ! total air density [kg/m^3]

  ! If false, use current earth-cartesian velocities to compute the w, v, and t
  ! donor point locations, which saves some computation, memory, and communication.
  ! If true, compute half-forward earth-cartesian velocities, which is more
  ! exact for the time differencing scheme:

  logical, parameter :: strict_wvt_donorpoint = .false.

Contains

!===============================================================================

  subroutine alloc_basic(mza,mva,mwa)
    use misc_coms, only: rinit, rinit8
    implicit none

    integer, intent(in) :: mza,mva,mwa

!   Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
!   and initialize allocated arrays to zero

    allocate (vmp(mza,mva)) ; vmp = rinit
    allocate (vmc(mza,mva)) ; vmc = rinit
    allocate (vc (mza,mva)) ; vc  = rinit

    if (strict_wvt_donorpoint) then
       ! needed for half-forward earth-cartesian velocities:
       allocate (vp (mza,mva)) ; vp  = rinit
    endif

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

  end subroutine filltab_basic

End Module mem_basic
