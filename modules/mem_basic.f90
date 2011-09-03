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

Module mem_basic

  use consts_coms, only: r8
  implicit none

  private :: r8

  real, allocatable, target :: ump  (:,:) ! past U horiz momentum [kg/(m^2 s)]
  real, allocatable, target :: umc  (:,:) ! current U horiz momentum [kg/(m^2 s)]
  real, allocatable, target :: uc   (:,:) ! current U horiz velocity [m/s]
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

  real(r8), allocatable, target :: press(:,:) ! air pressure [Pa]
  real(r8), allocatable, target :: rho  (:,:) ! total air density [kg/m^3]

Contains

!===============================================================================

  subroutine alloc_basic(meshtype,mza,mua,mva,mwa)
    use misc_coms, only: rinit, rinit8
    implicit none

    integer, intent(in) :: meshtype,mza,mua,mva,mwa

!   Allocate basic memory needed for 'INITIAL' or 'HISTORY' runs
!   and initialize allocated arrays to zero

    if (meshtype == 1) then

       allocate (ump(mza,mua)) ; ump = rinit
       allocate (umc(mza,mua)) ; umc = rinit
      
    elseif (meshtype == 2) then
   
       allocate (vmp(mza,mva)) ; vmp = rinit
       allocate (vmc(mza,mva)) ; vmc = rinit
       allocate (vp (mza,mva)) ; vp  = rinit

    endif
   
    allocate (rho  (mza,mwa)) ; rho   = rinit8
    allocate (press(mza,mwa)) ; press = rinit8

    allocate (uc   (mza,mua)) ; uc    = rinit
    allocate (vc   (mza,mva)) ; vc    = rinit
    allocate (wmc  (mza,mwa)) ; wmc   = rinit
    allocate (wc   (mza,mwa)) ; wc    = rinit
    allocate (thil (mza,mwa)) ; thil  = rinit
    allocate (theta(mza,mwa)) ; theta = rinit
    allocate (sh_w (mza,mwa)) ; sh_w  = rinit
    allocate (sh_v (mza,mwa)) ; sh_v  = rinit

  end subroutine alloc_basic

!===============================================================================

  subroutine dealloc_basic()
    implicit none

    if (allocated(ump))   deallocate (ump)
    if (allocated(umc))   deallocate (umc)
    if (allocated(uc))    deallocate (uc)
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

  end subroutine dealloc_basic

!===============================================================================

  subroutine filltab_basic()
    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(ump)) then
       call increment_vtable('UMP', 'AU')
       vtab_r(num_var)%rvar2_p => ump
    endif

    if (allocated(umc)) then
       call increment_vtable('UMC', 'AU')
       vtab_r(num_var)%rvar2_p => umc
    endif
    
    if (allocated(uc)) then
       call increment_vtable('UC',  'AU')
       vtab_r(num_var)%rvar2_p => uc
    endif

    if (allocated(vmp)) then
       call increment_vtable('VMP', 'AU')
       vtab_r(num_var)%rvar2_p => vmp
    endif

    if (allocated(vmc)) then
       call increment_vtable('VMC', 'AU')
       vtab_r(num_var)%rvar2_p => vmc
    endif

    if (allocated(vp)) then
       call increment_vtable('VP',  'AU')
       vtab_r(num_var)%rvar2_p => vp
    endif

    if (allocated(vc)) then
       call increment_vtable('VC',  'AU')
       vtab_r(num_var)%rvar2_p => vc
    endif

    if (allocated(wmc)) then
       call increment_vtable('WMC', 'AW')
       vtab_r(num_var)%rvar2_p => wmc
    endif

    if (allocated(wc)) then
       call increment_vtable('WC',  'AW')
       vtab_r(num_var)%rvar2_p => wc
    endif

    if (allocated(sh_w)) then
       call increment_vtable('SH_W', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_w
    endif

    if (allocated(sh_v)) then
       call increment_vtable('SH_V', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => sh_v
    endif

    if (allocated(thil)) then
       call increment_vtable('THIL',  'AW')
       vtab_r(num_var)%rvar2_p => thil
    endif

    if (allocated(theta)) then
       call increment_vtable('THETA', 'AW', mpt1=.true.)
       vtab_r(num_var)%rvar2_p => theta
    endif

    if (allocated(rho)) then
       call increment_vtable('RHO',   'AW')
       vtab_r(num_var)%dvar2_p => rho
    endif

    if (allocated(press)) then
       call increment_vtable('PRESS', 'AW')
       vtab_r(num_var)%dvar2_p => press
    endif

  end subroutine filltab_basic

End Module mem_basic
