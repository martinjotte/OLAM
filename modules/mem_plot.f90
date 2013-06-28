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

Module mem_plot

use consts_coms, only: r8

implicit none

real(r8) :: time8_prev1
real(r8) :: time8_prev2
real(r8) :: time8_prev3

real, allocatable :: accpmic_prev1(:)
real, allocatable :: accpmic_prev2(:)
real, allocatable :: accpmic_prev3(:)

real, allocatable :: accpcon_prev1(:)
real, allocatable :: accpcon_prev2(:)
real, allocatable :: accpcon_prev3(:)

real(r8), allocatable ::  press_init(:,:)
real(r8), allocatable ::    rho_init(:,:)
real,     allocatable ::  theta_init(:,:)
real,     allocatable ::     uc_init(:,:)
real,     allocatable ::     vc_init(:,:)
real,     allocatable :: addsc1_init(:,:)

Contains

!===============================================================================

  subroutine alloc_plot()

  use mem_grid,  only: mza, mua, mva, mwa
  use mem_basic, only: press, rho, theta, uc, vc
  use mem_addsc, only: addsc

  implicit none

  allocate ( accpmic_prev1(mwa))
  allocate ( accpmic_prev2(mwa))
  allocate ( accpmic_prev3(mwa))

  allocate ( accpcon_prev1(mwa))
  allocate ( accpcon_prev2(mwa))
  allocate ( accpcon_prev3(mwa))

  allocate ( press_init(mza,mwa))
  allocate (   rho_init(mza,mwa))
  allocate ( theta_init(mza,mwa))
  allocate (    uc_init(mza,mua))
  allocate (    vc_init(mza,mva))
! allocate (addsc1_init(mza,mwa))

   accpmic_prev1(:) = 0.
   accpmic_prev2(:) = 0.
   accpmic_prev3(:) = 0.

   accpcon_prev1(:) = 0.
   accpcon_prev2(:) = 0.
   accpcon_prev3(:) = 0.

  press_init(:,:) = press(:,:)
  rho_init  (:,:) = rho  (:,:)
  theta_init(:,:) = theta(:,:)

  if (allocated(uc)) uc_init(:,:) = uc(:,:)
  if (allocated(vc)) vc_init(:,:) = vc(:,:)
!   addsc1_init(:,:) = addsc(1)%sclp(:,:)

  time8_prev1 = 0.
  time8_prev2 = 0.
  time8_prev3 = 0.

  end subroutine alloc_plot

!===============================================================================

  subroutine copy_plot()

  use mem_cuparm, only: aconpr

  use mem_micro,  only: accpd, accpr, accpp, accps, accpa, accpg, accph
  use misc_coms, only: time8

  implicit none
 
! Shift previous values (#3's are oldest)

  time8_prev3 = time8_prev2
  time8_prev2 = time8_prev1

  accpmic_prev3(:) = accpmic_prev2(:)
  accpmic_prev2(:) = accpmic_prev1(:)

  accpcon_prev3(:) = accpcon_prev2(:)
  accpcon_prev2(:) = accpcon_prev1(:)

! Fill most recent previous values

  time8_prev1 = time8

  accpmic_prev1(:) = 0.
                     
  if (allocated(accpd)) accpmic_prev1(:) = accpmic_prev1(:) + accpd(:)
  if (allocated(accpr)) accpmic_prev1(:) = accpmic_prev1(:) + accpr(:)
  if (allocated(accpp)) accpmic_prev1(:) = accpmic_prev1(:) + accpp(:)
  if (allocated(accps)) accpmic_prev1(:) = accpmic_prev1(:) + accps(:)
  if (allocated(accpa)) accpmic_prev1(:) = accpmic_prev1(:) + accpa(:)
  if (allocated(accpg)) accpmic_prev1(:) = accpmic_prev1(:) + accpg(:)
  if (allocated(accph)) accpmic_prev1(:) = accpmic_prev1(:) + accph(:)

  if (allocated(aconpr)) accpcon_prev1(:) = aconpr(:)

  end subroutine copy_plot

End Module mem_plot
