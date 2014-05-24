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

real(r8) :: time8_prev0
real(r8) :: time8_prev1
real(r8) :: time8_prev2
real(r8) :: time8_prev3

real, allocatable :: accpmic_prev0(:)
real, allocatable :: accpmic_prev1(:)
real, allocatable :: accpmic_prev2(:)
real, allocatable :: accpmic_prev3(:)

real, allocatable :: accpcon_prev0(:)
real, allocatable :: accpcon_prev1(:)
real, allocatable :: accpcon_prev2(:)
real, allocatable :: accpcon_prev3(:)

real, allocatable :: rshort_accum_prev0(:)
real, allocatable :: rshort_accum_prev1(:)

real, allocatable :: rshortup_accum_prev0(:)
real, allocatable :: rshortup_accum_prev1(:)

real, allocatable :: rlong_accum_prev0(:)
real, allocatable :: rlong_accum_prev1(:)

real, allocatable :: rlongup_accum_prev0(:)
real, allocatable :: rlongup_accum_prev1(:)

real, allocatable :: rshort_top_accum_prev0(:)
real, allocatable :: rshort_top_accum_prev1(:)

real, allocatable :: rshortup_top_accum_prev0(:)
real, allocatable :: rshortup_top_accum_prev1(:)

real, allocatable :: rlongup_top_accum_prev0(:)
real, allocatable :: rlongup_top_accum_prev1(:)

real, allocatable :: sflux_t_accum_prev0(:)
real, allocatable :: sflux_t_accum_prev1(:)

real, allocatable :: sflux_r_accum_prev0(:)
real, allocatable :: sflux_r_accum_prev1(:)

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

  allocate ( accpmic_prev0(mwa))
  allocate ( accpmic_prev1(mwa))
  allocate ( accpmic_prev2(mwa))
  allocate ( accpmic_prev3(mwa))

  allocate ( accpcon_prev0(mwa))
  allocate ( accpcon_prev1(mwa))
  allocate ( accpcon_prev2(mwa))
  allocate ( accpcon_prev3(mwa))

  allocate ( rshort_accum_prev0(mwa))
  allocate ( rshort_accum_prev1(mwa))

  allocate ( rshortup_accum_prev0(mwa))
  allocate ( rshortup_accum_prev1(mwa))

  allocate ( rlong_accum_prev0(mwa))
  allocate ( rlong_accum_prev1(mwa))

  allocate ( rlongup_accum_prev0(mwa))
  allocate ( rlongup_accum_prev1(mwa))

  allocate ( rshort_top_accum_prev0(mwa))
  allocate ( rshort_top_accum_prev1(mwa))

  allocate ( rshortup_top_accum_prev0(mwa))
  allocate ( rshortup_top_accum_prev1(mwa))

  allocate ( rlongup_top_accum_prev0(mwa))
  allocate ( rlongup_top_accum_prev1(mwa))

  allocate ( sflux_t_accum_prev0(mwa))
  allocate ( sflux_t_accum_prev1(mwa))

  allocate ( sflux_r_accum_prev0(mwa))
  allocate ( sflux_r_accum_prev1(mwa))

  allocate ( press_init(mza,mwa))
  allocate (   rho_init(mza,mwa))
  allocate ( theta_init(mza,mwa))
  allocate (    uc_init(mza,mua))
  allocate (    vc_init(mza,mva))
! allocate (addsc1_init(mza,mwa))

  accpmic_prev0(:) = 0.
  accpmic_prev1(:) = 0.
  accpmic_prev2(:) = 0.
  accpmic_prev3(:) = 0.

  accpcon_prev0(:) = 0.
  accpcon_prev1(:) = 0.
  accpcon_prev2(:) = 0.
  accpcon_prev3(:) = 0.

  rshort_accum_prev0(:) = 0.
  rshort_accum_prev1(:) = 0.

  rshortup_accum_prev0(:) = 0.
  rshortup_accum_prev1(:) = 0.

  rlong_accum_prev0(:) = 0.
  rlong_accum_prev1(:) = 0.

  rlongup_accum_prev0(:) = 0.
  rlongup_accum_prev1(:) = 0.

  rshort_top_accum_prev0(:) = 0.
  rshort_top_accum_prev1(:) = 0.

  rshortup_top_accum_prev0(:) = 0.
  rshortup_top_accum_prev1(:) = 0.

  rlongup_top_accum_prev0(:) = 0.
  rlongup_top_accum_prev1(:) = 0.

  sflux_t_accum_prev0(:) = 0.
  sflux_t_accum_prev1(:) = 0.

  sflux_r_accum_prev0(:) = 0.
  sflux_r_accum_prev1(:) = 0.

  press_init(:,:) = press(:,:)
  rho_init  (:,:) = rho  (:,:)
  theta_init(:,:) = theta(:,:)

  if (allocated(uc)) uc_init(:,:) = uc(:,:)
  if (allocated(vc)) vc_init(:,:) = vc(:,:)
!   addsc1_init(:,:) = addsc(1)%sclp(:,:)

  time8_prev0 = 0.
  time8_prev1 = 0.
  time8_prev2 = 0.
  time8_prev3 = 0.

  end subroutine alloc_plot

!===============================================================================

  subroutine copy_plot(iplt_file)

  use mem_cuparm, only: aconpr

  use mem_micro,      only: accpd, accpr, accpp, accps, accpa, accpg, accph
  use misc_coms,      only: time8
  use mem_turb,       only: sflux_t, sflux_r
  use mem_flux_accum, only: rshort_accum, rshortup_accum, &
                            rlong_accum, rlongup_accum, &
                            rshort_top_accum, rshortup_top_accum, &
                            rlongup_top_accum, &
                            sflux_t_accum, sflux_r_accum

  implicit none

  integer, intent(in) :: iplt_file ! File index for 'PLOTONLY' runs; 0 otherwise

! The following IF/ENDIF statements (but not the lines between them) should be
! commented out for most applications.  Uncommenting them provides the ability
! to sum precipitation values from multiple history files into a single
! accpmic_prev0 or accpcon_prev0 value, before later differencing between 
! prev0, prev1, prev2, and prev3 groupings.  For example, the first application
! of this capability was to sum the accumulated precipitation of 'Simulation A'
! each March 1 over 5 consecutive years and store the sum in the prev0 arrays,
! copy the prev0 sum into prev1 and then sum the March 1 precipitation totals
! for 'Simulation B' into prev0, copy the forgoing into prev2 and prev1,
! respectively, and sum June 1 totals for 'Simulation A' into prev0, and finaly
! copy the foregoing into prev3, prev2, and prev1, respectively, and sum June 1
! totals from 'Simulation B' into prev0.  Oplot_lib then plotted differences of
! the above 4 totals to get the difference in springtime precipitation
! from A to B.

!  if (iplt_file ==  6 .or. &
!      iplt_file == 11 .or. &
!      iplt_file == 16 .or. &
!      iplt_file == 21 .or. &
!      iplt_file == 26 .or. &
!      iplt_file == 31 .or. &
!      iplt_file == 36 .or. &
!      iplt_file == 41 .or. &
!      iplt_file == 46) then

! Shift previous values (#3's are oldest)

     time8_prev3 = time8_prev2
     time8_prev2 = time8_prev1
     time8_prev1 = time8_prev0

     accpmic_prev3(:) = accpmic_prev2(:)
     accpmic_prev2(:) = accpmic_prev1(:)
     accpmic_prev1(:) = accpmic_prev0(:)

     accpcon_prev3(:) = accpcon_prev2(:)
     accpcon_prev2(:) = accpcon_prev1(:)
     accpcon_prev1(:) = accpcon_prev0(:)

           rshort_accum_prev1(:) =       rshort_accum_prev0(:)
         rshortup_accum_prev1(:) =     rshortup_accum_prev0(:)
            rlong_accum_prev1(:) =        rlong_accum_prev0(:)
          rlongup_accum_prev1(:) =      rlongup_accum_prev0(:)
       rshort_top_accum_prev1(:) =   rshort_top_accum_prev0(:)
     rshortup_top_accum_prev1(:) = rshortup_top_accum_prev0(:)
      rlongup_top_accum_prev1(:) =  rlongup_top_accum_prev0(:)
          sflux_t_accum_prev1(:) =      sflux_t_accum_prev0(:)
          sflux_r_accum_prev1(:) =      sflux_r_accum_prev0(:)

! Fill most recent previous values

     time8_prev0 = 0.

     accpmic_prev0(:) = 0.
     accpcon_prev0(:) = 0.

           rshort_accum_prev0(:) = 0.
         rshortup_accum_prev0(:) = 0.
            rlong_accum_prev0(:) = 0.
          rlongup_accum_prev0(:) = 0.
       rshort_top_accum_prev0(:) = 0.
     rshortup_top_accum_prev0(:) = 0.
      rlongup_top_accum_prev0(:) = 0.
          sflux_t_accum_prev0(:) = 0.
          sflux_r_accum_prev0(:) = 0.

!  endif

  time8_prev0 = time8_prev0 + time8

  if (allocated(accpd)) accpmic_prev0(:) = accpmic_prev0(:) + real(accpd(:))
  if (allocated(accpr)) accpmic_prev0(:) = accpmic_prev0(:) + real(accpr(:))
  if (allocated(accpp)) accpmic_prev0(:) = accpmic_prev0(:) + real(accpp(:))
  if (allocated(accps)) accpmic_prev0(:) = accpmic_prev0(:) + real(accps(:))
  if (allocated(accpa)) accpmic_prev0(:) = accpmic_prev0(:) + real(accpa(:))
  if (allocated(accpg)) accpmic_prev0(:) = accpmic_prev0(:) + real(accpg(:))
  if (allocated(accph)) accpmic_prev0(:) = accpmic_prev0(:) + real(accph(:))

  if (allocated(aconpr)) accpcon_prev0(:) = accpcon_prev0(:) + real(aconpr(:))

  if (allocated(rshort_accum)) rshort_accum_prev0(:) = &
                               rshort_accum_prev0(:) + &
                          real(rshort_accum(:))

  if (allocated(rshortup_accum)) rshortup_accum_prev0(:) = &
                                 rshortup_accum_prev0(:) + &
                            real(rshortup_accum(:))

  if (allocated(rlong_accum)) rlong_accum_prev0(:) = &
                              rlong_accum_prev0(:) + &
                         real(rlong_accum(:))

  if (allocated(rlongup_accum)) rlongup_accum_prev0(:) = &
                                rlongup_accum_prev0(:) + &
                           real(rlongup_accum(:))

  if (allocated(rshort_top_accum)) rshort_top_accum_prev0(:) = &
                                   rshort_top_accum_prev0(:) + &
                              real(rshort_top_accum(:))

  if (allocated(rshortup_top_accum)) rshortup_top_accum_prev0(:) = &
                                     rshortup_top_accum_prev0(:) + &
                                real(rshortup_top_accum(:))

  if (allocated(rlongup_top_accum)) rlongup_top_accum_prev0(:) = &
                                    rlongup_top_accum_prev0(:) + &
                               real(rlongup_top_accum(:))

  if (allocated(sflux_t_accum)) sflux_t_accum_prev0(:) = &
                                sflux_t_accum_prev0(:) + &
                           real(sflux_t_accum(:))

  if (allocated(sflux_r_accum)) sflux_r_accum_prev0(:) = &
                                sflux_r_accum_prev0(:) + &
                           real(sflux_r_accum(:))

  end subroutine copy_plot

End Module mem_plot
