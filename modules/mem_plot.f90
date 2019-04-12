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

real, allocatable :: rshort_clr_accum_prev0(:)
real, allocatable :: rshort_clr_accum_prev1(:)

real, allocatable :: rshortup_clr_accum_prev0(:)
real, allocatable :: rshortup_clr_accum_prev1(:)

real, allocatable :: rlong_clr_accum_prev0(:)
real, allocatable :: rlong_clr_accum_prev1(:)

real, allocatable :: rlongup_clr_accum_prev0(:)
real, allocatable :: rlongup_clr_accum_prev1(:)

real, allocatable :: rshort_top_clr_accum_prev0(:)
real, allocatable :: rshort_top_clr_accum_prev1(:)

real, allocatable :: rshortup_top_clr_accum_prev0(:)
real, allocatable :: rshortup_top_clr_accum_prev1(:)

real, allocatable :: rlongup_top_clr_accum_prev0(:)
real, allocatable :: rlongup_top_clr_accum_prev1(:)

real, allocatable :: sfluxt_accum_prev0(:)
real, allocatable :: sfluxt_accum_prev1(:)

real, allocatable :: sfluxr_accum_prev0(:)
real, allocatable :: sfluxr_accum_prev1(:)

real, allocatable :: vc_accum_prev0(:,:)
real, allocatable :: vc_accum_prev1(:,:)

real, allocatable :: wc_accum_prev0(:,:)
real, allocatable :: wc_accum_prev1(:,:)

real, allocatable :: press_accum_prev0(:,:)
real, allocatable :: press_accum_prev1(:,:)

real, allocatable :: tair_accum_prev0(:,:)
real, allocatable :: tair_accum_prev1(:,:)

real, allocatable :: sh_v_accum_prev0(:,:)
real, allocatable :: sh_v_accum_prev1(:,:)

real, allocatable :: latheat_liq_accum_prev0(:,:)
real, allocatable :: latheat_liq_accum_prev1(:,:)

real, allocatable :: latheat_ice_accum_prev0(:,:)
real, allocatable :: latheat_ice_accum_prev1(:,:)

real, allocatable :: vels_l_accum_prev0(:)
real, allocatable :: vels_l_accum_prev1(:)

real, allocatable :: airtemp_l_accum_prev0(:)
real, allocatable :: airtemp_l_accum_prev1(:)
real, allocatable :: airtemp_l_accum_prev2(:)
real, allocatable :: airtemp_l_accum_prev3(:)

real, allocatable :: airshv_l_accum_prev0(:)
real, allocatable :: airshv_l_accum_prev1(:)

real, allocatable :: cantemp_l_accum_prev0(:)
real, allocatable :: cantemp_l_accum_prev1(:)
real, allocatable :: cantemp_l_accum_prev2(:)
real, allocatable :: cantemp_l_accum_prev3(:)

real, allocatable :: canshv_l_accum_prev0(:)
real, allocatable :: canshv_l_accum_prev1(:)

real, allocatable :: skintemp_l_accum_prev0(:)
real, allocatable :: skintemp_l_accum_prev1(:)
real, allocatable :: skintemp_l_accum_prev2(:)
real, allocatable :: skintemp_l_accum_prev3(:)

real, allocatable :: sfluxt_l_accum_prev0(:)
real, allocatable :: sfluxt_l_accum_prev1(:)

real, allocatable :: sfluxr_l_accum_prev0(:)
real, allocatable :: sfluxr_l_accum_prev1(:)

real, allocatable :: wxfer1_l_accum_prev0(:)
real, allocatable :: wxfer1_l_accum_prev1(:)

real, allocatable :: vels_s_accum_prev0(:)
real, allocatable :: vels_s_accum_prev1(:)

real, allocatable :: airtemp_s_accum_prev0(:)
real, allocatable :: airtemp_s_accum_prev1(:)
real, allocatable :: airtemp_s_accum_prev2(:)
real, allocatable :: airtemp_s_accum_prev3(:)

real, allocatable :: airshv_s_accum_prev0(:)
real, allocatable :: airshv_s_accum_prev1(:)

real, allocatable :: cantemp_s_accum_prev0(:)
real, allocatable :: cantemp_s_accum_prev1(:)
real, allocatable :: cantemp_s_accum_prev2(:)
real, allocatable :: cantemp_s_accum_prev3(:)

real, allocatable :: canshv_s_accum_prev0(:)
real, allocatable :: canshv_s_accum_prev1(:)

real, allocatable :: skintemp_s_accum_prev0(:)
real, allocatable :: skintemp_s_accum_prev1(:)
real, allocatable :: skintemp_s_accum_prev2(:)
real, allocatable :: skintemp_s_accum_prev3(:)

real, allocatable :: sfluxt_s_accum_prev0(:)
real, allocatable :: sfluxt_s_accum_prev1(:)

real, allocatable :: sfluxr_s_accum_prev0(:)
real, allocatable :: sfluxr_s_accum_prev1(:)

real(r8), allocatable ::  press_init(:,:)
real(r8), allocatable ::    rho_init(:,:)
real,     allocatable ::  theta_init(:,:)
real,     allocatable ::     vc_init(:,:)
real,     allocatable :: addsc1_init(:,:)

Contains

!===============================================================================

  subroutine alloc_plot()

  use mem_grid,  only: mza, mva, mwa
  use mem_basic, only: press, rho, theta, vc
  use mem_addsc, only: addsc
  use leaf_coms, only: mwl
  use sea_coms,  only: mws

  implicit none

  allocate ( accpmic_prev0(mwa)) ; accpmic_prev0(:) = 0.
  allocate ( accpmic_prev1(mwa)) ; accpmic_prev1(:) = 0.
  allocate ( accpmic_prev2(mwa)) ; accpmic_prev2(:) = 0.
  allocate ( accpmic_prev3(mwa)) ; accpmic_prev3(:) = 0.

  allocate ( accpcon_prev0(mwa)) ; accpcon_prev0(:) = 0.
  allocate ( accpcon_prev1(mwa)) ; accpcon_prev1(:) = 0.
  allocate ( accpcon_prev2(mwa)) ; accpcon_prev2(:) = 0.
  allocate ( accpcon_prev3(mwa)) ; accpcon_prev3(:) = 0.

  allocate ( rshort_accum_prev0(mwa)) ; rshort_accum_prev0(:) = 0.
  allocate ( rshort_accum_prev1(mwa)) ; rshort_accum_prev1(:) = 0.

  allocate ( rshortup_accum_prev0(mwa)) ; rshortup_accum_prev0(:) = 0.
  allocate ( rshortup_accum_prev1(mwa)) ; rshortup_accum_prev1(:) = 0.

  allocate ( rlong_accum_prev0(mwa)) ; rlong_accum_prev0(:) = 0.
  allocate ( rlong_accum_prev1(mwa)) ; rlong_accum_prev1(:) = 0.

  allocate ( rlongup_accum_prev0(mwa)) ; rlongup_accum_prev0(:) = 0.
  allocate ( rlongup_accum_prev1(mwa)) ; rlongup_accum_prev1(:) = 0.

  allocate ( rshort_top_accum_prev0(mwa)) ; rshort_top_accum_prev0(:) = 0.
  allocate ( rshort_top_accum_prev1(mwa)) ; rshort_top_accum_prev1(:) = 0.

  allocate ( rshortup_top_accum_prev0(mwa)) ; rshortup_top_accum_prev0(:) = 0.
  allocate ( rshortup_top_accum_prev1(mwa)) ; rshortup_top_accum_prev1(:) = 0.

  allocate ( rlongup_top_accum_prev0(mwa)) ; rlongup_top_accum_prev0(:) = 0.
  allocate ( rlongup_top_accum_prev1(mwa)) ; rlongup_top_accum_prev1(:) = 0.

  allocate ( rshort_clr_accum_prev0(mwa)) ; rshort_clr_accum_prev0(:) = 0.
  allocate ( rshort_clr_accum_prev1(mwa)) ; rshort_clr_accum_prev1(:) = 0.

  allocate ( rshortup_clr_accum_prev0(mwa)) ; rshortup_clr_accum_prev0(:) = 0.
  allocate ( rshortup_clr_accum_prev1(mwa)) ; rshortup_clr_accum_prev1(:) = 0.

  allocate ( rlong_clr_accum_prev0(mwa)) ; rlong_clr_accum_prev0(:) = 0.
  allocate ( rlong_clr_accum_prev1(mwa)) ; rlong_clr_accum_prev1(:) = 0.

  allocate ( rlongup_clr_accum_prev0(mwa)) ; rlongup_clr_accum_prev0(:) = 0.
  allocate ( rlongup_clr_accum_prev1(mwa)) ; rlongup_clr_accum_prev1(:) = 0.

  allocate ( rshort_top_clr_accum_prev0(mwa)) ; rshort_top_clr_accum_prev0(:) = 0.
  allocate ( rshort_top_clr_accum_prev1(mwa)) ; rshort_top_clr_accum_prev1(:) = 0.

  allocate ( rshortup_top_clr_accum_prev0(mwa)) ; rshortup_top_clr_accum_prev0(:) = 0.
  allocate ( rshortup_top_clr_accum_prev1(mwa)) ; rshortup_top_clr_accum_prev1(:) = 0.

  allocate ( rlongup_top_clr_accum_prev0(mwa)) ; rlongup_top_clr_accum_prev0(:) = 0.
  allocate ( rlongup_top_clr_accum_prev1(mwa)) ; rlongup_top_clr_accum_prev1(:) = 0.

  allocate ( sfluxt_accum_prev0(mwa)) ; sfluxt_accum_prev0(:) = 0.
  allocate ( sfluxt_accum_prev1(mwa)) ; sfluxt_accum_prev1(:) = 0.

  allocate ( sfluxr_accum_prev0(mwa)) ; sfluxr_accum_prev0(:) = 0.
  allocate ( sfluxr_accum_prev1(mwa)) ; sfluxr_accum_prev1(:) = 0.

  allocate ( vc_accum_prev0(mza,mva)) ; vc_accum_prev0(:,:) = 0.
  allocate ( vc_accum_prev1(mza,mva)) ; vc_accum_prev1(:,:) = 0.

  allocate ( wc_accum_prev0(mza,mwa)) ; wc_accum_prev0(:,:) = 0.
  allocate ( wc_accum_prev1(mza,mwa)) ; wc_accum_prev1(:,:) = 0.

  allocate ( press_accum_prev0(mza,mwa)) ; press_accum_prev0(:,:) = 0.
  allocate ( press_accum_prev1(mza,mwa)) ; press_accum_prev1(:,:) = 0.

  allocate ( tair_accum_prev0(mza,mwa)) ; tair_accum_prev0(:,:) = 0.
  allocate ( tair_accum_prev1(mza,mwa)) ; tair_accum_prev1(:,:) = 0.

  allocate ( sh_v_accum_prev0(mza,mwa)) ; sh_v_accum_prev0(:,:) = 0.
  allocate ( sh_v_accum_prev1(mza,mwa)) ; sh_v_accum_prev1(:,:) = 0.

  allocate ( latheat_liq_accum_prev0(mza,mwa)) ; latheat_liq_accum_prev0(:,:) = 0.
  allocate ( latheat_liq_accum_prev1(mza,mwa)) ; latheat_liq_accum_prev1(:,:) = 0.

  allocate ( latheat_ice_accum_prev0(mza,mwa)) ; latheat_ice_accum_prev0(:,:) = 0.
  allocate ( latheat_ice_accum_prev1(mza,mwa)) ; latheat_ice_accum_prev1(:,:) = 0.

  allocate ( vels_l_accum_prev0(mwl)) ; vels_l_accum_prev0(:) = 0.
  allocate ( vels_l_accum_prev1(mwl)) ; vels_l_accum_prev1(:) = 0.

  allocate ( airtemp_l_accum_prev0(mwl)) ; airtemp_l_accum_prev0(:) = 0.
  allocate ( airtemp_l_accum_prev1(mwl)) ; airtemp_l_accum_prev1(:) = 0.
  allocate ( airtemp_l_accum_prev2(mwl)) ; airtemp_l_accum_prev2(:) = 0.
  allocate ( airtemp_l_accum_prev3(mwl)) ; airtemp_l_accum_prev3(:) = 0.

  allocate ( airshv_l_accum_prev0(mwl)) ; airshv_l_accum_prev0(:) = 0.
  allocate ( airshv_l_accum_prev1(mwl)) ; airshv_l_accum_prev1(:) = 0.

  allocate ( cantemp_l_accum_prev0(mwl)) ; cantemp_l_accum_prev0(:) = 0.
  allocate ( cantemp_l_accum_prev1(mwl)) ; cantemp_l_accum_prev1(:) = 0.
  allocate ( cantemp_l_accum_prev2(mwl)) ; cantemp_l_accum_prev2(:) = 0.
  allocate ( cantemp_l_accum_prev3(mwl)) ; cantemp_l_accum_prev3(:) = 0.

  allocate ( canshv_l_accum_prev0(mwl)) ; canshv_l_accum_prev0(:) = 0.
  allocate ( canshv_l_accum_prev1(mwl)) ; canshv_l_accum_prev1(:) = 0.

  allocate ( skintemp_l_accum_prev0(mwl)) ; skintemp_l_accum_prev0(:) = 0.
  allocate ( skintemp_l_accum_prev1(mwl)) ; skintemp_l_accum_prev1(:) = 0.
  allocate ( skintemp_l_accum_prev2(mwl)) ; skintemp_l_accum_prev2(:) = 0.
  allocate ( skintemp_l_accum_prev3(mwl)) ; skintemp_l_accum_prev3(:) = 0.

  allocate ( sfluxt_l_accum_prev0(mwl)) ; sfluxt_l_accum_prev0(:) = 0.
  allocate ( sfluxt_l_accum_prev1(mwl)) ; sfluxt_l_accum_prev1(:) = 0.

  allocate ( sfluxr_l_accum_prev0(mwl)) ; sfluxr_l_accum_prev0(:) = 0.
  allocate ( sfluxr_l_accum_prev1(mwl)) ; sfluxr_l_accum_prev1(:) = 0.

  allocate ( wxfer1_l_accum_prev0(mwl)) ; wxfer1_l_accum_prev0(:) = 0.
  allocate ( wxfer1_l_accum_prev1(mwl)) ; wxfer1_l_accum_prev1(:) = 0.

  allocate ( vels_s_accum_prev0(mws)) ; vels_s_accum_prev0(:) = 0.
  allocate ( vels_s_accum_prev1(mws)) ; vels_s_accum_prev1(:) = 0.

  allocate ( airtemp_s_accum_prev0(mws)) ; airtemp_s_accum_prev0(:) = 0.
  allocate ( airtemp_s_accum_prev1(mws)) ; airtemp_s_accum_prev1(:) = 0.
  allocate ( airtemp_s_accum_prev2(mws)) ; airtemp_s_accum_prev2(:) = 0.
  allocate ( airtemp_s_accum_prev3(mws)) ; airtemp_s_accum_prev3(:) = 0.

  allocate ( airshv_s_accum_prev0(mws)) ; airshv_s_accum_prev0(:) = 0.
  allocate ( airshv_s_accum_prev1(mws)) ; airshv_s_accum_prev1(:) = 0.

  allocate ( cantemp_s_accum_prev0(mws)) ; cantemp_s_accum_prev0(:) = 0.
  allocate ( cantemp_s_accum_prev1(mws)) ; cantemp_s_accum_prev1(:) = 0.
  allocate ( cantemp_s_accum_prev2(mws)) ; cantemp_s_accum_prev2(:) = 0.
  allocate ( cantemp_s_accum_prev3(mws)) ; cantemp_s_accum_prev3(:) = 0.

  allocate ( canshv_s_accum_prev0(mws)) ; canshv_s_accum_prev0(:) = 0.
  allocate ( canshv_s_accum_prev1(mws)) ; canshv_s_accum_prev1(:) = 0.

  allocate ( skintemp_s_accum_prev0(mws)) ; skintemp_s_accum_prev0(:) = 0.
  allocate ( skintemp_s_accum_prev1(mws)) ; skintemp_s_accum_prev1(:) = 0.
  allocate ( skintemp_s_accum_prev2(mws)) ; skintemp_s_accum_prev2(:) = 0.
  allocate ( skintemp_s_accum_prev3(mws)) ; skintemp_s_accum_prev3(:) = 0.

  allocate ( sfluxt_s_accum_prev0(mws)) ; sfluxt_s_accum_prev0(:) = 0.
  allocate ( sfluxt_s_accum_prev1(mws)) ; sfluxt_s_accum_prev1(:) = 0.

  allocate ( sfluxr_s_accum_prev0(mws)) ; sfluxr_s_accum_prev0(:) = 0.
  allocate ( sfluxr_s_accum_prev1(mws)) ; sfluxr_s_accum_prev1(:) = 0.

  allocate ( press_init(mza,mwa)) ; press_init(:,:) = press(:,:)
  allocate (   rho_init(mza,mwa)) ; rho_init  (:,:) = rho  (:,:)
  allocate ( theta_init(mza,mwa)) ; theta_init(:,:) = theta(:,:)
  allocate (    vc_init(mza,mva)) ;    vc_init(:,:) =    vc(:,:)

! allocate (addsc1_init(mza,mwa)) ; addsc1_init(:,:) = addsc(1)%sclp(:,:)

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
  use mem_turb,       only: sfluxt, sfluxr
  use mem_flux_accum, only:    rshort_accum,         rshortup_accum, &
                                rlong_accum,          rlongup_accum, &
                           rshort_top_accum,     rshortup_top_accum, &
                          rlongup_top_accum,                         &
                           rshort_clr_accum,     rshortup_clr_accum, &
                            rlong_clr_accum,      rlongup_clr_accum, &
                       rshort_top_clr_accum, rshortup_top_clr_accum, &
                      rlongup_top_clr_accum,                         &
                               sfluxt_accum,   sfluxr_accum, &
                                   vc_accum,       wc_accum, &
                                press_accum,     tair_accum, &
                                 sh_v_accum,                 &
                          latheat_liq_accum, latheat_ice_accum, &
                               vels_l_accum,                 &
                            airtemp_l_accum, airshv_l_accum, &
                            cantemp_l_accum, canshv_l_accum, &
                           skintemp_l_accum,                 &
                             sfluxt_l_accum, sfluxr_l_accum, &
                             wxfer1_l_accum,                 &
                               vels_s_accum,                 &
                            airtemp_s_accum, airshv_s_accum, &
                            cantemp_s_accum, canshv_s_accum, &
                           skintemp_s_accum,                 &
                             sfluxt_s_accum, sfluxr_s_accum
  use mem_ijtabs, only: jtab_w, jtw_prog
  use mem_grid,   only: mza, lpw
  use mem_basic,  only: theta

  implicit none

  integer, intent(in) :: iplt_file ! File index for 'PLOTONLY' runs; 0 otherwise

  integer :: k,j,iw

! The following IF/ENDIF statements (but not the lines between them) should be
! commented out for most applications.  Uncommenting them provides the ability
! to sum precipitation values from multiple history files into a single
! accpmic_prev0 or accpcon_prev0 value, before later differencing between 
! prev0, prev1, prev2, and prev3 groupings.  For example, the second application
! of this capability was to sum the accumulated precipitation of 'Simulation A'
! each March 1 over 5 consecutive years and store the sum in the prev0 arrays,
! copy the prev0 sum into prev1 and then sum the June 1 precipitation totals
! for 'Simulation A' into prev0, copy the forgoing into prev2 and prev1,
! respectively, and sum March 1 totals for 'Simulation B' into prev0, and finaly
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

           rshort_clr_accum_prev1(:) =       rshort_clr_accum_prev0(:)
         rshortup_clr_accum_prev1(:) =     rshortup_clr_accum_prev0(:)
            rlong_clr_accum_prev1(:) =        rlong_clr_accum_prev0(:)
          rlongup_clr_accum_prev1(:) =      rlongup_clr_accum_prev0(:)
       rshort_top_clr_accum_prev1(:) =   rshort_top_clr_accum_prev0(:)
     rshortup_top_clr_accum_prev1(:) = rshortup_top_clr_accum_prev0(:)
      rlongup_top_clr_accum_prev1(:) =  rlongup_top_clr_accum_prev0(:)

           sfluxt_accum_prev1(:) =        sfluxt_accum_prev0(:)
           sfluxr_accum_prev1(:) =        sfluxr_accum_prev0(:)
             vc_accum_prev1(:,:) =          vc_accum_prev0(:,:)
             wc_accum_prev1(:,:) =          wc_accum_prev0(:,:)
          press_accum_prev1(:,:) =       press_accum_prev0(:,:)
           tair_accum_prev1(:,:) =        tair_accum_prev0(:,:)
           sh_v_accum_prev1(:,:) =        sh_v_accum_prev0(:,:)
    latheat_liq_accum_prev1(:,:) = latheat_liq_accum_prev0(:,:)
    latheat_ice_accum_prev1(:,:) = latheat_ice_accum_prev0(:,:)

       airtemp_l_accum_prev3(:) =  airtemp_l_accum_prev2(:)
       airtemp_l_accum_prev2(:) =  airtemp_l_accum_prev1(:)
       airtemp_l_accum_prev1(:) =  airtemp_l_accum_prev0(:)

       cantemp_l_accum_prev3(:) =  cantemp_l_accum_prev2(:)
       cantemp_l_accum_prev2(:) =  cantemp_l_accum_prev1(:)
       cantemp_l_accum_prev1(:) =  cantemp_l_accum_prev0(:)

      skintemp_l_accum_prev3(:) = skintemp_l_accum_prev2(:)
      skintemp_l_accum_prev2(:) = skintemp_l_accum_prev1(:)
      skintemp_l_accum_prev1(:) = skintemp_l_accum_prev0(:)

          vels_l_accum_prev1(:) =     vels_l_accum_prev0(:)
        airshv_l_accum_prev1(:) =   airshv_l_accum_prev0(:)
        canshv_l_accum_prev1(:) =   canshv_l_accum_prev0(:)
        sfluxt_l_accum_prev1(:) =   sfluxt_l_accum_prev0(:)
        sfluxr_l_accum_prev1(:) =   sfluxr_l_accum_prev0(:)
        wxfer1_l_accum_prev1(:) =   wxfer1_l_accum_prev0(:)

       airtemp_s_accum_prev3(:) =  airtemp_s_accum_prev2(:)
       airtemp_s_accum_prev2(:) =  airtemp_s_accum_prev1(:)
       airtemp_s_accum_prev1(:) =  airtemp_s_accum_prev0(:)

       cantemp_s_accum_prev3(:) =  cantemp_s_accum_prev2(:)
       cantemp_s_accum_prev2(:) =  cantemp_s_accum_prev1(:)
       cantemp_s_accum_prev1(:) =  cantemp_s_accum_prev0(:)

      skintemp_s_accum_prev3(:) = skintemp_s_accum_prev2(:)
      skintemp_s_accum_prev2(:) = skintemp_s_accum_prev1(:)
      skintemp_s_accum_prev1(:) = skintemp_s_accum_prev0(:)

          vels_s_accum_prev1(:) =     vels_s_accum_prev0(:)
        airshv_s_accum_prev1(:) =   airshv_s_accum_prev0(:)
        canshv_s_accum_prev1(:) =   canshv_s_accum_prev0(:)
        sfluxt_s_accum_prev1(:) =   sfluxt_s_accum_prev0(:)
        sfluxr_s_accum_prev1(:) =   sfluxr_s_accum_prev0(:)

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

           rshort_clr_accum_prev0(:) = 0.
         rshortup_clr_accum_prev0(:) = 0.
            rlong_clr_accum_prev0(:) = 0.
          rlongup_clr_accum_prev0(:) = 0.
       rshort_top_clr_accum_prev0(:) = 0.
     rshortup_top_clr_accum_prev0(:) = 0.
      rlongup_top_clr_accum_prev0(:) = 0.

           sfluxt_accum_prev0(:) = 0.
           sfluxr_accum_prev0(:) = 0.
             vc_accum_prev0(:,:) = 0.
             wc_accum_prev0(:,:) = 0.
          press_accum_prev0(:,:) = 0.
           tair_accum_prev0(:,:) = 0.
           sh_v_accum_prev0(:,:) = 0.
    latheat_liq_accum_prev0(:,:) = 0.
    latheat_ice_accum_prev0(:,:) = 0.

          vels_l_accum_prev0(:) = 0.
       airtemp_l_accum_prev0(:) = 0.
        airshv_l_accum_prev0(:) = 0.
       cantemp_l_accum_prev0(:) = 0.
        canshv_l_accum_prev0(:) = 0.
      skintemp_l_accum_prev0(:) = 0.
        sfluxt_l_accum_prev0(:) = 0.
        sfluxr_l_accum_prev0(:) = 0.
        wxfer1_l_accum_prev0(:) = 0.

          vels_s_accum_prev0(:) = 0.
       airtemp_s_accum_prev0(:) = 0.
        airshv_s_accum_prev0(:) = 0.
       cantemp_s_accum_prev0(:) = 0.
        canshv_s_accum_prev0(:) = 0.
      skintemp_s_accum_prev0(:) = 0.
        sfluxt_s_accum_prev0(:) = 0.
        sfluxr_s_accum_prev0(:) = 0.

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

  if (allocated(rshort_clr_accum)) rshort_clr_accum_prev0(:) = &
                                   rshort_clr_accum_prev0(:) + &
                              real(rshort_clr_accum(:))

  if (allocated(rshortup_clr_accum)) rshortup_clr_accum_prev0(:) = &
                                     rshortup_clr_accum_prev0(:) + &
                                real(rshortup_clr_accum(:))

  if (allocated(rlong_clr_accum)) rlong_clr_accum_prev0(:) = &
                                  rlong_clr_accum_prev0(:) + &
                             real(rlong_clr_accum(:))

  if (allocated(rlongup_clr_accum)) rlongup_clr_accum_prev0(:) = &
                                    rlongup_clr_accum_prev0(:) + &
                               real(rlongup_clr_accum(:))

  if (allocated(rshort_top_clr_accum)) rshort_top_clr_accum_prev0(:) = &
                                       rshort_top_clr_accum_prev0(:) + &
                                  real(rshort_top_clr_accum(:))

  if (allocated(rshortup_top_clr_accum)) rshortup_top_clr_accum_prev0(:) = &
                                         rshortup_top_clr_accum_prev0(:) + &
                                    real(rshortup_top_clr_accum(:))

  if (allocated(rlongup_top_clr_accum)) rlongup_top_clr_accum_prev0(:) = &
                                        rlongup_top_clr_accum_prev0(:) + &
                                   real(rlongup_top_clr_accum(:))

  if (allocated(sfluxt_accum)) sfluxt_accum_prev0(:) = &
                               sfluxt_accum_prev0(:) + &
                          real(sfluxt_accum(:))

  if (allocated(sfluxr_accum)) sfluxr_accum_prev0(:) = &
                               sfluxr_accum_prev0(:) + &
                          real(sfluxr_accum(:))

  if (allocated(vc_accum)) vc_accum_prev0(:,:) = &
                           vc_accum_prev0(:,:) + &
                      real(vc_accum(:,:))

  if (allocated(wc_accum)) wc_accum_prev0(:,:) = &
                           wc_accum_prev0(:,:) + &
                      real(wc_accum(:,:))

  if (allocated(press_accum)) press_accum_prev0(:,:) = &
                              press_accum_prev0(:,:) + &
                         real(press_accum(:,:))

  if (allocated(tair_accum)) tair_accum_prev0(:,:) = &
                             tair_accum_prev0(:,:) + &
                        real(tair_accum(:,:))

  if (allocated(sh_v_accum)) sh_v_accum_prev0(:,:) = &
                             sh_v_accum_prev0(:,:) + &
                        real(sh_v_accum(:,:))

  if (allocated(latheat_liq_accum)) latheat_liq_accum_prev0(:,:) = &
                                    latheat_liq_accum_prev0(:,:) + &
                               real(latheat_liq_accum(:,:))

  if (allocated(latheat_ice_accum)) latheat_ice_accum_prev0(:,:) = &
                                    latheat_ice_accum_prev0(:,:) + &
                               real(latheat_ice_accum(:,:))

  if (allocated(vels_l_accum))  vels_l_accum_prev0(:) = &
                                vels_l_accum_prev0(:) + &
                           real(vels_l_accum(:))

  if (allocated(airtemp_l_accum))  airtemp_l_accum_prev0(:) = &
                                   airtemp_l_accum_prev0(:) + &
                              real(airtemp_l_accum(:))

  if (allocated(airshv_l_accum))  airshv_l_accum_prev0(:) = &
                                  airshv_l_accum_prev0(:) + &
                             real(airshv_l_accum(:))

  if (allocated(cantemp_l_accum))  cantemp_l_accum_prev0(:) = &
                                   cantemp_l_accum_prev0(:) + &
                              real(cantemp_l_accum(:))

  if (allocated(canshv_l_accum))  canshv_l_accum_prev0(:) = &
                                  canshv_l_accum_prev0(:) + &
                             real(canshv_l_accum(:))

  if (allocated(skintemp_l_accum))  skintemp_l_accum_prev0(:) = &
                                    skintemp_l_accum_prev0(:) + &
                               real(skintemp_l_accum(:))

  if (allocated(sfluxt_l_accum))  sfluxt_l_accum_prev0(:) = &
                                  sfluxt_l_accum_prev0(:) + &
                             real(sfluxt_l_accum(:))

  if (allocated(sfluxr_l_accum))  sfluxr_l_accum_prev0(:) = &
                                  sfluxr_l_accum_prev0(:) + &
                             real(sfluxr_l_accum(:))

  if (allocated(wxfer1_l_accum))  wxfer1_l_accum_prev0(:) = &
                                  wxfer1_l_accum_prev0(:) + &
                             real(wxfer1_l_accum(:))

  if (allocated(vels_s_accum))  vels_s_accum_prev0(:) = &
                                vels_s_accum_prev0(:) + &
                           real(vels_s_accum(:))

  if (allocated(airtemp_s_accum))  airtemp_s_accum_prev0(:) = &
                                   airtemp_s_accum_prev0(:) + &
                              real(airtemp_s_accum(:))

  if (allocated(airshv_s_accum))  airshv_s_accum_prev0(:) = &
                                  airshv_s_accum_prev0(:) + &
                             real(airshv_s_accum(:))

  if (allocated(cantemp_s_accum))  cantemp_s_accum_prev0(:) = &
                                   cantemp_s_accum_prev0(:) + &
                              real(cantemp_s_accum(:))

  if (allocated(canshv_s_accum))  canshv_s_accum_prev0(:) = &
                                  canshv_s_accum_prev0(:) + &
                             real(canshv_s_accum(:))

  if (allocated(skintemp_s_accum))  skintemp_s_accum_prev0(:) = &
                                    skintemp_s_accum_prev0(:) + &
                               real(skintemp_s_accum(:))

  if (allocated(sfluxt_s_accum))  sfluxt_s_accum_prev0(:) = &
                                  sfluxt_s_accum_prev0(:) + &
                             real(sfluxt_s_accum(:))

  if (allocated(sfluxr_s_accum))  sfluxr_s_accum_prev0(:) = &
                                  sfluxr_s_accum_prev0(:) + &
                             real(sfluxr_s_accum(:))

  end subroutine copy_plot

End Module mem_plot
