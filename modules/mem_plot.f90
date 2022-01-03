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

  private :: r8

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

  real, allocatable :: vc_accum_prev0(:,:)
  real, allocatable :: vc_accum_prev1(:,:)

  real, allocatable :: wc_accum_prev0(:,:)
  real, allocatable :: wc_accum_prev1(:,:)

  real, allocatable :: press_accum_prev0(:,:)
  real, allocatable :: press_accum_prev1(:,:)

  real, allocatable :: tair_accum_prev0(:,:)
  real, allocatable :: tair_accum_prev1(:,:)

  real, allocatable :: rr_v_accum_prev0(:,:)
  real, allocatable :: rr_v_accum_prev1(:,:)

  real, allocatable :: latheat_liq_accum_prev0(:,:)
  real, allocatable :: latheat_liq_accum_prev1(:,:)

  real, allocatable :: latheat_ice_accum_prev0(:,:)
  real, allocatable :: latheat_ice_accum_prev1(:,:)

  real, allocatable :: vels_accum_prev0(:)
  real, allocatable :: vels_accum_prev1(:)
  real, allocatable :: vels_accum_prev2(:)
  real, allocatable :: vels_accum_prev3(:)

  real, allocatable :: airtemp_accum_prev0(:)
  real, allocatable :: airtemp_accum_prev1(:)
  real, allocatable :: airtemp_accum_prev2(:)
  real, allocatable :: airtemp_accum_prev3(:)

  real, allocatable :: airrrv_accum_prev0(:)
  real, allocatable :: airrrv_accum_prev1(:)
  real, allocatable :: airrrv_accum_prev2(:)
  real, allocatable :: airrrv_accum_prev3(:)

  real, allocatable :: cantemp_accum_prev0(:)
  real, allocatable :: cantemp_accum_prev1(:)
  real, allocatable :: cantemp_accum_prev2(:)
  real, allocatable :: cantemp_accum_prev3(:)

  real, allocatable :: canrrv_accum_prev0(:)
  real, allocatable :: canrrv_accum_prev1(:)
  real, allocatable :: canrrv_accum_prev2(:)
  real, allocatable :: canrrv_accum_prev3(:)

  real, allocatable :: skintemp_accum_prev0(:)
  real, allocatable :: skintemp_accum_prev1(:)
  real, allocatable :: skintemp_accum_prev2(:)
  real, allocatable :: skintemp_accum_prev3(:)

  real, allocatable :: sfluxt_accum_prev0(:)
  real, allocatable :: sfluxt_accum_prev1(:)
  real, allocatable :: sfluxt_accum_prev2(:)
  real, allocatable :: sfluxt_accum_prev3(:)

  real, allocatable :: sfluxr_accum_prev0(:) ! fast can nud
  real, allocatable :: sfluxr_accum_prev1(:) ! fast can nud
  real, allocatable :: sfluxr_accum_prev2(:)
  real, allocatable :: sfluxr_accum_prev3(:)

  real, allocatable :: pcp_accum_prev0(:) ! fast can nud
  real, allocatable :: pcp_accum_prev1(:) ! fast can nud
  real, allocatable :: pcp_accum_prev2(:) ! fast can nud
  real, allocatable :: pcp_accum_prev3(:) ! fast can nud

  real, allocatable :: runoff_accum_prev0(:)
  real, allocatable :: runoff_accum_prev1(:)
  real, allocatable :: runoff_accum_prev2(:)
  real, allocatable :: runoff_accum_prev3(:)

  real, allocatable :: sfctemp_accum_prev0(:) ! fast can nud
  real, allocatable :: sfctemp_accum_prev1(:) ! fast can nud

  real, allocatable :: fracliq_accum_prev0(:) ! fast can nud
  real, allocatable :: fracliq_accum_prev1(:) ! fast can nud

  real, allocatable :: wxferi_accum_prev0(:)
  real, allocatable :: wxferi_accum_prev1(:)
  real, allocatable :: wxferi_accum_prev2(:)
  real, allocatable :: wxferi_accum_prev3(:)

  real, allocatable :: wxferp_accum_prev0(:)
  real, allocatable :: wxferp_accum_prev1(:)
  real, allocatable :: wxferp_accum_prev2(:)
  real, allocatable :: wxferp_accum_prev3(:)

  real, allocatable :: wxfer1_accum_prev0(:)
  real, allocatable :: wxfer1_accum_prev1(:)
  real, allocatable :: wxfer1_accum_prev2(:)
  real, allocatable :: wxfer1_accum_prev3(:)

  real, allocatable :: airtemp_dmin_accum_prev0(:)
  real, allocatable :: airtemp_dmin_accum_prev1(:)
  real, allocatable :: airtemp_dmin_accum_prev2(:)
  real, allocatable :: airtemp_dmin_accum_prev3(:)

  real, allocatable :: airtemp_dmax_accum_prev0(:)
  real, allocatable :: airtemp_dmax_accum_prev1(:)
  real, allocatable :: airtemp_dmax_accum_prev2(:)
  real, allocatable :: airtemp_dmax_accum_prev3(:)

  real, allocatable :: cantemp_dmin_accum_prev0(:)
  real, allocatable :: cantemp_dmin_accum_prev1(:)
  real, allocatable :: cantemp_dmin_accum_prev2(:)
  real, allocatable :: cantemp_dmin_accum_prev3(:)

  real, allocatable :: cantemp_dmax_accum_prev0(:)
  real, allocatable :: cantemp_dmax_accum_prev1(:)
  real, allocatable :: cantemp_dmax_accum_prev2(:)
  real, allocatable :: cantemp_dmax_accum_prev3(:)

  real, allocatable :: vegtemp_dmin_accum_prev0(:)
  real, allocatable :: vegtemp_dmin_accum_prev1(:)
  real, allocatable :: vegtemp_dmin_accum_prev2(:)
  real, allocatable :: vegtemp_dmin_accum_prev3(:)

  real, allocatable :: vegtemp_dmax_accum_prev0(:)
  real, allocatable :: vegtemp_dmax_accum_prev1(:)
  real, allocatable :: vegtemp_dmax_accum_prev2(:)
  real, allocatable :: vegtemp_dmax_accum_prev3(:)

  real, allocatable :: soiltemp_dmin_accum_prev0(:)
  real, allocatable :: soiltemp_dmin_accum_prev1(:)
  real, allocatable :: soiltemp_dmin_accum_prev2(:)
  real, allocatable :: soiltemp_dmin_accum_prev3(:)

  real, allocatable :: soiltemp_dmax_accum_prev0(:)
  real, allocatable :: soiltemp_dmax_accum_prev1(:)
  real, allocatable :: soiltemp_dmax_accum_prev2(:)
  real, allocatable :: soiltemp_dmax_accum_prev3(:)

!------------------------------------------------------------------------------
! Exception for mem_plot: Not an accumulated quantity but an instantaneous one.
! Perhaps in future add to mem_flux_accum arrays.

  real, allocatable :: soil_water_tot_prev0(:)
  real, allocatable :: soil_water_tot_prev1(:)

  real, allocatable :: head_wtab_prev0(:)
  real, allocatable :: head_wtab_prev1(:)
!------------------------------------------------------------------------------

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
!   use mem_addsc, only: addsc
    use mem_sfcg,  only: mwsfc
    use mem_land,  only: mland

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

    allocate ( vc_accum_prev0(mza,mva)) ; vc_accum_prev0(:,:) = 0.
    allocate ( vc_accum_prev1(mza,mva)) ; vc_accum_prev1(:,:) = 0.

    allocate ( wc_accum_prev0(mza,mwa)) ; wc_accum_prev0(:,:) = 0.
    allocate ( wc_accum_prev1(mza,mwa)) ; wc_accum_prev1(:,:) = 0.

    allocate ( press_accum_prev0(mza,mwa)) ; press_accum_prev0(:,:) = 0.
    allocate ( press_accum_prev1(mza,mwa)) ; press_accum_prev1(:,:) = 0.

    allocate ( tair_accum_prev0(mza,mwa)) ; tair_accum_prev0(:,:) = 0.
    allocate ( tair_accum_prev1(mza,mwa)) ; tair_accum_prev1(:,:) = 0.

    allocate ( rr_v_accum_prev0(mza,mwa)) ; rr_v_accum_prev0(:,:) = 0.
    allocate ( rr_v_accum_prev1(mza,mwa)) ; rr_v_accum_prev1(:,:) = 0.

    allocate ( latheat_liq_accum_prev0(mza,mwa)) ; latheat_liq_accum_prev0(:,:) = 0.
    allocate ( latheat_liq_accum_prev1(mza,mwa)) ; latheat_liq_accum_prev1(:,:) = 0.

    allocate ( latheat_ice_accum_prev0(mza,mwa)) ; latheat_ice_accum_prev0(:,:) = 0.
    allocate ( latheat_ice_accum_prev1(mza,mwa)) ; latheat_ice_accum_prev1(:,:) = 0.

    allocate ( vels_accum_prev0(mwsfc)) ; vels_accum_prev0(:) = 0.
    allocate ( vels_accum_prev1(mwsfc)) ; vels_accum_prev1(:) = 0.
    allocate ( vels_accum_prev2(mwsfc)) ; vels_accum_prev2(:) = 0.
    allocate ( vels_accum_prev3(mwsfc)) ; vels_accum_prev3(:) = 0.

    allocate ( airtemp_accum_prev0(mwsfc)) ; airtemp_accum_prev0(:) = 0.
    allocate ( airtemp_accum_prev1(mwsfc)) ; airtemp_accum_prev1(:) = 0.
    allocate ( airtemp_accum_prev2(mwsfc)) ; airtemp_accum_prev2(:) = 0.
    allocate ( airtemp_accum_prev3(mwsfc)) ; airtemp_accum_prev3(:) = 0.

    allocate ( airrrv_accum_prev0(mwsfc)) ; airrrv_accum_prev0(:) = 0.
    allocate ( airrrv_accum_prev1(mwsfc)) ; airrrv_accum_prev1(:) = 0.
    allocate ( airrrv_accum_prev2(mwsfc)) ; airrrv_accum_prev2(:) = 0.
    allocate ( airrrv_accum_prev3(mwsfc)) ; airrrv_accum_prev3(:) = 0.

    allocate ( cantemp_accum_prev0(mwsfc)) ; cantemp_accum_prev0(:) = 0.
    allocate ( cantemp_accum_prev1(mwsfc)) ; cantemp_accum_prev1(:) = 0.
    allocate ( cantemp_accum_prev2(mwsfc)) ; cantemp_accum_prev2(:) = 0.
    allocate ( cantemp_accum_prev3(mwsfc)) ; cantemp_accum_prev3(:) = 0.

    allocate ( canrrv_accum_prev0(mwsfc)) ; canrrv_accum_prev0(:) = 0.
    allocate ( canrrv_accum_prev1(mwsfc)) ; canrrv_accum_prev1(:) = 0.
    allocate ( canrrv_accum_prev2(mwsfc)) ; canrrv_accum_prev2(:) = 0.
    allocate ( canrrv_accum_prev3(mwsfc)) ; canrrv_accum_prev3(:) = 0.

    allocate ( skintemp_accum_prev0(mwsfc)) ; skintemp_accum_prev0(:) = 0.
    allocate ( skintemp_accum_prev1(mwsfc)) ; skintemp_accum_prev1(:) = 0.
    allocate ( skintemp_accum_prev2(mwsfc)) ; skintemp_accum_prev2(:) = 0.
    allocate ( skintemp_accum_prev3(mwsfc)) ; skintemp_accum_prev3(:) = 0.

    allocate ( sfluxt_accum_prev0(mwsfc)) ; sfluxt_accum_prev0(:) = 0.
    allocate ( sfluxt_accum_prev1(mwsfc)) ; sfluxt_accum_prev1(:) = 0.
    allocate ( sfluxt_accum_prev2(mwsfc)) ; sfluxt_accum_prev2(:) = 0.
    allocate ( sfluxt_accum_prev3(mwsfc)) ; sfluxt_accum_prev3(:) = 0.

    allocate ( sfluxr_accum_prev0(mwsfc)) ; sfluxr_accum_prev0(:) = 0. ! fast can nud
    allocate ( sfluxr_accum_prev1(mwsfc)) ; sfluxr_accum_prev1(:) = 0. ! fast can nud
    allocate ( sfluxr_accum_prev2(mwsfc)) ; sfluxr_accum_prev2(:) = 0.
    allocate ( sfluxr_accum_prev3(mwsfc)) ; sfluxr_accum_prev3(:) = 0.

    allocate ( pcp_accum_prev0(mwsfc)) ; pcp_accum_prev0(:) = 0. ! fast can nud
    allocate ( pcp_accum_prev1(mwsfc)) ; pcp_accum_prev1(:) = 0. ! fast can nud
    allocate ( pcp_accum_prev2(mwsfc)) ; pcp_accum_prev2(:) = 0. ! fast can nud
    allocate ( pcp_accum_prev3(mwsfc)) ; pcp_accum_prev3(:) = 0. ! fast can nud

    allocate ( runoff_accum_prev0(mwsfc)) ; runoff_accum_prev0(:) = 0.
    allocate ( runoff_accum_prev1(mwsfc)) ; runoff_accum_prev1(:) = 0.
    allocate ( runoff_accum_prev2(mwsfc)) ; runoff_accum_prev2(:) = 0.
    allocate ( runoff_accum_prev3(mwsfc)) ; runoff_accum_prev3(:) = 0.

    allocate ( sfctemp_accum_prev0(mwsfc)) ; sfctemp_accum_prev0(:) = 0. ! fast can nud
    allocate ( sfctemp_accum_prev1(mwsfc)) ; sfctemp_accum_prev1(:) = 0. ! fast can nud

    allocate ( fracliq_accum_prev0(mwsfc)) ; fracliq_accum_prev0(:) = 0. ! fast can nud
    allocate ( fracliq_accum_prev1(mwsfc)) ; fracliq_accum_prev1(:) = 0. ! fast can nud

    allocate ( wxferi_accum_prev0(mland)) ; wxferi_accum_prev0(:) = 0.
    allocate ( wxferi_accum_prev1(mland)) ; wxferi_accum_prev1(:) = 0.
    allocate ( wxferi_accum_prev2(mland)) ; wxferi_accum_prev2(:) = 0.
    allocate ( wxferi_accum_prev3(mland)) ; wxferi_accum_prev3(:) = 0.

    allocate ( wxferp_accum_prev0(mland)) ; wxferp_accum_prev0(:) = 0.
    allocate ( wxferp_accum_prev1(mland)) ; wxferp_accum_prev1(:) = 0.
    allocate ( wxferp_accum_prev2(mland)) ; wxferp_accum_prev2(:) = 0.
    allocate ( wxferp_accum_prev3(mland)) ; wxferp_accum_prev3(:) = 0.

    allocate ( wxfer1_accum_prev0(mland)) ; wxfer1_accum_prev0(:) = 0.
    allocate ( wxfer1_accum_prev1(mland)) ; wxfer1_accum_prev1(:) = 0.
    allocate ( wxfer1_accum_prev2(mland)) ; wxfer1_accum_prev2(:) = 0.
    allocate ( wxfer1_accum_prev3(mland)) ; wxfer1_accum_prev3(:) = 0.

    allocate ( airtemp_dmin_accum_prev0(mwsfc)) ; airtemp_dmin_accum_prev0(:) = 0.
    allocate ( airtemp_dmin_accum_prev1(mwsfc)) ; airtemp_dmin_accum_prev1(:) = 0.
    allocate ( airtemp_dmin_accum_prev2(mwsfc)) ; airtemp_dmin_accum_prev2(:) = 0.
    allocate ( airtemp_dmin_accum_prev3(mwsfc)) ; airtemp_dmin_accum_prev3(:) = 0.

    allocate ( airtemp_dmax_accum_prev0(mwsfc)) ; airtemp_dmax_accum_prev0(:) = 0.
    allocate ( airtemp_dmax_accum_prev1(mwsfc)) ; airtemp_dmax_accum_prev1(:) = 0.
    allocate ( airtemp_dmax_accum_prev2(mwsfc)) ; airtemp_dmax_accum_prev2(:) = 0.
    allocate ( airtemp_dmax_accum_prev3(mwsfc)) ; airtemp_dmax_accum_prev3(:) = 0.

    allocate ( cantemp_dmin_accum_prev0(mwsfc)) ; cantemp_dmin_accum_prev0(:) = 0.
    allocate ( cantemp_dmin_accum_prev1(mwsfc)) ; cantemp_dmin_accum_prev1(:) = 0.
    allocate ( cantemp_dmin_accum_prev2(mwsfc)) ; cantemp_dmin_accum_prev2(:) = 0.
    allocate ( cantemp_dmin_accum_prev3(mwsfc)) ; cantemp_dmin_accum_prev3(:) = 0.

    allocate ( cantemp_dmax_accum_prev0(mwsfc)) ; cantemp_dmax_accum_prev0(:) = 0.
    allocate ( cantemp_dmax_accum_prev1(mwsfc)) ; cantemp_dmax_accum_prev1(:) = 0.
    allocate ( cantemp_dmax_accum_prev2(mwsfc)) ; cantemp_dmax_accum_prev2(:) = 0.
    allocate ( cantemp_dmax_accum_prev3(mwsfc)) ; cantemp_dmax_accum_prev3(:) = 0.

    allocate ( vegtemp_dmin_accum_prev0(mland)) ; vegtemp_dmin_accum_prev0(:) = 0.
    allocate ( vegtemp_dmin_accum_prev1(mland)) ; vegtemp_dmin_accum_prev1(:) = 0.
    allocate ( vegtemp_dmin_accum_prev2(mland)) ; vegtemp_dmin_accum_prev2(:) = 0.
    allocate ( vegtemp_dmin_accum_prev3(mland)) ; vegtemp_dmin_accum_prev3(:) = 0.

    allocate ( vegtemp_dmax_accum_prev0(mland)) ; vegtemp_dmax_accum_prev0(:) = 0.
    allocate ( vegtemp_dmax_accum_prev1(mland)) ; vegtemp_dmax_accum_prev1(:) = 0.
    allocate ( vegtemp_dmax_accum_prev2(mland)) ; vegtemp_dmax_accum_prev2(:) = 0.
    allocate ( vegtemp_dmax_accum_prev3(mland)) ; vegtemp_dmax_accum_prev3(:) = 0.

    allocate ( soiltemp_dmin_accum_prev0(mland)) ; soiltemp_dmin_accum_prev0(:) = 0.
    allocate ( soiltemp_dmin_accum_prev1(mland)) ; soiltemp_dmin_accum_prev1(:) = 0.
    allocate ( soiltemp_dmin_accum_prev2(mland)) ; soiltemp_dmin_accum_prev2(:) = 0.
    allocate ( soiltemp_dmin_accum_prev3(mland)) ; soiltemp_dmin_accum_prev3(:) = 0.

    allocate ( soiltemp_dmax_accum_prev0(mland)) ; soiltemp_dmax_accum_prev0(:) = 0.
    allocate ( soiltemp_dmax_accum_prev1(mland)) ; soiltemp_dmax_accum_prev1(:) = 0.
    allocate ( soiltemp_dmax_accum_prev2(mland)) ; soiltemp_dmax_accum_prev2(:) = 0.
    allocate ( soiltemp_dmax_accum_prev3(mland)) ; soiltemp_dmax_accum_prev3(:) = 0.

!------------------------------------------------------------------------------
! Exception for mem_plot: Not an accumulated quantity but an instantaneous one.
! Perhaps in future add to mem_flux_accum arrays.

    allocate ( soil_water_tot_prev0(mland)) ; soil_water_tot_prev0(:) = 0.
    allocate ( soil_water_tot_prev1(mland)) ; soil_water_tot_prev1(:) = 0.

    allocate ( head_wtab_prev0(mwsfc)) ; head_wtab_prev0(:) = 0.
    allocate ( head_wtab_prev1(mwsfc)) ; head_wtab_prev1(:) = 0.
!------------------------------------------------------------------------------

    allocate ( press_init(mza,mwa)) ; press_init(:,:) = press(:,:)
    allocate (   rho_init(mza,mwa)) ; rho_init  (:,:) = rho  (:,:)
    allocate ( theta_init(mza,mwa)) ; theta_init(:,:) = theta(:,:)
    allocate (    vc_init(mza,mva)) ; vc_init(:,:) = vc(:,:)
!   allocate (addsc1_init(mza,mwa)) ; addsc1_init(:,:) = addsc(1)%sclp(:,:)

    time8_prev0 = 0._r8
    time8_prev1 = 0._r8
    time8_prev2 = 0._r8
    time8_prev3 = 0._r8

  end subroutine alloc_plot

!===============================================================================

  subroutine copy_plot(iplt_file)

    use mem_cuparm,     only: aconpr
    use mem_micro,      only: accpd, accpr, accpp, accps, accpa, accpg, accph
    use misc_coms,      only: time8
    use mem_sfcg,       only: sfcg, mwsfc, itab_wsfc
    use mem_land,       only: land, mland, omland, nzg, dslz, slzt
    use leaf4_soil,     only: soil_wat2pot
    use oname_coms,     only: nl
    use mem_para,       only: myrank
    use mem_flux_accum, only: rshort_accum,         rshortup_accum, &
                               rlong_accum,          rlongup_accum, &
                          rshort_top_accum,     rshortup_top_accum, &
                         rlongup_top_accum,                         &
                          rshort_clr_accum,     rshortup_clr_accum, &
                           rlong_clr_accum,      rlongup_clr_accum, &
                      rshort_top_clr_accum, rshortup_top_clr_accum, &
                     rlongup_top_clr_accum,                         &
                              sfluxt_accum,           sfluxr_accum, & ! ---/nud
                                 pcp_accum,           runoff_accum, & ! nud/---
                             sfctemp_accum,          fracliq_accum, & ! nud/nud
                                  vc_accum,               wc_accum, &
                               press_accum,             tair_accum, &
                                rr_v_accum,             vels_accum, &
                         latheat_liq_accum,      latheat_ice_accum, &
                             airtemp_accum,           airrrv_accum, &
                             cantemp_accum,           canrrv_accum, &
                            skintemp_accum,           wxferi_accum, &
                              wxferp_accum,           wxfer1_accum, &
                        airtemp_dmin_accum,     airtemp_dmax_accum, &
                        cantemp_dmin_accum,     cantemp_dmax_accum, &
                        vegtemp_dmin_accum,     vegtemp_dmax_accum, &
                       soiltemp_dmin_accum,    soiltemp_dmax_accum

    implicit none

    integer, intent(in) :: iplt_file ! File index for 'PLOTONLY' runs; 0 otherwise
    integer :: k, iwsfc, iland, klev

    real :: head(nzg)
    real :: psi, psi_slope

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

 ! if (iplt_file ==  4 .or. &
 !     iplt_file ==  7 .or. &
 !     iplt_file == 10 .or. &
 !     iplt_file == 13 .or. &
 !     iplt_file == 16 .or. &
 !     iplt_file == 19 .or. &
 !     iplt_file == 22 .or. &
 !     iplt_file == 25 .or. &
 !     iplt_file == 28 .or. &
 !     iplt_file == 31 .or. &
 !     iplt_file == 34 .or. &
 !     iplt_file == 37) then

!   if (iplt_file == 121) then

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

             vc_accum_prev1(:,:) =          vc_accum_prev0(:,:)
             wc_accum_prev1(:,:) =          wc_accum_prev0(:,:)
          press_accum_prev1(:,:) =       press_accum_prev0(:,:)
           tair_accum_prev1(:,:) =        tair_accum_prev0(:,:)
           rr_v_accum_prev1(:,:) =        rr_v_accum_prev0(:,:)
    latheat_liq_accum_prev1(:,:) = latheat_liq_accum_prev0(:,:)
    latheat_ice_accum_prev1(:,:) = latheat_ice_accum_prev0(:,:)

    airtemp_accum_prev3(:) = airtemp_accum_prev2(:)
    airtemp_accum_prev2(:) = airtemp_accum_prev1(:)
    airtemp_accum_prev1(:) = airtemp_accum_prev0(:)

    cantemp_accum_prev3(:) = cantemp_accum_prev2(:)
    cantemp_accum_prev2(:) = cantemp_accum_prev1(:)
    cantemp_accum_prev1(:) = cantemp_accum_prev0(:)

    skintemp_accum_prev3(:) = skintemp_accum_prev2(:)
    skintemp_accum_prev2(:) = skintemp_accum_prev1(:)
    skintemp_accum_prev1(:) = skintemp_accum_prev0(:)

    vels_accum_prev3(:) = vels_accum_prev2(:)
    vels_accum_prev2(:) = vels_accum_prev1(:)
    vels_accum_prev1(:) = vels_accum_prev0(:)

    airrrv_accum_prev3(:) = airrrv_accum_prev2(:)
    airrrv_accum_prev2(:) = airrrv_accum_prev1(:)
    airrrv_accum_prev1(:) = airrrv_accum_prev0(:)

    canrrv_accum_prev3(:) = canrrv_accum_prev2(:)
    canrrv_accum_prev2(:) = canrrv_accum_prev1(:)
    canrrv_accum_prev1(:) = canrrv_accum_prev0(:)

    sfluxt_accum_prev3(:) = sfluxt_accum_prev2(:)
    sfluxt_accum_prev2(:) = sfluxt_accum_prev1(:)
    sfluxt_accum_prev1(:) = sfluxt_accum_prev0(:)

    sfluxr_accum_prev3(:) = sfluxr_accum_prev2(:)
    sfluxr_accum_prev2(:) = sfluxr_accum_prev1(:)
    sfluxr_accum_prev1(:) = sfluxr_accum_prev0(:) ! fast can nud

    pcp_accum_prev3(:) = pcp_accum_prev2(:) ! fast can nud
    pcp_accum_prev2(:) = pcp_accum_prev1(:) ! fast can nud
    pcp_accum_prev1(:) = pcp_accum_prev0(:) ! fast can nud

    runoff_accum_prev3(:) = runoff_accum_prev2(:)
    runoff_accum_prev2(:) = runoff_accum_prev1(:)
    runoff_accum_prev1(:) = runoff_accum_prev0(:)

    sfctemp_accum_prev1(:) = sfctemp_accum_prev0(:) ! fast can nud
    fracliq_accum_prev1(:) = fracliq_accum_prev0(:) ! fast can nud

    wxferi_accum_prev3(:) = wxferi_accum_prev2(:)
    wxferi_accum_prev2(:) = wxferi_accum_prev1(:)
    wxferi_accum_prev1(:) = wxferi_accum_prev0(:)

    wxferp_accum_prev3(:) = wxferp_accum_prev2(:)
    wxferp_accum_prev2(:) = wxferp_accum_prev1(:)
    wxferp_accum_prev1(:) = wxferp_accum_prev0(:)

    wxfer1_accum_prev3(:) = wxfer1_accum_prev2(:)
    wxfer1_accum_prev2(:) = wxfer1_accum_prev1(:)
    wxfer1_accum_prev1(:) = wxfer1_accum_prev0(:)

    airtemp_dmin_accum_prev3(:) = airtemp_dmin_accum_prev2(:)
    airtemp_dmin_accum_prev2(:) = airtemp_dmin_accum_prev1(:)
    airtemp_dmin_accum_prev1(:) = airtemp_dmin_accum_prev0(:)

    airtemp_dmax_accum_prev3(:) = airtemp_dmax_accum_prev2(:)
    airtemp_dmax_accum_prev2(:) = airtemp_dmax_accum_prev1(:)
    airtemp_dmax_accum_prev1(:) = airtemp_dmax_accum_prev0(:)

    cantemp_dmin_accum_prev3(:) = cantemp_dmin_accum_prev2(:)
    cantemp_dmin_accum_prev2(:) = cantemp_dmin_accum_prev1(:)
    cantemp_dmin_accum_prev1(:) = cantemp_dmin_accum_prev0(:)

    cantemp_dmax_accum_prev3(:) = cantemp_dmax_accum_prev2(:)
    cantemp_dmax_accum_prev2(:) = cantemp_dmax_accum_prev1(:)
    cantemp_dmax_accum_prev1(:) = cantemp_dmax_accum_prev0(:)

    vegtemp_dmin_accum_prev3(:) = vegtemp_dmin_accum_prev2(:)
    vegtemp_dmin_accum_prev2(:) = vegtemp_dmin_accum_prev1(:)
    vegtemp_dmin_accum_prev1(:) = vegtemp_dmin_accum_prev0(:)

    vegtemp_dmax_accum_prev3(:) = vegtemp_dmax_accum_prev2(:)
    vegtemp_dmax_accum_prev2(:) = vegtemp_dmax_accum_prev1(:)
    vegtemp_dmax_accum_prev1(:) = vegtemp_dmax_accum_prev0(:)

    soiltemp_dmin_accum_prev3(:) = soiltemp_dmin_accum_prev2(:)
    soiltemp_dmin_accum_prev2(:) = soiltemp_dmin_accum_prev1(:)
    soiltemp_dmin_accum_prev1(:) = soiltemp_dmin_accum_prev0(:)

    soiltemp_dmax_accum_prev3(:) = soiltemp_dmax_accum_prev2(:)
    soiltemp_dmax_accum_prev2(:) = soiltemp_dmax_accum_prev1(:)
    soiltemp_dmax_accum_prev1(:) = soiltemp_dmax_accum_prev0(:)

!------------------------------------------------------------------------------
! Exception for mem_plot: Not an accumulated quantity but an instantaneous one.
! Perhaps in future add to mem_flux_accum arrays.

    soil_water_tot_prev1(:) = soil_water_tot_prev0(:)

    head_wtab_prev1(:) = head_wtab_prev0(:)
!------------------------------------------------------------------------------

! Fill most recent previous values

    time8_prev0 = 0._r8

    accpmic_prev0(:) = 0.
    accpcon_prev0(:) = 0.

          rshort_accum_prev0(:) = 0.
        rshortup_accum_prev0(:) = 0.
           rlong_accum_prev0(:) = 0.
         rlongup_accum_prev0(:) = 0.
      rshort_top_accum_prev0(:) = 0.
    rshortup_top_accum_prev0(:) = 0.
     rlongup_top_accum_prev0(:) = 0.

             vc_accum_prev0(:,:) = 0.
             wc_accum_prev0(:,:) = 0.
          press_accum_prev0(:,:) = 0.
           tair_accum_prev0(:,:) = 0.
           rr_v_accum_prev0(:,:) = 0.
    latheat_liq_accum_prev0(:,:) = 0.
    latheat_ice_accum_prev0(:,:) = 0.

        vels_accum_prev0(:) = 0.
     airtemp_accum_prev0(:) = 0.
      airrrv_accum_prev0(:) = 0.
     cantemp_accum_prev0(:) = 0.
      canrrv_accum_prev0(:) = 0.
    skintemp_accum_prev0(:) = 0.
      sfluxt_accum_prev0(:) = 0.
      sfluxr_accum_prev0(:) = 0. ! fast can nud
         pcp_accum_prev0(:) = 0. ! fast can nud
      runoff_accum_prev0(:) = 0.
     sfctemp_accum_prev0(:) = 0. ! fast can nud
     fracliq_accum_prev0(:) = 0. ! fast can nud
      wxferi_accum_prev0(:) = 0.
      wxferp_accum_prev0(:) = 0.
      wxfer1_accum_prev0(:) = 0.

     airtemp_dmin_accum_prev0(:) = 0.
     airtemp_dmax_accum_prev0(:) = 0.
     cantemp_dmin_accum_prev0(:) = 0.
     cantemp_dmax_accum_prev0(:) = 0.
     vegtemp_dmin_accum_prev0(:) = 0.
     vegtemp_dmax_accum_prev0(:) = 0.
    soiltemp_dmin_accum_prev0(:) = 0.
    soiltemp_dmax_accum_prev0(:) = 0.

!------------------------------------------------------------------------------
! Exception for mem_plot: Not an accumulated quantity but an instantaneous one.
! Perhaps in future add to mem_flux_accum arrays.

    soil_water_tot_prev0(:) = 0.
    head_wtab_prev0     (:) = 0.
!------------------------------------------------------------------------------

 !   endif

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

    if (allocated(rr_v_accum)) rr_v_accum_prev0(:,:) = &
                               rr_v_accum_prev0(:,:) + &
                          real(rr_v_accum(:,:))

    if (allocated(latheat_liq_accum)) latheat_liq_accum_prev0(:,:) = &
                                      latheat_liq_accum_prev0(:,:) + &
                                 real(latheat_liq_accum(:,:))

    if (allocated(latheat_ice_accum)) latheat_ice_accum_prev0(:,:) = &
                                      latheat_ice_accum_prev0(:,:) + &
                                 real(latheat_ice_accum(:,:))

    if (allocated(vels_accum))  vels_accum_prev0(:) = &
                                vels_accum_prev0(:) + &
                           real(vels_accum(:))

    if (allocated(airtemp_accum))  airtemp_accum_prev0(:) = &
                                   airtemp_accum_prev0(:) + &
                              real(airtemp_accum(:))

    if (allocated(airrrv_accum))  airrrv_accum_prev0(:) = &
                                  airrrv_accum_prev0(:) + &
                             real(airrrv_accum(:))

    if (allocated(cantemp_accum))  cantemp_accum_prev0(:) = &
                                   cantemp_accum_prev0(:) + &
                              real(cantemp_accum(:))

    if (allocated(canrrv_accum))  canrrv_accum_prev0(:) = &
                                  canrrv_accum_prev0(:) + &
                             real(canrrv_accum(:))

    if (allocated(skintemp_accum))  skintemp_accum_prev0(:) = &
                                    skintemp_accum_prev0(:) + &
                               real(skintemp_accum(:))

    if (allocated(sfluxt_accum))  sfluxt_accum_prev0(:) = &
                                  sfluxt_accum_prev0(:) + &
                             real(sfluxt_accum(:))

    if (allocated(sfluxr_accum))  sfluxr_accum_prev0(:) = & ! fast can nud
                                  sfluxr_accum_prev0(:) + &
                             real(sfluxr_accum(:))

    if (allocated(pcp_accum))  pcp_accum_prev0(:) = & ! fast can nud
                               pcp_accum_prev0(:) + &
                          real(pcp_accum(:))

    if (allocated(runoff_accum))  runoff_accum_prev0(:) = &
                                  runoff_accum_prev0(:) + &
                             real(runoff_accum(:))

    if (allocated(sfctemp_accum))  sfctemp_accum_prev0(:) = & ! fast can nud
                                   sfctemp_accum_prev0(:) + &
                              real(sfctemp_accum(:))

    if (allocated(fracliq_accum))  fracliq_accum_prev0(:) = & ! fast can nud
                                   fracliq_accum_prev0(:) + &
                              real(fracliq_accum(:))

    if (allocated(wxferi_accum))  wxferi_accum_prev0(:) = &
                                  wxferi_accum_prev0(:) + &
                             real(wxferi_accum(:))

    if (allocated(wxferp_accum))  wxferp_accum_prev0(:) = &
                                  wxferp_accum_prev0(:) + &
                             real(wxferp_accum(:))

    if (allocated(wxfer1_accum))  wxfer1_accum_prev0(:) = &
                                  wxfer1_accum_prev0(:) + &
                             real(wxfer1_accum(:))

    if (allocated(airtemp_dmin_accum))  airtemp_dmin_accum_prev0(:) = &
                                        airtemp_dmin_accum_prev0(:) + &
                                   real(airtemp_dmin_accum(:))

    if (allocated(airtemp_dmax_accum))  airtemp_dmax_accum_prev0(:) = &
                                        airtemp_dmax_accum_prev0(:) + &
                                   real(airtemp_dmax_accum(:))

    if (allocated(cantemp_dmin_accum))  cantemp_dmin_accum_prev0(:) = &
                                        cantemp_dmin_accum_prev0(:) + &
                                   real(cantemp_dmin_accum(:))

    if (allocated(cantemp_dmax_accum))  cantemp_dmax_accum_prev0(:) = &
                                        cantemp_dmax_accum_prev0(:) + &
                                   real(cantemp_dmax_accum(:))

    if (allocated(vegtemp_dmin_accum))  vegtemp_dmin_accum_prev0(:) = &
                                        vegtemp_dmin_accum_prev0(:) + &
                                   real(vegtemp_dmin_accum(:))

    if (allocated(vegtemp_dmax_accum))  vegtemp_dmax_accum_prev0(:) = &
                                        vegtemp_dmax_accum_prev0(:) + &
                                   real(vegtemp_dmax_accum(:))

    if (allocated(soiltemp_dmin_accum)) soiltemp_dmin_accum_prev0(:) = &
                                        soiltemp_dmin_accum_prev0(:) + &
                                   real(soiltemp_dmin_accum(:))

    if (allocated(soiltemp_dmax_accum)) soiltemp_dmax_accum_prev0(:) = &
                                        soiltemp_dmax_accum_prev0(:) + &
                                   real(soiltemp_dmax_accum(:))

!------------------------------------------------------------------------------
! Exception for mem_plot: Not an accumulated quantity but an instantaneous one.
! Perhaps in future add to mem_flux_accum arrays.

    do iland = 2, mland
       iwsfc = iland + omland
       if (itab_wsfc(iwsfc)%irank == myrank) then
          do k = 1,nzg
             soil_water_tot_prev0(iland) = soil_water_tot_prev0(iland) &
                                         + land%soil_water(k,iland) * dslz(k)
          enddo
       endif
    enddo

    do iwsfc = 2, mwsfc
       if (itab_wsfc(iwsfc)%irank /= myrank) cycle

       if (sfcg%leaf_class(iwsfc) < 2) then

          head_wtab_prev0(iwsfc) = head_wtab_prev0(iwsfc) &
                                 + sfcg%head1(iwsfc)
       else

          iland = iwsfc - omland

          do klev = nzg,1,-1
             call soil_wat2pot(klev, iland, land%soil_water(klev,iland), &
                  land%wresid_vg(klev,iland), land%wsat_vg(klev,iland), &
                  land%alpha_vg(klev,iland), land%en_vg(klev,iland), psi, psi_slope)

             ! Trial algorithm: Get head_wtab from highest saturated soil level

             if (psi > 1.e-2) then
                head(klev) = psi + slzt(klev)
                head_wtab_prev0(iwsfc) = head_wtab_prev0(iwsfc) &
                                       + head(klev)
                exit
             else
                head(klev) = psi + slzt(klev)
                if (klev == 1)  head_wtab_prev0(iwsfc) = head_wtab_prev0(iwsfc) &
                                                       + head(klev)
             endif
          enddo

       endif
    enddo
!------------------------------------------------------------------------------

  end subroutine copy_plot

End Module mem_plot
