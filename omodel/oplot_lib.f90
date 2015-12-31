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
subroutine oplot_lib(kk,ii,infotyp,fldname0,wtbot,wttop,fldval,notavail)

use mem_ijtabs,  only: itab_w, itab_v, itab_m, itabg_w, &
                       jtab_w, jtab_v, jtab_m, &
                       jtw_prog, jtv_prog, jtm_vadj
use mem_basic,   only: vmc, wmc, vmp, vc, wc, rho, press, &
                       thil, theta, tair, sh_w, sh_v, vxe, vye, vze
use mem_cuparm,  only: conprr, aconpr, qwcon

use mem_grid,    only: mza, mva, mwa, lpm, lpv, lpw, lsw, &
                       zm, zt, dnu, dnv, arw0, arm0, arv, arw, volt, &
                       volti, xem, yem, zem, &
                       xev, yev, zev, xew, yew, zew, &
                       topm, topw, glatm, glonm, glatv, glonv, glatw, glonw, &
                       unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz

use mem_leaf,    only: land, itab_wl

use mem_sea,     only: sea, itab_ws

use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                       con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                       con_ccn, con_gccn, con_ifn, &
                       pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                       accpd, accpr, accpp, accps, accpa, accpg, accph, &
                       cldnum
                      
use mem_radiate, only: fthrd_lw, fthrd_sw, rshort, rlong, rlongup, albedt, &
                       rshort_top, rshortup_top, rlongup_top
use mem_addsc,   only: addsc
use mem_tend,    only: vmt, wmt
use mem_para,    only: myrank
use mem_turb,    only: vkm_sfc, sfluxt, sfluxr, pblh, hkm, hkh, ustar, wstar
use mem_nudge,   only: rho_obs, theta_obs, shw_obs, uzonal_obs, umerid_obs, &
                       rho_sim, theta_sim, shw_sim, uzonal_sim, umerid_sim

use misc_coms,   only: io6, pr01d, dn01d, th01d, time8, isubdomain, &
                       naddsc, mdomain
use oplot_coms,  only: op
use consts_coms, only: p00i, rocp, erad, piu180, cp, alvl, grav, omega2
use leaf_coms,   only: slcpd, nzg, slmstsi_ch, slmstsi_vg, slz, dslz, mwl, dt_leaf
use sea_coms,    only: mws
use mem_flux_accum, only:     rshort_accum,         rshortup_accum, &
                               rlong_accum,          rlongup_accum, &
                          rshort_top_accum,     rshortup_top_accum, &
                         rlongup_top_accum,                         &
                          rshort_clr_accum,     rshortup_clr_accum, &
                           rlong_clr_accum,      rlongup_clr_accum, &
                      rshort_top_clr_accum, rshortup_top_clr_accum, &
                     rlongup_top_clr_accum,                         &
                              sfluxt_accum,           sfluxr_accum, &
                                  vc_accum,               wc_accum, &
                               press_accum,             tair_accum, &
                                sh_v_accum,                         &
                              vels_l_accum,                         &
                           airtemp_l_accum,         airshv_l_accum, &
                           cantemp_l_accum,         canshv_l_accum, &
                          skintemp_l_accum,                         &
                            sfluxt_l_accum,         sfluxr_l_accum, &
                            wxfer1_l_accum,                         &
                              vels_s_accum,                         &
                           airtemp_s_accum,         airshv_s_accum, &
                           cantemp_s_accum,         canshv_s_accum, &
                          skintemp_s_accum,                         &
                            sfluxt_s_accum,         sfluxr_s_accum

  use mem_average_vars, only: &
     npoints_mavg, npoints_davg, npoints_avg24, nz_avg, &
     press_mavg, rho_mavg, tempk_mavg, sh_v_mavg, sh_w_mavg, wc_mavg, &
     vxe_mavg, vye_mavg, vze_mavg, &
     rshort_mavg, rshort_top_mavg, rshortup_mavg, rshortup_top_mavg, &
     rlong_mavg, rlongup_mavg, rlongup_top_mavg, &
     latflux_mavg, sensflux_mavg, windspeed_mavg, &
     accpmic_mtot, accpcon_mtot, &
     press_avg24, tempk_avg24, shum_avg24, vxe_avg24, vye_avg24, vze_avg24, &
     rshort_avg24, rshort_top_avg24, rshortup_avg24, rshortup_top_avg24, &
     rlong_avg24, rlongup_avg24, rlongup_top_avg24, &
     latflux_avg24, sensflux_avg24, &
     accpmic_tot24, accpcon_tot24, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     airtempk_l_davg, airtempk_l_dmin, airtempk_l_dmax, &
     cantempk_l_davg, cantempk_l_dmin, cantempk_l_dmax, &
       vegtempk_davg,   vegtempk_dmin,   vegtempk_dmax, &
      soiltempk_davg,  soiltempk_dmin,  soiltempk_dmax, &
       sfluxt_l_davg,   sfluxr_l_davg,                  &
     airtempk_s_davg, airtempk_s_dmin, airtempk_s_dmax, &
     cantempk_s_davg, cantempk_s_dmin, cantempk_s_dmax, &
       sfluxt_s_davg,   sfluxr_s_davg

use oname_coms, only: nl

use mem_swtc5_refsoln_cubic

use mem_plot, only: &
      press_init,      rho_init,    theta_init, vc_init, addsc1_init, &
     time8_prev0,   time8_prev1,   time8_prev2,   time8_prev3, &
   accpmic_prev0, accpmic_prev1, accpmic_prev2, accpmic_prev3, &
   accpcon_prev0, accpcon_prev1, accpcon_prev2, accpcon_prev3, &
         rshort_accum_prev0,       rshort_accum_prev1, &
       rshortup_accum_prev0,     rshortup_accum_prev1, &
          rlong_accum_prev0,        rlong_accum_prev1, &
        rlongup_accum_prev0,      rlongup_accum_prev1, &
     rshort_top_accum_prev0,   rshort_top_accum_prev1, &
   rshortup_top_accum_prev0, rshortup_top_accum_prev1, &
    rlongup_top_accum_prev0,  rlongup_top_accum_prev1, &
         rshort_clr_accum_prev0,       rshort_clr_accum_prev1, &
       rshortup_clr_accum_prev0,     rshortup_clr_accum_prev1, &
          rlong_clr_accum_prev0,        rlong_clr_accum_prev1, &
        rlongup_clr_accum_prev0,      rlongup_clr_accum_prev1, &
     rshort_top_clr_accum_prev0,   rshort_top_clr_accum_prev1, &
   rshortup_top_clr_accum_prev0, rshortup_top_clr_accum_prev1, &
    rlongup_top_clr_accum_prev0,  rlongup_top_clr_accum_prev1, &
         sfluxt_accum_prev0,       sfluxt_accum_prev1, &
         sfluxr_accum_prev0,       sfluxr_accum_prev1, &
             vc_accum_prev0,           vc_accum_prev1, &
             wc_accum_prev0,           wc_accum_prev1, &
          press_accum_prev0,        press_accum_prev1, &
           tair_accum_prev0,         tair_accum_prev1, &
           sh_v_accum_prev0,         sh_v_accum_prev1, &
         vels_l_accum_prev0,       vels_l_accum_prev1, &
      airtemp_l_accum_prev0,    airtemp_l_accum_prev1,  airtemp_l_accum_prev2,  airtemp_l_accum_prev3, &
       airshv_l_accum_prev0,     airshv_l_accum_prev1, &
      cantemp_l_accum_prev0,    cantemp_l_accum_prev1,  cantemp_l_accum_prev2,  cantemp_l_accum_prev3, &
       canshv_l_accum_prev0,     canshv_l_accum_prev1, &
     skintemp_l_accum_prev0,   skintemp_l_accum_prev1, skintemp_l_accum_prev2, skintemp_l_accum_prev3, &
       sfluxt_l_accum_prev0,     sfluxt_l_accum_prev1, &
       sfluxr_l_accum_prev0,     sfluxr_l_accum_prev1, &
       wxfer1_l_accum_prev0,     wxfer1_l_accum_prev1, &
         vels_s_accum_prev0,       vels_s_accum_prev1, &
      airtemp_s_accum_prev0,    airtemp_s_accum_prev1,  airtemp_s_accum_prev2,  airtemp_s_accum_prev3, &
       airshv_s_accum_prev0,     airshv_s_accum_prev1, &
      cantemp_s_accum_prev0,    cantemp_s_accum_prev1,  cantemp_s_accum_prev2,  cantemp_s_accum_prev3, &
       canshv_s_accum_prev0,     canshv_s_accum_prev1, &
     skintemp_s_accum_prev0,   skintemp_s_accum_prev1, skintemp_s_accum_prev2, skintemp_s_accum_prev3, &
       sfluxt_s_accum_prev0,     sfluxt_s_accum_prev1, &
       sfluxr_s_accum_prev0,     sfluxr_s_accum_prev1

implicit none

integer, intent(in) :: kk,ii
character(len=*), intent(in) :: infotyp,fldname0
real, intent(out) :: fldval
real, intent(in) :: wtbot, wttop

integer, intent(out) :: notavail  ! 0 - variable is available
                                  ! 1 - variable is below ground
                                  ! 2 - variable is above model top
                                  ! 3 - variable is not available in this run
                                  ! 4 - variable is not available in current grid cell
                                  !     (e.g., no sfcwater fracliq when no sfcwater)

integer :: klev,nls,jv,iv,iw,kw,kp,k,i
real :: raxis,u,v,farv2,rpolyi
real :: vx, vy, vz, vxc, vyc, vzc
real :: tempk,fracliq
real :: vc_change, wc_change
integer :: iw1,iw2,iwl,iws
integer :: npoly,j
integer :: lenstr, ic, ifield
integer, save :: indp, icase
integer :: iland, jland, isea, jsea
real :: area_sum

real :: vcc
real :: vcc_init, vx_init, vy_init, vz_init, u_init, v_init
real :: accpboth_prev0, accpboth_prev1, accpboth_prev2, accpboth_prev3
real :: denom, fldval1, fldval2

integer :: isf, ilf
integer, save :: icall = 0
real, save, allocatable :: aux(:)

real :: zanal_swtc5, zanal0_swtc5

integer, parameter :: nfields = 503
character(len=40) :: fldlib(4,nfields)
character(len=40), save :: fldname

!  fldname     stagpt/dimens     field description & units                   field #
!-----------------------------------------------------------------------------------
! ATMOSPHERE - 3D

data fldlib(1:4,  1:34)/ &
 'VMC'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  1
 'WMC'           ,'W3' ,'W MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'            ,& !p  2
 'VMP'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  3
 'VC'            ,'V3' ,'V-NORMAL VELOCITY',' (m s:S2:-1  )'                ,& !p  4
 'WC'            ,'W3' ,'W VELOCITY',' (m s:S2:-1  )'                       ,& !p  5
 'RHO'           ,'T3' ,'AIR DENSITY',' (kg m:S2:-3  )'                     ,& !p  6
 'PRESS'         ,'T3' ,'PRESSURE',' (hPa)'                                 ,& !   7
 'THIL'          ,'T3' ,'ICE-LIQUID THETA',' (K)'                           ,& !p  8
 'THETA'         ,'T3' ,'THETA',' (K)'                                      ,& !p  9
 'AIRTEMPK'      ,'T3' ,'AIR TEMP',' (K)'                                   ,& !p 10
 'AIRTEMPC'      ,'T3' ,'AIR TEMP',' (C)'                                   ,& !p 11
 'SH_W'          ,'T3' ,'TOTAL WATER SPEC DENSITY',' (g kg:S2:-1  )'        ,& !p 12
 'SH_V'          ,'T3' ,'WATER VAPOR SPEC DENSITY',' (g kg:S2:-1  )'        ,& !p 13
 'SH_C'          ,'T3' ,'CLOUDWATER SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 14
 'SH_D'          ,'T3' ,'DRIZZLE SPEC DENSITY',' (g kg:S2:-1  )'            ,& !p 15
 'SH_R'          ,'T3' ,'RAIN SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 16
 'SH_P'          ,'T3' ,'PRISTINE ICE SPEC DENSITY',' (g kg:S2:-1  )'       ,& !p 17
 'SH_S'          ,'T3' ,'SNOW SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 18
 'SH_A'          ,'T3' ,'AGGREGATES SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 19
 'SH_G'          ,'T3' ,'GRAUPEL SPEC DENSITY',' (g kg:S2:-1  )'            ,& !p 20
 'SH_H'          ,'T3' ,'HAIL SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 21
 'SH_CP'         ,'T3' ,'CLOUD + PRIST ICE SPEC DENSITY',' (g kg:S2:-1  )'  ,& !p 22
 'SH_TOTCOND'    ,'T3' ,'CONDENSATE SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 23
 'CON_C'         ,'T3' ,'CLOUD DROPLET NUMBER CONCEN',' (# mg:S2:-1  )'     ,& !p 24
 'CON_D'         ,'T3' ,'DRIZZLE NUMBER CONCEN',' (# mg:S2:-1  )'           ,& !p 25
 'CON_R'         ,'T3' ,'RAIN NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 26
 'CON_P'         ,'T3' ,'PRISTINE ICE NUMBER CONCEN',' (# kg:S2:-1  )'      ,& !p 27
 'CON_S'         ,'T3' ,'SNOW NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 28
 'CON_A'         ,'T3' ,'AGGREGATES NUMBER CONCEN',' (# kg:S2:-1  )'        ,& !p 29
 'CON_G'         ,'T3' ,'GRAUPEL NUMBER CONCEN',' (# kg:S2:-1  )'           ,& !p 30
 'CON_H'         ,'T3' ,'HAIL NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 31
 'CON_CCN'       ,'T3' ,'CCN NUMBER CONCEN',' (# mg:S2:-1  )'               ,& !p 32
 'CON_GCCN'      ,'T3' ,'GCCN NUMBER CONCEN',' (# mg:S2:-1  )'              ,& !p 33
 'CON_IFN'       ,'T3' ,'IFN NUMBER CONCEN',' (# kg:S2:-1  )'                / !p 34

data fldlib(1:4, 35:56)/ &
 'HKM'           ,'T3' ,'HORIZ TURB MOMENTUM K',' (N s m:S2:-2  )'          ,& !p 35
 'FTHRD'         ,'T3' ,'RADIATIVE THETA TENDENCY',' (K s:S2:-1  )'         ,& !p 36
 'SPEEDW'        ,'T3' ,'WIND SPEED AT W',' (m s:S2:-1  )'                  ,& !p 37
 'AZIMW'         ,'T3' ,'WIND AZIMUTH AT W',' (deg)'                        ,& !p 38
 'ZONAL_WINDW'   ,'T3' ,'ZONAL WIND AT W',' (m s:S2:-1  )'                  ,& !p 39
 'MERID_WINDW'   ,'T3' ,'MERIDIONAL WIND AT W',' (m s:S2:-1  )'             ,& !p 40
 'RVORTZM'       ,'P3' ,'REL VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 41
 'TVORTZM'       ,'P3' ,'TOT VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 42
 'RVORTZM_P'     ,'P3' ,'REL VERT VORTICITY PERT AT M',' (s:S2:-1  )'       ,& !p 43
 'DIVERG'        ,'T3' ,'HORIZONTAL DIVERGENCE',' (s:S2:-1  )'              ,& !p 44
 'VMASSFLUX'     ,'V3' ,'GRID CELL V-FACE MASS FLUX',' (kg s:S2:-1  )'      ,& !  45
 'VC_P'          ,'V3' ,'NORMAL WIND PERT AT V',' (m s:S2:-1  )'            ,& !p 46
 'PRESS_P'       ,'T3' ,'PRESSURE PERT',' (hPa)'                            ,& !  47
 'RHO_P'         ,'T3' ,'DENSITY PERT',' (kg m:S2:-3  )'                    ,& !  48
 'THETA_P'       ,'T3' ,'THETA PERT',' (K)'                                 ,& !  49
 'AIRTEMPK_P'    ,'T3' ,'AIR TEMP PERT',' (K)'                              ,& !p 50
 'VMT'           ,'V3' ,'V-NORM MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'  ,& !  51
 'WMT'           ,'W3' ,'W MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'       ,& !  52
 'ADDSC'         ,'T3' ,'ADDED SCALAR AMOUNT PER KG AIR',' '                ,& !p 53
 'ADDSCP'        ,'T3' ,'SCALAR PERTURBATION',' ( )'                        ,& !  54
 'ZPLEV'         ,'T3' ,'HEIGHT OF CONST P SFC',' (m)'                      ,& !p 55
 'QWCON'         ,'T3' ,'CUPARM CONDENSATE SPEC DENSITY',' (g kg:S2:-1  )'   / !p 56

! ATMOSPHERE - 2D

data fldlib(1:4, 62:95)/ &
 'RSHORT_TOP'    ,'T2' ,'TOP DOWN SHORTWV FLX',' (W m:S2:-2  )'             ,& !  62
 'RSHORTUP_TOP'  ,'T2' ,'TOP UP SHORTWV FLX',' (W m:S2:-2  )'               ,& !  63
 'RLONGUP_TOP'   ,'T2' ,'TOP UP LONGWV FLX',' (W m:S2:-2  )'                ,& !  64

! ATMOSPHERE SURFACE - 2D

 'RSHORT'        ,'T2' ,'SFC DOWN SHORTWV FLX',' (W m:S2:-2  )'             ,& !  65
 'RSHORTUP'      ,'T2' ,'SFC UP SHORTWV FLX',' (W m:S2:-2  )'               ,& !  66
 'RLONG'         ,'T2' ,'SFC DOWN LONGWV FLX',' (W m:S2:-2  )'              ,& !  67
 'RLONGUP'       ,'T2' ,'SFC UP LONGWV FLX',' (W m:S2:-2  )'                ,& !  68
 'ALBEDT'        ,'T2' ,'NET GRID COLUMN SFC ALBEDO',' ( )'                 ,& !  69
 'VKM_SFC'       ,'T2' ,'SFC TURB K FOR MOMENTUM',' (N s m:S2:-2  )'        ,& !  70
 'USTAR'         ,'T2' ,'SFC FRICTION VELOCITY',' (m s:S2:-1  )'            ,& !  71
 'SENSFLUX'      ,'T2' ,'ATM SFC SENSIBLE HEAT FLUX',' (W m:S2:-2  )'       ,& !  72
 'VAPFLUX'       ,'T2' ,'ATM SFC VAPOR FLUX',' (kg m:S2:-2   s:S2:-1  )'    ,& !  73
 'LATFLUX'       ,'T2' ,'ATM SFC LATENT HEAT FLUX',' (W m:S2:-2  )'         ,& !  74
 'PCPRD'         ,'T2' ,'DRIZZLE PRECIP RATE',' (kg m:S2:-2   h:S2:-1  )'   ,& !  75
 'PCPRR'         ,'T2' ,'RAIN PRECIP RATE',' (kg m:S2:-2   h:S2:-1  )'      ,& !  76
 'PCPRP'         ,'T2' ,'PRIST ICE PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    ,& !  77
 'PCPRS'         ,'T2' ,'SNOW PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  78
 'PCPRA'         ,'T2' ,'AGGREGATES PCP RATE',' (kg m:S2:-2   h:S2:-1  )'   ,& !  79
 'PCPRG'         ,'T2' ,'GRAUPEL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'      ,& !  80
 'PCPRH'         ,'T2' ,'HAIL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  81
 'PCPRMIC'       ,'T2' ,'MICROPHYS PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    ,& !  82
 'PCPRCON'       ,'T2' ,'CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  83
 'PCPRBOTH'      ,'T2' ,'MICRO + CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )' ,& !  84
 'ACCPD'         ,'T2' ,'ACCUM DRIZZLE',' (kg m:S2:-2  )'                   ,& !  85
 'ACCPR'         ,'T2' ,'ACCUM RAIN',' (kg m:S2:-2  )'                      ,& !  86
 'ACCPP'         ,'T2' ,'ACCUM PRIST ICE',' (kg m:S2:-2  )'                 ,& !  87
 'ACCPS'         ,'T2' ,'ACCUM SNOW',' (kg m:S2:-2  )'                      ,& !  88
 'ACCPA'         ,'T2' ,'ACCUM AGGREGATES',' (kg m:S2:-2  )'                ,& !  89
 'ACCPG'         ,'T2' ,'ACCUM GRAUPEL',' (kg m:S2:-2  )'                   ,& !  90
 'ACCPH'         ,'T2' ,'ACCUM HAIL',' (kg m:S2:-2  )'                      ,& !  91
 'ACCPMIC'       ,'T2' ,'ACCUM MICPHYS PCP',' (kg m:S2:-2  )'               ,& !  92
 'ACCPCON'       ,'T2' ,'ACCUM CONV PCP',' (kg m:S2:-2  )'                  ,& !  93
 'ACCPBOTH'      ,'T2' ,'ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )'        ,& !  94
 'WSTAR'         ,'T2' ,'PBL CONVECTIVE VELOCITY',' (m s:S2:-1  )'           / !  95

! ATMOSPHERE DIF2 fields (3D & 2D)

data fldlib(1:4,101:132)/ &
 'ZONAL_WINDW_DIF2' ,'T3' ,'ATM ZONAL VELOCITY DIF2',' (m s:S2:-1  )'       ,& ! 101
 'MERID_WINDW_DIF2' ,'T3' ,'ATM MERID VELOCITY DIF2',' (m s:S2:-1  )'       ,& ! 102
 'WC_DIF2'          ,'T3' ,'ATM VERT VELOCITY DIF2',' (m s:S2:-1  )'        ,& ! 103
 'PRESS_DIF2'       ,'T3' ,'ATM PRESSURE DIF2',' (hPa)'                     ,& ! 104
 'AIRTEMPK_DIF2'    ,'T3' ,'ATM TEMP DIF2',' (K)'                           ,& ! 105
 'SH_V_DIF2'        ,'T3' ,'ATM VAP SPECIFIC DENSITY DIF2',' (g kg:S2:-1  )',& ! 106
 'PCPMIC_DIF2'      ,'T2' ,'MICPHYS PRECIP DIFF2',' (mm/day)'               ,& ! 107
 'PCPCON_DIF2'      ,'T2' ,'CONV PRECIP DIFF2',' (mm/day)'                  ,& ! 108
 'PCPBOTH_DIF2'     ,'T2' ,'MICPHYS + CONV PRECIP DIFF2',' (mm/day)'        ,& ! 109
 'RSHORT_DIF2'      ,'T2' ,'SFC DOWN SHORTWV FLX DIF2',' (W m:S2:-2  )'     ,& ! 110
 'RSHORTUP_DIF2'    ,'T2' ,'SFC UP SHORTWV FLX DIF2',' (W m:S2:-2  )'       ,& ! 111
 'RLONG_DIF2'       ,'T2' ,'SFC DOWN LONGWV FLX DIF2',' (W m:S2:-2  )'      ,& ! 112
 'RLONGUP_DIF2'     ,'T2' ,'SFC UP LONGWV FLX DIF2',' (W m:S2:-2  )'        ,& ! 113
 'RSHORT_TOP_DIF2'  ,'T2' ,'TOP DOWN SHORTWV FLX DIF2',' (W m:S2:-2  )'     ,& ! 114
 'RSHORTUP_TOP_DIF2','T2' ,'TOP UP SHORTWV FLX DIF2',' (W m:S2:-2  )'       ,& ! 115
 'RLONGUP_TOP_DIF2' ,'T2' ,'TOP UP LONGWV FLX DIF2',' (W m:S2:-2  )'        ,& ! 116
 'SENSFLUX_DIF2'    ,'T2' ,'ATM SFC SENS HEAT FLUX DIF2',' (W m:S2:-2  )'   ,& ! 117
 'LATFLUX_DIF2'     ,'T2' ,'ATM SFC LAT HEAT FLUX DIF2',' (W m:S2:-2  )'    ,& ! 118
 'VAPFLUX_DIF2'     ,'T2' ,'ATM SFC VAP FLUX DIF2',' (kg m:S2:-2   s:S2:-1  )' ,& ! 119
 'RSHORT_CLR_DIF2'      ,'T2' ,'SFC DOWN SHORTWV FLX CLR DIF2',' (W m:S2:-2  )',& ! 120
 'RSHORTUP_CLR_DIF2'    ,'T2' ,'SFC UP SHORTWV FLX CLR DIF2',' (W m:S2:-2  )'  ,& ! 121
 'RLONG_CLR_DIF2'       ,'T2' ,'SFC DOWN LONGWV FLX CLR DIF2',' (W m:S2:-2  )' ,& ! 122
 'RLONGUP_CLR_DIF2'     ,'T2' ,'SFC UP LONGWV FLX CLR DIF2',' (W m:S2:-2  )'   ,& ! 123
 'RSHORT_TOP_CLR_DIF2'  ,'T2' ,'TOP DOWN SHORTWV FLX CLR DIF2',' (W m:S2:-2  )',& ! 124
 'RSHORTUP_TOP_CLR_DIF2','T2' ,'TOP UP SHORTWV FLX CLR DIF2',' (W m:S2:-2  )'  ,& ! 125
 'RLONGUP_TOP_CLR_DIF2' ,'T2' ,'TOP UP LONGWV FLX CLR DIF2',' (W m:S2:-2  )'   ,& ! 126

! ATMOSPHERE DIF4 fields (2D)

 'PCPMIC_DIF4'    ,'T2' ,'MICPHYS PRECIP DIFF4',' (mm/day)'                  ,& ! 127
 'PCPCON_DIF4'    ,'T2' ,'CONV PRECIP DIFF4',' (mm/day)'                     ,& ! 128
 'PCPBOTH_DIF4'   ,'T2' ,'MICPHYS + CONV PRECIP DIFF4',' (mm/day)'           ,& ! 129
 'PCPMIC_REL4'    ,'T2' ,'MICPHYS PRECIP RELATIVE DIFF4',' '                 ,& ! 130
 'PCPCON_REL4'    ,'T2' ,'CONV PRECIP RELATIVE DIFF4',' '                    ,& ! 131
 'PCPBOTH_REL4'   ,'T2' ,'MICPHYS + CONV PRECIP RELATIVE DIFF4',' '           / ! 132

! LAND_CELLS - 3D

data fldlib(1:4,151:183)/ &
 'SOIL_TEXT'     ,'L3G','SOIL TEXTURAL CLASS',' ( )'                        ,& ! 151
 'SOIL_ENERGY'   ,'L3G','SOIL ENERGY',' (J cm:S2:-3  )'                     ,& ! 152
 'SOIL_TEMPK'    ,'L3G','SOIL TEMP',' (K)'                                  ,& ! 153
 'SOIL_FRACLIQ'  ,'L3G','LIQUID FRACTION OF SOIL WATER',' ( )'              ,& ! 154
 'SOIL_WATER'    ,'L3G','SOIL WATER CONTENT',' ( )'                         ,& ! 155
 'SFWAT_MASS'    ,'L3S','SFCWATER MASS',' (kg m:S2:-2  )'                   ,& ! 156
 'SFWAT_ENERGY'  ,'L3S','SFCWATER ENERGY',' (J g:S2:-1  )'                  ,& ! 157
 'SFWAT_TEMPK'   ,'L3S','SFCWATER TEMP',' (K)'                              ,& ! 158
 'SFWAT_FRACLIQ' ,'L3S','SFCWATER LIQUID FRACTION',' ( )'                   ,& ! 159
 'SFWAT_DEPTH'   ,'L3S','SFCWATER DEPTH',' (m)'                             ,& ! 160

! LAND_CELLS - 2D

 'NLEV_SFWAT'       ,'L2' ,'NUMBER OF SFCWATER LAYERS',' ( )'                ,& ! 161
 'VEG_NDVIC'        ,'L2' ,'VEGETATION NDVI',' ( )'                          ,& ! 162
 'VEG_TEMPC'        ,'L2' ,'VEGETATION TEMP',' (C)'                          ,& ! 163
 'VEG_TEMPK'        ,'L2' ,'VEGETATION TEMP',' (K)'                          ,& ! 164
 'VEG_WATER'        ,'L2' ,'VEGETATION SFC WATER ',' (kg m:S2:-2  )'         ,& ! 165
 'STOM_RESIST'      ,'L2' ,'STOMATAL RESISTANCE',' (s m:S2:-1  )'            ,& ! 166
 'SFCWATER_TOT'     ,'L2' ,'TOTAL SFCWATER MASS',' (kg m:S2:-2  )'           ,& ! 167
 'SFCWATER_TOP_TEMP','L2' ,'SFCWATER TOPLAYER TEMP',' (K)'                   ,& ! 168
 'SOIL_TOP_TEMP'    ,'L2' ,'SOIL TOPLAYER TEMP',' (K)'                       ,& ! 169
 'GROUND_SHV'       ,'L2' ,'EQUIL VAP SPEC DENSITY OF SOIL',' (g kg:S2:-1  )',& ! 170
 'SOIL_DEPTH'       ,'L2' ,'SOIL DEPTH',' (m)'                               ,& ! 171
 'SOIL_WATER_TOT'   ,'L2' ,'TOTAL SOIL WATER',' (m)'                         ,& ! 172
 'HEAD0'            ,'L2' ,'HEAD0',' (m)'                                    ,& ! 173

! LAND_CELLS - DIF2 fields

 'L_WXFER1_DIF2'    ,'L2' ,'L SOIL BOTTOM WATER FLUX DIF2',' (mm/day)'       ,& ! 174

! LAND_CELLS - ATM averages

 'AL_SFCWATER_TOT'  ,'T2' ,'AL TOTAL SFCWATER MASS',' (kg m:S2:-2  )'        ,& ! 175
 'AL_SOIL_WATER_TOT','T2' ,'AL TOTAL SOIL WATER',' (m)'                      ,& ! 176

! LAND_CELLS - ATM averages of DIF2 fields

 'AL_WXFER1_DIF2'   ,'T2' ,'AL SOIL BOTTOM WATER FLUX DIF2',' (mm/day)'      ,& ! 177

! SEA_CELLS - 2D

 'SEATP'         ,'S2' ,'SEA SFC TEMP (PAST DATA)',' (K)'                    ,& ! 178
 'SEATF'         ,'S2' ,'SEA SFC TEMP (FUTURE DATA)',' (K)'                  ,& ! 179
 'SEATC'         ,'S2' ,'SEA SFC TEMP (CURRENT)',' (K)'                      ,& ! 180
 'SEAICEP'       ,'S2' ,'SEAICE FRACTION (PAST DATA)',' ( )'                 ,& ! 181
 'SEAICEF'       ,'S2' ,'SEAICE FRACTION (FUTURE DATA)',' ( )'               ,& ! 182
 'SEAICEC'       ,'S2' ,'SEAICE FRACTION (CURRENT)',' ( )'                    / ! 183

! LAND AND SEA CELLS - 2D

data fldlib(1:4,201:224)/ &
 'LEAF_CLASS'       ,'B2' ,'LEAF CLASS',' ( )'                             ,& ! 201
 'LS_IW'            ,'B2' ,'LS ATM IW INDEX',' '                           ,& ! 202
 'LS_KW'            ,'B2' ,'LS ATM KW INDEX',' '                           ,& ! 203
 'LS_ARF_IW'        ,'B2' ,'LS IW AREA FRACTION',' '                       ,& ! 204
 'LS_ARF_KW'        ,'B2' ,'LS KW AREA FRACTION',' '                       ,& ! 205
 'LS_AREA'          ,'B2' ,'LS CELL AREA',' (m:S2:2  )'                    ,& ! 206
 'LS_TOPW'          ,'B2' ,'LS CELL TOPW',' (m)'                           ,& ! 207

 'LS_ROUGH'         ,'B2' ,'LS NET ROUGHNESS HEIGHT',' (m)'                ,& ! 208
 'LS_VELS'          ,'B2' ,'LS WIND SPEED',' (m s:S2:-1  )'                ,& ! 209
 'LS_AIRTEMPK'      ,'B2' ,'LS ATM TEMP',' (K)'                            ,& ! 210
 'LS_AIRSHV'        ,'B2' ,'LS ATM SHV',' (g kg:S2:-1  )'                  ,& ! 211
 'LS_CANTEMPK'      ,'B2' ,'LS CANOPY AIR TEMP',' (K)'                     ,& ! 212
 'LS_CANSHV'        ,'B2' ,'LS CANOPY VAPOR SPEC DENSITY',' (g kg:S2:-1  )',& ! 213
 'LS_SKINTEMPK'     ,'B2' ,'LS VEG/GROUND/SFCWATER/SEA TEMP',' (K)'        ,& ! 214
 'LS_GSS_SSH'       ,'B2' ,'LS SFC SAT VAP SPEC DENS',' (g kg:S2:-1  )'    ,& ! 215

 'LS_SENSFLUX'      ,'B2','LS SENS HEAT FLX',' (W m:S2:-2  )'              ,& ! 216
 'LS_LATFLUX'       ,'B2','LS LAT HEAT FLX',' (W m:S2:-2  )'               ,& ! 217
 'LS_VAPFLUX'       ,'B2','LS VAP FLX',' (kg m:S2:-2   s:S2:-1  )'         ,& ! 218

 'LS_RSHORT'        ,'B2','LS DOWN SW FLX',' (W m:S2:-2  )'                ,& ! 219
 'LS_RLONG'         ,'B2','LS DOWN LW FLX',' (W m:S2:-2  )'                ,& ! 220
 'LS_RLONGUP'       ,'B2','LS UP LW FLX',' (W m:S2:-2  )'                  ,& ! 221
 'LS_RLONG_ALBEDO'  ,'B2','LS NET SFC LW ALBEDO',' ( )'                    ,& ! 222
 'LS_ALBEDO_BEAM'   ,'B2','LS NET SFC BEAM ALBEDO',' ( )'                  ,& ! 223
 'LS_ALBEDO_DIFFUSE','B2','LS NET SFC DIFFUSE ALBEDO',' ( )'                / ! 224

! LAND AND SEA CELLS - DIF2 fields

data fldlib(1:4,241:249)/ &
 'LS_VELS_DIF2'      ,'B2' ,'LS WIND SPEED DIF2',' (m s:S2:-1  )'           ,& ! 241
 'LS_AIRTEMPK_DIF2'  ,'B2' ,'LS ATM TEMP DIF2',' (K)'                       ,& ! 242
 'LS_AIRSHV_DIF2'    ,'B2' ,'LS ATM SHV DIF2',' (g kg:S2:-1  )'             ,& ! 243
 'LS_CANTEMPK_DIF2'  ,'B2' ,'LS CANOPY AIR TEMP DIF2',' (K)'                ,& ! 244
 'LS_CANSHV_DIF2'    ,'B2' ,'LS CANOPY VAP SPEC DENS DIF2',' (g kg:S2:-1  )',& ! 245
 'LS_SKINTEMPK_DIF2' ,'B2' ,'LS VEG/GROUND/SFCWATER/SEA TEMP DIF2',' (K)'   ,& ! 246
 'LS_SENSFLUX_DIF2'  ,'B2' ,'LS SENS HEAT FLUX DIF2',' (W m:S2:-2  )'       ,& ! 247
 'LS_LATFLUX_DIF2'   ,'B2' ,'LS LAT HEAT FLUX DIF2',' (W m:S2:-2  )'        ,& ! 248
 'LS_VAPFLUX_DIF2'   ,'B2' ,'LS VAP FLUX DIF2',' (kg m:S2:-2   s:S2:-1  )'   / ! 249

! LAND AND SEA CELLS - ATM averages

data fldlib(1:4,261:269)/ &
 'ALS_VELS'         ,'T2' ,'ALS WIND SPEED',' (m s:S2:-1  )'                ,& ! 261
 'ALS_AIRTEMPK'     ,'T2' ,'ALS ATM TEMP',' (K)'                            ,& ! 262
 'ALS_AIRSHV'       ,'T2' ,'ALS ATM SHV',' (g kg:S2:-1  )'                  ,& ! 263
 'ALS_CANTEMPK'     ,'T2' ,'ALS CANOPY AIR TEMP',' (K)'                     ,& ! 264
 'ALS_CANSHV'       ,'T2' ,'ALS CANOPY VAP SPEC DENS',' (g kg:S2:-1  )'     ,& ! 265
 'ALS_SKINTEMPK'    ,'T2' ,'ALS VEG/GROUND/SFCWATER/SEA TEMP',' (K)'        ,& ! 266
 'ALS_SENSFLUX'     ,'T2' ,'ALS SENS HEAT FLUX',' (W m:S2:-2  )'            ,& ! 267
 'ALS_LATFLUX'      ,'T2' ,'ALS LAT HEAT FLUX',' (W m:S2:-2  )'             ,& ! 268
 'ALS_VAPFLUX'      ,'T2' ,'ALS VAP FLUX',' (kg m:S2:-2   s:S2:-1  )'        / ! 269

! LAND AND SEA CELLS - ATM averages of DIF2 fields

data fldlib(1:4,281:289)/ &
 'ALS_VELS_DIF2'     ,'T2' ,'ALS WIND SPEED DIF2',' (m s:S2:-1  )'           ,& ! 281
 'ALS_AIRTEMPK_DIF2' ,'T2' ,'ALS ATM TEMP DIF2',' (K)'                       ,& ! 282
 'ALS_AIRSHV_DIF2'   ,'T2' ,'ALS ATM SHV DIF2',' (g kg:S2:-1  )'             ,& ! 283
 'ALS_CANTEMPK_DIF2' ,'T2' ,'ALS CANOPY AIR TEMP DIF2',' (K)'                ,& ! 284
 'ALS_CANSHV_DIF2'   ,'T2' ,'ALS CANOPY VAP SPEC DENS DIF2',' (g kg:S2:-1  )',& ! 285
 'ALS_SKINTEMPK_DIF2','T2' ,'ALS VEG/GROUND/SFCWATER/SEA TEMP DIF2',' (K)'   ,& ! 286
 'ALS_SENSFLUX_DIF2' ,'T2' ,'ALS SENS HEAT FLUX DIF2',' (W m:S2:-2  )'       ,& ! 287
 'ALS_LATFLUX_DIF2'  ,'T2' ,'ALS LAT HEAT FLUX DIF2',' (W m:S2:-2  )'        ,& ! 288
 'ALS_VAPFLUX_DIF2'  ,'T2' ,'ALS VAP FLUX DIF2',' (kg m:S2:-2   s:S2:-1  )'   / ! 289

! LAND AND SEA CELLS - ATM averages of DIF4 fields

data fldlib(1:4,290:292)/ &
 'ALS_AIRTEMPK_DIF4' ,'T2' ,'ALS ATM TEMP DIF4',' (K)'                       ,& ! 290
 'ALS_CANTEMPK_DIF4' ,'T2' ,'ALS CANOPY AIR TEMP DIF4',' (K)'                ,& ! 291
 'ALS_SKINTEMPK_DIF4','T2' ,'ALS VEG/GROUND/SFCWATER/SEA TEMP DIF4',' (K)'    / ! 292

! GRID GEOMETRY - 3D

data fldlib(1:4,301:303)/ &
 'ARV'           ,'V3' ,'AREA OF GRID CELL V-FACE',' (m:S2:2  )'            ,& ! 301
 'ARW'           ,'W3' ,'AREA OF GRID CELL W-FACE',' (m:S2:2  )'            ,& ! 302
 'VOLT'          ,'T3' ,'GRID T-CELL VOLUME',' (m:S2:3  )'                   / ! 303

! GRID GEOMETRY - 2D

data fldlib(1:4,311:336)/ &
 'TOPM'          ,'M2' ,'TOPOGRAPHY HEIGHT',' (m)'                          ,& ! 311
 'TOPW'          ,'W2' ,'TOPOGRAPHY HEIGHT AT W',' (m)'                     ,& ! 312
 'GLATM'         ,'M2' ,'LATITUDE AT M',' (deg)'                            ,& ! 313
 'GLONM'         ,'M2' ,'LONGITUDE AT M',' (deg)'                           ,& ! 314
 'GLATV'         ,'V2' ,'LATITUDE AT V',' (deg)'                            ,& ! 315
 'GLONV'         ,'V2' ,'LONGITUDE AT V',' (deg)'                           ,& ! 316
 'GLATW'         ,'T2' ,'LATITUDE',' (deg)'                                 ,& ! 317
 'GLONW'         ,'T2' ,'LONGITUDE',' (deg)'                                ,& ! 318
 'LPM'           ,'M2' ,'LOWEST PREDICTED M LEVEL',' ( )'                   ,& ! 319
 'LPV'           ,'V2' ,'LOWEST PREDICTED V LEVEL',' ( )'                   ,& ! 320
 'LCV'           ,'V2' ,'LOWEST ACTIVE V CONTROL VOL',' ( )'                ,& ! 321
 'LPW'           ,'W2' ,'LOWEST PREDICTED W LEVEL',' ( )'                   ,& ! 322
 'LSW'           ,'W2' ,'NUMBER OF SFC W LEVELS',' ( )'                     ,& ! 323
 'XEM'           ,'M2' ,'EARTH-X COORD OF M POINT',' ( )'                   ,& ! 324
 'YEM'           ,'M2' ,'EARTH-Y COORD OF M POINT',' ( )'                   ,& ! 325
 'ZEM'           ,'M2' ,'EARTH-Z COORD OF M POINT',' ( )'                   ,& ! 326
 'XEV'           ,'V2' ,'EARTH-X COORD OF V POINT',' ( )'                   ,& ! 327
 'YEV'           ,'V2' ,'EARTH-Y COORD OF V POINT',' ( )'                   ,& ! 328
 'ZEV'           ,'V2' ,'EARTH-Z COORD OF V POINT',' ( )'                   ,& ! 329
 'XEW'           ,'W2' ,'EARTH-X COORD OF W POINT',' ( )'                   ,& ! 330
 'YEW'           ,'W2' ,'EARTH-Y COORD OF W POINT',' ( )'                   ,& ! 331
 'ZEW'           ,'W2' ,'EARTH-Z COORD OF W POINT',' ( )'                   ,& ! 332
 'DNU'           ,'V2' ,'DNU',' (m)'                                        ,& ! 333
 'DNV'           ,'V2' ,'DNV',' (m)'                                        ,& ! 334
 'ARM0'          ,'W2' ,'SFC AREA OF M CELL',' (m:S2:2  )'                  ,& ! 335
 'ARW0'          ,'W2' ,'SFC AREA OF W CELL',' (m:S2:2  )'                   / ! 336

! ITAB_M MEMBERS - 2D

data fldlib(1:4,341:348)/ &
 'ITAB_M_NPOLY'  ,'M2' ,'ITAB_M_NPOLY',' ( )'                               ,& ! 341
 'ITAB_M_IMGLOBE','M2' ,'ITAB_M_IMGLOBE',' ( )'                             ,& ! 342
 'ITAB_M_MRLM'   ,'M2' ,'ITAB_M_MRLM',' ( )'                                ,& ! 343
 'ITAB_M_MRLM_OR','M2' ,'ITAB_M_MRLM_ORIG',' ( )'                           ,& ! 344
 'ITAB_M_MROW'   ,'M2' ,'ITAB_M_MROW',' ( )'                                ,& ! 345
 'ITAB_M_NGR'    ,'M2' ,'ITAB_M_NGR',' ( )'                                 ,& ! 346
 'ITAB_M_IV'     ,'M2' ,'ITAB_M_IV',' ( )'                                  ,& ! 347
 'ITAB_M_IW'     ,'M2' ,'ITAB_M_IW',' ( )'                                   / ! 348

! ITAB_V MEMBERS - 2D

data fldlib(1:4,351:358)/  &
 'ITAB_V_IVP'    ,'V2' ,'ITAB_V_IVP',' ( )'                                 ,& ! 351
 'ITAB_V_IRANK'  ,'V2' ,'ITAB_V_IRANK',' ( )'                               ,& ! 352
 'ITAB_V_IVGLOBE','V2' ,'ITAB_V_IVGLOBE',' ( )'                             ,& ! 353
 'ITAB_V_MRLV'   ,'V2' ,'ITAB_V_MRLV',' ( )'                                ,& ! 354
 'ITAB_V_IM'     ,'V2' ,'ITAB_V_IM',' ( )'                                  ,& ! 355
 'ITAB_V_IV'     ,'V2' ,'ITAB_V_IV',' ( )'                                  ,& ! 356
 'ITAB_V_IW'     ,'V2' ,'ITAB_V_IW',' ( )'                                  ,& ! 357
 'ITAB_V_FARW'   ,'V2' ,'ITAB_V_FARW',' ( )'                                 / ! 358

! ITAB_W MEMBERS - 2D

data fldlib(1:4,361:375)/ &

 'ITAB_W_NPOLY'  ,'W2' ,'ITAB_W_NPOLY',' ( )'                               ,& ! 361
 'ITAB_W_IWP'    ,'W2' ,'ITAB_W_IWP',' ( )'                                 ,& ! 362
 'ITAB_W_IRANK'  ,'W2' ,'ITAB_W_IRANK',' ( )'                               ,& ! 363
 'ITAB_W_IWGLOBE','W2' ,'ITAB_W_IWGLOBE',' ( )'                             ,& ! 364
 'ITAB_W_MRLW'   ,'W2' ,'ITAB_W_MRLW',' ( )'                                ,& ! 365
 'ITAB_W_MRLW_OR','W2' ,'ITAB_W_MRLW_ORIG',' ( )'                           ,& ! 366
 'ITAB_W_NGR'    ,'W2' ,'ITAB_W_NGR',' ( )'                                 ,& ! 367
 'ITAB_W_IM'     ,'W2' ,'ITAB_W_IM',' ( )'                                  ,& ! 368
 'ITAB_W_IV'     ,'W2' ,'ITAB_W_IV',' ( )'                                  ,& ! 369
 'ITAB_W_IW'     ,'W2' ,'ITAB_W_IW',' ( )'                                  ,& ! 370
 'ITAB_W_DIRV'   ,'W2' ,'ITAB_W_DIRV',' ( )'                                ,& ! 371
 'ITAB_W_FARM'   ,'W2' ,'ITAB_W_FARM',' ( )'                                ,& ! 372
 'ITAB_W_FARV'   ,'W2' ,'ITAB_W_FARV',' ( )'                                ,& ! 373
 'ITAB_W_IWNUD'  ,'W2' ,'ITAB_W_IWNUD',' ( )'                               ,& ! 374
 'ITAB_W_FNUD'   ,'W2' ,'ITAB_W_FNUD',' ( )'                                 / ! 375

! Monthly and Daily averaged fields; daily min & max fields - 2D

data fldlib(1:4,401:454)/ &

 'PRESS_MAVG'    ,'T2','MONTH-AVG PRESSURE',' (hPa)'                        ,& ! 401
 'RHO_MAVG'      ,'T2','MONTH-AVG AIR DENSITY',' (kg m:S2:-3  )'            ,& ! 402
 'TEMPK_MAVG'    ,'T2','MONTH-AVG TEMPERATURE',' (K)'                       ,& ! 403
 'SH_V_MAVG'     ,'T2','MONTH-AVG VAPOR SPEC DENSITY',' (g kg:S2:-1  )'     ,& ! 404
 'SH_W_MAVG'     ,'T2','MONTH-AVG TOT WATER SPEC DENSITY',' (g kg:S2:-1  )' ,& ! 405
 'WC_MAVG'       ,'T2','MONTH-AVG W VELOCITY',' (m s:S2:-1  )'              ,& ! 406
 'ZONAL_WINDW_MAVG' ,'T2','MONTH-AVG ZONAL WIND',' (m s:S2:-1  )'           ,& ! 407
 'MERID_WINDW_MAVG' ,'T2','MONTH-AVG MERID WIND',' (m s:S2:-1  )'           ,& ! 408
 'RSHORT_MAVG'      ,'T2','MONTH-AVG SFC DOWNWARD S/W FLX',' (W m:S2:-2  )' ,& ! 409
 'RSHORT_TOP_MAVG'  ,'T2','MONTH-AVG TOP DOWNWARD S/W FLX',' (W m:S2:-2  )' ,& ! 410
 'RSHORTUP_MAVG'    ,'T2','MONTH-AVG SFC UPWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 411
 'RSHORTUP_TOP_MAVG','T2','MONTH-AVG TOP UPWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 412
 'RLONG_MAVG'       ,'T2','MONTH-AVG SFC DOWNWARD L/W FLX',' (W m:S2:-2  )' ,& ! 413
 'RLONGUP_MAVG'     ,'T2','MONTH-AVG SFC UPWARD L/W FLX',' (W m:S2:-2  )'   ,& ! 414
 'RLONGUP_TOP_MAVG' ,'T2','MONTH-AVG TOP UPWARD L/W FLX',' (W m:S2:-2  )'   ,& ! 415
 'LATFLUX_MAVG'     ,'T2','MONTH-AVG SFC LATENT HEAT FLX',' (W m:S2:-2  )'  ,& ! 416
 'SENSFLUX_MAVG'    ,'T2','MONTH-AVG SFC SENSIBLE HEAT FLX',' (W m:S2:-2  )',& ! 417
 'WINDSPEED_MAVG'   ,'T2','MONTH-AVG SFC WIND SPEED',' (m s:S2:-1  )'       ,& ! 418
 'ACCPMIC_MTOT'     ,'T2','MONTH-ACCUM MICPHYS PRECIP',' (kg m:S2:-2  )'    ,& ! 419
 'ACCPCON_MTOT'     ,'T2','MONTH-ACCUM CUPARM PRECIP',' (kg m:S2:-2  )'     ,& ! 420
 'ACCPBOTH_MTOT'    ,'T2','MONTH-ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )',& ! 421
 'PRESS_DAVG'       ,'T2','DAY-AVG SFC PRESSURE',' (hPa)'                   ,& ! 422
 'ZONAL_WINDW_DAVG' ,'T2','DAY-AVG SFC ZONAL WIND',' (m s:S2:-1  )'         ,& ! 423
 'MERID_WINDW_DAVG' ,'T2','DAY-AVG SFC MERID WIND',' (m s:S2:-1  )'         ,& ! 424
 'RSHORT_DAVG'      ,'T2','DAY-AVG SFC DOWNWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 425
 'TEMPK_DAVG'       ,'T2','DAY-AVG SFC TEMPERATURE',' (K)'                  ,& ! 426
 'TEMPK_DMIN'       ,'T2','DAY-MIN SFC TEMPERATURE',' (K)'                  ,& ! 427
 'TEMPK_DMAX'       ,'T2','DAY-MAX SFC TEMPERATURE',' (K)'                  ,& ! 428
 'ACCPMIC_DTOT'     ,'T2','DAY-ACCUM MICPHYS PRECIP',' (kg m:S2:-2  )'      ,& ! 429
 'ACCPCON_DTOT'     ,'T2','DAY-ACCUM CUPARM PRECIP',' (kg m:S2:-2  )'       ,& ! 430
 'ACCPBOTH_DTOT'    ,'T2','DAY-ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )'  ,& ! 431
 'PRESS_UL_DAVG'      ,'T2','DAY-AVG UL PRESSURE',' (hPa)'                  ,& ! 432
 'ZONAL_WINDW_UL_DAVG','T2','DAY-AVG UL ZONAL WIND',' (m s:S2:-1  )'        ,& ! 433
 'MERID_WINDW_UL_DAVG','T2','DAY-AVG UL MERID WIND',' (m s:S2:-1  )'        ,& ! 434
 'CANTEMPK_DAVG'    ,'B2','DAY-AVG CANOPY TEMPERATURE',' (K)'               ,& ! 435
 'CANTEMPK_DMIN'    ,'B2','DAY-MIN CANOPY TEMPERATURE',' (K)'               ,& ! 436
 'CANTEMPK_DMAX'    ,'B2','DAY-MAX CANOPY TEMPERATURE',' (K)'               ,& ! 437
 'VEGTEMPK_DAVG'    ,'L2','DAY-AVG VEG TEMPERATURE',' (K)'                  ,& ! 438
 'VEGTEMPK_DMIN'    ,'L2','DAY-MIN VEG TEMPERATURE',' (K)'                  ,& ! 439
 'VEGTEMPK_DMAX'    ,'L2','DAY-MAX VEG TEMPERATURE',' (K)'                  ,& ! 440
 'SOILTEMPK_DAVG'   ,'L2','DAY-AVG SOIL TEMPERATURE',' (K)'                 ,& ! 441
 'SOILTEMPK_DMIN'   ,'L2','DAY-MIN SOIL TEMPERATURE',' (K)'                 ,& ! 442
 'SOILTEMPK_DMAX'   ,'L2','DAY-MAX SOIL TEMPERATURE',' (K)'                 ,& ! 443
 'LS_AIRTEMPK_DAVG' ,'B2','L/S DAY-AVG ATM TEMP',' (K)'                     ,& ! 444
 'LS_AIRTEMPK_DMIN' ,'B2','L/S DAY-MIN ATM TEMP',' (K)'                     ,& ! 445
 'LS_AIRTEMPK_DMAX' ,'B2','L/S DAY-MAX ATM TEMP',' (K)'                     ,& ! 446
 'LS_CANTEMPK_DAVG' ,'B2','L/S DAY-AVG CAN TEMP',' (K)'                     ,& ! 447
 'LS_CANTEMPK_DMIN' ,'B2','L/S DAY-MIN CAN TEMP',' (K)'                     ,& ! 448
 'LS_CANTEMPK_DMAX' ,'B2','L/S DAY-MAX CAN TEMP',' (K)'                     ,& ! 449
 'LS_SENSFLUX_DAVG' ,'B2','L/S DAY-AVG SENS HEAT FLUX',' (W m:S2:-2  )'     ,& ! 450
 'LS_LATFLUX_DAVG'  ,'B2','L/S DAY-AVG LAT HEAT FLUX',' (W m:S2:-2  )'      ,& ! 451
 'SENSFLUX_DAVG'    ,'T2','DAY-AVG SENS HEAT FLUX',' (W m:S2:-2  )'         ,& ! 452
 'LATFLUX_DAVG'     ,'T2','DAY-AVG LAT HEAT FLUX',' (W m:S2:-2  )'          ,& ! 453
 'PRESS_AVG24'      ,'T2','MONTH-AVG HOURLY PRESSURE',' (hPa)'               / ! 454

! Miscellaneous and new additions

data fldlib(1:4,471:493)/ &
 'RHO_OBS'       ,'T3' ,'NUDGING OBS AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 471
 'THETA_OBS'     ,'T3' ,'NUDGING OBS THETA',' (K)'                          ,& ! 472
 'SHW_OBS'       ,'T3' ,'NUDGING OBS VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 473
 'UZONAL_OBS'    ,'T3' ,'NUDGING OBS ZONAL WIND',' (m s:S2:-1  )'           ,& ! 474
 'UMERID_OBS'    ,'T3' ,'NUDGING OBS MERID WIND',' (m s:S2:-1  )'           ,& ! 475
 'RHO_SIM'       ,'T3' ,'NUDGING SIM AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 476
 'THETA_SIM'     ,'T3' ,'NUDGING SIM THETA',' (K)'                          ,& ! 477
 'SHW_SIM'       ,'T3' ,'NUDGING SIM VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 478
 'UZONAL_SIM'    ,'T3' ,'NUDGING SIM ZONAL WIND',' (m s:S2:-1  )'           ,& ! 479
 'UMERID_SIM'    ,'T3' ,'NUDGING SIM MERID WIND',' (m s:S2:-1  )'           ,& ! 480
 'RHO_OBS_SIM'   ,'T3' ,'NUDGING DIF AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 481
 'THETA_OBS_SIM' ,'T3' ,'NUDGING DIF THETA',' (K)'                          ,& ! 482
 'SHW_OBS_SIM'   ,'T3' ,'NUDGING DIF VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 483
 'UZONAL_OBS_SIM','T3' ,'NUDGING DIF ZONAL WIND',' (m s:S2:-1  )'           ,& ! 484
 'UMERID_OBS_SIM','T3' ,'NUDGING DIF MERID WIND',' (m s:S2:-1  )'           ,& ! 485
 'VXE'           ,'T3' ,'EARTH CARTESIAN X WIND',' (m s:S2:-1  )'           ,& ! 486
 'VYE'           ,'T3' ,'EARTH CARTESIAN Y WIND',' (m s:S2:-1  )'           ,& ! 487
 'VZE'           ,'T3' ,'EARTH CARTESIAN Z WIND',' (m s:S2:-1  )'           ,& ! 488
 'PBLH'          ,'T2' ,'PBL HEIGHT',' (m)'                                 ,& ! 489
 'HKH'           ,'T3' ,'EDDY DIFFUSIVITY',' (m:S2:2 s:S2:-1  )'            ,& ! 490
 'SHW_HCONV'     ,'T2' ,'TOTAL WATER HORIZ CONV',' (kg m:S2:-2   s:S2:-1  )',& ! 491
 'SHV_HCONV'     ,'T2' ,'WATER VAPOR HORIZ CONV',' (kg m:S2:-2   s:S2:-1  )',& ! 492
 'CLDNUM'        ,'T2' ,'CLOUD # CONCEN (GEOG)',' (# mg:S2:-1  )'            / ! 493

! External fields

data fldlib(1:4,501:503)/ &
 'VORTP'         ,'P3' ,'VORTP',' (s:S2:-1  )'                              ,& ! 501
 'VORTN'         ,'N3' ,'VORTN',' (s:S2:-1  )'                              ,& ! 502
 'RKE'           ,'T3' ,'RKE',' (s:S2:-1  )'                                 / ! 503

if (icall /= 1) then
   icall = 1
   allocate (aux(mwa))
endif

k = kk
i = ii

if (infotyp == 'UNITS') then

! Print field name to be plotted

   if (myrank == 0) write(io6,*) 'oplib ',trim(fldname0)

! Default values

   fldname = trim(fldname0)
   op%label = ' '
   op%units = ' '
   indp = 0

! Search through fldname0 for occurrence of '(' character
 
   lenstr = len_trim(fldname)

   do ic = 1,lenstr

      if (fldname(ic:ic) == '(') then
         fldname = trim(fldname0(1:ic-1))

         if (fldname0(lenstr:lenstr) /= ')' .or. lenstr - ic < 2) go to 1000

         if (fldname0(lenstr-1:lenstr-1) == '1')   indp = indp + 1
         if (fldname0(lenstr-1:lenstr-1) == '2')   indp = indp + 2
         if (fldname0(lenstr-1:lenstr-1) == '3')   indp = indp + 3
         if (fldname0(lenstr-1:lenstr-1) == '4')   indp = indp + 4
         if (fldname0(lenstr-1:lenstr-1) == '5')   indp = indp + 5
         if (fldname0(lenstr-1:lenstr-1) == '6')   indp = indp + 6
         if (fldname0(lenstr-1:lenstr-1) == '7')   indp = indp + 7
         if (fldname0(lenstr-1:lenstr-1) == '8')   indp = indp + 8
         if (fldname0(lenstr-1:lenstr-1) == '9')   indp = indp + 9

         if (lenstr - ic > 2) then

            if (fldname0(lenstr-2:lenstr-2) == '1') indp = indp + 10
            if (fldname0(lenstr-2:lenstr-2) == '2') indp = indp + 20
            if (fldname0(lenstr-2:lenstr-2) == '3') indp = indp + 30
            if (fldname0(lenstr-2:lenstr-2) == '4') indp = indp + 40
            if (fldname0(lenstr-2:lenstr-2) == '5') indp = indp + 50
            if (fldname0(lenstr-2:lenstr-2) == '6') indp = indp + 60
            if (fldname0(lenstr-2:lenstr-2) == '7') indp = indp + 70
            if (fldname0(lenstr-2:lenstr-2) == '8') indp = indp + 80
            if (fldname0(lenstr-2:lenstr-2) == '9') indp = indp + 90

         endif

         exit
      endif

   enddo

! Search through FLDLIB data list for match to FLDNAME

   do ifield = 1,nfields
      icase = ifield

      if (trim(fldname) == fldlib(1,ifield)) exit

      if (ifield == nfields) then
         write(io6,*) 'Plot field ',trim(fldname),' not available in oplot_lib.'
         go to 1000
      endif
   enddo

   op%stagpt = fldlib(2,icase)(1:1)
   op%dimens = fldlib(2,icase)(2:3)
   op%label  = fldlib(3,icase)
   if (indp > 0) op%label = trim(fldname0)  
   op%units  = fldlib(4,icase)

   op%fldval_min = 1.e30
   op%fldval_max = -1.e30
   op%fldvalv_min = 1.e30
   op%fldvalv_max = -1.e30

   if (op%stagpt == 'T' .or. op%stagpt == 'W') then
      i = jtab_w(jtw_prog)%iw(1)
      k = lpw(i)
   elseif (op%stagpt == 'V') then
      i = jtab_v(jtv_prog)%iv(1)
      k = lpv(i)
   elseif (op%stagpt == 'M' .or. op%stagpt == 'P') then
      i = jtab_m(jtm_vadj)%im(1)
      k = lpm(i)
   elseif (op%stagpt == 'L' .or. op%stagpt == 'S' .or. op%stagpt == 'B') then
      i = 2
      k = 1
   endif
endif

if (infotyp == 'VALUV') then
   op%stagpt = 'V'
   op%dimens = '3'
   icase = 4
endif

if (op%stagpt == 'L' .and. mwl < 2) then
   notavail = 3
   return
endif

if (op%stagpt == 'S' .and. mws < 2) then
   notavail = 3
   return
endif

notavail = 0
kp = min(k+1,mza)

! Execute IF block below even when infotyp == 'UNITS'
! in order to check whether current plot field is available in this model run.
! For this type of call to this subroutine, (k,i) are passed in as (1,2).

!-----------------------------------------
! ATMOSPHERE - 3D
!-----------------------------------------

select case(icase)

case(1) ! 'VMC'

   if (.not. allocated(vmc)) go to 1000

   fldval = wtbot * vmc(k ,i) &
          + wttop * vmc(kp,i)

case(2) ! 'WMC'

! Need to re-examine use of k-1 when k = lpw

   fldval = wtbot * wmc(k ,i) &
          + wttop * wmc(kp,i)

case(3) ! 'VMP'

   if (.not. allocated(vmp)) go to 1000

   fldval = wtbot * vmp(k ,i) &
          + wttop * vmp(kp,i)

case(4) ! 'VC'

   if (.not. allocated(vc)) go to 1000

   if (k >= lpv(i)) then

      fldval = wtbot * vc(k ,i) &
             + wttop * vc(kp,i)

   else

      iw1 = itab_v(i)%iw(1)
      iw2 = itab_v(i)%iw(2)

      if (k >= lpw(iw1)) then
         fldval = vnx(i) * vxe(k,iw1) &
                + vny(i) * vye(k,iw1) &
                + vnz(i) * vze(k,iw1)
      elseif (k >= lpw(iw2)) then
         fldval = vnx(i) * vxe(k,iw2) &
                + vny(i) * vye(k,iw2) &
                + vnz(i) * vze(k,iw2)
      endif

   endif

case(5) ! 'WC'

   fldval = wtbot * wc(k ,i) &
          + wttop * wc(kp,i)

case(6) ! 'RHO'

   fldval = wtbot * rho(k ,i) &
          + wttop * rho(kp,i)

case(7) ! 'PRESS'

   fldval = press(k,i) * .01

   if (nl%test_case == 2 .or. nl%test_case == 5) then
      fldval = press(k,i)
   endif

case(8) ! 'THIL'

   fldval = wtbot * thil(k ,i) &
          + wttop * thil(kp,i)

case(9) ! 'THETA'

   fldval = wtbot * theta(k ,i) &
          + wttop * theta(kp,i)

case(10) ! 'AIRTEMPK'

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i)

case(11) ! 'AIRTEMPC'

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i) - 273.15

case(12) ! 'SH_W'

   fldval = (wtbot * sh_w(k ,i) &
          +  wttop * sh_w(kp,i)) * 1.e3

case(13) ! 'SH_V'

   fldval = (wtbot * sh_v(k ,i) &
          +  wttop * sh_v(kp,i)) * 1.e3

case(14) ! 'SH_C'

   if (.not. allocated(sh_c)) go to 1000

   fldval = (wtbot * sh_c(k ,i) &
          +  wttop * sh_c(kp,i)) * 1.e3

case(15) ! 'SH_D'

   if (.not. allocated(sh_d)) go to 1000

   fldval = (wtbot * sh_d(k ,i) &
          +  wttop * sh_d(kp,i)) * 1.e3

case(16) ! 'SH_R'

   if (.not. allocated(sh_r)) go to 1000

   fldval = (wtbot * sh_r(k ,i) &
          +  wttop * sh_r(kp,i)) * 1.e3

case(17) ! 'SH_P'

   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbot * sh_p(k ,i) &
          +  wttop * sh_p(kp,i)) * 1.e3

case(18) ! 'SH_S'

   if (.not. allocated(sh_s)) go to 1000

   fldval = (wtbot * sh_s(k ,i) &
          +  wttop * sh_s(kp,i)) * 1.e3

case(19) ! 'SH_A'

   if (.not. allocated(sh_a)) go to 1000

   fldval = (wtbot * sh_a(k ,i) &
          +  wttop * sh_a(kp,i)) * 1.e3

case(20) ! 'SH_G'

   if (.not. allocated(sh_g)) go to 1000

   fldval = (wtbot * sh_g(k ,i) &
          +  wttop * sh_g(kp,i)) * 1.e3

case(21) ! 'SH_H'

   if (.not. allocated(sh_h)) go to 1000

   fldval = (wtbot * sh_h(k ,i) &
          +  wttop * sh_h(kp,i)) * 1.e3

case(22) ! 'SH_CP'

   if (.not. allocated(sh_c)) go to 1000
   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbot * (sh_c(k ,i) + sh_p(k ,i)) &
          +  wttop * (sh_c(kp,i) + sh_p(kp,i))) * 1.e3

case(23) ! 'SH_TOTCOND'

   fldval = (wtbot * (sh_w(k ,i) - sh_v(k ,i)) &
          +  wttop * (sh_w(kp,i) - sh_v(kp,i))) * 1.e3

case(24) ! 'CON_C'

   if (.not. allocated(con_c)) go to 1000

   fldval = (wtbot * con_c(k ,i) &
          +  wttop * con_c(kp,i)) * 1.e-6

case(25) ! 'CON_D'

   if (.not. allocated(con_d)) go to 1000

   fldval = (wtbot * con_d(k ,i) &
          +  wttop * con_d(kp,i)) * 1.e-6

case(26) ! 'CON_R'

   if (.not. allocated(con_r)) go to 1000

   fldval = wtbot * con_r(k ,i) &
          + wttop * con_r(kp,i)

case(27) ! 'CON_P'

   if (.not. allocated(con_p)) go to 1000

   fldval = wtbot * con_p(k ,i) &
          + wttop * con_p(kp,i)

case(28) ! 'CON_S'

   if (.not. allocated(con_s)) go to 1000

   fldval = wtbot * con_s(k ,i) &
          + wttop * con_s(kp,i)

case(29) ! 'CON_A'

   if (.not. allocated(con_a)) go to 1000

   fldval = wtbot * con_a(k ,i) &
          + wttop * con_a(kp,i)

case(30) ! 'CON_G'

   if (.not. allocated(con_g)) go to 1000

   fldval = wtbot * con_g(k ,i) &
          + wttop * con_g(kp,i)

case(31) ! 'CON_H'

   if (.not. allocated(con_h)) go to 1000

   fldval = wtbot * con_h(k ,i) &
          + wttop * con_h(kp,i)

case(32) ! 'CON_CCN'

   if (.not. allocated(con_ccn)) go to 1000

   fldval = (wtbot * con_ccn(k ,i) &
          +  wttop * con_ccn(kp,i)) * 1.e-6

case(33) ! 'CON_GCCN'

   if (.not. allocated(con_gccn)) go to 1000

   fldval = (wtbot * con_gccn(k ,i) &
          +  wttop * con_gccn(kp,i)) * 1.e-6

case(34) ! 'CON_IFN'

   if (.not. allocated(con_ifn)) go to 1000

   fldval = wtbot * con_ifn(k ,i) &
          + wttop * con_ifn(kp,i)

case(35) ! 'HKM'

   fldval = wtbot * hkm(k ,i) &
          + wttop * hkm(kp,i)

case(36) ! 'FTHRD'

   if (.not. allocated(fthrd_sw)) goto 1000
   if (.not. allocated(fthrd_lw)) goto 1000

   fldval = wtbot * (fthrd_sw(k ,i) + fthrd_lw(k ,i)) &
          + wttop * (fthrd_sw(kp,i) + fthrd_lw(kp,i))

case(37:40) ! 'SPEEDW','AZIMW','ZONAL_WINDW','MERID_WINDW'

   npoly = itab_w(i)%npoly
   rpolyi = 1. / real(npoly)

   vx = wtbot * vxe(k ,i) &
      + wttop * vxe(kp,i)

   vy = wtbot * vye(k ,i) &
      + wttop * vye(kp,i)

   vz = wtbot * vze(k ,i) &
      + wttop * vze(kp,i)

   if (trim(fldname) == 'SPEEDW') then
      fldval = sqrt(vx**2 + vy**2 + vz**2)
   else
      
      if (mdomain < 2) then

         raxis = sqrt(xew(i)**2 + yew(i)**2)  ! dist from earth axis

         if (raxis > 1.e3) then
            u = (vy * xew(i) - vx * yew(i)) / raxis
            v = vz * raxis / erad &
              - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 
         else
            u = 0.
            v = 0.
         endif

      else
         u = vx
         v = vy
      endif

      if (trim(fldname) == 'AZIMW') then
         fldval = mod(450. - piu180 * atan2(v,u),360.)
      elseif (trim(fldname) == 'ZONAL_WINDW') then
         fldval = u
      elseif (trim(fldname) == 'MERID_WINDW') then
         fldval = v
      endif
   endif

case(41:43) ! 'RVORTZM','TVORTZM','RVORTZM_P'

   fldval = 0.

   do j = 1,itab_m(i)%npoly


      iv = itab_m(i)%iv(j)

      if (k >= lpv(iv)) then

         vcc = wtbot * vc(k ,iv) &
             + wttop * vc(kp,iv)

      else

         iw1 = itab_v(iv)%iw(1)
         iw2 = itab_v(iv)%iw(2)

         if (k >= lpw(iw1)) then
            vcc = vnx(iv) * vxe(k,iw1) &
                + vny(iv) * vye(k,iw1) &
                + vnz(iv) * vze(k,iw1)
         elseif (k >= lpw(iw2)) then
            vcc = vnx(iv) * vxe(k,iw2) &
                + vny(iv) * vye(k,iw2) &
                + vnz(iv) * vze(k,iw2)
         else
            vcc = 0.
         endif

      endif

      if (fldname == 'RVORTZM_P') then

         vcc_init = wtbot * vc_init(k ,iv) &
                  + wttop * vc_init(kp,iv)

         vcc = vcc - vcc_init

      endif

! Now reconstruct total wind vector projected onto vector from IW1 to IW2

      if (i == itab_v(iv)%im(2)) then
         fldval = fldval + dnv(iv) * vcc
      else
         fldval = fldval - dnv(iv) * vcc
      endif

   enddo

   fldval = fldval / arm0(i)

! For shallow water test case 2, subtract initial vorticity to get perturbation

   if (trim(fldname) == 'RVORTZM' .and. nl%test_case == 2) then
      fldval = fldval - omega2 * zem(i) / (12. * erad)
   endif

   if (trim(fldname) == 'TVORTZM') then
      fldval = fldval + omega2 * zem(i) / erad  ! add earth vorticity at M point
   endif

case(44) ! 'DIVERG'

   fldval = 0.

   npoly = itab_w(i)%npoly

   if (itab_w(i)%iwp == i) then

      do j = 1,npoly

         iv = itab_w(i)%iv(j)

         fldval = fldval &

                + wtbot * vmc(k,iv) * (-itab_w(i)%dirv(j)) * arv(k,iv) &
                        / (volt(k,i) * rho(k,i)) &

                + wttop * vmc(kp,iv) * (-itab_w(i)%dirv(j)) * arv(kp,iv) &
                        / (volt(kp,i) * rho(kp,i))
      enddo
   endif

case(45) ! 'VMASSFLUX'

   fldval = vmc(k,i) * arv(k,i)

case(46) ! 'VC_P'

   fldval = vc(k,i) - vc_init(k,i)

case(47) ! 'PRESS_P'

   fldval = press(k,i) - press_init(k,i)

case(48) ! 'RHO_P'

   fldval = wtbot * (rho(k ,i) - rho_init(k ,i)) &
          + wttop * (rho(kp,i) - rho_init(kp,i))

! For shallow water test case 5, define RHO_P from reference fields

   if (nl%test_case == 5) then
      call npr_bicubics(zanal00_swtc5,glatw(i),glonw(i),zanal0_swtc5)
      call npr_bicubics(zanal15_swtc5,glatw(i),glonw(i),zanal_swtc5)

      fldval = (rho(k,i) - rho_init(k,i)) - (zanal_swtc5 - zanal0_swtc5)
   endif

case(49) ! 'THETA_P'

   fldval = wtbot * (theta(k ,i) - theta_init(k ,i)) &
          + wttop * (theta(kp,i) - theta_init(kp,i))

case(50) ! 'AIRTEMPK_P'

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i) &
          - wtbot * theta_init(k ,i) * (press_init(k ,i) * p00i) ** rocp &
          - wttop * theta_init(kp,i) * (press_init(kp,i) * p00i) ** rocp

case(51) ! 'VMT'

   fldval = vmt(k,i)

case(52) ! 'WMT'

   fldval = wmt(k,i)

case(53) ! 'ADDSC'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000

   fldval = wtbot * addsc(indp)%sclp(k ,i) &
          + wttop * addsc(indp)%sclp(kp,i)

case(54) ! 'ADDSC_P'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000
   if (.not. allocated(addsc1_init)) go to 1000

   fldval = addsc(indp)%sclp(k,i) - addsc1_init(k,i) 
   
!-----------------------------------------
! ATMOSPHERE - 2D
!-----------------------------------------

case(55) ! 'ZPLEV'

   fldval = wtbot * zt(k ) &
          + wttop * zt(kp)

case(56) ! 'QWCON'

   fldval = (wtbot * qwcon(k ,i) &
          +  wttop * qwcon(kp,i)) * 1.e3

case(62) ! 'RSHORT_TOP'

   if (.not. allocated(rshort_top)) go to 1000

   fldval = rshort_top(i)

case(63) ! 'RSHORTUP_TOP'

   if (.not. allocated(rshortup_top)) go to 1000

   fldval = rshortup_top(i)

case(64) ! 'RLONGUP_TOP'

   if (.not. allocated(rlongup_top)) go to 1000

   fldval = rlongup_top(i)

!-----------------------------------------
! ATMOSPHERE SURFACE - 2D
!-----------------------------------------

case(65) ! 'RSHORT'

   if (.not. allocated(rshort)) go to 1000

   fldval = rshort(i)

case(66) ! 'RSHORTUP'

   if (.not. allocated(rshort)) go to 1000
   if (.not. allocated(albedt)) go to 1000

   fldval = rshort(i) * albedt(i)

case(67) ! 'RLONG'

   if (.not. allocated(rlong)) go to 1000

   fldval = rlong(i)

case(68) ! 'RLONGUP'

   if (.not. allocated(rlongup)) go to 1000

   fldval = rlongup(i)

case(69) ! 'ALBEDT'

   if (.not. allocated(albedt)) go to 1000

   fldval = albedt(i)

case(70) ! 'VKM_SFC'

   fldval = vkm_sfc(1,i)

case(71) ! 'USTAR'

   fldval = ustar(i)

case(72) ! 'SENSFLUX'

   fldval = sfluxt(i) * cp
   
case(73) ! 'VAPFLUX'

   fldval = sfluxr(i)

case(74) ! 'LATFLUX'

   fldval = sfluxr(i) * alvl

case(75) ! 'PCPRD'

   if (.not. allocated(pcprd)) go to 1000

   fldval = pcprd(i) * 3600.

case(76) ! 'PCPRR'

   if (.not. allocated(pcprr)) go to 1000

   fldval = pcprr(i) * 3600.

case(77) ! 'PCPRP'

   if (.not. allocated(pcprp)) go to 1000

   fldval = pcprp(i)  * 3600.   

case(78) ! 'PCPRS'

   if (.not. allocated(pcprs)) go to 1000

   fldval = pcprs(i) * 3600.

case(79) ! 'PCPRA'

   if (.not. allocated(pcpra)) go to 1000

   fldval = pcpra(i) * 3600.

case(80) ! 'PCPRG'

   if (.not. allocated(pcprg)) go to 1000

   fldval = pcprg(i) * 3600.

case(81) ! 'PCPRH'

   if (.not. allocated(pcprh)) go to 1000

   fldval = pcprh(i) * 3600.

case(82) ! 'PCPRMIC'

   fldval = 0.

   if (allocated(pcprd)) fldval = fldval + pcprd(i)
   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)

   fldval = fldval * 3600.

case(83) ! 'PCPRCON'

   if (.not. allocated(conprr)) go to 1000

   fldval = conprr(i) * 3600.

case(84) ! 'PCPRBOTH'

   fldval = 0.

   if (allocated(pcprd)) fldval = fldval + pcprd(i)
   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)
   if (allocated(conprr)) fldval = fldval + conprr(i)

   fldval = fldval * 3600.

case(85) ! 'ACCPD'

   if (.not. allocated(accpd)) go to 1000

   fldval = real(accpd(i))

case(86) ! 'ACCPR'

   if (.not. allocated(accpr)) go to 1000

   fldval = real(accpr(i))

case(87) ! 'ACCPP'

   if (.not. allocated(accpp)) go to 1000

   fldval = real(accpp(i))

case(88) ! 'ACCPS'

   if (.not. allocated(accps)) go to 1000

   fldval = real(accps(i))

case(89) ! 'ACCPA'

   if (.not. allocated(accpa)) go to 1000

   fldval = real(accpa(i))

case(90) ! 'ACCPG'

   if (.not. allocated(accpg)) go to 1000

   fldval = real(accpg(i))

case(91) ! 'ACCPH'

   if (.not. allocated(accph)) go to 1000

   fldval = real(accph(i))

case(92) ! 'ACCPMIC'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + real(accpd(i))
   if (allocated(accpr)) fldval = fldval + real(accpr(i))
   if (allocated(accpp)) fldval = fldval + real(accpp(i))
   if (allocated(accps)) fldval = fldval + real(accps(i))
   if (allocated(accpa)) fldval = fldval + real(accpa(i))
   if (allocated(accpg)) fldval = fldval + real(accpg(i))
   if (allocated(accph)) fldval = fldval + real(accph(i))

case(93) ! 'ACCPCON'

   if (.not. allocated(aconpr)) go to 1000

   fldval = real(aconpr(i))

case(94) ! 'ACCPBOTH'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + real(accpd(i))
   if (allocated(accpr)) fldval = fldval + real(accpr(i))
   if (allocated(accpp)) fldval = fldval + real(accpp(i))
   if (allocated(accps)) fldval = fldval + real(accps(i))
   if (allocated(accpa)) fldval = fldval + real(accpa(i))
   if (allocated(accpg)) fldval = fldval + real(accpg(i))
   if (allocated(accph)) fldval = fldval + real(accph(i))
   if (allocated(aconpr)) fldval = fldval + real(aconpr(i))

case(95) ! 'WSTAR'

   fldval = wstar(i)

!-----------------------------------------
! ATMOSPHERE DIF2 fields - (3D & 2D)
!-----------------------------------------

case(101:102) ! 'ZONAL_WINDW_DIF2', 'MERID_WINDW_DIF2'

   if (.not. allocated(vc_accum)) go to 1000

   npoly = itab_w(i)%npoly

   wc_change = 0.5 * (wc_accum_prev0(k-1,i) - wc_accum_prev1(k-1,i) &
                    + wc_accum_prev0(k  ,i) - wc_accum_prev1(k  ,i))

   vx = wc_change * wnx(i)
   vy = wc_change * wny(i)
   vz = wc_change * wnz(i)

   do jv = 1, npoly
      iv = itab_w(i)%iv(jv)

      vc_change = vc_accum_prev0(k,iv) - vc_accum_prev1(k,iv)

      vx = vx + itab_w(i)%ecvec_vx(jv) * vc_change
      vy = vy + itab_w(i)%ecvec_vy(jv) * vc_change
      vz = vz + itab_w(i)%ecvec_vz(jv) * vc_change
   enddo

   if (mdomain < 2) then

      raxis = sqrt(xew(i)**2 + yew(i)**2)  ! dist from earth axis

      if (raxis > 1.e3) then
         u = (vy * xew(i) - vx * yew(i)) / raxis
         v = vz * raxis / erad &
           - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 
      else
         u = 0.
         v = 0.
      endif

   else
      u = vx
      v = vy
   endif

   if (trim(fldname) == 'ZONAL_WINDW_DIF2') then
      fldval = u
   elseif (trim(fldname) == 'MERID_WINDW_DIF2') then
      fldval = v
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(103) ! 'WC_DIF2'

   if (.not. allocated(wc_accum)) go to 1000

   fldval = wc_accum_prev0(k,i) - wc_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(104) ! 'PRESS_DIF2'

   if (.not. allocated(press_accum)) go to 1000

   fldval = (press_accum_prev0(k,i) - press_accum_prev1(k,i)) * .01

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(105) ! 'AIRTEMPK_DIF2'

   if (.not. allocated(tair_accum)) go to 1000

   fldval = tair_accum_prev0(k,i) - tair_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(106) ! 'SH_V_DIF2'

   if (.not. allocated(sh_v_accum)) go to 1000

   fldval = sh_v_accum_prev0(k,i) - sh_v_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 1.e3 / (time8_prev0 - time8_prev1)
   endif

case(107) ! 'PCPMIC_DIF2'

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   fldval = accpmic_prev0(i) - accpmic_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(108) ! 'PCPCON_DIF2'

   if (.not. allocated(aconpr)) go to 1000

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   fldval = accpcon_prev0(i) - accpcon_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(109) ! 'PCPBOTH_DIF2'

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   accpboth_prev0 = accpmic_prev0(i) + accpcon_prev0(i)
   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)

   fldval = accpboth_prev0 - accpboth_prev1

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(110) ! 'RSHORT_DIF2'

   if (.not. allocated(rshort_accum)) go to 1000

   fldval = rshort_accum_prev0(i) - rshort_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(111) ! 'RSHORTUP_DIF2'

   if (.not. allocated(rshortup_accum)) go to 1000

   fldval = rshortup_accum_prev0(i) - rshortup_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(112) ! 'RLONG_DIF2'

   if (.not. allocated(rlong_accum)) go to 1000

   fldval = rlong_accum_prev0(i) - rlong_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(113) ! 'RLONGUP_DIF2'

   if (.not. allocated(rlongup_accum)) go to 1000

   fldval = rlongup_accum_prev0(i) - rlongup_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(114) ! 'RSHORT_TOP_DIF2'

   if (.not. allocated(rshort_top_accum)) go to 1000

   fldval = rshort_top_accum_prev0(i) - rshort_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(115) ! 'RSHORTUP_TOP_DIF2'

   if (.not. allocated(rshortup_top_accum)) go to 1000

   fldval = rshortup_top_accum_prev0(i) - rshortup_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(116) ! 'RLONGUP_TOP_DIF2'

   if (.not. allocated(rlongup_top_accum)) go to 1000

   fldval = rlongup_top_accum_prev0(i) - rlongup_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(117) ! 'SENSFLUX_DIF2'

   if (.not. allocated(sfluxt_accum)) go to 1000

   fldval = (sfluxt_accum_prev0(i) - sfluxt_accum_prev1(i)) * cp

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(118) ! 'LATFLUX_DIF2'

   if (.not. allocated(sfluxr_accum)) go to 1000

   fldval = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)) * alvl

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(119) ! 'VAPFLUX_DIF2'

   if (.not. allocated(sfluxr_accum)) go to 1000

   fldval = sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(120) ! 'RSHORT_CLR_DIF2'

   if (.not. allocated(rshort_clr_accum)) go to 1000

   fldval = rshort_clr_accum_prev0(i) - rshort_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(121) ! 'RSHORTUP_CLR_DIF2'

   if (.not. allocated(rshortup_clr_accum)) go to 1000

   fldval = rshortup_clr_accum_prev0(i) - rshortup_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(122) ! 'RLONG_CLR_DIF2'

   if (.not. allocated(rlong_clr_accum)) go to 1000

   fldval = rlong_clr_accum_prev0(i) - rlong_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(123) ! 'RLONGUP_CLR_DIF2'

   if (.not. allocated(rlongup_clr_accum)) go to 1000

   fldval = rlongup_clr_accum_prev0(i) - rlongup_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(124) ! 'RSHORT_TOP_CLR_DIF2'

   if (.not. allocated(rshort_top_clr_accum)) go to 1000

   fldval = rshort_top_clr_accum_prev0(i) - rshort_top_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(125) ! 'RSHORTUP_TOP_CLR_DIF2'

   if (.not. allocated(rshortup_top_clr_accum)) go to 1000

   fldval = rshortup_top_clr_accum_prev0(i) - rshortup_top_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(126) ! 'RLONGUP_TOP_CLR_DIF2'

   if (.not. allocated(rlongup_top_clr_accum)) go to 1000

   fldval = rlongup_top_clr_accum_prev0(i) - rlongup_top_clr_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(127) ! 'PCPMIC_DIF4'

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 second

   fldval1 = accpmic_prev0(i) - accpmic_prev1(i)
   fldval2 = accpmic_prev2(i) - accpmic_prev3(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

!   fldval =  accpmic_prev0(i) - accpmic_prev1(i) &
!          - (accpmic_prev2(i) - accpmic_prev3(i))

!   if (abs(time8_prev0 - time8_prev2) > .99) then
!      fldval = fldval * 86400. / (time8_prev0 - time8_prev2)
!   endif

case(128) ! 'PCPCON_DIF4'

   if (.not. allocated(aconpr)) go to 1000

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 second

   fldval1 = accpcon_prev0(i) - accpcon_prev1(i)
   fldval2 = accpcon_prev2(i) - accpcon_prev3(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

!   fldval =  accpcon_prev0(i) - accpcon_prev1(i) &
!          - (accpcon_prev2(i) - accpcon_prev3(i))

!   if (abs(time8_prev0 - time8_prev2) > .99) then
!      fldval = fldval * 86400. / (time8_prev0 - time8_prev2)
!   endif

case(129) ! 'PCPBOTH_DIF4'

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 second

   accpboth_prev0 = accpmic_prev0(i) + accpcon_prev0(i)
   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)
   accpboth_prev2 = accpmic_prev2(i) + accpcon_prev2(i)
   accpboth_prev3 = accpmic_prev3(i) + accpcon_prev3(i)

   fldval1 = accpboth_prev0 - accpboth_prev1
   fldval2 = accpboth_prev2 - accpboth_prev3

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

!   fldval =  accpboth_prev0 - accpboth_prev1 &
!          - (accpboth_prev2 - accpboth_prev3)

!   if (abs(time8_prev0 - time8_prev2) > .99) then
!      fldval = fldval * 86400. / (time8_prev0 - time8_prev2)
!   endif

case(130) ! 'PCPMIC_REL4'

! Compute relative differences involving 3 previously-stored fields
! [fldval=runB_end; prev1=runB_beg; prev2=runA_end; prev3=runA_beg]

   denom  = (accpmic_prev0(i) - accpmic_prev1(i)) &
          + (accpmic_prev2(i) - accpmic_prev3(i))

   fldval = (accpmic_prev0(i) - accpmic_prev1(i)) &
          - (accpmic_prev2(i) - accpmic_prev3(i))

   if (abs(denom) > 1.e-3) then
      fldval = max(-1.0,min(1.0,fldval / denom)) ! min/max in case of trunc error
   else
      fldval = 0.
   endif

case(131) ! 'PCPCON_REL4'

   if (.not. allocated(aconpr)) go to 1000

! Compute relative differences involving 3 previously-stored fields

   denom  = (accpcon_prev0(i) - accpcon_prev1(i)) &
          + (accpcon_prev2(i) - accpcon_prev3(i))

   fldval = (accpcon_prev0(i) - accpcon_prev1(i)) &
          - (accpcon_prev2(i) - accpcon_prev3(i))

   if (abs(denom) > 1.e-3) then
      fldval = max(-1.0,min(1.0,fldval / denom)) ! min/max in case of trunc error
   else
      fldval = 0.
   endif

case(132) ! 'PCPBOTH_REL4'

! Compute relative differences involving 3 previously-stored fields

   accpboth_prev0 = accpmic_prev0(i) + accpcon_prev0(i)
   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)
   accpboth_prev2 = accpmic_prev2(i) + accpcon_prev2(i)
   accpboth_prev3 = accpmic_prev3(i) + accpcon_prev3(i)

   denom  = (accpboth_prev0 - accpboth_prev1) &
          + (accpboth_prev2 - accpboth_prev3)

   fldval = (accpboth_prev0 - accpboth_prev1) &
          - (accpboth_prev2 - accpboth_prev3)

   if (abs(denom) > 1.e-3) then
      fldval = max(-1.0,min(1.0,fldval / denom)) ! min/max in case of trunc error
   else
      fldval = 0.
   endif

!-----------------------------------------
! LAND CELLS - 3D
!-----------------------------------------

case(151) ! 'SOIL_TEXT'

   fldval = real(land%ntext_soil(k,i))

case(152) ! 'SOIL_ENERGY'

   fldval = land%soil_energy(k,i) * 1.e-6

case(153) ! 'SOIL_TEMPK'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = tempk

case(154) ! 'SOIL_FRACLIQ'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = fracliq

case(155) ! 'SOIL_WATER'

   if (land%flag_vg(i)) then
      fldval = land%soil_water(k,i) * slmstsi_vg(land%ntext_soil(k,i))
   else
      fldval = land%soil_water(k,i) * slmstsi_ch(land%ntext_soil(k,i))
   endif

case(156) ! 'SFWAT_MASS'

   fldval = land%sfcwater_mass(k,i)

case(157) ! 'SFWAT_ENERGY'

   if (land%nlev_sfcwater(i) == 0) then
      notavail = 4
   else
      fldval = land%sfcwater_energy(k,i) * 1.e-3
   endif

case(158) ! 'SFWAT_TEMPK'

   if (land%nlev_sfcwater(i) == 0) then
      notavail = 4
   else
      call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
      fldval = tempk
   endif

case(159) ! 'SFWAT_FRACLIQ'

   if (land%nlev_sfcwater(i) == 0) then
      notavail = 4
   else
      call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
      fldval = fracliq
   endif

case(160) ! 'SFWAT_DEPTH'

   fldval = 0.
   do klev = 1,land%nlev_sfcwater(i)
      fldval = fldval + land%sfcwater_depth(klev,i)
   enddo

!-----------------------------------------
! LAND CELLS - 2D
!-----------------------------------------

case(161) ! 'NLEV_SFWAT'

   fldval = real(land%nlev_sfcwater(i))

case(162) ! 'VEG_NDVIC'

   fldval = land%veg_ndvic(i)

case(163) ! 'VEG_TEMPC'

   fldval = land%veg_temp(i) - 273.15

case(164) ! 'VEG_TEMPK'

   fldval = land%veg_temp(i)

case(165) ! 'VEG_WATER'

   fldval = land%veg_water(i)

case(166) ! 'STOM_RESIST'

   fldval = land%stom_resist(i)

case(167) ! 'SFCWATER_TOT'

   fldval = sum(land%sfcwater_mass(:,i))

case(168) ! 'SFCWATER_TOP_TEMP'

   if (land%nlev_sfcwater(i) == 0) then
      notavail = 4
   else
      nls = land%nlev_sfcwater(i)
      call qtk(land%sfcwater_energy(nls,i),tempk,fracliq)
      fldval = tempk
   endif

case(169) ! 'SOIL_TOP_TEMP'

   call qwtk(land%soil_energy(nzg,i)       &
            ,land%soil_water(nzg,i)*1.e3   &
            ,slcpd(land%ntext_soil(nzg,i)) &
            ,tempk, fracliq)
   fldval = tempk

case(170) ! 'GROUND_SHV'

   fldval = land%ground_shv(i) * 1.e3

case(171) ! 'SOIL_DEPTH'

   fldval = -slz(1)

case(172) ! 'SOIL_WATER_TOT'

   fldval = sum(land%soil_water(:,i) * dslz(:))

case(173) ! 'HEAD0'

   fldval = land%head0(i)

!-----------------------------------------
! LAND CELLS - DIF2 fields
!-----------------------------------------

case(174) ! 'L_WXFER1_DIF2'

   if (.not. allocated(wxfer1_l_accum)) go to 1000
   fldval = (wxfer1_l_accum_prev0(i) - wxfer1_l_accum_prev1(i)) * 1.e3

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

!-----------------------------------------
! LAND CELLS - ATM averages
!-----------------------------------------

case(175:176) ! 'AL_SFCWATER_TOT', 'AL_SOIL_WATER_TOT'

   fldval = 0.
   area_sum = 0.

   do jland = 1,itab_w(i)%nland
      iland = itab_w(i)%iland(jland)

      if (trim(fldname) == 'AL_SFCWATER_TOT') then
         fldval = fldval + sum(land%sfcwater_mass(:,iland)) * land%area(iland)
      elseif (trim(fldname) == 'AL_SOIL_WATER_TOT') then
         fldval = fldval + sum(land%soil_water(:,iland) * dslz(:)) * land%area(iland)
      endif

      area_sum = area_sum + land%area(iland)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

!-----------------------------------------
! LAND CELLS - ATM averages of DIF2 fields
!-----------------------------------------

case(177) ! 'AL_WXFER1_DIF2'

   if (.not. allocated(wxfer1_l_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do jland = 1,itab_w(i)%nland
      iland = itab_w(i)%iland(jland)

      fldval = fldval + (wxfer1_l_accum_prev0(iland) - wxfer1_l_accum_prev1(iland)) &
             * land%area(iland)

      area_sum = area_sum + land%area(iland)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum * 1.e3 ! convert to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

!-----------------------------------------
! SEA CELLS
!-----------------------------------------

case(178) ! 'SEATP'

   fldval = sea%seatp(i)   

case(179) ! 'SEATF'

   fldval = sea%seatf(i)

case(180) ! 'SEATC'

   fldval = sea%seatc(i)

case(181) ! 'SEAICEP'

   fldval = sea%seaicep(i)   

case(182) ! 'SEAICEF'

   fldval = sea%seaicef(i)

case(183) ! 'SEAICEC'

   fldval = sea%seaicec(i)

!-----------------------------------------
! LAND AND SEA CELLS - 2D
!-----------------------------------------

case(201) ! 'LEAF_CLASS'

   if (op%stagpt == 'S') then
      fldval = real(sea%leaf_class(i))
   elseif (op%stagpt == 'L') then
      fldval = real(land%leaf_class(i))
   endif

case(202) ! 'LS_IW'

   if (op%stagpt == 'S') then
      fldval = real(itab_ws(i)%iw)
   elseif (op%stagpt == 'L') then
      fldval = real(itab_wl(i)%iw)
   endif

case(203) ! 'LS_KW'

   if (op%stagpt == 'S') then
      fldval = real(itab_ws(i)%kw)
   elseif (op%stagpt == 'L') then
      fldval = real(itab_wl(i)%kw)
   endif

case(204) ! 'LS_ARF_IW'

   if (op%stagpt == 'S') then
      fldval = real(itab_ws(i)%arf_iw)
   elseif (op%stagpt == 'L') then
      fldval = real(itab_wl(i)%arf_iw)
   endif

case(205) ! 'LS_ARF_KW'

   if (op%stagpt == 'S') then
      fldval = real(itab_ws(i)%arf_kw)
   elseif (op%stagpt == 'L') then
      fldval = real(itab_wl(i)%arf_kw)
   endif

case(206) ! 'LS_AREA'

   if (op%stagpt == 'S') then
      fldval = sea%area(i)
   elseif (op%stagpt == 'L') then
      fldval = land%area(i)
   endif

case(207) ! 'LS_TOPW'

   if (op%stagpt == 'S') then
      fldval = sea%topw(i)
   elseif (op%stagpt == 'L') then
      fldval = land%topw(i)
   endif

case(208) ! 'LS_ROUGH'

   if (op%stagpt == 'S') then
      fldval = sea%rough(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rough(i)
   endif

case(209) ! 'LS_VELS'

   if (op%stagpt == 'S') then
      fldval = sea%vels(i)
   elseif (op%stagpt == 'L') then
      fldval = land%vels(i)
   endif

case(210) ! 'LS_AIRTEMPK'

   if (op%stagpt == 'S') then
      fldval = sea%airtemp(i)
   elseif (op%stagpt == 'L') then
      fldval = land%airtemp(i)
   endif

case(211) ! 'LS_AIRSHV'

   if (op%stagpt == 'S') then
      fldval = sea%airshv(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%airshv(i) * 1.e3
   endif

case(212) ! 'LS_CANTEMPK'

   if (op%stagpt == 'S') then
      fldval = sea%cantemp(i)
   elseif (op%stagpt == 'L') then
      fldval = land%cantemp(i)
   endif

case(213) ! 'LS_CANSHV'

   if (op%stagpt == 'S') then
      fldval = sea%canshv(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%canshv(i) * 1.e3
   endif

case(214) ! 'LS_SKINTEMPK'

   if (op%stagpt == 'S') then
      nls = sea%nlev_seaice(i)
      
      if (nls > 0) then
         fldval = (1.0 - sea%seaicec(i)) * sea%seatc(i) &
                +        sea%seaicec(i)  * sea%seaice_tempk(nls,i)
      else
         fldval = sea%seatc(i)
      endif

   elseif (op%stagpt == 'L') then
      nls = land%nlev_sfcwater(i)

      if (nls > 0) then
         call qtk(land%sfcwater_energy(nls,i),tempk,fracliq)
      else
         call qwtk(land%soil_energy(nzg,i)       &
                  ,land%soil_water(nzg,i)*1.e3   &
                  ,slcpd(land%ntext_soil(nzg,i)) &
                  ,tempk, fracliq)
      endif

      fldval = (1. - land%vf(i)) * tempk &
                   + land%vf(i)  * land%veg_temp(i)
   endif

case(215) ! 'LS_GSS_SSH'

   if (op%stagpt == 'S') then
      fldval = sea%surface_ssh(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%surface_ssh(i) * 1.e3
   endif

case(216) ! 'LS_SENSFLUX'

   if (op%stagpt == 'S') then
      fldval = sea%sfluxt(i) * cp
   elseif (op%stagpt == 'L') then
      fldval = land%sfluxt(i) * cp
   endif

case(217) ! 'LS_LATFLUX'

   if (op%stagpt == 'S') then
      fldval = sea%sfluxr(i) * alvl
   elseif (op%stagpt == 'L') then
      fldval = land%sfluxr(i) * alvl
   endif

case(218) ! 'LS_VAPFLUX'

   if (op%stagpt == 'S') then
      fldval = sea%sfluxr(i)
   elseif (op%stagpt == 'L') then
      fldval = land%sfluxr(i)
   endif

case(219) ! 'LS_RSHORT'

   if (op%stagpt == 'S') then
      fldval = sea%rshort(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort(i)
   endif

case(220) ! 'LS_RLONG'

   if (op%stagpt == 'S') then
      fldval = sea%rlong(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong(i)
   endif

case(221) ! 'LS_RLONGUP'

   if (op%stagpt == 'S') then
      fldval = sea%rlongup(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlongup(i)
   endif

case(222) ! 'LS_RLONG_ALBEDO'

   if (op%stagpt == 'S') then
      fldval = sea%rlong_albedo(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong_albedo(i)
   endif

case(223) ! 'LS_ALBEDO_BEAM'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_beam(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_beam(i)
   endif

case(224) ! 'LS_ALBEDO_DIFFUSE'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_diffuse(i)
   endif

!-----------------------------------------
! LAND AND SEA CELLS - DIF2 fields
!-----------------------------------------

case(241) ! 'LS_VELS_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(vels_l_accum)) go to 1000
      fldval = (vels_l_accum_prev0(i) - vels_l_accum_prev1(i))
   elseif (op%stagpt == 'S') then
      if (.not. allocated(vels_s_accum)) go to 1000
      fldval = (vels_s_accum_prev0(i) - vels_s_accum_prev1(i))
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(242) ! 'LS_AIRTEMPK_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(airtemp_l_accum)) go to 1000
      fldval = (airtemp_l_accum_prev0(i) - airtemp_l_accum_prev1(i))
   elseif (op%stagpt == 'S') then
      if (.not. allocated(airtemp_s_accum)) go to 1000
      fldval = (airtemp_s_accum_prev0(i) - airtemp_s_accum_prev1(i))
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(243) ! 'LS_AIRSHV_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(airshv_l_accum)) go to 1000
      fldval = (airshv_l_accum_prev0(i) - airshv_l_accum_prev1(i)) * 1.e3
   elseif (op%stagpt == 'S') then
      if (.not. allocated(airshv_s_accum)) go to 1000
      fldval = (airshv_s_accum_prev0(i) - airshv_s_accum_prev1(i)) * 1.e3
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(244) ! 'LS_CANTEMPK_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantemp_l_accum)) go to 1000
      fldval = (cantemp_l_accum_prev0(i) - cantemp_l_accum_prev1(i))
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantemp_s_accum)) go to 1000
      fldval = (cantemp_s_accum_prev0(i) - cantemp_s_accum_prev1(i))
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(245) ! 'LS_CANSHV_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(canshv_l_accum)) go to 1000
      fldval = (canshv_l_accum_prev0(i) - canshv_l_accum_prev1(i)) * 1.e3
   elseif (op%stagpt == 'S') then
      if (.not. allocated(canshv_s_accum)) go to 1000
      fldval = (canshv_s_accum_prev0(i) - canshv_s_accum_prev1(i)) * 1.e3
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(246) ! 'LS_SKINTEMPK_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(skintemp_l_accum)) go to 1000
      fldval = (skintemp_l_accum_prev0(i) - skintemp_l_accum_prev1(i))
   elseif (op%stagpt == 'S') then
      if (.not. allocated(skintemp_s_accum)) go to 1000
      fldval = (skintemp_s_accum_prev0(i) - skintemp_s_accum_prev1(i))
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(247) ! 'LS_SENSFLUX_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxt_l_accum)) go to 1000
      fldval = (sfluxt_l_accum_prev0(i) - sfluxt_l_accum_prev1(i)) * cp
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxt_s_accum)) go to 1000
      fldval = (sfluxt_s_accum_prev0(i) - sfluxt_s_accum_prev1(i)) * cp
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(248) ! 'LS_LATFLUX_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxr_l_accum)) go to 1000
      fldval = (sfluxr_l_accum_prev0(i) - sfluxr_l_accum_prev1(i)) * alvl
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxr_s_accum)) go to 1000
      fldval = (sfluxr_s_accum_prev0(i) - sfluxr_s_accum_prev1(i)) * alvl
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(249) ! 'LS_VAPFLUX_DIF2'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxr_l_accum)) go to 1000
      fldval = (sfluxr_l_accum_prev0(i) - sfluxr_l_accum_prev1(i))
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxr_s_accum)) go to 1000
      fldval = (sfluxr_s_accum_prev0(i) - sfluxr_s_accum_prev1(i))
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

!-----------------------------------------
! LAND AND SEA CELLS - ATM averages
!-----------------------------------------

case(261:269) ! 'ALS_VELS',       'ALS_AIRTEMPK', 'ALS_AIRSHV',
              ! 'ALS_CANTEMPK',   'ALS_CANSHV',
              ! 'ALS_SENSFLUX',   'ALS_LATFLUX',
              ! 'ALS_VAPFLUX',    'ALS_SKINTEMPK'

   fldval = 0.
   area_sum = 0.

   do jland = 1,itab_w(i)%nland
      iland = itab_w(i)%iland(jland)

      if     (trim(fldname) == 'ALS_VELS'    ) then
         fldval = fldval + land%vels   (iland) * land%area(iland)
      elseif (trim(fldname) == 'ALS_AIRTEMPK') then
         fldval = fldval + land%airtemp(iland) * land%area(iland) 
      elseif (trim(fldname) == 'ALS_AIRSHV'  ) then
         fldval = fldval + land%airshv (iland) * land%area(iland) * 1.e3
      elseif (trim(fldname) == 'ALS_CANTEMPK') then
         fldval = fldval + land%cantemp(iland) * land%area(iland)
      elseif (trim(fldname) == 'ALS_CANSHV'  ) then
         fldval = fldval + land%canshv (iland) * land%area(iland) * 1.e3
      elseif (trim(fldname) == 'ALS_SENSFLUX') then
         fldval = fldval + land%sfluxt (iland) * land%area(iland) * cp
      elseif (trim(fldname) == 'ALS_LATFLUX' ) then
         fldval = fldval + land%sfluxr (iland) * land%area(iland) * alvl
      elseif (trim(fldname) == 'ALS_VAPFLUX' ) then
         fldval = fldval + land%sfluxr (iland) * land%area(iland)
      elseif (trim(fldname) == 'ALS_SKINTEMPK') then
         nls = land%nlev_sfcwater(iland)

         if (nls > 0) then
            call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)
         else
            call qwtk(land%soil_energy(nzg,iland)       &
                     ,land%soil_water(nzg,iland)*1.e3   &
                     ,slcpd(land%ntext_soil(nzg,iland)) &
                     ,tempk, fracliq)
         endif

         fldval = fldval + ((1. - land%vf(iland)) * tempk &
                                + land%vf(iland)  * land%veg_temp(iland)) &
                * land%area(iland)
      else
         cycle
      endif

      area_sum = area_sum + land%area(iland)
   enddo

   do jsea = 1,itab_w(i)%nsea
      isea = itab_w(i)%isea(jsea)

      if     (trim(fldname) == 'ALS_VELS   '  ) then
         fldval = fldval + sea%vels   (isea) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_AIRTEMPK' ) then
         fldval = fldval + sea%airtemp(isea) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_AIRSHV'   ) then
         fldval = fldval + sea%airshv (isea) * sea%area(isea) * 1.e3
      elseif (trim(fldname) == 'ALS_CANTEMPK' ) then
         fldval = fldval + sea%cantemp(isea) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_CANSHV'   ) then
         fldval = fldval + sea%canshv (isea) * sea%area(isea) * 1.e3
      elseif (trim(fldname) == 'ALS_SENSFLUX' ) then
         fldval = fldval + sea%sfluxt (isea) * sea%area(isea) * cp
      elseif (trim(fldname) == 'ALS_LATFLUX'  ) then
         fldval = fldval + sea%sfluxr (isea) * sea%area(isea) * alvl
      elseif (trim(fldname) == 'ALS_VAPFLUX'  ) then
         fldval = fldval + sea%sfluxr (isea) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_SKINTEMPK') then
         nls = sea%nlev_seaice(isea)

         if (nls > 0) then
            fldval = fldval &
                   + ((1.0 - sea%seaicec(isea)) * sea%seatc(isea) &
                           + sea%seaicec(isea)  * sea%seaice_tempk(nls,isea)) &
                   * sea%area(isea)
         else
            fldval = fldval + sea%seatc(isea) * sea%area(isea)
         endif

      else
         cycle
      endif

      area_sum = area_sum + sea%area(isea)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

!-------------------------------------------------
! LAND AND SEA CELLS - ATM averages of DIF2 fields
!-------------------------------------------------

case(281:289) ! 'ALS_VELS_DIF2',       'ALS_AIRTEMPK_DIF2', 'ALS_AIRSHV_DIF2',
              ! 'ALS_CANTEMPK_DIF2',   'ALS_CANSHV_DIF2',
              ! 'ALS_SENSFLUX_DIF2',   'ALS_LATFLUX_DIF2',
              ! 'ALS_VAPFLUX_DIF2',    'ALS_SKINTEMPK_DIF2'

   fldval = 0.
   area_sum = 0.

   do jland = 1,itab_w(i)%nland
      iland = itab_w(i)%iland(jland)

      if     (trim(fldname) == 'ALS_VELS_DIF2'    ) then
         if (.not. allocated(vels_l_accum)) go to 1000
         fldval = fldval + (vels_l_accum_prev0(iland) - vels_l_accum_prev1(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_AIRTEMPK_DIF2') then
         if (.not. allocated(airtemp_l_accum)) go to 1000
         fldval = fldval + (airtemp_l_accum_prev0(iland) - airtemp_l_accum_prev1(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_AIRSHV_DIF2'  ) then
         if (.not. allocated(airshv_l_accum)) go to 1000
         fldval = fldval + (airshv_l_accum_prev0(iland) - airshv_l_accum_prev1(iland)) * land%area(iland) * 1.e3
      elseif (trim(fldname) == 'ALS_CANTEMPK_DIF2') then
         if (.not. allocated(cantemp_l_accum)) go to 1000
         fldval = fldval + (cantemp_l_accum_prev0(iland) - cantemp_l_accum_prev1(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_CANSHV_DIF2'  ) then
         if (.not. allocated(canshv_l_accum)) go to 1000
         fldval = fldval + (canshv_l_accum_prev0(iland) - canshv_l_accum_prev1(iland)) * land%area(iland) * 1.e3
      elseif (trim(fldname) == 'ALS_SKINTEMPK_DIF2') then
         if (.not. allocated(skintemp_l_accum)) go to 1000
         fldval = fldval + (skintemp_l_accum_prev0(iland) - skintemp_l_accum_prev1(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_SENSFLUX_DIF2') then
         if (.not. allocated(sfluxt_l_accum)) go to 1000
         fldval = fldval + (sfluxt_l_accum_prev0(iland) - sfluxt_l_accum_prev1(iland)) * land%area(iland) * cp
      elseif (trim(fldname) == 'ALS_LATFLUX_DIF2' ) then
         if (.not. allocated(sfluxr_l_accum)) go to 1000
         fldval = fldval + (sfluxr_l_accum_prev0(iland) - sfluxr_l_accum_prev1(iland)) * land%area(iland) * alvl
      elseif (trim(fldname) == 'ALS_VAPFLUX_DIF2' ) then
         if (.not. allocated(sfluxr_l_accum)) go to 1000
         fldval = fldval + (sfluxr_l_accum_prev0(iland) - sfluxr_l_accum_prev1(iland)) * land%area(iland)
      else
         cycle
      endif

      area_sum = area_sum + land%area(iland)
   enddo

   do jsea = 1,itab_w(i)%nsea
      isea = itab_w(i)%isea(jsea)

      if     (trim(fldname) == 'ALS_VELS_DIF2'     ) then
         if (.not. allocated(vels_s_accum)) go to 1000
         fldval = fldval + (vels_s_accum_prev0(isea) - vels_s_accum_prev1(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_AIRTEMPK_DIF2' ) then
         if (.not. allocated(airtemp_s_accum)) go to 1000
         fldval = fldval + (airtemp_s_accum_prev0(isea) - airtemp_s_accum_prev1(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_AIRSHV_DIF2'   ) then
         if (.not. allocated(airshv_s_accum)) go to 1000
         fldval = fldval + (airshv_s_accum_prev0(isea) - airshv_s_accum_prev1(isea)) * sea%area(isea) * 1.e3
      elseif (trim(fldname) == 'ALS_CANTEMPK_DIF2' ) then
         if (.not. allocated(cantemp_s_accum)) go to 1000
         fldval = fldval + (cantemp_s_accum_prev0(isea) - cantemp_s_accum_prev1(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_CANSHV_DIF2'   ) then
         if (.not. allocated(canshv_s_accum)) go to 1000
         fldval = fldval + (canshv_s_accum_prev0(isea) - canshv_s_accum_prev1(isea)) * sea%area(isea) * 1.e3
      elseif (trim(fldname) == 'ALS_SKINTEMPK_DIF2') then
         if (.not. allocated(skintemp_s_accum)) go to 1000
         fldval = fldval + (skintemp_s_accum_prev0(isea) - skintemp_s_accum_prev1(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_SENSFLUX_DIF2' ) then
         if (.not. allocated(sfluxt_s_accum)) go to 1000
         fldval = fldval + (sfluxt_s_accum_prev0(isea) - sfluxt_s_accum_prev1(isea)) * sea%area(isea) * cp
      elseif (trim(fldname) == 'ALS_LATFLUX_DIF2'  ) then
         if (.not. allocated(sfluxr_s_accum)) go to 1000
         fldval = fldval + (sfluxr_s_accum_prev0(isea) - sfluxr_s_accum_prev1(isea)) * sea%area(isea) * alvl
      elseif (trim(fldname) == 'ALS_VAPFLUX_DIF2'  ) then
         if (.not. allocated(sfluxr_s_accum)) go to 1000
         fldval = fldval + (sfluxr_s_accum_prev0(isea) - sfluxr_s_accum_prev1(isea)) * sea%area(isea)
      else
         cycle
      endif

      area_sum = area_sum + sea%area(isea)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

!-------------------------------------------------
! LAND AND SEA CELLS - ATM averages of DIF4 fields
!-------------------------------------------------

case(290:292) ! 'ALS_AIRTEMPK_DIF4', 'ALS_CANTEMPK_DIF4', 'ALS_SKINTEMPK_DIF4'

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do jland = 1,itab_w(i)%nland
      iland = itab_w(i)%iland(jland)

      if     (trim(fldname) == 'ALS_AIRTEMPK_DIF4') then
         if (.not. allocated(airtemp_l_accum)) go to 1000
         fldval1 = fldval1 + (airtemp_l_accum_prev0(iland) - airtemp_l_accum_prev1(iland)) * land%area(iland)
         fldval2 = fldval2 + (airtemp_l_accum_prev2(iland) - airtemp_l_accum_prev3(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_CANTEMPK_DIF4') then
         if (.not. allocated(cantemp_l_accum)) go to 1000
         fldval1 = fldval1 + (cantemp_l_accum_prev0(iland) - cantemp_l_accum_prev1(iland)) * land%area(iland)
         fldval2 = fldval2 + (cantemp_l_accum_prev2(iland) - cantemp_l_accum_prev3(iland)) * land%area(iland)
      elseif (trim(fldname) == 'ALS_SKINTEMPK_DIF4') then
         if (.not. allocated(skintemp_l_accum)) go to 1000
         fldval1 = fldval1 + (skintemp_l_accum_prev0(iland) - skintemp_l_accum_prev1(iland)) * land%area(iland)
         fldval2 = fldval2 + (skintemp_l_accum_prev2(iland) - skintemp_l_accum_prev3(iland)) * land%area(iland)
      else
         cycle
      endif

      area_sum = area_sum + land%area(iland)
   enddo

   do jsea = 1,itab_w(i)%nsea
      isea = itab_w(i)%isea(jsea)

      if     (trim(fldname) == 'ALS_AIRTEMPK_DIF4' ) then
         if (.not. allocated(airtemp_s_accum)) go to 1000
         fldval1 = fldval1 + (airtemp_s_accum_prev0(isea) - airtemp_s_accum_prev1(isea)) * sea%area(isea)
         fldval2 = fldval2 + (airtemp_s_accum_prev2(isea) - airtemp_s_accum_prev3(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_CANTEMPK_DIF4' ) then
         if (.not. allocated(cantemp_s_accum)) go to 1000
         fldval1 = fldval1 + (cantemp_s_accum_prev0(isea) - cantemp_s_accum_prev1(isea)) * sea%area(isea)
         fldval2 = fldval2 + (cantemp_s_accum_prev2(isea) - cantemp_s_accum_prev3(isea)) * sea%area(isea)
      elseif (trim(fldname) == 'ALS_SKINTEMPK_DIF4') then
         if (.not. allocated(skintemp_s_accum)) go to 1000
         fldval1 = fldval1 + (skintemp_s_accum_prev0(isea) - skintemp_s_accum_prev1(isea)) * sea%area(isea)
         fldval2 = fldval2 + (skintemp_s_accum_prev2(isea) - skintemp_s_accum_prev3(isea)) * sea%area(isea)
      else
         cycle
      endif

      area_sum = area_sum + sea%area(isea)
   enddo

   if (area_sum > 1.e-3) then
      fldval1 = fldval1 / area_sum
      fldval2 = fldval2 / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

! GRID GEOMETRY - 3D

case(301) ! 'ARV'
   if (.not. allocated(arv)) go to 1000
   fldval = arv(k,i)
case(302) ! 'ARW'
   fldval = arw(k,i)
case(303) ! 'VOLT'
   fldval = volt(k,i)

! GRID GEOMETRY - 2D

case(311) ! 'TOPM'
   fldval = topm(i)
case(312) ! 'TOPW'
   fldval = topw(i)
case(313) ! 'GLATM'
   fldval = glatm(i)
case(314) ! 'GLONM'
   fldval = glonm(i)
case(315) ! 'GLATV'
   if (.not. allocated(glatv)) go to 1000
   fldval = glatv(i)
case(316) ! 'GLONV'
   if (.not. allocated(glonv)) go to 1000
   fldval = glonv(i)
case(317) ! 'GLATW'
   fldval = glatw(i)
case(318) ! 'GLONW'
   fldval = glonw(i)
case(319) ! 'LPM'
   fldval = real(lpm(i))
case(320) ! 'LPV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(321) ! 'LCV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(322) ! 'LPW'
   fldval = real(lpw(i))
case(323) ! 'LSW'
   fldval = real(lsw(i))
case(324) ! 'XEM'
   fldval = xem(i)
case(325) ! 'YEM'
   fldval = yem(i)
case(326) ! 'ZEM'
   fldval = zem(i)
case(327) ! 'XEV'
   if (.not. allocated(xev)) go to 1000
   fldval = xev(i)
case(328) ! 'YEV'
   if (.not. allocated(yev)) go to 1000
   fldval = yev(i)
case(329) ! 'ZEV'
   if (.not. allocated(zev)) go to 1000
   fldval = zev(i)
case(330) ! 'XEW'
   fldval = xew(i)
case(331) ! 'YEW'
   fldval = yew(i)
case(332) ! 'ZEW'
   fldval = zew(i)
case(333) ! 'DNU'
   fldval = dnu(i)
case(334) ! 'DNV'
   fldval = dnv(i)
case(335) ! 'ARM0'
   fldval = arm0(i)
case(336) ! 'ARW0'
   fldval = arw0(i)

! ITAB_M MEMBERS

case(341) ! 'ITAB_M_NPOLY'
   fldval = itab_m(i)%npoly
case(342) ! 'ITAB_M_IMGLOBE'
   fldval = itab_m(i)%imglobe
case(343) ! 'ITAB_M_MRLMR'
   fldval = itab_m(i)%mrlm
case(344) ! 'ITAB_M_MRLM_OR'
   fldval = itab_m(i)%mrlm_orig
case(345) ! 'ITAB_M_MROW'
   fldval = itab_m(i)%mrow
case(346) ! 'ITAB_M_NGR'
   fldval = itab_m(i)%ngr
case(347) ! 'ITAB_M_IV'
   fldval = itab_m(i)%iv(indp)
case(348) ! 'ITAB_M_IW'
   fldval = itab_m(i)%iw(indp)

! ITAB_V MEMBERS

case(351) ! 'ITAB_V_IVP'
   fldval = itab_v(i)%ivp
case(352) ! 'ITAB_V_IRANK'
   fldval = itab_v(i)%irank
!   fldval = itabg_v(i)%irank
case(353) ! 'ITAB_V_IVGLOBE'
   fldval = itab_v(i)%ivglobe
case(354) ! 'ITAB_V_MRLV'
   fldval = itab_v(i)%mrlv
case(355) ! 'ITAB_V_IM'
   fldval = itab_v(i)%im(indp)
case(356) ! 'ITAB_V_IV'
   fldval = itab_v(i)%iv(indp)
case(357) ! 'ITAB_V_IW'
   fldval = itab_v(i)%iw(indp)
case(358) ! 'ITAB_V_FARW'
   fldval = itab_v(i)%farw(indp)

! ITAB_W MEMBERS

case(361) ! 'ITAB_W_NPOLY'
   fldval = itab_w(i)%npoly
case(362) ! 'ITAB_W_IWP'
   fldval = itab_w(i)%iwp
case(363) ! 'ITAB_W_IRANK'
   fldval = itab_w(i)%irank
!   fldval = itabg_w(i)%irank
case(364) ! 'ITAB_W_IWGLOBE'
   fldval = itab_w(i)%iwglobe
case(365) ! 'ITAB_W_MRLW'
   fldval = itab_w(i)%mrlw
case(366) ! 'ITAB_W_MRLW_OR'
   fldval = itab_w(i)%mrlw_orig
case(367) ! 'ITAB_W_NGR'
   fldval = itab_w(i)%ngr
case(368) ! 'ITAB_W_IM'
   fldval = itab_w(i)%im(indp)
case(369) ! 'ITAB_W_IV'
   fldval = itab_w(i)%iv(indp)
case(370) ! 'ITAB_W_IW'
   fldval = itab_w(i)%iw(indp)
case(371) ! 'ITAB_W_DIRV'
   fldval = itab_w(i)%dirv(indp)
case(372) ! 'ITAB_W_FARM'
   fldval = itab_w(i)%farm(indp)
case(373) ! 'ITAB_W_FARV'
   fldval = itab_w(i)%farv(indp)
case(374) ! 'ITAB_W_IWNUD'
   fldval = itab_w(i)%iwnud(indp)
case(375) ! 'ITAB_W_FNUD'
   fldval = itab_w(i)%fnud(indp)

! Time-averaged fields

case(401) ! 'PRESS_MAVG' [indp = 1:nz_avg is the selected vertical level]

   if (.not. allocated(press_mavg)) go to 1000

   fldval = press_mavg(indp,i) * .01

case(402) !  'RHO_MAVG'

   if (.not. allocated(rho_mavg)) go to 1000

   fldval = rho_mavg(indp,i)

case(403) !  'TEMPK_MAVG'

   if (.not. allocated(tempk_mavg)) go to 1000

   fldval = tempk_mavg(indp,i)

case(404) !  'SH_V_MAVG'

   if (.not. allocated(sh_v_mavg)) go to 1000

   fldval = sh_v_mavg(indp,i) * 1.e3

case(405) !  'SH_W_MAVG'

   if (.not. allocated(sh_w_mavg)) go to 1000

   fldval = sh_w_mavg(indp,i) * 1.e3

case(406) !  'WC_MAVG'

   if (.not. allocated(wc_mavg)) go to 1000

   fldval = wc_mavg(indp,i)

case(407:408) !  'ZONAL_WINDW_MAVG','MERID_WINDW_MAVG'

   if (.not. allocated(vxe_mavg)) go to 1000
   if (.not. allocated(vye_mavg)) go to 1000
   if (.not. allocated(vze_mavg)) go to 1000

   npoly = itab_w(i)%npoly
   rpolyi = 1. / real(npoly)

   vx = vxe_mavg(indp,i)
   vy = vye_mavg(indp,i)
   vz = vze_mavg(indp,i)

   if (mdomain < 2) then

      raxis = sqrt(xew(i) ** 2 + yew(i) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
         u = (vy * xew(i) - vx * yew(i)) / raxis
         v = vz * raxis / erad &
           - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 
      else
         u = 0.
         v = 0.
      endif

   else
      u = vx
      v = vy
   endif

   if (trim(fldname) == 'ZONAL_WINDW_MAVG') then
      fldval = u
   elseif (trim(fldname) == 'MERID_WINDW_MAVG') then
      fldval = v
   endif

case(409) ! 'RSHORT_MAVG'

   if (.not. allocated(rshort_mavg)) go to 1000

   fldval = rshort_mavg(i)

case(410) ! 'RSHORT_TOP_MAVG'

   if (.not. allocated(rshort_top_mavg)) go to 1000

   fldval = rshort_top_mavg(i)

case(411) ! 'RSHORTUP_MAVG'

   if (.not. allocated(rshortup_mavg)) go to 1000

   fldval = rshortup_mavg(i)

case(412) ! 'RSHORTUP_TOP_MAVG'

   if (.not. allocated(rshortup_top_mavg)) go to 1000

   fldval = rshortup_top_mavg(i)

case(413) ! 'RLONG_MAVG'

   if (.not. allocated(rlong_mavg)) go to 1000

   fldval = rlong_mavg(i)

case(414) ! 'RLONGUP_MAVG'

   if (.not. allocated(rlongup_mavg)) go to 1000

   fldval = rlongup_mavg(i)

case(415) ! 'RLONGUP_TOP_MAVG'

   if (.not. allocated(rlongup_top_mavg)) go to 1000

   fldval = rlongup_top_mavg(i)

case(416) ! 'LATFLUX_MAVG'

   if (.not. allocated(latflux_mavg)) go to 1000

   fldval = latflux_mavg(i)

case(417) ! 'SENSFLUX_MAVG'

   if (.not. allocated(sensflux_mavg)) go to 1000

   fldval = sensflux_mavg(i)

case(418) ! 'WINDSPEED_MAVG'

   if (.not. allocated(windspeed_mavg)) go to 1000

   fldval = windspeed_mavg(i)

case(419) ! 'ACCPMIC_MTOT'

   if (.not. allocated(accpmic_mtot)) go to 1000

   fldval = accpmic_mtot(i)

case(420) ! 'ACCPCON_MTOT'

   if (.not. allocated(accpcon_mtot)) go to 1000

   fldval = accpcon_mtot(i)

case(421) ! 'ACCPBOTH_MTOT'

   if (.not. allocated(accpmic_mtot)) go to 1000
   if (.not. allocated(accpcon_mtot)) go to 1000

   fldval = accpmic_mtot(i) + accpcon_mtot(i)

case(422) ! 'PRESS_DAVG'

   if (.not. allocated(press_davg)) go to 1000

   fldval = press_davg(i) * .01

case(423:424) ! 'ZONAL_WINDW_DAVG','MERID_WINDW_DAVG'

   if (.not. allocated(vxe_davg)) go to 1000
   if (.not. allocated(vye_davg)) go to 1000
   if (.not. allocated(vze_davg)) go to 1000

   npoly = itab_w(i)%npoly
   rpolyi = 1. / real(npoly)

   vx = vxe_davg(i)
   vy = vye_davg(i)
   vz = vze_davg(i)

   if (mdomain < 2) then

      raxis = sqrt(xew(i) ** 2 + yew(i) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
         u = (vy * xew(i) - vx * yew(i)) / raxis
         v = vz * raxis / erad &
           - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 
      else
         u = 0.
         v = 0.
      endif

   else
      u = vx
      v = vy
   endif

   if (trim(fldname) == 'ZONAL_WINDW_DAVG') then
      fldval = u
   elseif (trim(fldname) == 'MERID_WINDW_DAVG') then
      fldval = v
   endif

case(425) ! 'RSHORT_DAVG'

   if (.not. allocated(rshort_davg)) go to 1000

   fldval = rshort_davg(i)

case(426) ! 'TEMPK_DAVG'

   if (.not. allocated(tempk_davg)) go to 1000

   fldval = tempk_davg(i)

case(427) ! 'TEMPK_DMIN'

   if (.not. allocated(tempk_dmin)) go to 1000

   fldval = tempk_dmin(i)

case(428) ! 'TEMPK_DMAX'

   if (.not. allocated(tempk_dmax)) go to 1000

   fldval = tempk_dmax(i)

case(429) ! 'ACCPMIC_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000

   fldval = accpmic_dtot(i)

case(430) ! 'ACCPCON_DTOT'

   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpcon_dtot(i)

case(431) ! 'ACCPBOTH_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000
   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpmic_dtot(i) + accpcon_dtot(i)

case(432) ! 'PRESS_UL_DAVG'

   if (.not. allocated(press_ul_davg)) go to 1000

   fldval = press_ul_davg(i) * .01

case(433:434) ! 'ZONAL_WINDW_UL_DAVG','MERID_WINDW_UL_DAVG'

   if (.not. allocated(vxe_ul_davg)) go to 1000
   if (.not. allocated(vye_ul_davg)) go to 1000
   if (.not. allocated(vze_ul_davg)) go to 1000

   npoly = itab_w(i)%npoly
   rpolyi = 1. / real(npoly)

   vx = vxe_ul_davg(i)
   vy = vye_ul_davg(i)
   vz = vze_ul_davg(i)

   if (mdomain < 2) then

      raxis = sqrt(xew(i) ** 2 + yew(i) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
         u = (vy * xew(i) - vx * yew(i)) / raxis
         v = vz * raxis / erad &
           - (vx * xew(i) + vy * yew(i)) * zew(i) / (raxis * erad) 
      else
         u = 0.
         v = 0.
      endif

   else
      u = vx
      v = vy
   endif

   if (trim(fldname) == 'ZONAL_WINDW_UL_DAVG') then
      fldval = u
   elseif (trim(fldname) == 'MERID_WINDW_UL_DAVG') then
      fldval = v
   endif

case(435) ! 'CANTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_davg)) go to 1000
      fldval = cantempk_l_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_davg)) go to 1000
      fldval = cantempk_s_davg(i)
   endif

case(436) ! 'CANTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_dmin)) go to 1000
      fldval = cantempk_l_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_dmin)) go to 1000
      fldval = cantempk_s_dmin(i)
   endif

case(437) ! 'CANTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_dmax)) go to 1000
      fldval = cantempk_l_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_dmax)) go to 1000
      fldval = cantempk_s_dmax(i)
   endif

case(438) ! 'VEGTEMPK_DAVG'

   if (.not. allocated(vegtempk_davg)) go to 1000

   fldval = vegtempk_davg(i)

case(439) ! 'VEGTEMPK_DMIN'

   if (.not. allocated(vegtempk_dmin)) go to 1000

   fldval = vegtempk_dmin(i)

case(440) ! 'VEGTEMPK_DMAX'

   if (.not. allocated(vegtempk_dmax)) go to 1000

   fldval = vegtempk_dmax(i)

case(441) ! 'SOILTEMPK_DAVG'

   if (.not. allocated(soiltempk_davg)) go to 1000

   fldval = soiltempk_davg(i)

case(442) ! 'SOILTEMPK_DMIN'

   if (.not. allocated(soiltempk_dmin)) go to 1000

   fldval = soiltempk_dmin(i)

case(443) ! 'SOILTEMPK_DMAX'

   if (.not. allocated(soiltempk_dmax)) go to 1000

   fldval = soiltempk_dmax(i)

case(444) ! 'LS_AIRTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(airtempk_l_davg)) go to 1000
      fldval = airtempk_l_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(airtempk_s_davg)) go to 1000
      fldval = airtempk_s_davg(i)
   endif

case(445) ! 'LS_AIRTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(airtempk_l_dmin)) go to 1000
      fldval = airtempk_l_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(airtempk_s_dmin)) go to 1000
      fldval = airtempk_s_dmin(i)
   endif

case(446) ! 'LS_AIRTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(airtempk_l_dmax)) go to 1000
      fldval = airtempk_l_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(airtempk_s_dmax)) go to 1000
      fldval = airtempk_s_dmax(i)
   endif

case(447) ! 'LS_CANTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_davg)) go to 1000
      fldval = cantempk_l_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_davg)) go to 1000
      fldval = cantempk_s_davg(i)
   endif

case(448) ! 'LS_CANTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_dmin)) go to 1000
      fldval = cantempk_l_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_dmin)) go to 1000
      fldval = cantempk_s_dmin(i)
   endif

case(449) ! 'LS_CANTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_l_dmax)) go to 1000
      fldval = cantempk_l_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_s_dmax)) go to 1000
      fldval = cantempk_s_dmax(i)
   endif

case(450) ! 'LS_SENSFLUX_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxt_l_davg)) go to 1000
      fldval = sfluxt_l_davg(i) * cp
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxt_s_davg)) go to 1000
      fldval = sfluxt_s_davg(i) * cp
   endif

case(451) ! 'LS_LATFLUX_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxr_l_davg)) go to 1000
      fldval = sfluxr_l_davg(i) * alvl
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxr_s_davg)) go to 1000
      fldval = sfluxr_s_davg(i) * alvl
   endif

case(452) ! 'SENSFLUX_DAVG'

   if (infotyp == 'UNITS') then
      aux(:) = 0.
      do iws = 2,mws
         iw = itab_ws(iws)%iw
         aux(iw) = aux(iw) + itab_ws(iws)%arf_iw * sfluxt_s_davg(isf) * cp
      enddo
      do iwl = 2,mwl
         iw = itab_wl(iwl)%iw
         aux(iw) = aux(iw) + itab_wl(iwl)%arf_iw * sfluxt_l_davg(ilf) * cp
      enddo
   endif

   fldval = aux(i)

case(453) ! 'LATFLUX_DAVG'

   if (infotyp == 'UNITS') then
      aux(:) = 0.
      do iws = 2,mws
         iw = itab_ws(iws)%iw
         aux(iw) = aux(iw) + itab_ws(iws)%arf_iw * sfluxr_s_davg(isf) * alvl
      enddo
      do iwl = 2,mwl
         iw = itab_wl(iwl)%iw
         aux(iw) = aux(iw) + itab_wl(iwl)%arf_iw * sfluxr_l_davg(ilf) * alvl
      enddo
   endif

   fldval = aux(i)

case(454) ! 'PRESS_AVG24' [indp = 1:24 is the hour of day to be plotted]

   if (.not. allocated(press_avg24)) go to 1000

   fldval = press_avg24(indp,i)

! Miscellaneous and new additions

case(471) ! 'RHO_OBS'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1))

case(472) ! 'THETA_OBS'

   if (.not. allocated(theta_obs)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i)

case(473) ! 'SHW_OBS'

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3

case(474) ! 'UZONAL_OBS'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1))

case(475) ! 'UMERID_OBS'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1))

case(476) ! 'RHO_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_sim(k,itab_w(i)%iwnud(1))

case(477) ! 'THETA_SIM'

   if (.not. allocated(theta_sim)) go to 1000

!   fldval = theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_sim(k,i)

case(478) ! 'SHW_SIM'

   if (.not. allocated(shw_sim)) go to 1000

   fldval = shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(479) ! 'UZONAL_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_sim(k,itab_w(i)%iwnud(1))

case(480) ! 'UMERID_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_sim(k,itab_w(i)%iwnud(1))

case(481) ! 'RHO_OBS_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1)) &
          - rho_sim(k,itab_w(i)%iwnud(1))

case(482) ! 'THETA_OBS_SIM'

   if (.not. allocated(theta_sim)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1)) &
!          - theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i) &
          - theta_sim(k,i)

case(483) ! 'SHW_OBS_SIM'

   if (.not. allocated(shw_sim)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3 &
          - shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(484) ! 'UZONAL_OBS_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1)) &
          - uzonal_sim(k,itab_w(i)%iwnud(1))

case(485) ! 'UMERID_OBS_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1)) &
          - umerid_sim(k,itab_w(i)%iwnud(1))

case(486) ! 'VXE'

   if (.not. allocated(vxe)) go to 1000

   fldval = wtbot * vxe(k ,i) &
          + wttop * vxe(kp,i)

case(487) ! 'VYE'

   if (.not. allocated(vye)) go to 1000

   fldval = wtbot * vye(k ,i) &
          + wttop * vye(kp,i)

case(488) ! 'VZE'

   if (.not. allocated(vze)) go to 1000

   fldval = wtbot * vze(k ,i) &
          + wttop * vze(kp,i)

case(489) ! 'PBLH'

   if (.not. allocated(pblh)) go to 1000

   fldval = pblh(i)

case(490) ! 'HKH'

   if (.not. allocated(hkh)) go to 1000

   fldval = wtbot * hkh(k ,i) / rho(k ,i) &
          + wttop * hkh(kp,i) / rho(kp,i)

case(491) ! 'SHW_HCONV'

   fldval = 0.

   npoly = itab_w(i)%npoly

   do klev = lpw(i),mza
      do j = 1,npoly

         iv = itab_w(i)%iv(j)
         iw = itab_w(i)%iw(j)

         fldval = fldval &
                + vmc(klev,iv) * itab_w(i)%dirv(j) * arv(klev,iv) &
                * 0.5 * (sh_w(klev,i) + sh_w(klev,iw))
      enddo
   enddo

   fldval = fldval / arw0(i)

case(492) ! 'SHV_HCONV'

   fldval = 0.

   npoly = itab_w(i)%npoly

   do klev = lpw(i),mza
      do j = 1,npoly

         iv = itab_w(i)%iv(j)
         iw = itab_w(i)%iw(j)

         fldval = fldval &
                + vmc(klev,iv) * itab_w(i)%dirv(j) * arv(klev,iv) &
                * 0.5 * (sh_v(klev,i) + sh_v(klev,iw))
      enddo
   enddo

   fldval = fldval / arw0(i)

case(493) ! 'CLDNUM'

   fldval = cldnum(i) * 1.e-6

case default

! Optional plot of external field

   if (allocated(op%extfld)) then
      fldval = op%extfld(k,i)
   else
      go to 1000
   endif

end select

if (infotyp == 'VALUV') then
   if (op%fldvalv_min > fldval) op%fldvalv_min = fldval
   if (op%fldvalv_max < fldval) op%fldvalv_max = fldval
endif

! normal RETURN

return

! RETURN for "Plot field UNAVAILABLE"

1000 continue

notavail = 3

return
end subroutine oplot_lib
