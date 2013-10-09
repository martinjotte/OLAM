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
subroutine oplot_lib(k,i,infotyp,fldname0,wtbot,wttop,fldval,notavail)

use mem_ijtabs,  only: itab_w, itab_u, itab_v, itab_m, itabg_u, itabg_w
use mem_basic,   only: umc, vmc, wmc, ump, vmp, uc, vc, wc, rho, press, &
                       thil, theta, tair, sh_w, sh_v, vxe, vye, vze
use mem_cuparm,  only: conprr, aconpr

use mem_grid,    only: mza, mua, mva, mwa, lpm, lpu, lcu, lpv, lpw, lsw, &
                       zm, zt, dnu, dnv, arw0, arm0, aru, arv, arw, volt, &
                       volti, volui, volvi, volwi, xem, yem, zem, &
                       xeu, yeu, zeu, xev, yev, zev, xew, yew, zew, &
                       topm, topw, glatm, glonm, glatu, glonu, glatv, glonv, &
                       glatw, glonw, unx, uny, unz, vnx, vny, vnz, dnu, dnv

use mem_leaf,    only: land, itabg_wl

use mem_sea,     only: sea, itabg_ws

use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                       con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                       con_ccn, con_gccn, con_ifn, &
                       pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                       accpd, accpr, accpp, accps, accpa, accpg, accph
                      
use mem_radiate, only: fthrd_lw, fthrd_sw, rshort, rlong, rlongup, albedt, &
                       rshort_top, rshortup_top, rlongup_top
use mem_addsc,   only: addsc
use mem_tend,    only: umt, vmt, wmt
use mem_turb,    only: vkm, vkm_sfc, sflux_w, sflux_t, sflux_r, pblh, hkm
use mem_nudge,   only: rho_obs, theta_obs, shw_obs, uzonal_obs, umerid_obs, &
                       rho_sim, theta_sim, shw_sim, uzonal_sim, umerid_sim

use misc_coms,   only: io6, pr01d, dn01d, th01d, time8, iparallel, meshtype, &
                       naddsc, mdomain
use oplot_coms,  only: op
use consts_coms, only: p00i, rocp, erad, piu180, cp, alvl, grav, omega2
use leaf_coms,   only: slcpd, nzg, slmsts, slz, mwl, dt_leaf
use mem_sflux,   only: landflux, seaflux, mlandflux, mseaflux
use sea_coms,    only: mws
use mem_timeavg, only: rshort_avg, rshortup_avg, rlong_avg, rlongup_avg, &
                       rshort_top_avg, rshortup_top_avg, rlongup_top_avg, &
                       sflux_t_avg, sflux_r_avg

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
     canltempk_davg, canltempk_dmin, canltempk_dmax, &
      vegtempk_davg,  vegtempk_dmin,  vegtempk_dmax, &
     soiltempk_davg, soiltempk_dmin, soiltempk_dmax, &
     canstempk_davg, canstempk_dmin, canstempk_dmax, &
        tempk_lf_davg,    tempk_lf_dmin,    tempk_lf_dmax, &
     cantempk_lf_davg, cantempk_lf_dmin, cantempk_lf_dmax, &
       sfluxt_lf_davg,   sfluxr_lf_davg,                   &
        tempk_sf_davg,    tempk_sf_dmin,    tempk_sf_dmax, &
     cantempk_sf_davg, cantempk_sf_dmin, cantempk_sf_dmax, &
       sfluxt_sf_davg,   sfluxr_sf_davg

use oname_coms,  only: nl

use mem_swtc5_refsoln_cubic

use mem_plot,    only:   time8_prev1,   time8_prev2,   time8_prev3, &
                       accpmic_prev1, accpmic_prev2, accpmic_prev3, &
                       accpcon_prev1, accpcon_prev2, accpcon_prev3, &
                          press_init,      rho_init,    theta_init, &
                             uc_init,       vc_init,   addsc1_init

implicit none

integer, intent(in) :: k,i
character(len=*), intent(in) :: infotyp,fldname0
real, intent(out) :: fldval
real, intent(in) :: wtbot, wttop

integer, intent(out) :: notavail  ! 0 - variable is available
                                  ! 1 - variable is below ground
                                  ! 2 - variable is above model top
                                  ! 3 - variable is not available in this run

integer :: klev,nls,iv,iw,kw
real :: raxis,u,v,ucint,vcint,farv2,rpolyi
real :: vx, vy, vz, vxc, vyc, vzc
real :: tempk,fracliq
real :: contrib
real :: tuu1, tuu2, tuu3, tuu4
integer :: iw1,iw2,iwl,iws
integer :: iu, iu1, iu2, iu3, iu4
integer :: npoly,j
integer :: lenstr, ic, ifield
integer, save :: indp, icase

real :: ucc,vcc
real :: ucc_init, vcc_init, vx_init, vy_init, vz_init, u_init, v_init
real :: accpboth_prev1, accpboth_prev2, accpboth_prev3
real :: denom

integer :: isf, ilf
integer, save :: icall = 0
real, save, allocatable :: aux(:)

real :: zanal_swtc5, zanal0_swtc5

integer, parameter :: nfields = 365
character(len=40) :: fldlib(4,nfields)
character(len=40), save :: fldname

!  fldname     stagpt/dimens     field description & units                   field #
!-----------------------------------------------------------------------------------
! ATMOSPHERE - 3D

data fldlib(1:4,  1:37)/ &
 'UMC'           ,'V3' ,'U-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  1
 'VMC'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  2
 'WMC'           ,'W3' ,'W MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'            ,& !p  3
 'UMP'           ,'V3' ,'U-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  4
 'VMP'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  5
 'UC'            ,'V3' ,'U-NORMAL VELOCITY',' (m s:S2:-1  )'                ,& !p  6
 'VC'            ,'V3' ,'V-NORMAL VELOCITY',' (m s:S2:-1  )'                ,& !p  7
 'WC'            ,'W3' ,'W VELOCITY',' (m s:S2:-1  )'                       ,& !p  8
 'RHO'           ,'T3' ,'AIR DENSITY',' (kg m:S2:-3  )'                     ,& !p  9
 'PRESS'         ,'T3' ,'PRESSURE',' (hPa)'                                 ,& !  10
 'THIL'          ,'T3' ,'ICE-LIQUID THETA',' (K)'                           ,& !p 11
 'THETA'         ,'T3' ,'THETA',' (K)'                                      ,& !p 12
 'AIRTEMPK'      ,'T3' ,'AIR TEMP',' (K)'                                   ,& !p 13
 'AIRTEMPC'      ,'T3' ,'AIR TEMP',' (C)'                                   ,& !p 14
 'SH_W'          ,'T3' ,'TOTAL WATER SPEC DENSITY',' (g kg:S2:-1  )'        ,& !p 15
 'SH_V'          ,'T3' ,'WATER VAPOR SPEC DENSITY',' (g kg:S2:-1  )'        ,& !p 16
 'SH_C'          ,'T3' ,'CLOUDWATER SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 17
 'SH_D'          ,'T3' ,'DRIZZLE SPEC DENSITY',' (g kg:S2:-1  )'            ,& !p 18
 'SH_R'          ,'T3' ,'RAIN SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 19
 'SH_P'          ,'T3' ,'PRISTINE ICE SPEC DENSITY',' (g kg:S2:-1  )'       ,& !p 20
 'SH_S'          ,'T3' ,'SNOW SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 21
 'SH_A'          ,'T3' ,'AGGREGATES SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 22
 'SH_G'          ,'T3' ,'GRAUPEL SPEC DENSITY',' (g kg:S2:-1  )'            ,& !p 23
 'SH_H'          ,'T3' ,'HAIL SPEC DENSITY',' (g kg:S2:-1  )'               ,& !p 24
 'SH_CP'         ,'T3' ,'CLOUD + PRIST ICE SPEC DENSITY',' (g kg:S2:-1  )'  ,& !p 25
 'SH_TOTCOND'    ,'T3' ,'CONDENSATE SPEC DENSITY',' (g kg:S2:-1  )'         ,& !p 26
 'CON_C'         ,'T3' ,'CLOUD DROPLET NUMBER CONCEN',' (# mg:S2:-1  )'     ,& !p 27
 'CON_D'         ,'T3' ,'DRIZZLE NUMBER CONCEN',' (# mg:S2:-1  )'           ,& !p 28
 'CON_R'         ,'T3' ,'RAIN NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 29
 'CON_P'         ,'T3' ,'PRISTINE ICE NUMBER CONCEN',' (# kg:S2:-1  )'      ,& !p 30
 'CON_S'         ,'T3' ,'SNOW NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 31
 'CON_A'         ,'T3' ,'AGGREGATES NUMBER CONCEN',' (# kg:S2:-1  )'        ,& !p 32
 'CON_G'         ,'T3' ,'GRAUPEL NUMBER CONCEN',' (# kg:S2:-1  )'           ,& !p 33
 'CON_H'         ,'T3' ,'HAIL NUMBER CONCEN',' (# kg:S2:-1  )'              ,& !p 34
 'CON_CCN'       ,'T3' ,'CCN NUMBER CONCEN',' (# mg:S2:-1  )'               ,& !p 35
 'CON_GCCN'      ,'T3' ,'GCCN NUMBER CONCEN',' (# mg:S2:-1  )'              ,& !p 36
 'CON_IFN'       ,'T3' ,'IFN NUMBER CONCEN',' (# kg:S2:-1  )'                / !p 37

data fldlib(1:4, 38:67)/ &
 'VKM'           ,'W3' ,'VERT TURB MOMENTUM K',' (N s m:S2:-2  )'           ,& !p 38
 'FTHRD'         ,'T3' ,'RADIATIVE THETA TENDENCY',' (K s:S2:-1  )'         ,& !p 39
 'SPEEDV'        ,'V3' ,'WIND SPEED AT V',' (m s:S2:-1  )'                  ,& !p 40
 'AZIMV'         ,'V3' ,'WIND AZIMUTH AT V',' (deg)'                        ,& !p 41
 'ZONAL_WINDV'   ,'V3' ,'ZONAL WIND AT V',' (m s:S2:-1  )'                  ,& !p 42
 'MERID_WINDV'   ,'V3' ,'MERIDIONAL WIND AT V',' (m s:S2:-1  )'             ,& !p 43
 'ZONAL_WINDV_P' ,'V3' ,'ZONAL WIND PERT AT V',' (m s:S2:-1  )'             ,& !  44
 'MERID_WINDV_P' ,'V3' ,'MERIDIONAL WIND PERT V',' (m s:S2:-1  )'           ,& !  45
 'SPEEDW'        ,'T3' ,'WIND SPEED AT W',' (m s:S2:-1  )'                  ,& !p 46
 'AZIMW'         ,'T3' ,'WIND AZIMUTH AT W',' (deg)'                        ,& !p 47
 'ZONAL_WINDW'   ,'T3' ,'ZONAL WIND AT W',' (m s:S2:-1  )'                  ,& !p 48
 'MERID_WINDW'   ,'T3' ,'MERIDIONAL WIND AT W',' (m s:S2:-1  )'             ,& !p 49
 'RVORTZM'       ,'P3' ,'REL VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 50
 'TVORTZM'       ,'P3' ,'TOT VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 51
 'RVORTZM_P'     ,'P3' ,'REL VERT VORTICITY PERT AT M',' (s:S2:-1  )'       ,& !p 52
 'DIVERG'        ,'T3' ,'HORIZONTAL DIVERGENCE',' (s:S2:-1  )'              ,& !p 53
 'UMASSFLUX'     ,'V3' ,'GRID CELL U-FACE MASS FLUX',' (kg s:S2:-1  )'      ,& !  54
 'VMASSFLUX'     ,'V3' ,'GRID CELL V-FACE MASS FLUX',' (kg s:S2:-1  )'      ,& !  55
 'UC_P'          ,'V3' ,'NORMAL WIND PERT AT U',' (m s:S2:-1  )'            ,& !p 56
 'VC_P'          ,'V3' ,'NORMAL WIND PERT AT V',' (m s:S2:-1  )'            ,& !p 57
 'PRESS_P'       ,'T3' ,'PRESSURE PERT',' (hPa)'                            ,& !  58
 'RHO_P'         ,'T3' ,'DENSITY PERT',' (kg m:S2:-3  )'                    ,& !  59
 'THETA_P'       ,'T3' ,'THETA PERT',' (K)'                                 ,& !  60
 'AIRTEMPK_P'    ,'T3' ,'AIR TEMP PERT',' (K)'                              ,& !p 61
 'UMT'           ,'V3' ,'U-NORM MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'  ,& !  62
 'VMT'           ,'V3' ,'V-NORM MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'  ,& !  63
 'WMT'           ,'W3' ,'W MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'       ,& !  64
 'ADDSC'         ,'T3' ,'ADDED SCALAR AMOUNT PER KG AIR',' '                ,& !p 65
 'ADDSCP'        ,'T3' ,'SCALAR PERTURBATION',' ( )'                        ,& !  66
 'ZPLEV'         ,'T3' ,'HEIGHT OF CONST P SFC',' (m)'                       / !p 67

! ATMOSPHERE - 2D

data fldlib(1:4, 68:109)/ &
 'RSHORT_TOP'    ,'T2' ,'TOP DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'         ,& !  68
 'RSHORTUP_TOP'  ,'T2' ,'TOP UPWARD SHORTWV FLX',' (W m:S2:-2  )'           ,& !  69
 'RLONGUP_TOP'   ,'T2' ,'TOP UPWARD LONGWV FLX',' (W m:S2:-2  )'            ,& !  70

! ATMOSPHERE SURFACE - 2D

 'RSHORT'        ,'T2' ,'SFC DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'         ,& !  71
 'RSHORTUP'      ,'T2' ,'SFC UPWARD SHORTWV FLX',' (W m:S2:-2  )'           ,& !  72
 'RLONG'         ,'T2' ,'SFC DOWNWARD LONGWV FLX',' (W m:S2:-2  )'          ,& !  73
 'RLONGUP'       ,'T2' ,'SFC UPWARD LONGWV FLX',' (W m:S2:-2  )'            ,& !  74
 'ALBEDT'        ,'T2' ,'NET GRID COLUMN SFC ALBEDO',' ( )'                 ,& !  75
 'VKM_SFC'       ,'T2' ,'SFC TURB K FOR MOMENTUM',' (N s m:S2:-2  )'        ,& !  76
 'SFLUX_W'       ,'T2' ,'SFC W MOMENTUM FLUX',' (N m:S2:-2  )'              ,& !  77
 'SENSFLUX'      ,'T2' ,'ATM SFC SENSIBLE HEAT FLUX',' (W m:S2:-2  )'       ,& !  78
 'VAPFLUX'       ,'T2' ,'ATM SFC VAPOR FLUX',' (kg m:S2:-2   s:S2:-1  )'    ,& !  79
 'LATFLUX'       ,'T2' ,'ATM SFC LATENT HEAT FLUX',' (W m:S2:-2  )'         ,& !  80
 'PCPRD'         ,'T2' ,'DRIZZLE PRECIP RATE',' (kg m:S2:-2   h:S2:-1  )'   ,& !  81
 'PCPRR'         ,'T2' ,'RAIN PRECIP RATE',' (kg m:S2:-2   h:S2:-1  )'      ,& !  82
 'PCPRP'         ,'T2' ,'PRIST ICE PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    ,& !  83
 'PCPRS'         ,'T2' ,'SNOW PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  84
 'PCPRA'         ,'T2' ,'AGGREGATES PCP RATE',' (kg m:S2:-2   h:S2:-1  )'   ,& !  85
 'PCPRG'         ,'T2' ,'GRAUPEL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'      ,& !  86
 'PCPRH'         ,'T2' ,'HAIL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  87
 'PCPRMIC'       ,'T2' ,'MICROPHYS PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    ,& !  88
 'PCPRCON'       ,'T2' ,'CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  89
 'PCPRBOTH'      ,'T2' ,'MICRO + CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )' ,& !  90
 'ACCPD'         ,'T2' ,'ACCUM DRIZZLE',' (kg m:S2:-2  )'                   ,& !  91
 'ACCPR'         ,'T2' ,'ACCUM RAIN',' (kg m:S2:-2  )'                      ,& !  92
 'ACCPP'         ,'T2' ,'ACCUM PRIST ICE',' (kg m:S2:-2  )'                 ,& !  93
 'ACCPS'         ,'T2' ,'ACCUM SNOW',' (kg m:S2:-2  )'                      ,& !  94
 'ACCPA'         ,'T2' ,'ACCUM AGGREGATES',' (kg m:S2:-2  )'                ,& !  95
 'ACCPG'         ,'T2' ,'ACCUM GRAUPEL',' (kg m:S2:-2  )'                   ,& !  96
 'ACCPH'         ,'T2' ,'ACCUM HAIL',' (kg m:S2:-2  )'                      ,& !  97
 'ACCPMIC'       ,'T2' ,'ACCUM MICPHYS PCP',' (kg m:S2:-2  )'               ,& !  98
 'ACCPCON'       ,'T2' ,'ACCUM CONV PCP',' (kg m:S2:-2  )'                  ,& !  99
 'ACCPBOTH'      ,'T2' ,'ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )'        ,& ! 100
 'PCPDIF2MIC'    ,'T2' ,'MICPHYS PRECIP DIFF2',' (mm/day)'                  ,& ! 101
 'PCPDIF2CON'    ,'T2' ,'CONV PRECIP DIFF2',' (mm/day)'                     ,& ! 102
 'PCPDIF2BOTH'   ,'T2' ,'MICPHYS + CONV PRECIP DIFF2',' (mm/day)'           ,& ! 103
 'PCPDIF4MIC'    ,'T2' ,'MICPHYS PRECIP DIFF4',' (mm/day)'                  ,& ! 104
 'PCPDIF4CON'    ,'T2' ,'CONV PRECIP DIFF4',' (mm/day)'                     ,& ! 105
 'PCPDIF4BOTH'   ,'T2' ,'MICPHYS + CONV PRECIP DIFF4',' (mm/day)'           ,& ! 106
 'PCPREL4MIC'    ,'T2' ,'MICPHYS PRECIP RELATIVE DIFF4',' '                 ,& ! 107
 'PCPREL4CON'    ,'T2' ,'CONV PRECIP RELATIVE DIFF4',' '                    ,& ! 108
 'PCPREL4BOTH'   ,'T2' ,'MICPHYS + CONV PRECIP RELATIVE DIFF4',' '           / ! 109

! LAND_CELLS - 3D

data fldlib(1:4,110:132)/ &
 'SOIL_TEXT'     ,'L3G','SOIL TEXTURAL CLASS',' ( )'                        ,& ! 110
 'SOIL_ENERGY'   ,'L3G','SOIL ENERGY',' (J cm:S2:-3  )'                     ,& ! 111
 'SOIL_TEMPC'    ,'L3G','SOIL TEMP',' (C)'                                  ,& ! 112
 'SOIL_FRACLIQ'  ,'L3G','LIQUID FRACTION OF SOIL WATER',' ( )'              ,& ! 113
 'SOIL_WATER'    ,'L3G','SOIL WATER CONTENT',' ( )'                         ,& ! 114
 'SFWAT_MASS'    ,'L3S','SFCWATER MASS',' (kg m:S2:-2  )'                   ,& ! 115
 'SFWAT_ENERGY'  ,'L3S','SFCWATER ENERGY',' (J g:S2:-1  )'                  ,& ! 116
 'SFWAT_TEMPC'   ,'L3S','SFCWATER TEMP',' (C)'                              ,& ! 117
 'SFWAT_FRACLIQ' ,'L3S','SFCWATER LIQUID FRACTION',' ( )'                   ,& ! 118
 'SFWAT_DEPTH'   ,'L3S','SFCWATER DEPTH',' (m)'                             ,& ! 119

! LAND_CELLS - 2D

 'NLEV_SFWAT'    ,'L2' ,'NUMBER OF SFCWATER LAYERS',' ( )'                  ,& ! 120
 'VEG_NDVIC'     ,'L2' ,'VEGETATION NDVI',' ( )'                            ,& ! 121
 'VEG_TEMPC'     ,'L2' ,'VEGETATION TEMP',' (C)'                            ,& ! 122
 'VEG_WATER'     ,'L2' ,'VEGETATION SFC WATER ',' (kg m:S2:-2  )'           ,& ! 123
 'STOM_RESIST'   ,'L2' ,'STOMATAL RESISTANCE',' (s m:S2:-1  )'              ,& ! 124
 'GROUND_SHV'    ,'L2' ,'EQUIL VAP SPEC DENSITY OF SOIL',' (g kg:S2:-1  )'  ,& ! 125
 'SOIL_DEPTH'    ,'L2' ,'SOIL DEPTH',' (m)'                                 ,& ! 126

! SEA_CELLS - 2D

 'SEATP'         ,'S2' ,'SEA SFC TEMP (PAST DATA)',' (K)'                   ,& ! 127
 'SEATF'         ,'S2' ,'SEA SFC TEMP (FUTURE DATA)',' (K)'                 ,& ! 128
 'SEATC'         ,'S2' ,'SEA SFC TEMP (CURRENT)',' (K)'                     ,& ! 129
 'SEAICEP'       ,'S2' ,'SEAICE FRACTION (PAST DATA)',' ( )'                ,& ! 130
 'SEAICEF'       ,'S2' ,'SEAICE FRACTION (FUTURE DATA)',' ( )'              ,& ! 131
 'SEAICEC'       ,'S2' ,'SEAICE FRACTION (CURRENT)',' ( )'                   / ! 132

! LAND AND SEA CELLS - 2D

data fldlib(1:4,133:149)/ &
 'LEAF_CLASS'    ,'B2' ,'LEAF CLASS',' ( )'                                 ,& ! 133
 'AREA'          ,'B2' ,'LAND/SEA CELL AREA',' (m:S2:2  )'                  ,& ! 134
 'ROUGH'         ,'B2' ,'NET ROUGHNESS HEIGHT',' (m)'                       ,& ! 135
 'CAN_TEMPC'     ,'B2' ,'CANOPY AIR TEMP',' (C)'                            ,& ! 136
 'CAN_SHV'       ,'B2' ,'CANOPY VAPOR SPEC DENSITY',' (g kg:S2:-1  )'       ,& ! 137
 'SFC_TEMPC'     ,'B2' ,'SOIL/SFCWATER/SEA SFC TEMP',' (C)'                 ,& ! 138
 'SFC_SSH'       ,'B2' ,'LAND/SEA SFC SAT VAP SPEC DENS',' (g kg:S2:-1  )'  ,& ! 139

 'SENSFLUX_LS'      ,'B2','L/S CAN TOP SENS HEAT FLX',' (W m:S2:-2  )'      ,& ! 140
 'VAPFLUX_LS'       ,'B2','L/S CAN TOP VAP FLX',' (kg m:S2:-2   s:S2:-1  )' ,& ! 141
 'LATFLUX_LS'       ,'B2','L/S CAN TOP LAT HEAT FLX',' (W m:S2:-2  )'       ,& ! 142
 'RSHORT_LS'        ,'B2','L/S CAN TOP DOWN SW FLX',' (W m:S2:-2  )'        ,& ! 143
 'RSHORT_DIFFUSE_LS','B2','L/S CAN TOP DOWN DIFFUSE SW FLX',' (W m:S2:-2  )',& ! 144
 'RLONG_LS'         ,'B2','L/S CAN TOP DOWN LW FLX',' (W m:S2:-2  )'        ,& ! 145
 'RLONGUP_LS'       ,'B2','L/S CAN TOP UP LW FLX',' (W m:S2:-2  )'          ,& ! 146
 'RLONG_ALBEDO_LS'  ,'B2','L/S NET SFC LW ALBEDO',' ( )'                    ,& ! 147
 'ALBEDO_BEAM_LS'   ,'B2','L/S NET SFC BEAM ALBEDO',' ( )'                  ,& ! 148
 'ALBEDO_DIFFUSE_LS','B2','L/S NET SFC DIFFUSE ALBEDO',' ( )'                / ! 149

! FLUX CELLS - 2D

data fldlib(1:4,150:162)/  &
 'FCELL_ILSF'    ,'F2' ,'FLUXCELL INDEX',' '                                ,& ! 150
 'FCELL_IWLS'    ,'F2' ,'FLUXCELL LAND/SEA INDEX',' '                       ,& ! 151
 'FCELL_IW'      ,'F2' ,'FLUXCELL ATM IW INDEX',' '                         ,& ! 152
 'FCELL_KW'      ,'F2' ,'FLUXCELL ATM KW INDEX',' '                         ,& ! 153
 'FCELL_AREA'    ,'F2' ,'FLUXCELL AREA',' (m:S2:2  )'                       ,& ! 154
 'FCELL_ARFATM'  ,'F2' ,'FLUXCELL ATM AREA FRACTION',' '                    ,& ! 155
 'FCELL_ARFLS'   ,'F2' ,'FLUXCELL LAND/SEA AREA FRACTION',' '               ,& ! 156
 'FCELL_SENS'    ,'F2' ,'FLUXCELL SENS HEAT FLUX',' (W m:S2:-2  )'          ,& ! 157
 'FCELL_VAP'     ,'F2' ,'FLUXCELL VAP FLUX',' (kg m:S2:-2   s:S2:-1  )'     ,& ! 158
 'FCELL_LAT'     ,'F2' ,'FLUXCELL LAT HEAT FLUX',' (W m:S2:-2  )'           ,& ! 159
 'FCELL_AIRTEMPC','F2' ,'FLUXCELL AIR TEMP',' (C)'                          ,& ! 160
 'FCELL_AIRTEMPK','F2' ,'FLUXCELL AIR TEMP',' (K)'                          ,& ! 161
 'FCELL_CANTEMPC','F2' ,'FLUXCELL CAN TEMP',' (C)'                           / ! 162

! GRID GEOMETRY - 3D

data fldlib(1:4,163:169)/ &
 'ARU'           ,'V3' ,'AREA OF GRID CELL U-FACE',' (m:S2:2  )'            ,& ! 163
 'ARV'           ,'V3' ,'AREA OF GRID CELL V-FACE',' (m:S2:2  )'            ,& ! 164
 'ARW'           ,'W3' ,'AREA OF GRID CELL W-FACE',' (m:S2:2  )'            ,& ! 165
 'VOLT'          ,'T3' ,'GRID T-CELL VOLUME',' (m:S2:3  )'                  ,& ! 166
 'VOLU'          ,'V3' ,'GRID U-CELL VOLUME',' (m:S2:3  )'                  ,& ! 167
 'VOLV'          ,'V3' ,'GRID V-CELL VOLUME',' (m:S2:3  )'                  ,& ! 168
 'VOLW'          ,'W3' ,'GRID W-CELL VOLUME',' (m:S2:3  )'                   / ! 169

! GRID GEOMETRY - 2D

data fldlib(1:4,170:202)/ &
 'TOPM'          ,'M2' ,'TOPOGRAPHY HEIGHT',' (m)'                          ,& ! 170
 'TOPW'          ,'W2' ,'TOPOGRAPHY HEIGHT AT W',' (m)'                     ,& ! 171
 'GLATM'         ,'M2' ,'LATITUDE AT M',' (deg)'                            ,& ! 172
 'GLONM'         ,'M2' ,'LONGITUDE AT M',' (deg)'                           ,& ! 173
 'GLATU'         ,'V2' ,'LATITUDE AT U',' (deg)'                            ,& ! 174
 'GLONU'         ,'V2' ,'LONGITUDE AT U',' (deg)'                           ,& ! 175
 'GLATV'         ,'V2' ,'LATITUDE AT V',' (deg)'                            ,& ! 176
 'GLONV'         ,'V2' ,'LONGITUDE AT V',' (deg)'                           ,& ! 177
 'GLATW'         ,'T2' ,'LATITUDE',' (deg)'                                 ,& ! 178
 'GLONW'         ,'T2' ,'LONGITUDE',' (deg)'                                ,& ! 179
 'LPM'           ,'M2' ,'LOWEST PREDICTED M LEVEL',' ( )'                   ,& ! 180
 'LPU'           ,'V2' ,'LOWEST PREDICTED U LEVEL',' ( )'                   ,& ! 181
 'LCU'           ,'V2' ,'LOWEST ACTIVE U CONTROL VOL',' ( )'                ,& ! 182
 'LPV'           ,'V2' ,'LOWEST PREDICTED V LEVEL',' ( )'                   ,& ! 183
 'LCV'           ,'V2' ,'LOWEST ACTIVE V CONTROL VOL',' ( )'                ,& ! 184
 'LPW'           ,'W2' ,'LOWEST PREDICTED W LEVEL',' ( )'                   ,& ! 185
 'LSW'           ,'W2' ,'NUMBER OF SFC W LEVELS',' ( )'                     ,& ! 186
 'XEM'           ,'M2' ,'EARTH-X COORD OF M POINT',' ( )'                   ,& ! 187
 'YEM'           ,'M2' ,'EARTH-Y COORD OF M POINT',' ( )'                   ,& ! 188
 'ZEM'           ,'M2' ,'EARTH-Z COORD OF M POINT',' ( )'                   ,& ! 189
 'XEU'           ,'V2' ,'EARTH-X COORD OF U POINT',' ( )'                   ,& ! 190
 'YEU'           ,'V2' ,'EARTH-Y COORD OF U POINT',' ( )'                   ,& ! 191
 'ZEU'           ,'V2' ,'EARTH-Z COORD OF U POINT',' ( )'                   ,& ! 192
 'XEV'           ,'V2' ,'EARTH-X COORD OF V POINT',' ( )'                   ,& ! 193
 'YEV'           ,'V2' ,'EARTH-Y COORD OF V POINT',' ( )'                   ,& ! 194
 'ZEV'           ,'V2' ,'EARTH-Z COORD OF V POINT',' ( )'                   ,& ! 195
 'XEW'           ,'W2' ,'EARTH-X COORD OF W POINT',' ( )'                   ,& ! 196
 'YEW'           ,'W2' ,'EARTH-Y COORD OF W POINT',' ( )'                   ,& ! 197
 'ZEW'           ,'W2' ,'EARTH-Z COORD OF W POINT',' ( )'                   ,& ! 198
 'DNU'           ,'V2' ,'DNU',' (m)'                                        ,& ! 199
 'DNV'           ,'V2' ,'DNV',' (m)'                                        ,& ! 200
 'ARM0'          ,'W2' ,'SFC AREA OF M CELL',' (m:S2:2  )'                  ,& ! 201
 'ARW0'          ,'W2' ,'SFC AREA OF W CELL',' (m:S2:2  )'                   / ! 202

! ITAB_M MEMBERS - 2D

data fldlib(1:4,203:238)/ &
 'ITAB_M_NPOLY'  ,'M2' ,'ITAB_M_NPOLY',' ( )'                               ,& ! 203
 'ITAB_M_IMGLOBE','M2' ,'ITAB_M_IMGLOBE',' ( )'                             ,& ! 204
 'ITAB_M_MRLM'   ,'V2' ,'ITAB_M_MRLM',' ( )'                                ,& ! 205
 'ITAB_M_MRLM_OR','V2' ,'ITAB_M_MRLM_ORIG',' ( )'                           ,& ! 206
 'ITAB_M_MROW'   ,'M2' ,'ITAB_M_MROW',' ( )'                                ,& ! 207
 'ITAB_M_MROWH'  ,'M2' ,'ITAB_M_MROWH',' ( )'                               ,& ! 208
 'ITAB_M_IU'     ,'M2' ,'ITAB_M_IU',' ( )'                                  ,& ! 209
 'ITAB_M_IV'     ,'M2' ,'ITAB_M_IV',' ( )'                                  ,& ! 210
 'ITAB_M_IW'     ,'M2' ,'ITAB_M_IW',' ( )'                                  ,& ! 211
 'ITAB_M_FMW'    ,'M2' ,'ITAB_M_FMW',' ( )'                                 ,& ! 212

! ITAB_U MEMBERS - 2D

 'ITAB_U_IUP'    ,'V2' ,'ITAB_U_IUP',' ( )'                                 ,& ! 213
 'ITAB_U_IRANK'  ,'V2' ,'ITAB_U_IRANK',' ( )'                               ,& ! 214
 'ITAB_U_IUGLOBE','V2' ,'ITAB_U_IUGLOBE',' ( )'                             ,& ! 215
 'ITAB_U_MRLU'   ,'V2' ,'ITAB_U_MRLU',' ( )'                                ,& ! 216
 'ITAB_U_IM'     ,'V2' ,'ITAB_U_IM',' ( )'                                  ,& ! 217
 'ITAB_U_IU'     ,'V2' ,'ITAB_U_IU',' ( )'                                  ,& ! 218
 'ITAB_U_IW'     ,'V2' ,'ITAB_U_IW',' ( )'                                  ,& ! 219
 'ITAB_U_DIRU'   ,'V2' ,'ITAB_U_DIRU',' ( )'                                ,& ! 220
 'ITAB_U_FUU'    ,'V2' ,'ITAB_U_FUU',' ( )'                                 ,& ! 221
 'ITAB_U_FUW'    ,'V2' ,'ITAB_U_FUW',' ( )'                                 ,& ! 222
 'ITAB_U_TUU'    ,'V2' ,'ITAB_U_TUU',' ( )'                                 ,& ! 223
 'ITAB_U_GUW'    ,'V2' ,'ITAB_U_GUW',' ( )'                                 ,& ! 224
 'ITAB_U_PGC12'  ,'V2' ,'ITAB_U_PGC12',' (m:S2:-1  )'                       ,& ! 225
 'ITAB_U_PGC12b' ,'V2' ,'ITAB_U_PGC12b',' (m:S2:-1  )'                      ,& ! 226
 'ITAB_U_PGC12c' ,'V2' ,'ITAB_U_PGC12c',' (m:S2:-1  )'                      ,& ! 227
 'ITAB_U_PGC12d' ,'V2' ,'ITAB_U_PGC12d',' (m:S2:-1  )'                      ,& ! 228
 'ITAB_U_PGC45'  ,'V2' ,'ITAB_U_PGC45',' (m:S2:-1  )'                       ,& ! 229
 'ITAB_U_PGC45b' ,'V2' ,'ITAB_U_PGC45b',' (m:S2:-1  )'                      ,& ! 230
 'ITAB_U_PGC63'  ,'V2' ,'ITAB_U_PGC63',' (m:S2:-1  )'                       ,& ! 231
 'ITAB_U_PGC63c' ,'V2' ,'ITAB_U_PGC63c',' (m:S2:-1  )'                      ,& ! 232
 'ITAB_U_VXU_U'  ,'V2' ,'ITAB_U_VXU_U',' ( )'                               ,& ! 233
 'ITAB_U_VYU_U'  ,'V2' ,'ITAB_U_VYU_U',' ( )'                               ,& ! 234
 'ITAB_U_VXW_U'  ,'V2' ,'ITAB_U_VXW_U',' ( )'                               ,& ! 235
 'ITAB_U_VYW_U'  ,'V2' ,'ITAB_U_VYW_U',' ( )'                               ,& ! 236
 'ITAB_U_CROSSMM','V2' ,'ITAB_U_CROSSMM',' ( )'                             ,& ! 237
 'ITAB_U_CROSSWW','V2' ,'ITAB_U_CROSSWW',' ( )'                              / ! 238

! ITAB_V MEMBERS - 2D

data fldlib(1:4,239:249)/  &
 'ITAB_V_IVP'    ,'V2' ,'ITAB_V_IVP',' ( )'                                 ,& ! 239
 'ITAB_V_IRANK'  ,'V2' ,'ITAB_V_IRANK',' ( )'                               ,& ! 240
 'ITAB_V_IVGLOBE','V2' ,'ITAB_V_IVGLOBE',' ( )'                             ,& ! 241
 'ITAB_V_MRLV'   ,'V2' ,'ITAB_V_MRLV',' ( )'                                ,& ! 242
 'ITAB_V_IM'     ,'V2' ,'ITAB_V_IM',' ( )'                                  ,& ! 243
 'ITAB_V_IV'     ,'V2' ,'ITAB_V_IV',' ( )'                                  ,& ! 244
 'ITAB_V_IW'     ,'V2' ,'ITAB_V_IW',' ( )'                                  ,& ! 245
 'ITAB_V_FVV'    ,'V2' ,'ITAB_V_FVV',' ( )'                                 ,& ! 246
 'ITAB_V_FVW'    ,'V2' ,'ITAB_V_FVW',' ( )'                                 ,& ! 247
 'ITAB_V_FUV'    ,'V2' ,'ITAB_V_FUV',' ( )'                                 ,& ! 248
 'ITAB_V_FARW'   ,'V2' ,'ITAB_V_FARW',' ( )'                                 / ! 249

! ITAB_W MEMBERS - 2D

data fldlib(1:4,250:277)/ &

 'ITAB_W_NPOLY'  ,'W2' ,'ITAB_W_NPOLY',' ( )'                               ,& ! 250
 'ITAB_W_IWP'    ,'W2' ,'ITAB_W_IWP',' ( )'                                 ,& ! 251
 'ITAB_W_IRANK'  ,'W2' ,'ITAB_W_IRANK',' ( )'                               ,& ! 252
 'ITAB_W_IWGLOBE','W2' ,'ITAB_W_IWGLOBE',' ( )'                             ,& ! 253
 'ITAB_W_MRLW'   ,'W2' ,'ITAB_W_MRLW',' ( )'                                ,& ! 254
 'ITAB_W_MROW'   ,'W2' ,'ITAB_W_MROW',' ( )'                                ,& ! 255
 'ITAB_W_MROWH'  ,'W2' ,'ITAB_W_MROWH',' ( )'                               ,& ! 256
 'ITAB_W_IM'     ,'W2' ,'ITAB_W_IM',' ( )'                                  ,& ! 257
 'ITAB_W_IU'     ,'W2' ,'ITAB_W_IU',' ( )'                                  ,& ! 258
 'ITAB_W_IV'     ,'W2' ,'ITAB_W_IV',' ( )'                                  ,& ! 259
 'ITAB_W_IW'     ,'W2' ,'ITAB_W_IW',' ( )'                                  ,& ! 260
 'ITAB_W_DIRU'   ,'W2' ,'ITAB_W_DIRU',' ( )'                                ,& ! 261
 'ITAB_W_DIRV'   ,'W2' ,'ITAB_W_DIRV',' ( )'                                ,& ! 262
 'ITAB_W_FWV'    ,'W2' ,'ITAB_W_FWV',' ( )'                                 ,& ! 263
 'ITAB_W_FWW'    ,'W2' ,'ITAB_W_FWW',' ( )'                                 ,& ! 264
 'ITAB_W_FWU'    ,'W2' ,'ITAB_W_FWU',' ( )'                                 ,& ! 265
 'ITAB_W_VXU'    ,'W2' ,'ITAB_W_VXU',' ( )'                                 ,& ! 266
 'ITAB_W_VYU'    ,'W2' ,'ITAB_W_VYU',' ( )'                                 ,& ! 267
 'ITAB_W_VZU'    ,'W2' ,'ITAB_W_VZU',' ( )'                                 ,& ! 268
 'ITAB_W_VXW'    ,'W2' ,'ITAB_W_VXW',' ( )'                                 ,& ! 269
 'ITAB_W_VYW'    ,'W2' ,'ITAB_W_VYW',' ( )'                                 ,& ! 270
 'ITAB_W_VZW'    ,'W2' ,'ITAB_W_VZW',' ( )'                                 ,& ! 271
 'ITAB_W_VXU_W'  ,'W2' ,'ITAB_W_VXU_W',' ( )'                               ,& ! 272
 'ITAB_W_VYU_W'  ,'W2' ,'ITAB_W_VYU_W',' ( )'                               ,& ! 273
 'ITAB_W_FARM'   ,'W2' ,'ITAB_W_FARM',' ( )'                                ,& ! 274
 'ITAB_W_FARV'   ,'W2' ,'ITAB_W_FARV',' ( )'                                ,& ! 275
 'ITAB_W_IWNUD'  ,'W2' ,'ITAB_W_IWNUD',' ( )'                               ,& ! 276
 'ITAB_W_FNUD'   ,'W2' ,'ITAB_W_FNUD',' ( )'                                 / ! 277

! TIME-AVERAGED FIELDS - 2D

data fldlib(1:4,278:331)/ &

 'PRESS_MAVG'    ,'T2','MONTH-AVG PRESSURE',' (hPa)'                        ,& ! 278
 'RHO_MAVG'      ,'T2','MONTH-AVG AIR DENSITY',' (kg m:S2:-3  )'            ,& ! 279
 'TEMPK_MAVG'    ,'T2','MONTH-AVG TEMPERATURE',' (K)'                       ,& ! 280
 'SH_V_MAVG'     ,'T2','MONTH-AVG VAPOR SPEC DENSITY',' (g kg:S2:-1  )'     ,& ! 281
 'SH_W_MAVG'     ,'T2','MONTH-AVG TOT WATER SPEC DENSITY',' (g kg:S2:-1  )' ,& ! 282
 'WC_MAVG'       ,'T2','MONTH-AVG W VELOCITY',' (m s:S2:-1  )'              ,& ! 283
 'ZONAL_WINDW_MAVG' ,'T2','MONTH-AVG ZONAL WIND',' (m s:S2:-1  )'           ,& ! 284
 'MERID_WINDW_MAVG' ,'T2','MONTH-AVG MERID WIND',' (m s:S2:-1  )'           ,& ! 285
 'RSHORT_MAVG'      ,'T2','MONTH-AVG SFC DOWNWARD S/W FLX',' (W m:S2:-2  )' ,& ! 286
 'RSHORT_TOP_MAVG'  ,'T2','MONTH-AVG TOP DOWNWARD S/W FLX',' (W m:S2:-2  )' ,& ! 287
 'RSHORTUP_MAVG'    ,'T2','MONTH-AVG SFC UPWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 288
 'RSHORTUP_TOP_MAVG','T2','MONTH-AVG TOP UPWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 289
 'RLONG_MAVG'       ,'T2','MONTH-AVG SFC DOWNWARD L/W FLX',' (W m:S2:-2  )' ,& ! 290
 'RLONGUP_MAVG'     ,'T2','MONTH-AVG SFC UPWARD L/W FLX',' (W m:S2:-2  )'   ,& ! 291
 'RLONGUP_TOP_MAVG' ,'T2','MONTH-AVG TOP UPWARD L/W FLX',' (W m:S2:-2  )'   ,& ! 292
 'LATFLUX_MAVG'     ,'T2','MONTH-AVG SFC LATENT HEAT FLX',' (W m:S2:-2  )'  ,& ! 293
 'SENSFLUX_MAVG'    ,'T2','MONTH-AVG SFC SENSIBLE HEAT FLX',' (W m:S2:-2  )',& ! 294
 'WINDSPEED_MAVG'   ,'T2','MONTH-AVG SFC WIND SPEED',' (m s:S2:-1  )'       ,& ! 295
 'ACCPMIC_MTOT'     ,'T2','MONTH-ACCUM MICPHYS PRECIP',' (kg m:S2:-2  )'    ,& ! 296
 'ACCPCON_MTOT'     ,'T2','MONTH-ACCUM CUPARM PRECIP',' (kg m:S2:-2  )'     ,& ! 297
 'ACCPBOTH_MTOT'    ,'T2','MONTH-ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )',& ! 298
 'PRESS_DAVG'       ,'T2','DAY-AVG SFC PRESSURE',' (hPa)'                   ,& ! 299
 'ZONAL_WINDW_DAVG' ,'T2','DAY-AVG SFC ZONAL WIND',' (m s:S2:-1  )'         ,& ! 300
 'MERID_WINDW_DAVG' ,'T2','DAY-AVG SFC MERID WIND',' (m s:S2:-1  )'         ,& ! 301
 'RSHORT_DAVG'      ,'T2','DAY-AVG SFC DOWNWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 302
 'TEMPK_DAVG'       ,'T2','DAY-AVG SFC TEMPERATURE',' (K)'                  ,& ! 303
 'TEMPK_DMIN'       ,'T2','DAY-MIN SFC TEMPERATURE',' (K)'                  ,& ! 304
 'TEMPK_DMAX'       ,'T2','DAY-MAX SFC TEMPERATURE',' (K)'                  ,& ! 305
 'ACCPMIC_DTOT'     ,'T2','DAY-ACCUM MICPHYS PRECIP',' (kg m:S2:-2  )'      ,& ! 306
 'ACCPCON_DTOT'     ,'T2','DAY-ACCUM CUPARM PRECIP',' (kg m:S2:-2  )'       ,& ! 307
 'ACCPBOTH_DTOT'    ,'T2','DAY-ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )'  ,& ! 308
 'PRESS_UL_DAVG'      ,'T2','DAY-AVG UL PRESSURE',' (hPa)'                  ,& ! 309
 'ZONAL_WINDW_UL_DAVG','T2','DAY-AVG UL ZONAL WIND',' (m s:S2:-1  )'        ,& ! 310
 'MERID_WINDW_UL_DAVG','T2','DAY-AVG UL MERID WIND',' (m s:S2:-1  )'        ,& ! 311
 'CANTEMPK_DAVG'    ,'B2','DAY-AVG CANOPY TEMPERATURE',' (K)'               ,& ! 312
 'CANTEMPK_DMIN'    ,'B2','DAY-MIN CANOPY TEMPERATURE',' (K)'               ,& ! 313
 'CANTEMPK_DMAX'    ,'B2','DAY-MAX CANOPY TEMPERATURE',' (K)'               ,& ! 314
 'VEGTEMPK_DAVG'    ,'L2','DAY-AVG VEG TEMPERATURE',' (K)'                  ,& ! 315
 'VEGTEMPK_DMIN'    ,'L2','DAY-MIN VEG TEMPERATURE',' (K)'                  ,& ! 316
 'VEGTEMPK_DMAX'    ,'L2','DAY-MAX VEG TEMPERATURE',' (K)'                  ,& ! 317
 'SOILTEMPK_DAVG'   ,'L2','DAY-AVG SOIL TEMPERATURE',' (K)'                 ,& ! 318
 'SOILTEMPK_DMIN'   ,'L2','DAY-MIN SOIL TEMPERATURE',' (K)'                 ,& ! 319
 'SOILTEMPK_DMAX'   ,'L2','DAY-MAX SOIL TEMPERATURE',' (K)'                 ,& ! 320
 'FCELL_AIRTEMPK_DAVG','F2','DAY-AVG FLUXCELL ATM TEMP',' (K)'                ,& ! 321
 'FCELL_AIRTEMPK_DMIN','F2','DAY-MIN FLUXCELL ATM TEMP',' (K)'                ,& ! 322
 'FCELL_AIRTEMPK_DMAX','F2','DAY-MAX FLUXCELL ATM TEMP',' (K)'                ,& ! 323
 'FCELL_CANTEMPK_DAVG','F2','DAY-AVG FLUXCELL CAN TEMP',' (K)'                ,& ! 324
 'FCELL_CANTEMPK_DMIN','F2','DAY-MIN FLUXCELL CAN TEMP',' (K)'                ,& ! 325
 'FCELL_CANTEMPK_DMAX','F2','DAY-MAX FLUXCELL CAN TEMP',' (K)'                ,& ! 326
 'FCELL_SENSFLUX_DAVG','F2','DAY-AVG FLUXCELL SENS HEAT FLUX',' (W m:S2:-2  )',& ! 327
 'FCELL_LATFLUX_DAVG' ,'F2','DAY-AVG FLUXCELL LAT HEAT FLUX',' (W m:S2:-2  )' ,& ! 328
 'FTOA_SENSFLUX_DAVG' ,'T2','DAY-AVG FTOA SENS HEAT FLUX',' (W m:S2:-2  )'    ,& ! 329
 'FTOA_LATFLUX_DAVG'  ,'T2','DAY-AVG FTOA LAT HEAT FLUX',' (W m:S2:-2  )'     ,& ! 330
 'PRESS_AVG24'      ,'T2','MONTH-AVG HOURLY PRESSURE',' (hPa)'               / ! 331

data fldlib(1:4,332:341)/ &
 'RSHORT_AVG'      ,'T2' ,'AVG SFC DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'   ,& ! 332
 'RSHORTUP_AVG'    ,'T2' ,'AVG SFC UPWARD SHORTWV FLX',' (W m:S2:-2  )'     ,& ! 333
 'RLONG_AVG'       ,'T2' ,'AVG SFC DOWNWARD LONGWV FLX',' (W m:S2:-2  )'    ,& ! 334
 'RLONGUP_AVG'     ,'T2' ,'AVG SFC UPWARD LONGWV FLX',' (W m:S2:-2  )'      ,& ! 335
 'RSHORT_TOP_AVG'  ,'T2' ,'AVG TOP DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'   ,& ! 336
 'RSHORTUP_TOP_AVG','T2' ,'AVG TOP UPWARD SHORTWV FLX',' (W m:S2:-2  )'     ,& ! 337
 'RLONGUP_TOP_AVG' ,'T2' ,'AVG TOP UPWARD LONGWV FLX',' (W m:S2:-2  )'      ,& ! 338
 'SENSFLUX_AVG'    ,'T2' ,'AVG ATM SFC SENS HEAT FLUX',' (N m:S2:-2  )'     ,& ! 339
 'LATFLUX_AVG'     ,'T2' ,'AVG ATM SFC LAT HEAT FLUX',' (W m:S2:-2  )'      ,& ! 340
 'VAPFLUX_AVG'     ,'T2' ,'AVG ATM SFC VAP FLUX',' (kg m:S2:-2   s:S2:-1  )' / ! 341

! Miscellaneous and new additions

data fldlib(1:4,342:362)/ &
 'RHO_OBS'       ,'T3' ,'NUDGING OBS AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 342
 'THETA_OBS'     ,'T3' ,'NUDGING OBS THETA',' (K)'                          ,& ! 343
 'SHW_OBS'       ,'T3' ,'NUDGING OBS VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 344
 'UZONAL_OBS'    ,'T3' ,'NUDGING OBS ZONAL WIND',' (m s:S2:-1  )'           ,& ! 345
 'UMERID_OBS'    ,'T3' ,'NUDGING OBS MERID WIND',' (m s:S2:-1  )'           ,& ! 346
 'RHO_SIM'       ,'T3' ,'NUDGING SIM AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 347
 'THETA_SIM'     ,'T3' ,'NUDGING SIM THETA',' (K)'                          ,& ! 348
 'SHW_SIM'       ,'T3' ,'NUDGING SIM VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 349
 'UZONAL_SIM'    ,'T3' ,'NUDGING SIM ZONAL WIND',' (m s:S2:-1  )'           ,& ! 350
 'UMERID_SIM'    ,'T3' ,'NUDGING SIM MERID WIND',' (m s:S2:-1  )'           ,& ! 351
 'RHO_OBS_SIM'   ,'T3' ,'NUDGING DIF AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 352
 'THETA_OBS_SIM' ,'T3' ,'NUDGING DIF THETA',' (K)'                          ,& ! 353
 'SHW_OBS_SIM'   ,'T3' ,'NUDGING DIF VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 354
 'UZONAL_OBS_SIM','T3' ,'NUDGING DIF ZONAL WIND',' (m s:S2:-1  )'           ,& ! 355
 'UMERID_OBS_SIM','T3' ,'NUDGING DIF MERID WIND',' (m s:S2:-1  )'           ,& ! 356
 'VXE'           ,'T3' ,'EARTH CARTESIAN X WIND',' (m s:S2:-1  )'           ,& ! 357
 'VYE'           ,'T3' ,'EARTH CARTESIAN Y WIND',' (m s:S2:-1  )'           ,& ! 358
 'VZE'           ,'T3' ,'EARTH CARTESIAN Z WIND',' (m s:S2:-1  )'           ,& ! 359
 'PBLH'          ,'T2' ,'PBL HEIGHT',' (m)'                                 ,& ! 360
 'HKM'           ,'T3' ,'EDDY DIFFUSIVITY',' (m:S2:2 s:S2:-1  )'            ,& ! 361
 'HEAD0'         ,'L2' ,'HEAD0',' (m)'                                       / ! 362

! External fields

data fldlib(1:4,363:365)/ &
 'VORTP'         ,'P3' ,'VORTP',' (s:S2:-1  )'                              ,& ! 363
 'VORTN'         ,'N3' ,'VORTN',' (s:S2:-1  )'                              ,& ! 364
 'RKE'           ,'T3' ,'RKE',' (s:S2:-1  )'                                 / ! 365

if (icall /= 1) then
   icall = 1
   allocate (aux(mwa))
endif

if (infotyp == 'UNITS') then

! Print field name to be plotted

   write(io6,*) 'oplib ',trim(fldname0)

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

   if (fldname(1:6) == 'ITAB_V' .and. meshtype == 1) then
      write(io6,*) 'Plot field ',trim(fldname),' not available for hexagonal mesh.'
      go to 1000
   elseif (fldname(1:6) == 'ITAB_U' .and. meshtype == 2) then
      write(io6,*) 'Plot field ',trim(fldname),' not available for triangular mesh.'
      go to 1000
   endif

   op%stagpt = fldlib(2,icase)(1:1)
   op%dimens = fldlib(2,icase)(2:3)
   op%label  = fldlib(3,icase)
   if (indp > 0) op%label = trim(fldname0)  
   op%units  = fldlib(4,icase)

   op%fldval_min = 1.e30
   op%fldval_max = -1.e30
   op%fldvalv_min = 1.e30
   op%fldvalv_max = -1.e30

endif

if (infotyp == 'VALUV') then
   op%stagpt = 'V'
   op%dimens = '3'
   if (trim(fldname0) == 'UC') then
      icase = 6
   else
      icase = 7
   endif
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

! Execute IF block below even when infotyp == 'UNITS'
! in order to check whether current plot field is available in this model run.
! For this type of call to this subroutine, (k,i) are passed in as (1,2).

! ATMOSPHERE - 3D

select case(icase)

case(1) ! 'UMC'

   if (.not. allocated(umc)) go to 1000

   fldval = wtbot * umc(k  ,i) &
          + wttop * umc(k+1,i)

case(2) ! 'VMC'

   if (.not. allocated(vmc)) go to 1000

   fldval = wtbot * vmc(k  ,i) &
          + wttop * vmc(k+1,i)

case(3) ! 'WMC'

! Need to re-examine use of k-1 when k = lpw

   fldval = wtbot * wmc(k  ,i) &
          + wttop * wmc(k+1,i)

case(4) ! 'UMP'

   if (.not. allocated(ump)) go to 1000

   fldval = wtbot * ump(k  ,i) &
          + wttop * ump(k+1,i)

case(5) ! 'VMP'

   if (.not. allocated(vmp)) go to 1000

   fldval = wtbot * vmp(k  ,i) &
          + wttop * vmp(k+1,i)

case(6) ! 'UC'

   if (.not. allocated(uc)) go to 1000

   fldval = wtbot * uc(k  ,i) &
          + wttop * uc(k+1,i)

case(7) ! 'VC'

   if (.not. allocated(vc)) go to 1000

   fldval = wtbot * vc(k  ,i) &
          + wttop * vc(k+1,i)

case(8) ! 'WC'

   fldval = wtbot * wc(k  ,i) &
          + wttop * wc(k+1,i)

case(9) ! 'RHO'

   fldval = wtbot * rho(k  ,i) &
          + wttop * rho(k+1,i)

case(10) ! 'PRESS'

   fldval = press(k,i) * .01

   if (nl%test_case == 2 .or. nl%test_case == 5) then
      fldval = press(k,i)
   endif

case(11) ! 'THIL'

   fldval = wtbot * thil(k  ,i) &
          + wttop * thil(k+1,i)

case(12) ! 'THETA'

   fldval = wtbot * theta(k  ,i) &
          + wttop * theta(k+1,i)

case(13) ! 'AIRTEMPK'

   fldval = wtbot * tair(k  ,i) &
          + wttop * tair(k+1,i)

case(14) ! 'AIRTEMPC'

   fldval = wtbot * tair(k  ,i) &
          + wttop * tair(k+1,i) - 273.15

case(15) ! 'SH_W'

   fldval = (wtbot * sh_w(k  ,i) &
          +  wttop * sh_w(k+1,i)) * 1.e3

case(16) ! 'SH_V'

   fldval = (wtbot * sh_v(k  ,i) &
          +  wttop * sh_v(k+1,i)) * 1.e3

case(17) ! 'SH_C'

   if (.not. allocated(sh_c)) go to 1000

   fldval = (wtbot * sh_c(k  ,i) &
          +  wttop * sh_c(k+1,i)) * 1.e3

case(18) ! 'SH_D'

   if (.not. allocated(sh_d)) go to 1000

   fldval = (wtbot * sh_d(k  ,i) &
          +  wttop * sh_d(k+1,i)) * 1.e3

case(19) ! 'SH_R'

   if (.not. allocated(sh_r)) go to 1000

   fldval = (wtbot * sh_r(k  ,i) &
          +  wttop * sh_r(k+1,i)) * 1.e3

case(20) ! 'SH_P'

   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbot * sh_p(k  ,i) &
          +  wttop * sh_p(k+1,i)) * 1.e3

case(21) ! 'SH_S'

   if (.not. allocated(sh_s)) go to 1000

   fldval = (wtbot * sh_s(k  ,i) &
          +  wttop * sh_s(k+1,i)) * 1.e3

case(22) ! 'SH_A'

   if (.not. allocated(sh_a)) go to 1000

   fldval = (wtbot * sh_a(k  ,i) &
          +  wttop * sh_a(k+1,i)) * 1.e3

case(23) ! 'SH_G'

   if (.not. allocated(sh_g)) go to 1000

   fldval = (wtbot * sh_g(k  ,i) &
          +  wttop * sh_g(k+1,i)) * 1.e3

case(24) ! 'SH_H'

   if (.not. allocated(sh_h)) go to 1000

   fldval = (wtbot * sh_h(k  ,i) &
          +  wttop * sh_h(k+1,i)) * 1.e3

case(25) ! 'SH_CP'

   if (.not. allocated(sh_c)) go to 1000
   if (.not. allocated(sh_p)) go to 1000

   fldval = (wtbot * (sh_c(k  ,i) + sh_p(k  ,i)) &
          +  wttop * (sh_c(k+1,i) + sh_p(k+1,i))) * 1.e3

case(26) ! 'SH_TOTCOND'

   fldval = (wtbot * (sh_w(k  ,i) - sh_v(k  ,i)) &
          +  wttop * (sh_w(k+1,i) - sh_v(k+1,i))) * 1.e3

case(27) ! 'CON_C'

   if (.not. allocated(con_c)) go to 1000

   fldval = (wtbot * con_c(k  ,i) &
          +  wttop * con_c(k+1,i)) * 1.e-6

case(28) ! 'CON_D'

   if (.not. allocated(con_d)) go to 1000

   fldval = (wtbot * con_d(k  ,i) &
          +  wttop * con_d(k+1,i)) * 1.e-6

case(29) ! 'CON_R'

   if (.not. allocated(con_r)) go to 1000

   fldval = wtbot * con_r(k  ,i) &
          + wttop * con_r(k+1,i)

case(30) ! 'CON_P'

   if (.not. allocated(con_p)) go to 1000

   fldval = wtbot * con_p(k  ,i) &
          + wttop * con_p(k+1,i)

case(31) ! 'CON_S'

   if (.not. allocated(con_s)) go to 1000

   fldval = wtbot * con_s(k  ,i) &
          + wttop * con_s(k+1,i)

case(32) ! 'CON_A'

   if (.not. allocated(con_a)) go to 1000

   fldval = wtbot * con_a(k  ,i) &
          + wttop * con_a(k+1,i)

case(33) ! 'CON_G'

   if (.not. allocated(con_g)) go to 1000

   fldval = wtbot * con_g(k  ,i) &
          + wttop * con_g(k+1,i)

case(34) ! 'CON_H'

   if (.not. allocated(con_h)) go to 1000

   fldval = wtbot * con_h(k  ,i) &
          + wttop * con_h(k+1,i)

case(35) ! 'CON_CCN'

   if (.not. allocated(con_ccn)) go to 1000

   fldval = (wtbot * con_ccn(k  ,i) &
          +  wttop * con_ccn(k+1,i)) * 1.e-6

case(36) ! 'CON_GCCN'

   if (.not. allocated(con_gccn)) go to 1000

   fldval = (wtbot * con_gccn(k  ,i) &
          +  wttop * con_gccn(k+1,i)) * 1.e-6

case(37) ! 'CON_IFN'

   if (.not. allocated(con_ifn)) go to 1000

   fldval = wtbot * con_ifn(k  ,i) &
          + wttop * con_ifn(k+1,i)

case(38) ! 'VKM'

   fldval = wtbot * vkm(k  ,i) &
          + wttop * vkm(k+1,i)

case(39) ! 'FTHRD'

   if (.not. allocated(fthrd_sw)) goto 1000
   if (.not. allocated(fthrd_lw)) goto 1000

   fldval = wtbot * (fthrd_sw(k  ,i) + fthrd_sw(k  ,i)) &
          + wttop * (fthrd_sw(k+1,i) + fthrd_lw(k+1,i))

case(40:45) ! 'SPEEDV','AZIMV','ZONAL_WINDV','MERID_WINDV','ZONAL_WINDV_P','MERID_WINDV_P'

   if (meshtype == 1) then
      iw1 = itab_u(i)%iw(1)
      iw2 = itab_u(i)%iw(2)
   else
      iw1 = itab_v(i)%iw(1)
      iw2 = itab_v(i)%iw(2)
   endif

   vx = wtbot * 0.5 * (vxe(k  ,iw1) + vxe(k  ,iw2)) &
      + wttop * 0.5 * (vxe(k+1,iw1) + vxe(k+1,iw2))

   vy = wtbot * 0.5 * (vye(k  ,iw1) + vye(k  ,iw2)) &
      + wttop * 0.5 * (vye(k+1,iw1) + vye(k+1,iw2))

   vz = wtbot * 0.5 * (vze(k  ,iw1) + vze(k  ,iw2)) &
      + wttop * 0.5 * (vze(k+1,iw1) + vze(k+1,iw2))

   if (trim(fldname) == 'SPEEDV') then

      fldval = sqrt(vx** 2 + vy** 2 + vz** 2)

   else

      u = 0.
      v = 0.

      if (meshtype == 1) then

         raxis = sqrt(xeu(i) ** 2 + yeu(i) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            u = (vy * xeu(i) - vx * yeu(i)) / raxis
            v = vz * raxis / erad &
              - (vx * xeu(i) + vy * yeu(i)) * zeu(i) / (raxis * erad) 
         else
            u = 0.0
            v = 0.0
         endif

      elseif (mdomain < 2) then

         raxis = sqrt(xev(i) ** 2 + yev(i) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            u = (vy * xev(i) - vx * yev(i)) / raxis
            v = vz * raxis / erad &
              - (vx * xev(i) + vy * yev(i)) * zev(i) / (raxis * erad)
         else
            u = 0.0
            v = 0.0
         endif

      else
         
         u = vx
         v = vy

      endif

      if (trim(fldname) == 'AZIMV'           ) then
         fldval = mod(450. - piu180 * atan2(v,u),360.)
      elseif (trim(fldname) == 'ZONAL_WINDV' ) then
         fldval = u
      elseif (trim(fldname) == 'MERID_WINDV' ) then
         fldval = v
      else

! Tangential component of initial wind at V point

         ucc_init = wtbot * uc_init(k  ,i) &
                  + wttop * uc_init(k+1,i)

         vcc_init = wtbot * vc_init(k  ,i) &
                  + wttop * vc_init(k+1,i)

         vx_init = unx(i) * ucc_init + vnx(i) * vcc_init
         vy_init = uny(i) * ucc_init + vny(i) * vcc_init
         vz_init = unz(i) * ucc_init + vnz(i) * vcc_init

         u_init = 0.
         v_init = 0.

         if (meshtype == 1) then

            raxis = sqrt(xeu(i) ** 2 + yeu(i) ** 2)  ! dist from earth axis

            if (raxis > 1.e3) then
               u_init = (vy_init * xeu(i) - vx_init * yeu(i)) / raxis
               v_init = vz_init * raxis / erad &
                  - (vx_init * xeu(i) + vy_init * yeu(i)) * zeu(i) / (raxis * erad) 
            else
               u_init = 0.0
               v_init = 0.0
            endif

         elseif (mdomain < 2) then

            raxis = sqrt(xev(i) ** 2 + yev(i) ** 2)  ! dist from earth axis

            if (raxis > 1.e3) then
               u_init = (vy_init * xev(i) - vx_init * yev(i)) / raxis
               v_init = vz_init * raxis / erad &
                  - (vx_init * xev(i) + vy_init * yev(i)) * zev(i) / (raxis * erad) 
            else
               u_init = 0.0
               v_init = 0.0
            endif

         else

            u_init = vx_init
            v_init = vy_init

         endif

         if     (trim(fldname) == 'ZONAL_WINDV_P') then
            fldval = u - u_init
         elseif (trim(fldname) == 'MERID_WINDV_P') then
            fldval = v - v_init
         endif

      endif

   endif

case(46:49) ! 'SPEEDW','AZIMW','ZONAL_WINDW','MERID_WINDW'

   npoly = itab_w(i)%npoly
   rpolyi = 1. / real(npoly)

   vx = wtbot * vxe(k  ,i) &
      + wttop * vxe(k+1,i)

   vy = wtbot * vye(k  ,i) &
      + wttop * vye(k+1,i)

   vz = wtbot * vze(k  ,i) &
      + wttop * vze(k+1,i)

   if (trim(fldname) == 'SPEEDW') then
      fldval = sqrt(vx** 2 + vy** 2 + vz** 2)
   else
      
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

      if (trim(fldname) == 'AZIMW') then
         fldval = mod(450. - piu180 * atan2(v,u),360.)
      elseif (trim(fldname) == 'ZONAL_WINDW') then
         fldval = u
      elseif (trim(fldname) == 'MERID_WINDW') then
         fldval = v
      endif
   endif

case(50:52) ! 'RVORTZM','RVORTZM_P','TVORTZM'

   fldval = 0.

   do j = 1,itab_m(i)%npoly

      if (meshtype == 1) then

         iu = itab_m(i)%iu(j)

         iu1 = itab_u(iu)%iu(1)
         iu2 = itab_u(iu)%iu(2)
         iu3 = itab_u(iu)%iu(3)
         iu4 = itab_u(iu)%iu(4)

         iw1 = itab_u(iu)%iw(1)
         iw2 = itab_u(iu)%iw(2)

         tuu1 = itab_u(iu)%tuu(1)
         tuu2 = itab_u(iu)%tuu(2)
         tuu3 = itab_u(iu)%tuu(3)
         tuu4 = itab_u(iu)%tuu(4)

         ucc = wtbot * uc(k  ,iu) &
             + wttop * uc(k+1,iu)

         vcc = wtbot * (uc(k  ,iu1) * tuu1  &
                      + uc(k  ,iu2) * tuu2  &
                      + uc(k  ,iu3) * tuu3  &
                      + uc(k  ,iu4) * tuu4) &
             + wttop * (uc(k+1,iu1) * tuu1  &
                      + uc(k+1,iu2) * tuu2  &
                      + uc(k+1,iu3) * tuu3  &
                      + uc(k+1,iu4) * tuu4)
      else

         iv = itab_m(i)%iv(j)

         vcc = wtbot * vc(k  ,iv) &
             + wttop * vc(k+1,iv)

      endif
      
      if (fldname == 'RVORTZM_P') then

         if (meshtype == 1) then

            ucc_init = wtbot * uc_init(k  ,iu) &
                     + wttop * uc_init(k+1,iu)

            vcc_init = wtbot * (uc_init(k  ,iu1) * tuu1  &
                              + uc_init(k  ,iu2) * tuu2  &
                              + uc_init(k  ,iu3) * tuu3  &
                              + uc_init(k  ,iu4) * tuu4) &
                     + wttop * (uc_init(k+1,iu1) * tuu1  &
                              + uc_init(k+1,iu2) * tuu2  &
                              + uc_init(k+1,iu3) * tuu3  &
                              + uc_init(k+1,iu4) * tuu4)
            
            ucc = ucc - ucc_init
            vcc = vcc - vcc_init

         else

            vcc_init = wtbot * vc_init(k  ,iv) &
                     + wttop * vc_init(k+1,iv)

            vcc = vcc - vcc_init

         endif

      endif

! Now reconstruct total wind vector projected onto vector from IW1 to IW2

      if (meshtype == 1) then

         contrib = ucc * (unx(iu) * (xew(iw2) - xew(iw1))  &
                       +  uny(iu) * (yew(iw2) - yew(iw1))  &
                       +  unz(iu) * (zew(iw2) - zew(iw1))) &
                 + vcc * (vnx(iu) * (xew(iw2) - xew(iw1))  &
                       +  vny(iu) * (yew(iw2) - yew(iw1))  &
                       +  vnz(iu) * (zew(iw2) - zew(iw1)))

         if (i == itab_u(iu)%im(1)) then
            fldval = fldval + contrib
         else
            fldval = fldval - contrib
         endif

      else

         if (i == itab_v(iv)%im(2)) then
            fldval = fldval + dnv(iv) * vcc
         else
            fldval = fldval - dnv(iv) * vcc
         endif

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

case(53) ! 'DIVERG'

   fldval = 0.

   npoly = itab_w(i)%npoly

   if (itab_w(i)%iwp == i) then

   do j = 1,npoly

      if (meshtype == 1) then

         iv = itab_w(i)%iu(j)

         fldval = fldval &

                + wtbot * umc(k,iv) * (-itab_w(i)%diru(j)) * aru(k,iv) &
                        / (volt(k,i) * rho(k,i)) &

                + wttop * umc(k+1,iv) * (-itab_w(i)%diru(j)) * aru(k+1,iv) &
                        / (volt(k+1,i) * rho(k+1,i))

      else

         iv = itab_w(i)%iv(j)

         fldval = fldval &

                + wtbot * vmc(k,iv) * (-itab_w(i)%dirv(j)) * arv(k,iv) &
                        / (volt(k,i) * rho(k,i)) &

                + wttop * vmc(k+1,iv) * (-itab_w(i)%dirv(j)) * arv(k+1,iv) &
                        / (volt(k+1,i) * rho(k+1,i))
      endif

   enddo
   endif

case(54) ! 'UMASSFLUX'

   if (meshtype == 2) goto 1000
   fldval = umc(k,i) * aru(k,i)

case(55) ! 'VMASSFLUX'

   if (meshtype == 1) goto 1000
   fldval = vmc(k,i) * arv(k,i)

case(56) ! 'UC_P'

   if (meshtype == 2) goto 1000
   fldval = uc(k,i) - uc_init(k,i)

case(57) ! 'VC_P'

   if (meshtype == 1) goto 1000
   fldval = vc(k,i) - vc_init(k,i)

case(58) ! 'PRESS_P'

   fldval = press(k,i) - press_init(k,i)

case(59) ! 'RHO_P'

   fldval = wtbot * (rho(k  ,i) - rho_init(k  ,i)) &
          + wttop * (rho(k+1,i) - rho_init(k+1,i))

! For shallow water test case 5, define RHO_P from reference fields

   if (nl%test_case == 5) then
      call npr_bicubics(zanal00_swtc5,glatw(i),glonw(i),zanal0_swtc5)
      call npr_bicubics(zanal15_swtc5,glatw(i),glonw(i),zanal_swtc5)

      fldval = (rho(k,i) - rho_init(k,i)) - (zanal_swtc5 - zanal0_swtc5)
   endif

case(60) ! 'THETA_P'

   fldval = wtbot * (theta(k  ,i) - theta_init(k  ,i)) &
          + wttop * (theta(k+1,i) - theta_init(k+1,i))

case(61) ! 'AIRTEMPK_P'

   fldval = wtbot * tair(k  ,i) &
          + wttop * tair(k+1,i) &
          - wtbot * theta_init(k  ,i) * (press_init(k  ,i) * p00i) ** rocp &
          - wttop * theta_init(k+1,i) * (press_init(k+1,i) * p00i) ** rocp

case(62) ! 'UMT'

   if (meshtype == 2) goto 1000
   fldval = umt(k,i) * volui(k,i)

case(63) ! 'VMT'

   if (meshtype == 1) goto 1000
   fldval = vmt(k,i) * volvi(k,i)

case(64) ! 'WMT'

   fldval = wmt(k,i) * volwi(k,i)

case(65) ! 'ADDSC'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000

   fldval = wtbot * addsc(indp)%sclp(k  ,i) &
          + wttop * addsc(indp)%sclp(k+1,i)

case(66) ! 'ADDSC_P'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000
   if (.not. allocated(addsc1_init)) go to 1000

   fldval = addsc(indp)%sclp(k,i) - addsc1_init(k,i) 
   
! ATMOSPHERE - 2D

case(67) ! 'ZPLEV'

   fldval = wtbot * zt(k  ) &
          + wttop * zt(k+1)

case(68) ! 'RSHORT_TOP'

   if (.not. allocated(rshort_top)) go to 1000

   fldval = rshort_top(i)

case(69) ! 'RSHORTUP_TOP'

   if (.not. allocated(rshortup_top)) go to 1000

   fldval = rshortup_top(i)

case(70) ! 'RLONGUP_TOP'

   if (.not. allocated(rlongup_top)) go to 1000

   fldval = rlongup_top(i)

! ATMOSPHERE SURFACE - 2D

case(71) ! 'RSHORT'

   if (.not. allocated(rshort)) go to 1000

   fldval = rshort(i)

case(72) ! 'RSHORTUP'

   if (.not. allocated(rshort)) go to 1000
   if (.not. allocated(albedt)) go to 1000

   fldval = rshort(i) * albedt(i)

case(73) ! 'RLONG'

   if (.not. allocated(rlong)) go to 1000

   fldval = rlong(i)

case(74) ! 'RLONGUP'

   if (.not. allocated(rlongup)) go to 1000

   fldval = rlongup(i)

case(75) ! 'ALBEDT'

   if (.not. allocated(albedt)) go to 1000

   fldval = albedt(i)

case(76) ! 'VKM_SFC'

   fldval = vkm_sfc(1,i)

case(77) ! 'SFLUX_W'

   fldval = sflux_w(i)

case(78) ! 'SENSFLUX'

   fldval = sflux_t(i) * cp
   
case(79) ! 'VAPFLUX'

   fldval = sflux_r(i)

case(80) ! 'LATFLUX'

   fldval = sflux_r(i) * alvl

case(81) ! 'PCPRD'

   if (.not. allocated(pcprd)) go to 1000

   fldval = pcprd(i) * 3600.

case(82) ! 'PCPRR'

   if (.not. allocated(pcprr)) go to 1000

   fldval = pcprr(i) * 3600.

case(83) ! 'PCPRP'

   if (.not. allocated(pcprp)) go to 1000

   fldval = pcprp(i)  * 3600.   

case(84) ! 'PCPRS'

   if (.not. allocated(pcprs)) go to 1000

   fldval = pcprs(i) * 3600.

case(85) ! 'PCPRA'

   if (.not. allocated(pcpra)) go to 1000

   fldval = pcpra(i) * 3600.

case(86) ! 'PCPRG'

   if (.not. allocated(pcprg)) go to 1000

   fldval = pcprg(i) * 3600.

case(87) ! 'PCPRH'

   if (.not. allocated(pcprh)) go to 1000

   fldval = pcprh(i) * 3600.

case(88) ! 'PCPRMIC'

   fldval = 0.

   if (allocated(pcprd)) fldval = fldval + pcprd(i)
   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)

   fldval = fldval * 3600.

case(89) ! 'PCPRCON'

   if (.not. allocated(conprr)) go to 1000

   fldval = conprr(i) * 3600.

case(90) ! 'PCPRBOTH'

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

case(91) ! 'ACCPD'

   if (.not. allocated(accpd)) go to 1000

   fldval = accpd(i)

case(92) ! 'ACCPR'

   if (.not. allocated(accpr)) go to 1000

   fldval = accpr(i)

case(93) ! 'ACCPP'

   if (.not. allocated(accpp)) go to 1000

   fldval = accpp(i)

case(94) ! 'ACCPS'

   if (.not. allocated(accps)) go to 1000

   fldval = accps(i)

case(95) ! 'ACCPA'

   if (.not. allocated(accpa)) go to 1000

   fldval = accpa(i)

case(96) ! 'ACCPG'

   if (.not. allocated(accpg)) go to 1000

   fldval = accpg(i)

case(97) ! 'ACCPH'

   if (.not. allocated(accph)) go to 1000

   fldval = accph(i)

case(98) ! 'ACCPMIC'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

case(99) ! 'ACCPCON'

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

case(100) ! 'ACCPBOTH'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

case(101) ! 'PCPDIF2MIC'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 hour

   fldval = fldval - accpmic_prev1(i)
   if (time8 - time8_prev1 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev1)
   endif

case(102) ! 'PCPDIF2CON'

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 hour

   fldval = fldval - accpcon_prev1(i)
   if (time8 - time8_prev1 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev1)
   endif

case(103) ! 'PCPDIF2BOTH'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 hour

   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)

   fldval = fldval - accpboth_prev1
   if (time8 - time8_prev1 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev1)
   endif

case(104) ! 'PCPDIF4MIC'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 hour

   fldval = fldval - accpmic_prev1(i) - (accpmic_prev2(i) - accpmic_prev3(i))
   if (time8 - time8_prev2 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev2)
   endif

case(105) ! 'PCPDIF4CON'

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 hour

   fldval = fldval - accpcon_prev1(i) - (accpcon_prev2(i) - accpcon_prev3(i))
   if (time8 - time8_prev2 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev2)
   endif

case(106) ! 'PCPDIF4BOTH'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 hour

   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)
   accpboth_prev2 = accpmic_prev2(i) + accpcon_prev2(i)
   accpboth_prev3 = accpmic_prev3(i) + accpcon_prev3(i)

   fldval = fldval - accpboth_prev1 - (accpboth_prev2 - accpboth_prev3)
   if (time8 - time8_prev2 > 3599.9) then
      fldval = fldval * 86400. / (time8 - time8_prev2)
   endif

case(107) ! 'PCPREL4MIC'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

! Compute relative differences involving 3 previously-stored fields

   fldval = fldval - accpmic_prev1(i) - (accpmic_prev2(i) - accpmic_prev3(i))
   denom = fldval + accpmic_prev1(i) - (accpmic_prev2(i) + accpmic_prev3(i))
   if (abs(denom) > 1.e-3) then
      fldval = fldval / denom
   else
      fldval = 0.
   endif

case(108) ! 'PCPREL4CON'

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

! Compute relative differences involving 3 previously-stored fields

   fldval = fldval - accpcon_prev1(i) - (accpcon_prev2(i) - accpcon_prev3(i))
   denom = fldval + accpcon_prev1(i) - (accpcon_prev2(i) + accpcon_prev3(i))
   if (abs(denom) > 1.e-3) then
      fldval = fldval / denom
   else
      fldval = 0.
   endif

case(109) ! 'PCPREL4BOTH'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

! Compute relative differences involving 3 previously-stored fields

   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)
   accpboth_prev2 = accpmic_prev2(i) + accpcon_prev2(i)
   accpboth_prev3 = accpmic_prev3(i) + accpcon_prev3(i)

   fldval = fldval - accpboth_prev1 - (accpboth_prev2 - accpboth_prev3)
   denom = fldval + accpboth_prev1 - (accpboth_prev2 + accpboth_prev3)
   if (abs(denom) > 1.e-3) then
      fldval = fldval / denom
   else
      fldval = 0.
   endif

! LAND_CELLS - 3D

case(110) ! 'SOIL_TEXT'

   fldval = real(land%ntext_soil(k,i))

case(111) ! 'SOIL_ENERGY'

   fldval = land%soil_energy(k,i) * 1.e-6

case(112) ! 'SOIL_TEMPC'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = tempk - 273.15

case(113) ! 'SOIL_FRACLIQ'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = fracliq

case(114) ! 'SOIL_WATER'

   fldval = land%soil_water(k,i) / slmsts(land%ntext_soil(k,i))

case(115) ! 'SFWAT_MASS'

   fldval = land%sfcwater_mass(k,i)

case(116) ! 'SFWAT_ENERGY'

   fldval = land%sfcwater_energy(k,i) * 1.e-3

case(117) ! 'SFWAT_TEMPC'

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = tempk - 273.15

case(118) ! 'SFWAT_FRACLIQ'

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = fracliq

case(119) ! 'SFWAT_DEPTH'

   fldval = 0.
   do klev = 1,land%nlev_sfcwater(i)
      fldval = fldval + land%sfcwater_depth(klev,i)
   enddo

! LAND_CELLS - 2D

case(120) ! 'NLEV_SFWAT'

   fldval = real(land%nlev_sfcwater(i))

case(121) ! 'VEG_NDVIC'

   fldval = land%veg_ndvic(i)

case(122) ! 'VEG_TEMPC'

   fldval = land%veg_temp(i) - 273.15

case(123) ! 'VEG_WATER'

   fldval = land%veg_water(i)

case(124) ! 'STOM_RESIST'

   fldval = land%stom_resist(i)

case(125) ! 'GROUND_SHV'

   fldval = land%ground_shv(i) * 1.e3

case(126) ! 'SOIL_DEPTH'

   fldval = -slz(1)

! SEA_CELLS - 2D

case(127) ! 'SEATP'

   fldval = sea%seatp(i)   

case(128) ! 'SEATF'

   fldval = sea%seatf(i)

case(129) ! 'SEATC'

   fldval = sea%seatc(i)

case(130) ! 'SEAICEP'

   fldval = sea%seaicep(i)   

case(131) ! 'SEAICEF'

   fldval = sea%seaicef(i)

case(132) ! 'SEAICEC'

   fldval = sea%seaicec(i)

! LAND AND SEA CELLS - 2D

case(133) ! 'LEAF_CLASS'

   if (op%stagpt == 'S') then
      fldval = real(sea%leaf_class(i))
   elseif (op%stagpt == 'L') then
      fldval = real(land%leaf_class(i))
   endif

case(134) ! 'AREA'

   if (op%stagpt == 'S') then
      fldval = sea%area(i)
   elseif (op%stagpt == 'L') then
      fldval = land%area(i)
   endif

case(135) ! 'ROUGH'

   if (op%stagpt == 'S') then
      fldval = sea%rough(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rough(i)
   endif

case(136) ! 'CAN_TEMPC'

   if (op%stagpt == 'S') then
      fldval = sea%can_temp(i) - 273.15
   elseif (op%stagpt == 'L') then
      fldval = land%can_temp(i) - 273.15
   endif

case(137) ! 'CAN_SHV'

   if (op%stagpt == 'S') then
      fldval = sea%can_shv(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%can_shv(i) * 1.e3
   endif

case(138) ! 'SFC_TEMPC'

   if (op%stagpt == 'S') then
      nls = sea%nlev_seaice(i)
      
      if (nls > 0) then
         fldval = (1.0 - sea%seaicec(i)) * (sea%seatc(i) - 273.15) &
                +        sea%seaicec(i)  * (sea%seaice_tempk(nls,i) - 273.15)
      else
         fldval = sea%seatc(i) - 273.15
      endif

   elseif (op%stagpt == 'L') then
      nls = land%nlev_sfcwater(i)

      if (nls > 0) then
         call qtk(land%sfcwater_energy(nls,i),tempk,fracliq)
         fldval = tempk - 273.15
      else
         call qwtk(land%soil_energy(nzg,i)       &
                  ,land%soil_water(nzg,i)*1.e3   &
                  ,slcpd(land%ntext_soil(nzg,i)) &
                  ,tempk, fracliq)
         fldval = tempk - 273.15
      endif
   endif

case(139) ! 'SFC_SSH'

   if (op%stagpt == 'S') then
      fldval = sea%surface_ssh(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%surface_ssh(i) * 1.e3
   endif

case(140) ! 'SENSFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_tsav(i) * cp / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_tsav(i) * cp / dt_leaf
   endif

case(141) ! 'VAPFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_rsav(i) / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_rsav(i) / dt_leaf
   endif

case(142) ! 'LATFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_rsav(i) * alvl / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_rsav(i) * alvl / dt_leaf
   endif

case(143) ! 'RSHORT_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rshort(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort(i)
   endif

case(144) ! 'RSHORT_DIFFUSE_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rshort_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort_diffuse(i)
   endif

case(145) ! 'RLONG_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlong(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong(i)
   endif

case(146) ! 'RLONGUP_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlongup(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlongup(i)
   endif

case(147) ! 'RLONG_ALBEDO_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlong_albedo(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong_albedo(i)
   endif

case(148) ! 'ALBEDO_BEAM_LS'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_beam(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_beam(i)
   endif

case(149) ! 'ALBEDO_DIFFUSE_LS'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_diffuse(i)
   endif

! FLUX CELLS

case(150) ! 'FCELL_ILSF'

   fldval = real(i)

case(151) ! 'FCELL_IWLS'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iwls)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iwls)
   endif

case(152) ! 'FCELL_IW'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iw)
   endif

case(153) ! 'FCELL_KW'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%kw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%kw)
   endif

case(154) ! FCELL_AREA'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%area
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%area
   endif

case(155) ! 'FCELL_ARFATM'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_atm
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_atm
   endif

case(156) ! 'FCELL_ARFLS'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_sfc
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_sfc
   endif

case(157) ! 'FCELL_SENS'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxt * cp
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxt * cp
   endif

case(158) ! 'FCELL_VAP'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr
   endif

case(159) ! 'FCELL_LAT'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr * alvl
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr * alvl
   endif

case(160) ! 'FCELL_AIRTEMPC'

   if (op%stagpt == 'S') then
      kw = seaflux(i)%kw
      iw = seaflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = tair(kw,iw) - 273.15
   elseif (op%stagpt == 'L') then
      kw = landflux(i)%kw
      iw = landflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = tair(kw,iw) - 273.15
   endif

case(161) ! 'FCELL_AIRTEMPK'

   if (op%stagpt == 'S') then
      kw = seaflux(i)%kw
      iw = seaflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = tair(kw,iw)
   elseif (op%stagpt == 'L') then
      kw = landflux(i)%kw
      iw = landflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = tair(kw,iw)
   endif

case(162) ! 'FCELL_CANTEMPC'

   if (op%stagpt == 'S') then
      iws = real(seaflux(i)%iwls)
      if (iparallel == 1) then
         iws = itabg_ws(iws)%iws_myrank
      endif

      fldval = sea%can_temp(iws) - 273.15
   elseif (op%stagpt == 'L') then
      iwl = real(landflux(i)%iwls)
      if (iparallel == 1) then
         iwl = itabg_wl(iwl)%iwl_myrank
      endif

      fldval = land%can_temp(iwl) - 273.15
   endif

! GRID GEOMETRY - 3D

case(163) ! 'ARU'
   if (.not. allocated(aru)) go to 1000
   fldval = aru(k,i)
case(164) ! 'ARV'
   if (.not. allocated(arv)) go to 1000
   fldval = arv(k,i)
case(165) ! 'ARW'
   fldval = arw(k,i)
case(166) ! 'VOLT'
   fldval = volt(k,i)
case(167) ! 'VOLU'
   if (.not. allocated(volui)) go to 1000
   fldval = 1. / volui(k,i)
case(168) ! 'VOLV'
   if (.not. allocated(volvi)) go to 1000
   fldval = 1. / volvi(k,i)
case(169) ! 'VOLW'
   fldval = 1. / volwi(k,i)

! GRID GEOMETRY - 2D

case(170) ! 'TOPM'
   fldval = topm(i)
case(171) ! 'TOPW'
   fldval = topw(i)
case(172) ! 'GLATM'
   fldval = glatm(i)
case(173) ! 'GLONM'
   fldval = glonm(i)
case(174) ! 'GLATU'
   if (.not. allocated(glatu)) go to 1000
   fldval = glatu(i)
case(175) ! 'GLONU'
   if (.not. allocated(glonu)) go to 1000
   fldval = glonu(i)
case(176) ! 'GLATV'
   if (.not. allocated(glatv)) go to 1000
   fldval = glatv(i)
case(177) ! 'GLONV'
   if (.not. allocated(glonv)) go to 1000
   fldval = glonv(i)
case(178) ! 'GLATW'
   fldval = glatw(i)
case(179) ! 'GLONW'
   fldval = glonw(i)
case(180) ! 'LPM'
   fldval = real(lpm(i))
case(181) ! 'LPU'
   if (.not. allocated(lpu)) go to 1000
   fldval = real(lpu(i))
case(182) ! 'LCU'
   if (.not. allocated(lcu)) go to 1000
   fldval = real(lcu(i))
case(183) ! 'LPV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(184) ! 'LCV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(185) ! 'LPW'
   fldval = real(lpw(i))
case(186) ! 'LSW'
   fldval = real(lsw(i))
case(187) ! 'XEM'
   fldval = xem(i)
case(188) ! 'YEM'
   fldval = yem(i)
case(189) ! 'ZEM'
   fldval = zem(i)
case(190) ! 'XEU'
   if (.not. allocated(xeu)) go to 1000
   fldval = xeu(i)
case(191) ! 'YEU'
   if (.not. allocated(yeu)) go to 1000
   fldval = yeu(i)
case(192) ! 'ZEU'
   if (.not. allocated(zeu)) go to 1000
   fldval = zeu(i)
case(193) ! 'XEV'
   if (.not. allocated(xev)) go to 1000
   fldval = xev(i)
case(194) ! 'YEV'
   if (.not. allocated(yev)) go to 1000
   fldval = yev(i)
case(195) ! 'ZEV'
   if (.not. allocated(zev)) go to 1000
   fldval = zev(i)
case(196) ! 'XEW'
   fldval = xew(i)
case(197) ! 'YEW'
   fldval = yew(i)
case(198) ! 'ZEW'
   fldval = zew(i)
case(199) ! 'DNU'
   fldval = dnu(i)
case(200) ! 'DNV'
   fldval = dnv(i)
case(201) ! 'ARM0'
   fldval = arm0(i)
case(202) ! 'ARW0'
   fldval = arw0(i)

! ITAB_M MEMBERS

case(203) ! 'ITAB_M_NPOLY'
   fldval = itab_m(i)%npoly
case(204) ! 'ITAB_M_IMGLOBE'
   fldval = itab_m(i)%imglobe
case(205) ! 'ITAB_M_MRLMR'
   fldval = itab_m(i)%mrlm
case(206) ! 'ITAB_M_MRLM_OR'
   fldval = itab_m(i)%mrlm_orig
case(207) ! 'ITAB_M_MROW'
   fldval = itab_m(i)%mrow
case(208) ! 'ITAB_M_MROWH'
   fldval = itab_m(i)%mrowh
case(209) ! 'ITAB_M_IU'
   fldval = itab_m(i)%iu(indp)
case(210) ! 'ITAB_M_IV'
   fldval = itab_m(i)%iv(indp)
case(211) ! 'ITAB_M_IW'
   fldval = itab_m(i)%iw(indp)
case(212) ! 'ITAB_M_FMW'
   go to 1000 ! currently not used
   ! fldval = itab_m(i)%fmw(indp)

! ITAB_U MEMBERS

case(213) ! 'ITAB_U_IUP'
   fldval = itab_u(i)%iup
case(214) ! 'ITAB_U_IRANK'
   fldval = itab_u(i)%irank
!   fldval = itabg_u(i)%irank
case(215) ! 'ITAB_U_IUGLOBE'
   fldval = itab_u(i)%iuglobe
case(216) ! 'ITAB_U_MRLU'
   fldval = itab_u(i)%mrlu
case(217) ! 'ITAB_U_IM'
   fldval = itab_u(i)%im(indp)
case(218) ! 'ITAB_U_IU'
   fldval = itab_u(i)%iu(indp)
case(219) ! 'ITAB_U_IW'
   fldval = itab_u(i)%iw(indp)
case(220) ! 'ITAB_U_DIRU'
   fldval = itab_u(i)%diru(indp)
case(221) ! 'ITAB_U_FUU'
   fldval = itab_u(i)%fuu(indp)
case(222) ! 'ITAB_U_FUW'
   fldval = itab_u(i)%fuw(indp)
case(223) ! 'ITAB_U_TUU'
   fldval = itab_u(i)%tuu(indp)
case(224) ! 'ITAB_U_GUW'
   fldval = itab_u(i)%guw(indp)
case(225) ! 'ITAB_U_PGC12'
   fldval = itab_u(i)%pgc12
case(226) ! 'ITAB_U_PGC12b'
   fldval = itab_u(i)%pgc12b
case(227) ! 'ITAB_U_PGC12c'
   fldval = itab_u(i)%pgc12c
case(228) ! 'ITAB_U_PGC12d'
   fldval = itab_u(i)%pgc12d
case(229) ! 'ITAB_U_PGC45'
   fldval = itab_u(i)%pgc45
case(230) ! 'ITAB_U_PGC45b'
   fldval = itab_u(i)%pgc45b
case(231) ! 'ITAB_U_PGC63'
   fldval = itab_u(i)%pgc63
case(232) ! 'ITAB_U_PGC63c'
   fldval = itab_u(i)%pgc63c
case(233) ! 'ITAB_U_VXU_U'
   fldval = itab_u(i)%vxu_u(indp)
case(234) ! 'ITAB_U_VYU_U'
   fldval = itab_u(i)%vyu_u(indp)
case(235) ! 'ITAB_U_VXW_U''
   fldval = itab_u(i)%vxw_u(indp)
case(236) ! 'ITAB_U_VYW_U'
   fldval = itab_u(i)%vyw_u(indp)
case(237) ! 'ITAB_U_CROSSMM'
   fldval = itab_u(i)%crossmm
case(238) ! 'ITAB_U_CROSSWW'
   fldval = itab_u(i)%crossww

! ITAB_V MEMBERS

case(239) ! 'ITAB_V_IVP'
   fldval = itab_v(i)%ivp
case(240) ! 'ITAB_V_IRANK'
   fldval = itab_v(i)%irank
!   fldval = itabg_v(i)%irank
case(241) ! 'ITAB_V_IVGLOBE'
   fldval = itab_v(i)%ivglobe
case(242) ! 'ITAB_V_MRLV'
   fldval = itab_v(i)%mrlv
case(243) ! 'ITAB_V_IM'
   fldval = itab_v(i)%im(indp)
case(244) ! 'ITAB_V_IV'
   fldval = itab_v(i)%iv(indp)
case(245) ! 'ITAB_V_IW'
   fldval = itab_v(i)%iw(indp)
case(246) ! 'ITAB_V_FVV'
   go to 1000 ! currently not used
   ! fldval = itab_v(i)%fvv(indp)
case(247) ! 'ITAB_V_FVW'
   go to 1000 ! currently not used
   ! fldval = itab_v(i)%fvw(indp)
case(248) ! 'ITAB_V_FUV'
   go to 1000 ! currently not used
   ! fldval = itab_v(i)%fuv(indp)
case(249) ! 'ITAB_V_FARW'
   fldval = itab_v(i)%farw(indp)

! ITAB_W MEMBERS

case(250) ! 'ITAB_W_NPOLY'
   fldval = itab_w(i)%npoly
case(251) ! 'ITAB_W_IWP'
   fldval = itab_w(i)%iwp
case(252) ! 'ITAB_W_IRANK'
   fldval = itab_w(i)%irank
!   fldval = itabg_w(i)%irank
case(253) ! 'ITAB_W_IWGLOBE'
   fldval = itab_w(i)%iwglobe
case(254) ! 'ITAB_W_MRLW'
   fldval = itab_w(i)%mrlw
case(255) ! 'ITAB_W_MROW'
   fldval = itab_w(i)%mrow
case(256) ! 'ITAB_W_MROWH'
   fldval = itab_w(i)%mrowh
case(257) ! 'ITAB_W_IM'
   fldval = itab_w(i)%im(indp)
case(258) ! 'ITAB_W_IU'
   fldval = itab_w(i)%iu(indp)
case(259) ! 'ITAB_W_IV'
   fldval = itab_w(i)%iv(indp)
case(260) ! 'ITAB_W_IW'
   fldval = itab_w(i)%iw(indp)
case(261) ! 'ITAB_W_DIRU'
   fldval = itab_w(i)%diru(indp)
case(262) ! 'ITAB_W_DIRV'
   fldval = itab_w(i)%dirv(indp)
case(263) ! 'ITAB_W_FWV'
   fldval = itab_w(i)%fwv(indp)
case(264) ! 'ITAB_W_FWW'
   fldval = itab_w(i)%fww(indp)
case(265) ! 'ITAB_W_FWU'
   fldval = itab_w(i)%fwu(indp)
case(266) ! 'ITAB_W_VXU'
   fldval = itab_w(i)%vxu(indp)
case(267) ! 'ITAB_W_VYU'
   fldval = itab_w(i)%vyu(indp)
case(268) ! 'ITAB_W_VZU'
   fldval = itab_w(i)%vzu(indp)
case(269) ! 'ITAB_W_VXW'
   fldval = itab_w(i)%vxw
case(270) ! 'ITAB_W_VYW'
   fldval = itab_w(i)%vyw
case(271) ! 'ITAB_W_VZW'
   fldval = itab_w(i)%vzw
case(272) ! 'ITAB_W_VXU_W'
   fldval = itab_w(i)%vxu_w(indp)
case(273) ! 'ITAB_W_VYU_W'
   fldval = itab_w(i)%vyu_w(indp)
case(274) ! 'ITAB_W_FARM'
   fldval = itab_w(i)%farm(indp)
case(275) ! 'ITAB_W_FARV'
   fldval = itab_w(i)%farv(indp)
case(276) ! 'ITAB_W_IWNUD'
   fldval = itab_w(i)%iwnud(indp)
case(277) ! 'ITAB_W_FNUD'
   fldval = itab_w(i)%fnud(indp)

! Time-averaged fields

case(278) ! 'PRESS_MAVG' [indp = 1:nz_avg is the selected vertical level]

   if (.not. allocated(press_mavg)) go to 1000

   fldval = press_mavg(indp,i) * .01

case(279) !  'RHO_MAVG'

   if (.not. allocated(rho_mavg)) go to 1000

   fldval = rho_mavg(indp,i)

case(280) !  'TEMPK_MAVG'

   if (.not. allocated(tempk_mavg)) go to 1000

   fldval = tempk_mavg(indp,i)

case(281) !  'SH_V_MAVG'

   if (.not. allocated(sh_v_mavg)) go to 1000

   fldval = sh_v_mavg(indp,i) * 1.e3

case(282) !  'SH_W_MAVG'

   if (.not. allocated(sh_w_mavg)) go to 1000

   fldval = sh_w_mavg(indp,i) * 1.e3

case(283) !  'WC_MAVG'

   if (.not. allocated(wc_mavg)) go to 1000

   fldval = wc_mavg(indp,i)

case(284:285) !  'ZONAL_WINDW_MAVG','MERID_WINDW_MAVG'

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

case(286) ! 'RSHORT_MAVG'

   if (.not. allocated(rshort_mavg)) go to 1000

   fldval = rshort_mavg(i)

case(287) ! 'RSHORT_TOP_MAVG'

   if (.not. allocated(rshort_top_mavg)) go to 1000

   fldval = rshort_top_mavg(i)

case(288) ! 'RSHORTUP_MAVG'

   if (.not. allocated(rshortup_mavg)) go to 1000

   fldval = rshortup_mavg(i)

case(289) ! 'RSHORTUP_TOP_MAVG'

   if (.not. allocated(rshortup_top_mavg)) go to 1000

   fldval = rshortup_top_mavg(i)

case(290) ! 'RLONG_MAVG'

   if (.not. allocated(rlong_mavg)) go to 1000

   fldval = rlong_mavg(i)

case(291) ! 'RLONGUP_MAVG'

   if (.not. allocated(rlongup_mavg)) go to 1000

   fldval = rlongup_mavg(i)

case(292) ! 'RLONGUP_TOP_MAVG'

   if (.not. allocated(rlongup_top_mavg)) go to 1000

   fldval = rlongup_top_mavg(i)

case(293) ! 'LATFLUX_MAVG'

   if (.not. allocated(latflux_mavg)) go to 1000

   fldval = latflux_mavg(i)

case(294) ! 'SENSFLUX_MAVG'

   if (.not. allocated(sensflux_mavg)) go to 1000

   fldval = sensflux_mavg(i)

case(295) ! 'WINDSPEED_MAVG'

   if (.not. allocated(windspeed_mavg)) go to 1000

   fldval = windspeed_mavg(i)

case(296) ! 'ACCPMIC_MTOT'

   if (.not. allocated(accpmic_mtot)) go to 1000

   fldval = accpmic_mtot(i)

case(297) ! 'ACCPCON_MTOT'

   if (.not. allocated(accpcon_mtot)) go to 1000

   fldval = accpcon_mtot(i)

case(298) ! 'ACCPBOTH_MTOT'

   if (.not. allocated(accpmic_mtot)) go to 1000
   if (.not. allocated(accpcon_mtot)) go to 1000

   fldval = accpmic_mtot(i) + accpcon_mtot(i)

case(299) ! 'PRESS_DAVG'

   if (.not. allocated(press_davg)) go to 1000

   fldval = press_davg(i) * .01

case(300:301) ! 'ZONAL_WINDW_DAVG','MERID_WINDW_DAVG'

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

case(302) ! 'RSHORT_DAVG'

   if (.not. allocated(rshort_davg)) go to 1000

   fldval = rshort_davg(i)

case(303) ! 'TEMPK_DAVG'

   if (.not. allocated(tempk_davg)) go to 1000

   fldval = tempk_davg(i)

case(304) ! 'TEMPK_DMIN'

   if (.not. allocated(tempk_dmin)) go to 1000

   fldval = tempk_dmin(i)

case(305) ! 'TEMPK_DMAX'

   if (.not. allocated(tempk_dmax)) go to 1000

   fldval = tempk_dmax(i)

case(306) ! 'ACCPMIC_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000

   fldval = accpmic_dtot(i)

case(307) ! 'ACCPCON_DTOT'

   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpcon_dtot(i)

case(308) ! 'ACCPBOTH_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000
   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpmic_dtot(i) + accpcon_dtot(i)

case(309) ! 'PRESS_UL_DAVG'

   if (.not. allocated(press_ul_davg)) go to 1000

   fldval = press_ul_davg(i) * .01

case(310:311) ! 'ZONAL_WINDW_UL_DAVG','MERID_WINDW_UL_DAVG'

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

case(312) ! 'CANTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(canltempk_davg)) go to 1000
      fldval = canltempk_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(canstempk_davg)) go to 1000
      fldval = canstempk_davg(i)
   endif

case(313) ! 'CANTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(canltempk_dmin)) go to 1000
      fldval = canltempk_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(canstempk_dmin)) go to 1000
      fldval = canstempk_dmin(i)
   endif

case(314) ! 'CANTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(canltempk_dmax)) go to 1000
      fldval = canltempk_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(canstempk_dmax)) go to 1000
      fldval = canstempk_dmax(i)
   endif

case(315) ! 'VEGTEMPK_DAVG'

   if (.not. allocated(vegtempk_davg)) go to 1000

   fldval = vegtempk_davg(i)

case(316) ! 'VEGTEMPK_DMIN'

   if (.not. allocated(vegtempk_dmin)) go to 1000

   fldval = vegtempk_dmin(i)

case(317) ! 'VEGTEMPK_DMAX'

   if (.not. allocated(vegtempk_dmax)) go to 1000

   fldval = vegtempk_dmax(i)

case(318) ! 'SOILTEMPK_DAVG'

   if (.not. allocated(soiltempk_davg)) go to 1000

   fldval = soiltempk_davg(i)

case(319) ! 'SOILTEMPK_DMIN'

   if (.not. allocated(soiltempk_dmin)) go to 1000

   fldval = soiltempk_dmin(i)

case(320) ! 'SOILTEMPK_DMAX'

   if (.not. allocated(soiltempk_dmax)) go to 1000

   fldval = soiltempk_dmax(i)

case(321) ! 'FCELL_AIRTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(tempk_lf_davg)) go to 1000
      fldval = tempk_lf_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(tempk_sf_davg)) go to 1000
      fldval = tempk_sf_davg(i)
   endif

case(322) ! 'FCELL_AIRTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(tempk_lf_dmin)) go to 1000
      fldval = tempk_lf_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(tempk_sf_dmin)) go to 1000
      fldval = tempk_sf_dmin(i)
   endif

case(323) ! 'FCELL_AIRTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(tempk_lf_dmax)) go to 1000
      fldval = tempk_lf_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(tempk_sf_dmax)) go to 1000
      fldval = tempk_sf_dmax(i)
   endif

case(324) ! 'FCELL_CANTEMPK_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_lf_davg)) go to 1000
      fldval = cantempk_lf_davg(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_sf_davg)) go to 1000
      fldval = cantempk_sf_davg(i)
   endif

case(325) ! 'FCELL_CANTEMPK_DMIN'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_lf_dmin)) go to 1000
      fldval = cantempk_lf_dmin(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_sf_dmin)) go to 1000
      fldval = cantempk_sf_dmin(i)
   endif

case(326) ! 'FCELL_CANTEMPK_DMAX'

   if (op%stagpt == 'L') then
      if (.not. allocated(cantempk_lf_dmax)) go to 1000
      fldval = cantempk_lf_dmax(i)
   elseif (op%stagpt == 'S') then
      if (.not. allocated(cantempk_sf_dmax)) go to 1000
      fldval = cantempk_sf_dmax(i)
   endif

case(327) ! 'FCELL_SENSFLUX_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxt_lf_davg)) go to 1000
      fldval = sfluxt_lf_davg(i) * cp
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxt_sf_davg)) go to 1000
      fldval = sfluxt_sf_davg(i) * cp
   endif

case(328) ! 'FCELL_LATFLUX_DAVG'

   if (op%stagpt == 'L') then
      if (.not. allocated(sfluxr_lf_davg)) go to 1000
      fldval = sfluxr_lf_davg(i) * alvl
   elseif (op%stagpt == 'S') then
      if (.not. allocated(sfluxr_sf_davg)) go to 1000
      fldval = sfluxr_sf_davg(i) * alvl
   endif

case(329) ! 'FTOA_SENSFLUX_DAVG'

   if (infotyp == 'UNITS') then
      aux(:) = 0.
      do isf = 2,mseaflux
         iw = seaflux(isf)%iw
         aux(iw) = aux(iw) + seaflux(isf)%arf_atm * sfluxt_sf_davg(isf) * cp
      enddo
      do ilf = 2,mlandflux
         iw = landflux(ilf)%iw
         aux(iw) = aux(iw) + landflux(ilf)%arf_atm * sfluxt_lf_davg(ilf) * cp
      enddo
   endif

   fldval = aux(i)

case(330) ! 'FTOA_LATFLUX_DAVG'

   if (infotyp == 'UNITS') then
      aux(:) = 0.
      do isf = 2,mseaflux
         iw = seaflux(isf)%iw
         aux(iw) = aux(iw) + seaflux(isf)%arf_atm * sfluxr_sf_davg(isf) * alvl
      enddo
      do ilf = 2,mlandflux
         iw = landflux(ilf)%iw
         aux(iw) = aux(iw) + landflux(ilf)%arf_atm * sfluxr_lf_davg(ilf) * alvl
      enddo
   endif

   fldval = aux(i)

case(331) ! 'PRESS_AVG24' [indp = 1:24 is the hour of day to be plotted]

   if (.not. allocated(press_avg24)) go to 1000

   fldval = press_avg24(indp,i)

case(332) ! 'RSHORT_AVG'

   if (.not. allocated(rshort_avg)) go to 1000

   fldval = rshort_avg(i)

case(333) ! 'RSHORTUP_AVG'

   if (.not. allocated(rshortup_avg)) go to 1000

   fldval = rshortup_avg(i)

case(334) ! 'RLONG_AVG'

   if (.not. allocated(rlong_avg)) go to 1000

   fldval = rlong_avg(i)

case(335) ! 'RLONGUP_AVG'

   if (.not. allocated(rlongup_avg)) go to 1000

   fldval = rlongup_avg(i)

case(336) ! 'RSHORT_TOP_AVG'

   if (.not. allocated(rshort_top_avg)) go to 1000

   fldval = rshort_top_avg(i)

case(337) ! 'RSHORTUP_TOP_AVG'

   if (.not. allocated(rshortup_top_avg)) go to 1000

   fldval = rshortup_top_avg(i)

case(338) ! 'RLONGUP_TOP_AVG'

   if (.not. allocated(rlongup_top_avg)) go to 1000

   fldval = rlongup_top_avg(i)

case(339) ! 'SENSFLUX_AVG'

   if (.not. allocated(sflux_t_avg)) go to 1000

   fldval = sflux_t_avg(i) * cp

case(340) ! 'LATFLUX_AVG'

   if (.not. allocated(sflux_r_avg)) go to 1000

   fldval = sflux_r_avg(i) * alvl

case(341) ! 'VAPFLUX_AVG'

   if (.not. allocated(sflux_r_avg)) go to 1000

   fldval = sflux_r_avg(i)

! Miscellaneous and new additions

case(342) ! 'RHO_OBS'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1))

case(343) ! 'THETA_OBS'

   if (.not. allocated(theta_obs)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i)

case(344) ! 'SHW_OBS'

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3

case(345) ! 'UZONAL_OBS'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1))

case(346) ! 'UMERID_OBS'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1))

case(347) ! 'RHO_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_sim(k,itab_w(i)%iwnud(1))

case(348) ! 'THETA_SIM'

   if (.not. allocated(theta_sim)) go to 1000

!   fldval = theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_sim(k,i)

case(349) ! 'SHW_SIM'

   if (.not. allocated(shw_sim)) go to 1000

   fldval = shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(350) ! 'UZONAL_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_sim(k,itab_w(i)%iwnud(1))

case(351) ! 'UMERID_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_sim(k,itab_w(i)%iwnud(1))

case(352) ! 'RHO_OBS_SIM'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1)) &
          - rho_sim(k,itab_w(i)%iwnud(1))

case(353) ! 'THETA_OBS_SIM'

   if (.not. allocated(theta_obs)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1)) &
!          - theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i) &
          - theta_sim(k,i)

case(354) ! 'SHW_OBS_SIM'

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3 &
          - shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(355) ! 'UZONAL_OBS_SIM'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1)) &
          - uzonal_sim(k,itab_w(i)%iwnud(1))

case(356) ! 'UMERID_OBS_SIM'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1)) &
          - umerid_sim(k,itab_w(i)%iwnud(1))

case(357) ! 'VXE'

   if (.not. allocated(vxe)) go to 1000

   fldval = wtbot * vxe(k  ,i) &
          + wttop * vxe(k+1,i)

case(358) ! 'VYE'

   if (.not. allocated(vye)) go to 1000

   fldval = wtbot * vye(k  ,i) &
          + wttop * vye(k+1,i)

case(359) ! 'VZE'

   if (.not. allocated(vze)) go to 1000

   fldval = wtbot * vze(k  ,i) &
          + wttop * vze(k+1,i)

case(360) ! 'PBLH'

   if (.not. allocated(pblh)) go to 1000

   fldval = pblh(i)

case(361) ! 'HKM'

   if (.not. allocated(hkm)) go to 1000

   fldval = wtbot * hkm(k  ,i) / rho(k  ,i) &
          + wttop * hkm(k+1,i) / rho(k+1,i)

case(362) ! 'HEAD0'

   fldval = land%head0(i)

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
