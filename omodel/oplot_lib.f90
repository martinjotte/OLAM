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
                       thil, theta, sh_w, sh_v
use mem_cuparm,  only: conprr, aconpr

use mem_grid,    only: mza, mua, mva, mwa, lpm, lpu, lcu, lpv, lcv, lpw, lsw, &
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
                      
use mem_radiate, only: fthrd, rshort, rlong, rlongup, albedt, &
                       rshort_top, rshortup_top, rlongup_top
use mem_addsc,   only: addsc
use mem_tend,    only: umt, vmt, wmt
use mem_turb,    only: vkm, vkm_sfc, sflux_w, sflux_t, sflux_r
use mem_nudge,   only: rho_obs, theta_obs, shw_obs, uzonal_obs, umerid_obs, &
                       rho_sim, theta_sim, shw_sim, uzonal_sim, umerid_sim

use misc_coms,   only: io6, pr01d, dn01d, th01d, time8, iparallel, meshtype, naddsc
use oplot_coms,  only: op
use consts_coms, only: p00, rocp, erad, piu180, cp, alvl, grav, omega2
use leaf_coms,   only: slcpd, nzg, slmsts, slz, mwl, dt_leaf
use mem_sflux,   only: landflux, seaflux
use sea_coms,    only: mws
use mem_timeavg, only: rshort_avg, rshortup_avg, rlong_avg, rlongup_avg, &
                       rshort_top_avg, rshortup_top_avg, rlongup_top_avg, &
                       sflux_t_avg, sflux_r_avg

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
real :: vx,vy,vz,raxis,u,v,vxe,vye,vze,ucint,vcint,farv2,rpolyi
real :: tempk,fracliq
real :: contrib

integer :: iw1,iw2,iwl,iws
integer :: npoly,j
integer :: lenstr, ic, ifield
integer, save :: indp, icase

real :: ucc,vcc
real :: ucc_init, vcc_init, vx_init, vy_init, vz_init, u_init, v_init

! Stored initial values for perturbation calculations

integer, save :: ncall = 0
real(kind=8), save, allocatable :: press_init(:,:)
real(kind=8), save, allocatable :: rho_init(:,:)
real,         save, allocatable :: theta_init(:,:)
real,         save, allocatable :: uc_init(:,:)
real,         save, allocatable :: vc_init(:,:)
real,         save, allocatable :: addsc1_init(:,:)

integer, parameter :: nfields = 296
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

data fldlib(1:4, 68:100)/ &
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
 'PCPRTOT'       ,'T2' ,'MICROPHYS PCP RATE',' (kg m:S2:-2   h:S2:-1  )'    ,& !  88
 'CONPRR'        ,'T2' ,'CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )'         ,& !  89
 'TOTPRR'        ,'T2' ,'TOTAL PCP RATE',' (kg m:S2:-2   h:S2:-1  )'        ,& !  90
 'ACCPD'         ,'T2' ,'ACCUM DRIZZLE',' (kg m:S2:-2  )'                   ,& !  91
 'ACCPR'         ,'T2' ,'ACCUM RAIN',' (kg m:S2:-2  )'                      ,& !  92
 'ACCPP'         ,'T2' ,'ACCUM PRIST ICE',' (kg m:S2:-2  )'                 ,& !  93
 'ACCPS'         ,'T2' ,'ACCUM SNOW',' (kg m:S2:-2  )'                      ,& !  94
 'ACCPA'         ,'T2' ,'ACCUM AGGREGATES',' (kg m:S2:-2  )'                ,& !  95
 'ACCPG'         ,'T2' ,'ACCUM GRAUPEL',' (kg m:S2:-2  )'                   ,& !  96
 'ACCPH'         ,'T2' ,'ACCUM HAIL',' (kg m:S2:-2  )'                      ,& !  97
 'ACCPTOT'       ,'T2' ,'ACCUM PRECIP',' (kg m:S2:-2  )'                    ,& !  98
 'ACONPR'        ,'T2' ,'ACCUM CONV PCP',' (kg m:S2:-2  )'                  ,& !  99
 'ATOTPR'        ,'T2' ,'TOTAL ACCUM PCP',' (kg m:S2:-2  )'                  / ! 100

! LAND_CELLS - 3D

data fldlib(1:4,101:123)/ &
 'SOIL_TEXT'     ,'L3G','SOIL TEXTURAL CLASS',' ( )'                        ,& ! 101
 'SOIL_ENERGY'   ,'L3G','SOIL ENERGY',' (J cm:S2:-3  )'                     ,& ! 102
 'SOIL_TEMPC'    ,'L3G','SOIL TEMP',' (C)'                                  ,& ! 103
 'SOIL_FRACLIQ'  ,'L3G','LIQUID FRACTION OF SOIL WATER',' ( )'              ,& ! 104
 'SOIL_WATER'    ,'L3G','SOIL WATER CONTENT',' ( )'                         ,& ! 105
 'SFWAT_MASS'    ,'L3S','SFCWATER MASS',' (kg m:S2:-2  )'                   ,& ! 106
 'SFWAT_ENERGY'  ,'L3S','SFCWATER ENERGY',' (J g:S2:-1  )'                  ,& ! 107
 'SFWAT_TEMPC'   ,'L3S','SFCWATER TEMP',' (C)'                              ,& ! 108
 'SFWAT_FRACLIQ' ,'L3S','SFCWATER LIQUID FRACTION',' ( )'                   ,& ! 109
 'SFWAT_DEPTH'   ,'L3S','SFCWATER DEPTH',' (m)'                             ,& ! 110

! LAND_CELLS - 2D

 'NLEV_SFWAT'    ,'L2' ,'NUMBER OF SFCWATER LAYERS',' ( )'                  ,& ! 111
 'VEG_NDVIC'     ,'L2' ,'VEGETATION NDVI',' ( )'                            ,& ! 112
 'VEG_TEMPC'     ,'L2' ,'VEGETATION TEMP',' (C)'                            ,& ! 113
 'VEG_WATER'     ,'L2' ,'VEGETATION SFC WATER ',' (kg m:S2:-2  )'           ,& ! 114
 'STOM_RESIST'   ,'L2' ,'STOMATAL RESISTANCE',' (s m:S2:-1  )'              ,& ! 115
 'GROUND_SHV'    ,'L2' ,'EQUIL VAP SPEC DENSITY OF SOIL',' (g kg:S2:-1  )'  ,& ! 116
 'SOIL_DEPTH'    ,'L2' ,'SOIL DEPTH',' (m)'                                 ,& ! 117

! SEA_CELLS - 2D

 'SEATP'         ,'S2' ,'SEA SFC TEMP (PAST DATA)',' (K)'                   ,& ! 118
 'SEATF'         ,'S2' ,'SEA SFC TEMP (FUTURE DATA)',' (K)'                 ,& ! 119
 'SEATC'         ,'S2' ,'SEA SFC TEMP (CURRENT)',' (K)'                     ,& ! 120
 'SEAICEP'       ,'S2' ,'SEAICE FRACTION (PAST DATA)',' ( )'                ,& ! 121
 'SEAICEF'       ,'S2' ,'SEAICE FRACTION (FUTURE DATA)',' ( )'              ,& ! 122
 'SEAICEC'       ,'S2' ,'SEAICE FRACTION (CURRENT)',' ( )'                   / ! 123

! LAND AND SEA CELLS - 2D

data fldlib(1:4,124:140)/ &
 'LEAF_CLASS'    ,'B2' ,'LEAF CLASS',' ( )'                                 ,& ! 124
 'AREA'          ,'B2' ,'LAND/SEA CELL AREA',' (m:S2:2  )'                  ,& ! 125
 'ROUGH'         ,'B2' ,'NET ROUGHNESS HEIGHT',' (m)'                       ,& ! 126
 'CAN_TEMPC'     ,'B2' ,'CANOPY AIR TEMP',' (C)'                            ,& ! 127
 'CAN_SHV'       ,'B2' ,'CANOPY VAPOR SPEC DENSITY',' (g kg:S2:-1  )'       ,& ! 128
 'SFC_TEMPC'     ,'B2' ,'SOIL/SFCWATER/SEA SFC TEMP',' (C)'                 ,& ! 129
 'SFC_SSH'       ,'B2' ,'LAND/SEA SFC SAT VAP SPEC DENS',' (g kg:S2:-1  )'  ,& ! 130

 'SENSFLUX_LS'      ,'B2','L/S CAN TOP SENS HEAT FLX',' (W m:S2:-2  )'      ,& ! 131
 'VAPFLUX_LS'       ,'B2','L/S CAN TOP VAP FLX',' (kg m:S2:-2   s:S2:-1  )' ,& ! 132
 'LATFLUX_LS'       ,'B2','L/S CAN TOP LAT HEAT FLX',' (W m:S2:-2  )'       ,& ! 133
 'RSHORT_LS'        ,'B2','L/S CAN TOP DOWN SW FLX',' (W m:S2:-2  )'        ,& ! 134
 'RSHORT_DIFFUSE_LS','B2','L/S CAN TOP DOWN DIFFUSE SW FLX',' (W m:S2:-2  )',& ! 135
 'RLONG_LS'         ,'B2','L/S CAN TOP DOWN LW FLX',' (W m:S2:-2  )'        ,& ! 136
 'RLONGUP_LS'       ,'B2','L/S CAN TOP UP LW FLX',' (W m:S2:-2  )'          ,& ! 137
 'RLONG_ALBEDO_LS'  ,'B2','L/S NET SFC LW ALBEDO',' ( )'                    ,& ! 138
 'ALBEDO_BEAM_LS'   ,'B2','L/S NET SFC BEAM ALBEDO',' ( )'                  ,& ! 139
 'ALBEDO_DIFFUSE_LS','B2','L/S NET SFC DIFFUSE ALBEDO',' ( )'                / ! 140

! FLUX CELLS - 2D

data fldlib(1:4,141:153)/  &
 'FCELL_ILSF'    ,'F2' ,'FLUXCELL INDEX',' '                                ,& ! 141
 'FCELL_IWLS'    ,'F2' ,'FLUXCELL LAND/SEA INDEX',' '                       ,& ! 142
 'FCELL_IW'      ,'F2' ,'FLUXCELL ATM IW INDEX',' '                         ,& ! 143
 'FCELL_KW'      ,'F2' ,'FLUXCELL ATM KW INDEX',' '                         ,& ! 144
 'FCELL_AREA'    ,'F2' ,'FLUXCELL AREA',' (m:S2:2  )'                       ,& ! 145
 'FCELL_ARFATM'  ,'F2' ,'FLUXCELL ATM AREA FRACTION',' '                    ,& ! 146
 'FCELL_ARFLS'   ,'F2' ,'FLUXCELL LAND/SEA AREA FRACTION',' '               ,& ! 147
 'FCELL_SENS'    ,'F2' ,'FLUXCELL SENS HEAT FLUX',' (W m:S2:-2  )'          ,& ! 148
 'FCELL_VAP'     ,'F2' ,'FLUXCELL VAP FLUX',' (kg m:S2:-2   s:S2:-1  )'     ,& ! 149
 'FCELL_LAT'     ,'F2' ,'FLUXCELL LAT HEAT FLUX',' (W m:S2:-2  )'           ,& ! 150
 'FCELL_AIRTEMPC','F2' ,'FLUXCELL AIR TEMP',' (C)'                          ,& ! 151
 'FCELL_AIRTEMPK','F2' ,'FLUXCELL AIR TEMP',' (K)'                          ,& ! 152
 'FCELL_CANTEMPC','F2' ,'FLUXCELL CAN TEMP',' (C)'                           / ! 153

! GRID GEOMETRY - 3D

data fldlib(1:4,154:160)/ &
 'ARU'           ,'V3' ,'AREA OF GRID CELL U-FACE',' (m:S2:2  )'            ,& ! 154
 'ARV'           ,'V3' ,'AREA OF GRID CELL V-FACE',' (m:S2:2  )'            ,& ! 155
 'ARW'           ,'W3' ,'AREA OF GRID CELL W-FACE',' (m:S2:2  )'            ,& ! 156
 'VOLT'          ,'T3' ,'GRID T-CELL VOLUME',' (m:S2:3  )'                  ,& ! 157
 'VOLU'          ,'V3' ,'GRID U-CELL VOLUME',' (m:S2:3  )'                  ,& ! 158
 'VOLV'          ,'V3' ,'GRID V-CELL VOLUME',' (m:S2:3  )'                  ,& ! 159
 'VOLW'          ,'W3' ,'GRID W-CELL VOLUME',' (m:S2:3  )'                   / ! 160

! GRID GEOMETRY - 2D

data fldlib(1:4,161:193)/ &
 'TOPM'          ,'M2' ,'TOPOGRAPHY HEIGHT',' (m)'                          ,& ! 161
 'TOPW'          ,'W2' ,'TOPOGRAPHY HEIGHT AT W',' (m)'                     ,& ! 162
 'GLATM'         ,'M2' ,'LATITUDE AT M',' (deg)'                            ,& ! 163
 'GLONM'         ,'M2' ,'LONGITUDE AT M',' (deg)'                           ,& ! 164
 'GLATU'         ,'V2' ,'LATITUDE AT U',' (deg)'                            ,& ! 165
 'GLONU'         ,'V2' ,'LONGITUDE AT U',' (deg)'                           ,& ! 166
 'GLATV'         ,'V2' ,'LATITUDE AT V',' (deg)'                            ,& ! 167
 'GLONV'         ,'V2' ,'LONGITUDE AT V',' (deg)'                           ,& ! 168
 'GLATW'         ,'T2' ,'LATITUDE',' (deg)'                                 ,& ! 169
 'GLONW'         ,'T2' ,'LONGITUDE',' (deg)'                                ,& ! 170
 'LPM'           ,'M2' ,'LOWEST PREDICTED M LEVEL',' ( )'                   ,& ! 171
 'LPU'           ,'V2' ,'LOWEST PREDICTED U LEVEL',' ( )'                   ,& ! 172
 'LCU'           ,'V2' ,'LOWEST ACTIVE U CONTROL VOL',' ( )'                ,& ! 173
 'LPV'           ,'V2' ,'LOWEST PREDICTED V LEVEL',' ( )'                   ,& ! 174
 'LCV'           ,'V2' ,'LOWEST ACTIVE V CONTROL VOL',' ( )'                ,& ! 175
 'LPW'           ,'W2' ,'LOWEST PREDICTED W LEVEL',' ( )'                   ,& ! 176
 'LSW'           ,'W2' ,'NUMBER OF SFC W LEVELS',' ( )'                     ,& ! 177
 'XEM'           ,'M2' ,'EARTH-X COORD OF M POINT',' ( )'                   ,& ! 178
 'YEM'           ,'M2' ,'EARTH-Y COORD OF M POINT',' ( )'                   ,& ! 179
 'ZEM'           ,'M2' ,'EARTH-Z COORD OF M POINT',' ( )'                   ,& ! 180
 'XEU'           ,'V2' ,'EARTH-X COORD OF U POINT',' ( )'                   ,& ! 181
 'YEU'           ,'V2' ,'EARTH-Y COORD OF U POINT',' ( )'                   ,& ! 182
 'ZEU'           ,'V2' ,'EARTH-Z COORD OF U POINT',' ( )'                   ,& ! 183
 'XEV'           ,'V2' ,'EARTH-X COORD OF V POINT',' ( )'                   ,& ! 184
 'YEV'           ,'V2' ,'EARTH-Y COORD OF V POINT',' ( )'                   ,& ! 185
 'ZEV'           ,'V2' ,'EARTH-Z COORD OF V POINT',' ( )'                   ,& ! 186
 'XEW'           ,'W2' ,'EARTH-X COORD OF W POINT',' ( )'                   ,& ! 187
 'YEW'           ,'W2' ,'EARTH-Y COORD OF W POINT',' ( )'                   ,& ! 188
 'ZEW'           ,'W2' ,'EARTH-Z COORD OF W POINT',' ( )'                   ,& ! 189
 'DNU'           ,'V2' ,'DNU',' (m)'                                        ,& ! 190
 'DNV'           ,'V2' ,'DNV',' (m)'                                        ,& ! 191
 'ARM0'          ,'W2' ,'SFC AREA OF M CELL',' (m:S2:2  )'                  ,& ! 192
 'ARW0'          ,'W2' ,'SFC AREA OF W CELL',' (m:S2:2  )'                   / ! 193

! ITAB_M MEMBERS - 2D

data fldlib(1:4,194:229)/ &
 'ITAB_M_NPOLY'  ,'M2' ,'ITAB_M_NPOLY',' ( )'                               ,& ! 194
 'ITAB_M_IMGLOBE','M2' ,'ITAB_M_IMGLOBE',' ( )'                             ,& ! 195
 'ITAB_M_MRLM'   ,'V2' ,'ITAB_M_MRLM',' ( )'                                ,& ! 196
 'ITAB_M_MRLM_OR','V2' ,'ITAB_M_MRLM_ORIG',' ( )'                           ,& ! 197
 'ITAB_M_MROW'   ,'M2' ,'ITAB_M_MROW',' ( )'                                ,& ! 198
 'ITAB_M_MROWH'  ,'M2' ,'ITAB_M_MROWH',' ( )'                               ,& ! 199
 'ITAB_M_IU'     ,'M2' ,'ITAB_M_IU',' ( )'                                  ,& ! 200
 'ITAB_M_IV'     ,'M2' ,'ITAB_M_IV',' ( )'                                  ,& ! 201
 'ITAB_M_IW'     ,'M2' ,'ITAB_M_IW',' ( )'                                  ,& ! 202
 'ITAB_M_FMW'    ,'M2' ,'ITAB_M_FMW',' ( )'                                 ,& ! 203

! ITAB_U MEMBERS - 2D

 'ITAB_U_IUP'    ,'V2' ,'ITAB_U_IUP',' ( )'                                 ,& ! 204
 'ITAB_U_IRANK'  ,'V2' ,'ITAB_U_IRANK',' ( )'                               ,& ! 205
 'ITAB_U_IUGLOBE','V2' ,'ITAB_U_IUGLOBE',' ( )'                             ,& ! 206
 'ITAB_U_MRLU'   ,'V2' ,'ITAB_U_MRLU',' ( )'                                ,& ! 207
 'ITAB_U_IM'     ,'V2' ,'ITAB_U_IM',' ( )'                                  ,& ! 208
 'ITAB_U_IU'     ,'V2' ,'ITAB_U_IU',' ( )'                                  ,& ! 209
 'ITAB_U_IW'     ,'V2' ,'ITAB_U_IW',' ( )'                                  ,& ! 210
 'ITAB_U_DIRU'   ,'V2' ,'ITAB_U_DIRU',' ( )'                                ,& ! 211
 'ITAB_U_FUU'    ,'V2' ,'ITAB_U_FUU',' ( )'                                 ,& ! 212
 'ITAB_U_FUW'    ,'V2' ,'ITAB_U_FUW',' ( )'                                 ,& ! 213
 'ITAB_U_TUU'    ,'V2' ,'ITAB_U_TUU',' ( )'                                 ,& ! 214
 'ITAB_U_GUW'    ,'V2' ,'ITAB_U_GUW',' ( )'                                 ,& ! 215
 'ITAB_U_PGC12'  ,'V2' ,'ITAB_U_PGC12',' (m:S2:-1  )'                       ,& ! 216
 'ITAB_U_PGC12b' ,'V2' ,'ITAB_U_PGC12b',' (m:S2:-1  )'                      ,& ! 217
 'ITAB_U_PGC12c' ,'V2' ,'ITAB_U_PGC12c',' (m:S2:-1  )'                      ,& ! 218
 'ITAB_U_PGC12d' ,'V2' ,'ITAB_U_PGC12d',' (m:S2:-1  )'                      ,& ! 219
 'ITAB_U_PGC45'  ,'V2' ,'ITAB_U_PGC45',' (m:S2:-1  )'                       ,& ! 220
 'ITAB_U_PGC45b' ,'V2' ,'ITAB_U_PGC45b',' (m:S2:-1  )'                      ,& ! 221
 'ITAB_U_PGC63'  ,'V2' ,'ITAB_U_PGC63',' (m:S2:-1  )'                       ,& ! 222
 'ITAB_U_PGC63c' ,'V2' ,'ITAB_U_PGC63c',' (m:S2:-1  )'                      ,& ! 223
 'ITAB_U_VXU_U'  ,'V2' ,'ITAB_U_VXU_U',' ( )'                               ,& ! 224
 'ITAB_U_VYU_U'  ,'V2' ,'ITAB_U_VYU_U',' ( )'                               ,& ! 225
 'ITAB_U_VXW_U'  ,'V2' ,'ITAB_U_VXW_U',' ( )'                               ,& ! 226
 'ITAB_U_VYW_U'  ,'V2' ,'ITAB_U_VYW_U',' ( )'                               ,& ! 227
 'ITAB_U_CROSSMM','V2' ,'ITAB_U_CROSSMM',' ( )'                             ,& ! 228
 'ITAB_U_CROSSWW','V2' ,'ITAB_U_CROSSWW',' ( )'                              / ! 229

! ITAB_V MEMBERS - 2D

data fldlib(1:4,230:240)/  &
 'ITAB_V_IVP'    ,'V2' ,'ITAB_V_IVP',' ( )'                                 ,& ! 230
 'ITAB_V_IRANK'  ,'V2' ,'ITAB_V_IRANK',' ( )'                               ,& ! 231
 'ITAB_V_IVGLOBE','V2' ,'ITAB_V_IVGLOBE',' ( )'                             ,& ! 232
 'ITAB_V_MRLV'   ,'V2' ,'ITAB_V_MRLV',' ( )'                                ,& ! 233
 'ITAB_V_IM'     ,'V2' ,'ITAB_V_IM',' ( )'                                  ,& ! 234
 'ITAB_V_IV'     ,'V2' ,'ITAB_V_IV',' ( )'                                  ,& ! 235
 'ITAB_V_IW'     ,'V2' ,'ITAB_V_IW',' ( )'                                  ,& ! 236
 'ITAB_V_FVV'    ,'V2' ,'ITAB_V_FVV',' ( )'                                 ,& ! 237
 'ITAB_V_FVW'    ,'V2' ,'ITAB_V_FVW',' ( )'                                 ,& ! 238
 'ITAB_V_FUV'    ,'V2' ,'ITAB_V_FUV',' ( )'                                 ,& ! 239
 'ITAB_V_FARW'   ,'V2' ,'ITAB_V_FARW',' ( )'                                 / ! 240

! ITAB_W MEMBERS - 2D

data fldlib(1:4,241:268)/ &

 'ITAB_W_NPOLY'  ,'W2' ,'ITAB_W_NPOLY',' ( )'                               ,& ! 241
 'ITAB_W_IWP'    ,'W2' ,'ITAB_W_IWP',' ( )'                                 ,& ! 242
 'ITAB_W_IRANK'  ,'W2' ,'ITAB_W_IRANK',' ( )'                               ,& ! 243
 'ITAB_W_IWGLOBE','W2' ,'ITAB_W_IWGLOBE',' ( )'                             ,& ! 244
 'ITAB_W_MRLW'   ,'W2' ,'ITAB_W_MRLW',' ( )'                                ,& ! 245
 'ITAB_W_MROW'   ,'W2' ,'ITAB_W_MROW',' ( )'                                ,& ! 246
 'ITAB_W_MROWH'  ,'W2' ,'ITAB_W_MROWH',' ( )'                               ,& ! 247
 'ITAB_W_IM'     ,'W2' ,'ITAB_W_IM',' ( )'                                  ,& ! 248
 'ITAB_W_IU'     ,'W2' ,'ITAB_W_IU',' ( )'                                  ,& ! 249
 'ITAB_W_IV'     ,'W2' ,'ITAB_W_IV',' ( )'                                  ,& ! 250
 'ITAB_W_IW'     ,'W2' ,'ITAB_W_IW',' ( )'                                  ,& ! 251
 'ITAB_W_DIRU'   ,'W2' ,'ITAB_W_DIRU',' ( )'                                ,& ! 252
 'ITAB_W_DIRV'   ,'W2' ,'ITAB_W_DIRV',' ( )'                                ,& ! 253
 'ITAB_W_FWV'    ,'W2' ,'ITAB_W_FWV',' ( )'                                 ,& ! 254
 'ITAB_W_FWW'    ,'W2' ,'ITAB_W_FWW',' ( )'                                 ,& ! 255
 'ITAB_W_FWU'    ,'W2' ,'ITAB_W_FWU',' ( )'                                 ,& ! 256
 'ITAB_W_VXU'    ,'W2' ,'ITAB_W_VXU',' ( )'                                 ,& ! 257
 'ITAB_W_VYU'    ,'W2' ,'ITAB_W_VYU',' ( )'                                 ,& ! 258
 'ITAB_W_VZU'    ,'W2' ,'ITAB_W_VZU',' ( )'                                 ,& ! 259
 'ITAB_W_VXW'    ,'W2' ,'ITAB_W_VXW',' ( )'                                 ,& ! 260
 'ITAB_W_VYW'    ,'W2' ,'ITAB_W_VYW',' ( )'                                 ,& ! 261
 'ITAB_W_VZW'    ,'W2' ,'ITAB_W_VZW',' ( )'                                 ,& ! 262
 'ITAB_W_VXU_W'  ,'W2' ,'ITAB_W_VXU_W',' ( )'                               ,& ! 263
 'ITAB_W_VYU_W'  ,'W2' ,'ITAB_W_VYU_W',' ( )'                               ,& ! 264
 'ITAB_W_FARM'   ,'W2' ,'ITAB_W_FARM',' ( )'                                ,& ! 265
 'ITAB_W_FARV'   ,'W2' ,'ITAB_W_FARV',' ( )'                                ,& ! 266
 'ITAB_W_IWNUD'  ,'W2' ,'ITAB_W_IWNUD',' ( )'                               ,& ! 267
 'ITAB_W_FNUD'   ,'W2' ,'ITAB_W_FNUD',' ( )'                                 / ! 268

! Time-averaged fields

data fldlib(1:4,269:278)/ &
 'RSHORT_AVG'      ,'T2' ,'AVG SFC DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'   ,& ! 269
 'RSHORTUP_AVG'    ,'T2' ,'AVG SFC UPWARD SHORTWV FLX',' (W m:S2:-2  )'     ,& ! 270
 'RLONG_AVG'       ,'T2' ,'AVG SFC DOWNWARD LONGWV FLX',' (W m:S2:-2  )'    ,& ! 271
 'RLONGUP_AVG'     ,'T2' ,'AVG SFC UPWARD LONGWV FLX',' (W m:S2:-2  )'      ,& ! 272
 'RSHORT_TOP_AVG'  ,'T2' ,'AVG TOP DOWNWARD SHORTWV FLX',' (W m:S2:-2  )'   ,& ! 273
 'RSHORTUP_TOP_AVG','T2' ,'AVG TOP UPWARD SHORTWV FLX',' (W m:S2:-2  )'     ,& ! 274
 'RLONGUP_TOP_AVG' ,'T2' ,'AVG TOP UPWARD LONGWV FLX',' (W m:S2:-2  )'      ,& ! 275
 'SENSFLUX_AVG'    ,'T2' ,'AVG ATM SFC SENS HEAT FLUX',' (N m:S2:-2  )'     ,& ! 276
 'LATFLUX_AVG'     ,'T2' ,'AVG ATM SFC LAT HEAT FLUX',' (W m:S2:-2  )'      ,& ! 277
 'VAPFLUX_AVG'     ,'T2' ,'AVG ATM SFC VAP FLUX',' (kg m:S2:-2   s:S2:-1  )' / ! 278

! Miscellaneous and new additions

data fldlib(1:4,279:293)/ &
 'RHO_OBS'       ,'T3' ,'NUDGING OBS AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 279
 'THETA_OBS'     ,'T3' ,'NUDGING OBS THETA',' (K)'                          ,& ! 280
 'SHW_OBS'       ,'T3' ,'NUDGING OBS VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 281
 'UZONAL_OBS'    ,'T3' ,'NUDGING OBS ZONAL WIND',' (m s:S2:-1  )'           ,& ! 282
 'UMERID_OBS'    ,'T3' ,'NUDGING OBS MERID WIND',' (m s:S2:-1  )'           ,& ! 283
 'RHO_SIM'       ,'T3' ,'NUDGING SIM AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 284
 'THETA_SIM'     ,'T3' ,'NUDGING SIM THETA',' (K)'                          ,& ! 285
 'SHW_SIM'       ,'T3' ,'NUDGING SIM VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 286
 'UZONAL_SIM'    ,'T3' ,'NUDGING SIM ZONAL WIND',' (m s:S2:-1  )'           ,& ! 287
 'UMERID_SIM'    ,'T3' ,'NUDGING SIM MERID WIND',' (m s:S2:-1  )'           ,& ! 288
 'RHO_OBS_SIM'   ,'T3' ,'NUDGING DIF AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 289
 'THETA_OBS_SIM' ,'T3' ,'NUDGING DIF THETA',' (K)'                          ,& ! 290
 'SHW_OBS_SIM'   ,'T3' ,'NUDGING DIF VAPOR SPEC DENSITY',' (g kg:S2:-1  )'  ,& ! 291
 'UZONAL_OBS_SIM','T3' ,'NUDGING DIF ZONAL WIND',' (m s:S2:-1  )'           ,& ! 292
 'UMERID_OBS_SIM','T3' ,'NUDGING DIF MERID WIND',' (m s:S2:-1  )'            / ! 293

! External fields

data fldlib(1:4,294:296)/ &
 'VORTP'         ,'P3' ,'VORTP',' (s:S2:-1  )'                              ,& ! 294
 'VORTN'         ,'N3' ,'VORTN',' (s:S2:-1  )'                              ,& ! 295
 'RKE'           ,'T3' ,'RKE',' (s:S2:-1  )'                                 / ! 296

if (ncall /= 10) then
   ncall = 10
   allocate (press_init(mza,mwa))
   allocate (rho_init  (mza,mwa))
   allocate (theta_init(mza,mwa))
   allocate (uc_init   (mza,mua))
   allocate (vc_init   (mza,mua))
!   allocate (addsc1_init(mza,mwa))

   press_init(:,:) = press(:,:)
   rho_init  (:,:) = rho  (:,:)
   theta_init(:,:) = theta(:,:)
   uc_init   (:,:) = uc   (:,:)
   vc_init   (:,:) = vc   (:,:)
!   addsc1_init(:,:) = addsc(1)%sclp(:,:)

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

   fldval = wtbot * umc(k  ,i) &
          + wttop * umc(k+1,i)

case(2) ! 'VMC'

   fldval = wtbot * vmc(k  ,i) &
          + wttop * vmc(k+1,i)

case(3) ! 'WMC'

! Need to re-examine use of k-1 when k = lpw

   fldval = wtbot * wmc(k  ,i) &
          + wttop * wmc(k+1,i)

case(4) ! 'UMP'

   fldval = wtbot * ump(k  ,i) &
          + wttop * ump(k+1,i)

case(5) ! 'VMP'

   fldval = wtbot * vmp(k  ,i) &
          + wttop * vmp(k+1,i)

case(6) ! 'UC'

   fldval = wtbot * uc(k  ,i) &
          + wttop * uc(k+1,i)

case(7) ! 'VC'

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

case(11) ! 'THIL'

   fldval = wtbot * thil(k  ,i) &
          + wttop * thil(k+1,i)

case(12) ! 'THETA'

   fldval = wtbot * theta(k  ,i) &
          + wttop * theta(k+1,i)

case(13) ! 'AIRTEMPK'

   fldval = wtbot * theta(k  ,i) * (press(k  ,i) / p00) ** rocp &
          + wttop * theta(k+1,i) * (press(k+1,i) / p00) ** rocp

case(14) ! 'AIRTEMPC'

   fldval = wtbot * theta(k  ,i) * (press(k  ,i) / p00) ** rocp &
          + wttop * theta(k+1,i) * (press(k+1,i) / p00) ** rocp - 273.15

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

   if (.not. allocated(fthrd)) go to 1000

   fldval = wtbot * fthrd(k  ,i) &
          + wttop * fthrd(k+1,i)

case(40:45) ! 'SPEEDV','AZIMV','ZONAL_WINDV','MERID_WINDV','ZONAL_WINDV_P','MERID_WINDV_P'

   ucc = wtbot * uc(k  ,i) &
       + wttop * uc(k+1,i)

   vcc = wtbot * vc(k  ,i) &
       + wttop * vc(k+1,i)

   vx = unx(i) * ucc + vnx(i) * vcc
   vy = uny(i) * ucc + vny(i) * vcc
   vz = unz(i) * ucc + vnz(i) * vcc

   if (trim(fldname) == 'SPEEDV') then
      fldval = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
   else

      u = 0.
      v = 0.

      if (meshtype == 1) then

         raxis = sqrt(xeu(i) ** 2 + yeu(i) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            u = (vy * xeu(i) - vx * yeu(i)) / raxis
            v = vz * raxis / erad &
              - (vx * xeu(i) + vy * yeu(i)) * zeu(i) / (raxis * erad) 
         endif

      else

         raxis = sqrt(xev(i) ** 2 + yev(i) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            u = (vy * xev(i) - vx * yev(i)) / raxis
            v = vz * raxis / erad &
              - (vx * xev(i) + vy * yev(i)) * zev(i) / (raxis * erad)
         endif

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
            endif

         else

            raxis = sqrt(xev(i) ** 2 + yev(i) ** 2)  ! dist from earth axis

            if (raxis > 1.e3) then
               u_init = (vy_init * xev(i) - vx_init * yev(i)) / raxis
               v_init = vz_init * raxis / erad &
                  - (vx_init * xev(i) + vy_init * yev(i)) * zev(i) / (raxis * erad) 
            endif

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

   vxe = 0.
   vye = 0.
   vze = 0.

   do j = 1,npoly

      if (meshtype == 1) then

         iv = itab_w(i)%iu(j)
         
         ucint = wtbot * uc(k,iv) + wttop * uc(k+1,iv) 
         vcint = wtbot * vc(k,iv) + wttop * vc(k+1,iv) 

         vxe = vxe + rpolyi * (unx(iv) * ucint + vnx(iv) * vcint)
         vye = vye + rpolyi * (uny(iv) * ucint + vny(iv) * vcint)
         vze = vze + rpolyi * (unz(iv) * ucint + vnz(iv) * vcint)

      else

         iv = itab_w(i)%iv(j)
         farv2 = 2. * itab_w(i)%farv(j)

         vcint = wtbot * vc(k,iv) + wttop * vc(k+1,iv) 
         
         vxe = vxe + farv2 * vnx(iv) * vcint
         vye = vye + farv2 * vny(iv) * vcint
         vze = vze + farv2 * vnz(iv) * vcint

      endif

   enddo

   if (trim(fldname) == 'SPEEDW') then
      fldval = sqrt(vxe ** 2 + vye ** 2 + vze ** 2)
   else
      raxis = sqrt(xew(i) ** 2 + yew(i) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
         u = (vye * xew(i) - vxe * yew(i)) / raxis
         v = vze * raxis / erad &
           - (vxe * xew(i) + vye * yew(i)) * zew(i) / (raxis * erad) 
      else
         u = 0.
         v = 0.
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
         iv = itab_m(i)%iu(j)

         iw1 = itab_u(iv)%iw(1)
         iw2 = itab_u(iv)%iw(2)
      else
         iv = itab_m(i)%iv(j)
      endif

      ucc = wtbot * uc(k  ,iv) &
          + wttop * uc(k+1,iv)

      vcc = wtbot * vc(k  ,iv) &
          + wttop * vc(k+1,iv)

      if (trim(fldname) == 'RVORTZM_P') then

         ucc_init = wtbot * uc_init(k  ,iv) &
                  + wttop * uc_init(k+1,iv)

         vcc_init = wtbot * vc_init(k  ,iv) &
                  + wttop * vc_init(k+1,iv)

         ucc = ucc - ucc_init
         vcc = vcc - vcc_init

      endif

! Now reconstruct total wind vector projected onto vector from IW1 to IW2

      if (meshtype == 1) then

         contrib = ucc * (unx(iv) * (xew(iw2) - xew(iw1))  &
                       +  uny(iv) * (yew(iw2) - yew(iw1))  &
                       +  unz(iv) * (zew(iw2) - zew(iw1))) &
                 + vcc * (vnx(iv) * (xew(iw2) - xew(iw1))  &
                       +  vny(iv) * (yew(iw2) - yew(iw1))  &
                       +  vnz(iv) * (zew(iw2) - zew(iw1)))

         if (i == itab_u(iv)%im(1)) then
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

   if (trim(fldname) == 'TVORTZM') then

      fldval = fldval + omega2 * zem(i) / erad  ! add earth vorticity at M point

   endif

case(53) ! 'DIVERG'

   fldval = 0.

   npoly = itab_w(i)%npoly

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

case(54) ! 'UMASSFLUX'

   fldval = umc(k,i) * aru(k,i)

case(55) ! 'VMASSFLUX'

   fldval = vmc(k,i) * arv(k,i)

case(56) ! 'UC_P'

   fldval = uc(k,i) - uc_init(k,i)

case(57) ! 'VC_P'

   fldval = vc(k,i) - vc_init(k,i)

case(58) ! 'PRESS_P'

   fldval = press(k,i) - press_init(k,i)

case(59) ! 'RHO_P'

   fldval = wtbot * (rho(k  ,i) - rho_init(k  ,i)) &
          + wttop * (rho(k+1,i) - rho_init(k+1,i))

case(60) ! 'THETA_P'

   fldval = wtbot * (theta(k  ,i) - theta_init(k  ,i)) &
          + wttop * (theta(k+1,i) - theta_init(k+1,i))

case(61) ! 'AIRTEMPK_P'

   fldval = wtbot * theta(k  ,i) * (press(k  ,i) / p00) ** rocp &
          + wttop * theta(k+1,i) * (press(k+1,i) / p00) ** rocp &
          - wtbot * theta_init(k  ,i) * (press_init(k  ,i) / p00) ** rocp &
          - wttop * theta_init(k+1,i) * (press_init(k+1,i) / p00) ** rocp

case(62) ! 'UMT'

   fldval = umt(k,i) * volui(k,i)

case(63) ! 'VMT'

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

   fldval = vkm_sfc(i)

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

case(88) ! 'PCPRTOT'

   fldval = 0.

   if (allocated(pcprd)) fldval = fldval + pcprd(i)
   if (allocated(pcprr)) fldval = fldval + pcprr(i)
   if (allocated(pcprp)) fldval = fldval + pcprp(i)
   if (allocated(pcprs)) fldval = fldval + pcprs(i)
   if (allocated(pcpra)) fldval = fldval + pcpra(i)
   if (allocated(pcprg)) fldval = fldval + pcprg(i)
   if (allocated(pcprh)) fldval = fldval + pcprh(i)

   fldval = fldval * 3600.

case(89) ! 'CONPRR'

   if (.not. allocated(conprr)) go to 1000

   fldval = conprr(i) * 3600.

case(90) ! 'TOTPRR'

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

case(98) ! 'ACCPTOT'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)

case(99) ! 'ACONPR'

   if (.not. allocated(aconpr)) go to 1000

   fldval = aconpr(i)

case(100) ! 'ATOTPR'

   fldval = 0.

   if (allocated(accpd)) fldval = fldval + accpd(i)
   if (allocated(accpr)) fldval = fldval + accpr(i)
   if (allocated(accpp)) fldval = fldval + accpp(i)
   if (allocated(accps)) fldval = fldval + accps(i)
   if (allocated(accpa)) fldval = fldval + accpa(i)
   if (allocated(accpg)) fldval = fldval + accpg(i)
   if (allocated(accph)) fldval = fldval + accph(i)
   if (allocated(aconpr)) fldval = fldval + aconpr(i)

! LAND_CELLS - 3D

case(101) ! 'SOIL_TEXT'

   fldval = real(land%ntext_soil(k,i))

case(102) ! 'SOIL_ENERGY'

   fldval = land%soil_energy(k,i) * 1.e-6

case(103) ! 'SOIL_TEMPC'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = tempk - 273.15

case(104) ! 'SOIL_FRACLIQ'

   call qwtk(land%soil_energy(k,i)       &
            ,land%soil_water(k,i)*1.e3   &
            ,slcpd(land%ntext_soil(k,i)) &
            ,tempk, fracliq)
   fldval = fracliq

case(105) ! 'SOIL_WATER'

   fldval = land%soil_water(k,i)

case(106) ! 'SFWAT_MASS'

   fldval = land%sfcwater_mass(k,i)

case(107) ! 'SFWAT_ENERGY'

   fldval = land%sfcwater_energy(k,i) * 1.e-3

case(108) ! 'SFWAT_TEMPC'

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = tempk - 273.15

case(109) ! 'SFWAT_FRACLIQ'

   call qtk(land%sfcwater_energy(k,i),tempk,fracliq)
   fldval = fracliq

case(110) ! 'SFWAT_DEPTH'

   fldval = 0.
   do klev = 1,land%nlev_sfcwater(i)
      fldval = fldval + land%sfcwater_depth(klev,i)
   enddo

! LAND_CELLS - 2D

case(111) ! 'NLEV_SFWAT'

   fldval = real(land%nlev_sfcwater(i))

case(112) ! 'VEG_NDVIC'

   fldval = land%veg_ndvic(i)

case(113) ! 'VEG_TEMPC'

   fldval = land%veg_temp(i) - 273.15

case(114) ! 'VEG_WATER'

   fldval = land%veg_water(i)

case(115) ! 'STOM_RESIST'

   fldval = land%stom_resist(i)

case(116) ! 'GROUND_SHV'

   fldval = land%ground_shv(i) * 1.e3

case(117) ! 'SOIL_DEPTH'

   fldval = -slz(land%lsl(i))

! SEA_CELLS - 2D

case(118) ! 'SEATP'

   fldval = sea%seatp(i)   

case(119) ! 'SEATF'

   fldval = sea%seatf(i)

case(120) ! 'SEATC'

   fldval = sea%seatc(i)

case(121) ! 'SEAICEP'

   fldval = sea%seaicep(i)   

case(122) ! 'SEAICEF'

   fldval = sea%seaicef(i)

case(123) ! 'SEAICEC'

   fldval = sea%seaicec(i)

! LAND AND SEA CELLS - 2D

case(124) ! 'LEAF_CLASS'

   if (op%stagpt == 'S') then
      fldval = real(sea%leaf_class(i))
   elseif (op%stagpt == 'L') then
      fldval = real(land%leaf_class(i))
   endif

case(125) ! 'AREA'

   if (op%stagpt == 'S') then
      fldval = sea%area(i)
   elseif (op%stagpt == 'L') then
      fldval = land%area(i)
   endif

case(126) ! 'ROUGH'

   if (op%stagpt == 'S') then
      fldval = sea%rough(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rough(i)
   endif

case(127) ! 'CAN_TEMPC'

   if (op%stagpt == 'S') then
      fldval = sea%can_temp(i) - 273.15
   elseif (op%stagpt == 'L') then
      fldval = land%can_temp(i) - 273.15
   endif

case(128) ! 'CAN_SHV'

   if (op%stagpt == 'S') then
      fldval = sea%can_shv(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%can_shv(i) * 1.e3
   endif

case(129) ! 'SFC_TEMPC'

   if (op%stagpt == 'S') then
      fldval = sea%seatc(i) - 273.15
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

case(130) ! 'SFC_SSH'

   if (op%stagpt == 'S') then
      fldval = sea%surface_ssh(i) * 1.e3
   elseif (op%stagpt == 'L') then
      fldval = land%surface_ssh(i) * 1.e3
   endif

case(131) ! 'SENSFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_tsav(i) * cp / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_tsav(i) * cp / dt_leaf
   endif

case(132) ! 'VAPFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_rsav(i) / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_rsav(i) / dt_leaf
   endif

case(133) ! 'LATFLUX_LS'

   if (op%stagpt == 'S') then
      fldval = sea%sxfer_rsav(i) * alvl / dt_leaf
   elseif (op%stagpt == 'L') then
      fldval = land%sxfer_rsav(i) * alvl / dt_leaf
   endif

case(134) ! 'RSHORT_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rshort(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort(i)
   endif

case(135) ! 'RSHORT_DIFFUSE_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rshort_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rshort_diffuse(i)
   endif

case(136) ! 'RLONG_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlong(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong(i)
   endif

case(137) ! 'RLONGUP_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlongup(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlongup(i)
   endif

case(138) ! 'RLONG_ALBEDO_LS'

   if (op%stagpt == 'S') then
      fldval = sea%rlong_albedo(i)
   elseif (op%stagpt == 'L') then
      fldval = land%rlong_albedo(i)
   endif

case(139) ! 'ALBEDO_BEAM_LS'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_beam(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_beam(i)
   endif

case(140) ! 'ALBEDO_DIFFUSE_LS'

   if (op%stagpt == 'S') then
      fldval = sea%albedo_diffuse(i)
   elseif (op%stagpt == 'L') then
      fldval = land%albedo_diffuse(i)
   endif

! FLUX CELLS

case(141) ! 'FCELL_ILSF'

   fldval = real(i)

case(142) ! 'FCELL_IWLS'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iwls)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iwls)
   endif

case(143) ! 'FCELL_IW'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%iw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%iw)
   endif

case(144) ! 'FCELL_KW'

   if (op%stagpt == 'S') then
      fldval = real(seaflux(i)%kw)
   elseif (op%stagpt == 'L') then
      fldval = real(landflux(i)%kw)
   endif

case(145) ! FCELL_AREA'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%area
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%area
   endif

case(146) ! 'FCELL_ARFATM'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_atm
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_atm
   endif

case(147) ! 'FCELL_ARFLS'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%arf_sfc
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%arf_sfc
   endif

case(148) ! 'FCELL_SENS'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxt * cp
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxt * cp
   endif

case(149) ! 'FCELL_VAP'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr
   endif

case(150) ! 'FCELL_LAT'

   if (op%stagpt == 'S') then
      fldval = seaflux(i)%sfluxr * alvl
   elseif (op%stagpt == 'L') then
      fldval = landflux(i)%sfluxr * alvl
   endif

case(151) ! 'FCELL_AIRTEMPC'

   if (op%stagpt == 'S') then
      kw = seaflux(i)%kw
      iw = seaflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp - 273.15
   elseif (op%stagpt == 'L') then
      kw = landflux(i)%kw
      iw = landflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp - 273.15
   endif

case(152) ! 'FCELL_AIRTEMPK'

   if (op%stagpt == 'S') then
      kw = seaflux(i)%kw
      iw = seaflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp
   elseif (op%stagpt == 'L') then
      kw = landflux(i)%kw
      iw = landflux(i)%iw
      if (iparallel == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif
      fldval = theta(kw,iw) * (press(kw,iw) / p00) ** rocp
   endif

case(153) ! 'FCELL_CANTEMPC'

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

case(154) ! 'ARU'
   fldval = aru(k,i)
case(155) ! 'ARV'
   if (.not. allocated(arv)) go to 1000
   fldval = arv(k,i)
case(156) ! 'ARW'
   fldval = arw(k,i)
case(157) ! 'VOLT'
   fldval = volt(k,i)
case(158) ! 'VOLU'
   if (.not. allocated(volui)) go to 1000
   fldval = 1. / volui(k,i)
case(159) ! 'VOLV'
   if (.not. allocated(volvi)) go to 1000
   fldval = 1. / volvi(k,i)
case(160) ! 'VOLW'
   fldval = 1. / volwi(k,i)

! GRID GEOMETRY - 2D

case(161) ! 'TOPM'
   fldval = topm(i)
case(162) ! 'TOPW'
   fldval = topw(i)
case(163) ! 'GLATM'
   fldval = glatm(i)
case(164) ! 'GLONM'
   fldval = glonm(i)
case(165) ! 'GLATU'
   if (.not. allocated(glatu)) go to 1000
   fldval = glatu(i)
case(166) ! 'GLONU'
   if (.not. allocated(glonu)) go to 1000
   fldval = glonu(i)
case(167) ! 'GLATV'
   if (.not. allocated(glatv)) go to 1000
   fldval = glatv(i)
case(168) ! 'GLONV'
   if (.not. allocated(glonv)) go to 1000
   fldval = glonv(i)
case(169) ! 'GLATW'
   fldval = glatw(i)
case(170) ! 'GLONW'
   fldval = glonw(i)
case(171) ! 'LPM'
   fldval = real(lpm(i))
case(172) ! 'LPU'
   if (.not. allocated(lpu)) go to 1000
   fldval = real(lpu(i))
case(173) ! 'LCU'
   if (.not. allocated(lcu)) go to 1000
   fldval = real(lcu(i))
case(174) ! 'LPV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(175) ! 'LCV'
   if (.not. allocated(lcv)) go to 1000
   fldval = real(lcv(i))
case(176) ! 'LPW'
   fldval = real(lpw(i))
case(177) ! 'LSW'
   fldval = real(lsw(i))
case(178) ! 'XEM'
   fldval = xem(i)
case(179) ! 'YEM'
   fldval = yem(i)
case(180) ! 'ZEM'
   fldval = zem(i)
case(181) ! 'XEU'
   if (.not. allocated(xeu)) go to 1000
   fldval = xeu(i)
case(182) ! 'YEU'
   if (.not. allocated(yeu)) go to 1000
   fldval = yeu(i)
case(183) ! 'ZEU'
   if (.not. allocated(zeu)) go to 1000
   fldval = zeu(i)
case(184) ! 'XEV'
   if (.not. allocated(xev)) go to 1000
   fldval = xev(i)
case(185) ! 'YEV'
   if (.not. allocated(yev)) go to 1000
   fldval = yev(i)
case(186) ! 'ZEV'
   if (.not. allocated(zev)) go to 1000
   fldval = zev(i)
case(187) ! 'XEW'
   fldval = xew(i)
case(188) ! 'YEW'
   fldval = yew(i)
case(189) ! 'ZEW'
   fldval = zew(i)
case(190) ! 'DNU'
   fldval = dnu(i)
case(191) ! 'DNV'
   fldval = dnv(i)
case(192) ! 'ARM0'
   fldval = arm0(i)
case(193) ! 'ARW0'
   fldval = arw0(i)

! ITAB_M MEMBERS

case(194) ! 'ITAB_M_NPOLY'
   fldval = itab_m(i)%npoly
case(195) ! 'ITAB_M_IMGLOBE'
   fldval = itab_m(i)%imglobe
case(196) ! 'ITAB_M_MRLMR'
   fldval = itab_m(i)%mrlm
case(197) ! 'ITAB_M_MRLM_OR'
   fldval = itab_m(i)%mrlm_orig
case(198) ! 'ITAB_M_MROW'
   fldval = itab_m(i)%mrow
case(199) ! 'ITAB_M_MROWH'
   fldval = itab_m(i)%mrowh
case(200) ! 'ITAB_M_IU'
   fldval = itab_m(i)%iu(indp)
case(201) ! 'ITAB_M_IV'
   fldval = itab_m(i)%iv(indp)
case(202) ! 'ITAB_M_IW'
   fldval = itab_m(i)%iw(indp)
case(203) ! 'ITAB_M_FMW'
   fldval = itab_m(i)%fmw(indp)

! ITAB_U MEMBERS

case(204) ! 'ITAB_U_IUP'
   fldval = itab_u(i)%iup
case(205) ! 'ITAB_U_IRANK'
   fldval = itab_u(i)%irank
!   fldval = itabg_u(i)%irank
case(206) ! 'ITAB_U_IUGLOBE'
   fldval = itab_u(i)%iuglobe
case(207) ! 'ITAB_U_MRLU'
   fldval = itab_u(i)%mrlu
case(208) ! 'ITAB_U_IM'
   fldval = itab_u(i)%im(indp)
case(209) ! 'ITAB_U_IU'
   fldval = itab_u(i)%iu(indp)
case(210) ! 'ITAB_U_IW'
   fldval = itab_u(i)%iw(indp)
case(211) ! 'ITAB_U_DIRU'
   fldval = itab_u(i)%diru(indp)
case(212) ! 'ITAB_U_FUU'
   fldval = itab_u(i)%fuu(indp)
case(213) ! 'ITAB_U_FUW'
   fldval = itab_u(i)%fuw(indp)
case(214) ! 'ITAB_U_TUU'
   fldval = itab_u(i)%tuu(indp)
case(215) ! 'ITAB_U_GUW'
   fldval = itab_u(i)%guw(indp)
case(216) ! 'ITAB_U_PGC12'
   fldval = itab_u(i)%pgc12
case(217) ! 'ITAB_U_PGC12b'
   fldval = itab_u(i)%pgc12b
case(218) ! 'ITAB_U_PGC12c'
   fldval = itab_u(i)%pgc12c
case(219) ! 'ITAB_U_PGC12d'
   fldval = itab_u(i)%pgc12d
case(220) ! 'ITAB_U_PGC45'
   fldval = itab_u(i)%pgc45
case(221) ! 'ITAB_U_PGC45b'
   fldval = itab_u(i)%pgc45b
case(222) ! 'ITAB_U_PGC63'
   fldval = itab_u(i)%pgc63
case(223) ! 'ITAB_U_PGC63c'
   fldval = itab_u(i)%pgc63c
case(224) ! 'ITAB_U_VXU_U'
   fldval = itab_u(i)%vxu_u(indp)
case(225) ! 'ITAB_U_VYU_U'
   fldval = itab_u(i)%vyu_u(indp)
case(226) ! 'ITAB_U_VXW_U''
   fldval = itab_u(i)%vxw_u(indp)
case(227) ! 'ITAB_U_VYW_U'
   fldval = itab_u(i)%vyw_u(indp)
case(228) ! 'ITAB_U_CROSSMM'
   fldval = itab_u(i)%crossmm
case(229) ! 'ITAB_U_CROSSWW'
   fldval = itab_u(i)%crossww

! ITAB_V MEMBERS

case(230) ! 'ITAB_V_IVP'
   fldval = itab_v(i)%ivp
case(231) ! 'ITAB_V_IRANK'
   fldval = itab_v(i)%irank
!   fldval = itabg_v(i)%irank
case(232) ! 'ITAB_V_IVGLOBE'
   fldval = itab_v(i)%ivglobe
case(233) ! 'ITAB_V_MRLV'
   fldval = itab_v(i)%mrlv
case(234) ! 'ITAB_V_IM'
   fldval = itab_v(i)%im(indp)
case(235) ! 'ITAB_V_IV'
   fldval = itab_v(i)%iv(indp)
case(236) ! 'ITAB_V_IW'
   fldval = itab_v(i)%iw(indp)
case(237) ! 'ITAB_V_FVV'
   fldval = itab_v(i)%fvv(indp)
case(238) ! 'ITAB_V_FVW'
   fldval = itab_v(i)%fvw(indp)
case(239) ! 'ITAB_V_FUV'
   fldval = itab_v(i)%fuv(indp)
case(240) ! 'ITAB_V_FARW'
   fldval = itab_v(i)%farw(indp)

! ITAB_W MEMBERS

case(241) ! 'ITAB_W_NPOLY'
   fldval = itab_w(i)%npoly
case(242) ! 'ITAB_W_IWP'
   fldval = itab_w(i)%iwp
case(243) ! 'ITAB_W_IRANK'
   fldval = itab_w(i)%irank
!   fldval = itabg_w(i)%irank
case(244) ! 'ITAB_W_IWGLOBE'
   fldval = itab_w(i)%iwglobe
case(245) ! 'ITAB_W_MRLW'
   fldval = itab_w(i)%mrlw
case(246) ! 'ITAB_W_MROW'
   fldval = itab_w(i)%mrow
case(247) ! 'ITAB_W_MROWH'
   fldval = itab_w(i)%mrowh
case(248) ! 'ITAB_W_IM'
   fldval = itab_w(i)%im(indp)
case(249) ! 'ITAB_W_IU'
   fldval = itab_w(i)%iu(indp)
case(250) ! 'ITAB_W_IV'
   fldval = itab_w(i)%iv(indp)
case(251) ! 'ITAB_W_IW'
   fldval = itab_w(i)%iw(indp)
case(252) ! 'ITAB_W_DIRU'
   fldval = itab_w(i)%diru(indp)
case(253) ! 'ITAB_W_DIRV'
   fldval = itab_w(i)%dirv(indp)
case(254) ! 'ITAB_W_FWV'
   fldval = itab_w(i)%fwv(indp)
case(255) ! 'ITAB_W_FWW'
   fldval = itab_w(i)%fww(indp)
case(256) ! 'ITAB_W_FWU'
   fldval = itab_w(i)%fwu(indp)
case(257) ! 'ITAB_W_VXU'
   fldval = itab_w(i)%vxu(indp)
case(258) ! 'ITAB_W_VYU'
   fldval = itab_w(i)%vyu(indp)
case(259) ! 'ITAB_W_VZU'
   fldval = itab_w(i)%vzu(indp)
case(260) ! 'ITAB_W_VXW'
   fldval = itab_w(i)%vxw
case(261) ! 'ITAB_W_VYW'
   fldval = itab_w(i)%vyw
case(262) ! 'ITAB_W_VZW'
   fldval = itab_w(i)%vzw
case(263) ! 'ITAB_W_VXU_W'
   fldval = itab_w(i)%vxu_w(indp)
case(264) ! 'ITAB_W_VYU_W'
   fldval = itab_w(i)%vyu_w(indp)
case(265) ! 'ITAB_W_FARM'
   fldval = itab_w(i)%farm(indp)
case(266) ! 'ITAB_W_FARV'
   fldval = itab_w(i)%farv(indp)
case(267) ! 'ITAB_W_IWNUD'
   fldval = itab_w(i)%iwnud(indp)
case(268) ! 'ITAB_W_FNUD'
   fldval = itab_w(i)%fnud(indp)

! Time-averaged fields

case(269) ! 'RSHORT_AVG'

   if (.not. allocated(rshort_avg)) go to 1000

   fldval = rshort_avg(i)

case(270) ! 'RSHORTUP_AVG'

   if (.not. allocated(rshortup_avg)) go to 1000

   fldval = rshortup_avg(i)

case(271) ! 'RLONG_AVG'

   if (.not. allocated(rlong_avg)) go to 1000

   fldval = rlong_avg(i)

case(272) ! 'RLONGUP_AVG'

   if (.not. allocated(rlongup_avg)) go to 1000

   fldval = rlongup_avg(i)

case(273) ! 'RSHORT_TOP_AVG'

   if (.not. allocated(rshort_top_avg)) go to 1000

   fldval = rshort_top_avg(i)

case(274) ! 'RSHORTUP_TOP_AVG'

   if (.not. allocated(rshortup_top_avg)) go to 1000

   fldval = rshortup_top_avg(i)

case(275) ! 'RLONGUP_TOP_AVG'

   if (.not. allocated(rlongup_top_avg)) go to 1000

   fldval = rlongup_top_avg(i)

case(276) ! 'SENSFLUX_AVG'

   if (.not. allocated(sflux_t_avg)) go to 1000

   fldval = sflux_t_avg(i) * cp

case(277) ! 'LATFLUX_AVG'

   if (.not. allocated(sflux_r_avg)) go to 1000

   fldval = sflux_r_avg(i) * alvl

case(278) ! 'VAPFLUX_AVG'

   if (.not. allocated(sflux_r_avg)) go to 1000

   fldval = sflux_r_avg(i)

! Miscellaneous and new additions

case(279) ! 'RHO_OBS'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1))

case(280) ! 'THETA_OBS'

   if (.not. allocated(theta_obs)) go to 1000

   fldval = theta_obs(k,itab_w(i)%iwnud(1))

case(281) ! 'SHW_OBS'

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3

case(282) ! 'UZONAL_OBS'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1))

case(283) ! 'UMERID_OBS'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1))

case(284) ! 'RHO_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_sim(k,itab_w(i)%iwnud(1))

case(285) ! 'THETA_SIM'

   if (.not. allocated(theta_sim)) go to 1000

   fldval = theta_sim(k,itab_w(i)%iwnud(1))

case(286) ! 'SHW_SIM'

   if (.not. allocated(shw_sim)) go to 1000

   fldval = shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(287) ! 'UZONAL_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_sim(k,itab_w(i)%iwnud(1))

case(288) ! 'UMERID_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_sim(k,itab_w(i)%iwnud(1))

case(289) ! 'RHO_OBS_SIM'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1)) &
          - rho_sim(k,itab_w(i)%iwnud(1))

case(290) ! 'THETA_OBS_SIM'

   if (.not. allocated(theta_obs)) go to 1000

   fldval = theta_obs(k,itab_w(i)%iwnud(1)) &
          - theta_sim(k,itab_w(i)%iwnud(1))

case(291) ! 'SHW_OBS_SIM'

   if (.not. allocated(shw_obs)) go to 1000

   fldval = shw_obs(k,itab_w(i)%iwnud(1)) * 1.e3 &
          - shw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(292) ! 'UZONAL_OBS_SIM'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1)) &
          - uzonal_sim(k,itab_w(i)%iwnud(1))

case(293) ! 'UMERID_OBS_SIM'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1)) &
          - umerid_sim(k,itab_w(i)%iwnud(1))

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

