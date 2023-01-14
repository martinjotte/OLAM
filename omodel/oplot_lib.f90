subroutine oplot_lib(kk,ii,infotyp,fldname0,wtbot,wttop,fldval,notavail)

use mem_ijtabs,  only: itab_w, itab_v, itab_m, &
                       jtab_w, jtab_v, jtab_m, &
                       jtw_prog, jtv_prog, jtm_vadj
use mem_sfcg,    only: itab_wsfc, itab_msfc, itab_vsfc, sfcg
use mem_land,    only: land, mland, omland, nzg, slz, dslz, slzt
use mem_lake,    only: lake, mlake, omlake
use mem_sea,     only: sea,  msea,  omsea
use pom2k1d,     only: pom
use leaf4_soil,  only: soil_wat2pot

use mem_basic,   only: vmc, wmc, vc, wc, rho, press, &
                       thil, theta, tair, rr_w, rr_v, vxe, vye, vze
use mem_cuparm,  only: conprr, aconpr, qwcon, cbmf

use mem_grid,    only: mza, lpm, lpv, lpw, lsw, &
                       zt, dzt, dnu, dnv, arw0, arm0, arv, arw, volt, &
                       xem, yem, zem, &
                       xev, yev, zev, xew, yew, zew, &
                       topm, topw, glatm, glonm, glatv, glonv, glatw, glonw, &
                       vnx, vny, vnz, wnx, wny, wnz, &
                       wnxo2, wnyo2, wnzo2, dzt_bot

use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                       con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                       q2, q6, q7, ccntyp, con_gccn, con_ifn, &
                       pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                       accpd, accpr, accpp, accps, accpa, accpg, accph, &
                       cldnum

use micro_coms,  only: rxmin
use mem_co2,     only: rr_co2, co2_sh2ppm

use ccnbin_coms, only: nccntyp

use mem_radiate, only: fthrd_lw, fthrd_sw, rshort, rlong, rlongup, albedt, &
                       rshort_top, rshortup_top, rlongup_top
use mem_addsc,   only: addsc
use mem_tend,    only: vmxet, vmyet, vmzet
use mem_para,    only: myrank
use mem_turb,    only: vkm_sfc, sfluxt, sfluxr, pblh, vkm, vkh, ustar, wstar
use mem_nudge,   only: rho_obs, theta_obs, rrw_obs, uzonal_obs, umerid_obs, &
                       rho_sim, theta_sim, rrw_sim, uzonal_sim, umerid_sim
use therm_lib,   only: qtk, qwtk, rhovsl_inv
use misc_coms,   only: io6, naddsc, mdomain
use oplot_coms,  only: op
use consts_coms, only: p00i, rocp, erad, eradi, piu180, cp, alvl, grav, omega2, eps_virt
use mem_flux_accum, only:     rshort_accum,         rshortup_accum, &
                               rlong_accum,          rlongup_accum, &
                          rshort_top_accum,     rshortup_top_accum, &
                         rlongup_top_accum,                         &
                                  vc_accum,               wc_accum, &
                               press_accum,             tair_accum, &
                                rr_v_accum,      latheat_liq_accum, &
                         latheat_ice_accum,             vels_accum, &
                             airtemp_accum,           airrrv_accum, &
                             cantemp_accum,           canrrv_accum, &
                              sfluxt_accum,           sfluxr_accum, & ! ---/nud
                                 pcp_accum,           runoff_accum, & ! nud/---
                             sfctemp_accum,          fracliq_accum, & ! nud/nud
                            skintemp_accum,           wxferi_accum, &
                              wxferp_accum,           wxfer1_accum, &
                        airtemp_dmin_accum,     airtemp_dmax_accum, &
                        cantemp_dmin_accum,     cantemp_dmax_accum, &
                        vegtemp_dmin_accum,     vegtemp_dmax_accum, &
                       soiltemp_dmin_accum,    soiltemp_dmax_accum

  use mem_average_vars, only: nz_avg, &
         press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
         tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
      press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
      airtempk_davg,  airtempk_dmin,  airtempk_dmax, &
      cantempk_davg,  cantempk_dmin,  cantempk_dmax, &
      vegtempk_davg,  vegtempk_dmin,  vegtempk_dmax, &
     soiltempk_davg, soiltempk_dmin, soiltempk_dmax, &
        sfluxt_davg,    sfluxr_davg

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
             vc_accum_prev0,           vc_accum_prev1, &
             wc_accum_prev0,           wc_accum_prev1, &
          press_accum_prev0,        press_accum_prev1, &
           tair_accum_prev0,         tair_accum_prev1, &
           rr_v_accum_prev0,         rr_v_accum_prev1, &
    latheat_liq_accum_prev0,  latheat_liq_accum_prev1, &
    latheat_ice_accum_prev0,  latheat_ice_accum_prev1, &
           vels_accum_prev0,         vels_accum_prev1,     vels_accum_prev2,     vels_accum_prev3, &
        airtemp_accum_prev0,      airtemp_accum_prev1,  airtemp_accum_prev2,  airtemp_accum_prev3, &
         airrrv_accum_prev0,       airrrv_accum_prev1,   airrrv_accum_prev2,   airrrv_accum_prev3, &
        cantemp_accum_prev0,      cantemp_accum_prev1,  cantemp_accum_prev2,  cantemp_accum_prev3, &
         canrrv_accum_prev0,       canrrv_accum_prev1,   canrrv_accum_prev2,   canrrv_accum_prev3, &
       skintemp_accum_prev0,     skintemp_accum_prev1, skintemp_accum_prev2, skintemp_accum_prev3, &
         sfluxt_accum_prev0,       sfluxt_accum_prev1,   sfluxt_accum_prev2,   sfluxt_accum_prev3, &
         sfluxr_accum_prev0,       sfluxr_accum_prev1,   sfluxr_accum_prev2,   sfluxr_accum_prev3, & ! fast can nud
            pcp_accum_prev0,          pcp_accum_prev1,      pcp_accum_prev2,      pcp_accum_prev3, & ! fast can nud
        sfctemp_accum_prev0,      sfctemp_accum_prev1,                                             & ! fast can nud
        fracliq_accum_prev0,      fracliq_accum_prev1,                                             & ! fast can nud
         runoff_accum_prev0,       runoff_accum_prev1,   runoff_accum_prev2,   runoff_accum_prev3, &
         wxferi_accum_prev0,       wxferi_accum_prev1,   wxferi_accum_prev2,   wxferi_accum_prev3, &
         wxferp_accum_prev0,       wxferp_accum_prev1,   wxferp_accum_prev2,   wxferp_accum_prev3, &
         wxfer1_accum_prev0,       wxfer1_accum_prev1,   wxfer1_accum_prev2,   wxfer1_accum_prev3, &
       soil_water_tot_prev0,     soil_water_tot_prev1,                                             &
            head_wtab_prev0,          head_wtab_prev1,                                             &

   airtemp_dmin_accum_prev0,  airtemp_dmin_accum_prev1,  airtemp_dmin_accum_prev2,  airtemp_dmin_accum_prev3, &
   airtemp_dmax_accum_prev0,  airtemp_dmax_accum_prev1,  airtemp_dmax_accum_prev2,  airtemp_dmax_accum_prev3, &
   cantemp_dmin_accum_prev0,  cantemp_dmin_accum_prev1,  cantemp_dmin_accum_prev2,  cantemp_dmin_accum_prev3, &
   cantemp_dmax_accum_prev0,  cantemp_dmax_accum_prev1,  cantemp_dmax_accum_prev2,  cantemp_dmax_accum_prev3, &
   vegtemp_dmin_accum_prev0,  vegtemp_dmin_accum_prev1,  vegtemp_dmin_accum_prev2,  vegtemp_dmin_accum_prev3, &
   vegtemp_dmax_accum_prev0,  vegtemp_dmax_accum_prev1,  vegtemp_dmax_accum_prev2,  vegtemp_dmax_accum_prev3, &
  soiltemp_dmin_accum_prev0, soiltemp_dmin_accum_prev1, soiltemp_dmin_accum_prev2, soiltemp_dmin_accum_prev3, &
  soiltemp_dmax_accum_prev0, soiltemp_dmax_accum_prev1, soiltemp_dmax_accum_prev2, soiltemp_dmax_accum_prev3

use mem_sfcnud, only: sfcwat_nud, sfctemp_nud, fracliq_nud ! fast can nud

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

integer :: klev,nls,jv,im,iv,iw,kp,k,i
real :: raxis,u,v,rpolyi,raxisi
real :: vx, vy, vz
real :: tempk, fracliq
real :: vc_change, wc_change
integer :: iw1,iw2,iland,ilake,isea,iwsfc
integer :: npoly,j
integer :: lenstr, ic, ifield
integer, save :: indp, icase
integer :: jasfc
real :: area_sum
real :: zobs, press_zobs, exner_zobs, wind_zobs, theta_zobs, rrv_zobs
real :: canexner, cantheta, canthetav, airthetav, tstar, rstar, ufree

real :: vcc
real :: vcc_init
real :: accpboth_prev0, accpboth_prev1, accpboth_prev2, accpboth_prev3
real :: denom, fldval1, fldval2

real :: head(nzg)
real :: psi, psi_slope

real :: zanal_swtc5, zanal0_swtc5

real, parameter :: onethird = 1./3.

integer, parameter :: nfields = 803
character(len=40) :: fldlib(4,nfields)
character(len=40), save :: fldname

!  fldname     stagpt/dimens     field description & units                   field #
!-----------------------------------------------------------------------------------
! ATMOSPHERE - 3D

data fldlib(1:4,  1:36)/ &
 'VMC'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  1
 'WMC'           ,'W3' ,'W MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'            ,& !p  2
 'VMP'           ,'V3' ,'V-NORMAL MOMENTUM',' (kg m:S2:-2   s:S2:-1  )'     ,& !p  3
 'VC'            ,'V3' ,'V-NORMAL VELOCITY',' (m s:S2:-1  )'                ,& !p  4
 'WC'            ,'W3' ,'W',' (m s:S2:-1  )'                                ,& !p  5
 'RHO'           ,'T3' ,'AIR DENSITY',' (kg m:S2:-3  )'                     ,& !p  6
 'PRESS'         ,'T3' ,'PRESSURE',' (hPa)'                                 ,& !   7
 'THIL'          ,'T3' ,'ICE-LIQUID THETA',' (K)'                           ,& !p  8
 'THETA'         ,'T3' ,'THETA',' (K)'                                      ,& !p  9
 'AIRTEMPK'      ,'T3' ,'AIR TEMP',' (K)'                                   ,& !p 10
 'AIRTEMPC'      ,'T3' ,'AIR TEMP',' (C)'                                   ,& !p 11
 'RR_W'          ,'T3' ,'TOTAL WATER',' (g kg:S2:-1  )'                     ,& !p 12
 'RR_V'          ,'T3' ,'WATER VAPOR',' (g kg:S2:-1  )'                     ,& !p 13
 'RR_C'          ,'T3' ,'CLOUD',' (g kg:S2:-1  )'                           ,& !p 14
 'RR_D'          ,'T3' ,'DRIZZLE',' (g kg:S2:-1  )'                         ,& !p 15
 'RR_R'          ,'T3' ,'RAIN',' (g kg:S2:-1  )'                            ,& !p 16
 'RR_P'          ,'T3' ,'PRIS ICE',' (g kg:S2:-1  )'                        ,& !p 17
 'RR_S'          ,'T3' ,'SNOW',' (g kg:S2:-1  )'                            ,& !p 18
 'RR_A'          ,'T3' ,'AGGREGATES',' (g kg:S2:-1  )'                      ,& !p 19
 'RR_G'          ,'T3' ,'GRAUPEL',' (g kg:S2:-1  )'                         ,& !p 20
 'RR_H'          ,'T3' ,'HAIL',' (g kg:S2:-1  )'                            ,& !p 21
 'RR_CP'         ,'T3' ,'CLOUD + PRIS ICE',' (g kg:S2:-1  )'                ,& !p 22
 'RR_TOTLIQ'     ,'T3' ,'LIQUID',' (g kg:S2:-1  )'                          ,& !p 23
 'RR_TOTICE'     ,'T3' ,'ICE',' (g kg:S2:-1  )'                             ,& !p 24
 'RR_TOTCOND'    ,'T3' ,'CONDENSATE',' (g kg:S2:-1  )'                      ,& !p 25
 'CON_C'         ,'T3' ,'CLOUD NUM',' (# mg:S2:-1  )'                       ,& !p 26
 'CON_D'         ,'T3' ,'DRIZZLE NUM',' (# g:S2:-1  )'                      ,& !p 27
 'CON_R'         ,'T3' ,'RAIN NUM',' (# kg:S2:-1  )'                        ,& !p 28
 'CON_P'         ,'T3' ,'PRIS ICE NUM',' (# mg:S2:-1  )'                    ,& !p 29
 'CON_S'         ,'T3' ,'SNOW NUM',' (# kg:S2:-1  )'                        ,& !p 30
 'CON_A'         ,'T3' ,'AGGREGATES NUM',' (# kg:S2:-1  )'                  ,& !p 31
 'CON_G'         ,'T3' ,'GRAUPEL NUM',' (# kg:S2:-1  )'                     ,& !p 32
 'CON_H'         ,'T3' ,'HAIL NUM',' (# kg:S2:-1  )'                        ,& !p 33
 'CON_CCN'       ,'T3' ,'CCN NUM',' (# mg:S2:-1  )'                         ,& !p 34
 'CON_GCCN'      ,'T3' ,'GCCN NUM',' (# g:S2:-1  )'                         ,& !p 35
 'CON_IFN'       ,'T3' ,'IFN NUM',' (# g:S2:-1  )'                           / !p 36

data fldlib(1:4, 37:61)/ &
 'VKM'           ,'T3' ,'VERT TURB MOMENTUM K',' (N s m:S2:-2  )'           ,& !p 37
 'FTHRD'         ,'T3' ,'RADIATIVE THETA TENDENCY',' (K s:S2:-1  )'         ,& !p 38
 'SPEEDW'        ,'T3' ,'WIND SPEED AT W',' (m s:S2:-1  )'                  ,& !p 39
 'AZIMW'         ,'T3' ,'WIND AZIMUTH AT W',' (deg)'                        ,& !p 40
 'ZONAL_WINDW'   ,'T3' ,'ZONAL WIND AT W',' (m s:S2:-1  )'                  ,& !p 41
 'MERID_WINDW'   ,'T3' ,'MERIDIONAL WIND AT W',' (m s:S2:-1  )'             ,& !p 42
 'RVORTZM'       ,'P3' ,'REL VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 43
 'TVORTZM'       ,'P3' ,'TOT VERT VORTICITY AT M',' (s:S2:-1  )'            ,& !p 44
 'RVORTZM_P'     ,'P3' ,'REL VERT VORTICITY PERT AT M',' (s:S2:-1  )'       ,& !p 45
 'DIVERG'        ,'T3' ,'HORIZONTAL DIVERGENCE',' (s:S2:-1  )'              ,& !p 46
 'VMASSFLUX'     ,'V3' ,'GRID CELL V-FACE MASS FLUX',' (kg s:S2:-1  )'      ,& !  47
 'VC_P'          ,'V3' ,'NORMAL WIND PERT AT V',' (m s:S2:-1  )'            ,& !p 48
 'PRESS_P'       ,'T3' ,'PRESSURE PERT',' (hPa)'                            ,& !  49
 'RHO_P'         ,'T3' ,'DENSITY PERT',' (kg m:S2:-3  )'                    ,& !  50
 'THETA_P'       ,'T3' ,'THETA PERT',' (K)'                                 ,& !  51
 'AIRTEMPK_P'    ,'T3' ,'AIR TEMP PERT',' (K)'                              ,& !p 52
 'VMT'           ,'V3' ,'V-NORM MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'  ,& !  53
 'WMT'           ,'W3' ,'W MOMENTUM TEND',' (kg m:S2:-2   s:S2:-2  )'       ,& !  54
 'ADDSC'         ,'T3' ,'ADDED SCALAR AMOUNT PER KG AIR','( )'              ,& !p 55
 'ADDSCP'        ,'T3' ,'SCALAR PERTURBATION',' ( )'                        ,& !  56
 'ZPLEV'         ,'T3' ,'HEIGHT OF CONST P SFC',' (m)'                      ,& !p 57
 'QWCON'         ,'T3' ,'CUPARM CONDENSATE MIX RATIO',' (g kg:S2:-1  )'     ,& !p 58
 'CO2CON'        ,'T3' ,'CO2 CONCENTRATION',' (ppmV dry air)'               ,& !p 59
 'CO2PERT'       ,'T3' ,'CO2 CHANGE',' (ppmV dry air)'                      ,& !p 60
 'RH_LIQ'        ,'T3' ,'RH_LIQUID',' (%)'                                   / !p 61

! ATMOSPHERE - 2D

data fldlib(1:4, 62:98)/ &
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
 'PCPRBOTH'      ,'T2' ,'MIC+CONV PCP RATE',' (kg m:S2:-2   h:S2:-1  )'     ,& !  84
 'ACCPD'         ,'T2' ,'ACCUM DRIZZLE',' (kg m:S2:-2  )'                   ,& !  85
 'ACCPR'         ,'T2' ,'ACCUM RAIN',' (kg m:S2:-2  )'                      ,& !  86
 'ACCPP'         ,'T2' ,'ACCUM PRIST ICE',' (kg m:S2:-2  )'                 ,& !  87
 'ACCPS'         ,'T2' ,'ACCUM SNOW',' (kg m:S2:-2  )'                      ,& !  88
 'ACCPA'         ,'T2' ,'ACCUM AGGREGATES',' (kg m:S2:-2  )'                ,& !  89
 'ACCPG'         ,'T2' ,'ACCUM GRAUPEL',' (kg m:S2:-2  )'                   ,& !  90
 'ACCPH'         ,'T2' ,'ACCUM HAIL',' (kg m:S2:-2  )'                      ,& !  91
 'ACCPMIC'       ,'T2' ,'ACCUM MIC PCP',' (kg m:S2:-2  )'                   ,& !  92
 'ACCPCON'       ,'T2' ,'ACCUM CONV PCP',' (kg m:S2:-2  )'                  ,& !  93
 'ACCPBOTH'      ,'T2' ,'ACCUM MIC+CONV PCP',' (kg m:S2:-2  )'              ,& !  94
 'WSTAR'         ,'T2' ,'PBL CONVECTIVE VELOCITY',' (m s:S2:-1  )'          ,& !  95
 'PSFC'          ,'T2' ,'SURFACE PRESSURE',' (hPa)'                         ,& !  96
 'PMSL'          ,'T2' ,'SEA LEVEL PRESSURE',' (hPa)'                       ,& !  97
 'CBMF'          ,'T2' ,'CUPARM CLOUD BASE MASS FLUX','(kg m:S2:-2 s:S2:-1 )'/ !  98

! ATMOSPHERE DIF2 fields (3D & 2D)

data fldlib(1:4,101:118)/ &
 'ZONAL_WINDW_DIF2' ,'T3' ,'ATM ZONAL VELOCITY DIF2',' (m s:S2:-1  )'       ,& ! 101
 'MERID_WINDW_DIF2' ,'T3' ,'ATM MERID VELOCITY DIF2',' (m s:S2:-1  )'       ,& ! 102
 'WC_DIF2'          ,'T3' ,'ATM VERT VELOCITY DIF2',' (m s:S2:-1  )'        ,& ! 103
 'PRESS_DIF2'       ,'T3' ,'ATM PRESSURE DIF2',' (hPa)'                     ,& ! 104
 'AIRTEMPK_DIF2'    ,'T3' ,'ATM TEMP DIF2',' (K)'                           ,& ! 105
 'RR_V_DIF2'        ,'T3' ,'ATM VAP MIX RATIO DIF2',' (g kg:S2:-1  )'       ,& ! 106
 'LATHEAT_LIQ_DIF2' ,'T3' ,'ATM LATENT HEAT LIQ DIF2',' (K hr:S2:-1  )'     ,& ! 107
 'LATHEAT_ICE_DIF2' ,'T3' ,'ATM LATENT HEAT ICE DIF2',' (K hr:S2:-1  )'     ,& ! 108
 'PCPMIC_DIF2'      ,'T2' ,'MICPHYS PRECIP DIF2',' (mm/day)'                ,& ! 109
 'PCPCON_DIF2'      ,'T2' ,'CONV PRECIP DIF2',' (mm/day)'                   ,& ! 110
 'PCPBOTH_DIF2'     ,'T2' ,'MICPHYS + CONV PRECIP DIF2',' (mm/day)'         ,& ! 111
 'RSHORT_DIF2'      ,'T2' ,'SFC DOWN SHORTWV FLX DIF2',' (W m:S2:-2  )'     ,& ! 112
 'RSHORTUP_DIF2'    ,'T2' ,'SFC UP SHORTWV FLX DIF2',' (W m:S2:-2  )'       ,& ! 113
 'RLONG_DIF2'       ,'T2' ,'SFC DOWN LONGWV FLX DIF2',' (W m:S2:-2  )'      ,& ! 114
 'RLONGUP_DIF2'     ,'T2' ,'SFC UP LONGWV FLX DIF2',' (W m:S2:-2  )'        ,& ! 115
 'RSHORT_TOP_DIF2'  ,'T2' ,'TOP DOWN SHORTWV FLX DIF2',' (W m:S2:-2  )'     ,& ! 116
 'RSHORTUP_TOP_DIF2','T2' ,'TOP UP SHORTWV FLX DIF2',' (W m:S2:-2  )'       ,& ! 117
 'RLONGUP_TOP_DIF2' ,'T2' ,'TOP UP LONGWV FLX DIF2',' (W m:S2:-2  )'         / ! 118

! ATMOSPHERE DIF4 fields (2D)

data fldlib(1:4,129:134)/ &
 'PCPMIC_DIF4'    ,'T2' ,'MICPHYS PRECIP DIFF4',' (mm/day)'                  ,& ! 129
 'PCPCON_DIF4'    ,'T2' ,'CONV PRECIP DIFF4',' (mm/day)'                     ,& ! 130
 'PCPBOTH_DIF4'   ,'T2' ,'MICPHYS + CONV PRECIP DIFF4',' (mm/day)'           ,& ! 131
 'PCPMIC_REL4'    ,'T2' ,'MICPHYS PRECIP RELATIVE DIFF4',' '                 ,& ! 132
 'PCPCON_REL4'    ,'T2' ,'CONV PRECIP RELATIVE DIFF4',' '                    ,& ! 133
 'PCPBOTH_REL4'   ,'T2' ,'MICPHYS + CONV PRECIP RELATIVE DIFF4',' '           / ! 134

! LAND_CELLS - 3D

data fldlib(1:4,181:203)/ &
 'SAND'              ,'L3G','SOIL SAND FRACTION',' ( )'                     ,& ! 181
 'CLAY'              ,'L3G','SOIL CLAY FRACTION',' ( )'                     ,& ! 182
 'SILT'              ,'L3G','SOIL SILT FRACTION',' ( )'                     ,& ! 183
 'ORGAN'             ,'L3G','SOIL ORGANIC SPECIFIC DENSITY',' (g/kg )'      ,& ! 184
 'BULKDENS_DRYSOIL'  ,'L3G','SOIL DRY BULK DENSITY',' (kg m:S2:-3  )'       ,& ! 185
 'PH_SOIL'           ,'L3G','SOIL PH',' ( )'                                ,& ! 186
 'CEC_SOIL'          ,'L3G','SOIL CATION EXCHANGE CAPACITY',' (cmol+/kg )'  ,& ! 187
 'WRESID_VG'         ,'L3G','SOIL RESIDUAL WATER CONTENT (VG)',' ( )'       ,& ! 188
 'WSAT_VG'           ,'L3G','SOIL SATURATION WATER CONTENT (VG)',' ( )'     ,& ! 189
 'KSAT_VG'           ,'L3G','SOIL SATURATION HYD CONDUCTIVITY (VG)',' (m/s)',& ! 190
 'ALPHA_VG'          ,'L3G','SOIL ALPHA_VG',' (1/m)'                        ,& ! 191
 'EN_VG'             ,'L3G','SOIL EN_VG',' ( )'                             ,& ! 192
 'LAMBDA_VG'         ,'L3G','SOIL LAMBDA_VG',' ( )'                         ,& ! 193
 'SPECIFHEAT_DRYSOIL','L3G','SOIL DRY SPECIFIC HEAT',' ( )'                 ,& ! 194
 'SOIL_ENERGY'       ,'L3G','SOIL ENERGY',' (J cm:S2:-3  )'                 ,& ! 195
 'SOIL_TEMPK'        ,'L3G','SOIL TEMP',' (K)'                              ,& ! 196
 'SOIL_FRACLIQ'      ,'L3G','LIQUID FRACTION OF SOIL WATER',' ( )'          ,& ! 197
 'SOIL_WATER'        ,'L3G','SOIL WATER CONTENT',' ( )'                     ,& ! 198
 'SFWAT_MASS'        ,'L3S','SFCWATER MASS',' (kg m:S2:-2  )'               ,& ! 199
 'SFWAT_ENERGY'      ,'L3S','SFCWATER ENERGY',' (J g:S2:-1  )'              ,& ! 200
 'SFWAT_TEMPK'       ,'L3S','SFCWATER TEMP',' (K)'                          ,& ! 201
 'SFWAT_FRACLIQ'     ,'L3S','SFCWATER LIQUID FRACTION',' ( )'               ,& ! 202
 'SFWAT_DEPTH'       ,'L3S','SFCWATER DEPTH',' (m)'                          / ! 203

! LAND_CELLS - 2D

data fldlib(1:4,205:224)/ &
 'USDA_TEXT'        ,'L2' ,'USDA SOIL TEXTURAL CLASS (2D)',' ( )'            ,& ! 205
 'Z_BEDROCK'        ,'L2' ,'Z_BEDROCK',' (m)'                                ,& ! 206
 'GPP'              ,'L2' ,'GROSS PRIMARY PRODUCTION OF CARBON',' ()'        ,& ! 207
 'GLHYMPS_KSAT'     ,'L2' ,'GLHYMPS KSAT',' (m/s)'                           ,& ! 208
 'GLHYMPS_KSAT_PFR' ,'L2' ,'GLHYMPS KSAT WITH PERMAFROST',' (m/s)'           ,& ! 209
 'GLHYMPS_POROS'    ,'L2' ,'GLHYMPS POROSITY',' ( )'                         ,& ! 208
 'NLEV_SFWAT'       ,'L2' ,'NUMBER OF SFCWATER LAYERS',' ( )'                ,& ! 211
 'VEG_NDVIC'        ,'L2' ,'VEGETATION NDVI',' ( )'                          ,& ! 212
 'VEG_TEMPC'        ,'L2' ,'VEGETATION TEMP',' (C)'                          ,& ! 213
 'VEG_TEMPK'        ,'L2' ,'VEGETATION TEMP',' (K)'                          ,& ! 214
 'VEG_WATER'        ,'L2' ,'VEGETATION SFC WATER ',' (kg m:S2:-2  )'         ,& ! 215
 'STOM_RESIST'      ,'L2' ,'STOMATAL RESISTANCE',' (s m:S2:-1  )'            ,& ! 216
 'SFCWATER_TOT'     ,'L2' ,'TOTAL SFCWATER MASS',' (kg m:S2:-2  )'           ,& ! 217
 'SFCWATER_TOP_TEMP','L2' ,'SFCWATER TOPLAYER TEMP',' (K)'                   ,& ! 218
 'SOIL_TOP_TEMP'    ,'L2' ,'SOIL TOPLAYER TEMP',' (K)'                       ,& ! 219
 'GROUND_RRV'       ,'L2' ,'EQUIL VAP MIX RATIO OF SOIL',' (g kg:S2:-1  )'   ,& ! 220
 'SOIL_DEPTH'       ,'L2' ,'SOIL DEPTH',' (m)'                               ,& ! 221
 'SOIL_WATER_TOT'   ,'L2' ,'TOTAL SOIL WATER',' (m)'                         ,& ! 222
 'HEAD0'            ,'L2' ,'HEAD0',' (m)'                                    ,& ! 223
 'SLOPE_FACT'       ,'L2' ,'SUBGRID OROGRAPHY SLOPE FACTOR',' ( )'            / ! 224


! LAND_CELLS - DIF2 fields

data fldlib(1:4,230:250)/ &
 'WXFERI_DIF2'         ,'L2' ,'INFILTRATION FLUX DIF2',' (mm/day)'       ,& ! 230
 'WXFERP_DIF2'         ,'L2' ,'PERCOLATION FLUX DIF2',' (mm/day)'        ,& ! 231
 'WXFER1_DIF2'         ,'L2' ,'SOIL BOTTOM WATER FLUX DIF2',' (mm/day)'  ,& ! 232
 'SOIL_WATER_TOT_DIF2' ,'L2' ,'TOTAL SOIL WATER DIF2',' (m)'             ,& ! 233
 'VEGTEMPK_DMIN_DIF2'  ,'L2' ,'VEG DMIN TEMP DIF2',' (K)'                ,& ! 234
 'VEGTEMPK_DMAX_DIF2'  ,'L2' ,'VEG DMAX TEMP DIF2',' (K)'                ,& ! 235
 'VEGTEMPK_DSPAN_DIF2' ,'L2' ,'VEG DSPAN TEMP DIF2',' (K)'               ,& ! 236
 'SOILTEMPK_DMIN_DIF2' ,'L2' ,'SOIL DMIN TEMP DIF2',' (K)'               ,& ! 237
 'SOILTEMPK_DMAX_DIF2' ,'L2' ,'SOIL DMAX TEMP DIF2',' (K)'               ,& ! 238
 'SOILTEMPK_DSPAN_DIF2','L2' ,'SOIL DSPAN TEMP DIF2',' (K)'              ,& ! 239
 'WXFERIF_DIF2'        ,'L2' ,'INFILTRATION FRACTION DIF2',' ()'         ,& ! 240

! LAND_CELLS - DIF4 fields

 'WXFERI_DIF4'         ,'L2' ,'INFILTRATION FLUX DIF4',' (mm/day)'       ,& ! 241
 'WXFERP_DIF4'         ,'L2' ,'PERCOLATION FLUX DIF4',' (mm/day)'        ,& ! 242
 'WXFER1_DIF4'         ,'L2' ,'SOIL BOTTOM WATER FLUX DIF4',' (mm/day)'  ,& ! 243
 'VEGTEMPK_DMIN_DIF4'  ,'L2' ,'VEG DMIN TEMP DIF4',' (K)'                ,& ! 244
 'VEGTEMPK_DMAX_DIF4'  ,'L2' ,'VEG DMAX TEMP DIF4',' (K)'                ,& ! 245
 'VEGTEMPK_DSPAN_DIF4' ,'L2' ,'VEG DSPAN TEMP DIF4',' (K)'               ,& ! 246
 'SOILTEMPK_DMIN_DIF4' ,'L2' ,'SOIL DMIN TEMP DIF4',' (K)'               ,& ! 247
 'SOILTEMPK_DMAX_DIF4' ,'L2' ,'SOIL DMAX TEMP DIF4',' (K)'               ,& ! 248
 'SOILTEMPK_DSPAN_DIF4','L2' ,'SOIL DSPAN TEMP DIF4',' (K)'              ,& ! 249
 'WXFERIF_DIF4'        ,'L2' ,'INFILTRATION FRACTION DIF4',' ()'          / ! 250

! LAND_CELLS - ATM averages

data fldlib(1:4,251:280)/ &
 'AL_SFCWATER_TOT'        ,'T2' ,'AL TOTAL SFCWATER MASS',' (kg m:S2:-2  )'  ,& ! 251
 'AL_SOIL_WATER_TOT'      ,'T2' ,'AL TOTAL SOIL WATER',' (m)'                ,& ! 252

! LAND_CELLS - ATM averages of DIF2 fields

 'AL_WXFERI_DIF2'         ,'T2' ,'AL INFILTRATION FLUX DIF2',' (mm/day)'      ,& ! 253
 'AL_WXFERP_DIF2'         ,'T2' ,'AL PERCOLATION FLUX DIF2',' (mm/day)'       ,& ! 254
 'AL_WXFER1_DIF2'         ,'T2' ,'AL SOIL BOTTOM WATER FLUX DIF2',' (mm/day)' ,& ! 255
 'AL_SOIL_WATER_TOT_DIF2' ,'T2' ,'AL TOTAL SOIL WATER DIF2',' (m)'            ,& ! 256
 'AL_VEGTEMPK_DMIN_DIF2'  ,'T2' ,'AL VEG DMIN TEMP DIF2',' (K)'               ,& ! 257
 'AL_VEGTEMPK_DMAX_DIF2'  ,'T2' ,'AL VEG DMAX TEMP DIF2',' (K)'               ,& ! 258
 'AL_VEGTEMPK_DSPAN_DIF2' ,'T2' ,'AL VEG DSPAN TEMP DIF2',' (K)'              ,& ! 259
 'AL_SOILTEMPK_DMIN_DIF2' ,'T2' ,'AL SOIL DMIN TEMP DIF2',' (K)'              ,& ! 260
 'AL_SOILTEMPK_DMAX_DIF2' ,'T2' ,'AL SOIL DMAX TEMP DIF2',' (K)'              ,& ! 261
 'AL_SOILTEMPK_DSPAN_DIF2','T2' ,'AL SOIL DSPAN TEMP DIF2',' (K)'             ,& ! 262

! LAND_CELLS - ATM averages of DIF4 fields

 'AL_WXFERI_DIF4'         ,'T2' ,'AL INFILTRATION FLUX DIF4',' (mm/day)'      ,& ! 263
 'AL_WXFERP_DIF4'         ,'T2' ,'AL PERCOLATION FLUX DIF4',' (mm/day)'       ,& ! 264
 'AL_WXFER1_DIF4'         ,'T2' ,'AL SOIL BOTTOM WATER FLUX DIF4',' (mm/day)' ,& ! 265
 'AL_VEGTEMPK_DMIN_DIF4'  ,'T2' ,'AL VEG DMIN TEMP DIF4',' (K)'               ,& ! 266
 'AL_VEGTEMPK_DMAX_DIF4'  ,'T2' ,'AL VEG DMAX TEMP DIF4',' (K)'               ,& ! 267
 'AL_VEGTEMPK_DSPAN_DIF4' ,'T2' ,'AL VEG DSPAN TEMP DIF4',' (K)'              ,& ! 268
 'AL_SOILTEMPK_DMIN_DIF4' ,'T2' ,'AL SOIL DMIN TEMP DIF4',' (K)'              ,& ! 269
 'AL_SOILTEMPK_DMAX_DIF4' ,'T2' ,'AL SOIL DMAX TEMP DIF4',' (K)'              ,& ! 270
 'AL_SOILTEMPK_DSPAN_DIF4','T2' ,'AL SOIL DSPAN TEMP DIF4',' (K)'             ,& ! 271

! SEA_CELLS - 2D

 'SEATP'         ,'S2' ,'SEA SFC TEMP (PAST DATA)',' (K)'                    ,& ! 272
 'SEATF'         ,'S2' ,'SEA SFC TEMP (FUTURE DATA)',' (K)'                  ,& ! 273
 'SEATC'         ,'S2' ,'SEA SFC TEMP (CURRENT)',' (K)'                      ,& ! 274
 'SEAICEP'       ,'S2' ,'SEAICE FRACTION (PAST DATA)',' ( )'                 ,& ! 275
 'SEAICEF'       ,'S2' ,'SEAICE FRACTION (FUTURE DATA)',' ( )'               ,& ! 276
 'SEAICEC'       ,'S2' ,'SEAICE FRACTION (CURRENT)',' ( )'                   ,& ! 277
 'WDEPTH'        ,'S2' ,'SEA WDEPTH',' (m)'                                  ,& ! 278
 'POM_KBA'       ,'S2' ,'POM LEVELS',' (#)'                                  ,& ! 279
 'POM_TEMPSFC'   ,'S2' ,'POM SURFACE TEMP',' (K)'                             / ! 280

! SFC GRID CELLS - 2D

data fldlib(1:4,297:332)/ &
 'LEAF_CLASS'         ,'C2' ,'LEAF CLASS',' ( )'                               ,& ! 297
 'SFCG_AREA'          ,'C2' ,'SFCG CELL AREA',' (m:S2:2  )'                    ,& ! 298
 'SFCG_GLATW'         ,'C2' ,'SFCG CELL LATITUDE',' (deg)'                     ,& ! 299
 'SFCG_GLONW'         ,'C2' ,'SFCG CELL LONGITUDE',' (deg)'                    ,& ! 300
 'SFCG_TOPW'          ,'C2' ,'SFCG CELL TOPW',' (m)'                           ,& ! 301
 'SFCG_ROUGH'         ,'C2' ,'SFCG NET ROUGHNESS HEIGHT',' (m)'                ,& ! 302
 'SFCG_VELS'          ,'C2' ,'SFCG WIND SPEED',' (m s:S2:-1  )'                ,& ! 303
 'SFCG_PRSS'          ,'C2' ,'SFCG PRESSURE',' (hPa)'                          ,& ! 304
 'SFCG_RHOS'          ,'C2' ,'SFCG DENSITY',' (kg m:S2:-3  )'                  ,& ! 305
 'SFCG_AIRTEMPK'      ,'C2' ,'SFCG ATM TEMP',' (K)'                            ,& ! 306
 'SFCG_AIRRRV'        ,'C2' ,'SFCG ATM RRV',' (g kg:S2:-1  )'                  ,& ! 307
 'SFCG_CANTEMPK'      ,'C2' ,'SFCG CANOPY AIR TEMP',' (K)'                     ,& ! 308
 'SFCG_CANRRV'        ,'C2' ,'SFCG CANOPY VAPOR MIX RATIO',' (g kg:S2:-1  )'   ,& ! 309
 'SFCG_SKINTEMPK'     ,'C2' ,'SFCG SKIN TEMP',' (K)'                           ,& ! 310
 'SFCG_GSS_SRRV'      ,'C2' ,'SFCG SAT VAP MIX RATIO',' (g kg:S2:-1  )'        ,& ! 311
 'HEAD1'              ,'C2' ,'WATER SFC HEAD',' (m)'                           ,& ! 312
 'HEAD_WTAB'          ,'C2' ,'HEAD AT WATER TABLE',' (m)'                      ,& ! 313
 'SFCG_SENSFLUX'      ,'C2','SFCG SENS HEAT FLX',' (W m:S2:-2  )'              ,& ! 314
 'SFCG_LATFLUX'       ,'C2','SFCG LAT HEAT FLX',' (W m:S2:-2  )'               ,& ! 315
 'SFCG_VAPFLUX'       ,'C2','SFCG VAP FLX',' (kg m:S2:-2   s:S2:-1  )'         ,& ! 316
 'SFCG_SPEED10M'      ,'C2','SFCG 10M WIND SPEED',' (m s:S2:-1  )'             ,& ! 317
 'SFCG_SPEED2M'       ,'C2','SFCG 2M WIND SPEED',' (m s:S2:-1  )'              ,& ! 318
 'SFCG_TEMPK2M'       ,'C2','SFCG 2M TEMP',' (K)'                              ,& ! 319
 'SFCG_RVAP2M'        ,'C2','SFCG 2M RRV',' (g kg:S2:-1  )'                    ,& ! 320
 'SFCG_RSHORT'        ,'C2','SFCG DOWN SW FLX',' (W m:S2:-2  )'                ,& ! 321
 'SFCG_RLONG'         ,'C2','SFCG DOWN LW FLX',' (W m:S2:-2  )'                ,& ! 322
 'SFCG_RLONGUP'       ,'C2','SFCG UP LW FLX',' (W m:S2:-2  )'                  ,& ! 323
 'SFCG_RLONG_ALBEDO'  ,'C2','SFCG NET LW ALBEDO',' ( )'                        ,& ! 324
 'SFCG_ALBEDO_BEAM'   ,'C2','SFCG NET BEAM ALBEDO',' ( )'                      ,& ! 325
 'SFCG_ALBEDO_DIFFUSE','C2','SFCG NET DIFFUSE ALBEDO',' ( )'                   ,& ! 326
 'SFCG_BATHYM'        ,'C2','SFCG CELL BATHYMETRY',' (m)'                      ,& ! 327
 'SFCG_PCPG'          ,'C2','SFCG PCPG',' (kg m:S2:-2  )'                      ,& ! 328
 'SFCG_VC'            ,'B2','SFCG VC',' (m s:S2:-1  )'                         ,& ! 329
 'WAT_DEPTH'          ,'C2','WAT_DEPTH',' (m)'                                 ,& ! 330
 'WAT_TEMPK'          ,'C2','WAT_TEMP',' (K)'                                  ,& ! 331
 'HEAD1_MSL'          ,'C2','WATER SFC HEAD REL TO MSL',' (m)'                  / ! 332

!SFCG_SWM_ACTIVE

! Special SFC GRID nudging fields

data fldlib(1:4,335:337)/ &
 'SFCWAT_NUD'         ,'C2' ,'SFCWAT_NUD',' (mm/day)'                          ,& ! 335 ! fast can nud
 'SFCTEMP_NUD'        ,'C2' ,'SFCTEMP_NUD',' (K)'                              ,& ! 336 ! fast can nud
 'FRACLIQ_NUD'        ,'C2' ,'FRACLIQ_NUD',' ( )'                               / ! 337 ! fast can nud

! SFC GRID CELLS - DIF2 fields

data fldlib(1:4,341:361)/ &
 'SFCG_VELS_DIF2'      ,'C2' ,'SFCG WIND SPEED DIF2',' (m s:S2:-1  )'           ,& ! 341
 'SFCG_AIRTEMPK_DIF2'  ,'C2' ,'SFCG ATM TEMP DIF2',' (K)'                       ,& ! 342
 'SFCG_AIRRRV_DIF2'    ,'C2' ,'SFCG ATM RRV DIF2',' (g kg:S2:-1  )'             ,& ! 343
 'SFCG_CANTEMPK_DIF2'  ,'C2' ,'SFCG CANOPY AIR TEMP DIF2',' (K)'                ,& ! 344
 'SFCG_CANRRV_DIF2'    ,'C2' ,'SFCG CANOPY VAP MIX RATIO DIF2',' (g kg:S2:-1  )',& ! 345
 'SFCG_SKINTEMPK_DIF2' ,'C2' ,'SFCG VEG/GROUND/SFCWATER/SEA TEMP DIF2',' (K)'   ,& ! 346
 'SFCG_SENSFLUX_DIF2'  ,'C2' ,'SFCG SENS HEAT FLUX DIF2',' (W m:S2:-2  )'       ,& ! 347
 'SFCG_LATFLUX_DIF2'   ,'C2' ,'SFCG LAT HEAT FLUX DIF2',' (W m:S2:-2  )'        ,& ! 348
 'SFCG_VAPFLUX_DIF2'   ,'C2' ,'SFCG VAP FLUX DIF2',' (mm/day)'                  ,& ! 349 ! fast can nud
 'SFCG_PCP_DIF2'       ,'C2' ,'SFCG PCP DIF2',' (mm/day)'                       ,& ! 350 ! fast can nud
 'SFCG_RUNOFF_DIF2'    ,'C2' ,'SFCG RUNOFF DIF2',' (mm/day)'                    ,& ! 351
 'SFCG_SFCTEMP_DIF2'   ,'C2' ,'SFCG SFCTEMP DIF2',' (k)'                        ,& ! 352 ! fast can nud
 'SFCG_FRACLIQ_DIF2'   ,'C2' ,'SFCG FRACLIQ DIF2',' ( )'                        ,& ! 353 ! fast can nud
 'SFCG_AIRTEMPK_DMIN_DIF2' ,'C2' ,'SFCG ATM DMIN TEMP DIF2',' (K)'              ,& ! 354
 'SFCG_AIRTEMPK_DMAX_DIF2' ,'C2' ,'SFCG ATM DMAX TEMP DIF2',' (K)'              ,& ! 355
 'SFCG_AIRTEMPK_DSPAN_DIF2','C2' ,'SFCG ATM DSPAN TEMP DIF2',' (K)'             ,& ! 356
 'SFCG_CANTEMPK_DMIN_DIF2' ,'C2' ,'SFCG CAN DMIN TEMP DIF2',' (K)'              ,& ! 357
 'SFCG_CANTEMPK_DMAX_DIF2' ,'C2' ,'SFCG CAN DMAX TEMP DIF2',' (K)'              ,& ! 358
 'SFCG_CANTEMPK_DSPAN_DIF2','C2' ,'SFCG CAN DSPAN TEMP DIF2',' (K)'             ,& ! 359
 'SFCG_VAPFLUXF_DIF2'      ,'L2' ,'SFCG VAP FLUX FRACTION DIF2',' ()'           ,& ! 360 ! fast can nud
 'HEAD_WTAB_DIF2'          ,'C2' ,'HEAD AT WATER TABLE DIF2',' (m)'              / ! 361

! SFC GRID CELLS - DIF4 fields

data fldlib(1:4,362:372)/ &
 'SFCG_CANTEMPK_DIF4'  ,'C2' ,'SFCG CANOPY AIR TEMP DIF4',' (K)'                ,& ! 362
 'SFCG_CANRRV_DIF4'    ,'C2' ,'SFCG CANOPY VAP MIX RATIO DIF4',' (g kg:S2:-1  )',& ! 363
 'SFCG_SKINTEMPK_DIF4' ,'C2' ,'SFCG VEG/GROUND/SFCWATER/SEA TEMP DIF4',' (K)'   ,& ! 364
 'SFCG_SENSFLUX_DIF4'  ,'C2' ,'SFCG SENS HEAT FLUX DIF4',' (W m:S2:-2  )'       ,& ! 365
 'SFCG_LATFLUX_DIF4'   ,'C2' ,'SFCG LAT HEAT FLUX DIF4',' (W m:S2:-2  )'        ,& ! 366
 'SFCG_VAPFLUX_DIF4'   ,'C2' ,'SFCG VAP FLUX DIF4',' (kg m:S2:-2   s:S2:-1  )'  ,& ! 367
 'SFCG_RUNOFF_DIF4'    ,'C2' ,'SFCG RUNOFF DIF4',' (mm/day)'                    ,& ! 368
 'SFCG_CANTEMPK_DMIN_DIF4' ,'C2' ,'SFCG CAN DMIN TEMP DIF4',' (K)'              ,& ! 369
 'SFCG_CANTEMPK_DMAX_DIF4' ,'C2' ,'SFCG CAN DMAX TEMP DIF4',' (K)'              ,& ! 370
 'SFCG_CANTEMPK_DSPAN_DIF4','C2' ,'SFCG CAN DSPAN TEMP DIF4',' (K)'             ,& ! 371
 'SFCG_VAPFLUXF_DIF4'   ,'L2' ,'SFCG VAP FLUX FRACTION DIF4',' ()'               / ! 372

! SFC GRID CELLS - ATM averages

data fldlib(1:4,401:409)/ &
 'ASFCG_VELS'         ,'T2' ,'ASFCG WIND SPEED',' (m s:S2:-1  )'                ,& ! 401
 'ASFCG_AIRTEMPK'     ,'T2' ,'ASFCG ATM TEMP',' (K)'                            ,& ! 402
 'ASFCG_AIRRRV'       ,'T2' ,'ASFCG ATM VAP MIX RATIO',' (g kg:S2:-1  )'        ,& ! 403
 'ASFCG_CANTEMPK'     ,'T2' ,'ASFCG CANOPY AIR TEMP',' (K)'                     ,& ! 404
 'ASFCG_CANRRV'       ,'T2' ,'ASFCG CANOPY VAP MIX RATIO',' (g kg:S2:-1  )'     ,& ! 405
 'ASFCG_SKINTEMPK'    ,'T2' ,'ASFCG SKIN TEMP',' (K)'                           ,& ! 406
 'ASFCG_SENSFLUX'     ,'T2' ,'ASFCG SENS HEAT FLUX',' (W m:S2:-2  )'            ,& ! 407
 'ASFCG_LATFLUX'      ,'T2' ,'ASFCG LAT HEAT FLUX',' (W m:S2:-2  )'             ,& ! 408
 'ASFCG_VAPFLUX'      ,'T2' ,'ASFCG VAP FLUX',' (kg m:S2:-2   s:S2:-1  )'        / ! 409

! SFC GRID CELLS - ATM averages of DIF2 fields

data fldlib(1:4,421:435)/ &
 'ASFCG_VELS_DIF2'     ,'T2' ,'ASFCG WIND SPEED DIF2',' (m s:S2:-1  )'           ,& ! 421
 'ASFCG_AIRTEMPK_DIF2' ,'T2' ,'ASFCG ATM TEMP DIF2',' (K)'                       ,& ! 422
 'ASFCG_AIRRRV_DIF2'   ,'T2' ,'ASFCG ATM RRV DIF2',' (g kg:S2:-1  )'             ,& ! 423
 'ASFCG_CANTEMPK_DIF2' ,'T2' ,'ASFCG CANOPY AIR TEMP DIF2',' (K)'                ,& ! 424
 'ASFCG_CANRRV_DIF2'   ,'T2' ,'ASFCG CANOPY VAP MIX RATIO DIF2',' (g kg:S2:-1  )',& ! 425
 'ASFCG_SKINTEMPK_DIF2','T2' ,'ASFCG VEG/GROUND/SFCWATER/SEA TEMP DIF2',' (K)'   ,& ! 426
 'ASFCG_SENSFLUX_DIF2' ,'T2' ,'ASFCG SENS HEAT FLUX DIF2',' (W m:S2:-2  )'       ,& ! 427
 'ASFCG_LATFLUX_DIF2'  ,'T2' ,'ASFCG LAT HEAT FLUX DIF2',' (W m:S2:-2  )'        ,& ! 428
 'ASFCG_VAPFLUX_DIF2'  ,'T2' ,'ASFCG VAP FLUX DIF2',' (kg m:S2:-2   s:S2:-1  )'  ,& ! 429
 'ASFCG_AIRTEMPK_DMIN_DIF2' ,'T2' ,'ASFCG ATM DMIN TEMP DIF2',' (K)'             ,& ! 430
 'ASFCG_AIRTEMPK_DMAX_DIF2' ,'T2' ,'ASFCG ATM DMAX TEMP DIF2',' (K)'             ,& ! 431
 'ASFCG_AIRTEMPK_DSPAN_DIF2','T2' ,'ASFCG ATM DSPAN TEMP DIF2',' (K)'            ,& ! 432
 'ASFCG_CANTEMPK_DMIN_DIF2' ,'T2' ,'ASFCG CAN DMIN TEMP DIF2',' (K)'             ,& ! 433
 'ASFCG_CANTEMPK_DMAX_DIF2' ,'T2' ,'ASFCG CAN DMAX TEMP DIF2',' (K)'             ,& ! 434
 'ASFCG_CANTEMPK_DSPAN_DIF2','T2' ,'ASFCG CAN DSPAN TEMP DIF2',' (K)'             / ! 435

! SFC GRID CELLS - ATM averages of DIF4 fields

data fldlib(1:4,451:465)/ &
 'ASFCG_VELS_DIF4'     ,'T2' ,'ASFCG WIND SPEED DIF4',' (m s:S2:-1  )'           ,& ! 451
 'ASFCG_AIRTEMPK_DIF4' ,'T2' ,'ASFCG ATM TEMP DIF4',' (K)'                       ,& ! 452
 'ASFCG_AIRRRV_DIF4'   ,'T2' ,'ASFCG ATM MIX RATIO DIF4',' (g kg:S2:-1  )'       ,& ! 453
 'ASFCG_CANTEMPK_DIF4' ,'T2' ,'ASFCG CANOPY AIR TEMP DIF4',' (K)'                ,& ! 454
 'ASFCG_CANRRV_DIF4'   ,'T2' ,'ASFCG CANOPY VAP MIX RATIO DIF4',' (g kg:S2:-1  )',& ! 455
 'ASFCG_SKINTEMPK_DIF4','T2' ,'ASFCG VEG/GROUND/SFCWATER/SEA TEMP DIF4',' (K)'   ,& ! 456
 'ASFCG_SENSFLUX_DIF4' ,'T2' ,'ASFCG SENS HEAT FLUX DIF4',' (W m:S2:-2  )'       ,& ! 457
 'ASFCG_LATFLUX_DIF4'  ,'T2' ,'ASFCG LAT HEAT FLUX DIF4',' (W m:S2:-2  )'        ,& ! 458
 'ASFCG_VAPFLUX_DIF4'  ,'T2' ,'ASFCG VAP FLUX DIF4',' (kg m:S2:-2   s:S2:-1  )'  ,& ! 459
 'ASFCG_AIRTEMPK_DMIN_DIF4' ,'T2' ,'ASFCG ATM DMIN TEMP DIF4',' (K)'             ,& ! 460
 'ASFCG_AIRTEMPK_DMAX_DIF4' ,'T2' ,'ASFCG ATM DMAX TEMP DIF4',' (K)'             ,& ! 461
 'ASFCG_AIRTEMPK_DSPAN_DIF4','T2' ,'ASFCG ATM DSPAN TEMP DIF4',' (K)'            ,& ! 462
 'ASFCG_CANTEMPK_DMIN_DIF4' ,'T2' ,'ASFCG CAN DMIN TEMP DIF4',' (K)'             ,& ! 463
 'ASFCG_CANTEMPK_DMAX_DIF4' ,'T2' ,'ASFCG CAN DMAX TEMP DIF4',' (K)'             ,& ! 464
 'ASFCG_CANTEMPK_DSPAN_DIF4','T2' ,'ASFCG CAN DSPAN TEMP DIF4',' (K)'             / ! 465

! GRID GEOMETRY - 3D

data fldlib(1:4,501:503)/ &
 'ARV'           ,'V3' ,'AREA OF GRID CELL V-FACE',' (m:S2:2  )'            ,& ! 501
 'ARW'           ,'W3' ,'AREA OF GRID CELL W-FACE',' (m:S2:2  )'            ,& ! 502
 'VOLT'          ,'T3' ,'GRID T-CELL VOLUME',' (m:S2:3  )'                   / ! 503

! GRID GEOMETRY - 2D

data fldlib(1:4,511:537)/ &
 'TOPM'          ,'M2' ,'TOPOGRAPHY HEIGHT',' (m)'                          ,& ! 511
 'TOPW'          ,'W2' ,'TOPOGRAPHY HEIGHT AT W',' (m)'                     ,& ! 512
 'GLATM'         ,'M2' ,'LATITUDE AT M',' (deg)'                            ,& ! 513
 'GLONM'         ,'M2' ,'LONGITUDE AT M',' (deg)'                           ,& ! 514
 'GLATV'         ,'V2' ,'LATITUDE AT V',' (deg)'                            ,& ! 515
 'GLONV'         ,'V2' ,'LONGITUDE AT V',' (deg)'                           ,& ! 516
 'GLATW'         ,'T2' ,'LATITUDE',' (deg)'                                 ,& ! 517
 'GLONW'         ,'T2' ,'LONGITUDE',' (deg)'                                ,& ! 518
 'LPM'           ,'M2' ,'LOWEST PREDICTED M LEVEL',' ( )'                   ,& ! 519
 'LPV'           ,'V2' ,'LOWEST PREDICTED V LEVEL',' ( )'                   ,& ! 520
 'LCV'           ,'V2' ,'LOWEST ACTIVE V CONTROL VOL',' ( )'                ,& ! 521
 'LPW'           ,'W2' ,'LOWEST PREDICTED W LEVEL',' ( )'                   ,& ! 522
 'LSW'           ,'W2' ,'NUMBER OF SFC W LEVELS',' ( )'                     ,& ! 523
 'XEM'           ,'M2' ,'EARTH-X COORD OF M POINT',' ( )'                   ,& ! 524
 'YEM'           ,'M2' ,'EARTH-Y COORD OF M POINT',' ( )'                   ,& ! 525
 'ZEM'           ,'M2' ,'EARTH-Z COORD OF M POINT',' ( )'                   ,& ! 526
 'XEV'           ,'V2' ,'EARTH-X COORD OF V POINT',' ( )'                   ,& ! 527
 'YEV'           ,'V2' ,'EARTH-Y COORD OF V POINT',' ( )'                   ,& ! 528
 'ZEV'           ,'V2' ,'EARTH-Z COORD OF V POINT',' ( )'                   ,& ! 529
 'XEW'           ,'W2' ,'EARTH-X COORD OF W POINT',' ( )'                   ,& ! 530
 'YEW'           ,'W2' ,'EARTH-Y COORD OF W POINT',' ( )'                   ,& ! 531
 'ZEW'           ,'W2' ,'EARTH-Z COORD OF W POINT',' ( )'                   ,& ! 532
 'DNU'           ,'V2' ,'DNU',' (m)'                                        ,& ! 533
 'DNV'           ,'V2' ,'DNV',' (m)'                                        ,& ! 534
 'ARM0'          ,'W2' ,'SFC AREA OF M CELL',' (m:S2:2  )'                  ,& ! 535
 'ARW0'          ,'W2' ,'SFC AREA OF W CELL',' (m:S2:2  )'                  ,& ! 536
 'TOPMW'         ,'H2' ,'TOPOGRAPHY HEIGHT',' (m)'                           / ! 537

! ITAB_M MEMBERS - 2D

data fldlib(1:4,541:550)/ &
 'ITAB_M_NPOLY'  ,'M2' ,'ITAB_M_NPOLY',' ( )'                               ,& ! 541
 'ITAB_M_IMGLOBE','M2' ,'ITAB_M_IMGLOBE',' ( )'                             ,& ! 542
 'ITAB_M_MRLM'   ,'M2' ,'ITAB_M_MRLM',' ( )'                                ,& ! 543
 'ITAB_M_MRLM_OR','M2' ,'ITAB_M_MRLM_ORIG',' ( )'                           ,& ! 544
 'ITAB_M_MROW'   ,'M2' ,'ITAB_M_MROW',' ( )'                                ,& ! 545
 'ITAB_M_NGR'    ,'M2' ,'ITAB_M_NGR',' ( )'                                 ,& ! 546
 'ITAB_M_IV'     ,'M2' ,'ITAB_M_IV',' ( )'                                  ,& ! 547
 'ITAB_M_IW'     ,'M2' ,'ITAB_M_IW',' ( )'                                  ,& ! 548
 'ITAB_M_IM'     ,'M2' ,'ITAB_M_IM',' ( )'                                  ,& ! 549
 'ITAB_M_IRANK'  ,'M2' ,'ITAB_M_IRANK',' ( )'                                / ! 550

! ITAB_V MEMBERS - 2D

data fldlib(1:4,551:557)/  &
 'ITAB_V_IVP'    ,'V2' ,'ITAB_V_IVP',' ( )'                                 ,& ! 551
 'ITAB_V_IRANK'  ,'V2' ,'ITAB_V_IRANK',' ( )'                               ,& ! 552
 'ITAB_V_IVGLOBE','V2' ,'ITAB_V_IVGLOBE',' ( )'                             ,& ! 553
 'ITAB_V_MRLV'   ,'V2' ,'ITAB_V_MRLV',' ( )'                                ,& ! 554
 'ITAB_V_FARW'   ,'V2' ,'ITAB_V_FARW',' ( )'                                ,& ! 555
 'ITAB_V_IM'     ,'V2' ,'ITAB_V_IM',' ( )'                                  ,& ! 556
 'ITAB_V_IW'     ,'V2' ,'ITAB_V_IW',' ( )'                                   / ! 557
!'ITAB_V_IV'     ,'V2' ,'ITAB_V_IV',' ( )'                                  ,& ! 558

! ITAB_W MEMBERS - 2D

data fldlib(1:4,561:584)/ &

 'ITAB_W_NPOLY'     ,'W2' ,'ITAB_W_NPOLY',' ( )'                            ,& ! 561
 'ITAB_W_IWP'       ,'W2' ,'ITAB_W_IWP',' ( )'                              ,& ! 562
 'ITAB_W_IRANK'     ,'W2' ,'ITAB_W_IRANK',' ( )'                            ,& ! 563
 'ITAB_W_IWGLOBE'   ,'W2' ,'ITAB_W_IWGLOBE',' ( )'                          ,& ! 564
 'ITAB_W_MRLW'      ,'W2' ,'ITAB_W_MRLW',' ( )'                             ,& ! 565
 'ITAB_W_MRLW_OR'   ,'W2' ,'ITAB_W_MRLW_ORIG',' ( )'                        ,& ! 566
 'ITAB_W_NGR'       ,'W2' ,'ITAB_W_NGR',' ( )'                              ,& ! 567
 'ITAB_W_IM'        ,'W2' ,'ITAB_W_IM',' ( )'                               ,& ! 568
 'ITAB_W_IV'        ,'W2' ,'ITAB_W_IV',' ( )'                               ,& ! 569
 'ITAB_W_IW'        ,'W2' ,'ITAB_W_IW',' ( )'                               ,& ! 570
 'ITAB_W_DIRV'      ,'W2' ,'ITAB_W_DIRV',' ( )'                             ,& ! 571
 'ITAB_W_FARM'      ,'W2' ,'ITAB_W_FARM',' ( )'                             ,& ! 572
 'ITAB_W_FARV'      ,'W2' ,'ITAB_W_FARV',' ( )'                             ,& ! 573
 'ITAB_W_IWNUD'     ,'W2' ,'ITAB_W_IWNUD',' ( )'                            ,& ! 574
 'ITAB_W_FNUD'      ,'W2' ,'ITAB_W_FNUD',' ( )'                             ,& ! 575
 'ITAB_W_JLAND1'    ,'W2' ,'ITAB_W_JLAND1',' ( )'                           ,& ! 576
 'ITAB_W_JLAND2'    ,'W2' ,'ITAB_W_JLAND2',' ( )'                           ,& ! 577
 'ITAB_W_JLAKE1'    ,'W2' ,'ITAB_W_JLAKE1',' ( )'                           ,& ! 578
 'ITAB_W_JLAKE2'    ,'W2' ,'ITAB_W_JLAKE2',' ( )'                           ,& ! 579
 'ITAB_W_JSEA1'     ,'W2' ,'ITAB_W_JSEA1',' ( )'                            ,& ! 580
 'ITAB_W_JSEA2'     ,'W2' ,'ITAB_W_JSEA2',' ( )'                            ,& ! 581
 'ITAB_W_JSFC2'     ,'W2' ,'ITAB_W_J2SFC',' ( )'                            ,& ! 582
 'ITAB_W_IWSFC'     ,'W2' ,'ITAB_W_IWSFC',' ( )'                            ,& ! 583
 'ITAB_W_JASFC'     ,'W2' ,'ITAB_W_JASFC',' ( )'                             / ! 584

! ITAB_WSFC MEMBERS - 2D

data fldlib(1:4,602:614)/ &

 'ITAB_WSFC_IWGLOBE'   ,'C2' ,'ITAB_WSFC_IWGLOBE',' ( )'                      ,& ! 602
 'ITAB_WSFC_IRANK'     ,'C2' ,'ITAB_WSFC_IRANK',' ( )'                        ,& ! 603
 'ITAB_WSFC_NWATM'     ,'C2' ,'ITAB_WSFC_NWATM',' ( )'                        ,& ! 604
 'ITAB_WSFC_IWATM'     ,'C2' ,'ITAB_WSFC_IWATM',' ( )'                        ,& ! 605
 'ITAB_WSFC_KWATM'     ,'C2' ,'ITAB_WSFC_KWATM',' ( )'                        ,& ! 606
 'ITAB_WSFC_ARC'       ,'C2' ,'ITAB_WSFC_ARC',' ( )'                          ,& ! 607
 'ITAB_WSFC_ARCOARSFC' ,'C2' ,'ITAB_WSFC_ARCOARSFC',' ( )'                    ,& ! 608
 'ITAB_WSFC_ARCOARIW'  ,'C2' ,'ITAB_WSFC_ARCOARIW',' ( )'                     ,& ! 609
 'ITAB_WSFC_ARCOARKW'  ,'C2' ,'ITAB_WSFC_ARCOARKW',' ( )'                     ,& ! 610
 'ITAB_WSFC_NPOLY'     ,'C2' ,'ITAB_WSFC_NPOLY',' ( )'                        ,& ! 611
 'ITAB_WSFC_IMN'       ,'C2' ,'ITAB_WSFC_IMN',' ( )'                          ,& ! 612
 'ITAB_WSFC_IVN'       ,'C2' ,'ITAB_WSFC_IVN',' ( )'                          ,& ! 613
 'ITAB_WSFC_IWN'       ,'C2' ,'ITAB_WSFC_IWN',' ( )'                           / ! 614

! Monthly and Daily averaged fields; daily min & max fields - 2D

data fldlib(1:4,722:753)/ &

 'PRESS_DAVG'       ,'T2','DAY-AVG SFC PRESSURE',' (hPa)'                   ,& ! 722
 'ZONAL_WINDW_DAVG' ,'T2','DAY-AVG SFC ZONAL WIND',' (m s:S2:-1  )'         ,& ! 723
 'MERID_WINDW_DAVG' ,'T2','DAY-AVG SFC MERID WIND',' (m s:S2:-1  )'         ,& ! 724
 'RSHORT_DAVG'      ,'T2','DAY-AVG SFC DOWNWARD S/W FLX',' (W m:S2:-2  )'   ,& ! 725
 'TEMPK_DAVG'       ,'T2','DAY-AVG SFC TEMPERATURE',' (K)'                  ,& ! 726
 'TEMPK_DMIN'       ,'T2','DAY-MIN SFC TEMPERATURE',' (K)'                  ,& ! 727
 'TEMPK_DMAX'       ,'T2','DAY-MAX SFC TEMPERATURE',' (K)'                  ,& ! 728
 'ACCPMIC_DTOT'     ,'T2','DAY-ACCUM MICPHYS PRECIP',' (kg m:S2:-2  )'      ,& ! 729
 'ACCPCON_DTOT'     ,'T2','DAY-ACCUM CUPARM PRECIP',' (kg m:S2:-2  )'       ,& ! 730
 'ACCPBOTH_DTOT'    ,'T2','DAY-ACCUM MICPHYS + CONV PCP',' (kg m:S2:-2  )'  ,& ! 731
 'PRESS_UL_DAVG'      ,'T2','DAY-AVG UL PRESSURE',' (hPa)'                  ,& ! 732
 'ZONAL_WINDW_UL_DAVG','T2','DAY-AVG UL ZONAL WIND',' (m s:S2:-1  )'        ,& ! 733
 'MERID_WINDW_UL_DAVG','T2','DAY-AVG UL MERID WIND',' (m s:S2:-1  )'        ,& ! 734
 'CANTEMPK_DAVG'    ,'C2','DAY-AVG CANOPY TEMPERATURE',' (K)'               ,& ! 735
 'CANTEMPK_DMIN'    ,'C2','DAY-MIN CANOPY TEMPERATURE',' (K)'               ,& ! 736
 'CANTEMPK_DMAX'    ,'C2','DAY-MAX CANOPY TEMPERATURE',' (K)'               ,& ! 737
 'VEGTEMPK_DAVG'    ,'L2','DAY-AVG VEG TEMPERATURE',' (K)'                  ,& ! 738
 'VEGTEMPK_DMIN'    ,'L2','DAY-MIN VEG TEMPERATURE',' (K)'                  ,& ! 739
 'VEGTEMPK_DMAX'    ,'L2','DAY-MAX VEG TEMPERATURE',' (K)'                  ,& ! 740
 'SOILTEMPK_DAVG'   ,'L2','DAY-AVG SOIL TEMPERATURE',' (K)'                 ,& ! 741
 'SOILTEMPK_DMIN'   ,'L2','DAY-MIN SOIL TEMPERATURE',' (K)'                 ,& ! 742
 'SOILTEMPK_DMAX'   ,'L2','DAY-MAX SOIL TEMPERATURE',' (K)'                 ,& ! 743
 'SFCG_AIRTEMPK_DAVG' ,'C2','SFCG DAY-AVG ATM TEMP',' (K)'                  ,& ! 744
 'SFCG_AIRTEMPK_DMIN' ,'C2','SFCG DAY-MIN ATM TEMP',' (K)'                  ,& ! 745
 'SFCG_AIRTEMPK_DMAX' ,'C2','SFCG DAY-MAX ATM TEMP',' (K)'                  ,& ! 746
 'SFCG_CANTEMPK_DAVG' ,'C2','SFCG DAY-AVG CAN TEMP',' (K)'                  ,& ! 747
 'SFCG_CANTEMPK_DMIN' ,'C2','SFCG DAY-MIN CAN TEMP',' (K)'                  ,& ! 748
 'SFCG_CANTEMPK_DMAX' ,'C2','SFCG DAY-MAX CAN TEMP',' (K)'                  ,& ! 749
 'SFCG_SENSFLUX_DAVG' ,'C2','SFCG DAY-AVG SENS HEAT FLUX',' (W m:S2:-2  )'  ,& ! 750
 'SFCG_LATFLUX_DAVG'  ,'C2','SFCG DAY-AVG LAT HEAT FLUX',' (W m:S2:-2  )'   ,& ! 751
 'SENSFLUX_DAVG'    ,'T2','DAY-AVG SENS HEAT FLUX',' (W m:S2:-2  )'         ,& ! 752
 'LATFLUX_DAVG'     ,'T2','DAY-AVG LAT HEAT FLUX',' (W m:S2:-2  )'           / ! 753

! Miscellaneous and new additions

data fldlib(1:4,771:796)/ &
 'RHO_OBS'       ,'T3' ,'NUDGING OBS AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 771
 'THETA_OBS'     ,'T3' ,'NUDGING OBS THETA',' (K)'                          ,& ! 772
 'RRW_OBS'       ,'T3' ,'NUDGING OBS VAPOR MIX RATIO',' (g kg:S2:-1  )'     ,& ! 773
 'UZONAL_OBS'    ,'T3' ,'NUDGING OBS ZONAL WIND',' (m s:S2:-1  )'           ,& ! 774
 'UMERID_OBS'    ,'T3' ,'NUDGING OBS MERID WIND',' (m s:S2:-1  )'           ,& ! 775
 'RHO_SIM'       ,'T3' ,'NUDGING SIM AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 776
 'THETA_SIM'     ,'T3' ,'NUDGING SIM THETA',' (K)'                          ,& ! 777
 'RRW_SIM'       ,'T3' ,'NUDGING SIM VAPOR MIX RATIO',' (g kg:S2:-1  )'     ,& ! 778
 'UZONAL_SIM'    ,'T3' ,'NUDGING SIM ZONAL WIND',' (m s:S2:-1  )'           ,& ! 779
 'UMERID_SIM'    ,'T3' ,'NUDGING SIM MERID WIND',' (m s:S2:-1  )'           ,& ! 780
 'RHO_OBS_SIM'   ,'T3' ,'NUDGING DIF AIR DENSITY',' (kg m:S2:-3  )'         ,& ! 781
 'THETA_OBS_SIM' ,'T3' ,'NUDGING DIF THETA',' (K)'                          ,& ! 782
 'RRW_OBS_SIM'   ,'T3' ,'NUDGING DIF VAPOR MIX RATIO',' (g kg:S2:-1  )'     ,& ! 783
 'UZONAL_OBS_SIM','T3' ,'NUDGING DIF ZONAL WIND',' (m s:S2:-1  )'           ,& ! 784
 'UMERID_OBS_SIM','T3' ,'NUDGING DIF MERID WIND',' (m s:S2:-1  )'           ,& ! 785
 'VXE'           ,'T3' ,'EARTH CARTESIAN X WIND',' (m s:S2:-1  )'           ,& ! 786
 'VYE'           ,'T3' ,'EARTH CARTESIAN Y WIND',' (m s:S2:-1  )'           ,& ! 787
 'VZE'           ,'T3' ,'EARTH CARTESIAN Z WIND',' (m s:S2:-1  )'           ,& ! 788
 'PBLH'          ,'T2' ,'PBL HEIGHT',' (m)'                                 ,& ! 789
 'VKH'           ,'T3' ,'EDDY DIFFUSIVITY',' (m:S2:2 s:S2:-1  )'            ,& ! 790
 'RRW_HCONV'     ,'T2' ,'TOTAL WATER HORIZ CONV',' (kg m:S2:-2   s:S2:-1  )',& ! 791
 'RRV_HCONV'     ,'T2' ,'WATER VAPOR HORIZ CONV',' (kg m:S2:-2   s:S2:-1  )',& ! 792
 'CLDNUM'        ,'T2' ,'CLOUD # CONCEN (GEOG)',' (# mg:S2:-1  )'           ,& ! 793
 'ADDSC1_ZINT'   ,'T2' ,'VERTICAL INTEGRAL OF ADDED SCALAR 1',' ( )'        ,& ! 794
 'ADDSC2_ZINT'   ,'T2' ,'VERTICAL INTEGRAL OF ADDED SCALAR 2',' ( )'        ,& ! 795
 'ADDSC1P2_ZINT' ,'T2' ,'VERTICAL INTEGRAL OF ADDED SCALARS 1+2',' ( )'      / ! 796

! External fields

data fldlib(1:4,801:803)/ &
 'VORTP'         ,'P3' ,'VORTP',' (s:S2:-1  )'                              ,& ! 801
 'VORTN'         ,'N3' ,'VORTN',' (s:S2:-1  )'                              ,& ! 802
 'RKE'           ,'T3' ,'RKE',' (s:S2:-1  )'                                 / ! 803

if (len_trim(fldname0) >= 5) then
   if ( fldname0(1:5) == 'CHEM_' .or. &
        fldname0(1:5) == 'chem_' ) then
      call oplot_chem_lib(kk,ii,infotyp,fldname0,wtbot,wttop,fldval,notavail)
      return
   endif
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
   elseif (op%stagpt == 'L' .or. op%stagpt == 'R' .or. op%stagpt == 'S' .or. &
           op%stagpt == 'C' .or. op%stagpt == 'E') then
      i = 2
   endif
endif

if (infotyp == 'VALUV') then
   op%stagpt = 'V'
   op%dimens = '3'
   icase = 4
endif

if (op%stagpt == 'L' .and. mland < 2) then
   notavail = 3
   return
endif

if (op%stagpt == 'R' .and. mlake < 2) then
   notavail = 3
   return
endif

if (op%stagpt == 'S' .and. msea < 2) then
   notavail = 3
   return
endif

notavail = 0
kp = min(k+1,mza)

! Try these to prevent out-of-bounds access when infotyp == 'UNITS'

iland = max(2,i-omland)
ilake = max(2,i-omlake)
isea  = max(2,i-omsea)

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

   if (.not. allocated(wmc)) go to 1000

   fldval = wtbot * wmc(k ,i) &
          + wttop * wmc(kp,i)

case(3) ! 'VMP'  ! obsolete -- placeholder for now in data statement above

   if (.not. allocated(vmc)) go to 1000

   fldval = wtbot * vmc(k ,i) &
          + wttop * vmc(kp,i)

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

   if (.not. allocated(wc)) go to 1000

   fldval = wtbot * wc(k ,i) &
          + wttop * wc(kp,i)

case(6) ! 'RHO'

   if (.not. allocated(rho)) go to 1000

   fldval = wtbot * rho(k ,i) &
          + wttop * rho(kp,i)

case(7) ! 'PRESS'

   if (.not. allocated(press)) go to 1000

   fldval = press(k,i) * .01

   if (nl%test_case == 2 .or. nl%test_case == 5) then
      fldval = press(k,i)
   endif

case(8) ! 'THIL'

   if (.not. allocated(thil)) go to 1000

   fldval = wtbot * thil(k ,i) &
          + wttop * thil(kp,i)

case(9) ! 'THETA'

   if (.not. allocated(theta)) go to 1000

   fldval = wtbot * theta(k ,i) &
          + wttop * theta(kp,i)

case(10) ! 'AIRTEMPK'

   if (.not. allocated(tair)) go to 1000

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i)

case(11) ! 'AIRTEMPC'

   if (.not. allocated(tair)) go to 1000

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i) - 273.15

case(12) ! 'RR_W'

   if (.not. allocated(rr_w)) go to 1000

   fldval = (wtbot * rr_w(k ,i) &
          +  wttop * rr_w(kp,i)) * 1.e3

case(13) ! 'RR_V'

   if (.not. allocated(rr_v)) go to 1000

   fldval = (wtbot * rr_v(k ,i) &
          +  wttop * rr_v(kp,i)) * 1.e3

case(14) ! 'RR_C'

   if (.not. allocated(rr_c)) go to 1000

   fldval = (wtbot * rr_c(k ,i) &
          +  wttop * rr_c(kp,i)) * 1.e3

case(15) ! 'RR_D'

   if (.not. allocated(rr_d)) go to 1000

   fldval = (wtbot * rr_d(k ,i) &
          +  wttop * rr_d(kp,i)) * 1.e3

case(16) ! 'RR_R'

   if (.not. allocated(rr_r)) go to 1000

   fldval = (wtbot * rr_r(k ,i) &
          +  wttop * rr_r(kp,i)) * 1.e3

case(17) ! 'RR_P'

   if (.not. allocated(rr_p)) go to 1000

   fldval = (wtbot * rr_p(k ,i) &
          +  wttop * rr_p(kp,i)) * 1.e3

case(18) ! 'RR_S'

   if (.not. allocated(rr_s)) go to 1000

   fldval = (wtbot * rr_s(k ,i) &
          +  wttop * rr_s(kp,i)) * 1.e3

case(19) ! 'RR_A'

   if (.not. allocated(rr_a)) go to 1000

   fldval = (wtbot * rr_a(k ,i) &
          +  wttop * rr_a(kp,i)) * 1.e3

case(20) ! 'RR_G'

   if (.not. allocated(rr_g)) go to 1000

   fldval = (wtbot * rr_g(k ,i) &
          +  wttop * rr_g(kp,i)) * 1.e3

case(21) ! 'RR_H'

   if (.not. allocated(rr_h)) go to 1000

   fldval = (wtbot * rr_h(k ,i) &
          +  wttop * rr_h(kp,i)) * 1.e3

case(22) ! 'RR_CP'

   if (.not. allocated(rr_c)) go to 1000
   if (.not. allocated(rr_p)) go to 1000

   fldval = (wtbot * (rr_c(k ,i) + rr_p(k ,i)) &
          +  wttop * (rr_c(kp,i) + rr_p(kp,i))) * 1.e3

case(23) ! 'RR_TOTLIQ'

   fldval = 0.

   if (allocated(rr_c)) fldval = fldval + (wtbot * rr_c(k ,i) &
                                        +  wttop * rr_c(kp,i)) * 1.e3

   if (allocated(rr_d)) fldval = fldval + (wtbot * rr_d(k ,i) &
                                        +  wttop * rr_d(kp,i)) * 1.e3

   if (allocated(rr_r)) fldval = fldval + (wtbot * rr_r(k ,i) &
                                        +  wttop * rr_r(kp,i)) * 1.e3

   if (allocated(rr_g)) then
      if (rr_g(k,i) > rxmin(6)) then
         call qtk(q6(k,i)/rr_g(k,i),tempk,fracliq)
         fldval = fldval + wtbot * rr_g(k,i) * fracliq * 1.e3
      endif

      if (rr_g(kp,i) > rxmin(6)) then
         call qtk(q6(kp,i)/rr_g(kp,i),tempk,fracliq)
         fldval = fldval + wttop * rr_g(kp,i) * fracliq * 1.e3
      endif
   endif

   if (allocated(rr_h)) then
      if (rr_h(k,i) > rxmin(7)) then
         call qtk(q7(k,i)/rr_h(k,i),tempk,fracliq)
         fldval = fldval + wtbot * rr_h(k,i) * fracliq * 1.e3
      endif

      if (rr_h(kp,i) > rxmin(7)) then
         call qtk(q7(kp,i)/rr_h(kp,i),tempk,fracliq)
         fldval = fldval + wttop * rr_h(kp,i) * fracliq * 1.e3
      endif
   endif

case(24) ! 'RR_TOTICE'

   fldval = 0.

   if (allocated(rr_p)) fldval = fldval + (wtbot * rr_p(k ,i) &
                                        +  wttop * rr_p(kp,i)) * 1.e3

   if (allocated(rr_s)) fldval = fldval + (wtbot * rr_s(k ,i) &
                                        +  wttop * rr_s(kp,i)) * 1.e3

   if (allocated(rr_a)) fldval = fldval + (wtbot * rr_a(k ,i) &
                                        +  wttop * rr_a(kp,i)) * 1.e3

   if (allocated(rr_g)) then
      if (rr_g(k,i) > rxmin(6)) then
         call qtk(q6(k,i)/rr_g(k,i),tempk,fracliq)
         fldval = fldval + wtbot * rr_g(k,i) * (1.0 -fracliq) * 1.e3
      endif

      if (rr_g(kp,i) > rxmin(6)) then
         call qtk(q6(kp,i)/rr_g(kp,i),tempk,fracliq)
         fldval = fldval + wttop * rr_g(kp,i) * (1.0 -fracliq) * 1.e3
      endif
   endif

   if (allocated(rr_h)) then
      if (rr_h(k,i) > rxmin(7)) then
         call qtk(q7(k,i)/rr_h(k,i),tempk,fracliq)
         fldval = fldval + wtbot * rr_h(k,i) * (1.0 -fracliq) * 1.e3
      endif

      if (rr_h(kp,i) > rxmin(7)) then
         call qtk(q7(kp,i)/rr_h(kp,i),tempk,fracliq)
         fldval = fldval + wttop * rr_h(kp,i) * (1.0 -fracliq) * 1.e3
      endif
   endif

case(25) ! 'RR_TOTCOND'

!   fldval = (wtbot * (rr_w(k ,i) - rr_v(k ,i)) &
!          +  wttop * (rr_w(kp,i) - rr_v(kp,i))) * 1.e3

   fldval = 0.0

   if (allocated(rr_c)) fldval = fldval + (wtbot * rr_c(k ,i) &
                                        +  wttop * rr_c(kp,i)) * 1.e3

   if (allocated(rr_r)) fldval = fldval + (wtbot * rr_r(k ,i) &
                                        +  wttop * rr_r(kp,i)) * 1.e3

   if (allocated(rr_d)) fldval = fldval + (wtbot * rr_d(k ,i) &
                                        +  wttop * rr_d(kp,i)) * 1.e3

   if (allocated(rr_p)) fldval = fldval + (wtbot * rr_p(k ,i) &
                                        +  wttop * rr_p(kp,i)) * 1.e3

   if (allocated(rr_s)) fldval = fldval + (wtbot * rr_s(k ,i) &
                                        +  wttop * rr_s(kp,i)) * 1.e3

   if (allocated(rr_a)) fldval = fldval + (wtbot * rr_a(k ,i) &
                                        +  wttop * rr_a(kp,i)) * 1.e3

   if (allocated(rr_g)) fldval = fldval + (wtbot * rr_g(k ,i) &
                                        +  wttop * rr_g(kp,i)) * 1.e3

   if (allocated(rr_h)) fldval = fldval + (wtbot * rr_h(k ,i) &
                                        +  wttop * rr_h(kp,i)) * 1.e3

case(26) ! 'CON_C'

   if (.not. allocated(con_c)) go to 1000

   fldval = (wtbot * con_c(k ,i) &
          +  wttop * con_c(kp,i)) * 1.e-6

case(27) ! 'CON_D'

   if (.not. allocated(con_d)) go to 1000

   fldval = (wtbot * con_d(k ,i) &
          +  wttop * con_d(kp,i)) * 1.e-3

case(28) ! 'CON_R'

   if (.not. allocated(con_r)) go to 1000

   fldval = wtbot * con_r(k ,i) &
          + wttop * con_r(kp,i)

case(29) ! 'CON_P'

   if (.not. allocated(con_p)) go to 1000

   fldval = (wtbot * con_p(k ,i) &
          +  wttop * con_p(kp,i)) * 1.e-6

case(30) ! 'CON_S'

   if (.not. allocated(con_s)) go to 1000

   fldval = wtbot * con_s(k ,i) &
          + wttop * con_s(kp,i)

case(31) ! 'CON_A'

   if (.not. allocated(con_a)) go to 1000

   fldval = wtbot * con_a(k ,i) &
          + wttop * con_a(kp,i)

case(32) ! 'CON_G'

   if (.not. allocated(con_g)) go to 1000

   fldval = wtbot * con_g(k ,i) &
          + wttop * con_g(kp,i)

case(33) ! 'CON_H'

   if (.not. allocated(con_h)) go to 1000

   fldval = wtbot * con_h(k ,i) &
          + wttop * con_h(kp,i)

case(34) ! 'CON_CCN'

   if (indp > nccntyp) go to 1000

   if (.not. allocated(ccntyp(indp)%con_ccn)) go to 1000

   if (indp > 0) then

      ! Standard plot of single ccntyp species

      fldval = (wtbot * ccntyp(indp)%con_ccn(k ,i) &
             +  wttop * ccntyp(indp)%con_ccn(kp,i)) * 1.e-6

   else

      ! Special summation over all ccntype species

      fldval = 0.

      do j = 1,nccntyp
         fldval = fldval + (wtbot * ccntyp(j)%con_ccn(k ,i) &
                         +  wttop * ccntyp(j)%con_ccn(kp,i)) * 1.e-6
      enddo

   endif

case(35) ! 'CON_GCCN'

   if (.not. allocated(con_gccn)) go to 1000

   fldval = (wtbot * con_gccn(k ,i) &
          +  wttop * con_gccn(kp,i)) * 1.e-3

case(36) ! 'CON_IFN'

   if (.not. allocated(con_ifn)) go to 1000

   fldval = (wtbot * con_ifn(k ,i) &
          +  wttop * con_ifn(kp,i)) * 1.e-3

case(37) ! 'VKM'

   if (.not. allocated(vkm)) go to 1000

   fldval = wtbot * vkm(k ,i) &
          + wttop * vkm(kp,i)

case(38) ! 'FTHRD'

   if (.not. allocated(fthrd_sw)) goto 1000
   if (.not. allocated(fthrd_lw)) goto 1000

   fldval = wtbot * (fthrd_sw(k ,i) + fthrd_lw(k ,i)) &
          + wttop * (fthrd_sw(kp,i) + fthrd_lw(kp,i))

case(39:42) ! 'SPEEDW','AZIMW','ZONAL_WINDW','MERID_WINDW'

   if (.not. allocated(vxe)) go to 1000

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

case(43:45) ! 'RVORTZM','TVORTZM','RVORTZM_P'

   if (.not. allocated(vc)) go to 1000

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

case(46) ! 'DIVERG'

   if (.not. allocated(vmc)) go to 1000

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

case(47) ! 'VMASSFLUX'

   if (.not. allocated(vmc)) go to 1000

   fldval = vmc(k,i) * arv(k,i)

case(48) ! 'VC_P'

   if (.not. allocated(vc)) go to 1000

   fldval = vc(k,i) - vc_init(k,i)

case(49) ! 'PRESS_P'

   if (.not. allocated(press)) go to 1000

   fldval = press(k,i) - press_init(k,i)

case(50) ! 'RHO_P'

   if (.not. allocated(rho)) go to 1000

   fldval = wtbot * (rho(k ,i) - rho_init(k ,i)) &
          + wttop * (rho(kp,i) - rho_init(kp,i))

! For shallow water test case 5, define RHO_P from reference fields

   if (nl%test_case == 5) then
      call npr_bicubics(zanal00_swtc5,glatw(i),glonw(i),zanal0_swtc5)
      call npr_bicubics(zanal15_swtc5,glatw(i),glonw(i),zanal_swtc5)

      fldval = (rho(k,i) - rho_init(k,i)) - (zanal_swtc5 - zanal0_swtc5)
   endif

case(51) ! 'THETA_P'

   if (.not. allocated(theta)) go to 1000

   fldval = wtbot * (theta(k ,i) - theta_init(k ,i)) &
          + wttop * (theta(kp,i) - theta_init(kp,i))

case(52) ! 'AIRTEMPK_P'

   if (.not. allocated(tair)) go to 1000

   fldval = wtbot * tair(k ,i) &
          + wttop * tair(kp,i) &
          - wtbot * theta_init(k ,i) * (press_init(k ,i) * p00i) ** rocp &
          - wttop * theta_init(kp,i) * (press_init(kp,i) * p00i) ** rocp

case(53) ! 'VMT'

   if (.not. allocated(vmxet)) go to 1000

   iw1 = itab_v(i)%iw(1)
   iw2 = itab_v(i)%iw(2)

   fldval = .5 * (vnx(i) * (vmxet(k,iw1) + vmxet(k,iw2)) &
               +  vny(i) * (vmyet(k,iw1) + vmyet(k,iw2)) &
               +  vnz(i) * (vmzet(k,iw1) + vmzet(k,iw2)))


case(54) ! 'WMT'

   if (.not. allocated(vmxet)) go to 1000

   fldval = ( wnxo2(i) * (vmxet(k,i) + vmxet(k+1,i)) &
            + wnyo2(i) * (vmyet(k,i) + vmyet(k+1,i)) &
            + wnzo2(i) * (vmzet(k,i) + vmzet(k+1,i)) )

case(55) ! 'ADDSC'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000

   fldval = wtbot * addsc(indp)%sclp(k ,i) &
          + wttop * addsc(indp)%sclp(kp,i)

case(56) ! 'ADDSC_P'

   if (indp > naddsc) go to 1000
   if (.not. allocated(addsc(indp)%sclp)) go to 1000
   if (.not. allocated(addsc1_init)) go to 1000

   fldval = addsc(indp)%sclp(k,i) - addsc1_init(k,i)

case(57) ! 'ZPLEV'

   fldval = wtbot * zt(k ) &
          + wttop * zt(kp)

case(58) ! 'QWCON'

   if (.not. allocated(qwcon)) go to 1000

   fldval = (wtbot * qwcon(k ,i) &
          +  wttop * qwcon(kp,i)) * 1.e3

case(59:60) ! 'CO2CON', 'CO2PERT'

   if (.not. allocated(rr_co2)) go to 1000

   fldval = (wtbot * rr_co2(k ,i) &
          +  wttop * rr_co2(kp,i) ) * co2_sh2ppm

   if (icase == 60) fldval = fldval - nl%co2_ppmv_init

case(61) ! 'RH_LIQ'

   if (.not. allocated(rr_v)) go to 1000

   fldval = (wtbot * rr_v(k ,i) * real(rho(k ,i)) * rhovsl_inv(tair(k ,i)-273.15) &
          +  wttop * rr_v(kp,i) * real(rho(kp,i)) * rhovsl_inv(tair(kp,i)-273.15)) * 1.e2

!-----------------------------------------
! ATMOSPHERE - 2D
!-----------------------------------------

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

   if (.not. allocated(vkm_sfc)) go to 1000

   fldval = vkm_sfc(i)

case(71) ! 'USTAR'

   if (.not. allocated(ustar)) go to 1000

   fldval = ustar(i)

case(72) ! 'SENSFLUX'

   if (.not. allocated(sfluxt)) go to 1000

   fldval = sfluxt(i) * cp

case(73) ! 'VAPFLUX'

   if (.not. allocated(sfluxr)) go to 1000

   fldval = sfluxr(i)

case(74) ! 'LATFLUX'

   if (.not. allocated(sfluxr)) go to 1000

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

   if (.not. allocated(wstar)) go to 1000

   fldval = wstar(i)

case(96) ! 'PSFC'

   if (.not. allocated(press)) go to 1000

   fldval = (press(k,i) + (zt(k) - topw(i)) * rho(k,i) * grav) * .01  ! hydrostatic eqn.

case(97) ! 'PMSL'

   if (.not. allocated(press)) go to 1000

   fldval = press(k,i) * (1. - .0065 * zt(k) / (tair(k,i) + .0065 * zt(k)))**(-5.257) * .01  ! hydrostatic eqn.

case(98) ! 'CBMF'

   if (.not. allocated(cbmf)) go to 1000

   fldval = cbmf(i)

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
         raxisi = 1.0 / raxis

         u = (vy * xew(i) - vx * yew(i)) * raxisi
         v = vz * raxis * eradi &
           - (vx * xew(i) + vy * yew(i)) * zew(i) * raxisi * eradi
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

case(106) ! 'RR_V_DIF2'

   if (.not. allocated(rr_v_accum)) go to 1000

   fldval = rr_v_accum_prev0(k,i) - rr_v_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 1.e3 / (time8_prev0 - time8_prev1)
   endif

case(107) ! 'LATHEAT_LIQ_DIF2'

   if (.not. allocated(latheat_liq_accum)) go to 1000

   fldval = latheat_liq_accum_prev0(k,i) - latheat_liq_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 3600. / (time8_prev0 - time8_prev1)
   endif

case(108) ! 'LATHEAT_ICE_DIF2'

   if (.not. allocated(latheat_ice_accum)) go to 1000

   fldval = latheat_ice_accum_prev0(k,i) - latheat_ice_accum_prev1(k,i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 3600. / (time8_prev0 - time8_prev1)
   endif

case(109) ! 'PCPMIC_DIF2'

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   if (.not. allocated(accpmic_prev1)) go to 1000

   fldval = accpmic_prev0(i) - accpmic_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(110) ! 'PCPCON_DIF2'

   if (.not. allocated(aconpr)) go to 1000

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   fldval = accpcon_prev0(i) - accpcon_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(111) ! 'PCPBOTH_DIF2'

! Compute difference involving 1 previously-stored field
! Convert to mm/day if time interval >= 1 second

   if (.not. allocated(accpmic_prev1)) go to 1000

   accpboth_prev0 = accpmic_prev0(i) + accpcon_prev0(i)
   accpboth_prev1 = accpmic_prev1(i) + accpcon_prev1(i)

   fldval = accpboth_prev0 - accpboth_prev1

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(112) ! 'RSHORT_DIF2'

   if (.not. allocated(rshort_accum)) go to 1000

   fldval = rshort_accum_prev0(i) - rshort_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(113) ! 'RSHORTUP_DIF2'

   if (.not. allocated(rshortup_accum)) go to 1000

   fldval = rshortup_accum_prev0(i) - rshortup_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(114) ! 'RLONG_DIF2'

   if (.not. allocated(rlong_accum)) go to 1000

   fldval = rlong_accum_prev0(i) - rlong_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(115) ! 'RLONGUP_DIF2'

   if (.not. allocated(rlongup_accum)) go to 1000

   fldval = rlongup_accum_prev0(i) - rlongup_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(116) ! 'RSHORT_TOP_DIF2'

   if (.not. allocated(rshort_top_accum)) go to 1000

   fldval = rshort_top_accum_prev0(i) - rshort_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(117) ! 'RSHORTUP_TOP_DIF2'

   if (.not. allocated(rshortup_top_accum)) go to 1000

   fldval = rshortup_top_accum_prev0(i) - rshortup_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(118) ! 'RLONGUP_TOP_DIF2'

   if (.not. allocated(rlongup_top_accum)) go to 1000

   fldval = rlongup_top_accum_prev0(i) - rlongup_top_accum_prev1(i)

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(129) ! 'PCPMIC_DIF4'

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 second

   if (.not. allocated(accpmic_prev1)) go to 1000

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

case(130) ! 'PCPCON_DIF4'

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

case(131) ! 'PCPBOTH_DIF4'

! Compute differences involving 3 previously-stored fields
! Convert to mm/day if time interval >= 1 second

   if (.not. allocated(accpmic_prev1)) go to 1000

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

case(132) ! 'PCPMIC_REL4'

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

case(133) ! 'PCPCON_REL4'

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

case(134) ! 'PCPBOTH_REL4'

! Compute relative differences involving 3 previously-stored fields

   if (.not. allocated(accpmic_prev1)) go to 1000

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

case(181) ! 'SAND'

   fldval = real(land%sand(k,iland))

case(182) ! 'CLAY'

   fldval = real(land%clay(k,iland))

case(183) ! 'SILT'

   fldval = real(land%silt(k,iland))

case(184) ! 'ORGAN'

   fldval = real(land%organ(k,iland)) * 1000. ! converts from kg/kg to g/kg

case(185) ! 'BULKDENS_DRYSOIL'

   fldval = real(land%bulkdens_drysoil(k,iland))

case(186) ! 'PH_SOIL'

   fldval = real(land%pH_soil(k,iland))

case(187) ! 'CEC_SOIL'

   fldval = real(land%cec_soil(k,iland))

case(188) ! 'WRESID_VG'

   fldval = real(land%wresid_vg(k,iland))

case(189) ! 'WSAT_VG'

   fldval = real(land%wsat_vg(k,iland))

case(190) ! 'KSAT_VG'

   fldval = real(land%ksat_vg(k,iland))

case(191) ! 'ALPHA_VG'

   fldval = real(land%alpha_vg(k,iland))

case(192) ! 'EN_VG'

   fldval = real(land%en_vg(k,iland))

case(193) ! 'LAMBDA_VG'

   fldval = real(land%lambda_vg(k,iland))

case(194) ! 'SPECIFHEAT_DRYSOIL'

   fldval = real(land%specifheat_drysoil(k,iland))

case(195) ! 'SOIL_ENERGY'

   fldval = land%soil_energy(k,iland) * 1.e-6

case(196) ! 'SOIL_TEMPK'

   call qwtk(land%soil_energy(k,iland),        &
             land%soil_water(k,iland)*1.e3,    &
             land%specifheat_drysoil(k,iland), &
             tempk, fracliq)
   fldval = tempk

case(197) ! 'SOIL_FRACLIQ'

   call qwtk(land%soil_energy(k,iland),        &
             land%soil_water(k,iland)*1.e3,    &
             land%specifheat_drysoil(k,iland), &
             tempk, fracliq)
   fldval = fracliq

case(198) ! 'SOIL_WATER'

   fldval = land%soil_water(k,iland) / land%wsat_vg(k,iland)

case(199) ! 'SFWAT_MASS'

   fldval = land%sfcwater_mass(k,iland)

case(200) ! 'SFWAT_ENERGY'

   if (land%nlev_sfcwater(iland) == 0) then
      notavail = 4
   else
      fldval = land%sfcwater_energy(k,iland) * 1.e-3
   endif

case(201) ! 'SFWAT_TEMPK'

   if (land%nlev_sfcwater(iland) == 0) then
      notavail = 4
   else
      call qtk(land%sfcwater_energy(k,iland),tempk,fracliq)
      fldval = tempk
   endif

case(202) ! 'SFWAT_FRACLIQ'

   if (land%nlev_sfcwater(iland) == 0) then
      notavail = 4
   else
      call qtk(land%sfcwater_energy(k,iland),tempk,fracliq)
      fldval = fracliq
   endif

case(203) ! 'SFWAT_DEPTH'

   fldval = 0.
   do klev = 1,land%nlev_sfcwater(iland)
      fldval = fldval + land%sfcwater_depth(klev,iland)
   enddo

!-----------------------------------------
! LAND CELLS - 2D
!-----------------------------------------

! LAND_CELLS - 2D

case(205) ! 'USDA_TEXT'

   fldval = real(land%usdatext(iland))

case(206) ! 'Z_BEDROCK'

   fldval = real(land%z_bedrock(iland))

case(207) ! 'GPP'

   fldval = real(land%gpp(iland))

case(208) ! 'GLHYMPS_KSAT'

   fldval = real(land%glhymps_ksat(iland))

case(209) ! 'GLHYMPS_KSAT_PFR'

   fldval = real(land%glhymps_ksat_pfr(iland))

case(210) ! 'GLHYMPS_POROS'

   fldval = real(land%glhymps_poros(iland))

case(211) ! 'NLEV_SFWAT'

   fldval = real(land%nlev_sfcwater(iland))

case(212) ! 'VEG_NDVIC'

   fldval = land%veg_ndvic(iland)

case(213) ! 'VEG_TEMPC'

   fldval = land%veg_temp(iland) - 273.15

case(214) ! 'VEG_TEMPK'

   fldval = land%veg_temp(iland)

case(215) ! 'VEG_WATER'

   fldval = land%veg_water(iland)

case(216) ! 'STOM_RESIST'

   fldval = land%stom_resist(iland)

case(217) ! 'SFCWATER_TOT'

   fldval = sum(land%sfcwater_mass(:,iland))

case(218) ! 'SFCWATER_TOP_TEMP'


   if (land%nlev_sfcwater(iland) == 0) then
      notavail = 4
   else
      nls = land%nlev_sfcwater(iland)
      call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)
      fldval = tempk
   endif

case(219) ! 'SOIL_TOP_TEMP'

   call qwtk(land%soil_energy(nzg,iland),        &
             land%soil_water(nzg,iland)*1.e3,    &
             land%specifheat_drysoil(nzg,iland), &
             tempk, fracliq)
   fldval = tempk

case(220) ! 'GROUND_RRV'

   fldval = land%ground_rrv(iland) * 1.e3

case(221) ! 'SOIL_DEPTH'

   fldval = -slz(1)

case(222) ! 'SOIL_WATER_TOT'

   fldval = sum(land%soil_water(:,iland) * dslz(:))

case(223) ! 'HEAD0'

   fldval = land%head0(iland)

case(224) ! 'SLOPE_FACT'

   fldval = land%slope_fact(iland)

!-----------------------------------------
! LAND CELLS - DIF2 fields
!-----------------------------------------

case(230) ! 'WXFERI_DIF2'

   ! NOTE: wxferi_accum fluxes are positive upward; fldval flux (infiltration)
   !       is defined positive downward

   if (.not. allocated(wxferi_accum)) go to 1000
   fldval = -(wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(231) ! 'WXFERP_DIF2'

   ! NOTE: wxferp_accum fluxes are positive upward; fldval flux (percolation)
   !       is defined positive downward

   if (.not. allocated(wxferp_accum)) go to 1000
   fldval = -(wxferp_accum_prev0(iland) - wxferp_accum_prev1(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(232) ! 'WXFER1_DIF2'

   ! NOTE: wxfer1_accum fluxes are positive upward; fldval flux (at soil bottom)
   !       is defined positive downward

   if (.not. allocated(wxfer1_accum)) go to 1000
   fldval = -(wxfer1_accum_prev0(iland) - wxfer1_accum_prev1(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(233) ! 'SOIL_WATER_TOT_DIF2'

   ! Not an accumulated quantity

   fldval = soil_water_tot_prev0(iland) - soil_water_tot_prev1(iland)

case(234) ! 'VEGTEMPK_DMIN_DIF2'

   if (.not. allocated(vegtemp_dmin_accum)) go to 1000
   fldval = (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(235) ! 'VEGTEMPK_DMAX_DIF2'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000
   fldval = (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(236) ! 'VEGTEMPK_DSPAN_DIF2'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000
   fldval = (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland)) &
          - (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(237) ! 'SOILTEMPK_DMIN_DIF2'

   if (.not. allocated(soiltemp_dmin_accum)) go to 1000
   fldval = (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(238) ! 'SOILTEMPK_DMAX_DIF2'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000
   fldval = (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(239) ! 'SOILTEMPK_DSPAN_DIF2'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000
   fldval = (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland)) &
          - (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(240) ! 'WXFERIF_DIF2'

   ! NOTE: wxferi_accum fluxes are positive upward; fldval flux (infiltration)
   !       is defined positive downward

   if (.not. allocated(wxferi_accum)) go to 1000
   fldval = -(wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) * 1.e3 & ! convert from m to mm
          / (pcp_accum_prev0(i) - pcp_accum_prev1(i))

!-----------------------------------------
! LAND CELLS - DIF4 fields
!-----------------------------------------

case(241) ! 'WXFERI_DIF4'

   ! NOTE: wxferi_accum fluxes are positive upward; fldval fluxes (infiltration)
   !       are defined positive downward

   if (.not. allocated(wxferi_accum)) go to 1000
   fldval1 = -(wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) * 1.e3 ! convert from m to mm
   fldval2 = -(wxferi_accum_prev2(iland) - wxferi_accum_prev3(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(242) ! 'WXFERP_DIF4'

   ! NOTE: wxferp_accum fluxes are positive upward; fldval fluxes (percolation)
   !       are defined positive downward

   if (.not. allocated(wxferp_accum)) go to 1000
   fldval1 = -(wxferp_accum_prev0(iland) - wxferp_accum_prev1(iland)) * 1.e3 ! convert from m to mm
   fldval2 = -(wxferp_accum_prev2(iland) - wxferp_accum_prev3(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(243) ! 'WXFER1_DIF4'

   ! NOTE: wxfer1_accum fluxes are positive upward; fldval fluxes (at soil bottom)
   !       are defined positive downward

   if (.not. allocated(wxfer1_accum)) go to 1000
   fldval1 = -(wxfer1_accum_prev0(iland) - wxfer1_accum_prev1(iland)) * 1.e3 ! convert from m to mm
   fldval2 = -(wxfer1_accum_prev2(iland) - wxfer1_accum_prev3(iland)) * 1.e3 ! convert from m to mm

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(244) ! 'VEGTEMPK_DMIN_DIF4'

   if (.not. allocated(vegtemp_dmin_accum)) go to 1000
   fldval1 = (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))
   fldval2 = (vegtemp_dmin_accum_prev2(iland) - vegtemp_dmin_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(245) ! 'VEGTEMPK_DMAX_DIF4'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000
   fldval1 = (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland))
   fldval2 = (vegtemp_dmax_accum_prev2(iland) - vegtemp_dmax_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(246) ! 'VEGTEMPK_DSPAN_DIF4'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000
   fldval1 = (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland)) &
           - (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))
   fldval2 = (vegtemp_dmax_accum_prev2(iland) - vegtemp_dmax_accum_prev3(iland)) &
           - (vegtemp_dmin_accum_prev2(iland) - vegtemp_dmin_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(247) ! 'SOILTEMPK_DMIN_DIF4'

   if (.not. allocated(soiltemp_dmin_accum)) go to 1000
   fldval1 = (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))
   fldval2 = (soiltemp_dmin_accum_prev2(iland) - soiltemp_dmin_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(248) ! 'SOILTEMPK_DMAX_DIF4'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000
   fldval1 = (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland))
   fldval2 = (soiltemp_dmax_accum_prev2(iland) - soiltemp_dmax_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(249) ! 'SOILTEMPK_DSPAN_DIF4'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000
   fldval1 = (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland)) &
           - (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))
   fldval2 = (soiltemp_dmax_accum_prev2(iland) - soiltemp_dmax_accum_prev3(iland)) &
           - (soiltemp_dmin_accum_prev2(iland) - soiltemp_dmin_accum_prev3(iland))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(250) ! 'WXFERIF_DIF4'

   ! NOTE: wxferi_accum fluxes are positive upward; fldval fluxes (infiltration)
   !       are defined positive downward

   if (.not. allocated(wxferi_accum)) go to 1000
   fldval1 = -(wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) * 1.e3 & ! convert from m to mm
           / (pcp_accum_prev0(i) - pcp_accum_prev1(i))
   fldval2 = -(wxferi_accum_prev2(iland) - wxferi_accum_prev3(iland)) * 1.e3 & ! convert from m to mm
           / (pcp_accum_prev2(i) - pcp_accum_prev3(i))

   fldval = fldval1 - fldval2

!--------------------------
! LAND CELLS - ATM averages
!--------------------------

case(251:252) ! 'AL_SFCWATER_TOT', 'AL_SOIL_WATER_TOT'

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      if (trim(fldname) == 'AL_SFCWATER_TOT') then
         fldval = fldval + sum(land%sfcwater_mass(:,iland)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'AL_SOIL_WATER_TOT') then
         fldval = fldval + sum(land%soil_water(:,iland) * dslz(:)) * sfcg%area(iwsfc)
      endif

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

!-----------------------------------------
! LAND CELLS - ATM averages of DIF2 fields
!-----------------------------------------

case(253) ! 'AL_WXFERI_DIF2'

   if (.not. allocated(wxferi_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval - (wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(254) ! 'AL_WXFERP_DIF2'

   if (.not. allocated(wxferp_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval - (wxferp_accum_prev0(iland) - wxferp_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(255) ! 'AL_WXFER1_DIF2'

   if (.not. allocated(wxfer1_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval - (wxfer1_accum_prev0(iland) - wxfer1_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

case(256) ! 'AL_SOIL_WATER_TOT_DIF2'

   ! Not an accumulated quantity

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + (soil_water_tot_prev0(iland) - soil_water_tot_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

case(257) ! 'AL_VEGTEMP_DMIN_DIF2'

   if (.not. allocated(vegtemp_dmin_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(258) ! 'AL_VEGTEMP_DMAX_DIF2'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(259) ! 'AL_VEGTEMP_DSPAN_DIF2'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + ((vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland))  &
                       - (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(260) ! 'AL_SOILTEMP_DMIN_DIF2'

   if (.not. allocated(soiltemp_dmin_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(261) ! 'AL_SOILTEMP_DMAX_DIF2'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland)) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(262) ! 'AL_SOILTEMP_DSPAN_DIF2'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval = fldval + ((soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland))  &
                       - (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))) &
             * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

!-----------------------------------------
! LAND CELLS - ATM averages of DIF4 fields
!-----------------------------------------

case(263) ! 'AL_WXFERI_DIF4'

   if (.not. allocated(wxferi_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 - (wxferi_accum_prev0(iland) - wxferi_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 - (wxferi_accum_prev2(iland) - wxferi_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval1 = fldval1 / area_sum * 1.e3 ! convert from m to mm
      fldval2 = fldval2 / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(264) ! 'AL_WXFERP_DIF4'

   if (.not. allocated(wxferp_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 - (wxferp_accum_prev0(iland) - wxferp_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 - (wxferp_accum_prev2(iland) - wxferp_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval1 = fldval1 / area_sum * 1.e3 ! convert from m to mm
      fldval2 = fldval2 / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(265) ! 'AL_WXFER1_DIF4'

   if (.not. allocated(wxfer1_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 - (wxfer1_accum_prev0(iland) - wxfer1_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 - (wxfer1_accum_prev2(iland) - wxfer1_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval1 = fldval1 / area_sum * 1.e3 ! convert from m to mm
      fldval2 = fldval2 / area_sum * 1.e3 ! convert from m to mm
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1) ! convert from mm to mm/day
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3) ! convert from mm to mm/day
   endif

   fldval = fldval1 - fldval2

case(266) ! 'AL_VEGTEMP_DMIN_DIF4'

   if (.not. allocated(vegtemp_dmin_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 + (vegtemp_dmin_accum_prev2(iland) - vegtemp_dmin_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(267) ! 'AL_VEGTEMP_DMAX_DIF4'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + (vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 + (vegtemp_dmax_accum_prev2(iland) - vegtemp_dmax_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(268) ! 'AL_VEGTEMP_DSPAN_DIF4'

   if (.not. allocated(vegtemp_dmax_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + ((vegtemp_dmax_accum_prev0(iland) - vegtemp_dmax_accum_prev1(iland))  &
                         - (vegtemp_dmin_accum_prev0(iland) - vegtemp_dmin_accum_prev1(iland))) * sfcg%area(iwsfc)
      fldval2 = fldval2 + ((vegtemp_dmax_accum_prev2(iland) - vegtemp_dmax_accum_prev3(iland))  &
                         - (vegtemp_dmin_accum_prev2(iland) - vegtemp_dmin_accum_prev3(iland))) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(269) ! 'AL_SOILTEMP_DMIN_DIF4'

   if (.not. allocated(soiltemp_dmin_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 + (soiltemp_dmin_accum_prev2(iland) - soiltemp_dmin_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(270) ! 'AL_SOILTEMP_DMAX_DIF4'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + (soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland)) * sfcg%area(iwsfc)
      fldval2 = fldval2 + (soiltemp_dmax_accum_prev2(iland) - soiltemp_dmax_accum_prev3(iland)) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(271) ! 'AL_SOILTEMP_DSPAN_DIF4'

   if (.not. allocated(soiltemp_dmax_accum)) go to 1000

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)
      iland = iwsfc - omland

      fldval1 = fldval1 + ((soiltemp_dmax_accum_prev0(iland) - soiltemp_dmax_accum_prev1(iland))  &
                         - (soiltemp_dmin_accum_prev0(iland) - soiltemp_dmin_accum_prev1(iland))) * sfcg%area(iwsfc)
      fldval2 = fldval2 + ((soiltemp_dmax_accum_prev2(iland) - soiltemp_dmax_accum_prev3(iland))  &
                         - (soiltemp_dmin_accum_prev2(iland) - soiltemp_dmin_accum_prev3(iland))) * sfcg%area(iwsfc)

      area_sum = area_sum + sfcg%area(iwsfc)
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

!-----------------------------------------
! SEA CELLS
!-----------------------------------------

case(272) ! 'SEATP'

   fldval = sea%seatp(isea)

case(273) ! 'SEATF'

   fldval = sea%seatf(isea)

case(274) ! 'SEATC'

   fldval = sea%seatc(isea)

case(275) ! 'SEAICEP'

   fldval = sea%seaicep(isea)

case(276) ! 'SEAICEF'

   fldval = sea%seaicef(isea)

case(277) ! 'SEAICEC'

   fldval = sea%seaicec(isea)

case(278) ! 'WDEPTH'

   fldval = sea%wdepth(isea)

case(279) ! 'POM_KBA'

   fldval = pom%kba(isea)

case(280) ! 'POM_TEMPSFC'

   fldval = pom%potmp(1,isea)

!-----------------------------------------
! SFC GRID CELLS - 2D
!-----------------------------------------

case(297) ! 'LEAF_CLASS'

   fldval = real(sfcg%leaf_class(i))

case(298) ! 'SFCG_AREA'

   fldval = sfcg%area(i)

case(299) ! 'SFCG_GLATW'

   fldval = sfcg%glatw(i)

case(300) ! 'SFCG_GLONW'

   fldval = sfcg%glonw(i)

case(301) ! 'SFCG_TOPW'

   fldval = sfcg%topw(i)

case(302) ! 'SFCG_ROUGH'

   fldval = sfcg%rough(i)

case(303) ! 'SFCG_VELS'

   fldval = sfcg%vels(i)

case(304) ! 'SFCG_PRSS'

   fldval = sfcg%prss(i)

case(305) ! 'SFCG_RHOS'

   fldval = sfcg%rhos(i)

case(306) ! 'SFCG_AIRTEMPK'

   fldval = sfcg%airtemp(i)

case(307) ! 'SFCG_AIRRRV'

   fldval = sfcg%airrrv(i) * 1.e3

case(308) ! 'SFCG_CANTEMPK'

   fldval = sfcg%cantemp(i)

case(309) ! 'SFCG_CANRRV'

   fldval = sfcg%canrrv(i) * 1.e3

case(310) ! 'SFCG_SKINTEMPK'

   if (sfcg%leaf_class(i) == 0) then
      isea = i - omsea
      nls = sea%nlev_seaice(isea)

      if (nls > 0) then
         fldval = (1.0 - sea%seaicec(isea)) * sea%seatc(isea) &
                +        sea%seaicec(isea)  * sea%seaice_tempk(nls,isea)
      else
         fldval = sea%seatc(isea)
      endif

   elseif (sfcg%leaf_class(i) == 1) then
      ilake = i - omlake

      call qtk(lake%lake_energy(ilake),tempk,fracliq)

      fldval = tempk

   elseif (sfcg%leaf_class(i) >= 2) then
      iland = i - omland
      nls = land%nlev_sfcwater(iland)

      if (nls > 0) then
         call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)
      else
         call qwtk(land%soil_energy(nzg,iland),        &
                   land%soil_water(nzg,iland)*1.e3,    &
                   land%specifheat_drysoil(nzg,iland), &
                   tempk, fracliq)
      endif

      fldval = (1. - land%vf(iland)) * tempk &
                   + land%vf(iland)  * land%veg_temp(iland)
   endif

case(311) ! 'SFCG_GSS_SRRV'

   if (sfcg%leaf_class(i) == 0) then
      isea = i - omsea
      fldval = sea%surface_srrv(isea) * 1.e3
   elseif (sfcg%leaf_class(i) == 1) then
      ilake = i - omlake
      fldval = lake%surface_srrv(ilake) * 1.e3
   elseif (sfcg%leaf_class(i) >= 2) then
      iland = i - omland
      fldval = land%surface_srrv(iland) * 1.e3
   endif

case(312) ! 'HEAD1'

   fldval = sfcg%head1(i)

case(313) ! 'HEAD_WTAB'

   if (sfcg%leaf_class(i) < 2) then

      fldval = sfcg%head1(i)

   else

      iland = i - omland

      do klev = nzg,1,-1
         call soil_wat2pot(klev, iland, land%soil_water(klev,iland), &
              land%wresid_vg(klev,iland), land%wsat_vg(klev,iland), &
              land%alpha_vg(klev,iland), land%en_vg(klev,iland), psi, psi_slope)

         ! Trial algorithm: Get head_wtab from highest saturated soil level

         if (psi > 1.e-2) then
            head(klev) = psi + slzt(klev)
            fldval = head(klev)
            exit
         else
            head(klev) = psi + slzt(klev)
            if (klev == 1) fldval = head(klev)
         endif
      enddo

   endif

case(314) ! 'SFCG_SENSFLUX'

   fldval = sfcg%sfluxt(i) * cp

case(315) ! 'SFCG_LATFLUX'

   fldval = sfcg%sfluxr(i) * alvl

case(316) ! 'SFCG_VAPFLUX'

   fldval = sfcg%sfluxr(i)

case(317:320) ! 'SFCG_SPEED10M', 'SFCG_SPEED2M', 'SFCG_TEMPK2M', 'SFCG_RVAP2M'

   if (trim(fldname) == 'SFCG_SPEED10M') then
      zobs = 10.
   else
      zobs = 2.
   endif

   press_zobs = sfcg%prss(i) - zobs * sfcg%rhos(i) ! hydrostatic eqn.
   exner_zobs = (press_zobs * p00i) ** rocp

   canexner = (sfcg%prss(i) * p00i) ** rocp
   cantheta  = sfcg%cantemp(i) / canexner
   canthetav = cantheta         * (1.0 + eps_virt * sfcg%canrrv(i))
   airthetav = sfcg%airtheta(i) * (1.0 + eps_virt * sfcg%airrrv(i))

   tstar = -sfcg%sfluxt(i) / (sfcg%ustar(i) * sfcg%rhos(i))
   rstar = -sfcg%sfluxr(i) / (sfcg%ustar(i) * sfcg%rhos(i))

   ufree = (grav * sfcg%dzt_bot(i) * max(sfcg%wthv(i),0.0) / airthetav) ** onethird

   call sfclyr_profile (sfcg%vels(i), sfcg%ustar(i), tstar, rstar, &
                        sfcg%dzt_bot(i), sfcg%rough(i), ufree, &
                        cantheta, canthetav, sfcg%canrrv(i), airthetav, &
                        zobs, wind_zobs, theta_zobs, rrv_zobs)

   if (trim(fldname) == 'SFCG_SPEED10M' .or. trim(fldname) == 'SFCG_SPEED2M') then
      fldval = wind_zobs
   elseif (trim(fldname) == 'SFCG_TEMPK2M') then
      fldval = theta_zobs * exner_zobs
   else
      fldval = rrv_zobs * 1.e3 ! converting from kg/kg to g/kg
   endif

case(321) ! 'SFCG_RSHORT'

   fldval = sfcg%rshort(i)

case(322) ! 'SFCG_RLONG'

   fldval = sfcg%rlong(i)

case(323) ! 'SFCG_RLONGUP'

   fldval = sfcg%rlongup(i)

case(324) ! 'SFCG_RLONG_ALBEDO'

   fldval = sfcg%rlong_albedo(i)

case(325) ! 'SFCG_ALBEDO_BEAM'

   fldval = sfcg%albedo_beam(i)

case(326) ! 'SFCG_ALBEDO_DIFFUSE'

   fldval = sfcg%albedo_diffuse(i)

case(327) ! 'SFCG_BATHYM'

   fldval = sfcg%bathym(i)

case(328) ! 'SFCG_PCPG'

   fldval = sfcg%pcpg(i)

case(329) ! 'SFCG_VC'

   fldval = sfcg%vc(i)

case(330) ! 'WAT_DEPTH'

   if (sfcg%leaf_class(i) == 0) then
      isea = i - omsea
      fldval = sea%wdepth(isea)
   elseif (sfcg%leaf_class(i) == 1) then
      ilake = i - omlake
      fldval = lake%depth(ilake)
   elseif (sfcg%leaf_class(i) >= 2) then
      iland = i - omland
      fldval = 0.
      do klev = 1,land%nlev_sfcwater(iland)
         fldval = fldval + land%sfcwater_depth(klev,iland)
      enddo
   endif

case(331) ! 'WAT_TEMPK'

   if (sfcg%leaf_class(i) == 0) then
      isea = i - omsea
      fldval = sea%seatc(isea)
   elseif (sfcg%leaf_class(i) == 1) then
      ilake = i - omlake
      call qtk(lake%lake_energy(ilake),tempk,fracliq)
      fldval = tempk
   elseif (sfcg%leaf_class(i) >= 2) then
      iland = i - omland
      nls = land%nlev_sfcwater(iland)
      if (nls > 0) then
         call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)
         fldval = tempk
      else
         call qwtk(land%soil_energy(nzg,iland),        &
                   land%soil_water(nzg,iland)*1.e3,    &
                   land%specifheat_drysoil(nzg,iland), &
                   tempk, fracliq)
      endif

      fldval = tempk
   endif

case(332) ! 'HEAD1_MSL'

   fldval = sfcg%head1(i) + sfcg%topw(i)

case(335) ! 'SFCWAT_NUD'

   if (.not. allocated(sfcwat_nud)) go to 1000
   fldval = sfcwat_nud(i) * 86400.

case(336) ! 'SFCTEMP_NUD'

   if (.not. allocated(sfctemp_nud)) go to 1000
   fldval = sfctemp_nud(i)

case(337) ! 'FRACLIQ_NUD'

   if (.not. allocated(fracliq_nud)) go to 1000
   fldval = fracliq_nud(i)

!-----------------------------------------
! SFC GRID CELLS - DIF2 fields
!-----------------------------------------

case(341) ! 'SFCG_VELS_DIF2'

   if (.not. allocated(vels_accum)) go to 1000
   fldval = (vels_accum_prev0(i) - vels_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(342) ! 'SFCG_AIRTEMPK_DIF2'

   if (.not. allocated(airtemp_accum)) go to 1000
   fldval = (airtemp_accum_prev0(i) - airtemp_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(343) ! 'SFCG_AIRRRV_DIF2'

   if (.not. allocated(airrrv_accum)) go to 1000
   fldval = (airrrv_accum_prev0(i) - airrrv_accum_prev1(i)) * 1.e3

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(344) ! 'SFCG_CANTEMPK_DIF2'

   if (.not. allocated(cantemp_accum)) go to 1000
   fldval = (cantemp_accum_prev0(i) - cantemp_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(345) ! 'SFCG_CANRRV_DIF2'

   if (.not. allocated(canrrv_accum)) go to 1000
   fldval = (canrrv_accum_prev0(i) - canrrv_accum_prev1(i)) * 1.e3

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(346) ! 'SFCG_SKINTEMPK_DIF2'

   if (.not. allocated(skintemp_accum)) go to 1000
   fldval = (skintemp_accum_prev0(i) - skintemp_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(347) ! 'SFCG_SENSFLUX_DIF2'

   if (.not. allocated(sfluxt_accum)) go to 1000
   fldval = (sfluxt_accum_prev0(i) - sfluxt_accum_prev1(i)) * cp

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(348) ! 'SFCG_LATFLUX_DIF2'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)) * alvl

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(349) ! 'SFCG_VAPFLUX_DIF2'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

   fldval = fldval * 86400. ! Convert from [kg/(m^2 s)] to [mm/day] for plotting

case(350) ! 'SFCG_PCP_DIF2'

   if (.not. allocated(pcp_accum)) go to 1000
   fldval = (pcp_accum_prev0(i) - pcp_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

   fldval = fldval * 86400. ! Convert from [kg/(m^2 s)] to [mm/day] for plotting

case(351) ! 'SFCG_RUNOFF_DIF2'

   if (.not. allocated(runoff_accum)) go to 1000
   fldval = (runoff_accum_prev0(i) - runoff_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval * 86400. / (time8_prev0 - time8_prev1)
   endif

case(352) ! 'SFCG_SFCTEMP_DIF2'

   if (.not. allocated(sfctemp_accum)) go to 1000
   fldval = (sfctemp_accum_prev0(i) - sfctemp_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(353) ! 'SFCG_FRACLIQ_DIF2'

   if (.not. allocated(fracliq_accum)) go to 1000
   fldval = (fracliq_accum_prev0(i) - fracliq_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(354) ! 'SFCG_AIRTEMPK_DMIN_DIF2'

   if (.not. allocated(airtemp_dmin_accum)) go to 1000
   fldval = (airtemp_dmin_accum_prev0(i) - airtemp_dmin_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(355) ! 'SFCG_AIRTEMPK_DMAX_DIF2'

   if (.not. allocated(airtemp_dmax_accum)) go to 1000
   fldval = (airtemp_dmax_accum_prev0(i) - airtemp_dmax_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(356) ! 'SFCG_AIRTEMPK_DSPAN_DIF2'

   if (.not. allocated(airtemp_dmax_accum)) go to 1000
   fldval = (airtemp_dmax_accum_prev0(i) - airtemp_dmax_accum_prev1(i)) &
          - (airtemp_dmin_accum_prev0(i) - airtemp_dmin_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(357) ! 'SFCG_CANTEMPK_DMIN_DIF2'

   if (.not. allocated(cantemp_dmin_accum)) go to 1000
   fldval = (cantemp_dmin_accum_prev0(i) - cantemp_dmin_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(358) ! 'SFCG_CANTEMPK_DMAX_DIF2'

   if (.not. allocated(cantemp_dmax_accum)) go to 1000
   fldval = (cantemp_dmax_accum_prev0(i) - cantemp_dmax_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(359) ! 'SFCG_CANTEMPK_DSPAN_DIF2'

   if (.not. allocated(cantemp_dmax_accum)) go to 1000
   fldval = (cantemp_dmax_accum_prev0(i) - cantemp_dmax_accum_prev1(i)) &
          - (cantemp_dmin_accum_prev0(i) - cantemp_dmin_accum_prev1(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

case(360) ! 'SFCG_VAPFLUXF_DIF2'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)) / (pcp_accum_prev0(i) - pcp_accum_prev1(i))

case(361) ! 'HEAD_WTAB_DIF2'

   ! Not an accumulated quantity

   fldval = (head_wtab_prev0(i) - head_wtab_prev1(i)) / 120.

!-----------------------------------------
! SFC GRID CELLS - DIF4 fields
!-----------------------------------------

case(362) ! 'SFCG_CANTEMPK_DIF4'

   if (.not. allocated(cantemp_accum)) go to 1000
   fldval1 = (cantemp_accum_prev0(i) - cantemp_accum_prev1(i))
   fldval2 = (cantemp_accum_prev2(i) - cantemp_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(363) ! 'SFCG_CANRRV_DIF4'

   if (.not. allocated(canrrv_accum)) go to 1000
   fldval1 = (canrrv_accum_prev0(i) - canrrv_accum_prev1(i)) * 1.e3
   fldval2 = (canrrv_accum_prev2(i) - canrrv_accum_prev3(i)) * 1.e3

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(364) ! 'SFCG_SKINTEMPK_DIF4'

   if (.not. allocated(skintemp_accum)) go to 1000
   fldval1 = (skintemp_accum_prev0(i) - skintemp_accum_prev1(i))
   fldval2 = (skintemp_accum_prev2(i) - skintemp_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(365) ! 'SFCG_SENSFLUX_DIF4'

   if (.not. allocated(sfluxt_accum)) go to 1000
   fldval1 = (sfluxt_accum_prev0(i) - sfluxt_accum_prev1(i)) * cp
   fldval2 = (sfluxt_accum_prev2(i) - sfluxt_accum_prev3(i)) * cp

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(366) ! 'SFCG_LATFLUX_DIF4'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval1 = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)) * alvl
   fldval2 = (sfluxr_accum_prev2(i) - sfluxr_accum_prev3(i)) * alvl

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(367) ! 'SFCG_VAPFLUX_DIF4'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval1 = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i))
   fldval2 = (sfluxr_accum_prev2(i) - sfluxr_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(368) ! 'SFCG_RUNOFF_DIF4'

   if (.not. allocated(runoff_accum)) go to 1000
   fldval1 = (runoff_accum_prev0(i) - runoff_accum_prev1(i))
   fldval2 = (runoff_accum_prev2(i) - runoff_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 * 86400. / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 * 86400. / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(369) ! 'SFCG_CANTEMPK_DMIN_DIF4'

   if (.not. allocated(cantemp_dmin_accum)) go to 1000
   fldval1 = (cantemp_dmin_accum_prev0(i) - cantemp_dmin_accum_prev1(i))
   fldval2 = (cantemp_dmin_accum_prev2(i) - cantemp_dmin_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(370) ! 'SFCG_CANTEMPK_DMAX_DIF4'

   if (.not. allocated(cantemp_dmax_accum)) go to 1000
   fldval1 = (cantemp_dmax_accum_prev0(i) - cantemp_dmax_accum_prev1(i))
   fldval2 = (cantemp_dmax_accum_prev2(i) - cantemp_dmax_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(371) ! 'SFCG_CANTEMPK_DSPAN_DIF4'

    if (.not. allocated(cantemp_dmax_accum)) go to 1000
    fldval1 = (cantemp_dmax_accum_prev0(i) - cantemp_dmax_accum_prev1(i)) &
            - (cantemp_dmin_accum_prev0(i) - cantemp_dmin_accum_prev1(i))
    fldval2 = (cantemp_dmax_accum_prev2(i) - cantemp_dmax_accum_prev3(i)) &
            - (cantemp_dmin_accum_prev2(i) - cantemp_dmin_accum_prev3(i))

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval1 = fldval1 / (time8_prev0 - time8_prev1)
   endif

   if (abs(time8_prev2 - time8_prev3) > .99) then
      fldval2 = fldval2 / (time8_prev2 - time8_prev3)
   endif

   fldval = fldval1 - fldval2

case(372) ! 'SFCG_VAPFLUXF_DIF4'

   if (.not. allocated(sfluxr_accum)) go to 1000
   fldval1 = (sfluxr_accum_prev0(i) - sfluxr_accum_prev1(i)) / (pcp_accum_prev0(i) - pcp_accum_prev1(i))
   fldval2 = (sfluxr_accum_prev2(i) - sfluxr_accum_prev3(i)) / (pcp_accum_prev2(i) - pcp_accum_prev3(i))

   fldval = fldval1 - fldval2

!-----------------------------------------
! SFC GRID CELLS - ATM averages
!-----------------------------------------

case(401:409) ! 'ASFCG_VELS',       'ASFCG_AIRTEMPK', 'ASFCG_AIRRRV',
              ! 'ASFCG_CANTEMPK',   'ASFCG_CANRRV',
              ! 'ASFCG_SENSFLUX',   'ASFCG_LATFLUX',
              ! 'ASFCG_VAPFLUX',    'ASFCG_SKINTEMPK'

   fldval = 0.
   area_sum = 0.

   do j = 1,itab_w(i)%jsfc2
      iwsfc = itab_w(i)%iwsfc(j)

      if     (trim(fldname) == 'ASFCG_VELS'    ) then
         fldval = fldval + sfcg%vels   (iwsfc) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK') then
         fldval = fldval + sfcg%airtemp(iwsfc) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRRRV'  ) then
         fldval = fldval + sfcg%airrrv (iwsfc) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_CANTEMPK') then
         fldval = fldval + sfcg%cantemp(iwsfc) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANRRV'  ) then
         fldval = fldval + sfcg%canrrv (iwsfc) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_SENSFLUX') then
         fldval = fldval + sfcg%sfluxt (iwsfc) * sfcg%area(iwsfc) * cp
      elseif (trim(fldname) == 'ASFCG_LATFLUX' ) then
         fldval = fldval + sfcg%sfluxr (iwsfc) * sfcg%area(iwsfc) * alvl
      elseif (trim(fldname) == 'ASFCG_VAPFLUX' ) then
         fldval = fldval + sfcg%sfluxr (iwsfc) * sfcg%area(iwsfc)

      elseif (trim(fldname) == 'ASFCG_SKINTEMPK') then

         if (sfcg%leaf_class(iwsfc) >= 2) then

            iland = iwsfc - omland
            nls = land%nlev_sfcwater(iland)

            if (nls > 0) then
               call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)
            else
               call qwtk(land%soil_energy(nzg,iland),        &
                         land%soil_water(nzg,iland)*1.e3,    &
                         land%specifheat_drysoil(nzg,iland), &
                         tempk, fracliq)
            endif

            fldval = fldval + ((1. - land%vf(iland)) * tempk &
                                   + land%vf(iland)  * land%veg_temp(iland)) &
                   * sfcg%area(iwsfc)

         elseif (sfcg%leaf_class(iwsfc) == 1) then

            ilake = iwsfc - omlake

            call qtk(lake%lake_energy(ilake),tempk,fracliq)
            fldval = fldval + tempk * sfcg%area(iwsfc)

         elseif (sfcg%leaf_class(iwsfc) == 0) then

            isea = iwsfc - omsea
            nls = sea%nlev_seaice(isea)

            if (nls > 0) then
               fldval = fldval &
                      + ((1.0 - sea%seaicec(isea)) * sea%seatc(isea) &
                              + sea%seaicec(isea)  * sea%seaice_tempk(nls,isea)) &
                      * sfcg%area(iwsfc)
            else
               fldval = fldval + sea%seatc(isea) * sfcg%area(iwsfc)
            endif

         endif

      else
         cycle
      endif

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

!-------------------------------------------------
! SFC GRID CELLS - ATM averages of DIF2 fields
!-------------------------------------------------

case(421:435) ! 'ASFCG_VELS_DIF2',          'ASFCG_AIRTEMPK_DIF2',      'ASFCG_AIRRRV_DIF2'
              ! 'ASFCG_CANTEMPK_DIF2',      'ASFCG_CANRRV_DIF2',        'ASFCG_SKINTEMPK_DIF2'
              ! 'ASFCG_SENSFLUX_DIF2',      'ASFCG_LATFLUX_DIF2',       'ASFCG_VAPFLUX_DIF2'
              ! 'ASFCG_AIRTEMPK_DMIN_DIF2', 'ASFCG_AIRTEMPK_DMAX_DIF2', 'ASFCG_AIRTEMPK_DSPAN_DIF2'
              ! 'ASFCG_CANTEMPK_DMIN_DIF2', 'ASFCG_CANTEMPK_DMAX_DIF2', 'ASFCG_CANTEMPK_DSPAN_DIF2'

   fldval = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)

      if     (trim(fldname) == 'ASFCG_VELS_DIF2'    ) then
         if (.not. allocated(vels_accum)) go to 1000
         fldval = fldval +  (vels_accum_prev0(iwsfc)  &
                           - vels_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DIF2') then
         if (.not. allocated(airtemp_accum)) go to 1000
         fldval = fldval +  (airtemp_accum_prev0(iwsfc)  &
                           - airtemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRRRV_DIF2'  ) then
         if (.not. allocated(airrrv_accum)) go to 1000
         fldval = fldval +  (airrrv_accum_prev0(iwsfc)  &
                           - airrrv_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DIF2') then
         if (.not. allocated(cantemp_accum)) go to 1000
         fldval = fldval +  (cantemp_accum_prev0(iwsfc)  &
                           - cantemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANRRV_DIF2'  ) then
         if (.not. allocated(canrrv_accum)) go to 1000
         fldval = fldval +  (canrrv_accum_prev0(iwsfc)  &
                           - canrrv_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_SKINTEMPK_DIF2') then
         if (.not. allocated(skintemp_accum)) go to 1000
         fldval = fldval +  (skintemp_accum_prev0(iwsfc)  &
                           - skintemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_SENSFLUX_DIF2') then
         if (.not. allocated(sfluxt_accum)) go to 1000
         fldval = fldval +  (sfluxt_accum_prev0(iwsfc)  &
                           - sfluxt_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * cp
      elseif (trim(fldname) == 'ASFCG_LATFLUX_DIF2' ) then
         if (.not. allocated(sfluxr_accum)) go to 1000
         fldval = fldval +  (sfluxr_accum_prev0(iwsfc)  &
                           - sfluxr_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * alvl
      elseif (trim(fldname) == 'ASFCG_VAPFLUX_DIF2' ) then
         if (.not. allocated(sfluxr_accum)) go to 1000
         fldval = fldval +  (sfluxr_accum_prev0(iwsfc)  &
                           - sfluxr_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DMIN_DIF2') then
         if (.not. allocated(airtemp_dmin_accum)) go to 1000
         fldval = fldval +  (airtemp_dmin_accum_prev0(iwsfc)  &
                           - airtemp_dmin_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DMAX_DIF2') then
         if (.not. allocated(airtemp_dmax_accum)) go to 1000
         fldval = fldval +  (airtemp_dmax_accum_prev0(iwsfc)  &
                           - airtemp_dmax_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DSPAN_DIF2') then
         if (.not. allocated(airtemp_dmax_accum)) go to 1000
         fldval = fldval + ((airtemp_dmax_accum_prev0(iwsfc)   &
                          -  airtemp_dmin_accum_prev0(iwsfc))  &
                          - (airtemp_dmax_accum_prev1(iwsfc)   &
                          -  airtemp_dmin_accum_prev1(iwsfc))) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DMIN_DIF2') then
         if (.not. allocated(cantemp_dmin_accum)) go to 1000
         fldval = fldval +  (cantemp_dmin_accum_prev0(iwsfc)  &
                           - cantemp_dmin_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DMAX_DIF2') then
         if (.not. allocated(cantemp_dmax_accum)) go to 1000
         fldval = fldval +  (cantemp_dmax_accum_prev0(iwsfc)  &
                           - cantemp_dmax_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DSPAN_DIF2') then
         if (.not. allocated(cantemp_dmax_accum)) go to 1000
         fldval = fldval + ((cantemp_dmax_accum_prev0(iwsfc)   &
                          -  cantemp_dmin_accum_prev0(iwsfc))  &
                          - (cantemp_dmax_accum_prev1(iwsfc)   &
                          -  cantemp_dmin_accum_prev1(iwsfc))) * sfcg%area(iwsfc)
      else
         cycle
      endif

      area_sum = area_sum + sfcg%area(iwsfc)
   enddo

   if (area_sum > 1.e-3) then
      fldval = fldval / area_sum
   endif

   if (abs(time8_prev0 - time8_prev1) > .99) then
      fldval = fldval / (time8_prev0 - time8_prev1)
   endif

!-------------------------------------------------
! SFC GRID CELLS - ATM averages of DIF4 fields
!-------------------------------------------------

case(451:465) ! 'ASFCG_VELS_DIF4',          'ASFCG_AIRTEMPK_DIF4',      'ASFCG_AIRRRV_DIF4'
              ! 'ASFCG_CANTEMPK_DIF4',      'ASFCG_CANRRV_DIF4',        'ASFCG_SKINTEMPK_DIF4'
              ! 'ASFCG_SENSFLUX_DIF4',      'ASFCG_LATFLUX_DIF4',       'ASFCG_VAPFLUX_DIF4'
              ! 'ASFCG_AIRTEMPK_DMIN_DIF4', 'ASFCG_AIRTEMPK_DMAX_DIF4', 'ASFCG_AIRTEMPK_DSPAN_DIF4'
              ! 'ASFCG_CANTEMPK_DMIN_DIF4', 'ASFCG_CANTEMPK_DMAX_DIF4', 'ASFCG_CANTEMPK_DSPAN_DIF4'

   fldval1 = 0.
   fldval2 = 0.
   area_sum = 0.

   do j = itab_w(i)%jland1,itab_w(i)%jland2
      iwsfc = itab_w(i)%iwsfc(j)

      if     (trim(fldname) == 'ASFCG_VELS_DIF4'    ) then
         if ( .not. allocated(vels_accum)) go to 1000
         fldval1 = fldval1 + (vels_accum_prev0(iwsfc)  &
                            - vels_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (vels_accum_prev2(iwsfc)  &
                            - vels_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DIF4') then
         if ( .not. allocated(airtemp_accum)) go to 1000
         fldval1 = fldval1 + (airtemp_accum_prev0(iwsfc)  &
                            - airtemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (airtemp_accum_prev2(iwsfc) &
                            - airtemp_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRRRV_DIF4'  ) then
         if ( .not. allocated(airrrv_accum)) go to 1000
         fldval1 = fldval1 + (airrrv_accum_prev0(iwsfc) &
                            - airrrv_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * 1.e3
         fldval2 = fldval2 + (airrrv_accum_prev2(iwsfc)  &
                            - airrrv_accum_prev3(iwsfc)) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DIF4') then
         if ( .not. allocated(cantemp_accum)) go to 1000
         fldval1 = fldval1 + (cantemp_accum_prev0(iwsfc) &
                            - cantemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (cantemp_accum_prev2(iwsfc)  &
                            - cantemp_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANRRV_DIF4'  ) then
         if ( .not. allocated(canrrv_accum)) go to 1000
         fldval1 = fldval1 + (canrrv_accum_prev0(iwsfc)  &
                            - canrrv_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * 1.e3
         fldval2 = fldval2 + (canrrv_accum_prev2(iwsfc)  &
                            - canrrv_accum_prev3(iwsfc)) * sfcg%area(iwsfc) * 1.e3
      elseif (trim(fldname) == 'ASFCG_SKINTEMPK_DIF4') then
         if ( .not. allocated(skintemp_accum)) go to 1000
         fldval1 = fldval1 + (skintemp_accum_prev0(iwsfc)  &
                            - skintemp_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (skintemp_accum_prev2(iwsfc)  &
                            - skintemp_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_SENSFLUX_DIF4') then
         if ( .not. allocated(sfluxt_accum)) go to 1000
         fldval1 = fldval1 + (sfluxt_accum_prev0(iwsfc)  &
                            - sfluxt_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * cp
         fldval2 = fldval2 + (sfluxt_accum_prev2(iwsfc)  &
                            - sfluxt_accum_prev3(iwsfc)) * sfcg%area(iwsfc) * cp
      elseif (trim(fldname) == 'ASFCG_LATFLUX_DIF4' ) then
         if ( .not. allocated(sfluxr_accum)) go to 1000
         fldval1 = fldval1 + (sfluxr_accum_prev0(iwsfc)  &
                            - sfluxr_accum_prev1(iwsfc)) * sfcg%area(iwsfc) * alvl
         fldval2 = fldval2 + (sfluxr_accum_prev2(iwsfc)  &
                            - sfluxr_accum_prev3(iwsfc)) * sfcg%area(iwsfc) * alvl
      elseif (trim(fldname) == 'ASFCG_VAPFLUX_DIF4' ) then
         if ( .not. allocated(sfluxr_accum)) go to 1000
         fldval1 = fldval1 + (sfluxr_accum_prev0(iwsfc)  &
                            - sfluxr_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (sfluxr_accum_prev2(iwsfc)  &
                            - sfluxr_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DMIN_DIF4') then
         if ( .not. allocated(airtemp_dmin_accum)) go to 1000
         fldval1 = fldval1 + (airtemp_dmin_accum_prev0(iwsfc)  &
                            - airtemp_dmin_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (airtemp_dmin_accum_prev2(iwsfc)  &
                            - airtemp_dmin_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DMAX_DIF4') then
         if ( .not. allocated(airtemp_dmax_accum)) go to 1000
         fldval1 = fldval1 + (airtemp_dmax_accum_prev0(iwsfc)  &
                            - airtemp_dmax_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (airtemp_dmax_accum_prev2(iwsfc)  &
                            - airtemp_dmax_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_AIRTEMPK_DSPAN_DIF4') then
         if ( .not. allocated(airtemp_dmax_accum)) go to 1000
         fldval1 = fldval1 + ((airtemp_dmax_accum_prev0(iwsfc)   &
                           -   airtemp_dmax_accum_prev1(iwsfc))  &
                           -  (airtemp_dmin_accum_prev0(iwsfc)   &
                           -   airtemp_dmin_accum_prev1(iwsfc))) * sfcg%area(iwsfc)
         fldval2 = fldval2 + ((airtemp_dmax_accum_prev2(iwsfc)   &
                           -   airtemp_dmax_accum_prev3(iwsfc))  &
                           -  (airtemp_dmin_accum_prev2(iwsfc)   &
                           -   airtemp_dmin_accum_prev3(iwsfc))) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DMIN_DIF4') then
         if ( .not. allocated(cantemp_dmin_accum)) go to 1000
         fldval1 = fldval1 + (cantemp_dmin_accum_prev0(iwsfc)  &
                            - cantemp_dmin_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (cantemp_dmin_accum_prev2(iwsfc)  &
                            - cantemp_dmin_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DMAX_DIF4') then
         if ( .not. allocated(cantemp_dmax_accum)) go to 1000
         fldval1 = fldval1 + (cantemp_dmax_accum_prev0(iwsfc)  &
                            - cantemp_dmax_accum_prev1(iwsfc)) * sfcg%area(iwsfc)
         fldval2 = fldval2 + (cantemp_dmax_accum_prev2(iwsfc)  &
                            - cantemp_dmax_accum_prev3(iwsfc)) * sfcg%area(iwsfc)
      elseif (trim(fldname) == 'ASFCG_CANTEMPK_DSPAN_DIF4') then
         if ( .not. allocated(cantemp_dmax_accum)) go to 1000
         fldval1 = fldval1 + ((cantemp_dmax_accum_prev0(iwsfc)   &
                           -   cantemp_dmax_accum_prev1(iwsfc))  &
                           -  (cantemp_dmin_accum_prev0(iwsfc)   &
                           -   cantemp_dmin_accum_prev1(iwsfc))) * sfcg%area(iwsfc)
         fldval2 = fldval2 + ((cantemp_dmax_accum_prev2(iwsfc)   &
                           -   cantemp_dmax_accum_prev3(iwsfc))  &
                           -  (cantemp_dmin_accum_prev2(iwsfc)   &
                           -   cantemp_dmin_accum_prev3(iwsfc))) * sfcg%area(iwsfc)
      else
         cycle
      endif

      area_sum = area_sum + sfcg%area(iwsfc)
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

case(501) ! 'ARV'
   if (.not. allocated(arv)) go to 1000
   fldval = arv(k,i)
case(502) ! 'ARW'
   fldval = arw(k,i)
case(503) ! 'VOLT'
   fldval = volt(k,i)

! GRID GEOMETRY - 2D

case(511) ! 'TOPM'
   fldval = topm(i)
case(512) ! 'TOPW'
   fldval = topw(i)
case(513) ! 'GLATM'
   fldval = glatm(i)
case(514) ! 'GLONM'
   fldval = glonm(i)
case(515) ! 'GLATV'
   if (.not. allocated(glatv)) go to 1000
   fldval = glatv(i)
case(516) ! 'GLONV'
   if (.not. allocated(glonv)) go to 1000
   fldval = glonv(i)
case(517) ! 'GLATW'
   fldval = glatw(i)
case(518) ! 'GLONW'
   fldval = glonw(i)
case(519) ! 'LPM'
   fldval = real(lpm(i))
case(520) ! 'LPV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(521) ! 'LCV'
   if (.not. allocated(lpv)) go to 1000
   fldval = real(lpv(i))
case(522) ! 'LPW'
   fldval = real(lpw(i))
case(523) ! 'LSW'
   fldval = real(lsw(i))
case(524) ! 'XEM'
   fldval = xem(i)
case(525) ! 'YEM'
   fldval = yem(i)
case(526) ! 'ZEM'
   fldval = zem(i)
case(527) ! 'XEV'
   if (.not. allocated(xev)) go to 1000
   fldval = xev(i)
case(528) ! 'YEV'
   if (.not. allocated(yev)) go to 1000
   fldval = yev(i)
case(529) ! 'ZEV'
   if (.not. allocated(zev)) go to 1000
   fldval = zev(i)
case(530) ! 'XEW'
   fldval = xew(i)
case(531) ! 'YEW'
   fldval = yew(i)
case(532) ! 'ZEW'
   fldval = zew(i)
case(533) ! 'DNU'
   fldval = dnu(i)
case(534) ! 'DNV'
   fldval = dnv(i)
case(535) ! 'ARM0'
   fldval = arm0(i)
case(536) ! 'ARW0'
   fldval = arw0(i)
case(537) ! 'TOPMW'
   ! fldval filled directly in subroutine contslab_topmw

! ITAB_M MEMBERS

case(541) ! 'ITAB_M_NPOLY'
   fldval = itab_m(i)%npoly
case(542) ! 'ITAB_M_IMGLOBE'
   fldval = itab_m(i)%imglobe
case(543) ! 'ITAB_M_MRLM'
   fldval = itab_m(i)%mrlm
case(544) ! 'ITAB_M_MRLM_OR'
   fldval = itab_m(i)%mrlm_orig
case(545) ! 'ITAB_M_MROW'
   fldval = itab_m(i)%mrow
case(546) ! 'ITAB_M_NGR'
   fldval = itab_m(i)%ngr
case(547) ! 'ITAB_M_IV'
   fldval = itab_m(i)%iv(indp)
case(548) ! 'ITAB_M_IW'
   fldval = itab_m(i)%iw(indp)
case(549) ! 'ITAB_M_IM'
   fldval = itab_m(i)%im(indp)
case(550) ! 'ITAB_M_IRANK'
   fldval = itab_m(i)%irank

! ITAB_V MEMBERS

case(551) ! 'ITAB_V_IVP'
   fldval = itab_v(i)%ivp
case(552) ! 'ITAB_V_IRANK'
   fldval = itab_v(i)%irank
!   fldval = itabg_v(i)%irank
case(553) ! 'ITAB_V_IVGLOBE'
   fldval = itab_v(i)%ivglobe
case(554) ! 'ITAB_V_MRLV'
   fldval = itab_v(i)%mrlv
case(555) ! 'ITAB_V_FARW'
   fldval = itab_v(i)%farw(indp)
case(556) ! 'ITAB_V_IM'
   fldval = itab_v(i)%im(indp)
case(557) ! 'ITAB_V_IW'
   fldval = itab_v(i)%iw(indp)
!case(558) ! 'ITAB_V_IV'
!   fldval = itab_v(i)%iv(indp)

! ITAB_W MEMBERS

case(561) ! 'ITAB_W_NPOLY'
   fldval = itab_w(i)%npoly
case(562) ! 'ITAB_W_IWP'
   fldval = itab_w(i)%iwp
case(563) ! 'ITAB_W_IRANK'
   fldval = itab_w(i)%irank
case(564) ! 'ITAB_W_IWGLOBE'
   fldval = itab_w(i)%iwglobe
case(565) ! 'ITAB_W_MRLW'
   fldval = itab_w(i)%mrlw
case(566) ! 'ITAB_W_MRLW_OR'
   fldval = itab_w(i)%mrlw_orig
case(567) ! 'ITAB_W_NGR'
   fldval = itab_w(i)%ngr
case(568) ! 'ITAB_W_IM'
   fldval = itab_w(i)%im(indp)
case(569) ! 'ITAB_W_IV'
   fldval = itab_w(i)%iv(indp)
case(570) ! 'ITAB_W_IW'
   fldval = itab_w(i)%iw(indp)
case(571) ! 'ITAB_W_DIRV'
   fldval = itab_w(i)%dirv(indp)
case(572) ! 'ITAB_W_FARM'
   fldval = itab_w(i)%farm(indp)
case(573) ! 'ITAB_W_FARV'
   fldval = itab_w(i)%farv(indp)
case(574) ! 'ITAB_W_IWNUD'
   fldval = itab_w(i)%iwnud(indp)
case(575) ! 'ITAB_W_FNUD'
   fldval = itab_w(i)%fnud(indp)
case(576) ! 'ITAB_W_JLAND1'
   fldval = itab_w(i)%jland1
case(577) ! 'ITAB_W_JLAND2'
   fldval = itab_w(i)%jland2
case(578) ! 'ITAB_W_JLAKE1'
   fldval = itab_w(i)%jlake1
case(579) ! 'ITAB_W_JLAKE2'
   fldval = itab_w(i)%jlake2
case(580) ! 'ITAB_W_JSEA1'
   fldval = itab_w(i)%jsea1
case(581) ! 'ITAB_W_JSEA2'
   fldval = itab_w(i)%jsea2
case(582) ! 'ITAB_W_JSFC2'
   fldval = itab_w(i)%jsfc2
case(583) ! 'ITAB_W_IWSFC'
   fldval = itab_w(i)%iwsfc(indp)
case(584) ! 'ITAB_W_JASFC'
   fldval = itab_w(i)%jasfc(indp)

! ITAB_WSFC MEMBERS

case(602) ! 'ITAB_WSFC_IWGLOBE'
   fldval = itab_wsfc(i)%iwglobe
case(603) ! 'ITAB_WSFC_IRANK'
   fldval = itab_wsfc(i)%irank
case(604) ! 'ITAB_WSFC_NWATM'
   fldval = itab_wsfc(i)%nwatm
case(605) ! 'ITAB_WSFC_IWATM'
   iw = itab_wsfc(i)%iwatm(indp) ! local index
   fldval = itab_w(iw)%iwglobe   ! global index
case(606) ! 'ITAB_WSFC_KWATM'
   fldval = itab_wsfc(i)%kwatm(indp)
case(607) ! 'ITAB_WSFC_ARC'
   fldval = itab_wsfc(i)%arc(indp)
case(608) ! 'ITAB_WSFC_ARCOARSFC'
   fldval = itab_wsfc(i)%arcoarsfc(indp)
case(609) ! 'ITAB_WSFC_ARCOARIW'
   fldval = itab_wsfc(i)%arcoariw(indp)
case(610) ! 'ITAB_WSFC_ARCOARKW'
   fldval = itab_wsfc(i)%arcoarkw(indp)
case(611) ! 'ITAB_WSFC_NPOLY'
   fldval = itab_wsfc(i)%npoly
case(612) ! 'ITAB_WSFC_IMN'
   im = itab_wsfc(i)%imn(indp)
   if (im >= 2) then
      fldval = itab_msfc(im)%imglobe
   else
      notavail = 4
   endif
case(613) ! 'ITAB_WSFC_IVN'
   iv = itab_wsfc(i)%ivn(indp)
   if (iv >= 2) then
      fldval = itab_vsfc(iv)%ivglobe
   else
      notavail = 4
   endif
case(614) ! 'ITAB_WSFC_IWN'
   iw = itab_wsfc(i)%iwn(indp)
   if (iw >= 2) then
      fldval = itab_wsfc(iw)%iwglobe
   else
      notavail = 4
   endif

case(722) ! 'PRESS_DAVG'

   if (.not. allocated(press_davg)) go to 1000

   fldval = press_davg(i) * .01

case(723:724) ! 'ZONAL_WINDW_DAVG','MERID_WINDW_DAVG'

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
         raxisi = 1.0 / raxis

         u = (vy * xew(i) - vx * yew(i)) * raxisi
         v = vz * raxis * eradi &
           - (vx * xew(i) + vy * yew(i)) * zew(i) * raxisi * eradi
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

case(725) ! 'RSHORT_DAVG'

   if (.not. allocated(rshort_davg)) go to 1000

   fldval = rshort_davg(i)

case(726) ! 'TEMPK_DAVG'

   if (.not. allocated(tempk_davg)) go to 1000

   fldval = tempk_davg(i)

case(727) ! 'TEMPK_DMIN'

   if (.not. allocated(tempk_dmin)) go to 1000

   fldval = tempk_dmin(i)

case(728) ! 'TEMPK_DMAX'

   if (.not. allocated(tempk_dmax)) go to 1000

   fldval = tempk_dmax(i)

case(729) ! 'ACCPMIC_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000

   fldval = accpmic_dtot(i)

case(730) ! 'ACCPCON_DTOT'

   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpcon_dtot(i)

case(731) ! 'ACCPBOTH_DTOT'

   if (.not. allocated(accpmic_dtot)) go to 1000
   if (.not. allocated(accpcon_dtot)) go to 1000

   fldval = accpmic_dtot(i) + accpcon_dtot(i)

case(732) ! 'PRESS_UL_DAVG'

   if (.not. allocated(press_ul_davg)) go to 1000

   fldval = press_ul_davg(i) * .01

case(733:734) ! 'ZONAL_WINDW_UL_DAVG','MERID_WINDW_UL_DAVG'

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
         raxisi = 1.0 / raxis

         u = (vy * xew(i) - vx * yew(i)) * raxisi
         v = vz * raxis * eradi &
           - (vx * xew(i) + vy * yew(i)) * zew(i) * raxisi * eradi
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

case(735) ! 'CANTEMPK_DAVG'

   if (.not. allocated(cantempk_davg)) go to 1000
   fldval = cantempk_davg(i)

case(736) ! 'CANTEMPK_DMIN'

   if (.not. allocated(cantempk_dmin)) go to 1000
   fldval = cantempk_dmin(i)

case(737) ! 'CANTEMPK_DMAX'

   if (.not. allocated(cantempk_dmax)) go to 1000
   fldval = cantempk_dmax(i)

case(738) ! 'VEGTEMPK_DAVG'

   if (.not. allocated(vegtempk_davg)) go to 1000
   fldval = vegtempk_davg(iland)

case(739) ! 'VEGTEMPK_DMIN'

   if (.not. allocated(vegtempk_dmin)) go to 1000
   fldval = vegtempk_dmin(iland)

case(740) ! 'VEGTEMPK_DMAX'

   if (.not. allocated(vegtempk_dmax)) go to 1000
   fldval = vegtempk_dmax(iland)

case(741) ! 'SOILTEMPK_DAVG'

   if (.not. allocated(soiltempk_davg)) go to 1000
   fldval = soiltempk_davg(iland)

case(742) ! 'SOILTEMPK_DMIN'

   if (.not. allocated(soiltempk_dmin)) go to 1000
   fldval = soiltempk_dmin(iland)

case(743) ! 'SOILTEMPK_DMAX'

   if (.not. allocated(soiltempk_dmax)) go to 1000
   fldval = soiltempk_dmax(iland)

case(744) ! 'SFCG_AIRTEMPK_DAVG'

      if (.not. allocated(airtempk_davg)) go to 1000
      fldval = airtempk_davg(i)

case(745) ! 'SFCG_AIRTEMPK_DMIN'

      if (.not. allocated(airtempk_dmin)) go to 1000
      fldval = airtempk_dmin(i)

case(746) ! 'SFCG_AIRTEMPK_DMAX'

      if (.not. allocated(airtempk_dmax)) go to 1000
      fldval = airtempk_dmax(i)

case(747) ! 'SFCG_CANTEMPK_DAVG'

      if (.not. allocated(cantempk_davg)) go to 1000
      fldval = cantempk_davg(i)

case(748) ! 'SFCG_CANTEMPK_DMIN'

      if (.not. allocated(cantempk_dmin)) go to 1000
      fldval = cantempk_dmin(i)

case(749) ! 'SFCG_CANTEMPK_DMAX'

      if (.not. allocated(cantempk_dmax)) go to 1000
      fldval = cantempk_dmax(i)

case(750) ! 'SFCG_SENSFLUX_DAVG'

      if (.not. allocated(sfluxt_davg)) go to 1000
      fldval = sfluxt_davg(i) * cp

case(751) ! 'SFCG_LATFLUX_DAVG'

      if (.not. allocated(sfluxr_davg)) go to 1000
      fldval = sfluxr_davg(i) * alvl

case(752) ! 'SENSFLUX_DAVG'

   if (.not. allocated(sfluxt_davg)) go to 1000
   fldval = 0.0

   do j = 1,itab_w(i)%jsfc2
      iwsfc = itab_w(i)%iwsfc(j)
      jasfc = itab_w(i)%jasfc(j)
      fldval = fldval + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfluxt_davg(iwsfc) * cp
   enddo

case(753) ! 'LATFLUX_DAVG'

   if (.not. allocated(sfluxr_davg)) go to 1000
   fldval = 0.0

   do j = 1,itab_w(i)%jsfc2
      iwsfc = itab_w(i)%iwsfc(j)
      jasfc = itab_w(i)%jasfc(j)
      fldval = fldval + itab_wsfc(iwsfc)%arcoariw(jasfc) * sfluxr_davg(iwsfc) * alvl
   enddo

case(771) ! 'RHO_OBS'

   if (.not. allocated(rho_obs)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1))

case(772) ! 'THETA_OBS'

   if (.not. allocated(theta_obs)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i)

case(773) ! 'RRW_OBS'

   if (.not. allocated(rrw_obs)) go to 1000

   fldval = rrw_obs(k,itab_w(i)%iwnud(1)) * 1.e3

case(774) ! 'UZONAL_OBS'

   if (.not. allocated(uzonal_obs)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1))

case(775) ! 'UMERID_OBS'

   if (.not. allocated(umerid_obs)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1))

case(776) ! 'RHO_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_sim(k,itab_w(i)%iwnud(1))

case(777) ! 'THETA_SIM'

   if (.not. allocated(theta_sim)) go to 1000

!   fldval = theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_sim(k,i)

case(778) ! 'RRW_SIM'

   if (.not. allocated(rrw_sim)) go to 1000

   fldval = rrw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(779) ! 'UZONAL_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_sim(k,itab_w(i)%iwnud(1))

case(780) ! 'UMERID_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_sim(k,itab_w(i)%iwnud(1))

case(781) ! 'RHO_OBS_SIM'

   if (.not. allocated(rho_sim)) go to 1000

   fldval = rho_obs(k,itab_w(i)%iwnud(1)) &
          - rho_sim(k,itab_w(i)%iwnud(1))

case(782) ! 'THETA_OBS_SIM'

   if (.not. allocated(theta_sim)) go to 1000

!   fldval = theta_obs(k,itab_w(i)%iwnud(1)) &
!          - theta_sim(k,itab_w(i)%iwnud(1))
   fldval = theta_obs(k,i) &
          - theta_sim(k,i)

case(783) ! 'RRW_OBS_SIM'

   if (.not. allocated(rrw_sim)) go to 1000

   fldval = rrw_obs(k,itab_w(i)%iwnud(1)) * 1.e3 &
          - rrw_sim(k,itab_w(i)%iwnud(1)) * 1.e3

case(784) ! 'UZONAL_OBS_SIM'

   if (.not. allocated(uzonal_sim)) go to 1000

   fldval = uzonal_obs(k,itab_w(i)%iwnud(1)) &
          - uzonal_sim(k,itab_w(i)%iwnud(1))

case(785) ! 'UMERID_OBS_SIM'

   if (.not. allocated(umerid_sim)) go to 1000

   fldval = umerid_obs(k,itab_w(i)%iwnud(1)) &
          - umerid_sim(k,itab_w(i)%iwnud(1))

case(786) ! 'VXE'

   if (.not. allocated(vxe)) go to 1000

   fldval = wtbot * vxe(k ,i) &
          + wttop * vxe(kp,i)

case(787) ! 'VYE'

   if (.not. allocated(vye)) go to 1000

   fldval = wtbot * vye(k ,i) &
          + wttop * vye(kp,i)

case(788) ! 'VZE'

   if (.not. allocated(vze)) go to 1000

   fldval = wtbot * vze(k ,i) &
          + wttop * vze(kp,i)

case(789) ! 'PBLH'

   if (.not. allocated(pblh)) go to 1000

   fldval = pblh(i)

case(790) ! 'VKH'

   if (.not. allocated(vkh)) go to 1000

   fldval = wtbot * vkh(k ,i) / rho(k ,i) &
          + wttop * vkh(kp,i) / rho(kp,i)

case(791) ! 'RRW_HCONV'

   fldval = 0.

   npoly = itab_w(i)%npoly

   do klev = lpw(i),mza
      do j = 1,npoly

         iv = itab_w(i)%iv(j)
         iw = itab_w(i)%iw(j)

         fldval = fldval &
                + vmc(klev,iv) * itab_w(i)%dirv(j) * arv(klev,iv) &
                * 0.5 * (rr_w(klev,i) + rr_w(klev,iw))
      enddo
   enddo

   fldval = fldval / arw0(i)

case(792) ! 'RRV_HCONV'

   fldval = 0.

   npoly = itab_w(i)%npoly

   do klev = lpw(i),mza
      do j = 1,npoly

         iv = itab_w(i)%iv(j)
         iw = itab_w(i)%iw(j)

         fldval = fldval &
                + vmc(klev,iv) * itab_w(i)%dirv(j) * arv(klev,iv) &
                * 0.5 * (rr_v(klev,i) + rr_v(klev,iw))
      enddo
   enddo

   fldval = fldval / arw0(i)

case(793) ! 'CLDNUM'

   fldval = cldnum(i) * 1.e-6

case(794) ! 'ADDSC1_ZINT'

   if (naddsc < 1) go to 1000
   if (.not. allocated(addsc(1)%sclp)) go to 1000

   ! Vertical integral of addsc1

   fldval = 0.
   denom = 0.

   do k = lpw(i),mza-1
      fldval = fldval + addsc(1)%sclp(k,i) * dzt(k) * rho(k,i)
      denom = denom + dzt(k) * rho(k,i)
   enddo

   fldval = fldval / denom

   ! special for dcmip 2016:

   fldval = fldval / 4.0e-6

case(795) ! 'ADDSC2_ZINT'

   if (naddsc < 2) go to 1000
   if (.not. allocated(addsc(2)%sclp)) go to 1000

   ! Vertical integral of addsc2

   fldval = 0.
   denom = 0.

   do k = lpw(i),mza-1
      fldval = fldval + addsc(2)%sclp(k,i) * dzt(k) * rho(k,i)
      denom = denom + dzt(k) * rho(k,i)
   enddo

   fldval = fldval / denom

   ! special for dcmip 2016:

   fldval = fldval / 4.0e-6


case(796) ! 'ADDSC1P2_ZINT'

   if (naddsc < 2) go to 1000
   if (.not. allocated(addsc(2)%sclp)) go to 1000

   ! Vertical integral of (addsc1 + 2 * addsc2) - SHOULD WE USE RHO WEIGHTING FOR DCMIP 2016?

   fldval = 0.
   denom = 0.

   do k = lpw(i),mza-1
      fldval = fldval + (addsc(1)%sclp(k,i) + 2. * addsc(2)%sclp(k,i)) * dzt(k) ! * rho(k,i)
      denom = denom + dzt(k) ! * rho(k,i)
   enddo

   fldval = fldval / denom

   ! special for dcmip 2016:

   fldval = fldval / 4.0e-6

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

end subroutine oplot_lib
