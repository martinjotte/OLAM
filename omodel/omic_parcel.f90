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
subroutine parcel_env(iw)

! This subroutine specifies environmental parameters (pressure, etc.) as
! a function of time for parcel simulations

  use mem_basic,   only: thil, press, wc, rho, sh_w, sh_v, theta

  use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                         q2, q6, q7

  use misc_coms,   only: time8, timmax8, dtlm
  use micro_coms,  only: cfmasi, pwmasi, rxmin
  use consts_coms, only: cvocp, p00k, rdry, rvap, p00i, rocp, cpi
  use oname_coms,  only: nl

  implicit none

  integer, intent(in) :: iw
  integer :: k

  real :: thildif, addrain, addhail, exner, tairk

  k = 2

  if (nl%test_case == 901) then

! NOTES:
! This is first test case from Walko et al. (2000).  A parcel initially at a
! temperature of 14 deg C and relative humidity of 80% ascends at
! approximately 6 m/s for 2000 s.  During this time, the parcel cools and
! saturates, producing cloud water and later pristine ice crystals, the only
! two condensate categories allowed in this simulation.  (The drizzle category
! was subsequently added to microphysics, and drizzle can be included in this
! test.)  Collisions and sedimentation are turned off (by IF statements in
! subroutine micphys in omic_driv.f90) in order to examine conservation
! properties of the parcel. No hydrometeor collisons occur (save for
! cloud-cloud and cloud-drizzle for cases where drizzle is included) because
! of the imposed absence of most hydrometeor categories.

! OLAMIN microphysics settings are: ICLOUD = 4, IDRIZ = 5, IRAIN = 0,
! IPRIS = 5, ISNOW = 0, IAGGR = 0, IGRAUP = 0, IHAIL = 0,
! CCNPARM = 1000.E6, GCCNPARM = 0.
! Sounding levels are:
!                 0.0,   15.0,   80.0,   0.0,   0.0,  ! P1
!              1000.0,    8.5,   80.0,   0.0,   0.0,  ! P1

     press(k,iw) = press(k,iw) * .995

     rho(k,iw) = press(k,iw) ** cvocp * p00k &
               / (theta(k,iw) * (1. - sh_c(k,iw)) &
               * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))

  elseif (nl%test_case == 902) then

! NOTES:
! This is second test case from Walko et al. (2000).  The parcel is initially
! at 18 deg C and 30% relative humidity, and remains at constant height.  On
! the first timestep, 6 g/kg of rain is added suddenly, as if falling into the
! parcel from above.  The rain is not permitted to fall out of the parcel
! (prevented by an IF statement in subroutine micphys in omic_driv.f90), but
! it does keep its ventilation coefficient as if falling at its natural
! velocity.  The rain is given an initial temperature of 0 deg C as if it were
! just shed or melted from hail.  Hail with a mean-mass diameter of 2 mm and
! mixing ratio 6 g/kg is introduced into the parcel at 1500 s, and is prevented
! from falling out of the parcel (by the IF statement in subroutine micphys;
! collisions between hydrometeors are also turned off by another IF statement
! in subroutine micphys.)  The internal energy of hail is initially 0 J/kg,
! representing a temperature of 0 deg C and 1000% ice content.

! OLAMIN microphysics settings are: ICLOUD = 4, IDRIZ = 5, IRAIN = 2,
! IPRIS = 0, ISNOW = 0, IAGGR = 0, IGRAUP = 0, IHAIL = 2, 
! CCNPARM = 1000.E6, GCCNPARM = 10., HPARM = 2.0e-3
! Sounding levels are:
!               0.0,   20.0,   30.0,   0.0,   0.0,  ! P2
!            1000.0,   13.5,   30.0,   0.0,   0.0,  ! P2

     exner = (press(2,iw) * p00i) ** rocp  ! defined WITHOUT CP factor
     tairk = theta(2,iw) * exner

     addrain = .006

     if (real(time8) < .1) then
        sh_w(k,iw) = sh_w(k,iw) + addrain
        sh_r(k,iw) = sh_r(k,iw) + addrain
        q2(k,iw) = 334000. ! + 24. * 4186.

        thildif = - thil(k,iw) * thil(k,iw) &
                * (2820. * addrain - cpi * (q2(k,iw) * sh_r(k,iw) ))  &
                / (max(tairk, 253.) * theta(k,iw))

        thil(k,iw) = thil(k,iw) + thildif
     endif

     addhail = .006

     if (real(time8) >= 1499. .and. real(time8) <= 1505.) then
        sh_w(k,iw) = sh_w(k,iw) + addhail
        sh_h(k,iw) = sh_h(k,iw) + addhail
        q7(k,iw) = 0.

        thildif = - thil(k,iw) * thil(k,iw) &
           * (2820. * addhail - cpi * (q7(k,iw) * sh_h(k,iw) )) &
           / (max(tairk, 253.) * theta(k,iw))
 
        thil(k,iw) = thil(k,iw) + thildif
     endif

  endif

end subroutine parcel_env

!===========================================================================

subroutine parcel_plot(mza0,iw0,ncat,dtli0,jhcat,rx,cx,emb,qx,tx,vap, &
    press0,thil0,theta0,tairc,rhovslair,rhovsiair,rhov,rhoi,rhoa,rhow, &
    rnuc_vc,rnuc_vd,rnuc_cp_hom,rnuc_dp_hom,rnuc_vp_haze,rnuc_vp_immers, &
    cnuc_vc,cnuc_vd,cnuc_cp_hom,cnuc_dp_hom,cnuc_vp_haze,cnuc_vp_immers, &
    rpsxfer,epsxfer, &
    r1118,r8882,r1112,r1818,r1212,r8282,r3335,r4445,r3435,r3445, &
    r3535,r3636,r3737,r4545,r4646,r4747,r5656,r5757,r6767,r1413, &
    r1416,r1414,r1446,r1513,r1516,r1515,r1556,r1613,r1616,r8483, &
    r8486,r8484,r8446,r8583,r8586,r8585,r8556,r8683,r8686,r1713, &
    r1717,r8783,r8787,r2332,r2327,r2323,r2337,r2442,r2427,r2424, &
    r2447,r2552,r2527,r2525,r2557,r2662,r2627,r2626,r2667,r2772, &
    r2727,r0000, &
    e1111,e1118,e8888,e8882,e1112,e1811,e1211,e8288,e2222,e5555, &
    e6666,e7777,e3333,e3335,e4444,e4445,e3433,e3444,e3435,e3445, &
    e3533,e3633,e3733,e4544,e4644,e4744,e5655,e5755,e6766,e1413, &
    e1411,e1446,e1513,e1511,e1556,e1613,e1611,e8483,e8488,e8446, &
    e8583,e8588,e8556,e8683,e8688,e1713,e1711,e8783,e8788,e2322, &
    e2327,e2333,e2337,e2422,e2427,e2444,e2447,e2522,e2527,e2555, &
    e2557,e2622,e2627,e2666,e2667,e2722,e2727,e2777,e0000)

  use consts_coms, only: r8
  use misc_coms,   only: time8, timmax8, dtlm
  use micro_coms,  only: cfmasi, pwmasi, rxmin

  implicit none

  integer, intent(in) :: mza0, iw0, ncat

  integer, intent(in) :: jhcat(mza0,ncat)

  real, intent(in) :: dtli0

  real, intent(in) :: rx (mza0,ncat)
  real, intent(in) :: cx (mza0,ncat)
  real, intent(in) :: emb(mza0,ncat)
  real, intent(in) :: qx (mza0,ncat)
  real, intent(in) :: tx (mza0,ncat)
  real, intent(in) :: vap(mza0,ncat)

  real, intent(in) :: press0   (mza0)
  real, intent(in) :: thil0    (mza0)
  real, intent(in) :: theta0   (mza0)
  real, intent(in) :: tairc    (mza0)
  real, intent(in) :: rhovslair(mza0)
  real, intent(in) :: rhovsiair(mza0)
  real, intent(in) :: rhov     (mza0)
  real, intent(in) :: rhoi     (mza0)

  real(r8), intent(in) :: rhoa(mza0)
  real(r8), intent(in) :: rhow(mza0)

  real, intent(in) :: rnuc_vc       (mza0)
  real, intent(in) :: rnuc_vd       (mza0)
  real, intent(in) :: rnuc_cp_hom   (mza0)
  real, intent(in) :: rnuc_dp_hom   (mza0)
  real, intent(in) :: rnuc_vp_haze  (mza0)
  real, intent(in) :: rnuc_vp_immers(mza0)

  real, intent(in) :: cnuc_vc       (mza0)
  real, intent(in) :: cnuc_vd       (mza0)
  real, intent(in) :: cnuc_cp_hom   (mza0)
  real, intent(in) :: cnuc_dp_hom   (mza0)
  real, intent(in) :: cnuc_vp_haze  (mza0)
  real, intent(in) :: cnuc_vp_immers(mza0)

  real, intent(in) :: rpsxfer(mza0)
  real, intent(in) :: epsxfer(mza0)

  real, intent(in) :: &
  r1118(mza0,2),r8882(mza0,2),r1112(mza0,2),r1818(mza0,2),r1212(mza0,2), &
  r8282(mza0,2),r3335(mza0,2),r4445(mza0,2),r3435(mza0,2),r3445(mza0,2), &
  r3535(mza0,2),r3636(mza0,2),r3737(mza0,2),r4545(mza0,2),r4646(mza0,2), &
  r4747(mza0,2),r5656(mza0,2),r5757(mza0,2),r6767(mza0,2),r1413(mza0,2), &
  r1416(mza0,2),r1414(mza0,2),r1446(mza0,2),r1513(mza0,2),r1516(mza0,2), &
  r1515(mza0,2),r1556(mza0,2),r1613(mza0,2),r1616(mza0,2),r8483(mza0,2), &
  r8486(mza0,2),r8484(mza0,2),r8446(mza0,2),r8583(mza0,2),r8586(mza0,2), &
  r8585(mza0,2),r8556(mza0,2),r8683(mza0,2),r8686(mza0,2),r1713(mza0,2), &
  r1717(mza0,2),r8783(mza0,2),r8787(mza0,2),r2332(mza0,2),r2327(mza0,2), &
  r2323(mza0,2),r2337(mza0,2),r2442(mza0,2),r2427(mza0,2),r2424(mza0,2), &
  r2447(mza0,2),r2552(mza0,2),r2527(mza0,2),r2525(mza0,2),r2557(mza0,2), &
  r2662(mza0,2),r2627(mza0,2),r2626(mza0,2),r2667(mza0,2),r2772(mza0,2), &
  r2727(mza0,2),r0000(mza0,2)

  real, intent(in) :: &
  e1111(mza0),e1118(mza0),e8888(mza0),e8882(mza0),e1112(mza0), &
  e1811(mza0),e1211(mza0),e8288(mza0),e2222(mza0),e5555(mza0), &
  e6666(mza0),e7777(mza0),e3333(mza0),e3335(mza0),e4444(mza0), &
  e4445(mza0),e3433(mza0),e3444(mza0),e3435(mza0),e3445(mza0), &
  e3533(mza0),e3633(mza0),e3733(mza0),e4544(mza0),e4644(mza0), &
  e4744(mza0),e5655(mza0),e5755(mza0),e6766(mza0),e1413(mza0), &
  e1411(mza0),e1446(mza0),e1513(mza0),e1511(mza0),e1556(mza0), &
  e1613(mza0),e1611(mza0),e8483(mza0),e8488(mza0),e8446(mza0), &
  e8583(mza0),e8588(mza0),e8556(mza0),e8683(mza0),e8688(mza0), &
  e1713(mza0),e1711(mza0),e8783(mza0),e8788(mza0),e2322(mza0), &
  e2327(mza0),e2333(mza0),e2337(mza0),e2422(mza0),e2427(mza0), &
  e2444(mza0),e2447(mza0),e2522(mza0),e2527(mza0),e2555(mza0), &
  e2557(mza0),e2622(mza0),e2627(mza0),e2666(mza0),e2667(mza0), &
  e2722(mza0),e2727(mza0),e2777(mza0),e0000(mza0)

  real, save, allocatable :: press0_s   (:)
  real, save, allocatable :: thil0_s    (:)
  real, save, allocatable :: theta0_s   (:)
  real, save, allocatable :: tairc_s    (:)
  real, save, allocatable :: rhoa_s     (:)
  real, save, allocatable :: rhoi_s     (:)
  real, save, allocatable :: rhow_s     (:)
  real, save, allocatable :: rhov_s     (:)
  real, save, allocatable :: rhovslair_s(:)
  real, save, allocatable :: rhovsiair_s(:)
  real, save, allocatable :: rhl_s      (:)
  real, save, allocatable :: rhi_s      (:)
  real, save, allocatable :: hcat_s     (:)

  real, save, allocatable :: rnuc_vc_s  (:)
  real, save, allocatable :: rnuc_vd_s  (:)
  real, save, allocatable :: cnuc_vc_s  (:)
  real, save, allocatable :: cnuc_vd_s  (:)

  real, save, allocatable :: rnuc_cp_hom_s (:)
  real, save, allocatable :: rnuc_dp_hom_s (:)
  real, save, allocatable :: cnuc_cp_hom_s (:)
  real, save, allocatable :: cnuc_dp_hom_s (:)

  real, save, allocatable :: rnuc_vp_haze_s  (:)
  real, save, allocatable :: rnuc_vp_immers_s(:)
  real, save, allocatable :: cnuc_vp_haze_s  (:)
  real, save, allocatable :: cnuc_vp_immers_s(:)

  real, save, allocatable :: rpsxfer_s(:)
  real, save, allocatable :: epsxfer_s(:)

  real, save, allocatable :: rx_s (:,:)
  real, save, allocatable :: cx_s (:,:)
  real, save, allocatable :: emb_s(:,:)
  real, save, allocatable :: dmb_s(:,:)
  real, save, allocatable :: qx_s (:,:)
  real, save, allocatable :: tx_s (:,:)
  real, save, allocatable :: txa_s(:,:)
  real, save, allocatable :: vap_s(:,:)

  real, save, allocatable :: &
  r1118_s(:,:),r8882_s(:,:),r1112_s(:,:),r1818_s(:,:),r1212_s(:,:), &
  r8282_s(:,:),r3335_s(:,:),r4445_s(:,:),r3435_s(:,:),r3445_s(:,:), &
  r3535_s(:,:),r3636_s(:,:),r3737_s(:,:),r4545_s(:,:),r4646_s(:,:), &
  r4747_s(:,:),r5656_s(:,:),r5757_s(:,:),r6767_s(:,:),r1413_s(:,:), &
  r1416_s(:,:),r1414_s(:,:),r1446_s(:,:),r1513_s(:,:),r1516_s(:,:), &
  r1515_s(:,:),r1556_s(:,:),r1613_s(:,:),r1616_s(:,:),r8483_s(:,:), &
  r8486_s(:,:),r8484_s(:,:),r8446_s(:,:),r8583_s(:,:),r8586_s(:,:), &
  r8585_s(:,:),r8556_s(:,:),r8683_s(:,:),r8686_s(:,:),r1713_s(:,:), &
  r1717_s(:,:),r8783_s(:,:),r8787_s(:,:),r2332_s(:,:),r2327_s(:,:), &
  r2323_s(:,:),r2337_s(:,:),r2442_s(:,:),r2427_s(:,:),r2424_s(:,:), &
  r2447_s(:,:),r2552_s(:,:),r2527_s(:,:),r2525_s(:,:),r2557_s(:,:), &
  r2662_s(:,:),r2627_s(:,:),r2626_s(:,:),r2667_s(:,:),r2772_s(:,:), &
  r2727_s(:,:),r0000_s(:,:)

  real, save, allocatable :: &
  e1111_s(:),e1118_s(:),e8888_s(:),e8882_s(:),e1112_s(:), &
  e1811_s(:),e1211_s(:),e8288_s(:),e2222_s(:),e5555_s(:), &
  e6666_s(:),e7777_s(:),e3333_s(:),e3335_s(:),e4444_s(:), &
  e4445_s(:),e3433_s(:),e3444_s(:),e3435_s(:),e3445_s(:), &
  e3533_s(:),e3633_s(:),e3733_s(:),e4544_s(:),e4644_s(:), &
  e4744_s(:),e5655_s(:),e5755_s(:),e6766_s(:),e1413_s(:), &
  e1411_s(:),e1446_s(:),e1513_s(:),e1511_s(:),e1556_s(:), &
  e1613_s(:),e1611_s(:),e8483_s(:),e8488_s(:),e8446_s(:), &
  e8583_s(:),e8588_s(:),e8556_s(:),e8683_s(:),e8688_s(:), &
  e1713_s(:),e1711_s(:),e8783_s(:),e8788_s(:),e2322_s(:), &
  e2327_s(:),e2333_s(:),e2337_s(:),e2422_s(:),e2427_s(:), &
  e2444_s(:),e2447_s(:),e2522_s(:),e2527_s(:),e2555_s(:), &
  e2557_s(:),e2622_s(:),e2627_s(:),e2666_s(:),e2667_s(:), &
  e2722_s(:),e2727_s(:),e2777_s(:),e0000_s(:)

  integer, save :: ncnt = 0, ntim = 0
  integer :: k

  real :: time

  time = real(time8) + dtlm(1) ! (micphys plot called just before time8 update)

  if (ntim == 0) then
     ntim = nint(real(timmax8) / dtlm(1))

     allocate (press0_s   (ntim))
     allocate (thil0_s    (ntim))
     allocate (theta0_s   (ntim))
     allocate (tairc_s    (ntim))
     allocate (rhoa_s     (ntim))
     allocate (rhoi_s     (ntim))
     allocate (rhow_s     (ntim))
     allocate (rhov_s     (ntim))
     allocate (rhovslair_s(ntim))
     allocate (rhovsiair_s(ntim))
     allocate (rhl_s      (ntim))
     allocate (rhi_s      (ntim))
     allocate (hcat_s     (ntim))

     allocate (rnuc_vc_s       (ntim))
     allocate (rnuc_vd_s       (ntim))
     allocate (rnuc_cp_hom_s   (ntim))
     allocate (rnuc_dp_hom_s   (ntim))
     allocate (rnuc_vp_haze_s  (ntim))
     allocate (rnuc_vp_immers_s(ntim))

     allocate (cnuc_vc_s       (ntim))
     allocate (cnuc_vd_s       (ntim))
     allocate (cnuc_cp_hom_s   (ntim))
     allocate (cnuc_dp_hom_s   (ntim))
     allocate (cnuc_vp_haze_s  (ntim))
     allocate (cnuc_vp_immers_s(ntim))

     allocate (rpsxfer_s(ntim))
     allocate (epsxfer_s(ntim))

     allocate (rx_s  (ntim,ncat))
     allocate (cx_s  (ntim,ncat))
     allocate (emb_s (ntim,ncat))
     allocate (dmb_s (ntim,ncat))
     allocate (qx_s  (ntim,ncat))
     allocate (tx_s  (ntim,ncat))
     allocate (txa_s (ntim,ncat))
     allocate (vap_s (ntim,ncat))

  allocate ( &
  r1118_s(ntim,2),r8882_s(ntim,2),r1112_s(ntim,2),r1818_s(ntim,2),r1212_s(ntim,2), &
  r8282_s(ntim,2),r3335_s(ntim,2),r4445_s(ntim,2),r3435_s(ntim,2),r3445_s(ntim,2), &
  r3535_s(ntim,2),r3636_s(ntim,2),r3737_s(ntim,2),r4545_s(ntim,2),r4646_s(ntim,2), &
  r4747_s(ntim,2),r5656_s(ntim,2),r5757_s(ntim,2),r6767_s(ntim,2),r1413_s(ntim,2), &
  r1416_s(ntim,2),r1414_s(ntim,2),r1446_s(ntim,2),r1513_s(ntim,2),r1516_s(ntim,2), &
  r1515_s(ntim,2),r1556_s(ntim,2),r1613_s(ntim,2),r1616_s(ntim,2),r8483_s(ntim,2), &
  r8486_s(ntim,2),r8484_s(ntim,2),r8446_s(ntim,2),r8583_s(ntim,2),r8586_s(ntim,2), &
  r8585_s(ntim,2),r8556_s(ntim,2),r8683_s(ntim,2),r8686_s(ntim,2),r1713_s(ntim,2), &
  r1717_s(ntim,2),r8783_s(ntim,2),r8787_s(ntim,2),r2332_s(ntim,2),r2327_s(ntim,2), &
  r2323_s(ntim,2),r2337_s(ntim,2),r2442_s(ntim,2),r2427_s(ntim,2),r2424_s(ntim,2), &
  r2447_s(ntim,2),r2552_s(ntim,2),r2527_s(ntim,2),r2525_s(ntim,2),r2557_s(ntim,2), &
  r2662_s(ntim,2),r2627_s(ntim,2),r2626_s(ntim,2),r2667_s(ntim,2),r2772_s(ntim,2), &
  r2727_s(ntim,2),r0000_s(ntim,2))

  allocate ( &
  e1111_s(ntim),e1118_s(ntim),e8888_s(ntim),e8882_s(ntim),e1112_s(ntim), &
  e1811_s(ntim),e1211_s(ntim),e8288_s(ntim),e2222_s(ntim),e5555_s(ntim), &
  e6666_s(ntim),e7777_s(ntim),e3333_s(ntim),e3335_s(ntim),e4444_s(ntim), &
  e4445_s(ntim),e3433_s(ntim),e3444_s(ntim),e3435_s(ntim),e3445_s(ntim), &
  e3533_s(ntim),e3633_s(ntim),e3733_s(ntim),e4544_s(ntim),e4644_s(ntim), &
  e4744_s(ntim),e5655_s(ntim),e5755_s(ntim),e6766_s(ntim),e1413_s(ntim), &
  e1411_s(ntim),e1446_s(ntim),e1513_s(ntim),e1511_s(ntim),e1556_s(ntim), &
  e1613_s(ntim),e1611_s(ntim),e8483_s(ntim),e8488_s(ntim),e8446_s(ntim), &
  e8583_s(ntim),e8588_s(ntim),e8556_s(ntim),e8683_s(ntim),e8688_s(ntim), &
  e1713_s(ntim),e1711_s(ntim),e8783_s(ntim),e8788_s(ntim),e2322_s(ntim), &
  e2327_s(ntim),e2333_s(ntim),e2337_s(ntim),e2422_s(ntim),e2427_s(ntim), &
  e2444_s(ntim),e2447_s(ntim),e2522_s(ntim),e2527_s(ntim),e2555_s(ntim), &
  e2557_s(ntim),e2622_s(ntim),e2627_s(ntim),e2666_s(ntim),e2667_s(ntim), &
  e2722_s(ntim),e2727_s(ntim),e2777_s(ntim),e0000_s(ntim))

  endif

  ncnt = ncnt + 1

  k = 2

! Store quantities to be plotted in time-series arrays

  press0_s   (ncnt) = press0(k)    * 1.e-2 ! converts to hPa
  thil0_s    (ncnt) = thil0(k)
  theta0_s   (ncnt) = theta0(k)
  tairc_s    (ncnt) = tairc(k)
  rhoa_s     (ncnt) = rhoa(k)
  rhow_s     (ncnt) = rhow(k)      * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rhov_s     (ncnt) = rhov(k)      * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rhovslair_s(ncnt) = rhovslair(k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rhovsiair_s(ncnt) = rhovsiair(k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rhl_s      (ncnt) = rhov(k) / rhovslair(k) * 1.e2 ! converts to %
  rhi_s      (ncnt) = rhov(k) / rhovsiair(k) * 1.e2

  if (tairc(k) > 0.) then
           rhi_s(ncnt) =       rhl_s(ncnt)
     rhovsiair_s(ncnt) = rhovslair_s(ncnt)
  endif

  rnuc_vc_s       (ncnt) = rnuc_vc       (k) * rhoi(k) * 1.e6 ! converts to mg/kg from kg/m^3
  rnuc_vd_s       (ncnt) = rnuc_vd       (k) * rhoi(k) * 1.e6 ! converts to mg/kg from kg/m^3
  rnuc_cp_hom_s   (ncnt) = rnuc_cp_hom   (k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rnuc_dp_hom_s   (ncnt) = rnuc_dp_hom   (k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rnuc_vp_haze_s  (ncnt) = rnuc_vp_haze  (k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  rnuc_vp_immers_s(ncnt) = rnuc_vp_immers(k) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3

  cnuc_vc_s       (ncnt) = cnuc_vc       (k) * rhoi(k) * 1.e-6 ! converts to #/mg from #/m^3
  cnuc_vd_s       (ncnt) = cnuc_vd       (k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  cnuc_cp_hom_s   (ncnt) = cnuc_cp_hom   (k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  cnuc_dp_hom_s   (ncnt) = cnuc_dp_hom   (k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  cnuc_vp_haze_s  (ncnt) = cnuc_vp_haze  (k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  cnuc_vp_immers_s(ncnt) = cnuc_vp_immers(k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3

  rpsxfer_s(ncnt) = rpsxfer(k) * rhoi(k) * 1.e3  ! converts to g/kg from kg/m^3
  epsxfer_s(ncnt) = epsxfer(k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3

  rx_s(ncnt,1) = rx(k,1) * rhoi(k) * 1.e3 ! converts to g/kg from kg/kg
  rx_s(ncnt,2) = rx(k,2) * rhoi(k) * 1.e3
  rx_s(ncnt,3) = rx(k,3) * rhoi(k) * 1.e3
  rx_s(ncnt,4) = rx(k,4) * rhoi(k) * 1.e3
  rx_s(ncnt,5) = rx(k,5) * rhoi(k) * 1.e3
  rx_s(ncnt,6) = rx(k,6) * rhoi(k) * 1.e3
  rx_s(ncnt,7) = rx(k,7) * rhoi(k) * 1.e3
  rx_s(ncnt,8) = rx(k,8) * rhoi(k) * 1.e3

  cx_s(ncnt,1) = cx(k,1) * rhoi(k) * 1.e-6 ! converts to #/mg from #/m^3
  cx_s(ncnt,2) = cx(k,2) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  cx_s(ncnt,3) = cx(k,3) * rhoi(k) * 1.e-3
  cx_s(ncnt,4) = cx(k,4) * rhoi(k) * 1.e-3
  cx_s(ncnt,5) = cx(k,5) * rhoi(k) * 1.e-3
  cx_s(ncnt,6) = cx(k,6) * rhoi(k) * 1.e-3
  cx_s(ncnt,7) = cx(k,7) * rhoi(k) * 1.e-3
  cx_s(ncnt,8) = cx(k,8) * rhoi(k) * 1.e-3

  qx_s(ncnt,1) = qx(k,1) * 1.e-3 ! converts to J/g from J/kg
  qx_s(ncnt,2) = qx(k,2) * 1.e-3
  qx_s(ncnt,3) = qx(k,3) * 1.e-3
  qx_s(ncnt,4) = qx(k,4) * 1.e-3
  qx_s(ncnt,5) = qx(k,5) * 1.e-3
  qx_s(ncnt,6) = qx(k,6) * 1.e-3
  qx_s(ncnt,7) = qx(k,7) * 1.e-3
  qx_s(ncnt,8) = qx(k,8) * 1.e-3

  tx_s(ncnt,1) = tx(k,1)
  tx_s(ncnt,2) = tx(k,2)
  tx_s(ncnt,3) = tx(k,3)
  tx_s(ncnt,4) = tx(k,4)
  tx_s(ncnt,5) = tx(k,5)
  tx_s(ncnt,6) = tx(k,6)
  tx_s(ncnt,7) = tx(k,7)
  tx_s(ncnt,8) = tx(k,8)

  emb_s(ncnt,1) = emb(k,1) * 1.e12 ! converts to  nano g from kg
  emb_s(ncnt,2) = emb(k,2) * 1.e6  ! converts to      mg from kg
  emb_s(ncnt,3) = emb(k,3) * 1.e9  ! converts to micro g from kg
  emb_s(ncnt,4) = emb(k,4) * 1.e9  ! converts to micro g from kg
  emb_s(ncnt,5) = emb(k,5) * 1.e9  ! converts to micro g from kg
  emb_s(ncnt,6) = emb(k,6) * 1.e6  ! converts to      mg from kg
  emb_s(ncnt,7) = emb(k,7) * 1.e3  ! converts to       g from kg
  emb_s(ncnt,8) = emb(k,8) * 1.e9  ! converts to micro g from kg

  dmb_s(ncnt,1) = (emb(k,1) * cfmasi(1)) ** pwmasi(1) * 1.e6 ! converts to microns from m
  dmb_s(ncnt,2) = (emb(k,2) * cfmasi(2)) ** pwmasi(2) * 1.e3 ! converts to mm from m
  dmb_s(ncnt,3) = (emb(k,3) * cfmasi(jhcat(k,3))) ** pwmasi(jhcat(k,3)) * 1.e3 ! to mm from m
  dmb_s(ncnt,4) = (emb(k,4) * cfmasi(jhcat(k,4))) ** pwmasi(jhcat(k,4)) * 1.e3 ! to mm from m
  dmb_s(ncnt,5) = (emb(k,5) * cfmasi(5)) ** pwmasi(5) * 1.e3 ! converts to mm from m
  dmb_s(ncnt,6) = (emb(k,6) * cfmasi(6)) ** pwmasi(6) * 1.e3 ! converts to mm from m
  dmb_s(ncnt,7) = (emb(k,7) * cfmasi(7)) ** pwmasi(7) * 1.e3 ! converts to mm from m
  dmb_s(ncnt,8) = (emb(k,8) * cfmasi(8)) ** pwmasi(8) * 1.e6 ! converts to microns from m

  vap_s(ncnt,1) = vap(k,1) * rhoi(k) * 1.e3 ! converts to g/kg from kg/m^3
  vap_s(ncnt,2) = vap(k,2) * rhoi(k) * 1.e3
  vap_s(ncnt,3) = vap(k,3) * rhoi(k) * 1.e3
  vap_s(ncnt,4) = vap(k,4) * rhoi(k) * 1.e3
  vap_s(ncnt,5) = vap(k,5) * rhoi(k) * 1.e3
  vap_s(ncnt,6) = vap(k,6) * rhoi(k) * 1.e3
  vap_s(ncnt,7) = vap(k,7) * rhoi(k) * 1.e3
  vap_s(ncnt,8) = vap(k,8) * rhoi(k) * 1.e3

  txa_s(ncnt,1) = tx(k,1) - tairc(k)
  txa_s(ncnt,2) = tx(k,2) - tairc(k)
  txa_s(ncnt,3) = tx(k,3) - tairc(k)
  txa_s(ncnt,4) = tx(k,4) - tairc(k)
  txa_s(ncnt,5) = tx(k,5) - tairc(k)
  txa_s(ncnt,6) = tx(k,6) - tairc(k)
  txa_s(ncnt,7) = tx(k,7) - tairc(k)
  txa_s(ncnt,8) = tx(k,8) - tairc(k)

  if (rx(k,1) < rxmin(1)) txa_s(ncnt,1) = 0.
  if (rx(k,2) < rxmin(2)) txa_s(ncnt,2) = 0.
  if (rx(k,3) < rxmin(3)) txa_s(ncnt,3) = 0.
  if (rx(k,4) < rxmin(4)) txa_s(ncnt,4) = 0.
  if (rx(k,5) < rxmin(5)) txa_s(ncnt,5) = 0.
  if (rx(k,6) < rxmin(6)) txa_s(ncnt,6) = 0.
  if (rx(k,7) < rxmin(7)) txa_s(ncnt,7) = 0.
  if (rx(k,8) < rxmin(8)) txa_s(ncnt,8) = 0.

  if (jhcat(k,3) == 3) then
     hcat_s(ncnt) = 1
  else
     hcat_s(ncnt) = jhcat(k,3) - 7
  endif

  r1118_s(ncnt,:) = r1118(k,:) * rhoi(k) * 1.e3  ! converts to g/kg from kg/m^3
  r8882_s(ncnt,:) = r8882(k,:) * rhoi(k) * 1.e3
  r1112_s(ncnt,:) = r1112(k,:) * rhoi(k) * 1.e3
  r1818_s(ncnt,:) = r1818(k,:) * rhoi(k) * 1.e3
  r1212_s(ncnt,:) = r1212(k,:) * rhoi(k) * 1.e3
  r8282_s(ncnt,:) = r8282(k,:) * rhoi(k) * 1.e3
  r3335_s(ncnt,:) = r3335(k,:) * rhoi(k) * 1.e3
  r4445_s(ncnt,:) = r4445(k,:) * rhoi(k) * 1.e3
  r3435_s(ncnt,:) = r3435(k,:) * rhoi(k) * 1.e3
  r3445_s(ncnt,:) = r3445(k,:) * rhoi(k) * 1.e3
  r3535_s(ncnt,:) = r3535(k,:) * rhoi(k) * 1.e3
  r3636_s(ncnt,:) = r3636(k,:) * rhoi(k) * 1.e3
  r3737_s(ncnt,:) = r3737(k,:) * rhoi(k) * 1.e3
  r4545_s(ncnt,:) = r4545(k,:) * rhoi(k) * 1.e3
  r4646_s(ncnt,:) = r4646(k,:) * rhoi(k) * 1.e3
  r4747_s(ncnt,:) = r4747(k,:) * rhoi(k) * 1.e3
  r5656_s(ncnt,:) = r5656(k,:) * rhoi(k) * 1.e3
  r5757_s(ncnt,:) = r5757(k,:) * rhoi(k) * 1.e3
  r6767_s(ncnt,:) = r6767(k,:) * rhoi(k) * 1.e3
  r1413_s(ncnt,:) = r1413(k,:) * rhoi(k) * 1.e3
  r1416_s(ncnt,:) = r1416(k,:) * rhoi(k) * 1.e3
  r1414_s(ncnt,:) = r1414(k,:) * rhoi(k) * 1.e3
  r1446_s(ncnt,:) = r1446(k,:) * rhoi(k) * 1.e3
  r1513_s(ncnt,:) = r1513(k,:) * rhoi(k) * 1.e3
  r1516_s(ncnt,:) = r1516(k,:) * rhoi(k) * 1.e3
  r1515_s(ncnt,:) = r1515(k,:) * rhoi(k) * 1.e3
  r1556_s(ncnt,:) = r1556(k,:) * rhoi(k) * 1.e3
  r1613_s(ncnt,:) = r1613(k,:) * rhoi(k) * 1.e3
  r1616_s(ncnt,:) = r1616(k,:) * rhoi(k) * 1.e3
  r8483_s(ncnt,:) = r8483(k,:) * rhoi(k) * 1.e3
  r8486_s(ncnt,:) = r8486(k,:) * rhoi(k) * 1.e3
  r8484_s(ncnt,:) = r8484(k,:) * rhoi(k) * 1.e3
  r8446_s(ncnt,:) = r8446(k,:) * rhoi(k) * 1.e3
  r8583_s(ncnt,:) = r8583(k,:) * rhoi(k) * 1.e3
  r8586_s(ncnt,:) = r8586(k,:) * rhoi(k) * 1.e3
  r8585_s(ncnt,:) = r8585(k,:) * rhoi(k) * 1.e3
  r8556_s(ncnt,:) = r8556(k,:) * rhoi(k) * 1.e3
  r8683_s(ncnt,:) = r8683(k,:) * rhoi(k) * 1.e3
  r8686_s(ncnt,:) = r8686(k,:) * rhoi(k) * 1.e3
  r1713_s(ncnt,:) = r1713(k,:) * rhoi(k) * 1.e3
  r1717_s(ncnt,:) = r1717(k,:) * rhoi(k) * 1.e3
  r8783_s(ncnt,:) = r8783(k,:) * rhoi(k) * 1.e3
  r8787_s(ncnt,:) = r8787(k,:) * rhoi(k) * 1.e3
  r2332_s(ncnt,:) = r2332(k,:) * rhoi(k) * 1.e3
  r2327_s(ncnt,:) = r2327(k,:) * rhoi(k) * 1.e3
  r2323_s(ncnt,:) = r2323(k,:) * rhoi(k) * 1.e3
  r2337_s(ncnt,:) = r2337(k,:) * rhoi(k) * 1.e3
  r2442_s(ncnt,:) = r2442(k,:) * rhoi(k) * 1.e3
  r2427_s(ncnt,:) = r2427(k,:) * rhoi(k) * 1.e3
  r2424_s(ncnt,:) = r2424(k,:) * rhoi(k) * 1.e3
  r2447_s(ncnt,:) = r2447(k,:) * rhoi(k) * 1.e3
  r2552_s(ncnt,:) = r2552(k,:) * rhoi(k) * 1.e3
  r2527_s(ncnt,:) = r2527(k,:) * rhoi(k) * 1.e3
  r2525_s(ncnt,:) = r2525(k,:) * rhoi(k) * 1.e3
  r2557_s(ncnt,:) = r2557(k,:) * rhoi(k) * 1.e3
  r2662_s(ncnt,:) = r2662(k,:) * rhoi(k) * 1.e3
  r2627_s(ncnt,:) = r2627(k,:) * rhoi(k) * 1.e3
  r2626_s(ncnt,:) = r2626(k,:) * rhoi(k) * 1.e3
  r2667_s(ncnt,:) = r2667(k,:) * rhoi(k) * 1.e3
  r2772_s(ncnt,:) = r2772(k,:) * rhoi(k) * 1.e3
  r2727_s(ncnt,:) = r2727(k,:) * rhoi(k) * 1.e3
  r0000_s(ncnt,:) = r0000(k,:) * rhoi(k) * 1.e3

  e1111_s(ncnt) = e1111(k) * rhoi(k) * 1.e-3 ! converts to #/g from #/m^3
  e1118_s(ncnt) = e1118(k) * rhoi(k) * 1.e-3
  e8888_s(ncnt) = e8888(k) * rhoi(k) * 1.e-3
  e8882_s(ncnt) = e8882(k) * rhoi(k) * 1.e-3
  e1112_s(ncnt) = e1112(k) * rhoi(k) * 1.e-3
  e1811_s(ncnt) = e1811(k) * rhoi(k) * 1.e-3
  e1211_s(ncnt) = e1211(k) * rhoi(k) * 1.e-3
  e8288_s(ncnt) = e8288(k) * rhoi(k) * 1.e-3
  e2222_s(ncnt) = e2222(k) * rhoi(k) * 1.e-3
  e5555_s(ncnt) = e5555(k) * rhoi(k) * 1.e-3
  e6666_s(ncnt) = e6666(k) * rhoi(k) * 1.e-3
  e7777_s(ncnt) = e7777(k) * rhoi(k) * 1.e-3
  e3333_s(ncnt) = e3333(k) * rhoi(k) * 1.e-3
  e3335_s(ncnt) = e3335(k) * rhoi(k) * 1.e-3
  e4444_s(ncnt) = e4444(k) * rhoi(k) * 1.e-3
  e4445_s(ncnt) = e4445(k) * rhoi(k) * 1.e-3
  e3433_s(ncnt) = e3433(k) * rhoi(k) * 1.e-3
  e3444_s(ncnt) = e3444(k) * rhoi(k) * 1.e-3
  e3435_s(ncnt) = e3435(k) * rhoi(k) * 1.e-3
  e3445_s(ncnt) = e3445(k) * rhoi(k) * 1.e-3
  e3533_s(ncnt) = e3533(k) * rhoi(k) * 1.e-3
  e3633_s(ncnt) = e3633(k) * rhoi(k) * 1.e-3
  e3733_s(ncnt) = e3733(k) * rhoi(k) * 1.e-3
  e4544_s(ncnt) = e4544(k) * rhoi(k) * 1.e-3
  e4644_s(ncnt) = e4644(k) * rhoi(k) * 1.e-3
  e4744_s(ncnt) = e4744(k) * rhoi(k) * 1.e-3
  e5655_s(ncnt) = e5655(k) * rhoi(k) * 1.e-3
  e5755_s(ncnt) = e5755(k) * rhoi(k) * 1.e-3
  e6766_s(ncnt) = e6766(k) * rhoi(k) * 1.e-3
  e1413_s(ncnt) = e1413(k) * rhoi(k) * 1.e-3
  e1411_s(ncnt) = e1411(k) * rhoi(k) * 1.e-3
  e1446_s(ncnt) = e1446(k) * rhoi(k) * 1.e-3
  e1513_s(ncnt) = e1513(k) * rhoi(k) * 1.e-3
  e1511_s(ncnt) = e1511(k) * rhoi(k) * 1.e-3
  e1556_s(ncnt) = e1556(k) * rhoi(k) * 1.e-3
  e1613_s(ncnt) = e1613(k) * rhoi(k) * 1.e-3
  e1611_s(ncnt) = e1611(k) * rhoi(k) * 1.e-3
  e8483_s(ncnt) = e8483(k) * rhoi(k) * 1.e-3
  e8488_s(ncnt) = e8488(k) * rhoi(k) * 1.e-3
  e8446_s(ncnt) = e8446(k) * rhoi(k) * 1.e-3
  e8583_s(ncnt) = e8583(k) * rhoi(k) * 1.e-3
  e8588_s(ncnt) = e8588(k) * rhoi(k) * 1.e-3
  e8556_s(ncnt) = e8556(k) * rhoi(k) * 1.e-3
  e8683_s(ncnt) = e8683(k) * rhoi(k) * 1.e-3
  e8688_s(ncnt) = e8688(k) * rhoi(k) * 1.e-3
  e1713_s(ncnt) = e1713(k) * rhoi(k) * 1.e-3
  e1711_s(ncnt) = e1711(k) * rhoi(k) * 1.e-3
  e8783_s(ncnt) = e8783(k) * rhoi(k) * 1.e-3
  e8788_s(ncnt) = e8788(k) * rhoi(k) * 1.e-3
  e2322_s(ncnt) = e2322(k) * rhoi(k) * 1.e-3
  e2327_s(ncnt) = e2327(k) * rhoi(k) * 1.e-3
  e2333_s(ncnt) = e2333(k) * rhoi(k) * 1.e-3
  e2337_s(ncnt) = e2337(k) * rhoi(k) * 1.e-3
  e2422_s(ncnt) = e2422(k) * rhoi(k) * 1.e-3
  e2427_s(ncnt) = e2427(k) * rhoi(k) * 1.e-3
  e2444_s(ncnt) = e2444(k) * rhoi(k) * 1.e-3
  e2447_s(ncnt) = e2447(k) * rhoi(k) * 1.e-3
  e2522_s(ncnt) = e2522(k) * rhoi(k) * 1.e-3
  e2527_s(ncnt) = e2527(k) * rhoi(k) * 1.e-3
  e2555_s(ncnt) = e2555(k) * rhoi(k) * 1.e-3
  e2557_s(ncnt) = e2557(k) * rhoi(k) * 1.e-3
  e2622_s(ncnt) = e2622(k) * rhoi(k) * 1.e-3
  e2627_s(ncnt) = e2627(k) * rhoi(k) * 1.e-3
  e2666_s(ncnt) = e2666(k) * rhoi(k) * 1.e-3
  e2667_s(ncnt) = e2667(k) * rhoi(k) * 1.e-3
  e2722_s(ncnt) = e2722(k) * rhoi(k) * 1.e-3
  e2727_s(ncnt) = e2727(k) * rhoi(k) * 1.e-3
  e2777_s(ncnt) = e2777(k) * rhoi(k) * 1.e-3
  e0000_s(ncnt) = e0000(k) * rhoi(k) * 1.e-3

  if (ncnt == ntim) then
     call o_reopnwk()

     call o_sflush()
     call o_gsplci(10)
     call o_gstxci(10)

     call o_set(0.,1.,0.,1.,0.,1.,0.,1.,1)

! WALKO ET AL. 2000 FIRST PLOT

     call plotback()

     call plotpar(1,ntim,time,'CPV MIX RAT','(g/kg)',0.,10.,1.)

     call bdline (1,ntim,  rx_s(1:ntim,1))
     call bdline (2,ntim,  rx_s(1:ntim,3))
     call bdline (3,ntim,rhov_s(1:ntim  ))

     call plotpar(2,ntim,time,'CP VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim, vap_s(1:ntim,1))
     call bdline (2,ntim, vap_s(1:ntim,3))

     call plotpar(3,ntim,time,'CP INT ENERGY','(J/g)',-200.,500.,100.)

     call bdline (1,ntim,  qx_s(1:ntim,1))
     call bdline (2,ntim,  qx_s(1:ntim,3))

     call plotpar(4,ntim,time,'CP TEMP DIFF','(K)',-.1,1.,.1)

     call bdline (1,ntim, txa_s(1:ntim,1))
     call bdline (2,ntim, txa_s(1:ntim,3))

     call plotpar(5,ntim,time,'LIQ & ICE R.H.','(%)',50.,150.,10.)

     call bdline (1,ntim, rhl_s(1:ntim  ))
     call bdline (2,ntim, rhi_s(1:ntim  ))

     call o_frame()

! WALKO ET AL. 2000 SECOND PLOT

     call plotback()

     call plotpar(1,ntim,time,'RHV MIX RAT','(g/kg)',0.,10.,1.)

     call bdline (1,ntim,  rx_s(1:ntim,2))
     call bdline (2,ntim,  rx_s(1:ntim,7))
     call bdline (3,ntim,rhov_s(1:ntim  ))

     call plotpar(2,ntim,time,'RH VAP FLUX','(g/kg)',-.02,.01,.004)

     call bdline (1,ntim, vap_s(1:ntim,2))
     call bdline (2,ntim, vap_s(1:ntim,7))

     call plotpar(3,ntim,time,'RH INT ENERGY','(J/g)',0.,500.,100.)

     call bdline (1,ntim,  qx_s(1:ntim,2))
     call bdline (2,ntim,  qx_s(1:ntim,7))

     call plotpar(4,ntim,time,'R & AIR TEMP','(:0557:C)',8.,20.,2.)

     call bdline (1,ntim,   tx_s(1:ntim,2))
     call bdline (2,ntim,tairc_s(1:ntim))

     call plotpar(5,ntim,time,'LIQ R.H.','(%)',30.,110.,10.)

     call bdline (1,ntim, rhl_s(1:ntim  ))

     call o_frame()

! CLOUD PLOT

     call plotback()

     call plotpar(1,ntim,time,'C MIX RATIO','(g/kg)',0.,2.,.2)
     call bdline (1,ntim,rx_s(1:ntim,1))

     call plotpar(2,ntim,time,'C NUM CONC','(#/mg)',0.,2000.,200.)
     call bdline (1,ntim,cx_s(1:ntim,1))

     call plotpar(3,ntim,time,'C EMB','(ng)',0.,10.,1.)
     call bdline (1,ntim,emb_s(1:ntim,1))

     call plotpar(4,ntim,time,'C DMB','(:0714:m)',0.,50.,10.)
     call bdline (1,ntim,dmb_s(1:ntim,1))

     call plotpar(5,ntim,time,'C VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,1))

     call o_frame()

! DRIZZLE PLOT

     call plotback()

     call plotpar(1,ntim,time,'D MIX RATIO','(g/kg)',0.,.2,.02)
     call bdline (1,ntim,rx_s(1:ntim,8))

     call plotpar(2,ntim,time,'DR NUM CONC','(#/g)',0.,500.,50.)
     call bdline (1,ntim,cx_s(1:ntim,8))

     call plotpar(3,ntim,time,'D EMB','(:0714:g)',0.,1.,.1) ! :0714: is Greek "mu"
     call bdline (1,ntim, emb_s(1:ntim,8))

     call plotpar(4,ntim,time,'D DMB','(:0714:m)',0.,200.,20.)
     call bdline (1,ntim,dmb_s(1:ntim,8))

     call plotpar(5,ntim,time,'D VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,8))

     call o_frame()

! RAIN PLOT

     call plotback()

     call plotpar(1,ntim,time,'R MIX RATIO','(g/kg)',0.,1.,.1)
     call bdline (1,ntim,rx_s(1:ntim,2))

     call plotpar(2,ntim,time,'R NUM CONC','(#/g)',0.,500.,50.)
     call bdline (1,ntim,cx_s(1:ntim,2))

     call plotpar(3,ntim,time,'R EMB','(mg)',0.,10.,1.)
     call bdline (1,ntim,emb_s(1:ntim,2))

     call plotpar(4,ntim,time,'R DMB','(mm)',0.,10.,1.)
     call bdline (1,ntim,dmb_s(1:ntim,2))

     call plotpar(5,ntim,time,'R VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,2))

     call o_frame()

! PRISTINE ICE PLOT

     call plotback()

     call plotpar(1,ntim,time,'P MIX RATIO','(g/kg)',0.,.1,.01)
     call bdline (1,ntim,rx_s(1:ntim,3))

     call plotpar(2,ntim,time,'P NUM CONC','(#/g)',0.,500.,50.)
     call bdline (1,ntim,cx_s(1:ntim,3))

     call plotpar(3,ntim,time,'P EMB','(:0714:g)',0.,1.,.1)
     call bdline (1,ntim,emb_s(1:ntim,3))

     call plotpar(4,ntim,time,'P DMB','(mm)',0.,1.,.1)
     call bdline (1,ntim,dmb_s(1:ntim,3))

     call plotpar(5,ntim,time,'P VAP & XFER FLUX','(g/kg)',-.01,.01,.002)
     call bdline (1,ntim,vap_s(1:ntim,3))
     call bdline (2,ntim,rpsxfer_s(1:ntim))

     call o_frame()

! SNOW PLOT

     call plotback()

     call plotpar(1,ntim,time,'S MIX RATIO','(g/kg)',0.,.5,.1)
     call bdline (1,ntim,rx_s(1:ntim,4))

     call plotpar(2,ntim,time,'S NUM CONC','(#/g)',0.,500.,100.)
     call bdline (1,ntim,cx_s(1:ntim,4))

     call plotpar(3,ntim,time,'S EMB','(:0714:g)',0.,10.,1.)
     call bdline (1,ntim, emb_s(1:ntim,4))

     call plotpar(4,ntim,time,'S DMB','(mm)',0.,5.,1.)
     call bdline (1,ntim,dmb_s(1:ntim,4))

     call plotpar(5,ntim,time,'S VAP & XFER FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,4))
     call bdline (2,ntim,rpsxfer_s(1:ntim))

     call o_frame()

! AGGREGATES PLOT

     call plotback()

     call plotpar(1,ntim,time,'A MIX RATIO','(g/kg)',0.,10.,1.)
     call bdline (1,ntim,rx_s(1:ntim,5))

     call plotpar(2,ntim,time,'A NUM CONC','(#/g)',0.,500.,100.)
     call bdline (1,ntim,cx_s(1:ntim,5))

     call plotpar(3,ntim,time,'A EMB','(:0714:g)',0.,1000.,100.)
     call bdline (1,ntim, emb_s(1:ntim,5))

     call plotpar(4,ntim,time,'A DMB','(mm)',0.,10.,1.)
     call bdline (1,ntim,dmb_s(1:ntim,5))

     call plotpar(5,ntim,time,'A VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,5))

     call o_frame()

! GRAUPEL PLOT

     call plotback()

     call plotpar(1,ntim,time,'G MIX RATIO','(g/kg)',0.,1.,.1)
     call bdline (1,ntim,rx_s(1:ntim,6))

     call plotpar(2,ntim,time,'G NUM CONC','(#/g)',0.,10.,1.)
     call bdline (1,ntim,cx_s(1:ntim,6))

     call plotpar(3,ntim,time,'G EMB','(mg)',0.,10.,1.)
     call bdline (1,ntim, emb_s(1:ntim,6))

     call plotpar(4,ntim,time,'G DMB','(mm)',0.,10.,1.)
     call bdline (1,ntim,dmb_s(1:ntim,6))

     call plotpar(5,ntim,time,'G VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,6))

     call o_frame()

! HAIL PLOT

     call plotback()

     call plotpar(1,ntim,time,'H MIX RATIO','(g/kg)',0.,1.,.1)
     call bdline (1,ntim,rx_s(1:ntim,7))

     call plotpar(2,ntim,time,'H NUM CONC','(#/g)',0.,10.,1.)
     call bdline (1,ntim,cx_s(1:ntim,7))

     call plotpar(3,ntim,time,'H EMB','(g)',0.,10.,1.)
     call bdline (1,ntim, emb_s(1:ntim,7))

     call plotpar(4,ntim,time,'H DMB','(mm)',0.,50.,5.)
     call bdline (1,ntim,dmb_s(1:ntim,7))

     call plotpar(5,ntim,time,'H VAP FLUX','(g/kg)',-.05,.05,.01)
     call bdline (1,ntim,vap_s(1:ntim,7))

     call o_frame()

! PLOT 1

     call plotback()

     call plotpar(1,ntim,time,'C MIX RATIO','(g/kg)',0.,2.,.2)

     call bdline (1,ntim,rx_s(1:ntim,1))

     call plotpar(2,ntim,time,'D MIX RATIO','(g/kg)',0.,.2,.02)

     call bdline (1,ntim,rx_s(1:ntim,8))

     call plotpar(3,ntim,time,'R MIX RATIO','(g/kg)',0.,1.,.1)

     call bdline (1,ntim,rx_s(1:ntim,2))

     call plotpar(4,ntim,time,'PSA MIX RATIO','(g/kg)',0.,1.,.1)

     call bdline (1,ntim,rx_s(1:ntim,3))
     call bdline (2,ntim,rx_s(1:ntim,4))
     call bdline (3,ntim,rx_s(1:ntim,5))

     call plotpar(5,ntim,time,'GH MIX RATIO','(g/kg)',0.,1.,.1)

     call bdline (1,ntim,rx_s(1:ntim,6))
     call bdline (2,ntim,rx_s(1:ntim,7))

     call o_frame()

! PLOT 2

     call plotback()

     call plotpar(1,ntim,time,'C NUM CONC','(#/mg)',0.,2000.,200.)

     call bdline (1,ntim,cx_s(1:ntim,1))

     call plotpar(2,ntim,time,'DR NUM CONC','(#/g)',0.,500.,50.)

     call bdline (1,ntim,cx_s(1:ntim,8))
     call bdline (2,ntim,cx_s(1:ntim,2))

     call plotpar(3,ntim,time,'P NUM CONC','(#/g)',0.,5000.,500.)

     call bdline (1,ntim,cx_s(1:ntim,3))

     call plotpar(4,ntim,time,'SA NUM CONC','(#/g)',0.,1000.,100.)

     call bdline (1,ntim,cx_s(1:ntim,4))
     call bdline (2,ntim,cx_s(1:ntim,5))

     call plotpar(5,ntim,time,'GH NUM CONC','(#/g)',0.,10.,1.)

     call bdline (1,ntim,cx_s(1:ntim,6))
     call bdline (2,ntim,cx_s(1:ntim,7))

     call o_frame()

! PLOT 3

     call plotback()

     call plotpar(1,ntim,time,'C EMB','(ng)',0.,10.,1.)

     call bdline (1,ntim, emb_s(1:ntim,1))

     call plotpar(2,ntim,time,'D EMB','(:0714:g)',0.,1.,.1) ! :0714: is Greek "mu"

     call bdline (1,ntim, emb_s(1:ntim,8))

     call plotpar(3,ntim,time,'RG EMB','(mg)',0.,10.,1.)

     call bdline (1,ntim, emb_s(1:ntim,2))
     call bdline (2,ntim, emb_s(1:ntim,6))

     call plotpar(4,ntim,time,'PSA EMB','(mg)',0.,1.,.1)

     call bdline (1,ntim, emb_s(1:ntim,3))
     call bdline (2,ntim, emb_s(1:ntim,4))
     call bdline (3,ntim, emb_s(1:ntim,5))

     call plotpar(5,ntim,time,'H EMB','(g)',0.,10.,1.)

     call bdline (1,ntim, emb_s(1:ntim,7))

     call o_frame()

! PLOT 4

     call plotback()

     call plotpar(1,ntim,time,'C DMB','(:0714:m)',0.,50.,10.)

     call bdline (1,ntim,dmb_s(1:ntim,1))

     call plotpar(2,ntim,time,'D DMB','(:0714:m)',0.,200.,20.)

     call bdline (1,ntim,dmb_s(1:ntim,8))

     call plotpar(3,ntim,time,'RG DMB','(mm)',0.,10.,1.)

     call bdline (1,ntim,dmb_s(1:ntim,2))
     call bdline (2,ntim,dmb_s(1:ntim,6))

     call plotpar(4,ntim,time,'PSA DMB','(mm)',0.,10.,1.)

     call bdline (1,ntim,dmb_s(1:ntim,3))
     call bdline (2,ntim,dmb_s(1:ntim,4))
     call bdline (3,ntim,dmb_s(1:ntim,5))

     call plotpar(5,ntim,time,'H DMB','(mm)',0.,50.,5.)

     call bdline (1,ntim,dmb_s(1:ntim,7))

     call o_frame()

! PLOT 5

     call plotback()

     call plotpar(1,ntim,time,'C VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim,vap_s(1:ntim,1))

     call plotpar(2,ntim,time,'D VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim,vap_s(1:ntim,8))

     call plotpar(3,ntim,time,'R VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim,vap_s(1:ntim,2))

     call plotpar(4,ntim,time,'PSA VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim,vap_s(1:ntim,3))
     call bdline (2,ntim,vap_s(1:ntim,4))
     call bdline (3,ntim,vap_s(1:ntim,5))

     call plotpar(5,ntim,time,'GH VAP FLUX','(g/kg)',-.05,.05,.01)

     call bdline (1,ntim,vap_s(1:ntim,6))
     call bdline (2,ntim,vap_s(1:ntim,7))

     call o_frame()

! PLOT 6

     call plotback()

     call plotpar(1,ntim,time,'CDR INT ENERGY','(J/g)',-200.,500.,100.)

     call bdline (1,ntim,qx_s(1:ntim,1))
     call bdline (2,ntim,qx_s(1:ntim,8))
     call bdline (3,ntim,qx_s(1:ntim,2))

     call plotpar(2,ntim,time,'PSA INT ENERGY','(J/g)',-200.,500.,100.)

     call bdline (1,ntim,qx_s(1:ntim,3))
     call bdline (2,ntim,qx_s(1:ntim,4))
     call bdline (3,ntim,qx_s(1:ntim,5))

     call plotpar(3,ntim,time,'GH INT ENERGY','(J/g)',-200.,500.,100.)

     call bdline (1,ntim,qx_s(1:ntim,6))
     call bdline (2,ntim,qx_s(1:ntim,7))

     call o_frame()

! PLOT 7

     call plotback()

     call plotpar(1,ntim,time,'PRESSURE','(kPa)',0.,1100.,100.)

     call bdline (1,ntim,press0_s(1:ntim))

     call plotpar(2,ntim,time,'THIL & THETA','(K)',250.,350.,10.)

     call bdline (1,ntim, thil0_s(1:ntim))
     call bdline (2,ntim,theta0_s(1:ntim))

     call plotpar(3,ntim,time,'TAIRC','(:0557:C)',-60.,50.,10.) ! :0557: is 'degree' symbol

     call bdline (1,ntim,tairc_s(1:ntim))

     call plotpar(4,ntim,time,'RHOA','(kg/m^3)',0.,1.3,.1)

     call bdline (1,ntim,rhoa_s(1:ntim))

     call o_frame()

! PLOT 8

     call plotback()

     call plotpar(1,ntim,time,'RHOW & RHOV','(g/kg)',0.,20.,2.)

     call bdline (1,ntim,rhow_s(1:ntim))
     call bdline (2,ntim,rhov_s(1:ntim))

     call plotpar(2,ntim,time,'RSATL & RSATI','(g/kg)',0.,20.,2.)

     call bdline (1,ntim,rhovslair_s(1:ntim))
     call bdline (2,ntim,rhovsiair_s(1:ntim))

     call plotpar(3,ntim,time,'RHL & RHI','(%)',0.,150.,10.)

     call bdline (1,ntim,rhl_s(1:ntim))
     call bdline (2,ntim,rhi_s(1:ntim))

     call plotpar(4,ntim,time,'HCAT','(#)',0.,5.,1.)

     call bdline (1,ntim,hcat_s(1:ntim))

     call o_frame()

! PLOT 9

     call plotback()

     call plotpar(1,ntim,time,'CDR TEMP DIFF','(K)',-1.,1.,.2)

     call bdline (1,ntim,txa_s(1:ntim,1))
     call bdline (2,ntim,txa_s(1:ntim,8))
     call bdline (3,ntim,txa_s(1:ntim,2))

     call plotpar(2,ntim,time,'PSA TEMP DIFF','(K)',-1.,1.0,.2)

     call bdline (1,ntim,txa_s(1:ntim,3))
     call bdline (2,ntim,txa_s(1:ntim,4))
     call bdline (3,ntim,txa_s(1:ntim,5))

     call plotpar(3,ntim,time,'GH TEMP DIFF','(K)',-1.,1.0,.2)

     call bdline (1,ntim,txa_s(1:ntim,7))
     call bdline (2,ntim,txa_s(1:ntim,8))

     call o_frame()

! PLOT 10

     call plotback()

     call plotpar(1,ntim,time,'RNUC_VC','(mg/kg)',0.,5.,1.)

     call bdline (1,ntim,rnuc_vc_s(1:ntim))

     call plotpar(2,ntim,time,'RNUC_VD','(mg/kg)',0.,.05,.01)

     call bdline (1,ntim,rnuc_vd_s(1:ntim))

     call plotpar(3,ntim,time,'CNUC_VC','(#/mg)',0.,1000.,100.)

     call bdline (1,ntim,cnuc_vc_s(1:ntim))

     call plotpar(4,ntim,time,'CNUC_VD','(#/g)',0.,.2,.02)

     call bdline (1,ntim,cnuc_vd_s(1:ntim))

     call o_frame()

! PLOT 11

     call plotback()

     call plotpar(1,ntim,time,'RNUC_CP_HOM','(g/kg)',0.,1.,.2)

     call bdline (1,ntim,rnuc_cp_hom_s(1:ntim))

     call plotpar(2,ntim,time,'RNUC_DP_HOM','(g/kg)',0.,1.,.2)

     call bdline (1,ntim,rnuc_dp_hom_s(1:ntim))

     call plotpar(3,ntim,time,'RNUC_VP_immers','(g/kg)',0.,1.,.2)

     call bdline (1,ntim,rnuc_vp_immers_s(1:ntim))

     call plotpar(4,ntim,time,'RNUC_VP_HAZE','(g/kg)',0.,1.,.2)

     call bdline (1,ntim,rnuc_vp_haze_s(1:ntim))

     call o_frame()

! PLOT 12

     call plotback()

     call plotpar(1,ntim,time,'CNUC_CP_HOM','(#/g)',0.,1.,.2)

     call bdline (1,ntim,cnuc_cp_hom_s(1:ntim))

     call plotpar(2,ntim,time,'CNUC_DP_HOM','(#/g)',0.,1.,.2)

     call bdline (1,ntim,cnuc_dp_hom_s(1:ntim))

     call plotpar(3,ntim,time,'CNUC_VP_immers','(#/g)',0.,1.,.2)

     call bdline (1,ntim,cnuc_vp_immers_s(1:ntim))

     call plotpar(4,ntim,time,'CNUC_VP_HAZE','(#/g)',0.,1.,.2)

     call bdline (1,ntim,cnuc_vp_haze_s(1:ntim))

     call o_frame()

! R PLOT 1

     call plotback()

     call plotpar(1,ntim,time,'R1118','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1118_s(1:ntim,1))

     call plotpar(2,ntim,time,'R8882','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8882_s(1:ntim,1))

     call plotpar(3,ntim,time,'R1112','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1112_s(1:ntim,1))

     call plotpar(4,ntim,time,'R1818','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1818_s(1:ntim,1))

     call plotpar(5,ntim,time,'R1212','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1212_s(1:ntim,1))

     call o_frame()

! R PLOT 2

     call plotback()

     call plotpar(1,ntim,time,'R8282','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8282_s(1:ntim,1))

     call plotpar(2,ntim,time,'R3335','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3335_s(1:ntim,1))

     call plotpar(3,ntim,time,'R4445','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r4445_s(1:ntim,1))

     call plotpar(4,ntim,time,'R3435','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3435_s(1:ntim,1))

     call plotpar(5,ntim,time,'R3445','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3445_s(1:ntim,1))

     call o_frame()

! R PLOT 3

     call plotback()

     call plotpar(1,ntim,time,'R3535','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3535_s(1:ntim,1))

     call plotpar(2,ntim,time,'R3636','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3636_s(1:ntim,1))

     call plotpar(3,ntim,time,'R3737','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r3737_s(1:ntim,1))

     call plotpar(4,ntim,time,'R4545','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r4545_s(1:ntim,1))

     call plotpar(5,ntim,time,'R4646','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r4646_s(1:ntim,1))

     call o_frame()

! R PLOT 4

     call plotback()

     call plotpar(1,ntim,time,'R4747','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r4747_s(1:ntim,1))

     call plotpar(2,ntim,time,'R5656','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r5656_s(1:ntim,1))

     call plotpar(3,ntim,time,'R5757','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r5757_s(1:ntim,1))

     call plotpar(4,ntim,time,'R6767','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r6767_s(1:ntim,1))

     call plotpar(5,ntim,time,'R1413','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1413_s(1:ntim,1))

     call o_frame()

! R PLOT 5

     call plotback()

     call plotpar(1,ntim,time,'R1416','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1416_s(1:ntim,1))

     call plotpar(2,ntim,time,'R1414','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1414_s(1:ntim,1))

     call plotpar(3,ntim,time,'R1446','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1446_s(1:ntim,1))

     call plotpar(4,ntim,time,'R1513','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1513_s(1:ntim,1))

     call plotpar(5,ntim,time,'R1516','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1516_s(1:ntim,1))

     call o_frame()

! R PLOT 6

     call plotback()

     call plotpar(1,ntim,time,'R1515','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1515_s(1:ntim,1))

     call plotpar(2,ntim,time,'R1556','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1556_s(1:ntim,1))

     call plotpar(3,ntim,time,'R1613','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1613_s(1:ntim,1))

     call plotpar(4,ntim,time,'R1616','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1616_s(1:ntim,1))

     call plotpar(5,ntim,time,'R8483','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8483_s(1:ntim,1))

     call o_frame()

! R PLOT 7

     call plotback()

     call plotpar(1,ntim,time,'R8486','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8486_s(1:ntim,1))

     call plotpar(2,ntim,time,'R8484','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8484_s(1:ntim,1))

     call plotpar(3,ntim,time,'R8446','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8446_s(1:ntim,1))

     call plotpar(4,ntim,time,'R8583','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8583_s(1:ntim,1))

     call plotpar(5,ntim,time,'R8586','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8586_s(1:ntim,1))

     call o_frame()

! R PLOT 8

     call plotback()

     call plotpar(1,ntim,time,'R8585','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8585_s(1:ntim,1))

     call plotpar(2,ntim,time,'R8556','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8556_s(1:ntim,1))

     call plotpar(3,ntim,time,'R8683','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8683_s(1:ntim,1))

     call plotpar(4,ntim,time,'R8686','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8686_s(1:ntim,1))

     call plotpar(5,ntim,time,'R1713','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1713_s(1:ntim,1))

     call o_frame()

! R PLOT 9

     call plotback()

     call plotpar(1,ntim,time,'R1717','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r1717_s(1:ntim,1))

     call plotpar(2,ntim,time,'R8783','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8783_s(1:ntim,1))

     call plotpar(3,ntim,time,'R8787','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r8787_s(1:ntim,1))

     call plotpar(4,ntim,time,'R2332','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2332_s(1:ntim,1))

     call plotpar(5,ntim,time,'R2327','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2327_s(1:ntim,1))

     call o_frame()

! R PLOT 10

     call plotback()

     call plotpar(1,ntim,time,'R2323','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2323_s(1:ntim,1))

     call plotpar(2,ntim,time,'R2337','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2337_s(1:ntim,1))

     call plotpar(3,ntim,time,'R2442','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2442_s(1:ntim,1))

     call plotpar(4,ntim,time,'R2427','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2427_s(1:ntim,1))

     call plotpar(5,ntim,time,'R2424','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2424_s(1:ntim,1))

     call o_frame()

! R PLOT 11

     call plotback()

     call plotpar(1,ntim,time,'R2447','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2447_s(1:ntim,1))

     call plotpar(2,ntim,time,'R2552','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2552_s(1:ntim,1))

     call plotpar(3,ntim,time,'R2527','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2527_s(1:ntim,1))

     call plotpar(4,ntim,time,'R2525','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2525_s(1:ntim,1))

     call plotpar(5,ntim,time,'R2557','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2557_s(1:ntim,1))

     call o_frame()

! R PLOT 12

     call plotback()

     call plotpar(1,ntim,time,'R2662','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2662_s(1:ntim,1))

     call plotpar(2,ntim,time,'R2627','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2627_s(1:ntim,1))

     call plotpar(3,ntim,time,'R2626','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2626_s(1:ntim,1))

     call plotpar(4,ntim,time,'R2667','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2667_s(1:ntim,1))

     call plotpar(5,ntim,time,'R2772','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2772_s(1:ntim,1))

     call o_frame()

! R PLOT 13

     call plotback()

     call plotpar(1,ntim,time,'R2727','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r2727_s(1:ntim,1))

     call plotpar(2,ntim,time,'R0000','(g/kg)',0.,.1,.01)

     call bdline(1,ntim,r0000_s(1:ntim,1))

     call o_frame()

! E PLOT 1

     call plotback()

     call plotpar(1,ntim,time,'E1111','(#/g)',0.,5000.,500.)

     call bdline(1,ntim,e1111_s(1:ntim))

     call plotpar(2,ntim,time,'E1118','(#/g)',0.,5000.,500.)

     call bdline(1,ntim,e1118_s(1:ntim))

     call plotpar(3,ntim,time,'E8888','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8888_s(1:ntim))

     call plotpar(4,ntim,time,'E8882','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8882_s(1:ntim))

     call plotpar(5,ntim,time,'E1112','(#/g)',0.,5000.,500.)

     call bdline(1,ntim,e1112_s(1:ntim))

     call o_frame()

! E PLOT 2

     call plotback()

     call plotpar(1,ntim,time,'E1811','(#/g)',0.,5000.,500.)

     call bdline(1,ntim,e1811_s(1:ntim))

     call plotpar(2,ntim,time,'E1211','(#/g)',0.,5000.,500.)

     call bdline(1,ntim,e1211_s(1:ntim))

     call plotpar(3,ntim,time,'E8288','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8288_s(1:ntim))

     call plotpar(4,ntim,time,'E2222','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2222_s(1:ntim))

     call plotpar(5,ntim,time,'E5555','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e5555_s(1:ntim))

     call o_frame()

! E PLOT 3

     call plotback()

     call plotpar(1,ntim,time,'E6666','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e6666_s(1:ntim))

     call plotpar(2,ntim,time,'E7777','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e7777_s(1:ntim))

     call plotpar(3,ntim,time,'E3333','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3333_s(1:ntim))

     call plotpar(4,ntim,time,'E3335','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3335_s(1:ntim))

     call plotpar(5,ntim,time,'E4444','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e4444_s(1:ntim))

     call o_frame()

! E PLOT 4

     call plotback()

     call plotpar(1,ntim,time,'E4445','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e4445_s(1:ntim))

     call plotpar(2,ntim,time,'E3433','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3433_s(1:ntim))

     call plotpar(3,ntim,time,'E3444','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3444_s(1:ntim))

     call plotpar(4,ntim,time,'E3435','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3435_s(1:ntim))

     call plotpar(5,ntim,time,'E3445','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3445_s(1:ntim))

     call o_frame()

! E PLOT 5

     call plotback()

     call plotpar(1,ntim,time,'E3533','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3533_s(1:ntim))

     call plotpar(2,ntim,time,'E3633','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3633_s(1:ntim))

     call plotpar(3,ntim,time,'E3733','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e3733_s(1:ntim))

     call plotpar(4,ntim,time,'E4544','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e4544_s(1:ntim))

     call plotpar(5,ntim,time,'E4644','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e4644_s(1:ntim))

     call o_frame()

! E PLOT 6

     call plotback()

     call plotpar(1,ntim,time,'E4744','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e4744_s(1:ntim))

     call plotpar(2,ntim,time,'E5655','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e5655_s(1:ntim))

     call plotpar(3,ntim,time,'E5755','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e5755_s(1:ntim))

     call plotpar(4,ntim,time,'E6766','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e6766_s(1:ntim))

     call plotpar(5,ntim,time,'E1413','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1413_s(1:ntim))

     call o_frame()

! E PLOT 7

     call plotback()

     call plotpar(1,ntim,time,'E1411','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1411_s(1:ntim))

     call plotpar(2,ntim,time,'E1446','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1446_s(1:ntim))

     call plotpar(3,ntim,time,'E1513','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1513_s(1:ntim))

     call plotpar(4,ntim,time,'E1511','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1511_s(1:ntim))

     call plotpar(5,ntim,time,'E1556','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1556_s(1:ntim))

     call o_frame()

! E PLOT 8

     call plotback()

     call plotpar(1,ntim,time,'E1613','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1613_s(1:ntim))

     call plotpar(2,ntim,time,'E1611','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1611_s(1:ntim))

     call plotpar(3,ntim,time,'E8483','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8483_s(1:ntim))

     call plotpar(4,ntim,time,'E8488','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8488_s(1:ntim))

     call plotpar(5,ntim,time,'E8446','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8446_s(1:ntim))

     call o_frame()

! E PLOT 9

     call plotback()

     call plotpar(1,ntim,time,'E8583','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8583_s(1:ntim))

     call plotpar(2,ntim,time,'E8588','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8588_s(1:ntim))

     call plotpar(3,ntim,time,'E8556','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8556_s(1:ntim))

     call plotpar(4,ntim,time,'E8683','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8683_s(1:ntim))

     call plotpar(5,ntim,time,'E8688','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8688_s(1:ntim))

     call o_frame()

! E PLOT 10

     call plotback()

     call plotpar(1,ntim,time,'E1713','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1713_s(1:ntim))

     call plotpar(2,ntim,time,'E1711','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e1711_s(1:ntim))

     call plotpar(3,ntim,time,'E8783','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8783_s(1:ntim))

     call plotpar(4,ntim,time,'E8788','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e8788_s(1:ntim))

     call plotpar(5,ntim,time,'E2322','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2322_s(1:ntim))

     call o_frame()

! E PLOT 11

     call plotback()

     call plotpar(1,ntim,time,'E2327','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2327_s(1:ntim))

     call plotpar(2,ntim,time,'E2333','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2333_s(1:ntim))

     call plotpar(3,ntim,time,'E2337','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2337_s(1:ntim))

     call plotpar(4,ntim,time,'E2422','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2422_s(1:ntim))

     call plotpar(5,ntim,time,'E2427','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2427_s(1:ntim))

     call o_frame()

! E PLOT 12

     call plotback()

     call plotpar(1,ntim,time,'E2444','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2444_s(1:ntim))

     call plotpar(2,ntim,time,'E2447','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2447_s(1:ntim))

     call plotpar(3,ntim,time,'E2522','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2522_s(1:ntim))

     call plotpar(4,ntim,time,'E2527','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2527_s(1:ntim))

     call plotpar(5,ntim,time,'E2555','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2555_s(1:ntim))

     call o_frame()

! E PLOT 13

     call plotback()

     call plotpar(1,ntim,time,'E2557','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2557_s(1:ntim))

     call plotpar(2,ntim,time,'E2622','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2622_s(1:ntim))

     call plotpar(3,ntim,time,'E2627','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2627_s(1:ntim))

     call plotpar(4,ntim,time,'E2666','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2666_s(1:ntim))

     call plotpar(5,ntim,time,'E2667','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2667_s(1:ntim))

     call o_frame()

! E PLOT 14

     call plotback()

     call plotpar(1,ntim,time,'E2722','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2722_s(1:ntim))

     call plotpar(2,ntim,time,'E2727','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2727_s(1:ntim))

     call plotpar(3,ntim,time,'E2777','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e2777_s(1:ntim))

     call plotpar(4,ntim,time,'E0000','(#/g)',0.,1000.,100.)

     call bdline(1,ntim,e0000_s(1:ntim))

     call o_frame()

! END OF PLOTTING

     call o_clswk()
     stop 'parcel_done'

  endif

  return
end subroutine parcel_plot

!========================================================================

subroutine plotpar(ipanel,ntim,time,ylab1,ylab2,qmin,qmax,qinc)

  implicit none

  integer, intent(in) :: ipanel, ntim

  real, intent(in) :: time, qmin, qmax, qinc

  character(len=*), intent(in) :: ylab1, ylab2

  real, parameter :: sizelab = 11.
  
  real :: ybot, ymid, ytop, c1, x, y, p, dt_tick

  integer :: ntick, i, ii, nyint, j, ipcnt

  CHARACTER :: TITLE*88
  character(len=10)  :: numbr,numbr2

  real :: pmin, pmax, pinc, ybotp, yspc
  common / omicp / pmin, pmax, pinc, ybotp, yspc

  pmin = qmin
  pmax = qmax
  pinc = qinc

  ybot  = .06  + .19 * real(ipanel - 1)
  ytop  = ybot + .17
  ymid  = .5 * (ybot + ytop)

  ybotp = .075 + .19 * real(ipanel - 1)
  nyint = nint((pmax - pmin) / pinc)
  yspc = .15 / real(nyint)

  dt_tick = 100.
!  dt_tick = 10.
  ntick = nint(time / dt_tick)
  c1 = .76/real(ntick)

 ! panel borders

  call o_frstpt(.15,ybot)
  call o_vector(.95,ybot)
  call o_vector(.95,ytop)
  call o_vector(.15,ytop)
  call o_vector(.15,ybot)

! x label

  if (ipanel == 1) then
     call o_plchhq(.55,.01,'TIME (S)',sizelab,0.,0.)
  endif

! y labels

  call o_plchhq(.040,ymid,ylab1,sizelab,90.,0.)
  call o_plchhq(.065,ymid,ylab2,sizelab,90.,0.)

! short and long x-axis ticks and tick labels

  do i = 0,ntick

     x = .17 + c1 * real(i)
     call o_frstpt(x,ybot)
     call o_vector(x,ybot + .005)

     if (mod(i,5) == 0) then
        call o_vector(x,ybot + .010)

        if (ipanel == 1) then
           write (numbr,'(i5)') nint(i * dt_tick)
           call o_plchhq(x,.04,trim(adjustl(numbr)),sizelab, 0.,0.)
        endif
     endif

  enddo

! short and long y-axis ticks and tick labels

  do j = 0,nyint
     y = ybotp + real(j) * yspc

     call o_frstpt(.15 ,y)
     call o_vector(.155,y)
     call o_frstpt(.95 ,y)
     call o_vector(.945,y)

     p = pmin + real(j) * pinc
     ipcnt = nint(p / pinc)

     if (mod(ipcnt,5) == 0) then
        call o_frstpt(.15,y)
        call o_vector(.16,y)
        call o_frstpt(.95,y)
        call o_vector(.94,y)

        if (pinc * 5. >= .999) then
           write (numbr,'(i6)') nint(p)
        elseif (pinc * 5. >= .0999) then
           write (numbr,'(f5.1)') p
        elseif (pinc * 5. >= .00999) then
           write (numbr,'(f5.2)') p
        elseif (pinc * 5. >= .000999) then
           write (numbr,'(f6.3)') p
        else
           write (numbr,'(f7.4)') p
        endif
   
        call o_plchhq(.14,y,trim(adjustl(numbr)),sizelab, 0.,1.)

     endif
  enddo

  return
end subroutine plotpar

!========================================================================

subroutine bdline(ltype,ntim,pfld)

  implicit none

  integer, intent(in) :: ltype, ntim
  real, intent(in) :: pfld(ntim)

  real :: xpos(ntim), ypos(ntim)

  real :: dashlen, dashspc, xloc, yloc, remdist
  real :: xdist_nxtpt, ydist_nxtpt, dist_nxtpt
  real :: cosang, sinang

  integer :: ii, idash

  real :: pmin, pmax, pinc, ybotp, yspc
  common / omicp / pmin, pmax, pinc, ybotp, yspc

! x,y points in plot coordinates

  do ii = 1,ntim
     xpos(ii) = .17 + .76 * real(ii) / float(ntim)
     ypos(ii) = ybotp + yspc * (pfld(ii) - pmin) / pinc
  enddo

  if (ltype == 2) then
     dashlen = .002
     dashspc = .004
  elseif (ltype == 3) then
     dashlen = .015
     dashspc = .008
  endif

  if (ltype == 1) then

     call o_frstpt(xpos(1),ypos(1))
     do ii = 2,ntim
        call o_vector(xpos(ii),ypos(ii))
     enddo

  else

     call o_frstpt(xpos(1),ypos(1))
     xloc = xpos(1)
     yloc = ypos(1)
     remdist = dashlen
     idash = 1
     ii = 2

 47  continue

     xdist_nxtpt = xpos(ii) - xloc
     ydist_nxtpt = ypos(ii) - yloc

     dist_nxtpt = sqrt(xdist_nxtpt**2 + ydist_nxtpt**2)

     if (remdist < dist_nxtpt .and. idash == 1) then
        cosang = xdist_nxtpt / dist_nxtpt
        sinang = ydist_nxtpt / dist_nxtpt
        xloc = xloc + remdist * cosang
        yloc = yloc + remdist * sinang
        call o_vector(xloc,yloc)
        remdist = dashspc
        idash = 0
        go to 47
     elseif (remdist < dist_nxtpt .and. idash == 0) then
        cosang = xdist_nxtpt / dist_nxtpt
        sinang = ydist_nxtpt / dist_nxtpt
        xloc = xloc + remdist * cosang
        yloc = yloc + remdist * sinang
        call o_frstpt(xloc,yloc)
        remdist = dashlen
        idash = 1
        go to 47
     else
        remdist = remdist - dist_nxtpt
        xloc = xpos(ii)
        yloc = ypos(ii)
        ii = ii + 1
        if (ii > ntim) return
        go to 47
     endif

  endif

  return
end subroutine bdline


