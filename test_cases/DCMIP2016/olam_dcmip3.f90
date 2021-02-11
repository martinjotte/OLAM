subroutine olam_dcmip_init()

use mem_basic,   only: vc, vmc, wc, wmc, rho, thil, theta, tair, &
                       press, rr_w, rr_v
use mem_ijtabs,  only: jtab_v, jtab_w, itab_v, jtv_init, jtw_init
use mem_grid,    only: mza, mva, mwa, zt, dzt, xev, yev, zev, vnx, vny, vnz, lpw, &
                       glatw, glonw, lpv, zm, glatv, glonv, gdz_belo8, gdz_abov8
use mem_addsc,   only: addsc
use misc_coms,   only: io6, iparallel
use consts_coms, only: p00, p00i, p00kord, rocp, alvlocp, cvocp, eps_vapi, &
                       xscale, pio180, erad, r8
use mem_micro,   only: rr_c, con_c, cldnum
use micro_coms,  only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, &
   test1_advection_hadley, &
   test1_advection_orography, &
   test2_steady_state_mountain, &
   test2_schaer_mountain, &
   test3_gravity_wave

use baroclinic_wave, only: baroclinic_wave_test

use test2xdeep, only: deep_mountain_wave

use dcmip_initial_conditions_test_4, only: &
   test4_baroclinic_wave

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                        mpi_send_v, mpi_recv_v

use obnd,         only: lbcopy_v, lbcopy_w

use test41alt, only: baroclinic_instability_alt

use tropical_cyclone, only: tropical_cyclone_test

use supercell_testm, only: supercell_init, supercell_test

use terminator, only: initial_value_terminator, ctend1, ctend2

use oname_coms, only: nl
use therm_lib,  only: rhovsl

implicit none

real(8) :: time0

real(8) :: latrad,lonrad,latdeg,londeg
real(8) :: p,t,phis,ps,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr
real(8) :: hyam, hybm, gc, xscal, cl, cl2

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0, rm0, rt0, thetav

real(8) :: rho_tot(mza), pkhyd

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,rcoords,cfv,shear,moist,iter,ka,kbc, &
           subtest, deep, pertt, pert

logical :: hybrid_eta

real :: exner, temp, ccn

! In case standard OLAM hydrostatic balance is carried out, choose the middle
! level of the model (as counted by vertical index) to be an internal pressure
! boundary condition.  Ideally, this will be near the middle of the atmosphere
! in terms of mass.

kbc = mza / 2

zcoords = 1
rcoords = 1
time0   = 0.0d0
xscal  = xscale

hybrid_eta = .false.

cfv = 0
shear = 0

mrl = 1

! Allocate chemical tendency arrays for baroclinic wave test

if (nl%test_case == 110 .or. &
    nl%test_case == 111 .or. &
    nl%test_case == 112 .or. &
    nl%test_case == 113 .or. &
    nl%test_case == 114) then

   allocate (ctend1(mza,mwa), ctend2(mza,mwa))
endif

! Pre-initialization for supercell test case

if (nl%test_case == 131) call supercell_init()

! Horizontal loop over model T points

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------
   latrad = pio180 * glatw(iw)
   lonrad = pio180 * glonw(iw)
   latdeg =          glatw(iw)
   londeg =          glonw(iw)

   if (glonw(iw) < 0.) lonrad = pio180 * (glonw(iw) + 360.)
   if (glonw(iw) < 0.) londeg =          (glonw(iw) + 360.)

   ka = lpw(iw)

   do k = ka,mza
      p = press(k,iw)
      zm0 = zm(k)
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==============================================================================
! DCMIP-2012 TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==============================================================================

         call test1_advection_deformation(lonrad,latrad,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2,q3,q4,time0)

         call test1_advection_deformation(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2
         if (nl%naddsc >= 3) addsc(3)%sclp(k,iw) = q3
         if (nl%naddsc >= 4) addsc(4)%sclp(k,iw) = q4
         if (nl%naddsc >= 5) addsc(5)%sclp(k,iw) = 1.0d0

      elseif (nl%test_case == 12) then

!==============================================================================
! DCMIP-2012 TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==============================================================================

         call test1_advection_hadley(lonrad,latrad,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,time0)

         call test1_advection_hadley(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1

      elseif (nl%test_case == 13) then

!==============================================================================
! DCMIP-2012 TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS
! IN THE PRESENCE OF OROGRAPHY
!==============================================================================

         call test1_advection_orography(lonrad,latrad,p,zm0,zcoords, &
            cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2,q3,q4)

         call test1_advection_orography(lonrad,latrad,p,zt0,zcoords, &
            cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2
         if (nl%naddsc >= 3) addsc(3)%sclp(k,iw) = q3
         if (nl%naddsc >= 4) addsc(4)%sclp(k,iw) = q4

! For horizontal (not terrain-following) model levels

         wm0 = 0.0d0
         wt0 = 0.0d0

      elseif (nl%test_case == 200 .or. &
              nl%test_case == 201) then

!==============================================================================
! DCMIP-2012 Test 2-0:  Steady-State Atmosphere at Rest
! in the Presence of Orography
!==============================================================================

         call test2_steady_state_mountain(lonrad,latrad,p,zm0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wm0, &
            t,phis,ps,rhom0,q)

         call test2_steady_state_mountain(lonrad,latrad,p,zt0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 21 .or. &
              nl%test_case == 22) then

!==============================================================================
! DCMIP-2012 Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over
! a Schaer-type Mountain
!==============================================================================

         if (nl%test_case == 21) then
            shear = 0
            subtest = 1
         else  ! nl%test_case == 22
            shear = 1
            subtest = 2
         endif

         rm0 = 6371220. / xscale + zm(k)
         rt0 = 6371220. / xscale + zt(k)

         call deep_mountain_wave(subtest,Xscal,lonrad,latrad,p,rm0,rcoords, &
                                 u0,v0,wm0,t,phis,ps,rhom0)

         call deep_mountain_wave(subtest,Xscal,lonrad,latrad,p,rt0,rcoords, &
                                 u0,v0,wt0,t,phis,ps,rhot0)

         q = 0.0d0

! The following 2 calls (commented out) are for shallow atmosphere

!         call test2_schaer_mountain(lonrad,latrad,p,zm0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wm0, &
!            t,phis,ps,rhom0,q)

!         call test2_schaer_mountain(lonrad,latrad,p,zt0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wt0, &
!            t,phis,ps,rhot0,q)

     elseif (nl%test_case == 31) then

!==============================================================================
! DCMIP-2012 Test 3-1 - GRAVITY WAVES
!==============================================================================

         call test3_gravity_wave(lonrad,latrad,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q)

         call test3_gravity_wave(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 410 .or. &
              nl%test_case == 411 .or. &
              nl%test_case == 412 .or. &
              nl%test_case == 413) then

!==============================================================================
! DCMIP-2012 Test 41x - DRY BAROCLINIC INSTABILITY
!==============================================================================

! New code for deep atmosphere case:

         deep = 1

         call baroclinic_instability_alt(deep,Xscal,lonrad,latrad,p,zm0,zcoords, &
                                         u0,v0,wm0, &
                                         t,phis,ps,rhom0,q1,q2)

         call baroclinic_instability_alt(deep,Xscal,lonrad,latrad,p,zt0,zcoords, &
                                         u0,v0,wt0, &
                                         t,phis,ps,rhot0,q1,q2)

         q = 0.0d0

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2

      elseif (nl%test_case == 42 .or. &
              nl%test_case == 43) then

!==============================================================================
! DCMIP-2012 Tests 42,43 - MOIST BAROCLINIC INSTABILITY
!==============================================================================

         moist = 1

         call test4_baroclinic_wave(moist,Xscal,lonrad,latrad,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2)
 
         call test4_baroclinic_wave(moist,Xscal,lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2

      elseif (nl%test_case == 110 .or. &
              nl%test_case == 111 .or. &
              nl%test_case == 112 .or. &
              nl%test_case == 113 .or. &
              nl%test_case == 114) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-1 - Dry or Moist Baroclinic Wave 
!==============================================================================

         deep = 1  ! 0=no, 1=yes

         if     (nl%test_case == 110) then  ! Dry case
            pertt = 0 ! exponential type perturbation
            moist = 0 ! 0=no, 1=yes
         elseif (nl%test_case == 111) then
            pertt = 0 ! exponential type perturbation
            moist = 1 ! 0=no, 1=yes
         elseif (nl%test_case == 112) then
            pertt = 0 ! exponential type perturbation
            moist = 1 ! 0=no, 1=yes
         elseif (nl%test_case == 113) then
            pertt = 1 ! streamfunction type perturbation)
            moist = 1 ! 0=no, 1=yes
         else ! (nl%test_case == 114)
            pertt = 1 ! streamfunction type perturbation)
            moist = 1 ! 0=no, 1=yes
         endif

         call baroclinic_wave_test(deep,moist,pertt,Xscal,lonrad,latrad,p,zt0,zcoords, &
                                   u0,v0,t,thetav,phis,ps,rhot0,q)

         wm0 = 0.
         rhom0 = rhot0 ! Exact value not needed for this case

         call initial_value_Terminator( latdeg, londeg, cl, cl2 )

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = real(cl)
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = real(cl2)

      elseif (nl%test_case == 121 .or. &
              nl%test_case == 122) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-2 - Tropical Cyclone
!==============================================================================

         call tropical_cyclone_test(lonrad,latrad,p,zt0,zcoords, &
                                    u0,v0,t,thetav,phis,ps,rhot0,q)

         ! q is specific humidity. Convert to mixing ratio
         q = q / (1. - q)

         wm0 = 0.
         rhom0 = rhot0 ! Exact value not needed for this case

      elseif (nl%test_case == 131) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-3 - Supercell 
!==============================================================================

         pert = 1

         call supercell_test(lonrad,latrad,p,zt0,zcoords,u0,v0,t,thetav,ps,rhot0,q,pert)

         wm0 = 0.
         rhom0 = rhot0 ! Exact value not needed for this case

      endif

      wc(k,iw) = wm0
      wmc(k,iw) = wm0 * rhom0

!  Assume: q is mixing ratio, rhot0 include dry air and water vapor

      rho(k,iw) = rhot0 / (1. + q)
      press(k,iw) = p
      tair (k,iw) = t
      theta(k,iw) = t * (p00 / p) ** rocp
      thil(k,iw) = theta(k,iw)
      rr_w(k,iw) = q
      rr_v(k,iw) = q

   enddo

   ! Iterative hydrostatic balance procedure
   ! Procedure holds THETA and RR_W constant, after THETA is computed above

!   go to 100  ! BYPASS THE HYDROSTATIC BALANCE PROCEDURE

   do iter = 1,100

!  Compute density for all levels

      do k = ka, mza

         if (miclevel == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00kord / theta(k,iw)
            rho_tot(k) = rho(k,iw)
         elseif (miclevel == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                        ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )
            rho_tot(k) = rho(k,iw) * (1. + rr_v(k,iw))
         else
            exner = (real(press(k,iw)) * p00i) ** rocp   ! Defined WITHOUT CP factor
            temp  = exner * theta(k,iw)

            rr_c(k,iw) = max(0., rr_w(k,iw) - rhovsl(temp-273.15) / real(rho(k,iw)))
            rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

            rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                        ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )

            rho_tot(k) = rho(k,iw) * (1. + rr_w(k,iw))
        endif

      enddo

! Integrate hydrostatic equation upward and downward from kbc level.
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
! during iteration.  Use weighting to damp oscillations

      do k = kbc+1,mza
         pkhyd = press(k-1,iw) &
               - gdz_belo8(k-1) * rho_tot(k-1) - gdz_abov8(k-1) * rho_tot(k)
         press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)
      enddo

      do k = kbc-1,ka,-1
         pkhyd = press(k+1,iw) &
               + gdz_belo8(k) * rho_tot(k) + gdz_abov8(k) * rho_tot(k+1)
         press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)
      enddo

   enddo

   do k = ka, mza
      tair(k,iw) = theta(k,iw) * (real(press(k,iw)) * p00i) ** rocp
      if (miclevel <= 1) then
         thil(k,iw) = theta(k,iw)
      else
         thil(k,iw) = theta(k,iw) / (1. + alvlocp * rr_c(k,iw) / &
                                         ((1.0 + rr_v(k,iw)) * max(tair(k,iw),253.)))
      endif
   enddo

   thil(1:ka-1,iw) = thil(ka,iw)

! If there is cloud condensate, initialize con_c if prognosed

   if (miclevel == 3 .and. jnmb(1) == 5) then
      if (ccnparm > 1.e6) then
         ccn = ccnparm
      else
         ccn = cldnum(iw)
      endif

      do k = ka, mza
         if (rr_c(k,iw) > rxmin(1)) then
            con_c(k,iw) = ccn * real(rho(k,iw)) * zfactor_ccn(k)
         else
            con_c(k,iw) = 0.0
         endif
      enddo
   endif

enddo

100 continue

! MPI parallel send/recv of W group

if (iparallel == 1) then
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=rho, d2=press)

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   lonrad = pio180 * glonv(iv)
   latrad = pio180 * glatv(iv)

   if (glonv(iv) < 0.) lonrad = pio180 * (glonv(iv) + 360.)

   do k = lpv(iv),mza
      p = 1.0d5
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==============================================================================
! DCMIP-2012 TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==============================================================================

         call test1_advection_deformation(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

      elseif (nl%test_case == 12) then

!==============================================================================
! DCMIP-2012 TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==============================================================================

         call test1_advection_hadley(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

      elseif (nl%test_case == 13) then

!==============================================================================
! DCMIP-2012 TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS
! IN THE PRESENCE OF OROGRAPHY
!==============================================================================

         call test1_advection_orography(lonrad,latrad,p,zt0,zcoords, &
            cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4)

! Alteration for OLAM: Set velocity at mountain heights = 0.

            if (zt0 < 2000.d0) then
               u0 = 0.0d0
               v0 = 0.0d0
            endif

      elseif (nl%test_case == 200 .or. &
              nl%test_case == 201) then

!==============================================================================
! DCMIP-2012 Test 2-0:  Steady-State Atmosphere at Rest
! in the Presence of Orography
!==============================================================================

         call test2_steady_state_mountain(lonrad,latrad,p,zt0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 21 .or. &
              nl%test_case == 22) then

!==============================================================================
! DCMIP-2012 Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves
! over a Schaer-type Mountain
!==============================================================================

         if (nl%test_case == 21) then
            shear = 0
            subtest = 1
         else  ! nl%test_case == 22
            shear = 1
            subtest = 2
         endif

         rt0 = 6371220. / xscale + zt(k)

         call deep_mountain_wave(subtest,Xscal,lonrad,latrad,p,rt0,rcoords, &
                                 u0,v0,wt0,t,phis,ps,rhot0)

! The following call (commented out) is for shallow atmosphere approx.

!         call test2_schaer_mountain(lonrad,latrad,p,zt0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wt0, &
!            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 31) then

!==============================================================================
! DCMIP-2012 Test 3-1 - GRAVITY WAVES
!==============================================================================

         call test3_gravity_wave(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 410 .or. &
              nl%test_case == 411 .or. &
              nl%test_case == 412 .or. &
              nl%test_case == 413) then

!==============================================================================
! DCMIP-2012 Test 41x - DRY BAROCLINIC INSTABILITY
!==============================================================================

! New code for deep atmosphere case:

         deep = 1

         call baroclinic_instability_alt(deep,Xscal,lonrad,latrad,p,zt0,zcoords, &
                                         u0,v0,wt0, &
                                         t,phis,ps,rhot0,q1,q2)

      elseif (nl%test_case == 42 .or. &
              nl%test_case == 43) then

!==============================================================================
! DCMIP-2012 Tests 42,43 - MOIST BAROCLINIC INSTABILITY
!==============================================================================

         moist = 1

         call test4_baroclinic_wave(moist,Xscal,lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2)

      elseif (nl%test_case == 110 .or. &
              nl%test_case == 111 .or. &
              nl%test_case == 112 .or. &
              nl%test_case == 113 .or. &
              nl%test_case == 114) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-1 - Dry or Moist Baroclinic Wave 
!==============================================================================

         deep = 1  ! 0=no, 1=yes

         if     (nl%test_case == 110) then  ! Dry case
            pertt = 0 ! exponential type perturbation
            moist = 0 ! 0=no, 1=yes
         elseif (nl%test_case == 111) then
            pertt = 0 ! exponential type perturbation
            moist = 1 ! 0=no, 1=yes
         elseif (nl%test_case == 112) then
            pertt = 0 ! exponential type perturbation
            moist = 1 ! 0=no, 1=yes
         elseif (nl%test_case == 113) then
            pertt = 1 ! streamfunction type perturbation)
            moist = 1 ! 0=no, 1=yes
         else ! (nl%test_case == 114)
            pertt = 1 ! streamfunction type perturbation)
            moist = 1 ! 0=no, 1=yes
         endif

         call baroclinic_wave_test(deep,moist,pertt,Xscal,lonrad,latrad,p,zt0,zcoords, &
                                   u0,v0,t,thetav,phis,ps,rhot0,q)

      elseif (nl%test_case == 121 .or. &
              nl%test_case == 122) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-2 - Tropical Cyclone 
!==============================================================================

         call tropical_cyclone_test(lonrad,latrad,p,zt0,zcoords, &
                                    u0,v0,t,thetav,phis,ps,rhot0,q)

      elseif (nl%test_case == 131) then

!==============================================================================
! DCMIP-2016 TEST CASE 1-3 - Supercell 
!==============================================================================

         pert = 1

         call supercell_test(lonrad,latrad,p,zt0,zcoords,u0,v0,t,thetav,ps,rhot0,q,pert)

      endif

      raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

      if (raxis > 1.e3) then
         uv01dr = -v0 * zev(iv) / erad  ! radially outward from axis

         uv01dx = (-u0 * yev(iv) + uv01dr * xev(iv)) / raxis 
         uv01dy = ( u0 * xev(iv) + uv01dr * yev(iv)) / raxis 
         uv01dz =   v0 * raxis / erad 

         vc(k,iv) = uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv)
      else
         vc(k,iv) = 0.
      endif

      vmc(k,iv) = vc(k,iv) * rhot0

   enddo

   ! For below-ground points, set VMC and VC to zero

   ka = lpv(iv)
   vmc(1:ka-1,iv) = 0.
   vc (1:ka-1,iv) = 0.
enddo

! MPI parallel send/recv of V group

if (iparallel == 1) then
   call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
endif

! LBC copy of VMC, VC

call lbcopy_v(1, vmc=vmc, vc=vc)

end subroutine olam_dcmip_init

!==========================================================================

subroutine dcmip_save_initfields()

use mem_basic,  only: vxe, vye, vze, thil, rr_w
use mem_grid,   only: mza, mwa

use supercell_testm, only: vxe_init, vye_init, vze_init, thil_init, rr_w_init

implicit none

  ! Allocate OLAM arrays for storing initial velocity and thil

  allocate ( vxe_init(mza,mwa))
  allocate ( vye_init(mza,mwa))
  allocate ( vze_init(mza,mwa))
  allocate (thil_init(mza,mwa))
  allocate (rr_w_init(mza,mwa))

  ! Store initial values of velocity and thil

   vxe_init(:,:) =  vxe(:,:)
   vye_init(:,:) =  vye(:,:)
   vze_init(:,:) =  vze(:,:)
  thil_init(:,:) = thil(:,:)
  rr_w_init(:,:) = rr_w(:,:)

end subroutine dcmip_save_initfields

!==========================================================================

subroutine olam_dcmip_prescribedflow(time0,vmsca,wmsca,rho_old)

use mem_basic,  only: vc, vmc, wc, wmc, rho, thil, theta, press, rr_w, rr_v, &
                      vxe, vye, vze
use mem_ijtabs, only: itab_v, jtab_v, jtab_w, jtw_prog, itab_w
use mem_grid,   only: mza, mva, mwa, zt, xev, yev, zev, vnx, vny, vnz, lpw, &
                      glatw, glonw, lpv, zm, glatv, glonv, xew, yew, zew, dzt, &
                      arw, arv, volti
use mem_addsc,  only: addsc
use consts_coms,only: p00, rocp, erad, pio180, p00i
use misc_coms,  only: dtlong, dtsm, iparallel
use vel_t3d,    only: diagvel_t3d
use oname_coms, only: nl

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
use obnd,         only: lbcopy_w

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, test1_advection_deformation_w, &
   test1_advection_deformation_uv, test1_advection_hadley

implicit none

real(8), intent(in) :: time0

real, intent(out) :: vmsca(mza,mva), wmsca(mza,mwa)

!real(8), intent(out) :: rho_old(mza,mwa) ! density at beginning of timestep [kg/m^3]
real, intent(out) :: rho_old(mza,mwa) ! density at beginning of timestep [kg/m^3]

real(8) :: lonrad,latrad,p,z,t,phis,ps,rho0,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0, dts

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,test,ka,kb,jv

real(8) :: pv(mza), zmv(mza), wcv(mza), uz(mza), vz(mza), rhoz(mza)

real    :: hflux_rho(mza), dirv

zcoords = 1
mrl     = 1

!----------------------------------------------------------------------
!$omp parallel private(pv,zmv,wcv,uz,vz,rhoz,hflux_rho)
!$omp do private(lonrad,latrad,ka,kb,k,p,zm0,u0,v0,wm0, &
!$omp            t,phis,ps,rhom0,q,q1)
do iw = 2,mwa  ! do for all points in subdomain
!----------------------------------------------------------------------
   lonrad = pio180 * glonw(iw)
   latrad = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lonrad = pio180 * (glonw(iw) + 360.)

   ka = lpw(iw)

   wc     (mza,iw) = 0.
   wmc    (mza,iw) = 0.
   wmsca  (mza,iw) = 0.
   rho_old(mza,iw) = rho(mza,iw)
   rho_old(ka ,iw) = rho(ka ,iw)

   if (nl%test_case == 11) then

!==============================================================================
! DCMIP-2012 TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==============================================================================

      zmv = zm
      kb  = mza-1

      call test1_advection_deformation_w(lonrad,latrad,ka,kb,pv,zmv,zcoords, &
                                         wcv,rhoz,time0)

      do k = ka, mza-1
         wc   (k,iw) = wcv(k)
         wmc  (k,iw) = wcv(k) * rhoz(k)
         wmsca(k,iw) = wmc(k,iw) * arw(k,iw)
      enddo

   elseif (nl%test_case == 12) then

!==============================================================================
! DCMIP-2012 TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==============================================================================

      do k = ka, mza-1
         zm0 = zm(k)

         call test1_advection_hadley(lonrad,latrad,p,zm0,zcoords, &
            u0,v0,wm0,t,phis,ps,rhom0,q,q1,time0)

         wc   (k,iw) = wm0
         wmc  (k,iw) = wm0 * rhom0
         wmsca(k,iw) = wmc(k,iw) * arw(k,iw)
      enddo

   endif

   do k = ka, mza
      rho_old(k,iw) = rho(k,iw)
   enddo

enddo
!$omp end do

!----------------------------------------------------------------------
!$omp do private(iw1,iw2,lonrad,latrad,kb,k,p,zt0,     &
!$omp            u0,v0,wt0,t,phis,ps,rhot0,q,q1,raxis, &
!$omp            uv01dr,uv01dx,uv01dy,uv01dz)
do iv = 2,mva
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   lonrad = pio180 * glonv(iv)
   latrad = pio180 * glatv(iv)

   kb = lpv(iv)

   if (glonv(iv) < 0.) lonrad = pio180 * (glonv(iv) + 360.)

   if (nl%test_case == 11) then

!==============================================================================
! DCMIP-2012 TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==============================================================================

      zmv = zt

      call test1_advection_deformation_uv(lonrad,latrad,kb,mza,pv,zmv, &
                                          zcoords,uz,vz,rhoz,time0)

   elseif (nl%test_case == 12) then

!==============================================================================
! DCMIP-2012 TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==============================================================================

      do k = lpv(iv), mza
         zt0 = zt(k)

         call test1_advection_hadley(lonrad,latrad,p,zt0,zcoords, &
            u0,v0,wt0,t,phis,ps,rhot0,q,q1,time0)

         uz  (k) = u0
         vz  (k) = v0
         rhoz(k) = rhot0
      enddo

   endif

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

   if (raxis > 1.e3) then

      do k = kb, mza
         uv01dr = -vz(k) * zev(iv) / erad  ! radially outward from axis

         uv01dx = (-uz(k) * yev(iv) + uv01dr * xev(iv)) / raxis
         uv01dy = ( uz(k) * xev(iv) + uv01dr * yev(iv)) / raxis
         uv01dz =   vz(k) * raxis / erad

         vc(k,iv) = uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv)
      enddo

   else

      do k = kb, mza
         vc(k,iv) = 0.
      enddo

   endif

   do k = kb, mza
      vmc  (k,iv) = vc (k,iv) * real( rhoz(k) )
      vmsca(k,iv) = vmc(k,iv) * arv(k,iv)
   enddo

enddo
!$omp end do

! Even though the idealized conditions have no divergence, when applied to the
! grid slight divergences are produced. Scalar advection requires density
! be consistent with the scalar mixing ratio, so here we advect density too.

!$omp do private(iw,k,jv,iv,dirv,dts)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

   hflux_rho = 0.0

   do jv = 1, itab_w(iw)%npoly
      iv   = itab_w(iw)%iv(jv)
      dirv = itab_w(iw)%dirv(jv)

      ! Loop over T levels
      do k = lpv(iv), mza
         hflux_rho(k) = hflux_rho(k) + dirv * vmsca(k,iv)
      enddo
   enddo

   dts = dtsm(itab_w(iw)%mrlw)

   do k = lpw(iw), mza
      rho(k,iw) = rho(k,iw) + dts * volti(k,iw) &
                            * (hflux_rho(k) + wmsca(k-1,iw) - wmsca(k,iw))
   enddo

   do k = 1, lpw(iw)-1
      rho(k,iw) = rho(lpw(iw),iw)
   enddo

end do
!$omp end do
!$omp end parallel

! Compute Earth-Cartesian velocities

call diagvel_t3d(1)

if (iparallel == 1) call mpi_send_w(mrl, dvara1=rho)

! Make sure no physics tendencies have been added to tracers

if (nl%naddsc >= 1) addsc(1)%sclt = 0.
if (nl%naddsc >= 2) addsc(2)%sclt = 0.
if (nl%naddsc >= 3) addsc(3)%sclt = 0.
if (nl%naddsc >= 4) addsc(4)%sclt = 0.
if (nl%naddsc >= 5) addsc(5)%sclt = 0.

if (iparallel == 1) call mpi_recv_w(mrl, dvara1=rho)
call lbcopy_w(1, d1=rho)

end subroutine olam_dcmip_prescribedflow

!==============================================================================

subroutine olam_dcmip2016_phys()

! Subroutine olam_dcmip2016_phys is called from subroutine timestep and in turn
! calls subroutine dcmip2016_physics that provides tendencies or updates to model
! prognostic fields for DCMIP2016 test cases.

use mem_basic,   only: vc, vmc, wc, wmc, rho, thil, theta, press, &
                       rr_w, rr_v, vxe, vye, vze
use mem_ijtabs,  only: itab_v, jtab_v, jtab_w, jtw_prog, jtv_prog
use mem_grid,    only: mza, mva, mwa, zt, xev, yev, zev, vnx, vny, vnz, lpw, &
                       glatw, glonw, lpv, zm, glatv, glonv, xew, yew, zew, dzt, &
                       gravm, gravt
use mem_addsc,   only: addsc
use consts_coms, only: r8, p00, rocp, erad, eradi, pio180, p00i
use misc_coms,   only: dtlong, time8p, dtlm, io6, iparallel
use mem_micro,   only: pcprr, accpr, rr_c, rr_r
use mem_tend,    only: thilt, rr_wt, rr_ct, rr_rt, vmxet, vmyet, vmzet
use oname_coms,  only: nl

implicit none

real(r8) :: rhodrycol(mza)
real(r8) ::  exnercol(mza)
real(r8) ::      tcol(mza)
real(r8) ::      pcol(mza)
real(r8) ::      zcol(mza)
real(r8) ::     zicol(mza+1)
real(r8) ::  thetacol(mza), thetacol0(mza)
real(r8) ::     qvcol(mza),    qvcol0(mza)
real(r8) ::     qccol(mza),    qccol0(mza)
real(r8) ::     qrcol(mza),    qrcol0(mza)
real(r8) ::      ucol(mza),     ucol0(mza)
real(r8) ::      vcol(mza),     vcol0(mza)

real(r8) :: dudt, dvdt, dvrdt, raxis, qold, qnew
real(r8) :: lonrad,latrad,dtime
real(r8) :: precl, vx, vy, vz, za

integer :: mrl,j,iw,k,iv,iw1,iw2,nz,ka
integer :: test,pbl_type,prec_type

mrl = 1
nz = mza
dtime = dtlong

! Define flags for dcmip2016_physics subroutine

!    pbl_type = 0 ! Reed-Jablonowski PBL
!    pbl_type = 1 ! Modified Bryan PBL
!    pbl_type = 2 ! Do not use surface fluxes or PBL
!   prec_type = 0 ! Default Kessler physics
!   prec_type = 1 ! Reed-Jablonowski microphysics

if     (nl%test_case == 110) then  ! Dry case
        test = 1 ! Baroclinic wave test
    pbl_type = 2 ! Do not use surface fluxes or PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 111) then
        test = 1 ! Baroclinic wave test
    pbl_type = 2 ! Do not use surface fluxes or PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 112) then 
        test = 1 ! Baroclinic wave test
    pbl_type = 0 ! Reed-Jablonowski PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 113) then 
        test = 1 ! Baroclinic wave est
    pbl_type = 2 ! Do not use surface fluxes or PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 114) then 
        test = 1 ! Baroclinic wave test
    pbl_type = 0 ! Reed-Jablonowski PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 121) then 
        test = 2 ! Tropical cyclone test
    pbl_type = 0 ! Reed-Jablonowski PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 122) then 
        test = 2 ! Tropical cyclone test
    pbl_type = 1 ! Modified Bryan PBL
   prec_type = 0 ! Default Kessler physics
elseif (nl%test_case == 131) then
        test = 3 ! Supercell test
    pbl_type = 2 ! Do not use surface fluxes or PBL
   prec_type = 0 ! Default Kessler physics
endif
 
nz = mza

do k = 2,mza
  zcol(k) = zt(k)
  zicol(k) = zm(k-1) ! Different indexing from zm
enddo
zcol(1) = zt(1)
zicol(mza+1) = zm(mza)

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   lonrad = pio180 * glonw(iw)
   latrad = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lonrad = pio180 * (glonw(iw) + 360.)

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Fill vertical columns in preparation for call to subroutine simple_physics

   ka = lpw(iw)
   do k = ka,mza

      vx = vxe(k,iw)
      vy = vye(k,iw)
      vz = vze(k,iw)

      if (raxis > 1.e3) then
         ucol(k) = (vy * xew(iw) - vx * yew(iw)) / raxis
         vcol(k) = vz * raxis * eradi &
           - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad)
      else
         ucol(k) = 0.
         vcol(k) = 0.
      endif

      ! Apply fixer to prevent negative moisture values before calling DCMIP physics

      rr_v(k,iw) = max(rr_v(k,iw), 0.0)
      rr_c(k,iw) = max(rr_c(k,iw), 0.0)
      rr_r(k,iw) = max(rr_r(k,iw), 0.0)

       exnercol(k) = (real(press(k,iw)) * p00i) ** rocp      ! Defined WITHOUT CP factor
      rhodrycol(k) =    rho(k,iw)
       thetacol(k) =   thil(k,iw) ! THIL = THETA for DCMIP 2016 cases; THIL is prognostic
           pcol(k) =  press(k,iw)
          qvcol(k) =   rr_v(k,iw)
          qccol(k) =   rr_c(k,iw)
          qrcol(k) =   rr_r(k,iw)

      ! Save a copy of quantities that will change in physics update

      thetacol0(k) = thetacol(k)
         qvcol0(k) =    qvcol(k)
         qccol0(k) =    qccol(k)
         qrcol0(k) =    qrcol(k)
          ucol0(k) =     ucol(k)
          vcol0(k) =     vcol(k)

   enddo

   precl = 0.

   call DCMIP2016_PHYSICS(test, ucol, vcol, pcol, qvcol, qccol, qrcol, &
                          rhodrycol, thetacol, exnercol, zcol, zicol, &
                          dtime, latrad, nz, ka, iw, precl, pbl_type, prec_type)

! Apply physics updates to model tendency arrays for prognostic variables

   do k = ka,mza
      qold = qvcol0(k) + qccol0(k) + qrcol0(k)
      qnew = qvcol (k) + qccol (k) + qrcol (k)

      ! Convert mixing ratio changes to density-weighted tendencies used by OLAM

      rr_wt(k,iw) = rr_wt(k,iw) + (qnew - qold)                * rhodrycol(k) / dtime
      rr_ct(k,iw) = rr_ct(k,iw) + (qccol(k) - qccol0(k))       * rhodrycol(k) / dtime
      rr_rt(k,iw) = rr_rt(k,iw) + (qrcol(k) - qrcol0(k))       * rhodrycol(k) / dtime
      thilt(k,iw) = thilt(k,iw) + (thetacol(k) - thetacol0(k)) * rhodrycol(k) / dtime

      ! Convert velocity change to momentum tendencies used by OLAM

      if (raxis > 1.e3) then
         dudt = (ucol(k) - ucol0(k)) / dtime ! Zonal velocity tendency
         dvdt = (vcol(k) - vcol0(k)) / dtime ! Meridional velocity tendency

         dvrdt = -dvdt * zew(iw) * eradi  ! radially outward from axis

         vmxet(k,iw) = vmxet(k,iw) + rho(k,iw) * (-dudt * yew(iw) + dvrdt * xew(iw)) / raxis
         vmyet(k,iw) = vmyet(k,iw) + rho(k,iw) * ( dudt * xew(iw) + dvrdt * yew(iw)) / raxis
         vmzet(k,iw) = vmzet(k,iw) + rho(k,iw) *   dvdt * raxis * eradi
      endif
   enddo

! Convert precip rate from m/s to mm/s (or equivalent kg/m^2/s.
! Add single-timestep contribution to accumulated precipitation

   pcprr(iw) = precl * 1000.
   accpr(iw) = accpr(iw) + dtime * pcprr(iw)

enddo

end subroutine olam_dcmip2016_phys

!===============================================================================

subroutine thermo_dcmip()

use mem_ijtabs, only: jtab_w, istp, mrl_endl, jtw_prog
use mem_grid,   only: lpw, mza
use consts_coms,only: p00i, rocp
use mem_basic,  only: theta, thil, tair, rr_v, rr_w, press
use mem_micro,  only: rr_c, rr_r

implicit none

integer iw,j,mrl,k

! Horizontal loop over W/T points

!-------------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw,k)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!-------------------------------------------------------------------------

   do k = lpw(iw),mza
      theta(k,iw) = thil(k,iw)
      tair (k,iw) = thil(k,iw) * (real(press(k,iw)) * p00i) ** rocp
      rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw) - rr_r(k,iw)
   enddo

enddo
!$omp end parallel do
endif

end subroutine thermo_dcmip

!==============================================================================

subroutine olam_dcmip_terminator()

! Subroutine olam_dcmip_terminator is called from subroutine timestep and in
! turn calls subroutine terminator that provides updates to chemistry
! tendencies for DCMIP test cases.

use mem_grid,   only: mza, mwa, lpw, glatw, glonw
use mem_basic,  only: rho
use mem_ijtabs, only: jtab_w, jtw_prog
use mem_addsc,  only: addsc
use misc_coms,  only: time8p, dtlong
use terminator, only: ctend1, ctend2, tendency_terminator
use oname_coms, only: nl

implicit none

integer :: iw, k, ka, mrl, j

real(8) :: dt, latdeg, londeg, cl, cl2, cl_f, cl2_f

mrl = 1

! terminator coupling timestep

dt = 900._8
!dt = 1800._8
!dt = real(dtlong,8)

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   latdeg = glatw(iw)
   londeg = glonw(iw)

   ka = lpw(iw)

   do k = 2,mza

      ! Update chemistry tendencies only at specified time interval

      if (mod(time8p,dt) < dtlong) then
         cl  = real(addsc(1)%sclp(k,iw),8)
         cl2 = real(addsc(2)%sclp(k,iw),8)

         call tendency_Terminator(latdeg, londeg, cl, cl2, dt, cl_f, cl2_f)

         ctend1(k,iw) = cl_f
         ctend2(k,iw) = cl2_f
      endif

      ! Copy chemistry tendencies to main model arrays on each subroutine call

      if (nl%naddsc >= 2) then
         addsc(1)%sclt(k,iw) = addsc(1)%sclt(k,iw) + real(ctend1(k,iw) * rho(k,iw))
         addsc(2)%sclt(k,iw) = addsc(2)%sclt(k,iw) + real(ctend2(k,iw) * rho(k,iw))
      endif

   enddo

enddo

end subroutine olam_dcmip_terminator
