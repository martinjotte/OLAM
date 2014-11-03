subroutine olam_dcmip_init()

use mem_basic,  only: vc, vmc, vp, vmp, wc, wmc, rho, thil, theta, tair, &
                      press, sh_w, sh_v
use mem_ijtabs, only: itab_v, itab_w, jtab_v, jtab_w
use mem_grid,   only: mza, mva, mwa, zt, dzt, xev, yev, zev, vnx, vny, vnz, lpw, &
                      glatw, glonw, lpv, zm, glatv, glonv 
use mem_addsc,  only: addsc
use misc_coms,   only: io6, iparallel
use consts_coms, only: p00, rocp, erad, pio180, cvocp, p00k, rdry, rvap, &
                       alvlocp, gravo2, xscale
use mem_micro,   only: sh_c
use micro_coms,  only: level

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, &
   test1_advection_hadley, &
   test1_advection_orography, &
   test2_steady_state_mountain, &
   test2_schaer_mountain, &
   test3_gravity_wave

use test2xdeep, only: deep_mountain_wave

use dcmip_initial_conditions_test_4, only: &
   test4_baroclinic_wave

use dcmip_initial_conditions_test_5, only: &
   test5_tropical_cyclone

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                        mpi_send_v, mpi_recv_v 

use test41alt, only: baroclinic_instability_alt

use oname_coms, only: nl

implicit none

real(8) :: time0

real(8) :: lon,lat,p,t,phis,ps,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr
real(8) :: hyam, hybm, gc, xscal

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0, rm0, rt0

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,rcoords,cfv,shear,moist,iter,ka, &
           subtest, deep

logical :: hybrid_eta

real :: exner,temp,rhovs

real, external :: rhovsl

zcoords = 1
rcoords = 1
time0   = 0.0d0
xscal  = xscale

hybrid_eta = .false.

cfv = 0
shear = 0

mrl = 1

!----------------------------------------------------------------------
!do j = 1,jtab_w(14)%jend(1); iw = jtab_w(14)%iw(j)
do iw = 2,mwa  ! do for all points in subdomain
!----------------------------------------------------------------------
   lon = pio180 * glonw(iw)
   lat = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lon = pio180 * (glonw(iw) + 360.)

   do k = lpw(iw),mza
      p = press(k,iw)
      zm0 = zm(k)
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

         call test1_advection_deformation(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2,q3,q4,time0)

         call test1_advection_deformation(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2
         if (nl%naddsc >= 3) addsc(3)%sclp(k,iw) = q3
         if (nl%naddsc >= 4) addsc(4)%sclp(k,iw) = q4
         if (nl%naddsc >= 5) addsc(5)%sclp(k,iw) = 1.0d0

      elseif (nl%test_case == 12) then

!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

         call test1_advection_hadley(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,time0)

         call test1_advection_hadley(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1

      elseif (nl%test_case == 13) then

!==========================================================================================
! TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS IN THE PRESENCE OF OROGRAPHY
!==========================================================================================

         call test1_advection_orography(lon,lat,p,zm0,zcoords, &
            cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2,q3,q4)

         call test1_advection_orography(lon,lat,p,zt0,zcoords, &
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

!=========================================================================
! Test 2-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================

         call test2_steady_state_mountain(lon,lat,p,zm0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wm0, &
            t,phis,ps,rhom0,q)

         call test2_steady_state_mountain(lon,lat,p,zt0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 21 .or. &
              nl%test_case == 22) then

!=====================================================================================
! Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over a Schaer-type Mountain
!=====================================================================================

         if (nl%test_case == 21) then
            shear = 0
            subtest = 1
         else  ! nl%test_case == 22
            shear = 1
            subtest = 2
         endif

         rm0 = 6371220. / xscale + zm(k)
         rt0 = 6371220. / xscale + zt(k)

         call deep_mountain_wave(subtest,Xscal,lon,lat,p,rm0,rcoords, &
                                 u0,v0,wm0,t,phis,ps,rhom0)

         call deep_mountain_wave(subtest,Xscal,lon,lat,p,rt0,rcoords, &
                                 u0,v0,wt0,t,phis,ps,rhot0)

         q = 0.0d0

! The following 2 calls (commented out) are for shallow atmosphere

!         call test2_schaer_mountain(lon,lat,p,zm0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wm0, &
!            t,phis,ps,rhom0,q)

!         call test2_schaer_mountain(lon,lat,p,zt0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wt0, &
!            t,phis,ps,rhot0,q)

     elseif (nl%test_case == 31) then

!==========================
! Test 3-1 - GRAVITY WAVES
!==========================

         call test3_gravity_wave(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q)

         call test3_gravity_wave(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 410 .or. &
              nl%test_case == 411 .or. &
              nl%test_case == 412 .or. &
              nl%test_case == 413) then

!==========================
! Test 41x - DRY BAROCLINIC INSTABILITY
!==========================

! New code for deep atmosphere case:

         deep = 1

         call baroclinic_instability_alt(deep,Xscal,lon,lat,p,zm0,zcoords, &
                                         u0,v0,wm0, &
                                         t,phis,ps,rhom0,q1,q2)

         call baroclinic_instability_alt(deep,Xscal,lon,lat,p,zt0,zcoords, &
                                         u0,v0,wt0, &
                                         t,phis,ps,rhot0,q1,q2)

         q = 0.0d0

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2

      elseif (nl%test_case == 42 .or. &
              nl%test_case == 43) then

!==========================
! Tests 42,43 - MOIST BAROCLINIC INSTABILITY
!==========================

         moist = 1

         call test4_baroclinic_wave(moist,Xscal,lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2)
 
         call test4_baroclinic_wave(moist,Xscal,lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2)

         if (nl%naddsc >= 1) addsc(1)%sclp(k,iw) = q1
         if (nl%naddsc >= 2) addsc(2)%sclp(k,iw) = q2

      elseif (nl%test_case == 51 .or. &
              nl%test_case == 52) then

!==========================================================================================
! TEST CASE 5 - Tropical Cyclone 
!==========================================================================================

         call test5_tropical_cyclone(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q)

         call test5_tropical_cyclone(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      endif

      wc(k,iw) = wm0
      wmc(k,iw) = wm0 * rhom0

      rho(k,iw) = rhot0
      press(k,iw) = p
      tair (k,iw) = t
      theta(k,iw) = t * (p00 / p) ** rocp
      thil(k,iw) = theta(k,iw)
      sh_w(k,iw) = q
      sh_v(k,iw) = q

   enddo

! Iterative hydrostatic balance procedure
! (Is this necessary for any test cases???)
! Procedure holds THETA constant, after THETA is computed above from t and p

go to 100

   do iter = 1,100
   
      do k = 1,mza
      
         if (level == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (level == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) / p00) ** rocp   ! Defined WITHOUT CP factor
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)
            sh_c(k,iw) = max(0.,sh_w(k,iw)-rhovs/real(rho(k,iw)))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)
! As for iteration in subroutine satadjst, use (0.3,0.7) weighting to damp iteration
!            rho(k,iw) = .5 * rho(k,iw)  &
!                      + .5 * press(k,iw) ** cvocp * p00k  &
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (1. - sh_c(k,iw))  &
               * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
            thil(k,iw) = theta(k,iw)  &
               / (1. + alvlocp * sh_c(k,iw) / max(temp,253.))
         endif
         
         if (k >= 2) then
! Impose minimum value of 1 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations
             press(k,iw) = .05d0 * press(k,iw)  &
                         + .95d0 * max(1.d0,press(k-1,iw)  &
                - gravo2 * (rho(k-1,iw) * dzt(k-1) + rho(k,iw) * dzt(k)))
         endif
         
      enddo
   enddo

   do k = 1, mza
      tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
   enddo

100 continue

enddo

! MPI parallel send/recv of W group

if (iparallel == 1) then
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil))

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil))
endif

!----------------------------------------------------------------------
do iv = 2,mva
!----------------------------------------------------------------------
   lon = pio180 * glonv(iv)
   lat = pio180 * glatv(iv)

   if (glonv(iv) < 0.) lon = pio180 * (glonv(iv) + 360.)

   do k = lpv(iv),mza
      p = 1.0d5
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

         call test1_advection_deformation(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

      elseif (nl%test_case == 12) then

!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

         call test1_advection_hadley(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

      elseif (nl%test_case == 13) then

!==========================================================================================
! TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS IN THE PRESENCE OF OROGRAPHY
!==========================================================================================

         call test1_advection_orography(lon,lat,p,zt0,zcoords, &
            cfv,hybrid_eta,hyam,hybm,gc,u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4)

! Alteration for OLAM: Set velocity at mountain heights = 0.

            if (zt0 < 2000.d0) then
               u0 = 0.0d0
               v0 = 0.0d0
            endif

      elseif (nl%test_case == 200 .or. &
              nl%test_case == 201) then

!=========================================================================
! Test 2-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================

         call test2_steady_state_mountain(lon,lat,p,zt0,zcoords, &
            hybrid_eta,hyam,hybm,u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 21 .or. &
              nl%test_case == 22) then

!=====================================================================================
! Tests 2-1 and 2-2:  Non-hydrostatic Mountain Waves over a Schaer-type Mountain
!=====================================================================================

         if (nl%test_case == 21) then
            shear = 0
            subtest = 1
         else  ! nl%test_case == 22
            shear = 1
            subtest = 2
         endif

         rt0 = 6371220. / xscale + zt(k)

         call deep_mountain_wave(subtest,Xscal,lon,lat,p,rt0,rcoords, &
                                 u0,v0,wt0,t,phis,ps,rhot0)

if (iv == 1000) print*, 'init22 ',u0,v0

! The following call (commented out) is for shallow atmosphere approx.

!         call test2_schaer_mountain(lon,lat,p,zt0,zcoords, &
!            hybrid_eta,hyam,hybm,shear,u0,v0,wt0, &
!            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 31) then

!==========================
! Test 3-1 - GRAVITY WAVES
!==========================

         call test3_gravity_wave(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

      elseif (nl%test_case == 410 .or. &
              nl%test_case == 411 .or. &
              nl%test_case == 412 .or. &
              nl%test_case == 413) then

!==========================
! Test 41x - DRY BAROCLINIC INSTABILITY
!==========================

! New code for deep atmosphere case:

         deep = 1

         call baroclinic_instability_alt(deep,Xscal,lon,lat,p,zt0,zcoords, &
                                         u0,v0,wt0, &
                                         t,phis,ps,rhot0,q1,q2)

      elseif (nl%test_case == 42 .or. &
              nl%test_case == 43) then

!==========================
! Tests 42,43 - MOIST BAROCLINIC INSTABILITY
!==========================

         call test4_baroclinic_wave(moist,Xscal,lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2)

      elseif (nl%test_case == 51 .or. &
              nl%test_case == 52) then

!==========================================================================================
! TEST CASE 5 - Tropical Cyclone 
!==========================================================================================

         call test5_tropical_cyclone(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q)

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
enddo

! For below-ground points, set VC to LCV value.

do iv = 2,mva
   ka = lpv(iv)
   vc(1:ka-1,iv) = vc(ka,iv)
enddo

! MPI parallel send/recv of V group

if (iparallel == 1) then
   call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
endif

! Set VMP and VP

if (allocated(vmp)) vmp(:,:) = vmc(:,:)
if (allocated(vp )) vp (:,:) = vc (:,:)

return
end subroutine olam_dcmip_init

!==========================================================================

subroutine olam_dcmip_vel(time0,vmsc,wmsc,rho_old)

use mem_basic,  only: vc, vmc, wc, wmc, rho, thil, theta, press, sh_w, sh_v, &
                      vxe, vye, vze
use mem_ijtabs, only: itab_v, itab_w, jtab_v, jtab_w
use mem_grid,   only: mza, mva, mwa, zt, xev, yev, zev, vnx, vny, vnz, lpw, &
                      glatw, glonw, lpv, zm, glatv, glonv, xew, yew, zew, dzt
use mem_addsc,  only: addsc
use consts_coms, only: p00, rocp, erad, pio180, p00i, grav
use misc_coms,  only: dtlong
use mem_micro,  only: pcpgr

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, &
   test1_advection_hadley

use oname_coms, only: nl

implicit none

real(8), intent(in) :: time0

real, intent(out) :: vmsc(mza,mva),wmsc(mza,mwa)

real(8), intent(out) :: rho_old(mza,mwa) ! density at beginning of timestep [kg/m^3]

real(8) :: lon,lat,p,z,t,phis,ps,rho0,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,test

zcoords = 1.0d0

!----------------------------------------------------------------------
do iw = 2,mwa  ! do for all points in subdomain
!----------------------------------------------------------------------
   lon = pio180 * glonw(iw)
   lat = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lon = pio180 * (glonw(iw) + 360.)

   do k = lpw(iw),mza
      p = press(k,iw)
      zm0 = zm(k)
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

         call test1_advection_deformation(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,q2,q3,q4,time0)

         call test1_advection_deformation(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

         wc(k,iw) = wm0
         wmc(k,iw) = wm0 * rhom0
         rho_old(k,iw) = rho(k,iw)

      elseif (nl%test_case == 12) then

!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

         call test1_advection_hadley(lon,lat,p,zm0,zcoords, &
            u0,v0,wm0, &
            t,phis,ps,rhom0,q,q1,time0)

         call test1_advection_hadley(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

         wc(k,iw) = wm0
         wmc(k,iw) = wm0 * rhom0
         rho_old(k,iw) = rho(k,iw)

      endif

      wmsc(k,iw) = wmc(k,iw)

   enddo

enddo

!----------------------------------------------------------------------
do iv = 2,mva
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   lon = pio180 * glonv(iv)
   lat = pio180 * glatv(iv)

   if (glonv(iv) < 0.) lon = pio180 * (glonv(iv) + 360.)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

   do k = lpv(iv),mza
      p = 1.0d5
      zt0 = zt(k)

      if (nl%test_case == 11) then

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

         call test1_advection_deformation(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,q2,q3,q4,time0)

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

      elseif (nl%test_case == 12) then

!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

         call test1_advection_hadley(lon,lat,p,zt0,zcoords, &
            u0,v0,wt0, &
            t,phis,ps,rhot0,q,q1,time0)

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

      endif

      vmsc(k,iv) = vmc(k,iv)

   enddo
enddo

return
end subroutine olam_dcmip_vel

!==========================================================================

subroutine olam_dcmip_simplephys(rhot)

use mem_basic,  only: vc, vmc, wc, wmc, rho, thil, theta, press, sh_w, sh_v, &
                      vxe, vye, vze
use mem_ijtabs, only: itab_v, itab_w, jtab_v, jtab_w
use mem_grid,   only: mza, mva, mwa, zt, xev, yev, zev, vnx, vny, vnz, lpw, &
                      glatw, glonw, lpv, zm, glatv, glonv, xew, yew, zew, dzt
use mem_addsc,  only: addsc
use consts_coms, only: p00, rocp, erad, pio180, p00i, grav
use misc_coms,  only: dtlong
use mem_micro,  only: pcpgr

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, &
   test1_advection_hadley

use oname_coms, only: nl

implicit none

real, intent(inout) :: rhot(mza,mwa)

real(8) :: lon,lat,p,z,t,phis,ps,rho0,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0, dtime

real(8) :: tcol(mza), qcol(mza), ucol(mza), vcol(mza), u3d(mza,mwa), v3d(mza,mwa)
real(8) :: pmidcol(mza), pintcol(mza), pdelcol(mza), rpdelcol(mza)
real(8) :: precl, vx, vy, vz, exner, za, dtheta

real :: thlprt(mza),tprt(mza),qprt(mza),uprt(mza),vprt(mza)

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,test,pver

zcoords = 1.0d0
dtime = dtlong

 pver = mza - 1

 call olam_dcmip_simplephys0(pver,rhot)

 RETURN

!----------------------------------------------------------------------
do iw = 2,mwa  ! do for all points in subdomain
!----------------------------------------------------------------------
   lon = pio180 * glonw(iw)
   lat = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lon = pio180 * (glonw(iw) + 360.)

   if (nl%test_case == 51) test = 0
   if (nl%test_case == 43) test = 1
   if (nl%test_case == 42) test = 2

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   do k = 2,mza

      vx = vxe(k,iw)
      vy = vye(k,iw)
      vz = vze(k,iw)

      if (raxis > 1.e3) then
         u0 = (vy * xew(iw) - vx * yew(iw)) / raxis
         v0 = vz * raxis / erad &
           - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 
      else
         u0 = 0.
         v0 = 0.
      endif

      exner = (p00i * press(k,iw)) ** rocp

      tcol(k) = theta(k,iw) * exner
      qcol(k) = sh_w(k,iw)
      ucol(k) = u0
      vcol(k) = v0
      pmidcol(k) = press(k,iw)
      pintcol(k) = press(k,iw) - rho(k,iw) * grav * (zm(k) - zt(k))
      pdelcol(k) = rho(k,iw) * grav * dzt(k)
      rpdelcol(k) = 1. / pdelcol(k)

   enddo

   pintcol(1) = press(2,iw) + rho(2,iw) * grav * (zt(2) - zm(1))
   ps = pintcol(1)
   za = zt(2)

   call simple_physics (mza, 2, za, dtime, lat, tcol, qcol, ucol, vcol, &
      pmidcol, pintcol, pdelcol, rpdelcol, ps, precl, test, iw)

! Apply changes to model prognostic fields

   do k = 2,mza

      exner = (p00i * press(k,iw)) ** rocp

      dtheta = tcol(k) / exner - theta(k,iw)
      thil(k,iw) = thil(k,iw) + dtheta
      rho(k,iw) = rho(k,iw) + rho(k,iw) * (qcol(k) - sh_w(k,iw))
      sh_w(k,iw) = qcol(k)
      u3d(k,iw) = ucol(k)
      v3d(k,iw) = vcol(k)

      ! Should we update T, THETA here too??
      ! tair (k,iw) = tcol(k)
      ! theta(k,iw) = tcol(k) / exner

   enddo

   pcpgr(iw) = precl

enddo

!----------------------------------------------------------------------
do iv = 2,mva
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   lon = pio180 * glonv(iv)
   lat = pio180 * glatv(iv)

   if (glonv(iv) < 0.) lon = pio180 * (glonv(iv) + 360.)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

   do k = lpv(iv),mza

      u0 = .5 * (u3d(k,iw1) + u3d(k,iw2))
      v0 = .5 * (v3d(k,iw1) + v3d(k,iw2))

      if (raxis > 1.e3) then
         uv01dr = -v0 * zev(iv) / erad  ! radially outward from axis

         uv01dx = (-u0 * yev(iv) + uv01dr * xev(iv)) / raxis 
         uv01dy = ( u0 * xev(iv) + uv01dr * yev(iv)) / raxis 
         uv01dz =   v0 * raxis / erad 

         vc(k,iv) = uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv)
      else
         vc(k,iv) = 0.
      endif

      vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

   enddo
enddo

return
end subroutine olam_dcmip_simplephys

!==========================================================================

subroutine olam_dcmip_simplephys0(pver,rhot)

use mem_basic,  only: vc, vmc, wc, wmc, rho, thil, theta, press, sh_w, sh_v, &
                      vxe, vye, vze
use mem_ijtabs, only: itab_v, itab_w, jtab_v, jtab_w
use mem_grid,   only: mza, mva, mwa, zt, xev, yev, zev, vnx, vny, vnz, lpw, &
                      glatw, glonw, lpv, zm, glatv, glonv, xew, yew, zew, dzt
use mem_addsc,  only: addsc
use consts_coms, only: p00, rocp, erad, pio180, p00i, grav
use misc_coms,  only: dtlong
use mem_micro,  only: pcpgr

use mem_tend,   only: thilt, sh_wt, vmt

use dcmip_initial_conditions_test_1_2_3, only: &
   test1_advection_deformation, &
   test1_advection_hadley

use oname_coms, only: nl

implicit none

integer, intent(in) :: pver

real, intent(inout) :: rhot(mza,mwa)

real(8) :: lon,lat,p,z,t,phis,ps,rho0,q,q1,q2,q3,q4
real(8) :: raxis, u01d, v01d, uv01dx, uv01dy, uv01dz, uv01dr

real(8) :: zm0, zt0, rhom0, rhot0, u0, v0, wm0, wt0, dtime

real(8) :: tcol(pver), qcol(pver), ucol(pver), vcol(pver), dudt3d(pver,mwa), dvdt3d(pver,mwa)
real(8) :: pmidcol(pver), pintcol(pver), pdelcol(pver), rpdelcol(pver)
real(8) :: precl, vx, vy, vz, exner, za, dtheta, dudt0, dvdt0

! Physics Tendency Arrays
real(8) :: dtdt(pver)             ! Temperature tendency 
real(8) :: dqdt(pver)             ! Specific humidity tendency
real(8) :: dudt(pver)             ! Zonal wind tendency
real(8) :: dvdt(pver)             ! Meridional wind tendency

real :: thlprt(pver),tprt(pver),qprt(pver),uprt(pver),vprt(pver)

integer :: mrl,j,iw,k,iv,iw1,iw2,zcoords,test,kver

zcoords = 1.0d0
dtime = dtlong

!----------------------------------------------------------------------
do iw = 2,mwa  ! do for all points in subdomain
!----------------------------------------------------------------------
   lon = pio180 * glonw(iw)
   lat = pio180 * glatw(iw)

   if (glonw(iw) < 0.) lon = pio180 * (glonw(iw) + 360.)

   if (nl%test_case == 51) test = 0
   if (nl%test_case == 43) test = 1
   if (nl%test_case == 42) test = 2

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

   do k = 2,mza

      kver = mza + 1 - k

      vx = vxe(k,iw)
      vy = vye(k,iw)
      vz = vze(k,iw)

      if (raxis > 1.e3) then
         u0 = (vy * xew(iw) - vx * yew(iw)) / raxis
         v0 = vz * raxis / erad &
           - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 
      else
         u0 = 0.
         v0 = 0.
      endif

      exner = (p00i * press(k,iw)) ** rocp

      tcol(kver) = theta(k,iw) * exner
      qcol(kver) = sh_w(k,iw)
      ucol(kver) = u0
      vcol(kver) = v0
      pmidcol(kver) = press(k,iw)
      pintcol(kver) = press(k,iw) - rho(k,iw) * grav * (zm(k) - zt(k))
      pdelcol(kver) = rho(k,iw) * grav * dzt(k)
      rpdelcol(kver) = 1. / pdelcol(kver)

   enddo

   pintcol(pver+1) = press(2,iw) + rho(2,iw) * grav * (zt(2) - zm(1))
   ps = pintcol(pver+1)
   za = zt(2)

   call simple_physics0 (pver, dtime, lat, tcol, qcol, ucol, vcol, &
      pmidcol, pintcol, pdelcol, rpdelcol, ps, precl, test, iw, &
      dtdt, dqdt, dudt, dvdt)

! Apply changes to model prognostic fields

   do k = 2,mza

      kver = mza + 1 - k

      exner = (p00i * press(k,iw)) ** rocp

      thilt(k,iw) = thilt(k,iw) + rho(k,iw) * dtdt(kver) / exner
      sh_wt(k,iw) = rho(k,iw) * dqdt(kver)
      rhot(k,iw) = rhot(k,iw) + rho(k,iw) * dqdt(kver)

      rho(k,iw) = rho(k,iw) + rho(k,iw) * (qcol(kver) - sh_w(k,iw))
      sh_w(k,iw) = qcol(kver)

      if (nl%test_case > 42) then
         dudt3d(k,iw) = dudt(kver)
         dvdt3d(k,iw) = dvdt(kver)
      endif

   enddo

   pcpgr(iw) = precl

enddo

if (nl%test_case == 42) return

!----------------------------------------------------------------------
do iv = 2,mva
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   lon = pio180 * glonv(iv)
   lat = pio180 * glatv(iv)

   if (glonv(iv) < 0.) lon = pio180 * (glonv(iv) + 360.)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

   do k = lpv(iv),mza

      rho0  = .5 * (rho(k,iw1) + rho(k,iw2))

      dudt0 = .5 * (dudt3d(k,iw1) + dudt3d(k,iw2))
      dvdt0 = .5 * (dvdt3d(k,iw1) + dvdt3d(k,iw2))

      if (raxis > 1.e3) then
         uv01dr = -dvdt0 * zev(iv) / erad  ! radially outward from axis

         uv01dx = (-dudt0 * yev(iv) + uv01dr * xev(iv)) / raxis 
         uv01dy = ( dudt0 * xev(iv) + uv01dr * yev(iv)) / raxis 
         uv01dz =   dvdt0 * raxis / erad 

         vmt(k,iv) = vmt(k,iv) * rho0 * (uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv))
      endif

   enddo
enddo

return
end subroutine olam_dcmip_simplephys0

