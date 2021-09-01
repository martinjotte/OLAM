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
module hcane_rz

implicit none

! The purpose of module hcane_rz is to restore the intensity and eyewall
! diameter of a hurricane that is under-resolved in the initial conditions of
! an OLAM simulation.  This situation commonly arises when the only gridded
! initialization dataset is of low resolution (for a hurricane), such as the
! CFSR global gridded reanalysis.  The procedure consists of one or more
! dynamic initialization cycles in which (1) OLAM is first initialized by
! interpolation from the CFS (or other) gridded dataset, (2) an axisymmetric 
! tangential perturbation wind field is added to the under-resolved vortex such
! that the sum is close to in-situ observations of tangential winds and the
! eyewall diameter, (3) OLAM adds a compensating potential temperature
! perturbation field that approximately maintains hydrostatic and gradient wind
! balance with the strengthened tangential winds, (4) OLAM is integrated
! forward in time one or more hours to allow the hurricane to more completely
! adjust to the added perturbation fields and to develop internal dynamic, 
! thermodynamic, and moisture structures, and (5) the simulated hurricane is
! remapped from its final simulated location back to the grid cells at its
! initial location in preparation for the next cycle.
!
! The number of cycles performed, the temporal length of each cycle, and the
! target eyewall radius and tangential wind speed are all user-specified
! parameters.  Plots of the vortex are made at the beginning and end of
! each cycle to allow inspection of the progress and quality of the
! initialization.

  ! radius_ax is an array of radial distances of the radial-height grid on which
  ! axisymmetric hurricane fields are diagnosed and modified by perturbations
  ! defined by the user.  The radial spacing between these values needs to be at
  ! least twice the grid spacing of the OLAM hexagonal grid where the hurricane
  ! is located.

  integer, parameter :: nz = 60, nr = 28 ! Number of vertical, radial points in vortex profiles

  !real :: radius_ax(nr) = (/ &
  !   0.e3,  10.e3,  20.e3,  30.e3,  40.e3,  50.e3,  60.e3,  70.e3,  80.e3,  90.e3, &
  ! 100.e3, 120.e3, 140.e3, 160.e3, 180.e3, 200.e3, 250.e3, 300.e3, 350.e3, 400.e3 /)

  real :: radius_ax(nr) = (/ &
     0.e3,   5.e3,  10.e3,  15.e3,  20.e3,  25.e3,  30.e3,  35.e3,  40.e3,  45.e3, &
    50.e3,  55.e3,  60.e3,  70.e3,  80.e3,  90.e3, 100.e3, 120.e3, 140.e3, 160.e3, &
   180.e3, 200.e3, 240.e3, 280.e3, 320.e3, 360.e3, 400.e3, 440.e3 /)

  integer :: nzz

  integer :: ncycle_hurrinit
  integer :: icycle_hurrinit

  real :: hlat0, hlon0, hlat, hlon
  real :: hlat_hist = 0., hlon_hist = 0.

  real :: circ_avg(nr)

  real :: vtan_eyw    ! Target maximum tangential wind speed 
  real :: rad_eyw     ! radial index (of radius_ax array) where vtan_eyw applies
  real :: rad_env     ! radial index (of radius_ax array) where perturbation ceases

  real :: rexpon_delc   ! Radial power law of perturbation magnitude
  real :: zmax_delv     ! Maximum height of perturbation vortex
  real :: zexpon_delv   ! Vertical power law of perturbation magnitude

  ! The following 4 parameters control how a relocated vortex is blended in with
  ! model initial fields at the beginning of the second and subsequent dynamic
  ! initialization cycles.  At locations inside radius rad1 and at heights below
  ! z1, 100% of the relocated fields are used (implying 0% of the model initial
  ! fields).  At locations outside radius rad2 and/or at heights above z2, none of
  ! the relocated vortex is used, and the model initial fields remain unchanged.
  ! Between these limits, interpolation weights vary linearly with height
  ! and/or radius.  LIMIT RAD2_BLEND TO ABOUT [radius_ax(nr) - 20.e3] OR LESS.

  real :: rad1_blend    ! Inner radius for relocation blending weights
  real :: rad2_blend    ! Outer radius for relocation blending weights
  real :: z1_blend      ! Lower height for relocation blending weights
  real :: z2_blend      ! Upper height for relocation blending weights

  integer :: ir_eyw, ir_env

  ! Axisymmmetric vortex profile arrays

  real ::  thil_ax(nz,nr) ! ice-liquid potential temperature (K)
  real :: theta_ax(nz,nr) ! potential temperature (K)
  real ::  tair_ax(nz,nr) ! temperature (K)
  real ::   rrw_ax(nz,nr) ! total water mixing ratio (kg/kg_dryair)
  real ::   rrv_ax(nz,nr) ! vapor mixing ratio (kg/kg_dryair)
  real ::  vtan_ax(nz,nr) ! tangential wind (m/s)
  real ::  vrad_ax(nz,nr) ! radial wind (m/s)
  real ::     w_ax(nz,nr) ! vertical wind (m/s)
  real :: ssliq_ax(nz,nr) ! supsat wrt liquid (%)
  real :: thdif_ax(nz,nr) ! theta difference from supsat (K)

  ! Relocation fields for I/O

  integer, parameter :: nfld = 28    ! number of fields transferred

  integer :: nout_dim ! Horiz # of W points filled by interpolation (array dim)
  integer :: nout     ! Horiz # of W points filled by interpolation (actual count)

  integer, parameter :: nips = 200, njps = 200 ! size of PS arrays

  integer :: nwps(nips,njps)    ! # of iw indices in each ips,jps cell
  integer :: iwps(nips,njps,60) ! iw indices in each ips,jps cell (max of 60)

  integer, allocatable :: iwout(:)

  real, allocatable :: reloc_field(:,:,:)

  real :: reh0, xeh0, yeh0, zeh0

Contains

!===============================================================================

  subroutine hurricane_init()

  use mem_grid,    only: mza, mwa, zt
  use consts_coms, only: pio180, erad

  implicit none

  integer :: k, ir

  ! This subroutine is called once to initialize some quantities

  ! Count number of model levels that are below the chosen maximum height
  ! for the dynamic initialization procedure (e.g. 20 km).

  do k = 2,mza
     nzz = k
     if (zt(k) > 20000.) exit
  enddo

  ir_eyw = 2
  ir_env = 3

  do ir = 1,nr
     if (abs(radius_ax(ir) - rad_eyw) < abs(radius_ax(ir_eyw) - rad_eyw)) &
        ir_eyw = ir
     if (abs(radius_ax(ir) - rad_env) < abs(radius_ax(ir_env) - rad_env)) &
        ir_env = ir
  enddo

  print*, 'ir_eyw, ir_env ',ir_eyw,ir_env

  ! Find "earth" coordinates of hurricane center initial location

  reh0 = erad * cos(hlat0 * pio180)  ! distance from earth axis
  zeh0 = erad * sin(hlat0 * pio180)
  xeh0 = reh0 * cos(hlon0 * pio180)
  yeh0 = reh0 * sin(hlon0 * pio180)

  end subroutine hurricane_init

!==================================================================================

  subroutine vortex_center_diagnose()

  ! This subroutine diagnoses the latitude & longitude of the MSL pressure minimum 
  ! inside a tropical cyclone.  It requires a nearby starting location (hlat,hlon)
  ! to be provided, either from a user specification of the known cyclone location
  ! or from a recent location previously diagnosed from this subroutine.

  ! The algorithm of subroutine vortex_center_diagnose is not to simply search for
  ! the lowest MSL pressure in the vicinity, but rather to take a weighted average 
  ! over multiple grid cells having pressure close to that lowest value.  This yields
  ! a location that is more steady in time, less subject to turbulent fluctuations,
  ! and more tied to the circulation at the scale of the cyclone core region.

  ! NOTE: THIS SUBROUTINE IS NOT MPI-COMPATIBLE; IT ASSUMES THAT ALL POINTS IN THE
  ! INNER REGION OF THE CYLONE ARE CONTAINED ON ONE AND THE SAME COMPUTER PROCESS.  

  use mem_ijtabs, only: 
  use mem_basic, only: press, tair
  use misc_coms
  use mem_grid, only: mwa, zt, lpw, xew, yew, zew, arw0, glatw, glonw
  use consts_coms, only: pio180, erad

  implicit none

  integer :: iw, ka

  real :: reh, xeh, yeh, zeh
  real :: dist, weight
  real :: rlon, rlat

  real :: area_tot, pmsl_min, pmsl_avg, pmsl_thresh
  real :: xew_avg, yew_avg, zew_avg, weight_sum
  real :: pmsl(mwa)

  write(6,'(a,2f10.3)') 'vortex_center_diagnose BEGIN ',hlat,hlon

  ! Find "earth" coordinates of first-guess position (hlat,hlon)

  reh = erad * cos(hlat * pio180)  ! distance from earth center
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)
  zeh = erad * sin(hlat * pio180)

  ! Initialize quantities

  area_tot = 0.
  pmsl_min = 2.e5
  pmsl_avg = 0.

  xew_avg = 0.
  yew_avg = 0.
  zew_avg = 0.

  weight_sum = 0.

  ! Horizontal loop over all W points

  do iw = 2,mwa

     ! Skip current W point if its location is far from first guess position
     ! (This is rough check that eliminates most points)

     if (abs(glatw(iw) - hlat) > 5.) cycle
     if (abs(glonw(iw) - hlon) > 5.) cycle

     ! Distance of current W point to first guess position

     dist = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

     ! Skip current W point if its location is far from first guess position
     ! (This is finer check)

     if (dist > 100.e3) cycle

     ka = lpw(iw)

     pmsl(iw) = press(ka,iw) &
              * (1. - .0065 * zt(ka) / (tair(ka,iw) + .0065 * zt(ka)))**(-5.257)

     ! Determine lowest value of pmsl within search area

     if (pmsl_min > pmsl(iw)) then
         pmsl_min = pmsl(iw)
     endif

     ! Sum grid cell area and (area * pressure) product at each model level within search area

     area_tot = area_tot + arw0(iw)
     pmsl_avg = pmsl_avg + arw0(iw) * pmsl(iw)

  enddo  ! iw

  ! Compute average pressure for any k level with cells above ground
  ! Compute threshold pressure at 80% of the range from avg to min

  pmsl_avg = pmsl_avg / area_tot
  pmsl_thresh = pmsl_avg + .80 * (pmsl_min - pmsl_avg)

  ! Horizontal loop over all W points

  do iw = 2,mwa

     ! Skip current W point if its location is far from first guess position
     ! (This is rough check that eliminates most points)

     if (abs(glatw(iw) - hlat) > 5.) cycle
     if (abs(glonw(iw) - hlon) > 5.) cycle

     ! Distance of current W point to first guess position

     dist = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

     ! Skip current W point if its location is far from first guess position
     ! (This is finer check)

     if (dist > 100.e3) cycle

     ! If pmsl < threshold value, sum area-pressure weighted grid cell location

     if (pmsl(iw) < pmsl_thresh) then

        weight = arw0(iw) * (pmsl_thresh - pmsl(iw)) ! ** exponent  (if desired)

        xew_avg = xew_avg + weight * xew(iw)
        yew_avg = yew_avg + weight * yew(iw)
        zew_avg = zew_avg + weight * zew(iw)

        weight_sum = weight_sum + weight

     endif

  enddo  ! iw

  ! Compute mean location

  xew_avg = xew_avg / weight_sum
  yew_avg = yew_avg / weight_sum
  zew_avg = zew_avg / weight_sum

  ! Transform mean location to lat/lon coordinates

  call e_ec(xew_avg,yew_avg,zew_avg,rlon,rlat)

  hlat = rlat
  hlon = rlon

  hlat_hist = rlat
  hlon_hist = rlon

  write(6,'(a,2f10.3,f12.2)') &
     'vortex_center_diagnose END: hlat, hlon, pmsl_min   ',hlat, hlon, pmsl_min

  end subroutine vortex_center_diagnose

!==================================================================================

  subroutine vortex_azim_avg()

  ! This subroutine computes axisymmetric azimuthal averages, centered around the
  ! latitude & longitude of the minimum PMSL pressure in a tropical cyclone, of
  ! cyclone radial, azimuthal, and vertical wind components and a few scalar
  ! fields as a function of radius and height.  Azimuthal averages are computed
  ! on each model level up to a specified maximum height and over intervals
  ! between radial distance values defined in the radius_ax array.

  use mem_ijtabs,  only: jtw_init, jtab_w
  use mem_basic,   only: thil, theta, tair, rho, rr_w, rr_v, wc, vxe, vye, vze
  use mem_grid,    only: xew, yew, zew, lpw
  use consts_coms, only: erad, pio180, alvlocp
  use therm_lib,   only: rhovsl

  implicit none

  integer :: iw, j, k, irad
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan

  real :: rad,wrad1,wrad2
  real :: vtan_ax0, vrad_ax0, ss_liq, ss_thetadif

  real :: weight_t(nz,nr)   ! weight array for T points

  ! Find "earth" coordinates of hurricane center

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

  ! Initialize axixymmetric arrays to zero prior to summation

   thil_ax(:,:) = 0.
  theta_ax(:,:) = 0.
   tair_ax(:,:) = 0.
    rrw_ax(:,:) = 0.
    rrv_ax(:,:) = 0.
   vtan_ax(:,:) = 0.
   vrad_ax(:,:) = 0.
      w_ax(:,:) = 0.
  ssliq_ax(:,:) = 0.
  thdif_ax(:,:) = 0.

  weight_t(:,:) = 0.

  ! Loop over all W points

  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)  ! jend(1) = hardwired for mrl 1

     ! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

     ! Skip axisymmetric averaging for all points outside specified radius

     if (rad >= radius_ax(nr) - 1.) cycle

     ! Determine interpolation point in radial dimension

     irad = 1
     do while (rad > radius_ax(irad+1))
        irad = irad + 1
     enddo

     if (irad < 1 .or. irad+1 > nr) then
        print*, 'irad out of bounds1 ', irad, nr
     endif

     wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
     wrad1 = 1. - wrad2

     ! Unit normal vector components from hurricane center to current IW point

     wnxrad = (xew(iw) - xeh) / rad
     wnyrad = (yew(iw) - yeh) / rad
     wnzrad = (zew(iw) - zeh) / rad

     ! Unit vector components in direction of tangential vortex wind

     wnxtan = wnyh * wnzrad - wnzh * wnyrad
     wnytan = wnzh * wnxrad - wnxh * wnzrad
     wnztan = wnxh * wnyrad - wnyh * wnxrad

     ! Vertical loop over T levels

     do k = lpw(iw),nzz

        weight_t(k,irad)   = weight_t(k,irad)   + wrad1
        weight_t(k,irad+1) = weight_t(k,irad+1) + wrad2

        ! Diagnose axisymmetric component of model's own vortex

         thil_ax(k,irad)   =  thil_ax(k,irad)   + wrad1 *  thil(k,iw)
         thil_ax(k,irad+1) =  thil_ax(k,irad+1) + wrad2 *  thil(k,iw)

        theta_ax(k,irad)   = theta_ax(k,irad)   + wrad1 * theta(k,iw)
        theta_ax(k,irad+1) = theta_ax(k,irad+1) + wrad2 * theta(k,iw)

         tair_ax(k,irad)   =  tair_ax(k,irad)   + wrad1 *  tair(k,iw)
         tair_ax(k,irad+1) =  tair_ax(k,irad+1) + wrad2 *  tair(k,iw)

          rrw_ax(k,irad)   =   rrw_ax(k,irad)   + wrad1 *  rr_w(k,iw)
          rrw_ax(k,irad+1) =   rrw_ax(k,irad+1) + wrad2 *  rr_w(k,iw)

          rrv_ax(k,irad)   =   rrv_ax(k,irad)   + wrad1 *  rr_v(k,iw)
          rrv_ax(k,irad+1) =   rrv_ax(k,irad+1) + wrad2 *  rr_v(k,iw)

            w_ax(k,irad)   =     w_ax(k,irad)   + wrad1 *    wc(k,iw)
            w_ax(k,irad+1) =     w_ax(k,irad+1) + wrad2 *    wc(k,iw)

        ss_liq = max(0.,(rr_v(k,iw) * real(rho(k,iw)) / rhovsl(tair(k,iw)-273.15) - 1.0) * 1.e2)

        ssliq_ax(k,irad)   = ssliq_ax(k,irad)   + wrad1 * ss_liq
        ssliq_ax(k,irad+1) = ssliq_ax(k,irad+1) + wrad2 * ss_liq

        ss_thetadif = min(0.,(rhovsl(tair(k,iw)-273.15) - rr_v(k,iw) * real(rho(k,iw)))) &
                    * alvlocp / real(rho(k,iw))

        thdif_ax(k,irad)   = thdif_ax(k,irad)   + wrad1 * ss_thetadif
        thdif_ax(k,irad+1) = thdif_ax(k,irad+1) + wrad2 * ss_thetadif

        ! Diagnose axisymmetric component of model's own vortex

        vtan_ax0 = vxe(k,iw) * wnxtan + vye(k,iw) * wnytan + vze(k,iw) * wnztan
        vrad_ax0 = vxe(k,iw) * wnxrad + vye(k,iw) * wnyrad + vze(k,iw) * wnzrad

        vtan_ax(k,irad)   = vtan_ax(k,irad)   + wrad1 * vtan_ax0
        vtan_ax(k,irad+1) = vtan_ax(k,irad+1) + wrad2 * vtan_ax0

        vrad_ax(k,irad)   = vrad_ax(k,irad)   + wrad1 * vrad_ax0
        vrad_ax(k,irad+1) = vrad_ax(k,irad+1) + wrad2 * vrad_ax0

     enddo

  enddo

  ! Convert sums to averages

  do irad = nr,1,-1

     do k = 2,nzz

        if (weight_t(k,irad) < 1.e-6) then

           if (irad == nr) stop 'stop irad T'

            thil_ax(k,irad) =  thil_ax(k,irad+1)
           theta_ax(k,irad) = theta_ax(k,irad+1)
            tair_ax(k,irad) =  tair_ax(k,irad+1)
             rrw_ax(k,irad) =   rrw_ax(k,irad+1)
             rrv_ax(k,irad) =   rrv_ax(k,irad+1)
               w_ax(k,irad) =     w_ax(k,irad+1)
           ssliq_ax(k,irad) = ssliq_ax(k,irad+1)
           thdif_ax(k,irad) = thdif_ax(k,irad+1)

           if (irad > 1) then

              vtan_ax(k,irad) = vtan_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

              vrad_ax(k,irad) = vrad_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

           endif

        else

            thil_ax(k,irad) =  thil_ax(k,irad) / weight_t(k,irad)
           theta_ax(k,irad) = theta_ax(k,irad) / weight_t(k,irad)
            tair_ax(k,irad) =  tair_ax(k,irad) / weight_t(k,irad)
             rrw_ax(k,irad) =   rrw_ax(k,irad) / weight_t(k,irad)
             rrv_ax(k,irad) =   rrv_ax(k,irad) / weight_t(k,irad)
               w_ax(k,irad) =     w_ax(k,irad) / weight_t(k,irad)
           ssliq_ax(k,irad) = ssliq_ax(k,irad) / weight_t(k,irad)
           thdif_ax(k,irad) = thdif_ax(k,irad) / weight_t(k,irad)

           if (irad > 1) then

              vtan_ax(k,irad) = vtan_ax(k,irad) / weight_t(k,irad)
              vrad_ax(k,irad) = vrad_ax(k,irad) / weight_t(k,irad)

           endif

        endif

        vtan_ax(k,1) = 0.
        vrad_ax(k,1) = 0.

     enddo

     write(6,'(a,i5,3f10.1)') 'irad, radius_ax, vtan_ax, vrad_ax ', &
        irad, radius_ax(irad), vtan_ax(5,irad), vrad_ax(5,irad)

  enddo

  end subroutine vortex_azim_avg

!==============================================================================

  subroutine vortex_add_pert()

  ! This subroutine adds perturbation fields to a model initial state that act to
  ! restore the approximate intensity of a tropical cyclone that is poorly resolved in
  ! the standard analysis used for the initial fields (e.g., the GFS analysis).

  ! The primary perturbation field is the tangential wind component of an
  ! axysymmetric vortex whose peak value and vertical and radial profiles are
  ! determined by 6 user-specified parameters.  A secondary axisymmetric
  ! perturbation field of potential temperature is computed in this subroutine
  ! that, when added to the background potential temperature field will maintain 
  ! approximate gradient wind balance with the enhanced tangential winds.  A
  ! hydrostatic balance procedure is carried out following the addition of these
  ! perturbations in order to adjust pressure and density to the new potential
  ! temperature field.

  ! This subroutine is called at the beginning of the first dynamic initialization
  ! cycle.  Optionally, it may also be called at the beginning of subsequent cycles,
  ! although in such cases the added perturbation should be relatively weak.  The
  ! preferred application of the dynamic initialization procedure is to choose
  ! (perhaps by trial and error) a perturbation on the first cycle that will 
  ! eliminate the need for any perturbations to be added on subsequent cycles.

  use mem_ijtabs,   only: jtw_init, jtab_w, jtv_init, jtab_v, itab_v
  use mem_basic,    only: thil, theta, tair, rho, press, rr_v, rr_w, &
                          wc, wmc, vc, vmc
  use misc_coms,    only: iparallel
  use mem_grid,     only: mza, lpw, lpv, xew, yew, zew, xev, yev, zev, &
                          zt, gdz_belo8, gdz_abov8, dzim, vnx, vny, vnz
  use consts_coms,  only: r8, erad, pio180, pi2, omega, omega2, grav, cvocp, &
                          cp, p00kord, eps_vapi, p00i, rocp, alvlocp, eps_vapi
  use therm_lib,    only: rhovsl
  use mem_micro,    only: rr_c, con_c, cldnum
  use micro_coms,   only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
  use vel_t3d,      only: diagvel_t3d
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_v, mpi_recv_v

  use obnd,         only: lbcopy_v, lbcopy_w

  use plotcolors,   only: make_colortable

  ! Define initial perturbation using analytical functions

  implicit none

  integer :: iw,j,k,iv,iw1,iw2,iter
  integer :: ir
  real :: exner, temp, ccn
  real :: wnxh,wnyh,wnzh
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan

  real :: rad

  real :: wt_vert

  real :: radax_eyw, radax_env, countk
  real :: circ_eyw, circ_env_avg
  real :: wrad1, wrad2

  real :: circ_targ  (nr)
  real :: vtan_targ  (nr)
  real :: vtan_ax_avg(nr)

  real :: vtan_pert (nz,nr)
  real :: exner_pert(nz,nr)
  real :: theta_pert(nz,nr)

  real :: delrad, radmid, thetamid, vtanmid, vpertmid, vtotmid

  real :: zmin_circav, zmax_circav
  real :: fcor, omeglat
  real :: exner_up, exner_dn

  real(r8) :: pkhyd, rho_tot(mza)

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh0 / erad
  wnyh = yeh0 / erad
  wnzh = zeh0 / erad

  radax_eyw = radius_ax(ir_eyw)
  radax_env = radius_ax(ir_env)

  fcor    = omega2 * sin(hlat * pio180)
  omeglat = omega  * sin(hlat * pio180)

  ! Average vtan_ax over selected vertical levels

  zmin_circav = 200.
  zmax_circav = 1000.
  countk = 0.
  vtan_ax_avg(:) = 0.

  do k = 2,mza
     if (zt(k) >= zmin_circav .and. zt(k) <= zmax_circav) then
        countk = countk + 1.
        do ir = 1,nr
           vtan_ax_avg(ir) = vtan_ax_avg(ir) + vtan_ax(k,ir)
        enddo
     endif
  enddo
  vtan_ax_avg(1:nr) = vtan_ax_avg(1:nr) / countk

  ! Compute target radial profile of circulation between eyewall and 
  ! "environment radius", and from that, a target radial profile of vtan

  circ_eyw = radax_eyw * vtan_eyw
  circ_env_avg = radax_env * (omeglat * radax_env + vtan_ax_avg(ir_env))

  do ir = ir_eyw, ir_env
     circ_targ(ir) = circ_eyw + (circ_env_avg - circ_eyw) &
        * ((radius_ax(ir) - radax_eyw) / (radax_env - radax_eyw)) ** rexpon_delc

     vtan_targ(ir) = circ_targ(ir) / radius_ax(ir) - omeglat * radius_ax(ir)

     write(6,'(a,i5,3f12.1)') 'ir, radius_ax, circ_targ, vtan_targ ', &
        ir, 1.e-3 * radius_ax(ir), 1.e-3 * circ_targ(ir), vtan_targ(ir)

  enddo
  print*, ' '

  ! Compute target radial profile of vtan inside eyewall radius

  do ir = 1,ir_eyw-1
     vtan_targ(ir) = vtan_eyw * radius_ax(ir) / radax_eyw
  enddo

  ! Evaluate perturbation tangential velocity at radius_ax points)

  vtan_pert(:,:) = 0.
  theta_pert(:,:) = 0.
  exner_pert(:,:) = 0.

  do ir = ir_env-1,1,-1
     do k = 2,nzz

        ! Vertical dependence of tangential velocity perturbation

        if (zt(k) < zmax_delv) then
           wt_vert = 1. - (zt(k) / zmax_delv) ** zexpon_delv
        else
           wt_vert = 0.
        endif

        ! Compute vtan_pert with vertical weighting

        vtan_pert(k,ir) = max(0.,(vtan_targ(ir) - vtan_ax_avg(ir)) * wt_vert)
     enddo

     ! Compute reduction in gradient-wind Exner function due to addition of vtan_pert

     delrad   =        radius_ax(ir+1) - radius_ax(ir)
     radmid   = 0.5 * (radius_ax(ir+1) + radius_ax(ir))

     do k = 2,nzz

        thetamid = 0.5 * (theta_ax(k,ir+1) + theta_ax(k,ir))
        vtanmid  = 0.5 * (vtan_ax (k,ir+1) + vtan_ax (k,ir))
        vpertmid = 0.5 * (vtan_pert(k,ir+1) + vtan_pert(k,ir))
        vtotmid  = vtanmid + vpertmid

        exner_pert(k,ir) = exner_pert(k,ir+1) - delrad / thetamid &
           * ((vtotmid**2 - vtanmid**2) / radmid + fcor * (vtotmid - vtanmid))

     enddo

     ! Compute increase in theta to hydrostatically balance perturbation
     ! Exner function; this is gradient wind relationship

     do k = nzz-1,3,-1
        exner_up = cp * tair_ax(k+1,ir) / theta_ax(k+1,ir) + exner_pert(k+1,ir)
        exner_dn = cp * tair_ax(k-1,ir) / theta_ax(k-1,ir) + exner_pert(k-1,ir)

        ! Impose zero minimum so that theta_pert does not become negative.
        ! (Near-surface tangential winds are slowed by surface drag and therefore
        ! do not obey gradient wind relation used in getting exner_pert.)

        theta_pert(k,ir) = max(0.,-theta_ax(k,ir) &
           * (exner_pert(k+1,ir) - exner_pert(k-1,ir)) / (exner_up - exner_dn))
     enddo
     theta_pert(2,ir) = theta_pert(3,ir)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f7.2)') 'vtan_ax ',k,(vtan_ax(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f7.2)') 'vtan_pert ',k,(vtan_pert(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f7.2)') 'vtan_tot ',k,(vtan_ax (k,ir) &
                                           + vtan_pert(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f8.2)') 'theta_ax ',k,(theta_ax(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f8.2)') 'tair_ax ',k,(tair_ax(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz,2,-1
     write(6,'(a,i5,30f7.2)') 'exner_pert ',k,(exner_pert(k,ir),ir = 1,ir_env-1)
  enddo

  print*, ' '
  do k = nzz-1,2,-1
     write(6,'(a,i5,30f7.2)') 'theta_pert ',k,(theta_pert(k,ir),ir = 1,ir_env-1)
  enddo
  print*, ' '

  call make_colortable(506,'lin',  0.0, 30.0,2.0,.1)
  call o_reopnwk()
  call plotback()
  call vortex_rzplot('3','vtan_pert '   ,' (m/s)' , vtan_pert  ,1. ,0 ,506)
  call vortex_rzplot('4','theta_pert '  ,' (K)'   , theta_pert ,1. ,0 ,506)
  call o_frame()
  call o_clswk()

  ! Add perturbation to thermodynamic fields

  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)  ! jend(1) = hardwired for mrl 1

     ! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh0)**2 + (yew(iw)-yeh0)**2 + (zew(iw)-zeh0)**2)

     ! Skip hurricane assimilation for all points outside limiting perturbation radius

     if (rad >= radius_ax(nr) - 1) cycle

     ! Determine interpolation point in radial dimension

     ir = 1
     do while (rad > radius_ax(ir+1))
        ir = ir + 1
     enddo

     if (ir < 1 .or. ir+1 > nr) then
        print*, 'ir out of bounds2 ', ir, nr
     endif

     wrad2 = (rad - radius_ax(ir)) / (radius_ax(ir+1) - radius_ax(ir))
     wrad1 = 1. - wrad2

     ! Add theta perturbation to model thil and theta arrays

     do k = lpw(iw),nzz

        ! Skip hurricane assimilation for all points above limiting height

        if (zt(k) >= zmax_delv) exit

         thil(k,iw) =  thil(k,iw) + wrad1 * theta_pert(k,ir) + wrad2 * theta_pert(k,ir+1)
        theta(k,iw) = theta(k,iw) + wrad1 * theta_pert(k,ir) + wrad2 * theta_pert(k,ir+1)

     enddo

     ! Carry out iterative hydrostatic balance procedure

     ! Use nzz level as top boundary condition for hydrostatic integration:
     ! Profile remains unchanged at and above this level 

     if (miclevel == 0) then
        rho_tot(nzz) = rho(nzz,iw)
     else
        rho_tot(nzz) = rho(nzz,iw) * (1. + rr_w(k,iw))
     endif

     do iter = 1,100

        do k = nzz-1,1,-1

           !  Compute density

           if (miclevel == 0) then
              rho(k,iw) = press(k,iw) ** cvocp * p00kord / theta(k,iw)
              rho_tot(k) = rho(k,iw)
           elseif (miclevel == 1) then
              rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                          ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )
              rho_tot(k) = rho(k,iw) * (1. + rr_v(k,iw))
           else
!WW              exner = (real(press(k,iw)) * p00i) ** rocp   ! Defined WITHOUT CP factor
!WW              temp  = exner * theta(k,iw)

!WW              rr_c(k,iw) = max(0., rr_w(k,iw) - rhovsl(temp-273.15) / real(rho(k,iw)))
!WW              rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

              rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                          ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )

              rho_tot(k) = rho(k,iw) * (1. + rr_w(k,iw))
           endif

           ! Hydrostatically integrate downward using weighting to damp oscillations

           pkhyd = press(k+1,iw) &
                 + gdz_belo8(k) * rho_tot(k) + gdz_abov8(k) * rho_tot(k+1)
           press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)

        enddo
     enddo

     ! Vertical loop over T levels

     do k = lpw(iw), nzz
        tair(k,iw) = theta(k,iw) * (real(press(k,iw)) * p00i) ** rocp

        if (miclevel <= 1) then
           thil(k,iw) = theta(k,iw)
        else
           thil(k,iw) = theta(k,iw) / (1. + alvlocp * rr_c(k,iw) / &
                                      ((1.0 + rr_v(k,iw)) * max(tair(k,iw),253.)))
        endif
     enddo

     do k = 1, lpw(iw)-1
        thil(k,iw) = thil(lpw(iw),iw)
     enddo

   ! If there is condensate, initialize con_c if prognosed

     if (miclevel == 3 .and. jnmb(1) == 5) then
        if (ccnparm > 1.e6) then
           ccn = ccnparm
        else
           ccn = cldnum(iw)
        endif

        do k = lpw(iw), mza
           if (rr_c(k,iw) > rxmin(1)) then
              con_c(k,iw) = ccn * real(rho(k,iw)) * zfactor_ccn(k)
           else
              con_c(k,iw) = 0.0
           endif
        enddo
     endif

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_w(1, dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)

     call mpi_recv_w(1, dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)
  endif

  ! LBC copy (THETA and TAIR will be copied later with the scalars)

  call lbcopy_w(1, a1=thil, d1=press, d2=rho)

  do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)  ! jend(1) = hardwired for mrl 1
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
 
     ! Distance of this point from eye center

     rad = sqrt((xev(iv)-xeh0)**2 + (yev(iv)-yeh0)**2 + (zev(iv)-zeh0)**2)

     ! Skip hurricane assimilation for all points outside specified radius

     if (rad >= radax_env) cycle

     ! Determine interpolation point in radial dimension

     ir = 1
     do while (rad > radius_ax(ir+1))
        ir = ir + 1
     enddo

     wrad2 = (rad - radius_ax(ir)) / (radius_ax(ir+1) - radius_ax(ir))
     wrad1 = 1. - wrad2
   
     ! Unit normal vector components from hurricane center to current IV point

     vnxrad = (xev(iv) - xeh0) / rad
     vnyrad = (yev(iv) - yeh0) / rad
     vnzrad = (zev(iv) - zeh0) / rad

     ! Unit vector components in direction of tangential vortex wind

     vnxtan = wnyh * vnzrad - wnzh * vnyrad
     vnytan = wnzh * vnxrad - wnxh * vnzrad
     vnztan = wnxh * vnyrad - wnyh * vnxrad

     ! Vertical loop over T levels

     do k = lpv(iv),nzz

        ! Skip hurricane assimilation for all points above limiting height

        if (zt(k) >= zmax_delv) exit

        ! Add enhanced vortex to winds interpolated from NCEP reanalysis

        vc(k,iv) = vc(k,iv) + (wrad1 * vtan_pert(k,ir) + wrad2 * vtan_pert(k,ir+1)) &
           * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan)

        vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

     enddo

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_v(1, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(1, rvara1=vmc, rvara2=vc)
  endif

  ! LBC copy of VMC, VC

  call lbcopy_v(1, vmc=vmc, vc=vc)

  ! Re-diagnose earth-relative velocities

  call diagvel_t3d(1)

  end subroutine vortex_add_pert

!==================================================================================

  subroutine vortex_reloc3d()

  use mem_basic,   only: vc, wc, thil, theta, rr_w, rr_v, rho, press, &
                         vxe, vye, vze
  use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                         q2, q6, q7, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use mem_ijtabs,  only: itab_m, itab_w, jtab_w, jtw_init
  use misc_coms,   only: io6
  use mem_grid,    only: mza, mma, mwa, lpw, xem, yem, zem, xew, yew, zew
  use consts_coms, only: erad,piu180,pio180
  use max_dims,    only: pathlen

  implicit none

  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan
  real :: vtan,vrad

  ! PS index transfer

  integer :: iwiflag(mwa)

  integer :: im,iw,j,k,kr,jnext,iwnext,ips,jps,npoly,ipt,jpt,lpt
  integer :: ifld, iwi

  real :: rad, rad0, rpolyi

  real :: x(3),y(3),z(3)
  real :: xw(7),yw(7)

  real :: a(mza,nfld),b(mza,nfld),c(mza,nfld)
  real :: field(mza,7,nfld),field_avg(mza,nfld)

  real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
  real :: xwi,ywi

  logical, save :: first_call = .TRUE.

  if (first_call) then
     first_call = .FALSE.

     nout_dim = 0
     nwps(:,:) = 0

     ! Loop over all W points for recording iw indices in iwps array

     do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)  ! jend(1) = hardwired for mrl 1

        ! Distance of this IW point from initial eye center

        rad0 = sqrt((xew(iw)-xeh0)**2 + (yew(iw)-yeh0)**2 + (zew(iw)-zeh0)**2)

        ! Skip hurricane assimilation for all points outside specified radius

        if (rad0 >= rad2_blend + 50.e3) cycle

        ! Transform current W point to PS coordinates tangent at hurricane center

        call e_ps(xew(iw),yew(iw),zew(iw),hlat0,hlon0,xw(1),yw(1))

        ! Get ips,jps indices on PS grid (assumed to be 5 km mesh that is 1000 km wide)

        ips = nint((xw(1) + 500.e3) / 5.e3)
        jps = nint((yw(1) + 500.e3) / 5.e3)

        nout_dim = nout_dim + 1
        nwps(ips,jps) = nwps(ips,jps) + 1

        if (nwps(ips,jps) > 60) then
           print*, 'nwps at ',ips,',',jps,' exceeds 60'
           stop 'stop nwps '
        endif

        ! Store iw index in ips,jps element of iwps array to mark its location

        iwps(ips,jps,nwps(ips,jps)) = iw

     enddo

     write(6,'(/,a,3i8)') 'vortex_reloc3d: nout_dim ',nout_dim,mza,nfld

     allocate (iwout(nout_dim))               ; iwout (:) = 0
     allocate (reloc_field(mza,nout_dim,nfld)); reloc_field(:,:,:) = 0.

  endif

  nout = 0

  ! Find "earth" coordinates of hurricane center current location

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  write(6,'(/,a,2f10.3)') 'vortex_reloc3d: hlat,hlon ',hlat,hlon

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

  iwiflag(:) = 0

  ! Loop over M points for interpolating W points

  do im = 2,mma

     ! Distance of this IM point from current eye center

     rad = sqrt((xem(im)-xeh)**2 + (yem(im)-yeh)**2 + (zem(im)-zeh)**2)
   
     ! Skip hurricane assimilation for all points outside specified radius

     if (rad >= rad2_blend + 30.e3) cycle

     ! Transform M point to PS coordinates tangent at current hurricane center

     call e_ps(xem(im),yem(im),zem(im),hlat,hlon,x(1),y(1))

     npoly = itab_m(im)%npoly
     rpolyi = 1. / real(npoly)

     ! Initialize field average

     field_avg(1:mza,1:nfld) = 0.

     ! Loop over all W points that surround current M point

     do j = 1,npoly

        ! Current W point index

        iw = itab_m(im)%iw(j)

        ! Diagnose tangential and radial velocity components at T points

        ! Unit normal vector components from hurricane center to current IW point

        wnxrad = (xew(iw) - xeh) / rad
        wnyrad = (yew(iw) - yeh) / rad
        wnzrad = (zew(iw) - zeh) / rad

        ! Unit vector components in direction of tangential vortex wind

        wnxtan = wnyh * wnzrad - wnzh * wnyrad
        wnytan = wnzh * wnxrad - wnxh * wnzrad
        wnztan = wnxh * wnyrad - wnyh * wnxrad

        ! Transform W point to PS coordinates tangent at current hurricane center

        call e_ps(xew(iw),yew(iw),zew(iw),hlat,hlon,xw(j),yw(j))

        ! Vertical loop over T levels

        do k = 2,mza
           kr = max(k,lpw(iw))  ! In case TC is too close to land and has underground points

           vtan = vxe(kr,iw) * wnxtan + vye(kr,iw) * wnytan + vze(kr,iw) * wnztan
           vrad = vxe(kr,iw) * wnxrad + vye(kr,iw) * wnyrad + vze(kr,iw) * wnzrad

           field(k,j, 1) =  vtan
           field(k,j, 2) =  vrad
           field(k,j, 3) =    wc(kr,iw)
           field(k,j, 4) =  thil(kr,iw)
           field(k,j, 5) = theta(kr,iw)
           field(k,j, 6) =  rr_w(kr,iw)
           field(k,j, 7) =  rr_v(kr,iw)
           if (allocated(rr_c))  field(k,j, 8) =  rr_c(kr,iw)
           if (allocated(rr_d))  field(k,j, 9) =  rr_d(kr,iw)
           if (allocated(rr_r))  field(k,j,10) =  rr_r(kr,iw)
           if (allocated(rr_p))  field(k,j,11) =  rr_p(kr,iw)
           if (allocated(rr_s))  field(k,j,12) =  rr_s(kr,iw)
           if (allocated(rr_a))  field(k,j,13) =  rr_a(kr,iw)
           if (allocated(rr_g))  field(k,j,14) =  rr_g(kr,iw)
           if (allocated(rr_h))  field(k,j,15) =  rr_h(kr,iw)
           if (allocated(con_c)) field(k,j,16) = con_c(kr,iw)
           if (allocated(con_d)) field(k,j,17) = con_d(kr,iw)
           if (allocated(con_r)) field(k,j,18) = con_r(kr,iw)
           if (allocated(con_p)) field(k,j,19) = con_p(kr,iw)
           if (allocated(con_s)) field(k,j,20) = con_s(kr,iw)
           if (allocated(con_a)) field(k,j,21) = con_a(kr,iw)
           if (allocated(con_g)) field(k,j,22) = con_g(kr,iw)
           if (allocated(con_h)) field(k,j,23) = con_h(kr,iw)
           if (allocated(q2))    field(k,j,24) =    q2(kr,iw)
           if (allocated(q6))    field(k,j,25) =    q6(kr,iw)
           if (allocated(q7))    field(k,j,26) =    q7(kr,iw)
           field(k,j,27) =   rho(kr,iw)
           field(k,j,28) = press(kr,iw)

           field_avg(k,1:nfld) = field_avg(k,1:nfld) + field(k,j,1:nfld) * rpolyi
        enddo

     enddo

     ! Loop over all W points that surround current M point and fill field values

     do j = 1,npoly
        jnext = j + 1
        if (j == npoly) jnext = 1

        iw     = itab_m(im)%iw(j)
        iwnext = itab_m(im)%iw(jnext)

        x(2) = xw(j)
        y(2) = yw(j)

        x(3) = xw(jnext)
        y(3) = yw(jnext)

        ! Loop over vertical levels

        do k = 2,mza

           ! Loop over fields

           do ifld = 1,nfld

              z(1) = field_avg(k,ifld)
              z(2) = field(k,j,ifld)
              z(3) = field(k,jnext,ifld)

              ! Evaluate interpolation coefficients for current trio of points

              call matrix_3x3(1.  , x(1), y(1), &
                              1.  , x(2), y(2), &
                              1.  , x(3), y(3), &
                              z(1), z(2), z(3), &
                              a(k,ifld), b(k,ifld), c(k,ifld))

           enddo

        enddo

        ! Set up some triangle-check coefficients

        v0x = x(2) - x(1)
        v0y = y(2) - y(1)

        v1x = x(3) - x(1)
        v1y = y(3) - y(1)

        dot00 = v0x * v0x + v0y * v0y
        dot01 = v0x * v1x + v0y * v1y
        dot11 = v1x * v1x + v1y * v1y

        denomi = 1. / (dot00 * dot11 - dot01 * dot01)

        ! Get ips,jps indices on PS grid (assumed to be 5 km mesh that is 1000 km wide)
        ! for current M point relative to current hurricane location

        ips = nint((x(1) + 500.e3) / 5.e3)
        jps = nint((y(1) + 500.e3) / 5.e3)

        ! From the ips,jps indices, look up nearby W points in iwps table that
        ! have similar PS coordinates relative to initial hurricane location

        do jpt = jps-1,jps+1
           do ipt = ips-1,ips+1
              do lpt = 1,nwps(ipt,jpt)

                 iwi = iwps(ipt,jpt,lpt)

                 if (iwi < 2 .or. iwi > mwa) &
                    write(6,'(a,6i8)') 'iwi out of bounds ',iwi,mwa,ipt,jpt,lpt,ips

                 ! Distance of this IWI point from initial eye center

                 rad0 = sqrt((xew(iwi)-xeh0)**2 + (yew(iwi)-yeh0)**2 + (zew(iwi)-zeh0)**2)
   
                 ! Skip interpolation for all points outside specified radius

                 if (rad0 >= rad2_blend + 20.e3) cycle

                 ! If this point has already been interpolated, cycle

                 if (iwiflag(iwi) > 0) cycle

                 ! Transform IWI point to PS coordinates tangent at initial hurricane center

                 call e_ps(xew(iwi),yew(iwi),zew(iwi),hlat0,hlon0,xwi,ywi)

                 ! Set up remaining triangle_check coefficients

                 v2x = xwi - x(1)
                 v2y = ywi - y(1)

                 dot02 = v0x * v2x + v0y * v2y
                 dot12 = v1x * v2x + v1y * v2y

                 u = (dot11 * dot02 - dot01 * dot12) * denomi
                 v = (dot00 * dot12 - dot01 * dot02) * denomi

                 ! Check if current qx,qy point is inside or very near current triangle

                 if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

                    ! Point is inside or very near triangle; loop over vertical levels

                    iwiflag(iwi) = 1
                    nout = nout + 1

                    if (nout > nout_dim) then
                       print*, 'In subroutine vortex_reloc3D, nout > nout_dim. '
                       print*, 'nout, nout_dim = ',nout,nout_dim
                       stop 'stop: nout '
                    endif

                    iwout(nout) = iwi

                    do k = 2,mza

                       ! Interpolate to output field point

                       reloc_field(k,nout,1:nfld) = a(k,1:nfld)       &
                                                  + b(k,1:nfld) * xwi &
                                                  + c(k,1:nfld) * ywi    

                    enddo  ! k

                 endif  ! q point inside triangle

              enddo  ! lpt

           enddo   ! ipt

        enddo  ! jpt

     enddo   ! j

   9 continue

  enddo   ! im

  end subroutine vortex_reloc3d

!=======================================================================================

  subroutine vortex_relocated()

  ! This subroutine replaces a model initial state that contains a poorly-resolved
  ! tropical cyclone from a GFS or other analysis with a cyclone that was
  ! simulated in a previous model run and was relocated to the position 
  ! inferred in the initial GFS fields.

  use mem_basic,   only: vc, vmc, wc, wmc, thil, theta, tair, &
                         rr_w, rr_v, rho, press, vxe, vye, vze
  use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                         q2, q6, q7, cldnum, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use micro_coms,  only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
  use therm_lib,   only: rhovsl
  use mem_ijtabs,  only: itab_m, itab_v, itab_w, jtab_v, jtv_init
  use misc_coms,   only: io6, iparallel
  use mem_grid,    only: mza, mma, mwa, lpw, xev, yev, zev, xew, yew, zew, &
                         vnx, vny, vnz, zt, gdz_belo8, gdz_abov8
  use consts_coms, only: erad, piu180, pio180, grav, rvap, rdry, alvlocp, &
                         cvocp, rocp, p00k, p00i, r8, eps_vapi, p00kord
  use max_dims,    only: pathlen
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_v, mpi_recv_v 
  use obnd,         only: lbcopy_v, lbcopy_w
  use vel_t3d,      only: diagvel_t3d

  implicit none

  integer :: iout

  integer :: iw,j,k,iv,iw1,iw2,iter
  real :: exner, temp, ccn
  real :: wnxh,wnyh,wnzh
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan
  real :: wt1, wt2, wt2c, vtant, vradt

  real :: rad0, vreloc

  real(r8) :: pkhyd, rho_tot(mza)

  integer :: mrl

  real :: vtan(mza,mwa),vrad(mza,mwa)

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh0 / erad
  wnyh = yeh0 / erad
  wnzh = zeh0 / erad

  ! Horizontal loop over all active W points in file data 

!  print*, 'rld0 : nout ',nout

  do iout = 1,nout

     iw = iwout(iout)

     ! Distance of this IW point from initial eye center

     rad0 = sqrt((xew(iw)-xeh0)**2 + (yew(iw)-yeh0)**2 + (zew(iw)-zeh0)**2)

     ! Define radial weight coefficients, skipping hurricane assimilation for
     ! all points outside specified radius

     if (rad0 > rad2_blend + 10.e3) cycle

     if (rad0 > rad2_blend) then
        wt1 = 0.
     elseif (rad0 < rad1_blend) then
        wt1 = 1.
     else
        wt1 = (rad2_blend - rad0) / (rad2_blend - rad1_blend)
     endif

     ! Vertical loop over T levels

     do k = 2,mza

        ! Apply height weight coefficient

        if (zt(k) < z1_blend) then
           wt2 = wt1
        elseif (zt(k) > z2_blend) then
           wt2 = 0.
        else
           wt2 = wt1 * (z2_blend - zt(k)) / (z2_blend - z1_blend)
        endif

        wt2c = 1. - wt2

        ! Copy file data to OLAM grid, with weighting

         vtan(k,iw) =                      reloc_field(k,iout, 1)
         vrad(k,iw) =                      reloc_field(k,iout, 2)
           wc(k,iw) =    wc(k,iw) * wt2c + reloc_field(k,iout, 3) * wt2 ! applies at W level
         thil(k,iw) =  thil(k,iw) * wt2c + reloc_field(k,iout, 4) * wt2
        theta(k,iw) = theta(k,iw) * wt2c + reloc_field(k,iout, 5) * wt2
         rr_w(k,iw) =  rr_w(k,iw) * wt2c + reloc_field(k,iout, 6) * wt2
         rr_v(k,iw) =  rr_v(k,iw) * wt2c + reloc_field(k,iout, 7) * wt2
        if (allocated(rr_c)) &
         rr_c(k,iw) =  rr_c(k,iw) * wt2c + reloc_field(k,iout, 8) * wt2
        if (allocated(rr_d)) &
         rr_d(k,iw) =  rr_d(k,iw) * wt2c + reloc_field(k,iout, 9) * wt2
        if (allocated(rr_r)) &
         rr_r(k,iw) =  rr_r(k,iw) * wt2c + reloc_field(k,iout,10) * wt2
        if (allocated(rr_p)) &
         rr_p(k,iw) =  rr_p(k,iw) * wt2c + reloc_field(k,iout,11) * wt2
        if (allocated(rr_s)) &
         rr_s(k,iw) =  rr_s(k,iw) * wt2c + reloc_field(k,iout,12) * wt2
        if (allocated(rr_a)) &
         rr_a(k,iw) =  rr_a(k,iw) * wt2c + reloc_field(k,iout,13) * wt2
        if (allocated(rr_g)) &
         rr_g(k,iw) =  rr_g(k,iw) * wt2c + reloc_field(k,iout,14) * wt2
        if (allocated(rr_h)) &
         rr_h(k,iw) =  rr_h(k,iw) * wt2c + reloc_field(k,iout,15) * wt2
        if (allocated(con_c)) &
        con_c(k,iw) = con_c(k,iw) * wt2c + reloc_field(k,iout,16) * wt2
        if (allocated(con_d)) &
        con_d(k,iw) = con_d(k,iw) * wt2c + reloc_field(k,iout,17) * wt2
        if (allocated(con_r)) &
        con_r(k,iw) = con_r(k,iw) * wt2c + reloc_field(k,iout,18) * wt2
        if (allocated(con_p)) &
        con_p(k,iw) = con_p(k,iw) * wt2c + reloc_field(k,iout,19) * wt2
        if (allocated(con_s)) &
        con_s(k,iw) = con_s(k,iw) * wt2c + reloc_field(k,iout,20) * wt2
        if (allocated(con_a)) &
        con_a(k,iw) = con_a(k,iw) * wt2c + reloc_field(k,iout,21) * wt2
        if (allocated(con_g)) &
        con_g(k,iw) = con_g(k,iw) * wt2c + reloc_field(k,iout,22) * wt2
        if (allocated(con_h)) &
        con_h(k,iw) = con_h(k,iw) * wt2c + reloc_field(k,iout,23) * wt2
        if (allocated(q2)) &
           q2(k,iw) =    q2(k,iw) * wt2c + reloc_field(k,iout,24) * wt2
        if (allocated(q6)) &
           q6(k,iw) =    q6(k,iw) * wt2c + reloc_field(k,iout,25) * wt2
        if (allocated(q7)) &
           q7(k,iw) =    q7(k,iw) * wt2c + reloc_field(k,iout,26) * wt2
          rho(k,iw) =   rho(k,iw) * wt2c + reloc_field(k,iout,27) * wt2
        press(k,iw) = press(k,iw) * wt2c + reloc_field(k,iout,28) * wt2

     enddo

     ! Carry out iterative hydrostatic balance procedure

     if (miclevel == 0) then
        rho_tot(nzz) = rho(nzz,iw)
     else
        rho_tot(nzz) = rho(nzz,iw) * (1. + rr_w(nzz,iw))
     endif

     do iter = 1,100

        ! Use nzz level as top boundary condition for hydrostatic integration:
        ! Profile remains unchanged at and above this level 

        do k = nzz-1,1,-1

           !  Compute density

           if (miclevel == 0) then
              rho(k,iw) = press(k,iw) ** cvocp * p00kord / theta(k,iw)
              rho_tot(k) = rho(k,iw)
           elseif (miclevel == 1) then
              rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                          ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )
              rho_tot(k) = rho(k,iw) * (1. + rr_v(k,iw))
           else
!WW              exner = (real(press(k,iw)) * p00i) ** rocp   ! Defined WITHOUT CP factor
!WW              temp  = exner * theta(k,iw)

!WW              rr_c(k,iw) = max(0., rr_w(k,iw) - rhovsl(temp-273.15) / real(rho(k,iw)))
!WW              rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

              rho(k,iw) = press(k,iw) ** cvocp * p00kord / &
                          ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )

              rho_tot(k) = rho(k,iw) * (1. + rr_w(k,iw))

           endif

           ! Hydrostatically integrate downward using weighting to damp oscillations

           pkhyd = press(k+1,iw) &
                 + gdz_belo8(k) * rho_tot(k) + gdz_abov8(k) * rho_tot(k+1)
           press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)

        enddo
     enddo

     ! Vertical loop over T levels

     do k = lpw(iw),nzz
        wmc(k,iw) = wc(k,iw) * .5 * (rho(k,iw) + rho(k+1,iw))
        tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp

        if (miclevel <= 1) then
           thil(k,iw) = theta(k,iw)
        else
           thil(k,iw) = theta(k,iw) / (1. + alvlocp * rr_c(k,iw) / &
                                      ((1.0 + rr_v(k,iw)) * max(tair(k,iw),253.)))
        endif
     enddo

   ! If there is condensate, initialize con_c if prognosed

     if (miclevel == 3 .and. jnmb(1) == 5) then
        if (ccnparm > 1.e6) then
           ccn = ccnparm
        else
           ccn = cldnum(iw)
        endif

        do k = lpw(iw), mza
           if (rr_c(k,iw) > rxmin(1)) then
              con_c(k,iw) = ccn * real(rho(k,iw)) * zfactor_ccn(k)
           else
              con_c(k,iw) = 0.0
           endif
        enddo
     endif

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     mrl = 1
     call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)

     call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)
  endif

  ! LBC copy (THETA and TAIR will be copied later with the scalars)

  call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

  ! Initialize VC field

  do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)  ! jend(1) = hardwired for mrl 1
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
 
     ! Distance of this point from eye center

     rad0 = sqrt((xev(iv)-xeh0)**2 + (yev(iv)-yeh0)**2 + (zev(iv)-zeh0)**2)

     ! Skip hurricane assimilation for all points outside specified radius

     if (rad0 > rad2_blend) cycle

     ! Unit normal vector components from hurricane center to current IV point

     vnxrad = (xev(iv) - xeh0) / rad0
     vnyrad = (yev(iv) - yeh0) / rad0
     vnzrad = (zev(iv) - zeh0) / rad0

     ! Unit vector components in direction of tangential vortex wind

     vnxtan = wnyh * vnzrad - wnzh * vnyrad
     vnytan = wnzh * vnxrad - wnxh * vnzrad
     vnztan = wnxh * vnyrad - wnyh * vnxrad

     ! Define radial weight coefficients

     if (rad0 < rad1_blend) then
        wt1 = 1.
     elseif (rad0 > rad2_blend) then
        wt1 = 0.
     else
        wt1 = (rad2_blend - rad0) / (rad2_blend - rad1_blend)
     endif

     ! Vertical loop over T levels

     do k = 2,mza

        ! Apply height weight coefficient

        if (zt(k) < z1_blend) then
           wt2 = wt1
        elseif (zt(k) > z2_blend) then
           wt2 = 0.
        else
           wt2 = wt1 * (z2_blend - zt(k)) / (z2_blend - z1_blend)
        endif

        wt2c = 1. - wt2

        ! Average radial and tangential velocity components to V point

        vtant = .5 * (vtan(k,iw1) + vtan(k,iw2))
        vradt = .5 * (vrad(k,iw1) + vrad(k,iw2))

        vreloc = vtant * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan) &
               + vradt * (vnx(iv) * vnxrad + vny(iv) * vnyrad + vnz(iv) * vnzrad)

        vc(k,iv) = vc(k,iv) * wt2c + vreloc * wt2
        vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

     enddo

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_v(1, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(1, rvara1=vmc, rvara2=vc)
  endif

  ! LBC copy of VMC, VC

  call lbcopy_v(1, vmc=vmc, vc=vc)

  ! Re-diagnose earth-relative velocities

  call diagvel_t3d(1)

  end subroutine vortex_relocated

!==================================================================================

  subroutine vortex_rzplot1()

  ! This subroutine calls vortex_rzplot to contour plot azimuthal averages of the TC vortex

  use mem_para,   only: myrank
  use plotcolors, only: make_colortable

  implicit none

  real :: cond_ax(nz,nr)

  if (myrank /= 0) return

  print*, 'beginning vortex_diagnose'

  cond_ax(:,:) = rrw_ax(:,:) - rrv_ax(:,:)

  call make_colortable(500,'lin',280.0,500.0,5.0,1.)
  call make_colortable(501,'lin',  0.0, 50.0,5.0,1.)
  call make_colortable(502,'lin',-20.0, 20.0,2.0,1.)
  call make_colortable(503,'lin', -2.0,  5.0,0.2,1.)
  call make_colortable(504,'lin',250.0,310.0,2.0,1.)
  call make_colortable(505,'lin',  0.0, 30.0,1.0,.1)

  call o_reopnwk()

  !                  panel label    units      field   factor lbc ctab

  call plotback()
  call vortex_rzplot('3','theta1 '  ,' (K)'   , theta_ax ,1.   ,1 ,500)
  call vortex_rzplot('4','vtan1 '   ,' (m/s)' ,  vtan_ax ,1.   ,0 ,501)
  call vortex_rzplot('1','vrad1 '   ,' (m/s)' ,  vrad_ax ,1.   ,0 ,502)
  call vortex_rzplot('2','w1 '      ,' (m/s)' ,     w_ax ,1.   ,0 ,503)
  call o_frame()

  call plotback()
  call vortex_rzplot('3','tair1 '   ,' (K)'   ,  tair_ax ,1.   ,1 ,504)
  call vortex_rzplot('4','rrw1 '    ,' (g/kg)',   rrw_ax ,1.e3 ,1 ,505)
  call vortex_rzplot('1','rrv1 '    ,' (g/kg)',   rrv_ax ,1.e3 ,1 ,505)
  call vortex_rzplot('2','cond1 '   ,' (g/kg)',  cond_ax ,1.e3 ,1 ,505)
  call o_frame()

!  call vortex_rzplot('0','ssliq1 '  ,' (%)'   , ssliq_ax ,1.   ,0 ,443)
!  call vortex_rzplot('0','ss_thdif ',' (K)'   , thdif_ax ,1.   ,0 ,464)

  call o_clswk()

  end subroutine vortex_rzplot1

!==================================================================================

  subroutine vortex_rzplot(panel, label, units, fieldin, factor, lbc, colortab)

  ! This subroutine is a wrapper for plotting radius-height arrays of 
  ! dimension (nzz,nr) that contain azimuthal averages of the TC vortex.
  ! It calls oplot_zxy2 to carry out the actual plot.

  use misc_coms
  use mem_grid,   only: zm, zt
  use consts_coms
  use max_dims,   only: pathlen
  use mem_para,   only: myrank

  implicit none

  character(*), intent(in) :: panel, label, units
  real, intent(in) :: fieldin(nz,nr)
  real, intent(in) :: factor
  integer, intent(in) :: lbc, colortab

  integer :: ifill = 1

  real :: radius(nr), height(nzz), field(nzz,nr)

  real, save :: aspect = 1.0
  real, save :: scalelab = .014

  if (myrank /= 0) return

  radius(1:nr) = radius_ax(1:nr) * 1.e-3

  height(1:nzz) = zt(1:nzz) * 1.e-3
  height(1) = zm(1) * 1.e-3

  field(2:nzz,1:nr) = fieldin(2:nzz,1:nr) * factor
  field(1,1:nr) = field(2,1:nr)
  if (lbc == 0) field(1,1:nr) = 0.
  
  call oplot_zxy2(trim(panel),'N','a','c',aspect,scalelab,                 &
                  trim(label),trim(units), nr, nzz, radius, height, & 
                  'radius (km)','height (km)', &
                  field,    & ! 11N
                  colortab,ifill,0.,radius(nr),20.,5,               &
                  0.0,height(nzz),1.0,5  )
 
  end subroutine vortex_rzplot

end module hcane_rz
