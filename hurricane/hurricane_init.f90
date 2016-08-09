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

!------------------------------------------------------------------------------
! The purpose of module hcane_rz is to restore the intensity and eyewall
! diameter of a hurricane that is under-resolved in the initial conditions
! of an OLAM simulation.  This situation commonly arises when the only
! gridded initialization dataset is of low resolution (for a hurricane), such
! as the CFS global gridded reanalysis.  The procedure consists of one or more
! dynamic initialization cycles in which (1) OLAM is first initialized by
! interpolation from the CFS (or other) gridded dataset, (2) the user adds a
! perturbation axisymmetric tangential wind field to the under-resolved
! vortex such that the sum is close to in-situ observations of tangential winds
! and the eyewall diameter, (3) OLAM adds a compensating potential temperature
! perturbation field that approximately maintains hydrostatic and gradient wind
! balance with the strengthened tangential winds, (4) OLAM is integrated forward
! in time a few hours to allow the hurricane to more completely adjust to the
! added perturbation fields and to develop internal dynamic, thermodynamic, and
! moisture structures, and (5) the simulated hurricane is remapped from its
! final simulated location back to the grid cells at its initial location in
! preparation for the next cycle, but the result is only written to a file.
!
! Each cycle is performed as a separate OLAM run, and the user should compare
! plots from each run against in-situ hurricane measurements in order to
! determine how to construct or modify the perturbation to be added in the 
! subsequent cycle.  For each of these runs, or to run simulations without
! dynamic initialization, the parameter INIT_HURR_STEP must be set as follows:
!
! 0 = No hurricane enhancement or tracking; use for most OLAM simulations.
! 1 = For FIRST CYCLE of dynamic initialization procedure.
!     OLAM initial fields are interpolated from CFSR (or other) dataset.
!     An axisymmetric perturbation tangential wind field is added by the user.   
!     OLAM runs short (e.g. 6 hrs) forward integration while tracking hurricane location every timestep.
!     Axisymmetric hurricane fields are diagnosed and plotted at hourly intervals.
!     3D hurricane fields are remapped at hourly intervals from grid cells at present
!     location to grid cells at initial location, but are only written to files.
! 2 = For SUBSEQUENT CYCLES of dynamic initialization procedure.
!     OLAM initial fields are interpolated from CFSR (or other) dataset.
!     Remapped 3D wind fields from previous cycle are read from files and replace 
!     fields interpolated from CFSR where the hurricane is located.
!     An axisymmetric perturbation tangential wind field is added by the user.   
!     OLAM runs short (e.g. 6 hrs) forward integration while tracking hurricane location every timestep.
!     Axisymmetric hurricane fields are diagnosed and plotted at hourly intervals.
!     3D hurricane fields are remapped at hourly intervals from grid cells at present
!     location to grid cells at initial location, but are only written to files.
! 3 = For HURRICANE SIMULATION after completion of dynamic initialization cycles.
!     OLAM initial fields are interpolated from CFSR (or other) dataset.
!     Remapped 3D wind fields from previous cycle are read from files and replace 
!     fields interpolated from CFSR where the hurricane is located.
!     OLAM simulation proceeds normally with no hurricane tracking.
!
! Comments in subroutine hurricane_init (below) provide additional details.
!------------------------------------------------------------------------------

integer, parameter :: init_hurr_step = 0

! Set lat/lon coords of eye center (correct observed location)

! Frances location on 1 Sept 2004 at 00 UTC
real, parameter :: hcentlat = 20.6, hcentlon = -66.3

! Frances location on 2 Sept 2004 at 00 UTC
! real, parameter :: hcentlat = 22.3, hcentlon = -71.4

integer, parameter :: nz = 60, nr = 28 ! Number of vertical, radial points in vortex profiles
integer :: nzz

real :: hlat0, hlon0, hlat, hlon

real :: circ_avg(nr)

! radius_ax is an array of radial distances of the radial-height grid on which
! axisymmetric hurricane fields are diagnosed and modified by perturbations
! defined by the user.  The radial spacing between these values needs to be at
! least twice the grid spacing of the OLAM hexagonal grid where the hurricane
! is located.

real :: radius_ax(nr) = (/ &
   0.e3,   5.e3,  10.e3,  15.e3,  20.e3,  25.e3,  30.e3,  35.e3,  40.e3,  45.e3, &
  50.e3,  55.e3,  60.e3,  70.e3,  80.e3,  90.e3, 100.e3, 120.e3, 140.e3, 160.e3, &
 180.e3, 200.e3, 250.e3, 300.e3, 350.e3, 400.e3, 450.e3, 500.e3/)

! Axisymmmetric vortex profile arrays

  real ::  thil_ax1(nz,nr)
  real :: theta_ax1(nz,nr)
  real ::   shw_ax1(nz,nr)
  real ::   shv_ax1(nz,nr)
  real ::  vtan_ax1(nz,nr)
  real ::  vrad_ax1(nz,nr)
  real ::     w_ax1(nz,nr)

! Relocation fields for I/O

  integer, parameter :: nfld = 28   ! number of fields transferred
  integer, parameter :: nwi = 300000 ! Horizontal # of W points filled 
                                    ! by interpolation
Contains

!===============================================================================

  subroutine hurricane_init()

  use mem_grid,  only: mza, zt
  use misc_coms, only: iparallel, runtype

  implicit none

  integer, save :: newcall = 0

  integer :: k

! If this is first call to hurricane_init for this model run, count number
! of model levels that are below the chosen maximum height for the dynamic
! initialization procedure (e.g. 20 km).

  if (newcall /= 1) then
     newcall = 1

     do k = 2,mza
        nzz = k
        if (zt(k) > 20000.) exit
     enddo
  endif

  if (runtype /= 'INITIAL') return

  ! Each cycle of a hurricane dynamic initialization procedure is performed 
  ! as a separate run of the model, with RUNTYPE set to 'INITIAL'.  At this
  ! stage of each run, OLAM fields have been interpolated from CFSR or other
  ! gridded data, and they contain that dataset's low resolution representation
  ! of the hurricane, but no additional enhancement has yet been added.

  if (init_hurr_step == 1) then

     ! For the first cycle of the dynamic initialization procedure, 
     ! init_hurr_step is set to 1.

     ! hcentlat and hcentlon are observed coordinates of hurricane that are set
     ! by the user in a parameter statement above; copy these to hlat and hlon.

     hlat = hcentlat
     hlon = hcentlon

     ! Diagnose vortex center location as it is represented on the OLAM grid
     ! after initial interpolation from CFSR or other dataset.  This usually
     ! gives (slighly) different values of hlat and hlon.

     call vortex_center_diagnose()

     ! Copy hlat and hlon to hlat0 and hlon0, which will permanently save
     ! initial hurricane location based on CFSR data interpolated to OLAM grid.

     hlat0 = hlat
     hlon0 = hlon

  elseif (init_hurr_step == 2 .or. init_hurr_step == 3) then

     ! For the second and subsequent dynamic initialization cycles, 
     ! init_hurr_step is set to 2, while for the beginning of the model
     ! simulation after all cycles are completed, init_hurr_step is set to 3.

     ! Read interpolated 3D hurricane fields from one of the files that were 
     ! written during the forward integration of the previous dynamic
     ! initialization cycle, and relocate the fields to the initial position
     ! of the hurricane, specified by hlat0 and hlon0 which are also read from
     ! the file.  The user must first decide which of the files to read from
     ! and must modify subroutine hurricane_init2C so that it reads that file.

     call hurricane_init2C()

     ! Set hlat and hlon to the initial hurricane location.

     hlat = hlat0
     hlon = hlon0

  endif  

  if (init_hurr_step == 1 .or. init_hurr_step == 2) then
  
     ! For either the first or subsequent dynamic initialization cycles,
     ! perform azimuthal average of tangential, radial, and vertical velocity
     ! components, potential temperature, and moisture and plot each in
     ! radial-height cross section.

     call vortex_diagnose()

     ! Add axisymmetric tangential wind field to increase vortex intensity
     ! and/or to change its radius of maximum wind.  Subroutine
     ! hurricane_init1B requires user modification to specify the desired
     ! radial-height profile of the added perturbation, and the subroutine
     ! then adds a potential temperature perturbation that balances the
     ! increased tangential winds through hydrostatic and gradient wind
     ! balances.

     call hurricane_init1B()

     ! Following the addition of the perturbations, re-diagnose and plot
     ! the azimuthally-averaged fields.

     call vortex_diagnose()

     ! Based on what is seen in the plotted azimuthally-averaged fields and
     ! how they compare to in-situ measurements of the hurricane, the user
     ! may want to re-run the current cycle of the dynamic initialization
     ! procedure with a different wind perturbation profile in subroutine 
     ! hurricane_init1B.  If the user sets timmax = 0., the model will stop
     ! at this point, allowing this comparison and decision to be made.

     ! Otherwise, if timmax is set to a positive value (typically 6 hours), 
     ! the model will proceed next with the forward integration phase of
     ! the current dynamic initialization cycle, during which it will track
     ! the hurricane location and, at regular intervals (usually chosen to be
     ! hourly), output hurricane fields.  See comments in subroutine
     ! hurricane_track for more explanation of this process.

  endif

  end subroutine hurricane_init

!==================================================================================

  subroutine vortex_center_diagnose()

  use mem_ijtabs
  use mem_basic
  use misc_coms
  use mem_grid
  use consts_coms

  implicit none

  integer :: iw,k

  real :: reh,xeh,yeh,zeh
  real :: dist,area_tot,weight
  real :: rlon,rlat

  real :: press_min(mza),press_avg(mza),press_thresh(mza)
  real :: xew_avg(mza),yew_avg(mza),zew_avg(mza),weight_sum(mza)

  print*, 'vortex_center_diagnose BEGIN ',hlat,hlon

! Find "earth" coordinates of first-guess position (hlat,hlon)

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth center
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

! Initialize quantities

  area_tot = 0.

  press_min(1:mza) = 2.e5
  press_avg(1:mza) = 0.

  xew_avg(1:mza) = 0.
  yew_avg(1:mza) = 0.
  zew_avg(1:mza) = 0.

  weight_sum(1:mza) = 0.

! Horizontal loop over all W points

  do iw = 2,mwa

! Skip current W point if its location is far from first guess position
! (This is rough check that eliminates moist points)

     if (abs(glatw(iw) - hlat) > 10.) cycle
     if (abs(glonw(iw) - hlon) > 10.) cycle

! Distance of current W point to first guess position

     dist = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

! Skip current W point if its location is far from first guess position
! (This is finer check)

     if (dist > 200.e3) cycle

! Sum grid cell area within search area

     area_tot = area_tot + arw0(iw)

! Vertical loop over T levels

     do k = lpw(iw),mza

! Do not apply algorithm above threshold height

        if (zt(k) > 10000.) exit

! Determine lowest pressure at each model level within search area

        if (real(press(k,iw)) < press_min(k)) then
           press_min(k) = real(press(k,iw))
        endif

! Sum (area * pressure) product at each model level within search area

        press_avg(k) = press_avg(k) + arw0(iw) * press(k,iw)

     enddo  ! k

  enddo  ! iw

! Vertical loop over T levels

  do k = 2,mza

! Do not apply algorithm above threshold height

     if (zt(k) > 10000.) exit

! Compute average pressure

     press_avg(k) = press_avg(k) / area_tot

! Compute threshold pressure at 80% of the range from avg to min

     press_thresh(k) = press_avg(k) + .80 * (press_min(k) - press_avg(k))

  enddo

! Horizontal loop over all W points

  do iw = 2,mwa

! Skip current W point if its location is far from first guess position
! (This is rough check that eliminates moist points)

     if (abs(glatw(iw) - hlat) > 10.) cycle
     if (abs(glonw(iw) - hlon) > 10.) cycle

! Distance of current W point to first guess position

     dist = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

! Skip current W point if its location is far from first guess position
! (This is finer check)

     if (dist > 200.e3) cycle

! Vertical loop over T levels

     do k = lpw(iw),mza

! Do not apply algorithm above threshold height

        if (zt(k) > 10000.) exit

! If pressure < threshold value, sum area-pressure weighted grid cell location

        if (press(k,iw) < press_thresh(k)) then

           weight = arw0(iw) * (press_thresh(k) - press(k,iw)) ! ** exponent  (if desired)

           xew_avg(k) = xew_avg(k) + weight * xew(iw)
           yew_avg(k) = yew_avg(k) + weight * yew(iw)
           zew_avg(k) = zew_avg(k) + weight * zew(iw)

           weight_sum(k) = weight_sum(k) + weight

        endif

     enddo  ! k
   
  enddo  ! iw

! Vertical loop over T levels

  do k = 2,mza

! Do not apply algorithm above threshold height

    if (zt(k) > 10000.) exit

! Compute mean location

     xew_avg(k) = xew_avg(k) / weight_sum(k)
     yew_avg(k) = yew_avg(k) / weight_sum(k)
     zew_avg(k) = zew_avg(k) / weight_sum(k)

! Transform mean location to lat/lon coordinates

     call e_ec(xew_avg(k),yew_avg(k),zew_avg(k),rlon,rlat)

     if (k == 2) then
        hlat = rlat
        hlon = rlon
     endif

  enddo

! Also need to check for lpw(iw) > 2

  print*, 'vortex_center_diagnose END ',hlat,hlon,press_min(2)

  return
  end subroutine vortex_center_diagnose

!==================================================================================

  subroutine vortex_diagnose()

  use mem_ijtabs
  use mem_basic
  use misc_coms
  use mem_grid
  use consts_coms
  use max_dims,   only: pathlen
  use mem_para,   only: myrank

  implicit none

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2
  integer :: irad,ir

  integer, save :: ipert = 0

  real :: rad,rrad,wrad1,wrad2,vtan_pert0,shw_pert0
  real :: zero = 0.

  real :: circ(nr)

  ! Only works on node 0 until cplot is parallelized
  if (myrank /= 0) return

  call o_reopnwk()

  call vortex_diagnose0h(thil_ax1, theta_ax1, shw_ax1, shv_ax1, vtan_ax1, &
                         vrad_ax1,     w_ax1)

!  call plotback; call cplot(nz,nr, thil_ax1 , 58,'thil1'  ); call o_frame
  call plotback; call cplot(nz,nr,theta_ax1 , 58,'theta1' ); call o_frame
!  call plotback; call cplot(nz,nr,  shw_ax1 ,  5,'shw1'   ); call o_frame
!  call plotback; call cplot(nz,nr,  shv_ax1 ,  5,'shv1'   ); call o_frame
  call plotback; call cplot(nz,nr, vtan_ax1 ,159,'vtan1'  ); call o_frame
!  call plotback; call cplot(nz,nr, vrad_ax1 ,109,'vrad1'  ); call o_frame
!  call plotback; call cplot(nz,nr,    w_ax1 ,108,'w1'     ); call o_frame

  call o_clswk()

  return
  end subroutine vortex_diagnose

!==================================================================================

  subroutine vortex_diagnose0h(thil_ax, theta_ax, shw_ax, shv_ax, vtan_ax, &
                               vrad_ax,     w_ax)

  use mem_ijtabs
  use mem_basic
  use mem_micro
  use misc_coms
  use mem_grid
  use consts_coms

  implicit none

  real, intent(out) ::  thil_ax(nz,nr) ! ice-liquid potential temperature (K)
  real, intent(out) :: theta_ax(nz,nr) ! potential temperature (K)
  real, intent(out) ::   shw_ax(nz,nr) ! total water specific density (kg/kg)
  real, intent(out) ::   shv_ax(nz,nr) ! vapor specific density (kg/kg)
  real, intent(out) ::  vtan_ax(nz,nr) ! tangential wind (m/s)
  real, intent(out) ::  vrad_ax(nz,nr) ! radial wind (m/s)
  real, intent(out) ::     w_ax(nz,nr) ! vertical wind (m/s)

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan

  real :: rad,rrad,wrad1,wrad2
  real :: vtan_ax0,vrad_ax0

  integer :: irad,ir

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
    shw_ax(:,:) = 0.
    shv_ax(:,:) = 0.
   vtan_ax(:,:) = 0.
   vrad_ax(:,:) = 0.
      w_ax(:,:) = 0.

  weight_t(:,:) = 0.

! Loop over all W points

!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------

! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)
   
! Skip hurricane assimilation for all points outside specified radius

     if (rad >= radius_ax(nr) - 1.) cycle

! Determine interpolation point in radial dimension

     irad = 1
     do while (rad > radius_ax(irad+1))
        irad = irad + 1
     enddo

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

          shw_ax(k,irad)   =   shw_ax(k,irad)   + wrad1 *  sh_w(k,iw)
          shw_ax(k,irad+1) =   shw_ax(k,irad+1) + wrad2 *  sh_w(k,iw)

          shv_ax(k,irad)   =   shv_ax(k,irad)   + wrad1 *  sh_v(k,iw)
          shv_ax(k,irad+1) =   shv_ax(k,irad+1) + wrad2 *  sh_v(k,iw)

            w_ax(k,irad)   =     w_ax(k,irad)   + wrad1 *    wc(k,iw)
            w_ax(k,irad+1) =     w_ax(k,irad+1) + wrad2 *    wc(k,iw)

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
             shw_ax(k,irad) =   shw_ax(k,irad+1)
             shv_ax(k,irad) =   shv_ax(k,irad+1)
               w_ax(k,irad) =     w_ax(k,irad+1)

           if (irad > 1) then

              vtan_ax(k,irad) = vtan_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

              vrad_ax(k,irad) = vrad_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

           endif

        else

            thil_ax(k,irad) =  thil_ax(k,irad) / weight_t(k,irad)
           theta_ax(k,irad) = theta_ax(k,irad) / weight_t(k,irad)
             shw_ax(k,irad) =   shw_ax(k,irad) / weight_t(k,irad)
             shv_ax(k,irad) =   shv_ax(k,irad) / weight_t(k,irad)
               w_ax(k,irad) =     w_ax(k,irad) / weight_t(k,irad)

           if (irad > 1) then

              vtan_ax(k,irad) = vtan_ax(k,irad) / weight_t(k,irad)
              vrad_ax(k,irad) = vrad_ax(k,irad) / weight_t(k,irad)

           endif

        endif

     enddo

  enddo

  vtan_ax(:,1) = 0.
  vrad_ax(:,1) = 0.

  return
  end subroutine vortex_diagnose0h

!==============================================================================

  subroutine hurricane_init1B()

! This subroutine adds perturbation fields to a model initial state that act to
! restore the approximate intensity of a tropical cyclone that is poorly resolved in
! the standard analysis used for the initial fields (e.g., the GFS analysis).

! The added perturbation fields consist of potential temperature, water vapor
! specific density, and cyclonic winds.  The added perturbation fields are 
! axisymmetric, defined as a function of radius and height only.  The perturbation 
! characteristics are controlled by a set of parameters, whose values are specified
! below, and which should be modified by the user according to the particular
! initialization case.

! This subroutine is called only at the beginning of the first of several dynamic
! initialization cycles, an is only intended to get the vortex intensity and size
! somewhere in the ballpark of observations.  Subsequent cycles do not use this
! subroutine but instead use results from the previous cycle.

  use mem_ijtabs
  use mem_basic
  use misc_coms
  use mem_grid
  use consts_coms
  use vel_t3d,      only: diagvel_t3d
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_v, mpi_recv_v 

  use obnd,         only: lbcopy_v, lbcopy_w

! Define initial perturbation using analytical functions

  implicit none

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2,kbc,iter
  integer :: ir
  real :: pkhyd,exner,deltheta
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: unxrad,unyrad,unzrad,unxtan,unytan,unztan
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan

  real :: rad,rrad,vtan_pert0,shw_pert0,zdif

! TANGENTIAL VELOCITY PERTURBATION - assumed to be maximum at (r = r_vmax, z = 0)

  real,    parameter :: delv_eyw = 45.
  real,    parameter :: zmax_delv = 15000.
  real,    parameter :: zexpon_delv = 1.5
  integer, parameter :: ir_eyw = 7     ! radial index in axisymmetric profiles
  integer, parameter :: ir_env = 17    ! radial index in axisymmetric profiles

! real,    parameter :: rexpon_delc = 1.0
  real,    parameter :: rexpon_delc = 0.5

  real :: vpert_rad, vpert_vert

  real :: rad_eyw, rad_env
  real :: circ_eyw, circ_env, dcirc, circ0
  real :: wrad1, wrad2

  real :: del_circ(nr), circ(nr)

  real ::   theta_totw(nz,nr)
  real ::     vtan_tot(nz,nr)

  integer :: irad, mrl
  real :: theta_tot0, theta_tot1, theta_pert0, dbdz, dri
  real :: b(nz)

  print*, 'IN HURRICANE_INIT1 ',hlat,hlon

! Find "earth" coordinates of hurricane center

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

  circ_avg(:) = 0.

  do k = 10,14
     do ir = 1,nr
        circ0 = pi2 * radius_ax(ir) &
              * (omega * sin(hlat * pio180) * radius_ax(ir) &
              + vtan_ax1(k,ir))

! Cyclone circulation in Reanalysis, vertically averaged over 5 levels

        circ_avg(ir) = circ_avg(ir) + .2 * circ0
     enddo
  enddo

  do ir = 1,nr
     write(6,'(a,i5,2f10.1)') 'circ1_avg ',ir, &
        1.e-3*radius_ax(ir),1.e-6*circ_avg(ir)
  enddo

! Compute radial profile of circulation perturbation between eyewall and env.

  rad_eyw = radius_ax(ir_eyw)
  rad_env = radius_ax(ir_env)

  circ_eyw = circ_avg(ir_eyw) + pi2 * rad_eyw * delv_eyw
  circ_env = circ_avg(ir_env)

  del_circ(:) = 0.

  do ir = ir_eyw,ir_env

     circ(ir) = circ_eyw + (circ_env - circ_eyw) &
              * ((radius_ax(ir) - rad_eyw) / (rad_env - rad_eyw)) ** rexpon_delc 
      
     del_circ(ir) = max(0., circ(ir) - circ_avg(ir))
      
     print*, 'hi2 ',ir,1.e-6*circ(ir),1.e-6*circ_avg(ir),1.e-6*del_circ(ir)

  enddo   

! Average theta_ax1 vertically to W levels; result is theta_totw, which
! will be used to compute thermal wind for gradient wind balance 

  do ir = 1,nr
     do k = 2,nz-2
        theta_totw(k,ir) = .5 * (theta_ax1(k,ir) + theta_ax1(k+1,ir))
     enddo
     theta_totw(1   ,ir) = theta_ax1(2   ,ir)
     theta_totw(nz-1,ir) = theta_ax1(nz-1,ir)
  enddo

! Evaluate perturbation tangential velocity (at mid-points between radius_ax points)

  do ir = ir_env-1,1,-1

     rad = .5 * (radius_ax(ir+1) + radius_ax(ir))
     dri = 1. / (radius_ax(ir+1) - radius_ax(ir))

     if (rad <= rad_eyw) then
        vpert_rad = delv_eyw * rad / rad_eyw
     else
        vpert_rad = del_circ(ir) / (pi2 * rad)
     endif

! Vertical loop over T levels

     do k = 2,nzz

! Vertical dependence of tangential velocity perturbation

        if (zt(k) < zmax_delv) then
           vpert_vert = 1. - (zt(k) / zmax_delv) ** zexpon_delv
        else
           vpert_vert = 0.
        endif

! Combine vertical and radial parts of perturbation

        vtan_pert0 = vpert_rad * vpert_vert

! Add enhanced vortex to winds interpolated from reanalysis

        vtan_tot(k,ir) = .5 * (vtan_ax1(k,ir) + vtan_ax1(k,ir+1)) &
                       + vtan_pert0

        b(k) = vtan_tot(k,ir)**2 / rad &
             + omega2 * sin(hlat * pio180) * vtan_tot(k,ir)
     enddo

! Vertical loop over W levels; Evaluate thermal wind theta_totw

     do k = 2,nzz-1

        dbdz = min(0.,(b(k+1) - b(k)) * dzim(k)) ! min function to omit effect of sfc drag

        theta_totw(k,ir) = (- .5 * dbdz * theta_totw(k,ir+1) &
                           + grav * dri * theta_totw(k,ir+1)) &
                         / (.5 * dbdz + grav * dri)

     enddo
     theta_totw(1,ir) = theta_totw(2,ir) 
  enddo

! Add perturbation to thermodynamic fields

!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------

! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)
   
! Skip hurricane assimilation for all points outside limiting perturbation radius

     if (rad >= radius_ax(nr) - 1) cycle

! Determine interpolation point in radial dimension

     irad = 1
     do while (rad > radius_ax(irad+1))
        irad = irad + 1
     enddo

     wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
     wrad1 = 1. - wrad2
   
! Vertical loop over T levels

     do k = lpw(iw),nzz

! Skip hurricane assimilation for all points above limiting height

        if (zt(k) >= zmax_delv) exit

! Interpolate perturbation table values to current grid cell

        theta_tot0 = .5 * (theta_totw(k,irad  ) + theta_totw(k-1,irad  ))
        theta_tot1 = .5 * (theta_totw(k,irad+1) + theta_totw(k-1,irad+1))

        theta_pert0 = wrad1 * (theta_tot0 - theta_ax1(k,irad  )) &
                    + wrad2 * (theta_tot1 - theta_ax1(k,irad+1))

! Add theta and shw perturbations to model fields

         thil(k,iw) =  thil(k,iw) + theta_pert0
        theta(k,iw) = theta(k,iw) + theta_pert0

     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrostatically balance modified initial profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use nzz level as top boundary condition for hydrostatic integration:
! Profile remains unchanged at and above this level 

! Carry out iterative hydrostatic balance procedure

     do iter = 1,100

        do k = nzz-1,1,-1

!  Compute density - (ASSUME MICPHYS LEVEL = 1 FOR NOW)

           rho(k,iw) = press(k,iw) ** cvocp * p00k  &
              / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))

! Hydrostatically integrate downward using weighting to damp oscillations

           pkhyd = press(k+1,iw)  &
                 + grav * (rho(k+1,iw) * dzt_bot(k+1) + rho(k,iw) * dzt_top(k))
           press(k,iw) = .05 * press(k,iw) + .95 * pkhyd

        enddo
     enddo

! Vertical loop over T levels

     do k = lpw(iw), nzz
        tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
     enddo

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

  call lbcopy_w(1, a1=thil, d1=press, d2=rho)

!----------------------------------------------------------------------
  do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)  ! jend(1) = hardwired for mrl 1
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
 
! Distance of this point from eye center

     rad = sqrt((xev(iv)-xeh)**2 + (yev(iv)-yeh)**2 + (zev(iv)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

     if (rad >= rad_env) cycle

! Determine interpolation point in radial dimension

     ir = 1
     do while (rad > radius_ax(ir+1))
        ir = ir + 1
     enddo

     wrad2 = (rad - radius_ax(ir)) / (radius_ax(ir+1) - radius_ax(ir))
     wrad1 = 1. - wrad2
   
! Interpolate circulation perturbation table values to rad

     dcirc = wrad1 * del_circ(ir  ) &
           + wrad2 * del_circ(ir+1)

! Evaluate perturbation tangential velocity

     if (rad <= rad_eyw) then
        vpert_rad = delv_eyw * rad / rad_eyw
     else
        vpert_rad = dcirc / (pi2 * rad)
     endif

! Unit normal vector components from hurricane center to current IV point

     vnxrad = (xev(iv) - xeh) / rad
     vnyrad = (yev(iv) - yeh) / rad
     vnzrad = (zev(iv) - zeh) / rad

! Unit vector components in direction of tangential vortex wind

     vnxtan = wnyh * vnzrad - wnzh * vnyrad
     vnytan = wnzh * vnxrad - wnxh * vnzrad
     vnztan = wnxh * vnyrad - wnyh * vnxrad

! Vertical loop over T levels

     do k = lpv(iv),nzz

! Skip hurricane assimilation for all points above limiting height

        if (zt(k) >= zmax_delv) exit

! Vertical dependence of tangential velocity perturbation

        vpert_vert = 1. - (zt(k) / zmax_delv) ** zexpon_delv

! Combine vertical and radial parts of perturbation

        vtan_pert0 = vpert_rad * vpert_vert

! Add enhanced vortex to winds interpolated from NCEP reanalysis

        vc(k,iv) = vc(k,iv) + vtan_pert0  &
           * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan)

        vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

     enddo

  enddo

  mrl = 1

! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif

! LBC copy of VMC, VC

  call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

  if (allocated(vmp)) vmp(:,:) = vmc(:,:)
  if (allocated(vp )) vp (:,:) = vc (:,:)

! Re-diagnose earth-relative velocities

  call diagvel_t3d(mrl)

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
  endif

  call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)

  return
  end subroutine hurricane_init1B

!==================================================================================

  subroutine vortex_reloc3d()

  use mem_basic,   only: vc, wc, thil, theta, sh_w, sh_v, rho, press, &
                         vxe, vye, vze
  use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                         q2, q6, q7, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use mem_ijtabs,  only: itab_m, itab_w, jtab_w
  use misc_coms,   only: io6
  use mem_grid,    only: mza, mma, mwa, lpw, xem, yem, zem, xew, yew, zew
  use consts_coms, only: erad,piu180,pio180
  use max_dims,   only: pathlen

  implicit none

  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan
  real :: vtan,vrad
  real :: reloc_field(mza,nwi,nfld)

! PS index transfer

  integer, parameter :: nips = 200, njps = 200 ! size of PS arrays

  integer :: nwps(nips,njps)   ! # of iw indices in each ips,jps cell
  integer :: iwps(nips,njps,30) ! iw indices in each ips,jps cell (max of 10)

  integer :: iwiflag(mwa), iwout(mwa)

  integer :: im,iw,i,j,k,jnext,iwnext,ips,jps,npoly,ipt,jpt,lpt
  integer :: ifld, nout, iwi

  integer, save :: ipert = 0

  real :: rad, rpolyi

  real :: x(3),y(3),z(3)
  real :: xw(7),yw(7)

  real :: a(mza,nfld),b(mza,nfld),c(mza,nfld)
  real :: field(mza,7,nfld),field_avg(mza,nfld)

  real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
  real :: xwi,ywi

  character(pathlen) :: fname
  logical :: exans

  ipert = ipert + 1

     if (ipert == 1) then
        fname = 'frances_out1e'
     elseif (ipert == 2) then
        fname = 'frances_out2e'
     elseif (ipert == 3) then
        fname = 'frances_out3e'
     elseif (ipert == 4) then
        fname = 'frances_out4e'
     elseif (ipert == 5) then
        fname = 'frances_out5e'
     elseif (ipert == 6) then
        fname = 'frances_out6e'
     endif

  nwps(:,:) = 0

! Find "earth" coordinates of hurricane center initial location

  zeh = erad * sin(hlat0 * pio180)
  reh = erad * cos(hlat0 * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon0 * pio180)
  yeh = reh  * sin(hlon0 * pio180)

print*, 'hlat0,hlon0 ',hlat0,hlon0,xeh,yeh,zeh

! Loop over all W points

!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------

! Distance of this IW point from initial eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)
   
! Skip hurricane assimilation for all points outside specified radius

     if (rad >= 450.e3) cycle

! Transform current W point to PS coordinates tangent at hurricane center

     call e_ps(xew(iw),yew(iw),zew(iw),hlat0,hlon0,xw(1),yw(1))

! Get ips,jps indices on PS grid (assumed to be 5 km mesh that is 1000 km wide)

     ips = nint((xw(1) + 500.e3) / 5.e3)
     jps = nint((yw(1) + 500.e3) / 5.e3)

     nwps(ips,jps) = nwps(ips,jps) + 1

     if (nwps(ips,jps) > 30) then
        print*, 'nwps at ',ips,',',jps,' exceeds 30'
        stop 'stop nwps '
     endif

! Store iw index in ips,jps element of iwps array to mark its location

     iwps(ips,jps,nwps(ips,jps)) = iw

! write(6,'(a,5i7,3f10.1)') 'nwps1 ',j,iw,ips,jps,nwps(ips,jps),rad,xw(1)/1000.,yw(1)/1000.

  enddo

  nout = 0

! Find "earth" coordinates of hurricane center current location

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

print*, 'hlat,hlon ',hlat,hlon,xeh,yeh,zeh

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

     if (rad >= 420.e3) cycle

! Transform M point to PS coordinates tangent at current hurricane center

     call e_ps(xem(im),yem(im),zem(im),hlat,hlon,x(1),y(1))

     npoly = itab_m(im)%npoly
     rpolyi = 1. / real(npoly)

! Initialize field average and iwiflag

     field_avg(1:mza,1:nfld) = 0.

! Loop over all W points that surround current M point

     do j = 1,npoly

! Current W point index

        iw = itab_m(im)%iw(j)

!---------------------------------------------------------------
! Diagnose tangential and radial velocity components at T points
!---------------------------------------------------------------

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
           vtan = vxe(k,iw) * wnxtan + vye(k,iw) * wnytan + vze(k,iw) * wnztan
           vrad = vxe(k,iw) * wnxrad + vye(k,iw) * wnyrad + vze(k,iw) * wnzrad

           field(k,j, 1) =  vtan
           field(k,j, 2) =  vrad
           field(k,j, 3) =    wc(k,iw)
           field(k,j, 4) =  thil(k,iw)
           field(k,j, 5) = theta(k,iw)
           field(k,j, 6) =  sh_w(k,iw)
           field(k,j, 7) =  sh_v(k,iw)
           if (allocated(sh_c))  field(k,j, 8) =  sh_c(k,iw)
           if (allocated(sh_d))  field(k,j, 9) =  sh_d(k,iw)
           if (allocated(sh_r))  field(k,j,10) =  sh_r(k,iw)
           if (allocated(sh_p))  field(k,j,11) =  sh_p(k,iw)
           if (allocated(sh_s))  field(k,j,12) =  sh_s(k,iw)
           if (allocated(sh_a))  field(k,j,13) =  sh_a(k,iw)
           if (allocated(sh_g))  field(k,j,14) =  sh_g(k,iw)
           if (allocated(sh_h))  field(k,j,15) =  sh_h(k,iw)
           if (allocated(con_c)) field(k,j,16) = con_c(k,iw)
           if (allocated(con_d)) field(k,j,17) = con_d(k,iw)
           if (allocated(con_r)) field(k,j,18) = con_r(k,iw)
           if (allocated(con_p)) field(k,j,19) = con_p(k,iw)
           if (allocated(con_s)) field(k,j,20) = con_s(k,iw)
           if (allocated(con_a)) field(k,j,21) = con_a(k,iw)
           if (allocated(con_g)) field(k,j,22) = con_g(k,iw)
           if (allocated(con_h)) field(k,j,23) = con_h(k,iw)
           if (allocated(q2))    field(k,j,24) =    q2(k,iw)
           if (allocated(q6))    field(k,j,25) =    q6(k,iw)
           if (allocated(q7))    field(k,j,26) =    q7(k,iw)
           field(k,j,27) =   rho(k,iw)
           field(k,j,28) = press(k,iw)

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

                    iwout(nout) = iwi

                    do k = 2,mza

! Interpolate to output field point

                       reloc_field(k,nout,1:nfld) = a(k,1:nfld)       &
                                                  + b(k,1:nfld) * xwi &
                                                  + c(k,1:nfld) * ywi


!  if (k == 2) then
!     write(6,'(a,2i7,8f10.1)') 'fo1 ',nout,iwi,reloc_field(k,nout,1), &
!                                a(k,1),b(k,1),c(k,1),xwi/1000.,ywi/1000.
!  endif



                    enddo  ! k

                 endif  ! q point inside triangle

              enddo  ! lpt

           enddo   ! ipt

        enddo  ! jpt

     enddo   ! j

9    continue

  enddo   ! im

! Write output fields to file

  inquire(file=trim(fname),exist=exans)

  if (.not. exans) then
     open(32,file=trim(fname),status='new',form='unformatted')
     write(32) hlat0,hlon0,nout,iwout,reloc_field
     close(32)
  endif

  return
  end subroutine vortex_reloc3d

!=======================================================================================

  subroutine hurricane_init2C()

! This subroutine replaces a model initial state with a poorly-resolved
! tropical cyclone from a GFS or other analysis with a cyclone that was
! simulated in a previous model run and was relocated to the position 
! inferred in the initial GFS fields.

!! The perturbation fields that are added in this subroutine are generated 
!! by a dynamic initialization cycle, performed by previous integrations of OLAM.
!! The cycle should be repeated whenever major changes are made to the model
!! configuration, such as a change in grid resolution, or for a new vortex
!! initialization time and location.

  use mem_basic,   only: vc, vp, vmc, vmp, wc, wmc, thil, theta, tair, &
                         sh_w, sh_v, rho, press, vxe, vye, vze
  use mem_micro,   only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                         q2, q6, q7, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use mem_ijtabs,  only: itab_m, itab_v, itab_w, jtab_v, jtab_w
  use misc_coms,   only: io6, iparallel
  use mem_grid,    only: mza, mma, mwa, lpw, xev, yev, zev, xew, yew, zew, &
                         vnx, vny, vnz, zt, dzt_top, dzt_bot
  use consts_coms, only: erad, piu180, pio180, grav, rvap, rdry, &
                         cvocp, rocp, p00k, p00i
  use max_dims,    only: pathlen
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_v, mpi_recv_v 
  use obnd,         only: lbcopy_v, lbcopy_w
  use vel_t3d,      only: diagvel_t3d

  implicit none

  integer :: iwout(mwa)
  integer :: nout, iout

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2,kbc,iter
  real :: pkhyd,exner,deltheta
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: unxrad,unyrad,unzrad,unxtan,unytan,unztan
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan
  real :: wt1, wt2, wt2c, vtan0, vrad0

  real :: rad

  character(10) :: cstr
  integer :: istr, mrl
  real :: rstr

  real :: reloc_field(mza,nwi,nfld)

  real :: vtan(mza,mwa),vrad(mza,mwa)

  character(pathlen) :: fname
  logical :: exans

  fname = 'frances_out6e'

! Read relocated vortex data

  inquire(file=trim(fname),exist=exans)

  print*, 'exans ',exans

  if (.not. exans) then
     print*, 'file ',trim(fname),' does not exit '
     stop 'stop hfile '
  else
     open(32,file=trim(fname),status='old',form='unformatted')
     read(32) hlat0,hlon0,nout,iwout,reloc_field
     close(32)
  endif

  print*, 'Init2c - hlat0,hlon0 = ',hlat0,hlon0

! Find "earth" coordinates of hurricane center

  zeh = erad * sin(hlat0 * pio180)
  reh = erad * cos(hlat0 * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon0 * pio180)
  yeh = reh  * sin(hlon0 * pio180)

! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

! Horizontal loop over all active W points in file data 

  do iout = 1,nout

     iw = iwout(iout)

! Distance of this IW point from initial eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

     if (rad > 400.e3) cycle

! Define radial weight coefficients

     if (rad < 300.e3) then
        wt1 = 1.
     elseif (rad > 400.e3) then
        wt1 = 0.
     else
        wt1 = (400.e3 - rad) / 100.e3
     endif

! Vertical loop over T levels

     do k = 2,mza

! Apply height weight coefficient

        if (zt(k) < 13000.) then
           wt2 = wt1
        elseif (zt(k) > 16000.) then
           wt2 = 0.
        else
           wt2 = wt1 * (16000. - zt(k)) / 3000.
        endif

        wt2c = 1. - wt2

! Copy file data to OLAM grid, with weighting

         vtan(k,iw) =                      reloc_field(k,iout, 1)
         vrad(k,iw) =                      reloc_field(k,iout, 2)
           wc(k,iw) =    wc(k,iw) * wt2c + reloc_field(k,iout, 3) * wt2 ! applies at W level
         thil(k,iw) =  thil(k,iw) * wt2c + reloc_field(k,iout, 4) * wt2
        theta(k,iw) = theta(k,iw) * wt2c + reloc_field(k,iout, 5) * wt2
         sh_w(k,iw) =  sh_w(k,iw) * wt2c + reloc_field(k,iout, 6) * wt2
         sh_v(k,iw) =  sh_v(k,iw) * wt2c + reloc_field(k,iout, 7) * wt2
        if (allocated(sh_c)) &
         sh_c(k,iw) =  sh_c(k,iw) * wt2c + reloc_field(k,iout, 8) * wt2
        if (allocated(sh_d)) &
         sh_d(k,iw) =  sh_d(k,iw) * wt2c + reloc_field(k,iout, 9) * wt2
        if (allocated(sh_r)) &
         sh_r(k,iw) =  sh_r(k,iw) * wt2c + reloc_field(k,iout,10) * wt2
        if (allocated(sh_p)) &
         sh_p(k,iw) =  sh_p(k,iw) * wt2c + reloc_field(k,iout,11) * wt2
        if (allocated(sh_s)) &
         sh_s(k,iw) =  sh_s(k,iw) * wt2c + reloc_field(k,iout,12) * wt2
        if (allocated(sh_a)) &
         sh_a(k,iw) =  sh_a(k,iw) * wt2c + reloc_field(k,iout,13) * wt2
        if (allocated(sh_g)) &
         sh_g(k,iw) =  sh_g(k,iw) * wt2c + reloc_field(k,iout,14) * wt2
        if (allocated(sh_h)) &
         sh_h(k,iw) =  sh_h(k,iw) * wt2c + reloc_field(k,iout,15) * wt2
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hydrostatically balance modified initial profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use nzz level as top boundary condition for hydrostatic integration:
! Profile remains unchanged at and above this level 

! Carry out iterative hydrostatic balance procedure

     do iter = 1,100

        do k = nzz-1,1,-1

!  Compute density - (ASSUME MICPHYS LEVEL = 1 FOR NOW)

           rho(k,iw) = press(k,iw) ** cvocp * p00k  &
              / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))

! Hydrostatically integrate downward using weighting to damp oscillations

           pkhyd = press(k+1,iw)  &
              + grav * (rho(k+1,iw) * dzt_bot(k+1) + rho(k,iw) * dzt_top(k))
           press(k,iw) = .05 * press(k,iw) + .95 * pkhyd

        enddo
     enddo

! Vertical loop over T levels

     do k = lpw(iw),nzz
        wmc(k,iw) = wc(k,iw) * .5 * (rho(k,iw) + rho(k+1,iw))
        tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
     enddo

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

!----------------------------------------------------------------------
  do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)  ! jend(1) = hardwired for mrl 1
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
 
! Distance of this point from eye center

     rad = sqrt((xev(iv)-xeh)**2 + (yev(iv)-yeh)**2 + (zev(iv)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

     if (rad > 400.e3) cycle

! Unit normal vector components from hurricane center to current IV point

     vnxrad = (xev(iv) - xeh) / rad
     vnyrad = (yev(iv) - yeh) / rad
     vnzrad = (zev(iv) - zeh) / rad

! Unit vector components in direction of tangential vortex wind

     vnxtan = wnyh * vnzrad - wnzh * vnyrad
     vnytan = wnzh * vnxrad - wnxh * vnzrad
     vnztan = wnxh * vnyrad - wnyh * vnxrad

! Define radial weight coefficients

     if (rad < 300.e3) then
        wt1 = 1.
     elseif (rad > 400.e3) then
        wt1 = 0.
     else
        wt1 = (400.e3 - rad) / 100.e3
     endif

! Vertical loop over T levels

     do k = 2,mza

! Apply height weight coefficient

        if (zt(k) < 13000.) then
           wt2 = wt1
        elseif (zt(k) > 16000.) then
           wt2 = 0.
        else
           wt2 = wt1 * (16000. - zt(k)) / 3000.
        endif

        wt2c = 1. - wt2

! Average radial and tangential velocity components to V point

        vtan0 = .5 * (vtan(k,iw1) + vtan(k,iw2))
        vrad0 = .5 * (vrad(k,iw1) + vrad(k,iw2))

        vc(k,iv) = vc(k,iv) * wt2c &
                 + vtan0 * wt2 &
                 * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan) &
                 + vrad0 &
                 * (vnx(iv) * vnxrad + vny(iv) * vnyrad + vnz(iv) * vnzrad)

        vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

     enddo

  enddo

  mrl = 1

! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
  endif

! LBC copy of VMC, VC

  call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

  if (allocated(vmp)) vmp(:,:) = vmc(:,:)
  if (allocated(vp))  vp (:,:) = vc (:,:)

! Re-diagnose earth-relative velocities

  call diagvel_t3d(mrl)

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
  endif

  call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)

  return
  end subroutine hurricane_init2C

!==============================================================================

  subroutine cplot(nz,nr,slab,itab,fldname)

  use oplot_coms, only: op
  use plotcolors, only: clrtab
  use mem_grid,   only: zm, zt, lpw
  use misc_coms,  only: time8

  implicit none

  integer, intent(in) :: nz,nr
  real, intent(in) :: slab(nz,nr)
  integer, intent(in) :: itab

  character(len=*), intent(in) :: fldname

  real :: hcpn(4),vcpn(4),fldvals(4),dst(6),ind(8)
  integer :: iasf(18)

  integer :: k,km,i,ival,icolor,ibox,ln,nnn
  real :: bsize,yinc

  data iasf / 18*1 /
  character*8 :: number,numbr
  character*60 :: title

! turn off the clipping indicator.

  call o_gsclip (0)

! set all the gks aspect source flags to "individual".

  call o_gsasf (iasf)

! force solid fill.

  call o_gsfais (1)

  call o_set(.01,.87,.01,.96,0.,radius_ax(nr),0.,zm(nzz),1)

  op%xmin = 0.
  op%xmax = radius_ax(nr)
  op%ymin = 0.
  op%ymax = zm(nzz)

  call o_sfseti ('type of fill',0)

  do k = 2,nzz
     km = max(k-1,2)

     do i = 2,nr

        hcpn(1) = radius_ax(i-1)
        hcpn(2) = radius_ax(i)
        hcpn(3) = hcpn(2)
        hcpn(4) = hcpn(1)

        vcpn(1) = zt(km)
        if (k == 2) vcpn(1) = zm(k-1)
        vcpn(3) = zt(k)
        vcpn(2) = vcpn(1)
        vcpn(4) = vcpn(3)

        fldvals(1) = slab(km,i-1)
        fldvals(2) = slab(km,i  )
        fldvals(3) = slab(k ,i  )
        fldvals(4) = slab(k ,i-1)

! Convert atmospheric moisture values from kg/kg to g/kg

        if (fldname(1:2) == 'sh') then
           fldvals(1:4) = fldvals(1:4) * 1.e3
        endif

        call contpolyg(itab,1,4,hcpn,vcpn,fldvals)

     enddo
  enddo

  call o_sflush

! draw a color bar for the plot.

  call o_set (0.,1.,0.,1.,0.,1.,0.,1.,1)

  call o_gsplci(8)
  call o_gstxci(8)

  hcpn(1) = .88
  hcpn(2) = .91
  hcpn(3) = .91
  hcpn(4) = .88

  bsize = -.63

  yinc = .93 / float(clrtab(itab)%nvals)

  do ibox = 1,clrtab(itab)%nvals

     vcpn(1) = .01 + float(ibox-1) * yinc
     vcpn(2) = vcpn(1)
     vcpn(3) = vcpn(2) + yinc
     vcpn(4) = vcpn(3)
      
     write (number,'(f7.1)') clrtab(itab)%vals(ibox)

     numbr = adjustl(number)

     if(numbr(len_trim(numbr):len_trim(numbr)) == '.') then
        ln = len_trim(numbr) - 1
     else
        ln = len_trim(numbr)
     endif

     if (ibox < clrtab(itab)%nvals) then
        call o_plchlq (hcpn(2)+.01,vcpn(3),numbr(1:ln),bsize,0.,-1.)
     endif

     call o_sfsgfa (hcpn,vcpn,4,clrtab(itab)%ipal(ibox))
!    call fillpolyg (4,hcpn,vcpn,ipal(ibox))

  enddo

! draw a labelbar for the plot

  call o_sflush()

  call o_gsplci(10)
  call o_gstxci(10)
  call o_sflush()

  bsize = .015

  write(title, '(a,i6,a)') trim(fldname),int(real(time8)),' sec'

  call o_plchhq(.05,.97,trim(title),bsize,0.,-1.)

  call o_sflush()

  return
  end subroutine cplot

end module hcane_rz
