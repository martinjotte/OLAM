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
module hcane_rz

implicit none

!----------------------------------------------------------------
! SET INIT_HURR_STEP HERE TO PERFORM THE FOLLOWING:
!
! 0 = no hurricane initialization; no tracking
! 1 = initialize for first dynamic initialization cycle, track hurricane
! 2 = initialize for subsequent dynamic initialization cycles; track hurricane
! 3 = initialize for simulation; do not track hurricane
!----------------------------------------------------------------

integer, parameter :: init_hurr_step = 0

! Set lat/lon coords of eye center (correct observed location)

! Frances location on 1 Sept 2004 at 00 UTC
real, parameter :: hcentlat = 20.6, hcentlon = -66.3

! Frances location on 2 Sept 2004 at 00 UTC
! real, parameter :: hcentlat = 22.3, hcentlon = -71.4

integer, parameter :: nz = 60, nr = 28 ! Number of vertical, radial points in vortex profiles
integer :: nzz

real :: hlat, hlon

real :: circ_avg(nr)

real :: radius_ax(nr) = (/ &
   0.e3,   5.e3,  10.e3,  15.e3,  20.e3,  25.e3,  30.e3,  35.e3,  40.e3,  45.e3, &
  50.e3,  55.e3,  60.e3,  70.e3,  80.e3,  90.e3, 100.e3, 120.e3, 140.e3, 160.e3, &
 180.e3, 200.e3, 250.e3, 300.e3, 350.e3, 400.e3, 450.e3, 500.e3/)

! Axisymmmetric vortex profile arrays

  real ::  thil_ax1(nz,nr),  thil_ax2(nz,nr),  thil_ax12(nz,nr)
  real :: theta_ax1(nz,nr), theta_ax2(nz,nr), theta_ax12(nz,nr)
  real ::   shw_ax1(nz,nr),   shw_ax2(nz,nr),   shw_ax12(nz,nr)
  real ::   shv_ax1(nz,nr),   shv_ax2(nz,nr),   shv_ax12(nz,nr)
  real ::  vtan_ax1(nz,nr),  vtan_ax2(nz,nr),  vtan_ax12(nz,nr)
  real ::  vrad_ax1(nz,nr),  vrad_ax2(nz,nr),  vrad_ax12(nz,nr)
  real ::     w_ax1(nz,nr),     w_ax2(nz,nr),     w_ax12(nz,nr)
  real ::   shc_ax1(nz,nr),   shc_ax2(nz,nr),   shc_ax12(nz,nr)
  real ::   shd_ax1(nz,nr),   shd_ax2(nz,nr),   shd_ax12(nz,nr)
  real ::   shr_ax1(nz,nr),   shr_ax2(nz,nr),   shr_ax12(nz,nr)
  real ::   shp_ax1(nz,nr),   shp_ax2(nz,nr),   shp_ax12(nz,nr)
  real ::   shs_ax1(nz,nr),   shs_ax2(nz,nr),   shs_ax12(nz,nr)
  real ::   sha_ax1(nz,nr),   sha_ax2(nz,nr),   sha_ax12(nz,nr)
  real ::   shg_ax1(nz,nr),   shg_ax2(nz,nr),   shg_ax12(nz,nr)
  real ::   shh_ax1(nz,nr),   shh_ax2(nz,nr),   shh_ax12(nz,nr)

Contains

!===============================================================================

  subroutine hurricane_init()

  use mem_grid,  only: mza, zt
  use mem_basic, only: ump, umc, vmp, vmc
  use misc_coms, only: iparallel, meshtype

  implicit none

  integer, save :: newcall = 0

  integer :: k

! If this is first call to hurricane_init, initialize hurricane location
! and find nzz level 

  if (newcall /= 1) then
     newcall = 1

     hlat = hcentlat
     hlon = hcentlon

     do k = 2,mza-1
        nzz = k
        if (zt(k) > 25000.) exit
     enddo
  endif

  if (init_hurr_step == 1) then

! Diagnose vortex as it exists in initial (CFSR) data prior to adding pert.

     call vortex_center_diagnose()

     call vortex_diagnose(1)

! Add first-guess perturbation for first initialization cycle

     call hurricane_init1B()

     call vortex_diagnose(2)

  elseif (init_hurr_step == 2) then

! Add perturbations for second and subsequent initialization cycles

     call vortex_center_diagnose()

     call hurricane_init2B()

     call vortex_diagnose(2)

  elseif (init_hurr_step == 3) then

! Add perturbations for actual simulation

     call vortex_center_diagnose()

     call hurricane_init2B()

  endif

  end subroutine hurricane_init

!==================================================================================

  subroutine hurricane_track()

  use misc_coms, only: meshtype

  implicit none

! Get updated hurricane center location

  call vortex_center_diagnose()

  call vortex_diagnose(2)

  end subroutine hurricane_track

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

     do k = lpw(iw),mza-1

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

  do k = 2,mza-1

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

     do k = lpw(iw),mza-1

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

  do k = 2,mza-1

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

  print*, 'vortex_center_diagnose END ',hlat,hlon

  return
  end subroutine vortex_center_diagnose

!==================================================================================

  subroutine vortex_diagnose(ifield)

  use mem_ijtabs
  use mem_basic
  use misc_coms
  use mem_grid
  use consts_coms
  use max_dims,   only: pathlen

  implicit none

  integer, intent(in) :: ifield

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2
  integer :: irad,ir

  integer, save :: ipert = 0

  real :: rad,rrad,wrad1,wrad2,vtan_pert0,shw_pert0
  real :: zero = 0.

  real :: circ(nr)

  character(pathlen) :: fname
  logical :: exans

  fname = 'frances_diagnosis_1sept00z'

  11 format(/,15x  ,a,28f6.0,/)
  12 format(i3,f8.0,a,28f6.1)
  13 format(/,4x   ,a,28f6.1)

  call o_reopnwk()

  if (ifield == 1) then

     call vortex_diagnose0(thil_ax1, theta_ax1, shw_ax1, shv_ax1, vtan_ax1, &
                           vrad_ax1,     w_ax1, shc_ax1, shd_ax1,  shr_ax1, &
                            shp_ax1,   shs_ax1, sha_ax1, shg_ax1,  shh_ax1)

     inquire(file=trim(fname),exist=exans)

     if (.not. exans) then

        open(32,file=trim(fname),status='new',form='unformatted')
        write(32) thil_ax1, theta_ax1, shw_ax1, shv_ax1, vtan_ax1, &
                  vrad_ax1,     w_ax1, shc_ax1, shd_ax1,  shr_ax1, &
                   shp_ax1,   shs_ax1, sha_ax1, shg_ax1,  shh_ax1
        close(32)
     endif

  elseif (ifield == 2) then

     call vortex_diagnose0(thil_ax2, theta_ax2, shw_ax2, shv_ax2, vtan_ax2, &
                           vrad_ax2,     w_ax2, shc_ax2, shd_ax2,  shr_ax2, &
                            shp_ax2,   shs_ax2, sha_ax2, shg_ax2,  shh_ax2)

     open(32,file=trim(fname),status='old',form='unformatted')
     read(32) thil_ax1, theta_ax1, shw_ax1, shv_ax1, vtan_ax1, &
              vrad_ax1,     w_ax1, shc_ax1, shd_ax1,  shr_ax1, &
               shp_ax1,   shs_ax1, sha_ax1, shg_ax1,  shh_ax1
     close(32)

      thil_ax12(:,:) =  thil_ax2(:,:) -  thil_ax1(:,:)
     theta_ax12(:,:) = theta_ax2(:,:) - theta_ax1(:,:)
       shw_ax12(:,:) =   shw_ax2(:,:) -   shw_ax1(:,:)
       shv_ax12(:,:) =   shv_ax2(:,:) -   shv_ax1(:,:)
      vtan_ax12(:,:) =  vtan_ax2(:,:) -  vtan_ax1(:,:)
      vrad_ax12(:,:) =  vrad_ax2(:,:) -  vrad_ax1(:,:)
         w_ax12(:,:) =     w_ax2(:,:) -     w_ax1(:,:)
       shc_ax12(:,:) =   shc_ax2(:,:) -   shc_ax1(:,:)
       shd_ax12(:,:) =   shd_ax2(:,:) -   shd_ax1(:,:)
       shr_ax12(:,:) =   shr_ax2(:,:) -   shr_ax1(:,:)
       shp_ax12(:,:) =   shp_ax2(:,:) -   shp_ax1(:,:)
       shs_ax12(:,:) =   shs_ax2(:,:) -   shs_ax1(:,:)
       sha_ax12(:,:) =   sha_ax2(:,:) -   sha_ax1(:,:)
       shg_ax12(:,:) =   shg_ax2(:,:) -   shg_ax1(:,:)
       shh_ax12(:,:) =   shh_ax2(:,:) -   shh_ax1(:,:)

      thil_ax12(nzz,:) = 0.
     theta_ax12(nzz,:) = 0.
       shw_ax12(nzz,:) = 0.
       shv_ax12(nzz,:) = 0.
      vtan_ax12(nzz,:) = 0.
      vrad_ax12(nzz,:) = 0.
         w_ax12(nzz,:) = 0.
       shc_ax12(nzz,:) = 0.
       shd_ax12(nzz,:) = 0.
       shr_ax12(nzz,:) = 0.
       shp_ax12(nzz,:) = 0.
       shs_ax12(nzz,:) = 0.
       sha_ax12(nzz,:) = 0.
       shg_ax12(nzz,:) = 0.
       shh_ax12(nzz,:) = 0.

      thil_ax12(:,nr) = 0.
     theta_ax12(:,nr) = 0.
       shw_ax12(:,nr) = 0.
       shv_ax12(:,nr) = 0.
      vtan_ax12(:,nr) = 0.
      vrad_ax12(:,nr) = 0.
         w_ax12(:,nr) = 0.
       shc_ax12(:,nr) = 0.
       shd_ax12(:,nr) = 0.
       shr_ax12(:,nr) = 0.
       shp_ax12(:,nr) = 0.
       shs_ax12(:,nr) = 0.
       sha_ax12(:,nr) = 0.
       shg_ax12(:,nr) = 0.
       shh_ax12(:,nr) = 0.

!---------------------------------------------------------------------
GO TO 15
!---------------------------------------------------------------------

     if (ipert == 0) then
        open(31,file='frances_cyc1_pert0',status='new',form='formatted')
     elseif (ipert == 1) then
        open(31,file='frances_cyc1_pert1',status='new',form='formatted')
     elseif (ipert == 2) then
        open(31,file='frances_cyc1_pert2',status='new',form='formatted')
     elseif (ipert == 3) then
        open(31,file='frances_cyc1_pert3',status='new',form='formatted')
     elseif (ipert == 4) then
        open(31,file='frances_cyc1_pert4',status='new',form='formatted')
     elseif (ipert == 5) then
        open(31,file='frances_cyc1_pert5',status='new',form='formatted')
     endif

     ipert = ipert + 1

     write(31) thil_ax12, theta_ax12, shw_ax12, shv_ax12, vtan_ax12, &
               vrad_ax12,     w_ax12, shc_ax12, shd_ax12,  shr_ax12, &
                shp_ax12,   shs_ax12, sha_ax12, shg_ax12,  shh_ax12

     close(31)

!------------------------------------------------------------------------
15 CONTINUE
!------------------------------------------------------------------------

!     call plotback; call cplot(nz,nr, thil_ax1 , 58,'thil1'  ); call o_frame
!     call plotback; call cplot(nz,nr, thil_ax2 , 58,'thil2'  ); call o_frame
!     call plotback; call cplot(nz,nr, thil_ax12,130,'thil12' ); call o_frame

!     call plotback; call cplot(nz,nr,theta_ax1 , 58,'theta1' ); call o_frame
     call plotback; call cplot(nz,nr,theta_ax2 , 58,'theta2' ); call o_frame
!     call plotback; call cplot(nz,nr,theta_ax12,130,'theta12'); call o_frame

!     call plotback; call cplot(nz,nr,  shw_ax1 ,  5,'shw1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shw_ax2 ,  5,'shw2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shw_ax12,108,'shw12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shv_ax1 ,  5,'shv1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shv_ax2 ,  5,'shv2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shv_ax12,108,'shv12'  ); call o_frame

!     call plotback; call cplot(nz,nr, vtan_ax1 ,159,'vtan1'  ); call o_frame
     call plotback; call cplot(nz,nr, vtan_ax2 ,159,'vtan2'  ); call o_frame
!     call plotback; call cplot(nz,nr, vtan_ax12,159,'vtan12' ); call o_frame

!     call plotback; call cplot(nz,nr, vrad_ax1 ,109,'vrad1'  ); call o_frame
!     call plotback; call cplot(nz,nr, vrad_ax2 ,109,'vrad2'  ); call o_frame
!     call plotback; call cplot(nz,nr, vrad_ax12,109,'vrad12' ); call o_frame

!     call plotback; call cplot(nz,nr,    w_ax1 ,108,'w1'     ); call o_frame
!     call plotback; call cplot(nz,nr,    w_ax2 ,108,'w2'     ); call o_frame
!     call plotback; call cplot(nz,nr,    w_ax12,108,'w12'    ); call o_frame

!     call plotback; call cplot(nz,nr,  shc_ax1 ,  3,'shc1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shc_ax2 ,  5,'shc2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shc_ax12,108,'shc12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shd_ax1 ,  3,'shd1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shd_ax2 ,  3,'shd2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shd_ax12,108,'shd12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shr_ax1 ,  3,'shr1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shr_ax2 ,  3,'shr2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shr_ax12,108,'shr12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shp_ax1 ,  3,'shp1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shp_ax2 ,  5,'shp2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shp_ax12,108,'shp12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shs_ax1 ,  3,'shs1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shs_ax2 ,  3,'shs2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shs_ax12,108,'shs12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  sha_ax1 ,  3,'sha1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  sha_ax2 ,  3,'sha2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  sha_ax12,108,'sha12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shg_ax1 ,  3,'shg1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shg_ax2 ,  3,'shg2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shg_ax12,108,'shg12'  ); call o_frame

!     call plotback; call cplot(nz,nr,  shh_ax1 ,  3,'shh1'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shh_ax2 ,  3,'shh2'   ); call o_frame
!     call plotback; call cplot(nz,nr,  shh_ax12,108,'shh12'  ); call o_frame

  endif

  call o_clswk()

  return
  end subroutine vortex_diagnose

!==================================================================================

  subroutine vortex_diagnose0(thil_ax, theta_ax, shw_ax, shv_ax, vtan_ax, &
                              vrad_ax,     w_ax, shc_ax, shd_ax,  shr_ax, &
                               shp_ax,   shs_ax, sha_ax, shg_ax,  shh_ax)

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
  real, intent(out) ::   shc_ax(nz,nr) ! cloud water water specific density (kg/kg)
  real, intent(out) ::   shd_ax(nz,nr) ! drizzle water specific density (kg/kg)
  real, intent(out) ::   shr_ax(nz,nr) ! rain water specific density (kg/kg)
  real, intent(out) ::   shp_ax(nz,nr) ! pristine ice water specific density (kg/kg)
  real, intent(out) ::   shs_ax(nz,nr) ! snow water specific density (kg/kg)
  real, intent(out) ::   sha_ax(nz,nr) ! aggregates water specific density (kg/kg)
  real, intent(out) ::   shg_ax(nz,nr) ! graupel water specific density (kg/kg)
  real, intent(out) ::   shh_ax(nz,nr) ! hail water specific density (kg/kg)

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: unxrad,unyrad,unzrad,unxtan,unytan,unztan
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan

  real :: rad,rrad,wrad1,wrad2
  real :: vtan_ax0,vrad_ax0

  integer :: irad,ir

  real :: weight_t(nz,nr)   ! weight array for T points
  real :: weight_u(nz,nr)   ! weight array for U points
  real :: weight_v(nz,nr)   ! weight array for V points

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
    shc_ax(:,:) = 0.
    shd_ax(:,:) = 0.
    shr_ax(:,:) = 0.
    shp_ax(:,:) = 0.
    shs_ax(:,:) = 0.
    sha_ax(:,:) = 0.
    shg_ax(:,:) = 0.
    shh_ax(:,:) = 0.

  weight_t(:,:) = 0.

  call psub()
!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------
  call qsub('W',iw)

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

          shc_ax(k,irad)   =   shc_ax(k,irad)   + wrad1 *  sh_c(k,iw)
          shc_ax(k,irad+1) =   shc_ax(k,irad+1) + wrad2 *  sh_c(k,iw)

!          shd_ax(k,irad)   =   shd_ax(k,irad)   + wrad1 *  sh_d(k,iw)
!          shd_ax(k,irad+1) =   shd_ax(k,irad+1) + wrad2 *  sh_d(k,iw)

!          shr_ax(k,irad)   =   shr_ax(k,irad)   + wrad1 *  sh_r(k,iw)
!          shr_ax(k,irad+1) =   shr_ax(k,irad+1) + wrad2 *  sh_r(k,iw)

          shp_ax(k,irad)   =   shp_ax(k,irad)   + wrad1 *  sh_p(k,iw)
          shp_ax(k,irad+1) =   shp_ax(k,irad+1) + wrad2 *  sh_p(k,iw)

!          shs_ax(k,irad)   =   shs_ax(k,irad)   + wrad1 *  sh_s(k,iw)
!          shs_ax(k,irad+1) =   shs_ax(k,irad+1) + wrad2 *  sh_s(k,iw)

!          sha_ax(k,irad)   =   sha_ax(k,irad)   + wrad1 *  sh_a(k,iw)
!          sha_ax(k,irad+1) =   sha_ax(k,irad+1) + wrad2 *  sh_a(k,iw)

!          shg_ax(k,irad)   =   shg_ax(k,irad)   + wrad1 *  sh_g(k,iw)
!          shg_ax(k,irad+1) =   shg_ax(k,irad+1) + wrad2 *  sh_g(k,iw)

!          shh_ax(k,irad)   =   shh_ax(k,irad)   + wrad1 *  sh_h(k,iw)
!          shh_ax(k,irad+1) =   shh_ax(k,irad+1) + wrad2 *  sh_h(k,iw)

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
             shc_ax(k,irad) =   shc_ax(k,irad+1)
!             shd_ax(k,irad) =   shd_ax(k,irad+1)
!             shr_ax(k,irad) =   shr_ax(k,irad+1)
!             shp_ax(k,irad) =   shp_ax(k,irad+1)
!             shs_ax(k,irad) =   shs_ax(k,irad+1)
!             sha_ax(k,irad) =   sha_ax(k,irad+1)
!             shg_ax(k,irad) =   shg_ax(k,irad+1)
!             shh_ax(k,irad) =   shh_ax(k,irad+1)

        else

            thil_ax(k,irad) =  thil_ax(k,irad) / weight_t(k,irad)
           theta_ax(k,irad) = theta_ax(k,irad) / weight_t(k,irad)
             shw_ax(k,irad) =   shw_ax(k,irad) / weight_t(k,irad)
             shv_ax(k,irad) =   shv_ax(k,irad) / weight_t(k,irad)
               w_ax(k,irad) =     w_ax(k,irad) / weight_t(k,irad)
             shc_ax(k,irad) =   shc_ax(k,irad) / weight_t(k,irad)
!             shd_ax(k,irad) =   shd_ax(k,irad) / weight_t(k,irad)
!             shr_ax(k,irad) =   shr_ax(k,irad) / weight_t(k,irad)
             shp_ax(k,irad) =   shp_ax(k,irad) / weight_t(k,irad)
!             shs_ax(k,irad) =   shs_ax(k,irad) / weight_t(k,irad)
!             sha_ax(k,irad) =   sha_ax(k,irad) / weight_t(k,irad)
!             shg_ax(k,irad) =   shg_ax(k,irad) / weight_t(k,irad)
!             shh_ax(k,irad) =   shh_ax(k,irad) / weight_t(k,irad)

        endif

     enddo

  enddo

! If meshtype = 1, diagnose UC field

  if (meshtype == 1) then

     weight_u(:,:) = 0.

!----------------------------------------------------------------------
     do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------

! Distance of this point from eye center

        rad = sqrt((xeu(iu)-xeh)**2 + (yeu(iu)-yeh)**2 + (zeu(iu)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

        if (rad >= radius_ax(nr) - 1.) cycle

! Determine interpolation point in radial dimension

        irad = 1
        do while (rad > radius_ax(irad+1))
           irad = irad + 1
        enddo

        wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
        wrad1 = 1. - wrad2
   
! Unit normal vector components from hurricane center to current IU point

        unxrad = (xeu(iu) - xeh) / rad
        unyrad = (yeu(iu) - yeh) / rad
        unzrad = (zeu(iu) - zeh) / rad

! Unit vector components in direction of tangential vortex wind

        unxtan = wnyh * unzrad - wnzh * unyrad
        unytan = wnzh * unxrad - wnxh * unzrad
        unztan = wnxh * unyrad - wnyh * unxrad

! Vertical loop over T levels

        do k = lpu(iu),nzz

           weight_u(k,irad)   = weight_u(k,irad)   + wrad1
           weight_u(k,irad+1) = weight_u(k,irad+1) + wrad2

! Diagnose axisymmetric component of model's own vortex

           vtan_ax0 = uc(k,iu) &
                    * (unx(iu) * unxtan + uny(iu) * unytan + unz(iu) * unztan) &
                    + vc(k,iu) &
                    * (vnx(iu) * unxtan + vny(iu) * unytan + vnz(iu) * unztan)

           vrad_ax0 = uc(k,iu) &
                    * (unx(iu) * unxrad + uny(iu) * unyrad + unz(iu) * unzrad) &
                    + vc(k,iu) &
                    * (vnx(iu) * unxrad + vny(iu) * unyrad + vnz(iu) * unzrad)

           vtan_ax(k,irad)   = vtan_ax(k,irad)   + wrad1 * vtan_ax0
           vtan_ax(k,irad+1) = vtan_ax(k,irad+1) + wrad2 * vtan_ax0

           vrad_ax(k,irad)   = vrad_ax(k,irad)   + wrad1 * vrad_ax0
           vrad_ax(k,irad+1) = vrad_ax(k,irad+1) + wrad2 * vrad_ax0

        enddo

     enddo

! Convert sums to averages

     do irad = nr,2,-1

        do k = 2,nzz

           if (weight_u(k,irad) < 1.e-6) then

              if (irad == nr) stop 'stop irad U'

              vtan_ax(k,irad) = vtan_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

              vrad_ax(k,irad) = vrad_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

           else

              vtan_ax(k,irad) = vtan_ax(k,irad) / weight_u(k,irad)
              vrad_ax(k,irad) = vrad_ax(k,irad) / weight_u(k,irad)

           endif

        enddo

     enddo

     vtan_ax(:,1) = 0.
     vrad_ax(:,1) = 0.

  else  ! If meshtype = 2, diagnose VC field

     weight_v(:,:) = 0.

!----------------------------------------------------------------------
     do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
 
! Distance of this point from eye center

        rad = sqrt((xev(iv)-xeh)**2 + (yev(iv)-yeh)**2 + (zev(iv)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

        if (rad >= radius_ax(nr) - 1.) cycle

! Determine interpolation point in radial dimension

        irad = 1
        do while (rad > radius_ax(irad+1))
           irad = irad + 1
        enddo

        wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
        wrad1 = 1. - wrad2
   
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

           weight_v(k,irad)   = weight_v(k,irad)   + wrad1
           weight_v(k,irad+1) = weight_v(k,irad+1) + wrad2

! Diagnose axisymmetric component of model's own vortex

           vtan_ax0 = vc(k,iv) &
                    * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan) &
                    + uc(k,iv) &
                    * (unx(iv) * vnxtan + uny(iv) * vnytan + unz(iv) * vnztan)

           vrad_ax0 = vc(k,iv) &
                    * (vnx(iv) * vnxrad + vny(iv) * vnyrad + vnz(iv) * vnzrad) &
                    + uc(k,iv) &
                    * (unx(iv) * vnxrad + uny(iv) * vnyrad + unz(iv) * vnzrad)

           vtan_ax(k,irad)   = vtan_ax(k,irad)   + wrad1 * vtan_ax0
           vtan_ax(k,irad+1) = vtan_ax(k,irad+1) + wrad2 * vtan_ax0

           vrad_ax(k,irad)   = vrad_ax(k,irad)   + wrad1 * vrad_ax0
           vrad_ax(k,irad+1) = vrad_ax(k,irad+1) + wrad2 * vrad_ax0

        enddo

     enddo

! Convert sums to averages

     do irad = nr,2,-1

        do k = 2,nzz

           if (weight_v(k,irad) < 1.e-6) then

              if (irad == nr) stop 'stop irad V'

              vtan_ax(k,irad) = vtan_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

              vrad_ax(k,irad) = vrad_ax(k,irad+1) &
                              * radius_ax(irad) / radius_ax(irad+1)

           else

              vtan_ax(k,irad) = vtan_ax(k,irad) / weight_v(k,irad)
              vrad_ax(k,irad) = vrad_ax(k,irad) / weight_v(k,irad)

           endif

        enddo

     enddo

     vtan_ax(:,1) = 0.
     vrad_ax(:,1) = 0.

  endif

  return
  end subroutine vortex_diagnose0

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
  use massflux,  only: diagnose_uc, diagnose_vc

  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_u, mpi_recv_u, &
                          mpi_send_v, mpi_recv_v 

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
!  real,    parameter :: rexpon_delc = 1.0
  real,    parameter :: rexpon_delc = 0.5

  real :: vpert_rad, vpert_vert

  real :: rad_eyw, rad_env
  real :: circ_eyw, circ_env, dcirc, circ0
  real :: wrad1, wrad2

  real :: del_circ(nr), circ(nr)

  real ::   theta_totw(nz,nr)
  real ::     vtan_tot(nz,nr)

  integer :: irad
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

! Read diagnosed hurricane initial tangential wind from file

  open(32,file='frances_diagnosis_1sept00z',status='old',form='unformatted')
  read(32) thil_ax1, theta_ax1, shw_ax1, shv_ax1, vtan_ax1 !, &
!          vrad_ax1,     w_ax1, shc_ax1, shd_ax1,  shr_ax1, &
!           shp_ax1,   shs_ax1, sha_ax1, shg_ax1,  shh_ax1
  close(32)

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

  call psub()
!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------
  call qsub('W',iw)

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

!if (iw == 9485) then
!   print*, ' '
!   print*, 'x1 ',nzz,k,thil(k,iw),theta_pert0
!   print*, 'x2 ',theta(k,iw),sh_w(k,iw),sh_v(k,iw)
!   print*, 'x3 ',theta_tot0,theta_tot1,wrad1,wrad2
!   print*, 'x4 ',irad,theta_totw(k,irad),theta_totw(k-1,irad)
!   print*, 'x5 ',theta_totw(k,irad+1),theta_totw(k-1,irad+1)
!   print*, 'x6 ',theta_ax1(k,irad),theta_ax1(k,irad+1)
!endif

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
                 + gravo2 * (rho(k+1,iw) * dzt(k+1) + rho(k,iw) * dzt(k))
           press(k,iw) = .05 * press(k,iw) + .95 * pkhyd

        enddo
     enddo

!!     write(6,'(a,i8,2f10.0)') 'press ',iw,press(2,iw),rad

  enddo
  call rsub('W_frances_a',7)

! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_w('I')  ! Send W group
     call mpi_recv_w('I')  ! Recv W group
  endif

! If meshtype = 1, initialize UC field

  if (meshtype == 1) then

     call psub()
!----------------------------------------------------------------------
     do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
     call qsub('U',iu)

! Distance of this point from eye center

        rad = sqrt((xeu(iu)-xeh)**2 + (yeu(iu)-yeh)**2 + (zeu(iu)-zeh)**2)

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

! Unit normal vector components from hurricane center to current IU point

        unxrad = (xeu(iu) - xeh) / rad
        unyrad = (yeu(iu) - yeh) / rad
        unzrad = (zeu(iu) - zeh) / rad

! Unit vector components in direction of tangential vortex wind

        unxtan = wnyh * unzrad - wnzh * unyrad
        unytan = wnzh * unxrad - wnxh * unzrad
        unztan = wnxh * unyrad - wnyh * unxrad

! Vertical loop over T levels

        do k = lpu(iu),nzz

! Skip hurricane assimilation for all points above limiting height

           if (zt(k) >= zmax_delv) exit

! Vertical dependence of tangential velocity perturbation

           vpert_vert = 1. - (zt(k) / zmax_delv) ** zexpon_delv

! Combine vertical and radial parts of perturbation

           vtan_pert0 = vpert_rad * vpert_vert

! Add enhanced vortex to winds interpolated from NCEP reanalysis

           uc(k,iu) = uc(k,iu) + vtan_pert0  &
              * (unx(iu) * unxtan + uny(iu) * unytan + unz(iu) * unztan)

           umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
           ump(k,iu) = umc(k,iu)

        enddo

     enddo
     call rsub('U_frances_a',7)

! If using MPI, perform parallel send/recv

     if (iparallel == 1) then
        call mpi_send_u('I')
        call mpi_recv_u('I')
     endif

! Re-diagnose V velocity component

     call diagnose_vc()

  else  ! If meshtype = 2, initialize VC field

     call psub()
!----------------------------------------------------------------------
     do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
     call qsub('V',iv)
 
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
           vmp(k,iv) = vmc(k,iv)

        enddo

     enddo
     call rsub('V_frances_a',7)

! If using MPI, perform parallel send/recv

     if (iparallel == 1) then
        call mpi_send_v('I')
        call mpi_recv_v('I')
     endif

! Re-diagnose U velocity component

     call diagnose_uc()

  endif

  return
  end subroutine hurricane_init1B

!=======================================================================================

  subroutine hurricane_init2B()

! This subroutine adds perturbation fields to a model initial state that act to
! capture the full intensity of a tropical cyclone that is poorly resolved in
! the standard analysis used for the initial fields (e.g., the GFS analysis).

! The added perturbation fields currently consist of potential temperature, 
! water vapor specific density, and cyclonic winds.  Additional fields may be added
! in the future, including microphysics condensate species and radial-vertical
! circulation.

! The added perturbation fields currently are axisymmetric, defined as a function
! of radius and height only.  They can be made asymmetric with some code 
! development, for example by including azimuthal wavenumbers.

! The perturbation fields that are added in this subroutine are generated 
! by a dynamic initialization cycle, performed by previous integrations of OLAM.
! The cycle should be repeated whenever major changes are made to the model
! configuration, such as a change in grid resolution, or for a new vortex
! initialization time and location.

  use mem_ijtabs
  use mem_basic
  use mem_micro
  use misc_coms
  use mem_grid
  use consts_coms
  use massflux,  only: diagnose_uc, diagnose_vc

  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_u, mpi_recv_u, &
                          mpi_send_v, mpi_recv_v 

  implicit none

  integer :: iw,i,j,k,ngr,iu,iv,iw1,iw2,kbc,iter
  real :: pkhyd,exner,deltheta
  real :: zeh,reh,xeh,yeh
  real :: wnxh,wnyh,wnzh
  real :: unxrad,unyrad,unzrad,unxtan,unytan,unztan
  real :: vnxrad,vnyrad,vnzrad,vnxtan,vnytan,vnztan

  real :: rad,rrad,wrad1,wrad2

  real :: thil_ax, theta_ax, shw_ax, shv_ax, vtan_ax
  real :: vrad_ax,     w_ax, shc_ax, shd_ax,  shr_ax
  real ::  shp_ax,   shs_ax, sha_ax, shg_ax,  shh_ax

  integer :: irad,iz,ir

  character*10 :: cstr
  integer :: istr
  real :: rstr

  real :: rad_pert(nr)

!  11 format(/,a,28f6.0,/)
!  12 format(i3,f8.0,a,28f6.1)

  11 format(/,15x  ,a,28f6.0,/)
  12 format(i3,f8.0,a,28f6.1)

!  11 format(/,3x,f8.0,a,30f6.0)
!  12 format(/,i3,f8.0,a,30f6.1)

! Find "earth" coordinates of hurricane center

  zeh = erad * sin(hlat * pio180)
  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

! Thermodynamic fields

  open(31,file='frances_cyc1_pert3',status='old',form='formatted')

  read(31) thil_ax12, theta_ax12, shw_ax12, shv_ax12, vtan_ax12, &
           vrad_ax12,     w_ax12, shc_ax12, shd_ax12,  shr_ax12, &
            shp_ax12,   shs_ax12, sha_ax12, shg_ax12,  shh_ax12

  close(31)

  call psub()
!----------------------------------------------------------------------
  do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)  ! jend(1) = hardwired for mrl 1
!----------------------------------------------------------------------
  call qsub('W',iw)

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
   
! Vertical loop over T levels

     do k = lpw(iw),nzz

! Interpolate perturbation table values to current grid cell

         thil_ax = wrad1 *  thil_ax12(k,irad) + wrad2 *  thil_ax12(k,irad+1)
        theta_ax = wrad1 * theta_ax12(k,irad) + wrad2 * theta_ax12(k,irad+1)
          shw_ax = wrad1 *   shw_ax12(k,irad) + wrad2 *   shw_ax12(k,irad+1)
          shv_ax = wrad1 *   shv_ax12(k,irad) + wrad2 *   shv_ax12(k,irad+1)
            w_ax = wrad1 *     w_ax12(k,irad) + wrad2 *     w_ax12(k,irad+1)
          shc_ax = wrad1 *   shc_ax12(k,irad) + wrad2 *   shc_ax12(k,irad+1)
          shd_ax = wrad1 *   shd_ax12(k,irad) + wrad2 *   shd_ax12(k,irad+1)
          shr_ax = wrad1 *   shr_ax12(k,irad) + wrad2 *   shr_ax12(k,irad+1)
          shp_ax = wrad1 *   shp_ax12(k,irad) + wrad2 *   shp_ax12(k,irad+1)
          shs_ax = wrad1 *   shs_ax12(k,irad) + wrad2 *   shs_ax12(k,irad+1)
          sha_ax = wrad1 *   sha_ax12(k,irad) + wrad2 *   sha_ax12(k,irad+1)
          shg_ax = wrad1 *   shg_ax12(k,irad) + wrad2 *   shg_ax12(k,irad+1)
          shh_ax = wrad1 *   shh_ax12(k,irad) + wrad2 *   shh_ax12(k,irad+1)

!---------------------------------------------------------------------
!        if (thil_ax > 0.) thil_ax = thil_ax * 1.2
!---------------------------------------------------------------------

! Add theta and shw perturbations to model fields

         thil(k,iw) =  thil(k,iw) +  thil_ax
        theta(k,iw) = theta(k,iw) + theta_ax
         sh_w(k,iw) =  sh_w(k,iw) +   shw_ax
         sh_v(k,iw) =  sh_v(k,iw) +   shv_ax
           wc(k,iw) =    wc(k,iw) +     w_ax
         sh_c(k,iw) =  sh_c(k,iw) +   shc_ax
         sh_d(k,iw) =  sh_d(k,iw) +   shd_ax
         sh_r(k,iw) =  sh_r(k,iw) +   shr_ax
         sh_p(k,iw) =  sh_p(k,iw) +   shp_ax
         sh_s(k,iw) =  sh_s(k,iw) +   shs_ax
         sh_a(k,iw) =  sh_a(k,iw) +   sha_ax
         sh_g(k,iw) =  sh_g(k,iw) +   shg_ax
         sh_h(k,iw) =  sh_h(k,iw) +   shh_ax

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
              + gravo2 * (rho(k+1,iw) * dzt(k+1) + rho(k,iw) * dzt(k))
           press(k,iw) = .05 * press(k,iw) + .95 * pkhyd

        enddo
     enddo

! Vertical loop over T levels

     do k = lpw(iw),nzz
        wmc(k,iw) = wc(k,iw) * .5 * (rho(k,iw) + rho(k+1,iw))
     enddo

  enddo
  call rsub('W_frances_a',7)

! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_w('I')  ! Send W group
     call mpi_recv_w('I')  ! Recv W group
  endif

! If meshtype = 1, initialize UC field

  if (meshtype == 1) then

     call psub()
!----------------------------------------------------------------------
     do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
     call qsub('U',iu)

! Distance of this point from eye center

        rad = sqrt((xeu(iu)-xeh)**2 + (yeu(iu)-yeh)**2 + (zeu(iu)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

        if (rad >= radius_ax(nr) - 1.) cycle

! Determine interpolation point in radial dimension

        irad = 1
        do while (rad > radius_ax(irad+1))
           irad = irad + 1
        enddo

        wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
        wrad1 = 1. - wrad2
   
! Unit normal vector components from hurricane center to current IU point

        unxrad = (xeu(iu) - xeh) / rad
        unyrad = (yeu(iu) - yeh) / rad
        unzrad = (zeu(iu) - zeh) / rad

! Unit vector components in direction of tangential vortex wind

        unxtan = wnyh * unzrad - wnzh * unyrad
        unytan = wnzh * unxrad - wnxh * unzrad
        unztan = wnxh * unyrad - wnyh * unxrad

! Vertical loop over T levels

        do k = lpu(iu),nzz

! Interpolate perturbation table values to current grid cell

           vtan_ax = wrad1 * vtan_ax12(k,irad  )  &
                   + wrad2 * vtan_ax12(k,irad+1)

           vrad_ax = wrad1 * vrad_ax12(k,irad  )  &
                   + wrad2 * vrad_ax12(k,irad+1)

!---------------------------------------------------------------------
!           if (vtan_pert0 > 0.) vtan_pert0 = vtan_pert0 * 1.6
!---------------------------------------------------------------------

! Add enhanced vortex to winds interpolated from NCEP reanalysis

           uc(k,iu) = uc(k,iu) &
              + vtan_ax &
              * (unx(iu) * unxtan + uny(iu) * unytan + unz(iu) * unztan) &
              + vrad_ax &
              * (unx(iu) * unxrad + uny(iu) * unyrad + unz(iu) * unzrad)

           umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
           ump(k,iu) = umc(k,iu)

        enddo

     enddo
     call rsub('U_frances_a',7)

! If using MPI, perform parallel send/recv

     if (iparallel == 1) then
        call mpi_send_u('I')
        call mpi_recv_u('I')
     endif

! Re-diagnose V velocity component

     call diagnose_vc()

  else  ! If meshtype = 2, initialize VC field

     call psub()
!----------------------------------------------------------------------
     do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)  ! jend(1) = hardwired for mrl 1
        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
     call qsub('V',iv)
 
! Distance of this point from eye center

        rad = sqrt((xev(iv)-xeh)**2 + (yev(iv)-yeh)**2 + (zev(iv)-zeh)**2)

! Skip hurricane assimilation for all points outside specified radius

        if (rad >= radius_ax(nr) - 1.) cycle

! Determine interpolation point in radial dimension

        irad = 1
        do while (rad > radius_ax(irad+1))
           irad = irad + 1
        enddo

        wrad2 = (rad - radius_ax(irad)) / (radius_ax(irad+1) - radius_ax(irad))
        wrad1 = 1. - wrad2
   
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

! Interpolate perturbation table values to current grid cell

           vtan_ax = wrad1 * vtan_ax12(k,irad  ) &
                   + wrad2 * vtan_ax12(k,irad+1)

           vrad_ax = wrad1 * vrad_ax12(k,irad  ) &
                   + wrad2 * vrad_ax12(k,irad+1)

!---------------------------------------------------------------------
!           if (vtan_pert0 > 0.) vtan_pert0 = vtan_pert0 * 1.2
!---------------------------------------------------------------------

! Add enhanced vortex to winds interpolated from NCEP reanalysis

           vc(k,iv) = vc(k,iv) &
              + vtan_ax &
              * (vnx(iv) * vnxtan + vny(iv) * vnytan + vnz(iv) * vnztan) &
              + vrad_ax &
              * (vnx(iv) * vnxrad + vny(iv) * vnyrad + vnz(iv) * vnzrad)

           vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
           vmp(k,iv) = vmc(k,iv)
        enddo

     enddo
     call rsub('V_frances_a',7)

! If using MPI, perform parallel send/recv

     if (iparallel == 1) then
        call mpi_send_v('I')
        call mpi_recv_v('I')
     endif

! Re-diagnose U velocity component

     call diagnose_uc()

  endif

  return
  end subroutine hurricane_init2B

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

  integer :: k,km,i,ival,icolor,ibox,ln
  real :: bsize,yinc

  data iasf / 18*1 /
  character*8 :: number,numbr
  character*60 :: title

! turn off the clipping indicator.

  call gsclip (0)

! set all the gks aspect source flags to "individual".

  call gsasf (iasf)

! force solid fill.

  call gsfais (1)

  call set(.01,.87,.01,.96,0.,radius_ax(nr),0.,zm(nzz),1)

  op%xmin = 0.
  op%xmax = radius_ax(nr)
  op%ymin = 0.
  op%ymax = zm(nzz)

  call sfseti ('type of fill',0)

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

  call sflush

! draw a color bar for the plot.

  call set (0.,1.,0.,1.,0.,1.,0.,1.,1)

  call gsplci(8)
  call gstxci(8)

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
        call plchlq (hcpn(2)+.01,vcpn(3),numbr(1:ln),bsize,0.,-1.)
     endif

     call sfsgfa (hcpn,vcpn,4,dst,6,ind,8,clrtab(itab)%ipal(ibox))
!     call fillpolyg (4,hcpn,vcpn,ipal(ibox))

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

