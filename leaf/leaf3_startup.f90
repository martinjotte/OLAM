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
subroutine leaf3_startup()

use leaf_coms, only: nzg, nzs, landusefile, mml, mul, mwl, &
                     alloc_leafcol, isoildepthflg, ndviflg
use mem_leaf,  only: alloc_leaf, filltab_ED, filltab_leaf, land
use misc_coms, only: io6, runtype

implicit none

integer :: iwl

! Subroutine LEAF3_STARTUP allocates and initializes some leaf3 arrays.

! THIS SUBROUTINE DOES NOT INITIALIZE soil/veg/canopy temperature and moisture
! values, which depend on atmospheric conditions.

!-------------------------------------------------------------------------------
! STEP 1: Allocate leafcol arrays
!-------------------------------------------------------------------------------

call alloc_leafcol()

!-------------------------------------------------------------------------------
! STEP 2: Initialize inherent soil and vegetation properties
!-------------------------------------------------------------------------------

call sfcdata()

call load_ed_ecosystem_params()

!-------------------------------------------------------------------------------
! STEP 3: Allocate main LEAF arrays
!-------------------------------------------------------------------------------

call alloc_leaf(mwl, nzg, nzs)
call filltab_leaf()

!-------------------------------------------------------------------------------
! STEP 4: Read soil depth info and initialize ED sites
!-------------------------------------------------------------------------------

if (isoildepthflg == 1) call leaf_soil_depth_read()
call initialize_ed_sites()

! Here, specify what land-grid variables to write out

call filltab_ED()

!-------------------------------------------------------------------------------
! STEP 5: Fill ndvi values
!-------------------------------------------------------------------------------

if (ndviflg == 2) then

! Default initialization of NDVI

   do iwl = 2,mwl
      land%veg_ndvip(iwl) = .5
      land%veg_ndvic(iwl) = .5
   enddo

elseif (runtype /= 'PLOTONLY' .and. runtype /= 'PARCOMBINE') then

! Make inventory of ndvi files and initialize current/past time
! Not needed for a PLOTONLY run.

   call ndvi_database_read(0)

! Initialize future ndvi file time

   call ndvi_database_read(1)

endif

return
end subroutine leaf3_startup

!==========================================================================

subroutine sfcdata()

use leaf_coms, only: slz, nstyp, nvtyp, nzg, refdepth, slcons0, fhydraul,  &
                     slcons1, slpots, slmsts, slbs, slcons, slcpd, slden,  &
                     xsand, xclay, xorgan, xrobulk, soilcond0, soilcond1,  &
                     soilcond2, emisg, soilcp, albv_green, albv_brown,     &
                     emisv, sr_max, tai_max, sai, veg_clump, veg_frac,     &
                     veg_ht, dead_frac, rcmin, glai_max, dfpardsr,         &
                     fpar_max, fpar_min, sr_min, root, kroot,              &
                     dt_leaf, dslz, dslzo2, dslzi,                         &
                     dslzidt, slzt, dslzt, dslzti, dslztidt
use misc_coms, only: io6

implicit none

integer :: k
integer :: nnn

real :: soilparms(7,nstyp)
real :: soilparms2(7,nstyp)

real :: bioparms(12,0:nvtyp)

!  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke,
!                        1981; Pielke, 1984; Tremback & Kessler, 1985)
!
!  slpots  - sat moisture potential [m]
!  slmsts  - sat volumetric moisture content (soil porosity) [m^3_wat/m^3_tot]
!  slbs    - b exponent [dimensionless]
!  slcons  - saturation soil hydraulic conductivity [m/s]
!  slcons0 - surface value for slcons [m/s]
!  slcpd   - dry soil volumetric heat capacity [J/(m^3 K)]
!  slden   - dry soil density [kg/m3]

data soilparms/  &
!-----------------------------------------------------------------------------
!slpots        slbs          slcons0         slden       USDA SOIL CLASS
!      slmsts         slcons          slcpd              # AND NAME
!-----------------------------------------------------------------------------
-.121, .395,  4.05, .18e-3, .50e-3, 1465.e3, 1600.,  & !  1 sand
-.090, .410,  4.38, .16e-3, .60e-3, 1407.e3, 1600.,  & !  2 loamy sand
-.218, .435,  4.9 , .34e-4, .77e-3, 1344.e3, 1600.,  & !  3 sandy loam
-.786, .485,  5.3 , .72e-5, .11e-4, 1273.e3, 1600.,  & !  4 silt loam
-.478, .451,  5.39, .69e-5, .22e-2, 1214.e3, 1600.,  & !  5 loam
-.299, .420,  7.12, .63e-5, .15e-2, 1177.e3, 1600.,  & !  6 sandy clay loam
-.356, .477,  7.75, .17e-5, .11e-3, 1319.e3, 1600.,  & !  7 silty clay loam
-.630, .476,  8.52, .24e-5, .22e-2, 1227.e3, 1600.,  & !  8 clay loam
-.153, .426, 10.4 , .22e-5, .22e-5, 1177.e3, 1600.,  & !  9 sandy clay
-.490, .492, 10.4 , .10e-5, .10e-5, 1151.e3, 1600.,  & ! 10 silty clay
-.405, .482, 11.4 , .13e-5, .13e-5, 1088.e3, 1600.,  & ! 11 clay 
-.356, .863,  7.75, .80e-5, .80e-5,  874.e3,  300./    ! 12 peat

data soilparms2/  &
!-----------------------------------------------------------------------------
!xsand      xorgan        soilcond0    soilcond2    USDA SOIL CLASS
!      xclay      xrobulk       soilcond1           # AND NAME
!-----------------------------------------------------------------------------
.97,  .03,  .00,  1200.,  .30,  4.80,  -2.70,  & !  1 sand
.92,  .07,  .01,  1250.,  .30,  4.66,  -2.60,  & !  2 loamy sand
.80,  .18,  .02,  1300.,  .29,  4.27,  -2.31,  & !  3 sandy loam
.57,  .40,  .03,  1400.,  .27,  3.47,  -1.74,  & !  4 silt loam
.60,  .35,  .05,  1350.,  .28,  3.63,  -1.85,  & !  5 loam
.65,  .31,  .04,  1350.,  .28,  3.78,  -1.96,  & !  6 sandy clay loam
.35,  .59,  .06,  1500.,  .26,  2.73,  -1.20,  & !  7 silty clay loam
.48,  .45,  .07,  1450.,  .27,  3.23,  -1.56,  & !  8 clay loam
.50,  .42,  .08,  1450.,  .27,  3.32,  -1.63,  & !  9 sandy clay
.30,  .61,  .09,  1650.,  .25,  2.58,  -1.09,  & ! 10 silty clay
.25,  .65,  .10,  1700.,  .25,  2.40,  -0.96,  & ! 11 clay
.20,  .20,  .60,   500.,  .06,  0.46,   0.00/    ! 12 peat

!         LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

data bioparms/  &
!-----------------------------------------------------------------------------
!albv_green     sr_max         veg_clump       rootdep             LEAF-3 CLASS #
!     albv_brown     tai_max        veg_frac        dead_frac      AND DESCRIPTION
!          emisv          sai            veg_ht         rcmin
!-----------------------------------------------------------------------------
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  0  Ocean
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  1  Lakes, rivers, streams
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  2  Ice cap/glacier
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  3  Desert, bare soil
 .14, .24, .97, 5.4, 8.0, 1.0, 1.0, .80, 20.0, 5.0, .0, 500., & !  4  Evergreen needleleaf tree
 .14, .24, .95, 5.4, 8.0, 1.0, 1.0, .80, 22.0, 5.0, .0, 500., & !  5  Deciduous needleleaf tree
 .20, .24, .95, 6.2, 7.0, 1.0,  .0, .80, 22.0, 5.0, .0, 500., & !  6  Deciduous broadleaf tree
 .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 5.0, .0, 500., & !  7  Evergreen broadleaf tree
 .21, .43, .96, 5.1, 4.0, 1.0,  .0, .75,   .3,  .7, .7, 100., & !  8  Short grass
 .24, .43, .96, 5.1, 5.0, 1.0,  .0, .80,  1.2, 1.0, .7, 100., & !  9  Tall grass
 .24, .24, .96, 5.1, 1.0,  .2, 1.0, .20,   .7, 1.0, .0, 500., & ! 10  Semi-desert
 .20, .24, .95, 5.1, 4.5,  .5, 1.0, .60,   .2, 1.0, .0,  50., & ! 11  Tundra
 .14, .24, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 12  Evergreen shrub
 .20, .28, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 13  Deciduous shrub
 .16, .24, .96, 6.2, 7.0, 1.0,  .5, .80, 22.0, 5.0, .0, 500., & ! 14  Mixed woodland
 .22, .40, .95, 5.1, 5.0,  .5,  .0, .85,  1.0, 1.0, .0, 100., & ! 15  Crop/mixed farming, C3 grassland
 .18, .40, .95, 5.1, 5.0,  .5,  .0, .80,  1.1, 1.0, .0, 500., & ! 16  Irrigated crop
 .12, .43, .98, 5.1, 7.0, 1.0,  .0, .80,  1.6, 1.0, .0, 500., & ! 17  Bog or marsh
 .20, .36, .96, 5.1, 6.0, 1.0,  .0, .80,  7.0, 3.5, .0, 100., & ! 18  Wooded grassland 
 .20, .36, .90, 5.1, 3.6, 1.0,  .0, .74,  6.0, 2.5, .0, 500., & ! 19  Urban and built up
 .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 5.0, .0, 500./   ! 20  Wetland evergreen broadleaf tree
!.16, .24, .96, 5.1, 2.0, 1.5, 1.0, .10, 20.0, 1.5, .0, 500./   ! 21  Very urban

! Soil vertical grid spacing arrays (some with timestep info)

slz(nzg+1) = 0.

do k = 1,nzg
   dslz   (k) = slz(k+1) - slz(k)
   dslzo2 (k) = .5 * dslz(k)
   dslzi  (k) = 1. / dslz(k)
   dslzidt(k) = dslzi(k) * dt_leaf
   slzt   (k) = .5 * (slz(k) + slz(k+1))
enddo

do k = 2,nzg
   dslzt   (k) = slzt(k) - slzt(k-1)
   dslzti  (k) = 1. / dslzt(k)
   dslztidt(k) = dslzti(k) * dt_leaf
enddo

! Soil constants

refdepth = -2.0

do nnn = 1,nstyp
   slcons0 (nnn) = soilparms(5,nnn)
   fhydraul(nnn) = log (soilparms(4,nnn) / soilparms(5,nnn)) / refdepth

   do k = 1,nzg
      slcons1(k,nnn) = soilparms(4,nnn)     ! ORIGINAL form - const with depth
!     slcons1(k,nnn) = soilparms(5,nnn)  &  ! TOPMODEL form - large at surface
!        * exp(slz(k) * fhydraul(nnn))      !    and exp decrease with depth
   enddo

   slpots   (nnn) = soilparms(1,nnn)
   slmsts   (nnn) = soilparms(2,nnn)
   slbs     (nnn) = soilparms(3,nnn)
   slcons   (nnn) = soilparms(4,nnn)
   slcpd    (nnn) = soilparms(6,nnn)
   slden    (nnn) = soilparms(7,nnn)

   xsand    (nnn) = soilparms2(1,nnn)
   xclay    (nnn) = soilparms2(2,nnn)
   xorgan   (nnn) = soilparms2(3,nnn)
   xrobulk  (nnn) = soilparms2(4,nnn)
   soilcond0(nnn) = soilparms2(5,nnn)
   soilcond1(nnn) = soilparms2(6,nnn)
   soilcond2(nnn) = soilparms2(7,nnn)

   emisg (nnn) = .98
   soilcp(nnn) = 0.1 - 0.07 * xsand(nnn)
enddo

do nnn = 1,nvtyp
   albv_green(nnn) = bioparms(1,nnn)
   albv_brown(nnn) = bioparms(2,nnn)
   emisv     (nnn) = bioparms(3,nnn)
   sr_max    (nnn) = bioparms(4,nnn)
   tai_max   (nnn) = bioparms(5,nnn)
   sai       (nnn) = bioparms(6,nnn)
   veg_clump (nnn) = bioparms(7,nnn)
   veg_frac  (nnn) = bioparms(8,nnn)
   veg_ht    (nnn) = bioparms(9,nnn)
   dead_frac (nnn) = bioparms(11,nnn)
   rcmin     (nnn) = bioparms(12,nnn)
   glai_max  (nnn) = tai_max(nnn) - sai(nnn)
   dfpardsr  (nnn) = (fpar_max - fpar_min) / (sr_max(nnn) - sr_min)
   kroot     (nnn) = nzg

   do k = nzg-1,1,-1
      if (slz(k+1) > -bioparms(10,nnn)) kroot(nnn) = k
   enddo

!  root    (1,nnn) = 0.    ! not used
enddo

return
end subroutine sfcdata

!=======================================================================

subroutine leaf_soil_depth_read()

  use leaf_coms,   only: mwl, nzg, slz, soildepth_db
  use mem_leaf,    only: land
  use consts_coms, only: piu180, erad
  use misc_coms,   only: io6
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info

  implicit none

  integer :: nio, njo, nperdeg
  integer :: io1, io2, jo1, jo2
  integer :: iwl
  integer :: layer_index
  integer :: k, ndims, idims(2)
  integer :: nlon, nlat
  real    :: wio1, wio2, wjo1, wjo2
  real    :: rio, rjo
  real    :: glat, glon
  real    :: soild
  real    :: offpix
  logical :: l1

  real, allocatable :: soil_depth(:,:)

  inquire(file=soildepth_db, exist=l1)
  if(.not.l1)then
     write(io6,*) 'You have ISOILDEPTHFLG = 1, which means you read a file'
     write(io6,*) 'containing soil depths.  However, the '
     write(io6,*) 'file you specified for SOILDEPTH_DB,'
     write(io6,*)
     write(io6,*) trim(soildepth_db)
     write(io6,*)
     write(io6,*) 'does not exist.'
     stop
  endif

  write(io6,*) 'Reading soil depth database '//trim(soildepth_db)
  call shdf5_open(soildepth_db,'R')

  call shdf5_info('soildepth',ndims,idims)
  nlon = idims(1)
  nlat = idims(2)
  
  allocate(soil_depth(nlon,nlat))

  call shdf5_irec(ndims,idims,'soildepth',rvara=soil_depth)
  call shdf5_close()

  offpix = 0.
  nperdeg = nlon/360
  if (mod(nlon,nperdeg) == 2) offpix = .5

  ! Loop over land sites

  do iwl = 2, mwl
     glon = max(-179.999,min(179.999,land%glonw(iwl)))
  
     rio = 1. + (glon            + 180.) * nperdeg + offpix
     rjo = 1. + (land%glatw(iwl) +  90.) * nperdeg + offpix

     io1 = int(rio)
     jo1 = int(rjo)
         
     wio2 = rio - float(io1)
     wjo2 = rjo - float(jo1)
           
     wio1 = 1. - wio2
     wjo1 = 1. - wjo2

     io2 = io1 + 1
     jo2 = jo1 + 1

     soild = &
          wio1 * (wjo1 * soil_depth(io1,jo1) + wjo2 * soil_depth(io1,jo2)) + &
          wio2 * (wjo1 * soil_depth(io2,jo1) + wjo2 * soil_depth(io2,jo2))

     ! We require at least TWO soil layers
     layer_index = nzg - 1
     do k = nzg-2,1,-1
        if (slz(k+1) > -soild) layer_index = k
     enddo

     land%lsl(iwl) = layer_index
  enddo

  deallocate(soil_depth)

  return
end subroutine leaf_soil_depth_read
