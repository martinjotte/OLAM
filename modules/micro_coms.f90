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
Module micro_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

!  The product [(nthz-1)  * dthz ] must equal 25.0
!  The product [(nrhhz-1) * drhhz] must equal 0.18
!  The product [(ntc-1)   * dtc  ] must equal 20.0
!  The product [(ndnc-1)  * ddnc ] must equal 20.e-6

integer, parameter :: maxgrds0 = 20  ! used for homog freezing table - change later

integer, parameter :: &
    ncat   = 8    & ! # of hydrometeor categories without distinguishing ice habits
   ,nhcat  = 16   & ! # of hydrometeor categories including ice habits

   ,nthz   = 26   & ! # of temp values spanning haze nucleation table
   ,nrhhz  = 10   & ! # of R.H. values spanning haze nucleation table
   ,ngam   = 3000 & ! # of values in incomplete gamma function tables
   ,ninc   = 201  & ! # of liquid fraction values spanning shed and melt tables
   ,ndns   = 15   & ! # of diameter values spanning shed table
   ,ntc    = 21   & ! # of temperature values spanning homogeneous freezing table
   ,ndnc   = 11   & ! # of diameter values spanning homogeneous freezing table

   ,neff   = 8    & ! 2nd array dim for eff (# of types of coalescence efficiency) 

   ,npairc = 106  & ! # of pairs of species in number collection table
   ,npairx = 104  & ! # of pairs of species in x mass collection table
   ,npairy = 55   & ! # of pairs of species in y mass collection table
   ,nembc  = 20     ! # of diam values spanning main collection table (in 2 dims)

real, parameter :: & 
    dtc   = 1.     & ! temperature increment in homogeneous nucleation table
   ,ddnc  = 2.e-6  & ! characteristic diam increment in homogeneous nucleation table
   ,dthz  = 1.     & ! temperature increment in haze nucleation table
   ,drhhz = .02      ! R.H. increment in haze nucleation table

! dn increment in incomplete gamma tables, and its inverse

real, parameter :: ddn_ngam  = 5.e-6
real, parameter :: ddni_ngam = 1. / ddn_ngam

integer :: miclevel  ! microphysics complexity level read from namelist
integer :: icloud ! cloud control flag read from namelist
integer :: idriz  ! drizzle control flag read from namelist
integer :: irain  ! rain control flag read from namelist
integer :: ipris  ! pristine ice control flag read from namelist
integer :: isnow  ! snow control flag read from namelist
integer :: iaggr  ! aggregates control flag read from namelist
integer :: igraup ! graupel control flag read from namelist
integer :: ihail  ! hail control flag read from namelist
   
integer :: iccn   ! CCN flag read from namelist
integer :: igccn  ! GCCN flag read from namelist
integer :: iifn   ! IFN flag read from namelist

integer :: mza0  ! micphys copy of number of model vertical levels

integer :: jnmb(ncat)  ! category control parameter for each category

! table to map lhcat to lcat
integer :: lcat_lhcat(nhcat) = (/1,2,3,4,5,6,7,8,3,3,3,3,4,4,4,4/)

integer :: jhabtab(31,100,2)

real :: rparm ! rain parameter read from namelist
real :: sparm ! snow parameter read from namelist
real :: aparm ! aggregates parameter read from namelist
real :: gparm ! graupel parameter read from namelist
real :: hparm ! hail parameter read from namelist

real :: ccnparm  ! CCN concentration [#/kg] read from namelist
real :: gccnparm ! GCCN concentration [#/kg] read from namelist
real :: ifnparm  ! IFN concentration [#/kg] read from namelist

real :: rictmin ! minimum diameter-index for collection table lookup
real :: rictmax ! maximum diameter-index for collection table lookup
real :: dps     ! cutoff diameter (microns) between pristine ice and snow
real :: dps2    ! square of dps

real :: sedtime0 ! min generalized timestep considered in sedimentation table
real :: sedtime1 ! max generalized timestep considered in sedimentation table

real :: dmb0     (ncat) ! minimum mean-mass diameter for each category
real :: dmb1     (ncat) ! maximum mean-mass diameter for each category
real :: emb0     (ncat) ! minimum mean mass for each category
real :: emb1     (ncat) ! maximum mean mass for each category
real :: emb0i    (ncat) ! 1 / (minimum mean mass for each category)
real :: emb1i    (ncat) ! 1 / (maximum mean mass for each category)
real :: emb2     (ncat) ! mean hydrometeor mass for each category for case when jnmb(lcat) = 2
real :: gnu      (ncat) ! gamma distribution width parameter for each category
real :: emb0log  (ncat) ! log of emb0
real :: emb1log  (ncat) ! log of emb1
real :: dict     (ncat) ! log-mass increment of each category in collection tables
real :: rxmin    (ncat) ! minimum bulk density for each category
   
real :: shapefac (nhcat) ! shape factor for each category
real :: cfmas    (nhcat) ! mass power law coefficient for each category
real :: cfmasi   (nhcat) ! inverse of cfmas
real :: pwmas    (nhcat) ! mass power law exponent for each category
real :: pwmasi   (nhcat) ! inverse of pwmas
real :: cfvt     (nhcat) ! fall velocity power law coefficient for each category
real :: pwvt     (nhcat) ! fall velocity power law exponent for each category
real :: dpsmi    (nhcat) ! (1/mass) of 125 micron diameter pristine ice or snow
real :: cfemb0   (nhcat) ! constant dependent on gamma function and mass power law
real :: cfen0    (nhcat) ! constant dependent on gamma function and mass power law
real :: pwemb0   (nhcat) ! constant dependent on gamma function and mass power law
real :: pwen0    (nhcat) ! constant dependent on gamma function and mass power law
real :: vtfac    (nhcat) ! constant dependent on gamma function and mass and fall power laws
real :: frefac1  (nhcat) ! ventilation factor coefficient for constant term
real :: frefac2  (nhcat) ! ventilation factor coefficient for linear term
real :: cfmasft  (nhcat) ! constant dependent on gamma function and mass power law
real :: dnfac    (nhcat) ! constant dependent on gamma function and mass power law 
real :: sipfac   (nhcat) ! constant dependent on gamma function and mass power law
real :: ch3      (nhcat) ! constant dependent on gamma function and mass power law
real :: cdp1     (nhcat) ! constant dependent on mass and fall power laws
real :: pwvtmasi (nhcat) ! constant dependent on mass and fall power laws
real :: reffcof  (nhcat) ! constant dependent on gamma function and mass power law
                         ! for computing effective diameter
real :: dmncof   (nhcat) ! constant dependent on gamma function and mass power law
                         ! for computing mean diameter
real :: dmncofi  (nhcat) ! 1. / dmncof

real :: coltabc(nembc,nembc,npairc) ! collection table for bulk number
real :: coltabx(nembc,nembc,npairx) ! collection table for x bulk density
real :: coltaby(nembc,nembc,npairy) ! collection table for y bulk density

real :: driz_gammq(nembc)

real :: frachz(nrhhz,nthz)       ! Haze nucleation table
real :: fracc(ndnc,ntc,maxgrds0) ! Homogeneous freezing table

real :: gamm(4)          ! complete gamma func of nu, nu+1
real :: gamn1(4)         ! complete gamma func of nu, nu+1

real :: gam(ngam,3)      ! incomplete gamma func
real :: gaminc(ngam,2)   ! incomplete gamma func
real :: gamsip13(2,ngam) ! incomplete gamma func
real :: gamsip24(2,ngam) ! incomplete gamma func

real :: rmlttab(ninc)        ! melting table for bulk density
real :: enmlttab(ninc,nhcat) ! melting table for bulk number
real :: shedtab(ninc,ndns)   ! shedding table

real :: sc(2),sk(2),sl(2) ! specific heats, latent heats
real :: sj(ncat)          ! flag for rain, graupel, hail

real, allocatable :: zfactor_ccn(:)  ! height-based CCN concentration weight factor (when ICCN = 1)
real, allocatable :: zfactor_gccn(:) ! height-based GCCN concentration weight factor (when IGCCN = 1)
real, allocatable :: zfactor_ifn(:)  ! height-based IFN concentration weight factor (when IIFN = 1)

! Sedimentation table section

real :: ch2      (nhcat) ! sedimentation table increment for each category
real :: dispemb0 (nhcat) ! minimum vertical displacement in sedim table for each category
real :: dispemb0i(nhcat) ! inverse of dispemb0
real :: dispemb1 (nhcat) ! maximum vertical displacement in sedim table for each category
real :: vf_max

integer, parameter :: nembfall = 80 ! # of mass values spanning sedimentation table

real, allocatable :: pcpfillc(:,:,:,:) ! sedim table for bulk number
real, allocatable :: pcpfillr(:,:,:,:) ! sedim table for bulk density

integer, allocatable :: kfall(:,:,:)

integer, private :: iefcx, iefcy, ix, iy

real :: ipair(nhcat,nhcat,6)

!         ihx ihy    ipc ipx2 ipy2 ipx ipy ieff

data ipair( 1, 1,:) /   1,  0,  0,  1,  0,  1 /
data ipair( 1, 2,:) /   2,  0,  0,  2,  0,  1 /
data ipair( 2, 2,:) /   3,  0,  0,  0,  0,  1 /
data ipair( 8, 2,:) /   4,  0,  0,  3,  0,  1 /
data ipair( 2, 3,:) /   5,  0,  0,  4,  1, 10 /
data ipair( 3, 3,:) /   6,  0,  0,  5,  0, 10 /
data ipair( 1, 4,:) /   7,  0,  0,  6,  2,  2 /
data ipair( 2, 4,:) /   8,  0,  0,  7,  3, 10 /
data ipair( 3, 4,:) /   9,  0,  0,  8,  4, 10 /
data ipair( 4, 4,:) /  10,  0,  0,  9,  0, 10 /
data ipair( 8, 4,:) /  11,  0,  0, 10,  5,  2 /
data ipair( 9, 4,:) /  12,  0,  0, 11,  6, 10 /
data ipair(10, 4,:) /  13,  0,  0, 12,  7, 10 /
data ipair(11, 4,:) /  14,  0,  0, 13,  8, 10 /
data ipair(12, 4,:) /  15,  0,  0, 14,  9, 10 /
data ipair( 1, 5,:) /  16,  0,  0, 15, 10,  4 /
data ipair( 2, 5,:) /  17,  0,  0, 16, 11, 10 /
data ipair( 3, 5,:) /  18,  0,  0, 17,  0, 10 /
data ipair( 4, 5,:) /  19,  0,  0, 18,  0, 10 /
data ipair( 5, 5,:) /  20,  0,  0,  0,  0, 10 /
data ipair( 8, 5,:) /  21,  0,  0, 19, 12,  4 /
data ipair( 9, 5,:) /  22,  0,  0, 20,  0, 10 /
data ipair(10, 5,:) /  23,  0,  0, 21,  0, 10 /
data ipair(11, 5,:) /  24,  0,  0, 22,  0, 10 /
data ipair(12, 5,:) /  25,  0,  0, 23,  0, 10 /
data ipair(13, 5,:) /  26,  0,  0, 24,  0, 10 /
data ipair(14, 5,:) /  27,  0,  0, 25,  0, 10 /
data ipair(15, 5,:) /  28,  0,  0, 26,  0, 10 /
data ipair(16, 5,:) /  29,  0,  0, 27,  0, 10 /
data ipair( 1, 6,:) /  30,  0,  0, 28, 13,  5 /
data ipair( 2, 6,:) /  31,  0,  0, 29, 14, 10 /
data ipair( 3, 6,:) /  32,  0,  0, 30,  0, 10 /
data ipair( 4, 6,:) /  33,  0,  0, 31,  0, 10 /
data ipair( 5, 6,:) /  34,  0,  0, 32,  0, 10 /
data ipair( 6, 6,:) /  35,  0,  0,  0,  0, 10 /
data ipair( 8, 6,:) /  36,  0,  0, 33, 15,  5 /
data ipair( 9, 6,:) /  37,  0,  0, 34,  0, 10 /
data ipair(10, 6,:) /  38,  0,  0, 35,  0, 10 /
data ipair(11, 6,:) /  39,  0,  0, 36,  0, 10 /
data ipair(12, 6,:) /  40,  0,  0, 37,  0, 10 /
data ipair(13, 6,:) /  41,  0,  0, 38,  0, 10 /
data ipair(14, 6,:) /  42,  0,  0, 39,  0, 10 /
data ipair(15, 6,:) /  43,  0,  0, 40,  0, 10 /
data ipair(16, 6,:) /  44,  0,  0, 41,  0, 10 /
data ipair( 1, 7,:) /  45,  0,  0, 42, 16,  6 /
data ipair( 2, 7,:) /  46,  0,  0, 43, 17, 10 /
data ipair( 3, 7,:) /  47,  0,  0, 44,  0, 10 /
data ipair( 4, 7,:) /  48,  0,  0, 45,  0, 10 /
data ipair( 5, 7,:) /  49,  0,  0, 46,  0, 10 /
data ipair( 6, 7,:) /  50,  0,  0, 47,  0, 10 /
data ipair( 7, 7,:) /  51,  0,  0,  0,  0, 10 /
data ipair( 8, 7,:) /  52,  0,  0, 48, 18,  6 /
data ipair( 9, 7,:) /  53,  0,  0, 49,  0, 10 /
data ipair(10, 7,:) /  54,  0,  0, 50,  0, 10 /
data ipair(11, 7,:) /  55,  0,  0, 51,  0, 10 /
data ipair(12, 7,:) /  56,  0,  0, 52,  0, 10 /
data ipair(13, 7,:) /  57,  0,  0, 53,  0, 10 /
data ipair(14, 7,:) /  58,  0,  0, 54,  0, 10 /
data ipair(15, 7,:) /  59,  0,  0, 55,  0, 10 /
data ipair(16, 7,:) /  60,  0,  0, 56,  0, 10 /
data ipair( 1, 8,:) /  61,104, 55, 57,  0,  1 /
data ipair( 8, 8,:) /  62,  0,  0, 58,  0,  1 /
data ipair( 2, 9,:) /  63,  0,  0, 59, 19, 10 /
data ipair( 9, 9,:) /  64,  0,  0, 60,  0, 10 /
data ipair( 2,10,:) /  65,  0,  0, 61, 20, 10 /
data ipair(10,10,:) /  66,  0,  0, 62,  0, 10 /
data ipair( 2,11,:) /  67,  0,  0, 63, 21, 10 /
data ipair(11,11,:) /  68,  0,  0, 64,  0, 10 /
data ipair( 2,12,:) /  69,  0,  0, 65, 22, 10 /
data ipair(12,12,:) /  70,  0,  0, 66,  0, 10 /
data ipair( 1,13,:) /  71,  0,  0, 67, 23,  4 /
data ipair( 2,13,:) /  72,  0,  0, 68, 24, 10 /
data ipair( 3,13,:) /  73,  0,  0, 69, 25, 10 /
data ipair( 8,13,:) /  74,  0,  0, 70, 26,  4 /
data ipair( 9,13,:) /  75,  0,  0, 71, 27, 10 /
data ipair(10,13,:) /  76,  0,  0, 72, 28, 10 /
data ipair(11,13,:) /  77,  0,  0, 73, 29, 10 /
data ipair(12,13,:) /  78,  0,  0, 74, 30, 10 /
data ipair(13,13,:) /  79,  0,  0, 75,  0, 10 /
data ipair( 1,14,:) /  80,  0,  0, 76, 31,  3 /
data ipair( 2,14,:) /  81,  0,  0, 77, 32, 10 /
data ipair( 3,14,:) /  82,  0,  0, 78, 33, 10 /
data ipair( 8,14,:) /  83,  0,  0, 79, 34,  3 /
data ipair( 9,14,:) /  84,  0,  0, 80, 35, 10 /
data ipair(10,14,:) /  85,  0,  0, 81, 36, 10 /
data ipair(11,14,:) /  86,  0,  0, 82, 37, 10 /
data ipair(12,14,:) /  87,  0,  0, 83, 38, 10 /
data ipair(14,14,:) /  88,  0,  0, 84,  0, 10 /
data ipair( 1,15,:) /  89,  0,  0, 85, 39,  2 /
data ipair( 2,15,:) /  90,  0,  0, 86, 40, 10 /
data ipair( 3,15,:) /  91,  0,  0, 87, 41, 10 /
data ipair( 8,15,:) /  92,  0,  0, 88, 42,  2 /
data ipair( 9,15,:) /  93,  0,  0, 89, 43, 10 /
data ipair(10,15,:) /  94,  0,  0, 90, 44, 10 /
data ipair(11,15,:) /  95,  0,  0, 91, 45, 10 /
data ipair(12,15,:) /  96,  0,  0, 92, 46, 10 /
data ipair(15,15,:) /  97,  0,  0, 93,  0, 10 /
data ipair( 1,16,:) /  98,  0,  0, 94, 47,  3 /
data ipair( 2,16,:) /  99,  0,  0, 95, 48, 10 /
data ipair( 3,16,:) / 100,  0,  0, 96, 49, 10 /
data ipair( 8,16,:) / 101,  0,  0, 97, 50,  3 /
data ipair( 9,16,:) / 102,  0,  0, 98, 51, 10 /
data ipair(10,16,:) / 103,  0,  0, 99, 52, 10 /
data ipair(11,16,:) / 104,  0,  0,100, 53, 10 /
data ipair(12,16,:) / 105,  0,  0,101, 54, 10 /
data ipair(16,16,:) / 106,  0,  0,102,  0, 10 /

! Define parameters for collision efficiency table diameters

integer, parameter :: nefcx(6) = (/21,16,16,16,16,16/)
integer, parameter :: nefcy(6) = (/12, 8, 8, 8,14,14/)

real :: defcx(21,6)
real :: defcy(14,6)
real :: efctab(21,14,6)

! Define 3D collision efficiency table (READ IN REVERSE ORDER)
! (Feingold power laws for RAIN  based on mass from kernelsolver code)

data (defcx(iefcx,1),iefcx=1,nefcx(1)) / & ! 21 values (size ratio)
  .00, .05, .10, .15, .20, .25, .30, .35, .40, .45, &
  .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.00/

data (defcy(iefcy,1),iefcy=1,nefcy(1)) / & ! 12 values [m]
  .00e-3, .02e-3, .04e-3, .06e-3, .08e-3, .10e-3, &
  .12e-3, .14e-3, .20e-3, .30e-3, .40e-3, .60e-3/

data ((efctab(ix,iy,1),iy=1,nefcy(1)),ix=1,nefcx(1)) / &                ! ratio
! 0    20     40     60     80   100    120  140   200  300  400  600 mm
 0.0, .00,   .00,   .00,   .00,  .00,   .00, .00,  .00, .00, .00, .00 & !  .00
,0.0, .0001, .0001, .0001, .001, .005,  .05, .20,  .50, .77, .87, .97 & !  .05
,0.0, .0001, .0001,  .002,  .07,  .40,  .43, .58,  .79, .93, .96, 1.0 & !  .10
,0.0, .0001,  .005,   .02,  .28,  .60,  .64, .75,  .91, .97, .98, 1.0 & !  .15
,0.0, .014,   .016,   .04,  .50,  .70,  .77, .84,  .88, .95, 1.0, 1.0 & !  .20
,0.0, .017,   .022,  .085,  .62,  .78,  .84, .88,  .95, 1.0, 1.0, 1.0 & !  .25
,0.0, .019,    .03,   .17,  .68,  .83,  .87, .90,  1.0, 1.0, 1.0, 1.0 & !  .30
,0.0, .022,   .043,   .27,  .74,  .86,  .89, .92,  1.0, 1.0, 1.0, 1.0 & !  .35
,0.0, .027,   .052,   .40,  .78,  .88,  .90, .94,  1.0, 1.0, 1.0, 1.0 & !  .40
,0.0, .030,   .064,   .50,  .80,  .90,  .91, .95,  1.0, 1.0, 1.0, 1.0 & !  .45
,0.0, .033,   .072,   .55,  .80,  .90,  .91, .95,  1.0, 1.0, 1.0, 1.0 & !  .50
,0.0, .035,   .079,   .58,  .80,  .90,  .91, .95,  1.0, 1.0, 1.0, 1.0 & !  .55
,0.0, .037,   .082,   .59,  .78,  .90,  .91, .95,  1.0, 1.0, 1.0, 1.0 & !  .60
,0.0, .038,   .080,   .58,  .77,  .89,  .91, .95,  1.0, 1.0, 1.0, 1.0 & !  .65
,0.0, .038,   .076,   .54,  .76,  .88,  .92, .95,  1.0, 1.0, 1.0, 1.0 & !  .70
,0.0, .037,   .067,   .51,  .77,  .88,  .93, .97,  1.0, 1.0, 1.0, 1.0 & !  .75
,0.0, .036,   .057,   .49,  .77,  .89,  .95, 1.0,  1.0, 1.0, 1.0, 1.0 & !  .80
,0.0, .035,   .048,   .47,  .78,  .92,  1.0, 1.02, 1.0, 1.0, 1.0, 1.0 & !  .85
,0.0, .032,   .040,   .45,  .79, 1.01, 1.03, 1.04, 1.0, 1.0, 1.0, 1.0 & !  .90
,0.0, .029,   .033,   .47,  .95, 1.30, 1.70, 2.30, 1.0, 1.0, 1.0, 1.0 & !  .95
,0.0, .027,   .027,   .52, 1.40, 2.30, 3.00, 4.00, 1.0, 1.0, 1.0, 1.0 / ! 1.00

! TABLE 2: Cloud or Drizzle + Snow Columns (Wang & Ji 2000)
! (standard RAMS power law relations for COLUMNS)

data (defcx(iefcx,2),iefcx=1,nefcx(2)) / & ! 16 values [m]
    5.e-6, 10.e-6, 15.e-6, 20.e-6, 25.e-6, 30.e-6, 35.e-6, 40.e-6, &
   45.e-6, 50.e-6, 55.e-6, 60.e-6, 65.e-6, 70.e-6, 75.e-6, 80.e-6/

data (defcy(iefcy,2),iefcy=1,nefcy(2)) / & ! 8 values [m]
   .0671e-3, .0933e-3,  .1126e-3,  .1383e-3, & 
   .2374e-3, .5149e-3, 1.0670e-3, 2.4400e-3/

data ((efctab(ix,iy,2),iy=1,nefcy(2)),ix=1,nefcx(2)) / &
! .07   .09   .11   .14   .24   .51 1.07 2.44 mm
  .015, .015, .015, .015, .015, .03, .07, .10, & ! 5 micron
  .03,  .05,  .08,  .12,  .23,  .35, .40, .46, & ! 10
  .18,  .30,  .34,  .39,  .47,  .58, .64, .68, & ! 15
  .29,  .42,  .46,  .52,  .61,  .71, .75, .79, & ! 20
  .29,  .49,  .54,  .61,  .70,  .78, .83, .86, & ! 25
  .24,  .50,  .56,  .62,  .72,  .82, .86, .89, & ! 30
  .00,  .49,  .57,  .63,  .73,  .83, .88, .90, & ! 35
  .00,  .47,  .56,  .64,  .74,  .84, .89, .92, & ! 40
  .00,  .40,  .53,  .62,  .74,  .84, .90, .93, & ! 45
  .00,  .18,  .45,  .59,  .73,  .84, .91, .94, & ! 50
  .00,  .00,  .33,  .55,  .73,  .84, .91, .94, & ! 55
  .00,  .00,  .10,  .48,  .72,  .83, .91, .94, & ! 60
  .00,  .00,  .00,  .33,  .70,  .83, .91, .94, & ! 65
  .00,  .00,  .00,  .33,  .68,  .82, .91, .94, & ! 70
  .00,  .00,  .00,  .00,  .65,  .82, .91, .94, & ! 75
  .00,  .00,  .00,  .00,  .58,  .82, .90, .94  / ! 80

! TABLE 3: Cloud or Drizzle + Snow Broad-Branched Crystals (Wang & Ji 2000)
! (Mitchell (1996) power law relations for BROAD BRANCHED CRYSTALS

data (defcx(iefcx,3),iefcx=1,nefcx(3)) / & ! 16 values [m]
    5.e-6, 10.e-6, 15.e-6, 20.e-6, 25.e-6, 30.e-6, 35.e-6, 40.e-6, &
   45.e-6, 50.e-6, 55.e-6, 60.e-6, 65.e-6, 70.e-6, 75.e-6, 80.e-6/

data (defcy(iefcy,3),iefcy=1,nefcy(3)) / & ! 8 values [m]
    .20e-3,  .25e-3,  .70e-3, 1.00e-3,   &
   1.50e-3, 2.00e-3, 2.50e-3, 3.10e-3/ 

data ((efctab(ix,iy,3),iy=1,nefcy(3)),ix=1,nefcx(3)) / &
! .20   .25   .70  1.00  1.50  2.00 2.50 3.10 mm 
  .00,  .00,  .04,  .05,  .06,  .08, .08, .10, & ! 5 micron
  .00,  .00,  .05,  .08,  .10,  .15, .14, .14, & ! 10
  .00,  .01,  .11,  .20,  .20,  .29, .30, .27, & ! 15
  .00,  .20,  .30,  .39,  .44,  .46, .49, .44, & ! 20
  .00,  .35,  .44,  .52,  .55,  .57, .59, .54, & ! 25
  .00,  .49,  .58,  .65,  .66,  .67, .69, .63, & ! 30
  .00,  .50,  .63,  .69,  .71,  .72, .73, .69, & ! 35
  .00,  .51,  .68,  .74,  .75,  .76, .77, .74, & ! 40
  .00,  .29,  .70,  .77,  .76,  .78, .79, .76, & ! 45
  .00,  .08,  .68,  .79,  .77,  .80, .82, .78, & ! 50
  .00,  .04,  .50,  .79,  .77,  .81, .83, .81, & ! 55
  .00,  .00,  .08,  .80,  .78,  .83, .85, .83, & ! 60
  .00,  .00,  .04,  .79,  .77,  .84, .86, .86, & ! 65
  .00,  .00,  .00,  .79,  .77,  .84, .86, .86, & ! 70
  .00,  .00,  .00,  .73,  .76,  .83, .86, .86, & ! 75
  .00,  .00,  .00,  .68,  .75,  .83, .86, .86  / ! 80

! TABLE 4: Cloud or Drizzle + Ice Hexagonal Plates (Wang & Ji 2000)
! (standard RAMS power law relations for AGGREGATES)

data (defcx(iefcx,4),iefcx=1,nefcx(4)) / & ! 16 values [m]
    5.e-6, 10.e-6, 15.e-6, 20.e-6, 25.e-6, 30.e-6, 35.e-6, 40.e-6, &
   45.e-6, 50.e-6, 55.e-6, 60.e-6, 65.e-6, 70.e-6, 75.e-6, 80.e-6/

data (defcy(iefcy,4),iefcy=1,nefcy(4)) / & ! 8 values [m]
   .1600e-3,  .2266e-3,  .5066e-3,  .7164e-3, & 
   .9476e-3, 1.2400e-3, 1.5000e-3, 1.7000e-3/ 

data ((efctab(ix,iy,4),iy=1,nefcy(4)),ix=1,nefcx(4)) / &
! .16   .23   .51   .72   .95  1.24 1.50 1.70 mm
  .00,  .00,  .03,  .04,  .05,  .06, .07, .08, & ! 5 micron
  .00,  .00,  .08,  .18,  .25,  .29, .30, .31, & ! 10
  .00,  .00,  .40,  .46,  .52,  .53, .54, .55, & ! 15
  .15,  .35,  .57,  .63,  .67,  .68, .69, .70, & ! 20
  .39,  .49,  .70,  .71,  .74,  .77, .78, .79, & ! 25
  .40,  .54,  .75,  .80,  .82,  .83, .84, .84, & ! 30
  .38,  .56,  .78,  .83,  .85,  .86, .87, .87, & ! 35
  .00,  .55,  .82,  .87,  .88,  .89, .90, .90, & ! 40
  .00,  .49,  .85,  .88,  .89,  .90, .91, .91, & ! 45
  .00,  .38,  .85,  .89,  .90,  .91, .92, .92, & ! 50
  .00,  .00,  .86,  .90,  .91,  .92, .93, .93, & ! 55
  .00,  .00,  .85,  .91,  .92,  .93, .94, .94, & ! 60
  .00,  .00,  .83,  .91,  .92,  .93, .94, .94, & ! 65
  .00,  .00,  .82,  .90,  .92,  .93, .94, .94, & ! 70
  .00,  .00,  .80,  .89,  .92,  .93, .94, .94, & ! 75
  .00,  .00,  .75,  .88,  .92,  .93, .94, .94  / ! 80

! TABLE 5: Cloud or Drizzle + Graupel (Cober & List 1993)
! (standard RAMS power law relations for GRAUPEL)

data (defcx(iefcx,5),iefcx=1,nefcx(5)) / & ! 16 values [m]
    5.e-6, 10.e-6, 15.e-6, 20.e-6, 25.e-6, 30.e-6, 35.e-6, 40.e-6, &
   45.e-6, 50.e-6, 55.e-6, 60.e-6, 65.e-6, 70.e-6, 75.e-6, 80.e-6/

data (defcy(iefcy,5),iefcy=1,nefcy(5)) / & ! 14 values [m]
   .08e-3, .16e-3, .24e-3, .32e-3, .40e-3,  .48e-3,  .56e-3, &
   .64e-3, .72e-3, .80e-3, .88e-3, .96e-3, 1.04e-3, 1.12e-3/

data ((efctab(ix,iy,5),iy=1,nefcy(5)),ix=1,nefcx(5)) /                 &
! .08  .16  .24  .32  .40  .48  .56  .64  .72  .80  .88  .96 1.04 1.12 mm
  .34, .26, .21, .18, .15, .13, .11, .09, .08, .07, .06, .05, .04, .03, & ! 5 micron
  .67, .59, .54, .51, .48, .46, .44, .42, .41, .40, .39, .38, .37, .36, & ! 10
  .87, .78, .74, .70, .68, .65, .63, .62, .60, .59, .58, .57, .56, .55, & ! 15
  1.0, .92, .87, .84, .81, .79, .77, .76, .74, .73, .72, .71, .70, .69, & ! 20
  1.0, 1.0, .98, .94, .92, .90, .88, .86, .85, .83, .82, .81, .80, .79, & ! 25
  1.0, 1.0, 1.0, 1.0, 1.0, .98, .96, .95, .93, .92, .91, .90, .89, .88, & ! 30
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .99, .98, .97, .96, .95, & ! 35
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 40
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 45
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 50
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 55
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 60
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 65
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 70
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 75
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  / ! 80

! TABLE 6: Cloud or Drizzle + Hail (Greenan & List 1995)
! (standard RAMS power law relations for HAIL)

data (defcx(iefcx,6),iefcx=1,nefcx(6)) / & ! 16 values [m]
    5.e-6, 10.e-6, 15.e-6, 20.e-6, 25.e-6, 30.e-6, 35.e-6, 40.e-6, &
   45.e-6, 50.e-6, 55.e-6, 60.e-6, 65.e-6, 70.e-6, 75.e-6, 80.e-6/

data (defcy(iefcy,6),iefcy=1,nefcy(6)) / & ! 14 values [m]
    .1e-3,  .2e-3,  .4e-3,  .6e-3,  .9e-3, 1.2e-3,  1.6e-3, &
   2.0e-3, 2.5e-3, 3.0e-3, 4.0e-3, 6.0e-3, 9.0e-3, 12.0e-3/ 
 
data ((efctab(ix,iy,6),iy=1,nefcy(6)),ix=1,nefcx(6)) /                 &
! .1   .2   .4   .6   .9  1.2  1.6  2.0  2.5  3.0  4.0  6.0  9.0  12.0 mm
 .68, .65, .61, .59, .57, .56, .55, .54, .53, .52, .51, .50, .48, .47, & ! 5 micron
 .84, .79, .75, .73, .71, .69, .68, .67, .65, .65, .63, .61, .59, .58, & ! 10
 .95, .90, .85, .82, .80, .78, .77, .75, .74, .73, .71, .69, .67, .66, & ! 15
 1.0, .98, .93, .90, .87, .85, .83, .82, .81, .80, .78, .76, .73, .72, & ! 20
 1.0, 1.0, .99, .96, .93, .91, .89, .88, .86, .85, .83, .81, .78, .77, & ! 25
 1.0, 1.0, 1.0, 1.0, .99, .96, .94, .93, .91, .90, .88, .85, .83, .81, & ! 30
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .99, .97, .96, .94, .92, .89, .87, .85, & ! 35
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .99, .98, .96, .93, .90, .88, & ! 40
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .99, .97, .94, .92, & ! 45
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .97, .95, & ! 50
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, .99, .97, & ! 55
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 60
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 65
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 70
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & ! 75
 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/   ! 80

Contains

!===============================================================================

   subroutine init_nuc_zfactors(mza,zt)

   implicit none

   integer, intent(in) :: mza
   real, intent(in) :: zt(mza)
   
   allocate (zfactor_ccn (mza)) ! height-based CCN concentration weight factor (when ICCN = 1)
   allocate (zfactor_gccn(mza)) ! height-based GCCN concentration weight factor (when IGCCN = 1)
   allocate (zfactor_ifn (mza)) ! height-based IFN concentration weight factor (when IIFN = 1)

   ! Initialize height-based multipliers for concentrations [#/kg_air] of
   ! nucleating aerosols, used when ICCN, IGCCN, and/or IIFN are set to 1  

   ! The following profile decays exponentially with given scale height

   zfactor_ccn (:) = exp(min(0.,-zt(:) / 8000.))
   zfactor_gccn(:) = exp(min(0.,-zt(:) / 8000.)) 
   zfactor_ifn (:) = exp(min(0.,-zt(:) / 8000.)) 

   end subroutine init_nuc_zfactors

End Module micro_coms

