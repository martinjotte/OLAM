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
Module leaf_coms

!-------------------------------------------------------------------------
! This module defines memory for parameters, tables, and other quantities
! that are initialized only once at the beginning of a model run for each
! compute-NODE process.  Afterward, during timesteps, these quantities are
! only read and not written, and may be safely treated as shared memory.
!--------------------------------------------------------------------------

use max_dims,    only: nzgmax, maxndvifiles, pathlen
use consts_coms, only: r8
implicit none

private :: nzgmax, maxndvifiles, pathlen, r8
 
integer, parameter :: nstyp = 12    ! total # of soil textural classes
integer, parameter :: nvtyp = 20    ! total # of leaf ("veg") classes

! Parameters for subroutine vegndvi

real, parameter :: sr_min     = 1.081  ! Minimum simple ratio
real, parameter :: fpar_min   = .001   ! Minimum fpar
real, parameter :: fpar_max   = .950   ! Maximum fpar
real, parameter :: soil_rough = .05    ! soil roughness height
real, parameter :: snow_rough = .01    ! snowcover roughness height

character(pathlen) :: landusefile
character(pathlen) :: veg_database
character(pathlen) :: soil_database
character(pathlen) :: ndvi_database
character(pathlen) :: soilstate_db
character(pathlen) :: soildepth_db

integer            :: nndvifiles, indvifile
character(pathlen) :: fnames_ndvi  (maxndvifiles)
character(14)      :: ctotdate_ndvi(maxndvifiles)
real(r8)           :: s1900_ndvi(maxndvifiles)

integer :: ivegflg
integer :: isoilflg
integer :: ndviflg
integer :: isfcl
integer :: iupdndvi

integer :: isoilstateinit
integer :: isoildepthflg

integer :: nzg      ! # of soil layers
integer :: nzs      ! maximum # of snowcover layers

integer :: nml      ! Total # of land cell M pts in model domain
integer :: nul      ! Total # of land cell U pts in model domain
integer :: nwl      ! Total # of land cell W pts in model domain

integer :: mml      ! Total # of land cell M pts in model parallel sub_domain
integer :: mul      ! Total # of land cell U pts in model parallel sub_domain
integer :: mwl      ! Total # of land cell W pts in model parallel sub-domain

integer :: nslcon   ! soil textural class (optionally used OLAMIN parameter)
integer :: nvgcon   ! leaf class (optionally used OLAMIN parameter)
integer :: mrl_leaf ! mesh refinement level at which leaf is updated

real :: dt_leaf     ! leaf timestep [s]
real :: refdepth    ! reference depth for computing fhydraul [m]
   
real :: slden    (nstyp) ! dry soil density [kg/m3]
real :: slcpd    (nstyp) ! dry soil volumetric heat capacity [J/(m^3 K)]
real :: slbs     (nstyp) ! b exponent [dimensionless]
real :: slcons   (nstyp) ! sat soil hydraulic conductivity [m/s]
real :: slmsts   (nstyp) ! sat volumetric moist content (porosity) [m^3_wat/m^3_tot]
real :: slpots   (nstyp) ! sat moisture potential [m]
real :: soilcp   (nstyp) ! minimum soil moisture [m^3_wat/m^3_tot]
real :: emisg    (nstyp) ! soil infrared emissivity
real :: slcons0  (nstyp) ! surface value for slcons [m/s]
real :: fhydraul (nstyp) ! vertically varying hydraulic conductivity factor
real :: xsand    (nstyp) ! soil fractional sand content
real :: xclay    (nstyp) ! soil fractional clay content
real :: xorgan   (nstyp) ! soil fractional organic content
real :: xrobulk  (nstyp) ! soil bulk density [kg/m^3]
real :: soilcond0(nstyp) ! New soil cond (8/17/00)
real :: soilcond1(nstyp) ! New soil cond (8/17/00)
real :: soilcond2(nstyp) ! New soil cond (8/17/00)

integer :: kroot  (nvtyp)  ! k index of soil layer of lowest roots

real :: albv_green(nvtyp)  ! green vegetation albedo
real :: albv_brown(nvtyp)  ! brown vegetation albedo
real :: emisv     (nvtyp)  ! vegetation infrared emissivity
real :: sr_max    (nvtyp)  ! maximum simple ratio (for using NDVI)
real :: tai_max   (nvtyp)  ! maximum vegetation total area index
real :: sai       (nvtyp)  ! vegetation stem area index
real :: veg_clump (nvtyp)  ! vegetation clumping factor
real :: veg_frac  (nvtyp)  ! vegetation fractional coverage
real :: veg_ht    (nvtyp)  ! vegetation height [m]
real :: glai_max  (nvtyp)  ! vegation maximum green leaf area index
real :: dead_frac (nvtyp)  ! vegetation dead-material fraction
real :: rcmin     (nvtyp)  ! vegetation minimum stomatal resistance [s/m]
real :: dfpardsr  (nvtyp)  ! rate of change of fpar with simple ratio (for using NDVI)
   
real :: slz   (nzgmax)  ! Depth (neg height value) of bottom of each soil layer [m]
real :: slmstr(nzgmax)  ! Initial soil moisture content (range 0. to 1.)

real, save, allocatable :: root   (:,:) ! not used
real, save, allocatable :: slcons1(:,:) ! z-dependent soil sat hydraul cond [m/s]

real, save, allocatable :: dslz     (:) ! soil layer thickness at T pt [m]
real, save, allocatable :: dslzo2   (:) ! HALF soil layer thickness at T pt [m]
real, save, allocatable :: dslzi    (:) ! inverse soil layer thickness at T pt [1/m]
real, save, allocatable :: dslzidt  (:) ! dtll / soil layer thickness at T pt [s/m]
real, save, allocatable :: slzt     (:) ! soil depth at T pt [m]
real, save, allocatable :: dslzt    (:) ! soil layer thickness at M pt [m]
real, save, allocatable :: dslzti   (:) ! inverse soil layer thickness at M pt [1/m]
real, save, allocatable :: dslztidt (:) ! dtll / soil layer thickness at M pt [s/m]

Contains

   subroutine alloc_leafcol()

   implicit none

! Allocate leaf column arrays

   allocate (root   (nzg,nvtyp))
   allocate (slcons1(nzg,nstyp))

   allocate (dslz     (nzg))
   allocate (dslzo2   (nzg))
   allocate (dslzi    (nzg))
   allocate (dslzidt  (nzg))
   allocate (slzt     (nzg))
   allocate (dslzt    (nzg))
   allocate (dslzti   (nzg))
   allocate (dslztidt (nzg))

   return
   end subroutine alloc_leafcol

End Module leaf_coms

