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
  integer, parameter :: nvtyp = 21    ! total # of leaf ("veg") classes

  ! Parameters for subroutine vegndvi

  real, parameter :: sr_min     = 1.081  ! Minimum simple ratio
  real, parameter :: fpar_min   =  .001  ! Minimum fpar
  real, parameter :: fpar_max   =  .950  ! Maximum fpar
  real, parameter :: soil_rough =  .050  ! soil roughness height
  real, parameter :: snow_rough =  .010  ! snowcover roughness height

  ! Parameters for soil...

  ! At the higher end of the range of soil water content in both the van
  ! Genuchten and Clapp & Hornberger models for soil water potential (head),
  ! water potential is modified to represent a linear increase with water
  ! content so that a positive confinement pressure can be represented.  The
  ! constant gradient is defined using two values of water fraction,
  ! wfrac_high1 and wfrac_high2, along with corresponding head values for each.

  real, parameter :: wfrac_high1 =  0.95 ! Soil water fraction values at high
  real, parameter :: wfrac_high2 =  1.05 ! and low ends of normal range used to

  ! At the lower end of the range of soil water content in the van Genuchten
  ! model for soil water potential (head), water potential becomes singular.
  ! To avoid computational problems associated with the singularity and nearby
  ! steepness of the potential curve, water potential is modified to represent
  ! a linear increase with water.  The constant gradient is defined using two
  ! values of water fraction, wfrac_low1 and wfrac_low2, along with
  ! corresponding head values for each.

  real, parameter :: wfrac_low1  =  0.009 ! define zones of constant potential
  real, parameter :: wfrac_low2  =  0.010 ! (or head) gradient with respect to
                                         ! soil water content
  
  character(pathlen) :: landusefile
  character(pathlen) :: veg_database
  character(pathlen) :: soil_database
  character(pathlen) :: ndvi_database
  character(pathlen) :: soilstate_db
  character(pathlen) :: soildepth_db
  character(pathlen) :: watertab_db

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
  integer :: iwatertabflg

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
   
  real :: slpden      (nstyp) ! soil particle density [kg/m3]
  real :: slcpd       (nstyp) ! dry soil volumetric heat capacity [J/(m^3 K)]
  real :: slbs        (nstyp) ! b exponent [dimensionless]
  real :: slcons      (nstyp) ! sat soil hydraulic conductivity [m/s]
  real :: slmsts      (nstyp) ! sat volumetric moist content (porosity) [m^3_wat/m^3_tot]
  real :: slmstsh0    (nstyp) ! nzg volumetric moist content at head = 0 [m^3_wat/m^3_tot]
  real :: slpots      (nstyp) ! sat moisture potential [m]
  real :: soilcp      (nstyp) ! minimum soil moisture [m^3_wat/m^3_tot]
  real :: emisg       (nstyp) ! soil infrared emissivity
  real :: slcons0     (nstyp) ! surface value for slcons [m/s]
  real :: xsand       (nstyp) ! soil fractional sand content
  real :: xclay       (nstyp) ! soil fractional clay content
  real :: xorgan      (nstyp) ! soil fractional organic content
  real :: robulk      (nstyp) ! soil dry bulk density [kg/m^3]
  real :: soilwilt    (nstyp) ! vol moist content at wilting point [m^3_wat/m^3_tot]
  real :: soilcond0   (nstyp) ! New soil cond (8/17/00)
  real :: soilcond1   (nstyp) ! New soil cond (8/17/00)
  real :: soilcond2   (nstyp) ! New soil cond (8/17/00)
  real :: slmstsi     (nstyp) ! 1 / porosity
  real :: headp_high  (nstyp) ! Gradient of soil water potential at high water content
  real :: slpott_high1(nstyp) ! soil water potential [m] at wfrac_high1
  real :: slpott_high2(nstyp) ! soil water potential [m] at wfrac_high2

! Parameters for van Genuchten model (an alternative to Clapp & Hornberger) 
  real ::       slmsts_vg(nstyp) ! sat volumetric moist content (porosity) [m^3_wat/m^3_tot]
  real ::     slmstsh0_vg(nstyp) ! nzg volumetric moist content at head = 0 [m^3_wat/m^3_tot]
  real ::       soilcp_vg(nstyp) ! minimum soil moisture [m^3_wat/m^3_tot]
  real ::       slcons_vg(nstyp) ! sat soil hydraulic conductivity [m/s]
  real ::        alpha_vg(nstyp) ! alpha parameter [1/m]
  real ::       alphai_vg(nstyp) ! 1/alpha parameter [m]; ROUGHLY equivalent to slpots
  real ::           en_vg(nstyp) ! n parameter [ ]
  real ::          eni_vg(nstyp) ! 1/n parameter [ ]
  real ::           em_vg(nstyp) ! m parameter [ ]
  real ::          emi_vg(nstyp) ! 1/m parameter [ ]
  real ::   headp_high_vg(nstyp) ! Gradient of soil water potential at high water content
  real ::    headp_low_vg(nstyp) ! Gradient of soil water potential at low water content
  real :: slpott_high1_vg(nstyp) ! soil water potential [m] at wfrac_high1
  real :: slpott_high2_vg(nstyp) ! soil water potential [m] at wfrac_high2
  real ::  slpott_low1_vg(nstyp) ! soil water potential [m] at wfrac_low1
  real ::  slpott_low2_vg(nstyp) ! soil water potential [m] at wfrac_low2

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

  real, save, allocatable :: slreso2   (:,:) ! z-dependent soil sat hydraul resist [s]
  real, save, allocatable :: slreso2_vg(:,:) ! z-dependent soil sat hydraul resist [s]
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

    allocate (slreso2   (nzg,nstyp))
    allocate (slreso2_vg(nzg,nstyp))
    allocate (dslz    (nzg))
    allocate (dslzo2  (nzg))
    allocate (dslzi   (nzg))
    allocate (dslzidt (nzg))
    allocate (slzt    (nzg))
    allocate (dslzt   (nzg))
    allocate (dslzti  (nzg))
    allocate (dslztidt(nzg))

  end subroutine alloc_leafcol

End Module leaf_coms

