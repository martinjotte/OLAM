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
Module nuclei_coms

  !-----------------------------------------------------
  ! PHYSICAL PARAMETERS USED IN DUST SOURCE CALCULATIONS
  !-----------------------------------------------------

  !>>> dust 
  !Set up the bins/apportionment (Ginoux et al. 2001)
  !r1=0.15, r2=0.265, r3=0.471, r4=0.838, r5=1.5, r6=2.65, r7=4.71 microns
  !For the 4 small clay size bins, each bin contains the following fractions:
  ! 1=0.9%, 2=08.1%, 3=23.4%, 4=67.6%
  !However, the 4 small clay bins represent only 1/10th the total dust silt,
  ! so the parameter "sp(1:4)" represents relative abundance
  !The three large classes are weight equally at 30%.

  integer, parameter :: nbin = 7
  real, parameter :: rhoa = 1.25e-3, grav = 981., fact = 1.153e-16, pi = 3.141593

  real :: rad(nbin) ! dust particle radius [cm] for each bin
  real ::  sp(nbin) ! relative abundance [fractional #] of each bin
  real :: den(nbin) ! density [g/cm^3] of particles in each bin

  data rad /0.15e-4, 0.265e-4, 0.471e-4, 0.838e-4, 1.5e-4, 2.65e-4, 4.71e-4/
  data sp  / 0.0009,   0.0081,   0.0234,   0.0676,    0.3,     0.3,     0.3/
  data den /    2.5,      2.5,      2.5,      2.5,   2.65,    2.65,    2.65/

  logical :: nodust(0:21) ! True if no dust can be generated from a landuse class

  real :: amassi(nbin) ! Inverse mass of single particle for each bin
  real :: uth(nbin)    ! Threshold velocity [cm/s] from Marticorena & Bermagetti 1995

  ! Soil textural classes are: 1 sand, 2 loamy sand, 3 sandy loam, 4 silt loam,
  ! 5 loam, 6 sandy clay loam, 7 silty clay loam, 8 clay loam, 9 sandy clay, 
  ! 10 silty clay, 11 clay, 12 peat, 13 silt.  The assumed clay content of each
  ! is given in clayfrac.

  real :: clayfrac(13)
  data clayfrac /  0.0,  5.0, 10.0, 12.0, 18.0, & !  1 -  5
                  28.0, 33.0, 33.0, 42.0, 48.0, & !  6 - 10
                  70.0,  0.0,  0.0 /              ! 11 - 13

  ! wprime is used in parameterization by Fecan et al. and is based on clay
  ! percentage using the formula wprime = 0.0014(clayfrac)^2 + 0.17(clayfrac).

  real :: wprime(13)

  !---------------------------------------------------------------------
  ! PHYSICAL PARAMETERS USED IN NUCLEI DRY & WET DEPOSITION CALCULATIONS
  !---------------------------------------------------------------------

  !******************************************************************************
  !* salt_dust_sources.f90 computes  sea salt (2 modes), dust (2 modes)         *
  !* HERE SMALL MODE OF SEASALT -----> CCN                                      *
  !* HERE LARGE MODE OF SEASALT -----> GCCN                                     *
  !*                                                                            *
  !* For radiation purposes the small mode is a GCCN (large mode here)          *
  !* For radiation purposes the large mode is supergiant CCN (non-existent here)*
  !                                                                             *
  !* SO, the small "radiative" mode is updated as the large "source" mode       *
  !* AND the large "radiative" mode is NOT updated here                         *
  !*									        *
  !* Added: Seigel (9-23-10)						        *
  !******************************************************************************

  ! Constants from RAMS subroutine calvar

  real, parameter ::        pi1 = 3.1415926  ! pi
  real, parameter ::        pi3 = 3. * pi1   ! 3 * pi
  real, parameter ::      boltz = 1.3807e-23 ! Boltzmann constant
  real, parameter ::     t_suth = 110.4      ! Sutherland temperature [K]
  real, parameter ::      t_ref = 293.15     ! Reference temperature [K]
  real, parameter ::      p_ref = 101325.    ! Reference Pressure [Pa]
  real, parameter ::    dva_ref = 1.8203e-5  ! Dynamic viscosity of air at t_ref,p_ref [kg/(m s)]
  real, parameter :: airmfp_ref = 0.0651e-6  ! Air mean free path at t_ref,p_ref [m]
                                             ! Seinfeld and Pandis (2006) eq.(9.7)
  real, parameter :: e0 = 3.

  integer, parameter :: nddc = 15 ! number of dry-deposition landuse categories
  integer, parameter :: nt = 6    ! viscosity of water at different temperature 

  real :: z0    (nddc) ! roughness length
  real :: acoll (nddc) ! radius of collectors
  real :: alpha (nddc) ! property of vegtype
  real :: ggamma(nddc) ! property of vegtype  

  integer, parameter :: nvtyp = 21 ! should correspond to nvtyp in leaf_coms

  integer :: jsfcinx(0:nvtyp)

  ! Map 21 leaf classes into 15 "use types" for which data values 
  ! (dimensioned to nddc) are defined below

  data jsfcinx(0:nvtyp) /0,0,12,8,1,3,4,2,6,6,8,9,10,10,10,7,7,11,6,15,1,6/

  data acoll(1:nddc) /2.0,5.0,2.0, 5.0,5.0,2.0,2.0,0.0,0.0, &
                      10.0,10.0,0.0,0.0,0.0,10.0/

  data alpha(1:nddc) /1.0,  0.6,  1.1,  0.8,   0.8,   1.2, 1.2, 50.0, &
                     50.0,  1.3,  2.0, 50.0, 100.0, 100.0, 1.5/

  data ggamma(1:nddc) /0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54, &
                       0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.56/ 

  ! Dynamic viscosity of liquid water [kg/(m s)] as a function of temperature [K]

  ! /       0.,       5.,      10.,      20.,      30.,      40. /
  ! / 1.787e-3, 1.519e-3, 1.307e-3, 1.002e-3, 0.798e-3, 0.653e-3 /

  Contains

!===============================================================================

  subroutine dust_src_init()

     implicit none

     integer :: ibin, nts
     real :: diam

     ! Set logical flags that indicate which landuse classes produce no dust 

     nodust( 0) = .true.  ! Ocean
     nodust( 1) = .true.  ! Lakes, rivers, streams
     nodust( 2) = .true.  ! Ice cap/glacier
     nodust( 3) = .false. ! Desert, bare soil
     nodust( 4) = .true.  ! Evergreen needleleaf tree
     nodust( 5) = .true.  ! Deciduous needleleaf tree
     nodust( 6) = .true.  ! Deciduous broadleaf tree
     nodust( 7) = .true.  ! Evergreen broadleaf tree
     nodust( 8) = .false. ! Short grass
     nodust( 9) = .false. ! Tall grass
     nodust(10) = .false. ! Semi-desert
     nodust(11) = .true.  ! Tundra
     nodust(12) = .true.  ! Evergreen shrub
     nodust(13) = .true.  ! Deciduous shrub
     nodust(14) = .true.  ! Mixed woodland
     nodust(15) = .false. ! Crop/mixed farming, C3 grassland
     nodust(16) = .true.  ! Irrigated crop
     nodust(17) = .true.  ! Bog or marsh
     nodust(18) = .false. ! Wooded grassland 
     nodust(19) = .true.  ! Urban and built up
     nodust(20) = .true.  ! Wetland evergreen broadleaf tree
     nodust(21) = .true.  ! Deforested (Amazon - Medvigy)

     ! Loop over dust bins and for each compute particle inverse mass and
     ! threshold wind speed

     do ibin = 1,nbin
        amassi(ibin) = 1. / (den(ibin) * (4./3.) * pi * rad(ibin)**3)

        if (ibin <= 4) then  ! For first 4 bins (sub-micron), assume rad = 0.75 micron
           diam = 1.5e-4
        else
           diam = 2. * rad(ibin)
        endif

        uth(ibin) = 0.13 * sqrt(den(ibin) * grav * diam / rhoa) &
                  * sqrt(1. + 0.006 / (den(ibin) * grav * diam**2.5)) &
                  / sqrt(1.928 * (1331. * diam**1.56 + 0.38)**0.092 - 1.)
     enddo

     ! Compute wprime value for each soil textural class

     do nts = 1,13
        wprime(nts) = 0.0014 * clayfrac(nts)**2 + 0.17 * clayfrac(nts)
     enddo

  end subroutine dust_src_init

end module nuclei_coms

