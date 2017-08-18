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

subroutine leaf4_startup()

  use leaf_coms, only: nzg, nzs, landusefile, mml, mul, mwl, &
                       alloc_leafcol, ndviflg, iupdndvi
  use mem_leaf,  only: alloc_leaf, filltab_leaf, land
  use misc_coms, only: io6, runtype

  implicit none

  integer :: iwl

  ! Subroutine LEAF4_STARTUP allocates and initializes some leaf4 arrays.

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

  !-------------------------------------------------------------------------------
  ! STEP 3: Allocate main LEAF arrays
  !-------------------------------------------------------------------------------

  call alloc_leaf(mwl, nzg, nzs)
  call filltab_leaf()

    !-------------------------------------------------------------------------------
  ! STEP 4: Fill ndvi values
  !-------------------------------------------------------------------------------

if ( (iupdndvi /= 1) .and. &
     (runtype == 'HISTORY' .or. runtype == 'HISTADDGRID') ) then

! Do nothing if we are restarting and keeping NDVI constant.
! It will be read in from the history file. However, we still
! need to define values for leaf_init_atm before history read.
     do iwl = 2,mwl
        land%veg_ndvip(iwl) = .5
        land%veg_ndvif(iwl) = .5
        land%veg_ndvic(iwl) = .5
     enddo

elseif (ndviflg == 2) then

  ! Default initialization of NDVI

     do iwl = 2,mwl
        land%veg_ndvip(iwl) = .5
        land%veg_ndvif(iwl) = .5
        land%veg_ndvic(iwl) = .5
     enddo

  elseif (runtype == 'INITIAL' .or. &
          runtype == 'HISTORY' .or. &
          runtype == 'HISTADDGRID') then

  ! Make inventory of ndvi files and initialize current/past time
  ! Not needed for a PLOTONLY run.

     call ndvi_database_read(0)

  ! Initialize future ndvi file time

     call ndvi_database_read(1)

  endif


end subroutine leaf4_startup

!==========================================================================

subroutine sfcdata()

  use leaf_coms, only: slz, nstyp, nvtyp, nzg, slbs_ch, slcpd, slpden,     &
                       slreso2_ch, slreso2_vg, slpots_ch,                  &
                       slmsts_ch, slmsts_vg, slcons_ch, slcons_vg,         &
                       xsand, xclay, xorgan, soilcond0, soilcond1,         &
                       soilcond2, emisg, soilcp_ch, soilcp_vg,             &
                       albv_green, albv_brown,                             &
                       emisv, sr_max, tai_max, sai, veg_clump, veg_frac,   &
                       veg_ht, dead_frac, rcmin, glai_max, dfpardsr,       &
                       fpar_max, fpar_min, sr_min, kroot, soilwilt,        &
                       dt_leaf, dslz, dslzo2, dslzi,                       &
                       dslzidt, slzt, dslzt, dslzti, dslztidt,             &
                       wfrac_high1, wfrac_high2, wfrac_low1, wfrac_low2,   &
                       slpott_high1_ch, slpott_high1_vg,                   &
                       slpott_high2_ch, slpott_high2_vg,                   &
                       slmstsh0_ch, slmstsh0_vg, robulk_ch, robulk_vg,     &
                       headp_high_ch, headp_high_vg,                       &
                       slpott_low1_vg, slpott_low2_vg, headp_low_vg,       &
                       alpha_vg, alphai_vg, en_vg, eni_vg, em_vg, emi_vg,  &
                       slmstsi_ch, slmstsi_vg, slbsi_ch,                   &
                       headp_highi_ch, headp_highi_vg,  headp_lowi_vg,     &
                       slmsts_mscp_vg, slmsts_mscpi_vg


  use misc_coms, only: io6

  implicit none

  integer :: k
  integer :: nnn

  real :: fraclow1, fraclow2
  real :: wfiemi1, wfiemi2
  real :: psilow1, psilow2, psizero
  real :: extrap

!  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke,
!                        1981; Pielke, 1984; Tremback & Kessler, 1985)
!
!  slpots_ch  - sat moisture potential [m]
!  slmsts_ch  - sat volumetric moisture content (soil porosity) [m^3_wat/m^3_tot]
!  slbs_ch    - b exponent [dimensionless]
!  slcons_ch  - saturation soil hydraulic conductivity [m/s]
!  slcons0_ch - surface value for slcons [m/s]
!  slcpd      - dry soil volumetric heat capacity [J/(m^3 K)]
!  slpden     - dry soil density [kg/m3]

  real, parameter :: soilparms1(7,nstyp) = reshape( (/ &
      !------------------------------------------------------------------------
      !slpots_ch    slbs_ch      slcons0_ch        slpden      USDA SOIL CLASS
      !      slmsts_ch    slcons_ch         slcpd               # AND NAME
      !------------------------------------------------------------------------
       -.121, .395,  4.05, .18e-3, .50e-3, 1465.e3, 2650.,  & !  1 sand
       -.090, .410,  4.38, .16e-3, .60e-3, 1407.e3, 2650.,  & !  2 loamy sand
       -.218, .435,  4.9 , .34e-4, .77e-3, 1344.e3, 2650.,  & !  3 sandy loam
       -.786, .485,  5.3 , .72e-5, .11e-4, 1273.e3, 2650.,  & !  4 silt loam
       -.478, .451,  5.39, .69e-5, .22e-2, 1214.e3, 2650.,  & !  5 loam
       -.299, .420,  7.12, .63e-5, .15e-2, 1177.e3, 2650.,  & !  6 sandy clay loam
       -.356, .477,  7.75, .17e-5, .11e-3, 1319.e3, 2650.,  & !  7 silty clay loam
       -.630, .476,  8.52, .24e-5, .22e-2, 1227.e3, 2650.,  & !  8 clay loam
       -.153, .426, 10.4 , .22e-5, .22e-5, 1177.e3, 2650.,  & !  9 sandy clay
       -.490, .492, 10.4 , .10e-5, .10e-5, 1151.e3, 2650.,  & ! 10 silty clay
       -.405, .482, 11.4 , .13e-5, .13e-5, 1088.e3, 2650.,  & ! 11 clay 
       -.356, .863,  7.75, .80e-5, .80e-5,  874.e3,  500./),& ! 12 peat
       (/7,nstyp/) )

  real, parameter :: soilparms2(7,nstyp) = reshape( (/ &
      !-----------------------------------------------------------------
      !xsand      xorgan       soilcond0    soilcond2    USDA SOIL CLASS
      !      xclay      slwilt       soilcond1           # AND NAME
      !-----------------------------------------------------------------
       .97,  .03,  .00,  .070,  .30,  4.80,  -2.70,  & !  1 sand
       .92,  .07,  .01,  .075,  .30,  4.66,  -2.60,  & !  2 loamy sand
       .80,  .18,  .02,  .114,  .29,  4.27,  -2.31,  & !  3 sandy loam
       .57,  .40,  .03,  .179,  .27,  3.47,  -1.74,  & !  4 silt loam
       .60,  .35,  .05,  .155,  .28,  3.63,  -1.85,  & !  5 loam
       .65,  .31,  .04,  .175,  .28,  3.78,  -1.96,  & !  6 sandy clay loam
       .35,  .59,  .06,  .218,  .26,  2.73,  -1.20,  & !  7 silty clay loam
       .48,  .45,  .07,  .250,  .27,  3.23,  -1.56,  & !  8 clay loam
       .50,  .42,  .08,  .219,  .27,  3.32,  -1.63,  & !  9 sandy clay
       .30,  .61,  .09,  .283,  .25,  2.58,  -1.09,  & ! 10 silty clay
       .25,  .65,  .10,  .286,  .25,  2.40,  -0.96,  & ! 11 clay
       .20,  .20,  .60,  .200,  .06,  0.46,   0.00/),& ! 12 peat
       (/7,nstyp/) )

!  Soil Characteristics (see van Genuchten, 1980; Carsel and Parrish, 1988)
!
!  slmsts_vg  - sat volumetric moisture content (soil porosity) [m^3_wat/m^3_tot]
!  soilcp_vg  - minimum soil moisture [m^3_wat/m^3_tot]
!  slcons_vg  - saturation soil hydraulic conductivity [m/s]
!   alpha_vg  - alpha parameter [1/m]
!  alphai_vg  - 1/alpha_vg [m]; ROUGHLY equivalent to slpots
!      en_vg  - n parameter [ ]; equal to (1/slbs + 1) 
!     eni_vg  - 1/n [ ]; equal to (slbs/(1+slbs))
!      em_vg  - m parameter [ ]; could be independent of en_vg, but instead set
!               to (1/(slbs+1))
!     emi_vg  - 1/m [ ]; could be independent of en_vg, but instead set to (slbs+1)

  real, parameter :: soilparms3(5,nstyp) = reshape( (/ &
  !------------------------------------------------------------------------
  ! slmsts_vg  soilcp_vg  slcons_vg  alpha_vg  en_vg       USDA SOIL CLASS
  !                         (m/s)     (1/m)                  # AND NAME
  !------------------------------------------------------------------------
         
       .43,      .045,    .825e-4,    14.5,    2.68,  & !  1 sand
       .41,      .057,    .405e-4,    12.4,    2.28,  & !  2 loamy sand
       .41,      .065,    .123e-4,     7.5,    1.89,  & !  3 sandy loam
       .45,      .067,    .125e-5,     2.0,    1.41,  & !  4 silt loam
       .43,      .078,    .289e-5,     3.6,    1.56,  & !  5 loam
       .39,      .100,    .364e-5,     5.9,    1.48,  & !  6 sandy clay loam
       .43,      .089,    .194e-6,     1.0,    1.23,  & !  7 silty clay loam
       .41,      .095,    .722e-6,     1.9,    1.31,  & !  8 clay loam
       .38,      .100,    .333e-6,     2.7,    1.23,  & !  9 sandy clay
       .36,      .070,    .556e-7,      .5,    1.09,  & ! 10 silty clay
       .38,      .068,    .556e-6,      .8,    1.09,  & ! 11 clay 
       .46,      .034,    .694e-6,     1.6,    1.37/),& ! 12 silt
      (/5,nstyp/) )

! LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

  real, parameter :: bioparms(12,0:nvtyp) = reshape( (/ &
      !-----------------------------------------------------------------------------
      !albv_green     sr_max         veg_clump       rootdep             LEAF-3 CLASS #
      !     albv_brown     tai_max        veg_frac        dead_frac      AND DESCRIPTION
      !          emisv          sai            veg_ht         rcmin
      !-----------------------------------------------------------------------------
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0.,  & !  0  Ocean
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0.,  & !  1  Lakes, rivers, streams
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0.,  & !  2  Ice cap/glacier
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0.,  & !  3  Desert, bare soil
       .14, .24, .97, 5.4, 8.0, 1.0, 1.0, .80, 20.0, 5.0, .0, 500.,  & !  4  Evergreen needleleaf tree
       .14, .24, .95, 5.4, 8.0, 1.0, 1.0, .80, 22.0, 5.0, .0, 500.,  & !  5  Deciduous needleleaf tree
       .20, .24, .95, 6.2, 7.0, 1.0,  .0, .80, 22.0, 5.0, .0, 500.,  & !  6  Deciduous broadleaf tree
       .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 5.0, .0, 500.,  & !  7  Evergreen broadleaf tree
       .21, .43, .96, 5.1, 4.0, 1.0,  .0, .75,   .3,  .7, .7, 100.,  & !  8  Short grass
       .24, .43, .96, 5.1, 5.0, 1.0,  .0, .80,  1.2, 1.0, .7, 100.,  & !  9  Tall grass
       .24, .24, .96, 5.1, 1.0,  .2, 1.0, .20,   .7, 1.0, .0, 500.,  & ! 10  Semi-desert
       .20, .24, .95, 5.1, 4.5,  .5, 1.0, .60,   .2, 1.0, .0,  50.,  & ! 11  Tundra
       .14, .24, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500.,  & ! 12  Evergreen shrub
       .20, .28, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500.,  & ! 13  Deciduous shrub
       .16, .24, .96, 6.2, 7.0, 1.0,  .5, .80, 22.0, 5.0, .0, 500.,  & ! 14  Mixed woodland
       .22, .40, .95, 5.1, 5.0,  .5,  .0, .85,  1.0, 1.0, .0, 100.,  & ! 15  Crop/mixed farming, C3 grassland
       .18, .40, .95, 5.1, 5.0,  .5,  .0, .80,  1.1, 1.0, .0, 500.,  & ! 16  Irrigated crop
       .12, .43, .98, 5.1, 7.0, 1.0,  .0, .80,  1.6, 1.0, .0, 500.,  & ! 17  Bog or marsh
       .20, .36, .96, 5.1, 6.0, 1.0,  .0, .80,  7.0, 3.5, .0, 100.,  & ! 18  Wooded grassland 
       .20, .36, .90, 5.1, 3.6, 1.0,  .0, .74,  6.0, 2.5, .0, 500.,  & ! 19  Urban and built up
       .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 5.0, .0, 500.,  & ! 20  Wetland evergreen broadleaf tree
       .18, .18, .96, 5.1, 3.0,  .5,  .0, .80, 0.32, 1.0, .0, 100./),& ! 21  Deforested (Amazon - Medvigy)
      !.16, .24, .96, 5.1, 2.0, 1.5, 1.0, .10, 20.0, 1.5, .0, 500./),& ! 21  Very urban
       (/12,nvtyp+1/) )

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

  do nnn = 1,nstyp

     ! Soil parameters independent of soil model

     slcpd    (nnn) = soilparms1(6,nnn)
     slpden   (nnn) = soilparms1(7,nnn)
     xsand    (nnn) = soilparms2(1,nnn)
     xclay    (nnn) = soilparms2(2,nnn)
     xorgan   (nnn) = soilparms2(3,nnn)
     soilwilt (nnn) = soilparms2(4,nnn)
     soilcond0(nnn) = soilparms2(5,nnn)
     soilcond1(nnn) = soilparms2(6,nnn)
     soilcond2(nnn) = soilparms2(7,nnn)
     emisg    (nnn) = .98

     ! Clapp & Hornberger soil model parameters

     slpots_ch      (nnn) = soilparms1(1,nnn)
     slmsts_ch      (nnn) = soilparms1(2,nnn)
     slbs_ch        (nnn) = soilparms1(3,nnn)
     slcons_ch      (nnn) = soilparms1(4,nnn)

     robulk_ch      (nnn) = slpden(nnn) * (1.0 - slmsts_ch(nnn))
     soilcp_ch      (nnn) = 0.1 - 0.07 * xsand(nnn)
     slbsi_ch       (nnn) = 1.0 / slbs_ch(nnn)
     slmstsi_ch     (nnn) = 1.0 / slmsts_ch(nnn)
     slpott_high1_ch(nnn) = slpots_ch(nnn) * (1./wfrac_high1) ** slbs_ch(nnn)
     slpott_high2_ch(nnn) = 10.

     headp_high_ch  (nnn) = (slpott_high2_ch(nnn) - slpott_high1_ch(nnn)) &
                          / ((wfrac_high2 - wfrac_high1) * slmsts_ch(nnn))
     headp_highi_ch (nnn) = 1.0 / headp_high_ch(nnn)

     slmstsh0_ch    (nnn) = slmsts_ch(nnn) * (wfrac_high1 &
                          + (wfrac_high2 - wfrac_high1) * (-slpott_high1_ch(nnn) - slz(nzg)) &
                          / (slpott_high2_ch(nnn) - slpott_high1_ch(nnn)))

     ! van Genuchten soil model parameters

           slmsts_vg(nnn) =  soilparms3(1,nnn)
           soilcp_vg(nnn) =  soilparms3(2,nnn)
           slcons_vg(nnn) =  soilparms3(3,nnn)
            alpha_vg(nnn) = -soilparms3(4,nnn)
               en_vg(nnn) =  soilparms3(5,nnn)
          slmstsi_vg(nnn) = 1. / slmsts_vg(nnn)
           alphai_vg(nnn) = 1. / alpha_vg(nnn)
              eni_vg(nnn) = 1. / en_vg(nnn)
               em_vg(nnn) = 1. - eni_vg(nnn)
              emi_vg(nnn) = 1. / em_vg(nnn)
           robulk_vg(nnn) = slpden(nnn) * (1.0 - slmsts_vg(nnn))

     slpott_high1_vg(nnn) = alphai_vg(nnn) &
                          * (wfrac_high1**(-emi_vg(nnn)) - 1.)**eni_vg(nnn)
     slpott_high2_vg(nnn) = 10.

       headp_high_vg(nnn) = (slpott_high2_vg(nnn) - slpott_high1_vg(nnn)) &
                          / ((wfrac_high2 - wfrac_high1) * (slmsts_vg(nnn) - soilcp_vg(nnn)))
      headp_highi_vg(nnn) = 1.0 / headp_high_vg(nnn)

         slmstsh0_vg(nnn) = slmsts_vg(nnn) * (wfrac_high1 &
                          + (wfrac_high2 - wfrac_high1) * (-slpott_high1_vg(nnn) - slz(nzg)) &
                          / (slpott_high2_vg(nnn) - slpott_high1_vg(nnn)))

      slmsts_mscp_vg(nnn) = slmsts_vg(nnn) - soilcp_vg(nnn)
     slmsts_mscpi_vg(nnn) = 1.0 / slmsts_mscp_vg(nnn)

  ! Construct constant head gradient at lower limit of van Genuchten soil moisture
  ! curve to avoid singularity.

                   wfiemi1 = wfrac_low1**(-emi_vg(nnn))
                   wfiemi2 = wfrac_low2**(-emi_vg(nnn))
       slpott_low1_vg(nnn) = alphai_vg(nnn) * (wfiemi1 - 1.)**eni_vg(nnn)
       slpott_low2_vg(nnn) = alphai_vg(nnn) * (wfiemi2 - 1.)**eni_vg(nnn)

         headp_low_vg(nnn) = (slpott_low2_vg(nnn) - slpott_low1_vg(nnn)) &
                           / ((wfrac_low2 - wfrac_low1) * (slmsts_vg(nnn) - soilcp_vg(nnn)))

  ! Increase headp_low_vg as necessary so that head is at or below -1.e5 m
  ! at wfrac = 0. 

         headp_low_vg(nnn) = max(headp_low_vg(nnn), &
                             (slpott_low2_vg(nnn) + 1.e5) / wfrac_low2)

        headp_lowi_vg(nnn) = 1.0 / headp_low_vg(nnn)

     do k = 1,nzg
        slreso2_ch(k,nnn) = dslzo2(k) / slcons_ch(nnn)
        slreso2_vg(k,nnn) = dslzo2(k) / slcons_vg(nnn)
     enddo
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
  enddo

end subroutine sfcdata
