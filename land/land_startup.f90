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

subroutine land_startup()

  use leaf_coms, only: nzs, ndviflg, iupdndvi, isoilflg
  use mem_land,  only: alloc_land, filltab_land, land, mland, nzg, slzt, omland
  use misc_coms, only: io6, runtype, isubdomain
  use mem_sfcg,  only: itab_wsfc
  use mem_para,  only: myrank

  implicit none

  real, parameter :: specifheat_bedrock = 2.2e6 ! Specific heat of bedrock [J/(m^3 K)]
  integer :: iland, k, iwsfc

  ! Initialize some land & vegetation properties

  call land_parms()

  ! Allocate time-dependent LEAF arrays and add to history file I/O table

  call alloc_land(mland, nzg, nzs)
  call filltab_land()

  ! Fill ndvi values

  if ( (iupdndvi /= 1) .and. &
       (runtype == 'HISTORY' .or. runtype == 'HISTADDGRID') ) then

     ! Do nothing if we are restarting and keeping NDVI constant.
     ! It will be read in from the history file. However, we still
     ! need to define values for leaf_init_atm before history read.

     do iland = 2,mland
        land%veg_ndvip(iland) = .5
        land%veg_ndvif(iland) = .5
        land%veg_ndvic(iland) = .5
     enddo

  elseif (ndviflg == 2) then

     ! Default initialization of NDVI

     do iland = 2,mland
        land%veg_ndvip(iland) = .5
        land%veg_ndvif(iland) = .5
        land%veg_ndvic(iland) = .5
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

  ! Apply pedotransfer functions to obtain soil hydraulic properties and
  ! specific heat.  Loop over all land points, INCLUDING THOSE THAT ARE NOT
  ! PRIMARY ON THIS SUBDOMAIN. 

  do iland = 2,mland
     iwsfc = iland + omland

     do k = 1,nzg

        ! If isoilflg = 1, meaning that SoilGrids and GLHYMPS datasets are used,
        ! define a level in the soil above which SoilGrids composition and PTFs
        ! are used and below which GLHYMPS permeability and porosity are used.

        ! z_bedrock, the height of the upper boundary of bedrock, is used for
        ! this transition grid level when bedrock is near the surface.  However,
        ! where bedrock begins much deeper, the transition level is set above
        ! z_bedrock because SoilGrids data is defined only in the top 2 m, while
        ! GLHYMPS applies to roughly the top 100 m.  For now, we choose the 
        ! transition level to be no greater than 10 meters below the surface. 

        if (isoilflg == 1 .and. slzt(k) < max(-10.0, land%z_bedrock(iland))) then

           ! ISOILFLG = 1 and this soil grid level is below transition level;
           ! therefore, assign wsat and ksat based on glhymps dataset.

           ! We choose to set wresid_vg to 0.2 * wsat_vg, which is an average
           ! ratio for soils.  We also choose to set alpha_vg and en_vg to
           ! values given for USDA soil textural class of silt loam, which is
           ! an average soil.  lambda_vg is set to 0.5, which is the value used
           ! in some, but not all, pedotransfer functions.  There is no
           ! particular justification for these choices other than the need to
           ! provide something for the water retention and hydraulic conductivity
           ! curves.  Nevertheless, we expect that fairly moist conditions will
           ! prevail at these depths, in which case alpha_vg and en_vg are of
           ! minor importance.

           ! pH_soil is not used for physical computations in bedrock layers,
           ! leaving the array free for other uses.  We set pH_soil to -10.0
           ! in bedrock to be used as an identifier that a soil grid level is
           ! filled with bedrock.  

           land%wresid_vg         (k,iland) = 0.2 * land%glhymps_poros(iland)
           land%wsat_vg           (k,iland) =       land%glhymps_poros(iland)
           land%ksat_vg           (k,iland) =       land%glhymps_ksat (iland)
           land%alpha_vg          (k,iland) = -2.0
           land%en_vg             (k,iland) = 1.4
           land%lambda_vg         (k,iland) = 0.5
           land%pH_soil           (k,iland) = -10.0
           land%specifheat_drysoil(k,iland) = specifheat_bedrock

        else

           ! Either ISOILFLG /= 1 or this soil grid level is above transition level;
           ! therefore, apply soil pedotransfer functions

           call soil_ptf(iland,k,                          &
                         slzt                   (k),       &
                         land%usdatext            (iland), &
                         land%sand              (k,iland), &
                         land%clay              (k,iland), &
                         land%silt              (k,iland), &
                         land%organ             (k,iland), &
                         land%bulkdens_drysoil  (k,iland), &
                         land%cec_soil          (k,iland), &
                         land%pH_soil           (k,iland), &
                         land%wresid_vg         (k,iland), &
                         land%wsat_vg           (k,iland), &
                         land%alpha_vg          (k,iland), &
                         land%en_vg             (k,iland), &
                         land%ksat_vg           (k,iland), &
                         land%lambda_vg         (k,iland), &
                         land%specifheat_drysoil(k,iland)  )

        endif

     enddo
  enddo

end subroutine land_startup

!==========================================================================

subroutine land_parms()

  use leaf_coms,  only: nvtyp, albv_green, albv_brown, emisv,         &
                        sr_max, tai_max, sai, veg_clump, veg_frac,    &
                        veg_ht, dead_frac, rcmin, glai_max, dfpardsr, &
                        fpar_max, fpar_min, sr_min, z_root, kroot,    &
                        snowmin_expl, wcap_min, wcap_vmin, dt_leaf

  use mem_land,   only: nzg, slz, kperc, hptimi

  use misc_coms,  only: io6
  use oname_coms, only: nl

  implicit none

  integer :: k
  integer :: iveg

  real :: zdiff

  ! LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

  real, parameter :: bioparms(12,0:nvtyp) = reshape( (/ &
      !-----------------------------------------------------------------------------
      !albv_green     sr_max         veg_clump       rootdepth           LEAF-3 CLASS #
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

  do iveg = 1,nvtyp
     albv_green(iveg) =  bioparms( 1,iveg)
     albv_brown(iveg) =  bioparms( 2,iveg)
     emisv     (iveg) =  bioparms( 3,iveg)
     sr_max    (iveg) =  bioparms( 4,iveg)
     tai_max   (iveg) =  bioparms( 5,iveg)
     sai       (iveg) =  bioparms( 6,iveg)
     veg_clump (iveg) =  bioparms( 7,iveg)
     veg_frac  (iveg) =  bioparms( 8,iveg)
     veg_ht    (iveg) =  bioparms( 9,iveg)
     z_root    (iveg) = -bioparms(10,iveg)
     dead_frac (iveg) =  bioparms(11,iveg)
     rcmin     (iveg) =  bioparms(12,iveg)
     glai_max  (iveg) = tai_max(iveg) - sai(iveg)
     dfpardsr  (iveg) = (fpar_max - fpar_min) / (sr_max(iveg) - sr_min)
     kroot     (iveg) = nzg

     do k = nzg-1,1,-1
        if (slz(k+1) > -bioparms(10,iveg)) kroot(iveg) = k
     enddo
  enddo

  ! Minimum snow/water mass

  wcap_min  = dt_leaf * 1.e-7   !  1.e-7 is 1 mm pcp/dew in 115 days
  wcap_vmin = dt_leaf * 1.e-8

  ! Check whether running surface model as stand-alone with long timestep

  if (nl%igw_spinup /= 1) then

     ! Standard run with ATM coupling and short timestep

     snowmin_expl = max(10.0, 0.04 * dt_leaf)
     hptimi = 1.0 / 1800.

     ! Choose percolation level to be slz level that is closest to -3.0 m in height

     zdiff = 1.e4

     do k = 1,nzg
        if (zdiff > abs(slz(k) + 3.)) then
            zdiff = abs(slz(k) + 3.)
            kperc = k
        endif
     enddo

     write(io6,'(a,i5,f8.2)') 'Percolation level at kperc,slz(kperc) ',kperc,slz(kperc)

  else

     ! Surface stand-alone run with long timestep

     snowmin_expl = 300.
     hptimi = 1.0 / (86400. * 5.)

  endif

end subroutine land_parms

!==========================================================================

subroutine soil_ptf(iland, k, slzt, usdatext, sand, clay, silt,   &
                    organ0, bulkdens_drysoil, cec_soil, pH_soil,  &
                    wresid_vg, wsat_vg, alpha_vg, en_vg, ksat_vg, &
                    lambda_vg, specifheat_drysoil                 )

  use leaf_coms, only: isoilptf

  ! Computes static (time-independent) thermal and hydraulic parameters for
  ! soil layers, but not for bedrock layers.

  implicit none

  integer, intent(in)  :: iland, k, usdatext
  real,    intent(in)  :: slzt, sand, clay, silt, organ0
  real,    intent(in)  :: bulkdens_drysoil, cec_soil, pH_soil
  real,    intent(out) :: wresid_vg, wsat_vg, alpha_vg, en_vg, ksat_vg, &
                          lambda_vg, specifheat_drysoil

  ! Soil composition and van Genuchten hydraulic parameters based on USDA
  ! soil textural class (van Genuchten, 1980; Carsel and Parrish, 1988)

  ! Soil compositoin parameters are not applied here; they are included in the following
  ! parameter statement only to keep it identical with subroutine usda_composition.

  !   wsat_vg - sat volumetric moisture content (soil porosity) [m^3_wat/m^3_tot]
  ! wresid_vg - minimum soil moisture [m^3_wat/m^3_tot]
  !   ksat_vg - saturation soil hydraulic conductivity [m/s]
  !  alpha_vg - alpha parameter [1/m]; Alpha value in model is NEGATIVE; a coef of SWP, not tension
  !     en_vg - n parameter [ ]
  ! lambda_vg - lambda parameter [ ]

  real, parameter :: soilparms4(10,12) = reshape( (/ &
  !------------------------------------------------------------------------------------------------------
  ! sand   clay   silt  organ  wresid_vg wsat_vg  ksat_vg  alpha_vg  en_vg  slwilt      USDA SOIL CLASS
  !                                                 (m/s)    (1/m)                       # AND NAME
  !------------------------------------------------------------------------------------------------------
    .92,   .03,   .05,   .00,   .045,    .43,    .825e-4,    14.5,   2.68,  .070,   & !  1 sand
    .82,   .06,   .12,   .01,   .057,    .41,    .405e-4,    12.4,   2.28,  .075,   & !  2 loamy sand
    .58,   .10,   .32,   .02,   .065,    .41,    .123e-4,     7.5,   1.89,  .114,   & !  3 sandy loam
    .17,   .13,   .70,   .03,   .067,    .45,    .125e-5,     2.0,   1.41,  .179,   & !  4 silt loam
    .43,   .18,   .39,   .05,   .078,    .43,    .289e-5,     3.6,   1.56,  .155,   & !  5 loam
    .58,   .27,   .15,   .04,   .100,    .39,    .364e-5,     5.9,   1.48,  .175,   & !  6 sandy clay loam
    .10,   .34,   .56,   .06,   .089,    .43,    .194e-6,     1.0,   1.23,  .218,   & !  7 silty clay loam
    .32,   .34,   .34,   .07,   .095,    .41,    .722e-6,     1.9,   1.31,  .250,   & !  8 clay loam
    .52,   .42,   .06,   .08,   .100,    .38,    .333e-6,     2.7,   1.23,  .219,   & !  9 sandy clay
    .06,   .47,   .47,   .09,   .070,    .36,    .556e-7,      .5,   1.09,  .283,   & ! 10 silty clay
    .22,   .58,   .20,   .10,   .068,    .38,    .556e-6,      .8,   1.09,  .286,   & ! 11 clay 
    .10,   .06,   .84,   .05,   .034,    .46,    .694e-6,     1.6,   1.37,  .200/), & ! 12 silt
     (/10,12/) )

  real, parameter :: zero = 0.

  ! Vereecken PTF hydraulic parameter coefficients

  ! Sand, clay, and organ coefficient multipliers convert from fraction (kg/kg) to percentage.
  ! Bulkdens coefficient multiplier converts from kg/m^3 to g/cm^3.

  real, parameter :: a01 =  0.05            ! wresid intercept (Vereecken value 0.00, but need something > 0)
  real, parameter :: b01 =  zero    * 1.e2  ! wresid clay
  real, parameter :: c01 =  zero    * 1.e2  ! wresid organ

  real, parameter :: a02 =  0.6355          ! wsat intercept
  real, parameter :: b02 =  0.0013  * 1.e2  ! wsat clay
  real, parameter :: c02 = -0.1631  * 1.e-3 ! wsat bulkdens_drysoil

  real, parameter :: a03 = -4.3003          ! alpha intercept
  real, parameter :: b03 = -0.0097  * 1.e2  ! alpha clay 
  real, parameter :: c03 =  0.0138  * 1.e2  ! alpha sand
  real, parameter :: d03 =  zero    * 1.e-3 ! alpha bulkdens_drysoil
  real, parameter :: e03 = -0.0992  * 1.e2  ! alpha organ

  real, parameter :: a04 = -1.0846          ! en intercept
  real, parameter :: b04 = -0.0236  * 1.e2  ! en clay
  real, parameter :: c04 = -0.0085  * 1.e2  ! en sand
  real, parameter :: d04 =  0.0001  * 1.e4  ! en sand**2

  real, parameter :: a05 =  1.9582          ! ksat intercept
  real, parameter :: b05 =  zero    * 1.e2  ! ksat clay
  real, parameter :: c05 =  0.0308  * 1.e2  ! ksat sand
  real, parameter :: d05 = -0.6142  * 1.e-3 ! ksat bulkdens_drysoil
  real, parameter :: e05 = -0.1566  * 1.e2  ! ksat organ

  real, parameter :: a06 = -1.8642          ! lambda intercept
  real, parameter :: b06 = -0.1317  * 1.e2  ! lambda clay
  real, parameter :: c06 =  0.0067  * 1.e2  ! lambda sand
  real, parameter :: d06 =  zero    * 1.e-3 ! lambda bulkdens_drysoil
  real, parameter :: e06 =  zero    * 1.e2  ! lambda organ
 
  ! de Boer PTF hydraulic parameters

  ! Sand and clay multipliers convert from fraction to percentage.
  ! Organ multiplier converts from kg/kg to percentage.
  ! Bulkdens multiplier converts from kg/m^3 to g/cm^3.

  real, parameter :: a11a = 0.041, a11b = 0.179 ! wresid intercept

  real, parameter :: a12 =  0.83080           ! wsat intercept
  real, parameter :: b12 =  0.0002728 * 1.e2  ! wsat clay
  real, parameter :: c12 =  0.000187  * 1.e2  ! wsat silt 
  real, parameter :: d12 = -0.28217   * 1.e-3 ! wsat bulkdens_drysoil

  real, parameter :: a13 = -0.43348           ! alpha intercept
  real, parameter :: b13 = -0.01581   * 1.e2  ! alpha clay
  real, parameter :: c13 = -0.01207   * 1.e2  ! alpha silt
  real, parameter :: d13 = -0.41729   * 1.e-3 ! alpha bulkdens_drysoil
  real, parameter :: e13 = -0.04762   * 1.e2  ! alpha organ
  real, parameter :: f13 =  0.21810           ! alpha topsoil

  real, parameter :: a14 =  0.22236           ! en intercept
  real, parameter :: b14 = -0.005306  * 1.e2  ! en clay
  real, parameter :: c14 = -0.003084  * 1.e2  ! en silt
  real, parameter :: d14 = -0.30189   * 1.e-3 ! en bulkdens_drysoil
  real, parameter :: e14 = -0.01072   * 1.e2  ! en organ
  real, parameter :: f14 = -0.05558           ! en topsoil

  real, parameter :: a15 =  0.40220           ! ksat intercept
  real, parameter :: b15 = -0.02329   * 1.e2  ! ksat clay
  real, parameter :: c15 = -0.01265   * 1.e2  ! ksat silt
  real, parameter :: d15 = -0.01038           ! ksat cec_soil
  real, parameter :: e15 =  0.26122           ! ksat pH_soil
  real, parameter :: f15 =  0.44565           ! ksat topsoil

  ! Thermal parameters

  real, parameter :: specifheat_organic = 2.5e6 ! [J/(m^3 K)]
  real :: mineral
  real :: specifheat_mineral
  real :: topsoil
  real :: soilwilt

  real :: organ

  organ = organ0

!S  ! Limiter for organic content when Vereecken PTFs are used
!S  if (isoilptf == 1) organ = min(0.10,organ)
! Vereecken also suggested a possible need to limit clay content to 50% or so.

  mineral = 1. - organ

  if (isoilptf == 1) then

     ! Vereecken PTFs with zero coefficients

     ! wresid_vg =               a01 + b01 * clay + c01 * organ
     ! wsat_vg   =               a02 + b02 * clay + c02 * bulkdens_drysoil
     ! alpha_vg  =  -100. *  exp(a03 + b03 * clay + c03 * sand + d03 * bulkdens_drysoil + e03 * organ)
     ! en_vg     =      1. + exp(a04 + b04 * clay + c04 * sand + d04 * sand**2                       )
     ! ksat_vg   = .116e-6 * exp(a05 + b05 * clay + c05 * sand + d05 * bulkdens_drysoil + e05 * organ)
     ! lambda_vg =               a06 + b06 * clay + c06 * sand + d06 * bulkdens_drysoil + e06 * organ

     ! Vereecken PTFs without zero coefficients (alpha_vg converted from 1/cm to 1/m and multiplied
     ! by -1 to use with negative matric potential; ksat_vg converted from cm/day to m/s)

     wresid_vg =                 a01
     wsat_vg   =                 a02 + b02 * clay + c02 * bulkdens_drysoil
     alpha_vg  =     -100. * exp(a03 + b03 * clay + c03 * sand                          + e03 * organ)
     en_vg     =        1. + exp(a04 + b04 * clay + c04 * sand + d04 * sand**2                       )
     ksat_vg   = .11574e-6 * exp(a05              + c05 * sand + d05 * bulkdens_drysoil + e05 * organ)
     lambda_vg =                 a06 + b06 * clay + c06 * sand

  elseif (isoilptf == 2) then

     ! de Boer PTFs (alpha_vg converted from 1/cm to 1/m and multiplied by -1
     ! to use with negative matric potential; ksat_vg converted from cm/day to m/s)

     if (slzt < -0.30) then
        topsoil = 0. ! Topsoil flag
     else
        topsoil = 1. ! Topsoil flag
     endif

     if (sand > 2.0) then
        wresid_vg = a11a
     else
        wresid_vg = a11b
     endif

     wsat_vg   =                   a12 + b12 * clay + c12 * silt + d12 * bulkdens_drysoil
     alpha_vg  =     -100. * 10.**(a13 + b13 * clay + c13 * silt + d13 * bulkdens_drysoil + e13 * organ   + f13 * topsoil)
     en_vg     =        1. + 10.**(a14 + b14 * clay + c14 * silt + d14 * bulkdens_drysoil + e14 * organ   + f14 * topsoil)
     ksat_vg   = .11574e-6 * 10.**(a15 + b15 * clay + c15 * silt + d15 * cec_soil         + e15 * pH_soil + f15 * topsoil)
     lambda_vg = 0.5  

  elseif (isoilptf == 3) then

     ! van Genuchten parameters assigned from USDA soil textural class
     ! (alpha_vg multiplied by -1 to use with negative matric potential)

     wresid_vg =  soilparms4(5,usdatext)
     wsat_vg   =  soilparms4(6,usdatext)
     ksat_vg   =  soilparms4(7,usdatext)
     alpha_vg  = -soilparms4(8,usdatext)
     en_vg     =  soilparms4(9,usdatext)
     lambda_vg = 0.5  

     soilwilt  =  soilparms4(10,usdatext)

  endif

  ! Specific heat for dry soil material

  specifheat_mineral = (2.128e6 * sand + 2.385e6 * clay) / (sand + clay)
  specifheat_drysoil = (1. - wsat_vg) * (mineral * specifheat_mineral &
                                       + organ   * specifheat_organic)

end subroutine soil_ptf

