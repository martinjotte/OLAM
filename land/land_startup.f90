subroutine land_startup()

  use leaf_coms,   only: nzs, ndviflg, iupdndvi, isoilflg, specifheat_bedrock
  use mem_land,    only: alloc_land2, filltab_land, land, mland, nzg, slzt, omland
  use misc_coms,   only: runtype
  use mem_sfcg,    only: sfcg
  use consts_coms, only: cice
  use leaf4_soil,  only: soil_pot2wat, calc_wfrac_low

  implicit none

  integer :: iland, k, iwsfc
  real    :: psi

  ! Initialize some land & vegetation properties

  call land_parms()

  ! Allocate time-dependent LEAF arrays and add to history file I/O table

  call alloc_land2 ()
  call filltab_land()

  ! Fill ndvi values

  if ( (iupdndvi /= 1) .and. &
       (runtype == 'HISTORY' .or. runtype == 'HISTREGRID') ) then

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
          runtype == 'HISTREGRID') then

     ! Make inventory of ndvi files and initialize current/past time
     ! Not needed for a PLOTONLY run.

     call ndvi_database_read(0)

     ! Initialize future ndvi file time

     call ndvi_database_read(1)

  endif

  ! Apply pedotransfer functions to obtain soil hydraulic properties and
  ! specific heat.  Loop over all land points, INCLUDING THOSE THAT ARE NOT
  ! PRIMARY ON THIS SUBDOMAIN.

  !$omp parallel do private(iwsfc,k,psi)
  do iland = 2,mland
     iwsfc = iland + omland

     if (sfcg%leaf_class(iwsfc) == 2) then

        ! If the landuse type of this land cell is glacier/firn/ice cap,
        ! set some physical properties related to heat and water content
        ! that are more appropriate for firn.

        do k = 1, nzg
           land%wsat_vg           (k,iland) = 0.1 ! 10% porosity assumed (limits soil water content)
           land%wresid_vg         (k,iland) = land%wsat_vg(k,iland) * 0.2 ! wresid must be < wsat
           land%ksat_vg           (k,iland) = 0.  ! prevents water fluxes
           land%alpha_vg          (k,iland) = -2.0
           land%en_vg             (k,iland) = 1.4
           land%lambda_vg         (k,iland) = 0.5
           land%k_bedrock           (iland) = 0

           land%specifheat_drysoil(k,iland) = cice * 600. ! Assumes firn density of 600 kg/m^3
           land%bulkdens_drysoil  (k,iland) = 600.

           ! Default values for firn; will not be used in the code
           land%sand              (k,iland) = 0.
           land%clay              (k,iland) = 0.
           land%silt              (k,iland) = 1.
           land%organ             (k,iland) = 0.
           land%ph_soil           (k,iland) = 7.
           land%cec_soil          (k,iland) = 20.

        enddo

     else

        ! z_bedrock, the height of the upper boundary of bedrock, is used for
        ! the transition grid level when bedrock is near the surface.  However,
        ! where bedrock begins much deeper, the transition level is set above
        ! z_bedrock because SoilGrids data is defined only in the top 2 m, while
        ! GLHYMPS applies to roughly the top 100 m.  For now, we choose the
        ! transition level to be no greater than 10 meters below the surface.

        do k = nzg, 1, -1
           if (slzt(k) < max(-10.0, land%z_bedrock(iland))) exit
        enddo
        land%k_bedrock(iland) = k

        ! Since SoilGrids data is only applicable to the top 2 meters,
        ! GLHYMPS permeability and porosity will be used below the
        ! transition level to set wsat and ksat.
        !
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

        do k = 1, land%k_bedrock(iland)

           land%wresid_vg         (k,iland) = 0.2 * land%glhymps_poros(iland)
           land%wsat_vg           (k,iland) =       land%glhymps_poros(iland)
           land%ksat_vg           (k,iland) =       land%glhymps_ksat (iland)
           land%alpha_vg          (k,iland) = -2.0
           land%en_vg             (k,iland) = 1.4
           land%lambda_vg         (k,iland) = 0.5

           land%specifheat_drysoil(k,iland) = specifheat_bedrock
           land%bulkdens_drysoil  (k,iland) = 2700. * (1. - land%glhymps_poros(iland))

           ! Default values for bedrock; will not be used in the code
           land%sand              (k,iland) = 0.
           land%clay              (k,iland) = 0.
           land%silt              (k,iland) = 1.
           land%organ             (k,iland) = 0.
           land%ph_soil           (k,iland) = 7.
           land%cec_soil          (k,iland) = 20.

        enddo

        ! Above transition level apply soil pedotransfer functions

        do k = land%k_bedrock(iland)+1, nzg

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

        enddo

     endif

     ! Compute wfrac_low, the minimum soil water fraction beyond which we use
     ! a linear fit to the van Genuchten soil water retention curve. Wfrac_low
     ! is the water fraction that gives d(head)/d(vol_wat_fraction) equal to
     ! psi_slope_low, which is set to a finite value for numerical stability.

     call calc_wfrac_low(nzg, land%wfrac_low(:,iland), land%wsat_vg (:,iland), &
                              land%wresid_vg(:,iland), land%en_vg   (:,iland), &
                              land%alpha_vg (:,iland)  )

     ! Estimate top-layer soil field capacity as the soil moisture
     ! that gives 220 cm of head for sand and 340 cm for silt/clay

     psi = -2.2*land%sand(nzg,iland) - 3.4*(1.-land%sand(nzg,iland))
     call soil_pot2wat( psi, land%wresid_vg(nzg,iland), land%wsat_vg(nzg,iland), &
                             land%alpha_vg (nzg,iland), land%en_vg  (nzg,iland), &
                             land%wfrac_low(nzg,iland), land%soilfldcap (iland)  )

     ! Estimate top-layer soil wilting point as the soil moisture
     ! that gives 1500 kPa (153m) of hydraulic head

     psi = -153.0
     call soil_pot2wat( psi, land%wresid_vg(nzg,iland), land%wsat_vg(nzg,iland), &
                             land%alpha_vg (nzg,iland), land%en_vg  (nzg,iland), &
                             land%wfrac_low(nzg,iland), land%soilwilt   (iland)  )
  enddo

end subroutine land_startup

!==========================================================================

subroutine land_parms()

  use leaf_coms,  only: nvtyp, albv_green, albv_brown, emisv,         &
                        sr_max, tai_max, sai, veg_clump, veg_frac,    &
                        veg_ht, dead_frac, rcmin, glai_max, dfpardsr, &
                        fpar_max, fpar_min, sr_min, z_root, kroot,    &
                        snowmin_expl, wcap_min, wcap_vmin, dt_leaf

  use mem_land,   only: nzg, slz, kperc

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

  endif

end subroutine land_parms

!==========================================================================

subroutine soil_ptf(iland, k, slzt, usdatext, sand, clay, silt,   &
                    organ0, bulkdens_drysoil, cec_soil, pH_soil,  &
                    wresid_vg, wsat_vg, alpha_vg, en_vg, ksat_vg, &
                    lambda_vg, specifheat_drysoil                 )

  use leaf_coms, only: specifheat_coarse, specifheat_fine, specifheat_organic, &
                       isoilptf, isoilflg

  ! Computes static (time-independent) thermal and hydraulic parameters for
  ! soil layers, but not for bedrock layers.

  implicit none

  integer, intent(in)  :: iland, k, usdatext
  real,    intent(in)  :: slzt, sand, clay, silt, organ0
  real,    intent(in)  :: bulkdens_drysoil, cec_soil, pH_soil
  real,    intent(out) :: wresid_vg, wsat_vg, alpha_vg, en_vg, ksat_vg, &
                          lambda_vg, specifheat_drysoil

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

  real :: mineral
  real :: specifheat_mineral
  real :: topsoil

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

     if (sand > 0.02) then
        wresid_vg = a11a
     else
        wresid_vg = a11b
     endif

     wsat_vg   =                   a12 + b12 * clay + c12 * silt + d12 * bulkdens_drysoil
     alpha_vg  =     -100. * 10.**(a13 + b13 * clay + c13 * silt + d13 * bulkdens_drysoil + e13 * organ   + f13 * topsoil)
     en_vg     =        1. + 10.**(a14 + b14 * clay + c14 * silt + d14 * bulkdens_drysoil + e14 * organ   + f14 * topsoil)
     ksat_vg   = .11574e-6 * 10.**(a15 + b15 * clay + c15 * silt + d15 * cec_soil         + e15 * pH_soil + f15 * topsoil)
     lambda_vg = 0.5

  endif

  ! Specific heat for dry soil material

  specifheat_mineral = specifheat_coarse * sand + specifheat_fine * (1. - sand)
  specifheat_drysoil = (1. - wsat_vg) * ( mineral * specifheat_mineral &
                                        + organ   * specifheat_organic )

end subroutine soil_ptf
