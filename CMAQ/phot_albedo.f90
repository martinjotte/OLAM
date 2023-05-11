subroutine phot_albedo(iw, coszens, currhr_lst, julian_day, jyfreq, &
                       albedo_dir, albedo_dif)

  use phot_mod
  use csqy_data
  use mem_ijtabs, only: itab_w
  use mem_sfcg,   only: itab_wsfc, sfcg
  use mem_sea,    only: sea, omsea
  use mem_land,   only: land, omland
  use mem_grid,   only: glatw, nsw_max, lpw
  use therm_lib,  only: qtk

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: coszens
  real,    intent( in) :: currhr_lst
  real,    intent( in) :: julian_day
  real,    intent( in) :: jyfreq
  real,    intent(out) :: albedo_dir(nsw_max,nwl)
  real,    intent(out) :: albedo_dif(nsw_max,nwl)

  integer, parameter :: lsea  = 17
  integer, parameter :: lice  = 20
  integer, parameter :: lsnow = 19
  integer            :: lland

   real    :: mscale
   real    :: seas_modulate
   real    :: zfactor
   real    :: sfactor
   real    :: tempk, fracliq, snowfac
   integer :: jsfc, iwsfc, jasfc, isea, iland
   integer :: ka, kw, ks
   real    :: arcoarkw
   real    :: albedo_sea_dir(nwl)
   real    :: albedo_sea_dif(nwl)
   real    :: albedo_ice_dir(nwl)
   real    :: albedo_ice_dif(nwl)
   real    :: albedo_land_dir(nwl)
   real    :: albedo_land_dif(nwl)
   real    :: albedo_snow_dir(nwl)
   real    :: albedo_snow_dif(nwl)


!     CMAQ albedo category             LEAF category
!----------------------------------------------------------------
!   1 EVERGREEN NEEDLE FOREST           0  Ocean
!   2 EVERGREEN BROADLEAF FOREST        1  Lakes, rivers, streams
!   3 DECIDUOUS NEEDLE FOREST           2  Ice cap/glacier
!   4 DECIDUOUS BROADLEAF FOREST        3  Desert, bare soil
!   5 MIXED FOREST                      4  Evergreen needleleaf tree
!   6 CLOSED SHRUBS                     5  Deciduous needleleaf tree
!   7 OPEN / SHRUBS                     6  Deciduous broadleaf tree
!   8 WOODY SAVANNA                     7  Evergreen broadleaf tree
!   9 SAVANNA                           8  Short grass
!  10 GRASSLAND                         9  Tall grass
!  11 PERMANENT WETLANDS               10  Semi-desert
!  12 CROPLAND                         11  Tundra
!  13 URBAN                            12  Evergreen shrub
!  14 CROP MOSAIC                      13  Deciduous shrub
!  15 PERMANENT SNOW                   14  Mixed woodland
!  16 BARREN / DESSERT                 15  Crop/mixed farming, C3 grassland
!  17 OCEAN WATER                      16  Irrigated crop
!  18 TUNDRA                           17  Bog or marsh
!  19 FRESH SNOW                       18  Wooded grassland
!  20 SEA ICE                          19  Urban and built up
!                                      20  Wetland evergreen broadleaf tree
!                                      21  Deforested (Amazon - Medvigy)

   integer, parameter :: leaf2cmaq(0:21) = (/ 17, 17, 15, 16,  1,  3,  4, 2, 10, 10, 7, &
                                              18,  6,  6,  5, 12, 12, 11, 8, 13,  2, 7 /)

  ! Determine seasonal and snow corrections to surface albedo
  ! convert julian into time of year for grid cell
  ! seasonal adjustment has an 11 day phase delay in the solar cycle

  if (glatw(iw) .ge. 0.0 ) then
     SEAS_MODULATE = COS( JYFREQ * ( JULIAN_DAY + CURRHR_LST / 24.0 + 11.0 ) )
  ELSE
     SEAS_MODULATE = COS( JYFREQ * ( JULIAN_DAY + CURRHR_LST / 24.0 + 11.0 ) + PI )
  END IF

  IF ( SEAS_MODULATE .GE. 0.0 ) THEN
     MSCALE = 0.5 * ( 1.0 + SQRT( SEAS_MODULATE ) )
  ELSE
     MSCALE = 0.5 * ( 1.0 - SQRT( abs(SEAS_MODULATE) ) )
  END IF

  albedo_dif(:,1:nwl) = 0.0
  albedo_dir(:,1:nwl) = 0.0

  ka = lpw(iw)

  ! Check for sea area beneath this atmospheric grid column

  if (itab_w(iw)%jsea2 > 0) then

     ! Loop over sea cells beneath this atmospheric grid column

     do jsfc = itab_w(iw)%jsea1, itab_w(iw)%jsea2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        arcoarkw = itab_wsfc(iwsfc)%arcoarkw(jasfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)

        isea = iwsfc - omsea

        ! water

        zfactor = MAX( 0.8, ( 1.0 + ZENITH_COEFF_REF( lsea ) )  &
                            / ( 1.0 + 2.0 * COSZENS * ZENITH_COEFF_REF( lsea ) ) )

        sfactor = 1.0 / (1.0 + MSCALE * (SEASON_COEFF_REF( lsea ) - 1.0))

        albedo_sea_dif(1:nwl) = sfactor * SPECTRAL_ALBEDO_REF(1:nwl, lsea)
        albedo_sea_dir(1:nwl) = zfactor * albedo_sea_dif(1:nwl)

        ! ice

        if (sea%nlev_seaice(isea) > 0) then

           zfactor = MAX( 0.8, ( 1.0 + ZENITH_COEFF_REF( lice ) )  &
                               / ( 1.0 + 2.0 * COSZENS * ZENITH_COEFF_REF( lice ) ) )

           sfactor = 1.0 / (1.0 + MSCALE * (SEASON_COEFF_REF( lice ) - 1.0))

           albedo_ice_dif(1:nwl) = sfactor * SPECTRAL_ALBEDO_REF(1:nwl, lice)
           albedo_ice_dir(1:nwl) = zfactor * albedo_ice_dif(1:nwl)

           albedo_ice_dif(1:nwl) = min(albedo_ice_dif(1:nwl), SPECTRAL_ALBEDO_REF(1:nwl,lsnow))
           albedo_ice_dir(1:nwl) = min(albedo_ice_dir(1:nwl), SPECTRAL_ALBEDO_REF(1:nwl,lsnow))

           albedo_sea_dif(1:nwl) = (1.0 - sea%seaicec(isea)) * albedo_sea_dif(1:nwl) &
                                 +        sea%seaicec(isea)  * albedo_ice_dif(1:nwl)

           albedo_sea_dir(1:nwl) = (1.0 - sea%seaicec(isea)) * albedo_sea_dir(1:nwl) &
                                 +        sea%seaicec(isea)  * albedo_ice_dir(1:nwl)
        endif

        ks = kw - ka + 1

        albedo_dif(ks,1:nwl) = albedo_dif(ks,1:nwl) + arcoarkw * albedo_sea_dif(1:nwl)
        albedo_dir(ks,1:nwl) = albedo_dir(ks,1:nwl) + arcoarkw * albedo_sea_dir(1:nwl)

     enddo
  endif

  ! Check for land area beneath this atmospheric grid column

  if (itab_w(iw)%jland2 > 0) then

     ! Loop over land cells beneath this atmospheric grid column

     do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        arcoarkw = itab_wsfc(iwsfc)%arcoarkw(jasfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)

        iland = iwsfc - omland

        lland = leaf2cmaq(sfcg%leaf_class(iwsfc))

        zfactor = MAX( 0.8, ( 1.0 + ZENITH_COEFF_REF( lland ) )  &
                            / ( 1.0 + 2.0 * COSZENS * ZENITH_COEFF_REF( lland ) ) )

        sfactor = 1.0 / (1.0 + MSCALE * (SEASON_COEFF_REF( lland ) - 1.0))

        albedo_land_dif(1:nwl) = sfactor * SPECTRAL_ALBEDO_REF(1:nwl, lland)
        albedo_land_dir(1:nwl) = zfactor * albedo_land_dif(1:nwl)

        if (land%nlev_sfcwater(iland) > 0) then
           call qtk(land%sfcwater_energy(land%nlev_sfcwater(iland),iland),tempk,fracliq)
           snowfac = (1.0 - fracliq) * land%snowfac(iland)
        else
           snowfac = 0.0
        endif

        if (snowfac >= 0.01) then

           zfactor = MAX( 0.8, ( 1.0 + ZENITH_COEFF_REF( lsnow ) )  &
                               / ( 1.0 + 2.0 * COSZENS * ZENITH_COEFF_REF( lsnow ) ) )

           sfactor = 1.0 / (1.0 + MSCALE * (SEASON_COEFF_REF( lsnow ) - 1.0))

           albedo_snow_dif(1:nwl) = sfactor * SPECTRAL_ALBEDO_REF(1:nwl, lsnow)
           albedo_snow_dir(1:nwl) = zfactor * albedo_snow_dif(1:nwl)

           albedo_snow_dif(1:nwl) = min(albedo_snow_dif(1:nwl), SPECTRAL_ALBEDO_REF(1:nwl,lsnow))
           albedo_snow_dir(1:nwl) = min(albedo_snow_dir(1:nwl), SPECTRAL_ALBEDO_REF(1:nwl,lsnow))

           albedo_land_dif(1:nwl) = (1.0 - snowfac) * albedo_land_dif(1:nwl) &
                                    +        snowfac  * albedo_snow_dif(1:nwl)

           albedo_land_dir(1:nwl) = (1.0 - snowfac) * albedo_land_dir(1:nwl) &
                                   +        snowfac  * albedo_snow_dir(1:nwl)
        endif

        ks = kw - ka + 1

        albedo_dif(ks,1:nwl) = albedo_dif(ks,1:nwl) + arcoarkw * albedo_land_dif(1:nwl)
        albedo_dir(ks,1:nwl) = albedo_dir(ks,1:nwl) + arcoarkw * albedo_land_dir(1:nwl)

     enddo
  endif

end subroutine phot_albedo
