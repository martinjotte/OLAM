subroutine seacells(isea, timefac_sst, timefac_seaice)

  use sea_coms,    only: nzi, fssat0
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, omsea
  use misc_coms,   only: iparallel
  use consts_coms, only: grav, p00i, rocp, eps_virt, cpi, erad, eradi, cliq1000
  use mem_para,    only: myrank
  use oname_coms,  only: nl
  use pom2k1d,     only: pom, rhoref, pom_column
  use umwm_module, only: umwm

  implicit none

  integer, intent(in) :: isea
  real,    intent(in) :: timefac_sst    ! fraction of elapsed time from past to future SST obs
  real,    intent(in) :: timefac_seaice ! fraction of elapsed time from past to future SEA ICE obs

  ! Local variables

  integer :: iwsfc, ispeed10
  logical :: ispray_active, iumwm_active, use_umwm_roughness

  real :: canexneri, cantheta, canthetav
  real :: airthetav, wstar
  real :: usti, zw, zn1, zn2
  real :: raxis, windu, windv, cdtop, wusurf, wvsurf, wtsurf, wssurf, swrad, hfluxsea
  real :: sea_spray1_temp, sea_spray2_temp, fssat
  real :: vels0, richnum, a2fm, a2fh, speed10

  real, parameter :: z0fac = 0.011 / grav  ! factor for Charnock roughness height
  real, parameter :: ozo   = 1.59e-5       ! base roughness height in HWRF
  real, parameter :: onethird = 1./3.
  real, parameter :: ubmin = .1            ! lower bound on wind speed [m/s]
  real, parameter :: z10   = 10.           ! Height for evaluating 10m wind speed [m]

  ! Define variable roughness length from Breivik et al. 2022, Figure 3b,
  ! Ekofisk: ALT, which is based on winds at z = 10 m and spans wind speeds
  ! from 1 to 31 m/s.  At wind speeds below 15 m/s, drag coefficient is too
  ! close to zero to be easily discernible in the figure, so the much more
  ! sensitive Charnock parameter in Figure 3c is used instead to infer the
  ! roughness length.  Then, extrapolate downward-trending
  ! curve to 35 m/s, where the roughness length equals 0.0030 and hold
  ! roughness constant at 0.0030 for higher wind speeds.

  real, parameter :: rough_br(35) = &
     [.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0002,.0004, &
      .0006,.0008,.0011,.0015,.0021,.0030,.0045,.0060,.0075,.0090, &
      .0120,.0150,.0170,.0180,.0250,.0275,.0285,.0320,.0310,.0260, &
      .0205,.0150,.0095,.0040,.0030 ]

  ! Define variable roughness length from Breivik et al. 2022, Figure 3b,
  ! Ekofisk: ALT, which is based on winds at z = 10 m and spans wind speeds
  ! from 1 to 31 m/s.  At wind speeds below 15 m/s, drag coefficient is too
  ! close to zero to be easily discernible in the figure, so the much more
  ! sensitive Charnock parameter in Figure 3c is used instead to infer the
  ! roughness length.  Then, extrapolate downward-trending
  ! curve to 33 m/s, where the roughness length equals 0.0100 and hold
  ! roughness constant at 0.0100 for higher wind speeds.

  real, parameter :: rough_br2(35) = &
     [.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0002,.0004, &
      .0006,.0008,.0011,.0015,.0021,.0030,.0045,.0060,.0075,.0090, &
      .0120,.0150,.0170,.0180,.0250,.0275,.0285,.0320,.0310,.0260, &
      .0205,.0150,.0100,.0100,.0100 ]

  ! Define variable roughness length from Breivik et al. 2022, Figure 3b,
  ! Ekofisk: ALT, which is based on winds at z = 10 m and spans wind speeds
  ! from 1 to 31 m/s.  At wind speeds below 15 m/s, drag coefficient is too
  ! close to zero to be easily discernible in the figure, so the much more
  ! sensitive Charnock parameter in Figure 3c is used instead to infer the
  ! roughness length.  Increase roughness value at 23, 24, and 27 m/s to remove
  ! dips in the curve.  At wind speeds above 28 m/s, which is where roughness
  ! reaches a maximum value of .0320, fix roughness at that same value.

  real, parameter :: rough_br3(35) = &
     [.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0002,.0004, &
      .0006,.0008,.0011,.0015,.0021,.0030,.0045,.0060,.0075,.0090, &
      .0120,.0150,.0180,.0215,.0250,.0275,.0300,.0320,.0320,.0320, &
      .0320,.0320,.0320,.0320,.0320 ]

  iwsfc = isea + omsea

  ! Skip this cell if running in parallel and cell rank is not MYRANK

  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  ! Zero runoff for sea cells

  sfcg%runoff(iwsfc) = 0.

  ! Update seaice fraction

  sea%seaicec(isea) = sea%seaicep(isea) &
                    + timefac_seaice * (sea%seaicef(isea) - sea%seaicep(isea))

  ! Update SEATC if isea cell is not pom_active

  if (.not. sea%pom_active(isea)) then

     sea%seatc(isea) = sea%seatp(isea) &
                     + timefac_sst * (sea%seatf(isea) - sea%seatp(isea))
  endif

  ! Prepare to evaluate surface layer exchange parameters for open water areas

  airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))
  canexneri = 1. / sfcg%canexner(iwsfc)

  cantheta  = sea%sea_cantemp(isea) * canexneri
  canthetav = cantheta * (1.0 + eps_virt * sea%sea_canrrv(isea))

  wstar = (grav * sfcg%pblh(iwsfc) * max(sea%sea_wthv(isea),0.0) / airthetav) ** onethird

  ! If the roughness length from the UMWM will be used

  iumwm_active = .false.
  if ( allocated(umwm%iactive) ) iumwm_active = umwm%iactive(isea)

  use_umwm_roughness = .false.
  if ( iumwm_active) use_umwm_roughness = ( nl%use_umwm_roughness == 1  .and.        &
                                            umwm%wspd(isea) > nl%umwm_wind_threshold )

  ! sea water saturation reduction factor

  fssat = fssat0

  if (sea%pom_active(isea) .and. nl%sea_salinity_effect == 1) then
     fssat = 1.0 - (.02 / 35.) * pom%salin(1,isea)
  endif

  ! Evaluate sea roughness height

  if ( use_umwm_roughness ) then

     sea%sea_rough(isea) = max(0.0001,exp(umwm%alogzo(isea)))

  elseif (nl%iroughsea == 1) then

     ! Charnock (1955):
     ! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnock (1955)

     ! Davis et al. (2008) originally used in HWRF
     ! rough = 10. * exp(-10. / ustar ** .333333)
     ! rough = max(.125e-6, min(2.85e-3,rough))

     ! 2012 HWRF scheme; interpolates between the Charnock scheme at low wind
     ! and the Davis et al. curve fit at high wind speeds

     usti  = 1.0 / sea%sea_ustar(isea)
     zw    = min( (sea%sea_ustar(isea)/1.06)**onethird, 1.0 )
     zn1   = z0fac * sea%sea_ustar(isea) * sea%sea_ustar(isea) + ozo
     zn2   = 10. * exp(-9.5 * usti**onethird) + 1.65e-6 * usti

     sea%sea_rough(isea) = min( (1.0-zw) * zn1 + zw * zn2, 2.85e-3 )

  elseif (nl%iroughsea == 2) then

     ! Charnock (1955) with 0.018 parameter:
     sea%sea_rough(isea) = max(0.018 / grav * sea%sea_ustar(isea) ** 2,.0001)

  elseif (nl%iroughsea == 3) then

     ! Variable roughness length from Breivik et al. (2022)

     ! First, evaluate wind at 10 m based on current roughness, ustar,
     ! ATM, and canopy values

     vels0 = max(ubmin, sqrt(sfcg%vels(iwsfc)*sfcg%vels(iwsfc) + wstar*wstar))

     richnum = 2.0 * grav * sfcg%dzt_bot(iwsfc) * (airthetav - canthetav)  &
             / ( (airthetav + canthetav) * vels0 * vels0 )

     ! Get the Louis drag coefficients defined for wind at 10 m height
     call a2fmfh(richnum, z10, sea%sea_rough(isea), a2fm, a2fh)

     speed10 = sea%sea_ustar(isea) / sqrt(a2fm)
     ispeed10 = int(speed10)

     if (ispeed10 == 0) then
        sea%sea_rough(isea) = rough_br(1)
     elseif (ispeed10 >= 35) then
        sea%sea_rough(isea) = rough_br(35)
     else
        sea%sea_rough(isea) = rough_br(ispeed10)  &
             + (speed10 - real(ispeed10)) * (rough_br(ispeed10+1) - rough_br(ispeed10))
     endif

  elseif (nl%iroughsea == 4) then

     ! Variable roughness length from Breivik et al. (2022) [version 2]

     ! First, evaluate wind at 10 m based on current roughness, ustar,
     ! ATM, and canopy values

     vels0 = max(ubmin, sqrt(sfcg%vels(iwsfc)*sfcg%vels(iwsfc) + wstar*wstar))

     richnum = 2.0 * grav * sfcg%dzt_bot(iwsfc) * (airthetav - canthetav)  &
             / ( (airthetav + canthetav) * vels0 * vels0 )

     ! Get the Louis drag coefficients defined for wind at 10 m height
     call a2fmfh(richnum, z10, sea%sea_rough(isea), a2fm, a2fh)

     speed10 = sea%sea_ustar(isea) / sqrt(a2fm)
     ispeed10 = int(speed10)

     if (ispeed10 == 0) then
        sea%sea_rough(isea) = rough_br2(1)
     elseif (ispeed10 >= 35) then
        sea%sea_rough(isea) = rough_br2(35)
     else
        sea%sea_rough(isea) = rough_br2(ispeed10)  &
             + (speed10 - real(ispeed10)) * (rough_br2(ispeed10+1) - rough_br2(ispeed10))
     endif

  elseif (nl%iroughsea == 5) then

     ! Variable roughness length from Breivik et al. (2022) [version 3]

     ! First, evaluate wind at 10 m based on current roughness, ustar,
     ! ATM, and canopy values

     vels0 = max(ubmin, sqrt(sfcg%vels(iwsfc)*sfcg%vels(iwsfc) + wstar*wstar))

     richnum = 2.0 * grav * sfcg%dzt_bot(iwsfc) * (airthetav - canthetav)  &
             / ( (airthetav + canthetav) * vels0 * vels0 )

     ! Get the Louis drag coefficients defined for wind at 10 m height
     call a2fmfh(richnum, z10, sea%sea_rough(isea), a2fm, a2fh)

     speed10 = sea%sea_ustar(isea) / sqrt(a2fm)
     ispeed10 = int(speed10)

     if (ispeed10 == 0) then
        sea%sea_rough(isea) = rough_br3(1)
     elseif (ispeed10 >= 35) then
        sea%sea_rough(isea) = rough_br3(35)
     else
        sea%sea_rough(isea) = rough_br3(ispeed10)  &
             + (speed10 - real(ispeed10)) * (rough_br3(ispeed10+1) - rough_br3(ispeed10))
     endif

  else
     print*, 'invalid iroughsea value '
     stop 'stop: iroughsea '
  endif

  call stars(sfcg%dzt_bot  (iwsfc), &
             sea%sea_rough  (isea), &
             sfcg%vels     (iwsfc), &
             sfcg%rhos     (iwsfc), &
             wstar                , &
             airthetav            , &
             canthetav            , &
             sea%sea_vkmsfc (isea), &
             sea%sea_vkhsfc (isea), &
             sea%sea_ustar  (isea), &
             sea%sea_ggaer  (isea)  )

! When we have CO2:
!  sea%sea_sfluxc(isea) = sfcg%rhos(iwsfc) * sea%sea_ggaer(isea) &
!                       * (sea%sea_co2(isea) - air_co2)

  sea_spray1_temp = sea%seatc(isea)
  sea_spray2_temp = sea%seatc(isea)
  ispray_active   = .false.

  if (allocated(sea%spray_active)) ispray_active = sea%spray_active(isea)

  if (ispray_active) then
     if (allocated(sea%spraytemp )) sea_spray1_temp = sea%spraytemp (isea)
     if (allocated(sea%spray2temp)) sea_spray2_temp = sea%spray2temp(isea)
  endif

  ! Update SEA fields

  call seacell_2(isea, iwsfc,             &
                 sfcg%rhos      (iwsfc),  &
                 sea%sea_ustar   (isea),  &
                 sea%sea_vkhsfc  (isea),  &
                 sfcg%can_depth (iwsfc),  &
                 sea%seatc       (isea),  &
                 sfcg%vels      (iwsfc),  &
                 wstar                 ,  &
                 sfcg%prss      (iwsfc),  &
                 sfcg%glatw     (iwsfc),  &
                 sfcg%glonw     (iwsfc),  &
                 sfcg%airtheta  (iwsfc),  &
                 sfcg%airrrv    (iwsfc),  &
                 sfcg%canexner  (iwsfc),  &
                 sea%sea_cantemp (isea),  &
                 sea%sea_canrrv  (isea),  &
                 sea%sea_bcantemp(isea),  &
                 sea%sea_bcanrrv (isea),  &
                 sea_spray1_temp,         &
                 sea_spray2_temp,         &
                 sea%sea_sfluxt  (isea),  &
                 sea%sea_sfluxr  (isea),  &
                 sea%sea_rough   (isea),  &
                 hfluxsea,                &
                 fssat,                   &
                 ispray_active            )

  if (allocated(sea%spraytemp ))   sea%spraytemp   (isea) = sea_spray1_temp
  if (allocated(sea%spray2temp))   sea%spray2temp  (isea) = sea_spray2_temp
  if (allocated(sea%spray_active)) sea%spray_active(isea) = ispray_active

  ! New calculation of wthv for sfluxt units [W m^-2]

  sea%sea_wthv(isea) = ( sea%sea_sfluxt(isea) * cpi * canexneri * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
                     + sea%sea_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

  ! Update sea ice based on seaice fraction

  call prep_seaice(sea%seatc              (isea), &
                   sea%seaicec            (isea), &
                   sea%sea_cantemp        (isea), &
                   sea%ice_cantemp        (isea), &
                   sea%seaice_energy(1:nzi,isea), &
                   sea%seaice_tempk (1:nzi,isea), &
                   sea%nlev_seaice        (isea), &
                   sea%ice_albedo         (isea), &
                   sea%ice_rlongup        (isea), &
                   sea%ice_net_rshort     (isea), &
                   sea%ice_net_rlong      (isea), &
                   sfcg%rshort           (iwsfc), &
                   sfcg%rlong            (iwsfc), &
                   sea%ice_rough          (isea), &
                   sea%sea_canrrv         (isea), &
                   sea%ice_canrrv         (isea), &
                   sea%sea_ustar          (isea), &
                   sea%ice_ustar          (isea), &
                   sea%sea_ggaer          (isea), &
                   sea%ice_ggaer          (isea), &
                   sea%sea_wthv           (isea), &
                   sea%ice_wthv           (isea)  )

  ! Check if sea ice is present

  if (sea%nlev_seaice(isea) == 0) then

     ! If no sea ice is present in this cell, zero out fluxes over ice
     ! and set cell flux to open water part

     sfcg%vkmsfc(iwsfc) = sea%sea_vkmsfc(isea)
     sfcg%ustar (iwsfc) = sea%sea_ustar (isea)
     sfcg%ggaer (iwsfc) = sea%sea_ggaer (isea)
     sfcg%sfluxt(iwsfc) = sea%sea_sfluxt(isea)
     sfcg%sfluxr(iwsfc) = sea%sea_sfluxr(isea)
!    sfcg%sfluxc(iwsfc) = sea%sea_sfluxc(isea)
     sfcg%wthv  (iwsfc) = sea%sea_wthv  (isea)

     sea%ice_vkmsfc(isea) = 0.0
     sea%ice_ustar (isea) = 0.0
     sea%ice_ggaer (isea) = 0.0
     sea%ice_sfluxt(isea) = 0.0
     sea%ice_sfluxr(isea) = 0.0
!    sea%ice_sfluxc(isea) = 0.0

     sfcg%rough      (iwsfc) = sea%sea_rough   (isea)
     sfcg%cantemp    (iwsfc) = sea%sea_cantemp (isea)
     sfcg%canrrv     (iwsfc) = sea%sea_canrrv  (isea)

  else

     ! If sea ice is present in this cell, compute turbulent fluxes over ice

     cantheta  = sea%ice_cantemp(isea) * canexneri
     canthetav = cantheta * (1.0 + eps_virt * sea%ice_canrrv(isea))

     wstar = (grav * sfcg%pblh(iwsfc) * max(sea%ice_wthv(isea),0.0) / airthetav) ** onethird

     call stars(sfcg%dzt_bot  (iwsfc), &
                sea%ice_rough  (isea), &
                sfcg%vels     (iwsfc), &
                sfcg%rhos     (iwsfc), &
                wstar                , &
                airthetav            , &
                canthetav            , &
                sea%ice_vkmsfc (isea), &
                sea%ice_vkhsfc (isea), &
                sea%ice_ustar  (isea), &
                sea%ice_ggaer  (isea)  )

! When we have CO2:
!    sea%ice_sfluxc(isea) = sfcg%rhos(iwsfc) * sea%ice_ggaer(isea) &
!                         * (sea%ice_co2(isea) - airco2)

     ! Compute seaice canopy fluxes

     call seaice(isea, iwsfc,                   &
                 sea%nlev_seaice        (isea), &
                 sfcg%rhos             (iwsfc), &
                 sea%ice_ustar          (isea), &
                 sea%ice_vkhsfc         (isea), &
                 sfcg%can_depth        (iwsfc), &
                 sea%ice_net_rshort     (isea), &
                 sea%ice_net_rlong      (isea), &
                 sfcg%glatw            (iwsfc), &
                 sfcg%glonw            (iwsfc), &
                 sfcg%airtheta         (iwsfc), &
                 sfcg%airrrv           (iwsfc), &
                 sfcg%canexner         (iwsfc), &
                 sea%seaice_energy(1:nzi,isea), &
                 sea%seaice_tempk (1:nzi,isea), &
                 sea%ice_cantemp        (isea), &
                 sea%ice_canrrv         (isea), &
                 sea%ice_sfluxt         (isea), &
                 sea%ice_sfluxr         (isea)  )

! Original calculation of wthv for sfluxt units [kg_dry K m^-2 s^-1]

!     sea%ice_wthv(isea) = ( sea%ice_sfluxt(isea) * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
!        + sea%ice_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

! New calculation of wthv for sfluxt units [W m^-2]

     sea%ice_wthv(isea) = ( sea%ice_sfluxt(isea) * cpi * canexneri * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
          + sea%ice_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

     ! Combine sea and ice values based on ice fraction:

     sfcg%rough    (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_rough  (isea) + &
                                    sea%seaicec(isea)  * sea%ice_rough  (isea)

     sfcg%cantemp  (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_cantemp(isea) + &
                                    sea%seaicec(isea)  * sea%ice_cantemp(isea)

     sfcg%canrrv   (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_canrrv (isea) + &
                                    sea%seaicec(isea)  * sea%ice_canrrv (isea)

     sfcg%vkmsfc(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_vkmsfc(isea) &
                               + sea%seaicec(isea)  * sea%ice_vkmsfc(isea)

     sfcg%ustar(iwsfc)  = (1.0 - sea%seaicec(isea)) * sea%sea_ustar(isea) &
                               + sea%seaicec(isea)  * sea%ice_ustar(isea)

     sfcg%ggaer(iwsfc)  = (1.0 - sea%seaicec(isea)) * sea%sea_ggaer(isea) &
                               + sea%seaicec(isea)  * sea%ice_ggaer(isea)

     sfcg%sfluxt(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxt(isea) &
                               + sea%seaicec(isea)  * sea%ice_sfluxt(isea)

     sfcg%sfluxr(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxr(isea) &
                               + sea%seaicec(isea)  * sea%ice_sfluxr(isea)

     sfcg%wthv(iwsfc)   = (1.0 - sea%seaicec(isea)) * sea%sea_wthv(isea) &
                               + sea%seaicec(isea)  * sea%ice_wthv(isea)

!    sfcg%sfluxc(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_sfluxc(isea) &
!                              + sea%seaicec(isea)  * sea%ice_sfluxc(isea)

  endif

  ! Update POM1D vertical column variables if this sea cell is pom_active

  if (sea%pom_active(isea)) then

     ! Eastward and northward surface wind components for sea cell

     raxis = sqrt(sfcg%xew(iwsfc) ** 2 + sfcg%yew(iwsfc) ** 2)  ! dist from earth axis

     if (raxis > 1.e-3) then
        windu = (sea%windye(isea) * sfcg%xew(iwsfc) &
               - sea%windxe(isea) * sfcg%yew(iwsfc)) / raxis

        windv =  sea%windze(isea) * raxis * eradi    &
              - (sea%windxe(isea) * sfcg%xew(iwsfc)  &
               + sea%windye(isea) * sfcg%yew(iwsfc)) &
              * sfcg%zew(iwsfc) / (raxis * erad)
     else
        windu = 0.
        windv = 0.
     endif

     ! CDTOP is the drag coefficient between wind and water at the top water surface
     ! and is based on vkmsfc computed in subroutine stars for surface wind stress.

     cdtop = sfcg%vkmsfc(iwsfc) / (sfcg%dzt_bot(iwsfc) * rhoref)

     wusurf = cdtop * (pom%ub(1,isea) - windu)
     wvsurf = cdtop * (pom%vb(1,isea) - windv)
     wtsurf = hfluxsea / rhoref + (sfcg%rlong(iwsfc) - sfcg%rlongup(iwsfc)) / cliq1000
     wssurf = 0.
     swrad = sfcg%rshort(iwsfc) * (1. - sfcg%albedo_beam(iwsfc)) / cliq1000

     call pom_column(isea, sea%pom_kba(isea), wusurf, wvsurf, wtsurf, wssurf, swrad)

     sea%seatc(isea) = pom%potmp(1,isea) + 273.15

  endif

end subroutine seacells

!===============================================================================

subroutine seacell_2( isea, iwsfc, rhos, ustar, vkhsfc, can_depth, seatc, vels, &
                      wstar, prss, glatw, glonw, airtheta, airrrv, canexner, &
                      cantemp, canrrv, bcantemp, bcanrrv, spraytemp, spray2temp, &
                      sfluxt, sfluxr, rough, hfluxsea, fssat, spray_active )

  use mem_sfcg,    only: sfcg
  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav, alvl, cliq, r8, p00i, rocp, eps_virt
  use therm_lib,   only: rhovsl
  use matrix,      only: matrix8_NxN
  use leaf4_canopy,only: sing_print
  use oname_coms,  only: nl
  use umwm_module, only: umwm, swh

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
  integer, intent(in)    :: iwsfc        ! current sfcg cell index
  real,    intent(in)    :: rhos         ! air density [kg_dryair/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: vkhsfc       ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: seatc        ! current sea temp (obs time) [K]
  real,    intent(in)    :: vels         ! wind speed at top of sfc layer [m/s]
  real,    intent(in)    :: wstar        ! PBL convective velocity scale [m/s]
  real,    intent(in)    :: prss         ! canopy air pressure [Pa]
  real,    intent(in)    :: glatw        ! Latitude of land cell 'center' [deg]
  real,    intent(in)    :: glonw        ! Longitude of land cell 'center' [deg]
  real,    intent(in)    :: airtheta     ! atm potential temp [K]
  real,    intent(in)    :: airrrv       ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: canexner     ! canopy Exner function []
  real,    intent(inout) :: cantemp      ! can_air temp [K]
  real,    intent(inout) :: canrrv       ! can_air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(inout) :: bcantemp     ! bcan_air temp [K]
  real,    intent(inout) :: bcanrrv      ! bcan_air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(inout) :: spraytemp    ! seaspray temperature [K]
  real,    intent(inout) :: spray2temp   ! seaspray2 temperature [K]
  real,    intent(inout) :: sfluxt       ! can_air to atm sens heat flux [W m^-2]
  real,    intent(inout) :: sfluxr       ! can_air to atm vap flux [kg_vap m^-2 s^-1]
  real,    intent(in)    :: rough        ! sea cell roughess height [m]
  real,    intent(out)   :: hfluxsea     ! heat flux from sea to can_air [kg K/(m^2 s)]
  real,    intent(in)    :: fssat        ! sea water saturation reduction factor []
  logical, intent(inout) :: spray_active ! was sea spray layer active previous timestep

  ! Local parameters

  real, parameter :: one3    = 1./ 3.
  real, parameter :: fcn     = 0.75      ! Crank-Nicolson future time weight
  real, parameter :: ubmin   = .1        ! lower bound on wind speed
  real, parameter :: ubminsq = ubmin**2

  ! Local variables

  real :: rdi            ! sea surface to can_air conductance [m/s]
  real :: hxfersc        ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca        ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca        ! vapor xfer from can_air to atm this step [kg/m^2]

  real :: hxfercb   !    ! heat xfer from can_air to bcanopy this step [J/m^2]
  real :: wxfercb   !    ! vapor xfer from can_air to bcanopy this step [kg/m^2]
  real :: hxferba   !    ! heat xfer from bcanopy to atm this step [J/m^2]
  real :: wxferba   !    ! vapor xfer from bcanopy to atm this step [kg/m^2]

  real :: wxfersc        ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]
  real :: hxfersd        ! heat xfer from sea surface to seaspray drops this step [J/m^2]
  real :: hxferdb   !    ! heat xfer from seaspray drops to bcan_air this step [J/m^2]
  real :: wxferdb   !    ! vapor xfer from seaspray drops to bcan_air this step [kg_vap/m^2]
  real :: hxferse        ! heat xfer from sea surface to seaspray2 drops this step [J/m^2]
  real :: hxfereb   !    ! heat xfer from seaspray2 drops to bcan_air this step [J/m^2]
  real :: wxfereb   !    ! vapor xfer from seaspray2 drops to bcan_air this step [kg_vap/m^2]
  real :: sfc_rhovs      ! sat vapor density at sea surface temp [kg_vap/m^3]
  real :: spray_rhovs    ! sat vapor density at seaspray temp [kg_vap/m^3]
  real :: spray_rhovsp   ! sat vapor density derivative at seaspray temp [kg_vap/(K m^3)]
  real :: spray2_rhovs   ! sat vapor density at seaspray2 temp [kg_vap/m^3]
  real :: spray2_rhovsp  ! sat vapor density derivative at seaspray2 temp [kg_vap/(K m^3)]

  real :: can_rhov       ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair         ! Canopy air mass [kg/m^2]
  real :: hcapcan        ! Canopy air heat capacity [J/(m^2 K)]
  real :: bcan_rhov      ! Bcanopy air water vapor density [kg_vap/m^3]
  real :: bcanair        ! Bcanopy air mass [kg/m^2]
  real :: hcapbcan       ! Bcanopy air heat capacity [J/(m^2 K)]
  real :: hcapspray      ! Heat capacity of sea spray [J/(m^2 K)]
  real :: hcapspray2     ! Heat capacity of sea spray2 [J/(m^2 K)]

  real :: canairi        ! Inverse of canair
  real :: hcapcani       ! Inverse of hcapcan
  real :: bcanairi       ! Inverse of bcanair
  real :: hcapbcani      ! Inverse of hcapbcan
  real :: hcapsprayi     ! Inverse of hcapspray
  real :: hcapspray2i    ! Inverse of hcapspray2
  real :: zspray         ! height of top of spray layer [m]
  real :: zspray0        ! height of top of spray layer [m]
  real :: zsprayo2       ! middle of middle of spray layer [m]

  real    :: spray_massflux  ! Seaspray mass flux from sea [kg m^-2 s^-1]
  real    :: spray_mass      ! Seaspray steady-state mass in canopy [kg m^-2]
  real    :: spray_suodt     ! Seaspray steady-state vapor transfer coef [m s^-1]
  real    :: spray2_massflux ! Seaspray2 mass flux from sea [kg m^-2 s^-1]
  real    :: spray2_mass     ! Seaspray2 steady-state mass in canopy [kg m^-2]
  real    :: spray2_suodt    ! Seaspray2 steady-state vapor transfer coef [m s^-1]
  real    :: wt1, wt2        ! Interpolation weights [ ]
  integer :: ind             ! Interpolation index

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
  real(r8) :: h1, h2,     h4, h5, h6, h7, h8, h9, h10, h11
  real(r8) :: y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12

  real(r8) :: aa4(4,4), xx4(4), yy4(4)         ! 4x4 matrix equation terms
  real(r8) :: aa12(12,12), xx12(12), yy12(12)  ! 12x12 matrix equation terms

  logical :: sing

  ! Local variables used to evaluate wind speed at 10 m height

  real :: z10, cantheta, canthetav, airthetav
  real :: wind_z10, theta_z10, rrv_z10
  real :: press_bcan, exner_bcan, bcantheta, a2fm, a2fh

  real :: rs, rsi, ra, rb, rbi, richnum, vels0

  ! Seaspray is represented in the "sea canopy" in an analogous manner to water
  ! on the surface of vegetation in the land canopy.  HOWEVER, SINCE SEASPRAY
  ! DROPLETS ARE INJECTED TO HEIGHTS WELL ABOVE THE ROUGHNESS HEIGHT, A SECONDARY
  ! "CANOPY" LAYER, WHICH IS CALLED THE 'BCANOPY', IS INTRODUCED TO RECEIVE THE
  ! SEASPRAY DROPLETS.  TEMPERATURE AND VAPOR MIXING RATIO ARE PROGNOSED IN THE
  ! BCANOPY.  Heat and vapor are exchanged between seaspray and bcanopy air within
  ! an implicit solver that also incorporates heat and vapor exchange between the
  ! sea surface and canopy air, between canopy and bcanopy air, and between
  ! bcanopy air and the free atmosphere.  Three properties of seaspray are required:
  ! (1) the steady-state mass of seaspray droplets [kg m^-2] in the canopy air, (2)
  ! the steady-state vapor transfer coefficient [m s^-1] between seaspray and bcanopy
  ! air, and (3) the mass flux [kg m^-2 s^-1] of seaspray into the bcanopy air from
  ! the sea surface.  (By multiplying the mass flux by the temperature difference
  ! between new and returning sea spray droplets and by the specific heat of liquid
  ! water, we get the net energy flux from the sea surface to seaspray.)

  ! These three properties are tabulated below as a function of wind speed at 10 m
  ! height.  The tabulated values were generated in a separate program (ssgf.f90),
  ! which processed information presented in Figs 3i and 4a of Barr et al. (2023).
  ! ssgf.f90 numerically integrates over the (seastate-based) seaspray droplet size
  ! spectra in Fig. 4a to obtain all three properties.

  ! Vapor transfer coefficient suodt between sea spray and bcanopy air corresponds,
  ! after multiplication by the model timestep, to the quantity (su) in OLAM
  ! microphysics that governs vapor flux between liquid droplets and air.
  ! suodt is computed from the same formula, except for the following differences:
  !
  ! (1) su in microphysics is dimensionless, while (suodt * dt_sea) has dimensions
  !     of [m] representing a vertical integration over "bcanopy" depth.
  !
  ! (2) su in microphysics represents an integral over individual droplet size
  !     categories (cloud, drizzle, rain) assuming a gamma size distribution for
  !     each, whereas suodt is obtained by numerically integrating over the empirical
  !     seaspray droplet spectra in Fig. 4a.

  ! 10m height wind speed
  real ::   ssgf_windz10 (5) = [ 15.,    20.,    30.,    40.,    50. ] ! [m s^-1]

  ! 20-4000 MICRONS (band 1)
  real ::   ssgf_massflux(5) = [  0., .00007, .00135, .00702, .02769 ] ! [kg m^-2 s^-1]
  real ::   ssgf_mass    (5) = [  0., .00110, .01270, .04760, .12738 ] ! [kg m^-2]
  real ::   ssgf_suodt   (5) = [  0.,  0.059,  0.389,  1.086,  1.917 ] ! [m s^-1]

  ! 20-100 MICRONS (band 2)
  real ::  Assgf_massflux(5) = [  0., .00004, .00023, .00060, .00094 ] ! [kg m^-2 s^-1]
  real ::  Assgf_mass    (5) = [  0., .00089, .00665, .01918, .03373 ] ! [kg m^-2]
  real ::  Assgf_suodt   (5) = [  0.,  0.057,  0.339,  0.904,  1.481 ] ! [m s^-1]

  ! 20-200 MICRONS (band 3)
  real ::  Bssgf_massflux(5) = [  0., .00007, .00069, .00189, .00332 ] ! [kg m^-2 s^-1]
  real ::  Bssgf_mass    (5) = [  0., .00108, .01021, .03078, .05793 ] ! [kg m^-2]
  real ::  Bssgf_suodt   (5) = [  0.,  0.059,  0.378,  1.030,  1.744 ] ! [m s^-1]

  ! significant wave height
  real ::  sig_wavehght  (5) = [ 4.0,    6.0,    8.0,   10.0,   11.0 ] ! [m]

  integer :: band

  band = 1
  if (nl%iseasprayflg == 2) band = 35

  ! If doing seaspray and wind speed is at least seaspray_vmin, compute wind speed at 10 m level
  ! (wind_z10).  (Note that SFLUXT and SFLUXR are not actually used in sfclyr_profile
  ! to compute wind_z10, so it does not matter that they may represent bcanopy-to-atm
  ! fluxes.)

  if (nl%iseasprayflg > 0 .and. vels > nl%seaspray_vmin) then
     z10 = 10.

     cantheta  = cantemp / canexner
     canthetav = cantheta             * (1.0 + eps_virt * canrrv)
     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * airrrv)

     call sfclyr_profile (vels, rhos, canexner, ustar, sfluxt, sfluxr, &
                          sfcg%dzt_bot(iwsfc), rough, wstar, &
                          cantheta, canthetav, canrrv, airthetav, &
                          z10, wind_z10, theta_z10, rrv_z10)
  else
     wind_z10 = 0.
  endif

  ! Evaluate surface saturation vapor density and mixing ratio of sea surface

  sfc_rhovs = fssat * rhovsl(seatc-273.15)

  ! rdi = ustar/5 is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  ! Canopy air quantities

  can_rhov  = canrrv * rhos
  canair    = rhos * can_depth * 0.01
  canairi   = 1. / canair
  hcapcan   = cp * canair
  hcapcani  = 1. / hcapcan

  ! Set up and solve a linear system of equations that use trapezoidal-implicit
  ! differencing.  The solution of the system consists of turbulent heat and
  ! water vapor fluxes between canopy air, the ocean surface, sea spray droplets
  ! (if present), and the free atmosphere, and the consequent changes to water
  ! and temperature of canopy air and temperature of sea spray.

  if (wind_z10 <= nl%seaspray_vmin) then ! Case with NO SEA SPRAY (two uncoupled equations)

     spray_active = .false.

     ! Set up and solve 4x4 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

     a5  = dt_sea * rdi       ! sfc sublayer vap xfer coef
     a6  = cp * rhos * a5     ! sfc sublayer heat xfer coef
     a9  = dt_sea * vkhsfc / sfcg%dzt_bot(iwsfc)
     a10 = cp * a9

     h4 = fcn * rhos * canairi  ! = fcn / can_depth
     h7 = fcn * hcapcani
     h8 = fcn * canairi

     y2  = sfc_rhovs - can_rhov
     y5  = seatc     - cantemp
     y9  = canrrv    - airrrv
     y10 = cantemp   - canexner * airtheta

     aa4(1,1) = 1._r8 + a5 * h4
     aa4(1,2) = 0._r8
     aa4(1,3) =       - a5 * h4
     aa4(1,4) = 0._r8
     yy4(1)   =         a5 * y2   ! WSC row

     aa4(2,1) = 0._r8
     aa4(2,2) = 1._r8 + a6 * h7
     aa4(2,3) = 0._r8
     aa4(2,4) =       - a6 * h7
     yy4(2)   =         a6 * y5   ! HSC row

     aa4(3,1) =       - a9 * h8
     aa4(3,2) = 0._r8
     aa4(3,3) = 1._r8 + a9 * h8
     aa4(3,4) = 0._r8
     yy4(3)   =         a9 * y9   ! WCA row

     aa4(4,1) = 0._r8
     aa4(4,2) =       - a10 * h7
     aa4(4,3) = 0._r8
     aa4(4,4) = 1._r8 + a10 * h7
     yy4(4)   =         a10 * y10 ! HCA row

     call matrix8_NxN(4,aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'sea3',4,aa4,yy4,glatw,glonw)

     wxfersc = xx4(1)
     hxfersc = xx4(2)
     wxferca = xx4(3)
     hxferca = xx4(4)

     cantemp = cantemp + (hxfersc - hxferca) * hcapcani
     canrrv  = canrrv  + (wxfersc - wxferca) * canairi

     sfluxt = hxferca / dt_sea
     sfluxr = wxferca / dt_sea

     hfluxsea = hxfersc / (cp * dt_sea)

     if (nl%iseasprayflg > 0) then
        spraytemp  = seatc
        spray2temp = seatc
     endif

  else                                        ! Case WITH SEA SPRAY (12x12 matrix equation)

     ! For a given value of 10 m wind speed (windz10), interpolate from ssgf tables
     ! to get seaspray mass flux and steady-state values for mass and vapor transfer
     ! coefficient of seaspray in the canopy.

     ind = max(1,min(4,int(0.1 * wind_z10)))
     wt2 = (wind_z10 - ssgf_windz10(ind)) / (ssgf_windz10(ind+1) - ssgf_windz10(ind))
     wt1 = 1. - wt2

     if (band == 1) then
        ! 20-4000 MICRONS
        spray_massflux = wt1 * ssgf_massflux(ind) + wt2 * ssgf_massflux(ind+1)
        spray_mass     = wt1 * ssgf_mass    (ind) + wt2 * ssgf_mass    (ind+1)
        spray_suodt    = wt1 * ssgf_suodt   (ind) + wt2 * ssgf_suodt   (ind+1)

        spray2_massflux = 0.
        spray2_mass     = spray_mass
        spray2_suodt    = 0.

     elseif (band == 11) then  ! Don't use with iseasprayflg=1; spray2 is not saved
        ! 20-4000 MICRONS
        spray2_massflux = wt1 * ssgf_massflux(ind) + wt2 * ssgf_massflux(ind+1)
        spray2_mass     = wt1 * ssgf_mass    (ind) + wt2 * ssgf_mass    (ind+1)
        spray2_suodt    = wt1 * ssgf_suodt   (ind) + wt2 * ssgf_suodt   (ind+1)

        spray_massflux = 0.
        spray_mass     = spray2_mass
        spray_suodt    = 0.

     elseif (band == 24) then
        ! 20-100 MICRONS
        spray_massflux = wt1 * Assgf_massflux(ind) + wt2 * Assgf_massflux(ind+1)
        spray_mass     = wt1 * Assgf_mass    (ind) + wt2 * Assgf_mass    (ind+1)
        spray_suodt    = wt1 * Assgf_suodt   (ind) + wt2 * Assgf_suodt   (ind+1)

        ! 100-4000 MICRONS
        spray2_massflux = wt1 * (ssgf_massflux(ind)   - Assgf_massflux(ind  )) &
                        + wt2 * (ssgf_massflux(ind+1) - Assgf_massflux(ind+1))
        spray2_mass     = wt1 * (ssgf_mass    (ind)   - Assgf_mass    (ind  )) &
                        + wt2 * (ssgf_mass    (ind+1) - Assgf_mass    (ind+1))
        spray2_suodt    = wt1 * (ssgf_suodt   (ind)   - Assgf_suodt   (ind  )) &
                        + wt2 * (ssgf_suodt   (ind+1) - Assgf_suodt   (ind+1))

     elseif (band == 35) then
        ! 20-200 MICRONS
        spray_massflux = wt1 * Bssgf_massflux(ind) + wt2 * Bssgf_massflux(ind+1)
        spray_mass     = wt1 * Bssgf_mass    (ind) + wt2 * Bssgf_mass    (ind+1)
        spray_suodt    = wt1 * Bssgf_suodt   (ind) + wt2 * Bssgf_suodt   (ind+1)

        ! 200-4000 MICRONS
        spray2_massflux = wt1 * (ssgf_massflux(ind)   - Bssgf_massflux(ind  )) &
                        + wt2 * (ssgf_massflux(ind+1) - Bssgf_massflux(ind+1))
        spray2_mass     = wt1 * (ssgf_mass    (ind)   - Bssgf_mass    (ind  )) &
                        + wt2 * (ssgf_mass    (ind+1) - Bssgf_mass    (ind+1))
        spray2_suodt    = wt1 * (ssgf_suodt   (ind)   - Bssgf_suodt   (ind  )) &
                        + wt2 * (ssgf_suodt   (ind+1) - Bssgf_suodt   (ind+1))

     else
        stop 'ssgf - band value not valid'
     endif

     zspray0 = wt1 * sig_wavehght(ind) + wt2 * sig_wavehght(ind+1)
     zspray  = zspray0

     if (nl%umwmflg == 1 .and. nl%use_umwm_swh == 1) then
        if (umwm%iactive(isea)) then
           zspray      = max(3.0, swh(isea))
           spray_mass  = spray_mass  * (zspray / zspray0)
           spray2_mass = spray2_mass * (zspray / zspray0)
        endif
     endif

     zsprayo2 = 0.5 * zspray
     zsprayo2 = min(zsprayo2, .9*sfcg%dzt_bot(iwsfc)) ! limit middle of spray layer to no higher than
                                                      ! middle of lowest model layer for now

     press_bcan = prss - grav * zsprayo2 * rhos ! hydrostatic eqn.
     exner_bcan = (press_bcan * p00i) ** rocp

     vels0 = max(ubmin, sqrt(vels*vels + wstar*wstar))

     ! Surface layer bulk Richardson number computed from surface canopy
     richnum = 2.0 * grav * sfcg%dzt_bot(iwsfc) * (airthetav - canthetav)  &
             / ( (airthetav + canthetav) * vels0 * vels0 )

     ! Louis drag coefficients at middle of spray layer
     call a2fmfh(richnum, zsprayo2, rough, a2fm, a2fh)

     ! Convert heat/tracer drag coefficient to conductivity from can to bcan
     rsi = ustar * a2fh / sqrt(a2fm)

     rs = 1. / rsi                            ! can to bcan resistance
     ra = rhos * sfcg%dzt_bot(iwsfc) / vkhsfc ! can to atm resistance

     ! subtract the can-to-bcan resistance from the can-to-atm resistance
     ! to get the bcan-to-atm resistance

     rb  = ra - rs     ! bcan to atm resistance
     rbi = 1. / rb     ! bcan to atm conductivity

     if (spray_active) then

        bcantheta  = bcantemp / exner_bcan

     else

        spray_active = .true.

        ! Initialize bcan temp/vapor to give the same initial flux as from canopy
        bcantheta = airtheta + rb/ra * (cantheta - airtheta)
        bcanrrv   = airrrv   + rb/ra * (canrrv   - airrrv)

        bcantemp  = bcantheta * exner_bcan

     endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Test: method to turn off spray effect; must keep nonzero spray_mass
!  spray_massflux = 0.
!  spray_suodt = 0.

!  spray2_massflux = 0.
!  spray2_suodt = 0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     bcan_rhov = bcanrrv * rhos
     bcanair   = rhos * zspray
     bcanairi  = 1. / canair
     hcapbcan  = cp * canair
     hcapbcani = 1. / hcapcan

     spray_rhovs   = fssat * rhovsl(spraytemp        - 273.15)
     spray_rhovsp  = fssat * rhovsl(spraytemp  + 1.0 - 273.15) - spray_rhovs
     spray2_rhovs  = fssat * rhovsl(spray2temp       - 273.15)
     spray2_rhovsp = fssat * rhovsl(spray2temp + 1.0 - 273.15) - spray2_rhovs

     hcapspray   = spray_mass * cliq
     hcapsprayi  = 1. / hcapspray
     hcapspray2  = spray2_mass * cliq
     hcapspray2i = 1. / hcapspray2

     a1  = spray_suodt * dt_sea
     a2  = cp * rhos * a1
     a3  = spray2_suodt * dt_sea
     a4  = cp * rhos * a3
     a5  = dt_sea * rdi       ! sfc sublayer vap xfer coef
     a6  = cp * rhos * a5     ! sfc sublayer heat xfer coef
     a7  = spray_massflux * cliq * dt_sea
     a8  = spray2_massflux * cliq * dt_sea
     a9  = dt_sea * rbi * rhos  ! bcan to atm humidity mixing ratio xfer coef
     a10 = cp * a9              ! bcan to atm heat xfer coef
     a11 = dt_sea * rsi * rhos  ! can to bcan vap xfer coef
     a12 = cp * a11              ! can to bcan heat xfer coef

     h1  = fcn * hcapsprayi  * spray_rhovsp
     h2  = fcn * hcapspray2i * spray2_rhovsp
     h4  = fcn * rhos * canairi  ! = fcn / can_depth
     h5  = fcn * hcapsprayi
     h6  = fcn * hcapspray2i
     h7  = fcn * hcapcani
     h8  = fcn * canairi
     h9  = fcn * rhos * bcanairi ! = fcn / bcan_depth
     h10 = fcn * hcapbcani
     h11 = fcn * bcanairi

     y1  = spray_rhovs  - bcan_rhov
     y2  = sfc_rhovs    - can_rhov
     y4  = spraytemp    - bcantemp
     y5  = seatc        - cantemp
     y3  = seatc        - spraytemp
     y6  = spray2_rhovs - bcan_rhov
     y8  = spray2temp   - bcantemp
     y7  = seatc        - spray2temp
     y9  = bcanrrv      - airrrv
     y10 = bcantemp     - canexner * airtheta
     y11 = canrrv       - bcanrrv
     y12 = cantemp      - bcantemp

     ! Set up and solve 12x12 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

     aa12(1,1)  = 1._r8 + a1 * (h1 * alvl + h9)
     aa12(1,2)  = 0._r8
     aa12(1,3)  =         a1 * h1
     aa12(1,4)  = 0._r8
     aa12(1,5)  =       - a1 * h1
     aa12(1,6)  =         a1 * h9
     aa12(1,7)  = 0._r8
     aa12(1,8)  = 0._r8
     aa12(1,9)  =       - a1 * h9
     aa12(1,10) = 0._r8
     aa12(1,11) =         a1 * h9
     aa12(1,12) = 0._r8
     yy12(1)    =         a1 * y1  ! WDB row

     aa12(2,1)  = 0._r8
     aa12(2,2)  = 1._r8 + a5 * h4
     aa12(2,3)  = 0._r8
     aa12(2,4)  = 0._r8
     aa12(2,5)  = 0._r8
     aa12(2,6)  = 0._r8
     aa12(2,7)  = 0._r8
     aa12(2,8)  = 0._r8
     aa12(2,9)  = 0._r8
     aa12(2,10) = 0._r8
     aa12(2,11) =       - a5 * h4
     aa12(2,12) = 0._r8
     yy12(2)    =         a5 * y2  ! WSC row

     aa12(3,1)  =         a2 * h5 * alvl
     aa12(3,2)  = 0._r8
     aa12(3,3)  = 1._r8 + a2 * (h5 + h10)
     aa12(3,4)  = 0._r8
     aa12(3,5)  =       - a2 * h5
     aa12(3,6)  = 0._r8
     aa12(3,7)  =         a2 * h10
     aa12(3,8)  = 0._r8
     aa12(3,9)  = 0._r8
     aa12(3,10) =       - a2 * h10
     aa12(3,11) = 0._r8
     aa12(3,12) =         a2 * h10
     yy12(3)    =         a2 * y4  ! HDB row

     aa12(4,1)  = 0._r8
     aa12(4,2)  = 0._r8
     aa12(4,3)  = 0._r8
     aa12(4,4)  = 1._r8 + a6 * h7
     aa12(4,5)  = 0._r8
     aa12(4,6)  = 0._r8
     aa12(4,7)  = 0._r8
     aa12(4,8)  = 0._r8
     aa12(4,9)  = 0._r8
     aa12(4,10) = 0._r8
     aa12(4,11) = 0._r8
     aa12(4,12) =       - a6 * h7
     yy12(4)    =         a6 * y5  ! HSC row

     aa12(5,1)  =       - a7 * h5 * alvl
     aa12(5,2)  = 0._r8
     aa12(5,3)  =       - a7 * h5
     aa12(5,4)  = 0._r8
     aa12(5,5)  = 1._r8 + a7 * h5
     aa12(5,6)  = 0._r8
     aa12(5,7)  = 0._r8
     aa12(5,8)  = 0._r8
     aa12(5,9)  = 0._r8
     aa12(5,10) = 0._r8
     aa12(5,11) = 0._r8
     aa12(5,12) = 0._r8
     yy12(5)    =         a7 * y3  ! HSD row

     aa12(6,1)  =         a3 * h9
     aa12(6,2)  = 0._r8
     aa12(6,3)  = 0._r8
     aa12(6,4)  = 0._r8
     aa12(6,5)  = 0._r8
     aa12(6,6)  = 1._r8 + a3 * (h2 * alvl + h9)
     aa12(6,7)  =         a3 * h2
     aa12(6,8)  =       - a3 * h2
     aa12(6,9)  =       - a3 * h9
     aa12(6,10) = 0._r8
     aa12(6,11) =         a3 * h9
     aa12(6,12) = 0._r8
     yy12(6)    =         a3 * y6  ! WEB row

     aa12(7,1)  = 0._r8
     aa12(7,2)  = 0._r8
     aa12(7,3)  =         a4 * h10
     aa12(7,4)  = 0._r8
     aa12(7,5)  = 0._r8
     aa12(7,6)  =         a4 * h6 * alvl
     aa12(7,7)  = 1._r8 + a4 * (h6 + h10)
     aa12(7,8)  =       - a4 * h6
     aa12(7,9)  = 0._r8
     aa12(7,10) =       - a4 * h10
     aa12(7,11) = 0._r8
     aa12(7,12) =         a4 * h10
     yy12(7)    =         a4 * y8  ! HEB row

     aa12(8,1)  = 0._r8
     aa12(8,2)  = 0._r8
     aa12(8,3)  = 0._r8
     aa12(8,4)  = 0._r8
     aa12(8,5)  = 0._r8
     aa12(8,6)  =       - a8 * h6 * alvl
     aa12(8,7)  =       - a8 * h6
     aa12(8,8)  = 1._r8 + a8 * h6
     aa12(8,9)  = 0._r8
     aa12(8,10) = 0._r8
     aa12(8,11) = 0._r8
     aa12(8,12) = 0._r8
     yy12(8)    =         a8 * y7  ! HSE row

     aa12(9,1)  =       - a9 * h11
     aa12(9,2)  = 0._r8
     aa12(9,3)  = 0._r8
     aa12(9,4)  = 0._r8
     aa12(9,5)  = 0._r8
     aa12(9,6)  =       - a9 * h11
     aa12(9,7)  = 0._r8
     aa12(9,8)  = 0._r8
     aa12(9,9)  = 1._r8 + a9 * h11
     aa12(9,10) = 0._r8
     aa12(9,11) =       - a9 * h11
     aa12(9,12) = 0._r8
     yy12(9)    =         a9 * y9  ! WBA row

     aa12(10,1)  = 0._r8
     aa12(10,2)  = 0._r8
     aa12(10,3)  =       - a10 * h10
     aa12(10,4)  = 0._r8
     aa12(10,5)  = 0._r8
     aa12(10,6)  = 0._r8
     aa12(10,7)  =       - a10 * h10
     aa12(10,8)  = 0._r8
     aa12(10,9)  = 0._r8
     aa12(10,10) = 1._r8 + a10 * h10
     aa12(10,11) = 0._r8
     aa12(10,12) =       - a10 * h10
     yy12(10)    =         a10 * y10  ! HBA row

     aa12(11,1)  =         a11 * h11
     aa12(11,2)  =       - a11 * h8
     aa12(11,3)  = 0._r8
     aa12(11,4)  = 0._r8
     aa12(11,5)  = 0._r8
     aa12(11,6)  =         a11 * h11
     aa12(11,7)  = 0._r8
     aa12(11,8)  = 0._r8
     aa12(11,9)  =       - a11 * h11
     aa12(11,10) = 0._r8
     aa12(11,11) = 1._r8 + a11 * (h8 + h11)
     aa12(11,12) = 0._r8
     yy12(11)    =         a11 * y11    ! WCB row

     aa12(12,1)  = 0._r8
     aa12(12,2)  = 0._r8
     aa12(12,3)  =         a12 * h10
     aa12(12,4)  =       - a12 * h7
     aa12(12,5)  = 0._r8
     aa12(12,6)  = 0._r8
     aa12(12,7)  =         a12 * h10
     aa12(12,8)  = 0._r8
     aa12(12,9)  = 0._r8
     aa12(12,10) =       - a12 * h10
     aa12(12,11) = 0._r8
     aa12(12,12) = 1._r8 + a12 * (h7 + h10)
     yy12(12)    =         a12 * y12   ! HCB row

     call matrix8_NxN(12,aa12,yy12,xx12,sing); if (sing) call sing_print(iwsfc,'sea4',12,aa12,yy12,glatw,glonw)

     wxferdb = xx12(1)
     wxfersc = xx12(2)
     hxferdb = xx12(3)
     hxfersc = xx12(4)
     hxfersd = xx12(5)
     wxfereb = xx12(6)
     hxfereb = xx12(7)
     hxferse = xx12(8)
     wxferba = xx12(9)
     hxferba = xx12(10)
     wxfercb = xx12(11)
     hxfercb = xx12(12)

     cantemp    = cantemp    + (hxfersc - hxfercb) * hcapcani
     canrrv     = canrrv     + (wxfersc - wxfercb) * canairi

     bcantemp   = bcantemp   + (hxferdb + hxfereb + hxfercb - hxferba) * hcapbcani
     bcanrrv    = bcanrrv    + (wxferdb + wxfereb + wxfercb - wxferba) * bcanairi

     spraytemp  = spraytemp  + (hxfersd - hxferdb - wxferdb * alvl)    * hcapsprayi
     spray2temp = spray2temp + (hxferse - hxfereb - wxfereb * alvl)    * hcapspray2i

     sfluxt = hxferba / dt_sea
     sfluxr = wxferba / dt_sea

     hfluxsea = (hxfersc + hxferdb) / (cp * dt_sea)

  endif

end subroutine seacell_2
