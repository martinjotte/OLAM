subroutine seacells(isea, timefac_sst, timefac_seaice)

  use sea_coms,    only: nzi
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, omsea
  use misc_coms,   only: iparallel
  use consts_coms, only: grav, p00i, rocp, eps_virt, cpi, erad, eradi, cliq1000
  use mem_para,    only: myrank
  use oname_coms,  only: nl
  use pom2k1d,     only: pom, rhoref, pom_column

  implicit none

  integer, intent(in) :: isea
  real,    intent(in) :: timefac_sst    ! fraction of elapsed time from past to future SST obs
  real,    intent(in) :: timefac_seaice ! fraction of elapsed time from past to future SEA ICE obs

  ! Local variables

  integer :: iwsfc

  real :: canexneri, cantheta, canthetav
  real :: airthetav, ufree
  real :: usti, zw, zn1, zn2
  real :: raxis, windu, windv, cdtop, wusurf, wvsurf, wtsurf, wssurf, swrad, hfluxsea
  real :: sea_spray1_temp, sea_spray2_temp

  real, parameter :: z0fac = 0.011 / grav  ! factor for Charnok roughness height
  real, parameter :: ozo   = 1.59e-5       ! base roughness height in HWRF
  real, parameter :: onethird = 1./3.

  iwsfc = isea + omsea

  ! Skip this cell if running in parallel and cell rank is not MYRANK

  if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) return

  ! Zero runoff for sea cells

  sfcg%runoff(iwsfc) = 0.

  ! Update seaice fraction

  sea%seaicec(isea) = sea%seaicep(isea) &
                    + timefac_seaice * (sea%seaicef(isea) - sea%seaicep(isea))

  ! Update SEATC if isea cell is not pom_active or is is swm_active

  if (.not. sea%pom_active(isea) .or. sfcg%swm_active(iwsfc)) then

     sea%seatc(isea) = sea%seatp(isea) &
                     + timefac_sst * (sea%seatf(isea) - sea%seatp(isea))
  endif

  ! Evaluate sea roughness height

  ! Charnok (1955):
  ! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)

  ! Davis et al. (2008) originally used in HWRF
  ! rough = 10. * exp(-10. / ustar ** .333333)
  ! rough = max(.125e-6, min(2.85e-3,rough))

  ! 2012 HWRF scheme; interpolates between the Charnok scheme at low wind
  ! and the Davis et al. curve fit at high wind speeds

  usti  = 1.0 / sea%sea_ustar(isea)
  zw    = min( (sea%sea_ustar(isea)/1.06)**onethird, 1.0 )
  zn1   = z0fac * sea%sea_ustar(isea) * sea%sea_ustar(isea) + ozo
  zn2   = 10. * exp(-9.5 * usti**onethird) + 1.65e-6 * usti

  sea%sea_rough(isea) = min( (1.0-zw) * zn1 + zw * zn2, 2.85e-3 )

  ! Evaluate surface layer exchange coefficients vkmsfc and vkhsfc for open water areas

  airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))
  canexneri = 1. / sfcg%canexner(iwsfc)

  cantheta  = sea%sea_cantemp(isea) * canexneri
  canthetav = cantheta * (1.0 + eps_virt * sea%sea_canrrv(isea))

  ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sea%sea_wthv(isea),0.0) / airthetav) ** onethird

  call stars(sfcg%dzt_bot  (iwsfc), &
             sea%sea_rough  (isea), &
             sfcg%vels     (iwsfc), &
             sfcg%rhos     (iwsfc), &
             ufree                , &
             airthetav            , &
             canthetav            , &
             sea%sea_vkmsfc (isea), &
             sea%sea_vkhsfc (isea), &
             sea%sea_ustar  (isea), &
             sea%sea_ggaer  (isea)  )

! When we have CO2:
!  sea%sea_sfluxc(isea) = sfcg%rhos(iwsfc) * sea%sea_ggaer(isea) &
!                       * (sea%sea_co2(isea) - air_co2)

  if (allocated(sea%spraytemp )) sea_spray1_temp = sea%spraytemp (isea)
  if (allocated(sea%spray2temp)) sea_spray2_temp = sea%spray2temp(isea)

  ! Update SEA fields

  if (nl%iseasprayflg /= 1) then

     call seacell_2(isea, iwsfc,             &
                    sfcg%rhos      (iwsfc),  &
                    sea%sea_ustar   (isea),  &
                    sea%sea_vkhsfc  (isea),  &
                    sfcg%can_depth (iwsfc),  &
                    sea%seatc       (isea),  &
                    sfcg%vels      (iwsfc),  &
                    sfcg%prss      (iwsfc),  &
                    sea%sea_wthv    (isea),  &
                    sfcg%glatw     (iwsfc),  &
                    sfcg%glonw     (iwsfc),  &
                    sfcg%airtheta  (iwsfc),  &
                    sfcg%airrrv    (iwsfc),  &
                    sfcg%canexner  (iwsfc),  &
                    sea%sea_cantemp (isea),  &
                    sea%sea_canrrv  (isea),  &
                    sea_spray1_temp,         &
                    sea_spray2_temp,         &
                    sea%sea_sfluxt  (isea),  &
                    sea%sea_sfluxr  (isea),  &
                    sea%sea_sfc_srrv(isea),  &
                    sea%sea_rough   (isea),  &
                    hfluxsea                 )

  else

     call seacell_1(isea, iwsfc,             &
                    sfcg%rhos      (iwsfc),  &
                    sea%sea_ustar   (isea),  &
                    sea%sea_vkhsfc  (isea),  &
                    sfcg%can_depth (iwsfc),  &
                    sea%seatc       (isea),  &
                    sfcg%vels      (iwsfc),  &
                    sfcg%prss      (iwsfc),  &
                    sea%sea_wthv    (isea),  &
                    sfcg%glatw     (iwsfc),  &
                    sfcg%glonw     (iwsfc),  &
                    sfcg%airtheta  (iwsfc),  &
                    sfcg%airrrv    (iwsfc),  &
                    sfcg%canexner  (iwsfc),  &
                    sea%sea_cantemp (isea),  &
                    sea%sea_canrrv  (isea),  &
                    sea_spray1_temp,         &
                    sea%sea_sfluxt  (isea),  &
                    sea%sea_sfluxr  (isea),  &
                    sea%sea_sfc_srrv(isea),  &
                    sea%sea_rough   (isea),  &
                    hfluxsea                 )

  endif

  ! New calculation of wthv for sfluxt units [W m^-2]

  sea%sea_wthv(isea) = ( sea%sea_sfluxt(isea) * cpi * canexneri * (1.0 + eps_virt * sfcg%airrrv(iwsfc)) &
                     + sea%sea_sfluxr(isea) * eps_virt * sfcg%airtheta(iwsfc) ) / sfcg%rhos(iwsfc)

  ! Update POM1D vertical column variables if this sea cell is pom_active and is not swm_active

  if (sea%pom_active(isea) .and. .not. sfcg%swm_active(iwsfc)) then

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

     call pom_column(isea, pom%kba(isea), wusurf, wvsurf, wtsurf, wssurf, swrad)

  endif


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
     sea%surface_srrv (isea) = sea%sea_sfc_srrv(isea)

  else

     ! If sea ice is present in this cell, compute turbulent fluxes over ice

     cantheta  = sea%ice_cantemp(isea) * canexneri
     canthetav = cantheta * (1.0 + eps_virt * sea%ice_canrrv(isea))

     ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sea%ice_wthv(isea),0.0) / airthetav) ** onethird

     call stars(sfcg%dzt_bot  (iwsfc), &
                sea%ice_rough  (isea), &
                sfcg%vels     (iwsfc), &
                sfcg%rhos     (iwsfc), &
                ufree                , &
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
                 sea%sea_sfluxt         (isea), &
                 sea%sea_sfluxr         (isea), &
                 sea%ice_sfc_srrv       (isea)  )

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

     sea%surface_srrv(isea) = (1.0 - sea%seaicec(isea)) * sea%sea_sfc_srrv(isea) + &
                                     sea%seaicec(isea)  * sea%ice_sfc_srrv(isea)

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

end subroutine seacells

!===============================================================================

subroutine seacell_1( isea, iwsfc, rhos, ustar, vkhsfc, can_depth, seatc, &
                      vels, prss, wthv, glatw, glonw, airtheta, airrrv,  &
                      canexner, cantemp, canrrv, spraytemp, sfluxt, sfluxr, &
                      surface_srrv, rough, hfluxsea )

  use mem_sfcg,    only: sfcg
  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav, alvl, cliq, r8, p00i, rocp, eps_virt
  use therm_lib,   only: rhovsl
  use matrix,      only: matrix8_NxN
  use leaf4_canopy,only: sing_print
  use oname_coms,  only: nl

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
  integer, intent(in)    :: iwsfc        ! current sfcg cell index
  real,    intent(in)    :: rhos         ! air density [kg_dryair/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: vkhsfc       ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: seatc        ! current sea temp (obs time) [K]
  real,    intent(in)    :: vels         ! wind speed at top of sfc layer [m/s]
  real,    intent(in)    :: prss         ! canopy air pressure [Pa]
  real,    intent(in)    :: wthv         ! surface buoyancy flux [K m/s]
  real,    intent(in)    :: glatw        ! Latitude of sea cell 'center' [deg]
  real,    intent(in)    :: glonw        ! Longitude of sea cell 'center' [deg]
  real,    intent(in)    :: airtheta     ! atm potential temp [K]
  real,    intent(in)    :: airrrv       ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: canexner     ! canopy Exner function []
  real,    intent(inout) :: cantemp      ! can_air temp [K]
  real,    intent(inout) :: canrrv       ! can_air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(inout) :: spraytemp    ! seaspray temperature [K]
  real,    intent(inout) :: sfluxt       ! can_air to atm heat flux [W m^-2]
  real,    intent(inout) :: sfluxr       ! can_air to atm vapor flux [kg_vap m^-2 s^-1]
  real,    intent(out)   :: surface_srrv ! sea surface sat mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: rough        ! sea cell roughess height [m]
  real,    intent(out)   :: hfluxsea     ! heat flux from sea surface to can_air [kg K/(m^2 s)]

  ! Local parameters

  real, parameter :: one3  = 1./ 3.
  real, parameter :: fcn = 0.75            ! Crank-Nicolson future time weight

  ! Local variables

  real :: rdi            ! sea surface to can_air conductance [m/s]
  real :: hxfersc        ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca        ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca        ! vapor xfer from can_air to atm this step [kg/m^2]
  real :: wxfersc        ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]
  real :: hxfersd        ! heat xfer from sea surface to seaspray drops this step [J/m^2]
  real :: hxferdc        ! heat xfer from seaspray drops to can_air this step [J/m^2]
  real :: wxferdc        ! vapor xfer from seaspray drops to can_air this step [kg_vap/m^2]
  real :: sfc_rhovs      ! sat vapor density at sea surface temp [kg_vap/m^3]
  real :: spray_rhovs    ! sat vapor density at seaspray temp [kg_vap/m^3]
  real :: spray_rhovsp   ! sat vapor density derivative at seaspray temp [kg_vap/(K m^3)]
  real :: can_rhov       ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair         ! Canopy air mass [kg/m^2]
  real :: hcapcan        ! Canopy air heat capacity [J/(m^2 K)]
  real :: hcapspray      ! Heat capacity of sea spray [J/(m^2 K)]
  real :: canairi        ! Inverse of canair
  real :: hcapcani       ! Inverse of hcapcan
  real :: hcapsprayi     ! Inverse of hcapspray

  real    :: spray_massflux ! Seaspray mass flux from sea [kg m^-2 s^-1]
  real    :: spray_mass     ! Seaspray steady-state mass in canopy [kg m^-2]
  real    :: spray_suodt    ! Seaspray steady-state vapor transfer coef [m s^-1]
  real    :: wt1, wt2       ! Interpolation weights [ ]
  integer :: ind            ! Interpolation index

  real(r8) :: a1, a2, a5, a6, a7, a9, a10
  real(r8) :: h1, h4, h5, h7, h8
  real(r8) :: y1, y2, y4, y5, y3, y9, y10

  real(r8) :: aa4(4,4), xx4(4), yy4(4)  ! 4x4 matrix equation terms
  real(r8) :: aa7(7,7), xx7(7), yy7(7)  ! 7x7 matrix equation terms

  logical :: sing

  ! Local variables used to evaluate wind speed at 10 m height

  real :: z10, press_z10, exner_z10, cantheta, canthetav, airthetav
  real :: ufree, wind_z10, theta_z10, rrv_z10

  ! Seaspray is represented in the "sea canopy" in an analogous manner to water
  ! on the surface of vegetation in the land canopy.  Heat and vapor are exchanged
  ! between seaspray and canopy air within an implicit solver that also incorporates
  ! heat and vapor exchange between the sea surface and canopy air, and between
  ! canopy air and the free atmosphere.  Three properties of seaspray are required:
  ! (1) the steady-state mass of seaspray droplets [kg m^-2] in the canopy air, (2)
  ! the steady-state vapor transfer coefficient [m s^-1] between seaspray and canopy
  ! air, and (3) the mass flux [kg m^-2 s^-1] of seaspray into the canopy air from
  ! the sea surface.  (By multiplying the mass flux by the temperature difference
  ! between new and returning sea spray droplets and by the specific heat of liquid
  ! water, we get the net energy flux from the sea surface to seaspray.) 

  ! These three properties are tabulated below as a function of wind speed at 10 m
  ! height.  The tabulated values were generated in a separate program (ssgf.f90),
  ! which processed information presented in Figs 3i and 4a of Barr et al. (2023).
  ! ssgf.f90 numerically integrates over the (seastate-based) seaspray droplet size
  ! spectra in Fig. 4a to obtain all three properties.

  ! Vapor transfer coefficient suodt between sea spray and canopy air corresponds,
  ! after multiplication by the model timestep, to the quantity (su) in OLAM
  ! microphysics divided by that governs vapor flux between liquid droplets and air.
  ! suodt is computed from the same formula, except for the following differences:
  !
  ! (1) su in microphysics is dimensionless, while (suodt * dt_sea) has dimensions
  !     of [m] representing a vertical integration over "canopy" depth.
  !
  ! (2) su in microphysics represents an integral over individual droplet size
  !     categories (cloud, drizzle, rain) assuming a gamma size distribution for
  !     each, whereas suodt is obtained by numerically integrating over the empirical
  !     seaspray droplet spectra in Fig. 4a.

  ! 20-4000 MICRONS (band 1)
  real ::  ssgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  ssgf_massflux(5) = (/  0., .00007, .00135, .00702, .02769 /) ! [kg m^-2 s^-1]
  real ::  ssgf_mass    (5) = (/  0., .00110, .01270, .04760, .12738 /) ! [kg m^-2]
  real ::  ssgf_suodt   (5) = (/  0.,  0.059,  0.389,  1.086,  1.917 /) ! [m s^-1]

  ! 20-100 MICRONS (band 2)
  real ::  Assgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  Assgf_massflux(5) = (/  0., .00004, .00023, .00060, .00094 /) ! [kg m^-2 s^-1]
  real ::  Assgf_mass    (5) = (/  0., .00089, .00665, .01918, .03373 /) ! [kg m^-2]
  real ::  Assgf_suodt   (5) = (/  0.,  0.057,  0.339,  0.904,  1.481 /) ! [m s^-1]

  ! 20-200 MICRONS (band 3)
  real ::  Bssgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  Bssgf_massflux(5) = (/  0., .00007, .00069, .00189, .00332 /) ! [kg m^-2 s^-1]
  real ::  Bssgf_mass    (5) = (/  0., .00108, .01021, .03078, .05793 /) ! [kg m^-2]
  real ::  Bssgf_suodt   (5) = (/  0.,  0.059,  0.378,  1.030,  1.744 /) ! [m s^-1]

  ! 100-4000 MICRONS (band 4)

  ! 200-4000 MICRONS (band 5)

  integer, parameter :: band = 1

  ! If doing seaspray and wind speed is at least 15 m/s, compute wind speed at 10 m level
  ! (wind_z10) USING SFLUXT AND SFLUXR VALUES EVALUATED ON THE PREVIOUS TIMESTEP.
  ! Otherwise, set wind_z10 to zero as a flag to not compute seaspray flux.

  if (nl%iseasprayflg > 0 .and. vels > nl%seaspray_vmin) then
     z10 = 10.

     press_z10 = prss - grav * z10 * rhos ! hydrostatic eqn.
     exner_z10 = (press_z10 * p00i) ** rocp

     cantheta  = cantemp / canexner
     canthetav = cantheta             * (1.0 + eps_virt * canrrv)
     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * airrrv)

     ufree = (grav * sfcg%dzt_bot(iwsfc) * max(wthv,0.0) / airthetav) ** one3

     call sfclyr_profile (vels, rhos, canexner, ustar, sfluxt, sfluxr, &
                          sfcg%dzt_bot(iwsfc), rough, ufree, &
                          cantheta, canthetav, canrrv, airthetav, &
                          z10, wind_z10, theta_z10, rrv_z10)
  else
     wind_z10 = 0.
  endif

  ! Evaluate surface saturation vapor density and mixing ratio of sea surface

  sfc_rhovs    = rhovsl(seatc-273.15)
  surface_srrv = sfc_rhovs / rhos

  ! rdi = ustar/5 is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  ! Canopy air quantities

  can_rhov = canrrv * rhos
  canair   = rhos * can_depth
  canairi  = 1. / canair
  hcapcan  = cp * canair
  hcapcani = 1. / hcapcan

  ! Set up and solve a linear system of equations that use trapezoidal-implicit
  ! differencing.  The solution of the system consists of turbulent heat and
  ! water vapor fluxes between canopy air, the ocean surface, sea spray droplets
  ! (if present), and the free atmosphere, and the consequent changes to water
  ! and temperature of canopy air and temperature of sea spray.

  a5  = dt_sea * rdi   ! sfc sublayer vap xfer coef
  a6  = cp * rhos * a5 ! sfc sublayer heat xfer coef
  a9  = dt_sea * vkhsfc / sfcg%dzt_bot(iwsfc)
  a10 = cp * a9

  h4 = fcn * rhos * canairi  ! = fcn / can_depth
  h7 = fcn * hcapcani
  h8 = fcn * canairi

  y2  = sfc_rhovs - can_rhov
  y5  = seatc     - cantemp
  y9  = canrrv    - airrrv
  y10 = cantemp   - canexner * airtheta

  if (wind_z10 < nl%seaspray_vmin) then ! Case with NO SEA SPRAY (two uncoupled equations)

     ! Set up and solve 4x4 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

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

     call matrix8_NxN(4,aa4,yy4,xx4,sing); if (sing) call sing_print(iwsfc,'sea1',4,aa4,yy4,glatw,glonw)

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
        spraytemp = seatc
     endif

  else                                        ! Case WITH SEA SPRAY (7x7 matrix equation)

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

     elseif (band == 2) then
        ! 20-100 MICRONS
        spray_massflux = wt1 * Assgf_massflux(ind) + wt2 * Assgf_massflux(ind+1)
        spray_mass     = wt1 * Assgf_mass    (ind) + wt2 * Assgf_mass    (ind+1)
        spray_suodt    = wt1 * Assgf_suodt   (ind) + wt2 * Assgf_suodt   (ind+1)

     elseif (band == 3) then
        ! 20-200 MICRONS
        spray_massflux = wt1 * Bssgf_massflux(ind) + wt2 * Bssgf_massflux(ind+1)
        spray_mass     = wt1 * Bssgf_mass    (ind) + wt2 * Bssgf_mass    (ind+1)
        spray_suodt    = wt1 * Bssgf_suodt   (ind) + wt2 * Bssgf_suodt   (ind+1)

     elseif (band == 4) then
        ! 100-4000 MICRONS
        spray_massflux = wt1 * (ssgf_massflux(ind)   - Assgf_massflux(ind  )) &
                       + wt2 * (ssgf_massflux(ind+1) - Assgf_massflux(ind+1))
        spray_mass     = wt1 * (ssgf_mass    (ind)   - Assgf_mass    (ind  )) &
                       + wt2 * (ssgf_mass    (ind+1) - Assgf_mass    (ind+1))
        spray_suodt    = wt1 * (ssgf_suodt   (ind)   - Assgf_suodt   (ind  )) &
                       + wt2 * (ssgf_suodt   (ind+1) - Assgf_suodt   (ind+1))

     elseif (band == 5) then
        ! 200-4000 MICRONS
        spray_massflux = wt1 * (ssgf_massflux(ind)   - Bssgf_massflux(ind  )) &
                       + wt2 * (ssgf_massflux(ind+1) - Bssgf_massflux(ind+1))
        spray_mass     = wt1 * (ssgf_mass    (ind)   - Bssgf_mass    (ind  )) &
                       + wt2 * (ssgf_mass    (ind+1) - Bssgf_mass    (ind+1))
        spray_suodt    = wt1 * (ssgf_suodt   (ind)   - Bssgf_suodt   (ind  )) &
                       + wt2 * (ssgf_suodt   (ind+1) - Bssgf_suodt   (ind+1))

     else
        stop 'ssgf - band value not valid'
     endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Test: method to turn off spray effect; must keep nonzero spray_mass
!  spray_massflux = 0.
!  spray_suodt = 0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     a1 = spray_suodt * dt_sea
     a2 = cp * rhos * a1
     a7 = spray_massflux * cliq * dt_sea
     a9  = vkhsfc * dt_sea / sfcg%dzt_bot(iwsfc)
     a10 = cp * a9

     spray_rhovs  = rhovsl(spraytemp       - 273.15)
     spray_rhovsp = rhovsl(spraytemp + 1.0 - 273.15) - spray_rhovs

     hcapspray = spray_mass * cliq
     hcapsprayi = 1. / hcapspray

     h1 = fcn * hcapsprayi * spray_rhovsp
     h5 = fcn * hcapsprayi

     y1 = spray_rhovs - can_rhov
     y4 = spraytemp   - cantemp
     y3 = seatc       - spraytemp
     y9  = canrrv     - airrrv
     y10 = cantemp    - canexner * airtheta

     ! Set up and solve 7x7 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

     aa7(1,1) = 1._r8 + a1 * (h1 * alvl + h4)
     aa7(1,2) =         a1 * h4
     aa7(1,3) =         a1 * h1
     aa7(1,4) = 0._r8
     aa7(1,5) =       - a1 * h1
     aa7(1,6) =       - a1 * h4
     aa7(1,7) = 0._r8
     yy7(1)   =         a1 * y1  ! WDC row

     aa7(2,1) =         a5 * h4
     aa7(2,2) = 1._r8 + a5 * h4
     aa7(2,3) = 0._r8
     aa7(2,4) = 0._r8
     aa7(2,5) = 0._r8
     aa7(2,6) =       - a5 * h4
     aa7(2,7) = 0._r8
     yy7(2)   =         a5 * y2  ! WSC row

     aa7(3,1) =         a2 * h5 * alvl
     aa7(3,2) = 0._r8
     aa7(3,3) = 1._r8 + a2 * (h5 + h7)
     aa7(3,4) =         a2 * h7
     aa7(3,5) =       - a2 * h5
     aa7(3,6) = 0._r8
     aa7(3,7) =       - a2 * h7
     yy7(3)   =         a2 * y4  ! HDC row

     aa7(4,1) = 0._r8
     aa7(4,2) = 0._r8
     aa7(4,3) =         a6 * h7
     aa7(4,4) = 1._r8 + a6 * h7
     aa7(4,5) = 0._r8
     aa7(4,6) = 0._r8
     aa7(4,7) =       - a6 * h7
     yy7(4)   =         a6 * y5  ! HSC row

     aa7(5,1) =       - a7 * h5 * alvl
     aa7(5,2) = 0._r8
     aa7(5,3) =       - a7 * h5
     aa7(5,4) = 0._r8
     aa7(5,5) = 1._r8 + a7 * h5
     aa7(5,6) = 0._r8
     aa7(5,7) = 0._r8
     yy7(5)   =         a7 * y3  ! HSD row

     aa7(6,1) =       - a9 * h8
     aa7(6,2) =       - a9 * h8
     aa7(6,3) = 0._r8
     aa7(6,4) = 0._r8
     aa7(6,5) = 0._r8
     aa7(6,6) = 1._r8 + a9 * h8
     aa7(6,7) = 0._r8
     yy7(6)   =         a9 * y9  ! WCA row

     aa7(7,1) = 0._r8
     aa7(7,2) = 0._r8
     aa7(7,3) =       - a10 * h7
     aa7(7,4) =       - a10 * h7
     aa7(7,5) = 0._r8
     aa7(7,6) = 0._r8
     aa7(7,7) = 1._r8 + a10 * h7
     yy7(7)   =         a10 * y10  ! HCA row

     call matrix8_NxN(7,aa7,yy7,xx7,sing); if (sing) call sing_print(iwsfc,'sea2',7,aa7,yy7,glatw,glonw)

     wxferdc = xx7(1)
     wxfersc = xx7(2)
     hxferdc = xx7(3)
     hxfersc = xx7(4)
     hxfersd = xx7(5)
     wxferca = xx7(6)
     hxferca = xx7(7)

     cantemp   = cantemp   + (hxferdc + hxfersc - hxferca)        * hcapcani
     canrrv    = canrrv    + (wxferdc + wxfersc - wxferca)        * canairi
     spraytemp = spraytemp + (hxfersd - hxferdc - wxferdc * alvl) * hcapsprayi

     sfluxt = hxferca / dt_sea
     sfluxr = wxferca / dt_sea

     hfluxsea = (hxfersc + hxferdc) / (cp * dt_sea)

  endif

end subroutine seacell_1

!===============================================================================

subroutine seacell_2( isea, iwsfc, rhos, ustar, vkhsfc, can_depth, seatc, &
                      vels, prss, wthv, glatw, glonw, airtheta, airrrv, &
                      canexner, cantemp, canrrv, spraytemp, spray2temp, sfluxt, sfluxr, &
                      surface_srrv, rough, hfluxsea )

  use mem_sfcg,    only: sfcg
  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav, alvl, cliq, r8, p00i, rocp, eps_virt
  use therm_lib,   only: rhovsl
  use matrix,      only: matrix8_NxN
  use leaf4_canopy,only: sing_print
  use oname_coms,  only: nl

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
  integer, intent(in)    :: iwsfc        ! current sfcg cell index
  real,    intent(in)    :: rhos         ! air density [kg_dryair/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: vkhsfc       ! can_air to atm heat & vapor transfer coef [kg_dryair m^-1 s^-1]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: seatc        ! current sea temp (obs time) [K]
  real,    intent(in)    :: vels         ! wind speed at top of sfc layer [m/s]
  real,    intent(in)    :: prss         ! canopy air pressure [Pa]
  real,    intent(in)    :: wthv         ! surface buoyancy flux [K m/s]
  real,    intent(in)    :: glatw        ! Latitude of land cell 'center' [deg]
  real,    intent(in)    :: glonw        ! Longitude of land cell 'center' [deg]
  real,    intent(in)    :: airtheta     ! atm potential temp [K]
  real,    intent(in)    :: airrrv       ! atm vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: canexner     ! canopy Exner function []
  real,    intent(inout) :: cantemp      ! can_air temp [K]
  real,    intent(inout) :: canrrv       ! can_air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(inout) :: spraytemp    ! seaspray temperature [K]
  real,    intent(inout) :: spray2temp   ! seaspray2 temperature [K]
  real,    intent(inout) :: sfluxt       ! can_air to atm sens heat flux [W m^-2]
  real,    intent(inout) :: sfluxr       ! can_air to atm vap flux [kg_vap m^-2 s^-1]
  real,    intent(out)   :: surface_srrv ! sea surface sat mixing ratio [kg_vap/kg_dryair]
  real,    intent(in)    :: rough        ! sea cell roughess height [m]
  real,    intent(out)   :: hfluxsea     ! heat flux from sea to can_air [kg K/(m^2 s)]

  ! Local parameters

  real, parameter :: one3  = 1./ 3.
  real, parameter :: fcn = 0.75            ! Crank-Nicolson future time weight

  ! Local variables

  real :: rdi            ! sea surface to can_air conductance [m/s]
  real :: hxfersc        ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca        ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxferca        ! vapor xfer from can_air to atm this step [kg/m^2]
  real :: wxfersc        ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]
  real :: hxfersd        ! heat xfer from sea surface to seaspray drops this step [J/m^2]
  real :: hxferdc        ! heat xfer from seaspray drops to can_air this step [J/m^2]
  real :: wxferdc        ! vapor xfer from seaspray drops to can_air this step [kg_vap/m^2]
  real :: hxferse        ! heat xfer from sea surface to seaspray2 drops this step [J/m^2]
  real :: hxferec        ! heat xfer from seaspray2 drops to can_air this step [J/m^2]
  real :: wxferec        ! vapor xfer from seaspray2 drops to can_air this step [kg_vap/m^2]
  real :: sfc_rhovs      ! sat vapor density at sea surface temp [kg_vap/m^3]
  real :: spray_rhovs    ! sat vapor density at seaspray temp [kg_vap/m^3]
  real :: spray_rhovsp   ! sat vapor density derivative at seaspray temp [kg_vap/(K m^3)]
  real :: spray2_rhovs   ! sat vapor density at seaspray2 temp [kg_vap/m^3]
  real :: spray2_rhovsp  ! sat vapor density derivative at seaspray2 temp [kg_vap/(K m^3)]
  real :: can_rhov       ! Canopy air water vapor density [kg_vap/m^3]
  real :: canair         ! Canopy air mass [kg/m^2]
  real :: hcapcan        ! Canopy air heat capacity [J/(m^2 K)]
  real :: hcapspray      ! Heat capacity of sea spray [J/(m^2 K)]
  real :: hcapspray2     ! Heat capacity of sea spray2 [J/(m^2 K)]
  real :: canairi        ! Inverse of canair
  real :: hcapcani       ! Inverse of hcapcan
  real :: hcapsprayi     ! Inverse of hcapspray
  real :: hcapspray2i    ! Inverse of hcapspray2

  real    :: spray_massflux  ! Seaspray mass flux from sea [kg m^-2 s^-1]
  real    :: spray_mass      ! Seaspray steady-state mass in canopy [kg m^-2]
  real    :: spray_suodt     ! Seaspray steady-state vapor transfer coef [m s^-1]
  real    :: spray2_massflux ! Seaspray2 mass flux from sea [kg m^-2 s^-1]
  real    :: spray2_mass     ! Seaspray2 steady-state mass in canopy [kg m^-2]
  real    :: spray2_suodt    ! Seaspray2 steady-state vapor transfer coef [m s^-1]
  real    :: wt1, wt2        ! Interpolation weights [ ]
  integer :: ind             ! Interpolation index

  real(r8) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
  real(r8) :: h1, h2, h4, h5, h6, h7, h8
  real(r8) :: y1, y2, y3, y4, y5, y6, y7, y8, y9, y10

  real(r8) :: aa4(4,4), xx4(4), yy4(4)         ! 4x4 matrix equation terms
  real(r8) :: aa10(10,10), xx10(10), yy10(10)  ! 10x10 matrix equation terms

  logical :: sing

  ! Local variables used to evaluate wind speed at 10 m height

  real :: z10, press_z10, exner_z10, cantheta, canthetav, airthetav
  real :: ufree, wind_z10, theta_z10, rrv_z10

  ! Seaspray is represented in the "sea canopy" in an analogous manner to water
  ! on the surface of vegetation in the land canopy.  Heat and vapor are exchanged
  ! between seaspray and canopy air within an implicit solver that also incorporates
  ! heat and vapor exchange between the sea surface and canopy air, and between
  ! canopy air and the free atmosphere.  Three properties of seaspray are required:
  ! (1) the steady-state mass of seaspray droplets [kg m^-2] in the canopy air, (2)
  ! the steady-state vapor transfer coefficient [m s^-1] between seaspray and canopy
  ! air, and (3) the mass flux [kg m^-2 s^-1] of seaspray into the canopy air from
  ! the sea surface.  (By multiplying the mass flux by the temperature difference
  ! between new and returning sea spray droplets and by the specific heat of liquid
  ! water, we get the net energy flux from the sea surface to seaspray.)

  ! These three properties are tabulated below as a function of wind speed at 10 m
  ! height.  The tabulated values were generated in a separate program (ssgf.f90),
  ! which processed information presented in Figs 3i and 4a of Barr et al. (2023).
  ! ssgf.f90 numerically integrates over the (seastate-based) seaspray droplet size
  ! spectra in Fig. 4a to obtain all three properties.

  ! Vapor transfer coefficient suodt between sea spray and canopy air corresponds,
  ! after multiplication by the model timestep, to the quantity (su) in OLAM
  ! microphysics divided by that governs vapor flux between liquid droplets and air.
  ! suodt is computed from the same formula, except for the following differences:
  !
  ! (1) su in microphysics is dimensionless, while (suodt * dt_sea) has dimensions
  !     of [m] representing a vertical integration over "canopy" depth.
  !
  ! (2) su in microphysics represents an integral over individual droplet size
  !     categories (cloud, drizzle, rain) assuming a gamma size distribution for
  !     each, whereas suodt is obtained by numerically integrating over the empirical
  !     seaspray droplet spectra in Fig. 4a.

  ! 20-4000 MICRONS (band 1)
  real ::  ssgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  ssgf_massflux(5) = (/  0., .00007, .00135, .00702, .02769 /) ! [kg m^-2 s^-1]
  real ::  ssgf_mass    (5) = (/  0., .00110, .01270, .04760, .12738 /) ! [kg m^-2]
  real ::  ssgf_suodt   (5) = (/  0.,  0.059,  0.389,  1.086,  1.917 /) ! [m s^-1]

  ! 20-100 MICRONS (band 2)
  real ::  Assgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  Assgf_massflux(5) = (/  0., .00004, .00023, .00060, .00094 /) ! [kg m^-2 s^-1]
  real ::  Assgf_mass    (5) = (/  0., .00089, .00665, .01918, .03373 /) ! [kg m^-2]
  real ::  Assgf_suodt   (5) = (/  0.,  0.057,  0.339,  0.904,  1.481 /) ! [m s^-1]

  ! 20-200 MICRONS (band 3)
  real ::  Bssgf_windz10 (5) = (/ 15.,    20.,    30.,    40.,    50. /) ! [m s^-1]
  real ::  Bssgf_massflux(5) = (/  0., .00007, .00069, .00189, .00332 /) ! [kg m^-2 s^-1]
  real ::  Bssgf_mass    (5) = (/  0., .00108, .01021, .03078, .05793 /) ! [kg m^-2]
  real ::  Bssgf_suodt   (5) = (/  0.,  0.059,  0.378,  1.030,  1.744 /) ! [m s^-1]

  ! 100-4000 MICRONS (band 4)

  ! 200-4000 MICRONS (band 5)

  integer, parameter :: band = 35

  ! If doing seaspray and wind speed is at least 15 m/s, compute wind speed at 10 m level
  ! (wind_z10) USING SFLUXT AND SFLUXR VALUES EVALUATED ON THE PREVIOUS TIMESTEP.
  ! Otherwise, set wind_z10 to zero as a flag to not compute seaspray flux.

  if (nl%iseasprayflg > 0 .and. vels > nl%seaspray_vmin) then
     z10 = 10.

     press_z10 = prss - grav * z10 * rhos ! hydrostatic eqn.
     exner_z10 = (press_z10 * p00i) ** rocp

     cantheta  = cantemp / canexner
     canthetav = cantheta             * (1.0 + eps_virt * canrrv)
     airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * airrrv)

     ufree = (grav * sfcg%dzt_bot(iwsfc) * max(wthv,0.0) / airthetav) ** one3

     call sfclyr_profile (vels, rhos, canexner, ustar, sfluxt, sfluxr, &
                          sfcg%dzt_bot(iwsfc), rough, ufree, &
                          cantheta, canthetav, canrrv, airthetav, &
                          z10, wind_z10, theta_z10, rrv_z10)
  else
     wind_z10 = 0.
  endif

  ! Evaluate surface saturation vapor density and mixing ratio of sea surface

  sfc_rhovs    = rhovsl(seatc-273.15)
  surface_srrv = sfc_rhovs / rhos

  ! rdi = ustar/5 is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  ! Canopy air quantities

  can_rhov = canrrv * rhos
  canair   = rhos * can_depth
  canairi  = 1. / canair
  hcapcan  = cp * canair
  hcapcani = 1. / hcapcan

  ! Set up and solve a linear system of equations that use trapezoidal-implicit
  ! differencing.  The solution of the system consists of turbulent heat and
  ! water vapor fluxes between canopy air, the ocean surface, sea spray droplets
  ! (if present), and the free atmosphere, and the consequent changes to water
  ! and temperature of canopy air and temperature of sea spray.

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

  if (wind_z10 < nl%seaspray_vmin) then ! Case with NO SEA SPRAY (two uncoupled equations)

     ! Set up and solve 4x4 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

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

  else                                        ! Case WITH SEA SPRAY (10x10 matrix equation)

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

     elseif (band == 11) then
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Test: method to turn off spray effect; must keep nonzero spray_mass
!  spray_massflux = 0.
!  spray_suodt = 0.

!  spray2_massflux = 0.
!  spray2_suodt = 0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     a1  = spray_suodt * dt_sea
     a2  = cp * rhos * a1
     a7  = spray_massflux * cliq * dt_sea
     a3  = spray2_suodt * dt_sea
     a4  = cp * rhos * a3
     a8  = spray2_massflux * cliq * dt_sea

     spray_rhovs   = rhovsl(spraytemp        - 273.15)
     spray_rhovsp  = rhovsl(spraytemp  + 1.0 - 273.15) - spray_rhovs
     spray2_rhovs  = rhovsl(spray2temp       - 273.15)
     spray2_rhovsp = rhovsl(spray2temp + 1.0 - 273.15) - spray2_rhovs

     hcapspray   = spray_mass * cliq
     hcapsprayi  = 1. / hcapspray
     hcapspray2  = spray2_mass * cliq
     hcapspray2i = 1. / hcapspray2

     h1 = fcn * hcapsprayi * spray_rhovsp
     h5 = fcn * hcapsprayi
     h2 = fcn * hcapspray2i * spray2_rhovsp
     h6 = fcn * hcapspray2i

     y1  = spray_rhovs  - can_rhov
     y4  = spraytemp    - cantemp
     y3  = seatc        - spraytemp
     y6  = spray2_rhovs - can_rhov
     y8  = spray2temp   - cantemp
     y7  = seatc        - spray2temp

     ! Set up and solve 10x10 matrix equation (trapezoidal implicit method) to balance
     ! vapor and heat fluxes between canopy air, sea spray, and the ocean surface.

     aa10(1,1)  = 1._r8 + a1 * (h1 * alvl + h4)
     aa10(1,2)  =         a1 * h4
     aa10(1,3)  =         a1 * h1
     aa10(1,4)  = 0._r8
     aa10(1,5)  =       - a1 * h1
     aa10(1,6)  =         a1 * h4
     aa10(1,7)  = 0._r8
     aa10(1,8)  = 0._r8
     aa10(1,9)  =       - a1 * h4
     aa10(1,10) = 0._r8
     yy10(1)    =         a1 * y1  ! WDC row

     aa10(2,1)  =         a5 * h4
     aa10(2,2)  = 1._r8 + a5 * h4
     aa10(2,3)  = 0._r8
     aa10(2,4)  = 0._r8
     aa10(2,5)  = 0._r8
     aa10(2,6)  =         a5 * h4
     aa10(2,7)  = 0._r8
     aa10(2,8)  = 0._r8
     aa10(2,9)  =       - a5 * h4
     aa10(2,10) = 0._r8
     yy10(2)    =         a5 * y2  ! WSC row

     aa10(3,1)  =         a2 * h5 * alvl
     aa10(3,2)  = 0._r8
     aa10(3,3)  = 1._r8 + a2 * (h5 + h7)
     aa10(3,4)  =         a2 * h7
     aa10(3,5)  =       - a2 * h5
     aa10(3,6)  = 0._r8
     aa10(3,7)  =         a2 * h7
     aa10(3,8)  = 0._r8
     aa10(3,9)  = 0._r8
     aa10(3,10) =       - a2 * h7
     yy10(3)    =         a2 * y4  ! HDC row

     aa10(4,1)  = 0._r8
     aa10(4,2)  = 0._r8
     aa10(4,3)  =         a6 * h7
     aa10(4,4)  = 1._r8 + a6 * h7
     aa10(4,5)  = 0._r8
     aa10(4,6)  = 0._r8
     aa10(4,7)  =         a6 * h7
     aa10(4,8)  = 0._r8
     aa10(4,9)  = 0._r8
     aa10(4,10) =       - a6 * h7
     yy10(4)    =         a6 * y5  ! HSC row

     aa10(5,1)  =       - a7 * h5 * alvl
     aa10(5,2)  = 0._r8
     aa10(5,3)  =       - a7 * h5
     aa10(5,4)  = 0._r8
     aa10(5,5)  = 1._r8 + a7 * h5
     aa10(5,6)  = 0._r8
     aa10(5,7)  = 0._r8
     aa10(5,8)  = 0._r8
     aa10(5,9)  = 0._r8
     aa10(5,10) = 0._r8
     yy10(5)    =         a7 * y3  ! HSD row

     aa10(6,1)  =         a3 * h4
     aa10(6,2)  =         a3 * h4
     aa10(6,3)  = 0._r8
     aa10(6,4)  = 0._r8
     aa10(6,5)  = 0._r8
     aa10(6,6)  = 1._r8 + a3 * (h2 * alvl + h4)
     aa10(6,7)  =         a3 * h2
     aa10(6,8)  =       - a3 * h2
     aa10(6,9)  =       - a3 * h4
     aa10(6,10) = 0._r8
     yy10(6)    =         a3 * y6  ! WEC row

     aa10(7,1)  = 0._r8
     aa10(7,2)  = 0._r8
     aa10(7,3)  =         a4 * h7
     aa10(7,4)  =         a4 * h7
     aa10(7,5)  = 0._r8
     aa10(7,6)  =         a4 * h6 * alvl
     aa10(7,7)  = 1._r8 + a4 * (h6 + h7)
     aa10(7,8)  =       - a4 * h6
     aa10(7,9)  = 0._r8
     aa10(7,10) =       - a4 * h7
     yy10(7)    =         a4 * y8  ! HEC row

     aa10(8,1)  = 0._r8
     aa10(8,2)  = 0._r8
     aa10(8,3)  = 0._r8
     aa10(8,4)  = 0._r8
     aa10(8,5)  = 0._r8
     aa10(8,6)  =       - a8 * h6 * alvl
     aa10(8,7)  =       - a8 * h6
     aa10(8,8)  = 1._r8 + a8 * h6
     aa10(8,9)  = 0._r8
     aa10(8,10) = 0._r8
     yy10(8)    =         a8 * y7  ! HSE row

     aa10(9,1)  =       - a9 * h8
     aa10(9,2)  =       - a9 * h8
     aa10(9,3)  = 0._r8
     aa10(9,4)  = 0._r8
     aa10(9,5)  = 0._r8
     aa10(9,6)  =       - a9 * h8
     aa10(9,7)  = 0._r8
     aa10(9,8)  = 0._r8
     aa10(9,9)  = 1._r8 + a9 * h8
     aa10(9,10) = 0._r8
     yy10(9)    =         a9 * y9  ! WCA row

     aa10(10,1)  = 0._r8
     aa10(10,2)  = 0._r8
     aa10(10,3)  =       - a10 * h7
     aa10(10,4)  =       - a10 * h7
     aa10(10,5)  = 0._r8
     aa10(10,6)  = 0._r8
     aa10(10,7)  =       - a10  * h7
     aa10(10,8)  = 0._r8
     aa10(10,9)  = 0._r8
     aa10(10,10) = 1._r8 + a10 * h7
     yy10(10)    =         a10 * y10  ! HCA row

     call matrix8_NxN(10,aa10,yy10,xx10,sing); if (sing) call sing_print(iwsfc,'sea4',10,aa10,yy10,glatw,glonw)

     wxferdc = xx10(1)
     wxfersc = xx10(2)
     hxferdc = xx10(3)
     hxfersc = xx10(4)
     hxfersd = xx10(5)
     wxferec = xx10(6)
     hxferec = xx10(7)
     hxferse = xx10(8)
     wxferca = xx10(9)
     hxferca = xx10(10)

     cantemp    = cantemp    + (hxferdc + hxferec + hxfersc - hxferca) * hcapcani
     canrrv     = canrrv     + (wxferdc + wxferec + wxfersc - wxferca) * canairi
     spraytemp  = spraytemp  + (hxfersd - hxferdc - wxferdc * alvl)    * hcapsprayi
     spray2temp = spray2temp + (hxferse - hxferec - wxferec * alvl)    * hcapspray2i

     sfluxt = hxferca / dt_sea
     sfluxr = wxferca / dt_sea

     hfluxsea = (hxfersc + hxferdc) / (cp * dt_sea)

  endif

end subroutine seacell_2
