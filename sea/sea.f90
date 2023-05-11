subroutine seacells()

  use sea_coms,    only: nzi, iupdsst, s1900_sst, isstfile, nsstfiles,  &
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, msea, omsea
  use misc_coms,   only: s1900_sim, iparallel
  use mem_para,    only: myrank

  implicit none

! Local variables

  integer :: isea      ! sea cell loop counter
  integer :: iwsfc
  real    :: timefac_sst   ! fraction of elapsed time from past to future SST obs
  real    :: timefac_seaice   ! fraction of elapsed time from past to future SEA ICE obs

! Time interpolation factors for updating SST and SEA ICE

  timefac_sst    = 0.0
  timefac_seaice = 0.0

  if (iupdsst == 1 .and. nsstfiles > 1) then
     timefac_sst = (s1900_sim           - s1900_sst(isstfile-1))  &
                 / (s1900_sst(isstfile) - s1900_sst(isstfile-1))
  endif

  if (iupdseaice == 1 .and. nseaicefiles > 1) then
     timefac_seaice = (s1900_sim                 - s1900_seaice(iseaicefile-1)) &
                    / (s1900_seaice(iseaicefile) - s1900_seaice(iseaicefile-1))
  endif

! Loop over ALL SEA CELLS

  !$omp parallel do private (iwsfc)
  do isea = 2, msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

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

     ! Update SEA fields

     call seacell(isea, iwsfc,             &
                  sfcg%rhos      (iwsfc),  &
                  sea%sea_ustar   (isea),  &
                  sea%sea_sxfer_t (isea),  &
                  sea%sea_sxfer_r (isea),  &
                  sfcg%can_depth (iwsfc),  &
                  sea%seatc       (isea),  &
                  sea%sea_cantemp (isea),  &
                  sea%sea_canrrv  (isea),  &
                  sea%sea_sfc_srrv(isea),  &
                  sea%sea_rough   (isea)   )

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
                      sea%ice_wthv           (isea), &
                      sea%ice_sxfer_t        (isea), &
                      sea%ice_sxfer_r        (isea)  )

! If ice exists, compute seaice canopy fluxes

     if (sea%nlev_seaice(isea) > 0) then

        call seaice(sea%seaice_energy(1:nzi,isea), &
                    sea%seaice_tempk (1:nzi,isea), &
                    sea%nlev_seaice        (isea), &
                    sea%ice_net_rshort     (isea), &
                    sea%ice_net_rlong      (isea), &
                    sfcg%rhos             (iwsfc), &
                    sea%ice_ustar          (isea), &
                    sfcg%can_depth        (iwsfc), &
                    sea%ice_cantemp        (isea), &
                    sea%ice_canrrv         (isea), &
                    sea%ice_sfc_srrv       (isea), &
                    sea%ice_sxfer_t        (isea), &
                    sea%ice_sxfer_r        (isea)  )

        sfcg%rough    (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_rough  (isea) + &
                                       sea%seaicec(isea)  * sea%ice_rough  (isea)

        sfcg%cantemp  (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_cantemp(isea) + &
                                       sea%seaicec(isea)  * sea%ice_cantemp(isea)

        sfcg%canrrv   (iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_canrrv (isea) + &
                                       sea%seaicec(isea)  * sea%ice_canrrv (isea)

        sea%surface_srrv(isea) = (1.0 - sea%seaicec(isea)) * sea%sea_sfc_srrv(isea) + &
                                        sea%seaicec(isea)  * sea%ice_sfc_srrv(isea)

     else

        sfcg%rough      (iwsfc) = sea%sea_rough   (isea)
        sfcg%cantemp    (iwsfc) = sea%sea_cantemp (isea)
        sfcg%canrrv     (iwsfc) = sea%sea_canrrv  (isea)
        sea%surface_srrv (isea) = sea%sea_sfc_srrv(isea)

     endif

  enddo
  !$omp end parallel do

end subroutine seacells

!===============================================================================

subroutine seacell( isea, iwsfc, rhos, ustar, sxfer_t, sxfer_r, can_depth, &
                    seatc, cantemp, canrrv, surface_srrv, rough     )

  use mem_sfcg,    only: sfcg
  use mem_sea,     only: sea
  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav, cliq1000, erad, eradi
  use therm_lib,   only: rhovsl
  use pom2k1d,     only: pom, rhoref, pom_column

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
  integer, intent(in)    :: iwsfc        ! current sfcg cell index
  real,    intent(in)    :: rhos         ! air density [kg/m^3]
  real,    intent(in)    :: ustar        ! friction velocity [m/s]
  real,    intent(in)    :: sxfer_t      ! can_air to atm heat xfer this step [kg_air K/m^2]
  real,    intent(in)    :: sxfer_r      ! can_air to atm vapor xfer this step (kg_vap/m^2]
  real,    intent(in)    :: can_depth    ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: seatc        ! current sea temp (obs time) [K]
  real,    intent(inout) :: cantemp      ! "canopy" air temp [K]
  real,    intent(inout) :: canrrv       ! "canopy" air vapor mixing ratio [kg_vap/kg_dryair]
  real,    intent(out)   :: surface_srrv ! sea surface sat mixing ratio [kg_vap/kg_dryair]
  real,    intent(out)   :: rough        ! sea cell roughess height [m]

! Local parameter

  real, parameter :: z0fac = 0.011 / grav  ! factor for Charnok roughness height
  real, parameter :: ozo   = 1.59e-5       ! base roughness height in HWRF
  real, parameter :: one3  = 1./ 3.

! Local variables

  real :: rdi     ! canopy conductance [m/s]
  real :: hfluxgc ! heat flux from sea surface to can_air [kg K/(m^2 s)]
  real :: hxfergc ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxfergc ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]

  real :: zn1, zn2, zw, usti
  real :: raxis, windu, windv, cdtop, wusurf, wvsurf, wtsurf, wssurf, swrad

! Evaluate surface saturation specific humidity

  surface_srrv = rhovsl(seatc-273.15) / rhos

! Update temperature and vapor specific humidity of "canopy" from
! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  hfluxgc = rhos * rdi * (seatc - cantemp)

  hxfergc = dt_sea * cp * hfluxgc
  wxfergc = dt_sea * rhos * rdi * (surface_srrv - canrrv)

  hxferca = cp * sxfer_t  ! sxfer_t and sxfer_r already incorporate dt_sea

  cantemp = cantemp + (hxfergc - hxferca) / (can_depth * rhos * cp)
  canrrv  = canrrv  + (wxfergc - sxfer_r) / (can_depth * rhos)

! Evaluate sea roughness height

! Charnok (1955):
! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)

! Davis et al. (2008) originally used in HWRF
! rough = 10. * exp(-10. / ustar ** .333333)
! rough = max(.125e-6, min(2.85e-3,rough))

! 2012 HWRF scheme; interpolates between the Charnok scheme at low wind
! and the Davis et al. curve fit at high wind speeds

  usti  = 1.0 / ustar
  zw    = min( (ustar/1.06)**one3, 1.0 )
  zn1   = z0fac * ustar * ustar + ozo
  zn2   = 10. * exp(-9.5 * usti**one3) + 1.65e-6 * usti
  rough = (1.0-zw) * zn1 + zw * zn2
  rough = min( rough, 2.85e-3)

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
     wtsurf = hfluxgc / rhoref + (sfcg%rlong(iwsfc) - sfcg%rlongup(iwsfc)) / cliq1000
     wssurf = 0.
     swrad = sfcg%rshort(iwsfc) * (1. - sfcg%albedo_beam(iwsfc)) / cliq1000

     call pom_column(isea, pom%kba(isea), wusurf, wvsurf, wtsurf, wssurf, swrad)

!TMP     sea%seatc(isea) = pom%potmp(1,isea)

  endif

end subroutine seacell
