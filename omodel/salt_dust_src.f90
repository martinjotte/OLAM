subroutine sea_spray()

  ! Compute sea spray number concentration fluxes, applied as tendencies in
  ! the grid level adjacent to the sea surface.  Based on O'Dowd, C. D.,
  ! Smith, M. H., and Jennings, S. G., (1993): Submicron aerosol, radon and
  ! soot carbon characteristics over the North East Atlantic. J. Geophys.
  ! Res., 98, 1132-1136.  Updated for O'Dowd (1997, 1999).

  use mem_grid,    only: dzt_bot, volti
  use mem_basic,   only: rho
  use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use micro_coms,  only: igccn
  use ccnbin_coms, only: isalt
  use mem_micro,   only: ccntyp, con_gccn
  use mem_tend,    only: con_gccnt
  use mem_sea,     only: sea, omsea

  implicit none

  real, parameter :: timescale_salt = 3600.0 ! relaxation time [s]
  real, parameter :: depthscale_salt = 100.  ! assumed depth scale filled by fluxes [m]
  real, parameter :: dssotss = timescale_salt / depthscale_salt
  integer :: j, iw, isea, kw, jsfc, iwsfc, jasfc
  real :: vels10, film_diag, jet_diag, film_flux, jet_flux

  ! Return if neither ccn nor gccn are prognosed

  if (isalt < 1 .and. igccn < 2) return

  ! Loop over prognosed atmospheric grid columns

  !$omp parallel do private(iw, jsfc, iwsfc, jasfc, kw, isea, vels10, &
  !$omp                     film_diag, jet_diag, film_flux, jet_flux)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Check for sea area beneath this atmospheric grid column; cycle if none

     if (itab_w(iw)%jsea2 == 0) cycle

     ! Loop over sea cells beneath this atmospheric grid column

     do jsfc = itab_w(iw)%jsea1, itab_w(iw)%jsea2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)

        isea = iwsfc - omsea

        ! Skip if not seawater
        if (sfcg%leaf_class(iwsfc) /= 0) cycle

        ! Diagnose wind speed at 10 m height, and limit to a maximum of 20 m/s.

        vels10 = min(20., sfcg%vels(iwsfc) * log(10.         / sea%sea_rough(isea)) &
                                           / log(dzt_bot(kw) / sea%sea_rough(isea)))

        ! Diagnose equilibrium sea salt concentrations [#/m^3] for 10 m wind speed

        film_diag = 10.**(0.0950 * vels10 + 6.2830) ! Film mode diameter 0.1 micron
        jet_diag  = 10.**(0.0422 * vels10 + 5.7122) ! Jet mode diameter 1.0 micron
        !spm_diag = 10.**(0.0690 * vels10 + 0.1900)

        ! If model surface concentrations of sea salt are less than diagnosed
        ! values, compute surface fluxes [#/(m^2 s)] for film and jet modes
        ! that would bring model concentrations up to diagnosed concentrations
        ! in a time scale of timescale_salt, assuming that the model deficit
        ! applies and must be filled over a depth scale of depthscale_salt.
        ! Fixing the depth scale in this way makes the fluxes independent of
        ! the chosen vertical grid spacing.

        ! Convert surface fluxes to concentration tendency [#/(m^3 s)] in
        ! atmospheric cells that are adjacent to the sea surface.  Tendency is
        ! obtained by multiplying fluxes by surface area and dividing by grid
        ! cell volume.

        if (isalt > 0) then
           film_flux = dssotss * (film_diag - ccntyp(isalt)%con_ccn(kw,iw) * rho(kw,iw))

           if (film_flux > 0.) then
                ccntyp(isalt)%con_ccnt(kw,iw) &
              = ccntyp(isalt)%con_ccnt(kw,iw) &
              + film_flux *  itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw) * (1.0 - sea%seaicec(isea))
           endif
        endif

        if (igccn == 2) then
           jet_flux = dssotss * ( jet_diag - con_gccn(kw,iw) * rho(kw,iw))

           if (jet_flux > 0.) then
                con_gccnt(kw,iw) &
              = con_gccnt(kw,iw) &
              + jet_flux * itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw) * (1.0 - sea%seaicec(isea))
           endif
        endif

     enddo  ! jws/isea

  enddo  ! iw
  !$omp end parallel do

end subroutine sea_spray

!===============================================================================

subroutine dust_src()

  use mem_grid,    only: dzt_bot, volti
  use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use ccnbin_coms, only: idust1, idust2, idust3, idust4
  use mem_micro,   only: ccntyp
  use mem_land,    only: land, omland, nzg
  use nuclei_coms, only: nbin, fact, sp, nodust, amassi, uth

  implicit none

  integer :: j, iw, iland, kw
  integer :: leaf_class, ibin
  integer :: jsfc, iwsfc, jasfc
  real :: vels10, vels10_cm
  real :: gwet, wetfac, flux_m
  real :: wprime
  real :: flux_n(nbin)

  ! Return if no ccn dust categories are prognosed

  if (idust1 < 1 .and. idust2 < 1 .and. idust3 < 1 .and. idust4 < 1) return

  ! Loop over prognosed atmospheric grid columns

  !$omp parallel private(flux_n)
  !$omp do private(iw,jsfc,iwsfc,jasfc,kw,iland,leaf_class,gwet, &
  !$omp            vels10,wprime,wetfac,vels10_cm,flux_m)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Check for land area beneath this atmospheric grid column; cycle if none

     if (itab_w(iw)%jland2 == 0) cycle

     ! Loop over land cells beneath this atmospheric grid column

     do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
        iwsfc = itab_w(iw)%iwsfc(jsfc)
        jasfc = itab_w(iw)%jasfc(jsfc)

        kw = itab_wsfc(iwsfc)%kwatm(jasfc)

        iland = iwsfc - omland

        leaf_class = sfcg%leaf_class(iwsfc)

        ! If leaf_class is not type that can generate dust, skip this land cell

        if (nodust(leaf_class)) cycle

        ! Diagnose fractional soil wetness, and cycle if >= 0.5

        gwet = land%soil_water(nzg,iland) / land%wsat_vg(nzg,iland)

        if (gwet >= 0.5) cycle

        ! Diagnose wind speed at 10 m height, and limit to a maximum of 20 m/s.

        vels10 = min(20., sfcg%vels(iwsfc) * log(10.         / sfcg%rough(iwsfc)) &
                                           / log(dzt_bot(kw) / sfcg%rough(iwsfc)))

        vels10_cm = 100. * vels10

        ! Diagnose wprime, which is used in parameterization by
        ! Fecan et al. and is based on clay percentage using the formula
        ! wprime = 0.0014(%clay)^2 + 0.17(%clay).

        wprime = (14. * land%clay(nzg,iland) + 17.) * land%clay(nzg,iland)

        if (gwet * 100. > wprime) then
           wetfac = sqrt(1. + 1.21 * (gwet * 100. - wprime)**0.68)
        else
           wetfac = 1.
        endif

        do ibin = 1,nbin

           !> Mass flux [g cm^-2 s^-1] from Ginoux et al. 2001 (Eq. 2)

           flux_m = fact * sp(ibin) * (vels10_cm**2.) * (vels10_cm - uth(ibin) * wetfac)

           !> Convert mass flux to number flux [#/(m^2 s)]

           flux_n(ibin) = max(0.,flux_m * 1.e4 * amassi(ibin))
        enddo

        ! Convert surface fluxes to concentration tendency [#/(m^3 s)] in
        ! atmospheric cells that are adjacent to the land surface.  Tendency is
        ! obtained by multiplying fluxes by surface area and dividing by grid
        ! cell volume.

        if (idust1 > 0) ccntyp(idust1)%con_ccnt(kw,iw) &
                      = ccntyp(idust1)%con_ccnt(kw,iw) &
                      + (flux_n(1) + flux_n(2) + flux_n(3) + flux_n(4)) &
                      * itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw)

        if (idust2 > 0) ccntyp(idust2)%con_ccnt(kw,iw) &
                      = ccntyp(idust2)%con_ccnt(kw,iw) &
                      + flux_n(5) * itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw)

        if (idust3 > 0) ccntyp(idust3)%con_ccnt(kw,iw) &
                      = ccntyp(idust3)%con_ccnt(kw,iw) &
                      + flux_n(6) * itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw)

        if (idust4 > 0) ccntyp(idust4)%con_ccnt(kw,iw) &
                      = ccntyp(idust4)%con_ccnt(kw,iw) &
                      + flux_n(7) * itab_wsfc(iwsfc)%arc(jasfc) * volti(kw,iw)

     enddo  ! jwl/iland

  enddo  ! iw
  !$omp end do
  !$omp end parallel

end subroutine dust_src
