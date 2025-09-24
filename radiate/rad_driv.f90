module raddriv

contains

!============================================================================

subroutine radiate()

  use mem_tend,     only: thilt
  use mem_ijtabs,   only: jtab_w, itab_w, mrl_begl, istp, jtw_prog
  use mem_sfcg,     only: itab_wsfc, sfcg, mwsfc
  use mem_lake,     only: lake, mlake, omlake
  use mem_land,     only: land, mland, omland, nzg, nzs
  use mem_sea,      only: sea, msea, omsea
  use sea_coms,     only: nzi
  use leaf4_surface,only: sfcrad_land, sfcrad_prep, sfcrad_rlongup
  use mem_radiate,  only: sunx, suny, sunz, cosz, cosz_rad, nadd_rad,        &
                          rlongup, rlong_albedo, albedt, albedt_beam,        &
                          albedt_diffuse, fthrd_sw, rshort, fthrd_lw,        &
                          rshort_top, rshortup_top, rshort_diffuse, dlong,   &
                          par, par_diffuse, uva, uvb, uvc, pbl_cld_forc,     &
                          ppfd, ppfd_diffuse, rlong_ks, rshort_ks,           &
                          rshort_diffuse_ks, ppfd_ks, ppfd_diffuse_ks, cosz_min
  use mem_basic,    only: theta, tair
  use consts_coms,  only: stefan, pio180, eradi, r8
  use misc_coms,    only: io6, time8p, time_istp8, radfrq, ilwrtyp, do_chem, &
                          iswrtyp, dtlong, iparallel, mstp, runtype
  use mem_grid,     only: lpw, mza, wnx, wny, wnz, lsw, nsw_max
  use therm_lib,    only: qtk
  use mem_para,     only: myrank
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_turb,     only: frac_land
  use rrtmg_driver, only: rrtmg_raddriv, update_rrtmg_lw_intermediate

  implicit none

  integer :: j
  integer :: iw
  integer :: k
  integer :: nsfc, nwatm
  integer :: isea, iland, ilake, iwsfc, jsfc, jasfc
  integer :: ka, ks, kw
  integer :: koff
  integer :: nrad, nrad_dt
  real    :: sfc_cosz, alpha, Zsurf, DZsurf
  real    :: wt, wti, wtk
  real    :: tempk, fracliq
  real(r8):: radtime

  real    :: albedt_ks        (nsw_max)
  real    :: albedt_diffuse_ks(nsw_max)
  real    :: rlongup_ks       (nsw_max)
  real    :: rlong_albedo_ks  (nsw_max)

  ! Check whether it is time to update radiative fluxes and heating rates

  if ((istp == 1 .and. mod(time8p, radfrq) < dtlong) .or. &
      (istp == 1 .and. mstp == 0 .and. runtype == 'HISTREGRID')) then

     ! Print message that radiative transfer is being computed

     write(io6, '(A,f0.2,A)') &
          ' Radiation tendencies updated at ', time_istp8/3600., &
          ' hrs into simulation'

     ! number of timesteps before radiation is called again

     nrad_dt = int( (radfrq - mod(time8p, radfrq)) / dtlong ) + 1
     radtime = nrad_dt * dtlong * 0.5_r8

     ! Compute components of unit vector pointing to sun
     ! Compute coefficient for solar constant to represent varying earth-sun distance.

     call sunloc( doffset=radtime )

     ! Loop over all SEA cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary.

     !$omp parallel
     !$omp do private(iwsfc, nwatm, sfc_cosz)
     do isea = 2, msea
        iwsfc = isea + omsea

        ! We can skip a sea cell that is not attached to a primary atm cell
        if (iparallel == 1) then
           nwatm = itab_wsfc(iwsfc)%nwatm
           if ( all( itab_w( itab_wsfc(iwsfc)%iwatm(1:nwatm) )%irank /= myrank ) ) cycle
        endif

        ! Get surface radiative properties (albedos and rlongup) for each sea cell.
        ! Compute solar zenith angle for sea cells

        sfc_cosz = ( sfcg%xew(iwsfc) * sunx &
                   + sfcg%yew(iwsfc) * suny &
                   + sfcg%zew(iwsfc) * sunz ) * eradi

        ! Water albedo from Atwater and Bell (1981).

        if (sfc_cosz > .03) then
           sea%sea_albedo(isea) = min(max(-.0139 + .0467 * tan(acos(sfc_cosz)),.03),.999)
        else
           sea%sea_albedo(isea) = 0.999
        endif

        sea%sea_rlongup(isea) = stefan * sea%seatc(isea) ** 4

        if (sea%nlev_seaice(isea) > 0) then

           ! Get seaice albedo and upward longwave

           call sfcrad_seaice_1( sea%ice_rlongup(isea),       &
                                 sea%ice_albedo(isea),        &
                                 sea%nlev_seaice(isea),       &
                                 sea%ice_cantemp(isea),       &
                                 sea%seaice_tempk(1:nzi,isea) )

           ! Average ice and water components based on seaice fraction

           sfcg%rlongup(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_rlongup(isea) + &
                                        sea%seaicec(isea)  * sea%ice_rlongup(isea)

           sfcg%albedo_beam(iwsfc) = (1.0 - sea%seaicec(isea)) * sea%sea_albedo(isea) + &
                                            sea%seaicec(isea)  * sea%ice_albedo(isea)

        else

           sea%ice_rlongup(isea) = 0.0
           sea%ice_albedo (isea) = 0.0

           sfcg%rlongup    (iwsfc) = sea%sea_rlongup(isea)
           sfcg%albedo_beam(iwsfc) = sea%sea_albedo (isea)

        endif

        sfcg%rlong_albedo  (iwsfc) = 0.0  ! [water longwave albedo assumed to be zero]
        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)
        sfcg%cosz_rad      (iwsfc) = max(sfc_cosz, 0.)

     enddo
     !$omp end do nowait

     ! Loop over all LAKE cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary.

     !$omp do private (iwsfc, nwatm, sfc_cosz, tempk, fracliq)
     do ilake = 2, mlake
        iwsfc = ilake + omlake

        ! We can skip a lake cell that is not attached to a primary atm cell
        if (iparallel == 1) then
           nwatm = itab_wsfc(iwsfc)%nwatm
           if ( all( itab_w( itab_wsfc(iwsfc)%iwatm(1:nwatm) )%irank /= myrank ) ) cycle
        endif

        ! Get surface radiative properties (albedos and rlongup) for each sea cell.
        ! Compute solar zenith angle for sea cells

        sfc_cosz = ( sfcg%xew(iwsfc) * sunx &
                   + sfcg%yew(iwsfc) * suny &
                   + sfcg%zew(iwsfc) * sunz ) * eradi

        ! Water albedo from Atwater and Bell (1981).

        if (sfc_cosz > .03) then
           sfcg%albedo_beam(iwsfc) = min(max(-.0139 + .0467 * tan(acos(sfc_cosz)),.03),.999)
        else
           sfcg%albedo_beam(iwsfc) = 0.999
        endif

        call qtk(lake%lake_energy(ilake),tempk,fracliq)

        sfcg%rlongup       (iwsfc) = stefan * tempk ** 4
        sfcg%rlong_albedo  (iwsfc)  = 0.0  ! [water longwave albedo assumed to be zero]
        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)
        sfcg%cosz_rad      (iwsfc) = max(sfc_cosz, 0.)

     enddo
     !$omp end do nowait

     ! Loop over all LAND cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary.

     !$omp do private(iwsfc,nwatm,sfc_cosz)
     do iland = 2, mland
        iwsfc = iland + omland

        ! We can skip a land cell that is not attached to a primary atm cell
        if (iparallel == 1) then
           nwatm = itab_wsfc(iwsfc)%nwatm
           if ( all( itab_w( itab_wsfc(iwsfc)%iwatm(1:nwatm) )%irank /= myrank ) ) cycle
        endif

        sfc_cosz = ( sfcg%xew(iwsfc) * sunx &
                   + sfcg%yew(iwsfc) * suny &
                   + sfcg%zew(iwsfc) * sunz ) * eradi

        call sfcrad_prep(iland, iwsfc,           &
             sfcg%leaf_class        (    iwsfc), &
             sfcg%wnx               (    iwsfc), &
             sfcg%wny               (    iwsfc), &
             sfcg%wnz               (    iwsfc), &
             land%skncomp           (    iland), &
             land%sfcwater_mass     (:,  iland), &
             land%sfcwater_energy   (:,  iland), &
             land%sfcwater_depth    (:,  iland), &
             land%veg_temp          (    iland), &
             land%veg_fracarea      (    iland), &
             land%veg_height        (    iland), &
             land%veg_albedo        (    iland), &
             land%soil_energy       (nzg,iland), &
             land%soil_water        (nzg,iland), &
             land%wresid_vg         (nzg,iland), &
             land%wsat_vg           (nzg,iland), &
             land%specifheat_drysoil(nzg,iland), &
             land%sand              (nzg,iland), &
             land%snowfac           (    iland), &
             land%vf                (    iland), &
             land%cosz              (    iland), &
             sfcg%rlongup           (    iwsfc), &
             sfcg%rlong_albedo      (    iwsfc), &
             sfcg%albedo_beam       (    iwsfc), &
             land%gnd_albedo        (    iwsfc), &
             land%gnd_emiss         (    iwsfc), &
             land%slong             (    iwsfc), &
             land%vlong             (    iwsfc)  )

           ! For LEAF land cells, there is no distinction between beam and
           ! diffuse radiation.

        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)
        sfcg%cosz_rad      (iwsfc) = max(sfc_cosz, 0.)

     enddo
     !$omp end do nowait
     !$omp end parallel

     ! Loop over all radiative IW grid columns that are primary in this subdomain

     !$omp parallel private(rlongup_ks,rlong_albedo_ks,albedt_ks,albedt_diffuse_ks)
     !$omp do private (iw, ka, nsfc, jsfc, iwsfc, jasfc, wti, wtk, &
     !$omp                      kw, ks, koff, nrad) &
     !$omp             schedule(guided)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ka   = lpw(iw)
        nsfc = lsw(iw)

        ! Compute solar zenith angle for atmosphere cells corresponding to
        ! the radiation time

        cosz_rad(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

        ! Zero out fields that are summed from SFC grid cells

        rlongup       (iw) = 0.0
        rlong_albedo  (iw) = 0.0
        albedt_beam   (iw) = 0.0
        albedt_diffuse(iw) = 0.0

        rlongup_ks       (1:nsfc) = 0.0
        rlong_albedo_ks  (1:nsfc) = 0.0
        albedt_ks        (1:nsfc) = 0.0
        albedt_diffuse_ks(1:nsfc) = 0.0

        ! Zero out solar radiation fields since solar radiation calls are skipped at night

        rshort        (iw) = 0.
        rshort_diffuse(iw) = 0.
        rshort_top    (iw) = 0.
        rshortup_top  (iw) = 0.

!        rshort_clr      (iw) = 0.
!        rshortup_clr    (iw) = 0.
!        rshort_top_clr  (iw) = 0.
!        rshortup_top_clr(iw) = 0.

        fthrd_sw(ka:mza,iw) = 0.

        par(iw) = 0.
        par_diffuse(iw) = 0.
        ppfd(iw) = 0.
        ppfd_diffuse(iw) = 0.
        uva(iw) = 0.
        uvb(iw) = 0.
        uvc(iw) = 0.

        pbl_cld_forc(iw) = 0.

        rlong_ks         (1:nsfc,iw) = 0.
        rshort_ks        (1:nsfc,iw) = 0.
        rshort_diffuse_ks(1:nsfc,iw) = 0.
!       par_ks           (1:nsfc,iw) = 0.
!       par_diffuse_ks   (1:nsfc,iw) = 0.
        ppfd_ks          (1:nsfc,iw) = 0.
        ppfd_diffuse_ks  (1:nsfc,iw) = 0.

        ! Loop over SFC grid cells beneath this ATM grid column and sum
        ! SFC grid radiative properties to this ATM column

        do jsfc = 1, itab_w(iw)%jsfc2
           iwsfc = itab_w(iw)%iwsfc(jsfc)
           jasfc = itab_w(iw)%jasfc(jsfc)

           wti = itab_wsfc(iwsfc)%arcoariw(jasfc)
           wtk = itab_wsfc(iwsfc)%arcoarkw(jasfc)

           kw  = itab_wsfc(iwsfc)%kwatm(jasfc)
           ks  = kw - ka + 1

           rlongup_ks       (ks) = rlongup_ks       (ks) + wtk * sfcg%rlongup       (iwsfc)
           rlong_albedo_ks  (ks) = rlong_albedo_ks  (ks) + wtk * sfcg%rlong_albedo  (iwsfc)
           albedt_ks        (ks) = albedt_ks        (ks) + wtk * sfcg%albedo_beam   (iwsfc)
           albedt_diffuse_ks(ks) = albedt_diffuse_ks(ks) + wtk * sfcg%albedo_diffuse(iwsfc)

           rlongup       (iw) = rlongup       (iw) + wti * sfcg%rlongup       (iwsfc)
           rlong_albedo  (iw) = rlong_albedo  (iw) + wti * sfcg%rlong_albedo  (iwsfc)
           albedt_beam   (iw) = albedt_beam   (iw) + wti * sfcg%albedo_beam   (iwsfc)
           albedt_diffuse(iw) = albedt_diffuse(iw) + wti * sfcg%albedo_diffuse(iwsfc)
        enddo

        ! If level exists with no land/lake/sea cell overlaps, just copy from level below.

        do ks = 2, nsfc
           if (rlongup_ks(ks) < 1.e-7) then
              rlongup_ks       (ks) = rlongup_ks       (ks-1)
              rlong_albedo_ks  (ks) = rlong_albedo_ks  (ks-1)
              albedt_ks        (ks) = albedt_ks        (ks-1)
              albedt_diffuse_ks(ks) = albedt_diffuse_ks(ks-1)
           endif
        enddo

        ! Set total surface albedo to surface direct albedo
        ! (LEAF doesn't differentiate between diffuse and direct albedo)

        albedt(iw) = albedt_beam(iw)

        ! Do RRTMg radiation if specified

        if (ilwrtyp > 0 .or. (iswrtyp > 0 .and. cosz_rad(iw) > cosz_min)) then

           ! K index offset for radiation column arrays

           koff = ka - 1
           nrad = mza - koff + nadd_rad

           call rrtmg_raddriv( iw, ka, nrad, koff, nsfc, &
                             rlongup_ks, rlong_albedo_ks, albedt_ks, albedt_diffuse_ks )

        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

     ! MPI send/recv of surface downward radiative fluxes

     if (iparallel == 1) then

        if (do_chem == 1) then
           call mpi_send_w(svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks, &
                           svara4=ppfd_ks, svara5=ppfd_diffuse_ks, r1dvara1=rshort_top)
           call mpi_recv_w(svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks, &
                           svara4=ppfd_ks, svara5=ppfd_diffuse_ks, r1dvara1=rshort_top)
        else
           call mpi_send_w(svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks, r1dvara1=rshort_top)
           call mpi_recv_w(svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks, r1dvara1=rshort_top)
        endif

     endif

     ! Compute components of unit vector pointing to sun at CURRENT time

     call sunloc( do_mclat=.false. )

     ! Loop over all SFC grid cells

     !$omp parallel
     !$omp do private (j, iw, kw, ka, ks, wt, iland, sfc_cosz, alpha, Zsurf, DZsurf)
     do iwsfc = 2,mwsfc

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        ! Prepare SFC grid downward shortwave and longwave fluxes for summation from ATM columns

        sfcg%rlong         (iwsfc) = 0.
        sfcg%rshort_rad    (iwsfc) = 0.
        sfcg%rshort_dif_rad(iwsfc) = 0.
        sfcg%rshort_toa_rad(iwsfc) = 0.

        if (do_chem == 1 .and. sfcg%leaf_class(iwsfc) >= 2) then
           iland = iwsfc - omland
           land%ppfd        (iland) = 0.
           land%ppfd_diffuse(iland) = 0.
        endif

        ! Loop over all ATM grid cells that are coupled to this SFC grid cell
        ! and sum ATM radiative fluxes to this SFCG cell

        do j = 1,itab_wsfc(iwsfc)%nwatm
           iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
           kw = itab_wsfc(iwsfc)%kwatm(j)
           ka = lpw(iw)
           ks = kw - ka + 1

           wt = itab_wsfc(iwsfc)%arcoarsfc(j)

           sfcg%rlong         (iwsfc) = sfcg%rlong         (iwsfc) + wt * rlong_ks         (ks,iw)
           sfcg%rshort_rad    (iwsfc) = sfcg%rshort_rad    (iwsfc) + wt * rshort_ks        (ks,iw)
           sfcg%rshort_dif_rad(iwsfc) = sfcg%rshort_dif_rad(iwsfc) + wt * rshort_diffuse_ks(ks,iw)
           sfcg%rshort_toa_rad(iwsfc) = sfcg%rshort_toa_rad(iwsfc) + wt * rshort_top          (iw)

           ! par and ppfd are members of land% and not sfcg%

           if (do_chem == 1 .and. sfcg%leaf_class(iwsfc) >= 2) then
            ! land%par         (iland) = land%par         (iland) + wt * par_ks         (ks,iw)
            ! land%par_diffuse (iland) = land%par_diffuse (iland) + wt * par_diffuse_ks (ks,iw)
              land%ppfd        (iland) = land%ppfd        (iland) + wt * ppfd_ks        (ks,iw)
              land%ppfd_diffuse(iland) = land%ppfd_diffuse(iland) + wt * ppfd_diffuse_ks(ks,iw)
           endif

        enddo

        ! scale solar radiation to current time

        if (sfcg%rshort_toa_rad(iwsfc) > .1) then

           sfc_cosz = ( sfcg%xew(iwsfc) * sunx &
                      + sfcg%yew(iwsfc) * suny &
                      + sfcg%zew(iwsfc) * sunz ) * eradi

           alpha  = max( sfc_cosz, cosz_min ) / max( sfcg%cosz_rad(iwsfc), cosz_min )

           Zsurf  = sfcg%rshort_rad(iwsfc) - sfcg%rshort_dif_rad(iwsfc)
           DZsurf = ( Zsurf / sfcg%rshort_toa_rad(iwsfc) )**(1./alpha) * sfcg%rshort_toa_rad(iwsfc) - Zsurf

           sfcg%rshort        (iwsfc) = alpha * (sfcg%rshort_rad    (iwsfc) + 0.5 * DZsurf)
           sfcg%rshort_diffuse(iwsfc) = alpha * (sfcg%rshort_dif_rad(iwsfc) - 0.5 * DZsurf)

        else

           sfcg%rshort        (iwsfc) = sfcg%rshort_rad    (iwsfc)
           sfcg%rshort_diffuse(iwsfc) = sfcg%rshort_dif_rad(iwsfc)

        endif

     enddo
     !$omp end do

     ! Loop over all SEA cells to compute radiative fluxes for all
     ! seaice components, given that rshort and rlong are now updated.

    !$omp do private(iwsfc)
     do isea = 2, msea
        iwsfc = isea + omsea

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        call sfcrad_seaice_2( sea%ice_net_rshort(isea), &
                              sea%ice_net_rlong (isea), &
                              sea%nlev_seaice   (isea), &
                              sfcg%rshort      (iwsfc), &
                              sfcg%rlong       (iwsfc), &
                              sea%ice_rlongup   (isea), &
                              sea%ice_albedo    (isea)  )

     enddo
     !$omp end do

     ! Loop over all LAND cells to compute radiative fluxes
     ! for all cell components, given that rshort and rlong are now updated.

     !$omp do private (iwsfc)
     do iland = 2,mland
        iwsfc = iland + omland

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        call sfcrad_land(iland, iwsfc, &
               sfcg%leaf_class(iwsfc), &
               sfcg%rshort    (iwsfc), &
               sfcg%rlong     (iwsfc), &
               land%slong     (iland), &
               land%vlong     (iland), &
               land%gnd_albedo(iland), &
               land%gnd_emiss (iland), &
               land%vf        (iland), &
               land%veg_albedo(iland), &
               land%rshort_s  (iland), &
               land%rlong_s   (iland), &
               land%rshort_v  (iland), &
               land%rlong_v   (iland)  )

     enddo
     !$omp end do nowait
     !$omp end parallel


  elseif (mrl_begl(istp) > 0) then

     ! Compute components of unit vector pointing to sun

     call sunloc( do_mclat=.false. )

     ! Update outgoing longwave from land cells due to change in skin temperature

     !$omp parallel
     !$omp do private(iwsfc)
     do iland = 2, mland
        iwsfc = iland + omland

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        call sfcrad_rlongup(iland, iwsfc, &
             sfcg%leaf_class        (    iwsfc), &
             land%skncomp           (    iland), &
             land%sfcwater_mass     (:,  iland), &
             land%sfcwater_energy   (:,  iland), &
             land%sfcwater_depth    (:,  iland), &
             land%veg_temp          (    iland), &
             land%soil_energy       (nzg,iland), &
             land%soil_water        (nzg,iland), &
             land%specifheat_drysoil(nzg,iland), &
             land%gnd_emiss         (    iland), &
             land%vf                (    iland), &
             sfcg%rlongup           (    iwsfc), &
             land%slong             (    iland), &
             land%vlong             (    iland)  )

     enddo
     !$omp end do

     !$omp do private(iw,ks,k)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Initialize rlong_ks to current downward longwave flux
        do ks = 1, lsw(iw)
           k = lpw(iw) + ks - 2
           rlong_ks(ks,iw) = dlong(k,iw)
        enddo

        if (frac_land(iw) > .01) then
           call update_rrtmg_lw_intermediate(iw, lpw(iw), lsw(iw))
        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

     if (iparallel == 1) then
        call mpi_send_w(svara1=rlong_ks)
        call mpi_recv_w(svara1=rlong_ks)
     endif

     ! Loop over all SFC grid cells

     !$omp parallel
     !$omp do private(j, iw, kw, wt, ks, sfc_cosz, alpha, Zsurf, DZsurf)
     do iwsfc = 2, mwsfc

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        ! Prepare SFC grid downward longwave fluxes for summation from ATM columns

        sfcg%rlong(iwsfc) = 0.

        ! Loop over all ATM grid cells that are coupled to this SFC grid cell
        ! and sum ATM radiative fluxes to this SFCG cell

        do j = 1,itab_wsfc(iwsfc)%nwatm
           iw = itab_wsfc(iwsfc)%iwatm    (j)  ! local index
           kw = itab_wsfc(iwsfc)%kwatm    (j)
           wt = itab_wsfc(iwsfc)%arcoarsfc(j)
           ks = kw - lpw(iw) + 1

           sfcg%rlong(iwsfc) = sfcg%rlong(iwsfc) + wt * rlong_ks(ks,iw)
        enddo

        ! scale solar radiation to current time

        if (sfcg%rshort_toa_rad(iwsfc) > .1) then

           sfc_cosz = ( sfcg%xew(iwsfc) * sunx &
                      + sfcg%yew(iwsfc) * suny &
                      + sfcg%zew(iwsfc) * sunz ) * eradi

           alpha  = max( sfc_cosz, cosz_min ) / max( sfcg%cosz_rad(iwsfc), cosz_min )

           Zsurf  = sfcg%rshort_rad(iwsfc) - sfcg%rshort_dif_rad(iwsfc)
           DZsurf = ( Zsurf / sfcg%rshort_toa_rad(iwsfc) )**(1./alpha) * sfcg%rshort_toa_rad(iwsfc) - Zsurf

           sfcg%rshort        (iwsfc) = alpha * (sfcg%rshort_rad    (iwsfc) + 0.5 * DZsurf)
           sfcg%rshort_diffuse(iwsfc) = alpha * (sfcg%rshort_dif_rad(iwsfc) - 0.5 * DZsurf)

        else

           sfcg%rshort        (iwsfc) = sfcg%rshort_rad    (iwsfc)
           sfcg%rshort_diffuse(iwsfc) = sfcg%rshort_dif_rad(iwsfc)

        endif

     enddo
     !$omp end do

     ! Loop over all SEA cells to compute radiative fluxes for all
     ! seaice components, given that rshort and rlong are now updated.

    !$omp do private(iwsfc)
     do isea = 2, msea
        iwsfc = isea + omsea

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        call sfcrad_seaice_2( sea%ice_net_rshort(isea), &
                              sea%ice_net_rlong (isea), &
                              sea%nlev_seaice   (isea), &
                              sfcg%rshort      (iwsfc), &
                              sfcg%rlong       (iwsfc), &
                              sea%ice_rlongup   (isea), &
                              sea%ice_albedo    (isea)  )

     enddo
     !$omp end do

     ! Loop over all LAND cells to compute radiative fluxes
     ! for all cell components, given that rshort and rlong are now updated.

     !$omp parallel do private(iwsfc,j,iw,kw,wt,ka,ks)
     do iland = 2, mland
        iwsfc = iland + omland

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK

        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        ! Update land canopy radiation budget

        call sfcrad_land(iland, iwsfc, &
               sfcg%leaf_class(iwsfc), &
               sfcg%rshort    (iwsfc), &
               sfcg%rlong     (iwsfc), &
               land%slong     (iland), &
               land%vlong     (iland), &
               land%gnd_albedo(iland), &
               land%gnd_emiss (iland), &
               land%vf        (iland), &
               land%veg_albedo(iland), &
               land%rshort_s  (iland), &
               land%rlong_s   (iland), &
               land%rshort_v  (iland), &
               land%rlong_v   (iland)  )
     enddo
     !$omp end parallel do

  endif

  ! Apply radiation tendencies to THILT at start of each long timestep

  if (mrl_begl(istp) > 0) then

     !$omp parallel do private(iw,k,cosz,alpha)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Current solar zenith angle for atm cells
        cosz(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

        alpha = max( cosz(iw), cosz_min ) / max( cosz_rad(iw), cosz_min )

        do k = lpw(iw), mza

           thilt(k,iw) = thilt(k,iw) &
                       + ( alpha * fthrd_sw(k,iw) + fthrd_lw(k,iw) ) * theta(k,iw) / tair(k,iw)
        enddo

     enddo
     !$omp end parallel do

  endif

end subroutine radiate

!============================================================================

subroutine sunloc( doffset, do_mclat )

  use misc_coms,   only: imonth1, idate1, iyear1, itime1, time_istp8
  use consts_coms, only: pi2, pio180, r8
  use mem_radiate, only: jday, solfac, sunx, suny, sunz
  use mem_mclat,   only: mclat_spline

  implicit none

  real(r8), optional, intent(in) :: doffset  ! time offset from current time (s)
  logical,  optional, intent(in) :: do_mclat ! switch to turn off preparing Mclatchy soundings

  integer :: outyear  ! current simulation year
  integer :: outmonth ! current simulation month
  integer :: outdate  ! current simulation date
  integer :: outhour  ! current simulation hour/min/sec (6 digits)
  integer :: ihour2   ! current simulation hour
  integer :: imin2    ! current simulation min
  integer :: isec2    ! current simulation sec

  real :: t1             ! 2 pi times fraction of year elapsed
  real :: t2             ! 2 pi times fraction of year elapsed with offset
  real :: declin         ! solar declination angle [deg]
  real :: eqn_of_time    ! equation of time solution [s]
  real :: d0             ! coefficient for solfac computation
  real :: d02            ! coefficient for solfac computation
  real :: utc_sec        ! seconds elapsed in current simulation day (UTC)
  real :: sun_longitude  ! longitude where sun is at zenith

  real(r8) :: stime
  logical  :: do_mclat_spline

  integer, external :: julday

  ! Handle optional arguments

  stime = time_istp8
  if (present(doffset)) stime = stime + doffset

  do_mclat_spline = .true.
  if (present(do_mclat)) do_mclat_spline = do_mclat

  ! Find current simulation date/time by adding elapsed simulation time to
  ! initial simulation date/time

  call date_add_to8( iyear1, imonth1, idate1, itime1*100, stime, 's', &
                     outyear, outmonth, outdate, outhour )

  ! Find current Julian day

  jday = julday(outmonth,outdate,outyear)

  ! Solfac is a multiplier of the solar constant to correct for Earth's
  ! varying distance to the sun.

  d0 = pi2 * real(jday-1) / 365.
  d02 = d0 * 2.
  solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
                    + 0.000719 * cos(d02) + 0.000077 * sin(d02)

  ! Declin is the solar latitude in degrees

  t1 = pi2 * real(jday) / 366.

  declin = .322003                 &
         - 22.971  * cos(t1)       &
         - .357898 * cos(t1 * 2.)  &
         - .14398  * cos(t1 * 3.)  &
         + 3.94638 * sin(t1)       &
         + .019334 * sin(t1 * 2.)  &
         + .05928  * sin(t1 * 3.)

  t2 = (279.134 + .985647 * real(jday)) * pio180

  ! The equation of time gives the number of seconds by which sundial time
  ! leads clock time

  eqn_of_time = 5.0323                  &
              - 100.976 * sin(t2)       &
              + 595.275 * sin(t2 * 2.)  &
              + 3.6858  * sin(t2 * 3.)  &
              - 12.47   * sin(t2 * 4.)  &
              - 430.847 * cos(t2)       &
              + 12.5024 * cos(t2 * 2.)  &
              + 18.25   * cos(t2 * 3.)

  ! Find the longitude where the sun is at zenith

  ihour2 = outhour / 10000
  imin2  = (outhour - 10000 * ihour2) / 100
  isec2  = outhour - 10000 * ihour2 - 100 * imin2

  utc_sec = real(ihour2) * 3600. + real(imin2) * 60. + real(isec2)

  sun_longitude = 180. - 360. * (utc_sec + eqn_of_time) / 86400.

  sunx = cos(declin * pio180) * cos(sun_longitude * pio180)
  suny = cos(declin * pio180) * sin(sun_longitude * pio180)
  sunz = sin(declin * pio180)

  ! Interpolate Mclatchy soundings between summer and winter values, and prepare
  ! spline coefficients for interpolation by latitude.

  if (do_mclat_spline) then
     call mclat_spline(jday)
  endif

end subroutine sunloc

!============================================================================

subroutine radinit()

  use mem_radiate,   only: maxadd_rad, nadd_rad, zmrad, mcica_seed, &
                           sunx, suny, sunz, cosz
  use mem_grid,      only: mza, zm, glatw, glonw, nwa, wnx, wny, wnz
  use misc_coms,     only: iswrtyp, ilwrtyp
  use consts_coms,   only: cp
  use rrtmg_sw_init, only: rrtmg_sw_ini
  use rrtmg_lw_init, only: rrtmg_lw_ini
  use rrtmg_cloud,   only: rsw_cld_optics_init, rlw_cloud_optics_init
  use clouds_gno,    only: gno_lookup_init
  use mem_ijtabs,    only: jtab_w, jtw_prog, itab_w

  implicit none

  real :: deltaz
  integer :: j, iw

  ! Compute NADD_RAD, the number of radiation levels to be added above the top
  ! model prognostic level.  (Added levels will be filled elsewhere with data
  ! from Mclatchy soundings.)

  ! (3/13/2013) Following recent recommendations that at least one level should
  ! always be added, we choose here to always add at least 5 levels AND to always
  ! carry the radiation up to at least 45 km.

  ! Estimate a reasonable value for the height increment between added levels.
  ! Make it approximately the increment between the two highest model levels,
  ! but no less than allowed by maxadd_rad, the maximum number of added levels.

  zmrad = 45.e3

  deltaz = max( zm(mza) - zm(mza-1), (zmrad - zm(mza)) / real(maxadd_rad) )

  zmrad = max(zmrad, zm(mza) + 5. * deltaz)

  nadd_rad = nint( (zmrad - zm(mza)) / deltaz )
  nadd_rad = max(nadd_rad,1)
  nadd_rad = min(nadd_rad,maxadd_rad)

  ! Initialize RRTMG s/w scheme

  if (iswrtyp > 0) then
     call rrtmg_sw_ini(cp)
     call rsw_cld_optics_init()
  endif

  ! Initialize RRTMG l/w scheme

  if (ilwrtyp > 0) then
     call rrtmg_lw_ini(cp)
     call rlw_cloud_optics_init()
  endif

  ! Initialize solar zenith angle information

  call sunloc()

  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Seed the random number generator for RRTMg's cloud overlap scheme.
     ! This should only be done at initial time, but if we restart a run
     ! where the seeds weren't saved we need to re-seed again

     if (all(mcica_seed(1:4,iw) == 0)) then
        mcica_seed(1,iw) = 1 + itab_w(iw)%iwglobe
        mcica_seed(2,iw) = 1 + itab_w(iw)%iwglobe + nwa
        mcica_seed(3,iw) = nint( (glatw(iw) +   1) * 1.e4 )
        mcica_seed(4,iw) = nint( (glonw(iw) + 181) * 1.e4 )
     endif

     ! Compute initial solar zenith angle for atmosphere cells

     cosz(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

  enddo

  ! Read in the convective cloud fraction lookup table

  call gno_lookup_init()

end subroutine radinit

!============================================================================

end module raddriv
