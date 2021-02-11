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
subroutine radiate()

  use mem_tend,    only: thilt
  use mem_ijtabs,  only: jtab_w, itab_w, itabg_w, mrl_begl, istp, jtw_prog
  use mem_sfcg,    only: itab_wsfc, sfcg, mwsfc
  use mem_lake,    only: lake, mlake, omlake
  use mem_land,    only: land, mland, omland, nzg
  use mem_sea,     only: sea, msea, omsea
  use sea_coms,    only: nzi
  use leaf4_surface,only: sfcrad_land
  use mem_radiate, only: solfac, sunx, suny, sunz, cosz, nadd_rad,          &
                         rlongup, rlong_albedo, albedt, albedt_beam,        &
                         albedt_diffuse, fthrd_sw, rshort, rlong, fthrd_lw, &
                         rshort_top, rshortup_top, rshort_diffuse,          &
                         rshort_clr, rshortup_clr,                          &
                         rshort_top_clr, rshortup_top_clr,                  &
                         par, par_diffuse, uva, uvb, uvc, pbl_cld_forc,     &
                         ppfd, ppfd_diffuse, rlong_ks, rshort_ks,           &
                         rshort_diffuse_ks, ppfd_ks, ppfd_diffuse_ks
  use mem_basic,   only: thil, theta, tair
  use consts_coms, only: stefan, pio180, eradi, r8
  use misc_coms,   only: io6, time8p, time_istp8, radfrq, ilwrtyp, &
                         iswrtyp, dtlong, iparallel, isubdomain, mstp, runtype
  use mem_grid,    only: lpw, mza, wnx, wny, wnz, lsw, nsw_max
  use mem_turb,    only: frac_sfc
  use therm_lib,   only: qtk
  use mem_para,    only: myrank
  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w

  implicit none

  integer :: j
  integer :: iw
  integer :: k
  integer :: nsfc, nzw
  integer :: isea, iland, ilake, iwsfc, jsfc, jasfc
  integer :: ka, ks, kw
  integer :: koff
  integer :: nrad
  integer :: mrl
  real    :: sea_cosz, lake_cosz
  real    :: wdepth
  real    :: wt, wti, wtk
  real    :: tempk, fracliq

  real    :: albedt_ks        (nsw_max)
  real    :: albedt_diffuse_ks(nsw_max)
  real    :: rlongup_ks       (nsw_max)
  real    :: rlong_albedo_ks  (nsw_max)

  ! Check whether it is time to update radiative fluxes and heating rates

  if ((istp == 1 .and. mod(time8p, radfrq) < dtlong) .or. &
      (istp == 1 .and. mstp == 0 .and. runtype == 'HISTADDGRID')) then

     ! Print message that radiative transfer is being computed

     write(io6, '(A,f0.2,A)') &
          ' Radiation tendencies updated at ', time_istp8/3600., &
          ' hrs into simulation'

     ! Compute components of unit vector pointing to sun
     ! Compute coefficient for solar constant to represent varying earth-sun distance.

     call sunloc()

     ! Loop over all SEA cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary. 

     !$omp parallel
     !$omp do private (iwsfc, sea_cosz)
     do isea = 2,msea
        iwsfc = isea + omsea

        ! Get surface radiative properties (albedos and rlongup) for each sea cell.
        ! Compute solar zenith angle for sea cells

        sea_cosz = (sfcg%xew(iwsfc) * sunx  &
                 +  sfcg%yew(iwsfc) * suny  &
                 +  sfcg%zew(iwsfc) * sunz) * eradi

        ! Water albedo from Atwater and Bell (1981).

        if (sea_cosz > .03) then
           sea%sea_albedo(isea) = min(max(-.0139 + .0467 * tan(acos(sea_cosz)),.03),.999)
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

        sfcg%rlong_albedo(iwsfc)   = 0.0  ! [water longwave albedo assumed to be zero]
        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)
     enddo
     !$omp end do nowait

     ! Loop over all LAKE cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary. 

     !$omp do private (iwsfc, lake_cosz, tempk, fracliq)
     do ilake = 2,mlake
        iwsfc = ilake + omlake

        ! Get surface radiative properties (albedos and rlongup) for each sea cell. 
        ! Compute solar zenith angle for sea cells

        lake_cosz = (sfcg%xew(iwsfc) * sunx  &
                  +  sfcg%yew(iwsfc) * suny  &
                  +  sfcg%zew(iwsfc) * sunz) * eradi

        ! Water albedo from Atwater and Bell (1981).

        if (lake_cosz > .03) then
           sfcg%albedo_beam(iwsfc) = min(max(-.0139 + .0467 * tan(acos(lake_cosz)),.03),.999)
        else
           sfcg%albedo_beam(iwsfc) = 0.999
        endif

        call qtk(lake%lake_energy(ilake),tempk,fracliq)

        sfcg%rlongup       (iwsfc) = stefan * tempk ** 4
        sfcg%rlong_albedo  (iwsfc)  = 0.0  ! [water longwave albedo assumed to be zero]
        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)
     enddo
     !$omp end do

     ! Loop over all LAND cells in subdomain, EVEN THOSE THAT ARE NOT PRIMARY,
     ! so that all surface albedos and upward longwave radiative fluxes are
     ! available beneath all ATM columns that are primary. 

     !$omp do private (iwsfc, nzw, wdepth)
     do iland = 2,mland
        iwsfc = iland + omland

        ! Get surface radiative properties (albedos and rlongup) for each land cell.

        sfcg%rshort(iwsfc) = 0.
        sfcg%rlong (iwsfc) = 0.

        nzw = max(land%nlev_sfcwater(iland), 1)

        if (land%nlev_sfcwater(iland) > 0) then
           wdepth = sum(land%sfcwater_depth(1:land%nlev_sfcwater(iland),iland))
        else
           wdepth = 0.
        endif

! CHANGES POSSIBLY MADE TO THIS CALL, PENDING QUESTION ON DELETING RSHORT_S

        call sfcrad_land(iland,                  &
           sfcg%leaf_class        (    iwsfc), &
           sfcg%xew               (    iwsfc), &
           sfcg%yew               (    iwsfc), &
           sfcg%zew               (    iwsfc), &
           sfcg%wnx               (    iwsfc), &
           sfcg%wny               (    iwsfc), &
           sfcg%wnz               (    iwsfc), &
           sfcg%rshort            (    iwsfc), &
           sfcg%rlong             (    iwsfc), &
           sfcg%rlongup           (    iwsfc), &
           sfcg%rlong_albedo      (    iwsfc), &
           sfcg%albedo_beam       (    iwsfc), &
           land%sfcwater_energy   (nzw,iland), &
           land%sfcwater_depth    (nzw,iland), &
           land%rshort_s          (nzw,iland), &
           land%rshort_g          (    iland), &
           land%rshort_v          (    iland), &
           land%rlong_g           (    iland), &
           land%rlong_s           (    iland), &
           land%rlong_v           (    iland), &
           land%nlev_sfcwater     (    iland), &
           land%veg_temp          (    iland), &
           land%veg_fracarea      (    iland), &
           land%veg_height        (    iland), &
           land%veg_albedo        (    iland), &
           land%snowfac           (    iland), &
           land%vf                (    iland), &
           land%cosz              (    iland), &
           land%soil_energy       (nzg,iland), &
           land%soil_water        (nzg,iland), &
           land%wresid_vg         (nzg,iland), &
           land%wsat_vg           (nzg,iland), &
           land%specifheat_drysoil(nzg,iland), &
           land%sand              (nzg,iland)  )

           ! For LEAF land cells, there is no distinction between beam and 
           ! diffuse radiation.

        sfcg%albedo_diffuse(iwsfc) = sfcg%albedo_beam(iwsfc)

     enddo
     !$omp end do nowait
     !$omp end parallel

     ! Loop over all radiative IW grid columns that are primary in this subdomain

     !$omp parallel private(rlongup_ks,rlong_albedo_ks,albedt_ks,albedt_diffuse_ks)
     !$omp do private (iw, ka, nsfc, jsfc, iwsfc, jasfc, wti, wtk, &
     !$omp                      kw, ks, koff, nrad) &
     !$omp             schedule(guided)
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j) ! jend(1) = hardw for mrl=1

        ka   = lpw(iw)
        nsfc = lsw(iw)

        ! Compute solar zenith angle for atmosphere cells

        cosz(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

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

        ! If level exists with no land/lake/sea cells, just copy from level below.

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

        if (ilwrtyp > 0 .or. (iswrtyp > 0 .and. cosz(iw) > 0.03)) then

           ! K index offset for radiation column arrays

           koff = ka - 1
           nrad = mza - koff + nadd_rad

           call rrtmg_raddriv( iw, ka, nrad, koff, nsfc, &
                             rlongup_ks, rlong_albedo_ks, albedt_ks, albedt_diffuse_ks )

        endif

     enddo
     !$omp end do
     !$omp end parallel

     ! MPI send/recv of surface downward radiative fluxes

     if (iparallel == 1) then
        call mpi_send_w(1, svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks)
        call mpi_recv_w(1, svara1=rlong_ks, svara2=rshort_ks, svara3=rshort_diffuse_ks)
     endif

     ! Loop over all SFC grid cells

     !$omp parallel
     !$omp do private (j, iw, kw, ka, ks, wt, iland)
     do iwsfc = 2,mwsfc

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        ! Prepare SFC grid downward shortwave and longwave fluxes for summation from ATM columns

        sfcg%rlong         (iwsfc) = 0.
        sfcg%rshort        (iwsfc) = 0.
        sfcg%rshort_diffuse(iwsfc) = 0.
        sfcg%rshort_clr    (iwsfc) = 0.

        if (sfcg%leaf_class(iwsfc) >= 2) then
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
           sfcg%rshort        (iwsfc) = sfcg%rshort        (iwsfc) + wt * rshort_ks        (ks,iw)
           sfcg%rshort_diffuse(iwsfc) = sfcg%rshort_diffuse(iwsfc) + wt * rshort_diffuse_ks(ks,iw)

           ! par and ppfd are members of land% and not sfcg%

           if (sfcg%leaf_class(iwsfc) >= 2) then

            ! land%par         (iland) = land%par         (iland) + wt * par_ks         (ks,iw)
            ! land%par_diffuse (iland) = land%par_diffuse (iland) + wt * par_diffuse_ks (ks,iw)
              land%ppfd        (iland) = land%ppfd        (iland) + wt * ppfd_ks        (ks,iw)
              land%ppfd_diffuse(iland) = land%ppfd_diffuse(iland) + wt * ppfd_diffuse_ks(ks,iw)

           endif

        enddo

     enddo   
     !$omp end do

     ! Loop over all SEA cells to compute radiative fluxes for all 
     ! seaice components, given that rshort and rlong are now updated.

    !$omp do private(iwsfc)
     do isea = 2, msea
        iwsfc = isea + omsea

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

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

     !$omp do private (iwsfc,nzw,wdepth)
     do iland = 2,mland
        iwsfc = iland + omland

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        nzw = max(land%nlev_sfcwater(iland), 1)

        if (land%nlev_sfcwater(iland) > 0) then
           wdepth = sum(land%sfcwater_depth(1:land%nlev_sfcwater(iland),iland))
        else
           wdepth = 0.
        endif

! CHANGES POSSIBLY MADE TO THIS CALL, PENDING QUESTION ON DELETING RSHORT_S

        call sfcrad_land(iland,                  &
           sfcg%leaf_class        (    iwsfc), &
           sfcg%xew               (    iwsfc), &
           sfcg%yew               (    iwsfc), &
           sfcg%zew               (    iwsfc), &
           sfcg%wnx               (    iwsfc), &
           sfcg%wny               (    iwsfc), &
           sfcg%wnz               (    iwsfc), &
           sfcg%rshort            (    iwsfc), &
           sfcg%rlong             (    iwsfc), &
           sfcg%rlongup           (    iwsfc), &
           sfcg%rlong_albedo      (    iwsfc), &
           sfcg%albedo_beam       (    iwsfc), &
           land%sfcwater_energy   (nzw,iland), &
           land%sfcwater_depth    (nzw,iland), &
           land%rshort_s          (nzw,iland), &
           land%rshort_g          (    iland), &
           land%rshort_v          (    iland), &
           land%rlong_g           (    iland), &
           land%rlong_s           (    iland), &
           land%rlong_v           (    iland), &
           land%nlev_sfcwater     (    iland), &
           land%veg_temp          (    iland), &
           land%veg_fracarea      (    iland), &
           land%veg_height        (    iland), &
           land%veg_albedo        (    iland), &
           land%snowfac           (    iland), &
           land%vf                (    iland), &
           land%cosz              (    iland), &
           land%soil_energy       (nzg,iland), &
           land%soil_water        (nzg,iland), &
           land%wresid_vg         (nzg,iland), &
           land%wsat_vg           (nzg,iland), &
           land%specifheat_drysoil(nzg,iland), &
           land%sand              (nzg,iland)  )

     enddo
     !$omp end do nowait
     !$omp end parallel

  endif

  ! Apply radiation tendencies in FTHRD to THILT

  mrl = mrl_begl(istp)
  if (mrl > 0) then
  !$omp parallel do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza

      if (tair(k,iw) > 253.) then
         thilt(k,iw) = thilt(k,iw) &
                     + (fthrd_sw(k,iw) + fthrd_lw(k,iw)) * theta(k,iw) / tair(k,iw)
      else
         thilt(k,iw) = thilt(k,iw) &
                     + (fthrd_sw(k,iw) + fthrd_lw(k,iw)) * thil(k,iw) / tair(k,iw)
      endif

     enddo

  enddo
  !$omp end parallel do
  endif

end subroutine radiate

!============================================================================

subroutine sunloc()

  use misc_coms,   only: io6, imonth1, idate1, iyear1, itime1, time_istp8,  &
                         iswrtyp, ilwrtyp

  use consts_coms, only: pi2, pio180

  use mem_radiate, only: jday, solfac, sunx, suny, sunz

  use mem_mclat,   only: mclat_spline

  implicit none

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

  integer, external :: julday

  ! Find current simulation date/time by adding elapsed simulation time to
  ! initial simulation date/time

  call date_add_to8(iyear1,imonth1,idate1,itime1*100  &
     ,time_istp8,'s',outyear,outmonth,outdate,outhour)

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

  call mclat_spline(jday)

end subroutine sunloc

!============================================================================

subroutine radinit()

  use mem_radiate,   only: maxadd_rad, nadd_rad, zmrad, mcica_seed, &
                           sunx, suny, sunz, cosz
  use mem_grid,      only: mza, zm, glatw, glonw, nwa, wnx, wny, wnz
  use misc_coms,     only: io6, iswrtyp, ilwrtyp
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

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

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
