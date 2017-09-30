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
use mem_ijtabs,  only: jtab_w, itabg_w, mrl_begl, istp, jtw_prog, itab_w
use mem_leaf,    only: land, itab_wl
use mem_sea,     only: sea, itab_ws
use leaf_coms,   only: nzg, nzs, mwl
use sea_coms,    only: mws, nzi
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
                       iswrtyp, dtlong, isubdomain, mstp, runtype
use mem_grid,    only: lpw, mza, wnx, wny, wnz, lsw, nsw_max
use mem_turb,    only: frac_sfc
use ed_misc_coms,only: ed2_active

implicit none

integer :: j
integer :: iw
integer :: k
integer :: nsfc
integer :: iws, jws
integer :: iwl, jwl
integer :: ka, ks, kw
integer :: koff
integer :: nrad
integer :: mrl
real    :: arf_iw, arf_kw
real    :: sea_cosz

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

! Loop over all SEA cells

   !$omp parallel do private (sea_cosz)
   do iws = 2,mws

! Get surface radiative properties (albedos and rlongup) for each sea cell. 
! Compute solar zenith angle for sea cells

      sea_cosz = (sea%xew(iws) * sunx  &
               +  sea%yew(iws) * suny  &
               +  sea%zew(iws) * sunz) * eradi

      ! Water albedo from Atwater and Bell (1981).

      if (sea_cosz > .03) then
         sea%sea_albedo(iws) = min(max(-.0139 + .0467 * tan(acos(sea_cosz)),.03),.999)
      else
         sea%sea_albedo(iws) = 0.999
      endif

      sea%sea_rlongup(iws) = stefan * sea%seatc(iws) ** 4

      if (sea%nlev_seaice(iws) > 0) then

         ! Get seaice albedo and upward longwave

         call sfcrad_seaice_1( sea%ice_rlongup(iws),       &
                               sea%ice_albedo(iws),        &
                               sea%nlev_seaice(iws),       &
                               sea%ice_cantemp(iws),       &
                               sea%seaice_tempk(1:nzi,iws) )

         ! Average ice and water components based on seaice fraction

         sea%rlongup(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_rlongup(iws) + &
                                   sea%seaicec(iws)  * sea%ice_rlongup(iws)

         sea%albedo_beam(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_albedo(iws) + &
                                       sea%seaicec(iws)  * sea%ice_albedo(iws)

      else
         
         sea%ice_rlongup(iws) = 0.0
         sea%ice_albedo (iws) = 0.0

         sea%rlongup    (iws) = sea%sea_rlongup(iws)
         sea%albedo_beam(iws) = sea%sea_albedo(iws)

      endif

      sea%rlong_albedo(iws)   = 0.0  ! [water longwave albedo assumed to be zero]
      sea%albedo_diffuse(iws) = sea%albedo_beam(iws)

   enddo
   !$omp end parallel do

! Loop over all LAND cells

   !$omp parallel do
   do iwl = 2,mwl

! Get surface radiative properties (albedos and rlongup) for each land cell.

      if (land%ed_flag(iwl) == 0) then

! THIS IS ONLY FOR CELLS NOT RUNNING ED.

         land%rshort(iwl) = 0.
         land%rlong (iwl) = 0.
         
         call sfcrad_land(iwl,                                               &
            land%leaf_class    (      iwl), land%ntext_soil     (  nzg,iwl), &
            land%nlev_sfcwater (      iwl), land%sfcwater_energy(1:nzs,iwl), &
            land%sfcwater_depth(1:nzs,iwl), land%soil_energy    (  nzg,iwl), &
            land%soil_water    (  nzg,iwl), land%veg_temp       (      iwl), &
            land%veg_fracarea  (      iwl), land%veg_height     (      iwl), &
            land%veg_albedo    (      iwl), land%rshort         (      iwl), &
            land%rlong         (      iwl), land%rshort_s       (1:nzs,iwl), &
            land%rshort_g      (      iwl), land%rshort_v       (      iwl), &
            land%rlong_g       (      iwl), land%rlong_s        (      iwl), &
            land%rlong_v       (      iwl), land%rlongup        (      iwl), &
            land%rlong_albedo  (      iwl), land%albedo_beam    (      iwl), &
            land%snowfac       (      iwl), land%vf             (      iwl), &
            land%cosz          (      iwl), land%xew            (      iwl), &
            land%yew           (      iwl), land%zew            (      iwl), &
            land%wnx           (      iwl), land%wny            (      iwl), &
            land%wnz           (      iwl), land%flag_vg        (      iwl))

! For LEAF land cells, there is no distinction between beam and 
! diffuse radiation.

         land%albedo_diffuse(iwl) = land%albedo_beam(iwl)

      endif

   enddo
   !$omp end parallel do

! Do radiation for all ED grid cells.

#ifdef USE_ED2
   if (ed2_active == 1) then
      call ed_rad_wrapper(1)
   endif
#endif

! Loop over all radiative IW grid columns

!----------------------------------------------------------------------
   !$omp parallel do private (iw,ka,koff,nrad,nsfc,rlongup_ks,rlong_albedo,&
   !$omp                      albedt_ks,albedt_diffuse_ks)
   do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j) ! jend(1) = hardw for  mrl=1
!----------------------------------------------------------------------

      ka   = lpw(iw)
      nsfc = lsw(iw)

! Compute solar zenith angle for atmosphere cells

      cosz(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

! Zero out fields that are summed from land/sea cells

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

      rshort_clr      (iw) = 0.
      rshortup_clr    (iw) = 0.
      rshort_top_clr  (iw) = 0.
      rshortup_top_clr(iw) = 0.

      fthrd_sw(ka:mza,iw) = 0.

      par(iw) = 0.
      par_diffuse(iw) = 0.
      ppfd(iw) = 0.
      ppfd_diffuse(iw) = 0.
      uva(iw) = 0.
      uvb(iw) = 0.
      uvc(iw) = 0.
      pbl_cld_forc(iw) = 0.

! Sum sea cell contributions to this iw column

      do jws = 1, itab_w(iw)%nsea
         iws = itab_w(iw)%isea(jws)

         arf_iw = itab_ws(iws)%arf_iw     ! sea cell area frac of atm IW col
         arf_kw = itab_ws(iws)%arf_kw     ! sea cell area frac of atm (KW,IW)
                                          !    cell contact with surface
         kw = itab_ws(iws)%kw
         ks = kw - ka + 1

         rlongup_ks       (ks) = rlongup_ks       (ks) + arf_kw * sea%rlongup(iws)
         rlong_albedo_ks  (ks) = rlong_albedo_ks  (ks) + arf_kw * sea%rlong_albedo(iws)
         albedt_ks        (ks) = albedt_ks        (ks) + arf_kw * sea%albedo_beam(iws)
         albedt_diffuse_ks(ks) = albedt_diffuse_ks(ks) + arf_kw * sea%albedo_diffuse(iws)

         rlongup       (iw) = rlongup       (iw) + arf_iw * sea%rlongup(iws)
         rlong_albedo  (iw) = rlong_albedo  (iw) + arf_iw * sea%rlong_albedo(iws)
         albedt_beam   (iw) = albedt_beam   (iw) + arf_iw * sea%albedo_beam(iws)
         albedt_diffuse(iw) = albedt_diffuse(iw) + arf_iw * sea%albedo_diffuse(iws)
      enddo

! Sum land cell contributions to this iw column

      do jwl = 1, itab_w(iw)%nland
         iwl = itab_w(iw)%iland(jwl)

         arf_iw = itab_wl(iwl)%arf_iw      ! land cell area frac of atm IW col
         arf_kw = itab_wl(iwl)%arf_kw      ! land cell area frac of atm (KW,IW)
                                           !    cell contact with surface
         kw = itab_wl(iwl)%kw
         ks = kw - ka + 1

         rlongup_ks       (ks) = rlongup_ks       (ks) + arf_kw * land%rlongup(iwl)
         rlong_albedo_ks  (ks) = rlong_albedo_ks  (ks) + arf_kw * land%rlong_albedo(iwl)
         albedt_ks        (ks) = albedt_ks        (ks) + arf_kw * land%albedo_beam(iwl)
         albedt_diffuse_ks(ks) = albedt_diffuse_ks(ks) + arf_kw * land%albedo_diffuse(iwl)

         rlongup       (iw) = rlongup       (iw) + arf_iw * land%rlongup(iwl)
         rlong_albedo  (iw) = rlong_albedo  (iw) + arf_iw * land%rlong_albedo(iwl)
         albedt_beam   (iw) = albedt_beam   (iw) + arf_iw * land%albedo_beam(iwl)
         albedt_diffuse(iw) = albedt_diffuse(iw) + arf_iw * land%albedo_diffuse(iwl)
      enddo

! If no land/sea cells are at this level, just copy from level below.

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
!$omp end parallel do

! If running leaf, loop over SEA cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to sea cells

!$omp parallel do private (iw,ks)
   do iws = 2,mws

      iw = itab_ws(iws)%iw   ! global index

! If run is parallel, get local rank indices

      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      ks = itab_ws(iws)%kw - lpw(iw) + 1

      sea%rlong         (iws) = rlong_ks         (ks,iw)
      sea%rshort        (iws) = rshort_ks        (ks,iw)
      sea%rshort_diffuse(iws) = rshort_diffuse_ks(ks,iw)

   enddo
!$omp end parallel do

! If running leaf, loop over LAND cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to land cells

!$omp parallel do private (iw,ks)
   do iwl = 2,mwl

      iw = itab_wl(iwl)%iw   ! global index

! If run is parallel, get local rank indices

      if (isubdomain == 1) then
         iw = itabg_w(iw)%iw_myrank
      endif

      ks = itab_wl(iwl)%kw - lpw(iw) + 1

      land%rlong         (iwl) = rlong_ks         (ks,iw)
      land%rshort        (iwl) = rshort_ks        (ks,iw)
      land%rshort_diffuse(iwl) = rshort_diffuse_ks(ks,iw)
!     land%par           (iwl) = par_ks           (ks,iw)
!     land%par_diffuse   (iwl) = par_diffuse_ks   (ks,iw)
      land%ppfd          (iwl) = ppfd_ks          (ks,iw)
      land%ppfd_diffuse  (iwl) = ppfd_diffuse_ks  (ks,iw)

   enddo   
!$omp end parallel do

! Loop over all SEA cells to compute radiative fluxes for all 
! seaice components, given that rshort and rlong are now updated.

   !$omp parallel do
   do iws = 2, mws

      call sfcrad_seaice_2( sea%ice_net_rshort(iws), &
                            sea%ice_net_rlong (iws), &
                            sea%nlev_seaice   (iws), &
                            sea%rshort        (iws), &
                            sea%rlong         (iws), &
                            sea%ice_rlongup   (iws), &
                            sea%ice_albedo    (iws)  )

   enddo
   !$omp end parallel do

! If running leaf, loop over all LAND cells to compute radiative fluxes 
! for all cell components, given that rshort and rlong are now updated.

!$omp parallel do
   do iwl = 2,mwl

      if (land%ed_flag(iwl) == 0) then

! THIS IS ONLY FOR CELLS NOT RUNNING ED.

         call sfcrad_land(iwl,                                               &
            land%leaf_class    (      iwl), land%ntext_soil     (  nzg,iwl), &
            land%nlev_sfcwater (      iwl), land%sfcwater_energy(1:nzs,iwl), &
            land%sfcwater_depth(1:nzs,iwl), land%soil_energy    (  nzg,iwl), &
            land%soil_water    (  nzg,iwl), land%veg_temp       (      iwl), &
            land%veg_fracarea  (      iwl), land%veg_height     (      iwl), &
            land%veg_albedo    (      iwl), land%rshort         (      iwl), &
            land%rlong         (      iwl), land%rshort_s       (1:nzs,iwl), &
            land%rshort_g      (      iwl), land%rshort_v       (      iwl), &
            land%rlong_g       (      iwl), land%rlong_s        (      iwl), &
            land%rlong_v       (      iwl), land%rlongup        (      iwl), &
            land%rlong_albedo  (      iwl), land%albedo_beam    (      iwl), &
            land%snowfac       (      iwl), land%vf             (      iwl), &
            land%cosz          (      iwl), land%xew            (      iwl), &
            land%yew           (      iwl), land%zew            (      iwl), &
            land%wnx           (      iwl), land%wny            (      iwl), &
            land%wnz           (      iwl), land%flag_vg        (      iwl))

      endif

   enddo
!$omp end parallel do

! Do radiation for all ED grid cells.

#ifdef USE_ED2
   if (ed2_active == 1) then
      call ed_rad_wrapper(2)
   endif
#endif

endif

! Apply radiation tendencies in FTHRD to THILT

!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

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

! Update ED output variables

#ifdef USE_ED2
if (ed2_active == 1) then
   call ed_rad_wrapper(3)
endif
#endif

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

  use mem_radiate,   only: maxadd_rad, nadd_rad, zmrad, mcica_seed
  use mem_grid,      only: mza, zm, glatw, glonw, nwa
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

! Seed the random number generator for RRTMg's cloud overlap scheme.
! This should only be done at initial time, but if we restart a run
! where the seeds weren't saved we need to re-seed again

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
     if (all(mcica_seed(1:4,iw) == 0)) then
        mcica_seed(1,iw) = 1 + itab_w(iw)%iwglobe
        mcica_seed(2,iw) = 1 + itab_w(iw)%iwglobe + nwa
        mcica_seed(3,iw) = nint( (glatw(iw) +   1) * 1.e4 )
        mcica_seed(4,iw) = nint( (glonw(iw) + 181) * 1.e4 )
     endif
  enddo

! Read in the convective cloud fraction lookup table

  call gno_lookup_init()

end subroutine radinit
