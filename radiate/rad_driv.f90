!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine radiate()

use mem_tend,    only: thilt
use mem_ijtabs,  only: jtab_w, itabg_w, mrl_begl, istp
use mem_leaf,    only: land, itabg_wl, itab_wl, first_site
use mem_sea,     only: sea, itabg_ws, itab_ws
use leaf_coms,   only: nzg, nzs, mwl
use sea_coms,    only: mws
use mem_radiate, only: solfac, sunx, suny, sunz, cosz, nadd_rad,    &
                       rlongup, rlong_albedo, albedt, albedt_beam,  &
                       albedt_diffuse, fthrd, rshort, rlong, fthrd_lw,  &
                       rshort_top, rshortup_top, rshort_diffuse
use mem_basic,   only: rho
use micro_coms,  only: level
use consts_coms, only: stefan, pio180, eradi
use misc_coms,   only: io6, time_istp8, radfrq, itime1, ilwrtyp, iswrtyp,   &
                       dtlong, current_time, iparallel
use mem_grid,    only: lpw, mza, mwa
use mem_sflux,   only: mseaflux, seaflux, jseaflux,  &
                       mlandflux, landflux, jlandflux
use mem_grid,    only: wnx, wny, wnz
use mem_para,    only: myrank

use ed_structure_defs

!$ use omp_lib

implicit none

integer :: j
integer :: iw
integer :: k
integer :: isf
integer :: ilf
integer :: iws
integer :: iwl
integer :: ka
integer :: koff
integer :: nrad
integer :: mrl

real :: water_albedo
real :: arf_atm
real :: arf_land
real :: arf_sea
real :: flux
real :: sea_cosz

type(site), pointer :: ed_site
type(patch), pointer :: ed_patch

integer, external :: julday
integer :: jday
real :: rlong_previous(mwa)

! Check whether it is time to update radiative fluxes and heating rates

if (istp == 1 .and. mod(time_istp8 + .001d0,dble(radfrq)) < dtlong) then

! Print message that radiative transfer is being computed

   write(io6, '(A,f0.2,A)') &
        ' Radiation tendencies updated at ', time_istp8/3600., &
        ' hrs into simulation'

! Compute components of unit vector pointing to sun
! Compute coefficient for solar constant to represent varying earth-sun distance.

   call sunloc()

! Loop over all radiative IW grid cells where radiation may be done

   call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,k)
   do j = 1,jtab_w(12)%jend(1); iw = jtab_w(12)%iw(j)! jend(1) = hardw for mrl = 1
!----------------------------------------------------------------------
   call qsub('W',iw)

! Compute solar zenith angle for atmosphere cells

      cosz(iw) = wnx(iw) * sunx + wny(iw) * suny + wnz(iw) * sunz

! Zero out rlongup, albedt, rshort, rlong, and fthrd prior to summing over 
! land/sea flux cells

      rlong_previous(iw) = rlong(iw)
      rlong         (iw) = 0.
      rlongup       (iw) = 0.
      rlong_albedo  (iw) = 0.
      albedt        (iw) = 0.
      albedt_beam   (iw) = 0.
      rshort        (iw) = 0.
      rshort_diffuse(iw) = 0.
      albedt_diffuse(iw) = 0.
      rshort_top    (iw) = 0.
      rshortup_top  (iw) = 0.
      
      do k = lpw(iw),mza-1
         fthrd(k,iw) = 0.
         fthrd_lw(k,iw) = 0.
      enddo

   enddo
!$omp end parallel do
   call rsub('Wa',12)

! If running leaf3, loop over all SEA cells.

!$omp parallel do private (sea_cosz,water_albedo)
   do iws = 2,mws

! Skip IWS cell if running in parallel and primary rank of IWS /= MYRANK

      if (iparallel == 1) then
         if (itab_ws(iws)%irank /= myrank) cycle
      endif

! Zero out sea cell downward radiative fluxes prior to summation 
! over flux cells.

      sea%rlong(iws) = 0.                           
      sea%rshort(iws) = 0.                           
      sea%rshort_diffuse(iws) = 0.

! Get surface radiative properties (albedos and rlongup) for each sea cell. 
! Compute solar zenith angle for sea cells

      sea_cosz = (sea%xew(iws) * sunx  &
               +  sea%yew(iws) * suny  &
               +  sea%zew(iws) * sunz) * eradi

      water_albedo = .999

      if (nint(sea%seaicec(iws)) == 0) then 

! Water albedo from Atwater and Bell (1981).

         if (sea_cosz > .03) water_albedo =  &
              min(max(-.0139 + .0467 * tan(acos(sea_cosz)),.03),.999)

      else
         ! In the Los Alamos sea ice model, visible albedo is 0.78 and 
         ! NIR albedo is 0.36 for temperatures below -1C.  Here, we assume
         ! equal parts visible and NIR, yielding an albedo of 0.57.  This
         ! albedo decreases to 0.5 as the temperature increases to 0C.

         if(sea%seatc(iws) < 272.15)then
            water_albedo = 0.72  ! Dave modification
!            water_albedo = 0.57
         else
            water_albedo = 0.57 - 0.07 * (min(273.15,sea%seatc(iws)) - 272.15) 
         endif

      endif

      sea%rlongup(iws) = stefan * sea%seatc(iws) ** 4

      sea%rlong_albedo(iws) = 0.  ! [water longwave albedo assumed to be zero]
      
      sea%albedo_beam(iws)    = water_albedo 
      sea%albedo_diffuse(iws) = sea%albedo_beam(iws)

   enddo
!$omp end parallel do

! Do parallel send of SEA albedos and rlongup

   if (iparallel == 1) then
      call mpi_send_ws('R')
   endif

! If running leaf3, loop over all LAND cells.

!$omp parallel do
   do iwl = 2,mwl

! Skip IWL cell if running in parallel and primary rank of IWL /= MYRANK

      if (iparallel == 1) then
         if (itab_wl(iwl)%irank /= myrank) cycle
      endif

! Zero out land cell downward radiative fluxes prior to summation 
! over flux cells.

      land%rlong(iwl) = 0.                           
      land%rshort(iwl) = 0.                           
      land%rshort_diffuse(iwl) = 0.

! Get surface radiative properties (albedos and rlongup) for each land cell.

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
            land%wnz           (      iwl)                                 )

! For LEAF3 land cells, there is no distinction between beam and 
! diffuse radiation.

         land%albedo_diffuse(iwl) = land%albedo_beam(iwl)

      endif

   enddo
!$omp end parallel do

! Set the pointer for the first ED site

   ed_site => first_site

   do while(associated(ed_site))

! Loop over subgrid-scale patches

      ed_patch => ed_site%oldest_patch
      do while(associated(ed_patch))

! Get the unnormalized radiative transfer information

         call sfcrad_ed(land%cosz(ed_site%iland), ed_patch,   &
              ed_patch%cohort_count)

         ed_patch => ed_patch%younger
      enddo

! Copy results over to iwl (iland)

      call ed2land_radiation(ed_site)

! At this point, we know the beam and diffuse albedo, but since they
! are in general different, we do not know the net albedo.

      ed_site => ed_site%next_site

   enddo

! Do parallel recv of SEA albedos and rlongup

   if (iparallel == 1) then  
      call mpi_recv_ws('R')
   endif

! Do parallel send of LAND albedos and rlongup

   if (iparallel == 1) then  
      call mpi_send_wl('R')
   endif

! If running leaf3, loop over all SEAFLUX cells to get mean surface radiative
! properties for each IW grid cell.

!$omp parallel do private (isf,iw,iws,arf_atm)
   do j = 1,jseaflux(1)%jend(1)
      isf = jseaflux(1)%iseaflux(j)
      iw  = seaflux(isf)%iw        ! global index
      iws = seaflux(isf)%iwls      ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
         iws = itabg_ws(iws)%iws_myrank
      endif

      arf_atm = seaflux(isf)%arf_atm

      rlongup       (iw) = rlongup       (iw) + arf_atm * sea%rlongup(iws)       
      rlong_albedo  (iw) = rlong_albedo  (iw) + arf_atm * sea%rlong_albedo(iws)  
      albedt_beam   (iw) = albedt_beam   (iw) + arf_atm * sea%albedo_beam(iws)   
      albedt_diffuse(iw) = albedt_diffuse(iw) + arf_atm * sea%albedo_diffuse(iws)

   enddo
!$omp end parallel do

! Do parallel recv of LAND albedos and rlongup

   if (iparallel == 1) then  
      call mpi_recv_wl('R')
   endif

! If running leaf3, loop over all LANDFLUX cells to get mean surface radiative
! properties for each IW grid cell.

!$omp parallel do private (ilf,iw,iwl,arf_atm)
   do j = 1,jlandflux(1)%jend(1)
      ilf = jlandflux(1)%ilandflux(j)
      iw  = landflux(ilf)%iw         ! global index
      iwl = landflux(ilf)%iwls       ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
         iwl = itabg_wl(iwl)%iwl_myrank
      endif

      arf_atm = landflux(ilf)%arf_atm

      rlongup       (iw) = rlongup(iw)         &
                         + arf_atm * land%rlongup(iwl)
      rlong_albedo  (iw) = rlong_albedo(iw)    &
                         + arf_atm * land%rlong_albedo(iwl)

      albedt_beam   (iw) = albedt_beam(iw)     &
                         + arf_atm * land%albedo_beam(iwl)
      albedt_diffuse(iw) = albedt_diffuse(iw)  &
                         + arf_atm * land%albedo_diffuse(iwl)

   enddo
!$omp end parallel do

! Loop over all radiative IW grid columns

   call psub()
!----------------------------------------------------------------------
!$omp parallel do private (iw,ka,koff,nrad)
   do j = 1,jtab_w(12)%jend(1); iw = jtab_w(12)%iw(j) ! jend(1) = hardw for  mrl=1
!----------------------------------------------------------------------
   call qsub('W',iw)

! K index of lowest predicted level

      ka = lpw(iw)
      
! K index offset for radiation column arrays

      koff = ka - 2
   
! Do Mahrer-Pielke and/or Chen-Cotton radiation if specified

      if (ilwrtyp == 1 .or. ilwrtyp == 2 .or.  &
          iswrtyp == 1 .or. iswrtyp == 2) then  

         call ccmp_raddriv(iw,ka,koff)
      endif

! Do Harrington radiation if specified

      if (ilwrtyp == 3 .or. iswrtyp == 3) then
         nrad = mza - 1 - koff + nadd_rad

         call harr_raddriv(iw,ka,nrad,koff,rlong_previous(iw))
      endif

   enddo
!$omp end parallel do
   call rsub('Wb',12)

! If running leaf3, loop over SEAFLUX cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to sea cells

!$omp parallel do private (isf,iw,iws,arf_sea)
   do j = 1,jseaflux(1)%jend(1)
      isf = jseaflux(1)%iseaflux(j)
      iw  = seaflux(isf)%iw         ! global index
      iws = seaflux(isf)%iwls       ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
         iws = itabg_ws(iws)%iws_myrank
      endif

      arf_sea = seaflux(isf)%arf_sfc

      seaflux(isf)%rlong          = arf_sea * rlong(iw)
      seaflux(isf)%rshort         = arf_sea * rshort(iw)
      seaflux(isf)%rshort_diffuse = arf_sea * rshort_diffuse(iw)

   enddo
!$omp end parallel do

! Do parallel send of atm radiative fluxes to sea

   if (iparallel == 1) then  
      mrl = 1
      call mpi_send_wsf('R',mrl)
   endif

! If running leaf3, loop over LANDFLUX cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to land cells

!$omp parallel do private (ilf,iw,iwl,arf_land)
   do j = 1,jlandflux(1)%jend(1)
      ilf = jlandflux(1)%ilandflux(j)
      iw  = landflux(ilf)%iw         ! global index
      iwl = landflux(ilf)%iwls       ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iw  = itabg_w(iw)%iw_myrank
         iwl = itabg_wl(iwl)%iwl_myrank
      endif

      arf_land = landflux(ilf)%arf_sfc

      landflux(ilf)%rlong          = arf_land * rlong(iw)
      landflux(ilf)%rshort         = arf_land * rshort(iw)
      landflux(ilf)%rshort_diffuse = arf_land * rshort_diffuse(iw)

   enddo   
!$omp end parallel do

! Do parallel recv of atm radiative fluxes to sea

   if (iparallel == 1) then
      mrl = 1
      call mpi_recv_wsf('R',mrl)
   endif

! Do parallel send of atm radiative fluxes to land

   if (iparallel == 1) then  
      mrl = 1
      call mpi_send_wlf('R',mrl)
   endif

! If running leaf3, loop over SEAFLUX cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to sea cells

!$omp parallel do private (isf,iws)
   do j = 1,jseaflux(2)%jend(1)
      isf = jseaflux(2)%iseaflux(j)
      iws = seaflux(isf)%iwls       ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iws = itabg_ws(iws)%iws_myrank
      endif

      sea%rlong(iws)          = sea%rlong(iws)          + seaflux(isf)%rlong
      sea%rshort(iws)         = sea%rshort(iws)         + seaflux(isf)%rshort
      sea%rshort_diffuse(iws) = sea%rshort_diffuse(iws)  &
                              + seaflux(isf)%rshort_diffuse

   enddo
!$omp end parallel do

! Do parallel recv of atm radiative fluxes to land

   if (iparallel == 1) then  
      mrl = 1
      call mpi_recv_wlf('R',mrl)
   endif

! If running leaf3, loop over LANDFLUX cells to transfer downward surface
! shortwave and longwave fluxes from IW atmospheric column to land cells

!$omp parallel do private (ilf,iwl,arf_land)
   do j = 1,jlandflux(2)%jend(1)
      ilf = jlandflux(2)%ilandflux(j)
      iwl = landflux(ilf)%iwls       ! global index

! If run is parallel, get local rank indices

      if (iparallel == 1) then
         iwl = itabg_wl(iwl)%iwl_myrank
      endif

      arf_land = landflux(ilf)%arf_sfc

      land%rlong         (iwl) = land%rlong         (iwl) + landflux(ilf)%rlong
      land%rshort        (iwl) = land%rshort        (iwl) + landflux(ilf)%rshort
      land%rshort_diffuse(iwl) = land%rshort_diffuse(iwl)  &
                               + landflux(ilf)%rshort_diffuse
   enddo
!$omp end parallel do

! If running leaf3, loop over all LAND cells to compute radiative fluxes 
! for all cell components, given that rshort and rlong are now updated.

!$omp parallel do
   do iwl = 2,mwl

! If current IWL land cell is not prognosed on this rank, skip to next cell

      if (iparallel == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

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
            land%wnz           (      iwl)                                   )

      endif

   enddo
!$omp end parallel do

! Set the pointer for the first ED site

   ed_site => first_site

   do while(associated(ed_site))

! Loop over subgrid-scale patches

      ed_patch => ed_site%oldest_patch

      do while(associated(ed_patch))
         call scale_ed_radiation(ed_patch)
         ed_patch => ed_patch%younger
      enddo
      
      ed_site => ed_site%next_site

   enddo

endif

! Apply radiation tendencies in FTHRD to THILT

call psub()
!----------------------------------------------------------------------
mrl = mrl_begl(istp)
if (mrl > 0) then
!$omp parallel do private(iw,k)
do j = 1,jtab_w(12)%jend(mrl); iw = jtab_w(12)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   do k = lpw(iw),mza-1
      thilt(k,iw) = thilt(k,iw) + rho(k,iw) * fthrd(k,iw)
   enddo

enddo      
!$omp end parallel do
endif
call rsub('Wc',12)

return
end subroutine radiate

!============================================================================

subroutine sunloc()

use misc_coms,   only: io6, imonth1, idate1, iyear1, itime1, time_istp8,  &
                       iswrtyp, ilwrtyp

use consts_coms, only: pi2, pio180

use mem_radiate, only: jday, solfac, sunx, suny, sunz

use mem_mclat,   only: mclat_spline
use mem_harr,    only: nsolb, solar0, solar1

implicit none

integer :: outyear  ! current simulation year
integer :: outmonth ! current simulation month
integer :: outdate  ! current simulation date
integer :: outhour  ! current simulation hour/min/sec (6 digits)
integer :: ihour2   ! current simulation hour
integer :: imin2    ! current simulation min
integer :: isec2    ! current simulation sec
integer :: is       ! counter over solar bands in Harrington radiation

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

! Check whether Harrington or Fu-Liou shortwave or longwave radiation is used

if (iswrtyp >= 3 .or. ilwrtyp >= 3) then

! Adjust solar fluxes at top of atmosphere for current Earth-Sun distance
! for Harrington shortwave radiation

   do is = 1,nsolb
      solar1(is) = solar0(is) * solfac
   enddo

! Interpolate Mclatchy soundings between summer and winter values, and prepare
! spline coefficients for interpolation by latitude.

   call mclat_spline(jday)

endif

return
end subroutine sunloc
