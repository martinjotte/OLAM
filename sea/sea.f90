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
subroutine seacells()

  use sea_coms,    only: nzi, iupdsst, s1900_sst, isstfile, nsstfiles,  &
                         iupdseaice, s1900_seaice, iseaicefile, nseaicefiles
  use mem_ijtabs,  only: itabg_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_sea,     only: sea, msea, omsea
  use misc_coms,   only: io6, time8, s1900_sim, isubdomain
  use mem_para,    only: myrank
  use mem_basic,   only: rho, press, theta, tair, vxe, vye, vze, rr_v
  use mem_micro,   only: rr_c
  use consts_coms, only: grav

  implicit none

! Local variables

  integer :: isea      ! sea cell loop counter
  integer :: iw, kw, iwsfc, j
  real    :: timefac_sst   ! fraction of elapsed time from past to future SST obs
  real    :: timefac_seaice   ! fraction of elapsed time from past to future SEA ICE obs
  real    :: psfc, vels

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
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

! Zero runoff for sea cells

     sfcg%runoff(iwsfc) = 0.

! Update SEATC and seaice fraction

     sea%seatc(isea) = sea%seatp(isea) &
                     + timefac_sst * (sea%seatf(isea) - sea%seatp(isea))

     sea%seaicec(isea) = sea%seaicep(isea) &
                       + timefac_seaice * (sea%seaicef(isea) - sea%seaicep(isea))

! Update SEA fields

     call seacell(isea                  ,  &
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

subroutine seacell( isea, rhos, ustar, sxfer_t, sxfer_r, can_depth, &
                    seatc, cantemp, canrrv, surface_srrv, rough     )

  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav
  use misc_coms,   only: io6
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in)    :: isea         ! current sea cell index
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
  real :: hxfergc ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxfergc ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]

  real :: zn1, zn2, zw, usti

! Evaluate surface saturation specific humidity

  surface_srrv = rhovsl(seatc-273.15) / rhos

! Update temperature and vapor specific humidity of "canopy" from
! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  hxfergc = dt_sea * cp * rhos * rdi * (seatc        - cantemp)
  wxfergc = dt_sea *      rhos * rdi * (surface_srrv - canrrv)

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

end subroutine seacell
