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

  use sea_coms,  only: mws, nzi, iupdsst, s1900_sst, isstfile, nsstfiles,  &
                       iupdseaice, s1900_seaice, iseaicefile, nseaicefiles

  use mem_sea,   only: sea, itab_ws
  use misc_coms, only: io6, time8, s1900_sim, isubdomain
  use mem_para,  only: myrank
!$use omp_lib

  implicit none

! Local variables

  integer :: iws      ! sea cell loop counter
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

!$omp parallel do
  do iws = 2, mws

! Update SEATC and seaice fraction

     sea%seatc(iws) = sea%seatp(iws) + &
                      timefac_sst * (sea%seatf(iws) - sea%seatp(iws))

     sea%seaicec(iws) = sea%seaicep(iws) + &
                        timefac_seaice * (sea%seaicef(iws) - sea%seaicep(iws))

! Update SEA fields

     call seacell(iws                 ,  &
                  sea%rhos       (iws),  &
                  sea%sea_ustar  (iws),  &
                  sea%sea_ggaer  (iws),  &
                  sea%sea_sxfer_t(iws),  &
                  sea%sea_sxfer_r(iws),  &
                  sea%can_depth  (iws),  &
                  sea%seatc      (iws),  &
                  sea%sea_cantemp(iws),  &
                  sea%sea_canshv (iws),  &
                  sea%sea_sfc_ssh(iws),  &
                  sea%sea_rough  (iws)   )

! Update sea ice based on seaice fraction

     call prep_seaice(sea%seatc              (iws), &      
                      sea%seaicec            (iws), &
                      sea%sea_cantemp        (iws), &
                      sea%ice_cantemp        (iws), &
                      sea%seaice_energy(1:nzi,iws), &
                      sea%seaice_tempk (1:nzi,iws), &
                      sea%nlev_seaice        (iws), &
                      sea%ice_albedo         (iws), &
                      sea%ice_rlongup        (iws), &
                      sea%ice_net_rshort     (iws), &
                      sea%ice_net_rlong      (iws), &
                      sea%rshort             (iws), &
                      sea%rlong              (iws), &
                      sea%ice_rough          (iws), &
                      sea%sea_canshv         (iws), &
                      sea%ice_canshv         (iws), &
                      sea%sea_ustar          (iws), &
                      sea%ice_ustar          (iws), &
                      sea%sea_ggaer          (iws), &
                      sea%ice_ggaer          (iws), &
                      sea%ice_sxfer_t        (iws), &
                      sea%ice_sxfer_r        (iws)  )

! If ice exists, compute seaice canopy fluxes

     if (sea%nlev_seaice(iws) > 0) then

        call seaice(sea%seaice_energy(1:nzi,iws), &
                    sea%seaice_tempk (1:nzi,iws), &
                    sea%nlev_seaice        (iws), &
                    sea%ice_net_rshort     (iws), &
                    sea%ice_net_rlong      (iws), &
                    sea%rhos               (iws), &
                    sea%ice_ustar          (iws), &
                    sea%ice_ggaer          (iws), &
                    sea%can_depth          (iws), &
                    sea%ice_cantemp        (iws), &
                    sea%ice_canshv         (iws), &
                    sea%ice_sfc_ssh        (iws), &
                    sea%ice_sxfer_t        (iws), &
                    sea%ice_sxfer_r        (iws)  )

        sea%rough      (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_rough  (iws) + &
                                      sea%seaicec(iws)  * sea%ice_rough  (iws)

        sea%cantemp    (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_cantemp(iws) + &
                                      sea%seaicec(iws)  * sea%ice_cantemp(iws)

        sea%canshv     (iws) = (1.0 - sea%seaicec(iws)) * sea%sea_canshv (iws) + &
                                      sea%seaicec(iws)  * sea%ice_canshv (iws)

        sea%surface_ssh(iws) = (1.0 - sea%seaicec(iws)) * sea%sea_sfc_ssh(iws) + &
                                      sea%seaicec(iws)  * sea%ice_sfc_ssh(iws)
     else

        sea%rough      (iws) = sea%sea_rough  (iws)
        sea%cantemp    (iws) = sea%sea_cantemp (iws)
        sea%canshv     (iws) = sea%sea_canshv  (iws)
        sea%surface_ssh(iws) = sea%sea_sfc_ssh(iws)

     endif

! Zero out SEA%SXFER_T(IWS) and SEA%SXFER_R(IWS) now that they have 
! been applied to the canopy

   sea%sxfer_t(iws) = 0.
   sea%sxfer_r(iws) = 0.
   sea%sxfer_c(iws) = 0.

   sea%sea_sxfer_t(iws) = 0.
   sea%sea_sxfer_r(iws) = 0.

   sea%ice_sxfer_t(iws) = 0.
   sea%ice_sxfer_r(iws) = 0.
                   
enddo
!$omp end parallel do

return
end subroutine seacells

!===============================================================================

subroutine seacell( iws, rhos, ustar, ggaer, sxfer_t, sxfer_r, can_depth, &
                    seatc, cantemp, canshv, surface_ssh, rough           )

  use sea_coms,    only: dt_sea
  use consts_coms, only: cp, grav
  use misc_coms,   only: io6

  implicit none

  integer, intent(in)    :: iws         ! current sea cell index
  real,    intent(in)    :: rhos        ! air density [kg/m^3]
  real,    intent(in)    :: ustar       ! friction velocity [m/s]
  real,    intent(in)    :: ggaer       ! sfc. aerodynamic conductance [m/s]
  real,    intent(in)    :: sxfer_t     ! can_air to atm heat xfer this step [kg_air K/m^2]
  real,    intent(in)    :: sxfer_r     ! can_air to atm vapor xfer this step (kg_vap/m^2]
  real,    intent(in)    :: can_depth   ! "canopy" depth for heat and vap capacity [m]
  real,    intent(in)    :: seatc       ! current sea temp (obs time) [K]
  real,    intent(inout) :: cantemp     ! "canopy" air temp [K]
  real,    intent(inout) :: canshv      ! "canopy" air vapor spec hum [kg_vap/kg_air]
  real,    intent(out)   :: surface_ssh ! sea surface sat spec hum [kg_vap/kg_air] 
  real,    intent(out)   :: rough       ! sea cell roughess height [m] 
  
! Local parameter

  real, parameter :: z0fac_water = .016 / grav  ! factor for Charnok roughness height
  real, parameter :: ozo = 1.59e-5              ! base roughness height in HWRF

! Local variables

  real :: rdi     ! canopy conductance [m/s]
  real :: hxfergc ! heat xfer from sea surface to can_air this step [J/m^2]
  real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]
  real :: wxfergc ! vapor xfer from sea surface to can_air this step [kg_vap/m^2]

  real :: zn1, zn2, zw
  real, external :: rhovsil  ! function to compute saturation vapor density

! Evaluate surface saturation specific humidity

  surface_ssh = rhovsil(seatc-273.15) / rhos

! Update temperature and vapor specific humidity of "canopy" from
! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity derived from Garratt (1992), or
! we can use the bare surface aerodynamic conductance ggaer computed
! from the surface layer similarity relationships

  rdi = .2 * ustar
! rdi = ggaer

  hxfergc = dt_sea * cp * rhos * rdi * (seatc - cantemp)
  wxfergc = dt_sea *      rhos * rdi * (surface_ssh - canshv)

  hxferca = cp * sxfer_t  ! sxfer_t and sxfer_r already incorporate dt_sea

  cantemp = cantemp + (hxfergc - hxferca) / (can_depth * rhos * cp)
  canshv  = canshv  + (wxfergc - sxfer_r) / (can_depth * rhos)             

! Evaluate sea roughness height

! Charnok (1955):
! rough = max(z0fac_water * ustar ** 2,.0001)  ! Charnok (1955)

! Davis et al. (2008) originally used in HWRF
! rough = 10. * exp(-10. / ustar ** .333333)
! rough = max(.125e-6, min(2.85e-3,rough))

! 2012 HWRF scheme; interpolates between the Charnok scheme at low wind
! and the Davis et al. curve fit at high wind speeds
  zw    = min( (ustar/1.06)**0.3, 1.0)
  zn1   = 0.011 * ustar * ustar /grav + ozo
  zn2   = 10. * exp(-9.5 * ustar**(-.3333333)) + 1.65e-6 / ustar
  rough = (1.0-zw) * zn1 + zw * zn2
  rough = min( rough, 2.85e-3)

end subroutine seacell
