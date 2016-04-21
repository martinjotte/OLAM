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

subroutine prep_seaice( sst, seaice, sea_cantemp, ice_cantemp,        &
                        seaice_energy, seaice_tempk, nlev_seaice,     &
                        ice_albedo, ice_rlongup, rshort_i, rlong_i,   &
                        rshort, rlong, rough, sea_canshv, ice_canshv, &
                        sea_ustar, ice_ustar, sea_ggaer, ice_ggaer,   &
                        sea_wthv, ice_wthv, sxfer_t, sxfer_r          )

  use consts_coms, only: cice, t00
  use sea_coms,    only: nzi, t00sea

  implicit none

  real,    intent(in)    :: sst
  real,    intent(in)    :: seaice
  real,    intent(in)    :: sea_cantemp
  real,    intent(out)   :: ice_cantemp
  real,    intent(out)   :: seaice_energy(nzi) ! seaice layer energy [J/kg]
  real,    intent(out)   :: seaice_tempk (nzi) ! seaice layer temperature [K]
  integer, intent(inout) :: nlev_seaice
  real,    intent(out)   :: ice_albedo
  real,    intent(out)   :: ice_rlongup
  real,    intent(out)   :: rshort_i
  real,    intent(out)   :: rlong_i
  real,    intent(in)    :: rshort
  real,    intent(in)    :: rlong
  real,    intent(out)   :: rough  ! sea cell roughess height [m] 
  real,    intent(in)    :: sea_canshv
  real,    intent(out)   :: ice_canshv
  real,    intent(in)    :: sea_ustar
  real,    intent(out)   :: ice_ustar
  real,    intent(in)    :: sea_ggaer
  real,    intent(out)   :: ice_ggaer
  real,    intent(in)    :: sea_wthv
  real,    intent(out)   :: ice_wthv
  real,    intent(out)   :: sxfer_t
  real,    intent(out)   :: sxfer_r
  integer                :: k
  
! First check whether sea ice is present

  if (seaice < 0.01) then

     ! If no sea ice present, nullify quantities and return

     nlev_seaice      = 0
     seaice_tempk (:) = 0.0
     seaice_energy(:) = 0.0
     ice_cantemp      = 0.0
     ice_canshv       = 0.0
     ice_albedo       = 0.0
     ice_rlongup      = 0.0
     rshort_i         = 0.0
     rlong_i          = 0.0
     rough            = 0.0
     ice_ustar        = 0.0
     ice_ggaer        = 0.0
     ice_wthv         = 0.0
     sxfer_t          = 0.0
     sxfer_r          = 0.0

     return

  else

     ! Always set bottom layer temperature equal to the SST

     seaice_tempk (1) = min(sst, t00sea)
     seaice_energy(1) = cice * (seaice_tempk(1) - t00sea)
     
     ! Ice roughness height

     rough = 5.0e-4 ! [m]; taken from Los Alamos sea ice model

  endif

! Seaice is present. If there was no seaice the previous timestep,
! initialize the ice temperature and energy based on the SST and air
! temperature, and initialize the ice canopy temperature, vapor, and
! ustar to the sea canopy values

  if (nlev_seaice == 0) then

     nlev_seaice = nzi
     ice_cantemp = sea_cantemp
     ice_canshv  = sea_canshv
     ice_ustar   = sea_ustar
     ice_ggaer   = sea_ggaer
     ice_wthv    = sea_wthv
     sxfer_t     = 0.0
     sxfer_r     = 0.0

     ! Initialize top layer to the canopy temperature

     seaice_tempk (nlev_seaice) = min(ice_cantemp, t00sea)
     seaice_energy(nlev_seaice) = cice * (seaice_tempk(nlev_seaice) - t00sea)

     ! Interpolate other layers between the top and bottom temperature

     do k = 2, nlev_seaice-1

        seaice_tempk(k) = real(k-1)/real(nlev_seaice-1) * seaice_tempk(nlev_seaice) &
                        + real(nlev_seaice-k)/real(nlev_seaice-1) * seaice_tempk(1)

        seaice_energy(k) = cice * (seaice_tempk(k) - t00sea)

     enddo

     ! Estimate ice albedo and outgoing longwave in case we have
     ! created seaice between calls to the radiation scheme

     call sfcrad_seaice_1( ice_rlongup,        &
                           ice_albedo,         &
                           nlev_seaice,        &
                           ice_cantemp,        &
                           seaice_tempk(1:nzi) )

     ! Estimate net seaice longwave and shortwave radiative fluxes in case
     ! seaice has been created between calls to the radiation scheme

     call sfcrad_seaice_2( rshort_i,    &
                           rlong_i,     &
                           nlev_seaice, &
                           rshort,      &
                           rlong,       &
                           ice_rlongup, &
                           ice_albedo   )

  endif

end subroutine prep_seaice
  




subroutine seaice( seaice_energy, seaice_tempk, nlev_seaice,      &
                   rshort_i, rlong_i, rhos, ustar, can_depth,     &
                   cantemp, canshv, surface_ssh, sxfer_t, sxfer_r )

  use consts_coms, only: alvi, cice, t00, cp, alli
  use sea_coms,    only: dt_sea, t00sea, nzi

  implicit none

  real,    intent(inout) :: seaice_energy(nzi) ! seaice layer energy [J/kg]
  real,    intent(inout) :: seaice_tempk (nzi) ! seaice layer temperature [K]

  integer, intent(in)    :: nlev_seaice ! number of seaice levels
  real,    intent(in)    :: rshort_i    ! s/w net rad flux to seaice [W/m^2]
  real,    intent(in)    :: rlong_i     ! l/w net rad flux to seaice [W/m^2]
  real,    intent(in)    :: rhos        ! air density [kg/m^3]
  real,    intent(in)    :: ustar       ! friction velocity [m/s]
  real,    intent(in)    :: can_depth   ! "canopy" depth for heat and vap capacity [m]
  real,    intent(inout) :: cantemp     ! ice "canopy" air temp [K]
  real,    intent(inout) :: canshv      ! ice "canopy" air vapor spec hum [kg_vap/kg_air]
  real,    intent(out)   :: surface_ssh ! ice surface sat spec hum [kg_vap/kg_air] 
  real,    intent(in)    :: sxfer_t     ! can_air to atm heat xfer this step [kg_air K/m^2]
  real,    intent(in)    :: sxfer_r     ! can_air to atm vapor xfer this step (kg_vap/m^2]
 
! Local parameters

  real, parameter :: iceden   = 900.0   ! Density of seaice layers [kg/m^3]
  real, parameter :: icethick = 0.5     ! Thickness of seaice layers [m]
  real, parameter :: icek     = 2.0     ! Thermal conductivity of seaice [W/(K m)]

  real, parameter :: icemass  = icethick * iceden ! mass of unit area of ice [kg/m^2]
  real, parameter :: icer     = icethick / icek   ! ice thermal resistivity [K m^2/W]

! Local variables

  real :: energy_per_m2(nzi) ! seaice energy [J/m^2]
  real :: hxfers(nzi+1)      ! seaice heat xfer [J/m2] 

  real :: fracliq ! fraction of water in liquid phase
  real :: rdi     ! canopy conductance [m/s]
  real :: hxferic ! seaice-to-can_air heat xfer this step [J/m^2]
  real :: wxferic ! seaice-to-can_air vap xfer this step [kg_vap/m^2]
  real :: hxferca ! heat xfer from can_air to atm this step [J/m^2]

  real, external :: rhovsil ! function to compute saturation vapor density
  integer :: k

  if (nlev_seaice == 0) return

! Compute sfcwater energy per m^2   

  do k = 1, nlev_seaice
     energy_per_m2(k) = seaice_energy(k) * icemass
  enddo

! Evaluate surface saturation specific humidity

  surface_ssh = rhovsil(seaice_tempk(nlev_seaice)-t00) / rhos

! Update temperature and vapor specific humidity of "canopy" from
! divergence of xfers with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity derived from Garratt (1992)

  rdi = .2 * ustar

  hxferic = dt_sea * cp * rhos * rdi * (seaice_tempk(nlev_seaice) - cantemp)
  wxferic = dt_sea *      rhos * rdi * (surface_ssh - canshv)

  hxferca = cp * sxfer_t  ! sxfer_t and sxfer_r already incorporate dt_sea

  cantemp = cantemp + (hxferic - hxferca) / (can_depth * rhos * cp)
  canshv  = canshv  + (wxferic - sxfer_r) / (can_depth * rhos)             

! Zero out sfcwater internal heat transfer array at top and bottom surfaces.
! Energy transfer at top is applied separately.
  
  hxfers(1)             = 0.0
  hxfers(nlev_seaice+1) = 0.0

! Compute internal seaice energy xfer if at least two layers exist [J/m2]
   
  do k = 2, nlev_seaice
     hxfers(k) = dt_sea * (seaice_tempk(k-1) - seaice_tempk(k)) / icer
  enddo

! Add contributions to seaice energy from internal transfers of heat

  do k = 2, nlev_seaice
     energy_per_m2(k) = energy_per_m2(k) + hxfers(k) - hxfers(k+1)
  enddo

! Apply long- and short-wave radiative transfer and seaice-to-canopy 
! sensible and latent heat transfer to top seaice layer

  energy_per_m2(nlev_seaice) = energy_per_m2(nlev_seaice)  &
       + dt_sea * (rshort_i + rlong_i) - hxferic - wxferic * alvi

! Compute new seaice_energy

  do k = 2, nlev_seaice

     seaice_energy(k) = energy_per_m2(k) / icemass

     ! Bound seaice_energy so that the fraction of liquid water in each layer
     ! is never larger than 50%. For now, we trust that ice should exist when
     ! the input seaice data indicates that ice is present, and this keeps the
     ! seaice temperature at or below its freezing point

     seaice_energy(k) = min( seaice_energy(k), 0.5 * alli )

  enddo

! Update seaice temperature
   
  do k = 2, nlev_seaice
     call qtk_sea( seaice_energy(k), seaice_tempk(k), fracliq )
  enddo

end subroutine seaice



subroutine sfcrad_seaice_1( rlongup, albedo, nlev_seaice, cantemp, seaice_tempk )

  use sea_coms,    only: nzi, emi
  use consts_coms, only: stefan
  implicit none

  real,    intent(out) :: rlongup   ! seaice outgoing L/W radiation[W/m^2]
  real,    intent(out) :: albedo    ! seaice albedo                [0-1]

  integer, intent(in)  :: nlev_seaice       ! number of seaice levels
  real,    intent(in)  :: cantemp           ! seaice canopy temperature [K]
  real,    intent(in)  :: seaice_tempk(nzi) ! seaice layer temperatures [K]

  if (nlev_seaice > 0) then

     rlongup = emi * stefan * seaice_tempk(nlev_seaice) ** 4

     ! In the Los Alamos sea ice model, visible albedo is 0.78 and 
     ! NIR albedo is 0.36 for temperatures below -1C.  Here, we assume
     ! equal parts visible and NIR, yielding an albedo of 0.57.  This
     ! albedo decreases to 0.5 as the temperature increases to 0C.

     if (cantemp < 272.15) then
        albedo = 0.72  ! Dave modification
     !  albedo = 0.57
     else
        albedo = 0.72 - 0.22 * ( min( 273.15, cantemp) - 272.15 )
     endif

  else

     rlongup = 0.0
     albedo  = 0.0

  end if

end subroutine sfcrad_seaice_1
     


subroutine sfcrad_seaice_2( rshort_i, rlong_i, nlev_seaice, &
                            rshort, rlong, rlongup, albedo  )

  use sea_coms, only: nzi, emi
  implicit none

  real,    intent(out) :: rshort_i    ! net S/W radiation at ice surface [W/m^2]
  real,    intent(out) :: rlong_i     ! net L/W radiation at ice surface [W/m^2]

  integer, intent(in)  :: nlev_seaice ! number of seaice levls
  real,    intent(in)  :: rshort      ! incoming S/W radiation       [W/m^2]
  real,    intent(in)  :: rlong       ! incoming L/W radiation       [W/m^2]
  real,    intent(in)  :: rlongup     ! seaice outgoing L/W radiation[W/m^2]
  real,    intent(in)  :: albedo      ! seaice albedo                [0-1]

  if (nlev_seaice > 0) then

     rshort_i = (1.0 - albedo) * rshort
     rlong_i  = rlong * emi - rlongup

  else

     rshort_i = 0.0
     rlong_i  = 0.0

  end if

end subroutine sfcrad_seaice_2
