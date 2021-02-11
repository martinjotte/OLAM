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
Module mem_sea

   implicit none

   integer :: nzsea = 0 ! Max number of vertical levels in sea/ocean model

   integer :: nsea      ! Total # of sea cell W pts in model global domain
   integer :: msea      ! Total # of sea cell W pts in model parallel sub-domain

   integer :: onsea     ! Offset # of sea cell W pts in model global domain
   integer :: omsea     ! Offset # of sea cell W pts in model parallel sub-domain
                        ! (isea + omsea = iwsfc)

! SEA GRID TABLES

   Type itab_sea_vars
      integer :: iwglobe = 1
   End type

   type (itab_sea_vars), allocatable, target :: itab_sea(:)

! SEA MODEL VARIABLES

   Type sea_vars

! Canopy to atmosphere turbulent flux quantities

      real, allocatable :: sea_ustar  (:) ! friction velocity [m/s]
      real, allocatable :: ice_ustar  (:) ! friction velocity [m/s]

      real, allocatable :: sea_vkmsfc (:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: ice_vkmsfc (:) ! surface drag coefficient [kg/(m s)]

      real, allocatable :: sea_sfluxt (:)
      real, allocatable :: ice_sfluxt (:)

      real, allocatable :: sea_sfluxr (:)
      real, allocatable :: ice_sfluxr (:)

      real, allocatable :: sea_sfluxc (:)
      real, allocatable :: ice_sfluxc (:)

      real, allocatable :: sea_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: ice_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]

      real, allocatable :: sea_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: ice_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]

      real, allocatable :: sea_sxfer_c(:) ! can_air-to-atm CO2 xfer [ppm/m^2]
      real, allocatable :: ice_sxfer_c(:) ! can_air-to-atm CO2 xfer [ppm/m^2]

      real, allocatable :: sea_ggaer  (:) ! surface aerodynamic conductance over water [m/s]
      real, allocatable :: ice_ggaer  (:) ! surface aerodynamic conductance over ice [m/s]

      real, allocatable :: sea_wthv   (:) ! surface buoyancy flux over water [K m/s]
      real, allocatable :: ice_wthv   (:) ! surface buoyancy flux over ice [K m/s]

! Radiative flux quantities

      real, allocatable :: sea_albedo    (:) ! water s/w albedo [0-1]
      real, allocatable :: ice_albedo    (:) ! seaice s/w albedo [0-1]

      real, allocatable :: sea_rlongup   (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: ice_rlongup   (:) ! upward can-top l/w flux [W/m^2]

      real, allocatable :: ice_net_rlong (:)
      real, allocatable :: ice_net_rshort(:)

! Canopy and surface quantities:

      real, allocatable :: sea_cantemp   (:) ! "canopy" air temperature [K]
      real, allocatable :: ice_cantemp   (:) ! "canopy" air temperature [K]

      real, allocatable :: sea_canrrv    (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: ice_canrrv    (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]

      integer, allocatable :: nlev_seaice(:) ! number of seaice layers in current cell

      real, allocatable :: seatp         (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf         (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc         (:) ! current sea temperature [K]
      real, allocatable :: seaicep       (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef       (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec       (:) ! seaice fraction [0-1]

      real, allocatable :: surface_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: sea_sfc_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]
      real, allocatable :: ice_sfc_srrv  (:) ! sea surface sat vapor mixing ratio [kg_vap/kg_dryair]

      real, allocatable :: sea_rough     (:) ! water surface roughness height [m]
      real, allocatable :: ice_rough     (:) ! water surface roughness height [m]

      real, allocatable :: seaice_energy(:,:)
      real, allocatable :: seaice_tempk (:,:)

   End Type sea_vars

   type (sea_vars) :: sea

Contains

!=========================================================================

   subroutine alloc_sea(msea)

     use misc_coms, only: rinit
     use sea_coms,  only: nzi
     use mem_co2,   only: co2flag

     implicit none

     integer, intent(in) :: msea

!    Allocate and initialize sea arrays

     allocate (sea%sea_ustar     (msea)) ; sea%sea_ustar      = rinit
     allocate (sea%ice_ustar     (msea)) ; sea%ice_ustar      = rinit

     allocate (sea%sea_vkmsfc    (msea)) ; sea%sea_vkmsfc     = rinit
     allocate (sea%ice_vkmsfc    (msea)) ; sea%ice_vkmsfc     = rinit

     allocate (sea%sea_sfluxt    (msea)) ; sea%sea_sfluxt     = rinit
     allocate (sea%ice_sfluxt    (msea)) ; sea%ice_sfluxt     = rinit

     allocate (sea%sea_sfluxr    (msea)) ; sea%sea_sfluxr     = rinit
     allocate (sea%ice_sfluxr    (msea)) ; sea%ice_sfluxr     = rinit

     allocate (sea%sea_sxfer_t   (msea)) ; sea%sea_sxfer_t    = 0.0
     allocate (sea%ice_sxfer_t   (msea)) ; sea%ice_sxfer_t    = 0.0

     allocate (sea%sea_sxfer_r   (msea)) ; sea%sea_sxfer_r    = 0.0
     allocate (sea%ice_sxfer_r   (msea)) ; sea%ice_sxfer_r    = 0.0

     if (co2flag /= 0) then
        allocate (sea%sea_sfluxc    (msea)) ; sea%sea_sfluxc     = rinit
        allocate (sea%ice_sfluxc    (msea)) ; sea%ice_sfluxc     = rinit

        allocate (sea%sea_sxfer_c   (msea)) ; sea%sea_sxfer_c    = 0.0
        allocate (sea%ice_sxfer_c   (msea)) ; sea%ice_sxfer_c    = 0.0
     endif

     allocate (sea%sea_ggaer     (msea)) ; sea%sea_ggaer      = rinit
     allocate (sea%ice_ggaer     (msea)) ; sea%ice_ggaer      = rinit

     allocate (sea%sea_wthv      (msea)) ; sea%sea_wthv       = rinit
     allocate (sea%ice_wthv      (msea)) ; sea%ice_wthv       = rinit

     allocate (sea%sea_albedo    (msea)) ; sea%sea_albedo     = 0.0
     allocate (sea%ice_albedo    (msea)) ; sea%ice_albedo     = 0.0

     allocate (sea%sea_rlongup   (msea)) ; sea%sea_rlongup    = 0.0
     allocate (sea%ice_rlongup   (msea)) ; sea%ice_rlongup    = 0.0

     allocate (sea%ice_net_rlong (msea)) ; sea%ice_net_rlong  = 0.0
     allocate (sea%ice_net_rshort(msea)) ; sea%ice_net_rshort = 0.0

     allocate (sea%sea_cantemp   (msea)) ; sea%sea_cantemp    = rinit
     allocate (sea%ice_cantemp   (msea)) ; sea%ice_cantemp    = rinit

     allocate (sea%sea_canrrv    (msea)) ; sea%sea_canrrv     = rinit
     allocate (sea%ice_canrrv    (msea)) ; sea%ice_canrrv     = rinit

     allocate (sea%nlev_seaice   (msea)) ; sea%nlev_seaice    = 0

     allocate (sea%seatp         (msea)) ; sea%seatp          = rinit
     allocate (sea%seatf         (msea)) ; sea%seatf          = rinit
     allocate (sea%seatc         (msea)) ; sea%seatc          = rinit
     allocate (sea%seaicep       (msea)) ; sea%seaicep        = rinit
     allocate (sea%seaicef       (msea)) ; sea%seaicef        = rinit
     allocate (sea%seaicec       (msea)) ; sea%seaicec        = rinit

     allocate (sea%surface_srrv  (msea)) ; sea%surface_srrv   = rinit
     allocate (sea%sea_sfc_srrv  (msea)) ; sea%sea_sfc_srrv   = rinit
     allocate (sea%ice_sfc_srrv  (msea)) ; sea%ice_sfc_srrv   = rinit

     allocate (sea%sea_rough     (msea)) ; sea%sea_rough      = rinit
     allocate (sea%ice_rough     (msea)) ; sea%ice_rough      = rinit

     allocate (sea%seaice_energy(nzi,msea)) ; sea%seaice_energy = rinit
     allocate (sea%seaice_tempk (nzi,msea)) ; sea%seaice_tempk  = rinit

   end subroutine alloc_sea

!=========================================================================

   subroutine filltab_sea()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(sea%sea_ustar)) call increment_vtable('SEA%SEA_USTAR', 'SW', rvar1=sea%sea_ustar)
     if (allocated(sea%ice_ustar)) call increment_vtable('SEA%ICE_USTAR', 'SW', rvar1=sea%ice_ustar)

     if (allocated(sea%sea_vkmsfc)) call increment_vtable('SEA%SEA_VKMSFC', 'SW', rvar1=sea%sea_vkmsfc)
     if (allocated(sea%ice_vkmsfc)) call increment_vtable('SEA%ICE_VKMSFC', 'SW', rvar1=sea%ice_vkmsfc)

     if (allocated(sea%sea_sfluxt)) call increment_vtable('SEA%SEA_SFLUXT', 'SW', rvar1=sea%sea_sfluxt)
     if (allocated(sea%ice_sfluxt)) call increment_vtable('SEA%ICE_SFLUXT', 'SW', rvar1=sea%ice_sfluxt)

     if (allocated(sea%sea_sfluxr)) call increment_vtable('SEA%SEA_SFLUXR', 'SW', rvar1=sea%sea_sfluxr)
     if (allocated(sea%ice_sfluxr)) call increment_vtable('SEA%ICE_SFLUXR', 'SW', rvar1=sea%ice_sfluxr)

     if (allocated(sea%sea_sfluxc)) call increment_vtable('SEA%SEA_SFLUXC', 'SW', rvar1=sea%sea_sfluxc)
     if (allocated(sea%ice_sfluxc)) call increment_vtable('SEA%ICE_SFLUXC', 'SW', rvar1=sea%ice_sfluxc)

     if (allocated(sea%sea_sxfer_t)) call increment_vtable('SEA%SEA_SXFER_T', 'SW', rvar1=sea%sea_sxfer_t)
     if (allocated(sea%ice_sxfer_t)) call increment_vtable('SEA%ICE_SXFER_T', 'SW', rvar1=sea%ice_sxfer_t)

     if (allocated(sea%sea_sxfer_r)) call increment_vtable('SEA%SEA_SXFER_R', 'SW', rvar1=sea%sea_sxfer_r)
     if (allocated(sea%ice_sxfer_r)) call increment_vtable('SEA%ICE_SXFER_R', 'SW', rvar1=sea%ice_sxfer_r)

     if (allocated(sea%sea_sxfer_c)) call increment_vtable('SEA%SEA_SXFER_C', 'SW', rvar1=sea%sea_sxfer_c)
     if (allocated(sea%ice_sxfer_c)) call increment_vtable('SEA%ICE_SXFER_C', 'SW', rvar1=sea%ice_sxfer_c)

     if (allocated(sea%sea_wthv)) call increment_vtable('SEA%SEA_WTHV', 'SW', rvar1=sea%sea_wthv)
     if (allocated(sea%ice_wthv)) call increment_vtable('SEA%ICE_WTHV', 'SW', rvar1=sea%ice_wthv)

     if (allocated(sea%sea_albedo)) call increment_vtable('SEA%SEA_ALBEDO', 'SW', rvar1=sea%sea_albedo)
     if (allocated(sea%ice_albedo)) call increment_vtable('SEA%ICE_ALBEDO', 'SW', rvar1=sea%ice_albedo)

     if (allocated(sea%sea_rlongup)) call increment_vtable('SEA%SEA_RLONGUP', 'SW', rvar1=sea%sea_rlongup)
     if (allocated(sea%ice_rlongup)) call increment_vtable('SEA%ICE_RLONGUP', 'SW', rvar1=sea%ice_rlongup)

      if (allocated(sea%ice_net_rlong)) call increment_vtable('SEA%ICE_NET_RLONG', 'SW', rvar1=sea%ice_net_rlong)

      if (allocated(sea%ice_net_rshort)) call increment_vtable('SEA%ICE_NET_RSHORT', 'SW', rvar1=sea%ice_net_rshort)

     if (allocated(sea%seatc)) call increment_vtable('SEA%SEATC', 'SW', rvar1=sea%seatc)

     if (allocated(sea%sea_cantemp)) call increment_vtable('SEA%SEA_CANTEMP', 'SW', rvar1=sea%sea_cantemp)
     if (allocated(sea%ice_cantemp)) call increment_vtable('SEA%ICE_CANTEMP', 'SW', rvar1=sea%ice_cantemp)

     if (allocated(sea%sea_canrrv)) call increment_vtable('SEA%SEA_CANRRV', 'SW', rvar1=sea%sea_canrrv)
     if (allocated(sea%ice_canrrv)) call increment_vtable('SEA%ICE_CANRRV', 'SW', rvar1=sea%ice_canrrv)

     if (allocated(sea%nlev_seaice)) call increment_vtable('SEA%NLEV_SEAICE', 'SW', ivar1=sea%nlev_seaice)

     if (allocated(sea%seaicec)) call increment_vtable('SEA%SEAICEC', 'SW', rvar1=sea%seaicec)

     if (allocated(sea%surface_srrv)) call increment_vtable('SEA%SURFACE_SRRV', 'SW', rvar1=sea%surface_srrv)
     if (allocated(sea%sea_sfc_srrv)) call increment_vtable('SEA%SEA_SFC_SRRV', 'SW', rvar1=sea%sea_sfc_srrv)
     if (allocated(sea%ice_sfc_srrv)) call increment_vtable('SEA%ICE_SFC_SRRV', 'SW', rvar1=sea%ice_sfc_srrv)

     if (allocated(sea%sea_rough)) call increment_vtable('SEA%SEA_ROUGH', 'SW', rvar1=sea%sea_rough)
     if (allocated(sea%ice_rough)) call increment_vtable('SEA%ICE_ROUGH', 'SW', rvar1=sea%ice_rough)

     if (allocated(sea%seaice_energy)) call increment_vtable('SEA%SEAICE_ENERGY', 'SW', rvar2=sea%seaice_energy)

     if (allocated(sea%seaice_tempk)) call increment_vtable('SEA%SEAICE_TEMPK', 'SW', rvar2=sea%seaice_tempk)

   end subroutine filltab_sea

End Module mem_sea
