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

   use mem_mksfc, only: itab_wls_vars

   implicit none

   private :: itab_wls_vars

! SEA SURFACE GRID TABLES

   type (itab_wls_vars), allocatable, target :: itab_ws(:)

! DERIVED TYPES TO HOLD GLOBAL GRID INDICES FOR A PARALLEL RUN

   Type itabg_ws_vars
      integer :: iws_myrank = -1
      integer :: irank      = -1
   End Type itabg_ws_vars

   type (itabg_ws_vars), allocatable, target :: itabg_ws(:)

! SEA SURFACE MODEL VARIABLES

   Type sea_vars

! W-point grid variables stored on SEAFILE:

      real, allocatable :: area (:) ! cell surface area [m^2]
      real, allocatable :: glatw(:) ! latitude of sea cell W points
      real, allocatable :: glonw(:) ! longitude of sea cell W points
      real, allocatable :: xew  (:) ! earth x coord of sea cell W points
      real, allocatable :: yew  (:) ! earth y coord of sea cell W points
      real, allocatable :: zew  (:) ! earth z coord of sea cell W points
      real, allocatable :: topw (:) ! topographic height of sea cell W points

      integer, allocatable :: leaf_class (:) ! sea cell leaf class

! Atmospheric near-surface properties

      real, allocatable :: rhos   (:) ! air density [kg_air/m^3]
      real, allocatable :: vels   (:) ! wind speed [m/s]
      real, allocatable :: prss   (:) ! air pressure [Pa]
      real, allocatable :: airtemp(:) ! air temperature [K]
      real, allocatable :: airshv (:) ! air specific humidity[kg_vap/kg_air]

! Canopy to atmosphere turbulent flux quantities

      real, allocatable :: ustar      (:) ! friction velocity [m/s]
      real, allocatable :: sea_ustar  (:) ! friction velocity [m/s]
      real, allocatable :: ice_ustar  (:) ! friction velocity [m/s]

      real, allocatable :: vkmsfc     (:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: sea_vkmsfc (:) ! surface drag coefficient [kg/(m s)]
      real, allocatable :: ice_vkmsfc (:) ! surface drag coefficient [kg/(m s)]

      real, allocatable :: sfluxt     (:)
      real, allocatable :: sea_sfluxt (:)
      real, allocatable :: ice_sfluxt (:)

      real, allocatable :: sfluxr     (:)
      real, allocatable :: sea_sfluxr (:)
      real, allocatable :: ice_sfluxr (:)

      real, allocatable :: sfluxc     (:)

      real, allocatable :: sxfer_t    (:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: sea_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: ice_sxfer_t(:) ! can_air-to-atm heat xfer [kg_air K/m^2]

      real, allocatable :: sxfer_r    (:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: sea_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: ice_sxfer_r(:) ! can_air-to-atm vapor xfer [kg_vap/m^2]

      real, allocatable :: sxfer_c    (:) ! can_air-to-atm CO2 xfer [ppm/m^2]

      real, allocatable ::     ggaer  (:) ! surface aerodynamic conductance [m/s]
      real, allocatable :: sea_ggaer  (:) ! surface aerodynamic conductance over water [m/s]
      real, allocatable :: ice_ggaer  (:) ! surface aerodynamic conductance over ice [m/s]

! Radiative flux quantities

      real, allocatable :: albedo_beam   (:) ! water s/w beam albedo [0-1]
      real, allocatable :: albedo_diffuse(:) ! water s/w diffuse albedo [0-1]
      real, allocatable :: sea_albedo    (:) ! water s/w albedo [0-1]
      real, allocatable :: ice_albedo    (:) ! seaice s/w albedo [0-1]

      real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
      real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
      real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
      real, allocatable :: rlong_albedo  (:) ! surface l/w lbedo [0-1]

      real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: sea_rlongup   (:) ! upward can-top l/w flux [W/m^2]
      real, allocatable :: ice_rlongup   (:) ! upward can-top l/w flux [W/m^2]
      
      real, allocatable :: ice_net_rlong (:)
      real, allocatable :: ice_net_rshort(:)

! Precipitation fluxes

      real, allocatable :: pcpg (:) ! new pcp amount this timestep [kg/m^2]
      real, allocatable :: qpcpg(:) ! new pcp energy this timestep [J/m^2]
      real, allocatable :: dpcpg(:) ! new pcp depth this timestep [m]

! Canopy and surface quantities:

      real, allocatable :: can_depth     (:) ! "canopy" depth for heat & vap capacity [m]
      real, allocatable :: cantemp       (:) ! "canopy" air temperature [K]
      real, allocatable :: sea_cantemp   (:) ! "canopy" air temperature [K]
      real, allocatable :: ice_cantemp   (:) ! "canopy" air temperature [K]

      real, allocatable :: canshv        (:) ! "canopy" vapor spec hum [kg_vap/kg_air]
      real, allocatable :: sea_canshv    (:) ! "canopy" vapor spec hum [kg_vap/kg_air]
      real, allocatable :: ice_canshv    (:) ! "canopy" vapor spec hum [kg_vap/kg_air]

      integer, allocatable :: nlev_seaice(:) ! number of seaice layers in current cell

      real, allocatable :: seatp         (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf         (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc         (:) ! current sea temperature [K]
      real, allocatable :: seaicep       (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef       (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec       (:) ! seaice fraction [0-1]

      real, allocatable :: surface_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]
      real, allocatable :: sea_sfc_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]
      real, allocatable :: ice_sfc_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]

      real, allocatable :: rough         (:) ! water surface roughness height [m]
      real, allocatable :: sea_rough     (:) ! water surface roughness height [m]
      real, allocatable :: ice_rough     (:) ! water surface roughness height [m]

      real, allocatable :: seaice_energy(:,:)
      real, allocatable :: seaice_tempk (:,:)

   End Type sea_vars

   type (sea_vars), target :: sea

Contains

!=========================================================================

   subroutine alloc_sea_grid(mws)

     use misc_coms, only: rinit
     use sea_coms,  only: iseagrid

     implicit none

     integer, intent(in) :: mws

!    Allocate sea grid tables

     allocate (itab_ws(mws))

!    Allocate and initialize sea grid information from mksfc

     allocate (sea%area      (mws)) ; sea%area       = rinit
     allocate (sea%glatw     (mws)) ; sea%glatw      = rinit
     allocate (sea%glonw     (mws)) ; sea%glonw      = rinit
     allocate (sea%xew       (mws)) ; sea%xew        = rinit
     allocate (sea%yew       (mws)) ; sea%yew        = rinit
     allocate (sea%zew       (mws)) ; sea%zew        = rinit
     allocate (sea%topw      (mws)) ; sea%topw       = rinit

     allocate (sea%leaf_class(mws)) ; sea%leaf_class = 0

   end subroutine alloc_sea_grid

!=========================================================================

   subroutine alloc_sea(mws)

     use misc_coms, only: rinit
     use sea_coms,  only: nzi

     implicit none

     integer, intent(in) :: mws

!    Allocate and initialize sea arrays

     allocate (sea%rhos          (mws)) ; sea%rhos           = rinit
     allocate (sea%vels          (mws)) ; sea%vels           = rinit
     allocate (sea%prss          (mws)) ; sea%prss           = rinit
     allocate (sea%airtemp       (mws)) ; sea%airtemp        = rinit
     allocate (sea%airshv        (mws)) ; sea%airshv         = rinit

     allocate (sea%ustar         (mws)) ; sea%ustar          = rinit
     allocate (sea%sea_ustar     (mws)) ; sea%sea_ustar      = rinit
     allocate (sea%ice_ustar     (mws)) ; sea%ice_ustar      = rinit

     allocate (sea%vkmsfc        (mws)) ; sea%vkmsfc         = rinit
     allocate (sea%sea_vkmsfc    (mws)) ; sea%sea_vkmsfc     = rinit
     allocate (sea%ice_vkmsfc    (mws)) ; sea%ice_vkmsfc     = rinit

     allocate (sea%sfluxt        (mws)) ; sea%sfluxt         = rinit
     allocate (sea%sea_sfluxt    (mws)) ; sea%sea_sfluxt     = rinit
     allocate (sea%ice_sfluxt    (mws)) ; sea%ice_sfluxt     = rinit

     allocate (sea%sfluxr        (mws)) ; sea%sfluxr         = rinit
     allocate (sea%sea_sfluxr    (mws)) ; sea%sea_sfluxr     = rinit
     allocate (sea%ice_sfluxr    (mws)) ; sea%ice_sfluxr     = rinit

     allocate (sea%sfluxc        (mws)) ; sea%sfluxc         = rinit

     allocate (sea%sxfer_t       (mws)) ; sea%sxfer_t        = 0.0
     allocate (sea%sea_sxfer_t   (mws)) ; sea%sea_sxfer_t    = 0.0
     allocate (sea%ice_sxfer_t   (mws)) ; sea%ice_sxfer_t    = 0.0

     allocate (sea%sxfer_r       (mws)) ; sea%sxfer_r        = 0.0
     allocate (sea%sea_sxfer_r   (mws)) ; sea%sea_sxfer_r    = 0.0
     allocate (sea%ice_sxfer_r   (mws)) ; sea%ice_sxfer_r    = 0.0

     allocate (sea%sxfer_c       (mws)) ; sea%sxfer_c        = 0.0

     allocate (sea%ggaer         (mws)) ; sea%ggaer          = rinit
     allocate (sea%sea_ggaer     (mws)) ; sea%sea_ggaer      = rinit
     allocate (sea%ice_ggaer     (mws)) ; sea%ice_ggaer      = rinit

     allocate (sea%albedo_beam   (mws)) ; sea%albedo_beam    = 0.0
     allocate (sea%albedo_diffuse(mws)) ; sea%albedo_diffuse = 0.0
     allocate (sea%sea_albedo    (mws)) ; sea%sea_albedo     = 0.0
     allocate (sea%ice_albedo    (mws)) ; sea%ice_albedo     = 0.0

     allocate (sea%rshort        (mws)) ; sea%rshort         = 0.0
     allocate (sea%rshort_diffuse(mws)) ; sea%rshort_diffuse = 0.0
     allocate (sea%rlong         (mws)) ; sea%rlong          = 0.0
     allocate (sea%rlong_albedo  (mws)) ; sea%rlong_albedo   = 0.0

     allocate (sea%rlongup       (mws)) ; sea%rlongup        = 0.0
     allocate (sea%sea_rlongup   (mws)) ; sea%sea_rlongup    = 0.0
     allocate (sea%ice_rlongup   (mws)) ; sea%ice_rlongup    = 0.0

     allocate (sea%ice_net_rlong (mws)) ; sea%ice_net_rlong  = 0.0
     allocate (sea%ice_net_rshort(mws)) ; sea%ice_net_rshort = 0.0

     allocate (sea%pcpg          (mws)) ; sea%pcpg           = 0.0
     allocate (sea%qpcpg         (mws)) ; sea%qpcpg          = 0.0
     allocate (sea%dpcpg         (mws)) ; sea%dpcpg          = 0.0

     allocate (sea%can_depth     (mws)) ; sea%can_depth      = rinit

     allocate (sea%cantemp       (mws)) ; sea%cantemp        = rinit
     allocate (sea%sea_cantemp   (mws)) ; sea%sea_cantemp    = rinit
     allocate (sea%ice_cantemp   (mws)) ; sea%ice_cantemp    = rinit

     allocate (sea%canshv        (mws)) ; sea%canshv         = rinit
     allocate (sea%sea_canshv    (mws)) ; sea%sea_canshv     = rinit
     allocate (sea%ice_canshv    (mws)) ; sea%ice_canshv     = rinit

     allocate (sea%nlev_seaice   (mws)) ; sea%nlev_seaice    = 0

     allocate (sea%seatp         (mws)) ; sea%seatp          = rinit
     allocate (sea%seatf         (mws)) ; sea%seatf          = rinit
     allocate (sea%seatc         (mws)) ; sea%seatc          = rinit
     allocate (sea%seaicep       (mws)) ; sea%seaicep        = rinit
     allocate (sea%seaicef       (mws)) ; sea%seaicef        = rinit
     allocate (sea%seaicec       (mws)) ; sea%seaicec        = rinit

     allocate (sea%surface_ssh   (mws)) ; sea%surface_ssh    = rinit
     allocate (sea%sea_sfc_ssh   (mws)) ; sea%sea_sfc_ssh    = rinit
     allocate (sea%ice_sfc_ssh   (mws)) ; sea%ice_sfc_ssh    = rinit

     allocate (sea%rough         (mws)) ; sea%rough          = rinit
     allocate (sea%sea_rough     (mws)) ; sea%sea_rough      = rinit
     allocate (sea%ice_rough     (mws)) ; sea%ice_rough      = rinit

     allocate (sea%seaice_energy(nzi,mws)) ; sea%seaice_energy = rinit
     allocate (sea%seaice_tempk (nzi,mws)) ; sea%seaice_tempk  = rinit

   end subroutine alloc_sea

!=========================================================================

   subroutine filltab_sea()

     use var_tables, only: increment_vtable
     implicit none

     if (allocated(sea%rhos)) call increment_vtable('SEA%RHOS', 'SW', rvar1=sea%rhos)

     if (allocated(sea%vels)) call increment_vtable('SEA%VELS', 'SW', rvar1=sea%vels)

     if (allocated(sea%prss)) call increment_vtable('SEA%PRSS', 'SW', rvar1=sea%prss)

     if (allocated(sea%airtemp)) call increment_vtable('SEA%AIRTEMP', 'SW', rvar1=sea%airtemp)

     if (allocated(sea%airshv)) call increment_vtable('SEA%AIRSHV', 'SW', rvar1=sea%airshv)

     if (allocated(sea%ustar)) call increment_vtable('SEA%USTAR', 'SW', rvar1=sea%ustar)

     if (allocated(sea%sea_ustar)) call increment_vtable('SEA%SEA_USTAR', 'SW', rvar1=sea%sea_ustar)

     if (allocated(sea%ice_ustar)) call increment_vtable('SEA%ICE_USTAR', 'SW', rvar1=sea%ice_ustar)

     if (allocated(sea%vkmsfc)) call increment_vtable('SEA%VKMSFC', 'SW', rvar1=sea%vkmsfc)

     if (allocated(sea%sea_vkmsfc)) call increment_vtable('SEA%SEA_VKMSFC', 'SW', rvar1=sea%sea_vkmsfc)

     if (allocated(sea%ice_vkmsfc)) call increment_vtable('SEA%ICE_VKMSFC', 'SW', rvar1=sea%ice_vkmsfc)

     if (allocated(sea%sfluxt)) call increment_vtable('SEA%SFLUXT', 'SW', rvar1=sea%sfluxt)

     if (allocated(sea%sea_sfluxt)) call increment_vtable('SEA%SEA_SFLUXT', 'SW', rvar1=sea%sea_sfluxt)

     if (allocated(sea%ice_sfluxt)) call increment_vtable('SEA%ICE_SFLUXT', 'SW', rvar1=sea%ice_sfluxt)

     if (allocated(sea%sfluxr)) call increment_vtable('SEA%SFLUXR', 'SW', rvar1=sea%sfluxr)

     if (allocated(sea%sea_sfluxr)) call increment_vtable('SEA%SEA_SFLUXR', 'SW', rvar1=sea%sea_sfluxr)

     if (allocated(sea%ice_sfluxr)) call increment_vtable('SEA%ICE_SFLUXR', 'SW', rvar1=sea%ice_sfluxr)

     if (allocated(sea%sfluxc)) call increment_vtable('SEA%SFLUXC', 'SW', rvar1=sea%sfluxc)

     if (allocated(sea%sxfer_t)) call increment_vtable('SEA%SXFER_T', 'SW', rvar1=sea%sxfer_t)

     if (allocated(sea%sea_sxfer_t)) call increment_vtable('SEA%SEA_SXFER_T', 'SW', rvar1=sea%sea_sxfer_t)

     if (allocated(sea%ice_sxfer_t)) call increment_vtable('SEA%ICE_SXFER_T', 'SW', rvar1=sea%ice_sxfer_t)

     if (allocated(sea%sxfer_r)) call increment_vtable('SEA%SXFER_R', 'SW', rvar1=sea%sxfer_r)

     if (allocated(sea%sea_sxfer_r)) call increment_vtable('SEA%SEA_SXFER_R', 'SW', rvar1=sea%sea_sxfer_r)

     if (allocated(sea%ice_sxfer_r)) call increment_vtable('SEA%ICE_SXFER_R', 'SW', rvar1=sea%ice_sxfer_r)

     if (allocated(sea%sxfer_c)) call increment_vtable('SEA%SXFER_C', 'SW', rvar1=sea%sxfer_c)

     if (allocated(sea%albedo_beam)) call increment_vtable('SEA%ALBEDO_BEAM', 'SW', rvar1=sea%albedo_beam)

     if (allocated(sea%albedo_diffuse)) call increment_vtable('SEA%ALBEDO_DIFFUSE', 'SW', rvar1=sea%albedo_diffuse)

     if (allocated(sea%sea_albedo)) call increment_vtable('SEA%SEA_ALBEDO', 'SW', rvar1=sea%sea_albedo)

      if (allocated(sea%ice_albedo)) call increment_vtable('SEA%ICE_ALBEDO', 'SW', rvar1=sea%ice_albedo)

     if (allocated(sea%rshort)) call increment_vtable('SEA%RSHORT', 'SW', rvar1=sea%rshort)

     if (allocated(sea%rshort_diffuse)) call increment_vtable('SEA%RSHORT_DIFFUSE', 'SW', rvar1=sea%rshort_diffuse)

     if (allocated(sea%rlong)) call increment_vtable('SEA%RLONG', 'SW', rvar1=sea%rlong)

     if (allocated(sea%rlong_albedo)) call increment_vtable('SEA%RLONG_ALBEDO', 'SW', rvar1=sea%rlong_albedo)

     if (allocated(sea%rlongup)) call increment_vtable('SEA%RLONGUP', 'SW', rvar1=sea%rlongup)

     if (allocated(sea%sea_rlongup)) call increment_vtable('SEA%SEA_RLONGUP', 'SW', rvar1=sea%sea_rlongup)

     if (allocated(sea%ice_rlongup)) call increment_vtable('SEA%ICE_RLONGUP', 'SW', rvar1=sea%ice_rlongup)

      if (allocated(sea%ice_net_rlong)) call increment_vtable('SEA%ICE_NET_RLONG', 'SW', rvar1=sea%ice_net_rlong)

      if (allocated(sea%ice_net_rshort)) call increment_vtable('SEA%ICE_NET_RSHORT', 'SW', rvar1=sea%ice_net_rshort)

     if (allocated(sea%pcpg)) call increment_vtable('SEA%PCPG', 'SW', rvar1=sea%pcpg)

     if (allocated(sea%qpcpg)) call increment_vtable('SEA%QPCPG', 'SW', rvar1=sea%qpcpg)

     if (allocated(sea%dpcpg)) call increment_vtable('SEA%DPCPG', 'SW', rvar1=sea%dpcpg)

     if (allocated(sea%can_depth)) call increment_vtable('SEA%CAN_DEPTH', 'SW', rvar1=sea%can_depth)

     if (allocated(sea%seatc)) call increment_vtable('SEA%SEATC', 'SW', rvar1=sea%seatc)

     if (allocated(sea%cantemp)) call increment_vtable('SEA%CANTEMP', 'SW', rvar1=sea%cantemp)

     if (allocated(sea%sea_cantemp)) call increment_vtable('SEA%SEA_CANTEMP', 'SW', rvar1=sea%sea_cantemp)

     if (allocated(sea%ice_cantemp)) call increment_vtable('SEA%ICE_CANTEMP', 'SW', rvar1=sea%ice_cantemp)

     if (allocated(sea%canshv)) call increment_vtable('SEA%CANSHV', 'SW', rvar1=sea%canshv)

     if (allocated(sea%sea_canshv)) call increment_vtable('SEA%SEA_CANSHV', 'SW', rvar1=sea%sea_canshv)

     if (allocated(sea%ice_canshv)) call increment_vtable('SEA%ICE_CANSHV', 'SW', rvar1=sea%ice_canshv)

     if (allocated(sea%nlev_seaice)) call increment_vtable('SEA%NLEV_SEAICE', 'SW', ivar1=sea%nlev_seaice)

     if (allocated(sea%seaicec)) call increment_vtable('SEA%SEAICEC', 'SW', rvar1=sea%seaicec)

     if (allocated(sea%surface_ssh)) call increment_vtable('SEA%SURFACE_SSH', 'SW', rvar1=sea%surface_ssh)

     if (allocated(sea%sea_sfc_ssh)) call increment_vtable('SEA%SEA_SFC_SSH', 'SW', rvar1=sea%sea_sfc_ssh)

     if (allocated(sea%ice_sfc_ssh)) call increment_vtable('SEA%ICE_SFC_SSH', 'SW', rvar1=sea%ice_sfc_ssh)

     if (allocated(sea%rough)) call increment_vtable('SEA%ROUGH', 'SW', rvar1=sea%rough)

     if (allocated(sea%sea_rough)) call increment_vtable('SEA%SEA_ROUGH', 'SW', rvar1=sea%sea_rough)

     if (allocated(sea%ice_rough)) call increment_vtable('SEA%ICE_ROUGH', 'SW', rvar1=sea%ice_rough)

     if (allocated(sea%seaice_energy)) call increment_vtable('SEA%SEAICE_ENERGY', 'SW', rvar2=sea%seaice_energy)
    
     if (allocated(sea%seaice_tempk)) call increment_vtable('SEA%SEAICE_TEMPK', 'SW', rvar2=sea%seaice_tempk)
    
   end subroutine filltab_sea

End Module mem_sea
