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

      real, allocatable :: sxfer_tsav (:) ! saved previous value of sxfer_t
      real, allocatable :: sxfer_rsav (:) ! saved previous value of sxter_r
      real, allocatable :: sxfer_csav (:) ! saved previous value of sxter_c

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

     allocate (sea%sxfer_tsav    (mws)) ; sea%sxfer_tsav     = 0.0
     allocate (sea%sxfer_rsav    (mws)) ; sea%sxfer_rsav     = 0.0
     allocate (sea%sxfer_csav    (mws)) ; sea%sxfer_csav     = 0.0

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

     use var_tables, only: vtab_r, num_var, increment_vtable
     implicit none

     if (allocated(sea%rhos)) then
        call increment_vtable('SEA%RHOS', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rhos
     endif

     if (allocated(sea%vels)) then
        call increment_vtable('SEA%VELS', 'SW')
        vtab_r(num_var)%rvar1_p => sea%vels
     endif

     if (allocated(sea%prss)) then
        call increment_vtable('SEA%PRSS', 'SW')
        vtab_r(num_var)%rvar1_p => sea%prss
     endif

     if (allocated(sea%airtemp)) then
        call increment_vtable('SEA%AIRTEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%airtemp
     endif

     if (allocated(sea%airshv)) then
        call increment_vtable('SEA%AIRSHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%airshv
     endif

     if (allocated(sea%ustar)) then
        call increment_vtable('SEA%USTAR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ustar
     endif

     if (allocated(sea%sea_ustar)) then
        call increment_vtable('SEA%SEA_USTAR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_ustar
     endif

     if (allocated(sea%ice_ustar)) then
        call increment_vtable('SEA%ICE_USTAR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_ustar
     endif

     if (allocated(sea%vkmsfc)) then
        call increment_vtable('SEA%VKMSFC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%vkmsfc
     endif

     if (allocated(sea%sea_vkmsfc)) then
        call increment_vtable('SEA%SEA_VKMSFC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_vkmsfc
     endif

     if (allocated(sea%ice_vkmsfc)) then
        call increment_vtable('SEA%ICE_VKMSFC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_vkmsfc
     endif

     if (allocated(sea%sfluxt)) then
        call increment_vtable('SEA%SFLUXT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sfluxt
     endif

     if (allocated(sea%sea_sfluxt)) then
        call increment_vtable('SEA%SEA_SFLUXT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_sfluxt
     endif

     if (allocated(sea%ice_sfluxt)) then
        call increment_vtable('SEA%ICE_SFLUXT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_sfluxt
     endif

     if (allocated(sea%sfluxr)) then
        call increment_vtable('SEA%SFLUXR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sfluxr
     endif

     if (allocated(sea%sea_sfluxr)) then
        call increment_vtable('SEA%SEA_SFLUXR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_sfluxr
     endif

     if (allocated(sea%ice_sfluxr)) then
        call increment_vtable('SEA%ICE_SFLUXR', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_sfluxr
     endif

     if (allocated(sea%sfluxc)) then
        call increment_vtable('SEA%SFLUXC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sfluxc
     endif

     if (allocated(sea%sxfer_t)) then
        call increment_vtable('SEA%SXFER_T', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sxfer_t
     endif

     if (allocated(sea%sea_sxfer_t)) then
        call increment_vtable('SEA%SEA_SXFER_T', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_sxfer_t
     endif

     if (allocated(sea%ice_sxfer_t)) then
        call increment_vtable('SEA%ICE_SXFER_T', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_sxfer_t
     endif

     if (allocated(sea%sxfer_r)) then
        call increment_vtable('SEA%SXFER_R', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sxfer_r
     endif

     if (allocated(sea%sea_sxfer_r)) then
        call increment_vtable('SEA%SEA_SXFER_R', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_sxfer_r
     endif

     if (allocated(sea%ice_sxfer_r)) then
        call increment_vtable('SEA%ICE_SXFER_R', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_sxfer_r
     endif

     if (allocated(sea%sxfer_c)) then
        call increment_vtable('SEA%SXFER_C', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sxfer_c
     endif

     if (allocated(sea%albedo_beam)) then
        call increment_vtable('SEA%ALBEDO_BEAM', 'SW')
        vtab_r(num_var)%rvar1_p => sea%albedo_beam
     endif

     if (allocated(sea%albedo_diffuse)) then
        call increment_vtable('SEA%ALBEDO_DIFFUSE', 'SW')
        vtab_r(num_var)%rvar1_p => sea%albedo_diffuse
     endif

     if (allocated(sea%sea_albedo)) then
        call increment_vtable('SEA%SEA_ALBEDO', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_albedo
     endif

      if (allocated(sea%ice_albedo)) then
        call increment_vtable('SEA%ICE_ALBEDO', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_albedo
     endif

     if (allocated(sea%rshort)) then
        call increment_vtable('SEA%RSHORT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rshort
     endif

     if (allocated(sea%rshort_diffuse)) then
        call increment_vtable('SEA%RSHORT_DIFFUSE', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rshort_diffuse
     endif

     if (allocated(sea%rlong)) then
        call increment_vtable('SEA%RLONG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rlong
     endif

     if (allocated(sea%rlong_albedo)) then
        call increment_vtable('SEA%RLONG_ALBEDO', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rlong_albedo
     endif

     if (allocated(sea%rlongup)) then
        call increment_vtable('SEA%RLONGUP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rlongup
     endif

     if (allocated(sea%sea_rlongup)) then
        call increment_vtable('SEA%SEA_RLONGUP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_rlongup
     endif

     if (allocated(sea%ice_rlongup)) then
        call increment_vtable('SEA%ICE_RLONGUP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_rlongup
     endif

      if (allocated(sea%ice_net_rlong)) then
        call increment_vtable('SEA%ICE_NET_RLONG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_net_rlong
     endif

      if (allocated(sea%ice_net_rshort)) then
        call increment_vtable('SEA%ICE_NET_RSHORT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_net_rshort
     endif

     if (allocated(sea%pcpg)) then
        call increment_vtable('SEA%PCPG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%pcpg
     endif

     if (allocated(sea%qpcpg)) then
        call increment_vtable('SEA%QPCPG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%qpcpg
     endif

     if (allocated(sea%dpcpg)) then
        call increment_vtable('SEA%DPCPG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%dpcpg
     endif

     if (allocated(sea%can_depth)) then
        call increment_vtable('SEA%CAN_DEPTH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%can_depth
     endif

     if (allocated(sea%seatc)) then

     if (allocated(sea%cantemp)) then
        call increment_vtable('SEA%CANTEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%cantemp
     endif

     if (allocated(sea%sea_cantemp)) then
        call increment_vtable('SEA%SEA_CANTEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_cantemp
     endif

     if (allocated(sea%ice_cantemp)) then
        call increment_vtable('SEA%ICE_CANTEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_cantemp
     endif

     if (allocated(sea%canshv)) then
        call increment_vtable('SEA%CANSHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%canshv
     endif

     if (allocated(sea%sea_canshv)) then
        call increment_vtable('SEA%SEA_CANSHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_canshv
     endif

     if (allocated(sea%ice_canshv)) then
        call increment_vtable('SEA%ICE_CANSHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_canshv
     endif

     if (allocated(sea%nlev_seaice)) then
        call increment_vtable('SEA%NLEV_SEAICE', 'SW')
        vtab_r(num_var)%ivar1_p => sea%nlev_seaice
     endif

        call increment_vtable('SEA%SEATC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seatc
     endif

     if (allocated(sea%seaicec)) then
        call increment_vtable('SEA%SEAICEC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seaicec
     endif

     if (allocated(sea%surface_ssh)) then
        call increment_vtable('SEA%SURFACE_SSH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%surface_ssh
     endif

     if (allocated(sea%sea_sfc_ssh)) then
        call increment_vtable('SEA%SEA_SFC_SSH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_sfc_ssh
     endif

     if (allocated(sea%ice_sfc_ssh)) then
        call increment_vtable('SEA%ICE_SFC_SSH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_sfc_ssh
     endif

     if (allocated(sea%rough)) then
        call increment_vtable('SEA%ROUGH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rough
     endif

     if (allocated(sea%sea_rough)) then
        call increment_vtable('SEA%SEA_ROUGH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%sea_rough
     endif

     if (allocated(sea%ice_rough)) then
        call increment_vtable('SEA%ICE_ROUGH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_rough
     endif

     if (allocated(sea%seaice_energy)) then
        call increment_vtable('SEA%SEAICE_ENERGY', 'SW')
        vtab_r(num_var)%rvar2_p => sea%seaice_energy
     endif
    
     if (allocated(sea%seaice_tempk)) then
        call increment_vtable('SEA%SEAICE_TEMPK', 'SW')
        vtab_r(num_var)%rvar2_p => sea%seaice_tempk
     endif
    
   end subroutine filltab_sea

End Module mem_sea
