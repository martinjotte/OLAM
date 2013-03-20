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
Module mem_sea

   use max_dims,  only: maxremote
   use mem_mksfc, only: itab_mls_vars, itab_uls_vars, itab_wls_vars
   implicit none

   private :: maxremote, itab_mls_vars, itab_uls_vars, itab_wls_vars

! SEA SURFACE GRID TABLES

   type (itab_mls_vars), allocatable, target :: itab_ms(:)
   type (itab_uls_vars), allocatable, target :: itab_us(:)
   type (itab_wls_vars), allocatable, target :: itab_ws(:)

   type (itab_uls_vars), allocatable, target :: itab_us_pd(:)
   type (itab_wls_vars), allocatable, target :: itab_ws_pd(:)

! DERIVED TYPES TO HOLD GLOBAL GRID INDICES FOR A PARALLEL RUN

   Type itabg_ms_vars
      integer :: ims_myrank = -1
      integer :: irank      = -1
   End Type itabg_ms_vars

   Type itabg_us_vars
      integer :: ius_myrank = -1
      integer :: irank      = -1
   End Type itabg_us_vars

   Type itabg_ws_vars
      integer :: iws_myrank = -1
      integer :: irank      = -1
   End Type itabg_ws_vars

   type (itabg_ms_vars), allocatable, target :: itabg_ms(:)
   type (itabg_us_vars), allocatable, target :: itabg_us(:)
   type (itabg_ws_vars), allocatable, target :: itabg_ws(:)

! DERIVED TYPE TO HOLD MPI TABLES FOR A PARALLEL RUN

   Type jtab_ws_mpi_vars
      integer, allocatable :: iws(:)
      integer, allocatable :: jend(:)
   End Type jtab_ws_mpi_vars

   type (jtab_ws_mpi_vars) :: jtab_ws_mpi(maxremote)

! SEA SURFACE MODEL VARIABLES

   Type sea_vars
      integer, allocatable :: leaf_class (:) ! sea cell leaf class
      integer, allocatable :: nlev_seaice(:) ! number of seaice layers in current cell

      real, allocatable :: area          (:) ! sea cell surface area [m^2]
      real, allocatable :: xew           (:) ! earth x coord of sea cell W points
      real, allocatable :: yew           (:) ! earth y coord of sea cell W points
      real, allocatable :: zew           (:) ! earth z coord of sea cell W points
      real, allocatable :: glatw         (:) ! latitude of sea cell W points
      real, allocatable :: glonw         (:) ! longitude of sea cell W points

      real, allocatable :: xem           (:)  ! earth x coord of sea cell M points
      real, allocatable :: yem           (:)  ! earth y coord of sea cell M points
      real, allocatable :: zem           (:)  ! earth z coord of sea cell M points
      real, allocatable :: zm            (:)  ! topography height at sea cell M points
      real, allocatable :: glatm         (:)  ! latitude of sea cell M points
      real, allocatable :: glonm         (:)  ! longitude of sea cell M points

      real, allocatable :: rhos          (:) ! air density [kg/m^3]

      real, allocatable :: ustar         (:) ! friction velocity [m/s]
      real, allocatable :: sea_ustar     (:) ! friction velocity [m/s]
      real, allocatable :: ice_ustar     (:) ! friction velocity [m/s]

      real, allocatable :: sxfer_t       (:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: sea_sxfer_t   (:) ! can_air-to-atm heat xfer [kg_air K/m^2]
      real, allocatable :: ice_sxfer_t   (:) ! can_air-to-atm heat xfer [kg_air K/m^2]

      real, allocatable :: sxfer_r       (:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: sea_sxfer_r   (:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
      real, allocatable :: ice_sxfer_r   (:) ! can_air-to-atm vapor xfer [kg_vap/m^2]

      real, allocatable :: sxfer_c       (:) ! can_air-to-atm CO2 xfer [ppm/m^2]

      real, allocatable :: sxfer_tsav    (:) ! saved previous value of sxfer_t
      real, allocatable :: sxfer_rsav    (:) ! saved previous value of sxter_r
      real, allocatable :: sxfer_csav    (:) ! saved previous value of sxter_c

      real, allocatable :: can_depth     (:) ! "canopy" depth for heat & vap capacity [m]
      real, allocatable :: seatp         (:) ! past sea temperature (obs time) [K]
      real, allocatable :: seatf         (:) ! future sea temperature (obs time) [K]
      real, allocatable :: seatc         (:) ! current sea temperature [K]
      real, allocatable :: seaicep       (:) ! past seaice fraction (obs time) [0-1]
      real, allocatable :: seaicef       (:) ! future seaice fraction (obs time) [0-1]
      real, allocatable :: seaicec       (:) ! seaice fraction [0-1]

      real, allocatable :: can_temp      (:) ! "canopy" air temperature [K]
      real, allocatable :: seacan_temp   (:) ! "canopy" air temperature [K]
      real, allocatable :: icecan_temp   (:) ! "canopy" air temperature [K]

      real, allocatable :: can_shv       (:) ! "canopy" vapor spec hum [kg_vap/kg_air]
      real, allocatable :: seacan_shv    (:) ! "canopy" vapor spec hum [kg_vap/kg_air]
      real, allocatable :: icecan_shv    (:) ! "canopy" vapor spec hum [kg_vap/kg_air]

      real, allocatable :: surface_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]
      real, allocatable :: sea_sfc_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]
      real, allocatable :: ice_sfc_ssh   (:) ! sea surface sat spec hum [kg_vap/kg_air]

      real, allocatable :: rough         (:) ! water surface roughness height [m]
      real, allocatable :: sea_rough     (:) ! water surface roughness height [m]
      real, allocatable :: ice_rough     (:) ! water surface roughness height [m]

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

      real, allocatable :: seaice_energy(:,:)
      real, allocatable :: seaice_tempk (:,:)

      real, allocatable :: pcpg          (:) ! new pcp amount this timestep [kg/m^2]
      real, allocatable :: qpcpg         (:) ! new pcp energy this timestep [J/m^2]
      real, allocatable :: dpcpg         (:) ! new pcp depth this timestep [m]

      real, allocatable ::     ggaer     (:) ! surface aerodynamic conductance [m/s]
      real, allocatable :: sea_ggaer     (:) ! surface aerodynamic conductance over water [m/s]
      real, allocatable :: ice_ggaer     (:) ! surface aerodynamic conductance over ice [m/s]

   End Type sea_vars

   type (sea_vars), target :: sea

   Type sea_pd_vars
      real, allocatable :: xem           (:)  ! earth x coord of sea cell M points
      real, allocatable :: yem           (:)  ! earth y coord of sea cell M points
      real, allocatable :: zem           (:)  ! earth z coord of sea cell M points
   End Type sea_pd_vars

   type (sea_vars), target :: sea_pd

Contains

!=========================================================================

   subroutine alloc_sea_grid(mms,mus,mws)
     use misc_coms, only: rinit
     implicit none

     integer, intent(in) :: mms, mus, mws

!    Allocate sea grid tables

     allocate (itab_ms(mms))
     allocate (itab_us(mus))
     allocate (itab_ws(mws))

!    Allocate and initialize sea grid information from mksfc

     allocate (sea%leaf_class(mws)) ; sea%leaf_class = 0
     allocate (sea%area      (mws)) ; sea%area       = rinit
     allocate (sea%xew       (mws)) ; sea%xew        = rinit
     allocate (sea%yew       (mws)) ; sea%yew        = rinit
     allocate (sea%zew       (mws)) ; sea%zew        = rinit
     allocate (sea%glatw     (mws)) ; sea%glatw      = rinit
     allocate (sea%glonw     (mws)) ; sea%glonw      = rinit

     allocate (sea%xem       (mms)) ; sea%xem        = rinit
     allocate (sea%yem       (mms)) ; sea%yem        = rinit
     allocate (sea%zem       (mms)) ; sea%zem        = rinit
     allocate (sea%zm        (mms)) ; sea%zm         = rinit
     allocate (sea%glatm     (mms)) ; sea%glatm      = rinit
     allocate (sea%glonm     (mms)) ; sea%glonm      = rinit

   end subroutine alloc_sea_grid

!=========================================================================

   subroutine alloc_sea_grid_pd(mms,mus,mws)
     use misc_coms, only: rinit
     implicit none

     integer, intent(in) :: mms, mus, mws

!    Allocate sea grid tables

     allocate (itab_us_pd(mus))
     allocate (itab_ws_pd(mws))

     allocate (sea_pd%xem       (mms)) ; sea_pd%xem        = rinit
     allocate (sea_pd%yem       (mms)) ; sea_pd%yem        = rinit
     allocate (sea_pd%zem       (mms)) ; sea_pd%zem        = rinit

   end subroutine alloc_sea_grid_pd

!=========================================================================

   subroutine alloc_sea(mws)

     use misc_coms, only: rinit
     use sea_coms,  only: nzi
     implicit none

     integer, intent(in) :: mws

!    Allocate and initialize sea arrays

     allocate (sea%nlev_seaice   (mws)) ; sea%nlev_seaice    = 0

     allocate (sea%rhos          (mws)) ; sea%rhos           = rinit

     allocate (sea%ustar         (mws)) ; sea%ustar          = rinit
     allocate (sea%sea_ustar     (mws)) ; sea%sea_ustar      = rinit
     allocate (sea%ice_ustar     (mws)) ; sea%ice_ustar      = rinit

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

     allocate (sea%can_depth     (mws)) ; sea%can_depth      = rinit
     allocate (sea%seatp         (mws)) ; sea%seatp          = rinit
     allocate (sea%seatf         (mws)) ; sea%seatf          = rinit
     allocate (sea%seatc         (mws)) ; sea%seatc          = rinit
     allocate (sea%seaicep       (mws)) ; sea%seaicep        = rinit
     allocate (sea%seaicef       (mws)) ; sea%seaicef        = rinit
     allocate (sea%seaicec       (mws)) ; sea%seaicec        = rinit

     allocate (sea%can_temp      (mws)) ; sea%can_temp       = rinit
     allocate (sea%seacan_temp   (mws)) ; sea%seacan_temp    = rinit
     allocate (sea%icecan_temp   (mws)) ; sea%icecan_temp    = rinit

     allocate (sea%can_shv       (mws)) ; sea%can_shv        = rinit
     allocate (sea%seacan_shv    (mws)) ; sea%seacan_shv     = rinit
     allocate (sea%icecan_shv    (mws)) ; sea%icecan_shv     = rinit

     allocate (sea%surface_ssh   (mws)) ; sea%surface_ssh    = rinit
     allocate (sea%sea_sfc_ssh   (mws)) ; sea%sea_sfc_ssh    = rinit
     allocate (sea%ice_sfc_ssh   (mws)) ; sea%ice_sfc_ssh    = rinit

     allocate (sea%rough         (mws)) ; sea%rough          = rinit
     allocate (sea%sea_rough     (mws)) ; sea%sea_rough      = rinit
     allocate (sea%ice_rough     (mws)) ; sea%ice_rough      = rinit

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

     allocate (sea%seaice_energy(nzi,mws)) ; sea%seaice_energy = rinit
     allocate (sea%seaice_tempk (nzi,mws)) ; sea%seaice_tempk  = rinit

     allocate (sea%ggaer         (mws)) ; sea%ggaer          = rinit
     allocate (sea%sea_ggaer     (mws)) ; sea%sea_ggaer      = rinit
     allocate (sea%ice_ggaer     (mws)) ; sea%ice_ggaer      = rinit

   end subroutine alloc_sea

!=========================================================================

   subroutine filltab_sea()

     use var_tables, only: vtab_r, num_var, increment_vtable
     use misc_coms,  only: iparallel, ipar_out, runtype
     implicit none

     if ( (runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE') .or. &
          (iparallel == 1 .and. ipar_out == 0) ) then

        ! THESE ONLY NEED TO BE WRITTEN TO HISTORY FILE FOR PARALLEL RUNS
        ! WRITING TO MULTIPLE FILES, AND READ FOR PLOTONLY OR PARCOMBINE RUNS

        if (allocated(itab_ms)) then

           call increment_vtable('IMGLOBE_S', 'SM', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ms%imglobe
        endif

        if (allocated(itab_us)) then

           call increment_vtable('IUGLOBE_S', 'SU', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_us%iuglobe

           call increment_vtable('IRANKU_S',  'SU', noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_us%irank

        endif

        if (allocated(itab_ws)) then 

           call increment_vtable('IWGLOBE_S', 'SW',  noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ws%iwglobe

           call increment_vtable('IRANKW_S',  'SW',  noread=.true.)
           vtab_r(num_var)%ivar1_p => itab_ws%irank

        endif

     endif

     if (allocated(sea%nlev_seaice)) then
        call increment_vtable('SEA%NLEV_SEAICE', 'SW')
        vtab_r(num_var)%ivar1_p => sea%nlev_seaice
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

     if (allocated(sea%can_depth)) then
        call increment_vtable('SEA%CAN_DEPTH', 'SW')
        vtab_r(num_var)%rvar1_p => sea%can_depth
     endif

     if (allocated(sea%seatc)) then
        call increment_vtable('SEA%SEATC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seatc
     endif

     if (allocated(sea%seaicec)) then
        call increment_vtable('SEA%SEAICEC', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seaicec
     endif

     if (allocated(sea%can_temp)) then
        call increment_vtable('SEA%CAN_TEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%can_temp
     endif

     if (allocated(sea%seacan_temp)) then
        call increment_vtable('SEA%SEACAN_TEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seacan_temp
     endif

     if (allocated(sea%icecan_temp)) then
        call increment_vtable('SEA%ICECAN_TEMP', 'SW')
        vtab_r(num_var)%rvar1_p => sea%icecan_temp
     endif

     if (allocated(sea%can_shv)) then
        call increment_vtable('SEA%CAN_SHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%can_shv
     endif

     if (allocated(sea%seacan_shv)) then
        call increment_vtable('SEA%SEACAN_SHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%seacan_shv
     endif

     if (allocated(sea%icecan_shv)) then
        call increment_vtable('SEA%ICECAN_SHV', 'SW')
        vtab_r(num_var)%rvar1_p => sea%icecan_shv
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

     if (allocated(sea%rlong_albedo)) then
        call increment_vtable('SEA%RLONG_ALBEDO', 'SW')
        vtab_r(num_var)%rvar1_p => sea%rlong_albedo
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

      if (allocated(sea%ice_net_rlong)) then
        call increment_vtable('SEA%ICE_NET_RLONG', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_net_rlong
     endif

      if (allocated(sea%ice_net_rshort)) then
        call increment_vtable('SEA%ICE_NET_RSHORT', 'SW')
        vtab_r(num_var)%rvar1_p => sea%ice_net_rshort
     endif

     if (allocated(sea%seaice_energy)) then
        call increment_vtable('SEA%SEAICE_ENERGY', 'SW')
        vtab_r(num_var)%rvar2_p => sea%seaice_energy
     endif
    
     if (allocated(sea%seaice_tempk)) then
        call increment_vtable('SEA%SEAICE_TEMPK', 'SW')
        vtab_r(num_var)%rvar2_p => sea%seaice_tempk
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

   end subroutine filltab_sea

!===============================================================================

   subroutine fill_jsea()

   use mem_ijtabs, only: mrls

   use misc_coms,  only: io6, iparallel

   use mem_para,   only: mgroupsize, myrank,  &
!                        send_us, recv_us,  &
                         send_ws, recv_ws,  &
                         send_wsf, recv_wsf,  &
!                        nsends_us, nrecvs_us,  &
                         nsends_ws, nrecvs_ws,  &
                         nsends_wsf, nrecvs_wsf

   use sea_coms,   only: mms, mus, mws, nms, nus, nws

   implicit none
   
   integer :: jsend, iws, jend

! Allocate and zero-fill JTAB_WS_MPI%JEND

   do jsend = 1,maxremote
      allocate (jtab_ws_mpi(jsend)%jend(mrls))
                jtab_ws_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Return if run is not parallel (jtab not needed)

   if (iparallel == 0) return
   
! Compute and store JTAB_WS_MPI%JEND(1)

   do jsend = 1,nsends_ws(1)
      jtab_ws_mpi(jsend)%jend(1) = 0
      do iws = 2,mws
         if (itab_ws(iws)%send(jsend)) then
            jtab_ws_mpi(jsend)%jend(1) = jtab_ws_mpi(jsend)%jend(1) + 1
         endif
      enddo
      jtab_ws_mpi(jsend)%jend(1) = max(1,jtab_ws_mpi(jsend)%jend(1))
   enddo

! Allocate and zero-fill JTAB_WS_MPI%IWS

   do jsend = 1,nsends_ws(1)
      jend = jtab_ws_mpi(jsend)%jend(1)
      allocate (jtab_ws_mpi(jsend)%iws(jend))
                jtab_ws_mpi(jsend)%iws(1:jend) = 0
   enddo

! Initialize JTAB_WS_MPI%JEND counters to zero

   do jsend = 1,nsends_ws(1)
      jtab_ws_mpi(jsend)%jend(1:mrls) = 0
   enddo

! Compute JTAB_WS_MPI%IWS (only fill for MRL = 1)

   do iws = 2,mws
      do jsend = 1,nsends_ws(1)

         if (itab_ws(iws)%send(jsend)) then
            jtab_ws_mpi(jsend)%jend(1) = jtab_ws_mpi(jsend)%jend(1) + 1
            jtab_ws_mpi(jsend)%iws(jtab_ws_mpi(jsend)%jend(1)) = iws
         endif

      enddo
   enddo

   return
   end subroutine fill_jsea

End Module mem_sea
