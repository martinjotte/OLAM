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
Module mem_sfcg
  
  use max_dims,   only: maxnlspoly, maxgrds, maxngrdll, pathlen, maxremote
  use mem_ijtabs, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                        nest_ud_vars, nest_wd_vars

  implicit none

  character(pathlen) :: sfcgfile

! HEX-GRID INFORMATION FOR LAND/SEA CELLS (ALL GLOBAL INDICES)

  integer :: nmsfc, mmsfc
  integer :: nvsfc, mvsfc
  integer :: nwsfc, mwsfc

  Type itab_msfc_vars
     integer :: imglobe = 1
     integer :: ivn(3) = 1       ! array of V neighbors of this M pt
     integer :: iwn(3) = 1       ! array of W neighbors of this M pt
  End type itab_msfc_vars

  Type itab_vsfc_vars
     integer :: ivglobe = 1
     integer :: imn(2)
     integer :: iwn(2)
  End type itab_vsfc_vars

  Type itab_wsfc_vars
     logical :: send(maxremote) = .false.
     integer :: iwglobe = 1  ! global sfcg index of this WSFC pt (in parallel run)
     integer :: irank = -1   ! rank of parallel process at this WSFC pt
     integer :: ivoronoi = 0 ! Force to be Voronoi cell if >= 2

     integer :: npoly = 0

     integer :: imn(7)
     integer :: ivn(7)
     integer :: iwn(7)

     ! Dimension of (8) is estimate of max possible number of atm cells that could
     ! overlap this sfc cell

     integer :: nwatm        = 0  ! number of atm cells coupled to this sfc cell
     integer :: iwatm    (8) = 0  ! local-rank atm index of atm cell coupled to this sfc cell
     integer :: kwatm    (8) = 0  ! vertical index of atm cell coupled to this sfc cell
     real    :: arc      (8) = 0. ! coupling area between this sfc cell and 1 atm cell [m^2]
     real    :: arcoarsfc(8) = 0. ! arc over (divided by) area of sfc cell
     real    :: arcoariw (8) = 0. ! arc over area of atm cell
     real    :: arcoarkw (8) = 0. ! arc over area of atm cell that contacts sfc at kw level
  End type itab_wsfc_vars

  type (itab_msfc_vars), allocatable, target :: itab_msfc(:)
  type (itab_vsfc_vars), allocatable, target :: itab_vsfc(:)
  type (itab_wsfc_vars), allocatable, target :: itab_wsfc(:)

  Type itab_msfc_pd_vars
     integer :: ivn(3)
     integer :: iwn(3)
  End type itab_msfc_pd_vars

  Type itab_vsfc_pd_vars
     integer :: imn(2)
     integer :: iwn(2)
  End type itab_vsfc_pd_vars

  Type itab_wsfc_pd_vars
     integer :: ivoronoi = 0 ! Force Voronoi cell if >= 2

     integer :: npoly = 0

     integer :: imn(7)
     integer :: ivn(7)
     integer :: iwn(7)

     integer :: nwatm    = 0
     integer :: iwatm(8) = 0  ! dimension is estimate of max possible

     integer :: leaf_class = -1
  End type itab_wsfc_pd_vars

  type (itab_msfc_pd_vars), allocatable, target :: itab_msfc_pd(:)
  type (itab_vsfc_pd_vars), allocatable, target :: itab_vsfc_pd(:)
  type (itab_wsfc_pd_vars), allocatable, target :: itab_wsfc_pd(:)

  Type itabg_msfc_vars            ! Global data structure for MSFC pts
     integer :: imsfc_myrank = -1 ! local (parallel subdomain) index of a MSFC pt
     integer :: irank = -1        ! rank of parallel process at a MSFC pt
  End Type itabg_msfc_vars

  Type itabg_vsfc_vars            ! Global data structure for VSFC pts
     integer :: ivsfc_myrank = -1 ! local (parallel subdomain) index of a VSFC pt
     integer :: irank = -1        ! rank of parallel process at a VSFC pt
  End Type itabg_vsfc_vars

  Type itabg_wsfc_vars            ! Global data structure for WSFC pts
     integer :: iwsfc_myrank = -1 ! local (parallel subdomain) WSFC index 
     integer :: iland_myrank = -1 ! local (parallel subdomain) LAND index 
     integer :: ilake_myrank = -1 ! local (parallel subdomain) LAKE index 
     integer :: isea_myrank  = -1 ! local (parallel subdomain) SEA index 
     integer :: irank = -1        ! rank of parallel process at a WSFC pt
  End Type itabg_wsfc_vars

  type (itabg_msfc_vars), allocatable, target :: itabg_msfc(:)
  type (itabg_vsfc_vars), allocatable, target :: itabg_vsfc(:)
  type (itabg_wsfc_vars), allocatable, target :: itabg_wsfc(:)

  ! DERIVED TYPE TO HOLD MPI TABLES FOR A PARALLEL RUN

  Type jtab_wsfc_mpi_vars
     integer, allocatable :: iwsfc(:)
     integer, allocatable :: jend(:)
  End Type jtab_wsfc_mpi_vars

  type (jtab_wsfc_mpi_vars) :: jtab_wsfc_mpi(maxremote) 

  Type surface_grid_vars

     ! Surface grid geometry and location

     real, allocatable :: xem  (:) ! earth x coord of sfc M points
     real, allocatable :: yem  (:) ! earth y coord of sfc M points
     real, allocatable :: zem  (:) ! earth z coord of sfc M points
     real, allocatable :: glatm(:) ! latitude of sfc cell M points
     real, allocatable :: glonm(:) ! longitude of sfc cell M points
     real, allocatable :: topm (:) ! topographic height of sfc W points

     real, allocatable :: xev (:) !
     real, allocatable :: yev (:) !
     real, allocatable :: zev (:) !
     real, allocatable :: dnu (:) ! 
     real, allocatable :: dniu(:) ! 
     real, allocatable :: dnv (:) ! 
     real, allocatable :: dniv(:) ! 

     real, allocatable :: area    (:) ! cell surface area [m^2]
     real, allocatable :: xew     (:) ! earth x coord of sfc W points
     real, allocatable :: yew     (:) ! earth y coord of sfc W points
     real, allocatable :: zew     (:) ! earth z coord of sfc W points
     real, allocatable :: glatw   (:) ! latitude of sfc cell W points
     real, allocatable :: glonw   (:) ! longitude of sfc cell W points
     real, allocatable :: topw    (:) ! topographic height of sfc W points
     real, allocatable :: wnx     (:) ! norm unit vector x comp of sfc cells
     real, allocatable :: wny     (:) ! norm unit vector y comp of sfc cells
     real, allocatable :: wnz     (:) ! norm unit vector z comp of sfc cells
     real, allocatable :: dzt_bot (:) ! surface similarity grid-height 

     ! Surface type (land/vegetation, lake, or sea)

     integer, allocatable :: leaf_class(:) ! leaf ("vegetation") class
     integer, allocatable :: ioge      (:) ! integer array for storing database data

     ! Atmospheric near-surface properties

     real, allocatable :: vels    (:) ! wind speed [m/s]
     real, allocatable :: prss    (:) ! surface air pressure [Pa]
     real, allocatable :: rhos    (:) ! air density [kg_air/m^3]
     real, allocatable :: airtemp (:) ! air temperature [K]
     real, allocatable :: airtheta(:) ! air potential temperature [K]
     real, allocatable :: airrrv  (:) ! air mixing ratio humidity [kg_vap/kg_dryair]
     real, allocatable :: airco2  (:) ! air mixing ratio of CO2 [kg_co2/kg_dryair]

     ! Canopy to atmosphere turbulent flux quantities

     real, allocatable :: ustar   (:) ! friction velocity [m/s]
     real, allocatable :: vkmsfc  (:) ! surface drag coefficient [kg/(m s)]
     real, allocatable :: sfluxt  (:)
     real, allocatable :: sfluxr  (:)
     real, allocatable :: sfluxc  (:)
     real, allocatable :: sxfer_t (:) ! can_air-to-atm heat xfer [kg_air K/m^2]
     real, allocatable :: sxfer_r (:) ! can_air-to-atm vapor xfer [kg_vap/m^2]
     real, allocatable :: sxfer_c (:) ! can_air-to-atm CO2 xfer [ppm/m^2]
     real, allocatable :: ggaer   (:) ! surface aerodynamic conductance [m/s]
     real, allocatable :: zeta    (:) ! surface z / M-O length [ ]
     real, allocatable :: wthv    (:) ! surface buoyancy flux [K m/s]

     ! Radiative flux quantities

     real, allocatable :: albedo_beam   (:) ! surface s/w beam albedo [0-1]
     real, allocatable :: albedo_diffuse(:) ! surface s/w diffuse albedo [0-1]
     real, allocatable :: rshort        (:) ! downward can-top s/w flux [W/m^2]
     real, allocatable :: rshort_diffuse(:) ! downward diffuse can-top s/w flux [W/m2]
     real, allocatable :: rshort_clr    (:) ! downward can-top s/w clr flux [W/m^2]
     real, allocatable :: rlong         (:) ! downward can-top l/w flux [W/m^2]
     real, allocatable :: rlong_albedo  (:) ! surface l/w lbedo [0-1]
     real, allocatable :: rlongup       (:) ! upward can-top l/w flux [W/m^2]

     ! Precipitation fluxes and runoff

     real, allocatable :: pcpg  (:) ! new pcp amount this timestep [kg/m^2]
     real, allocatable :: qpcpg (:) ! new pcp energy this timestep [J/m^2]
     real, allocatable :: dpcpg (:) ! new pcp depth this timestep [m]
     real, allocatable :: runoff(:) ! new runoff mass this timestep [kg/m^2]

     ! Canopy and surface quantities:

     real, allocatable :: can_depth(:) ! "canopy" depth for heat & vap capacity [m]
     real, allocatable :: cantemp  (:) ! "canopy" air temperature [K]
     real, allocatable :: canrrv   (:) ! "canopy" vapor mixing ratio [kg_vap/kg_dryair]
     real, allocatable :: rough    (:) ! roughness height [m]
     real, allocatable :: head1    (:) ! water surface hydraulic head (rel to topo) [m]

  End type surface_grid_vars

  type (surface_grid_vars), target :: sfcg

! TRI-GRID INFORMATION FOR INDEPENDENT REFINING OF SFC GRID

  integer :: nsfcgrids

  integer :: nsfcgrid_root

  integer :: nsfcgrdll(maxgrds)
  real    :: sfcgrdrad(maxgrds,maxngrdll)
  real    :: sfcgrdlat(maxgrds,maxngrdll)
  real    :: sfcgrdlon(maxgrds,maxngrdll)
  
  integer :: nmd
  integer :: nud
  integer :: nwd

  real, allocatable :: xemd(:), yemd(:), zemd(:)

  type (itab_md_vars), allocatable :: itab_md(:)
  type (itab_ud_vars), allocatable :: itab_ud(:)
  type (itab_wd_vars), allocatable :: itab_wd(:)

  type (itab_md_vars), allocatable :: ltab_md(:)
  type (itab_ud_vars), allocatable :: ltab_ud(:)
  type (itab_wd_vars), allocatable :: ltab_wd(:)

  type (nest_ud_vars), allocatable :: nest_ud(:)
  type (nest_wd_vars), allocatable :: nest_wd(:)
   
Contains

!=========================================================================

  subroutine alloc_sfcgrid1(mmsfc0, mvsfc0, mwsfc0)

  use misc_coms, only: rinit, runtype

  implicit none

  integer, intent(in) :: mmsfc0, mvsfc0, mwsfc0

  ! Allocate surface grid tables

  allocate (itab_msfc(mmsfc0))
  allocate (itab_vsfc(mvsfc0))
  allocate (itab_wsfc(mwsfc0))

  ! Allocate sfcg member arrays and initialize them with default values

  allocate (sfcg%xem  (mmsfc0)) ; sfcg%xem   = rinit
  allocate (sfcg%yem  (mmsfc0)) ; sfcg%yem   = rinit
  allocate (sfcg%zem  (mmsfc0)) ; sfcg%zem   = rinit
  allocate (sfcg%glatm(mmsfc0)) ; sfcg%glatm = rinit
  allocate (sfcg%glonm(mmsfc0)) ; sfcg%glonm = rinit
  allocate (sfcg%topm (mmsfc0)) ; sfcg%topm  = rinit

  allocate (sfcg%xev  (mvsfc0)) ; sfcg%xev   = rinit
  allocate (sfcg%yev  (mvsfc0)) ; sfcg%yev   = rinit
  allocate (sfcg%zev  (mvsfc0)) ; sfcg%zev   = rinit
  allocate (sfcg%dnu  (mvsfc0)) ; sfcg%dnu   = rinit
  allocate (sfcg%dniu (mvsfc0)) ; sfcg%dniu  = rinit
  allocate (sfcg%dnv  (mvsfc0)) ; sfcg%dnv   = rinit
  allocate (sfcg%dniv (mvsfc0)) ; sfcg%dniv  = rinit

  allocate (sfcg%area    (mwsfc0)) ; sfcg%area     = rinit
  allocate (sfcg%glatw   (mwsfc0)) ; sfcg%glatw    = rinit
  allocate (sfcg%glonw   (mwsfc0)) ; sfcg%glonw    = rinit
  allocate (sfcg%xew     (mwsfc0)) ; sfcg%xew      = rinit
  allocate (sfcg%yew     (mwsfc0)) ; sfcg%yew      = rinit
  allocate (sfcg%zew     (mwsfc0)) ; sfcg%zew      = rinit
  allocate (sfcg%topw    (mwsfc0)) ; sfcg%topw     = rinit
  allocate (sfcg%wnx     (mwsfc0)) ; sfcg%wnx      = rinit
  allocate (sfcg%wny     (mwsfc0)) ; sfcg%wny      = rinit
  allocate (sfcg%wnz     (mwsfc0)) ; sfcg%wnz      = rinit
  allocate (sfcg%dzt_bot (mwsfc0)) ; sfcg%dzt_bot  = rinit

  allocate (sfcg%leaf_class(mwsfc0)) ; sfcg%leaf_class = 0
  allocate (sfcg%ioge      (mwsfc0)) ; sfcg%ioge       = 0

  end subroutine alloc_sfcgrid1

!=========================================================================

  subroutine alloc_sfcgrid2(mwsfc)

     use misc_coms, only: rinit
     use mem_co2,   only: co2flag

     implicit none

     integer, intent(in) :: mwsfc

!    Allocate and initialize surface grid arrays

     allocate (sfcg%vels          (mwsfc)) ; sfcg%vels           = rinit
     allocate (sfcg%prss          (mwsfc)) ; sfcg%prss           = rinit
     allocate (sfcg%rhos          (mwsfc)) ; sfcg%rhos           = rinit
     allocate (sfcg%airtemp       (mwsfc)) ; sfcg%airtemp        = rinit
     allocate (sfcg%airtheta      (mwsfc)) ; sfcg%airtheta       = rinit
     allocate (sfcg%airrrv        (mwsfc)) ; sfcg%airrrv         = rinit
     allocate (sfcg%airco2        (mwsfc)) ; sfcg%airco2         = rinit
     allocate (sfcg%ustar         (mwsfc)) ; sfcg%ustar          = rinit
     allocate (sfcg%vkmsfc        (mwsfc)) ; sfcg%vkmsfc         = rinit
     allocate (sfcg%sfluxt        (mwsfc)) ; sfcg%sfluxt         = rinit
     allocate (sfcg%sfluxr        (mwsfc)) ; sfcg%sfluxr         = rinit
     allocate (sfcg%sxfer_t       (mwsfc)) ; sfcg%sxfer_t        = 0.0
     allocate (sfcg%sxfer_r       (mwsfc)) ; sfcg%sxfer_r        = 0.0
     if (co2flag /= 0) then
        allocate (sfcg%sfluxc        (mwsfc)) ; sfcg%sfluxc         = rinit
        allocate (sfcg%sxfer_c       (mwsfc)) ; sfcg%sxfer_c        = 0.0
     endif
     allocate (sfcg%ggaer         (mwsfc)) ; sfcg%ggaer          = rinit
     allocate (sfcg%wthv          (mwsfc)) ; sfcg%wthv           = rinit
     allocate (sfcg%albedo_beam   (mwsfc)) ; sfcg%albedo_beam    = 0.0
     allocate (sfcg%albedo_diffuse(mwsfc)) ; sfcg%albedo_diffuse = 0.0
     allocate (sfcg%rshort        (mwsfc)) ; sfcg%rshort         = 0.0
     allocate (sfcg%rshort_diffuse(mwsfc)) ; sfcg%rshort_diffuse = 0.0
     allocate (sfcg%rshort_clr    (mwsfc)) ; sfcg%rshort_clr     = 0.0
     allocate (sfcg%rlong         (mwsfc)) ; sfcg%rlong          = 0.0
     allocate (sfcg%rlong_albedo  (mwsfc)) ; sfcg%rlong_albedo   = 0.0
     allocate (sfcg%rlongup       (mwsfc)) ; sfcg%rlongup        = 0.0
     allocate (sfcg%pcpg          (mwsfc)) ; sfcg%pcpg           = 0.0
     allocate (sfcg%qpcpg         (mwsfc)) ; sfcg%qpcpg          = 0.0
     allocate (sfcg%dpcpg         (mwsfc)) ; sfcg%dpcpg          = 0.0
     allocate (sfcg%runoff        (mwsfc)) ; sfcg%runoff         = 0.0
     allocate (sfcg%can_depth     (mwsfc)) ; sfcg%can_depth      = rinit
     allocate (sfcg%cantemp       (mwsfc)) ; sfcg%cantemp        = rinit
     allocate (sfcg%canrrv        (mwsfc)) ; sfcg%canrrv         = rinit
     allocate (sfcg%rough         (mwsfc)) ; sfcg%rough          = rinit
     allocate (sfcg%head1         (mwsfc)) ; sfcg%head1          = rinit

  end subroutine alloc_sfcgrid2

!=========================================================================

  subroutine resize_sfcgrid(mmsfc1, mvsfc1, mwsfc1, iwnew)

  use misc_coms, only: rinit, runtype

  implicit none

  integer, intent(in) :: mmsfc1, mvsfc1, mwsfc1
  integer, optional, intent(in) :: iwnew(mwsfc1)

  integer :: mmsfc0, mvsfc0, mwsfc0, mmsfc2, mvsfc2, mwsfc2

  real,    allocatable :: rscr(:,:)
  integer, allocatable :: iscr(:,:)

  ! Resize surface grid member arrays

  mmsfc0 = size(sfcg%xem)
  mvsfc0 = size(sfcg%xev)
  mwsfc0 = size(sfcg%xew)

  mmsfc2 = min(mmsfc0,mmsfc1)
  mvsfc2 = min(mvsfc0,mvsfc1)
  mwsfc2 = min(mwsfc0,mwsfc1)

  allocate(rscr(mmsfc2,6))

  rscr(1:mmsfc2,1) = sfcg%xem  (1:mmsfc2)
  rscr(1:mmsfc2,2) = sfcg%yem  (1:mmsfc2)
  rscr(1:mmsfc2,3) = sfcg%zem  (1:mmsfc2)
  rscr(1:mmsfc2,4) = sfcg%glatm(1:mmsfc2)
  rscr(1:mmsfc2,5) = sfcg%glonm(1:mmsfc2)
  rscr(1:mmsfc2,6) = sfcg%topm (1:mmsfc2)

  deallocate (sfcg%xem)   ; allocate (sfcg%xem  (mmsfc1))
  deallocate (sfcg%yem)   ; allocate (sfcg%yem  (mmsfc1))
  deallocate (sfcg%zem)   ; allocate (sfcg%zem  (mmsfc1))
  deallocate (sfcg%glatm) ; allocate (sfcg%glatm(mmsfc1))
  deallocate (sfcg%glonm) ; allocate (sfcg%glonm(mmsfc1))
  deallocate (sfcg%topm)  ; allocate (sfcg%topm (mmsfc1))

  sfcg%xem  (1:mmsfc2) = rscr(1:mmsfc2,1) 
  sfcg%yem  (1:mmsfc2) = rscr(1:mmsfc2,2) 
  sfcg%zem  (1:mmsfc2) = rscr(1:mmsfc2,3) 
  sfcg%glatm(1:mmsfc2) = rscr(1:mmsfc2,4) 
  sfcg%glonm(1:mmsfc2) = rscr(1:mmsfc2,5) 
  sfcg%topm (1:mmsfc2) = rscr(1:mmsfc2,6) 

  deallocate(rscr) ; allocate(rscr(mvsfc2,7))

  rscr(1:mvsfc2,1) = sfcg%xev (1:mvsfc2)
  rscr(1:mvsfc2,2) = sfcg%yev (1:mvsfc2)
  rscr(1:mvsfc2,3) = sfcg%zev (1:mvsfc2)
  rscr(1:mvsfc2,4) = sfcg%dnu (1:mvsfc2)
  rscr(1:mvsfc2,5) = sfcg%dniu(1:mvsfc2)
  rscr(1:mvsfc2,6) = sfcg%dnv (1:mvsfc2)
  rscr(1:mvsfc2,7) = sfcg%dniv(1:mvsfc2)

  deallocate (sfcg%xev)  ; allocate (sfcg%xev(mvsfc1))
  deallocate (sfcg%yev)  ; allocate (sfcg%yev(mvsfc1))
  deallocate (sfcg%zev)  ; allocate (sfcg%zev(mvsfc1))
  deallocate (sfcg%dnu)  ; allocate (sfcg%dnu(mvsfc1))
  deallocate (sfcg%dniu) ; allocate (sfcg%dniu(mvsfc1))
  deallocate (sfcg%dnv)  ; allocate (sfcg%dnv(mvsfc1))
  deallocate (sfcg%dniv) ; allocate (sfcg%dniv(mvsfc1))

  sfcg%xev (1:mvsfc2) = rscr(1:mvsfc2,1) 
  sfcg%yev (1:mvsfc2) = rscr(1:mvsfc2,2) 
  sfcg%zev (1:mvsfc2) = rscr(1:mvsfc2,3) 
  sfcg%dnu (1:mvsfc2) = rscr(1:mvsfc2,4) 
  sfcg%dniu(1:mvsfc2) = rscr(1:mvsfc2,5) 
  sfcg%dnv (1:mvsfc2) = rscr(1:mvsfc2,6) 
  sfcg%dniv(1:mvsfc2) = rscr(1:mvsfc2,7) 

  deallocate(rscr) ; allocate(rscr(mwsfc2,11), iscr(mwsfc2,2))

  rscr(1:mwsfc2, 1) = sfcg%area      (1:mwsfc2)
  rscr(1:mwsfc2, 2) = sfcg%glatw     (1:mwsfc2)
  rscr(1:mwsfc2, 3) = sfcg%glonw     (1:mwsfc2)
  rscr(1:mwsfc2, 4) = sfcg%xew       (1:mwsfc2)
  rscr(1:mwsfc2, 5) = sfcg%yew       (1:mwsfc2)
  rscr(1:mwsfc2, 6) = sfcg%zew       (1:mwsfc2)
  rscr(1:mwsfc2, 7) = sfcg%topw      (1:mwsfc2)
  rscr(1:mwsfc2, 8) = sfcg%wnx       (1:mwsfc2)
  rscr(1:mwsfc2, 9) = sfcg%wny       (1:mwsfc2)
  rscr(1:mwsfc2,10) = sfcg%wnz       (1:mwsfc2)
  rscr(1:mwsfc2,11) = sfcg%dzt_bot   (1:mwsfc2)
  iscr(1:mwsfc2, 1) = sfcg%leaf_class(1:mwsfc2)
  iscr(1:mwsfc2, 2) = sfcg%ioge      (1:mwsfc2)

  deallocate (sfcg%area)       ; allocate (sfcg%area      (mwsfc1))
  deallocate (sfcg%glatw)      ; allocate (sfcg%glatw     (mwsfc1))
  deallocate (sfcg%glonw)      ; allocate (sfcg%glonw     (mwsfc1))
  deallocate (sfcg%xew)        ; allocate (sfcg%xew       (mwsfc1))
  deallocate (sfcg%yew)        ; allocate (sfcg%yew       (mwsfc1))
  deallocate (sfcg%zew)        ; allocate (sfcg%zew       (mwsfc1))
  deallocate (sfcg%topw)       ; allocate (sfcg%topw      (mwsfc1))
  deallocate (sfcg%wnx)        ; allocate (sfcg%wnx       (mwsfc1))
  deallocate (sfcg%wny)        ; allocate (sfcg%wny       (mwsfc1))
  deallocate (sfcg%wnz)        ; allocate (sfcg%wnz       (mwsfc1))
  deallocate (sfcg%dzt_bot)    ; allocate (sfcg%dzt_bot   (mwsfc1))
  deallocate (sfcg%leaf_class) ; allocate (sfcg%leaf_class(mwsfc1))
  deallocate (sfcg%ioge)       ; allocate (sfcg%ioge      (mwsfc1))

  if (present(iwnew)) then
     sfcg%area      (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 1)
     sfcg%glatw     (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 2) 
     sfcg%glonw     (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 3)
     sfcg%xew       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 4)
     sfcg%yew       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 5)
     sfcg%zew       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 6)
     sfcg%topw      (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 7)
     sfcg%wnx       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 8)
     sfcg%wny       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2, 9)
     sfcg%wnz       (iwnew(1:mwsfc2)) = rscr(1:mwsfc2,10)
     sfcg%dzt_bot   (iwnew(1:mwsfc2)) = rscr(1:mwsfc2,11)
     sfcg%leaf_class(iwnew(1:mwsfc2)) = iscr(1:mwsfc2, 1)
     sfcg%ioge      (iwnew(1:mwsfc2)) = iscr(1:mwsfc2, 2)
 else
     sfcg%area      (1:mwsfc2) = rscr(1:mwsfc2, 1)
     sfcg%glatw     (1:mwsfc2) = rscr(1:mwsfc2, 2) 
     sfcg%glonw     (1:mwsfc2) = rscr(1:mwsfc2, 3)
     sfcg%xew       (1:mwsfc2) = rscr(1:mwsfc2, 4)
     sfcg%yew       (1:mwsfc2) = rscr(1:mwsfc2, 5)
     sfcg%zew       (1:mwsfc2) = rscr(1:mwsfc2, 6)
     sfcg%topw      (1:mwsfc2) = rscr(1:mwsfc2, 7)
     sfcg%wnx       (1:mwsfc2) = rscr(1:mwsfc2, 8)
     sfcg%wny       (1:mwsfc2) = rscr(1:mwsfc2, 9)
     sfcg%wnz       (1:mwsfc2) = rscr(1:mwsfc2,10)
     sfcg%dzt_bot   (1:mwsfc2) = rscr(1:mwsfc2,11)
     sfcg%leaf_class(1:mwsfc2) = iscr(1:mwsfc2, 1)
     sfcg%ioge      (1:mwsfc2) = iscr(1:mwsfc2, 2)
  endif

  deallocate(rscr, iscr)

  end subroutine resize_sfcgrid

!=========================================================================

   subroutine filltab_sfcg()

     use var_tables, only: increment_vtable

     implicit none

     if (allocated(sfcg%vels))           call increment_vtable('SFCG%VELS',           'CW', rvar1=sfcg%vels)
     if (allocated(sfcg%prss))           call increment_vtable('SFCG%PRSS',           'CW', rvar1=sfcg%prss)
     if (allocated(sfcg%rhos))           call increment_vtable('SFCG%RHOS',           'CW', rvar1=sfcg%rhos)
     if (allocated(sfcg%airtemp))        call increment_vtable('SFCG%AIRTEMP',        'CW', rvar1=sfcg%airtemp)
     if (allocated(sfcg%airtheta))       call increment_vtable('SFCG%AIRTHETA',       'CW', rvar1=sfcg%airtheta)
     if (allocated(sfcg%airrrv))         call increment_vtable('SFCG%AIRRRV',         'CW', rvar1=sfcg%airrrv)
     if (allocated(sfcg%airco2))         call increment_vtable('SFCG%AIRCO2',         'CW', rvar1=sfcg%airco2)
     if (allocated(sfcg%ustar))          call increment_vtable('SFCG%USTAR',          'CW', rvar1=sfcg%ustar)
     if (allocated(sfcg%vkmsfc))         call increment_vtable('SFCG%VKMSFC',         'CW', rvar1=sfcg%vkmsfc)
     if (allocated(sfcg%sfluxt))         call increment_vtable('SFCG%SFLUXT',         'CW', rvar1=sfcg%sfluxt)
     if (allocated(sfcg%sfluxr))         call increment_vtable('SFCG%SFLUXR',         'CW', rvar1=sfcg%sfluxr)
     if (allocated(sfcg%sfluxc))         call increment_vtable('SFCG%SFLUXC',         'CW', rvar1=sfcg%sfluxc)
     if (allocated(sfcg%sxfer_t))        call increment_vtable('SFCG%SXFER_T',        'CW', rvar1=sfcg%sxfer_t)
     if (allocated(sfcg%sxfer_r))        call increment_vtable('SFCG%SXFER_R',        'CW', rvar1=sfcg%sxfer_r)
     if (allocated(sfcg%sxfer_c))        call increment_vtable('SFCG%SXFER_C',        'CW', rvar1=sfcg%sxfer_c)
     if (allocated(sfcg%ggaer))          call increment_vtable('SFCG%GGAER',          'CW', rvar1=sfcg%ggaer)
     if (allocated(sfcg%zeta))           call increment_vtable('SFCG%ZETA',           'CW', rvar1=sfcg%zeta)
     if (allocated(sfcg%wthv))           call increment_vtable('SFCG%WTHV',           'CW', rvar1=sfcg%wthv)
     if (allocated(sfcg%albedo_beam))    call increment_vtable('SFCG%ALBEDO_BEAM',    'CW', rvar1=sfcg%albedo_beam)
     if (allocated(sfcg%albedo_diffuse)) call increment_vtable('SFCG%ALBEDO_DIFFUSE', 'CW', rvar1=sfcg%albedo_diffuse)
     if (allocated(sfcg%rshort))         call increment_vtable('SFCG%RSHORT',         'CW', rvar1=sfcg%rshort)
     if (allocated(sfcg%rshort_diffuse)) call increment_vtable('SFCG%RSHORT_DIFFUSE', 'CW', rvar1=sfcg%rshort_diffuse)
     if (allocated(sfcg%rshort_clr))     call increment_vtable('SFCG%RSHORT_CLR',     'CW', rvar1=sfcg%rshort_clr)
     if (allocated(sfcg%rlong))          call increment_vtable('SFCG%RLONG',          'CW', rvar1=sfcg%rlong)
     if (allocated(sfcg%rlong_albedo))   call increment_vtable('SFCG%RLONG_ALBEDO',   'CW', rvar1=sfcg%rlong_albedo)
     if (allocated(sfcg%rlongup))        call increment_vtable('SFCG%RLONGUP',        'CW', rvar1=sfcg%rlongup)
     if (allocated(sfcg%pcpg))           call increment_vtable('SFCG%PCPG',           'CW', rvar1=sfcg%pcpg)
     if (allocated(sfcg%qpcpg))          call increment_vtable('SFCG%QPCPG',          'CW', rvar1=sfcg%qpcpg)
     if (allocated(sfcg%dpcpg))          call increment_vtable('SFCG%DPCPG',          'CW', rvar1=sfcg%dpcpg)
     if (allocated(sfcg%runoff))         call increment_vtable('SFCG%RUNOFF',         'CW', rvar1=sfcg%runoff)
     if (allocated(sfcg%can_depth))      call increment_vtable('SFCG%CAN_DEPTH',      'CW', rvar1=sfcg%can_depth)
     if (allocated(sfcg%cantemp))        call increment_vtable('SFCG%CANTEMP',        'CW', rvar1=sfcg%cantemp)
     if (allocated(sfcg%canrrv))         call increment_vtable('SFCG%CANRRV',         'CW', rvar1=sfcg%canrrv)
     if (allocated(sfcg%rough))          call increment_vtable('SFCG%ROUGH',          'CW', rvar1=sfcg%rough)
     if (allocated(sfcg%head1))          call increment_vtable('SFCG%HEAD1',          'CW', rvar1=sfcg%head1)
    
   end subroutine filltab_sfcg

!=========================================================================

   subroutine sfcg_avgatm()

     use mem_basic,   only: rho, press, theta, tair, rr_v, vxe, vye, vze
     use mem_micro,   only: rr_c
     use mem_co2,     only: rr_co2
     use misc_coms,   only: isubdomain
     use consts_coms, only: grav
     use mem_ijtabs,  only: itabg_w
     use mem_grid,    only: gdz_abov
     use mem_para,    only: myrank

     implicit none

     integer :: iwsfc, j, iw, kw
     real :: vels, psfc

     ! Average atmospheric quantities to each SFC grid cell location

     ! Loop over all SFC grid cells

     !$omp parallel do private(j,iw,kw,vels,psfc)
     do iwsfc = 2, mwsfc

        ! Skip this SFC grid cell if running in parallel and cell rank is not MYRANK
        if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        sfcg%vels    (iwsfc) = 0.
        sfcg%prss    (iwsfc) = 0.
        sfcg%rhos    (iwsfc) = 0.
        sfcg%airtemp (iwsfc) = 0.
        sfcg%airtheta(iwsfc) = 0.
        sfcg%airrrv  (iwsfc) = 0.
        sfcg%airco2  (iwsfc) = 0.

        ! Loop over all ATM grid cells that couple to this SFC grid cell

        do j = 1,itab_wsfc(iwsfc)%nwatm
           iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
           kw = itab_wsfc(iwsfc)%kwatm(j)

           vels = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
           psfc = press(kw,iw) + gdz_abov(kw) * rho(kw,iw) * grav  ! hydrostatic eqn.

           sfcg%vels    (iwsfc) = sfcg%vels    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * vels
           sfcg%prss    (iwsfc) = sfcg%prss    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * psfc
           sfcg%rhos    (iwsfc) = sfcg%rhos    (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rho  (kw,iw)
           sfcg%airtemp (iwsfc) = sfcg%airtemp (iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * tair (kw,iw)
           sfcg%airtheta(iwsfc) = sfcg%airtheta(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * theta(kw,iw)

           if (allocated(rr_c)) then
              sfcg%airrrv(iwsfc) = sfcg%airrrv(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * (rr_v(kw,iw) + rr_c(kw,iw))
           else
              sfcg%airrrv(iwsfc) = sfcg%airrrv(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rr_v(kw,iw)
           endif

           if (allocated(rr_co2)) then
              sfcg%airco2(iwsfc) = sfcg%airco2(iwsfc) + itab_wsfc(iwsfc)%arcoarsfc(j) * rr_co2(kw,iw)
           endif
        enddo

     enddo
     !$omp end parallel do

  end subroutine sfcg_avgatm

!===============================================================================

  subroutine fill_jwsfc()

  use mem_ijtabs, only: mrls

  use misc_coms,  only: io6, iparallel
   
  use mem_para,   only: mgroupsize, myrank,  &
                        send_wsfc, recv_wsfc,  &
                        nsends_wsfc, nrecvs_wsfc

  implicit none
   
  integer :: jsend, iwsfc, jend

  ! Allocate and zero-fill JTAB_WSFC_MPI%JEND

  do jsend = 1,maxremote
     allocate (jtab_wsfc_mpi(jsend)%jend(mrls))
               jtab_wsfc_mpi(jsend)%jend(1:mrls) = 0
  enddo

  ! Return if run is not parallel (jtab not needed)

  if (iparallel == 0) return

  ! Compute and store JTAB_WSFC_MPI%JEND(1)

  do jsend = 1,nsends_wsfc(1)
     jtab_wsfc_mpi(jsend)%jend(1) = 0
     do iwsfc = 2,mwsfc
        if (itab_wsfc(iwsfc)%send(jsend)) then
           jtab_wsfc_mpi(jsend)%jend(1) = jtab_wsfc_mpi(jsend)%jend(1) + 1
        endif
     enddo
     jtab_wsfc_mpi(jsend)%jend(1) = max(1,jtab_wsfc_mpi(jsend)%jend(1))
  enddo

  ! Allocate and zero-fill JTAB_WSFC_MPI%IWSFC

  do jsend = 1,nsends_wsfc(1)
     jend = jtab_wsfc_mpi(jsend)%jend(1)
     allocate (jtab_wsfc_mpi(jsend)%iwsfc(jend))
               jtab_wsfc_mpi(jsend)%iwsfc(1:jend) = 0
  enddo

  ! Initialize JTAB_WSFC_MPI%JEND counters to zero

  do jsend = 1,nsends_wsfc(1)
     jtab_wsfc_mpi(jsend)%jend(1:mrls) = 0
  enddo

  ! Compute JTAB_WSFC_MPI%IWSFC (only fill for MRL = 1)

  do iwsfc = 2,mwsfc
     do jsend = 1,nsends_wsfc(1)

        if (itab_wsfc(iwsfc)%send(jsend)) then
           jtab_wsfc_mpi(jsend)%jend(1) = jtab_wsfc_mpi(jsend)%jend(1) + 1
           jtab_wsfc_mpi(jsend)%iwsfc(jtab_wsfc_mpi(jsend)%jend(1)) = iwsfc
        endif

     enddo
  enddo

  end subroutine fill_jwsfc

End Module mem_sfcg
