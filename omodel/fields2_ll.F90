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
subroutine fields2_ll()

! This subroutine is a template, intended for user modification, for
! interpolating selected fields from the OLAM grid to a structured lat-lon
! grid of either limited or global area.  

! The following tasks are performed:

! (1) Allocate arrays for interpolated fields
! (2) Allocate scratch arrays
! (3) Copy each selected model field to scratch array, or compute it from
!     model fields if required
! (4) Call subroutine to interpolate each selected field to lat-lon array
! (5) Write lat-lon arrays to hdf5 file
! (6) Plot lat-lon arrays
! (7) Plot longitudinal averages of some lat-lon arrays

! New fields may be added following the examples in this file, and existing
! fields, or selected processes such as hdf5 output and plotting, may be
! disabled by commenting out lines of code.

  use mem_ijtabs,  only: itab_w, itabg_w
  use mem_basic,   only: vc, wc, rho, press, theta, sh_w, sh_v, &
                         vxe, vye, vze, tair
  use mem_grid,    only: mza, mwa, lpw, xew, yew, zew, topw, zt, &
                         wnx, wny, wnz, dzt
  use mem_leaf,    only: land
  use leaf_coms,   only: dslz
  use mem_sea,     only: sea
  use misc_coms,   only: io6, current_time, hfilepref, iclobber, iparallel, &
                         mdomain, time8
  use consts_coms, only: p00, rocp, piu180, erad, eradi, pio180, cp, alvl, &
                         rvap, r8
  use hdf5_utils,  only: shdf5_open, shdf5_orec, shdf5_orec_ll, shdf5_close, &
                         shdf5_write_global_attribute
  use max_dims,    only: pathlen
  use oname_coms,  only: nl
  use mem_para,    only: myrank, mgroupsize

#ifdef OLAM_MPI
  use mpi
#endif

! NOTE: The fields used from the MEM_MICRO, MEM_CUPARM, and MEM_FLUX_ACCUM
! modules are time-integrated fluxes of mass or energy and have units of
! [kg/m^2] or [J/m^2].  Their values written to any history file represent flux
! integrals from the beginning of a simulation to the time of the history file
! write.  Thus, time averages of a flux over any time interval can be computed
! as the difference between the flux integrals in two different history files
! divided by the difference in simulation times when the files were written.

  use mem_micro,  only: accpd, accpr, accpp, accps, accpa, accpg, accph

  use mem_cuparm, only: aconpr

  use mem_flux_accum, only:   rshort_accum,         rshortup_accum, &
                               rlong_accum,          rlongup_accum, &
                          rshort_top_accum,     rshortup_top_accum, &
                         rlongup_top_accum,                         &
                          rshort_clr_accum,     rshortup_clr_accum, &
                           rlong_clr_accum,      rlongup_clr_accum, &
                      rshort_top_clr_accum, rshortup_top_clr_accum, &
                     rlongup_top_clr_accum,                         &
                              sfluxt_accum,           sfluxr_accum, &
                                  vc_accum,               wc_accum, &
                               press_accum,             tair_accum, &
                                sh_v_accum,                         &
                              vels_l_accum,                         &
                           airtemp_l_accum,         airshv_l_accum, &
                           cantemp_l_accum,         canshv_l_accum, &
                          skintemp_l_accum,                         &
                            sfluxt_l_accum,         sfluxr_l_accum, &
                            wxfer1_l_accum,                         &
                              vels_s_accum,                         &
                           airtemp_s_accum,         airshv_s_accum, &
                           cantemp_s_accum,         canshv_s_accum, &
                          skintemp_s_accum,                         &
                            sfluxt_s_accum,         sfluxr_s_accum

  implicit none

!--------------------------------------------------------------------------------
! THE USER MUST SPECIFY THE FOLLOWING PARAMETERS THAT DEFINE THE OUTPUT
! LATITUDE-LONGITUDE ARRAYS. LONGITUDE VALUES ARE DEFINED TO BE IN THE
! RANGE (-180,180], WHICH DOES NOT PRECLUDE THAT THE LIMITED-AREA
! LATITUDE-LONGITUDE GRID CROSS THE 180 DEGREE MERIDIAN.
!--------------------------------------------------------------------------------

! Global domain:
  real   , parameter :: beglon = -180.   ! minimum (westernmost) longitude (deg)
  real   , parameter :: endlon =  180.   ! maximum (easternmost) longitude (deg)
  real   , parameter :: beglat =  -90.   ! minimum (southernmost) latitude (deg)
  real   , parameter :: endlat =   90.   ! maximum (northernmost) latitude (deg)
  integer, parameter :: rf = 1           ! horiz resolution factor (pts per deg)

! Regional domain:
! real   , parameter :: beglon = -92.   ! minimum (westernmost) longitude (deg)
! real   , parameter :: endlon = -72.   ! maximum (easternmost) longitude (deg)
! real   , parameter :: beglat =  23.   ! minimum (southernmost) latitude (deg)
! real   , parameter :: endlat =  35.   ! maximum (northernmost) latitude (deg)
! integer, parameter :: rf = 30         ! horiz resolution factor (pts per deg)

! DCMIP Global domain:
! real, parameter :: beglon = -179.5 ! minimum (westernmost) longitude (deg)
! real, parameter :: endlon =  179.5 ! maximum (easternmost) longitude (deg)
! real, parameter :: beglat =  -89.5 ! minimum (southernmost) latitude (deg)
! real, parameter :: endlat =   89.5 ! maximum (northernmost) latitude (deg)
! integer, parameter :: rf = 1       ! horiz resolution factor (pts per deg)

! NLON and NLAT are the number of longitude and latitude points in the lat-lon
! arrays that get filled by interpolation from the native OLAM grid.  They may
! be computed from values of BEGLON, ENDLON, BEGLAT, ENDLAT, and RF, or they
! may be specified directly.

  integer, save :: nlon, nlat

! Missing data indicators (mostly for underground points)

  integer,  parameter :: imissing = -999
  real,     parameter :: rmissing = -999.9
  real(r8), parameter :: dmissing = -999.9_r8

!--------------------------------------------------------------------------------
! THE FOLLOWING ARRAYS WILL CONTAIN THE FIELDS THAT ARE DEFINED ON THE 
! LATITUDE-LONGITUDE GRID.  EXCEPT FOR THE 1D VERTICAL ARRAYS, THEY ARE
! DIMENSIONED USING THE ABOVE NLON,NLAT PARAMETERS.  THE USER SHOULD ADD
! NEW ARRAYS AS REQUIRED.
!--------------------------------------------------------------------------------

!----------
! 1D FIELDS
!----------

  real, allocatable, save :: alat(:)  ! latitudes of grid points (deg)
  real, allocatable, save :: alon(:)  ! longitudes of grid points (deg)
  real, allocatable, save :: zlev(:)  ! heights above sea level of grid points (m)
  real, allocatable       :: value(:) ! longitudinal average of a field
  real, allocatable       :: arr_ll(:,:)

!----------
! 2D FIELDS
!----------
  
  real, allocatable :: topo_ll              (:) ! topography height (m)
  real, allocatable :: u_lpw_ll             (:) ! lpw zonal wind component (m/s)
  real, allocatable :: v_lpw_ll             (:) ! lpw meridional wind component (m/s)
  real, allocatable :: t_lpw_ll             (:) ! lpw air temperature (K)
  real, allocatable :: r_lpw_ll             (:) ! lpw water vapor specific density (kg/kg)
  real, allocatable :: pvap_lpw_ll          (:) ! lpw vapor pressure (Pa)
  real, allocatable :: slp_ll               (:) ! sea level pressure (Pa)

  real, allocatable :: accpmic_ll           (:) ! accum microphysics precip (kg/m^2)
  real, allocatable :: accpcon_ll           (:) ! accum convective precip (kg/m^2)
  real, allocatable :: vapflux_accum_ll     (:) ! accum sfc vapor (kg/m^2)
  real, allocatable :: sensflux_accum_ll    (:) ! accum sfc sensible heat (J/m^2)
  real, allocatable :: latflux_accum_ll     (:) ! accum sfc latent heat (J/m^2)
  real, allocatable :: rshort_accum_ll      (:) ! accum sfc downward s/w rad (J/m^2)
  real, allocatable :: rshortup_accum_ll    (:) ! accum sfc upward s/w rad (J/m^2)
  real, allocatable :: rlong_accum_ll       (:) ! accum sfc downward l/w rad (J/m^2)
  real, allocatable :: rlongup_accum_ll     (:) ! accum sfc upward l/w rad (J/m^2)
  real, allocatable :: rshort_top_accum_ll  (:) ! accum TOA downward s/w rad (J/m^2)
  real, allocatable :: rshortup_top_accum_ll(:) ! accum TOA upward s/w rad (J/m^2)
  real, allocatable :: rlongup_top_accum_ll (:) ! accum TOA upward l/w rad (J/m^2)

  real, allocatable :: rshort_clr_accum_ll      (:) ! accum sfc downward s/w rad (J/m^2)
  real, allocatable :: rshortup_clr_accum_ll    (:) ! accum sfc upward s/w rad (J/m^2)
  real, allocatable :: rlong_clr_accum_ll       (:) ! accum sfc downward l/w rad (J/m^2)
  real, allocatable :: rlongup_clr_accum_ll     (:) ! accum sfc upward l/w rad (J/m^2)
  real, allocatable :: rshort_top_clr_accum_ll  (:) ! accum TOA downward s/w rad (J/m^2)
  real, allocatable :: rshortup_top_clr_accum_ll(:) ! accum TOA upward s/w rad (J/m^2)
  real, allocatable :: rlongup_top_clr_accum_ll (:) ! accum TOA upward l/w rad (J/m^2)

  real, allocatable :: als_vels_accum_ll     (:) ! als wind speed accum (m)
  real, allocatable :: als_airtempk_accum_ll (:) ! als atm temp accum (K s)
  real, allocatable :: als_airshv_accum_ll   (:) ! als atm shv accum (kg s/kg)
  real, allocatable :: als_cantempk_accum_ll (:) ! als canopy air temp accum (K s)
  real, allocatable :: als_canshv_accum_ll   (:) ! als canopy air vap spec dens accum (kg s/kg)
  real, allocatable :: als_skintempk_accum_ll(:) ! als canopy air temp accum (K s)
  real, allocatable :: als_sensflux_accum_ll (:) ! als sens heat flux accum (J/m^2)
  real, allocatable :: als_latflux_accum_ll  (:) ! als lat heat flux accum (J/m^2)
  real, allocatable :: als_vapflux_accum_ll  (:) ! als vap flux accum (kg/m^2)
  real, allocatable :: al_wxfer1_accum_ll    (:) ! al soil bottom water flux accum (m)

  real, allocatable :: al_sfcwater_tot_ll  (:) ! al total sfcwater mass (kg/m^2)
  real, allocatable :: al_soil_water_tot_ll(:) ! al total soil water (m)

  real, allocatable, save :: pcpmic_dif2_ll      (:)
  real, allocatable, save :: pcpcon_dif2_ll      (:)
  real, allocatable, save :: pcpboth_dif2_ll     (:)
  real, allocatable, save :: rshort_dif2_ll      (:)
  real, allocatable, save :: rshortup_dif2_ll    (:)
  real, allocatable, save :: rlong_dif2_ll       (:)
  real, allocatable, save :: rlongup_dif2_ll     (:)
  real, allocatable, save :: rshort_top_dif2_ll  (:)
  real, allocatable, save :: rshortup_top_dif2_ll(:)
  real, allocatable, save :: rlongup_top_dif2_ll (:)
  real, allocatable, save :: sensflux_dif2_ll    (:)
  real, allocatable, save :: latflux_dif2_ll     (:)

!----------
! 3D FIELDS
!----------

  real, allocatable :: u_ll (:,:) ! zonal wind component (m/s)
  real, allocatable :: v_ll (:,:) ! meridional wind component (m/s)
  real, allocatable :: w_ll (:,:) ! vertical wind component (m/s)
  real, allocatable :: t_ll (:,:) ! air temperature (K)
  real, allocatable :: r_ll (:,:) ! water vapor specific density (kg/kg)
  real, allocatable :: p_ll (:,:) ! air pressure (Pa)

  real, allocatable :: u_accum_ll (:,:) ! zonal wind component accum (m/s)
  real, allocatable :: v_accum_ll (:,:) ! meridional wind component accum (m/s)
  real, allocatable :: w_accum_ll (:,:) ! vertical wind component accum (m/s)
  real, allocatable :: t_accum_ll (:,:) ! air temperature accum (K)
  real, allocatable :: r_accum_ll (:,:) ! water vapor specific density accum (kg/kg)
  real, allocatable :: p_accum_ll (:,:) ! air pressure accum (Pa)

  integer       :: k, iw, ilat, ilon, n, kb, ier, np
  integer       :: iland, jland, isea, jsea, iv, jv, npoly, klev
  integer       :: ndims, idims(3)
  character(30) :: dimnames(3)

  real :: raxis, raxisi, area_land_sum, area_sea_sum
  real :: vx, vy, vz

  real :: scr1a(mwa), scr1b(mwa), scr1c(mwa), scr1d(mwa)
  real :: scr1e(mwa), scr1f(mwa), scr1g(mwa), scr1h(mwa)
  real :: scr1i(mwa), scr1j(mwa), scr1k(mwa), scr1l(mwa)

  real :: scr2a(mza,mwa), scr2b(mza,mwa)
  real, allocatable :: scr2_ll(:), scr3_ll(:,:)

  real :: timefac
  real :: aspect, scalelab, ymin, ymax, yinc, alatinc

! SEIGEL 2013 - Added for ll interp writeout
  character(pathlen) :: hnamel
  character(20)      :: ofrq

  integer,  save :: npts
  integer,  save :: init = 0
  real,     save :: dlon, dlat
  real(r8), save :: time8_prev

  integer, allocatable :: iws_l (:)
  integer, allocatable :: iws_gl(:,:)
  integer, allocatable :: iws_ll(:,:,:)
  real,    allocatable :: wts_ll(:,:,:)

  integer, allocatable, save :: iws_loc(:,:), lls_loc(:)
  real,    allocatable, save :: wts_loc(:,:)

  ! These switches control which lat/lon interpolations are performed

  logical, parameter :: dosfc   = .true.
  logical, parameter :: doaccum = .false.
  logical, parameter :: dodifs  = .false.
  logical, parameter :: do3d    = .false.
  logical, parameter :: dopress = .false. ! TODO!

!------------------------------------------------------
! 3D FIELDS - interpolation to pressure levels (TODO!)
!------------------------------------------------------

  integer, parameter :: npress = 8
  real :: plev(npress) = (/ 1000., 925., 850., 700., 500., 300., 200., 100. /)

  if (nl%ioutput_latlon /= 1 .and. nl%latlonplot /= 1) return

  if (.not. (dosfc .or. doaccum .or. dodifs .or. do3d .or. dopress)) return

  if (init == 0) then

     init = 1

     if (myrank == 0) write(io6,'(/,a)') "Initializing overlaps for lat/lon outputs..."

     ! Compute latitude and longitude of output grid points (assuming uniform spacing)

     nlat = nint((endlat - beglat) * real(rf)) + 1  ! # of lat values
     if (nlat < 2) stop 'stop: nlat < 2 in fields2_ll '
     dlat = (endlat  - beglat) / real(nlat-1)

     allocate(alat(nlat))

     do ilat = 1, nlat
        alat(ilat) = beglat + dlat * real(ilat-1)
     enddo

     if (endlon > beglon) then
        nlon = nint((endlon - beglon) * real(rf)) + 1  ! # of lon values
        if (nlon < 2) stop 'stop: nlon < 2 in fields2_ll '
        dlon = (endlon - beglon) / real(nlon-1)
     else
        nlon = nint((endlon + 360. - beglon) * rf) + 1  ! # of lon values
        if (nlon < 2) stop 'stop: nlon < 2 in fields2_ll '
        dlon = (endlon + 360. - beglon) / real(nlon-1)
     endif

     allocate(alon(nlon))

     do ilon = 1, nlon
        alon(ilon) = beglon + dlon * real(ilon-1)
        if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.
     enddo

     allocate(iws_ll(nlon,nlat,3)) ; iws_ll = 1
     allocate(wts_ll(nlon,nlat,3)) ; wts_ll = 0.0

     ! Allocate and fill zlev with model grid levels

     allocate(zlev(mza-1))

     do k = 2, mza
        zlev(k-1) = zt(k)
     enddo

     ! Find the 3 IW points and weights for interpolation to each lat/lon point

     call find_3iws_ll(nlon,nlat,alon,alat,iws_ll,wts_ll)

     ! Compute number of lat/lon points

     if (iparallel == 0) then

        npts = nlon*nlat

     else

     ! In parallel run, subroutine find_3iws_ll may have duplicated a few
     ! lat/lon interpolation points across different nodes.  Sort this out
     ! in order to eliminate duplicates.

        allocate (iws_l(nlon),iws_gl(nlon,mgroupsize))

        do ilat = 1,nlat

           iws_l(:) = itab_w( iws_ll(:,ilat,1) )%iwglobe

#ifdef OLAM_MPI
           call MPI_Allgather(iws_l, nlon, MPI_INTEGER, iws_gl, nlon, MPI_INTEGER, MPI_COMM_WORLD, ier)
#endif

           do ilon = 1, nlon

              ! In case of multiple matches in parallel, select the cell with
              ! the lowest global rank to match the single-processor result

              if ( count(iws_gl(ilon,:)>1) > 1 ) then
                 np = minloc(iws_gl(ilon,:), mask=iws_gl(ilon,:)>1, dim=1) - 1
                 if (myrank /= np) then
                    iws_ll(ilon,ilat,:) = 1
                    wts_ll(ilon,ilat,:) = 0.0
                 endif
              endif

           enddo
        enddo

        npts = count(iws_ll(:,:,1) > 1)

     endif

     ! Store node-local copies of the interpolation IW indices and weights

     allocate(iws_loc(npts,3)); iws_loc = 0  ! local iw points for interpolating each lat/lon point
     allocate(wts_loc(npts,3)); wts_loc = 0. ! local weights for interpolating each lat/lon point
     allocate(lls_loc(npts))  ; lls_loc = 0  ! lat/lon indices on this node

     n = 0
     do ilat = 1, nlat
        do ilon = 1, nlon
           if (iws_ll(ilon,ilat,1) > 1) then
              n = n + 1
              iws_loc(n,1:3) = iws_ll(ilon,ilat,1:3) 
              wts_loc(n,1:3) = wts_ll(ilon,ilat,1:3) 
              lls_loc(n)     = ilon + (ilat-1)*nlon
           endif
        enddo
     enddo

     deallocate(iws_ll)
     deallocate(wts_ll)

     if (dodifs) then
        allocate( pcpmic_dif2_ll      (npts) ) ; pcpmic_dif2_ll       = 0.
        allocate( pcpcon_dif2_ll      (npts) ) ; pcpcon_dif2_ll       = 0.
        allocate( pcpboth_dif2_ll     (npts) ) ; pcpboth_dif2_ll      = 0.
        allocate( rshort_dif2_ll      (npts) ) ; rshort_dif2_ll       = 0.
        allocate( rshortup_dif2_ll    (npts) ) ; rshortup_dif2_ll     = 0.
        allocate( rlong_dif2_ll       (npts) ) ; rlong_dif2_ll        = 0.
        allocate( rlongup_dif2_ll     (npts) ) ; rlongup_dif2_ll      = 0.
        allocate( rshort_top_dif2_ll  (npts) ) ; rshort_top_dif2_ll   = 0.
        allocate( rshortup_top_dif2_ll(npts) ) ; rshortup_top_dif2_ll = 0.
        allocate( rlongup_top_dif2_ll (npts) ) ; rlongup_top_dif2_ll  = 0.
        allocate( sensflux_dif2_ll    (npts) ) ; sensflux_dif2_ll     = 0.
        allocate( latflux_dif2_ll     (npts) ) ; latflux_dif2_ll      = 0.
     endif

  endif

! Initialize latitude-longitude arrays to zero prior to interpolation.
! For certain applications, zero should be replaced with "missing value".

  if (dosfc) then
     allocate( topo_ll              (npts) ) ; topo_ll               = rmissing
     allocate( u_lpw_ll             (npts) ) ; u_lpw_ll              = rmissing
     allocate( v_lpw_ll             (npts) ) ; v_lpw_ll              = rmissing
     allocate( t_lpw_ll             (npts) ) ; t_lpw_ll              = rmissing
     allocate( r_lpw_ll             (npts) ) ; r_lpw_ll              = rmissing
     allocate( pvap_lpw_ll          (npts) ) ; pvap_lpw_ll           = rmissing
     allocate( slp_ll               (npts) ) ; slp_ll                = rmissing
  endif

  if (doaccum .or. dodifs) then
     allocate( accpmic_ll           (npts) ) ; accpmic_ll            = rmissing
     allocate( accpcon_ll           (npts) ) ; accpcon_ll            = rmissing
     allocate( vapflux_accum_ll     (npts) ) ; vapflux_accum_ll      = rmissing
     allocate( sensflux_accum_ll    (npts) ) ; sensflux_accum_ll     = rmissing
     allocate( latflux_accum_ll     (npts) ) ; latflux_accum_ll      = rmissing
     allocate( rshort_accum_ll      (npts) ) ; rshort_accum_ll       = rmissing
     allocate( rshortup_accum_ll    (npts) ) ; rshortup_accum_ll     = rmissing
     allocate( rlong_accum_ll       (npts) ) ; rlong_accum_ll        = rmissing
     allocate( rlongup_accum_ll     (npts) ) ; rlongup_accum_ll      = rmissing
     allocate( rshort_top_accum_ll  (npts) ) ; rshort_top_accum_ll   = rmissing
     allocate( rshortup_top_accum_ll(npts) ) ; rshortup_top_accum_ll = rmissing
     allocate( rlongup_top_accum_ll (npts) ) ; rlongup_top_accum_ll  = rmissing
  endif

  if (doaccum) then
     allocate( rshort_clr_accum_ll      (npts) ) ; rshort_clr_accum_ll       = rmissing
     allocate( rshortup_clr_accum_ll    (npts) ) ; rshortup_clr_accum_ll     = rmissing
     allocate( rlong_clr_accum_ll       (npts) ) ; rlong_clr_accum_ll        = rmissing
     allocate( rlongup_clr_accum_ll     (npts) ) ; rlongup_clr_accum_ll      = rmissing
     allocate( rshort_top_clr_accum_ll  (npts) ) ; rshort_top_clr_accum_ll   = rmissing
     allocate( rshortup_top_clr_accum_ll(npts) ) ; rshortup_top_clr_accum_ll = rmissing
     allocate( rlongup_top_clr_accum_ll (npts) ) ; rlongup_top_clr_accum_ll  = rmissing

     allocate( al_sfcwater_tot_ll  (npts) ) ; al_sfcwater_tot_ll   = rmissing
     allocate( al_soil_water_tot_ll(npts) ) ; al_soil_water_tot_ll = rmissing

     allocate( als_vels_accum_ll     (npts) ) ; als_vels_accum_ll      = rmissing
     allocate( als_airtempk_accum_ll (npts) ) ; als_airtempk_accum_ll  = rmissing
     allocate( als_airshv_accum_ll   (npts) ) ; als_airshv_accum_ll    = rmissing
     allocate( als_cantempk_accum_ll (npts) ) ; als_cantempk_accum_ll  = rmissing
     allocate( als_canshv_accum_ll   (npts) ) ; als_canshv_accum_ll    = rmissing
     allocate( als_skintempk_accum_ll(npts) ) ; als_skintempk_accum_ll = rmissing
     allocate( als_sensflux_accum_ll (npts) ) ; als_sensflux_accum_ll  = rmissing
     allocate( als_latflux_accum_ll  (npts) ) ; als_latflux_accum_ll   = rmissing
     allocate( als_vapflux_accum_ll  (npts) ) ; als_vapflux_accum_ll   = rmissing
     allocate( al_wxfer1_accum_ll    (npts) ) ; al_wxfer1_accum_ll     = rmissing
  endif

  if (doaccum .and. do3d) then
     allocate( u_accum_ll(npts,mza-1) ) ; u_accum_ll = rmissing
     allocate( v_accum_ll(npts,mza-1) ) ; v_accum_ll = rmissing
     allocate( w_accum_ll(npts,mza-1) ) ; w_accum_ll = rmissing
     allocate( t_accum_ll(npts,mza-1) ) ; t_accum_ll = rmissing
     allocate( r_accum_ll(npts,mza-1) ) ; r_accum_ll = rmissing
     allocate( p_accum_ll(npts,mza-1) ) ; p_accum_ll = rmissing
  endif

  if (do3d) then
     allocate( u_ll(npts,mza-1) ) ; u_ll = rmissing
     allocate( v_ll(npts,mza-1) ) ; v_ll = rmissing
     allocate( w_ll(npts,mza-1) ) ; w_ll = rmissing
     allocate( t_ll(npts,mza-1) ) ; t_ll = rmissing
     allocate( r_ll(npts,mza-1) ) ; r_ll = rmissing
     allocate( p_ll(npts,mza-1) ) ; p_ll = rmissing
  endif

  ! scratch arrays
  allocate( scr2_ll(npts) )
  allocate( scr3_ll(npts,mza-1) )

  if (myrank == 0) write(io6,'(/,a)') "Interpolating fields to lat/lon grid..."

!------------------------------------------------------------
! Interpolate topography height
!------------------------------------------------------------

  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,topw,topo_ll)
  
!------------------------------------------------------------
! Compute zonal and meridional wind components on OLAM grid
! and copy their values at lowest prognosed model level to
! separate arrays
!------------------------------------------------------------

  do iw = 2, mwa
     kb = lpw(iw)

     if (mdomain < 2) then
        raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

        do k = kb, mza
           if (raxis > 1.e3) then
              scr2a(k,iw) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) / raxis
              scr2b(k,iw) = vze(k,iw) * raxis * eradi &
                          - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) / (raxis * erad)
           else
              scr2a(k,iw) = 0.
              scr2b(k,iw) = 0.
           endif
        enddo

     else

        do k = kb, mza
           scr2a(k,iw) = vxe(k,iw)
           scr2b(k,iw) = vye(k,iw)
        enddo

     endif

     scr1a(iw) = scr2a(kb,iw)
     scr1b(iw) = scr2b(kb,iw)

  enddo

!------------------------------------------------------------
! Interpolate zonal and meridional wind components, and their
! near-surface values
!------------------------------------------------------------

  if (do3d)  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,u_ll)
  if (do3d)  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,v_ll)

  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,u_lpw_ll)
  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,v_lpw_ll)

!------------------------------------------------------------
! Compute temperature on OLAM grid and copy its value at 
! lowest prognosed model level to separate array.
! Compute sea level pressure based on values at lowest
! prognosed model level. 
!------------------------------------------------------------

     do iw = 2, mwa
        kb = lpw(iw)
        scr1a(iw) = tair(kb,iw)
        scr1b(iw) = press(kb,iw) &
                  * (1. - .0065 * zt(kb) / (tair(kb,iw) + .0065 * zt(kb)))**(-5.257)
     enddo

!------------------------------------------------------------
! Interpolate temperature, its surface value, and sea level pressure
!------------------------------------------------------------

  if (do3d)  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,tair,t_ll)
  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,t_lpw_ll)
  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,slp_ll)

!------------------------------------------------------------
! Copy vapor specific density at lowest prognosed model level
! on OLAM grid to separate array.
! Compute vertical velocity at OLAM T point.
! Compute vapor pressure at lowest prognosed model level on OLAM grid.
!------------------------------------------------------------

  do iw = 2,mwa
     kb = lpw(iw)

     do k = kb,mza
        scr2a(k,iw) = 0.5 * (wc(k,iw) + wc(k-1,iw))
     enddo

     scr1a(iw) = sh_v(kb,iw)
     scr1b(iw) = sh_v(kb,iw) * rho(kb,iw) * rvap * tair(kb,iw)
  enddo

!------------------------------------------------------------
! Interpolate vapor specific density ratio, its surface value,
! and the surface value of vapor pressure
!------------------------------------------------------------

  if (do3d)  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,w_ll)
  if (do3d)  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,sh_v,r_ll)
  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,r_lpw_ll)
  if (dosfc) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,pvap_lpw_ll)

!------------------------------------------------------------
! Interpolate atmospheric pressure
!------------------------------------------------------------

! Copy pressure to real array

  scr2a = 0.0
  do iw = 2,mwa
     scr2a(lpw(iw):mza,iw) = real(press(lpw(iw):mza,iw))
  enddo

  if (do3d) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,p_ll)

!------------------------------------------------------------
! Compute total (resolved + parameterized) accumulated
! precipitation on OLAM grid
!------------------------------------------------------------

  if (doaccum .or. dodifs) then

     do iw = 2, mwa
        scr1a(iw) = 0.

        if (allocated(accpd))  scr1a(iw) = scr1a(iw) + real(accpd(iw))
        if (allocated(accpr))  scr1a(iw) = scr1a(iw) + real(accpr(iw))
        if (allocated(accpp))  scr1a(iw) = scr1a(iw) + real(accpp(iw))
        if (allocated(accps))  scr1a(iw) = scr1a(iw) + real(accps(iw))
        if (allocated(accpa))  scr1a(iw) = scr1a(iw) + real(accpa(iw))
        if (allocated(accpg))  scr1a(iw) = scr1a(iw) + real(accpg(iw))
        if (allocated(accph))  scr1a(iw) = scr1a(iw) + real(accph(iw))
     enddo

!-------------------------------------------------------------
! Compute and interpolate accumulated precipitation and fluxes
!-------------------------------------------------------------

     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,accpmic_ll)

     if (allocated(aconpr)) then
        scr1a(:) = real(aconpr(:))
        call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,accpcon_ll)
     endif

     scr1a(:) = real(sfluxr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,vapflux_accum_ll)

     scr1a(:) = real(sfluxt_accum(:)) * cp
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,sensflux_accum_ll)

     scr1a(:) = real(sfluxr_accum(:)) * alvl
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,latflux_accum_ll)

     scr1a(:) = real(rshort_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshort_accum_ll)

     scr1a(:) = real(rshortup_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshortup_accum_ll)

     scr1a(:) = real(rlong_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlong_accum_ll)

     scr1a(:) = real(rlongup_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlongup_accum_ll)

     scr1a(:) = real(rshort_top_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshort_top_accum_ll)

     scr1a(:) = real(rshortup_top_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshortup_top_accum_ll)

     scr1a(:) = real(rlongup_top_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlongup_top_accum_ll)

  endif

  if (doaccum) then

     scr1a(:) = real(rshort_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshort_clr_accum_ll)

     scr1a(:) = real(rshortup_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshortup_clr_accum_ll)

     scr1a(:) = real(rlong_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlong_clr_accum_ll)

     scr1a(:) = real(rlongup_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlongup_clr_accum_ll)

     scr1a(:) = real(rshort_top_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshort_top_clr_accum_ll)

     scr1a(:) = real(rshortup_top_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rshortup_top_clr_accum_ll)

     scr1a(:) = real(rlongup_top_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlongup_top_clr_accum_ll)

     scr1a(:) = real(rlongup_top_clr_accum(:))
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,rlongup_top_clr_accum_ll)

     do iw = 2, mwa
        scr1a(iw) = 0.
        scr1b(iw) = 0.
        scr1c(iw) = 0.
        scr1d(iw) = 0.
        scr1e(iw) = 0.
        scr1f(iw) = 0.
        scr1g(iw) = 0.
        scr1h(iw) = 0.
        scr1i(iw) = 0.
        scr1j(iw) = 0.
        scr1k(iw) = 0.
        scr1l(iw) = 0.

        area_land_sum = 0.
        area_sea_sum  = 0.

        do jland = 1,itab_w(iw)%nland
           iland = itab_w(iw)%iland(jland)

           scr1a(iw) = scr1a(iw) + real(    vels_l_accum(iland)) * land%area(iland)
           scr1b(iw) = scr1b(iw) + real( airtemp_l_accum(iland)) * land%area(iland)
           scr1c(iw) = scr1c(iw) + real(  airshv_l_accum(iland)) * land%area(iland)
           scr1d(iw) = scr1d(iw) + real( cantemp_l_accum(iland)) * land%area(iland)
           scr1e(iw) = scr1e(iw) + real(  canshv_l_accum(iland)) * land%area(iland)
           scr1f(iw) = scr1f(iw) + real(skintemp_l_accum(iland)) * land%area(iland)
           scr1g(iw) = scr1g(iw) + real(  sfluxt_l_accum(iland)) * land%area(iland) * cp
           scr1h(iw) = scr1h(iw) + real(  sfluxr_l_accum(iland)) * land%area(iland) * alvl
           scr1i(iw) = scr1i(iw) + real(  sfluxr_l_accum(iland)) * land%area(iland)
           scr1j(iw) = scr1j(iw) + real(  wxfer1_l_accum(iland)) * land%area(iland)
           scr1k(iw) = scr1k(iw) &
                     + sum(land%sfcwater_mass(:,iland))          * land%area(iland)
           scr1l(iw) = scr1l(iw) &
                     + sum(land%soil_water(:,iland) * dslz(:))   * land%area(iland)

           area_land_sum = area_land_sum + land%area(iland)
        enddo

        do jsea = 1,itab_w(iw)%nsea
           isea = itab_w(iw)%isea(jsea)

           scr1a(iw) = scr1a(iw) + real(    vels_s_accum(isea)) * sea%area(isea)
           scr1b(iw) = scr1b(iw) + real( airtemp_s_accum(isea)) * sea%area(isea)
           scr1c(iw) = scr1c(iw) + real(  airshv_s_accum(isea)) * sea%area(isea)
           scr1d(iw) = scr1d(iw) + real( cantemp_s_accum(isea)) * sea%area(isea)
           scr1e(iw) = scr1e(iw) + real(  canshv_s_accum(isea)) * sea%area(isea)
           scr1f(iw) = scr1f(iw) + real(skintemp_s_accum(isea)) * sea%area(isea)
           scr1g(iw) = scr1g(iw) + real(  sfluxt_s_accum(isea)) * sea%area(isea) * cp
           scr1h(iw) = scr1h(iw) + real(  sfluxr_s_accum(isea)) * sea%area(isea) * alvl
           scr1i(iw) = scr1i(iw) + real(  sfluxr_s_accum(isea)) * sea%area(isea)

           area_sea_sum = area_sea_sum + sea%area(isea)
        enddo

        if (area_land_sum + area_sea_sum > 1.e0) then
           scr1a(iw) = scr1a(iw) / (area_land_sum + area_sea_sum)
           scr1b(iw) = scr1b(iw) / (area_land_sum + area_sea_sum)
           scr1c(iw) = scr1c(iw) / (area_land_sum + area_sea_sum)
           scr1d(iw) = scr1d(iw) / (area_land_sum + area_sea_sum)
           scr1e(iw) = scr1e(iw) / (area_land_sum + area_sea_sum)
           scr1f(iw) = scr1f(iw) / (area_land_sum + area_sea_sum)
           scr1g(iw) = scr1g(iw) / (area_land_sum + area_sea_sum)
           scr1h(iw) = scr1h(iw) / (area_land_sum + area_sea_sum)
           scr1i(iw) = scr1i(iw) / (area_land_sum + area_sea_sum)
        endif

        if (area_land_sum > 1.e0) then
           scr1j(iw) = scr1j(iw) /  area_land_sum
           scr1k(iw) = scr1k(iw) /  area_land_sum
           scr1l(iw) = scr1l(iw) /  area_land_sum
        endif

     enddo

     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,als_vels_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,als_airtempk_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1c,als_airshv_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1d,als_cantempk_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1e,als_canshv_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1f,als_skintempk_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1g,als_sensflux_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1h,als_latflux_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1i,als_vapflux_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1j,al_wxfer1_accum_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1k,al_sfcwater_tot_ll)
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1l,al_soil_water_tot_ll)

     if (do3d) then

        do iw = 2, mwa
           npoly = itab_w(iw)%npoly

           do k = lpw(iw), mza

              vx = 0.5 * real(wc_accum(k-1,iw) + wc_accum(k,iw)) * wnx(iw)
              vy = 0.5 * real(wc_accum(k-1,iw) + wc_accum(k,iw)) * wny(iw)
              vz = 0.5 * real(wc_accum(k-1,iw) + wc_accum(k,iw)) * wnz(iw)

              do jv = 1, npoly
                 iv = itab_w(iw)%iv(jv)

                 vx = vx + itab_w(iw)%ecvec_vx(jv) * real(vc_accum(k,iv))
                 vy = vy + itab_w(iw)%ecvec_vy(jv) * real(vc_accum(k,iv))
                 vz = vz + itab_w(iw)%ecvec_vz(jv) * real(vc_accum(k,iv))
              enddo

              if (mdomain < 2) then

                 raxis = sqrt(xew(iw)**2 + yew(iw)**2)  ! dist from earth axis

                 if (raxis > 1.e3) then
                    scr2a(k,iw) = (vy * xew(iw) - vx * yew(iw)) / raxis
                    scr2b(k,iw) = vz * raxis / erad &
                       - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 
                 else
                    scr2a(k,iw) = 0.
                    scr2b(k,iw) = 0.
                 endif

              else
                 scr2a(k,iw) = vx
                 scr2b(k,iw) = vy
              endif

           enddo
        enddo

        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,u_accum_ll)
        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,v_accum_ll)

        do iw = 2,mwa
           kb = lpw(iw)

           do k = kb,mza
              scr2a(k,iw) = 0.5 * real(wc_accum(k,iw) + wc_accum(k-1,iw))
              scr2b(k,iw) = real(press_accum(k,iw))
           enddo
        enddo

        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,w_accum_ll)
        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,p_accum_ll)

        do iw = 2,mwa
           kb = lpw(iw)

           do k = kb,mza
              scr2a(k,iw) = real(tair_accum(k,iw))
              scr2b(k,iw) = real(sh_v_accum(k,iw))
           enddo
        enddo

        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,t_accum_ll)
        call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,r_accum_ll)

     endif
  endif

! HDF5 write

  if (nl%ioutput_latlon == 1) then

     if (myrank == 0) write(io6,'(/,a)') "Writing lat/lon fields to disk..."

     call makefnam(hnamel, hfilepref, current_time, 'LL', '$', 'h5')
     call shdf5_open(hnamel,'W',iclobber) 

! Write any global attributes to file

     write(ofrq,'(I0,A1)') nint(nl%frqlatlon), 's'

     call shdf5_write_global_attribute("Conventions",    cvalue = "CF-1.0")
     call shdf5_write_global_attribute("model",          cvalue = "OLAM")
     call shdf5_write_global_attribute("time_frequency", cvalue = trim(ofrq))
     call shdf5_write_global_attribute("grid",           cvalue = "hexagonal")

! Write coordinate variables to disk (lat, lon, height)

     ndims    = 1
     idims(2) = 1
     idims(3) = 1

     idims(1) = nlon

     CALL shdf5_orec(ndims, idims, 'lon', rvar1=alon, isdim=.true., &
                     long_name = "longitude",                       &
                     standard_name = "longitude",                   &
                     units = "degrees_east"                         )

     idims(1) = nlat

     CALL shdf5_orec(ndims, idims, 'lat', rvar1=alat, isdim=.true., &
                     long_name = "latitude",                        &
                     standard_name = "latitude",                    &
                     units = "degrees_north"                        )

     idims(1) = mza-1

     CALL shdf5_orec(ndims, idims, 'z', rvar1=zlev, isdim=.true., &
                    long_name = "height above mean sea level",    &
                     standard_name = "altitude",                  &
                     units = "m",                                 &
                     positive = "up"                              )

  ! If we ever output data on pressure levels:

     if (dopress) then

        idims(1) = npress

        CALL shdf5_orec(ndims, idims, 'pres', rvar1=plev, isdim=.true., &
                        long_name = "pressure",                         &
                        standard_name = "air_pressure",                 &
                        units = "hPa",                                  &
                        positive = "down"                               )
     endif

! Now write lat/lon interpolated variables to disk.  
! THESE WRITES NEED THE ROUTINE SHDF5_OREC_LL TO WORK IN PARALLEL!!

     ndims    = 2
     idims(1) = nlon
     idims(2) = nlat
     dimnames(1) = 'lon'
     dimnames(2) = 'lat'

! Surface Quantities

     if (dosfc) then

        CALL shdf5_orec_ll(ndims, idims, 'TOPO_LL', rvar1=topo_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                     &
                           long_name = "topography height",                         &
                           standard_name = "surface_altitude",                      &
                           units = "m",                                             &
                           rmissing = rmissing                                      )

        CALL shdf5_orec_ll(ndims, idims, 'U_LPW_LL', rvar1=u_lpw_ll, gpoints=lls_loc,       &
                           dimnames = dimnames,                                             &
                           long_name = "eastward wind at lowest model layer above surface", &
                           standard_name = "eastward_wind",                                 &
                           units = "m s-1",                                                 &
                           rmissing = rmissing                                              )

        CALL shdf5_orec_ll(ndims, idims, 'V_LPW_LL', rvar1=v_lpw_ll, gpoints=lls_loc,        &
                           dimnames = dimnames,                                              &
                           long_name = "northward wind at lowest model layer above surface", &
                           standard_name = "northward_wind",                                 &
                           units = "m s-1",                                                  &
                           rmissing = rmissing                                               )

        CALL shdf5_orec_ll(ndims, idims, 'T_LPW_LL', rvar1=t_lpw_ll, gpoints=lls_loc,         &
                           dimnames = dimnames,                                               &
                           long_name = "air temperature at lowest model layer above surface", &
                           standard_name = "surface_air_temperature",                         &
                           units = "K",                                                       &
                           rmissing = rmissing                                                )

        CALL shdf5_orec_ll(ndims, idims, 'R_LPW_LL', rvar1=r_lpw_ll, gpoints=lls_loc,                      &
                           dimnames = dimnames,                                                            &
                           long_name = "water vapor specific density at lowest model layer above surface", &
                           standard_name = "vapor_specific_density",                                       &
                           units = "kg kg-1",                                                              &
                           rmissing = rmissing                                                             )

        CALL shdf5_orec_ll(ndims, idims, 'SLP_LL', rvar1=slp_ll, gpoints=lls_loc,    &
                           dimnames = dimnames,                                      &
                           long_name = "surface pressure reduced to mean sea level", &
                           standard_name = "surface_air_pressure_at_sea_level",      &
                           units = "Pa",                                             &
                           rmissing = rmissing                                       )

        CALL shdf5_orec_ll(ndims, idims, 'PVAP_LPW_LL', rvar1=pvap_lpw_ll, gpoints=lls_loc,        &
                           dimnames = dimnames,                                                    &
                           long_name = "water vapor pressure at lowest model layer above surface", &
                           standard_name = "water_vapor_partial_pressure_in_air",                  &
                           units = "Pa",                                                           &
                           rmissing = rmissing                                                     )
     endif

  ! Accumulated precipitation/water vapor at surface

     if (doaccum) then

        CALL shdf5_orec_ll(ndims, idims, 'ACCPMIC_LL', rvar1=accpmic_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated resolved precipitation",              &
                           standard_name = "large_scale_precipitation_amount",            &
                           units = "kg m-2",                                              &
                           rmissing = rmissing                                            )

        CALL shdf5_orec_ll(ndims, idims, 'ACCPCON_LL', rvar1=accpcon_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated convective precipitation",            &
                           standard_name = "convective_precipitation_amount",             &
                           units = "kg m-2",                                              &
                           rmissing = rmissing                                            )

        CALL shdf5_orec_ll(ndims, idims, 'VAPFLUX_ACCUM_LL', rvar1=vapflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                       &
                           long_name = "Accumulated water vapor at surface",                          &
                           standard_name = "integral_of_surface_upward_water_vapor_flux_wrt_time",    &
                           units = "kg m-2",                                                          &
                           rmissing = rmissing                                                        )

! Accumulated sensible/latent heat fluxes at surface

        CALL shdf5_orec_ll(ndims, idims, 'SENSFLUX_ACCUM_LL', rvar1=sensflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                         &
                           long_name = "Accumulated upward sensible heat flux at surface",              &
                           standard_name = "integral_of_surface_upward_sensible_heat_flux_wrt_time",    &
                           units = "J m-2",                                                             &
                           rmissing = rmissing                                                          )

        CALL shdf5_orec_ll(ndims, idims, 'LATFLUX_ACCUM_LL', rvar1=latflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                        &
                           long_name = "Accumulated upward latent heat flux at surface",               &
                           standard_name = "integral_of_surface_upward_latent_heat_flux_wrt_time",     &
                           units = "J m-2",                                                            &
                           rmissing = rmissing                                                         )

! Accumulated longwave and shortwave radiative fluxes at surface

        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_ACCUM_LL', rvar1=rshort_accum_ll, gpoints=lls_loc,   &
                           dimnames = dimnames,                                                       &
                           long_name = "Accumulated downwelling shortwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_shortwave_flux_wrt_time", &
                           units = "J m-2",                                                           &
                           rmissing = rmissing                                                        )

        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_ACCUM_LL', rvar1=rshortup_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                         &
                           long_name = "Accumulated upwelling shortwave flux at surface",               &
                           standard_name = "integral_of_surface_upwelling_shortwave_flux_wrt_time",     &
                           units = "J m-2",                                                             &
                           rmissing = rmissing                                                          )
                    
        CALL shdf5_orec_ll(ndims, idims, 'RLONG_ACCUM_LL', rvar1=rlong_accum_ll, gpoints=lls_loc,    &
                           dimnames = dimnames,                                                      &
                           long_name = "Accumulated downwelling longwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_longwave_flux_wrt_time", &
                           units = "J m-2",                                                          &
                           rmissing = rmissing                                                       )

        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_ACCUM_LL', rvar1=rlongup_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                       &
                           long_name = "Accumulated upwelling longwave flux at surface",              &
                           standard_name = "integral_of_surface_upwelling_longwave_flux_wrt_time",    &
                           units = "J m-2",                                                           &
                           rmissing = rmissing                                                        )

! Accumulated longwave and shortwave radiative fluxes at top-of-atmosphere
! Note: at TOA "incoming" and "outgoing" are used in place of "downwelling" and "upwelling"

        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_TOP_ACCUM_LL', rvar1=rshort_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                             &
                           long_name = "Accumulated incoming shortwave flux at TOA",                        &
                           standard_name = "integral_of_toa_incoming_shortwave_flux_wrt_time",              &
                           units = "J m-2",                                                                 &
                           rmissing = rmissing                                                              )

        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_TOP_ACCUM_LL', rvar1=rshortup_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                 &
                           long_name = "Accumulated outgoing shortwave flux at TOA",                            &
                           standard_name = "integral_of_toa_outgoing_shortwave_flux_wrt_time",                  &
                           units = "J m-2",                                                                     &
                           rmissing = rmissing                                                                  )

        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_TOP_ACCUM_LL' , rvar1=rlongup_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "Accumulated outgoing longwave flux at TOA",                            &
                           standard_name = "integral_of_toa_outgoing_longwave_flux_wrt_time",                  &
                           units = "J m-2",                                                                    &
                           rmissing = rmissing                                                                 )

! Accumulated longwave and shortwave clear-air radiative fluxes at surface

        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_CLR_ACCUM_LL', rvar1=rshort_clr_accum_ll, gpoints=lls_loc,     &
                           dimnames = dimnames,                                                                 &
                           long_name = "Accumulated downwelling clear-air shortwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_clear_air_shortwave_flux_wrt_time", &
                           units = "J m-2",                                                                     &
                           rmissing = rmissing                                                                  )

        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_CLR_ACCUM_LL', rvar1=rshortup_clr_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                 &
                           long_name = "Accumulated upwelling clear-air shortwave flux at surface",             &
                           standard_name = "integral_of_surface_upwelling_clear_air_shortwave_flux_wrt_time",   &
                           units = "J m-2",                                                                     &
                           rmissing = rmissing                                                                  )
                    
        CALL shdf5_orec_ll(ndims, idims, 'RLONG_CLR_ACCUM_LL', rvar1=rlong_clr_accum_ll, gpoints=lls_loc,      &
                           dimnames = dimnames,                                                                &
                           long_name = "Accumulated downwelling clear-air longwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_clear_air_longwave_flux_wrt_time", &
                           units = "J m-2",                                                                    &
                           rmissing = rmissing                                                                 )

        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_CLR_ACCUM_LL', rvar1=rlongup_clr_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                               &
                           long_name = "Accumulated upwelling clear-air longwave flux at surface",            &
                           standard_name = "integral_of_surface_upwelling_clear_air_longwave_flux_wrt_time",  &
                           units = "J m-2",                                                                   &
                           rmissing = rmissing                                                                )

! Accumulated longwave and shortwave clear-air radiative fluxes at top-of-atmosphere
! Note: at TOA "incoming" and "outgoing" are used in place of "downwelling" and "upwelling"

        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_TOP_CLR_ACCUM_LL', rvar1=rshort_top_clr_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                     &
                           long_name = "Accumulated incoming clear-air shortwave flux at TOA",                      &
                           standard_name = "integral_of_toa_incoming_clear_air_shortwave_flux_wrt_time",            &
                           units = "J m-2",                                                                         &
                           rmissing = rmissing                                                                      )

        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_TOP_CLR_ACCUM_LL', rvar1=rshortup_top_clr_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                         &
                           long_name = "Accumulated outgoing clear-air shortwave flux at TOA",                          &
                           standard_name = "integral_of_toa_outgoing_clear_air_shortwave_flux_wrt_time",                &
                           units = "J m-2",                                                                             &
                           rmissing = rmissing                                                                          )

        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_TOP_CLR_ACCUM_LL' , rvar1=rlongup_top_clr_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                        &
                           long_name = "Accumulated outgoing clear-air longwave flux at TOA",                          &
                           standard_name = "integral_of_toa_outgoing_clear_air_longwave_flux_wrt_time",                &
                           units = "J m-2",                                                                            &
                           rmissing = rmissing                                                                         )

! 'ALS' quantities are area-weighted averages of land & sea cell quantities over a single atmosphere column

        CALL shdf5_orec_ll(ndims, idims, 'ALS_VELS_ACCUM_LL' , rvar1=als_vels_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                          &
                           long_name = "ALS average of accumulated atmosphere surface wind speed",       &
                           standard_name = "als_average_of_integral_of_atm_wind_speed_wrt_time",         &
                           units = "m",                                                                  &
                           rmissing = rmissing                                                           )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_AIRTEMPK_ACCUM_LL' , rvar1=als_airtempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                  &
                           long_name = "ALS average of accumulated atmosphere temperature",                      &
                           standard_name = "als_average_of_integral_of_atm_temperature_wrt_time",                &
                           units = "K s",                                                                        &
                           rmissing = rmissing                                                                   )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_AIRSHV_ACCUM_LL' , rvar1=als_airshv_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                              &
                           long_name = "ALS average of accumulated atmosphere vapor specific density",       &
                           standard_name = "als_average_of_integral_of_atm_vapor_specific_density_wrt_time", &
                           units = "kg s kg-1",                                                              &
                           rmissing = rmissing                                                               )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_CANTEMPK_ACCUM_LL' , rvar1=als_cantempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                  &
                           long_name = "ALS average of accumulated canopy air temperature",                      &
                           standard_name = "als_average_of_integral_of_canopy_air_temperature_wrt_time",         &
                           units = "K s",                                                                        &
                           rmissing = rmissing                                                                   )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_CANSHV_ACCUM_LL' , rvar1=als_canshv_accum_ll, gpoints=lls_loc,         &
                           dimnames = dimnames,                                                                     &
                           long_name = "ALS average of accumulated canopy air vapor specific density",              &
                           standard_name = "als_average_of_integral_of_canopy_air_vapor_specific_density_wrt_time", &
                           units = "kg s kg-1",                                                                     &
                           rmissing = rmissing                                                                      )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_SKINTEMPK_ACCUM_LL' , rvar1=als_skintempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                    &
                           long_name = "ALS average of accumulated skin temperature",                              &
                           standard_name = "als_average_of_integral_of_skin_temperature_wrt_time",                 &
                           units = "J m-2",                                                                        &
                           rmissing = rmissing                                                                     )
   
        CALL shdf5_orec_ll(ndims, idims, 'ALS_SENSFLUX_ACCUM_LL' , rvar1=als_sensflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                  &
                           long_name = "ALS average of accumulated sensible heat flux",                          &
                           standard_name = "als_average_of_integral_of_sensible_heat_flux_wrt_time",             &
                           units = "J m-2",                                                                      &
                           rmissing = rmissing                                                                   )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_LATFLUX_ACCUM_LL' , rvar1=als_latflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "ALS average of accumulated latent heat flux",                          &
                           standard_name = "als_average_of_integral_of_latent_heat_flux_wrt_time",             &
                           units = "J m-2",                                                                    &
                           rmissing = rmissing                                                                 )

        CALL shdf5_orec_ll(ndims, idims, 'ALS_VAPFLUX_ACCUM_LL' , rvar1=als_vapflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "ALS average of accumulated vapor flux",                                &
                           standard_name = "als_average_of_integral_of_vapor_flux_wrt_time",                   &
                           units = "kg m-2",                                                                   &
                           rmissing = rmissing                                                                 )

! 'AL' quantities are area-weighted averages of land cell quantities over a single atmosphere column

        CALL shdf5_orec_ll(ndims, idims, 'AL_WXFER1_ACCUM_LL' , rvar1=al_wxfer1_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of accumulated soil bottom water flux",                 &
                           standard_name = "al_average_of_integral_of_soil_bottom_water_flux_wrt_time",    &
                           units = "m",                                                                    &
                           rmissing = rmissing                                                             )

        CALL shdf5_orec_ll(ndims, idims, 'AL_SFCWATER_TOT_LL' , rvar1=al_sfcwater_tot_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of total surface water",                                &
                           standard_name = "al_average_of_total_surface_water",                            &
                           units = "kg m-2",                                                               &
                           rmissing = rmissing                                                             )

        CALL shdf5_orec_ll(ndims, idims, 'AL_SOIL_WATER_TOT_LL' , rvar1=al_soil_water_tot_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "AL average of total soil water",                                       &
                           standard_name = "al_average_of_total_soil_water",                                   &
                           units = "m",                                                                        &
                           rmissing = rmissing                                                                 )

! 3D accumulated atmospheric fields

        ndims = 3
        idims(1) = nlon
        idims(2) = nlat
        idims(3) = mza-1
        dimnames(1) = 'lon'
        dimnames(2) = 'lat'
        dimnames(3) = 'z'

        if (do3d) then

           CALL shdf5_orec_ll(ndims, idims, 'U_ACCUM_LL', rvar2=u_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated eastward wind",                       &
                              standard_name = "integral_of_eastward_wind_wrt_time",          &
                              units = "m",                                                   &
                              rmissing=rmissing                                              )

           CALL shdf5_orec_ll(ndims, idims, 'V_ACCUM_LL', rvar2=v_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated northward wind",                      &
                              standard_name = "integral_of_northward_wind_wrt_time",         &
                              units = "m",                                                   &
                              rmissing=rmissing                                              )

           CALL shdf5_orec_ll(ndims, idims, 'W_ACCUM_LL', rvar2=w_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated upward wind",                         &
                              standard_name = "integral_of_upward_wind_wrt_time",            &
                              units = "m",                                                   &
                              rmissing=rmissing                                              )

           CALL shdf5_orec_ll(ndims, idims, 'T_ACCUM_LL', rvar2=t_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated air temperature",                     &
                              standard_name = "integral_of_air_temperature_wrt_time",        &
                              units = "K s",                                                 &
                              rmissing=rmissing                                              )

           CALL shdf5_orec_ll(ndims, idims, 'R_ACCUM_LL', rvar2=r_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated water vapor specific density ratio",  &
                              standard_name = "integral_of_vapor_specific_density_wrt_time", &
                              units = "kg s kg-1",                                           &
                              rmissing=rmissing                                              )

           CALL shdf5_orec_ll(ndims, idims, 'P_ACCUM_LL', rvar2=p_accum_ll, gpoints=lls_loc, &
                              dimnames = dimnames,                                           &
                              long_name = "Accumulated air pressure",                        &
                              standard_name = "integral_of_air_pressure_wrt_time",           &
                              units = "Pa s",                                                &
                              rmissing=rmissing                                              )

        endif

     endif


! 3D atmospheric fields

     ndims = 3
     idims(1) = nlon
     idims(2) = nlat
     idims(3) = mza-1
     dimnames(1) = 'lon'
     dimnames(2) = 'lat'
     dimnames(3) = 'z'

     if (do3d) then

        CALL shdf5_orec_ll(ndims, idims, 'U_LL', rvar2=u_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "eastward wind",                       &
                           standard_name = "eastward_wind",                   &
                           units = "m s-1",                                   &
                           rmissing=rmissing                                  )

        CALL shdf5_orec_ll(ndims, idims, 'V_LL', rvar2=v_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "northward wind",                      &
                           standard_name = "northward_wind",                  &
                           units = "m s-1",                                   &
                           rmissing=rmissing                                  )

        CALL shdf5_orec_ll(ndims, idims, 'W_LL', rvar2=w_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "upward wind",                         &
                           standard_name = "upward_wind",                     &
                           units = "m s-1",                                   &
                           rmissing=rmissing                                  )

        CALL shdf5_orec_ll(ndims, idims, 'T_LL', rvar2=t_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "air temperature",                     &
                           standard_name = "air_temperature",                 &
                           units = "K",                                       &
                           rmissing=rmissing                                  )

        CALL shdf5_orec_ll(ndims, idims, 'R_LL', rvar2=r_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "water vapor specific density ratio",            &
                           standard_name = "vapor_specific_density",           &
                           units = "kg kg-1",                                 &
                           rmissing=rmissing                                  )

        CALL shdf5_orec_ll(ndims, idims, 'P_LL', rvar2=p_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "air pressure",                        &
                           standard_name = "air_pressure",                    &
                           units = "Pa",                                      &
                           rmissing=rmissing                                  )

     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TODO: output variables on pressure levels???
!
! ndims = 3
! idims(1) = nlon
! idims(2) = nlat
! idims(3) = npres
! dimnames(1) = 'lon'
! dimnames(2) = 'lat'
! dimnames(3) = 'pres'
!
! if (dopress) then
!     ....
! endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SEIGEL 2013 - close interpolation write
     call shdf5_close()

  endif

! Return if we do not want to plot the interpolated fields

  if (nl%latlonplot == 0) return

! Reopen the current graphics output workstation if it is closed

  if (myrank == 0) call o_reopnwk()
  if (myrank == 0) write(io6,'(/,a)') "Plotting interpolated lat/lon fields..."

  if (dosfc) then
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,topo_ll,402, &
                      "Topography Height", "(m)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,u_lpw_ll,113, &
                      "Zonal wind at lpw", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,v_lpw_ll,113, &
                      "Meridional wind at lpw", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,t_lpw_ll,4, &
                      "Air temperature at lpw", "(K)", rmissing)

     ! Preserve missing values
     where( abs(r_lpw_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) = r_lpw_ll(:) * 1000.
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,5, &
                      "Water vapor specific density at lpw", "(g/kg)", rmissing)

     ! Preserve missing values
     where( abs(slp_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) =  slp_ll(:) * .01 
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,17, &
                      "Sea level pressure at suface", "(mb)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pvap_lpw_ll,35, &
                      "Water vapor pressure at lpw", "(mb)", rmissing)
  endif

  if (doaccum) then
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpmic_ll,204, &
                      "Accum microphysics precip", "(kg/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpcon_ll,204, &
                      "Accum convective precip", "(kg/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,vapflux_accum_ll,304, &
                      "Accum surface vapor", "(kg/m^2)", rmissing)

     scr2_ll(:) = sensflux_accum_ll(:) * 1.e-6 
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,304, &
                      "Accum surface sensible heat", "(MJ/m^2)", rmissing)

     scr2_ll(:) = latflux_accum_ll(:) * 1.e-6 
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,304, &
                      "Accum surface latent heat", "(MJ/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_accum_ll,204, &
                      "Accum surface downward s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_accum_ll,204, &
                      "Accum surface upward s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlong_accum_ll,204, &
                      "Accum surface downward l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_accum_ll,204, &
                      "Accum surface upward l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_top_accum_ll,204, &
                      "Accum TOA downward s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_top_accum_ll,204, &
                      "Accum TOA upward s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_top_accum_ll,204, &
                      "Accum TOA upward l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_clr_accum_ll,204, &
                      "Accum surface downward clear-air s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_clr_accum_ll,204, &
                      "Accum surface upward clear-air s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlong_clr_accum_ll,204, &
                      "Accum surface downward clear-air l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_clr_accum_ll,204, &
                      "Accum surface upward clear-air l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_top_clr_accum_ll,204, &
                      "Accum TOA downward clear-air s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_top_clr_accum_ll,204, &
                      "Accum TOA upward clear-air s/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_top_clr_accum_ll,204, &
                      "Accum TOA upward clear-air l/w radiation", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_vels_accum_ll,204, &
                      "als accum wind speed", "(m)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_airtempk_accum_ll,204, &
                      "als accum atmosphere temperature", "(K s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_airshv_accum_ll,204, &
                      "als accum atmosphere vapor specific density", "(kg s / kg)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_cantempk_accum_ll,204, &
                      "als accum canopy air temperature", "(K s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_canshv_accum_ll,204, &
                      "als accum canopy air vapor specific density", "(kg s / kg)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_skintempk_accum_ll,204, &
                      "als accum skin temperature", "(K s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_sensflux_accum_ll,304, &
                      "als accum sensible heat flux", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_latflux_accum_ll,304, &
                      "als accum latent heat flux", "(J/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,als_vapflux_accum_ll,304, &
                      "als accum vapor flux", "(kg/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_wxfer1_accum_ll,304, &
                      "als accum soil bottom water flux", "(m)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_sfcwater_tot_ll,200, &
                      "als total surface water", "(kg/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_soil_water_tot_ll,3, &
                      "als total soil water", "(m)", rmissing)

!------------------------------------------------------------
! Contour plot or tile plot 2D slices of 3D lat-lon fields
!------------------------------------------------------------

     if (do3d) then

        ! klev = vertical level to plot in 3D fields

        klev = min(23,mza)

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,u_accum_ll,304, &
                         "Accumulated zonal wind", "(m)", rmissing)

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,v_accum_ll,304, &
                         "Accumulated meridional wind", "(m)", rmissing)

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,w_accum_ll,304, &
                         "Accumulated vertical wind", "(m)", rmissing)

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,t_accum_ll,204, &
                         "Accumulated temperature", "(K s)", rmissing)

        ! Preserve missing values
        where( abs(r_ll(:,:) - rmissing) > 1.e-7 )
           scr3_ll(:,:) = r_accum_ll(:,:) * 1000.
        elsewhere
           scr3_ll(:,:) = rmissing
        endwhere

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,204, &
             "Accumulated water vapor specific density", "(g s/kg)", rmissing)

        ! Preserve missing values
        where( abs(p_ll(:,:) - rmissing) > 1.e-7 )
           scr3_ll(:,:) = p_accum_ll(:,:) * 0.01
        elsewhere
           scr3_ll(:,:) = rmissing
        endwhere

        call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,204, &
             "Accumulated air pressure", "(hPa s)", rmissing)

     endif

  endif

!------------------------------------------------------------
! Contour plot or tile plot 2D slices of 3D lat-lon fields
!------------------------------------------------------------

  if (do3d) then

     ! klev = vertical level to plot in 3D fields

     klev = min(23,mza)

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,u_ll,113, &
                      "Zonal wind", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,v_ll,113, &
                      "Meridional wind", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,w_ll,125, &
                      "Vertical wind", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,t_ll,4, &
                      "Temperature", "(K)", rmissing)

     ! Preserve missing values
     where( abs(r_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = r_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,5, &
          "Water vapor specific density", "(g/kg)", rmissing)

     ! Preserve missing values
     where( abs(p_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = p_ll(:,:) * 0.01
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,415, &
          "Air Pressure", "(mb)", rmissing)

  endif

  if (dodifs) then

!---------------------------------------------------------------------------
! Contour plot 1D arrays that are longitudinal averages of 2D lat-lon fields
! This will only work for serial runs that have the entire lat/lon memory
!---------------------------------------------------------------------------

     if (iparallel == 1 .or. npts /= nlon*nlat) then
        ! For a parallel run, close graphics and return to program
        if (myrank == 0) call o_clswk()
        return
     endif

     timefac = 1.0
     if (abs(time8 - time8_prev) > .99) then
        timefac = 1.0 / real(time8 - time8_prev)
     endif

           pcpmic_dif2_ll(:) = (           accpmic_ll(:) -       pcpmic_dif2_ll(:)) * timefac * 86400.
           pcpcon_dif2_ll(:) = (           accpcon_ll(:) -       pcpcon_dif2_ll(:)) * timefac * 86400.
          pcpboth_dif2_ll(:) =         pcpmic_dif2_ll(:) +       pcpcon_dif2_ll(:)
           rshort_dif2_ll(:) = (      rshort_accum_ll(:) -       rshort_dif2_ll(:)) * timefac
         rshortup_dif2_ll(:) = (    rshortup_accum_ll(:) -     rshortup_dif2_ll(:)) * timefac
            rlong_dif2_ll(:) = (       rlong_accum_ll(:) -        rlong_dif2_ll(:)) * timefac
          rlongup_dif2_ll(:) = (     rlongup_accum_ll(:) -      rlongup_dif2_ll(:)) * timefac
       rshort_top_dif2_ll(:) = (  rshort_top_accum_ll(:) -   rshort_top_dif2_ll(:)) * timefac
     rshortup_top_dif2_ll(:) = (rshortup_top_accum_ll(:) - rshortup_top_dif2_ll(:)) * timefac
      rlongup_top_dif2_ll(:) = ( rlongup_top_accum_ll(:) -  rlongup_top_dif2_ll(:)) * timefac
         sensflux_dif2_ll(:) = (    sensflux_accum_ll(:) -     sensflux_dif2_ll(:)) * timefac
          latflux_dif2_ll(:) = (     latflux_accum_ll(:) -      latflux_dif2_ll(:)) * timefac

     aspect = .7
     scalelab = .012
     alatinc = 10.

     call plotback()

     allocate(value(nlat))
     allocate(arr_ll(nlon,nlat))

!---------------------------------------------------------------

     ymin = -1.
     ymax = 10.
     yinc = 1.

     arr_ll = reshape( pcpmic_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,8,0,   &
                    nlat, alat, value,             &
                    'latitude','PCPDIF2 (mm/day)', &
                    alat(1),alat(nlat),alatinc,3,  &
                    ymin, ymax, yinc, 5            )

     arr_ll = reshape(pcpcon_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,1,0,  &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape(pcpboth_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,11,0, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

!---------------------------------------------------------------

     ymin = -10.
     ymax = 500.
     yinc = 10.

     arr_ll = reshape(rshort_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,8,0,  &
                    nlat, alat, value,            &
                    'latitude','SFC RAD (W/m^2)', &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape(rshortup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape(rlong_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,1,0,  &
                    nlat, alat, value,            &
                    'latitude','',                &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape(rlongup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

! Plot reference lines

     value(1:nlat) = 100.
     call oplot_xy2('2','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 200.
     call oplot_xy2('2','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 300.
     call oplot_xy2('2','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 400.
     call oplot_xy2('2','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

!---------------------------------------------------------------

     arr_ll = reshape(rshort_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,8,0,   &
                    nlat, alat, value,             &
                    'latitude','TOA RAD (W/m^2) ', &
                    alat(1),alat(nlat),alatinc,3,  &
                    ymin, ymax, yinc, 5            )

     arr_ll = reshape(rshortup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape(rlongup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

! Plot reference lines

     value(1:nlat) = 100.
     call oplot_xy2('3','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 200.
     call oplot_xy2('3','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 300.
     call oplot_xy2('3','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     value(1:nlat) = 400.
     call oplot_xy2('3','N',aspect,scalelab,10,0, &
                    nlat, alat, value,            &
                    ' ',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

!---------------------------------------------------------------

     ymin = -50.
     ymax = 200.
     yinc = 10.

     arr_ll = reshape(sensflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N',aspect,scalelab,17,0,      &
                    nlat, alat, value,                 &
                    'latitude','SFC S/L FLUX (W/m^2)', &
                    alat(1),alat(nlat),alatinc,3,      &
                    ymin, ymax, yinc, 5                )

     arr_ll = reshape(latflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N',aspect,scalelab,9,0,  &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

!---------------------------------------------------------------

     call o_frame()

     call o_clswk()

           pcpmic_dif2_ll(:) =            accpmic_ll(:)
           pcpcon_dif2_ll(:) =            accpcon_ll(:)
           rshort_dif2_ll(:) =       rshort_accum_ll(:)
         rshortup_dif2_ll(:) =     rshortup_accum_ll(:)
            rlong_dif2_ll(:) =        rlong_accum_ll(:)
          rlongup_dif2_ll(:) =      rlongup_accum_ll(:)
       rshort_top_dif2_ll(:) =   rshort_top_accum_ll(:)
     rshortup_top_dif2_ll(:) = rshortup_top_accum_ll(:)
      rlongup_top_dif2_ll(:) =  rlongup_top_accum_ll(:)
         sensflux_dif2_ll(:) =     sensflux_accum_ll(:)
          latflux_dif2_ll(:) =      latflux_accum_ll(:)

                  time8_prev = time8

  endif

end subroutine fields2_ll

!==========================================================================

subroutine tileplot_ll(nlon,nlat,nlev,ilev,alon,alat,npts,lls_loc,fld,itab, &
                       fldname, units, rmissing)

  use plotcolors, only: clrtab
  use oplot_coms, only: op
  use misc_coms,  only: io6, iparallel
  use mem_para,   only: myrank, mgroupsize, nbytes_int, nbytes_real

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer, intent(in) :: nlon,nlat,nlev,ilev,npts,itab
  integer, intent(in) :: lls_loc(npts)
  real,    intent(in) :: alon(nlon),alat(nlat)
  real,    intent(in) :: fld(npts,nlev)
  real,    intent(in) :: rmissing

  character(*), intent(in) :: fldname
  character(*), intent(in) :: units

  real :: glon(nlon)
  real :: htpn(4),vtpn(4),fldval
  integer :: ilat, ilon, ival, icolor

  integer, allocatable :: buffer(:)
  integer :: nu, ier, buffsize, ipos, base, n, nn, j
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

  real,         parameter :: aspect =  0.
  character(1), parameter :: proj   = 'L'
  character(1), parameter :: cbar   = 'c'
  character(1), parameter :: panel  = 'n'

  real, save :: scalelab = .014

  real    :: dum1(1),dum2(1), bsize

  ! Set up communication buffers for parallel run

#ifdef OLAM_MPI
  if (iparallel == 1) then
     base = 8 * nbytes_real + nbytes_int

     if (myrank > 0) then
        buffsize = npts * base
        allocate( buffer(buffsize) )
     endif

     ipos = 0
     nu   = 0
  endif
#endif

  ! Set plot bounds

  glon(:) = alon(:)

  if (alon(nlon) < alon(1)) then
     do ilon = 1,nlon
        if (glon(ilon) < alon(1)) glon(ilon) = glon(ilon) + 360.
     enddo
  endif

  op%xmin = 2. * glon(1) - glon(2)
  op%xmax = 2. * glon(nlon) - glon(nlon-1)
  op%ymin = 2. * alat(1) - alat(2)
  op%ymax = 2. * alat(nlat) - alat(nlat-1)

  ! Scale plot window and draw title and frame

  if (myrank == 0) call plotback()

  call oplot_panel(panel, cbar, aspect, proj)

  if (myrank == 0) then

     call o_gsplci(10)
     call o_gstxci(10)
     call o_sflush()

     bsize = .016 * (op%hp2 - op%hp1)
     call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)
     call o_plchhq(0.5, op%fnamey, trim(trim(fldname)//' '//trim(units)), bsize, 0., 0.)

     call o_set(op%h1,op%h2,op%v1,op%v2,op%xmin,op%xmax,op%ymin,op%ymax,1)

     call oplot_xy2(panel, cbar, 0., scalelab, 10, 0,     &
                    1, dum1,dum2,                         &
                    'LONGITUDE (deg)', 'LATITUDE (deg)',  &
                    op%xmin, op%xmax, 10., 3,             &
                    op%ymin, op%ymax, 10., 3              )
  endif

  ! Tile plot each lat/lon point
 
  do nn = 1, npts

     ! 1d global lat/lon index

     n = lls_loc(nn)            

     ! convert 1d lat/lon index to 2d

     ilat = (n-1)/nlon + 1
     ilon = mod(n-1,nlon) + 1

     ! set field value and skip missing data

     fldval = fld(nn,ilev)

     if (abs(fldval - rmissing) < 1.e-7) cycle
     
     ! get lat/lon cell boundaries

     if (ilat == 1) then
        vtpn(1) = 1.5 * alat(1) - .5 * alat(2)
     else
        vtpn(1) = .5 * (alat(ilat) + alat(ilat-1))
     endif

     if (ilat == nlat) then
        vtpn(3) = 1.5 * alat(nlat) - .5 * alat(nlat-1)
     else
        vtpn(3) = .5 * (alat(ilat) + alat(ilat+1))
     endif

     vtpn(2) = vtpn(1)
     vtpn(4) = vtpn(3)

     if (ilon == 1) then
        htpn(1) = 1.5 * glon(1) - .5 * glon(2)
     else
        htpn(1) = .5 * (glon(ilon) + glon(ilon-1))
     endif

     if (ilon == nlon) then
        htpn(2) = 1.5 * glon(nlon) - .5 * glon(nlon-1)
     else
        htpn(2) = .5 * (glon(ilon) + glon(ilon+1))
     endif

     htpn(3) = htpn(2)
     htpn(4) = htpn(1)

     ! Set field value and extract contour color from color table

     fldval = fld(nn,ilev)

     ival = 1
     do while (fldval > clrtab(itab)%vals(ival) .and.  &
                 ival < clrtab(itab)%nvals             )
        ival = ival + 1
     enddo
     icolor = clrtab(itab)%ipal(ival)

     ! Plot the current pixel, or send to node 0 to plot

     if (myrank == 0) then

        call fillpolyg(4,htpn,vtpn,icolor)

     else

#ifdef OLAM_MPI
        nu = nu + 1
        call MPI_Pack(htpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(vtpn,   4, MPI_REAL,    buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
        call MPI_Pack(icolor, 1, MPI_INTEGER, buffer, buffsize, ipos, MPI_COMM_WORLD, ier)
#endif

     endif

  enddo  ! local lat/lon loop

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

     if (myrank > 0 .and. nu > 0) then
        call MPI_Send(buffer, ipos, MPI_PACKED, 0, itag, MPI_COMM_WORLD, ier)
     endif

     if (myrank == 0) then

        buffsize = maxval(nus(2:mgroupsize)) * base
        allocate( buffer( buffsize ) )

        do n = 2, mgroupsize

           if (nus(n) > 0) then

              call MPI_Recv( buffer, buffsize, MPI_PACKED, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

              ipos = 0

              do j = 1, nus(n)
               
                 call MPI_Unpack(buffer, buffsize, ipos, htpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, vtpn,   4, MPI_REAL,    MPI_COMM_WORLD, ier)
                 call MPI_Unpack(buffer, buffsize, ipos, icolor, 1, MPI_INTEGER, MPI_COMM_WORLD, ier)

                 call fillpolyg(4,htpn,vtpn,icolor)

              enddo

           endif
        enddo
     endif
        
     deallocate(buffer)
  endif
#endif

  if (myrank == 0) then
     call mkmap_ll()
     if (cbar == 'c') then
        ! trick trunc_polyg since our user/fractional grid mapping is different
        op%xmin = 0.0
        op%ymin = 0.0
        op%xmax = 1.0
        op%ymax = 1.0
        call plot_colorbar(itab)
     endif
     call o_frame()
  endif

end subroutine tileplot_ll

!==========================================================================

subroutine contplot_ll(nlon,nlat,nlev,ilev,alon,alat,fld,itab)

  use oplot_coms, only: op
  implicit none

  integer, intent(in) :: nlon,nlat,nlev,ilev,itab
  real, intent(in) :: alon(nlon),alat(nlat)
  real, intent(in) :: fld(nlon,nlat,nlev)

  real :: htpn(4),vtpn(4),fldvals(4)
  real :: glon(nlon)
  integer :: ilat, ilon

  call plotback()

  glon(:) = alon(:)

  if (alon(nlon) < alon(1)) then
     do ilon = 1,nlon
        if (glon(ilon) < alon(1)) glon(ilon) = glon(ilon) + 360.
     enddo
  endif

  op%xmin = 2. * glon(1) - glon(2)
  op%xmax = 2. * glon(nlon) - glon(nlon-1)
  op%ymin = 2. * alat(1) - alat(2)
  op%ymax = 2. * alat(nlat) - alat(nlat-1)

! The following call to o_set uses full plot window, which will distort shape
! unless op%xmax - op%xmin = op%ymax - op%ymin

  call o_set(0.,1.,0.,1.,op%xmin,op%xmax,op%ymin,op%ymax,1)

  do ilat = 1,nlat-1
     vtpn(1) = alat(ilat)
     vtpn(2) = vtpn(1)
     vtpn(3) = alat(ilat+1)
     vtpn(4) = vtpn(3)

     do ilon = 1,nlon-1
        htpn(1) = glon(ilon)
        htpn(2) = glon(ilon+1)
        htpn(3) = htpn(2)
        htpn(4) = htpn(1)

! Set field value and extract contour color from color table

        fldvals(1) = fld(ilon  ,ilat  ,ilev)
        fldvals(2) = fld(ilon+1,ilat  ,ilev)
        fldvals(3) = fld(ilon+1,ilat+1,ilev)
        fldvals(4) = fld(ilon  ,ilat+1,ilev)

        call contpolyg(itab,1,4,htpn,vtpn,fldvals)
        call contpolyg(itab,0,4,htpn,vtpn,fldvals)
     enddo
  enddo

  call o_frame()

end subroutine contplot_ll

!==========================================================================

subroutine mkmap_ll()

  use oplot_coms,  only: op
  use mem_para,    only: myrank

  implicit none
  
  if (myrank /= 0) return

  ! Set plot color (black)

  call o_gsplci(10)
  call o_gsfaci(10)
  call o_gstxci(10)
  call o_sflush()

  call o_mapint()
  call o_mappos(op%h1,op%h2,op%v1,op%v2)

  call o_mapsti('GR',20)    ! For 20-degree spacing of plotted parallels and meridians
  call o_mapsti('DA',65535) ! To plot parallels and meridians with solid lines
! call o_mapstc('OU','PS')  ! To plot bnds of continents, countries, and US  states
! call o_mapstc('OU','NO')  ! To plot NO geographic information
  call o_mapstc('OU','CO')  ! To plot continental outlines
! call o_mapstc('OU','US')  ! To plot US state outlines
! call o_mapstc('OU','PO')  ! To plot continental outlines + international outlines

  call o_maproj('CE', 0., 0., 0.)  ! Force plat to be zero for CE projection
  call o_mapset('LI', op%xmin, op%xmax, op%ymin, op%ymax)

  call o_mapint()      ! Initialize the above parameters
! call mapgrd()
  call o_maplot()
  call o_sflush()

end subroutine mkmap_ll
