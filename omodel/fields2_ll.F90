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

  use mem_ijtabs,  only: itab_w, itab_v, itabg_m
  use mem_basic,   only: vc, wc, rho, press, theta, sh_w, sh_v, &
                         vxe, vye, vze, tair

  use mem_grid,    only: mza, mva, mwa, lpv, lpw, &
                         xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                         topw, glatm, glonm, glatw, glonw, zm, zt, &
                         vnx, vny, vnz, dzt

  use misc_coms,   only: io6, current_time, hfilepref, iclobber, iparallel

  use consts_coms, only: p00, rocp, piu180, erad, eradi, pio180, cp, alvl, rvap, r8

  use hdf5_utils,  only: shdf5_open, shdf5_orec, shdf5_orec_ll, shdf5_close

  use max_dims,    only: pathlen
  use oname_coms,  only: nl
  use mem_para,    only: myrank
  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w

! NOTE: The fields used from the MEM_MICRO, MEM_CUPARM, and MEM_FLUX_ACCUM
! modules are time-integrated fluxes of mass or energy and have units of
! [kg/m^2] or [J/m^2].  Their values written to any history file represent flux
! integrals from the beginning of a simulation to the time of the history file
! write.  Thus, time averages of a flux over any time interval can be computed
! as the difference between the flux integrals in two different history files
! divided by the difference in simulation times when the files were written.

  use mem_micro,  only: accpd, accpr, accpp, accps, accpa, accpg, accph

  use mem_cuparm, only: aconpr

  use mem_flux_accum, only: rshort_accum, rshortup_accum, rlong_accum, rlongup_accum, &
                            rshort_top_accum, rshortup_top_accum, rlongup_top_accum, &
                            sfluxt_accum, sfluxr_accum

! NOTE: The fields used from the MEM_PLOT module are time-integrated fluxes of
! mass or energy and have units of [kg/m^2] or [J/m^2].  They are filled from
! corresponding values from the MEM_MICRO, MEM_CUPARM, and MEM_FLUX_ACCUM
! modules, either as values at a single time or as sums of values read from
! history files at multiple times.  MEM_PLOT field names ending with 'prev0'
! are the most recent (in time) single or summed values, while those ending
! with 'prev1' are single or summed values from an earlier time.  Their
! difference divided by the difference in single or summed times gives the 
! mean flux rate over the differential time period.

  use mem_plot, only: time8_prev0, accpmic_prev0, accpcon_prev0,        &
                      time8_prev1, accpmic_prev1, accpcon_prev1,        &
                      rshort_accum_prev0,     rshortup_accum_prev0,     &
                      rshort_accum_prev1,     rshortup_accum_prev1,     &
                      rlong_accum_prev0,      rlongup_accum_prev0,      &
                      rlong_accum_prev1,      rlongup_accum_prev1,      &
                      rshort_top_accum_prev0, rshortup_top_accum_prev0, &
                      rshort_top_accum_prev1, rshortup_top_accum_prev1, &
                      rlongup_top_accum_prev0,                          &
                      rlongup_top_accum_prev1,                          &
                      sfluxt_accum_prev0,     sfluxr_accum_prev0,       &
                      sfluxt_accum_prev1,     sfluxr_accum_prev1

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

!Regional domain:
!  real   , parameter :: beglon = -92.   ! minimum (westernmost) longitude (deg)
!  real   , parameter :: endlon = -72.   ! maximum (easternmost) longitude (deg)
!  real   , parameter :: beglat =  23.   ! minimum (southernmost) latitude (deg)
!  real   , parameter :: endlat =  35.   ! maximum (northernmost) latitude (deg)
!  integer, parameter :: rf = 30         ! horiz resolution factor (pts per deg)

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
  real, allocatable :: u_sfc_ll             (:) ! surface zonal wind component (m/s)
  real, allocatable :: v_sfc_ll             (:) ! surface meridional wind component (m/s)
  real, allocatable :: t_sfc_ll             (:) ! surface air temperature (K)
  real, allocatable :: r_sfc_ll             (:) ! surface water vapor mixing ratio (kg/kg)
  real, allocatable :: pvap_sfc_ll          (:) ! surface vapor pressure (Pa)
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

  real, allocatable :: pcpdif2mic_ll        (:) ! microphysics precip diff2 (mm/day)
  real, allocatable :: pcpdif2con_ll        (:) ! convective precip diff2 (mm/day)
  real, allocatable :: pcpdif2both_ll       (:) ! total precip diff2 (mm/day)
  real, allocatable :: vapflux_dif2_ll      (:) ! sfc vapor flux diff2 (kg/(m^2 s))
  real, allocatable :: sensflux_dif2_ll     (:) ! sfc sensible heat flux diff2 (W/m^2)
  real, allocatable :: latflux_dif2_ll      (:) ! sfc latent heat flux diff2 (W/m^2)
  real, allocatable :: rshort_dif2_ll       (:) ! sfc downward s/w rad flux diff2 (W/m^2)
  real, allocatable :: rshortup_dif2_ll     (:) ! sfc upward s/w rad flux diff2 (W/m^2)
  real, allocatable :: rlong_dif2_ll        (:) ! sfc downward l/w rad flux diff2 (W/m^2)
  real, allocatable :: rlongup_dif2_ll      (:) ! sfc upward l/w rad flux diff2 (W/m^2)
  real, allocatable :: rshort_top_dif2_ll   (:) ! TOA downward s/w rad flux diff2 (W/m^2)
  real, allocatable :: rshortup_top_dif2_ll (:) ! TOA upward s/w rad flux diff2 (W/m^2)
  real, allocatable :: rlongup_top_dif2_ll  (:) ! TOA upward l/w rad flux diff2 (W/m^2)

!----------
! 3D FIELDS
!----------

  real, allocatable :: u_ll (:,:) ! zonal wind component (m/s)
  real, allocatable :: v_ll (:,:) ! meridional wind component (m/s)
  real, allocatable :: t_ll (:,:) ! air temperature (K)
  real, allocatable :: r_ll (:,:) ! water vapor mixing ratio (kg/kg)
  real, allocatable :: p_ll (:,:) ! air pressure (Pa)

  integer       :: k,iw,ilat,ilon,n,kb
  integer       :: ndims, idims(3)
  character(30) :: dimnames(3)

  real :: raxis,raxisi

  real :: scr1a(mwa), scr1b(mwa)
  real :: scr2a(mza,mwa), scr2b(mza,mwa)
  real, allocatable :: scr2_ll(:), scr3_ll(:,:)

  real :: timefac
  real :: aspect, scalelab, ymin, ymax, yinc, alatinc

! SEIGEL 2013 - Added for ll interp writeout
  character(pathlen) :: hnamel

  real,    save :: dlon, dlat
  integer, save :: npts
  integer, save :: init = 0

  integer, allocatable, save :: ims_loc(:), lls_loc(:), ijs_loc(:)

  integer, allocatable :: ims(:,:)
  integer, allocatable :: ijs(:,:)

  ! These switches control whick lat/lon interpolations are performed

  logical, parameter :: dosfc   = .true.
  logical, parameter :: doaccum = .false.
  logical, parameter :: dodifs  = .false.
  logical, parameter :: do3d    = .false.
  logical, parameter :: dopress = .false. ! TODO!

!------------------------------------------------------
! 3D FIELDS - interpolation to pressure levels (TODO!)
!------------------------------------------------------

  integer, parameter :: npress = 7
  real :: plev(npress) = (/ 1000., 925., 850., 700., 500., 250., 100. /)

!------------------------------------------------------------------------------

  if (nl%ioutput_latlon /= 1) return

!------------------------------------------------------------------------------

  if (.not. (dosfc .or. doaccum .or. dodifs .or. do3d .or. dopress)) return

  if (init == 0) then

     init = 1

     if (myrank == 0) write(io6,'(/,a)') "Initializing lat/lon overlaps..."

     ! Compute latitude and longitude of output grid points (assuming uniform spacing)

     nlon = nint(endlon - beglon) * rf + 1  ! # of lon values
     nlat = nint(endlat - beglat) * rf + 1  ! # of lat values

     dlat = (endlat  - beglat) / real(nlat-1)

     if (endlon >= beglon) then
        dlon = (endlon - beglon) / real(nlon-1)
     else
        dlon = (endlon + 360. - beglon) / real(nlon-1)
     endif
     
     ! for a global domain, don't include the redundant final cyclic point

     if (endlon - 360.0 < beglon + 1.e-6) nlon = nlon - 1

     allocate(alat(nlat))
     allocate(alon(nlon))

     do ilat = 1, nlat
        alat(ilat) = beglat + dlat * real(ilat-1)
     enddo

     do ilon = 1, nlon
        alon(ilon) = beglon + dlon * real(ilon-1)
        if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.
     enddo

     allocate(ims(nlon,nlat))
     allocate(ijs(nlon,nlat))

     ! Allocate and fill zlev with model grid levels

     allocate(zlev(mza-1))

     do k = 2, mza
        zlev(k-1) = zt(k)
     enddo
     
     ! Find the M point that corresponds to each lat/lon point and store
     ! the triangle sector that each lat/lon point is in
     
     call find_closest_m_ll(nlon,nlat,alon,alat,ims,ijs)

     ! Compute number of lat/lon points on current node

     if (iparallel == 0) then

        npts = nlon*nlat

     else

        npts = 0
        do ilat = 1, nlat
           do ilon = 1, nlon
              if (itabg_m(ims(ilon,ilat))%irank == myrank) npts = npts + 1
           enddo
        enddo

     endif
     
     ! Store the lat/lon index and M point of each lat/lon cell on current node

     allocate(ims_loc(npts)) ! local im location of each lat/lon point
     allocate(lls_loc(npts)) ! lat/lon indices on this node
     allocate(ijs_loc(npts)) ! the triangle around the M point that this lat/lon point is in

     n = 0
     do ilat = 1, nlat
        do ilon = 1, nlon
           if (itabg_m(ims(ilon,ilat))%irank == myrank) then
              n = n + 1
              ims_loc(n) = ims(ilon,ilat) 
              ! convert 2d lat/lon index to 1d
              lls_loc(n) = ilon + (ilat-1)*nlon
              ijs_loc(n) = ijs(ilon,ilat)
           endif
        enddo
     enddo

     deallocate(ims)
     deallocate(ijs)

  endif

! Initialize latitude-longitude arrays to zero prior to interpolation.
! For certain applications, zero should be replaced with "missing value".

  if (dosfc) then
     allocate( topo_ll              (npts) ) ; topo_ll               = rmissing
     allocate( u_sfc_ll             (npts) ) ; u_sfc_ll              = rmissing
     allocate( v_sfc_ll             (npts) ) ; v_sfc_ll              = rmissing
     allocate( t_sfc_ll             (npts) ) ; t_sfc_ll              = rmissing
     allocate( r_sfc_ll             (npts) ) ; r_sfc_ll              = rmissing
     allocate( pvap_sfc_ll          (npts) ) ; pvap_sfc_ll           = rmissing
     allocate( slp_ll               (npts) ) ; slp_ll                = rmissing
  endif

  if (doaccum) then
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

  if (dodifs) then
     allocate( pcpdif2mic_ll        (npts) ) ; pcpdif2mic_ll         = rmissing
     allocate( pcpdif2con_ll        (npts) ) ; pcpdif2con_ll         = rmissing
     allocate( pcpdif2both_ll       (npts) ) ; pcpdif2both_ll        = rmissing
     allocate( vapflux_dif2_ll      (npts) ) ; vapflux_dif2_ll       = rmissing
     allocate( sensflux_dif2_ll     (npts) ) ; sensflux_dif2_ll      = rmissing
     allocate( latflux_dif2_ll      (npts) ) ; latflux_dif2_ll       = rmissing
     allocate( rshort_dif2_ll       (npts) ) ; rshort_dif2_ll        = rmissing
     allocate( rshortup_dif2_ll     (npts) ) ; rshortup_dif2_ll      = rmissing
     allocate( rlong_dif2_ll        (npts) ) ; rlong_dif2_ll         = rmissing
     allocate( rlongup_dif2_ll      (npts) ) ; rlongup_dif2_ll       = rmissing
     allocate( rshort_top_dif2_ll   (npts) ) ; rshort_top_dif2_ll    = rmissing
     allocate( rshortup_top_dif2_ll (npts) ) ; rshortup_top_dif2_ll  = rmissing
     allocate( rlongup_top_dif2_ll  (npts) ) ; rlongup_top_dif2_ll   = rmissing
  endif

  if (do3d) then
     allocate( u_ll(npts,mza-1) ) ; u_ll = rmissing
     allocate( v_ll(npts,mza-1) ) ; v_ll = rmissing
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

  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,topw,topo_ll)
  
!------------------------------------------------------------
! Compute zonal and meridional wind components on OLAM grid
! and copy their values at lowest prognosed model level to
! separate arrays
!------------------------------------------------------------

  do iw = 2, mwa
     kb = lpw(iw)

     raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

     if (raxis > 1.e3) then
        raxisi = 1. / raxis

        do k = kb, mza
           scr2a(k,iw) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
           scr2b(k,iw) = vze(k,iw) * raxis * eradi &
                       - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi
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

  if (do3d)  call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,mza,mza-1,alon,alat,scr2a,u_ll)
  if (do3d)  call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,mza,mza-1,alon,alat,scr2b,v_ll)

  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,u_sfc_ll)
  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1b,v_sfc_ll)

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

  if (do3d)  call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,mza,mza-1,alon,alat,tair,t_ll)
  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,t_sfc_ll)
  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1b,slp_ll)

!------------------------------------------------------------
! Compute vapor mixing ratio and copy its value at lowest
! prognosed model level on OLAM grid to separate array.
! Compute vapor pressure at lowest prognosed model level on OLAM grid.
!------------------------------------------------------------

  do iw = 2,mwa
     kb = lpw(iw)

     do k = kb,mza
        scr2a(k,iw) = sh_v(k,iw) / (1. - sh_w(k,iw))
     enddo

     scr1a(iw) = scr2a(kb,iw)
     scr1b(iw) = sh_v(kb,iw) * rho(kb,iw) * rvap * tair(kb,iw)
  enddo

!------------------------------------------------------------
! Interpolate vapor mixing ratio, its surface value, and
! the surface value of vapor pressure
!------------------------------------------------------------

  if (do3d)  call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,mza,mza-1,alon,alat,scr2a,r_ll)
  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,r_sfc_ll)
  if (dosfc) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1b,pvap_sfc_ll)

!------------------------------------------------------------
! Interpolate atmospheric pressure
!------------------------------------------------------------

! Copy pressure to real array

  scr2a(:,:) = real(press(:,:))

  if (do3d) call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,mza,mza-1,alon,alat,scr2a,p_ll)

!------------------------------------------------------------
! Compute total (resolved + parameterized) accumulated
! precipitation on OLAM grid
!------------------------------------------------------------

  if (doaccum) then

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

     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,accpmic_ll)

     if (allocated(aconpr)) then
        scr1a(:) = real(aconpr(:))
        call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,accpcon_ll)
     endif

     scr1a(:) = real(sfluxr_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,vapflux_accum_ll)

     scr1a(:) = real(sfluxt_accum(:)) * cp
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,sensflux_accum_ll)

     scr1a(:) = real(sfluxr_accum(:)) * alvl
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,latflux_accum_ll)

     scr1a(:) = real(rshort_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshort_accum_ll)

     scr1a(:) = real(rshortup_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshortup_accum_ll)

     scr1a(:) = real(rlong_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlong_accum_ll)

     scr1a(:) = real(rlongup_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlongup_accum_ll)

     scr1a(:) = real(rshort_top_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshort_top_accum_ll)

     scr1a(:) = real(rshortup_top_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshortup_top_accum_ll)

     scr1a(:) = real(rlongup_top_accum(:))
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlongup_top_accum_ll)

  endif

!---------------------------------------------------------------------
! Compute and interpolate MEM_PLOT MEAN precipitation rates and fluxes
!---------------------------------------------------------------------

  if (dodifs) then

     timefac = 1.0
     if (abs(time8_prev0 - time8_prev1) > .99) then
        timefac = 1.0 / (time8_prev0 - time8_prev1)
     endif

     scr1a(:) = (accpmic_prev0(:) - accpmic_prev1(:)) * timefac * 86400.
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,pcpdif2mic_ll)

     scr1a(:) = (accpcon_prev0(:) - accpcon_prev1(:)) * timefac * 86400.
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,pcpdif2con_ll)

     scr1a(:) = (accpmic_prev0(:) - accpmic_prev1(:) &
              + accpcon_prev0(:) - accpcon_prev1(:)) * timefac * 86400.
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,pcpdif2both_ll)

     scr1a(:) = (sfluxr_accum_prev0(:) - sfluxr_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,vapflux_dif2_ll)

     scr1a(:) = (sfluxt_accum_prev0(:) - sfluxt_accum_prev1(:)) * timefac * cp
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,sensflux_dif2_ll)

     scr1a(:) = (sfluxr_accum_prev0(:) - sfluxr_accum_prev1(:)) * timefac * alvl
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,latflux_dif2_ll)

     scr1a(:) = (rshort_accum_prev0(:) - rshort_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshort_dif2_ll)

     scr1a(:) = (rshortup_accum_prev0(:) - rshortup_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshortup_dif2_ll)

     scr1a(:) = (rlong_accum_prev0(:) - rlong_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlong_dif2_ll)

     scr1a(:) = (rlongup_accum_prev0(:) - rlongup_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlongup_dif2_ll)

     scr1a(:) = (rshort_top_accum_prev0(:) - rshort_top_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshort_top_dif2_ll)

     scr1a(:) = (rshortup_top_accum_prev0(:) - rshortup_top_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rshortup_top_dif2_ll)

     scr1a(:) = (rlongup_top_accum_prev0(:) - rlongup_top_accum_prev1(:)) * timefac
     call interp_htw_ll(nlon,nlat,npts,ims_loc,lls_loc,ijs_loc,1,1,alon,alat,scr1a,rlongup_top_dif2_ll)

  endif

! HDF5 write

  if (myrank == 0) write(io6,'(/,a)') "Writing lat/lon fields to disk..."

  call makefnam(hnamel, hfilepref, current_time, 'LL', '$', 'h5')
  call shdf5_open(hnamel,'W',iclobber) 

! First write coordinate variables to disk (lat, lon, height)

  ndims    = 1
  idims(2) = 1
  idims(3) = 1

  idims(1) = nlon

  CALL shdf5_orec(ndims, idims, 'lon', rvara=alon, isdim=.true., &
                  long_name = "longitude",                       &
                  standard_name = "longitude",                   &
                  units = "degrees_north"                        )

  idims(1) = nlat

  CALL shdf5_orec(ndims, idims, 'lat', rvara=alat, isdim=.true., &
                  long_name = "latitude",                        &
                  standard_name = "latitude",                    &
                  units = "degrees_east"                         )

  idims(1) = mza-1

  CALL shdf5_orec(ndims, idims, 'z', rvara=zlev, isdim=.true., &
                  long_name = "height above mean sea level",   &
                  standard_name = "altitude",                  &
                  units = "m",                                 &
                  positive = "up"                              )

  ! If we ever output data on pressure levels:

  if (dopress) then

     idims(1) = npress

     CALL shdf5_orec(ndims, idims, 'pres', rvara=plev, isdim=.true., &
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

     CALL shdf5_orec_ll(ndims, idims, 'TOPO_LL', rvara=topo_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                     &
                        long_name = "topography height",                         &
                        standard_name = "surface_altitude",                      &
                        units = "m",                                             &
                        rmissing = rmissing                                      )

     CALL shdf5_orec_ll(ndims, idims, 'U_SFC_LL', rvara=u_sfc_ll, gpoints=lls_loc,       &
                        dimnames = dimnames,                                             &
                        long_name = "eastward wind at lowest model layer above surface", &
                        standard_name = "eastward_wind",                                 &
                        units = "m s-1",                                                 &
                        rmissing = rmissing                                              )

     CALL shdf5_orec_ll(ndims, idims, 'V_SFC_LL', rvara=v_sfc_ll, gpoints=lls_loc,        &
                        dimnames = dimnames,                                              &
                        long_name = "northward wind at lowest model layer above surface", &
                        standard_name = "northward_wind",                                 &
                        units = "m s-1",                                                  &
                        rmissing = rmissing                                               )

     CALL shdf5_orec_ll(ndims, idims, 'T_SFC_LL', rvara=t_sfc_ll, gpoints=lls_loc,         &
                        dimnames = dimnames,                                               &
                        long_name = "air temperature at lowest model layer above surface", &
                        standard_name = "surface_air_temperature",                         &
                        units = "K",                                                       &
                        rmissing = rmissing                                                )

     CALL shdf5_orec_ll(ndims, idims, 'R_SFC_LL', rvara=r_sfc_ll, gpoints=lls_loc,                  &
                        dimnames = dimnames,                                                        &
                        long_name = "water vapor mixing ratio at lowest model layer above surface", &
                        standard_name = "humidity_mixing_ratio",                                    &
                        units = "kg kg-1",                                                          &
                        rmissing = rmissing                                                         )

     CALL shdf5_orec_ll(ndims, idims, 'SLP_LL', rvara=slp_ll, gpoints=lls_loc,    &
                        dimnames = dimnames,                                      &
                        long_name = "surface pressure reduced to mean sea level", &
                        standard_name = "surface_air_pressure_at_sea_level",      &
                        units = "Pa",                                             &
                        rmissing = rmissing                                       )

     CALL shdf5_orec_ll(ndims, idims, 'PVAP_SFC_LL', rvara=pvap_sfc_ll, gpoints=lls_loc,        &
                        dimnames = dimnames,                                                    &
                        long_name = "water vapor pressure at lowest model layer above surface", &
                        standard_name = "water_vapor_partial_pressure_in_air",                  &
                        units = "Pa",                                                           &
                        rmissing = rmissing                                                     )
  endif

  ! Accumulated precipitation/water vapor at surface

  if (doaccum) then

     CALL shdf5_orec_ll(ndims, idims, 'ACCPMIC_LL', rvara=accpmic_ll, gpoints=lls_loc, &
                     dimnames = dimnames,                                           &
                     long_name = "Accumulated resolved precipitation",              &
                     standard_name = "large_scale_precipitation_amount",            &
                     units = "kg m-2",                                              &
                     rmissing = rmissing                                            )

     CALL shdf5_orec_ll(ndims, idims, 'ACCPCON_LL', rvara=accpcon_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                           &
                        long_name = "Accumulated convective precipitation",            &
                        standard_name = "convective_precipitation_amount",             &
                        units = "kg m-2",                                              &
                        rmissing = rmissing                                            )

     CALL shdf5_orec_ll(ndims, idims, 'VAPFLUX_ACCUM_LL', rvara=vapflux_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                       &
                        long_name = "Accumulated water vapor at surface",                          &
                        standard_name = "integral_of_surface_upward_water_vapor_flux_wrt_time",    &
                        units = "kg m-2",                                                          &
                        rmissing = rmissing                                                        )

! Accumulated sensible/latent heat fluxes at surface

     CALL shdf5_orec_ll(ndims, idims, 'SENSFLUX_ACCUM_LL', rvara=sensflux_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                         &
                        long_name = "Accumulated upward sensible heat flux at surface",              &
                        standard_name = "integral_of_surface_upward_sensible_heat_flux_wrt_time",    &
                        units = "J m-2",                                                             &
                        rmissing = rmissing                                                          )

     CALL shdf5_orec_ll(ndims, idims, 'LATFLUX_ACCUM_LL', rvara=latflux_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                        &
                        long_name = "Accumulated upward latent heat flux at surface",               &
                        standard_name = "integral_of_surface_upward_latent_heat_flux_wrt_time",     &
                        units = "J m-2",                                                            &
                        rmissing = rmissing                                                         )

! Accumulated longwave and shortwave radiative fluxes at surface

     CALL shdf5_orec_ll(ndims, idims, 'RSHORT_ACCUM_LL', rvara=rshort_accum_ll, gpoints=lls_loc,   &
                        dimnames = dimnames,                                                       &
                        long_name = "Accumulated downwelling shortwave flux at surface",           &
                        standard_name = "integral_of_surface_downwelling_shortwave_flux_wrt_time", &
                        units = "J m-2",                                                           &
                        rmissing = rmissing                                                        )

     CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_ACCUM_LL', rvara=rshortup_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                         &
                        long_name = "Accumulated upwelling shortwave flux at surface",               &
                        standard_name = "integral_of_surface_upwelling_shortwave_flux_wrt_time",     &
                        units = "J m-2",                                                             &
                        rmissing = rmissing                                                          )
                    
     CALL shdf5_orec_ll(ndims, idims, 'RLONG_ACCUM_LL', rvara=rlong_accum_ll, gpoints=lls_loc,    &
                        dimnames = dimnames,                                                      &
                        long_name = "Accumulated downwelling longwave flux at surface",           &
                        standard_name = "integral_of_surface_downwelling_longwave_flux_wrt_time", &
                        units = "J m-2",                                                          &
                        rmissing = rmissing                                                       )

     CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_ACCUM_LL', rvara=rlongup_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                       &
                        long_name = "Accumulated upwelling longwave flux at surface",              &
                        standard_name = "integral_of_surface_upwelling_longwave_flux_wrt_time",    &
                        units = "J m-2",                                                           &
                        rmissing = rmissing                                                        )

! Accumulated longwave and shortwave radiative fluxes at top-of-atmosphere
! Note: at TOA "incoming" and "outgoing" are used in place of "downwelling" and "upwelling"

     CALL shdf5_orec_ll(ndims, idims, 'RSHORT_TOP_ACCUM_LL', rvara=rshort_top_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                             &
                        long_name = "Accumulated incoming shortwave flux at TOA",                        &
                        standard_name = "integral_of_toa_incoming_shortwave_flux_wrt_time",              &
                        units = "J m-2",                                                                 &
                        rmissing = rmissing                                                              )

     CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_TOP_ACCUM_LL', rvara=rshortup_top_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                                 &
                        long_name = "Accumulated outgoing shortwave flux at TOA",                            &
                        standard_name = "integral_of_toa_outgoing_shortwave_flux_wrt_time",                  &
                        units = "J m-2",                                                                     &
                        rmissing = rmissing                                                                  )

     CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_TOP_ACCUM_LL' , rvara=rlongup_top_accum_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                                &
                        long_name = "Accumulated outgoing longwave flux at TOA",                            &
                        standard_name = "integral_of_toa_outgoing_longwave_flux_wrt_time",                  &
                        units = "J m-2",                                                                    &
                        rmissing = rmissing                                                                 )

  endif

! Time-averaged precipitation at surface

  if (dodifs) then

     CALL shdf5_orec_ll(ndims, idims, 'PCPDIF2MIC_LL', rvara=pcpdif2mic_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                 &
                        long_name = "Time-averaged resolved precipitation rate",             &
                        standard_name = "large_scale_precipitation_rate",                    &
                        units = "mm day-1",                                                  &
                        cell_methods = "time: mean",                                         &
                        rmissing = rmissing                                                  )

     CALL shdf5_orec_ll(ndims, idims, 'PCPDIF2CON_LL', rvara=pcpdif2con_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                 &
                        long_name = "Time-averaged convective precipitation rate",           &
                        standard_name = "convective_precipitation_rate",                     &
                        units = "mm day-1",                                                  &
                        cell_methods = "time: mean",                                         &
                        rmissing = rmissing                                                  )

! Time-averaged fluxes at surface

     CALL shdf5_orec_ll(ndims, idims, 'VAPFLUX_DIF2_LL', rvara=vapflux_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                     &
                        long_name = "Time-averaged water vapor flux at surface",                 &
                        standard_name = "surface_upward_water_vapor_flux",                       &
                        units = "kg m-2 s-1",                                                    &
                        cell_methods = "time: mean",                                             &
                        rmissing = rmissing                                                      )

     CALL shdf5_orec_ll(ndims, idims, 'SENSFLUX_DIF2_LL', rvara=sensflux_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                       &
                        long_name = "Time-averaged upward sensible heat flux at surface",          &
                        standard_name = "surface_upward_sensible_heat_flux",                       &
                        units = "W m-2",                                                           &
                        cell_methods = "time: mean",                                               &
                        rmissing = rmissing                                                        )

     CALL shdf5_orec_ll(ndims, idims, 'LATFLUX_DIF2_LL', rvara=latflux_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                     &
                        long_name = "Time-averaged upward latent heat flux at surface",          &
                        standard_name = "surface_upward_latent_heat_flux",                       &
                        units = "W m-2",                                                         &
                        cell_methods = "time: mean",                                             &
                        rmissing = rmissing                                                      )

! Time-averaged longwave and shortwave radiative fluxes at surface

     CALL shdf5_orec_ll(ndims, idims, 'RSHORT_DIF2_LL', rvara=rshort_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                   &
                        long_name = "Time-averaged downwelling shortwave flux at surface",     &
                        standard_name = "surface_downwelling_shortwave_flux",                  &
                        units = "W m-2",                                                       &
                        cell_methods = "time: mean",                                           &
                        rmissing = rmissing                                                    )

     CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_DIF2_LL', rvara=rshortup_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                       &
                        long_name = "Time-averaged upwelling shortwave flux at surface",           &
                        standard_name = "surface_upwelling_shortwave_flux",                        &
                        units = "W m-2",                                                           &
                        cell_methods = "time: mean",                                               &
                        rmissing = rmissing                                                        )

     CALL shdf5_orec_ll(ndims, idims, 'RLONG_DIF2_LL', rvara=rlong_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                 &
                        long_name = "Time-averaged downwelling longwave flux at surface",    &
                        standard_name = "surface_downwelling_longwave_flux",                 &
                        units = "W m-2",                                                     &
                        cell_methods = "time: mean",                                         &
                        rmissing = rmissing                                                  )

     CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_DIF2_LL', rvara=rlongup_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                     &
                        long_name = "Time-averaged upwelling longwave flux at surface",          &
                        standard_name = "surface_upwelling_longwave_flux",                       &
                        units = "W m-2",                                                         &
                        cell_methods = "time: mean",                                             &
                        rmissing = rmissing                                                      )

! Time-averaged longwave and shortwave radiative fluxes at top-of-atmosphere
! Note: at TOA "incoming" and "outgoing" are used in place of "downwelling" and "upwelling"


     CALL shdf5_orec_ll(ndims, idims, 'RSHORT_TOP_DIF2_LL', rvara=rshort_top_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                           &
                        long_name = "Time-averaged incoming shortwave flux at TOA",                    &
                        standard_name = "toa_incoming_shortwave_flux",                                 &
                        units = "W m-2",                                                               &
                        cell_methods = "time: mean",                                                   &
                        rmissing = rmissing                                                            )

     CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_TOP_DIF2_LL', rvara=rshortup_top_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                               &
                        long_name = "Time-averaged outgoing shortwave flux at TOA",                        &
                        standard_name = "toa_outgoing_shortwave_flux",                                     &
                        units = "W m-2",                                                                   &
                        cell_methods = "time: mean",                                                       &
                        rmissing = rmissing                                                                )


     CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_TOP_DIF2_LL', rvara=rlongup_top_dif2_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                                                             &
                        long_name = "Time-averaged outgoing longwave flux at TOA",                       &
                        standard_name = "toa_outgoing_longwave_flux",                                    &
                        units = "W m-2",                                                                 &
                        cell_methods = "time: mean",                                                     &
                        rmissing = rmissing                                                              )

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

     CALL shdf5_orec_ll(ndims, idims, 'U_LL', rvara=u_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                               &
                        long_name = "eastward wind",                       &
                        standard_name = "eastward_wind",                   &
                        units = "m s-1",                                   &
                        rmissing=rmissing                                  )

     CALL shdf5_orec_ll(ndims, idims, 'V_LL', rvara=v_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                               &
                        long_name = "northward wind",                      &
                        standard_name = "northward_wind",                  &
                        units = "m s-1",                                   &
                        rmissing=rmissing                                  )

     CALL shdf5_orec_ll(ndims, idims, 'T_LL', rvara=t_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                               &
                        long_name = "air temperature",                     &
                        standard_name = "air_temperature",                 &
                        units = "K",                                       &
                        rmissing=rmissing                                  )

     CALL shdf5_orec_ll(ndims, idims, 'R_LL', rvara=r_ll, gpoints=lls_loc, &
                        dimnames = dimnames,                               &
                        long_name = "water vapor mixing ratio",            &
                        standard_name = "humidity_mixing_ratio",           &
                        units = "kg kg-1",                                 &
                        rmissing=rmissing                                  )

     CALL shdf5_orec_ll(ndims, idims, 'P_LL', rvara=p_ll, gpoints=lls_loc, &
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

! Return if we do not want to plot the interpolated fields

  if (nl%latlonplot == 0) return

!------------------------------------------------------------
! Contour plot or tile plot 2D lat-lon fields
! In parallel, contplot does not work and will default
! to tile plots
!------------------------------------------------------------

! Reopen the current graphics output workstation if it is closed

  if (myrank == 0) call o_reopnwk()
  if (myrank == 0) write(io6,'(/,a)') "Plotting interpolated lat/lon fields..."

  if (dosfc) then
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,topo_ll,402, &
                      "Topography Height", "(m)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,u_sfc_ll,113, &
                      "Zonal wind at surface", "(m/s)", rmissing)
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,v_sfc_ll,113, &
                      "Meridional wind at surface", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,t_sfc_ll,2, &
                      "Air temperature at suface", "(K)", rmissing)

     ! Preserve missing values
     where( abs(r_sfc_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) = r_sfc_ll(:) * 1000.
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,5, &
                      "Water vapor mixing ratio at suface", "(g/kg)", rmissing)

     ! Preserve missing values
     where( abs(slp_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) =  slp_ll(:) * .01 
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,17, &
                      "Sea level pressure at suface", "(mb)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pvap_sfc_ll,35, &
                      "Water vapor pressure at surface", "(mb)", rmissing)
  endif

  if (doaccum) then
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpmic_ll,54, &
                      "Accum microphysics precip", "(kg/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpcon_ll,54, &
                      "Accum convective precip", "(kg/m^2)", rmissing)

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

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,vapflux_accum_ll,304, &
                      "Accum surface vapor", "(kg/m^2)", rmissing)

     scr2_ll(:) = sensflux_accum_ll(:) * 1.e-6 
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,304, &
                      "Accum surface sensible heat", "(MJ/m^2)", rmissing)

     scr2_ll(:) = latflux_accum_ll(:) * 1.e-6 
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,304, &
                      "Accum surface latent heat", "(MJ/m^2)", rmissing)
  endif

  if (dodifs) then
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pcpdif2mic_ll,418, &
                      "Microphysics precip diff2", "(mm/day)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pcpdif2con_ll,418, &
                      "Convective precip diff2", "(mm/day)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pcpdif2both_ll,418, &
                      "Total precip diff2", "(mm/day)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_dif2_ll,63, &
                      "Surface downward s/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_dif2_ll,63, &
                      "Surface upward s/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlong_dif2_ll,63, &
                      "Surface downward l/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_dif2_ll,63, &
                      "Surface upward s/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_top_dif2_ll,63, &
                      "TOA downward s/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_top_dif2_ll,63, &
                      "TOA upward s/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_top_dif2_ll,63, &
                      "TOA upward l/w radiation diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,vapflux_dif2_ll,128, &
                      "Surface vapor flux diff2", "(kg/(m^2 s))", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,sensflux_dif2_ll,163, &
                      "Surface sensible heat flux diff2", "(W/m^2)", rmissing)

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,latflux_dif2_ll,163, &
                      "Surface latent heat flux diff2", "(W/m^2)", rmissing)

  endif

!------------------------------------------------------------
! Contour plot or tile plot 2D slices of 3D lat-lon fields
!------------------------------------------------------------

  if (do3d) then
     call tileplot_ll(nlon,nlat,mza-1,1,alon,alat,npts,lls_loc,u_ll,113, &
                      "Zonal wind at lowest level", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,mza-1,1,alon,alat,npts,lls_loc,v_ll,113, &
                      "Meridional wind at lowest level", "(m/s)", rmissing)

     call tileplot_ll(nlon,nlat,mza-1,1,alon,alat,npts,lls_loc,t_ll,2, &
                      "Temperature at lowest level", "(K)", rmissing)

     ! Preserve missing values
     where( abs(r_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = r_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,1,alon,alat,npts,lls_loc,scr3_ll,5, &
          "Water vapor mixing ratio at lowest level", "(g/kg)", rmissing)

     ! Preserve missing values
     where( abs(p_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = p_ll(:,:) * 0.01
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,1,alon,alat,npts,lls_loc,scr3_ll,8, &
          "Air Pressure at lowest level", "(mb)", rmissing)

  endif

!---------------------------------------------------------------------------
! Contour plot 1D arrays that are longitudinal averages of 2D lat-lon fields
! This will only work for serial runs that have the entire lat/lon memory
!---------------------------------------------------------------------------

  if (iparallel == 1 .or. npts /= nlon*nlat) then
     ! For a parallel run, close graphics and return to program
     if (myrank == 0) call o_clswk()
     return
  endif

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

  if (dodifs) then

     arr_ll = reshape( pcpdif2mic_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,8,0,   &
                    nlat, alat, value,             &
                    'latitude','PCPDIF2 (mm/day)', &
                    alat(1),alat(nlat),alatinc,3,  &
                    ymin, ymax, yinc, 5            )

     arr_ll = reshape( pcpdif2con_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,1,0,  &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape( pcpdif2both_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('1','N',aspect,scalelab,11,0, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

!---------------------------------------------------------------

  ymin = -10.
  ymax = 500.
  yinc = 10.

  if (dodifs) then

     arr_ll = reshape( rshort_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,8,0,  &
                  nlat, alat, value,            &
                  'latitude','SFC RAD (W/m^2)', &
                  alat(1),alat(nlat),alatinc,3, &
                  ymin, ymax, yinc, 5           )

     arr_ll = reshape( rshortup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape( rlong_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,1,0,  &
                    nlat, alat, value,            &
                    'latitude','',                &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape( rlongup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif
!---------------------------------------------------------------

  if (dodifs) then

     arr_ll = reshape( rshort_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,8,0,   &
                    nlat, alat, value,             &
                    'latitude','TOA RAD (W/m^2) ', &
                    alat(1),alat(nlat),alatinc,3,  &
                    ymin, ymax, yinc, 5            )

     arr_ll = reshape( rshortup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )

     arr_ll = reshape( rlongup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif
!---------------------------------------------------------------

  ymin = -50.
  ymax = 200.
  yinc = 10.

  if (dodifs) then

     arr_ll = reshape( sensflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N',aspect,scalelab,17,0,      &
                    nlat, alat, value,                 &
                    'latitude','SFC S/L FLUX (W/m^2)', &
                    alat(1),alat(nlat),alatinc,3,      &
                    ymin, ymax, yinc, 5                )

     arr_ll = reshape( latflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N',aspect,scalelab,9,0,  &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif
!---------------------------------------------------------------

  call o_frame()
  call o_clswk()

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

     call oplot_xy2(panel, cbar, 0., .016, 10, 0,         &
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

subroutine contplot_ll(nlon,nlat,nlev,ilev,alon,alat,npts,lls_loc,fld,itab)
  use plotcolors, only: clrtab
  use oplot_coms, only: op
  use misc_coms,  only: iparallel, io6

  implicit none

  integer, intent(in) :: nlon,nlat,nlev,ilev,npts,itab
  integer, intent(in) :: lls_loc(npts)
  real,    intent(in) :: alon(nlon),alat(nlat)
  real,    intent(in) :: fld(npts,nlev)

  real :: htpn(4),vtpn(4),fldvals(4)
  real :: glon(nlon)
  integer :: ilat, ilon
  integer :: ll1, ll2, ll3, ll4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This routine does not work in parallel unless we implement MPI communication
  ! of lat/lon fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (iparallel == 1) then
     write(io6,*) "Contour plotting of lat/lon fields does not work in parallel."
     return
  endif

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

        ! converd 2d lat/lon indices to 1d

        ll1 = ilon     + (ilat-1) * nlon
        ll2 = ilon + 1 + (ilat-1) * nlon
        ll3 = ilon + 1 + (ilat  ) * nlon
        ll4 = ilon     + (ilat  ) * nlon

        ! Set field value and extract contour color from color table
        
        ! This currently only works if the node has data for 
        ! every lat/lon point (i.e., serial run)!

        fldvals(1) = fld(ll1, ilev)
        fldvals(2) = fld(ll2, ilev)
        fldvals(3) = fld(ll3, ilev)
        fldvals(4) = fld(ll4, ilev)

        call contpolyg(itab,1,4,htpn,vtpn,fldvals)
        call contpolyg(itab,0,4,htpn,vtpn,fldvals)
     enddo
  enddo

  call o_frame()
end subroutine contplot_ll

!==========================================================================

subroutine mkmap_ll()

  use oplot_coms,  only: op
  use consts_coms, only: erad, erad2
  use misc_coms,   only: io6
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
