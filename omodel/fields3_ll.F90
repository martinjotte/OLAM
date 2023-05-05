subroutine fields3_ll()

! This subroutine interpolates selected fields from the OLAM grid to a
! structured lat-lon grid of either limited or global area.

! The following tasks are performed:

! (1) Allocate arrays for interpolated fields
! (2) Allocate scratch arrays
! (3) Copy each selected model field to scratch array, or compute it from
!     model fields if required
! (4) Call subroutine to interpolate each selected field to lat-lon array
! (5) Write lat-lon arrays to hdf5 file
! (6) Plot lat-lon arrays
! (7) Plot longitudinal averages of some lat-lon arrays

! New fields may be added following the examples in this file.

  use mem_ijtabs,  only: itab_w, itabg_w, jtab_w, jtw_prog
  use mem_basic,   only: wc, rho, press, theta, rr_w, rr_v, tair, ue, ve
  use mem_grid,    only: mza, mwa, lpw, topw, zt, dzt, &
                         vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_land,    only: land, omland, dslz
  use leaf_coms,   only: isfcl
  use misc_coms,   only: io6, current_time, hfilepref, iclobber, iparallel, &
                         mdomain, time8, nxp, ngrids
  use consts_coms, only: p00, p00i, rocp, cp, alvl, rvap, r8, grav, eps_virt
  use hdf5_utils,  only: shdf5_open, shdf5_orec, shdf5_orec_ll, shdf5_close, &
                         shdf5_write_global_attribute
  use max_dims,    only: pathlen, maxlatlon
  use mem_addsc,   only: addsc
  use oname_coms,  only: nl
  use mem_para,    only: myrank, mgroupsize

#ifdef OLAM_MPI
  use mpi
#endif

! NOTE: Fields from the MEM_MICRO (acc fields only), MEM_CUPARM, and MEM_FLUX_ACCUM
! modules are time-integrated fluxes of mass or energy and have units of
! [kg/m^2] or [J/m^2].  Their values written to any history file represent flux
! integrals from the beginning of a simulation to the time of the history file
! write.  Thus, time averages of a flux over any time interval can be computed
! as the difference between the flux integrals in two different history files
! divided by the difference in simulation times when the files were written.

  use mem_micro,  only: accpd, accpr, accpp, accps, accpa, accpg, accph, &
                        pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                        rr_c, rr_r

  use mem_cuparm, only: aconpr

  use mem_flux_accum, only:   rshort_accum,         rshortup_accum, &
                               rlong_accum,          rlongup_accum, &
                          rshort_top_accum,     rshortup_top_accum, &
                         rlongup_top_accum,                         &
                                  vc_accum,               wc_accum, &
                               press_accum,             tair_accum, &
                                rr_v_accum,                         &
  ! SFC grid variables

                                vels_accum,                         &
                             airtemp_accum,           airrrv_accum, &
                             cantemp_accum,           canrrv_accum, &
                            skintemp_accum,                         &
                              sfluxt_accum,           sfluxr_accum, &
                              wxferi_accum,           wxferp_accum, &
                              wxfer1_accum

  implicit none

!--------------------------------------------------------------------------------
! THE USER MUST SPECIFY THE FOLLOWING PARAMETERS THAT DEFINE THE OUTPUT
! LATITUDE-LONGITUDE ARRAYS. LONGITUDE VALUES ARE DEFINED TO BE IN THE
! RANGE (-180,180], WHICH DOES NOT PRECLUDE THAT THE LIMITED-AREA
! LATITUDE-LONGITUDE GRID CROSS THE 180 DEGREE MERIDIAN.
!--------------------------------------------------------------------------------

! NLON and NLAT are the number of longitude and latitude points in the lat-lon
! arrays that get filled by interpolation from the native OLAM grid.

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

!----------------------------
! ATM 3D INSTANTANEOUS FIELDS
!----------------------------

  real, allocatable ::  u_ll(:,:) ! zonal wind component (m/s)
  real, allocatable ::  v_ll(:,:) ! meridional wind component (m/s)
  real, allocatable ::  w_ll(:,:) ! vertical wind component (m/s)
  real, allocatable ::  t_ll(:,:) ! air temperature (K)
  real, allocatable :: th_ll(:,:) ! air potential temperature (K)
  real, allocatable :: rv_ll(:,:) ! water vapor mixing ratio (kg/kg_dryair)
  real, allocatable :: rc_ll(:,:) ! cloud water mixing ratio (kg/kg_dryair)
  real, allocatable :: rr_ll(:,:) ! rain mixing ratio (kg/kg_dryair)
  real, allocatable ::  p_ll(:,:) ! air pressure (Pa)
  real, allocatable :: q1_ll(:,:) ! added scalar #1
  real, allocatable :: q2_ll(:,:) ! added scalar #2
  real, allocatable :: q3_ll(:,:) ! added scalar #3
  real, allocatable :: q4_ll(:,:) ! added scalar #4
  real, allocatable :: q5_ll(:,:) ! added scalar #5

!----------------------------
! ATM 2D INSTANTANEOUS FIELDS
!----------------------------

  real, allocatable :: u_lpw_ll   (:) ! lpw zonal wind component (m/s)
  real, allocatable :: v_lpw_ll   (:) ! lpw meridional wind component (m/s)
  real, allocatable :: t_lpw_ll   (:) ! lpw air temperature (K)
  real, allocatable :: rv_lpw_ll  (:) ! lpw water vapor mixing ratio (kg/kg_dryair)
  real, allocatable :: pvap_lpw_ll(:) ! lpw vapor pressure (Pa)
  real, allocatable :: slp_ll     (:) ! sea level pressure (Pa)
  real, allocatable :: pcprmic_ll (:) ! microphysics precip rate (kg/(m^2 s))
  real, allocatable :: topo_ll    (:) ! topography height (m)
  real, allocatable :: q1zint_ll  (:) ! vertically integrated scalar #1
  real, allocatable :: q2zint_ll  (:) ! vertically integrated scalar #2
  real, allocatable :: qyzint_ll  (:) ! vertically integrated (scalar #1 + 2 * scalar #2)

!-------------------------------
! ATM 3D TIME-ACCUMULATED FIELDS
!-------------------------------

  real, allocatable ::  u_accum_ll(:,:) ! zonal wind component accum (m/s)
  real, allocatable ::  v_accum_ll(:,:) ! meridional wind component accum (m/s)
  real, allocatable ::  w_accum_ll(:,:) ! vertical wind component accum (m/s)
  real, allocatable ::  t_accum_ll(:,:) ! air temperature accum (K)
  real, allocatable :: rv_accum_ll(:,:) ! water vapor mixing ratio accum (kg/kg_dryair)
  real, allocatable ::  p_accum_ll(:,:) ! air pressure accum (Pa)

!-------------------------------
! ATM 2D TIME-ACCUMULATED FIELDS
!-------------------------------

  real, allocatable :: rshort_accum_ll      (:) ! accum sfc downward s/w rad (J/m^2)
  real, allocatable :: rshortup_accum_ll    (:) ! accum sfc upward s/w rad (J/m^2)
  real, allocatable :: rlong_accum_ll       (:) ! accum sfc downward l/w rad (J/m^2)
  real, allocatable :: rlongup_accum_ll     (:) ! accum sfc upward l/w rad (J/m^2)
  real, allocatable :: rshort_top_accum_ll  (:) ! accum TOA downward s/w rad (J/m^2)
  real, allocatable :: rshortup_top_accum_ll(:) ! accum TOA upward s/w rad (J/m^2)
  real, allocatable :: rlongup_top_accum_ll (:) ! accum TOA upward l/w rad (J/m^2)
  real, allocatable :: accpmic_ll           (:) ! accum microphysics precip (kg/m^2)
  real, allocatable :: accpcon_ll           (:) ! accum convective precip (kg/m^2)

!-----------------------------------------------
! SURFACE-GRID FIELDS AVERAGED TO ATM GRID CELLS
!-----------------------------------------------

  real, allocatable :: asfc_u10m_ll           (:) ! asfc 10m u wind (m/s)
  real, allocatable :: asfc_v10m_ll           (:) ! asfc 10m v wind (m/s)
  real, allocatable :: asfc_speed10m_ll       (:) ! asfc 10m wind speed (m/s)
  real, allocatable :: asfc_u2m_ll            (:) ! asfc 2m u wind (m/s)
  real, allocatable :: asfc_v2m_ll            (:) ! asfc 2m v wind (m/s)
  real, allocatable :: asfc_speed2m_ll        (:) ! asfc 2m wind speed (m/s)
  real, allocatable :: asfc_temp2m_ll         (:) ! asfc 2m temp (K)
  real, allocatable :: asfc_rvap2m_ll         (:) ! asfc 2m vapor mixing ratio (kg/kg_dryair)
  real, allocatable :: asfc_vels_accum_ll     (:) ! asfc wind speed accum (m)
  real, allocatable :: asfc_airtempk_accum_ll (:) ! asfc atm temp accum (K s)
  real, allocatable :: asfc_airrv_accum_ll    (:) ! asfc atm vapor mixing ratio accum (kg s/kg_dryair)
  real, allocatable :: asfc_cantempk_accum_ll (:) ! asfc canopy air temp accum (K s)
  real, allocatable :: asfc_canrv_accum_ll    (:) ! asfc canopy air vap mixing ratio accum (kg s/kg_dryair)
  real, allocatable :: asfc_skintempk_accum_ll(:) ! asfc canopy air temp accum (K s)
  real, allocatable :: asfc_sensflux_accum_ll (:) ! asfc sens heat flux accum (J/m^2)
  real, allocatable :: asfc_latflux_accum_ll  (:) ! asfc lat heat flux accum (J/m^2)
  real, allocatable :: asfc_vapflux_accum_ll  (:) ! asfc vap flux accum (kg/m^2)
  real, allocatable :: al_wxferi_accum_ll     (:) ! al infiltration flux accum (m)
  real, allocatable :: al_wxferp_accum_ll     (:) ! al percolation flux accum (m)
  real, allocatable :: al_wxfer1_accum_ll     (:) ! al soil bottom water flux accum (m)
  real, allocatable :: al_sfcwater_tot_ll     (:) ! al total sfcwater mass (kg/m^2)
  real, allocatable :: al_soil_water_tot_ll   (:) ! al total soil water (m)

!----------------------------------------------------------------
! ATM 2D TIME-ACCUMULATED FIELDS, DIFFERENCES OVER TIME INTERVALS
!----------------------------------------------------------------

  real, allocatable, save :: pcpmic_dif2_ll       (:) ! Differences between consecutive calls
  real, allocatable, save :: pcpcon_dif2_ll       (:) ! to this subroutine of various accumulated
  real, allocatable, save :: pcpboth_dif2_ll      (:) ! quantities.  Division by intervening time
  real, allocatable, save :: rshort_dif2_ll       (:) ! interval gives mean rates of accumulation
  real, allocatable, save :: rshortup_dif2_ll     (:) ! over that interval.  Often, this subroutine
  real, allocatable, save :: rlong_dif2_ll        (:) ! is called in a PLOTONLY run, once for each
  real, allocatable, save :: rlongup_dif2_ll      (:) ! history file that is read.
  real, allocatable, save :: rshort_top_dif2_ll   (:)
  real, allocatable, save :: rshortup_top_dif2_ll (:)
  real, allocatable, save :: rlongup_top_dif2_ll  (:)
  real, allocatable, save :: asfc_sensflux_dif2_ll(:)
  real, allocatable, save :: asfc_latflux_dif2_ll (:)

!------------------------
! GRID AND SCRATCH ARRAYS
!------------------------

  real, allocatable, save ::  alat(:) ! latitudes of grid points (deg)
  real, allocatable, save ::  alon(:) ! longitudes of grid points (deg)
  real, allocatable, save ::  zlev(:) ! heights above sea level of grid points (m)
  real, allocatable       :: value(:) ! longitudinal average of a field
  real, allocatable       :: arr_ll(:,:)

  integer       :: k, iw, ilat, ilon, n, kb, ier, np, j, i, iwsfc, jasfc, kwatm
  integer       :: iland, iv, jv, npoly, klev
  integer       :: ndims, idims(3)
  character(30) :: dimnames(3), fld_ll

  real :: area_land_sum
  real :: vx, vy, vz
  real :: zobs, press_zobs, exner_zobs, wind_zobs, theta_zobs, rrv_zobs
  real :: canexner, cantheta, canthetav, airthetav, tstar, rstar, ufree
  real :: frac, uwind, vwind

  real :: scr1a(mwa), scr1b(mwa), scr1c(mwa), scr1d(mwa)
  real :: scr1e(mwa), scr1f(mwa), scr1g(mwa), scr1h(mwa)
  real :: scr1i(mwa), scr1j(mwa), scr1k(mwa), scr1l(mwa)
  real :: scr1m(mwa), scr1n(mwa), scr1o(mwa), scr1p(mwa)
  real :: scr1q(mwa), scr1r(mwa), scr1s(mwa), scr1t(mwa)
  real :: scr1u(mwa), scr1v(mwa)

  real :: scr2a(mza,mwa), scr2b(mza,mwa)
  real, allocatable :: scr2_ll(:), scr3_ll(:,:)

  real :: timedifi
  real :: aspect, scalelab, ymin, ymax, yinc, alatinc

! SEIGEL 2013 - Added for ll interp writeout
  character(pathlen) :: hnamel
  character(20)      :: ofrq
  character(3)       :: levs

  integer,  save :: npts
  real,     save :: dlon, dlat
  real(r8), save :: time8_prev
  logical,  save :: first_call = .true.

  integer, allocatable :: iws_l (:)
  integer, allocatable :: iws_gl(:,:)
  integer, allocatable :: iws_ll(:,:,:)
  real,    allocatable :: wts_ll(:,:,:)

  integer, allocatable, save :: iws_loc(:,:), lls_loc(:)
  real,    allocatable, save :: wts_loc(:,:)

  ! These switches control which lat/lon interpolations are performed

  logical, parameter :: dopress = .false. ! TODO!
  logical, parameter :: latplot = .false.

!------------------------------------------------------
! 3D FIELDS - interpolation to pressure levels (TODO!)
!------------------------------------------------------

  real,    parameter :: onethird = 1./3.
  integer, parameter :: npress = 8
  real :: plev(npress) = (/ 1000., 925., 850., 700., 500., 300., 200., 100. /)

  if (nl%ioutput_latlon /= 1 .and. nl%latlonplot /= 1) return

  if (first_call) then

     if (myrank == 0) write(io6,'(/,a)') "Initializing overlaps for lat/lon outputs..."

     ! Compute latitude and longitude of output grid points (assuming uniform spacing)

     nlat = nint((nl%endlat - nl%beglat) * real(nl%nlatlonpd)) + 1  ! # of lat values
     if (nlat < 2) stop 'stop: nlat < 2 in fields3_ll '
     dlat = (nl%endlat  - nl%beglat) / real(nlat-1)

     allocate(alat(nlat))

     do ilat = 1, nlat
        alat(ilat) = nl%beglat + dlat * real(ilat-1)
     enddo

     if (nl%endlon > nl%beglon) then
        nlon = nint((nl%endlon - nl%beglon) * real(nl%nlatlonpd)) + 1  ! # of lon values
        if (nlon < 2) stop 'stop: nlon < 2 in fields3_ll '
        dlon = (nl%endlon - nl%beglon) / real(nlon-1)
     else
        nlon = nint((nl%endlon + 360. - nl%beglon) * nl%nlatlonpd) + 1  ! # of lon values
        if (nlon < 2) stop 'stop: nlon < 2 in fields3_ll '
        dlon = (nl%endlon + 360. - nl%beglon) / real(nlon-1)
     endif

     allocate(alon(nlon))

     do ilon = 1, nlon
        alon(ilon) = nl%beglon + dlon * real(ilon-1)
        if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.
     enddo

     ! Print lat-lon grid information

     if (myrank == 0) then
        write(io6,'(/,a)')          'fields3_ll lat-lon grid information '
        write(io6,'(/,a,3f9.3,i5)') 'beglat,endlat,dlat,nlat ',nl%beglat,nl%endlat,dlat,nlat
        write(io6,'(/,a,3f9.3,i5)') 'beglon,endlon,dlon,nlon ',nl%beglon,nl%endlon,dlon,nlon
     endif

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

  endif

  ! Selectively allocate _ll arrays according to specification in OLAMIN namelist

  do j = 1,maxlatlon
     fld_ll = nl%latlon_vars(j)

     select case(fld_ll)

     case('')
        exit
     case('U_LL')
        allocate(  u_ll(npts,mza-1) ) ;  u_ll = rmissing
     case('V_LL')
        allocate(  v_ll(npts,mza-1) ) ;  v_ll = rmissing
     case('W_LL')
        allocate(  w_ll(npts,mza-1) ) ;  w_ll = rmissing
     case('T_LL')
        allocate(  t_ll(npts,mza-1) ) ;  t_ll = rmissing
     case('TH_LL')
        allocate( th_ll(npts,mza-1) ) ; th_ll = rmissing
     case('RV_LL')
        allocate( rv_ll(npts,mza-1) ) ; rv_ll = rmissing
     case('RC_LL')
        allocate( rc_ll(npts,mza-1) ) ; rc_ll = rmissing
     case('RR_LL')
        allocate( rr_ll(npts,mza-1) ) ; rr_ll = rmissing
     case('P_LL')
        allocate(  p_ll(npts,mza-1) ) ;  p_ll = rmissing
     case('Q1_LL')
        if (nl%naddsc >= 1) then
           allocate( q1_ll(npts,mza-1) ) ; q1_ll = rmissing
        endif
     case('Q2_LL')
        if (nl%naddsc >= 2) then
           allocate( q2_ll(npts,mza-1) ) ; q2_ll = rmissing
        endif
     case('Q3_LL')
        if (nl%naddsc >= 3) then
           allocate( q3_ll(npts,mza-1) ) ; q3_ll = rmissing
        endif
     case('Q4_LL')
        if (nl%naddsc >= 4) then
           allocate( q4_ll(npts,mza-1) ) ; q4_ll = rmissing
        endif
     case('Q5_LL')
        if (nl%naddsc >= 5) then
           allocate( q5_ll(npts,mza-1) ) ; q5_ll = rmissing
        endif

     case('U_LPW_LL')
        allocate( u_lpw_ll(npts) ) ; u_lpw_ll = rmissing
     case('V_LPW_LL')
        allocate( v_lpw_ll(npts) ) ; v_lpw_ll = rmissing
     case('T_LPW_LL')
        allocate( t_lpw_ll(npts) ) ; t_lpw_ll = rmissing
     case('RV_LPW_LL')
        allocate( rv_lpw_ll(npts) ) ; rv_lpw_ll = rmissing
     case('PVAP_LPW_LL')
        allocate( pvap_lpw_ll(npts) ) ; pvap_lpw_ll = rmissing
     case('SLP_LL')
        allocate( slp_ll(npts) ) ; slp_ll = rmissing
     case('PCPRMIC_LL')
        allocate( pcprmic_ll           (npts) ) ; pcprmic_ll            = rmissing
     case('TOPO_LL')
        allocate(  topo_ll(npts) ) ; topo_ll = rmissing

     case('U_ACCUM_LL')
        allocate(  u_accum_ll(npts,mza-1) ) ;  u_accum_ll = rmissing
     case('V_ACCUM_LL')
        allocate(  v_accum_ll(npts,mza-1) ) ;  v_accum_ll = rmissing
     case('W_ACCUM_LL')
        allocate(  w_accum_ll(npts,mza-1) ) ;  w_accum_ll = rmissing
     case('T_ACCUM_LL')
        allocate(  t_accum_ll(npts,mza-1) ) ;  t_accum_ll = rmissing
     case('RV_ACCUM_LL')
        allocate( rv_accum_ll(npts,mza-1) ) ; rv_accum_ll = rmissing
     case('P_ACCUM_LL')
        allocate(  p_accum_ll(npts,mza-1) ) ;  p_accum_ll = rmissing

     case('RSHORT_ACCUM_LL')
        allocate( rshort_accum_ll      (npts) ) ; rshort_accum_ll       = rmissing
     case('RSHORTUP_ACCUM_LL')
        allocate( rshortup_accum_ll    (npts) ) ; rshortup_accum_ll     = rmissing
     case('RLONG_ACCUM_LL')
        allocate( rlong_accum_ll       (npts) ) ; rlong_accum_ll        = rmissing
     case('RLONGUP_ACCUM_LL')
        allocate( rlongup_accum_ll     (npts) ) ; rlongup_accum_ll      = rmissing
     case('RSHORT_TOP_ACCUM_LL')
        allocate( rshort_top_accum_ll  (npts) ) ; rshort_top_accum_ll   = rmissing
     case('RSHORTUP_TOP_ACCUM_LL')
        allocate( rshortup_top_accum_ll(npts) ) ; rshortup_top_accum_ll = rmissing
     case('RLONGUP_TOP_ACCUM_LL')
        allocate( rlongup_top_accum_ll (npts) ) ; rlongup_top_accum_ll  = rmissing
     case('ACCPMIC_LL')
        allocate( accpmic_ll           (npts) ) ; accpmic_ll            = rmissing
     case('ACCPCON_LL')
        allocate( accpcon_ll           (npts) ) ; accpcon_ll            = rmissing

     case('ASFC_U10M_LL')
        allocate( asfc_u10m_ll           (npts) ) ; asfc_u10m_ll            = rmissing
     case('ASFC_V10M_LL')
        allocate( asfc_v10m_ll           (npts) ) ; asfc_v10m_ll            = rmissing
     case('ASFC_SPEED10M_LL')
        allocate( asfc_speed10m_ll       (npts) ) ; asfc_speed10m_ll        = rmissing
     case('ASFC_U2M_LL')
        allocate( asfc_u2m_ll            (npts) ) ; asfc_u2m_ll             = rmissing
     case('ASFC_V2M_LL')
        allocate( asfc_v2m_ll            (npts) ) ; asfc_v2m_ll             = rmissing
     case('ASFC_SPEED2M_LL')
        allocate( asfc_speed2m_ll        (npts) ) ; asfc_speed2m_ll         = rmissing
     case('ASFC_TEMP2M_LL')
        allocate( asfc_temp2m_ll         (npts) ) ; asfc_temp2m_ll          = rmissing
     case('ASFC_RVAP2M_LL')
        allocate( asfc_rvap2m_ll         (npts) ) ; asfc_rvap2m_ll          = rmissing
     case('ASFC_VELS_ACCUM_LL')
        allocate( asfc_vels_accum_ll     (npts) ) ; asfc_vels_accum_ll      = rmissing
     case('ASFC_AIRTEMPK_ACCUM_LL')
        allocate( asfc_airtempk_accum_ll (npts) ) ; asfc_airtempk_accum_ll  = rmissing
     case('ASFC_AIRRV_ACCUM_LL')
        allocate( asfc_airrv_accum_ll    (npts) ) ; asfc_airrv_accum_ll     = rmissing
     case('ASFC_CANTEMPK_ACCUM_LL')
        allocate( asfc_cantempk_accum_ll (npts) ) ; asfc_cantempk_accum_ll  = rmissing
     case('ASFC_CANRV_ACCUM_LL')
        allocate( asfc_canrv_accum_ll    (npts) ) ; asfc_canrv_accum_ll     = rmissing
     case('ASFC_SKINTEMPK_ACCUM_LL')
        allocate( asfc_skintempk_accum_ll(npts) ) ; asfc_skintempk_accum_ll = rmissing
     case('ASFC_SENSFLUX_ACCUM_LL')
        allocate( asfc_sensflux_accum_ll (npts) ) ; asfc_sensflux_accum_ll  = rmissing
     case('ASFC_LATFLUX_ACCUM_LL')
        allocate( asfc_latflux_accum_ll  (npts) ) ; asfc_latflux_accum_ll   = rmissing
     case('ASFC_VAPFLUX_ACCUM_LL')
        allocate( asfc_vapflux_accum_ll  (npts) ) ; asfc_vapflux_accum_ll   = rmissing
     case('AL_WXFERI_ACCUM_LL')
        allocate( al_wxferi_accum_ll     (npts) ) ; al_wxferi_accum_ll      = rmissing
     case('AL_WXFERP_ACCUM_LL')
        allocate( al_wxferp_accum_ll     (npts) ) ; al_wxferp_accum_ll      = rmissing
     case('AL_WXFER1_ACCUM_LL')
        allocate( al_wxfer1_accum_ll     (npts) ) ; al_wxfer1_accum_ll      = rmissing
     case('AL_SFCWATER_TOT_LL')
        allocate( al_sfcwater_tot_ll     (npts) ) ; al_sfcwater_tot_ll      = rmissing
     case('AL_SOIL_WATER_TOT_LL')
        allocate( al_soil_water_tot_ll   (npts) ) ; al_soil_water_tot_ll    = rmissing
     end select
  enddo

  ! 2nd round of allocations, which depend on allocations in the 1st round

  do j = 1,maxlatlon
     fld_ll = nl%latlon_vars(j)

     select case(fld_ll)

     case('')
        exit
     case('Q1ZINT_LL')
        if (allocated(q1_ll)) then
           allocate(q1zint_ll(npts))               ; q1zint_ll = rmissing
        endif
     case('Q2ZINT_LL')
        if (allocated(q2_ll)) then
           allocate(q2zint_ll(npts))               ; q2zint_ll = rmissing
        endif
     case('QYZINT_LL')
        if (allocated(q1_ll) .and. allocated(q2_ll)) then
           allocate(qyzint_ll(npts))               ; qyzint_ll = rmissing
        endif
     case('PCPMIC_DIF2_LL')
        if (allocated(accpmic_ll) .and. first_call) then
           allocate( pcpmic_dif2_ll(npts) )       ; pcpmic_dif2_ll = 0.
        endif
     case('PCPCON_DIF2_LL')
        if (allocated(accpcon_ll) .and. first_call) then
           allocate( pcpcon_dif2_ll(npts) )       ; pcpcon_dif2_ll = 0.
        endif
     case('PCPBOTH_DIF2_LL')
        if (allocated(pcpmic_dif2_ll) .and. first_call .and. &
            allocated(pcpcon_dif2_ll)) then
           allocate( pcpboth_dif2_ll(npts) )      ; pcpboth_dif2_ll = 0.
        endif
     case('RSHORT_DIF2_LL')
        if (allocated(rshort_accum_ll) .and. first_call) then
           allocate( rshort_dif2_ll(npts) )       ; rshort_dif2_ll = 0.
        endif
     case('RSHORTUP_DIF2_LL')
        if (allocated(rshortup_accum_ll) .and. first_call) then
           allocate( rshortup_dif2_ll(npts) )     ; rshortup_dif2_ll = 0.
        endif
     case('RLONG_DIF2_LL')
        if (allocated(rlong_accum_ll) .and. first_call) then
           allocate( rlong_dif2_ll(npts) )        ; rlong_dif2_ll = 0.
        endif
     case('RLONGUP_DIF2_LL')
        if (allocated(rlongup_accum_ll) .and. first_call) then
           allocate( rlongup_dif2_ll(npts) )      ; rlongup_dif2_ll = 0.
        endif
     case('RSHORT_TOP_DIF2_LL')
        if (allocated(rshort_top_accum_ll) .and. first_call) then
           allocate( rshort_top_dif2_ll(npts) )   ; rshort_top_dif2_ll = 0.
        endif
     case('RSHORTUP_TOP_DIF2_LL')
        if (allocated(rshortup_top_accum_ll) .and. first_call) then
           allocate( rshortup_top_dif2_ll(npts) ) ; rshortup_top_dif2_ll = 0.
        endif
     case('RLONGUP_TOP_DIF2_LL')
        if (allocated(rlongup_top_accum_ll) .and. first_call) then
           allocate( rlongup_top_dif2_ll(npts) )  ; rlongup_top_dif2_ll = 0.
        endif
     case('ASFC_SENSFLUX_DIF2_LL')
        if (allocated(asfc_sensflux_accum_ll) .and. first_call) then
           allocate( asfc_sensflux_dif2_ll(npts) ); asfc_sensflux_dif2_ll = 0.
        endif
     case('ASFC_LATFLUX_DIF2_LL')
        if (allocated(asfc_latflux_accum_ll) .and. first_call) then
           allocate( asfc_latflux_dif2_ll(npts) ) ; asfc_latflux_dif2_ll = 0.
        endif
     end select
  enddo

  first_call = .false.

  ! scratch arrays
  allocate( scr2_ll(npts) )
  allocate( scr3_ll(npts,mza-1) )

  if (myrank == 0) write(io6,'(/,a)') "Interpolating fields to lat/lon grid..."

  !------------------------------------------------------------------------
  ! Copy some atmospheric quantities to scratch arrays before interpolation
  !------------------------------------------------------------------------

  do iw = 2, mwa
     kb = lpw(iw)

     scr1a(iw) = ue(kb,iw)
     scr1b(iw) = ve(kb,iw)
     scr1c(iw) = tair(kb,iw)
     scr1d(iw) = press(kb,iw) &
               * (1. - .0065 * zt(kb) / (tair(kb,iw) + .0065 * zt(kb)))**(-5.257)
     scr1e(iw) = rr_v(kb,iw)
     scr1f(iw) = rr_v(kb,iw) * rho(kb,iw) * rvap * tair(kb,iw)

     do k = kb,mza
        scr2a(k,iw) = 0.5 * (wc(k,iw) + wc(k-1,iw))
        scr2b(k,iw) = real(press(k,iw))
     enddo
  enddo

  !--------------------------------------------------------------------------------
  ! Interpolate 3D atmospheric fields
  !--------------------------------------------------------------------------------

  if (allocated( u_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,ue           , u_ll)
  if (allocated( v_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,ve           , v_ll)
  if (allocated( w_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a        , w_ll)
  if (allocated( t_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,tair         , t_ll)
  if (allocated(th_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,theta        ,th_ll)
  if (allocated(rv_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,rr_v         ,rv_ll)
  if (allocated(rc_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,rr_c         ,rc_ll)
  if (allocated(rr_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,rr_r         ,rr_ll)
  if (allocated( p_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b        , p_ll)
  if (allocated(q1_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,addsc(1)%sclp,q1_ll)
  if (allocated(q2_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,addsc(2)%sclp,q2_ll)
  if (allocated(q3_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,addsc(3)%sclp,q3_ll)
  if (allocated(q4_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,addsc(4)%sclp,q4_ll)
  if (allocated(q5_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,addsc(5)%sclp,q5_ll)

  !------------------------------------------------------------------------------------------
  ! Interpolate 2D atmospheric fields
  !------------------------------------------------------------------------------------------

  if (allocated(   u_lpw_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,   u_lpw_ll)
  if (allocated(   v_lpw_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,   v_lpw_ll)
  if (allocated(   t_lpw_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1c,   t_lpw_ll)
  if (allocated(  rv_lpw_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1e,  rv_lpw_ll)
  if (allocated(pvap_lpw_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1f,pvap_lpw_ll)
  if (allocated(     slp_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1d,     slp_ll)
  if (allocated(    topo_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,topw ,    topo_ll)

  !----------------------------------------------------------------
  ! Compute vertical averages of added scalars and interpolate them
  !----------------------------------------------------------------

  if (allocated(q1zint_ll)) then
     do iw = 2, mwa
        scr1a(iw) = 0.
        scr1b(iw) = 0.
        do k = lpw(iw),mza
           scr1a(iw) = scr1a(iw) + addsc(1)%sclp(k,iw) * dzt(k) * rho(k,iw)
           scr1b(iw) = scr1b(iw) + dzt(k) * rho(k,iw)
        enddo
        scr1a(iw) = scr1a(iw) / scr1b(iw)
     enddo
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,q1zint_ll)
  endif

  if (allocated(q2zint_ll)) then
     do iw = 2, mwa
        scr1a(iw) = 0.
        scr1b(iw) = 0.
        do k = lpw(iw),mza
           scr1a(iw) = scr1a(iw) + addsc(2)%sclp(k,iw) * dzt(k) * rho(k,iw)
           scr1b(iw) = scr1b(iw) + dzt(k) * rho(k,iw)
        enddo
        scr1a(iw) = scr1a(iw) / scr1b(iw)
     enddo
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,q2zint_ll)
  endif

  if (allocated(qyzint_ll)) then
     do iw = 2, mwa
        scr1a(iw) = 0.
        scr1b(iw) = 0.
        do k = lpw(iw),mza
           scr1a(iw) = scr1a(iw) + (addsc(1)%sclp(k,iw) + 2. * addsc(2)%sclp(k,iw)) * dzt(k) * rho(k,iw)
           scr1b(iw) = scr1b(iw) + dzt(k) * rho(k,iw)
        enddo
        scr1a(iw) = scr1a(iw) / scr1b(iw)
     enddo
     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,qyzint_ll)
  endif

  !--------------------------------------------------------------------------------
  ! Interpolate 3D atmospheric time-accumulated fields
  !--------------------------------------------------------------------------------

  if (allocated(u_accum_ll) .or. allocated(v_accum_ll)) then
     do iw = 2, mwa
        npoly = itab_w(iw)%npoly

        do k = lpw(iw), mza

           vx = 0.0
           vy = 0.0
           vz = 0.0

           do jv = 1, npoly
              iv = itab_w(iw)%iv(jv)

              vx = vx + itab_w(iw)%ecvec_vx(jv) * real(vc_accum(k,iv))
              vy = vy + itab_w(iw)%ecvec_vy(jv) * real(vc_accum(k,iv))
              vz = vz + itab_w(iw)%ecvec_vz(jv) * real(vc_accum(k,iv))
           enddo

           scr2a(k,iw) = vx * vxn_ew(iw) + vy * vyn_ew(iw)
           scr2b(k,iw) = vx * vxn_ns(iw) + vy * vyn_ns(iw) + vz * vzn_ns(iw)

        enddo
     enddo

     if (allocated(u_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,u_accum_ll)
     if (allocated(v_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,v_accum_ll)
  endif

  if (allocated(w_accum_ll) .or. allocated(p_accum_ll)) then
     do iw = 2,mwa
        kb = lpw(iw)

        do k = kb,mza
           scr2a(k,iw) = 0.5 * real(wc_accum(k,iw) + wc_accum(k-1,iw))
           scr2b(k,iw) = real(press_accum(k,iw))
        enddo
     enddo

     if (allocated(w_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,w_accum_ll)
     if (allocated(p_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,p_accum_ll)
  endif

  if (allocated(t_accum_ll) .or. allocated(rv_accum_ll)) then
     do iw = 2,mwa
        kb = lpw(iw)

        do k = kb,mza
           scr2a(k,iw) = real(tair_accum(k,iw))
           scr2b(k,iw) = real(rr_v_accum(k,iw))
        enddo
     enddo

     if (allocated(t_accum_ll))  call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2a,t_accum_ll)
     if (allocated(rv_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,mza,mza-1,scr2b,rv_accum_ll)
  endif

  !------------------------------------------------------------------------
  ! Sum precipitation rates and accumulations over all microphysics species
  !------------------------------------------------------------------------

  scr1a(:) = 0.
  scr1b(:) = 0.
  scr1c(:) = 0.
  scr1d(:) = 0.
  scr1e(:) = 0.
  scr1f(:) = 0.
  scr1g(:) = 0.
  scr1h(:) = 0.
  scr1i(:) = 0.
  scr1j(:) = 0.

  if (allocated(pcprd))   scr1a(:) = scr1a(:) + real(pcprd(:))
  if (allocated(pcprr))   scr1a(:) = scr1a(:) + real(pcprr(:))
  if (allocated(pcprp))   scr1a(:) = scr1a(:) + real(pcprp(:))
  if (allocated(pcprs))   scr1a(:) = scr1a(:) + real(pcprs(:))
  if (allocated(pcpra))   scr1a(:) = scr1a(:) + real(pcpra(:))
  if (allocated(pcprg))   scr1a(:) = scr1a(:) + real(pcprg(:))
  if (allocated(pcprh))   scr1a(:) = scr1a(:) + real(pcprh(:))

  if (allocated(accpd))   scr1b(:) = scr1b(:) + real(accpd(:))
  if (allocated(accpr))   scr1b(:) = scr1b(:) + real(accpr(:))
  if (allocated(accpp))   scr1b(:) = scr1b(:) + real(accpp(:))
  if (allocated(accps))   scr1b(:) = scr1b(:) + real(accps(:))
  if (allocated(accpa))   scr1b(:) = scr1b(:) + real(accpa(:))
  if (allocated(accpg))   scr1b(:) = scr1b(:) + real(accpg(:))
  if (allocated(accph))   scr1b(:) = scr1b(:) + real(accph(:))

  if (allocated(aconpr))             scr1c(:) = real(            aconpr(:))
  if (allocated(rshort_accum))       scr1d(:) = real(      rshort_accum(:))
  if (allocated(rshortup_accum))     scr1e(:) = real(    rshortup_accum(:))
  if (allocated(rlong_accum))        scr1f(:) = real(       rlong_accum(:))
  if (allocated(rlongup_accum))      scr1g(:) = real(     rlongup_accum(:))
  if (allocated(rshort_top_accum))   scr1h(:) = real(  rshort_top_accum(:))
  if (allocated(rshortup_top_accum)) scr1i(:) = real(rshortup_top_accum(:))
  if (allocated(rlongup_top_accum))  scr1j(:) = real( rlongup_top_accum(:))

  ! For DCMIP 2016, convert precip rate from mm/s to m/s and accum precip from mm to m

  !DCMIP scr1a(:) = scr1a(:) * 1.e-3 ! converting precip rate to m/s for DCMIP 2016
  !DCMIP scr1b(:) = scr1b(:) * 1.e-3 ! converting accum precip to m for DCMIP 2016

  !-----------------------------------------------------------------------
  ! Interpolate accumulated precipitation, turbulent, and radiation fluxes
  !-----------------------------------------------------------------------

  if (allocated(           pcprmic_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,           pcprmic_ll)
  if (allocated(           accpmic_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,           accpmic_ll)
  if (allocated(           accpcon_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1c,           accpcon_ll)
  if (allocated(      rshort_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1d,      rshort_accum_ll)
  if (allocated(    rshortup_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1e,    rshortup_accum_ll)
  if (allocated(       rlong_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1f,       rlong_accum_ll)
  if (allocated(     rlongup_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1g,     rlongup_accum_ll)
  if (allocated(  rshort_top_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1h,  rshort_top_accum_ll)
  if (allocated(rshortup_top_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1i,rshortup_top_accum_ll)
  if (allocated( rlongup_top_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1j, rlongup_top_accum_ll)

  !---------------------------------------------------------------------------
  ! If surface model is activated, average accumulated surface-grid quantities
  ! over atmosphere grid cells and interpolate the averages
  !---------------------------------------------------------------------------

  if (isfcl == 1) then

     do n = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(n)
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
        scr1m(iw) = 0.
        scr1n(iw) = 0.
        scr1o(iw) = 0.
        scr1p(iw) = 0.
        scr1q(iw) = 0.
        scr1r(iw) = 0.
        scr1s(iw) = 0.
        scr1t(iw) = 0.
        scr1u(iw) = 0.
        scr1v(iw) = 0.

        ! Averages over SFC grid cells

        do j = 1,itab_w(iw)%jsfc2
           iwsfc = itab_w(iw)%iwsfc(j)
           jasfc = itab_w(iw)%jasfc(j)
           kwatm = itab_wsfc(iwsfc)%kwatm(jasfc)

           scr1a(iw) = scr1a(iw) + real(    vels_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1b(iw) = scr1b(iw) + real( airtemp_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1c(iw) = scr1c(iw) + real(  airrrv_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1d(iw) = scr1d(iw) + real( cantemp_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1e(iw) = scr1e(iw) + real(  canrrv_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1f(iw) = scr1f(iw) + real(skintemp_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)
           scr1g(iw) = scr1g(iw) + real(  sfluxt_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc) * cp
           scr1h(iw) = scr1h(iw) + real(  sfluxr_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc) * alvl
           scr1i(iw) = scr1i(iw) + real(  sfluxr_accum(iwsfc)) * itab_wsfc(iwsfc)%arcoariw(jasfc)

           ! Evaluate surface-layer quantities at 10 m and 2 m

           do k = 0,1
              zobs = real(10 - 8 * k)

              press_zobs = sfcg%prss(iwsfc) - zobs * sfcg%rhos(iwsfc) ! hydrostatic eqn.
              exner_zobs = (press_zobs * p00i) ** rocp

              canexner = (sfcg%prss(iwsfc) * p00i) ** rocp
              cantheta  = sfcg%cantemp(iwsfc) / canexner
              canthetav = cantheta         * (1.0 + eps_virt * sfcg%canrrv(iwsfc))
              airthetav = sfcg%airtheta(iwsfc) * (1.0 + eps_virt * sfcg%airrrv(iwsfc))

              tstar = -sfcg%sfluxt(iwsfc) / (sfcg%ustar(iwsfc) * sfcg%rhos(iwsfc))
              rstar = -sfcg%sfluxr(iwsfc) / (sfcg%ustar(iwsfc) * sfcg%rhos(iwsfc))

              ufree = (grav * sfcg%dzt_bot(iwsfc) * max(sfcg%wthv(iwsfc),0.0) / airthetav) ** onethird

              call sfclyr_profile (sfcg%vels(iwsfc), sfcg%ustar(iwsfc), tstar, rstar, &
                                   sfcg%dzt_bot(iwsfc), sfcg%rough(iwsfc), ufree, &
                                   cantheta, canthetav, sfcg%canrrv(iwsfc), airthetav, &
                                   zobs, wind_zobs, theta_zobs, rrv_zobs)

              frac = wind_zobs / max(1.e-3,sfcg%vels(iwsfc))

              uwind = ue(kwatm,iw) * frac
              vwind = ve(kwatm,iw) * frac

              if (k == 0) then
                 scr1s(iw) = scr1s(iw) + uwind                   * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1t(iw) = scr1t(iw) + vwind                   * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1o(iw) = scr1o(iw) + wind_zobs               * itab_wsfc(iwsfc)%arcoariw(jasfc)
              else
                 scr1u(iw) = scr1u(iw) + uwind                   * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1v(iw) = scr1v(iw) + vwind                   * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1p(iw) = scr1p(iw) + wind_zobs               * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1q(iw) = scr1q(iw) + theta_zobs * exner_zobs * itab_wsfc(iwsfc)%arcoariw(jasfc)
                 scr1r(iw) = scr1r(iw) + rrv_zobs                * itab_wsfc(iwsfc)%arcoariw(jasfc)
              endif
           enddo

        enddo

        ! Averages over LAND cells

        area_land_sum = 0.

        do j = itab_w(iw)%jland1,itab_w(iw)%jland2
           iwsfc = itab_w(iw)%iwsfc(j)
           iland = iwsfc - omland

           scr1j(iw) = scr1j(iw) + real(  wxferi_accum(iland)) * sfcg%area(iwsfc)
           scr1k(iw) = scr1k(iw) + real(  wxferp_accum(iland)) * sfcg%area(iwsfc)
           scr1l(iw) = scr1l(iw) + real(  wxfer1_accum(iland)) * sfcg%area(iwsfc)
           scr1m(iw) = scr1m(iw) &
                     + sum(land%sfcwater_mass(:,iland))        * sfcg%area(iwsfc)
           scr1n(iw) = scr1n(iw) &
                     + sum(land%soil_water(:,iland) * dslz(:)) * sfcg%area(iwsfc)

           area_land_sum = area_land_sum + sfcg%area(iwsfc)
        enddo

        if (area_land_sum > 1.e0) then
           scr1j(iw) = scr1j(iw) / area_land_sum
           scr1k(iw) = scr1k(iw) / area_land_sum
           scr1l(iw) = scr1l(iw) / area_land_sum
           scr1m(iw) = scr1m(iw) / area_land_sum
           scr1n(iw) = scr1n(iw) / area_land_sum
        endif

     enddo

     if (allocated(asfc_u10m_ll))            call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1s,asfc_u10m_ll)
     if (allocated(asfc_v10m_ll))            call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1t,asfc_v10m_ll)
     if (allocated(asfc_speed10m_ll))        call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1o,asfc_speed10m_ll)
     if (allocated(asfc_u2m_ll))             call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1u,asfc_u2m_ll)
     if (allocated(asfc_v2m_ll))             call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1v,asfc_v2m_ll)
     if (allocated(asfc_speed2m_ll))         call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1p,asfc_speed2m_ll)
     if (allocated(asfc_temp2m_ll))          call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1q,asfc_temp2m_ll)
     if (allocated(asfc_rvap2m_ll))          call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1r,asfc_rvap2m_ll)
     if (allocated(asfc_vels_accum_ll))      call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1a,asfc_vels_accum_ll)
     if (allocated(asfc_airtempk_accum_ll))  call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1b,asfc_airtempk_accum_ll)
     if (allocated(asfc_airrv_accum_ll))     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1c,asfc_airrv_accum_ll)
     if (allocated(asfc_cantempk_accum_ll))  call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1d,asfc_cantempk_accum_ll)
     if (allocated(asfc_canrv_accum_ll))     call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1e,asfc_canrv_accum_ll)
     if (allocated(asfc_skintempk_accum_ll)) call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1f,asfc_skintempk_accum_ll)
     if (allocated(asfc_sensflux_accum_ll))  call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1g,asfc_sensflux_accum_ll)
     if (allocated(asfc_latflux_accum_ll))   call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1h,asfc_latflux_accum_ll)
     if (allocated(asfc_vapflux_accum_ll))   call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1i,asfc_vapflux_accum_ll)
     if (allocated(al_wxferi_accum_ll))      call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1j,al_wxferi_accum_ll)
     if (allocated(al_wxferp_accum_ll))      call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1k,al_wxferp_accum_ll)
     if (allocated(al_wxfer1_accum_ll))      call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1l,al_wxfer1_accum_ll)
     if (allocated(al_sfcwater_tot_ll))      call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1m,al_sfcwater_tot_ll)
     if (allocated(al_soil_water_tot_ll))    call interp_htw_ll(npts,iws_loc,wts_loc,1,1,scr1n,al_soil_water_tot_ll)

  endif

  ! Compute difference fields

  timedifi = 1.0
  if (abs(time8 - time8_prev) > .99) then
     timedifi = 1.0 / real(time8 - time8_prev)
  endif

  if (allocated(pcpmic_dif2_ll)) &
     pcpmic_dif2_ll(:) =                  (accpmic_ll(:) -       pcpmic_dif2_ll(:)) * timedifi
  if (allocated(pcpcon_dif2_ll)) &
     pcpcon_dif2_ll(:) =                  (accpcon_ll(:) -       pcpcon_dif2_ll(:)) * timedifi

  !--------------------------------------------------------------------------------------
  !DCMIP if (allocated(pcpmic_dif2_ll)) &
  !DCMIP    pcpmic_dif2_ll(:) = pcpmic_dif2_ll(:) * 1.e-3 ! converting to m/s for DCMIP 2016
  !DCMIP if (allocated(pcpcon_dif2_ll)) &
  !DCMIP    pcpcon_dif2_ll(:) = pcpcon_dif2_ll(:) * 1.e-3 ! converting to m/s for DCMIP 2016
  !--------------------------------------------------------------------------------------
 
  if (allocated(pcpboth_dif2_ll)) &
     pcpboth_dif2_ll(:) =                pcpmic_dif2_ll(:) +        pcpcon_dif2_ll(:)
  if (allocated(rshort_dif2_ll)) &
     rshort_dif2_ll(:) =               (rshort_accum_ll(:) -        rshort_dif2_ll(:)) * timedifi
  if (allocated(rshortup_dif2_ll)) &
     rshortup_dif2_ll(:) =           (rshortup_accum_ll(:) -      rshortup_dif2_ll(:)) * timedifi
  if (allocated(rlong_dif2_ll)) &
     rlong_dif2_ll(:) =                 (rlong_accum_ll(:) -         rlong_dif2_ll(:)) * timedifi
  if (allocated(rlongup_dif2_ll)) &
     rlongup_dif2_ll(:) =             (rlongup_accum_ll(:) -       rlongup_dif2_ll(:)) * timedifi
  if (allocated(rshort_top_dif2_ll)) &
     rshort_top_dif2_ll(:) =       (rshort_top_accum_ll(:) -    rshort_top_dif2_ll(:)) * timedifi
  if (allocated(rshortup_top_dif2_ll)) &
     rshortup_top_dif2_ll(:) =   (rshortup_top_accum_ll(:) -  rshortup_top_dif2_ll(:)) * timedifi
  if (allocated(rlongup_top_dif2_ll)) &
     rlongup_top_dif2_ll(:) =     (rlongup_top_accum_ll(:) -   rlongup_top_dif2_ll(:)) * timedifi
  if (allocated(asfc_sensflux_dif2_ll)) &
     asfc_sensflux_dif2_ll(:) = (asfc_sensflux_accum_ll(:) - asfc_sensflux_dif2_ll(:)) * timedifi
  if (allocated(asfc_latflux_dif2_ll)) &
     asfc_latflux_dif2_ll(:) =   (asfc_latflux_accum_ll(:) -  asfc_latflux_dif2_ll(:)) * timedifi

  ! HDF5 write

  if (nl%ioutput_latlon == 1) then

     if (myrank == 0) write(io6,'(/,a)') "Writing lat/lon fields to file..."

!!  model.experiment_id.horizontal_resolution.levels.grid.equation.description.variable.nc

!!  olam.163.r25.L30.voronoi.nonhydro.variable_220_25km.PS.nc

     call makefnam(hnamel, hfilepref, current_time, 'LL', '$', 'h5')
     call shdf5_open(hnamel,'W',iclobber,trypario=.true.)

     ! Write any global attributes to file

     write(ofrq,'(I0,A1)') nint(nl%frqlatlon), 's'

     call shdf5_write_global_attribute("Conventions",    cvalue = "CF-1.0")
     call shdf5_write_global_attribute("model_id",       cvalue = "OLAM")
     call shdf5_write_global_attribute("time_frequency", cvalue = trim(ofrq))
     call shdf5_write_global_attribute("grid",           cvalue = "hexagonal")

     !-------------------------------------------------------------------------------------
     ! Write more global attributes if specifically for DCMIP 2016

     if (nl%test_case >= 110 .and. nl%test_case <= 131) then

        write (levs, '(a1,i2)') 'L',mza-1

        call shdf5_write_global_attribute("project_id",            cvalue = "DCMIP2016")
        call shdf5_write_global_attribute("institute_id",          cvalue = "University of Miami")
        call shdf5_write_global_attribute("modeling_realm",        cvalue = "atmos")

        if     (nl%test_case == 110) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_1-dry") 
           call shdf5_write_global_attribute("frequency",          cvalue = "1day")
           call shdf5_write_global_attribute("mdt",                cvalue = "900s")
        elseif (nl%test_case == 111) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_1-preciponly") 
           call shdf5_write_global_attribute("frequency",          cvalue = "1day")
           call shdf5_write_global_attribute("mdt",                cvalue = "900s")
        elseif (nl%test_case == 112) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_1-rjpbl") 
           call shdf5_write_global_attribute("frequency",          cvalue = "1day") 
           call shdf5_write_global_attribute("mdt",                cvalue = "900s")
        elseif (nl%test_case == 121) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_2-rjpbl") 
           call shdf5_write_global_attribute("frequency",          cvalue = "6hr") 
        elseif (nl%test_case == 122) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_2-bryanpbl") 
           call shdf5_write_global_attribute("frequency",          cvalue = "6hr") 
        elseif (nl%test_case == 131) then
           call shdf5_write_global_attribute("experiment_id",      cvalue = "2016_3") 
           call shdf5_write_global_attribute("frequency",          cvalue = "300s") 
        endif

        if     (nxp < 40 .and. ngrids == 1) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r200")
           call shdf5_write_global_attribute("description",           cvalue = "")
        elseif (nxp < 80 .and. ngrids == 1) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r100")
           call shdf5_write_global_attribute("description",           cvalue = "")
        elseif (nxp < 160 .and. ngrids == 1) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r50")
           call shdf5_write_global_attribute("description",           cvalue = "")
        elseif (nxp < 40 .and. ngrids == 2) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r100")
           call shdf5_write_global_attribute("description",           cvalue = "variable_200_100km")
        elseif (nxp < 80 .and. ngrids == 2) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r50")
           call shdf5_write_global_attribute("description",           cvalue = "variable_100_50km")
        elseif (nxp < 160 .and. ngrids == 2) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r25")
           call shdf5_write_global_attribute("description",           cvalue = "variable_50_25km")
        elseif (nxp < 40 .and. ngrids == 3) then
          call shdf5_write_global_attribute("horizontal resolution", cvalue = "r50")
           call shdf5_write_global_attribute("description",           cvalue = "variable_200_50km")
        elseif (nxp < 80 .and. ngrids == 3) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r25")
           call shdf5_write_global_attribute("description",           cvalue = "variable_100_25km")
        elseif (nxp < 160 .and. ngrids == 3) then
           call shdf5_write_global_attribute("horizontal resolution", cvalue = "r12")
           call shdf5_write_global_attribute("description",           cvalue = "variable_50_12km")
        endif

     endif
     !-------------------------------------------------------------------------------------

     ! Write coordinate variables to file (lat, lon, height)

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
                     long_name = "height above mean sea level",   &
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

     ! Now write lat/lon interpolated variables to file.
     ! THESE WRITES NEED THE ROUTINE SHDF5_OREC_LL TO WORK IN PARALLEL!!

     ndims    = 2
     idims(1) = nlon
     idims(2) = nlat
     dimnames(1) = 'lon'
     dimnames(2) = 'lat'

     ! Surface Quantities

     if (allocated(topo_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'TOPO_LL', rvar1=topo_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                     &
                           long_name = "topography height",                         &
                           standard_name = "surface_altitude",                      &
                           units = "m",                                             &
                           rmissing = rmissing,                                     &
                           stagpt = "LL"                                            )

     if (allocated(u_lpw_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'U_LPW_LL', rvar1=u_lpw_ll, gpoints=lls_loc,       &
                           dimnames = dimnames,                                             &
                           long_name = "eastward wind at lowest model layer above surface", &
                           standard_name = "eastward_wind",                                 &
                           units = "m s-1",                                                 &
                           rmissing = rmissing,                                             &
                           stagpt = "LL"                                                    )

     if (allocated(v_lpw_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'V_LPW_LL', rvar1=v_lpw_ll, gpoints=lls_loc,        &
                           dimnames = dimnames,                                              &
                           long_name = "northward wind at lowest model layer above surface", &
                           standard_name = "northward_wind",                                 &
                           units = "m s-1",                                                  &
                           rmissing = rmissing,                                              &
                           stagpt = "LL"                                                     )

     if (allocated(t_lpw_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'T_LPW_LL', rvar1=t_lpw_ll, gpoints=lls_loc,         &
                           dimnames = dimnames,                                               &
                           long_name = "air temperature at lowest model layer above surface", &
                           standard_name = "surface_air_temperature",                         &
                           units = "K",                                                       &
                           rmissing = rmissing,                                               &
                           stagpt = "LL"                                                      )

     if (allocated(rv_lpw_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RV_LPW_LL', rvar1=rv_lpw_ll, gpoints=lls_loc,                &
                           dimnames = dimnames,                                                        &
                           long_name = "water vapor mixing ratio at lowest model layer above surface", &
                           standard_name = "vapor_mixing_ratio",                                       &
                           units = "kg kg-1",                                                          &
                           rmissing = rmissing,                                                        &
                           stagpt = "LL"                                                               )

     if (allocated(slp_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'PS',     rvar1=slp_ll, gpoints=lls_loc,    &
        CALL shdf5_orec_ll(ndims, idims, 'SLP_LL', rvar1=slp_ll, gpoints=lls_loc,    &
                           dimnames = dimnames,                                      &
                           long_name = "surface pressure reduced to mean sea level", &
                           standard_name = "surface_air_pressure_at_sea_level",      &
                           units = "Pa",                                             &
                           rmissing = rmissing,                                      &
                           stagpt = "LL"                                             )

     if (allocated(pvap_lpw_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'PVAP_LPW_LL', rvar1=pvap_lpw_ll, gpoints=lls_loc,        &
                           dimnames = dimnames,                                                    &
                           long_name = "water vapor pressure at lowest model layer above surface", &
                           standard_name = "water_vapor_partial_pressure_in_air",                  &
                           units = "Pa",                                                           &
                           rmissing = rmissing,                                                    &
                           stagpt = "LL"                                                           )

     ! Precipitation rate at surface

     if (allocated(pcprmic_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'PRECL',      rvar1=pcprmic_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'PCPRMIC_LL', rvar1=pcprmic_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Resolved precipitation rate",                     &
                           standard_name = "large_scale_precipitation_rate",              &
!DCMIP                     units = "m s-1",                                               &
                           units = "kg m-2 s-1",                                          &
                           rmissing = rmissing,                                           &
                           stagpt = "LL"                                                  )

     if (allocated(q1zint_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'qCl',    rvar1=q1zint_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                      &
                           long_name = "Vertically integrated tracer qCl",           &
                           standard_name = "Added Scalar 1",                         &
                           units = "kg/kg",                                          &
                           rmissing=rmissing,                                        &
                           stagpt = "LL"                                             )

     if (allocated(q2zint_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'qCl2',  rvar1=q2zint_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                     &
                           long_name = "Vertically integrated tracer qCl2",         &
                           standard_name = "Added Scalar 2",                        &
                           units = "kg/kg",                                         &
                           rmissing=rmissing,                                       &
                           stagpt = "LL"                                            )

     if (allocated(qyzint_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'qCly',  rvar1=qyzint_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                     &
                           long_name = "Vertically integrated tracer sum qCly",     &
                           standard_name = "Added Scalar 1&2 sum",                  &
                           units = "kg/kg",                                         &
                           rmissing=rmissing,                                       &
                           stagpt = "LL"                                            )

     ! Accumulated precipitation flux at surface

     if (allocated(accpmic_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'PRECT',   rvar1=accpmic_ll, gpoints=lls_loc,    &
        CALL shdf5_orec_ll(ndims, idims, 'ACCPMIC_LL', rvar1=accpmic_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated resolved precipitation",              &
                           standard_name = "large_scale_precipitation_amount",            &
!DCMIP                     units = "m",                                                   &
                           units = "kg m-2",                                              &
                           rmissing = rmissing,                                           &
                           stagpt = "LL"                                                  )

     if (allocated(accpcon_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ACCPCON_LL', rvar1=accpcon_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated convective precipitation",            &
                           standard_name = "convective_precipitation_amount",             &
                           units = "kg m-2",                                              &
                           rmissing = rmissing,                                           &
                           stagpt = "LL"                                                  )

     ! Accumulated longwave and shortwave radiative fluxes at surface

     if (allocated(rshort_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_ACCUM_LL', rvar1=rshort_accum_ll, gpoints=lls_loc,   &
                           dimnames = dimnames,                                                       &
                           long_name = "Accumulated downwelling shortwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_shortwave_flux_wrt_time", &
                           units = "J m-2",                                                           &
                           rmissing = rmissing,                                                       &
                           stagpt = "LL"                                                              )

     if (allocated(rshortup_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_ACCUM_LL', rvar1=rshortup_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                         &
                           long_name = "Accumulated upwelling shortwave flux at surface",               &
                           standard_name = "integral_of_surface_upwelling_shortwave_flux_wrt_time",     &
                           units = "J m-2",                                                             &
                           rmissing = rmissing,                                                         &
                           stagpt = "LL"                                                                )

     if (allocated(rlong_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RLONG_ACCUM_LL', rvar1=rlong_accum_ll, gpoints=lls_loc,    &
                           dimnames = dimnames,                                                      &
                           long_name = "Accumulated downwelling longwave flux at surface",           &
                           standard_name = "integral_of_surface_downwelling_longwave_flux_wrt_time", &
                           units = "J m-2",                                                          &
                           rmissing = rmissing,                                                      &
                           stagpt = "LL"                                                             )

     if (allocated(rlongup_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_ACCUM_LL', rvar1=rlongup_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                       &
                           long_name = "Accumulated upwelling longwave flux at surface",              &
                           standard_name = "integral_of_surface_upwelling_longwave_flux_wrt_time",    &
                           units = "J m-2",                                                           &
                           rmissing = rmissing,                                                       &
                           stagpt = "LL"                                                              )

     ! Accumulated longwave and shortwave radiative fluxes at top-of-atmosphere
     ! Note: at TOA "incoming" and "outgoing" are used in place of "downwelling" and "upwelling"

     if (allocated(rshort_top_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RSHORT_TOP_ACCUM_LL', rvar1=rshort_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                             &
                           long_name = "Accumulated incoming shortwave flux at TOA",                        &
                           standard_name = "integral_of_toa_incoming_shortwave_flux_wrt_time",              &
                           units = "J m-2",                                                                 &
                           rmissing = rmissing,                                                             &
                           stagpt = "LL"                                                                    )

     if (allocated(rshortup_top_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RSHORTUP_TOP_ACCUM_LL', rvar1=rshortup_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                 &
                           long_name = "Accumulated outgoing shortwave flux at TOA",                            &
                           standard_name = "integral_of_toa_outgoing_shortwave_flux_wrt_time",                  &
                           units = "J m-2",                                                                     &
                           rmissing = rmissing,                                                                 &
                           stagpt = "LL"                                                                        )

     if (allocated(rlongup_top_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RLONGUP_TOP_ACCUM_LL' , rvar1=rlongup_top_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "Accumulated outgoing longwave flux at TOA",                            &
                           standard_name = "integral_of_toa_outgoing_longwave_flux_wrt_time",                  &
                           units = "J m-2",                                                                    &
                           rmissing = rmissing,                                                                &
                           stagpt = "LL"                                                                       )

     ! 'ASFC' quantities are area-weighted averages of sfcgrid cell quantities over a single atmosphere column

     if (allocated(asfc_u10m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_U10M_LL' , rvar1=asfc_u10m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 10m u wind_component",                             &
                           standard_name = "asfc_average_of_10m_u_wind_component",                         &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_v10m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_V10M_LL' , rvar1=asfc_v10m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 10m v wind component",                             &
                           standard_name = "asfc_average_of_10m_v_wind_component",                         &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_speed10m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_SPEED10M_LL' , rvar1=asfc_speed10m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 10m wind speed",                                   &
                           standard_name = "asfc_average_of_10m_wind_speed",                               &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_u2m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_U2M_LL' , rvar1=asfc_u2m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 2m wind u component",                                    &
                           standard_name = "asfc_average_of_2m_wind_u_component",                                &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_v2m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_V2M_LL' , rvar1=asfc_v2m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 2m wind v component",                                    &
                           standard_name = "asfc_average_of_2m_wind_v_component",                                &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_speed2m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_SPEED2M_LL' , rvar1=asfc_speed2m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 2m wind speed",                                    &
                           standard_name = "asfc_average_of_2m_wind_speed",                                &
                           units = "m s-1",                                                                &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_temp2m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_TEMP2M_LL' , rvar1=asfc_temp2m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 2m temperature",                                   &
                           standard_name = "asfc_average_of_2m_temperature",                               &
                           units = "K",                                                                    &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_rvap2m_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_RVAP2M_LL' , rvar1=asfc_rvap2m_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of 2m vapor mixing ratio",                            &
                           standard_name = "asfc_average_of_2m_vapor_mixing_ratio",                        &
                           units = "kg kg-1",                                                              &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_vels_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_VELS_ACCUM_LL' , rvar1=asfc_vels_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "ASFC average of accumulated atmosphere surface wind speed",        &
                           standard_name = "asfc_average_of_integral_of_atm_wind_speed_wrt_time",          &
                           units = "m",                                                                    &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(asfc_airtempk_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_AIRTEMPK_ACCUM_LL' , rvar1=asfc_airtempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                    &
                           long_name = "ASFC average of accumulated atmosphere temperature",                       &
                           standard_name = "asfc_average_of_integral_of_atm_temperature_wrt_time",                 &
                           units = "K s",                                                                          &
                           rmissing = rmissing,                                                                    &
                           stagpt = "LL"                                                                           )

     if (allocated(asfc_airrv_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_AIRRV_ACCUM_LL' , rvar1=asfc_airrv_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                              &
                           long_name = "ASFC average of accumulated atmosphere vapor mixing ratio",          &
                           standard_name = "asfc_average_of_integral_of_atm_vapor_mixing_ratio_wrt_time",    &
                           units = "kg s kg-1",                                                              &
                           rmissing = rmissing,                                                              &
                           stagpt = "LL"                                                                     )

     if (allocated(asfc_cantempk_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_CANTEMPK_ACCUM_LL' , rvar1=asfc_cantempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                    &
                           long_name = "ASFC average of accumulated canopy air temperature",                       &
                           standard_name = "asfc_average_of_integral_of_canopy_air_temperature_wrt_time",          &
                           units = "K s",                                                                          &
                           rmissing = rmissing,                                                                    &
                           stagpt = "LL"                                                                           )

     if (allocated(asfc_canrv_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_CANRV_ACCUM_LL' , rvar1=asfc_canrv_accum_ll, gpoints=lls_loc,     &
                           dimnames = dimnames,                                                                  &
                           long_name = "ASFC average of accumulated canopy air vapor mixing ratio",              &
                           standard_name = "asfc_average_of_integral_of_canopy_air_vapor_mixing_ratio_wrt_time", &
                           units = "kg s kg-1",                                                                  &
                           rmissing = rmissing,                                                                  &
                           stagpt = "LL"                                                                         )

     if (allocated(asfc_skintempk_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_SKINTEMPK_ACCUM_LL' , rvar1=asfc_skintempk_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                      &
                           long_name = "ASFC average of accumulated skin temperature",                               &
                           standard_name = "asfc_average_of_integral_of_skin_temperature_wrt_time",                  &
                           units = "J m-2",                                                                          &
                           rmissing = rmissing,                                                                      &
                           stagpt = "LL"                                                                             )

     if (allocated(asfc_sensflux_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_SENSFLUX_ACCUM_LL' , rvar1=asfc_sensflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                    &
                           long_name = "ASFC average of accumulated sensible heat flux",                           &
                           standard_name = "asfc_average_of_integral_of_sensible_heat_flux_wrt_time",              &
                           units = "J m-2",                                                                        &
                           rmissing = rmissing,                                                                    &
                           stagpt = "LL"                                                                           )

     if (allocated(asfc_latflux_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_LATFLUX_ACCUM_LL' , rvar1=asfc_latflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                  &
                           long_name = "ASFC average of accumulated latent heat flux",                           &
                           standard_name = "asfc_average_of_integral_of_latent_heat_flux_wrt_time",              &
                           units = "J m-2",                                                                      &
                           rmissing = rmissing,                                                                  &
                           stagpt = "LL"                                                                         )

     if (allocated(asfc_vapflux_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'ASFC_VAPFLUX_ACCUM_LL' , rvar1=asfc_vapflux_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                  &
                           long_name = "ASFC average of accumulated vapor flux",                                 &
                           standard_name = "asfc_average_of_integral_of_vapor_flux_wrt_time",                    &
                           units = "kg m-2",                                                                     &
                           rmissing = rmissing,                                                                  &
                           stagpt = "LL"                                                                         )

     ! 'AL' quantities are area-weighted averages of land cell quantities over a single atmosphere column

     if (allocated(al_wxferi_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'AL_WXFERI_ACCUM_LL' , rvar1=al_wxferi_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of accumulated infiltration flux",                      &
                           standard_name = "al_average_of_integral_of_infiltration_flux_wrt_time",         &
                           units = "m",                                                                    &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(al_wxferp_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'AL_WXFERP_ACCUM_LL' , rvar1=al_wxferp_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of accumulated percolation flux",                       &
                           standard_name = "al_average_of_integral_of_percolation_flux_wrt_time",          &
                           units = "m",                                                                    &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(al_wxfer1_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'AL_WXFER1_ACCUM_LL' , rvar1=al_wxfer1_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of accumulated soil bottom water flux",                 &
                           standard_name = "al_average_of_integral_of_soil_bottom_water_flux_wrt_time",    &
                           units = "m",                                                                    &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(al_sfcwater_tot_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'AL_SFCWATER_TOT_LL' , rvar1=al_sfcwater_tot_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                            &
                           long_name = "AL average of total surface water",                                &
                           standard_name = "al_average_of_total_surface_water",                            &
                           units = "kg m-2",                                                               &
                           rmissing = rmissing,                                                            &
                           stagpt = "LL"                                                                   )

     if (allocated(al_soil_water_tot_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'AL_SOIL_WATER_TOT_LL' , rvar1=al_soil_water_tot_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                                                &
                           long_name = "AL average of total soil water",                                       &
                           standard_name = "al_average_of_total_soil_water",                                   &
                           units = "m",                                                                        &
                           rmissing = rmissing,                                                                &
                           stagpt = "LL"                                                                       )

     ! 3D accumulated atmospheric fields

     ndims = 3
     idims(1) = nlon
     idims(2) = nlat
     idims(3) = mza-1
     dimnames(1) = 'lon'
     dimnames(2) = 'lat'
     dimnames(3) = 'z'

     if (allocated(u_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'U_ACCUM_LL', rvar2=u_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated eastward wind",                       &
                           standard_name = "integral_of_eastward_wind_wrt_time",          &
                           units = "m",                                                   &
                           rmissing=rmissing,                                             &
                           stagpt = "LL"                                                  )

     if (allocated(v_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'V_ACCUM_LL', rvar2=v_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated northward wind",                      &
                           standard_name = "integral_of_northward_wind_wrt_time",         &
                           units = "m",                                                   &
                           rmissing=rmissing,                                             &
                           stagpt = "LL"                                                  )

     if (allocated(w_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'W_ACCUM_LL', rvar2=w_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated upward wind",                         &
                           standard_name = "integral_of_upward_wind_wrt_time",            &
                           units = "m",                                                   &
                           rmissing=rmissing,                                             &
                           stagpt = "LL"                                                  )

     if (allocated(t_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'T_ACCUM_LL', rvar2=t_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated air temperature",                     &
                           standard_name = "integral_of_air_temperature_wrt_time",        &
                           units = "K s",                                                 &
                           rmissing=rmissing,                                             &
                           stagpt = "LL"                                                  )

     if (allocated(rv_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'RV_ACCUM_LL', rvar2=rv_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                             &
                           long_name = "Accumulated water vapor mixing ratio",              &
                           standard_name = "integral_of_vapor_mixing_ratio_wrt_time",       &
                           units = "kg s kg-1",                                             &
                           rmissing=rmissing,                                               &
                           stagpt = "LL"                                                    )

     if (allocated(p_accum_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'P_ACCUM_LL', rvar2=p_accum_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                           &
                           long_name = "Accumulated air pressure",                        &
                           standard_name = "integral_of_air_pressure_wrt_time",           &
                           units = "Pa s",                                                &
                           rmissing=rmissing,                                             &
                           stagpt = "LL"                                                  )

     ! 3D atmospheric fields

     ndims = 3
     idims(1) = nlon
     idims(2) = nlat
     idims(3) = mza-1
     dimnames(1) = 'lon'
     dimnames(2) = 'lat'
     dimnames(3) = 'z'

     if (allocated(u_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'U',    rvar2=u_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'U_LL', rvar2=u_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "eastward wind",                       &
                           standard_name = "eastward_wind",                   &
                           units = "m s-1",                                   &
                           rmissing=rmissing,                                 &
                           stagpt = "LL"                                      )

     if (allocated(v_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'V',    rvar2=v_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'V_LL', rvar2=v_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "northward wind",                      &
                           standard_name = "northward_wind",                  &
                           units = "m s-1",                                   &
                           rmissing=rmissing,                                 &
                           stagpt = "LL"                                      )

     if (allocated(w_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'W',    rvar2=w_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'W_LL', rvar2=w_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "upward wind",                         &
                           standard_name = "upward_wind",                     &
                           units = "m s-1",                                   &
                           rmissing=rmissing,                                 &
                           stagpt = "LL"                                      )

     if (allocated(th_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'THETA',    rvar2=th_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'THETA_LL', rvar2=th_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                       &
                           long_name = "air potential temperature",                   &
                           standard_name = "air_potential_temperature",               &
                           units = "K",                                               &
                           rmissing=rmissing,                                         &
                           stagpt = "LL"                                              )

     if (allocated(t_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'T_LL', rvar2=t_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "air temperature",                     &
                           standard_name = "air_temperature",                 &
                           units = "K",                                       &
                           rmissing=rmissing,                                 &
                           stagpt = "LL"                                      )

     if (allocated(rv_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'Qv',    rvar2=rv_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'RV_LL', rvar2=rv_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "water vapor mixing ratio",              &
                           standard_name = "vapor_mixing_ratio",                &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(rc_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'Qc',    rvar2=rc_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'RC_LL', rvar2=rc_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "cloud water mixing ratio ratio",        &
                           standard_name = "cloud_water_mixing_ratio",          &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(rr_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'Qr',    rvar2=rr_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'RR_LL', rvar2=rr_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "rain mixing ratio ratio",               &
                           standard_name = "rain_mixing_ratio",                 &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(p_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'P',    rvar2=p_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'P_LL', rvar2=p_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                               &
                           long_name = "air pressure",                        &
                           standard_name = "air_pressure",                    &
                           units = "Pa",                                      &
                           rmissing=rmissing,                                 &
                           stagpt = "LL"                                      )

     if (allocated(q1_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'Q1',    rvar2=q1_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'Q1_LL', rvar2=q1_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "Added Scalar 1",                        &
                           standard_name = "Added Scalar 1",                    &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(q2_ll)) &
!DCMIP  CALL shdf5_orec_ll(ndims, idims, 'Q2',    rvar2=q2_ll, gpoints=lls_loc, &
        CALL shdf5_orec_ll(ndims, idims, 'Q2_LL', rvar2=q2_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "Added Scalar 2",                        &
                           standard_name = "Added Scalar 2",                    &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(q3_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'Q3_LL', rvar2=q3_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "Added Scalar 3",                        &
                           standard_name = "Added Scalar 3",                    &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(q4_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'Q4_LL', rvar2=q4_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "Added Scalar 4",                        &
                           standard_name = "Added Scalar 4",                    &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

     if (allocated(q5_ll)) &
        CALL shdf5_orec_ll(ndims, idims, 'Q5_LL', rvar2=q5_ll, gpoints=lls_loc, &
                           dimnames = dimnames,                                 &
                           long_name = "Added Scalar 5",                        &
                           standard_name = "Added Scalar 5",                    &
                           units = "kg kg-1",                                   &
                           rmissing=rmissing,                                   &
                           stagpt = "LL"                                        )

! 'DIF2_LL' fields are not written to hdf5 files; they are only plotted

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

  !------------------------------------------------------------
  ! Contour plot or tile plot 2D slices of 3D lat-lon fields
  !------------------------------------------------------------

  ! klev = vertical level to plot in 3D fields

  ! klev = min(23,mza)
  klev = 2

  if (allocated(u_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,u_ll,113, &
                      "Zonal wind", "(m/s)", rmissing)

  if (allocated(v_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,v_ll,113, &
                      "Meridional wind", "(m/s)", rmissing)

  if (allocated(w_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,w_ll,125, &
                      "Vertical wind", "(m/s)", rmissing)

  if (allocated(t_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,t_ll,4, &
                      "Temperature", "(K)", rmissing)

  if (allocated(th_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,th_ll,4, &
                      "Potential Temperature", "(K)", rmissing)

  if (allocated(rv_ll)) then
     ! Preserve missing values
     where( abs(rv_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = rv_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,5, &
                      "Water vapor mixing ratio", "(g/kg)", rmissing)
  endif

  if (allocated(rc_ll)) then
     ! Preserve missing values
     where( abs(rc_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = rc_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,5, &
                      "Cloud water mixing ratio", "(g/kg)", rmissing)
  endif

  if (allocated(rr_ll)) then
     ! Preserve missing values
     where( abs(rr_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = rr_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,5, &
                      "Rain mixing ratio", "(g/kg)", rmissing)
  endif

  if (allocated(p_ll)) then
     ! Preserve missing values
     where( abs(p_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = p_ll(:,:) * 0.01
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,17, &
                      "Air Pressure", "(mb)", rmissing)
  endif

  if (allocated(q1_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,q1_ll,5, &
                      "Added scalar 1", "(kg/kg)", rmissing)

  if (allocated(q2_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,q2_ll,5, &
                      "Added scalar 2", "(kg/kg)", rmissing)

  if (allocated(q3_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,q3_ll,5, &
                      "Added scalar 3", "(kg/kg)", rmissing)

  if (allocated(q4_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,q4_ll,5, &
                      "Added scalar 4", "(kg/kg)", rmissing)

  if (allocated(q5_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,q5_ll,5, &
                      "Added scalar 5", "(kg/kg)", rmissing)

  if (allocated(u_accum_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,u_accum_ll,304, &
                      "Accumulated zonal wind", "(m)", rmissing)

  if (allocated(v_accum_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,v_accum_ll,304, &
                      "Accumulated meridional wind", "(m)", rmissing)

  if (allocated(w_accum_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,w_accum_ll,304, &
                      "Accumulated vertical wind", "(m)", rmissing)

  if (allocated(t_accum_ll)) &
     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,t_accum_ll,204, &
                      "Accumulated temperature", "(K s)", rmissing)

  if (allocated(rv_accum_ll)) then
     ! Preserve missing values
     where( abs(rv_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = rv_accum_ll(:,:) * 1000.
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,204, &
                      "Accumulated water vapor mixing ratio", "(g s/kg)", rmissing)
  endif

  if (allocated(p_accum_ll)) then
     ! Preserve missing values
     where( abs(p_ll(:,:) - rmissing) > 1.e-7 )
        scr3_ll(:,:) = p_accum_ll(:,:) * 0.01
     elsewhere
        scr3_ll(:,:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,mza-1,klev,alon,alat,npts,lls_loc,scr3_ll,204, &
                      "Accumulated air pressure", "(hPa s)", rmissing)
  endif

  !------------------------------------------------------------
  ! Contour plot or tile plot 2D lat-lon fields
  !------------------------------------------------------------

  if (allocated(u_lpw_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,u_lpw_ll,113, &
                      "Zonal wind at lpw", "(m/s)", rmissing)

  if (allocated(v_lpw_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,v_lpw_ll,113, &
                      "Meridional wind at lpw", "(m/s)", rmissing)

  if (allocated(t_lpw_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,t_lpw_ll,4, &
                      "Air temperature at lpw", "(K)", rmissing)

  if (allocated(rv_lpw_ll)) then
     ! Preserve missing values
     where( abs(rv_lpw_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) = rv_lpw_ll(:) * 1000.
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,5, &
                      "Water vapor mixing ratio at lpw", "(g/kg)", rmissing)
  endif

  if (allocated(slp_ll)) then
     ! Preserve missing values
     where( abs(slp_ll(:) - rmissing) > 1.e-7 )
        scr2_ll(:) = slp_ll(:) * .01 
     elsewhere
        scr2_ll(:) = rmissing
     endwhere

     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,17, &
                      "Sea level pressure", "(mb)", rmissing)
  endif

  if (allocated(pvap_lpw_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pvap_lpw_ll,35, &
                      "Water vapor pressure at lpw", "(mb)", rmissing)

  if (allocated(pcprmic_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pcprmic_ll,204, &
                      "Microphysics precip rate", "(kg/m^2/s)", rmissing)

  if (allocated(topo_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,topo_ll,476, &
                      "Topography Height", "(m)", rmissing)

  if (allocated(q1zint_ll)) then
     scr2_ll(:) = q1zint_ll(:) / 4.e-6
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,170, &
                      "q1zint", " ", rmissing)
  endif

  if (allocated(q2zint_ll)) then
     scr2_ll(:) = q2zint_ll(:) / 4.e-6
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,175, &
                      "q2zint", " ", rmissing)
  endif

  if (allocated(qyzint_ll)) then
     scr2_ll(:) = qyzint_ll(:) / 4.e-6
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,scr2_ll,176, &
                      "qyzint", " ", rmissing)
  endif

  if (allocated(rshort_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_accum_ll,204, &
                      "Accum surface downward s/w radiation", "(J/m^2)", rmissing)

  if (allocated(rshortup_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_accum_ll,204, &
                      "Accum surface upward s/w radiation", "(J/m^2)", rmissing)

  if (allocated(rlong_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlong_accum_ll,204, &
                      "Accum surface downward l/w radiation", "(J/m^2)", rmissing)

  if (allocated(rlongup_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_accum_ll,204, &
                      "Accum surface upward l/w radiation", "(J/m^2)", rmissing)

  if (allocated(rshort_top_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshort_top_accum_ll,204, &
                      "Accum TOA downward s/w radiation", "(J/m^2)", rmissing)

  if (allocated(rshortup_top_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rshortup_top_accum_ll,204, &
                      "Accum TOA upward s/w radiation", "(J/m^2)", rmissing)

  if (allocated(rlongup_top_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,rlongup_top_accum_ll,204, &
                      "Accum TOA upward l/w radiation", "(J/m^2)", rmissing)

  if (allocated(accpmic_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpmic_ll,441, &
                      "Accum microphysics precip", "(kg/m^2)", rmissing)

  if (allocated(accpcon_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,accpcon_ll,441, &
                      "Accum convective precip", "(kg/m^2)", rmissing)

  if (allocated(asfc_u10m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_u10m_ll,79, &
                      "asfc 10m wind u component", "(m s-1)", rmissing)

  if (allocated(asfc_v10m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_v10m_ll,79, &
                      "asfc 10m wind v component", "(m s-1)", rmissing)

  if (allocated(asfc_speed10m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_speed10m_ll,79, &
                      "asfc 10m wind speed", "(m s-1)", rmissing)

  if (allocated(asfc_u2m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_u2m_ll,79, &
                      "asfc 2m wind u component", "(m s-1)", rmissing)

  if (allocated(asfc_v2m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_v2m_ll,79, &
                      "asfc 2m wind v component", "(m s-1)", rmissing)

  if (allocated(asfc_speed2m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_speed2m_ll,79, &
                      "asfc 2m wind speed", "(m s-1)", rmissing)

  if (allocated(asfc_temp2m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_temp2m_ll,4, &
                      "asfc 2m temperature", "(K)", rmissing)

  if (allocated(asfc_rvap2m_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_rvap2m_ll,84, &
                      "asfc 2m vapor mixing ratio", "(kg kg-1)", rmissing)

  if (allocated(asfc_vels_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_vels_accum_ll,113, &
                      "asfc accum wind speed", "(m)", rmissing)

  if (allocated(asfc_airtempk_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_airtempk_accum_ll,204, &
                      "asfc accum atmosphere temperature", "(K s)", rmissing)

  if (allocated(asfc_airrv_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_airrv_accum_ll,204, &
                      "asfc accum atmosphere vapor mixing ratio", "(kg s / kg)", rmissing)

  if (allocated(asfc_cantempk_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_cantempk_accum_ll,204, &
                      "asfc accum canopy air temperature", "(K s)", rmissing)

  if (allocated(asfc_canrv_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_canrv_accum_ll,204, &
                      "asfc accum canopy air vapor mixing ratio", "(kg s / kg)", rmissing)

  if (allocated(asfc_skintempk_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_skintempk_accum_ll,204, &
                      "asfc accum skin temperature", "(K s)", rmissing)

  if (allocated(asfc_sensflux_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_sensflux_accum_ll,304, &
                      "asfc accum sensible heat flux", "(J/m^2)", rmissing)

  if (allocated(asfc_latflux_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_latflux_accum_ll,304, &
                      "asfc accum latent heat flux", "(J/m^2)", rmissing)

  if (allocated(asfc_vapflux_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,asfc_vapflux_accum_ll,304, &
                      "asfc accum vapor flux", "(kg/m^2)", rmissing)

  if (allocated(al_wxferi_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_wxferi_accum_ll,304, &
                      "al accum infiltration flux", "(m)", rmissing)

  if (allocated(al_wxferp_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_wxferp_accum_ll,304, &
                      "al accum percolation flux", "(m)", rmissing)

  if (allocated(al_wxfer1_accum_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_wxfer1_accum_ll,304, &
                      "al accum soil bottom water flux", "(m)", rmissing)

  if (allocated(al_sfcwater_tot_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_sfcwater_tot_ll,200, &
                      "al total surface water", "(kg/m^2)", rmissing)

  if (allocated(al_soil_water_tot_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,al_soil_water_tot_ll,3, &
                      "al total soil water", "(m)", rmissing)

  if (allocated(pcpmic_dif2_ll)) &
     call tileplot_ll(nlon,nlat,1,1,alon,alat,npts,lls_loc,pcpmic_dif2_ll,204, &
                      "Time-interval accum microphysics precip", "(kg/m^2)", rmissing)

  !---------------------------------------------------------------------------
  ! Contour plot 1D arrays that are longitudinal averages of 2D lat-lon fields
  ! This will only work for serial runs that have the entire lat/lon memory
  !---------------------------------------------------------------------------

  if (.not. latplot) return

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

  arr_ll = reshape( pcpmic_dif2_ll, (/nlon,nlat/) )
  value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

  call oplot_xy2('1','N','a','N',aspect,scalelab,8,0,   &
                 nlat, alat, value,             &
                 'latitude','PCPDIF2 (mm/day)', &
                 alat(1),alat(nlat),alatinc,3,  &
                 ymin, ymax, yinc, 5            )

  arr_ll = reshape(pcpcon_dif2_ll, (/nlon,nlat/) )
  value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

  call oplot_xy2('1','N','a','N',aspect,scalelab,1,0,  &
                 nlat, alat, value,            &
                 'latitude',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  arr_ll = reshape(pcpboth_dif2_ll, (/nlon,nlat/) )
  value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

  call oplot_xy2('1','N','a','N',aspect,scalelab,11,0, &
                 nlat, alat, value,            &
                 'latitude',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  !---------------------------------------------------------------

  ymin = -10.
  ymax = 500.
  yinc = 10.

  if (allocated(rshort_dif2_ll)) then
     arr_ll = reshape(rshort_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N','a','N',aspect,scalelab,8,0,  &
                    nlat, alat, value,            &
                    'latitude','SFC RAD (W/m^2)', &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  if (allocated(rshortup_dif2_ll)) then
     arr_ll = reshape(rshortup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N','a','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  if (allocated(rlong_dif2_ll)) then
     arr_ll = reshape(rlong_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N','a','N',aspect,scalelab,1,0,  &
                    nlat, alat, value,            &
                    'latitude','',                &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  if (allocated(rlongup_dif2_ll)) then
     arr_ll = reshape(rlongup_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('2','N','a','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  ! Plot reference lines

  value(1:nlat) = 100.
  call oplot_xy2('2','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 200.
  call oplot_xy2('2','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 300.
  call oplot_xy2('2','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 400.
  call oplot_xy2('2','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',               &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  !---------------------------------------------------------------

  if (allocated(rshort_top_dif2_ll)) then
     arr_ll = reshape(rshort_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N','a','N',aspect,scalelab,8,0,   &
                    nlat, alat, value,             &
                    'latitude','TOA RAD (W/m^2) ', &
                    alat(1),alat(nlat),alatinc,3,  &
                    ymin, ymax, yinc, 5            )
  endif

  if (allocated(rshortup_top_dif2_ll)) then
     arr_ll = reshape(rshortup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N','a','N',aspect,scalelab,8,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  if (allocated(rlongup_top_dif2_ll)) then
     arr_ll = reshape(rlongup_top_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('3','N','a','N',aspect,scalelab,1,40, &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  ! Plot reference lines

  value(1:nlat) = 100.
  call oplot_xy2('3','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',                      &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 200.
  call oplot_xy2('3','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',                      &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 300.
  call oplot_xy2('3','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',                      &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  value(1:nlat) = 400.
  call oplot_xy2('3','N','a','N',aspect,scalelab,10,0, &
                 nlat, alat, value,            &
                 ' ',' ',                      &
                 alat(1),alat(nlat),alatinc,3, &
                 ymin, ymax, yinc, 5           )

  !---------------------------------------------------------------

  ymin = -50.
  ymax = 200.
  yinc = 10.

  if (allocated(asfc_sensflux_dif2_ll)) then
     arr_ll = reshape(asfc_sensflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N','a','N',aspect,scalelab,17,0,      &
                    nlat, alat, value,                 &
                    'latitude','ASFC S/L FLUX (W/m^2)', &
                    alat(1),alat(nlat),alatinc,3,      &
                    ymin, ymax, yinc, 5                )
  endif

  if (allocated(asfc_latflux_dif2_ll)) then
     arr_ll = reshape(asfc_latflux_dif2_ll, (/nlon,nlat/) )
     value(1:nlat) = sum(arr_ll(:,1:nlat),1) / real(nlon)

     call oplot_xy2('4','N','a','N',aspect,scalelab,9,0,  &
                    nlat, alat, value,            &
                    'latitude',' ',               &
                    alat(1),alat(nlat),alatinc,3, &
                    ymin, ymax, yinc, 5           )
  endif

  !---------------------------------------------------------------

  call o_frame()

  call o_clswk()

  if (allocated(       pcpmic_dif2_ll))        pcpmic_dif2_ll(:) =             accpmic_ll(:)
  if (allocated(       pcpcon_dif2_ll))        pcpcon_dif2_ll(:) =             accpcon_ll(:)
  if (allocated(       rshort_dif2_ll))        rshort_dif2_ll(:) =        rshort_accum_ll(:)
  if (allocated(     rshortup_dif2_ll))      rshortup_dif2_ll(:) =      rshortup_accum_ll(:)
  if (allocated(        rlong_dif2_ll))         rlong_dif2_ll(:) =         rlong_accum_ll(:)
  if (allocated(      rlongup_dif2_ll))       rlongup_dif2_ll(:) =       rlongup_accum_ll(:)
  if (allocated(   rshort_top_dif2_ll))    rshort_top_dif2_ll(:) =    rshort_top_accum_ll(:)
  if (allocated( rshortup_top_dif2_ll))  rshortup_top_dif2_ll(:) =  rshortup_top_accum_ll(:)
  if (allocated(  rlongup_top_dif2_ll))   rlongup_top_dif2_ll(:) =   rlongup_top_accum_ll(:)
  if (allocated(asfc_sensflux_dif2_ll)) asfc_sensflux_dif2_ll(:) = asfc_sensflux_accum_ll(:)
  if (allocated( asfc_latflux_dif2_ll))  asfc_latflux_dif2_ll(:) =  asfc_latflux_accum_ll(:)

  time8_prev = time8

end subroutine fields3_ll

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

  real,         parameter :: aspect   =  0.
  character(1), parameter :: proj     = 'L'
  character(1), parameter :: cbar     = 'c'
  character(1), parameter :: panel    = 'n'
  character(1), parameter :: frameoff = 'n'
  character(1), parameter :: pltborder= 'a'

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

  call oplot_panel(panel, frameoff, pltborder, cbar, aspect, proj)

  if (myrank == 0) then

     call o_gsplci(10)
     call o_gstxci(10)
     call o_sflush()

     bsize = .016 * (op%hp2 - op%hp1)
     call o_set (op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)
     call o_plchhq(0.5, op%fnamey, trim(trim(fldname)//' '//trim(units)), bsize, 0., 0.)

     call o_set(op%h1,op%h2,op%v1,op%v2,op%xmin,op%xmax,op%ymin,op%ymax,1)

     call oplot_xy2(panel, frameoff, pltborder, cbar, 0., scalelab, 10, 0,     &
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
