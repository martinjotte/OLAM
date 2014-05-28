subroutine fields_ll()

! This subroutine is a template, intended for user modification, for interpolating
! selected fields from the OLAM grid to a structured latitude-longitude grid
! of either limited or global area.  

! The following tasks are performed:

! (1) Allocate arrays for interpolated fields
! (2) Allocate scratch arrays
! (3) Copy each selected model field to scratch array, or compute it from
!     model fields if required
! (4) Call subroutine to interpolate each selected field to latitude-longitude array

use mem_ijtabs,  only: itab_w, itab_u, itab_v
use mem_basic,   only: uc, vc, wc, rho, press, theta, sh_w, sh_v, &
                       vxe, vye, vze

use mem_grid,    only: mza, mua, mva, mwa, lpu, lpv, lpw, &
                       xem, yem, zem, xeu, yeu, zeu, &
                       xev, yev, zev, xew, yew, zew, &
                       topw, glatm, glonm, glatw, glonw, zm, zt, &
                       unx, uny, unz, vnx, vny, vnz, dzt

use misc_coms,   only: io6, meshtype, current_time, hfilepref, iclobber

use consts_coms, only: p00, rocp, piu180, erad, eradi, pio180, cp, alvl, rvap, r8

use hdf5_utils,  only: shdf5_open, shdf5_orec, shdf5_close

use max_dims,    only: pathlen

! NOTE: The fields used from the MEM_MICRO, MEM_CUPARM, and MEM_FLUX_ACCUM
! modules are time-integrated fluxes of mass or energy and have units of
! [kg/m^2] or [J/m^2].  Their values written to any history file represent flux
! integrals from the beginning of a simulation to the time of the history file
! write.  Thus, time averages of a flux over any time interval can be computed
! as the difference between the flux integrals in two different history files
! divided by the difference in simulation times when the files were written.

use mem_micro,   only: accpd, accpr, accpp, accps, accpa, accpg, accph

use mem_cuparm,  only: aconpr

use mem_flux_accum, only: rshort_accum, rshortup_accum, rlong_accum, rlongup_accum, &
                          rshort_top_accum, rshortup_top_accum, rlongup_top_accum, &
                          sflux_t_accum, sflux_r_accum

! NOTE: The fields used from the MEM_TIMEAVG module are time-averaged surface
! or top-of-atmosphere fluxes.  OLAM currently evaluates these as averages over
! the time interval between consecutive history-file writes, and it therefore 
! resets the averages to zero after each history write.  Consequently, for these
! fields to have their correct meaning when interpolated to the latitude-longitude
! grid, the interpolation should be done only at times when history file writes 
! are performed.  Both the accumulation and resetting of averages are performed
! in subroutine accum_timeavg in the source code file mem_timeavg.f90.

use mem_timeavg, only: rshort_avg, rshortup_avg, rlong_avg, rlongup_avg, &
                       rshort_top_avg, rshortup_top_avg, rlongup_top_avg, &
                       sflux_t_avg, sflux_r_avg

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
!real   , parameter :: beglon = -92.   ! minimum (westernmost) longitude (deg)
!real   , parameter :: endlon = -72.   ! maximum (easternmost) longitude (deg)
!real   , parameter :: beglat =  23.   ! minimum (southernmost) latitude (deg)
!real   , parameter :: endlat =  35.   ! maximum (northernmost) latitude (deg)
!integer, parameter :: rf = 30         ! horiz resolution factor (pts per deg)

!real   , parameter :: beglon = -91.   ! minimum (westernmost) longitude (deg)
!real   , parameter :: endlon = -71.   ! maximum (easternmost) longitude (deg)
!real   , parameter :: beglat =  17.   ! minimum (southernmost) latitude (deg)
!real   , parameter :: endlat =  37.   ! maximum (northernmost) latitude (deg)
!integer, parameter :: rf = 30         ! horiz resolution factor (pts per deg)

! NLON and NLAT are the number of longitude and latitude points in the lat-lon
! arrays that get filled by interpolation from the native OLAM grid.  They may
! be computed from values of BEGLON, ENDLON, BEGLAT, ENDLAT, and RF, or they
! may be specified directly.

integer, parameter :: nlon = nint(endlon - beglon) * rf + 1  ! # of lon values
integer, parameter :: nlat = nint(endlat - beglat) * rf + 1  ! # of lat values

!integer, parameter :: nlon = 641   ! # of lon values
!integer, parameter :: nlat = 641   ! # of lat values

!--------------------------------------------------------------------------------
! THE FOLLOWING ARRAYS WILL CONTAIN THE FIELDS THAT ARE DEFINED ON THE 
! LATITUDE-LONGITUDE GRID.  EXCEPT FOR THE 1D VERTICAL ARRAYS, THEY ARE
! DIMENSIONED USING THE ABOVE NLON,NLAT PARAMETERS.  THE USER SHOULD ADD
! NEW ARRAYS AS REQUIRED.
!--------------------------------------------------------------------------------

!----------
! 1D FIELDS
!----------

real :: alat(nlat)  ! latitudes of grid points (deg)
real :: alon(nlon)  ! longitudes of grid points (deg)
real :: zlev(mza-2) ! heights above sea level of grid points (m)

!-------------------------------------
! 2D FIELDS (TA means 'time averaged')
!-------------------------------------

real :: topo_ll            (nlon,nlat) ! topography height (m)
real :: u_sfc_ll           (nlon,nlat) ! surface zonal wind component (m/s)
real :: v_sfc_ll           (nlon,nlat) ! surface meridional wind component (m/s)
real :: t_sfc_ll           (nlon,nlat) ! surface air temperature (K)
real :: r_sfc_ll           (nlon,nlat) ! surface water vapor mixing ratio (kg/kg)
real :: pvap_sfc_ll        (nlon,nlat) ! surface vapor pressure (Pa)
real :: slp_ll             (nlon,nlat) ! sea level pressure (Pa)

real :: accpmic_ll           (nlon,nlat) ! accum microphysics precip (kg/m^2)
real :: accpcon_ll           (nlon,nlat) ! accum convective precip (kg/m^2)
real :: vapflux_accum_ll     (nlon,nlat) ! accum sfc vapor (kg/m^2)
real :: sensflux_accum_ll    (nlon,nlat) ! accum sfc sensible heat (J/m^2)
real :: latflux_accum_ll     (nlon,nlat) ! accum sfc latent heat (J/m^2)
real :: rshort_accum_ll      (nlon,nlat) ! accum sfc downward s/w rad (J/m^2)
real :: rshortup_accum_ll    (nlon,nlat) ! accum sfc upward s/w rad (J/m^2)
real :: rlong_accum_ll       (nlon,nlat) ! accum sfc downward l/w rad (J/m^2)
real :: rlongup_accum_ll     (nlon,nlat) ! accum sfc upward l/w rad (J/m^2)
real :: rshort_top_accum_ll  (nlon,nlat) ! accum TOA downward s/w rad (J/m^2)
real :: rshortup_top_accum_ll(nlon,nlat) ! accum TOA upward s/w rad (J/m^2)
real :: rlongup_top_accum_ll (nlon,nlat) ! accum TOA upward l/w rad (J/m^2)

real :: vapflux_avg_ll     (nlon,nlat) ! avg sfc vapor flux (kg/m^2/s)
real :: sensflux_avg_ll    (nlon,nlat) ! avg sfc sensible heat flux (W/m^2)
real :: latflux_avg_ll     (nlon,nlat) ! avg sfc latent heat flux (W/m^2)
real :: rshort_avg_ll      (nlon,nlat) ! avg sfc downward s/w rad flux (W/m^2)
real :: rshortup_avg_ll    (nlon,nlat) ! avg sfc upward s/w rad flux (W/m^2)
real :: rlong_avg_ll       (nlon,nlat) ! avg sfc downward l/w rad flux (W/m^2)
real :: rlongup_avg_ll     (nlon,nlat) ! avg sfc upward l/w rad flux (W/m^2)
real :: rshort_top_avg_ll  (nlon,nlat) ! avg TOA downward s/w rad flux (W/m^2)
real :: rshortup_top_avg_ll(nlon,nlat) ! avg TOA upward s/w rad flux (W/m^2)
real :: rlongup_top_avg_ll (nlon,nlat) ! avg TOA upward l/w rad flux (W/m^2)

!----------
! 3D FIELDS
!----------

real :: u_ll (nlon,nlat,mza-2) ! zonal wind component (m/s)
real :: v_ll (nlon,nlat,mza-2) ! meridional wind component (m/s)
real :: t_ll (nlon,nlat,mza-2) ! air temperature (K)
real :: r_ll (nlon,nlat,mza-2) ! water vapor mixing ratio (kg/kg)
real :: p_ll (nlon,nlat,mza-2) ! air pressure (Pa)

integer :: k,iw,ilat,ilon
integer :: npoly,kb

integer :: ndims, idims(3)

real :: raxis,raxisi
real :: dlon,dlat

real :: scr1a(mwa), scr1b(mwa)
real :: scr2a(mza,mwa), scr2b(mza,mwa)
real :: scr2_ll(nlon,nlat), scr3_ll(nlon,nlat,mza-2)
real :: uzonal(mza,mua), umerid(mza,mua)

! SEIGEL 2013 - Added for ll interp writeout
character(pathlen) :: hnamel

! Initialize latitude-longitude arrays to zero prior to interpolation.
! For certain applications, zero should be replaced with "missing value".

topo_ll            (:,:) = 0.
u_sfc_ll           (:,:) = 0.
v_sfc_ll           (:,:) = 0.
t_sfc_ll           (:,:) = 0.
r_sfc_ll           (:,:) = 0.
pvap_sfc_ll        (:,:) = 0.
slp_ll             (:,:) = 0.

accpmic_ll           (:,:) = 0.
accpcon_ll           (:,:) = 0.
vapflux_accum_ll     (:,:) = 0.
sensflux_accum_ll    (:,:) = 0.
latflux_accum_ll     (:,:) = 0.
rshort_accum_ll      (:,:) = 0.
rshortup_accum_ll    (:,:) = 0.
rlong_accum_ll       (:,:) = 0.
rlongup_accum_ll     (:,:) = 0.
rshort_top_accum_ll  (:,:) = 0.
rshortup_top_accum_ll(:,:) = 0.
rlongup_top_accum_ll (:,:) = 0.

vapflux_avg_ll     (:,:) = 0.
sensflux_avg_ll    (:,:) = 0.
latflux_avg_ll     (:,:) = 0.
rshort_avg_ll      (:,:) = 0.
rshortup_avg_ll    (:,:) = 0.
rlong_avg_ll       (:,:) = 0.
rlongup_avg_ll     (:,:) = 0.
rshort_top_avg_ll  (:,:) = 0.
rshortup_top_avg_ll(:,:) = 0.
rlongup_top_avg_ll (:,:) = 0.

u_ll(:,:,:) = 0.
v_ll(:,:,:) = 0.
t_ll(:,:,:) = 0.
r_ll(:,:,:) = 0.
p_ll(:,:,:) = 0.

! Compute latitude and longitude of output grid points (assuming uniform spacing)

dlat = (endlat  - beglat) / real(nlat-1)

do ilat = 1,nlat
   alat(ilat) = beglat + dlat * real(ilat-1)
enddo

if (endlon >= beglon) then
   dlon = (endlon - beglon) / real(nlon-1)
else
   dlon = (endlon + 360. - beglon) / real(nlon-1)
endif

do ilon = 1,nlon
   alon(ilon) = beglon + dlon * real(ilon-1)
   if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.
enddo

! Fill zlev with model grid levels

do k = 2,mza-1
   zlev(k-1) = zt(k)
enddo

!------------------------------------------------------------
! Interpolate topography height
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,topw,topo_ll)

!------------------------------------------------------------
! Compute zonal and meridional wind components on OLAM grid
! and copy their values at lowest prognosed model level to
! separate arrays
!------------------------------------------------------------

do iw = 2,mwa

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

   uzonal(:,iw) = 0.
   umerid(:,iw) = 0.

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb,mza-1
         scr2a(k,iw) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) * raxisi
         scr2b(k,iw) = vze(k,iw) * raxis * eradi &
            - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) * zew(iw) * raxisi * eradi
      enddo

   endif

   scr1a(iw) = scr2a(kb,iw)
   scr1b(iw) = scr2b(kb,iw)

enddo

!------------------------------------------------------------
! Interpolate zonal and meridional wind components, and their
! near-surface values
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr2a,u_ll)
call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr2b,v_ll)

call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,u_sfc_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1b,v_sfc_ll)

!------------------------------------------------------------
! Compute temperature on OLAM grid and copy its value at 
! lowest prognosed model level to separate array.
! Compute sea level pressure based on values at lowest
! prognosed model level. 
!------------------------------------------------------------

do iw = 2,mwa
   kb = lpw(iw)

   do k = kb,mza-1
      scr2a(k,iw) = theta(k,iw) * (press(k,iw) / p00) ** rocp
   enddo

   scr1a(iw) = scr2a(kb,iw)

   scr1b(iw) = press(kb,iw) &
             * (1. - .0065 * zt(kb) / (scr2a(kb,iw) + .0065 * zt(kb)))**(-5.257)
enddo

!------------------------------------------------------------
! Interpolate temperature, its surface value, and sea level pressure
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr2a,t_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,t_sfc_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1b,slp_ll)

!------------------------------------------------------------
! Compute vapor mixing ratio and copy its value at lowest
! prognosed model level on OLAM grid to separate array.
! Compute vapor pressure at lowest prognosed model level on OLAM grid.
!------------------------------------------------------------

do iw = 2,mwa
   kb = lpw(iw)

   do k = kb,mza-1
      scr2a(k,iw) = sh_v(k,iw) / (1. - sh_w(k,iw))
   enddo

   scr1a(iw) = scr2a(kb,iw)
   scr1b(iw) = sh_v(kb,iw) * rho(kb,iw) * rvap &
             * theta(kb,iw) * (press(kb,iw) / p00) ** rocp

enddo

!------------------------------------------------------------
! Interpolate vapor mixing ratio, its surface value, and
! the surface value of vapor pressure
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr2a,r_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,r_sfc_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1b,pvap_sfc_ll)

!------------------------------------------------------------
! Interpolate atmospheric pressure
!------------------------------------------------------------

! Copy pressure to real array

scr2a(:,:) = real(press(:,:))

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr2a,p_ll)

!------------------------------------------------------------
! Compute total (resolved + parameterized) accumulated
! precipitation on OLAM grid
!------------------------------------------------------------

do iw = 2,mwa
   scr1a(iw) = 0.

   if (allocated(accpd))  scr1a(iw) = scr1a(iw) + real(accpd(iw))
   if (allocated(accpr))  scr1a(iw) = scr1a(iw) + real(accpr(iw))
   if (allocated(accpp))  scr1a(iw) = scr1a(iw) + real(accpp(iw))
   if (allocated(accps))  scr1a(iw) = scr1a(iw) + real(accps(iw))
   if (allocated(accpa))  scr1a(iw) = scr1a(iw) + real(accpa(iw))
   if (allocated(accpg))  scr1a(iw) = scr1a(iw) + real(accpg(iw))
   if (allocated(accph))  scr1a(iw) = scr1a(iw) + real(accph(iw))
enddo

!------------------------------------------------------------
! Interpolate accumulated microphysics and convective precipitation
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,accpmic_ll)

if (allocated(aconpr)) then
   scr1a(:) = real(aconpr(:))
   call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,accpcon_ll)
endif

!------------------------------------------------------------
! Compute (as necessary) and interpolate time-integrated surface fluxes
!------------------------------------------------------------

scr1a(:) = real(sflux_r_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,vapflux_accum_ll)

scr1a(:) = real(sflux_t_accum(:)) * cp
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,sensflux_accum_ll)

scr1a(:) = real(sflux_r_accum(:)) * alvl
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,latflux_accum_ll)

scr1a(:) = real(rshort_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rshort_accum_ll)

scr1a(:) = real(rshortup_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rshortup_accum_ll)

scr1a(:) = real(rlong_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rlong_accum_ll)

scr1a(:) = real(rlongup_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rlongup_accum_ll)

!------------------------------------------------------------
! Interpolate time-averaged top-of-atmosphere radiative fluxes
!------------------------------------------------------------

scr1a(:) = real(rshort_top_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rshort_top_accum_ll)

scr1a(:) = real(rshortup_top_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rshortup_top_accum_ll)

scr1a(:) = real(rlongup_top_accum(:))
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,rlongup_top_accum_ll)

!------------------------------------------------------------
! Compute (as necessary) and interpolate time-averaged surface fluxes
!------------------------------------------------------------

scr1a(:) = sflux_t_avg(:) * cp
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,sensflux_avg_ll)

scr1a(:) = sflux_r_avg(:) * alvl
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,latflux_avg_ll)

call interp_htw_ll(nlon,nlat,1,1,alon,alat,sflux_r_avg,vapflux_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshort_avg,rshort_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshortup_avg,rshortup_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlong_avg,rlong_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlongup_avg,rlongup_avg_ll)

!------------------------------------------------------------
! Interpolate time-averaged top-of-atmosphere radiative fluxes
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshort_top_avg,rshort_top_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshortup_top_avg,rshortup_top_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlongup_top_avg,rlongup_top_avg_ll)

! HDF5 write

  ! SEIGEL 2013 - add write for LL interpolation  
  write(io6,*) 'GOT HERE'
  call makefnam(hnamel, hfilepref, current_time, 'LL', '$', 'h5')
  call shdf5_open(hnamel,'W',iclobber) 

  ndims = 1
  idims(2) = 1
  idims(3) = 1

  idims(1) = nlat

  CALL shdf5_orec(ndims, idims, 'ALAT', rvara=alat)

  idims(1) = nlon

  CALL shdf5_orec(ndims, idims, 'ALON', rvara=alon)

  idims(1) = mza-2

  CALL shdf5_orec(ndims, idims, 'ZLEV', rvara=zlev)

  ndims = 2
  !idims(1) = nlat
  !idims(2) = nlon
  idims(1) = nlon
  idims(2) = nlat

  CALL shdf5_orec(ndims, idims, 'TOPO_LL'            , rvara=topo_ll)
  CALL shdf5_orec(ndims, idims, 'U_SFC_LL'           , rvara=u_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'V_SFC_LL'           , rvara=v_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'T_SFC_LL'           , rvara=t_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'R_SFC_LL'           , rvara=r_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'SLP_LL'             , rvara=slp_ll)
  CALL shdf5_orec(ndims, idims, 'PVAP_SFC_LL'        , rvara=pvap_sfc_ll)

  CALL shdf5_orec(ndims, idims, 'ACCPMIC_LL'           , rvara=accpmic_ll)
  CALL shdf5_orec(ndims, idims, 'ACCPCON_LL'           , rvara=accpcon_ll)
  CALL shdf5_orec(ndims, idims, 'VAPFLUX_ACCUM_LL'     , rvara=vapflux_accum_ll)
  CALL shdf5_orec(ndims, idims, 'SENSFLUX_ACCUM_LL'    , rvara=sensflux_accum_ll)
  CALL shdf5_orec(ndims, idims, 'LATFLUX_ACCUM_LL'     , rvara=latflux_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_ACCUM_LL'      , rvara=rshort_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_ACCUM_LL'    , rvara=rshortup_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RLONG_ACCUM_LL'       , rvara=rlong_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_ACCUM_LL'     , rvara=rlongup_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_TOP_ACCUM_LL'  , rvara=rshort_top_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_TOP_ACCUM_LL', rvara=rshortup_top_accum_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_TOP_ACCUM_LL' , rvara=rlongup_top_accum_ll)

  CALL shdf5_orec(ndims, idims, 'VAPFLUX_AVG_LL'     , rvara=vapflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'SENSFLUX_AVG_LL'    , rvara=sensflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'LATFLUX_AVG_LL'     , rvara=latflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_AVG_LL'      , rvara=rshort_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_AVG_LL'    , rvara=rshortup_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONG_AVG_LL'       , rvara=rlong_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_AVG_LL'     , rvara=rlongup_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_TOP_AVG_LL'  , rvara=rshort_top_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_TOP_AVG_LL', rvara=rshortup_top_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_TOP_AVG_LL' , rvara=rlongup_top_avg_ll)

  ndims = 3
  !idims(1) = nlat
  !idims(2) = nlon
  idims(1) = nlon
  idims(2) = nlat
  idims(3) = mza-2

  CALL shdf5_orec(ndims, idims, 'U_LL', rvara=u_ll)
  CALL shdf5_orec(ndims, idims, 'V_LL', rvara=v_ll)
  CALL shdf5_orec(ndims, idims, 'T_LL', rvara=t_ll)
  CALL shdf5_orec(ndims, idims, 'R_LL', rvara=r_ll)
  CALL shdf5_orec(ndims, idims, 'P_LL', rvara=p_ll)

  ! SEIGEL 2013 - close interpolation write
  call shdf5_close()

!------------------------------------------------------------
! Plot test
!------------------------------------------------------------

! Reopen the current graphics output workstation if it is closed

 call o_reopnwk()

 call contplot_ll(nlon,nlat,1,1,alon,alat,topo_ll,402)
 call contplot_ll(nlon,nlat,1,1,alon,alat,u_sfc_ll,113)
 call contplot_ll(nlon,nlat,1,1,alon,alat,v_sfc_ll,113)
 call contplot_ll(nlon,nlat,1,1,alon,alat,t_sfc_ll,2)
 scr2_ll(:,:) = r_sfc_ll(:,:) * 1000. 
 call contplot_ll(nlon,nlat,1,1,alon,alat,scr2_ll,5)
 scr2_ll(:,:) = slp_ll(:,:) * .01 
 call contplot_ll(nlon,nlat,1,1,alon,alat,scr2_ll,17)
 call contplot_ll(nlon,nlat,1,1,alon,alat,pvap_sfc_ll,35)

 call contplot_ll(nlon,nlat,1,1,alon,alat,accpmic_ll,54)
 call contplot_ll(nlon,nlat,1,1,alon,alat,accpcon_ll,54)

 call contplot_ll(nlon,nlat,1,1,alon,alat,vapflux_accum_ll,306)
 call contplot_ll(nlon,nlat,1,1,alon,alat,sensflux_accum_ll,163)
 call contplot_ll(nlon,nlat,1,1,alon,alat,latflux_accum_ll,163)

 call contplot_ll(nlon,nlat,1,1,alon,alat,rshort_accum_ll,63)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rshortup_accum_ll,63)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlong_accum_ll,63)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlongup_accum_ll,63)

 call contplot_ll(nlon,nlat,1,1,alon,alat,rshort_top_accum_ll,63)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rshortup_top_accum_ll,63)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlongup_top_accum_ll,63)

 call contplot_ll(nlon,nlat,1,1,alon,alat,vapflux_avg_ll,306)
 call contplot_ll(nlon,nlat,1,1,alon,alat,sensflux_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,latflux_avg_ll,102)

 call contplot_ll(nlon,nlat,1,1,alon,alat,rshort_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rshortup_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlong_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlongup_avg_ll,102)

 call contplot_ll(nlon,nlat,1,1,alon,alat,rshort_top_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rshortup_top_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,rlongup_top_avg_ll,102)

 call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,u_ll,113)
 call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,v_ll,113)
 call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,t_ll,2)
 scr3_ll(:,:,:) = r_ll(:,:,:) * 1000. 
 call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,scr3_ll,5)
 scr3_ll(:,:,:) = p_ll(:,:,:) * .01 
 call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,scr3_ll,8)

 call o_clswk()

return
end subroutine fields_ll

!==========================================================================

subroutine tileplot_ll(nlon,nlat,nlev,ilev,alon,alat,fld,itab)

use plotcolors,  only: clrtab
use oplot_coms, only: op

implicit none

integer, intent(in) :: nlon,nlat,nlev,ilev,itab
real, intent(in) :: alon(nlon),alat(nlat)
real, intent(in) :: fld(nlon,nlat,nlev)

real :: glon(nlon)
real :: htpn(4),vtpn(4),fldval
integer :: ilat, ilon, ival, icolor

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

do ilat = 1,nlat

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

   do ilon = 1,nlon

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

      fldval = fld(ilon,ilat,ilev)

      ival = 1
      do while (fldval > clrtab(itab)%vals(ival) .and.  &
                  ival < clrtab(itab)%nvals             )
         ival = ival + 1
      enddo
      icolor = clrtab(itab)%ipal(ival)
      
      call fillpolyg(4,htpn,vtpn,icolor)

   enddo
enddo

call o_frame()

return
end subroutine tileplot_ll

!==========================================================================

subroutine contplot_ll(nlon,nlat,nlev,ilev,alon,alat,fld,itab)

use plotcolors, only: clrtab
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

return
end subroutine contplot_ll
