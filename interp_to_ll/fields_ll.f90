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

!--------------------------------------------------------------------------------
! NOTE: Some of the fields included in this template are time-averaged surface
! or top-of-atmosphere fluxes.  OLAM currently evaluates these as averages over
! the time interval between consecutive history-file writes, and it therefore 
! resets the averages to zero after each history write.  Consequently, for these
! fields to have their correct meaning when interpolated to the latitude-longitude
! grid, the interpolation should be done only at times when history file writes 
! are performed.  Both the accumulation and resetting of averages are performed
! in subroutine accum_timeavg in the source code file mem_timeavg.f90.
!--------------------------------------------------------------------------------

use mem_ijtabs,  only: itab_w, itab_u, itab_v
use mem_basic,   only: uc, vc, wc, rho, press, theta, sh_w, sh_v

use mem_grid,    only: mza, mua, mva, mwa, lpu, lpv, lpw, &
                       xem, yem, zem, xeu, yeu, zeu, &
                       xev, yev, zev, xew, yew, zew, &
                       topw, glatm, glonm, glatw, glonw, zm, zt, &
                       unx, uny, unz, vnx, vny, vnz, dzt

use misc_coms,   only: io6, meshtype

use mem_micro,   only: accpd, accpr, accpp, accps, accpa, accpg, accph

use mem_cuparm,  only: aconpr

use consts_coms, only: p00, rocp, piu180, erad, eradi, pio180, cp, alvl, rvap

use mem_timeavg, only: rshort_avg, rshortup_avg, rlong_avg, rlongup_avg, &
                       rshort_top_avg, rshortup_top_avg, rlongup_top_avg, &
                       sflux_t_avg, sflux_r_avg

use hdf5_utils,  only: shdf5_orec

implicit none

!--------------------------------------------------------------------------------
! THE USER MUST SPECIFY THE FOLLOWING PARAMETERS THAT DEFINE THE OUTPUT
! LATITUDE-LONGITUDE ARRAYS. LONGITUDE VALUES ARE DEFINED TO BE IN THE
! RANGE (-180,180], WHICH DOES NOT PRECLUDE THAT THE LIMITED-AREA
! LATITUDE-LONGITUDE GRID CROSS THE 180 DEGREE MERIDIAN.
!--------------------------------------------------------------------------------

integer, parameter :: nlon = 51   ! number of longitude values
integer, parameter :: nlat = 51   ! number of latitude values

real, parameter :: beglon = -130. ! minimum (westernmost) longitude (deg)
real, parameter :: endlon = -80.  ! maximum (easternmost) longitude (deg)
real, parameter :: beglat = 10.   ! minimum (southernmost) latitude (deg)
real, parameter :: endlat = 60.   ! maximum (northernmost) latitude (deg)

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
real :: atotpr_ll          (nlon,nlat) ! accumulated precipitation (kg/m2)
real :: sensflux_avg_ll    (nlon,nlat) ! TA surface sensible heat flux (W/m2)
real :: latflux_avg_ll     (nlon,nlat) ! TA surface latent heat flux (W/m2)
real :: vapflux_avg_ll     (nlon,nlat) ! TA surface vapor flux (kg/m2/s)
real :: rshort_avg_ll      (nlon,nlat) ! TA surface downward s/w rad flux (W/m2)
real :: rshortup_avg_ll    (nlon,nlat) ! TA surface upward s/w rad flux (W/m2)
real :: rlong_avg_ll       (nlon,nlat) ! TA surface downward l/w rad flux (W/m2)
real :: rlongup_avg_ll     (nlon,nlat) ! TA surface upward l/w rad flux (W/m2)
real :: rshort_top_avg_ll  (nlon,nlat) ! TA top of atm downward s/w rad flux (W/m2)
real :: rshortup_top_avg_ll(nlon,nlat) ! TA top of atm upward s/w rad flux (W/m2)
real :: rlongup_top_avg_ll (nlon,nlat) ! TA top of atm upward l/w rad flux (W/m2)

!----------
! 3D FIELDS
!----------

real :: u_ll (nlon,nlat,mza-2) ! zonal wind component (m/s)
real :: v_ll (nlon,nlat,mza-2) ! meridional wind component (m/s)
real :: t_ll (nlon,nlat,mza-2) ! air temperature (K)
real :: r_ll (nlon,nlat,mza-2) ! water vapor mixing ratio (kg/kg)
real :: p_ll (nlon,nlat,mza-2) ! air pressure (Pa)

integer :: k,iv,iw,lpuv,ilat,ilon
integer :: npoly,kb,jv

integer :: ndims, idims(3)

real :: vx,vy,vz,raxis,raxisi,u,v
real :: dlon,dlat
real :: xeuv,yeuv,zeuv,farv2

real :: scr1a(mwa), scr1b(mwa)
real :: scr2a(mza,mwa), scr2b(mza,mwa)
real :: scr2_ll(nlon,nlat), scr3_ll(nlon,nlat,mza-2)

real :: uzonal(mza,mua), umerid(mza,mua)

real :: vxe(mza),vye(mza),vze(mza)
real :: alat1,alat2

! Initialize latitude-longitude arrays to zero prior to interpolation.
! For certain applications, zero should be replaced with "missing value".

topo_ll            (:,:) = 0.
u_sfc_ll           (:,:) = 0.
v_sfc_ll           (:,:) = 0.
t_sfc_ll           (:,:) = 0.
r_sfc_ll           (:,:) = 0.
pvap_sfc_ll        (:,:) = 0.
slp_ll             (:,:) = 0.
atotpr_ll          (:,:) = 0.
sensflux_avg_ll    (:,:) = 0.
latflux_avg_ll     (:,:) = 0.
vapflux_avg_ll     (:,:) = 0.
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
enddo

if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.

! Fill zlev with model grid levels

do k = 2,mza
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

   vxe(:) = 0.
   vye(:) = 0.
   vze(:) = 0.

   uzonal(:,iw) = 0.
   umerid(:,iw) = 0.

   do jv = 1,npoly

      if (meshtype == 1) then
      
         iv = itab_w(iw)%iu(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + itab_w(iw)%vxu(jv) * uc(k,iv)
            vye(k) = vye(k) + itab_w(iw)%vyu(jv) * uc(k,iv)
            vze(k) = vze(k) + itab_w(iw)%vzu(jv) * uc(k,iv)
         enddo

      else

         iv = itab_w(iw)%iv(jv)
         farv2 = 2. * itab_w(iw)%farv(jv)

         do k = kb,mza-1
            vxe(k) = vxe(k) + farv2 * vc(k,iv) * vnx(iv)
            vye(k) = vye(k) + farv2 * vc(k,iv) * vny(iv)
            vze(k) = vze(k) + farv2 * vc(k,iv) * vnz(iv)
         enddo

      endif

   enddo

   raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

! Evaluate zonal and meridional wind components from model

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = kb,mza-1
         scr2a(k,iw) = (vye(k) * xew(iw) - vxe(k) * yew(iw)) * raxisi
         scr2b(k,iw) = vze(k) * raxis * eradi &
            - (vxe(k) * xew(iw) + vye(k) * yew(iw)) * zew(iw) * raxisi * eradi
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

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,press,p_ll)

!------------------------------------------------------------
! Compute total (resolved + parameterized) accumulated
! precipitation on OLAM grid
!------------------------------------------------------------

do iw = 2,mwa
   scr1a(iw) = 0.

   if (allocated(accpd))  scr1a(iw) = scr1a(iw) + accpd(iw)
   if (allocated(accpr))  scr1a(iw) = scr1a(iw) + accpr(iw)
   if (allocated(accpp))  scr1a(iw) = scr1a(iw) + accpp(iw)
   if (allocated(accps))  scr1a(iw) = scr1a(iw) + accps(iw)
   if (allocated(accpa))  scr1a(iw) = scr1a(iw) + accpa(iw)
   if (allocated(accpg))  scr1a(iw) = scr1a(iw) + accpg(iw)
   if (allocated(accph))  scr1a(iw) = scr1a(iw) + accph(iw)
   if (allocated(aconpr)) scr1a(iw) = scr1a(iw) + aconpr(iw)
enddo

!------------------------------------------------------------
! Interpolate total accumulated precipitation
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1a,atotpr_ll)

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
  idims(1) = nlat
  idims(2) = nlon

  CALL shdf5_orec(ndims, idims, 'TOPO_LL'            , rvara=topo_ll)
  CALL shdf5_orec(ndims, idims, 'U_SFc_LL'           , rvara=u_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'V_SFC_LL'           , rvara=v_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'T_SFC_LL'           , rvara=t_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'R_SFC_LL'           , rvara=r_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'SLP_LL'             , rvara=slp_ll)
  CALL shdf5_orec(ndims, idims, 'PVAP_SFC_LL'        , rvara=pvap_sfc_ll)
  CALL shdf5_orec(ndims, idims, 'ATOTPR_LL'          , rvara=atotpr_ll)
  CALL shdf5_orec(ndims, idims, 'SENSFLUX_AVG_LL'    , rvara=sensflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'LATFLUX_AVG_LL'     , rvara=latflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'VAPFLUX_AVG_LL'     , rvara=vapflux_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_AVG_LL'      , rvara=rshort_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_AVG_LL'    , rvara=rshortup_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONG_AVG_LL'       , rvara=rlong_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_AVG_LL'     , rvara=rlongup_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORT_TOP_AVG_LL'  , rvara=rshort_top_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RSHORTUP_TOP_AVG_LL', rvara=rshortup_top_avg_ll)
  CALL shdf5_orec(ndims, idims, 'RLONGUP_TOP_AVG_LL' , RVARA=RLONGUP_TOP_AVG_LL)

  ndims = 3
  idims(1) = nlat
  idims(2) = nlon
  idims(3) = mza-2

  CALL shdf5_orec(ndims, idims, 'U_LL', rvara=u_ll)
  CALL shdf5_orec(ndims, idims, 'V_LL', rvara=v_ll)
  CALL shdf5_orec(ndims, idims, 'T_LL', rvara=t_ll)
  CALL shdf5_orec(ndims, idims, 'R_LL', rvara=r_ll)
  CALL shdf5_orec(ndims, idims, 'P_LL', rvara=p_ll)

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
 call contplot_ll(nlon,nlat,1,1,alon,alat,atotpr_ll,54)

 call contplot_ll(nlon,nlat,1,1,alon,alat,sensflux_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,latflux_avg_ll,102)
 call contplot_ll(nlon,nlat,1,1,alon,alat,vapflux_avg_ll,306)

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
