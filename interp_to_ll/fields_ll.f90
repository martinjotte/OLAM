subroutine fields_ll()

! This subroutine is a template, intended for user modification, for interpolating
! selected fields from the OLAM grid to a structured latitude-longitude grid
! of either limited or global area.  

! The following tasks are performed:

! (1) Allocate arrays for interpolated fields
! (2) Allocate scratch array
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
use mem_basic,   only: uc, vc, wc, rho, press, theta

use mem_grid,    only: mza, mua, mva, mwa, lpu, lpv, lpw, &
                       xem, yem, zem, xeu, yeu, zeu, &
                       xev, yev, zev, xew, yew, zew, &
                       topw, glatm, glonm, glatw, glonw, zm, zt, &
                       unx, uny, unz, vnx, vny, vnz, dzt

use misc_coms,   only: io6, meshtype

use mem_micro,   only: accpd, accpr, accpp, accps, accpa, accpg, accph
                      
use mem_cuparm,  only: aconpr

use consts_coms, only: p00, rocp, piu180, erad, pio180, cp, alvl

use mem_timeavg, only: rshort_avg, rshortup_avg, rlong_avg, rlongup_avg, &
                       rshort_top_avg, rshortup_top_avg, rlongup_top_avg, &
                       sflux_t_avg, sflux_r_avg

implicit none

!--------------------------------------------------------------------------------
! THE USER MUST SPECIFY THE FOLLOWING PARAMETERS THAT DEFINE THE OUTPUT
! LATITUDE-LONGITUDE ARRAYS; LONGITUDES ARE ASSUMED TO BE IN THE RANGE (-180,180].
!--------------------------------------------------------------------------------

integer, parameter :: nlon = 151   ! number of longitude values
integer, parameter :: nlat = 151   ! number of latitude values

real, parameter :: beglon = -85.64 ! minimum (westernmost) longitude (deg)
real, parameter :: endlon = -64.36  ! maximum (easternmost) longitude (deg)
real, parameter :: beglat = 10.   ! minimum (southernmost) latitude (deg)
real, parameter :: endlat = 30.   ! maximum (northernmost) latitude (deg)

!--------------------------------------------------------------------------------
! THE FOLLOWING ARRAYS ARE DIMENSIONED USING THE ABOVE NLON,NLAT PARAMETERS,
! AND WILL CONTAIN THE INTERPOLATED FIELDS.  THE USER SHOULD ADD NEW ARRAYS
! AS REQUIRED.
!--------------------------------------------------------------------------------

real :: t_ll (nlon,nlat,mza-2)
real :: u_ll (nlon,nlat,mza-2)
real :: v_ll (nlon,nlat,mza-2)

real :: topo_ll(nlon,nlat)

real :: atotpr_ll(nlon,nlat)

real :: sensflux_avg_ll    (nlon,nlat)
real :: latflux_avg_ll     (nlon,nlat)
real :: vapflux_avg_ll     (nlon,nlat)

real :: rshort_avg_ll      (nlon,nlat)
real :: rshortup_avg_ll    (nlon,nlat)
real :: rlong_avg_ll       (nlon,nlat)
real :: rlongup_avg_ll     (nlon,nlat)

real :: rshort_top_avg_ll  (nlon,nlat)
real :: rshortup_top_avg_ll(nlon,nlat)
real :: rlongup_top_avg_ll (nlon,nlat)

real :: alon(nlon)
real :: alat(nlat)
real :: zlev(mza-2)

integer :: k,iv,iw,lpuv,ilat,ilon
real :: vx,vy,vz,raxis,u,v
real :: dlon,dlat
real :: xeuv,yeuv,zeuv

real :: scr(mza,mwa), scr1(mwa)

real :: uzonal(mza,mua), umerid(mza,mua)

real :: alat1,alat2,dnorm,dnorm1,dnorm2,anorm1,anorm2,ubar

! Initialize latitude-longitude arrays to zero prior to interpolation.
! For certain applications, zero should be replaced with "missing value".

t_ll(:,:,:) = 0.
u_ll(:,:,:) = 0.
v_ll(:,:,:) = 0.

topo_ll(:,:) = 0.

atotpr_ll(:,:) = 0.

sensflux_avg_ll(:,:) = 0.
latflux_avg_ll(:,:) = 0.
vapflux_avg_ll(:,:) = 0.

rshort_avg_ll(:,:) = 0.
rshortup_avg_ll(:,:) = 0.
rlong_avg_ll(:,:) = 0.
rlongup_avg_ll(:,:) = 0.

rshort_top_avg_ll(:,:) = 0.
rshortup_top_avg_ll(:,:) = 0.
rlongup_top_avg_ll(:,:) = 0.

! Compute longitude and latitude of output grid points (assuming uniform spacing)

if (endlon >= beglon) then
   dlon = (endlon - beglon) / real(nlon-1)
else
   dlon = (endlon + 360. - beglon) / real(nlon-1)
endif

dlat = (endlat  - beglat) / real(nlat-1)

do ilon = 1,nlon
   alon(ilon) = beglon + dlon * real(ilon-1)
enddo

if (alon(ilon) > 180.) alon(ilon) = alon(ilon) - 360.

do ilat = 1,nlat
   alat(ilat) = beglat + dlat * real(ilat-1)
enddo

!------------------------------------------------------------
! Temperature
!------------------------------------------------------------

scr(:,:) = 0.

do iw = 2,mwa
   do k = lpw(iw),mza-1
      scr(k,iw) = theta(k,iw) * (press(k,iw) / p00) ** rocp
   enddo
enddo

call interp_htw_ll(nlon,nlat,mza,mza-2,alon,alat,scr,t_ll)

!------------------------------------------------------------
! UZONAL and UMERID wind components at U/V point
!------------------------------------------------------------

do iv = 2,mva
   uzonal(:,iv) = 0.
   umerid(:,iv) = 0.
   if (meshtype == 1) then
      lpuv = lpu(iv)
      xeuv = xeu(iv)
      yeuv = yeu(iv)
      zeuv = zeu(iv)
   else
      lpuv = lpv(iv)
      xeuv = xev(iv)
      yeuv = yev(iv)
      zeuv = zev(iv)
   endif

   raxis = sqrt(xeuv ** 2 + yeuv ** 2) ! dist from earth axis

   if (raxis > 1.e3) then
      do k = lpuv,mza-1

         vx = unx(iv) * uc(k,iv) + vnx(iv) * vc(k,iv)
         vy = uny(iv) * uc(k,iv) + vny(iv) * vc(k,iv)
         vz = unz(iv) * uc(k,iv) + vnz(iv) * vc(k,iv)

         uzonal(k,iv) = (vy * xeuv - vx * yeuv) / raxis
         umerid(k,iv) = vz * raxis / erad  &
            - (vx * xeuv + vy * yeuv) * zeuv / (raxis * erad)

      enddo
   endif

enddo

call interp_hvn_ll(nlon,nlat,mza,mza-2,alon,alat,uzonal,u_ll)
call interp_hvn_ll(nlon,nlat,mza,mza-2,alon,alat,umerid,v_ll)

!------------------------------------------------------------
! TOPO (topography height)
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,topw,topo_ll)

!------------------------------------------------------------
! Total (resolved + parameterized) accumulated surface precipitation
!------------------------------------------------------------

do iw = 2,mwa
   scr1(iw) = 0.

   if (allocated(accpd))  scr1(iw) = scr1(iw) + accpd(iw)
   if (allocated(accpr))  scr1(iw) = scr1(iw) + accpr(iw)
   if (allocated(accpp))  scr1(iw) = scr1(iw) + accpp(iw)
   if (allocated(accps))  scr1(iw) = scr1(iw) + accps(iw)
   if (allocated(accpa))  scr1(iw) = scr1(iw) + accpa(iw)
   if (allocated(accpg))  scr1(iw) = scr1(iw) + accpg(iw)
   if (allocated(accph))  scr1(iw) = scr1(iw) + accph(iw)
   if (allocated(aconpr)) scr1(iw) = scr1(iw) + aconpr(iw)
enddo

call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1,atotpr_ll)

!------------------------------------------------------------
! Time-averaged surface fluxes
!------------------------------------------------------------

scr1(:) = sflux_t_avg(:) * cp

call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1,sensflux_avg_ll)

scr1(:) = sflux_r_avg(:) * alvl
call interp_htw_ll(nlon,nlat,1,1,alon,alat,scr1,latflux_avg_ll)

call interp_htw_ll(nlon,nlat,1,1,alon,alat,sflux_r_avg,vapflux_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshort_avg,rshort_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshortup_avg,rshortup_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlong_avg,rlong_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlongup_avg,rlongup_avg_ll)

!------------------------------------------------------------
! Time-averaged top-of-atmosphere radiative fluxes
!------------------------------------------------------------

call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshort_top_avg,rshort_top_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rshortup_top_avg,rshortup_top_avg_ll)
call interp_htw_ll(nlon,nlat,1,1,alon,alat,rlongup_top_avg,rlongup_top_avg_ll)

!------------------------------------------------------------
! Plot test
!------------------------------------------------------------

! Reopen the current graphics output workstation if it is closed

call o_reopnwk()

call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,t_ll,2)
call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,u_ll,113)
call tileplot_ll(nlon,nlat,mza-2,1,alon,alat,v_ll,113)

call contplot_ll(nlon,nlat,mza-2,1,alon,alat,t_ll,2)
call contplot_ll(nlon,nlat,mza-2,1,alon,alat,u_ll,113)
call contplot_ll(nlon,nlat,mza-2,1,alon,alat,v_ll,113)

call contplot_ll(nlon,nlat,1,1,alon,alat,topo_ll,402)
call contplot_ll(nlon,nlat,1,1,alon,alat,atotpr_ll,5)

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
