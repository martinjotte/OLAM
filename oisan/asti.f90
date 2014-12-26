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
subroutine isnstage(p_u, p_v, p_t, p_z, p_r, &
                    p_topo, p_prsfc, p_tsfc, p_shsfc, &
                    o_rho, o_theta, o_shv, o_uzonal, o_umerid, o_vc)

use max_dims,   only: maxpr
use isan_coms,  only: nprz, npry, nprx, nprz_rh, &
                      xswlat, xswlon, gdatdx, gdatdy, &
                      npd, kzonoff, levpr, lzon_bot, ipoffset
use mem_grid,   only: glatw, glonw, mza, mwa, mva, lpv, lpw, &
                      xev, yev, zev, unx, uny, unz, vnx, vny, vnz
use mem_ijtabs, only: jtab_v, jtab_w, itab_v, itab_w, jtv_init, jtw_init
use mem_zonavg, only: zonp_vect, zont, zonz, zonr, zonu
use consts_coms,only: r8, eradi, rocp, p00i, cp
use misc_coms,  only: io6, iparallel

use olam_mpi_atm, only: mpi_send_w, mpi_send_v, mpi_recv_w, mpi_recv_v
use obnd,         only: lbcopy_v

implicit none

real, intent(in) :: p_u(nprx+4,npry+4,nprz)
real, intent(in) :: p_v(nprx+4,npry+4,nprz)
real, intent(in) :: p_t(nprx+4,npry+4,nprz)
real, intent(in) :: p_z(nprx+4,npry+4,nprz)
real, intent(in) :: p_r(nprx+4,npry+4,nprz)

real, intent(in) :: p_topo (nprx+4,npry+4)
real, intent(in) :: p_prsfc(nprx+4,npry+4)
real, intent(in) :: p_tsfc (nprx+4,npry+4)
real, intent(in) :: p_shsfc(nprx+4,npry+4)

real(r8), intent(inout) :: o_rho   (mza,mwa)
real,     intent(inout) :: o_theta (mza,mwa)
real,     intent(inout) :: o_shv   (mza,mwa)
real,     intent(inout) :: o_uzonal(mza,mwa)
real,     intent(inout) :: o_umerid(mza,mwa)
real,     intent(inout) :: o_vc    (mza,mva)

real :: pcol_p    (maxpr+2)
real :: pcol_temp (maxpr+2)
real :: pcol_z    (maxpr+2)
real :: pcol_u    (maxpr+2)
real :: pcol_v    (maxpr+2)
real :: pcol_rt   (maxpr+2)
real :: pcol_r    (maxpr+2)
real :: pcol_exner(maxpr+2)

real :: pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc

character(3) :: csuff
character(8) :: rot_type

integer :: ngrd,lv,lf,k,levp,j,iw,ilat,iu,iv,iuv,mrl,ka
integer :: iw1,iw2,nlevs,kstrt

real :: wt2,grx,gry,rlat,qlatu,qlonu,dummy,ug,vg
real :: qlatv,qlonv,cosuv,qlatuv,qlonuv,sinuv,uvgx,uvgy,uvgz,uvgr,raxis,raxisi

real :: r_interp(22)

! Determine index of lowest ZONAVG pressure level that is at least 1/2
! ZONAVG pressure level higher than highest input pressure data level 
! (i.e., maximum zonp_vect value that is less than 82.5% of levpr(nprz),
! which is in hPa)

lzon_bot = nint(19. - 6. * alog10(real(levpr(nprz)))) + 1

! Fill column array of pressure level pressure values in MKS

do k = 1,nprz
   pcol_p(k+2) = levpr(k) * 100.
enddo
pcol_p(2) = 110000.   ! Phony underground level
pcol_p(1) = 120000.   ! Phony underground level

! Fill upper ppd levels from ZONAVG data

k = nprz + 2                ! highest ppd level filled so far
kzonoff = k + 1 - lzon_bot  ! k offset for copying pzon levels to ppd
npd = 22 + kzonoff          ! Total number of pressure levels

do levp = lzon_bot,22
   k = levp + kzonoff
   pcol_p(k) = zonp_vect(levp)
enddo

! Fill column array of exner function values

do k = 1, npd
   pcol_exner(k) = cp * (pcol_p(k) * p00i)**rocp
enddo

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

! fractional x/y indices in pressure data arrays at current iw point location 

   gry = (glatw(iw) - xswlat) / gdatdy + 3.
   grx = (glonw(iw) - xswlon) / gdatdx + 1. + real(ipoffset) 

! Horizontally interpolate gridded pressure-level data to column
! at location of current W point

   do k = 1,nprz
      call gdtost(p_t(:,:,k),nprx+4,npry+4,grx,gry,pcol_temp(k+2)) ! temp
      call gdtost(p_z(:,:,k),nprx+4,npry+4,grx,gry,pcol_z(k+2))    ! geop ht
      call gdtost(p_u(:,:,k),nprx+4,npry+4,grx,gry,pcol_u(k+2))    ! uzonal
      call gdtost(p_v(:,:,k),nprx+4,npry+4,grx,gry,pcol_v(k+2))    ! umerid
   enddo

   do k=1,nprz_rh
      call gdtost(p_r(:,:,k),nprx+4,npry+4,grx,gry,pcol_rt(k+2))    ! s.h.
   enddo

! Horizontally interpolate surface gridded data to column
! at location of current W point

   if (ihydsfc == 1) then
      call gdtost(p_topo (:,:),nprx+4,npry+4,grx,gry,pcol_topo )
      call gdtost(p_prsfc(:,:),nprx+4,npry+4,grx,gry,pcol_prsfc)
      call gdtost(p_tsfc (:,:),nprx+4,npry+4,grx,gry,pcol_tsfc )
      call gdtost(p_shsfc(:,:),nprx+4,npry+4,grx,gry,pcol_shsfc)
   endif
 
   pcol_prsfc = pcol_prsfc * 100.
   pcol_exnersfc = cp * (pcol_prsfc * p00i)**rocp

! If water vapor stops at a lower level, fill pcol_rt from levels nprz_rh+3
! to nprz+2 with values interpolated from the ZONAVG data

   rlat = .4 * (glatw(iw) + 93.75)
   ilat = int(rlat)
   wt2 = rlat - float(ilat)

   if (nprz_rh < nprz) then
      nlevs = nprz - nprz_rh
      kstrt = nprz_rh + 3
      r_interp(:) = (1. - wt2) * zonr(ilat,:) + wt2 * zonr(ilat+1,:)
      call pintrp_ee(22, r_interp, zonp_vect, nlevs, pcol_rt(kstrt), pcol_p(kstrt))
   endif

! Linearly interpolate zonavg arrays by latitude to current IW column 
! and K level. WIND IS ALL U AND NO V.

   do levp = lzon_bot,22
      k = levp + kzonoff
      pcol_temp(k) = (1. - wt2) * zont(ilat,levp) + wt2 * zont(ilat+1,levp)
      pcol_z(k)    = (1. - wt2) * zonz(ilat,levp) + wt2 * zonz(ilat+1,levp)
      pcol_rt(k)   = (1. - wt2) * zonr(ilat,levp) + wt2 * zonr(ilat+1,levp)
      pcol_u(k)    = (1. - wt2) * zonu(ilat,levp) + wt2 * zonu(ilat+1,levp)
      pcol_v(k)    = 0.
   enddo

! Fill phony underground values for PPD levels 1100 mb and 1200 mb

   pcol_u(1:2)  = pcol_u(3)
   pcol_v(1:2)  = pcol_v(3)
   pcol_rt(1:2) = pcol_rt(3)
   pcol_temp(2) = pcol_temp(3) + 5.3  ! Uses approx std lapse rate
   pcol_temp(1) = pcol_temp(2) + 4.9  ! Uses approx std lapse rate

! Vertically interpolate current column to model grid and 
! perform iterative hydrostatic balance 

   call vterpp_s(iw,o_rho,o_theta,o_shv,o_uzonal,o_umerid, &
                 pcol_p, pcol_temp, pcol_z, pcol_u, pcol_v, &
                 pcol_rt, pcol_r, pcol_exner, &
                 pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc)

enddo

if (iparallel == 1) then
   mrl = 1
   call mpi_send_w(mrl, rvara1=o_uzonal, rvara2=o_umerid)
   call mpi_recv_w(mrl, rvara1=o_uzonal, rvara2=o_umerid)
endif

! Initialize V wind component

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   if (iw1 < 2) iw1 = iw2
   if (iw2 < 2) iw2 = iw1

   ka = lpv(iv)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

! Average winds to V point and rotate at V point

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = ka,mza
         ug = .5 * (o_uzonal(k,iw1) + o_uzonal(k,iw2))
         vg = .5 * (o_umerid(k,iw1) + o_umerid(k,iw2))

         uvgr = -vg * zev(iv) * eradi  ! radially outward from axis

         uvgx = (-ug * yev(iv) + uvgr * xev(iv)) * raxisi 
         uvgy = ( ug * xev(iv) + uvgr * yev(iv)) * raxisi 
         uvgz =   vg * raxis * eradi 

         o_vc(k,iv) = uvgx * vnx(iv) + uvgy * vny(iv) + uvgz * vnz(iv)
      enddo

   else
      o_vc(:,iv) = 0.
   endif

enddo

! MPI parallel send/recv of o_vc

if (iparallel == 1) then
   mrl = 1
   call mpi_send_v(mrl, rvara1=o_vc)
   call mpi_recv_v(mrl, rvara1=o_vc)
endif

! LBC copy of o_uvc

call lbcopy_v(1, vc=o_vc)

end subroutine isnstage

!===============================================================================

subroutine vterpp_s(iw,o_rho,o_theta,o_shv,o_uzonal,o_umerid, &
                 pcol_p, pcol_temp, pcol_z, pcol_u, pcol_v, &
                 pcol_rt, pcol_r, pcol_exner, &
                 pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc)

use max_dims,    only: maxpr
use isan_coms,   only: npd, ihydsfc
use consts_coms, only: r8, grav2, grav, cvocp, p00k, rdry, rvap, p00, &
                       cp, rocp, gravi, eps_virt, eps_vap
use mem_grid,    only: mwa, mza, lpw, dzt_top, dzt_bot, zt
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw

real(r8), intent(inout) :: o_rho   (mza,mwa)
real,     intent(inout) :: o_theta (mza,mwa)
real,     intent(inout) :: o_shv   (mza,mwa)
real,     intent(inout) :: o_uzonal(mza,mwa)
real,     intent(inout) :: o_umerid(mza,mwa)

real(r8) :: o_press(mza)  ! automatic array
real(r8) :: pressnew, pkhyd

real :: vctr1(mza)  ! automatic array
real :: vctr2(mza)  ! automatic array
real :: vctr3(mza)  ! automatic array
real :: vctr4(mza)  ! automatic array
real :: vctr5(mza)  ! automatic array

real, intent(in)    :: pcol_p    (maxpr+2)
real, intent(in)    :: pcol_temp (maxpr+2)
real, intent(inout) :: pcol_z    (maxpr+2)
real, intent(in)    :: pcol_u    (maxpr+2)
real, intent(in)    :: pcol_v    (maxpr+2)
real, intent(inout) :: pcol_rt   (maxpr+2)
real, intent(in)    :: pcol_r    (maxpr+2)
real, intent(in)    :: pcol_exner(maxpr+2)

real, intent(in) :: pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc

real :: pcol_thet (maxpr+2)
real :: pcol_thv  (maxpr+2)

real :: pcol_thetsfc, pcol_thvsfc

integer :: i,j,k,mcnt,kl,kpbc,klo,khi,kbc,levp,kother,iter,ka
real :: cpo2g, z, extrap, vapor_press
real, external :: eslf

kpbc = 0

do k = 1,npd

!! NOTE:
!! HUMIDITY ALREADY CONVERTED TO SPECIFIC HUMIDITY

!!! Compute ambient vapor pressure based on relative humidity (pcol_r)
!!! and saturation vapor pressure (eslf)
!!
!!   vapor_press = pcol_r(k) * eslf(pcol_temp(k)-273.15)
!!
!!! Do not allow vapor pressure to exceed ambient pressure
!!
!!   vapor_press = min(pcol_p(k),vapor_press)
!!
!!! Compute specific humidity (pcol_rt) from vapor pressure and ambient pressure
!!
!!   pcol_rt(k) = eps_vap * vapor_press  &
!!              / (pcol_p(k) + vapor_press * (eps_vap - 1.))

! Compute potential temperatures

   pcol_thet(k) = pcol_temp(k) * cp / pcol_exner(k)                ! dry theta
   pcol_thv (k) = pcol_thet(k) * (1. + eps_virt * pcol_rt(k))   ! virtual theta

! Find highest k level for which pressure is greater than surface pressure

   if (ihydsfc == 1 .and. pcol_p(k) > pcol_prsfc) kpbc = k

enddo

! If geopotential heights of pressure levels are not available, compute them by
! hydrostatic integration using topography and surface pressure (if available)
! as the boundary condition.

if (ihydsfc == 1) then

   pcol_thetsfc = pcol_tsfc * cp / pcol_exnersfc                ! dry theta
   pcol_thvsfc  = pcol_thetsfc * (1. + eps_virt * pcol_shsfc)   ! virtual theta

   pcol_z(kpbc) = pcol_topo + .5 * (pcol_thvsfc + pcol_thv(kpbc)) &
             * (pcol_exnersfc - pcol_exner(kpbc)) * gravi

   pcol_z(kpbc+1) = pcol_topo + .5 * (pcol_thvsfc + pcol_thv(kpbc+1)) &
             * (pcol_exnersfc - pcol_exner(kpbc+1)) * gravi

   do k = kpbc+2,npd
      pcol_z(k) = pcol_z(k-1) + .5 * (pcol_thv(k-1) + pcol_thv(k)) &
                * (pcol_exner(k-1) - pcol_exner(k)) * gravi
   enddo

   do k = kpbc-1,1,-1
      pcol_z(k) = pcol_z(k+1) + .5 * (pcol_thv(k+1) + pcol_thv(k)) &
                * (pcol_exner(k+1) - pcol_exner(k)) * gravi
   enddo

else

! If geopotential heights of pressure levels ARE available, use hydrostatic
! integration only downward from k = 3 level to get heights of phony
! underground pressure levels

   pcol_z(2) = pcol_z(3) + .5 * (pcol_thv(3) + pcol_thv(2)) &
             * (pcol_exner(3) - pcol_exner(2)) * gravi

   pcol_z(1) = pcol_z(2) + .5 * (pcol_thv(2) + pcol_thv(1)) &
             * (pcol_exner(2) - pcol_exner(1)) * gravi

endif

!!cpo2g = cp / grav2
!!
!!do levp = lzon_bot,22
!!   k = levp + kzonoff
!!   pcol_z(k) = pcol_z(k-1) + cpo2g * (pcol_exner(k-1) - pcol_exner(k))   &
!!      * (pcol_thv(k-1) + pcol_thv(k))
!!enddo

! Now all "pcol" arrays are filled; vertically interpolate them to model grid

call hintrp_cc(npd,pcol_p   ,pcol_z,mza,vctr1,zt)  ! pressure
call hintrp_cc(npd,pcol_thet,pcol_z,mza,vctr2,zt)  ! theta
call hintrp_cc(npd,pcol_rt  ,pcol_z,mza,vctr3,zt)  ! vapor specific humidity
call hintrp_cc(npd,pcol_u   ,pcol_z,mza,vctr4,zt)  ! zonal wind
call hintrp_cc(npd,pcol_v   ,pcol_z,mza,vctr5,zt)  ! merid wind

ka = lpw(iw)

do k = ka,mza
   o_press (k)    = vctr1(k)
   o_theta (k,iw) = vctr2(k)
   o_shv   (k,iw) = max(1.e-8,vctr3(k))
   o_uzonal(k,iw) = vctr4(k)
   o_umerid(k,iw) = vctr5(k)
enddo

! Choose as an internal pressure boundary condition the pcol_p pressure level 
! at or below (in elevation) the 49900 Pa surface.  Find the k index of this level.

kpbc = npd
do while (pcol_p(kpbc) < 49900.)
   kpbc = kpbc - 1
enddo

! Determine which two model zt levels bracket pcol_z(kpbc) in this column

khi = 2
do while (zt(khi) < pcol_z(kpbc))
   khi = khi + 1
enddo
klo = khi - 1

! Determine whether zt(klo) or zt(khi) is closer to pcol_z(kpbc).  The closer
! one will undergo direct pressure adjustment to satisfy the internal b.c.

if (zt(khi) - pcol_z(kpbc) < pcol_z(kpbc) - zt(klo)) then
   kbc = khi
   kother = klo
else
   kbc = klo
   kother = khi
endif

extrap = (zt(kbc) - zt(kother)) / (pcol_z(kpbc) - zt(kother))

! Carry out iterative hydrostatic balance procedure

do iter = 1,100

! Adjust pressure at k = kbc.  Use temporal weighting for damping

   pressnew = o_press(kother) * (pcol_p(kpbc) / o_press(kother)) ** extrap
   o_press(kbc) = .1 * o_press(kbc) + .9 * pressnew

! Compute density for all levels (assume micphys level = 1)

   do k = ka,mza
      o_rho(k,iw) = o_press(k) ** cvocp * p00k &
         / (o_theta(k,iw) * (rdry * (1. - o_shv(k,iw)) + rvap * o_shv(k,iw)))
   enddo

! Integrate hydrostatic equation upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations

   do k = kbc+1,mza
      pkhyd = o_press(k-1) &
            - grav * (o_rho(k-1,iw) * dzt_top(k-1) + o_rho(k,iw) * dzt_bot(k))
      o_press(k) = .05 * o_press(k) + .95 * max(.1,pkhyd)
   enddo

   do k = kbc-1,ka,-1
      pkhyd = o_press(k+1) &
            + grav * (o_rho(k+1,iw) * dzt_bot(k+1) + o_rho(k,iw) * dzt_top(k))
      o_press(k) = .05 * o_press(k) + .95 * max(.1,pkhyd)
   enddo

enddo

end subroutine vterpp_s
