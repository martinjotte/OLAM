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
subroutine isnstage(p_u,p_v,p_t,p_z,p_r, &
                    o_rho, o_theta, o_shv, o_uzonal, o_umerid, o_uvc)

use max_dims,   only: maxpr
use isan_coms,  only: nprz, npry, nprx, nprz_rh, pcol_v, &
                      pcol_u, pcol_rt, pcol_z, pcol_temp, gdatdx, gdatdy, &
                      npd, kzonoff, levpr, lzon_bot, pcol_p
use mem_grid,   only: glatw, glonw, mza, mwa, mva, &
                      xeu, yeu, zeu, xev, yev, zev, &
                      unx, uny, unz, vnx, vny, vnz
use mem_ijtabs, only: jtab_u, jtab_v, jtab_w, itab_u, itab_v, itab_w
use mem_zonavg, only: zonp_vect, zont, zonz, zonr, zonu
use consts_coms, only: eradi
use misc_coms,  only: io6, meshtype

implicit none

real, intent(in) :: p_u(nprx+3,npry+2,nprz)
real, intent(in) :: p_v(nprx+3,npry+2,nprz)
real, intent(in) :: p_t(nprx+3,npry+2,nprz)
real, intent(in) :: p_z(nprx+3,npry+2,nprz)
real, intent(in) :: p_r(nprx+3,npry+2,nprz)

real(kind=8), intent(out) :: o_rho   (mza,mwa)
real,         intent(out) :: o_theta (mza,mwa)
real,         intent(out) :: o_shv   (mza,mwa)
real,         intent(out) :: o_uzonal(mza,mwa)
real,         intent(out) :: o_umerid(mza,mwa)
real,         intent(out) :: o_uvc   (mza,mva)

character(3) :: csuff
character(8) :: rot_type

integer :: ngrd,lv,lf,k,levp,j,iw,ilat,iu,iv,iuv  &
   ,iw1,iw2,nlevs,kstrt

real :: wt2,grx,gry,rlat,qlatu,qlonu,dummy,ug,vg  &
   ,qlatv,qlonv,cosuv,qlatuv,qlonuv,sinuv,uvgx,uvgy,uvgz,uvgr,raxis,raxisi

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

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(7)%jend(1); iw = jtab_w(7)%iw(j)
!---------------------------------------------------------------------
call qsub('W',iw)

! fractional x/y indices in pressure data arrays at current iw point location 

   gry = (glatw(iw) +  90.) / gdatdy + 2.
   grx = (glonw(iw) + 180.) / gdatdx + 2.

! Horizontally interpolate gridded pressure-level data to column
! at location of current W point

   do k = 1,nprz
      call gdtost(p_t(1,1,k),nprx+3,npry+2,grx,gry,pcol_temp(k+2)) ! temp
      call gdtost(p_z(1,1,k),nprx+3,npry+2,grx,gry,pcol_z(k+2))    ! geop ht
      call gdtost(p_u(1,1,k),nprx+3,npry+2,grx,gry,pcol_u(k+2))    ! uzonal
      call gdtost(p_v(1,1,k),nprx+3,npry+2,grx,gry,pcol_v(k+2))    ! umerid
   enddo

   do k=1,nprz_rh
      call gdtost(p_r(1,1,k),nprx+3,npry+2,grx,gry,pcol_rt(k+2))    ! s.h.
   enddo

! surface gridded data

!   call gdtost(p_sft(1,1),nprx+3,npry+2,grx,gry,sfc_temp(iw))
!   call gdtost(p_sst(1,1),nprx+3,npry+2,grx,gry,seatp(iw))
!   call gdtost(p_snow(1,1),nprx+3,npry+2,grx,gry,snow(iw))
 
   rlat = .4 * (glatw(iw) + 93.75)
   ilat = int(rlat)
   wt2 = rlat - float(ilat)

! If water vapor stops at a lower level, fill pcol_rt from levels nprz_rh+3
! to nprz+2 with values interpolated from the ZONAVG data

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

! Vertically interpolate current column to model grid and 
! perform iterative hydrostatic balance 

   call vterpp_s(iw,o_rho,o_theta,o_shv,o_uzonal,o_umerid)

enddo
call rsub('Wa',7)

if (meshtype == 1) then

! If triangular mesh, initialize U wind component

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(7)%jend(1); iu = jtab_u(7)%iu(j)
      iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      if (iw1 < 2) iw1 = iw2
      if (iw2 < 2) iw2 = iw1

      raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

! Average winds to U point and rotate at U point

      if (raxis > 1.e3) then
         raxisi = 1. / raxis

         do k = 1,mza
            ug = .5 * (o_uzonal(k,iw1) + o_uzonal(k,iw2))
            vg = .5 * (o_umerid(k,iw1) + o_umerid(k,iw2))

            uvgr = -vg * zeu(iu) * eradi  ! radially outward from axis

            uvgx = (-ug * yeu(iu) + uvgr * xeu(iu)) * raxisi
            uvgy = ( ug * xeu(iu) + uvgr * yeu(iu)) * raxisi 
            uvgz =   vg * raxis * eradi 

            o_uvc(k,iu) = uvgx * unx(iu) + uvgy * uny(iu) + uvgz * unz(iu)
         enddo
      else
         o_uvc(:,iu) = 0.
      endif

   enddo
   call rsub('Ua',7)

else

! If using hexagonal mesh, initialize V wind component

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_v(7)%jend(1); iv = jtab_v(7)%iv(j)
      iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   call qsub('V',iv)

      if (iw1 < 2) iw1 = iw2
      if (iw2 < 2) iw2 = iw1

      raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

! Average winds to V point and rotate at V point

      if (raxis > 1.e3) then
         raxisi = 1. / raxis

         do k = 1,mza
            ug = .5 * (o_uzonal(k,iw1) + o_uzonal(k,iw2))
            vg = .5 * (o_umerid(k,iw1) + o_umerid(k,iw2))

            uvgr = -vg * zev(iv) * eradi  ! radially outward from axis

            uvgx = (-ug * yev(iv) + uvgr * xev(iv)) * raxisi 
            uvgy = ( ug * xev(iv) + uvgr * yev(iv)) * raxisi 
            uvgz =   vg * raxis * eradi 

            o_uvc(k,iv) = uvgx * vnx(iv) + uvgy * vny(iv) + uvgz * vnz(iv)
         enddo

      else
         o_uvc(:,iv) = 0.
      endif

   enddo
   call rsub('Va',7)

endif

end subroutine isnstage

!===============================================================================

subroutine vterpp_s(iw,o_rho,o_theta,o_shv,o_uzonal,o_umerid)

use isan_coms,   only: pcol_z, pcol_thv, pcol_rt, pcol_thet, pcol_pi,  &
                       pcol_p, pcol_pk, npd, pcol_temp, pcol_r, pcol_v,  &
                       pcol_u
use consts_coms, only: grav2, gravo2, cvocp, p00k, rdry, rvap, p00,  &
                       cp, rocp, grav, eps_virt, eps_vap
use mem_grid,    only: mwa, mza, dzt, zt
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw

real(kind=8), intent(out) :: o_rho   (mza,mwa)
real,         intent(out) :: o_theta (mza,mwa)
real,         intent(out) :: o_shv   (mza,mwa)
real,         intent(out) :: o_uzonal(mza,mwa)
real,         intent(out) :: o_umerid(mza,mwa)

real(kind=8) :: o_press(mza)  ! automatic array

real :: vctr1(mza)  ! automatic array
real :: vctr2(mza)  ! automatic array
real :: vctr3(mza)  ! automatic array
real :: vctr4(mza)  ! automatic array
real :: vctr5(mza)  ! automatic array

integer :: i,j,k,mcnt,kl,kpbc,klo,khi,kbc,levp,kother,iter
real :: pbc,cpo2g,piocp,z,extrap,pressnew,pkhyd, vapor_press
real, external :: eslf

! Fill phony underground values for PPD levels 1100 mb and 1200 mb

pcol_u(1:2)  = pcol_u(3)
pcol_v(1:2)  = pcol_v(3)
pcol_rt(1:2) = pcol_rt(3)
pcol_temp(2) = pcol_temp(3) + 5.3  ! Uses approx std lapse rate
pcol_temp(1) = pcol_temp(2) + 4.9  ! Uses approx std lapse rate

do k = 1,npd
   pcol_pk(k) = pcol_p(k)**rocp
   piocp = (pcol_p(k)/p00)**rocp
   pcol_pi(k) = cp * piocp
   pcol_thet(k) = pcol_temp(k) / piocp                  ! dry theta

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

! Compute virtual potential temperature

   pcol_thv(k) = pcol_thet(k) * (1. + eps_virt * pcol_rt(k))
enddo

! Use hydrostatic integration to get heights of phony underground press levels

pcol_z(2) = pcol_z(3) + .5 * (pcol_thv(3) + pcol_thv(2))   &
          * (pcol_pi(3) - pcol_pi(2)) / grav

pcol_z(1) = pcol_z(2) + .5 * (pcol_thv(2) + pcol_thv(1))   &
          * (pcol_pi(2) - pcol_pi(1)) / grav

cpo2g = cp / grav2

!!do levp = lzon_bot,22
!!   k = levp + kzonoff
!!   pcol_z(k) = pcol_z(k-1) + cpo2g * (pcol_pi(k-1) - pcol_pi(k))   &
!!      * (pcol_thv(k-1) + pcol_thv(k))
!!enddo

! Now all "pcol" arrays are filled; vertically interpolate them to model grid

call hintrp_cc(npd,pcol_p   ,pcol_z,mza,vctr1,zt)  ! pressure
call hintrp_cc(npd,pcol_thet,pcol_z,mza,vctr2,zt)  ! theta
call hintrp_cc(npd,pcol_rt  ,pcol_z,mza,vctr3,zt)  ! vapor specific humidity
call hintrp_cc(npd,pcol_u   ,pcol_z,mza,vctr4,zt)  ! zonal wind
call hintrp_cc(npd,pcol_v   ,pcol_z,mza,vctr5,zt)  ! merid wind

do k = 1,mza
   o_press (k)    = vctr1(k)
   o_theta (k,iw) = vctr2(k)
   o_shv   (k,iw) = max(1.e-8,vctr3(k))
   o_uzonal(k,iw) = vctr4(k)
   o_umerid(k,iw) = vctr5(k)
enddo

! Hydrostatically balance fields on model grid

! Choose as an internal pressure boundary condition the pcol_p pressure level 
! at or below (in elevation) the 599 mb surface.  Find the k index of this level.

pbc = 59900.
kpbc = npd
do while (pcol_p(kpbc) < pbc)
   kpbc = kpbc - 1
enddo

! Determine which two model zt levels bracket pcol_z(kpbc) in this column

khi = 2
do while (zt(khi) < pcol_z(kpbc))
   khi = khi + 1
enddo
klo = khi - 1

! Determine whether zt(klo) or zt(khi) is closer to zd(kpbc).  The closer
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

! Adjust pressure at k = kbc.  Use weighting for damping

   pressnew = o_press(kother) * (pcol_p(kpbc) / o_press(kother)) ** extrap
   o_press(kbc) = .1 * o_press(kbc) + .9 * pressnew

!  Compute density for all levels

   do k = 1,mza
      o_rho(k,iw) = o_press(k) ** cvocp * p00k  &
         / (o_theta(k,iw) * (rdry * (1. - o_shv(k,iw)) + rvap * o_shv(k,iw)))
   enddo

! Hydrostatically integrate upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations

   do k = kbc+1,mza
      pkhyd = o_press(k-1)  &
         - gravo2 * (o_rho(k-1,iw) * dzt(k-1) + o_rho(k,iw) * dzt(k))
      o_press(k) = .05 * o_press(k) + .95 * max(.1,pkhyd)
   enddo

   do k = kbc-1,1,-1
      pkhyd = o_press(k+1)  &
         + gravo2 * (o_rho(k+1,iw) * dzt(k+1) + o_rho(k,iw) * dzt(k))
      o_press(k) = .05 * o_press(k) + .95 * max(.1,pkhyd)
   enddo

enddo

end subroutine vterpp_s
