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
subroutine isnstage(p_u, p_v, p_t, p_z, p_r, p_o, &
                    p_topo, p_prsfc, p_tsfc, p_shsfc, &
                    o_press, o_rho, o_theta, o_rrw, o_uzonal, o_umerid, o_ozone)

use max_dims,     only: maxpr
use isan_coms,    only: nprz, npry, nprx, nprz_rh, haso3, nbot_o3, &
                        xswlat, xswlon, gdatdx, gdatdy, glat, &
                        npd, kzonoff, levpr, lzon_bot, ipoffset, inproj
use mem_grid,     only: glatw, glonw, mza, mwa
use mem_ijtabs,   only: jtab_w, jtw_init
use mem_zonavg,   only: zonp_vect, zont, zonz, zonr, zonu, zono
use consts_coms,  only: r8, rocp, p00i, cp
use misc_coms,    only: iparallel
use isan_coms,    only: ihydsfc
use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

implicit none

real, intent(in) :: p_u(nprx+4,npry+4,nprz)
real, intent(in) :: p_v(nprx+4,npry+4,nprz)
real, intent(in) :: p_t(nprx+4,npry+4,nprz)
real, intent(in) :: p_z(nprx+4,npry+4,nprz)
real, intent(in) :: p_r(nprx+4,npry+4,nprz)
real, intent(in) :: p_o(nprx+4,npry+4,nprz)

real, intent(in) :: p_topo (nprx+4,npry+4)
real, intent(in) :: p_prsfc(nprx+4,npry+4)
real, intent(in) :: p_tsfc (nprx+4,npry+4)
real, intent(in) :: p_shsfc(nprx+4,npry+4)

real(r8), intent(inout) :: o_press (mza,mwa)
real(r8), intent(inout) :: o_rho   (mza,mwa)
real,     intent(inout) :: o_theta (mza,mwa)
real,     intent(inout) :: o_rrw   (mza,mwa)
real,     intent(inout) :: o_uzonal(mza,mwa)
real,     intent(inout) :: o_umerid(mza,mwa)
real,     intent(inout) :: o_ozone (mza,mwa)

real :: pcol_p    (maxpr+2)
real :: pcol_temp (maxpr+2)
real :: pcol_z    (maxpr+2)
real :: pcol_u    (maxpr+2)
real :: pcol_v    (maxpr+2)
real :: pcol_rt   (maxpr+2)
real :: pcol_exner(maxpr+2)
real :: pcol_o3   (maxpr+2)

real :: plat(npry+4)

real :: pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc

integer :: k,levp,j,iw,ilat,mrl
integer :: nlevs,kstrt

real :: wt2,grx,gry,rlat
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

! For INPROJ = 2, where latitudes are specified, fill expanded latitude array

if (inproj == 2) then

   do ilat = 1,npry
      plat(ilat+2) = glat(ilat)
   enddo

   plat(2) = plat(3) - (plat(4) - plat(3))
   plat(1) = plat(2) - (plat(3) - plat(2))

   plat(npry+3) = plat(npry+2) + (plat(npry+2) - plat(npry+1))
   plat(npry+4) = plat(npry+3) + (plat(npry+3) - plat(npry+2))

endif

!----------------------------------------------------------------------
!$omp parallel do private(iw,gry,grx,ilat,k,pcol_temp,pcol_z,pcol_u,pcol_v,&
!$omp                     pcol_rt,pcol_o3,pcol_topo,pcol_prsfc,pcol_tsfc,&
!$omp                     pcol_shsfc,pcol_exnersfc,rlat,wt2,nlevs,kstrt,&
!$omp                     r_interp,levp)
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

! fractional x/y indices in pressure data arrays at current iw point location

   if (inproj == 1) then

      gry = (glatw(iw) - xswlat) / gdatdy + 3.
      grx = (glonw(iw) - xswlon) / gdatdx + 1. + real(ipoffset)

   elseif (inproj == 2) then

      ! estimate latitude index assuming uniform spacing of plat

      ilat = 2 + npry * int((glatw(iw) - plat(2)) / 180.)

      ! find correct latitude index

      if (plat(ilat) > glatw(iw)) then
         do while(plat(ilat) > glatw(iw))
            ilat = ilat - 1
         enddo
      elseif (plat(ilat+1) < glatw(iw)) then
         do while(plat(ilat+1) < glatw(iw))
            ilat = ilat + 1
         enddo
      endif

      gry = (glatw(iw) - plat(ilat)) / (plat(ilat+1) - plat(ilat)) + real(ilat)
      grx = (glonw(iw) - xswlon) / gdatdx + 1. + real(ipoffset)

   endif

! Horizontally interpolate gridded pressure-level data to column
! at location of current W point

   do k = 1,nprz
      call gdtost(p_t(:,:,k),nprx+4,npry+4,grx,gry,pcol_temp(k+2)) ! temp
      call gdtost(p_z(:,:,k),nprx+4,npry+4,grx,gry,pcol_z(k+2))    ! geop ht
      call gdtost(p_u(:,:,k),nprx+4,npry+4,grx,gry,pcol_u(k+2))    ! uzonal
      call gdtost(p_v(:,:,k),nprx+4,npry+4,grx,gry,pcol_v(k+2))    ! umerid
   enddo

   do k=1,nprz_rh
      call gdtost(p_r(:,:,k),nprx+4,npry+4,grx,gry,pcol_rt(k+2))   ! mix. ratio
   enddo

   if (haso3) then
      do k = nbot_o3, nprz
         call gdtost(p_o(:,:,k),nprx+4,npry+4,grx,gry,pcol_o3(k+2)) ! ozone
      enddo
   endif

! Horizontally interpolate surface gridded data to column
! at location of current W point

   if (ihydsfc == 1) then
      call gdtost(p_topo (:,:),nprx+4,npry+4,grx,gry,pcol_topo )
      call gdtost(p_prsfc(:,:),nprx+4,npry+4,grx,gry,pcol_prsfc)
      call gdtost(p_tsfc (:,:),nprx+4,npry+4,grx,gry,pcol_tsfc )
      call gdtost(p_shsfc(:,:),nprx+4,npry+4,grx,gry,pcol_shsfc)

      pcol_prsfc = pcol_prsfc * 100.
      pcol_exnersfc = cp * (pcol_prsfc * p00i)**rocp
   endif

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

   if (haso3) then

      if (nbot_o3 /= 1) then

         ! If ozone is not reported at lower levels, fill pcol_o3 from levels
         ! 3 to nbot_o3 + 1 with values interpolated from the ZONAVG data

         nlevs = nbot_o3 - 1
         kstrt = 3
         r_interp(:) = (1. - wt2) * zono(ilat,:) + wt2 * zono(ilat+1,:)
         call pintrp_ee(22, r_interp, zonp_vect, nlevs, pcol_o3(kstrt), pcol_p(kstrt))
      endif

   else

      ! If ozone is not reported at all, fill pcol_o3 from levels 3
      ! to nbrz+2 with values interpolated from the ZONAVG data

      nlevs = nprz
      kstrt = 3
      r_interp(:) = (1. - wt2) * zono(ilat,:) + wt2 * zono(ilat+1,:)
      call pintrp_ee(22, r_interp, zonp_vect, nlevs, pcol_o3(kstrt), pcol_p(kstrt))

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
      pcol_o3(k)   = (1. - wt2) * zono(ilat,levp) + wt2 * zono(ilat+1,levp)
   enddo

! Fill phony underground values for PPD levels 1100 mb and 1200 mb

   pcol_u(1:2)  = pcol_u(3)
   pcol_v(1:2)  = pcol_v(3)
   pcol_rt(1:2) = pcol_rt(3)
   pcol_temp(2) = pcol_temp(3) + 5.3  ! Uses approx std lapse rate
   pcol_temp(1) = pcol_temp(2) + 4.9  ! Uses approx std lapse rate
   pcol_o3(1:2) = pcol_o3(3)

! Vertically interpolate current column to model grid and
! perform iterative hydrostatic balance

   call vterpp_s(iw,o_press,o_rho,o_theta,o_rrw,o_uzonal,o_umerid,o_ozone, &
                 pcol_p, pcol_temp, pcol_z, pcol_u, pcol_v, &
                 pcol_rt, pcol_exner, pcol_o3, &
                 pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc)

enddo
!$omp end parallel do

if (iparallel == 1) then
   mrl = 1
   call mpi_send_w(mrl, rvara1=o_uzonal, rvara2=o_umerid)
   call mpi_recv_w(mrl, rvara1=o_uzonal, rvara2=o_umerid)
endif

end subroutine isnstage

!===============================================================================

subroutine vterpp_s(iw,o_press,o_rho,o_theta,o_rrw,o_uzonal,o_umerid,o_ozone, &
                 pcol_p, pcol_temp, pcol_z, pcol_u, pcol_v, &
                 pcol_rt, pcol_exner, pcol_o3, &
                 pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc)

use max_dims,    only: maxpr
use isan_coms,   only: npd, ihydsfc
use consts_coms, only: r8, grav2, grav, cvocp, p00kord, rvap, p00, p00i, &
                       t00, cp, rocp, gravi, eps_virt, eps_vapi
use mem_grid,    only: mwa, mza, lpw, zt, gdz_belo8, gdz_abov8
use micro_coms,  only: miclevel
use therm_lib,   only: rhovsl

implicit none

integer,  intent(in)    :: iw

real(r8), intent(inout) :: o_rho   (mza,mwa)
real(r8), intent(inout) :: o_press (mza,mwa)
real,     intent(inout) :: o_theta (mza,mwa)
real,     intent(inout) :: o_rrw   (mza,mwa)
real,     intent(inout) :: o_uzonal(mza,mwa)
real,     intent(inout) :: o_umerid(mza,mwa)
real,     intent(inout) :: o_ozone (mza,mwa)

real, intent(in)    :: pcol_p    (maxpr+2)
real, intent(in)    :: pcol_temp (maxpr+2)
real, intent(inout) :: pcol_z    (maxpr+2)
real, intent(in)    :: pcol_u    (maxpr+2)
real, intent(in)    :: pcol_v    (maxpr+2)
real, intent(inout) :: pcol_rt   (maxpr+2)
real, intent(in)    :: pcol_exner(maxpr+2)
real, intent(in)    :: pcol_o3   (maxpr+2)

real, intent(in) :: pcol_topo, pcol_prsfc, pcol_tsfc, pcol_shsfc, pcol_exnersfc

real(r8) :: rho_tot(mza)  ! automatic array
real(r8) :: pressnew, pkhyd

real :: vctr1(mza)  ! automatic array
real :: vctr2(mza)  ! automatic array
real :: vctr3(mza)  ! automatic array
real :: vctr4(mza)  ! automatic array
real :: vctr5(mza)  ! automatic array
real :: vctr6(mza)  ! automatic array

real, parameter :: mwair  = 28.9628             ! molecular weight of air
real, parameter :: mwo3   = 48.0                ! molecular weight of ozone
real, parameter :: cnvto3 = mwair / mwo3 * 1.e6 ! ozone mixing ratio to ppmV

real :: pcol_thet(maxpr+2)
real :: pcol_thv (maxpr+2)

real :: pcol_thetsfc, pcol_thvsfc

integer :: k,kpbc,klo,khi,kbc,kother,iter,ka
real    :: extrap, exner, tairc, cond, rrv

kpbc = 0

do k = 1,npd

! Compute potential temperatures

   pcol_thet(k) = pcol_temp(k) * cp / pcol_exner(k)             ! dry theta
   pcol_thv (k) = pcol_thet(k) * (1. + eps_virt * pcol_rt(k))   ! virtual theta

! Find highest k level for which pressure is greater than surface pressure

   if (ihydsfc == 1) then
      if (pcol_p(k) > pcol_prsfc) kpbc = k
   endif

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

! Now all "pcol" arrays are filled; vertically interpolate them to model grid

call hintrp_cc(npd,pcol_p   ,pcol_z,mza,vctr1,zt)  ! pressure
call hintrp_cc(npd,pcol_thet,pcol_z,mza,vctr2,zt)  ! theta
call hintrp_cc(npd,pcol_rt  ,pcol_z,mza,vctr3,zt)  ! water mixing ratio
call hintrp_cc(npd,pcol_u   ,pcol_z,mza,vctr4,zt)  ! zonal wind
call hintrp_cc(npd,pcol_v   ,pcol_z,mza,vctr5,zt)  ! merid wind
call hintrp_cc(npd,pcol_o3  ,pcol_z,mza,vctr6,zt)  ! ozone

ka = lpw(iw)

do k = ka,mza
   o_press (k,iw) = vctr1(k)
   o_theta (k,iw) = vctr2(k)
   o_rrw   (k,iw) = max(1.e-8,vctr3(k))
   o_uzonal(k,iw) = vctr4(k)
   o_umerid(k,iw) = vctr5(k)
   o_ozone (k,iw) = max( 1.e-30, vctr6(k)*cnvto3) ! mix ratio to ppmV
enddo

if (miclevel > 2) then
   do k = ka, mza
      o_rho(k,iw) = o_press(k,iw) ** cvocp * p00kord / &
                    ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
   enddo
endif

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

! Carry out iterative hydrostatic balance procedure keeping theta constant

do iter = 1,100

! Adjust pressure at k = kbc.  Use temporal weighting for damping

   pressnew = o_press(kother,iw) * (pcol_p(kpbc) / o_press(kother,iw)) ** extrap
   o_press(kbc,iw) = .1_r8 * o_press(kbc,iw) + .9_r8 * pressnew

! Compute density for all levels

   do k = ka, mza

      if (miclevel == 0) then
         o_rho(k,iw) = o_press(k,iw) ** cvocp * p00kord / o_theta(k,iw)
         rho_tot(k) = o_rho(k,iw)
      elseif (miclevel == 1) then
         o_rho(k,iw) = o_press(k,iw) ** cvocp * p00kord / &
              ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
         rho_tot(k) = o_rho(k,iw) * (1. + o_rrw(k,iw))
      else
         exner = (real(o_press(k,iw)) * p00i) ** rocp
         tairc = exner * o_theta(k,iw) - t00
         cond  = max(0., o_rrw(k,iw) - rhovsl(tairc) / real(o_rho(k,iw)))
         rrv   = o_rrw(k,iw) - cond

         o_rho(k,iw) = o_press(k,iw) ** cvocp * p00kord / &
              ( o_theta(k,iw) * (1.0 + eps_vapi * rrv) )
         rho_tot(k) = o_rho(k,iw) * (1. + o_rrw(k,iw))
      endif

   enddo

! Integrate hydrostatic equation upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
! during iteration.  Use weighting to damp oscillations

   do k = kbc+1,mza
      pkhyd = o_press(k-1,iw) &
            - gdz_belo8(k-1) * rho_tot(k-1) - gdz_abov8(k-1) * rho_tot(k)
      o_press(k,iw) = .05_r8 * o_press(k,iw) + .95_r8 * max(.1_r8,pkhyd)
   enddo

   do k = kbc-1,ka,-1
      pkhyd = o_press(k+1,iw) &
            + gdz_belo8(k) * rho_tot(k) + gdz_abov8(k) * rho_tot(k+1)
      o_press(k,iw) = .05_r8 * o_press(k,iw) + .95_r8 * max(.1_r8,pkhyd)
   enddo

enddo

end subroutine vterpp_s
