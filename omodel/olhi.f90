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
subroutine fldslhi()

use mem_basic,   only: theta, thil, rho, press, rr_w, &
                       wc, wmc, vc, vp, vmc, vmp, ue, ve
use mem_ijtabs,  only: jtab_w, jtab_v, itab_v, &
                       jtv_init, jtw_init
!use consts_coms, only: p00, p00i, rocp, cvocp, p00kord, rdry, alvlocp, &
!                       eps_vapi, r8
use mem_grid,    only: mza, lpv, vcn_ew, zm, zt
use mem_zonavg,  only: zonavg_init
use misc_coms,   only: io6, iparallel, idate1, imonth1, iyear1
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w,  &
                       mpi_send_v, mpi_recv_v
use obnd,        only: lbcopy_v, lbcopy_w

implicit none

integer :: j,iw,k,ka,iv,iw1,iw2
real    :: ug

! Fill zonavg arrays for initialization time

call zonavg_init(idate1,imonth1,iyear1)

! Fill and interpolate thermodynamic fields

!----------------------------------------------------------------------
!$omp parallel do private(iw)
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!----------------------------------------------------------------------

   call interp_thermo_olhi(iw)

enddo
!$omp end parallel do

! Parallel send and LBC copy (THETA and TAIR will be copied later with the scalars)

if (iparallel == 1) then
   call mpi_send_w(1, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil, &
                   rvara4=ue, rvara5=ve)

   call mpi_recv_w(1, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil, &
                   rvara4=ue, rvara5=ve)
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, a4=ue, a5=ve, d1=press, d2=rho)

! Initialize VMC, VC

!----------------------------------------------------------------------
!$omp parallel do private(iv,iw1,iw2,ka,k,ug)
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   ka = lpv(iv)

   ! Average winds to V point and rotate at V point (assumed to be global simulation)

   do k = ka, mza
      ug = .5 * (ue(k,iw1) + ue(k,iw2))
      vc( k,iv) = ug * vcn_ew(iv)
      vmc(k,iv) = vc(k,iv) * .5 * real(rho(k,iw1) + rho(k,iw2))
   enddo

   ! For below-ground points, set VC to 0

   vc (1:ka-1,iv) = 0.0
   vmc(1:ka-1,iv) = 0.0

enddo
!$omp end parallel do

! MPI parallel send/recv of V group

if (iparallel == 1) then
   call mpi_send_v(1, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(1, rvara1=vmc, rvara2=vc)
endif

! LBC copy of VMC, VC

call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

if (allocated(vmp)) vmp(:,:) = vmc(:,:)
if (allocated(vp )) vp (:,:) = vc (:,:)

! print out initial state from 1st jtw_init column

iw = jtab_w(jtw_init)%iw(1)

write(io6,*) ' '
write(io6,*) '========================================================================='
write(io6,*) '                    OLAM INITIAL STATE COLUMN (lhi)'
write(io6,*) '========================================================================='
write(io6,*) '   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   rr_w(g/kg)'
write(io6,*) '========================================================================='
write(io6,*) ' '

do k = mza,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),rr_w(k,iw)*1.e3
enddo

write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

end subroutine fldslhi



subroutine interp_thermo_olhi(iw)

  use mem_basic,   only: theta, thil, tair, rho, press, rr_w, rr_v,  &
                         wc, wmc, ue, ve
  use mem_micro,   only: rr_c, con_c, cldnum
  use micro_coms,  only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
  use consts_coms, only: p00, p00i, rocp, cvocp, p00kord, rdry, alvlocp, &
                         eps_vapi, r8
  use mem_grid,    only: mza, lpw, zt, gdz_belo8, gdz_abov8, glatw
  use mem_zonavg,  only: zonp_vect, zonz, zonu, zont, zonr
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in) :: iw

  integer :: k,ka,iter,ilat
  integer :: kpbc, kbc, kother, khi, klo

  real :: temp,exner,ccn
  real :: rlat,wt2
  real :: extrap

  real :: vctr1(mza), vctr2(mza), vctr3(mza), vctr4(mza)

  real(r8) :: pressnew, pkhyd, rho_tot(mza), dp(mza)

  real :: zonz_vect(22)
  real :: zont_vect(22)
  real :: zonu_vect(22)
  real :: zonr_vect(22)

! Choose as an internal pressure boundary condition the pressure level
! zonp_vect(3), whose pressure is 46415.89 Pa.

  kpbc = 3

! Linearly interpolate zonavg arrays by latitude to current IW column
! and K level

  ka = lpw(iw)

  rlat = .4 * (glatw(iw) + 93.75)
  ilat = int(rlat)
  wt2 = rlat - float(ilat)

  zonz_vect(1:22) = (1. - wt2) * zonz(ilat,1:22) + wt2 * zonz(ilat+1,1:22)
  zont_vect(1:22) = (1. - wt2) * zont(ilat,1:22) + wt2 * zont(ilat+1,1:22)
  zonu_vect(1:22) = (1. - wt2) * zonu(ilat,1:22) + wt2 * zonu(ilat+1,1:22)
  zonr_vect(1:22) = (1. - wt2) * zonr(ilat,1:22) + wt2 * zonr(ilat+1,1:22)

! Interpolate zonavg vector arrays in height to model levels

! Fill pressure, theta, air density, and vapor density arrays from zonavg
! vector arrays prior to iteration

  call hintrp_ee(22,zonp_vect,zonz_vect,mza,vctr1,zt) ! pressure
  call hintrp_ee(22,zont_vect,zonz_vect,mza,vctr2,zt) ! temp
  call hintrp_ee(22,zonr_vect,zonz_vect,mza,vctr3,zt) ! specific humidity
  call hintrp_ee(22,zonu_vect,zonz_vect,mza,vctr4,zt) ! uzonal

  press(1:mza,iw) = vctr1(1:mza)
  theta(1:mza,iw) = vctr2(1:mza) * (p00 / vctr1(1:mza)) ** rocp  ! temp to theta
  thil (1:mza,iw) = theta(1:mza,iw)
  rr_v (1:mza,iw) = vctr3(1:mza)  ! mixing ratio
  rr_w (1:mza,iw) = rr_v(1:mza,iw)
  ue   (1:mza,iw) = vctr4(1:mza)
  ve   (1:mza,iw) = 0.
  wc   (1:mza,iw) = 0.
  wmc  (1:mza,iw) = 0.

  if (miclevel == 0) then
     rho(1:mza,iw) = press(1:mza,iw) ** cvocp * p00kord / theta(1:mza,iw)
  else
     rho(1:mza,iw) = press(1:mza,iw) ** cvocp * p00kord &
                   / (theta(1:mza,iw) * (1.0 + eps_vapi * rr_w(1:mza,iw)))
  endif

! Make sure kpbc is above the surface

  k = ka + 1
  if (zonp_vect(kpbc) > vctr1(k)) then
     do while (zonp_vect(kpbc) > vctr1(k))
        kpbc = kpbc + 1
     enddo
  endif

! Determine which two model zt levels bracket zonz_vect(kpbc) in this column

  khi = ka + 1
  do while (zt(khi) < zonz_vect(kpbc))
     khi = khi + 1
  enddo
  klo = khi - 1

! Determine whether zt(klo) or zt(khi) is closer to zonz_vect(kpbc).  The closer
! one will undergo direct pressure adjustment to satisfy the internal b.c.

  if (zt(khi) - zonz_vect(kpbc) < zonz_vect(kpbc) - zt(klo)) then
     kbc = khi
     kother = klo
  else
     kbc = klo
     kother = khi
  endif

  extrap = (zt(kbc) - zt(kother)) / (zonz_vect(kpbc) - zt(kother))

! Iterative hydrostatic integration

  do iter = 1, 100

! Adjust pressure at k = kbc.  Use temporal weighting for damping

     pressnew = press(kother,iw) * (zonp_vect(kpbc) / press(kother,iw)) ** extrap
     dp(kbc) = pressnew - press(kbc,iw)
     press(kbc,iw) = .1_r8 * press(kbc,iw) + .9_r8 * pressnew

!  Compute density for all levels

     do k = ka,mza

! Try this: hold Mclatchy temp (vctr2) constant during iterations

        theta(k,iw) = vctr2(k) * (p00 / real(press(k,iw))) ** rocp

        if (miclevel == 0) then
           rho (k,iw) = press(k,iw) ** cvocp * p00kord / theta(k,iw)
           rho_tot(k) = rho(k,iw)
        elseif (miclevel == 1) then
           rho (k,iw) = press(k,iw) ** cvocp * p00kord / &
                ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )
           rho_tot(k) = rho(k,iw) * (1. + rr_v(k,iw))
        else
           exner = (real(press(k,iw)) * p00i) ** rocp  ! Defined WITHOUT CP factor
           temp = exner * theta(k,iw)

           rr_c(k,iw) = max(0., rr_w(k,iw) - rhovsl(temp-273.15) / real(rho(k,iw)))
           rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

           rho (k,iw) = press(k,iw) ** cvocp * p00kord / &
                ( theta(k,iw) * (1.0 + eps_vapi * rr_v(k,iw)) )

           rho_tot(k) = rho(k,iw) * (1. + rr_v(k,iw))
        endif

     enddo

! Integrate hydrostatic equation upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
! during iteration.  Use weighting to damp oscillations

     do k = kbc+1,mza
        pkhyd = press(k-1,iw) &
              - gdz_belo8(k-1) * rho_tot(k-1) - gdz_abov8(k-1) * rho_tot(k)
        dp(k) = pkhyd - press(k,iw)
        press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)
     enddo

     do k = kbc-1,ka,-1
        pkhyd = press(k+1,iw) &
              + gdz_belo8(k) * rho_tot(k) + gdz_abov8(k) * rho_tot(k+1)
        dp(k) = pkhyd - press(k,iw)
        press(k,iw) = .05_r8 * press(k,iw) + .95_r8 * max(.1_r8, pkhyd)
     enddo

! Exit if pressure has converged after at least 15 iterations

     if (iter > 15) then
        if (maxval( abs(dp(ka:mza)) / press(ka:mza,iw) ) < 1.d-12) exit
     endif

  enddo

  do k = ka, mza
     tair(k,iw) = theta(k,iw) * (real(press(k,iw)) * p00i) ** rocp
     if (miclevel <= 1) then
        thil(k,iw) = theta(k,iw)
     else
        thil(k,iw) = theta(k,iw) / (1. + alvlocp * rr_c(k,iw) / &
                                       ((1.0 + rr_c(k,iw)) * max(temp,253.)))
     endif
  enddo

  do k = 1, ka-1
     thil(k,iw) = thil(ka,iw)
  enddo

  ! If there is condensate, initialize con_c if prognosed

  if (miclevel == 3 .and. jnmb(1) == 5) then
     if (ccnparm > 1.e6) then
        ccn = ccnparm
     else
        ccn = cldnum(iw)
     endif

     do k = ka, mza
        if (rr_c(k,iw) > rxmin(1)) then
           con_c(k,iw) = ccn * real(rho(k,iw)) * zfactor_ccn(k)
        else
           con_c(k,iw) = 0.0
        endif
     enddo
  endif

end subroutine interp_thermo_olhi
