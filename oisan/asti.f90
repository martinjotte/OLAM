
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
subroutine isnstage(iaction)

  use mem_basic,    only: vmc, vc, thil, rr_w, rr_v, vxe, vye, vze, tair, &
                          wmc, wc, theta, rho, press, ue, ve, vc, vmc
  use mem_grid,     only: mza, lpv, lpw, zm, zt, vnxo2, vnyo2, vnzo2, &
                          vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns
  use mem_ijtabs,   only: jtab_w, jtw_init, jtab_v, itab_v, jtv_init
  use consts_coms,  only: rocp, p00i, alvlocp, t00
  use misc_coms,    only: iparallel, runtype, io6, i_o3
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use obnd,         only: lbcopy_v, lbcopy_w
  use mem_micro,    only: rr_c, con_c, cldnum
  use micro_coms,   only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
  use isan_coms,    only: o_rho, o_press, o_theta, o_rrw, o_uzonal, o_umerid, &
                          o_ozone
  use var_tables,   only: scalar_tab
  use therm_lib,    only: rhovsl

  implicit none

  integer, intent(in) :: iaction

  integer :: j, iw, k, ka, iv, iw1, iw2
  real    :: cond, tt, ccn

  write(io6,*)
  write(io6,*) "Performing iterative hydrostatic balancing."

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     ! Perform iterative hydrostatic balance and copy to main model arrays
     call vterpp_s(iw)

  enddo
  !$omp end parallel do

  ! Skip setting model arrays if we are only nudging and not initialising

  if (iaction == 1 .or. runtype /= 'INITIAL') return

  ! Copy to model arrays

  !$omp parallel do private(iw,ka,k,cond,tt,ccn)
  do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)

     ka = lpw(iw)

     wc (1:mza,iw) = 0.0
     wmc(1:mza,iw) = 0.0

     do k = ka, mza

        rho  (k,iw) = o_rho   (k,iw)
        rr_w (k,iw) = o_rrw   (k,iw)
        theta(k,iw) = o_theta (k,iw)
        press(k,iw) = o_press (k,iw)
        tair (k,iw) = o_theta (k,iw) * (o_press(k,iw) * p00i) ** rocp
        rr_v (k,iw) = o_rrw   (k,iw)
        rr_c (k,iw) = 0.0
        thil (k,iw) = theta   (k,iw)
        ue   (k,iw) = o_uzonal(k,iw)
        ve   (k,iw) = o_umerid(k,iw)

        vxe  (k,iw) = vxn_ew(iw) * o_uzonal(k,iw) + vxn_ns(iw) * o_umerid(k,iw)
        vye  (k,iw) = vyn_ew(iw) * o_uzonal(k,iw) + vyn_ns(iw) * o_umerid(k,iw)
        vze  (k,iw) =                               vzn_ns(iw) * o_umerid(k,iw)

     enddo

     if (miclevel >= 2) then

        do k = ka, mza
           cond = rr_w(k,iw) * o_rho(k,iw) - rhovsl(tair(k,iw)-t00)

           if (cond > rxmin(1)) then
              rr_c(k,iw) = cond / o_rho(k,iw)
              rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

              tt         = max(tair(k,iw),253.) * (1.0 + rr_v(k,iw))
              thil(k,iw) = theta(k,iw) * tt / (tt + alvlocp * rr_c(k,iw))
           endif

        enddo
     endif

     do k = 1, ka-1
        thil (k,iw) = thil (ka,iw)
        theta(k,iw) = theta(ka,iw)
        rho  (k,iw) = rho  (ka,iw)
     enddo

     ! Initialize any variable in the var_table that is named ozone.
     ! Assumed units are ppmv

     if (i_o3 > 0) then
        scalar_tab(i_o3)%var_p(ka:mza,iw) = o_ozone(ka:mza,iw)
        scalar_tab(i_o3)%var_p(1:ka-1,iw) = o_ozone(ka,iw)
     endif

     ! If there is initial cloud water and we are prognosing cloud number,
     ! initialize con_c based on default CCN. This can be overwriten later
     ! if we read in geos-chem or CMAQ CCN.

     if (miclevel == 3 .and. jnmb(1) == 5) then
        if (ccnparm > 1.e6) then
           ccn = ccnparm
        else
           ccn = cldnum(iw)
        endif

        do k = ka, mza
           if (rr_c(k,iw) * o_rho(k,iw) >= rxmin(1)) then
              con_c(k,iw) = ccn * o_rho(k,iw) * zfactor_ccn(k)
           else
              con_c(k,iw) = 0.0
           endif
        enddo
     endif

  enddo
  !$omp end parallel do

  ! LBC copy (THETA and TAIR and OZONE will be copied later with the scalars)

  if (iparallel == 1) then
     call mpi_send_w(1, dvara1=press, dvara2=rho, &
                     rvara2=wc, rvara3=wmc, rvara4=thil, &
                     rvara5=ue, rvara6=ve, &
                     rvara7=vxe, rvara8=vye, rvara9=vze)

     call mpi_recv_w(1, dvara1=press, dvara2=rho, &
                     rvara2=wc, rvara3=wmc, rvara4=thil, &
                     rvara5=ue, rvara6=ve, &
                     rvara7=vxe, rvara8=vye, rvara9=vze)
  endif

  call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, a4=ue, a5=ve, a6=vxe, a7=vye, &
                a8=vze, d1=press, d2=rho)

  ! Initialize VMC, VC

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vc(k,iv) = vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                 + vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                 + vnzo2(iv) * (vze(k,iw1) + vze(k,iw2))

        vmc(k,iv) = vc(k,iv) * 0.5 * real(rho(k,iw1) + rho(k,iw2))
     enddo

     ! For below-ground points, set VC to 0

     vc (1:lpv(iv)-1,iv) = 0.
     vmc(1:lpv(iv)-1,iv) = 0.

  enddo
  !$omp end parallel do

  ! MPI parallel send/recv of V group

  if (iparallel == 1) then
     call mpi_send_v(1, rvara1=vmc, rvara2=vc)
     call mpi_recv_v(1, rvara1=vmc, rvara2=vc)
  endif

  ! LBC copy of VMC, VC

  call lbcopy_v(1, vmc=vmc, vc=vc)

  ! Print out initial state from 1st jtw_init column

  iw = jtab_w(jtw_init)%iw(1)

  write(io6,*)' '
  write(io6,*)'========================================================================'
  write(io6,*)'                    OLAM INITIAL STATE COLUMN (vari)'
  write(io6,*)'========================================================================'
  write(io6,*)'   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   rr_w(g/kg)'
  write(io6,*)'========================================================================'
  write(io6,*)' '

  do k = mza,lpw(iw),-1
     write(io6, '(f10.2,1x,9(''-------''))') zm(k)
     write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)') &
          k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),rr_w(k,iw)*1.e3
  enddo

  write(io6, '(f10.2,1x,9(''-------''))') zm(1)
  write(io6,*)' '

end subroutine isnstage

!===============================================================================

subroutine vterpp_s(iw)

  use isan_coms,   only: o_rho, o_press, o_theta, o_rrw, npd, pcol_p, pcol_z
  use consts_coms, only: grav2, grav, cvocp, p00kord, rvap, p00, p00i, &
                         t00, cp, rocp, gravi, eps_virt, eps_vapi
  use mem_grid,    only: mza, lpw, zt, gravm, volt, volwi, dzm, &
                         gdz_belo, gdz_abov
  use misc_coms,   only: nrk_wrtv
  use micro_coms,  only: miclevel
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in) :: iw

  real    :: rho_tot(mza), v4(mza), fact(mza), dp(mza)
  integer :: k, kpbc, klo, khi, kbc, kother, iter, ka
  real    :: extrap, exner, tairc, cond, rrv
  real    :: x0, x1, pkhyd

  ka = lpw(iw)

  ! Choose as an internal pressure boundary condition the pcol_p pressure level
  ! at or below (in elevation) the 49900 Pa surface.  Find the k index of this level.

  kpbc = npd
  do while (pcol_p(kpbc) < 49900.)
     kpbc = kpbc - 1
  enddo

  ! Make sure we are above the surface though!

  if (pcol_z(kpbc,iw) < zt(ka)) then
     do while (pcol_z(kpbc,iw) < zt(ka) .and. kpbc < npd-1)
        kpbc = kpbc + 1
     enddo
  endif

  ! Determine which two model zt levels bracket pcol_z(kpbc) in this column

  khi = ka + 1
  do while (zt(khi) < pcol_z(kpbc,iw) .and. khi < mza-1)
     khi = khi + 1
  enddo
  klo = khi - 1

  ! Determine whether zt(klo) or zt(khi) is closer to pcol_z(kpbc).  The closer
  ! one will undergo direct pressure adjustment to satisfy the internal b.c.

  if (zt(khi) - pcol_z(kpbc,iw) < pcol_z(kpbc,iw) - zt(klo)) then
     kbc = khi
     kother = klo
  else
     kbc = klo
     kother = khi
  endif

  extrap = (zt(kbc) - zt(kother)) / (pcol_z(kpbc,iw) - zt(kother))

  if (nrk_wrtv > 1) then
     do k = ka, mza
        v4  (k) = real(volt(k,iw))
        fact(k) = gravm(k) * volwi(k,iw) * dzm(k)
     enddo
  endif

  if (miclevel > 1) then
     do k = ka, mza
        o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / &
                      ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
     enddo
  endif

  ! Carry out iterative hydrostatic balance procedure keeping Theta constant

  do iter = 1, 100

     x0 = 0.02 * real(min(iter,20))
     x1 = 1.0 - x0

     ! Adjust pressure at k = kbc.  Use temporal weighting for damping

     pkhyd = o_press(kother,iw) * ( pcol_p(kpbc) / o_press(kother,iw) )**extrap

     dp     (kbc)    = abs( pkhyd - o_press(kbc,iw) )
     o_press(kbc,iw) = x0 * o_press(kbc,iw) + x1 * pkhyd

     ! Compute density for all levels

     if (miclevel == 0) then

        ! No influence of water on thermodynamics
        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / o_theta(k,iw)
           rho_tot(k)  = o_rho(k,iw)
        enddo

     elseif (miclevel == 1) then

        ! Only assume water vapor effects on thermodynamics (no condensate)
        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / &
                         ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
           rho_tot(k)  = o_rho(k,iw) * (1. + o_rrw(k,iw))
        enddo

     else

        ! Water vapor and condensate allowed
        do k = ka, mza
           exner = (o_press(k,iw) * p00i) ** rocp
           tairc = exner * o_theta(k,iw) - t00
           cond  = max(0., o_rrw(k,iw) - rhovsl(tairc) / o_rho(k,iw))
           rrv   = o_rrw(k,iw) - cond

           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / &
                         ( o_theta(k,iw) * (1.0 + eps_vapi * rrv) )
           rho_tot(k)  = o_rho(k,iw) * (1. + o_rrw(k,iw))
        enddo

     endif

     ! Integrate hydrostatic equation upward and downward from kbc level
     ! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
     ! during iteration.  Use weighting to damp oscillations

     if (nrk_wrtv == 1) then

        ! Bob's single-step time integration scheme height-averages density
        ! in the hydrostatic balance

        do k = kbc+1, mza
           pkhyd = o_press(k-1,iw) &
                 - ( gdz_belo(k-1) * rho_tot(k-1) + gdz_abov(k-1) * rho_tot(k) )

           dp     (k)    = abs( pkhyd - o_press(k,iw) )
           o_press(k,iw) = x0 * o_press(k,iw) + x1 * max(.1, pkhyd)
        enddo

        do k = kbc-1, ka, -1
           pkhyd = o_press(k+1,iw) &
                 + ( gdz_belo(k) * rho_tot(k) + gdz_abov(k) * rho_tot(k+1) )

           dp     (k)    = abs( pkhyd - o_press(k,iw) )
           o_press(k,iw) = x0 * o_press(k,iw) + x1 * pkhyd
        enddo

     else

        ! Runge-Kutta time integration scheme volume-averages density
        ! in the hydrostatic balance

        do k = kbc+1, mza
           pkhyd = o_press(k-1,iw) &
                 - ( v4(k-1) * rho_tot(k-1) + v4(k) * rho_tot(k) ) * fact(k-1)

           dp     (k)    = abs( pkhyd - o_press(k,iw) )
           o_press(k,iw) = x0 * o_press(k,iw) + x1 * max(.1, pkhyd)
        enddo

        do k = kbc-1, ka, -1
           pkhyd = o_press(k+1,iw) &
                 + ( v4(k) * rho_tot(k) + v4(k+1) * rho_tot(k+1) ) * fact(k)

           dp     (k)    = abs( pkhyd - o_press(k,iw) )
           o_press(k,iw) = x0 * o_press(k,iw) + x1 * pkhyd
        enddo

     endif

     ! Exit if pressure has converged after at least 8 iterations

     if (iter >= 8) then
        if ( maxval( dp(ka:mza) / o_press(ka:mza,iw) ) < 1.e-6 .or. &
             maxval( dp(ka:mza) ) < 1.e-2 ) exit
     endif

  enddo

  ! Return if we are performing a nudging prep and not initialising model arays

end subroutine vterpp_s
