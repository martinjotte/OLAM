subroutine vterpp()

  use mem_ijtabs,  only: jtab_w, jtw_init
  implicit none

  integer :: j, iw

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Perform iterative hydrostatic balance on analysis arrays
     call vterpp_s(iw)

  enddo
  !$omp end parallel do

end subroutine vterpp



subroutine isnstage()

  use mem_basic,    only: vmc, vc, thil, rr_w, rr_v, vxe, vye, vze, tair, &
                          wmc, wc, theta, rho, press, ue, ve, vc, vmc
  use mem_grid,     only: mza, lpv, lpw, lve2, zm, zt, vcn_ew, vcn_ns, &
                          vxn_ew, vxn_ns, vyn_ew, vyn_ns, vzn_ns
  use mem_ijtabs,   only: jtab_w, jtw_init, jtab_v, itab_v, jtv_init
  use consts_coms,  only: rocp, p00i, alvlocp, t00
  use misc_coms,    only: iparallel, io6, i_o3
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use obnd,         only: lbcopy_v, lbcopy_w
  use mem_micro,    only: rr_c, con_c, cldnum
  use micro_coms,   only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
  use isan_coms,    only: o_rho, o_press, o_theta, o_rrw, o_uzonal, o_umerid, &
                          o_ozone
  use var_tables,   only: scalar_tab
  use therm_lib,    only: rhovsl

  implicit none

  integer :: j, iw, k, ka, iv, iw1, iw2
  real    :: cond, tt, ccn

  ! Copy to model arrays

  !$omp parallel do private(iw,ka,k,cond,tt,ccn)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ka = lpw(iw)

     wc (1:mza,iw) = 0.0
     wmc(1:mza,iw) = 0.0

     do k = ka, mza

        rho  (k,iw) = o_rho   (k,iw)
        rr_w (k,iw) = o_rrw   (k,iw)
        theta(k,iw) = o_theta (k,iw)
        thil (k,iw) = o_theta (k,iw)
        press(k,iw) = o_press (k,iw)
        tair (k,iw) = o_theta (k,iw) * (o_press(k,iw) * p00i) ** rocp
        rr_v (k,iw) = o_rrw   (k,iw)
        ue   (k,iw) = o_uzonal(k,iw)
        ve   (k,iw) = o_umerid(k,iw)

     enddo

     ! Earth-cartesian velocities only needed for initializing "underground"
     ! shaved-cell V faces for the new prognostic shaved-cell method.

     do k = ka, ka + lve2(iw) - 1

        vxe(k,iw) = vxn_ew(iw) * o_uzonal(k,iw) + vxn_ns(iw) * o_umerid(k,iw)
        vye(k,iw) = vyn_ew(iw) * o_uzonal(k,iw) + vyn_ns(iw) * o_umerid(k,iw)
        vze(k,iw) =                               vzn_ns(iw) * o_umerid(k,iw)

     enddo

     if (miclevel >= 2) then

        do k = ka, mza
           cond = rr_w(k,iw) * o_rho(k,iw) - rhovsl(tair(k,iw)-t00)

           if (cond > rxmin(1)) then
              rr_c(k,iw) = cond / o_rho(k,iw)
              rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)

              ! THIIL is defined using mixing ratios in Tripoli and Cotton
              tt         = max(tair(k,iw),253.)
              thil(k,iw) = theta(k,iw) * tt / (tt + alvlocp * rr_c(k,iw))
           endif

        enddo
     endif

     do k = 1, ka-1
        thil (k,iw) = thil (ka,iw)
        theta(k,iw) = theta(ka,iw)
        rho  (k,iw) = rho  (ka,iw)
        ue   (k,iw) = ue   (ka,iw)
        ve   (k,iw) = ve   (ka,iw)
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
           endif
        enddo
     endif

  enddo
  !$omp end parallel do

  ! LBC copy (THETA and TAIR and OZONE will be copied later with the scalars)

  if (iparallel == 1) then
     call mpi_send_w(dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil, rvara4=ue, rvara5=ve)

     call mpi_recv_w(dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil, rvara4=ue, rvara5=ve)
  endif

  call lbcopy_w(a1=wc, a2=wmc, a3=thil, a4=ue, a5=ve, d1=press, d2=rho)

  ! Initialize VMC, VC

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_init)%jend; iv = jtab_v(jtv_init)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        vc (k,iv) = 0.5 * ( (ue(k,iw1) + ue(k,iw2)) * vcn_ew(iv) &
                          + (ve(k,iw2) + ve(k,iw2)) * vcn_ns(iv) )
        vmc(k,iv) = vc(k,iv) * 0.5 * real(rho(k,iw1) + rho(k,iw2))
     enddo

     ! For below-ground points, set VC to 0

     vc (1:lpv(iv)-1,iv) = 0.
     vmc(1:lpv(iv)-1,iv) = 0.

  enddo
  !$omp end parallel do

  ! MPI parallel send/recv of V group

  if (iparallel == 1) then
     call mpi_send_v(rvara1=vmc, rvara2=vc)
     call mpi_recv_v(rvara1=vmc, rvara2=vc)
  endif

  ! LBC copy of VMC, VC

  call lbcopy_v(vmc=vmc, vc=vc)

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

  use isan_coms,   only: o_rho, o_press, o_theta, o_rrw, pbc, z_pbc
  use consts_coms, only: cvocp, p00kord, p00i, t00, rocp, eps_vapi
  use mem_grid,    only: mza, lpw, zt, gdz_belo, gdz_abov
  use micro_coms,  only: miclevel
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in) :: iw

  real    :: rho_tot(mza), pkhyd(mza), delp(mza)
  integer :: k, kbc, iter, ka, khi, klo, kother
  real    :: exner, tairc, cond, rrv, x0, x1, extrap

  ! Find olam levels that bracket height z_pbc

  do khi = 2, mza-1
     if (zt(khi) >= z_pbc(iw)) exit
  enddo
  klo = khi - 1

  ! Determine whether zt(klo) or zt(khi) is closer to pcol_z(kpbc).  The closer
  ! one will undergo direct pressure adjustment to satisfy the internal b.c.

  if (zt(khi) - z_pbc(iw) < z_pbc(iw) - zt(klo)) then
     kbc = khi
     kother = klo
  else
     kbc = klo
     kother = khi
  endif

  ka = min(kother,lpw(iw))

  ! Logarithmic interpolation factor for pressure boundary condition

  extrap = (zt(kbc) - zt(kother)) / (z_pbc(iw) - zt(kother))

  ! Carry out iterative hydrostatic balance procedure keeping theta constant

  do iter = 1, 100

     ! Slowly ramp up adjustment weighting factor

     x0 = 0.02 * real(min(iter,20))
     x1 = 1.0 - x0

     if (miclevel == 0) then

        ! No influence of water on thermodynamics
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw)
        enddo

     else

        ! Air density includes water
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw) + o_rho(k,iw) * o_rrw(k,iw)
        enddo

     endif

     ! Adjust pressure at k = kbc.  Use temporal weighting for damping

     pkhyd(kbc) = o_press(kother,iw) * ( pbc / o_press(kother,iw) )**extrap

     o_press(kbc,iw) = x0 * o_press(kbc,iw) + x1 * pkhyd(kbc)

     ! Integrate hydrostatic equation upward and downward from kbc level
     ! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
     ! during iteration.  Use weighting to damp oscillations

     do k = kbc+1, mza
        pkhyd(k) = o_press(k-1,iw) &
                 - ( gdz_belo(k-1) * rho_tot(k-1) + gdz_abov(k-1) * rho_tot(k) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * max(.1,pkhyd(k))
     enddo

     do k = kbc-1, ka, -1
        pkhyd(k) = o_press(k+1,iw) &
                 + ( gdz_belo(k) * rho_tot(k) + gdz_abov(k) * rho_tot(k+1) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * pkhyd(k)
     enddo

     ! Compute density for all levels using new pressure

     if (miclevel == 0) then

        ! No influence of water on thermodynamics
        do k = ka, mza
           o_rho(k,iw) = rho_tot(k)
        enddo

     elseif (miclevel == 1) then

        ! Only assume water vapor effects on thermodynamics (no condensate)
        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / &
                         ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
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
        enddo

     endif

     ! Exit if pressure has converged after at least 8 iterations

     if (iter >= 8) then

        do k = ka, mza
           delp(k) = abs( pkhyd(k) - o_press(k,iw) )
        enddo

        if ( maxval( delp(ka:mza) / o_press(ka:mza,iw) ) < 1.e-6 .or. &
             maxval( delp(ka:mza) ) < 1.e-2 ) exit

     endif

  enddo

end subroutine vterpp_s
