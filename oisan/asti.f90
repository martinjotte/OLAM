subroutine vterpp()

  use mem_ijtabs,  only: jtab_w, jtw_init
  implicit none

  logical, parameter :: tconst = .false.  ! true:  keep tair  constant during adjustent
                                          ! false: keep theta constant during adjustent
                                          ! maybe use namelist parameter?
  integer :: j, iw

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Perform iterative hydrostatic balance on analysis arrays

     if (tconst) then
        call vterpp_tair(iw)
     else
        call vterpp_theta(iw)
     endif

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

              ! THIL is defined using mixing ratios in Tripoli and Cotton
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

subroutine vterpp_tair(iw)

  use isan_coms,   only: o_rho, o_press, o_theta, o_rrw, pbc, z_pbc
  use consts_coms, only: p00kord, p00i, t00, rocp, eps_vapi, rdry, gordry, pi1
  use mem_grid,    only: mza, lpw, zt, zm, gdz_belo, gdz_abov
  use micro_coms,  only: miclevel
  use therm_lib,   only: rhovsl

  implicit none

  ! TCONST: hold temperature constant during hydrostatic balancing.
  ! This preserves the temperature at each (geopotential) height.

  integer, intent(in)  :: iw

  real    :: rho_tot(mza), pold(mza), delp(mza), tair(mza), rhov_sat(mza)
  integer :: k, kbc, iter, ka
  real    :: exner, x0, x1, pkhyd, rrv

  ! Find olam layer that contains height z_pbc

  do kbc = 2, mza-1
     if (zm(kbc) >= z_pbc(iw)) exit
  enddo

  ka = min(kbc,lpw(iw))

  ! Compute air temperature as it will be held constant

  do k = ka, mza
     exner   = (o_press(k,iw) * p00i) ** rocp
     tair(k) = exner * o_theta(k,iw)
  enddo

  ! Compute water vapor density at saturation

  if (miclevel >= 2) then
     do k = ka, mza
        rhov_sat(k) = rhovsl( tair(k)-t00 )
     enddo
  endif

  ! Compute pressure closest to the z_pbc level from hypsometric equation
  ! This pressure will be held constant through the iterations.

  if (miclevel == 0) then

     o_press(kbc,iw) = pbc * exp( gordry * (z_pbc(iw) - zt(kbc)) / tair(kbc) )

     o_rho(kbc,iw) = o_press(kbc,iw) / ( tair(kbc) * rdry )

  else

     o_press(kbc,iw) = pbc * exp( gordry * (z_pbc(iw) - zt(kbc)) * (1. + o_rrw(kbc,iw))   &
                              / (tair(kbc) * (1.0 + eps_vapi * o_rrw(kbc,iw))) )

     o_rho(kbc,iw) = o_press(kbc,iw) &
                   / ( tair(kbc) * rdry * (1.0 + eps_vapi * o_rrw(kbc,iw)) )
  endif

  ! Check if there is condensate at kbc level, and adjust P and RHO if necessary

  if (miclevel >= 2) then
     if ( o_rrw(kbc,iw) * o_rho(kbc,iw) > rhov_sat(kbc) ) then

        do iter = 1, 10

           rrv = o_rrw(kbc,iw) - max(0., o_rrw(kbc,iw) - rhov_sat(kbc) / o_rho(kbc,iw))

           pkhyd = pbc * exp( gordry * (z_pbc(iw) - zt(kbc)) * (1. + o_rrw(kbc,iw)) &
                            / ( tair(kbc) * (1.0 + eps_vapi * rrv) ) )

           o_press(kbc,iw) = .25 * o_press(kbc,iw) + .75 * pkhyd

           o_rho(kbc,iw) = o_press(kbc,iw) &
                         / ( tair(kbc) * rdry * (1.0 + eps_vapi * rrv) )
        enddo

     endif
  endif

  ! Carry out iterative hydrostatic balance procedure keeping tair constant

  do iter = 1, 100

     pold(ka:mza) = o_press(ka:mza,iw)

     ! Vary adjustment weighting factor around a mean of .85
     ! (helps to contain oscillations)

     x1 = 0.85 + 0.125 * sin( pi1 * (real(iter-1) / 6. - 0.5) )
     x0 = 1.0 - x1

     if (miclevel == 0) then

        ! No influence of water on thermodynamics
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw)
        enddo

     else

        ! Air density includes any water
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw) + o_rho(k,iw) * o_rrw(k,iw)
        enddo

     endif

     ! Integrate hydrostatic equation upward and downward from kbc level
     ! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
     ! during iteration.  Use weighting to damp oscillations

     do k = kbc+1, mza
        pkhyd = o_press(k-1,iw) &
                 - ( gdz_belo(k-1) * rho_tot(k-1) + gdz_abov(k-1) * rho_tot(k) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * max(.1,pkhyd)
     enddo

     do k = kbc-1, ka, -1
        pkhyd = o_press(k+1,iw) &
                 + ( gdz_belo(k) * rho_tot(k) + gdz_abov(k) * rho_tot(k+1) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * pkhyd
     enddo

     ! Compute density for all levels using new pressure

     if (miclevel == 0) then

        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw) / ( tair(k) * rdry )
        enddo

     elseif (miclevel == 1) then

        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw) &
                       / ( tair(k) * rdry * (1.0 + eps_vapi * o_rrw(k,iw)) )
        enddo

     else  ! (miclevel == 2,3)

        do k = ka, mza
           rrv = o_rrw(k,iw) - max(0., o_rrw(k,iw) - rhov_sat(k) / o_rho(k,iw))

           o_rho(k,iw) = o_press(k,iw) &
                       / ( tair(k) * rdry * (1.0 + eps_vapi * rrv) )
        enddo

     endif

     ! Exit if pressure has converged

     if (iter >= 5) then

        do k = ka, mza
           delp(k) = abs( pold(k) - o_press(k,iw) )
        enddo

        if ( maxval( delp(ka:mza) / o_press(ka:mza,iw) ) < 2.e-6 .or. &
             maxval( delp(ka:mza) ) < 1.e-2 ) exit

     endif

  enddo

  ! Recompute theta

  do k = ka, mza
     exner = (o_press(k,iw) * p00i) ** rocp
     o_theta(k,iw) = tair(k) / exner
  enddo

end subroutine vterpp_tair

!===============================================================================

subroutine vterpp_theta(iw)

  use isan_coms,   only: o_rho, o_press, o_theta, o_rrw, pbc, z_pbc
  use consts_coms, only: p00kord, p00i, t00, cvocp, rocp, eps_vapi, cpor, gocp, pi1
  use mem_grid,    only: mza, lpw, zt, zm, gdz_belo, gdz_abov
  use micro_coms,  only: miclevel
  use therm_lib,   only: rhovsl

  implicit none

  integer, intent(in)  :: iw

  ! THETACONST: hold theta constant during hydrostatic balancing.
  ! This preserves the temperature at each pressure level.

  real    :: rho_tot(mza), pold(mza), delp(mza)
  integer :: k, kbc, iter, ka
  real    :: exner, tairc, x0, x1, pkhyd, rrv

  ! Find olam layer that contains height z_pbc

  do kbc = 2, mza-1
     if (zm(kbc) >= z_pbc(iw)) exit
  enddo

  ka = min(kbc,lpw(iw))

  ! Compute pressure closest to the z_pbc level from hypsometric equation
  ! This will be held constant through the iterations.

  exner = (pbc * p00i)**rocp

  if (miclevel == 0) then

     o_press(kbc,iw) = pbc * (1. + gocp * (z_pbc(iw) - zt(kbc)) &
                                 / (o_theta(kbc,iw) * exner) )**cpor

     o_rho(kbc,iw) = o_press(kbc,iw)**cvocp * p00kord / o_theta(kbc,iw)

  else

     o_press(kbc,iw) = pbc * (1. + gocp * (z_pbc(iw) - zt(kbc)) * (1. + o_rrw(kbc,iw)) &
                                 / (o_theta(kbc,iw) * exner * (1.0 + eps_vapi * o_rrw(kbc,iw))) )**cpor

     o_rho(kbc,iw) = o_press(kbc,iw)**cvocp * p00kord &
                   / ( o_theta(kbc,iw) * (1.0 + eps_vapi * o_rrw(kbc,iw)) )
  endif

  ! Check if there is condensate at kbc level, and adjust P and RHO if necessary

  if (miclevel >= 2) then

     tairc = o_theta(kbc,iw) * (o_press(kbc,iw) * p00i)**rocp - t00
     if ( o_rrw(kbc,iw) * o_rho(kbc,iw) > rhovsl(tairc) ) then

        do iter = 1, 10

           tairc = o_theta(kbc,iw) * (o_press(kbc,iw) * p00i)**rocp - t00

           rrv = o_rrw(kbc,iw) - max(0., o_rrw(kbc,iw) - rhovsl(tairc) / o_rho(kbc,iw))

           pkhyd = pbc * (1. + gocp * (z_pbc(iw) - zt(kbc)) * (1. + o_rrw(kbc,iw)) &
                             / (o_theta(kbc,iw) * exner * (1.0 + eps_vapi * rrv)) )**cpor

           o_press(kbc,iw) = .25 * o_press(kbc,iw) + .75 * pkhyd

           o_rho(kbc,iw) = o_press(kbc,iw)**cvocp * p00kord &
                         / ( o_theta(kbc,iw) * (1.0 + eps_vapi * rrv) )
        enddo

     endif
  endif

  ! Carry out iterative hydrostatic balance procedure keeping theta or tair constant

  do iter = 1, 100

     pold(ka:mza) = o_press(ka:mza,iw)

     ! Vary adjustment weighting factor around a mean of .85
     ! (helps to contain oscillations)

     x1 = 0.85 + 0.125 * sin( pi1 * (real(iter-1) / 6. - 0.5) )
     x0 = 1.0 - x1

     if (miclevel == 0) then

        ! No influence of water on thermodynamics
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw)
        enddo

     else

        ! Air density includes any water
        do k = ka, mza
           rho_tot(k) = o_rho(k,iw) + o_rho(k,iw) * o_rrw(k,iw)
        enddo

     endif

     ! Integrate hydrostatic equation upward and downward from kbc level
     ! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
     ! during iteration.  Use weighting to damp oscillations

     do k = kbc+1, mza
        pkhyd = o_press(k-1,iw) &
                 - ( gdz_belo(k-1) * rho_tot(k-1) + gdz_abov(k-1) * rho_tot(k) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * max(.1,pkhyd)
     enddo

     do k = kbc-1, ka, -1
        pkhyd = o_press(k+1,iw) &
                 + ( gdz_belo(k) * rho_tot(k) + gdz_abov(k) * rho_tot(k+1) )

        o_press(k,iw) = x0 * o_press(k,iw) + x1 * pkhyd
     enddo

     ! Compute density for all levels using new pressure

     if (miclevel == 0) then

        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / o_theta(k,iw)
        enddo

     elseif (miclevel == 1) then

        do k = ka, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord &
                       / ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
        enddo

     else  ! (miclevel == 2,3)

        do k = ka, mza
           tairc = o_theta(k,iw) * (o_press(k,iw) * p00i)**rocp - t00

           rrv = o_rrw(k,iw) - max(0., o_rrw(k,iw) - rhovsl(tairc) / o_rho(k,iw))

           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord &
                       / ( o_theta(k,iw) * (1.0 + eps_vapi * rrv) )
        enddo

     endif



     if (iter >= 5) then

        do k = ka, mza
           delp(k) = abs( pold(k) - o_press(k,iw) )
        enddo

        if ( maxval( delp(ka:mza) / o_press(ka:mza,iw) ) < 2.e-6 .or. &
             maxval( delp(ka:mza) ) < 1.e-2 ) exit

     endif

  enddo

end subroutine vterpp_theta
