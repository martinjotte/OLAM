subroutine swtc_init()

  use mem_basic,   only: press, rho, vc, vmc
  use mem_ijtabs,  only: jtab_w, jtab_v, itab_v, jtv_init, jtw_init
  use misc_coms,   only: iparallel
  use consts_coms, only: gravo2, gravi, erad, pio180, pi2, omega, pi1
  use mem_grid,    only: vnx, vny, glatw, glonw, xev, yev, zev
  use oname_coms,  only: nl
  use olam_mpi_atm,only: mpi_send_w, mpi_recv_w,  &
                         mpi_send_v, mpi_recv_v
  use obnd,        only: lbcopy_w, lbcopy_v

  implicit none

  integer :: j,iw,iv,iw1,iw2,im1,im2

  real :: u0_swtc, uv01dx, uv01dy, raxis
  real :: xseg, yseg, zseg, glats, glons, useg, vseg
  real :: rad0_swtc, rad_swtc, topo_swtc

! This subroutine initializes prognostic variables for shallow-water
! test cases 1, 2, and 5.  A global spherical domain (mdomain = 1) is assumed.

! Equatorial wind speeds for SWTC 1 & 2 and for 5

  if (nl%test_case == 1 .or. nl%test_case == 2) then
     u0_swtc = pi2 * erad / (12. * 86400.)
  else
     u0_swtc  = 20.
  endif

!----------------------------------------------------------------------
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

     if (nl%test_case == 1 .or. nl%test_case == 2) then

! Initial height and mass fields for SWTC 1 & 2

        press(2,iw) = 29400. - (erad * omega * u0_swtc + .5 * u0_swtc ** 2) &
                    * (sin(glatw(iw) * pio180)) ** 2

        press(2,iw) = press(2,iw) * gravi

        rho(2,iw) = press(2,iw)

     else

! Initial height and mass fields for SWTC 5

        rad0_swtc = pi1 / 9.

        rad_swtc = sqrt((glonw(iw) * pio180 + 0.5 * pi1)**2 &
                 + (glatw(iw) * pio180 - pi1 / 6.) ** 2)

        topo_swtc = max(0., 2000. * (1. - rad_swtc / rad0_swtc))

        press(2,iw) = 58408. - (erad * omega * u0_swtc + .5 * u0_swtc ** 2) &
                    * (sin(glatw(iw) * pio180)) ** 2

        press(2,iw) = press(2,iw) * gravi

        rho(2,iw) = press(2,iw) - topo_swtc

     endif

  enddo

  if (iparallel == 1) then
     call mpi_send_w(dvara1=press, dvara2=rho)
     call mpi_recv_w(dvara1=press, dvara2=rho)
  endif

! LBC copy

  call lbcopy_w(d1=press, d2=rho)

! Initialize VMC, VC

!----------------------------------------------------------------------
  do j = 1,jtab_v(jtv_init)%jend; iv = jtab_v(jtv_init)%iv(j)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

     if (iw1 == 1) iw1 = iw2
     if (iw2 == 1) iw2 = iw1

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

! V point coordinates and normal vector components

     xseg = xev(iv)
     yseg = yev(iv)
     zseg = zev(iv)

     raxis = sqrt(xseg ** 2 + yseg ** 2)

     if (raxis > 1.e3) then

        glats = atan2(zseg,raxis)
        glons = atan2(yseg,xseg)

        useg = u0_swtc * cos(glats)
        vseg = 0.

        uv01dx = -useg * yseg / raxis
        uv01dy =  useg * xseg / raxis

        vc(2,iv) = uv01dx * vnx(iv) + uv01dy * vny(iv)

     else

        vc(2,iv) = 0.

     endif

     vmc(2,iv) = vc(2,iv) * .5 * (rho(2,iw1) + rho(2,iw2))

! For below-ground points, set VC to LPV value.

     vc(1,iv) = vc(2,iv)

  enddo

! MPI parallel send/recv of V group

  if (iparallel == 1) then
     call mpi_send_v(rvara1=vmc, rvara2=vc)
     call mpi_recv_v(rvara1=vmc, rvara2=vc)
  endif

  call lbcopy_v(vmc=vmc, vc=vc)

end subroutine swtc_init
