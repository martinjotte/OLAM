module vel_t3d

contains

!===============================================================================

subroutine diagvel_t3d()

  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use misc_coms,    only: iparallel
  use mem_ijtabs,   only: jtab_w, jtw_prog
  use mem_basic,    only: vxe, vye, vze

  implicit none

  integer :: j, iw

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call vel_t3d_hex(iw)

  enddo
  !$omp end parallel do

! Parallel send-recieve of Earth Cartesian velocities

  if (iparallel == 1) then
     call mpi_send_w(rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(rvara1=vxe, rvara2=vye, rvara3=vze)
  endif
  call lbcopy_w(a1=vxe, a2=vye, a3=vze)

end subroutine diagvel_t3d

!===============================================================================

subroutine vel_t3d_hex(iw)

  use mem_ijtabs, only: itab_w
  use mem_basic,  only: vc, wc, vxe, vye, vze
  use mem_grid,   only: mza, lpw, lpv, vnx, vny, vnz, wnxo2, wnyo2, wnzo2, &
                        lve2, nve2_max
  implicit none

  integer, intent(in) :: iw
  integer             :: npoly,ka,k,jv,iv,ksw,kbv
  real                :: wst, vcwall
  real                :: vxe1(nve2_max), vye1(nve2_max), vze1(nve2_max)

  npoly = itab_w(iw)%npoly
  ka    = lpw(iw)

  ! Store prognosed velocities in cells with closed V faces

  do ksw = 1, lve2(iw)
     k = ksw + lpw(iw) - 1
     vxe1(ksw) = vxe(k,iw)
     vye1(ksw) = vye(k,iw)
     vze1(ksw) = vze(k,iw)
  enddo

  ! W contributions

  do k = ka, mza
     wst       = wc(k-1,iw) + wc(k,iw)
     vxe(k,iw) = wst * wnxo2(iw)
     vye(k,iw) = wst * wnyo2(iw)
     vze(k,iw) = wst * wnzo2(iw)
  enddo

  ! VC contribution

  do jv = 1, npoly
     iv  = itab_w(iw)%iv(jv)
     kbv = lpv(iv)

     ! Vertical loop over V levels that are above ground
     do k = kbv, mza
        vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vc(k,iv)
        vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vc(k,iv)
        vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vc(k,iv)
     enddo

     ! Vertical loop over all closed V faces
     do k = ka, kbv-1
        ksw = k - ka + 1

        ! Use projected vcwall from prognosed cell-centered velocities
        vcwall = vnx(iv) * vxe1(ksw) &
               + vny(iv) * vye1(ksw) &
               + vnz(iv) * vze1(ksw)

        vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vcwall
        vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vcwall
        vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vcwall
     enddo

  enddo

  ! Set underground points to values at first level

  vxe(1:ka-1,iw) = vxe(ka,iw)
  vye(1:ka-1,iw) = vye(ka,iw)
  vze(1:ka-1,iw) = vze(ka,iw)

end subroutine vel_t3d_hex

!===============================================================================

subroutine diag_uzonal_umerid()

  use mem_grid,  only: mwa
  use misc_coms, only: mdomain

  implicit none

  integer :: iw

  if (mdomain > 1) return

  !$omp parallel do private(iw)
  do iw = 2, mwa
     call diag_uzonal_umerid_iw(iw)
  enddo
  !$omp end parallel do

end subroutine diag_uzonal_umerid

!===============================================================================

subroutine diag_uzonal_umerid_iw(iw)

  use mem_basic, only: vxe, vye, vze, ue, ve
  use mem_grid,  only: lpw, mza, vxn_ew, vyn_ew, vxn_ns, vyn_ns, vzn_ns

  implicit none

  integer, intent(in) :: iw
  integer             :: k

  ! Reconstruct UZONAL and UMERID from VXE, VYE, VZE

  do k = lpw(iw), mza
     ue(k,iw) = vxe(k,iw) * vxn_ew(iw) + vye(k,iw) * vyn_ew(iw)
     ve(k,iw) = vxe(k,iw) * vxn_ns(iw) + vye(k,iw) * vyn_ns(iw) &
              + vze(k,iw) * vzn_ns(iw)
  enddo

  do k = 2, lpw(iw) - 1
     ue(k,iw) = ue(lpw(iw),iw)
     ve(k,iw) = ve(lpw(iw),iw)
  enddo

end subroutine diag_uzonal_umerid_iw

!===============================================================================

end module vel_t3d
