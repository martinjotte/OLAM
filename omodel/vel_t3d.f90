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

module vel_t3d

  real, allocatable :: exm(:,:)
  real, allocatable :: eym(:,:)
  real, allocatable :: ezm(:,:)

  real, allocatable :: vxe1(:,:)
  real, allocatable :: vye1(:,:)
  real, allocatable :: vze1(:,:)

  integer, parameter :: icut_vel = 1

contains

!===============================================================================

subroutine init_velt3d()

  use mem_grid,   only: mwa, nve2_max, lve2, lpw, lpv, vnx, vny, vnz, wnx, wny, wnz
  use mem_ijtabs, only: jtab_w, jtw_prog, itab_w
  use misc_coms,  only: rinit

  implicit none

  integer :: j, iw, iv, jv, k, ks

  real :: exsum(nve2_max)
  real :: eysum(nve2_max)
  real :: ezsum(nve2_max)

  real :: extot, eytot, eztot

  if (icut_vel == 1) then

     allocate (vxe1(nve2_max,mwa)) ; vxe1 = rinit
     allocate (vye1(nve2_max,mwa)) ; vye1 = rinit
     allocate (vze1(nve2_max,mwa)) ; vze1 = rinit

  else

     allocate(exm(nve2_max,mwa)) ; exm = 1.0
     allocate(eym(nve2_max,mwa)) ; eym = 1.0
     allocate(ezm(nve2_max,mwa)) ; ezm = 1.0

     !$omp parallel private(exsum,eysum,ezsum)
     !$omp do private(iw,jv,iv,k,ks,extot,eytot,eztot) schedule(guided)
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

        if (lve2(iw) > 0) then

           exsum(1:lve2(iw)) = 0.0
           eysum(1:lve2(iw)) = 0.0
           ezsum(1:lve2(iw)) = 0.0

           do jv = 1, itab_w(iw)%npoly
              iv = itab_w(iw)%iv(jv)

              do k = lpv(iv), lpw(iw) + lve2(iw) - 1
                 ks = k - lpw(iw) + 1
                 exsum(ks) = exsum(ks) + itab_w(iw)%ecvec_vx(jv) * vnx(iv)
                 eysum(ks) = eysum(ks) + itab_w(iw)%ecvec_vy(jv) * vny(iv)
                 ezsum(ks) = ezsum(ks) + itab_w(iw)%ecvec_vz(jv) * vnz(iv)
              enddo
           enddo

           extot = 1. - wnx(iw)**2
           eytot = 1. - wny(iw)**2
           eztot = 1. - wnz(iw)**2

           do ks = 1, lve2(iw)
              exm(ks,iw) = min(4., extot / max(exsum(ks), 1.e-7))
              eym(ks,iw) = min(4., eytot / max(eysum(ks), 1.e-7))
              ezm(ks,iw) = min(4., eztot / max(ezsum(ks), 1.e-7))
           enddo

        endif

     enddo
     !$omp end do
     !$omp end parallel

  endif

end subroutine init_velt3d

!===============================================================================

subroutine diagvel_t3d(mrl)

  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use misc_coms,    only: iparallel
  use mem_ijtabs,   only: jtab_w, jtw_prog
  use mem_basic,    only: vxe, vye, vze

  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iw

  if (mrl == 0) return

  !$omp parallel do private(iw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     call vel_t3d_hex(iw)

  enddo
  !$omp end parallel do

! Parallel send-recieve of Earth Cartesian velocities

  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
     call mpi_recv_w(mrl, rvara1=vxe, rvara2=vye, rvara3=vze)
  endif
  call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)

end subroutine diagvel_t3d

!===============================================================================

subroutine vel_t3d_hex(iw)

  use mem_ijtabs, only: itab_w
  use mem_basic,  only: vc, wc, vxe, vye, vze
  use mem_grid,   only: mza, lpw, lve2, lpv, vnx, vny, vnz, wnxo2, wnyo2, wnzo2

  implicit none

  integer, intent(in) :: iw
  integer             :: npoly,ka,k,jv,iv,ksw,kbv
  real                :: wst, vcwall

  npoly = itab_w(iw)%npoly
  ka    = lpw(iw)

  vxe(:,iw) = 0.0
  vye(:,iw) = 0.0
  vze(:,iw) = 0.0

  ! Diagnose 3D earth-velocity vector at T points; VC contribution

  do jv = 1, npoly
     iv  = itab_w(iw)%iv(jv)
     kbv = lpv(iv)

     ! Vertical loop over V levels that are above ground
     do k = kbv, mza
        vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vc(k,iv)
        vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vc(k,iv)
        vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vc(k,iv)
     enddo

     if (icut_vel == 1) then

        ! Vertical loop over all closed V faces
        do k = ka, kbv-1

           ! Use projected vcwall for closed faces

           ksw = k - ka + 1

           vcwall = vnx(iv) * vxe1(ksw,iw) &
                  + vny(iv) * vye1(ksw,iw) &
                  + vnz(iv) * vze1(ksw,iw)

           vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vcwall
           vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vcwall
           vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vcwall

        enddo

     endif

  enddo

  if (icut_vel /= 1) then

     ! Underground velocity contribution: extrapolate from prognosed vc
     do ksw = 1, lve2(iw)
        k = ka + ksw - 1
        vxe(k,iw) = vxe(k,iw) * exm(ksw,iw)
        vye(k,iw) = vye(k,iw) * eym(ksw,iw)
        vze(k,iw) = vze(k,iw) * ezm(ksw,iw)
     enddo

  endif

  ! W contributions

  do k = ka, mza
     wst = wc(k-1,iw) + wc(k,iw)
     vxe(k,iw) = vxe(k,iw) + wst * wnxo2(iw)
     vye(k,iw) = vye(k,iw) + wst * wnyo2(iw)
     vze(k,iw) = vze(k,iw) + wst * wnzo2(iw)
  enddo

  vxe(1:ka-1,iw) = vxe(ka,iw)
  vye(1:ka-1,iw) = vye(ka,iw)
  vze(1:ka-1,iw) = vze(ka,iw)

end subroutine vel_t3d_hex

!===============================================================================

subroutine diagvel_t3d_init(mrl)

  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
  use mem_grid,   only: lpw, lve2, lpv
  use mem_basic,  only: vxe, vye, vze

  implicit none

  integer, intent(in) :: mrl

  integer :: j,iw,ka,k,ksw

  if (mrl == 0 .or. icut_vel /= 1) return

  ! Horizontal loop over W columns

  !----------------------------------------------------------------------
  !$omp parallel do private(iw,ka,k,ksw)
  do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
  !----------------------------------------------------------------------

     if (lve2(iw) > 0) then

        ka = lpw(iw)
        do ksw = 1, lve2(iw)
           k = ksw + ka - 1

           vxe1(ksw,iw) = vxe(k,iw)
           vye1(ksw,iw) = vye(k,iw)
           vze1(ksw,iw) = vze(k,iw)
        enddo

     endif

  enddo
  !$omp end parallel do

end subroutine diagvel_t3d_init

!===============================================================================

subroutine diag_uzonal_umerid(mrl)

  use mem_grid, only: mwa
  implicit none

  integer, intent(in) :: mrl
  integer             :: iw

  if (mrl == 0) return

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


end module vel_t3d
