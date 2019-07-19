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

contains

!===============================================================================

subroutine diag_earth_vels(mrl, vxesc, vyesc, vzesc)

  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use mem_ijtabs,   only: jtab_w, itab_w, jtw_prog, jtw_wadj, mrl_ends, mrl_endl
  use mem_basic,    only: vxe, vye, vze
  use misc_coms,    only: iparallel
  use mem_grid,     only: mza, mwa

  implicit none

  integer, intent(in)    :: mrl
  real,    intent(inout) :: vxesc(mza,mwa), vyesc(mza,mwa), vzesc(mza,mwa)
  integer                :: j, iw, k

  if (mrl == 0) return

  !$omp parallel do private(iw,k)
  do j = 1,jtab_w(jtw_wadj)%jend(mrl); iw = jtab_w(jtw_wadj)%iw(j)

     do k = 2, mza
        vxesc(k,iw) = vxesc(k,iw) + 0.5 * vxe(k,iw)
        vyesc(k,iw) = vyesc(k,iw) + 0.5 * vye(k,iw)
        vzesc(k,iw) = vzesc(k,iw) + 0.5 * vze(k,iw)
     enddo

  enddo
  !$omp end parallel do

  call diagvel_t3d(mrl)

  !$omp parallel do private(iw,k)
  do j = 1,jtab_w(jtw_wadj)%jend(mrl); iw = jtab_w(jtw_wadj)%iw(j)

     do k = 2, mza
        vxesc(k,iw) = vxesc(k,iw) + 0.5 * vxe(k,iw)
        vyesc(k,iw) = vyesc(k,iw) + 0.5 * vye(k,iw)
        vzesc(k,iw) = vzesc(k,iw) + 0.5 * vze(k,iw)
     enddo

  enddo
  !$omp end parallel do

end subroutine diag_earth_vels

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
  use mem_basic,  only: vc, wc, vxe, vye, vze, vxe2, vye2, vze2
  use mem_grid,   only: mza, lpw, lve2, lpv, wnx, wny, wnz

  implicit none

  integer, intent(in) :: iw
  integer             :: npoly,ka,k,jv,iv,ksw,kbv
  real                :: wst

  npoly = itab_w(iw)%npoly
  ka = lpw(iw)

! Vertical loop over T levels

  do k = ka, mza

! Diagnose 3D earth-velocity vector at T points; W contribution first

     wst = 0.5 * (wc(k-1,iw) + wc(k,iw))

     vxe(k,iw) = wst * wnx(iw)
     vye(k,iw) = wst * wny(iw)
     vze(k,iw) = wst * wnz(iw)

  enddo

! Effective contribution from submerged V faces

  if (lve2(iw) > 0) then
     do ksw = 1,lve2(iw)
        k = ka + ksw - 1

        vxe(k,iw) = vxe(k,iw) + vxe2(ksw,iw)
        vye(k,iw) = vye(k,iw) + vye2(ksw,iw)
        vze(k,iw) = vze(k,iw) + vze2(ksw,iw)
     enddo
  endif

! Loop over V neighbors of this W cell

  do jv = 1, npoly

     iv  = itab_w(iw)%iv(jv)
     kbv = lpv(iv)

! Vertical loop over V levels that are above ground

     do k = kbv, mza

! Diagnose 3D earth-velocity vector at T points; VC contribution

        vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vc(k,iv)
        vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vc(k,iv)
        vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vc(k,iv)
     enddo

  enddo

  vxe(1:ka-1,iw) = vxe(ka,iw)
  vye(1:ka-1,iw) = vye(ka,iw)
  vze(1:ka-1,iw) = vze(ka,iw)

end subroutine vel_t3d_hex

!===============================================================================

subroutine diagvel_t3d_init(mrl)

use mem_basic,  only: vc, vxe2, vye2, vze2
use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
use mem_grid,   only: lpw, lve2, lpv

implicit none

integer, intent(in) :: mrl

integer :: j,iw,npoly,ka,k,jv,iv,ksw,kbv

if (mrl == 0) return

! Horizontal loop over W columns

!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,ka,k,jv,iv,kbv,ksw)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   npoly = itab_w(iw)%npoly
   ka = lpw(iw)

   vxe2(:,iw) = 0.
   vye2(:,iw) = 0.
   vze2(:,iw) = 0.

   if (lve2(iw) > 0) then

! Loop over adjacent V faces

      do jv = 1, npoly

         iv  = itab_w(iw)%iv(jv)
         kbv = lpv(iv)

! Check if any V faces are below ground

         if (ka < kbv) then
            do k = ka, kbv-1
               ksw = k - ka + 1

! Project INITIAL VC from below-ground V faces back to (vxe2, vye2, vze2)

               vxe2(ksw,iw) = vxe2(ksw,iw) + itab_w(iw)%ecvec_vx(jv) * vc(kbv,iv)
               vye2(ksw,iw) = vye2(ksw,iw) + itab_w(iw)%ecvec_vy(jv) * vc(kbv,iv)
               vze2(ksw,iw) = vze2(ksw,iw) + itab_w(iw)%ecvec_vz(jv) * vc(kbv,iv)
            enddo
         endif

      enddo
   endif

enddo
!$omp end parallel do

end subroutine diagvel_t3d_init

end module vel_t3d
