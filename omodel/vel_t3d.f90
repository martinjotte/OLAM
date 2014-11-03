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
subroutine diagvel_t3d(mrl)

  use mem_basic, only: vc, wc, vxe, vye, vze
  implicit none

  integer, intent(in) :: mrl

  call vel_t3d_hex(mrl, vc, wc, vxe, vye, vze)

  return
end subroutine diagvel_t3d

!===============================================================================

subroutine vel_t3d_hex(mrl, vs, ws, vxe, vye, vze)

use mem_ijtabs, only: jtab_w, itab_v, itab_w, jtw_prog
use mem_grid,   only: mza, mva, mwa, lpw, lpv, vnx, vny, vnz, wnx, wny, wnz
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mrl

real, intent(in)  :: vs(mza,mva)
real, intent(in)  :: ws(mza,mwa)

real, intent(out) :: vxe(mza,mwa)
real, intent(out) :: vye(mza,mwa)
real, intent(out) :: vze(mza,mwa)

integer :: j,iw,npoly,kb,k,jv,iv
real    :: farv2, wst

if (mrl == 0) return

! Horizontal loop over W columns

!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,kb,k,jv,iv,wst) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

! Vertical loop over T levels

   do k = kb, mza

! Diagnose 3D earth-velocity vector at T points; W contribution first

      wst = 0.5 * (ws(k-1,iw) + ws(k,iw))

      vxe(k,iw) = wst * wnx(iw)
      vye(k,iw) = wst * wny(iw)
      vze(k,iw) = wst * wnz(iw)

   enddo

! Loop over V neighbors of this W cell

   do jv = 1, npoly

      iv = itab_w(iw)%iv(jv)

! Vertical loop over T levels

      do k = kb, mza

! Diagnose 3D earth-velocity vector at T points; VC contribution

         vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_vx(jv) * vs(k,iv)
         vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_vy(jv) * vs(k,iv)
         vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_vz(jv) * vs(k,iv)

      enddo
      
   enddo
   
   vxe(1:kb-1,iw) = vxe(kb,iw)
   vye(1:kb-1,iw) = vye(kb,iw)
   vze(1:kb-1,iw) = vze(kb,iw)

enddo
!$omp end parallel do

return
end subroutine vel_t3d_hex
