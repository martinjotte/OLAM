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
Module obnd

Contains

subroutine trsets(mrl)

use var_tables, only: num_scalar, scalar_tab
use mem_basic,  only: theta, tair, rr_v
use micro_coms, only: miclevel

implicit none

integer, intent(in) :: mrl
integer             :: n

! SCALARS - Lateral, top, and bottom boundary conditions.

if (mrl < 1) return

! Prognostic scalars

do n = 1, num_scalar
   call latsett(mrl,scalar_tab(n)%var_p)
   call botset (mrl,scalar_tab(n)%var_p)
enddo

! THETA, TAIR, and RR_V

call latsett(mrl,theta)
call botset (mrl,theta)

call latsett(mrl,tair)
call botset (mrl,tair)

if (miclevel >= 1) then
   call latsett(mrl,rr_v)
   call botset (mrl,rr_v)
endif

end subroutine trsets

!===============================================================================

subroutine latsett(mrl,sclr)

use mem_ijtabs, only: jtab_w, itab_w, jtw_lbcp
use mem_grid,   only: mza, mwa, lpw

implicit none

integer, intent(in)    :: mrl
real,    intent(inout) :: sclr(mza,mwa)

integer :: j,iw,iwp,k

! LBC for scalars (usually cyclic)

if (mrl < 1) return
if (jtab_w(jtw_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
   do k = lpw(iw),mza
      sclr(k,iw) = sclr(k,iwp)
   enddo
enddo

end subroutine latsett

!===============================================================================

subroutine botset(mrl,sclr)

use mem_ijtabs, only: jtab_w, jtw_prog
use mem_grid,   only: mza, mwa, lpw

implicit none

integer, intent(in)    :: mrl
real,    intent(inout) :: sclr(mza,mwa)

integer :: iw,j,k,ka

! Top/bottom boundary condition for scalars: zero-gradient

if (mrl < 1) return

!----------------------------------------------------------------------
!$omp parallel do private (iw,ka,k)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   ka = lpw(iw)
   do k = ka-1,1,-1
      sclr(k,iw) = sclr(ka,iw)
   enddo

enddo
!$omp end parallel do

end subroutine botset

!===============================================================================

subroutine set_bottom()

  use mem_ijtabs, only: jtab_w, jtw_prog
  use mem_grid,   only: lpw
  use var_tables, only: nvar_par, vtab_r, nptonv
  use mem_basic,  only: thil

  implicit none

  integer :: j, iw, ka, n, k

  !$omp parallel do private (iw,ka,n,k)
  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     ka = lpw(iw)

     do n = 1, nvar_par
        do k = 1, ka-1
           vtab_r(nptonv(n))%rvar2_p(k,iw) = vtab_r(nptonv(n))%rvar2_p(ka,iw)
        enddo
     enddo

     do k = 1, ka-1
        thil(k,iw) = thil(ka,iw)
     enddo

  enddo
  !$omp end parallel do

end subroutine set_bottom

!===============================================================================

subroutine lbcopy_m(mrl, a1)

use mem_ijtabs, only: jtab_m, itab_m, jtm_lbcp
use mem_grid,   only: mza, mma

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1(mza,mma)

integer :: j,im,imp

! Lateral boundary copy (usually cyclic)

if (mrl < 1) return
if (jtab_m(jtm_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
do j = 1,jtab_m(jtm_lbcp)%jend(mrl); im = jtab_m(jtm_lbcp)%im(j)
   imp = itab_m(im)%imp
!----------------------------------------------------------------------

   if (present(a1)) a1(:,im) = a1(:,imp) 

enddo

end subroutine lbcopy_m

!===============================================================================

subroutine lbcopy_v(mrl, vmc, vc, iv1, iv2, iv3)

use mem_ijtabs, only: jtab_v, itab_v, jtv_lbcp
use mem_grid,   only: mza, mva

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: vmc(mza,mva)
real, optional, intent(inout) :: vc (mza,mva)

integer, optional, intent(inout) :: iv1(mva)
integer, optional, intent(inout) :: iv2(mva)
integer, optional, intent(inout) :: iv3(mva)

integer :: j,iv,ivp

! Lateral boundary copy (usually cyclic)

if (mrl < 1) return
if (jtab_v(jtv_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
!$omp parallel
!$omp do private(iv,ivp)
do j = 1,jtab_v(jtv_lbcp)%jend(mrl); iv = jtab_v(jtv_lbcp)%iv(j)
   ivp = itab_v(iv)%ivp
!----------------------------------------------------------------------

   if (present(vmc)) vmc(:,iv) = vmc(:,ivp)
   if (present(vc))  vc (:,iv) = vc (:,ivp)

   if (present(iv1)) iv1(iv) = iv1(ivp)
   if (present(iv2)) iv2(iv) = iv2(ivp)
   if (present(iv3)) iv3(iv) = iv3(ivp)
enddo
!$omp end do
!$omp end parallel

end subroutine lbcopy_v

!===============================================================================

subroutine lbcopy_w(mrl, a1,  a2,  a3,  a4,  a5,  a6,  a7,  a8,  a9,  a10,   &
                         a11, a12, a13, a14, a15, a16, a17, a18, a19, a20,   &
                         d1,  d2,  s1,  s2,  v1,  v2,  v3,  v4, iv1, iv2, iv3)

use mem_ijtabs,  only: jtab_w, itab_w, jtw_lbcp
use mem_grid,    only: mza, mwa
use consts_coms, only: r8

implicit none

integer, intent(in) :: mrl

! For real 2D arrays dimensioned to mza
real, optional, intent(inout) :: a1 (mza,mwa)
real, optional, intent(inout) :: a2 (mza,mwa)
real, optional, intent(inout) :: a3 (mza,mwa)
real, optional, intent(inout) :: a4 (mza,mwa)
real, optional, intent(inout) :: a5 (mza,mwa)
real, optional, intent(inout) :: a6 (mza,mwa)
real, optional, intent(inout) :: a7 (mza,mwa)
real, optional, intent(inout) :: a8 (mza,mwa)
real, optional, intent(inout) :: a9 (mza,mwa)
real, optional, intent(inout) :: a10(mza,mwa)
real, optional, intent(inout) :: a11(mza,mwa)
real, optional, intent(inout) :: a12(mza,mwa)
real, optional, intent(inout) :: a13(mza,mwa)
real, optional, intent(inout) :: a14(mza,mwa)
real, optional, intent(inout) :: a15(mza,mwa)
real, optional, intent(inout) :: a16(mza,mwa)
real, optional, intent(inout) :: a17(mza,mwa)
real, optional, intent(inout) :: a18(mza,mwa)
real, optional, intent(inout) :: a19(mza,mwa)
real, optional, intent(inout) :: a20(mza,mwa)

! For real*8 2D arrays dimensioned to mza
real(8), optional, intent(inout) :: d1(mza,mwa)
real(8), optional, intent(inout) :: d2(mza,mwa)

! 1D real vectors
real, optional, intent(inout) :: v1(mwa)
real, optional, intent(inout) :: v2(mwa)
real, optional, intent(inout) :: v3(mwa)
real, optional, intent(inout) :: v4(mwa)

! 1D integer vectors
integer, optional, intent(inout) :: iv1(mwa)
integer, optional, intent(inout) :: iv2(mwa)
integer, optional, intent(inout) :: iv3(mwa)

! For real 2D arrays where the size of the first dimension is not mza
real, optional, contiguous, intent(inout) :: s1(:,:)
real, optional, contiguous, intent(inout) :: s2(:,:)

integer :: j,iw,iwp

! Lateral boundary copy (usually cyclic)

if (mrl < 1) return
if (jtab_w(jtw_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
!$omp parallel
!$omp do private(iw,iwp)
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

   if (present(a1 )) a1 (:,iw) = a1 (:,iwp)
   if (present(a2 )) a2 (:,iw) = a2 (:,iwp)
   if (present(a3 )) a3 (:,iw) = a3 (:,iwp)
   if (present(a4 )) a4 (:,iw) = a4 (:,iwp)
   if (present(a5 )) a5 (:,iw) = a5 (:,iwp)
   if (present(a6 )) a6 (:,iw) = a6 (:,iwp)
   if (present(a7 )) a7 (:,iw) = a7 (:,iwp)
   if (present(a8 )) a8 (:,iw) = a8 (:,iwp)
   if (present(a9 )) a9 (:,iw) = a9 (:,iwp)
   if (present(a10)) a10(:,iw) = a10(:,iwp)
   if (present(a11)) a11(:,iw) = a11(:,iwp)
   if (present(a12)) a12(:,iw) = a12(:,iwp)
   if (present(a13)) a13(:,iw) = a13(:,iwp)
   if (present(a14)) a14(:,iw) = a14(:,iwp)
   if (present(a15)) a15(:,iw) = a15(:,iwp)
   if (present(a16)) a16(:,iw) = a16(:,iwp)
   if (present(a17)) a17(:,iw) = a17(:,iwp)
   if (present(a18)) a18(:,iw) = a18(:,iwp)
   if (present(a19)) a19(:,iw) = a19(:,iwp)
   if (present(a20)) a20(:,iw) = a20(:,iwp)

   if (present(v1 )) v1   (iw) = v1   (iwp)
   if (present(v2 )) v2   (iw) = v2   (iwp)
   if (present(v3 )) v3   (iw) = v3   (iwp)
   if (present(v4 )) v4   (iw) = v4   (iwp)

   if (present(s1 )) s1 (:,iw) = s1 (:,iwp)
   if (present(s2 )) s2 (:,iw) = s2 (:,iwp)

   if (present(d1 )) d1 (:,iw) = d1 (:,iwp)
   if (present(d2 )) d2 (:,iw) = d2 (:,iwp)

   if (present(iv1)) iv1  (iw) = iv1  (iwp)
   if (present(iv2)) iv2  (iw) = iv2  (iwp)
   if (present(iv3)) iv3  (iw) = iv3  (iwp)

enddo
!$omp end do
!$omp end parallel

end subroutine lbcopy_w

!===============================================================================

subroutine lbcopy_wex(mrl, n1, a1)

use mem_ijtabs, only: jtab_w, itab_w, jtw_lbcp
use mem_grid,   only: mza, mwa

implicit none

integer, intent(in)    :: mrl, n1
real,    intent(inout) :: a1(n1,mwa)
integer                :: j, iw, iwp

! Lateral boundary copy (usually cyclic)

if (mrl < 1) return
if (jtab_w(jtw_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
!$omp parallel
!$omp do private(iw,iwp)
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

   a1(:,iw) = a1(:,iwp)

enddo
!$omp end do
!$omp end parallel

end subroutine lbcopy_wex

!===============================================================================

subroutine lbcopy_w1d(mrl, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, d1)

use mem_ijtabs,  only: jtab_w, itab_w, jtw_lbcp
use mem_grid,    only: mwa
use consts_coms, only: r8

implicit none

integer, intent(in) :: mrl

real, optional, intent(inout) :: a1 (mwa)
real, optional, intent(inout) :: a2 (mwa)
real, optional, intent(inout) :: a3 (mwa)
real, optional, intent(inout) :: a4 (mwa)
real, optional, intent(inout) :: a5 (mwa)
real, optional, intent(inout) :: a6 (mwa)
real, optional, intent(inout) :: a7 (mwa)
real, optional, intent(inout) :: a8 (mwa)
real, optional, intent(inout) :: a9 (mwa)
real, optional, intent(inout) :: a10(mwa)
real, optional, intent(inout) :: a11(mwa)
real, optional, intent(inout) :: a12(mwa)

real(r8), optional, intent(inout) :: d1(mwa)

integer :: j,iw,iwp

! Lateral boundary copy (usually cyclic)

if (mrl < 1) return
if (jtab_w(jtw_lbcp)%jend(mrl) < 1) return

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_lbcp)%jend(mrl); iw = jtab_w(jtw_lbcp)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------

   if (present(a1 )) a1 (iw) = a1 (iwp) 
   if (present(a2 )) a2 (iw) = a2 (iwp) 
   if (present(a3 )) a3 (iw) = a3 (iwp) 
   if (present(a4 )) a4 (iw) = a4 (iwp) 
   if (present(a5 )) a5 (iw) = a5 (iwp) 
   if (present(a6 )) a6 (iw) = a6 (iwp) 
   if (present(a7 )) a7 (iw) = a7 (iwp) 
   if (present(a8 )) a8 (iw) = a8 (iwp) 
   if (present(a9 )) a9 (iw) = a9 (iwp) 
   if (present(a10)) a10(iw) = a10(iwp) 
   if (present(a11)) a11(iw) = a11(iwp) 
   if (present(a12)) a12(iw) = a12(iwp) 
   if (present(d1 )) d1 (iw) = d1 (iwp) 

enddo

end subroutine lbcopy_w1d

End Module obnd
