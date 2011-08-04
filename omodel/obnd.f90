!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine trsets()

use var_tables, only: num_scalar, scalar_tab
use mem_basic,  only: theta, sh_v
use micro_coms, only: level
use misc_coms,  only: io6

implicit none

integer :: n,j,iw,iwp,k

! SCALARS - Lateral, top, and bottom boundary conditions.

! Prognostic scalars other than thil

do n = 1,num_scalar

!      call cyclic(nzp,nxp,nyp,scalarp,'T')
   call latsett(scalar_tab(n)%var_p)
   call topbott(scalar_tab(n)%var_p)

enddo

! THETA and SH_V

call latsett(theta)
call topbott(theta)

if (level >= 1) then
   call latsett(sh_v)
   call topbott(sh_v)
endif

return
end subroutine trsets

!===============================================================================

subroutine latsett(sclr)

use mem_ijtabs, only: jtab_w, itab_w, istp, mrl_endl
use mem_grid,   only: mza, mwa, lpw
use misc_coms,  only: io6

implicit none

real, intent(inout) :: sclr(mza,mwa)

integer :: j,iw,iwp,k,mrl

! if (ipara == 0) then
!    call cyclic (m1,m2,m3,vp,'V')
! endif

! All boundaries - zero gradient condition

!!----------------------------------------------------------------------
!do j = 1,mw_zglbc(mrl); iw = jw_zglbc(j,mrl); iwp = iwa(3,iw)
!!----------------------------------------------------------------------

!   do k = lpw(iw),mza-1
!      sclr(k,iw) = sclr(k,iwp)
!   enddo
!enddo

! Neumann BC for scalars

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
do j = 1,jtab_w(31)%jend(mrl); iw = jtab_w(31)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
      sclr(k,iw) = sclr(k,iwp)
   enddo
enddo
endif
call rsub('W',31)

! Dirichlet BC for scalars

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
do j = 1,jtab_w(32)%jend(mrl); iw = jtab_w(32)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza-1
!      sclr(k,iw) = (remains unchanged)
   enddo
enddo
endif
call rsub('W',32)

return
end subroutine latsett

!===============================================================================

subroutine topbott(sclr)

use mem_ijtabs, only: jtab_w, istp, mrl_endl
use mem_grid,   only: mza, mwa, lpw, dzm, dzim

!$ use omp_lib

implicit none

real, intent(inout) :: sclr(mza,mwa)

integer :: iw,j,k,ka,mrl
real :: dzmr

dzmr = dzm(mza-1) * dzim(mza-2)

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)

if (mrl > 0) then
!$omp parallel do private (iw,ka,k)
do j = 1,jtab_w(33)%jend(mrl); iw = jtab_w(33)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   sclr(mza,iw) = max(0.,  &
      sclr(mza-1,iw) + dzmr * (sclr(mza-1,iw) - sclr(mza-2,iw)))

   ka = lpw(iw)
   do k = ka-1,1,-1
      sclr(k,iw) = sclr(ka,iw)
   enddo
   
enddo
!$omp end parallel do
endif
call rsub('W',33)

return
end subroutine topbott

