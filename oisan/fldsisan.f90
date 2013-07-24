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
subroutine fldsisan(o_rho, o_theta, o_shv, o_uvc, o_ozone)

use mem_basic,   only: umc, ump, uc, vmc, vmp, vc, vp, thil, sh_w, sh_v, &
                       wmc, wc, theta, tair, rho, press
use mem_grid,    only: mza, mua, mva, mwa, lcu, lpv, lpw, zm, zt
use misc_coms,   only: io6, deltax, iparallel, runtype, meshtype, do_chem
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_u, jtab_v, jtab_w, itab_u, itab_v, itab_w, &
                       jtu_init, jtv_init, jtw_init
use consts_coms, only: pc1, rdry, rvap, cpocv, rocp, p00i

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                        mpi_send_u, mpi_recv_u, &
                        mpi_send_v, mpi_recv_v 

use obnd,         only: lbcopy_u, lbcopy_v, lbcopy_w

use cgrid_defn,   only: ns_o3
use var_tables,   only: scalar_tab

implicit none

real(kind=8), intent(in) :: o_rho   (mza,mwa)
real,         intent(in) :: o_theta (mza,mwa)
real,         intent(in) :: o_shv   (mza,mwa)
real,         intent(in) :: o_uvc   (mza,mva)
real,         intent(in) :: o_ozone (mza,mwa)

integer :: j,iw,k,ka,iu,iv,iw1,iw2,iup,ivp

real :: rcloud,temp,rvls,uu,vv
real :: alph_p

! If initializing the model, fill the main model arrays
! and initialize related arrays

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------
call qsub('W',iw)
   do k = 1,mza

      rho  (k,iw) = o_rho(k,iw)
      theta(k,iw) = o_theta(k,iw)
      sh_v (k,iw) = o_shv(k,iw)

      thil(k,iw) = theta(k,iw)
      sh_w(k,iw) = sh_v(k,iw)   ! no condensate considered here
      if (level > 1) then
         sh_c(k,iw) = 0.        ! no condensate considered here
      endif

      if (do_chem) then
         scalar_tab(ns_o3)%var_p(k,iw) = o_ozone(k,iw)
      endif

! alph_p is like alpha_press except that the (theta/thil) factor is excluded
! here (it's equal to 1)

      alph_p = pc1 * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)) ** cpocv

      press(k,iw) = alph_p * (rho(k,iw) * thil(k,iw)) ** cpocv

      tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp

      wc(k,iw) = 0.
      wmc(k,iw) = 0.

   enddo
enddo
call rsub('Wb',7)

! LBC copy (THETA and TAIR will be copied later with the scalars)

if (iparallel == 1) then
   call mpi_send_w('I')  ! Send W group
   call mpi_recv_w('I')  ! Recv W group
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

if (meshtype == 1) then

! If triangular mesh, initialize UMC, UC

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(jtu_init)%jend(1); iu = jtab_u(jtu_init)%iu(j)
      iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      do k = 1,mza
         uc(k,iu) = o_uvc(k,iu)
         umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
      enddo

! For below-ground points, set UC to LCU value.

      ka = lcu(iu)
      uc(1:ka-1,iu) = uc(ka,iu)

   enddo
   
   call rsub('Ub',7)

! MPI parallel send/recv of U group

   if (iparallel == 1) then
      call mpi_send_u('I')
      call mpi_recv_u('I')
   endif

! LBC copy of UMC, UC

   call lbcopy_u(1, a1=umc, a2=uc)

! Set UMP to UMC

   ump(:,:) = umc(:,:)

else

! If using hexagonal mesh, initialize VMC, VC

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
      iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   call qsub('V',iv)

      do k = 1,mza
         vc(k,iv) = o_uvc(k,iv)
         vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
      enddo

! For below-ground points, set VC to LPV value.

      ka = lpv(iv)
      vc(1:ka-1,iv) = vc(ka,iv)

   enddo
   
   call rsub('Vb',7)

! MPI parallel send/recv of V group

   if (iparallel == 1) then
      call mpi_send_v('I')
      call mpi_recv_v('I')
   endif

! LBC copy of VMC, VC

   call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

   if (allocated(vmp)) vmp(:,:) = vmc(:,:)
   if (allocated(vp )) vp (:,:) = vc (:,:)

endif

! Print out initial state from 1st jtw_init column

iw = jtab_w(jtw_init)%iw(1)

write(io6,*)' '
write(io6,*)'========================================================================'
write(io6,*)'                    OLAM INITIAL STATE COLUMN (vari)'
write(io6,*)'========================================================================'
write(io6,*)'   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*)'========================================================================'
write(io6,*)' '

do k = mza-1,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)') &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
enddo

write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*)' '

return
end subroutine fldsisan
