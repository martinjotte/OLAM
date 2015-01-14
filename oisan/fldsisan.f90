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
subroutine fldsisan(o_rho, o_theta, o_shv, o_vc, o_ozone)

use mem_basic,   only: vmc, vmp, vc, vp, thil, sh_w, sh_v, &
                       wmc, wc, theta, tair, rho, press
use mem_grid,    only: mza, mva, mwa, lpv, lpw, zm, zt
use misc_coms,   only: io6, deltax, iparallel, runtype
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_v, jtab_w, itab_v, itab_w, jtv_init, jtw_init
use consts_coms, only: r8, pc1, rdry, rvap, cpocv, rocp, p00i
use var_tables,  only: num_var, vtab_r
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

use obnd,         only: lbcopy_v, lbcopy_w

implicit none

real(r8), intent(in) :: o_rho   (mza,mwa)
real,     intent(in) :: o_theta (mza,mwa)
real,     intent(in) :: o_shv   (mza,mwa)
real,     intent(in) :: o_vc    (mza,mva)
real,     intent(in) :: o_ozone (mza,mwa)

integer :: j,iw,k,ka,iv,iw1,iw2,mrl,n
real    :: alph_p

! If initializing the model, fill the main model arrays
! and initialize related arrays

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

   ka = lpw(iw)

   do k = ka,mza

      rho  (k,iw) = o_rho(k,iw)
      theta(k,iw) = o_theta(k,iw)
      sh_v (k,iw) = o_shv(k,iw)

      thil(k,iw) = theta(k,iw)
      sh_w(k,iw) = sh_v(k,iw)   ! no condensate considered here
      if (level > 1) then
         sh_c(k,iw) = 0.        ! no condensate considered here
      endif

! alph_p is like alpha_press except that the (theta/thil) factor is excluded
! here (it's equal to 1)

      alph_p = pc1 * (((1. - sh_w(k,iw)) * rdry + sh_v(k,iw) * rvap)) ** cpocv

      press(k,iw) = alph_p * (rho(k,iw) * thil(k,iw)) ** cpocv

      tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp

      wc(k,iw) = 0.
      wmc(k,iw) = 0.

   enddo

   ! Initialize any variable in the var_table that is named ozone.
   ! Assumed units are ppmv

   do n = 1, num_var
      if ( (vtab_r(n)%name == 'O3'   ) .or. (vtab_r(n)%name == 'o3'   ) .or. &
           (vtab_r(n)%name == 'OZONE') .or. (vtab_r(n)%name == 'ozone') ) then

         if ( associated( vtab_r(n)%rvar2_p ) ) then
            vtab_r(n)%rvar2_p(ka:mza,iw) = o_ozone(ka:mza,iw)
         endif
      endif
   enddo

enddo

! LBC copy (THETA and TAIR and OZONE will be copied later with the scalars)

if (iparallel == 1) then
   mrl = 1
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

! Initialize VMC, VC

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   ka = lpv(iv)

   do k = ka,mza
      vc(k,iv) = o_vc(k,iv)
      vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
   enddo

! For below-ground points, set VC to LPV value.

   vc(1:ka-1,iv) = vc(ka,iv)

enddo
   
! MPI parallel send/recv of V group

if (iparallel == 1) then
   mrl = 1
   call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
endif

! LBC copy of VMC, VC

call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

if (allocated(vmp)) vmp(:,:) = vmc(:,:)
if (allocated(vp )) vp (:,:) = vc (:,:)

! Print out initial state from 1st jtw_init column

iw = jtab_w(jtw_init)%iw(1)

write(io6,*)' '
write(io6,*)'========================================================================'
write(io6,*)'                    OLAM INITIAL STATE COLUMN (vari)'
write(io6,*)'========================================================================'
write(io6,*)'   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*)'========================================================================'
write(io6,*)' '

do k = mza,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)') &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
enddo

write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*)' '

return
end subroutine fldsisan
