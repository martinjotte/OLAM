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
subroutine fldsisan(o_press, o_rho, o_theta, o_rrw, o_uzonal, o_umerid, o_ozone)

use mem_basic,   only: vmc, vmp, vc, vp, thil, rr_w, rr_v, &
                       wmc, wc, theta, tair, rho, press
use mem_grid,    only: mza, mwa, lpv, lpw, zm, zt, xev, yev, zev, &
                       vnx, vny, vnz
use misc_coms,   only: io6, iparallel
use mem_micro,   only: rr_c, con_c, cldnum
use micro_coms,  only: miclevel, ccnparm, jnmb, rxmin, zfactor_ccn
use mem_ijtabs,  only: jtab_v, jtab_w, itab_v, jtv_init, jtw_init
use consts_coms, only: r8, rocp, t00, p00i, alvlocp, eradi
use var_tables,  only: num_var, vtab_r
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
use obnd,        only: lbcopy_v, lbcopy_w
use therm_lib,   only: rhovsl

implicit none

real(r8), intent(in) :: o_press (mza,mwa)
real(r8), intent(in) :: o_rho   (mza,mwa)
real,     intent(in) :: o_theta (mza,mwa)
real,     intent(in) :: o_rrw   (mza,mwa)
real,     intent(in) :: o_uzonal(mza,mwa)
real,     intent(in) :: o_umerid(mza,mwa)
real,     intent(in) :: o_ozone (mza,mwa)

integer :: j,iw,k,ka,iv,iw1,iw2,mrl,n
real    :: ccn, cond, sh_c
real    :: raxis,raxisi,ug,vg,uvgr,uvgx,uvgy,uvgz

! If initializing the model, fill the main model arrays
! and initialize related arrays

!----------------------------------------------------------------------

!$omp parallel do private(iw,ka,k,cond,sh_c,n,ccn)
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

   ka = lpw(iw)

   wc (1:mza,iw) = 0.0
   wmc(1:mza,iw) = 0.0

   do k = ka, mza

      rho  (k,iw) = o_rho  (k,iw)
      rr_w (k,iw) = o_rrw  (k,iw)
      theta(k,iw) = o_theta(k,iw)
      press(k,iw) = o_press(k,iw)
      tair (k,iw) = o_theta(k,iw) * (real(press(k,iw)) * p00i) ** rocp
      rr_v (k,iw) = o_rrw(k,iw)
      thil (k,iw) = theta(k,iw)


      if (miclevel >= 2) then
         cond = rr_w(k,iw) * real(rho(k,iw)) - rhovsl(tair(k,iw)-t00)

         if (cond > rxmin(1)) then
            rr_c(k,iw) = cond / real(rho(k,iw))
            rr_v(k,iw) = rr_w(k,iw) - rr_c(k,iw)
            sh_c       = rr_c(k,iw) / (1.0 + rr_v(k,iw))
            thil(k,iw) = theta(k,iw) / (1. + alvlocp * sh_c / max(tair(k,iw),253.))
         else
            rr_c(k,iw) = 0.0
         endif
      endif

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

   do k = 1, ka-1
      thil(k,iw) = thil(ka,iw)
   enddo

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
         if (rr_c(k,iw) * real(rho(k,iw)) > rxmin(1)) then
            con_c(k,iw) = ccn * real(rho(k,iw)) * zfactor_ccn(k)
         else
            con_c(k,iw) = 0.0
         endif
      enddo
   endif

enddo
!$omp end parallel do

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
!$omp parallel do private(iv,iw1,iw2,ka,raxis,raxisi,k,ug,vg,uvgr,uvgx,uvgy,uvgz)
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   ka = lpv(iv)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

! Average winds to V point and rotate at V point

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = ka, mza
         ug = .5 * (o_uzonal(k,iw1) + o_uzonal(k,iw2))
         vg = .5 * (o_umerid(k,iw1) + o_umerid(k,iw2))

         uvgr = -vg * zev(iv) * eradi  ! radially outward from axis

         uvgx = (-ug * yev(iv) + uvgr * xev(iv)) * raxisi
         uvgy = ( ug * xev(iv) + uvgr * yev(iv)) * raxisi
         uvgz =   vg * raxis * eradi

         vc (k,iv) = uvgx * vnx(iv) + uvgy * vny(iv) + uvgz * vnz(iv)
         vmc(k,iv) = vc(k,iv) * 0.5 * real(rho(k,iw1) + rho(k,iw2))
      enddo
   else
      vc (ka,iv) = 0.
      vmc(ka,iw) = 0.
   endif

! For below-ground points, set VC to 0

   vc (1:ka-1,iv) = 0.
   vmc(1:ka-1,iv) = 0.

enddo
!$omp end parallel do

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

end subroutine fldsisan
