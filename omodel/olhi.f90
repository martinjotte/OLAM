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
subroutine fldslhi()

use mem_basic,   only: theta, thil, tair, rho, press, sh_w, sh_v,  &
                       wc, wmc, vc, vp, vmc, vmp
use mem_micro,   only: sh_c
use micro_coms,  only: miclevel
use mem_ijtabs,  only: jtab_w, jtab_v, itab_v, &
                       jtv_init, jtw_init
use consts_coms, only: p00, p00i, rocp, cvocp, p00k, rdry, rvap, alvl, cp, grav
use mem_grid,    only: mza, mva, mwa, lpv, lpw, zt, zm, dzt_top, dzt_bot, &
                       xev, yev, zev, unx, uny, vnx, vny, glatw, arv
use mem_zonavg,  only: zonz_vect, zonu_vect, zont_vect, zonr_vect,  &
                       zonp_vect, zonz, zonu, zont, zonr, zonavg_init
use misc_coms,   only: io6, iparallel, idate1, imonth1, iyear1
use therm_lib,   only: rhovsl

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w,  &
                        mpi_send_v, mpi_recv_v 
    
use obnd,         only: lbcopy_v, lbcopy_w

implicit none   

integer :: j,iw,k,ka,iu,iv,iter,iw1,iw2,iup,ivp,llat,mrl
integer :: im,iplev,ilat,im1,im2,ilatn,ilats
integer :: kpbc, kbc, kother, khi, klo

real :: rcloud,temp,rvls,exner
real :: ugx,ugy,raxis,raxisi,pnorth,psouth,ug,ug1,ug2,wt1
real :: alat,dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,rlat,wt2,cpo2g,rhovs
real :: extrap, pkhyd, pressnew, qhydm

character(len=1) :: line

real :: vctr1(mza), vctr2(mza), vctr3(mza)
real :: uzonal(mza,mwa)

! Choose as an internal pressure boundary condition the pressure level
! zonp_vect(3), whose pressure is 46415.89 Pa.

kpbc = 3

! Fill zonavg arrays for initialization time

call zonavg_init(idate1,imonth1,iyear1)

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!----------------------------------------------------------------------

! Linearly interpolate zonavg arrays by latitude to current IW column 
! and K level

   ka = lpw(iw)

   rlat = .4 * (glatw(iw) + 93.75)
   ilat = int(rlat)
   wt2 = rlat - float(ilat)
   
   zonz_vect(1:22) = (1. - wt2) * zonz(ilat,1:22) + wt2 * zonz(ilat+1,1:22)
   zont_vect(1:22) = (1. - wt2) * zont(ilat,1:22) + wt2 * zont(ilat+1,1:22)
   zonu_vect(1:22) = (1. - wt2) * zonu(ilat,1:22) + wt2 * zonu(ilat+1,1:22)
   zonr_vect(1:22) = (1. - wt2) * zonr(ilat,1:22) + wt2 * zonr(ilat+1,1:22)
   
! Interpolate zonavg vector arrays in height to model levels

! Fill pressure, theta, air density, and vapor density arrays from zonavg
! vector arrays prior to iteration

   call hintrp_ee(22,zonp_vect,zonz_vect,mza,vctr1,zt) ! pressure
   call hintrp_ee(22,zont_vect,zonz_vect,mza,vctr2,zt) ! temp
   call hintrp_ee(22,zonr_vect,zonz_vect,mza,vctr3,zt) ! specific humidity
   call hintrp_ee(22,zonu_vect,zonz_vect,mza,uzonal(1,iw),zt) ! uzonal

   press(1:mza,iw) = vctr1(1:mza)
   theta(1:mza,iw) = vctr2(1:mza) * (p00 / vctr1(1:mza)) ** rocp  ! temp to theta
   thil(1:mza,iw) = theta(1:mza,iw)
   rho(1:mza,iw)  = press(1:mza,iw) ** cvocp * p00k / (rdry * theta(1:mza,iw))
   sh_v(1:mza,iw) = vctr3(1:mza)  ! spec hum
   sh_w(1:mza,iw) = sh_v(1:mza,iw)
   wc(1:mza,iw)   = 0.
   wmc(1:mza,iw)  = 0.

! Determine which two model zt levels bracket zonz_vect(kpbc) in this column

   khi = 2
   do while (zt(khi) < zonz_vect(kpbc))
      khi = khi + 1
   enddo
   klo = khi - 1

! Determine whether zt(klo) or zt(khi) is closer to zonz_vect(kpbc).  The closer
! one will undergo direct pressure adjustment to satisfy the internal b.c.

   if (zt(khi) - zonz_vect(kpbc) < zonz_vect(kpbc) - zt(klo)) then
      kbc = khi
      kother = klo
   else
      kbc = klo
      kother = khi
   endif

   extrap = (zt(kbc) - zt(kother)) / (zonz_vect(kpbc) - zt(kother))

! Iterative hydrostatic integration

   do iter = 1,100
   
! Adjust pressure at k = kbc.  Use temporal weighting for damping

      pressnew = press(kother,iw) * (zonp_vect(kpbc) / press(kother,iw)) ** extrap
      press(kbc,iw) = .1 * press(kbc,iw) + .9 * pressnew

!  Compute density for all levels

      do k = ka,mza

! Try this: hold Mclatchy temp (vctr2) constant during iterations
      
         theta(k,iw) = vctr2(k) * (p00 / press(k,iw)) ** rocp
         thil(k,iw) = theta(k,iw)

         if (miclevel == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (miclevel == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) * p00i) ** rocp  ! Defined WITHOUT CP factor
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)

            sh_c(k,iw) = sh_w(k,iw) - rhovs/real(rho(k,iw))
            sh_c(k,iw) = max(0.,sh_c(k,iw))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)

            rho(k,iw) = press(k,iw) ** cvocp * p00k &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))

            qhydm = alvl * sh_c(k,iw)
            thil(k,iw) = theta(k,iw) / (1. + qhydm / (cp * max(temp,253.)))
         endif

      enddo

! Integrate hydrostatic equation upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations

      do k = kbc+1,mza
         pkhyd = press(k-1,iw) &
               - grav * (rho(k-1,iw) * dzt_top(k-1) + rho(k,iw) * dzt_bot(k))
         press(k,iw) = .05 * press(k,iw) + .95 * max(.1,pkhyd)
      enddo

      do k = kbc-1,ka,-1
         pkhyd = press(k+1,iw) &
               + grav * (rho(k+1,iw) * dzt_bot(k+1) + rho(k,iw) * dzt_top(k))
         press(k,iw) = .05 * press(k,iw) + .95 * max(.1,pkhyd)
      enddo

   enddo

   do k = ka, mza
      tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
   enddo

   do k = 1, ka-1
      thil(k,iw) = thil(ka,iw)
   enddo

enddo

! LBC copy (THETA and TAIR will be copied later with the scalars)

mrl = 1

if (iparallel == 1) then
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)

   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc, rvara2=wmc, rvara3=thil)

   call mpi_send_w(mrl, rvara1=uzonal)
   call mpi_recv_w(mrl, rvara1=uzonal)
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

! Initialize VMC, VC

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   ka = lpv(iv)

   raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

! Average winds to V point and rotate at V point (assumed to be global simulation)

   if (raxis > 1.e3) then
      raxisi = 1. / raxis

      do k = ka,mza
         ug = .5 * (uzonal(k,iw1) + uzonal(k,iw2))

         ugx = -ug * yev(iv) * raxisi 
         ugy =  ug * xev(iv) * raxisi

         vc(k,iv) = ugx * vnx(iv) + ugy * vny(iv)
      enddo

   else
      vc(ka:mza,iv) = 0.
   endif

   vmc(ka:mza,iv) = vc(ka:mza,iv) * .5 * (rho(ka:mza,iw1) + rho(ka:mza,iw2))

! For below-ground points, set VC to 0.

   vc (1:ka-1,iv) = 0.0
   vmc(1:ka-1,iv) = 0.0

enddo

! MPI parallel send/recv of V group

if (iparallel == 1) then
   call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
endif

! LBC copy of VMC, VC

call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

if (allocated(vmp)) vmp(:,:) = vmc(:,:)
if (allocated(vp )) vp (:,:) = vc (:,:)

! print out initial state from 1st jtw_init column

iw = jtab_w(jtw_init)%iw(1)

write(io6,*) ' '
write(io6,*) '========================================================================='
write(io6,*) '                    OLAM INITIAL STATE COLUMN (lhi)'
write(io6,*) '========================================================================='
write(io6,*) '   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*) '========================================================================='
write(io6,*) ' '

do k = mza,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
enddo
   
write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

return
end subroutine fldslhi
