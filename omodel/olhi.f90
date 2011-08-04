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
subroutine fldslhi()

use mem_basic,   only: theta, thil, rho, press, sh_w, sh_v,  &
                       wc, wmc, uc, umc, ump, vc, vmc, vmp
use mem_micro,   only: sh_c
use micro_coms,  only: level
use mem_ijtabs,  only: jtab_w, jtab_u, itab_u, jtab_v, itab_v
use consts_coms, only: p00, rocp, cvocp, p00k, rdry, rvap, alvlocp, gravo2
use mem_grid,    only: mza, mua, mva, mwa, lcu, lcv, zt, dzt, zm, &
                       xeu, yeu, zeu, xev, yev, zev, unx, uny, vnx, vny, &
                       glatw, aru, arv
use mem_zonavg,  only: zonz_vect, zonu_vect, zont_vect, zonr_vect,  &
                       zonp_vect, zonz, zonu, zont, zonr
use misc_coms,   only: io6, iparallel, meshtype
use massflux,    only: diagnose_uc, diagnose_vc

use olam_mpi_atm, only: mpi_send_w, mpi_recv_w,  &
                        mpi_send_u, mpi_recv_u,  &
                        mpi_send_v, mpi_recv_v 
    
implicit none   

integer :: j,iw,k,ka,iu,iv,iter,iw1,iw2,iup,ivp,llat  &
   ,im,iplev,ilat,im1,im2,ilatn,ilats
real :: rcloud,temp,rvls,exner  &
   ,ugx,ugy,raxis,pnorth,psouth,ug,ug1,ug2,wt1,pkhyd  &
   ,alat,dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,rlat,wt2,cpo2g,rhovs
character(len=1) :: line

real, dimension(mza) :: vctr1,vctr2,vctr3  ! automatic arrays
real, dimension(mza,mwa) :: uzonal  ! automatic array

real, external :: rhovsl

call zonavg_init()  ! Read in 'ZONAVG_CLIMATE' data and fill zonavg arrarys

call psub()
!----------------------------------------------------------------------
do j = 1,jtab_w(8)%jend(1); iw = jtab_w(8)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

! Linearly interpolate zonavg arrays by latitude to current IW column 
! and K level

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

! Iterative hydrostatic integration

   do iter = 1,100
   
      do k = 1,mza

! Try this: hold Mclatchy temp (vctr2) constant during iterations
      
         theta(k,iw) = vctr2(k) * (p00 / press(k,iw)) ** rocp
         thil(k,iw) = theta(k,iw)

         if (level == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (level == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) / p00) ** rocp   ! Defined WITHOUT CP factor
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)
            sh_c(k,iw) = max(0.,sh_w(k,iw)-rhovs/real(rho(k,iw)))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)

! SPECIAL: Remove any initial cloud water by setting sh_w to sh_v
!            sh_w(k,iw) = sh_v(k,iw)
! END SPECIAL

! As is done for iteration in sub satadjst, use (0.3,0.7) weights to damp oscil.
!            rho(k,iw) = .5 * rho(k,iw)  &
!                      + .5 * press(k,iw) ** cvocp * p00k  &
            rho(k,iw) = press(k,iw) ** cvocp * p00k  &
               / (theta(k,iw) * (1. - sh_c(k,iw))  &
               * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
            thil(k,iw) = theta(k,iw)  &
               / (1. + alvlocp * sh_c(k,iw) / max(temp,253.))
         endif

         if (k >= 2) then
            pkhyd = press(k-1,iw)  &
               - gravo2 * (rho(k-1,iw) * dzt(k-1) + rho(k,iw) * dzt(k))
! Impose minimum value of .01 Pa to avoid overshoot to negative values 
! during iteration.  Use weighting to damp oscillations
            press(k,iw) = .15 * press(k,iw) + .85 * max(.01,pkhyd)

         endif

      enddo
      
   enddo

enddo
call rsub('Wb',8)

if (iparallel == 1) then
   call mpi_send_w('I')  ! Send W group
   call mpi_recv_w('I')  ! Recv W group
endif

if (meshtype == 1) then

! For triangle mesh, initialize UMC, UC

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_u(8)%jend(1); iu = jtab_u(8)%iu(j)
      iw1 = itab_u(iu)%iw(1); iw2 = itab_u(iu)%iw(2)
!----------------------------------------------------------------------
   call qsub('U',iu)

      if (iw1 == 1) iw1 = iw2
      if (iw2 == 1) iw2 = iw1

! Average winds to U point and rotate at U point (assumed to be global simulation)

      do k = 1,mza
         ug = .5 * (uzonal(k,iw1) + uzonal(k,iw2))

         raxis = sqrt(xeu(iu) ** 2 + yeu(iu) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            ugx = -ug * yeu(iu) / raxis 
            ugy =  ug * xeu(iu) / raxis 

            uc(k,iu) = ugx * unx(iu) + ugy * uny(iu)
         else
            uc(k,iu) = 0.
         endif

         umc(k,iu) = uc(k,iu) * .5 * (rho(k,iw1) + rho(k,iw2))
      enddo

   enddo
   call rsub('Ub',8)

! Set UMC = 0 wherever ARU = 0.

   call psub()
!----------------------------------------------------------------------
   do iu = 2,mua
!----------------------------------------------------------------------
   call qsub('U',iu)
!x      do k = 1,mza-1
!x         if (aru(k,iu) < 1.e-9) then
!x            umc(k,iu) = 0.
!x         endif
!x      enddo

! For below-ground points, set UC to LCU value.

      ka = lcu(iu)
      uc(1:ka-1,iu) = uc(ka,iu)

   enddo
   call rsub('Ub',0)

! MPI parallel send/recv of U group

      if (iparallel == 1) then
         call mpi_send_u('I')
         call mpi_recv_u('I')
      endif

! Set UMP to UMC

   ump(:,:) = umc(:,:)

! Diagnose VC

   call diagnose_vc()

else

! For hexagon grid, initialize VMC, VC

   call psub()
!----------------------------------------------------------------------
   do j = 1,jtab_v(8)%jend(1); iv = jtab_v(8)%iv(j)
      iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------
   call qsub('V',iv)

      if (iw1 == 1) iw1 = iw2
      if (iw2 == 1) iw2 = iw1

! Average winds to V point and rotate at V point (assumed to be global simulation)

      do k = 1,mza
         ug = .5 * (uzonal(k,iw1) + uzonal(k,iw2))

         raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            ugx = -ug * yev(iv) / raxis 
            ugy =  ug * xev(iv) / raxis 

            vc(k,iv) = ugx * vnx(iv) + ugy * vny(iv)
         else
            vc(k,iv) = 0.
         endif
         
         vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
      enddo

   enddo
   call rsub('Vb',8)

! Set VMC = 0 wherever ARV = 0.

   call psub()
!----------------------------------------------------------------------
   do iv = 2,mva
!----------------------------------------------------------------------
   call qsub('V',iv)
!x      do k = 1,mza-1
!x         if (arv(k,iv) < 1.e-9) then
!x            vmc(k,iv) = 0.
!x         endif
!x      enddo

! For below-ground points, set VC to LCV value.

      ka = lcv(iv)
      vc(1:ka-1,iv) = vc(ka,iv)

   enddo
   call rsub('Vb',0)

! MPI parallel send/recv of V group

   if (iparallel == 1) then
      call mpi_send_v('I')  ! Send V group
      call mpi_recv_v('I')  ! Recv V group
   endif

! Set VMP to VMC

   vmp(:,:) = vmc(:,:)

! Diagnose UC

   call diagnose_uc()

endif

! print out initial state column from column 2

iw = 2

write(io6,*) ' '
write(io6,*) '========================================================================='
write(io6,*) '                    OLAM INITIAL STATE COLUMN (lhi)'
write(io6,*) '========================================================================='
write(io6,*) '   zm(m)    k     zt(m)   press(Pa)   rho(kg/m3)   theta(K)   sh_w(g/kg)'
write(io6,*) '========================================================================='
write(io6,*) ' '

do k = mza-1,2,-1
   write(io6, '(f10.2,1x,9(''-------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f11.2,f11.4,f12.2,f11.4)')  &
       k,zt(k),press(k,iw),rho(k,iw),theta(k,iw),sh_w(k,iw)*1.e3
enddo
   
write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

deallocate(zont,zonu,zonz,zonr  &
   ,zonp_vect,zont_vect,zonu_vect,zonz_vect,zonr_vect)

return
end subroutine fldslhi

!===============================================================================

subroutine zonavg_init()

use misc_coms,   only: io6, idate1, imonth1, zonclim                    
use consts_coms, only: pio180, cp, rocp, eps_virt, omega, grav2, dlat
use mem_zonavg,  only: zonz, zonu, zont, zonr, zonp_vect,  &
                       alloc_zonavg, fill_zonavg

implicit none   

integer :: iplev,ilat,ilatn,ilats
real :: dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,cpo2g

! Read in 'ZONAVG_CLIMATE' dataset

call alloc_zonavg()
call fill_zonavg(zonclim,idate1,imonth1)

! In order to obtain zonally averaged water vapor mixing ratio, interpolate
! Mclatchy soundings in time to mclat array and obtain coefficients for 
! latitudinal spline interpolation

call fill_zonr_mclat()

! Compute zonal height field from hydrostatic and geostrophic balances

zonz(37,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)
zonz(38,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)

! Use geostrophic balance eqn to get 1000 mb heights at all latitudes

dyo2g = 2.5 * dlat / grav2  ! spacing (m) between zonavg values divided by 2g

do ilatn = 39,74
   ilats = 75 - ilatn 

! Coriolis parameter at midpoint between zonavg values (north and south hemispheres)
   
   fcorn = 2. * omega * sin(2.5 * (ilatn-38) * pio180)
   fcors = 2. * omega * sin(2.5 * (ilats-37) * pio180)

! Compute and apply delta zp using midpoint average of two zonu values

   zonz(ilatn,1) = zonz(ilatn-1,1)   &
      - fcorn * dyo2g * (zonu(ilatn-1,1) + zonu(ilatn,1))
   zonz(ilats,1) = zonz(ilats+1,1)  &
      - fcorn * dyo2g * (zonu(ilats+1,1) + zonu(ilats,1))
enddo

! Use hydrostatic equation to get all heights above 1000 mb surface

cpo2g = cp / grav2

do ilat = 1,74
   do iplev = 2,22
      pilo = (1.e-5 * zonp_vect(iplev-1)) ** rocp
      pihi = (1.e-5 * zonp_vect(iplev)) ** rocp
      thetavlo = zont(ilat,iplev-1) * (1. + eps_virt * zonr(ilat,iplev-1))  &
               / pilo
      thetavhi = zont(ilat,iplev) * (1. + eps_virt * zonr(ilat,iplev))  &
               / pihi
      zonz(ilat,iplev) = zonz(ilat,iplev-1) + cpo2g * (pilo - pihi)  &
         * (thetavlo + thetavhi)

!write(io6,3302) ilat,iplev,zonp_vect(iplev),pihi,thetavhi,zont(ilat,iplev)  &
!    ,zonr(ilat,iplev),zonz(ilat,iplev)
!3302 format('zonz1 ',2i6,12e15.6)
         
   enddo
enddo

return
end subroutine zonavg_init

!===============================================================================

subroutine fill_zonr_mclat()

use misc_coms,   only: io6, imonth1, idate1, iyear1
use mem_mclat,   only: slat, mclat, ypp_mclat, mcol, mclat_spline
use mem_zonavg,  only: zonr_vect, zonp_vect, zonr
use mem_radiate, only: jday

implicit none   

integer :: ilat,lv
real :: alat
integer, external :: julday

! Subroutine fill_zonr_mclat fills the ZONR array with specific humidity values 
! by means of interpolation from the McLatchy sounding.  The ZONR array is 2-D, 
! defined at 22 vertical pressure levels and 74 latitudes (2.5 degree spacing),
! and is required to supplement the zonu and zont arrays which are filled with
! monthly climatological data from the U.K. Met office.  

! Get current Julian day

jday = julday(imonth1,idate1,iyear1)

! Compute spline coefficients in preparation for interpolation

call mclat_spline(jday)  

! Loop over latitude columns of zonavg data

do ilat = 1,74
   alat = -93.75 + 2.5 * ilat

! Loop over vertical levels in McLatchy sounding

   do lv = 1,33

! Spline-interpolate (by latitude) Mclatchy pressure, vapor density, and 
! total density to current latitude
 
      call spline2(13,slat,mclat(1,lv,2),ypp_mclat(1,lv,2),alat,mcol(lv,2))
      call spline2(13,slat,mclat(1,lv,4),ypp_mclat(1,lv,4),alat,mcol(lv,4))
      call spline2(13,slat,mclat(1,lv,6),ypp_mclat(1,lv,6),alat,mcol(lv,6))

! Compute specific humidity from quotient of vapor and total density

      mcol(lv,4) = mcol(lv,4) / mcol(lv,6)
   enddo

! Vertically interpolate Mclatchy vapor specific humidity BY PRESSURE to zonavg 
! data levels

   call pintrp_ee(33,mcol(1,4),mcol(1,2),22,zonr_vect,zonp_vect)
   zonr(ilat,1:22) = zonr_vect(1:22)
enddo

return
end subroutine fill_zonr_mclat



