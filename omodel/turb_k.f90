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
subroutine turb_k()

use mem_turb,   only: hkm, vkm, vkh
use mem_ijtabs, only: istp, jtab_w, itab_w, mrl_endl
use mem_grid,   only: dzim, lpw, mza
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer :: j,iw,k,iwp,mrl

real :: dzimo2(mza)   ! automatic array
real :: dzim2 (mza)   ! automatic array

do k = 2,mza-2
   dzimo2(k) = .5 * dzim(k)
   dzim2(k) = dzim(k) * dzim(k)
enddo

! Horizontal loop over W/T points

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw)
do j = 1,jtab_w(34)%jend(mrl); iw = jtab_w(34)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call turb_k0(iw,dzimo2,dzim2)
   
enddo
!$omp end parallel do
endif
call rsub('W',34)

! Grid 1 lateral boundary

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw,iwp,k)
do j = 1,jtab_w(35)%jend(mrl); iw = jtab_w(35)%iw(j)
   iwp = itab_w(iw)%iwp
!----------------------------------------------------------------------
call qsub('W',iw)
   do k = lpw(iw),mza
      hkm(k,iw) = hkm(k,iwp)
      vkm(k,iw) = vkm(k,iwp)
      vkh(k,iw) = vkh(k,iwp)
   enddo
enddo
!$omp end parallel do
endif
call rsub('W',35)

return
end subroutine turb_k

!===============================================================================

subroutine turb_k0(iw,dzimo2,dzim2)

use mem_turb,    only: vkm, vkh, hkm
use misc_coms,   only: io6, idiffk, csx, csz, zkhkm, akmin, meshtype
use mem_ijtabs,  only: itab_w
use mem_basic,   only: uc, vc, rho
use consts_coms, only: vonk
use mem_grid,    only: mza, lpw, arw0, dzm, zm, lpu, lpv,  &
                       unx, uny, unz, vnx, vny, vnz

implicit none

integer, intent(in) :: iw

real, intent(in) :: dzimo2(mza)
real, intent(in) :: dzim2 (mza)

integer :: k,mrlw,ka,npoly,j,iv
real :: enfl,richnum,ambda,ambda2,hill_term,richnum_term  &
   ,scalen_asympt,scalen_vert,scalen_horiz,vkz2,dudz,dvdz,dwdz,volfrac  &
   ,bkmin,zi,sbf,polyi

real :: bvfreq2(mza)
real :: strain2(mza)
real :: vx(mza),vy(mza),vz(mza)

!real, parameter :: richcof = 1.0     ! Lilly-Richardson number coefficient
!real, parameter :: rchmax = 1.0 + 9.0 * richcof
real, parameter :: rchmax = 3.  ! Test with new asympt vert scale length
real, parameter :: rmin = -100.
real, parameter :: rmax = 1. / 3.  !  1. / zkhkm(mrl)

mrlw = itab_w(iw)%mrlw

! Compute Brunt Vaisala frequency squared: bvfreq2(k)

call bruvais (iw,bvfreq2)
   
! Compute vertical strain rate squared (based on du/dz and dv/dz): strain2
! Transform uc into earth components for each iw since coefficients are
! available (Similar to leaf3_eachcall)

npoly = itab_w(iw)%npoly
polyi = 1. / real(npoly)

ka = lpw(iw)

vx(ka:mza-1) = 0.
vy(ka:mza-1) = 0.
vz(ka:mza-1) = 0.
   
! Loop over U/V neighbors of W

do j = 1,npoly

   if (meshtype == 1) then
      iv = itab_w(iw)%iu(j)
   else
      iv = itab_w(iw)%iv(j)
   endif

! Vertical loop over T levels

   do k = ka,mza-2

! Sum over V neighbors the 3 EARTH components of horizontal wind   
   
      vx(k) = vx(k) + uc(k,iv) * unx(iv) + vc(k,iv) * vnx(iv)
      vy(k) = vy(k) + uc(k,iv) * uny(iv) + vc(k,iv) * vny(iv)
      vz(k) = vz(k) + uc(k,iv) * unz(iv) + vc(k,iv) * vnz(iv)
   enddo

enddo

vx(ka:mza-1) = vx(ka:mza-1) * polyi
vy(ka:mza-1) = vy(ka:mza-1) * polyi
vz(ka:mza-1) = vz(ka:mza-1) * polyi

do k = ka,mza-2
   strain2(k) = dzim2(k) * ((vx(k+1) - vx(k))**2  &
                          + (vy(k+1) - vy(k))**2  &
                          + (vz(k+1) - vz(k))**2)

enddo
      
! Select turbulence scheme

if (idiffk(mrlw) == 2) then

! This is Smagorinsky-Lilly-Hill parameterization

   scalen_horiz = csx(mrlw) * sqrt(arw0(iw))  ! change this later?
   scalen_asympt = csx(mrlw) * 300.

   do k = ka,mza-2

!  Compute Richardson number and Lilly Richardson-number term

      richnum = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-15)))
      richnum_term = min(rchmax,sqrt(max(0.,(1.-zkhkm(mrlw)*richnum))))

! Compute vertical and net scale lengths: scalen_vert & ambda

      scalen_vert = csz(mrlw) * dzm(k)
      ambda = max(scalen_vert,min(scalen_asympt,scalen_horiz))
      ambda2 = ambda ** 2
            
      vkz2 = (vonk * (zm(k) - zm(ka-1))) ** 2  ! could move into vctr* outside iw loop 

      if (bvfreq2(k) < -1.e-12) then
         hill_term = sqrt(-bvfreq2(k))
      else
         hill_term = 0.
      endif

      vkm(k,iw) = .5 * (rho(k,iw) + rho(k+1,iw))   & ! density factor
                * vkz2 * ambda2 / (vkz2 + ambda2)  & ! lengthscale^2 factor
                * (sqrt(strain2(k)) + hill_term)   & ! strain rate + Hill term
                * richnum_term                       ! Lilly Richnum term

!!!!!!!!!!!!!!!!!!!!!!!!hs3
!     vkm(k,iw) = 0.   !hs3
!!!!!!!!!!!!!!!!!!!!!!!!hs3

      vkh(k,iw) = vkm(k,iw) * zkhkm(mrlw)  ! Need to change this factor later

   enddo
      
! Zero values for top and bottom boundaries   
   
   vkm(ka-1,iw) = 0.
   vkh(ka-1,iw) = 0.
   vkm(mza-1,iw) = 0.
   vkh(mza-1,iw) = 0.

! Horizontal diffusion coefficient (current version)

   bkmin = akmin(1) * .075 * arw0(iw) ** .666667
! akmin hardwired for "grid 1" value for now

   do k = ka,mza-1

!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      hkm(k,iw) = max(vkm(k-1,iw),vkm(k,iw),bkmin * real(rho(k,iw)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
!    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! elseif (idiffk(mrlw) == 1) then
   
! This is Mellor and Yamada parameterization

! elseif (idiffk(mrlw) == 4) then
   
! This is Deardorff parameterization
! Deardorff (needs 3d strain rate here and cross terms in veltend_long)
   
elseif (idiffk(mrlw) == 7) then

! This is boundary layer parameterization based on Taylor hypothesis
! and originally programmed by Haroldo Fraga de Campos Velho for RAMS

   call kcgcv2(iw,mrlw,bvfreq2,strain2)

! Horizontal diffusion coefficient (current version)

   bkmin = akmin(1) * .075 * arw0(iw) ** .666667
! akmin hardwired for "grid 1" value for now

   do k = ka,mza-1
      hkm(k,iw) = max(vkm(k-1,iw),vkm(k,iw),bkmin * real(rho(k,iw)))
   enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
!    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   
else
   
! For IDIFFK not equal to any of the above options, set diffusion
! coefficients to zero
   
   do k = ka-1,mza-1
      hkm(k,iw) = 0.
      vkm(k,iw) = 0.
      vkh(k,iw) = 0.
   enddo

endif

return
end subroutine turb_k0

!===============================================================================

subroutine bruvais(iw,bvfreq2)

use mem_basic,   only: theta, sh_v
use consts_coms, only: alvl, grav2, cp, rvap
use mem_grid,    only: mza, lpw, dzim
use micro_coms,  only: level
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw
real, intent(out) :: bvfreq2(mza)

real, parameter :: c2 = alvl * alvl / (cp * rvap)
                  
integer :: k
real :: exnrf,tempw,tempw2,thetaw,rvlsw

real, dimension(mza) :: thetav  ! automatic array

real, external :: rhovsl

! Calculate brunt-vaisalla frequency squared (en2) at W point

if (level == 0) then
   
   do k = lpw(iw),mza-2
      bvfreq2(k) = grav2 * dzim(k)  &
         * (theta(k+1,iw) - theta(k,iw)) / (theta(k+1,iw) + theta(k,iw)) 
   enddo
      
!!!elseif (level == 1) then BYPASS SATURATED CONDITION DUE TO APPARENT PROBLEMS
else
   
   do k = lpw(iw),mza-1
      thetav(k) = theta(k,iw) * (1. + .61 * sh_v(k,iw))
   enddo

   do k = lpw(iw),mza-2
      bvfreq2(k) = grav2 * dzim(k)  &
         * (thetav(k+1) - thetav(k)) / (thetav(k+1) + thetav(k)) 
   enddo

!else    ! level >= 2
   
!   do k = lpw(iw),mza-1
!      vctr1(k) = theta(k,iw) * (1. + .61 * sh_v(k,iw))  ! theta_v
!      vctr2(k) = sh_w(k,iw) - sh_v(k,iw)                ! condensate spec. hum.
!      exnrf = (press(k,iw) / p00) ** rocp               ! Exner function
!      vctr3(k) = theta(k,iw) * exnrf                    ! temperature (K)
!      vctr4(k) = rhovsl(vctr3(k)-273.15)                ! sat vapor density over liq

!   enddo

!   do k = lpw(iw),mza-2
      
!      if (sh_c(k,iw) > 0.) then

!         tempw = .5 * (vctr3(k) + vctr3(k+1))
!         tempw2 = tempw ** 2              ! multiplying RAMS num & denom by this
!         rvlsw = .5 * (vctr4(k) + vctr4(k+1))
!         thetaw = .5 * (theta(k,iw) + theta(k+1,iw))
      
!         vctr10(k) = grav * dzim(k) * (  &
!            (tempw2 + alvlor * rvlsw * tempw) / (tempw2 + c2 * rvlsw)  &
!            * ((theta(k+1,iw) - theta(k,iw)) / thetaw  &
!            +  alvlocp / tempw * (vctr4(k+1) - vctr4(k)))  &
!            - (sh_w(k+1,iw) - sh_w(k,iw)) )
            
!      else  
          
!         vctr10(k) = grav2 * dzim(k)  &
!            * (vctr1(k+1) - vctr1(k)) / (vctr1(k+1) + vctr1(k))
               
!   enddo
      
endif

return
end subroutine bruvais

!===============================================================================

!     subroutine friclyr(m1,m2,m3,scr1,ustarl,tstarl,ustarw,tstarw
!    +    ,pctlnd,theta,rtgt)
!     include 'rcommons.h'
!     dimension scr1(m1,m2,m3),ustarl(m2,m3),ustarw(m2,m3)
!    +   ,tstarl(m2,m3),tstarw(m2,m3),pctlnd(m2,m3),theta(m1,m2,m3)
!    +   ,rtgt(m2,m3)
!
!     do j=1,n3
!        do i=1,n2
!           ust=ustarl(i,j)*pctlnd(i,j)+ustarw(i,j)*(1.-pctlnd(i,j))
!           tst=tstarl(i,j)*pctlnd(i,j)+tstarw(i,j)*(1.-pctlnd(i,j))
!           if(abs(tst).lt.1.e-10)tst=1.e-10
!           obuklen=theta(2,i,j)*ust*ust/(3.92*tst)
!           do k=2,m1-1
!              zst=zt(k)*rtgt(i,j)
!              xi=zst/obuklen
!              if(xi.le.0.)then
!                 phim=(1.-15.*xi)**-0.25
!              else
!                 phim=1.+4.7*xi
!              endif
!              scr1(k,i,j)=min(scr1(k,i,j),0.4*ust*zst/phim)
!           enddo
!        enddo
!     enddo
!     return
!     end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!===============================================================================

subroutine kcgcv2(iw,mrlw,bvfreq2,strain2)

!################################################################
!#                                                              #
!# PURPOSE: Compute vertical turbulent exchange coefficient     #
!#          based on Taylor theory                              #
!#                                                              #
!# Original programmer: Haroldo Fraga de Campos Velho           #
!#                                                              #
!#            Permanent Address:                                #
!#                                                              #
!#                LAC-INPE                                      #
!#                P.O. Box 515                                  #
!#                12.201-970 - Sao Jose dos Campos (SP)         #
!#                Brasil                                        #
!#                                                              #
!#                Fax:   +55  012  345-6375                     #
!#                Phone: +55  012  345-6356                     #
!#                E-mail: haroldo@lac.inpe.br                   #
!#                Home-page: http://www.lac.inpe.br/~haroldo/   #
!#                                                              #
!#                                                              #
!# Date: April 15, 2000                                         #
!#       Last alteration: June 19, 2006                         #
!#                                                              #
!# Modified for OLAM November 4, 2006                           #
!#                                                              #
!# REFERENCES:                                                  #
!#                                                              #
!# 1. G.A. Degrazia,  H.F. Campos Velho,  J.C. Carvalho (1997): #
!#      "Nonlocal  Exchange  Coefficients  for the  Convective  #
!#      Boundary  Layer  Derived  from  Spectral  Properties",  #
!#      Beitrage  zur  Physik der Atmosphare, Vol. 70,  No. 1,  #
!#      pp. 57-64.                                              #
!#                                                              #
!# 2. G.A. Degrazia,  O.L.L. Moraes  (1992): "A Model for Eddy  #
!#      Diffusivity  in  a  Stable Boundary Layer", Boundary-   #
!#      Layer Meteorology, Vol. 58, pp. 205-124.                #
!#                                                              #
!# 3. G.A. Degrazia, D. Anfossi,  J.C. Carvalho,  T. Tirabassi, #
!#      H.F. Campos Velho (2000): "Turbulence Parameterization  #
!#      for PBL Dispersion Models in All Stability Conditions", #
!#      Atmospheric Environment (to appear)                     #
!#                                                              #
!################################################################
!#                                                              #
!# OUTPUT VARIABLES:                                            #
!#   VKM, VKH: Eddy diffusivity coefficients                    #
!#                                                              #
!################################################################

use mem_basic,   only: theta, sh_v, rho
use mem_turb,    only: vkm, vkh, hkm, sflux_t, sflux_r, ustar
use consts_coms, only: vonk, grav
use misc_coms,   only: io6, idiffk, zkhkm, akmin
use mem_grid,    only: mza, lpw, zm, zt

implicit none

! Input variables

integer, intent(in) :: iw        ! horizontal index of host model (OLAM)
integer, intent(in) :: mrlw      ! mesh refinement level for IW column  

real, intent(in) :: bvfreq2(mza) ! Brunt-Vaisala frequency squared
real, intent(in) :: strain2(mza) ! Vertical strain rate squared 

! Local variables

integer :: ka         ! vertical k index of lowest predicted W or T level
integer :: k          ! vertical k index of any model level
integer :: kzi        ! vertical k index of PBL top

real :: lamb          ! local Monin-Obukhov length, see reference 2.
real :: q             ! stability function, Eq. 2.17 in reference 1.
real :: corr          ! correction factor
real :: kn            ! vertical eddy diffusivity from wind shear
real :: sbf           ! surface buoyancy flux
real :: z             ! height of grid level above surface
real :: zl            ! Monin Obukhov length
real :: zi            ! PBL height above ground level
real :: zozi          ! z / zi
real :: wstar         ! convective velocity scale
real :: aux1          ! scratch variable
real :: aux2          ! scratch variable
real :: rnu           !
real :: den           !
real :: thvirt        ! virtual potential temperature 
real :: richnum       ! Richardson number
real :: richnum_term  ! Lilly Richardson-number coefficient

! 1 - Stable parameterization (see ref. 2)

real, parameter :: alp1 = 1.5
real, parameter :: alp2 = 1.
real, parameter :: auxl = 1.5 * alp1 - alp2
real, parameter :: fc = 1.e-4      ! Coriolis parameter (see ref. 3)
real, parameter :: rmin = -100.
real, parameter :: rmax = 1. / 3.  !  1. / zkhkm(mrl)

ka = lpw(iw)

! Compute surface virtual potential temperature
     
thvirt = theta(ka,iw) * (1. + .61 * sh_v(ka,iw))

! Compute surface buoyancy flux

sbf = grav / thvirt  &      
    * (sflux_t(iw) * (1. + .61 * sh_v(ka,iw))  &
    +  sflux_r(iw) * .61 * theta(ka,iw))

! Compute Monin Obukhov height

zl = 500. 
if (sbf /= 0.) zl = -ustar(iw) ** 3 / (vonk * sbf)

! Compute boundary layer depth

if (zl <= 500.0 .and. zl > 0) then
      
! Stable case (from Zilitinkevich 1972)      

   zi = 2.4e3 * ustar(iw) ** 1.5

else
      
! Unstable case - using gradient Richardon number...

! Initialize boundary layer height (minimum value is 1/2 deltaz)

   zi = zt(ka) - zm(ka-1)

! Vertical loop over W levels

   do k = ka,mza-2
      
!  Compute Richardson number and Lilly Richardson-number term

      richnum = max(rmin,min(rmax,bvfreq2(k) / max(strain2(k),1.e-15)))
      richnum_term = min(10.,sqrt(max(0.,(1.-zkhkm(mrlw)*richnum))))

! Find the lowest W(k) level where richnum_term < 1.e-5 (the lowest neutral or 
! stable layer).  Take zt(k) to be top of unstable layer.

      if (richnum_term <= 1.e-5) then
         zi = zt(k) - zm(ka-1)
         exit
      endif

   enddo
      
endif

! Compute convective velocity scale

wstar = (max(0., sbf * zi)) ** .3333333
   
! Find the highest model level zt(k) that is at or below boundary layer
! height h.

kzi = ka
do while (zm(kzi) < zm(ka-1) + zi)
   kzi = kzi + 1
enddo

! For model W levels at and above zm(kzi), set vertical turbulent mixing 
! coefficient to zero. (Lowest level of vkm = 0 is at height zm(kzi)-zm(ka-1), 
! which is 1/2 deltaz above zi at height zt(kzi).)

! Vertical loop over W levels

do k = kzi,mza-1
   vkm(k,iw) = 0.
   vkh(k,iw) = 0.
enddo
      
! For model levels below k = kzi, set vertical turbulent mixing 
! coefficient according to Taylor theory.

! Vertical loop over W levels

do k = ka,kzi-1

! Height above ground level of current model level

   z = zt(k) - zm(ka-1)

! Nondimensional vertical coordinate

   zozi = min(1.,z/zi)

! Contribution from wind shear to turbulence (see Eq(28) in ref. 3)

   aux1 = .4 * (1 - zozi) ** .85 * ustar(iw) * z
   aux2 = (1. + 15. * (fc * z / ustar(iw))) ** 1.333333
   kn = aux1 / aux2

   if (abs(zl) > 500.) then

! Nearly-neutral case

      vkm(k,iw) = kn

   elseif (zl > 0.) then

! Stable case

! 1.2 - Calculate the local Monin-Obukhov length (see Eq(3) in ref. 2)

      lamb = abs(zl) * (1. - zozi) ** auxl

! 1.3 - Stable eddy diffusivity (see Eq(29) in ref. 3)

      rnu = .4 * (1. - zozi) ** .75
      den = 1. + 3.7 * z / lamb
      vkm(k,iw) = rnu / den * ustar(iw) * z

   else

! Convective case

! Correction factor (see Eq. (14) in ref. 3)

      corr = sqrt(.01 * zi / (-zl)) ! See Eq. 14 4

! Stability function: Eq. (2.17) in reference 1.

      q = 1. - exp(-4. * zozi) - .0003 * exp(8. * zozi)

! Eddy diffusivity: Eq. (26) in reference 3.        
! Eddy diffusivity exponent = 1.333333: Eq. (2.18) in reference 1.

      vkm(k,iw) = .16 * wstar * zi * corr * q ** 1.333333

! Shear effect contribution (see Eq. (15), ref. 3)

      vkm(k,iw) = vkm(k,iw) + kn
      vkh(k,iw) = vkm(k,iw) * zkhkm(mrlw)  ! Need to change this factor later

   endif

enddo

! Zero values for bottom boundary 
   
vkm(ka-1,iw) = 0.
vkh(ka-1,iw) = 0.

return
end subroutine kcgcv2


