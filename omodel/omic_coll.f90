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
subroutine effxy(lpw0,k1,k2,rx,qr,emb,tx,eff)

use micro_coms, only: mza0, ncat, jnmb, rxmin, neff, cfmasi, pwmasi
use misc_coms,  only: io6

implicit none

integer, intent(in) :: lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in) :: rx (mza0,ncat)
real, intent(in) :: qr (mza0,ncat)
real, intent(in) :: emb(mza0,ncat)
real, intent(in) :: tx (mza0,ncat)

real, intent(out) :: eff(mza0,neff)

real, parameter :: al10 = log(10.)
real, parameter :: ef0 = -0.7 * al10
real, parameter :: ef1 = .035 * al10

integer         :: k
real            :: dmb

! This subroutine sets COALLESCENCE EFFICIENCIES for all hydrometeor collisions.
! Some of these depend on hydrometeor temperatures, while others are constant.

! COLLISION efficiencies are not considered here, but are all set when the
! collection tables are generated.

! 1 = cc,cd,cr,   cs,ca,cg,ch,
!        dd,dr,   ds,da,dg,dh,
!              rp,rs,ra,rg,rh

eff(lpw0:mza0,1) = 1.0

! 2 = pp,ps,pa

if (jnmb(5) >= 1) then

   do k = k1(3),k2(3)

      if (abs(tx(k,3) + 14.) <= 2.) then
         eff(k,2) = 1.4
      else
         eff(k,2) = min(.2, exp( ef0 + ef1 * tx(k,3) ))
      endif

   enddo

! 3 = ss,sa

   do k = k1(4),k2(4)

      if (abs(tx(k,4) + 14.) <= 2.) then
         eff(k,3) = 1.4
      else
         eff(k,3) = min(.2, exp( ef0 + ef1 * tx(k,4) ))
      endif

   enddo

! 4 = aa

   do k = k1(5),k2(5)

      if (abs(tx(k,5) + 14.) <= 2.) then
         eff(k,4) = 1.4
      elseif (tx(k,5) >= -1.) then
         eff(k,4) = 1.
      else
         eff(k,4) = min(.2, exp( ef0 + ef1 * tx(k,5) ))
      endif

   enddo

endif

! 5 = pg,sg,ag,gg,gh

if (jnmb(6) >= 1) then

   do k = k1(6),k2(6)

      if (qr(k,6) > 0.) then
         eff(k,5) = 1.0
      else
         eff(k,5) = min(.2, exp( ef0 + ef1 * tx(k,6) ))
      endif

   enddo

endif

! 6 = ph,sh,ah

if (jnmb(7) >= 1) then

   do k = k1(7),k2(7)

      if (qr(k,7) > 0.) then
         eff(k,6) = 1.0
      else
         eff(k,6) = min(.2, exp( ef0 + ef1 * tx(k,7) ))
      endif

! 7 = hh (experimental; should be tested and improved; large hail probably
!         should not coallesce)

      eff(k,7) = max(0., .1 + .005 * tx(k,7))

   enddo

endif

! 8 = rr (rain-rain collisional breakup at large sizes)

if (jnmb(2) >= 1) then

   do k = k1(2),k2(2)

      ! RAMS version

    !  if (emb(k,2) < .113e-6) then
    !     eff(k,8) = 1.0
    !  elseif (emb(k,2) > .158e-5) then
    !     eff(k,8) = -5.0
    !  else
    !     eff(k,8) = 2. - exp(.1326e7 * (emb(k,2) - .113e-6))
    !  endif

    ! Version designed to limit max droplet size to about 4 mm

      if (emb(k,2) < .524e-6) then      ! This is 1 mm diameter
         eff(k,8) = 1.0
      elseif (emb(k,2) > 33.536e-6) then  ! This is 4 mm diameter
         eff(k,8) = -1.0
      else
         dmb = (cfmasi(2) * emb(k,2)) ** pwmasi(2)
         eff(k,8) = 1.0 - 666.667 * (dmb - 1.e-3)
      endif

   enddo

endif

end subroutine effxy

!===============================================================================

subroutine cols(mx,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,exxxx)

use micro_coms, only: mza0, ncat, rxmin, ipair, coltabc, neff
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: exxxx(mza0)

integer :: ipc,k

real :: tabc,colc

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle

   ipc = ipair(jhcat(k,mx),jhcat(k,mx),1)

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = colfac(k) * eff(k,meff) * cx(k,mx)**2 * exp(tabc)

   exxxx(k) = min(0.5 * cx(k,mx),colc)

enddo

end subroutine cols

!===============================================================================

subroutine col1188(mx,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxxxz,exxxx,exxxz)

use micro_coms, only: mza0, ncat, rxmin, ipair, jnmb, neff, &
                      coltabc, coltabx
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxxxz(mza0,2)
real, intent(inout) :: exxxx(mza0)
real, intent(inout) :: exxxz(mza0)

integer :: k,ipc,ipc2,ipx,indx

real :: c1,c2,tabc,tabc2,tabx,embxzi,colc

if (mx == 1) then         ! Cloud-cloud (mx = 1 and mz = 8)
   embxzi = 1. / 15.e-12  ! Droplet mass for transfer to drizzle
else                      ! Drizzle-drizzle (mx = 8 and mz = 2)
   embxzi = 1. / 4.e-9    ! Droplet mass for transfer to rain
endif

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle

   ipc  = ipair(jhcat(k,mx),jhcat(k,mx),1)
   ipx  = ipair(jhcat(k,mx),jhcat(k,mx),4)

   c1 = eff(k,meff) * colfac(k) * cx(k,mx) ** 2
   c2 = 2. * c1

! Interpolate from coltabx

   tabx =      wct1(k,mx) * wct1(k,mx) * coltabx(ict1(k,mx),ict1(k,mx),ipx) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabx(ict1(k,mx),ict2(k,mx),ipx) &
        +      wct2(k,mx) * wct2(k,mx) * coltabx(ict2(k,mx),ict2(k,mx),ipx)

! Hydrometeor mass transfer

   rxxxz(k,1) = min(rx(k,mx), c2 * exp(tabx)) ! c2 since both tabx droplets go to z
   rxxxz(k,2) = rxxxz(k,1) * qx(k,mx)

   if (jnmb(mz) < 5) cycle

! Hydrometeor number transfer

   exxxz(k) = rxxxz(k,1) * embxzi

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = min(0.5 * cx(k,mx), c1 * exp(tabc)) ! c1 so colc represents total number loss

! Hydrometeor number loss

   exxxx(k) = colc

enddo

end subroutine col1188

!===============================================================================

subroutine col1882(mx,my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxyxy,rxyxz,rxyyz,exyxx,exyyz)

use micro_coms, only: mza0, ncat, rxmin, ipair, jnmb, neff, &
                      coltabc, coltabx, coltaby, driz_gammq, emb1
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxyxy(mza0,2)
real, intent(inout) :: rxyxz(mza0,2)
real, intent(inout) :: rxyyz(mza0,2)
real, intent(inout) :: exyxx(mza0)
real, intent(inout) :: exyyz(mza0)

integer :: k,ipc,ipx,ipx2,ipy2

real :: c1,tabc,tabc2,tabx,tabx2,taby2,gammq

real, parameter :: embxzi = 1. / 4.e-9 ! Droplet mass for transfer to rain

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

! FOR CLOUD-DRIZZLE COLLISIONS THAT MAKE LARGER DRIZZLE

   ipc  = ipair(jhcat(k,mx),jhcat(k,my),1)
   ipx  = ipair(jhcat(k,mx),jhcat(k,my),4)

! FOR CLOUD-DRIZZLE COLLISIONS THAT MOVE DRIZZLE TO RAIN

   ipx2 = ipair(jhcat(k,mx),jhcat(k,my),2)
   ipy2 = ipair(jhcat(k,mx),jhcat(k,my),3)

   c1 = eff(k,meff) * colfac(k) * cx(k,mx) * cx(k,my)

! FOR CLOUD-DRIZZLE COLLISIONS THAT MAKE LARGER DRIZZLE

! Interpolate from coltabx

   tabx = wct1(k,mx) * wct1(k,my) * coltabx (ict1(k,mx),ict1(k,my),ipx) &
        + wct2(k,mx) * wct1(k,my) * coltabx (ict2(k,mx),ict1(k,my),ipx) &
        + wct1(k,mx) * wct2(k,my) * coltabx (ict1(k,mx),ict2(k,my),ipx) &
        + wct2(k,mx) * wct2(k,my) * coltabx (ict2(k,mx),ict2(k,my),ipx)

! Hydrometeor mass transfer

   rxyxy(k,1) = min(rx(k,mx), c1 * exp(tabx))
   rxyxy(k,2) = rxyxy(k,1) * qx(k,mx)

! Interpolate from coltabc

   tabc = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc) &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc) &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc) &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

! Hydrometeor number loss (This includes situation when mz=rain)

   exyxx(k) = min(0.5 * cx(k,mx), c1 * exp(tabc))

! FOR CLOUD-DRIZZLE COLLISIONS THAT MOVE DRIZZLE TO RAIN

! Interpolate from coltabx

   tabx2 = wct1(k,mx) * wct1(k,my) * coltabx (ict1(k,mx),ict1(k,my),ipx2) &
         + wct2(k,mx) * wct1(k,my) * coltabx (ict2(k,mx),ict1(k,my),ipx2) &
         + wct1(k,mx) * wct2(k,my) * coltabx (ict1(k,mx),ict2(k,my),ipx2) &
         + wct2(k,mx) * wct2(k,my) * coltabx (ict2(k,mx),ict2(k,my),ipx2)

! Hydrometeor mass transfer

   rxyxz(k,1) = min(rx(k,mx), c1 * exp(tabx2))
   rxyxz(k,2) = rxyxz(k,1) * qx(k,mx)

! Interpolate from coltaby

   taby2 = wct1(k,my) * wct1(k,mx) * coltaby (ict1(k,my),ict1(k,mx),ipy2) &
         + wct2(k,my) * wct1(k,mx) * coltaby (ict2(k,my),ict1(k,mx),ipy2) &
         + wct1(k,my) * wct2(k,mx) * coltaby (ict1(k,my),ict2(k,mx),ipy2) &
         + wct2(k,my) * wct2(k,mx) * coltaby (ict2(k,my),ict2(k,mx),ipy2)

! Hydrometeor mass transfer

   gammq = wct1(k,my) * driz_gammq(ict1(k,my)) &
         + wct2(k,my) * driz_gammq(ict2(k,my))

   rxyyz(k,1) = min(rx(k,my)*gammq, c1 * exp(taby2))
   rxyyz(k,2) = rxyyz(k,1) * qx(k,my)

! Hydrometeor number transfer

   exyyz(k) = rxyyz(k,1) * embxzi

enddo

end subroutine col1882

!===============================================================================

subroutine col3344(mx,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxxxz,exxxx,exxxz)

use micro_coms, only: mza0, ncat, rxmin, ipair, jnmb, neff,  &
                      coltabc, coltabx
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac2(mza0)

real, intent(inout) :: rxxxz(mza0,2)
real, intent(inout) :: exxxx(mza0)
real, intent(inout) :: exxxz(mza0)

integer :: k,ipc,ipx

real :: c1,tabc,tabx,colc,colx

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle

   ipc = ipair(jhcat(k,mx),jhcat(k,mx),1)
   ipx = ipair(jhcat(k,mx),jhcat(k,mx),4)

   c1 = eff(k,meff) * colfac2(k) * cx(k,mx) ** 2

! Interpolate from coltabx

   tabx =      wct1(k,mx) * wct1(k,mx) * coltabx(ict1(k,mx),ict1(k,mx),ipx) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabx(ict1(k,mx),ict2(k,mx),ipx) &
        +      wct2(k,mx) * wct2(k,mx) * coltabx(ict2(k,mx),ict2(k,mx),ipx)

   colx = min(rx(k,mx), c1 * exp(tabx))

   rxxxz(k,1) = colx
   rxxxz(k,2) = colx * qx(k,mx)

   if (jnmb(mz) < 5) cycle

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = min(0.5 * cx(k,mx), c1 * exp(tabc))

   exxxx(k) = colc
   exxxz(k) = colc

enddo

end subroutine col3344

!===============================================================================

subroutine col3443(meff,j1,j2,iw0, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
   rxyxz,rxyyz,exyxx,exyyy,exyxz,exyyz)

use micro_coms, only: mza0, ncat, rxmin, neff, ipair, &
                      coltabc, coltabx, coltaby
use misc_coms,  only: io6

implicit none

integer, intent(in) :: meff

integer, intent(in) :: j1
integer, intent(in) :: j2
integer, intent(in) :: iw0

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxyxz(mza0,2)
real, intent(inout) :: rxyyz(mza0,2)
real, intent(inout) :: exyxx(mza0)
real, intent(inout) :: exyyy(mza0)
real, intent(inout) :: exyxz(mza0)
real, intent(inout) :: exyyz(mza0)

integer :: k,jhcatx,jhcaty,ipx,ipy,ipc

real :: c1,tabc,tabx,taby,colc,colx,coly

do k = j1,j2

   if (rx(k,3) < rxmin(3) .or. rx(k,4) < rxmin(4)) cycle

   jhcatx = jhcat(k,3)
   jhcaty = jhcat(k,4)

   ipc = ipair(jhcatx,jhcaty,1)
   ipx = ipair(jhcatx,jhcaty,4)
   ipy = ipair(jhcatx,jhcaty,5)

   c1 = eff(k,meff) * colfac(k) * cx(k,3) * cx(k,4)

! Interpolate from coltabx

   tabx = wct1(k,3) * wct1(k,4) * coltabx (ict1(k,3),ict1(k,4),ipx) &
        + wct2(k,3) * wct1(k,4) * coltabx (ict2(k,3),ict1(k,4),ipx) &
        + wct1(k,3) * wct2(k,4) * coltabx (ict1(k,3),ict2(k,4),ipx) &
        + wct2(k,3) * wct2(k,4) * coltabx (ict2(k,3),ict2(k,4),ipx)

   colx = min(rx(k,3), c1 * exp(tabx))

! Interpolate from coltaby

   taby = wct1(k,4) * wct1(k,3) * coltaby (ict1(k,4),ict1(k,3),ipy) &
        + wct2(k,4) * wct1(k,3) * coltaby (ict2(k,4),ict1(k,3),ipy) &
        + wct1(k,4) * wct2(k,3) * coltaby (ict1(k,4),ict2(k,3),ipy) &
        + wct2(k,4) * wct2(k,3) * coltaby (ict2(k,4),ict2(k,3),ipy)

   coly = min(rx(k,4), c1 * exp(taby))

   rxyxz(k,1) = colx
   rxyxz(k,2) = colx * qx(k,3)

   rxyyz(k,1) = coly
   rxyyz(k,2) = coly * qx(k,4)

! Interpolate from coltabc

   tabc = wct1(k,3) * wct1(k,4) * coltabc (ict1(k,3),ict1(k,4),ipc) &
        + wct2(k,3) * wct1(k,4) * coltabc (ict2(k,3),ict1(k,4),ipc) &
        + wct1(k,3) * wct2(k,4) * coltabc (ict1(k,3),ict2(k,4),ipc) &
        + wct2(k,3) * wct2(k,4) * coltabc (ict2(k,3),ict2(k,4),ipc)

   colc = c1 * exp(tabc)

   if (cx(k,3) > cx(k,4)) then
      exyyz(k) = min(cx(k,4),colc)
      exyxx(k) = min(cx(k,3),colc)
   else
      exyxz(k) = min(cx(k,3),colc)
      exyyy(k) = min(cx(k,4),colc)
   endif

! also loss for aerosol

enddo

end subroutine col3443

!===============================================================================

subroutine col1(mx,my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxyxz,exyxx)

use micro_coms, only: mza0, ncat, neff, rxmin, ipair, jnmb, &
                      coltabc, coltabx
use misc_coms,  only: io6

implicit none

integer, intent(in) :: mx
integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)
real, intent(inout) :: rxyxz(mza0,2)
real, intent(inout) :: exyxx(mza0)

integer :: k,ipx,ipc

real :: c1,tabc,tabx,colc,colx

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

   ipc = ipair(jhcat(k,mx),jhcat(k,my),1)
   ipx = ipair(jhcat(k,mx),jhcat(k,my),4)

   c1 = eff(k,meff) * colfac(k) * cx(k,mx) * cx(k,my)

! Interpolate from coltabx

   tabx = wct1(k,mx) * wct1(k,my) * coltabx (ict1(k,mx),ict1(k,my),ipx) &
        + wct2(k,mx) * wct1(k,my) * coltabx (ict2(k,mx),ict1(k,my),ipx) &
        + wct1(k,mx) * wct2(k,my) * coltabx (ict1(k,mx),ict2(k,my),ipx) &
        + wct2(k,mx) * wct2(k,my) * coltabx (ict2(k,mx),ict2(k,my),ipx)

   colx = min(rx(k,mx), c1 * exp(tabx))

   rxyxz(k,1) = colx
   rxyxz(k,2) = colx * qx(k,mx)

   if (jnmb(mx) < 5) cycle

! Interpolate from coltabc

   tabc = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc) &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc) &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc) &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

   colc = c1 * exp(tabc)

   exyxx(k) = min(colc,cx(k,mx))

! also loss for aerosol

enddo

end subroutine col1

!===============================================================================

subroutine col2(mx,my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
   rxyx3,rxyxz,rxyxy,rxyyz,exyx3,exyxx,exyyz)

use micro_coms, only: mza0, ncat, rxmin, ipair, coltabx, coltaby, coltabc, &
                      sipfac, pwmasi, emb1, gamsip13, gamsip24, emb0, neff, &
                      jnmb, ngam
use misc_coms,  only: io6
use therm_lib,  only: qtc

implicit none

integer, intent(in) :: mx
integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: emb (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxyx3(mza0,2)
real, intent(inout) :: rxyxz(mza0,2)
real, intent(inout) :: rxyxy(mza0,2)
real, intent(inout) :: rxyyz(mza0,2)
real, intent(inout) :: exyx3(mza0)
real, intent(inout) :: exyxx(mza0)
real, intent(inout) :: exyyz(mza0)

real, intent(in) :: dtl0

integer :: k,jhcatx,jhcaty,ipx,ipy,ipc,it,lx

real :: c1,c2,tabc,tabx,taby,colc,colx,coly,colc0,colqrx,colqry, &
        rcoal,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice, &
        area,cn13,cn24,sip,rsip,qrsip,rfinlz,xtoz

real :: alpha(16)
real :: beta(16)

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00.,00.,1.0,1.0,1.0,1.0,00.,00.,00.,00.,00.,1.0,1.0,1.0,1.0/
data beta  /00.,00.,00.,0.5,0.5,0.0,0.0,00.,00.,00.,00.,00.,0.5,0.5,0.5,0.5/

lx = 1
if (mx == 8) lx = 2

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)

   ipc = ipair(jhcatx,jhcaty,1)
   ipx = ipair(jhcatx,jhcaty,4)
   ipy = ipair(jhcatx,jhcaty,5)

   c2 = colfac(k) * cx(k,mx) * cx(k,my)
   c1 = eff(k,meff) * c2

! Interpolate from coltabx

   tabx = wct1(k,mx) * wct1(k,my) * coltabx (ict1(k,mx),ict1(k,my),ipx) &
        + wct2(k,mx) * wct1(k,my) * coltabx (ict2(k,mx),ict1(k,my),ipx) &
        + wct1(k,mx) * wct2(k,my) * coltabx (ict1(k,mx),ict2(k,my),ipx) &
        + wct2(k,mx) * wct2(k,my) * coltabx (ict2(k,mx),ict2(k,my),ipx)

   colx = min(rx(k,mx), c1 * exp(tabx))

! Interpolate from coltaby

   taby = wct1(k,my) * wct1(k,mx) * coltaby (ict1(k,my),ict1(k,mx),ipy) &
        + wct2(k,my) * wct1(k,mx) * coltaby (ict2(k,my),ict1(k,mx),ipy) &
        + wct1(k,my) * wct2(k,mx) * coltaby (ict1(k,my),ict2(k,mx),ipy) &
        + wct2(k,my) * wct2(k,mx) * coltaby (ict2(k,my),ict2(k,mx),ipy)

   coly = min(rx(k,my), c1 * exp(taby))

! Interpolate from coltabc

   tabc = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc) &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc) &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc) &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

   colc0 = c2 * exp(tabc)
   colc = colc0 * eff(k,meff)

   colqrx = colx * qx(k,mx)
   colqry = coly * qx(k,my)

   rcoal = colx + coly
   qrcoal = colqrx + colqry
   qcoal = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

! Secondary ice production: cn24 is the number fraction of collected cloud
! droplets larger than 24 microns and is obtained from an incomplete gamma
! function table. cn13 is the fraction of collected cloud droplets
! smaller than 13 microns. "area" is cross section area of collecting ice
! per m^3 of atmospheric volume.

! Saleeby(6/3/02): Hallett-Mossop is done for both cloud droplet modes, though
! contribution from the large droplet mode is minimal compared to the small
! droplet mode. Ice splintering is only done if number concentration is
! prognostic for at least one of the two hydrometeor species involved.

   if (tcoal > -8.0 .and. tcoal < -3.0 .and. &
      (jnmb(mx) == 5 .or. jnmb(my) == 5)) then  ! Steve, why this condition?

      area = cx(k,my) * sipfac(jhcaty) * emb(k,my) ** (2.*pwmasi(jhcaty)) ! remd rhoa
      it = max(1,nint(emb(k,mx) / emb1(mx) * ngam))

      cn13 = colc * gamsip13(lx,it) / (area * dtl0)
      cn24 = min(cx(k,mx),colc0) * gamsip24(lx,it)

      sip = 9.1e-10 * cn24 * cn13 ** .93           ! has units of #/m^3
      if (tcoal < -5.) then
         sip = 0.33333 * (tcoal + 8.) * sip
      else
         sip = -0.5 * (tcoal + 3.) * sip
      endif

      rsip = sip * emb0(3)           ! has units of kg/m^3
      qrsip = qcoal * rsip           ! has units of J/m^3

      rcoal = rcoal - rsip
      qrcoal = qrcoal - qrsip

      exyx3(k) = sip
      rxyx3(k,1) = rsip
      rxyx3(k,2) = qrsip

   endif

! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
! fracliq may be a little larger than fracx due to collected
! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
! or 1.2, or need ALPHA itself to be at least 1.0.

   rfinlz = min(rcoal,alpha(jhcaty) * coalliq + beta(jhcaty) * colx)

   xtoz = min(colx,rfinlz)

   if (my /= mz) then
      rxyxz(k,1) = xtoz
      rxyxy(k,1) = colx - xtoz
      rxyyz(k,1) = rfinlz - xtoz

      rxyxz(k,2) = qx(k,mx) * xtoz
      rxyxy(k,2) = qx(k,mx) * (colx - xtoz)
      rxyyz(k,2) = qx(k,my) * (rfinlz - xtoz)

      exyyz(k) = (rfinlz - xtoz) * min(colc,cx(k,my)) / max(1.e-20,coly)
   else
      rxyxz(k,1) = colx
      rxyxz(k,2) = qx(k,mx) * colx
   endif

   exyxx(k) = min(colc,cx(k,mx))

! also include loss of aerosol

enddo

end subroutine col2

!===============================================================================

subroutine col3(my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
   r2yy2,r2y2z,r2y2y,r2yyz,e2y22,e2y2z,e2yyy,e2yyz)

use micro_coms, only: mza0, ncat, rxmin, ipair, jnmb, neff, &
                      coltabx, coltaby, coltabc
use misc_coms,  only: io6
use therm_lib,  only: qtc

implicit none

integer, intent(in) :: my
integer, intent(in) :: mz
integer, intent(in) :: meff
integer, intent(in) :: j1
integer, intent(in) :: j2

integer, intent(in) :: jhcat(mza0,ncat)
integer, intent(in) :: ict1 (mza0,ncat)
integer, intent(in) :: ict2 (mza0,ncat)

real, intent(in) :: wct1(mza0,ncat)
real, intent(in) :: wct2(mza0,ncat)
real, intent(in) :: rx  (mza0,ncat)
real, intent(in) :: cx  (mza0,ncat)
real, intent(in) :: qx  (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: r2yy2(mza0,2)
real, intent(inout) :: r2y2z(mza0,2)
real, intent(inout) :: r2y2y(mza0,2)
real, intent(inout) :: r2yyz(mza0,2)
real, intent(inout) :: e2y22(mza0)
real, intent(inout) :: e2y2z(mza0)
real, intent(inout) :: e2yyy(mza0)
real, intent(inout) :: e2yyz(mza0)

integer :: k,ipx,ipy,ipc,jhcaty

real :: c1,tabc,tabx,taby,colc,colx,coly,colcx,colcy,colqrx,colqry, &
        coalnum,rcoal,qrcoal,qcoal,fracliq,coalliq,coalice, &
        xtoz,rfinlz,tcoal,cfinlz

real :: alpha(16)
real :: beta(16)

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00., 1., 1., 1., 1., 1.,00., 1., 1., 1., 1., 1., 1., 1., 1./
data beta  /00.,00., 2., 2., 2., 1., 0.,00., 2., 2., 2., 2., 2., 2., 2., 2./

do k = j1,j2

   if (rx(k,2) < rxmin(2) .or. rx(k,my) < rxmin(my)) cycle

   jhcaty = jhcat(k,my)

   ipc = ipair(jhcat(k,2),jhcaty,1)
   ipx = ipair(jhcat(k,2),jhcaty,4)
   ipy = ipair(jhcat(k,2),jhcaty,5)

   c1 = eff(k,meff) * colfac(k) * cx(k,2) * cx(k,my)

! Interpolate from coltabx

   tabx = wct1(k,2) * wct1(k,my) * coltabx(ict1(k,2),ict1(k,my),ipx) &
        + wct2(k,2) * wct1(k,my) * coltabx(ict2(k,2),ict1(k,my),ipx) &
        + wct1(k,2) * wct2(k,my) * coltabx(ict1(k,2),ict2(k,my),ipx) &
        + wct2(k,2) * wct2(k,my) * coltabx(ict2(k,2),ict2(k,my),ipx)

   colx = min(rx(k,2), c1 * exp(tabx))

! Interpolate from coltaby

   taby = wct1(k,my) * wct1(k,2) * coltaby(ict1(k,my),ict1(k,2),ipy) &
        + wct2(k,my) * wct1(k,2) * coltaby(ict2(k,my),ict1(k,2),ipy) &
        + wct1(k,my) * wct2(k,2) * coltaby(ict1(k,my),ict2(k,2),ipy) &
        + wct2(k,my) * wct2(k,2) * coltaby(ict2(k,my),ict2(k,2),ipy)

   coly = min(rx(k,my), c1 * exp(taby))

   if (jnmb(2) == 5) then

! Interpolate from coltabc

      tabc = wct1(k,2) * wct1(k,my) * coltabc(ict1(k,2),ict1(k,my),ipc) &
           + wct2(k,2) * wct1(k,my) * coltabc(ict2(k,2),ict1(k,my),ipc) &
           + wct1(k,2) * wct2(k,my) * coltabc(ict1(k,2),ict2(k,my),ipc) &
           + wct2(k,2) * wct2(k,my) * coltabc(ict2(k,2),ict2(k,my),ipc)

      colc = c1 * exp(tabc)

      colcx = min(cx(k,2),colc)
      colcy = min(cx(k,my),colc)

      coalnum = min(colcx,colcy)

   endif

   colqrx = colx * qx(k,2)
   colqry = coly * qx(k,my)

   rcoal  = colx + coly
   qrcoal = colqrx + colqry
   qcoal  = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

   if (fracliq >= .99) then

      r2yy2(k,1) = coly
      r2yy2(k,2) = colqry

      if (jnmb(2) == 5) e2yyy(k) = colcy

   else

      rfinlz = min(rcoal, alpha(jhcaty) * coalliq + beta(jhcaty) * colx)

      xtoz = min(colx,rfinlz)

      if (my /= mz) then
         r2y2z(k,1) = xtoz
         r2y2y(k,1) = colx - xtoz
         r2yyz(k,1) = rfinlz - xtoz

! NEED TO USE QCOAL TO TRANSFER Q?

         r2y2z(k,2) = qx(k,2) * xtoz
         r2y2y(k,2) = qx(k,2) * (colx - xtoz)
         r2yyz(k,2) = qx(k,my) * (rfinlz - xtoz)
      else
         r2y2z(k,1) = colx
         r2y2z(k,2) = qx(k,2) * colx
      endif

      if (jnmb(2) == 5) then

         if (my == mz) then

            e2y22(k) = colcx

         elseif (colcy >= colcx) then

            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)

            e2y2z(k) = cfinlz
            e2y22(k) = colcx - cfinlz
            e2yyy(k) = colcy

         else

            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)

            e2yyz(k) = cfinlz
            e2y22(k) = colcx
            e2yyy(k) = colcy - cfinlz
         endif

      endif

   endif

enddo

! also include loss of aerosol

end subroutine col3

!===============================================================================

subroutine colxfers(iw0,k1,k2,rx,cx,qr, &
   r1118,r8882,r1112,r1818,r1212,r8282,r3335,r4445,r3435,r3445, &
   r3535,r3636,r3737,r4545,r4646,r4747,r5656,r5757,r6767,r1413, &
   r1416,r1414,r1446,r1513,r1516,r1515,r1556,r1613,r1616,r8483, &
   r8486,r8484,r8446,r8583,r8586,r8585,r8556,r8683,r8686,r1713, &
   r1717,r8783,r8787,r2332,r2327,r2323,r2337,r2442,r2427,r2424, &
   r2447,r2552,r2527,r2525,r2557,r2662,r2627,r2626,r2667,r2772, &
   r2727,r0000,r1812,r1882, &
   e1111,e1118,e8888,e8882,e1112,e1811,e1211,e8288,e2222,e5555, &
   e6666,e7777,e3333,e3335,e4444,e4445,e3433,e3444,e3435,e3445, &
   e3533,e3633,e3733,e4544,e4644,e4744,e5655,e5755,e6766,e1413, &
   e1411,e1446,e1513,e1511,e1556,e1613,e1611,e8483,e8488,e8446, &
   e8583,e8588,e8556,e8683,e8688,e1713,e1711,e8783,e8788,e2322, &
   e2327,e2333,e2337,e2422,e2427,e2444,e2447,e2522,e2527,e2555, &
   e2557,e2622,e2627,e2666,e2667,e2722,e2727,e2777,e0000,e1882, &
   con_ccnx)

use micro_coms,  only: mza0, ncat, jnmb, iccn, rxmin
use misc_coms,   only: io6
use ccnbin_coms, only: nccntyp, nbins, relcon_bin, ihyg, iccntyp

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(inout) :: rx(mza0,ncat)
real, intent(inout) :: cx(mza0,ncat)
real, intent(inout) :: qr(mza0,ncat)

real, intent(inout) :: &
  r1118(mza0,2),r8882(mza0,2),r1112(mza0,2),r1818(mza0,2),r1212(mza0,2), &
  r8282(mza0,2),r3335(mza0,2),r4445(mza0,2),r3435(mza0,2),r3445(mza0,2), &
  r3535(mza0,2),r3636(mza0,2),r3737(mza0,2),r4545(mza0,2),r4646(mza0,2), &
  r4747(mza0,2),r5656(mza0,2),r5757(mza0,2),r6767(mza0,2),r1413(mza0,2), &
  r1416(mza0,2),r1414(mza0,2),r1446(mza0,2),r1513(mza0,2),r1516(mza0,2), &
  r1515(mza0,2),r1556(mza0,2),r1613(mza0,2),r1616(mza0,2),r8483(mza0,2), &
  r8486(mza0,2),r8484(mza0,2),r8446(mza0,2),r8583(mza0,2),r8586(mza0,2), &
  r8585(mza0,2),r8556(mza0,2),r8683(mza0,2),r8686(mza0,2),r1713(mza0,2), &
  r1717(mza0,2),r8783(mza0,2),r8787(mza0,2),r2332(mza0,2),r2327(mza0,2), &
  r2323(mza0,2),r2337(mza0,2),r2442(mza0,2),r2427(mza0,2),r2424(mza0,2), &
  r2447(mza0,2),r2552(mza0,2),r2527(mza0,2),r2525(mza0,2),r2557(mza0,2), &
  r2662(mza0,2),r2627(mza0,2),r2626(mza0,2),r2667(mza0,2),r2772(mza0,2), &
  r2727(mza0,2),r0000(mza0,2),r1812(mza0,2),r1882(mza0,2)

real, intent(inout) :: &
  e1111(mza0),e1118(mza0),e8888(mza0),e8882(mza0),e1112(mza0), &
  e1811(mza0),e1211(mza0),e8288(mza0),e2222(mza0),e5555(mza0), &
  e6666(mza0),e7777(mza0),e3333(mza0),e3335(mza0),e4444(mza0), &
  e4445(mza0),e3433(mza0),e3444(mza0),e3435(mza0),e3445(mza0), &
  e3533(mza0),e3633(mza0),e3733(mza0),e4544(mza0),e4644(mza0), &
  e4744(mza0),e5655(mza0),e5755(mza0),e6766(mza0),e1413(mza0), &
  e1411(mza0),e1446(mza0),e1513(mza0),e1511(mza0),e1556(mza0), &
  e1613(mza0),e1611(mza0),e8483(mza0),e8488(mza0),e8446(mza0), &
  e8583(mza0),e8588(mza0),e8556(mza0),e8683(mza0),e8688(mza0), &
  e1713(mza0),e1711(mza0),e8783(mza0),e8788(mza0),e2322(mza0), &
  e2327(mza0),e2333(mza0),e2337(mza0),e2422(mza0),e2427(mza0), &
  e2444(mza0),e2447(mza0),e2522(mza0),e2527(mza0),e2555(mza0), &
  e2557(mza0),e2622(mza0),e2627(mza0),e2666(mza0),e2667(mza0), &
  e2722(mza0),e2727(mza0),e2777(mza0),e0000(mza0),e1882(mza0)

real, intent(inout) :: con_ccnx(mza0,nccntyp)

integer :: lcat,jcat,kd1,kd2,k
real :: fac

real ::  rloss(mza0,ncat)
real :: qrloss(mza0,ncat)
real :: enloss(mza0,ncat)

integer :: ibin, jbin, jc
real :: con_ccny(nccntyp)
real :: ccnloss, tot, con_bin

! Limit rxfer and enxfer values so no hydrometeor category gets over-depleted.
! Limit the losses to each category without accounting for potential 
! compensating gains.

! All enxfer and rxfer values are nonnegative.

 rloss = 0.
qrloss = 0.
enloss = 0.

! CLOUD loss adjustments (kg/m^3 and #/m^3)

if (jnmb(1) >= 1) then

   do k = k1(1),k2(1)
      if (rx(k,1) < rxmin(1)) cycle

      rloss(k,1) = r1118(k,1) + r1112(k,1) + r1818(k,1) + r1212(k,1) &
                 + r1413(k,1) + r1416(k,1) + r1414(k,1) + r1513(k,1) &
                 + r1516(k,1) + r1515(k,1) + r1613(k,1) + r1616(k,1) &
                 + r1713(k,1) + r1717(k,1) + r1812(k,1)

      qrloss(k,1) = r1118(k,2) + r1112(k,2) + r1818(k,2) + r1212(k,2) &
                  + r1413(k,2) + r1416(k,2) + r1414(k,2) + r1513(k,2) &
                  + r1516(k,2) + r1515(k,2) + r1613(k,2) + r1616(k,2) &
                  + r1713(k,2) + r1717(k,2) + r1812(k,2)

      if (rloss(k,1) > rx(k,1)) then
         fac = rx(k,1) / max(1.e-20,rloss(k,1))

          rloss(k,1) =     rx(k,1)
         qrloss(k,1) = qrloss(k,1) * fac

         r1118(k,:) = r1118(k,:) * fac
         r1112(k,:) = r1112(k,:) * fac
         r1818(k,:) = r1818(k,:) * fac
         r1212(k,:) = r1212(k,:) * fac
         r1413(k,:) = r1413(k,:) * fac
         r1416(k,:) = r1416(k,:) * fac
         r1414(k,:) = r1414(k,:) * fac
         r1513(k,:) = r1513(k,:) * fac
         r1516(k,:) = r1516(k,:) * fac
         r1515(k,:) = r1515(k,:) * fac
         r1613(k,:) = r1613(k,:) * fac
         r1616(k,:) = r1616(k,:) * fac
         r1713(k,:) = r1713(k,:) * fac
         r1717(k,:) = r1717(k,:) * fac
         r1812(k,:) = r1812(k,:) * fac
      endif

      enloss(k,1) = e1111(k) + e1118(k) + e1112(k) + e1811(k) &
                  + e1211(k) + e1413(k) + e1411(k) + e1513(k) &
                  + e1511(k) + e1613(k) + e1611(k) + e1713(k) &
                  + e1711(k)

      if (enloss(k,1) > cx(k,1)) then
         fac = cx(k,1) / max(1.e-20,enloss(k,1))

         enloss(k,1) = cx(k,1)

         e1118(k) = e1118(k) * fac
         e1112(k) = e1112(k) * fac
         e1413(k) = e1413(k) * fac
         e1513(k) = e1513(k) * fac
         e1613(k) = e1613(k) * fac
         e1713(k) = e1713(k) * fac
      endif
   enddo

endif

! RAIN loss adjustments (kg/m^3 and #/m^3)

if (jnmb(2) >= 1) then

   do k = k1(2),k2(2)
      if (rx(k,2) < rxmin(2)) cycle

      rloss(k,2) = r2327(k,1) + r2323(k,1) + r2427(k,1) + r2424(k,1) &
                 + r2527(k,1) + r2525(k,1) + r2627(k,1) + r2626(k,1) &
                 + r2727(k,1)

      qrloss(k,2) = r2327(k,2) + r2323(k,2) + r2427(k,2) + r2424(k,2) &
                  + r2527(k,2) + r2525(k,2) + r2627(k,2) + r2626(k,2) &
                  + r2727(k,2)

      if (rloss(k,2) > rx(k,2)) then
         fac = rx(k,2) / max(1.e-20,rloss(k,2))

          rloss(k,2) =     rx(k,2)
         qrloss(k,2) = qrloss(k,2) * fac

         r2327(k,:) = r2327(k,:) * fac
         r2323(k,:) = r2323(k,:) * fac
         r2427(k,:) = r2427(k,:) * fac
         r2424(k,:) = r2424(k,:) * fac
         r2527(k,:) = r2527(k,:) * fac
         r2525(k,:) = r2525(k,:) * fac
         r2627(k,:) = r2627(k,:) * fac
         r2626(k,:) = r2626(k,:) * fac
         r2727(k,:) = r2727(k,:) * fac
      endif

      enloss(k,2) = e2222(k) + e2322(k) + e2327(k) + e2422(k) &
                  + e2427(k) + e2522(k) + e2527(k) + e2622(k) &
                  + e2627(k) + e2722(k) + e2727(k)

      if (enloss(k,2) > cx(k,2)) then
         fac = cx(k,2) / max(1.e-20,enloss(k,2))

         enloss(k,2) = cx(k,2)

         e2327(k) = e2327(k) * fac
         e2427(k) = e2427(k) * fac
         e2527(k) = e2527(k) * fac
         e2627(k) = e2627(k) * fac
         e2727(k) = e2727(k) * fac
      endif
   enddo

endif

! PRISTINE ICE loss adjustments (kg/m^3 and #/m^3)

if (jnmb(3) >= 1) then

   do k = k1(3),k2(3)
      if (rx(k,3) < rxmin(3)) cycle

      rloss(k,3) = r3335(k,1) + r3435(k,1) + r3535(k,1) + r3636(k,1) &
                 + r3737(k,1) + r2332(k,1) + r2337(k,1)

      qrloss(k,3) = r3335(k,2) + r3435(k,2) + r3535(k,2) + r3636(k,2) &
                  + r3737(k,2) + r2332(k,2) + r2337(k,2)

      if (rloss(k,3) > rx(k,3)) then
         fac = rx(k,3) / max(1.e-20,rloss(k,3))

          rloss(k,3) =     rx(k,3)
         qrloss(k,3) = qrloss(k,3) * fac

         r3335(k,:) = r3335(k,:) * fac
         r3435(k,:) = r3435(k,:) * fac
         r3535(k,:) = r3535(k,:) * fac
         r3636(k,:) = r3636(k,:) * fac
         r3737(k,:) = r3737(k,:) * fac
         r2332(k,:) = r2332(k,:) * fac
         r2337(k,:) = r2337(k,:) * fac
      endif

      enloss(k,3) = e3333(k) + e3335(k) + e3433(k) + e3435(k) &
                  + e3533(k) + e3633(k) + e3733(k) + e2333(k) &
                  + e2337(k)

      if (enloss(k,3) > cx(k,3)) then
         fac = cx(k,3) / max(1.e-20,enloss(k,3))

         enloss(k,3) = cx(k,3)

         e3335(k) = e3335(k) * fac
         e3435(k) = e3435(k) * fac
         e2337(k) = e2337(k) * fac
      endif
   enddo

endif

! SNOW loss adjustments (kg/m^3 and #/m^3)

if (jnmb(4) >= 1) then

   do k = k1(4),k2(4)
      if (rx(k,4) < rxmin(4)) cycle

      rloss(k,4) = r4445(k,1) + r3445(k,1) + r4545(k,1) + r4646(k,1) &
                 + r4747(k,1) + r1446(k,1) + r8446(k,1) + r2442(k,1) &
                 + r2447(k,1)

      qrloss(k,4) = r4445(k,2) + r3445(k,2) + r4545(k,2) + r4646(k,2) &
                  + r4747(k,2) + r1446(k,2) + r8446(k,2) + r2442(k,2) &
                  + r2447(k,2)

      if (rloss(k,4) > rx(k,4)) then
         fac = rx(k,4) / max(1.e-20,rloss(k,4))

          rloss(k,4) =     rx(k,4)
         qrloss(k,4) = qrloss(k,4) * fac

         r4445(k,:) = r4445(k,:) * fac
         r3445(k,:) = r3445(k,:) * fac
         r4545(k,:) = r4545(k,:) * fac
         r4646(k,:) = r4646(k,:) * fac
         r4747(k,:) = r4747(k,:) * fac
         r1446(k,:) = r1446(k,:) * fac
         r8446(k,:) = r8446(k,:) * fac
         r2442(k,:) = r2442(k,:) * fac
         r2447(k,:) = r2447(k,:) * fac
      endif

      enloss(k,4) = e4444(k) + e4445(k) + e3444(k) + e3445(k) &
                  + e4544(k) + e4644(k) + e4744(k) + e2444(k) &
                  + e2447(k) + e1446(k) + e8446(k)

      if (enloss(k,4) > cx(k,4)) then
         fac = cx(k,4) / max(1.e-20,enloss(k,4))

         enloss(k,4) = cx(k,4)

         e4445(k) = e4445(k) * fac
         e3445(k) = e3445(k) * fac
         e2447(k) = e2447(k) * fac
         e1446(k) = e1446(k) * fac
         e8446(k) = e8446(k) * fac
      endif
   enddo

endif

! AGGREGATES loss adjustments (kg/m^3 and #/m^3)

if (jnmb(5) >= 1) then

   do k = k1(5),k2(5)
      if (rx(k,5) < rxmin(5)) cycle

      rloss(k,5) = r5656(k,1) + r5757(k,1) + r1556(k,1) + r8556(k,1) &
                 + r2552(k,1) + r2557(k,1)

      qrloss(k,5) = r5656(k,2) + r5757(k,2) + r1556(k,2) + r8556(k,2) &
                  + r2552(k,2) + r2557(k,2)

      if (rloss(k,5) > rx(k,5)) then
         fac = rx(k,5) / max(1.e-20,rloss(k,5))

          rloss(k,5) =     rx(k,5)
         qrloss(k,5) = qrloss(k,5) * fac

         r5656(k,:) = r5656(k,:) * fac
         r5757(k,:) = r5757(k,:) * fac
         r1556(k,:) = r1556(k,:) * fac
         r8556(k,:) = r8556(k,:) * fac
         r2552(k,:) = r2552(k,:) * fac
         r2557(k,:) = r2557(k,:) * fac
      endif

      enloss(k,5) = e5555(k) + e5655(k) + e5755(k) + e2555(k) &
                  + e2557(k) + e1556(k) + e8556(k)

      if (enloss(k,5) > cx(k,5)) then
         fac = cx(k,5) / max(1.e-20,enloss(k,5))

         enloss(k,5) = cx(k,5)

         e2557(k) = e2557(k) * fac
         e1556(k) = e1556(k) * fac
         e8556(k) = e8556(k) * fac
      endif
   enddo

endif

! GRAUPEL loss adjustments (kg/m^3 and #/m^3)

if (jnmb(6) >= 1) then

   do k = k1(6),k2(6)
      if (rx(k,6) < rxmin(6)) cycle

      rloss(k,6) = r6767(k,1) + r2662(k,1) + r2667(k,1)

      qrloss(k,6) = r6767(k,2) + r2662(k,2) + r2667(k,2)

      if (rloss(k,6) > rx(k,6)) then
         fac = rx(k,6) / max(1.e-20,rloss(k,6))

          rloss(k,6) =     rx(k,6)
         qrloss(k,6) = qrloss(k,6) * fac

         r6767(k,:) = r6767(k,:) * fac
         r2662(k,:) = r2662(k,:) * fac
         r2667(k,:) = r2667(k,:) * fac
      endif

      enloss(k,6) = e6666(k) + e6766(k) + e2666(k) + e2667(k)

      if (enloss(k,6) > cx(k,6)) then
         fac = cx(k,6) / max(1.e-20,enloss(k,6))

         enloss(k,6) = cx(k,6)

         e2667(k) = e2667(k) * fac
      endif
   enddo

endif

! HAIL loss adjustments (kg/m^3 and #/m^3)

if (jnmb(7) >= 1) then

   do k = k1(7),k2(7)
      if (rx(k,7) < rxmin(7)) cycle

      rloss(k,7) = r2772(k,1)

      qrloss(k,7) = r2772(k,2)

      if (rloss(k,7) > rx(k,7)) then
         fac = rx(k,7) / max(1.e-20,rloss(k,7))

          rloss(k,7) =     rx(k,7)
         qrloss(k,7) = qrloss(k,7) * fac

         r2772(k,:) = r2772(k,:) * fac
      endif

      enloss(k,7) = e7777(k) + e2777(k)

      if (enloss(k,7) > cx(k,7)) then
         fac = cx(k,7) / max(1.e-20,enloss(k,7))

         enloss(k,7) = cx(k,7)

! no enxfer adjustments needed
      endif
   enddo

endif

! DRIZZLE loss adjustments (kg/m^3 and #/m^3)

if (jnmb(8) >= 1) then

   do k = k1(8),k2(8)
      if (rx(k,8) < rxmin(8)) cycle

      rloss(k,8) = r8882(k,1) + r8282(k,1) + r8483(k,1) + r8486(k,1) &
                 + r8484(k,1) + r8583(k,1) + r8586(k,1) + r8585(k,1) &
                 + r8683(k,1) + r8686(k,1) + r8783(k,1) + r8787(k,1) + r1882(k,1)

      qrloss(k,8) = r8882(k,2) + r8282(k,2) + r8483(k,2) + r8486(k,2) &
                  + r8484(k,2) + r8583(k,2) + r8586(k,2) + r8585(k,2) &
                  + r8683(k,2) + r8686(k,2) + r8783(k,2) + r8787(k,2) + r1882(k,2)

      if (rloss(k,8) > rx(k,8)) then
         fac = rx(k,8) / max(1.e-20,rloss(k,8))

          rloss(k,8) =     rx(k,8)
         qrloss(k,8) = qrloss(k,8) * fac

         r8882(k,:) = r8882(k,:) * fac
         r8282(k,:) = r8282(k,:) * fac
         r8483(k,:) = r8483(k,:) * fac
         r8486(k,:) = r8486(k,:) * fac
         r8484(k,:) = r8484(k,:) * fac
         r8583(k,:) = r8583(k,:) * fac
         r8586(k,:) = r8586(k,:) * fac
         r8585(k,:) = r8585(k,:) * fac
         r8683(k,:) = r8683(k,:) * fac
         r8686(k,:) = r8686(k,:) * fac
         r8783(k,:) = r8783(k,:) * fac
         r8787(k,:) = r8787(k,:) * fac
         r1882(k,:) = r1882(k,:) * fac
      endif

      enloss(k,8) = e8888(k) + e8882(k) + e8288(k) + e8483(k) &
                  + e8488(k) + e8583(k) + e8588(k) + e8683(k) &
                  + e8688(k) + e8783(k) + e8788(k) + e1882(k)

 ! NEW

      if (enloss(k,8) > cx(k,8)) then
         fac = cx(k,8) / max(1.e-20,enloss(k,8))

         enloss(k,8) = cx(k,8)

         e8882(k) = e8882(k) * fac
         e8483(k) = e8483(k) * fac
         e8583(k) = e8583(k) * fac
         e8683(k) = e8683(k) * fac
         e8783(k) = e8783(k) * fac
         e1882(k) = e1882(k) * fac
      endif
   enddo

endif

! Transfer bulk density, energy, and number between hydrometeor categories 
! as a result of hydrometeor collisions.

! CLOUD as donor

if (jnmb(1) >= 1) then

   do k = k1(1),k2(1)
      if (rloss(k,1) < 1.e-20) cycle

      rx(k,1) = rx(k,1) -  rloss(k,1)
      qr(k,1) = qr(k,1) - qrloss(k,1)
      cx(k,1) = cx(k,1) - enloss(k,1)

      rx(k,2) = rx(k,2) + r1112(k,1) + r1212(k,1) + r1812(k,1)
      rx(k,3) = rx(k,3) + r1413(k,1) + r1513(k,1) + r1613(k,1) + r1713(k,1)
      rx(k,4) = rx(k,4) + r1414(k,1)
      rx(k,5) = rx(k,5) + r1515(k,1)
      rx(k,6) = rx(k,6) + r1416(k,1) + r1516(k,1) + r1616(k,1)
      rx(k,7) = rx(k,7) + r1717(k,1)
      rx(k,8) = rx(k,8) + r1118(k,1) + r1818(k,1)

      qr(k,2) = qr(k,2) + r1112(k,2) + r1212(k,2) + r1812(k,2)
      qr(k,3) = qr(k,3) + r1413(k,2) + r1513(k,2) + r1613(k,2) + r1713(k,2)
      qr(k,4) = qr(k,4) + r1414(k,2)
      qr(k,5) = qr(k,5) + r1515(k,2)
      qr(k,6) = qr(k,6) + r1416(k,2) + r1516(k,2) + r1616(k,2)
      qr(k,7) = qr(k,7) + r1717(k,2)
      qr(k,8) = qr(k,8) + r1118(k,2) + r1818(k,2)

      cx(k,2) = cx(k,2) + e1112(k)
      cx(k,3) = cx(k,3) + e1413(k) + e1513(k) + e1613(k) + e1713(k)
      cx(k,8) = cx(k,8) + e1118(k)

   enddo

endif

! RAIN as donor

if (jnmb(2) >= 1) then

   do k = k1(2),k2(2)
      if (rloss(k,2) < 1.e-20) cycle

      rx(k,2) = rx(k,2) -  rloss(k,2)
      qr(k,2) = qr(k,2) - qrloss(k,2)
      cx(k,2) = cx(k,2) - enloss(k,2)

      rx(k,3) = rx(k,3) + r2323(k,1)
      rx(k,4) = rx(k,4) + r2424(k,1)
      rx(k,5) = rx(k,5) + r2525(k,1)
      rx(k,6) = rx(k,6) + r2626(k,1) 
      rx(k,7) = rx(k,7) + r2327(k,1) + r2427(k,1) + r2527(k,1) &
                        + r2627(k,1) + r2727(k,1)

      qr(k,3) = qr(k,3) + r2323(k,2)
      qr(k,4) = qr(k,4) + r2424(k,2)
      qr(k,5) = qr(k,5) + r2525(k,2)
      qr(k,6) = qr(k,6) + r2626(k,2) 
      qr(k,7) = qr(k,7) + r2327(k,2) + r2427(k,2) + r2527(k,2) &
                        + r2627(k,2) + r2727(k,2)

      cx(k,7) = cx(k,7) + e2327(k) + e2427(k) + e2527(k) &
                        + e2627(k) + e2727(k)

   enddo

endif

! PRISTINE ICE as donor

if (jnmb(3) >= 1) then

   do k = k1(3),k2(3)
      if (rloss(k,3) < 1.e-20) cycle

      rx(k,3) = rx(k,3) -  rloss(k,3)
      qr(k,3) = qr(k,3) - qrloss(k,3)
      cx(k,3) = cx(k,3) - enloss(k,3)

      rx(k,2) = rx(k,2) + r2332(k,1)
      rx(k,5) = rx(k,5) + r3335(k,1) + r3435(k,1) + r3535(k,1)
      rx(k,6) = rx(k,6) + r3636(k,1) 
      rx(k,7) = rx(k,7) + r3737(k,1) + r2337(k,1)

      qr(k,2) = qr(k,2) + r2332(k,2)
      qr(k,5) = qr(k,5) + r3335(k,2) + r3435(k,2) + r3535(k,2)
      qr(k,6) = qr(k,6) + r3636(k,2) 
      qr(k,7) = qr(k,7) + r3737(k,2) + r2337(k,2)

      cx(k,5) = cx(k,5) + e3335(k) + e3435(k) 
      cx(k,7) = cx(k,7) + e2337(k)
   enddo

endif

! SNOW as donor

if (jnmb(4) >= 1) then

   do k = k1(4),k2(4)
      if (rloss(k,4) < 1.e-20) cycle

      rx(k,4) = rx(k,4) -  rloss(k,4)
      qr(k,4) = qr(k,4) - qrloss(k,4)
      cx(k,4) = cx(k,4) - enloss(k,4)

      rx(k,2) = rx(k,2) + r2442(k,1) 
      rx(k,5) = rx(k,5) + r4445(k,1) + r3445(k,1) + r4545(k,1)
      rx(k,6) = rx(k,6) + r4646(k,1) + r1446(k,1) + r8446(k,1)
      rx(k,7) = rx(k,7) + r4747(k,1) + r2447(k,1)

      qr(k,2) = qr(k,2) + r2442(k,2)
      qr(k,5) = qr(k,5) + r4445(k,2) + r3445(k,2) + r4545(k,2)
      qr(k,6) = qr(k,6) + r4646(k,2) + r1446(k,2) + r8446(k,2)
      qr(k,7) = qr(k,7) + r4747(k,2) + r2447(k,2)

      cx(k,5) = cx(k,5) + e4445(k) + e3445(k)
      cx(k,6) = cx(k,6) + e1446(k) + e8446(k)
      cx(k,7) = cx(k,7) + e2447(k)
   enddo

endif

! AGGREGATES as donor

if (jnmb(5) >= 1) then

   do k = k1(5),k2(5)
      if (rloss(k,5) < 1.e-20) cycle

      rx(k,5) = rx(k,5) -  rloss(k,5)
      qr(k,5) = qr(k,5) - qrloss(k,5)
      cx(k,5) = cx(k,5) - enloss(k,5)

      rx(k,2) = rx(k,2) + r2552(k,1) 
      rx(k,6) = rx(k,6) + r5656(k,1) + r1556(k,1) + r8556(k,1) 
      rx(k,7) = rx(k,7) + r5757(k,1) + r2557(k,1)

      qr(k,2) = qr(k,2) + r2552(k,2) 
      qr(k,6) = qr(k,6) + r5656(k,2) + r1556(k,2) + r8556(k,2) 
      qr(k,7) = qr(k,7) + r5757(k,2) + r2557(k,2)

      cx(k,6) = cx(k,6) + e1556(k) + e8556(k)
      cx(k,7) = cx(k,7) + e2557(k)
   enddo

endif

! GRAUPEL as donor

if (jnmb(6) >= 1) then

   do k = k1(6),k2(6)
      if (rloss(k,6) < 1.e-20) cycle

      rx(k,6) = rx(k,6) -  rloss(k,6)
      qr(k,6) = qr(k,6) - qrloss(k,6)
      cx(k,6) = cx(k,6) - enloss(k,6)

      rx(k,2) = rx(k,2) + r2662(k,1)
      rx(k,7) = rx(k,7) + r6767(k,1) + r2667(k,1)

      qr(k,2) = qr(k,2) + r2662(k,2)
      qr(k,7) = qr(k,7) + r6767(k,2) + r2667(k,2)

      cx(k,7) = cx(k,7) + e2667(k)
   enddo

endif

! HAIL as donor

if (jnmb(7) >= 1) then

   do k = k1(7),k2(7)
      if (rloss(k,7) < 1.e-20) cycle

      rx(k,7) = rx(k,7) -  rloss(k,7)
      qr(k,7) = qr(k,7) - qrloss(k,7)
      cx(k,7) = cx(k,7) - enloss(k,7)

      rx(k,2) = rx(k,2) + r2772(k,1)

      qr(k,2) = qr(k,2) + r2772(k,2)

! no cx recipients
   enddo

endif

! DRIZZLE as donor

if (jnmb(8) >= 1) then

   do k = k1(8),k2(8)
      if (rloss(k,8) < 1.e-20) cycle

      rx(k,8) = rx(k,8) -  rloss(k,8)
      qr(k,8) = qr(k,8) - qrloss(k,8)
      cx(k,8) = cx(k,8) - enloss(k,8)

      rx(k,2) = rx(k,2) + r8882(k,1) + r8282(k,1) + r1882(k,1)
      rx(k,3) = rx(k,3) + r8483(k,1) + r8583(k,1) + r8683(k,1) + r8783(k,1)
      rx(k,4) = rx(k,4) + r8484(k,1)
      rx(k,5) = rx(k,5) + r8585(k,1) 
      rx(k,6) = rx(k,6) + r8486(k,1) + r8586(k,1) + r8686(k,1) 
      rx(k,7) = rx(k,7) + r8787(k,1)

      qr(k,2) = qr(k,2) + r8882(k,2) + r8282(k,2) + r1882(k,2)
      qr(k,3) = qr(k,3) + r8483(k,2) + r8583(k,2) + r8683(k,2) + r8783(k,2)
      qr(k,4) = qr(k,4) + r8484(k,2)
      qr(k,5) = qr(k,5) + r8585(k,2) 
      qr(k,6) = qr(k,6) + r8486(k,2) + r8586(k,2) + r8686(k,2)
      qr(k,7) = qr(k,7) + r8787(k,2)

      cx(k,2) = cx(k,2) + e8882(k) + e1882(k)
      cx(k,3) = cx(k,3) + e8483(k) + e8583(k) + e8683(k) + e8783(k)

   enddo

endif

! Scavenge CCN in accordance with collisional losses of cloud number and
! pristine ice number.  Loop through all bins, using ihyg array to reorder bins
! from lowest to highest critical supersaturation (as used for activation).
! Scavenge CCN in this order until the correct number to scavenge is reached. 

! NOTE: GCCN are scavenged immediately when nucleated to drizzle.
! Pure IFN are not currently scavenged in the model, but CCN that have
! IFN properties are automatically scavenged here.

if (iccn >= 2) then

   do k = k1(11),k2(11)

      ! Define ccnloss as sum of cloud droplet and pristine ice losses due to
      ! collisions.  Perhaps it is better to be more selective in the types of
      ! losses by summing individual contributions like e1411(k), e3533(k), etc.

      ccnloss = enloss(k,1) + enloss(k,3)

      if (ccnloss < 1.e-20) cycle

      ! Save a copy of CCN concentrations before scavenging

      con_ccny(1:nccntyp) = con_ccnx(k,1:nccntyp)

      tot = 0.

      do ibin = 1,nbins
         jbin = ihyg(ibin)
         jc = iccntyp(jbin)
         con_bin = relcon_bin(jbin) * con_ccny(jc)
         tot = tot + con_bin

         if (tot < ccnloss) then
            con_ccnx(k,jc) = con_ccnx(k,jc) - con_bin
         else
            con_ccnx(k,jc) = con_ccnx(k,jc) - con_bin + (tot - ccnloss)
            exit
         endif
      enddo

   enddo

endif

end subroutine colxfers
