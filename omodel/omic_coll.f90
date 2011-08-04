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
subroutine effxy(lpw0,k1,k2,rx,qr,emb,tx,eff)

use micro_coms, only: mza0, ncat, jnmb, rxmin, neff
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

integer :: k

real :: dmr

! This subroutine sets COALLESCENCE EFFICIENCIES for all hydrometeor collisions.
! Some of these depend on hydrometeor temperatures, while others are constant.

! COLLISION efficiencies are not considered here, but are all set when the
! collection tables are generated.

! 1 = cc,cd,cr,   cs,ca,cg,ch,
!        dd,dr,   ds,da,dg,dh,
!           rr,rp,rs,ra,rg,rh

eff(lpw0:mza0,1) = 1.0

! 2 = pp,ps,pa

if (jnmb(5) >= 1) then

   do k = k1(3),k2(3)

      if (abs(tx(k,3) + 14.) <= 2.) then
         eff(k,2) = 1.4
      else
         eff(k,2) = min(.2,10. ** (.035 * tx(k,3) - .7))
      endif

   enddo

! 3 = ss,sa

   do k = k1(4),k2(4)

      if (abs(tx(k,4) + 14.) <= 2.) then
         eff(k,3) = 1.4
      else
         eff(k,3) = min(.2,10. ** (.035 * tx(k,4) - .7))
      endif

   enddo

! 4 = aa

   do k = k1(5),k2(5)

      if (abs(tx(k,5) + 14.) <= 2.) then
         eff(k,4) = 1.4
      elseif (tx(k,5) >= -1.) then
         eff(k,4) = 1.
      else
         eff(k,4) = min(.2,10. ** (.035 * tx(k,5) - .7))
      endif

   enddo

endif

! 5 = pg,sg,ag,gg,gh

if (jnmb(6) >= 1) then

   do k = k1(6),k2(6)

      if (qr(k,6) > 0.) then
         eff(k,5) = 1.0
      else
         eff(k,5) = min(.2,10. ** (.035 * tx(k,6) - .7))
      endif

   enddo

endif

! 6 = ph,sh,ah

if (jnmb(7) >= 1) then

   do k = k1(7),k2(7)

      if (qr(k,7) > 0.) then
         eff(k,6) = 1.0
      else
         eff(k,6) = min(.2,10. ** (.035 * tx(k,7) - .7))
      endif

! 7 = hh (experimental; should be tested and improved; large hail probably
!         should not coallesce)

      eff(k,7) = max(0.,.1 + .005 * tx(k,7))

   enddo

endif

return
end subroutine effxy

!===============================================================================

subroutine cols(mx,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairc, coltabc, neff
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

real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: ipc,k

real :: tabc,colc

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle

   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = colfac(k) * eff(k,meff) * (cx(k,mx) ** 2) * (10. ** (-tabc))

   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colc)

enddo

return
end subroutine cols

!===============================================================================

subroutine col1188(mx,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, jnmb, neff,  &
                      coltabc, coltabr
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

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,ipr,ipc,ipz

real :: c1,tabc,tabrx,tabcz,colc,colrx,colcz

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle
   
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
   ipr = ipairr(jhcat(k,mx),jhcat(k,mx))
   ipz = ipairr(jhcat(k,mz),jhcat(k,mx)) ! X to Z number xfer from coltabr

   c1 = eff(k,meff) * colfac(k) * cx(k,mx) ** 2

! Interpolate from coltabr

   tabrx =      wct1(k,mx) * wct1(k,mx) * coltabr(ict1(k,mx),ict1(k,mx),ipr) &
         + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ipr) &
         +      wct2(k,mx) * wct2(k,mx) * coltabr(ict2(k,mx),ict2(k,mx),ipr)

   colrx = min(rx(k,mx),c1 * 10. ** (-tabrx))

! Hydrometeor mass transfer

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colrx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colrx * qx(k,mx)

   if (jnmb(mz) < 5) cycle

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = min(0.5 * cx(k,mx),c1 * 10. ** (-tabc))

! Hydrometeor number loss

   enxfer(k,mx,mx) = enxfer(k,mx,mx) + colc

! Interpolate number from special section of coltabr

   tabcz =      wct1(k,mx) * wct1(k,mx) * coltabr(ict1(k,mx),ict1(k,mx),ipz) &
         + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ipz) &
         +      wct2(k,mx) * wct2(k,mx) * coltabr(ict2(k,mx),ict2(k,mx),ipz)

   colcz = min(0.5 * cx(k,mx),c1 * 10. ** (-tabcz))

! Hydrometeor number transfer

   enxfer(k,mx,mz) = enxfer(k,mx,mz) + colcz

enddo

return
end subroutine col1188

!===============================================================================

subroutine col3344(mx,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, jnmb, neff,  &
                      coltabc, coltabr
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

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,ipr,ipc

real :: c1,tabc,tabrx,colc,colrx

do k = j1,j2

   if (rx(k,mx) < rxmin(mx)) cycle
   
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
   ipr = ipairr(jhcat(k,mx),jhcat(k,mx))

   c1 = eff(k,meff) * colfac2(k) * cx(k,mx) ** 2

! Interpolate from coltabr

   tabrx =      wct1(k,mx) * wct1(k,mx) * coltabr(ict1(k,mx),ict1(k,mx),ipr) &
         + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ipr) &
         +      wct2(k,mx) * wct2(k,mx) * coltabr(ict2(k,mx),ict2(k,mx),ipr)

   colrx = min(rx(k,mx),c1 * 10. ** (-tabrx))

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colrx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colrx * qx(k,mx)

   if (jnmb(mz) < 5) cycle

! Interpolate from coltabc

   tabc =      wct1(k,mx) * wct1(k,mx) * coltabc(ict1(k,mx),ict1(k,mx),ipc) &
        + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc) &
        +      wct2(k,mx) * wct2(k,mx) * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colc = min(0.5 * cx(k,mx),c1 * 10. ** (-tabc))

   enxfer(k,mx,mz) = enxfer(k,mx,mz) + colc
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + colc

enddo

return
end subroutine col3344

!===============================================================================

subroutine col3443(meff,j1,j2,iw0, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, neff, ipairc, ipairr, coltabc, coltabr
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

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,jhcatx,jhcaty,iprxy,ipryx,ipc

real :: c1,tabc,tabrx,tabry,colc,colrx,colry

do k = j1,j2

   if (rx(k,3) < rxmin(3) .or. rx(k,4) < rxmin(4)) cycle

   jhcatx = jhcat(k,3)
   jhcaty = jhcat(k,4)

   ipc   = ipairc(jhcatx,jhcaty)
   iprxy = ipairr(jhcatx,jhcaty)
   ipryx = ipairr(jhcaty,jhcatx)

   c1 = eff(k,meff) * colfac(k) * cx(k,3) * cx(k,4)

! Interpolate from coltabr

   tabrx = wct1(k,3) * wct1(k,4) * coltabr (ict1(k,3),ict1(k,4),iprxy) &
         + wct2(k,3) * wct1(k,4) * coltabr (ict2(k,3),ict1(k,4),iprxy) &
         + wct1(k,3) * wct2(k,4) * coltabr (ict1(k,3),ict2(k,4),iprxy) &
         + wct2(k,3) * wct2(k,4) * coltabr (ict2(k,3),ict2(k,4),iprxy)

   colrx = min(rx(k,3),c1 * 10. ** (-tabrx))

! Interpolate from coltabr

   tabry = wct1(k,4) * wct1(k,3) * coltabr (ict1(k,4),ict1(k,3),ipryx) &
         + wct2(k,4) * wct1(k,3) * coltabr (ict2(k,4),ict1(k,3),ipryx) &
         + wct1(k,4) * wct2(k,3) * coltabr (ict1(k,4),ict2(k,3),ipryx) &
         + wct2(k,4) * wct2(k,3) * coltabr (ict2(k,4),ict2(k,3),ipryx)

   colry = min(rx(k,4),c1 * 10. ** (-tabry))

   rxfer(k,3,5) = rxfer(k,3,5) + colrx
   qrxfer(k,3,5) = qrxfer(k,3,5) + colrx * qx(k,3)

   rxfer(k,4,5) = rxfer(k,4,5) + colry
   qrxfer(k,4,5) = qrxfer(k,4,5) + colry * qx(k,4)

! Interpolate from coltabc

   tabc = wct1(k,3) * wct1(k,4) * coltabc (ict1(k,3),ict1(k,4),ipc) &
        + wct2(k,3) * wct1(k,4) * coltabc (ict2(k,3),ict1(k,4),ipc) &
        + wct1(k,3) * wct2(k,4) * coltabc (ict1(k,3),ict2(k,4),ipc) &
        + wct2(k,3) * wct2(k,4) * coltabc (ict2(k,3),ict2(k,4),ipc)

   colc = c1 * 10. ** (-tabc)

   if (cx(k,3) > cx(k,4)) then
      enxfer(k,4,5) = enxfer(k,4,5) + min(cx(k,4),colc)
      enxfer(k,3,3) = enxfer(k,3,3) + min(cx(k,3),colc)
   else
      enxfer(k,3,5) = enxfer(k,3,5) + min(cx(k,3),colc)
      enxfer(k,4,4) = enxfer(k,4,4) + min(cx(k,4),colc)
   endif

! also loss for aerosol

enddo

return
end subroutine col3443

!===============================================================================

subroutine col1(mx,my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, neff, rxmin, ipairr, ipairc, jnmb, &
                      coltabc, coltabr
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

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,iprxy,ipc

real :: c1,tabc,tabrx,colc,colrx

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

   ipc   = ipairc(jhcat(k,mx),jhcat(k,my))
   iprxy = ipairr(jhcat(k,mx),jhcat(k,my))

   c1 = eff(k,meff) * colfac(k) * cx(k,mx) * cx(k,my)

! Interpolate from coltabr

   tabrx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),iprxy) &
         + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),iprxy) &
         + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),iprxy) &
         + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),iprxy)

   colrx = min(rx(k,mx),c1 * 10. ** (-tabrx))

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colrx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colrx * qx(k,mx)

   if (jnmb(mx) < 5) cycle

! Interpolate from coltabc

   tabc = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc) &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc) &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc) &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

   colc = c1 * 10. ** (-tabc)

   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colc,cx(k,mx))

! also loss for aerosol

enddo

return
end subroutine col1

!===============================================================================

subroutine col2(mx,my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
   rxfer,qrxfer,enxfer,dtl0)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, coltabr, coltabc, &
                      sipfac, pwmasi, emb1, gamsip13, gamsip24, emb0, neff, jnmb
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
real, intent(in) :: emb (mza0,ncat)
real, intent(in) :: eff (mza0,neff)

real, intent(in) :: colfac(mza0)

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

real, intent(in) :: dtl0

integer :: k,jhcatx,jhcaty,iprxy,ipryx,ipc,it

real :: c1,c2,tabc,tabrx,tabry,colc,colrx,colry,colc0,colqrx,colqry, &
        rcoal,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice, &
        area,cn13,cn24,sip,rsip,qrsip,rfinlz,xtoz

real :: alpha(16)
real :: beta(16)

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00.,00.,1.0,1.0,1.0,1.0,00.,00.,00.,00.,00.,1.0,1.0,1.0,1.0/
data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,00.,1.2,1.1,1.1,1.3/

do k = j1,j2

   if (rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle

   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)

   ipc   = ipairc(jhcatx,jhcaty)
   iprxy = ipairr(jhcatx,jhcaty)
   ipryx = ipairr(jhcaty,jhcatx)

   c2 = colfac(k) * cx(k,mx) * cx(k,my)
   c1 = eff(k,meff) * c2

! Interpolate from coltabr

   tabrx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),iprxy) &
         + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),iprxy) &
         + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),iprxy) &
         + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),iprxy)

   colrx = min(rx(k,mx),c1 * 10. ** (-tabrx))

! Interpolate from coltabr

   tabry = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipryx) &
         + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipryx) &
         + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipryx) &
         + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipryx)

   colry = min(rx(k,my),c1 * 10. ** (-tabry))

! Interpolate from coltabc

   tabc = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc) &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc) &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc) &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

   colc0 = c2 * 10. ** (-tabc)
   colc = colc0 * eff(k,meff)

   colqrx = colrx * qx(k,mx)
   colqry = colry * qx(k,my)

   rcoal = colrx + colry
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
! prognostic for at least one of the two hydromet species involved. This is
! specified above in the calculations for "colc".

   if (tcoal > -8.0 .and. tcoal < -3.0 .and. &
      (jnmb(mx) >= 5 .or. jnmb(my) >= 5)) then  ! Steve, why this condition?

      area = cx(k,my) * sipfac(jhcaty) * emb(k,my) ** (2.*pwmasi(jhcaty)) ! remd rhoa
      it = nint(emb(k,mx) / emb1(mx) * 5000.)

      if (mx == 1) then
         cn13 = colc * gamsip13(1,it) / (area * dtl0)
         cn24 = min(cx(k,mx),colc0) * gamsip24(1,it)
      elseif (mx == 8) then
         cn13 = colc * gamsip13(2,it) / (area * dtl0)
         cn24 = min(cx(k,mx),colc0) * gamsip24(2,it)
      endif

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

      enxfer(k,mx,3) = enxfer(k,mx,3) + sip
      rxfer (k,mx,3) = rxfer (k,mx,3) + rsip
      qrxfer(k,mx,3) = qrxfer(k,mx,3) + qrsip

   endif

! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
! fracliq may be a little larger than fracx due to collected
! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
! or 1.2, or need ALPHA itself to be at least 1.0.

   rfinlz = min(rcoal,alpha(jhcaty) * coalliq + beta(jhcaty) * colrx)

   xtoz = min(colrx,rfinlz)

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
   rxfer(k,mx,my) = rxfer(k,mx,my) + colrx - xtoz
   if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
   qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (colrx - xtoz)
   if (my /= mz) qrxfer(k,my,mz) = qrxfer(k,my,mz) + qx(k,my) * (rfinlz - xtoz)

   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colc,cx(k,mx))
   if (my /= mz) enxfer(k,my,mz) = enxfer(k,my,mz)  &
      + (rfinlz - xtoz) * min(colc,cx(k,my)) / max(1.e-20,colry)
      
! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y

! also include loss of aerosol

enddo

return
end subroutine col2

!===============================================================================

subroutine col3(my,mz,meff,j1,j2, &
   jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, rxmin, ipairr, ipairc, jnmb, neff, &
                      coltabr, coltabc
use misc_coms,  only: io6

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

real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: k,iprxy,ipryx,ipc,jhcaty

real :: c1,tabc,tabrx,tabry,colc,colrx,colry,colcx,colcy,colqrx,colqry, &
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

   ipc   = ipairc(jhcat(k,2),jhcaty)
   iprxy = ipairr(jhcat(k,2),jhcaty)
   ipryx = ipairr(jhcaty,jhcat(k,2))

   c1 = eff(k,meff) * colfac(k) * cx(k,2) * cx(k,my)

! Interpolate from coltabr

   tabrx = wct1(k,2) * wct1(k,my) * coltabr(ict1(k,2),ict1(k,my),iprxy) &
         + wct2(k,2) * wct1(k,my) * coltabr(ict2(k,2),ict1(k,my),iprxy) &
         + wct1(k,2) * wct2(k,my) * coltabr(ict1(k,2),ict2(k,my),iprxy) &
         + wct2(k,2) * wct2(k,my) * coltabr(ict2(k,2),ict2(k,my),iprxy)

   colrx = min(rx(k,2),c1 * 10. ** (-tabrx))

! Interpolate from coltabr

   tabry = wct1(k,my) * wct1(k,2) * coltabr(ict1(k,my),ict1(k,2),ipryx) &
         + wct2(k,my) * wct1(k,2) * coltabr(ict2(k,my),ict1(k,2),ipryx) &
         + wct1(k,my) * wct2(k,2) * coltabr(ict1(k,my),ict2(k,2),ipryx) &
         + wct2(k,my) * wct2(k,2) * coltabr(ict2(k,my),ict2(k,2),ipryx)

   colry = min(rx(k,my),c1 * 10. ** (-tabry))

   if (jnmb(2) == 5) then

! Interpolate from coltabc

      tabc = wct1(k,2) * wct1(k,my) * coltabc(ict1(k,2),ict1(k,my),ipc) &
           + wct2(k,2) * wct1(k,my) * coltabc(ict2(k,2),ict1(k,my),ipc) &
           + wct1(k,2) * wct2(k,my) * coltabc(ict1(k,2),ict2(k,my),ipc) &
           + wct2(k,2) * wct2(k,my) * coltabc(ict2(k,2),ict2(k,my),ipc)

      colc = c1 * 10. ** (-tabc)

      colcx = min(cx(k,2),colc)
      colcy = min(cx(k,my),colc)

      coalnum = min(colcx,colcy)

   endif

   colqrx = colrx * qx(k,2)
   colqry = colry * qx(k,my)

   rcoal  = colrx + colry
   qrcoal = colqrx + colqry
   qcoal  = qrcoal / (1.e-20 + rcoal)

   call qtc(qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

   if (fracliq >= .99) then

      rxfer(k,my,2) = rxfer(k,my,2) + colry
      qrxfer(k,my,2) = qrxfer(k,my,2) + colqry

      if (jnmb(2) == 5) enxfer(k,my,my) = enxfer(k,my,my) + colcy

   else

      rfinlz = min(rcoal, alpha(jhcaty) * coalliq + beta(jhcaty) * colrx)

      xtoz = min(colrx,rfinlz)

      rxfer(k,2,mz) = rxfer(k,2,mz) + xtoz
      rxfer(k,2,my) = rxfer(k,2,my) + colrx - xtoz

      if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

! NEED TO USE QCOAL TO TRANSFER Q?

      qrxfer(k,2,mz) = qrxfer(k,2,mz) + qx(k,2) * xtoz
      qrxfer(k,2,my) = qrxfer(k,2,my) + qx(k,2) * (colrx - xtoz)

      if (my /= mz) &
         qrxfer(k,my,mz) = qrxfer(k,my,mz) + qx(k,my) * (rfinlz - xtoz)

      if (jnmb(2) == 5) then

         if (my == mz) then

            enxfer(k,2,2) = enxfer(k,2,2) + colcx

         elseif (colcy >= colcx) then

            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)

            enxfer(k,2,mz) = enxfer(k,2,mz) + cfinlz
            enxfer(k,2,2) = enxfer(k,2,2) + colcx - cfinlz
            enxfer(k,my,my) = enxfer(k,my,my) + colcy

         else

            cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)

            enxfer(k,my,mz) = enxfer(k,my,mz) + cfinlz
            enxfer(k,2,2) = enxfer(k,2,2) + colcx
            enxfer(k,my,my) = enxfer(k,my,my) + colcy - cfinlz

         endif

      endif

   endif

enddo

! also include loss of aerosol

return
end subroutine col3

!===============================================================================

subroutine colxfers(k1,k2,rx,cx,qr,rxfer,qrxfer,enxfer)

use micro_coms, only: mza0, ncat, jnmb
use misc_coms,  only: io6

implicit none

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(inout) :: rx    (mza0,ncat)
real, intent(inout) :: cx    (mza0,ncat)
real, intent(inout) :: qr    (mza0,ncat)
real, intent(inout) :: rxfer (mza0,ncat,ncat)
real, intent(inout) :: qrxfer(mza0,ncat,ncat)
real, intent(inout) :: enxfer(mza0,ncat,ncat)

integer :: lcat,jcat,kd1,kd2,k

! Automatic arrays

real :: rloss (mza0)
real :: enloss(mza0)

! Limit rxfer and enxfer values so no hydrometeor category gets over-depleted.
! Limit the losses to each category without accounting for potential 
! compensating gains.

! All enxfer and rxfer values are nonnegative.

! rxfer(i,j) and qrxfer(i,j) are zero for i = j.

! Loop over donor categories

do lcat = 1,ncat

   if (jnmb(lcat) >= 1) then

      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do k = kd1,kd2
         rloss(k) = 0.
         enloss(k) = 0.
      enddo

! Loop over recipient categories

      do jcat = 1,ncat

         if (jnmb(jcat) >= 1) then

! Sum losses from donor category

            do k = kd1,kd2
               rloss (k) = rloss (k) + rxfer (k,lcat,jcat) ! kg/m^3
               enloss(k) = enloss(k) + enxfer(k,lcat,jcat) ! #/m^3
            enddo

         endif

      enddo

! Convert rloss and enloss to scaling factors that will reduce, if necessary,
! rxfer and enxfer values in order to prevent over-depletion of donor category.

      do k = kd1,kd2
         rloss(k) = min(1.,rx(k,lcat) / max(1.e-20,rloss(k)))
         enloss(k) = min(1.,cx(k,lcat) / max(1.e-10,enloss(k)))
      enddo

! Apply possible reduction to rxfer, qrxfer, and enxfer arrays.

! Loop over recipient categories

      do jcat = 1,ncat

         if (jnmb(jcat) >= 1) then

            do k = kd1,kd2
               rxfer (k,lcat,jcat) = rxfer (k,lcat,jcat) * rloss (k)
               qrxfer(k,lcat,jcat) = qrxfer(k,lcat,jcat) * rloss (k)
               enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat) * enloss(k)
            enddo

         endif

      enddo

   endif

enddo

! Transfer bulk density, energy, and number between hydrometeor categories 
! as a result of hydrometeor collisions.

! Loop over donor categories

do lcat = 1,ncat

   if (jnmb(lcat) >= 1) then

      kd1 = k1(lcat)
      kd2 = k2(lcat)

! Loop over recipient categories

      do jcat = 1,ncat

         if (jnmb(jcat) >= 1) then

            do k = kd1,kd2
               rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
               rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)

               qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
               qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)

               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
               cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
            enddo

         endif

      enddo

   endif

enddo

return
end subroutine colxfers


