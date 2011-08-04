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
subroutine micro()

use mem_ijtabs, only: jtab_w, istp, mrl_endl
use micro_coms, only: level
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer :: j,iw,mrl

if (level /= 3) return

call psub()
!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw)
do j = 1,jtab_w(30)%jend(mrl); iw = jtab_w(30)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   call micphys(iw)

enddo
!$omp end parallel do
endif
call rsub('W',30)

return
end subroutine micro

!===============================================================================

subroutine micphys(iw)

use micro_coms,  only: mza0, ncat, cfvt, jnmb, emb2, rxmin, neff, nhcat, emb1
use consts_coms, only: pi4
use misc_coms,   only: io6, dtlm
use mem_ijtabs,  only: itab_w
use mem_grid,    only: lpw,glatw,glonw

implicit none

integer, intent(in) :: iw

integer :: k,jflag,jcat,lcat,j1,j2

integer :: j12,j14,j15,j16,j17,j18,j23,j24,j25,j26,j27
integer :: j35,j36,j37,j45,j46,j47,j56,j57,j67,j82,j84,j85,j86,j87
integer :: k12,k14,k15,k16,k17,k18,k23,k24,k25,k26,k27
integer :: k35,k36,k37,k45,k46,k47,k56,k57,k67,k82,k84,k85,k86,k87

integer :: lpw0,iw0,mrl0

real :: pcpg0,qpcpg0,dpcpg0,dtl0,dtli0

real :: pi4dt  ! delta_t * pi * 4

real :: frac

! Automatic arrays

integer :: ict1 (mza0,ncat)
integer :: ict2 (mza0,ncat)
integer :: jhcat(mza0,ncat)

real :: wct1(mza0,ncat)
real :: wct2(mza0,ncat)
real :: rx  (mza0,ncat)
real :: cx  (mza0,ncat)
real :: qx  (mza0,ncat)
real :: qr  (mza0,ncat)
real :: emb (mza0,ncat)
real :: tx  (mza0,ncat)
real :: vap (mza0,ncat)
real :: sb  (mza0,ncat)
real :: sd  (mza0,ncat)
real :: se  (mza0,ncat)
real :: sf  (mza0,ncat)
real :: sg  (mza0,ncat)
real :: sh  (mza0,ncat)
real :: sm  (mza0,ncat)
real :: ss  (mza0,ncat)
real :: su  (mza0,ncat)
real :: sw  (mza0,ncat)
real :: sy  (mza0,ncat)
real :: sz  (mza0,ncat)

real :: eff(mza0,neff)

real :: colfac (mza0)
real :: colfac2(mza0)
real :: dsed_thil(mza0)

real :: rxfer (mza0,ncat,ncat)
real :: qrxfer(mza0,ncat,ncat)
real :: enxfer(mza0,ncat,ncat)

integer :: k1(11)
integer :: k2(11)
integer :: k3(11)

real :: thil0 (mza0)
real :: theta0(mza0)
real :: press0(mza0)
real :: exner0(mza0)
real :: wc0   (mza0)

real :: tair     (mza0)
real :: tairc    (mza0)
real :: tairstrc (mza0)
real :: rhovstr  (mza0)
real :: rhov     (mza0)
real :: rhoi     (mza0)
real :: rhovslair(mza0)
real :: rhovsiair(mza0)
real :: thrmcon  (mza0)
real :: vapdif   (mza0)
real :: dynvisc  (mza0)
real :: rdynvsci (mza0)
real :: denfac   (mza0)
real :: sumuy    (mza0)
real :: sumuz    (mza0)
real :: sumvr    (mza0)
real :: con_ccnx (mza0)
real :: con_gccnx(mza0)
real :: con_ifnx (mza0)
real :: totcond  (mza0)

real(kind=8) :: rhoa(mza0)
real(kind=8) :: rhow(mza0)

real :: tref     (mza0,2)
real :: rhovsref (mza0,2)
real :: rhovsrefp(mza0,2)

real :: sa(mza0,9)

real :: pcprx(ncat)
real :: accpx(ncat)

real :: ch1(nhcat)  ! delta_t * fall speed coefficient (automatic array)

! Set constants for this column

iw0  = iw
lpw0 = lpw(iw0)
mrl0 = itab_w(iw)%mrlw
dtl0 = dtlm(mrl0)

dtli0 = 1. / dtl0
pi4dt = pi4 * dtl0

ch1(2:nhcat) = dtl0 * cfvt(2:nhcat)

! Zero out microphysics scratch arrays for the present i,j column

ict1 (:,:) = 0
ict2 (:,:) = 0
jhcat(:,:) = 0

wct1(:,:) = 0.
wct2(:,:) = 0.
rx  (:,:) = 0.
cx  (:,:) = 0.
qx  (:,:) = 0.
qr  (:,:) = 0.
emb (:,:) = 0.
tx  (:,:) = 0.
vap (:,:) = 0.
sb  (:,:) = 0.
sd  (:,:) = 0.
se  (:,:) = 0.
sf  (:,:) = 0.
sg  (:,:) = 0.
sh  (:,:) = 0.
sm  (:,:) = 0.
ss  (:,:) = 0.
su  (:,:) = 0.
sw  (:,:) = 0.
sy  (:,:) = 0.
sz  (:,:) = 0.

eff(:,:) = 0.

colfac (:) = 0.
colfac2(:) = 0.
dsed_thil(:) = 0.

rxfer (:,:,:) = 0.
qrxfer(:,:,:) = 0.
enxfer(:,:,:) = 0.

k1(:) = 0
k2(:) = 0
k3(:) = 0

thil0 (:) = 0.
theta0(:) = 0.
press0(:) = 0.
exner0(:) = 0.
wc0   (:) = 0.

tair     (:) = 0.
tairc    (:) = 0.
tairstrc (:) = 0.
rhovstr  (:) = 0.
rhov     (:) = 0.
rhoi     (:) = 0.
rhovslair(:) = 0.
rhovsiair(:) = 0.
thrmcon  (:) = 0.
vapdif   (:) = 0.
dynvisc  (:) = 0.
rdynvsci (:) = 0.
denfac   (:) = 0.
sumuy    (:) = 0.
sumuz    (:) = 0.
sumvr    (:) = 0.
con_ccnx (:) = 0.
con_gccnx(:) = 0.
con_ifnx (:) = 0.
totcond  (:) = 0.

rhoa(:) = 0.
rhow(:) = 0.

tref     (:,:) = 0.
rhovsref (:,:) = 0.
rhovsrefp(:,:) = 0.

sa(:,:) = 0.

pcprx(:) = 0.
accpx(:) = 0.

! Copy hydrometeor bulk mass and number concentration from main model arrays
! to microphysics column arrays

call mic_copy(iw0,lpw0,thil0,press0,wc0,rhoa,rhow,rhoi, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx,qr)

! Loop over all vertical levels

do k = lpw0,mza0

! Compute total condensate in k level

   totcond(k) = 1.001 * (rx(k,1) + rx(k,2) + rx(k,3) + rx(k,4) &
                       + rx(k,5) + rx(k,6) + rx(k,7) + rx(k,8))

! If total water exceeds condensate, no corrections are necessary

   if (real(rhow(k)) > totcond(k)) cycle

! If total water density is negative, increase to zero and increase total air density
! rhoa by same amount to preserve dry-air content

   if (real(rhow(k)) < 0.) then
      rhoa(k) = rhoa(k) - rhow(k)
      rhow(k) = 0.
   endif

! Adjust condensate amounts downward if their sum exceeds rhow 

   if (totcond(k) > real(rhow(k))) then

      frac = rhow(k) / totcond(k)

      do lcat = 1,ncat
         rx(k,lcat) = rx(k,lcat) * frac
         cx(k,lcat) = cx(k,lcat) * frac
         qr(k,lcat) = qr(k,lcat) * frac
      enddo

   endif
enddo

! Find minimum and maximum model level in this column for each species

do lcat = 1,ncat

! Find new k2(lcat)

   k = mza0
   do while (k > lpw0 .and. rx(k,lcat) < rxmin(lcat))
      k = k - 1
   enddo
   k2(lcat) = k

! Find new k1(lcat)

   k = lpw0
   do while (k <= k2(lcat) .and. rx(k,lcat) < rxmin(lcat))
      k = k + 1
   enddo
   k1(lcat) = k

enddo

! Save initial k2 values in k3 for copyback

k3(1) = k2(1)
k3(3) = k2(3)
k3(8) = k2(8)

! Min/max heights for any liquid, any ice, and any of either

k1(9)  = min(k1(1),k1(2),k1(8))
k2(9)  = max(k2(1),k2(2),k2(8))
k1(10) = min(k1(3),k1(4),k1(5),k1(6),k1(7))
k2(10) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
k1(11) = min(k1(9),k1(10))
k2(11) = max(k2(9),k2(10))

call thrmstr(iw0,lpw0,k1,k2, &
   press0,thil0,rhow,rhoi,exner0,tair,theta0,rhov,rhovstr,tairstrc, &
   rx,qx,sa)

call each_column(lpw0,iw0,k1,k2,dtl0,                          &
   jhcat,press0,tair,tairc,tairstrc,rhovstr,rhoa,rhov,rhoi,    &
   rhovslair,rhovsiair,thrmcon,vapdif,dynvisc,rdynvsci,denfac, &
   colfac,colfac2,sumuy,sumuz,sumvr,                           &
   tx,sh,sm,sa,tref,rhovsref,rhovsrefp)


! USEFUL SAMPLE PRINTS FOR EXAMINING A SINGLE IW0 COLUMN BEFORE
! MAIN MICROPHYSICS COMPUTATIONS ARE DONE

!print*, 'micphys1 ',iw0,glatw(iw0),glonw(iw0)

!if (iw0 == 2443) then
!   print*, ' '
!   do k = lpw0,mza0
!      write(6,'(a,i6,10e12.3)') 'rx ',k,(rx(k,j1),j1=1,8)
!   enddo

!   print*, ' '
!   do k = lpw0,mza0
!      write(6,'(a,i6,10e12.3)') 'cx ',k,(cx(k,j1),j1=1,8)
!   enddo

!   print*, ' '
!   do k = lpw0,mza0
!      write(6,'(a,i6,10e12.3)') 'qr ',k,(qr(k,j1),j1=1,8)
!   enddo

!   print*, ' '
!   do k = lpw0,mza0
!      write(6,'(a,i6,10e12.3)') 'env ',k,press0(k),theta0(k),wc0(k),tair(k), &
!                                tairc(k),rhoa(k),rhow(k),rhov(k)
!   enddo
!endif

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

jflag = 1

if (k1(1) <= k2(1)) &
   call enemb(1,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(2) <= k2(2)) &
   call enemb(2,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(3) <= k2(3)) &
   call enemb(3,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(4) <= k2(4)) &
   call enemb(4,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(5) <= k2(5)) &
   call enemb(5,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(6) <= k2(6)) &
   call enemb(6,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(7) <= k2(7)) &
   call enemb(7,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

! Set up matrix for heat/vapor diffusion computation

if (k1(1) <= k2(1)) &
   call diffprep(iw0,1,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(2) <= k2(2)) &
   call diffprep(iw0,2,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(3) <= k2(3)) &
   call diffprep(iw0,3,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(4) <= k2(4)) &
   call diffprep(iw0,4,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(5) <= k2(5)) &
   call diffprep(iw0,5,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(6) <= k2(6)) &
   call diffprep(iw0,6,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(7) <= k2(7)) &
   call diffprep(iw0,7,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

if (k1(8) <= k2(8)) &
   call diffprep(iw0,8,k1,k2, &
      pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
      rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

! Implicit matrix solution of atmospheric vapor density

call vapdiff(k1(11),k2(11),rhov,rhovstr,sumuy,sumuz)

! Vapor flux applied to each category.  Do not change the order of these

if (k1(1) <= k2(1)) &
   call vapflux(iw0,1,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(8) <= k2(8)) &
   call vapflux(iw0,8,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(3) <= k2(3)) &
   call vapflux(iw0,3,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(4) <= k2(4)) &
   call vapflux(iw0,4,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(5) <= k2(5)) &
   call vapflux(iw0,5,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(2) <= k2(2)) &
   call vapflux(iw0,2,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(6) <= k2(6)) &
   call vapflux(iw0,6,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

if (k1(7) <= k2(7)) &
   call vapflux(iw0,7,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

! Conversion between pristine ice and snow due to vapor flux

if (jnmb(4) >= 1) then
   k1(3) = min(k1(3),k1(4))
   k2(3) = max(k2(3),k2(4))
   k1(4) = k1(3)
   k2(4) = k2(3)
   if (k1(3) <= k2(3)) call psxfer(iw0,k1(3),k2(3),jhcat,vap,rx,cx,qx,qr)
endif

! Diagnose new air temperature following heat and vapor fluxes

call newtemp(k1(11),k2(11), &
   tairstrc,rhoi,rhovstr,rhov,exner0,tairc,tair,theta0,rhovslair,rhovsiair,sa)

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

jflag = 2

if (k1(1) <= k2(1)) &
   call enemb(1,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(2) <= k2(2)) &
   call enemb(2,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(3) <= k2(3)) &
   call enemb(3,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(4) <= k2(4)) &
   call enemb(4,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(5) <= k2(5)) &
   call enemb(5,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(6) <= k2(6)) &
   call enemb(6,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(7) <= k2(7)) &
   call enemb(7,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

! Determine coalescence efficiencies

call effxy(lpw0,k1,k2,rx,qr,emb,tx,eff)

! Liquid water collisions

j18 = max(k1(1),k1(8)); k18 = min(k2(1),k2(8))
j82 = max(k1(8),k1(2)); k82 = min(k2(8),k2(2))
j12 = max(k1(1),k1(2)); k12 = min(k2(1),k2(2))

! Self-collection of cloud and drizzle

if (jnmb(8) >= 1 .and. k1(1) <= k2(1)) &
   call col1188(1,8,1,k1(1),k2(1),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

if (jnmb(2) >= 1 .and. k1(8) <= k2(8)) &
   call col1188(8,2,1,k1(8),k2(8),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

! If drizzle is NOT active but rain is

if (jnmb(8) == 0 .and. k1(1) <= k2(1)) &
   call col1188(1,2,1,k1(1),k2(1),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

! Collisions between different liquid species

if (j18 <= k18) &
   call col1(1,8,8,1,j18,k18, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j12 <= k12) &
   call col1(1,2,2,1,j12,k12, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j82 <= k82) &
   call col1(8,2,2,1,j82,k82, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)
   
! Self collection of rain, aggregates, graupel, and hail: number change only

if (jnmb(2) == 5 .and. k1(2) <= k2(2)) &
   call cols(2,1,k1(2),k2(2),jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

if (jnmb(5) == 5 .and. k1(5) <= k2(5)) &
   call cols(5,4,k1(5),k2(5),jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

if (jnmb(6) == 5 .and. k1(6) <= k2(6)) &
   call cols(6,5,k1(6),k2(6),jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

if (jnmb(7) == 5 .and. k1(7) <= k2(7)) &
   call cols(7,7,k1(7),k2(7),jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,enxfer)

! Self collection of pristine ice and of snow

if (jnmb(5) >= 1 .and. k1(3) <= k2(3)) &
   call col3344(3,5,2,k1(3),k2(3),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

if (jnmb(5) >= 1 .and. k1(4) <= k2(4)) &
   call col3344(4,5,3,k1(4),k2(4),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,rxfer,qrxfer,enxfer)

! Collection between pristine ice and snow (k1,k2 now same for lcats 3 & 4)

if (jnmb(5) >= 1 .and. k1(3) <= k2(3)) &
   call col3443(2,k1(3),k2(3),iw0,     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

! Other ice-ice collisions

j35 = max(k1(3),k1(5)); k35 = min(k2(3),k2(5))
j36 = max(k1(3),k1(6)); k36 = min(k2(3),k2(6))
j37 = max(k1(3),k1(7)); k37 = min(k2(3),k2(7))
j45 = max(k1(4),k1(5)); k45 = min(k2(4),k2(5))
j46 = max(k1(4),k1(6)); k46 = min(k2(4),k2(6))
j47 = max(k1(4),k1(7)); k47 = min(k2(4),k2(7))
j56 = max(k1(5),k1(6)); k56 = min(k2(5),k2(6))
j57 = max(k1(5),k1(7)); k57 = min(k2(5),k2(7))
j67 = max(k1(6),k1(7)); k67 = min(k2(6),k2(7))

if (j35 <= k35) &
   call col1(3,5,5,2,j35,k35, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j36 <= k36) &
   call col1(3,6,6,5,j36,k36, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j37 <= k37) &
   call col1(3,7,7,6,j37,k37, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j45 <= k45) &
   call col1(4,5,5,3,j45,k45, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j46 <= k46) &
   call col1(4,6,6,5,j46,k46, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j47 <= k47) &
   call col1(4,7,7,6,j47,k47, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j56 <= k56) &
   call col1(5,6,6,5,j56,k56, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j57 <= k57) &
   call col1(5,7,7,6,j57,k57, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

if (j67 <= k67) &
   call col1(6,7,7,5,j67,k67, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

! Ice-cloud and ice-drizzle collisions with graupel by-product

if (jnmb(6) >= 1) then
   j14 = max(k1(1),k1(4)); k14 = min(k2(1),k2(4))
   j15 = max(k1(1),k1(5)); k15 = min(k2(1),k2(5))
   j16 = max(k1(1),k1(6)); k16 = min(k2(1),k2(6))

   j84 = max(k1(8),k1(4)); k84 = min(k2(8),k2(4))
   j85 = max(k1(8),k1(5)); k85 = min(k2(8),k2(5))
   j86 = max(k1(8),k1(6)); k86 = min(k2(8),k2(6))
   
   if (j14 <= k14) &
      call col2(1,4,6,1,j14,k14, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j15 <= k15) &
      call col2(1,5,6,1,j15,k15, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j16 <= k16) &
      call col2(1,6,6,1,j16,k16, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j84 <= k84) &
      call col2(8,4,6,1,j84,k84, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j85 <= k85) &
      call col2(8,5,6,1,j85,k85, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j86 <= k86) &
      call col2(8,6,6,1,j86,k86, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

endif

! Hail-cloud and hail-drizzle collisions (hail by-product)

if (jnmb(7) >= 1) then

   j17 = max(k1(1),k1(7)); k17 = min(k2(1),k2(7))
   j87 = max(k1(8),k1(7)); k87 = min(k2(8),k2(7))
   
   if (j17 <= k17) &
      call col2(1,7,7,1,j17,k17, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

   if (j87 <= k87) &
      call col2(8,7,7,1,j87,k87, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac, &
         rxfer,qrxfer,enxfer,dtl0)

! Ice-rain collisions (hail byproduct)

   j23 = max(k1(2),k1(3)); k23= min(k2(2),k2(3))
   j24 = max(k1(2),k1(4)); k24= min(k2(2),k2(4))
   j25 = max(k1(2),k1(5)); k25= min(k2(2),k2(5))
   j26 = max(k1(2),k1(6)); k26= min(k2(2),k2(6))
   j27 = max(k1(2),k1(7)); k27= min(k2(2),k2(7))

   if (j23 <= k23) &
      call col3(3,7,1,j23,k23, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

   if (j24 <= k24) &
      call col3(4,7,1,j24,k24, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

   if (j25 <= k25) &
      call col3(5,7,1,j25,k25, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

   if (j26 <= k26) &
      call col3(6,7,1,j26,k26, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)

   if (j27 <= k27) &
      call col3(7,7,1,j27,k27, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,rxfer,qrxfer,enxfer)
   
endif

! Apply transfers of bulk mass, energy, and number from all collisions

call colxfers(k1,k2,rx,cx,qr,rxfer,qrxfer,enxfer)

! Do not change order of the following x02 calls

if (jnmb(3) >= 1) &
   call x02(iw0,3,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(1) >= 1) &
   call x02(iw0,1,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(8) >= 1) &
   call x02(iw0,8,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(4) >= 1) &
   call x02(iw0,4,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(5) >= 1) &
   call x02(iw0,5,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(6) >= 1) &
   call x02(iw0,6,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(7) >= 1) &
   call x02(iw0,7,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

if (jnmb(2) >= 1) &
   call x02(iw0,2,k1,k2, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

! Nucleation of cloud droplets

if (jnmb(1) >= 1) &
   call cldnuc(iw0,lpw0,dtli0,rx,cx,con_ccnx,con_gccnx,rhov,rhoi,rhoa, &
               tair,tairc,wc0,rhovslair)

! Rediagnose k2(1) and k3(1) because of possible new cloud nucleation

k = mza0
do while (k > lpw0 .and. rx(k,1) < rxmin(1))
   k = k - 1
enddo
k2(1) = k
k3(1) = max(k2(1),k3(1))

! Rediagnose k1(1) because of possible new cloud nucleation

k = lpw0
do while (k <= k2(1) .and. rx(k,1) < rxmin(1))
   k = k + 1
enddo
k1(1) = k

! Rediagnose k2(8) and k3(8) because of possible new cloud nucleation

k = mza0
do while (k > lpw0 .and. rx(k,8) < rxmin(8))
   k = k - 1
enddo
k2(8) = k
k3(8) = max(k2(8),k3(8))

! Rediagnose k1(8) because of possible new cloud nucleation

k = lpw0
do while (k <= k2(8) .and. rx(k,8) < rxmin(8))
   k = k + 1
enddo
k1(8) = k

! Re-diagnose cloud droplet and drizzle mean mass

jflag = 1

if (jnmb(1) >= 3 .and. k1(1) <= k2(1)) &
   call enemb(1,jflag,k1,k2,ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (jnmb(8) >= 3 .and. k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2,ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

! Nucleation of ice crystals

if (jnmb(3) >= 1) &
   call icenuc(k1,k2,lpw0,mrl0,iw0, &
      rx,cx,qr,emb,vap,tx,rhov,rhoa,press0,dynvisc,thrmcon, &
      tair,tairc,rhovslair,rhovsiair,con_ccnx,con_ifnx,dtl0)

! Rediagnose k2(3) and k3(3) because of possible new ice nucleation

k = mza0
do while (k > lpw0 .and. rx(k,3) < rxmin(3))
   k = k - 1
enddo
k2(3) = k
k3(3) = max(k2(3),k3(3))

! Rediagnose k1(3) because of possible new ice nucleation

k = lpw0
do while (k <= k2(3) .and. rx(k,3) < rxmin(3))
   k = k + 1
enddo
k1(3) = k

! Re-diagnose pristine ice and cloud droplet mean mass

jflag = 1

if (k2(3) >= k1(3)) &
   call enemb(3,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k2(1) >= k1(1)) &
   call enemb(1,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

if (k2(8) >= k1(8)) &
   call enemb(8,jflag,k1,k2, &
      ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

! Zero out surface precip arrays (for leaf, will need to coordinate w/new leaf4)

pcpg0  = 0.
qpcpg0 = 0.
dpcpg0 = 0.

accpx(1:ncat) = 0.
pcprx(1:ncat) = 0.

! Zero out array for accumulating sedim changes to thil 

dsed_thil(lpw0:mza0) = 0.

! Compute sedimentation for all 7 precipitating categories

if (k1(2) <= k2(2)) &
   call sedim(2,iw0,lpw0,k1,k2,.001, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(3) <= k2(3)) &
   call sedim(3,iw0,lpw0,k1,k2,.010, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(4) <= k2(4))&
   call sedim(4,iw0,lpw0,k1,k2,.010, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(5) <= k2(5)) &
   call sedim(5,iw0,lpw0,k1,k2,.010,  &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(6) <= k2(6)) &
   call sedim(6,iw0,lpw0,k1,k2,.003, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(7) <= k2(7)) &
   call sedim(7,iw0,lpw0,k1,k2,.001, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

if (k1(8) <= k2(8)) &
   call sedim(8,iw0,lpw0,k1,k2,.001, &
      dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
      rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

! Apply change to thil from sedim of all categories

thil0(lpw0:mza0) = thil0(lpw0:mza0) + dsed_thil(lpw0:mza0)

! Copy hydrometeor bulk mass and number concentration and surface precipitation
! from microphysics column arrays to main model arrays

call mic_copyback(iw0,lpw0,k1,k2,k3, &
   pcpg0,qpcpg0,dpcpg0,dtli0,accpx,pcprx,thil0,theta0,rhoa,rhow,rhov,rhoi, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx)

return
end subroutine micphys

!===============================================================================

subroutine mic_copy(iw0,lpw0,thil0,press0,wc0,rhoa,rhow,rhoi, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx,qr)

use micro_coms, only: mza0, ncat, jnmb, rxmin
use mem_basic,  only: thil, press, wc, rho, sh_w, sh_v

use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                      q2, q6, q7, con_ccn, con_gccn, con_ifn, &
                      con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h

use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

real, intent(out) :: thil0 (mza0)
real, intent(out) :: press0(mza0)
real, intent(out) :: wc0   (mza0)
real, intent(out) :: rhoi  (mza0)
real, intent(out) :: con_ccnx (mza0)
real, intent(out) :: con_gccnx(mza0)
real, intent(out) :: con_ifnx (mza0)

real(kind=8), intent(out) :: rhoa(mza0)
real(kind=8), intent(out) :: rhow(mza0)

real, intent(inout) :: rx(mza0,ncat)
real, intent(inout) :: cx(mza0,ncat)
real, intent(inout) :: qx(mza0,ncat)
real, intent(inout) :: qr(mza0,ncat)

integer :: k

! copy atmospheric variables to micphys column vectors

do k = lpw0,mza0
   thil0 (k) = thil (k,iw0)
   press0(k) = press(k,iw0)
   wc0   (k) = wc   (k,iw0)
   rhoa  (k) = rho  (k,iw0)
   rhow  (k) = sh_w (k,iw0) * rho(k,iw0)

   rhoi  (k) = 1. / rhoa(k)
enddo

! Cloud water

if (jnmb(1) >= 5) then

   do k = lpw0,mza0

! Make sure cloud water number concentration is positive-definite

      con_c(k,iw0) = max(con_c(k,iw0), 0.0)

! If cloud water bulk density is sufficiently abundant, copy to rx
! and copy cloud water number concentration to cx.
! Otherwise, make sure cloud water bulk density is positive-definite
! and set cloud water number concentration to zero.

      if (sh_c(k,iw0) >= rxmin(1)) then
         rx(k,1) = sh_c (k,iw0) * rhoa(k)
         cx(k,1) = con_c(k,iw0) * rhoa(k)
      else
         sh_c (k,iw0) = max(sh_c(k,iw0), 0.0)
         con_c(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(1) >= 1) then

   do k = lpw0,mza0

! If cloud water bulk density is sufficiently abundant, copy to rx.
! Otherwise, if cloud water bulk density is negative, set to zero.

      if (sh_c(k,iw0) >= rxmin(1)) then
         rx(k,1) = sh_c(k,iw0) * rhoa(k)
      elseif (sh_c(k,iw0) < 0.) then
         sh_c(k,iw0) = 0.
      endif

   enddo

endif

! Rain

if (jnmb(2) >= 5) then

   do k = lpw0,mza0

! Make sure rain number concentration is positive-definite

      con_r(k,iw0) = max(con_r(k,iw0), 0.0)

! If rain bulk density is sufficiently abundant, copy to rx
! and copy rain number concentration to cx.  Fill qx and compute qr.
! Otherwise, make sure rain bulk density is positive-definite
! and set rain number concentration to zero.

      if (sh_r(k,iw0) >= rxmin(2)) then
         rx(k,2) = sh_r (k,iw0) * rhoa(k)
         cx(k,2) = con_r(k,iw0) * rhoa(k)
         qx(k,2) = q2(k,iw0)
         qr(k,2) = qx(k,2) * rx(k,2)
      else
         sh_r (k,iw0) = max(sh_r(k,iw0), 0.0)
         con_r(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(2) >= 1) then

   do k = lpw0,mza0

! If rain bulk density is sufficiently abundant, copy to rx.
! Fill qx and compute qr.
! Otherwise, if rain bulk density is negative, set to zero.

      if (sh_r(k,iw0) >= rxmin(2)) then
         rx(k,2) = sh_r(k,iw0) * rhoa(k)
         qx(k,2) = q2(k,iw0)
         qr(k,2) = qx(k,2) * rx(k,2)
      elseif (sh_r(k,iw0) < 0.) then
         sh_r(k,iw0) = 0.
      endif

   enddo

endif

! Pristine ice

if (jnmb(3) >= 5) then

   do k = lpw0,mza0

! Make sure pristine ice number concentration is positive-definite

      con_p(k,iw0) = max(con_p(k,iw0), 0.0)

! If pristine ice bulk density is sufficiently abundant, copy to rx
! and copy pristine ice number concentration to cx.
! Otherwise, make sure pristine ice bulk density is positive-definite
! and set pristine ice number concentration to zero.

      if (sh_p(k,iw0) >= rxmin(3)) then
         rx(k,3) = sh_p (k,iw0) * rhoa(k)
         cx(k,3) = con_p(k,iw0) * rhoa(k)
      else
         sh_p (k,iw0) = max(sh_p(k,iw0), 0.0)
         con_p(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(3) >= 1) then

   do k = lpw0,mza0

! If pristine ice bulk density is sufficiently abundant, copy to rx.
! Otherwise, if pristine ice bulk density is negative, set to zero.

      if (sh_p(k,iw0) >= rxmin(3)) then
         rx(k,3) = sh_p(k,iw0) * rhoa(k)
      elseif (sh_p(k,iw0) < 0.) then
         sh_p(k,iw0) = 0.
      endif

   enddo

endif

! Snow

if (jnmb(4) >= 5) then

   do k = lpw0,mza0

! Make sure snow number concentration is positive-definite

      con_s(k,iw0) = max(con_s(k,iw0), 0.0)

! If snow bulk density is sufficiently abundant, copy to rx
! and copy snow number concentration to cx.
! Otherwise, make sure snow bulk density is positive-definite
! and set snow number concentration to zero.

      if (sh_s(k,iw0) >= rxmin(4)) then
         rx(k,4) = sh_s (k,iw0) * rhoa(k)
         cx(k,4) = con_s(k,iw0) * rhoa(k)
      else
         sh_s (k,iw0) = max(sh_s(k,iw0), 0.0)
         con_s(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(4) >= 1) then

   do k = lpw0,mza0

! If snow bulk density is sufficiently abundant, copy to rx.
! Otherwise, if snow bulk density is negative, set to zero.

      if (sh_s(k,iw0) >= rxmin(4)) then
         rx(k,4) = sh_s(k,iw0) * rhoa(k)
      elseif (sh_s(k,iw0) < 0.) then
         sh_s(k,iw0) = 0.
      endif

   enddo

endif

! Aggregates

if (jnmb(5) >= 5) then

   do k = lpw0,mza0

! Make sure aggregates number concentration is positive-definite

      con_a(k,iw0) = max(con_a(k,iw0), 0.0)

! If aggregates bulk density is sufficiently abundant, copy to rx
! and copy aggregates number concentration to cx.
! Otherwise, make sure aggregates bulk density is positive-definite
! and set aggregates number concentration to zero.

      if (sh_a(k,iw0) >= rxmin(5)) then
         rx(k,5) = sh_a (k,iw0) * rhoa(k)
         cx(k,5) = con_a(k,iw0) * rhoa(k)
      else
         sh_a (k,iw0) = max(sh_a(k,iw0), 0.0)
         con_a(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(5) >= 1) then

   do k = lpw0,mza0

! If aggregates bulk density is sufficiently abundant, copy to rx.
! Otherwise, if aggregates bulk density is negative, set to zero.

      if (sh_a(k,iw0) >= rxmin(5)) then
         rx(k,5) = sh_a(k,iw0) * rhoa(k)
      elseif (sh_a(k,iw0) < 0.) then
         sh_a(k,iw0) = 0.
      endif

   enddo

endif

! Graupel

if (jnmb(6) >= 5) then

   do k = lpw0,mza0

! Make sure graupel number concentration is positive-definite

      con_g(k,iw0) = max(con_g(k,iw0), 0.0)

! If graupel bulk density is sufficiently abundant, copy to rx
! and copy graupel number concentration to cx.  Fill qx and compute qr. 
! Otherwise, make sure graupel bulk density is positive-definite
! and set graupel number concentration to zero.

      if (sh_g(k,iw0) >= rxmin(6)) then
         rx(k,6) = sh_g (k,iw0) * rhoa(k)
         cx(k,6) = con_g(k,iw0) * rhoa(k)
         qx(k,6) = q6(k,iw0)
         qr(k,6) = qx(k,6) * rx(k,6)
      else
         sh_g (k,iw0) = max(sh_g(k,iw0), 0.0)
         con_g(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(6) >= 1) then

   do k = lpw0,mza0

! If graupel bulk density is sufficiently abundant, copy to rx.
! Fill qx and compute qr.
! Otherwise, if graupel bulk density is negative, set to zero.

      if (sh_g(k,iw0) >= rxmin(6)) then
         rx(k,6) = sh_g(k,iw0) * rhoa(k)
         qx(k,6) = q6(k,iw0)
         qr(k,6) = qx(k,6) * rx(k,6)
      elseif (sh_g(k,iw0) < 0.) then
         sh_g(k,iw0) = 0.
      endif

   enddo

endif

! Hail

if (jnmb(7) >= 5) then

   do k = lpw0,mza0

! Make sure hail number concentration is positive-definite

      con_h(k,iw0) = max(con_h(k,iw0), 0.0)

! If hail bulk density is sufficiently abundant, copy to rx
! and copy hail number concentration to cx.  Fill qx and compute qr. 
! Otherwise, make sure hail bulk density is positive-definite
! and set hail number concentration to zero.

      if (sh_h(k,iw0) >= rxmin(7)) then
         rx(k,7) = sh_h (k,iw0) * rhoa(k)
         cx(k,7) = con_h(k,iw0) * rhoa(k)
         qx(k,7) = q7(k,iw0)
         qr(k,7) = qx(k,7) * rx(k,7)
      else
         sh_h (k,iw0) = max(sh_h(k,iw0), 0.0)
         con_h(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(7) >= 1) then

   do k = lpw0,mza0

! If hail bulk density is sufficiently abundant, copy to rx.
! Fill qx and compute qr.
! Otherwise, if hail bulk density is negative, set to zero.

      if (sh_h(k,iw0) >= rxmin(7)) then
         rx(k,7) = sh_h(k,iw0) * rhoa(k)
         qx(k,7) = q7(k,iw0)
         qr(k,7) = qx(k,7) * rx(k,7)
      elseif (sh_h(k,iw0) < 0.) then
         sh_h(k,iw0) = 0.
      endif

   enddo

endif

! Drizzle

if (jnmb(8) >= 5) then

   do k = lpw0,mza0

! Make sure drizzle number concentration is positive-definite

      con_d(k,iw0) = max(con_d(k,iw0), 0.0)

! If drizzle bulk density is sufficiently abundant, copy to rx
! and copy drizzle number concentration to cx.
! Otherwise, make sure drizzle bulk density is positive-definite
! and set drizle number concentration to zero.

      if (sh_d(k,iw0) >= rxmin(8)) then
         rx(k,8) = sh_d (k,iw0) * rhoa(k)
         cx(k,8) = con_d(k,iw0) * rhoa(k)
      else
         sh_d (k,iw0) = max(sh_d(k,iw0), 0.0)
         con_d(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(8) >= 1) then

   do k = lpw0,mza0

! If drizzle bulk density is sufficiently abundant, copy to rx.
! Otherwise, if drizzle bulk density is negative, set to zero.

      if (sh_d(k,iw0) >= rxmin(8)) then
         rx(k,8) = sh_d(k,iw0) * rhoa(k)
      elseif (sh_d(k,iw0) < 0.) then
         sh_d(k,iw0) = 0.
      endif

   enddo

endif

! Copy CCN number concentration to con_ccnx 

if (jnmb(1) == 7) then
   do k = lpw0,mza0
      con_ccnx(k) = con_ccn(k,iw0) * rhoa(k)
   enddo
endif

! Copy GCCN number concentration to con_gccnx 

if (jnmb(8) == 7) then
   do k = lpw0,mza0
      con_gccnx(k) = con_gccn(k,iw0) * rhoa(k)
   enddo
endif

! Copy IFN number concentration to con_ifnx 

if (jnmb(3) == 7) then
   do k = lpw0,mza0
      con_ifnx(k) = con_ifn(k,iw0) * rhoa(k)
   enddo
endif

return
end subroutine mic_copy

!===============================================================================

subroutine mic_copyback(iw0,lpw0,k1,k2,k3, &
   pcpg0,qpcpg0,dpcpg0,dtli0,accpx,pcprx,thil0,theta0,rhoa,rhow,rhov,rhoi, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx)

use micro_coms, only: mza0, ncat, jnmb
use mem_basic,  only: thil, theta, rho, sh_w, sh_v
use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                      q2, q6, q7, pcpgr, qpcpgr, dpcpgr, &
                      con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                      accpd, accpr, accpp, accps, accpa, accpg, accph, &
                      pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                      con_ccn, con_gccn, con_ifn

use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)
integer, intent(in) :: k3(11)

real, intent(in) :: pcpg0
real, intent(in) :: qpcpg0
real, intent(in) :: dpcpg0
real, intent(in) :: dtli0

real, intent(in) :: accpx(ncat)
real, intent(in) :: pcprx(ncat)

real, intent(in)  :: thil0(mza0)
real, intent(in)  :: theta0(mza0)
real, intent(in)  :: rhov(mza0)
real, intent(out) :: rhoi(mza0)

real(kind=8), intent(in) :: rhoa(mza0)
real(kind=8), intent(in) :: rhow(mza0)

real, intent(in) :: con_ccnx (mza0)
real, intent(in) :: con_gccnx(mza0)
real, intent(in) :: con_ifnx (mza0)

real, intent(in) :: rx(mza0,ncat)
real, intent(in) :: cx(mza0,ncat)
real, intent(in) :: qx(mza0,ncat)

integer :: k,kk

! Convert surface precipitation from amount to rate

pcpgr (iw0) = pcpg0 * dtli0
qpcpgr(iw0) = qpcpg0 * dtli0
dpcpgr(iw0) = dpcpg0 * dtli0

do k = lpw0,mza0
   thil(k,iw0)  = thil0(k)
   theta(k,iw0) = theta0(k)
   rho(k,iw0)   = rhoa(k)
   rhoi(k)      = 1. / rhoa(k)
   sh_w(k,iw0)  = rhow(k) * rhoi(k)
   sh_v(k,iw0)  = min(rhov(k) * rhoi(k), sh_w(k,iw0))
enddo

! Copy cloud water bulk density back to main array

if (jnmb(1) >= 1) then
   do k = lpw0,k3(1)
      sh_c(k,iw0) = rx(k,1) * rhoi(k)
   enddo
endif

! Copy rain bulk density, internal energy, and surface precip back to main arrays

if (jnmb(2) >= 1) then
   accpr(iw0) = accpr(iw0) + accpx(2)
   pcprr(iw0) = pcprx(2)

   do k = lpw0,k2(11)
      sh_r(k,iw0) = rx(k,2) * rhoi(k)
      q2  (k,iw0) = qx(k,2)
   enddo
endif

! Copy pristine ice bulk density and surface precip back to main arrays

if (jnmb(3) >= 1) then
   accpp(iw0) = accpp(iw0) + accpx(3)
   pcprp(iw0) = pcprx(3)

   do k = lpw0,k3(3)
      sh_p(k,iw0) = rx(k,3) * rhoi(k)
   enddo
endif

! Copy snow bulk density and surface precip back to main arrays

if (jnmb(4) >= 1) then
   accps(iw0) = accps(iw0) + accpx(4)
   pcprs(iw0) = pcprx(4)

   do k = lpw0,k2(11)
      sh_s(k,iw0) = rx(k,4) * rhoi(k)
   enddo
endif

! Copy aggregates bulk density and surface precip back to main arrays

if (jnmb(5) >= 1) then
   accpa(iw0) = accpa(iw0) + accpx(5)
   pcpra(iw0) = pcprx(5)

   do k = lpw0,k2(11)
      sh_a(k,iw0) = rx(k,5) * rhoi(k)
   enddo
endif

! Copy graupel bulk density, internal energy, and surface precip back to main arrays

if (jnmb(6) >= 1) then
   accpg(iw0) = accpg(iw0) + accpx(6)
   pcprg(iw0) = pcprx(6)
   do k = lpw0,k2(11)
      sh_g(k,iw0) = rx(k,6) * rhoi(k)
      q6  (k,iw0) = qx(k,6)
   enddo
endif

! Copy hail bulk density, internal energy, and surface precip back to main arrays

if (jnmb(7) >= 1) then
   accph(iw0) = accph(iw0) + accpx(7)
   pcprh(iw0) = pcprx(7)
   do k = lpw0,k2(11)
      sh_h(k,iw0) = rx(k,7) * rhoi(k)
      q7  (k,iw0) = qx(k,7)
   enddo
endif

! Copy drizzle bulk density and surface precip back to main arrays

if (jnmb(8) >= 1) then
   accpd(iw0) = accpd(iw0) + accpx(8)
   pcprd(iw0) = pcprx(8)

   do k = lpw0,k2(11)
      sh_d(k,iw0) = rx(k,8) * rhoi(k)
   enddo
endif

! Copy cloud water number concentration back to main array

if (jnmb(1) >= 5) then
   do k = lpw0,k2(11)
      con_c(k,iw0) = cx(k,1) * rhoi(k)
   enddo
endif

! Copy rain number concentration back to main array

if (jnmb(2) == 5) then
   do k = lpw0,k2(11)
      con_r(k,iw0) = cx(k,2) * rhoi(k)
   enddo
endif

! Copy pristine ice number concentration back to main array

if (jnmb(3) >= 5) then
   do k = lpw0,k2(11)
      con_p(k,iw0) = cx(k,3) * rhoi(k)
   enddo
endif

! Copy snow number concentration back to main array

if (jnmb(4) == 5) then
   do k = lpw0,k2(11)
      con_s(k,iw0) = cx(k,4) * rhoi(k)
   enddo
endif

! Copy aggregates number concentration back to main array

if (jnmb(5) == 5) then
   do k = lpw0,k2(11)
      con_a(k,iw0) = cx(k,5) * rhoi(k)
   enddo
endif

! Copy graupel number concentration back to main array

if (jnmb(6) == 5) then
   do k = lpw0,k2(11)
      con_g(k,iw0) = cx(k,6) * rhoi(k)
   enddo
endif

! Copy hail number concentration back to main array

if (jnmb(7) == 5) then
   do k = lpw0,k2(11)
      con_h(k,iw0) = cx(k,7) * rhoi(k)
   enddo
endif

! Copy drizzle number concentration back to main array

if (jnmb(8) >= 5) then
   do k = lpw0,k2(11)
      con_d(k,iw0) = cx(k,8) * rhoi(k)
   enddo
endif

! Copy CCN number concentration back to main array

if (jnmb(1) == 7) then
   do k = lpw0,mza0
      con_ccn(k,iw0) = con_ccnx(k) * rhoi(k)
   enddo
endif

! Copy GCCN number concentration back to main array

if (jnmb(8) == 7) then
   do k = lpw0,mza0
      con_gccn(k,iw0) = con_gccnx(k) * rhoi(k)
   enddo
endif

! Copy IFN number concentration back to main array

if (jnmb(3) == 7) then
   do k = lpw0,mza0
      con_ifn(k,iw0) = con_ifnx(k) * rhoi(k)
   enddo
endif

return
end subroutine mic_copyback

