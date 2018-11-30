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
subroutine micro()

use mem_ijtabs, only: jtab_w, istp, mrl_endl, jtw_prog
use micro_coms, only: miclevel
use misc_coms,  only: io6, time8
use oname_coms, only: nl
use mem_grid,   only: mwa, glatw, glonw

implicit none

integer :: j,iw,mrl,nbc

! The nbincall array stores the number of grid levels where subroutine ccnbin
! is called for each IW.  It is for informational (printing) purposes only and
! can be removed from the model once the performance of ccnbin is better
! understood.

integer :: nbincall(mwa)

if (miclevel /= 3) return

! Section for parcel model simulation (only one grid cell is prognosed)

if (nl%test_case >= 901 .and. nl%test_case <= 949) then

! Specify IW of grid cell to prognose (k=2 assumed for parcel cell)

   iw = 3

   call parcel_env(iw)
   call micphys(iw,nbincall(iw))

   RETURN
endif

!----------------------------------------------------------------------
mrl = mrl_endl(istp)
if (mrl > 0) then
!$omp parallel do private (iw) schedule(guided)
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

   call micphys(iw,nbincall(iw))

enddo
!$omp end parallel do
endif

!nbc = 0
!do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!   nbc = nbc + nbincall(iw)
!enddo
!
!print*, 'nbincall ',nbc

end subroutine micro

!===============================================================================

subroutine micphys(iw,nbincall)

use micro_coms,  only: mza0, ncat, nhcat, neff, cfvt, jnmb, emb2, rxmin, emb1
use ccnbin_coms, only: nccntyp, nnuc
use consts_coms, only: r8, pi4
use misc_coms,   only: io6, dtlm, time_istp8, timmax8
use mem_ijtabs,  only: itab_w
use mem_grid,    only: lpw,glatw,glonw
use oname_coms,  only: nl

implicit none

integer, intent(in) :: iw
integer, intent(inout) :: nbincall

integer :: k,jflag,jcat,lcat,j1,j2

integer :: j12,j14,j15,j16,j17,j18,j23,j24,j25,j26,j27
integer :: j35,j36,j37,j45,j46,j47,j56,j57,j67,j82,j84,j85,j86,j87
integer :: k12,k14,k15,k16,k17,k18,k23,k24,k25,k26,k27
integer :: k35,k36,k37,k45,k46,k47,k56,k57,k67,k82,k84,k85,k86,k87

integer :: lpw0,iw0,mrl0,kend

real :: dtl0,dtli0

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
real :: dmb (mza0,ncat)
real :: pcpvel  (mza0,ncat)
real :: pcpfluxc(mza0,ncat)
real :: pcpfluxr(mza0,ncat)
real :: pcpfluxq(mza0,ncat)
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

integer :: k1(11)
integer :: k2(11)
integer :: k3(11)

real :: voa   (mza0)
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
real :: con_ccnx (mza0,nccntyp)
real :: con_gccnx(mza0)
real :: con_ifnx (mza0)
real :: totcond  (mza0)

! New arrays (2/25/2012) to retain mass and number amounts of each
! cloud and ice nucleation type (for plotting)

real :: rnuc_vc       (mza0)
real :: rnuc_vd       (mza0)
real :: rnuc_cp_hom   (mza0)
real :: rnuc_dp_hom   (mza0)
real :: rnuc_vp_haze  (mza0)
real :: rnuc_vp_immers(mza0)

real :: cnuc_vc       (mza0)
real :: cnuc_vd       (mza0)
real :: cnuc_cp_hom   (mza0)
real :: cnuc_dp_hom   (mza0)
real :: cnuc_vp_haze  (mza0)
real :: cnuc_vp_immers(mza0)

! New transfer arrays for subroutine psxfer

real :: rpsxfer(mza0)
real :: epsxfer(mza0)

! Replacement of rxfer and enxfer arrays (2/25/2012) to retain transfer
! information (for plotting) from each separate collection pair interaction.

! The new arrays have names 'r' for mass and heat transfer, and 'e' for number
! transfer, followed by 4 digts.  The first 2 digits represent the 2 colliding
! categories, the 3rd digit represents the donor category, and the 4th digit
! represents the recipient category.  For example, the consequence of a
! collision between rain and pristine ice is a mass transfer from rain to hail
! that is denoted as 'r2327', which uses species numbers: 1=cloud, 2=rain, 
! 3=pristine_ice, 4=snow, 5=aggregates, 6=graupel, 7=hail, 8=drizzle (=cloud2).

real :: &
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

real :: &
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

real(r8) :: rhoa(mza0)
real(r8) :: rhow(mza0)

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
dmb (:,:) = 0.
pcpvel  (:,:) = 0.
pcpfluxc(:,:) = 0.
pcpfluxr(:,:) = 0.
pcpfluxq(:,:) = 0.
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

k1(:) = 0
k2(:) = 0
k3(:) = 0

voa   (:) = 0.
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
con_ccnx (:,:) = 0.
con_gccnx(:) = 0.
con_ifnx (:) = 0.
totcond  (:) = 0.

 rnuc_vc       (:) = 0.
 rnuc_vd       (:) = 0.
 rnuc_cp_hom   (:) = 0.
 rnuc_dp_hom   (:) = 0.
 rnuc_vp_haze  (:) = 0.
 rnuc_vp_immers(:) = 0.

 cnuc_vc       (:) = 0.
 cnuc_vd       (:) = 0.
 cnuc_cp_hom   (:) = 0.
 cnuc_dp_hom   (:) = 0.
 cnuc_vp_haze  (:) = 0.
 cnuc_vp_immers(:) = 0.

 rpsxfer(:) = 0.
 epsxfer(:) = 0.

 r1118(:,:) = 0.;r8882(:,:) = 0.;r1112(:,:) = 0.;r1818(:,:) = 0.;r1212(:,:) = 0.
 r8282(:,:) = 0.;r3335(:,:) = 0.;r4445(:,:) = 0.;r3435(:,:) = 0.;r3445(:,:) = 0.
 r3535(:,:) = 0.;r3636(:,:) = 0.;r3737(:,:) = 0.;r4545(:,:) = 0.;r4646(:,:) = 0.
 r4747(:,:) = 0.;r5656(:,:) = 0.;r5757(:,:) = 0.;r6767(:,:) = 0.;r1413(:,:) = 0.
 r1416(:,:) = 0.;r1414(:,:) = 0.;r1446(:,:) = 0.;r1513(:,:) = 0.;r1516(:,:) = 0.
 r1515(:,:) = 0.;r1556(:,:) = 0.;r1613(:,:) = 0.;r1616(:,:) = 0.;r8483(:,:) = 0.
 r8486(:,:) = 0.;r8484(:,:) = 0.;r8446(:,:) = 0.;r8583(:,:) = 0.;r8586(:,:) = 0.
 r8585(:,:) = 0.;r8556(:,:) = 0.;r8683(:,:) = 0.;r8686(:,:) = 0.;r1713(:,:) = 0.
 r1717(:,:) = 0.;r8783(:,:) = 0.;r8787(:,:) = 0.;r2332(:,:) = 0.;r2327(:,:) = 0.
 r2323(:,:) = 0.;r2337(:,:) = 0.;r2442(:,:) = 0.;r2427(:,:) = 0.;r2424(:,:) = 0.
 r2447(:,:) = 0.;r2552(:,:) = 0.;r2527(:,:) = 0.;r2525(:,:) = 0.;r2557(:,:) = 0.
 r2662(:,:) = 0.;r2627(:,:) = 0.;r2626(:,:) = 0.;r2667(:,:) = 0.;r2772(:,:) = 0.
 r2727(:,:) = 0.;r0000(:,:) = 0.;r1812(:,:) = 0.;r1882(:,:) = 0.

 e1111(:) = 0.;e1118(:) = 0.;e8888(:) = 0.;e8882(:) = 0.;e1112(:) = 0.
 e1811(:) = 0.;e1211(:) = 0.;e8288(:) = 0.;e2222(:) = 0.;e5555(:) = 0.
 e6666(:) = 0.;e7777(:) = 0.;e3333(:) = 0.;e3335(:) = 0.;e4444(:) = 0.
 e4445(:) = 0.;e3433(:) = 0.;e3444(:) = 0.;e3435(:) = 0.;e3445(:) = 0.
 e3533(:) = 0.;e3633(:) = 0.;e3733(:) = 0.;e4544(:) = 0.;e4644(:) = 0.
 e4744(:) = 0.;e5655(:) = 0.;e5755(:) = 0.;e6766(:) = 0.;e1413(:) = 0.
 e1411(:) = 0.;e1446(:) = 0.;e1513(:) = 0.;e1511(:) = 0.;e1556(:) = 0.
 e1613(:) = 0.;e1611(:) = 0.;e8483(:) = 0.;e8488(:) = 0.;e8446(:) = 0.
 e8583(:) = 0.;e8588(:) = 0.;e8556(:) = 0.;e8683(:) = 0.;e8688(:) = 0.
 e1713(:) = 0.;e1711(:) = 0.;e8783(:) = 0.;e8788(:) = 0.;e2322(:) = 0.
 e2327(:) = 0.;e2333(:) = 0.;e2337(:) = 0.;e2422(:) = 0.;e2427(:) = 0.
 e2444(:) = 0.;e2447(:) = 0.;e2522(:) = 0.;e2527(:) = 0.;e2555(:) = 0.
 e2557(:) = 0.;e2622(:) = 0.;e2627(:) = 0.;e2666(:) = 0.;e2667(:) = 0.
 e2722(:) = 0.;e2727(:) = 0.;e2777(:) = 0.;e0000(:) = 0.;e1882(:) = 0.

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

 call mic_copy(iw0,lpw0,voa,thil0,press0,wc0,rhoa,rhow,rhoi, &
   theta0,exner0,rhov, &
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
   do while (k >= lpw0 .and. rx(k,lcat) < rxmin(lcat))
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
   press0,thil0,rhow,rhoi,exner0,tair,theta0,rhov,rhovstr,tairstrc,rx,qx,sa)

call each_column(lpw0,iw0,k1,k2,dtl0,                          &
   jhcat,press0,tair,tairc,tairstrc,rhovstr,rhoa,rhov,rhoi,    &
   rhovslair,rhovsiair,thrmcon,vapdif,dynvisc,rdynvsci,denfac, &
   colfac,colfac2,sumuy,sumuz,sumvr,                           &
   tx,sh,sm,sa,tref,rhovsref,rhovsrefp)

! USEFUL SAMPLE PRINTS FOR EXAMINING A SINGLE IW0 COLUMN BEFORE
! MAIN MICROPHYSICS COMPUTATIONS ARE DONE

!if (time_istp8 > 17300.) then
!   print*, 'micphys1 ',iw0,glatw(iw0),glonw(iw0)

!if (iw0 == 84702) then
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
!      write(6,'(a,i6,10e12.3)') 'qx ',k,(qx(k,j1),j1=1,8)
!   enddo

!   print*, ' '
!   do k = lpw0,mza0
!      write(6,'(a,i6,10e12.3)') 'con_ccnx ',k,(con_ccnx(k,j1),j1=1,nccntyp)
!   enddo

!   print*, ' '
!   do k = lpw0,mza0
!      print*, ' '
!      print*, 'micphys1 ',iw0,glatw(iw0),glonw(iw0)

!      write(6,'(a,i6,10e12.3)') 'env ',k,press0(k),theta0(k),wc0(k),tair(k), &
!                                tairc(k),rhoa(k),rhow(k),rhov(k)
!   enddo
!endif

!endif

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

jflag = 1

if (k1(1) <= k2(1)) &
   call enemb(1,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(2) <= k2(2)) &
   call enemb(2,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(3) <= k2(3)) &
   call enemb(3,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(4) <= k2(4)) &
   call enemb(4,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(5) <= k2(5)) &
   call enemb(5,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(6) <= k2(6)) &
   call enemb(6,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(7) <= k2(7)) &
   call enemb(7,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

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

 call vapdiff(iw0,k1(11),k2(11),rhov,rhovstr,sumuy,sumuz)

! Vapor flux applied to each category.  Do not change the order of these

if (k1(1) <= k2(1)) &
   call vapflux(iw0,1,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(8) <= k2(8)) &
   call vapflux(iw0,8,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(3) <= k2(3)) &
   call vapflux(iw0,3,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(4) <= k2(4)) &
   call vapflux(iw0,4,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(5) <= k2(5)) &
   call vapflux(iw0,5,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(2) <= k2(2)) &
   call vapflux(iw0,2,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(6) <= k2(6)) &
   call vapflux(iw0,6,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

if (k1(7) <= k2(7)) &
   call vapflux(iw0,7,k1,k2, &
      jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
      rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr)

! Conversion between pristine ice and snow due to vapor flux

if (k1(3) <= k2(3) .and. jnmb(4) >= 1) then
   call psxfer(iw0,k1(3),k2(3),vap,rpsxfer,epsxfer,rx,cx,qx,qr)
   k1(4) = min(k1(3),k1(4))
   k2(4) = max(k2(3),k2(4))
endif

! Diagnose new air temperature following heat and vapor fluxes

call newtemp(k1(11),k2(11), &
   tairstrc,rhoi,rhovstr,rhov,exner0,tairc,tair,theta0,rhovslair,rhovsiair,sa)

! No diagnosis or collisions if running parcel test 901 or 902

if (nl%test_case == 901 .or. nl%test_case == 902) go to 1411

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

jflag = 2

if (k1(1) <= k2(1)) &
   call enemb(1,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(2) <= k2(2)) &
   call enemb(2,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(3) <= k2(3)) &
   call enemb(3,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(4) <= k2(4)) &
   call enemb(4,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(5) <= k2(5)) &
   call enemb(5,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(6) <= k2(6)) &
   call enemb(6,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(7) <= k2(7)) &
   call enemb(7,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

! Determine coalescence efficiencies

call effxy(lpw0,k1,k2,rx,qr,emb,tx,eff)

! Liquid water collisions

j18 = max(k1(1),k1(8)); k18 = min(k2(1),k2(8))
j82 = max(k1(8),k1(2)); k82 = min(k2(8),k2(2))
j12 = max(k1(1),k1(2)); k12 = min(k2(1),k2(2))

! Self-collection of cloud droplets - transfer to drizzle

if (jnmb(8) >= 1 .and. k1(1) <= k2(1)) &
   call col1188(1,8,1,k1(1),k2(1), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
      r1118,e1111,e1118)

! Collision between cloud and drizzle - transfer to drizzle and rain

if (j18 <= k18) &
   call col1882(1,8,2,1,j18,k18, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
      r1818,r1812,r1882,e1811,e1882)

! Self-collection of drizzle - transfer to rain

if (jnmb(2) >= 1 .and. k1(8) <= k2(8)) &
   call col1188(8,2,1,k1(8),k2(8), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2, &
      r8882,e8888,e8882)

! Collisions of cloud or drizzle with rain

if (j12 <= k12) &
   call col1(1,2,2,1,j12,k12, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r1212,e1211)

if (j82 <= k82) &
   call col1(8,2,2,1,j82,k82, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r8282,e8288)

! Self collection of rain, aggregates, graupel, and hail: number change only

if (jnmb(2) == 5 .and. k1(2) <= k2(2)) &
   call cols(2,8,k1(2),k2(2), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,e2222)

if (jnmb(5) == 5 .and. k1(5) <= k2(5)) &
   call cols(5,4,k1(5),k2(5), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,e5555)

if (jnmb(6) == 5 .and. k1(6) <= k2(6)) &
   call cols(6,5,k1(6),k2(6), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,e6666)

if (jnmb(7) == 5 .and. k1(7) <= k2(7)) &
   call cols(7,7,k1(7),k2(7), &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,eff,colfac,e7777)

! Self collection of pristine ice and of snow

if (jnmb(5) >= 1 .and. k1(3) <= k2(3)) &
   call col3344(3,5,2,k1(3),k2(3),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,r3335,e3333,e3335)

if (jnmb(5) >= 1 .and. k1(4) <= k2(4)) &
   call col3344(4,5,3,k1(4),k2(4),     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac2,r4445,e4444,e4445)

! Collection between pristine ice and snow (k1,k2 now same for lcats 3 & 4)

if (jnmb(5) >= 1 .and. k1(3) <= k2(3)) &
   call col3443(2,k1(3),k2(3),iw0,     &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
      r3435,r3445,e3433,e3444,e3435,e3445)

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
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r3535,e3533)

if (j36 <= k36) &
   call col1(3,6,6,5,j36,k36, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r3636,e3633)

if (j37 <= k37) &
   call col1(3,7,7,6,j37,k37, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r3737,e3733)

if (j45 <= k45) &
   call col1(4,5,5,3,j45,k45, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r4545,e4544)

if (j46 <= k46) &
   call col1(4,6,6,5,j46,k46, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r4646,e4644)

if (j47 <= k47) &
   call col1(4,7,7,6,j47,k47, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r4747,e4744)

if (j56 <= k56) &
   call col1(5,6,6,5,j56,k56, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r5656,e5655)

if (j57 <= k57) &
   call col1(5,7,7,6,j57,k57, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r5757,e5755)

if (j67 <= k67) &
   call col1(6,7,7,5,j67,k67, &
      jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac,r6767,e6766)

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
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r1413,r1416,r1414,r1446,e1413,e1411,e1446)

   if (j15 <= k15) &
      call col2(1,5,6,1,j15,k15, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r1513,r1516,r1515,r1556,e1513,e1511,e1556)

   if (j16 <= k16) &
      call col2(1,6,6,1,j16,k16, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r1613,r1616,r0000,r0000,e1613,e1611,e0000)

   if (j84 <= k84) &
      call col2(8,4,6,1,j84,k84, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r8483,r8486,r8484,r8446,e8483,e8488,e8446)

   if (j85 <= k85) &
      call col2(8,5,6,1,j85,k85, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r8583,r8586,r8585,r8556,e8583,e8588,e8556)

   if (j86 <= k86) &
      call col2(8,6,6,1,j86,k86, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r8683,r8686,r0000,r0000,e8683,e8688,e0000)

endif

! Hail-cloud and hail-drizzle collisions (hail by-product)

if (jnmb(7) >= 1) then

   j17 = max(k1(1),k1(7)); k17 = min(k2(1),k2(7))
   j87 = max(k1(8),k1(7)); k87 = min(k2(8),k2(7))
   
   if (j17 <= k17) &
      call col2(1,7,7,1,j17,k17, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r1713,r1717,r0000,r0000,e1713,e1711,e0000)

   if (j87 <= k87) &
      call col2(8,7,7,1,j87,k87, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,emb,eff,colfac,dtl0, &
         r8783,r8787,r0000,r0000,e8783,e8788,e0000)

! Ice-rain collisions (hail byproduct)

   j23 = max(k1(2),k1(3)); k23= min(k2(2),k2(3))
   j24 = max(k1(2),k1(4)); k24= min(k2(2),k2(4))
   j25 = max(k1(2),k1(5)); k25= min(k2(2),k2(5))
   j26 = max(k1(2),k1(6)); k26= min(k2(2),k2(6))
   j27 = max(k1(2),k1(7)); k27= min(k2(2),k2(7))

   if (j23 <= k23) &
      call col3(3,7,1,j23,k23, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
         r2332,r2327,r2323,r2337,e2322,e2327,e2333,e2337)

   if (j24 <= k24) &
      call col3(4,7,1,j24,k24, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
         r2442,r2427,r2424,r2447,e2422,e2427,e2444,e2447)

   if (j25 <= k25) &
      call col3(5,7,1,j25,k25, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
         r2552,r2527,r2525,r2557,e2522,e2527,e2555,e2557)

   if (j26 <= k26) &
      call col3(6,7,1,j26,k26, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
         r2662,r2627,r2626,r2667,e2622,e2627,e2666,e2667)

   if (j27 <= k27) &
      call col3(7,7,1,j27,k27, &
         jhcat,ict1,ict2,wct1,wct2,rx,cx,qx,eff,colfac, &
         r2772,r2727,r0000,r0000,e2722,e0000,e2777,e0000)

endif

! Apply transfers of bulk mass, energy, and number from all collisions

 call colxfers(iw0,k1,k2,rx,cx,qr, &
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

1411 continue

! Nucleation of cloud droplets

if (jnmb(1) >= 1) &
   call cldnuc(iw0,lpw0,dtli0,nbincall, &
               rx,cx,qr,qx,con_ccnx,con_gccnx,rhov,rhoi,rhoa,press0, &
               tair,tairc,wc0,rhovslair,rnuc_vc,rnuc_vd,cnuc_vc,cnuc_vd)

! Rediagnose k2(1) and k3(1) because of possible new cloud nucleation

k = mza0
do while (k >= lpw0 .and. rx(k,1) < rxmin(1))
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
do while (k >= lpw0 .and. rx(k,8) < rxmin(8))
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
   call enemb(1,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

if (jnmb(8) >= 3 .and. k1(8) <= k2(8)) &
   call enemb(8,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

! Nucleation of ice crystals

if (jnmb(3) >= 1) &
   call icenuc(k1,k2,lpw0,mrl0,iw0, &
      rx,cx,qr,qx,emb,vap,tx,rhov,rhoa,press0,dynvisc,thrmcon, &
      tair,tairc,rhovslair,rhovsiair,con_ccnx,con_ifnx,dtl0, &
      rnuc_cp_hom,rnuc_dp_hom,rnuc_vp_haze,rnuc_vp_immers, &
      cnuc_cp_hom,cnuc_dp_hom,cnuc_vp_haze,cnuc_vp_immers)

! Rediagnose k2(3) and k3(3) because of possible new ice nucleation

k = mza0
do while (k >= lpw0 .and. rx(k,3) < rxmin(3))
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

! Max height for any liquid, any ice, and any of either

k2(9)  = max(k2(1),k2(2),k2(8))
k2(10) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
k2(11) = max(k2(9),k2(10))

! Do not change order of the following x02 calls

if (jnmb(3) >= 1) &
   call x02(iw0,lpw0,3,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(1) >= 1) &
   call x02(iw0,lpw0,1,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(8) >= 1) &
   call x02(iw0,lpw0,8,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(4) >= 1) &
   call x02(iw0,lpw0,4,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(5) >= 1) &
   call x02(iw0,lpw0,5,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(6) >= 1) &
   call x02(iw0,lpw0,6,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(7) >= 1) &
   call x02(iw0,lpw0,7,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

if (jnmb(2) >= 1) &
   call x02(iw0,lpw0,2,k1,k2,con_ccnx, &
      jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

  accpx(1:ncat) = 0.
  pcprx(1:ncat) = 0.

! No sedimentation, dry nuclei deposition, or precipitation scavenging of
! nuclei if running parcel tests 901-999

if (nl%test_case >= 901 .and. nl%test_case <= 999) go to 1412

! Compute sedimentation for all 7 precipitating categories

  dsed_thil(lpw0:mza0) = 0.

  call sedim2(iw0,lpw0,k1,k2,jhcat,dtl0, &
   voa, denfac, tair, thil0, theta0, dsed_thil, rhoi, rhoa, rhow, &
   cx, rx, qx, qr, emb, dmb, pcpvel, pcpfluxc, pcpfluxr, pcpfluxq, accpx, pcprx)

  ! Apply change to thil from sedim of all categories

  thil0(lpw0:mza0) = thil0(lpw0:mza0) + dsed_thil(lpw0:mza0) / rhoa(lpw0:mza0)

  if (nnuc > 0) &
     call nuclei_deposition(iw0, k1, k2, dtl0, voa, rhoa, press0, tair, &
        tairc, dynvisc, rhov, rhovslair, dmb, pcpvel, pcpfluxr, con_ccnx, &
        con_gccnx, con_ifnx)

1412 continue

! If running parcel test cases 901-999, call parcel plot to store and
! (on last timestep) plot parcel fields

if (nl%test_case >= 901 .and. nl%test_case <= 949) then
   k = 2
   kend = 2 ! Parcel simulations

   call parcel_plot(k,kend,mza0,iw0,ncat,dtli0,jhcat,rx,cx,emb,qx,tx,vap, &
      con_ccnx, &
      press0,thil0,theta0,tairc,rhovslair,rhovsiair,rhov,rhoi,rhoa,rhow, &
      rnuc_vc,rnuc_vd,rnuc_cp_hom,rnuc_dp_hom,rnuc_vp_haze,rnuc_vp_immers, &
      cnuc_vc,cnuc_vd,cnuc_cp_hom,cnuc_dp_hom,cnuc_vp_haze,cnuc_vp_immers, &
      rpsxfer,epsxfer, &
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
      e2557,e2622,e2627,e2666,e2667,e2722,e2727,e2777,e0000,e1882)

elseif (nl%test_case >= 950 .and. nl%test_case <= 999 .and. iw0 == 50) then

   if (mod(real(time_istp8),3600.) < dtl0) then

 !    kend = 18 ! Column simulation: cutoff near T = -30C (simulation/grid dependent)
      kend = 39 ! Column simulation: cutoff near T = -30C (simulation/grid dependent)

      do k = 2,kend

         call parcel_plot(k,kend,mza0,iw0,ncat,dtli0,jhcat,rx,cx,emb,qx,tx,vap, &
            con_ccnx, &
            press0,thil0,theta0,tairc,rhovslair,rhovsiair,rhov,rhoi,rhoa,rhow, &
            rnuc_vc,rnuc_vd,rnuc_cp_hom,rnuc_dp_hom,rnuc_vp_haze,rnuc_vp_immers, &
            cnuc_vc,cnuc_vd,cnuc_cp_hom,cnuc_dp_hom,cnuc_vp_haze,cnuc_vp_immers, &
            rpsxfer,epsxfer, &
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
            e2557,e2622,e2627,e2666,e2667,e2722,e2727,e2777,e0000,e1882)

      enddo

   endif

endif

! Copy hydrometeor bulk mass and number concentration and surface precipitation
! from microphysics column arrays to main model arrays

 call mic_copyback(iw0,lpw0,k1,k2,k3, &
    dtli0,accpx,pcprx,thil0,theta0,tair,rhoa,rhow,rhov, &
    con_ccnx,con_gccnx,con_ifnx,rx,cx,qx,qr,exner0)

end subroutine micphys

!===============================================================================

subroutine mic_copy(iw0,lpw0,voa,thil0,press0,wc0,rhoa,rhow,rhoi, &
   theta0,exner0,rhov, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx,qr)

use micro_coms, only: mza0, ncat, iccn, igccn, iifn, jnmb, rxmin, &
                      ccnparm, gccnparm, ifnparm, &
                      zfactor_ccn, zfactor_gccn, zfactor_ifn

use ccnbin_coms, only: nccntyp

use mem_grid,   only: zfacm2, zfacim2, arw0, arw, volt

use mem_basic,  only: thil, press, wc, rho, sh_w, sh_v, theta

use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                      q2, q6, q7, ccntyp, con_gccn, con_ifn, &
                      con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                      cldnum

use misc_coms,  only: io6, time_istp8
use consts_coms,only: r8, p00i, rocp, cpi
use therm_lib,  only: qtc

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

real, intent(inout) :: voa   (mza0)
real, intent(inout) :: thil0 (mza0)
real, intent(inout) :: press0(mza0)
real, intent(inout) :: wc0   (mza0)
real, intent(inout) :: rhoi  (mza0)
real, intent(inout) :: theta0(mza0)
real, intent(inout) :: exner0(mza0)
real, intent(inout) :: rhov  (mza0)
real, intent(inout) :: con_ccnx (mza0,nccntyp)
real, intent(inout) :: con_gccnx(mza0)
real, intent(inout) :: con_ifnx (mza0)

real(r8), intent(inout) :: rhoa(mza0)
real(r8), intent(inout) :: rhow(mza0)

real, intent(inout) :: rx(mza0,ncat)
real, intent(inout) :: cx(mza0,ncat)
real, intent(inout) :: qx(mza0,ncat)
real, intent(inout) :: qr(mza0,ncat)

integer :: k, ic

! Ratio of grid cell volume to top horizontal area arw projected onto W(k-1) level

voa(mza0) = volt(mza0,iw0) / (arw0(iw0) * zfacm2(mza0-1))
do k = lpw0,mza0-1
   voa(k) = volt(k,iw0) / (arw(k,iw0) * zfacim2(k) * zfacm2(k-1))
enddo

! copy atmospheric variables to micphys column vectors

do k = lpw0,mza0
   thil0 (k) = thil (k,iw0)
   theta0(k) = theta(k,iw0)
   press0(k) = press(k,iw0)
   wc0   (k) = wc   (k,iw0)
   rhoa  (k) = rho  (k,iw0)
   rhow  (k) = sh_w (k,iw0) * rho(k,iw0)
   rhov  (k) = sh_v (k,iw0) * rho(k,iw0)

   rhoi  (k) = 1. / rhoa(k)

   exner0(k) = (press0(k) * p00i) ** rocp  ! defined WITHOUT CP factor
enddo

! Cloud water

if (jnmb(1) == 5) then

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
         sh_c (k,iw0) = 0.0
         con_c(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(1) >= 1) then

   do k = lpw0,mza0

! If cloud water bulk density is sufficiently abundant, copy to rx.
! Otherwise, if cloud water bulk density is negative, set to zero.

      if (sh_c(k,iw0) >= rxmin(1)) then
         rx(k,1) = sh_c(k,iw0) * rhoa(k)
      else
         sh_c(k,iw0) = 0.
      endif

   enddo

endif

! Rain

if (jnmb(2) == 5) then

   do k = lpw0,mza0

! Make sure rain number concentration is positive-definite

      con_r(k,iw0) = max(con_r(k,iw0), 0.0)

! If rain bulk density is sufficiently abundant, copy to rx and copy rain
! number concentration to cx, compute qx and impose limits on it in case
! advection of very dry air has resulted in a large magnitude of q2/sh_r, and
! compute qr.  Otherwise, set rain bulk density, rain number concentration,
! and q2 to zero.

      if (sh_r(k,iw0) >= rxmin(2)) then
         rx(k,2) = sh_r (k,iw0) * rhoa(k)
         cx(k,2) = con_r(k,iw0) * rhoa(k)
         qx(k,2) = max(-20000.,min(500000.,q2(k,iw0) / sh_r(k,iw0))) ! Limits -10C to 40C
         qr(k,2) = qx(k,2) * rx(k,2)
      else
         sh_r (k,iw0) = 0.0
         con_r(k,iw0) = 0.0
         q2   (k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(2) >= 1) then

   do k = lpw0,mza0

! If rain bulk density is sufficiently abundant, copy to rx, compute qx and
! impose limits on it in case advection of very dry air has resulted in a large
! magnitude of q2/sh_r, and compute qr.  Otherwise, set rain bulk density and
! q2 to zero.

      if (sh_r(k,iw0) >= rxmin(2)) then
         rx(k,2) = sh_r(k,iw0) * rhoa(k)
         qx(k,2) = max(-20000.,min(500000.,q2(k,iw0) / sh_r(k,iw0))) ! Limits -10C to 40C
         qr(k,2) = qx(k,2) * rx(k,2)
      else
         sh_r(k,iw0) = 0.
         q2  (k,iw0) = 0.
      endif

   enddo

endif

! Pristine ice

if (jnmb(3) == 5) then

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
         sh_p (k,iw0) = 0.0
         con_p(k,iw0) = 0.0
      endif

   enddo

endif

! Snow

if (jnmb(4) == 5) then

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
         sh_s (k,iw0) = 0.0
         con_s(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(4) >= 1) then

   do k = lpw0,mza0

! If snow bulk density is sufficiently abundant, copy to rx.
! Otherwise, if snow bulk density is negative, set to zero.

      if (sh_s(k,iw0) >= rxmin(4)) then
         rx(k,4) = sh_s(k,iw0) * rhoa(k)
      else
         sh_s(k,iw0) = 0.
      endif

   enddo

endif

! Aggregates

if (jnmb(5) == 5) then

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
         sh_a (k,iw0) = 0.0
         con_a(k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(5) >= 1) then

   do k = lpw0,mza0

! If aggregates bulk density is sufficiently abundant, copy to rx.
! Otherwise, if aggregates bulk density is negative, set to zero.

      if (sh_a(k,iw0) >= rxmin(5)) then
         rx(k,5) = sh_a(k,iw0) * rhoa(k)
      else
         sh_a(k,iw0) = 0.
      endif

   enddo

endif

! Graupel

if (jnmb(6) == 5) then

   do k = lpw0,mza0

! Make sure graupel number concentration is positive-definite

      con_g(k,iw0) = max(con_g(k,iw0), 0.0)

! If graupel bulk density is sufficiently abundant, copy to rx and copy graupel
! number concentration to cx, compute qx and impose limits on it in case
! advection of very dry air has resulted in a large magnitude of q6/sh_g, and
! compute qr.  Otherwise, set graupel bulk density, graupel number concentration,
! and q6 to zero.

      if (sh_g(k,iw0) >= rxmin(6)) then
         rx(k,6) = sh_g (k,iw0) * rhoa(k)
         cx(k,6) = con_g(k,iw0) * rhoa(k)
         qx(k,6) = max(-100000.,min(334000.,q6(k,iw0) / sh_g(k,iw0))) ! Limits -50C to 0C
         qr(k,6) = qx(k,6) * rx(k,6)
      else
         sh_g (k,iw0) = 0.0
         con_g(k,iw0) = 0.0
         q6   (k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(6) >= 1) then

   do k = lpw0,mza0

! If graupel bulk density is sufficiently abundant, copy to rx, compute qx and
! impose limits on it in case advection of very dry air has resulted in a large
! magnitude of q6/sh_g, and compute qr.  Otherwise, set graupel bulk density and
! q6 to zero.

      if (sh_g(k,iw0) >= rxmin(6)) then
         rx(k,6) = sh_g(k,iw0) * rhoa(k)
         qx(k,6) = max(-100000.,min(334000.,q6(k,iw0) / sh_g(k,iw0))) ! Limits -50C to 0C
         qr(k,6) = qx(k,6) * rx(k,6)
      else
         sh_g(k,iw0) = 0.
         q6  (k,iw0) = 0.
      endif

   enddo

endif

! Hail

if (jnmb(7) == 5) then

   do k = lpw0,mza0

! Make sure hail number concentration is positive-definite

      con_h(k,iw0) = max(con_h(k,iw0), 0.0)

! If hail bulk density is sufficiently abundant, copy to rx and copy hail
! number concentration to cx, compute qx and impose limits on it in case
! advection of very dry air has resulted in a large magnitude of q7/sh_h, and
! compute qr.  Otherwise, set hail bulk density, hail number concentration,
! and q7 to zero.

      if (sh_h(k,iw0) >= rxmin(7)) then
         rx(k,7) = sh_h (k,iw0) * rhoa(k)
         cx(k,7) = con_h(k,iw0) * rhoa(k)
         qx(k,7) = max(-100000.,min(334000.,q7(k,iw0) / sh_h(k,iw0))) ! Limits -50C to 0C
         qr(k,7) = qx(k,7) * rx(k,7)
      else
         sh_h (k,iw0) = 0.0
         con_h(k,iw0) = 0.0
         q7   (k,iw0) = 0.0
      endif

   enddo

elseif (jnmb(7) >= 1) then

   do k = lpw0,mza0

! If hail bulk density is sufficiently abundant, copy to rx, compute qx and
! impose limits on it in case advection of very dry air has resulted in a large
! magnitude of q7/sh_h, and compute qr.  Otherwise, set hail bulk density and
! q7 to zero.

      if (sh_h(k,iw0) >= rxmin(7)) then
         rx(k,7) = sh_h(k,iw0) * rhoa(k)
         qx(k,7) = max(-100000.,min(334000.,q7(k,iw0) / sh_h(k,iw0))) ! Limits -50C to 0C
         qr(k,7) = qx(k,7) * rx(k,7)
      else
         sh_h(k,iw0) = 0.
         q7  (k,iw0) = 0.
      endif

   enddo

endif

! Drizzle

if (jnmb(8) == 5) then

   do k = lpw0,mza0

! Make sure drizzle number concentration is positive-definite

      con_d(k,iw0) = max(con_d(k,iw0), 0.0)

! If drizzle bulk density is sufficiently abundant, copy to rx
! and copy drizzle number concentration to cx.
! Otherwise, make sure drizzle bulk density is positive-definite
! and set drizzle number concentration to zero.

      if (sh_d(k,iw0) >= rxmin(8)) then
         rx(k,8) = sh_d (k,iw0) * rhoa(k)
         cx(k,8) = con_d(k,iw0) * rhoa(k)
      else
         sh_d (k,iw0) = 0.0
         con_d(k,iw0) = 0.0
      endif

   enddo

endif

! Fill column CCN values [#/m^3] 

if (iccn == 1) then
   if (ccnparm > 1.e6) then
      do k = lpw0,mza0
         con_ccnx(k,1) = ccnparm * zfactor_ccn(k) * rhoa(k)
      enddo
   else
      do k = lpw0,mza0
         con_ccnx(k,1) = cldnum(iw0) * zfactor_ccn(k) * rhoa(k)
      enddo
   endif
else
   do ic = 1,nccntyp
      do k = lpw0,mza0
         ccntyp(ic)%con_ccn(k,iw0) = max(0.,ccntyp(ic)%con_ccn(k,iw0))
         con_ccnx(k,ic) = ccntyp(ic)%con_ccn(k,iw0) * rhoa(k)
      enddo
   enddo
endif

! Fill column GCCN values [#/m^3] 

if (igccn == 1) then
   do k = lpw0,mza0
      con_gccnx(k) = gccnparm * zfactor_gccn(k) * rhoa(k)
   enddo
else
   do k = lpw0,mza0
      con_gccn(k,iw0) = max(0.,con_gccn(k,iw0))
      con_gccnx(k) = con_gccn(k,iw0) * rhoa(k)
   enddo
endif

! Fill column IFN values [#/m^3] 

if (iifn == 1) then
   do k = lpw0,mza0
      con_ifnx(k) = ifnparm * zfactor_ifn(k) * rhoa(k)
   enddo
else
   do k = lpw0,mza0
      con_ifn(k,iw0) = max(0.,con_ifn(k,iw0))
      con_ifnx(k) = con_ifn(k,iw0) * rhoa(k)
   enddo
endif

end subroutine mic_copy

!===============================================================================

subroutine mic_copyback(iw0,lpw0,k1,k2,k3, &
   dtli0,accpx,pcprx,thil0,theta0,tair0,rhoa,rhow,rhov, &
   con_ccnx,con_gccnx,con_ifnx,rx,cx,qx,qr,exner0)

use micro_coms, only: mza0, ncat, jnmb, iccn, igccn, iifn, rxmin
use ccnbin_coms, only: nccntyp
use mem_basic,  only: thil, theta, tair, rho, sh_w, sh_v, wmc, wc
use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                      q2, q6, q7, &
                      con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                      accpd, accpr, accpp, accps, accpa, accpg, accph, &
                      pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                      ccntyp, con_gccn, con_ifn
use misc_coms,  only: io6
use consts_coms,only: r8, p00i, rocp, alvl, alvi, cpi4, cp253i
use mem_grid,   only: zwgt_bot, zwgt_top
use mem_tend,   only: num_omic
use var_tables, only: num_scalar, scalar_tab
use therm_lib,  only: qtc

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)
integer, intent(in) :: k3(11)

real, intent(in) :: dtli0

real, intent(in) :: accpx(ncat)
real, intent(in) :: pcprx(ncat)

real, intent(in)  :: thil0(mza0)
real, intent(in)  :: theta0(mza0)
real, intent(in)  :: tair0(mza0)
real, intent(in)  :: rhov(mza0)

real(r8), intent(in) :: rhoa(mza0)
real(r8), intent(in) :: rhow(mza0)

real, intent(in) :: con_ccnx (mza0,nccntyp)
real, intent(in) :: con_gccnx(mza0)
real, intent(in) :: con_ifnx (mza0)

real, intent(in) :: rx(mza0,ncat)
real, intent(in) :: cx(mza0,ncat)
real, intent(in) :: qx(mza0,ncat)
real, intent(in) :: qr(mza0,ncat)

real, intent(in) :: exner0(mza0)

real     :: rfact(mza0)
real(r8) :: rhoi(mza0)
integer  :: k, n, ic
real     :: til, qhydm, tc, fracliq

real :: rhoice(mza0), rholiq(mza0)

! Copy base thermodynamic variables

do k = lpw0,mza0
   thil(k,iw0)  = thil0(k)
   theta(k,iw0) = theta0(k)
   tair(k,iw0)  = tair0(k)
   rhoi(k)      = 1.0_r8 / rhoa(k)
   rfact(k)     = rhoi(k) * rho(k,iw0)
   rho(k,iw0)   = rhoa(k)
   sh_w(k,iw0)  = rhow(k) * rhoi(k)
   sh_v(k,iw0)  = min( real(rhov(k) * rhoi(k)), sh_w(k,iw0))
   rholiq(k)    = 0.
   rhoice(k)    = 0.
enddo

! Adjust scalar specific densities to conserve mass when air density changes

do n = num_omic+1, num_scalar
   !dir$ ivdep
   do k = lpw0, mza0
      scalar_tab(n)%var_p(k,iw0) = scalar_tab(n)%var_p(k,iw0) * rfact(k)
   enddo
enddo

! Conserve vertical momentum when air density changes

do k = lpw0, mza0-1
   wmc(k,iw0) = wc(k,iw0) * (zwgt_bot(k) * rho(k,iw0) + zwgt_top(k) * rho(k+1,iw0))
enddo

! Copy cloud water bulk density back to main array

if (jnmb(1) >= 1) then
   do k = lpw0,k3(1)
      sh_c(k,iw0) = rx(k,1) * rhoi(k)
      if (rx(k,1) > rxmin(1)) rholiq(k) = rholiq(k) + rx(k,1)
   enddo
endif

! Copy rain bulk density, internal energy, and surface precip back to main arrays

if (jnmb(2) >= 1) then
   accpr(iw0) = accpr(iw0) + real(accpx(2),r8)
   pcprr(iw0) = pcprx(2)

   do k = lpw0,k2(11)
      sh_r(k,iw0) = rx(k,2) * rhoi(k)
      q2  (k,iw0) = qr(k,2) * rhoi(k)
      if (rx(k,2) > rxmin(2)) rholiq(k) = rholiq(k) + rx(k,2)
   enddo
endif

! Copy pristine ice bulk density and surface precip back to main arrays

if (jnmb(3) >= 1) then
   accpp(iw0) = accpp(iw0) + real(accpx(3),r8)
   pcprp(iw0) = pcprx(3)

   do k = lpw0,k3(3)
      sh_p(k,iw0) = rx(k,3) * rhoi(k)
      if (rx(k,3) > rxmin(3)) rhoice(k) = rhoice(k) + rx(k,3)
   enddo
endif

! Copy snow bulk density and surface precip back to main arrays

if (jnmb(4) >= 1) then
   accps(iw0) = accps(iw0) + real(accpx(4),r8)
   pcprs(iw0) = pcprx(4)

   do k = lpw0,k2(11)
      sh_s(k,iw0) = rx(k,4) * rhoi(k)
      if (rx(k,4) > rxmin(4)) rhoice(k) = rhoice(k) + rx(k,4)
   enddo
endif

! Copy aggregates bulk density and surface precip back to main arrays

if (jnmb(5) >= 1) then
   accpa(iw0) = accpa(iw0) + real(accpx(5),r8)
   pcpra(iw0) = pcprx(5)

   do k = lpw0,k2(11)
      sh_a(k,iw0) = rx(k,5) * rhoi(k)
      if (rx(k,5) > rxmin(5)) rhoice(k) = rhoice(k) + rx(k,5)
   enddo
endif

! Copy graupel bulk density, internal energy, and surface precip back to main arrays

if (jnmb(6) >= 1) then
   accpg(iw0) = accpg(iw0) + real(accpx(6),r8)
   pcprg(iw0) = pcprx(6)
   do k = lpw0,k2(11)
      sh_g(k,iw0) = rx(k,6) * rhoi(k)
      q6  (k,iw0) = qr(k,6) * rhoi(k)
      if (rx(k,6) > rxmin(6)) then
         call qtc(qx(k,6),tc,fracliq)
         rholiq(k) = rholiq(k) + rx(k,6) * fracliq
         rhoice(k) = rhoice(k) + rx(k,6) * (1. - fracliq)
      endif
   enddo
endif

! Copy hail bulk density, internal energy, and surface precip back to main arrays

if (jnmb(7) >= 1) then
   accph(iw0) = accph(iw0) + real(accpx(7),r8)
   pcprh(iw0) = pcprx(7)
   do k = lpw0,k2(11)
      sh_h(k,iw0) = rx(k,7) * rhoi(k)
      q7  (k,iw0) = qr(k,7) * rhoi(k)
      if (rx(k,7) > rxmin(7)) then
         call qtc(qx(k,7),tc,fracliq)
         rholiq(k) = rholiq(k) + rx(k,7) * fracliq
         rhoice(k) = rhoice(k) + rx(k,7) * (1. - fracliq)
      endif
   enddo
endif

! Copy drizzle bulk density and surface precip back to main arrays

if (jnmb(8) >= 1) then
   accpd(iw0) = accpd(iw0) + real(accpx(8),r8)
   pcprd(iw0) = pcprx(8)

   do k = lpw0,k3(8)
      sh_d(k,iw0) = rx(k,8) * rhoi(k)
      if (rx(k,8) > rxmin(8)) rholiq(k) = rholiq(k) + rx(k,8)
   enddo
endif

! Copy cloud water number concentration back to main array

if (jnmb(1) == 5) then
   do k = lpw0,k3(1)
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

if (jnmb(3) == 5) then
   do k = lpw0,k3(3)
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

if (jnmb(8) == 5) then
   do k = lpw0,k3(8)
      con_d(k,iw0) = cx(k,8) * rhoi(k)
   enddo
endif

! Copy CCN number concentration back to main array

if (iccn >= 2) then
   do ic = 1,nccntyp
      do k = lpw0,mza0
         ccntyp(ic)%con_ccn(k,iw0) = con_ccnx(k,ic) * rhoi(k)
      enddo
   enddo
endif

! Copy GCCN number concentration back to main array

if (igccn == 2) then
   do k = lpw0,mza0
      con_gccn(k,iw0) = con_gccnx(k) * rhoi(k)
   enddo
endif

! Copy IFN number concentration back to main array

if (iifn == 2) then
   do k = lpw0,mza0
      con_ifn(k,iw0) = con_ifnx(k) * rhoi(k)
   enddo
endif

! Rediagnose tair and theta

do k = lpw0, max( k3(1), k3(3), k3(8), k2(11) )

   til = thil0(k) * exner0(k)
   qhydm = alvl * rholiq(k) + alvi * rhoice(k)

   if (tair0(k) > 253.) then
      tair(k,iw0) = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm * rhoi(k))))
   else
      tair(k,iw0) = til * (1. + qhydm * rhoi(k) * cp253i)
   endif

   theta(k,iw0) = tair(k,iw0) / exner0(k)

enddo

end subroutine mic_copyback
