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
subroutine cuparm_kuo(iw)

use mem_tend,    only: thilt, sh_wt
use mem_cuparm,  only: thsrc, rtsrc, conprr
use mem_basic,   only: uc, vc, wc, theta, press, rho, sh_v
use mem_grid,    only: mza, mwa, lpu, lpv, lpw, zm, zt, volt,  &
                       unx, uny, unz, vnx, vny, vnz
use misc_coms,   only: io6, initial, time_istp8, confrq, dtlong, itime1, &
                       meshtype, wcldbs
use mem_ijtabs,  only: itab_w, istp
use kuo_coms,    only: icprtfl, icpltfl, wconmin, contim, igo, cprecip
use consts_coms, only: p00i, rocp, cp, grav, alvl

implicit none

integer, intent(in) :: iw

! automatic arrays

real, dimension(mza) :: zzcon, zcon, uecon, vecon, wecon, wcon, thtcon, hzcon
real, dimension(mza) :: prcon, dncon, rvcon, picon, tmpcon, ftcon, frcon

integer :: k,iv,ka,k2,kc,kv,npoly,j
integer :: icpcnt=0,iprtfrq,iwqmax,kqmax
real    :: dthmax,piocp,polyi

!        FLAG TO CONTROL PRINTOUT
!          ICPRTFL=0 - NO PRINTOUT
!                  1 - BRIEF INFO ON UP/DOWN DRAFT AND ITERATIONS
!                  2 - 1 PLUS MODEL TENDENCIES
!                  3 - 2 PLUS FINAL CONVECTIVE STRUCTURE
!                  4 - 3 PLUS UP/DOWN DRAFT AND ENVIRONMENT

!-------------------------------------------------------------------------

! Set flags/counter for print/plot options

icprtfl = 0
iprtfrq = 8
icpltfl = 0
icpcnt = icpcnt + 1

if (mod(icpcnt-iprtfrq+1,iprtfrq) == 0) then
   icprtfl = 1
endif

! Do modified Kuo convective cumulus parameterization

ka = lpw(iw)

! Compute one more than the number of prognostic levels above ground 
! for this IW column

k2 = mza + 1 - ka 

! Zero out surface convective precipitation rate 

conprr(iw) = 0.

! Zero out 3 EARTH wind components

uecon(1:mza) = 0.
vecon(1:mza) = 0.
wecon(1:mza) = 0.

! Number of polygon edges for current IW

npoly = itab_w(iw)%npoly
polyi = 1. / real(npoly)

! Vertical loop over prognostic T levels

do k = ka,mza-1

! Zero out 3D cumulus parameterization heating and moistening rates

   thsrc(k,iw) = 0.
   rtsrc(k,iw) = 0.

! K index for above-ground convective column of model levels

   kc = k - ka + 2 

! Compute heights AGL using zm at bottom of lpw cell as ground level

   zzcon(kc) = zm(k) - zm(ka-1)
   zcon(kc)  = zt(k) - zm(ka-1)

   wcon(kc)   = wc(k,iw)
   thtcon(kc) = theta(k,iw)
   prcon(kc)  = press(k,iw)
   dncon(kc)  = rho(k,iw)
   rvcon(kc)  = sh_v(k,iw)

   piocp      = (prcon(kc) * p00i)**rocp
   picon(kc)  = cp * piocp
   tmpcon(kc) = thtcon(kc) * piocp
   hzcon(kc)  = cp * tmpcon(kc) + grav * zcon(kc) + alvl * rvcon(kc)

! Loop over U/V neighbors of W
! Avoid underground levels for IV point   

   do j = 1,npoly
   
      if (meshtype == 1) then
         iv = itab_w(iw)%iu(j)
         kv = max(k,lpu(iv))
      else
         iv = itab_w(iw)%iv(j)
         kv = max(k,lpv(iv))
      endif

! Sum over V neighbors the 3 EARTH components of horizontal wind   
   
      uecon(kc) = uecon(kc) + uc(kv,iv) * unx(iv) + vc(kv,iv) * vnx(iv)
      vecon(kc) = vecon(kc) + uc(kv,iv) * uny(iv) + vc(kv,iv) * vny(iv)
      wecon(kc) = wecon(kc) + uc(kv,iv) * unz(iv) + vc(kv,iv) * vnz(iv)
   enddo

! Get average EARTH wind components

   uecon(kc) = uecon(kc) * polyi
   vecon(kc) = vecon(kc) * polyi
   wecon(kc) = wecon(kc) * polyi

enddo

! Fill values for kc = 1 in case any get used

zzcon(1) = 0. 
zcon(1)  = 0.

uecon(1)  = uecon(2)
vecon(1)  = vecon(2)
wecon(1)  = wecon(2)

wcon(1)   = wcon(2)
thtcon(1) = thtcon(2)
prcon(1)  = prcon(2)
dncon(1)  = dncon(2)
rvcon(1)  = rvcon(2) 

picon(1)  = picon(2)
tmpcon(1) = tmpcon(2)
hzcon(1)  = hzcon(2)

wconmin = wcldbs
contim = confrq

call cu_environ(k2, zcon, zzcon, tmpcon, rvcon, uecon, vecon, wecon, &
                wcon, dncon, thtcon, picon, hzcon)

if (igo /= 0) call kuocp(iw)

if (igo /= 0) then

   call cp2mod(k2,ftcon,frcon,zzcon,dncon,picon)

! Bob: multiply ftcon and frcon by rho(k,iw) so that 
!      thsrc units are [kg_air K / (m^3 s)] and
!      rtsrc units are [kg_wat / (m^3 s)].

   do k = ka,mza-1
      kc = k - ka + 2 
      thsrc(k,iw) = ftcon(kc) * rho(k,iw)
      rtsrc(k,iw) = frcon(kc) * rho(k,iw)
   enddo

   conprr(iw) = cprecip

   if (icprtfl > 0) then
      write(io6, '(A,I0,A,F0.3)') ' CONVECTION AT IW=',iw,' TIME=',time_istp8
   endif

endif

return
end subroutine cuparm_kuo

!===============================================================================

subroutine cu_environ(k2, zcon, zzcon, tmpcon, rvcon, uecon, vecon, wecon, &
                      wcon, dncon, thtcon, picon, hzcon)

use kuo_coms,    only: dzlow, dzhigh, zmid, cdzmin, igo, wconmin, zc, &
                       ze, kmt, uepe, vepe, wepe, wpe, the, rve, pke, &
                       thve, te, pe, rhoe, thee, kcon, tlcl, plcl,    &
                       dzlcl, klcl, theu, ketl, nkp

use consts_coms, only: cp, cpi, grav, grav2, alvl, alvlocp, cpor, p00, &
                       rdry, eps_virt
use misc_coms,   only: io6

implicit none

integer, intent(in) :: k2
real,    intent(in) :: zcon(k2),  zzcon(k2),  tmpcon(k2), rvcon(k2)
real,    intent(in) :: uecon(k2), vecon(k2),  wecon(k2),  wcon(k2)
real,    intent(in) :: dncon(k2), thtcon(k2), picon(k2),  hzcon(k2)

real :: wcpmax,themax,tlll,plll,rlll,zlll,dzlll,dzdd,abe,thdu,tdu,rdsu,znz
integer :: k,nkmid

!   Basic constants

dzlow  = 200.
dzhigh = 500.
zmid   = 3000.
cdzmin = 3000.

!   Check for conditional instability and any upward motion
!      greater than WCONMIN under ZMID

igo = 0
do k = 1,k2-1
   if (hzcon(k) > hzcon(k+1)) then
      igo = 1
      exit
   endif
enddo

if (igo == 0) return

igo = 0
wcpmax = -1.e10

do k = 1,k2
   if (zcon(k) > zmid) exit
   wcpmax = max(wcpmax,wcon(k))
enddo

if (wcpmax > 0. .and. wcpmax > wconmin) igo = 1
if (igo == 0) return

!   INTERPOLATE MODEL SOUNDING (ENVIRONMENT) TO HIGHER RESOLUTION GRID

nkmid = zmid / dzlow + 1
zc(1) = 0.

do k = 2,nkmid
   zc(k) = zc(k-1) + dzlow
enddo

do k = nkmid+1,nkp
   zc(k) = zc(k-1) + dzhigh
enddo

ze(1) = 0.
do k = 2,nkp
   ze(k) = (zc(k) + zc(k-1)) * .5
enddo

!   FIND MODEL TOP ON CONVECTIVE GRID

znz = zcon(k2)

do k = nkp,1,-1
   if (ze(k) < znz) go to 13
enddo
stop ' envir stop 12'
13 continue
kmt = k

!   DO ACTUAL INTERPOLATION

call hintrp_cc(k2, uecon, zcon,kmt,uepe,ze)
call hintrp_cc(k2, vecon, zcon,kmt,vepe,ze)
call hintrp_cc(k2, wecon, zcon,kmt,wepe,ze)
call hintrp_cc(k2,  wcon,zzcon,kmt,wpe, ze)
call hintrp_cc(k2,thtcon, zcon,kmt,the, ze)
call hintrp_cc(k2, rvcon, zcon,kmt,rve, ze)

do k = 1,kmt
   rve(k) = max(rve(k),1.e-8)
enddo

!   COMPUTE THETA V, THETA E, AND GET PRESSURE PROFILE

do k = 1,kmt
   thve(k) = the(k) * (1. + eps_virt * rve(k))
enddo

pke(1) = picon(1)
do k = 2,kmt
   pke(k) = pke(k-1) - grav2 * (ze(k) - ze(k-1)) / (thve(k) + thve(k-1))
enddo

do k = 1,kmt
   te(k) = the(k) * pke(k) * cpi
   pe(k) = (pke(k) * cpi)**cpor * p00
   rhoe(k) = pe(k) / (rdry * te(k) * (1. + eps_virt * rve(k)))
enddo

do k = 1,kmt
   call thetae(pe(k),te(k),rve(k),thee(k))
enddo

!   FIND THE MAIN SOURCE LEVEL OF THE UPDRAFT.

!   FIRST TEST - ANY INVERSION BELOW 1.2 KM

do k = 3,nkmid
   if (te(k) > te(k-1) .and. te(k) > te(k+1) .and. ze(k) <= 1200.) then
      kcon = k
      go to 77
   endif
enddo

!   IF THERE ISN'T AN INVERSION, USE THE LEVEL OF HIGHEST THETA E .

themax = 0.
do k = 2,nkmid
   if (thee(k) > themax) then
      themax = thee(k)
      kcon = k
   endif
enddo

!   FIND THE LCL OF A LAYER AVERAGE AROUND THE SOURCE LEVEL

77 continue

tlll = (te(kcon) + te(kcon+1) + te(kcon-1)) / 3.
plll = pe(kcon)
rlll = (rve(kcon) + rve(kcon+1) + rve(kcon-1)) / 3.
zlll = ze(kcon)

call lcl(tlll,plll,rlll,tlcl,plcl,dzlcl)

!   FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL

dzlll = 1.e20
do k = 1,kmt
   dzdd = abs(ze(k) - (zlll + dzlcl))
   if (dzdd < dzlll) then
      dzlll = dzdd
      klcl = k
   endif
enddo
klcl = max(klcl,3)

!   IF THERE IS NOT UPWARD MOTION AT THE LCL, NO CONVECTION
!     (MUST BE GREATER THAN WCONMIN )

if (wpe(klcl) < 0. .or. wpe(klcl) < wconmin) then
   igo = 0
   return
endif

!   LOCATE EQUILIBRIUM TEMPERATURE LEVEL OF AN UNENTRAINED PARCEL.
!   COMPUTE INITIAL ABE.  IF ABE IS LESS THAN 0, NO CONVECTION.

theu(klcl) = the(kcon) * exp(alvlocp * rve(kcon) / tlcl)

do k = klcl,kmt
   if (theu(klcl) <= thve(k)) go to 66
enddo

write(io6,*) 'convection above model top:',klcl,theu(klcl),thve(k)

ketl = kmt - 2
!cccccc      stop 65

66 continue

ketl = k

if (ze(ketl) - ze(klcl) < cdzmin) then
   igo = 0
   return
endif

abe = 0.
do k = klcl,ketl
   call the2t(theu(klcl),pe(k),thdu,tdu,rdsu)
   abe = abe + (thdu * (1. + eps_virt * rdsu) - thve(k)) / thve(k)  &
       * (zc(k) - zc(k-1))
enddo

if (abe <= 0.) then
   igo = 0
   return
endif

!     if(icprtfl > 0)then
!       write(io6, 899)
! 899   format(///,' * convection is activated * ')
!     endif

return
end subroutine cu_environ

!===============================================================================

subroutine kuocp(iw)

use kuo_coms,    only: nkp, ftcone, frcone, rhoe, rve, wpe, supply, igo, theu, &
                       the, tlcl, klcl, kmt, pe, thu, tu, rsu, klfc, ketl,  &
                       kct, ze, cdzmin, zc, thee, thd, kcon, wtd, envshr, &
                       uepe, vepe, wepe, preff, vheat, vmois, vmdry, pke,  &
                       qvct1, contim, cprecip, icprtfl, thcone
use consts_coms, only: alvl, cp
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw

integer :: k,idownd,klfs,kdiv,kdet,kover,kcoolh,kheat
real :: supplyw,anegl,apos,anegh,dddt,dzdiv,wtlfs,wtlcl,wtdiv,wtgnd &
       ,bkuo,zdetr,dzdet,vhint,vmint,vdint,avgmin,avtdiff,overmax,factr &
       ,heatmx,coolhi,c1 

!     Downdraft flag - 0 - no downdrafts
!                      1 - simple downdraft model

idownd = 1

do k = 1,nkp
   ftcone(k) = 0.
   frcone(k) = 0.
enddo

!     Compute vertical moisture convergence into the cloud layer.
!       Vertical flux out cloud top is assumed small.

supplyw = rhoe(klcl) * rve(klcl) * (wpe(klcl) + wpe(klcl-1)) * .5

supply = supplyw

!if (iw == 100) write(io6,*) 'supply1 ',supply,wpe(klcl),rve(klcl),rhoe(klcl)


if (supply <= 0.) then
   igo = 0
   return
endif

!     This is the cloud model.  Updraft is constant THETA e and
!       saturated with respect to water.  There is no ice.
!       Cloud top is one level above ETL.
!
!     THETA e of the updraft

theu(klcl) = the(kcon) * exp(alvl * rve(kcon) / (cp * tlcl))

!     Equilibrium Temperature Level of the source level air.

igo = 0
do k = klcl,kmt
   call the2t(theu(klcl),pe(k),thu(k),tu(k),rsu(k))
   if (thu(k) > the(k) .and. igo == 0) then
      igo = 1
      klfc = k
   endif
   if (thu(k) <= the(k) .and. igo == 1) go to 66
enddo

if (igo == 0) return

write(io6,*) ' Convection beyond model top - THup, THenv ',THU(KMT),THE(KMT)

k = kmt - 1

66 continue

ketl = min(k,kmt)
kct = min(ketl+1,kmt)

call the2t(theu(klcl),pe(kct),thu(kct),tu(kct),rsu(kct))

do k = 1,klcl-1
   thu(k) = the(k)
enddo

!     If the cloud is not at least CDZMIN deep or cloud top is
!       under 500 mb, no convection.

if (ze(ketl) - ze(klfc) < cdzmin .or. pe(kct) > 50000.) then
   igo = 0
   return
endif

!     Require the positive area be 50% greater than the negative
!       area below the LFC and  5% greater in total.

anegl = 0.

do k = klcl,klfc-1
   anegl = anegl + (thu(k) - the(k)) * (zc(k) - zc(k-1))
enddo

apos = 0.

do k = klfc,ketl-1
   apos = apos + (thu(k) - the(k)) * (zc(k) - zc(k-1))
enddo

anegh = 0.

do k = ketl,kct
   anegh = anegh + (thu(k) - the(k)) * (zc(k) - zc(k-1))
enddo

if (apos < abs(anegl) * 1.5 .or. apos < abs(anegl + anegh) * 1.05) then
   igo = 0
   return
endif

if (idownd == 1) then

!     The downdraft model - starts at THETA e minimum (LFS).
!         Downdraft is 2 degrees colder than
!         environment at cloud base increasing to 5 degrees
!         colder at the ground.


!         Find LFS as THETA e minimum

   do k = kct,2,-1
      if (thee(k) < thee(k+1) .and. thee(k) < thee(k-1)) go to 11
   enddo

   k = 2
   11 continue
   klfs = k

   if (klfs <= klcl) klfs = klcl + 1
   thd(klfs) = the(klfs)

!    Limit dd deficit at the ground to the maximum of positive
!      temperature difference of updraft if less than 2.5 degrees.

   dddt = 0.

   do k = klcl,kct
      dddt = max(dddt,thu(k)-the(k))
   enddo

   if (dddt > 2.5) dddt = 5.

   thd(2) = the(2) - dddt
   thd(klcl) = the(klcl) - dddt * .2

   do k = klcl,klfs
      thd(k) = thd(klcl) + (thd(klfs) - thd(klcl)) / (ze(klfs) - ze(klcl))  &
             * (ze(k) - ze(klcl))
   enddo

   do k = 3,klcl-1
      thd(k) = thd(2) + (thd(klcl) - thd(2)) / (ze(klcl) - ze(2))  &
             * (ze(k) - ze(2))
   enddo

!     Now we need to weight the downdraft relative to the updraft.
!       Assume that the dd weight is zero at the LFS, 1/2 of
!       updraft at cloud base, and equal to the updraft at cloud
!       base at the ground.

   dzdiv = 1e20
   do k = 1,kmt
      if (abs(ze(k) - 800.) < dzdiv) then
         kdiv = k
         dzdiv = abs(ze(k) - 800.)
      endif
   enddo

   kdiv = max(min(klcl,kdiv),2)
   if (kdiv == klcl) kdiv = klcl - 1

   do k = 1,nkp
      wtd(k) = 0.
   enddo

   wtlfs = 0.
   wtlcl = .1
   wtdiv = .2
   wtgnd = 1.

   do k = klcl+1,klfs
      wtd(k) = wtlcl + (wtlfs - wtlcl) / (ze(klfs) - ze(klcl))  &
             * (ze(k) - ze(klcl))
   enddo

   do k = kdiv,klcl
      wtd(k) = wtdiv + (wtlcl - wtdiv) / (ze(klcl) - ze(kdiv))  &
             * (ze(k) - ze(kdiv))
   enddo

   do k = 2,kdiv-1
      wtd(k) = wtgnd + (wtdiv - wtgnd) / (ze(kdiv) - ze(2))  &
             * (ze(k) - ze(2))
   enddo

else

   do k = 1,nkp
      wtd(k) = 0.
   enddo

   do k = 2,klcl-1
      thu(k) = the(k)
   enddo 

endif

!     Compute infamous b parameter.  Use Fritsch/Chappell's
!       precipitation efficiency.

envshr = sqrt((uepe(kct) - uepe(klfc))**2  &
       +      (vepe(kct) - vepe(klfc))**2  &
       +      (wepe(kct) - wepe(klfc))**2) &
       / (ze(kct) - ze(klfc)) * 1.e3

if (envshr > 1.35) then
   preff = 1.591 - .639 * envshr + .0953 * envshr**2 - .00496 * envshr**3
else
   preff = .9
endif

bkuo = 1. - preff

!     Vertical profiles of convective heating and moistening

do k = 2,kmt
   vheat(k) = 0.
   vmois(k) = 0.
   vmdry(k) = 0.
enddo

!     Find the weighted THETA to use for the convection.

do k = 2,kct
   thcone(k) = wtd(k) * thd(k) + (1. - wtd(k)) * thu(k)
enddo

!     Heating profile is difference between convective THETAs and environment.

do k = 2,kct
   vheat(k) = thcone(k) - the(k)
enddo

!     Moisture profile is difference between vapor's of updraft and
!       environment in the cloud layer.  Below cloud base, air is
!       dried by SUPPLY.  Downdrafts are assumed to have no effect
!       on this.

zdetr = .66667 * ze(kct)
dzdet = 1000000.

do k = klcl,kct
   if (abs(ze(k) - zdetr) < dzdet) then
      dzdet = abs(ze(k) - zdetr)
      kdet = k
   endif
enddo

do k = kdet,kct
   vmois(k) = 1.
enddo

!   do k=klcl,kct
!     vmois(k)=rsu(k)-rve(k)
!   enddo

do k = 2,klcl-1
   vmdry(k) = rve(k)
enddo

vhint = 0.
vmint = 0.
vdint = 0.
do k = 2,kmt
   vhint = vhint + vheat(k) * (zc(k) - zc(k-1))
   vmint = vmint + vmois(k) * (zc(k) - zc(k-1))
   vdint = vdint + vmdry(k) * (zc(k) - zc(k-1))
enddo

!     If VHINT is less than 0, there is more negative area than
!       positive area.  No convection allowed.

if (vhint <= 0.) then
   igo = 0
   return
endif

!     Also require that there is a minimum average
!       temperature difference between the updraft and environment
!       from the LFC to the ETL.  This eliminates the cases where
!       VHINT is very small and the heating and cooling rates get
!       astronomically large.

avgmin = .10
avtdiff = 0.

do k = klfc,ketl-1
   avtdiff = avtdiff + (thcone(k) - the(k))
enddo

avtdiff = avtdiff / max(1,ketl-klfc)
if (avtdiff < avgmin) then
   igo = 0
   return
endif

!     Heating and moistening rates

3100 continue

do k = 2,kmt
   ftcone(k) = alvl * preff * supply * vheat(k) / (pke(k) * rhoe(k) * vhint)
enddo

do k = klcl,kct
   frcone(k) = bkuo * supply * vmois(k) / (rhoe(k) * vmint)
enddo

do k = 2,klcl-1
   frcone(k) = -supply * vmdry(k) / (rhoe(k) * vdint)
enddo

do k = klfc,ketl-1
   qvct1(k) = the(k) + contim * ftcone(k)
enddo

overmax = 0.

do k = klfc,ketl-1
   if (qvct1(k) - thu(k) > overmax) then
      overmax = (qvct1(k) - thu(k)) / (ftcone(k) * contim)
      kover = k
   endif
enddo

if (overmax > 0.) then
   factr = 1. - overmax
   supply = factr * supply
!bob        if(icprtfl.ge.1)write(io6,*) ' reducing supply ',kover,factr
!bob     +   ,qvct1(kover),thu(kover)
   go to 3100
endif

cprecip = preff * supply

!if (iw == 100) write(io6,*) 'cprecip1 ',cprecip,preff,supply


if (icprtfl > 0) then
!bob        write(io6,*)' ----------------------------------------------------'
!bob        write(io6,898) ZE(KCON),ZE(KLCL),ZE(KLFC),ZE(KETL),ZE(KCT)
!bob  898   FORMAT(' CLOUD LEVELS - SOURCE,LCL,LFC,ETL,TOP(KM) ',-3P,5F5.1)
!bobC       write(io6,896) SUPPLY/SUPPLYW
! 896   FORMAT(' SUPPLIES ' ,F8.4)
!bob        write(io6,897) THEU(KLCL),RSU(KLCL)*1E3
897   FORMAT(' CLOUD PROPERTIES -  THETA E, RS AT LCL',F6.1,F8.2)
   coolhi = 100000.
   heatmx = -10000.
   kcoolh = 0
   kheat = 0
   do k = ketl,kct
      if (ftcone(k) < coolhi) then
         coolhi = ftcone(k)
         kcoolh = k
      endif
   enddo
   do k = klcl,kct
      if (ftcone(k) > heatmx) then
         heatmx = ftcone(k)
         kheat = k
      endif
   enddo

   C1 = 86400.
!bob        write(io6,905) PE(KHEAT),HEATMX*C1,PE(KCOOLH),COOLHI*C1
905   FORMAT(' MAX-MIN HEATING- P,(K/DAY)', 2(-2PF8.1,0PF7.1))
!bob        write(io6,906) PREFF,PREFF*SUPPLY*3600.*10.
!bob     +           ,(WPE(KLCL)+WPE(KLCL-1))*.5
!bob     +           ,RVE(KLCL)*1E3
906   FORMAT(' PRECIPITATION - EFFICIENCY,RATE(MM/HR)',F5.2,F6.2,  &
  '  LCL W,RV',F8.4,F7.2)
ENDIF
!
IF (ICPRTFL >= 2) THEN
!bob      write(io6,95) (K,ZE(K),PE(K),TE(K),THE(K),THEE(K),RVE(K),UCON(K)
!bob     +   ,VCON(K),WPE(K),K=1,KMT)
95 FORMAT(//' ENVIRONMENT-K,Z,P,TE,THE,THEE,RVE,UP,VP,WP'/,  &
   (I3, -2P,F8.1,-3P,F8.2,0P,3F7.2,3P,F6.2, -2P,3F8.2))
!bob      write(io6,96) (K,THE(K),THU(K)   ,FTCON(K)*86400.,RVE(K),RSU(K),
!bob     +       FRCON(K)*86400.,VHEAT(K)*(ZC(K)-ZC(K-1))/VHINT,
!bob     +                       VMOIS(K)*(ZC(K)-ZC(K-1))/VMINT,
!bob     +                       VMDRY(K)*(ZC(K)-ZC(K-1))/VDINT,
!bob     +  K=2,KMT)
96 FORMAT(//' HEATING/MOISTENING'  &
 ,'-K,THE,THU,THSRC,RVE,RSU,RTSRC,HEAT%,MOIST%,DRY%'/,  &
   (I3,0P,3F7.1,3P,3F7.1,0P,3F6.2))
ENDIF
!
!      IF(ICPLTFL.EQ.1)THEN
!      FTCON(1)=FTCON(2)
!      FRCON(1)=FRCON(2)
!      FTMAX=-2000.
!      FRMAX=-2000.
!      FTMIN=2000.
!      FRMIN=2000.
!      DO 400 K=1,KMT
!      QVCT1(K)=FTCON(K)*86400.
!      QVCT2(K)=FRCON(K)*86400.
!      FTMAX=MAX(FTMAX,QVCT1(K))
!      FTMIN=MIN(FTMIN,QVCT1(K))
!      FRMAX=MAX(FRMAX,QVCT2(K))
!      FRMIN=MIN(FRMIN,QVCT2(K))
!      QVCT3(K)=ZE(K)*1E-5
!  400 CONTINUE
!
!        CALL DISPLA(2,1,1)
!        CALL SET (.05,.55,.15,.9,FTMIN-20.,FTMAX+20.
!     +            ,QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'dTH/dt$','height (km)$',1,4,1,'$$$$')
!c       CALL LINE(0.,0.,0.,QVCT3(KMT))
!       CALL EZXY(QVCT1,QVCT3,KMT,'Heating rates (K/day) $')
!C
!        CALL DISPLA(2,1,1)
!        CALL SET (.55,.99,.15,.9,FRMIN-.1,FRMAX+.1
!     +            ,QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'drT/dt$',' $',1,4,0,0)
!        CALL LINE(0.,0.,0.,QVCT3(KMT))
!        CALL EZXY(QVCT2,QVCT3,KMT,'Moistening rates (g/g/day)$')
!        CALL FRAME
!
!C
!      DO 410 K=2,KMT
!      THGRID=THE(K)+FTCON(K)*CONTIM
!      QVCT2(K)=RVE(K)+FRCON(K)*CONTIM
!      QVCT4(K)=RVE(K)
!      TGRID=THGRID*(PE(K)/P00)**ROCP
!      RVVV=RS(PE(K),TGRID)
!      QVCT1(K)=THGRID+AKLV*MAX(0.,QVCT2(K)-RVVV)
!  410 CONTINUE
!      QVCT1(1)=QVCT1(2)
!      QVCT2(1)=QVCT2(2)
!      QVCT4(1)=QVCT4(2)
!        CALL DISPLA(2,1,1)
!        CALL SET (.05,.55,.15,.8,290.,QVCT1(KMT),QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'THETA$','height (km)$',1,4,1,'$$$$')
!        CALL EZXY(THE,QVCT3,KMT,'INITIAL-SOLID  AFTER-DASHED$')
!        CALL ANOTAT(CHAR(0),CHAR(0),4,4,1,'$''$''')
!        CALL EZXY(QVCT1,QVCT3,KMT,CHAR(0))
!
!        CALL DISPLA(2,1,1)
!        CALL SET (.55,.97,.15,.8,  0.,QVCT4(2),QVCT3(1),QVCT3(KMT),1)
!        CALL ANOTAT( 'MIXING RATIO$',' $',1,4,1,'$$$$')
!        CALL EZXY(QVCT4,QVCT3,KMT,'INITIAL-SOLID  AFTER-DASHED$')
!        CALL ANOTAT(CHAR(0),' $',4,4,1,'$''$''')
!        CALL EZXY(QVCT2,QVCT3,KMT,CHAR(0))
!        CALL FRAME
!
!      ENDIF
!-----------------------------------------------------------------------
return
end subroutine kuocp

!===============================================================================

subroutine cp2mod(k2,ftcon,frcon,zzcon,dncon,picon)

use kuo_coms,    only: kmt, qvct1, rhoe, ftcone, pke, qvct2, frcone, qvct3, &
                       zc, qvct4
use consts_coms, only: alvl
use misc_coms,   only: io6

implicit none

integer, intent(in)  :: k2
real,    intent(in)  :: zzcon(k2), dncon(k2), picon(k2)
real,    intent(out) :: ftcon(k2), frcon(k2)

integer :: k
real    :: tftc,tftm,tfrc,tfrm,ftres,frres

real :: vctr5(k2)  ! automatic array
real :: vctr6(k2)  ! automatic array

!    Compute integrated heating and moistening tendencies

do k = 2,kmt
   qvct1(k) = rhoe(k) * ftcone(k) * pke(k)
   qvct2(k) = rhoe(k) * alvl * frcone(k)
   qvct3(k) = (zc(k) - zc(k-1)) * qvct1(k)
   qvct4(k) = (zc(k) - zc(k-1)) * qvct2(k)
enddo

tftc = sum(qvct3(2:kmt))
tfrc = sum(qvct4(2:kmt))

!     Transfer tendencies to model grid

call vertmap3(qvct1,zc,kmt,vctr5,zzcon,k2)
call vertmap3(qvct2,zc,kmt,vctr6,zzcon,k2)

do k = 2,k2
   vctr5(k) = vctr5(k) * (zzcon(k) - zzcon(k-1))
   vctr6(k) = vctr6(k) * (zzcon(k) - zzcon(k-1))
enddo

!     Make sure the transfer from the convective grid to the model
!       grid happened correctly.

tftm = sum(vctr5(2:k2))
tfrm = sum(vctr6(2:k2))

ftres = tftm - tftc
frres = tfrm - tfrc

if (abs(ftres) > .01 * abs(tftc)) then
   write(io6,*) ' energy error in grid tranfser in convective param.'
   write(io6,*) ' tftm,tftc ',tftm,tftc
endif

!     Change energy tendencies to temperature and mixing ratio tendencies.

do k = 2,k2
   ftcon(k) = vctr5(k) / ((zzcon(k) - zzcon(k-1)) * dncon(k) * picon(k))
   frcon(k) = vctr6(k) / ((zzcon(k) - zzcon(k-1)) * dncon(k) * alvl)
enddo

return
end subroutine cp2mod

!===============================================================================

subroutine vertmap3(datin,zin,n3in,datout,zout,n3out)

! Given two arrays of vertical height values, ZIN and ZOUT, and data values
! DATIN defined in layers between ZIN heights, subroutine vertmap3 interpolates
! DATIN values to DATOUT values defined in layers between ZOUT heights over
! the range of heights common to both.  The vertical integral is preserved 
! over this range.  

implicit none

integer, intent(in) :: n3in
integer, intent(in) :: n3out

real, intent(in)  :: datin(n3in)
real, intent(in)  :: zin(n3in)
real, intent(in)  :: zout(n3out)
real, intent(out) :: datout(n3out)

integer :: kin,kout
real :: overlap

! Set datout values to zero prior to summation

do kout = 1,n3out
   datout(kout) = 0.
enddo

! Begin with level 2 of both datin and datout

kin = 2
kout = 2

do while (kin <= n3in .and. kout <= n3out) 
   overlap = max(0.,min(zin(kin),zout(kout))-max(zin(kin-1),zout(kout-1)))
   datout(kout) = datout(kout) + overlap * datin(kin)
   if (zin(kin) < zout(kout)) then
      kin = kin + 1
   else
      kout = kout + 1
   endif
enddo

! Normalize datout values by thickness of datout layer

do kout = 2,n3out
   datout(kout) = datout(kout) / (zout(kout) - zout(kout-1))
enddo

return
end subroutine vertmap3

!===============================================================================

subroutine lcl(t0,pp0,r0,tlcl,plcl,dzlcl)

use consts_coms, only: p00, p00k, p00ki, cp, cpi, grav, rocp, cpor, cpog, &
                       eps_virt
use misc_coms,  only: io6

implicit none

real, intent(in) :: t0,pp0,r0
real, intent(out) :: tlcl,plcl,dzlcl

integer :: nitt,ip
real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs,weight
real, external :: td, rslf

ip = 0
11 continue

if(ip == 0)then
   weight = 1.0
else
   weight = 0.3
endif

plcl = pp0
tlcl = t0
p0k = pp0**rocp
pi0i = p0k * p00ki * cp
ttth0= t0 * p00k / p0k
ttd = td(pp0,r0)
dz = weight * cpog * (t0 - ttd)

if (dz <= 0.) then
   dzlcl = 0.
   return
endif

do nitt = 1,50
   pki = pi0i - grav * dz / (ttth0 * (1. + eps_virt * r0))
   pppi = (pki * cpi)**cpor * p00
   ti = ttth0 * pki * cpi
   rvs = rslf(pppi,ti)
   if (abs(rvs - r0) < .00003) go to 110
   ttd = td(pppi,r0)
   dz = dz + weight * cpog * (ti - ttd)
enddo

write(io6,*) 'no converge in LCL:',t0,pp0,r0
ip = ip + 1
if (ip == 1) go to 11
stop 'LCL no convergence'

110 continue
plcl = pppi
tlcl = ti
dzlcl = dz

return
end subroutine lcl
