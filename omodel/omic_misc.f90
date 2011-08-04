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

!  FOR ICNFLG=3, DEBATING WHETHER TO KEEP IT #/M4 OR CHANGE
!  PARM TO #/KG/M.  NEED TO DEFINE AVMIPSA, ETC. FOR ALL CATEGORIES.
!  MAY WANT TO DEFINE C1 TOO AND RENAME IT.
!  IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections
!  need to consider collection efficiencies for different habits?
!  collection efficiency for hail too high.  big hail should not
!  coallesce.

subroutine each_column(lpw0,iw0,k1,k2,dtl0,                    &
   jhcat,press0,tair,tairc,tairstrc,rhovstr,rhoa,rhov,rhoi,    &
   rhovslair,rhovsiair,thrmcon,vapdif,dynvisc,rdynvsci,denfac, &
   colfac,colfac2,sumuy,sumuz,sumvr,                           &
   tx,sh,sm,sa,tref,rhovsref,rhovsrefp)

use micro_coms,  only: mza0, ncat, jhabtab
use consts_coms, only: alvl, alvi
use misc_coms,   only: io6

implicit none

integer, intent(in) :: lpw0
integer, intent(in) :: iw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in) :: dtl0

integer, intent(out) :: jhcat(mza0,ncat)

real, intent(in)  :: press0   (mza0)
real, intent(in ) :: tair     (mza0)
real, intent(out) :: tairc    (mza0)
real, intent(in)  :: tairstrc (mza0)
real, intent(in)  :: rhovstr  (mza0)
real, intent(in)  :: rhov     (mza0)
real, intent(in)  :: rhoi     (mza0)
real, intent(out) :: rhovslair(mza0)
real, intent(out) :: rhovsiair(mza0)
real, intent(out) :: thrmcon  (mza0)
real, intent(out) :: vapdif   (mza0)
real, intent(out) :: dynvisc  (mza0)
real, intent(out) :: rdynvsci (mza0)
real, intent(out) :: denfac   (mza0)
real, intent(out) :: colfac   (mza0)
real, intent(out) :: colfac2  (mza0)
real, intent(out) :: sumuy    (mza0)
real, intent(out) :: sumuz    (mza0)
real, intent(out) :: sumvr    (mza0)

real(kind=8), intent(in) :: rhoa(mza0)

real, intent(out) :: tx(mza0,ncat)
real, intent(out) :: sh(mza0,ncat)
real, intent(out) :: sm(mza0,ncat)

real, intent(inout) :: sa(mza0,9)

real, intent(out) :: tref     (mza0,2)
real, intent(out) :: rhovsref (mza0,2)
real, intent(out) :: rhovsrefp(mza0,2)

integer :: k,nt,ns
real :: ck1,ck2,ck3,rhovsref1,rhovsref2,relhum,colf
real, external :: rhovsl,rhovsi

data ck1,ck2,ck3/-4.818544e-3,1.407892e-4,-1.249986e-7/

colf  = .785 * dtl0

! Loop over all atmospheric levels for this column

do k = lpw0,mza0
   tairc(k) = tair(k) - 273.15
   tx(k,1) = tairc(k)
   thrmcon(k) = ck1 + (ck2 + ck3 * tair(k)) * tair(k)
   dynvisc(k) = .1718e-4 + .49e-7 * tairc(k)
   denfac(k) = sqrt(rhoi(k))

   rhovslair(k) = rhovsl(tairc(k))
   rhovsiair(k) = rhovsi(tairc(k))

! Diagnose habit of pristine ice and snow

   nt = max(1,min(31,-nint(tairc(k))))
   relhum = min(1.,rhov(k) / rhovslair(k))
   ns = max(1,nint(100. * relhum))

   jhcat(k,1) = 1
   jhcat(k,2) = 2
   jhcat(k,3) = jhabtab(nt,ns,1)
   jhcat(k,4) = jhabtab(nt,ns,2)
   jhcat(k,5) = 5
   jhcat(k,6) = 6
   jhcat(k,7) = 7
   jhcat(k,8) = 8
   
enddo

! Loop over the range of levels with pre-existing condensate

do k = k1(11),k2(11)
   vapdif(k) = 2.14 * (tair(k) / 273.15) ** 1.94 / press0(k)
   rdynvsci(k) = sqrt(1. / dynvisc(k))

   colfac(k)  = colf * denfac(k)
   colfac2(k) = 2. * colfac(k)

   tref(k,1) = tairc(k) - min(25.,700. * (rhovslair(k) - rhov(k)) * rhoi(k))

   sa(k,2) = thrmcon(k) * sa(k,1)  ! stays the same
   sa(k,3) = thrmcon(k) * (tairstrc(k) * rhoa(k) + sa(k,1) * rhovstr(k))  
    ! sa3 picks up factor of rhoa (rhoa had to be inserted in first term)

   sumuy(k) = 0.
   sumuz(k) = 0.
   sumvr(k) = 0.
enddo

! Loop over the range of levels with pre-existing liquid

do k = k1(9),k2(9)

! Compute rhovsrefp by centered finite difference over range of 1 K

   rhovsref(k,1)  = rhovsl(tref(k,1)     )
   rhovsref2      = rhovsl(tref(k,1) + .5)
   rhovsref1      = rhovsl(tref(k,1) - .5)
   rhovsrefp(k,1) = rhovsref2 - rhovsref1

   sa(k,4) = rhovsrefp(k,1) * tref(k,1) - rhovsref(k,1)  ! x rhoa factor
   sa(k,6) = alvl * rhovsrefp(k,1)                       ! x rhoa factor
   sa(k,8) = alvl * sa(k,4)                              ! x rhoa factor

   sh(k,1) = 0.
   sh(k,2) = 1.
   sh(k,8) = 0.

   sm(k,1) = 1.
   sm(k,2) = 1.
   sm(k,8) = 1.
enddo

! Loop over the range of levels with pre-existing ice

do k = k1(10),k2(10)
   tref(k,2)    = min(0.,tref(k,1))

! Compute rhovsrefp by onc-sided finite difference over range of 1 K

   rhovsref(k,2)  = rhovsi(tref(k,2)     )
   rhovsref1      = rhovsi(tref(k,2) - 1.)
   rhovsrefp(k,2) = rhovsref(k,2) - rhovsref1

   sa(k,5) = rhovsrefp(k,2) * tref(k,2) - rhovsref(k,2)  ! x rhoa factor
   sa(k,7) = alvi * rhovsrefp(k,2)                       ! x rhoa factor
   sa(k,9) = alvi * sa(k,5)                              ! x rhoa factor

   sh(k,3) = 0.
   sh(k,4) = 0.
   sh(k,5) = 0.
   sh(k,6) = 1.
   sh(k,7) = 1.
enddo

return
end subroutine each_column

!===============================================================================

subroutine enemb(lcat,jflag,k1,k2, &
   ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

use micro_coms, only: mza0, ncat, jnmb, emb2, cfemb0, pwemb0, cfen0, pwen0, &
                      parm, emb0, emb1, rxmin, enmlttab, dict, emb0log, &
                      rictmin, rictmax
use misc_coms,  only: io6

implicit none

integer, intent(in) :: lcat
integer, intent(in) :: jflag

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

integer, intent(out) :: ict1(mza0,ncat)
integer, intent(out) :: ict2(mza0,ncat)

real, intent(out) :: wct1(mza0,ncat)
real, intent(out) :: wct2(mza0,ncat)

real, intent(in)    :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(out)   :: emb (mza0,ncat)
real, intent(in)    :: vap (mza0,ncat)
real, intent(in)    :: rhoi(mza0)

real(kind=8), intent(in) :: rhoa(mza0)

integer :: k,lhcat
real :: embi,parmi,fracmass,cxloss
real :: rict,rictmm

if (jnmb(lcat) == 2) then

   do k = k1(lcat),k2(lcat)
      emb(k,lcat) = emb2(lcat)
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo

elseif (jnmb(lcat) == 4) then

   parmi = 1. / parm(lcat)
   do k = k1(lcat),k2(lcat)
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo

elseif (jnmb(lcat) >= 5) then

   do k = k1(lcat),k2(lcat)
      emb(k,lcat) = &
         max(emb0(lcat),min(emb1(lcat),rx(k,lcat) / max(1.e-12,cx(k,lcat))))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo

endif

if (jflag == 2) then
   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) < rxmin(lcat)) cycle

      rict = dict(lcat) * (log(emb(k,lcat)) - emb0log(lcat)) + 1.
      rictmm = max(rictmin,min(rictmax,rict))
      ict1(k,lcat) = int(rictmm)
      ict2(k,lcat) = ict1(k,lcat) + 1
      wct2(k,lcat) = rictmm - float(ict1(k,lcat))
      wct1(k,lcat) = 1.0 - wct2(k,lcat)

   enddo
endif

return
end subroutine enemb

!===============================================================================

subroutine x02(iw0,lcat,k1,k2, &
   jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qx,qr,vap,rhoa,rhoi)

use micro_coms,  only: mza0, ncat, rxmin, enmlttab, dnfac, pwmasi, gnu, shedtab
use consts_coms, only: alli
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lcat

integer, intent(inout) :: k1(11)
integer, intent(inout) :: k2(11)

integer, intent(in) :: jhcat(mza0,ncat)

integer, intent(out) :: ict1(mza0,ncat)
integer, intent(out) :: ict2(mza0,ncat)

real, intent(out) :: wct1(mza0,ncat)
real, intent(out) :: wct2(mza0,ncat)

real, intent(inout) :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(inout) :: qx  (mza0,ncat)
real, intent(inout) :: qr  (mza0,ncat)
real, intent(inout) :: emb (mza0,ncat)
real, intent(in)    :: vap (mza0,ncat)
real, intent(in)    :: rhoi(mza0)

real(kind=8), intent(in) :: rhoa(mza0)

integer :: k,jflag,lhcat,inc,idns
real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed, &
        qrmltshed,fracmloss,dn

real, parameter :: shedmass = 5.236e-7

! Now, any category may have mass anywhere in the range k1(11),k2(11).  
! Rediagnose k1(lcat) and k2(lcat).

! Find new k2(lcat)

k = k2(11)
do while (k > k1(11) .and. rx(k,lcat) < rxmin(lcat))
   k = k - 1
enddo
k2(lcat) = k

! Find new k1(lcat)

k = k1(11)
do while (k <= k2(lcat) .and. rx(k,lcat) < rxmin(lcat))
   k = k + 1
enddo
k1(lcat) = k

! Return if there is no significant bulk mass for lcat category in this column

if (k1(lcat) > k2(lcat)) return

! Diagnose bulk mean mass and/or number concentration for lcat

if (lcat /= 1 .and. lcat /= 3 .and. lcat /= 8) then

   jflag = 1
   call enemb(lcat,jflag,k1,k2  &
      ,ict1,ict2,wct1,wct2,rx,cx,emb,vap,rhoa,rhoi)

endif

! PRISTINE ICE category

if (lcat == 3) then

! Steve: Allow pristine ice to melt to cloud1 since we assume that smaller
! particles will melt first. Perhaps need a way in the future to treat
! melting of larger pristine ice into cloud2.

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rinv

         call qtc(qx(k,lcat),tcoal,fracliq)

! If PRISTINE ICE is melting, transfer liquid to cloud category

         if (fracliq > 1.e-6) then

            rmelt = rx(k,lcat) * fracliq
            cmelt = cx(k,lcat) * fracliq

            rx(k,lcat) = rx(k,lcat) - rmelt
            cx(k,lcat) = cx(k,lcat) - cmelt
            qx(k,lcat) = 0. ! (not needed)
            qr(k,lcat) = 0.

            rx(k,1) = rx(k,1) + rmelt
            cx(k,1) = cx(k,1) + cmelt

         endif

      endif

   enddo

! SNOW and AGGREGATES categories

elseif (lcat == 4 .or. lcat == 5) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rinv

         call qtc(qx(k,lcat),tcoal,fracliq)

! If SNOW or AGGREGATES are melting, transfer liquid plus some ice 
! to GRAUPEL category

         if (fracliq > 1.e-6) then

            rmelt = rx(k,lcat) * fracliq
            ricetor6 = min(rmelt, rx(k,lcat) - rmelt)

! change this??? move to rain instead ??? look at melting decisions in col2

            rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
            qx(k,lcat) = 0. ! (not needed)
            qr(k,lcat) = 0.
            
            rx(k,6) = rx(k,6) + rmelt + ricetor6
            qr(k,6) = qr(k,6) + rmelt * alli

! keep the above the same with ricetor6
! meyers - use sa melt table here? yes
!
            fracmloss = (rmelt + ricetor6) * rinv
            closs = enmlttab(int(200. * fracmloss) + 1,jhcat(k,lcat)) * cx(k,lcat)

            cx(k,lcat) = cx(k,lcat) - closs
            cx(k,6) = cx(k,6) + closs

         endif

      endif

   enddo

! GRAUPEL category

elseif (lcat == 6) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rxinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rxinv

         call qtc(qx(k,lcat),tcoal,fracliq)

! If GRAUPEL is more than 95% melted, transfer all to RAIN category

         if (fracliq > 0.95) then

            rx(k,2) = rx(k,2) + rx(k,6)
            qr(k,2) = qr(k,2) + rx(k,6) * alli
            cx(k,2) = cx(k,2) + cx(k,6)

            rx(k,6) = 0.
            qr(k,6) = 0.
            cx(k,6) = 0.
            qx(k,6) = 0. ! (not needed)

         endif

      endif

   enddo

! HAIL category

elseif (lcat == 7) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rxinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rxinv

!c          qx(k,lcat) = max(-50.,qx(k,lcat))

         call qtc(qx(k,lcat),tcoal,fracliq)

! If HAIL is more than 95% melted, transfer all to RAIN category

         if (fracliq > 0.95) then

            rx(k,2) = rx(k,2) + rx(k,7)
            qr(k,2) = qr(k,2) + rx(k,7) * alli
            cx(k,2) = cx(k,2) + cx(k,7)

            rx(k,7) = 0.
            qr(k,7) = 0.
            cx(k,7) = 0.
            qx(k,7) = 0. ! (not needed)
         
! Otherwise, if HAIL is more than 30% liquid, compute amount to be shed to rain

!  take out following IF statement?

         elseif (fracliq > 0.3) then

            lhcat = jhcat(k,lcat)
            inc = nint(200. * fracliq) + 1
            dn = dnfac(lhcat) * emb(k,lcat) ** pwmasi(lhcat)
            idns = max(1,nint(1.e3 * dn * gnu(lcat)))
            rshed = rx(k,lcat) * shedtab(inc,idns)
!cc               rmltshed = rx(k,lcat) * rmlttab(inc) + rshed
            rmltshed = rshed
            qrmltshed = rmltshed * alli

            rx(k,2) = rx(k,2) + rmltshed
            qr(k,2) = qr(k,2) + qrmltshed

            rx(k,lcat) = rx(k,lcat) - rmltshed
            qr(k,lcat) = qr(k,lcat) - qrmltshed
!               closs = cx(k,lcat) * enmlttab(inc,lhcat)
!               cx(k,lcat) = cx(k,lcat) - closs
!               cx(k,2) = cx(k,2) + closs + rshed / shedmass
            cx(k,2) = cx(k,2) + rshed / shedmass

! The following part not needed

            if (rx(k,7) > rxmin(7)) then
               qx(k,7) = qr(k,7) / rx(k,7)
            else
               qx(k,7) = 0.
            endif
         endif

      endif

   enddo

! RAIN category

elseif (lcat == 2) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rxinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rxinv

! limit rain temperature to range of about (-40C,40C)

         qx(k,lcat) = max(.5 * alli,min(1.5 * alli,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rx(k,lcat)

      endif

   enddo

endif

return
end subroutine x02

!===============================================================================

subroutine sedim(lcat,iw0,lpw0,k1,k2,alphasfc, &
   dtli0,accpx,pcprx,pcpg0,qpcpg0,dpcpg0,jhcat, &
   rx,cx,qx,qr,emb,thil0,theta0,tair,denfac,rhoi,rhoa,rhow,dsed_thil,ch1)

use micro_coms,  only: mza0, ncat, rxmin, cfmasi, ch3, ch2, dispemb0i, &
                       nembfall, pcpfillc, pcpfillr, dztf, nhcat, maxkfall, &
                       emb1
use consts_coms, only: cpi
use misc_coms,   only: io6

implicit none

integer, intent(in) :: lcat
integer, intent(in) :: iw0
integer, intent(in) :: lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in)    :: alphasfc
real, intent(in)    :: dtli0
real, intent(out)   :: accpx(ncat)
real, intent(out)   :: pcprx(ncat)
real, intent(inout) :: pcpg0
real, intent(inout) :: qpcpg0
real, intent(inout) :: dpcpg0

integer, intent(in) :: jhcat(mza0,ncat)

real, intent(inout) :: rx (mza0,ncat)
real, intent(inout) :: cx (mza0,ncat)
real, intent(inout) :: qx (mza0,ncat)
real, intent(in)    :: qr (mza0,ncat)
real, intent(in)    :: emb(mza0,ncat)

real, intent(in)    :: thil0    (mza0)
real, intent(in)    :: theta0   (mza0)
real, intent(in)    :: tair     (mza0)
real, intent(in)    :: denfac   (mza0)
real, intent(in)    :: rhoi     (mza0)
real, intent(inout) :: dsed_thil(mza0)

real(kind=8), intent(inout) :: rhoa(mza0)
real(kind=8), intent(inout) :: rhow(mza0)

real, intent(in) :: ch1(nhcat)

real, parameter :: iplaws = 0

integer :: k,lhcat,iemb,iemb2,kkf,kk
real :: dispemb,riemb,wt2,rsfc,qrsfc

! automatic arrays

real :: rfall (2-maxkfall:mza0)
real :: cfall (2-maxkfall:mza0)
real :: qrfall(2-maxkfall:mza0)

! Zero out any "fall" cells that might accumulate precipitation

cfall (lpw0+1-maxkfall:k2(lcat)) = 0.
rfall (lpw0+1-maxkfall:k2(lcat)) = 0.
qrfall(lpw0+1-maxkfall:k2(lcat)) = 0.

! Loop over potential donor cells

do k = k1(lcat),k2(lcat)
   lhcat = jhcat(k,lcat)

! Jump to end of loop if current cell has no hydrometeor mass

   if (rx(k,lcat) < rxmin(lcat)) cycle

! Compute displacement over one timestep and sedimentation table index

! Here determine which set of powerlaws to use: the original
!  ones in RAMS or the ones from R.Carver adapted from Mitchell 1996.
! The Mitchell power laws are not based at sea level so we adjust the
!  density factor based at 0.7 kg/m3 instead of 1.0 kg/m3.

   if (iplaws == 0) then

      dispemb = ch1(lhcat)  &
         * (emb(k,lcat) * cfmasi(lhcat)) ** ch3(lhcat) * denfac(k)

   else

      dispemb = ch1(lhcat)  &
         * (emb(k,lcat) * cfmasi(lhcat)) ** ch3(lhcat) * (0.7 * rhoi(k))**.362

   endif

   riemb = 1. + ch2(lhcat) * log10(dispemb * dispemb0i(lhcat))

!Bob (10/24/00):  Now, limiting iemb to max of nembfall

   iemb = min(nint(riemb),nembfall)

! Loop over receptor cells

   do kkf = 1,maxkfall
      kk = k + 1 - kkf   ! receptor ("fall") cell index

! Accumulate density of hydrometeor number, mass, and energy in receptor cell kk

      cfall (kk) = cfall (kk) + cx(k,lcat) * pcpfillc(k,kkf,iemb,lhcat)
      rfall (kk) = rfall (kk) + rx(k,lcat) * pcpfillr(k,kkf,iemb,lhcat)
      qrfall(kk) = qrfall(kk) + qr(k,lcat) * pcpfillr(k,kkf,iemb,lhcat)
   
   enddo

enddo

! Copy accumulated precip in "below-ground" cells to surface precip

rsfc = 0.
qrsfc = 0.

do k = lpw0+1-maxkfall,lpw0-1

   rsfc  = rsfc  + rfall(k)  * dztf(k)
   qrsfc = qrsfc + qrfall(k) * dztf(k)
   
enddo

accpx(lcat) = rsfc
pcprx(lcat) = rsfc * dtli0

pcpg0  = pcpg0  + rsfc
qpcpg0 = qpcpg0 + qrsfc
dpcpg0 = dpcpg0 + rsfc * alphasfc

! Compute change in rhoa, rhow, and thil

do k = lpw0,k2(lcat)

   rhoa(k) = rhoa(k) + rfall(k) - rx(k,lcat)
   rhow(k) = rhow(k) + rfall(k) - rx(k,lcat)
   
   dsed_thil(k) = dsed_thil(k) - thil0(k) * thil0(k)  &
      * (2820. * (rfall(k) - rx(k,lcat))  &
      - cpi * (qrfall(k) - qr(k,lcat)))  &
      / (max(tair(k), 253.) * theta0(k))

! Transfer "fall" amounts to category arrays

   rx(k,lcat) = rfall(k)
   cx(k,lcat) = cfall(k)
   qx(k,lcat) = qrfall(k) / max(rxmin(lcat),rfall(k))

   if (rx(k,lcat) < rxmin(lcat)) then
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
   endif

enddo

return
end subroutine sedim

