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

!  collection efficiency for hail too high.  big hail should not
!  coallesce.

subroutine each_column(lpw0,iw0,k1,k2,dtl0,                    &
   jhcat,press0,tair,tairc,tairstrc,rhovstr,rhoa,rhov,rhoi,    &
   rhovslair,rhovsiair,thrmcon,vapdif,dynvisc,rdynvsci,denfac, &
   colfac,colfac2,sumuy,sumuz,sumvr,                           &
   tx,sh,sm,sa,tref,rhovsref,rhovsrefp)

use micro_coms,  only: mza0, ncat, jhabtab
use consts_coms, only: r8, alvl, alvi
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

real(r8), intent(in) :: rhoa(mza0)

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
   dynvisc(k) = .1718e-4 + .49e-7 * tairc(k) ! Units are [kg/(m^2 s)]
   denfac(k) = sqrt(rhoi(k))
   colfac(k)  = colf * denfac(k)
   colfac2(k) = 2. * colfac(k)

   rhovslair(k) = rhovsl(tairc(k))
   rhovsiair(k) = rhovsi(tairc(k))

   ! New arrays for dry/wet deposition and sedim2

   

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

end subroutine each_column

!===============================================================================

subroutine enemb(lcat,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

use micro_coms,  only: mza0, ncat, jnmb, emb2, emb0, emb1, rxmin, &
                       enmlttab, dict, emb0log, rictmin, rictmax
use ccnbin_coms, only: nccntyp
use misc_coms,   only: io6
use consts_coms, only: r8

implicit none

integer, intent(in) :: lcat
integer, intent(in) :: jflag

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in) :: con_ccnx(mza0,nccntyp)

integer, intent(out) :: ict1(mza0,ncat)
integer, intent(out) :: ict2(mza0,ncat)

real, intent(out) :: wct1(mza0,ncat)
real, intent(out) :: wct2(mza0,ncat)

real, intent(in)    :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(out)   :: emb (mza0,ncat)
real, intent(in)    :: vap (mza0,ncat)

integer :: k,lhcat
real :: rict,rictmm

if (jnmb(lcat) == 2) then

   do k = k1(lcat),k2(lcat)
      emb(k,lcat) = emb2(lcat)
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo

elseif (jnmb(lcat) == 4) then ! As of version 5.0.0, can only apply to cloud

   do k = k1(lcat), k2(lcat)
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat), rx(k,lcat) / max(1.e-12,con_ccnx(k,1))))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo

elseif (jnmb(lcat) == 5) then

   do k = k1(lcat), k2(lcat)
      emb(k,lcat) = max( emb0(lcat),                                            &
                         min( emb1(lcat), rx(k,lcat) / max(1.e-12,cx(k,lcat)) ) )
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

end subroutine enemb

!===============================================================================

subroutine x02(iw0,lpw0,lcat,k1,k2,con_ccnx, &
   jhcat,ict1,ict2,wct1,wct2,rx,emb,cx,qr,qx,tx,vap)

use micro_coms,  only: mza0, ncat, rxmin, enmlttab, dnfac, pwmasi, &
                       gnu, shedtab
use ccnbin_coms, only: nccntyp
use consts_coms, only: r8, alli
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0
integer, intent(in) :: lcat

integer, intent(inout) :: k1(11)
integer, intent(inout) :: k2(11)

real, intent(in) :: con_ccnx(mza0,nccntyp)

integer, intent(in) :: jhcat(mza0,ncat)

integer, intent(out) :: ict1(mza0,ncat)
integer, intent(out) :: ict2(mza0,ncat)

real, intent(out) :: wct1(mza0,ncat)
real, intent(out) :: wct2(mza0,ncat)

real, intent(inout) :: rx  (mza0,ncat)
real, intent(inout) :: cx  (mza0,ncat)
real, intent(inout) :: qr  (mza0,ncat)
real, intent(inout) :: qx  (mza0,ncat)
real, intent(inout) :: tx  (mza0,ncat)
real, intent(inout) :: emb (mza0,ncat)
real, intent(in)    :: vap (mza0,ncat)

integer :: k,lhcat,inc,idns
real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,ricetor6,rshed,rmltshed, &
        qrmltshed,fracmloss,dn

integer, parameter :: jflag = 1
real, parameter :: shedmass = 5.236e-7

! Collection can change the vertical range over which a hydrometeor
! category is present, so rediagnose k1(lcat) and k2(lcat).

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

! Return if there is no significant bulk mass for lcat category in this column

if (k1(lcat) > k2(lcat)) return

! Diagnose bulk mean mass and/or number concentration for lcat

 call enemb(lcat,jflag,k1,k2,con_ccnx,ict1,ict2,wct1,wct2,rx,cx,emb,vap)

! CLOUD, RAIN, and DRIZZLE categories

if (lcat == 1 .or. lcat == 2 .or. lcat == 8) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rxinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rxinv

! limit cloud, rain, drizzle temperature to range of about (-40C,40C)

         qx(k,lcat) = max(.5 * alli,min(1.5 * alli,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rx(k,lcat)

      endif

   enddo

! PRISTINE ICE category

elseif (lcat == 3) then

! Steve: Allow pristine ice to melt to cloud1 since we assume that smaller
! particles will melt first. Perhaps need a way in the future to treat
! melting of larger pristine ice into cloud2.

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin(lcat)) then

         rinv = 1. / rx(k,lcat)
         qx(k,lcat) = qr(k,lcat) * rinv

         call qtc(qx(k,lcat),tx(k,lcat),fracliq)

! If PRISTINE ICE is melting, transfer liquid to cloud category

         if (fracliq > 1.e-6) then

            rmelt = rx(k,lcat) * fracliq
            cmelt = cx(k,lcat) * fracliq

            rx(k,lcat) = rx(k,lcat) - rmelt
            cx(k,lcat) = cx(k,lcat) - cmelt
            qr(k,lcat) = 0.
            qx(k,lcat) = 0.
            tx(k,lcat) = 0.

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

         call qtc(qx(k,lcat),tx(k,lcat),fracliq)

! If SNOW or AGGREGATES are melting, transfer liquid plus some ice 
! to GRAUPEL category

         if (fracliq > 1.e-6) then

            rmelt = rx(k,lcat) * fracliq
            ricetor6 = min(rmelt, rx(k,lcat) - rmelt)

! change this??? move to rain instead ??? look at melting decisions in col2

            rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
            qr(k,lcat) = 0.
            qx(k,lcat) = 0.
            tx(k,lcat) = 0.
            
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

         call qtc(qx(k,lcat),tx(k,lcat),fracliq)

! If GRAUPEL is more than 95% melted, transfer all to RAIN category

         if (fracliq > 0.95) then

            rx(k,2) = rx(k,2) + rx(k,6)
            qr(k,2) = qr(k,2) + rx(k,6) * alli
            cx(k,2) = cx(k,2) + cx(k,6)

            rx(k,6) = 0.
            cx(k,6) = 0.
            qr(k,6) = 0.
            qx(k,6) = 0.
            tx(k,6) = 0.

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

         call qtc(qx(k,lcat),tx(k,lcat),fracliq)

! If HAIL is more than 95% melted, transfer all to RAIN category

         if (fracliq > 0.95) then

            rx(k,2) = rx(k,2) + rx(k,7)
            qr(k,2) = qr(k,2) + rx(k,7) * alli
            cx(k,2) = cx(k,2) + cx(k,7)

            rx(k,7) = 0.
            cx(k,7) = 0.
            qr(k,7) = 0.
            qx(k,7) = 0.
            tx(k,7) = 0.
         
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

            if (rx(k,7) > rxmin(7)) then
               qx(k,7) = qr(k,7) / rx(k,7)
               call qtc(qx(k,lcat),tx(k,lcat),fracliq)
            else
               qx(k,7) = 0.
               tx(k,7) = 0.
            endif
         endif

      endif

   enddo

endif

end subroutine x02

!===============================================================================

subroutine sedim2(iw0,lpw0,k1,k2,jhcat,dtl0, &
   voa, denfac, tair, thil0, theta0, dsed_thil, rhoi, rhoa, rhow, &
   cx, rx, qx, qr, emb, dmb, pcpvel, pcpfluxc, pcpfluxr, pcpfluxq, accpx, pcprx)

use micro_coms,  only: mza0, ncat, rxmin, cfmasi, pwmasi, cfvt, pwvt, jnmb, emb1
use consts_coms, only: r8, cpi, alviocp
use misc_coms,   only: io6
use mem_grid,    only: zm, dzt, dzit, zfacm2, zfacim2, arw0, arw, volti
use mem_ijtabs,  only: itab_w
use mem_sea,     only: sea, itab_ws
use mem_leaf,    only: land, itab_wl

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

integer, intent(in) :: jhcat(mza0,ncat)

real, intent(in) :: dtl0

real, intent(in)    :: voa   (mza0)
real, intent(in)    :: denfac(mza0)
real, intent(in)    :: tair     (mza0)
real, intent(in)    :: thil0    (mza0)
real, intent(in)    :: theta0   (mza0)
real, intent(inout) :: dsed_thil(mza0)
real, intent(in)    :: rhoi     (mza0)

real(r8), intent(inout) :: rhoa(mza0)
real(r8), intent(inout) :: rhow(mza0)

real, intent(inout) :: cx      (mza0,ncat)
real, intent(inout) :: rx      (mza0,ncat)
real, intent(inout) :: qx      (mza0,ncat)
real, intent(inout) :: qr      (mza0,ncat) 
real, intent(inout) :: emb     (mza0,ncat)
real, intent(inout) :: dmb     (mza0,ncat)
real, intent(inout) :: pcpvel  (mza0,ncat)
real, intent(inout) :: pcpfluxc(mza0,ncat)
real, intent(inout) :: pcpfluxr(mza0,ncat)
real, intent(inout) :: pcpfluxq(mza0,ncat)

real, intent(inout) :: accpx(ncat)
real, intent(inout) :: pcprx(ncat)

integer, parameter :: iplaws = 0

real, parameter :: alphasfc(ncat) = (/.001,.001,.010,.010,.010,.003,.001,.001/)

integer :: lcat, lhcat, k, kk, kw, nsea, nland, jws, iws, jwl, iwl

real :: sourcec, sourcer, sourceq
real :: zbotnew, areascale, fracwkk

real :: cxnew(mza0,ncat)
real :: rxnew(mza0,ncat)
real :: qrnew(mza0,ncat) 

real    :: dispemb,riemb,rsfc,qrsfc

  ! Loop over precipitation categories

  do lcat = 2,ncat
     if (k1(lcat) > k2(lcat)) cycle

     ! Loop over precipation source levels

     do k = k1(lcat),k2(lcat)
        lhcat = jhcat(k,lcat)

        ! Precipitation number, mass, and energy in source grid cell per m^2 of
        ! source grid cell BOTTOM horiz area

        sourcec = cx(k,lcat) * voa(k)
        sourcer = rx(k,lcat) * voa(k)
        sourceq = qr(k,lcat) * voa(k)

        ! Diameter of mean-mass hydrometeor

        dmb(k,lcat) = (emb(k,lcat) * cfmasi(lhcat))**pwmasi(lhcat)

        ! Fall velocity [m/s] of mean-mass hydrometeor, limited to max of 15 m/s

        ! Here determine which set of powerlaws to use: the original
        !  ones in RAMS or the ones from R.Carver adapted from Mitchell 1996.
        ! The Mitchell power laws are not based at sea level so we adjust the
        !  density factor based at 0.7 kg/m3 instead of 1.0 kg/m3.

        if (iplaws == 0) then
           pcpvel(k,lcat) = min(15., cfvt(lhcat) * dmb(k,lcat) ** pwvt(lhcat) * denfac(k))
        else
           pcpvel(k,lcat) = min(15., cfvt(lhcat) * dmb(k,lcat) ** pwvt(lhcat) * (0.7 * rhoi(k))**.362)
        endif

        ! New bottom height of source-cell precipitation after fall for 1 time step

        zbotnew = zm(k-1) - dtl0 * pcpvel(k,lcat)

        ! Loop over W levels that source grid cell precipitation can reach this time step

        do kk = k-1,lpw0-1,-1
           if (zm(kk) <= zbotnew) exit

           ! Horizontal area scale factor for current w(kk) level

           areascale = zfacim2(k-1) * zfacm2(kk)

           ! Fraction of source grid cell precipitation that crosses w(kk) level

           fracwkk = min(1.0, (zm(kk) - zbotnew) * dzit(k))

           ! Add source cell contribution to precipitation number flux [#/m^2],
           ! mass flux [kg/m^2], and internal energy flux [J/m^2] across w(kk) level

           pcpfluxc(kk,lcat) = pcpfluxc(kk,lcat) + sourcec * fracwkk * areascale
           pcpfluxr(kk,lcat) = pcpfluxr(kk,lcat) + sourcer * fracwkk * areascale
           pcpfluxq(kk,lcat) = pcpfluxq(kk,lcat) + sourceq * fracwkk * areascale
        enddo

     enddo ! k

     ! Apply number and mass transfers to obtain new hydrometeor concentrations
     ! Method should be positive definite, but check this.

     do k = k2(lcat),lpw0,-1
        cxnew(k,lcat) = cx(k,lcat) + volti(k,iw0) &
           * (pcpfluxc(k,lcat) * arw(k,iw0) - pcpfluxc(k-1,lcat) * arw(k-1,iw0))
        rxnew(k,lcat) = rx(k,lcat) + volti(k,iw0) &
           * (pcpfluxr(k,lcat) * arw(k,iw0) - pcpfluxr(k-1,lcat) * arw(k-1,iw0))
        qrnew(k,lcat) = qr(k,lcat) + volti(k,iw0) &
           * (pcpfluxq(k,lcat) * arw(k,iw0) - pcpfluxq(k-1,lcat) * arw(k-1,iw0))
     enddo

  enddo ! lcat

  ! Check for sea area beneath this atmospheric grid column

  nsea = itab_w(iw0)%nsea
  if (nsea > 0) then

     ! Loop over sea cells beneath this atmospheric grid column

     do jws = 1,nsea
        iws = itab_w(iw0)%isea(jws)
        kw = itab_ws(iws)%kw

        ! Loop over precipitation categories

        do lcat = 2,ncat
           if (k1(lcat) > k2(lcat)) cycle
           cxnew(kw,lcat) = cxnew(kw,lcat) - volti(kw,iw0) * pcpfluxc(kw-1,lcat) * sea%area(iws)
           rxnew(kw,lcat) = rxnew(kw,lcat) - volti(kw,iw0) * pcpfluxr(kw-1,lcat) * sea%area(iws)
           qrnew(kw,lcat) = qrnew(kw,lcat) - volti(kw,iw0) * pcpfluxq(kw-1,lcat) * sea%area(iws)

           pcprx(lcat) = pcprx(lcat) + pcpfluxr(kw-1,lcat) * sea%area(iws) / (arw0(iw0) * dtl0)
           accpx(lcat) = accpx(lcat) + pcpfluxr(kw-1,lcat) * sea%area(iws) / arw0(iw0)

           sea%pcpg (iws) = sea%pcpg (iws) + pcpfluxr(kw-1,lcat)
           sea%qpcpg(iws) = sea%qpcpg(iws) + pcpfluxq(kw-1,lcat)
           sea%dpcpg(iws) = sea%dpcpg(iws) + pcpfluxr(kw-1,lcat) * alphasfc(lcat)
        enddo

     enddo
  endif

  ! Check for land area beneath this atmospheric grid column

  nland = itab_w(iw0)%nland
  if (nland > 0) then

     ! Loop over land cells beneath this atmospheric grid column

     do jwl = 1,nland
        iwl = itab_w(iw0)%iland(jwl)
        kw = itab_wl(iwl)%kw

        ! Loop over precipitation categories

        do lcat = 2,ncat
           if (k1(lcat) > k2(lcat)) cycle
           cxnew(kw,lcat) = cxnew(kw,lcat) - volti(kw,iw0) * pcpfluxc(kw-1,lcat) * land%area(iwl)
           rxnew(kw,lcat) = rxnew(kw,lcat) - volti(kw,iw0) * pcpfluxr(kw-1,lcat) * land%area(iwl)
           qrnew(kw,lcat) = qrnew(kw,lcat) - volti(kw,iw0) * pcpfluxq(kw-1,lcat) * land%area(iwl)

           pcprx(lcat) = pcprx(lcat) + pcpfluxr(kw-1,lcat) * land%area(iwl) / (arw0(iw0) * dtl0)
           accpx(lcat) = accpx(lcat) + pcpfluxr(kw-1,lcat) * land%area(iwl) / arw0(iw0)

           land%pcpg (iwl) = land%pcpg (iwl) + pcpfluxr(kw-1,lcat)
           land%qpcpg(iwl) = land%qpcpg(iwl) + pcpfluxq(kw-1,lcat)
           land%dpcpg(iwl) = land%dpcpg(iwl) + pcpfluxr(kw-1,lcat) * alphasfc(lcat)
        enddo

     enddo
  endif

  ! Loop over precipitation categories

  do lcat = 2,ncat
     if (k1(lcat) > k2(lcat)) cycle

     ! Loop over precipation source levels

     do k = lpw0,k2(lcat)

        if (rxnew(k,lcat) < rxmin(lcat)) then
           rxnew(k,lcat) = 0.
           cxnew(k,lcat) = 0.
           qrnew(k,lcat) = 0.
        elseif (jnmb(lcat) == 5) then
           cxnew(k,lcat) = max(cxnew(k,lcat), rxnew(k,lcat) / emb1(lcat))
        endif

        rhoa(k) = rhoa(k) + rxnew(k,lcat) - rx(k,lcat)
        rhow(k) = rhow(k) + rxnew(k,lcat) - rx(k,lcat)

        dsed_thil(k) = dsed_thil(k) - thil0(k) * thil0(k)  &
           * (alviocp * (rxnew(k,lcat) - rx(k,lcat))  &
           - cpi * (qrnew(k,lcat) - qr(k,lcat)))  &
           / (max(tair(k), 253.) * theta0(k))

        ! Transfer "new" amounts to category arrays

        rx(k,lcat) = rxnew(k,lcat)
        cx(k,lcat) = cxnew(k,lcat)
        qx(k,lcat) = qrnew(k,lcat) / max(rxmin(lcat),rxnew(k,lcat))

     enddo

  enddo ! lcat

end subroutine sedim2

