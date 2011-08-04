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
subroutine thrmstr(iw0,lpw0,k1,k2, &
   press0,thil0,rhow,rhoi,exner0,tair,theta0,rhov,rhovstr,tairstrc, &
   rx,qx,sa)

use micro_coms,  only: mza0, ncat
use consts_coms, only: p00i, rocp, alvl, alvi, cpi4, cpi, cp253i
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw0,lpw0

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in)  :: press0  (mza0)
real, intent(in)  :: thil0   (mza0)
real, intent(in)  :: rhoi    (mza0)
real, intent(out) :: exner0  (mza0)
real, intent(out) :: tair    (mza0)
real, intent(out) :: theta0  (mza0)
real, intent(out) :: rhov    (mza0)
real, intent(out) :: rhovstr (mza0)
real, intent(out) :: tairstrc(mza0)

real(kind=8), intent(in) :: rhow(mza0)

real, intent(in) :: rx(mza0,ncat)
real, intent(in) :: qx(mza0,ncat)

real, intent(out) :: sa(mza0,9)

integer :: k,lcat
real :: fracliq,tcoal,tairstr

! automatic arrays

real :: rholiq(mza0)
real :: rhoice(mza0)
real :: qhydm (mza0)
real :: til   (mza0)

! Loop over whole column

do k = lpw0,mza0
   exner0(k) = (press0(k) * p00i) ** rocp  ! defined WITHOUT CP factor
   theta0(k) = thil0(k)
   tair(k) = theta0(k) * exner0(k)
   rhov(k) = rhow(k)
enddo

! Loop over levels that may have any type of condensate

do k = k1(11),k2(11)
   til(k) = thil0(k) * exner0(k)
   rholiq(k) = 0.
   rhoice(k) = 0.
enddo

! Loop over levels that may have cloud

do k = k1(1),k2(1)
   rholiq(k) = rholiq(k) + rx(k,1)
enddo

! Loop over levels that may have rain

do k = k1(2),k2(2)
   rholiq(k) = rholiq(k) + rx(k,2)
enddo

! Loop over levels that may have drizzle

do k = k1(8),k2(8)
   rholiq(k) = rholiq(k) + rx(k,8)
enddo

! Loop over levels that may have pristine ice

do k = k1(3),k2(3)
   rhoice(k) = rhoice(k) + rx(k,3)
enddo

! Loop over levels that may have snow

do k = k1(4),k2(4)
   rhoice(k) = rhoice(k) + rx(k,4)
enddo

! Loop over levels that may have aggregates

do k = k1(5),k2(5)
   rhoice(k) = rhoice(k) + rx(k,5)
enddo

! Loop over levels that may have graupel 
! (qtc diagnoses graupel temp and liquid fraction)

do k = k1(6),k2(6)
   call qtc(qx(k,6),tcoal,fracliq)
   rholiq(k) = rholiq(k) + rx(k,6) * fracliq
   rhoice(k) = rhoice(k) + rx(k,6) * (1. - fracliq)
enddo

! Loop over levels that may have hail
! (qtc diagnoses graupel temp and liquid fraction)

do k = k1(7),k2(7)
   call qtc(qx(k,7),tcoal,fracliq)
   rholiq(k) = rholiq(k) + rx(k,7) * fracliq
   rhoice(k) = rhoice(k) + rx(k,7) * (1. - fracliq)
enddo

! Loop over levels that may have any type of condensate

do k = k1(11),k2(11)
   qhydm(k) = alvl * rholiq(k) + alvi * rhoice(k)
   rhovstr(k) = rhow(k) - rholiq(k) - rhoice(k)
   sa(k,1) = til(k) * qhydm(k) / (1.e-12 + rholiq(k) + rhoice(k)) ! stays the same

   if (tair(k) > 253.) then

!ORIG   tairstr = .5 * (til(k) + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
! Change in tairstr computation since qhydm is now J/m^3 instead of J/kg:    

      tairstr = .5 &
         * (til(k) + sqrt(til(k) * (til(k) + cpi4 * qhydm(k) * rhoi(k))))
      sa(k,1) = sa(k,1) * cpi / (2. * tairstr - til(k)) ! stays the same

   else

!ORIG   tairstr = til(k) * (1. + qhydm(k) * cp253i)
! Change in tairstr computation since qhydm is now J/m^3 instead of J/kg:    

      tairstr = til(k) * (1. + qhydm(k) * rhoi(k) * cp253i)
      sa(k,1) = sa(k,1) * cp253i ! stays the same

   endif
   tairstrc(k) = tairstr - 273.15
     
enddo

return
end subroutine thrmstr

!===============================================================================

subroutine diffprep(iw0,lcat,k1,k2, &
   pi4dt,jhcat,sa,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz, &
   rx,cx,qr,emb,rhoa,rhov,rhovsrefp,rdynvsci,vapdif,thrmcon,sumuy,sumuz)

use micro_coms, only: rxmin, frefac1, pwmasi, frefac2, cdp1, sl, sj, sc, sk, &
                      mza0, ncat
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0,lcat

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

real, intent(in) :: pi4dt

integer, intent(in) :: jhcat(mza0,ncat)

real, intent(in) :: sa(mza0,9)

real, intent(out)   :: sb(mza0,ncat)
real, intent(out)   :: sd(mza0,ncat)
real, intent(out)   :: se(mza0,ncat)
real, intent(out)   :: sf(mza0,ncat)
real, intent(out)   :: sg(mza0,ncat)
real, intent(inout) :: sh(mza0,ncat)
real, intent(inout) :: sm(mza0,ncat)
real, intent(out)   :: ss(mza0,ncat)
real, intent(out)   :: su(mza0,ncat)
real, intent(out)   :: sw(mza0,ncat)
real, intent(out)   :: sy(mza0,ncat)
real, intent(out)   :: sz(mza0,ncat)

real, intent(in) :: rx (mza0,ncat)
real, intent(in) :: cx (mza0,ncat)
real, intent(in) :: qr (mza0,ncat)
real, intent(in) :: emb(mza0,ncat)

real, intent(in) :: rhov     (mza0)
real, intent(in) :: rhovsrefp(mza0,2)
real, intent(in) :: rdynvsci (mza0)
real, intent(in) :: vapdif   (mza0)
real, intent(in) :: thrmcon  (mza0)

real(kind=8), intent(in) :: rhoa(mza0)

real, intent(inout) :: sumuy(mza0)
real, intent(inout) :: sumuz(mza0)

integer :: k,mynum,if1,if4,if6,if8,lhcat
real :: fre,scdei

real, dimension(mza0) :: ttest ! automatic array

!CODE BASED ON WALKO ET AL 2000
!EFFICIENT COMPUTATION OF VAPOR AND HEAT DIFFUSION BETWEEN HYDROMETEORS
!IN A NUMERICAL MODEL

! DETERMINES WHETHER TO CALCULATE USING LIQUID OR ICE "SA" ARRAYS

if (lcat <= 2 .or. lcat == 8) then
   if1 = 1
   if4 = 4
   if6 = 6
   if8 = 8
else
   if1 = 2
   if4 = 5
   if6 = 7
   if8 = 9
endif

do k = k1(lcat),k2(lcat)
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) < rxmin(lcat)) cycle

   fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat) &
      + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

   sb(k,lcat) = cx(k,lcat) * fre * pi4dt  ! stays the same (rhoa factor removed)
   su(k,lcat) = vapdif(k) * sb(k,lcat)    ! stays the same
   sd(k,lcat) = sh(k,lcat) * rx(k,lcat)   ! x rhoa
   se(k,lcat) = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k) * rhoa(K)
    ! se picked up rhoa factor (rhoa factor had to be inserted in second term)
   sf(k,lcat) = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)  ! stays the same
   sg(k,lcat) = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3) &
              + sj(lcat) * qr(k,lcat) ! x rhoa
!     + lambda_j [Joules/m^3 added by radiative heating this timestep]
   scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat)) ! x (1/rhoa)
   ss(k,lcat) = sf(k,lcat) * scdei   ! x (1/rhoa)
   sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei  ! stays the same
   ttest(k) = ss(k,lcat) * rhov(k) + sw(k,lcat)  ! stays the same

! FOR ALL ICE HYDROS, "SM" IS 1
! IF PRISTINE,SNOW,AGG TTEST >= ZERO, "SH" IS 1

   if (lcat >= 3 .and. lcat <= 5) then
      if (ttest(k) >= 0.) then
         sm(k,lcat) = 0.
         sh(k,lcat) = 1.
         sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
         scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
         ss(k,lcat) = sf(k,lcat) * scdei
         sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      else
         sm(k,lcat) = 1.
      endif
   endif

! FOR MIXED-PHASE HYDROMETEORS, "SM" IS 0

   if (lcat == 6 .or. lcat == 7) then
      if (ttest(k) >= 0.) then
         sm(k,lcat) = 0.
      else
         sm(k,lcat) = 1.
      endif
   endif

   sy(k,lcat) = rhovsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
   sz(k,lcat) = 1. - rhovsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
   sumuy(k) = sumuy(k) + su(k,lcat) * sy(k,lcat)
   sumuz(k) = sumuz(k) + su(k,lcat) * sz(k,lcat)

enddo

return
end subroutine diffprep

!===============================================================================

subroutine vapdiff(j1,j2,rhov,rhovstr,sumuy,sumuz)

use micro_coms, only: mza0
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

real, intent(out) :: rhov   (mza0)
real, intent(in)  :: rhovstr(mza0)
real, intent(in)  :: sumuy  (mza0)
real, intent(in)  :: sumuz  (mza0)

integer :: k

do k = j1,j2
   rhov(k) = (rhovstr(k) + sumuy(k)) / (1.0 + sumuz(k))
enddo

return
end subroutine vapdiff

!===============================================================================

subroutine vapflux(iw0,lcat,k1,k2, &
   jhcat,sa,sd,se,sf,sg,sm,ss,su,sw,sy,sz,rx,cx,qx,qr,tx,vap, &
   rhovsrefp,rhov,rhovstr,sumuy,sumuz,sumvr,con_ccnx,con_gccnx)

use micro_coms, only: mza0, ncat, rxmin, sc, sk, jnmb, enmlttab, iccnlev
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0,lcat

integer, intent(in) :: k1(11)
integer, intent(in) :: k2(11)

integer, intent(in) :: jhcat(mza0,ncat)

real, intent(in) :: sa(mza0,9)

real, intent(in) :: sd(mza0,ncat)
real, intent(in) :: se(mza0,ncat)
real, intent(in) :: sf(mza0,ncat)
real, intent(in) :: sg(mza0,ncat)
real, intent(in) :: sm(mza0,ncat)
real, intent(in) :: ss(mza0,ncat)
real, intent(in) :: su(mza0,ncat)
real, intent(in) :: sw(mza0,ncat)
real, intent(in) :: sy(mza0,ncat)
real, intent(in) :: sz(mza0,ncat)

real, intent(inout) :: rx (mza0,ncat)
real, intent(inout) :: cx (mza0,ncat)
real, intent(out)   :: qx (mza0,ncat)
real, intent(out)   :: qr (mza0,ncat)
real, intent(out)   :: tx (mza0,ncat)
real, intent(out)   :: vap(mza0,ncat)

real, intent(in)    :: rhovsrefp(mza0,2)
real, intent(inout) :: rhov     (mza0)
real, intent(in)    :: rhovstr  (mza0)
real, intent(inout) :: sumuy    (mza0)
real, intent(inout) :: sumuz    (mza0)
real, intent(inout) :: sumvr    (mza0)
real, intent(inout) :: con_ccnx (mza0)
real, intent(inout) :: con_gccnx(mza0)

integer :: k,if1,if4
real :: rxx

!-----------------------------------------------------------------------
! Variables for aerosol restoration upon droplet evaporation

integer :: rgccn1,rgb,ic,ct
real :: cnnum,cnmass,jrg1,jrg2,fracmass,cxloss,rg

real :: rg_ccn(14) ! Given median radius of CCN spectra for mass calculation

data rg_ccn / 0.01e-4,0.015e-4,0.02e-4,0.03e-4,0.04e-4,0.06e-4,0.08e-4 &
             ,0.12e-4,0.160e-4,0.24e-4,0.32e-4,0.48e-4,0.64e-4,0.96e-4 /
!-----------------------------------------------------------------------

! UPDATES MIXING RATIO AND HEAT DUE TO FLUX OF VAPOR
! ALSO LINKED TO WALKO ET AL 2000

if (lcat <= 2 .or. lcat == 8) then
   if1 = 1
   if4 = 4
else
   if1 = 2
   if4 = 5
endif

do k = k1(lcat),k2(lcat)

   if (rx(k,lcat) < rxmin(lcat)) cycle

   tx(k,lcat) = (ss(k,lcat) * rhov(k) + sw(k,lcat)) * sm(k,lcat)
   vap(k,lcat) = su(k,lcat) * (rhov(k) + sa(k,if4) - rhovsrefp(k,if1) * tx(k,lcat))

! Do this section if vapor transfer does NOT deplete all of LCAT category

   if (vap(k,lcat) > -rx(k,lcat)) then

      rxx = rx(k,lcat) + vap(k,lcat)

      if (sm(k,lcat) > .5) then
         qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
         qr(k,lcat) = qx(k,lcat) * rxx
      else
         qx(k,lcat) = (rhov(k) * sf(k,lcat) + sg(k,lcat) &
                    - tx(k,lcat) * se(k,lcat)) / sd(k,lcat)
         qx(k,lcat) = min(350000.,max(-100000.,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rxx
      endif

   endif

! Do this section if vapor transfer DOES deplete all of LCAT category

! Also do this section if LCAT is pristine ice and it totally melts:
! (evaporate it too).

   if ((vap(k,lcat) <= -rx(k,lcat)) .or. &
       (lcat == 3 .and. qx(k,lcat) > 330000.)) then

      sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
      sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
      sumvr(k) = sumvr(k) + rx(k,lcat)

      rhov(k) = (rhovstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k))

      vap(k,lcat) = - rx(k,lcat)
      tx(k,lcat) = 0.
      rx(k,lcat) = 0.
      qx(k,lcat) = 0.
      qr(k,lcat) = 0.

! If we are doing nucleation scavenging, return ccn or gccn when complete
! evaporation occurs

      if (iccnlev == 1) then
         if (lcat == 1 .or. (lcat == 8 .and. jnmb(lcat) < 5)) &
            con_ccnx(k) = con_ccnx(k) + cx(k,lcat)

         if (lcat == 8 .and. jnmb(lcat) >= 5) &
            con_gccnx(k) = con_gccnx(k) + cx(k,lcat)
      endif

      cx(k,lcat) = 0.

   else

! If we are doing nucleation scavenging, return ccn or gccn when partial
! evaporation occurs

      if (vap(k,lcat) < 0.) then
         fracmass = min(1.,-vap(k,lcat) / rx(k,lcat))
         cxloss = cx(k,lcat) * enmlttab( int(200.*fracmass)+1, jhcat(k,lcat) )
         cx(k,lcat) = cx(k,lcat) - cxloss

         if (iccnlev == 1) then
            if (lcat == 1 .or. (lcat == 8 .and. jnmb(lcat) < 5)) &
               con_ccnx(k) = con_ccnx(k) + cxloss

            if (lcat == 8 .and. jnmb(lcat) >= 5) &
               con_gccnx(k) = con_gccnx(k) + cxloss
         endif
      endif

      rx(k,lcat) = rxx

   endif

enddo

return
end subroutine vapflux

!===============================================================================

subroutine psxfer(iw0,j1,j2,jhcat,vap,rx,cx,qx,qr)

use micro_coms, only: rxmin, dnfac, pwmasi, gam, dps2, dps, gnu, gamn1, &
                      pwmas, dpsmi, mza0, ncat, emb0, emb1
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw0,j1,j2

integer, intent(in) :: jhcat (mza0,ncat)

real, intent(in)    :: vap(mza0,ncat)
real, intent(inout) :: rx (mza0,ncat)
real, intent(inout) :: cx (mza0,ncat)
real, intent(in)    :: qx (mza0,ncat)
real, intent(inout) :: qr (mza0,ncat)

integer :: k,lhcat,it
real :: embx,dn,xlim,dvap,dqr,dnum

! Variables for Carrio adjustment to number transfer between snow and pristine ice

real :: old_m,old_c,old_r,prelim_m,delta_c, delta_r


do k = j1,j2

   if (vap(k,3) > 0.) then

      lhcat = jhcat(k,3)
      embx = max(emb0(3),min(emb1(3),rx(k,3) / max(1.e-9,cx(k,3))))

      dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
      it = nint(dn * 1.e6)

      xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.) &
         / (gamn1(3) * pwmas(lhcat) * dn ** 2)

      dvap = min(rx(k,3),vap(k,3) * (xlim + gam(it,1) / gamn1(3))) ! x rhoa
      dqr = dvap * qx(k,3)                                         ! x rhoa
      dnum = dvap * min(dpsmi(lhcat),1./embx)                      ! x rhoa

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

   elseif (vap(k,4) < 0. .and. rx(k,4) >= 0.) then

      lhcat = jhcat(k,4)
      embx = max(emb0(4),min(emb1(4),rx(k,4) / max(1.e-9,cx(k,4))))

      dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
      it = nint(dn * 1.e6)

      xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.) &
         / (gamn1(4) * pwmas(lhcat) * dn ** 2)

      dvap = max(-rx(k,4),vap(k,4) * xlim)     ! x rhoa
      dqr = dvap * qx(k,4)                     ! x rhoa
      dnum = dvap * max(dpsmi(lhcat),1./embx)  ! x rhoa

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

! Label/check this for rho factor

      !Carrio 2003: Better xfer pristine to snow to prevent subroutine
      !enemb from artificially creating pristine number concentration
      !when a bounds problem appears with the pristine ice cutoff size.
      !Compare preliminary calculation pristine mean mass with maximum.

!      delta_r = 0.0
      prelim_m = rx(k,3) / max(cx(k,3), 1.e-9)

      if (prelim_m > emb1(3)) then
          old_c = cx(k,3) + dnum
!          old_m = (rx(k,3) + dvap) / (cx(k,3) + dnum)
!          old_r = rx(k,3) + dvap
          delta_r = rx(k,3) - old_c * emb1(3)
          delta_c = delta_r / emb1(3) ! delta_c only adds to snow
          rx(k,3) = rx(k,3) - delta_r
          cx(k,3) = old_c
          rx(k,4) = rx(k,4) + delta_r
          cx(k,4) = cx(k,4) + delta_c
      endif

   endif

enddo

return
end subroutine psxfer

!===============================================================================

subroutine newtemp(j1,j2, &
   tairstrc,rhoi,rhovstr,rhov,exner0,tairc,tair,theta0,rhovslair,rhovsiair,sa)

use micro_coms, only: mza0
use misc_coms,  only: io6

implicit none

integer, intent(in) :: j1
integer, intent(in) :: j2

real, intent(in)  :: tairstrc (mza0)
real, intent(in)  :: rhoi     (mza0)
real, intent(in)  :: rhovstr  (mza0)
real, intent(in)  :: rhov     (mza0)
real, intent(in)  :: exner0   (mza0)
real, intent(out) :: tairc    (mza0)
real, intent(out) :: tair     (mza0)
real, intent(out) :: theta0   (mza0)
real, intent(out) :: rhovslair(mza0)
real, intent(out) :: rhovsiair(mza0)

real, intent(in) :: sa (mza0,9)

real, external :: rhovsl,rhovsi

integer :: k

do k = j1,j2
   tairc(k) = tairstrc(k) + sa(k,1) * rhoi(k) * (rhovstr(k) - rhov(k)) ! rhoi inserted
   tair(k)  = tairc(k) + 273.15
   theta0(k) = tair(k) / exner0(k)
   rhovslair(k) = rhovsl(tairc(k))
   rhovsiair(k) = rhovsi(tairc(k))
enddo

return
end subroutine newtemp
