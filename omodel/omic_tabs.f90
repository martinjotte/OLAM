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
subroutine haznuc()

use micro_coms,  only: nthz, dthz, nrhhz, drhhz, frachz, cnparm
use misc_coms,   only: io6
use consts_coms, only: r8

implicit none

integer  :: ithz,irhhz,k
real     :: denccn,gnuccn,dnccn,rhhz,c1hz,c2hz,c3hz,bhz,thz
real(r8) :: dum, dccn, dm, sum, y, ddccn

real, external :: gammln

!  Haze nucleation table

denccn = 1.769       !Density of ammonium sulfate
gnuccn = 2.          !Saleeby(02-21-2007) Originally = 1.
dnccn  = 2. * cnparm !Saleeby(02-21-2007) Originally = .075e-4 cm
ddccn  = .005d-4     !Increment [cm] of CCN diameter in numerical integration
                     !   over gamma size distribution

do ithz = 1,nthz
   thz = -60. + dthz * real(ithz - 1)
   do irhhz = 1,nrhhz
      rhhz = 0.82 + drhhz * real(irhhz - 1)
      c1hz = (3.14159 * denccn / 6.) ** (-0.3333333)
      c2hz = -14.65 - 1.045 * thz
      c3hz = -492.35 - 8.34 * thz - 0.0608 * thz ** 2
      bhz = min(38., max(-38., c2hz + c3hz * (1. - rhhz)))
      dm = c1hz * 10. ** (-bhz / 6.)

      sum  = 0.0_r8
      dccn = 0.0_r8
      do k = 1, 200
         dccn = dccn + ddccn
         y = dccn / dnccn
         dum = min(100.0_r8, (dccn / dm) ** 6)
         sum = sum + y ** (gnuccn - 1.) * exp(-y) * (1.0_r8 - exp(-dum))
      enddo
      frachz(irhhz,ithz) = sum * ddccn / (exp(gammln(gnuccn)) * dnccn)

   enddo
enddo

return
end subroutine haznuc

!===============================================================================

subroutine homfrzcl(dtl,mrl)

use micro_coms,  only: ntc, dtc, ndnc, ddnc, fracc, gnu
use misc_coms,   only: io6
use consts_coms, only: r8

implicit none

integer, intent(in) :: mrl
real,    intent(in) :: dtl

integer  :: itc, k, idnc
real     :: tc, y
real(r8) :: sum, dc, ddc, ajlso, v1, dnc, dum1, dum2
real, external :: gammln

!  Make table for homogeneous freezing of cloud droplets
!  Uses characteristic diameter instead of true diameter

! ndnc = 11
! ddnc = 2.e-6

ddc = 0.5d-6

do itc = 1, ntc

   tc = -50. + dtc * real(itc-1)
   y = -(606.3952+tc*(52.6611+tc*(1.7439+tc*(.0265+tc*1.536e-4))))
   ajlso = 1.d6 * 10. ** y

   do idnc = 1, ndnc

      dnc = ddnc * real(idnc,r8)
      sum = 0.0_r8
      dc  = 0.0_r8

      do k = 1, 2000
         dc = dc + ddc
         v1 = 0.523599_r8 * dc ** 3
         dum1 = dc / dnc
         dum2 = min(100.0_r8, ajlso * v1 * dtl)
         sum = sum + dum1 ** (gnu(1) - 1.) * exp(-dum1) * (1.0_r8 - exp(-dum2))
      enddo

      fracc(idnc,itc,mrl) = sum * ddc / (exp(gammln(gnu(1))) * dnc)
   enddo

enddo

return
end subroutine homfrzcl

!===============================================================================

subroutine mksedim_tab(m1,zm,dzt,dzit)

use micro_coms, only: sedtime0, sedtime1, mza0, zmf, dztf, dzitf, maxkfall,  &
                      nhcat, lcat_lhcat, dispemb0, cfvt, emb0, cfmasi,  &
                      pwvt, pwmasi, dispemb0i, dispemb1, ch2, nembfall, gnu,  &
                      pwmas, pcpfillc, pcpfillr, emb1
use misc_coms,  only: io6

implicit none

integer, intent(in) :: m1
real, intent(in) :: zm(m1),dzt(m1),dzit(m1)

integer, parameter :: nbin=50
integer :: iembs,lcat,lhcat,k,kkf,ibin,kk,jbin

real :: dmbodn,diam0,diam1,fac1,fac3,sumc,sumr,diam,fac2,fac4  &
       ,disp,ztopnew,zbotnew,fallin,delzsfc,dispemb,dispmax,dispmx
real, external :: gammln,gammp

real :: cbin   (nbin)
real :: rbin   (nbin)
real :: reldisp(nbin)

! To cover the possible range over all timesteps, define sedtime0 and sedtime1
! here as 0.1 seconds and 3000 seconds.  The former is supposed to be
! less than 0.7 of the shortest timestep on any grid (sqrt(rhoi) never exceeds
! 0.7) and the latter is the longest timestep expected to ever be used (600
! seconds) times a factor of 5 for the largest value of sqrt(rhoi).

sedtime0 = .1
sedtime1 = 3000.
dispmax = 500.

! Fill zmf and ztf arrays to easily handle "underground" points

do k = 1,mza0
   zmf(k) = zm(k)
   dztf(k) = dzt(k)
   dzitf(k) = dzit(k)
enddo

do k = 0,2-maxkfall,-1
   zmf(k) = 2. * zmf(k+1) - zmf(k+2)
   dztf(k+1) = zmf(k+1) - zmf(k)
   dzitf(k+1) = 1. / dztf(k+1)
enddo

! Loop over hydrometeor categories

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   dispemb0(lhcat) = sedtime0 * cfvt(lhcat)  &
      * (emb0(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))
      
   dispemb0i(lhcat) = 1. / dispemb0(lhcat)

   dispemb1(lhcat) = sedtime1 * cfvt(lhcat)  &
      * (emb1(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

!Bob (10/24/00):  Limit dispemb1 to a maximum of dispmax

   if (dispemb1(lhcat) > dispmax) dispemb1(lhcat) = dispmax

   ch2(lhcat) = real(nembfall-1) / log10(dispemb1(lhcat) * dispemb0i(lhcat))

! Loop over bins, filling them with fractional number, fractional mass,
! and displacement quotient relative to emb.

   dmbodn = (exp(gammln(gnu(lcat) + pwmas(lhcat))  &
      - gammln(gnu(lcat)))) ** pwmasi(lhcat)
   diam0 = 0.06 * dmbodn
   diam1 = 1.0 * dmbodn
   fac1 = gammp(gnu(lcat),diam0)
   fac3 = gammp(gnu(lcat) + pwmas(lhcat),diam0)
   sumc = 0.
   sumr = 0.

   do jbin = 1,nbin

      diam = diam0 * (diam1 / diam0) ** (real(jbin)/real(nbin))
      fac2 = gammp(gnu(lcat),diam)
      fac4 = gammp(gnu(lcat) + pwmas(lhcat),diam)
      cbin(jbin) = fac2 - fac1
      rbin(jbin) = fac4 - fac3
      fac1 = fac2
      fac3 = fac4
      sumc = sumc + cbin(jbin)
      sumr = sumr + rbin(jbin)
      reldisp(jbin) = diam ** pwvt(lhcat)

   enddo

   do jbin = 1,nbin
      cbin(jbin) = cbin(jbin) / sumc
      rbin(jbin) = rbin(jbin) / sumr
   enddo

! Loop over displacement distance for size emb.

   do iembs = 1,nembfall
      dispemb = dispemb0(lhcat)  &
         * (dispemb1(lhcat) * dispemb0i(lhcat))  &
         ** (real(iembs-1) / real(nembfall-1))

! Zero out concentration and mass fill arrays and surface precip array
! before accumulation.

      do k = 1,mza0
         do kkf = 1,maxkfall
            pcpfillc(k,kkf,iembs,lhcat) = 0.
            pcpfillr(k,kkf,iembs,lhcat) = 0.
         enddo
      enddo

! Loop over vertical grid index.

      do k = 2,mza0

!Bob (10/24/00):  Limit disp distance to (maxkfall-1) levels

         dispmx = min(dispmax,zmf(k-1) - zmf(k-maxkfall))

! Loop over bins

         do ibin = 1,nbin
            disp = dispemb * reldisp(ibin)
            if (disp > dispmx) disp = dispmx

            ztopnew = zmf(k) - disp
            zbotnew = zmf(k-1) - disp

! Loop over grid cells that a parcel falls into, including the one it starts from.

            do kkf = 1,maxkfall

               kk = k + 1 - kkf
               if (zbotnew > zmf(kk)) go to 50

               if (ztopnew <= zmf(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzitf(kk) *  &
                     (min(zmf(kk),ztopnew) - max(zmf(kk-1),zbotnew))
               endif

               pcpfillc(k,kkf,iembs,lhcat) = pcpfillc(k,kkf,iembs,lhcat)  &
                  + fallin * cbin(ibin)

               pcpfillr(k,kkf,iembs,lhcat) = pcpfillr(k,kkf,iembs,lhcat)  &
                  + fallin * rbin(ibin)

            enddo

50          continue

         enddo

      enddo
   enddo
enddo

return
end subroutine mksedim_tab

!===============================================================================

subroutine tabmelt()

use micro_coms, only: nhcat, lcat_lhcat, gnu, rmlttab, enmlttab, ndns,  &
                      shedtab, cfmas, pwmas, cfvt, pwvt, shapefac, ninc, ncat
use misc_coms,  only: io6

implicit none

integer, parameter :: nbins=500

integer :: lhcat,lcat,ndns1,ibin,inc,iter,idns
real :: dn,gammaa,totfmg,totmass,vtx,fre,totqm,qmgoal,qmnow,totmdqdt,deltat  &
       ,pliqmass,picemass,critmass,vk

real, dimension(nbins) :: db,fmg,pmass,binmass,dqdt,q
real, dimension(ncat) :: dmean
real, external :: gammln

data vk/0.2123e-04/

! Choose some reasonable mean diameter for all hydrometeor categories to be 
! used for generating enmlttab.  Actual sizes are not important since enmlttab 
! only contains relative number vs. mass melting information.

data dmean/20.e-6,500.e-6,30.e-6,500.e-6,500.e-6,500.e-6,8000.e-6,60.e-6/

! Loop over all hydrometeor categories, but enmlttab values for cloud,
! drizzle, and rain are never used.

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   dn = dmean(lcat) / gnu(lcat)
   gammaa = exp(gammln(gnu(lcat)))

   rmlttab(1) = 0.0
   rmlttab(ninc) = 1.0
   enmlttab(1,lhcat) = 0.0
   enmlttab(ninc,lhcat) = 1.0

   ndns1 = 1
   if (lcat == 7) ndns1 = ndns

! Loop over multiple characteristic diameters only for hail category
! and single characteristic diameter for all other categories

   do idns = 1,ndns1
      shedtab(1,idns) = 0.0
      shedtab(ninc,idns) = 0.0

! Re-define characteristic diameter only for hail category

      if (ndns1 > 1) dn = 1.e-3 * real(idns) / gnu(lcat)

      totfmg = 0.
      totmass = 0.

      do ibin = 1,nbins
         db(ibin) = 0.02 * dn * (real(ibin) - 0.5)
         fmg(ibin) = (db(ibin) / dn) ** (gnu(lcat) - 1.)  &
            / (dn * gammaa) * exp(-db(ibin) / dn)
         totfmg = totfmg + fmg(ibin)
         q(ibin) = 0.
         pmass(ibin) = cfmas(lhcat) * db(ibin) ** pwmas(lhcat)
         binmass(ibin) = pmass(ibin) * fmg(ibin)
         totmass = totmass + binmass(ibin)
         vtx = cfvt(lhcat) * db(ibin) ** pwvt(lhcat)
         fre = (1.0 + 0.229 * sqrt(vtx * db(ibin) / vk))  &
            * shapefac(lhcat)
         dqdt(ibin) = db(ibin) ** (1. - pwmas(lhcat)) * fre
      enddo

      totqm = totmass * 80.

! Loop over inc value (representing total liquid fraction)

      do inc = 2,ninc-1
         qmgoal = totqm * real(inc-1) / real(ninc-1)

         do iter = 1,2
            qmnow = 0.
            totmdqdt = 0.
            do ibin = 1,nbins
               if (q(ibin) < 79.9999) then
                  totmdqdt = totmdqdt + binmass(ibin) * dqdt(ibin)
               endif
               qmnow = qmnow + q(ibin) * binmass(ibin)
            enddo
            deltat = max(0.,(qmgoal - qmnow) / totmdqdt)
            do ibin = 1,nbins
               q(ibin) = min(80.,q(ibin) + dqdt(ibin) * deltat)
            enddo
         enddo

! For idns = 7 (an intermediate value that only happens once over idns loop and
! only for hail), compute melted mixing ratio (rmlttab) from totally-melted bins.

         if (idns == 7) then
            rmlttab(inc) = 0.0

            do ibin = 1,nbins
               if (q(ibin) > 79.9) then
                  rmlttab(inc) = rmlttab(inc) + binmass(ibin)
               endif
            enddo

            rmlttab(inc) = rmlttab(inc) / totmass
         endif

! Compute melted number (enmlttab) from totally-melted bins.  This is done
! for only one value of idns for each category.

         if (idns == 7 .or. ndns1 == 1) then
            enmlttab(inc,lhcat) = 0.0

            do ibin = 1,nbins
               if (q(ibin) > 79.9) then
                  enmlttab(inc,lhcat) = enmlttab(inc,lhcat) + fmg(ibin)
               endif
            enddo

            enmlttab(inc,lhcat) = enmlttab(inc,lhcat) / totfmg
         endif

! For hail category only, compute shedded mixing ratio (shedtab) from
! partially-melted bins.  This is done for all values of idns.

         if (lcat == 7) then
            shedtab(inc,idns) = 0.0
!                  do ibin = kbin,nbins

            do ibin = 1,nbins
               if (q(ibin) <= 79.9) then
                  pliqmass = pmass(ibin) * q(ibin) / 80.
                  picemass = pmass(ibin) - pliqmass
                  critmass = .268e-3 + .1389 * picemass
                  shedtab(inc,idns) = shedtab(inc,idns)  &
                     + max(0.0, pliqmass - critmass) * fmg(ibin)
               endif
            enddo

            shedtab(inc,idns) = shedtab(inc,idns) / totmass
         endif

      enddo
   enddo
enddo
return
end subroutine tabmelt

!===============================================================================

subroutine mkcoltb_brute()

use micro_coms, only: nhcat, lcat_lhcat, gnu, pwmas, emb0, cfmasi, pwmasi,  &
                      emb1, ipairc, ipairr, nembc, cfvt, pwvt, cfmas,  &
                      coltabc, coltabr, jnmb
use misc_coms,  only: io6

implicit none

integer, parameter :: ndx=100,ndy=100

integer :: ihx,ix,ihy,iy,iemby,iembx,idx,idy,ipc,iprxy,ipryx

real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gxn,gyn,gym  &
       ,dnminy,dnmaxy,dny,dnx,bint,sum_num,sum_xmass,sum_ymass,vx,vy,dx,dy  &
       ,fgamx,emx,dx1,dx2,fgamy,emy,dy1,dy2,dyhi,dylo

real, external :: gammln, efc

! Loop over colliding category X

do ihx = 1,nhcat
   ix = lcat_lhcat(ihx)

   gxm = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
   dnminx = ((emb0(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dnmaxx = ((emb1(ix) * cfmasi(ihx)) * gxm) ** pwmasi(ihx)
   dxlo = .01 * dnminx
   dxhi = 10. * dnmaxx

! Loop over colliding category Y

   do ihy = 1,nhcat

! Check if colliding pair (X,Y) is to be considered.  

      ipc   = ipairc(ihx,ihy)
      iprxy = ipairr(ihx,ihy)
      ipryx = ipairr(ihy,ihx)
      
! Bob (11/11/2009): We only need to check ipairc because 
! ipairr(IHX,IHY) and ipairr(IHY,IHX) are now filled together.

      if (ipc == 0) cycle

      iy = lcat_lhcat(ihy)

      gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
      dnminy = ((emb0(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)
      dnmaxy = ((emb1(iy) * cfmasi(ihy)) * gym) ** pwmasi(ihy)
      dylo = .01 * dnminy
      dyhi = 100. * dnmaxy

! Loop over all mean-mass table values for category Y

      do iemby = 1,nembc
         dny = dnminy * (dnmaxy / dnminy) ** (real(iemby-1) / real(nembc-1))

! Loop over all mean-mass table values for category X

         do iembx = 1,nembc
            dnx = dnminx * (dnmaxx / dnminx) ** (real(iembx-1) / real(nembc-1))

            sum_num = 0.    ! Initialize integral number sum to zero
            sum_xmass = 0.  ! Initialize integral xmass sum to zero
            sum_ymass = 0.  ! Initialize integral ymass sum to zero
            
! Loop over spectrum of Y diameters for current mean-mass Y value

            do idy = 1,ndy
               dy1 = dylo * (dyhi / dylo) ** (real(idy-1) / real(ndy))
               dy2 = dylo * (dyhi / dylo) ** (real(idy) / real(ndy))
               dy = .5 * (dy1 + dy2)
               vy = cfvt(ihy) * dy ** pwvt(ihy)
               emy = cfmas(ihy) * dy ** pwmas(ihy)

               gyn = exp(gammln(gnu(iy)))
               fgamy = (dy / dny) ** (gnu(iy) - 1) * exp(-dy / dny)  &
                  / (gyn * dny)

! Loop over spectrum of X diameters for current mean-mass X value

               do idx = 1,ndx
                  dx1 = dxlo * (dxhi / dxlo) ** (real(idx-1) / real(ndx))
                  dx2 = dxlo * (dxhi / dxlo) ** (real(idx) / real(ndx))
                  dx = .5 * (dx1 + dx2)
                  vx = cfvt(ihx) * dx ** pwvt(ihx)
                  emx = cfmas(ihx) * dx ** pwmas(ihx)

                  gxn = exp(gammln(gnu(ix)))
                  fgamx = (dx / dnx) ** (gnu(ix) - 1) * exp(-dx / dnx)  &
                     / (gxn * dnx)

!------------------------------------------------------------------
! GENERAL COMMENTS FOR GENERALIZING THIS SUBROUTINE

! 1. EMX and EMY are masses of individual X and Y hydrometeors in input bins
! 2. BINT is number of X-Y collisions | per second
!                                     | per bin pair
!                                     | per unit concentration of X and Y (1/m^3)
! 3. In contrast to (deleted subroutine) SXY, this subroutine allows us to
!    distinguish between X and Y during collision.
! 4. Don't know at the time this table is made what fraction of X and Y collide
!    each second or the relative inputs of X and Y to coalesced hydrometeors.
!    This is only knowable during model timestep.
! 5. We CAN get some info, if it would be helpful, on sizes of X and Y in each
!    X-Y bin collision.  For example, would big drizzle and small ice go to
!    category Z at a different rate than big ice and small cloud/drizzle?
! 6. On the other hand, is the info in #5 useful?  Perhaps the collection
!    subroutines provide all the info needed (bulk mass, bulk number, and
!    hence, mean hydrometeor mass) colliding each timestep) in order to decide
!    how much to send to category Z.
! 7. FOR DRIZZLE ACTIVATED:
!      For CLOUD-CLOUD, determine how much goes to DRIZZLE.
!         Use ipairc(8,1) for change in DRIZZLE number?
!      for DRIZZLE-DRIZZLE, determine how much goes to RAIN.
!         Use ipairc(2,8) for change in RAIN number?
!    FOR NO DRIZZLE:
!      For CLOUD-CLOUD, determine how much goes to RAIN.
!         Use ipairc(2,1) for change in RAIN number?
!------------------------------------------------------------------

! BINT is (integrand * del_dx * del_dy)

                  bint = (dx + dy) ** 2 * abs(vx - vy) * fgamx * fgamy  &
                     * efc(ihx,ihy,dx,dy) * (dy2 - dy1) * (dx2 - dx1)

! Every collision gets counted for loss of number for category X

                  sum_num = sum_num + bint

! CONDITIONAL SUMMATION FOR CLOUD-CLOUD AND FOR DRIZZLE-DRIZZLE COLLISIONS
! BASED ON MASS CUTOFF THRESHOLDS.  
!
! The assigned thresholds influence the rate at which mass and number are
! transferred from smaller to larger droplet species when collisions occur
! within the smaller species.  These thresholds should be considered subject
! to modification as model performance is evaluated (11/22/2009).
!
! For comparison, the mass cutoff thresholds used in subroutine SXY of the
! CSU RAMS microphysics model are given.  With the drizzle category activated,
! bins 1-13, 14-17, and 18-36 represent cloud, drizzle, and rain, respectively.
! With the drizzle category NOT activated, bins 1-15 and 16-36 represent
! cloud and rain, respectively.  The diameters and masses of the bins on
! either side of these thresholds are:
!
!          DIAMETER(m)   MASS(kg)
!--------------------------------
! BIN 13     50.0e-6      65.e-12
! BIN 14     63.0e-6     131.e-12
! BIN 15     79.4e-6     262.e-12
! BIN 16    100.0e-6     524.e-12
! BIN 17    126.0e-6    1047.e-12
! BIN 18    158.7e-6    2094.e-12

                  if (ipc == 1 .and. jnmb(8) == 0) then

! Case 1: Cloud-cloud where drizzle is NOT activated in this model run
!         (sum_ymass is actually a number concentration here, not a mass.)

                     if (emx + emy > 400.e-12) then

                        sum_xmass = sum_xmass + bint * emx
                        sum_ymass = sum_ymass + bint

                     endif

                  elseif (ipc == 1 .and. jnmb(8) > 0) then

! Case 2: Cloud-cloud where drizzle IS activated in this model run
!         (sum_ymass is actually a number concentration here, not a mass.)

                     if (emx + emy > 100.e-12) then

                        sum_xmass = sum_xmass + bint * emx
                        sum_ymass = sum_ymass + bint

                     endif

                  elseif (ipc == 62) then

! Case 3: Drizzle-drizzle
!         (sum_ymass is actually a number concentration here, not a mass.)

                     if (emx + emy > 1500.e-12) then

                        sum_xmass = sum_xmass + bint * emx
                        sum_ymass = sum_ymass + bint

                     endif

                  else

! Case 4: All other collisions - standard summation

                     sum_xmass = sum_xmass + bint * emx
                     sum_ymass = sum_ymass + bint * emy

                  endif

               enddo
            enddo

! sum_num and sum_xmass are the definite integral sums of number and mass
! for the current X and Y mean-mass diameter pair.  Enter these in tables.

!no log10   coltabc(iembx,iemby,ipc) = sum_num
            coltabc(iembx,iemby,ipc) = -log10(max(1.e-30,sum_num))

            if (iprxy > 0) then
!no log10      coltabr(iembx,iemby,iprxy) = sum_xmass
               coltabr(iembx,iemby,iprxy) = -log10(max(1.e-30,sum_xmass))
            endif

            if (ipryx > 0) then
!no log10      coltabr(iemby,iembx,ipryx) = sum_ymass
               coltabr(iemby,iembx,ipryx) = -log10(max(1.e-30,sum_ymass))
            endif

         enddo
      enddo
   enddo
enddo

return
end subroutine mkcoltb_brute

!===============================================================================

real function efc(ihx,ihy,dx,dy)

! Evaluate hydrometeor collision efficiency

use micro_coms, only: nhcat, ipairk, nefcx, nefcy, defcx, defcy, efctab

implicit none

integer, intent(in) :: ihx,ihy ! colliding species X and Y
real, intent(in) :: dx,dy      ! diameters of X and Y [m]

integer :: knum ! collision efficiency table number
integer :: ix,iy
integer :: iefcx,iefcy

real :: dx0, dy0
real :: hx,wtx1,wtx2,wty1,wty2

knum = ipairk(ihx,ihy)

! For most collision pairs, set collision effiency to 1.0 and return

if (knum == 10) then
   efc = 1.0
   return
endif

! Collisions between only Cloud, Drizzle, and/or Rain
! (Special steps are taken because collision efficiency table is
! triangular over the two diameters.)

if (knum == 1) then

! Require that dy0 be the larger of the input diameters

   dy0 = max(dx,dy)
   dx0 = min(dx,dy)

! Limit dy0 to maximum diameter in table

   dy0 = min(dy0,.9999 * defcy(nefcy(knum),knum))

! Limit dx0 to a maximum of dy0 (in case dy0 was limited)

   dx0 = min(dx0,.9999 * dy0)
   
! Re-define dx0 as ratio of dx0 to dy0 since for knum = 1, X-dimension of
! efctab is this ratio

   dx0 = dx0 / dy0   

else

! Collisions between Ice and small droplets (Cloud or Drizzle)

! If dy is too small, set efc to zero and return

   if (dy < defcy(1,knum)) then
      efc = 0.
      return
   endif

! If dx is too small, set efc to zero and return

   if (dx < defcx(1,knum)) then
      efc = 0.
      return
   endif

! Limit dy0 to maximum diameter in table

   dy0 = min(dy,.9999 * defcy(nefcy(knum),knum))

! Limit dx0 to maximum diameter in table

   dx0 = min(dx,.9999 * defcx(nefcx(knum),knum))

endif

! Determine interpolation point in X-dimension

iefcx = 1
do while (dx0 > defcx(iefcx+1,knum))
   iefcx = iefcx + 1
enddo

wtx2 = (dx0 - defcx(iefcx,knum)) / (defcx(iefcx + 1,knum) - defcx(iefcx,knum))
wtx1 = 1. - wtx2

! Determine interpolation point in Y-dimension

iefcy = 1
do while (dy0 > defcy(iefcy+1,knum))
   iefcy = iefcy + 1
enddo

wty2 = (dy0 - defcy(iefcy,knum)) / (defcy(iefcy + 1,knum) - defcy(iefcy,knum))
wty1 = 1. - wty2

! Interpolate from table

efc = wtx1 * wty1 * efctab(iefcx  ,iefcy  ,knum) &
    + wtx2 * wty1 * efctab(iefcx+1,iefcy  ,knum) &
    + wtx1 * wty2 * efctab(iefcx  ,iefcy+1,knum) &
    + wtx2 * wty2 * efctab(iefcx+1,iefcy+1,knum)

return
end function efc

!===============================================================================

subroutine tabhab()

use micro_coms, only: jhabtab
use misc_coms,  only: io6

implicit none

integer, parameter :: nhab=0
integer :: it,is

if (nhab ==  0) write(io6,*) 'VARIABLE HABIT PREDICTION'
if (nhab ==  3) write(io6,*) 'ASSUMED HABIT IS COLUMNS'
if (nhab ==  8) write(io6,*) 'ASSUMED HABIT IS HEX PLATES'
if (nhab ==  9) write(io6,*) 'ASSUMED HABIT IS DENDRITES'
if (nhab == 10) write(io6,*) 'ASSUMED HABIT IS NEEDLES'
if (nhab == 11) write(io6,*) 'ASSUMED HABIT IS ROSETTES'
!c    if (nhab .eq.  x) write(io6,*) 'ASSUMED HABIT IS SPHERES'

! nt is temp, ns = satur (liq)

do it = 1,31
   do is = 1,100
      if (nhab == 0) then
         if (it >= 0 .and. it <= 2) then
            if (is <= 95) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            endif
         else if(it > 2 .and. it <= 4) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            endif
         else if(it > 4 .and. it <= 6) then
            if (is < 85) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         else if(it > 6 .and. it <= 9) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         else if(it > 9 .and. it <= 22) then
            if (is < 90) then
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         elseif(it > 22 .and. it <= 30) then
            if (is < 80) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         elseif(it > 30) then
            if (is < 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 12
               jhabtab(it,is,2) = 16
            endif
         endif
      else
         jhabtab(it,is,1) = nhab
         jhabtab(it,is,2) = nhab + 4
         if (nhab == 3) jhabtab(it,is,2) = 4
      endif
   enddo
enddo
return
end subroutine tabhab
