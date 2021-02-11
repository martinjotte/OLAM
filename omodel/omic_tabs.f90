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
subroutine haznuc()

use micro_coms,  only: nthz, dthz, nrhhz, drhhz, frachz
use consts_coms, only: r8

implicit none

integer  :: ithz,irhhz,k
real     :: rhhz,c2hz,c3hz,bhz,thz
real(r8) :: dum, dm, sum, y, dccn

real, external :: gammln

!  Haze nucleation table

real,     parameter :: denccn = 1.769       !Density of ammonium sulfate
real,     parameter :: gnuccn = 2.          !Saleeby(02-21-2007) Originally = 1.
real,     parameter :: dnccn  = 2. * .04e-6 !units are [m]; Saleeby(02-21-2007) Originally = .075e-4 [cm] = .075e-6 [m]
real(r8), parameter :: ddccn  = .005e-4     !Increment [cm] of CCN diameter in numerical integration
                                            !   over gamma size distribution
real,     parameter :: c1hz   = (3.14159 * denccn / 6.) ** (-0.33333333)

do ithz = 1,nthz
   thz = -60. + dthz * real(ithz - 1)
   do irhhz = 1,nrhhz
      rhhz = 0.82 + drhhz * real(irhhz - 1)
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

end subroutine haznuc

!===============================================================================

subroutine homfrzcl(dtl,mrl)

use micro_coms,  only: ntc, dtc, ndnc, ddnc, fracc, gnu
use consts_coms, only: r8

implicit none

integer,  intent(in) :: mrl
real(r8), intent(in) :: dtl

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

end subroutine homfrzcl

!===============================================================================

subroutine tabmelt()

use micro_coms, only: nhcat, lcat_lhcat, gnu, rmlttab, enmlttab, ndns,  &
                      shedtab, cfmas, pwmas, cfvt, pwvt, shapefac, ninc, ncat

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

end subroutine tabmelt

!===============================================================================

subroutine mkcoltb_brute()

use micro_coms, only: nhcat, lcat_lhcat, nembc, gnu, &
                      pwmas, pwmasi, dmb0, dmb1, cfvt, pwvt, cfmas, &
                      ipair, coltabc, coltabx, coltaby, driz_gammq
use consts_coms, only: r8

implicit none

integer, parameter :: ndx=50,ndy=50

integer :: ihx,ix,ihy,iy,iemby,iembx,idx,idy
integer :: ipc,ipx,ipy,ipx2,ipy2

real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gxn,gyn,gym
real :: dnminy,dnmaxy,dny,dnx,bint
real :: sum_num, sum_xmass, sum_ymass, sum_xmass2, sum_ymass2
real :: vx,vy,dx,dy
real :: emx,dx1,dx2,emy,dy1,dy2,dyhi,dylo

real :: emby, dmby, embx, dmbx, sumddx, sumddy

real, external :: gammln, efc

real(r8) :: fgamx, fgamy, dxodnx, dyodny, dxodnxeg1, dyodnyeg1, expdxodnx, expdyodny

! Initialize drizzle size incomplete gamma function

driz_gammq(:) = 0.

! Loop over colliding category X

do ihx = 1,nhcat
   ix = lcat_lhcat(ihx)

   gxm = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
   dnminx = dmb0(ix) * gxm ** pwmasi(ihx)
   dnmaxx = dmb1(ix) * gxm ** pwmasi(ihx)

! Loop over colliding category Y

   do ihy = 1,nhcat

! Get number-concentration collection table number for this colliding pair (X,Y)

      ipc  = ipair(ihx,ihy,1)
      ipx2 = ipair(ihx,ihy,2)
      ipy2 = ipair(ihx,ihy,3)
      ipx  = ipair(ihx,ihy,4)
      ipy  = ipair(ihx,ihy,5)

! If collection table number is zero, this interaction is not considered.
! Bob (11/11/2009): We only need to check ipairc because ipairrx
! and ipairry are now filled together.

      if (ipc == 0) cycle

print*, 'ihx,ihy,ipc ',ihx,ihy,ipc

      iy = lcat_lhcat(ihy)

      gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
      dnminy = dmb0(iy) * gym ** pwmasi(ihy)
      dnmaxy = dmb1(iy) * gym ** pwmasi(ihy)

! Loop over all mean-mass table values for category Y

      do iemby = 1,nembc
         dmby = dmb0(iy) * (dmb1(iy) / dmb0(iy)) ** (real(iemby-1) / real(nembc-1))
         emby = cfmas(iy) * dmby ** pwmas(iy)
         dny = dmby * gym ** pwmasi(ihy)

! Changed Sept 2018: Integrate over narrower diameter range, and base range
! on dmb rather than dn.  Narrower range justified in order to avoid tails.

         dylo = .03 * dmby
         dyhi = 3. * dmby

! Loop over all mean-mass table values for category X

         do iembx = 1,nembc
            dmbx = dmb0(ix) * (dmb1(ix) / dmb0(ix)) ** (real(iembx-1) / real(nembc-1))
            dnx = dmbx * gxm ** pwmasi(ihx)

            dxlo = .03 * dmbx
            dxhi = 3. * dmbx

            sum_num    = 0. ! Initialize integral number sum to zero
            sum_xmass  = 0. ! Initialize integral xmass  sum to zero
            sum_ymass  = 0. ! Initialize integral ymass  sum to zero
            sum_xmass2 = 0. ! Initialize integral xmass2 sum to zero
            sum_ymass2 = 0. ! Initialize integral ymass2 sum to zero

! Loop over spectrum of Y diameters for current mean-mass Y value

            do idy = 1,ndy
               dy1 = dylo * (dyhi / dylo) ** (real(idy-1) / real(ndy))
               dy2 = dylo * (dyhi / dylo) ** (real(idy) / real(ndy))
               dy = .5 * (dy1 + dy2)
               vy = cfvt(ihy) * dy ** pwvt(ihy)

               ! FOR CLOUD, RAIN, AND DRIZZLE, REPLACE FALL SPEED POWER LAW
               ! WITH TABULATED FALL SPEEDS

               if (ihy == 1 .or. ihy == 2 .or. ihy == 8) call vterm_liq(dy,vy)

               emy = cfmas(ihy) * dy ** pwmas(ihy)

               gyn = exp(gammln(gnu(iy)))

               dyodny = real(dy / dny,r8)
               dyodnyeg1 = dyodny ** real(gnu(iy) - 1.0,r8)
               expdyodny = real(exp(-dyodny),r8)

               fgamy = dyodnyeg1 * expdyodny / (gyn * dny)

   if (idy == 1) then
      sumddy = 0.
   endif
   sumddy = sumddy + fgamy*(dy2-dy1)

! Loop over spectrum of X diameters for current mean-mass X value

               do idx = 1,ndx
                  dx1 = dxlo * (dxhi / dxlo) ** (real(idx-1) / real(ndx))
                  dx2 = dxlo * (dxhi / dxlo) ** (real(idx) / real(ndx))
                  dx = .5 * (dx1 + dx2)
                  vx = cfvt(ihx) * dx ** pwvt(ihx)

                  ! FOR CLOUD, RAIN, AND DRIZZLE, REPLACE FALL SPEED POWER LAW
                  ! WITH TABULATED FALL SPEEDS

                  if (ihx == 1 .or. ihx == 2 .or. ihx == 8) call vterm_liq(dx,vx)

                  emx = cfmas(ihx) * dx ** pwmas(ihx)

                  gxn = exp(gammln(gnu(ix)))

                  dxodnx = real(dx / dnx,r8)
                  dxodnxeg1 = dxodnx ** real(gnu(ix) - 1.0,r8)
                  expdxodnx = real(exp(-dxodnx),r8)

                  fgamx = dxodnxeg1 * expdxodnx / (gxn * dnx)

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
!      for DRIZZLE-DRIZZLE, determine how much goes to RAIN.
!    FOR NO DRIZZLE:
!      For CLOUD-CLOUD, determine how much goes to RAIN.
!------------------------------------------------------------------

! BINT is (integrand * del_dx * del_dy)

                  bint = (dx + dy) ** 2 * abs(vx - vy) * fgamx * fgamy  &
                     * efc(ihx,ihy,dx,dy) * (dy2 - dy1) * (dx2 - dx1)


! Every collision gets counted for loss of number for category X

                  sum_num = sum_num + bint

!if (ihx == 1 .and. ihy == 2 .and. &
!   (iembx == 1 .or. iembx == nembc) .and. &
!   (iemby == 1 .or. iemby == nembc)) then

!   if (idx == 1) then
!      sumddx = 0.

!      print*, ' '
!      write(6,'(120a)') '        idx  idy     dx          dy       dxodnx   dyodny', &
!           '     fgamx       fgamy     fgamx*ddy   fgamy*ddy', &
!           '     sumddx      sumddy       bint      sum_num'
!      print*, ' '
!   endif

!   sumddx = sumddx + fgamx*(dx2-dx1)
!   write(6,'(a,2i5,2f12.8,2f9.3,8e12.3)') 'omt0 ', &
!      idx,idy,dx,dy,dxodnx,dyodny,fgamx,fgamy, &
!      fgamx*(dx2-dx1),fgamy*(dy2-dy1),sumddx,sumddy,bint,sum_num
!endif

! CONDITIONAL SUMMATION FOR CLOUD-CLOUD, DRIZZLE-DRIZZLE, AND CLOUD-DRIZZLE
! COLLISIONS BASED ON MASS CUTOFF THRESHOLDS (for cases where ipc = 1, 61, or 62)  
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
! SXY      DIAMETER(m)   MASS(kg)
!--------------------------------
! BIN 10     25.0e-6       8.e-12
! BIN 11     31.5e-6      16.e-12
! BIN 12     39.7e-6      33.e-12
! BIN 13     50.0e-6      65.e-12
! BIN 14     63.0e-6     131.e-12
! BIN 15     79.4e-6     262.e-12
! BIN 16    100.0e-6     524.e-12
! BIN 17    126.0e-6    1047.e-12
! BIN 18    158.7e-6    2094.e-12
! BIN 19    200.0e-6    4188.e-12
! BIN 20    252.0e-6    8376.e-12
! BIN 21    320.0e-6     16.7e-9
! BIN 22    400.0e-6     33.4e-9
! BIN 23    500.0e-6     66.8e-9
! BIN 24    630.0e-6    133.6e-9

! The following 4 size thresholds for (emx + emy) are related to embxz in
! subroutines col1188 and col1882, but do not necessarily need to have the
! same values.

                  if (ipc == 1) then      ! Cloud-cloud interaction table

                     if (emx + emy > 33.e-12) then
                        sum_xmass = sum_xmass + bint * emx
                     endif

                  elseif (ipc == 62) then ! Drizzle-drizzle interaction table

                     if (emx + emy > 33.e-9) then
                        sum_xmass = sum_xmass + bint * emx
                     endif

                  elseif (ipc == 61) then ! Cloud-drizzle interaction table

                     if (emx + emy > 33.e-9) then  ! For transfer to rain
          !!!        if (emx + emy > 4.   ) then  ! For transfer to rain

                        sum_xmass2 = sum_xmass2 + bint * emx
                        sum_ymass2 = sum_ymass2 + bint * emy

                        if (iembx == 1 .and. idx == 1) then
                           driz_gammq(iemby) = driz_gammq(iemby) + fgamy * emy * (dy2 - dy1)
                        endif

                     else                         ! For transfer to drizzle

                        sum_xmass = sum_xmass + bint * emx
                        sum_ymass = sum_ymass + bint * emy

                     endif

                  else                    ! All other collisions - standard summation

                     sum_xmass = sum_xmass + bint * emx
                     sum_ymass = sum_ymass + bint * emy

                  endif

               enddo ! idx
            enddo ! idy

! sum_num, sum_xmass, sum_ymass, sum_xmass2, and sum_ymass2 are the definite
! integral sums of number and mass for the current X and Y mean-mass diameter
! pair.  Enter these in tables.

                           coltabc(iembx,iemby,ipc ) = log(max(1.e-32,sum_num))
             if (ipx  > 0) coltabx(iembx,iemby,ipx ) = log(max(1.e-32,sum_xmass))
             if (ipx2 > 0) coltabx(iembx,iemby,ipx2) = log(max(1.e-32,sum_xmass2))
             if (ipy  > 0) coltaby(iemby,iembx,ipy ) = log(max(1.e-32,sum_ymass))
             if (ipy2 > 0) coltaby(iemby,iembx,ipy2) = log(max(1.e-32,sum_ymass2))

         enddo  ! iembx
         
         if (ipc == 61) then ! Cloud-drizzle interaction table
            driz_gammq(iemby) = driz_gammq(iemby) / emby
         endif

      enddo  ! iemby

print*, ' '
print*, 'ihx,ihy ',ihx,ihy
print*, ' '
write(6,'(100a)') '             1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20'    
print*, ' '

do iemby = 1,20
   write(6,'(a,i5,20f5.1)') 'ipc  ',iemby,(coltabc(iembx,iemby,ipc ),iembx=1,20)
enddo

if (ipx > 0) then
   print*, ' '
   write(6,'(100a)') '             1    2    3    4    5   6    7    8    9   10   11   12   13   14   15   16   17   18   19   20'    
   print*, ' '

   do iemby = 1,20
      write(6,'(a,i5,20f5.1)') 'ipx  ',iemby,(coltabx(iembx,iemby,ipx ) ,iembx=1,20)
   enddo
endif

if (ipy > 0) then
   print*, ' '
   write(6,'(100a)') '             1    2    3    4   5    6    7   8    9   10   11   12   13   14   15   16   17   18   19   20'    
   print*, ' '

   do iemby = 1,20
      write(6,'(a,i5,20f5.1)') 'ipy  ',iemby,(coltaby(iemby,iembx,ipy ),iembx=1,20)
   enddo
endif

if (ipx2 > 0) then
   print*, ' '
   write(6,'(100a)') '             1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20'    
   print*, ' '

   do iemby = 1,20
      write(6,'(a,i5,20f5.1)') 'ipx2 ',iemby,(coltabx(iembx,iemby,ipx2) ,iembx=1,20)
   enddo
endif

if (ipy2 > 0) then
   print*, ' '
   write(6,'(100a)') '             1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20'    
   print*, ' '

   do iemby = 1,20
      write(6,'(a,i5,20f5.1)') 'ipy2  ',iemby,(coltaby(iemby,iembx,ipy2) ,iembx=1,20)
   enddo
endif

!if (ipc == 61) then ! Cloud-drizzle interaction table
   do iemby = 1,20
      write(6,'(a,i5,f10.6)') 'driz_gammq ',iemby,driz_gammq(iemby)
   enddo
!endif

   enddo  ! ihy
enddo  ! ihx

do iemby = 1,20
   write(6,'(a,i5,f10.6)') 'driz_gammq2 ',iemby,driz_gammq(iemby)
enddo

end subroutine mkcoltb_brute

!===============================================================================

subroutine vterm_liq(d,v)

! Based partly on Gunn & Kinzer (1949)

implicit none

real, intent(in)  :: d  ! droplet diameter (m)
real, intent(out) :: v  ! droplet terminal velocity (m/s)

integer :: i

real, parameter :: d0(10) = (/.000080, .000100, .000150,  .000200,  .000300, .000500, .000900, .001800, .003200, .005800/)
real, parameter :: v0(10) = (/   .192,     .27,     .48,      .72,      1.2,     2.0,     3.7,     6.1,     8.3,     9.2/)

if (d <= d0(1)) then

   v = .3e8 * d**2

elseif (d >= d0(10)) then

   v = 9.2

else

   i = 1

   do while(d > d0(i+1))
      i = i + 1
   enddo

   v = v0(i) + (v0(i+1) - v0(i)) * (d - d0(i)) / (d0(i+1) - d0(i))

endif

end subroutine vterm_liq

!===============================================================================

real function efc(ihx,ihy,dx,dy)

! Evaluate hydrometeor collision efficiency

use micro_coms, only: nhcat, ipair, nefcx, nefcy, defcx, defcy, efctab

implicit none

integer, intent(in) :: ihx,ihy ! colliding species X and Y
real, intent(in) :: dx,dy      ! diameters of X and Y [m]

integer :: knum ! collision efficiency table number
integer :: iefcx,iefcy

real :: dx0, dy0
real :: wtx1,wtx2,wty1,wty2

knum = ipair(ihx,ihy,6)

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

end subroutine tabhab

!===============================================================================

subroutine mksedim_tab()

use micro_coms, only: sedtime0, sedtime1, mza0, vf_max,  &
                      nhcat, lcat_lhcat, dispemb0, cfvt, emb0, cfmasi,  &
                      pwvt, pwmasi, dispemb0i, dispemb1, ch2, nembfall, gnu,  &
                      pwmas, pcpfillc, pcpfillr, emb1, kfall, ch3
use mem_grid,   only: zm, dzt, zfacm2, zfacim2
use misc_coms,  only: dtlm

implicit none

!integer, parameter :: nbin=50
integer, parameter :: nbin=200
integer :: iembs,lcat,lhcat,k,ibin,kk,jbin

real :: dmbodn,diam0,diam1,diam,areascale, &
        disp,zbotnew,fallin,dispemb,dispmax

real, external :: gammln,gammp

real(8) :: fac1, fac2, fac3, fac4, sumc, sumr, sumv

real :: cbin   (nbin)
real :: rbin   (nbin)
real :: reldisp(nbin)

allocate (pcpfillc(mza0,mza0,nembfall,nhcat))
allocate (pcpfillr(mza0,mza0,nembfall,nhcat))
allocate (kfall   (         mza0,nembfall,nhcat))

pcpfillc = 0.0
pcpfillr = 0.0
kfall    = 1

! To cover the possible range over all timesteps, define sedtime0 and sedtime1
! here as 0.1 seconds and 3000 seconds.  The former is supposed to be
! less than 0.7 of the shortest timestep on any grid (sqrt(rhoi) never exceeds
! 0.7) and the latter is the longest timestep expected to ever be used (600
! seconds) times a factor of 5 for the largest value of sqrt(rhoi).

sedtime0 = dtlm(1) * 0.7
sedtime1 = dtlm(1) * 4.

dispmax  = 1500.
vf_max   = 12.0

! Loop over hydrometeor categories

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   dispemb0(lhcat) = sedtime0 * cfvt(lhcat)  &
      * (emb0(lcat) * cfmasi(lhcat)) ** ch3(lhcat)

   dispemb0i(lhcat) = 1. / dispemb0(lhcat)

   dispemb1(lhcat) = sedtime1 * min(vf_max, cfvt(lhcat) &
      * (emb1(lcat) * cfmasi(lhcat)) ** ch3(lhcat) )

!Bob (10/24/00):  Limit dispemb1 to a maximum of dispmax

!   if (dispemb1(lhcat) > dispmax) dispemb1(lhcat) = dispmax

   ch2(lhcat) = real(nembfall-1) / log(dispemb1(lhcat) * dispemb0i(lhcat))


!!   ri0(lhcat) = 1.0 + ch2(lhcat) * log(dispemb0i(lhcat))
!!   ri1(lhcat) = ch2(lhcat) * log(vf_max*real(dtlm(1)))
!!   r12(lhcat) = ch2(lhcat) * log(cfvt(lhcat)*cfmasi(lhcat)**ch3(lhcat))

! Loop over bins, filling them with fractional number, fractional mass,
! and displacement quotient relative to emb.

   dmbodn = (exp(gammln(gnu(lcat) + pwmas(lhcat))  &
      - gammln(gnu(lcat)))) ** pwmasi(lhcat)

   diam0 = 0.06 * dmbodn
   diam1 = 4.0 * dmbodn

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

   sumv = 0.

   do jbin = 1, nbin
      cbin(jbin) = cbin(jbin) / sumc
      rbin(jbin) = rbin(jbin) / sumr
      sumv = sumv + rbin(jbin) * reldisp(jbin)
   enddo

   do jbin = 1, nbin
      reldisp(jbin) = reldisp(jbin) / sumv
   enddo

! Loop over displacement distance for size emb.

   do iembs = 1, nembfall

      dispemb = dispemb0(lhcat)  &
         * (dispemb1(lhcat) * dispemb0i(lhcat))  &
         ** (real(iembs-1) / real(nembfall-1))

! Loop over vertical grid index.

      do k = 2, mza0

! Loop over bins

         do ibin = 1,nbin
            disp = min( dispemb * reldisp(ibin), dispmax)

            ! New bottom height of source-cell precipitation after fall for 1 time step
            zbotnew = zm(k-1) - disp

 ! Loop over grid cells that a parcel falls into

            do kk = k-1, 1, -1

               if (zbotnew > zm(kk)) exit

               ! Horizontal area scale factor for current w(kk) level
               ! areascale = zfacim2(k-1) * zfacm2(kk)

               ! Depth of source grid cell precipitation that crosses w(kk) level
               fallin = min(dzt(k), zm(kk) - zbotnew)

               pcpfillc(kk,k,iembs,lhcat) = pcpfillc(kk,k,iembs,lhcat)  &
                  + fallin * cbin(ibin) !* areascale

               pcpfillr(kk,k,iembs,lhcat) = pcpfillr(kk,k,iembs,lhcat)  &
                  + fallin * rbin(ibin) !* areascale
            enddo

         enddo

         do kk = k-2, 1, -1
            if (pcpfillr(kk,k,iembs,lhcat) < .005 * dzt(k)) exit
         enddo

         kfall(k,iembs,lhcat) = k - kk - 1

      enddo
   enddo
enddo

end subroutine mksedim_tab
