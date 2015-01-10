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
subroutine micinit_fields()

use mem_basic,   only: rho

use mem_micro,   only: sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                       accpd, accpr, accpp, accps, accpa, accpg, accph, &
                       pcprd, pcprr, pcprp, pcprs, pcpra, pcprg, pcprh, &
                       con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h, &
                       con_ccn, con_gccn, con_ifn, q2, q6, q7, &
                       pcpgr, qpcpgr, dpcpgr, cldnum

use misc_coms,   only: io6, runtype
use mem_ijtabs,  only: jtab_w, jtw_init
use mem_grid,    only: mza, glatw, glonw
use consts_coms, only: r8

use micro_coms,  only: level, icloud, idriz, irain, ipris, isnow, iaggr, &
                       igraup, ihail, jnmb, cparm, dparm, pparm

implicit none

integer :: j, iw
integer :: ilat1, ilon1, ilat2, ilon2
real :: qlat, qlon, rlat, rlon, wlat1, wlon1, wlat2, wlon2

! Geographic map of cloud droplet concentration [*1.e7/m^3]

integer :: numcld(37,19)

data ((numcld(ilon1,ilat1),ilon1=1,37),ilat1=19,1,-1)/    &
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !9
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !8
  3, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 3, 3, 3, 3, 3, 8, 8,13,18,18,18,18,18,18,18,13,13,13, 8, 8, 8, 8, 8, 8, 3, & !7
  8, 8, 8, 8, 8,13,13,18,13, 8, 8, 8, 8, 8, 8, 8, 8,13,18,28,38,48,58,53,48,43,33,23,23,18,18,13,13, 8, 8, 8, 8, & !6
  8, 8, 8, 8, 8, 8,13,18,18,23,23,18,18,13, 8, 8,13,13,23,53,58,58,58,48,43,43,38,33,28,28,33,33,28,18,13, 8, 8, & !5
  8, 8, 8, 8, 8, 8,13,18,28,43,48,28,13,13, 8,13,13,18,28,38,43,43,38,28,28,33,33,23,33,58,93,58,38,18,13, 8, 8, & !4
  8, 8, 8, 8, 8,13,13,18,23,33,28,13,13,13, 8,13,13,13,18,23,38,33,23,23,33,33,28,23,33,83,83,23,18,13, 8, 8, 8, & !3
  8, 8, 8, 8, 8,13,13,18,23,13,13, 8, 8, 8, 8,13,13,13,18,23,28,28,33,33,28,23,33,23,33,43,23,13, 8, 8, 8, 8, 8, & !2
  8, 8, 8, 8, 8, 8, 8,13,13,13,13,13,13, 8, 8, 8,13,18,23,23,23,23,18,18,13,13,18,13,13,13,13, 8, 8, 8, 8, 8, 8, & !1
  8, 8, 8, 8, 8, 8,13,13,13,13,13, 8, 8, 8, 8, 8, 8, 8, 8,13,13,13,13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, & !0
  8, 8, 8, 8, 8, 8, 8, 8,13,13,13, 8, 8, 8, 8, 8, 8, 8, 8,13,18,18,13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, & !1
  3, 3, 3, 3, 3, 3, 3, 3, 8,13,13,23,13,13, 8, 3, 8, 8, 8,13,18,13, 8, 8, 3, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, & !2
  3, 3, 3, 3, 3, 3, 3, 3, 3, 8,13,18,13,13, 8, 3, 3, 3, 8, 8,18,18, 8, 8, 3, 3, 3, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, & !3
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 8,13,13, 8, 8, 3, 3, 3, 3, 3, 8, 8, 8, 8, 3, 3, 3, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, & !4
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 8, 8, 8, 8, 8, 8, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !5
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !6
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !7
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, & !8
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3  / !9

!18 17 16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18

! Initialize geographic-based cloud number concentration

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!----------------------------------------------------------------------

   qlon = max(-179.9999,min(179.9999,glonw(iw)))
   qlat = max( -89.9999,min( 89.9999,glatw(iw)))

   rlon = 0.1 * (qlon + 180.) + 1.
   rlat = 0.1 * (qlat +  90.) + 1.

   ilon1 = int(rlon)
   ilat1 = int(rlat)

   wlon2 = rlon - real(ilon1)
   wlat2 = rlat - real(ilat1)

   ilon2 = ilon1 + 1
   ilat2 = ilat1 + 1

   wlon1 = 1. - wlon2
   wlat1 = 1. - wlat2

! Interpolate from 4 surrounding values

   cldnum(iw) = wlon1 * wlat1 * real(numcld(ilon1,ilat1)) &
              + wlon2 * wlat1 * real(numcld(ilon2,ilat1)) &
              + wlon1 * wlat2 * real(numcld(ilon1,ilat2)) &
              + wlon2 * wlat2 * real(numcld(ilon2,ilat2))

   cldnum(iw) = cldnum(iw) * 1.e7

   if (cparm > 1.e6) cldnum(iw) = cparm

enddo

! Initialize 3D and 2D microphysics fields

if (runtype == 'INITIAL') then

!----------------------------------------------------------------------
   do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!----------------------------------------------------------------------

      if (allocated(sh_d))  then
         sh_d(1:mza,iw) = 0.
         accpd(iw) = 0._r8
         pcprd(iw) = 0.
      endif

      if (allocated(sh_r))  then
         sh_r(1:mza,iw) = 0.
         accpr(iw) = 0._r8
         pcprr(iw) = 0.
         q2(1:mza,iw) = 0.
      endif
   
      if (allocated(sh_p))  then
         sh_p(1:mza,iw) = 0.
         accpp(iw) = 0._r8
         pcprp(iw) = 0.
      endif

      if (allocated(sh_s))  then
         sh_s(1:mza,iw) = 0.
         accps(iw) = 0._r8
         pcprs(iw) = 0.
      endif

      if (allocated(sh_a))  then
         sh_a(1:mza,iw) = 0.
         accpa(iw) = 0._r8
         pcpra(iw) = 0.
      endif

      if (allocated(sh_g)) then
         sh_g(1:mza,iw) = 0.
         accpg(iw) = 0._r8
         pcprg(iw) = 0.
         q6(1:mza,iw) = 0.
      endif

      if (allocated(sh_h))  then
         sh_h(1:mza,iw) = 0.
         accph(iw) = 0._r8
         pcprh(iw) = 0.
         q7(1:mza,iw) = 0.
      endif

      if (allocated(con_c)) con_c(1:mza,iw) = 0.
      if (allocated(con_r)) con_r(1:mza,iw) = 0.
      if (allocated(con_p)) con_p(1:mza,iw) = 0.
      if (allocated(con_s)) con_s(1:mza,iw) = 0.
      if (allocated(con_a)) con_a(1:mza,iw) = 0.
      if (allocated(con_g)) con_g(1:mza,iw) = 0.
      if (allocated(con_h)) con_h(1:mza,iw) = 0.
      if (allocated(con_d)) con_d(1:mza,iw) = 0.

! Initialize CCN field if activated
      
      if (allocated(con_ccn)) then
         con_ccn(1:mza,iw) = cparm  ! Default specification

! Steve's example (changed to #/kg)

!        con_ccn(k,iw) = max(100.e6,cparm * (1. - zt(k) / 4000.))

      endif

! Initialize GCCN field if activated

      if (allocated(con_gccn)) then
         con_gccn(1:mza,iw) = dparm  ! Default specification

! Steve's example (changed to #/kg)

!        con_gccn(k,iw) = max(10.,dparm * (1. - zt(k) / 4000.))

      endif

! Initialize IFN field if activated

      if (allocated(con_ifn)) then
         con_ifn (1:mza,iw) = pparm * rho(1:mza,iw) ** 5.4  ! Default specification
      endif

! Initialize accumulated precipitation for surface model

      if (allocated(pcpgr)) then
         pcpgr(iw)  = 0.
         qpcpgr(iw) = 0.
         dpcpgr(iw) = 0.
      endif

   enddo

endif ! runtype == 'INITIAL'

end subroutine micinit_fields

!===============================================================================

subroutine micinit_tabs()

use mem_basic, only: rho

use misc_coms,   only: io6, dtlm, runtype
use mem_ijtabs,  only: jtab_w, mrls, jtw_init
use mem_grid,    only: mza, zm, dzt, dzit
use mem_para,    only: myrank

use micro_coms,  only: level, jnmb, nembc, &
                       cfmas, pwmas, cfvt, pwvt, &
                       npairx, npairy, npairc, coltabx, coltaby, coltabc, &
                       alloc_sedimtab

use hdf5_utils, only: shdf5_irec, shdf5_orec, shdf5_open, shdf5_close

implicit none

  integer :: j,iw,k,mrl,ndims,idims(3)

  character(len=80) :: coltabfile

  logical :: exans

  call micinit_gam()

  if (level < 3) return

  call haznuc()

  call tabmelt()

  call tabhab()

  call alloc_sedimtab(mza)
  call mksedim_tab(mza,zm,dzt,dzit)

  do mrl = 1,mrls
     call homfrzcl(dtlm(mrl),mrl)
  enddo

! Check if collection table file exists

  coltabfile = './COLTABFILE2'

  inquire(file=coltabfile, exist=exans)

  if (exans) then

! Collection table file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening collection table file ', trim(coltabfile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(trim(coltabfile),'R')

     ndims = 3
     idims(1) = nembc
     idims(2) = nembc
     idims(3) = npairc

     call shdf5_irec(ndims, idims, 'COLTABC', rvara=coltabc)

     idims(3) = npairx

     call shdf5_irec(ndims, idims, 'COLTABX', rvara=coltabx)

     idims(3) = npairy

     call shdf5_irec(ndims, idims, 'COLTABY', rvara=coltaby)

! Close the collection table file

     call shdf5_close()

  else

! Collection table file does not exist.

! Make collection table

     call mkcoltb_brute()

! Open, write, and close collection table file.

     if (myrank == 0) then

        call shdf5_open(trim(coltabfile),'W',0)

        ndims = 3
        idims(1) = nembc
        idims(2) = nembc
        idims(3) = npairc

        call shdf5_orec(ndims, idims, 'COLTABC', rvara=coltabc)

        idims(3) = npairx

        call shdf5_orec(ndims, idims, 'COLTABX', rvara=coltabx)

        idims(3) = npairy

        call shdf5_orec(ndims, idims, 'COLTABY', rvara=coltaby)

! Close the collection table file

        call shdf5_close()

     endif
  endif

end subroutine micinit_tabs

!===============================================================================

subroutine micinit_gam()

use micro_coms,  only: icloud, idriz, irain, ipris, isnow, iaggr, igraup, ihail, &
                       parm, nhcat, shapefac, cfmas, pwmas, cfvt, pwvt, &
                       ncat, emb0, emb1, gnu, rxmin, level, sl, sc, sj, &
                       cparm, dparm, rparm, pparm, sparm, aparm, gparm, hparm, &
                       sk, dps, dps2, rictmin, rictmax, nembc, lcat_lhcat, &
                       emb0log, emb1log, emb2, cfmasi, pwmasi, pwen0, &
                       pwemb0, ch3, cdp1, pwvtmasi, jnmb, cfemb0, cfen0, &
                       dnfac, vtfac, frefac1, frefac2, sipfac, cfmasft, &
                       dict, dpsmi, gamm, gamn1, ngam, gam, gaminc, &
                       gamsip13, gamsip24, ddn_ngam, reffcof, dmncof
                      
use consts_coms, only: alvl, alvi, alli
use misc_coms,   only: io6

implicit none

integer :: lhcat,lcat,igam
real :: c1,glg,glg1,glg2,glgm,glgc,glgmv,gym,flngi,dpsi,embsip,dnsip
real :: gammln,gammp,gammq

real :: dstprms(9,nhcat) = reshape( (/ &  ! Carver/Mitchell 1996 power laws
!----------------------------------------------------------------------
!shapefac  cfmas  pwmas   cfvt    pwvt    dmb0     dmb1  gnu  rxmin
!----------------------------------------------------------------------
    .5,     524.,    3.,  3173.,    2.,  2.e-6,  50.e-6,  9., 1.e-12, & !cloud
    .5,     524.,    3.,   144.,  .497,  .1e-3,   5.e-3,  2.,  1.e-9, & !rain
  .179,    110.8,  2.91,  1538.,  1.00, 15.e-6, 125.e-6,  2., 1.e-12, & !pris col
  .179, 2.739e-3,  1.74,   27.7,  .484,  .1e-3,  10.e-3,  2.,  1.e-9, & !snow col
    .5,     .496,   2.4,   16.1,  .416,  .1e-3,  10.e-3,  2.,  1.e-9, & !aggreg
    .5,     157.,    3.,   332.,  .786,  .1e-3,   5.e-3,  2.,  1.e-9, & !graupel
    .5,     471.,    3.,  152.1,  .497,  .8e-3,  10.e-3,  2.,  1.e-9, & !hail
    .5,     524.,    3., 1.26e7,  1.91, 65.e-6, 100.e-6,  4., 1.e-12, & !drizzle
  .429,    .8854,   2.5, 20801., 1.377,    00.,     00., 00.,    00., & !pris hex
 .3183,  .377e-2,    2.,   56.4,  .695,    00.,     00., 00.,    00., & !pris den
 .1803,  1.23e-3,   1.8, 1617.9,  .983,    00.,     00., 00.,    00., & !pris ndl
    .5,    .1001, 2.256,  6239.,  1.24,    00.,     00., 00.,    00., & !pris ros
  .429,    .8854,   2.5,  30.08,  .563,    00.,     00., 00.,    00., & !snow hex
 .3183,  .377e-2,    2.,   3.39,  .302,    00.,     00., 00.,    00., & !snow den
 .1803,  1.23e-3,   1.8,   44.6,  .522,    00.,     00., 00.,    00., & !snow ndl
    .5,    .1001, 2.256,  125.7,  .716,    00.,     00., 00.,    00.  & !snow ros
    /), (/ 9,nhcat /) )

! Initialize arrays based on microphysics namelist parameters

parm(1) = cparm
parm(2) = rparm
parm(3) = pparm
parm(4) = sparm
parm(5) = aparm
parm(6) = gparm
parm(7) = hparm
parm(8) = dparm

if (icloud <= 1) parm(1) = .3e9
if (irain  == 1) parm(2) = .1e-2
if (ipris  == 1) parm(3) = 100. ! obsolete
if (isnow  == 1) parm(4) = .1e-2
if (iaggr  == 1) parm(5) = .1e-2
if (igraup == 1) parm(6) = .1e-2
if (ihail  == 1) parm(7) = .3e-2
if (idriz  == 1) parm(8) = .1e6  !# per kg (mid-range avg from Feingold(99)

! Copy individual arrays from dstprms data table

do lhcat = 1,nhcat
   shapefac(lhcat) = dstprms(1,lhcat)
   cfmas   (lhcat) = dstprms(2,lhcat)
   pwmas   (lhcat) = dstprms(3,lhcat)
   cfvt    (lhcat) = dstprms(4,lhcat)
   pwvt    (lhcat) = dstprms(5,lhcat)
enddo

do lcat = 1,ncat
   emb0 (lcat) = cfmas(lcat) * dstprms(6,lcat) ** pwmas(lcat)
   emb1 (lcat) = cfmas(lcat) * dstprms(7,lcat) ** pwmas(lcat)
   emb2 (lcat) = cfmas(lcat) * parm(lcat) ** pwmas(lcat)
   gnu  (lcat) = dstprms(8,lcat)
   rxmin(lcat) = dstprms(9,lcat)
enddo

if (level < 3) RETURN

! Initialize constants for vapor diffusion

sl(1) = alvl
sl(2) = alvi
sc(1) = 4186.
sc(2) = 2093.  ! 2106 is correct value
sj(1) = 0
sj(2) = 1
sj(3) = 0
sj(4) = 0
sj(5) = 0
sj(6) = 1
sj(7) = 1
sj(8) = 0
sk(1) = alli
sk(2) = 0.

dps = 125.e-6
dps2 = dps ** 2
rictmin = 1.0001
rictmax = 0.9999 * float(nembc)

do lhcat = 1,nhcat
   lcat = lcat_lhcat(lhcat)

   emb0log(lcat) = log(emb0(lcat))
   emb1log(lcat) = log(emb1(lcat))

! Define coefficients [vtfac, frefac1, frefac2] used for terminal velocity
! and Reynolds number

   cfmasi(lhcat)   = 1. / cfmas(lhcat)
   pwmasi(lhcat)   = 1. / pwmas(lhcat)
   pwen0(lhcat)    = 1. / (pwmas(lhcat) + 1.)
   pwemb0(lhcat)   = pwmas(lhcat) / (pwmas(lhcat) + 1.)
   ch3(lhcat)      = pwvt(lhcat) * pwmasi(lhcat)
   cdp1(lhcat)     = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
   pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)

   c1 = 1.5 + .5 * pwvt(lhcat)
   glg = gammln(gnu(lcat))
   glg1 = gammln(gnu(lcat) + 1.)
   glg2 = gammln(gnu(lcat) + 2.)
   glgm = gammln(gnu(lcat) + pwmas(lhcat))
   glgc = gammln(gnu(lcat) + c1)
   glgmv = gammln(gnu(lcat) + pwmas(lhcat) + pwvt(lhcat))

   dnfac(lhcat) = (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

   vtfac(lhcat) = cfvt(lhcat) * exp(glgmv - glgm)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwvt(lhcat) * pwmasi(lhcat))

   frefac1(lhcat) = shapefac(lhcat) * exp(glg1 - glg)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac2(lhcat) = shapefac(lhcat) * 0.229 * sqrt(cfvt(lhcat))  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)  &
      * exp(glgc - glg)

   sipfac(lhcat) = .785 * exp(glg2 - glg)  &
      * (cfmasi(lhcat) * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

   cfmasft(lhcat) = cfmas(lhcat) * exp(glgm - glg)

   dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

   dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))

   reffcof(lhcat) = 0.5 * (gnu(lcat) + 2.) * dnfac(lhcat)

   dmncof(lhcat) = gnu(lcat) * dnfac(lhcat)

   if (lhcat <= 4) gamm(lhcat) = exp(glg)
   if (lhcat <= 4) gamn1(lhcat) = exp(glg1)

enddo

flngi = 1. / real(ngam)

do igam = 1,ngam
   dpsi = dps / (ddn_ngam * real(igam))

! gam1   :  the integral of the pristine distribution from dps to infty
! gam2   :  the integral of the snow dist. from 0 to dps
! gam3   :  values of the exponential exp(-dps/dn)

   gam(igam,1) = gammq(gnu(3) + 1., dpsi)
   gam(igam,2) = gammp(gnu(4) + 1., dpsi)
   gam(igam,3) = exp(-dpsi)

   GAMINC(igam,1) = GAMMQ(GNU(3),dpsi)
   GAMINC(igam,2) = GAMMP(GNU(4),dpsi)

! Secondary ice production arrays

   embsip = emb1(1) * real(igam) * flngi
   dnsip = dnfac(1) * embsip ** pwmasi(1)
   gamsip13(1,igam) = gammp(gnu(1),13.e-6/dnsip)
   gamsip24(1,igam) = gammq(gnu(1),24.e-6/dnsip)

   embsip = emb1(8) * real(igam) * flngi
   dnsip = dnfac(8) * embsip ** pwmasi(8)
   gamsip13(2,igam)= gammp(gnu(8),13.e-6/dnsip)
   gamsip24(2,igam)= gammq(gnu(8),24.e-6/dnsip)
enddo

return
end subroutine micinit_gam

!===============================================================================

subroutine jnmbinit()

use micro_coms, only: level, jnmb, icloud, idriz, irain, ipris, isnow, iaggr, &
                      igraup, ihail
use misc_coms,  only: io6

implicit none

if (level /= 3) then

   if (level <= 1) then
      jnmb(1) = 0
   else
      jnmb(1) = 4
   endif

   jnmb(2) = 0
   jnmb(3) = 0
   jnmb(4) = 0
   jnmb(5) = 0
   jnmb(6) = 0
   jnmb(7) = 0
   jnmb(8) = 0

else

   jnmb(1) = icloud
   jnmb(2) = irain
   jnmb(3) = ipris
   jnmb(4) = isnow
   jnmb(5) = iaggr
   jnmb(6) = igraup
   jnmb(7) = ihail
   jnmb(8) = idriz

   if (icloud == 1) jnmb(1) = 4
   if (irain  == 1) jnmb(2) = 2
   if (ipris  == 1) jnmb(3) = 5
   if (isnow  == 1) jnmb(4) = 2
   if (iaggr  == 1) jnmb(5) = 2
   if (igraup == 1) jnmb(6) = 2
   if (ihail  == 1) jnmb(7) = 2
   if (idriz  == 1) jnmb(8) = 4

   if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or. &
      igraup == 5 .or. ihail == 5) then

      if (irain  >= 1) jnmb(2) = 5
      if (isnow  >= 1) jnmb(4) = 5
      if (iaggr  >= 1) jnmb(5) = 5
      if (igraup >= 1) jnmb(6) = 5
      if (ihail  >= 1) jnmb(7) = 5

   endif

endif

return
end subroutine jnmbinit

