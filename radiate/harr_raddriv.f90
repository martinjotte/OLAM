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
subroutine harr_raddriv(iw,ka,nrad,koff,rlong_previous)

use mem_harr, only: mg, mb

use mem_grid, only: mza, zm, zt, glatw, glonw

use mem_basic, only: rho, press, theta, sh_v

use misc_coms, only: io6, iswrtyp, ilwrtyp, time8

use consts_coms, only: rocp, p00i, stefan, cp

use mem_radiate, only: rshort, rlong, fthrd, rlongup, cosz, albedt,   &
     rshort_top,rshortup_top,rlongup_top, fthrd_lw, albedt_beam,   &
     albedt_diffuse,rshort_diffuse,rlong_albedo

use micro_coms, only: ncat

!use radsave

implicit none

integer, intent(in) :: iw
integer, intent(in) :: ka
integer, intent(in) :: nrad
integer, intent(in) :: koff
real, intent(in) :: rlong_previous

! Local arrays

integer :: jhcat(mza,ncat)  ! hydrom category table with ice habits

real :: tairk(mza) ! air temperature [K]
real :: rhov (mza) ! vapor density [kg_vap/m^3]
real :: rx   (mza,ncat)  ! hydrom bulk spec dens [kg_hyd/kg_air]
real :: cx   (mza,ncat)  ! hydrom bulk number [num_hyd/kg_air]
real :: emb  (mza,ncat)  ! hydrom mean particle mass [kg/particle]

real :: rl (nrad) ! vapor density of all radiation levels (kg/m^3)
real :: dzl(nrad) ! delta-z (m) of all radiation levels
real :: dl (nrad) ! air density of all radiation levels (kg/m^3)
real :: pl (nrad) ! pressure (Pa)
real :: o3l(nrad) ! stores the calculated ozone profile (g/m^3)
real :: vp (nrad) ! vapor pressure (Pa)

real :: u(nrad,3) ! path-length for gases (H_2O, CO_2, O_3)  (Pa)

real :: tp  (nrad,mb) ! optical depth of hydrometeors (m^-1)
real :: omgp(nrad,mb) ! Single scatter albedo of hydrometeors
real :: gp  (nrad,mb) ! Asymmetry factor of hydrometeors

real :: zml  (nrad) ! heights of W points of all radiation levels (m)
real :: ztl  (nrad) ! heights of T points of all radiation levels (m)
real :: tl   (nrad) ! temperature (K)
real :: flxus(nrad) ! Total upwelling s/w flux (W/m^2)
real :: flxds(nrad) ! Total downwelling s/w flux (W/m^2)
real :: flxul(nrad) ! Total upwelling l/w flux (W/m^2)
real :: flxdl(nrad) ! Total downwelling l/w flux (W/m^2)

real :: fu   (nrad,6) ! upwelling fluxes for pseudo-bands (W/m^2)
real :: fd   (nrad,6) ! downwelling fluxes for pseudo-bands (W/m^2)

! Set activation flags for gases of importance:  
! Flag = 1: gas active;    Flag = 0: gas not active

integer, save :: ngass(mg)=(/1, 1, 1/)  ! Flags for (H2O, CO2, O3) for shortwave
integer, save :: ngast(mg)=(/1, 1, 1/)  ! Flags for (H2O, CO2, O3) for longwave

integer k,ib,ig,kk,ik,krad,mcat

real rmix
real dzl9,rvk0,rvk1
real :: flx_diff

! Copy surface and vertical-column values from model to radiation memory space
! In this loop, (k-koff) ranges from 2 to mza + 1 - lpw(iw)

do k = ka,mza-1
   krad = k - koff

   tairk(k) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
   rhov(k) = max(0.,sh_v(k,iw)) * rho(k,iw)

   dl(krad)  = rho  (k,iw)
   pl(krad)  = press(k,iw)
   tl(krad)  = tairk(k)
   rl(krad)  = rhov (k)
   zml(krad) = zm   (k)
   ztl(krad) = zt   (k)
enddo
 
! Fill surface values

zml(1) = zm(1+koff)
ztl(1) = zt(1+koff)
pl (1) = pl(2) + (zml(1) - ztl(2)) / (ztl(2) - ztl(3)) * (pl(2) - pl(3))

call rad_mclat(iw,nrad,koff,glatw(iw),dl,pl,rl,tl,o3l,zml,ztl,dzl)

tl(1) = sqrt(sqrt((rlongup(iw) + rlong_previous * rlong_albedo(iw))/ stefan))
dl(1) = dl(2)
rl(1) = rl(2)

! zero out scratch arrays

u (:,:) = 0.0
fu(:,:) = 0.0
fd(:,:) = 0.0

! Fill arrays rx, cx, and emb with hydrometeor properties

call cloudprep_rad(iw,ka,mcat,jhcat,tairk,rhov,rx,cx,emb)

! Fill hydrometeor optical property arrays [tp, omgp, gp]

call cloud_opt(iw,ka,nrad,koff,mcat,jhcat,cx,emb,tp,omgp,gp,real(time8))

! Get the path lengths for the various gases...

call path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl)

do k = 1,nrad
   if (rl(k) <   0. .or.  &
       dl(k) <   0. .or.  &
       pl(k) <   0. .or.  &
      o3l(k) <   0. .or.  &
       tl(k) < 100.) then
      write(io6,*) 'Temperature too low or negative value of'
      write(io6,*) 'density, vapor, pressure, or ozone'
      write(io6,*) 'before calling Harrington radiation'
      write(io6,*) 'at k,iw = ',k,iw,' lat/lon = ',glatw(iw),glonw(iw)
      write(io6,*) 'stopping model'
      write(io6,*) 'rad: k, rl(k), dl(k), pl(k), o3l(k), tl(k)'
      do kk=1,nrad
         write(io6,'(i3,5g15.6)') kk, rl(kk), dl(kk), pl(kk), o3l(kk), tl(kk)
      enddo
      stop 'stop: radiation call'
   endif

enddo

! Harrington shortwave scheme (valid only if cosz > .03)

if (iswrtyp == 3 .and. cosz(iw) > 0.03) then

   flxus(:) = 0.0
   flxds(:) = 0.0

   call harr_swrad(nrad,iw,albedt_beam(iw),albedt_diffuse(iw),cosz(iw),  &
        real(time8),u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxus,flxds,ngass,flx_diff)

   rshort(iw) = flxds(1)
   rshort_diffuse(iw) = flx_diff
   rshort_top(iw) = flxds(nrad)
   rshortup_top(iw) = flxus(nrad)
   albedt(iw) = albedt_beam(iw) + rshort_diffuse(iw) / rshort(iw) *  &
        (albedt_diffuse(iw) - albedt_beam(iw))

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw)  &
         + (flxds(krad) - flxds(krad-1) + flxus(krad-1) - flxus(krad))  &
         / (dl(krad) * dzl(krad) * cp)
   enddo

endif

! Harrington longwave scheme

if (ilwrtyp == 3) then

   flxul(:) = 0.0
   flxdl(:) = 0.0

   call harr_lwrad(nrad,u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxul,flxdl,ngast)

   rlong(iw) = flxdl(1)
   rlongup(iw) = flxul(1)
   rlongup_top(iw) = flxul(nrad)

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw)  &
         + (flxdl(krad) - flxdl(krad-1) + flxul(krad-1) - flxul(krad))  &
         / (dl(krad) * dzl(krad) * cp)
      fthrd_lw(k,iw) = fthrd_lw(k,iw)  &
         + (flxdl(krad) - flxdl(krad-1) + flxul(krad-1) - flxul(krad))  &
         / (dl(krad) * dzl(krad) * cp)
   enddo

endif

!      if (i == isave .and. j == jsave) then
!         downs = rshort
!         downl = rlong
!         do k = lpw,m1-1
!            kk = k - koff
!            ttendl(k) = (flxdl(kk) - flxdl(kk-1) + flxul(kk-1) - flxul(kk)) &
!               / (dl(kk) * dzl(kk) * cp)
!            ttends(k) = (flxds(kk) - flxds(kk-1) + flxus(kk-1) - flxus(kk)) &
!               / (dl(kk) * dzl(kk) * cp)
!         enddo
!         do k = lpw-1,m1-1
!            kk = k - koff
!            sfluxu(k) = flxus(kk)
!            sfluxd(k) = flxds(kk)
!            lfluxu(k) = flxul(kk)
!            lfluxd(k) = flxdl(kk)
!         enddo
!      endif

return
end subroutine harr_raddriv

!===============================================================================

subroutine cloudprep_rad(iw,ka,mcat,jhcat,tairk,rhov,rx,cx,emb)

! This subroutine was developed from parts of subroutine MIC_COPY in
! omic_driv.f90 and subroutines EACH_COLUMN and ENEMB in omic_misc.f90.

! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR, 
! NOT PER M^3 AS IN THE ORIGINAL MICROPHYSICS SUBROUTINES.

use micro_coms, only: ncat, jnmb, rxmin, jhabtab, level,  &
                      emb0, emb1, emb2, parm

use mem_micro,  only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, sh_d,  &
                      con_c, con_r, con_p, con_s, con_a, con_g, con_h, con_d

use mem_grid,   only: mza

implicit none

integer, intent(in) :: iw
integer, intent(in) :: ka

integer, intent(out) :: mcat  ! # of active hydrom categories (0,1, or 7)
integer, intent(out) :: jhcat(mza,ncat)  ! hydrom category table with ice habits

real, intent(in) :: tairk(mza) ! air temperature [K]
real, intent(in) :: rhov(mza)  ! water vapor density [kg_vap/m^3]

real, intent(out) :: rx (mza,ncat)  ! hydrom bulk spec dens [kg_hyd/kg_air]
real, intent(out) :: cx (mza,ncat)  ! hydrom bulk number [num_hyd/kg_air]
real, intent(out) :: emb(mza,ncat)  ! hydrom mean particle mass [kg/particle]

integer :: k
integer :: icat
integer :: ihcat
integer :: ns
integer :: nt

real :: rhovslair
real :: relhum
real :: tairc
real :: parmi

real, external :: rhovsl

! If level <= 1, there is no condensate of any type in this simulation.

if (level <= 1) then

! Set mcat to 0 and return

   mcat = 0
   return
endif

! If level = 2, cloud water is the only form of condensate that may exist. 

if (level == 2) then

! Set mcat to 1

   mcat = 1

! Set jnmb flag for cloud water to 1

   jnmb(1) = 1
   
! In OLAM, with level = 2, cloud number concentration is specified in cparm
! or parm(1).  Diagnose cloud droplet mean mass.

   parmi = 1. / parm(1)
   do k = ka,mza-1
      rx(k,1) = sh_c(k,iw)
      emb(k,1) = rx(k,1) * parmi
      cx(k,1) = parm(1)

      jhcat(k,1) = 1
   enddo

! Return

   return
   
endif

! If level = 3, up to 8 forms of condensate that may exist.

if (level == 3) then

! Set mcat to 8.

   mcat = 8

! Zero out microphysics scratch arrays for the present iw column

   rx(:,:) = 0.
   cx(:,:) = 0.

! Copy hydrometeor bulk mass and number concentration from main model arrays
! to microphysics column arrays rx and cx

! Cloud water

   if (jnmb(1) >= 1) then
      do k = ka,mza-1
         if (sh_c(k,iw) >= rxmin(1)) then

! If cloud bulk density is sufficiently abundant, copy to rx.

            rx(k,1) = sh_c(k,iw)

! If cloud water number concentration is prognosed, copy to cx.

            if (jnmb(1) >= 5) cx(k,1) = con_c(k,iw)

         endif
      enddo
   endif

! Rain

   if (jnmb(2) >= 1) then
      do k = ka,mza-1
         if (sh_r(k,iw) >= rxmin(2)) then

! If rain bulk density is sufficiently abundant, copy to rx,

            rx(k,2) = sh_r(k,iw)

! If rain water number concentration is prognosed, copy to cx.

            if (jnmb(2) >= 5) cx(k,2) = con_r(k,iw)

         endif
      enddo
   endif

! Pristine ice

   if (jnmb(3) >= 1) then
      do k = ka,mza-1
         if (sh_p(k,iw) >= rxmin(3)) then

! If pristine ice bulk density is sufficiently abundant, copy to rx.

            rx(k,3) = sh_p(k,iw)

! If pristine ice number concentration is prognosed, copy to cx.

            if (jnmb(3) >= 5) cx(k,3) = con_p(k,iw)

         endif
      enddo
   endif

! Snow

   if (jnmb(4) >= 1) then
      do k = ka,mza-1
         if (sh_s(k,iw) >= rxmin(4)) then

! If snow bulk density is sufficiently abundant, copy to rx.

            rx(k,4) = sh_s(k,iw)

! If snow number concentration is prognosed, copy to cx.

            if (jnmb(4) >= 5) cx(k,4) = con_s(k,iw)

         endif
      enddo
   endif

! Aggregates

   if (jnmb(5) >= 1) then
      do k = ka,mza-1
         if (sh_a(k,iw) >= rxmin(5)) then

! If aggregates bulk density is sufficiently abundant, copy to rx.

            rx(k,5) = sh_a(k,iw)

! If aggregates number concentration is prognosed, copy to cx.

            if (jnmb(5) >= 5) cx(k,5) = con_a(k,iw)

         endif
      enddo
   endif

! Graupel

   if (jnmb(6) >= 1) then
      do k = ka,mza-1
         if (sh_g(k,iw) >= rxmin(6)) then

! If graupel bulk density is sufficiently abundant, copy to rx,

            rx(k,6) = sh_g(k,iw)

! If graupel number concentration is prognosed, copy to cx.

            if (jnmb(6) >= 5) cx(k,6) = con_g(k,iw)

         endif
      enddo
   endif

! Hail

   if (jnmb(7) >= 1) then
      do k = ka,mza-1
         if (sh_h(k,iw) >= rxmin(7)) then

! If hail bulk density is sufficiently abundant, copy to rx,

            rx(k,7) = sh_h(k,iw)

! If hail number concentration is prognosed, copy to cx.

            if (jnmb(7) >= 5) cx(k,7) = con_h(k,iw)

         endif
      enddo
   endif

! Drizzle

   if (jnmb(8) >= 1) then
      do k = ka,mza-1
         if (sh_d(k,iw) >= rxmin(8)) then

! If drizzle bulk density is sufficiently abundant, copy to rx,

            rx(k,8) = sh_d(k,iw)

! If drizzle number concentration is prognosed, copy to cx.

            if (jnmb(8) >= 5) cx(k,8) = con_d(k,iw)

         endif
      enddo
   endif

! Diagnose pristine ice and snow habits from atmospheric temperature and humidity.
! This section of code copied or adapted from subroutines THRMSTR in omic_vap.f90 
! and EACH_COLUMN in omic_misc.f90

   do k = ka,mza-1
      tairc = tairk(k) - 273.15
      rhovslair = rhovsl(tairc)
      relhum = min(1.,rhov(k) / rhovslair)

      ns = max(1,nint(100. * relhum))
      nt = max(1,min(31,-nint(tairc)))

      jhcat(k,1) = 1
      jhcat(k,2) = 2
      jhcat(k,3) = jhabtab(nt,ns,1)
      jhcat(k,4) = jhabtab(nt,ns,2)
      jhcat(k,5) = 5
      jhcat(k,6) = 6
      jhcat(k,7) = 7
      jhcat(k,8) = 8
   enddo

! Loop over all hydrometeor categories

   do icat = 1,ncat

! Evaluate hydrometeor mean mass emb and concentration cx   
   
! This section of code was developed from subroutine enemb in omic_misc.f90
! by removing parts that are not needed for radiation calculations.  
! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR, 
! NOT PER M^3 AS IN THE ORIGINAL SUBROUTINE MIC_COPY.

      if (jnmb(icat) == 2) then

         do k = ka,mza-1
!!          ihcat = jhcat(k,icat)
!!          emb(k,icat) = emb2(ihcat)
            emb(k,icat) = emb2(icat)
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      elseif (jnmb(icat) == 4) then

         parmi = 1. / parm(icat)
         do k = ka,mza-1
            emb(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat) * parmi))
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      elseif (jnmb(icat) >= 5) then

         do k = ka,mza-1
            emb(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat)  &
                        / max(1.e-12,cx(k,icat))))
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      endif

   enddo

endif

return
end subroutine cloudprep_rad

!===============================================================================

subroutine cloud_opt(iw,ka,nrad,koff,mcat,jhcat,cx,emb,tp,omgp,gp,time)

use mem_harr, only: mb, nb, ocoef, bcoef, gcoef, nsolb

use micro_coms, only: ncat, jnmb, pwmasi, dnfac

use consts_coms, only: p00i, rocp

use mem_micro, only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h,  &
                     con_c, con_r, con_p, con_s, con_a, con_g, con_h

use mem_basic, only: rho

use mem_grid,  only: mza, dzt

use mem_radiate, only: rad_region

! computing properties of spherical liquid water and irregular ice
! using fits to adt theory
!
! ib .......... band number
! mb .......... maximum number of bands
! nb .......... total number of bands
! mza.......... number of vertical levels
! dzl ......... delta z in each level (m)
! dn .......... characteristic diameter (m)
! emb ......... mean hydrometeor mass (kg)
! cx .......... hydrometeor concentration (#/kg)
! ocoef ....... scattering albedo fit coefficients
! bcoef ....... extinction fit coefficients
! gcoef ....... asymmetry fit coefficients
! ncog ........ number of fit coefficients (omega and asym)
! ncb ......... number of fit coefficients (extinction)
! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                 numbers as a function of 15 microphysics category numbers
! rho ......... model air density (kg/m^3)
! dnfac ....... factor for computing dn from emb
! pwmasi ...... inverse of power used in mass power law

implicit none

integer, intent(in) :: iw
integer, intent(in) :: ka
integer, intent(in) :: nrad
integer, intent(in) :: koff
integer, intent(in) :: mcat

real, intent(in) :: time

integer, intent(in) :: jhcat(mza,ncat)

real, intent(in) :: cx(mza,ncat)
real, intent(in) :: emb(mza,ncat)

real, intent(out) :: tp(nrad,mb)   ! optical depth
real, intent(out) :: omgp(nrad,mb) ! scattering albedo
real, intent(out) :: gp(nrad,mb)   ! asymmetry parameter (assumes all 
                                     !   particles are spherical)

integer ib,iz,krc
integer icat,k,ihcat,krad

real :: dn  ! hydrometeor characteristic diameter [microns]
real :: ext
real :: om
real :: gg

! Use the following dn limiters (rather than those in microphysics) for 
! consistency with fit coefficients

real, parameter :: dnmin(8) = (/    1.,    10.,   1.,   125.,    10.,    10.,    10.,    1. /)
real, parameter :: dnmax(8) = (/ 1000., 10000., 125., 10000., 10000., 10000., 10000., 1000. /)

! Array kradcat maps RAMS/OLAM microphysics hydrometeor categories to those
! represented in Harrington radiation code according to the following numbering:

!     Harrington radiation code             Microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   rain                        2.  rain
!  3:   pristine ice columns        3.  pristine ice columns
!  4:   pristine ice rosettes       4.  snow columns
!  5:   pristine ice plates         5.  aggregates
!  6:   snow columns                6.  graupel
!  7:   snow rosettes               7.  hail
!  8:   snow plates                 8.  drizzle
!  9:   aggregates columns          9.  pristine ice hexagonal plates
!  10:  aggregates rosettes        10.  pristine ice dendrites
!  11:  aggregates plates          11.  pristine ice needles
!  12:  graupel                    12.  pristine ice rosettes
!  13:  hail                       13.  snow hexagonal plates
!                                  14.  snow dendrites
!                                  15.  snow needles
!                                  16.  snow rosettes

integer, save :: kradcat(16) = (/1,2,3,6,10,12,13,1,5,5,3,4,8,8,6,7/)

! Initialize arrays to zero prior to summation over any hydrometeor species

tp  (:,:) = 0.0
omgp(:,:) = 0.0
gp  (:,:) = 0.0

! Loop over active (mcat) hydrometeor categories

do icat = 1,mcat
   if (jnmb(icat) > 0) then

      do k = ka,mza-1
         krad = k - koff

         if (cx(k,icat) > 1.e-9) then

            ihcat = jhcat(k,icat)
            krc = kradcat(ihcat)
            dn = 1.e6 * dnfac(ihcat) * emb(k,icat) ** pwmasi(ihcat)
            dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns

            do ib = 1,nb

               ext = cx(k,icat) * rho(k,iw) * dzt(k)  &
                  * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)

               om = ocoef(1,ib,krc)  &
                  + ocoef(2,ib,krc) * exp(ocoef(3,ib,krc) * dn)  &
                  + ocoef(4,ib,krc) * exp(ocoef(5,ib,krc) * dn)

               gg = gcoef(1,ib,icat)  &
                  + gcoef(2,ib,icat) * exp(gcoef(3,ib,icat) * dn)  &
                  + gcoef(4,ib,icat) * exp(gcoef(5,ib,icat) * dn)

!--------------------------------------------------------------------
! THIS SECTION ADDED BASED ON OPTIMIZATION PROCEDURE

!               if (ib <= nsolb) then
!                  if (icat == 1 .or. icat == 3) then
!                     gg = gg * 1.143
!                     if (rad_region(iw) == 1) then
!                        gg = gg * 1.13
!                     elseif (rad_region(iw) == 2) then
!                        gg = gg * 1.10
!                     elseif (rad_region(iw) == 3) then
!                        gg = gg * 0.50
!                     endif
!                  endif
!               endif
!--------------------------------------------------------------------

               tp(krad,ib) = tp(krad,ib) + ext

               omgp(krad,ib) = omgp(krad,ib) + om * ext
               gp(krad,ib) = gp(krad,ib) + gg * om * ext

            enddo

         endif
      enddo

   endif
enddo

! Combine the optical properties....

do ib = 1,nb

   do k = ka,mza-1
      krad = k - koff

      if (tp(krad,ib) > 1.e-15 .and. omgp(krad,ib) > 1.e-15) then
         gp  (krad,ib) = min(0.9999, gp  (krad,ib) / omgp(krad,ib))
         omgp(krad,ib) = min(0.9999, omgp(krad,ib) / tp  (krad,ib))
      else
         omgp(krad,ib) = 0.0
         gp  (krad,ib) = 0.0
      endif

! Check for validity of opt values before calling radiation
!      if (tp(krad,ib) < 0) then
!         write(io6,*) 'tp(krad,ib) less than zero for krad,ib = ',krad,ib
!         write(io6,*) 'tp(krad,ib) = ',tp(krad,ib)
!         stop 'opt1'
!      endif
!      if (omgp(krad,ib) < 0. .or. omgp(krad,ib) > 1.) then
!         write(io6,*) 'omgp(krad,ib) out of range [0,1] for krad,ib = ',krad,ib
!         write(io6,*) 'omgp(krad,ib) = ',omgp(krad,ib)
!         stop 'opt2'
!      endif
!      if (gp(krad,ib) < 0. .or. gp(krad,ib) > 1.) then
!         write(io6,*) 'gp(krad,ib) out of range [0,1] for krad,ib = ',krad,ib
!         write(io6,*) 'gp(krad,ib) = ',gp(krad,ib)
!         stop 'opt3'
!      endif

   enddo
enddo

return
end subroutine cloud_opt

!===============================================================================

subroutine path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl)

! Get the path lengths for the various gases...

use consts_coms, only: grav, eps_virt

implicit none

integer :: nrad

real, intent(out) :: u  (nrad,3)
real, intent(out) :: vp (nrad)

real, intent(in)  :: rl (nrad)
real, intent(in)  :: dzl(nrad)
real, intent(in)  :: dl (nrad)
real, intent(in)  :: o3l(nrad)
real, intent(in)  :: pl (nrad)

real, parameter :: eps_rad = 1.e-15

real :: rvk0,rvk1,dzl9,rmix
integer :: k
real, parameter :: co2_mixing_ratio = 360.0e-6 * 44.011 / 28.966 ! [kg/kg]

u(1,1) = .5 * (rl(2) + rl(1)) * grav * dzl(1)
u(1,2) = .5 * (dl(2) + dl(1)) * co2_mixing_ratio * grav * dzl(1)
u(1,3) = o3l(1) * grav * dzl(1)

do k = 2,nrad
   dzl9   = grav * dzl(k)
   rmix = rl(k) / dl(k)
   vp(k)  = pl(k) * rmix / (eps_virt + rmix)
   u(k,1) = 0.5 * dzl9 * (rl(k) + rl(k-1))
   u(k,2) = 0.5 * dzl9 * (dl(k) + dl(k-1)) * co2_mixing_ratio
   u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
enddo

vp(1) = vp(2)

return
end subroutine path_lengths

!===============================================================================

subroutine define_rad_regions()

  !  This subroutine classifies each grid cell as:
  !  (1) Polar  (mreg = 8)
  !  (2) Temperate land (mreg = 2)
  !  (3) Temperate ocean (mreg = 7)
  !  (4) Tropical ocean (mreg = 4)
  !  (5) Tropical land (mreg = 5)

  !  However, the current radiation parameters have only a 
  !  latitude-dependence, not a land/sea dependence.  Thus, their are 
  !  only 3 radiation regions:
  !  (1)  rad_region = 1 corresponds to mreg = 4 and 5 (tropical)
  !  (2)  rad_region = 2 corresponds to mreg = 2 and 7 (temperate)
  !  (3)  rad_region = 3 corresponds to mreg = 8 (polar)

  use mem_grid, only: mwa, glatw, glonw
  
  use mem_radiate, only: rad_region

  implicit none

  integer :: iw, mreg
  real :: lat, lon
 
  do iw = 2, mwa
     lat = glatw(iw)
     lon = glonw(iw)
     if(lat > 70.0)then
        mreg = 8
     elseif(lat > 60.0)then
        if(lon < -160.0)then
           mreg = 8
        elseif(lon < -60.0)then
           mreg = 2
        elseif(lon < 10.0)then
           mreg = 8
        else
           mreg = 2
        endif
     elseif(lat > 50.0)then
        if(lon < -130.0)then
           mreg = 7
        elseif(lon < -60.0)then
           mreg = 2
        elseif(lon < -10.0)then
           mreg = 7
        elseif(lon < 0.0)then
           mreg = 2
        elseif(lon < 10.0)then
           mreg = 7
        elseif(lon < 140.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > 40.0)then
        if(lon < -120.0)then
           mreg = 7
        elseif(lon < -60.0)then
           mreg = 2
        elseif(lon < -10.0)then
           mreg = 7
        elseif(lon < 130.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > 30.0)then
        if(lon < -120.0)then
           mreg = 7
        elseif(lon < -80.0)then
           mreg = 2
        elseif(lon < -10.0)then
           mreg = 7
        elseif(lon < 120.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > 20.0)then
        if(lon < -110.0)then
           mreg = 7
        elseif(lon < -100.0)then
           mreg = 2
        elseif(lon < -10.0)then
           mreg = 7
        elseif(lon < 120.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > 10.0)then
        if(lon < -90.0)then
           mreg = 4
        elseif(lon < -80.0)then
           mreg = 5
        elseif(lon < -20.0)then
           mreg = 4
        elseif(lon < 50.0)then
           mreg = 5
        elseif(lon < 70.0)then
           mreg = 4
        elseif(lon < 80.0)then
           mreg = 5
        elseif(lon < 90.0)then
           mreg = 4
        elseif(lon < 130.0)then
           mreg = 5
        else
           mreg = 4
        endif
     elseif(lat > 0.0)then
        if(lon < -80.0)then
           mreg = 4
        elseif(lon < -50.0)then
           mreg = 5
        elseif(lon < -10.0)then
           mreg = 4
        elseif(lon < 50.0)then
           mreg = 5
        elseif(lon < 90.0)then
           mreg = 4
        elseif(lon < 130.0)then
           mreg = 5
        else
           mreg = 4
        endif
     elseif(lat > -10.0)then
        if(lon < -80.0)then
           mreg = 4
        elseif(lon < -30.0)then
           mreg = 5
        elseif(lon < 10.0)then
           mreg = 4
        elseif(lon < 40.0)then
           mreg = 5
        elseif(lon < 90.0)then
           mreg = 4
        elseif(lon < 140.0)then
           mreg = 5
        else
           mreg = 4
        endif
     elseif(lat > -20.0)then
        if(lon < -80.0)then
           mreg = 4
        elseif(lon < -40.0)then
           mreg = 5
        elseif(lon < 10.0)then
           mreg = 4
        elseif(lon < 40.0)then
           mreg = 5
        elseif(lon < 130.0)then
           mreg = 4
        elseif(lon < 150.0)then
           mreg = 5
        else
           mreg = 4
        endif
     elseif(lat > -30.0)then
        if(lon < -70.0)then
           mreg = 7
        elseif(lon < -50.0)then
           mreg = 2
        elseif(lon < 10.0)then
           mreg = 7
        elseif(lon < 30.0)then
           mreg = 2
        elseif(lon < 110.0)then
           mreg = 7
        elseif(lon < 150.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > -40.0)then
        if(lon < -70.0)then
           mreg = 7
        elseif(lon < -50.0)then
           mreg = 2
        elseif(lon < 20.0)then
           mreg = 7
        elseif(lon < 30.0)then
           mreg = 2
        elseif(lon < 120.0)then
           mreg = 7
        elseif(lon < 150.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > -50.0)then
        if(lon < -70.0)then
           mreg = 7
        elseif(lon < -60.0)then
           mreg = 2
        else
           mreg = 7
        endif
     elseif(lat > -60.0)then
        mreg = 7
     else
        mreg = 8
     endif     
     if(mreg == 4 .or. mreg == 5)rad_region(iw) = 1
     if(mreg == 2 .or. mreg == 7)rad_region(iw) = 2
     if(mreg == 8)rad_region(iw) = 3
  enddo

  return
end subroutine define_rad_regions

