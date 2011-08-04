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
subroutine ccmp_raddriv(iw,ka,koff)

use mem_basic,   only: press, theta, sh_v, sh_w, rho

use mem_radiate, only: solfac, sunx, suny, sunz, cosz,                       &
                       rlongup, rlong, albedt, albedt_beam, albedt_diffuse,  &
                       rshort, rshort_diffuse, fthrd, rshort, rlong_albedo,  &
                       rlongup_top

use misc_coms,   only: io6, ilwrtyp, iswrtyp, time8, itime1
use consts_coms, only: p00, rocp, stefan, solar, pio180, erad
use mem_grid,    only: dzim, dzit, glonw, mza
use micro_coms,  only: level, ipris
use mem_micro,   only: sh_c, sh_p

implicit none

integer, intent(in) :: iw
integer, intent(in) :: ka
integer, intent(in) :: koff

integer :: j
integer :: k
integer :: krad
integer :: nlev

real :: r1  ! Radial expansion factor due to topography at im1
real :: r2  ! Radial expansion factor due to topography at im2
real :: r3  ! Radial expansion factor due to topography at im3

real :: tnx  ! x-component of unit vector normal to mean terrain [m]
real :: tny  ! y-component of unit vector normal to mean terrain [m]
real :: tnz  ! z-component of unit vector normal to mean terrain [m]

real :: tcosz  ! cosine of angle between topography-normal and sun unit vectors

! Local arrays

real :: rvr   (mza+1) 
real :: rcr   (mza+1) 
real :: dn0r  (mza+1) 
real :: pird  (mza+1) 
real :: prd   (mza+1) 
real :: fthrl (mza+1) 
real :: dzmr  (mza+1) 
real :: dztr  (mza+1) 
real :: fthrs (mza+1) 
real :: temprd(mza+1) 
!  Since these schemes do not calculate beam and diffuse radiation separately,
!  we here assume a default diffuse fraction.
real, parameter :: diffuse_fraction = 0.5

! NLEV (number of radiative levels) includes predicted T levels (k=ka to k=mza-1)
! plus an extra bottom level as an auxiliary memory space.  
   
nlev = mza - ka + 1

! Loop over predicted model T levels

do k = ka,mza-1
   
! KRAD (radiative level) equals model level k when ka = 2 because of extra
! bottom level in radiation column vectors.   
   
   krad = k - koff

   pird(krad)   = (press(k,iw) / p00) ** rocp
   temprd(krad) = theta(k,iw) * pird(krad)
   rvr(krad)    = max(0.0, sh_v(k,iw))

! Convert the next 4 variables to cgs for now.

   prd(krad)  = press(k,iw) * 10.
   dn0r(krad) = rho(k,iw) * 1.e-3
   dzmr(krad) = dzim(k) * 1.e-2
   dztr(krad) = dzit(k) * 1.e-2

enddo

if (level >= 3 .and. ipris >= 1) then

! Fill rcr with cloud water mixing ratio plus a fraction of pristine ice
! mixing ratio.  CHEN-COTTON ASSUMES THIS CONDENSATE (I.E., RCR) HAS THE 
! RADIATIVE PROPERTIES OF CLOUD WATER.

   do k = ka,mza-1
      krad = k - koff
      rcr(krad) = max(0.0, sh_c(k,iw)+0.1*sh_p(k,iw) )
   enddo

elseif (level >= 2) then

! USE THE CLOUD WATER FROM THE MICROPHYSICS SCHEME

   do k = ka,mza-1
      krad = k - koff
      rcr(krad) = max(0.0, sh_c(k,iw))
   enddo

else

! NO CLOUDS ALLOWED FOR LEVEL=1

   do k = ka,mza-1
      krad = k - koff
      rcr(krad) = 0.0
   enddo

endif

pird(1) = pird(2)
rvr(1)  = rvr(2)
prd(1)  = prd(2)
dn0r(1) = dn0r(2)
dzmr(1) = dzmr(2) 
dztr(1) = dztr(2)
rcr(1)  = rcr(2)

temprd(1) = sqrt(sqrt(rlongup(iw) / stefan))
temprd(nlev+1) = temprd(nlev)

! Call the longwave parameterizations.

if (ilwrtyp == 2) then

   call lwradp(nlev,temprd,rvr,dn0r,dztr,pird,fthrl,rlong(iw))

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw) + fthrl(krad)
   enddo

! Convert the downward flux at the ground to SI.

   rlong(iw) = rlong(iw) * 1.e-3

elseif (ilwrtyp == 1) then

   rlongup(iw) = rlongup(iw) * 1.0e3   ! converts to cgs

   call lwradc(nlev+1,rvr,rcr,dn0r,temprd,prd,dztr,fthrl,rlong(iw)  &
      ,rlongup(iw),rlong_albedo(iw),rlongup_top(iw))
   rlongup(iw) = rlongup(iw) * 1.0e-3  ! converts to SI (mks)

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw) + fthrl(krad)
   enddo

! Convert the downward flux at the ground to SI.

   rlong(iw) = rlong(iw) * 1.e-3
   rlongup_top(iw) = rlongup_top(iw) * 1.0e-3

endif

! Call the shortwave parameterizations.

! Compute the weighted shortwave albedo.
albedt(iw) = albedt_beam(iw) + diffuse_fraction * (albedt_diffuse(iw) -   &
     albedt_beam(iw))

! The shortwave parameterizations are only valid if the cosine
!    of the zenith angle is greater than .03 .

if (iswrtyp == 2 .and. cosz(iw) > .03) then

   call shradp(nlev,rvr,dn0r,dzmr,pird,cosz(iw)  &
      ,albedt(iw),solar*1e3*solfac,fthrs,rshort(iw))

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw) + fthrs(krad)
   enddo

! Convert the downward flux at the ground to SI.

   rshort(iw) = rshort(iw) * 1.e-3 / (1. - albedt(iw))
   rshort_diffuse(iw) = rshort(iw) * diffuse_fraction
   
elseif (iswrtyp == 1 .and. cosz(iw) > .03) then

   call shradc(nlev+1,rvr,rcr,dn0r,dztr,prd  &
      ,albedt(iw),solar*1.e3*solfac,cosz(iw),fthrs,rshort(iw))

   do k = ka,mza-1
      krad = k - koff
      fthrd(k,iw) = fthrd(k,iw) + fthrs(krad)
   enddo

! Convert the downward flux at the ground to SI.

   rshort(iw) = rshort(iw) * 1.e-3 / (1. - albedt(iw))
   rshort_diffuse(iw) = rshort(iw) * diffuse_fraction
   
endif


return
end subroutine ccmp_raddriv

