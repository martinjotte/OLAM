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

subroutine cloudprep_rad(iw,ka,mcat,jhcat,rhov,rx,cx,emb)

! This subroutine was developed from parts of subroutine MIC_COPY in
! omic_driv.f90 and subroutines EACH_COLUMN and ENEMB in omic_misc.f90.

! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR, 
! NOT PER M^3 AS IN THE ORIGINAL MICROPHYSICS SUBROUTINES.

use micro_coms, only: ncat, jnmb, rxmin, jhabtab, miclevel,  &
                      emb0, emb1, emb2, zfactor_ccn

use mem_micro,  only: sh_c, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, sh_d,  &
                      con_c, con_r, con_p, con_s, con_a, con_g, con_h, con_d, &
                      cldnum

use mem_basic,  only: tair

use mem_grid,   only: mza

use therm_lib,  only: rhovsl

implicit none

integer, intent(in) :: iw
integer, intent(in) :: ka

integer, intent(out) :: mcat  ! # of active hydrom categories (0,1, or 7)
integer, intent(out) :: jhcat(mza,ncat)  ! hydrom category table with ice habits

real, intent(in) :: rhov(mza)  ! water vapor density [kg_vap/m^3]

real, intent(out) :: rx (mza,ncat)  ! hydrom bulk spec dens [kg_hyd/kg_air]
real, intent(out) :: cx (mza,ncat)  ! hydrom bulk number [num_hyd/kg_air]
real, intent(out) :: emb(mza,ncat)  ! hydrom mean particle mass [kg/particle]

real :: con_ccnx(mza)  ! Local array for CCN/cloud concentration [#/kg_air]

integer :: k
integer :: icat
integer :: ihcat
integer :: ns
integer :: nt

real :: rhovslair
real :: relhum
real :: tairc

! If miclevel <= 1, there is no condensate of any type in this simulation.

if (miclevel <= 1) then

! Set mcat to 0 and return

   mcat = 0
   return
endif

! If miclevel = 2, cloud water is the only form of condensate that may exist. 

if (miclevel == 2) then

! Set mcat to 1

   mcat = 1

! Set jnmb flag for cloud water to 1

   jnmb(1) = 1
   
! In OLAM, with miclevel = 2, cloud number concentration is specified in cldnum.
! Diagnose cloud droplet mean mass.

   do k = ka,mza
      rx(k,1) = sh_c(k,iw)
      cx(k,1) = cldnum(iw) * zfactor_ccn(k)
      emb(k,1) = rx(k,1) / cx(k,1)

      jhcat(k,1) = 1
   enddo

! Return

   return
   
endif

! If miclevel = 3, up to 8 forms of condensate that may exist.

if (miclevel == 3) then

! Set mcat to 8.

   mcat = 8

! Zero out microphysics scratch arrays for the present iw column

   rx(:,:) = 0.
   cx(:,:) = 0.

! Copy hydrometeor bulk mass and number concentration from main model arrays
! to microphysics column arrays rx and cx

! Cloud water

   if (jnmb(1) >= 1) then
      do k = ka,mza
         if (sh_c(k,iw) >= rxmin(1)) then

! If cloud bulk density is sufficiently abundant, copy to rx.

            rx(k,1) = sh_c(k,iw)

! If cloud water number concentration is prognosed, copy to cx.

            if (jnmb(1) == 5) cx(k,1) = con_c(k,iw)

         endif
      enddo
   endif

! Rain

   if (jnmb(2) >= 1) then
      do k = ka,mza
         if (sh_r(k,iw) >= rxmin(2)) then

! If rain bulk density is sufficiently abundant, copy to rx,

            rx(k,2) = sh_r(k,iw)

! If rain water number concentration is prognosed, copy to cx.

            if (jnmb(2) == 5) cx(k,2) = con_r(k,iw)

         endif
      enddo
   endif

! Pristine ice

   if (jnmb(3) >= 1) then
      do k = ka,mza
         if (sh_p(k,iw) >= rxmin(3)) then

! If pristine ice bulk density is sufficiently abundant, copy to rx.

            rx(k,3) = sh_p(k,iw)

! If pristine ice number concentration is prognosed, copy to cx.

            if (jnmb(3) == 5) cx(k,3) = con_p(k,iw)

         endif
      enddo
   endif

! Snow

   if (jnmb(4) >= 1) then
      do k = ka,mza
         if (sh_s(k,iw) >= rxmin(4)) then

! If snow bulk density is sufficiently abundant, copy to rx.

            rx(k,4) = sh_s(k,iw)

! If snow number concentration is prognosed, copy to cx.

            if (jnmb(4) == 5) cx(k,4) = con_s(k,iw)

         endif
      enddo
   endif

! Aggregates

   if (jnmb(5) >= 1) then
      do k = ka,mza
         if (sh_a(k,iw) >= rxmin(5)) then

! If aggregates bulk density is sufficiently abundant, copy to rx.

            rx(k,5) = sh_a(k,iw)

! If aggregates number concentration is prognosed, copy to cx.

            if (jnmb(5) == 5) cx(k,5) = con_a(k,iw)

         endif
      enddo
   endif

! Graupel

   if (jnmb(6) >= 1) then
      do k = ka,mza
         if (sh_g(k,iw) >= rxmin(6)) then

! If graupel bulk density is sufficiently abundant, copy to rx,

            rx(k,6) = sh_g(k,iw)

! If graupel number concentration is prognosed, copy to cx.

            if (jnmb(6) == 5) cx(k,6) = con_g(k,iw)

         endif
      enddo
   endif

! Hail

   if (jnmb(7) >= 1) then
      do k = ka,mza
         if (sh_h(k,iw) >= rxmin(7)) then

! If hail bulk density is sufficiently abundant, copy to rx,

            rx(k,7) = sh_h(k,iw)

! If hail number concentration is prognosed, copy to cx.

            if (jnmb(7) == 5) cx(k,7) = con_h(k,iw)

         endif
      enddo
   endif

! Drizzle

   if (jnmb(8) >= 1) then
      do k = ka,mza
         if (sh_d(k,iw) >= rxmin(8)) then

! If drizzle bulk density is sufficiently abundant, copy to rx,

            rx(k,8) = sh_d(k,iw)

! If drizzle number concentration is prognosed, copy to cx.

            if (jnmb(8) == 5) cx(k,8) = con_d(k,iw)

         endif
      enddo
   endif

! Diagnose pristine ice and snow habits from atmospheric temperature and humidity.
! This section of code copied or adapted from subroutines THRMSTR in omic_vap.f90 
! and EACH_COLUMN in omic_misc.f90

   do k = ka,mza
      tairc = tair(k,iw) - 273.15
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

         do k = ka,mza
            emb(k,icat) = emb2(icat)
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      elseif (jnmb(icat) == 4) then ! As of version 5.0.0, can only apply to cloud

         do k = ka,mza
            con_ccnx(k) = cldnum(iw) * zfactor_ccn(k)
            emb(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat) / con_ccnx(k)))
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      elseif (jnmb(icat) == 5) then

         do k = ka,mza
            emb(k,icat) = max(emb0(icat), &
                          min(emb1(icat),rx(k,icat) / max(1.e-12,cx(k,icat))))
            cx(k,icat) = rx(k,icat) / emb(k,icat)
         enddo

      endif

   enddo

endif

end subroutine cloudprep_rad
