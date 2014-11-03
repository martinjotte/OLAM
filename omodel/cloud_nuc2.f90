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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Bob notes 9/24/2014: Subroutine cldnuc2 begun as alternative to subroutine
! cldnud.  Cldnuc2 intended to interface with subroutine ccnbin, the new
! "on-line" code for nucleating CCN during a model run.  Suspended development
! until further considering Lagrangian parcel formulation, which is the 
! correct way to consider nucleation and environmental forcing rates.  Code
! below that relates a1, a1inv, vertical velocity, and supsat tendency are all
! central to this consideration.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


subroutine cldnuc2(iw0,lpw0,dtli0, &
                  rx,cx,qr,qx,con_ccnx,con_gccnx,rhov,rhoi,rhoa, &
                  tair,tairc,wc0,rhovslair,rnuc_vc,rnuc_vd,cnuc_vc,cnuc_vd)

use micro_coms,  only: jnmb, parm, emb0, emb1, cparm, mza0, ncat, &
                       iccnlev, dparm, cnparm, gnparm, rxmin
use misc_coms,   only: io6
use mem_grid,    only: zt
use consts_coms, only: rvap, grav, alvl, cp, cliq, alli, eps_vapi

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

real, intent(in) :: dtli0

real, intent(inout) :: rx(mza0,ncat)
real, intent(inout) :: cx(mza0,ncat)
real, intent(inout) :: qr(mza0,ncat)
real, intent(inout) :: qx(mza0,ncat)

real, intent(inout) :: con_ccnx (mza0)
real, intent(inout) :: con_gccnx(mza0)

real, intent(inout) :: rhov  (mza0)
real, intent(in) :: rhoi     (mza0)
real, intent(in) :: tair     (mza0)
real, intent(in) :: tairc    (mza0)
real, intent(in) :: wc0      (mza0)
real, intent(in) :: rhovslair(mza0)

real, intent(inout) :: cnuc_vc(mza0)
real, intent(inout) :: cnuc_vd(mza0)
real, intent(inout) :: rnuc_vc(mza0)
real, intent(inout) :: rnuc_vd(mza0)

real(kind=8), intent(in) :: rhoa(mza0)

integer :: k,jtemp,jw,jcon,iw,iconc

real :: rnuc,cnuc,excessrhov,rhocnew,tab,concen_tab,cxadd, &
   tairc_nuc,w_nuc,rjw,wtw2,wtw1,con_ccnk,con_gccnk,rjcon,wtcon2,wtcon1

real :: con_ccn_tab, con_gccn_tab, con_ccn_nuc, con_gccn_nuc, frac

integer :: ic,rgb,jrg,ct,ctc,maxct,sps

real :: grat,crat,wtrg1,wtrg2,gvap,cvap,rg,rg1,rg2
real :: vaprccn,vaprgccn,mult1,mult2,supsat,supsat_rate,a1inv,w_pseudo

! If NOT prognosing number concentration of cloud

if (jnmb(1) == 4) then

! Loop over all levels in column

   do k = lpw0,mza0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Initial and final environmental saturation values (actual final value will
! be less than satenv1 due to condensation onto CCN)

  satenv0 = rhov0 / (1.0e-3 * rhovsl(tempk0-273.15))
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Diagnose excess of vapor density over saturation; cycle if not > 0.

      excessrhov = rhov(k) - 1.00001 * rhovslair(k) ! x rhoa
      if (excessrhov <= 0.) cycle                   ! x rhoa

! parm(1) is specified cloud # per kg_air; use as upper bound on nucleation 

      cnuc = parm(1) * rhoa(k)   ! #_nucleated / m^3
      rnuc = cnuc * emb0(1)      ! kg_nucleated / m^3

      rnuc_vc(k) = min(rnuc,.5 * excessrhov)   ! x rhoa
      rx(k,1) = rx(k,1) + rnuc_vc(k)           ! x rhoa
      rhov(k) = rhov(k) - rnuc_vc(k)           ! x rhoa

      cnuc_vc(k) = min(cnuc,rx(k,1) / emb0(1)) ! x rhoa
      cx(k,1) = cnuc_vc(k)
      qr(k,1) = qr(k,1) + rnuc_vc(k) * (tairc(k) * cliq + alli) ! (x rhoa)

   enddo

! Saleeby(6/3/02) Prognosing number concentration of cloud and/or drizzle

elseif (jnmb(1) >= 5) then

! cloud number predicted from ccn field

! Loop over all levels in column

   do k = lpw0,mza0

! Diagnose excess of vapor density over saturation; cycle if not > 0.

      excessrhov = rhov(k) - 1.00001 * rhovslair(k) ! x rhoa
      if (excessrhov <= 0.) cycle                   ! x rhoa

      con_ccn_nuc = 0.

! CCN concentration
! (CON_CCNK HAS UNITS OF #/KG TO MATCH NUCLEATION PARAMETERIZATION)

      if (jnmb(1) == 5) then

! Constant in time and space CCN concentration

         con_ccnk = cparm       ! cparm units are #/kg 

      elseif (jnmb(1) == 6) then

! Horizontally-homogeneous, time-constant vertical profile of CCN
! (Example profile by Saleeby)

         con_ccnk = max(100.e6,cparm * (1. - zt(k) / 4000.))

      elseif (jnmb(1) == 7) then

! Prognostic CCN field

         con_ccnk = con_ccnx(k) * rhoi(k)   ! con_ccnx units are #/m^3

      endif  

! Do not allow con_ccnk to exceed upper limit

      if (con_ccnk > 9999.e6) con_ccnk = 9999.e6
      
! Nucleate only if CON_CCNK concentration exceeds lower threshold (10**7/kg_air)

      if (con_ccnk > 10.001e6) then


! VERTICAL VELOCITY indices and weights for CCN nucleation table
! (impose upper & lower bounds)
!---------------------------------------------------------------

! Inverse of A1 coefficient in Pruppacher and Klett Eqn. 13-31.

         a1inv = (rvap * tair(k)) / (grav * (alvl / (cp * tair(k)) - eps_vapi))

         a1 = grav * (alvl / (cp * tair(k)) - eps_vapi) / (rvap * tair(k))

         a1 = (grav / (rdry * tair(k))) * (alvl * eps_vap / (cp * tair(k)) - 1.) 

! Pseudo vertical velocity based on supersaturation and its assumed 
! production over 1 timestep

         supsat = rhov(k) / (1.00001 * rhovslair(k)) - 1.0

         supsat_rate = supsat * dtli0

         w_pseudo = a1inv * supsat_rate

! Make w_nuc the larger of wc0 and w_pseudo

         w_nuc = max(wc0(k),w_pseudo)

         if (w_nuc < .010001) then 
            w_nuc = .010001
         elseif (w_nuc > 99.99) then
            w_nuc = 99.99
         endif

! interpolate from CLDNUCTAB10 table

            tab = &    ! Fraction of CCN to nucleate

! Get number to nucleate

         con_ccn_tab = con_ccnk * tab * rhoa(k)   ! (#/m^3) (x rhoa)

! If you're NOT depleting aerosols, reduce nucleatable CCN concentration 
! by existing cloud concentration

! Bob's comment (11/16/2009): This should apply in 2 cases: Not depleting aerosols
! in nucleation or anywhere else, and not depleting aerosols in nucleation 
! because you're depleting them in collisions.  (I think it is better to never
! deplete aerosols in nucleation, but to deplete them in collisions instead.)

         con_ccn_nuc = con_ccn_tab
         if (iccnlev == 0) con_ccn_nuc = max(0.,con_ccn_tab - cx(k,1))

      endif ! con_ccnk > 0

! Vapor mass to nucleate to cloud droplets (subject to possible adjustment)

      vaprccn = con_ccn_nuc * emb0(1)

! If DRIZZLE number concentration is being predicted (which in the model
! also means that DRIZZLE nucleation from GCCN is active), excess vapor 
! must be shared between CLOUD and DRIZZLE nucleation.

      vaprgccn = 0.

      if (jnmb(8) >= 5) then

! GCCN concentration

         if (jnmb(8) == 5) then

! Constant in time and space GCCN concentration

            con_gccn_tab = dparm * rhoa(k)      ! dparm units are #/kg 

         elseif (jnmb(8) == 6) then

! Horizontally-homogeneous, time-constant vertical profile of GCCN
! (Example profile by Saleeby)

            con_gccn_tab = max(10.,dparm * (1. - zt(k) / 4000.)) * rhoa(k)

         elseif (jnmb(8) == 7) then

! Prognostic CCN field

            con_gccn_tab = con_gccnx(k)   ! con_gccnx units are #/m^3

         endif  

! If you're NOT depleting aerosols, reduce nucleatable GCCN concentration 
! by existing drizzle concentration

! Bob's comment (11/16/2009): This should apply in 2 cases: Not depleting aerosols
! in nucleation or anywhere else, and not depleting aerosols in nucleation 
! because you're depleting them in collisions.  (I think it is better to never
! deplete aerosols in nucleation, but to deplete them in collisions instead.)

         con_gccn_nuc = con_gccn_tab
         if (iccnlev == 0) con_gccn_nuc = max(0.,con_gccn_nuc - cx(k,8))

! Vapor mass to nucleate to cloud droplets (subject to possible adjustment)

         vaprgccn = con_gccn_nuc * emb0(8)

      endif

! If nucleation to cloud and/or drizzle is limited by available vapor,
! reduce nucleation accordingly

      if (vaprccn + vaprgccn > excessrhov) then

         frac = excessrhov / (vaprccn + vaprgccn)

         vaprccn  = vaprccn  * frac
         vaprgccn = vaprgccn * frac

         con_ccn_nuc  = vaprccn  / emb0(1)
         con_gccn_nuc = vaprgccn / emb0(8)
      endif

! Update cloud concentration and mass from nucleation

      if (vaprccn > 0.) then

         cnuc_vc(k) = con_ccn_nuc
         rnuc_vc(k) = vaprccn

         cx(k,1) = cx(k,1) + cnuc_vc(k) ! (#/m^3)  x rhoa
         rx(k,1) = rx(k,1) + rnuc_vc(k) ! (kg/m^3) x rhoa
         qr(k,1) = qr(k,1) + rnuc_vc(k) * (tairc(k) * cliq + alli) ! (x rhoa)

         if (rx(k,1) >= rxmin(1)) qx(k,1) = qr(k,1) / rx(k,1)

! Update CCN concentration if doing optional nucleation scavenging

         if (iccnlev == 1) con_ccnx(k) = con_ccnx(k) - con_ccn_nuc ! (#/m^3)

      endif ! (vaprccn > 0)

! Update drizzle concentration and mass from nucleation

      if (jnmb(8) >= 5 .and. vaprgccn > 0.) then

         cnuc_vd(k) = con_gccn_nuc
         rnuc_vd(k) = vaprgccn

         cx(k,8) = cx(k,8) + cnuc_vd(k) ! (#/m^3)  x rhoa
         rx(k,8) = rx(k,8) + rnuc_vd(k) ! (kg/m^3) x rhoa
         qr(k,8) = qr(k,8) + rnuc_vd(k) * (tairc(k) * cliq + alli) ! (x rhoa)

         if (rx(k,8) >= rxmin(8)) qx(k,8) = qr(k,8) / rx(k,8)

! Update GCCN concentration if doing optional nucleation scavenging

         if (iccnlev == 1) con_gccnx(k) = con_gccnx(k) - con_gccn_nuc ! (#/m^3)

      endif ! (jnmb(8) >= 5 .and. vaprgccn > 0.)

! If nucleatable CCN concentration is positive, carry out nucleation

   enddo ! Loop over k

endif ! (jnmb(1) >= 5)

return
end subroutine cldnuc2
