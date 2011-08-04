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
Module grell_shallow_coms

  real, parameter :: entr_rate = 0.2 / 666.0 ! gross entrainment rate
  real, parameter :: cap_maxs = 120.0
  integer, parameter :: maxens =3  !3  ensemble one on mbdt
  integer, parameter :: maxens2=1  !1 ensemble two on precip efficiency
  integer, parameter :: maxens3=10 !10  ensemble three done in cup_forcing
  integer, save :: ensdim = maxens * maxens2 * maxens3
  real, parameter :: cap_inc = 0.0
  real, parameter :: zkbmax = 4000.0  ! max height(m) above ground where shallow clouds can originate
!srf- ICOIC is used for choice a specific closure for shallow cumulus
! icoic = 0 -> ensemble (all closures)
! icoic = 1 -> Grell
! icoic = 4 -> Arakawa-Schubert
! icoic = 8 -> like Fritsch Chappel or Kain Fritsch
  integer, parameter :: icoic = 8

end Module grell_shallow_coms
