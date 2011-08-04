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
Module grell_deep_coms

  integer, parameter :: maxens = 3       !  ensemble one on cap_max
  integer, parameter :: maxens2 = 3     ! ensemble two on precip efficiency
  integer, parameter :: maxens3 = 16    ! ensemble three done in cup_forcing
  integer :: ensdim = maxens*maxens2*maxens3
  real, parameter :: entr_rate = 0.2 / 12000.0
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
  real, parameter :: cap_maxs = 100.0
!--- depth(m) over which downdraft detrains all its mass
  real, parameter :: z_detr=1250.
  real, parameter :: cap_max_increment = 50.
!--- height(m) above which no downdrafts are allowed to originate
  real, parameter :: zcutdown=3000.
!--- minimum depth (m), clouds must have
  real, parameter :: depth_min = 500.0
!--- max height(m) above ground where updraft air can originate
  real, parameter :: zkbmax=4000.
!---  max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!     base mass flux
  real, parameter :: edtmax = .95
  real, parameter :: edtmin = .2


End Module grell_deep_coms
