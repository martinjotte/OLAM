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
module c34constants
  implicit none

  type farqdata
     real :: D0
     real :: alpha
     real :: gamma
     real :: m
     real :: b
  end type farqdata

  type metdat
     real :: ea,ca,rn,tl,par,gbc,gbw,ta,el,compp,eta,gbci
  end type metdat

  type glim
     real :: sigma,rho,vm,tau,nu,k1,k2
  end type glim

  type solution
     real, dimension(2,2) :: es,ci,cs,a,gsw
     real :: eps
     integer :: ninterval
  end type solution

  Type stoma_data

     integer :: recalc=1
     real :: T_L
     real :: e_A
     real :: PAR
     real :: rb_factor
     real :: prss
     real :: phenology_factor
     real :: gsw_open
     integer :: ilimit
     
     real :: T_L_residual
     real :: e_a_residual
     real :: par_residual
     real :: rb_residual
     real :: prss_residual
     real :: leaf_residual
     real :: gsw_residual
     
  End Type stoma_data

end module c34constants
