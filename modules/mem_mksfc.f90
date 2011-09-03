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
Module mem_mksfc
  
  use max_dims, only: maxremote, maxnlspoly
  implicit none

  private :: maxremote, maxnlspoly

! DATA STRUCTURES TO HOLD THE MKSFC GRID INFORMATION FOR ALL LAND/SEA CELLS

  Type itab_mls_vars
     integer :: imglobe = 1
  End type itab_mls_vars


  Type itab_uls_vars
     integer :: im(2)   =  1
     integer :: iw(2)   =  1
     integer :: irank   = -1
     integer :: iuglobe =  1
  End type itab_uls_vars


  Type itab_wls_vars
     integer :: im(maxnlspoly)  =  1
     integer :: iu(maxnlspoly)  =  1
     integer :: irank           = -1
     integer :: iwglobe         =  1
     integer :: npoly           =  0
     logical :: send(maxremote) = .false.
  End type itab_wls_vars

End Module mem_mksfc
