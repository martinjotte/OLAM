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
Module mem_mksfc
  
  use max_dims, only: maxnlspoly

  implicit none

! DATA STRUCTURES TO HOLD THE MKSFC GRID INFORMATION FOR ALL LAND/SEA CELLS

  Type itab_wls_vars
     integer :: irank = -1
     integer :: iwglobe = 1  ! = iwlglobe or iwsglobe
     integer :: iw = 1
     integer :: kw = 1
     integer :: npoly = 0

     real :: arf_iw = 0.
     real :: arf_kw = 0.

     real :: xem(maxnlspoly) = 1
     real :: yem(maxnlspoly) = 1
     real :: zem(maxnlspoly) = 1
  End type itab_wls_vars

  Type landsea_grid_vars
      real :: area  = 0. ! cell surface area [m^2]
      real :: glatw = 0. ! latitude of sfc cell W points
      real :: glonw = 0. ! longitude of sfc cell W points
      real :: xew   = 0. ! earth x coord of sfc W points
      real :: yew   = 0. ! earth y coord of sfc W points
      real :: zew   = 0. ! earth z coord of sfc W points
      real :: topw  = 0. ! topographic height of sfc W points
      real :: wnx   = 0. ! norm unit vector x comp of sfc cells
      real :: wny   = 0. ! norm unit vector y comp of sfc cells
      real :: wnz   = 0. ! norm unit vector z comp of sfc cells

      integer :: leaf_class = 0 ! leaf ("vegetation") class

      integer :: wpt_sea  = 0 ! Counter of W pt for sea cell
      integer :: wpt_land = 0 ! Counter of W pt for land cell
      integer :: idatq    = 0 ! integer array for storing database data
  End type landsea_grid_vars

End Module mem_mksfc
