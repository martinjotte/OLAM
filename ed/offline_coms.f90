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
Module offline_coms

  integer, parameter :: max_ol_vars = 12
  integer :: nformats
  character(len=128), allocatable, dimension(:) :: ed_ol_names
  integer, allocatable, dimension(:) :: ed_ol_nlon
  integer, allocatable, dimension(:) :: ed_ol_nlat
  real, allocatable, dimension(:) :: ed_ol_dx
  real, allocatable, dimension(:) :: ed_ol_dy
  real, allocatable, dimension(:) :: ed_ol_xmin
  real, allocatable, dimension(:) :: ed_ol_ymin
  integer, allocatable, dimension(:) :: ed_ol_nv
  character(len=16), allocatable, dimension(:,:) :: ed_ol_vars
  real, allocatable, dimension(:,:) :: ed_ol_frq
  integer, allocatable, dimension(:,:) :: ed_ol_interp

  Type met_driv_data

     real, allocatable, dimension(:) :: nbdsf
     real, allocatable, dimension(:) :: nddsf
     real, allocatable, dimension(:) :: vbdsf
     real, allocatable, dimension(:) :: vddsf
     real, allocatable, dimension(:) :: prate
     real, allocatable, dimension(:) :: dlwrf
     real, allocatable, dimension(:) :: pres
     real, allocatable, dimension(:) :: hgt
     real, allocatable, dimension(:) :: ugrd
     real, allocatable, dimension(:) :: vgrd
     real, allocatable, dimension(:) :: sh
     real, allocatable, dimension(:) :: tmp

     real :: nir_beam
     real :: nir_diffuse
     real :: par_beam
     real :: par_diffuse
     real :: geoht
     real :: atm_tmp
     real :: atm_shv
     real :: rhos

  end Type met_driv_data

end Module offline_coms
