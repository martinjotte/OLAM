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
Module mem_ed

! Suppose you want to run ED only at a few scattered test sites.  These are called sites of interest (SOIs).  max_soi is the maximum allowable.
integer, parameter :: max_soi = 10
integer :: n_soi  ! n_soi is the actual number of SOIs from OLAMIN.
real, dimension(max_soi) :: soi_lat  ! latitude (-90:90) of each SOI
real, dimension(max_soi) :: soi_lon  ! longitude (-180:180) of each SOI

! Suppose you want to run ED for a few rectangular regions.  You can run for a global region, if you want.  You can also run several sub-global regions and several SOIs.  The maximum number of regions is max_ed_regions.
integer, parameter :: max_ed_regions = 10
integer :: n_ed_region  ! actual number of ED regions from OLAMIN.
real, dimension(max_ed_regions) :: ed_reg_latmin ! minimum latitude (-90:90) of each ED region. 
real, dimension(max_ed_regions) :: ed_reg_latmax ! maximum latitude (-90:90) of each ED region. 
real, dimension(max_ed_regions) :: ed_reg_lonmin ! minimum longitude (-180:180) of each ED region. 
real, dimension(max_ed_regions) :: ed_reg_lonmax ! maximum longitude (-180:180) of each ED region. 

integer, parameter :: n_pft = 11 ! number of plant functional types:
     ! 1   C4 grass
     ! 2   early successional broadleaf evergreen
     ! 3   mid successional broadleaf evergreen
     ! 4   late successional broadleaf evergreen
     ! 5   C3 grass
     ! 6   northern pines
     ! 7   southern pines
     ! 8   late successional conifers
     ! 9   early successional broadleaf deciduous
     ! 10  early successional broadleaf deciduous
     ! 11  early successional broadleaf deciduous

integer, parameter :: n_dbh = 11 ! Number of DBH bins for output quantities

end Module mem_ed
