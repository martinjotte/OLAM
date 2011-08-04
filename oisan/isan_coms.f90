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
Module isan_coms

use max_dims,    only: maxisfiles, maxpr, maxisdirs, pathlen
use consts_coms, only: r8

integer :: ioflgisz, ioflgvar, iszstage, ivrstage
integer :: iyear, imonth, idate, ihour, ipoffset
integer :: npd,lzon_bot,kzonoff,isdirs

real, dimension(maxpr+2) :: pcol_p, pcol_thet, pcol_pk, pcol_u, pcol_v, pcol_z
real, dimension(maxpr+2) :: pcol_r, pcol_pi, pcol_temp, pcol_rt, pcol_thv

integer :: nfgfiles, ifgfile

character(pathlen) :: iapr(maxisdirs)
character(pathlen) :: fnames_fg(maxisfiles)
character(pathlen) :: innpr
character(14)      :: ctotdate_fg(maxisfiles)

real(r8)  :: s1900_fg(maxisfiles)

! Pressure header variables:
integer :: marker, isversion, iyy, imm, idd, ihh, itinc, inproj, ivertcoord
real    :: xnelat, xnelon, cntlat, cntlon, secondlat

integer :: nprx, npry, nprz, nprz_rh
integer :: levpr(maxpr)

real    :: xswlon, xswlat, gdatdx, gdatdy
real    :: pnpr(maxpr)

End module isan_coms
