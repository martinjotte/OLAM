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
Module isan_coms

  use max_dims,    only: maxisdirs, pathlen
  use consts_coms, only: r8

  implicit none

  integer :: ioflgisz, ioflgvar, iszstage, ivrstage
  integer :: iyear, imonth, idate, ihour, ipoffset
  integer :: npd,lzon_bot,kzonoff

  integer :: nfgfiles = -1
  integer :: ifgfile  =  0

  character(pathlen) :: iapr(maxisdirs)
  character(pathlen) :: innpr

  character(pathlen), allocatable :: fnames_fg  (:)
  character(14),      allocatable :: ctotdate_fg(:)
  real(r8),           allocatable :: s1900_fg   (:)

  logical :: haso3
  integer :: nbot_o3

  ! Only used on MPI node 0:
  character(6) :: o3name = " "

! Pressure header variables:
  integer :: isversion, iyy, imm, idd, ihh, itinc, inproj, ivertcoord
  real    :: xnelat, xnelon, cntlat, secondlat

  integer :: nprx, npry, nprz
  integer, allocatable :: levpr(:)

  real    :: xswlon, xswlat, gdatdx, gdatdy
  real, allocatable :: pnpr(:)

  real, allocatable :: glat(:), plat(:)

  real, allocatable :: pcol_p(:)
  real, allocatable :: pcol_z(:,:)

  ! input variables interpolated to height at each OLAM column
  real, allocatable :: o_press (:,:)
  real, allocatable :: o_rho   (:,:)

  real, allocatable :: o_theta (:,:)
  real, allocatable :: o_rrw   (:,:)
  real, allocatable :: o_uzonal(:,:)
  real, allocatable :: o_umerid(:,:)
  real, allocatable :: o_ozone (:,:)

End module isan_coms
