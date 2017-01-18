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
Module misc_coms

use max_dims,    only: maxsndg, maxgrds, maxngrdll, pathlen
use consts_coms, only: r8
implicit none

private :: maxsndg, maxgrds, maxngrdll, pathlen, r8

type simtime
   integer  :: year
   integer  :: month
   integer  :: date
   real(r8) :: time
end type simtime

type(simtime) :: current_time ! current simulation time

character(pathlen) :: tmpdir = "/tmp"
character(64)      :: expnme
character(16)      :: runtype
character(1)       :: timeunit
character(pathlen) :: gridfile
character(pathlen) :: hfilin
character(pathlen) :: hfilepref
character(pathlen) :: zonclim
character(pathlen) :: topo_database

real     :: rinit  = 0.0
real(r8) :: rinit8 = 0.0_r8

logical :: debug_fp  = .false.
logical :: init_nans = .false.

integer :: io6
integer :: initial
integer :: ngrids
integer :: ngrids_old
integer :: nzp
integer :: ilwrtyp
integer :: iswrtyp
integer :: icfrac
integer :: mdomain
integer :: nsndg
integer :: iflag
integer :: iyear1
integer :: imonth1
integer :: idate1
integer :: ihour1
integer :: itime1
integer :: naddsc
integer :: ipsflg
integer :: itsflg
integer :: irtsflg
integer :: iusflg
integer :: ioutput
integer :: iclobber
integer :: itopoflg
integer :: nzpp
integer :: nscl
integer :: nxp
integer :: iparallel = 0
integer :: isubdomain = 0
integer :: ndz
integer :: mstp

integer :: idiffk   (maxgrds)
integer :: ndtrat   (maxgrds)
integer :: nacoust  (maxgrds)
integer :: nqparm   (maxgrds)

real(r8) :: dtlong
real(r8) :: frqstate
real(r8) :: radfrq
real(r8) :: confrq

real :: ubmin
real :: topref
real :: polelat
real :: polelon
real :: deltax
real :: hdz(10)
real :: dz(10)
real :: p_sfc
real :: cfracrh1
real :: cfracrh2
real :: cfraccup

integer :: ngrdll(maxgrds)

real :: grdrad(maxgrds)
real :: grdlat(maxgrds,maxngrdll)
real :: grdlon(maxgrds,maxngrdll)
real :: cflxy (maxgrds)
real :: cflz  (maxgrds)
real :: csz   (maxgrds)
real :: csx   (maxgrds)
real :: akmin (maxgrds)

real(r8) :: dtlm(maxgrds)
real(r8) :: dtsm(maxgrds)

real, allocatable :: u01d (:)
real, allocatable :: v01d (:)
real, allocatable :: pr01d(:)
real, allocatable :: th01d(:)
real, allocatable :: dn01d(:)
real, allocatable :: rt01d(:)

real :: us  (maxsndg)
real :: vs  (maxsndg)
real :: ts  (maxsndg)
real :: thds(maxsndg)
real :: ps  (maxsndg)
real :: hs  (maxsndg)
real :: rts (maxsndg)

real(r8) :: time8
real(r8) :: time8p
real(r8) :: time_istp8
real(r8) :: time_istp8p
real(r8) :: timmax8
real(r8) :: s1900_init
real(r8) :: s1900_sim
real(r8) :: time_prevhist = 0.0_r8
real(r8) :: time_bias ! A number small compared to the timestep

integer       :: do_chem   =  0

Contains

!===============================================================================

   subroutine alloc_misc(mza)

   implicit none
   
   integer, intent(in) :: mza

   allocate (u01d (mza))
   allocate (v01d (mza))
   allocate (pr01d(mza))
   allocate (th01d(mza))
   allocate (dn01d(mza))
   allocate (rt01d(mza))
   
   return
   end subroutine alloc_misc
   
End Module


