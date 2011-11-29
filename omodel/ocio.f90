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

subroutine commio(action)

! THIS ROUTINE READS OR WRITES NON-HORIZONTALLY VARYING FIELDS
! THAT ARE COMMON TO ALL PROCESSES. CALL THIS ROUTINE WITH
! 'ACTION' EQUALS 'READ' OR 'WRITE' TO INPUT OR OUTPUT THE FIELDS

use max_dims,   only: maxngrdll
use misc_coms,  only: io6, itime1, idate1, imonth1, iyear1, nxp, ngrids, &
                      ngrdll, grdrad, grdlat, grdlon, nzp, &
                      mdomain, meshtype, deltax, &
                      itopoflg, time8, ndz, hdz, dz, current_time
use leaf_coms,  only: nzg, nzs, slz, ivegflg, isfcl
use hdf5_utils, only: shdf5_orec, shdf5_irec, shdf5_io

implicit none
integer          :: ndims, idims(2), k, ihour
character(len=*) :: action

if (action /= "READ" .and. action /= "WRITE") then
   write(io6,*) 'Illegal action in routine commio'
   write(io6,*) 'action must be "READ" or "WRITE"'
   stop     'Stopping run'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GROUP1: NAMELIST PARAMETERS THAT WE DON'T WANT CHANGED ON HISTORY RESTART.
!         THESE SHOULD CORRESPOND WITH THE VARIABLES IN THE NOT_HISTORY 
!         SECTION OF THE SUBROUTINE COPY_NL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndims = 1
idims(1) = 1
call shdf5_io(action, ndims, idims, 'nl%itime1',  ivars=itime1)
call shdf5_io(action, ndims, idims, 'nl%idate1' , ivars=idate1)
call shdf5_io(action, ndims, idims, 'nl%imonth1', ivars=imonth1)
call shdf5_io(action, ndims, idims, 'nl%iyear1',  ivars=iyear1)
call shdf5_io(action, ndims, idims, 'nl%ngrids',  ivars=ngrids)
call shdf5_io(action, ndims, idims, 'nl%nzp',     ivars=nzp)
call shdf5_io(action, ndims, idims, 'nl%nzg',     ivars=nzg)
call shdf5_io(action, ndims, idims, 'nl%nzs',     ivars=nzs)
call shdf5_io(action, ndims, idims, 'nl%nxp',     ivars=nxp)
call shdf5_io(action, ndims, idims, 'nl%mdomain', ivars=mdomain)
call shdf5_io(action, ndims, idims, 'nl%meshtype',ivars=meshtype)
call shdf5_io(action, ndims, idims, 'nl%isfcl',   ivars=isfcl)
call shdf5_io(action, ndims, idims, 'nl%itopoflg',ivars=itopoflg)
call shdf5_io(action, ndims, idims, 'nl%ivegflg', ivars=ivegflg)
call shdf5_io(action, ndims, idims, 'nl%deltax',  rvars=deltax)

call shdf5_io(action, ndims, idims, 'nl%ndz',     ivars=ndz)

ndims = 1
idims(1) = ndz
call shdf5_io(action, ndims, idims, 'nl%hdz',    rvara=hdz)
call shdf5_io(action, ndims, idims, 'nl%dz',     rvara=dz)

ndims = 1
idims(1) = ngrids
call shdf5_io(action, ndims, idims, 'nl%ngrdll', ivara=ngrdll)
call shdf5_io(action, ndims, idims, 'nl%grdrad', rvara=grdrad)

ndims = 2
idims(1) = ngrids
idims(2) = maxngrdll

call shdf5_io(action, ndims, idims, 'nl%grdlat', rvara=grdlat(1:ngrids,:))
call shdf5_io(action, ndims, idims, 'nl%grdlon', rvara=grdlon(1:ngrids,:))

ndims=1
idims(1) = nzg
call shdf5_io(action, ndims, idims, 'nl%slz', rvara=slz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GROUP 2: NON-HORIZONTALLY-VARYING (COMMOM) MODEL VARIABLES
!           THAT ARE NOT PART OF THE VTABLES I/O MECHANISM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndims=1
idims(1) = 1
call shdf5_io(action, ndims, idims, 'time8',     dvars=time8)
call shdf5_io(action, ndims, idims, 'cur%year',  ivars=current_time%year)
call shdf5_io(action, ndims, idims, 'cur%month', ivars=current_time%month)
call shdf5_io(action, ndims, idims, 'cur%date',  ivars=current_time%date)
call shdf5_io(action, ndims, idims, 'cur%time',  rvars=current_time%time)

return
end subroutine commio
