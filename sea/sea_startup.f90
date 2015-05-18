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
subroutine sea_startup()

use sea_coms,  only: mms, mus, mws, isstflg, seatmp,  &
                     iupdsst, iseaiceflg, iupdseaice, seaice

use mem_sea,   only: sea, alloc_sea, filltab_sea

use misc_coms, only: io6, runtype

implicit none

integer :: iws

! Subroutine SEA_STARTUP allocates some sea arrays and initializes sst

! THIS SUBROUTINE DOES NOT INITIALIZE canopy temperature and moisture
! values, which depend on atmospheric conditions.

!-------------------------------------------------------------------------------
! STEP 1: Call alloc_sea and filltab_sea (sea grid arrays already allocated)
!-------------------------------------------------------------------------------

call alloc_sea(mws)
call filltab_sea()

!-------------------------------------------------------------------------------
! STEP 2a: Fill sst values
!-------------------------------------------------------------------------------

if (isstflg == 0) then

! Default initialization of SST

   do iws = 2,mws
      sea%seatp(iws) = seatmp
      sea%seatf(iws) = seatmp
   enddo

elseif (isstflg == 1) then

   if (runtype == 'INITIAL' .or. &
       runtype == 'HISTORY' .or. &
       runtype == 'HISTADDGRID') then

! Read standard SST database
! Not needed for a plotonly run

      write(io6,'(/,a)') 'calling sst_database_read(0)'
      call sst_database_read(0)

! Future SSTs are only needed if iupdsst == 1

      if (iupdsst == 1) then
         write(io6,'(/,a)') 'calling sst_database_read(1)'
         call sst_database_read(1)
      endif

   endif

elseif (isstflg == 2) then

   if (runtype == 'INITIAL' .or. &
       runtype == 'HISTORY' .or. &
       runtype == 'HISTADDGRID') then

! Read SST from degribbed analysis files
! Not needed for a plotonly run

      write(io6,'(/,a)') 'calling read_sst_analysis(0)'
      call read_sst_analysis(0)

! Future SSTs are only needed if iupdsst == 1

      if (iupdsst == 1) then
         write(io6,'(/,a)') 'calling read_sst_analysis(1)'
         call read_sst_analysis(1)
      endif

   endif

endif

!-------------------------------------------------------------------------------
! STEP 2b: Fill sea ice values
!-------------------------------------------------------------------------------

if (iseaiceflg == 0) then

! Default initialization of SEAICE

   do iws = 2,mws
      sea%seaicep(iws) = seaice
      sea%seaicef(iws) = seaice
   enddo

elseif (iseaiceflg == 1) then

   if (runtype == 'INITIAL' .or. &
       runtype == 'HISTORY' .or. &
       runtype == 'HISTADDGRID') then

! Read standard SEAICE database
! Not needed for a plotonly run

      write(io6,'(/,a)') 'calling seaice_database_read(0)'
      call seaice_database_read(0)

! Future SEAICE is only needed if iupdseaice == 1

      if (iupdseaice == 1) then
         write(io6,'(/,a)') 'calling seaice_database_read(1)'
         call seaice_database_read(1)
      endif
   
   endif

elseif (iseaiceflg == 2) then

   if (runtype == 'INITIAL' .or. &
       runtype == 'HISTORY' .or. &
       runtype == 'HISTADDGRID') then

! Read SEAICE from degribbed analysis files
! Not needed for a plotonly run

      write(io6,'(/,a)') 'calling read_seaice_analysis(0)'
      call read_seaice_analysis(0)

! Future SEAICE is only needed if iupdseaice == 1

      if (iupdseaice == 1) then
         write(io6,'(/,a)') 'calling read_seaice_analysis(1)'
         call read_seaice_analysis(1)
      endif

   endif

endif

end subroutine sea_startup
