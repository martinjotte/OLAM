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
subroutine sst_database_read(iaction)

use mem_sea,     only: sea, itab_ws

use sea_coms,    only: mms, mws, iupdsst, isstcyclic, nsstfiles,  &
                       fnames_sst, ctotdate_sst, s1900_sst,       &
                       isstfile, sst_database, isstflg

use misc_coms,   only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                       time8, runtype, s1900_init, s1900_sim

use consts_coms, only: erad, piu180, t00
use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
use max_dims,    only: maxsstfiles

implicit none

integer, intent(in) :: iaction

integer :: issty,isstm,isstd,issth
integer :: iyears, imonths, idates, ihours

real, allocatable :: dato(:,:)

integer :: nio, njo, nperdeg
integer :: io1, io2, jo1, jo2
integer :: nf
integer :: iws
integer :: ntimes, jtime
integer :: slen
integer :: ndims, idims(2)

real :: wio1, wio2, wjo1, wjo2
real :: offpix
real :: glat, glon
real :: rio, rjo

character(len=128) :: flnm
character(len=10)  :: sdate

logical :: nocall

! This subroutine is simpler than topm_database because it assumes that 
! each sst_database file covers the entire geographic area of the model.
! If this ever changes, this subroutine must be modified.

! Check type of call to sst_database_read

if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

   call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize sst cyclic flag to zero

   isstcyclic = 0
   nsstfiles  = 0

   ! new schema to avoid call systems
   IF (isstflg == 0) THEN
      flnm = TRIM(sst_database)
      nocall = .TRUE.
      isstflg = 1
   ELSE
      flnm = TRIM(sst_database)//'??????????.h5'
      nocall = .FALSE.
   ENDIF

   write(io6,*) 'Checking for sst database files'
   CALL OLAM_filelist(fnames_sst, maxsstfiles, flnm, nsstfiles, nocall)

   if (nsstfiles < 1) then
      write(io6,*) 'SST database files '//flnm//' were not found.'
      write(io6,*) 'Stopping run.'
      stop 'stop: no sea file'
   endif
   
   ntimes = nsstfiles

   do jtime=1,ntimes

      ! Assume SST file names should always end with YYYYMMDDHH.h5, and
      ! use the file name to infer the file date and time

      flnm  = fnames_sst(jtime)
      slen = len_trim(flnm)
      sdate = flnm(slen-12:slen-3)

      read(sdate,'(i4,i2,i2,i2)') issty, isstm, isstd, issth
      
      ! If file year is read as zero, sst data is expected to be cyclic
      ! over 1 year. Increment isstcyclic to indicate this and use current
      ! simulation year for sst database file times.

      if (issty == 0) then
         isstcyclic = isstcyclic + 1
         issty = iyears
      endif

      call date_make_big (issty,isstm,isstd,issth*100,ctotdate_sst(jtime))
      call date_abs_secs2(issty,isstm,isstd,issth*100,s1900_sst(jtime))
   enddo

   ! Make sure files are sorted by date
   call dintsort28(nsstfiles,ctotdate_sst,fnames_sst,s1900_sst)

! If sst cyclic flag > 0, check its value against ntimes and stop if they 
! are unequal.  If they are equal, reset sst cyclic flag to 1 and augment
! sst file arrays by 1 at each end.

   if (isstcyclic > 0) then
      if (isstcyclic /= ntimes) then
         write(io6,'(/,a)') 'Some but not all sst database files do not have'
         write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
         stop 'stop_sst_inv'
      endif
      isstcyclic = 1
      nsstfiles = ntimes + 2

! Shift sst data file names and times by one array element

      do jtime = ntimes,1,-1
         fnames_sst  (jtime+1) = fnames_sst  (jtime)
         ctotdate_sst(jtime+1) = ctotdate_sst(jtime)
         s1900_sst   (jtime+1) = s1900_sst   (jtime)
      enddo      

! Add new sst member at beginning of time sequence

      fnames_sst(1) = fnames_sst(ntimes+1)

      call date_unmake_big(issty,isstm,isstd,issth,ctotdate_sst(ntimes+1))
      call date_make_big(issty-1,isstm,isstd,issth,ctotdate_sst(1))
      call date_abs_secs2(issty-1,isstm,isstd,issth,s1900_sst(1))

! Add new sst member at end of time sequence

      fnames_sst(ntimes+2) = fnames_sst(2)

      call date_unmake_big(issty,isstm,isstd,issth,ctotdate_sst(2))
      call date_make_big(issty+1,isstm,isstd,issth,ctotdate_sst(ntimes+2))
      call date_abs_secs2(issty+1,isstm,isstd,issth,s1900_sst(ntimes+2))

   endif

! Loop over number of SST_DATABASE file times and search for the one that
! corresponds to current or most recent model time.

   isstfile = 0
   do nf = 1,nsstfiles

      write(io6,*) 'nsstf0 ',nf,s1900_sst(nf),' ',s1900_sim

      if (s1900_sst(nf) <= s1900_sim) then
         isstfile = nf
      endif
   enddo

   if (isstfile < 1) then
      write(io6,*) ' '
      write(io6,*) 'Unable to find previous or current sst file for current'
      write(io6,*) 'model time.  Stopping model.'
      stop 'stop: no current sst file'
   endif

elseif (iaction == 1) then

! Processing next sst file (only called with iaction = 1 if iupdsst = 1)

   isstfile = isstfile + 1
   
   if (isstfile > nsstfiles) then
      if(isstcyclic == 0)then
         write(io6,*) ' '
         write(io6,*) 'No future sst file is available for nudging '
         write(io6,*) 'Stopping model '
         stop 'stop: no future sst file'
      else
         isstfile = 3
         do jtime = 1, nsstfiles
            call date_unmake_big(issty,isstm,isstd,issth,ctotdate_sst(jtime))
            call date_make_big(issty+1,isstm,isstd,issth,ctotdate_sst(jtime))
            call date_abs_secs2(issty+1,isstm,isstd,issth,s1900_sst(jtime))
         enddo
      endif
   endif

   sea%seatp(:) = sea%seatf(:)   

endif

! Open and read sst_database file

write(io6,*) 'sst_database_read2 ', isstfile, trim(fnames_sst(isstfile))

call shdf5_open(fnames_sst(isstfile),'R')
call shdf5_info('sst',ndims,idims)
nio = idims(1)
njo = idims(2)
allocate(dato(nio,njo))
call shdf5_irec(ndims,idims,'sst',rvara=dato)
call shdf5_close()

offpix = 0.
nperdeg = nio/360
if (mod(nio,nperdeg) == 2) offpix = .5

! Fill sst array

do iws = 2,mws
   
   glat = asin(sea%zew(iws) / erad) * piu180
   glon = atan2(sea%yew(iws), sea%xew(iws)) * piu180

   glon = max(-179.999,min(179.999,glon))

   rio = 1. + (glon + 180.) * nperdeg + offpix
   rjo = 1. + (glat +  90.) * nperdeg + offpix

   io1 = int(rio)
   jo1 = int(rjo)
         
   wio2 = rio - float(io1)
   wjo2 = rjo - float(jo1)
           
   wio1 = 1. - wio2
   wjo1 = 1. - wjo2

   io2 = io1 + 1
   jo2 = jo1 + 1

   sea%seatf(iws) = t00  &
        + wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2))  &
        + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

enddo

deallocate(dato)

if (iaction == 0) then
   sea%seatp(:) = sea%seatf(:)
endif

return
end subroutine sst_database_read
