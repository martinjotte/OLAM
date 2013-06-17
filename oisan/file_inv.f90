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
subroutine isan_file_inv()

  use max_dims,  only: pathlen, maxisdirs
  use misc_coms, only: io6, iyear1, imonth1, idate1, itime1
  use isan_coms, only: iapr, fnames_fg, s1900_fg, ctotdate_fg, nfgfiles

  implicit none

  integer :: n, nf, lnf, lns, istat, nfsize
  integer :: iyear, imonth, idate, ihour
  logical :: exists

  character(pathlen), allocatable :: fnames_tmp(:)
  character(pathlen)              :: filename

  integer, parameter :: iun = 98
  integer, parameter :: incr = 500

! Go through first guess files and make inventory

  nfgfiles = 0
  filename = ""
  allocate( fnames_tmp(incr) )

! Loop through all of the isan file listings given in the namelist

  do n = 1, maxisdirs

     if (len_trim(iapr(n)) < 1) exit

     ! Open the current file listing and read through the list of files

     inquire(file=iapr(n), exist=exists)
     if (.not. exists) then
        write(*,*) "file_inv: Error opening file list " // trim(iapr(n))
        stop       "File listing first-guess fields cannot be found"
     endif
     
     open(unit=iun, file=iapr(n), status="OLD", action="READ")
   
   ! Loop until the end-of file

   do while (.true.)

      read(iun, '(a)', iostat=istat) filename

      ! iostat /= 0 indicates end-of-file or error so exit loop
      IF (istat /= 0) exit
      
      ! Skip over any blank lines in filelist
      IF (len_trim(filename) < 1) cycle

      nfgfiles = nfgfiles + 1
      nfsize   = size(fnames_tmp)

      ! allocate more space to hold the filenames if necessary
      if (nfgfiles > nfsize) then
         allocate( fnames_fg(nfsize+incr) )
         fnames_fg(1:nfsize) = fnames_tmp
         call move_alloc( fnames_fg, fnames_tmp )
      endif

      ! store the filenames
      fnames_tmp(nfgfiles) = filename

   enddo

   close(iun)

enddo

allocate( fnames_fg  (nfgfiles) )
allocate( ctotdate_fg(nfgfiles) )
allocate( s1900_fg   (nfgfiles) )

fnames_fg(1:nfgfiles) = fnames_tmp(1:nfgfiles)

deallocate( fnames_tmp )

do nf = 1, nfgfiles

   ! strip any directory prefix from the filename

   lnf = len_trim( fnames_fg(nf) )
   lns = index   ( fnames_fg(nf), "/", back=.true.)

   filename = fnames_fg(nf)(lns+1:lnf)

   ! Starting from the right, check for a "." in the file name, and assume
   ! this is the file suffix if it is within 10 characters of the end of the
   ! filename (the suffix should only be .h5 or .hdf5)

   lnf = len_trim( filename )
   lns = index   ( filename, ".", back=.true.)

   if (lns > 0 .and. lnf-lns < 10) then
      filename = filename(1:lns-1)
      lnf      = lns - 1
   endif

   ! Some degribbed filenames have a grid number (such as -g2) appended.
   ! Strip this out if present and if it is at the right

   lns = index( filename, "-g", back=.true.)

   if (lns > 0 .and. lnf-lns < 5) then
      filename = filename(1:lns-1)
      lnf      = lns - 1
   endif

   ! With any suffix removed, assume files have the form *2005-09-01-hhmm 
   ! or *2005-09-01-hhmmss so check if seconds are present by counting how many
   ! characters are present after the last '-' character.

   lns = index( filename, "-", back=.true.)

   if (lnf-lns == 4) then
      
      ! time is in hhmm
      read (filename(lnf-14:lnf), '(i4,1x,i2,1x,i2,1x,i4)') iyear, imonth, idate, ihour

   elseif (lnf-lns == 6) then

      ! time is in hhmmss
      read (filename(lnf-16:lnf-2), '(i4,1x,i2,1x,i2,1x,i4)') iyear, imonth, idate, ihour

   else

      write(*,*) "file_inv: Error opening file:"
      write(*,*) trim(fnames_fg(nf))
      write(*,*) "Date/time portion of filename should be YYYY-MM-DD-hhmm or YYYY-MM-DD-hhmmss"
      stop       "Invalid date/time in analysis filename."

   endif

   ! Note that subroutines date_make_big and date_abs_secs2 assume that
   ! the format of ihour is hhmmss

   call date_make_big (iyear,imonth,idate,ihour*100,ctotdate_fg(nf))
   call date_abs_secs2(iyear,imonth,idate,ihour*100,s1900_fg(nf))

enddo

call dintsort28(nfgfiles,ctotdate_fg,fnames_fg,s1900_fg)

write(io6,'(/A,I0,A)') " Found ", nfgfiles, " analysis files:"
do nf = 1, nfgfiles
   write(io6,*) nf, ctotdate_fg(nf), " ", trim(fnames_fg(nf)) !, s1900_fg(nf)
enddo

end subroutine isan_file_inv
