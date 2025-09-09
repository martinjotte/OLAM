subroutine isan_file_inv()

  use max_dims,  only: pathlen, maxisdirs
  use misc_coms, only: io6
  use isan_coms, only: iapr, fnames_fg, s1900_fg, ctotdate_fg, nfgfiles

  implicit none

  integer  :: n, nf, lnf, lns, istat, nfsize, islash
  integer  :: iyear, imonth, idate, ihour, nff, ndays
  logical  :: exists, ierr

  character(pathlen), allocatable :: fnames_tmp(:)
  character(pathlen)              :: filename
  character(pathlen)              :: dir_prefix

  integer, parameter :: iun = 98
  integer, parameter :: incr = 500
  integer, parameter :: days_in_month(12) = [ 31,28,31,30,31,30,31,31,30,31,30,31 ]

  ! Return if analysis file inventory was already computed
  ! TODO: seperate analysis files for atm/sea/seaice?

  if (allocated(fnames_fg)) return

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

     ! Store the seaice_database directory prefix

     islash = index( iapr(n), '/', back=.true. )

     if (islash > 0) then
        dir_prefix = iapr(n)(1:islash)
     else
        dir_prefix = ''
     endif

   ! Loop until the end-of file

   do while (.true.)

      read(iun, '(a)', iostat=istat) filename

      ! iostat /= 0 indicates end-of-file or error so exit loop
      IF (istat /= 0) exit

      ! Skip over any blank lines in filelist
      IF (len_trim(filename) < 1) cycle

      ! Remove leading blanks from filenames
      filename = adjustl(filename)

      ! If the filename is a relative path, make it relative
      ! to the text listing path

      if ( filename(1:1) /= '/' .and. filename(1:1) /= '~' ) then
         if ( islash > 0 ) filename = trim(dir_prefix) // trim(filename)
      endif

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

!fnames_fg(1:nfgfiles) = fnames_tmp(1:nfgfiles)

!deallocate( fnames_tmp )

nff = 0
do nf = 1, nfgfiles

   ! strip any directory prefix from the filename

   lnf = len_trim( fnames_tmp(nf) )
   lns = index   ( fnames_tmp(nf), "/", back=.true.)

   filename = fnames_tmp(nf)(lns+1:lnf)

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

      write(io6,*)
      write(io6,*) "file_inv: skipping file ", trim(fnames_tmp(nf))
      write(io6,*) "Date/time portion of filename should be YYYY-MM-DD-hhmm or YYYY-MM-DD-hhmmss"
      cycle

   endif

   ! Check that dates are (reasonably) valid

   ierr = .false.

   if (iyear < 1900)                ierr = .true.
   if (imonth < 1 .or. imonth > 12) ierr = .true.
   if (ihour < 0 .or. ihour > 2400) ierr = .true.

   if (.not. ierr) then
      ndays = days_in_month(imonth)

      if (imonth == 2) then
         if ( (mod(iyear,400) == 0) .or. &
              (mod(iyear,  4) == 0 .and. mod(iyear,100) /= 0) ) ndays = ndays + 1
      endif

      if (idate < 1 .or. idate > ndays) ierr = .true.
   endif

   if (ierr) then
      write(io6,*)
      write(io6,*) "file_inv: skipping file ", trim(fnames_tmp(nf))
      write(io6,'(A,i4.4,1x,i2.2,1x,i2.2,1x,i4.4)') " invalid date/time ", iyear, imonth, idate, ihour
      cycle
   endif

   nff = nff + 1

   call date_abs_secs2(iyear,imonth,idate,ihour*100,s1900_fg   (nff))
   call date_make_big (iyear,imonth,idate,ihour*100,ctotdate_fg(nff))

   fnames_fg(nff) = fnames_tmp(nf)

enddo

nfgfiles = nff

call dintsort28(nfgfiles,ctotdate_fg,fnames_fg,s1900_fg)

write(io6,'(/A,I0,A)') " Found ", nfgfiles, " analysis files:"
do nf = 1, nfgfiles
   write(io6,*) nf, ctotdate_fg(nf), " ", trim(fnames_fg(nf)) !, s1900_fg(nf)
enddo

end subroutine isan_file_inv
