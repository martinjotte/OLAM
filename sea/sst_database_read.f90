subroutine sst_database_read(iaction)

  use mem_sea,     only: sea, msea, omsea
  use mem_sfcg,    only: sfcg, itab_wsfc

  use sea_coms,    only: isstcyclic, nsstfiles, fnames_sst, &
                         ctotdate_sst, s1900_sst,  &
                         isstfile, sst_database, isstflg

  use misc_coms,   only: io6, s1900_sim, iparallel
  use consts_coms, only: t00
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use mem_para,    only: myrank

  implicit none

  integer, intent(in) :: iaction

  integer :: issty,isstm,isstd,issth
  integer :: iyears, imonths, idates, ihours

  real, allocatable :: dato(:,:)

  integer :: nio, njo
  integer :: nx, ny
  integer :: io1, io2, jo1, jo2
  integer :: nf, nfsize
  integer :: lns, istat
  integer :: isea, iwsfc
  integer :: ntimes, jtime
  integer :: slen
  integer :: ndims, idims(2)

  real :: wio1, wio2, wjo1, wjo2
  real :: xoffpix, yoffpix
  real :: glat, glon
  real :: rio, rjo
  real :: xperdeg, yperdeg

  character(pathlen), allocatable :: fnames_tmp(:)
  character(pathlen)              :: filename
  character(pathlen)              :: dir_prefix
  character(10)                   :: sdate

  integer, parameter :: iun = 98
  integer, parameter :: incr = 500

  logical :: exists

! Nothing to do here if isstflg is not 1

  if (isstflg /= 1) return

! This subroutine is simpler than land_database_read because it assumes that
! each sst_database file covers the entire geographic area of the model.
! If this ever changes, this subroutine must be modified.

! Check type of call to sst_database_read

  if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

     call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize sst cyclic flag to zero

     isstcyclic = 0
     nsstfiles  = 0

     ! Open the sst filelist and read through the files

     INQUIRE(file=sst_database, exist=exists)
     IF (.not. exists) THEN
        WRITE(*,*) 'sst_database_read: Error opening sst file list ' // trim(sst_database)
        STOP
     ENDIF

     ! Store the sst_database directory prefix

     lns = index( sst_database, '/', back=.true. )

     if (lns > 0) then
        dir_prefix = sst_database(1:lns)
     else
        dir_prefix = ''
     endif

     ! Loop until the end-of file and read the file listing

     OPEN( unit=iun, file=sst_database, status='OLD', action='READ' )

     nsstfiles = 0
     filename  = ''
     allocate( fnames_tmp(incr) )

     do while (.true.)

        read(iun, '(a)', iostat=istat) filename

        ! iostat /= 0 indicates end-of-file or error so exit loop
        IF (istat /= 0) exit

        ! Skip over any blank lines in filelist
        IF (len_trim(filename) < 1) cycle

        ! Remove leading blanks from filenames
        filename = adjustl(filename)

        ! If the filename is a relative path, make it relative
        ! to the sst_database path

        if ( filename(1:1) /= '/' .and. filename(1:1) /= '~' ) then
           if ( lns > 0 ) filename = trim(dir_prefix) // trim(filename)
        endif

        nsstfiles = nsstfiles + 1
        nfsize    = size(fnames_tmp)

        ! allocate more space to hold the filenames if necessary
        if (nsstfiles > nfsize) then
           allocate( fnames_sst(nfsize+incr) )
           fnames_sst(1:nfsize) = fnames_tmp
           call move_alloc( fnames_sst, fnames_tmp )
        endif

        ! store the filenames
        fnames_tmp(nsstfiles) = filename

     enddo

     close(iun)

     if (nsstfiles < 1) then
        write(io6,*) 'SST database files ' // trim(sst_database) // ' were not found.'
        write(io6,*) 'Stopping run.'
        stop 'stop: no sst files found'
     endif

     allocate( fnames_sst  (nsstfiles+2) )
     allocate( ctotdate_sst(nsstfiles+2) )
     allocate( s1900_sst   (nsstfiles+2) )

     fnames_sst(1:nsstfiles) = fnames_tmp(1:nsstfiles)

     deallocate( fnames_tmp )

     ntimes = nsstfiles

     do jtime = 1, ntimes

        ! Assume SST file names should always end with YYYYMMDDHH.h5, and
        ! use the file name to infer the file date and time

        filename = fnames_sst(jtime)
        slen     = len_trim(filename)
        sdate    = filename(slen-12:slen-3)

        read(sdate,'(i4,i2,i2,i2)') issty, isstm, isstd, issth

        ! Convert issth format from hh to hhmmss

        issth = issth * 10000

        ! If file year is read as zero, sst data is expected to be cyclic
        ! over 1 year. Increment isstcyclic to indicate this and use current
        ! simulation year for sst database file times.

        if (issty == 0) then
           isstcyclic = isstcyclic + 1
           issty = iyears
        endif

        call date_make_big (issty,isstm,isstd,issth,ctotdate_sst(jtime))
        call date_abs_secs2(issty,isstm,isstd,issth,s1900_sst(jtime))
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
     do nf = 1, nsstfiles

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

  call shdf5_open(fnames_sst(isstfile),'R',trypario=.true.)
  call shdf5_info('sst',ndims,idims)
  nio = idims(1)
  njo = idims(2)
  allocate(dato(nio,njo))
  call shdf5_irec(ndims,idims,'sst',rvar2=dato)
  call shdf5_close()

! If data size is even, assume that it is offset 1/2 delta lat-lon from
! the lower left point, and that we have added an extra cyclic row and
! column around the border. This is the format of the climatalogical
! and Reynolds SSTs

! If data size is odd, assume that it is "unstaggered" like grib datasets,
! so that the SW point is at (-180, -90) and the NW point is at (+180, +90)
! with no extra rows/columns

  if (mod(nio,2) == 0) then
     xoffpix = 0.5
     nx = nio - 2
  else
     xoffpix = 0.0
     nx = nio - 1
  endif

  if (mod(njo,2) == 0) then
     yoffpix = 0.5
     ny = njo - 2
  else
     yoffpix = 0.0
     ny = njo - 1
  endif

  xperdeg = real(nx) / 360.0
  yperdeg = real(ny) / 180.0

  ! Fill sst array

  do isea = 2, msea
     iwsfc = isea + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     glat = sfcg%glatw(iwsfc)
     glon = sfcg%glonw(iwsfc)

     glon = max(-179.999,min(179.999,glon))

     rio = 1. + (glon + 180.) * xperdeg + xoffpix
     rjo = 1. + (glat +  90.) * yperdeg + yoffpix

     io1 = int(rio)
     jo1 = int(rjo)

     wio2 = rio - real(io1)
     wjo2 = rjo - real(jo1)

     wio1 = 1. - wio2
     wjo1 = 1. - wjo2

     io2 = min(nio, io1 + 1)
     jo2 = min(njo, jo1 + 1)

     sea%seatf(isea) = t00  &
          + wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2))  &
          + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

  enddo

  deallocate(dato)

  if (iaction == 0) then
     sea%seatp(:) = sea%seatf(:)
  endif

end subroutine sst_database_read
