SUBROUTINE OLAM_filelist(fnames,nfnam,file_prefix,nfile,nocall)

  use mem_para,  only: myrank
  use misc_coms, only: io6, tmpdir
  use max_dims,  only: pathlen
  implicit none

  integer, intent(in)  :: nfnam
  integer, intent(out) :: nfile

  character(len=*), intent(in)  :: file_prefix
  character(len=*), intent(out) :: fnames(nfnam)

  LOGICAL, intent(in) :: nocall

  character(len=256) :: command
  character(pathlen) :: file, cdir, dirname
  character(pathlen) :: tmpname, fname, dirtail
  character(len=10)  :: zrank
  character(len=30)  :: tmpfile

  integer, parameter :: iun=98
  integer            :: iprelen,nc,nf,lndir,slen,istat
  logical            :: exists

  IF (nocall) THEN

     ! file_prefix is a list with all files, identical to the call system
     ! with the find command. So tmpname is file_prefix in this case.
     tmpname = file_prefix

  ELSE

     WRITE(io6,*)
     WRITE(io6,*) 'OLAM_filelist: Checking prefix: '//TRIM(file_prefix)

     nfile = 0
     iprelen = LEN_TRIM(file_prefix)

     lndir = SCAN(file_prefix, '/', back=.TRUE.)
     IF (lndir > 0) THEN
        dirname = file_prefix(1:lndir-1)
        fname =   file_prefix(lndir+1:iprelen)
        slen =    SCAN(dirname, '/', back=.TRUE.) + 1
        dirtail = dirname(slen:)
     ELSE
        dirname = "."
        dirtail = "."
        fname = file_prefix(1:iprelen)
     ENDIF

     ! Determine a unique temporary filename

     WRITE(zrank,'(I0)') myrank
     tmpfile = 'olam_' // trim(zrank) // '.XXXXXXXX'
     CALL mktempname(tmpfile)
     tmpname = trim(tmpdir) // '/' // trim(tmpfile)

     ! Make sure we can write to the tmp directory

     OPEN(unit=iun, file=tmpname, action='write', iostat=istat)
     CLOSE(iun, status='delete')

     IF (istat /= 0) THEN
        WRITE(*,*)
        WRITE(*,*) "Error writing to tmp directory " // trim(tmpdir)
        STOP       "OLAM_filelist: Please set $TMPDIR to the temporary directory"
     ENDIF

     ! The 'ls' command only takes a limited number of file names as
     ! arguments. Instead, we use the 'find' command to locate files and
     ! the '-prune' option to not recursively decend into sub-directories

     command = 'find '//TRIM(dirname)//'/. "(" -type d -a ! -name . -prune ")" &
                -o -name "'//TRIM(fname)//'" -print > '//TRIM(tmpname)
     CALL system(command)

  ENDIF

  ! Open the directory list and read through the files

  INQUIRE(file=tmpname, exist=exists)
  IF (.not. exists) THEN
     WRITE(*,*) 'OLAM_filelist: Error opening temporary OLAM_filelist'
     STOP 'OLAM_filelist: /tmp file error. Run again.'
  ENDIF

  OPEN(unit=iun, file=tmpname, status='old')

  nfile = 0
  DO nf = 1, nfnam+10
     READ(iun, '(a)', iostat=istat) file

     ! iostat /= 0 indicates end-of-file or error so exit loop
     IF (istat /= 0) exit

     ! Skip over any blank lines in filelist
     IF (len_trim(file) < 1) cycle

     nfile = nfile + 1

     IF (nfile > nfnam) THEN
        WRITE(*,*) 'OLAM_filelist: too many files of the form:'
        WRITE(*,*) trim(file_prefix)
        STOP       'Please increase the appropriate parameter in max_dims'
     ENDIF

     fnames(nfile) = file
  ENDDO

  if (nocall) then
     close(iun)
  else
     close(iun, status='delete')
  endif

END SUBROUTINE OLAM_filelist


!===============================================================================


subroutine mktempname(tempname)

! FORTRAN ROUTINE TO MIMIC THE BEHAVIOUR OF THE MKTEMP SYSTEM FUNCTION
! THE INPUT IS A STRING TEMPLATE REPRESENTING A FILE NAME SUCH AS
! /tmp/temp.XXXX.
! ON OUTPUT, THE X'S ARE REPLACED BY RANDOM COMBINATION OF LETTERS AND NUMBERS

  implicit none
  character(len=*), intent(inout) :: tempname

  integer :: date_time(8)
  integer :: seed_size, s, i
  integer, allocatable :: seed(:)
  real :: z

! IF tempname DOESN'T CONTAIN ANY X'S, DO NOTHING

  if (scan(trim(tempname),'X') .eq. 0) return

! INITIALIZE THE RANDOM NUMBER GENERATOR WITH A SEED BASED
! ON THE CURRENT DATA AND TIME

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 0

  call date_and_time(values=date_time)

  do i=1,min(8,seed_size)
     seed(i) = date_time(8-i+1)
  enddo

  call random_seed(put=seed)
  deallocate(seed)

! NOW LOOP OVER EACH LETTER IN tempname, REPLACING EACH X WITH A RANDOM
! LETTER OR NUMBER

  do i=1,len_trim(tempname)
     if (tempname(i:i).eq.'X') then

        call random_number(z) ! 0.0 < z < 1.0
        s = nint(z*61.0)      !   0 < s < 61

        if (s.le.9) then
           s = s+48           ! ASCII 0-9
        elseif (s.le.35) then
           s = s+55           ! ASCII A-Z
        else
           s = s+61           ! ASCII a-z
        endif

        ! convert ASCII -> character
        tempname(i:i) = achar(s)

     endif
  enddo

end subroutine mktempname
