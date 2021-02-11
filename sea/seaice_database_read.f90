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

subroutine seaice_database_read(iaction)

  use mem_sea,     only: sea, msea, omsea
  use mem_sfcg,     only: sfcg, itab_wsfc

  use sea_coms,    only: iupdseaice, iseaicecyclic, nseaicefiles,  &
                         fnames_seaice, ctotdate_seaice, s1900_seaice,       &
                         iseaicefile, seaice_database, iseaiceflg

  use misc_coms,   only: io6, iyear1, imonth1, idate1, itime1, timmax8,  &
                         time8, runtype, s1900_init, s1900_sim, isubdomain

  use consts_coms, only: erad, piu180
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use mem_para,    only: myrank

  implicit none

  integer, intent(in) :: iaction

  integer :: iseaicey,iseaicem,iseaiced,iseaiceh
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

! Nothing to do here if iseaiceflg is not 1

  if (iseaiceflg /= 1) return

! This subroutine is simpler than topm_database because it assumes that 
! each seaice_database file covers the entire geographic area of the model.
! If this ever changes, this subroutine must be modified.

! Check type of call to seaice_database_read

  if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

     call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize seaice cyclic flag to zero

     iseaicecyclic = 0
     nseaicefiles  = 0

     ! Open the seaice filelist and read through the files

     INQUIRE(file=seaice_database, exist=exists)
     IF (.not. exists) THEN
        WRITE(*,*) 'seaice_database_read: Error opening seaice file list ' // trim(seaice_database)
        STOP
     ENDIF

     ! Store the seaice_database directory prefix

     lns = index( seaice_database, '/', back=.true. )

     if (lns > 0) then
        dir_prefix = seaice_database(1:lns)
     else
        dir_prefix = ''
     endif

     ! Loop until the end-of file and read the file listing

     OPEN( unit=iun, file=seaice_database, status='OLD', action='READ' )

     nseaicefiles = 0
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
        ! to the seaice_database path

        if ( filename(1:1) /= '/' .and. filename(1:1) /= '~' ) then
           if ( lns > 0 ) filename = trim(dir_prefix) // trim(filename)
        endif

        nseaicefiles = nseaicefiles + 1
        nfsize    = size(fnames_tmp)
     
        ! allocate more space to hold the filenames if necessary
        if (nseaicefiles > nfsize) then
           allocate( fnames_seaice(nfsize+incr) )
           fnames_seaice(1:nfsize) = fnames_tmp
           call move_alloc( fnames_seaice, fnames_tmp )
        endif
      
        ! store the filenames
        fnames_tmp(nseaicefiles) = filename

     enddo

     close(iun)
 
     if (nseaicefiles < 1) then
        write(io6,*) 'SEAICE database files ' // trim(seaice_database) // ' were not found.'
        write(io6,*) 'Stopping run.'
        stop 'stop: no seaice files found'
     endif

     allocate( fnames_seaice  (nseaicefiles+2) )
     allocate( ctotdate_seaice(nseaicefiles+2) )
     allocate( s1900_seaice   (nseaicefiles+2) )

     fnames_seaice(1:nseaicefiles) = fnames_tmp(1:nseaicefiles)

     deallocate( fnames_tmp )

     ntimes = nseaicefiles

     do jtime = 1, ntimes

        ! Assume SEAICE file names should always end with YYYYMMDDHH.h5, and
        ! use the file name to infer the file date and time

        filename = fnames_seaice(jtime)
        slen     = len_trim(filename)
        sdate    = filename(slen-12:slen-3)

        read(sdate,'(i4,i2,i2,i2)') iseaicey, iseaicem, iseaiced, iseaiceh

        ! Convert iseaiceh format from hh to hhmmss

        iseaiceh = iseaiceh * 10000
      
        ! If file year is read as zero, seaice data is expected to be cyclic
        ! over 1 year. Increment iseaicecyclic to indicate this and use current
        ! simulation year for seaice database file times.

        if (iseaicey == 0) then
           iseaicecyclic = iseaicecyclic + 1
           iseaicey = iyears
        endif

        call date_make_big (iseaicey,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(jtime))
        call date_abs_secs2(iseaicey,iseaicem,iseaiced,iseaiceh,s1900_seaice(jtime))
     enddo

     ! Make sure files are sorted by date
     call dintsort28(nseaicefiles,ctotdate_seaice,fnames_seaice,s1900_seaice)

     ! If seaice cyclic flag > 0, check its value against ntimes and stop if they 
     ! are unequal.  If they are equal, reset seaice cyclic flag to 1 and augment
     ! seaice file arrays by 1 at each end.

     if (iseaicecyclic > 0) then

        if (iseaicecyclic /= ntimes) then
           write(io6,'(/,a)') 'Some but not all seaice database files do not have'
           write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
           stop 'stop_seaice_inv'
        endif

        iseaicecyclic = 1
        nseaicefiles = ntimes + 2
        
        ! Shift seaice data file names and times by one array element

        do jtime = ntimes,1,-1
           fnames_seaice  (jtime+1) = fnames_seaice  (jtime)
           ctotdate_seaice(jtime+1) = ctotdate_seaice(jtime)
           s1900_seaice   (jtime+1) = s1900_seaice   (jtime)
        enddo

        ! Add new seaice member at beginning of time sequence

        fnames_seaice(1) = fnames_seaice(ntimes+1)
      
        call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(ntimes+1))
        call date_make_big(iseaicey-1,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(1))
        call date_abs_secs2(iseaicey-1,iseaicem,iseaiced,iseaiceh,s1900_seaice(1))

        ! Add new seaice member at end of time sequence

        fnames_seaice(ntimes+2) = fnames_seaice(2)

        call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(2))
        call date_make_big(iseaicey+1,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(ntimes+2))
        call date_abs_secs2(iseaicey+1,iseaicem,iseaiced,iseaiceh,s1900_seaice(ntimes+2))

     endif

     ! Loop over number of SEAICE_DATABASE file times and search for the one that
     ! corresponds to current or most recent model time.

     iseaicefile = 0
     do nf = 1, nseaicefiles

        write(io6,*) 'nseaicef0 ',nf,s1900_seaice(nf),' ',s1900_sim

        if (s1900_seaice(nf) <= s1900_sim) then
           iseaicefile = nf
        endif
     enddo

     if (iseaicefile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find previous or current seaice file for current'
        write(io6,*) 'model time.  Stopping model.'
        stop 'stop: no current seaice file'
     endif

  elseif (iaction == 1) then

     ! Processing next seaice file (only called with iaction = 1 if iupdseaice = 1)

     iseaicefile = iseaicefile + 1
   
     if (iseaicefile > nseaicefiles) then
        if(iseaicecyclic == 0)then
           write(io6,*) ' '
           write(io6,*) 'No future seaice file is available for nudging '
           write(io6,*) 'Stopping model '
           stop 'stop: no future seaice file'
        else
           iseaicefile = 3
           do jtime = 1, nseaicefiles
              call date_unmake_big(iseaicey,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(jtime))
              call date_make_big(iseaicey+1,iseaicem,iseaiced,iseaiceh,ctotdate_seaice(jtime))
              call date_abs_secs2(iseaicey+1,iseaicem,iseaiced,iseaiceh,s1900_seaice(jtime))
           enddo
        endif
     endif

     sea%seaicep(:) = sea%seaicef(:)   

  endif

  ! Open and read seaice_database file

  write(io6,*) 'seaice_database_read2 ', iseaicefile, trim(fnames_seaice(iseaicefile))

  call shdf5_open(fnames_seaice(iseaicefile),'R')
  call shdf5_info('ice',ndims,idims)
  nio = idims(1)
  njo = idims(2)
  allocate(dato(nio,njo))
  call shdf5_irec(ndims,idims,'ice',rvar2=dato)
  call shdf5_close()

! If data size is even, assume that it is offset 1/2 delta lat-lon from
! the lower left point, and that we have added an extra cyclic row and
! column around the border. This is the format of the climatalogical
! and Reynolds SEAICEs

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

  ! Fill seaice array

  do isea = 2, msea
     iwsfc = isea + omsea
   
     ! Skip this cell if running in parallel and cell rank is not MYRANK
     if (isubdomain == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

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

     sea%seaicef(isea) =   &
          + wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2))  &
          + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

     sea%seaicef(isea) = max( min( sea%seaicef(isea), 1.0 ), 0.0 )

  enddo

  deallocate(dato)

  if (iaction == 0) then
     sea%seaicep(:) = sea%seaicef(:)
  endif

end subroutine seaice_database_read
