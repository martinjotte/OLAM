subroutine ndvi_database_read(iaction)

  use mem_land,  only: land, mland, omland
  use mem_sfcg,  only: sfcg, itab_wsfc
  use leaf_coms, only: fnames_ndvi, ctotdate_ndvi, s1900_ndvi,  &
                       indvifile, ndvi_database
  use misc_coms, only: io6, s1900_sim, iparallel
  use land_db,   only: land_database_read
  use mem_para,  only: myrank

  implicit none

  integer, intent(in) :: iaction

  integer :: indviy,indvim,indvid,indvih
  integer :: iyears, imonths, idates, ihours

  character(len=1)  :: dummy
  character(len=128):: flnm
  character(len=5)  :: filemonth

  integer :: nf
  integer :: iland, iwsfc
  integer :: ntimes, jtime
  logical :: there

  real :: datq(mland)

  integer, save :: indvicyclic, nndvifiles

  real, allocatable :: glatwl(:), glonwl(:)

! Check type of call to ndvi_database_read

  if (iaction == 0) then

! Convert current model time from s1900 to years, months, dates, hours

     call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

! Initialize ndvi cyclic flag to zero

     indvicyclic = 0

     flnm = trim(ndvi_database)//'HEADER'

     write(io6,*)  'Checking ndvi database header file ',trim(flnm)

     inquire(file=flnm,exist=there)

     if (.not. there) then
        write(io6,*) 'NDVI database file header was not found - stopping run'
        stop 'stop: no ndvi file'
     endif

! Open & read NDVI_DATABASE HEADER file; make inventory of file names & times
! Assumes data in header file is chronologically ordered oldest to newest.

     open(29, file=flnm, form='FORMATTED', status='OLD', action='READ')
     read(29,*) dummy
     read(29,*) ntimes

     nndvifiles = ntimes

     do jtime = 1,ntimes
        read(29,*) filemonth, indviy, indvim, indvid, indvih

! Convert indvih from hh to hhmmss format

        indvih = indvih * 10000

! If file year is read as zero, ndvi data is expected to be cyclic over 1 year.
! Increment indvicyclic to indicate this and use current simulation year for
! ndvi database file times.

        if (indviy == 0) then
           indvicyclic = indvicyclic + 1
           indviy = iyears
        endif

! Compute and store ndvi database file times in arrays

        fnames_ndvi(jtime) = trim(ndvi_database)//trim(filemonth)
        call date_make_big (indviy,indvim,indvid,indvih,ctotdate_ndvi(jtime))
        call date_abs_secs2(indviy,indvim,indvid,indvih,s1900_ndvi(jtime))

     enddo

     close(29)

! If ndvi cyclic flag > 0, check its value against ntimes and stop if they
! are unequal.  If they are equal, reset ndvi cyclic flag to 1 and augment
! ndvi file arrays by 1 at each end.

     if (indvicyclic > 0) then
        if (indvicyclic /= ntimes) then
           write(io6,'(/,a)') 'Some but not all ndvi database files do not have'
           write(io6,'(a)')   'year 0000, which is ambiguous.  Stopping model.'
           stop 'stop_ndvi_inv'
        endif
        indvicyclic = 1
        nndvifiles = ntimes + 2

! Shift ndvi data file names and times by one array element

        do jtime = ntimes,1,-1
           fnames_ndvi  (jtime+1) = fnames_ndvi  (jtime)
           ctotdate_ndvi(jtime+1) = ctotdate_ndvi(jtime)
           s1900_ndvi   (jtime+1) = s1900_ndvi   (jtime)
        enddo

! Add new ndvi member at beginning of time sequence

        fnames_ndvi(1) = fnames_ndvi(ntimes+1)

        call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(ntimes+1))
        call date_make_big(indviy-1,indvim,indvid,indvih,ctotdate_ndvi(1))
        call date_abs_secs2(indviy-1,indvim,indvid,indvih,s1900_ndvi(1))

! Add new ndvi member at end of time sequence

        fnames_ndvi(ntimes+2) = fnames_ndvi(2)

        call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(2))
        call date_make_big(indviy+1,indvim,indvid,indvih,ctotdate_ndvi(ntimes+2))
        call date_abs_secs2(indviy+1,indvim,indvid,indvih,s1900_ndvi(ntimes+2))

     endif

! Loop over number of NDVI_DATABASE file times and search for the one that
! corresponds to current or most recent model time.

     indvifile = 0
     do nf = 1,nndvifiles

        write(io6,*) 'nndvif0 ',nf,s1900_ndvi(nf),' ',s1900_sim

        if (s1900_ndvi(nf) <= s1900_sim) then
           indvifile = nf
        endif
     enddo

     if (indvifile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find previous or current ndvi file for current'
        write(io6,*) 'model time.  Stopping model '
        stop 'stop: no current ndvi file'
     endif

  elseif (iaction == 1) then

! Processing next ndvi file (only called with iaction = 1 if iupdndvi = 1)

     indvifile = indvifile + 1

     if (indvifile > nndvifiles)then
        if(indvicyclic == 0) then
           write(io6,*) ' '
           write(io6,*) 'No future ndvi file is available for nudging '
           write(io6,*) 'Stopping model '
           stop 'stop: no future ndvi file'
        else
           indvifile = 3
           do jtime = 1, nndvifiles
              call date_unmake_big(indviy,indvim,indvid,indvih,ctotdate_ndvi(jtime))
              call date_make_big(indviy+1,indvim,indvid,indvih,ctotdate_ndvi(jtime))
              call date_abs_secs2(indviy+1,indvim,indvid,indvih,s1900_ndvi(jtime))
           enddo
        endif
     endif

     land%veg_ndvip(:) = land%veg_ndvif(:)

  endif

  ! Allocate and fill temporary lat/lon arrays for land points only

  allocate (glatwl(mland), glonwl(mland)); glatwl(:) = 0.; glonwl(:) = 0.

  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     glatwl(iland) = sfcg%glatw(iwsfc)
     glonwl(iland) = sfcg%glonw(iwsfc)
  enddo

! Read and interpolate current ndvi file

  call land_database_read(mland,                  &
                          glatwl,                 &
                          glonwl,                 &
                          ndvi_database,          &
                          fnames_ndvi(indvifile), &
                          'ndvi',                 &
                          datq=datq               )

  deallocate (glatwl, glonwl)

! Loop over all land cells (already defined and filled with leaf_class)

  do iland = 2,mland
     iwsfc = iland + omland

     ! Skip this cell if running in parallel and cell rank is not MYRANK
     ! if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

     land%veg_ndvif(iland) = max(.05,datq(iland))
  enddo

  if (iaction == 0) then
     land%veg_ndvip(:) = land%veg_ndvif(:)
  endif

end subroutine ndvi_database_read
