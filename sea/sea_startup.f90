subroutine sea_startup()

  use sea_coms,  only: isstflg, seatmp,  &
                       iupdsst, iseaiceflg, iupdseaice, seaice

  use mem_sea,   only: msea, sea, alloc_sea2, filltab_sea

  use misc_coms, only: io6, runtype
  use perth,     only: ntides
  use sea_swm,   only: use_tides

  implicit none

  integer :: isea

  ! Subroutine SEA_STARTUP allocates some sea arrays and initializes sst

  ! THIS SUBROUTINE DOES NOT INITIALIZE canopy temperature and moisture
  ! values, which depend on atmospheric conditions.

  !-------------------------------------------------------------------------------
  ! STEP 1: Call alloc_sea and filltab_sea (sea grid arrays already allocated)
  !-------------------------------------------------------------------------------

  call alloc_sea2(ntides, use_tides)
  call filltab_sea()

  !-------------------------------------------------------------------------------
  ! STEP 2a: Fill sst values
  !-------------------------------------------------------------------------------

  if ( (iupdsst /= 1) .and. &
       (runtype == 'HISTORY' .or. runtype == 'HISTREGRID') ) then

  ! Do nothing if we are restarting and keeping SST constant.
  ! It will be read in from the history file

  elseif (isstflg == 0) then

     ! Default initialization of SST

     do isea = 2,msea
        sea%seatp(isea) = seatmp
        sea%seatf(isea) = seatmp
     enddo

  elseif (isstflg == 1) then

     if (runtype == 'INITIAL' .or. &
         runtype == 'HISTORY' .or. &
         runtype == 'HISTREGRID') then

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
         runtype == 'HISTREGRID') then

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

  if ( (iupdseaice /= 1) .and. &
       (runtype == 'HISTORY' .or. runtype == 'HISTREGRID') ) then

     ! Do nothing if we are restarting and keeping SEAICE constant.
     ! It will be read in from the history file.

  elseif (iseaiceflg == 0) then

     ! Default initialization of SEAICE

     do isea = 2,msea
        sea%seaicep(isea) = seaice
        sea%seaicef(isea) = seaice
     enddo

  elseif (iseaiceflg == 1) then

     if (runtype == 'INITIAL' .or. &
         runtype == 'HISTORY' .or. &
         runtype == 'HISTREGRID') then

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
         runtype == 'HISTREGRID') then

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
