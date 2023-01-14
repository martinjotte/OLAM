Module soilgrids_db

Contains

  subroutine soilgrids_read()

  ! The letter "q" represents any point in the grid stagger, as determined by
  ! the routine that calls this subroutine.

  use consts_coms, only: piu180
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec
  use misc_coms,   only: io6
  use max_dims,    only: pathlen
  use mem_sfcg,    only: sfcg
  use mem_land,    only: nzg, slzt, nland, onland, land
  use leaf_coms,   only: soilgrids_database

  implicit none

  integer :: qtable(nland)

  integer :: nio, njo
  integer :: niosh, njosh
  integer :: isoc, iwoc
  integer :: isocpt, isocpo
  integer :: iwocph, iwocpt, iwocpo
  integer :: ind
  integer :: ind1, ind2
  integer :: io1, io2, jo1, jo2
  integer :: ifiles, jfiles
  integer :: ifile, jfile
  integer :: io_full, jo_full
  integer :: iland, iwsfc
  integer :: nperdeg
  integer :: ndims, idims(2)
  integer :: idataset

  real :: offpix
  real :: qlat1, qlon1
  real :: wio1, wio2, wjo1, wjo2
  real :: rio_full, rjo_full

  real :: datc(nzg)

  integer(kind=4), allocatable :: idato (:,:,:)

  integer, allocatable :: numq    (:,:)
  integer, allocatable :: numqind1(:,:)
  integer, allocatable :: numqind2(:,:) ! (ifiles,jfiles)

  character(len=3)   :: title1
  character(len=4)   :: title2
  character(pathlen) :: fname, data_name, slev

  logical :: l1,l2

  integer :: ksg, k, islev, islev2
  integer :: ksg1(nzg),ksg2(nzg)
  real    :: wsg1(nzg),wsg2(nzg)

  ! Standard SoilGrids layer heights (negative of depth below surface)

  real, parameter :: sgz(7) = (/0.00,-0.05,-0.15,-0.30,-0.60,-1.00,-2.00/)

  ! Fill vertical interpolation indices and weights

  do k = 1,nzg
     if (slzt(k) < sgz(7)) then
        ksg2(k) = 7 ; wsg2(k) = 1.0
        ksg1(k) = 7 ; wsg1(k) = 0.0
     else
        ksg = 7
        do while(sgz(ksg) < slzt(k))
           ksg = ksg - 1
        enddo
        ksg2(k) = ksg;   wsg2(k) = (slzt(k) - sgz(ksg+1)) / (sgz(ksg) - sgz(ksg+1))
        ksg1(k) = ksg+1; wsg1(k) = 1.0 - wsg2(k)
     endif

     write(6,'(a,3i5,3f10.4)') 'ksg,wsg ',k,ksg1(k),ksg2(k),wsg1(k),wsg2(k),slzt(k)

  enddo

  ! Open, read, and close dataset header file

  fname = trim(soilgrids_database)//'HEADER'
  inquire(file=fname, exist=l1)

  if (.not. l1) then
     write(io6,*)
     write(io6,*) '==================================================='
     write(io6,*) '| Problem in soilgrids_read:'
     write(io6,*) '| Header file ', trim(fname)
     write(io6,*) '| not found!'
     write(io6,*) '==================================================='
     stop 'soilgrids_read'
  endif

  open(29, file=fname, form='FORMATTED', status='OLD', action='READ')
  read(29,*) nio, njo, nperdeg
  close(29)

  ! Compute number of pixels in a shift to adjacent file (niosh and njosh).
  ! Compute number of files in database that span all latitudes and
  ! longitudes on earth. [This algorithm will change when multiple resolutions
  ! of the SRTM data become available.]

  if (mod(nio,nperdeg) == 2) then
     offpix = .5
     niosh = nio - 2
     njosh = njo - 2
  else
     offpix = 0.
     niosh = nio - 1
     njosh = njo - 1
  endif

  ifiles = 360 * nperdeg / niosh
  jfiles = 180 * nperdeg / njosh

  ! Allocate 4 arrays.

  allocate (numq    (ifiles,jfiles))
  allocate (numqind1(ifiles,jfiles))
  allocate (numqind2(ifiles,jfiles))
  allocate (idato   (nio,njo,7))

  do jfile = 1,jfiles
     do ifile = 1,ifiles
        numq(ifile,jfile) = 0
     enddo
  enddo

  ! Loop over all geographic points that need data filled, determine which
  ! file to read for each one, and count number of points to be read from
  ! each file.

  do iland = 2,nland
     iwsfc = iland + onland

     qlat1 = max(-89.9999,min(89.9999,sfcg%glatw(iwsfc)))
     qlon1 = sfcg%glonw(iwsfc)

     if (qlon1 >=  180.) qlon1 = qlon1 - 360.
     if (qlon1 <  -180.) qlon1 = qlon1 + 360.
     qlon1 = max(-179.9999,min(179.9999,qlon1))

     rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
     rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

     io_full = int(rio_full)
     jo_full = int(rjo_full)

     ! If io_full and/or jo_full are at max value for database, decrease by 1 so
     ! that out-of-range file will not be sought.

     if (io_full == 360 * nperdeg) io_full = io_full - 1
     if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

     ifile = io_full / niosh + 1
     jfile = jo_full / njosh + 1

     numq(ifile,jfile) = numq(ifile,jfile) + 1

  enddo

  ! Set up array index values for qtable array

  ind = 1
  do jfile = 1,jfiles
     do ifile = 1,ifiles
        numqind1(ifile,jfile) = ind
        numqind2(ifile,jfile) = ind
        ind = ind + numq(ifile,jfile)
     enddo
  enddo

  ! Fill qtable array

  do iland = 2,nland
     iwsfc = iland + onland

     qlat1 = max(-89.9999,min(89.9999,sfcg%glatw(iwsfc)))
     qlon1 = sfcg%glonw(iwsfc)

     if (qlon1 >=  180.) qlon1 = qlon1 - 360.
     if (qlon1 <  -180.) qlon1 = qlon1 + 360.
     qlon1 = max(-179.9999,min(179.9999,qlon1))

     rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
     rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

     io_full = int(rio_full)
     jo_full = int(rjo_full)

     ! If io_full and/or jo_full are at max value for database, decrease by 1 so
     ! that out-of-range file will not be sought.

     if (io_full == 360 * nperdeg) io_full = io_full - 1
     if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

     ifile = io_full / niosh + 1
     jfile = jo_full / njosh + 1

     ind = numqind2(ifile,jfile)

     qtable(ind) = iland
     numqind2(ifile,jfile) = numqind2(ifile,jfile) + 1 ! Sums to ending index

  enddo

  ! Read files and extract data

  do jfile = 1,jfiles
     do ifile = 1,ifiles

        ind1 = numqind1(ifile,jfile)
        ind2 = numqind2(ifile,jfile)

        if (ind2 > ind1) then
           iwoc = (ifile - 1) * niosh / nperdeg - 180  ! SW longitude of current file
           isoc = (jfile - 1) * njosh / nperdeg -  90  ! SW latitude of current file

           ! Construct filename

           isocpt = abs(isoc) / 10
           isocpo = abs(isoc) - isocpt*10
           iwocph = abs(iwoc) / 100
           iwocpt = (abs(iwoc) - iwocph * 100) / 10
           iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10

           if (isoc >= 0) then
              write(title1,'(2i1,a1)') isocpt,isocpo,'N'
           else
              write(title1,'(2i1,a1)') isocpt,isocpo,'S'
           endif

           if (iwoc >= 0) then
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
           else
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
           endif

           ! Loop over Soilgrids datasets

           ! GLHYMPS porosity and permeability are read separately because
           ! they are represented on coarser (1/20 degree) pixels.  SRSTHICK
           ! dataset reserved for possible future use, and would be read here.

           do idataset = 1,9

              ! Set data_name

              if     (idataset == 1) then
                 data_name  = 'GPP'
              elseif (idataset == 2) then
                 data_name  = 'BDTICM'
              elseif (idataset == 3) then
                 data_name  = 'BLDFIE'
              elseif (idataset == 4) then
                 data_name  = 'CECSOL'
              elseif (idataset == 5) then
                 data_name  = 'CLYPPT'
              elseif (idataset == 6) then
                 data_name  = 'ORCDRC'
              elseif (idataset == 7) then
                 data_name  = 'PHIHOX'
              elseif (idataset == 8) then
                 data_name  = 'SLTPPT'
              elseif (idataset == 9) then
                 data_name  = 'SNDPPT'
              endif

              ! Do only one soilgrids level for GPP and BDTICM, and 7 levels
              ! for other SoilGrids datasets

              if (idataset <= 2) then
                 islev2 = 1
              else
                 islev2 = 7
              endif

              ! Loop over all soilgrids levels to process

              do islev = 1,islev2

                 ! Set soil level character string

                 if     (idataset == 1) then
                    slev = '2014_'
                 elseif (idataset == 2) then
                    slev = '_M_1km_'
                 elseif (islev == 1) then
                    slev = '_M_sl1_1km_'
                 elseif (islev == 2) then
                    slev = '_M_sl2_1km_'
                 elseif (islev == 3) then
                    slev = '_M_sl3_1km_'
                 elseif (islev == 4) then
                    slev = '_M_sl4_1km_'
                 elseif (islev == 5) then
                    slev = '_M_sl5_1km_'
                 elseif (islev == 6) then
                    slev = '_M_sl6_1km_'
                 elseif (islev == 7) then
                    slev = '_M_sl7_1km_'
                 endif

                 ! Set file name

                 fname = trim(soilgrids_database)//trim(data_name)//trim(slev)//title1//title2//'.h5'

                 inquire(file=fname,exist=l1,opened=l2)

                 ! Read file

                 if (l1) then
                    write(io6,*) 'getting file ',trim(fname)

                    call shdf5_open(trim(fname),'R')

                    ndims = 2
                    idims(1) = nio
                    idims(2) = njo

                    call shdf5_irec(ndims,idims,trim(data_name),ivar2=idato(:,:,islev))

                    call shdf5_close()
                 else
                    write(io6,*) 'In soilgrids_read, input file is missing'
                    write(io6,*) 'Filename = ',trim(fname)
                    write(io6,*) 'Stopping model run'
                    stop 'stop_soilgrids_input2'
                 endif

              enddo ! islev

              do ind = ind1,ind2-1

                 iland = qtable(ind)
                 iwsfc = iland + onland

                 qlat1 = max(-89.9999,min(89.9999,sfcg%glatw(iwsfc)))
                 qlon1 = sfcg%glonw(iwsfc)

                 if (qlon1 >=  180.) qlon1 = qlon1 - 360.
                 if (qlon1 <= -180.) qlon1 = qlon1 + 360.

                 rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
                 rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

                 io_full = int(rio_full)
                 jo_full = int(rjo_full)

                 ! If io_full and/or jo_full are at max value for database,
                 ! decrease by 1 so that out-of-range file will not be sought.

                 if (io_full == 360 * nperdeg) io_full = io_full - 1
                 if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

                 io1 = mod(io_full,niosh) + 1
                 jo1 = mod(jo_full,njosh) + 1

                 wio2 = rio_full - float(io_full) + offpix
                 wjo2 = rjo_full - float(jo_full) + offpix

                 ! At this point, io1, jo1, wio2, and wjo2 are correct if offpix = 0,
                 ! but need correction if offpix = .5

                 if (wio2 > 1.) then
                    wio2 = wio2 - 1.
                    io1 = io1 + 1
                 endif

                 if (wjo2 > 1.) then
                    wjo2 = wjo2 - 1.
                    jo1 = jo1 + 1
                 endif

                 ! Correction for offpix = .5 is now completed

                 io2 = io1 + 1
                 jo2 = jo1 + 1

                 wio1 = 1. - wio2
                 wjo1 = 1. - wjo2

                 ! Interpolate from 4 surrounding values unless any are "missing" (negative)

                 if (min(idato(io1,jo1,1), idato(io1,jo2,1), &
                         idato(io2,jo1,1), idato(io2,jo2,1)) >= 0) then

                    do islev = 1,islev2
                       datc(islev) &
                          = wio1 * (wjo1 * idato(io1,jo1,islev) + wjo2 * idato(io1,jo2,islev)) &
                          + wio2 * (wjo1 * idato(io2,jo1,islev) + wjo2 * idato(io2,jo2,islev))
                    enddo

                 ! Otherwise, fill from first neighbor found that is not missing

                 elseif (idato(io1,jo1,1) >= 0) then

                    do islev = 1,islev2
                       datc(islev) = idato(io1,jo1,islev)
                    enddo

                 elseif (idato(io2,jo2,1) >= 0) then

                    do islev = 1,islev2
                       datc(islev) = idato(io2,jo2,islev)
                    enddo

                 elseif (idato(io1,jo2,1) >= 0) then

                    do islev = 1,islev2
                       datc(islev) = idato(io1,jo2,islev)
                    enddo

                 elseif (idato(io2,jo1,1) >= 0) then

                    do islev = 1,islev2
                       datc(islev) = idato(io2,jo1,islev)
                    enddo

                 ! If all neighbor points have missing values, cycle past this
                 ! soil model column.  Column will thus retain values that were
                 ! already filled in subroutine usda_composition.

                 else

                    cycle

                 endif

                 ! Apply soilgrids data to OLAM-SOIL grid.  For soilgrids datasets
                 ! that are multi-leveled (covering the uppermost 2 m), interpolate
                 ! vertically to the OLAM-SOIL grid within that vertical range and
                 ! assign deepest ('sl7') dataset value to deeper OLAM-SOIL grid
                 ! levels.  (Some of these deeper values will be ignored, depending
                 ! on bedrock depth and model input parameters.)

                 if (data_name == 'GPP') then
                    land%gpp(iland) = datc(1)
                 elseif (data_name == 'BDTICM') then
                    land%z_bedrock(iland) = datc(1) * (-0.01) ! Converting from [cm] to [-m]
                 elseif (data_name == 'BLDFIE') then
                    do k = 1,nzg
                       land%bulkdens_drysoil(k,iland) = wsg2(k) * datc(ksg2(k)) &
                                                      + wsg1(k) * datc(ksg1(k)) ! [kg/m^3]
                    enddo

                 elseif (data_name == 'CECSOL') then
                    do k = 1,nzg
                       land%cec_soil(k,iland) = wsg2(k) * datc(ksg2(k)) &
                                              + wsg1(k) * datc(ksg1(k))         ! [cmol+/kg]
                    enddo
                 elseif (data_name == 'CLYPPT') then
                    do k = 1,nzg
                       land%clay(k,iland) = (wsg2(k) * datc(ksg2(k)) &          ! Converting from
                                          +  wsg1(k) * datc(ksg1(k))) * 0.01    ! [%] to [fraction]
                    enddo
                 elseif (data_name == 'ORCDRC') then
                    do k = 1,nzg
                       land%organ(k,iland) = (wsg2(k) * datc(ksg2(k)) &         !  Converting from
                                           +  wsg1(k) * datc(ksg1(k))) * 0.001  ! [g/kg] to [kg/kg]
                    enddo
                 elseif (data_name == 'PHIHOX') then
                    do k = 1,nzg
                       land%pH_soil(k,iland) = (wsg2(k) * datc(ksg2(k)) &       ! Converting from
                                             +  wsg1(k) * datc(ksg1(k))) * 0.1  ! [10 * pH] to [pH]
                    enddo
                 elseif (data_name == 'SLTPPT') then
                    do k = 1,nzg
                       land%silt(k,iland) = (wsg2(k) * datc(ksg2(k)) &          ! Converting from
                                          +  wsg1(k) * datc(ksg1(k))) * 0.01    ! [%] to [fraction]
                    enddo
                 elseif (data_name == 'SNDPPT') then
                    do k = 1,nzg
                       land%sand(k,iland) = (wsg2(k) * datc(ksg2(k)) &          ! Converting from
                                          +  wsg1(k) * datc(ksg1(k))) * 0.01    ! [%] to [fraction]
                    enddo
                 endif

              enddo    ! ind

           enddo    ! idataset

        endif    ! ind2 > ind1
     enddo    ! ifile
  enddo    ! jfile

  print*, 'soilgrids_read final deallocation '

  deallocate(numq,numqind1,numqind2,idato)

  end subroutine soilgrids_read

End module soilgrids_db
