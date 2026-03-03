subroutine pom_startup()

  use pom2k1d,     only: nzpom, dtl2, cbc, cbcmin, cbcmax, y, yy, &
                         vonk, z0b, pom, alloc_pom, filltab_pom

  use consts_coms, only: omega2, pio180

  use mem_sfcg,    only: sfcg

  use mem_sea,     only: sea, msea, omsea, npomzons

  use leaf_coms,   only: dt_leaf

  implicit none

  integer :: isea, iwsfc, k

  if (npomzons < 1) return

  ! Later, get these from OLAM simulation parameters

  dtl2 = dt_leaf * 2.0

  call alloc_pom(msea)
  call filltab_pom()

  do isea = 2,msea
     iwsfc = isea + omsea

     if (sea%pom_active(isea)) then
        pom%cor(isea) = omega2 * sin(sfcg%glatw(iwsfc) * pio180)
     endif
  enddo

  do k = 1,nzpom-1
     cbc(k) = (vonk / log((yy(k) - y(k+1)) / z0b))**2
     cbc(k) = max(cbcmin, cbc(k))

     ! If the following is invoked, then it is probable that the wrong
     ! choice of z0b or vertical spacing has been made:

     cbc(k) = min(cbcmax, cbc(k))
  enddo
  cbc(nzpom) = cbc(nzpom-1)

end subroutine pom_startup

!===============================================================================

subroutine pom_init()

  use mem_sea,  only: sea, msea, omsea, npomzons
  use pom2k1d,  only: nzpom, pom, y, yy, small
  use sea_coms, only: pom_idata

  implicit none

  integer :: k, isea, iwsfc
  real    :: tlen  ! turbulent length scale [m]
  real    :: seatempc, deeptempc

  if (npomzons < 1) return

  do isea = 2,msea
     iwsfc = isea + omsea

     ! Fill surface temperature of all POM cells, even if outside pom_active regions

     seatempc = sea%seatc(isea) - 273.15
        deeptempc = min(4.,seatempc)

     pom%potmp(1,isea) = seatempc

     if (sea%pom_active(isea)) then

        tlen = 0.1 * y( sea%pom_kba(isea) )

        pom%wubot(isea) = 0.
        pom%wvbot(isea) = 0.

        ! Default initialization: Use seatempc at surface, decrease linearly with depth
        ! to min(4,seatempc) at 1000 m depth, and use constant min(4,seatempc) below 1000 m.

        deeptempc = min(4.,seatempc)

        do k = 1,nzpom
           pom%potmp(k,isea) = seatempc + (deeptempc - seatempc) &
                             * min(1.0, (yy(1) - yy(k)) / 1000.)
           pom%salin(k,isea) = 35.
           pom%q2   (k,isea) = small
           pom%q2l  (k,isea) = pom%q2(k,isea) * tlen
           pom%u    (k,isea) = 0.
           pom%v    (k,isea) = 0.

           pom%potmpb(k,isea) = pom%potmp(k,isea)
           pom%salinb(k,isea) = pom%salin(k,isea)
           pom%q2b   (k,isea) = pom%q2   (k,isea)
           pom%q2lb  (k,isea) = pom%q2l  (k,isea)
           pom%ub    (k,isea) = pom%u    (k,isea)
           pom%vb    (k,isea) = pom%v    (k,isea)

           pom%kh    (k,isea) = tlen * sqrt(pom%q2b(k,isea))
           pom%km    (k,isea) = pom%kh(k,isea)
           pom%kq    (k,isea) = pom%kh(k,isea)
        enddo

     endif
  enddo

  ! If pom_idata > 0, override default initialization by reading in 3D temperature and salinity data 

  if (pom_idata == 1) then
     call rtofs_read()
  elseif (pom_idata == 2) then
     call hycom_read()
  elseif (pom_idata == 3) then
     call glorys_read()
  endif

end subroutine pom_init

!===============================================================================

subroutine rtofs_read()

  use pom2k1d,     only: nzpom, yy, pom
  use mem_sea,     only: sea, nsea, onsea
  use mem_sfcg,    only: sfcg
  use consts_coms, only: piu180, pio180, pio2
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use sea_coms,    only: pom_database
  use misc_coms,   only: io6

  implicit none

  integer, parameter :: nio = 742, njo = 1710, nko = 40

  real, allocatable :: dummy(:,:,:,:), dato(:,:,:,:), datocol(:,:)

  integer :: io, jo, io1, jo1, io2, jo2, k, k1, k2, kofs, idatotyp
  integer :: isea, iwsfc
  integer :: ndims, idims(4)

  real :: glat, glon, rio, rjo, wx1, wx2, wy1, wy2, w, wsum

  logical :: exists

  integer :: kofs1(nzpom),kofs2(nzpom)
  real    :: wofs1(nzpom),wofs2(nzpom)

  ! rtofs layer heights (negative of depth below surface)

  real, parameter :: zofs(40) = (/   0.,   -2.,   -4.,   -6.,   -8.,  -10.,  -12.,  -15., &
                                   -20.,  -25.,  -30.,  -35.,  -40.,  -45.,  -50.,  -60., &
                                   -70.,  -80.,  -90., -100., -125., -150., -200., -250., &
                                  -300., -350., -400., -500., -600., -700., -800., -900., &
                                 -1000.,-1250.,-1500.,-2000.,-2500.,-3000.,-4000.,-5000. /)

  ! This subroutine reads ocean temperature and salinity from the following rtofs file
  ! (or one of the same type at a different date/time) to be used by POM1D in OLAM:
  ! nomads.ncep.noaa.gov/pub/data/nccf/com/rtofs/prod/rtofs20211127/
  ! rtofs_glo_3dz_n006_6hrly_hvr_US_east.nc
  ! This class of rtofs files (containing "US_east" in their name) covers only the North
  ! Atlantic from Longitude 100 W to 40.72 W and from the equator northward.  Temperature
  ! and salinity are 4D arrays with dimensions (742,1710,40,1), which are hardwired herein.
  ! The rtofs grid projection is Mercator south of 47 N and bipolar north of 47 N.  For
  ! simplicity, this subroutine is hardwired to interpolate only from points south of 47 N.
  ! This subroutine will require modification if it is ever to be used to input rtofs data
  ! outside the lat/lon limits stated above.

  ! Check for existence of the rtofs file

  call shdf5_exists(trim(pom_database), exists)

  if (.not. exists) then
     write(io6,*) 'rtofs: error opening rtofs file ' // trim(pom_database)
     stop 'no rtofs file'
  endif

  ! Allocate rtofs input array and column array

  allocate (dummy(nio,njo,nko,1), dato(nio,njo,nko,2), datocol(nko,2))

  ! Open rtofs file and read 2 fields from it

  call shdf5_open(trim(pom_database), 'R')

  ndims = 4
  idims(1) = nio
  idims(2) = njo
  idims(3) = nko
  idims(4) = 1

  call shdf5_irec(ndims,idims,'temperature',rvar4=dummy)
  dato(:,:,:,1) = dummy(:,:,:,1)
  call shdf5_irec(ndims,idims,'salinity'   ,rvar4=dummy)
  dato(:,:,:,2) = dummy(:,:,:,1)

  call shdf5_close()

  ! Fill vertical interpolation indices and weights to interpolate vertically from
  ! rtofs levels to POM1D levels

  do k = 1,nzpom
     if (yy(k) < zofs(40)) then
        kofs2(k) = 40 ; wofs2(k) = 1.0
        kofs1(k) = 40 ; wofs1(k) = 0.0
     else
        kofs = 40
        do while(zofs(kofs) < yy(k))
           kofs = kofs - 1
        enddo
        kofs2(k) = kofs;   wofs2(k) = (yy(k) - zofs(kofs+1)) / (zofs(kofs) - zofs(kofs+1))
        kofs1(k) = kofs+1; wofs1(k) = 1.0 - wofs2(k)
     endif

     write(6,'(a,3i5,3f10.4)') 'kofs,wofs ',k,kofs1(k),kofs2(k),wofs1(k),wofs2(k),yy(k)
  enddo

  ! Loop over sea cells and select those where POM1D is active

  do isea = 2, nsea

     if (.not. sea%pom_active(isea)) cycle

     iwsfc = isea + onsea

     glat = sfcg%glatw(iwsfc)
     glon = sfcg%glonw(iwsfc)

     ! Exclude sea points that are outside rtofs (US_east) Mercator zone

     if (glat < 0.001 .or. glat > 46.999 .or. glon < -99.999 .or. glon > -40.721) cycle

     ! Westernmost rtofs data values are at -100.00 deg longitude and are spaced
     ! 12.5 values per degree of longitude

     rio = (glon + 100.) * 12.5 + 1.0

     ! rtofs latitudes are equally spaced in Mercator projection (up to 47 N latitude),
     ! with 12.5 values per degree of latitude at equator

     rjo = 12.5 * piu180 * log( tan(0.5 * (glat * pio180 + pio2)) ) + 1.0

     io1 = int(rio)
     jo1 = int(rjo)

     io2 = io1 + 1
     jo2 = jo1 + 1

     wx2 = rio - real(io1)
     wy2 = rjo - real(jo1)

     wx1 = 1. - wx2
     wy1 = 1. - wy2

     ! Loop over all vertical levls in RTOFS and interpolate temperature and salinity
     ! values horizontally to location of isea column.  Avoid interpolating missing values,
     ! but if all 4 surrounding values at a given level are missing, set result to missing.

     do k = 1,40
        do idatotyp = 1,2 ! 1 for temperature, 2 for salinity

           wsum = 0.0

           datocol(k,idatotyp) = 0.

           if (abs(dato(io1,jo1,k,idatotyp)) < 50.) then
              w    = wx1 * wy1
              wsum = w
              datocol(k,idatotyp) = w * dato(io1,jo1,k,idatotyp)
           endif

           if (abs(dato(io2,jo1,k,idatotyp)) < 50.) then
              w    = wx2 * wy1
              wsum = wsum + w
              datocol(k,idatotyp) = datocol(k,idatotyp) + w * dato(io2,jo1,k,idatotyp)
           endif

           if (abs(dato(io1,jo2,k,idatotyp)) < 50.) then
              w    = wx2 * wy1
              wsum = wsum + w
              datocol(k,idatotyp) = datocol(k,idatotyp) + w * dato(io1,jo2,k,idatotyp)
           endif

           if (abs(dato(io2,jo2,k,idatotyp)) < 50.) then
              w    = wx2 * wy1
              wsum = wsum + w
              datocol(k,idatotyp) = datocol(k,idatotyp) + w * dato(io2,jo2,k,idatotyp)
           endif

           if (wsum > 1.e-9) then
              datocol(k,idatotyp) = datocol(k,idatotyp) / wsum
           else
              datocol(k,idatotyp) = 1.e10 ! missing value
           endif

        enddo ! idatotyp
     enddo    ! k

     ! Loop over all active points in current POM1D column and vertically interpolate 
     ! temperature and salinity from horizontally-interpolated rtofs column 

     do k = 1,sea%pom_kba(isea)

        k1 = kofs1(k)
        k2 = kofs2(k)

        ! Interpolate if both upper and lower rtofs values are ok.  Assign from upper
        ! rtofs value if it alone is ok.  Otherwise, do not reassign POM values.

        if (abs(datocol(k2,1)) < 100.) then
           if (abs(datocol(k1,1)) < 100.) then
              pom%potmp(k,isea) = wofs2(k) * datocol(k2,1) &
                                + wofs1(k) * datocol(k1,1)
           else
              pom%potmp(k,isea) = datocol(k2,1)
           endif
           pom%potmpb(k,isea) = pom%potmp(k,isea)
        endif

        if (abs(datocol(k2,2)) < 100.) then
           if (abs(datocol(k1,2)) < 100.) then
              pom%salin(k,isea) = wofs2(k) * datocol(k2,2) &
                                + wofs1(k) * datocol(k1,2)
           else
              pom%salin(k,isea) = datocol(k2,2)
           endif
           pom%salinb(k,isea) = pom%salin(k,isea)
        endif

     enddo ! k

  enddo    ! isea

  deallocate(dato,datocol)

end subroutine rtofs_read

!===============================================================================

subroutine hycom_read()

  use pom2k1d,     only: nzpom, yy, pom
  use mem_sea,     only: sea, nsea, onsea
  use mem_sfcg,    only: sfcg
  use consts_coms, only: piu180, pio180, pio2
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use sea_coms,    only: pom_database
  use misc_coms,   only: io6

  implicit none

  integer, parameter :: nio = 4500, njo = 4251, nko = 40

  real, allocatable :: dummy(:,:,:,:), temp(:,:,:), salin(:,:,:), tempcol(:), salincol(:)

  integer :: io, jo, io1, jo1, io2, jo2, k, k1, k2, kofs, idatotyp
  integer :: isea, iwsfc
  integer :: ndims, idims(4)

  real :: glat, glon, rio, rjo, wx1, wx2, wy1, wy2, w, wsum

  character(pathlen) :: temp_fname, salin_fname

  logical :: exists

  integer :: kofs1(nzpom),kofs2(nzpom)
  real    :: wofs1(nzpom),wofs2(nzpom)

  ! rtofs layer heights (negative of depth below surface)

  real, parameter :: zofs(nko) = (/   0.,   -2.,   -4.,   -6.,   -8.,  -10.,  -12.,  -15., &
                                    -20.,  -25.,  -30.,  -35.,  -40.,  -45.,  -50.,  -60., &
                                    -70.,  -80.,  -90., -100., -125., -150., -200., -250., &
                                   -300., -350., -400., -500., -600., -700., -800., -900., &
                                  -1000.,-1250.,-1500.,-2000.,-2500.,-3000.,-4000.,-5000. /)

  ! This subroutine reads ocean temperature and salinity from the following hycom files
  ! (or others of the same type at a different date/time) to be used by POM1D in OLAM:
  ! ncss/hycom/org/thredds/ncss/grid/ESPC-D-V02/t3z/2024/dataset.html/t3z_2024.nc4
  ! and
  ! ncss/hycom/org/thredds/ncss/grid/ESPC-D-V02/s3z/2024/dataset.html/s3z_2024.nc4
  ! These hycom files use a global (north of 80 deg S) lat/lon grid with uniform 0.04 degree
  ! spacing in latitude and 0.08 degree spacing in longitude.  Temperature and salinity are
  ! 4D arrays with dimensions (4500,4251,40,1), which are hardwired herein.

  ! Check for existence of the hycom temperature and salinity files

  temp_fname = trim(pom_database)//'_t3z.nc4'
  call shdf5_exists(trim(temp_fname), exists)
  if (.not. exists) then
     write(io6,*) 'hycom temperature file not found '//trim(temp_fname)
     stop 'no hycom temp file'
  endif

  salin_fname = trim(pom_database)//'_s3z.nc4'
  call shdf5_exists(trim(salin_fname), exists)
  if (.not. exists) then
     write(io6,*) 'hycom salinity file not found '//trim(salin_fname)
     stop 'no hycom_salinity file'
  endif

  ! Allocate temp, salin, and column arrays

  allocate (dummy(nio,njo,nko,1), temp(nio+1,njo,nko), salin(nio+1,njo,nko), tempcol(nko), salincol(nko))

  ndims = 4
  idims(1) = nio
  idims(2) = njo
  idims(3) = nko
  idims(4) = 1

  ! Open hycom_temp_file, read temperature from it, and copy to temp array with longitude shift
  ! since hycom file longitudes range from 0 to 360.  Also apply scale factor and offset.

  call shdf5_open(trim(temp_fname), 'R')
  call shdf5_irec(ndims,idims,'water_temp',rvar4=dummy)
  call shdf5_close()
  temp(   1:2250,:,:) = dummy(2251:4500,:,:,1) * 0.001 + 20.
  temp(2251:4501,:,:) = dummy(   1:2251,:,:,1) * 0.001 + 20.

  ! Open hycom_salin_file, read salinity from it, and copy to salin array with longitude shift
  ! since hycom file longitudes range from 0 to 360.  Also apply scale factor and offset.

  call shdf5_open(trim(salin_fname), 'R')
  call shdf5_irec(ndims,idims,'salinity',rvar4=dummy)
  call shdf5_close()
  salin(   1:2250,:,:) = dummy(2251:4500,:,:,1) * 0.001 + 20.
  salin(2251:4501,:,:) = dummy(   1:2251,:,:,1) * 0.001 + 20.

  ! Fill vertical interpolation indices and weights to interpolate vertically from
  ! rtofs levels to POM1D levels

  do k = 1,nzpom
     if (yy(k) < zofs(nko)) then
        kofs2(k) = nko ; wofs2(k) = 1.0
        kofs1(k) = nko ; wofs1(k) = 0.0
     else
        kofs = nko
        do while(zofs(kofs) < yy(k))
           kofs = kofs - 1
        enddo
        kofs2(k) = kofs;   wofs2(k) = (yy(k) - zofs(kofs+1)) / (zofs(kofs) - zofs(kofs+1))
        kofs1(k) = kofs+1; wofs1(k) = 1.0 - wofs2(k)
     endif

     write(6,'(a,3i5,3f11.4)') 'kofs,wofs ',k,kofs1(k),kofs2(k),wofs1(k),wofs2(k),yy(k)
  enddo

  ! Loop over sea cells and select those where POM1D is active

  do isea = 2, nsea

     if (.not. sea%pom_active(isea)) cycle

     iwsfc = isea + onsea

     glat = max( -89.999,min( 89.999,sfcg%glatw(iwsfc)))
     glon = max(-179.999,min(179.999,sfcg%glonw(iwsfc)))

     ! Westernmost hycom data values are at -180.00 deg longitude and are spaced
     ! 12.5 values per degree of longitude

     rio = (glon + 180.) * 12.5 + 1.0

     ! Southernmost hycom data values are at -80.00 deg latitude and are spaced
     ! with 25.0 values per degree of latitude

     rjo = (glat + 80.) * 25.0 + 1.0

     io1 = int(rio)
     jo1 = int(rjo)

     io2 = io1 + 1
     jo2 = jo1 + 1

     wx2 = rio - real(io1)
     wy2 = rjo - real(jo1)

     wx1 = 1. - wx2
     wy1 = 1. - wy2

     ! Loop over all vertical levls in RTOFS and interpolate temperature and salinity values
     ! horizontally to location of isea column.  Avoid interpolating missing values, but if
     ! all 4 surrounding values at a given level are missing, set result to missing.  For both
     ! temperature and salinity, missing value is -10 (after applying scale factor and offset).

     do k = 1,nko
        wsum = 0.0
        tempcol(k) = 0.

        if (temp(io1,jo1,k) > -9.) then
           w    = wx1 * wy1
           wsum = w
           tempcol(k) = w * temp(io1,jo1,k)
        endif

        if (temp(io2,jo1,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io2,jo1,k)
        endif

        if (temp(io1,jo2,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io1,jo2,k)
        endif

        if (temp(io2,jo2,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io2,jo2,k)
        endif

        if (wsum > 1.e-9) then
           tempcol(k) = tempcol(k) / wsum
        else
           tempcol(k) = -10. ! missing value
        endif

        wsum = 0.0
        salincol(k) = 0.

        if (salin(io1,jo1,k) >= 0.) then
           w    = wx1 * wy1
           wsum = w
           salincol(k) = w * salin(io1,jo1,k)
        endif

        if (salin(io2,jo1,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io2,jo1,k)
        endif

        if (salin(io1,jo2,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io1,jo2,k)
        endif

        if (salin(io2,jo2,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io2,jo2,k)
        endif

        if (wsum > 1.e-9) then
           salincol(k) = salincol(k) / wsum
        else
           salincol(k) = -10. ! missing value
        endif

     enddo    ! k

     ! Loop over all active points in current POM1D column and vertically interpolate 
     ! temperature and salinity from horizontally-interpolated hycom column 

     do k = 1,sea%pom_kba(isea)

        k1 = kofs1(k)
        k2 = kofs2(k)

        ! Interpolate if both upper and lower rtofs values are not missing.  Assign from upper
        ! rtofs value if it alone is not missing.  Otherwise, do not reassign POM values.

        if (tempcol(k2) > -9.) then
           if (tempcol(k1) > -9.) then
              pom%potmp(k,isea) = wofs2(k) * tempcol(k2) &
                                + wofs1(k) * tempcol(k1)
           else
              pom%potmp(k,isea) = tempcol(k2)
           endif
           pom%potmpb(k,isea) = pom%potmp(k,isea)
        endif

        if (salincol(k2) >= 0.) then
           if (salincol(k1) >= 0.) then
              pom%salin(k,isea) = wofs2(k) * salincol(k2) &
                                + wofs1(k) * salincol(k1)
           else
              pom%salin(k,isea) = salincol(k2)
           endif
           pom%salinb(k,isea) = pom%salin(k,isea)
        endif

     enddo ! k

  enddo    ! isea

  deallocate(temp, salin, tempcol, salincol)

end subroutine hycom_read

!===============================================================================

subroutine glorys_read()

  use pom2k1d,     only: nzpom, yy, pom
  use mem_sea,     only: sea, nsea, onsea
  use mem_sfcg,    only: sfcg
  use consts_coms, only: piu180, pio180, pio2
  use hdf5_utils,  only: shdf5_exists, shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,    only: pathlen
  use sea_coms,    only: pom_database
  use misc_coms,   only: io6

  implicit none

  integer, parameter :: nio = 4320, njo = 2041, nko = 50

  integer, allocatable :: dummy(:,:,:,:)
  real,    allocatable :: temp(:,:,:), salin(:,:,:), tempcol(:), salincol(:)

  integer :: io, jo, io1, jo1, io2, jo2, k, k1, k2, kofs, idatotyp
  integer :: isea, iwsfc
  integer :: ndims, idims(4)

  real :: glat, glon, rio, rjo, wx1, wx2, wy1, wy2, w, wsum

  logical :: exists

  integer :: kofs1(nzpom),kofs2(nzpom)
  real    :: wofs1(nzpom),wofs2(nzpom)

  ! rtofs layer heights (negative of depth below surface)

  real, parameter :: zofs(nko) = (/    -0.49,    -1.54,    -2.65,    -3.82,    -5.08, &
                                       -6.44,    -7.93,    -9.57,   -11.40,   -13.47, &
                                      -15.81,   -18.50,   -21.60,   -25.21,   -29.44, &
                                      -34.43,   -40.34,   -47.37,   -55.76,   -65.81, &
                                      -77.85,   -92.33,  -109.73,  -130.67,  -155.85, &
                                     -186.13,  -222.48,  -266.04,  -318.13,  -380.21, &
                                     -453.94,  -541.09,  -643.57,  -763.33,  -902.34, &
                                    -1062.44, -1245.29, -1452.25, -1684.28, -1941.89, &
                                    -2225.08, -2533.34, -2865.7,  -3220.82, -3597.03, &
                                    -3992.48, -4405.22, -4833.29, -5274.78, -5727.92 /)

  ! This subroutine reads ocean temperature and salinity from the following Copernicus file
  ! (or others of the same type at a different date/time) to be used by POM1D in OLAM:
  ! mercatorglorys12v1_gl12_mean_20241007_R20241009.nc
  ! Downloaded from: https://data.marine.copernicus.eu
  ! These Copernicus files use a global lat/lon grid (from -80 lat to 90 lat) with uniform (1/12) degree
  ! spacing in latitude and longitude.  Temperature and salinity are
  ! 4D arrays with dimensions (4320,2041,50,1), which are hardwired herein.

  ! Check for existence of the Copernicus file

  call shdf5_exists(trim(pom_database), exists)
  if (.not. exists) then
     write(io6,*) 'Copernicus ocean file not found '//trim(pom_database)
     stop 'no Copernicus ocean file'
  endif

  ! Allocate temp, salin, and column arrays

  allocate (dummy(nio,njo,nko,1), temp(nio+1,njo,nko), salin(nio+1,njo,nko), tempcol(nko), salincol(nko))

  ndims = 4
  idims(1) = nio
  idims(2) = njo
  idims(3) = nko
  idims(4) = 1

  ! Open Copernicus ocean file, read temperature and salinity from it, and copy to
  ! temp and salin arrays.  Also apply scale factor and offset.

  call shdf5_open(trim(pom_database), 'R')

  call shdf5_irec(ndims,idims,'thetao',ivar4=dummy)
  temp(1:4320,:,:) = real(dummy(1:4320,:,:,1)) * .7324442e-3 + 21.  ! 0.000732444226741791 + 21.
  temp  (4321,:,:) = temp(1,:,:)

  call shdf5_irec(ndims,idims,'so',ivar4=dummy)
  salin(1:4320,:,:) = real(dummy(1:4320,:,:,1)) * .1525926e-2 - .1525926e-2  ! * 0.00152592547237873 - 0.00152592547237873
  salin  (4321,:,:) = salin(1,:,:)

  call shdf5_close()

  ! Fill vertical interpolation indices and weights to interpolate vertically from
  ! rtofs levels to POM1D levels

  do k = 1,nzpom
     if (yy(k) < zofs(nko)) then
        kofs2(k) = nko ; wofs2(k) = 1.0
        kofs1(k) = nko ; wofs1(k) = 0.0
     elseif (yy(k) >= zofs(1)) then
        kofs2(k) = 1 ; wofs2(k) = 1.0
        kofs1(k) = 1 ; wofs1(k) = 0.0
     else
        kofs = nko
        do while(zofs(kofs) < yy(k))
           kofs = kofs - 1
        enddo
        kofs2(k) = kofs;   wofs2(k) = (yy(k) - zofs(kofs+1)) / (zofs(kofs) - zofs(kofs+1))
        kofs1(k) = kofs+1; wofs1(k) = 1.0 - wofs2(k)
     endif

     write(6,'(a,3i5,3f11.4)') 'kofs,wofs ',k,kofs1(k),kofs2(k),wofs1(k),wofs2(k),yy(k)
  enddo

  ! Loop over sea cells and select those where POM1D is active

  do isea = 2, nsea

     if (.not. sea%pom_active(isea)) cycle

     iwsfc = isea + onsea

     glat = max( -89.999,min( 89.999,sfcg%glatw(iwsfc)))
     glon = max(-179.999,min(179.999,sfcg%glonw(iwsfc)))

     ! Westernmost Copernicus data values are at -180.00 deg longitude and are spaced
     ! 12.0 values per degree of longitude

     rio = (glon + 180.) * 12.0 + 1.0

     ! Southernmost Copernicus data values are at -80.00 deg latitude and are spaced
     ! with 12.0 values per degree of latitude

     rjo = (glat + 80.) * 12.0 + 1.0

     io1 = int(rio)
     jo1 = int(rjo)

     io2 = io1 + 1
     jo2 = jo1 + 1

     wx2 = rio - real(io1)
     wy2 = rjo - real(jo1)

     wx1 = 1. - wx2
     wy1 = 1. - wy2

     ! Loop over all vertical levls in Copernicus and interpolate temperature and salinity values
     ! horizontally to location of isea column.  Avoid interpolating missing values, but if
     ! all 4 surrounding values at a given level are missing, set result to missing.  For both
     ! temperature and salinity, missing value is -10 (after applying scale factor and offset).

     do k = 1,nko
        wsum = 0.0
        tempcol(k) = 0.

        if (temp(io1,jo1,k) > -9.) then
           w    = wx1 * wy1
           wsum = w
           tempcol(k) = w * temp(io1,jo1,k)
        endif

        if (temp(io2,jo1,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io2,jo1,k)
        endif

        if (temp(io1,jo2,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io1,jo2,k)
        endif

        if (temp(io2,jo2,k) > -9.) then
           w    = wx2 * wy1
           wsum = wsum + w
           tempcol(k) = tempcol(k) + w * temp(io2,jo2,k)
        endif

        if (wsum > 1.e-9) then
           tempcol(k) = tempcol(k) / wsum
        else
           tempcol(k) = -10. ! missing value
        endif

        wsum = 0.0
        salincol(k) = 0.

        if (salin(io1,jo1,k) >= 0.) then
           w    = wx1 * wy1
           wsum = w
           salincol(k) = w * salin(io1,jo1,k)
        endif

        if (salin(io2,jo1,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io2,jo1,k)
        endif

        if (salin(io1,jo2,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io1,jo2,k)
        endif

        if (salin(io2,jo2,k) >= 0.) then
           w    = wx2 * wy1
           wsum = wsum + w
           salincol(k) = salincol(k) + w * salin(io2,jo2,k)
        endif

        if (wsum > 1.e-9) then
           salincol(k) = salincol(k) / wsum
        else
           salincol(k) = -10. ! missing value
        endif

     enddo    ! k

     ! Loop over all active points in current POM1D column and vertically interpolate 
     ! temperature and salinity from horizontally-interpolated hycom column 

     do k = 1,sea%pom_kba(isea)

        k1 = kofs1(k)
        k2 = kofs2(k)

        ! Interpolate if both upper and lower rtofs values are not missing.  Assign from upper
        ! rtofs value if it alone is not missing.  Otherwise, do not reassign POM values.

        if (tempcol(k2) > -9.) then
           if (tempcol(k1) > -9.) then
              pom%potmp(k,isea) = wofs2(k) * tempcol(k2) &
                                + wofs1(k) * tempcol(k1)
           else
              pom%potmp(k,isea) = tempcol(k2)
           endif
           pom%potmpb(k,isea) = pom%potmp(k,isea)
        endif

        if (salincol(k2) >= 0.) then
           if (salincol(k1) >= 0.) then
              pom%salin(k,isea) = wofs2(k) * salincol(k2) &
                                + wofs1(k) * salincol(k1)
           else
              pom%salin(k,isea) = salincol(k2)
           endif
           pom%salinb(k,isea) = pom%salin(k,isea)
        endif

     enddo ! k

  enddo    ! isea

  deallocate(temp, salin, tempcol, salincol)

end subroutine glorys_read

!===============================================================================

subroutine plot_pom()

  use mem_sea, only: omsea, sea
  use pom2k1d, only: nzpom, pom, yy

  implicit none

  real, parameter :: aspect = .7
  real, parameter :: scalelab = .014

  real, allocatable :: vctr18(:), val(:,:)

  integer :: labincx, labincy, iwsfc, isea, iv, k, jpom, kb

  integer, parameter :: icolor(6) = [16, 109, 12, 11, 9, 8]
! orange, blue green, dark red, purple, dark green, dark blue

  real :: xmin, xmax, xinc
  real :: ymin, ymax, yinc

  integer, parameter :: ijpom(4) = [179785, 145400, 139184, 157375] ! selected POM1D columns to plot (iwsfc values)
! integer, parameter :: ijpom(4) = [120592,  86207,  79991,  98182] ! selected POM1D columns to plot
  character(len=1), parameter :: ip(4) = ['1','2','3','4']

  labincx = 5
  labincy = 5

  xmin = -11.
  xmax = 41.
  xinc = 2.

  ymin =  1.05 * yy(nzpom)
  ymax = -0.05 * yy(nzpom)
  yinc = 100.

  call o_reopnwk()

  call plotback()

  do jpom = 1,4
     iwsfc = ijpom(jpom)
     isea = iwsfc - omsea

     kb = sea%pom_kba(isea)

     allocate(vctr18(kb-1))
     allocate(val(kb-1,6))

     do k = 1,kb-1
        vctr18(k) = yy(k)

        val(k,1) = pom%potmp(k,isea)
        val(k,2) = pom%salin(k,isea)
        val(k,3) = pom%q2(k,isea) * 100.
        val(k,4) = pom%q2l(k,isea) * 100.
        val(k,5) = pom%u(k,isea) * 10.
        val(k,6) = pom%v(k,isea) * 10.
     enddo

     do iv = 1,6

        call oplot_xy2l(ip(jpom),'N','a','N',aspect,scalelab,icolor(iv),0, &
                       kb-1,  val(:,iv), &
                       vctr18, &
                       'vals','Z', &
                       xmin, xmax, xinc, labincx, &
                       ymin, ymax, yinc, labincy  )

     enddo

     deallocate(vctr18,val)

  enddo

  call o_frame()

  call o_clswk()

end subroutine plot_pom

!===============================================================================

subroutine oplot_xy2l(panel,frameoff,pltborder,colorbar0,aspect,scalelab,&
                      linecolor,ndashes,n,xval,yval,xlab,ylab, &
                      xmin,xmax,xinc,labincx,ymin,ymax,yinc,labincy)

  use oplot_coms, only: op

! This routine is a substitute for NCAR Graphics routine ezxy to allow
! control over fonts, labels, axis labels, line width, scaling, etc.
! Pass in a value of 1 for n to not plot (only draw frame and ticks)

  implicit none

  character(len=1), intent(in) :: panel,frameoff,pltborder,colorbar0
  integer, intent(in) :: n,labincx,labincy,linecolor,ndashes
  real, intent(in) :: aspect,scalelab,xmin,xmax,xinc,ymin,ymax,yinc
  real, intent(in) :: xval(n),yval(n)
  character(len=*), intent(in) :: xlab,ylab

  integer :: i,logy,itickvalq
  real :: dx,dy,tickval,sizelab,xlabx,ylaby,tickvalq
  character(len=20)  :: numbr,numbr2

  integer :: icyc
  real :: dashlen, dist, remain, step, xfac, yfac, asp2, x, y, eps

! Set plot color (black)

  call o_sflush()
  call o_gsplci(10)
  call o_gsfaci(10)
  call o_gstxci(10)
  call o_gslwsc(1.) ! line width

! Scale local working window (0,1,0,1)
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

  call oplot_panel(panel,frameoff,pltborder,colorbar0,aspect,'N')
  call o_set(op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Draw frame

  call o_frstpt(op%fx1,op%fy1)
  call o_vector(op%fx2,op%fy1)
  call o_vector(op%fx2,op%fy2)
  call o_vector(op%fx1,op%fy2)
  call o_vector(op%fx1,op%fy1)

! Specify font # and scale font size to designated plotter coordinates

  call o_sflush()
  call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
  call o_pcsetr('CL',1.)  ! set character line width to 1

  sizelab = scalelab * (op%hp2 - op%hp1)

! Write x axis label

  if (pltborder == 'a' .or. &
      panel     == 'N' .or. &
      panel     == '1' .or. &
      panel     == '2' .or. &
      panel     == '9') then

     xlabx = .5 * (op%fx1 + op%fx2)
     call o_plchhq(xlabx,op%xlaby,trim(xlab),sizelab, 0.,0.)
  endif

! Write y axis label

  if (pltborder == 'a' .or. &
       panel    == 'N' .or. &
       panel    == '1' .or. &
       panel    == '3' .or. &
       panel    == '5') then

     ylaby = .5 * (op%fy1 + op%fy2)
     call o_plchhq(op%ylabx,ylaby,trim(ylab),sizelab,90.,0.)
  endif

! Scale local working window (xmin,xmax,0.,1.)
! to plotter coordinates (op%h1,op%h2,op%vp1,op%vp2)

  call o_set(op%h1,op%h2,op%vp1,op%vp2,xmin,xmax,0.,1.,1)

! Plot and label X-axis ticks

  tickval = nint(xmin/xinc) * xinc
  if (tickval < xmin - .001 * xinc) tickval = tickval + xinc

  do while (tickval < xmax + .001 * xinc)

     if (mod(nint(tickval/xinc),labincx) == 0) then  ! Only for long ticks

        dy = .014

        ! Encode and plot current X tick label

        if (pltborder == 'a' .or. &
            panel     == 'N' .or. &
            panel     == '1' .or. &
            panel     == '2' .or. &
            panel     == '9') then

           if (xinc * labincx >= .999) then
              write (numbr,'(i6)') nint(tickval)
           elseif (xinc * labincx >= .0999) then
              write (numbr,'(f5.1)') tickval
           elseif (xinc * labincx >= .00999) then
              write (numbr,'(f5.2)') tickval
           elseif (xinc * labincx >= .000999) then
              write (numbr,'(f6.3)') tickval
           else
              write (numbr,'(f7.4)') tickval
           endif

           call o_plchhq(tickval,op%xtlaby,trim(adjustl(numbr)),sizelab,0.,0.)
        endif

     else             ! Only for short ticks
        dy = .007
     endif

! Plot current X tick

     call o_frstpt(tickval,op%fy1)
     call o_vector(tickval,op%fy1 + dy)
     call o_frstpt(tickval,op%fy2)
     call o_vector(tickval,op%fy2 - dy)

     tickval = tickval + xinc
  enddo

! Scale local working window (0.,1.,ymin,ymax)
! to plotter coordinates (op%hp1,op%hp2,op%v1,op%v2)

  call o_set(op%hp1,op%hp2,op%v1,op%v2,0.,1.,ymin,ymax,1)

! Plot and label Y-axis ticks

  tickval = nint(ymin/yinc) * yinc

  if ((tickval - ymax) / (ymin - ymax) > 1.001) tickval = tickval + yinc

  do while ((tickval - ymin) / (ymax - ymin) < 1.001)

     if (mod(nint(tickval/yinc),labincy) == 0) then  ! Only for long ticks

        dx = .014

        ! Encode and plot current Y tick label

        if (pltborder == 'a' .or. &
            panel     == 'N' .or. &
            panel     == '1' .or. &
            panel     == '3' .or. &
            panel     == '5') then

           ! Encode current Y tick label

           if (abs(yinc * labincy) >= .999) then
              write (numbr,'(i6)') nint(tickval)
           elseif (yinc * labincy >= .0999) then
              write (numbr,'(f5.1)') tickval
           elseif (yinc * labincy >= .00999) then
              write (numbr,'(f5.2)') tickval
           elseif (yinc * labincy >= .000999) then
              write (numbr,'(f6.3)') tickval
           else

              logy = int(log10(yinc * labincy)) - 2
              tickvalq = tickval * 10. ** (-logy)
              itickvalq = nint(tickvalq)

              ! If significand is at or above 10, reduce

              do while (abs(itickvalq) >= 10)
                 logy = logy + 1
                 tickvalq = tickvalq * .1
                 itickvalq = nint(tickvalq)
              enddo

              ! Use only one of the following 2 lines

            ! write (numbr,'(f4.1)') tickvalq   ! If real significand is required
              write (numbr,'(i2)') itickvalq    ! If integer significand is ok

              write (numbr2,'(i3)') logy
              numbr = trim(adjustl(numbr))
              numbr2 = trim(adjustl(numbr2))

              ! Determine whether significand, power of 10, or both are to be plotted

              if (itickvalq == 1) then
                 numbr = '10:S3:'//trim(numbr2)//'        '
              elseif (itickvalq == -1) then
                 numbr = '-10:S3:'//trim(numbr2)//'        '
              elseif (itickvalq /= 0) then
                 numbr = trim(adjustl(numbr))//'x'//'10:S3:'//trim(adjustl(numbr2))//'        '
              endif

           endif

           ! Plot Y tick label

         ! call o_plchhq(op%ytlabx,tickval,numbr(1:len_trim(numbr)),sizelab,0.,1.)
           call o_plchhq(op%ytlabx,tickval,trim(adjustl(numbr)),sizelab,0.,1.)

        endif ! frameoff/panel

     else                      ! Only for short ticks
        dx = .007
     endif

! Plot current Y tick

     call o_frstpt(op%fx1     ,tickval)
     call o_vector(op%fx1 + dx,tickval)
     call o_frstpt(op%fx2     ,tickval)
     call o_vector(op%fx2 - dx,tickval)

     tickval = tickval + yinc
  enddo

! Scale local working window (xmin,xmax,ymin,ymax)
!  to plotter coordinates (op%h1,op%h2,op%v1,op%v2)

  call o_set(op%h1,op%h2,op%v1,op%v2,xmin,xmax,ymin,ymax,1)

! Plot values

! Set plot color (linecolor)

  call o_sflush()
  call o_gsplci(linecolor)
  call o_gsfaci(linecolor)
  call o_gstxci(linecolor)
  call o_gslwsc(1.) ! line width

  if (panel == '3')then
     if     (linecolor == 16) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.70,'Temperature (deg C)',1.5*sizelab,0.,-1.)
     elseif (linecolor == 109) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.60,'Salinity (g/kg) ',1.5*sizelab,0.,-1.)
     elseif (linecolor == 12) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.50,'TKE*2 (m^2/s^2 x 100) ',1.5*sizelab,0.,-1.)
     elseif (linecolor == 11) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.40,'TKE*L*2 (m^3/s^2 x 100) ',1.5*sizelab,0.,-1.)
     elseif (linecolor == 9) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.30,'U (m/s x 10) ',1.5*sizelab,0.,-1.)
     elseif (linecolor == 8) then
        call o_plchhq(xmin,ymax+(ymax-ymin)*0.20,'V (m/s x 10) ',1.5*sizelab,0.,-1.)
     endif
  endif

  if (ndashes <= 0) then

     call o_frstpt(xval(1),yval(1))
     do i = 2,n
        call o_vector(xval(i),yval(i))
     enddo

  else

     eps = 1.e-6 * (xmax - xmin)

     xfac = (op%h2 - op%h1) / (xmax - xmin)
     yfac = (op%v2 - op%v1) / (ymax - ymin)

     asp2 = (yfac / xfac)**2

     dashlen = (xmax - xmin) / real(2 * ndashes)
     remain = dashlen

     x = xval(1)
     y = yval(1)

     call o_frstpt(x,y)

     icyc = 1
     i = 2

     do while (i < n)

        dist = sqrt((xval(i) - x)**2 + asp2 * (yval(i) - y)**2)

        if (remain > dist + eps) then

           step = dist

           x = xval(i)
           y = yval(i)

           if (icyc > 0) then
              call o_vector(x,y)
           else
              call o_frstpt(x,y)
           endif

           remain = remain - step
           i = i + 1

        else

           step = remain

           x = x + (xval(i) - x) * step / dist
           y = y + (yval(i) - y) * step / dist

           if (icyc > 0) then
              call o_vector(x,y)
           else
              call o_frstpt(x,y)
           endif

           remain = dashlen
           icyc = -icyc

        endif

     enddo

  endif

end subroutine oplot_xy2l

