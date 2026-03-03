module analysis_lib

contains

!===============================================================================

subroutine gdtost_ll_to_atm(field, nx, ny, swlon, swlat, dlon, dlat, &
                            inproj, r2d, d2d, i2d, b2d, plats, irev_ns, ilevel)

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8
  use               mem_ijtabs,      only: jtab_w, jtw_prog
  use               mem_grid,        only: mwa, glatw, glonw

  implicit none

  real,                  intent(out) :: field(mwa)
  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)
  integer,     optional, intent(in)  :: ilevel

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  integer :: j, iw

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call gdtost_ll( nx, ny, glonw(iw), glatw(iw), field(iw), swlon, swlat, dlon, dlat, &
                     inproj, r2d, d2d, i2d, b2d, plats, irev_ns, ilevel )

  enddo
  !$omp end parallel do

end subroutine gdtost_ll_to_atm

!===============================================================================

subroutine gdtost_ll_to_sfc( field, nx, ny, swlon, swlat, dlon, dlat, &
                             inproj, r2d, d2d, i2d, b2d, plats, irev_ns, ilevel, &
                             primary_only, iland, isea, ilake )

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

  use mem_sfcg,  only: mwsfc, sfcg, itab_wsfc
  use mem_land,  only: mland, omland
  use mem_sea,   only: msea, omsea
  use mem_lake,  only: mlake, omlake
  use misc_coms, only: iparallel
  use mem_para,  only: myrank

  implicit none

  real,      contiguous, intent(out) :: field(:)
  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)
  integer,     optional, intent(in)  :: ilevel

  logical,     optional, intent(in)  :: primary_only
  logical,     optional, intent(in)  :: iland
  logical,     optional, intent(in)  :: isea
  logical,     optional, intent(in)  :: ilake

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  integer :: is, iwsfc, iend
  logical :: prim_only

  prim_only = .true.
  if (present(primary_only)) prim_only = primary_only

  if ( count( [present(iland), present(isea), present(ilake)] ) > 1 ) then
     write(*,*) "Error in gdtost_ll_to_sfc:"
     write(*,*) "only one of iland, isea, or ilake can be specified"
     stop
  endif

  if (present(iland)) then
     iend = mland
     if (size(field) /= mland) then
        write(*,*) "Error in gdtost_ll_to_sfc:"
        write(*,*) "field must be dimensioned to mland if iland is specified."
        stop
     endif
  elseif (present(ilake)) then
     iend = mlake
     if (size(field) /= mlake) then
        write(*,*) "Error in gdtost_ll_to_sfc:"
        write(*,*) "field must be dimensioned to mlake if ilake is specified."
        stop
     endif
  elseif (present(isea)) then
     iend = msea
     if (size(field) /= msea) then
        write(*,*) "Error in gdtost_ll_to_sfc:"
        write(*,*) "field must be dimensioned to msea if isea is specified."
        stop
     endif
  else
     iend = mwsfc
     if (size(field) /= mwsfc) then
        write(*,*) "Error in gdtost_ll_to_sfc:"
        write(*,*) "field must be dimensioned to mwsfc if iland, ilake, or isea is not specified."
        stop
     endif
  endif

  !$omp parallel do private(iwsfc)
  do is = 2, iend

     iwsfc = is
     if (present(iland)) iwsfc = is + omland
     if (present(ilake)) iwsfc = is + omlake
     if (present(isea )) iwsfc = is + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (prim_only) then
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle
     endif

     call gdtost_ll( nx, ny, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), field(is), swlon, swlat, &
                     dlon, dlat, inproj, r2d, d2d, i2d, b2d, plats, irev_ns, ilevel )

  enddo
  !$omp end parallel do

end subroutine gdtost_ll_to_sfc

!===============================================================================

subroutine gdtost_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                     inproj, r2d, d2d, i2d, b2d, plats, irev_ns, ilevel)

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

! import, only: quad_avg
  implicit none

  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: glon, glat
  real,                  intent(out) :: staval

  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)
  integer,     optional, intent(in)  :: ilevel

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  !  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  !  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  !  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(NX,NY),WHERE
  !  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  !  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  !  AND STATION COLUMN.
  !
  !  ON INPUT, VALUES OF -999 AND SMALLER INDICATE MISSING DATA.
  !  ON OUTPUT, STAVAL OF -1.E30 INDICATES MISSING DATA.

  integer :: level

  level = 2
  if (present(ilevel)) level = ilevel

  if (level <= 0) then
     call gdtost_nearest_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                            inproj, r2d, d2d, i2d, b2d, plats, irev_ns)
  elseif (level == 1) then
     call gdtost_bilinear_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                             inproj, r2d, d2d, i2d, b2d, plats, irev_ns)
  else
     call gdtost_quadratic_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                              inproj, r2d, d2d, i2d, b2d, plats, irev_ns)
  endif

end subroutine gdtost_ll

!===============================================================================

subroutine gdtost_quadratic_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                               inproj, r2d, d2d, i2d, b2d, plats, irev_ns)

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

! import, only: quad_avg
  implicit none

  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: glon, glat
  real,                  intent(out) :: staval

  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  !  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  !  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  !  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(NX,NY),WHERE
  !  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  !  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  !  AND STATION COLUMN.
  !
  !  ON INPUT, VALUES OF -999 AND SMALLER INDICATE MISSING DATA.
  !  ON OUTPUT, STAVAL OF -1.E30 INDICATES MISSING DATA.

  integer        :: iy1,ix1,ii,i,jj,j,iin,iio,jjo
  real           :: yy,xx,stax,stay
  real           :: r(4), scr(4)

  i = count( [present(r2d), present(d2d), present(i2d), present(b2d)] )
  if (i /= 1) stop 'Must have one input array to subroutine gdtost'

  if (inproj == 2) then

     if (.not. present(plats)) stop 'gdtost: plats array must be input with inproj=2'

     if (glat < plats(1)) then

        iy1 = -1
        yy = 3.0 - (glat - plats(1)) / (plats(1) - plats(2))

     elseif (glat >= plats(ny)) then

        iy1 = ny - 1
        yy = 2.0 + (glat - plats(ny)) / (plats(ny) - plats(ny-1))

     else

        ! estimate latitude index assuming uniform spacing of plats

        iy1 = floor( real(ny) * (glat + 90.0) / 180. )
        iy1 = max(0, min(ny-2, iy1))

        ! find correct latitude index

        if (glat < plats(iy1+1)) then
           do while ( glat < plats(iy1+1) )
              iy1 = iy1 - 1
           enddo
        elseif (glat >= plats(iy1+2)) then
           do while ( glat >= plats(iy1+2) )
              iy1 = iy1 + 1
           enddo
        endif

        yy = 2. + (glat - plats(iy1+1)) / (plats(iy1+2) - plats(iy1+1))

     endif

  else

     stay = (glat - swlat) / dlat
     iy1 = floor(stay)
     yy = stay - real(iy1 - 2)

  endif

  stax = (glon - swlon) / dlon
  ix1 = floor(stax)
  xx = stax - real(ix1 - 2)

  do i = ix1, ix1+3

     ii  = i - ix1 + 1
     iin = modulo(i-1,nx) + 1

     ! Initialize data to missing values
     r(1:4) = -1.e30

     do j = iy1, iy1+3
        jj = j - iy1 + 1

        if (j < 1) then

           if (inproj==1 .and. abs(swlat + 90.) < .1 * dlat) then
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 - j
           else
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 1 - j
           endif

        elseif (j > ny) then

           if (inproj==1 .and. abs(swlat + 90.) < .1 * dlat) then
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 * ny - j
           else
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 * ny - j + 1
           endif

        else

           iio = iin
           jjo = j

        endif

        if (present(irev_ns)) then
           if (irev_ns) jjo = ny - jjo + 1
        endif

        if     (present(r2d)) then
           if (r2d(iio,jjo) > -999.)    r(jj) = r2d(iio,jjo)
        elseif (present(d2d)) then
           if (d2d(iio,jjo) > -999._r8) r(jj) = d2d(iio,jjo)
        elseif (present(i2d)) then
           if (i2d(iio,jjo) > -999)     r(jj) = i2d(iio,jjo)
        elseif (present(b2d)) then
           if (b2d(iio,jjo) > -127_i1)  r(jj) = b2d(iio,jjo)
        endif

     enddo

     scr(ii) = quad_avg(yy, 1., 2., 3., 4., r(1), r(2), r(3), r(4) )

  enddo

  staval = quad_avg(xx, 1., 2., 3., 4., scr(1),scr(2),scr(3),scr(4))

end subroutine gdtost_quadratic_ll

!===============================================================================

subroutine gdtost_bilinear_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                     inproj, r2d, d2d, i2d, b2d, plats, irev_ns)

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

! import, only: quad_avg
  implicit none

  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: glon, glat
  real,                  intent(out) :: staval

  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  !  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  !  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  !  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(NX,NY),WHERE
  !  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  !  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  !  AND STATION COLUMN.
  !
  !  ON INPUT, VALUES OF -999 AND SMALLER INDICATE MISSING DATA.
  !  ON OUTPUT, STAVAL OF -1.E30 INDICATES MISSING DATA.

  integer        :: iy1,ix1,ii,i,jj,j,iin,iio,jjo
  real           :: yy,xx,stax,stay
  real           :: r(4), scr(4)

  i = count( [present(r2d), present(d2d), present(i2d), present(b2d)] )
  if (i /= 1) stop 'Must have one input array to subroutine gdtost'

  if (inproj == 2) then

     if (.not. present(plats)) stop 'gdtost: plats array must be input with inproj=2'

     if (glat < plats(1)) then

        iy1 = 0
        yy = 1. + (glat - plats(1)) / (plats(2) - plats(1))

     elseif (glat >= plats(ny)) then

        iy1 = ny
        yy = (glat - plats(ny)) / (plats(ny) - plats(ny-1))

     else

        ! estimate latitude index assuming uniform spacing of plats

        iy1 = floor( real(ny) * (glat + 90.0) / 180. )
        iy1 = max(0, min(ny-1, iy1))

        ! find correct latitude index

        if (glat < plats(iy1)) then
           do while ( glat < plats(iy1) )
              iy1 = iy1 - 1
           enddo
        elseif (glat >= plats(iy1+1)) then
           do while ( glat >= plats(iy1+1) )
              iy1 = iy1 + 1
           enddo
        endif

        yy = (glat - plats(iy1)) / (plats(iy1+1) - plats(iy1))

     endif

  else

     stay = (glat - swlat) / dlat
     iy1 = floor(stay) + 1
     yy = stay - real(iy1 - 1)

  endif

  stax = (glon - swlon) / dlon
  ix1 = floor(stax) + 1
  xx = stax - real(ix1 - 1)

  do i = ix1, ix1+1

     ii  = i - ix1 + 1
     iin = modulo(i-1,nx) + 1

     ! Initialize data to missing values
     r(1:4) = -1.e30

     do j = iy1, iy1+1
        jj = j - iy1 + 1

        if (j < 1) then

           if (inproj==1 .and. abs(swlat + 90.) < .1 * dlat) then
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 - j
           else
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 1 - j
           endif

        elseif (j > ny) then

           if (inproj==1 .and. abs(swlat + 90.) < .1 * dlat) then
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 * ny - j
           else
              iio = mod(iin + nx/2 - 1, nx) + 1
              jjo = 2 * ny - j + 1
           endif

        else

           iio = iin
           jjo = j

        endif

        if (present(irev_ns)) then
           if (irev_ns) jjo = ny - jjo + 1
        endif

        if     (present(r2d)) then
           if (r2d(iio,jjo) > -999.)    r(jj) = r2d(iio,jjo)
        elseif (present(d2d)) then
           if (d2d(iio,jjo) > -999._r8) r(jj) = d2d(iio,jjo)
        elseif (present(i2d)) then
           if (i2d(iio,jjo) > -999)     r(jj) = i2d(iio,jjo)
        elseif (present(b2d)) then
           if (b2d(iio,jjo) > -127_i1)  r(jj) = b2d(iio,jjo)
        endif

     enddo

     if (r(1) <= -1.e30 .and. r(2) <= -1.e30) then
        scr(ii) = -1.e30
     elseif (r(1) <= -1.e30) then
        scr(ii) = r(2)
     elseif (r(2) <= -1.e30) then
        scr(ii) = r(1)
     else
        scr(ii) = r(1) + yy * (r(2) - r(1))
     endif

  enddo

  if (scr(1) <= -1.e30 .and. scr(2) <= -1.e30) then
     staval = -1.e30
  elseif (scr(1) <= -1.e30) then
     staval = scr(2)
  elseif (scr(2) <= -1.e30) then
     staval = scr(1)
  else
     staval = scr(1) + xx * (scr(2) - scr(1))
  endif

end subroutine gdtost_bilinear_ll

!===============================================================================

subroutine gdtost_nearest_ll(nx, ny, glon, glat, staval, swlon, swlat, dlon, dlat, &
                             inproj, r2d, d2d, i2d, b2d, plats, irev_ns)

  use, intrinsic :: iso_fortran_env, only: r8=>real64, i1=>int8

  implicit none

  integer,               intent(in)  :: nx, ny
  real,                  intent(in)  :: glon, glat
  real,                  intent(out) :: staval

  real,                  intent(in)  :: swlon, swlat ! lon and lat of lower left point
  real,                  intent(in)  :: dlon, dlat   ! lon and lat spacing
  integer,               intent(in)  :: inproj       ! 1 if lat-lon, 2 if gaussian
  real,        optional, intent(in)  :: plats(nx)    ! list of gaussian lats
  logical,     optional, intent(in)  :: irev_ns      ! for files with reversed lats (N->S)

  real,        optional, intent(in)  :: r2d(nx,ny)
  real(r8),    optional, intent(in)  :: d2d(nx,ny)
  integer,     optional, intent(in)  :: i2d(nx,ny)
  integer(i1), optional, intent(in)  :: b2d(nx,ny)

  !  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  !  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  !  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(NX,NY),WHERE
  !  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  !  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  !  AND STATION COLUMN.
  !
  !  ON INPUT, VALUES OF -999 AND SMALLER INDICATE MISSING DATA.
  !  ON OUTPUT, STAVAL OF -1.E30 INDICATES MISSING DATA.

  integer :: i, ix, iy
  real    :: stax, stay

  i = count( [present(r2d), present(d2d), present(i2d), present(b2d)] )
  if (i /= 1) stop 'Must have one input array to subroutine gdtost'

  if (inproj == 2) then

     if (.not. present(plats)) stop 'gdtost: plats array must be input with inproj=2'

     if (glat < 0.5*(plats(1)+plats(2))) then
        iy = 1
     elseif (glat > 0.5*(plats(ny)+plats(ny-1))) then
        iy = ny
     else
        iy = minloc( abs(glat - plats(:)), 1)
     endif

  else

     stay = (glat - swlat) / dlat
     iy   = max(1, min(ny, nint(stay) + 1))

  endif

  stax = (glon - swlon) / dlon
  ix   = modulo(nint(stax),nx) + 1

  if (present(irev_ns)) then
     if (irev_ns) iy = ny - iy + 1
  endif

  if     (present(r2d)) then
     staval = r2d(ix,iy)
  elseif (present(d2d)) then
     staval = d2d(ix,iy)
  elseif (present(i2d)) then
     staval = i2d(ix,iy)
  elseif (present(b2d)) then
     staval = b2d(ix,iy)
  endif

end subroutine gdtost_nearest_ll

!===============================================================================

subroutine gdtost(a,ix,iy,stax,stay,staval)

! import, only: quad_avg
  implicit none

  integer, intent(in)  :: ix,iy
  real,    intent(in)  :: a(ix,iy), stax, stay
  real,    intent(out) :: staval

  !  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  !  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  !  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
  !  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  !  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  !  AND STATION COLUMN.
  !
  !  ON INPUT, VALUES OF 1.E30 AND LARGER OR -999 AND SMALLER
  !  INDICATE MISSING DATA.
  !  ON OUTPUT, STAVAL OF -1.0E30 INDICATES MISSING DATA.

  integer :: iy1,iy2,ix1,ix2,ii,i,jj,j
  real    :: fiym2,fixm2,yy,xx
  real    :: r(4), scr(4)

  iy1 = int(stay) - 1
  iy2 = iy1 + 3

  ix1 = int(stax)-1
  ix2 = ix1+3

  fiym2=real(iy1 - 1)
  fixm2=real(ix1 - 1)

  yy = stay - fiym2
  xx = stax - fixm2

  do i  = ix1, ix2
     ii = i - ix1 + 1

     if (i < 1 .or. i > ix) then

        scr(ii) = -1.e30

     else

        do j  = iy1, iy2
           jj = j - iy1 + 1

           if (j < 1 .or. j > iy) then
              r(jj) = -1.e30
           elseif (a(i,j) <= -999.) then
              r(jj) = -1.e30
           else
              r(jj) = a(i,j)
           endif
        enddo

        scr(ii) = quad_avg(yy, 1., 2., 3., 4., r(1), r(2), r(3), r(4) )

     endif

  enddo

  staval = quad_avg(xx,  1., 2., 3., 4., scr(1),scr(2),scr(3),scr(4))

end subroutine gdtost

!===============================================================================

real function quad_avg ( x, x0, x1, x2, x3, y0, y1, y2, y3 )

! import, only: quad
  implicit none

  real, intent(in) :: x, x0, x1, x2, x3, y0, y1, y2, y3
  real             :: q1, q2

  !  Since there are three points required for a quadratic, we compute it twice
  !  (once with x0, x1, x2 and once with x1, x2, x3), and then average those.
  !  This will reduce overshoot. The "x" point is where we are interpolating to.

  q1 = quad( x, x0, x1, x2,     y0, y1, y2     )
  q2 = quad( x,     x1, x2, x3,     y1, y2, y3 )

  if (abs(q1) < 1.e30 .and. abs(q2) < 1.e30) then

     quad_avg = ( q1 * ( x2 - x  ) &
                + q2 * ( x  - x1 ) ) / ( x2 - x1 )

  elseif (abs(q1) < 1.e30) then

     quad_avg = q1

  elseif (abs(q2) < 1.e30) then

     quad_avg = q2

  else

     quad_avg = -1.e30

  endif

end function quad_avg

!===============================================================================

real function quad ( x, x0, x1, x2, y0, y1, y2 )

  implicit none

  real    :: x , x0, x1, x2, y0, y1, y2, ymax, ymin
  integer :: n
  logical :: mask(3)

  mask = abs( [y0, y1, y2] ) < 1.e30

  n = count( mask )

  if (n == 0) then
     quad = -1.e30
     return
  elseif (n == 1 .or. n == 2) then
     quad = sum( [y0, y1, y2], mask = mask ) / real(n)
     return
  endif

  !  Lagrange = sum     prod    ( x  - xj )
  !             i=0,n ( j=0,n    ---------  * yi )
  !                     j<>i    ( xi - xj )

  !  For a quadratic, in the above equation, we are setting n=2. Three points
  !  required for a quadratic, points x0, x1, x2 (hence n=2).

  quad = (x-x1)*(x-x2)*y0 / ( (x1-x0)*(x2-x0) ) &
       - (x-x0)*(x-x2)*y1 / ( (x1-x0)*(x2-x1) ) &
       + (x-x0)*(x-x1)*y2 / ( (x2-x0)*(x2-x1) )

  ! Should we bound the quadratic interpolants?

  ymax = max(y0, y1, y2)
  ymin = min(y0, y1, y2)

  quad = min( max( quad, ymin ), ymax )

end function quad

!===============================================================================

end module analysis_lib
