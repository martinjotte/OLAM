module netcdf_utils

  interface check_temperature_units
     module procedure check_temperature_units_2d
     module procedure check_temperature_units_3d
  end interface check_temperature_units

  interface check_fractional_units
     module procedure check_fractional_units_2d
     module procedure check_fractional_units_3d
  end interface check_fractional_units

contains

!============================================================================

subroutine get_timeunit_offset(tunits, toff8, error)

  use consts_coms, only: r8
  use token_mod,   only: tokenize_ws, tokenize

  implicit none

  character(*), intent(in) :: tunits
  real(r8), intent(out) :: toff8
  logical, intent(out) :: error

  character( 5) :: yymmdd(4)
  character( 5) :: hhmmss(4)
  character(20) :: tokens(5)
  character(20) :: tzstring

  integer :: iyear0, imnth0, idate0, ihour0, imins0, isecs0, ihhmmss0
  integer :: tzhr, tzmn, tzsc, tzsg, tzoff, ntok, itok, slen, zlen
  real    :: rsecs0

  ! CF-convention time units are of the form:
  ! "seconds since 1992-10-8 15:15:42.5 -6:00"
  ! "Seconds since 1970-01-01 00:00"

  error = .false.
  toff8 = 0.0_r8

  ! Split the time unit string by spaces

  call tokenize_ws(tunits,tokens,ntok)

  if (ntok < 3 .or. tokens(2) /= "since") then
     error = .true.
     return
  endif

  ! Parse the third substring for the date

  call tokenize(tokens(3), yymmdd, itok, '-')

  if (itok /= 3) then
     error = .true.
     return
  endif

  read(yymmdd(1),*) iyear0
  read(yymmdd(2),*) imnth0
  read(yymmdd(3),*) idate0

  ! Parse the optional fourth substring for the time (if present)

  ihour0 = 0
  imins0 = 0
  isecs0 = 0

  if (ntok > 3) then

     call tokenize(tokens(4), hhmmss, itok, ':')

     if (itok >= 1) read(hhmmss(1),*) ihour0
     if (itok >= 2) read(hhmmss(2),*) imins0

     ! seconds could be a decimal number
     if (itok >= 3) then
        read(hhmmss(3),*) rsecs0
        isecs0 = nint(rsecs0)
     endif

  endif

  ! Parse the optional fifth substring for time zone offset (if present)

  tzoff = 0

  if (ntok > 4) then

     slen = len_trim(tokens(5))

     if (slen > 0) then

        tzhr = 0
        tzmn = 0
        tzsc = 0
        tzsg = 1

        ! Check for leading minus sign, indicating time offset is west of UTC

        if (tokens(5)(1:1) == '-') then
           tzsg = 1
           tzstring = tokens(5)(2:slen)
        else
           tzsg = -1
           tzstring = tokens(5)(1:slen)
        endif

        call tokenize(tzstring, hhmmss, ntok, ':')

        if (ntok == 1) then

           ! string is hh or hhmm
           zlen = len_trim(hhmmss(1))

           if (zlen < 3) then

              ! string is hh
              read(hhmmss(1),*) tzhr

           elseif (zlen < 5) then

              ! string is hmm or hhmm
              read(hhmmss(1)(     1:zlen-2),*) tzhr
              read(hhmmss(1)(zlen-1:zlen  ),*) tzmn

           elseif (zlen < 7) then

              ! string is hmmss or hhmmss
              read(hhmmss(1)(     1:zlen-4),*) tzhr
              read(hhmmss(1)(zlen-3:zlen-2),*) tzmn
              read(hhmmss(1)(zlen-1:zlen  ),*) tzsc

           endif


        elseif (ntok == 2) then

           ! string is hh:mm
           read(hhmmss(1),*) tzhr
           read(hhmmss(2),*) tzmn

        elseif (ntok == 3) then

           ! string is hh:mm:ss
           read(hhmmss(1),*) tzhr
           read(hhmmss(2),*) tzmn
           read(hhmmss(3),*) tzsc

        endif

        tzoff = tzsg * (tzhr * 3600 + tzmn * 60 + tzsc)

     endif
  endif

  ! Compute number of seconds past 1 Jan 1900 00Z

  ihhmmss0 = ihour0 * 10000 + imins0 * 100 + isecs0

  call date_abs_secs2(iyear0, imnth0, idate0, ihhmmss0, toff8)

  ! Include the time zone offset

  toff8 = toff8 + real(tzoff,r8)

end subroutine get_timeunit_offset

!============================================================================

subroutine ncf_times_to_s1900(tstring, ntimes, dtimes, error)

  use consts_coms, only: r8
  implicit none

  character(*), intent(in   ) :: tstring
  integer,      intent(in   ) :: ntimes
  real(r8),     intent(inout) :: dtimes(ntimes)
  logical,      intent(  out) :: error

  real(r8)                    :: toff8

  ! Convert time units base to time offset from 0Z Jan 01 1900

  call get_timeunit_offset(tstring,toff8,error)

  ! Convert time to seconds since 0Z Jan 01 1900

  if ( any(tstring(1:1) == ['s','S']) ) then
     dtimes(1:ntimes) = toff8 + dtimes(1:ntimes)
  elseif ( any(tstring(1:1) == ['m','M']) ) then
     dtimes(1:ntimes) = toff8 + dtimes(1:ntimes) * 60._r8
  elseif ( any(tstring(1:1) == ['h','H']) ) then
     dtimes(1:ntimes) = toff8 + dtimes(1:ntimes) * 3600._r8
  elseif ( any(tstring(1:1) == ['d','D']) ) then
     dtimes(1:ntimes) = toff8 + dtimes(1:ntimes) * 86400._r8
  else
     error = .true.
  endif

end subroutine ncf_times_to_s1900

!============================================================================

subroutine ncftime_to_s1900(tstring, dtime, error)

  use consts_coms, only: r8
  implicit none

  character(*), intent(in   ) :: tstring
  real(r8),     intent(inout) :: dtime
  logical,      intent(  out) :: error
  real(r8)                    :: toff8

  ! Convert time units base to time offset from 0Z Jan 01 1900

  call get_timeunit_offset(tstring,toff8,error)

  ! Convert time to seconds since 0Z Jan 01 1900

  if ( tstring(1:1) == 's') then
     dtime = toff8 + dtime
  elseif ( tstring(1:1) == 'm') then
     dtime = toff8 + dtime * 60._r8
  elseif ( tstring(1:1) == 'h') then
     dtime = toff8 + dtime * 3600._r8
  elseif ( tstring(1:1) == 'd') then
     dtime = toff8 + dtime * 86400._r8
  else
     error = .true.
  endif

end subroutine ncftime_to_s1900

!============================================================================

subroutine check_temperature_units_2d(units, nx, ny, field)

  use string_lib,  only: lowercase, strip_char
  use consts_coms, only: t00

  implicit none

  character(*), intent(inout) :: units
  integer,      intent(in)    :: nx, ny
  real,         intent(inout) :: field(nx,ny)
  integer                     :: i, j

  if (len_trim(units) < 1) return

  call lowercase (units)
  call strip_char(units, ' ')

  if ( units == 'c'       .or. units == 'degreesc'       .or. &
       units == 'celsius' .or. units == 'degreescelsius' ) then

     !$omp parallel do private(i)
     do j = 1, ny
        do i = 1, nx
           if (field(i,j) > -999.) then
              field(i,j) = field(i,j) + t00
           endif
        enddo
     enddo
     !$omp end parallel do

  elseif ( units /= 'k'      .and. units /= 'degreesk'      .and. &
           units /= 'kelvin' .and. units /= 'degreeskelvin' ) then

     write(*,*) "Unknown temperature units ", trim(units)
     stop

  endif

end subroutine check_temperature_units_2d

!============================================================================

subroutine check_temperature_units_3d(units, nx, ny, nz, field)

  use string_lib,  only: lowercase, strip_char
  use consts_coms, only: t00

  implicit none

  character(*), intent(inout) :: units
  integer,      intent(in)    :: nx, ny, nz
  real,         intent(inout) :: field(nx,ny,nz)
  integer                     :: i, j, k

  if (len_trim(units) < 1) return

  call lowercase (units)
  call strip_char(units, ' ')

  if ( units == 'c'       .or. units == 'degreesc'       .or. &
       units == 'celsius' .or. units == 'degreescelsius' ) then

     !$omp parallel do private(j,i)
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              if (field(i,j,k) > -999.) then
                 field(i,j,k) = field(i,j,k) + t00
              endif
           enddo
        enddo
     enddo
     !$omp end parallel do

  elseif ( units /= 'k'      .and. units /= 'degreesk'      .and. &
           units /= 'kelvin' .and. units /= 'degreeskelvin' ) then

     write(*,*) "Unknown temperature units ", trim(units)
     stop

  endif

end subroutine check_temperature_units_3d

!============================================================================

subroutine check_fractional_units_2d(units, nx, ny, field)

  use string_lib,  only: lowercase, strip_chars

  implicit none

  character(*), intent(inout) :: units
  integer,      intent(in)    :: nx, ny
  real,         intent(inout) :: field(nx,ny)
  integer                     :: i, j

  if (len_trim(units) < 1) return

  call lowercase  (units)
  call strip_chars(units, [ ' ', '*', '.', '[', ']', '(', ')' ])

  if ( units == '%' .or. units == 'percent' ) then

     !$omp parallel do private(i)
     do j = 1, ny
        do i = 1, nx
           if (field(i,j) > -999.) then
              field(i,j) = 0.01 * field(i,j)
           endif
        enddo
     enddo
     !$omp end parallel do

  elseif ( units /= '0-1' .and. units /= 'fraction' ) then

     write(*,*) "Unknown fractional units ", trim(units)
     stop

  endif

end subroutine check_fractional_units_2d

!============================================================================

subroutine check_fractional_units_3d(units, nx, ny, nz, field)

  use string_lib,  only: lowercase, strip_chars

  implicit none

  character(*), intent(inout) :: units
  integer,      intent(in)    :: nx, ny, nz
  real,         intent(inout) :: field(nx,ny,nz)
  integer                     :: i, j, k

  if (len_trim(units) < 1) return

  call lowercase  (units)
  call strip_chars(units, [ ' ', '*', '.', '[', ']', '(', ')' ])

  if ( units == '%' .or. units == 'percent' ) then

     !$omp parallel do private(j,i)
     do k = 1, nz
        do j = 1, ny
           do i = 1, nx
              if (field(i,j,k) > -999.) then
                 field(i,j,k) = 0.01 * field(i,j,k)
              endif
           enddo
        enddo
     enddo
     !$omp end parallel do

  elseif ( units /= '0-1' .and. units /= 'fraction' ) then

     write(*,*) "Unknown fractional units ", trim(units)
     stop

  endif

end subroutine check_fractional_units_3d

!============================================================================

end module netcdf_utils
