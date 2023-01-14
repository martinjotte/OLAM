subroutine date_abs_secs(cdate,seconds)

use consts_coms, only: r8

implicit none

character(len=14), intent(in) :: cdate
real(r8), intent(out) :: seconds

! Compute number of seconds past 1 January 1900 12:00 am
! The format of hour1 is hhmmss

integer :: nyears,ndays
real(r8) :: s1,s2,s3,s4
integer, external :: julday

integer :: year1, month1, date1, hour1

call date_unmake_big(year1,month1,date1,hour1,cdate)

nyears = year1 - 1900
ndays = nyears * 365 + (nyears-1)/4 + julday(month1,date1,year1) - 1

! Integer arithmetic in ndays calculation correctly treats 1900 as
! non-leapyear and 2000 as leapyear.  Correct for non-leapyears
! 2100, 2200, and 2300 if necessary.

if (year1 > 2100) ndays = ndays - 1
if (year1 > 2200) ndays = ndays - 1
if (year1 > 2300) ndays = ndays - 1

s1 = real(ndays,r8) * 86400._r8
s2 = real(hour1/10000,r8) * 3600._r8
s3 = real(mod(hour1,10000) / 100,r8) * 60._r8
s4 = real(mod(hour1,100),r8)

seconds = s1 + s2 + s3 +s4

return
end subroutine date_abs_secs

!===============================================================================

subroutine date_abs_secs2(year1,month1,date1,hour1,seconds)

use consts_coms, only: r8

implicit none

integer, intent(in) :: year1,month1,date1,hour1
real(r8), intent(out) :: seconds

! Compute number of seconds past 1 January 1900 12:00 am
! The format of hour1 is hhmmss

integer :: nyears,ndays
real(r8) :: s1,s2,s3,s4
integer, external :: julday

nyears = year1 - 1900
ndays = nyears * 365 + (nyears-1)/4 + julday(month1,date1,year1) - 1

! Integer arithmetic in ndays calculation correctly treats 1900 as
! non-leapyear and 2000 as leapyear.  Correct for non-leapyears
! 2100, 2200, and 2300 if necessary.

if (year1 > 2100) ndays = ndays - 1
if (year1 > 2200) ndays = ndays - 1
if (year1 > 2300) ndays = ndays - 1

s1 = real(ndays,r8) * 86400._r8
s2 = real(hour1/10000,r8) * 3600._r8
s3 = real(mod(hour1,10000) / 100,r8) * 60._r8
s4 = real(mod(hour1,100),r8)

seconds = s1 + s2 + s3 + s4

return
end subroutine date_abs_secs2

!===============================================================================

subroutine date_subtract(cdate1,cdate2,tinc,tunits)

use consts_coms, only: r8

implicit none

character(len=14), intent(in) :: cdate1, cdate2
real, intent(out) :: tinc
character(len=1), intent(in) :: tunits

! add (or subracts) a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(r8) :: secs1,secs2
real :: ttinc

call date_abs_secs(cdate1,secs1)
call date_abs_secs(cdate2,secs2)

! convert time to requested unit

ttinc = secs2 - secs1
if (tunits == 's') tinc = ttinc
if (tunits == 'm') tinc = ttinc / 60.
if (tunits == 'h') tinc = ttinc / 3600.
if (tunits == 'd') tinc = ttinc / 86400.

return
end subroutine date_subtract

!===============================================================================

subroutine date_secs_ymdt(seconds,iyear1,imonth1,idate1,ihour1)

use consts_coms, only: r8

implicit none

real(r8), intent(in) :: seconds
integer, intent(out) :: iyear1,imonth1,idate1,ihour1

! compute real time given number of seconds past 1 January 1900 12:00 am

real(r8) :: s1, secs_in_year, secs_in_month
integer  :: iyear, ileap, imonth, ihour, iminute, isecond, days_in_month
integer  :: mondays(12) = [ 31,28,31,30,31,30,31,31,30,31,30,31 ]

! Get what year it is

s1 = seconds

do iyear = 1900,10000

   if ( (mod(iyear,400) == 0) .or. &
        (mod(iyear,  4) == 0 .and. mod(iyear,100) /= 0) ) then
      ileap = 1
   else
      ileap = 0
   endif

   secs_in_year = real((365+ileap)*86400,r8)

   if (s1 - secs_in_year < 0._r8) then
      iyear1 = iyear
      exit
   endif

   s1 = s1 - secs_in_year
enddo

! s1 is now number of secs into the year; get month

do imonth = 1,12
   days_in_month = mondays(imonth)
   if (imonth == 2) days_in_month = days_in_month + ileap

   secs_in_month = real(days_in_month*86400,r8)

   if (s1 - secs_in_month < 0._r8) then
      imonth1 = imonth
      exit
   endif

   s1 = s1 - secs_in_month
enddo

! s1 is now number of secs into the month
!   Get date and time

idate1 = int(s1 / 86400._r8)
s1 = s1 - real(idate1*86400,r8)
idate1 = idate1 + 1 ! Since date starts at 1

ihour = int(s1 / 3600._r8)
s1 = s1 - real(ihour*3600,r8)

iminute = int(s1 / 60._r8)
s1 = s1 - real(iminute*60,r8)

isecond = int(s1)
ihour1 = ihour * 10000 + iminute * 100 + isecond

return
end subroutine date_secs_ymdt

!===============================================================================

subroutine date_make_big(iyear,imonth,idate,ihour,cdate)

implicit none

integer, intent(in) :: iyear,imonth,idate,ihour
character(len=14), intent(out) ::  cdate

write(cdate(1:4),10) iyear
write(cdate(5:6),11) imonth
write(cdate(7:8),11) idate
write(cdate(9:14),12) ihour
10 format (i4.4)
11 format (i2.2)
12 format (i6.6)

return
end subroutine date_make_big

!===============================================================================

subroutine date_unmake_big(iyear,imonth,idate,ihour,cdate)

implicit none

integer, intent(out) :: iyear,imonth,idate,ihour
character(len=14), intent(in) :: cdate

read(cdate(1:4),10) iyear
read(cdate(5:6),11) imonth
read(cdate(7:8),11) idate
read(cdate(9:14),12) ihour
10 format (i4)
11 format (i2)
12 format (i6)

return
end subroutine date_unmake_big

!===============================================================================

subroutine dintsort(ni,chnums,cstr)

implicit none

integer :: ni
character(len=14) :: chnums(*)
character(len=*) :: cstr(*)

! sort an array of character strings by an associated character field

character(len=200) :: cscr
character(len=14) :: mini,nscr

integer :: n,nm,nmm

do n = 1,ni
   mini = '99999999999999'
   do nm = n,ni

      if (chnums(nm) < mini) then
         nmm = nm
         mini = chnums(nm)
      endif
   enddo

   nscr        = chnums(n)
   chnums(n)   = chnums(nmm)
   chnums(nmm) = nscr

   cscr      = cstr(n)
   cstr(n)   = cstr(nmm)
   cstr(nmm) = cscr
enddo

return
end subroutine dintsort

!===============================================================================

subroutine unique_dint(n1,ca1)

implicit none

integer, intent(inout) :: n1
character(len=14), intent(inout) :: ca1(*)

integer :: n,nt,nn

! reduce an array to get rid of duplicate entries

nt = n1
10 continue
do n = 2,nt
   if (ca1(n) == ca1(n-1)) then
      do nn = n,nt
         ca1(nn-1) = ca1(nn)
      enddo
      nt = nt - 1
      goto 10
   endif
enddo
n1 = nt

return
end subroutine unique_dint

!===============================================================================

integer function julday(imonth,iday,iyear)

implicit none

integer, intent(in) :: imonth,iday,iyear
integer             :: nfeb
logical             :: isleap

isleap = ( (mod(iyear,400) == 0) .or. &
           (mod(iyear,100) /= 0 .and. mod(iyear,4) == 0) )

if (isleap) then
   nfeb = 29
else
   nfeb = 28
endif

! compute the julian day from a normal date

julday= iday  &
      + min(1,max(0,imonth-1))*31   &
      + min(1,max(0,imonth-2))*nfeb &
      + min(1,max(0,imonth-3))*31   &
      + min(1,max(0,imonth-4))*30   &
      + min(1,max(0,imonth-5))*31   &
      + min(1,max(0,imonth-6))*30   &
      + min(1,max(0,imonth-7))*31   &
      + min(1,max(0,imonth-8))*31   &
      + min(1,max(0,imonth-9))*30   &
      + min(1,max(0,imonth-10))*31  &
      + min(1,max(0,imonth-11))*30  &
      + min(1,max(0,imonth-12))*31

return
end function julday

!===============================================================================

subroutine date_add_to_big8(cindate,tinc8,tunits,coutdate)

use consts_coms, only: r8

implicit none

character(len=14), intent(in) :: cindate
real(r8), intent(in) :: tinc8
character(len=1), intent(in) :: tunits
character(len=14), intent(out) :: coutdate

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(r8) :: ttinc,secs
integer :: inyear,inmonth,indate,inhour  &
          ,outyear,outmonth,outdate,outhour

! convert input time to seconds

ttinc = tinc8
if (tunits == 'm') ttinc = tinc8 * 60._r8
if (tunits == 'h') ttinc = tinc8 * 3600._r8
if (tunits == 'd') ttinc = tinc8 * 86400._r8

call date_unmake_big(inyear,inmonth,indate,inhour,cindate)
call date_abs_secs2(inyear,inmonth,indate,inhour,secs)

secs = secs + ttinc

call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)
call date_make_big(outyear,outmonth,outdate,outhour,coutdate)

return
end subroutine date_add_to_big8

!===============================================================================

subroutine date_add_to8(inyear,inmonth,indate,inhour,tinc8,tunits  &
                       ,outyear,outmonth,outdate,outhour)

use consts_coms, only: r8

implicit none

integer, intent(in) :: inyear,inmonth,indate,inhour
real(r8), intent(in) :: tinc8
character(len=1), intent(in) :: tunits
integer, intent(out) :: outyear,outmonth,outdate,outhour

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(r8) :: ttinc,secs

! convert input time to seconds

ttinc = tinc8
if (tunits == 'm') ttinc = tinc8 * 60._r8
if (tunits == 'h') ttinc = tinc8 * 3600._r8
if (tunits == 'd') ttinc = tinc8 * 86400._r8

call date_abs_secs2(inyear,inmonth,indate,inhour,secs)

secs = secs + ttinc

call date_secs_ymdt(secs,outyear,outmonth,outdate,outhour)

return
end subroutine date_add_to8

!===============================================================================

subroutine makefnam(fname, prefix, ctime, ftype, post, fmt)

  use misc_coms, only: simtime

  implicit none

! Creates standard timestamped filename from date/time

  character(len=*), intent(out) :: fname
  character(len=*), intent(in)  :: prefix, ftype, post, fmt
  type(simtime),    intent(in)  :: ctime

  character(len=40) :: dstring
  integer           :: ihhmmss
  integer, external :: timetohhmmss

  ihhmmss = timetohhmmss(ctime%time)

  write(dstring,'(3a,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)') '-', trim(ftype),  &
       '-', ctime%year, '-', ctime%month, '-', ctime%date, '-', ihhmmss

  fname = trim(prefix)//trim(dstring)
  if (post(1:1) /= '$') then
     fname = trim(fname)//'_'//trim(post)
  endif
  fname = trim(fname)//'.'//trim(fmt)

end subroutine makefnam

!===============================================================================

integer function timetohhmmss(time)

  use consts_coms, only: r8

  implicit none

  real(r8), intent(in) :: time

  integer  :: ihr, imn, isc
  real(r8) :: t

! Converts # of seconds since begining of day to hhmmss format

  t = time + .0001_r8    ! Adding small bias
  ihr = int(t/3600._r8)

  t = t - ihr*3600._r8
  imn = int(t/60._r8)

  t = t - imn*60._r8
  isc = int(t)

  timetohhmmss = ihr*10000 + imn*100 + isc
end function timetohhmmss

!===============================================================================

subroutine dintsort28(ni,chnums,cstr,array8)

use consts_coms, only: r8

implicit none

integer,           intent(in)    :: ni
character(len=14), intent(inout) :: chnums(ni)
character(len=*),  intent(inout) :: cstr  (ni)
real(r8),          intent(inout) :: array8(ni)

! sort an array of character strings by an associated character field

character(len=200) :: cscr
character(len=14)  :: mini,nscr
real(r8)           :: aa

integer :: n
integer :: nm
integer :: nmm

do n = 1,ni
   mini = '99999999999999'
   do nm = n,ni

      if (chnums(nm) < mini) then
         nmm = nm
         mini = chnums(nm)
      endif
   enddo

   nscr        = chnums(n)
   chnums(n)   = chnums(nmm)
   chnums(nmm) = nscr

   cscr      = cstr(n)
   cstr(n)   = cstr(nmm)
   cstr(nmm) = cscr

   aa          = array8(n)
   array8(n)   = array8(nmm)
   array8(nmm) = aa

enddo

return
end subroutine dintsort28

!===============================================================================

subroutine update_model_time(ctime,dtlong8)

  use misc_coms,   only: simtime
  use consts_coms, only: r8

  implicit none

  type(simtime), intent(inout) :: ctime
  real(r8),      intent(in)    :: dtlong8
  logical                      :: isleap

! Update the year, month, day, and time structure based on the current timestep

  isleap = ( (mod(ctime%year,400) == 0) .or. &
             (mod(ctime%year,100) /= 0 .and. mod(ctime%year,4) == 0) )

  ctime%time = ctime%time + dtlong8

  if (ctime%time + .0001_r8 >= 86400.0_r8) then  ! Using a small bias of .0001_r8

     ctime%time = ctime%time - 86400.0_r8
     ctime%date = ctime%date + 1
     if (ctime%date == 29 .and. ctime%month == 2 .and. .not. isleap) then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif (ctime%date == 30 .and. ctime%month == 2) then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif (ctime%date == 31) then
        if (ctime%month == 4 .or. ctime%month == 6 .or. ctime%month == 9  &
             .or. ctime%month == 11) then
           ctime%date = 1
           ctime%month = ctime%month + 1
        endif
     elseif (ctime%date == 32) then
        ctime%date = 1
        ctime%month = ctime%month + 1
        if(ctime%month == 13)then
           ctime%month = 1
           ctime%year = ctime%year + 1
        endif
     endif

  elseif (ctime%time + 0.0001_r8 < 0.0_r8) then  ! Using a small bias of .0001_r8

     ctime%time = ctime%time + 86400.0_r8
     ctime%date = ctime%date - 1
     if (ctime%date == 0) then
        ctime%month = ctime%month - 1
        if (ctime%month == 0) then
           ctime%month = 12
           ctime%year = ctime%year - 1
           ctime%date = 31
        else
           if (ctime%month == 2) then
              if (isleap) then
                 ctime%date = 29
              else
                 ctime%date = 28
              endif
           elseif (ctime%month == 4 .or. ctime%month == 6 .or.   &
                   ctime%month == 9 .or. ctime%month == 11) then
              ctime%date = 30
           else
              ctime%date = 31
           endif
        endif
     endif

  endif

end subroutine update_model_time

!===============================================================================

real function walltime(wstart)

  implicit none

  real, intent(in) :: wstart

  integer       :: ii
  integer, save :: previous = 0
  integer       :: imax
  real          :: rate
  real,    save :: adjustOverflow = 0.0

  call system_clock(count=ii, count_rate=rate, count_max=imax)

  if (ii < previous) then
    adjustOverflow = adjustOverflow + real(imax)/rate
  end if

  previous = ii

  walltime = adjustOverflow + real(ii)/rate - wstart

end function walltime

!===============================================================================

integer function day_of_week(m, d, y)
  implicit none

  integer, intent(in) :: d, m, y
  integer             :: j, k, mm, yy

  ! Given the month, day, and year, this function computes the
  ! day of week (Sun=1, Mon=2, .., Sat=7).

  mm = m
  yy = y
  if (mm <= 2) then
     mm = mm + 12
     yy = yy -  1
  end if
  j = yy / 100
  k = mod(yy, 100)
  day_of_week = mod(d + ((mm+1)*26)/10 + k + k/4 + j/4 + 5*j, 7)

end function day_of_week
