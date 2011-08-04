!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine date_abs_secs(cdate,seconds)

implicit none

character(len=14), intent(in) :: cdate
real(kind=8), intent(out) :: seconds

! compute number of seconds past 1 January 1900 12:00 am

real(kind=8) :: s1,s2,s3,s4
integer :: year1,month1,date1,hour1,iy,ndays
integer, external :: julday

call date_unmake_big(year1,month1,date1,hour1,cdate)

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1= dble(ndays) *86400.
s2= dble(hour1/10000)*3600.
s3= dble(mod(hour1,10000)/100)*60.
s4= dble(mod(hour1,100))
seconds= s1+s2+s3+s4

return
end subroutine date_abs_secs

!===============================================================================

subroutine date_abs_secs2(year1,month1,date1,hour1,seconds)

implicit none

integer, intent(in) :: year1,month1,date1,hour1
real(kind=8), intent(out) :: seconds

! compute number of seconds past 1 January 1900 12:00 am

integer :: iy,ndays
real(kind=8) :: s1,s2,s3,s4
integer, external :: julday

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1 = dble(ndays) * 86400.
s2 = dble(hour1/10000) * 3600.
s3 = dble(mod(hour1,10000)/100) * 60.
s4 = dble(mod(hour1,100))
seconds = s1 + s2 + s3 + s4

return
end subroutine date_abs_secs2

!===============================================================================

subroutine date_subtract(cdate1,cdate2,tinc,tunits)

implicit none

character(len=14), intent(in) :: cdate1, cdate2
real, intent(out) :: tinc
character(len=1), intent(in) :: tunits

! add (or subracts) a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

integer :: mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/
integer :: year1,month1,date1,hour1,year2,month2,date2,hour2
real(kind=8) :: secs1,secs2
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

implicit none

real(kind=8), intent(in) :: seconds
integer, intent(out) :: iyear1,imonth1,idate1,ihour1

! compute real time given number of seconds past 1 January 1900 12:00 am  

real(kind=8) :: s1
integer :: ny,nyr,ileap,nm,nd,ihr,imn,isc
integer :: mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/

! Get what year it is

s1=seconds
do ny=0,10000
   ileap=0
   if(mod(1900+ny,4) == 0) ileap=1
   s1=s1-(365.+ileap)*86400.
   if(s1 < 0.) then
      nyr=ny
      s1=s1+(365.+ileap)*86400.
      exit
   endif
enddo
iyear1=1900+nyr

! s1 is now number of secs into the year
!   Get month

do nm=1,12
   ileap=0
   if(mod(1900+ny,4) == 0 .and. nm == 2) ileap=1
   s1=s1-(mondays(nm)+ileap)*86400.
   if(s1 < 0.) then
      s1=s1+(mondays(nm)+ileap)*86400.
      exit
   endif
enddo
imonth1=nm

! s1 is now number of secs into the month
!   Get date and time

idate1=int(s1/86400.)
s1=s1-idate1*86400.
idate1=idate1+1 ! Since date starts at 1

ihr=int(s1/3600.)
s1=s1-ihr*3600.
imn=int(s1/60.)
s1=s1-imn*60.
isc=s1
ihour1=ihr*10000+imn*100+isc

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

integer :: imonth,iday,iyear

! compute the julian day from a normal date

julday= iday  &
      + min(1,max(0,imonth-1))*31  &
      + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4))))  &
      + min(1,max(0,imonth-3))*31  &
      + min(1,max(0,imonth-4))*30  &
      + min(1,max(0,imonth-5))*31  &
      + min(1,max(0,imonth-6))*30  &
      + min(1,max(0,imonth-7))*31  &
      + min(1,max(0,imonth-8))*31  &
      + min(1,max(0,imonth-9))*30  &
      + min(1,max(0,imonth-10))*31  &
      + min(1,max(0,imonth-11))*30  &
      + min(1,max(0,imonth-12))*31

return
end function julday

!===============================================================================

subroutine date_add_to_big8(cindate,tinc8,tunits,coutdate)

implicit none

character(len=14), intent(in) :: cindate
real(kind=8), intent(in) :: tinc8
character(len=1), intent(in) :: tunits
character(len=14), intent(out) :: coutdate

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(kind=8) :: ttinc,secs
integer :: inyear,inmonth,indate,inhour  &
          ,outyear,outmonth,outdate,outhour

! convert input time to seconds

ttinc = tinc8
if (tunits == 'm') ttinc = tinc8 * 60.
if (tunits == 'h') ttinc = tinc8 * 3600.
if (tunits == 'd') ttinc = tinc8 * 86400.

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
implicit none

integer, intent(in) :: inyear,inmonth,indate,inhour
real(kind=8), intent(in) :: tinc8
character(len=1), intent(in) :: tunits
integer, intent(out) :: outyear,outmonth,outdate,outhour

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(kind=8) :: ttinc,secs

! convert input time to seconds

ttinc = tinc8
if (tunits == 'm') ttinc = tinc8 * 60.
if (tunits == 'h') ttinc = tinc8 * 3600.
if (tunits == 'd') ttinc = tinc8 * 86400.

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
  implicit none
  real, intent(in) :: time

  integer          :: ihr, imn, isc
  real             :: t

! Converts # of seconds since begining of day to hhmmss format

  t = time
  ihr = int(t/3600.)

  t = t - ihr*3600.
  imn = int(t/60.)

  t = t - imn*60.
  isc = int(t)

  timetohhmmss = ihr*10000+imn*100+isc
end function timetohhmmss

!===============================================================================

subroutine dintsort28(ni,chnums,cstr,array8)

implicit none

integer,           intent(in)    :: ni
character(len=14), intent(inout) :: chnums(ni)
character(len=*),  intent(inout) :: cstr  (ni)
real(kind=8),      intent(inout) :: array8(ni)

! sort an array of character strings by an associated character field

character(len=200) :: cscr
character(len=14)  :: mini,nscr
real(kind=8)       :: aa

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

subroutine update_model_time(ctime,dtlong)
  use misc_coms, only: simtime
  implicit none

  type(simtime), intent(inout) :: ctime
  real,          intent(in)    :: dtlong
  logical                      :: isleap

! Update the year, month, day, and time structure based on the current timestep
  
  isleap = ( (mod(ctime%year,400) == 0) .or. &
             (mod(ctime%year,100) /= 0 .and. mod(ctime%year,4) == 0) )
  
  ctime%time = ctime%time + dtlong

  if (ctime%time.ge.86400.0) then

     ctime%time = ctime%time - 86400.0
     ctime%date = ctime%date + 1
     if (ctime%date.eq.29 .and. ctime%month.eq.2 .and. .not. isleap) then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif (ctime%date.eq.30 .and. ctime%month.eq.2) then
        ctime%date = 1
        ctime%month = ctime%month + 1
     elseif (ctime%date.eq.31) then
        if (ctime%month.eq.4 .or. ctime%month.eq.6 .or. ctime%month.eq.9  &
             .or. ctime%month.eq.11) then
           ctime%date = 1
           ctime%month = ctime%month + 1
        endif
     elseif (ctime%date.eq.32) then
        ctime%date = 1
        ctime%month = ctime%month + 1
        if(ctime%month.eq.13)then
           ctime%month = 1
           ctime%year = ctime%year + 1
        endif
     endif

  elseif (ctime%time.lt.0.0) then

     ctime%time = ctime%time + 86400.0
     ctime%date = ctime%date - 1
     if (ctime%date.eq.0) then
        ctime%month = ctime%month - 1
        if (ctime%month.eq.0) then
           ctime%month = 12
           ctime%year = ctime%year - 1
           ctime%date = 31
        else
           if (ctime%month.eq.2) then
              if (isleap) then
                 ctime%date = 29
              else
                 ctime%date = 28
              endif
           elseif (ctime%month.eq.4 .or. ctime%month.eq.6 .or.   &
                   ctime%month.eq.9 .or. ctime%month.eq.11) then
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
  integer :: ii,ir,im
  integer, save :: previous = 0
  integer, save :: adjustOverflow = 0.0

  call system_clock(count=ii,count_rate=ir,count_max=im)

  if (ii < previous) then
    adjustOverflow = adjustOverflow + real(im)/real(ir)
  end if

  previous = ii

  walltime=adjustOverflow + real(ii)/real(ir) - wstart

end function walltime
