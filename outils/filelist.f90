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
subroutine OLAM_filelist(fnames,nfnam,file_prefix,nfile)

use mem_para,  only: myrank
use misc_coms, only: io6, tmpdir
use max_dims,  only: pathlen
implicit none

integer, intent(in)  :: nfnam
integer, intent(out) :: nfile

character(len=*), intent(in)  :: file_prefix
character(len=*), intent(out) :: fnames(nfnam)

character(len=256) :: command
character(pathlen) :: file, cdir, dirname
character(pathlen) :: tmpname, fname, dirtail
character(len=10)  :: zrank

integer :: iprelen,nc,nf,iun,lndir,slen,istat
logical :: exists
      
write(io6,*) 'OLAM_filelist: Checking prefix: '//trim(file_prefix)

nfile = 0
iprelen = len_trim(file_prefix)

lndir = scan(file_prefix, '/', back=.true.)
if (lndir > 0) then
   dirname = file_prefix(1:lndir-1)
   fname =   file_prefix(lndir+1:iprelen)
   slen =    scan(dirname, '/', back=.true.) + 1
   dirtail = dirname(slen:)
else
   dirname = "."
   dirtail = "."
   fname = file_prefix(1:iprelen)
endif

! Determine a unique temporary filename

write(zrank,'(I0)') myrank
tmpname = trim(tmpdir) // '/olam' // trim(zrank) // '.XXXXXXXX'
call mktempname(tmpname)
 
! The 'ls' command only takes a limited number of file names as
! arguments. Instead, we use the 'find' command to locate files and
! the '-prune' option to not recursively decend into sub-directories

command = 'find '//trim(dirname)//'/. "(" -type d -a ! -name . -prune ")" -o -name "'//trim(fname)//'" -print > '//trim(tmpname)
call system(command)

! Open the directory list and read through the files
 
inquire(file=tmpname, exist=exists)
if (.not. exists) then
   write(*,*) 'OLAM_filelist: Error opening temporary OLAM_filelist'
   stop 'OLAM_filelist: /tmp file error. Run again.'
endif

iun=98
open(unit=iun,file=tmpname,status='old')

do nf=1,nfnam+1
   read(iun,'(a)', iostat=istat) file

   ! iostat /= 0 indicates end-of-file or error so exit loop
   if (istat /= 0) exit

   if (nf.gt.nfnam) then
      write(*,*) 'OLAM_filelist: too many files of the form:'
      write(*,*) trim(file_prefix)
      stop 'Please increase the appropriate parameter in max_dims'
   endif

   fnames(nf) = file
enddo
      
close(iun)
nfile=nf-1

command = '/bin/rm -f '//tmpname
call system(command)

return 
end subroutine OLAM_filelist


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
