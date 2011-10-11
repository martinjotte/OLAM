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
subroutine isan_file_inv()

use max_dims,  only: maxisfiles, pathlen
use misc_coms, only: io6, iyear1, imonth1, idate1, itime1
use isan_coms, only: iapr, isdirs, fnames_fg, s1900_fg, ctotdate_fg,  &
                     nfgfiles

implicit none

integer :: nc,nf,lnf,nn,ndates,isan_err_flag
integer :: nfg_tmp,n
integer :: iyear,imonth,idate,ihour

character(pathlen) :: fnames_tmp(maxisfiles)
character(15) :: suffix
logical :: nocall

! Go through first guess files and make inventory

nfgfiles = 0
fnames_fg = ''

! new schema to avoid call systems
IF (isdirs == -1) THEN
   isdirs = 1
   nocall = .TRUE.
ELSE
   suffix = "????-??-??-????"
   nocall = .FALSE.
ENDIF

do n = 1,isdirs 
   
   if (len_trim(iapr(n)) == 0) cycle
   if (iapr(n)(1:1) == char(0)) cycle

   iapr(n) = trim(adjustl(iapr(n)))

   if (n >= 2) then
      if (any(iapr(n) == iapr(1:n-1))) cycle
   endif

   nfg_tmp=-1
   IF (nocall) THEN
      CALL OLAM_filelist(fnames_tmp, maxisfiles, &
           TRIM(iapr(n)), nfg_tmp, nocall)
   ELSE
      CALL OLAM_filelist(fnames_tmp, maxisfiles, &
           TRIM(iapr(n))//TRIM(suffix), nfg_tmp, nocall)
   ENDIF
   nf = nfgfiles + 1
   nfgfiles = nfgfiles + nfg_tmp
   
   if (nfgfiles > maxisfiles) then
      write(io6,*) 'Too many first guess files'
      write(io6,*) 'Increase maxisfiles in maxdims.f90 if you need more first' 
      write(io6,*) 'guess files' 
      stop         'maxisfiles exceeded'
   endif

   fnames_fg(nf:nfgfiles) = fnames_tmp(1:nfg_tmp)

enddo

do nf = 1,nfgfiles

   lnf=len_trim(fnames_fg(nf))

! assume files have the form of dp-p2005-09-01-0000

   read (fnames_fg(nf)(lnf-14:lnf), '(i4,1x,i2,1x,i2,1x,i4)' ) &
        iyear, imonth, idate, ihour

! assume files have the form of dp-p2005-09-01-000000-g1.h5

!!   read (fnames_fg(nf)(lnf-22:lnf-8), '(i4,1x,i2,1x,i2,1x,i4)' ) &
!!        iyear, imonth, idate, ihour

   call date_make_big (iyear,imonth,idate,ihour*100,ctotdate_fg(nf))
   call date_abs_secs2(iyear,imonth,idate,ihour*100,s1900_fg(nf))

enddo

call dintsort28(nfgfiles,ctotdate_fg,fnames_fg,s1900_fg)

do nf = 1,nfgfiles

   write(io6,*) ' '
   write(io6,*) 'fgfiles0 ',nfgfiles,nf,ctotdate_fg(nf),s1900_fg(nf)
   write(io6,*) 'fgfiles1 ',fnames_fg(nf)
   
enddo

return
end subroutine isan_file_inv
