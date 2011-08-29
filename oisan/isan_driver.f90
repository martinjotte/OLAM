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
subroutine isan_driver(iaction)

use isan_coms,  only: innpr, ihour, idate, imonth, iyear, nfgfiles, ifgfile, &
                      ctotdate_fg, fnames_fg, s1900_fg

use misc_coms,  only: io6, initial, runtype, s1900_init, s1900_sim

use mem_zonavg, only: zonavg_init

implicit none

integer, intent(in) :: iaction

character(len=3) :: fform

integer :: i, k, ifileok, nf

! Check type of call to isan_driver

if (iaction == 0) then

! Model run is being started, either for initialization or history restart

! Do inventory of isan file names and times

   call ISAN_file_inv()

! Loop over isan files and search for the one that corresponds to current
! or most recent model time.

   ifgfile = 0
   do nf = 1,nfgfiles
      if (s1900_fg(nf) <= s1900_sim) then
         ifgfile = nf
      endif
   enddo

   if (ifgfile < 1) then
      write(io6,*) ' '
      write(io6,*) 'Unable to find isan file for current model time '
      write(io6,*) 'Stopping model '
      stop 'stop: no current isan file'
   elseif (runtype == 'INITIAL' .and. &
           s1900_fg(ifgfile) < s1900_init - 1800.d0) then
      write(io6,*) ' '
      write(io6,*) 'Available isan file is more than 1800 seconds '
      write(io6,*) 'prior to simulation initialization time. '
      write(io6,*) 'Stopping model '
      stop 'stop: initial isan file'
   endif

! Return if not initializing model but only setting ifgfile value

   if (runtype /= 'INITIAL') return

elseif (iaction == 1) then

! Processing next isan file (only called with iaction = 1 if nudflag = 1)

   ifgfile = ifgfile + 1

   if (ifgfile > nfgfiles) then
      write(io6,*) ' '
      write(io6,*) 'No future isan file is available for nudging '
      write(io6,*) 'Stopping model '
      stop 'stop: no future isan file'
   endif

endif

! Process current isan file
        
call date_unmake_big(iyear,imonth,idate,ihour,ctotdate_fg(ifgfile))
ihour = ihour / 100

write(io6,*)
write(io6,*) 'Reading ISAN file ifgfile = ',ifgfile
write(io6,*) fnames_fg(ifgfile)
write(io6,'(a,4i6)') ctotdate_fg(ifgfile),iyear,imonth,idate,ihour

innpr = fnames_fg(ifgfile)

! Fill zonavg arrays for current time

call zonavg_init(idate,imonth,iyear)

! Read header information from gridded pressure-level files for this file time.
! This information includes input data array dimensions.

call read_press_header(fform)

! Process isan data for this file time.

call isan_singletime(iaction,fform)

return
end subroutine isan_driver

!=======================================================================

subroutine isan_singletime(iaction,fform)

use isan_coms, only: nprx, npry, nprz
use mem_grid,  only: mza, mwa, mva
use misc_coms, only: io6

implicit none

integer, intent(in) :: iaction
character(len=3), intent(in) :: fform

! Automatic arrays

real :: p_u(nprx+3,npry+2,nprz)
real :: p_v(nprx+3,npry+2,nprz)
real :: p_t(nprx+3,npry+2,nprz)
real :: p_z(nprx+3,npry+2,nprz)
real :: p_r(nprx+3,npry+2,nprz)

real :: p_slp (nprx+3,npry+2)
real :: p_sfp (nprx+3,npry+2)
real :: p_sft (nprx+3,npry+2)
real :: p_snow(nprx+3,npry+2)
real :: p_sst (nprx+3,npry+2)

real(kind=8) :: o_rho   (mza,mwa)
real         :: o_theta (mza,mwa)
real         :: o_shv   (mza,mwa)
real         :: o_umzonal(mza,mwa)
real         :: o_ummerid(mza,mwa)
real         :: o_uvc    (mza,mva) ! uc or vc wind component

! Fill surface arrays with 'missing' values in case they are not present

p_slp (1:nprx+3,1:npry+2) = -999.0
p_sfp (1:nprx+3,1:npry+2) = -999.0
p_sft (1:nprx+3,1:npry+2) = -999.0
p_snow(1:nprx+3,1:npry+2) = -999.0
p_sst (1:nprx+3,1:npry+2) = -999.0

! Read in gridded pressure-level data and copy to isan arrays

call pressure_stage(fform,p_u, p_v, p_t, p_z, p_r, &
                    p_slp, p_sfp, p_sft, p_snow, p_sst)

! Add pressure-level data at higher levels from climatology and interpolate
! combined data to OLAM grid

call isnstage(p_u,p_v,p_t,p_z,p_r, &
              o_rho, o_theta, o_shv, o_umzonal, o_ummerid, o_uvc)

! Fill model initial fields and/or nudging obs values

call fldsisan(iaction, o_rho, o_theta, o_shv, o_umzonal, o_ummerid, o_uvc)

return
end subroutine isan_singletime
