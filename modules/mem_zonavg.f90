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
Module mem_zonavg

real, allocatable :: zont12(:,:,:)
real, allocatable :: zonu12(:,:,:)
real, allocatable :: zont(:,:)
real, allocatable :: zonu(:,:)
real, allocatable :: zonr(:,:)
real, allocatable :: zonz(:,:)
real, allocatable :: zonp_vect(:)
real, allocatable :: zonr_vect(:)
real, allocatable :: zonz_vect(:)
real, allocatable :: zont_vect(:)
real, allocatable :: zonu_vect(:)

Contains

   subroutine alloc_zonavg()

   implicit none
   
   allocate(zont(74,22),zonu(74,22),zonr(74,22),zonz(74,22)  &
      ,zonp_vect(22),zonr_vect(22),zonz_vect(22),zont_vect(22),zonu_vect(22))

   return
   end subroutine alloc_zonavg

!===============================================================================

   subroutine dealloc_zonavg()

   implicit none
   
   deallocate(zont,zonu,zonz,zonr  &
      ,zonp_vect,zont_vect,zonu_vect,zonz_vect,zonr_vect)

   return
   end subroutine dealloc_zonavg

!===============================================================================

   subroutine fill_zonavg(zonclim,idate1,imonth1)

   use misc_coms, only: io6

   implicit none

   integer, intent(in) :: idate1,imonth1
   character(len=*), intent(in) :: zonclim

   integer :: ilat,iplev,im,im1,im2
   real :: wt1,wt2
   logical :: fexists
   character(len=1) :: line
   integer, external :: julday
   
   allocate(zont12(74,22,12),zonu12(74,22,12))
   
   inquire(file=zonclim, exist=fexists)
   if (.not. fexists) then
      write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(io6,*) '!!!   Trying to open zonally-averaged climate file:'
      write(io6,*) '!!!   '//trim(zonclim)
      write(io6,*) '!!!   but it does not exist. Please set the correct'
      write(io6,*) '!!!   path for ZONCLIM in the OLAMIN file.         '
      write(io6,*) '!!!   The run is ended.                            '
      write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop 'in fill_zonavg'
   endif

   open(29, file=zonclim, form='FORMATTED', status='OLD', action='READ')

   do im = 1,12  ! 'month' index
      read(29,*) line
      do iplev = 22,1,-1
         read(29,'(72f7.2)') (zont12(ilat+1,iplev,im),ilat=1,72)
      enddo
   enddo

   do im = 1,12
      read(29,*) line
      do iplev = 22,1,-1
         read(29,'(72f7.2)') (zonu12(ilat+1,iplev,im),ilat=1,72)
      enddo
   enddo

   close(29)

! Array time indices and interpolation weights for bracketing months

   im1 = imonth1
   if (idate1 < 16) im1 = mod(imonth1+10,12) + 1
   im2 = mod(im1,12) + 1
   wt2 = .033333 * mod(idate1+15,31)
   wt1 = 1. - wt2

! Interpolate zonal T and u in time, fill in N/S boundary values,
! and fill pressure column array

   do iplev = 1,22
      do ilat = 2,73
         zont(ilat,iplev) = wt1 * zont12(ilat,iplev,im1)  &
                          + wt2 * zont12(ilat,iplev,im2)
         zonu(ilat,iplev) = wt1 * zonu12(ilat,iplev,im1)  &
                          + wt2 * zonu12(ilat,iplev,im2)
      enddo
      zont(1 ,iplev) =  zont(2 ,iplev)
      zont(74,iplev) =  zont(73,iplev)
      zonu(1 ,iplev) = -zonu(2 ,iplev)
      zonu(74,iplev) = -zonu(73,iplev)
      zonp_vect(iplev) = 10. ** (float(31-iplev)/6.) 
   
      write(io6,*) 'zonp_vect1 ',iplev,zonp_vect(iplev)
   
   enddo

   deallocate(zont12,zonu12)

   return
   end subroutine fill_zonavg

End Module mem_zonavg


