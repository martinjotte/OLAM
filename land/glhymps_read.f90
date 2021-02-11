!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

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

!===============================================================================

subroutine glhymps_read()

  use mem_land,    only: land, nland, onland
  use mem_sfcg,    only: sfcg
  use consts_coms, only: erad, piu180
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use misc_coms,   only: io6
  use max_dims,    only: pathlen
  use leaf_coms,   only: glhymps_database

  implicit none

  real, allocatable :: dato(:,:)

  integer :: nio, njo
  integer :: io, jo
  integer :: idataset, iland, iwsfc
  integer :: ndims, idims(2)

  real :: xoffpix, yoffpix
  real :: glat, glon
  real :: xperdeg, yperdeg

  character(pathlen) :: fname, fullname

  logical :: exists

  ! This subroutine is simpler than land_database_read because it assumes that 
  ! each glhymps database file covers the entire geographic area of the model.
  ! If this ever changes, this subroutine must be modified.

  ! Loop over the three glhymps datasets

  do idataset = 1,3

     if     (idataset == 1) then
        fname = 'permeability_no_permafrost_5km.h5'
     elseif (idataset == 2) then
        fname = 'permeability_permafrost_5km.h5'
     elseif (idataset == 3) then
        fname = 'porosity_5km.h5'
     endif

     ! Make full filename

     fullname = trim(glhymps_database)//trim(fname)

     ! Open the sst filelist and read through the files

     print*, 'glhymps_read reading ',trim(fullname)

     INQUIRE(file=fullname, exist=exists)
     IF (.not. exists) THEN
        WRITE(*,*) 'glhymps: Error opening glhymps file ' // trim(fullname)
        STOP
     ENDIF

     ! Open and read glhymps_database file

     call shdf5_open(trim(fullname),'R')
     call shdf5_info('Band1',ndims,idims)

     nio = idims(1)
     njo = idims(2)
     allocate(dato(nio,njo))

     call shdf5_irec(ndims,idims,'Band1',rvar2=dato)
     call shdf5_close()

     ! glhymps data files do NOT have an added perimeter, and data values are
     ! offset by 1/2 pixel from integer lat-lon values.

     xoffpix = 0.5
     yoffpix = 0.5

     xperdeg = real(nio) / 360.0
     yperdeg = real(njo) / 180.0

     ! Fill olam-soil array from nearest neighbor in dataset field

     do iland = 2, nland
        iwsfc = iland + onland
   
        glat = sfcg%glatw(iwsfc)
        glon = sfcg%glonw(iwsfc)

        glon = max(-179.999,min(179.999,glon))

        io = nint((glon + 180.) * xperdeg + xoffpix)
        jo = nint((glat +  90.) * yperdeg + yoffpix)

        ! Exclude "undefined" and "missing data values.
        ! For idataset = 1 and 2, dato = log10(permeability). Convert this to ksat; 
        ! ksat is 1.e7 * permeability

        if (iwsfc == 13328) then
           print*, 'glhymps_read1 ',idataset,iland,io,jo,dato(io,jo)
        endif

        if     (idataset == 1) then
           if (dato(io,jo) >= -20. .and. dato(io,jo) < -10.) then
              land%glhymps_ksat(iland) = 10.0**(dato(io,jo) + 7.0)
           else
              land%glhymps_ksat(iland) = 10.0**(-20.0 + 7.0) ! Require a nonzero minimum value
           endif
        elseif (idataset == 2) then
           if (dato(io,jo) >= -20. .and. dato(io,jo) < -10.) then
              land%glhymps_ksat_pfr(iland) = 10.0**(dato(io,jo) + 7.0)
           else
              land%glhymps_ksat_pfr(iland) = 10.0**(-20.0 + 7.0) ! Require a nonzero minimum value
           endif
        elseif (idataset == 3) then
           if (dato(io,jo) >= 0.01 .and. dato(io,jo) < 0.5) then
              land%glhymps_poros(iland) = dato(io,jo)
           else
              land%glhymps_poros(iland) = 0.01 ! Require a nonzero minimum value
           endif
        endif
     enddo

     deallocate(dato)

  enddo ! idataset

end subroutine glhymps_read
