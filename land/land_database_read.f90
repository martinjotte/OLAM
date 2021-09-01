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
Module land_db

Contains

  subroutine land_database_read(nqa,glatq,glonq,ofn,ofn2,iaction,idatq,datq)

  ! The letter "q" represents any point in the grid stagger, as determined by
  ! the routine that calls this subroutine.

  use consts_coms, only: erad, piu180
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec
  use misc_coms,   only: io6
  use max_dims,    only: pathlen

  implicit none

  integer, intent(in) :: nqa
  real,    intent(in) :: glatq(nqa), glonq(nqa)

  character(*), intent(in) :: ofn, ofn2
  character(*), intent(in) :: iaction

  integer, optional, intent(out) :: idatq(nqa)
  real,    optional, intent(out) ::  datq(nqa)

  integer :: qtable(nqa)

  integer :: nio, njo
  integer :: niosh, njosh
  integer :: isoc, iwoc
  integer :: isocpt, isocpo
  integer :: iwocph, iwocpt, iwocpo
  integer :: io, jo
  integer :: ind
  integer :: ind1, ind2
  integer :: io1, io2, jo1, jo2
  integer :: ifiles, jfiles
  integer :: ifile, jfile
  integer :: io_full, jo_full
  integer :: iq
  integer :: nperdeg
  integer :: ndims, idims(2)

  real :: offpix 
  real :: qlat1, qlon1
  real :: wio1, wio2, wjo1, wjo2
  real :: rio_full, rjo_full

  integer(kind=4), allocatable :: idato (:,:)
  real(kind=4), allocatable :: dato (:,:)

  integer, allocatable :: numq    (:,:)
  integer, allocatable :: numqind1(:,:)
  integer, allocatable :: numqind2(:,:) ! (ifiles,jfiles)

  character(len=3)   :: title1
  character(len=4)   :: title2
  character(pathlen) :: fname

  logical :: l1,l2

  ! Open, read, and close dataset header file

  fname = trim(ofn)//'HEADER'
  inquire(file=fname, exist=l1)

  if (.not. l1) then
     write(io6,*)
     write(io6,*) '==================================================='
     write(io6,*) '| Problem in land_database_read:'
     write(io6,*) '| Header file ', trim(fname)
     write(io6,*) '| not found!'
     write(io6,*) '==================================================='
     stop 'land_database_read'
  endif

  open(29, file=fname, form='FORMATTED', status='OLD', action='READ')
  read(29,*) nio, njo, nperdeg
  close(29)

  ! Compute number of pixels in a shift to adjacent file (niosh and njosh).
  ! Compute number of files in database that span all latitudes and
  ! longitudes on earth. [This algorithm will change when multiple resolutions
  ! of the SRTM data become available.]

  if (mod(nio,nperdeg) == 2) then
     offpix = .5
     niosh = nio - 2
     njosh = njo - 2
  else
     offpix = 0.
     niosh = nio - 1
     njosh = njo - 1
  endif

  ifiles = 360 * nperdeg / niosh
  jfiles = 180 * nperdeg / njosh

  ! Allocate 5 arrays.

  allocate (numq    (ifiles,jfiles))
  allocate (numqind1(ifiles,jfiles))
  allocate (numqind2(ifiles,jfiles))
  allocate (idato   (nio,njo))
  allocate (dato    (nio,njo))

  do jfile = 1,jfiles
     do ifile = 1,ifiles
        numq(ifile,jfile) = 0
     enddo
  enddo

  ! Loop over all geographic points that need data filled, determine which
  ! file to read for each one, and count number of points to be read from
  ! each file.

  do iq = 2,nqa

     qlat1 = max(-89.9999,min(89.9999,glatq(iq)))
     qlon1 = glonq(iq)

     if (qlon1 >=  180.) qlon1 = qlon1 - 360.
     if (qlon1 <  -180.) qlon1 = qlon1 + 360.
     qlon1 = max(-179.9999,min(179.9999,qlon1))

     rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
     rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

     io_full = int(rio_full)
     jo_full = int(rjo_full)

     ! If io_full and/or jo_full are at max value for database, decrease by 1 so
     ! that out-of-range file will not be sought.

     if (io_full == 360 * nperdeg) io_full = io_full - 1
     if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

     ifile = io_full / niosh + 1
     jfile = jo_full / njosh + 1

     numq(ifile,jfile) = numq(ifile,jfile) + 1

  enddo

  ! Set up array index values for qtable array

  ind = 1
  do jfile = 1,jfiles
     do ifile = 1,ifiles
        numqind1(ifile,jfile) = ind
        numqind2(ifile,jfile) = ind
        ind = ind + numq(ifile,jfile)
     enddo
  enddo

  ! Fill qtable array

  do iq = 2,nqa

     qlat1 = max(-89.9999,min(89.9999,glatq(iq)))
     qlon1 = glonq(iq)

     if (qlon1 >=  180.) qlon1 = qlon1 - 360.
     if (qlon1 <  -180.) qlon1 = qlon1 + 360.
     qlon1 = max(-179.9999,min(179.9999,qlon1))

     rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
     rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

     io_full = int(rio_full)
     jo_full = int(rjo_full)

     ! If io_full and/or jo_full are at max value for database, decrease by 1 so
     ! that out-of-range file will not be sought.

     if (io_full == 360 * nperdeg) io_full = io_full - 1
     if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

     ifile = io_full / niosh + 1
     jfile = jo_full / njosh + 1

     ind = numqind2(ifile,jfile)

     qtable(ind) = iq
     numqind2(ifile,jfile) = numqind2(ifile,jfile) + 1 ! Sums to ending index

  enddo

  ! Read files and extract data

  do jfile = 1,jfiles
     do ifile = 1,ifiles

        ind1 = numqind1(ifile,jfile)
        ind2 = numqind2(ifile,jfile)

        if (ind2 > ind1) then
           iwoc = (ifile - 1) * niosh / nperdeg - 180  ! SW longitude of current file
           isoc = (jfile - 1) * njosh / nperdeg -  90  ! SW latitude of current file

           ! Construct filename

           isocpt = abs(isoc) / 10
           isocpo = abs(isoc) - isocpt*10
           iwocph = abs(iwoc) / 100
           iwocpt = (abs(iwoc) - iwocph * 100) / 10
           iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10

           if (isoc >= 0) then
              write(title1,'(2i1,a1)') isocpt,isocpo,'N'
           else
              write(title1,'(2i1,a1)') isocpt,isocpo,'S'
           endif

           if (iwoc >= 0) then
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
           else
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
           endif

           fname = trim(ofn2)//title1//title2//'.h5'

!D Modification for reading Amazon deforestation OGE files

!D         if (trim(iaction) == 'leaf_class') then
!D            if (title1//title2 == '00N060W' .or. &
!D                title1//title2 == '00N090W' .or. &
!D                title1//title2 == '30S060W' .or. &
!D                title1//title2 == '30S090W') then

!D  fname = '../../olamdatah5/BAU_SoaresFilho/2050/OGE2_'//title1//title2//'.h5'

!D            endif
!D         endif

           inquire(file=fname,exist=l1,opened=l2)

           ! Read file

           if (l1) then
              write(io6,*) 'getting file ',trim(fname)

              call shdf5_open(fname,'R')

              ndims = 2
              idims(1) = nio
              idims(2) = njo

              if     (trim(iaction) == 'topo') then
                 call shdf5_irec(ndims,idims,'topo',rvar2=dato)
              elseif (trim(iaction) == 'leaf_class') then
                 call shdf5_irec(ndims,idims,'oge2',ivar2=idato)
              elseif (trim(iaction) == 'soil_text') then
                 call shdf5_irec(ndims,idims,'fao',ivar2=idato)
              elseif (trim(iaction) == 'ndvi') then
                 call shdf5_irec(ndims,idims,'ndvi',rvar2=dato)
              elseif (trim(iaction) == 'wtd') then
                 call shdf5_irec(ndims,idims,'WTD',rvar2=dato)
              elseif (iaction == 'orog') then
                 call shdf5_irec(ndims,idims,'orog',rvar2=dato)
              elseif (iaction == 'etopo1') then
                 call shdf5_irec(ndims,idims,'etopo1',rvar2=dato)
              else
                 write(io6,*) 'incorrect action specified in leaf_database'
                 write(io6,*) 'stopping run'
                 stop 'stop landuse_input1'
              endif

              call shdf5_close()
           else
              write(io6,*) 'In landuse_input, ',iaction,' file is missing'
              write(io6,*) 'Filename = ',trim(fname)
              write(io6,*) 'Stopping model run'
              stop 'stop_landuse_input2'
           endif

           !$omp parallel do private(iq,qlat1,qlon1,rio_full,rjo_full,io_full, &
           !$omp                     jo_full,io1,jo1,wio2,wjo2,io2,jo2,wio1,wjo1,io,jo)
           do ind = ind1,ind2-1
              iq = qtable(ind)

              qlat1 = max(-89.9999,min(89.9999,glatq(iq)))
              qlon1 = glonq(iq)

              if (qlon1 >=  180.) qlon1 = qlon1 - 360.
              if (qlon1 <= -180.) qlon1 = qlon1 + 360.

              rio_full = (qlon1 + 180.) * nperdeg ! must ignore pixel offset here
              rjo_full = (qlat1 +  90.) * nperdeg ! must ignore pixel offset here

              io_full = int(rio_full)
              jo_full = int(rjo_full)

              ! If io_full and/or jo_full are at max value for database,
              ! decrease by 1 so that out-of-range file will not be sought.

              if (io_full == 360 * nperdeg) io_full = io_full - 1
              if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

              io1 = mod(io_full,niosh) + 1
              jo1 = mod(jo_full,njosh) + 1

              wio2 = rio_full - float(io_full) + offpix
              wjo2 = rjo_full - float(jo_full) + offpix

              ! At this point, io1, jo1, wio2, and wjo2 are correct if offpix = 0,
              ! but need correction if offpix = .5

              if (wio2 > 1.) then
                 wio2 = wio2 - 1.
                 io1 = io1 + 1
              endif

              if (wjo2 > 1.) then
                 wjo2 = wjo2 - 1.
                 jo1 = jo1 + 1
              endif

              ! Correction for offpix = .5 is now completed

              io2 = io1 + 1
              jo2 = jo1 + 1

              wio1 = 1. - wio2
              wjo1 = 1. - wjo2

              if (iaction == 'topo' .or. &
                  iaction == 'ndvi' .or. &
                  iaction == 'orog' .or. &
                  iaction == 'etopo1') then

                 ! Interpolate from 4 surrounding values

                 datq(iq) &
                    = wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2)) &
                    + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

              elseif (trim(iaction) == 'wtd') then

                 ! Interpolate from 4 surrounding values unless any are "missing" (negative)

                 if (min(dato(io1,jo1), dato(io1,jo2), &
                         dato(io2,jo1), dato(io2,jo2)) > -1.e-6) then

                    datq(iq) &
                       = wio1 * (wjo1 * dato(io1,jo1) + wjo2 * dato(io1,jo2)) &
                       + wio2 * (wjo1 * dato(io2,jo1) + wjo2 * dato(io2,jo2))

                 else

                    ! If any values are missing, set datp to maximum of surrounding values

                    datq(iq) = max(dato(io1,jo1), dato(io1,jo2), &
                                   dato(io2,jo1), dato(io2,jo2))

                 endif

              elseif (trim(iaction) == 'leaf_class' .or. &
                      trim(iaction) == 'soil_text') then

                 ! Use nearest data point - do not interpolate

                 io = io2
                 jo = jo2
                 if (wio2 < .5) io = io1
                 if (wjo2 < .5) jo = jo1

                 idatq(iq) = idato(io,jo)

              endif    ! iaction

           enddo    ! ind
           !$omp end parallel do

       endif     ! ind2 > ind1
     enddo    ! ifile
  enddo    ! jfile

  deallocate(numq,numqind1,numqind2,idato,dato)

  end subroutine land_database_read

End module land_db
