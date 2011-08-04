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
subroutine topo_database_read(nqa,xeq,yeq,zeq,topq,arq0,glatq,glonq)

! The letter "q" represents any point in the grid stagger, as determined by
! the subroutine call statement.

use misc_coms,   only: io6, topo_database
use consts_coms, only: erad, piu180, pi1, pi2
use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec

implicit none

integer, intent(in) :: nqa

real, intent(in) :: xeq(nqa),yeq(nqa),zeq(nqa),arq0(nqa),glatq(nqa),glonq(nqa)

real, intent(out) :: topq(nqa)

integer, parameter :: nradp = 1
integer, parameter :: np = 3*nradp*(nradp-1)+1 ! # data values avgd per topm pt

integer, allocatable :: ptable(:)
integer, allocatable :: nump(:,:)
integer, allocatable :: numpind1(:,:)
integer, allocatable :: numpind2(:,:)

real, allocatable :: dato(:,:)
real, allocatable :: datp(:,:)
real, allocatable :: pqlat(:,:)
real, allocatable :: pqlon(:,:)

integer :: ndims,idims(2)
logical l1,l2

character(len=3)  :: title1
character(len=4)  :: title2
character(len=99) :: fname

integer :: iq,ip,jp,nio,njo,niosh,njosh,nperdeg,isoc,iwoc,isocpt,isocpo,  &
           iwocph,iwocpt,iwocpo,lb,ind,j2d,j1d,ind1,ind2,io1,jo1,io2,jo2, &
           ifiles,jfiles,ifile,jfile,ptab,ptab0,j,io_full,jo_full,iw,     &
           iradp,nazimp,iazimp

real :: offlat,offlon,pqlat1,pqlon1,dlato,dlono,wio2,wjo2,wio1,wjo1, &
        yp,xp,c1,c2,rio_full,rjo_full,radq,radp,azimpoff,azimp,offpix

write(io6,*) 'Starting topography database read'

! Check topography file prefix
   
lb = len_trim(topo_database)
if (lb <= 0) then
   write(io6,*)
   write(io6,*) '==================================================='
   write(io6,*) '|  Problem in topo_database: Input data prefix incorrect !'
   write(io6,*) '|  File prefix:',topo_database
   write(io6,*) '==================================================='
   stop 'topo_database_file'
endif

! Open, read, and close topography dataset header file

fname = trim(topo_database)//'HEADER'
inquire(file=fname, exist=l1)
if (.not. l1) then
   write(io6,*)
   write(io6,*) '==================================================='
   write(io6,*) '|  Problem in topo_database:'
   write(io6,*) '|  Header file ', trim(fname)
   write(io6,*) '|  not found!'
   write(io6,*) '==================================================='
   stop 'topo_database_file'
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
   
! Allocate 8 arrays

allocate (nump    (ifiles,jfiles)) ! # pts filled by file(ifile,jfile)
allocate (numpind1(ifiles,jfiles))
allocate (numpind2(ifiles,jfiles))
allocate (ptable  (np*nqa))

allocate (dato(nio,njo))
allocate (datp(np,nqa))
allocate (pqlat(np,nqa))
allocate (pqlon(np,nqa))

! Fill pqlat and pqlon arrays with offset latitudes and longitudes of all p points

do iq = 2,nqa

! Compute radius from current Q point for averaging topo data values

   radq = sqrt(arq0(iq) / pi1)

! Loop over array (array size options = 1,7,19,37) of pts within radius 
! RADQ - get lat/lon for each pt

   ip = 0

   do iradp = 1,nradp

      if (iradp == 1) then

         ip = ip + 1
         pqlat(ip,iq) = glatq(iq)
         pqlon(ip,iq) = glonq(iq)

      else

         radp = real(2 * iradp - 2) / real(2 * nradp - 1) * radq
         nazimp = 6 * (iradp - 1)
         azimpoff = .5 * mod(iradp,2)  ! offset to stagger azimp of alt. rows

         do iazimp = 1,nazimp
            azimp = (real(iazimp) + azimpoff) / real(nazimp) * pi2  
            xp = radp * cos(azimp)
            yp = radp * sin(azimp)
            ip = ip + 1
            call xy_ll(pqlat(ip,iq),pqlon(ip,iq),glatq(iq),glonq(iq),xp,yp)
         enddo

      endif
   enddo

enddo

do jfile = 1,jfiles
   do ifile = 1,ifiles
      nump(ifile,jfile) = 0
   enddo
enddo

! Get file index (ifile,jfile) within full dataset and count number of p 
! points (nump) that occur in each file

do iq = 2,nqa
   do ip = 1,np

      pqlat1 = max( -89.9999,min( 89.9999,pqlat(ip,iq)))
      pqlon1 = max(-179.9999,min(179.9999,pqlon(ip,iq)))

      rio_full = (pqlon1 + 180.) * nperdeg ! must ignore pixel offset here
      rjo_full = (pqlat1 +  90.) * nperdeg ! must ignore pixel offset here

      io_full = int(rio_full)
      jo_full = int(rjo_full)

! If io_full and/or jo_full are at max value for database, decrease by 1 so
! that out-of-range file will not be sought.

      if (io_full == 360 * nperdeg) io_full = io_full - 1
      if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

      ifile = io_full / niosh + 1
      jfile = jo_full / njosh + 1

      nump(ifile,jfile) = nump(ifile,jfile) + 1  ! Summation to # pts 
                                     ! filled by file (ifile,jfile) in dataset

   enddo
enddo

! Set up array index values for ptable array

ind = 1
do jfile = 1,jfiles
   do ifile = 1,ifiles
      numpind1(ifile,jfile) = ind   ! Beginning pt index for file (ifile,jfile)
      numpind2(ifile,jfile) = ind   ! For now, same as numpind1
      ind = ind + nump(ifile,jfile)
   enddo
enddo

! Fill ptable array

do iq = 2,nqa

   j1d = (iq - 1) * np
   do ip = 1,np

      pqlat1 = max( -89.9999,min( 89.9999,pqlat(ip,iq)))
      pqlon1 = max(-179.9999,min(179.9999,pqlon(ip,iq)))

      rio_full = (pqlon1 + 180.) * nperdeg ! must ignore pixel offset here
      rjo_full = (pqlat1 +  90.) * nperdeg ! must ignore pixel offset here

      io_full = int(rio_full)
      jo_full = int(rjo_full)

! If io_full and/or jo_full are at max value for database, decrease by 1 so
! that out-of-range file will not be sought.

      if (io_full == 360 * nperdeg) io_full = io_full - 1
      if (jo_full == 180 * nperdeg) jo_full = jo_full - 1

      ifile = io_full / niosh + 1
      jfile = jo_full / njosh + 1

      ind = numpind2(ifile,jfile)

      ptable(ind) = j1d + ip  ! point index in 1-d
      numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1 ! Sums to ending index
            
   enddo

enddo

! Read files and extract data

do jfile = 1,jfiles
   do ifile = 1,ifiles
   
      ind1 = numpind1(ifile,jfile)
      ind2 = numpind2(ifile,jfile)
   
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

         fname = trim(topo_database)//title1//title2//'.h5'
         inquire(file=fname, exist=l1, opened=l2)

! Read file

         if (l1) then
            write(io6,*) 'topo_database_read2 ',trim(fname)

            call shdf5_open(fname,'R')
            ndims=2 ; idims(1)=nio ; idims(2)=njo
            call shdf5_irec(ndims,idims,'topo',rvara=dato)
            call shdf5_close()
         else
            write(io6,*) 'Topography file is missing'
            write(io6,*) 'Topography filename = ',fname(1:lb)
            write(io6,*) 'Stopping model run'
            stop 'missing topo file'
         endif
         
         do ind = ind1,ind2-1

            ptab = ptable(ind)         
            ptab0 = ptab - 1

            iq = ptab0 / np + 1
            j1d = (iq - 1) * np

            ip = ptab - j1d

            pqlat1 = max( -89.9999,min( 89.9999,pqlat(ip,iq)))
            pqlon1 = max(-179.9999,min(179.9999,pqlon(ip,iq)))

            rio_full = (pqlon1 + 180.) * nperdeg ! must ignore pixel offset here
            rjo_full = (pqlat1 +  90.) * nperdeg ! must ignore pixel offset here

            io_full = int(rio_full)
            jo_full = int(rjo_full)

! If io_full and/or jo_full are at max value for database, decrease by 1 so
! that out-of-range file will not be sought.

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
            
! This is end of correction
            
            io2 = io1 + 1
            jo2 = jo1 + 1

            wio1 = 1. - wio2
            wjo1 = 1. - wjo2

            datp(ip,iq)                                                 &
               = wio1 * (wjo1 * dato(io1,jo1) +  wjo2 * dato(io1,jo2))  &
               + wio2 * (wjo1 * dato(io2,jo1) +  wjo2 * dato(io2,jo2))
               
         enddo
         
      endif
   enddo
enddo

! Average datp points together to get each topq value

c1 = 1. / float(np)

do iq = 2,nqa

   topq(iq) = 0.
   do ip = 1,np
      topq(iq) = topq(iq) + datp(ip,iq)
   enddo
   topq(iq) = topq(iq) * c1

enddo

deallocate(nump,numpind1,numpind2,ptable,dato,datp,pqlat,pqlon)

return
end subroutine topo_database_read
