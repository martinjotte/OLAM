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
subroutine read_press_header(fform)

use isan_coms,  only: iyear, nprz, levpr, ivertcoord, secondlat, cntlon, &
                      cntlat, xnelon, xnelat, itinc, inproj, gdatdy, &
                      gdatdx, xswlat, xswlon, npry, nprx, ihh, idd, imm, &
                      iyy, isversion, marker, innpr, imonth, idate, ihour, &
                      ipoffset
use misc_coms,  only: io6
use hdf5_utils, only: shdf5_open, shdf5_irec

implicit none

character(len=3), intent(inout) :: fform

character(len=16) :: ext

integer, external :: lastdot
integer :: lv,n,iunit
integer :: ndims, idims(2)

! Read the header of input pressure file.

write(io6,'(/,a,/)') 'Reading pressure gridded data header '//trim(innpr)

! Set flag for type of file:
!     "GDF" - text Ralph II file (none or .vfm extension)
!     "HD5"  - HDF5 GDF file ("h5" or "hdf5")
!                     2 or 3D field records (x,y) or (x,y,z/p)
!                       (1,1)   - (SW lon,lat) --
!                       (1,1,1) - (SW lon,lat,bottom z/p level)
!                     range of header info

! Find file name extension

ext = trim( innpr(lastdot(innpr)+1:) )
write(io6,*) 'ffffff-ext:-',trim(ext),'----',lastdot(innpr)

fform = "GDF"
if (trim(ext) == 'h5' .or. trim(ext) == 'hdf5' .or. &
    trim(ext) == 'H5' .or. trim(ext) == 'HDF5') &
    fform = 'HD5'

if (fform == 'GDF') then

   open(11,file=innpr)
   read(11,*) marker,isversion
   if(marker.ne.999999) isversion=1

   if (isversion == 1) then
      rewind 11
      read(11,*) iyy,imm,idd,ihh,nprz,nprx,npry,xswlon,xswlat,gdatdx,gdatdy
      read(11,*) (levpr(n),n=1,nprz)
      inproj = 1
      ihh = ihh * 100
   elseif (isversion == 2) then
      write(io6,*) 'doing RALPH 2 format'
      read(11,*) iyy,imm,idd,ihh,itinc,nprz,nprx,npry
      read(11,*) inproj,gdatdx,gdatdy,xswlat,xswlon &
                ,xnelat,xnelon,cntlat,cntlon,secondlat
      read(11,*) ivertcoord,(levpr(lv),lv=1,nprz)
   endif

elseif (fform == 'HD5') then

   call shdf5_open (trim(innpr), 'R')

   ndims = 1
   idims(1) = 1
   idims(2) = 1

   call shdf5_irec(ndims, idims, 'version',ivars=isversion)
   call shdf5_irec(ndims, idims, 'year'   ,ivars=iyy)
   call shdf5_irec(ndims, idims, 'month'  ,ivars=imm)
   call shdf5_irec(ndims, idims, 'day'    ,ivars=idd)
   call shdf5_irec(ndims, idims, 'hour'   ,ivars=ihh)
   call shdf5_irec(ndims, idims, 'ftime'  ,ivars=itinc)
   call shdf5_irec(ndims, idims, 'nx'     ,ivars=nprx)
   call shdf5_irec(ndims, idims, 'ny'     ,ivars=npry)
   call shdf5_irec(ndims, idims, 'nlev'   ,ivars=nprz)
   call shdf5_irec(ndims, idims, 'iproj'  ,ivars=inproj)
   call shdf5_irec(ndims, idims, 'vcoord' ,ivars=ivertcoord)
   call shdf5_irec(ndims, idims, 'swlat'  ,rvars=xswlat)
   call shdf5_irec(ndims, idims, 'swlon'  ,rvars=xswlon)
   call shdf5_irec(ndims, idims, 'nelat'  ,rvars=xnelat)
   call shdf5_irec(ndims, idims, 'nelon'  ,rvars=xnelon)
   call shdf5_irec(ndims, idims, 'dx'     ,rvars=gdatdx)
   call shdf5_irec(ndims, idims, 'dy'     ,rvars=gdatdy)
   call shdf5_irec(ndims, idims, 'reflat1',rvars=cntlat)
   call shdf5_irec(ndims, idims, 'reflat2',rvars=secondlat)

   idims(1) = nprz

   call shdf5_irec(ndims, idims, 'levels' ,ivara=levpr)

endif

write(io6,*) 'nprz1 ',nprz,nprx,npry

! Check for consistency between file parameters and namelist parameters

if (iyy /= iyear .or. imm /= imonth .or. idd /= idate .or. ihh /= ihour) then

   write(io6,*) 'Pressure file dates not the same as namelist!'
   write(io6,*) 'Year :',iyy,iyear
   write(io6,*) 'Month:',imm,imonth
   write(io6,*) 'Day  :',idd,idate
   write(io6,*) 'Hour :',ihh,ihour
   stop 'pr_dates'
endif

! Check pressure data domain size and location

if (inproj /= 1) then
   write(io6,*) 'You must input a lat-lon pressure grid '
   stop 'glob-no-press'
endif

! We make the requirement that a full global domain of pressure-level data
! be read in.  Check this here.  Following the convention for the NCEP/DOE
! Reanalysis2 data, assume that data exists at both latitudinal boundaries
! (-90. and 90. degrees) but that the longitudinal boundary is not repeated.
! If either is not the case, this check will stop execution.

if (abs(nprx      * gdatdx - 360.) > .1 .or. &
    abs ((npry-1) * gdatdy - 180.) > .1) then

    write(io6,*) 'Gridded pressure level data must have global coverage'
    write(io6,*) 'nprx,npry = ',nprx,npry
    write(io6,*) 'gdatdx,gdatdy = ',gdatdx,gdatdy
    write(io6,*) 'xswlat,xswlon= ',xswlat,xswlon
    stop 'astp - non-global domain in input pressure data'
endif

! Compute longitudinal offset index, which is the index in the expanded isan
! pressure arrays where the first input data point (at xswlon) is locatd.

ipoffset = (xswlon + 180.) / gdatdx + 1

end subroutine read_press_header

!===============================================================================

subroutine pressure_stage(fform,p_u,p_v,p_t,p_z,p_r, &
                          p_slp, p_sfp, p_sft, p_snow, p_sst)

use isan_coms,   only: pnpr, levpr, nprx, npry, nprz, nprz_rh
use consts_coms, only: rocp, p00, eps_vap
use misc_coms,   only: io6
use hdf5_utils,  only: shdf5_close

implicit none

character(len=*), intent(in) :: fform

real, intent(out) :: p_u(nprx+3,npry+2,nprz)
real, intent(out) :: p_v(nprx+3,npry+2,nprz)
real, intent(out) :: p_t(nprx+3,npry+2,nprz)
real, intent(out) :: p_z(nprx+3,npry+2,nprz)
real, intent(out) :: p_r(nprx+3,npry+2,nprz)

real, intent(out) :: p_slp (nprx+3,npry+2)
real, intent(out) :: p_sfp (nprx+3,npry+2)
real, intent(out) :: p_sft (nprx+3,npry+2)
real, intent(out) :: p_snow(nprx+3,npry+2)
real, intent(out) :: p_sst (nprx+3,npry+2)

real :: thmax,thmin,vapor_press
integer :: i,j,k,lv,n,iunit

real, external :: eslf

iunit = 11

write(io6, *) ''
write(io6, *) '*****************************************************'
write(io6, *) '     Access pressure level data'
write(io6, *) '*****************************************************'

do k = 1,nprz
   pnpr(k) = levpr(k) * 100.
enddo

! Call routine to fill pressure arrays from the chosen dataset.

call get_press (fform,iunit,p_u,p_v,p_t,p_z,p_r, &
                p_slp, p_sfp, p_sft, p_snow, p_sst)

!!!!!!!! Be careful !!!!!!!!!
!  Check input humidity variable p_r.  Assume that if the maximum of the field
!  is greater than 1.1 (which allows for some machine roundoff),
!  it is specific humidity (in g/kg) which needs to be converted to kg/kg,
!  else it is R.H. which needs to be converted to specific humidity

if (maxval(p_r(1:nprx+3,1:npry+2,1:nprz_rh)) > 1.1) then

   ! Convert specific humidity to kg/kg units

   write(io6, *) ''
   write(io6, *) 'Converting specific humidity to kg/kg'

   do k = 1,nprz
      do j = 1,npry+2
         do i = 1,nprx+3

            p_r(i,j,k) = .001 * p_r(i,j,k)

         enddo
      enddo
   enddo

else

   ! Convert R.H. to specific humidity

   write(io6, *) ''
   write(io6, *) 'Converting relative humidity to specific humidity'

   do k = 1,nprz
      do j = 1,npry+2
         do i = 1,nprx+3

            ! Compute ambient vapor pressure based on R.H.
            ! and saturation vapor pressure (eslf)

            vapor_press = p_r(i,j,k) * eslf(p_t(i,j,k)-273.15)

            ! Do not allow vapor pressure to exceed ambient pressure

            vapor_press = min(pnpr(k),vapor_press)

            ! Compute specific humidity from vapor press and ambient press

            p_r(i,j,k) = eps_vap * vapor_press &
                 / (pnpr(k) + vapor_press * (eps_vap - 1.))

         enddo
      enddo
   enddo

endif

! Print max-min theta at bottom and top levels

thmax = maxval(p_t(:,:,nprz)) * (p00/pnpr(nprz))**rocp
thmin = minval(p_t(:,:,1   )) * (p00/pnpr(1))**rocp

write(io6, *) ''
write(io6, "(' Minimum THETA at ',I4,' mb: ',F9.3)") levpr(1), thmin
write(io6, "(' Maximum THETA at ',I4,' mb: ',F9.3)") levpr(nprz), thmax

if (fform == 'GDF') then
   close(iunit)
elseif (fform == 'HD5') then
   call shdf5_close()
endif

write(io6,*) ''
write(io6,*) 'nprz2 ',nprz

end subroutine pressure_stage

!===============================================================================

subroutine get_press (fform,iunit,p_u,p_v,p_t,p_z,p_r, &
                      p_slp, p_sfp, p_sft, p_snow, p_sst)

use max_dims,   only: maxpr
use isan_coms,  only: nprz, npry, nprx, pnpr, iyear, imonth, idate, &
                      ihour, levpr, ipoffset, nprz_rh
use misc_coms,  only: io6
use hdf5_utils, only: shdf5_irec

implicit none

character(len=*), intent(in) :: fform
integer, intent(in) :: iunit

real, intent(out) :: p_u(nprx+3,npry+2,nprz)
real, intent(out) :: p_v(nprx+3,npry+2,nprz)
real, intent(out) :: p_t(nprx+3,npry+2,nprz)
real, intent(out) :: p_z(nprx+3,npry+2,nprz)
real, intent(out) :: p_r(nprx+3,npry+2,nprz)

real, intent(out) :: p_slp (nprx+3,npry+2)
real, intent(out) :: p_sfp (nprx+3,npry+2)
real, intent(out) :: p_sft (nprx+3,npry+2)
real, intent(out) :: p_snow(nprx+3,npry+2)
real, intent(out) :: p_sst (nprx+3,npry+2)

real :: as(nprx,npry)
real :: as3(nprx,npry,nprz)

integer :: i,j,k,nv,nvar,misstot,lv,n
integer :: ithere(maxpr,5),isfthere(5)
character(len=1) :: idat(5) = (/ 'T','R','U','V','H' /)
integer :: ndims, idims(3)

! Initialize with missing data flag

ithere   = -999
isfthere = -999

!  Read upper air fields

write(io6,*) ' '

if (fform == 'GDF') then

   do lv = 1,nprz

! Zonal wind component

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,ipoffset,as,p_u(1,1,lv))
      write(io6, '('' ==  Read UE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Meridional wind component

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,ipoffset,as,p_v(1,1,lv))
      write(io6, '('' ==  Read VE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Temperature

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,ipoffset,as,p_t(1,1,lv))
      write(io6, '('' ==  Read T on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear
      write(io6,*) 'prread1 ',as(5,5),nprx,npry,ipoffset,p_t(5,5,lv)

! Geopotential height

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,ipoffset,as,p_z(1,1,lv))
      write(io6, '('' ==  Read Z on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

! Relative humidity (or specific humidity)

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
      call prfill(nprx,npry,ipoffset,as,p_r(1,1,lv))
      write(io6, '('' ==  Read RH on P lev '',i4,'' at UTC '',i6.4,2i3,i5)') &
         levpr(lv),ihour,idate,imonth,iyear

   enddo

!  Read surface fields

   write(io6,*) ' '

! SLP

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_slp(1,1))
   write(io6, '('' ==  Read SLP at sfc at UTC '',i6.4,2i3,i5)') &
      ihour,idate,imonth,iyear

! SFP

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_sfp(1,1))
   write(io6, '('' ==  Read SFP at sfc at UTC '',i6.4,2i3,i5)') &
      ihour,idate,imonth,iyear

! SFT

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_sft(1,1))
   write(io6, '('' ==  Read SFT at sfc at UTC '',i6.4,2i3,i5)') &
      ihour,idate,imonth,iyear

! SNOW

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_snow(1,1))
   write(io6, '('' ==  Read SNOW at sfc at UTC '',i6.4,2i3,i5)') &
      ihour,idate,imonth,iyear

! SST at surface

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_sst(1,1))
   write(io6, '('' ==  Read SST at sfc at UTC '',i6.4,2i3,i5)') &
      ihour,idate,imonth,iyear

   write(io6,*) ''

   goto 71

   70 CONTINUE
   write(io6,*) 'Premature end of file or error in pressure input file!'
   write(io6,*) 'We''ll close our eyes and pretend it didn''t happen!'
   71 continue

elseif (fform == 'HD5') then

! Read 3D met vars

   ndims = 3
   idims(1) = nprx
   idims(2) = npry
   idims(3) = nprz

   call shdf5_irec(ndims, idims,'UP',rvara = as3)
   call prfill3(nprx,npry,nprz,ipoffset,as3,p_u)

   call shdf5_irec(ndims, idims,'VP',rvara = as3)
   call prfill3(nprx,npry,nprz,ipoffset,as3,p_v)

   call shdf5_irec(ndims, idims,'TEMP',rvara = as3)
   call prfill3(nprx,npry,nprz,ipoffset,as3,p_t)

   call shdf5_irec(ndims, idims,'GEO',rvara = as3)
   call prfill3(nprx,npry,nprz,ipoffset,as3,p_z)

   call shdf5_irec(ndims, idims,'RELHUM',rvara = as3)
   call prfill3(nprx,npry,nprz,ipoffset,as3,p_r)

!   if (ivertcoord == 3) call shdf5_irec('PRESS',rvara = p_p)

   print*,'uu max-min:', maxval(p_u), minval(p_u)
   print*,'vv max-min:', maxval(p_v), minval(p_v)
   print*,'tt max-min:', maxval(p_t), minval(p_t)
   print*,'zz max-min:', maxval(p_z), minval(p_z)
   print*,'rr max-min:', maxval(p_r), minval(p_r)
!   print*,'pp max-min:', maxval(p_p), minval(p_p)

   ! 3D scalar vars

!   if (num_scvars > 0) then
!      do nv = 1, num_scvars
!         call shdf5_irec(trim(in_scvars(nv)),rvara = p_scalar(1,1,1,nv))

!         print*,trim(in_scvars(nv)),' max-min:', maxval(p_scalar(:,:,:,nv)), &
!                      minval(p_scalar(:,:,:,nv))

!      enddo
!   endif

! Read 2D met vars

   ndims = 2
   idims(1) = nprx
   idims(2) = npry
   idims(3) = 1

   call shdf5_irec(ndims, idims, 'PRMSL'    ,rvara = p_slp)
   call shdf5_irec(ndims, idims, 'PRSFC'    ,rvara = p_sfp)
   call shdf5_irec(ndims, idims, 'TSFC'     ,rvara = p_sft)
   call shdf5_irec(ndims, idims, 'SNOW_MASS',rvara = p_snow)
   call shdf5_irec(ndims, idims, 'SST'      ,rvara = p_sst)

endif

! Special for RH:
! Reanalysis typically only reports RH up to 100mb, but still reports the other
! fields up to higher levels. Check for the highest level that reports RH:

do k = nprz,1,-1
   if (all(p_r(:,:,k) > -998.)) then
      nprz_rh = k
      exit
   endif
enddo

write(io6, *) '----------------------------------------------------'
write(io6, "(A,I0,A,I0,A)") ' Pressure-level data has ', nprz, &
     ' levels, goes up to ', levpr(nprz), ' mb.'
write(io6, "(A,I0,A,I0,A)") ' Water vapor has ', nprz_rh, &
     ' levels, is reported up to ', levpr(nprz_rh), ' mb.'
write(io6, *) ''

! Check for missing data

do k = 1,nprz
   ithere(k,1) = count(p_t(:,:,k) < -998.0)
   ithere(k,5) = count(p_r(:,:,k) < -998.0)
   ithere(k,2) = count(p_u(:,:,k) < -998.0)
   ithere(k,3) = count(p_v(:,:,k) < -998.0)
   ithere(k,4) = count(p_z(:,:,k) < -998.0)
enddo

where (ithere(:,:) > nprx*npry) ithere(:,:) = -1

isfthere(1) = count(p_slp (:,:) < -998.0)
isfthere(2) = count(p_sfp (:,:) < -998.0)
isfthere(3) = count(p_sft (:,:) < -998.0)
isfthere(4) = count(p_snow(:,:) < -998.0)
isfthere(5) = count(p_sst (:,:) < -998.0)

where (isfthere(:) > nprx*npry) isfthere(:) = -1

write(io6,*) '---------------------------------------------------'
write(io6,*) ' # of missing values per level (-1 = all missing): '
write(io6,*) '---------------------------------------------------'
write(io6, '(a,5(6x,a1))') '   P (mb)', (idat(n),n=1,5)

do k = 1,nprz
   write(io6, '(f10.2,t10,5(i7))') pnpr(k)/100., (ithere(k,n),n=1,5)
enddo

write(io6,*) ''
write(io6,*) '---------------------------------------------------'
write(io6,*) ' # of missing surface values (-1 = all missing):   '
write(io6,*) '---------------------------------------------------'
write(io6,*) '             SLP    SFP    SFT    SNOW   SST       '

write(io6, '(t10,5(i7))') (isfthere(n),n=1,5)

if (any(ithere(1:nprz,(/1,3,4,5/)) > 0) .or. any(ithere(1:nprz_rh,2) > 0)) then
   write(io6,*) ''
   write(io6,*) '-------------------------------------------------------'
   write(io6,*) 'WARNING - There appears to be missing data in the input'
   write(io6,*) '          pressure-level data files'
   write(io6,*) '-------------------------------------------------------'
endif

end subroutine get_press

!===============================================================================

subroutine prfill (nprx,npry,ipoffset,xx,dn)

implicit none

integer, intent(in) :: nprx,npry,ipoffset

real, intent(in)  :: xx(nprx,npry)
real, intent(out) :: dn(nprx+3,npry+2)

integer :: i,j,nprxo2,iin

nprxo2 = nprx / 2

! Copy pressure-level or surface data from xx to dn array, shifting by 1 row
! and column, and set any required missing values

do j = 1,npry
   do i = 1,nprx+3
      iin = mod(i+nprx-ipoffset-1,nprx)+1
      dn(i,j+1) = xx(iin,j)
   enddo
enddo

! Fill in N/S boundary points

do i = 1,nprxo2
   dn(i,1)      = dn(i+nprxo2,3)
   dn(i,npry+2) = dn(i+nprxo2,npry)
enddo

do i = nprxo2+1,nprx+3
   dn(i,1)      = dn(i-nprxo2,3)
   dn(i,npry+2) = dn(i-nprxo2,npry)
enddo

return
end subroutine prfill

!===============================================================================

subroutine prfill3 (nprx,npry,nprz,ipoffset,xx,dn)

implicit none

integer, intent(in) :: nprx,npry,nprz,ipoffset

real, intent(in)  :: xx(nprx,npry,nprz)
real, intent(out) :: dn(nprx+3,npry+2,nprz)

integer :: i,j,k,nprxo2,iin

nprxo2 = nprx / 2

! Copy pressure-level or surface data from xx to dn array, shifting by 1 row
! and column, and set any required missing values

do k = 1,nprz
   do j = 1,npry
      do i = 1,nprx+3
         iin = mod(i+nprx-ipoffset-1,nprx)+1
         dn(i,j+1,k) = xx(iin,j,k)
      enddo
   enddo

! Fill in N/S boundary points

   do i = 1,nprxo2
      dn(i,1,k)      = dn(i+nprxo2,3,k)
      dn(i,npry+2,k) = dn(i+nprxo2,npry,k)
   enddo

   do i = nprxo2+1,nprx+3
      dn(i,1,k)      = dn(i-nprxo2,3,k)
      dn(i,npry+2,k) = dn(i-nprxo2,npry,k)
   enddo
enddo

return
end subroutine prfill3

!===============================================================================

integer function lastdot(str)
implicit none
character(len=*) :: str
integer :: n,ln

! returns last . character position from a string
!     trailing blanks are ignored

ln=len_trim(str)
do n=ln,1,-1
   if(str(n:n) == '.') then
      lastdot = n
      return
   endif
enddo
lastdot = 0

return
end

