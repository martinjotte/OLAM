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
subroutine read_press_header()

use isan_coms, only: iyear, nprz, levpr, ivertcoord, secondlat, cntlon,  &
                     cntlat, xnelon, xnelat, itinc, inproj, gdatdy,  &
                     gdatdx, xswlat, xswlon, npry, nprx, ihh, idd, imm,  &
                     iyy, isversion, marker, innpr, imonth, idate, ihour,  &
                     ipoffset
use misc_coms, only: io6

implicit none

integer :: lv,n,iunit

! Read the header of input pressure file.

write(io6, *) ''
write(io6, *) 'Reading RALPH data header '//trim(innpr)
write(io6, *) ''

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
   read(11,*) inproj,gdatdx,gdatdy,xswlat,xswlon  &
             ,xnelat,xnelon,cntlat,cntlon,secondlat
   read(11,*) ivertcoord,(levpr(lv),lv=1,nprz)
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

if (abs(nprx      * gdatdx - 360.) > .1 .or.   &
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

subroutine pressure_stage(p_u,p_v,p_t,p_z,p_r  &
                         ,p_slp, p_sfp, p_sft, p_snow, p_sst)

use isan_coms,   only: pnpr, levpr, nprx, npry, nprz, nprz_rh
use consts_coms, only: rocp, p00, eps_vap
use misc_coms,   only: io6

implicit none

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
write(io6, *) '     Access RALPH format pressure level data'
write(io6, *) '*****************************************************'

do k = 1,nprz
   pnpr(k) = levpr(k) * 100.
enddo

! Call routine to fill pressure arrays from the chosen dataset.

call get_press (iunit,p_u,p_v,p_t,p_z,p_r  &
               ,p_slp, p_sfp, p_sft, p_snow, p_sst)

!!!!!!!! Be careful !!!!!!!!!
!  Check input humidity variable p_r.  Assume that if the maximum of the field 
!  is greater than 1.1 (which allows for some machine roundoff), 
!  it is specific humidity (in g/kg) which needs to be converted to kg/kg,
!  else it is R.H. which needs to be converted to specifit humidity

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

close(iunit)

write(io6,*) ''
write(io6,*) 'nprz2 ',nprz

end subroutine pressure_stage

!===============================================================================

subroutine get_press (iunit,p_u,p_v,p_t,p_z,p_r  &
                     ,p_slp, p_sfp, p_sft, p_snow, p_sst)

use max_dims,  only: maxpr
use isan_coms, only: nprz, npry, nprx, pnpr, iyear, imonth, idate,  &
                     ihour, levpr, ipoffset, nprz_rh
use misc_coms,   only: io6

implicit none

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

real :: as(nprx,npry)  ! automatic array

integer :: i,j,k,nv,nvar,misstot,lv,n
integer :: ithere(maxpr,5),isfthere(5)
character(len=1) :: idat(5) = (/ 'T','R','U','V','H' /)

! Initialize with missing data flag

ithere   = -999
isfthere = -999

!  Read upper air fields

write(io6,*) ' '

do lv = 1,nprz

! Zonal wind component

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_u(1,1,lv))
   write(io6, '('' ==  Read UE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)')  &
      levpr(lv),ihour,idate,imonth,iyear

! Meridional wind component

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_v(1,1,lv))
   write(io6, '('' ==  Read VE on P lev '',i4,'' at UTC '',i6.4,2i3,i5)')  &
      levpr(lv),ihour,idate,imonth,iyear

! Temperature

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_t(1,1,lv))
   write(io6, '('' ==  Read T on P lev '',i4,'' at UTC '',i6.4,2i3,i5)')  &
      levpr(lv),ihour,idate,imonth,iyear
   write(io6,*) 'prread1 ',as(5,5),nprx,npry,ipoffset,p_t(5,5,lv)

! Geopotential height

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_z(1,1,lv))
   write(io6, '('' ==  Read Z on P lev '',i4,'' at UTC '',i6.4,2i3,i5)')  &
      levpr(lv),ihour,idate,imonth,iyear

! Relative humidity (or specific humidity)

   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
   call prfill(nprx,npry,ipoffset,as,p_r(1,1,lv))
   write(io6, '('' ==  Read RH on P lev '',i4,'' at UTC '',i6.4,2i3,i5)')  &
      levpr(lv),ihour,idate,imonth,iyear

enddo

!  Read surface fields

write(io6,*) ' '

! SLP

read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
call prfill(nprx,npry,ipoffset,as,p_slp(1,1))
write(io6, '('' ==  Read SLP at sfc at UTC '',i6.4,2i3,i5)')  &
   ihour,idate,imonth,iyear

! SFP

read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
call prfill(nprx,npry,ipoffset,as,p_sfp(1,1))
write(io6, '('' ==  Read SFP at sfc at UTC '',i6.4,2i3,i5)')  &
   ihour,idate,imonth,iyear

! SFT

read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
call prfill(nprx,npry,ipoffset,as,p_sft(1,1))
write(io6, '('' ==  Read SFT at sfc at UTC '',i6.4,2i3,i5)')  &
   ihour,idate,imonth,iyear

! SNOW

read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
call prfill(nprx,npry,ipoffset,as,p_snow(1,1))
write(io6, '('' ==  Read SNOW at sfc at UTC '',i6.4,2i3,i5)')  &
   ihour,idate,imonth,iyear

! SST at surface

read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
call prfill(nprx,npry,ipoffset,as,p_sst(1,1))
write(io6, '('' ==  Read SST at sfc at UTC '',i6.4,2i3,i5)')  &
   ihour,idate,imonth,iyear

write(io6,*) ''

goto 71

70 CONTINUE
write(io6,*) 'Premature end of file or error in pressure input file!'
write(io6,*) 'We''ll close our eyes and pretend it didn''t happen!'
71 continue

! Special for RH:
! Reanalysis typically only reports RH up to 100mb, but still reports the other
! fields up to higher levels. Check for the highest level that reports RH:
do k=nprz, 1, -1
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


!!if (misstot > 0) then
!!   ! Let's see if we can get creative and make up data for the missing fields
!!
!!   call press_miss(nprx+3,npry+2,nprz,p_u,p_v,p_t,p_z,p_r  &
!!                  ,ithere,maxpr,levpr)
!!
!!  ! Check again for missing fields
!!
!!   misstot = 0
!!   do nv = 1,5
!!      do k = 1,nprz
!!         if (ithere(k,nv) == 0) misstot = misstot + 1
!!      enddo
!!   enddo
!!
!!   write(io6,*) '------------------------------------------------'
!!   write(io6,*) ' After missing parameters check: 0 = all missing'
!!   write(io6,*) '------------------------------------------------'
!!   write(io6, '(t20,5(a1,6x))') (idat(n),n=1,5)
!!
!!   do k = 1,nprz
!!      write(io6, '(f10.1,t14,5(i7))',pnpr(k),(ithere(k,n),n=1,5)
!!   enddo
!!
!!endif

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


!!!===============================================================================
!!
!!subroutine press_miss (n1,n2,n3,un,vn,tn,zn,rn,ithere,maxpr,levpr)
!!
!!implicit none
!!
!!integer, intent(in) :: n1,n2,n3,maxpr,levpr(*)
!!
!!real, intent(inout), dimension(n1,n2,n3) ::  un,vn,tn,zn,rn
!!
!!integer, intent(inout) :: ithere(maxpr,5)
!!
!!real :: prs(n3),prsln(n3)
!!integer :: it=1,ir=2,iu=3,iv=4,iz=5
!!integer :: ierr,i,j,k
!!
!!do k = 1,n3
!!   prs(k) = float(levpr(k))
!!   prsln(k) = log(prs(k))
!!enddo
!!
!!! first do moisture. since there are no physical relationships we can
!!!   use to help us, we will simply interpolate to a missing level. If
!!!   the top or bottom level is missing, we will fill it with the relative
!!!   humidity above or below.
!!
!!call pr_miss_fill (n1,n2,n3,rn,ithere(1,ir),ithere(1,ir)  &
!!     ,prsln,'rel hum',ierr)
!!if (ierr == 1) return
!!
!!! do temperature in a similar manner, but only if there are levels
!!!   where height is missing also
!!
!!call pr_miss_fill (n1,n2,n3,tn,ithere(1,it),ithere(1,iz)  &
!!     ,prsln,'temp',ierr)
!!if (ierr == 1) return
!!
!!
!!! check height at each level. It can be computed hydrostatically
!!!   from temp as long as temp is not missing at this level and temp and
!!!   height is available on a level above or below.
!!!   Make a downward sweep first, then an upward.
!!
!!
!!do k = n3-1,1,-1
!!   if (ithere(k,iz) == 0) then
!!      ! check for temp, z, above
!!      if (ithere(k+1,iz) == 1 .and. ithere(k+1,it) == 1 .and.  &
!!         ithere(k,it) == 1) then
!!         !ok, we can get this
!!         call pr_hystatic_z (n1*n2,zn(1,1,k),zn(1,1,k+1)  &
!!                            ,tn(1,1,k),tn(1,1,k+1),rn(1,1,k),rn(1,1,k+1)  &
!!                            ,prs(k),prs(k+1))
!!         ithere(k,iz) = 1
!!         write(io6,*) '-->Computing hydrostatic pressure at level:',k
!!      endif
!!   endif
!!enddo
!!
!!do k = 2,n3
!!   if (ithere(k,iz) == 0) then
!!      ! check for temp, z, below
!!      if(ithere(k-1,iz) == 1 .and. ithere(k-1,it) == 1 .and.  &
!!           ithere(k,it) == 1) then
!!         ! ok, we can get this
!!         call pr_hystatic_z (n1*n2,zn(1,1,k),zn(1,1,k-1)  &
!!                            ,tn(1,1,k),tn(1,1,k-1),rn(1,1,k),rn(1,1,k-1)  &
!!                            ,prs(k),prs(k-1))
!!         write(io6,*) '-->Computing hydrostatic pressure at level:',k
!!         ithere(k,iz) = 1 
!!      endif
!!   endif
!!enddo
!!
!!! try temperature in a similar manner. It can also be computed
!!!   from height as long as height is not missing at this level and temp and
!!!   height is available on a level above or below.
!!!   Note that vapor mixing ratio (used to compute virtual temperature)
!!!   is somewhat difficult to compute from relative humidity if temperature
!!!   is missing. Therefore, the vitual temperature factor at the
!!!   non-missing level is assumed at the missing level.
!!
!!do k = n3-1,1,-1
!!   if (ithere(k,it) == 0) then
!!      ! check for temp, z, above
!!      if (ithere(k+1,iz) == 1 .and. ithere(k+1,it) == 1 .and.  &
!!         ithere(k,iz) == 1) then
!!         ! ok, we can get this
!!         call pr_hystatic_t (n1*n2,zn(1,1,k),zn(1,1,k+1)  &
!!                           ,tn(1,1,k),tn(1,1,k+1),rn(1,1,k),rn(1,1,k+1)  &
!!                           ,prs(k),prs(k+1))
!!         ithere(k,it) = 1
!!         write(io6,*) '-->Computing hydrostatic temperature at level:',k
!!      endif
!!   endif
!!enddo
!!
!!do k = 2,n3
!!   if (ithere(k,it) == 0) then
!!      ! check for temp, z, below
!!      if (ithere(k-1,iz) == 1 .and. ithere(k-1,it) == 1 .and.  &
!!         ithere(k,iz) == 1) then
!!         ! ok, we can get this
!!         call pr_hystatic_t (n1*n2,zn(1,1,k),zn(1,1,k-1)  &
!!                            ,tn(1,1,k),tn(1,1,k-1),rn(1,1,k),rn(1,1,k-1)  &
!!                            ,prs(k),prs(k-1))
!!         ithere(k,it) = 1
!!         write(io6,*) '-->Computing hydrostatic temperature at level:',k
!!      endif
!!   endif
!!enddo
!!
!!! For the u and v components, do a straight interpolation again like we
!!!   did for rel humidity.
!!
!!call pr_miss_fill (n1,n2,n3,un,ithere(1,iu),ithere(1,iu),prsln,'u-comp',ierr)
!!if (ierr == 1) return
!!
!!call pr_miss_fill (n1,n2,n3,vn,ithere(1,iv),ithere(1,iv),prsln,'v-comp',ierr)
!!if (ierr == 1) return
!!
!!return
!!end
!!
!!!===============================================================================
!!
!!subroutine pr_hystatic_z(np,z1,z2,t1,t2,r1,r2,p1,p2)
!!
!!use consts_coms, only: eps_virt, grav, rdry, grav2
!!
!!implicit none
!!
!!integer, intent(in) :: np
!!
!!real, intent(in) :: z2(np),t2(np),r1(np),r2(np),p1,p2
!!real, intent(inout) :: z1(np),t1(np)
!!
!!integer :: n
!!real :: tv1,rslf,tv2,vtfact
!!
!!do n = 1,np
!!   if (z2(n) < 1.e20 .and. t1(n) < 1.e20 .and. t2(n) < 1.e20 .and.  &
!!       r1(n) < 1.e20 .and. r2(n) < 1.e20 ) then
!!
!!      tv1 = t1(n) * (1. + eps_virt * rslf(p1*100.,t1(n)) * r1(n))
!!      tv2 = t2(n) * (1. + eps_virt * rslf(p2*100.,t2(n)) * r2(n))
!!      z1(n) = z2(n) - rdry * .5 * (tv1 + tv2)  &
!!            * (log(p1 * 100.) - log(p2 * 100.)) / grav
!!
!!   else
!!      z1(n) = 1.e30
!!   endif
!!enddo
!!
!!return
!!
!!entry pr_hystatic_t(np,z1,z2,t1,t2,r1,r2,p1,p2)
!!
!!do n = 1,np
!!   if (t2(n) < 1.e20 .and. z1(n) < 1.e20 .and. z2(n) < 1.e20 .and.  &
!!       r1(n) < 1.e20 .and. r2(n) < 1.e20 ) then
!!
!!      vtfact = (1. + eps_virt * rslf(p2*100.,t2(n)) * r2(n))
!!
!!      t1(n) = -t2(n) - (grav2 * (z1(n) - z2(n))  &
!!            / (rdry * (log(p1 * 100.) - log(p2 * 100.))) ) / vtfact
!!
!!   else
!!      t1(n) = 1.e30
!!   endif
!!enddo
!!
!!return
!!end
!!
!!!===============================================================================
!!
!!subroutine pr_miss_fill (n1,n2,n3,a,ithere,ithere2,prsln,varname,ierr)
!!
!!! simple filling for missing pressure fields. Assuming there are no
!!!   physical relationships we can
!!!   use to help us, we will simply interpolate to a missing level. If
!!!   the top or bottom level is missing, we will fill it with the next
!!!   non-missing field above or below.
!!
!!implicit none
!!
!!integer, intent(in) :: n1,n2,n3
!!real, intent(inout) :: a(n1,n2,n3)
!!integer, intent(inout) :: ithere(n3)
!!integer, intent(in) :: ithere2(n3)
!!real, intent(in) :: prsln(n3)
!!character(len=*), intent(in) :: varname
!!
!!integer, intent(out) :: ierr
!!
!!
!!integer :: k,kk,ilev
!!
!!ierr = 0
!!k = 1
!!if (ithere(k) == 0 .and. ithere2(k) == 0) then
!!   ! find first non-missing level above
!!   ilev = 0
!!   do kk = 2,n3
!!      if (ithere(kk) == 1) then
!!         ilev = kk
!!         goto 10
!!      endif
!!   enddo
!!   write(io6,*) 'All data missing; fixing stopped for:',varname
!!   ierr = 1
!!   return
!!   10 continue
!!   a(1:n1,1:n2,k) = a(1:n1,1:n2,ilev)
!!   ithere(k) = 1
!!   write(io6,*) '-->Filling:',varname,' at level:',k,' from level:',ilev
!!endif
!!
!!k = n3
!!if (ithere(k) == 0 .and. ithere2(k) == 0) then
!!   ! find first non-missing level below
!!   ilev = 0
!!   do kk = n3-1,1,-1
!!      if (ithere(kk) == 1) then
!!         ilev = kk
!!         goto 11
!!      endif
!!   enddo
!!   write(io6,*) 'All data missing; fixing stopped for:',varname
!!   ierr = 1
!!   return
!!   11 continue
!!   a(1:n1,1:n2,k) = a(1:n1,1:n2,ilev)
!!   ithere(k) = 1
!!   write(io6,*) '-->Filling:',varname,' at level:',k,' from level:',ilev
!!endif
!!
!!! Now interpolate to internal missing levels
!!
!!do k = 2,n3-1
!!   if (ithere(k) == 0 .and. ithere2(k) == 0) then
!!      ! find first non-missing level above
!!      ilev = 0
!!      do kk = k+1,n3
!!         if (ithere(kk) == 1) then
!!            ilev = kk
!!            goto 12
!!         endif
!!      enddo
!!      stop 'pr_miss_fill'
!!      12 continue
!!      call pr_interp (n1*n2,a(1,1,k-1),a(1,1,k),a(1,1,ilev)  &
!!                     ,prsln(k-1),prsln(k),prsln(ilev))
!!      ithere(k) = 1
!!      write(io6,*) '-->Interpolating:',varname,' at level:',k  &
!!            ,' from levels:',k-1,' and:',ilev
!!
!!   endif
!!enddo
!!
!!return
!!end
!!
!!!===============================================================================
!!
!!subroutine pr_interp(np,a1,a2,a3,pln1,pln2,pln3)
!!
!!implicit none
!!
!!integer, intent(in) :: np
!!
!!real, intent(in) :: a1(np),a3(np),pln1,pln2,pln3
!!real, intent(out) :: a2(np)
!!
!!integer :: n
!!
!!do n = 1,np
!!   if (a3(n) < 1.e20 .and. a1(n) < 1.e20) then
!!      a2(n) = a1(n) +  (pln2 - pln1) * (a3(n) - a1(n)) / (pln3 - pln1)
!!   else
!!      a2(n) = 1.e30
!!   endif
!!enddo
!!
!!return
!!end subroutine pr_interp

