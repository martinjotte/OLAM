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

! Arrays zont12 and zonu12 will contain the entire ZONAVG

real :: zont12(74,22,12)
real :: zonu12(74,22,12)
real :: zont(74,22)
real :: zonu(74,22)
real :: zonr(74,22)
real :: zonz(74,22)
real :: zonp_vect(22)
real :: zonr_vect(22)
real :: zonz_vect(22)
real :: zont_vect(22)
real :: zonu_vect(22)

Contains

!===============================================================================

  subroutine zonavg_init(idate,imonth,iyear)

  use misc_coms,   only: io6, zonclim
  use consts_coms, only: pio180, cp, rocp, eps_virt, omega, grav2, dlat

  implicit none   

  integer, intent(in) :: idate,imonth,iyear

  integer :: iplev,ilat,ilatn,ilats
  real :: dyo2g,fcorn,fcors,pilo,pihi,thetavlo,thetavhi,cpo2g

! Read in 'ZONAVG_CLIMATE' dataset (if not already read in), and interpolate
! to the current time.  This fills array ZONU with zonal wind and array ZONT
! with temperature, both at 74 latitudes and 22 pressure levels, and fills
! array ZONP with the pressure values at those 22 levels.

  call fill_zonavg(zonclim,idate,imonth)

! In order to obtain zonally averaged water vapor mixing ratio, interpolate
! Mclatchy soundings in time to mclat array and obtain coefficients for
! latitudinal spline interpolation

  call fill_zonr_mclat(idate,imonth,iyear)

! Compute zonal height field from hydrostatic and geostrophic balances

  zonz(37,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)
  zonz(38,1) = 113. ! Std hgt of 1000 mb sfc in tropics (from Mclatchy sounding)

! Use geostrophic balance eqn to get 1000 mb heights at all latitudes

  dyo2g = 2.5 * dlat / grav2  ! spacing (m) between zonavg values divided by 2g

  do ilatn = 39,74
     ilats = 75 - ilatn 

! Coriolis parameter at midpoint between zonavg values (north and south hemispheres)
   
     fcorn = 2. * omega * sin(2.5 * (ilatn-38) * pio180)
     fcors = 2. * omega * sin(2.5 * (ilats-37) * pio180)

! Compute and apply delta zp using midpoint average of two zonu values

     zonz(ilatn,1) = zonz(ilatn-1,1) &
        - fcorn * dyo2g * (zonu(ilatn-1,1) + zonu(ilatn,1))
     zonz(ilats,1) = zonz(ilats+1,1) &
        - fcorn * dyo2g * (zonu(ilats+1,1) + zonu(ilats,1))
  enddo

! Use hydrostatic equation to get all heights above 1000 mb surface

  cpo2g = cp / grav2

  do ilat = 1,74
     do iplev = 2,22
        pilo = (1.e-5 * zonp_vect(iplev-1)) ** rocp
        pihi = (1.e-5 * zonp_vect(iplev)) ** rocp
        thetavlo = zont(ilat,iplev-1) * (1. + eps_virt * zonr(ilat,iplev-1)) &
                 / pilo
        thetavhi = zont(ilat,iplev) * (1. + eps_virt * zonr(ilat,iplev)) &
                 / pihi
        zonz(ilat,iplev) = zonz(ilat,iplev-1) + cpo2g * (pilo - pihi) &
           * (thetavlo + thetavhi)
     enddo
  enddo

  return
  end subroutine zonavg_init

!===============================================================================

  subroutine fill_zonavg(zonclim,idate,imonth)

  use misc_coms, only: io6

  implicit none

  character(len=*), intent(in) :: zonclim
  integer, intent(in) :: idate,imonth

  integer :: ilat,iplev,im,im1,im2
  integer, save :: iread = 0
  real :: wt1,wt2
  logical :: fexists
  character(len=1) :: line
  integer, external :: julday

! Read ZONCLIM dataset if it has not been read yet

  if (iread /= 1) then
     iread = 1
   
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
           READ(29,*) (zont12(ilat+1,iplev,im),ilat=1,72)
        enddo
     enddo

     do im = 1,12
        read(29,*) line
        do iplev = 22,1,-1
           READ(29,*) (zonu12(ilat+1,iplev,im),ilat=1,72)
        enddo
     enddo

     close(29)
  endif

! Array time indices and interpolation weights for bracketing months

! (Weight computation assumes that all months have 31 days, which results in
! modest jumps in value from the end of shorter months to the beginning of
! new months.  However, the method ensures that values are always interpolated
! between consecutive months, and never extrapolated.)

  if (idate < 16) then
     im1 = imonth - 1
     im2 = imonth
  else
     im1 = imonth
     im2 = imonth + 1
  endif
  
  if (im1 == 0) im1 = 12
  if (im2 == 13) im2 = 1

  wt2 = .033333 * mod(idate+15,31)
  wt1 = 1. - wt2

! Loop over all pressure levels and latitudes

  do iplev = 1,22
     do ilat = 2,73

! Interpolate zonal T and u in time

        zont(ilat,iplev) = wt1 * zont12(ilat,iplev,im1) &
                         + wt2 * zont12(ilat,iplev,im2)
        zonu(ilat,iplev) = wt1 * zonu12(ilat,iplev,im1) &
                         + wt2 * zonu12(ilat,iplev,im2)
     enddo

! Fill in N/S boundary values

     zont(1 ,iplev) =  zont(2 ,iplev)
     zont(74,iplev) =  zont(73,iplev)
     zonu(1 ,iplev) = -zonu(2 ,iplev)
     zonu(74,iplev) = -zonu(73,iplev)

! Fill pressure column array

     zonp_vect(iplev) = 10. ** (float(31-iplev)/6.)
  enddo

  return
  end subroutine fill_zonavg

!===============================================================================

  subroutine fill_zonr_mclat(idate,imonth,iyear)

  use misc_coms,   only: io6
  use mem_mclat,   only: slat, mclat, ypp_mclat, mcol, mclat_spline

  implicit none   

  integer, intent(in) :: idate,imonth,iyear

  integer :: ilat,lv,jday
  real :: alat
  integer, external :: julday

! Subroutine fill_zonr_mclat fills the ZONR array with specific humidity values
! by means of interpolation from the McLatchy sounding.  The ZONR array is 2-D,
! defined at 22 vertical pressure levels and 74 latitudes (2.5 degree spacing),
! and is required to supplement the zonu and zont arrays which are filled with
! monthly climatological data from the U.K. Met office.

! Get current Julian day

  jday = julday(imonth,idate,iyear)

! Compute spline coefficients in preparation for interpolation

  call mclat_spline(jday)  

! Loop over latitude columns of zonavg data

  do ilat = 1,74
     alat = -93.75 + 2.5 * ilat

! Loop over vertical levels in McLatchy sounding

     do lv = 1,33

! Spline-interpolate (by latitude) Mclatchy pressure, vapor density, and
! total density to current latitude

        call spline2(13,slat,mclat(1,lv,2),ypp_mclat(1,lv,2),alat,mcol(lv,2))
        call spline2(13,slat,mclat(1,lv,4),ypp_mclat(1,lv,4),alat,mcol(lv,4))
        call spline2(13,slat,mclat(1,lv,6),ypp_mclat(1,lv,6),alat,mcol(lv,6))

! Compute specific humidity from quotient of vapor and total density

        mcol(lv,4) = mcol(lv,4) / mcol(lv,6)
     enddo

! Vertically interpolate Mclatchy vapor specific humidity BY PRESSURE to zonavg
! data levels

     call pintrp_ee(33,mcol(1,4),mcol(1,2),22,zonr_vect,zonp_vect)
     zonr(ilat,1:22) = zonr_vect(1:22)
  enddo

  return
  end subroutine fill_zonr_mclat

End Module mem_zonavg


