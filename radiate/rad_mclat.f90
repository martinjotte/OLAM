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
subroutine rad_mclat(iw,nrad,koff,glat,dl,pl,rl,tl,o3l,zml,ztl,dzl,dpl)

use mem_mclat,   only: sslat, mclat, ypp_mclat, mcol
use mem_radiate, only: nadd_rad, zmrad
use consts_coms, only: gordry, grav, rdry
use mem_grid,    only: mza, zm
use misc_coms,   only: io6

implicit none

integer, intent(in) :: iw       ! current column index
integer, intent(in) :: nrad     ! # of vertical radiation levels
integer, intent(in) :: koff     ! offset between model and radiation levels

real, intent(in) :: glat

real, intent(inout) :: dl (nrad)
real, intent(inout) :: pl (nrad)
real, intent(inout) :: rl (nrad)
real, intent(inout) :: tl (nrad)
real, intent(inout) :: o3l(nrad)
real, intent(inout) :: zml(nrad)
real, intent(inout) :: ztl(nrad)
real, intent(inout) :: dzl(nrad)
real, intent(inout) :: dpl(nrad)

integer :: lv, lf, k, kadd
real    :: deltaz, tavg, ptop

! Before this subroutine is called, McClatchey soundings have been interpolated
! in time between summer and winter values, and spline coefficients for latitudinal
! interpolation have been pre-computed.
!
! The following call to spline2 completes the latitudinal spline interpolation 
! of a complete vertical column for a specific latitude.
! The result is in mcol(lv,lf).

do lf = 1,6  ! Loop over number of data types

   ! no need for pressure
   if (lf == 2) cycle                      

   ! only need o3 and z if nadd_rad = 0 or 1
   if (nadd_rad <= 1 .and. lf /= 5 .and. lf /= 1) cycle

   do lv = 1,33    ! Loop over number of vertical levels

      call spline2(13,sslat,mclat(1,lv,lf),ypp_mclat(1,lv,lf),glat,mcol(lv,lf))

   enddo
enddo

! Model values of dl, pl, tl, rl, zml, and ztl were filled in radiation driver
! from k = 1 to k = mza - 1 - koff

if (nadd_rad > 1) then

   if (zmrad < zm(mza-1)) then
      write(io6,*) 'Error - top of radiation grid is below the model grid'
      stop    'in rad_mclat'
   endif

   ! Compute heights of added levels for this column, excluding the
   ! top isothermal layer that extends to TOA

   deltaz = (zmrad - zm(mza-1)) / real(nadd_rad-1)

   do k = mza-koff, nrad-1
      zml(k) = zml(k-1) + deltaz
      ztl(k) = .5 * (zml(k) + zml(k-1))
   enddo

   ! Interpolate temperature and vapor density to the added levels,
   ! excluding the top isothermal layer

   kadd = mza - koff
   call hintrp_cc(33, mcol(1,3), mcol(1,1), nadd_rad-1, tl(kadd), ztl(kadd))
   call hintrp_cc(33, mcol(1,4), mcol(1,1), nadd_rad-1, rl(kadd), ztl(kadd))
   call hintrp_cc(33, mcol(1,6), mcol(1,1), nadd_rad-1, dl(kadd), ztl(kadd))

   ! Compute pressure of added levels by hydrostatic integration.
   
   do k = kadd, nrad-1
      tavg = 0.5 * (tl(k)+tl(k-1))
      pl(k) = pl(k-1) * exp( -gordry * (ztl(k) - ztl(k-1)) / tavg )
   enddo

endif

if (nadd_rad > 0) then

   ! Define profile values for extra layer from model top to top of atmosphere
   
   ptop = pl(nrad-1) * exp( -gordry * (zml(nrad-1) - ztl(nrad-1)) / tl(nrad-1) )

   tl(nrad) = tl(nrad-1)
   pl(nrad) = 0.5 * ptop
   dl(nrad) = pl(nrad) / (rdry * tl(nrad))
   rl(nrad) = rl(nrad-1) * dl(nrad) / dl(nrad-1)

   zml(nrad) = zml(nrad-1) + ptop / (dl(nrad) * grav)
   ztl(nrad) = 0.5 * (zml(nrad) + zml(nrad-1))

endif   

! Compute dzl values.

do k = 2,nrad
   dzl(k) = zml(k) - zml(k-1)
enddo
dzl(1) = zml(2)

! Compute dpl values (just estimate hydrostatically for now)

do k = 2, nrad-1
   dpl(k) = dl(k) * grav * dzl(k)
enddo
dpl(1)    = dpl(2)
dpl(nrad) = ptop

! Interpolate O3 from Mclatchy sounding to all levels in radiation column

call hintrp_cc(33, mcol(1,5), mcol(1,1), nrad, o3l, ztl)

return
end subroutine rad_mclat
