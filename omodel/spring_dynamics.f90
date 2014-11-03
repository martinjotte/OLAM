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
subroutine spring_dynamics(ngr,nconcave)

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd
use mem_grid,    only: nma, nua, nwa, xem, yem, zem, impent, mrows
use consts_coms, only: pi2, erad, erador5
use misc_coms,   only: io6, nxp, ngrids, mdomain, deltax
use oplot_coms,  only: op

implicit none

integer, intent(in) :: ngr,nconcave

integer :: niter
real, parameter :: relax = .04, beta = 1.2

! Automatic arrays

real :: dist (nua)
real :: dist0(nua)

real :: dx(nua)
real :: dy(nua)
real :: dz(nua)

real :: dirs(nma,7)

integer :: iu,iu1,iu2,iu3,iu4
integer :: im,im1,im2,im3,im4
integer :: iw,iw1,iw2,j
integer :: iter,ipent,mrow1,mrow2,npoly,ngrw,mrmax,mrmin

integer :: iskip
real :: xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2

real :: cosphi3, cosphi4
real :: coslim3, coslim4, ratio

real :: dist00,distm,expansion,frac_change
character(10) :: string

! Array MOVEM is used to flag selected M points that will be allowed to move
! in the spring dynamics procedure.  If MOVEM = 1 the point is allowed to move,
! and if MOVEM = 0 the point remains stationary.  This selective movement is
! optional for static mesh refinements, but it will be a requirement for 
! dynamically adaptive refinements.  

integer :: movem(nma)
integer :: moveall(ngrids)

write(io6,*) "In spring dynamics...."

!--------------------------------------------------------------
! SPECIAL CODE TO BE MODIFIED FOR DYNAMIC NESTS

if (mdomain == 5) then
   moveall(:) = 0
else
   moveall(:) = 1
endif

! END SPECIAL CODE
!--------------------------------------------------------------

! special
!RETURN
! end special

! For the case where nconcave = 3 and moveall(ngr) = 0, reset all movem values
! to 0, and then set those in border zone of NGR back to 1.  Specifically, 
! border zone consists of those M points that are adjacent to M points with
! MROW values of -3, -2, -1, or 1.  

if (ngr > 1 .and. nconcave == 3 .and. moveall(ngr) == 0) then

   movem(:) = 0
   niter    = 750

   do im = 2,nma
      npoly = itab_md(im)%npoly

      do j = 1,npoly
         iw = itab_md(im)%iw(j)

         ngrw = abs(itab_wd(iw)%mrowh) / 100

         if (ngrw /= ngr) cycle
         
         if (itab_wd(iw)%mrow == -3 .or. &
             itab_wd(iw)%mrow == -2 .or. &
             itab_wd(iw)%mrow == -1 .or. &
             itab_wd(iw)%mrow ==  1) then
            movem(im) = 1
            cycle ! once this point has been flagged, we can move on
         endif

      enddo
   enddo      

else

! The default is to set MOVEM flag = 1 for all points

   movem(:) = 1
   niter    = 2000

endif

! Compute mean length of coarse mesh U segments

if (mdomain < 2) then
   dist00 = beta * pi2 * erad / (5. * real(nxp))
else
   dist00 = deltax * sqrt(2.0) / sqrt(sqrt(3.0))
endif

! Compute target length of each triangle U segment
! Loop over all U points

!$omp parallel
!$omp do private (im1,im2,iw1,iw2,mrow1,mrow2,mrmax,mrmin)
do iu = 2, nua

   im1 = itab_ud(iu)%im(1)
   im2 = itab_ud(iu)%im(2)

   if (movem(im1) == 0 .and. movem(im2) == 0) cycle

! Compute target distance for any MRL value

   dist0(iu) = dist00 * 2.0 ** (1 - itab_ud(iu)%mrlu)

! Modified distance in MRL border zone

   if (ngr > 1) then

      iw1 = itab_ud(iu)%iw(1)
      iw2 = itab_ud(iu)%iw(2)

      mrow1 = itab_wd(iw1)%mrow
      mrow2 = itab_wd(iw2)%mrow

      if (nconcave == 3) then
      
         mrmax = max(mrow1,mrow2)
         mrmin = min(mrow1,mrow2)
      
         if     (mrmax == -2 .and. mrmin == -2) then
            dist0(iu) = dist0(iu) *  7. / 6.  !* .90
         elseif (mrmax == -1 .and. mrmin == -2) then
            dist0(iu) = dist0(iu) *  8. / 6.  !* .90
         elseif (mrmax == -1 .and. mrmin == -1) then
            dist0(iu) = dist0(iu) *  9. / 6.  !* .90
         elseif (mrmax == 1 .and. mrmin == -1) then
            dist0(iu) = dist0(iu) * 10. / 6.  !* .90
         elseif (mrmax == 1 .and. mrmin == 1) then
            dist0(iu) = dist0(iu) * 11. / 12. !* .90
         endif

      elseif (mrow1 > 0 .and. mrow1 <= mrows .and.  &
              mrow2 > 0 .and. mrow2 <= mrows) then

         dist0(iu) = dist0(iu) * (.5 + .25/real(mrows) * (mrow1 + mrow2 - 1))      

      endif ! nconcave

   endif ! ngr

enddo
!$omp end do

!$omp do private (npoly,j,iu)
do im = 2, nma

   if (movem(im) == 0) cycle

   npoly = itab_md(im)%npoly
   do j = 1, npoly

      iu = itab_md(im)%iu(j)
      if (itab_ud(iu)%im(2) == im) then
         dirs(im,j) =  relax
      else
         dirs(im,j) = -relax
      endif

   enddo
enddo
!$omp end do

! Main iteration loop 
do iter = 1, niter

! Compute length of each U segment once per iteration
! Loop over all U points

   !$omp do private (iu1,iu3,im1,im2,im3,im4)
   do iu = 2, nua

      iu1 = itab_ud(iu)%iu(1)
      iu3 = itab_ud(iu)%iu(3)

      im1 = itab_ud(iu)%im(1)
      im2 = itab_ud(iu)%im(2)

      if (itab_ud(iu1)%im(1) == im1) then
         im3 = itab_ud(iu1)%im(2)
      else
         im3 = itab_ud(iu1)%im(1)
      endif

      if (itab_ud(iu3)%im(1) == im1) then
         im4 = itab_ud(iu3)%im(2)
      else
         im4 = itab_ud(iu3)%im(1)
      endif

      if ( movem(im1) == 0 .and. movem(im2) == 0 .and. &
           movem(im3) == 0 .and. movem(im4) == 0 ) cycle

      dist(iu) = sqrt( (xem(im1) - xem(im2)) ** 2 &
                     + (yem(im1) - yem(im2)) ** 2 &
                     + (zem(im1) - zem(im2)) ** 2 )
   enddo
   !$omp end do
   
! Loop over all U points

   !$omp do private (im1,im2,iu1,iu2,iu3,iu4,cosphi3,cosphi4, &
   !$omp             ratio,distm,frac_change)
   do iu = 2, nua

      im1 = itab_ud(iu)%im(1)
      im2 = itab_ud(iu)%im(2)

      if (movem(im1) == 0 .and. movem(im2) == 0) cycle

! Adjustment of dist0 based on opposite angles of triangles

      iu1 = itab_ud(iu)%iu(1)
      iu2 = itab_ud(iu)%iu(2)
      iu3 = itab_ud(iu)%iu(3)
      iu4 = itab_ud(iu)%iu(4)

! Compute cosine of angles at IM3 and IM4

      cosphi3 = (dist(iu1)**2 + dist(iu2)**2 - dist(iu)**2) / (2. * dist(iu1) * dist(iu2))
      cosphi4 = (dist(iu3)**2 + dist(iu4)**2 - dist(iu)**2) / (2. * dist(iu3) * dist(iu4))

! Ratio of smaller cosine to limiting value of cos(72 deg)

      ratio = min(cosphi3,cosphi4) / .309
      
! Decrease dist0 if ratio < 1
! (NOTE: could just use .309 limit regardless of npoly)

      if (ratio < 1.) then
         distm = dist0(iu) * ratio
      else
         distm = dist0(iu)
      endif

! Fractional change to dist that would make it equal dist0

      frac_change = (distm - dist(iu)) / dist(iu)

! Compute components of displacement that gives dist0

      dx(iu) = (xem(im2) - xem(im1)) * frac_change
      dy(iu) = (yem(im2) - yem(im1)) * frac_change
      dz(iu) = (zem(im2) - zem(im1)) * frac_change

   enddo
   !$omp end do
   
! Loop over all M points

   !$omp do private(npoly,j,iu,expansion)
   do im = 2, nma
      
      if (movem(im) == 0) cycle

! For preventing either polar M point from moving:
!      if (im == impent(1 )) cycle
!      if (im == impent(12)) cycle

! For preventing all pentagonal points from moving:
!     if (any(im == impent(1:12)) cycle

! Apply the displacement components to each M point

      npoly = itab_md(im)%npoly
      do j = 1, npoly
         iu = itab_md(im)%iu(j)
         xem(im) = xem(im) + dirs(im,j) * dx(iu)
         yem(im) = yem(im) + dirs(im,j) * dy(iu)
         zem(im) = zem(im) + dirs(im,j) * dz(iu)
      enddo

! Push M point coordinates out to earth radius
      
      if (mdomain < 2) then
         expansion = erad / sqrt(xem(im) ** 2 + yem(im) ** 2 + zem(im) ** 2)
         xem(im) = xem(im) * expansion
         yem(im) = yem(im) * expansion
         zem(im) = zem(im) * expansion
      endif

   enddo
   !$omp end do

! Section for plotting grid at intermediate stages of spring dynamics adjustment

!plt   if (ngr > 1 .and. mod(iter,20) == 1) then
!plt
!plt! Plot grid lines
!plt
!plt      call o_reopnwk()
!plt      call plotback()
!plt
!plt      call oplot_set(1)
!plt 
!plt      do iu = 2,nua
!plt         im1 = itab_ud(iu)%im(1)
!plt         im2 = itab_ud(iu)%im(2)
!plt
!plt         call oplot_transform(1,xem(im1),yem(im1),zem(im1),xp1,yp1)
!plt         call oplot_transform(1,xem(im2),yem(im2),zem(im2),xp2,yp2)
!plt
!plt         call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)
!plt
!plt         if (iskip == 1) cycle
!plt
!plt         call o_frstpt (xq1,yq1)
!plt         call o_vector (xq2,yq2)
!plt      enddo
!plt
!plt! print mrow values
!plt
!plt      do iw = 2, nwa
!plt         im1 = itab_wd(iw)%im(1)
!plt         im2 = itab_wd(iw)%im(2)
!plt         im3 = itab_wd(iw)%im(3)
!plt
!plt         call oplot_transform(1, (xem(im1)+xem(im2)+xem(im3))/3., &
!plt                                 (yem(im1)+yem(im2)+yem(im3))/3., &
!plt                                 (zem(im1)+zem(im2)+zem(im3))/3., &
!plt                                 xp1, yp1                         )
!plt
!plt         if ( xp1 < op%xmin .or.  &
!plt              xp1 > op%xmax .or.  &
!plt              yp1 < op%ymin .or.  &
!plt              yp1 > op%ymax ) cycle
!plt
!plt         write(string,'(I0)') itab_wd(iw)%mrow
!plt         call o_plchlq (xp1,yp1,trim(adjustl(string)),0.01,0.,0.)
!plt      enddo
!plt
!plt      call o_frame()
!plt      call o_clswk()
!plt
!plt   endif ! mod(iter,*)

enddo ! iter

!$omp end parallel

return
end subroutine spring_dynamics
