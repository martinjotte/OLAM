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
subroutine cart_hex()

use mem_ijtabs, only: itab_md, itab_ud, itab_wd, mrls, alloc_itabsd, &
                      jtm_grid, jtu_grid, jtv_grid, jtw_grid, &
                      jtm_init, jtu_init, jtv_init, jtw_init, &
                      jtm_prog, jtu_prog, jtv_prog, jtw_prog, &
                      jtm_lbcp, jtu_lbcp, jtv_lbcp, jtw_lbcp, &
                      jtm_wadj, jtu_wadj, jtv_wadj, jtw_wadj, &
                      jtm_wstn, jtu_wstn, jtv_wstn, jtw_wstn, &
                      jtm_vadj, jtu_wall, jtv_wall, jtw_vadj

use mem_grid,   only: nma, nua, nwa, mma, mua, mwa, xem, yem, zem, alloc_xyzem
use misc_coms,  only: io6, nxp, deltax
use oplot_coms, only: op
use mem_para,   only: myrank

implicit none

integer :: i, j, im, iu, iw, ir, irm, irp, jw0, iup, iwp, iskip
integer :: im1, im2, im3, im4, iu1, iu2, iu3, iu4, iu5, iw1, iw2, iw3, iw4

real :: unit_dist, xstart, ystart, xm, ym, rxx, rxy, ryx, ryy

real :: xp1, yp1, xp2, yp2, xq1, yq1, xq2, yq2, psiz
real :: xu, yu, xw1, yw1, xw2, yw2

character(10) :: string

real :: jm1(  nxp+1,  nxp+1,3)
real :: ju1(  nxp+1,  nxp+1,3)
real :: ju2(  nxp+1,  nxp+1,3)
real :: ju3(  nxp+1,  nxp+1,3)
real :: jw1(0:nxp+1,  nxp+1,3)
real :: jw2(  nxp+1,0:nxp+1,3)

! Define triangles, edges, and vertices for 3D cartesian hexagonal domain.
! Radius of prognostic region of domain is nxp, and one boundary row
! is added for application of cyclic lateral boundary conditions.
! Cyclic repeat distance of domain is 2 * nxp.
! THIS CASE APPLIES FOR MDOMAIN = 5.

mrls = 1  ! Default value

! Fill temporary structured arrays with indices for all staggered 
! points in unstructured grid, and add up unstructured grid points.

nma = 1
nua = 1
nwa = 1

do ir = 1,3           ! rhombus counter
   do j = 1,nxp       ! row within rhombus
      do i = 1,nxp+1  ! column within rhombus row

         jm1(i,j,ir) = nma + 1

         ju1(i,j,ir) = nua + 1
         ju2(i,j,ir) = nua + 2
         ju3(i,j,ir) = nua + 3

         jw1(i,j,ir) = nwa + 1
         jw2(i,j,ir) = nwa + 2

         nma = nma + 1
         nua = nua + 3
         nwa = nwa + 2

      enddo

      jw1(0,j,ir) = nwa + 1

      nwa = nwa + 1

   enddo

   do i = 1,nxp+1

      jm1(i,nxp+1,ir) = nma + 1
      ju1(i,nxp+1,ir) = nua + 1
      jw2(i,    0,ir) = nwa + 1

      nma = nma + 1
      nua = nua + 1
      nwa = nwa + 1

   enddo
enddo

jw0 = nwa + 1

nwa = nwa + 1

! Allocate itab and earth-coordinate arrays

call alloc_itabsd(nma,nua,nwa)

call alloc_xyzem(nma)

! Assign imp, iup, iwp members of itab arrays equal to array indices
! by default.  Correct for boundary points later.

do im = 2,nma
   itab_md(im)%imp = im
   itab_md(im)%mrlm = 1
   itab_md(im)%mrlm_orig = 1
   call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_wstn, 0)
enddo

do iu = 2,nua
   itab_ud(iu)%iup = iu
   call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
enddo

do iw = 2,nwa
   itab_wd(iw)%iwp = iw
   itab_wd(iw)%mrlw = 1
   itab_wd(iw)%mrlw_orig = 1
   call wdloopf('f',iw, jtw_grid, jtw_vadj, 0, 0, 0, 0)
enddo

! DELTAX, whose value is defined in OLAMIN, is a measure of the
! horizontal grid spacing.  Here, it is defined as the square root
! of the area of a regular hexagonal grid cell on the plane.

! Based on the above definition of DELTAX, calculate unit_dist,
! which is the distance between the centers of adjacent hexagons, or
! equivalently, the length of an edge on the regular triangle grid.

unit_dist = sqrt(sqrt(4./3.)) * deltax

! Assign minimal neighbor index information from which complete stencil
! can (later) be filled.

! Define xstart and ystart as reference coordinates of rhombus corner within
! rhombus-relative coordinate system.  Coordinates of all triangle M points
! within each rhombus are computed relative to this reference, after
! which rotations are performed for rhombuses 2 and 3.

xstart = -real(nxp+1) * .5 * unit_dist
ystart = -(real(nxp) + 1./3.) * .5 * sqrt(3.) * unit_dist

do ir = 1,3
   irm = ir - 1
   if (ir == 1) irm = 3

   irp = ir + 1
   if (ir == 3) irp = 1

! Rotation matrix elements

   if (ir == 1) then
      rxx = 1.
      rxy = 0.
      ryx = 0.
      ryy = 1.
   elseif (ir == 2) then
      rxx = -.5
      rxy = -.5 * sqrt(3.)
      ryx =  .5 * sqrt(3.)
      ryy = -.5
   else
      rxx = -.5
      rxy =  .5 * sqrt(3.)
      ryx = -.5 * sqrt(3.)
      ryy = -.5
   endif

   do j = 1,nxp
      do i = 1,nxp+1

         im1 = jm1(i,j,ir)

! x,y coordinates in frame of current rhombus

         xm = xstart + (real(i-1) - .5 * real(j-1)) * unit_dist
         ym = ystart + real(j-1) * .5 * sqrt(3.) * unit_dist

! Rotation to universal "earth" frame (on planar surface)

         xem(im1) = rxx * xm + rxy * ym
         yem(im1) = ryx * xm + ryy * ym
         zem(im1) = 0.

         iu1 = ju1(i,j,ir)
         iu2 = ju2(i,j,ir)
         iu3 = ju3(i,j,ir)

         iw1 = jw1(i,j,ir)
         iw2 = jw2(i,j,ir)
         iw3 = jw2(i,j-1,ir)
         iw4 = jw1(i-1,j,ir)

         im3 = jm1(i,j+1,ir)
         iu5 = ju1(i,j+1,ir)

         if (i <= nxp) then
            im2 = jm1(i+1,j,  ir)
            im4 = jm1(i+1,j+1,ir)
            iu4 = ju3(i+1,j  ,ir)
         else
            im2 = jm1(j  ,nxp+1,irp)
            im4 = jm1(j+1,nxp+1,irp)
            iu4 = ju1(j  ,nxp+1,irp)
         endif

         if (ir == 1) then
            itab_ud(iu1)%im(1) = im1
            itab_ud(iu1)%im(2) = im2
            itab_ud(iu1)%iw(1) = iw3
            itab_ud(iu1)%iw(2) = iw1
         else
            itab_ud(iu1)%im(1) = im2
            itab_ud(iu1)%im(2) = im1
            itab_ud(iu1)%iw(1) = iw1
            itab_ud(iu1)%iw(2) = iw3
         endif
         
         if (ir == 1 .or. ir == 3) then
            itab_ud(iu2)%im(1) = im1
            itab_ud(iu2)%im(2) = im4
            itab_ud(iu2)%iw(1) = iw1
            itab_ud(iu2)%iw(2) = iw2
         else
            itab_ud(iu2)%im(1) = im4
            itab_ud(iu2)%im(2) = im1
            itab_ud(iu2)%iw(1) = iw2
            itab_ud(iu2)%iw(2) = iw1
         endif

         if (ir == 3) then
            itab_ud(iu3)%im(1) = im1
            itab_ud(iu3)%im(2) = im3
            itab_ud(iu3)%iw(1) = iw2
            itab_ud(iu3)%iw(2) = iw4
         else
            itab_ud(iu3)%im(1) = im3
            itab_ud(iu3)%im(2) = im1
            itab_ud(iu3)%iw(1) = iw4
            itab_ud(iu3)%iw(2) = iw2
         endif

         itab_wd(iw1)%iu(1) = iu1
         itab_wd(iw1)%iu(2) = iu4
         itab_wd(iw1)%iu(3) = iu2

         itab_wd(iw2)%iu(1) = iu2
         itab_wd(iw2)%iu(2) = iu5
         itab_wd(iw2)%iu(3) = iu3

         if (i == 1 .and. j == 1) then

            itab_md(im1)%imp = jm1(2,2,irp)
            call mdloopf('f',im1, jtm_lbcp, jtm_wadj, jtm_wstn, 0, 0, 0)

            if (ir == 2) then
               itab_ud(iu1)%iup = ju3(2,1,irm)
            else
               itab_ud(iu1)%iup = ju2(1,1,irp)
            endif
            call udloopf('f',iu1, jtu_lbcp, 0, 0, 0, 0, 0)

            if (ir == 3) then
               itab_ud(iu2)%iup = ju3(2,1,irp)
               call udloopf('f',iu2, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
            endif

            itab_ud(iu3)%iup = ju1(2,2,irp)
            call udloopf('f',iu3, jtu_lbcp, 0, 0, 0, 0, 0)

            itab_wd(iw3)%iwp = jw2(2,1,irm)
            itab_wd(iw3)%iu(1) = iu1
            call wdloopf('f',iw3, jtw_lbcp, 0, 0, 0, 0, 0)

            itab_wd(iw4)%iwp = jw1(2,2,irp)
            itab_wd(iw4)%iu(1) = iu3
            call wdloopf('f',iw4, jtw_lbcp, 0, 0, 0, 0, 0)

         elseif (i == 1) then

            itab_md(im1)%imp = jm1(j+1,2,irp)
            call mdloopf('f',im1, jtm_lbcp, jtm_wadj, jtm_wstn, 0, 0, 0)

            if (ir /= 2) then
               itab_ud(iu1)%iup = ju2(j,1,irp)
               call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
            endif

            if (ir == 3) then
               itab_ud(iu2)%iup = ju3(j+1,1,irp)
               call udloopf('f',iu2, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
            endif

            itab_ud(iu3)%iup = ju1(j+1,2,irp)
            call udloopf('f',iu3, jtu_lbcp, 0, 0, 0, 0, 0)

            itab_wd(iw4)%iwp = jw1(j+1,2,irp)
            itab_wd(iw4)%iu(1) = iu3
            call wdloopf('f',iw4, jtw_lbcp, 0, 0, 0, 0, 0)

         elseif (j == 1) then

            itab_md(im1)%imp = jm1(2,i,irm)
            call mdloopf('f',im1, jtm_lbcp, jtm_wadj, jtm_wstn, 0, 0, 0)
            
            call udloopf('f',iu1, jtu_lbcp, 0, 0, 0, 0, 0)

            if (i == nxp+1 .and. ir == 2) then
               itab_ud(iu1)%iup = ju1(1,nxp+1,ir)
            elseif (i == nxp+1) then
               itab_ud(iu1)%iup = ju2(nxp,1,irp)
            else
               itab_ud(iu1)%iup = ju3(2,i,irm)
            endif

            if (ir == 3) then
               itab_ud(iu2)%iup = ju1(1,i,irm)
               call udloopf('f',iu2, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
            endif

            if (ir /= 1) then
               itab_ud(iu3)%iup = ju2(1,i-1,irm)
               call udloopf('f',iu3, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
            endif

            if (i == nxp+1) then
               itab_wd(iw3)%iwp = jw2(nxp+1,1,irp)
            else
               itab_wd(iw3)%iwp = jw2(2,i,irm)
            endif
            itab_wd(iw3)%iu(1) = iu1
            call wdloopf('f',iw3, jtw_lbcp, 0, 0, 0, 0, 0)

         endif

      enddo
   enddo

   do i = 1,nxp+1

      im1 = jm1(i,nxp+1,ir)
      iu1 = ju1(i,nxp+1,ir)
      iw3 = jw2(i,nxp  ,ir)

! x,y coordinates in frame of current rhombus

      xm = xstart + (real(i-1) - .5 * real(nxp)) * unit_dist
      ym = ystart + real(nxp) * .5 * sqrt(3.) * unit_dist

! Rotation to universal "earth" frame (on planar surface)

      xem(im1) = rxx * xm + rxy * ym
      yem(im1) = ryx * xm + ryy * ym
      zem(im1) = 0.

      if (i <= nxp) then
         im2 = jm1(i+1,nxp+1,ir)
         iw1 = jw1(nxp+1,i,irm)
      else
         im2 = jm1(i,nxp+1,irp)
         iw1 = jw0

         itab_wd(iw1)%iu(ir) = iu1
      endif

      if (ir == 1) then
         itab_ud(iu1)%im(1) = im1
         itab_ud(iu1)%im(2) = im2
         itab_ud(iu1)%iw(1) = iw3
         itab_ud(iu1)%iw(2) = iw1
      else
         itab_ud(iu1)%im(1) = im2
         itab_ud(iu1)%im(2) = im1
         itab_ud(iu1)%iw(1) = iw1
         itab_ud(iu1)%iw(2) = iw3
      endif

      if (i == 1) then

         itab_md(im1)%imp = jm1(2,nxp+1,irm)
         call mdloopf('f',im1, jtm_lbcp, jtm_wadj, jtm_wstn, 0, 0, 0)

         if (ir /= 2) then
            itab_ud(iu1)%iup = ju2(nxp+1,1,irp)
            call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
         endif

      endif

   enddo

enddo

call tri_neighbors()

! Plot grid lines

if (.false.) then

   if (myrank /= 0) return

   call o_reopnwk()
   call plotback()

   mma = nma
   mua = nua
   mwa = nwa

   call oplot_set(1)

   psiz = .035 / real(nxp)

   do iu = 2,nua
      im1 = itab_ud(iu)%im(1)
      im2 = itab_ud(iu)%im(2)
      iw1 = itab_ud(iu)%iw(1)
      iw2 = itab_ud(iu)%iw(2)

      call oplot_transform(1,xem(im1),yem(im1),zem(im1),xp1,yp1)
      call oplot_transform(1,xem(im2),yem(im2),zem(im2),xp2,yp2)

      call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

      if (iskip == 1) cycle

      call o_frstpt (xq1,yq1)
      call o_vector (xq2,yq2)

      if ( xp1 < op%xmin .or.  &
           xp1 > op%xmax .or.  &
           yp1 < op%ymin .or.  &
           yp1 > op%ymax ) cycle

      write(string,'(I0)') im1
      call o_plchlq (xp1,yp1,trim(adjustl(string)),psiz,0.,0.)

      write(string,'(I0)') im2
      call o_plchlq (xp2,yp2,trim(adjustl(string)),psiz,0.,0.)

      ! Compute locations of U and W points based on M point locations

      xu = .5 * (xp1 + xp2)
      yu = .5 * (yp1 + yp2)

      xw1 = xu + sqrt(3.)/6. * (yp2 - yp1)
      yw1 = yu - sqrt(3.)/6. * (xp2 - xp1)
      xw2 = xu - sqrt(3.)/6. * (yp2 - yp1)
      yw2 = yu + sqrt(3.)/6. * (xp2 - xp1)

      write(string,'(I0)') iu
      call o_plchlq (xu,yu,trim(adjustl(string)),psiz,0.,0.)

      write(string,'(I0)') iw1
      call o_plchlq (xw1,yw1,trim(adjustl(string)),psiz,0.,0.)

      write(string,'(I0)') iw2
      call o_plchlq (xw2,yw2,trim(adjustl(string)),psiz,0.,0.)
   enddo

   call o_frame()
   call o_clswk()

endif

end subroutine cart_hex

