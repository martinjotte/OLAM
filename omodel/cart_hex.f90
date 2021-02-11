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

  integer :: jm1(  nxp+1,  nxp+1,3)
  integer :: ju1(  nxp+1,  nxp+1,3)
  integer :: ju2(  nxp+1,  nxp+1,3)
  integer :: ju3(  nxp+1,  nxp+1,3)
  integer :: jw1(0:nxp+1,  nxp+1,3)
  integer :: jw2(  nxp+1,0:nxp+1,3)

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
     itab_md(im)%ngr = 1
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
     itab_wd(iw)%ngr = 1
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

  call tri_neighbors(nma, nua, nwa, itab_md, itab_ud, itab_wd)

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

!===============================================================================

subroutine cart4_hex()

  use mem_ijtabs, only: itab_m, itab_v, itab_w, mrls, alloc_itabs, &
                        jtm_grid, jtv_grid, jtw_grid, &
                        jtm_init, jtv_init, jtw_init, &
                        jtm_prog, jtv_prog, jtw_prog, &
                        jtm_lbcp, jtv_lbcp, jtw_lbcp, &
                        jtm_wadj, jtv_wadj, jtw_wadj, &
                        jtm_wstn, jtv_wstn, jtw_wstn, &
                        jtm_vadj, jtv_wall, jtw_vadj

  use mem_grid,   only: nma, nua, nva, nwa, mma, mua, mva, mwa, &
                        xem, yem, zem, xew, yew, zew, &
                        alloc_xyzem, alloc_xyzew

  use misc_coms,  only: io6, nxp, deltax

  implicit none

  integer :: i, j, iw, im, iv
  integer :: jw(11), jaw(11), jzw(11)
  integer :: jm(29), jam(29), jzm(29)
  integer :: jv(28), jav(28), jzv(28)

  real :: centx, unit_dist

  character(10) :: string

  ! Constructs hexagon grid directly (not beginning from triangle grid), and
  ! therefore is not succeeded by calls to tri_neighbors, voronoi and pcvt.
  ! However, grid_geometry_hex IS called after this subroutine.

  ! Channel width = 4 allocated IW rows with cyclic repeat distance 2 rows wide
  ! THIS CASE APPLIES FOR MDOMAIN = 4 (cyclic x bnds & cyclic y bnds)

  ! Channel length is NXP, and cyclic repeat distance is NXP * UNIT_DIST
  ! (unit_dist defined below).  NXP MUST BE AT LEAST 2.

  if (nxp < 2) then
     print*, 'nxp must be at least 2 when mdomain = 4. '
     stop 'stop1 cart4_hex'
  endif

  mrls = 1  ! Default value

    ! Use nxp to count hexagons

  nwa =  4 * nxp + 7  ! Includes index=1 point
  nma =  8 * nxp + 11 ! Includes index=1 point
  nva = 10 * nxp + 10 ! Includes index=1 point
  nua = nva

  mma = nma
  mwa = nwa
  mva = nva
  mua = nua

! Allocate itab and earth-coordinate arrays

  call alloc_itabs(nma,nva,nwa,0)

  call alloc_xyzew(nwa)
  call alloc_xyzem(nma)

  zew(:) = 0.
  zem(:) = 0.

  do iw = 2,nwa
     itab_w(iw)%iwglobe   = iw
     itab_w(iw)%mrlw      = 1
     itab_w(iw)%mrlw_orig = 1
     itab_w(iw)%ngr = 1
  enddo

  do im = 2,nma
     itab_m(im)%imglobe   = im
     itab_m(im)%mrlm      = 1
     itab_m(im)%mrlm_orig = 1
     itab_m(im)%ngr = 1
  enddo

  do iv = 2,nva
     itab_v(iv)%ivglobe   = iv
     itab_v(iv)%mrlv      = 1
  enddo

  ! DELTAX, whose value is defined in OLAMIN, is a measure of the
  ! horizontal grid spacing.  Here, it is defined as the square root
  ! of the area of a regular hexagonal grid cell on the plane.

  ! Based on the above definition of DELTAX, calculate unit_dist,
  ! which is the distance between the centers of adjacent hexagons, or
  ! equivalently, the length of an edge on the regular triangle grid.

  unit_dist = sqrt(sqrt(4./3.)) * deltax  ! approx = 1.07457 * deltax

  ! Loop along channel length, covering a full transverse row at each step

  do i = 1,nxp

     ! Fill grid point index values in stencil for given i value

     ! JW stencil points used: 2,3,4,5,6,7,8,9,10,11
     ! JM stencil points used: 2,4,7,9,10,11,12,13,14,15,16,18,19,20,21,22,23,29
     ! JV stencil points used: 2,5,9,11,12,13,14,15,16,17,18,19,20,23,24,25,26,27,28

     do j = 1,11
        jw(j) =  4 * i + j -  4
     enddo

     do j = 1,29
        jm(j) =  8 * i + j - 12
     enddo

     do j = 1,28
        jv(j) = 10 * i + j - 16
     enddo

     ! Set unused stencil indices to 1

     jw(1) = 1

     jm(1 :5 :2) = 1
     jm(6 :8 :2) = 1
     jm(17     ) = 1
     jm(24:28  ) = 1

     jv( 1   ) = 1
     jv( 3: 4) = 1
     jv( 6: 8) = 1
     jv(10   ) = 1
     jv(21:22) = 1

     ! Modify M and V indices at ends of channel (no modification needed for W)

     if (i == 1) then

        jm(2) = 2
        jm(4) = 3
        jm(7) = 4

        jv(2) = 2
        jv(5) = 3
        jv(9) = 4

        jaw(:) = jw(:)
        jam(:) = jm(:)
        jav(:) = jv(:)

     elseif (i == nxp) then

        jm(18:23) = jm(18:23) - 1
        jm(29   ) = jm(29   ) - 6

        jv(23:28) = jv(23:28) - 2

        jzw(:) = jw(:)
        jzm(:) = jm(:)
        jzv(:) = jv(:)

     endif

     centx = (2 * i - nxp - 1) * 0.5 * unit_dist

     xew(jw( 3     )) = centx - 1.25 * unit_dist
     xew(jw( 2: 4:2)) = centx - 0.75 * unit_dist
     xew(jw( 5: 7:2)) = centx - 0.25 * unit_dist
     xew(jw( 6: 8:2)) = centx + 0.25 * unit_dist
     xew(jw( 9:11:2)) = centx + 0.75 * unit_dist
     xew(jw(10     )) = centx + 1.25 * unit_dist

     yew(jw(5 : 9:4)) =  2.25 * unit_dist / sqrt(3.)
     yew(jw(2 :10:4)) =  0.75 * unit_dist / sqrt(3.)
     yew(jw(3 :11:4)) = -0.75 * unit_dist / sqrt(3.)
     yew(jw(4 : 8:4)) = -2.25 * unit_dist / sqrt(3.)

     xem(jm( 4: 7: 3)) = centx - 1.25 * unit_dist
     xem(jm( 2      )) = centx - 0.75 * unit_dist
     xem(jm(13:14   )) = centx - 0.75 * unit_dist
     xem(jm(11:12   )) = centx - 0.25 * unit_dist
     xem(jm(15:16   )) = centx - 0.25 * unit_dist
     xem(jm( 9:10   )) = centx + 0.25 * unit_dist
     xem(jm(21:22   )) = centx + 0.25 * unit_dist
     xem(jm(19:20   )) = centx + 0.75 * unit_dist
     xem(jm(23      )) = centx + 0.75 * unit_dist
     xem(jm(18:29:11)) = centx + 1.25 * unit_dist

     yem(jm( 9     )) =  2.75 * unit_dist / sqrt(3.) 
     yem(jm( 2:18:8)) =  1.75 * unit_dist / sqrt(3.) 
     yem(jm(11:19:8)) =  1.25 * unit_dist / sqrt(3.) 
     yem(jm( 4:20:8)) =  0.25 * unit_dist / sqrt(3.) 
     yem(jm(13:29:8)) = -0.25 * unit_dist / sqrt(3.) 
     yem(jm(14:22:8)) = -1.25 * unit_dist / sqrt(3.) 
     yem(jm( 7:23:8)) = -1.75 * unit_dist / sqrt(3.)
     yem(jm(16     )) = -2.75 * unit_dist / sqrt(3.) 

     call wloopf('f',jw(5), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     call wloopf('f',jw(6), jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     call wloopf('f',jw(7), jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     call wloopf('f',jw(8), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)

     itab_w(jw(5))%iwp   = jw(7)

     itab_w(jw(6))%npoly = 6
     itab_w(jw(6))%iwp   = jw(6)

     itab_w(jw(6))%iw(1) = jw(11)
     itab_w(jw(6))%iw(2) = jw(10)
     itab_w(jw(6))%iw(3) = jw(9)
     itab_w(jw(6))%iw(4) = jw(5)
     itab_w(jw(6))%iw(5) = jw(2)
     itab_w(jw(6))%iw(6) = jw(7)

     itab_w(jw(6))%im(1) = jm(20)
     itab_w(jw(6))%im(2) = jm(19)
     itab_w(jw(6))%im(3) = jm(10)
     itab_w(jw(6))%im(4) = jm(11)
     itab_w(jw(6))%im(5) = jm(12)
     itab_w(jw(6))%im(6) = jm(21)

     itab_w(jw(6))%iv(1) = jv(26)
     itab_w(jw(6))%iv(2) = jv(24)
     itab_w(jw(6))%iv(3) = jv(12)
     itab_w(jw(6))%iv(4) = jv(13)
     itab_w(jw(6))%iv(5) = jv(14)
     itab_w(jw(6))%iv(6) = jv(15)

     itab_w(jw(6))%dirv(1:3) = -1.0
     itab_w(jw(6))%dirv(4:6) =  1.0

     itab_w(jw(7))%npoly = 6
     itab_w(jw(7))%iwp   = jw(7)

     itab_w(jw(7))%iw(1) = jw(8)
     itab_w(jw(7))%iw(2) = jw(11)
     itab_w(jw(7))%iw(3) = jw(6)
     itab_w(jw(7))%iw(4) = jw(2)
     itab_w(jw(7))%iw(5) = jw(3)
     itab_w(jw(7))%iw(6) = jw(4)

     itab_w(jw(7))%im(1) = jm(22)
     itab_w(jw(7))%im(2) = jm(21)
     itab_w(jw(7))%im(3) = jm(12)
     itab_w(jw(7))%im(4) = jm(13)
     itab_w(jw(7))%im(5) = jm(14)
     itab_w(jw(7))%im(6) = jm(15)

     itab_w(jw(7))%iv(1) = jv(19)
     itab_w(jw(7))%iv(2) = jv(27)
     itab_w(jw(7))%iv(3) = jv(15)
     itab_w(jw(7))%iv(4) = jv(16)
     itab_w(jw(7))%iv(5) = jv(17)
     itab_w(jw(7))%iv(6) = jv(18)

     itab_w(jw(7))%dirv(1:3) = -1.0
     itab_w(jw(7))%dirv(4:6) =  1.0

     itab_w(jw(8))%iwp = jw(6)

     call mloopf('f',jm(10), jtm_lbcp, jtm_grid, 0, 0, 0, 0)
     call mloopf('f',jm(11), jtm_grid, jtm_vadj, 0, 0, 0, 0)
     call mloopf('f',jm(12), jtm_grid, jtm_vadj, 0, 0, 0, 0)
     call mloopf('f',jm(13), jtm_grid, jtm_vadj, 0, 0, 0, 0)
     call mloopf('f',jm(14), jtm_grid, jtm_vadj, 0, 0, 0, 0)
     call mloopf('f',jm(15), jtm_grid, jtm_vadj, 0, 0, 0, 0)
     call mloopf('f',jm(16), jtm_lbcp, 0, 0, 0, 0, 0)

     itab_m(jm(9))%npoly = 1
     itab_m(jm(9))%imp   = jm(21) ! corrected later for jmz row
     itab_m(jm(9))%iv(1) = jv(11)

     itab_m(jm(10))%npoly = 3
     itab_m(jm(10))%imp   = jm(22) ! corrected later for jmz row
     itab_m(jm(10))%iw(1) = jw(6)
     itab_m(jm(10))%iw(2) = jw(9)
     itab_m(jm(10))%iw(3) = jw(5)
     itab_m(jm(10))%iv(1) = jv(11)
     itab_m(jm(10))%iv(2) = jv(13)
     itab_m(jm(10))%iv(3) = jv(12)

     itab_m(jm(11))%npoly = 3
     itab_m(jm(11))%imp   = jm(15)
     itab_m(jm(11))%iw(1) = jw(5)
     itab_m(jm(11))%iw(2) = jw(2)
     itab_m(jm(11))%iw(3) = jw(6)
     itab_m(jm(11))%iv(1) = jv(14)
     itab_m(jm(11))%iv(2) = jv(13)
     itab_m(jm(11))%iv(3) = jv(2)

     itab_m(jm(12))%npoly = 3
     itab_m(jm(12))%imp   = jm(12)
     itab_m(jm(12))%iw(1) = jw(7)
     itab_m(jm(12))%iw(2) = jw(6)
     itab_m(jm(12))%iw(3) = jw(2)
     itab_m(jm(12))%iv(1) = jv(14)
     itab_m(jm(12))%iv(2) = jv(16)
     itab_m(jm(12))%iv(3) = jv(15)

     itab_m(jm(13))%npoly = 3
     itab_m(jm(13))%imp   = jm(13)
     itab_m(jm(13))%iw(1) = jw(2)
     itab_m(jm(13))%iw(2) = jw(3)
     itab_m(jm(13))%iw(3) = jw(7)
     itab_m(jm(13))%iv(1) = jv(17)
     itab_m(jm(13))%iv(2) = jv(16)
     itab_m(jm(13))%iv(3) = jv(5)

     itab_m(jm(14))%npoly = 3
     itab_m(jm(14))%imp   = jm(14)
     itab_m(jm(14))%iw(1) = jw(4)
     itab_m(jm(14))%iw(2) = jw(7)
     itab_m(jm(14))%iw(3) = jw(3)
     itab_m(jm(14))%iv(1) = jv(17)
     itab_m(jm(14))%iv(2) = jv(9)
     itab_m(jm(14))%iv(3) = jv(18)

     itab_m(jm(15))%npoly = 3
     itab_m(jm(15))%imp   = jm(15)
     itab_m(jm(15))%iw(1) = jw(7)
     itab_m(jm(15))%iw(2) = jw(4)
     itab_m(jm(15))%iw(3) = jw(8)
     itab_m(jm(15))%iv(1) = jv(20)
     itab_m(jm(15))%iv(2) = jv(19)
     itab_m(jm(15))%iv(3) = jv(18)

     itab_m(jm(16))%npoly = 1
     itab_m(jm(16))%imp   = jm(12)
     itab_m(jm(16))%iv(1) = jv(20)

     call vloopf('f',jv(11), jtv_lbcp, 0, 0, 0, 0, 0)
     call vloopf('f',jv(12), jtv_grid, jtv_lbcp, jtv_wadj, jtv_wstn, 0, 0)
     call vloopf('f',jv(13), jtv_grid, jtv_lbcp, jtv_wadj, jtv_wstn, 0, 0)
     call vloopf('f',jv(14), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(15), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(16), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(17), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(18), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(19), jtv_grid, jtv_init, jtv_prog, jtv_wadj, jtv_wstn, 0)
     call vloopf('f',jv(20), jtv_lbcp, 0, 0, 0, 0, 0)

     itab_v(jv(11))%ivp   = jv(27) ! corrected later for jzv row
     itab_v(jv(11))%iw(1) = jw(5)
     itab_v(jv(11))%iw(2) = jw(9)
     itab_v(jv(11))%im(1) = jm(10)
     itab_v(jv(11))%im(2) = jm(9)

     itab_v(jv(12))%ivp   = jv(28) ! corrected later for jzv row
     itab_v(jv(12))%iw(1) = jw(6)
     itab_v(jv(12))%iw(2) = jw(9)
     itab_v(jv(12))%iw(3) = jw(10)
     itab_v(jv(12))%iw(4) = jw(5)
     itab_v(jv(12))%im(1) = jm(19)
     itab_v(jv(12))%im(2) = jm(10)
     itab_v(jv(12))%im(3) = jm(20)
     itab_v(jv(12))%im(4) = jm(18)
     itab_v(jv(12))%im(5) = jm(11)
     itab_v(jv(12))%im(6) = jm(9)
     itab_v(jv(12))%iv(1) = jv(24)
     itab_v(jv(12))%iv(2) = jv(23)
     itab_v(jv(12))%iv(3) = jv(13)
     itab_v(jv(12))%iv(4) = jv(11)

     itab_v(jv(13))%ivp   = jv(19)
     itab_v(jv(13))%iw(1) = jw(5)
     itab_v(jv(13))%iw(2) = jw(6)
     itab_v(jv(13))%iw(3) = jw(2)
     itab_v(jv(13))%iw(4) = jw(9)
     itab_v(jv(13))%im(1) = jm(11)
     itab_v(jv(13))%im(2) = jm(10)
     itab_v(jv(13))%im(3) = jm(2)
     itab_v(jv(13))%im(4) = jm(12)
     itab_v(jv(13))%im(5) = jm(9)
     itab_v(jv(13))%im(6) = jm(19)
     itab_v(jv(13))%iv(1) = jv(2)
     itab_v(jv(13))%iv(2) = jv(14)
     itab_v(jv(13))%iv(3) = jv(11)
     itab_v(jv(13))%iv(4) = jv(12)

     itab_v(jv(14))%ivp   = jv(14)
     itab_v(jv(14))%iw(1) = jw(2)
     itab_v(jv(14))%iw(2) = jw(6)
     itab_v(jv(14))%iw(3) = jw(7)
     itab_v(jv(14))%iw(4) = jw(5)
     itab_v(jv(14))%im(1) = jm(12)
     itab_v(jv(14))%im(2) = jm(11)
     itab_v(jv(14))%im(3) = jm(13)
     itab_v(jv(14))%im(4) = jm(21)
     itab_v(jv(14))%im(5) = jm(2)
     itab_v(jv(14))%im(6) = jm(10)
     itab_v(jv(14))%iv(1) = jv(16)
     itab_v(jv(14))%iv(2) = jv(15)
     itab_v(jv(14))%iv(3) = jv(2)
     itab_v(jv(14))%iv(4) = jv(13)

     itab_v(jv(15))%ivp   = jv(15)
     itab_v(jv(15))%iw(1) = jw(7)
     itab_v(jv(15))%iw(2) = jw(6)
     itab_v(jv(15))%iw(3) = jw(11)
     itab_v(jv(15))%iw(4) = jw(2)
     itab_v(jv(15))%im(1) = jm(21)
     itab_v(jv(15))%im(2) = jm(12)
     itab_v(jv(15))%im(3) = jm(22)
     itab_v(jv(15))%im(4) = jm(20)
     itab_v(jv(15))%im(5) = jm(13)
     itab_v(jv(15))%im(6) = jm(11)
     itab_v(jv(15))%iv(1) = jv(27)
     itab_v(jv(15))%iv(2) = jv(26)
     itab_v(jv(15))%iv(3) = jv(16)
     itab_v(jv(15))%iv(4) = jv(14)

     itab_v(jv(16))%ivp   = jv(16)
     itab_v(jv(16))%iw(1) = jw(2)
     itab_v(jv(16))%iw(2) = jw(7)
     itab_v(jv(16))%iw(3) = jw(3)
     itab_v(jv(16))%iw(4) = jw(6)
     itab_v(jv(16))%im(1) = jm(13)
     itab_v(jv(16))%im(2) = jm(12)
     itab_v(jv(16))%im(3) = jm(4)
     itab_v(jv(16))%im(4) = jm(14)
     itab_v(jv(16))%im(5) = jm(11)
     itab_v(jv(16))%im(6) = jm(21)
     itab_v(jv(16))%iv(1) = jv(5)
     itab_v(jv(16))%iv(2) = jv(17)
     itab_v(jv(16))%iv(3) = jv(14)
     itab_v(jv(16))%iv(4) = jv(15)

     itab_v(jv(17))%ivp   = jv(17)
     itab_v(jv(17))%iw(1) = jw(3)
     itab_v(jv(17))%iw(2) = jw(7)
     itab_v(jv(17))%iw(3) = jw(4)
     itab_v(jv(17))%iw(4) = jw(2)
     itab_v(jv(17))%im(1) = jm(14)
     itab_v(jv(17))%im(2) = jm(13)
     itab_v(jv(17))%im(3) = jm(7)
     itab_v(jv(17))%im(4) = jm(15)
     itab_v(jv(17))%im(5) = jm(4)
     itab_v(jv(17))%im(6) = jm(12)
     itab_v(jv(17))%iv(1) = jv(9)
     itab_v(jv(17))%iv(2) = jv(18)
     itab_v(jv(17))%iv(3) = jv(5)
     itab_v(jv(17))%iv(4) = jv(16)

     itab_v(jv(18))%ivp   = jv(18)
     itab_v(jv(18))%iw(1) = jw(4)
     itab_v(jv(18))%iw(2) = jw(7)
     itab_v(jv(18))%iw(3) = jw(8)
     itab_v(jv(18))%iw(4) = jw(3)
     itab_v(jv(18))%im(1) = jm(15)
     itab_v(jv(18))%im(2) = jm(14)
     itab_v(jv(18))%im(3) = jm(16)
     itab_v(jv(18))%im(4) = jm(22)
     itab_v(jv(18))%im(5) = jm(7)
     itab_v(jv(18))%im(6) = jm(13)
     itab_v(jv(18))%iv(1) = jv(20)
     itab_v(jv(18))%iv(2) = jv(19)
     itab_v(jv(18))%iv(3) = jv(9)
     itab_v(jv(18))%iv(4) = jv(17)

     itab_v(jv(19))%ivp   = jv(19)
     itab_v(jv(19))%iw(1) = jw(7)
     itab_v(jv(19))%iw(2) = jw(8)
     itab_v(jv(19))%iw(3) = jw(4)
     itab_v(jv(19))%iw(4) = jw(11)
     itab_v(jv(19))%im(1) = jm(15)
     itab_v(jv(19))%im(2) = jm(22)
     itab_v(jv(19))%im(3) = jm(14)
     itab_v(jv(19))%im(4) = jm(16)
     itab_v(jv(19))%im(5) = jm(21)
     itab_v(jv(19))%im(6) = jm(23)
     itab_v(jv(19))%iv(1) = jv(18)
     itab_v(jv(19))%iv(2) = jv(20)
     itab_v(jv(19))%iv(3) = jv(27)
     itab_v(jv(19))%iv(4) = jv(28)

     itab_v(jv(20))%ivp   = jv(14)
     itab_v(jv(20))%iw(1) = jw(4)
     itab_v(jv(20))%iw(2) = jw(8)
     itab_v(jv(20))%im(1) = jm(16)
     itab_v(jv(20))%im(2) = jm(15)

  enddo

  call wloopf('f',jaw( 2), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
  call wloopf('f',jaw( 3), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
  call wloopf('f',jaw( 4), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
  call wloopf('f',jzw( 9), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
  call wloopf('f',jzw(10), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
  call wloopf('f',jzw(11), jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)

  itab_w(jaw( 2))%iwp = jzw(6)
  itab_w(jaw( 3))%iwp = jzw(7)
  itab_w(jaw( 4))%iwp = jzw(6)
  itab_w(jzw( 9))%iwp = jaw(7)
  itab_w(jzw(10))%iwp = jaw(6)
  itab_w(jzw(11))%iwp = jaw(7)

  call mloopf('f',jam( 2), jtm_lbcp, 0, 0, 0, 0, 0)
  call mloopf('f',jam( 4), jtm_lbcp, 0, 0, 0, 0, 0)
  call mloopf('f',jam( 7), jtm_lbcp, 0, 0, 0, 0, 0)
  call mloopf('f',jzm(18), jtm_lbcp, 0, 0, 0, 0, 0)
  call mloopf('f',jzm(19), jtm_lbcp, jtm_grid, 0, 0, 0, 0)
  call mloopf('f',jzm(20), jtm_lbcp, jtm_grid, 0, 0, 0, 0)
  call mloopf('f',jzm(21), jtm_lbcp, jtm_grid, jtm_vadj, 0, 0, 0)
  call mloopf('f',jzm(22), jtm_lbcp, jtm_grid, jtm_vadj, 0, 0, 0)
  call mloopf('f',jzm(23), jtm_lbcp, 0, 0, 0, 0, 0)
  call mloopf('f',jzm(29), jtm_lbcp, 0, 0, 0, 0, 0)

  itab_m(jam(2))%npoly  = 1
  itab_m(jam(2))%imp    = jam(14)
  itab_m(jam(2))%iv(1)  = jav(2)

  itab_m(jam(4))%npoly  = 1
  itab_m(jam(4))%imp    = jzm(12)
  itab_m(jam(4))%iv(1)  = jav(5)

  itab_m(jam(7))%npoly  = 1
  itab_m(jam(7))%imp    = jzm(15)
  itab_m(jam(7))%iv(1)  = jav(9)

  itab_m(jzm(9))%imp    = jam(13)

  itab_m(jzm(10))%imp   = jam(14)

  itab_m(jzm(18))%npoly = 1
  itab_m(jzm(18))%imp   = jam(22)
  itab_m(jzm(18))%iv(1) = jzv(23)

  itab_m(jzm(19))%npoly = 3
  itab_m(jzm(19))%imp   = jam(15)
  itab_m(jzm(19))%iw(1) = jzw(9)
  itab_m(jzm(19))%iw(2) = jzw(6)
  itab_m(jzm(19))%iw(3) = jzw(10)
  itab_m(jzm(19))%iv(1) = jzv(24)
  itab_m(jzm(19))%iv(2) = jzv(23)
  itab_m(jzm(19))%iv(3) = jzv(12)

  itab_m(jzm(20))%npoly = 3
  itab_m(jzm(20))%imp   = jam(12)
  itab_m(jzm(20))%iw(1) = jzw(11)
  itab_m(jzm(20))%iw(2) = jzw(10)
  itab_m(jzm(20))%iw(3) = jzw(6)
  itab_m(jzm(20))%iv(1) = jzv(24)
  itab_m(jzm(20))%iv(2) = jzv(26)
  itab_m(jzm(20))%iv(3) = jzv(25)

  itab_m(jzm(21))%npoly = 3
  itab_m(jzm(21))%imp   = jam(13)
  itab_m(jzm(21))%iw(1) = jzw(6)
  itab_m(jzm(21))%iw(2) = jzw(7)
  itab_m(jzm(21))%iw(3) = jzw(11)
  itab_m(jzm(21))%iv(1) = jzv(27)
  itab_m(jzm(21))%iv(2) = jzv(26)
  itab_m(jzm(21))%iv(3) = jzv(15)

  itab_m(jzm(22))%npoly = 3
  itab_m(jzm(22))%imp   = jam(14)
  itab_m(jzm(22))%iw(1) = jzw(8)
  itab_m(jzm(22))%iw(2) = jzw(11)
  itab_m(jzm(22))%iw(3) = jzw(7)
  itab_m(jzm(22))%iv(1) = jzv(27)
  itab_m(jzm(22))%iv(2) = jzv(19)
  itab_m(jzm(22))%iv(3) = jzv(28)

  itab_m(jzm(23))%npoly = 1
  itab_m(jzm(23))%imp   = jam(15)
  itab_m(jzm(23))%iv(1) = jzv(28)

  itab_m(jzm(29))%npoly = 1
  itab_m(jzm(29))%imp   = jam(21)
  itab_m(jzm(29))%iv(1) = jzv(25)

  call vloopf('f',jav( 2), jtv_lbcp, 0, 0, 0, 0, 0)
  call vloopf('f',jav( 5), jtv_lbcp, 0, 0, 0, 0, 0)
  call vloopf('f',jav( 9), jtv_lbcp, 0, 0, 0, 0, 0)
  call vloopf('f',jzv(23), jtv_lbcp, 0, 0, 0, 0, 0)
  call vloopf('f',jzv(24), jtv_grid, jtv_lbcp, jtv_wadj, jtv_wstn, 0, 0)
  call vloopf('f',jzv(25), jtv_lbcp, 0, 0, 0, 0, 0)
  call vloopf('f',jzv(26), jtv_grid, jtv_lbcp, jtv_wadj, jtv_wstn, 0, 0)
  call vloopf('f',jzv(27), jtv_grid, jtv_lbcp, jtv_wadj, jtv_wstn, 0, 0)
  call vloopf('f',jzv(28), jtv_lbcp, 0, 0, 0, 0, 0)

  itab_v(jav(2))%ivp   = jav(18)
  itab_v(jav(2))%iw(1) = jaw(2)
  itab_v(jav(2))%iw(2) = jaw(5)
  itab_v(jav(2))%im(1) = jam(11)
  itab_v(jav(2))%im(2) = jam(2)

  itab_v(jav(5))%ivp   = jzv(15)
  itab_v(jav(5))%iw(1) = jaw(3)
  itab_v(jav(5))%iw(2) = jaw(2)
  itab_v(jav(5))%im(1) = jam(13)
  itab_v(jav(5))%im(2) = jam(4)

  itab_v(jav(9))%ivp   = jzv(19)
  itab_v(jav(9))%iw(1) = jaw(3)
  itab_v(jav(9))%iw(2) = jaw(4)
  itab_v(jav(9))%im(1) = jam(7)
  itab_v(jav(9))%im(2) = jam(14)

  itab_v(jzv(11))%ivp  = jav(17)

  itab_v(jzv(12))%ivp  = jav(18)

  itab_v(jzv(23))%ivp   = jav(19)
  itab_v(jzv(23))%iw(1) = jzw(9)
  itab_v(jzv(23))%iw(2) = jzw(10)
  itab_v(jzv(23))%im(1) = jzm(19)
  itab_v(jzv(23))%im(2) = jzm(18)

  itab_v(jzv(24))%ivp   = jav(14)
  itab_v(jzv(24))%iw(1) = jzw(6)
  itab_v(jzv(24))%iw(2) = jzw(10)
  itab_v(jzv(24))%iw(3) = jzw(11)
  itab_v(jzv(24))%iw(4) = jzw(9)
  itab_v(jzv(24))%im(1) = jzm(20)
  itab_v(jzv(24))%im(2) = jzm(19)
  itab_v(jzv(24))%im(3) = jzm(21)
  itab_v(jzv(24))%im(4) = jzm(29)
  itab_v(jzv(24))%im(5) = jzm(10)
  itab_v(jzv(24))%im(6) = jzm(18)
  itab_v(jzv(24))%iv(1) = jzv(26)
  itab_v(jzv(24))%iv(2) = jzv(25)
  itab_v(jzv(24))%iv(3) = jzv(12)
  itab_v(jzv(24))%iv(4) = jzv(23)

  itab_v(jzv(25))%ivp   = jav(15)
  itab_v(jzv(25))%iw(1) = jzw(11)
  itab_v(jzv(25))%iw(2) = jzw(10)
  itab_v(jzv(25))%im(1) = jzm(29)
  itab_v(jzv(25))%im(2) = jzm(20)

  itab_v(jzv(26))%ivp   = jav(16)
  itab_v(jzv(26))%iw(1) = jzw(6)
  itab_v(jzv(26))%iw(2) = jzw(11)
  itab_v(jzv(26))%iw(3) = jzw(7)
  itab_v(jzv(26))%iw(4) = jzw(10)
  itab_v(jzv(26))%im(1) = jzm(21)
  itab_v(jzv(26))%im(2) = jzm(20)
  itab_v(jzv(26))%im(3) = jzm(12)
  itab_v(jzv(26))%im(4) = jzm(22)
  itab_v(jzv(26))%im(5) = jzm(19)
  itab_v(jzv(26))%im(6) = jzm(29)
  itab_v(jzv(26))%iv(1) = jzv(15)
  itab_v(jzv(26))%iv(2) = jzv(27)
  itab_v(jzv(26))%iv(3) = jzv(24)
  itab_v(jzv(26))%iv(4) = jzv(25)

  itab_v(jzv(27))%ivp   = jav(17)
  itab_v(jzv(27))%iw(1) = jzw(7)
  itab_v(jzv(27))%iw(2) = jzw(11)
  itab_v(jzv(27))%iw(3) = jzw(8)
  itab_v(jzv(27))%iw(4) = jzw(6)
  itab_v(jzv(27))%im(1) = jzm(22)
  itab_v(jzv(27))%im(2) = jzm(21)
  itab_v(jzv(27))%im(3) = jzm(15)
  itab_v(jzv(27))%im(4) = jzm(23)
  itab_v(jzv(27))%im(5) = jzm(12)
  itab_v(jzv(27))%im(6) = jzm(20)
  itab_v(jzv(27))%iv(1) = jzv(19)
  itab_v(jzv(27))%iv(2) = jzv(28)
  itab_v(jzv(27))%iv(3) = jzv(15)
  itab_v(jzv(27))%iv(4) = jzv(26)

  itab_v(jzv(28))%ivp   = jav(18)
  itab_v(jzv(28))%iw(1) = jzw(8)
  itab_v(jzv(28))%iw(2) = jzw(11)
  itab_v(jzv(28))%im(1) = jzm(23)
  itab_v(jzv(28))%im(2) = jzm(22)

end subroutine cart4_hex
