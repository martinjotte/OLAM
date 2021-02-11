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
subroutine cartesian_2d()

  use mem_ijtabs, only: itab_md, itab_ud, itab_wd, mrls, alloc_itabsd, &
                        jtm_grid, jtu_grid, jtv_grid, jtw_grid, &
                        jtm_init, jtu_init, jtv_init, jtw_init, &
                        jtm_prog, jtu_prog, jtv_prog, jtw_prog, &
                        jtm_lbcp, jtu_lbcp, jtv_lbcp, jtw_lbcp, &
                        jtm_wadj, jtu_wadj, jtv_wadj, jtw_wadj, &
                        jtm_wstn, jtu_wstn, jtv_wstn, jtw_wstn, &
                        jtm_vadj, jtu_wall, jtv_wall, jtw_vadj

  use mem_grid,   only: nma, nua, nwa, xem, yem, zem, alloc_xyzem
  use misc_coms,  only: io6, nxp, mdomain, deltax

  implicit none

  integer :: i,im,iu,iw,im1,im2,im3,im4,iu1,iu2,iu3,iu4,iu5,iw0,iw1,iw2,iw3
  real :: unit_dist,diamond_centx

  ! Define triangles, edges, and vertices for 2D cartesian channel domain
  ! THIS CASE APPLIES FOR MDOMAIN = 3 (cyclic x bnds)

  mrls = 1  ! Default value

  ! Use nxp to count triangles

  nma = 2 * nxp + 2 + 1 ! ADDING 1 for reference point (index = 1)
  nua = 4 * nxp + 1 + 1 ! ADDING 1 for reference point (index = 1)
  nwa = 2 * nxp     + 1 ! ADDING 1 for reference point (index = 1)

  call alloc_itabsd(nma,nua,nwa)

  call alloc_xyzem(nma)

  do im = 2,nma
     itab_md(im)%imp = im
     call mdloopf('f',im, jtm_grid, jtm_vadj, 0, 0, 0, 0)
  enddo

  do i = 1,nxp

     im1 = 2 * i
     im2 = 2 * i + 1
     im3 = 2 * i + 2
     im4 = 2 * i + 3

     iu1 = 4 * i - 2
     iu2 = 4 * i - 1
     iu3 = 4 * i
     iu4 = 4 * i + 1
     iu5 = 4 * i + 2

     iw0 = 2 * i - 1
     iw1 = 2 * i
     iw2 = 2 * i + 1
     iw3 = 2 * i + 2

     if (i == 1  ) iw0 = 1
     if (i == nxp) iw3 = 1

     ! IU1

     itab_ud(iu1)%im(1) = im2
     itab_ud(iu1)%im(2) = im1
     itab_ud(iu1)%iw(1) = iw0
     itab_ud(iu1)%iw(2) = iw1

     if (i == 1) then
        itab_ud(iu1)%iup = 4*nxp-6
        call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu1)%iup = 6
        call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu1)%iup = iu1
        call udloopf('f',iu1, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU2

     itab_ud(iu2)%im(1) = im1
     itab_ud(iu2)%im(2) = im3
     itab_ud(iu2)%iw(1) = 1
     itab_ud(iu2)%iw(2) = iw1
     itab_ud(iu2)%iup = iu2

     call udloopf('f',iu2, jtu_wall, jtu_wadj, jtu_wstn, 0, 0, 0)

     ! IU3

     itab_ud(iu3)%im(1) = im2
     itab_ud(iu3)%im(2) = im3
     itab_ud(iu3)%iw(1) = iw1
     itab_ud(iu3)%iw(2) = iw2

     if (i == 1) then
        itab_ud(iu3)%iup = 4*nxp-4
        call udloopf('f',iu3, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu3)%iup = 8
        call udloopf('f',iu3, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu3)%iup = iu3
        call udloopf('f',iu3, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU4

     itab_ud(iu4)%im(1) = im2
     itab_ud(iu4)%im(2) = im4
     itab_ud(iu4)%iw(1) = iw2
     itab_ud(iu4)%iw(2) = 1
     itab_ud(iu4)%iup = iu4

     call udloopf('f',iu4, jtu_wall, jtu_wadj, jtu_wstn, 0, 0, 0)

     ! IU5

     if (i == nxp) then
        itab_ud(iu5)%im(1) = im4
        itab_ud(iu5)%im(2) = im3
        itab_ud(iu5)%iw(1) = iw2
        itab_ud(iu5)%iw(2) = iw3

        itab_ud(iu5)%iup = 10
        call udloopf('f',iu5, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IW1

     itab_wd(iw1)%iu(1) = iu1
     itab_wd(iw1)%iu(2) = iu2
     itab_wd(iw1)%iu(3) = iu3
     itab_wd(iw1)%mrlw = 1
     itab_wd(iw1)%mrlw_orig = 1
     itab_wd(iw1)%ngr = 1

     if (i == 1) then
        itab_wd(iw1)%iwp = 2*nxp-2
        call wdloopf('f',iw1, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw1)%iwp = 4
        call wdloopf('f',iw1, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw1)%iwp = iw1
        call wdloopf('f',iw1, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     ! IW2

     itab_wd(iw2)%iu(1) = iu3
     itab_wd(iw2)%iu(2) = iu5
     itab_wd(iw2)%iu(3) = iu4
     itab_wd(iw2)%mrlw = 1
     itab_wd(iw2)%mrlw_orig = 1
     itab_wd(iw2)%ngr = 1
     if (i == 1) then
        itab_wd(iw2)%iwp = 2*nxp-1
        call wdloopf('f',iw2, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw2)%iwp = 5
        call wdloopf('f',iw2, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw2)%iwp = iw2
        call wdloopf('f',iw2, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     unit_dist = .5 * sqrt(3.) * deltax ! This is 1/2 of triangle face width
     diamond_centx = (2 * i - 1 - nxp) * unit_dist

     xem(im1) = diamond_centx - 1.5 * unit_dist
     xem(im2) = diamond_centx -  .5 * unit_dist
     xem(im3) = diamond_centx +  .5 * unit_dist
     xem(im4) = diamond_centx + 1.5 * unit_dist

     yem(im1) = - .75 * deltax
     yem(im2) =   .75 * deltax
     yem(im3) = - .75 * deltax
     yem(im4) =   .75 * deltax

  enddo  ! end i loop

  zem(2:nma) = 0.

  call tri_neighbors(nma, nua, nwa, itab_md, itab_ud, itab_wd)

end subroutine cartesian_2d

!=============================================================================

subroutine cartesian_3d()

  use mem_ijtabs, only: itab_md, itab_ud, itab_wd, mrls, alloc_itabsd, &
                        jtm_grid, jtu_grid, jtw_grid, &
                        jtu_init, jtu_prog, jtu_wadj, jtu_wstn, jtu_lbcp, &
                        jtw_init, jtw_prog, jtw_wadj, jtw_wstn, jtw_lbcp

  use mem_grid,   only: nma, nua, nwa, xem, yem, zem, alloc_xyzem

  use misc_coms,  only: io6, nxp, mdomain, deltax

  implicit none

  integer :: i,im,iu,iw,iucyc,iwcyc
  integer :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10
  integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10,iu11,iu12,iu13,iu14,iu15
  integer :: iw1,iw2,iw3,iw4,iw5,iw6,iw7,iw8,iw9,iw10,iw11,iw12,iw13,iw14,iw15
  integer :: iu16,iu17,iw16

  real :: unit_dist,diamond_centx

  ! Define triangles, edges, and vertices for 3D cartesian channel domain
  ! Channel width = 4 allocated rows with cyclic repeat distance 2 rows wide
  ! THIS CASE APPLIES FOR MDOMAIN = 4 (cyclic x bnds & cyclic y bnds)

  mrls = 1  ! Default value

  ! Use nxp to count triangles

  nma =  5 * nxp + 5 + 1 ! ADDING 1 for reference point (index = 1)
  nua = 13 * nxp + 4 + 1 ! ADDING 1 for reference point (index = 1)
  nwa =  8 * nxp     + 1 ! ADDING 1 for reference point (index = 1)

  call alloc_itabsd(nma,nua,nwa)

  call alloc_xyzem(nma)

  do im = 2,nma
     itab_md(im)%imp = im
     call mdloopf('f',im, jtm_grid, 0, 0, 0, 0, 0)
  enddo

  do i = 1,nxp

     im1  = 5 * i - 3
     im2  = 5 * i - 2
     im3  = 5 * i - 1
     im4  = 5 * i
     im5  = 5 * i + 1
     im6  = 5 * i + 2
     im7  = 5 * i + 3
     im8  = 5 * i + 4
     im9  = 5 * i + 5
     im10 = 5 * i + 6

     iu1  = 13 * i - 11
     iu2  = 13 * i - 10
     iu3  = 13 * i - 9
     iu4  = 13 * i - 8
     iu5  = 13 * i - 7
     iu6  = 13 * i - 6
     iu7  = 13 * i - 5
     iu8  = 13 * i - 4
     iu9  = 13 * i - 3
     iu10 = 13 * i - 2
     iu11 = 13 * i - 1
     iu12 = 13 * i
     iu13 = 13 * i + 1
     iu14 = 13 * i + 2
     iu15 = 13 * i + 3
     iu16 = 13 * i + 4
     iu17 = 13 * i + 5

     iw1  = 8 * i - 10
     iw2  = 8 * i - 9
     iw3  = 8 * i - 8
     iw4  = 8 * i - 7
     iw5  = 8 * i - 6
     iw6  = 8 * i - 5
     iw7  = 8 * i - 4
     iw8  = 8 * i - 3
     iw9  = 8 * i - 2
     iw10 = 8 * i - 1
     iw11 = 8 * i
     iw12 = 8 * i + 1
     iw13 = 8 * i + 2
     iw14 = 8 * i + 3
     iw15 = 8 * i + 4
     iw16 = 8 * i + 5

     if (i == 1) then
        iw1 = 1
        iw2 = 1
        iw3 = 1
        iw4 = 1
     endif
   
     if (i == nxp) then
        iw13 = 1
        iw14 = 1
        iw15 = 1
        iw16 = 1
     endif

     ! Cyclic repeat distance along channel

     iucyc = 13 * (nxp - 2)
     iwcyc =  8 * (nxp - 2)

     ! IU1

     itab_ud(iu1)%im(1) = im4
     itab_ud(iu1)%im(2) = im1
     itab_ud(iu1)%iw(1) = iw1
     itab_ud(iu1)%iw(2) = iw5

     if (i == 1) then
        itab_ud(iu1)%iup = iu3+iucyc
        call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu1)%iup = iu3-iucyc
        call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu1)%iup = iu3
        call udloopf('f',iu1, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU2

     itab_ud(iu2)%im(1) = im2
     itab_ud(iu2)%im(2) = im4
     itab_ud(iu2)%iw(1) = iw2
     itab_ud(iu2)%iw(2) = iw6

     if (i == 1) then
        itab_ud(iu2)%iup = iu2+iucyc
        call udloopf('f',iu2, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu2)%iup = iu2-iucyc
        call udloopf('f',iu2, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu2)%iup = iu2
        call udloopf('f',iu2, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU3

     itab_ud(iu3)%im(1) = im5
     itab_ud(iu3)%im(2) = im2
     itab_ud(iu3)%iw(1) = iw3
     itab_ud(iu3)%iw(2) = iw7

     if (i == 1) then
        itab_ud(iu3)%iup = iu3+iucyc
        call udloopf('f',iu3, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu3)%iup = iu3-iucyc
        call udloopf('f',iu3, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu3)%iup = iu3
        call udloopf('f',iu3, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU4

     itab_ud(iu4)%im(1) = im3
     itab_ud(iu4)%im(2) = im5
     itab_ud(iu4)%iw(1) = iw4
     itab_ud(iu4)%iw(2) = iw8

     if (i == 1) then
        itab_ud(iu4)%iup = iu2+iucyc
        call udloopf('f',iu4, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu4)%iup = iu2-iucyc
        call udloopf('f',iu4, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu4)%iup = iu2
        call udloopf('f',iu4, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU5

     itab_ud(iu5)%im(1) = im1
     itab_ud(iu5)%im(2) = im6
     itab_ud(iu5)%iw(1) = 1
     itab_ud(iu5)%iw(2) = iw5

     if (i == 1) then
        itab_ud(iu5)%iup = iu6+iucyc
        call udloopf('f',iu5, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu5)%iup = iu6-iucyc
        call udloopf('f',iu5, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu5)%iup = iu6
        call udloopf('f',iu5, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU6

     itab_ud(iu6)%im(1) = im2
     itab_ud(iu6)%im(2) = im7
     itab_ud(iu6)%iw(1) = iw6
     itab_ud(iu6)%iw(2) = iw7

     if (i == 1) then
        itab_ud(iu6)%iup = iu6+iucyc
        call udloopf('f',iu6, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu6)%iup = iu6-iucyc
        call udloopf('f',iu6, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu6)%iup = iu6
        call udloopf('f',iu6, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU7

     itab_ud(iu7)%im(1) = im3
     itab_ud(iu7)%im(2) = im8
     itab_ud(iu7)%iw(1) = iw8
     itab_ud(iu7)%iw(2) = 1

     if (i == 1) then
        itab_ud(iu7)%iup = iu6+iucyc
        call udloopf('f',iu7, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu7)%iup = iu6-iucyc
        call udloopf('f',iu7, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu7)%iup = iu6
        call udloopf('f',iu7, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

  ! IU8

     itab_ud(iu8)%im(1) = im4
     itab_ud(iu8)%im(2) = im6
     itab_ud(iu8)%iw(1) = iw5
     itab_ud(iu8)%iw(2) = iw9

     if (i == 1) then
        itab_ud(iu8)%iup = iu10+iucyc
        call udloopf('f',iu8, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu8)%iup = iu10-iucyc
        call udloopf('f',iu8, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu8)%iup = iu10
        call udloopf('f',iu8, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU9

     itab_ud(iu9)%im(1) = im7
     itab_ud(iu9)%im(2) = im4
     itab_ud(iu9)%iw(1) = iw6
     itab_ud(iu9)%iw(2) = iw10

     if (i == 1) then
        itab_ud(iu9)%iup = iu9+iucyc
        call udloopf('f',iu9, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu9)%iup = iu9-iucyc
        call udloopf('f',iu9, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu9)%iup = iu9
        call udloopf('f',iu9, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif 

     ! IU10

     itab_ud(iu10)%im(1) = im5
     itab_ud(iu10)%im(2) = im7
     itab_ud(iu10)%iw(1) = iw7
     itab_ud(iu10)%iw(2) = iw11

     if (i == 1) then
        itab_ud(iu10)%iup = iu10+iucyc
        call udloopf('f',iu10, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu10)%iup = iu10-iucyc
        call udloopf('f',iu10, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu10)%iup = iu10
        call udloopf('f',iu10, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU11

     itab_ud(iu11)%im(1) = im8
     itab_ud(iu11)%im(2) = im5
     itab_ud(iu11)%iw(1) = iw8
     itab_ud(iu11)%iw(2) = iw12

     if (i == 1) then
        itab_ud(iu11)%iup = iu9+iucyc
        call udloopf('f',iu11, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu11)%iup = iu9-iucyc
        call udloopf('f',iu11, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu11)%iup = iu9
        call udloopf('f',iu11, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU12

     itab_ud(iu12)%im(1) = im4
     itab_ud(iu12)%im(2) = im9
     itab_ud(iu12)%iw(1) = iw9
     itab_ud(iu12)%iw(2) = iw10

     if (i == 1) then
        itab_ud(iu12)%iup = iu12+iucyc
        call udloopf('f',iu12, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu12)%iup = iu12-iucyc
        call udloopf('f',iu12, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu12)%iup = iu12
        call udloopf('f',iu12, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     endif

     ! IU13

     itab_ud(iu13)%im(1) = im5
     itab_ud(iu13)%im(2) = im10
     itab_ud(iu13)%iw(1) = iw11
     itab_ud(iu13)%iw(2) = iw12

     if (i == 1) then
        itab_ud(iu13)%iup = iu12+iucyc
        call udloopf('f',iu13, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_ud(iu13)%iup = iu12-iucyc
        call udloopf('f',iu13, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     else
        itab_ud(iu13)%iup = iu12
        call udloopf('f',iu13, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU14

     if (i == nxp) then
        itab_ud(iu14)%im(1) = im9
        itab_ud(iu14)%im(2) = im6
        itab_ud(iu14)%iw(1) = iw9
        itab_ud(iu14)%iw(2) = iw13

        itab_ud(iu14)%iup = iu16-iucyc
        call udloopf('f',iu14, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU15

     if (i == nxp) then
        itab_ud(iu15)%im(1) = im7
        itab_ud(iu15)%im(2) = im9
        itab_ud(iu15)%iw(1) = iw10
        itab_ud(iu15)%iw(2) = iw14

        itab_ud(iu15)%iup = iu15-iucyc
        call udloopf('f',iu15, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU16

     if (i == nxp) then
        itab_ud(iu16)%im(1) = im10
        itab_ud(iu16)%im(2) = im7
        itab_ud(iu16)%iw(1) = iw11
        itab_ud(iu16)%iw(2) = iw15

        itab_ud(iu16)%iup = iu13-iucyc
        call udloopf('f',iu16, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IU17

     if (i == nxp) then
        itab_ud(iu17)%im(1) = im8
        itab_ud(iu17)%im(2) = im10
        itab_ud(iu17)%iw(1) = iw12
        itab_ud(iu17)%iw(2) = iw16

        itab_ud(iu17)%iup = iu15-iucyc
        call udloopf('f',iu17, jtu_lbcp, jtu_wadj, jtu_wstn, 0, 0, 0)
     endif

     ! IW5

     itab_wd(iw5)%iu(1) = iu1
     itab_wd(iw5)%iu(2) = iu5
     itab_wd(iw5)%iu(3) = iu8
     itab_wd(iw5)%mrlw = 1
     itab_wd(iw5)%mrlw_orig = 1
     itab_wd(iw5)%ngr = 1

     if (i == 1) then
        itab_wd(iw5)%iwp = iw7+iwcyc
        call wdloopf('f',iw5, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw5)%iwp = iw7-iwcyc
        call wdloopf('f',iw5, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw5)%iwp = iw7
        call wdloopf('f',iw5, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     endif

     ! IW6

     itab_wd(iw6)%iu(1) = iu2
     itab_wd(iw6)%iu(2) = iu9
     itab_wd(iw6)%iu(3) = iu6
     itab_wd(iw6)%mrlw = 1
     itab_wd(iw6)%mrlw_orig = 1
     itab_wd(iw6)%ngr = 1

     if (i == 1) then
        itab_wd(iw6)%iwp = iw6+iwcyc
        call wdloopf('f',iw6, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw6)%iwp = iw6-iwcyc
        call wdloopf('f',iw6, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw6)%iwp = iw6
        call wdloopf('f',iw6, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     ! IW7

     itab_wd(iw7)%iu(1) = iu3
     itab_wd(iw7)%iu(2) = iu6
     itab_wd(iw7)%iu(3) = iu10
     itab_wd(iw7)%mrlw = 1
     itab_wd(iw7)%mrlw_orig = 1
     itab_wd(iw7)%ngr = 1

     if (i == 1) then
        itab_wd(iw7)%iwp = iw7+iwcyc
        call wdloopf('f',iw7, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw7)%iwp = iw7-iwcyc
        call wdloopf('f',iw7, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw7)%iwp = iw7
        call wdloopf('f',iw7, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     ! IW8

     itab_wd(iw8)%iu(1) = iu4
     itab_wd(iw8)%iu(2) = iu11
     itab_wd(iw8)%iu(3) = iu7
     itab_wd(iw8)%mrlw = 1
     itab_wd(iw8)%mrlw_orig = 1
     itab_wd(iw8)%ngr = 1

     if (i == 1) then
        itab_wd(iw8)%iwp = iw6+iwcyc
        call wdloopf('f',iw8, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw8)%iwp = iw6-iwcyc
        call wdloopf('f',iw8, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw8)%iwp = iw6
        call wdloopf('f',iw8, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     endif

     ! IW9

     itab_wd(iw9)%iu(1) = iu8
     itab_wd(iw9)%iu(2) = iu14
     itab_wd(iw9)%iu(3) = iu12
     itab_wd(iw9)%mrlw = 1
     itab_wd(iw9)%mrlw_orig = 1
     itab_wd(iw9)%ngr = 1

     if (i == 1) then
        itab_wd(iw9)%iwp = iw11+iwcyc
        call wdloopf('f',iw9, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw9)%iwp = iw11-iwcyc
        call wdloopf('f',iw9, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw9)%iwp = iw11
        call wdloopf('f',iw9, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     endif

     ! IW10

     itab_wd(iw10)%iu(1) = iu9
     itab_wd(iw10)%iu(2) = iu12
     itab_wd(iw10)%iu(3) = iu15
     itab_wd(iw10)%mrlw = 1
     itab_wd(iw10)%mrlw_orig = 1
     itab_wd(iw10)%ngr = 1

     if (i == 1) then
        itab_wd(iw10)%iwp = iw10+iwcyc
        call wdloopf('f',iw10, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw10)%iwp = iw10-iwcyc
        call wdloopf('f',iw10, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw10)%iwp = iw10
        call wdloopf('f',iw10, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     ! IW11

     itab_wd(iw11)%iu(1) = iu10
     itab_wd(iw11)%iu(2) = iu16
     itab_wd(iw11)%iu(3) = iu13
     itab_wd(iw11)%mrlw = 1
     itab_wd(iw11)%mrlw_orig = 1
     itab_wd(iw11)%ngr = 1

     if (i == 1) then
        itab_wd(iw11)%iwp = iw11+iwcyc
        call wdloopf('f',iw11, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw11)%iwp = iw11-iwcyc
        call wdloopf('f',iw11, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw11)%iwp = iw11
        call wdloopf('f',iw11, jtw_grid, jtw_init, jtw_prog, jtw_wadj, jtw_wstn, 0)
     endif

     ! IW12

     itab_wd(iw12)%iu(1) = iu11
     itab_wd(iw12)%iu(2) = iu13
     itab_wd(iw12)%iu(3) = iu17
     itab_wd(iw12)%mrlw = 1
     itab_wd(iw12)%mrlw_orig = 1
     itab_wd(iw12)%ngr = 1

     if (i == 1) then
        itab_wd(iw12)%iwp = iw10+iwcyc
        call wdloopf('f',iw12, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     elseif (i == nxp) then
        itab_wd(iw12)%iwp = iw10-iwcyc
        call wdloopf('f',iw12, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     else
        itab_wd(iw12)%iwp = iw10
        call wdloopf('f',iw12, jtw_lbcp, jtw_wadj, jtw_wstn, 0, 0, 0)
     endif

     unit_dist = .5 * sqrt(3.) * deltax ! This is 1/2 of triangle face width
     diamond_centx = (2 * i - 1 - nxp) * unit_dist

     xem(im1)  = diamond_centx - 1.5 * unit_dist
     xem(im2)  = xem(im1)
     xem(im3)  = xem(im1)
     xem(im4)  = diamond_centx -  .5 * unit_dist
     xem(im5)  = xem(im4)
     xem(im6)  = diamond_centx +  .5 * unit_dist
     xem(im7)  = xem(im6)
     xem(im8)  = xem(im6)
     xem(im9)  = diamond_centx + 1.5 * unit_dist
     xem(im10) = xem(im9)

     yem(im1)  = - 3.0 * deltax
     yem(im2)  =    .0 * deltax
     yem(im3)  =   3.0 * deltax
     yem(im4)  = - 1.5 * deltax
     yem(im5)  =   1.5 * deltax
     yem(im6)  = - 3.0 * deltax
     yem(im7)  =    .0 * deltax
     yem(im8)  =   3.0 * deltax
     yem(im9)  = - 1.5 * deltax
     yem(im10) =   1.5 * deltax

  enddo  ! end i loop

  zem(2:nma) = 0.

  call tri_neighbors(nma, nua, nwa, itab_md, itab_ud, itab_wd)

end subroutine cartesian_3d

